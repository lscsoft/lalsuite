#include "burstdso.h"
#include <lal/FrameStream.h>
#include <lalapps/lalapps.h>
#include "Translation.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <FrameL.h>
#include <unistd.h>

RCSID( "getFrameData.c" );
NRCSID( GETFRAMEDATAC, "getFrameData.c" );

/*****************************************/
/* acquire data */

/* loop over lines in fQuery, get data, add to datacond callchain */
/* Format is:
   frame_type Site [file] GPS_Start-GPS_End channel alias
   file is optional, and is a frame file that overrides LSCdataFind
   channel can be adc(H2:LSC-AS_Q), proc(H2:CAL-OLOOP_FAC), or proc(H2:CAL-RESPONSE!0!7000.0001!), for instance; same format as LDAS
   alias is a name to represent the data in the datacondAPI
*/
/*****************************************/

int getFrameData(char *fQuery, 
		 char *dataserver,
		 int *Nsymbols,
		 datacond_symbol_type **symbols) {

  /*
    Arguments:
    fQuery: string containing the frame queries
    dataserver: server for LSCdataFind, or frame cache file
    Nsymbols: number of datacondAPI symbols
    symbols: datacondAPI symbols
  */

  /* NOTE: symbols->name and symbols->s_user_data must be freed */

  static LALStatus      stat;  /* LAL status pointer */

  char *p0, *p1;   /* pointers for string manipulation */
  char *buf;       /* buffer for local copy of fQuery */

  buf = (char *)calloc(1+strlen(fQuery),sizeof(char));

  /* process line by line */

  /* 1st line: */
  p0 = fQuery;
  p1 = strchr(fQuery,'\n');
  if(p1) {
    memcpy(buf,p0,p1-p0+1);
  } else {
    strcpy(buf,p0);
  }

  /* loop as long as line contains stuff: */
  while(strlen(buf)) {

    char type[256], IFO[256], channel[256], Times[256], File[1024]; /* fields */
    char *alias = (char *)calloc(256, sizeof(char)); /* alias name */
    char *q, *T0, *T1;  /* for string manipulation */

    /* data of different types: */
    REAL4TimeSeries *inputR4;
    REAL8TimeSeries *inputR8;
    COMPLEX8TimeSeries *inputC8;
    COMPLEX16TimeSeries *inputC16;

    REAL4FrequencySeries *inputFR4;
    REAL8FrequencySeries *inputFR8;
    COMPLEX8FrequencySeries *inputFC8;
    COMPLEX16FrequencySeries *inputFC16;
    
    LALTYPECODE otype;   /* LAL type */
    ObjectType stype;    /* frame object type */

    LIGOTimeGPS endEpoch; 

    double start = 0.0, end = -1.0;  /* for name!lower!upper selection */

    if(buf[0]!='\n') {

#ifdef DEBUGBURST
      fprintf(stderr,"Processing: %s\n", buf);
#endif

      /* Examine line
	 Format: type IFO times channel alias
	 or type IFO File times channel alias
      */
      if(sscanf(buf,"%s\t%s\t%s\t%s\t%s\t%s",type,IFO,File,Times,channel,alias) != 6) {
	if(sscanf(buf,"%s\t%s\t%s\t%s\t%s",type,IFO,Times,channel,alias) != 5) {	
	  fprintf(stderr,"Malformed frame query: %s\n",buf);
	  return 1;
	}
      }
    
      /* split time interval in start and end times */
      q = strchr(Times,'-');
      if(!q) {
	fprintf(stderr,"Times: T0-T1\n");
	return 1;
      } 
      *q=0;
      T0 = Times;
      T1 = q+1;

      /* get data */
      {
	char tname[] = "/tmp/getFrameDataXXXXXX";   /* temporary filename */
	char cmd[1024];  /* shell command */
	char *path;      /* path */
	char *p;
	FrCache              *frameCache    = NULL;  /* LAL frame cache */
	FrStream             *stream        = NULL;  /* LAL frame stream */
	FrChanIn             ChannelIn;              /* LAL channel info */

	/* build channel type: adc, proc, sim */
	if((p=strstr(channel,"adc("))) {
	  ChannelIn.name = p + 4; /* channel name */
	  p=strchr(p,')');
	  *p=0;
	  ChannelIn.type = ADCDataChannel;  /* channel type */
	} else if((p=strstr(channel,"proc("))) {
	  ChannelIn.name = p + 5;
	  p=strchr(p,')');
	  *p=0;
	  ChannelIn.type = ProcDataChannel;
	} else if((p=strstr(channel,"sim("))) {
	  ChannelIn.name = p + 4;
	  p=strchr(p,')');
	  *p=0;
	  ChannelIn.type = SimDataChannel;
	} else {
	  fprintf(stderr,"Unsupported channel type\n");
	  return 1;
	}

	/* check for subset request with name!low!high notation */
	if((p=strchr(ChannelIn.name,'!'))) {
	  
	  /* Channel name contains a ! */

	  char *q;  /* for string manipulation */

	  /* parse lower bound */
	  *p=0;
	  p=p+1;
	  q=strchr(p,'!');
	  if(!q) {
	    fprintf(stderr,"Malformed channel name\n");
	    return 1;
	  }

	  *q=0;
	  start = atof(p);

	  /* parse upper bound */
	  p=q+1;
	  q=strchr(p,'!');
	  if(!q) {
	    fprintf(stderr,"Malformed channel name\n");
	    return 1;
	  }

	  *q=0;
	  end = atof(p);

	  /* check it makes sense */
	  if(end<=start) {
	    fprintf(stderr,"Malformed channel subset request\n");
	    return 1;
	  }
	}
	
	/* find the data */
	if(!strstr(dataserver,CACHEFILENAME)) {
	  /* use LSCdataFind since dataserver doesn't contain
	     the string defined in CACHEFILENAME */
	  {
	    int fid = mkstemp(tname); /* create unique temporary file */
	    if(fid==-1) {
	      fprintf(stderr,"Can't create temp file\n");
	      return 1;
	    }
	    close(fid);
	  }
	
	  /* get executable location from environment variable */
	  path = getenv("LSC_DATAGRID_CLIENT_LOCATION");
	  if(!path) {
	    fprintf(stderr,"Environment variable LSC_DATAGRID_CLIENT_LOCATION not set\n");
	    return 1;
	  }
      
	  /* construct LSCdataFind command */
	  /* send output to temporary file tname */
	  sprintf(cmd,"source %s/setup.sh; %s/ldg-client/bin/LSCdataFind --server %s --observatory %s --type %s --gps-start-time %s --gps-end-time %s --url-type file --lal-cache > %s", path, path, dataserver, IFO, type, T0, T1, tname);

	  /* run & check success of shell command */
	  if(system(cmd) == -1) {
	    fprintf(stderr,"system call failed\n");
	    perror("Error");
	    unlink(tname);
	    return 1;
	  }

	  /* get the data with LAL functions */
	  LAL_CALL( LALFrCacheImport( &stat, &frameCache, tname ), &stat );
	  LAL_CALL( LALFrCacheOpen( &stat, &stream, frameCache ), &stat );
	  LAL_CALL( LALDestroyFrCache( &stat, &frameCache ), &stat );
	
	  unlink(tname); /* remove temporary file */

	} else {
	  /* we have a cache file */
	  /* dataserver thus points to file name */

	  char tname[] = "/tmp/getFrameDataXXXXXX"; /* temporary file name */
	  char cmd[2048];  /* shell command */
	  int fid = mkstemp(tname); /* create and open temp file */
	  char gotIt = 0;

	  if(fid==-1) {
	    fprintf(stderr,"Can't create temp file\n");
	    return 1;
	  }

	  /* must "grep" for the right frame type */
	  {
	    FILE *in;
	    char buf[2048];

	    /* open cache file */
	    if((in=fopen(dataserver,"r"))==NULL) {
	      fprintf(stderr,"Can't open %s\n",dataserver);
	      return 1;
	    }
	    
	    /* loop over file; save in tempfile lines with
	       the right frame type */
	    while(fgets(buf,2048,in)) {
	      if(strstr(buf,type)) {
		gotIt = 1;
		write(fid,buf,strlen(buf));
	      }
	    }

	    /* close input file */
	    fclose(in);
	  }

	  /* close output file */
	  *buf = '\n';
	  write(fid,buf,1);
	  close(fid);

	  /* check for at least one match */
	  if(!gotIt) {
	    fprintf(stderr,"No frame file matches request for type %s\n",type);
	    return 1;
	  }

	  /* get the data with LAL functions */
	  LAL_CALL( LALFrCacheImport( &stat, &frameCache, tname ), &stat );
	  LAL_CALL( LALFrCacheOpen( &stat, &stream, frameCache ), &stat );
	  LAL_CALL( LALDestroyFrCache( &stat, &frameCache ), &stat );

	  /* remove temp file */
	  unlink(tname);
	}


	/* obtain the type of the series */
	LAL_CALL( LALFrGetSeriesType(&stat, &otype, &stype, &ChannelIn, stream), &stat);

	switch ( stype ) {   /* switch between time and frequency series */

	case TimeSeries:

	  switch ( otype ) {  /* switch between datatypes */

	  case LAL_S_TYPE_CODE:  /* REAL4 */

	    /* find data in frame: */
	    inputR4 = (REAL4TimeSeries *)calloc(1,sizeof(REAL4TimeSeries));
	    inputR4->data = NULL;

	    /* get series info */
	    LAL_CALL( LALFrGetREAL4TimeSeries( &stat, inputR4, &ChannelIn, stream),&stat); 

	    /* set requested start: */
	    inputR4->epoch.gpsSeconds     = (INT4)floor(atof(T0));
	    inputR4->epoch.gpsNanoSeconds = (INT4)(1e9*(atof(T0) - floor(atof(T0))));
	    
	    /* look for data in frame: */
	    LAL_CALL( LALFrSeek(&stat, &(inputR4->epoch), stream),&stat);

	    {
	      /* allocate memory: */

	      double t0 = (double)(inputR4->epoch.gpsSeconds) + 1e-9*(double)(inputR4->epoch.gpsNanoSeconds); /* start time */
	      INT4 length = (INT4)floor((atof(T1) - t0)/inputR4->deltaT); /* number of data points */

	      if(inputR4->data) { /* destroy whatever is there already */
		LAL_CALL(LALDestroyVector(&stat,&(inputR4->data)),&stat);
		inputR4->data = NULL;
	      }

	      /* allocate memory: */
	      LAL_CALL ( LALCreateVector( &stat, &(inputR4->data), length), &stat);
	    }
	    
	    /* get the data */
	    LAL_CALL( LALFrGetREAL4TimeSeries( &stat, inputR4, &ChannelIn, stream),&stat);
	    
	    /* datacond convention: f0 = 1/sampling rate */
	    inputR4->f0 = 1.0/inputR4->deltaT;

	    {
	      /* Calculate end */

	      int iend; /* end index */

	      iend = (int)floor((atof(T1) - (double)(inputR4->epoch.gpsSeconds) + 1e-9*(double)(inputR4->epoch.gpsNanoSeconds))/inputR4->deltaT);

	      iend = (iend>0) ? iend : 1; /* make sure have at least one datapoint */
	      inputR4->data->length = iend;
	    }

	    /* truncate and shift if necessary */
	    if(end>0) {
	      int ind, sind, eind;
	      sind = (int)floor(start / inputR4->deltaT); /* start index */
	      eind = (int)floor((start+end) / inputR4->deltaT); /* end index */
	      eind = (eind > inputR4->data->length) ? inputR4->data->length : eind;

	      /* shift */
	      for(ind=sind; ind<eind; ind++) {
		inputR4->data->data[ind-sind] = inputR4->data->data[ind];
	      }

	      /* reset length */
	      inputR4->data->length = eind - sind;

	    }

	    break; /* end of REAL4 case */
	    
	  case LAL_D_TYPE_CODE:
	    inputR8 = (REAL8TimeSeries *)calloc(1,sizeof(REAL8TimeSeries));
	    inputR8->data = NULL;
	    LAL_CALL( LALFrGetREAL8TimeSeries( &stat, inputR8, &ChannelIn, stream),&stat);
	    inputR8->epoch.gpsSeconds     = (INT4)floor(atof(T0));
	    inputR8->epoch.gpsNanoSeconds = (INT4)(1e9*(atof(T0) - floor(atof(T0))));
	    LAL_CALL( LALFrSeek(&stat, &(inputR8->epoch), stream),&stat);
	    /* allocate memory */
	    {
	      double t0 = (double)(inputR8->epoch.gpsSeconds) + 1e-9*(double)(inputR8->epoch.gpsNanoSeconds);
	      INT4 length = (INT4)floor((atof(T1) - t0)/inputR8->deltaT);

	      if(inputR8->data) {
		LAL_CALL(LALDDestroyVector(&stat,&(inputR8->data)),&stat);
		inputR8->data = NULL;
	      }

	      LAL_CALL ( LALDCreateVector( &stat, &(inputR8->data), length), &stat);
	    }
	    
	    /* get the data */
	    LAL_CALL( LALFrGetREAL8TimeSeries( &stat, inputR8, &ChannelIn, stream),&stat);
	    
	    /* datacond convention: f0 = 1/sampling rate */
	    inputR8->f0 = 1.0/inputR8->deltaT;
	    
	    /* Calculate end */
	    {
	      int iend;
	      iend = (int)floor((atof(T1) - (double)(inputR8->epoch.gpsSeconds) + 1e-9*(double)(inputR8->epoch.gpsNanoSeconds))/inputR8->deltaT);
	      iend = (iend>0) ? iend : 1;
	      inputR8->data->length = iend;
	    }

	    /* truncate and shift if necessary */
	    if(end>0) {
	      int ind, sind, eind;
	      sind = (int)floor(start / inputR8->deltaT);
	      eind = (int)floor((start+end) / inputR8->deltaT);
	      eind = (eind > inputR8->data->length) ? inputR8->data->length : eind;

	      for(ind=sind; ind<eind; ind++) {
		inputR8->data->data[ind-sind] = inputR8->data->data[ind];
	      }

	      inputR8->data->length = eind - sind;

	    }

	    break;
	    
	  case LAL_C_TYPE_CODE:
	    inputC8 = (COMPLEX8TimeSeries *)calloc(1,sizeof(COMPLEX8TimeSeries));
	    inputC8->data = NULL;
	    LAL_CALL( LALFrGetCOMPLEX8TimeSeries( &stat, inputC8, &ChannelIn, stream),&stat);
	    inputC8->epoch.gpsSeconds     = (INT4)floor(atof(T0));
	    inputC8->epoch.gpsNanoSeconds = (INT4)(1e9*(atof(T0) - floor(atof(T0))));
	    LAL_CALL( LALFrSeek(&stat, &(inputC8->epoch), stream),&stat);
	    /* allocate memory */
	    {
	      double t0 = (double)(inputC8->epoch.gpsSeconds) + 1e-9*(double)(inputC8->epoch.gpsNanoSeconds);
	      INT4 length = (INT4)floor((atof(T1) - t0)/inputC8->deltaT);

	      if(inputC8->data) {
		LAL_CALL(LALCDestroyVector(&stat,&(inputC8->data)),&stat);
		inputC8->data = NULL;
	      }

	      LAL_CALL ( LALCCreateVector( &stat, &(inputC8->data), length), &stat);
	    }
	    
	    /* get the data */
	    LAL_CALL( LALFrGetCOMPLEX8TimeSeries( &stat, inputC8, &ChannelIn, stream),&stat);
	    
	    /* datacond convention: f0 = 1/sampling rate */
	    inputC8->f0 = 1.0/inputC8->deltaT;
	    
	    /* Calculate end */
	    {
	      int iend;
	      iend = (int)floor((atof(T1) - (double)(inputC8->epoch.gpsSeconds) + 1e-9*(double)(inputC8->epoch.gpsNanoSeconds))/inputC8->deltaT);
	      iend = (iend>0) ? iend : 1;
	      inputC8->data->length = iend;
	    }

	    /**********************************************************/
	    /**********************************************************/
	    /**********************************************************/
	    /* To get same results as with LDAS: */
	    /**********************************************************/
	    /**********************************************************/
	    /**********************************************************/

	    if((strstr(ChannelIn.name,":CAL-CAV_FAC") ||
	       strstr(ChannelIn.name,":CAL-OLOOP_FAC")) &&
	       inputC8->epoch.gpsSeconds != (INT4)floor(atof(T0)))
	      {
		int i;
		for(i=0;i<inputC8->data->length-1;i++) {
		  inputC8->data->data[i] = inputC8->data->data[i+1];
		}
	      }

	    /**********************************************************/
	    /**********************************************************/
	    /**********************************************************/

	    /* truncate and shift if necessary */
	    if(end>0) {
	      int ind, sind, eind;
	      sind = (int)floor(start / inputC8->deltaT);
	      eind = (int)floor((start+end) / inputC8->deltaT);

	      eind = (eind > inputC8->data->length) ? inputC8->data->length : eind;

	      for(ind=sind; ind<eind; ind++) {
		inputC8->data->data[ind-sind] = inputC8->data->data[ind];
	      }

	      inputC8->data->length = eind - sind;

	    }

	    break;
	    
	  case LAL_Z_TYPE_CODE:

	    /***********************************************************/
	    /* Note: might not be valid; see LAL_C_TYPE_CODE */
	    /***********************************************************/

	    inputC16 = (COMPLEX16TimeSeries *)calloc(1,sizeof(COMPLEX16TimeSeries));
	    inputC16->data = NULL;
	    LAL_CALL( LALFrGetCOMPLEX16TimeSeries( &stat, inputC16, &ChannelIn, stream),&stat);
	    inputC16->epoch.gpsSeconds     = (INT4)floor(atof(T0));
	    inputC16->epoch.gpsNanoSeconds = (INT4)(1e9*(atof(T0) - floor(atof(T0))));
	    LAL_CALL( LALFrSeek(&stat, &(inputC16->epoch), stream),&stat);
	    /* allocate memory */
	    {
	      double t0 = (double)(inputC16->epoch.gpsSeconds) + 1e-9*(double)(inputC16->epoch.gpsNanoSeconds);
	      INT4 length = (INT4)floor((atof(T1) - t0)/inputC16->deltaT);

	      if(inputC16->data) {
		LAL_CALL(LALZDestroyVector(&stat,&(inputC16->data)),&stat);
		inputC16->data = NULL;
	      }

	      LAL_CALL ( LALZCreateVector( &stat, &(inputC16->data), length), &stat);
	    }
	    
	    /* get the data */
	    LAL_CALL( LALFrGetCOMPLEX16TimeSeries( &stat, inputC16, &ChannelIn, stream),&stat);
	    
	    /* datacond convention: f0 = 1/sampling rate */
	    inputC16->f0 = 1.0/inputC16->deltaT;
	  
	    /* Calculate end */
	    {
	      int iend;
	      iend = (int)floor((atof(T1) - (double)(inputC16->epoch.gpsSeconds) + 1e-9*(double)(inputC16->epoch.gpsNanoSeconds))/inputC16->deltaT);
	      iend = (iend>0) ? iend : 1;
	      inputC16->data->length = iend;
	    }

	    /* truncate and shift if necessary */
	    if(end>0) {
	      int ind, sind, eind;
	      sind = (int)floor(start / inputC16->deltaT);
	      eind = (int)floor((start+end) / inputC16->deltaT);
	      eind = (eind > inputC16->data->length) ? inputC16->data->length : eind;

	      for(ind=sind; ind<eind; ind++) {
		inputC16->data->data[ind-sind] = inputC16->data->data[ind];
	      }

	      inputC16->data->length = eind - sind;

	    }

	    break;
	  
	  default:
	    fprintf(stderr,"unsupported frame type\n");
	    return 1;
	  };
	  break;

	  /*************************************************************/

	case FrequencySeries:

	  switch ( otype ) {

	  case LAL_S_TYPE_CODE:

	    fprintf(stderr,"unsupported frame type\n");
	    return 1;

	    break;
	    
	  case LAL_D_TYPE_CODE:

	    fprintf(stderr,"unsupported frame type\n");
	    return 1;

	    break;
	    
	  case LAL_C_TYPE_CODE:
	    inputFC8 = (COMPLEX8FrequencySeries *)calloc(1,sizeof(COMPLEX8FrequencySeries));
	    inputFC8->data = NULL;
	    LAL_CALL( LALFrGetCOMPLEX8FrequencySeries( &stat, inputFC8, &ChannelIn, stream),&stat);

	    { 
	      /* end GPS time needed by translation function */

	      double t = atof(T1);
	      endEpoch.gpsSeconds = (INT4)floor(t);
	      endEpoch.gpsNanoSeconds = (INT4)(t - floor(t));
	    }

	    /* truncate and shift if necessary */
	    if(end>0) {
	      int ind, sind, eind;
	      sind = (int)floor(start / inputFC8->deltaF);
	      eind = (int)floor((start+end) / inputFC8->deltaF);
	      eind = (eind > inputFC8->data->length) ? inputFC8->data->length : eind;

	      for(ind=sind; ind<eind; ind++) {
		inputFC8->data->data[ind-sind] = inputFC8->data->data[ind];
	      }

	      inputFC8->data->length = eind - sind;

	    }

	    break;
	    
	  case LAL_Z_TYPE_CODE:

	    fprintf(stderr,"unsupported frame type\n");
	    return 1;

	    break;
	  
	  default:
	    fprintf(stderr,"unsupported frame type\n");
	    return 1;
	  };
	  break;

	default:
	  fprintf(stderr,"unsupported frame type\n");
	  return 1;
	};

	/* close the frame stream */
	LAL_CALL( LALFrClose( &stat, &stream ),&stat);

      }

      /* add to datacond symbols */
      *symbols = (datacond_symbol_type *)realloc(*symbols, (1 + *Nsymbols) * sizeof(datacond_symbol_type));
      (*symbols + *Nsymbols)->s_direction = DATACOND_SYMBOL_INPUT; /* input */
      (*symbols + *Nsymbols)->s_symbol_name = alias;
      (*symbols + *Nsymbols)->s_aux_data = NULL;
      
      switch ( stype ) {
      case TimeSeries:

	switch ( otype ) {
	case LAL_S_TYPE_CODE:
	  (*symbols + *Nsymbols)->s_translator = TranslateREAL4TimeSeries;
	  (*symbols + *Nsymbols)->s_user_data = inputR4;
	  break;
	case LAL_D_TYPE_CODE:
	  (*symbols + *Nsymbols)->s_translator = TranslateREAL8TimeSeries;
	  (*symbols + *Nsymbols)->s_user_data = inputR8;
	  break;
	case LAL_C_TYPE_CODE:
	  (*symbols + *Nsymbols)->s_translator = TranslateCOMPLEX8TimeSeries;
	  (*symbols + *Nsymbols)->s_user_data = inputC8;
	  break;
	case LAL_Z_TYPE_CODE:
	  (*symbols + *Nsymbols)->s_translator = TranslateCOMPLEX16TimeSeries;
	  (*symbols + *Nsymbols)->s_user_data = inputC16;
	  break;
	default:
	  fprintf(stderr,"unsupported frame type\n");
	  return 1;
	};
	break;

      case FrequencySeries:

	switch ( otype ) {
	case LAL_S_TYPE_CODE:
	  /* note: invalid */
	  (*symbols + *Nsymbols)->s_translator = TranslateREAL4FrequencySeries;
	  (*symbols + *Nsymbols)->s_user_data = inputFR4;
	  break;
	case LAL_D_TYPE_CODE:
	  /* note: invalid */
	  (*symbols + *Nsymbols)->s_translator = TranslateREAL8FrequencySeries;
	  (*symbols + *Nsymbols)->s_user_data = inputFR8;
	  break;
	case LAL_C_TYPE_CODE:
	  (*symbols + *Nsymbols)->s_translator = TranslateCOMPLEX8FrequencySeries;
	  (*symbols + *Nsymbols)->s_user_data = inputFC8;
	  (*symbols + *Nsymbols)->s_aux_data = &endEpoch; 
	  break;
	case LAL_Z_TYPE_CODE:
	  /* note: invalid */
	  (*symbols + *Nsymbols)->s_translator = TranslateCOMPLEX16FrequencySeries;
	  (*symbols + *Nsymbols)->s_user_data = inputFC16;
	  break;
	default:
	  fprintf(stderr,"unsupported frame type\n");
	  return 1;
	};
	break;

      default:
	fprintf(stderr,"unsupported frame type\n");
	return 1;
      };
	
      (*Nsymbols)++;
    }

    /* update line being processed */
    if(p1) {
      p0 = p1+1;
      p1 = strchr(p0,'\n');
      bzero(buf,strlen(buf));
      if(p1) {
	memcpy(buf,p0,p1-p0+1);
      } else {
	strcpy(buf,p0);
      }
    } else {
      break;
    }

  }

  free(buf);

  return 0;
}


/* code stolen from FrameSeries.c */
static FrVect * FrVectReadInfo( FrFile *iFile, FRULONG *pos )
{
  FrVect *v;
  unsigned short type;
  FRULONG localpos;
  if ( FrFileIGoToNextRecord( iFile ) != iFile->vectorType )
    return NULL;
  v = calloc( 1, sizeof( FrVect ) );
  if ( ! v )
  {
    iFile->error = FR_ERROR_MALLOC_FAILED;
    return NULL;
  }
  FrReadHeader( iFile, v );
  FrReadSChar( iFile, &v->name );
  FrReadShortU( iFile, &v->compress );
  if ( v->compress == 256 )
    v->compress = 0; /* we will swap bytes at reading time */
  FrReadShortU( iFile, &type );
  v->type = type;
  switch ( v->type )
  {
    case FR_VECT_4R:
      v->wSize = sizeof( float );
      break;
    case FR_VECT_8R:
      v->wSize = sizeof( double );
      break;
    case FR_VECT_C:
      v->wSize = sizeof( char );
      break;
    case FR_VECT_1U:
      v->wSize = sizeof( char );
      break;
    case FR_VECT_2S:
      v->wSize = sizeof( short );
      break;
    case FR_VECT_2U:
      v->wSize = sizeof( short );
      break;
    case FR_VECT_4S:
      v->wSize = sizeof( int );
      break;
    case FR_VECT_4U:
      v->wSize = sizeof( int );
      break;
    case FR_VECT_8S:
      v->wSize = sizeof( FRLONG );
      break;
    case FR_VECT_8U:
      v->wSize = sizeof( FRLONG );
      break;
    case FR_VECT_8C:
      v->wSize = 2 * sizeof( float );
      break;
    case FR_VECT_16C:
      v->wSize = 2 * sizeof( double );
      break;
    case FR_VECT_8H:
      v->wSize = 2 * sizeof( float );
      break;
    case FR_VECT_16H:
      v->wSize = 2 * sizeof( double );
      break;
    default:
      v->wSize = 0;
  }
  if ( iFile->fmtVersion > 5 )
  {
    FrReadLong( iFile, (FRLONG *)&v->nData );
    FrReadLong( iFile, (FRLONG *)&v->nBytes );
  }
  else
  {
    unsigned int nData, nBytes;
    FrReadIntU( iFile, &nData );
    FrReadIntU( iFile, &nBytes );
    v->nData  = nData;
    v->nBytes = nBytes;
  }
  v->space = v->nData;

  /* skip the data */
  localpos = FrIOTell( iFile->frfd );
  if ( pos )
    *pos = localpos;
  localpos += v->nBytes;
  FrIOSet( iFile->frfd, localpos );
  FrReadIntU( iFile, &v->nDim );
  FrReadVL( iFile, (FRLONG**)&v->nx, v->nDim );
  FrReadVD( iFile, &v->dx, v->nDim );
  FrReadVD( iFile, &v->startX, v->nDim );
  FrReadVQ( iFile, &v->unitX, v->nDim );
  FrReadSChar( iFile, &v->unitY );
  FrReadStruct( iFile, &v->next );
 
  if ( pos )
    *pos = localpos;
  return v;
}


void
LALFrGetSeriesType(
    LALStatus   *status,
    LALTYPECODE *output,
    ObjectType *objtype,
    FrChanIn    *chanin,
    FrStream    *stream
    )
{ 
  FrChanType chantype;
  INT4 type = -1;
  INT4 otype = 1; /* assume time series */
  FrTOCts    *ts   = NULL;
  FrProcData *proc = NULL;
  FrAdcData  *adc  = NULL;
  FrSimData  *sim  = NULL;
  INITSTATUS( status, "FUNC", GETFRAMEDATAC );  

  ASSERT( output, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( chanin, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( chanin->name, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );

  if ( stream->state & LAL_FR_ERR )
  {
    ABORT( status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR );
  }
  if ( stream->state & LAL_FR_END )
  {
    ABORT( status, FRAMESTREAMH_EDONE, FRAMESTREAMH_MSGEDONE );
  }

  if ( ! stream->file->toc )
  {
    if ( FrTOCReadFull( stream->file ) == NULL )
    {
      LALPrintError( "Error opening frame TOC %s\n", stream->file );
      stream->state |= LAL_FR_ERR | LAL_FR_TOC;
      ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
    }
  }

  /* scan adc channels */
  chantype = LAL_ADC_CHAN;
  ts = stream->file->toc->adc;
  while ( ts && strcmp( chanin->name, ts->name ) )
    ts = ts->next;

  if ( ! ts )
  {
    /* scan sim channels */
    chantype = LAL_SIM_CHAN;
    ts = stream->file->toc->sim;
    while ( ts && strcmp( chanin->name, ts->name ) )
      ts = ts->next;
  }

  if ( ! ts )
  {
    /* scan proc channels */
    chantype = LAL_PROC_CHAN;
    ts = stream->file->toc->proc;
    while ( ts && strcmp( chanin->name, ts->name ) )
      ts = ts->next;
  }

  if ( ts ) /* the channel was found */
  {
    FrTOCSetPos( stream->file, ts->position[0] );
    switch ( chantype )
    {
      case LAL_ADC_CHAN:
        if ( ( adc = FrAdcDataRead( stream->file ) ) )
        {
          adc->data = FrVectReadInfo( stream->file, NULL );
          type = adc->data ? adc->data->type : -1;
          FrAdcDataFree( adc );
        }
        break;
      case LAL_SIM_CHAN:
        if ( ( sim = FrSimDataRead( stream->file ) ) )
        {
          sim->data = FrVectReadInfo( stream->file, NULL );
          type = sim->data ? sim->data->type : -1;
          FrSimDataFree( sim );
        }
        break;
      case LAL_PROC_CHAN:
        if ( ( proc = FrProcDataRead( stream->file ) ) )
        {
          proc->data = FrVectReadInfo( stream->file, NULL );
	  otype = proc->data ? proc->type : -1;
          type = proc->data ? proc->data->type : -1;
          FrProcDataFree( proc );
        }
        break;
      default:
        type = -1;
    }
  } else {
    ABORT ( status, FRAMESTREAMH_ETYPE, "Channel not found" );
  }

  switch ( type )
  {
    case FR_VECT_C:
      *output = LAL_CHAR_TYPE_CODE;
      break;
    case FR_VECT_2S:
      *output = LAL_I2_TYPE_CODE;
      break;
    case FR_VECT_4S:
      *output = LAL_I4_TYPE_CODE;
      break;
    case FR_VECT_8S:
      *output = LAL_I8_TYPE_CODE;
      break;
    case FR_VECT_1U:
      *output = LAL_UCHAR_TYPE_CODE;
      break;
    case FR_VECT_2U:
      *output = LAL_U2_TYPE_CODE;
      break;
    case FR_VECT_4U:
      *output = LAL_U4_TYPE_CODE;
      break;
    case FR_VECT_8U:
      *output = LAL_U8_TYPE_CODE;
      break;
    case FR_VECT_4R:
      *output = LAL_S_TYPE_CODE;
      break;
    case FR_VECT_8R:
      *output = LAL_D_TYPE_CODE;
      break;
    case FR_VECT_C8:
      *output = LAL_C_TYPE_CODE;
      break;
    case FR_VECT_C16:
      *output = LAL_Z_TYPE_CODE;
      break;
    default:
      ABORT( status, FRAMESTREAMH_ETYPE, FRAMESTREAMH_MSGETYPE );      
  }

  switch ( otype ) {
  case 0:
    *objtype = Unknown;
    break;
  case 1:
    *objtype = TimeSeries;
    break;
  case 2:
    *objtype = FrequencySeries;
    break;
  case 3:
    *objtype = Other1D;
    break;
  case 4:
    *objtype = TimeFrequency;
    break;
  case 5:
    *objtype = Wavelets;
    break;
  case 6:
    *objtype = MultiDimensional;
    break;
  default:
    ABORT( status, FRAMESTREAMH_ETYPE, FRAMESTREAMH_MSGETYPE );      
  }

  RETURN( status );
}



