/*********************************************************************************/
/*                         Driver for calibration function                       */
/*                                                                               */
/*                                  X. Siemens                                   */
/*                                                                               */
/*                               UWM - August 2004                               */
/*********************************************************************************/

#include <config.h>
#if !defined HAVE_LIBGSL || !defined HAVE_LIBLALFRAME
#include <stdio.h>
int main(void) {fputs("disabled, no gsl or no lal frame library support.\n", stderr);return 1;}
#else

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <glob.h>
#include <errno.h>

#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/AVFactories.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <lal/Window.h>
#include <lal/Calibration.h>
#include <lal/LALConstants.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/ConfigFile.h>

extern char *optarg;
extern int optind, opterr, optopt;

#define MAXFORDER 100                     /* Maximum filter order */

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 100; } else ((void)0)

/***************************************************************************/

/* STRUCTURES */
struct CommandLineArgsTag {
  REAL8 f;                 /* Frequency of the calibration line */
  REAL8 To;                /* factors integration time */
  REAL8 G0Re;              /* Real part of open loop gain at cal line freq.*/
  REAL8 G0Im;              /* Imaginary part of open loop gain at cal line freq. */
  REAL8 D0Re;              /* Real part of digital filter at cal line freq.*/
  REAL8 D0Im;              /* Imaginary part of digital filter at cal line freq. */
  INT4 GPSStart;           /* Start GPS time for the segment to be calibrated */
  INT4 GPSEnd;             /* End GPS time for the segment to be calibrated */
  char *FrCacheFile;       /* Frame cache file for corresponding time */
  char *exc_chan;          /* excitation channel name */    
  char *darm_chan;         /* darm channel name */ 
  char *asq_chan;          /* asq channel name */
  char *filterfile;        /* file with filter coefficients */
} CommandLineArgs;


/***************************************************************************/

/* GLOBAL VARIABLES */

StrainIn InputData;
StrainOut OutputData;
INT4 duration;
LIGOTimeGPS gpsStartepoch;

static LALStatus status;
INT4 lalDebugLevel=3;
FrCache *framecache;                                           /* frame reading variables */
FrStream *framestream=NULL;

/***************************************************************************/

/* FUNCTION PROTOTYPES */
/* Reads the command line */
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);

/* Reads T seconds of AS_Q, EXC, and DARM_CTRL data */
int ReadData(struct CommandLineArgsTag CLA);

/* Reads the filters file */
int ReadFiltersFile(struct CommandLineArgsTag CLA);

/* Writes the filters file */
int WriteFiltersFile(struct CommandLineArgsTag CLA);

/* Frees the memory */
int FreeMem(void);                                        


/************************************* MAIN PROGRAM *************************************/

int main(int argc,char *argv[])
{
  int p;

  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;
  if (ReadData(CommandLineArgs)) return 2;
  if (ReadFiltersFile(CommandLineArgs)) return 3;
  
/*   LALComputeStrain(&status, &OutputData, &InputData); */
/*   TESTSTATUS( &status ); */

  if (WriteFiltersFile(CommandLineArgs)) return 4;

/*   for (p=0; p < (INT4) OutputData.h.data->length; p++) { */
/*      fprintf(stdout,"%e %1.17e\n",p*OutputData.h.deltaT,OutputData.h.data->data[p]);  */
/*   } */

  if(FreeMem()) return 8;

  return 0;
}

/************************************* MAIN PROGRAM ENDS *************************************/

/*  FUNCTIONS */

/*******************************************************************************/

int ReadData(struct CommandLineArgsTag CLA)
{

FrPos pos1;

static FrChanIn chanin_darm;
static FrChanIn chanin_asq;
static FrChanIn chanin_exc;


  chanin_asq.type  = ADCDataChannel;
  chanin_darm.type = ADCDataChannel;
  chanin_exc.type  = ADCDataChannel;

  chanin_asq.name  = CLA.asq_chan;
  chanin_darm.name = CLA.darm_chan;
  chanin_exc.name  = CLA.exc_chan; 

  /* create Frame cache, open frame stream and delete frame cache */
  LALFrCacheImport(&status,&framecache,CommandLineArgs.FrCacheFile);
  TESTSTATUS( &status );
  LALFrCacheOpen(&status,&framestream,framecache);
  TESTSTATUS( &status );
  LALDestroyFrCache(&status,&framecache);
  TESTSTATUS( &status );

  /* Get channel time step size by calling LALFrGetREAL4TimeSeries */
  LALFrSeek(&status,&gpsStartepoch,framestream);
  TESTSTATUS( &status );
  LALFrGetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.AS_Q,&chanin_asq,framestream);
  TESTSTATUS( &status );
  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.DARM,&chanin_darm,framestream);
  TESTSTATUS( &status );
  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.EXC,&chanin_exc,framestream);
  TESTSTATUS( &status );

  /* Allocate space for data vectors */
  LALSCreateVector(&status,&InputData.AS_Q.data,(UINT4)(duration/InputData.AS_Q.deltaT +0.5));
  TESTSTATUS( &status );
  LALSCreateVector(&status,&InputData.DARM.data,(UINT4)(duration/InputData.DARM.deltaT +0.5));
  TESTSTATUS( &status );
  LALSCreateVector(&status,&InputData.EXC.data,(UINT4)(duration/InputData.EXC.deltaT +0.5));
  TESTSTATUS( &status );

  /* Read in the data */
  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.AS_Q,&chanin_asq,framestream);
  TESTSTATUS( &status );
  
  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.DARM,&chanin_darm,framestream);
  TESTSTATUS( &status );

  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.EXC,&chanin_exc,framestream);
  TESTSTATUS( &status );

  LALFrClose(&status,&framestream);
  TESTSTATUS( &status );

  /* Set the rest of the input variables */
  InputData.Go.re=CLA.G0Re;
  InputData.Go.im=CLA.G0Im;
  InputData.Do.re=CLA.D0Re;
  InputData.Do.im=CLA.D0Im;
  InputData.f=CLA.f;
  InputData.To=CLA.To;

  /* Allocate output data */
  OutputData.h.epoch=InputData.AS_Q.epoch;
  OutputData.h.deltaT=InputData.AS_Q.deltaT;
  LALDCreateVector(&status,&OutputData.h.data,(UINT4)(duration/OutputData.h.deltaT +0.5));
  TESTSTATUS( &status );

  OutputData.alpha.deltaT=CLA.To;
  LALZCreateVector(&status,&OutputData.alpha.data,(UINT4)(duration/OutputData.alpha.deltaT +0.5));
  TESTSTATUS( &status );

  OutputData.beta.deltaT=CLA.To;
  LALZCreateVector(&status,&OutputData.beta.data,(UINT4)(duration/OutputData.beta.deltaT +0.5));
  TESTSTATUS( &status );

  return 0;
}



/*******************************************************************************/

int ReadFiltersFile(struct CommandLineArgsTag CLA)
{
  LALParsedDataFile *Filters =NULL;	/* pre-parsed contents of Filters-file */
  int numlines,fileGPS,i,j,n;
  CHAR *thisline;
  char gpsstr[3],sensingstr[7],usfstr[17], delaystr[5],servostr[5];
  char aastr[16],axstr[11],aystr[11];
  int Corders[MAXFORDER], Gorders[MAXFORDER], AAorders[MAXFORDER], AXorders[MAXFORDER], AYorders[MAXFORDER];

  LALParseDataFile (&status, &Filters, CLA.filterfile);
  TESTSTATUS( &status );

  numlines = Filters->lines->nTokens; /* how many lines of data */

  /* check that file is not empty */
  if (numlines == 0)
    {
      fprintf(stderr,"File %s has no contents!\n",CLA.filterfile);
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;
    }
  
  i=0; /* start at line 0 */
  /* get GPS time of file to later agrees with starting GPS */
  thisline = Filters->lines->tokens[i];	/* get line i */
  sscanf (thisline,"%" LAL_INT4_FORMAT " %s", &fileGPS, gpsstr);
  if ( strcmp(gpsstr, "GPS" ) ) 
    {
      fprintf(stderr,"ERROR: First line (%s) of file %s is not properly terminated by '%s' marker!\n\n", thisline, CLA.filterfile, "GPS");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;
    }

  /**-----------------------------------------------------------------------**/
  /* read sensing function info */
  i++; /*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  sscanf (thisline,"%s", sensingstr);
  if ( strcmp(sensingstr, "SENSING" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", thisline, CLA.filterfile, "SENSING");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;
    }

  i++; /*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  sscanf (thisline,"%" LAL_INT4_FORMAT " %s", &InputData.CinvUSF, usfstr);
  if ( strcmp(usfstr, "UPSAMPLING_FACTOR" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", thisline, CLA.filterfile, "UPSAMPLING_FACTOR");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;
    }
  /* FIXME: Check upsamplig factor USF, positive, and mod 2=0 */

  /* Read Delay */
  i++; /*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  sscanf (thisline,"%" LAL_INT4_FORMAT " %s", &InputData.CinvDelay, delaystr);
  if ( strcmp(delaystr, "DELAY" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", thisline, CLA.filterfile, "DELAY");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;
    }

  /* Read number of sensing filters and their orders */
  i++;/*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  InputData.NCinv=strtol(thisline, &thisline,10);
  /* FIXME: Check that InputData.NCinv < 100 */
  for(j=0; j < InputData.NCinv; j++)
    {
      Corders[j]=strtol(thisline, &thisline,10);
    }  
  if ( strcmp(thisline, " FILTERS_ORDERS" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", thisline, CLA.filterfile, "FILTERS_ORDERS");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;      
    }
   
  /* Allocate inverse sensing funtion filters */
  InputData.Cinv=(REAL8IIRFilter *)LALMalloc(sizeof(REAL8IIRFilter) * InputData.NCinv); 

  /* Allocate inverse sensing function filter */
  for(n = 0; n < InputData.NCinv; n++)
    {
      int l;
      InputData.Cinv[n].directCoef=NULL;
      InputData.Cinv[n].recursCoef=NULL;
      InputData.Cinv[n].history=NULL;

      LALDCreateVector(&status,&(InputData.Cinv[n].directCoef),Corders[n]);
      LALDCreateVector(&status,&(InputData.Cinv[n].recursCoef),Corders[n]);
      LALDCreateVector(&status,&(InputData.Cinv[n].history),Corders[n]-1);

      for(l=0;l<Corders[n];l++) InputData.Cinv[n].directCoef->data[l]=0.0;
      for(l=0;l<Corders[n];l++) InputData.Cinv[n].recursCoef->data[l]=0.0;
      for(l=0;l<Corders[n]-1;l++) InputData.Cinv[n].history->data[l]=0.0;
    }
  for(n = 0; n < InputData.NCinv; n++)
    {
      /* read direct coeffs */
      i++;/*advance one line */
      thisline = Filters->lines->tokens[i];	/* get line i */
      
      for (j=0; j < (int) InputData.Cinv[n].directCoef->length; j++)
	{
	  InputData.Cinv[n].directCoef->data[j]=strtod(thisline, &thisline);
	}      
      /* read recursive coeffs */
      i++;/*advance one line */
      thisline = Filters->lines->tokens[i];	/* get line i */
      for (j=0; j < (int) InputData.Cinv[n].recursCoef->length; j++)
	{
	  InputData.Cinv[n].recursCoef->data[j]=-strtod(thisline, &thisline);
	}
      /* read histories; if this is the correct file */
      i++;/*advance one line */
      thisline = Filters->lines->tokens[i];	/* get line i */
      if (fileGPS == OutputData.h.epoch.gpsSeconds)
	{
	  
	  for (j=0; j < (int) InputData.Cinv[n].history->length; j++)
	    {
	      InputData.Cinv[n].history->data[j]=strtod(thisline, &thisline);
	    }
	}
    }

  /**-----------------------------------------------------------------------**/
  /* Read in servo filters */
  i++; /*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  sscanf (thisline,"%s", servostr);
  if ( strcmp(servostr, "SERVO" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", thisline, CLA.filterfile, "SERVO");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;
    }
  /* Read number of sensing filters and their orders */
  i++;/*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  InputData.NG=strtol(thisline, &thisline,10);
  /* FIXME: Check that InputData.NCinv < 100 */
  for(j=0; j < InputData.NG; j++)
    {
      Gorders[j]=strtol(thisline, &thisline,10);
    }  
  if ( strcmp(thisline, " FILTERS_ORDERS" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", thisline, CLA.filterfile, "FILTERS_ORDERS");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;      
    }
   
  /* Allocate inverse sensing funtion filters */
  InputData.G=(REAL8IIRFilter *)LALMalloc(sizeof(REAL8IIRFilter) * InputData.NG); 

  /* Allocate inverse sensing function filter */
  for(n = 0; n < InputData.NG; n++)
    {
      int l;
      InputData.G[n].directCoef=NULL;
      InputData.G[n].recursCoef=NULL;
      InputData.G[n].history=NULL;

      LALDCreateVector(&status,&(InputData.G[n].directCoef),Gorders[n]);
      LALDCreateVector(&status,&(InputData.G[n].recursCoef),Gorders[n]);
      LALDCreateVector(&status,&(InputData.G[n].history),Gorders[n]-1);

      for(l=0;l<Gorders[n];l++) InputData.G[n].directCoef->data[l]=0.0;
      for(l=0;l<Gorders[n];l++) InputData.G[n].recursCoef->data[l]=0.0;
      for(l=0;l<Gorders[n]-1;l++) InputData.G[n].history->data[l]=0.0;
    }
  for(n = 0; n < InputData.NG; n++)
    {
      /* read direct coeffs */
      i++;/*advance one line */
      thisline = Filters->lines->tokens[i];	/* get line i */
      
      for (j=0; j < (int) InputData.G[n].directCoef->length; j++)
	{
	  InputData.G[n].directCoef->data[j]=strtod(thisline, &thisline);
	}      
      /* read recursive coeffs */
      i++;/*advance one line */
      thisline = Filters->lines->tokens[i];	/* get line i */
      for (j=0; j < (int) InputData.G[n].recursCoef->length; j++)
	{
	  InputData.G[n].recursCoef->data[j]=-strtod(thisline, &thisline);
	}
      /* read histories; if this is the correct file */
      i++;/*advance one line */
      thisline = Filters->lines->tokens[i];	/* get line i */
      if (fileGPS == OutputData.h.epoch.gpsSeconds)
	{
	  
	  for (j=0; j < (int) InputData.G[n].history->length; j++)
	    {
	      InputData.G[n].history->data[j]=strtod(thisline, &thisline);
	    }
	}
    }

  /**-----------------------------------------------------------------------**/
  /* Read in analog actuation filters */
  i++; /*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  sscanf (thisline,"%s", aastr);
  if ( strcmp(aastr, "ACTUATION_ANALOG" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", thisline, CLA.filterfile, "ACTUATION_ANALOG");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;
    }
  /* Read Delay */
  i++; /*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  sscanf (thisline,"%" LAL_INT4_FORMAT " %s", &InputData.AADelay, delaystr);
  if ( strcmp(delaystr, "DELAY" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", thisline, CLA.filterfile, "DELAY");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;
    }

  /* Read number of sensing filters and their orders */
  i++;/*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  InputData.NAA=strtol(thisline, &thisline,10);
  /* FIXME: Check that InputData.NCinv < 100 */
  for(j=0; j < InputData.NAA; j++)
    {
      AAorders[j]=strtol(thisline, &thisline,10);
    }  
  if ( strcmp(thisline, " FILTERS_ORDERS" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", thisline, CLA.filterfile, "FILTERS_ORDERS");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;      
    }
   
  /* Allocate inverse sensing funtion filters */
  InputData.AA=(REAL8IIRFilter *)LALMalloc(sizeof(REAL8IIRFilter) * InputData.NAA); 

  /* Allocate inverse sensing function filter */
  for(n = 0; n < InputData.NAA; n++)
    {
      int l;
      InputData.AA[n].directCoef=NULL;
      InputData.AA[n].recursCoef=NULL;
      InputData.AA[n].history=NULL;

      LALDCreateVector(&status,&(InputData.AA[n].directCoef),AAorders[n]);
      LALDCreateVector(&status,&(InputData.AA[n].recursCoef),AAorders[n]);
      LALDCreateVector(&status,&(InputData.AA[n].history),AAorders[n]-1);

      for(l=0;l<AAorders[n];l++) InputData.AA[n].directCoef->data[l]=0.0;
      for(l=0;l<AAorders[n];l++) InputData.AA[n].recursCoef->data[l]=0.0;
      for(l=0;l<AAorders[n]-1;l++) InputData.AA[n].history->data[l]=0.0;
    }
  for(n = 0; n < InputData.NAA; n++)
    {
      /* read direct coeffs */
      i++;/*advance one line */
      thisline = Filters->lines->tokens[i];	/* get line i */
      
      for (j=0; j < (int) InputData.AA[n].directCoef->length; j++)
	{
	  InputData.AA[n].directCoef->data[j]=strtod(thisline, &thisline);
	}      
      /* read recursive coeffs */
      i++;/*advance one line */
      thisline = Filters->lines->tokens[i];	/* get line i */
      for (j=0; j < (int) InputData.AA[n].recursCoef->length; j++)
	{
	  InputData.AA[n].recursCoef->data[j]=-strtod(thisline, &thisline);
	}
      /* read histories; if this is the correct file */
      i++;/*advance one line */
      thisline = Filters->lines->tokens[i];	/* get line i */
      if (fileGPS == OutputData.h.epoch.gpsSeconds)
	{
	  
	  for (j=0; j < (int) InputData.AA[n].history->length; j++)
	    {
	      InputData.AA[n].history->data[j]=strtod(thisline, &thisline);
	    }
	}
    }

  /**-----------------------------------------------------------------------**/
  /* Read in x-arm actuation filters */
  i++; /*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  sscanf (thisline,"%s", axstr);
  if ( strcmp(axstr, "ACTUATION_X" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", thisline, CLA.filterfile, "ACTUATION_X");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;
    }

  /* Read number of sensing filters and their orders */
  i++;/*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  InputData.NAX=strtol(thisline, &thisline,10);
  /* FIXME: Check that InputData.NCinv < 100 */
  for(j=0; j < InputData.NAX; j++)
    {
      AXorders[j]=strtol(thisline, &thisline,10);
    }  
  if ( strcmp(thisline, " FILTERS_ORDERS" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", thisline, CLA.filterfile, "FILTERS_ORDERS");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;      
    }
   
  /* Allocate inverse sensing funtion filters */
  InputData.AX=(REAL8IIRFilter *)LALMalloc(sizeof(REAL8IIRFilter) * InputData.NAX); 

  /* Allocate inverse sensing function filter */
  for(n = 0; n < InputData.NAX; n++)
    {
      int l;
      InputData.AX[n].directCoef=NULL;
      InputData.AX[n].recursCoef=NULL;
      InputData.AX[n].history=NULL;

      LALDCreateVector(&status,&(InputData.AX[n].directCoef),AXorders[n]);
      LALDCreateVector(&status,&(InputData.AX[n].recursCoef),AXorders[n]);
      LALDCreateVector(&status,&(InputData.AX[n].history),AXorders[n]-1);

      for(l=0;l<AXorders[n];l++) InputData.AX[n].directCoef->data[l]=0.0;
      for(l=0;l<AXorders[n];l++) InputData.AX[n].recursCoef->data[l]=0.0;
      for(l=0;l<AXorders[n]-1;l++) InputData.AX[n].history->data[l]=0.0;
    }
  for(n = 0; n < InputData.NAX; n++)
    {
      /* read direct coeffs */
      i++;/*advance one line */
      thisline = Filters->lines->tokens[i];	/* get line i */
      
      for (j=0; j < (int) InputData.AX[n].directCoef->length; j++)
	{
	  InputData.AX[n].directCoef->data[j]=strtod(thisline, &thisline);
	}      
      /* read recursive coeffs */
      i++;/*advance one line */
      thisline = Filters->lines->tokens[i];	/* get line i */
      for (j=0; j < (int) InputData.AX[n].recursCoef->length; j++)
	{
	  InputData.AX[n].recursCoef->data[j]=-strtod(thisline, &thisline);
	}
      /* read histories; if this is the correct file */
      i++;/*advance one line */
      thisline = Filters->lines->tokens[i];	/* get line i */
      if (fileGPS == OutputData.h.epoch.gpsSeconds)
	{
	  
	  for (j=0; j < (int) InputData.AX[n].history->length; j++)
	    {
	      InputData.AX[n].history->data[j]=strtod(thisline, &thisline);
	    }
	}
    }


  /**-----------------------------------------------------------------------**/
  /* Read in y-arm actuation filters */
  i++; /*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  sscanf (thisline,"%s", aystr);
  if ( strcmp(aystr, "ACTUATION_Y" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", thisline, CLA.filterfile, "ACTUATION_Y");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;
    }

  /* Read number of sensing filters and their orders */
  i++;/*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  InputData.NAY=strtol(thisline, &thisline,10);
  /* FIXME: Check that InputData.NCinv < 100 */
  for(j=0; j < InputData.NAY; j++)
    {
      AYorders[j]=strtol(thisline, &thisline,10);
    }  
  if ( strcmp(thisline, " FILTERS_ORDERS" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", thisline, CLA.filterfile, "FILTERS_ORDERS");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;      
    }
   
  /* Allocate inverse sensing funtion filters */
  InputData.AY=(REAL8IIRFilter *)LALMalloc(sizeof(REAL8IIRFilter) * InputData.NAY); 

  /* Allocate inverse sensing function filter */
  for(n = 0; n < InputData.NAY; n++)
    {
      int l;
      InputData.AY[n].directCoef=NULL;
      InputData.AY[n].recursCoef=NULL;
      InputData.AY[n].history=NULL;

      LALDCreateVector(&status,&(InputData.AY[n].directCoef),AYorders[n]);
      LALDCreateVector(&status,&(InputData.AY[n].recursCoef),AYorders[n]);
      LALDCreateVector(&status,&(InputData.AY[n].history),AYorders[n]-1);

      for(l=0;l<AYorders[n];l++) InputData.AY[n].directCoef->data[l]=0.0;
      for(l=0;l<AYorders[n];l++) InputData.AY[n].recursCoef->data[l]=0.0;
      for(l=0;l<AYorders[n]-1;l++) InputData.AY[n].history->data[l]=0.0;
    }
  for(n = 0; n < InputData.NAY; n++)
    {
      /* read direct coeffs */
      i++;/*advance one line */
      thisline = Filters->lines->tokens[i];	/* get line i */
      
      for (j=0; j < (int) InputData.AY[n].directCoef->length; j++)
	{
	  InputData.AY[n].directCoef->data[j]=strtod(thisline, &thisline);
	}      
      /* read recursive coeffs */
      i++;/*advance one line */
      thisline = Filters->lines->tokens[i];	/* get line i */
      for (j=0; j < (int) InputData.AY[n].recursCoef->length; j++)
	{
	  InputData.AY[n].recursCoef->data[j]=-strtod(thisline, &thisline);
	}
      /* read histories; if this is the correct file */
      i++;/*advance one line */
      thisline = Filters->lines->tokens[i];	/* get line i */
      if (fileGPS == OutputData.h.epoch.gpsSeconds)
	{
	  
	  for (j=0; j < (int) InputData.AY[n].history->length; j++)
	    {
	      InputData.AY[n].history->data[j]=strtod(thisline, &thisline);
	    }
	}
    }

  /**-----------------------------------------------------------------------**/
  LALDestroyParsedDataFile ( &status, &Filters );
  TESTSTATUS( &status );

  return 0;
}

/*******************************************************************************/

int WriteFiltersFile(struct CommandLineArgsTag CLA)
{
  int j,n;
  FILE *FiltersFile;

  FiltersFile=fopen(CLA.filterfile, "w");

  /* write GPS string */
  fprintf(FiltersFile,"%d GPS\n", (int)(OutputData.h.epoch.gpsSeconds
					+OutputData.h.data->length*OutputData.h.deltaT));
  /* sengsing section */
  fprintf(FiltersFile,"SENSING\n");
  fprintf(FiltersFile,"%d UPSAMPLING_FACTOR\n",InputData.CinvUSF);
  fprintf(FiltersFile,"%d DELAY\n",InputData.CinvDelay);

  fprintf(FiltersFile,"%d ",InputData.NCinv);
  n=0;
  while (n < InputData.NCinv) 
    { 
      fprintf(FiltersFile,"%d ",InputData.Cinv[n].directCoef->length);
      n++;
    }
  fprintf(FiltersFile,"FILTERS_ORDERS\n");
  
  for(n = 0; n < InputData.NCinv; n++)
    {
      /* write direct coeffs */
      for(j=0; j < (int)InputData.Cinv[n].directCoef->length; j++)
	{
	  fprintf(FiltersFile,"%1.16e ", InputData.Cinv[n].directCoef->data[j]);
	}
      fprintf(FiltersFile,"\n");
      
      /* write recursive coeffs */
      for (j=0; j < (int)InputData.Cinv[n].recursCoef->length; j++)
	{
	  fprintf(FiltersFile,"%1.16e ", -InputData.Cinv[n].recursCoef->data[j]);
	}
      fprintf(FiltersFile,"\n");

      /* write histories */
      for (j=0; j < (int)InputData.Cinv[n].history->length; j++)
	{
	  fprintf(FiltersFile,"%1.16e ", InputData.Cinv[n].history->data[j]);
	} 
      fprintf(FiltersFile,"\n");
    }
  
  /* servo section */
  fprintf(FiltersFile,"SERVO\n");
  fprintf(FiltersFile,"%d ",InputData.NG);
  n=0;
  while (n < InputData.NG) 
    { 
      fprintf(FiltersFile,"%d ",InputData.G[n].directCoef->length);
      n++;
    }
  fprintf(FiltersFile,"FILTERS_ORDERS\n");
  
  for(n = 0; n < InputData.NG; n++)
    {
      /* write direct coeffs */
      for(j=0; j < (int)InputData.G[n].directCoef->length; j++)
	{
	  fprintf(FiltersFile,"%1.16e ", InputData.G[n].directCoef->data[j]);
	}
      fprintf(FiltersFile,"\n");
      
      /* write recursive coeffs */
      for (j=0; j < (int)InputData.G[n].recursCoef->length; j++)
	{
	  fprintf(FiltersFile,"%1.16e ", -InputData.G[n].recursCoef->data[j]);
	}
      fprintf(FiltersFile,"\n");

      /* write histories */
      for (j=0; j < (int)InputData.G[n].history->length; j++)
	{
	  fprintf(FiltersFile,"%1.16e ", InputData.G[n].history->data[j]);
	} 
      fprintf(FiltersFile,"\n");
    }

  /* analog actuation section */
  fprintf(FiltersFile,"ACTUATION_ANALOG\n");
  fprintf(FiltersFile,"%d DELAY\n",InputData.AADelay);

  fprintf(FiltersFile,"%d ",InputData.NAA);
  n=0;
  while (n < InputData.NAA) 
    { 
      fprintf(FiltersFile,"%d ",InputData.AA[n].directCoef->length);
      n++;
    }
  fprintf(FiltersFile,"FILTERS_ORDERS\n");
  
  for(n = 0; n < InputData.NAA; n++)
    {
      /* write direct coeffs */
      for(j=0; j < (int)InputData.AA[n].directCoef->length; j++)
	{
	  fprintf(FiltersFile,"%1.16e ", InputData.AA[n].directCoef->data[j]);
	}
      fprintf(FiltersFile,"\n");
      
      /* write recursive coeffs */
      for (j=0; j < (int)InputData.AA[n].recursCoef->length; j++)
	{
	  fprintf(FiltersFile,"%1.16e ", -InputData.AA[n].recursCoef->data[j]);
	}
      fprintf(FiltersFile,"\n");

      /* write histories */
      for (j=0; j < (int)InputData.AA[n].history->length; j++)
	{
	  fprintf(FiltersFile,"%1.16e ", InputData.AA[n].history->data[j]);
	} 
      fprintf(FiltersFile,"\n");
    }
  
  /* x-arm actuation section */
  fprintf(FiltersFile,"ACTUATION_X\n");

  fprintf(FiltersFile,"%d ",InputData.NAX);
  n=0;
  while (n < InputData.NAX) 
    { 
      fprintf(FiltersFile,"%d ",InputData.AX[n].directCoef->length);
      n++;
    }
  fprintf(FiltersFile,"FILTERS_ORDERS\n");
  
  for(n = 0; n < InputData.NAX; n++)
    {
      /* write direct coeffs */
      for(j=0; j < (int)InputData.AX[n].directCoef->length; j++)
	{
	  fprintf(FiltersFile,"%1.16e ", InputData.AX[n].directCoef->data[j]);
	}
      fprintf(FiltersFile,"\n");
      
      /* write recursive coeffs */
      for (j=0; j < (int)InputData.AX[n].recursCoef->length; j++)
	{
	  fprintf(FiltersFile,"%1.16e ", -InputData.AX[n].recursCoef->data[j]);
	}
      fprintf(FiltersFile,"\n");

      /* write histories */
      for (j=0; j < (int)InputData.AX[n].history->length; j++)
	{
	  fprintf(FiltersFile,"%1.16e ", InputData.AX[n].history->data[j]);
	} 
      fprintf(FiltersFile,"\n");
    }
  
  /* y-arm actuation section */
  fprintf(FiltersFile,"ACTUATION_Y\n");

  fprintf(FiltersFile,"%d ",InputData.NAY);
  n=0;
  while (n < InputData.NAY) 
    { 
      fprintf(FiltersFile,"%d ",InputData.AY[n].directCoef->length);
      n++;
    }
  fprintf(FiltersFile,"FILTERS_ORDERS\n");
  
  for(n = 0; n < InputData.NAY; n++)
    {
      /* write direct coeffs */
      for(j=0; j < (int)InputData.AY[n].directCoef->length; j++)
	{
	  fprintf(FiltersFile,"%1.16e ", InputData.AY[n].directCoef->data[j]);
	}
      fprintf(FiltersFile,"\n");
      
      /* write recursive coeffs */
      for (j=0; j < (int)InputData.AY[n].recursCoef->length; j++)
	{
	  fprintf(FiltersFile,"%1.16e ", -InputData.AY[n].recursCoef->data[j]);
	}
      fprintf(FiltersFile,"\n");

      /* write histories */
      for (j=0; j < (int)InputData.AY[n].history->length; j++)
	{
	  fprintf(FiltersFile,"%1.16e ", InputData.AY[n].history->data[j]);
	} 
      fprintf(FiltersFile,"\n");
    }


  fclose(FiltersFile);  
  return 0;
}

/*******************************************************************************/

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA) 
{
  INT4 c, errflg = 0;
  optarg = NULL;
  
  /* Initialize default values */
  CLA->f=0.0;
  CLA->To=0.0;
  CLA->FrCacheFile=NULL;
  CLA->filterfile=NULL;
  CLA->exc_chan=NULL;
  CLA->darm_chan=NULL;
  CLA->asq_chan=NULL;
  CLA->GPSStart=0;
  CLA->GPSEnd=0;
  CLA->G0Re=0.0;
  CLA->G0Im=0.0;
  CLA->D0Re=0.0;
  CLA->D0Im=0.0;

  /* Scan through list of command line arguments */
  while (!errflg && ((c = getopt(argc, argv,"hf:C:A:E:D:s:e:i:j:k:l:t:F:"))!=-1))
    switch (c) {
    case 'f':
      /* calibration line frequency */
      CLA->f=atof(optarg);
      break;
    case 't':
      /* calibration line frequency */
      CLA->To=atof(optarg);
      break;
    case 'C':
      /* name of frame cache file */
      CLA->FrCacheFile=optarg;
      break;
    case 'F':
      /* name of frame cache file */
      CLA->filterfile=optarg;
      break;
    case 'E':
      /* name of excitation channel */
      CLA->exc_chan=optarg;
      break;    
    case 'A':
      /* name of as_q channel */
      CLA->asq_chan=optarg;
      break;    
    case 'D':
      /* name of darm channel */
      CLA->darm_chan=optarg;
      break;    
    case 's':
      /* calibration line frequency */
      CLA->GPSStart=atof(optarg);
      break;
    case 'e':
      /* calibration line frequency */
      CLA->GPSEnd=atof(optarg);
      break;
    case 'i':
      /* calibration line frequency */
      CLA->G0Re=atof(optarg);
      break;
    case 'j':
      /* calibration line frequency */
      CLA->G0Im=atof(optarg);
      break;
    case 'k':
      /* calibration line frequency */
      CLA->D0Re=atof(optarg);
      break;
    case 'l':
      /* calibration line frequency */
      CLA->D0Im=atof(optarg);
      break;
   case 'h':
      /* print usage/help message */
      fprintf(stdout,"All arguments are required. They are:\n");
      fprintf(stdout,"\t-f\tFLOAT\t Calibration line frequency in Hz.\n");
      fprintf(stdout,"\t-t\tFLOAT\t Factors integration time in seconds.\n");
      fprintf(stdout,"\t-i\tFLOAT\t Real part of the open loop gain at the calibration line frequency.\n");
      fprintf(stdout,"\t-j\tFLOAT\t Imaginary part of the open loop gain at the calibration line frequency.\n");
      fprintf(stdout,"\t-k\tFLOAT\t Real part of the digital filter at the calibration line frequency.\n");
      fprintf(stdout,"\t-l\tFLOAT\t Imaginary part of digital filter at the calibration line frequency.\n");
      fprintf(stdout,"\t-s\tINT\t GPS start time.\n");
      fprintf(stdout,"\t-e\tINT\t GPS end time.\n");
      fprintf(stdout,"\t-F\tSTRING\t Name of file containing filters and histories.\n");
      fprintf(stdout,"\t-C\tSTRING\t Name of frame cache file.\n");
      fprintf(stdout,"\t-A\tSTRING\t AS_Q channel name (eg, L1:LSC-AS_Q).\n");
      fprintf(stdout,"\t-E\tSTRING\t Excitation channel name (eg, L1:LSC-ETMX_EXC_DAQ)\n");
      fprintf(stdout,"\t-D\tSTRING\t Darm channel name (eg, L1:LSC-DARM_CTRL)\n");
      exit(0);
      break;
    default:
      /* unrecognized option */
      errflg++;
      fprintf(stderr,"Unrecognized option argument %c\n",c);
      exit(1);
      break;
    }
  if(CLA->f == 0.0)
    {
      fprintf(stderr,"No calibration line frequency specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }      
  if(CLA->To == 0.0)
    {
      fprintf(stderr,"No integration time for the factors specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }      
  if(CLA->G0Re == 0.0 )
    {
      fprintf(stderr,"No real part of open loop gain specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }
  if(CLA->G0Im == 0.0 )
    {
      fprintf(stderr,"No imaginary part of open loop gain specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }
  if(CLA->D0Re == 0.0 )
    {
      fprintf(stderr,"No real part of digital filter specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }
  if(CLA->D0Im == 0.0 )
    {
      fprintf(stderr,"No imaginary part of digital filter specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }
  if(CLA->GPSStart == 0)
    {
      fprintf(stderr,"No GPS start time specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }
  if(CLA->GPSEnd == 0)
    {
      fprintf(stderr,"No GPS end time specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }
  if(CLA->FrCacheFile == NULL)
    {
      fprintf(stderr,"No frame cache file specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }      
  if(CLA->filterfile == NULL)
    {
      fprintf(stderr,"No filter file specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }      
   if(CLA->exc_chan == NULL)
    {
      fprintf(stderr,"No excitation channel specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }      
   if(CLA->darm_chan == NULL)
    {
      fprintf(stderr,"No darm channel specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }      
   if(CLA->asq_chan == NULL)
    {
      fprintf(stderr,"No asq channel specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }      

   /* Set some global variables */
   duration = CLA->GPSEnd - CLA->GPSStart;
   gpsStartepoch.gpsSeconds = CLA->GPSStart;
   gpsStartepoch.gpsNanoSeconds = 0;

  return errflg;
}

/*******************************************************************************/

int FreeMem(void)
{
  int p;

  /* Free input data */
  LALSDestroyVector(&status,&InputData.AS_Q.data);
  TESTSTATUS( &status );
  LALSDestroyVector(&status,&InputData.DARM.data);
  TESTSTATUS( &status );
  LALSDestroyVector(&status,&InputData.EXC.data);
  TESTSTATUS( &status );

  /* Free filters */
  for(p=0;p<InputData.NCinv;p++){
    LALDDestroyVector(&status,&InputData.Cinv[p].directCoef);
    TESTSTATUS( &status );
    LALDDestroyVector(&status,&InputData.Cinv[p].recursCoef);
    TESTSTATUS( &status );
    LALDDestroyVector(&status,&InputData.Cinv[p].history);   
    TESTSTATUS( &status );
  }
  LALFree(InputData.Cinv);

  for(p=0;p<InputData.NG;p++){
    LALDDestroyVector(&status,&InputData.G[p].directCoef);
    TESTSTATUS( &status );
    LALDDestroyVector(&status,&InputData.G[p].recursCoef);
    TESTSTATUS( &status );
    LALDDestroyVector(&status,&InputData.G[p].history);   
    TESTSTATUS( &status );
  }
  LALFree(InputData.G);

  for(p=0;p<InputData.NAA;p++){
    LALDDestroyVector(&status,&InputData.AA[p].directCoef);
    TESTSTATUS( &status );
    LALDDestroyVector(&status,&InputData.AA[p].recursCoef);
    TESTSTATUS( &status );
    LALDDestroyVector(&status,&InputData.AA[p].history);   
    TESTSTATUS( &status );
  }
  LALFree(InputData.AA);

  for(p=0;p<InputData.NAX;p++){
    LALDDestroyVector(&status,&InputData.AX[p].directCoef);
    TESTSTATUS( &status );
    LALDDestroyVector(&status,&InputData.AX[p].recursCoef);
    TESTSTATUS( &status );
    LALDDestroyVector(&status,&InputData.AX[p].history);   
    TESTSTATUS( &status );
  }
  LALFree(InputData.AX);

  for(p=0;p<InputData.NAY;p++){
    LALDDestroyVector(&status,&InputData.AY[p].directCoef);
    TESTSTATUS( &status );
    LALDDestroyVector(&status,&InputData.AY[p].recursCoef);
    TESTSTATUS( &status );
    LALDDestroyVector(&status,&InputData.AY[p].history);   
    TESTSTATUS( &status );
  }
  LALFree(InputData.AY);

  /* Free output data */
  LALDDestroyVector(&status,&OutputData.h.data);
  TESTSTATUS( &status );
  LALZDestroyVector(&status,&OutputData.alpha.data);
  TESTSTATUS( &status );
  LALZDestroyVector(&status,&OutputData.beta.data);
  TESTSTATUS( &status );

  LALCheckMemoryLeaks();
 
  return 0;
}

/*******************************************************************************/
#endif
