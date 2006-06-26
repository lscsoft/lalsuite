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
#include <getopt.h>
#include <stdarg.h>

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
#include <lal/TimeSeries.h>

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
  REAL8 W0Re;              /* Real part of whitening filter at cal line freq.*/
  REAL8 W0Im;              /* Imaginary part of whitening filter at cal line freq. */
  INT4 GPSStart;           /* Start GPS time for the segment to be calibrated */
  INT4 GPSEnd;             /* End GPS time for the segment to be calibrated */
  INT4 testsensing;
  INT4 testactuation;
  char *FrCacheFile;       /* Frame cache file for corresponding time */
  char *exc_chan;          /* excitation channel name */    
  char *darm_chan;         /* darm channel name */ 
  char *darmerr_chan;      /* darm_err  channel name */ 
  char *asq_chan;          /* asq channel name */
  char *filterfile;        /* file with filter coefficients */
  char *frametype;
  char *strainchannel;
  char *datadir;
  char *renamefrom;
  char *renameto;
} CommandLineArgs;


/***************************************************************************/

/* GLOBAL VARIABLES */


StrainIn InputData;
StrainOut OutputData;
INT4 duration;
LIGOTimeGPS gpsStartepoch;

static LALStatus status;
INT4 lalDebugLevel=0;
FrCache *framecache;                                           /* frame reading variables */
FrStream *framestream=NULL;
char ifo[2];
char site[1];

/***************************************************************************/

/* FUNCTION PROTOTYPES */
/* Reads the command line */
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);

/* Reads T seconds of AS_Q, EXC, and DARM_CTRL data */
int ReadData(struct CommandLineArgsTag CLA);

/* Reads the filters file */
int ReadFiltersFile(struct CommandLineArgsTag CLA);

/* Frees the memory */
int FreeMem(void);                                        


/************************************* MAIN PROGRAM *************************************/

int main(int argc,char *argv[])
{
  FrOutPar opar;

  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;
  if (ReadData(CommandLineArgs)) return 2;
  if (ReadFiltersFile(CommandLineArgs)) return 3;
  
  LALComputeStrain(&status, &OutputData, &InputData);
  TESTSTATUS( &status );

  LALResizeREAL8TimeSeries(&status, &(OutputData.hR), (int)(InputData.wings/OutputData.hR.deltaT),
			   OutputData.hR.data->length-2*(UINT4)(InputData.wings/OutputData.hR.deltaT));
  TESTSTATUS( &status );
  LALResizeREAL8TimeSeries(&status, &(OutputData.hC), (int)(InputData.wings/OutputData.hC.deltaT),
			   OutputData.hC.data->length-2*(UINT4)(InputData.wings/OutputData.hC.deltaT));
  TESTSTATUS( &status );
  LALResizeREAL8TimeSeries(&status, &(OutputData.h), (int)(InputData.wings/OutputData.h.deltaT),
			   OutputData.h.data->length-2*(UINT4)(InputData.wings/OutputData.h.deltaT));
  TESTSTATUS( &status );

  strncpy( OutputData.hR.name, CommandLineArgs.strainchannel, sizeof( OutputData.hR.name ) );
  strncpy( OutputData.hC.name, CommandLineArgs.strainchannel, sizeof( OutputData.hC.name ) );
  strncpy( OutputData.h.name,  CommandLineArgs.strainchannel, sizeof( OutputData.h.name  ) );

  opar.source=CommandLineArgs.datadir;
  opar.description=CommandLineArgs.frametype;
  opar.type = ProcDataChannel;
  opar.nframes = 1;
  opar.frame = 0;
  opar.run = 4;

  if (CommandLineArgs.testsensing)
    {
      LALFrWriteREAL8TimeSeries( &status, &OutputData.hR, &opar );
      TESTSTATUS( &status );
    }else if(CommandLineArgs.testactuation){      
      LALFrWriteREAL8TimeSeries( &status, &OutputData.hC, &opar );
      TESTSTATUS( &status );
    }else{
      LALFrWriteREAL8TimeSeries( &status, &OutputData.h, &opar );
      TESTSTATUS( &status );
    }

  if(CommandLineArgs.renamefrom && CommandLineArgs.renameto)
    {
      rename(CommandLineArgs.renamefrom,CommandLineArgs.renameto);
    }
      
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
static FrChanIn chanin_darmerr;
static FrChanIn chanin_asq;
static FrChanIn chanin_exc;

  chanin_asq.type  = ADCDataChannel;
  chanin_darm.type = ADCDataChannel;
  chanin_darmerr.type = ADCDataChannel;
  chanin_exc.type  = ADCDataChannel;

  chanin_asq.name  = CLA.asq_chan;
  chanin_darm.name = CLA.darm_chan;
  chanin_darmerr.name = CLA.darmerr_chan;
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
  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.DARM_ERR,&chanin_darmerr,framestream);
  TESTSTATUS( &status );

  /* Allocate space for data vectors */
  LALSCreateVector(&status,&InputData.AS_Q.data,(UINT4)(duration/InputData.AS_Q.deltaT +0.5));
  TESTSTATUS( &status );
  LALSCreateVector(&status,&InputData.DARM.data,(UINT4)(duration/InputData.DARM.deltaT +0.5));
  TESTSTATUS( &status );
  LALSCreateVector(&status,&InputData.DARM_ERR.data,(UINT4)(duration/InputData.DARM_ERR.deltaT +0.5));
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
  LALFrGetREAL4TimeSeries(&status,&InputData.DARM_ERR,&chanin_darmerr,framestream);
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
  InputData.Wo.re=CLA.W0Re;
  InputData.Wo.im=CLA.W0Im;
  InputData.f=CLA.f;
  InputData.To=CLA.To;

  /* check input data epoch agrees with command line arguments */
  if ( InputData.AS_Q.epoch.gpsSeconds != CLA.GPSStart )
    {
      fprintf(stderr,"GPS start time of data (%d) does not agree with requested start time (%d). Exiting.", 
	      InputData.AS_Q.epoch.gpsSeconds, CLA.GPSStart);
      return 1;
    }

  /* Allocate output data */
  OutputData.h.epoch=InputData.AS_Q.epoch;
  OutputData.h.deltaT=InputData.AS_Q.deltaT;
  LALDCreateVector(&status,&OutputData.h.data,(UINT4)(duration/OutputData.h.deltaT +0.5));
  TESTSTATUS( &status );

  OutputData.hC.epoch=InputData.AS_Q.epoch;
  OutputData.hC.deltaT=InputData.AS_Q.deltaT;
  LALDCreateVector(&status,&OutputData.hC.data,(UINT4)(duration/OutputData.hC.deltaT +0.5));
  TESTSTATUS( &status );

  OutputData.hR.epoch=InputData.AS_Q.epoch;
  OutputData.hR.deltaT=InputData.AS_Q.deltaT;
  LALDCreateVector(&status,&OutputData.hR.data,(UINT4)(duration/OutputData.hR.deltaT +0.5));
  TESTSTATUS( &status );

  OutputData.alpha.deltaT=CLA.To;
  LALZCreateVector(&status,&OutputData.alpha.data,(UINT4)(duration/OutputData.alpha.deltaT +0.5));
  TESTSTATUS( &status );

  OutputData.alphabeta.deltaT=CLA.To;
  LALZCreateVector(&status,&OutputData.alphabeta.data,(UINT4)(duration/OutputData.alphabeta.deltaT +0.5));
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
  int numlines,i,n;
  CHAR *thisline;
  char sensingstr[7],usfstr[17], delaystr[5];
  char aastr[9], servostr[5];
  int NCinv, NA, ND, l;

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
  
  /**-----------------------------------------------------------------------**/
  /* read sensing function info */
  i=0; /*start with first line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  sscanf (thisline,"%s", sensingstr);
  if ( strcmp(sensingstr, "SENSING" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", 
	      thisline, CLA.filterfile, "SENSING");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;
    }

  i++; /*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  sscanf (thisline,"%" LAL_INT4_FORMAT " %s", &InputData.CinvUSF, usfstr);
  if ( strcmp(usfstr, "UPSAMPLING_FACTOR" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", 
	      thisline, CLA.filterfile, "UPSAMPLING_FACTOR");
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
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", 
	      thisline, CLA.filterfile, "DELAY");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;
    }

  /* Read number of sensing filters and their orders */
  i++;/*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  NCinv=strtol(thisline, &thisline,10);
  if ( strcmp(thisline, " FILTER_ORDER" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", 
	      thisline, CLA.filterfile, "FILTERS_ORDERS");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;      
    }
   
  /* Allocate inverse sensing funtion filters */
  InputData.Cinv=(REAL8IIRFilter *)LALMalloc(sizeof(REAL8IIRFilter)); 

  /* Allocate inverse sensing function filter */
  InputData.Cinv->directCoef=NULL;
  InputData.Cinv->recursCoef=NULL;
  InputData.Cinv->history=NULL;

  LALDCreateVector(&status,&(InputData.Cinv->directCoef),NCinv);
  LALDCreateVector(&status,&(InputData.Cinv->recursCoef),NCinv);
  LALDCreateVector(&status,&(InputData.Cinv->history),NCinv-1);

  for(l=0;l<NCinv;l++) InputData.Cinv->directCoef->data[l]=0.0;
  for(l=0;l<NCinv;l++) InputData.Cinv->recursCoef->data[l]=0.0;
  for(l=0;l<NCinv-1;l++) InputData.Cinv->history->data[l]=0.0;
  
  for(n = 0; n < NCinv; n++)
    {
      /* read direct coeffs */
      i++;/*advance one line */
      thisline = Filters->lines->tokens[i];	/* get line i */
      InputData.Cinv->directCoef->data[n]=strtod(thisline, &thisline);  
    }

  /**-----------------------------------------------------------------------**/
  /* Read in servo */
  i++; /*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  sscanf (thisline,"%s", servostr);
  if ( strcmp(servostr, "SERVO" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", 
	      thisline, CLA.filterfile, "SERVO");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;
    }
  /* Read number of sensing filters and their orders */
  i++;/*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  ND=strtol(thisline, &thisline,10);
  if ( strcmp(thisline, " FILTER_ORDER" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", 
	      thisline, CLA.filterfile, "FILTER_ORDER");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;      
    }
   
  /* Allocate inverse sensing funtion filters */
  InputData.D=(REAL8IIRFilter *)LALMalloc(sizeof(REAL8IIRFilter)); 

  /* Allocate inverse sensing function filter */
  InputData.D->directCoef=NULL;
  InputData.D->recursCoef=NULL;
  InputData.D->history=NULL;

  LALDCreateVector(&status,&(InputData.D->directCoef),ND);
  LALDCreateVector(&status,&(InputData.D->recursCoef),ND);
  LALDCreateVector(&status,&(InputData.D->history),ND-1);

  for(l=0;l<ND;l++) InputData.D->directCoef->data[l]=0.0;
  for(l=0;l<ND;l++) InputData.D->recursCoef->data[l]=0.0;
  for(l=0;l<ND-1;l++) InputData.D->history->data[l]=0.0;
  
  for(n = 0; n < ND; n++)
    {
      /* read direct coeffs */
      i++;/*advance one line */
      thisline = Filters->lines->tokens[i];	/* get line i */
      InputData.D->directCoef->data[n]=strtod(thisline, &thisline);
    }

  /**-----------------------------------------------------------------------**/
  /* Read in actuation */
  i++; /*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  sscanf (thisline,"%s", aastr);
  if ( strcmp(aastr, "ACTUATION" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", 
	      thisline, CLA.filterfile, "ACTUATION");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;
    }
  /* Read number of sensing filters and their orders */
  i++;/*advance one line */
  thisline = Filters->lines->tokens[i];	/* get line i */
  NA=strtol(thisline, &thisline,10);
  if ( strcmp(thisline, " FILTER_ORDER" ) ) 
    {
      fprintf(stderr,"ERROR: Line (%s) of file %s is not properly terminated by '%s' marker!\n\n", 
	      thisline, CLA.filterfile, "FILTER_ORDER");
      LALDestroyParsedDataFile ( &status, &Filters ); TESTSTATUS( &status );
      return 1;      
    }
   
  /* Allocate inverse sensing funtion filters */
  InputData.A=(REAL8IIRFilter *)LALMalloc(sizeof(REAL8IIRFilter)); 

  /* Allocate inverse sensing function filter */
  InputData.A->directCoef=NULL;
  InputData.A->recursCoef=NULL;
  InputData.A->history=NULL;

  LALDCreateVector(&status,&(InputData.A->directCoef),NA);
  LALDCreateVector(&status,&(InputData.A->recursCoef),NA);
  LALDCreateVector(&status,&(InputData.A->history),NA-1);

  for(l=0;l<NA;l++) InputData.A->directCoef->data[l]=0.0;
  for(l=0;l<NA;l++) InputData.A->recursCoef->data[l]=0.0;
  for(l=0;l<NA-1;l++) InputData.A->history->data[l]=0.0;
  
  for(n = 0; n < NA; n++)
    {
      /* read direct coeffs */
      i++;/*advance one line */
      thisline = Filters->lines->tokens[i];	/* get line i */
      InputData.A->directCoef->data[n]=strtod(thisline, &thisline);
    }

  /**-----------------------------------------------------------------------**/

  LALDestroyParsedDataFile ( &status, &Filters );
  TESTSTATUS( &status );

  return 0;
}

/*******************************************************************************/

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA) 
{
  INT4 errflg=0;
  struct option long_options[] = {
    {"cal-line-freq",       required_argument, NULL,           'f'},
    {"factors-time",        required_argument, NULL,           't'},
    {"filters-file",        required_argument, NULL,           'F'},
    {"frame-cache",         required_argument, NULL,           'C'},
    {"exc-channel",         required_argument, NULL,           'E'},
    {"asq-channel",         required_argument, NULL,           'A'},
    {"darm-channel",        required_argument, NULL,           'D'},
    {"darmerr-channel",     required_argument, NULL,           'R'},
    {"gps-start-time",      required_argument, NULL,           's'},
    {"gps-end-time",        required_argument, NULL,           'e'},
    {"olg-re",              required_argument, NULL,           'i'},
    {"olg-im",              required_argument, NULL,           'j'},
    {"servo-re",            required_argument, NULL,           'k'},
    {"servo-im",            required_argument, NULL,           'l'},
    {"whitener-re",         required_argument, NULL,           'm'},
    {"whitener-im",         required_argument, NULL,           'n'},
    {"wings",               required_argument, NULL,           'o'},
    {"test-sensing",        no_argument, NULL,                 'r'},
    {"test-actuation",      no_argument, NULL,                 'c'},
    {"delta",               no_argument, NULL,                 'd'},
    {"no-factors",          no_argument, NULL,                 'u'},
    {"td-fir",              no_argument, NULL,                 'x'},
    {"output-factors",      no_argument, NULL,                 'y'},
    {"require-histories",   required_argument, NULL,           'H'},
    {"frame-type",          required_argument, NULL,           'T'},
    {"strain-channel",      required_argument, NULL,           'S'},
    {"data-dir",            required_argument, NULL,           'z'},
    {"rename-from",         required_argument, NULL,           'a'},
    {"rename-to",           required_argument, NULL,           'w'},
    {"help",                no_argument, NULL,                 'h' },
    {0, 0, 0, 0}
  };
  char args[] = "hrcduxyf:C:A:E:D:R:F:s:e:i:j:k:l:m:n:t:o:H:T:S:z:";
  
  /* Initialize default values */
  CLA->f=0.0;
  CLA->To=0.0;
  CLA->FrCacheFile=NULL;
  CLA->filterfile=NULL;
  CLA->exc_chan=NULL;
  CLA->darm_chan=NULL;
  CLA->darmerr_chan=NULL;
  CLA->asq_chan=NULL;
  CLA->GPSStart=0;
  CLA->GPSEnd=0;
  CLA->G0Re=0.0;
  CLA->G0Im=0.0;
  CLA->D0Re=0.0;
  CLA->D0Im=0.0;
  CLA->W0Re=0.0;
  CLA->W0Im=0.0;
  CLA->testsensing=0;
  CLA->testactuation=0;
  CLA->frametype=NULL;
  CLA->strainchannel=NULL;
  CLA->datadir=NULL;
  CLA->renamefrom=NULL;
  CLA->renameto=NULL;

  InputData.delta=0;
  InputData.wings=0;
  InputData.usefactors=1;
  InputData.fftconv=1;
  InputData.outalphas=0;

  /* Scan through list of command line arguments */
  while ( 1 )
  {
    int option_index = 0; /* getopt_long stores long option here */
    int c;

    c = getopt_long_only( argc, argv, args, long_options, &option_index );
    if ( c == -1 ) /* end of options */
      break;

    switch ( c )
    {
    case 'f':
      /* calibration line frequency */
      CLA->f=atof(optarg);
      break;
    case 't':
      /* factors integration time */
      CLA->To=atof(optarg);
      break;
    case 'C':
      /* name of frame cache file */
      CLA->FrCacheFile=optarg;
      break;
    case 'F':
      /* name of filter cache file */
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
    case 'R':
      /* name of darm err channel */
      CLA->darmerr_chan=optarg;
      break;    
    case 's':
      /* GPS start */
      CLA->GPSStart=atof(optarg);
      break;
    case 'e':
      /* GPS end */
      CLA->GPSEnd=atof(optarg);
      break;
    case 'i':
      /* real part of OLG */
      CLA->G0Re=atof(optarg);
      break;
    case 'j':
      /* imaginary part of OLG */
      CLA->G0Im=atof(optarg);
      break;
    case 'k':
      /*  real part of servo*/
      CLA->D0Re=atof(optarg);
      break;
    case 'l':
      /*  imaginary part of servo */
      CLA->D0Im=atof(optarg);
      break;
    case 'm':
      /*  real part of servo*/
      CLA->W0Re=atof(optarg);
      break;
    case 'n':
      /*  imaginary part of servo */
      CLA->W0Im=atof(optarg);
      break;
    case 'd':
      /*  use unit impulse*/
      InputData.delta=1;
      break;
    case 'r':
      /*  output residual signal */
      CLA->testsensing=1; 
      break;
    case 'c':
      /* output control signal */
      CLA->testactuation=1;
      break;
    case 'o':
      /*  wing size (at each end have this many extra seconds) */
      InputData.wings=atoi(optarg);
      break;
    case 'u':
      /* don't use calibration factors in teh strain computation */
      InputData.usefactors=0;
      break;
    case 'x':
      /* use time domain FIR filtering instead of FFT convolution */
      InputData.fftconv=0;
      break;
    case 'y':
      /* output calibration factors (the code will output them into the residual signal
	 hence the testsensing=1) */
      InputData.outalphas=1;
      CLA->testsensing=1;       
      break;
    case 'T':
      CLA->frametype=optarg;       
      break;
    case 'S':
      CLA->strainchannel=optarg;       
      break;
    case 'z':
      CLA->datadir=optarg;       
      break;
    case 'a':
      CLA->renamefrom=optarg;       
      break;
    case 'w':
      CLA->renameto=optarg;       
      break;
    case 'h':
      /* print usage/help message */
      fprintf(stdout,"Arguments are:\n");
      fprintf(stdout,"\tcal-line-freq (-f)\tFLOAT\t Calibration line frequency in Hz.\n");
      fprintf(stdout,"\tfactors-time (-t)\tFLOAT\t Factors integration time in seconds.\n");
      fprintf(stdout,"\tolg-re (-i)\tFLOAT\t Real part of the open loop gain at the calibration line frequency.\n");
      fprintf(stdout,"\tolg-im (-j)\tFLOAT\t Imaginary part of the open loop gain at the calibration line frequency.\n");
      fprintf(stdout,"\tservo-re (-k)\tFLOAT\t Real part of the digital filter at the calibration line frequency.\n");
      fprintf(stdout,"\tservo-im (-l)\tFLOAT\t Imaginary part of digital filter at the calibration line frequency.\n");
      fprintf(stdout,"\twhitener-re (-m)\tFLOAT\t Real part of the whitening filter at the calibration line frequency.\n");
      fprintf(stdout,"\twhitener-im (-n)\tFLOAT\t Imaginary part of whitening filter at the calibration line frequency.\n");
      fprintf(stdout,"\tgps-start-time (-s)\tINT\t GPS start time.\n");
      fprintf(stdout,"\tgps-end-time (-e)\tINT\t GPS end time.\n");
      fprintf(stdout,"\tfilters-file (-F)\tSTRING\t Name of file containing filters and histories.\n");
      fprintf(stdout,"\tframe-cache (-C)\tSTRING\t Name of frame cache file.\n");
      fprintf(stdout,"\tasq-channel (-A)\tSTRING\t AS_Q channel name (eg, L1:LSC-AS_Q).\n");
      fprintf(stdout,"\texc-channel (-E)\tSTRING\t Excitation channel name (eg, L1:LSC-ETMX_EXC_DAQ)\n");
      fprintf(stdout,"\tdarm-channel (-D)\tSTRING\t Darm channel name (eg, L1:LSC-DARM_CTRL)\n");
      fprintf(stdout,"\tdarmerr-channel (-R)\tSTRING\t Darm ERR channel name (eg, L1:LSC-DARM_ERR)\n");
      fprintf(stdout,"\twings (-o)\tINTEGER\t Size of wings in seconds.\n");
      fprintf(stdout,"\ttest-sensing (-r)\tFLAG\t Output residual strain only.\n");
      fprintf(stdout,"\ttest-actuation (-c)\tFLAG\t Output control strain only.\n");
      fprintf(stdout,"\tdelta (-d)\tFLAG\t Use unit impulse.\n");
      fprintf(stdout,"\tno-factors (-u)\tFLAG\t Do not use factors in strain computation.\n");
      fprintf(stdout,"\ttd-fir (-x)\tFLAG\t Use time-domain FIR filtering (default is FFT convolution).\n");
      fprintf(stdout,"\toutput-factors (-y)\tFLAG\t Outputs upsampled calibration factors time series.\n");
      fprintf(stdout,"\tframe-type (-T)\tSTRING\t Frame type to be written (eg, H1_RDS_C01_LX)\n");
      fprintf(stdout,"\tstrain-channel (-S)\tSTRING\t Strain channel name in frame (eg, H1:LSC-STRAIN)\n");
      fprintf(stdout,"\tdata-dir (-z)\tSTRING\t Ouput frame to this directory (eg, /tmp/S4/H1/H). Don't forget the H or L at the end!\n");
      fprintf(stdout,"\thelp (-h)\tFLAG\t This message\n");    
      exit(0);
      break;
    default:
      /* unrecognized option */
      errflg++;
      fprintf(stderr,"Unrecognized option argument %c\n",c);
      exit(1);
      break;
    }
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
   if(CLA->darmerr_chan == NULL)
    {
      fprintf(stderr,"No darm err channel specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }      
   if(CLA->asq_chan == NULL)
    {
      fprintf(stderr,"No asq channel specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }      
   if(CLA->frametype == NULL)
    {
      fprintf(stderr,"No frame type specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }      
   if(CLA->strainchannel == NULL)
    {
      fprintf(stderr,"No strain channel specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }      
   if(CLA->datadir == NULL)
    {
      fprintf(stderr,"No data directory specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }      

   if ( (InputData.wings < 2) || (InputData.wings%2 != 0) )
     {
       fprintf(stderr,"Overlap %d is either not >=2, or not exactly divisible by two. Exiting.",InputData.wings);
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

  /* Free input data */
  LALSDestroyVector(&status,&InputData.AS_Q.data);
  TESTSTATUS( &status );
  LALSDestroyVector(&status,&InputData.DARM.data);
  TESTSTATUS( &status );
  LALSDestroyVector(&status,&InputData.DARM_ERR.data);
  TESTSTATUS( &status );
  LALSDestroyVector(&status,&InputData.EXC.data);
  TESTSTATUS( &status );

  /* Free filters */
  LALDDestroyVector(&status,&InputData.Cinv->directCoef);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&InputData.Cinv->recursCoef);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&InputData.Cinv->history);   
  TESTSTATUS( &status );
  LALFree(InputData.Cinv);

  LALDDestroyVector(&status,&InputData.A->directCoef);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&InputData.A->recursCoef);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&InputData.A->history);   
  TESTSTATUS( &status );
  LALFree(InputData.A);

  LALDDestroyVector(&status,&InputData.D->directCoef);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&InputData.D->recursCoef);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&InputData.D->history);   
  TESTSTATUS( &status );
  LALFree(InputData.D);

  /* Free output data */
  LALDDestroyVector(&status,&OutputData.h.data);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&OutputData.hC.data);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&OutputData.hR.data);
  TESTSTATUS( &status );
  LALZDestroyVector(&status,&OutputData.alpha.data);
  TESTSTATUS( &status );
  LALZDestroyVector(&status,&OutputData.beta.data);
  TESTSTATUS( &status );
  LALZDestroyVector(&status,&OutputData.alphabeta.data);
  TESTSTATUS( &status );

  LALCheckMemoryLeaks();
 
  return 0;
}

/*******************************************************************************/
#endif
