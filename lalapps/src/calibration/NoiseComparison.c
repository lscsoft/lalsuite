/*********************************************************************************/
/* Noise comparison code between h(t) and  frequency domain calibrated DARM_ERR  */
/*                                                                               */
/*                                  X. Siemens                                   */
/*                                                                               */
/*                Caltech -- December 2006                  */
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
#include <lal/TimeFreqFFT.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>

#include <FrameL.h>
#include <series.h>

extern char *optarg;
extern int optind, opterr, optopt;

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 100; } else ((void)0)
#define MAXLINESRS   500000     /* Maximum # of lines in a Response file */


/***************************************************************************/

/* STRUCTURES */

struct CommandLineArgsTag {
  REAL8 f;                 /* Starting Frequency */
  REAL8 b;                 /* band */    
  REAL8 t;                 /* Time interval to compute the FFT */
  char *hoftFrCacheFile;       /* Frame cache file for h(t) */
  char *derrFrCacheFile;       /* Frame cache file for DARM_ERR */
  char *calFrCacheFile;       /* calibration frame cache file */
  char *derr_chan;         /* DARM_ERR channel name */ 
  char *hoft_chan;         /* h(t) channel name */
  char *noisefile;         /* output file for the noise */
  char *ResponseFile;         /* response function file */
  INT4 GPSStart;           /* Start and end GPS times of the segment to be compared*/
  INT4 GPSEnd;
} CommandLineArgs;

typedef struct ResponseFunctionTag
{
  REAL4 Frequency[MAXLINESRS];
  REAL4 Magnitude[MAXLINESRS];
  REAL4 Phase[MAXLINESRS];
  REAL4 re[MAXLINESRS];
  REAL4 im[MAXLINESRS];
} Response;

/***************************************************************************/

/* GLOBAL VARIABLES */

static LALStatus status;
INT4 lalDebugLevel=3;

FrCache *hoftframecache;                                           /* frame reading variables */
FrStream *hoftframestream=NULL;

FrCache *derrframecache;                                           /* frame reading variables */
FrStream *derrframestream=NULL;

LIGOTimeGPS gpsepoch;

INT4 duration;

static REAL4TimeSeries derr;
static REAL8TimeSeries hoft;
static FrChanIn chanin_derr;
static FrChanIn chanin_hoft;
REAL4Vector *derrwin=NULL,*hoftwin=NULL;

LALWindowParams winparams;
FILE *fpout=NULL;
COMPLEX16Vector *ffthtData = NULL;
COMPLEX8Vector *fftderrData = NULL;
REAL8FFTPlan *fftPlanDouble=NULL;
REAL4FFTPlan *fftPlan=NULL;
FrPos derrpos;
FrPos hoftpos;

Response R0,R;

/***************************************************************************/

/* FUNCTION PROTOTYPES */
/* */
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);
/* */
int Initialise(struct CommandLineArgsTag CLA);
/* */
int ReadResponse(struct CommandLineArgsTag CLA);
/* */
int ComputeNoise(struct CommandLineArgsTag CLA);   
/* */
int Finalise(void);

/************************************* MAIN PROGRAM *************************************/

int main(int argc,char *argv[])
{
  int i, N; 

  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;

  if (ReadResponse(CommandLineArgs)) return 3;

  duration = CommandLineArgs.GPSEnd-CommandLineArgs.GPSStart;

  N = duration / CommandLineArgs.t ;

  gpsepoch.gpsSeconds=CommandLineArgs.GPSStart;
  gpsepoch.gpsNanoSeconds=0;

  if(Initialise(CommandLineArgs)) return 1;
  
  for(i=0;i<N;i++)
    {
      gpsepoch.gpsSeconds=CommandLineArgs.GPSStart+i*CommandLineArgs.t;
      gpsepoch.gpsNanoSeconds=0;

      if(ComputeNoise(CommandLineArgs)) return 2;
    }
  
  if(Finalise()) return 1;
 
  return 0;
}

/************************************* MAIN PROGRAM ENDS *************************************/


/*  FUNCTIONS */
/*******************************************************************************/
int Initialise(struct CommandLineArgsTag CLA)
{

  /* create Frame cache, open frame stream and delete frame cache */
  LALFrCacheImport(&status,&derrframecache,CommandLineArgs.derrFrCacheFile);
  TESTSTATUS( &status );
  LALFrCacheOpen(&status,&derrframestream,derrframecache);
  TESTSTATUS( &status );
  LALDestroyFrCache(&status,&derrframecache);
  TESTSTATUS( &status );

  LALFrCacheImport(&status,&hoftframecache,CommandLineArgs.hoftFrCacheFile);
  TESTSTATUS( &status );
  LALFrCacheOpen(&status,&hoftframestream,hoftframecache);
  TESTSTATUS( &status );
  LALDestroyFrCache(&status,&hoftframecache);
  TESTSTATUS( &status );

  chanin_derr.type  = ADCDataChannel;
  chanin_hoft.type = ProcDataChannel;

  chanin_hoft.name  = CLA.hoft_chan;
  chanin_derr.name = CLA.derr_chan;

  /* Get channel time step size by calling LALFrGetREAL4TimeSeries */
  /* DARM_ERR */
  LALFrSeek(&status,&gpsepoch,derrframestream);
  TESTSTATUS( &status );
  LALFrGetPos(&status,&derrpos,derrframestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&derr,&chanin_derr,derrframestream);
  TESTSTATUS( &status );
  /* h(t) */
  LALFrSeek(&status,&gpsepoch,hoftframestream);
  TESTSTATUS( &status );
  LALFrGetPos(&status,&hoftpos,hoftframestream);
  TESTSTATUS( &status );
  LALFrGetREAL8TimeSeries(&status,&hoft,&chanin_hoft,hoftframestream);
  TESTSTATUS( &status );

  /* Allocate space for data vectors */
  LALCreateVector(&status,&derr.data,(UINT4)(CLA.t/derr.deltaT +0.5));
  TESTSTATUS( &status );
  LALDCreateVector(&status,&hoft.data,(UINT4)(CLA.t/hoft.deltaT +0.5));
  TESTSTATUS( &status );

  /* Create Window vectors */
  LALCreateVector(&status,&derrwin,(UINT4)(CLA.t/derr.deltaT +0.5));
  TESTSTATUS( &status );
  LALCreateVector(&status,&hoftwin,(UINT4)(CLA.t/hoft.deltaT +0.5));
  TESTSTATUS( &status );

  winparams.type=Hann;
   
  /* make windows  */
  /*  DARM_ERR */
  winparams.length=(INT4)(CLA.t/derr.deltaT +0.5);
  LALWindow(&status,derrwin,&winparams);
  TESTSTATUS( &status );
  
  /* h(t) */
  winparams.length=(INT4)(CLA.t/hoft.deltaT +0.5);
  LALWindow(&status,hoftwin,&winparams);
  TESTSTATUS( &status );

  /* Open output file */
  fpout=fopen(CLA.noisefile,"w");
  if (fpout==NULL) 
    {
      fprintf(stderr,"Could not open %s!\n",CLA.noisefile);
      return 1;
    }
  setvbuf( fpout, NULL, _IONBF, 0 );  /* This stops buffering */
  
  fprintf(fpout,"%s  GPS start time      h(t) noise     Calibrated DARM_ERR noise      Ratio\n", "%");

  LALCreateForwardREAL8FFTPlan( &status, &fftPlanDouble, hoft.data->length, 0 );
  TESTSTATUS( &status );
  LALCreateForwardREAL4FFTPlan( &status, &fftPlan, derr.data->length, 0 );
  TESTSTATUS( &status );

  LALZCreateVector( &status, &ffthtData, hoft.data->length / 2 + 1 );
  TESTSTATUS( &status );  

  LALCCreateVector( &status, &fftderrData, hoft.data->length / 2 + 1 );
  TESTSTATUS( &status );  
  
  return 0;
}
/*******************************************************************************/

int ComputeNoise(struct CommandLineArgsTag CLA)
{

  INT4 k;
  REAL8 mean_Sh_derr;
  REAL8 mean_Sh_hoft;

  /* Fill data vectors with data */
  /* DARM_ERR */
  LALFrSeek(&status,&gpsepoch,derrframestream);
  TESTSTATUS( &status );
  LALFrGetPos(&status,&derrpos,derrframestream);
  TESTSTATUS( &status );

  LALFrSetPos(&status,&derrpos,derrframestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&derr,&chanin_derr,derrframestream);
  TESTSTATUS( &status );
  /* h(t) */
  LALFrSeek(&status,&gpsepoch,hoftframestream);
  TESTSTATUS( &status );
  LALFrGetPos(&status,&hoftpos,hoftframestream);
  TESTSTATUS( &status );

  LALFrSetPos(&status,&hoftpos,hoftframestream);
  TESTSTATUS( &status );
  LALFrGetREAL8TimeSeries(&status,&hoft,&chanin_hoft,hoftframestream);
  TESTSTATUS( &status );

  /* Window the data */
  for(k=0;k<(INT4)(CLA.t/derr.deltaT +0.5);k++)
    {
      derr.data->data[k] *= 2.0*derrwin->data[k];
    }
  for(k=0;k<(INT4)(CLA.t/hoft.deltaT +0.5);k++)
    {
      hoft.data->data[k] *= 2.0*hoftwin->data[k];
    }

  /* FFT the data */
  XLALREAL8ForwardFFT( ffthtData, hoft.data, fftPlanDouble );
  XLALREAL4ForwardFFT( fftderrData, derr.data, fftPlan );

  /* Compute the noise */
  mean_Sh_derr = 0.0;
  mean_Sh_hoft = 0.0;
  {
    int firstbin=(INT4)(CLA.f*CLA.t+0.5);
    int nsamples=(INT4)(CLA.b*CLA.t+0.5);
    
    for (k=0; k<nsamples; k++)
    {
      REAL4 dr= fftderrData->data[k+firstbin].re;
      REAL4 di= fftderrData->data[k+firstbin].im;
      REAL8 dpsq;
      REAL4 caldr,caldi;

      REAL4 hr= ffthtData->data[k+firstbin].re;
      REAL4 hi= ffthtData->data[k+firstbin].im;
      REAL8 hpsq;

      /* Apply response to DARM_ERR */
      caldr=dr*R.re[k]-di*R.im[k];
      caldi=di*R.re[k]+dr*R.im[k];
      
      dpsq=pow(caldr,2.0)+pow(caldi,2);
      hpsq=pow(hr,2.0)+pow(hi,2);

      mean_Sh_derr += dpsq;
      mean_Sh_hoft += hpsq;
      
    }
    mean_Sh_derr *= 2.0*derr.deltaT/(REAL4)derr.data->length;
    mean_Sh_hoft *= 2.0*hoft.deltaT/(REAL4)hoft.data->length;
  }


  fprintf(fpout, "    %d         %e          %e          %e\n", gpsepoch.gpsSeconds, 
	  sqrt(mean_Sh_derr), sqrt(mean_Sh_hoft), sqrt(mean_Sh_hoft)/sqrt(mean_Sh_derr));
  
  return 0;
}


/*******************************************************************************/

int ReadResponse(struct CommandLineArgsTag CLA)
{
  char line[256];

  INT4 i,j;
  FILE *fpR;
  int nsamples=(INT4)(CLA.b*CLA.t+0.5);
  REAL4 f,df;

 /* ------ Open and read Response File ------ */
 i=0;
 fpR=fopen(CLA.ResponseFile,"r");
 if (fpR==NULL) 
   {
     fprintf(stderr,"Weird... %s doesn't exist!\n",CLA.ResponseFile);
     return 1;
   }
 while(fgets(line,sizeof(line),fpR))
   {
     if(*line == '#') continue;
     if(*line == '%') continue;
     if (i >= MAXLINESRS)
       {
	 fprintf(stderr,"Too many lines in file %s! Exiting... \n", CLA.ResponseFile);
	 return 1;
       }
     sscanf(line,"%e %e %e",&R0.Frequency[i],&R0.Magnitude[i],&R0.Phase[i]);
     R0.re[i]=R0.Magnitude[i]*cos(R0.Phase[i]);
     R0.im[i]=R0.Magnitude[i]*sin(R0.Phase[i]);
     i++;
   }
 /* -- close Response file -- */
 fclose(fpR);     

 /* Now compute a response appropriate for the frequency spacing of our data */
 i=0;
 for (j=0; j < nsamples;j++)               
    {
      REAL4 a,b;

      f = CLA.f  + j/CLA.t;

      while (i < MAXLINESRS-1 && f > R0.Frequency[i]) i++; 

      if(i == MAXLINESRS-1 && f > R0.Frequency[i])
	{
	  fprintf(stderr,"No calibration info for frequency %f!\n",f);
	  return 1;
	}
      /* checks order */

      /* Since now Rraw.Frequency[i-1] < f =< Rraw.Frequency[i] ||*/
      /* can compute a weighted average of the raw frequencies at i-1 and i */

      /* check both bounds!!!!! */
      if(f < R0.Frequency[i-1] || f > R0.Frequency[i])
	{
	  fprintf(stderr,"Frequency %f in SFT does not lie between two lines in Response file!\n",f);
	  return 1;
	}

      /* If the frequencies are closely spaced this may create dangerous floating point errors */
      df=R0.Frequency[i]-R0.Frequency[i-1];

      a=(f-R0.Frequency[i-1])/df;
      if (a>1.0) a=1.0;
      b=1.0-a;

      R.Frequency[j]=f;
      R.Magnitude[j]=a*R0.Magnitude[i]+b*R0.Magnitude[i-1];
      R.Phase[j]=a*R0.Phase[i]+b*R0.Phase[i-1];

      R.re[j]=a*R0.re[i]+b*R0.re[i-1];
      R.im[j]=a*R0.im[i]+b*R0.im[i-1];      
    }




  return 0;
}

/*******************************************************************************/

int Finalise(void)
{

  LALDestroyVector(&status,&derr.data);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&hoft.data);
  TESTSTATUS( &status );

  LALDestroyVector(&status,&derrwin);
  TESTSTATUS( &status );
  LALDestroyVector(&status,&hoftwin);
  TESTSTATUS( &status );

  LALZDestroyVector( &status, &ffthtData);
  LALCDestroyVector( &status, &fftderrData);

  LALDestroyREAL8FFTPlan( &status, &fftPlanDouble );
  LALDestroyREAL4FFTPlan( &status, &fftPlan );

  LALFrClose(&status,&derrframestream);
  TESTSTATUS( &status );

  LALFrClose(&status,&hoftframestream);
  TESTSTATUS( &status );

  fclose(fpout);
  
  LALCheckMemoryLeaks();

  return 0;
}

/*******************************************************************************/

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA) 
{

  INT4 errflg=0;
  struct option long_options[] = {
    {"freq",                 required_argument, NULL,           'a'},
    {"time",                 required_argument, NULL,           'b'},
    {"band",                 required_argument, NULL,           'c'},
    {"hoft-cache",           required_argument, NULL,           'd'},
    {"derr-cache",           required_argument, NULL,           'e'},
    {"cal-cache",            required_argument, NULL,           'f'},
    {"derr-channel",         required_argument, NULL,           'g'},
    {"hoft-channel",         required_argument, NULL,           'i'},
    {"gps-start-time",       required_argument, NULL,           'j'},
    {"gps-end-time",         required_argument, NULL,           'k'},
    {"output-file",          required_argument, NULL,           'l'},
    {"response-file",        required_argument, NULL,           'm'},
    {"help",                 no_argument, NULL,                 'h'},
    {0, 0, 0, 0}
  };
  char args[] = "ha:b:c:d:e:f:g:i:j:k:l:";

  /* Initialize default values */
  CLA->f=0.0;
  CLA->b=0.0;
  CLA->t=0.0;
  CLA->hoftFrCacheFile=NULL;
  CLA->derrFrCacheFile=NULL;
  CLA->calFrCacheFile=NULL;
  CLA->derr_chan=NULL;
  CLA->hoft_chan=NULL;
  CLA->noisefile=NULL;
  CLA->ResponseFile=NULL;
  CLA->GPSStart = 0;
  CLA->GPSEnd = 0;

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
    case 'a':
      CLA->f=atof(optarg);
      break;
    case 'b':
      CLA->t=atof(optarg);
      break;
    case 'c':
      CLA->b=atof(optarg);
      break;
    case 'd':
      CLA->hoftFrCacheFile=optarg;
      break;
    case 'e':
      CLA->derrFrCacheFile=optarg;
      break;    
    case 'f':
      /* name of as_q channel */
      CLA->calFrCacheFile=optarg;
      break;    
    case 'g':
      /* name of darm channel */
      CLA->derr_chan=optarg;
      break;    
    case 'i':
      /* name of darm err channel */
      CLA->hoft_chan=optarg;
      break;    
    case 'j':
      /* GPS start */
      CLA->GPSStart=atoi(optarg);
      break;
    case 'k':
      /* GPS end */
      CLA->GPSEnd=atof(optarg);
      break;
    case 'l':
      /* real part of OLG */
      CLA->noisefile=optarg;
      break;
    case 'm':
      /* real part of OLG */
      CLA->ResponseFile=optarg;
      break;
    case 'h':
      /* print usage/help message */
      fprintf(stdout,"Arguments are:\n");
      fprintf(stdout,"\t freq (-a)\tFLOAT\t Start frequency (Hz).\n");
      fprintf(stdout,"\t time (-b)\tFLOAT\t Integration time (s).\n");
      fprintf(stdout,"\t band (-c)\tFLOAT\t Frequency band (Hz).\n");
      fprintf(stdout,"\t hoft-cache (-d)\tFLOAT\t h(t) cache file name.\n");
      fprintf(stdout,"\t derr-cache (-e)\tFLOAT\t DARM_ERR cache file name.\n");
      fprintf(stdout,"\t cal-cache (-f)\tFLOAT\t Calibration cache file name.\n");
      fprintf(stdout,"\t derr-channel (-g)\tFLOAT\t DARM_ERR channel name.\n");
      fprintf(stdout,"\t hoft-channel (-i)\tFLOAT\t h(t) channel name.\n");
      fprintf(stdout,"\t gps-start-time (-j)\tINT\t GPS start time.\n");
      fprintf(stdout,"\t gps-end-time (-k)\tINT\t GPS end time.\n");
      fprintf(stdout,"\t ouput-file (-l)\tSTRING\t Name of output file.\n");
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
      fprintf(stderr,"No start frequency specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }      

  if(CLA->b == 0.0)
    {
      fprintf(stderr,"No frequency band specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }      

  if(CLA->t == 0.0)
    {
      fprintf(stderr,"No integration time specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }      
  if(CLA->GPSStart == 0.0)
    {
      fprintf(stderr,"No GPS start time specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }      
  if(CLA->GPSEnd == 0.0)
    {
      fprintf(stderr,"No GPS end time specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }     
  if(CLA->hoftFrCacheFile==NULL)
    {
      fprintf(stderr,"No h(t) frame cache file specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->derrFrCacheFile==NULL)
    {
      fprintf(stderr,"No DARM_ERR frame cache file specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->calFrCacheFile==NULL)
    {
      fprintf(stderr,"No calibration frame cache file specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->derr_chan==NULL)
    {
      fprintf(stderr,"No DARM_ERR channel name specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->hoft_chan==NULL)
    {
      fprintf(stderr,"No h(t) channel name specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->noisefile==NULL)
    {
      fprintf(stderr,"No output file specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
 
  /* Some sanity checks */
  if (CLA->b < 1/CLA->t)
    {
      fprintf(stderr,"Frequency band requested (%e) is smaller than frequency resolution[1/integration time] (%e)\n", CLA->b, 1.0/CLA->t);
      return 1;
    }

  return errflg;
}

/*******************************************************************************/
#endif
