/*********************************************************************************/
/*                         Calibration factor generation code                    */
/*                                                                               */
/*                                  X. Siemens                                   */
/*                                                                               */
/*                 Albert Einstein Institute/UWM - October 2003                  */
/*********************************************************************************/

#include <config.h>
#if !defined HAVE_LIBLALFRAME
#include <stdio.h>
int main(void) {fputs("disabled, no lal frame library support.\n", stderr);return 1;}
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

extern char *optarg;
extern int optind, opterr, optopt;

#define MAXLINESEGS 10000                 /* Maximum number of science segments */

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 100; } else ((void)0)


/***************************************************************************/

/* STRUCTURES */
struct CommandLineArgsTag {
  REAL8 f;                 /* Frequency of the calibration line */
  REAL8 t;                 /* Time interval to calculate factors */
  REAL8 G0Re;              /* Real part of open loop gain at cal line freq.*/
  REAL8 G0Im;              /* Imaginary part of open loop gain at cal line freq. */
  REAL8 D0Re;              /* Real part of digital filter at cal line freq.*/
  REAL8 D0Im;              /* Imaginary part of digital filter at cal line freq. */
  char *FrCacheFile;       /* Frame cache file */
  char *SegmentsFile;      /* Text file with the segments */
  char *exc_chan;          /* excitation channel name */    
  char *darm_chan;         /* darm channel name */ 
  char *asq_chan;          /* asq channel name */
  char *alphafile;         /* file to store factors */
} CommandLineArgs;

typedef
struct SegmentListTag {
  INT4 gpsstart;          /* GPS Seconds of start of segment group */
  INT4 gpsend;            /* number of segments starting at tgps */
  INT4 seglength;         /* length of segment in seconds */       
} SegmentList;


typedef 
struct GlobalVariablesTag {
COMPLEX16 Rf0,Cf0,Af0;              /* Response, sensing and actuation function values at the frequency of the calibration line */
SegmentList SL[MAXLINESEGS];        /* Structure containing science segement info */
INT4 numsegs;                       /* number of science segments */
LIGOTimeGPS gpsepoch;               /* Global variables epoch and duration used to calculate the calibration factors */
INT4 duration;                      /* they are set every science segment to GPS start time and duration of segment */
} GlobalVariables;

/***************************************************************************/

/* GLOBAL VARIABLES */

static LALStatus status;
INT4 lalDebugLevel=3;
FrCache *framecache;                                           /* frame reading variables */
FrStream *framestream=NULL;

GlobalVariables GV;   /* A bunch of stuff is stored in here; mainly to protect it from accidents */

/***************************************************************************/

/* FUNCTION PROTOTYPES */
/* Reads the command line */
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);

/* Reads time segment, cache, response,sensing and actuation files */
int ReadFiles(struct CommandLineArgsTag CLA);             

/* Produces a table of calibration factors for a science segment; 
the factors are calculated every interval CLA.t */
int GetFactors(struct CommandLineArgsTag CLA);   

/************************************* MAIN PROGRAM *************************************/

int main(int argc,char *argv[])
{
int i; 


  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;
  if (ReadFiles(CommandLineArgs)) return 3;

  /* create Frame cache, open frame stream and delete frame cache */
  LALFrCacheImport(&status,&framecache,CommandLineArgs.FrCacheFile);
  TESTSTATUS( &status );
  LALFrCacheOpen(&status,&framestream,framecache);
  TESTSTATUS( &status );
  LALDestroyFrCache(&status,&framecache);
  TESTSTATUS( &status );

  for(i=0;i<GV.numsegs;i++)
    {
      GV.gpsepoch.gpsSeconds=GV.SL[i].gpsstart; /* Set global variable epoch */
      GV.gpsepoch.gpsNanoSeconds=0;
      GV.duration=GV.SL[i].seglength;

      if(GetFactors(CommandLineArgs)) return 5;
    }
  

  LALFrClose(&status,&framestream);
  TESTSTATUS( &status );

  LALCheckMemoryLeaks();
 
  return 0;
}

/************************************* MAIN PROGRAM ENDS *************************************/


/*  FUNCTIONS */
/*******************************************************************************/

int GetFactors(struct CommandLineArgsTag CLA)
{

FrPos pos1;

static REAL4TimeSeries darm;
static REAL4TimeSeries asq;
static REAL4TimeSeries exc;

static FrChanIn chanin_darm;
static FrChanIn chanin_asq;
static FrChanIn chanin_exc;

CalFactors factors;
UpdateFactorsParams params;

REAL4Vector *asqwin=NULL,*excwin=NULL,*darmwin=NULL;  /* windows */

LALWindowParams winparams;

REAL4 magexc,magasq,magdarm;


INT4 k,m;
LIGOTimeGPS localgpsepoch=GV.gpsepoch; /* Local variable epoch used to calculate the calibration factors */
long double gtime=(long double)(localgpsepoch.gpsSeconds+(long double)localgpsepoch.gpsNanoSeconds*1E-9);
 
FILE *fpAlpha=NULL;

  chanin_asq.type  = ADCDataChannel;
  chanin_darm.type = ADCDataChannel;
  chanin_exc.type  = ADCDataChannel;

  chanin_asq.name  = CLA.asq_chan;
  chanin_darm.name = CLA.darm_chan;
  chanin_exc.name  = CLA.exc_chan; 

  /* Get channel time step size by calling LALFrGetREAL4TimeSeries */
  LALFrSeek(&status,&localgpsepoch,framestream);
  TESTSTATUS( &status );
  LALFrGetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&asq,&chanin_asq,framestream);
  TESTSTATUS( &status );
  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&darm,&chanin_darm,framestream);
  TESTSTATUS( &status );
  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&exc,&chanin_exc,framestream);
  TESTSTATUS( &status );


  /* Allocate space for data vectors */
  LALCreateVector(&status,&asq.data,(UINT4)(CLA.t/asq.deltaT +0.5));
  TESTSTATUS( &status );
  LALCreateVector(&status,&darm.data,(UINT4)(CLA.t/darm.deltaT +0.5));
  TESTSTATUS( &status );
  LALCreateVector(&status,&exc.data,(UINT4)(CLA.t/exc.deltaT +0.5));
  TESTSTATUS( &status );


  /* Create Window vectors */
  LALCreateVector(&status,&asqwin,(UINT4)(CLA.t/asq.deltaT +0.5));
  TESTSTATUS( &status );
  LALCreateVector(&status,&darmwin,(UINT4)(CLA.t/darm.deltaT +0.5));
  TESTSTATUS( &status );
  LALCreateVector(&status,&excwin,(UINT4)(CLA.t/exc.deltaT +0.5));
  TESTSTATUS( &status );

  winparams.type=Hann;
   
  /* windows for time domain channels */
  /* asq */
  winparams.length=(INT4)(CLA.t/asq.deltaT +0.5);
  LALWindow(&status,asqwin,&winparams);
  TESTSTATUS( &status );
  
  /* darm */
  winparams.length=(INT4)(CLA.t/darm.deltaT +0.5);
  LALWindow(&status,darmwin,&winparams);
  TESTSTATUS( &status );

  /* exc */
  winparams.length=(INT4)(CLA.t/exc.deltaT +0.5);
  LALWindow(&status,excwin,&winparams);
  TESTSTATUS( &status );

  /* Open user input factors file */
  fpAlpha=fopen(CLA.alphafile,"w");
  if (fpAlpha==NULL) 
    {
      fprintf(stderr,"Could not open %s!\n",CLA.alphafile);
      return 1;
    }
  setvbuf( fpAlpha, NULL, _IONBF, 0 );  /* This stops buffering -- From Duncan Brown*/

  fprintf(fpAlpha,"%s Re(Go) = %e, Im(Go)= %e \n", "%",CLA.G0Re, CLA.G0Im);
  fprintf(fpAlpha,"%s Re(Do) = %e, Im(Do)= %e \n", "%",CLA.D0Re, CLA.D0Im);
  fprintf(fpAlpha,"%s  GPS Time          Re(a)    Im(a)     Re(b)    Im(b)    Re(ab)    Im(ab)   Re(AS_Q) Im(AS_Q) Re(DARM) Im(DARM) Re(EXC)  Im(EXC)\n", "%");

  for(m=0;m < (INT4)(GV.duration/CLA.t);m++)
    {
      /* Fill data vectors with data */
      LALFrSeek(&status,&localgpsepoch,framestream);
      TESTSTATUS( &status );
      LALFrGetPos(&status,&pos1,framestream);
      TESTSTATUS( &status );

      LALFrSetPos(&status,&pos1,framestream);
      TESTSTATUS( &status );
      LALFrGetREAL4TimeSeries(&status,&asq,&chanin_asq,framestream);
      TESTSTATUS( &status );

      LALFrSetPos(&status,&pos1,framestream);
      TESTSTATUS( &status );
      LALFrGetREAL4TimeSeries(&status,&darm,&chanin_darm,framestream);
      TESTSTATUS( &status );

      LALFrSetPos(&status,&pos1,framestream);
      TESTSTATUS( &status );
      LALFrGetREAL4TimeSeries(&status,&exc,&chanin_exc,framestream);
      TESTSTATUS( &status );

      /* Window the data */
      for(k=0;k<(INT4)(CLA.t/asq.deltaT +0.5);k++)
	{
	  asq.data->data[k] *= 2.0*asqwin->data[k];
	}
      for(k=0;k<(INT4)(CLA.t/darm.deltaT +0.5);k++)
	{
	  darm.data->data[k] *= 2.0*darmwin->data[k];
	}
      for(k=0;k<(INT4)(CLA.t/exc.deltaT +0.5);k++)
	{
	  exc.data->data[k] *= 2.0*excwin->data[k];
	}
      
      /* set params to call LALComputeCalibrationFactors */
      params.darmCtrl = &darm;
      params.asQ = &asq;
      params.exc = &exc;
      params.lineFrequency = CLA.f;
      params.openloop.re =  CLA.G0Re;
      params.openloop.im =  CLA.G0Im;
      params.digital.re = CLA.D0Re;
      params.digital.im = CLA.D0Im;
	  
      LALComputeCalibrationFactors(&status,&factors,&params);
      TESTSTATUS( &status );

      magexc=sqrt(factors.exc.re*factors.exc.re+factors.exc.im*factors.exc.im)*2/CLA.t;
      magasq=sqrt(factors.asq.re*factors.asq.re+factors.asq.im*factors.asq.im)*2/CLA.t;
      magdarm=sqrt(factors.darm.re*factors.darm.re+factors.darm.im*factors.darm.im)*2/CLA.t;

      fprintf(fpAlpha,"%18.9Lf %f %f %f %f %f %f %f %f %f %f %f %f \n",gtime,
	      factors.alpha.re,factors.alpha.im,
	      factors.beta.re,factors.beta.im,
	      factors.alphabeta.re,factors.alphabeta.im,
	      factors.asq.re*2/CLA.t,factors.asq.im*2/CLA.t,
	      factors.darm.re*2/CLA.t,factors.darm.im*2/CLA.t,
	      factors.exc.re*2/CLA.t,factors.exc.im*2/CLA.t);
 
      gtime += CLA.t;	
      localgpsepoch.gpsSeconds = (INT4)gtime;
      localgpsepoch.gpsNanoSeconds = (INT4)((gtime-(INT4)gtime)*1E+09);      
    }

  /* Clean up */
  LALDestroyVector(&status,&darm.data);
  TESTSTATUS( &status );
  LALDestroyVector(&status,&exc.data);
  TESTSTATUS( &status );
  LALDestroyVector(&status,&asq.data);
  TESTSTATUS( &status );

  LALDestroyVector(&status,&asqwin);
  TESTSTATUS( &status );
  LALDestroyVector(&status,&darmwin);
  TESTSTATUS( &status );
  LALDestroyVector(&status,&excwin);
  TESTSTATUS( &status );

  fclose(fpAlpha);
  
  return 0;
}

/*******************************************************************************/

int ReadFiles(struct CommandLineArgsTag CLA)
{
  char line[256];
  INT4 i;
  FILE *fpSeg;

 /* ------ Open and read Segment file ------ */
 i=0;
 fpSeg=fopen(CLA.SegmentsFile,"r");
 if (fpSeg==NULL) 
   {
     fprintf(stderr,"That's weird... %s doesn't exist!\n",CLA.SegmentsFile);
     return 1;
   }
 while(fgets(line,sizeof(line),fpSeg))
   {
     if(*line == '#') continue;
     if(*line == '%') continue;
     if (i > MAXLINESEGS-1)
       {
	 fprintf(stderr,"Too many lines in file %s! Exiting... \n", CLA.SegmentsFile);
	 return 1;
       }
     sscanf(line,"%d %d %d",&GV.SL[i].gpsstart,&GV.SL[i].gpsend,&GV.SL[i].seglength);
     i++;
   }
 GV.numsegs=i;
 fclose(fpSeg);     
 /* -- close Sensing file -- */

 return 0;
}

/*******************************************************************************/

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA) 
{
  INT4 c, errflg = 0;
  optarg = NULL;
  
  /* Initialize default values */
  CLA->f=0.0;
  CLA->t=0.0;
  CLA->G0Re=0.0;
  CLA->G0Im=0.0;
  CLA->D0Re=0.0;
  CLA->D0Im=0.0;
  CLA->FrCacheFile=NULL;
  CLA->SegmentsFile=NULL;
  CLA->alphafile=NULL;
  CLA->exc_chan=NULL;
  CLA->darm_chan=NULL;
  CLA->asq_chan=NULL;
  CLA->alphafile=NULL;

  /* Scan through list of command line arguments */
  while (!errflg && ((c = getopt(argc, argv,"hf:F:S:A:E:D:b:t:i:j:k:l:"))!=-1))
    switch (c) {
    case 'f':
      /* calibration line frequency */
      CLA->f=atof(optarg);
      break;
    case 't':
      /* calibration line frequency */
      CLA->t=atof(optarg);
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
    case 'F':
      /* name of frame cache file */
      CLA->FrCacheFile=optarg;
      break;
    case 'S':
      /* name of segments file */
      CLA->SegmentsFile=optarg;
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
    case 'b':
      /* name of darm channel */
      CLA->alphafile=optarg;
      break;    
   case 'h':
      /* print usage/help message */
      fprintf(stdout,"All arguments are required. They are:\n");
      fprintf(stdout,"\t-f\tFLOAT\t Calibration line frequency in Hz.\n");
      fprintf(stdout,"\t-t\tFLOAT\t Time interval to calculate factors in seconds (>0.0005).\n");
      fprintf(stdout,"\t-i\tFLOAT\t Real part of the open loop gain at the calibration line frequency.\n");
      fprintf(stdout,"\t-j\tFLOAT\t Imaginary part of the open loop gain at the calibration line frequency.\n");
      fprintf(stdout,"\t-k\tFLOAT\t Real part of the digital filter at the calibration line frequency.\n");
      fprintf(stdout,"\t-l\tFLOAT\t Imaginary part of digital filter at the calibration line frequency.\n");
      fprintf(stdout,"\t-F\tSTRING\t Name of frame cache file.\n");
      fprintf(stdout,"\t-S\tSTRING\t Name of segment list file.\n");
      fprintf(stdout,"\t-A\tSTRING\t AS_Q channel name (eg, L1:LSC-AS_Q).\n");
      fprintf(stdout,"\t-E\tSTRING\t Excitation channel name (eg, L1:LSC-ETMX_EXC_DAQ)\n");
      fprintf(stdout,"\t-D\tSTRING\t Darm channel name (eg, L1:LSC-DARM_CTRL)\n");
      fprintf(stdout,"\t-b\tSTRING\t Output file for calibration factors.\n");
      exit(0);
      break;
    default:
      /* unrecognized option */
      errflg++;
      fprintf(stderr,"Unrecognized option argument %c\n",c);
      exit(1);
      break;
    }

  if(CLA->f == 0)
    {
      fprintf(stderr,"No calibration line frequency specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }      
  if(CLA->t == 0)
    {
      fprintf(stderr,"No time interval to calculate factors.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }      
  if(CLA->t < 0.0005)
    {
      fprintf(stderr,"Time interval to calculate factors too small (<0.0005).\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }      
  if(CLA->G0Re == 0.0 )
    {
      fprintf(stderr,"No real part of open loop gain specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }      
  if(CLA->G0Im == 0.0 )
    {
      fprintf(stderr,"No imaginary part of open loop gain specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }      
  if(CLA->D0Re == 0.0 )
    {
      fprintf(stderr,"No real part of digital filter specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }      
  if(CLA->D0Im == 0.0 )
    {
      fprintf(stderr,"No imaginary part of digital filter specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }      
  if(CLA->alphafile == NULL)
    {
      fprintf(stderr,"No output factors file specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }      
  if(CLA->FrCacheFile == NULL)
    {
      fprintf(stderr,"No frame cache file specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }      
  if(CLA->SegmentsFile == NULL)
    {
      fprintf(stderr,"No segments file specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }      
   if(CLA->exc_chan == NULL)
    {
      fprintf(stderr,"No excitation channel specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }      
   if(CLA->darm_chan == NULL)
    {
      fprintf(stderr,"No darm channel specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }      
   if(CLA->asq_chan == NULL)
    {
      fprintf(stderr,"No asq channel specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }      

  return errflg;
}

/*******************************************************************************/
#endif
