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

#define MAXLINERS 76800                   /* Max lines read in Freq Response files */
#define MAXLINESEGS 10000                 /* Maximum number of science segments */

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 100; } else ((void)0)

/* run ./lalapps_ComputeFactors -b Factors.0 -f 973.3 -F S2-H1-Cache.0 -S NewSeg.0 -r R.txt -c C.txt -a A.txt -A H1:LSC-AS_Q -E H1:LSC-ETMX_EXC_DAQ -D H1:LSC-DARM_CTRL -m -0.493 */

/***************************************************************************/
/* Complex division routine -- used to calculate beta from alpha and alpha*beta */
static COMPLEX16 *cdiv( COMPLEX16 *pc, COMPLEX16 *pa, COMPLEX16 *pb )
{
  COMPLEX16 a = *pa;
  COMPLEX16 b = *pb;
  COMPLEX16 c;
  REAL8 rat;
  REAL8 den;
  if ( fabs( b.re ) > fabs( b.im ) )
  {
    rat = b.im / b.re;
    den = b.re + rat * b.im;
    c.re = ( a.re + rat * a.im ) / den;
    c.im = ( a.im - rat * a.re ) / den;
  }
  else
  {
    rat = b.re / b.im;
    den = b.im + rat * b.re;
    c.re = ( a.re * rat + a.im ) / den;
    c.im = ( a.im * rat - a.re ) / den;
  }
  *pc = c;
  return pc;
}

/***************************************************************************/

/* STRUCTURES */
struct CommandLineArgsTag {
  REAL8 f;                 /* Frequency of the calibration line */
  REAL8 t;                 /* Time interval to calculate factors */
  REAL8 mu;                /* Fraction of actuation seen by excitation channel */
  char *RFile;             /* Text file with the response funtion */
  char *CFile;             /* Text file with the sensing function */
  char *AFile;             /* Text file with the actuation function */
  char *FrCacheFile;       /* Frame cache file */
  char *SegmentsFile;      /* Text file with the segments */
  char *exc_chan;          /* excitation channel name */    
  char *darm_chan;         /* darm channel name */ 
  char *asq_chan;          /* asq channel name */
  char *alphafile;         /* asq channel name */
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

static REAL4TimeSeries gain;
static REAL4TimeSeries icm;


static FrChanIn chanin_darm;
static FrChanIn chanin_asq;
static FrChanIn chanin_exc;

static FrChanIn chanin_gain;
static FrChanIn chanin_icm;

CalFactors factors;
UpdateFactorsParams params;

REAL4Vector *asqwin=NULL,*excwin=NULL,*darmwin=NULL;  /* windows */

LALWindowParams winparams;

REAL4 g,i,magexc,magasq,magdarm;


COMPLEX16 be;
INT4 k,m;
LIGOTimeGPS localgpsepoch=GV.gpsepoch; /* Local variable epoch used to calculate the calibration factors */
long double gtime=(long double)(localgpsepoch.gpsSeconds+(long double)localgpsepoch.gpsNanoSeconds*1E-9);
 
FILE *fpAlpha=NULL;

  chanin_asq.type  = ADCDataChannel;
  chanin_darm.type = ADCDataChannel;
  chanin_exc.type  = ADCDataChannel;

  chanin_gain.type = ADCDataChannel;
  chanin_icm.type  = ADCDataChannel;

  chanin_asq.name  = CLA.asq_chan;
  chanin_darm.name = CLA.darm_chan;
  chanin_exc.name  = CLA.exc_chan; 

  chanin_gain.name = "L1:LSC-DARM_GAIN";
  chanin_icm.name  = "L1:LSC-ICMTRX_01";

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

/*   LALFrSetPos(&status,&pos1,framestream); */
/*   TESTSTATUS( &status ); */
/*   LALFrGetREAL4TimeSeries(&status,&gain,&chanin_gain,framestream); */
/*   TESTSTATUS( &status ); */
/*   LALFrSetPos(&status,&pos1,framestream); */
/*   TESTSTATUS( &status ); */
/*   LALFrGetREAL4TimeSeries(&status,&icm,&chanin_icm,framestream); */
/*   TESTSTATUS( &status ); */

  /* Allocate space for data vectors */
  LALCreateVector(&status,&asq.data,(UINT4)(CLA.t/asq.deltaT +0.5));
  TESTSTATUS( &status );
  LALCreateVector(&status,&darm.data,(UINT4)(CLA.t/darm.deltaT +0.5));
  TESTSTATUS( &status );
  LALCreateVector(&status,&exc.data,(UINT4)(CLA.t/exc.deltaT +0.5));
  TESTSTATUS( &status );

/*   LALCreateVector(&status,&gain.data,(UINT4)(CLA.t/gain.deltaT +0.5)); */
/*   TESTSTATUS( &status ); */
/*   LALCreateVector(&status,&icm.data,(UINT4)(CLA.t/icm.deltaT +0.5)); */
/*   TESTSTATUS( &status ); */

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
  fpAlpha=fopen(CLA.alphafile,"a");
  if (fpAlpha==NULL) 
    {
      fprintf(stderr,"Could not open %s!\n",CLA.alphafile);
      return 1;
    }
  setvbuf( fpAlpha, NULL, _IONBF, 0 );  /* This stops buffering -- From Duncan Brown*/

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

/*       LALFrSetPos(&status,&pos1,framestream); */
/*       TESTSTATUS( &status ); */
/*       LALFrGetREAL4TimeSeries(&status,&gain,&chanin_gain,framestream); */
/*       TESTSTATUS( &status ); */

/*       LALFrSetPos(&status,&pos1,framestream); */
/*       TESTSTATUS( &status ); */
/*       LALFrGetREAL4TimeSeries(&status,&icm,&chanin_icm,framestream); */
/*       TESTSTATUS( &status ); */

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
      
      /* Average gain and icm channels */
      g=0.0;
/*       for(k=0;k<(INT4)(CLA.t/gain.deltaT +0.5);k++) */
/* 	{ */
/* 	  g += gain.data->data[k]/((INT4)(CLA.t/gain.deltaT+0.5)); */
/* 	} */
      i=0.0;
/*       for(k=0;k<(INT4)(CLA.t/icm.deltaT +0.5);k++) */
/* 	{ */
/* 	  i += icm.data->data[k]/((INT4)(CLA.t/icm.deltaT +0.5)); */
/* 	} */

      /* set params to call LALComputeCalibrationFactors */
      params.darmCtrl = &darm;
      params.asQ = &asq;
      params.exc = &exc;
      params.lineFrequency = CLA.f;
      params.mu = CLA.mu;
      params.actuationFactor.re = GV.Af0.re ;
      params.actuationFactor.im = GV.Af0.im ;
      params.responseFactor = GV.Rf0;
      params.sensingFactor = GV.Cf0;
	  
      LALComputeCalibrationFactors(&status,&factors,&params);
      TESTSTATUS( &status );

      cdiv(&be,&factors.alphabeta,&factors.alpha);

      magexc=sqrt(factors.exc.re*factors.exc.re+factors.exc.im*factors.exc.im)*2/CLA.t;
      magasq=sqrt(factors.asq.re*factors.asq.re+factors.asq.im*factors.asq.im)*2/CLA.t;
      magdarm=sqrt(factors.darm.re*factors.darm.re+factors.darm.im*factors.darm.im)*2/CLA.t;

      fprintf(fpAlpha,"%18.9Lf %f %f %f %f %f %f %f %f %f %f \n",
	      gtime,factors.alpha.re,factors.alpha.im,be.re,be.im,
	      factors.alphabeta.re,factors.alphabeta.im,g*i,magexc,magasq,magdarm);

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
  FILE *fpS,*fpR,*fpA,*fpSeg;
  REAL8 Cmag,Cphase,Rmag,Rphase,Amag,Aphase,freq,x,y;
  static COMPLEX16FrequencySeries R0;
  static COMPLEX16FrequencySeries C0;
  static COMPLEX16FrequencySeries A0;

  /* Allocate space for response and sensing functions; just enough to read first 1200Hz */
  LALZCreateVector( &status, &R0.data, MAXLINERS);
  TESTSTATUS( &status );
  LALZCreateVector( &status, &C0.data, MAXLINERS);
  TESTSTATUS( &status );
  LALZCreateVector( &status, &A0.data, MAXLINERS);
  TESTSTATUS( &status );

  /* Fill in R0, C0 data */
  R0.f0=0.0;
  R0.deltaF=1.0/64.0;                   /*ACHTUNG: HARDWIRED !!*/
  C0.f0=0.0;
  C0.deltaF=1.0/64.0;                   /*ACHTUNG: HARDWIRED !!*/
  A0.f0=0.0;
  A0.deltaF=1.0/64.0;                   /*ACHTUNG: HARDWIRED !!*/

 /* This is kinda messy... Unfortunately there's no good way of doing this */ 
 /* ------ Open and read Sensing file ------ */
 i=0;
 fpS=fopen(CLA.CFile,"r");
 if (fpS==NULL) 
   {
     fprintf(stderr,"That's weird... %s doesn't exist!\n",CLA.CFile);
     return 1;
   }
 while(fgets(line,sizeof(line),fpS))
   {
     if(*line == '#') continue;
     if(*line == '%') continue;
     if (i > MAXLINERS-1)
       {
	 /* done reading file */
	 break;
       }
     sscanf(line,"%le %le %le",&freq,&Cmag,&Cphase);
     C0.data->data[i].re=Cmag*cos(Cphase);
     C0.data->data[i].im=Cmag*sin(Cphase);
     i++;
   }
 fclose(fpS);     
 /* -- close Sensing file -- */

 /* ------ Open and read Response file ------ */
 i=0;
 fpR=fopen(CLA.RFile,"r");
 if (fpR==NULL) 
   {
     fprintf(stderr,"That's weird... %s doesn't exist!\n",CLA.RFile);
     return 1;
   }
 while(fgets(line,sizeof(line),fpR))
   {
     if(*line == '#') continue;
     if(*line == '%') continue;
     if (i > MAXLINERS-1)
       {
	 /* done reading file */
	 break;
       }
     sscanf(line,"%le %le %le",&freq,&Rmag,&Rphase);
     R0.data->data[i].re=Rmag*cos(Rphase);
     R0.data->data[i].im=Rmag*sin(Rphase);
     i++;
   }
 fclose(fpR);     
 /* -- close Sensing file -- */

 /* ------ Open and read Response file ------ */
 i=0;
 fpA=fopen(CLA.AFile,"r");
 if (fpA==NULL) 
   {
     fprintf(stderr,"That's weird... %s doesn't exist!\n",CLA.AFile);
     return 1;
   }
 while(fgets(line,sizeof(line),fpA))
   {
     if(*line == '#') continue;
     if(*line == '%') continue;
     if (i > MAXLINERS-1)
       {
	 /* done reading file */
	 break;
       }
     sscanf(line,"%le %le %le",&freq,&Amag,&Aphase);
     A0.data->data[i].re=Amag*cos(Aphase);
     A0.data->data[i].im=Amag*sin(Aphase);
     i++;
   }
 fclose(fpA);     
 /* -- close Sensing file -- */

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

  /* compute global varibles GV.Cf0, GV.Af0 and GV.Rf0 at correct frequency */
  /* use linear interpolation */

  x = modf( CLA.f / R0.deltaF, &y );
  i = floor( y );
  GV.Rf0.re  = ( 1 - x ) * R0.data->data[i].re;
  GV.Rf0.re += x * R0.data->data[i].re;
  GV.Rf0.im  = ( 1 - x ) * R0.data->data[i].im;
  GV.Rf0.im += x * R0.data->data[i].im;
  x = modf( CLA.f / C0.deltaF, &y );
  i = floor( y );
  GV.Cf0.re  = ( 1 - x ) * C0.data->data[i].re;
  GV.Cf0.re += x * C0.data->data[i].re;
  GV.Cf0.im  = ( 1 - x ) * C0.data->data[i].im;
  GV.Cf0.im += x * C0.data->data[i].im;
  x = modf( CLA.f / A0.deltaF, &y );
  i = floor( y );
  GV.Af0.re  = ( 1 - x ) * A0.data->data[i].re;
  GV.Af0.re += x * A0.data->data[i].re;
  GV.Af0.im  = ( 1 - x ) * A0.data->data[i].im;
  GV.Af0.im += x * A0.data->data[i].im;

  /* clean up */
 LALZDestroyVector(&status,&R0.data);
 TESTSTATUS( &status );
 LALZDestroyVector(&status,&C0.data);
 TESTSTATUS( &status );
 LALZDestroyVector(&status,&A0.data);
 TESTSTATUS( &status );

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
  CLA->mu=0.0;  
  CLA->RFile=NULL;
  CLA->CFile=NULL;   
  CLA->AFile=NULL;   
  CLA->FrCacheFile=NULL;
  CLA->SegmentsFile=NULL;
  CLA->exc_chan=NULL;
  CLA->darm_chan=NULL;
  CLA->asq_chan=NULL;
  CLA->alphafile=NULL;

  /* Scan through list of command line arguments */
  while (!errflg && ((c = getopt(argc, argv,"hf:F:r:S:c:A:E:D:a:m:b:t:"))!=-1))
    switch (c) {
    case 'f':
      /* calibration line frequency */
      CLA->f=atof(optarg);
      break;
    case 't':
      /* calibration line frequency */
      CLA->t=atof(optarg);
      break;
    case 'm':
      /* darm to etmx output matrix value */
      CLA->mu=atof(optarg);
      break;
    case 'F':
      /* name of frame cache file */
      CLA->FrCacheFile=optarg;
      break;
    case 'r':
      /* name of response file */
      CLA->RFile=optarg;
      break;
    case 'c':
      /* name of sensing file */
      CLA->CFile=optarg;
      break;
    case 'a':
      /* name of actuation file */
      CLA->AFile=optarg;
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
      fprintf(stdout,"All arguments are required except -b. They are:\n");
      fprintf(stdout,"\t-f\tFLOAT\t Calibration line frequency in Hz.\n");
      fprintf(stdout,"\t-t\tFLOAT\t Time interval to calculate factors in seconds (>0.0005).\n");
      fprintf(stdout,"\t-m\tFLOAT\t Fraction of the actuation seen by excitation channel (=-1 if injected into darm (eg S3, or Gx/(Ky Gy- Kx Gx)...).\n");
      fprintf(stdout,"\t-F\tSTRING\t Name of frame cache file.\n");
      fprintf(stdout,"\t-r\tSTRING\t Name of response function file.\n");
      fprintf(stdout,"\t-c\tSTRING\t Name of sensing function file.\n");
      fprintf(stdout,"\t-a\tSTRING\t Name of actuation function file.\n");
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
 if(CLA->mu == 0)
    {
      fprintf(stderr,"No value of the output matrix to x arm specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }      
  if(CLA->CFile == NULL)
    {
      fprintf(stderr,"No sensing function file specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }      
  if(CLA->RFile == NULL)
    {
      fprintf(stderr,"No response function file specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }      
  if(CLA->AFile == NULL)
    {
      fprintf(stderr,"No actuation function file specified.\n");
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
