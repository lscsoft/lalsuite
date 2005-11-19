/*********************************************************************************/
/*                                                                               */
/* File: MakeSFTs.c                                                              */
/* Purpose: generate SFTs                                                        */
/* Origin: first written by Xavier Siemens, UWM - May 2005,                      */
/*         based on make_sfts.c by Bruce Allen                                   */
/* $Id$                     */
/*********************************************************************************/

/* REVISIONS: */
/* 11/02/05 gam; To save memory, change lalDebugLevel to 0; note when this is set to 3 that 3x the memory is used! */
/* 11/02/05 gam; To save memory, do not save the window function in memory; just recompute this for each SFT. */
/* 11/02/05 gam; To save memory, do not hold dataSingle.data and vtilde in memory at the same time. */
/* 11/03/05 gam; Add TRACKMEMUSE preprocessor flag and function for tracking memory usage; copied from make_sfts.c */
/* 11/19/05 gam; Reorganize code so that function appear in order called */
/* 11/19/05 gam; Add command line option to use single rather than double precision; will use LALDButterworthREAL4TimeSeries in single precision case. */
/* 11/19/05 gam; Rename vtilde fftDataDouble; add in fftDatasingle; rename fplan fftPlanDouble; add in fftPlanSingle */

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
#include <lal/TimeFreqFFT.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>

extern char *optarg;
extern int optind, opterr, optopt;

/* track memory usage under linux */
#define TRACKMEMUSE 0

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 100; } else ((void)0)

#define FMIN  (48.0)
#define DF    (2000.0)

/* This repeatedly tries to re-open a file.  This is useful on a
   cluster because automount may fail or time out.  This allows a
   number of tries with a bit of rest in between.
*/
FILE* tryopen(char *name, char *mode){
  int count=0;
  FILE *fp;
  
  while (!(fp=fopen(name, mode))){
    fprintf(stderr,"Unable to open file %s in mode %s.  Will retry...\n", name, mode);
    if (count++<10)
      sleep(10);
    else
      exit(3);
  }
  return fp;
}


/***************************************************************************/

/* STRUCTURES */
struct CommandLineArgsTag {
  REAL8 HPf;              /* High pass filtering frequency */
  INT4 T;                 /* SFT duration */
  INT4 GPSStart;
  INT4 GPSEnd;
  BOOLEAN htdata;          /* flag that indicates we're doing h(t) data */
  char *FrCacheFile;       /* Frame cache file */
  char *ChannelName;
  char *SFTpath;           /* path to SFT file location */
  BOOLEAN useSingle;       /* 11/19/05 gam; use single rather than double precision */
} CommandLineArgs;

struct headertag {
  REAL8 endian;
  INT4  gps_sec;
  INT4  gps_nsec;
  REAL8 tbase;
  INT4  firstfreqindex;
  INT4  nsamples;
} header;


/***************************************************************************/

/* GLOBAL VARIABLES */

static LALStatus status;
INT4 lalDebugLevel=0;        /* 11/02/05 gam; changed from 3 to 0. */
FrCache *framecache;         /* frame reading variables */
FrStream *framestream=NULL;

REAL8TimeSeries dataDouble;
REAL4TimeSeries dataSingle;
INT4 SegmentDuration;
LIGOTimeGPS gpsepoch;

/* REAL8Vector *window; */ /* 11/02/05 gam */

REAL8FFTPlan *fftPlanDouble;           /* fft plan and data container, double precision case */
COMPLEX16Vector *fftDataDouble = NULL;

REAL4FFTPlan *fftPlanSingle;           /* 11/19/05 gam; fft plan and data container, single precision case */
COMPLEX8Vector *fftDataSingle = NULL;

/***************************************************************************/

/* FUNCTION PROTOTYPES */
/* Reads the command line */
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);

/* Allocates space for data time series */
int AllocateData(struct CommandLineArgsTag CLA);

/* Reads data */
int ReadData(struct CommandLineArgsTag CLA);

/* High passes data */
int HighPass(struct CommandLineArgsTag CLA);

/* Windows data */
int WindowData(struct CommandLineArgsTag CLA);

/* create an SFT */
int CreateSFT(struct CommandLineArgsTag CLA);

/* writes out an SFT */
int WriteSFT(struct CommandLineArgsTag CLA);

/* Frees the memory */
int FreeMem(struct CommandLineArgsTag CLA);

#if TRACKMEMUSE
void printmemuse() {
   pid_t mypid=getpid();
   char commandline[256];
   fflush(NULL);
   sprintf(commandline,"cat /proc/%d/status | /bin/grep Vm | /usr/bin/fmt -140 -u", (int)mypid);
   system(commandline);
   fflush(NULL);
 }
#endif

/************************************* MAIN PROGRAM *************************************/

int main(int argc,char *argv[])
{
  int j;

  #if TRACKMEMUSE
    printf("Memory use at startup is:\n"); printmemuse();
  #endif

  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;
  SegmentDuration = CommandLineArgs.GPSEnd - CommandLineArgs.GPSStart ;

  /* create Frame cache, open frame stream and delete frame cache */
  LALFrCacheImport(&status,&framecache,CommandLineArgs.FrCacheFile);
  TESTSTATUS( &status );
  LALFrCacheOpen(&status,&framestream,framecache);
  TESTSTATUS( &status );
  LALDestroyFrCache(&status,&framecache);
  TESTSTATUS( &status );

  #if TRACKMEMUSE
    printf("Memory use after reading command line arguments and reading frame cache:\n"); printmemuse();
  #endif

  if( SegmentDuration < CommandLineArgs.T)
    {
      fprintf(stderr,"Cannot fit an SFT of duration %d between %d and %d\n",
	      CommandLineArgs.T, CommandLineArgs.GPSStart, CommandLineArgs.GPSEnd );
      return 0;;
    }

  gpsepoch.gpsSeconds = CommandLineArgs.GPSStart;
  gpsepoch.gpsNanoSeconds = 0;

  /* Allocates space for data */
  if (AllocateData(CommandLineArgs)) return 2;

  for(j=0; j < SegmentDuration/CommandLineArgs.T; j++) 
    {

      gpsepoch.gpsSeconds = CommandLineArgs.GPSStart + j*CommandLineArgs.T;
      gpsepoch.gpsNanoSeconds = 0;

      /* Reads T seconds of data */
      if (ReadData(CommandLineArgs)) return 3;

      /* High-pass data with Butterworth filter */
      if(HighPass(CommandLineArgs)) return 4;

      /* Window data */
      if(WindowData(CommandLineArgs)) return 5;

      /* create an SFT */
      if(CreateSFT(CommandLineArgs)) return 6;

      /* write out sft */
      if(WriteSFT(CommandLineArgs)) return 7;

    }

  if(FreeMem(CommandLineArgs)) return 8;

  #if TRACKMEMUSE
        printf("Memory use after freeing all allocated memory:\n"); printmemuse();
  #endif

  return 0;
}

/************************************* MAIN PROGRAM ENDS *************************************/

/*** FUNCTIONS ***/

/*******************************************************************************/
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA)
{
  INT4 errflg=0;
  struct option long_options[] = {
    {"high-pass-freq",       required_argument, NULL,          'f'},
    {"sft-duration",         required_argument, NULL,          't'},
    {"sft-write-path",       required_argument, NULL,          'p'},
    {"frame-cache",          required_argument, NULL,          'C'},
    {"channel-name",         required_argument, NULL,          'N'},
    {"gps-start-time",       required_argument, NULL,          's'},
    {"gps-end-time",         required_argument, NULL,          'e'},
    {"ht-data",              no_argument,       NULL,          'H'},
    {"use-single",           no_argument,       NULL,          'S'},    
    {"help",                 no_argument,       NULL,          'h'},
    {0, 0, 0, 0}
  };
  char args[] = "hHSf:t:C:N:s:e:p:";
  
  /* Initialize default values */
  CLA->HPf=-1.0;
  CLA->T=0.0;
  CLA->FrCacheFile=NULL;
  CLA->GPSStart=0;
  CLA->GPSEnd=0;
  CLA->ChannelName=NULL;
  CLA->SFTpath=NULL;
  CLA->htdata = 0;
  CLA->useSingle = 0; /* 11/19/05 gam; default is to use double precision, not single. */
  
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
    case 'H':
      /* high pass frequency */
      CLA->htdata=1;
      break;
    case 'S':
      /* use single precision */
      CLA->useSingle=1;
      break;
    case 'f':
      /* high pass frequency */
      CLA->HPf=atof(optarg);
      break;
    case 't':
      /* SFT time */
      CLA->T=atoi(optarg);
      break;
    case 'C':
      /* name of frame cache file */
      CLA->FrCacheFile=optarg;
      break;
    case 's':
      /* GPS start */
      CLA->GPSStart=atof(optarg);
      break;
    case 'e':
      /* GPS end */
      CLA->GPSEnd=atof(optarg);
      break;
    case 'N':
      CLA->ChannelName=optarg;       
      break;
    case 'p':
      CLA->SFTpath=optarg;       
      break;
    case 'h':
      /* print usage/help message */
      fprintf(stdout,"Arguments are:\n");
      fprintf(stdout,"\thigh-pass-freq (-f)\tFLOAT\t High pass filtering frequency in Hz.\n");
      fprintf(stdout,"\tsft-duration (-t)\tFLOAT\t SFT duration in seconds.\n");
      fprintf(stdout,"\tsft-write-path (-p)\tFLOAT\t Location of ouput SFTs.\n");
      fprintf(stdout,"\tframe-cache (-C)\tSTRING\t Path to frame cache file (including the filename).\n");
      fprintf(stdout,"\tgps-start-time (-s)\tINT\t GPS start time of segment.\n");
      fprintf(stdout,"\tgps-end-time (-e)\tINT\t GPS end time of segment.\n");
      fprintf(stdout,"\tchannel-name (-N)\tSTRING\t Name of channel to read within a frame.\n");
      fprintf(stdout,"\tht-data (-H)\t\tFLAG\t This is h(t) data (use Real8's; Procdata ).\n");
      fprintf(stdout,"\tuse-single (-S)\t\tFLAG\t Use single precision for window, plan, and fft; double precision filtering is always done.\n");
      fprintf(stdout,"\thelp (-h)\t\tFLAG\t This message\n");    
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

  if(CLA->HPf < 0 )
    {
      fprintf(stderr,"No high pass filtering frequency specified.\n"); 
      fprintf(stderr,"If you don't want to high pass filter set the frequency to 0.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }      
  if(CLA->T == 0.0)
    {
      fprintf(stderr,"No SFT duration specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }      
  if(CLA->GPSStart == 0)
    {
      fprintf(stderr,"No GPS start time specified.\n");
       fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->GPSEnd == 0)
    {
      fprintf(stderr,"No GPS end time specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->FrCacheFile == NULL)
    {
      fprintf(stderr,"No frame cache file specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }      
  if(CLA->ChannelName == NULL)
    {
      fprintf(stderr,"No data channel name specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }      
  if(CLA->SFTpath == NULL)
    {
      fprintf(stderr,"No output path specified for SFTs.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }      

  return errflg;
}
/*******************************************************************************/

/*******************************************************************************/
int AllocateData(struct CommandLineArgsTag CLA)
{
  static FrChanIn chanin;

  chanin.name  = CLA.ChannelName;

  /* These calls just return deltaT for the channel */
  if(CLA.htdata)
    {
      chanin.type  = ProcDataChannel;
      /* Get channel time step size by calling LALFrGetREAL8TimeSeries */
      LALFrSeek(&status,&gpsepoch,framestream);
      TESTSTATUS( &status );
      LALFrGetREAL8TimeSeries(&status,&dataDouble,&chanin,framestream);
      TESTSTATUS( &status );
      dataSingle.deltaT=dataDouble.deltaT;      
    }
  else
    {
      chanin.type  = ADCDataChannel;
      /* Get channel time step size by calling LALFrGetREAL4TimeSeries */
      LALFrSeek(&status,&gpsepoch,framestream);
      TESTSTATUS( &status );
      LALFrGetREAL4TimeSeries(&status,&dataSingle,&chanin,framestream);
      TESTSTATUS( &status );
      dataDouble.deltaT=dataSingle.deltaT;
    }

  /*  11/19/05 gam; will keep dataDouble or dataSingle in memory at all times, depending on whether using single or double precision */
  if(CLA.useSingle) {
      LALCreateVector(&status,&dataSingle.data,(UINT4)(CLA.T/dataSingle.deltaT +0.5));
      TESTSTATUS( &status );

      #if TRACKMEMUSE
         printf("Memory use after creating dataSingle and before calling LALCreateForwardRealFFTPlan.\n"); printmemuse();
      #endif

      LALCreateForwardRealFFTPlan( &status, &fftPlanSingle, dataSingle.data->length, 0 );
      TESTSTATUS( &status );

      #if TRACKMEMUSE
         printf("Memory use after creating dataSingle and after calling LALCreateForwardRealFFTPlan.\n"); printmemuse();
      #endif

  } else {  
      LALDCreateVector(&status,&dataDouble.data,(UINT4)(CLA.T/dataDouble.deltaT +0.5));
      TESTSTATUS( &status );

      #if TRACKMEMUSE
         printf("Memory use after creating dataDouble and before calling LALCreateForwardREAL8FFTPlan.\n"); printmemuse();
      #endif
      
      LALCreateForwardREAL8FFTPlan( &status, &fftPlanDouble, dataDouble.data->length, 0 );
      TESTSTATUS( &status );

      #if TRACKMEMUSE
            printf("Memory use after creating dataDouble and after calling LALCreateForwardREAL8FFTPlan.\n"); printmemuse();
      #endif
  }

  return 0;
}
/*******************************************************************************/

/*******************************************************************************/
int ReadData(struct CommandLineArgsTag CLA)
{
  static FrChanIn chanin;
  INT4 k;
  chanin.name  = CLA.ChannelName;
  LALFrSeek(&status,&gpsepoch,framestream);
  TESTSTATUS( &status );

  /* 11/19/05 gam */
  if(CLA.useSingle) {
    
    if(CLA.htdata)
    {

      #if TRACKMEMUSE
        printf("Memory use before creating dataDouble and calling LALFrGetREAL8TimeSeries.\n"); printmemuse();
      #endif

      LALDCreateVector(&status,&dataDouble.data,(UINT4)(CLA.T/dataDouble.deltaT +0.5));
      TESTSTATUS( &status );

      chanin.type  = ProcDataChannel;
      LALFrGetREAL8TimeSeries(&status,&dataDouble,&chanin,framestream);
      TESTSTATUS( &status );

      #if TRACKMEMUSE
        printf("Memory use after creating dataDouble and calling LALFrGetREAL8TimeSeries.\n"); printmemuse();
      #endif

      /*copy the data into the single precision timeseries */
      for (k = 0; k < dataSingle.data->length; k++) {
          dataSingle.data->data[k] = dataDouble.data->data[k];
      }

      LALDDestroyVector(&status,&dataDouble.data);
      TESTSTATUS( &status );

      #if TRACKMEMUSE
        printf("Memory use after destroying dataDouble.\n"); printmemuse();
      #endif

    }
    else
    {

      #if TRACKMEMUSE
        printf("Memory use before calling LALFrGetREAL4TimeSeries.\n"); printmemuse();
      #endif

      chanin.type  = ADCDataChannel;
      LALFrGetREAL4TimeSeries(&status,&dataSingle,&chanin,framestream);
      TESTSTATUS( &status );

      #if TRACKMEMUSE
        printf("Memory use after calling LALFrGetREAL4TimeSeries.\n"); printmemuse();
      #endif

    }

  } else {

    if(CLA.htdata)
    {
    
      #if TRACKMEMUSE
        printf("Memory use before calling LALFrGetREAL8TimeSeries.\n"); printmemuse();
      #endif
    
      chanin.type  = ProcDataChannel;
      LALFrGetREAL8TimeSeries(&status,&dataDouble,&chanin,framestream);
      TESTSTATUS( &status );

      #if TRACKMEMUSE
        printf("Memory use after calling LALFrGetREAL8TimeSeries.\n"); printmemuse();
      #endif

    }
    else
    {
      #if TRACKMEMUSE
        printf("Memory use before creating dataSingle and calling LALFrGetREAL4TimeSeries.\n"); printmemuse();
      #endif

      /* 11/02/05 gam; add next two lines */
      LALCreateVector(&status,&dataSingle.data,(UINT4)(CLA.T/dataSingle.deltaT +0.5));
      TESTSTATUS( &status );
      
      chanin.type  = ADCDataChannel;
      LALFrGetREAL4TimeSeries(&status,&dataSingle,&chanin,framestream);
      TESTSTATUS( &status );
 
      #if TRACKMEMUSE
        printf("Memory use after creating dataSingle and calling LALFrGetREAL4TimeSeries.\n"); printmemuse();
      #endif

      /*copy the data into the double precision timeseries */
      for (k = 0; k < (int)dataDouble.data->length; k++) {
          dataDouble.data->data[k] = dataSingle.data->data[k];
      }     

      LALDestroyVector(&status,&dataSingle.data);
      TESTSTATUS( &status );

      #if TRACKMEMUSE
        printf("Memory use after destroying dataSingle.\n"); printmemuse();
      #endif

    }

  }

  return 0;
}
/*******************************************************************************/

/*******************************************************************************/
int HighPass(struct CommandLineArgsTag CLA)
{
  PassBandParamStruc filterpar;

  filterpar.name  = "Butterworth High Pass";
  filterpar.nMax  = 10;
  filterpar.f2    = CLA.HPf;
  filterpar.a2    = 0.5;
  filterpar.f1    = -1.0;
  filterpar.a1    = -1.0;

  if (CLA.HPf > 0.0 )
    {
      if(CLA.useSingle) {
        #if TRACKMEMUSE
          printf("Memory use before calling LALDButterworthREAL4TimeSeries:\n"); printmemuse();
        #endif

        /* High pass the time series */
        LALDButterworthREAL4TimeSeries(&status,&dataSingle,&filterpar);
        TESTSTATUS( &status );

        #if TRACKMEMUSE
          printf("Memory use after calling LALDButterworthREAL4TimeSeries:\n"); printmemuse();
        #endif      
      } else {
        #if TRACKMEMUSE
          printf("Memory use before calling LALButterworthREAL8TimeSeries:\n"); printmemuse();
        #endif

        /* High pass the time series */
        LALButterworthREAL8TimeSeries(&status,&dataDouble,&filterpar);
        TESTSTATUS( &status );

        #if TRACKMEMUSE
          printf("Memory use after calling LALButterworthREAL8TimeSeries:\n"); printmemuse();
        #endif
      }
    }

  return 0;
}
/*******************************************************************************/

/*******************************************************************************/
int WindowData(struct CommandLineArgsTag CLA)
{
  REAL8 r = 0.001;
  INT4 k, N, kl, kh;
  /* This implementation of a Tukey window is describes 
     in the Matlab tukey window documentation */

  /* 11/19/05 gam */
  if(CLA.useSingle) {
    N=dataSingle.data->length;
    kl=r/2*(N-1)+1;
    kh=N-r/2*(N-1)+1;
    for(k = 1; k < kl; k++) 
    {
      dataSingle.data->data[k-1]=dataSingle.data->data[k-1]*0.5*( (REAL4)( 1.0 + cos(LAL_TWOPI/r*(k-1)/(N-1) - LAL_PI) ) );
    }
    for(k = kh; k <= N; k++) 
    {
      dataSingle.data->data[k-1]=dataSingle.data->data[k-1]*0.5*( (REAL4)( 1.0 + cos(LAL_TWOPI/r - LAL_TWOPI/r*(k-1)/(N-1) - LAL_PI) ) );
    }    
  } else {
    N=dataDouble.data->length;
    kl=r/2*(N-1)+1;
    kh=N-r/2*(N-1)+1;
    for(k = 1; k < kl; k++) 
    {
      dataDouble.data->data[k-1]=dataDouble.data->data[k-1]*0.5*( 1.0 + cos(LAL_TWOPI/r*(k-1)/(N-1) - LAL_PI) );
    }
    for(k = kh; k <= N; k++) 
    {
      dataDouble.data->data[k-1]=dataDouble.data->data[k-1]*0.5*( 1.0 + cos(LAL_TWOPI/r - LAL_TWOPI/r*(k-1)/(N-1) - LAL_PI) );
    }
  }

  return 0;
}
/*******************************************************************************/

/*******************************************************************************/
int CreateSFT(struct CommandLineArgsTag CLA)
{

  /* 11/19/05 gam */
  if(CLA.useSingle) {
      #if TRACKMEMUSE
        printf("Memory use before creating output vector fftDataSingle and calling LALForwardRealFFT:\n"); printmemuse();
      #endif

      /* 11/02/05 gam; allocate container for SFT data */
      LALCCreateVector( &status, &fftDataSingle, dataSingle.data->length / 2 + 1 );
      TESTSTATUS( &status );  

      /* compute sft */
      LALForwardRealFFT(&status, fftDataSingle, dataSingle.data, fftPlanSingle );
      TESTSTATUS( &status );      

      #if TRACKMEMUSE
        printf("Memory use after creating output vector fftDataSingle and calling LALForwardRealFFT:\n"); printmemuse();
      #endif  
  } else {
      #if TRACKMEMUSE
        printf("Memory use before creating output vector fftDataDouble and calling XLALREAL8ForwardFFT:\n"); printmemuse();
      #endif

      /* 11/02/05 gam; allocate container for SFT data */
      LALZCreateVector( &status, &fftDataDouble, dataDouble.data->length / 2 + 1 );
      TESTSTATUS( &status );  

      /* compute sft */
      XLALREAL8ForwardFFT( fftDataDouble, dataDouble.data, fftPlanDouble );

      #if TRACKMEMUSE
        printf("Memory use after creating output vector fftDataDouble and calling XLALREAL8ForwardFFT:\n"); printmemuse();
      #endif
  }
      return 0;
}
/*******************************************************************************/

/*******************************************************************************/
int WriteSFT(struct CommandLineArgsTag CLA)
{
  char sftname[256];
  char ifo[2];
  int firstbin=(INT4)(FMIN*CLA.T+0.5), k;
  FILE *fpsft;
  int errorcode1=0,errorcode2=0;
  char gpstime[10];
  REAL4 rpw,ipw;
  
  strncpy( ifo, CLA.ChannelName, 2 );
  strcpy( sftname, CLA.SFTpath );

  strcat(sftname,"/SFT_");
  strncat(sftname,ifo,2);
  sprintf(gpstime,".%09d",gpsepoch.gpsSeconds);
  strcat(sftname,gpstime);

  /* open SFT file for writing */
  fpsft=tryopen(sftname,"w");
      
  /* write header */
  header.endian=1.0;
  header.gps_sec=gpsepoch.gpsSeconds;
  header.gps_nsec=gpsepoch.gpsNanoSeconds;
  header.tbase=CLA.T;
  header.firstfreqindex=firstbin;
  header.nsamples=(INT4)(DF*CLA.T+0.5);
  if (1!=fwrite((void*)&header,sizeof(header),1,fpsft)){
    fprintf(stderr,"Error in writing header into file %s!\n",sftname);
    return 7;
  }

  /* 11/19/05 gam */
  if(CLA.useSingle) {
    /* Write SFT */
    for (k=0; k<header.nsamples; k++)
    {
      rpw=((REAL4)(((REAL8)DF)/(0.5*(REAL8)(1/dataSingle.deltaT))))
	* fftDataSingle->data[k+firstbin].re;
      ipw=((REAL4)(((REAL8)DF)/(0.5*(REAL8)(1/dataSingle.deltaT))))
	* fftDataSingle->data[k+firstbin].im;
      errorcode1=fwrite((void*)&rpw, sizeof(REAL4),1,fpsft);
      errorcode2=fwrite((void*)&ipw, sizeof(REAL4),1,fpsft);
    } 
    LALCDestroyVector( &status, &fftDataSingle );
    TESTSTATUS( &status );    
  } else {    
    /* Write SFT */
    for (k=0; k<header.nsamples; k++)
    {
      rpw=(((REAL8)DF)/(0.5*(REAL8)(1/dataDouble.deltaT))) 
	* fftDataDouble->data[k+firstbin].re;
      ipw=(((REAL8)DF)/(0.5*(REAL8)(1/dataDouble.deltaT))) 
	* fftDataDouble->data[k+firstbin].im;
      errorcode1=fwrite((void*)&rpw, sizeof(REAL4),1,fpsft);
      errorcode2=fwrite((void*)&ipw, sizeof(REAL4),1,fpsft);
    }
    /* 11/02/05 gam; deallocate container for SFT data */
    LALZDestroyVector( &status, &fftDataDouble );
    TESTSTATUS( &status );    
  }

  /* Check that there were no errors while writing SFTS */
  if (errorcode1-1 || errorcode2-1)
    {
      fprintf(stderr,"Error in writing data into SFT file %s!\n",sftname);
      return 8;
    }
  
  fclose(fpsft);


  return 0;
}
/*******************************************************************************/

/*******************************************************************************/
int FreeMem(struct CommandLineArgsTag CLA)
{

  LALFrClose(&status,&framestream);
  TESTSTATUS( &status );

  /* 11/19/05 gam */
  if(CLA.useSingle) {
    LALDestroyVector(&status,&dataSingle.data);
    TESTSTATUS( &status );
    LALDestroyRealFFTPlan( &status, &fftPlanSingle );
    TESTSTATUS( &status );
  } else {
    LALDDestroyVector(&status,&dataDouble.data);
    TESTSTATUS( &status );
    LALDestroyREAL8FFTPlan( &status, &fftPlanDouble );
    TESTSTATUS( &status );    
  }

  LALCheckMemoryLeaks();
 
  return 0;
}
/*******************************************************************************/

#endif
