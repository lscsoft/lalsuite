/*#include <FrameL.h> */
#include <stdio.h>
#include <unistd.h>
#include <math.h>

/* INT_MAX */
#include <limits.h>

/* Does a file exist */
#include <sys/types.h>
#include <sys/stat.h>

/* LAL stuff */
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/PrintFTSeries.h>
#include <lal/LALConstants.h>

/*#include <lal/FrameStream.h> */
#include <lal/RealFFT.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/PrintVector.h>

/* Frame headers */
#include <FrameL.h>
#ifndef FRLONG
#include "FrIO.h" /* corrected Version of FrIO.h */
#endif

typedef struct tagFrFileInfo {
  INT4  ind;
  CHAR *url;
  INT4  t0;
  INT4  dt;
} FrFileInfo;

/* Definition of FrStream */
struct tagFrStream {
  FrFileInfo     *filelist;
  UINT4           numfiles;
  UINT4           filenum;
  struct FrFile  *frfile;
  struct FrameH  *frame;
  LIGOTimeGPS     epoch;
  INT4            end;
  INT4            err;
  INT4            gap;
};

#define WINSTART 150
#define WINEND   200
#define WINLEN (WINEND-WINSTART)
/* with i ranging from 0 to WINLEN... */
#define WINFUN(i) ((sin(i*LAL_PI/(WINLEN)-LAL_PI_2)+1.0)/2.0)

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 666; } else ((void)0)

/* IFO sample rate in Hz */
#define SRATE 16384

/* If DOUBLEDATA not defined, use REAL4 in SFT files */
/* #define DOUBLEDATA */

INT4 lalDebugLevel = LALERROR | LALWARNING | LALINFO /* | LALTRACE */ ;

struct headertag {
  REAL8 endian;
  INT4  gps_sec;
  INT4  gps_nsec;
  REAL8 tbase;
  INT4  firstfreqindex;
  INT4  nsamples;
} header;

/* Minimum frequency (Hz, float) bandwidth (Hz, float) and baseline
   time (seconds,k int) for SFTs */
#define FMIN  (50.0)
#define DF    (2000.0)

/* maximum number of start times from input file */
#define MAXSTART 1000

/* This will hold the SFT and PLAN */
COMPLEX8Vector *fvec = NULL;
RealFFTPlan *pfwd = NULL;

/* This repeatedly tries to re-open a file.  This is useful on a
   cluster because it's possible that an automount may fail or time
   out.  This allows a number of tries with a bit of rest in between.
*/
FILE* tryopen(char *name, char *mode){
  int count=0;
  FILE *fp;

  while (!(fp=fopen(name, mode))){
    fprintf(stderr,"Unable to open file %s in mode %s\n", name, mode);
    fflush(stderr);
    if (count++<10)
      sleep(10);
    else
      exit(3);
  }
  return fp;
}



int main(int argc,char *argv[]){
  static LALStatus status;
  UINT4 npts;
  static REAL4TimeSeries chan;
  LIGOTimeGPS epoch;
  INT4 jobnum;
  INT4 starts[MAXSTART],nstarts=0,count=0;
  INT4 i,len2=0;
  char chname[256]; 
  INT4 tbase=INT_MAX,firstbin;
  PassBandParamStruc filterpar;
  REAL8 window[WINLEN];
  char framelist[256];
  int opencount=0;
  
/* FrameL frame file */
  FrFile *frfile;
 /* vector holding the frame data */
  FrVect *frvect;        

  printf("Normal startup\n");
  fflush(stdout);
  fprintf(stderr,"Normal startup\n");
  fflush(stderr);

  /* Check frame library version */
#if (0)
    if ( FRAMELIB_VERSION < 4.53 ) {
    fprintf(stderr,"FRAMELIB_VERSION=%f too early\n",FRAMELIB_VERSION);
    fflush(stderr);
    return 1;
    }
#endif 

    /* check command syntax.  If code compiled without CONDOR enabled,
       then it reads arguments from the command line.  If it's
       compiled with CONDOR defined, then it just takes a single
       jobnumber from stdin, and other inputs from command line.
    */
#ifndef CONDOR
    if (argc !=5 || (jobnum=atoi(argv[1]))<0 || jobnum>99999){
   int a;
    fprintf(stderr,"Syntax:\n\t%s N DIR1 DIR2 DETECTOR\nwhere 0<=N<=99999.\n",argv[0]);
    fprintf(stderr,"Files used are jobdata.N, jobtimes.N, where N has five digits\n");
    fprintf(stderr,"Input files are taken from DIR1\n");
    fprintf(stderr,"Output files are put into DIR2\n");
    fprintf(stderr,"DETECTOR is either H1, H2, or L1\n\n");
    fprintf(stderr,"There were argc=%d arguments:\n",argc);
    fflush(stderr);
    for (a=0;a<argc;a++)
	fprintf(stderr,"arg=%d: %s\n",a,argv[a]);
    fflush(stderr);
    return 2;
  }
#else
  /* condor versions get their jobnumber from stdin. All other */
  /* arguments are constant for a given job. */
  scanf("%d",&jobnum);
#endif 

  /* construct channel name */
  sprintf(chname,"%s:LSC-AS_Q",argv[4]);
  
  { /* read start times values and how many of them there are */
    char timelist[256];
    INT4 tbasei;
    FILE *stfp=NULL;
    
    /* construct name of file containing segment start times */
    sprintf(timelist,"%s/jobtimes.%05d",argv[2],jobnum);
    stfp=tryopen(timelist,"r");
    
    /* read file into an array */
    while (2==fscanf(stfp,"%d %d",starts+nstarts,&tbasei)){
      if (nstarts==0)
	tbase=tbasei;
      else if (tbasei!=tbase){
	fprintf(stderr,"File %s contains inconsistent SFT time baselines %d != %d\n",timelist,tbasei,tbase);
	fflush(stderr);
	return 3;
      }

      /* increment counter of lines read */
      nstarts++;

      /* and make sure that we are not yet maxed out */
      if (nstarts>=MAXSTART){
	fprintf(stderr,"More than MAXSTART=%d lines in file: %s.  Increase MAXSTART and recompile\n",
		MAXSTART,timelist);
	fflush(stderr);
	return 3;
      }
    }

    /* check that file contained at least some valid data */
    if (!nstarts){
      fprintf(stderr,"File %s didn't contain any valid lines!\n",timelist);
      fflush(stderr);
      return 3;
    }
    fclose(stfp);
  }
  
  /* Compute things that require the time baseline */
  firstbin=(INT4)(FMIN*tbase);
  npts = SRATE*tbase;
  len2=npts/2+1;

  /* init window function */
  for(i=0; i<WINLEN; i++)
    window[i]=WINFUN(i);

  /* set filtering parameters */
  filterpar.name = "Butterworth High Pass";
  filterpar.nMax  = 5;
  filterpar.f2 = 100.0;
  filterpar.a2 = 0.5;
  /* values that are 'not given' = out of range */
  filterpar.f1 = 0.0;
  filterpar.a1 = 0.0;

  /* Initialize frame library with correct file list */
  sprintf(framelist,"%s/jobdata.%05d.ffl",argv[2],jobnum);
  opencount=0;
  while (!(frfile = FrFileINew(framelist))){
    fprintf(stderr, "Couldnt open frame file %s\n", framelist);
    fflush(stderr);
    if (opencount++<10)
      sleep(10);
    else
      return(3);
  }
 
  /* create structure to store channel data  */
  chan.data = NULL;
  LALSCreateVector(&status, &chan.data, npts);
  TESTSTATUS(&status);

  /* Create vector to hold signal frequency series */
  LALCCreateVector(&status, &fvec, (UINT4)len2);
  TESTSTATUS(&status);
  
  /* Compute measured plan for FFTW */
  LALCreateForwardRealFFTPlan(&status, &pfwd, (UINT4)npts, 0);
  TESTSTATUS(&status);

  /* printf("Starting main loop to write SFTs\n"); */
  /* now start reading in data */
  for (count=0;count<nstarts;count++) {
    /*    LALTYPECODE typecode; */
    struct stat buff;
    char sftname[256];
    int filesize=(INT4)(DF*tbase);
    
    /* time of correct start */
    epoch.gpsSeconds=starts[count];
    epoch.gpsNanoSeconds=0;

    /* to check that an existing file has the correct size.... */
#ifdef DOUBLEDATA
    filesize *= 2*sizeof(double);
#else
    filesize *= 2*sizeof(float);
#endif
    filesize+=sizeof(header);  

    /* construct SFT name.  If file exists, and has correct size, just continue */
    sprintf(sftname,"%s/SFT_%s.%d",argv[3],argv[4],epoch.gpsSeconds);
    if (!stat(sftname, &buff) && buff.st_size==filesize){
      printf("File: %s already exists with correct size; skipping it.\n",sftname);
      fflush(stdout);
      continue;
    }
    
    /* read in correct data */
    frvect = FrFileIGetVAdc(frfile, chname, epoch.gpsSeconds, tbase, 0);
    if (frvect == NULL) {
      fprintf(stderr, "Data between times %d and %d not present\n",epoch.gpsSeconds,epoch.gpsSeconds+tbase);
      fflush(stderr);
      continue;
    }

    /* Check that data type is correct */
    if ( frvect->type != FR_VECT_4R ){
      fprintf(stderr, "Wrong data type (notFR_VECT_4R) found in frame!\n" );
      fflush(stderr);
      return(5);
    }

    /* check for gaps */
    if (frvect->next){
      fprintf(stderr, "Data between times %d and %d had a gap\n",epoch.gpsSeconds,epoch.gpsSeconds+tbase);
      fflush(stderr);
      FrVectFree(frvect); /* free space allocated by frame library */
      frvect=NULL;
      continue;
    }
    
    {
      /* We found all needed data -- output an SFT! */
      FILE *fpsft;
      unsigned int i;
      int k;

      /* set up LAL time series object sample interval, time epoch */
      strcpy(chan.name, chname);
      chan.deltaT = 1.0/SRATE;
      chan.epoch = epoch;
      chan.data->length = npts;

      /* copy data */
      for (i=0;i<npts;i++)
	chan.data->data[i]=frvect->dataF[i];

      /* free framevec -- no longer needed */
      FrVectFree(frvect); 
      frvect=NULL;
      
      /* apply butterworth filter */
      LALButterworthREAL4TimeSeries(&status, &chan, &filterpar);
      TESTSTATUS(&status);

      /* Turn on windowing */
#if (1)
      /* window data.  Off portion */
      for (i=0; i<WINSTART; i++) {
	chan.data->data[i] = 0.0;
	chan.data->data[chan.data->length - 1 - i] = 0.0;
      }
      /* window data, smooth turn-on portion */
      for (i=WINSTART; i<WINEND; i++) {
	double win=window[i - WINSTART];
	chan.data->data[i] *= win;
	chan.data->data[chan.data->length - 1 - i]  *= win;
      }
#endif

      /* open SFT file for writing */
      fpsft=tryopen(sftname,"w");
      
      /* start writing header */
      header.endian=1.0;
      header.gps_sec=epoch.gpsSeconds;
      header.gps_nsec=epoch.gpsNanoSeconds;
      header.tbase=tbase;
      header.firstfreqindex=firstbin;
      header.nsamples=(INT4)(DF*tbase);
      if (1!=fwrite((void*)&header,sizeof(header),1,fpsft)){
	fprintf(stderr,"Error in writing header into file %s!\n",sftname);
	fflush(stderr);
	return 7;
      }
      
      /* take forward FFT */
      LALForwardRealFFT(&status, fvec, chan.data, pfwd);
      TESTSTATUS(&status);
      
#ifdef PRINTFREQSERIES
      /* for debugging -- print freq series */
      LALCPrintVector(fvec);
      exit(0);
#endif    

      /* Write SFT */
      for (k=0; k<header.nsamples; k++){
#ifdef DOUBLEDATA
	REAL8 rpw=fvec->data[k+firstbin].re;
	REAL8 ipw=fvec->data[k+firstbin].im;
	INT4 errorcode1=fwrite((void*)&rpw, sizeof(REAL8),1,fpsft);
	INT4 errorcode2=fwrite((void*)&ipw, sizeof(REAL8),1,fpsft);  
#else 
	REAL4 rpw=fvec->data[k+firstbin].re;
	REAL4 ipw=fvec->data[k+firstbin].im;
	INT4 errorcode1=fwrite((void*)&rpw, sizeof(REAL4),1,fpsft);
  	INT4 errorcode2=fwrite((void*)&ipw, sizeof(REAL4),1,fpsft);
#endif
 	if (errorcode1-1 || errorcode2-1){
	  fprintf(stderr,"Error in writing data into SFT file %s!\n",sftname);
	  fflush(stderr);
	  return 8;
	}
      }
      fclose(fpsft);
    }
  }
  

  /* Clean up and exit */
  /* printf("Starting memory clean-up\n"); */

  FrFileIEnd(frfile);

  LALSDestroyVector( &status, &chan.data );
  TESTSTATUS( &status );
  
  LALCDestroyVector(&status, &fvec);
  TESTSTATUS( &status );
  
  LALDestroyRealFFTPlan(&status, &pfwd);
  TESTSTATUS( &status );
  
  LALCheckMemoryLeaks();

  printf("Normal exit\n");
  fflush(stdout);

  fprintf(stderr,"Normal exit\n");
  fflush(stderr);
  
  return 0;
}
