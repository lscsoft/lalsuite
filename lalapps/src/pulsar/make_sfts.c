#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdarg.h>
#include <sys/types.h>
#include <errno.h>

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
#include <lal/RealFFT.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/PrintVector.h>

/* Frame headers */
#include <FrameL.h>

/* Define the parameters to make the window */
#define WINSTART 150
#define WINEND   200
#define WINLEN (WINEND-WINSTART)
/* with i ranging from 0 to WINLEN... */
#define WINFUN(i) ((sin(i*LAL_PI/(WINLEN)-LAL_PI_2)+1.0)/2.0)

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 100; } else ((void)0)

/* IFO sample rate in Hz */
#define SRATE 16384

/* debug level for LAL */
INT4 lalDebugLevel = LALERROR | LALWARNING | LALINFO;

/* Header that defines the GEO SFT Data Format */
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
#define FMIN  (48.0)
#define DF    (2000.0)

/* maximum number of start times from input file */
#define MAXSTART 1000

/* This will hold the SFT and PLAN */
COMPLEX8Vector *fvec = NULL;
RealFFTPlan *pfwd = NULL;

/* Need prototype: this is POSIX not ANSI, sigh */
#include <unistd.h>
/* int gethostname(char *name, size_t len); */

/* This is an error handler that prints a bit of extra info */
void pout(char *fmt, ...)  
     __attribute__ ((format (printf, 1, 2)));

void pout(char *fmt, ...){
  va_list ap;
  char hostsname[256];
  pid_t processid;

  processid=getpid();

  /* not guaranteed null terminated! */
  gethostname(hostsname, 255);
  hostsname[255]='\0';

  /* printout process id and hostname into error message */
  fprintf(stderr,"%s [PID=%d] ", hostsname, processid);

  /* initialize variable argument list */
  va_start(ap,fmt);

  /* print out, and flush */
  vfprintf(stderr, fmt,ap);
  va_end(ap);
  fflush(stderr);
  return;
}


/* This repeatedly tries to re-open a file.  This is useful on a
   cluster because automount may fail or time out.  This allows a
   number of tries with a bit of rest in between.
*/
FILE* tryopen(char *name, char *mode){
  int count=0;
  FILE *fp;
  
  while (!(fp=fopen(name, mode))){
    pout("Unable to open file %s in mode %s.  Will retry...\n", name, mode);
    if (count++<10)
      sleep(10);
    else
      exit(3);
  }
  return fp;
}

int getenvval(const char *valuename){
  int retval;
  
  char *value=getenv(valuename);
  if (value)
    retval=atoi(value);
  else
    retval=0;

  printf("Environment variable %s = %d\n", valuename, retval);

  return retval;
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

  /* Control variable for run-time options */
  int doubledata=getenvval("DOUBLEDATA");
  int windowdata=getenvval("WINDOWDATA");
  int printfreqseries=getenvval("PRINTFREQSERIES");

  /* FrameL frame file */
  FrFile *frfile;
  /* vector holding the frame data */
  FrVect *frvect;        
  
  printf("Normal startup\n");
  fflush(stdout);
  
  /* check command syntax */
  if (argc !=5 || (jobnum=atoi(argv[1]))<0 || jobnum>99999){
    int a;
    pout("Syntax:\n\t%s N DIR1 DIR2 DETECTOR\nwhere 0<=N<=99999.\n",argv[0]);
    pout("Files used are jobdata.N, jobtimes.N, where N has five digits\n");
    pout("Input files are taken from DIR1\n");
    pout("Output files are put into DIR2\n");
    pout("DETECTOR is either H1, H2, or L1\n\n");
    pout("There were argc=%d arguments:\n",argc);
    for (a=0;a<argc;a++)
      pout("arg=%d: %s\n",a,argv[a]);
    return 2;
  }
  
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
	pout("File %s contains inconsistent SFT time baselines %d != %d\n",timelist,tbasei,tbase);
	return 3;
      }
      
      /* increment counter of lines read */
      nstarts++;
      
      /* and make sure that we are not yet maxed out */
      if (nstarts>=MAXSTART){
	pout("More than MAXSTART=%d lines in file: %s.  Increase MAXSTART and recompile\n",
		MAXSTART,timelist);
	return 3;
      }
    }
    
    /* check that file contained at least some valid data */
    if (!nstarts){
      pout("File %s didn't contain any valid lines!\n",timelist);
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
    pout( "Couldnt open frame file list %s\n", framelist);
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
  
  /* Main loop to write SFTs */
  for (count=0;count<nstarts;count++) {
    struct stat buff;
    char sftname[256];
    int filesize=(INT4)(DF*tbase);

    /* time of correct start */
    epoch.gpsSeconds=starts[count];
    epoch.gpsNanoSeconds=0;
    
    /* to check that an existing file has the correct size */
    if (doubledata)
      filesize *= 2*sizeof(double);
    else
      filesize *= 2*sizeof(float);
    filesize+=sizeof(header);
    
    /* construct SFT name.  If file exists, and has correct size, just continue */
    sprintf(sftname,"%s/SFT_%s.%d",argv[3],argv[4],epoch.gpsSeconds);
    if (!stat(sftname, &buff) && buff.st_size==filesize){
      printf("File: %s already exists with correct size; skipping it.\n",sftname);
      fflush(stdout);
      continue;
    }
    
    /* read in correct data */
    errno=0;
    frvect = FrFileIGetVAdc(frfile, chname, epoch.gpsSeconds, tbase, 0);
    
    /* This block will detect permission denied and other access errors */
    if (errno){
      char command[256];
      pout("System error when reading Frame data: %s\n", strerror(errno));
      pout("Following is output of ls -l on Frame file list:\n");
      sprintf(command, "cat %s | xargs ls -l 1>&2", framelist);
      system(command);
      return 3;
    }
    
    /* This block detects the case where the frame file list did not contain all data */
    if (frvect==NULL) {
      pout( "Data missing between times %d and %d\n",epoch.gpsSeconds,epoch.gpsSeconds+tbase);
      continue;
    }
    
    /* Check that data type is correct */
    if (frvect->type != FR_VECT_4R){
      pout( "Wrong data type (notFR_VECT_4R) found in frame!\n" );
      return(5);
    }
    
    /* check for gaps */
    if (frvect->next){
      pout( "Data between times %d and %d had a gap\n",epoch.gpsSeconds,epoch.gpsSeconds+tbase);
      /* free space allocated by frame library */
      FrVectFree(frvect);
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
      
      /* apply high-pass Butterworth filter */
      LALButterworthREAL4TimeSeries(&status, &chan, &filterpar);
      TESTSTATUS(&status);
      
      /* Turn on windowing */
      if (windowdata){
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
      }
      
      /* open SFT file for writing */
      fpsft=tryopen(sftname,"w");
      
      /* write header */
      header.endian=1.0;
      header.gps_sec=epoch.gpsSeconds;
      header.gps_nsec=epoch.gpsNanoSeconds;
      header.tbase=tbase;
      header.firstfreqindex=firstbin;
      header.nsamples=(INT4)(DF*tbase);
      if (1!=fwrite((void*)&header,sizeof(header),1,fpsft)){
	pout("Error in writing header into file %s!\n",sftname);
	return 7;
      }
      
      /* take forward FFT */
      LALForwardRealFFT(&status, fvec, chan.data, pfwd);
      TESTSTATUS(&status);
      
      if (printfreqseries){
	/* for debugging -- print freq series */
	LALCPrintVector(fvec);
	exit(0);
      }   
      
      /* Write SFT */
      for (k=0; k<header.nsamples; k++){
	int errorcode1,errorcode2;

	if (doubledata){
	  REAL8 rpw=fvec->data[k+firstbin].re;
	  REAL8 ipw=fvec->data[k+firstbin].im;
	  errorcode1=fwrite((void*)&rpw, sizeof(REAL8),1,fpsft);
	  errorcode2=fwrite((void*)&ipw, sizeof(REAL8),1,fpsft);  
	}
	else {
	  REAL4 rpw=fvec->data[k+firstbin].re;
	  REAL4 ipw=fvec->data[k+firstbin].im;
	  errorcode1=fwrite((void*)&rpw, sizeof(REAL4),1,fpsft);
	  errorcode2=fwrite((void*)&ipw, sizeof(REAL4),1,fpsft);
	}
	
	/* Check that there were no errors while writing SFTS */
 	if (errorcode1-1 || errorcode2-1){
	  pout("Error in writing data into SFT file %s!\n",sftname);
	  return 8;
	}
      }
      fclose(fpsft);
    }
  }
  

  /* Clean up and exit */
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
 
  return 0;
}
