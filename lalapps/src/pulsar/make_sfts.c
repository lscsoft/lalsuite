#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdarg.h>
#include <sys/types.h>
#include <errno.h>
#include <signal.h>

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
#define WINSTART 4096
#define WINEND   8192
#define WINLEN (WINEND-WINSTART)
/* with i ranging from 0 to WINLEN... */
#define WINFUN(i) ((sin(i*LAL_PI/(WINLEN)-LAL_PI_2)+1.0)/2.0)

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 100; } else ((void)0)

/* IFO sample rate in Hz */
#define SRATE 16384

/* should we print the first 50 elements of various series */
/* #define PRINT50 */
#define PRINT50 0

/* should we use doubles for the time-domain filtering */
#define TDDOUBLE 1

/* should we shift timestamps to compensate for timing errors */
#define TIMESHIFT  0

/* should we use the new 2*DF/SRATE normalization convention? */
#define NEWNORM 1

/* provide data to Teviet for working on LAL filtering */
#define PRINTTEV 0

/* check of changing knee frequency */
#define HIGHFREQ 0

/* check of insensitivity to LSB of gravity-wave channel */
#define QUANTIZATION_TEST 0

/* track memory usage under linux */
#define TRACKMEMUSE 0

/* make calibrated SFTs starting with calibrated data */
#define USETDCALDATA 1

/* Use high pass filtering.  This *should* NOT be needed when working
   in the time domain with S3 data, which is already aggresively
   high-pass filtered at 40 Hz.
*/
#define HIGHPASS_FILTER 1


/* debug level for LAL */
INT4 lalDebugLevel = LALERROR | LALWARNING | LALINFO | LALNMEMDBG;

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
#define MAXSTART 10000

/* This will hold the SFT and PLAN */
COMPLEX8Vector *fvec = NULL;
RealFFTPlan *pfwd = NULL;

/* Need prototype: this is POSIX not ANSI, sigh */
#include <config.h>
#include <unistd.h>
#ifndef HAVE_GETHOSTNAME_PROTOTYPE
int gethostname(char *name, int len);
#endif

/* This is an error handler that prints a bit of extra info */
#ifdef __GNUC__
void pout(char *fmt, ...) __attribute__ ((format (printf, 1, 2)));
#else
void pout(char *fmt, ...)  ;
#endif

void pout(char *fmt, ...){
  va_list ap;
  char hostsname[256];
  pid_t processid;

  processid=getpid();

  /* not guaranteed null terminated! */
  gethostname((char *)hostsname, (size_t)255);
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

  printf("%%Environment variable %s = %d\n", valuename, retval);

  return retval;
}

/* Compensate for timing glitches in S1 data.  Returns the number of
   microseconds that has to be added to the time recorded in the frame
   (LSC timestamp) to get the correct GPS time.  Thus a positive value
   indicates that the first datum in a frame file was actually
   recorded at a time later than the timestamp of the frame file.

   Data tables taken from:
   http://www.ligo.caltech.edu/~smarka/drafts/LIGO-T030070-00-D.pdf
*/
int deltatime(const char *instrument, int gpstime, int *valid){
  int (*data)[4];
  int i;
  
  /* corrections for Livingston L1 during S1.  It would probably make
     sense to discard the first few hours for which there IS no timing
     data, but since the first SFT that we make starts at time
     L1.714179401 there is no need. And the first calibrated SFT is at
     time 714299280 */
  int l1time[][4]={
    {714177080, 714855840, -121, 678760},
    {714855841, 715618813, -120, 762972},
    /* value to use if time range not found */
    {0,         INT_MAX,           -121,  0},
  };
  
  /* corrections for Hanford H1 during S1.  Again, it would probably
     make sense to discard the early times for which we have no timing
     data (about the first fifty hours of the run). In fact the first
     SFT is at time 714151661 which is BEFORE the first time listed
     below.  But in fact we don't have calibration information until
     time 714357188.  So this table is OK. */
  int h1time[][4]={
    {714335520, 715618813,   -95, 1283293},
    /* value to use if time range not found */
    {0,         INT_MAX,     -95,  0},
  };
  
  /* corrections for Hanford H2 during S1.  Here the pattern of timing
     is SO irregular that we discard any segments for which we have no
     data. */
  int h2time[][4]={
    {714256920, 714407160,  5463, 150240},
    {714407640, 714507180, -164, 99540},
    {714507300, 714587880,  141, 80580},
    {714697680, 714704820,  995, 7140},
    {714705000, 714786780, -162, 81780},
    {714786900, 715008660,  204, 221760},
    {715008780, 715026480,  448, 17700},
    {715026780, 715079460, -163, 52680},
    {715079580, 715101120,  142, 21540},
    {715101240, 715129260,  447, 28020},
    {715129500, 715165740, -162, 36240},
    {715166640, 715265820, -163, 99180},
    {715266180, 715274580, -102, 8400},
    {715274700, 715296960,  630, 22260},
    {715297440, 715298280, -163, 840},
    {715299420, 715305660, -163, 6240},
    {715306920, 715545120, -163, 238200},
    {715545240, 715556400,  814, 11160},
    {715556640, 715594020, -163, 37380},
    {715594140, 715618800,   81, 24660},
    /* discard data if time range not found */
    {-1,        0,           0,  0}
};

  /* select the correct instrument */
  if (!strcmp(instrument,"H1"))
    data=h1time;
  else if (!strcmp(instrument,"H2"))
    data=h2time;
  else if (!strcmp(instrument,"L1"))
    data=l1time;
  else {
    pout("Unrecognized detector: %s in deltatime(). Not H1, H2, or L1.\n",
	 instrument);
    exit(1);
  }
  
  /* search along list to see if we find the correct time range */
  for (i=0; data[i][0]>=0; i++)
    if (data[i][0]<=gpstime && gpstime<=data[i][1]){
      *valid=1;
      return data[i][2];
    }
  /* value we should use if time range NOT found */
  *valid=0;
  return 0;
}

  /* check a single value of the timing correction */
  void checkone(int gpstime,const char *instrument){
    int valid;
    int dt=deltatime(instrument, gpstime, &valid);
    printf("%s  time %d  correction %d%s\n", instrument, gpstime, dt, valid?"":" INVALID");
    return;
  }

/* check a number of bounary values of the timing correction */
void checktimingcorrections(){
  

  /* L1 checks */
  checkone(10, "L1");
  checkone(10000, "L1");
  checkone(715618799, "L1");
  printf("\n");

  /* H1 checks */
  checkone(10, "H1");
  checkone(10000, "H1");
  checkone(715618799, "H1");
  printf("\n");

  /* H2 checks */
  checkone(0, "H2");    
  checkone(10, "H2"); 
  checkone(10000, "H2"); 
  checkone(714256919, "H2");
  checkone(714256920, "H2");
  checkone(714256921, "H2");
  checkone(714407159, "H2"); 
  checkone(714407160, "H2");
  checkone(715594139, "H2");
  checkone(715594140, "H2"); 
  checkone(715594141, "H2");
  checkone(715618799, "H2"); 
  checkone(715618800, "H2"); 
  checkone(715618801, "H2"); 
}

/* Utility function for cyclically shifting an array "in place".
   Written for clarity and simplicity, not for efficiency.
   shift >= 0 :
               data[0]  (moves to)    -> data[shift]
               data[1]                -> data[shift+1]
	       data[length-1-shift]   -> data[length-1]
               data[length-1-shift+1] -> data[0]
               data[length-1-shift+2] -> data[1]
               ...
   shift < 0 :
               replace shift by length+shift and follow 
	       the rules above.               ...
 */
void shifter(float *data, int length, int shift){
  float *temp=NULL;
  int delta=shift>=0?shift:length+shift;
  int i;

  /* if no shift is required, we are done */
  if (shift==0)
    return;

  /* check that shift range seems reasonable */
  if (abs(shift)>length/8){
    pout("shifter(): shift amount %d seems too big/small for length %d array\n",
	 shift, length);
    exit(1);
  }
  
  /* allocate memory */
  if (!(temp=(float *)LALMalloc(sizeof(float)*length))){
    pout("Unable to allocate %d bytes of memory in shifter\n",
	 sizeof(float)*length);
    exit(1);
  }
  
  /* copy data */
  memcpy(temp, data, sizeof(float)*length);

  /* now do shift */  
  for (i=0; i<length; i++)
    data[(delta+i) % length]=temp[i];
  
  /* free memory and return */
  LALFree(temp);
  return;
}

#if PRINT50
void print50(float *array, char *name){
  int i;
  printf("\n%s\n", name);
  for (i=0; i<50; i++)
    printf("%02d %e\n", i, array[i+176]);
  fflush(stdout);
  return;
}

void print50double(double *array, char *name){
  int i;
  printf("\n%s\n", name);
  for (i=0; i<50; i++)
    printf("%02d %e\n", i, array[i+176]);
  fflush(stdout);
  return;
}
#endif

/* some signal handlers */
void sighandler(int sig){

  if (sig==SIGSEGV)
    fprintf(stderr,"Caught signal SIGSEGV\n");
  if (sig==SIGFPE)
    fprintf(stderr,"Caught signal SIGFPE\n");

  fflush(NULL);

  exit(128+sig);
}


/* routine to track memory usage in different categories */
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

#if TDDOUBLE
  /* for double-precision experiments */
  static REAL8TimeSeries chand;
#endif

  /* Control variable for run-time options */
  int doubledata=getenvval("DOUBLEDATA");
  int windowdata=getenvval("WINDOWDATA");
  int printfreqseries=getenvval("PRINTFREQSERIES");

  /* FrameL frame file */
  FrFile *frfile;
  /* vector holding the frame data */
  FrVect *frvect;        

  /* install signal handlers */
  if (signal(SIGSEGV, sighandler)==SIG_IGN)
    signal(SIGSEGV, SIG_IGN);

  if (signal(SIGFPE, sighandler)==SIG_IGN)
    signal(SIGFPE, SIG_IGN);

#if (0)
  /* check the timing correction code */
  checktimingcorrections();
  exit(0);
#endif

#if (0)
  /* check the shifter code */
  {
    float fun[30];
    int i;
    
    for (i=0; i<30; i++)
      fun[i]=(float) i;
    
    shifter(fun, 30, -2);
    
    for (i=0;i<30;i++)
      printf("%2d: %f\n", i, fun[i]);
    
  }
  exit(0);
#endif

  printf("%%Normal startup\n");
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

#if USETDCALDATA
  sprintf(chname,"%s:LSC-STRAIN",argv[4]);
#else
  sprintf(chname,"%s:LSC-AS_Q",argv[4]);
#endif

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
#if HIGHFREQ
  filterpar.f2 = 300.0;
#else
  filterpar.f2 = 40.0;
#endif
  filterpar.a2 = 0.5;
 
/* values that are 'not given' = out of range */
  filterpar.f1 = -1.0;
  filterpar.a1 = -1.0;
  
  /* Initialize frame library with correct file list */
  sprintf(framelist,"%s/jobdata.%05d.ffl",argv[2],jobnum);
  opencount=0;

#if TRACKMEMUSE
    printf("[ENTRY]      Memory usage at iteration %d is:\n", count);
    printmemuse();
#endif

  while (!(frfile = FrFileINew(framelist))){
    pout( "Couldnt open frame file list %s\n", framelist);
    if (opencount++<10)
      sleep(10);
    else
      return(3);
  }
  
#if TRACKMEMUSE
    printf("[After FFL]  Memory usage at iteration %d is:\n", count);
    printmemuse();
#endif

  /* Main loop to write SFTs */
  for (count=0;count<nstarts;count++) {
    struct stat buff;
    char sftname[256];
    int filesize=(INT4)(DF*tbase);
    int lsffl=0,verifyframes=0;

#if TRACKMEMUSE
    printf("[TOP]        Memory usage at iteration %d is:\n", count);
    printmemuse();
#endif

#if TIMESHIFT
    int microsec, valid=0;
#endif

    /* time of correct start */
    epoch.gpsSeconds=starts[count];
    epoch.gpsNanoSeconds=0;
    
#if TIMESHIFT
    /* compute timeshift needed */
    microsec=deltatime(argv[4], epoch.gpsSeconds, &valid);
    
    /* if there is no valid timeshift info, skip SFT */
    if (!valid) {
      printf("SFT at time %d not made.  No valid timing data\n", epoch.gpsSeconds);
      continue;
    }
#endif
    
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

#if USETDCALDATA
    frvect = FrFileIGetVProc(frfile, chname, epoch.gpsSeconds, tbase, 0);
#else 
    frvect = FrFileIGetVAdc(frfile, chname, epoch.gpsSeconds, tbase, 0);
#endif

#if TRACKMEMUSE
    printf("[After Fr]   Memory usage at iteration %d is:\n", count);
    printmemuse();
#endif
   
    /* This block will detect permission denied and other access errors */
    if (errno || frvect==NULL || frvect->next){
      int saveerrno=errno;

      /* if we haven't printed ls of the files, do so now (but just once!)*/
      if (!lsffl){
	char command[256];
	if (saveerrno)
	  pout("System error when reading Frame data: %s\n", strerror(saveerrno));
	pout("Following is output of ls -l on Frame file list:\n");
	sprintf(command, "cat %s | xargs ls -l 1>&2", framelist);
	system(command);
	lsffl=1;
      }

      if (saveerrno)
	return 3;
    }
    
    /* This block detects the case where the frame file list did not contain all data */
    if (frvect==NULL) {
      pout( "%s data missing between times %d and %d\n",argv[4], epoch.gpsSeconds,epoch.gpsSeconds+tbase);

      if (!verifyframes){
	char command[256];
	pout("Using FrCheck to validate frame data:\n");
	sprintf(command, "cat %s | xargs --maxlines=1 /home/ballen/projects/LIGO/src/v6r06/Linux-i686/FrCheck -d 1 -i 1>&2", framelist);
	system(command);
	verifyframes=1;
      }
      continue;
    }
    
    /* Check that data type is correct */

#if USETDCALDATA
    if (frvect->type != FR_VECT_8R){
      pout( "Wrong data type (notFR_VECT_8R) found in frame!\n" );
      return(5);
    }
#else
    if (frvect->type != FR_VECT_4R){
      pout( "Wrong data type (notFR_VECT_4R) found in frame!\n" );
      return(5);
    }
#endif

    /* check for gaps */
    if (frvect->next){
      pout( "%s data between times %d and %d had a gap\n",argv[4], epoch.gpsSeconds,epoch.gpsSeconds+tbase);
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
      
#if PRINT50
#if USETDCALDATA
      print50double(frvect->dataD, "DOUBLE_PRECISION_FRAMEDATA");
#else
      print50(frvect->dataF, "FRAMEDATA");
#endif
#endif
      
      /* create structure to store channel data  */
      /* NOTE: if we are using time-domain calibrated data, then we
	 are going to convert a double-precision quantity to a
	 single-precision one here */

      chan.data = NULL;
      LALSCreateVector(&status, &chan.data, npts);
      TESTSTATUS(&status);
      
#if TRACKMEMUSE
      printf("[After SVEC] Memory usage at iteration %d is:\n", count);
      printmemuse();
#endif
      /* set up LAL time series object sample interval, time epoch */
      strcpy(chan.name, chname);
      chan.deltaT = 1.0/SRATE;
      chan.epoch = epoch;
      chan.data->length = npts;
      
      /* copy data */
      for (i=0;i<npts;i++){
	chan.data->data[i]=frvect->dataD[i];

#if QUANTIZATION_TEST
	{
	  /* Note -- this mask MAY depend on big versus little endian.
	     Currently correct for x86 (little endian, format is
	     E...ESM...M, where E=exponent, S=sign, M=mantissa and MSB
	     is to the left.) */
	  const unsigned int mask = (unsigned int)(~(0x01));
	  union {
	    float f;
	    unsigned int i;
	  } maskfloat;
	  
	  /* zero LSB for testing quantization noise */
	  maskfloat.f=chan.data->data[i];
	  maskfloat.i &= mask;
	  chan.data->data[i]=maskfloat.f;
	}
#endif
	
#if NEWNORM
	/* Normalize data according to the GEO SFT specifications */
	chan.data->data[i]*=(((double)DF)/(0.5*(double)SRATE));
#endif
      }

/*       for (i=0; i<npts; i++) */
/* 	fprintf(stdout,"%e\n",chan.data->data[i]); */



#if PRINT50
      print50(chan.data->data, "LAL_TIMESERIES");
#endif

      /* free framevec -- no longer needed */
      FrVectFree(frvect); 
      frvect=NULL;

#if PRINTTEV
      for (i=0; i<npts; i++)
	fprintf(stderr, "%f\n", chan.data->data[i]);
      exit(0);
#endif 

#if TDDOUBLE
      /* create structure to store channel data */
      chand.data = NULL;
      LALDCreateVector(&status, &chand.data, npts);
      TESTSTATUS(&status);
      strcpy(chand.name, chname);
      chand.deltaT = 1.0/SRATE;
      chand.epoch = epoch;
      chand.data->length = npts;

#if TRACKMEMUSE
      printf("[After DVEC] Memory usage at iteration %d is:\n", count);
      printmemuse();
#endif

#if HIGHPASS_FILTER
      /* put into doubles */
      for (i=0; i<npts; i++)
	chand.data->data[i]=chan.data->data[i];
      
      /* filter */
      LALButterworthREAL8TimeSeries(&status, &chand, &filterpar);
      TESTSTATUS(&status);
      
      /* and copy back */
      for (i=0; i<npts; i++)
	chan.data->data[i]=chand.data->data[i];
#endif /* HIGHPASS_FILTER */


      /* to save memory, get rid of vector */
      LALDDestroyVector( &status, &chand.data );
      TESTSTATUS( &status );
      chand.data=NULL;
#else /* #if TDDOUBLE */
#if HIGHPASS_FILTER

      /* apply high-pass Butterworth filter */
      LALButterworthREAL4TimeSeries(&status, &chan, &filterpar);
      TESTSTATUS(&status);
#endif /* HIGHPASS_FILTER */
#endif /* TDDOUBLE */

#if TRACKMEMUSE
    printf("[After FILT] Memory usage at iteration %d is:\n", count);
    printmemuse();
#endif

#if PRINT50
      print50(chan.data->data, "LAL_FILTERED");
#endif

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

#if PRINT50
      print50(chan.data->data, "LAL_WINDOWED");
#endif

#if TIMESHIFT
      /* correct data for timing errors.  This is a bit of a hack.  We
	 can do it either just before or just after the FFT. This is
	 equivalent if the timing error equals an integer multiple of
	 the sample time.  Otherwise the correct way to do it is by
	 modifying the FFT output by multiplying by exp(2pi i f dt).
	 Here, we take the easy route, which gets it right to
	 sufficient accuracty (+- 30 us) for our purposes.  Note that
	 there are a number of approximations built in to either way
	 of doing it, all of which boil down to: it works if the
	 timing errors are very small compared to the SFT time.
	 
	 A positive value of deltatime() means that the correct time
	 stamp for the first datum recorded in the frame is LATER than
	 the frame time stamp.  A negative value of deltatime() means
	 that the correct time stamp for the first datum recorded in
	 the frame is EARLIER than the frame time stamp.
      */
      if (microsec){
	int nshift;
	double fshift=1.e-6*microsec*SRATE;
	
	if (fshift>=0.0)
	  nshift=(int)(fshift+0.5);
	else
	  nshift=(int)(fshift-0.5);
	
	if (nshift){
	  shifter(chan.data->data, npts, nshift);
	  printf("Shifting %s data sample at time %d forward %d samples\n",
		 argv[4], epoch.gpsSeconds, nshift);
	}
      }
#endif
      
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
      
#if PRINT50
      print50(chan.data->data, "SHIFTED");
#endif

      /* Create vector to hold signal frequency series.  To save memory,
	 do this only when needed */
      LALCCreateVector(&status, &fvec, (UINT4)len2);
      TESTSTATUS(&status);

#if TRACKMEMUSE
    printf("[After FVEC] Memory usage at iteration %d is:\n", count);
    printmemuse();
#endif
 
      /* Compute measured plan for FFTW.  To save memory, do only when needed */
      LALCreateForwardRealFFTPlan(&status, &pfwd, (UINT4)npts, 0);
      TESTSTATUS(&status);

#if TRACKMEMUSE
    printf("[After PLAN] Memory usage at iteration %d is:\n", count);
    printmemuse();
#endif

      /* take forward FFT */
      LALForwardRealFFT(&status, fvec, chan.data, pfwd);
      TESTSTATUS(&status);
      
      /* to save memory, destroy plan */
      LALDestroyRealFFTPlan(&status, &pfwd);
      TESTSTATUS( &status );
      pfwd=NULL;

#if PRINT50
      print50(chan.data->data, "FFT");
      exit(0);
#endif

      LALSDestroyVector( &status, &chan.data );
      TESTSTATUS( &status );
      chan.data=NULL;

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

      /* to save memory, get rid of vector when not needed any more */
      LALCDestroyVector(&status, &fvec);
      TESTSTATUS( &status );
      fvec=NULL;
    }
  }
  

  /* Clean up and exit */
  FrFileIEnd(frfile);
    
  LALCheckMemoryLeaks();
  
  printf("%%Normal exit\n");
  fflush(stdout);
 
  return 0;
}
