/*
*  Copyright (C) 2003-9 Bruce Allen, Peter Shawhan
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/* 
   multipulsar injection routine, written for E10/S3 by Bruce Allen,
   10/2003.  Calls to Signal Injection Library added by Peter
   Shawhan.  
   
   2005/02 - Duncan Brown renamed variable to avoid Condor conflict

   2005/02 - Reinhard Prix removed the actuation scaling options

   28 May 2009: renamed this code from 's3inject' to 'psinject' in lalsuite-GIT
*/

#define _GNU_SOURCE

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <errno.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <signal.h>
#include <math.h>
#include <sys/stat.h>
#ifdef ONLINE
#include "SIStr.h"
#endif

#include <FrameL.h>

/* include lalsuite git IDs for this code */
#include <lal/lalGitID.h>
#include <lalappsGitID.h>


#define MAXPULSARS 64
/* blocksize for gathering results from children */
#define BLOCKSIZE 16384
/* sample rate of output */
#define SRATE 16384

/* Give 'ident' something to recognize in the executable */
const char *cvsid="$Id$";

/* Flag to initiate graceful shutdown, set by sighandler */
volatile int shutdown_pulsar_injection = 0;
int npulsars=0;
float data[BLOCKSIZE];
float total[BLOCKSIZE];
float testval=1234.5;
char *directory;
char *channel=NULL; 
int  gpstime=0;
char *programname;
int show=0;
int do_axis=0;
int do_text=0;
int write_frames=0;
int secs_per_framefile=60;
long long count=0;
long blocks=0;
const char *ifo_name = NULL;
const char *actuation = NULL;

/* capacity, occupancy, and pointers to input buffers. Cleared to zero
   on starup because they are initialized data.  Units of buffer size
   and offsets are counted in floats!
*/
int bufflen[MAXPULSARS];
int readfrombuff[MAXPULSARS];
float *buffs[MAXPULSARS];

double calamp[3]={0.0, 0.0, 0.0};

/* 
   Calibration line frequencies. THE CODE ASSUMES THAT THESE ARE
   POSITIVE.  Order is Low, Medium, High frequency.  If you change
   these, you MUST choose a frequency that can be exactly represented
   as an IEEE754 floating point floats. For example
   52+1/4+1/32+1/64=52.296875
*/
double calfreq[3]={
  52.296875,         /* 52+1/4+1/32+1/64 Hz */
  166.6875,          /* 166+1/2+1/8+1/16 Hz */
  927.6875           /* 927+1/2+1/8+1/16 Hz */
};

/*
  Calibration line fiducial time.
*/
int tfiducial=751680013;

FILE *fp[MAXPULSARS];

/* forward declaration */
void syserror(int showerrno, const char *fmt, ...);
void sighandler(int sig);
void usage(FILE *filep);
int parseinput(int argc, char **argv);
void OutputVersion ( void );

/* signal handler to monitor the status of children catch SIGCHLD */
void sighandler(int sig){
  pid_t pid;
  int status;
  long long seconds=((long long)blocks)*((long long)BLOCKSIZE)/((long long)SRATE);
  int minutes=seconds/60;
  int hours=minutes/60;
  int days=hours/24;
  int secs=seconds % 60;
  
  syserror(0, "made %d days %d hours %d minutes %d seconds of data\n",
	   days, hours % 24, minutes % 60, secs);

  syserror (0, "Received signal %d\n", sig);
  
  if (sig==SIGTERM) {
    syserror(0, "received SIGTERM, initiating graceful shutdown...\n");
    shutdown_pulsar_injection = 1;
  }
  
  if (sig!=SIGCHLD) 
    return;
  
  if((pid=waitpid(-1, &status, WNOHANG | WUNTRACED))>0) {
    /* if child stopped, make log entry then return */
    if (WIFSTOPPED(status)) {
      syserror(0, "Subprocess [PID=%d] stopped because it caught signal %d [%s]\n",
	       (int)pid, WSTOPSIG(status), strsignal(WSTOPSIG(status)));
      return;
    }
    
    /* otherwise something more serious is wrong... */
    syserror(0, "Subprocess [PID=%d] is misbehaving.\n", (int)pid);
    if (WIFEXITED(status))
      syserror(0, "Subprocess [PID=%d] or shell did exit(%d)\n", (int)pid, WEXITSTATUS(status));
    
    if (WIFSIGNALED(status))
      syserror(0, "Subprocess [PID=%d] terminated because it caught signal %d [%s]\n", 
	       (int)pid, WTERMSIG(status), strsignal(WTERMSIG(status)));      
    exit(1);
  }
  else
    syserror(1, "waitpid() returned -1.  Call Bruce...\n");
  return; 
}

/* Like perror() but takes variable numbers of arguments and includes program name */
void syserror(int showerrno, const char *fmt, ...){
  char *thiserror = NULL;
  pid_t pid = getpid();
  time_t t = time(NULL);
  va_list ap;
  /* initialize variable argument list  */
  va_start(ap,fmt);
  /* print a standardized header with time and PID info */
  if (showerrno && errno && (thiserror=strerror(errno)))
    fprintf(stderr,"%s [PID=%d] %.24s: %s: ", programname, (int)pid, ctime(&t), thiserror);
  else
    fprintf(stderr,"%s [PID=%d] %.24s: ", programname, (int)pid, ctime(&t));
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fflush(stderr);
  return;
}

/* usage message */
void usage(FILE *filep){    
  fprintf(filep,
          "--------------------------------------------------------------------------------\n"
          "%s: %s\n" 
          "--------------------------------------------------------------------------------\n"
	  "Options are:\n"
	  "-h            THIS help message\n"
          "-v            VCS ID information\n"
	  "-n INT        Number of pulsars to simulate: 0, 1, 2, ..., %d\n"
	  "-d DIRECTORY  Directory containing command line files in.0, in.1, ...\n"
	  "              Default is: . (current working directory)\n"
	  "-e CHANNEL    Inject excitation into CHANNEL\n"
	  "-D            Turn on Debug output for routines in signal injection library\n"
	  "-G INTGPS     GPS time (integer seconds) that will be passed on command line\n"
	  "-T            Print human readable text rather than binary to stdout\n"
	  "-X            Include X-axis in human-readable text output\n"
	  "-s            Print commands that would be fork(2)ed, then exit(0)\n"
	  "              ---------- --------------------------------------------------\n"
          "-L DOUBLE     | Inject  calibration lines. Here L,M,H denote Low/Mid/High |\n"
          "-M DOUBLE     | frequency.  The DOUBLE values specifies the corresponding |\n"
	  "-H DOUBLE     | amplitude. If NOT given, the amplitude defaults to 0.0    |\n"
	  "              ------------- -----------------------------------------------\n"
	  "-p            Print the calibration line frequencies in Hz then exit(0)\n"
	  "-I STRING     Detector: LHO, LLO, GEO, VIRGO, TAMA, CIT, ROME [REQUIRED]\n"
	  "-A STRING     File containing detector actuation-function     [OPTIONAL]\n"
          "-F INT        Keep N frame files on disk.  If N==0 write all frames immediately.\n"
	  "-S INT        Number of 1-second frames per frame file (default 60).\n"
          "--------------------------------------------------------------------------------\n"
	  , programname, lalappsGitCommitID, MAXPULSARS
	  );
  return;
}


int parseinput(int argc, char **argv){
  
  int c;
  const char *optionlist="hL:M:H:n:d:e:DG:TXspI:A:F:vS:";
  opterr=0;
  
  /* set some defaults */
  directory = strdup(".");

  programname=argv[0];

  while (-1 != (c = getopt(argc, argv, optionlist))) {
    char *end;
    double tempamp;
    switch (c) {
    case 'v':
        OutputVersion();
        exit(0);
    case 'p':
      printf("The calibration line frequencies are:\n"
	     "  L: %18.14f Hz\n"
	     "  M: %18.14f Hz\n"
	     "  H: %18.14f Hz\n",
	     calfreq[0], calfreq[1], calfreq[2]);
      exit(0);
    case 'h':
      /* usage message */
      usage(stdout);
      exit(0);
    case 'n':
      /* number of pulsars to simulate */
      npulsars=atoi(optarg);
      if (npulsars<0 || npulsars>MAXPULSARS){
	syserror(0, "%s: Number of pulsars (-n argument = %d) must be non-negative and less than %d\n",
		argv[0], npulsars, MAXPULSARS+1);
	exit(1);
      }
      break;
    case 'L':
    case 'M':
    case 'H':
      /* calibration-line amplitude */
      tempamp=strtod(optarg, &end);
      if (*end){
	syserror(1, "-%c %s is invalid. -%c takes a double-precision amplitude.\n", c, optarg, c);
	exit(1);
      }
      /* assign amplitude to the correct line */
      if (c=='L')
	calamp[0]=tempamp;
      else if (c=='M')
	calamp[1]=tempamp;
      else
	calamp[2]=tempamp;
      break;
    case 'd':
      /* directory path */
      if (directory) free(directory);
      if (!(directory=strdup(optarg))){
	syserror(1, "Out of memory to duplicate -d directory path %s\n", optarg);
	exit(1);
      }
      break;
    case 'e':
#ifdef ONLINE
      /* Excitation channel into which to inject */
      if (!(channel=strdup(optarg))){
	syserror(1, "Out of memory to duplicate -e channel name %s\n", optarg);
	exit(1);
      }
#else
      syserror(0, "The -e option to enable online signal injection requires compilation with -DONLINE\n");
      exit(1);
#endif
      break;
    case 'D':
#ifdef ONLINE
      /* Turn on debugging for SIStr library routines */
      SIStr_debug++;
#else
      syserror(0, "The -D option to enable SIStr debugging requires compilation with -DONLINE\n");
      exit(1);
#endif
      break;
    case 'G':
      /* GPS time to pass as argument */
      gpstime=atoi(optarg);
      break;
    case 'T':
      /* enable text output */
      do_text=1;
      break;
    case 'X':
      /* include x axis in text output */
      do_axis=1;
      break;
    case 's':
      /* show commands rather than executing them */
      show=1;
      break;
    case 'I':
      ifo_name = optarg;
      break;
    case 'A':
      actuation = optarg;
      break;
    case 'F':
	{
	    int how_many = atoi(optarg);
	    if (how_many < 0) {
		syserror(0,"%s: fatal error, argument -F %d must be non-negative.\n", argv[0], how_many);
		exit(1);
	    }
	    write_frames = 1 + how_many;
	    break;
	}
    case 'S':
	secs_per_framefile = atoi(optarg);
        if (secs_per_framefile < 1) {
	    syserror(0,"%s: fatal error, argument -S %d must be at least 1 second.\n", argv[0], secs_per_framefile);
	    exit(1);
	}
        if (secs_per_framefile > 3600) {
	    syserror(0,"%s: caution, argument -S %d seconds is more than one hour!\n", argv[0], secs_per_framefile);
	}
	break;
    default:
      /* error case -- option not recognized */
      syserror(0,"%s: Option argument: -%c unrecognized or missing argument.\n"
	       "\t\tPlease use the '-h' option to print usage message\n"
	       ,argv[0], optopt);
      exit(1);
      break;
    }
  }

  /* sanity checks on command line arguments */
  if (do_axis && !do_text) {
    syserror(0, "The -X (axis) option only works with the -T (human readable) option\n");
    exit(1);
  }
  if (do_text && channel){
    syserror(0, "Can't use both '-T' and '-e' together\n");
    exit(1);
  }
  if (channel && (strlen(channel) <3 || channel[2] != ':' )) {
    syserror(0, "Excitation channel %s not of form CC:CCC...\n", channel );
    exit(1);
  }
  if ( ifo_name == NULL ) {
    syserror(0, "You must specify the IFO name (-I)\n");
    exit(1);
  }

  
#ifndef ONLINE
  if (channel) {
    syserror(0, "Can't do exicitations. Code not compiled with ONLINE defined\n");
    exit(1);
  }
#endif      
  
  return 0;
}

#define MAXLINE 1024
int main(int argc, char *argv[]){
  int i,j;
#ifdef ONLINE
  char *cwd;
  int iarg, status;
  SIStream sis;
  char info[MAXLINE];
  memset(&sis, '\0', sizeof(SIStream));
#endif

  parseinput(argc, argv);
  
  syserror(0, "Starting up\n");

  /* install signal handler to catch SIGCHLD. Note that solaris uses
     SysV rather than BSD semantics and doesn't automaticall restart system
     calls like fread and fwrite.  So we use sigaction() rather than
     signal() so that we can set SA_RESTART. */
  {
    struct sigaction sig;
    memset(&sig, '\0', sizeof(sig));
    sig.sa_flags=SA_RESTART;
    sig.sa_handler=sighandler;
    if (sigaction(SIGCHLD, &sig, NULL)){
      syserror(1, "Unable to install signal handler for messages about troubled children\n");
    }
    if (sigaction(SIGTERM, &sig, NULL)){
      syserror(1, "Unable to install signal handler for logging output rate data and terminating\n");
    }
    if (sigaction(SIGUSR1, &sig, NULL)){
      syserror(1, "Unable to install signal handler for logging output rate data\n");
    }
  }
  
  for (i=0; i<npulsars; i++){
    char command[MAXLINE];
    char filename[MAXLINE];
    FILE *fpc=NULL;
    int length=0;
    char *newlineloc=NULL;

    /* construct file name */
    if (snprintf(filename, MAXLINE, "%s/in.%d", directory, i)>MAXLINE-1){
      syserror(0, "%s: file name %s/in.%d has more than MAXLINE=%d characters\n",
	      programname, directory, i, MAXLINE);
      exit(1);
    }
    
    /* open file */
    if (!(fpc=fopen(filename, "r"))) {
      syserror(1, "Can't open file %s for reading\n", filename);
      exit(1);
    }

    /* read command line from file */
    if (!(fgets(command, MAXLINE, fpc))){
      syserror(1, "Command line file %s was empty!\n", filename);
      exit(1);
    }
    fclose(fpc);
    
    /* check that contents are not too large */
    length=strlen(command);
    if (length>=MAXLINE-1) {
      syserror(0, "Command line file %s has line >= to MAXLINE=%d characters!\n",
	       filename, MAXLINE);
      exit(1);
    }
    
    /* replace first NEWLINE to null terminate string */
    if ((newlineloc=index(command, '\n')))
	*newlineloc='\0';

    /* append additional arguments to command line */
    /* GPS starttime */
    length=strlen(command);
    if (snprintf(command+length, MAXLINE-length, " -G %d", gpstime)>MAXLINE-length-1){
      command[length]='\0';
      syserror(0, "%s: command line has >= MAXLINE=%d characters\n", programname, MAXLINE);
      exit(1);
    }
    /* IFO */
    length=strlen(command);
    if ( snprintf(command+length, MAXLINE-length, " -I %s", ifo_name) > MAXLINE-length-1 ) {
      command[length]='\0';
      syserror(0, "%s: command line has >= MAXLINE=%d characters\n", programname, MAXLINE);
      exit(1);
    }
    /* Actuation-function if given */
    if (actuation)
      {
	length=strlen(command);
	if ( snprintf(command+length, MAXLINE-length, " --actuation=%s", actuation) > MAXLINE-length-1 ) {
	  command[length]='\0';
	  syserror(0, "%s: command line has >= MAXLINE=%d characters\n", programname, MAXLINE);
	  exit(1);
	}
      }

    /* now either show the command or execute it */
    if (show)
      printf("[%02d] %s\n", i, command);
    else {
      errno=0;
      if (!(fp[i]=popen(command, "r")) || errno){
	syserror(1, "Unable to popen(3) %s\n", command);
	exit(1);
      }
    }
  }
  
  /* a useful option for debugging -- show the output */
  if (show) {
    if (!npulsars)
      printf("%s: Warning: n=0 so an infinite-length zero strength signal will be made!\n", argv[0]);
    exit(0);
  }

#if 0
  {
    pid_t pid;
    int status;
    /* wait a couple of seconds, then check that all processes are running happily! */
    sleep(2);
    pid=waitpid(-1, &status, WNOHANG | WUNTRACED);
    if (pid) {
      syserror(0,"Subprocess with PID=%d is misbehaving.\n", (int)pid);
      if (WIFEXITED(status))
	syserror(0, "Subprocess or shell did exit(%d)\n", WEXITSTATUS(status));
      
      if (WIFSIGNALED(status))
	syserror(0, "Subprocess terminated because it caught signal %d [%s]\n", 
		WTERMSIG(status), strsignal(WTERMSIG(status)));
      exit(1);
    }
  }
#endif
    
  /* processes opened, read data*/  
  for (i=0; i<npulsars; i++){
    if (fread(&testval, sizeof(float), 1, fp[i]) != 1) {
      syserror(1, "Could not read first float 1234.5 from %d'th signal source\n", i);
      exit(1);
    } else if (testval != 1234.5) {
      syserror(0, "First value (%f) from %d'th signal source was not 1234.5\n", testval,  i);
      exit(1);
    } else if (fread(bufflen+i, sizeof(int), 1, fp[i]) != 1) {
      syserror(1, "Could not read buffer size from %d'th signal source\n", i);
      exit(1);
    } else if (bufflen[i]<SRATE || bufflen[i]>SRATE*60) {
      syserror(0, "Bad buffer size %d floats from %d'th signal source (expect %d <= size <= %d)\n", bufflen[i], i, SRATE, 60*SRATE);
      exit(1);
    } else if (bufflen[i]% BLOCKSIZE) {
      syserror(0, "Bad buffer size %d floats from %d'th signal source NOT a multiple of BLOCKSIZE=%d\n", bufflen[i], i, BLOCKSIZE);
      exit(1);  
    } else if (!(buffs[i]=(float *)calloc(bufflen[i], sizeof(float)))) {
      syserror(1, "Can't allocate buffer of %d floats for %d'th signal source\n", bufflen[i], i);
      exit(1);
    }
    /* ensure that we read buffers on first pass */
    readfrombuff[i]=bufflen[i];
  }
  
#if 0
  /* are we writing frames? */
  if (write_frames) {
      /* put whatever is needed here to init frame lib */
  }
#endif

  /* are we calling the excitation engine directly? */
  if (channel) {

#ifdef ONLINE
    /* Report some information about this injection client */
    cwd = getcwd( NULL, 256 );
    if ( cwd ) {
      sprintf( info, "%s %s", argv[0], cwd );
    } else {
      sprintf( info, "%s unknown_directory", argv[0] );
    }
    free( cwd );
    SIStrAppInfo( info );
    
    /* Open the Signal Injection Stream */
    status = SIStrOpen( &sis, channel, 16384, (double) gpstime );
    if ( SIStr_debug ) { 
      syserror(0, "SIStrOpen() returned %d\n", status ); 
    }
    if ( status != SIStr_OK ) {
      syserror(0, "SIStrOpen() error opening SIStream: %s\n", SIStrErrorMsg(status) );
      exit(2);
    }
#endif
  }
  else 
    if (do_text)
    printf("1234.5\n");
  else {
    /* write out 1234.5 */
    if (1!=fwrite(&testval, sizeof(float), 1, stdout)){
      syserror(1, "Unable to output key value 1234.5\n");
      exit(1);
    }
  }
  
  /* now read data blocks unless a SIGTERM has set shutdown */
  while (!shutdown_pulsar_injection) {
    int num=0;
    int line;
    int tdelt=gpstime-tfiducial;

    /* clear block that will contain totals */
    for (j=0; j<BLOCKSIZE; j++)
      total[j]=0.0;
    
    /* if needed, insert calibration line(s) */
    for (line=0; line<3; line++){
      if (calamp[line] != 0.0) {
	/* normal int and double variables for integer/fractional
	   time.  In this and in the code that follows, _fra refers to
	   the fractional part [0,1) and _int refers to the integer
	   part. */

	int    dt_int;
	double dt_fra;

	/*  Avoid long longs in inner loop as they are slow. */
	long long t_rem=blocks, t_int=BLOCKSIZE;

	/* line amplitude and frequency (integer + fractional parts) */
	double f_fra  = calfreq[line];
	int    f_int  = (int)f_fra;
	f_fra -= f_int;
	
	/* integer and fractional time offsets of first sample */
	t_rem   *= t_int;
	t_int    = t_rem;
	t_int   /= SRATE;
	t_rem   -= t_int*SRATE;
	t_int   += tdelt;

	dt_int   = t_int;
	dt_fra   = t_rem;
	dt_fra  /= SRATE;

	for (j=0; j<BLOCKSIZE; j++) {
	  double cycles1, cycles2, cycles3;
	  double tlocal_fra  = dt_fra+(double)j/((double)SRATE);
	  int    tlocal_int  = (int)tlocal_fra;
	  tlocal_fra        -= tlocal_int;
	  tlocal_int        += t_int;
	  cycles1            = f_fra*tlocal_int;
	  cycles1           -= (int)cycles1;
	  cycles2            = tlocal_fra*f_int;
	  cycles2           -= (int)cycles2;
	  cycles3            = tlocal_fra*f_fra;
	  cycles3           -= (int)cycles3;

	  total[j]=calamp[line]*sin(2*M_PI*(cycles1+cycles2+cycles3));
	}
      }
    }
    
    /* loop over the different pulsars */
    for (i=0; i<npulsars; i++) {
      float *where;
      
      if (readfrombuff[i]==bufflen[i]){
	/* read data from the i'th signal */
	if (bufflen[i]!=(num=fread(buffs[i], sizeof(float), bufflen[i], fp[i]))){
	  syserror(1, "Only read %d floats (not %d) from %d'th signal source\n", num, bufflen[i], i);
	  exit(1);
	}
#ifdef DEBUGTIMING
	syserror(0, "just read %d seconds of data from %d'th signal\n", num/SRATE, i);
#endif
	readfrombuff[i]=0;
      }
      
      /* location of signal in buffer */
      where=buffs[i]+readfrombuff[i];
      
      /* add i'th pulsar to the total signal */
      for (j=0; j<BLOCKSIZE; j++)
	total[j]+=where[j];
      
      readfrombuff[i]+=BLOCKSIZE;
    }
    
    /* now output the total signal to frames */
    if (write_frames) {

	static FrFile *oFile = NULL;
	FrameH *frame;
	FrSimData *sim;

	int m, level = 0;
	double sampleRate = SRATE;
	long ndata = BLOCKSIZE;

	char framename[256];
	static int counter=0;
	struct stat statbuf;
	
	/* Create next frame file. This leads to a names like:
	   CW_Injection-921517800-60.gwf
	*/
	sprintf(framename, "CW_Injection");

	/* sprintf(framename, "CW_Injection_%d.gwf", counter + gpstime); */
	frame = FrameHNew(framename);
	if (!frame) {
	    syserror(1, "FrameNew failed (%s)", FrErrorGetHistory());
	    exit(1);
	}
	
	/* set up GPS time, sample interval, copy data */
	frame->GTimeS = gpstime + counter;
	frame->GTimeN = 0;
	frame->dt = ndata/sampleRate;  
	sim = FrSimDataNew(frame,"CW_simulated", sampleRate, ndata, -32);
	for (m=0; m < ndata; m++) {
	    sim->data->dataF[m] = total[m];
	}
    
	
	/* open file, write file, close file */
	if  (!(counter % secs_per_framefile)) {
	    if (oFile) {
		FrFileOEnd(oFile);
		oFile = NULL;
	    }

	    /* Does the user want us to keep a limited set of frames on disk? */
	    if (write_frames>1) {
		char listname[256];
		int watchtime = gpstime + secs_per_framefile*(counter/secs_per_framefile - write_frames + 1);
		sprintf(listname, "CW_Injection-%d-%d.gwf",  watchtime, secs_per_framefile);
		while (!stat(listname, &statbuf)) {
		    /* if enough files in place, sleep 0.1 seconds */
		    struct timespec rqtp;
		    rqtp.tv_sec = 0;
		    rqtp.tv_nsec = 100000000; 
		    nanosleep(&rqtp, NULL);
		}
	    }
	    oFile = FrFileONewM(framename, level, argv[0], secs_per_framefile);
	}

	if (!oFile) {
	    syserror(1, "Cannot open output file %s\n", framename);
	    exit(1);
	}
	if (FR_OK != FrameWrite(frame, oFile)) {
	    syserror(1, "Error during frame write\n"
		   "  Last errors are:\n%s", FrErrorGetHistory());
	    exit(1);
	}  
	/* FrFileOEnd(oFile); */
	FrameFree(frame);

	counter++;
    }

    /* now output the total signal... */
    else if (channel){
#ifdef ONLINE
      /* ... to the excitation engine ... */
      status = SIStrAppend( &sis, total, BLOCKSIZE, 1.0 );
      if ( SIStr_debug >= 2 )
	syserror(0, "SIStrAppend() returned %d\n", status );
      if ( status != SIStr_OK ) {
	syserror(0, "SIStrAppend() error adding data to stream: %s\n",
		  SIStrErrorMsg(status) );
	break;
      }
#endif
    }
    /* ... or as text ... */
    else if (do_text){
      if (do_axis){
	/* ... either as text with both x and y axes ... */
	long long x1=gpstime;
	long long E9=1000000000;
	x1*=E9;
	
	for (j=0; j<BLOCKSIZE; j++){
	  long long x2=count, x3;
	  x2 *= E9;
	  x2  /= (SRATE);
	  x2 += x1;
	  x3 =  x2;
	  x3 /= E9;
	  x2 -= x3*E9;
	  printf("%lld.%09lld %f\n", x3, x2, total[j]);
	  count++;
	}
      }
      else {
	/* ... or as y-axis text only ... */
	for (j=0; j<BLOCKSIZE; j++)
	  printf("%f\n", total[j]);
      }
    }
    else {
      /* ... or in raw binary form. */
      if (BLOCKSIZE!=(num=fwrite(total, sizeof(float), BLOCKSIZE, stdout))){
	syserror(1, "Only wrote %d values (not %d)\n", num, BLOCKSIZE);
	exit(1);
      }
#ifdef DEBUGTIMING
      syserror(0, "Just sent %d seconds of data to system\n", BLOCKSIZE/SRATE);
      sleep(BLOCKSIZE/SRATE);
#endif
    }

    /* increment counter of blocks sent out */
    blocks++; 
  }
  
  /* We'll be exiting, so clean up */
  if (channel) {
#ifdef ONLINE
    /* Close the signal injection stream */
    status = SIStrClose( &sis );
    if ( SIStr_debug )
      syserror(0, "SIStrClose returned %d\n", status );
    if ( status != SIStr_OK ) {
      syserror(0, "Error while closing SIStream: %s\n", SIStrErrorMsg(status) );
      exit(2);
    }
#endif
  }

#ifdef _LINUX  
  /* shut down signal handler for SIGCHLD */
  {
    struct sigaction sig;
    memset(&sig, '\0', sizeof(sig));
    sig.sa_flags=SA_RESTART | SA_NOCLDSTOP;
    sig.sa_handler=SIG_IGN;
    if (sigaction(SIGCHLD, &sig, NULL)){
      syserror(1, "Unable to install signal handler for exiting\n");
    }
    if (sigaction(SIGPIPE, &sig, NULL)){
      syserror(1, "Unable to install signal handler for exiting\n");
    }
  }

  for (i=0; i<npulsars; i++) {
    int status;

    __fpurge(fp[i]);
    status=pclose(fp[i]);
    if (status!=-1)
      syserror(0, "The %d'th signal generator did exit(%d)\n", i, (int)WEXITSTATUS(status));
    else {
      syserror(1, "Trouble shutting down the %d'th signal generator\n", i);
      exit(1);
    }
  }
#endif

  /* and exit cleanly */
  exit(0);
}


/** Simply output version information to stdout */
void
OutputVersion ( void )
{
  printf ( "%s\n", lalGitID );
  printf ( "%s\n", lalappsGitID );

  return;

} /* OutputVersion() */
