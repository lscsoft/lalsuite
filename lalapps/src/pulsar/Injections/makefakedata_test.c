 /*-----------------------------------------------------------------------
 *
 * File Name: makefakedata_v2.c
 *
 * Authors: Papa, M.A., 
 *
 * Revision: $Id$
 *
 * History:   Created by Papa July 2002
 *            Modified by X. Siemens to fix random seed bug
 *
 *            2003/10/03 modified by Bruce Allen for S2 pulsar
 *                       injections
 *           
 *-----------------------------------------------------------------------
 */

/*
 * 1.  An author and Id block
 */

/************************************ <lalVerbatim file="makefakedataCV">
Author: Papa, M.A. 
$Id$
************************************* </lalVerbatim> */


/* ************************************************ <lalLaTeX>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Program \ \texttt{makefakedata.c}}
\label{s:makefakedata.c}
Produces fake SFT data.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Usage}
\begin{verbatim}
makefakedata [-d debuglevel] [-o outfile] [-f f0] [-p alpha delta] [-I input dir]
\end{verbatim}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}

This program uses Teviet Creighton's LAL CW signal routines in order
to produce SFT files in the GEO binary format.


The signal parameters are specified in an input file. The default name
for the input file is In.data, however a file with another name can
also be used, as long as the name is specified with the -i command
line argument. The input file must contain the following information:
\begin{tabular}{|r|r|}
time-baseline of the SFT      &          (Tsft\_in\_sec)\\
number of SFTs to be produced &       (nTsft)\\
\end{tabular}

frequency of first bin of SFTs       (first\_SFT\_frequency\_in\_Hz)
band of produced SFT data            (SFT\_freq\_band\_in\_Hz)
standard deviation of noise 
    (for real and imag) SFT          (std\_of\_noise.When=0\_only\_signal\_present)
amplitude of plus polarization       (Aplus)
amplitude of cross polarization      (Across)
polarization angle                   (psi)
initial phase                        (phi0)
intrinsic emission frequency 
    at the beginning of observation   (f0)
position of source (eq. coordinates)  (latitude\_in\_degrees)
                      "               (longitude\_in\_degrees)
maximum spin-down order               (max\_spin-down\_param\_order)
name of time-stamps file              (name\_of\_time-stamps\_file)
The information in parenthesis above shows what appears in the In.data
input file as a comment to help you remember what the different
entries are.

A number of SFTs will be created. The names will be 
NAME.00001
NAME.00002
.....
and so on.
The default name for the SFT files is TEST\_SFT however a different
name can be specified using the -n command line argument.

How many SFTs will be created is specified in the input file mentioned
above (In.data is the default).

The time of the first datum of each SFT has to be specified. This is
done with a time-stamps file. The time-stamps file must have at least
as many rows as the number of SFTs that one wants to create and it has
two columns: one column for the gps seconds and one for the gps
nanoseconds of each time-stamp. The name of the time-stamps file is
specified in the input file (In.data default name).

If one uses the command line argument -t to specify a filename, say 
TIMEFILES, then a set of files:
TIMEFILES.00001
TIMEFILES.00002
.....
and so on is created. These contain the time series data that each SFT
is computed from. Note that Teviet's routines allow idealized
heterodyning that is used inthis code to reduce the amount of produced
data.

If sigma in the input file is negative, then the noise data is read in
from SFT files that are specified in the code.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Exit codes}
\vspace{0.1in}
\input{MAKEFAKEDATACErrorTable}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Uses}
\begin{verbatim}
LALGenerateTaylorCW()
LALSimulateCoherentGW()
LALMalloc()
LALSCreateVector()
LALSDestroyVector()
LALFree()
LALCheckMemoryLeaks()
\end{verbatim}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\subsubsection*{Notes}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{TestDriveHoughCV}}

********************************************   </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <getopt.h>

NRCSID (MAKEFAKEDATAC, "$Id$");


/* Error codes and messages */


/************** <lalErrTable file="MAKEFAKEDATACErrorTable"> */
#define MAKEFAKEDATAC_ENORM 0
#define MAKEFAKEDATAC_ESUB  1
#define MAKEFAKEDATAC_EARG  2
#define MAKEFAKEDATAC_EBAD  3
#define MAKEFAKEDATAC_EFILE 4

#define MAKEFAKEDATAC_MSGENORM "Normal exit"
#define MAKEFAKEDATAC_MSGESUB  "Subroutine failed"
#define MAKEFAKEDATAC_MSGEARG  "Error parsing arguments"
#define MAKEFAKEDATAC_MSGEBAD  "Bad argument values"
#define MAKEFAKEDATAC_MSGEFILE "Could not create output file"
/******************************************** </lalErrTable> */


/***************************************************/

#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <glob.h>
#include <string.h>


#include <lal/LALConfig.h>
/* #include <lal/LALStdlib.h> */
#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/SeqFactories.h>
#include <lal/Random.h>
#include <lal/RealFFT.h>
#include <lal/LALConstants.h>
#include <lal/DetResponse.h>
#include <lal/DetectorSite.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GenerateTaylorCW.h>
#include <lal/LALDatatypes.h>
#include <lal/LALBarycenter.h> 
#include <lal/LALInitBarycenter.h> 
#include <lal/LALVersion.h>
#include <lal/GenerateSpinOrbitCW.h>

/* Locations of the earth and sun ephemeris data */
#define EARTHDATA "earth00-04.dat"
#define SUNDATA "sun00-04.dat"

int mycalls=0;
int myclears=0;
void *ptr;

/* For debugging memory allocation */
#if(0)
#define LALMalloc(how)    ptr=LALMalloc(how);\
                          mycalls++,\
                          fprintf(stderr,"Allocating %d bytes with %d calls tp %p\n",(INT4)(how),mycalls,ptr)

#define LALFree(quantity) LALFree(ptr=quantity);\
                          myclears++;\
                          fprintf(stderr,"Clearing with %d calls to LALFree at %p\n",myclears,ptr)
#endif

/***************************/
/* Usage format string. */

#define USAGE "Usage: %s [-d debuglevel] [-o outfile] [-f f0] [-p alpha delta]\n"

/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, MAKEFAKEDATAC, statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              MAKEFAKEDATAC, (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( MAKEFAKEDATAC_ESUB, MAKEFAKEDATAC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return MAKEFAKEDATAC_ESUB;                                  \
  }                                                                  \
} while (0)
/******************************************************************/

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

/***************************/

/* Seconds added atthe end and at the beginning of Taylor series*/
#define LTT 1000
/* Maximum # of files in a directory  */
#define MAXFILES 40000        
/* Maximum # of characters of a SFT filename */
#define MAXFILENAMELENGTH 256   

/*Default input data file name*/
const char *inDataFilename="In.data";
char timestampsname[128];

/* non-null values mean to create the files! */
char *freqbasefilename=NULL;
char *timebasefilename=NULL;
char *noisedir;
char filelist[MAXFILES][MAXFILENAMELENGTH];
REAL8 GPStime=-1.0;
INT4 pulsar_defined_at_fiducial_SSB=0;
LIGOTimeGPS SSBpulsarparams;
INT4 do_windowing=0;
INT4 nomagic=0;
INT4 binaryoutput=0;
REAL8 xaxis=0.0;
INT4 doxaxis=0.0;
char *programname=NULL;
char *earthdata;
char *sundata;

/* timebaseline of SFT in sec, band SFT in Hz */
REAL8 Tsft,B,sigma;
/* smallest frequency in the band B */
REAL8 global_fmin;
/* How many SFTS we'll produce*/
INT4 nTsft; 

/* Do you want noise ?*/
INT2 inoise;

/*SCALE factor*/
REAL8 scale = 1.0;

/* Our detector*/
LALDetector Detector;

/* Stores earth/sun ephemeris data for barycentering */
EphemerisData *edat=NULL;
EarthState earth;
EmissionTime emit;

/* Stores detector location and other barycentering data 
 for the first and the follow-up demodulation */
BarycenterInput baryinput;

/* Store time stamps from SFT data */
LIGOTimeGPS *timestamps=NULL;
LIGOTimeGPS SSBfirst,SSBlast;

/*Time series for every SFT*/
REAL4TimeSeries *timeSeries = NULL;

/* Signal parameters to generate signal at source */
TaylorCWParamStruc genTayParams;
CoherentGW cgwOutput;
DetectorResponse cwDetector;

/*This will hold the SFT*/
COMPLEX8Vector *fvec = NULL;
COMPLEX8Vector *fvecn = NULL;

/*FFT plan*/
RealFFTPlan *pfwd = NULL;

INT4 lalDebugLevel = 0;

/* Prototypes for the functions defined in this file */
int read_commandline_and_file(LALStatus *, int argc, char *argv[]);
int set_default(void);
int read_timestamps(REAL8);
int write_modulated_amplitudes_file(LALStatus *);
int prepare_baryinput(LALStatus *);
int prepare_cwDetector(LALStatus *);
int prepare_timeSeries(void);
int compute_SSBtimes(LALStatus *);
int prepare_fvec(LALStatus *);
int make_and_add_time_domain_noise(LALStatus *);
int make_filelist(void);
int read_and_add_freq_domain_noise(LALStatus *, int iSFT);
int add(LALStatus *);
int write_SFTS(int iSFT);
int write_timeseries(int iSFT);
int cleanup(LALStatus *);
int freemem(LALStatus *);
int window_data(void);
int correct_phase(LALStatus *);
void syserror(const char *fmt, ...);
void error(const char *fmt, ...);
void compute_one_SSB(LALStatus* status, LIGOTimeGPS *ssbout, LIGOTimeGPS *gpsin);
INT4 myRound (REAL8 x); 
int parseR4(FILE *fp, const char* vname, REAL4 *data);
int parseR8(FILE *fp, const char* vname, REAL8 *data);
int parseI4(FILE *fp, const char* vname, INT4 *data);
void usage(FILE *fp);

extern void write_timeSeriesR4 (FILE *fp, const REAL4TimeSeries *series);
extern void write_timeSeriesR8 (FILE *fp, const REAL8TimeSeries *series);

/* Like perror() but takes variable numbers of arguments and includes
   program name*/
void syserror(const char *fmt, ...){
  char *thiserror=NULL;
  pid_t pid=getpid();
  va_list ap;
  /* initialize variable argument list */
  va_start(ap,fmt);
  if (errno && (thiserror=strerror(errno)))
    fprintf(stderr,"%s [PID=%d]: %s: ", programname, (int)pid, thiserror);
  else
    fprintf(stderr,"%s [PID=%d]: ", programname, (int)pid);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  return;
}

void error(const char *fmt, ...){
  pid_t pid=getpid();
  va_list ap;
  /* initialize variable argument list  */
  va_start(ap,fmt);
  fprintf(stderr,"%s [PID=%d]: ", programname, pid);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  return;
}

#if 0
void printclock(char *message){

  /* Define this to get crude profiling information */
  char datestring[64];
  time_t t=time(NULL);
  sprintf(datestring, "%s", ctime(&t));
  datestring[strlen(datestring)-1]='\0';
  error("%s : %s\n", datestring, message);
  return;
}
#else
#define printclock(a)
#endif

int main(int argc,char *argv[]) {
  
  static LALStatus status;
  int iSFT;
  SpinOrbitCWParamStruc spinorbit;
  REAL4 dfdt;
  
  programname=argv[0];
  
  /* check that LAL header and library versions are consistent */
  if (
      strcmp(lalVersion,LAL_VERSION) ||
      strcmp(lalConfigureArgs,LAL_CONFIGURE_ARGS) ||
      strcmp(lalConfigureDate,LAL_CONFIGURE_DATE) ||
      fabs(lalVersionMajor-LAL_VERSION_MAJOR)>1.e-3 ||
      fabs(lalVersionMinor-LAL_VERSION_MINOR)>1.e-3
      ) {
    error( "Mismatch between compile time header versions and run-time library version:\n");
    error( "LAL Version: %s\nMajor Version: %d\nMinor Version: %d\nConfig Args: %s\nBuild Date: %s\n",
	    lalVersion,lalVersionMajor,lalVersionMinor,lalConfigureArgs,lalConfigureDate);
    error( "LAL Headers: %s\nMajor Version: %d\nMinor Version: %d\nConfig Args: %s\nBuild Date: %s\n",
	    LAL_VERSION, LAL_VERSION_MAJOR, LAL_VERSION_MINOR, LAL_CONFIGURE_ARGS, LAL_CONFIGURE_DATE);
    exit(1);
  }
  
  printclock("Start");
  
  /* set genTayParams.f = NULL and 60s genTayParams interpolation interval */
  if (set_default())
    return 1;
  
  /* read command line arguments and input variables from input data file */
  if (read_commandline_and_file(&status,argc,argv))
    return 1;
  
  if (sigma < 0.0){
    /* make file list for noise data files */
    if (make_filelist())
      return 1;
  }
  
  /* Read in timestamps and place them in timestamps vec*/
  if (read_timestamps(GPStime))
    return 1;

  printclock("Making baryinput");  

  /* Fill-in edat and baryinput.  Uses alpha,delta from genTayParams. */
  if (prepare_baryinput(&status))
    return 1;
  
  printclock("Computing SSB times");

  /* compute and write SSB times in SSBtimestamps corresponding to
     detector GPS times */
  if (compute_SSBtimes(&status))
    return 1;
  
  /* Fill-in cwDetector */
  if (prepare_cwDetector(&status))
    return 1;

  /* Allocate space for time series that will hold SFT time chunks and
     fill-in all fields that are independent of the particular SFT */
  if (prepare_timeSeries())
    return 1;

  /* Allocate space for SFT vector and fill-in appropriate fields, if
     we will be making SFTs. Also, create plan for forward FFT.  */
  if (freqbasefilename && prepare_fvec(&status))
    return 1;
  
  /* If option active, write amplitude-modulation info file*/
#if(0)
  if (write_modulated_amplitudes_file(&status))
    return 1;
#endif

  printclock("Starting Taylor approximations");
  
  /*
  SSBpulsarparams == time at which pulsar parameters defined
  SSBfirst        == SSB time corresponding to first output sample
  */
  
  memset(&cgwOutput, 0, sizeof(CoherentGW));
  memset(&spinorbit, '\0', sizeof(spinorbit));
  
  /* The GENERATE routines work entirely in Barycentric Time. This is
     the fiducial time at which pulsar parameters are defined */
  spinorbit.spinEpoch=SSBpulsarparams;
  
  /* print out the SSB time corresponding to first detector GPS,
     useful for analysis and checking */
  /*
  error("\n"
	"    Starting time is GPS  = %d.%09d\n"
	"    Params defined at SSB = %d.%09d by file %s\n",
	timestamps[0].gpsSeconds, timestamps[0].gpsNanoSeconds,
	spinorbit.spinEpoch.gpsSeconds, spinorbit.spinEpoch.gpsNanoSeconds, inDataFilename);
  */
  /* These define the start and end of an interpolation table.  Set
     epoch at least the light travel time EARLIER than the SSB time of
     the first timestamp, so that the interpolation table goes beyond the
     times that we need at the begining */
  spinorbit.epoch    = SSBfirst;
  spinorbit.epoch.gpsSeconds -= 0.75*LTT;
  
  /* Set the interpolation table length large enough that the
     interpolation table extends beyond the end time. */
  spinorbit.length   = (1.5*LTT + Tsft + SSBlast.gpsSeconds -
			SSBfirst.gpsSeconds)/genTayParams.deltaT;
  
  /* These quantities come straight from user input */
  spinorbit.deltaT   = genTayParams.deltaT;
  spinorbit.aPlus    = genTayParams.aPlus;
  spinorbit.aCross   = genTayParams.aCross;
  spinorbit.phi0     = genTayParams.phi0;
  spinorbit.f0       = genTayParams.f0;
  spinorbit.position = genTayParams.position;
  spinorbit.psi      = genTayParams.psi;
  /* a copy of the pointer to fdot values */
  spinorbit.f        = genTayParams.f;
  
  /* the following input fields are NOT used when calling  LALGenerateSpinOrbitCW() */
  /* spinorbit.orbitEpoch */
  /* spinorbit.omega */
  /* spinorbit.oneMinuteEcc */
  /* spinorbit.angularSpeed */
  
  /* this is how we select to just call LALGenerateTaylorCW */
  spinorbit.rPeriNorm = 0.0; 

#if 0
  /*  Reproduce the error that was in makefakedata_v2 until around Oct
  7th 2003 when corrected by Bruce.  This should normally NOT be
  enabled. */
  error( 
	  "WARNING: reproducing error from < Oct 7, 2003!\n"
	  "This error can be seen in the time domain.  The first few\n"
	  "seconds of output are zero.  This is because the interpolation\n"
	  "table used within the code doesn't extend far enough to the\n"
	  "past and future.\n"
	  );
  spinorbit.length=(LTT+Tsft+timestamps[nTsft-1].gpsSeconds-timestamps[0].gpsSeconds)/genTayParams.deltaT;
  spinorbit.epoch =SSBfirst;
#endif
  
  /* Now correct the pulsar parameters, if needed, for the Barycentric
     epoch for which they are defined  */
  SUB(LALGenerateSpinOrbitCW(&status, &cgwOutput, &spinorbit), &status);
  dfdt=spinorbit.dfdt;
  
  /* check that output is OK */
  if (dfdt > 2.0) {
    error("Waveform sampling interval is too large:\n"
	  "\tmaximum df*dt = %f\n", dfdt );
    return 1;
  }

  if (lalDebugLevel >= 3)
    {
      FILE *fp;
      fp = fopen ("debug_phi_v2.dat", "w");
      write_timeSeriesR8 (fp, cgwOutput.phi);
      fclose (fp);
    }


  /* This is the main loop that produces output data */
  for (iSFT=0;iSFT<nTsft;iSFT++){

    /* This sets the time at which the output is given...*/
    timeSeries->epoch=timestamps[iSFT];
    
    /* Note that we DON'T update cwDetector Heterodyne Epoch. Teviet
    says: "You can set it to anything you like; it doesn't really
    matter as long as it's the same from one stretch of simulated data
    to the next."  See above for more remarks about this. */

    printclock("Starting simulate coherent");

    /* produce a time series simulation of a CW signal */
    SUB( LALSimulateCoherentGW(&status, timeSeries, &cgwOutput, &cwDetector), &status);

    if (lalDebugLevel >= 3)
      {  
	FILE *fp;
	CHAR fname[512];
	sprintf (fname, "Tseries_v2_%05d.dat", iSFT);
	fp = fopen (fname, "w");
	write_timeSeriesR4 (fp, timeSeries);
	fclose (fp);
      }

    /*if you want noise, make it and add to timeseries */
    if (sigma > 0.0 && make_and_add_time_domain_noise(&status))
      return 1;
    
    /* Write Time-domain strain file in the file specified by the path
       name. If no path is set, this means don't output in Time domain */
    if (write_timeseries(iSFT))
      return 1;

    /* if we've asked for freq-domain output... */
    if (freqbasefilename) {
      
      /* window data in the time domain before FFTing it */
      if (do_windowing && window_data())
	return 1;

      /* Perform FFTW-LAL Fast Fourier Transform */
      SUB(LALForwardRealFFT(&status, fvec, timeSeries->data, pfwd), &status);

#if 1
      /* correct phase */
      correct_phase(&status);
#endif
      
      /* if you want noise added in the FREQ domain only, read from files and add in */
      if (sigma < 0.0 && read_and_add_freq_domain_noise(&status, iSFT))
	return 1;
      
      /* Write the SFT file in the file that is specified by the path
	 name.  If no path name is set, this means don't output SFTs*/
      if (write_SFTS(iSFT))
	return 1;  
    }

  } /* end of loop over different SFTs */
  
  if (cleanup(&status))
    return 1;
  
  if (freemem(&status))
    return 1;

   LALCheckMemoryLeaks(); 

/*   INFO(MAKEFAKEDATAC_MSGENORM); */

/*   return MAKEFAKEDATAC_ENORM; */

  return 0;
}







/*cleans up before computing another SFT*/
int cleanup(LALStatus* status){

  LALSDestroyVectorSequence( status, &(cgwOutput.a->data));
  LALSDestroyVector( status, &( cgwOutput.f->data ) );
  LALDDestroyVector( status, &( cgwOutput.phi->data ) );
  
  LALFree(cgwOutput.a);
  LALFree(cgwOutput.f);
  LALFree(cgwOutput.phi);
  

  return 0;

}


/* This routine frees up all the memory */
int freemem(LALStatus* status){


  LALSDestroyVector( status, &(timeSeries->data) );
  LALFree(timeSeries);
  /* Clean up earth/sun Ephemeris tables */
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);

  /* Clean up cwDetector */
  LALCDestroyVector( status, &( cwDetector.transfer->data ) );
  LALFree(cwDetector.transfer);
/*   LALFree( cwDetector.ephemerides->ephemE->pos ); */
/*   LALFree( cwDetector.ephemerides->ephemS ); */
/*   LALFree( cwDetector.ephemerides ); */
/*   LALFree( Detector); */
/*   LALFree( cwDetector.site ); */

  /* Clean up timestamps */
  LALFree(timestamps); 

  /* Clean up FFTs of signal and of noise - each a complex8vector */
  if (fvec)
    LALCDestroyVector(status, &fvec);

  if (pfwd)
    LALDestroyRealFFTPlan(status, &pfwd);

  if (genTayParams.f)
    LALDDestroyVector( status, &(genTayParams.f) );

  if (fvecn)
    LALCDestroyVector(status, &fvecn);  

  return 0;
}


/*sets default values */
int set_default(void) {

  genTayParams.f=NULL;
  genTayParams.position.system= COORDINATESYSTEM_EQUATORIAL;
  genTayParams.deltaT=60;

  return 0;
}


void compute_one_SSB(LALStatus* status, LIGOTimeGPS *ssbout, LIGOTimeGPS *gpsin) 
{

  LALBarycenterEarth(status, &earth, gpsin, edat);
  baryinput.tgps = *gpsin;
  LALBarycenter(status, &emit, &baryinput, &earth);
  
  *ssbout = emit.te;
  return;
}

/* Computes SSB times corresponding to given GPS times on the earth */
int compute_SSBtimes(LALStatus* status) {
  compute_one_SSB(status, &SSBfirst, &timestamps[0]);
  compute_one_SSB(status, &SSBlast,  &timestamps[nTsft-1]);
  
  /* 
     if user has not defined the epoch at which the pulsar parameters
     are defined, then choose it to be the SSB time of the first
     timestamps:
  */

  if (!pulsar_defined_at_fiducial_SSB)
    SSBpulsarparams=SSBfirst;
  
  return 0;
}

#if 0
/*applies correct heterodyning from one SFT to the next*/
int correct_phase(LALStatus* status, int iSFT) {

  int i;
  REAL8 cosx,sinx,x;
  COMPLEX8 fvec1;

  x=2.0*LAL_PI*global_fmin*Tsft*iSFT;
  cosx=cos(x);
  sinx=sin(x);
  for (i = 0; i < fvec->length; ++i){
    fvec1=fvec->data[i];
    fvec->data[i].re=fvec1.re*cosx+fvec1.im*sinx;
    fvec->data[i].im=fvec1.im*cosx-fvec1.re*sinx;
  }
  
  return 0;
}

#endif



/* correct the heterodyned phase */
/*
 * "If the starting (which is also the heterodyning) frequency of the 
 * In.data file does not go through an integer number of cycles in the 
 * time corresponding to the gap, then the phase of the signal gets 
 * screwed up by the heterodyning procedure.", according to Xavier.
 *
 * This makes the resulting F stat values different from the correct 
 * ones by 20,000% in the worst case I have seen.  
 *
 * The correct_phase() solves this problem as far as F stat 
 * values are concerned. 
 *
 * Detail of validation:
 * Performing 60,000 MC experiments over differnt signal parameters 
 * (alpha,delta,phi0,psi,cosiota,fsignal, starting_frequency_in_In.data) 
 * for 10 hours-worth observation time with gaps for LLO, the resulting 
 * F stat values differ from the correct ones by up to 0.01% and for 99% 
 * of the time of the experiments less than 0.003%. 
 *
 * When there is no gap, the differences are less than 0.004% 
 * among 60,000 MCs. 
 * Differences between with this routine and without this routine are at most 
 * 6e-7% out of 120,000 MCs. 
 * Yousuke 24 Mar 2004.
 *
 */
int correct_phase(LALStatus* status) {

  UINT4 i;
  REAL8 cosx,sinx;
  COMPLEX8 fvec1;
  LALTimeInterval deltaGPS;
  LIGOTimeGPS gps1,gps2;
  REAL8 deltaT;

  gps1=timeSeries->epoch;
  gps2=cwDetector.heterodyneEpoch;

  LALDeltaGPS(status,&deltaGPS,&gps1,&gps2);

  LALIntervalToFloat(status,&deltaT,&deltaGPS);

  deltaT *= LAL_TWOPI*timeSeries->f0; 

  cosx=cos(deltaT);
  sinx=sin(deltaT);
  for (i = 0; i < fvec->length; ++i){
    fvec1=fvec->data[i];
    fvec->data[i].re=fvec1.re*cosx-fvec1.im*sinx;
    fvec->data[i].im=fvec1.im*cosx+fvec1.re*sinx;
  }
  
  return 0;
}


/* windows the data*/
int window_data(void){

  REAL4 *window,timedata1,timedata2;
  REAL4 frac;
  REAL8 x;
  INT4 nbinw,i;
  UINT4 data_length;

  frac=0.05; /*this is half the fraction of modified data points*/
  
  data_length=timeSeries->data->length;
  nbinw=0.5+data_length*frac;

  window=(REAL4 *)LALMalloc(nbinw*2*sizeof(REAL4)); 

  if (nbinw < 10) {
        error( "Not enough data points to window !\n");
	return 1;
  }
  
  for (i = 0; i < nbinw; ++i){
    x=-LAL_PI+i*(LAL_PI/nbinw);
    window[i]=(cos(x)+1)/2.0;
  }

  
  for (i = 0; i < nbinw; ++i){
    timedata1=timeSeries->data->data[i];
    timedata2=timeSeries->data->data[data_length-1-i];
    
    timeSeries->data->data[i]=timedata1*window[i];
    timeSeries->data->data[data_length-1-i]=timedata2*window[i];
  }

  

  LALFree(window);

  return 0;
}



/* Note (Bruce Allen).  If the pulsar parameters are being specified
   with the -S SSBpulsarparams option, then it may be necessary to put in a
   DIFFERENT sky position in the routine that follows, if the pulsar
   has a significant proper motion (nonzero velocity of pulsar) that
   carries it to a new sky position between the fiducial SSB time and
   time of the simulated data. */

/* Sets up edat and baryinput: reads ephemeris data files and fills-in
   baryinput fields */
int prepare_baryinput(LALStatus* status){
  
  /* Quantities computed for barycentering */
  edat=(EphemerisData *)LALMalloc(sizeof(EphemerisData));
  (*edat).ephiles.earthEphemeris = earthdata;
  (*edat).ephiles.sunEphemeris =   sundata;
  (*edat).leap=13; 
  
  /* Read in ephemerides */  
  LALInitBarycenter(status, edat);
  
  /* Getting detector coords from DetectorSite module of tools package */     
  baryinput.site=Detector;
  baryinput.site.location[0]=Detector.location[0]/LAL_C_SI;
  baryinput.site.location[1]=Detector.location[1]/LAL_C_SI;
  baryinput.site.location[2]=Detector.location[2]/LAL_C_SI;
  baryinput.alpha=genTayParams.position.longitude;
  baryinput.delta=genTayParams.position.latitude;
  baryinput.dInv=0.e0;

  return 0;
}

/* prepares cwDetector */
int prepare_cwDetector(LALStatus* status){

  memset(&cwDetector, 0, sizeof(DetectorResponse));
  /* The ephemerides */
  cwDetector.ephemerides = edat;
  /* Specifying the detector site (set above) */
  cwDetector.site = &Detector;  
  /* The transfer function.  
   * Note, this xfer function has only two points at it extends 
   between 0 and 16384 Hz. The routine that will generate the signal as 
   output from the detector on Earth will interpolate*/
  cwDetector.transfer = (COMPLEX8FrequencySeries *)LALMalloc(sizeof(COMPLEX8FrequencySeries));
  memset(cwDetector.transfer, 0, sizeof(COMPLEX8FrequencySeries));
  /* it does not change so just use first timestamp. Does not
   seem to matter whether SSBtimestamps or timestamps are used */
  cwDetector.transfer->epoch = timestamps[0]; 
  cwDetector.transfer->f0 = 0.0;
  cwDetector.transfer->deltaF = 16384.0;
  cwDetector.transfer->data = NULL;
  LALCCreateVector(status, &(cwDetector.transfer->data), 2);
  
  /* unit response function */
  cwDetector.transfer->data->data[0].re = 1.0;
  cwDetector.transfer->data->data[1].re = 1.0;
  cwDetector.transfer->data->data[0].im = 0.0;
  cwDetector.transfer->data->data[1].im = 0.0;
  
  /* 
     Note that we DON'T update cwDetector Heterodyne Epoch.  Teviet
     says: "You can set it to anything you like; it doesn't really
     matter as long as it's the same from one stretch of simulated
     data to the next."
     
     For this reason, we try and set it in a way that will be
     consistent if the user run makefakedata several times, to make
     different stretches of data from 'the same' source.  The idea is
     that, if they do this, they will be using the '-S' option to
     define an SSB time at which the pulsar parameters are defined.
     So we compute, at the given detector, the GPS time corresponding
     to that SSB time, and then use THAT GPS time to define the
     heterodyneEpoch.  If you don't like it, complain to Bruce.
  */

  if (!pulsar_defined_at_fiducial_SSB)
    cwDetector.heterodyneEpoch=timestamps[0];
  else 
    cwDetector.heterodyneEpoch=SSBpulsarparams;
#if 0 
{

    /* Find GPS detector time corresponding to SSBpulsarparams. To
       start root finding, use SSBpulsarparams as guess (not off by
       more than 400 secs! */
    
    LIGOTimeGPS SSBofguess, GPSguess=SSBpulsarparams;
    INT4 iterations, E9=1000000000;
    INT8 delta, guess;

    /* now find GPS time corresponding to SSBpulsarparams */
    for (iterations=0; iterations<100; iterations++) {

      /* find SSB time of guess */
      compute_one_SSB(status, &SSBofguess, &GPSguess);
      
      /* compute difference between that and what we want.  Be careful
	 with operations in INT4s. They will overflow if you are not
	 careful!  */
      delta  = SSBpulsarparams.gpsSeconds;
      delta -= SSBofguess.gpsSeconds;
      delta *= E9;
      delta += SSBpulsarparams.gpsNanoSeconds;
      delta -= SSBofguess.gpsNanoSeconds;
      
      /* break if we've converged */
      /*      if (delta>-2 && delta<2)*/
      if (delta == 0)	/* try to be ns-precise */
	break;

      /* use delta to make next guess.  Be careful of the order of
	 operations here. Remember that these quantities overflow
	 INT4s. */
      guess  = GPSguess.gpsSeconds;
      guess *= E9;
      guess += GPSguess.gpsNanoSeconds;
      guess += delta;

      /* from here on use delta as a temporary variable */
      GPSguess.gpsSeconds = delta = guess/E9;
      delta *= E9;
      guess -= delta;
      GPSguess.gpsNanoSeconds=guess;
    }

    /* check for convergence of root finder */
    if (iterations==100){
      error("Computation of GPS time for heterodyne epoch did not converge!\n");
      return 1;
    }

    /* Now that we've found the GPS time that corresponds to the SSB
       time at which the user specified the pulsar's parameters,
       please use that to heterodyne by. */
    cwDetector.heterodyneEpoch=GPSguess;
  }
#endif
  
  return 0;
}

/*Allocates space for timeseries */
int prepare_timeSeries(void) {

  timeSeries = (REAL4TimeSeries *)LALMalloc(sizeof(REAL4TimeSeries));
  timeSeries->data = (REAL4Vector *)LALMalloc(sizeof(REAL4Vector));
  
  timeSeries->data->length = B*Tsft*2;


  timeSeries->data->data = (REAL4 *)LALMalloc(timeSeries->data->length*sizeof(REAL4));
  timeSeries->deltaT=Tsft/timeSeries->data->length;

  timeSeries->f0=global_fmin;

  return 0;
}


/*Allocates space for timeseries */
int prepare_fvec(LALStatus* status) {

  INT4 len=timeSeries->data->length;
  INT4 len2=len/2+1;
  
  /* Create vector to hold signal frequency series */
  LALCCreateVector(status, &fvec, (UINT4)len2);
    
  /* Compute measured plan for FFTW */
  LALCreateForwardRealFFTPlan(status, &pfwd, (UINT4)len, 0);
  
  return 0;
}


int make_and_add_time_domain_noise(LALStatus* status) {


  REAL4Vector    *v1=NULL;
  RandomParams   *randpar=NULL;
  INT4            i,numPoints,seed,errorcode;
  FILE *devrandom;

  numPoints=timeSeries->data->length;


  LALCreateVector (status, &v1, numPoints);
  
  /*
   * Modified so as to not create random number parameters with seed
   * drawn from clock.  Seconds don't change fast enough and sft's
   * look alike.  We open /dev/urandom and read a 4 byte integer from
   * it and use that as our seed.  Note: /dev/random is slow after the
   * first, few accesses.
   */
  
  if (!(devrandom=fopen("/dev/urandom","r"))){
    syserror("Unable to open device /dev/urandom\n");
    return 1;
  }
  errorcode=fread((void*)&seed,sizeof(INT4),1,devrandom);
  if (errorcode!=1)
    {
    syserror( "Error reading /dev/random file!\n");
    return 1;
    }
  fclose(devrandom);
  
  LALCreateRandomParams (status, &randpar, seed);

  LALNormalDeviates(status, v1, randpar);

  for (i = 0; i < numPoints; ++i)
    timeSeries->data->data[i]+=sigma*v1->data[i];


  /* destroy randpar*/
  LALDestroyRandomParams (status, &randpar);
  
  /*   destroy v1 */
  LALDestroyVector (status, &v1);

  return 0;

}

/*reads timestamps file and fills-in timestamps vector*/
int read_timestamps(REAL8 startattime) {
  
  FILE *fp;
  int i,r;
 
  timestamps=(LIGOTimeGPS *)LALMalloc(nTsft*sizeof(LIGOTimeGPS)); 

  if (startattime==-1.0) {
    /* Read */
    if (!(fp=fopen(timestampsname,"r"))) {
      syserror("Unable to find timestampsname file %s\n",timestampsname);
      return 1;
    }
    
    for (i=0;i<nTsft;i++){
      if (2!=(r=fscanf(fp,"%d  %d\n", &timestamps[i].gpsSeconds, &timestamps[i].gpsNanoSeconds))){
	syserror("Unable to read datum from line # %d from file %s\n", i+1, timestampsname);
	return 1; 
      } 
    }  
    fclose(fp);
  }
  else {
    REAL8 time0 = startattime;
    REAL8 frac=0.0;
    
    /* set up array based on timestamps implied */    
    for (i=0; i<nTsft; i++){
      frac=time0 -(int)time0;
      timestamps[i].gpsSeconds=(int)time0;
      timestamps[i].gpsNanoSeconds=1000000000*frac;
      time0 += Tsft;
    }
  }
  return 0;
}

int write_modulated_amplitudes_file(LALStatus* status){
  FILE *fp;
  LALDetAMResponse  amresp;
  LALSource         source;
  LALDetAndSource   detectorandsource;
  LIGOTimeGPS       gps;
  LALGPSandAcc      gpsandacc;
  const char *filename="AmplMod.dat";
  int i;

  
  if (!(fp=fopen(filename,"w"))) {
    syserror("Unable to find file %s in write_modulated_amplitudes_file routine\n",filename);
    return 1;
  }

  detectorandsource.pDetector=&Detector;
  source.equatorialCoords=genTayParams.position;
  source.orientation=genTayParams.psi;  
  detectorandsource.pSource=&source;

  for (i=0;i<=nTsft;i++){
  
    gps.gpsSeconds=timestamps[i].gpsSeconds;
    gps.gpsNanoSeconds=timestamps[i].gpsNanoSeconds;
    gpsandacc.gps=gps;
    gpsandacc.accuracy=LALLEAPSEC_STRICT;

    LALComputeDetAMResponse(status, &amresp, &detectorandsource, &gpsandacc);
    fprintf(fp,"%f  %f\n",amresp.plus,amresp.cross);
    
  }
  
  fclose(fp);

  return 0;
}




int make_filelist(void) {

  UINT4 fileno=0;
  char command[256];
  glob_t globbuf;

  strcpy(command,noisedir);
  strcat(command,"*SFT*");
  
  globbuf.gl_offs = 1;
  glob(command, GLOB_ERR, NULL, &globbuf);

  /* read file names -- MUST NOT FORGET TO PUT ERROR CHECKING IN HERE !!!! */
  while (fileno< globbuf.gl_pathc) 
    {
      strcpy(filelist[fileno],globbuf.gl_pathv[fileno]);
      fileno++;
      if (fileno > MAXFILES)
	{
	  error("Too many files in directory! Exiting... \n");
	  return 1;
	}
    }
  globfree(&globbuf);


  return 0;

}





int read_and_add_freq_domain_noise(LALStatus* status, int iSFT) {

  FILE *fp;
  size_t errorcode;
  UINT4 i;
  REAL4 norm;

  struct headertag {
    REAL8 endian;
    INT4  gps_sec;
    INT4  gps_nsec;
    REAL8 tbase;
    INT4  firstfreqindex;
    INT4  nsamples;
  } header;

  /**********************/


  /* open FIRST file and get info from it*/
  if (!(fp=fopen(filelist[iSFT],"r"))){
    syserror("Unable to open the fisrt SFT file file %s\n", filelist[iSFT]);
    return 1;
  }
  /* read in the header from the file */
  errorcode=fread((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1) 
    {
      error("No header in data file %s\n",filelist[iSFT]);
      return 1;
    }

  /* check that data is correct endian order */
  if (header.endian!=1.0)
    {
      error("First object in file %s is not (double)1.0!\n",filelist[iSFT]);
      error("It could be a file format error (big/little\n");
      error("endian) or the file might be corrupted\n\n");
      return 2;
    }


  /* check frequency range */
  if ((INT4)(global_fmin*Tsft+0.5) < header.firstfreqindex){
    error("Frequency band of noise data out of range !\n");
    return 1;
  }

  if ( (INT4)((global_fmin*Tsft+0.5)+fvec->length-1) > (header.firstfreqindex + header.nsamples)){
    error("Frequency band of noise data out of range !\n");
    return 1;
  }

  /* seek to position */
  if (0 != fseek(fp, myRound((global_fmin*Tsft - header.firstfreqindex) * 4.0 * 2.0), SEEK_CUR)){
    syserror("file too short (could'n fssek to position\n");
    return 1;
  }


  /* calculate size of the data in bytes */
  /* datasize = (f2ind - f1ind) * head.dsize * 2;*/


  /* allocate storage space if needed */
  if (!fvecn) {
    INT4 len=timeSeries->data->length;
    INT4 len2=len/2+1;
    LALCCreateVector(status, &fvecn, (UINT4)len2);
  }

  /* read the data */
  if (1 != fread(fvecn->data,(fvecn->length)* sizeof(*fvecn->data), 1, fp)) {
      syserror("Could not read the data \n");
    return 1;
  }

  fclose(fp);

  norm=((REAL4)(fvec->length)/((REAL4)header.nsamples));

  for (i = 0; i < fvec->length; ++i) {
    fvec->data[i].re += norm * fvecn->data[i].re;
    fvec->data[i].im += norm * fvecn->data[i].im;
  }
  
  return 0;
}


/*-----------------------------------------------------------------*/
/* myRound function.                                              
   This function returns the nearest integer to x; 
   myRound(x) = sign(x) * {floor(|x|)+step(|x|-floor(|x|)-0.5)}
   step(x) = 1 if x>=0
             0 if x<0   
   Temporal replacement.  The C99 "round" is unavailable now...
   Yousuke, 07-Mar-2004 
*/
/*-----------------------------------------------------------------*/
INT4 myRound(REAL8 x)
{
  REAL8 sign=1.0;
  REAL8 roundedValue=0.0;
  REAL8 remainder=0.0;

  if(x<0) sign=-1.0;
  roundedValue= floor(sign*x);
  remainder=sign*x-roundedValue;
  if(remainder>=0.5) 
    roundedValue=roundedValue+1.0;
  roundedValue=sign*roundedValue;

  return (INT4) roundedValue;
}
/*-----------------------------------------------------------------*/
/*                                                                 */
/*-----------------------------------------------------------------*/







int write_SFTS(int iSFT){


  FILE *fp;
  REAL4 rpw,ipw;
  /* REAL8 fr; */
  char filename[256], filenumber[16];
  int errorcode;
  UINT4 i;

  struct headertag {
    REAL8 endian;
    INT4  gps_sec;
    INT4  gps_nsec;
    REAL8 tbase;
    INT4  firstfreqindex;
    INT4  nsamples;
  } header;

  /* Open SFT data file */
  strcpy(filename,freqbasefilename);
  sprintf(filenumber,".%05d",iSFT); /*ACHTUNG: used to be iSFT+1 - now starts from .00000*/
  strcat(filename,filenumber);
  fp=fopen(filename,"w");
  if (fp==NULL) {
    syserror("Unable to find the SFT data file %s\n",filename);
    return 1;
  }

  header.endian=1.0;
  header.gps_sec=timestamps[iSFT].gpsSeconds;
  header.gps_nsec=timestamps[iSFT].gpsNanoSeconds;
  header.tbase=Tsft;
  header.firstfreqindex=(INT4)(global_fmin*Tsft+0.5);
  header.nsamples=fvec->length;
  
  /* write header */
  errorcode=fwrite((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1){
    syserror( "Error in writing header into file!\n");
    return 1;
  }

  for (i=0;i<fvec->length;i++){

    rpw = fvec->data[i].re;
    ipw = fvec->data[i].im;

    errorcode=fwrite((void*)&rpw, sizeof(REAL4),1,fp);  
    if (errorcode!=1){
      syserror( "Error in writing data into SFT file!\n");
      return 1;
    }
    errorcode=fwrite((void*)&ipw, sizeof(REAL4),1,fp);  
    if (errorcode!=1){
      syserror( "Error in writing data into SFT file!\n");
      return 1;
    }
        
  }
  
  fclose(fp);
  return 0;  
  
}

/* This writes out the simulated dat in the time-domain. */
int write_timeseries(int iSFT){
  
  /* write binary form of the output to stdout.  This is useful for
     hardware pulsar injections at the sites */
  if (binaryoutput){
    REAL4 magic=1234.5;
    UINT4  length=timeSeries->data->length;
    REAL4 *datap=timeSeries->data->data;
    
    if (
	(!nomagic && 1 != fwrite(&magic, sizeof(magic), 1, stdout))
	||
	length!=fwrite(datap, sizeof(datap[0]), length, stdout)
	) {
      syserror( "Fatal error in writing binary data to stdout\n");
      exit(1);
    }

    /* only print the magic number at the start of the stream */
    nomagic=1;
  }
  
  /* write data to a file.  This is useful for debugging, and lots of
     other things too! */
  if (timebasefilename) {
    FILE *fp;
    REAL4 pw;
    char filename[256], filenumber[16];
    int i;
    
    struct headertag {
      REAL8 endian;
      INT4  gps_sec;
      INT4  gps_nsec;
      REAL8 tbase;
      INT4  firstfreqindex;
      INT4  nsamples;
    } header;
    
    /* Open time-domain data file */
    strcpy(filename,timebasefilename);
    sprintf(filenumber,".%05d",iSFT);
    strcat(filename,filenumber);
    fp=fopen(filename,"w");
    if (fp==NULL) {
      syserror("Unable to open file %s\n",filename);
      return 1;
    }
    header.endian=1.0;
    header.gps_sec=timestamps[iSFT].gpsSeconds;
    header.gps_nsec=timestamps[iSFT].gpsNanoSeconds;
    header.tbase=Tsft;
    header.firstfreqindex=(INT4)(global_fmin*Tsft);
    header.nsamples=timeSeries->data->length;
    
    /* write header */
#if(0)
    fprintf(fp, "%e\n",header.endian);
#endif
    
    /* now print data to a file in either one or two column format */
    if (!doxaxis) {
      for (i=0;i<header.nsamples;i++) {    
	pw=timeSeries->data->data[i];
	fprintf(fp,"%f\n",pw);
      }
    }
    else {
      for (i=0;i<header.nsamples;i++) {
	REAL8 ts=header.gps_sec-xaxis+header.tbase*i/header.nsamples;
	pw=timeSeries->data->data[i];
	fprintf(fp,"%f %f\n",ts, pw);
      }
    }
    
    fclose(fp);
  }
  
  return 0;  
}

int parseR4(FILE *fp, const char* vname, REAL4 *data){
  char junk[1024], junk2[1024];
  char command[1024];
  int r;
  
  memset(junk, '\0', 1024);
  memset(junk2,'\0', 1024);
  
  r=fscanf(fp, "%f%[\t ]%[^\012]", data, junk, junk2);
  if (r!=3)  {
    error("Unable to assign %s from file: %s\n"
	  "The entry must be of the form:\n"
	  "NUMBER TEXT\n"
	  "with white space in between. TEXT is NOT optional!\n",
	  vname, inDataFilename);
    sprintf(command, "cat %s 1>&2\n", inDataFilename);
    system(command);
    return 1;
  }
      return 0;
}

int parseR8(FILE *fp, const char* vname, REAL8 *data){
  char junk[1024], junk2[1024];
  char command[1024];
  int r;
  
  memset(junk, '\0', 1024);
  memset(junk2,'\0', 1024);
  
  r=fscanf(fp, "%lf%[\t ]%[^\n]", data, junk, junk2);
  if (r!=3)  {
    error("Unable to assign %s from file: %s\n"
	  "The entry must be of the form:\n"
	  "NUMBER TEXT\n"
	  "with white space in between. TEXT is NOT optional!\n",
	  vname, inDataFilename);
    sprintf(command, "cat %s 1>&2\n", inDataFilename);
    system(command);
    return 1;
  }
      return 0;
}
int parseI4(FILE *fp, const char* vname, INT4 *data){
  char junk[1024], junk2[1024];
  char command[1024];
  int r;
  
  memset(junk, '\0', 1024);
  memset(junk2,'\0', 1024);
  
  r=fscanf(fp, "%d%[\t ]%[^\n]", data, junk, junk2);
  if (r!=3)  {
    error("Unable to assign %s from file: %s\n"
	  "The entry must be of the form:\n"
	  "NUMBER TEXT\n"
	  "with white space in between. TEXT is NOT optional!\n",
	  vname, inDataFilename);
    sprintf(command, "cat %s 1>&2\n", inDataFilename);
    system(command);
    return 1;
  }
      return 0;
}



/* print usage/help message */
void usage(FILE *fp){
  fprintf(fp,"Recognized arguments [defaults] are:\n"
	  "-i Character String           Name of input parameter file               [In.data]\n"
	  "-n Character String           Basefilename of output SFT files           [Makes no SFTs]\n"
	  "-t Character String           Basefilename of output STRAIN files        [Makes no Strains]\n"
	  "-I Character String           LHO, LLO, VIRGO, GEO, TAMA, CIT, ROME      [No default]\n"
	  "-G Double-precision number    Detector GPS time to start data            [Use Timestamp files]\n"
	  "-S Double-precision number    SSB fiducial time at which pulsar defined  [Use Detector's first timestamp]\n"
	  "-X Double-precision number    Include time (minus arg) in STRAIN file    [No column of times]\n"
	  "-E Character String           Directory path for ephemeris files         [./ (current directory)]\n"
	  "-D Character String           Input noise dir                            [/sft/S2-LIGO/S2_H1_FunkyCal30MinSFTs/]\n"
	  "-w                            Window data in time domain before doing FFT\n"
	  "-b                            Output time-domain data in IEEE754 binary format\n"
	  "-m                            DON'T output 1234.5 before time-domain binary samples\n"
	  "-h                            Print this help/usage message\n"
	  );
  return;
}

int read_commandline_and_file(LALStatus* status, int argc,char *argv[]) {
  
  char dmp[128];
  int c, errflg = 0;
  int r,i,msp;
  UINT4 imin, nsamples;  
  FILE *fp;
  char *endptr;
  int detectorset=0;
  REAL8 temptime;
  
  /* scan through the list of arguments on the command line 
     and get the input data filename*/
  
  opterr=0;

  while (!errflg && ((c = getopt(argc, argv,":i:n:t:I:G:S:X:E:D:wbmhd:"))!=-1))
    switch (c) {
    case 'i':
      /* Name of input data file */
      inDataFilename=optarg;
      break;
    case 'n':
      /* Name of SFT data file */
      freqbasefilename=optarg;
      break;
    case 't':
      /* Name of TDD file */
      timebasefilename=optarg;
      break;
    case 'I':
      /* detector site */
      if (!strcmp(optarg,"LHO")) {
	Detector=lalCachedDetectors[LALDetectorIndexLHODIFF];
      } else if (!strcmp(optarg,"LLO")) {
	Detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
      } else if (!strcmp(optarg,"VIRGO")) {
	Detector=lalCachedDetectors[LALDetectorIndexVIRGODIFF];
      } else if (!strcmp(optarg,"GEO")) {
	Detector=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
      } else if (!strcmp(optarg,"TAMA")) {
	Detector=lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
      } else if (!strcmp(optarg,"CIT")) {
	Detector=lalCachedDetectors[LALDetectorIndexCIT40DIFF];
      } else if (!strcmp(optarg,"ROME")) {
	LALFrDetector detector_params;
	LALDetectorType bar;
	LALDetector Detector1;
	
	bar=LALDETECTORTYPE_CYLBAR;
	strcpy(detector_params.name,"NAUTILUS");
	detector_params.vertexLongitudeRadians=12.67*LAL_PI/180.0;
	detector_params.vertexLatitudeRadians=41.82*LAL_PI/180.0;
	detector_params.vertexElevation=300.0;
	detector_params.xArmAltitudeRadians=0.0;
	detector_params.xArmAzimuthRadians=44.0*LAL_PI/180.0;
	
	LALCreateDetector(status,&Detector1,&detector_params,bar);
	
	Detector=Detector1;
      } else {
	error( 
		"Invalid detector choice: -I %s\n"
		"Allowed choices are: LHO, LLO, VIRGO, GEO, TAMA, CIT, ROME\n", optarg);
	return 1;
      }
      detectorset=1;
      break;
    case 'G':
      /* GPS starting time (don't read from file) */
      GPStime=strtod(optarg, NULL);
      if (GPStime<0.0 || GPStime>1.e9){
	error( "GPS time argument to -G = %f must be non-negative and <= 1.e9\n", (double)GPStime);
	exit(1);
      }
      break;
    case 'S':
      /* SSB time at which pulsar parameters defined */
      pulsar_defined_at_fiducial_SSB=1;
      temptime=strtod(optarg, NULL);
      if (temptime<0.0 || temptime>1.e9){
	error( "SSB time argument to -S = %f must be non-negative and <= 1.e-9\n", (double)temptime);
	exit(1);
      }
      LALFloatToGPS(status, &SSBpulsarparams, &temptime);
      break;
    case 'X':
      /* include x axis in time strain output */
      doxaxis=1;
      xaxis=strtod(optarg, &endptr);
      if (!endptr || *endptr!='\0') {
	error( "Xaxis offset argument to -X = %s must be double\n", optarg);
	exit(1);
      }
      break;
    case 'E':    
      /* path to ephemeris files */
      if (!(earthdata=(char *)malloc((strlen(optarg)+strlen(EARTHDATA)+2))) ||
	  !(sundata=  (char *)malloc((strlen(optarg)+strlen(SUNDATA)+2)))) {
	syserror("No memory remaining to store filenames for Ephemeris data\n");
	exit(1);
      }
      /* construct data file names */
      sprintf(earthdata,"%s/%s", optarg, EARTHDATA);
      sprintf(sundata,  "%s/%s", optarg, SUNDATA);
      break;
    case 'D':
      if (!(noisedir=(char *)malloc((strlen(optarg)+2)))){
	syserror("No memory remaining to store input sft dir name\n");
	exit(1);
      }
      sprintf(noisedir, "%s", optarg);
      break; 
      /* input sft directory */
    case 'w':
      /* window data in time domain before FFTing it */
      do_windowing=1;
      break;
    case 'b':
      /* output data in binary format for on-line injection studies */
      binaryoutput=1;
      break;
    case 'm':
      /* DON'T output MAGIC 1234.5 as first sample when using -b flag above */
      nomagic=1;
      break;
    case 'h':
      usage(stdout);
      exit(0);
      break;
    case 'd':
      lalDebugLevel = atoi (optarg);
      break;
    default:
      
      /* unrecognized option */
      errflg++;
      if (c == '?')
	error("Unrecognized option argument -%c\n",optopt);
      else
	error("Option argument -%c requires an argument!\n",optopt);
      usage(stderr);
      exit(1); 
    }
  
  if (!detectorset) {
    error( "You must use the -I option to choose an IFO location\n");
    exit(1);
  }

  /* check that ephemeris files exist and can be read... */
  if (!(fp=fopen(earthdata, "r")) || fclose(fp)) {
    syserror("Unable to read ephemeris file %s\n", earthdata);
    exit(1);
  }
  if (!(fp=fopen(sundata, "r")) || fclose(fp)) {
    syserror("Unable to read ephemeris file %s\n", sundata);
    exit(1);
  }
  
  /* Open input data file */
  if (!(fp=fopen(inDataFilename,"r"))) {
    syserror("Unable to find the inDataFilename file %s\n",inDataFilename);
    return 1;
  }
  
  if (parseR8(fp, "SFT time baseline Tsft", &Tsft))
    return 1;

  if (parseI4(fp, "# of SFTs nTsft", &nTsft))
    return 1;

  if (parseR8(fp, "minimum frequency fmin", &global_fmin))
    return 1;

  if (parseR8(fp, "bandwidth B", &B))
    return 1;

  if (parseR8(fp, "noise variance sigma", &sigma))
    return 1;
  
  if (parseR4(fp, "Plus polarization amplitude aPlus", &genTayParams.aPlus))
    return 1;

  if (parseR4(fp, "Cross polarization amplitude aCross", &genTayParams.aCross))
    return 1;

  if (parseR4(fp, "polarization angle psi", &genTayParams.psi))
    return 1;

  if (parseR8(fp, "initial phase phi", &genTayParams.phi0))
    return 1;

  if (parseR8(fp, "frequency f0", &genTayParams.f0))
    return 1;

  if (parseR8(fp, "declination [radians] delta (lattitude) ", &genTayParams.position.latitude))
    return 1;

  if (parseR8(fp, "right ascension [radians] alpha (longitude)", &genTayParams.position.longitude))
    return 1;

  if (parseI4(fp, "max spin-down order msp", &msp))
    return 1;

  /* if there are spin down parameters read them in */
  if (msp > 0){
    LALDCreateVector(status, &(genTayParams.f), msp); 
    genTayParams.f->length=msp;
    for (i=0;i<msp;i++){
      REAL8 *fi=&(genTayParams.f->data[i]);
      char name[128];
      sprintf(name,"spin down parameter f_%d", i+1);
      if (parseR8(fp, name, fi))
	return 1;
      /*Redefine spin-down parameters to make them consistent with
	what Teviet's routine wants */
      *fi /= genTayParams.f0;
    }
  }
  
  /* timestamps file name */
  r=fscanf(fp, "%s %[^\n]",timestampsname,dmp);
  if (r==2 && GPStime != -1.0) {
    if (lalDebugLevel) error( "Since -G option given, ignoring timestamps file %s\n",
	   timestampsname);
  }
  else  if (r!=2 && GPStime==-1.0) {
    error("Unable to assign timestamps file name from %s\n",inDataFilename);
    return 1;   
  }
  
  fclose(fp);
  
  /* update global variables global_fmin, B */
  imin=floor(global_fmin*Tsft);
  global_fmin=imin/Tsft;
  /*idealized heterodyning in a band of amplitude B*/
  nsamples=2*ceil(B*Tsft);
  B=(nsamples/2)/Tsft;
  
    
  /*and return */
  return errflg;
}





/* File format

The total number of bytes in a file is
8+4+4*N
where N is  the number of timestamps. The format is:
1.0 (this entry is to be sure that the endian order is correct)
   FORMAT: REAL 8
Number N of timestamps
   FORMAT: INT4
Timestamps:
   FORMAT: INT4 SECONDS

*/


