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
makefakedata [-d debuglevel] [-o outfile] [-f f0] [-p alpha delta]
\end{verbatim}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}

This program uses Teviet's LAL CW signal routines in order to produce SFT
files in the GEO binary format.


The signal parameters are specified in an input file. The default name for the input file is In.data, however a file with another name can also be used, as long as the name is specified with the -i command line argument. The input file must contain the following information:
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
The information in parenthesis above shows what appears in the In.data input file as a comment to help you remember what the different entries are.

A number of SFTs will be created. The names will be 
NAME.00001
NAME.00002
.....
and so on.
The default name for the SFT files is TEST\_SFT however a different name can be specified using the -n command line argument.

How many SFTs will be created is specified in the input file mentioned above (In.data is the default).

The time of the first datum of each SFT has to be specified. This is done with a time-stamps file. The time-stamps file must have at least as many rows as the number of SFTs that one wants to create and it has two columns: one column for the gps seconds and one for the gps nanoseconds of each time-stamp. The name of the time-stamps file is specified in the input file (In.data default name).

If one uses the command line argument -t to specify a filename, say 
TIMEFILES, then a set of files:
TIMEFILES.00001
TIMEFILES.00002
.....
and so on
is created. These contain the time series data that each SFT is computed from. Note that Teviet's routines allow idealized heterodyning that is used inthis code to reduce the amount of produced data.

If sigma in the input file is negative, then the noise data is read in from 
SFT files that are specified in the code.

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

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALConfig.h>
#include <glob.h>
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
#include "lal/LALVersion.h"

/* Locations of the earth and sun ephemeris data */
#define EARTHDATA "./earth03.dat"
#define SUNDATA "./sun03.dat"

int mycalls=0;
int myclears=0;
void *ptr;

/* For debugging memory allocation */
#if(0)
#define LALMalloc(how)  ptr=LALMalloc(how);mycalls++,fprintf(stderr,"Allocating %d bytes with %d calls tp %p\n",(INT4)(how),mycalls,ptr)
#define LALFree(quantity) LALFree(ptr=quantity); myclears++; fprintf(stderr,"Clearing with %d calls to LALFree at %p\n",myclears,ptr)
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
#define MAXFILES 40000         // Maximum # of files in a directory 
#define MAXFILENAMELENGTH 256   // Maximum # of characters of a SFT filename

/*Default input data file name*/
char *inDataFilename="In.data";
char timestampsname[128];
char *basefilename="SFT";
char *timebasefilename="/scratch/papa/data4/TIME";
char *noisedir="/scratch/papa/AllSkyH1/";
char filelist[MAXFILES][MAXFILENAMELENGTH];

/* timebaseline of SFT in sec, band SFT in Hz */
REAL4 Tsft,B,sigma;
/* smallest frequency in the band B */
REAL8 fmin;
/* How many SFTS we'll produce*/
INT4 nTsft;

/* Do you want noise ?*/
INT2 inoise;

/*SCALE factor*/
REAL4 scale=1.E19;

/* Our detector*/
LALDetector Detector;
int siteinfo=0;

/* Stores earth/sun ephemeris data for barycentering */
EphemerisData *edat=NULL;
EarthState earth;
EmissionTime emit;

/* Stores detector location and other barycentering data 
 for the first and the follow-up demodulation */
BarycenterInput baryinput;

/* Store time stamps from SFT data */
LIGOTimeGPS *timestamps=NULL;
LIGOTimeGPS *SSBtimestamps=NULL;

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

INT4 lalDebugLevel=3;

/* Prototypes for the functions defined in this file */
int read_file(LALStatus *, int argc, char *argv[]);
int set_default(LALStatus *);
int read_timestamps(LALStatus *);
int write_modulated_amplitudes_file(LALStatus *);
int prepare_baryinput(LALStatus *);
int prepare_cwDetector(LALStatus *);
int prepare_timeSeries(LALStatus *);
int compute_SSBtimes(LALStatus *);
int prepare_fvec(LALStatus *);
int make_noise(LALStatus *);
int make_filelist(LALStatus *);
int read_noise(LALStatus *, int iSFT);
int add(LALStatus *);
int write_SFTS(LALStatus *, int iSFT);
int write_timefile(LALStatus *, int iSFT);
int cleanup(LALStatus *);
int freemem(LALStatus *);
int window_data(LALStatus *);
int correct_phase(LALStatus *, int iSFT);

int main(int argc,char *argv[]) {

  static LALStatus status;
  int iSFT;
  

printf("LAL Version: %s\nMajor Version: %d\nMinor Version: %d\nConfig Args: %s\nBuild Date: %s\n",
lalVersion,lalVersionMajor,lalVersionMinor,lalConfigureArgs,lalConfigureDate);

   
  /* Define Detector (to GEO) 
   set genTayParams.f = NULL; */
  if (set_default(&status))
    return 1;

  /* read input variables from input data file */
  if (read_file(&status,argc,argv))
    return 1;

    if (sigma < 0.0){
      /* make file list for noise data files */
      if (make_filelist(&status))
	return 1;
    }

  /* Read in timestamps and place them in timestamps vec*/
  if (read_timestamps(&status))
    return 1;

  /* Fill-in genTayParams*/
  if (prepare_baryinput(&status))
    return 1;

   /* compute and write SSB times in SSBtimestamps */
  if (compute_SSBtimes(&status))
    return 1;


  /* Fill-in cwDetector*/
  if (prepare_cwDetector(&status))
    return 1;

  /*Allocate space for time series that will hold SFT time chunks
   and fill-in all fields that are independent of the particular SFT*/
  if (prepare_timeSeries(&status))
    return 1;

  /*Allocate space for SFT vector and fill-in appropriate fields*/
  if (prepare_fvec(&status))
    return 1;

  /* If option active, write amplitude-modulation info file*/
#if(0)
  if (write_modulated_amplitudes_file(&status))
    return 1;
#endif

 

  memset(&cgwOutput, 0, sizeof(CoherentGW));
  
   /* complete filling-in the "input fields" of genTayparams */
  genTayParams.epoch.gpsSeconds = SSBtimestamps[0].gpsSeconds; 
  genTayParams.epoch.gpsNanoSeconds = SSBtimestamps[0].gpsNanoSeconds; 
  
  genTayParams.length=(LTT+Tsft+timestamps[nTsft-1].gpsSeconds-timestamps[0].gpsSeconds)/genTayParams.deltaT;

  /*Generate signal at the source for entire observation time */
  SUB( LALGenerateTaylorCW(&status, &cgwOutput, &genTayParams), &status);

  /******** check that the output is OK ********************************/
  if ( genTayParams.dfdt > 2.0 ) {
    /* LALSnprintf() can't seem to print floating-point formats.
       LALSnprintf( message, MSGLEN,
       "Waveform sampling interval is too large:\n"
       "\tmaximum df*dt = %f", params.dfdt );
    */
    printf( "Waveform sampling interval is too large:\n"
	    "\tmaximum df*dt = %f", genTayParams.dfdt );
    return 1;
  }
/***********************************************************************/


  for (iSFT=0;iSFT<nTsft;iSFT++){

    /* This matters a lot */
    timeSeries->epoch=timestamps[iSFT];
    
    SUB( LALSimulateCoherentGW(&status, timeSeries, &cgwOutput, &cwDetector), &status);

    /*if you want noise*/
    if (sigma > 0.0){
      /*produce noise and add it to timeSeries*/
      if (make_noise(&status))
	return 1;
    }

#if(0)
    /*window data*/
    if (window_data(&status))
      return 1;
#endif

    /* Perform FFTW-LAL Fast Fourier Transform */
#if(1)    
    SUB( LALForwardRealFFT(&status, fvec, timeSeries->data, pfwd), &status);
#endif


/*if you want noise*/
    if (sigma < 0.0 ){
      /*read data and add it to fvec*/
      if (read_noise(&status, iSFT))
	return 1;
    }

#if(0)
    if (timebasefilename != NULL){
      if (write_timefile(&status, iSFT))
	return 1;
    }
#endif

    /*Write the SFT file in the file that is specified by the path name*/
#if(1)    
    if (write_SFTS(&status, iSFT))
      return 1;
#endif

  }
  
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
  LALFree(SSBtimestamps);
 

  /* Clean up FFTs of signal and of noise - each a complex8vector */
  LALCDestroyVector(status, &fvec);
  LALDestroyRealFFTPlan(status, &pfwd);
  

  if ( genTayParams.f )
    LALDDestroyVector( status, &(genTayParams.f) );


  if (sigma > 0.0 || sigma < 0)
    LALCDestroyVector(status, &fvecn);
  

  return 0;
}


/*sets default values */
int set_default(LALStatus* status) {

  genTayParams.f=NULL;
  genTayParams.position.system= COORDINATESYSTEM_EQUATORIAL;
  genTayParams.deltaT=60;

  return 0;
}


/*sets default values */
int compute_SSBtimes(LALStatus* status) {
  
  int i;
  REAL8 floatTime;
  REAL4 Ts,Tns;
  LIGOTimeGPS ssb;
  
  for (i=0; i < nTsft; ++i){
    
    Ts=timestamps[i].gpsSeconds;
    Tns=timestamps[i].gpsNanoSeconds;

    LALBarycenterEarth(status, &earth, &timestamps[i], edat);
    LALBarycenter(status, &emit, &baryinput, &earth);
    floatTime= emit.deltaT + Ts + Tns*1.E-9;
    LALFloatToGPS(status,&ssb, &floatTime);
    SSBtimestamps[i]=ssb;
    
  }
  return 0;
  
}



/*applies correct heterodyning from one SFT to the next*/
int correct_phase(LALStatus* status, int iSFT) {

  int i;
  REAL8 cosx,sinx,x;
  COMPLEX8 fvec1;

  x=2.0*LAL_PI*fmin*Tsft*iSFT;
  cosx=cos(x);
  sinx=sin(x);
  for (i = 0; i < fvec->length; ++i){
    fvec1=fvec->data[i];
    fvec->data[i].re=fvec1.re*cosx+fvec1.im*sinx;
    fvec->data[i].im=fvec1.im*cosx-fvec1.re*sinx;
  }
  
  return 0;
}


/* windows the data*/
int window_data(LALStatus* status){

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
        printf("Not enough data points to window !\n");
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






/* reads ephemeris data files and fills-in baryinput fields*/
int prepare_baryinput(LALStatus* status){

  /* Quantities computed for barycentering */
  edat=(EphemerisData *)LALMalloc(sizeof(EphemerisData));
  (*edat).ephiles.earthEphemeris = "./earth03.dat";
  (*edat).ephiles.sunEphemeris = "./sun03.dat";
   (*edat).leap=13; 

  /* Read in ephemerides */  
  LALInitBarycenter(status, edat);
  /*Getting detector coords from DetectorSite module of tools package */   
  
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
  
  /*  cwDetector.heterodyneEpoch=(LIGOTimeGPS *)LALMalloc(sizeof(LIGOTimeGPS)); */
  /* SSBtimestamps or not, without heterodyning it does not seem to make a difference*/
  cwDetector.heterodyneEpoch.gpsSeconds=timestamps[0].gpsSeconds;
  cwDetector.heterodyneEpoch.gpsNanoSeconds=timestamps[0].gpsNanoSeconds;

  return 0;
}



/*Allocates space for timeseries */
int prepare_timeSeries(LALStatus* status) {

  timeSeries = (REAL4TimeSeries *)LALMalloc(sizeof(REAL4TimeSeries));
  timeSeries->data = (REAL4Vector *)LALMalloc(sizeof(REAL4Vector));
  
  timeSeries->data->length = B*Tsft*2;
  timeSeries->data->data = (REAL4 *)LALMalloc(timeSeries->data->length*sizeof(REAL4));
  timeSeries->deltaT=Tsft/timeSeries->data->length;

  timeSeries->f0=fmin;

  return 0;
}


/*Allocates space for timeseries */
int prepare_fvec(LALStatus* status) {

  INT4 len,len2;

  len=timeSeries->data->length;
  len2=(len/2)+1;
  
  /* Create vector to hold signal frequency series */
  LALCCreateVector(status, &fvec, (UINT4)len2);
  
  if (sigma > 0.0 || sigma < 0.0)
    LALCCreateVector(status, &fvecn, (UINT4)len2);
  
  /* Compute measured plan for FFTW */
  LALCreateForwardRealFFTPlan(status, &pfwd, (UINT4)len, 0);
  
  return 0;
}



int add(LALStatus* status) {

  int i;

  for (i = 0; i < fvec->length; ++i)
  {
    fvec->data[i].re=fvec->data[i].re+fvecn->data[i].re;
    fvec->data[i].im=fvec->data[i].im+fvecn->data[i].im;
  }


  return 0;
}




int make_noise(LALStatus* status) {


  REAL4Vector    *v1=NULL;
  RandomParams   *randpar=NULL;
  INT4            i,numPoints,seed,errorcode;
  FILE *devrandom;

  numPoints=timeSeries->data->length;


  LALCreateVector (status, &v1, numPoints);

  /*
   * Modified so as to not create random number parameters with seed drawn from clock.
   * Seconds don't change fast enough and sft's look alike.
   * We open /dev/urandom and read a 4 byte integer from it and use that 
as our seed.  Note: /dev/random is slow after the first, few accesses.
   */

  devrandom=fopen("/dev/urandom","r");
  errorcode=fread((void*)&seed,sizeof(INT4),1,devrandom);
  if (errorcode!=1)
    {
    printf("Error reading /dev/random file!\n");
    return 1;
    }
  fclose(devrandom);
  
  LALCreateRandomParams (status, &randpar, seed);

  LALNormalDeviates(status, v1, randpar);

  for (i = 0; i < numPoints; ++i)
  {
    timeSeries->data->data[i]=timeSeries->data->data[i]+sigma*v1->data[i];
  }

/* destroy randpar*/
  LALDestroyRandomParams (status, &randpar);
/*   destroy v1 */
  LALDestroyVector (status, &v1);

  return 0;

}




/*reads timestamps file and fills-in timestamps vector*/
int read_timestamps(LALStatus* status) {
  
  FILE *fp;
  int i,r;
 
  
  /*   %strcpy(filename,inDataFilename); */
  fp=fopen(timestampsname,"r");
  if (fp==NULL) {
    fprintf(stderr,"Unable to find file %s\n",timestampsname);
    return 1;
  }
  timestamps=(LIGOTimeGPS *)LALMalloc(nTsft*sizeof(LIGOTimeGPS)); 
  SSBtimestamps=(LIGOTimeGPS *)LALMalloc(nTsft*sizeof(LIGOTimeGPS)); 

  for (i=0;i<nTsft;i++){
    r=fscanf(fp,"%d  %d\n", &timestamps[i].gpsSeconds, &timestamps[i].gpsNanoSeconds);
    if ( r !=2 ) {
      fprintf(stderr,"Unable to read datum # %d\n",i);
      fprintf(stderr,"from file %s\n",timestampsname);
      return 1; 
    } 
  } 
  
  fclose(fp);
  return 0;
  
}




int write_modulated_amplitudes_file(LALStatus* status){
  FILE *fp;
  LALDetAMResponse  amresp;
  LALSource         source;
  LALDetAndSource   detectorandsource;
  LIGOTimeGPS       gps;
  char *filename="AmplMod.dat";
  int i;


  fp=fopen(filename,"w");
  if (fp==NULL) {
    fprintf(stderr,"Unable to find file %s\n",filename);
    return 1;
  }

  detectorandsource.pDetector=&Detector;
  source.equatorialCoords=genTayParams.position;
  source.orientation=genTayParams.psi;  
  detectorandsource.pSource=&source;

  for (i=0;i<=nTsft;i++){
  
    gps.gpsSeconds=timestamps[i].gpsSeconds;
    gps.gpsNanoSeconds=timestamps[i].gpsNanoSeconds;

    LALComputeDetAMResponse(status, &amresp, &detectorandsource, &gps);
    fprintf(fp,"%f  %f\n",amresp.plus,amresp.cross);
    
  }
  
  return 0;
}




int make_filelist(LALStatus* status) {

  INT4 fileno=0;
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
	  fprintf(stderr,"Too many files in directory! Exiting... \n");
	  return 1;
	}
    }
  globfree(&globbuf);


  return 0;

}





int read_noise(LALStatus* status, int iSFT) {

  FILE *fp;
  REAL4 norm;
  size_t errorcode;
  INT4 i;

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
  fp=fopen(filelist[iSFT],"r");
  /* read in the header from the file */
  errorcode=fread((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1) 
    {
      fprintf(stderr,"No header in data file %s\n",filelist[iSFT]);
      return 1;
    }

  /* check that data is correct endian order */
  if (header.endian!=1.0)
    {
      fprintf(stderr,"First object in file %s is not (double)1.0!\n",filelist[iSFT]);
      fprintf(stderr,"It could be a file format error (big/little\n");
      fprintf(stderr,"endian) or the file might be corrupted\n\n");
      return 2;
    }


  /* check frequency range */
  if ((fmin*Tsft < header.firstfreqindex) ||
      ((fmin*Tsft+fvec->length-1) > header.firstfreqindex + header.nsamples)){
    fprintf(stderr,"Frequency band of noise data out of range !\n");
    return 1;
  }

  /* seek to position */
  if (0 != fseek(fp, (fmin*Tsft - header.firstfreqindex) * 4.0 * 2, SEEK_CUR)){
    fprintf(stderr,"file too short (could'n fssek to position\n");
    return 1;
  }


  /* calculate size of the data in bytes */
  /* datasize = (f2ind - f1ind) * head.dsize * 2;*/

  /* read the data */
  if (1 != fread(fvecn->data,(fvecn->length-1)*2*4.0,1,fp)) {
      fprintf(stderr,"Could not read the data \n");
    return 1;
  }

  fclose(fp);

  norm=((REAL4)(fvec->length-1)*1.0/((REAL4)header.nsamples));

  for (i = 0; i < fvec->length; ++i)
  {
    fvec->data[i].re=fvec->data[i].re+scale*fvecn->data[i].re*norm;
    fvec->data[i].im=fvec->data[i].im+scale*fvecn->data[i].im*norm;
  }
  
  return 0;
}






int write_SFTS(LALStatus* status, int iSFT){


  FILE *fp;
  REAL4 rpw,ipw;
  REAL8 fr;
  char filename[256], filenumber[16];
  int i,errorcode;

  struct headertag {
    REAL8 endian;
    INT4  gps_sec;
    INT4  gps_nsec;
    REAL8 tbase;
    INT4  firstfreqindex;
    INT4  nsamples;
  } header;
  
  /* Open SFT data file */
  strcpy(filename,basefilename);
  sprintf(filenumber,".%05d",iSFT); /*ACHTUNG: used to be iSFT+1 - now starts from .00000*/
  strcat(filename,filenumber);
  fp=fopen(filename,"w");
  if (fp==NULL) {
    fprintf(stderr,"Unable to find file %s\n",filename);
    return 1;
  }

  header.endian=1.0;
  header.gps_sec=timestamps[iSFT].gpsSeconds;
  header.gps_nsec=timestamps[iSFT].gpsNanoSeconds;
  header.tbase=Tsft;
  header.firstfreqindex=(INT4)(fmin*Tsft+0.5);
  header.nsamples=fvec->length-1;
  
  /* write header */
  errorcode=fwrite((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1){
    printf("Error in writing header into file!\n");
    return 1;
  }

  for (i=0;i<fvec->length-1;i++){

    rpw=fvec->data[i].re;
    ipw=fvec->data[i].im;

    errorcode=fwrite((void*)&rpw, sizeof(REAL4),1,fp);  
    if (errorcode!=1){
      printf("Error in writing data into SFT file!\n");
      return 1;
    }
    errorcode=fwrite((void*)&ipw, sizeof(REAL4),1,fp);  
    if (errorcode!=1){
      printf("Error in writing data into SFT file!\n");
      return 1;
    }
        
  }
  
  fclose(fp);
  return 0;  
  
}




int write_timefile(LALStatus* status, int iSFT){

  FILE *fp;
  REAL4 pw;
  char filename[256], filenumber[16];
  int i,errorcode;

  struct headertag {
    REAL8 endian;
    INT4  gps_sec;
    INT4  gps_nsec;
    REAL8 tbase;
    INT4  firstfreqindex;
    INT4  nsamples;
  } header;
  
  /* Open SFT data file */
  strcpy(filename,timebasefilename);
  sprintf(filenumber,".%05d",iSFT+1);
  strcat(filename,filenumber);
  fp=fopen(filename,"w");
  if (fp==NULL) {
    fprintf(stderr,"Unable to open file %s\n",filename);
    return 1;
  }
  header.endian=1.0;
  header.gps_sec=timestamps[iSFT].gpsSeconds;
  header.gps_nsec=timestamps[iSFT].gpsNanoSeconds;
  header.tbase=Tsft;
  header.firstfreqindex=(INT4)(fmin*Tsft);
  header.nsamples=timeSeries->data->length;
  
  /* write header */
#if(0)
  fprintf(fp, "%e\n",header.endian);
#endif
  
  for (i=0;i<header.nsamples;i++){
    
    pw=timeSeries->data->data[i];
    fprintf(fp,"%f\n",pw);  
        
  }
  
  fclose(fp);
  return 0;  
  
}





/* options are:
   -i  name of input data file
   -n  name of SFT files
   -t  name of TDD (time-domain-data) files
   -h  print help
*/
int read_file(LALStatus* status, int argc,char *argv[]) {
  
  char filename[256], dmp[128];
  int c, errflg = 0;
  int r,i,msp;
  UINT4 imin, nsamples;  
  FILE *fp;
  extern char *optarg;
  
  /* scan through the list of arguments on the command line 
     and get the input data filename*/
  
  while (!errflg && ((c = getopt(argc, argv,"i:n:hb:t:o:I:"))!=-1))
    switch (c) {
      
    case 'i':
      /* Name of input data file */
      inDataFilename=optarg;
      break;
    case 'n':
      /* Name of SFT data file */
      basefilename=optarg;
      break;
    case 't':
      /* Name of TDD file */
      timebasefilename=optarg;
      break;
    case 'I':
      /* detector site*/
      siteinfo=atol(optarg);
      if (siteinfo > 3 ) {
	      fprintf(stderr,"Chosen detector: GEO600\n");
      }
      break;
    case 'h':
      /* print usage/help message */
      fprintf(stderr,"Arguments are:\n");
      fprintf(stderr,"-i\tCharacter String\t( Name of input parameter file)\n");
      fprintf(stderr,"-n\tCharacter String\t( Basefilename of output data files)\n");
      fprintf(stderr,"-I\tInteger number\t\t( IFO: default: GEO, 1 LLO, 2 LHO, 3 Roman Bar)\n");
      exit(0);
      break;
    default:
      /* unrecognized option */
      errflg++;
      fprintf(stderr,"Unrecognized option argument %c\n",c);
      exit(1);
      break;
    }

  Detector=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if (siteinfo == 1)
    Detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
  if (siteinfo == 2)
    Detector=lalCachedDetectors[LALDetectorIndexLHODIFF];
  if (siteinfo == 3)
    {
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
    }


  /* Open input data file */
  strcpy(filename,inDataFilename);
  fp=fopen(filename,"r");
  if (fp==NULL) {
    fprintf(stderr,"Unable to find file %s\n",filename);
    return 1;
  }

  /*Which detector*/
/*   r=fscanf(fp, "s\n",&detectorname); */
/*   if ( r !=1 ) { */
/*     fprintf(stderr,"Unable to assign f0 from %s\n",filename); */
/*     return 1;    */
/*   } */
  
  
  /*Tsft REAL8*/
  r=fscanf(fp, "%f %s\n",&Tsft,dmp);
  if ( r !=2 ) {
    fprintf(stderr,"Unable to assign Tsft from file: %s\n",filename);
    return 1;   
  }
  r=fscanf(fp, "%d %s\n",&nTsft,dmp);
  if ( r !=2 ) {
    fprintf(stderr,"Unable to assign how many SFTs from %s\n",filename);
    return 1;   
  }
  /*fmin REAL8*/
  r=fscanf(fp, "%lf %s\n",&fmin,dmp);
  if ( r !=2 ) {
    fprintf(stderr,"Unable to assign fmin from %s\n",filename);
    return 1;   
  }
  /*Band REAL4*/
  r=fscanf(fp, "%f %s\n",&B,dmp);
  if ( r !=2 ) {
    fprintf(stderr,"Unable to assign band width from %s\n",filename);
    return 1;   
  }
  /*sigma REAL4*/
  r=fscanf(fp, "%f %s\n",&sigma,dmp);
  if ( r !=2 ) {
    fprintf(stderr,"Unable to assign sigma from %s\n",filename);
    return 1;   
  }
  /*Ap REAL4*/
  r=fscanf(fp, "%f %s\n",&genTayParams.aPlus,dmp);
  if ( r !=2 ) {
	  char command[1024];
    fprintf(stderr,"Unable to assign A plus from %s\n",filename);
   sprintf(command, "cat %s 1>&2\n", filename);
    system(command);

    return 1;   
  }
  /*Ac REAL4*/
  r=fscanf(fp, "%f %s\n",&genTayParams.aCross,dmp);
  if ( r !=2 ) {
	  char command[1024];
    fprintf(stderr,"Read %d items: Unable to assign A cross from %s.\n", r, filename);
    if (r==1) fprintf(stderr,"Have A cross=%f\n", genTayParams.aCross);
    if (r==0) fprintf(stderr,"Last item read was %s\n", dmp);
    sprintf(command, "cat %s 1>&2\n", filename);
    system(command);
    return 1;   
  }
  /*psi REAL4*/
  r=fscanf(fp, "%f %s\n",&genTayParams.psi,dmp);
  if ( r !=2 ) {
	          char command[1024];
    fprintf(stderr,"Unable to assign psi from %s\n",filename);
    sprintf(command, "cat %s 1>&2\n", filename);
    system(command);
    
    return 1;   
  }
  /*phi0 REAL8*/
  r=fscanf(fp, "%lf %s\n",&genTayParams.phi0,dmp);
  if ( r !=2 ) {
	          char command[1024];
    fprintf(stderr,"Unable to assign phi0 from %s\n",filename);
   sprintf(command, "cat %s 1>&2\n", filename);
    system(command);

    return 1;   
  }
  /*f0 REAL8*/
  r=fscanf(fp, "%lf %s\n",&genTayParams.f0,dmp);
  if ( r !=2 ) {
    fprintf(stderr,"Unable to assign f0 from %s\n",filename);
    return 1;   
  }
  /*alpha - from input file in degrees. Converted into radians.*/
  r=fscanf(fp, "%lf %s\n",&genTayParams.position.latitude,dmp);
  if ( r !=2 ) {
    fprintf(stderr,"Unable to assign source latitude from %s\n",filename);
    return 1;   
  }
/*  genTayParams.position.latitude=genTayParams.position.latitude*LAL_TWOPI/360.0; */

  /*delta - from input file in degrees. Converted into radians.*/
  r=fscanf(fp, "%lf %s\n",&genTayParams.position.longitude,dmp);
  if ( r !=2 ) {
    fprintf(stderr,"Unable to assign source longitude from %s\n",filename);
    return 1;   
  }
/*  genTayParams.position.longitude=genTayParams.position.longitude*LAL_TWOPI/360.0; */

  /* max spin-down parameter order */
  r=fscanf(fp, "%d %s\n",&msp,dmp);
  if ( r !=2 ) {
    fprintf(stderr,"Unable to assign max spin-down order from %s\n",filename);
    return 1;   
  }
  /* if there are spin down parameters read them in */
  if (msp > 0){
    LALDCreateVector(status, &(genTayParams.f), msp); 
    genTayParams.f->length=msp;
    for (i=0;i<msp;i++){
      r=fscanf(fp, "%le %s\n",&genTayParams.f->data[i],dmp);
      if ( r !=2 ) {
	fprintf(stderr,"Unable to assign spin-down value %ih\n",i);
	fprintf(stderr,"from file %s\n",filename);
	return 1;   
      }/*endif it was read in OK*/
      /*Redefine spin-down parameters to make them consistent with 
	what Teviet's routine wants */
      genTayParams.f->data[i]=genTayParams.f->data[i]/genTayParams.f0;
    }/*endfor over different spindown orders*/
  }/*endif there were spindown values at all*/
  
  /* timestamps file name */
  r=fscanf(fp, "%s %s\n",timestampsname,dmp);
  if ( r !=2 ) {
    fprintf(stderr,"Unable to assign timestamps file name from %s\n",filename);
    return 1;   
  }
  
  /* update global variables fmin, B */
  imin=floor(fmin*Tsft);
  fmin=imin/Tsft;
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


