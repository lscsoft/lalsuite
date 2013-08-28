/*
*  Copyright (C) 2007 Chris Messenger, Reinhard Prix
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

/*-----------------------------------------------------------------------
 *
 * File Name: makebinaryfakedata_v1.c
 *
 * Authors: Papa, M.A.,
 *
 * History:   Created by Papa July 2002
 *            Modified by X. Siemens to fix random seed bug
 *            and modified by C. Messenger to inject binary signal
 *
 *-----------------------------------------------------------------------
 */

/**
 * \file
 * \ingroup pulsarApps
 * \author M.A. Papa
 * \brief
 * Produces fake SFT data.
 *
 * \heading{Usage}
 * \code
 * makefakedatabinary [-d debuglevel] [-o outfile] [-f f0] [-p alpha delta]
 * \endcode
 *
 * \heading{Description}
 *
 * This program uses Teviet's LAL CW signal routines in order to produce SFT
 * files in the GEO binary format.
 *
 * The signal parameters are specified in an input file. The default name for the input file is In.data, however a file with another name can also be used, as long as the name is specified with the -i command line argument. The input file must contain the following information:
 * <table><tr><td>
 * time-baseline of the SFT</td><td>(Tsft_in_sec)</td></tr>
 * <tr><td>number of SFTs to be produced</td><td>(nTsft)</td></tr>
 * </table>
 *
 * frequency of first bin of SFTs       (first_SFT_frequency_in_Hz)
 * band of produced SFT data            (SFT_freq_band_in_Hz)
 * standard deviation of noise
 * (for real and imag) SFT          (std_of_noise.When=0_only_signal_present)
 * amplitude of plus polarization       (Aplus)
 * amplitude of cross polarization      (Across)
 * polarization angle                   (psi)
 * initial phase                        (phi0)
 * intrinsic emission frequency
 * at the beginning of observation   (f0)
 * position of source (eq. coordinates)  (latitude_in_degrees)
 * (longitude_in_degrees)
 *
 * Orbital Semi-major axis(sec)
 * Orbital Period(sec)
 * Time of observed periapse passage of binary in SSB frame(GPS seconds)
 * Time of observed periapse passage of binary in SSB frame(GPS nanoseconds)
 * Argument of periapse(rad)
 * Orbital eccentricity
 * First spin down parameter(Hz/sec)
 *
 * maximum spin-down order               (max_spin-down_param_order)
 * name of time-stamps file              (name_of_time-stamps_file)
 * The information in parenthesis above shows what appears in the In.data input file as a comment to help you remember what the different entries are.
 *
 * A number of SFTs will be created. The names will be
 * NAME.00001
 * NAME.00002
 * .....
 * and so on.
 * The default name for the SFT files is TEST_SFT however a different name can be specified using the -n command line argument.
 *
 * How many SFTs will be created is specified in the input file mentioned above (In.data is the default).
 *
 * The time of the first datum of each SFT has to be specified. This is done with a time-stamps file. The time-stamps file must have at least as many rows as the number of SFTs that one wants to create and it has two columns: one column for the gps seconds and one for the gps nanoseconds of each time-stamp. The name of the time-stamps file is specified in the input file (In.data default name).
 *
 * If one uses the command line argument -t to specify a filename, say
 * TIMEFILES, then a set of files:
 * TIMEFILES.00001
 * TIMEFILES.00002
 * .....
 * and so on
 * is created. These contain the time series data that each SFT is computed from. Note that Teviet's routines allow idealized heterodyning that is used in this code to reduce the amount of produced data.
 *
 * If sigma in the input file is negative, then the noise data is read in from
 * SFT files that are specified in the code.
 *
 * \heading{Uses}
 * \code
 * LALGenerateTaylorCW()
 * LALPulsarSimulateCoherentGW()
 * LALMalloc()
 * LALSCreateVector()
 * LALSDestroyVector()
 * LALFree()
 * LALCheckMemoryLeaks()
 * \endcode
 *
 */

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALStdlib.h>

/**\name Error Codes */ /*@{*/
#define MAKEFAKEDATABINARYC_ENORM 0
#define MAKEFAKEDATABINARYC_ESUB  1
#define MAKEFAKEDATABINARYC_EARG  2
#define MAKEFAKEDATABINARYC_EBAD  3
#define MAKEFAKEDATABINARYC_EFILE 4

#define MAKEFAKEDATABINARYC_MSGENORM "Normal exit"
#define MAKEFAKEDATABINARYC_MSGESUB  "Subroutine failed"
#define MAKEFAKEDATABINARYC_MSGEARG  "Error parsing arguments"
#define MAKEFAKEDATABINARYC_MSGEBAD  "Bad argument values"
#define MAKEFAKEDATABINARYC_MSGEFILE "Could not create output file"
/*@}*/


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
#include <lal/PulsarSimulateCoherentGW.h>
#include <lal/GenerateTaylorCW.h>
#include <lal/LALDatatypes.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LALVersion.h>
#include <lal/GenerateSpinOrbitCW.h>

#include "getopt.h"

/* Locations of the earth and sun ephemeris data */
CHAR earthdata[]= "./earth00-04.dat";
CHAR sundata[] = "./sun00-04.dat";

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
    XLALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, "$Id$", statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    XLALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              "$Id$", (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( MAKEFAKEDATABINARYC_ESUB, MAKEFAKEDATABINARYC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return MAKEFAKEDATABINARYC_ESUB;                                  \
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
#define MAXFILES 40000         /* Maximum # of files in a directory */
#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename */

/*Default input data file name*/
const char *inDataFilename="In_bin.data";
char timestampsname[128];
const char *basefilename="SFT";
const char *timebasefilename="./scratch/TIME_e";
char *noisedir;        /* ="/gwave1/cm/GWsearch/data/GEOS1/"; */
char *programname=NULL;
char filelist[MAXFILES][MAXFILENAMELENGTH];

/* timebaseline of SFT in sec, band SFT in Hz */
REAL4 Tsft,B,sigma;
/* smallest frequency in the band B */
REAL8 f_min;
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
SpinOrbitCWParamStruc genTayParams;
PulsarCoherentGW cgwOutput;
PulsarDetectorResponse cwDetector;
COMPLEX8FrequencySeries *transferFunction;

/*This will hold the SFT*/
COMPLEX8Vector *fvec = NULL;
COMPLEX8Vector *fvecn = NULL;

/*FFT plan*/
RealFFTPlan *pfwd = NULL;

/*Calculate power in signal*/
REAL8 power;

/* some temporary orbital parameters */
REAL8 SemiMajorAxis;
LIGOTimeGPS TperiapseSSB;
REAL8 OrbitalEccentricity;
REAL8 ArgPeriapse;
REAL8 OrbitalPeriod;


/* Prototypes for the functions defined in this file */
int read_file(LALStatus *, int argc, char *argv[]);
int set_default(void);
int read_timestamps(void);
int write_modulated_amplitudes_file(LALStatus *);
int prepare_baryinput(LALStatus *);
int prepare_cwDetector(LALStatus *);
int prepare_timeSeries(void);
int compute_SSBtimes(LALStatus *);
int prepare_fvec(LALStatus *);
int make_noise(LALStatus *);
int make_filelist(void);
int read_noise(LALStatus *, int iSFT);
int write_SFTS(int iSFT);
int write_timefile(int iSFT);
int cleanup(LALStatus *);
int freemem(LALStatus *);
int window_data(void);
int correct_phase(void);
int compute_power(void);
int SetupSigGenParams(void);
static void TimeToFloat(REAL8 *, LIGOTimeGPS *);
static void FloatToTime(LIGOTimeGPS *, REAL8 *);

void syserror(const char *fmt, ...);
void error(const char *fmt, ...);

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

int main(int argc,char *argv[]) {

  static LALStatus status;
  int iSFT;

  programname=argv[0];

/* printf("LAL Version: %s\nMajor Version: %d\nMinor Version: %d\nConfig Args: %s\nBuild Date: %s\n", */
/* lalVersion,lalVersionMajor,lalVersionMinor,lalConfigureArgs,lalConfigureDate); */


  /* Define Detector (to GEO)
   set genTayParams.f = NULL; */
  if (set_default())
    return 1;

  /* read input variables from input data file */
  if (read_file(&status,argc,argv))
    return 1;

    if (sigma < 0.0){
      /* make file list for noise data files */
      if (make_filelist())
	return 1;
    }

  /* Read in timestamps and place them in timestamps vec*/
  if (read_timestamps())
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
  if (prepare_timeSeries())
    return 1;

  /*Allocate space for SFT vector and fill-in appropriate fields*/
  if (prepare_fvec(&status))
    return 1;

  /* If option active, write amplitude-modulation info file*/
#if(0)
  if (write_modulated_amplitudes_file(&status))
    return 1;
#endif



  memset(&cgwOutput, 0, sizeof(PulsarCoherentGW));

  /* complete filling-in the "input fields" of genTayparams */
  if (SetupSigGenParams())
    return 1;

  /*Generate signal at the source for entire observation time */
  SUB( LALGenerateSpinOrbitCW(&status, &cgwOutput, &genTayParams), &status);      /* changed for orbital motion */

  /******** check that the output is OK ********************************/
  if ( genTayParams.dfdt > 2.0 ) {
    /* snprintf() can't seem to print floating-point formats.
       snprintf( message, MSGLEN,
       "Waveform sampling interval is too large:\n"
       "\tmaximum df*dt = %f", params.dfdt );
    */
    printf( "Waveform sampling interval is too large:\n"
	    "\tmaximum df*dt = %f", genTayParams.dfdt );
    return 1;
  }
/***********************************************************************/

  power=0.0;

  for (iSFT=0;iSFT<nTsft;iSFT++){

    /* This matters a lot */
    timeSeries->epoch=timestamps[iSFT];

    SUB( LALPulsarSimulateCoherentGW(&status, timeSeries, &cgwOutput, &cwDetector), &status);

    /* lets calculate the power in the signal */
    /* if (compute_power()) */
    /*  return 1; */

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
      if (write_timefile(iSFT))
	return 1;
    }
#endif

    /*Write the SFT file in the file that is specified by the path name*/
#if(1)
    if (write_SFTS(iSFT))
      return 1;
#endif

  }

  /* fprintf(stdout,"%le\n",power*timeSeries->deltaT);*/

  if (cleanup(&status))
    return 1;

  if (freemem(&status))
    return 1;

   LALCheckMemoryLeaks();

/*   INFO(MAKEBINARYFAKEDATAC_MSGENORM); */

/*   return MAKEBINARYFAKEDATAC_ENORM; */

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
  LALCDestroyVector( status, &( transferFunction->data ) );
  LALFree(transferFunction);
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


  if (fvecn)
    LALCDestroyVector(status, &fvecn);


  return 0;
}


/*sets default values */
int set_default(void) {

  genTayParams.f=NULL;
  genTayParams.position.system= COORDINATESYSTEM_EQUATORIAL;
  genTayParams.deltaT=10;                      /* changed to accomodate greater frequency changes */

  return 0;
}


/*sets default values */
int compute_SSBtimes(LALStatus* status) {

  int i;
  REAL8 floatTime;
  REAL8 Ts,Tns;
  LIGOTimeGPS ssb;

  for (i=0; i < nTsft; ++i){

    Ts=timestamps[i].gpsSeconds;
    Tns=timestamps[i].gpsNanoSeconds;

    LALBarycenterEarth(status, &earth, &timestamps[i], edat);
    LALBarycenter(status, &emit, &baryinput, &earth);
    floatTime= emit.deltaT + Ts + Tns*1.E-9;
    XLALGPSSetREAL8(&ssb, floatTime);
    SSBtimestamps[i]=ssb;

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



/* reads ephemeris data files and fills-in baryinput fields*/
int prepare_baryinput(LALStatus* status){

  /* Quantities computed for barycentering */
  edat=(EphemerisData *)LALMalloc(sizeof(EphemerisData));
  (*edat).ephiles.earthEphemeris = earthdata;
  (*edat).ephiles.sunEphemeris = sundata;

  /* Read in ephemerides */
  LALInitBarycenter(status, edat);
  /*Getting detector coords from DetectorSite module of tools package */

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

  memset(&cwDetector, 0, sizeof(PulsarDetectorResponse));
  /* The ephemerides */
  cwDetector.ephemerides = edat;
  /* Specifying the detector site (set above) */
  cwDetector.site = &Detector;
  /* The transfer function.
   * Note, this xfer function has only two points at it extends
   between 0 and 16384 Hz. The routine that will generate the signal as
   output from the detector on Earth will interpolate*/
  transferFunction = (COMPLEX8FrequencySeries *)LALMalloc(sizeof(COMPLEX8FrequencySeries));
  memset(transferFunction, 0, sizeof(COMPLEX8FrequencySeries));
  /* it does not change so just use first timestamp. Does not
   seem to matter whether SSBtimestamps or timestamps are used */
  transferFunction->epoch = timestamps[0];
  transferFunction->f0 = 0.0;
  transferFunction->deltaF = 16384.0;
  transferFunction->data = NULL;
  LALCCreateVector(status, &(transferFunction->data), 2);

  /* unit response function */
  transferFunction->data->data[0].realf_FIXME = 1.0;
  transferFunction->data->data[1].realf_FIXME = 1.0;
  transferFunction->data->data[0].imagf_FIXME = 0.0;
  transferFunction->data->data[1].imagf_FIXME = 0.0;

  cwDetector.transfer = transferFunction;

  /*  cwDetector.heterodyneEpoch=(LIGOTimeGPS *)LALMalloc(sizeof(LIGOTimeGPS)); */
  /* SSBtimestamps or not, without heterodyning it does not seem to make a difference*/
  cwDetector.heterodyneEpoch.gpsSeconds=timestamps[0].gpsSeconds;
  cwDetector.heterodyneEpoch.gpsNanoSeconds=timestamps[0].gpsNanoSeconds;

  return 0;
}



/*Allocates space for timeseries */
int prepare_timeSeries(void) {

  timeSeries = (REAL4TimeSeries *)LALMalloc(sizeof(REAL4TimeSeries));
  timeSeries->data = (REAL4Vector *)LALMalloc(sizeof(REAL4Vector));

  timeSeries->data->length = B*Tsft*2;
  timeSeries->data->data = (REAL4 *)LALMalloc(timeSeries->data->length*sizeof(REAL4));
  timeSeries->deltaT=Tsft/timeSeries->data->length;

  timeSeries->f0=f_min;

  return 0;
}


/*Allocates space for timeseries */
int prepare_fvec(LALStatus* status) {

  INT4 len,len2;

  len=timeSeries->data->length;
  len2=(INT4)(len/2)+1;


  /* Create vector to hold signal frequency series */
  LALCCreateVector(status, &fvec, (UINT4)len2);

  /* if (sigma > 0.0 || sigma < 0.0) */
  /* LALCCreateVector(status, &fvecn, (UINT4)len2); */

  /* Compute measured plan for FFTW */
  LALCreateForwardRealFFTPlan(status, &pfwd, (UINT4)len, 0);

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

/* computes power in generated signal */
int compute_power(void) {

  unsigned int i;

  for (i = 0; i < timeSeries->data->length; i++)
    {
      power+=timeSeries->data->data[i]*timeSeries->data->data[i];
    }
  return 0;
}


  /*   %strcpy(filename,inDataFilename); */
/*  fp=fopen(timestampsname,"r");*/
/*  if (fp==NULL) {*/
  /*  fprintf(stderr,"Unable to find file %s\n",timestampsname);*/
  /*  return 1;*/
/*  } */
/*  timestamps=(LIGOTimeGPS *)LALMalloc(nTsft*sizeof(LIGOTimeGPS)); */
/*  SSBtimestamps=(LIGOTimeGPS *)LALMalloc(nTsft*sizeof(LIGOTimeGPS));*/

/*  for (i=0;i<nTsft;i++){*/
 /*   r=fscanf(fp,"%d  %d\n", &timestamps[i].gpsSeconds, &timestamps[i].gpsNanoSeconds);*/
  /*  if ( r !=2 ) {*/
   /*   fprintf(stderr,"Unable to read datum # %d\n",i);*/
    /*  fprintf(stderr,"from file %s\n",timestampsname);*/
    /*  return 1; */
  /*  } */
/*  } */

/*  fclose(fp);*/
/*  return 0;*/

/*} */


/*reads timestamps file and fills-in timestamps vector*/
int read_timestamps() {

  FILE *fp;
  int r;
  int i;

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
  const char *filename="AmplMod.dat";
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
    LALComputeDetAMResponse(status, &amresp, &detectorandsource, &timestamps[i]);
    fprintf(fp,"%f  %f\n",amresp.plus,amresp.cross);

  }

  return 0;
}




int make_filelist(void) {

  unsigned int filenum=0;
  char command[256];
  glob_t globbuf;

  strcpy(command,noisedir);
  strcat(command,"*SFT*");

  globbuf.gl_offs = 1;
  glob(command, GLOB_ERR, NULL, &globbuf);

  /* read file names -- MUST NOT FORGET TO PUT ERROR CHECKING IN HERE !!!! */
  while (filenum< globbuf.gl_pathc)
    {
      strcpy(filelist[filenum],globbuf.gl_pathv[filenum]);
      filenum++;
      if (filenum > MAXFILES)
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
  unsigned int i;

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
    syserror("Unable to open file %s\n", filelist[iSFT]);
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
  if ((f_min*Tsft < header.firstfreqindex) ||
      ((f_min*Tsft+fvec->length-1) > header.firstfreqindex + header.nsamples)){
    error("Frequency band of noise data out of range !\n");
    return 1;
  }

  /* seek to position */
  if (0 != fseek(fp, (UINT4)(f_min*Tsft - header.firstfreqindex)*2*sizeof(REAL4), SEEK_CUR)){
    syserror("file too short (could'n fssek to position\n");
    return 1;
  }


  /* calculate size of the data in bytes */
  /* datasize = (f2ind - f1ind) * head.dsize * 2;*/

  /* allocate storage space if needed */
  if (!fvecn) {
    INT4 len=timeSeries->data->length;
    INT4 len2=(INT4)(len/2)+1;
    LALCCreateVector(status, &fvecn, (UINT4)len2);
  }

  /* read the data */
  /* printf("fvac->length is %d\n",fvecn->length); */
  /* printf("header.nsamples is %d\n",header.nsamples); */
  if (1 != fread(fvecn->data,(UINT4)(fvecn->length-1)*sizeof(COMPLEX8),1,fp)) {
      syserror("Could not read the data \n");
    return 1;
  }

  fclose(fp);

  norm=((REAL4)(fvec->length-1)*1.0/((REAL4)header.nsamples));

  for (i = 0; i < fvec->length; ++i) {
    fvec->data[i].realf_FIXME += scale*crealf(fvecn->data[i])*norm;
    fvec->data[i].imagf_FIXME += scale*cimagf(fvecn->data[i])*norm;
  }

  return 0;
}








int write_SFTS(int iSFT){


  FILE *fp;
  REAL4 rpw,ipw;
  char filename[256], filenumber[16];
  int errorcode;
  unsigned int i;

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
  header.firstfreqindex=(INT4)(f_min*Tsft+0.5);
  header.nsamples=fvec->length-1;


  /* write header */
  errorcode=fwrite((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1){
    printf("Error in writing header into file!\n");
    return 1;
  }

  for (i=0;i<fvec->length-1;i++){

    rpw=crealf(fvec->data[i]);
    ipw=cimagf(fvec->data[i]);

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




int write_timefile(int iSFT){

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
  header.firstfreqindex=(INT4)(f_min*Tsft);
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

  /* scan through the list of arguments on the command line
     and get the input data filename*/

  while (!errflg && ((c = getopt(argc, argv,"i:n:hb:D:t:o:I:"))!=-1))
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
    case 'D':
      /* Name of directory where noise data is held */
      noisedir=optarg;
      break;
    case 'I':
      /* detector site*/
      siteinfo=atol(optarg);
      if (siteinfo > 2 ) {
	      fprintf(stderr,"Chosen detector: GEO600\n");
      }
      break;
    case 'h':
      /* print usage/help message */
      fprintf(stderr,"Arguments are:\n");
      fprintf(stderr,"-i\tCharacter String\t( Name of input parameter file)\n");
      fprintf(stderr,"-n\tCharacter String\t( Basefilename of output data files)\n");
      fprintf(stderr,"-I\tInteger number\t\t( IFO: default: GEO, 1 LLO, 2 LHO)\n");
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
  /*f_min REAL8*/
  r=fscanf(fp, "%lf %s\n",&f_min,dmp);
  if ( r !=2 ) {
    fprintf(stderr,"Unable to assign f_min from %s\n",filename);
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
   if ( system(command) ) error("\nsystem(%s) returned non-zero status!\n", command );
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
   if ( system(command) ) error("\nsystem(%s) returned non-zero status!\n", command );
    return 1;
  }
  /*psi REAL4*/
  r=fscanf(fp, "%f %s\n",&genTayParams.psi,dmp);
  if ( r !=2 ) {
	          char command[1024];
    fprintf(stderr,"Unable to assign psi from %s\n",filename);
    sprintf(command, "cat %s 1>&2\n", filename);
   if ( system(command) ) error("\nsystem(%s) returned non-zero status!\n", command );

    return 1;
  }
  /*phi0 REAL8*/
  r=fscanf(fp, "%lf %s\n",&genTayParams.phi0,dmp);
  if ( r !=2 ) {
	          char command[1024];
    fprintf(stderr,"Unable to assign phi0 from %s\n",filename);
   sprintf(command, "cat %s 1>&2\n", filename);
   if ( system(command) ) error("\nsystem(%s) returned non-zero status!\n", command );

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

  /*Orbital Semi-Major Axis REAL8*/
  r=fscanf(fp, "%lf %s\n",&SemiMajorAxis,dmp);
  if ( r !=2 ) {
	          char command[1024];
    fprintf(stderr,"Unable to assign orbital radius from %s\n",filename);
    sprintf(command, "cat %s 1>&2\n", filename);
    if ( system(command) ) error("\nsystem(%s) returned non-zero status!\n", command );


    return 1;
  }
  /*Orbital Period REAL8*/
  r=fscanf(fp, "%lf %s\n",&OrbitalPeriod,dmp);
  if ( r !=2 ) {
	          char command[1024];
    fprintf(stderr,"Unable to assign orbital angular velocity from %s\n",filename);
    sprintf(command, "cat %s 1>&2\n", filename);
    if ( system(command) ) error("\nsystem(%s) returned non-zero status!\n", command );

    return 1;
  }

  /* Time of observed periapse passage in the SSB frame (sec) INT8 */
  r=fscanf(fp, "%d %s\n",&TperiapseSSB.gpsSeconds,dmp);
  if ( r !=2 ) {
	          char command[1024];
    fprintf(stderr,"Unable to assign orbital initial phase from %s\n",filename);
    sprintf(command, "cat %s 1>&2\n", filename);
    if ( system(command) ) error("\nsystem(%s) returned non-zero status!\n", command );

    return 1;
  }

  /*Time of observed periapse passage in the SSB frame (nanosec) INT8*/
  r=fscanf(fp, "%d %s\n",&TperiapseSSB.gpsNanoSeconds,dmp);
  if ( r !=2 ) {
	          char command[1024];
    fprintf(stderr,"Unable to assign orbital initial phase from %s\n",filename);
    sprintf(command, "cat %s 1>&2\n", filename);
    if ( system(command) ) error("\nsystem(%s) returned non-zero status!\n", command );

    return 1;
  }

  /*Argument of periapse REAL8 */
  r=fscanf(fp, "%lf %s\n",&ArgPeriapse,dmp);
  if ( r !=2 ) {
	          char command[1024];
    fprintf(stderr,"Unable to assign orbital initial phase from %s\n",filename);
    sprintf(command, "cat %s 1>&2\n", filename);
    if ( system(command) ) error("\nsystem(%s) returned non-zero status!\n", command );

    return 1;
  }

  /*Orbital eccentricity REAL8*/
  r=fscanf(fp, "%lf %s\n",&OrbitalEccentricity,dmp);
  if ( r !=2 ) {
	          char command[1024];
    fprintf(stderr,"Unable to assign orbital initial phase from %s\n",filename);
    sprintf(command, "cat %s 1>&2\n", filename);
    if ( system(command) ) error("\nsystem(%s) returned non-zero status!\n", command );

    return 1;
  }

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

  /* update global variables f_min, B */
  imin=floor(f_min*Tsft);
  f_min=imin/Tsft;
  /*idealized heterodyning in a band of amplitude B*/
  nsamples=2*ceil(B*Tsft);
  B=(nsamples/2)/Tsft;


  /*and return */
  return errflg;
}


/*sets up input parameters for signal generation*/
int SetupSigGenParams(void){

  REAL8 Ecc;
  REAL8 OneMEcc;
  REAL8 OnePEcc;
  REAL8 TperiTrue;
  REAL8 TperiObs;
  REAL8 correction;
  LIGOTimeGPS TperiapseTrue;

  Ecc=OrbitalEccentricity;
  OneMEcc=1.0-Ecc;
  OnePEcc=1.0+Ecc;

  /* we need to convert the observed time of periapse passage to the "true" time of periapse passage */
  /* this correction is due to the light travel time from binary barycenter to the source at periapse */
  /* it is what Teviets codes require */
  /* IMPORTANT - note that we also define the spin down epoch as the same time as the orbit epoch for simplicity */

  /* printf("Input Tperiapse is %d.%d\n",TperiapseSSB.gpsSeconds,TperiapseSSB.gpsNanoSeconds); */
     /* printf("First SSB timestamp is %d.%d\n",SSBtimestamps[0].gpsSeconds,SSBtimestamps[0].gpsNanoSeconds); */
  correction=SemiMajorAxis*OneMEcc*sin(ArgPeriapse);
  /* printf("Correction is %lf\n",correction); */
  TimeToFloat(&TperiObs,&TperiapseSSB);
  /* printf("TperiObs is %lf\n",TperiObs); */
  TperiTrue=TperiObs-correction;
  /* printf("TperiTrue is %lf\n",TperiTrue); */
  FloatToTime(&TperiapseTrue,&TperiTrue);

  /* convert our input orbital parameters to required orbital parameters */
  genTayParams.oneMinusEcc=OneMEcc;
  genTayParams.rPeriNorm=SemiMajorAxis*OneMEcc;
  genTayParams.angularSpeed=(LAL_TWOPI/OrbitalPeriod)*sqrt(OnePEcc/(OneMEcc*OneMEcc*OneMEcc));
  genTayParams.omega=ArgPeriapse;

  /* set required epochs */
  genTayParams.epoch.gpsSeconds = SSBtimestamps[0].gpsSeconds;
  genTayParams.epoch.gpsNanoSeconds = SSBtimestamps[0].gpsNanoSeconds;
  genTayParams.spinEpoch.gpsSeconds = TperiapseTrue.gpsSeconds;   /* included for binary motion */
  genTayParams.spinEpoch.gpsNanoSeconds = TperiapseTrue.gpsNanoSeconds;  /* included for binary motion */
  genTayParams.orbitEpoch.gpsSeconds = TperiapseTrue.gpsSeconds;    /* included for binary motion */
  genTayParams.orbitEpoch.gpsNanoSeconds = TperiapseTrue.gpsNanoSeconds;  /* included for binary motion */
  /* printf("OrbitEpoch is %d.%d\n",genTayParams.orbitEpoch.gpsSeconds,genTayParams.orbitEpoch.gpsNanoSeconds); */


  /* define length of output */
  genTayParams.length=(LTT+Tsft+timestamps[nTsft-1].gpsSeconds-timestamps[0].gpsSeconds)/genTayParams.deltaT;
  /* printf("Orbit time 0 is %d %d \n",genTayParams.orbitEpoch.gpsSeconds,genTayParams.orbitEpoch.gpsNanoSeconds); */
  return 0;

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


/* Internal routines */
static void TimeToFloat(REAL8 *f, LIGOTimeGPS *tgps)
{
  INT4 x, y;

  x=tgps->gpsSeconds;
  y=tgps->gpsNanoSeconds;
  *f=(REAL8)x+(REAL8)y*1.e-9;
}

static void FloatToTime(LIGOTimeGPS *tgps, REAL8 *f)
{
  REAL8 temp0, temp2, temp3;
  REAL8 temp1, temp4;

  temp0 = floor(*f);     /* this is tgps.S */
  temp1 = (*f) * 1.e10;
  temp2 = fmod(temp1, 1.e10);
  temp3 = fmod(temp1, 1.e2);
  temp4 = (temp2-temp3) * 0.1;

  tgps->gpsSeconds = (INT4)temp0;
  tgps->gpsNanoSeconds = (INT4)temp4;
}
