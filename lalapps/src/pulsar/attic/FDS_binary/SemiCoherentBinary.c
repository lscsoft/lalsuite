/*  Copyright (C) 2010 Chris Messenger
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

/**
 * \author C.Messenger
 * \ingroup pulsarApps
 * \file
 * \brief
 * This code is designed to compute the Bayes factor for a semi-coherent analysis
 * of input SFT data specific to searching for continuous signals in a binary system.
 *
 * It generates likelihood samples from a coarse grid of templates placed on each SFT and
 * combines them using a fine binary template band.  The parameter space is integrated over
 * and a Bayes factor is produced.
 *
 */

/***********************************************************************************************/
/* includes */
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <gsl/gsl_interp.h>        /* needed for the gsl interpolation */
#include <gsl/gsl_spline.h>        /* needed for the gsl interpolation */
#include <gsl/gsl_rng.h>           /* for random number generation */ 
#include <gsl/gsl_randist.h>       /* for random number generation */ 
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf_log.h>        /* for log computation */
#include <lal/TimeSeries.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lal/SFTutils.h>
#include <lal/SFTfileIO.h>
#include <lal/ComplexFFT.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lalapps.h>

/***********************************************************************************************/
/* some global constants */

#define STRINGLENGTH 256              /* the length of general string */
#define APIDLENGTH 5                  /* the length of an APID string */
#define LONGSTRINGLENGTH 1024         /* the length of general string */
#define NFREQMAX 4                    /* the max dimensionality of the frequency derivitive grid */
#define NBINMAX 4                     /* the number of binary parameter dimensions */
#define NBINS 4                       /* the number of bins to add to each side of the fft for safety */
#define WINGS_FACTOR 2                /* the safety factor in reading extra frequency from SFTs */
#define PCU_AREA 0.13                 /* the collecting area of a single PCU in square metres */
#define DEFAULT_SOURCE "SCOX1"        /* the default source name */
#define AMPVECLENGTH 25               /* the fixed number of amplitude values to sample */
#define NLOGLUT 64                    /* the number of elements in the log LUT */
#define NBESSELLUT 256                 /* the number of elements in the bessel function LUT */

/***********************************************************************************************/
/* some useful macros */

#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

/***********************************************************************************************/
/* define internal structures */

/**
 * A single parameter prior pdf
 */
typedef struct { 
  REAL8Vector *logpriors;           /**< vector that stores the log of the prior pdf */
  REAL8 logdelta;                   /**< the log of the spacing */
  BOOLEAN gaussian;                 /**< are we using a Gaussian prior on this parameter */
} REAL8Priors;

/**
 * A vector of prior pdfs for many dimensions
 */
typedef struct { 
  REAL8Priors *data;                /**< points to the prior data */
  UINT4 ndim;                       /**< the dimensionality of the prior space */
} REAL8PriorsVector;

/**
 * A single parameter dimensions boundaries
 */
typedef struct { 
  REAL8 min;                        /**< the parameter space minimum */
  REAL8 max;                        /**< the parameter space maximium */
  REAL8 mid;                        /**< the parameter space mid point (where a Gaussian prior is centered) */
  REAL8 sig;                        /**< the one-sigma uncertainty on the parameter */
  REAL8 span;                       /**< the parameter space span */
  BOOLEAN gaussian;                 /**< are we using a Gaussian prior on this parameter */
  CHAR name[LALNameLength];         /**< string containing the name of the dimension */
} REAL8Dimension;

/**
 * A vector of parameter space boundary information
 */
typedef struct { 
  REAL8Dimension *data;             /**< the boundaries, span, etc for a single dimension */
  UINT4 ndim;                       /**< the number of dimensions */
} REAL8Space;

/**
 * Stores the gridding parameters for a single dimension
 */
typedef struct { 
  REAL8 min;                        /**< the starting points of the grid */
  REAL8 delta;                      /**< the grid spacings */
  REAL8 oneoverdelta;               /**< the inverse of the spacing */
  UINT4 length;                     /**< the number of templates in each dimension */
  CHAR name[LALNameLength];         /**< string containing the name of the dimension */
} Grid;

/**
 * Stores the current location in a hyper-cubic parameter space
 */
typedef struct { 
  REAL8 *x;                         /**< the location in parameter space */
  INT4 *idx;                        /**< the index of each dimension for this template */
  UINT4 ndim;                       /**< the dimension of the parameter space */
  UINT4 currentidx;                 /**< the current index value of the template */
} Template;

/**
 * Stores the gridding parameters for a hypercubic grid of templates
 */
typedef struct { 
  Grid *grid;                       /**< stores the parameters defining a single dimension */
  UINT4 ndim;                       /**< the number of dimensions */
  UINT4 *prod;                      /**< internal variable used to store the size of sub-dimensions */
  UINT4 max;                        /**< the maximum (total) number of templates */
  REAL8 mismatch;                   /**< the mismatch */
} GridParameters;

/**
 * Stores the parameters of an injection
 */
typedef struct { 
  Template temp;                    /**< stores the parameters of the signal */
  REAL8 amp;                        /**< the injected amplitude */
} InjectionParameters;

/**
 * contains information regarding the search parameter space
 */
typedef struct { 
  REAL8Space *space;                /**< stores the parameter space boundaries */
  REAL8Dimension *ampspace;         /**< the amplitude space */
  GridParameters *gridparams;       /**< stores the grid */
  Grid *ampgrid;                    /**< stores a 1D grid on amplitude */ 
  REAL8PriorsVector *priors;        /**< stores the priors on the paramaters */
  REAL8Priors *amppriors;           /**< the priors on the amplitude */
  LIGOTimeGPS epoch;                /**< the start time of the entire observation */
  REAL8 span;                       /**< the span of the entire observation */
  REAL8 tseg;                       /**< the coherent time */
  CHAR source[LALNameLength];       /**< the name of the source */
  InjectionParameters *inj;         /**< stores the injected signal parameters (if any) */
} ParameterSpace;

/*****************************************************************************************/

/**
 * Stores the gridding parameters for a hypercubic grid of templates
 */
typedef struct { 
  GridParameters **segment;         /**< stores the parameters defining a single dimension */
  UINT4 length;                     /**< the number of segments */
} GridParametersVector;

/**
 * Stores segment parameters
 */
typedef struct { 
  INT4Vector *npcus;                /**< a vector of PCUs */
  REAL8Vector *dt;                  /**< a vector of sampling times */
} SegmentParams;

/**
 * Stores parameters useful for the efficient calcualtion of the likelihood
 */
typedef struct { 
  REAL8 logsqrtP;                   /**< intermediate variable for the phase and amp marginalised likelihood calculation */
  REAL8 PQ;                         /**< intermediate variable for the phase and amp marginalised likelihood calculation */
  REAL8Vector *alphasqY;            /**< a vector of a variable computed for each amplitude value */
  REAL8Vector *alphaX;              /**< another vector of a variable computed for each amplitude value */
} LikelihoodParams;

/**
 * Stores parameters useful for the efficient calcualtion of the likelihood
 */
typedef struct { 
  LikelihoodParams *data;           /**< a vector of likelihood parameter structures */
  UINT4 length;                     /**< the length of the vector */
  REAL8Vector *logLratio_phase;     /**< a temporary storage for the fixed amplitude logL ratio for each alpha value */
  REAL8Vector *power;               /**< stores the power for a single template for all segments */
  REAL8Vector *alphasqsumY;         /**< the sum of Y over the segments multiplied by the amplitude */ 
  gsl_interp_accel *log_acc;        /**< gsl interpolation structure for log LUT */
  gsl_spline *log_spline;           /**< gsl interpolation structure for log LUT */
  gsl_interp_accel *logbesselI0_acc;   /**< gsl interpolation structure for bessel LUT */
  gsl_spline *logbesselI0_spline;      /**< gsl interpolation structure for bessel LUT */
} LikelihoodParamsVector;

/**
 * Stores the results of a Bayesian posterior integration (Bayes factor, evidence, posteriors, etc...)
 */
typedef struct { 
  REAL8 logBayesFactor_phaseamp;                /**< the log Bayes factor for phase and amplitude marginalised per segment */
  REAL8 logBayesFactor_phase;                   /**< the log Bayes factor for phase marginalised per segment */
  REAL8Vector *logBayesFactor_phaseamp_vector;  /**< the log Bayes factor for each segment individually */
 /*  REAL8Vector *logBayesFactor_phase_vector;     /\**< the log Bayes factor for fixed amplitude for each segment individually *\/ */
  REAL8Vector **logposteriors_phaseamp;         /**< the output log posteriors for phase and amplitude marginalised per segment */
  REAL8Vector **logposteriors_phase;            /**< the output log posteriors for phase marginalised per segment */
  REAL8Vector *logposterior_amp;                /**< the log posterior for the fixed amplitude parameter */
  GridParameters *gridparams;                   /**< the grid used for the marginalisation */
  Grid *ampgrid;                                /**< the grid used for amplitude marginalisation */
  UINT4 ndim;                                   /**< the dimensionality of the space */
  LIGOTimeGPS *epoch;                           /**< the epochs of each segment */
  UINT4 nsegments;                              /**< the number of segments used */
} BayesianProducts;

/**
 * Storage for the demodulated power from a single segment
 */
typedef struct { 
  REAL4Vector *data;                /**< pointer to the power data stored sequentially */
  REAL8 r;                          /**< the estimated noise background */
  UINT4 npcus;                      /**< the number of PCUs */
  REAL8 dt;                         /**< the original sampling time */
  LIGOTimeGPS epoch;                /**< the epoch of the segment */
  REAL8 duration;                   /**< the duration of the segment */
  GridParameters *gridparams;       /**< the grid on which the power was computed */
} REAL4DemodulatedPower;

/**
 * Storage for the demodulated power
 */
typedef struct { 
  REAL4DemodulatedPower **segment;  /**< pointer to a set of REAL4VectorArrays */
  UINT4 length;                     /**< the number of segments */
} REAL4DemodulatedPowerVector;

/**
 * An array of COMPLEX8TimeSeries
 */
typedef struct { 
  COMPLEX8TimeSeries **data;        /**< pointer to a set of COMPLEX8TimeSeries */
  UINT4 length;                     /**< the number of vectors */
} COMPLEX8TimeSeriesArray;

/**
 * A structure that stores user input variables
 */
typedef struct { 
  BOOLEAN help;		            /**< trigger output of help string */
  CHAR *sftbasename;                /**< basename of input SFT files */
  CHAR *outputdir;                  /**< the output directory */
  CHAR *source;                     /**< the name of the source */
  REAL8 freq;                       /**< the starting frequency */
  REAL8 freqband;                   /**< the search band width */
  REAL8 orbperiod;                  /**< the central orbital period value */
  REAL8 deltaorbperiod;             /**< the uncertainty on the orbital period */
  REAL8 asini;                      /**< the central orbital semi-major axis value */
  REAL8 deltaasini;                 /**< the uncertainty on the orbital semi-major axis value */
  REAL8 tasc;                       /**< the central time of ascension value */
  REAL8 deltatasc;                  /**< the uncertainty on the central time of ascension value */
  REAL8 nsig;                       /**< the width (in sigmas) of the Gaussian priors */
  REAL8 mismatch;                   /**< the grid mismatch */      
  REAL8 sigalpha;                   /**< the amplitude prior sigma */
  INT4 gpsstart;                    /**< the min GPS time to include */
  INT4 gpsend;                      /**< the max GPS time to include */
  BOOLEAN gaussianpriors;           /**< flag for using Gaussian priors on the orbital parameters */
  CHAR *obsid_pattern;              /**< the OBS ID substring */
  INT4 seed;                        /**< fix the random number generator seed */
  REAL8 inject_amplitude;           /**< the amplitude of the injected signal */
  BOOLEAN fixedamp;                 /**< use the fixed amplitude model */
  BOOLEAN version;	            /**< output version-info */
} UserInput_t;

/***********************************************************************************************/
/* global variables */
extern int vrbflg;	 	/**< defined in lalapps.c */

/* parameters for bessel function calculation (taken from Abramowitz and Stegun P.367) */
UINT4 LEN_BESSCO_HIGH = 9;
REAL8 BESSCO_HIGH[] = {0.39894228,0.01328592,0.00225319,-0.00157565,0.00916281,-0.02057706,0.02635537,-0.01647633,0.00392377};
UINT4 LEN_BESSCO_LOW = 7;
REAL8 BESSCO_LOW[] = {1.0,3.5156229,3.0899424,1.2067492,0.2659732,0.0360768,0.0045813};

/***********************************************************************************************/
/* define functions */
int main(int argc,char *argv[]);
int XLALReadUserVars(int argc,char *argv[],UserInput_t *uvar, CHAR **clargs);
int XLALDefineBinaryParameterSpace(REAL8Space **space,LIGOTimeGPS epoch, REAL8 span,UserInput_t *uvar); 
int XLALComputeAmplitudeParams(REAL8Dimension **ampspace,Grid **ampgrid,REAL8Priors **amppriors,REAL8 ampsigma); 
int XLALReadSFTs(SFTVector **sfts,SegmentParams **segparams,CHAR *sftbasename, REAL8 freq, REAL8 freqband, INT4 gpsstart, INT4 gpsend,CHAR *obsid_pattern);
int XLALComputeFreqGridParamsVector(GridParametersVector **freqgridparams,REAL8Space *pspace, SFTVector *sftvec, REAL8 mu);
int XLALComputeFreqGridParams(GridParameters **freqgridparams,REAL8Space *pspace, REAL8 tmid,REAL8 tsft, REAL8 mu);
int XLALSFTVectorToCOMPLEX8TimeSeriesArray(COMPLEX8TimeSeriesArray **dstimevec, SFTVector *sftvec);
int XLALComputeDemodulatedPower(REAL4DemodulatedPower **power,COMPLEX8TimeSeries *time,GridParameters *gridparams);
int XLALComputeDemodulatedPowerVector(REAL4DemodulatedPowerVector **power,COMPLEX8TimeSeriesArray *time,GridParametersVector *gridparams);
int XLALEstimateBackgroundFlux(REAL8Vector **background, SegmentParams *segparams, SFTVector *sftvec);
int XLALComputeBinaryGridParams(GridParameters **binarygridparams,REAL8Space *space,REAL8 T,REAL8 DT,REAL8 mu);
int XLALComputeBinaryPriors(REAL8PriorsVector **priors,REAL8Space *space,GridParameters *gridparams);
int XLALComputeBayesFactor(BayesianProducts **Bayes,REAL4DemodulatedPowerVector *power,ParameterSpace *pspace,REAL8 sigalpha);
int XLALGetNextBinaryTemplate(Template **temp,GridParameters *gridparams);
int XLALComputeBinaryFreqDerivitives(Template *fdots,Template *bintemp,REAL8 tmid);
REAL8 XLALComputePhaseAmpMargLogLRatio(REAL8 X,LikelihoodParams *Lparams);
REAL8 XLALComputePhaseAmpMargLogLRatioLUT(REAL8 X,LikelihoodParams *Lparams,gsl_interp_accel *logbesselI0_acc,gsl_spline *logbesselI0_spline);
int XLALSetupLikelihood(LikelihoodParamsVector **Lparamsvec,BayesianProducts **Bayes,REAL4DemodulatedPowerVector *power,GridParameters *binarygrid,Grid *ampgrid,REAL8 sigalpha);
REAL8 XLALLogBesselI0(REAL8 z);
REAL8 XLALLogSumExp(REAL8 logx,REAL8 logy);
REAL8 XLALLogSumExpLUT(REAL8 logx,REAL8 logy,gsl_interp_accel *log_acc,gsl_spline *log_spline);    
int XLALFreeParameterSpace(ParameterSpace *pspace);
int XLALFreeREAL4DemodulatedPowerVector(REAL4DemodulatedPowerVector *power);
int XLALFreeBayesianProducts(BayesianProducts *Bayes);
int XLALOutputBayesResults(CHAR *outputdir,BayesianProducts *Bayes,ParameterSpace *pspace,CHAR *clargs,CHAR *obsid_pattern);
int XLALAddBinarySignalToSFTVector(SFTVector **sftvec,ParameterSpace *pspace,REAL8 inject_amplitude,INT4 seed);
int XLALInitgslrand(gsl_rng **gslrnd,INT8 seed);
int XLALComputePhaseMargLogLRatio(REAL8Vector *logLratio,REAL8 X,LikelihoodParams *Lparams);
int XLALComputePhaseMargLogLRatioLUT(REAL8Vector *logLratio,REAL8 X,LikelihoodParams *Lparams,gsl_interp_accel *logbesselI0_acc,gsl_spline *logbesselI0_spline);
int XLALComputePhaseMargLogLRatioVectorLUT(REAL8Vector *logLratio,REAL8Vector *power,LikelihoodParamsVector *Lparamsvec,gsl_interp_accel *logbesselI0_acc,gsl_spline *logbesselI0_spline);

/***********************************************************************************************/
/* empty initializers */
UserInput_t empty_UserInput;
ParameterSpace empty_ParameterSpace;

/**
 * The main function of semicoherentbinary.c
 */
int main( int argc, char *argv[] )  {

  UserInput_t uvar = empty_UserInput;           /* user input variables */
  CHAR *clargs = NULL;                          /* store the command line args */
  SFTVector *sftvec = NULL;                     /* stores the input SFTs */
  SegmentParams *segparams = NULL;              /* stores the number of PCUs and sampling time per SFT */
  REAL8Vector *background = NULL;               /* estimates of the background for each SFT */
  ParameterSpace pspace = empty_ParameterSpace; /* the search parameter space */
  COMPLEX8TimeSeriesArray *dstimevec = NULL;    /* contains the downsampled inverse FFT'd SFTs */
  REAL4DemodulatedPowerVector *power = NULL;    /* contains the demodulated power for all SFTs */
  GridParametersVector *freqgridparams = NULL;  /* the coherent grid on the frequency derivitive parameter space */
  BayesianProducts *Bayes = NULL;               /* the Bayesian results */
  REAL8 fmin_read,fmax_read,fband_read;         /* the range of frequencies to be read from SFTs */
  UINT4 i;                                      /* counter */

  vrbflg = 1;	                        /* verbose error-messages */

  /* turn off default GSL error handler */
  gsl_set_error_handler_off();

  /* setup LAL debug level */
  LogSetLevel(lalDebugLevel);

  /* register and read all user-variables */
  if (XLALReadUserVars(argc,argv,&uvar,&clargs)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadUserVars() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : read in uservars\n",__func__);

  /* make crude but safe estimate of the bandwidth required for the source */
 /*  fmin_read = pspace.space->data[0].min - WINGS_FACTOR*pspace.space->data[0].min*pspace.space->data[1].max*pspace.space->data[3].max; */
/*   fmax_read = pspace.space->data[0].max + WINGS_FACTOR*pspace.space->data[0].max*pspace.space->data[1].max*pspace.space->data[3].max; */
  {
    REAL8 wings = LAL_TWOPI*(uvar.asini + uvar.deltaasini*uvar.nsig)/(uvar.orbperiod - uvar.deltaorbperiod*uvar.nsig);
    fmin_read = uvar.freq - WINGS_FACTOR*uvar.freq*wings;
    fmax_read = uvar.freq + uvar.freqband + WINGS_FACTOR*(uvar.freq + uvar.freqband)*wings;
    fband_read = fmax_read - fmin_read;
    LogPrintf(LOG_DEBUG,"%s : reading in SFT frequency band [%f -> %f]\n",__func__,fmin_read,fmax_read); 
  }
 
  /**********************************************************************************/
  /* READ THE SFT DATA */
  /**********************************************************************************/
  
  /* load in the SFTs - also fill in the segment parameters structure */
  if (XLALReadSFTs(&sftvec,&segparams,uvar.sftbasename,fmin_read,fband_read,uvar.gpsstart,uvar.gpsend,uvar.obsid_pattern)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadSFTs() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : read in SFTs\n",__func__);
  
  /* define SFT length and the start and span of the observations plus the definitive segment time */
  pspace.tseg = 1.0/sftvec->data[0].deltaF;
  memcpy(&(pspace.epoch),&(sftvec->data[0].epoch),sizeof(LIGOTimeGPS));
  pspace.span = XLALGPSDiff(&(sftvec->data[sftvec->length-1].epoch),&(sftvec->data[0].epoch)) + pspace.tseg;
  sprintf(pspace.source,"%s",uvar.source);
  LogPrintf(LOG_DEBUG,"%s : SFT length = %f seconds\n",__func__,pspace.tseg);
  LogPrintf(LOG_DEBUG,"%s : entire dataset starts at GPS time %d contains %d SFTS and spans %.0f seconds\n",__func__,pspace.epoch.gpsSeconds,sftvec->length,pspace.span);
  
  /**********************************************************************************/
  /* DEFINE THE BINARY PARAMETER SPACE */
  /**********************************************************************************/
  
  /* register and read all user-variables */
  if (XLALDefineBinaryParameterSpace(&(pspace.space),pspace.epoch,pspace.span,&uvar)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALDefineBinaryParameterSpace() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : defined binary parameter prior space\n",__func__);

  /**********************************************************************************/
  /* DEFINE THE AMPLITUDE PARAMETERS IF USING A FIXED AMPLITUDE SIGNAL MODEL */
  /**********************************************************************************/
  
  if (uvar.fixedamp) {
    if (XLALComputeAmplitudeParams(&(pspace.ampspace),&(pspace.ampgrid),&(pspace.amppriors),uvar.sigalpha)) {
      LogPrintf(LOG_CRITICAL,"%s : XLALComputeAmplitudeParams() failed with error = %d\n",__func__,xlalErrno);
      return 1;
    }
  }

  /**********************************************************************************/
  /* COMPUTE THE FINE GRID PARAMETERS */
  /**********************************************************************************/
  
  /* compute the fine grid on the binary parameters */
  if (XLALComputeBinaryGridParams(&(pspace.gridparams),pspace.space,pspace.span,pspace.tseg,uvar.mismatch)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALComputeBinaryGridParams() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : computed the binary parameter space grid\n",__func__);
 
  /**********************************************************************************/
  /* COMPUTE THE BINARY PARAMETER SPACE PRIORS */
  /**********************************************************************************/
  
  /* compute the priors on the binary parameters */
  if (XLALComputeBinaryPriors(&(pspace.priors),pspace.space,pspace.gridparams)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALComputeBinaryPriors() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : computed the binary parameter space priors\n",__func__);
  
  /**********************************************************************************/
  /* ADD A SIMULATED SIGNAL TO THE SFT DATA */
  /**********************************************************************************/

  if (XLALUserVarWasSet(&(uvar.inject_amplitude))) {
    if (XLALAddBinarySignalToSFTVector(&sftvec,&pspace,uvar.inject_amplitude,uvar.seed)) {
      LogPrintf(LOG_CRITICAL,"%s : XLALAddBinarySignalToSFTVector() failed with error = %d\n",__func__,xlalErrno);
      return 1;
    }
    LogPrintf(LOG_DEBUG,"%s : added a simulated signal to the SFTs\n",__func__);
  }

  /**********************************************************************************/
  /* ESTIMATE BACKGROUND NOISE FROM SFTS */
  /**********************************************************************************/
  
  /* compute the background noise using the sfts */
  if (XLALEstimateBackgroundFlux(&background,segparams,sftvec)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALEstimateBackgroundFlux() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : estimated the background noise from the SFTs\n",__func__);
 
  /**********************************************************************************/
  /* CONVERT ALL SFTS TO DOWNSAMPLED TIMESERIES */
  /**********************************************************************************/
  
  /* convert sfts to downsample dtimeseries */
  if (XLALSFTVectorToCOMPLEX8TimeSeriesArray(&dstimevec,sftvec)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALSFTVectorToCOMPLEX8TimeSeriesArray() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : converted SFTs to downsampled timeseries\n",__func__);

  /**********************************************************************************/
  /* COMPUTE THE COARSE GRID ON FREQUENCY DERIVITIVE */
  /**********************************************************************************/

  /* compute the grid parameters for all SFTs */
  if (XLALComputeFreqGridParamsVector(&freqgridparams,pspace.space,sftvec,uvar.mismatch)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALComputeFreqGridParams() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }

  /* free un-needed original SFT vector */
  XLALDestroySFTVector(sftvec);
  LogPrintf(LOG_DEBUG,"%s : Freed the SFT memory\n",__func__);

  /**********************************************************************************/
  /* COMPUTE THE STATISTICS ON THE COARSE GRID */
  /**********************************************************************************/

  /* compute the statistic on the grid */
  if (XLALComputeDemodulatedPowerVector(&power,dstimevec,freqgridparams)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALComputeDemodulatedPower() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : computed the demodulated power\n",__func__);
  
  /* fill in segment parameters */
  for (i=0;i<power->length;i++) {
    power->segment[i]->r = background->data[i];
    power->segment[i]->dt = segparams->dt->data[i];
    power->segment[i]->npcus = segparams->npcus->data[i];
    power->segment[i]->duration = pspace.tseg; 
  }

  /* free un-needed downsampled timeseries */
  for (i=0;i<dstimevec->length;i++) {
    XLALDestroyCOMPLEX8TimeSeries(dstimevec->data[i]);
  }
  XLALFree(dstimevec->data);
  XLALFree(dstimevec);
  LogPrintf(LOG_DEBUG,"%s : freed the downsampled timeseries memory\n",__func__);

  /* free frequency grid - the contents of each segment have been moved to the power structure and are freed later */
  XLALFree(freqgridparams->segment);
  XLALFree(freqgridparams);

  /* free the background estimate */
  XLALDestroyREAL8Vector(background);
  LogPrintf(LOG_DEBUG,"%s : freed the background noise estimate\n",__func__);

  /* free the segment params */
  XLALDestroyREAL8Vector(segparams->dt);
  XLALDestroyINT4Vector(segparams->npcus);
  XLALFree(segparams);
  LogPrintf(LOG_DEBUG,"%s : freed the segment parameters\n",__func__);

  /**********************************************************************************/
  /* INTEGRATE OVER THE FINE GRID TO COMPUTE THE BAYES FACTOR */
  /**********************************************************************************/

  /* compute the Bayes factor and the posterior distributions */
  if (XLALComputeBayesFactor(&Bayes,power,&pspace,uvar.sigalpha)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALComputeBayesFactor() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : computed the BayesFactor and posteriors\n",__func__);
 
  /**********************************************************************************/
  /* OUTPUT RESULTS TO FILE */
  /**********************************************************************************/

  if (XLALOutputBayesResults(uvar.outputdir,Bayes,&pspace,clargs,uvar.obsid_pattern)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALOutputBayesResults() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : output results to file.\n",__func__);
 
  /**********************************************************************************/
  /* CLEAN UP */
  /**********************************************************************************/

  /* clean up the parameter space */
  if (XLALFreeParameterSpace(&pspace)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALFreeParameterSpace() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : freed the parameter space\n",__func__);
 
  /* clean up the demodulated power */
  if (XLALFreeREAL4DemodulatedPowerVector(power)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALFreeREAL4DemodulatedPowerVector() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }  
  LogPrintf(LOG_DEBUG,"%s : freed the demodulated power\n",__func__);

  /* clean up the demodulated power */
  if (XLALFreeBayesianProducts(Bayes)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALFreeBayesianProducts() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }  
  LogPrintf(LOG_DEBUG,"%s : freed the Bayesian results\n",__func__);

  /* Free config-Variables and userInput stuff */
  XLALDestroyUserVars();
  XLALFree(clargs);

  /* did we forget anything ? */
  LALCheckMemoryLeaks();
  LogPrintf(LOG_DEBUG,"%s : successfully checked memory leaks.\n",__func__);

  LogPrintf(LOG_DEBUG,"%s : successfully completed.\n",__func__);
  return 0;
  
} /* end of main */

/**
 * Read in input user arguments
 */
int XLALReadUserVars(int argc,            /**< [in] the command line argument counter */ 
		     char *argv[],        /**< [in] the command line arguments */
		     UserInput_t *uvar,   /**< [out] the user input structure */
		     CHAR **clargs        /**< [out] the command line args string */
		     )
{
  CHAR *version_string;
  INT4 i;

  /* initialise user variables */
  uvar->sftbasename = NULL; 
  uvar->obsid_pattern = NULL;
  uvar->gpsstart = -1;
  uvar->gpsend = -1;
  uvar->mismatch = 0.2;
  uvar->nsig = 0;
  uvar->gaussianpriors = 0;
  uvar->seed = 0;

  /* initialise default source as SCOX1 */
  {
    UINT4 n = strlen(DEFAULT_SOURCE) + 1;
    uvar->source = XLALCalloc(n,sizeof(CHAR));
    snprintf(uvar->source,n,"%s",DEFAULT_SOURCE);
  }

  /* initialise all parameter space ranges to zero */
  uvar->freqband = 0;
  uvar->deltaorbperiod = 0.0;
  uvar->deltaasini = 0.0;
  uvar->deltatasc = 0.0;

  /* ---------- register all user-variables ---------- */
  XLALregBOOLUserStruct(help, 		        'h', UVAR_HELP,     "Print this message");
  XLALregSTRINGUserStruct(sftbasename, 	        'i', UVAR_REQUIRED, "The basename of the input SFT files"); 
  XLALregSTRINGUserStruct(outputdir, 	        'o', UVAR_REQUIRED, "The output directory name"); 
  XLALregSTRINGUserStruct(source, 	        'x', UVAR_OPTIONAL, "The source name (default SCOX1)"); 
  XLALregREALUserStruct(freq,                   'f', UVAR_REQUIRED, "The starting frequency (Hz)");
  XLALregREALUserStruct(freqband,   	        'b', UVAR_OPTIONAL, "The frequency band (Hz)");
  XLALregREALUserStruct(orbperiod,              'P', UVAR_REQUIRED, "The central orbital period value (sec)");
  XLALregREALUserStruct(deltaorbperiod,   	'p', UVAR_OPTIONAL, "The orbital period uncertainty (sec)");
  XLALregREALUserStruct(asini,                  'A', UVAR_REQUIRED, "The central orbital semi-major axis (sec)");
  XLALregREALUserStruct(deltaasini,       	'a', UVAR_OPTIONAL, "The orbital semi-major axis uncertainty (sec)");
  XLALregREALUserStruct(tasc,                   'T', UVAR_REQUIRED, "The central orbital time of ascension (GPS)");
  XLALregREALUserStruct(deltatasc,      	't', UVAR_OPTIONAL, "The orbital time of ascension uncertainty (GPS)");
  XLALregREALUserStruct(nsig,            	'g', UVAR_OPTIONAL, "The width of the Gaussian priors (in sigmas).  If unset flat 1-sig priors are used.");
  XLALregREALUserStruct(mismatch,        	'm', UVAR_OPTIONAL, "The grid mismatch (0->1)");
  XLALregREALUserStruct(sigalpha,        	'z', UVAR_REQUIRED, "The stdev of the zero mean amplitude prior in cnts/s/m^2");
  XLALregINTUserStruct(gpsstart,                's', UVAR_OPTIONAL, "The minimum start time (GPS sec)");
  XLALregINTUserStruct(gpsend,          	'e', UVAR_OPTIONAL, "The maximum end time (GPS sec)");
  XLALregSTRINGUserStruct(obsid_pattern,        'O', UVAR_OPTIONAL, "The observation ID substring to match"); 
  XLALregINTUserStruct(seed,                    'd', UVAR_SPECIAL,  "Fix the random number generator seed");
  XLALregREALUserStruct(inject_amplitude,     	'J', UVAR_SPECIAL,  "The amplitude of the injected signal (in cnts/s/m^2)");
  XLALregBOOLUserStruct(fixedamp,               'F', UVAR_OPTIONAL, "Use the fixed amplitude signal model");
  XLALregBOOLUserStruct(version,                'V', UVAR_SPECIAL,  "Output code version");

  /* do ALL cmdline and cfgfile handling */
  if (XLALUserVarReadAllInput(argc, argv)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALUserVarReadAllInput() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* if help was requested, we're done here */
  if (uvar->help) exit(0);

  if ((version_string = XLALGetVersionString(0)) == NULL) {
    XLALPrintError("XLALGetVersionString(0) failed.\n");
    exit(1);
  }
  
  if (uvar->version) {
    printf("%s\n",version_string);
    exit(0);
  }
  XLALFree(version_string);

  /* set priors flag if the sigma width has been set */
  if (XLALUserVarWasSet(&(uvar->nsig))) {
    uvar->gaussianpriors = 1;
    LogPrintf(LOG_DEBUG,"%s : using Gaussian priors on orbital parameters.\n",__func__);
  }
  else uvar->nsig = 1.0;

  if (uvar->nsig < 0.0) {
    LogPrintf(LOG_CRITICAL,"%s : the user input nsig must > 0.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* put clargs into string */
  *clargs = XLALCalloc(1,sizeof(CHAR));
  for (i=0;i<argc;i++) {
    INT4 len = 2 + strlen(argv[i]) + strlen(*clargs);
    *clargs = XLALRealloc(*clargs,len*sizeof(CHAR));
    strcat(*clargs,argv[i]);
    strcat(*clargs," ");
  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}

/**
 * Computes the binary parameter space boundaries given the user input args
 * For each search dimension we define the min, max, mid, and span of that dimension
 * plus we give each dimension a name, define whether it is to be given a flat or
 * Gaussian prior and specify the sigma of that prior.
 *
 */
int XLALDefineBinaryParameterSpace(REAL8Space **space,                 /**< [out] the parameter space  */
				   LIGOTimeGPS epoch,                  /**< [in] the observation start epoch */
				   REAL8 span,                         /**< [in] the observation span */
				   UserInput_t *uvar                   /**< [in] the user input variables */
				   )
{
  REAL8 midpoint;              /* the midpoint of the observation */
  REAL8 newtasc;               /* shifted value of tasc */
  REAL8 newdeltatasc;          /* updated uncertainty on tasc after shifting */

  /* validate input variables */
  if ((*space) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input REAL8Space boundary structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }  
  if (uvar == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input UserInput_t structure = NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (epoch.gpsSeconds < 0) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, observation epoch < 0.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (span < 0) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, observation span < 0.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  
  /* allocate memory for the parameter space */
  if ( ((*space) = XLALCalloc(1,sizeof(REAL8Space))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  (*space)->ndim = NBINMAX;
  if ( ((*space)->data = XLALCalloc((*space)->ndim,sizeof(REAL8Dimension))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }

  /* define observaton midpoint */
  midpoint = XLALGPSGetREAL8(&epoch) + 0.5*span;

  /* find the closest instance of ascension to the midpoint */
  /* the midpoint is defined in the detector frame and the time of */
  /* ascension is in the SSB frame but we only need to be roughly */
  /* correct in the number of orbits we shift by */
  {
    INT4 n = (INT4)floor(0.5 + (midpoint - uvar->tasc)/uvar->orbperiod);
    newtasc = uvar->tasc + n*uvar->orbperiod;
    newdeltatasc = sqrt(uvar->deltatasc*uvar->deltatasc + n*n*uvar->deltaorbperiod*uvar->deltaorbperiod);
    LogPrintf(LOG_DEBUG,"%s : shifted tasc by %d orbits, uncertainty equals %f sec\n",__func__,n,newdeltatasc);
  }
  
  /* this represents a hyper-cubic parameter space */
  /* we make sure that parameter ranges are consistent i.e asini > 0 etc.. */
  /* frequency */
  snprintf((*space)->data[0].name,LALNameLength,"nu");
  (*space)->data[0].min = uvar->freq;
  (*space)->data[0].max = uvar->freq + uvar->freqband;
  (*space)->data[0].mid = uvar->freq + 0.5*uvar->freqband;
  (*space)->data[0].sig = 0.5*uvar->freqband;
  (*space)->data[0].span = uvar->freqband;
  (*space)->data[0].gaussian = 0;

  /* asini */
  snprintf((*space)->data[1].name,LALNameLength,"asini");
  (*space)->data[1].min = uvar->asini - fabs(uvar->nsig*uvar->deltaasini) > 0.0 ? uvar->asini - fabs(uvar->nsig*uvar->deltaasini) : 0.0;
  (*space)->data[1].max = uvar->asini + fabs(uvar->nsig*uvar->deltaasini);
  (*space)->data[1].mid = uvar->asini;
  (*space)->data[1].sig = uvar->deltaasini;
  (*space)->data[1].span = (*space)->data[1].max - (*space)->data[1].min;
  (*space)->data[1].gaussian = uvar->gaussianpriors;
  
  /* tasc */
  snprintf((*space)->data[2].name,LALNameLength,"tasc");
  (*space)->data[2].min = newtasc - fabs(uvar->nsig*newdeltatasc);
  (*space)->data[2].max = newtasc + fabs(uvar->nsig*newdeltatasc);
  (*space)->data[2].mid = newtasc;
  (*space)->data[2].sig = newdeltatasc;
  (*space)->data[2].span = (*space)->data[2].max - (*space)->data[2].min;
  (*space)->data[2].gaussian = uvar->gaussianpriors;
  
  /* omega */
  snprintf((*space)->data[3].name,LALNameLength,"omega");
  (*space)->data[3].min = LAL_TWOPI/(uvar->orbperiod + fabs(uvar->nsig*uvar->deltaorbperiod));
  (*space)->data[3].max = LAL_TWOPI/(uvar->orbperiod - fabs(uvar->nsig*uvar->deltaorbperiod));
  (*space)->data[3].mid = LAL_TWOPI/uvar->orbperiod;
  (*space)->data[3].sig = (*space)->data[3].mid - (*space)->data[3].min;
  if ((*space)->data[3].max < 0.0) {
    LogPrintf(LOG_CRITICAL,"%s: max boundary on omega is < 0.  Exiting.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  (*space)->data[3].gaussian = uvar->gaussianpriors;
  (*space)->data[3].span = (*space)->data[3].max - (*space)->data[3].min;

  /* output boundaries to screen */
  if (uvar->gaussianpriors) LogPrintf(LOG_DEBUG,"%s : using Gaussian priors on the following %.2f sigma ranges (except nu)\n",__func__,uvar->nsig);
  else LogPrintf(LOG_DEBUG,"%s : using flat priors on the following ranges\n",__func__);
  LogPrintf(LOG_DEBUG,"%s : parameter space, %s = [%e -> %e]\n",__func__,(*space)->data[0].name,(*space)->data[0].min,(*space)->data[0].max);
  LogPrintf(LOG_DEBUG,"%s : parameter space, %s = [%e -> %e]\n",__func__,(*space)->data[1].name,(*space)->data[1].min,(*space)->data[1].max);
  LogPrintf(LOG_DEBUG,"%s : parameter space, %s = [%e -> %e]\n",__func__,(*space)->data[2].name,(*space)->data[2].min,(*space)->data[2].max);
  LogPrintf(LOG_DEBUG,"%s : parameter space, %s = [%e -> %e]\n",__func__,(*space)->data[3].name,(*space)->data[3].min,(*space)->data[3].max);
 
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}
/**
 * Compute the prior probability density functions on the search parameters
 * The priors are computed on the search grid and correctly normalised
 *
 */
int XLALComputeAmplitudeParams(REAL8Dimension **ampspace,        /**< [out] the amplitude parameter space */ 
			       Grid **ampgrid,                   /**< [out] the amplitude grid */
			       REAL8Priors **amppriors,          /**< [out] the amplitude priors */
			       REAL8 ampsigma                    /**< [in] the amplitude prior standard deviation */
			       )
{
  UINT4 j;

  /* validate input */
  if ((*ampspace) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input REAL8Dimension structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }  
  if ((*ampgrid) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input Grid structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }  
  if ((*amppriors) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input REAL8Priors structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }  
  if (ampsigma <= 0.0) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input ampsigma <= 0.0.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }  

  /* allocate memory for the amplitude space */
  if ( ((*ampspace) = XLALCalloc(1,sizeof(REAL8Dimension))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }

  /* define amplitude space */
  snprintf((*ampspace)->name,LALNameLength,"alpha");
  (*ampspace)->min = 0.0;
  (*ampspace)->max = 2.0*ampsigma;
  (*ampspace)->mid = 0.0;
  (*ampspace)->sig = ampsigma;
  (*ampspace)->span = (*ampspace)->max;
  (*ampspace)->gaussian = 1;

  /* allocate memory for the amplitude grid */
  if ( ((*ampgrid) = XLALCalloc(1,sizeof(Grid))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }

  /* define amplitude grid params */
  (*ampgrid)->delta = (*ampspace)->max/(REAL8)(AMPVECLENGTH - 1);
  (*ampgrid)->oneoverdelta = 1.0/(*ampgrid)->delta;
  (*ampgrid)->length = AMPVECLENGTH;
  (*ampgrid)->min = (*ampspace)->min;
  strncpy((*ampgrid)->name,(*ampspace)->name,LALNameLength*sizeof(CHAR));
  
  /* allocate memory for ampltude priors */
  if (((*amppriors) = (REAL8Priors *)XLALCalloc(1,sizeof(REAL8Priors))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (((*amppriors)->logpriors = XLALCreateREAL8Vector((*ampgrid)->length)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCrateREAL8Vector() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  
  /* compute prior function on the grid */
  {
    REAL8 x0 = (*ampspace)->mid;
    REAL8 sig = (*ampspace)->sig;
    REAL8 norm = (-0.5)*log(LAL_PI) + 0.5*LAL_LN2 - log(sig);
    LogPrintf(LOG_DEBUG,"%s : computing Gaussian priors for parameter %s\n",__func__,(*ampspace)->name);
    
    /* compute prior - with amplitude prior centered on zero we double the Gaussian profile */
    for (j=0;j<(*ampgrid)->length;j++) {
      REAL8 x = (*ampgrid)->min + j*(*ampgrid)->delta;
      (*amppriors)->logpriors->data[j] = norm - 0.5*pow((x-x0)/sig,2.0);
    }
    (*amppriors)->logdelta = log((*ampgrid)->delta);
    (*amppriors)->gaussian = (*ampspace)->gaussian;
  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/**
 * Compute the prior probability density functions on the search parameters
 * The priors are computed on the search grid and correctly normalised
 *
 */
int XLALComputeBinaryPriors(REAL8PriorsVector **priors,        /**< [out] the priors on each search parameter */
			    REAL8Space *space,                /**< [in] the parameter space */
			    GridParameters *gridparams        /**< [in] the grid on this parameter space */
			    )
{
  UINT4 i,j;                   /* counters */

  /* validate input */
  if ((*priors) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input REAL8PriorsVector structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }  
  if (space == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input REAL8Space structure = NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }  
  if (gridparams == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input GridParameters structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }  
  
  /* allocate memory for the priors */
  if (((*priors) = (REAL8PriorsVector*)XLALCalloc(1,sizeof(REAL8PriorsVector))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (((*priors)->data = (REAL8Priors *)XLALCalloc(gridparams->ndim,sizeof(REAL8Priors))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  for (i=0;i<gridparams->ndim;i++) {
    if (((*priors)->data[i].logpriors = XLALCreateREAL8Vector(gridparams->grid[i].length)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : XLALCrateREAL8Vector() failed with error = %d\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_ENOMEM);
    }
  }
  
  /* loop over each dimension */
  for (i=0;i<gridparams->ndim;i++) {
       
    /* if we're using Gaussian priors on this parameter */
    if (space->data[i].gaussian) {
      
      REAL8 x0 = space->data[i].mid;
      REAL8 sig = space->data[i].sig;
      LogPrintf(LOG_DEBUG,"%s : computing Gaussian priors for parameter %s\n",__func__,space->data[i].name);
      
      /* account for single template situations, i.e known parameters */ 
      if (gridparams->grid[i].length>1) {
	for (j=0;j<gridparams->grid[i].length;j++) {
	  REAL8 x = gridparams->grid[i].min + j*gridparams->grid[i].delta;
	  REAL8 norm = (-0.5)*log(LAL_TWOPI) - log(sig);
	  (*priors)->data[i].logpriors->data[j] = norm - 0.5*pow((x-x0)/sig,2.0);
	}
	(*priors)->data[i].logdelta = log(gridparams->grid[i].delta);
      }
      else {
	(*priors)->data[i].logpriors->data[0] = 0.0;
	(*priors)->data[i].logdelta = 0.0;
      }
      
    }
    /* otherwise we use flat priors */
    else {
      
      LogPrintf(LOG_DEBUG,"%s : computing Flat priors for parameter %s\n",__func__,space->data[i].name);
      
      /* set flat prior such that for a single template the prior has no effect once multiplied by deltax */
      /* account for single template situations, i.e known parameters */ 
      if (gridparams->grid[i].length>1) {
	for (j=0;j<gridparams->grid[i].length;j++) {
	  REAL8 flat = log(1.0/(gridparams->grid[i].delta*gridparams->grid[i].length));
	  (*priors)->data[i].logpriors->data[j] = flat;
	}
	(*priors)->data[i].logdelta = log(gridparams->grid[i].delta);
      }
      else {
	(*priors)->data[i].logpriors->data[0] = 0.0;
	(*priors)->data[i].logdelta = 0.0;
      }
      
    }
       
    /* record gaussian flag */
    (*priors)->data[i].gaussian = space->data[i].gaussian;

  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}

/**
 * Adds a simulated signal to the existing SFTs
 * The parameters of the injection are drawn from the binary parameter priors
 * with the exception of the amplitude.  We cannot do this exactly since the
 * signal is not additive but we approximate this by drawing the signal component
 * from a Poisson distribution and adding it to the existing noise.
 *
 */
int XLALAddBinarySignalToSFTVector(SFTVector **sftvec,           /**< [in/out] the input SFTs into which we add a signal */ 
				   ParameterSpace *pspace,       /**< [in] the parameter space */
				   REAL8 inject_amplitude,       /**< [in] the injection amplitude in cnts/s/m^2 */
				   INT4 seed                     /**< [in] the random number seed */
				   )
{
  COMPLEX8FFTPlan *plan = NULL;         /* FFT plan */
  COMPLEX8Vector *xt = NULL;            /* the downsampled timeseries */
  COMPLEX8Vector *xf = NULL;            /* the fft'd timeseries */
  REAL8 nu,a,tasc,omega,phi0;           /* randomly selected injection parameters */
  UINT4 i,j;                            /* counters */

  /* allocate memory for the injection parameters */
  if ((pspace->inj = XLALCalloc(1,sizeof(InjectionParameters))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if ((pspace->inj->temp.x = XLALCalloc(pspace->space->ndim,sizeof(REAL8))) == NULL)  {
    LogPrintf(LOG_CRITICAL,"%s: XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  pspace->inj->temp.ndim = pspace->space->ndim;

  /* randomly select values from the priors */
  {
    gsl_rng * r;
    REAL8PriorsVector *priors = pspace->priors;
    REAL8Space *space = pspace->space;

    /* initialise the random number generator */
    if (XLALInitgslrand(&r,seed)) {
      LogPrintf(LOG_CRITICAL,"%s: XLALinitgslrand() failed with error = %d\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_EFAULT);
    }
    
    if (priors->data[0].gaussian) nu = space->data[0].mid + gsl_ran_gaussian(r,space->data[0].sig);
    else nu = gsl_ran_flat(r,space->data[0].min,space->data[0].max);
    if (priors->data[0].gaussian) a = space->data[1].mid + gsl_ran_gaussian(r,space->data[1].sig);
    else a = gsl_ran_flat(r,space->data[1].min,space->data[1].max);
    if (priors->data[0].gaussian) tasc = space->data[2].mid + gsl_ran_gaussian(r,space->data[2].sig);
    else tasc = gsl_ran_flat(r,space->data[2].min,space->data[2].max);
    if (priors->data[0].gaussian) omega = space->data[3].mid + gsl_ran_gaussian(r,space->data[3].sig);
    else omega = gsl_ran_flat(r,space->data[3].min,space->data[3].max);

    LogPrintf(LOG_DEBUG,"%s: the injected signal parameters are :\n",__func__);
    LogPrintf(LOG_DEBUG,"%s: injected nu = %6.12f\n",__func__,nu);
    LogPrintf(LOG_DEBUG,"%s: injected asini = %6.12f\n",__func__,a);
    LogPrintf(LOG_DEBUG,"%s: injected tasc = %6.12f\n",__func__,tasc);
    LogPrintf(LOG_DEBUG,"%s: injected omega = %6.12e\n",__func__,omega);
     
    /* nuiseance parameters */
    phi0 = gsl_ran_flat(r,0,LAL_TWOPI);

    gsl_rng_free (r);
  }

  {
   
    /* define some parameters common to all SFTs */
    REAL8 tsft = 1.0/(*sftvec)->data[0].deltaF;
    UINT4 numBins = (*sftvec)->data[0].data->length;
    REAL8 f0 = (*sftvec)->data[0].f0;
    REAL8 deltaf = (*sftvec)->data[0].deltaF;
    REAL8 fhet = f0 + deltaf*floor((numBins+1)/2);
    REAL8 deltat = tsft/numBins;
    REAL8 tref = 0.0;

    /* allocate memory for the timeseries and the frequency domain output */
    if ( (xt = XLALCreateCOMPLEX8Vector(numBins)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCreateCOMPLEX8Vector() failed with error = %d\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    if ( (xf = XLALCreateCOMPLEX8Vector(numBins)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCreateCOMPLEX8Vector() failed with error = %d\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    
    /* make the reverse plan */
    if ((plan = XLALCreateForwardCOMPLEX8FFTPlan(numBins,1)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCreateCOMPLEX8FFTPlan() failed with error = %d\n",__func__,xlalErrno);
      return XLAL_EINVAL;
    }
    LogPrintf(LOG_DEBUG,"%s : created the COMPLEX8 FFT plan for signal injection\n",__func__);
    
    /* we generate a down-sampled time series for each SFT */
    for (i=0;i<(*sftvec)->length;i++) {
      
      /* initialise input and output */
      memset(xt->data,0.0,numBins*sizeof(COMPLEX8));
      memset(xf->data,0.0,numBins*sizeof(COMPLEX8));
      
      /* generate the heterodyned timeseries */
      for (j=0;j<numBins;j++) {
	REAL8 t = j*deltat + XLALGPSGetREAL8(&((*sftvec)->data[i].epoch));
	
	REAL8 phase = phi0 + LAL_TWOPI*((t-tref)*(nu-fhet) - nu*a*sin(omega*(t-tasc))) - LAL_PI/2.0;
	xt->data[j] = crectf( 0.5*inject_amplitude*cos(phase)*deltat, 0.5*inject_amplitude*sin(phase)*deltat );
      }
      
      /* perform fft */
      if (XLALCOMPLEX8VectorFFT(xf,xt,plan)) {
	LogPrintf(LOG_CRITICAL,"%s: XLALCOMPLEX8VectorFFT() failed with error = %d\n",__func__,xlalErrno);
	return XLAL_EINVAL;
      }
      LogPrintf(LOG_DEBUG,"%s : performed %d/%d FFT for signal injection\n",__func__,i+1,(*sftvec)->length);
    
      /* now add directly to the input sft - making sure to flip negative frequencies */
      for (j=0;j<(UINT4)floor((numBins+1)/2);j++) {
	(*sftvec)->data[i].data->data[j] += xf->data[j+(UINT4)floor((numBins)/2)];
      }
      for (j=(UINT4)floor((numBins+1)/2);j<numBins;j++) {
	(*sftvec)->data[i].data->data[j] += xf->data[j-(UINT4)floor((numBins+1)/2)];
      }
      
      /* TESTING */
     /*  { */
/* 	FILE *fp = NULL; */
/* 	CHAR name[256]; */
/* 	sprintf(name,"/home/chmess/temp/injtest_%d.txt",i); */
/* 	fp = fopen(name,"w"); */
/* 	for (j=0;j<(UINT4)floor((numBins+1)/2);j++) { */
/* 	  printf("i = %d j = %d\n",j,j+(UINT4)floor((numBins)/2)); */
/* 	  fprintf(fp,"%6.12f %6.12f %6.12f\n",f0+j*deltaf,xf->data[j+(UINT4)floor((numBins)/2)].re,xf->data[j+(UINT4)floor((numBins)/2)].im); */
/* 	} */
/* 	for (j=(UINT4)floor((numBins+1)/2);j<numBins;j++) { */
/* 	  printf("i = %d j = %d\n",j,j-(UINT4)floor((numBins+1)/2)); */
/* 	  fprintf(fp,"%6.12f %6.12f %6.12f\n",f0+j*deltaf,xf->data[j-(UINT4)floor((numBins+1)/2)].re,xf->data[j-(UINT4)floor((numBins+1)/2)].im); */
/* 	} */
/* 	fclose(fp); */
/*       } */
     
    }
    
  }

  /* fill in injection parameters */
  pspace->inj->amp = inject_amplitude;
  pspace->inj->temp.x[0] = nu;
  pspace->inj->temp.x[1] = a;
  pspace->inj->temp.x[2] = tasc;
  pspace->inj->temp.x[3] = omega;

  /* free memeory */
  XLALDestroyCOMPLEX8FFTPlan(plan);
  XLALDestroyCOMPLEX8Vector(xt);
  XLALDestroyCOMPLEX8Vector(xf);
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}

/**
 * Read in SFTs to an SFTVector
 */
int XLALReadSFTs(SFTVector **sftvec,        /**< [out] the input SFT data */
		 SegmentParams **segparams, /**< [out] the segment parameters (noise, sampling time, etc..) */
		 CHAR *sftbasename,         /**< [in] the SFT file basename to read in */
		 REAL8 freq,                /**< [in] the starting frequency to read in */
		 REAL8 freqband,            /**< [in] the bandwidth to read */
		 INT4 start,                /**< [in] the min GPS time of the input data */
		 INT4 end,                  /**< [in] the max GPS time of the input data*/
		 CHAR *obsid_pattern        /**< [in] the OBS-ID pattern */
  		 )
{
  static SFTConstraints constraints;
  SFTCatalog *catalog = NULL;
  SFTCatalog *newcat = NULL;
  INT4 sft_check_result = 0;
  REAL8 freqmin,freqmax;
  LIGOTimeGPS *dummy_gpsstart = NULL;
  LIGOTimeGPS *dummy_gpsend = NULL;
  LIGOTimeGPS gpsstart, gpsend;
  UINT4 i;                                    /* counters */
  INT4 count = 0;
  CHAR apid[APIDLENGTH];
  LALStatus status = blank_status;              /* for use wih non-XLAL functions */
  
  /* validate input variables */
  if (*sftvec != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input SFTVector structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }  
  if (*segparams != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input INT4Vector structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }  
  if (sftbasename == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input SFT basename string == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }  
  if (freqband < 0 ) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, frequency band must be > 0.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }  
  if ((start > 0) && (end > 0)) {
    if (start - end >= 0) {
      LogPrintf(LOG_CRITICAL,"%s: Invalid input, the start time %d >= %d end time.\n",__func__,start,end);
      XLAL_ERROR(XLAL_EINVAL);
    }
  }

  /* get sft catalog */
  /* if the input gps times are negative i.e. not set, then we pass null pointers to LALLALSFTDataFind */
  if ((start > 0) && (obsid_pattern == NULL)) {
    XLALGPSSetREAL8(&gpsstart,(REAL8)start);
    dummy_gpsstart = &gpsstart;
  }
  if ((end > 0) && (obsid_pattern == NULL)) {
    XLALGPSSetREAL8(&gpsend,(REAL8)end);
    dummy_gpsend = &gpsend;
  }  
  constraints.startTime = dummy_gpsstart;
  constraints.endTime = dummy_gpsend;
  LAL_CALL( LALSFTdataFind( &status, &catalog, sftbasename, &constraints), &status);
  LogPrintf(LOG_DEBUG,"%s : found %d SFTs\n",__func__,catalog->length);
  
  /* define actual frequency range to read in */
  freqmin = freq;
  freqmax = freqmin + freqband;

  /* allocate memory for the temporary catalog */
  if (obsid_pattern != NULL) {
    snprintf(apid,APIDLENGTH,"%s",obsid_pattern); 
    if ( (newcat = XLALCalloc(1,sizeof(SFTCatalog))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    if ( (newcat->data = (SFTDescriptor *)XLALCalloc(catalog->length,sizeof(SFTDescriptor))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    newcat->length = 0;
   
    /* manually restrict the catalog to the correct OBS ID pattern */
    /* we also check the frequency range manually because it's a lot quicker than the LAL functions */
    for (i=0;i<catalog->length;i++) {
      
      CHAR *c = NULL;
      CHAR *s_obsid,*e_obsid;
      CHAR obsid_string[STRINGLENGTH];
      REAL8 sft_fmin = catalog->data[i].header.f0;
      REAL8 sft_fmax = sft_fmin + catalog->data[i].numBins*catalog->data[i].header.deltaF;

      if ((c = strstr(catalog->data[i].comment,"Additional comment")) == NULL) {
	LogPrintf(LOG_CRITICAL,"%s: Error, couldn't find required header information in SFT files.\n",__func__);
	XLAL_ERROR(XLAL_EINVAL);
      }
      
      /* extract the OBS ID string from the comment field */
      s_obsid = strstr(c,"OBS_ID") + 9;
      e_obsid = strstr(s_obsid,"\n");
      snprintf(obsid_string,strlen(s_obsid) - strlen(e_obsid) + 1,"%s",s_obsid);

      /* if the obsid is not consistent with the requested pattern then we ignore this SFT */
      if (strstr(obsid_string,apid) != NULL) {
	if (!((sft_fmin>freqmax) || (sft_fmax<freqmin))) {
	  memcpy(&(newcat->data[newcat->length]),&(catalog->data[i]),sizeof(SFTDescriptor));
	  newcat->length++;
	}
      }
      
    }
    LogPrintf(LOG_DEBUG,"%s : found %d SFTs matching %s pattern\n",__func__,newcat->length,obsid_pattern);
  }
  else newcat = catalog;
  
  /* check CRC sums of SFTs */
  LAL_CALL ( LALCheckSFTCatalog ( &status, &sft_check_result, newcat ), &status );
  if (sft_check_result) {
    LogPrintf(LOG_CRITICAL,"%s : LALCheckSFTCatalogSFT() validity check failed with error = %d\n", sft_check_result);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : checked the SFTs\n",__func__);

  /* allocate memory for the output pcu vector (allocate full possible length) */
  if ( ((*segparams) = XLALCalloc(1,sizeof(SegmentParams))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }

  /* allocate memory for the output pcu and dt vectors (allocate full possible length) */
  if ( ((*segparams)->npcus = XLALCreateINT4Vector(catalog->length)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCreateINT4Vector() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if ( ((*segparams)->dt = XLALCreateREAL8Vector(catalog->length)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCreateREAL8Vector() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  
  /* load the SFT-vectors */
  LAL_CALL( LALLoadSFTs ( &status, sftvec, newcat, freqmin, freqmax ), &status);
  LogPrintf(LOG_DEBUG,"%s : loaded the sfts\n",__func__);

  /* associate each SFT with the corresponding number of PCUs, sampling time and OBS ID */
  for (i=0;i<(*sftvec)->length;i++) {

    /* find the corresponding catalog extry and extract the number of PCUs, sampling time and OBS ID */
    UINT4 j = 0;
    INT4 notfound = 1;
    while ((j<newcat->length) && notfound ) {
      
      /* if the epoch matches the catalog epoch then get info from catalog comment field */
      if (XLALGPSCmp(&(newcat->data[j].header.epoch),&((*sftvec)->data[i].epoch)) == 0) {
	CHAR *c = NULL;
	CHAR *s_pcus,*s_dt,*e_dt;
	CHAR npcus_string[STRINGLENGTH];
	CHAR dt_string[STRINGLENGTH];
	
	if ((c = strstr(newcat->data[j].comment,"Additional comment")) == NULL) {
	  LogPrintf(LOG_CRITICAL,"%s: Error, couldn't find required header information in SFT files.\n",__func__);
	  XLAL_ERROR(XLAL_EINVAL);
	}
	
	/* extract the relavent strings from the comment field */
	s_pcus = strstr(c,"NPCUS") + 8;	
	snprintf(npcus_string,2,"%s",s_pcus);
	s_dt = strstr(c,"DELTAT") + 9;
	e_dt = strstr(s_dt,",");
	snprintf(dt_string,strlen(s_dt) - strlen(e_dt) + 1,"%s",s_dt);
	(*segparams)->npcus->data[count] = (INT4)atoi(npcus_string);
	(*segparams)->dt->data[count] = (REAL8)atof(dt_string);
	
	/* unset the not found flag and increment the count */
	notfound = 0;
	count++;
      }
      /* increment the count over catalog entries */
      j++;
    }

  }
 
  /* resize the segment params vectors */
  if ( ((*segparams)->npcus = XLALResizeINT4Vector((*segparams)->npcus,count)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALResizeINT4Vector() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if ( ((*segparams)->dt = XLALResizeREAL8Vector((*segparams)->dt,count)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALResizeREAL8Vector() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
 
  /* we don't need the original catalog anymore -  also free the new catalog */
  LAL_CALL( LALDestroySFTCatalog( &status, &catalog ), &status);
  if (obsid_pattern) {
    XLALFree(newcat->data);
    XLALFree(newcat);
  }
  LogPrintf(LOG_DEBUG,"%s : destroyed the catalogue(s)\n",__func__);
  
  /* check if we found any SFTs */
  if ((*sftvec)->length == 0) {
    LogPrintf(LOG_CRITICAL,"%s : No SFTs found in specified frequency range.  Exiting.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/**
 * Inverse FFT all narrowband SFTs
 * In order to apply the frequency derivitive corrections we must work in the time domain
 * so here we convert all SFTs to the complex time domain.
 *
 */
int XLALSFTVectorToCOMPLEX8TimeSeriesArray(COMPLEX8TimeSeriesArray **dstimevec,      /**< [out] the downsampled timeseries */
					   SFTVector *sftvec                        /**< [in] the input SFT vector */
					   )
{
  INT4 i;                                /* counter */
  COMPLEX8FFTPlan *plan = NULL;         /* inverse FFT plan */
  UINT4 N;                               /* the length of the SFTs */

  /* validate input arguments */
  if ((*dstimevec) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output COMPLEX8TimeSeriesArray structure != NULL.\n",__func__);
    return XLAL_EINVAL;
  }
  if (sftvec == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input SFTVector structure == NULL.\n",__func__);
    return XLAL_EINVAL;
  }

  /* we check that all input SFTs are of identical length so we make a single plan */
  N = sftvec->data[0].data->length;
  for (i=0;i<(INT4)sftvec->length;i++) {
    if (sftvec->data[0].data->length != N) {
      LogPrintf(LOG_CRITICAL,"%s: Invalid input, input SFTs have different lengths.\n",__func__);
      return XLAL_EINVAL;
    }
  }
  LogPrintf(LOG_DEBUG,"%s : checked that all SFTs have length %d.\n",__func__,N);

  /* make the reverse plan */
  if ((plan = XLALCreateReverseCOMPLEX8FFTPlan(N,1)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCreateReverseCOMPLEX8FFTPlan() failed with error = %d\n",__func__,xlalErrno);
    return XLAL_EINVAL;
  }
  LogPrintf(LOG_DEBUG,"%s : created the inverse FFT plan\n",__func__);

  /* allocate memory for output */
  if (((*dstimevec) = (COMPLEX8TimeSeriesArray*)XLALCalloc(1,sizeof(COMPLEX8TimeSeriesArray))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for a COMPLEX8TimeSeriesArray structure\n",__func__,xlalErrno);
    return XLAL_ENOMEM;
  }
  if (((*dstimevec)->data = (COMPLEX8TimeSeries**)XLALCalloc(sftvec->length,sizeof(COMPLEX8TimeSeries *))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for a vector of COMPLEX8TimeSeries pointers\n",__func__,xlalErrno);
    return XLAL_ENOMEM;
  }
  (*dstimevec)->length = sftvec->length;
  LogPrintf(LOG_DEBUG,"%s : allocated memory for the output data structure\n",__func__);

  /* loop over each SFT */
  for (i=0;i<(INT4)sftvec->length;i++) {
  
    COMPLEX8Vector temp_output;
    REAL8 Tsft = 1.0/sftvec->data[i].deltaF;
    REAL8 deltaT = Tsft/N;

    /* point to input */
    COMPLEX8Vector temp_input;
    temp_input.length = sftvec->data[i].data->length;
    temp_input.data = (COMPLEX8*)sftvec->data[i].data->data;
   
    /* allocate output memory - create a COMPLEX8TimeSeries */
    if (((*dstimevec)->data[i] = XLALCreateCOMPLEX8TimeSeries("DS",&(sftvec->data[i].epoch),sftvec->data[i].f0,deltaT,&lalDimensionlessUnit,N)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCreateCOMPLEX8TimeSeries() failed to allocate memory for inverse FFT output.\n",__func__);
      return XLAL_ENOMEM;
    }
    LogPrintf(LOG_DEBUG,"%s : allocated memory for the %d/%d inverse FFT\n",__func__,i+1,sftvec->length);

    /* point temp output to timeseries */
    temp_output.length = N;
    temp_output.data = (COMPLEX8*)(*dstimevec)->data[i]->data->data;

    /* perform inverse FFT */
    if (XLALCOMPLEX8VectorFFT(&temp_output, &temp_input, plan)) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCOMPLEX16VectorFFT() failed with error = %d\n",__func__,xlalErrno);
      return XLAL_EINVAL;
    }
    LogPrintf(LOG_DEBUG,"%s : performed %d/%d inverse FFT\n",__func__,i+1,sftvec->length);

  }
  LogPrintf(LOG_DEBUG,"%s : performed inverse FFT on all %d SFTs\n",__func__,sftvec->length);

  /* free memeory */
  XLALDestroyCOMPLEX8FFTPlan(plan);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/**
 * Compute the gridding parameters on spin derivitives for all segments
 * This is simply a wrapper for the single segment function
 *
 */
int XLALComputeFreqGridParamsVector(GridParametersVector **freqgridparams,    /**< [out] the gridding parameters */
				    REAL8Space *space,                        /**< [in] the orbital parameter space */
				    SFTVector *sftvec,                        /**< [in] the input SFTs */
				    REAL8 mu                                  /**< [in] the required mismatch */ 
				    )
{
  UINT4 i;                              /* counter */

  /* validate input arguments */
  if ((*freqgridparams) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output GridParamsVector structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (space == NULL) {
     LogPrintf(LOG_CRITICAL,"%s: Invalid input, input REAL8Space structure == NULL.\n",__func__);
     XLAL_ERROR(XLAL_EINVAL);
  }
  if (sftvec == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input SFTVector structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if ((mu < 0)||(mu>1)) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input mismatch parameter, not in range 0 -> 1.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  
  /* allocate memory for each set of grid parameters */
  if (((*freqgridparams) = (GridParametersVector*)XLALCalloc(1,sizeof(GridParametersVector))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCalloc() falied with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (((*freqgridparams)->segment = (GridParameters**)XLALCalloc(sftvec->length,sizeof(GridParameters*))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for a COMPLEX8TimeSeriesArray structure\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  (*freqgridparams)->length = sftvec->length;

  /* loop over each SFT */
  for (i=0;i<sftvec->length;i++) {
    
    REAL8 t0 = XLALGPSGetREAL8(&(sftvec->data[i].epoch));
    REAL8 tsft = 1.0/sftvec->data[i].deltaF;
    REAL8 tmid = t0 + 0.5*tsft;

    if (XLALComputeFreqGridParams(&((*freqgridparams)->segment[i]),space,tmid,tsft,mu)) {
      LogPrintf(LOG_CRITICAL,"%s: XLALComputeFreqGridParams() failed with error = %d\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_EINVAL);
    }
    LogPrintf(LOG_DEBUG,"%s : computed frequency grid for SFT %d/%d\n",__func__,i+1,sftvec->length);
    
  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/**
 * Compute the gridding parameters on spin derivitives
 * The circular orbit binary phase model is phi = 2*pi*nu*( (t-tref) - a*sin( W*(t-tasc) )
 * from which we compute the min and maximum instantaneous spin derivitives.
 *
 */
int XLALComputeFreqGridParams(GridParameters **gridparams,              /**< [out] the gridding parameters */
			      REAL8Space *space,                        /**< [in] the orbital parameter space */
			      REAL8 tmid,                               /**< [in] the segment mid point */
			      REAL8 Tseg,                               /**< [in] the segment length */
			      REAL8 mu                                  /**< [in] the required mismatch */ 
			      )
{
  UINT4 i,j,k,l;                         /* counters */
  INT4 n;                                /* counter */
 /*  REAL8 nu,asini,omega,tasc;    */          /* temporary orbital parameters */
 /*  REAL8 fn[NFREQMAX];    */                 /* temporary values of spin derivitives */
  REAL8 fnmin[NFREQMAX],fnmax[NFREQMAX]; /* min and max values of spin derivitives */
  INT4 dim[NFREQMAX];                    /* flag indicating whether a dimension has width */
  INT4 ndim = -1;                        /* the number of spin derivitive dimensions required */      
  Template fdots;                        /* template for an instance of spin parameters */
  Template bintemp;                      /* template for instance of binary parameters */              
  UINT4 ngrid = 100;                     /* the number of grid points per omega and tasc to search */

  /* validate input arguments */
  if ((*gridparams) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output GridParams structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (space == NULL) {
     LogPrintf(LOG_CRITICAL,"%s: Invalid input, input ParameterSpace structure == NULL.\n",__func__);
     XLAL_ERROR(XLAL_EINVAL);
   }
   if (tmid < 0) {
     LogPrintf(LOG_CRITICAL,"%s: Invalid input, input GPS time < 0.\n",__func__);
     XLAL_ERROR(XLAL_EINVAL);
   }
   if (Tseg < 0) {
     LogPrintf(LOG_CRITICAL,"%s: Invalid input, input Tseg parameter < 0.\n",__func__);
     XLAL_ERROR(XLAL_EINVAL);
   }
   if ((mu < 0)||(mu>1)) {
     LogPrintf(LOG_CRITICAL,"%s: Invalid input, input mismatch parameter, not in range 0 -> 1.\n",__func__);
     XLAL_ERROR(XLAL_EINVAL);
   }

   /* allocte memory */
   if (((*gridparams) = (GridParameters*)XLALCalloc(1,sizeof(GridParameters))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }

   /* allocate memory for the fdots */
   if ((fdots.x = XLALCalloc(NFREQMAX,sizeof(REAL8))) == NULL) {
     LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
     XLAL_ERROR(XLAL_ENOMEM);
   }
   if ((bintemp.x = XLALCalloc(NBINMAX,sizeof(REAL8))) == NULL) {
     LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
     XLAL_ERROR(XLAL_ENOMEM);
   }
   fdots.ndim = NFREQMAX;
   bintemp.ndim = NBINMAX;
   
   /* initialise the min and max spin derivitives */
   for (n=0;n<NFREQMAX;n++) {
     fnmin[n] = 1e38;
     fnmax[n] = -1e38;
   }

   /* loop over each parameter in turn and compute the spin derivitives at the corners of the parameter space */
   for (i=0;i<2;i++) {     /* nu */
     if (i==0) bintemp.x[0] = space->data[0].min;
     else bintemp.x[0] = space->data[0].max;
     
     for (j=0;j<2;j++) {    /* a */
       if (j==0) bintemp.x[1] = space->data[1].min;
       else bintemp.x[1] = space->data[1].max;
       
       /* tasc and omega are the problematic ones so we'll perform a fine grid search over them */
       for (k=0;k<ngrid;k++) {   /* tasc */
	 bintemp.x[2] = space->data[2].min + k*(space->data[2].max-space->data[2].min)/(ngrid-1);
	 
	 for (l=0;l<ngrid;l++) {  /* omega */
	   bintemp.x[3] = space->data[3].min + l*(space->data[3].max-space->data[3].min)/(ngrid-1);
	   
	   if (XLALComputeBinaryFreqDerivitives(&fdots,&bintemp,tmid)) {
	     LogPrintf(LOG_CRITICAL,"%s : XLALComputeBinaryFreqDerivitives() failed with error = %d\n",__func__,xlalErrno);
	     XLAL_ERROR(XLAL_EFAULT);
	   }

	 /*   /\* the instantanous frequency is therefore f0 = nu - a*nu*W*cos(W*(t-tasc) ) *\/ */
/* 	   fn[0] = nu - asini*nu*omega*cos(omega*(tmid-tasc)); */
	   
/* 	   /\* the instantanous nth frequency derivitive is therefore fn = - a * nu * W^(n+1) * cos ( W*(t-tasc) + n*pi/2 ) *\/ */
/* 	   for (n=1;n<NFREQMAX;n++) { */
/* 	     fn[n] = (-1.0)*asini*nu*pow(omega,n+1)*cos( omega*(tmid-tasc) + 0.5*n*LAL_PI); */
/* 	   } */

	   /* find min and max values */
	   for (n=0;n<NFREQMAX;n++) {
	     if (fdots.x[n] < fnmin[n]) fnmin[n] = fdots.x[n];
	     if (fdots.x[n] > fnmax[n]) fnmax[n] = fdots.x[n];
	   }

	 }

       }

     }

   }
   for (n=0;n<NFREQMAX;n++) {
     LogPrintf(LOG_DEBUG,"%s : determined f%d range as [%6.12e -> %6.12e].\n",__func__,n,fnmin[n],fnmax[n]);
   }
   LogPrintf(LOG_DEBUG,"%s : midpoint epoch for this SFT is %6.12f\n",__func__,tmid);

   /* free templates */
   XLALFree(fdots.x);
   XLALFree(bintemp.x);

   /* compute the required dimensionality of the frequency derivitive grid */
   /* we check the width of a 1-D template across each dimension span */ 
   for (n=0;n<NFREQMAX;n++) {
     REAL8 gnn = pow(LAL_PI,2.0)*pow(Tseg,2*n+2)/(pow(2.0,2*n)*(2*n+3.0));
     REAL8 deltafn = 2.0*sqrt(mu/gnn);
     REAL8 span = fnmax[n] - fnmin[n];
     dim[n] = 0;
     if (span > deltafn) dim[n] = 1;
     LogPrintf(LOG_DEBUG,"%s : single template span for %d'th derivitive = %e.\n",__func__,n,deltafn);
   }
   n = NFREQMAX-1;
   while ( (n>=0) && (ndim == -1) ) {
     if (dim[n] > 0) ndim = n+1;
     n--;
   }
   if (ndim < 0) {
      LogPrintf(LOG_CRITICAL,"%s: dimensionality of frequency space < 0.  No templates required.\n",__func__);
      return XLAL_EINVAL;
   }
   LogPrintf(LOG_DEBUG,"%s : determined dimensionality of frequency space = %d.\n",__func__,ndim);

   /* allocate memory to the output */
   if ( ((*gridparams)->grid = XLALCalloc(ndim,sizeof(Grid))) == NULL) {
     LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for gridparams->grid.\n",__func__);
     return XLAL_ENOMEM;
   }
   (*gridparams)->ndim = ndim;
   LogPrintf(LOG_DEBUG,"%s : allocated memory for the output grid parameters.\n",__func__);

   /* Compute the grid spacing, grid start and span for each spin derivitive dimension */
   for (n=0;n<ndim;n++) {

     /* compute diagonal metric element and corresponding spacing */
     REAL8 gnn = pow(LAL_PI,2.0)*pow(Tseg,2*n+2)/(pow(2.0,2*n)*(2*n+3.0));
     REAL8 deltafn = 2.0*sqrt(mu/(ndim*gnn));

     /* compute number of grid points in this dimension and enforce a grid centered on the middle of the parameter space */
     INT4 length = (INT4)ceil((fnmax[n]-fnmin[n])/deltafn) + 2*NBINS;                      /* add bins at each end for safety */
     REAL8 minfn = 0.5*(fnmin[n]+fnmax[n]) - 0.5*(length-1)*deltafn;  
     
     (*gridparams)->grid[n].delta = deltafn;
     (*gridparams)->grid[n].oneoverdelta = 1.0/deltafn;
     (*gridparams)->grid[n].length = length;
     (*gridparams)->grid[n].min = minfn;
     snprintf((*gridparams)->grid[n].name,LALNameLength,"f%d",n);
     
     LogPrintf(LOG_DEBUG,"%s : %s -> [%e - %e] (%e) %d grid points.\n",__func__,(*gridparams)->grid[n].name,(*gridparams)->grid[n].min,
	       (*gridparams)->grid[n].min+((*gridparams)->grid[n].length-1)*(*gridparams)->grid[n].delta,
	       (*gridparams)->grid[n].delta,(*gridparams)->grid[n].length);
   }
   LogPrintf(LOG_DEBUG,"%s : computed output grid parameters.\n",__func__);

   /* compute some internally required parameters for the grid */
   if ( ((*gridparams)->prod = XLALCalloc(ndim,sizeof(UINT4))) == NULL) {
     LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for Template structure.\n",__func__);
     return XLAL_ENOMEM;
   }
   (*gridparams)->ndim = ndim;
   (*gridparams)->mismatch = mu;
   (*gridparams)->max = 1;
   for (k=0;k<(*gridparams)->ndim;k++) (*gridparams)->max *= (*gridparams)->grid[k].length;
   
   (*gridparams)->prod[0] = 1;
   for (k=1;k<(*gridparams)->ndim;k++) (*gridparams)->prod[k] = (*gridparams)->prod[k-1]*(*gridparams)->grid[k-1].length;
   
   LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
   return XLAL_SUCCESS;

 } 

/**
 * Compute the demodulated power for all downsampled timeseries
 * This function is simply a wrapper for XLALComputeDemodulatedPower
 *
 */
int XLALComputeDemodulatedPowerVector(REAL4DemodulatedPowerVector **power,     /**< [out] the spin derivitive demodulated power */ 
				      COMPLEX8TimeSeriesArray *dsdata,         /**< [in] the downsampled SFT data */
				      GridParametersVector *gridparams         /**< [in/out] the spin derivitive gridding parameters */
				      )
{
  UINT4 i;

  /* validate input */
  if ((*power) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output REAL4DemodulatedPowerVector structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (dsdata == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input COMPLEX8TimeSeriesArray structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (gridparams == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input GridParametersVector structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (dsdata->length != gridparams->length) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, length of downsampled data vector and grid parameters vector not equal.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* allocate memory */
  if ( ((*power) = (REAL4DemodulatedPowerVector*)XLALCalloc(1,sizeof(REAL4DemodulatedPowerVector))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for REAL4DemodulatedPowerVector structure.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if ( ((*power)->segment = (REAL4DemodulatedPower**)XLALCalloc(dsdata->length,sizeof(REAL4DemodulatedPower*))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for REAL4DemodulatedPower structure.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  (*power)->length = dsdata->length;

  /* loop over each segment */
  for (i=0;i<dsdata->length;i++) {
    
    if (XLALComputeDemodulatedPower(&((*power)->segment[i]),dsdata->data[i],gridparams->segment[i])) {
      LogPrintf(LOG_CRITICAL,"%s: XLALComputeDemodulatedPwer() failed with error = %d\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_EINVAL);
    }
    LogPrintf(LOG_DEBUG,"%s : computed demodulated power for SFT %d/%d\n",__func__,i+1,dsdata->length);

  }
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/**
 * Compute the demodulated power for a single SFT
 * This involves taking the downsampled SFT timeseries and multiplying by the
 * timeseries spin-derivitive templates in turn.  Then we inverse FFT the result
 * and square to obtain the power at a given set of freq derivitive parameters.
 *
 */
int XLALComputeDemodulatedPower(REAL4DemodulatedPower **power,     /**< [out] the spin derivitive demodulated power */ 
				COMPLEX8TimeSeries *dsdata,        /**< [in] the downsampled SFT data */
				GridParameters *gridparams             /**< [in/out] the spin derivitive gridding parameters */
				)
{
  COMPLEX8FFTPlan *plan = NULL;           /* plan for the inverse FFT */
  UINT4 i,j,k,n;                          /* counters */
  REAL8 freqoffset = 0.0;                 /* the offset between the desired starting freq and the closest bin */
  INT4 binoffset = 0;                     /* the offset in output frequency bins from the start of the sft */
  COMPLEX8Vector *temp_input = NULL;      /* the temporary input of the inverse FFT = data*template */
  COMPLEX8Vector *temp_output = NULL;     /* the temporary output of the inverse FFT */
  REAL8Vector *fn = NULL;                 /* used to store the current frequency derivitive values (0'th element not used) */
  REAL8 norm = 1.0;                       /* the normalisation factor required to return to original scaling */
  UINT4 Nvec = 1;                         /* the number of spin derivitives on which to perform FFT */
  UINT4 cnt = 0;                          /* indexes the output vector */

  /* validate input arguments */
  if ((*power) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output REAL4DemodulatedPower structure != NULL.\n",__func__);
    return XLAL_EINVAL;
  }
  if (dsdata == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input COMPLEX8TimeSeries structure == NULL.\n",__func__);
    return XLAL_EINVAL;
  }
  if (gridparams == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input GridParameters structure == NULL.\n",__func__);
    return XLAL_EINVAL;
  }
  
  /* allocate memory for the output structure */
  if ( (*power = XLALCalloc(1,sizeof(REAL4DemodulatedPower))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for REAL4DemodulatedPower structure.\n",__func__);
    return XLAL_ENOMEM;
  }

  /* allocate memory for sequentially stored power */
  if ( ((*power)->data = XLALCreateREAL4Vector(gridparams->max)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCreateREAL4Vector() failed with error = %d\n",__func__,xlalErrno);
    return XLAL_ENOMEM;
  }
   
  /* compute the number of spin-derivitive templates (not including freq) */
  for (n=1;n<gridparams->ndim;n++) {
    Nvec *= gridparams->grid[n].length;
  }
  LogPrintf(LOG_DEBUG,"%s : computed number of spin derivitives as = %d.\n",__func__,Nvec);

  /* compute timeseries parameters given the requested frequency resolution */
  {
    REAL8 Teff = gridparams->grid[0].oneoverdelta;
    UINT4 N = floor(0.5 + Teff/dsdata->deltaT);
    Teff = (REAL8)N*dsdata->deltaT;
    gridparams->grid[0].delta = 1.0/Teff;
    gridparams->grid[0].oneoverdelta = Teff;
    norm = pow(1.0/(REAL8)dsdata->data->length,2.0);
   /*  LogPrintf(LOG_DEBUG,"%s : length of FFT input = %d\n",__func__,N); */
/*     LogPrintf(LOG_DEBUG,"%s : computed effective length for inverse FFT as %f sec.\n",__func__,Teff); */
/*     LogPrintf(LOG_DEBUG,"%s : computed modified frequency resolution as %e Hz.\n",__func__,gridparams->grid[0].delta); */

    /* allocate memory for the temporary zero-padded input data */
    if ( (temp_input = XLALCreateCOMPLEX8Vector(N)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCreateCOMPLEX8Vector() failed with error = %d.\n",__func__,xlalErrno);
      return XLAL_ENOMEM;
    }
    
    /* allocate memory for the time domain zero-padded phase correction vector */
    if ( (temp_output = XLALCreateCOMPLEX8Vector(N)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCreateCOMPLEX8Vector() failed with error = %d\n",__func__,xlalErrno);
      return XLAL_ENOMEM;
    }
    
    /* create a forward complex fft plan */
    if ((plan = XLALCreateForwardCOMPLEX8FFTPlan(N,0)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCreateForwardCOMPLEX8FFTPlan() failed with error = %d\n",__func__,xlalErrno);
      return XLAL_ENOMEM;
    }
    
    /* create a vector to store the frequency derivitive values */
    if ((fn = XLALCreateREAL8Vector(gridparams->ndim)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCreateREAL8Vector() failed with error = %d\n",__func__,xlalErrno);
      return XLAL_ENOMEM;
    }

  }

  /* compute frequency offset needed for heterodyne so that desired frequency lies at an exact bin */
  {
    REAL8 newdf = gridparams->grid[0].delta;
    REAL8 closest;
    binoffset = (INT4)floor(0.5 + (gridparams->grid[0].min - dsdata->f0)/newdf);
    closest = dsdata->f0 + newdf*binoffset;
    freqoffset = gridparams->grid[0].min - closest;

    /* LogPrintf(LOG_DEBUG,"%s : requested start frequency = %e Hz -> closest frequency = %e Hz.\n",__func__,gridparams->grid[0].min,closest); */
/*     LogPrintf(LOG_DEBUG,"%s : frequency offset from closets bin = %e Hz\n",__func__,freqoffset); */
/*     LogPrintf(LOG_DEBUG,"%s : offset from first frequency = %d bins\n",__func__,binoffset); */
  }

  /* loop over each value of the spin derivitive grid (not including the frequency dimension itself) */
  for (i=0;i<Nvec;i++) {
    
    /* define current spin derivitive values - since the vectors are stored linearly we need to be able */
    /* to access the i'th element and know what spin derivitives it is for */
    UINT4 idx = i;
    for (j=gridparams->ndim-1;j>0;j--) {
      UINT4 prod = 1;
      for (k=j;k>0;k--) prod *= gridparams->grid[k].length;
      idx = idx - (UINT4)floor(idx/prod)*prod; 
      fn->data[j] = gridparams->grid[j].min + idx*gridparams->grid[j].delta;
      LogPrintf(LOG_DEBUG,"%s : for derivitive index %d -> f%d index = %d value = %e\n",__func__,i,j,idx,fn->data[j]);
    }
  
    /* initialise the input data */
    memset(temp_input->data,0.0,temp_input->length*sizeof(COMPLEX8));

    /* apply time domain phase correction - first loop over time and then spin derivitive */
    for (j=0;j<dsdata->data->length;j++) {
      
      /* compute phase correction including heterodyne to shift frequencies to match up with grid */
      REAL8 tn = j*dsdata->deltaT - 0.5*dsdata->deltaT*dsdata->data->length;
      REAL8 arg = (-1.0)*LAL_TWOPI*freqoffset*tn;
      REAL8 xr, xi;
      UINT4 fac = 1; 

      /* loop over each spin derivitive and add to phase contribution for current time sample */
      for (k=1;k<gridparams->ndim;k++) {
	tn *= tn;
	fac *= k+1;
	arg += (-1.0)*LAL_TWOPI*fn->data[k]*tn/fac;
      }

      /* compute real and imaginary parts of phase correction timeseries */
      xr = cos(arg);
      xi = sin(arg);
    
      /* multiply data by phase correction - leave the zero-padding */
      temp_input->data[j] = crectf( crealf(dsdata->data->data[j])*xr - cimagf(dsdata->data->data[j])*xi, crealf(dsdata->data->data[j])*xi + cimagf(dsdata->data->data[j])*xr );
				     
    }  

    /* FFT and square to get power */
    if (XLALCOMPLEX8VectorFFT(temp_output,temp_input,plan)) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCOMPLEX8VectorFFT() failed with error = %d\n",__func__,xlalErrno);
      return XLAL_ENOMEM;
    }
    /* LogPrintf(LOG_DEBUG,"%s : computed the FFT\n",__func__); */

    /* check that we will not overrun the output vector */
    if ( (binoffset < 0) || (binoffset + (INT4)gridparams->grid[0].length > (INT4)temp_output->length) ) {
      LogPrintf(LOG_CRITICAL,"%s: strange, required bins from demodulated power out of range of result.\n",__func__,xlalErrno);
      return XLAL_EFAULT;
    }
    
    /* fill in the actual output making sure that the frequency bins of interest are used and the result is normalised */
    for (j=0;j<gridparams->grid[0].length;j++) {
      k = j + binoffset;
      (*power)->data->data[cnt] = norm*(crealf(temp_output->data[k])*crealf(temp_output->data[k]) + cimagf(temp_output->data[k])*cimagf(temp_output->data[k])); 
      cnt++;
    }
    
  }

  /* point results gridparams pointer to the actual gridparams structure */
  (*power)->gridparams = gridparams;
   memcpy(&((*power)->epoch),&(dsdata->epoch),sizeof(LIGOTimeGPS));

  /* free memory */
  XLALDestroyCOMPLEX8Vector(temp_input);
  XLALDestroyCOMPLEX8Vector(temp_output);
  XLALDestroyREAL8Vector(fn);
  XLALDestroyCOMPLEX8FFTPlan(plan);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/**
 * Compute the background photon flux for each SFT
 * This uses the median of the power to obtain the quantity r (the background
 * photon flux).
 *
 */
int XLALEstimateBackgroundFlux(REAL8Vector **background,     /**< [out] the background flux estimate */
			       SegmentParams *segparams,     /**< [in] the number of operational PCUs for each SFT */ 
			       SFTVector *sftvec             /**< [in] the SFTs */
 			       )
{
  LALStatus status = blank_status;        /* for use wih non-XLAL functions */
  UINT4 i,j;                               /* counters */

  /* validate input arguments */
  if ((*background) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output REAL8Vector structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (segparams == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input INT4Vector structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (sftvec == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input SFTVector structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (segparams->npcus->length != sftvec->length) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, length of sftvector != length of npcus vector.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* allocate memory for background estaimte results */
  if (((*background) = XLALCreateREAL8Vector(sftvec->length)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCreateREAL8Vector() failed with error = %d.\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }

  /* loop over each SFT */
  for (i=0;i<sftvec->length;i++) {

    COMPLEX8Sequence *sft = sftvec->data[i].data;
    REAL8 *P = NULL;
    REAL8 medianbias;                       /* the bias from computing the median */
    REAL8 median;                           /* the median of the power */
    REAL8 T = 1.0/sftvec->data[i].deltaF;

    /* allocate temporary memory */
    if ((P = XLALCalloc(sft->length,sizeof(REAL8))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for gridparams->grid.\n",__func__);
      XLAL_ERROR(XLAL_ENOMEM);
    }

    /* loop over each element in the SFT and record the power */
    for (j=0;j<sft->length;j++) P[j] = (crealf(sft->data[j])*crealf(sft->data[j]) + cimagf(sft->data[j])*cimagf(sft->data[j]));

    /* sort the data */
    gsl_sort(P,1,sft->length);
 
    /* compute median */
    median = gsl_stats_median_from_sorted_data(P,1,sft->length);
   
    /* compute the median bias */
    LAL_CALL ( LALRngMedBias( &status, &medianbias, sft->length ), &status);
   
    /* record estimate */
    (*background)->data[i] = (segparams->npcus->data[i]*PCU_AREA/T)*median/medianbias;
    LogPrintf(LOG_DEBUG,"%s : Estimated the background for SFT starting at %d as %f cnts/s/m^2.\n",__func__,sftvec->data[i].epoch.gpsSeconds,(*background)->data[i]);
    /* free the power */
    XLALFree(P);
    
  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/**
 * Compute the grid on binary parameters based on the semi-coherent metric
 * We use this grid to perform the integration of the posterior and to ultimately
 * compute the Bayes factor.
 *
 */
int XLALComputeBinaryGridParams(GridParameters **binarygridparams,  /**< [out] the binary parameter grid */
				REAL8Space *space,                  /**< [in] the signal parameter space */
				REAL8 T,                            /**< [in] the duration of the observation */
				REAL8 DT,                           /**< [in] the length of the coherent segments */
				REAL8 mu                            /**< [in] the mismatch */
				)
{
  REAL8 gnn[NBINMAX];                    /* stores the diagonal metric elements */ 
  INT4 ndim = 0;                         /* the number of actual search dimensions */      
  INT4 n,k;                              /* counters */

  /* validate input arguments */
  if ((*binarygridparams) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output GridParameters structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (space == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input ParameterSpace structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (T < 0) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input T parameter < 0.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if ((DT < 0) || (DT > T)) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input DT parameter < 0 or < T.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if ( (mu < 0) || (mu>1) ) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input mismatch parameter, not in range 0 -> 1.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
 
  /* compute the semi-coherent binary metric diagonal elements */
  {
    REAL8 numax = space->data[0].max;
    REAL8 amax = space->data[1].max;
    REAL8 omegamax = space->data[3].max;
    gnn[0] = (pow(LAL_PI,2.0)/6.0)*pow(DT,2.0);                                  /* nu */
    gnn[1] = (pow(LAL_PI,2.0)/6.0)*pow(numax*omegamax*DT,2.0);                   /* a */
    gnn[2] = (pow(LAL_PI,2.0)/6.0)*pow(numax*amax*omegamax*omegamax*DT,2.0);     /* tasc */
    gnn[3] = (pow(LAL_PI,2.0)/6.0)*pow(numax*amax*omegamax*DT*T,2.0)/12.0;       /* W */
    
    /* add eccentricity parameters at some point */
    /*  gnn->data[4] = (pow(LAL_PI,2.0)/6.0)*pow(numax*amax*omegamax*DT,2.0);       /\* kappa *\/ */
    /*  gnn->data[5] = (pow(LAL_PI,2.0)/6.0)*pow(numax*amax*omegamax*DT,2.0);       /\* eta *\/ */
  }    

   /* allocate memory to the output */
  if ( ((*binarygridparams) = (GridParameters*)XLALCalloc(1,sizeof(GridParameters))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for gridparams->grid.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if ( ((*binarygridparams)->grid = (Grid*)XLALCalloc(NBINMAX,sizeof(Grid))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for gridparams->grid.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if ( ((*binarygridparams)->prod = (UINT4*)XLALCalloc(NBINMAX,sizeof(UINT4))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for Template structure.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  (*binarygridparams)->ndim = NBINMAX;
  LogPrintf(LOG_DEBUG,"%s : allocated memory for the output grid parameters.\n",__func__);
  
  /* we need to determine the true number of searchable dimensions */
  /* we check the width of a 1-D template across each dimension span */ 
  for (n=0;n<NBINMAX;n++) {
    REAL8 deltax = 2.0*sqrt(mu/gnn[n]);
    if (space->data[n].span > deltax) ndim++;
  }
  LogPrintf(LOG_DEBUG,"%s : determined true dimensionality of binary space = %d.\n",__func__,ndim);

  /* Compute the grid spacing, grid start and span for each spin derivitive dimension */
  for (n=0;n<NBINMAX;n++) {
    
    REAL8 deltax;
    UINT4 length;
    REAL8 xmin;

    /* only if we have a non-zero span */
    if (space->data[n].span > 0) {
      
      /* compute spacing for this parameter given the total number of true search dimensions and the mismatch */
      deltax = 2.0*sqrt(mu/(ndim*gnn[n]));
      
      /* compute number of grid points in this dimension and enforce a grid centered on the middle of the parameter space */
      length = MYMAX((UINT4)ceil((space->data[n].span)/deltax),1);
      xmin = 0.5*(space->data[n].min + space->data[n].max) - 0.5*(length-1)*deltax;
      
    }
    else {
      
      /* otherwise set the space boundaries accordingly */
      deltax = 0.0;
      length = 1;
      xmin = space->data[n].min;

    }
    
    (*binarygridparams)->grid[n].delta = deltax;
    (*binarygridparams)->grid[n].oneoverdelta = 1.0/deltax;
    (*binarygridparams)->grid[n].length = length;
    (*binarygridparams)->grid[n].min = xmin;
    strncpy((*binarygridparams)->grid[n].name,space->data[n].name,LALNameLength*sizeof(CHAR));

    LogPrintf(LOG_DEBUG,"%s : %s -> [%e - %e] (%e) %d grid points.\n",__func__,(*binarygridparams)->grid[n].name,(*binarygridparams)->grid[n].min,
	      (*binarygridparams)->grid[n].min+(*binarygridparams)->grid[n].length*(*binarygridparams)->grid[n].delta,
	      (*binarygridparams)->grid[n].delta,(*binarygridparams)->grid[n].length);
  }
  LogPrintf(LOG_DEBUG,"%s : computed output grid parameters.\n",__func__);

  /* compute some internally required parameters for the grid */  
  (*binarygridparams)->mismatch = mu;
  (*binarygridparams)->max = 1;
  for (k=0;k<(INT4)(*binarygridparams)->ndim;k++) (*binarygridparams)->max *= (*binarygridparams)->grid[k].length;
  
  (*binarygridparams)->prod[0] = 1;
  for (k=1;k<(INT4)(*binarygridparams)->ndim;k++) (*binarygridparams)->prod[k] = (*binarygridparams)->prod[k-1]*(*binarygridparams)->grid[k-1].length;
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/**
 * Compute the Bayes factor for the semi-coherent search
 * This function performs the integral over the binary parameter space on
 * the posterior probability distribution (likelihood*prior).
 *
 */
int XLALComputeBayesFactor(BayesianProducts **Bayes,                /**< [out] the Bayesian analysis products */
			   REAL4DemodulatedPowerVector *power,      /**< [in] the input data in the form of power */
			   ParameterSpace *pspace,                  /**< [in] the parameter space */
			   REAL8 sigalpha                           /**< [in] the signal amplitude prior sigma */
			   )
{  
  LikelihoodParamsVector *Lparamsvec = NULL;          /* stores parameters required for the likelihood calculation */
  Template *bintemp = NULL;                           /* the binary parameter space template */
  Template fdots;                                     /* the freq derivitive template for each segment */
  UINT4 i,j;                                          /* counters */
  REAL8 logBayesfactor = -1e200;                      /* the final BayesFactor result */
  REAL8 logBayesfactor_phase = -1e200;                /* the final BayesFactor result for the fixed amplitude signal model */
  REAL8PriorsVector *priors = pspace->priors;         /* shortcut pointer to priors */
  UINT4 percent = 0;                                  /* counter for status update */
  REAL8 thisdelta = 0.0;                              /* initialise the prior spacing contribution for this binary template */
  gsl_interp_accel *bess_acc = NULL;                  /* gsl interpolation structure for bessel LUT */
  gsl_spline *bess_spline = NULL;                     /* gsl interpolation structure for bessel LUT */
  gsl_interp_accel *log_acc = NULL;                  /* gsl interpolation structure for log LUT */
  gsl_spline *log_spline = NULL;                     /* gsl interpolation structure for log LUT */
 
  /* validate input parameters */
  if ((*Bayes) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output BayesianProducts structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (power == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input REAL4DemodulatedPowerVector structure = NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (pspace == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input GridParameters structure = NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (sigalpha < 0.0) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input sigalphs must be > 0.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  
  /* setup parameters for the likelihood computation */
  if (XLALSetupLikelihood(&Lparamsvec,Bayes,power,pspace->gridparams,pspace->ampgrid,sigalpha)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALSetupLikelihood() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_EFAULT);
  }
  bess_acc = Lparamsvec->logbesselI0_acc;
  bess_spline = Lparamsvec->logbesselI0_spline;
  log_acc = Lparamsvec->log_acc;
  log_spline = Lparamsvec->log_spline;

  /* allocate memory for the fdots */
  if ((fdots.x = XLALCalloc(power->segment[0]->gridparams->ndim,sizeof(REAL8))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  fdots.ndim = power->segment[0]->gridparams->ndim;

  /* compute binary priors grid spacing contribution for all binary templates */
  for (j=0;j<pspace->gridparams->ndim;j++) {
    thisdelta += priors->data[j].logdelta;
  }

  /* single loop over binary templates */
  while (XLALGetNextBinaryTemplate(&bintemp,pspace->gridparams)) {

    REAL8 thisprior = 0.0;                          /* initialise the prior for this binary template */
    REAL8 logLratiosum = 0.0;                       /* initialise likelihood ratio */
    REAL8 logLratiosum_phase = -1e200;              /* initialise likelihood ratio for the fixed amplitude signal model */

    /* initialise the temporary logLratio vector */
    memset(Lparamsvec->logLratio_phase->data,0.0,Lparamsvec->logLratio_phase->length*sizeof(REAL8));

    /* compute binary priors contribution for current n-dim template */
    for (j=0;j<pspace->gridparams->ndim;j++) {
      thisprior += priors->data[j].logpriors->data[bintemp->idx[j]];
    }
    
    /** loop over segments **********************************************************************************/
    for (i=0;i<power->length;i++) {
      
      REAL4DemodulatedPower *currentpower = power->segment[i];
      GridParameters *fdotgrid = power->segment[i]->gridparams;
      REAL8 tmid = XLALGPSGetREAL8(&(power->segment[i]->epoch)) + 0.5*pspace->tseg;
      LikelihoodParams Lparams = Lparamsvec->data[i];
      /* UINT4 idx = 0; */
      INT4 idx = 0;
      REAL8 logLratio = 0.0;
      
      /* compute instantaneous frequency derivitives corresponding to the current template for this segment */
      XLALComputeBinaryFreqDerivitives(&fdots,bintemp,tmid);

      /* find indices corresponding to the spin derivitive values for the segment power */
      for (j=0;j<fdots.ndim;j++) {
	UINT4 tempidx = 0.5 + (fdots.x[j] - fdotgrid->grid[j].min)*fdotgrid->grid[j].oneoverdelta;
	idx += tempidx*fdotgrid->prod[j];
      }

      /** DEBUGGING - remove me **/
      /* if ((idx>=(INT4)currentpower->data->length)||(idx<0)) { 
        printf("segment number %d start time = %d %d\n",i,power->segment[i]->epoch.gpsSeconds,power->segment[i]->epoch.gpsNanoSeconds);
        printf("idx = %d\n",idx);
	printf("binary template params (nu = %6.12f asini = %6.12f tasc = %6.12f W = %6.12f)\n",bintemp->x[0],bintemp->x[1],bintemp->x[2],bintemp->x[3]);
        printf("freq derivitives (f0 = %6.12f)\n",fdots.x[0]);
        printf("data boundaries (fstart = %6.12f fend = %6.12f delta = %6.12f length = %d)\n",power->segment[i]->gridparams->grid[0].min,power->segment[i]->gridparams->grid[0].min+power->segment[i]->gridparams->grid[0].delta*power->segment[i]->gridparams->grid[0].length,power->segment[i]->gridparams->grid[0].delta,power->segment[i]->gridparams->grid[0].length);
      } */
      
      /* define the power at this location in this segment */
      Lparamsvec->power->data[i] = currentpower->data->data[idx];

      /* compute the likelihood for this location given the power value */
      /* inside loop over segments we compute the product of likelihood ratios */
      /* this is the sum of log-likelihood ratios */
      /* logLratio = XLALComputePhaseAmpMargLogLRatio(X,&Lparams); */
      logLratio = XLALComputePhaseAmpMargLogLRatioLUT(Lparamsvec->power->data[i],&Lparams,bess_acc,bess_spline);
      logLratiosum += logLratio;

      /** individual SFT stuff *************************************************************/
      /* apply binary priors for the individual SFT Bayes factors - this is a multiplication of likelihoods OR a sum in log-likelihoods */
      logLratio += thisprior;
      
      /* record the log BayesFactor for each SFT */
      (*Bayes)->logBayesFactor_phaseamp_vector->data[i] = XLALLogSumExpLUT((*Bayes)->logBayesFactor_phaseamp_vector->data[i],logLratio,log_acc,log_spline);
      /*************************************************************************************/
      
    } /* end loop over segments */
    /*************************************************************************************/
    
    /* compute the Bayes factor for the fixed amplitude model by summing over the amplitudes */
    /* we specifically weight the first point by 0.5 since this is on the boundary of the parameter */
    /* space and for no signal is usually a large contribution */
    if (pspace->ampspace) {

      /* compute the log likelihood ratio at each value of alpha */
      XLALComputePhaseMargLogLRatioVectorLUT(Lparamsvec->logLratio_phase,Lparamsvec->power,Lparamsvec,bess_acc,bess_spline);

      /* HACK */
      REAL8 temp = Lparamsvec->logLratio_phase->data[0] + pspace->amppriors->logpriors->data[0] - LAL_LN2;
      logLratiosum_phase = XLALLogSumExpLUT(logLratiosum_phase,temp,log_acc,log_spline);
      
      /* integrate over amplitude */
      for (j=1;j<pspace->ampgrid->length;j++) {
	REAL8 temp2 = Lparamsvec->logLratio_phase->data[j] + pspace->amppriors->logpriors->data[j];
	logLratiosum_phase = XLALLogSumExpLUT(logLratiosum_phase,temp2,log_acc,log_spline);
      }
      logLratiosum_phase += pspace->amppriors->logdelta;
    }

    /* apply binary priors - this is a multiplication of likelihoods OR a sum in log-likelihoods */
    logLratiosum += thisprior;
    if (pspace->ampspace) logLratiosum_phase += thisprior;
 
    /* for this template we contribute to each posterior vector */
    /* we sum likelihood-ratios NOT log-likelihood-ratios */
    for (i=0;i<pspace->gridparams->ndim;i++) {
      REAL8 temp = (*Bayes)->logposteriors_phaseamp[i]->data[bintemp->idx[i]];
      (*Bayes)->logposteriors_phaseamp[i]->data[bintemp->idx[i]] = XLALLogSumExpLUT(temp,logLratiosum,log_acc,log_spline); 
    }

    /* for fixed amplitude and for this template we contribute to each posterior vector */
    /* we sum likelihood-ratios NOT log-likelihood-ratios */
    for (i=0;i<pspace->gridparams->ndim;i++) {
      REAL8 temp = (*Bayes)->logposteriors_phase[i]->data[bintemp->idx[i]];
      (*Bayes)->logposteriors_phase[i]->data[bintemp->idx[i]] = XLALLogSumExpLUT(temp,logLratiosum_phase,log_acc,log_spline); 
    }

    /* for fixed amplitude and for the amplitude parameter itself we construct the posterior */
    if (pspace->ampspace) {
      for (j=0;j<pspace->ampgrid->length;j++) {
	REAL8 temp = Lparamsvec->logLratio_phase->data[j] + pspace->amppriors->logpriors->data[j] + thisprior;
	(*Bayes)->logposterior_amp->data[j] = XLALLogSumExpLUT((*Bayes)->logposterior_amp->data[j],temp,log_acc,log_spline); 
      }
    }

    /* we also sum likelihood-ratios to compute the overall Bayes-factors */
    logBayesfactor = XLALLogSumExpLUT(logBayesfactor,logLratiosum,log_acc,log_spline);
    if (pspace->ampspace) logBayesfactor_phase = XLALLogSumExpLUT(logBayesfactor_phase,logLratiosum_phase,log_acc,log_spline);
    
    /* output status to screen */
    if ((UINT4)floor(0.5 + 100*bintemp->currentidx/pspace->gridparams->max) > percent) {
      percent = (UINT4)floor(0.5 + 100*bintemp->currentidx/pspace->gridparams->max);
      LogPrintf(LOG_DEBUG,"%s : completed %d%%\n",__func__,percent);
    }

  } /* end loop over templates */
  /*************************************************************************************/

  /* normalise the Bayesfactors with the binary grid spacing */
  logBayesfactor += thisdelta;
  for (j=0;j<power->length;j++) {
    (*Bayes)->logBayesFactor_phaseamp_vector->data[j] += thisdelta;
  }
  if (pspace->ampspace) {
    logBayesfactor_phase += thisdelta;
  }
  LogPrintf(LOG_DEBUG,"%s : computed log(B) = %e\n",__func__,logBayesfactor);
  LogPrintf(LOG_DEBUG,"%s : computed log(B) (fixed amp) = %e\n",__func__,logBayesfactor_phase);

  /* point the Bayesfactor results grid to the grid used and the result obtained */
  (*Bayes)->gridparams = pspace->gridparams;
  (*Bayes)->ampgrid = pspace->ampgrid;
  (*Bayes)->logBayesFactor_phaseamp = logBayesfactor;
  if (pspace->ampspace) (*Bayes)->logBayesFactor_phase = logBayesfactor_phase;

  /* free template memory */
  XLALFree(fdots.x);

  /* free likelihood params */
  for (i=0;i<power->length;i++) {
    XLALDestroyREAL8Vector(Lparamsvec->data[i].alphaX);
    XLALDestroyREAL8Vector(Lparamsvec->data[i].alphasqY);
  }
  XLALDestroyREAL8Vector(Lparamsvec->logLratio_phase);
  XLALDestroyREAL8Vector(Lparamsvec->power);
  XLALDestroyREAL8Vector(Lparamsvec->alphasqsumY);
  gsl_spline_free(Lparamsvec->log_spline);
  gsl_interp_accel_free(Lparamsvec->log_acc);
  gsl_spline_free(Lparamsvec->logbesselI0_spline);
  gsl_interp_accel_free(Lparamsvec->logbesselI0_acc);
  XLALFree(Lparamsvec->data);
  XLALFree(Lparamsvec);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/**
 * Compute some repeatedly used parameters in the likelihood computation
 * To make the likelihood computation efficient we compute some parameters before
 * cycling over templates.  We also compute the priors.
 *
 */
int XLALSetupLikelihood(LikelihoodParamsVector **Lparamsvec,       /**< [out] set of likelihood params for each segment */
			BayesianProducts **Bayes,                  /**< [out] the output products of the Bayesian search */
			REAL4DemodulatedPowerVector *power,        /**< [in] the data in the form of power */
			GridParameters *binarygrid,                /**< [in] the binary parameter grid */
			Grid *ampgrid,                             /**< [in] the amplitude grid (NULL if not used) */
			REAL8 sigalpha                             /**< [in] the amplitude sigma prior */
			)
{
  UINT4 i,j;                            /* counters */
  REAL8 maxpower = 0.0;                 /* initialise the maximum power in the input grid */
  REAL8 maxmodpower = 0.0;              /* initialise the maximum mod power in the input grid */ 
  REAL8 maxarg = 0.0;                   /* initialise the maximum bessel function argument */
  REAL8 sumY = 0.0;                     /* initialise the sum of Y */

  /* validate input */
  

  /* allocate memory for the likelihood params */
  if (((*Lparamsvec) = XLALCalloc(1,sizeof(LikelihoodParamsVector))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (((*Lparamsvec)->data = XLALCalloc(power->length,sizeof(LikelihoodParams))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }  
  if (ampgrid) {
    if (((*Lparamsvec)->power = XLALCreateREAL8Vector(power->length)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : XLALCrateREAL8Vector() failed with error = %d\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    if (((*Lparamsvec)->logLratio_phase = XLALCreateREAL8Vector(ampgrid->length)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : XLALCrateREAL8Vector() failed with error = %d\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_ENOMEM);
    }   
    if (((*Lparamsvec)->alphasqsumY = XLALCreateREAL8Vector(ampgrid->length)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : XLALCrateREAL8Vector() failed with error = %d\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_ENOMEM);
    }  
    for (i=0;i<power->length;i++) {
      if (((*Lparamsvec)->data[i].alphasqY = XLALCreateREAL8Vector(ampgrid->length)) == NULL) {
	LogPrintf(LOG_CRITICAL,"%s : XLALCrateREAL8Vector() failed with error = %d\n",__func__,xlalErrno);
	XLAL_ERROR(XLAL_ENOMEM);
      }
      if (((*Lparamsvec)->data[i].alphaX = XLALCreateREAL8Vector(ampgrid->length)) == NULL) {
	LogPrintf(LOG_CRITICAL,"%s : XLALCrateREAL8Vector() failed with error = %d\n",__func__,xlalErrno);
	XLAL_ERROR(XLAL_ENOMEM);
      }
    }
  }

  /* allocate memory for the Bayesian output products */
  if (((*Bayes) = XLALCalloc(1,sizeof(BayesianProducts))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (((*Bayes)->logBayesFactor_phaseamp_vector = XLALCreateREAL8Vector(power->length)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCrateREAL8Vector() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (((*Bayes)->epoch = XLALCalloc(power->length,sizeof(LIGOTimeGPS))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (((*Bayes)->logposteriors_phaseamp = XLALCalloc(binarygrid->ndim,sizeof(REAL8Vector *))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  
  for (i=0;i<binarygrid->ndim;i++) { 
    if (((*Bayes)->logposteriors_phaseamp[i] = XLALCreateREAL8Vector(binarygrid->grid[i].length)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : XLALCrateREAL8Vector() failed with error = %d\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    /* initialise results */
    for (j=0;j<binarygrid->grid[i].length;j++) (*Bayes)->logposteriors_phaseamp[i]->data[j] = -1e200;
  }
  
  /* allocate memory for the fixed amplitude results */
  if (ampgrid) {
  /*   if (((*Bayes)->logBayesFactor_phase_vector = XLALCreateREAL8Vector(power->length)) == NULL) { */
/*       LogPrintf(LOG_CRITICAL,"%s : XLALCrateREAL8Vector() failed with error = %d\n",__func__,xlalErrno); */
/*       XLAL_ERROR(XLAL_ENOMEM); */
/*     } */
    if (((*Bayes)->logposterior_amp = XLALCreateREAL8Vector(ampgrid->length)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : XLALCrateREAL8Vector() failed with error = %d\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    if (((*Bayes)->logposteriors_phase = XLALCalloc(binarygrid->ndim,sizeof(REAL8Vector *))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    for (i=0;i<binarygrid->ndim;i++) {
      if (((*Bayes)->logposteriors_phase[i] = XLALCreateREAL8Vector(binarygrid->grid[i].length)) == NULL) {
	LogPrintf(LOG_CRITICAL,"%s : XLALCrateREAL8Vector() failed with error = %d\n",__func__,xlalErrno);
	XLAL_ERROR(XLAL_ENOMEM);
      }
      /* initialise results */
      for (j=0;j<binarygrid->grid[i].length;j++) (*Bayes)->logposteriors_phase[i]->data[j] = -1e200;
    }
    for (j=0;j<ampgrid->length;j++) (*Bayes)->logposterior_amp->data[j] = -1e200;
  }
  
  (*Bayes)->gridparams = binarygrid;
  (*Bayes)->ndim = binarygrid->ndim;
  (*Bayes)->nsegments = power->length;

   /* initialise results - we add using XLALlogsumexp so we initialise to the log of a very low number */
  for (j=0;j<power->length;j++) {
    (*Bayes)->logBayesFactor_phaseamp_vector->data[j] = -1e200;
  /*   (*Bayes)->logBayesFactor_phase_vector->data[j] = -1e200; */
  }

  /**************************************************************************************************/
  /* parameters for phase and amplitude marginalisation */

  /* loop over each segment and compute parameters required for the likelihood computation */
  for (i=0;i<power->length;i++) {

    /* define commonly used individual sft quantities */
    REAL8 T = power->segment[i]->duration;
    REAL8 nV = (REAL8)power->segment[i]->npcus*PCU_AREA;
    REAL8 r = power->segment[i]->r;
   
    /* find max power value */
    for (j=0;j<power->segment[i]->data->length;j++) {
      if (power->segment[i]->data->data[j] > maxpower) maxpower = power->segment[i]->data->data[j];
    }
    maxmodpower = sqrt(maxpower);

    /* compute the parameters needed for the computation of the */
    /* phase and amplitude marginalised log-likelihood ratio */ 
    REAL8 Y = 0.25*T*nV/r;
    REAL8 X = nV/r;
    REAL8 P = 1.0/(2.0*Y*sigalpha*sigalpha + 1.0);
    (*Lparamsvec)->data[i].logsqrtP = 0.5*log(P);
    (*Lparamsvec)->data[i].PQ = P*0.25*sigalpha*sigalpha*X*X;
    sumY += Y;

    /* record max value of PQ*power for bessel LUT */
    if ((*Lparamsvec)->data[i].PQ*maxpower) maxarg = (*Lparamsvec)->data[i].PQ*maxpower;

    /* if we're dealing with an amplitude grid */
    if (ampgrid) {
      for (j=0;j<ampgrid->length;j++) {
	REAL8 alpha = ampgrid->min + ampgrid->delta*(REAL8)j;
	(*Lparamsvec)->data[i].alphaX->data[j] = alpha*X;
	(*Lparamsvec)->data[i].alphasqY->data[j] = alpha*alpha*Y;
 	/* LogPrintf(LOG_DEBUG,"%s : computed alphaX = %e alpha*alpha*Y = %e for SFT %d/%d\n",__func__,(*Lparamsvec)->data[i].alphaX->data[j],(*Lparamsvec)->data[i].alphasqY->data[j],i+1,power->length); */
	
	if ((*Lparamsvec)->data[i].alphaX->data[j]*maxmodpower) maxarg = (*Lparamsvec)->data[i].alphaX->data[j]*maxmodpower;
      }
    }
    LogPrintf(LOG_DEBUG,"%s : computed X = %e Y = %e P = %e PQ = %e for SFT %d/%d\n",__func__,X,Y,P,(*Lparamsvec)->data[i].PQ,i+1,power->length);
    
    /* record epoch */
    memcpy(&((*Bayes)->epoch[i]),&(power->segment[i]->epoch),sizeof(LIGOTimeGPS));
    
  }
  LogPrintf(LOG_DEBUG,"%s : found maximum bessel function argument = %f\n",__func__,maxarg);
  
  /* add sumY to all params structures */
  for (i=0;i<ampgrid->length;i++) {
    REAL8 alpha = ampgrid->min + ampgrid->delta*(REAL8)i;
    (*Lparamsvec)->alphasqsumY->data[i] = alpha*alpha*sumY;
  }

  /* log look-up-table */
  {
    REAL8 x[NLOGLUT];
    REAL8 y[NLOGLUT];
    (*Lparamsvec)->log_acc = gsl_interp_accel_alloc();
    (*Lparamsvec)->log_spline = gsl_spline_alloc(gsl_interp_cspline,NLOGLUT);    
    
    /* precompute discrete function */
    for (i=0;i<NLOGLUT;i++) {
      y[i] = 2.0*(REAL8)i/(REAL8)(NLOGLUT-1);
      x[i] = exp(y[i]);
    }
    gsl_spline_init((*Lparamsvec)->log_spline,x,y,NLOGLUT);  
  }    
 
  /* bessel look-up-table */
  {
    REAL8 x[NBESSELLUT];
    REAL8 y[NBESSELLUT];
    (*Lparamsvec)->logbesselI0_acc = gsl_interp_accel_alloc();
    (*Lparamsvec)->logbesselI0_spline = gsl_spline_alloc(gsl_interp_cspline,NBESSELLUT);
   
    /* precompute discrete function */
    for (i=0;i<NBESSELLUT;i++) {
      x[i] = i*(1.1*maxarg)/(NBESSELLUT-1);
      y[i] = XLALLogBesselI0(x[i]);
    }
    gsl_spline_init((*Lparamsvec)->logbesselI0_spline,x,y,NBESSELLUT);  
  }    

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}

/**
 * Compute the next binary template in the grid
 * the templates are generated sequentially from the grid parameters file.  The n-dimensional
 * virtual indices of each dimension are also generated
 *
 */
int XLALGetNextBinaryTemplate(Template **temp,                        /**< [out] the signal model template parameters */
			      GridParameters *gridparams              /**< [in] the parameter space grid params */
			      )
{  
  UINT4 idx;                             /* the index variable */ 
  INT4 j;                                /* counters */

  /* if the input template is null then we allocate memory and assume we are on the first template */
  if ((*temp) == NULL) {
    if ( ((*temp) = XLALCalloc(1,sizeof(Template))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for Template structure.\n",__func__);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    if ( ((*temp)->x = XLALCalloc(gridparams->ndim,sizeof(REAL8))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for Template structure.\n",__func__);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    if ( ((*temp)->idx = XLALCalloc(gridparams->ndim,sizeof(UINT4))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for Template structure.\n",__func__);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    (*temp)->currentidx = 0;
    (*temp)->ndim = gridparams->ndim;
  }
  else if ((*temp)->currentidx == gridparams->max) {
    
    /* free binary template memory */
    XLALFree((*temp)->x);
    XLALFree((*temp)->idx);
    XLALFree(*temp);
    
    LogPrintf(LOG_DEBUG,"%s: at last template.\n",__func__);
    return 0;
  }
    
  /* initialise index */
  idx = (*temp)->currentidx;

  /* loop over each dimension and obtain the index for that dimension (store both) */
  for (j=gridparams->ndim-1;j>=0;j--) {

    /* compute the index for the j'th dimension and compute the actual value */
    UINT4 q = idx/gridparams->prod[j];
    (*temp)->x[j] = gridparams->grid[j].min + q*gridparams->grid[j].delta;
    (*temp)->idx[j] = q;

    /* update the index variable for the next dimension */
    idx = idx - q*gridparams->prod[j];
    
  }
  
  /* update index */
  (*temp)->currentidx++;

  return 1;

}

/**
 * Compute the instantaneous frequency derivitives for a given binary template and segment
 */
int XLALComputeBinaryFreqDerivitives(Template *fdots,                        /**< [out] the frequency derivitives */
				     Template *bintemp,                      /**< [in] the binary template */
				     REAL8 tmid                              /**< [in] the midpoint time of the segment */
				     )
{  
  UINT4 n;                               /* counters */
  REAL8 nu = bintemp->x[0];              /* define nu */
  REAL8 asini = bintemp->x[1];           /* define asini */
  REAL8 tasc = bintemp->x[2];            /* define tasc */
  REAL8 omega = bintemp->x[3];           /* define omega */
  
  /* precompute repeated quantities */
  REAL8 nuasiniomega = nu*asini*omega;
  REAL8 orbphase = omega*(tmid-tasc);
  REAL8 omegan = 1;

  /* the instantanous frequency is therefore f0 = nu - a*nu*W*cos(W*(t-tasc) ) */
  fdots->x[0] = nu - nuasiniomega*cos(orbphase);

  /* the instantanous nth frequency derivitive is therefore fn = - a * nu * W^(n+1) * cos ( W*(t-tasc) + n*pi/2 ) */
  for (n=1;n<fdots->ndim;n++) {
    omegan *= omega;
    fdots->x[n] = (-1.0)*nuasiniomega*omegan*cos(orbphase + 0.5*n*LAL_PI);
  }
  
  return XLAL_SUCCESS;

}

/**
 * Compute the phase and amplitude marginalised log-likelihood for a signal in Poisson noise
 * This function computes (as efficiently as possible) the log-likelihood of obtaining a particular
 * power value given Poisson noise and marginalising over an unknown phase and amplitude.
 *
 */
REAL8 XLALComputePhaseAmpMargLogLRatio(REAL8 X,                       /**< [in] the Fourier power */ 
				       LikelihoodParams *Lparams      /**< [in] pre-computed parameters useful in the likelihood */
				       ) 
{
  REAL8 arg = Lparams->PQ*X;    /* define the argument to the bessel function */

  /* compute the log-likelihood ratio */
  return Lparams->logsqrtP + arg +  XLALLogBesselI0(arg);
  
}

/**
 * Compute the phase and amplitude marginalised log-likelihood for a signal in Poisson noise
 * This function computes (as efficiently as possible) the log-likelihood of obtaining a particular
 * power value given Poisson noise and marginalising over an unknown phase and amplitude.
 *
 */
REAL8 XLALComputePhaseAmpMargLogLRatioLUT(REAL8 X,                              /**< [in] the Fourier power */ 
					  LikelihoodParams *Lparams,            /**< [in] pre-computed parameters useful in the likelihood */
					  gsl_interp_accel *logbesselI0_acc,    /**< [in] gsl interpolation accellerator */
					  gsl_spline *logbesselI0_spline        /**< [in] gsl interpolation structure */
					  ) 
{
  REAL8 arg = Lparams->PQ*X;    /* define the argument to the bessel function */

  /* compute the log-likelihood ratio */
  return Lparams->logsqrtP + arg + gsl_spline_eval(logbesselI0_spline,arg,logbesselI0_acc);
  
}

/**
 * Compute the phase marginalised log-likelihood for a signal in Poisson noise assuming constant amplitude
 * This function computes (as efficiently as possible) the log-likelihood of obtaining a particular
 * power value given Poisson noise and marginalising over an unknown phase only for a vector of amplitudes.
 *
 */
int XLALComputePhaseMargLogLRatio(REAL8Vector *logLratio_phase,  /**< [out] the output log-likelihood ratio vector (result added to input) */
				  REAL8 X,                       /**< [in] the Fourier power */ 
				  LikelihoodParams *Lparams      /**< [in] pre-computed parameters useful in the likelihood */
				  ) 
{
  REAL8 modX = sqrt(X);         /* define part of the argument to the bessel function */
  UINT4 i;                      /* counter */

  /* compute the log-likelihood ratio */
  for (i=0;i<Lparams->alphasqY->length;i++) {
    logLratio_phase->data[i] += (-1.0)*Lparams->alphasqY->data[i] + XLALLogBesselI0(Lparams->alphaX->data[i]*modX);
  }

  return XLAL_SUCCESS;

}

/**
 * Compute the phase marginalised log-likelihood for a signal in Poisson noise assuming constant amplitude
 * This function computes (as efficiently as possible) the log-likelihood of obtaining a particular
 * power value given Poisson noise and marginalising over an unknown phase only for a vector of amplitudes.
 *
 */
int XLALComputePhaseMargLogLRatioLUT(REAL8Vector *logLratio_phase,         /**< [out] the output log-likelihood ratio vector (result added to input) */
				     REAL8 X,                              /**< [in] the Fourier power */ 
				     LikelihoodParams *Lparams,            /**< [in] pre-computed parameters useful in the likelihood */
				     gsl_interp_accel *logbesselI0_acc,    /**< [in] gsl interpolation accellerator */
				     gsl_spline *logbesselI0_spline        /**< [in] gsl interpolation structure */
				     ) 
{
  REAL8 modX = sqrt(X);         /* define part of the argument to the bessel function */
  UINT4 i;                      /* counter */

  /* compute the log-likelihood ratio */
  for (i=0;i<Lparams->alphasqY->length;i++) {
    logLratio_phase->data[i] += (-1.0)*Lparams->alphasqY->data[i] + gsl_spline_eval(logbesselI0_spline,Lparams->alphaX->data[i]*modX,logbesselI0_acc);
  }
  
  return XLAL_SUCCESS;

}

/**
 * Compute the phase marginalised log-likelihood for a signal in Poisson noise assuming constant amplitude
 * This function computes (as efficiently as possible) the log-likelihood of obtaining a particular set of
 * power values given Poisson noise and marginalising over an unknown phase only for a vector of amplitudes.
 *
 */
int XLALComputePhaseMargLogLRatioVectorLUT(REAL8Vector *logLratio_phase,         /**< [out] the output log-likelihood ratio vector (result added to input) */
					   REAL8Vector *power,                   /**< [in] the Fourier power for each segment */ 
					   LikelihoodParamsVector *Lparamsvec,   /**< [in] pre-computed parameters useful in the likelihood */
					   gsl_interp_accel *logbesselI0_acc,    /**< [in] gsl interpolation accellerator */
					   gsl_spline *logbesselI0_spline        /**< [in] gsl interpolation structure */
					   ) 
{

  UINT4 i,j;                    /* counters */
 
  /* compute the log-likelihood ratio - first loop over each power value */
  for (i=0;i<power->length;i++) {
    
    LikelihoodParams *Lparams = &(Lparamsvec->data[i]);

    /* loop over each amplitude - this way of ordering the loops takes advantage of the interpolation accelerator */
    for (j=0;j<Lparamsvec->alphasqsumY->length;j++) {
    
      REAL8 arg = Lparams->alphaX->data[j]*sqrt(power->data[i]);
      logLratio_phase->data[j] += gsl_spline_eval(logbesselI0_spline,arg,logbesselI0_acc);
    }
    
  }

  /* add component independent of the data */
  for (j=0;j<Lparamsvec->alphasqsumY->length;j++) logLratio_phase->data[j] += (-1.0)*Lparamsvec->alphasqsumY->data[j];

  return XLAL_SUCCESS;

}

/**
 * Output the results to file
 * We choose to output all results from a specific analysis to a single file
 *
 */
int XLALOutputBayesResults(CHAR *outputdir,            /**< [in] the output directory name */
			   BayesianProducts *Bayes,    /**< [in] the results structure */
			   ParameterSpace *pspace,     /**< [in] the parameter space */ 
			   CHAR *clargs,               /**< [in] the command line args */
			   CHAR *obsid_pattern         /**< [in] the obsid string */
			   )
{
  CHAR outputfile[LONGSTRINGLENGTH];    /* the output filename */
  time_t curtime = time(NULL);          /* get the current time */
  CHAR *time_string = NULL;             /* stores the current time */
  CHAR *version_string = NULL;           /* pointer to a string containing the git version information */
  FILE *fp = NULL;                      /* pointer to the output file */
  UINT4 i,j;                            /* counters */

  /* validate input */
  if (outputdir == NULL) { 
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output directory string == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (Bayes == NULL) { 
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, results BayesProducts structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (pspace == NULL) { 
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, ParameterSpace structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  
  /* define the output filename */
  /* the format we adopt is the following BayesianResults-<SOURCE>-<START>_<END>-<MIN_FREQ_INT>_<MIN_FREQ_mHZ>_ <MAX_FREQ_INT>_<MAX_FREQ_mHZ>.txt */
  {
    UINT4 min_freq_int = floor(pspace->space->data[0].min);
    UINT4 max_freq_int = floor(pspace->space->data[0].max);
    UINT4 min_freq_mhz = (UINT4)floor(0.5 + (pspace->space->data[0].min - (REAL8)min_freq_int)*1e3);
    UINT4 max_freq_mhz = (UINT4)floor(0.5 + (pspace->space->data[0].max - (REAL8)max_freq_int)*1e3);
    UINT4 end = (UINT4)ceil(XLALGPSGetREAL8(&(pspace->epoch)) + pspace->span);
    if (obsid_pattern == NULL) snprintf(outputfile,LONGSTRINGLENGTH,"%s/BayesianResults-%s-%d_%d-%04d_%03d_%04d_%03d.txt",
					outputdir,pspace->source,pspace->epoch.gpsSeconds,end,min_freq_int,min_freq_mhz,max_freq_int,max_freq_mhz); 
    else snprintf(outputfile,LONGSTRINGLENGTH,"%s/BayesianResults-%s-%s-%04d_%03d_%04d_%03d.txt",
		  outputdir,pspace->source,obsid_pattern,min_freq_int,min_freq_mhz,max_freq_int,max_freq_mhz);
  }
  LogPrintf(LOG_DEBUG,"%s : output %s\n",__func__,outputfile);

  /* open the output file */
  if ((fp = fopen(outputfile,"w")) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Error, failed to open file %s for writing.  Exiting.\n",__func__,outputfile);
    XLAL_ERROR(XLAL_EINVAL);
  }
  
  /* Convert time to local time representation */
  {
    struct tm *loctime = localtime(&curtime);
    CHAR *temp_time = asctime(loctime);    
    UINT4 n = strlen(temp_time);
    time_string = XLALCalloc(n,sizeof(CHAR));
    snprintf(time_string,n-1,"%s",temp_time);
  }
  
  /* get GIT version information */
  {
    CHAR *temp_version = XLALGetVersionString(0); 
    UINT4 n = strlen(temp_version);
    version_string = XLALCalloc(n,sizeof(CHAR));
    snprintf(version_string,n-1,"%s",temp_version);
    XLALFree(temp_version);
  }

  /* output header information */
  fprintf(fp,"%s \n",version_string);
  fprintf(fp,"%%%% command line args\t\t= %s\n",clargs);
  fprintf(fp,"%%%% filename\t\t\t\t= %s\n",outputfile);
  fprintf(fp,"%%%% date\t\t\t\t\t= %s\n",time_string);
  fprintf(fp,"%%%% start time (GPS sec)\t\t= %d\n",pspace->epoch.gpsSeconds);
  fprintf(fp,"%%%% observation span (sec)\t= %d\n",(UINT4)pspace->span);
  fprintf(fp,"%%%% coherent time (sec)\t\t= %d\n",(UINT4)pspace->tseg);
  fprintf(fp,"%%%% number of segments\t\t= %d\n",Bayes->nsegments);
  fprintf(fp,"%%%% number of dimensions\t= %d\n",Bayes->gridparams->ndim);
  if (pspace->ampspace) fprintf(fp,"%%%% amplitude dimension\t\t\t= 1\n");
  else fprintf(fp,"%%%% amplitude dimension\t\t\t= 0\n");
  fprintf(fp,"%%%% mismatch\t\t\t\t= %6.12f\n",Bayes->gridparams->mismatch);
  fprintf(fp,"%%%%\n");

  /* if an injection has been performed we output the injection parameters */
  if (pspace->inj) {
    fprintf(fp,"%%%% injection parameters --------------------------------------------------------------------------------\n");
    fprintf(fp,"%%%% inj_amp\t\t\t= %6.12f\n",pspace->inj->amp);
    fprintf(fp,"%%%% inj_nu\t\t\t= %6.12f\n",pspace->inj->temp.x[0]);
    fprintf(fp,"%%%% inj_asini\t\t\t= %6.12f\n",pspace->inj->temp.x[1]);
    fprintf(fp,"%%%% inj_tasc\t\t\t= %6.12f\n",pspace->inj->temp.x[2]);
    fprintf(fp,"%%%% inj_omega\t\t\t= %6.12e\n",pspace->inj->temp.x[3]);
    fprintf(fp,"%%%% -----------------------------------------------------------------------------------------------------\n");
    fprintf(fp,"%%%%\n");
  }

  /* output the main Bayes factor results */
  fprintf(fp,"%%%% log Bayes Factor (phase and amplitude marginalised per segment)\t= %6.12e\n",Bayes->logBayesFactor_phaseamp);
  fprintf(fp,"%%%% log Bayes Factor (phase marginalised per segment)\t\t\t= %6.12e\n",Bayes->logBayesFactor_phase);
  fprintf(fp,"%%%%\n");
  fprintf(fp,"%%%% GPS start\tGPS end\tlog Bayes Factor\n");
  fprintf(fp,"%%%%\n");

  /* output the Bayes factor for each segment */
  for (i=0;i<Bayes->nsegments;i++) fprintf(fp,"%d\t%d\t%6.12e\n",Bayes->epoch[i].gpsSeconds,
					   Bayes->epoch[i].gpsSeconds+(UINT4)pspace->tseg,
					   Bayes->logBayesFactor_phaseamp_vector->data[i]);
  
  fprintf(fp,"%%%%\n");
  
  /* output the amplitude posterior */
  if (pspace->ampspace) {
    fprintf(fp,"%%%% -------------------------------------------------------------------------------------------------------\n%%%%\n");
    fprintf(fp,"%%%% name_0\t= %s\n",Bayes->ampgrid->name);
    fprintf(fp,"%%%% min_0\t= %6.12e\n",pspace->ampspace->min);
    fprintf(fp,"%%%% max_0\t= %6.12e\n",pspace->ampspace->max);
    fprintf(fp,"%%%% sig_0\t= %6.12e\n",pspace->ampspace->sig);
    fprintf(fp,"%%%% start_0\t= %6.12e\n",Bayes->ampgrid->min);
    fprintf(fp,"%%%% delta_0\t= %6.12e\n",Bayes->ampgrid->delta);
    fprintf(fp,"%%%% length_0\t= %d\n",Bayes->ampgrid->length);
    if (pspace->amppriors->gaussian) fprintf(fp,"%%%% prior_0\t= GAUSSIAN\n");
    else fprintf(fp,"%%%% prior_0\t= FLAT\n"); 
    fprintf(fp,"%%%%\n%%%%\t%s\t\tlog_post(%s)\t\tnorm_post(%s)\tlog_post_fixedamp(%s)\t\tnorm_post_fixedamp(%s)\tnorm_prior(%s)\n%%%%\n",
	    Bayes->ampgrid->name,Bayes->ampgrid->name,
	    Bayes->ampgrid->name,Bayes->ampgrid->name,
	    Bayes->ampgrid->name,Bayes->ampgrid->name);
  
    /* output posteriors - we output un-normalised and normalised posteriors plus priors */
    {
      REAL8 sum = 0.0;
      REAL8 mx = Bayes->logposterior_amp->data[0];
      for (j=0;j<Bayes->ampgrid->length;j++) if (Bayes->logposterior_amp->data[j] > mx) mx = Bayes->logposterior_amp->data[j];
  
      /* compute normalising constant for the variable amplitude posteriors */
      for (j=0;j<Bayes->ampgrid->length;j++) {
	sum += exp(Bayes->logposterior_amp->data[j]-mx)*Bayes->ampgrid->delta;
      }
  
      /* output posteriors and priors to file */
      for (j=0;j<Bayes->ampgrid->length;j++) {
	REAL8 x = Bayes->ampgrid->min + j*Bayes->ampgrid->delta;
	REAL8 log_post = Bayes->logposterior_amp->data[j];
	REAL8 norm_post = exp(Bayes->logposterior_amp->data[j]-mx)/sum;
	REAL8 norm_prior = exp(pspace->amppriors->logpriors->data[j]);
	fprintf(fp,"%6.12e\t%6.12e\t%6.12e\t0.0\t0.0\t%6.12e\n",x,log_post,norm_post,norm_prior);
      }
  
    }

  }

  /* loop over each search dimension and output the grid parameters and posteriors */
  for (i=0;i<Bayes->gridparams->ndim;i++) {
    UINT4 idx = i;
    if (pspace->ampspace) idx = i+1;
    fprintf(fp,"%%%% -------------------------------------------------------------------------------------------------------\n%%%%\n");
    fprintf(fp,"%%%% name_%d\t= %s\n",idx,Bayes->gridparams->grid[i].name);
    fprintf(fp,"%%%% min_%d\t= %6.12e\n",idx,pspace->space->data[i].min);
    fprintf(fp,"%%%% max_%d\t= %6.12e\n",idx,pspace->space->data[i].max);
    fprintf(fp,"%%%% sig_%d\t= %6.12e\n",idx,pspace->space->data[i].sig);
    fprintf(fp,"%%%% start_%d\t= %6.12e\n",idx,Bayes->gridparams->grid[i].min);
    fprintf(fp,"%%%% delta_%d\t= %6.12e\n",idx,Bayes->gridparams->grid[i].delta);
    fprintf(fp,"%%%% length_%d\t= %d\n",idx,Bayes->gridparams->grid[i].length);
    if (pspace->priors->data[i].gaussian) fprintf(fp,"%%%% prior_%d\t= GAUSSIAN\n",idx);
    else fprintf(fp,"%%%% prior_%d\t= FLAT\n",idx); 
    fprintf(fp,"%%%%\n%%%%\t%s\t\tlog_post(%s)\t\tnorm_post(%s)\tlog_post_fixedamp(%s)\t\tnorm_post_fixedamp(%s)\tnorm_prior(%s)\n%%%%\n",
	    Bayes->gridparams->grid[i].name,Bayes->gridparams->grid[i].name,
	    Bayes->gridparams->grid[i].name,Bayes->gridparams->grid[i].name,
	    Bayes->gridparams->grid[i].name,Bayes->gridparams->grid[i].name);

    /* output posteriors - we output un-normalised and normalised posteriors plus priors */
    {
      REAL8 sum = 0.0;
      REAL8 sum_phase = 0.0;
      REAL8 mx = Bayes->logposteriors_phaseamp[i]->data[0];
      REAL8 mx_phase = 0.0;
      for (j=0;j<Bayes->gridparams->grid[i].length;j++) if (Bayes->logposteriors_phaseamp[i]->data[j] > mx) mx = Bayes->logposteriors_phaseamp[i]->data[j];
      
      /* compute normalising constant for the variable amplitude posteriors */
      for (j=0;j<Bayes->gridparams->grid[i].length;j++) {
	sum += exp(Bayes->logposteriors_phaseamp[i]->data[j]-mx)*Bayes->gridparams->grid[i].delta;
      }

      /* if we're also doing the fixed amplitude case */
      if (pspace->ampspace) {

	mx_phase = Bayes->logposteriors_phase[i]->data[0];
	for (j=0;j<Bayes->gridparams->grid[i].length;j++) if (Bayes->logposteriors_phase[i]->data[j] > mx_phase) mx_phase = Bayes->logposteriors_phase[i]->data[j];
	
	/* compute normalising constant for the fixed amplitude posteriors */
	for (j=0;j<Bayes->gridparams->grid[i].length;j++) {
	  sum_phase += exp(Bayes->logposteriors_phase[i]->data[j]-mx_phase)*Bayes->gridparams->grid[i].delta;
	}
      }
      
      /* output posteriors and priors to file */
      for (j=0;j<Bayes->gridparams->grid[i].length;j++) {
	REAL8 x = Bayes->gridparams->grid[i].min + j*Bayes->gridparams->grid[i].delta;
	REAL8 log_post = Bayes->logposteriors_phaseamp[i]->data[j];
	REAL8 norm_post = exp(Bayes->logposteriors_phaseamp[i]->data[j]-mx)/sum;
	REAL8 norm_prior = exp(pspace->priors->data[i].logpriors->data[j]);
	REAL8 log_post_phase = 0.0;
	REAL8 norm_post_phase = 0.0;
	if (pspace->ampspace) {
	  log_post_phase = Bayes->logposteriors_phase[i]->data[j];
	  norm_post_phase = exp(Bayes->logposteriors_phase[i]->data[j]-mx_phase)/sum_phase;
	}
	fprintf(fp,"%6.12e\t%6.12e\t%6.12e\t%6.12e\t%6.12e\t%6.12e\n",x,log_post,norm_post,log_post_phase,norm_post_phase,norm_prior);
      }
    
    }
  
  }

  /* close the file */
  fclose(fp);

  /* free memory */
  XLALFree(time_string);
  XLALFree(version_string);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/**
 * function to compute the log of the bessel function without computing the bessel function directly
 * We compute the log of the Bessel function I0(z) bypassing the computation of the
 * function and then taking the log.  This avoids numerical problems.  The expansion
 * used is taken from Abramowitz and Stegun P.378
 *
 */
REAL8 XLALLogBesselI0(REAL8 z            /**< [in] the argument of the Bessel function */
		      )
{
  UINT4 i;                /* counter */
  REAL8 tn = 1.0;         /* initialise the variable storing t^n */

  /* for large input args */
  if (z>3.75) {

     REAL8 invt = 3.75/z;
     REAL8 y = BESSCO_HIGH[0];
    
    /* compute expansion */
    for (i=1;i<LEN_BESSCO_HIGH;i++) {
      tn = tn*invt;
      y += BESSCO_HIGH[i]*tn;
    }
    
    /* compute log of bessel function */
    return gsl_sf_log(y) + z - 0.5*gsl_sf_log(z);
    
  }
  /* for small input args */
  else {
    
    REAL8 I0 = BESSCO_LOW[0];
    REAL8 tnsq = z*z/14.0625; 
    
    for (i=1;i<LEN_BESSCO_LOW;i++) {
      tn *= tnsq;
      I0 += BESSCO_LOW[i]*tn;
    }

    return gsl_sf_log(I0); 
    
  }
  
}

/* function to compute the log of the sum of the arguments of two logged quantities
 *
 * Eg. input log(x) and log(y) -> output log(x+y)
 *
 * If you do this by exponentiating first, then summing and then logging again you
 * can easily gat overflow errors.  We use a trick to avoid this.
 */
REAL8 XLALLogSumExp(REAL8 logx,      /**< [in] the log of x */  
		    REAL8 logy       /**< [in] the log of y */
		    )
{

  /* compute log(x + y) = logmax + log(1.0 + exp(logmin - logmax)) */
  /* this way the argument to the exponential is always negative */
  /* return logmax + log(1.0 + exp(logmin - logmax)); */
  if (logy>=logx) {
    return logy + log(1.0 + exp(logx - logy));
  }
  else { 
    return logx + log(1.0 + exp(logy - logx));
  }
  
}

/* function to compute the log of the sum of the arguments of two logged quantities
 *
 * Eg. input log(x) and log(y) -> output log(x+y)
 *
 * If you do this by exponentiating first, then summing and then logging again you
 * can easily gat overflow errors.  We use a trick to avoid this.
 */
REAL8 XLALLogSumExpLUT(REAL8 logx,                   /**< [in] the log of x */  
		       REAL8 logy,                   /**< [in] the log of y */
		       gsl_interp_accel *log_acc,    /**< [in] gsl interpolation accellerator */
		       gsl_spline *log_spline        /**< [in] gsl interpolation structure */
		       )
{
  
  /* compute log(x + y) = logmax + log(1.0 + exp(logmin - logmax)) */
  /* this way the argument to the exponential is always negative */
  /* return logmax + log(1.0 + exp(logmin - logmax)); */
  if (logy>=logx) {
    REAL8 arg = 1.0 + exp(logx - logy);
    return logy + gsl_spline_eval(log_spline,arg,log_acc);
  }
  else { 
    REAL8 arg = 1.0 + exp(logy - logx);
    return logx + gsl_spline_eval(log_spline,arg,log_acc);
  }
  
}

/**
 * Free the memory allocated within a ParameterSpace structure
 */
int XLALFreeParameterSpace(ParameterSpace *pspace            /**< [in] the parameter space to be freed */
			   )
{
  UINT4 i;                     /* counter */

  /* free parameter space */
  XLALFree(pspace->space->data);
  XLALFree(pspace->space);
  
  /* free prior params */
  for (i=0;i<pspace->gridparams->ndim;i++) {
    XLALDestroyREAL8Vector(pspace->priors->data[i].logpriors);
  }
  XLALFree(pspace->priors->data);
  XLALFree(pspace->priors);
  LogPrintf(LOG_DEBUG,"%s : freed the prior parameters\n",__func__);
  
  /* free binary grid params */
  XLALFree(pspace->gridparams->grid);
  XLALFree(pspace->gridparams->prod);
  XLALFree(pspace->gridparams);
  LogPrintf(LOG_DEBUG,"%s : freed the binary grid parameters\n",__func__);

  /* free the injection parameters if used */
  if (pspace->inj) {
    XLALFree(pspace->inj->temp.x);
    XLALFree(pspace->inj);
  }

  /* free amplitude parameters */
  if (pspace->ampspace) XLALFree(pspace->ampspace);
  if (pspace->ampgrid) XLALFree(pspace->ampgrid);
  if (pspace->amppriors) {
    XLALDestroyREAL8Vector(pspace->amppriors->logpriors);
    XLALFree(pspace->amppriors);
  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/**
 * Free the memory allocated within a REAL4DemodulatedPowerVector structure
 */
int XLALFreeREAL4DemodulatedPowerVector(REAL4DemodulatedPowerVector *power            /**< [in] the data to be freed */
					)
{
  UINT4 i;                     /* counter */

  /* free each segment */
  for (i=0;i<power->length;i++) {  

    /* if there is a non-null gridparams structure then free it aswell */
    if (power->segment[i]->gridparams) {
      XLALFree(power->segment[i]->gridparams->grid);
      XLALFree(power->segment[i]->gridparams->prod);
      XLALFree(power->segment[i]->gridparams);
    }
    XLALDestroyREAL4Vector(power->segment[i]->data);
    XLALFree(power->segment[i]);
  }
  XLALFree(power->segment);
  XLALFree(power);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/**
 * Free the memory allocated within a BayesianProducts structure
 */
int XLALFreeBayesianProducts(BayesianProducts *Bayes            /**< [in] the data to be freed */
			     )
{
  UINT4 i;                     /* counter */

  /* free results */
  XLALDestroyREAL8Vector(Bayes->logBayesFactor_phaseamp_vector); 
  for (i=0;i<Bayes->ndim;i++) {
    XLALDestroyREAL8Vector(Bayes->logposteriors_phaseamp[i]);
  }
  XLALFree(Bayes->logposteriors_phaseamp);
  XLALFree(Bayes->epoch);
 
  /* if using a fixed amplitude */
 /*  if (Bayes->logBayesFactor_phase_vector) XLALDestroyREAL8Vector(Bayes->logBayesFactor_phase_vector); */
  if (Bayes->logposteriors_phase) {
    for (i=0;i<Bayes->ndim;i++) {
      XLALDestroyREAL8Vector(Bayes->logposteriors_phase[i]);
    }
  }
  if (Bayes->logposteriors_phase) XLALFree(Bayes->logposteriors_phase);
  if (Bayes->logposterior_amp) XLALDestroyREAL8Vector(Bayes->logposterior_amp);
 
  XLALFree(Bayes);
  LogPrintf(LOG_DEBUG,"%s : freed the Bayesian results\n",__func__); 

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}

/**
 * this function initialises the gsl random number generation
 * If the input seed is zero then a random seed is drawn from
 * /dev/urandom.
 *
 */
int XLALInitgslrand(gsl_rng **gslrnd,     /**< [out] the gsl random number generator */
		    INT8 seed             /**< [in] the random number generator seed */
		    )
{  
  FILE *devrandom = NULL;      /* pointer to the /dev/urandom file */
  
  /* if the seed is 0 then we draw a random seed from /dev/urandom */
  if (seed == 0) {
    
    /* open /dev/urandom */
    if ((devrandom=fopen("/dev/urandom","r")) == NULL)  {
      LogPrintf(LOG_CRITICAL,"%s: Error, unable to open device /dev/random\n",__func__);
      XLAL_ERROR(XLAL_EINVAL);
    }
    
    /* read a random seed */
    if (fread((void*)&seed,sizeof(INT8),1,devrandom) != 1) {
      LogPrintf(LOG_CRITICAL,"%s: Error, unable to read /dev/random\n",__func__);
      XLAL_ERROR(XLAL_EINVAL);
    }
    fclose(devrandom);
    
  }
  
  /* setup gsl random number generation */   
  *gslrnd = gsl_rng_alloc(gsl_rng_taus2); 
  gsl_rng_set(*gslrnd,seed);
 
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}
