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

/** \author C.Messenger
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
#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <gsl/gsl_interp.h>        /* needed for the gsl interpolation */
#include <gsl/gsl_spline.h>        /* needed for the gsl interpolation */
#include <gsl/gsl_rng.h>           /* for random number generation */
#include <gsl/gsl_randist.h>       /* for random number generation */
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf_log.h>        /* for log computation */
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
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

/* #define STRINGLENGTH 256              /\* the length of general string *\/ */
#define MINBAND 0.1                    /* the minimum amnount of bandwidth to read in */
#define SAMPFREQ 2048                    /* the orginal sampling frequency */
#define LONGSTRINGLENGTH 1024         /* the length of general string */
#define NFREQMAX 4                    /* the max dimensionality of the frequency derivitive grid */
#define NBINMAX 4                        /* the number of binary parameter dimensions */ 
#define NBINS 4                       /* the number of bins to add to each side of the fft for safety */
#define WINGS_FACTOR 2                /* the safety factor in reading extra frequency from SFTs */

/***********************************************************************************************/
/* some useful macros */

#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

/***********************************************************************************************/
/* define internal structures */

/** A single parameter dimensions boundaries
 */
typedef struct {
  REAL8 min;                        /**< the parameter space minimum */
  REAL8 max;                        /**< the parameter space maximium */
  REAL8 span;                       /**< the parameter space span */
  CHAR name[LALNameLength];         /**< string containing the name of the dimension */
} REAL8Dimension;

/** A vector of parameter space boundary information
 */
typedef struct {
  REAL8Dimension *data;             /**< the boundaries, span, etc for a single dimension */
  UINT4 ndim;                       /**< the number of dimensions */
} REAL8Space;

/** Stores the gridding parameters for a single dimension
 */
typedef struct {
  REAL8 min;                        /**< the starting points of the grid */
  REAL8 delta;                      /**< the grid spacings */
  REAL8 oneoverdelta;               /**< the inverse of the spacing */
  UINT4 length;                     /**< the number of templates in each dimension */
  CHAR name[LALNameLength];         /**< string containing the name of the dimension */
} Grid;

/** Stores the current location in a hyper-cubic parameter space
 */
typedef struct {
  REAL8 *x;                         /**< the location in parameter space */
  INT4 *idx;                        /**< the index of each dimension for this template */
  UINT4 ndim;                       /**< the dimension of the parameter space */
  UINT4 currentidx;                 /**< the current index value of the template */
} Template;

/** Stores the gridding parameters for a hypercubic grid of templates
 */
typedef struct {
  Grid *grid;                       /**< stores the parameters defining a single dimension */
  UINT4 ndim;                       /**< the number of dimensions */
  UINT4 *prod;                      /**< internal variable used to store the size of sub-dimensions */
  UINT4 max;                        /**< the maximum (total) number of templates */
  REAL8 mismatch;                   /**< the mismatch */
  INT4 Nr;
} GridParameters;

/** Stores the parameters of an injection
 */
typedef struct {
  Template temp;                    /**< stores the parameters of the signal */
  REAL8 amp;                        /**< the injected amplitude */
  REAL8 background;
} InjectionParameters;

/** contains information regarding the search parameter space
 */
typedef struct { 
  REAL8Space *space;                /**< stores the parameter space boundaries */
/*   GridParameters *gridparams;       /\**< stores the grid *\/ */
  LIGOTimeGPS epoch;                /**< the start time of the entire observation */
  REAL8 span;                       /**< the span of the entire observation */
  REAL8 tseg;                       /**< the coherent time */
  InjectionParameters *inj;         /**< stores the injected signal parameters (if any) */
} ParameterSpace;

/*****************************************************************************************/

/** Stores the gridding parameters for a hypercubic grid of templates
 */
typedef struct {
  GridParameters **segment;         /**< stores the parameters defining a single dimension */
  UINT4 length;                     /**< the number of segments */
} GridParametersVector;

/** Storage for the demodulated power from a single segment
 */
typedef struct {
  REAL4Vector *data;                /**< pointer to the power data stored sequentially */
  /* REAL8 r;                          /\**< the estimated noise background *\/ */
  LIGOTimeGPS epoch;                /**< the epoch of the segment */
/*   REAL8 duration;                   /\**< the duration of the segment *\/ */
/*   REAL8 fband;                      /\**< the bandwidth of the data *\/ */
/*   REAL8 fmin;                       /\**< the starting frequency *\/ */
/*   GridParameters *gridparams;       /\**< the grid on which the power was computed *\/ */
} REAL4DemodulatedPower;

/** Storage for the demodulated power
 */
typedef struct {
  REAL4DemodulatedPower **segment;  /**< pointer to a set of REAL4VectorArrays */
  UINT4 length;                     /**< the number of segments */
} REAL4DemodulatedPowerVector;

/** An array of COMPLEX8TimeSeries
 */
typedef struct {
  COMPLEX8TimeSeries **data;        /**< pointer to a set of COMPLEX8TimeSeries */
  UINT4 length;                     /**< the number of vectors */
 /*  REAL8 duration;                   /\**< the duration of each segment *\/ */
/*   REAL8 fband;                      /\**< the bandwidth of the data *\/ */
/*   REAL8 fmin;                       /\**< the starting frequency *\/ */
/*   UINT4 N;                          /\**< the number of complex elements per timeseries *\/ */
/*   REAL8 dt;                         /\**< the sampling time *\/ */
} COMPLEX8TimeSeriesArray;

/** A structure that stores user input variables 
 */
typedef struct { 
  BOOLEAN help;		            /**< trigger output of help string */
  CHAR *sftbasename;                /**< basename of input SFT files */
  CHAR *outputdir;                  /**< the output directory */
  REAL8 freq;                       /**< the starting frequency */
  REAL8 freqband;                   /**< the search band width */
  REAL8 minorbperiod;               /**< the minimum orbital period value */
  REAL8 maxorbperiod;               /**< the maximum orbital period value*/
  REAL8 minasini;                   /**< the minimum orbital semi-major axis value */
  REAL8 maxasini;                   /**< the maximum orbital semi-major axis value */
  REAL8 tasc;                       /**< the best guess orbital time of ascension */
  REAL8 deltaorbphase;              /**< the orbital phase uncertainty (cycles) */
  REAL8 mismatch;                   /**< the grid mismatch */      
  REAL8 frac;                         /**< the number of results to output */
  INT4 gpsstart;                    /**< the min GPS time to include */
  INT4 gpsend;                      /**< the max GPS time to include */
  INT4 seed;                        /**< fix the random number generator seed */
  REAL8 coverage;                   /**< random template bank coverage */
/*   REAL8 inject_amplitude;           /\**< the fractional amplitude of the injected signal *\/ */
  CHAR *comment;
  BOOLEAN version;	            /**< output version-info */
} UserInput_t;

/***********************************************************************************************/
/* global variables */
extern int vrbflg;	 	/**< defined in lalapps.c */
gsl_rng * r;

/***********************************************************************************************/
/* define functions */
int main(int argc,char *argv[]);
int XLALReadUserVars(int argc,char *argv[],UserInput_t *uvar, CHAR **clargs);
int XLALDefineBinaryParameterSpace(REAL8Space **space,LIGOTimeGPS epoch, REAL8 span,UserInput_t *uvar);
int XLALReadSFTs(SFTVector **sfts,CHAR *sftbasename, REAL8 freq, REAL8 freqband, INT4 gpsstart, INT4 gpsend);
int XLALComputeFreqGridParamsVector(GridParametersVector **freqgridparams,REAL8Space *pspace, SFTVector *sftvec, REAL8 mu);
int XLALComputeFreqGridParams(GridParameters **freqgridparams,REAL8Space *pspace, REAL8 tmid,REAL8 tsft, REAL8 mu);
int XLALSFTVectorToCOMPLEX8TimeSeriesArray(COMPLEX8TimeSeriesArray **dstimevec, SFTVector *sftvec);
int XLALSFTToCOMPLEX8TimeSeries(COMPLEX8TimeSeries **ts, COMPLEX8FrequencySeries *sft,COMPLEX8FFTPlan **plan);
int XLALCOMPLEX8TimeSeriesToCOMPLEX8FrequencySeries(COMPLEX8FrequencySeries **fs,const COMPLEX8TimeSeries *ts,GridParameters **gridparams);
int XLALCOMPLEX8TimeSeriesArrayToDemodPowerVector(REAL4DemodulatedPowerVector **power,COMPLEX8TimeSeriesArray *time,GridParametersVector *gridparams);
int XLALApplyPhaseCorrection(COMPLEX8TimeSeries **outts, COMPLEX8TimeSeries *ints, Template *fn);
int XLALEstimateBackgroundFlux(REAL8Vector **background, SFTVector *sftvec);
int XLALComputeBinaryGridParams(GridParameters **binarygridparams,REAL8Space *space,REAL8 T,REAL8 DT,REAL8 mu,REAL8 coverage);
int XLALComputeSemiCoherentStat(FILE *fp,REAL4DemodulatedPowerVector *power,ParameterSpace *pspace,GridParametersVector *fgrid,GridParameters *bingrid,REAL8Vector *background,REAL8 top);
int XLALGetNextTemplate(Template **temp,GridParameters *gridparams, ParameterSpace *space);
int XLALGetNextRandomBinaryTemplate(Template **temp,GridParameters *gridparams, ParameterSpace *space);				    
int XLALComputeBinaryFreqDerivitives(Template *fdots,Template *bintemp,REAL8 tmid);
int XLALFreeParameterSpace(ParameterSpace *pspace);
int XLALFreeREAL4DemodulatedPowerVector(REAL4DemodulatedPowerVector *power);
int XLALOpenSemiCoherentResultsFile(FILE **fp,CHAR *outputdir,ParameterSpace *pspace,CHAR *clargs,UserInput_t *uvar);
/* int XLALAddBinarySignalToSFTVector(SFTVector **sftvec,ParameterSpace *pspace,REAL8Vector *background,REAL8 inject_amplitude, INT4 seed); */
int XLALReplaceSFTVectornoise(SFTVector **sftvec,REAL8 background,INT4 seed);
int XLALInitgslrand(gsl_rng **gslrnd,INT8 seed);

/***********************************************************************************************/
/* empty initializers */
UserInput_t empty_UserInput;
ParameterSpace empty_ParameterSpace;

/** The main function of semicoherentbinary.c
 *
 */
int main( int argc, char *argv[] )  {

  UserInput_t uvar = empty_UserInput;           /* user input variables */
  CHAR *clargs = NULL;                          /* store the command line args */
  SFTVector *sftvec = NULL;                     /* stores the input SFTs */
  REAL8Vector *background = NULL;                  /* estimates of the background for each SFT */
  ParameterSpace pspace = empty_ParameterSpace;    /* the search parameter space */ 
  COMPLEX8TimeSeriesArray *dstimevec = NULL;       /* contains the downsampled inverse FFT'd SFTs */
  REAL4DemodulatedPowerVector *dmpower = NULL;    /* contains the demodulated power for all SFTs */
  GridParametersVector *freqgridparams = NULL;  /* the coherent grid on the frequency derivitive parameter space */
  GridParameters *bingridparams = NULL;
 /*  REAL8Vector *SemiCo = NULL;                   /\* the semi-coherent statistic results *\/    */         
  REAL8 fmin_read,fmax_read,fband_read;         /* the range of frequencies to be read from SFTs */
  UINT4 i;                                      /* counter */
  FILE *fp = NULL;

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
  {
    REAL8 wings = LAL_TWOPI*uvar.maxasini/uvar.minorbperiod;
    fmin_read = MINBAND*floor((uvar.freq - WINGS_FACTOR*uvar.freq*wings)/MINBAND);
    fmax_read = MINBAND*ceil((uvar.freq + uvar.freqband + WINGS_FACTOR*(uvar.freq + uvar.freqband)*wings)/MINBAND);
    fband_read = fmax_read - fmin_read;
    LogPrintf(LOG_DEBUG,"%s : reading in SFT frequency band [%f -> %f]\n",__func__,fmin_read,fmax_read); 
  }
 
  /* initialise the random number generator */
  if (XLALInitgslrand(&r,uvar.seed)) {
    LogPrintf(LOG_CRITICAL,"%s: XLALinitgslrand() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_EFAULT);
  }

  /**********************************************************************************/
  /* READ THE SFT DATA */
  /**********************************************************************************/
  
  /* load in the SFTs - also fill in the segment parameters structure */
  if (XLALReadSFTs(&sftvec,uvar.sftbasename,fmin_read,fband_read,uvar.gpsstart,uvar.gpsend)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadSFTs() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : read in SFTs\n",__func__);
  
  /* define SFT length and the start and span of the observations plus the definitive segment time */
  pspace.tseg = 1.0/sftvec->data[0].deltaF;
  memcpy(&(pspace.epoch),&(sftvec->data[0].epoch),sizeof(LIGOTimeGPS));
  pspace.span = XLALGPSDiff(&(sftvec->data[sftvec->length-1].epoch),&(sftvec->data[0].epoch)) + pspace.tseg;
  LogPrintf(LOG_DEBUG,"%s : SFT length = %f seconds\n",__func__,pspace.tseg);
  LogPrintf(LOG_DEBUG,"%s : entire dataset starts at GPS time %d contains %d SFTS and spans %.0f seconds\n",__func__,pspace.epoch.gpsSeconds,sftvec->length,pspace.span);

  /**********************************************************************************/
  /* ESTIMATE BACKGROUND NOISE FROM SFTS */
  /**********************************************************************************/
  
  /* compute the background noise using the sfts */
  if (XLALEstimateBackgroundFlux(&background,sftvec)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALEstimateBackgroundFlux() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : estimated the background noise from the SFTs\n",__func__); 

  /**********************************************************************************/
  /* DEFINE THE BINARY PARAMETER SPACE */
  /**********************************************************************************/
  
  /*define the binary parameter space */
  if (XLALDefineBinaryParameterSpace(&(pspace.space),pspace.epoch,pspace.span,&uvar)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALDefineBinaryParameterSpace() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : defined binary parameter prior space\n",__func__);

  /**********************************************************************************/
  /* COMPUTE THE COARSE GRID ON FREQUENCY DERIVITIVES */
  /**********************************************************************************/

  /* compute the grid parameters for all SFTs */
  if (XLALComputeFreqGridParamsVector(&freqgridparams,pspace.space,sftvec,uvar.mismatch)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALComputeFreqGridParams() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }

  /**********************************************************************************/
  /* COMPUTE THE FINE GRID PARAMETERS */
  /**********************************************************************************/
  
  /* compute the fine grid on the binary parameters */
  if (XLALComputeBinaryGridParams(&bingridparams,pspace.space,pspace.span,pspace.tseg,uvar.mismatch,uvar.coverage)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALComputeBinaryGridParams() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : computed the binary parameter space grid\n",__func__);
  
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
  /* COMPUTE THE STATISTICS ON THE COARSE GRID */
  /**********************************************************************************/

  /* compute the demodulated power on the frequency derivitive grid */
  if (XLALCOMPLEX8TimeSeriesArrayToDemodPowerVector(&dmpower,dstimevec,freqgridparams)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCOMPLEX8TimeSeriesArrayToDemodPowerVector() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : computed the demodulated power\n",__func__);

  /**********************************************************************************/
  /* OPEN RESULTS FILE */
  /**********************************************************************************/
  
  if (XLALOpenSemiCoherentResultsFile(&fp,uvar.outputdir,&pspace,clargs,&uvar)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALOutputBayesResults() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : output results to file.\n",__func__);

   /**********************************************************************************/
  /* COMPUTE THE STATISTICS ON THE FINE GRID */
  /**********************************************************************************/

  /* compute the semi-coherent detection statistic on the fine grid */
  if (XLALComputeSemiCoherentStat(fp,dmpower,&pspace,freqgridparams,bingridparams,background,uvar.frac)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALComputeSemiCoherentStat() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  fclose(fp);
  LogPrintf(LOG_DEBUG,"%s : computed the semi-coherent statistic\n",__func__);

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
  if (XLALFreeREAL4DemodulatedPowerVector(dmpower)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALFreeREAL4DemodulatedPowerVector() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : freed the demodulated power\n",__func__);

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

  /* free un-needed original SFT vector */
  XLALDestroySFTVector(sftvec);
  LogPrintf(LOG_DEBUG,"%s : Freed the SFT memory\n",__func__);

  /* free semi-coherent results */
 /*  XLALDestroyREAL8Vector(SemiCo); */
/*   LogPrintf(LOG_DEBUG,"%s : Freed the semi-coherent results memory\n",__func__); */

  /* Free config-Variables and userInput stuff */
  XLALDestroyUserVars();
  XLALFree(clargs);

  /* did we forget anything ? */
  LALCheckMemoryLeaks();
  LogPrintf(LOG_DEBUG,"%s : successfully checked memory leaks.\n",__func__);

  LogPrintf(LOG_DEBUG,"%s : successfully completed.\n",__func__);
  return 0;
  
} /* end of main */

/*******************************************************************************/
/** Read in input user arguments
 *
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
  uvar->comment = NULL;
  uvar->gpsstart = -1;
  uvar->gpsend = -1;
  uvar->mismatch = 0.2;
  uvar->frac = 0.01;
  uvar->coverage = -1;
  uvar->seed = 1;

  /* initialise all parameter space ranges to zero */
  uvar->freqband = 0;
  uvar->minorbperiod = 0.0;
  uvar->maxorbperiod = 0.0;
  uvar->minasini = 0.0;
  uvar->maxasini = 0.0;
  uvar->tasc = -1.0;
  uvar->deltaorbphase = 2.0*M_PI;

  /* ---------- register all user-variables ---------- */
  XLALregBOOLUserStruct(help, 		        'h', UVAR_HELP,     "Print this message");
  XLALregSTRINGUserStruct(sftbasename, 	        'i', UVAR_REQUIRED, "The basename of the input SFT files"); 
  XLALregSTRINGUserStruct(outputdir, 	        'o', UVAR_REQUIRED, "The output directory name"); 
  XLALregSTRINGUserStruct(comment, 	        'C', UVAR_REQUIRED, "An analysis descriptor string"); 
  XLALregREALUserStruct(freq,                   'f', UVAR_REQUIRED, "The starting frequency (Hz)");
  XLALregREALUserStruct(freqband,   	        'b', UVAR_OPTIONAL, "The frequency band (Hz)");
  XLALregREALUserStruct(minorbperiod,           'p', UVAR_REQUIRED, "The minimum orbital period value (sec)");
  XLALregREALUserStruct(maxorbperiod,   	'P', UVAR_OPTIONAL, "The maximum orbital period value (sec)");
  XLALregREALUserStruct(minasini,               'a', UVAR_REQUIRED, "The minimum orbital semi-major axis (sec)");
  XLALregREALUserStruct(maxasini,       	'A', UVAR_OPTIONAL, "The maximum orbital semi-major axis (sec)");
  XLALregREALUserStruct(tasc,                   't', UVAR_REQUIRED, "The best guess orbital time of ascension (rads)");
  XLALregREALUserStruct(deltaorbphase,      	'T', UVAR_OPTIONAL, "The orbital phase uncertainty (cycles)");
  XLALregREALUserStruct(mismatch,        	'm', UVAR_OPTIONAL, "The grid mismatch (0->1)");
  XLALregREALUserStruct(coverage,        	'c', UVAR_OPTIONAL, "The random template coverage (0->1)");
  XLALregREALUserStruct(frac,                   'x', UVAR_OPTIONAL, "output this top fraction of results");
  XLALregINTUserStruct(seed,                    'X', UVAR_OPTIONAL, "The random number seed (0 = clock)");
  XLALregINTUserStruct(gpsstart,                's', UVAR_OPTIONAL, "The minimum start time (GPS sec)");
  XLALregINTUserStruct(gpsend,          	'e', UVAR_OPTIONAL, "The maximum end time (GPS sec)");
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

/*******************************************************************************/
/** Computes the binary parameter space boundaries given the user input args
 *
 * For each search dimension we define the min, max, mid, and span of that dimension
 * plus we give each dimension a name.
 *
 */
int XLALDefineBinaryParameterSpace(REAL8Space **space,                 /**< [out] the parameter space  */
				   LIGOTimeGPS epoch,                  /**< [in] the observation start epoch */
				   REAL8 span,                         /**< [in] the observation span */
				   UserInput_t *uvar                   /**< [in] the user input variables */
				   )
{
  REAL8 midpoint;              /* the midpoint of the observation */
 /*  REAL8 newtasc;               /\* shifted value of tasc *\/ */
/*   REAL8 newdeltatasc;          /\* updated uncertainty on tasc after shifting *\/ */
  REAL8 mintasc;               /* the minimum time of ascension */
  REAL8 maxtasc;               /* the maximum time of ascension */

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
  if (uvar->tasc>0) {    
    REAL8 meanperiod = 0.5*(uvar->maxorbperiod + uvar->minorbperiod);
    INT4 n = (INT4)floor(0.5 + (midpoint - uvar->tasc)/meanperiod);
    REAL8 newtasc = uvar->tasc + n*meanperiod;
    mintasc = newtasc - 0.5*uvar->maxorbperiod*uvar->deltaorbphase;
    maxtasc = newtasc + 0.5*uvar->maxorbperiod*uvar->deltaorbphase;
    LogPrintf(LOG_DEBUG,"%s : shifted tasc by %d orbits\n",__func__,n);
  }
  else {   /* we have no orbital info so we search a full period at the midpoint */
    mintasc = midpoint - 0.5*uvar->maxorbperiod;
    maxtasc = midpoint + 0.5*uvar->maxorbperiod;
  }
  
  /* this represents a hyper-cubic parameter space */
  /* we make sure that parameter ranges are consistent i.e asini > 0 etc.. */
  /* frequency */
  snprintf((*space)->data[0].name,LALNameLength,"nu");
  (*space)->data[0].min = uvar->freq;
  (*space)->data[0].max = uvar->freq + uvar->freqband;
  (*space)->data[0].span = (*space)->data[0].max - (*space)->data[0].min;

  /* asini */
  snprintf((*space)->data[1].name,LALNameLength,"asini");
  (*space)->data[1].min = uvar->minasini;
  (*space)->data[1].max = uvar->maxasini;
  (*space)->data[1].span = (*space)->data[1].max - (*space)->data[1].min;
  
  /* orbphase */
  snprintf((*space)->data[2].name,LALNameLength,"tasc");
  (*space)->data[2].min = mintasc;
  (*space)->data[2].max = maxtasc;
  (*space)->data[2].span = (*space)->data[2].max - (*space)->data[2].min;
  
  /* omega */
  snprintf((*space)->data[3].name,LALNameLength,"omega");
  (*space)->data[3].min = LAL_TWOPI/uvar->maxorbperiod;
  (*space)->data[3].max = LAL_TWOPI/uvar->minorbperiod;
  (*space)->data[3].span = (*space)->data[3].max - (*space)->data[3].min;

  /* output boundaries to screen */
  LogPrintf(LOG_DEBUG,"%s : using flat priors on the following ranges\n",__func__);
  LogPrintf(LOG_DEBUG,"%s : parameter space, %s = [%e -> %e]\n",__func__,(*space)->data[0].name,(*space)->data[0].min,(*space)->data[0].max);
  LogPrintf(LOG_DEBUG,"%s : parameter space, %s = [%e -> %e]\n",__func__,(*space)->data[1].name,(*space)->data[1].min,(*space)->data[1].max);
  LogPrintf(LOG_DEBUG,"%s : parameter space, %s = [%e -> %e]\n",__func__,(*space)->data[2].name,(*space)->data[2].min,(*space)->data[2].max);
  LogPrintf(LOG_DEBUG,"%s : parameter space, %s = [%e -> %e]\n",__func__,(*space)->data[3].name,(*space)->data[3].min,(*space)->data[3].max);
 
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}



/*******************************************************************************/
/** Read in SFTs to an SFTVector
 *
 */
int XLALReadSFTs(SFTVector **sftvec,        /**< [out] the input SFT data */
		 CHAR *sftbasename,         /**< [in] the SFT file basename to read in */
		 REAL8 freq,                /**< [in] the starting frequency to read in */
		 REAL8 freqband,            /**< [in] the bandwidth to read */
		 INT4 start,                /**< [in] the min GPS time of the input data */
		 INT4 end                   /**< [in] the max GPS time of the input data*/
  		 )
{
  static SFTConstraints constraints;
  SFTCatalog *catalog = NULL;
  INT4 sft_check_result = 0;
  REAL8 freqmin,freqmax;
  LIGOTimeGPS *dummy_gpsstart = NULL;
  LIGOTimeGPS *dummy_gpsend = NULL;
  LIGOTimeGPS gpsstart, gpsend;
  LALStatus status = blank_status;              /* for use wih non-XLAL functions */
  
  /* validate input variables */
  if (*sftvec != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input SFTVector structure != NULL.\n",__func__);
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
  if (start > 0) {
    XLALGPSSetREAL8(&gpsstart,(REAL8)start);
    dummy_gpsstart = &gpsstart;
  }
  if (end > 0) {
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
  
  /* check CRC sums of SFTs */
  LAL_CALL ( LALCheckSFTCatalog ( &status, &sft_check_result, catalog ), &status );
  if (sft_check_result) {
    LogPrintf(LOG_CRITICAL,"%s : LALCheckSFTCatalogSFT() validity check failed with error = %d\n", sft_check_result);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : checked the SFTs\n",__func__);
  
  /* load the SFT-vectors */
  LAL_CALL( LALLoadSFTs ( &status, sftvec, catalog, freqmin, freqmax ), &status);
  LogPrintf(LOG_DEBUG,"%s : loaded the sfts\n",__func__);
  
  /* we don't need the original catalog anymore */
  LAL_CALL( LALDestroySFTCatalog( &status, &catalog ), &status);
  LogPrintf(LOG_DEBUG,"%s : destroyed the catalogue(s)\n",__func__);
  
  /* check if we found any SFTs */
  if ((*sftvec)->length == 0) {
    LogPrintf(LOG_CRITICAL,"%s : No SFTs found in specified frequency range.  Exiting.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}

/*******************************************************************************/
/** Inverse FFT all narrowband SFTs
 *
 * In order to apply the frequency derivitive corrections we must work in the time domain
 * so here we convert all SFTs to the complex time domain.
 *
 */
int XLALSFTVectorToCOMPLEX8TimeSeriesArray(COMPLEX8TimeSeriesArray **dstimevec,      /**< [out] the downsampled timeseries */
					   SFTVector *sftvec                         /**< [in] the input SFT vector */
					   )
{
  INT4 i;                                /* counter */
  COMPLEX8FFTPlan *plan = NULL;          /* inverse FFT plan */

  /* validate input arguments */
  if ((*dstimevec) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output COMPLEX8TimeSeriesArray structure != NULL.\n",__func__);
    return XLAL_EINVAL;
  }
  if (sftvec == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input SFTVector structure == NULL.\n",__func__);
    return XLAL_EINVAL;
  }

  /* allocate memory for output */
  if (((*dstimevec) = (COMPLEX8TimeSeriesArray*)XLALCalloc(1,sizeof(COMPLEX8TimeSeriesArray))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for a COMPLEX8TimeSeriesArray structure\n",__func__,xlalErrno);
    return XLAL_ENOMEM;
  }
  if (((*dstimevec)->data = (COMPLEX8TimeSeries**)XLALCalloc(sftvec->length,sizeof(COMPLEX8TimeSeries *))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for a vector of COMPLEX8TimeSeries pointers\n",__func__,xlalErrno);
    return XLAL_ENOMEM;
  }
  (*dstimevec)->length = sftvec->length;                       /* the number of timeseries */  
  LogPrintf(LOG_DEBUG,"%s : allocated memory for the output data structure\n",__func__);

  /* loop over each SFT */
  for (i=0;i<(INT4)sftvec->length;i++) {
  
    /* convert each SFT to complex timeseries */
    if (XLALSFTToCOMPLEX8TimeSeries(&((*dstimevec)->data[i]),&(sftvec->data[i]),&plan)) {
      LogPrintf(LOG_CRITICAL,"%s: XLALSFTToCOMPLEX8TimeSeries() failed with error = %d\n",__func__,xlalErrno);
      return XLAL_EINVAL;
    }

  }
  LogPrintf(LOG_DEBUG,"%s : performed inverse FFT on all %d SFTs\n",__func__,sftvec->length);

  /* free memory */
  XLALDestroyCOMPLEX8FFTPlan(plan);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/*******************************************************************************/
/** Inverse FFT all narrowband SFTs
 *
 * In order to apply the frequency derivitive corrections we must work in the time domain
 * so here we convert all SFTs to the complex time domain.
 *
 */
int XLALSFTToCOMPLEX8TimeSeries(COMPLEX8TimeSeries **ts,           /**< [out] the downsampled timeseries */
				COMPLEX8FrequencySeries *sft,      /**< [in] the input SFT vector */
				COMPLEX8FFTPlan **plan             /**< [in/out] the FFT plan */ 
				)
{
  /* INT4 i;                                /\* counter *\/ */
 /*  COMPLEX8FFTPlan *plan = NULL;          /\* inverse FFT plan *\/ */
  UINT4 N;                               /* the length of the SFTs */

  /* validate input arguments */
  if ((*ts) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output COMPLEX8TimeSeries structure != NULL.\n",__func__);
    return XLAL_EINVAL;
  }
  if (sft == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input SFT structure == NULL.\n",__func__);
    return XLAL_EINVAL;
  }

  /* we check that all input SFTs are of identical length so we make a single plan */
  N = sft->data->length;

  /* make the reverse plan if not cached */
  if ((*plan)==NULL) {
    if (((*plan) = XLALCreateReverseCOMPLEX8FFTPlan(N,1)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCreateReverseCOMPLEX8FFTPlan() failed with error = %d\n",__func__,xlalErrno);
      return XLAL_EINVAL;
    }
    LogPrintf(LOG_DEBUG,"%s : created the inverse FFT plan\n",__func__);
  }

  COMPLEX8Vector temp_output;
  REAL8 deltaF = sft->deltaF;
  REAL8 Tsft = 1.0/sft->deltaF;
  REAL8 deltaT = Tsft/N;
  UINT4 j;
   
  /* allocate output memory - create a COMPLEX8TimeSeries */
  if (((*ts) = XLALCreateCOMPLEX8TimeSeries("DS",&(sft->epoch),sft->f0,deltaT,&lalDimensionlessUnit,N)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCreateCOMPLEX8TimeSeries() failed to allocate memory for inverse FFT output.\n",__func__);
    return XLAL_ENOMEM;
  }

  /* point to input */
  COMPLEX8Vector temp_input;
  temp_input.length = sft->data->length;
  temp_input.data = (COMPLEX8*)sft->data->data;
  
  /* point temp output to timeseries */
  temp_output.length = N;
  temp_output.data = (COMPLEX8*)(*ts)->data->data;

  /* perform inverse FFT */
  if (XLALCOMPLEX8VectorFFT(&temp_output, &temp_input,(*plan))) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCOMPLEX8VectorFFT() failed with error = %d\n",__func__,xlalErrno);
    return XLAL_EINVAL;
  }

  /* normalise outputs by multiplying by df */
  for (j=0;j<N;j++) {
    temp_output.data[j].realf_FIXME = temp_output.data[j].realf_FIXME*deltaF;
    temp_output.data[j].imagf_FIXME = temp_output.data[j].imagf_FIXME*deltaF;
    /* fprintf(stdout,"%.12f %.12f %.12f\n",j*deltaT,temp_output.data[j].realf_FIXME,temp_output.data[j].imagf_FIXME); */
  }
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/*******************************************************************************/
/** Compute the gridding parameters on spin derivitives for all segments
 *
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

/*******************************************************************************/
/** Compute the gridding parameters on spin derivitives
 *
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
  REAL8 fnmin[NFREQMAX],fnmax[NFREQMAX]; /* min and max values of spin derivitives */
  INT4 dim[NFREQMAX];                    /* flag indicating whether a dimension has width */
  INT4 ndim = -1;                        /* the number of spin derivitive dimensions required */
  Template fdots;                        /* template for an instance of spin parameters */
  Template bintemp;                      /* template for instance of binary parameters */
  UINT4 ngrid = 100;                     /* the number of grid points per omega and tasc to search for finding true fdot spans per sft */

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

/*******************************************************************************/
/** Compute the demodulated power for all downsampled timeseries
 *
 * This function is simply a wrapper for XLALCOMPLEX8TimeSeriesToCOMPLEX8FrequencySeries()
 *
 */
int XLALCOMPLEX8TimeSeriesArrayToDemodPowerVector(REAL4DemodulatedPowerVector **power,     /**< [out] the spin derivitive demodulated power */
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
    
    COMPLEX8TimeSeries *ts = dsdata->data[i];
    COMPLEX8TimeSeries *temp_ts = NULL;
    COMPLEX8FrequencySeries *fs = NULL;
    UINT4 j;
    UINT4 idx = 0;
    
    /* allocate memory for all templates in the frequency derivitive space for this segment */
    if ( ((*power)->segment[i] = XLALCalloc(1,sizeof(REAL4Vector*))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
      return XLAL_ENOMEM;
    }
    if ( ((*power)->segment[i]->data = XLALCreateREAL4Vector(gridparams->segment[i]->max)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCreateREAL4Vector() failed with error = %d\n",__func__,xlalErrno);
      return XLAL_ENOMEM;
    }
    
    /* allocate memory for temporary complex timeseries */
    if ( (temp_ts = XLALCreateCOMPLEX8TimeSeries("DS",&(ts->epoch),ts->f0,ts->deltaT,&lalDimensionlessUnit,ts->data->length)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: XLALCreateCOMPLEX8TimeSeries() failed to allocate memory for inverse FFT output.\n",__func__);
      return XLAL_ENOMEM;
    }

    /* redefine grid on spin derivitives only */
    GridParameters tempgrid;
    Template *spintemp = NULL;
    tempgrid.ndim = gridparams->segment[i]->ndim - 1;
    tempgrid.grid = XLALCalloc(tempgrid.ndim,sizeof(Grid));
    tempgrid.prod = XLALCalloc(tempgrid.ndim,sizeof(UINT4));
    tempgrid.max = 1;
    for (j=0;j<tempgrid.ndim;j++) {
      tempgrid.grid[j].min = gridparams->segment[i]->grid[j+1].min;
      tempgrid.grid[j].delta = gridparams->segment[i]->grid[j+1].delta;
      tempgrid.grid[j].length = gridparams->segment[i]->grid[j+1].length;
      tempgrid.max *= tempgrid.grid[j].length;
    }
    tempgrid.prod[0] = 1;
    for (j=1;j<tempgrid.ndim;j++) {
      tempgrid.prod[j] = tempgrid.prod[j-1]*tempgrid.grid[j-1].length;
    }
    
    /* loop over spin derivitive values - not over frequency */
    ParameterSpace *temp = NULL;
    while (XLALGetNextTemplate(&spintemp,&tempgrid,temp)) {

      fprintf(stdout,"current idx = %d\n",spintemp->currentidx);

      /* reinitilaise the reused complex timeseries */
      fs = NULL;
 
      /* apply phase correction to complex timeseries */
      if (XLALApplyPhaseCorrection(&temp_ts,ts,spintemp)) {
	LogPrintf(LOG_CRITICAL,"%s : XLALApplyPhseCorrection() failed with error = %d\n",__func__,xlalErrno);
	XLAL_ERROR(XLAL_EINVAL);
      }
      LogPrintf(LOG_DEBUG,"%s : applied phase correction for template index %d/%d on segment %d/%d\n",__func__,spintemp->currentidx,tempgrid.max,i,dsdata->length);

      /* convert to the complex frequency domain - on the frequency grid specified */
      if (XLALCOMPLEX8TimeSeriesToCOMPLEX8FrequencySeries(&fs,temp_ts,&(gridparams->segment[i]))) {
	LogPrintf(LOG_CRITICAL,"%s : XLALSFTCOMPLEX8TimeseriesToCOMPLEX8FrequencySeries() failed with error = %d\n",__func__,xlalErrno);
	XLAL_ERROR(XLAL_EINVAL);
      }
      LogPrintf(LOG_DEBUG,"%s : computed demodulated frequency series for SFT %d/%d\n",__func__,i+1,dsdata->length);
      
      /* compute power and store it */
      for (j=0;j<fs->data->length;j++) {
	(*power)->segment[i]->data->data[idx] = crealf(fs->data->data[j])*crealf(fs->data->data[j]) + cimagf(fs->data->data[j])*cimagf(fs->data->data[j]);
	idx++;
      }
      (*power)->segment[i]->epoch.gpsSeconds = ts->epoch.gpsSeconds;
      (*power)->segment[i]->epoch.gpsNanoSeconds = ts->epoch.gpsNanoSeconds;
    
      /* free memory */
      XLALDestroyCOMPLEX8FrequencySeries(fs);

    } /* end loop over spin derivitive templates */
    
    /* free memory */
    XLALDestroyCOMPLEX8TimeSeries(temp_ts);

  } /* end the loop over segments */

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}

/** Compute the demodulated power for a single SFT
 *
 * This involves taking the downsampled SFT timeseries and multiplying by the
 * timeseries spin-derivitive templates in turn.  Then we inverse FFT the result
 * and square to obtain the power at a given set of freq derivitive parameters.
 *
 */
int XLALCOMPLEX8TimeSeriesToCOMPLEX8FrequencySeries(COMPLEX8FrequencySeries **fs,     /**< [out] the over-resolved frequency series */
						    const COMPLEX8TimeSeries *ts,     /**< [in] the downsampled SFT data */
						    GridParameters **gridparams        /**< [in/out] the spin derivitive gridding parameters */
						    )
{
 
  COMPLEX8FFTPlan *plan = NULL;           /* plan for the inverse FFT */
  COMPLEX8Vector *temp_input = NULL;      /* the temporary input of the inverse FFT = data*template */
  COMPLEX8Vector *temp_output = NULL;     /* the temporary output of the inverse FFT = data*template */
  UINT4 j;
  
  /* validate input arguments */
  if ((*fs) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output COMPLEX8FrequencySeries structure != NULL.\n",__func__);
    return XLAL_EINVAL;
  }
  if (ts == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input COMPLEX8TimeSeries structure == NULL.\n",__func__);
    return XLAL_EINVAL;
  }
  if (gridparams == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input GridParameters structure == NULL.\n",__func__);
    return XLAL_EINVAL;
  }

  UINT4 N = ts->data->length;                         /* the original length of the time series */
  REAL8 deltaT = ts->deltaT;                          /* the fixed time sampling of the input data */
  REAL8 T = N*deltaT;                                 /* the intrinsic duration */
  REAL8 deltaF = 1.0/(ts->deltaT*ts->data->length);   /* the intrinsic deltaF */

  /* define new deltaF accounting for the fixed deltaT */
  REAL8 newdeltaF = (*gridparams)->grid[0].delta;
/*   REAL8 oldmaxf = (*gridparams)->grid[0].min + (*gridparams)->grid[0].length*(*gridparams)->grid[0].delta; */
  UINT4 newN = ceil((T*deltaF/newdeltaF)/deltaT);
  REAL8 newT = deltaT*newN;
  newdeltaF = 1/newT;
   
  /* allocate memory for the temporary zero-padded input data */
  if ( (temp_input = XLALCreateCOMPLEX8Vector(newN)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCreateCOMPLEX8Vector() failed with error = %d.\n",__func__,xlalErrno);
    return XLAL_ENOMEM;
  }
  /* allocate memory for the temporary output data */
  if ( (temp_output = XLALCreateCOMPLEX8Vector(newN)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCreateCOMPLEX8Vector() failed with error = %d.\n",__func__,xlalErrno);
    return XLAL_ENOMEM;
  }
 
  /* create a forward complex fft plan */
  if ((plan = XLALCreateForwardCOMPLEX8FFTPlan(newN,0)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCreateForwardCOMPLEX8FFTPlan() failed with error = %d\n",__func__,xlalErrno);
    return XLAL_ENOMEM;
  }
  
  /* initialise the input data with zeros */
  memset(temp_input->data,0.0,temp_input->length*sizeof(COMPLEX8));

  /* put the input data into the temporary input structure and normalise it with deltaT */
  for (j=0;j<N;j++) {
    temp_input->data[j].realf_FIXME = ts->data->data[j].realf_FIXME*deltaT;
    temp_input->data[j].imagf_FIXME = ts->data->data[j].imagf_FIXME*deltaT;
  }

  /* FFT the data */
  if (XLALCOMPLEX8VectorFFT(temp_output,temp_input,plan)) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCOMPLEX8VectorFFT() failed with error = %d\n",__func__,xlalErrno);
    return XLAL_ENOMEM;
  }
  LogPrintf(LOG_DEBUG,"%s : computed the FFT\n",__func__);

  /* and modify grid based on new resolution  - just shift the lowest frequency to match up with existing bins */
  /* also update the resolution BUT keep the same number of bins */
  INT4 binoffset = (INT4)floor(((*gridparams)->grid[0].min - ts->f0)/newdeltaF);
  if (binoffset<0) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, we have a negative binoffset.\n",__func__);
    return XLAL_EINVAL;
  }
  (*gridparams)->grid[0].min = ts->f0 + newdeltaF*binoffset;
  (*gridparams)->grid[0].delta = newdeltaF;

  /* allocate memory for output */
  if ( ((*fs) = XLALCreateCOMPLEX8FrequencySeries("FS",&(ts->epoch),(*gridparams)->grid[0].min,(*gridparams)->grid[0].delta,&lalDimensionlessUnit,(*gridparams)->grid[0].length)) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: XLALCreateREAL4Vector() failed with error = %d\n",__func__,xlalErrno);
    return XLAL_ENOMEM;
  }

  /*  extract desired frequencies */
  for (j=0;j<(*gridparams)->grid[0].length;j++) {
    INT4 k = j + binoffset;
    (*fs)->data->data[j].realf_FIXME = temp_output->data[k].realf_FIXME;
    (*fs)->data->data[j].imagf_FIXME = temp_output->data[k].imagf_FIXME;
  }

  /* free memory */
  XLALDestroyCOMPLEX8Vector(temp_input);
  XLALDestroyCOMPLEX8Vector(temp_output);
  XLALDestroyCOMPLEX8FFTPlan(plan);
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/** Compute the demodulated power for a single SFT
 *
 * This involves taking the downsampled SFT timeseries and multiplying by the
 * timeseries spin-derivitive templates in turn.  Then we inverse FFT the result
 * and square to obtain the power at a given set of freq derivitive parameters.
 *
 */
int XLALApplyPhaseCorrection(COMPLEX8TimeSeries **outts,            /**< [out] the output complex time series */
			     COMPLEX8TimeSeries *ints,          /**< [in] the input complex time series */
			     Template *fn                       /**< [in] the spin derivitives */
			     )
{
  
  UINT4 j,k;

  /* validate input arguments */
  if ((*outts) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output COMPLEX8Vector structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (ints == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input COMPLEX8 structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  
  /* apply time domain phase correction - first loop over time and then spin derivitive */
  for (j=0;j<ints->data->length;j++) {
    
    /* compute phase correction including heterodyne to shift frequencies to match up with grid */
    REAL8 tn = j*ints->deltaT - 0.5*ints->deltaT*ints->data->length;
    REAL8 arg = 0.0;
    UINT4 fac = 1;
    
    /* loop over each spin derivitive and add to phase contribution for current time sample */
    for (k=1;k<fn->ndim;k++) {
      tn *= tn;
      fac *= k+1;
      arg += (-1.0)*LAL_TWOPI*fn->x[k]*tn/fac;
    }
        
    /* compute real and imaginary parts of phase correction timeseries */
    REAL8 xr = cos(arg);
    REAL8 xi = sin(arg);
    
    /* multiply data by phase correction - leave the zero-padding */
    (*outts)->data->data[j].realf_FIXME = crealf(ints->data->data[j])*xr - cimagf(ints->data->data[j])*xi;
    (*outts)->data->data[j].imagf_FIXME = crealf(ints->data->data[j])*xi + cimagf(ints->data->data[j])*xr;
    
  }
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}



/** Compute the background photon flux for each SFT
 *
 * This computes the median of the power to obtain the average counts per second
 *
 */
int XLALEstimateBackgroundFlux(REAL8Vector **background,     /**< [out] the background flux estimate */
			       SFTVector *sftvec             /**< [in/out] the SFTs */
			       )
{
  LALStatus status = blank_status;        /* for use wih non-XLAL functions */
  UINT4 i,j;                               /* counters */

  /* validate input arguments */
  if ((*background) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output REAL8Vector structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (sftvec == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input SFTVector structure == NULL.\n",__func__);
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
  /*   REAL8 T = 1.0/sftvec->data[i].deltaF; */
    REAL8 R = 0.0;

    /* allocate temporary memory */
    if ((P = XLALCalloc(sft->length,sizeof(REAL8))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for gridparams->grid.\n",__func__);
      XLAL_ERROR(XLAL_ENOMEM);
    }

    /* loop over each element in the SFT and record the power */
    for (j=0;j<sft->length;j++) {
      P[j] = (crealf(sft->data[j])*crealf(sft->data[j]) + cimagf(sft->data[j])*cimagf(sft->data[j]));
      /* fprintf(stdout,"%.12f %.12f %.12f\n",sftvec->data[0].f0 + j*sftvec->data[i].deltaF,crealf(sft->data[j]),cimagf(sft->data[j])); */
    }
    
    /* sort the data */
    gsl_sort(P,1,sft->length);
 
    /* compute median */
    median = gsl_stats_median_from_sorted_data(P,1,sft->length);
   
    /* compute the median bias */
    LAL_CALL ( LALRngMedBias( &status, &medianbias, sft->length ), &status);
   
    /* record estimate - this is the expected mean of the power spectrum */
    R = (median/medianbias);
    (*background)->data[i] = R; /* T*R/pow((REAL8)SAMPFREQ,2.0); */
    LogPrintf(LOG_DEBUG,"%s : Estimated the background for SFT starting at %d as %.3e cnts/s.\n",__func__,sftvec->data[i].epoch.gpsSeconds,R);
    fprintf(stdout,"sctually using %f a norm\n",(*background)->data[i]);

    /* free the power */
    XLALFree(P);
    
  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/** Compute the grid on binary parameters based on the semi-coherent metric
 *
 * We use this grid to perform the integration of the posterior and to ultimately
 * compute the Bayes factor.
 *
 */
int XLALComputeBinaryGridParams(GridParameters **binarygridparams,  /**< [out] the binary parameter grid */
				REAL8Space *space,                  /**< [in] the signal parameter space */
				REAL8 T,                            /**< [in] the duration of the observation */
				REAL8 DT,                           /**< [in] the length of the coherent segments */
				REAL8 mu,                           /**< [in] the mismatch */
				REAL8 coverage
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
  
  /* if we've specified a random template bank */
  if (coverage>0) {

    REAL8 Vn = pow(LAL_PI,ndim/2.0)/gsl_sf_gamma(1.0+ndim/2.0);
  
    REAL8 G11 = pow(LAL_PI,2.0)*DT*DT/3;
    REAL8 G22 = pow(LAL_PI,2.0)*DT*DT/6;
    REAL8 G33 = pow(LAL_PI,2.0)*DT*DT/6;
    REAL8 G44 = pow(LAL_PI,2.0)*DT*DT*T*T/72;
    
    REAL8 dVr = sqrt(G11*G22*G33*G44);
    REAL8 Vsr = dVr*space->data[2].span*(1.0/60.0)*(pow(space->data[3].max,5.0)-pow(space->data[3].min,5.0))*(pow(space->data[0].max,4.0)-pow(space->data[0].min,4.0))*(pow(space->data[1].max,3.0)-pow(space->data[1].min,3.0));

    (*binarygridparams)->Nr = (INT4)ceil((1.0/Vn)*log(1.0/(1.0-coverage))*(pow(mu,-ndim/2.0))*Vsr);
    LogPrintf(LOG_DEBUG,"%s : computed the number of random binary templates to be %d.\n",__func__,(*binarygridparams)->Nr);
    LogPrintf(LOG_DEBUG,"%s : to br compared to the total number of cubic templates %d (%.6f).\n",__func__,(*binarygridparams)->max,(REAL8)(*binarygridparams)->max/(REAL8)(*binarygridparams)->Nr);
 
  }
  else (*binarygridparams)->Nr = -1;

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/** Compute the semi-coherent statistics
 *
 * This function computes the semi-coherent s.
 *
 */
int XLALComputeSemiCoherentStat(FILE *fp,                                /**< [in] the output file pointer */
				REAL4DemodulatedPowerVector *power,      /**< [in] the input data in the form of power */
				ParameterSpace *pspace,                  /**< [in] the parameter space */
				GridParametersVector *fgrid,
				GridParameters *bingrid,                 /**< [in] the grid parameters */
				REAL8Vector *background,
				REAL8 frac
				)
{

  Template *bintemp = NULL;                           /* the binary parameter space template */
  Template fdots;                                     /* the freq derivitive template for each segment */
  UINT4 i,j;                                          /* counters */
  UINT4 percent = 0;                                  /* counter for status update */
 
  /* validate input parameters */
 /*  if ((*SemiCo) != NULL) { */
/*     LogPrintf(LOG_CRITICAL,"%s: Invalid input, output BayesianProducts structure != NULL.\n",__func__); */
/*     XLAL_ERROR(XLAL_EINVAL); */
/*   } */
  if (power == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input REAL4DemodulatedPowerVector structure = NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (pspace == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input GridParameters structure = NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (fgrid == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input GridParametersVector structure = NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (bingrid == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input GridParameters structure = NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
 /*  /\* allocate memory for the results *\/ */
/*   if ((*SemiCo = XLALCreateREAL8Vector(bingrid->max)) == NULL) { */
/*     LogPrintf(LOG_CRITICAL,"%s : XLALCrateREAL8Vector() failed with error = %d\n",__func__,xlalErrno); */
/*     XLAL_ERROR(XLAL_ENOMEM); */
/*   } */

  /* allocate memory for the fdots */
  if ((fdots.x = XLALCalloc(fgrid->segment[0]->ndim,sizeof(REAL8))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  fdots.ndim = fgrid->segment[0]->ndim;

  /* define the chi-squared threshold */
  fprintf(stdout,"frac = %f\n",frac);
  REAL8 thr = gsl_cdf_chisq_Qinv(frac,2*power->length);
  LogPrintf(LOG_DEBUG,"%s : computed the threshold as %f\n",__func__,thr);

  int (*getnext)(Template **temp,GridParameters *gridparams, ParameterSpace *space);
  fprintf(stdout,"bingrid->max = %d\n",bingrid->max);
  INT4 newmax = bingrid->max;
  ParameterSpace *temppspace = NULL;
  fprintf(stdout,"newmax = %d\n",newmax);
  fprintf(stdout,"bingrid->Nr = %d\n",bingrid->Nr);
  if (bingrid->Nr>0) {
    fprintf(stdout,"inside\n");
    getnext = &XLALGetNextRandomBinaryTemplate;
    newmax = bingrid->Nr;
    temppspace = pspace;
  }
  else getnext = &XLALGetNextTemplate;
  fprintf(stdout,"max = %d Nr = %d\n",newmax,bingrid->Nr);
  /* single loop over binary templates */
  while (getnext(&bintemp,bingrid,temppspace)) {

    REAL8 logLratiosum = 0.0;                       /* initialise likelihood ratio */

    /** loop over segments **********************************************************************************/
    for (i=0;i<power->length;i++) {
      
      REAL4DemodulatedPower *currentpower = power->segment[i];
      GridParameters *fdotgrid = fgrid->segment[i];
      REAL8 tmid = XLALGPSGetREAL8(&(power->segment[i]->epoch)) + 0.5*pspace->tseg;
      REAL8 norm = background->data[i];
      INT4 idx = 0;
    
      /* compute instantaneous frequency derivitives corresponding to the current template for this segment */
      XLALComputeBinaryFreqDerivitives(&fdots,bintemp,tmid);
       
      /* find indices corresponding to the spin derivitive values for the segment power */
      for (j=0;j<fdots.ndim;j++) {
	UINT4 tempidx = 0.5 + (fdots.x[j] - fdotgrid->grid[j].min)*fdotgrid->grid[j].oneoverdelta;
	idx += tempidx*fdotgrid->prod[j];
      }
      
      /* define the power at this location in this segment */
      logLratiosum += currentpower->data->data[idx]/norm;
       
    } /* end loop over segments */
    /*************************************************************************************/
   
    /* output semi-coherent statistic for this template if it exceeds the threshold */
    logLratiosum *= 2.0;     /* make it a true chi-squared variable */
    if (logLratiosum>thr) {
      for (j=0;j<bingrid->ndim;j++) fprintf(fp,"%6.12f\t",bintemp->x[j]);
      fprintf(fp,"%6.12e\n",logLratiosum);
    }
  
    /* output status to screen */
    if (floor(100.0*(REAL8)bintemp->currentidx/(REAL8)newmax) > (REAL8)percent) {
      percent = (UINT4)floor(100*(REAL8)bintemp->currentidx/(REAL8)newmax);
      LogPrintf(LOG_DEBUG,"%s : completed %d%% (%d/%d)\n",__func__,percent,bintemp->currentidx,newmax);
    }
   
  } /* end loop over templates */
  /*************************************************************************************/

  /* free template memory */
  XLALFree(fdots.x);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}




/** Compute the next template in a grid
 *
 * the templates are generated sequentially from the grid parameters file.  The n-dimensional
 * virtual indices of each dimension are also generated
 *
 */
int XLALGetNextTemplate(Template **temp,                        /**< [out] the signal model template parameters */
			GridParameters *gridparams,              /**< [in] the parameter space grid params */
			ParameterSpace *space
			)
{
  UINT4 idx;                             /* the index variable */
  INT4 j;                                /* counters */

 
  /* check input */
  if (space != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input ParameterSpace structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

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
  else if ((*temp)->currentidx == gridparams->max - 1) {
    
    /* free binary template memory */
    XLALFree((*temp)->x);
    XLALFree((*temp)->idx);
    XLALFree(*temp);
    
    LogPrintf(LOG_DEBUG,"%s: at last template.\n",__func__);
    return 0;
  }
  else (*temp)->currentidx++;
    
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
 /*  (*temp)->currentidx++; */
    

  return 1;

}

/** Compute the next template in a grid
 *
 * the templates are generated sequentially from the grid parameters file.  The n-dimensional
 * virtual indices of each dimension are also generated
 *
 */
int XLALGetNextRandomBinaryTemplate(Template **temp,                        /**< [out] the signal model template parameters */
				    GridParameters *gridparams,             /**< [in] the parameter space grid params */
				    ParameterSpace *space
				    )
{
  
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
  else if ((*temp)->currentidx == (UINT4)gridparams->Nr - 1) {
    
    /* free binary template memory */
    XLALFree((*temp)->x);
    XLALFree((*temp)->idx);
    XLALFree(*temp);
    
    LogPrintf(LOG_DEBUG,"%s: at last template.\n",__func__);
    return 0;
  }
  else (*temp)->currentidx++;
    
  REAL8 n1 = space->space->data[0].min;
  REAL8 n2 = space->space->data[0].max;
  REAL8 a1 = space->space->data[1].min;
  REAL8 a2 = space->space->data[1].max;
  REAL8 t1 = space->space->data[2].min;
  REAL8 t2 = space->space->data[2].max;
  REAL8 O1 = space->space->data[3].min;
  REAL8 O2 = space->space->data[3].max;
  REAL8 nu;
  REAL8 Om;
  REAL8 a;
  REAL8 tasc;

  /* while we haven't selected a template in the space */
  INT4 flag = 0;
  while (!flag) {
    
    REAL8 temp1 = gsl_ran_flat(r,0,1);
    nu = n2*pow((pow(n1/n2,4.0) + (1.0 - pow(n1/n2,4.0))*temp1),1.0/4.0);
    REAL8 temp2 = gsl_ran_flat(r,0,1);
    a = a2*pow((pow(a1/a2,3.0) + (1.0 - pow(a1/a2,3.0))*temp2),1.0/3.0);
    REAL8 temp3 = gsl_ran_flat(r,0,1);
    tasc = t1 + (t2-t1)*temp3;
    REAL8 temp4 = gsl_ran_flat(r,0,1);
    Om = O2*pow((pow(O1/O2,5.0) + (1.0 - pow(O1/O2,5.0))*temp4),1.0/5.0);

    /* fprintf(stdout,"%f (%f %f) %f (%f %f) %f (%f %f) %f (%f %f)\n",nu,n1,n2,a,a1,a2,tasc,t1,t2,Om,O1,O2); */
    
    if ((nu>n1)&&(nu<n2)&&(a>a1)&&(a<a2)&&(tasc>t1)&&(tasc<t2)&&(Om>O1)&&(Om<O2)) flag = 1;

  }
  
  (*temp)->x[0] = nu;
  (*temp)->x[1] = a;
  (*temp)->x[2] = tasc;
  (*temp)->x[3] = Om;
   
  return 1;

}

/** Compute the instantaneous frequency derivitives for a given binary template and segment
 *
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

/** Output the results to file
 *
 * We choose to output all results from a specific analysis to a single file
 *
 */
int XLALOpenSemiCoherentResultsFile(FILE **fp,                   /**< [in] filepointer to output file */
				    CHAR *outputdir,            /**< [in] the output directory name */
				    ParameterSpace *pspace,     /**< [in] the parameter space */
				    CHAR *clargs,               /**< [in] the command line args */
				    UserInput_t *uvar
				    )
{
  CHAR outputfile[LONGSTRINGLENGTH];    /* the output filename */
/*   Template *bintemp = NULL;             /\* the binary parameter space template *\/   */
  time_t curtime = time(NULL);          /* get the current time */
  CHAR *time_string = NULL;             /* stores the current time */
  CHAR *version_string = NULL;          /* pointer to a string containing the git version information */
 /*  UINT4 j;                            /\* counters *\/ */
/*   UINT4 cnt = 0;                         */

  /* validate input */
  if (outputdir == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output directory string == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
 /*  if (SemiCo == NULL) { */
/*     LogPrintf(LOG_CRITICAL,"%s: Invalid input, results BayesProducts structure == NULL.\n",__func__); */
/*     XLAL_ERROR(XLAL_EINVAL); */
/*   } */
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
    snprintf(outputfile,LONGSTRINGLENGTH,"%s/SemiCoherentResults-%s-%d_%d-%04d_%03d_%04d_%03d.txt",
	     outputdir,uvar->comment,pspace->epoch.gpsSeconds,end,min_freq_int,min_freq_mhz,max_freq_int,max_freq_mhz);
  }
  LogPrintf(LOG_DEBUG,"%s : output %s\n",__func__,outputfile);

  /* open the output file */
  if (((*fp) = fopen(outputfile,"w")) == NULL) {
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
  fprintf((*fp),"%s \n",version_string);
  fprintf((*fp),"%%%% command line args\t\t= %s\n",clargs);
  fprintf((*fp),"%%%% filename\t\t\t\t= %s\n",outputfile);
  fprintf((*fp),"%%%% date\t\t\t\t\t= %s\n",time_string);
  fprintf((*fp),"%%%% start time (GPS sec)\t\t= %d\n",pspace->epoch.gpsSeconds);
  fprintf((*fp),"%%%% observation span (sec)\t= %d\n",(UINT4)pspace->span);
  fprintf((*fp),"%%%% coherent time (sec)\t\t= %d\n",(UINT4)pspace->tseg);
 /*  fprintf(fp,"%%%% number of segments\t\t= %d\n",Bayes->nsegments); */
/*   fprintf(fp,"%%%% number of dimensions\t= %d\n",Bayes->gridparams->ndim); */
/*   if (pspace->ampspace) fprintf(fp,"%%%% amplitude dimension\t\t\t= 1\n"); */
/*   else fprintf(fp,"%%%% amplitude dimension\t\t\t= 0\n"); */
/*   fprintf(fp,"%%%% mismatch\t\t\t\t= %6.12f\n",Bayes->gridparams->mismatch); */
  fprintf((*fp),"%%%%\n");

 /*  /\* sort the results and find the threshold corresponding to the top n'th value *\/ */
/*   size_t *p = NULL; */
/*   if ((p = LALCalloc(SemiCo->length,sizeof(size_t))) == NULL) { */
/*     LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for index structure.\n",__func__); */
/*     XLAL_ERROR(XLAL_ENOMEM); */
/*   } */
/*   gsl_sort_index(p,SemiCo->data,1,SemiCo->length); */
/*   REAL8 thr = SemiCo->data[p[SemiCo->length-top+1]]; */
/*   fprintf(stdout,"thresh = %.12f\n",thr); */

/*   /\* loop over each fine grid template and output the coordinates and the semi-coherent statistic *\/ */
/*   while (XLALGetNextTemplate(&bintemp,bingrid)) { */
  
/*     if (SemiCo->data[cnt]>=thr) { */
/*       for (j=0;j<bingrid->ndim;j++) fprintf(fp,"%6.12f\t",bintemp->x[j]); */
/*       fprintf(fp,"%6.12e\n",SemiCo->data[cnt]); */
/*     } */
/*     cnt++; */
  
/*   } */

/*   /\* close the file *\/ */
/*   fclose(fp); */

  /* free memory */
  XLALFree(time_string);
  XLALFree(version_string);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}



/** Free the memory allocated within a ParameterSpace structure
 *
 */
int XLALFreeParameterSpace(ParameterSpace *pspace            /**< [in] the parameter space to be freed */
			   )
{
 /*  UINT4 i;                     /\* counter *\/ */

  /* free parameter space */
  XLALFree(pspace->space->data);
  XLALFree(pspace->space);
  
  /* /\* free prior params *\/ */
/*   for (i=0;i<pspace->gridparams->ndim;i++) { */
/*     XLALDestroyREAL8Vector(pspace->priors->data[i].logpriors); */
/*   } */
/*   XLALFree(pspace->priors->data); */
/*   XLALFree(pspace->priors); */
/*   LogPrintf(LOG_DEBUG,"%s : freed the prior parameters\n",__func__); */
  
  /* free binary grid params */
 /*  XLALFree(pspace->gridparams->grid); */
/*   XLALFree(pspace->gridparams->prod); */
/*   XLALFree(pspace->gridparams); */
/*   LogPrintf(LOG_DEBUG,"%s : freed the binary grid parameters\n",__func__); */

  /* /\* free the injection parameters if used *\/ */
/*   if (pspace->inj) { */
/*     XLALFree(pspace->inj->temp.x); */
/*     XLALFree(pspace->inj); */
/*   } */

  /* /\* free amplitude parameters *\/ */
/*   if (pspace->ampspace) XLALFree(pspace->ampspace); */
/*   if (pspace->ampgrid) XLALFree(pspace->ampgrid); */
/*   if (pspace->amppriors) { */
/*     XLALDestroyREAL8Vector(pspace->amppriors->logpriors); */
/*     XLALFree(pspace->amppriors); */
/*   } */

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

/** Free the memory allocated within a REAL4DemodulatedPowerVector structure
 *
 */
int XLALFreeREAL4DemodulatedPowerVector(REAL4DemodulatedPowerVector *power            /**< [in] the data to be freed */
					)
{
  UINT4 i;                     /* counter */

  /* free each segment */
  for (i=0;i<power->length;i++) {

    /* if there is a non-null gridparams structure then free it aswell */
   /*  if (power->segment[i]->gridparams) { */
/*       XLALFree(power->segment[i]->gridparams->grid); */
/*       XLALFree(power->segment[i]->gridparams->prod); */
/*       XLALFree(power->segment[i]->gridparams); */
/*     } */
    XLALDestroyREAL4Vector(power->segment[i]->data);
    XLALFree(power->segment[i]);
  }
  XLALFree(power->segment);
  XLALFree(power);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}



/* /\** this function initialises the gsl random number generation */
/*  * */
/*  * If the input seed is zero then a random seed is drawn from */
/*  * /dev/urandom. */
/*  * */
/*  *\/ */
/* int XLALInitgslrand(gsl_rng **gslrnd,     /\**< [out] the gsl random number generator *\/ */
/* 		    INT8 seed             /\**< [in] the random number generator seed *\/ */
/* 		    ) */
/* { */
/*   FILE *devrandom = NULL;      /\* pointer to the /dev/urandom file *\/ */
  
/*   /\* if the seed is 0 then we draw a random seed from /dev/urandom *\/ */
/*   if (seed == 0) { */
    
/*     /\* open /dev/urandom *\/ */
/*     if ((devrandom=fopen("/dev/urandom","r")) == NULL)  { */
/*       LogPrintf(LOG_CRITICAL,"%s: Error, unable to open device /dev/random\n",__func__); */
/*       XLAL_ERROR(XLAL_EINVAL); */
/*     } */
    
/*     /\* read a random seed *\/ */
/*     if (fread((void*)&seed,sizeof(INT8),1,devrandom) != 1) { */
/*       LogPrintf(LOG_CRITICAL,"%s: Error, unable to read /dev/random\n",__func__); */
/*       XLAL_ERROR(XLAL_EINVAL); */
/*     } */
/*     fclose(devrandom); */
    
/*   } */
  
/*   /\* setup gsl random number generation *\/ */
/*   *gslrnd = gsl_rng_alloc(gsl_rng_taus2); */
/*   gsl_rng_set(*gslrnd,seed); */
 
/*   LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__); */
/*   return XLAL_SUCCESS; */
  
/* } */


/* /\* TESTING *\/ */
/*   if (0) */
/*     { */

/*       COMPLEX8TimeSeries *ts = NULL; */
/*       COMPLEX8FrequencySeries *fs = NULL; */
/*       COMPLEX8FFTPlan *plan = NULL; */
/*       FILE *fp = NULL; */
/*       UINT4 j; */

/*       /\* output SFT to text file *\/ */
/*       if ((fp = fopen("/Users/chrismessenger/temp/sft.txt","w"))==NULL) { */
/* 	LogPrintf(LOG_CRITICAL,"%s: Couldn't open file. Failed with error = %d\n",__func__,xlalErrno); */
/* 	return XLAL_EINVAL; */
/*       } */
/*       for (j=0;j<sftvec->data[0].data->length;j++) fprintf(fp,"%.12f %.12f %.12f\n",sftvec->data[0].f0 + j*sftvec->data[0].deltaF,crealf(sftvec->data[0].data->data[j]),cimagf(sftvec->data[0].data->data[j])); */
/*       fclose(fp); */
/*       LogPrintf(LOG_DEBUG,"%s : output first SFT\n",__func__); */

/*       /\* convert single SFT to complex timeseries *\/ */
/*       if (XLALSFTToCOMPLEX8TimeSeries(&ts,&(sftvec->data[0]),&plan)) { */
/* 	LogPrintf(LOG_CRITICAL,"%s : XLALSFTtoCOMPLEX8Timeseries() failed with error = %d\n",__func__,xlalErrno); */
/* 	return 1; */
/*       } */
/*       LogPrintf(LOG_DEBUG,"%s : converted SFT to complext timeseries\n",__func__); */

/*       /\* output timeseries to file *\/  */
/*       if ((fp = fopen("/Users/chrismessenger/temp/complexts.txt","w"))==NULL) { */
/* 	LogPrintf(LOG_CRITICAL,"%s: Couldn't open file. Failed with error = %d\n",__func__,xlalErrno); */
/* 	return XLAL_EINVAL; */
/*       } */
/*       for (j=0;j<ts->data->length;j++) fprintf(fp,"%.12f %.12f %.12f\n",j*ts->deltaT,crealf(ts->data->data[j]),cimagf(ts->data->data[j])); */
/*       fclose(fp); */
/*       LogPrintf(LOG_DEBUG,"%s : output first complex timeseries\n",__func__); */

/*       /\* compute over-resolved frequency series *\/ */
/*       if (XLALCOMPLEX8TimeSeriesToCOMPLEX8FrequencySeries(&fs,ts,&(freqgridparams->segment[0]))) { */
/* 	LogPrintf(LOG_CRITICAL,"%s : XLALSFTCOMPLEX8TimeseriesToCOMPLEX8FrequencySeries() failed with error = %d\n",__func__,xlalErrno); */
/* 	return 1; */
/*       } */
/*       LogPrintf(LOG_DEBUG,"%s : converted complext timeseries back to frequency domain\n",__func__); */

/*       /\* output timeseries to file *\/  */
/*       if ((fp = fopen("/Users/chrismessenger/temp/complexfs.txt","w"))==NULL) { */
/* 	LogPrintf(LOG_CRITICAL,"%s: Couldn't open file. Failed with error = %d\n",__func__,xlalErrno); */
/* 	return XLAL_EINVAL; */
/*       } */
/*       for (j=0;j<fs->data->length;j++) fprintf(fp,"%.12f %.12f %.12f\n",fs->f0 + j*fs->deltaF,crealf(fs->data->data[j]),cimagf(fs->data->data[j])); */
/*       fclose(fp); */
/*       LogPrintf(LOG_DEBUG,"%s : output converted frequency series\n",__func__); */
      
/*       XLALDestroyCOMPLEX8TimeSeries(ts); */
/*       XLALDestroyCOMPLEX8FrequencySeries(fs); */

/*       exit(0); */
/*     } */

/** this function initialises the gsl random number generation
 *
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
