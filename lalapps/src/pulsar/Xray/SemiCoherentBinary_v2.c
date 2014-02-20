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
#include <sys/stat.h>
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

#include "SemiCoherent.h"

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
  CHAR *comment;
  CHAR *tempdir;                    /**< a temporary directory for keeping the results */
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
int XLALComputeSemiCoherentStat(FILE *fp,REAL4DemodulatedPowerVector *power,ParameterSpace *pspace,GridParametersVector *fgrid,GridParameters *bingrid,REAL8Vector *background,REAL8 top);

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
  CHAR newnewtemp[LONGSTRINGLENGTH];
 /*  REAL8Vector *SemiCo = NULL;                   /\* the semi-coherent statistic results *\/    */
  REAL8 fmin_read,fmax_read,fband_read;         /* the range of frequencies to be read from SFTs */
  UINT4 i;                                      /* counter */
  FILE *sfp = NULL;
  FILE *cfp = NULL;

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

  /* initialise the random number generator */
  if (XLALInitgslrand(&r,0)) {
    LogPrintf(LOG_CRITICAL,"%s: XLALinitgslrand() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_EFAULT);
  }

  /* make temporary directory */
  if (uvar.tempdir) {

    /* initialise the random number generator  - use the clock */
    gsl_rng * q;
    if (XLALInitgslrand(&q,0)) {
      LogPrintf(LOG_CRITICAL,"%s: XLALinitgslrand() failed with error = %d\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_EFAULT);
    }

    CHAR newtemp[LONGSTRINGLENGTH];
    INT4 id = (INT4)(1e9*gsl_rng_uniform(q));
    sprintf(newtemp,"%s/%09d",uvar.tempdir,id); 
    if (mkdir(newtemp,0755)) {
      LogPrintf(LOG_DEBUG,"%s : Unable to make temporary directory %s.  Might be a problem.\n",__func__,newtemp);
    }
    sprintf(newnewtemp,"%s/%.3f-%.3f",newtemp,uvar.freq,uvar.freq+uvar.freqband);
    if (mkdir(newnewtemp,0755)) {
      LogPrintf(LOG_CRITICAL,"%s : Unable to make temporary directory %s\n",__func__,newnewtemp);
      return 1;
    }
  }

  /* initialise the random number generator */
  if (XLALInitgslrand(&r,uvar.seed)) {
    LogPrintf(LOG_CRITICAL,"%s: XLALinitgslrand() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_EFAULT);
  }

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
  /* OPEN INTERMEDIATE RESULTS FILE */
  /**********************************************************************************/
  
  if (XLALOpenSemiCoherentResultsFile(&cfp,newnewtemp,&pspace,clargs,&uvar,1)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALOpenCoherentResultsFile() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : opened coherent results file.\n",__func__);

   /**********************************************************************************/
  /* COMPUTE THE STATISTICS ON THE COARSE GRID */
  /**********************************************************************************/

  /* compute the demodulated power on the frequency derivitive grid */
  if (XLALCOMPLEX8TimeSeriesArrayToDemodPowerVector(&dmpower,dstimevec,freqgridparams,background,cfp)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCOMPLEX8TimeSeriesArrayToDemodPowerVector() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  fclose(cfp);
  LogPrintf(LOG_DEBUG,"%s : computed the demodulated power\n",__func__);

  /**********************************************************************************/
  /* OPEN RESULTS FILE */
  /**********************************************************************************/
  
  if (XLALOpenSemiCoherentResultsFile(&sfp,newnewtemp,&pspace,clargs,&uvar,0)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALOutputBayesResults() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : output results to file.\n",__func__);

   /**********************************************************************************/
  /* COMPUTE THE STATISTICS ON THE FINE GRID */
  /**********************************************************************************/

  /* compute the semi-coherent detection statistic on the fine grid */
  if (XLALComputeSemiCoherentStat(sfp,dmpower,&pspace,freqgridparams,bingridparams,background,uvar.frac)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALComputeSemiCoherentStat() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  fclose(sfp);
  LogPrintf(LOG_DEBUG,"%s : computed the semi-coherent statistic\n",__func__);

  /**********************************************************************************/
  /* CLEAN UP */
  /**********************************************************************************/

  /* move the temporary directory to the final location */
  if (uvar.tempdir) {
    CHAR newoutputdir[LONGSTRINGLENGTH];
    sprintf(newoutputdir,"%s/%.3f-%.3f",uvar.outputdir,uvar.freq,uvar.freq+uvar.freqband);
    if (rename(newnewtemp,newoutputdir)) {
      LogPrintf(LOG_CRITICAL,"%s : unable to move final results directory %s -> %s.  Exiting.\n",newnewtemp,newoutputdir);
      return 1;
    }
  }

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
  uvar->tempdir = NULL;

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
  XLALregSTRINGUserStruct(tempdir,              'z', UVAR_OPTIONAL, "A temporary directory");
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
  REAL8 thr = gsl_cdf_chisq_Qinv(frac,2*power->length);
  LogPrintf(LOG_DEBUG,"%s : computed the threshold as %f\n",__func__,thr);

  int (*getnext)(Template **temp,GridParameters *gridparams, ParameterSpace *space);
  INT4 newmax = bingrid->max;
  ParameterSpace *temppspace = NULL;
  if (bingrid->Nr>0) {
    getnext = &XLALGetNextRandomBinaryTemplate;
    newmax = bingrid->Nr;
    temppspace = pspace;
  }
  else getnext = &XLALGetNextTemplate;

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


