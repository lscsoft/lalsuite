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
 * \defgroup lalapps_pulsar_Xray X-ray Search Applications
 * \ingroup lalapps_pulsar_Apps
 */

/** \author C.Messenger
 * \ingroup lalapps_pulsar_Xray
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
#include "config.h"
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
/* define internal structures */

/** A structure that stores user input variables
 */
typedef struct {
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
  INT4 gpsstart;                    /**< the min GPS time to include */
  INT4 gpsend;                      /**< the max GPS time to include */
  INT4 seed;                        /**< fix the random number generator seed */
  REAL8 coverage;                   /**< random template bank coverage */
  INT4 blocksize;                  /**< the running median blocksize */
  INT4 ntoplist;                   /**< the number of results to record */
  INT4 tsft;			   /**< the length of the input sfts */
  CHAR *comment;
  CHAR *tempdir;                    /**< a temporary directory for keeping the results */
  BOOLEAN version;	            /**< output version-info */
} UserInput_t;

typedef struct {
  REAL8 *data;
  INT4 idx;
  INT4 n;
  REAL8 **params;
} toplist;

/***********************************************************************************************/
/* global variables */
extern int vrbflg;	 	/**< defined in lalapps.c */

/***********************************************************************************************/
/* define functions */
int main(int argc,char *argv[]);
int XLALReadUserVars(int argc,char *argv[],UserInput_t *uvar, CHAR **clargs);
int XLALComputeSemiCoherentStat(FILE *fp,REAL4DemodulatedPowerVector *power,ParameterSpace *pspace,GridParametersVector *fgrid,GridParameters *bingrid,INT4 ntoplist);
int XLALDefineBinaryParameterSpace(REAL8Space **, LIGOTimeGPS, REAL8, UserInput_t *);
int XLALOpenSemiCoherentResultsFile(FILE **,CHAR *,ParameterSpace *,CHAR *,UserInput_t *);
int XLALtoplist(REAL8 x,Template *params,toplist *TL);

/***********************************************************************************************/
/* empty initializers */
UserInput_t empty_UserInput;
ParameterSpace empty_ParameterSpace;
gsl_rng *r = NULL;

/** The main function of semicoherentbinary.c
 *
 */
int main( int argc, char *argv[] )  {

  UserInput_t uvar = empty_UserInput;           /* user input variables */
  CHAR *clargs = NULL;                          /* store the command line args */
  SFTVector *sftvec = NULL;                     /* stores the input SFTs */
  /* REAL4VectorArray *background = NULL;                  /\* running median estimates of the background for each SFT *\/ */
  ParameterSpace pspace = empty_ParameterSpace;    /* the search parameter space */
  COMPLEX8TimeSeriesArray *dstimevec = NULL;       /* contains the downsampled inverse FFT'd SFTs */
  REAL4DemodulatedPowerVector *dmpower = NULL;    /* contains the demodulated power for all SFTs */
  GridParametersVector *freqgridparams = NULL;  /* the coherent grid on the frequency derivitive parameter space */
  GridParameters *bingridparams = NULL;
  CHAR newnewtemp[LONGSTRINGLENGTH];
 /*  REAL8Vector *SemiCo = NULL;                   /\* the semi-coherent statistic results *\/    */
  REAL8 fmin_read,fmax_read,fband_read;         /* the range of frequencies to be read from SFTs */
  UINT4 i;                                      /* counters */
  FILE *sfp = NULL;
  /* FILE *cfp = NULL; */

  vrbflg = 0;	                        /* verbose error-messages */

  /* turn off default GSL error handler */
  gsl_set_error_handler_off();

  /* register and read all user-variables */
  if (XLALReadUserVars(argc,argv,&uvar,&clargs)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadUserVars() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : read in uservars\n",__func__);

  /* initialise the random number generator */
  if (XLALInitgslrand(&r,uvar.seed)) {
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

  /* make crude but safe estimate of the bandwidth required for the source - now includes running median wings */
  {
    REAL8 wings = LAL_TWOPI*uvar.maxasini/uvar.minorbperiod;
    fmin_read = MINBAND*floor((uvar.freq - WINGS_FACTOR*uvar.freq*wings - (REAL8)uvar.blocksize/(REAL8)uvar.tsft)/MINBAND);
    fmax_read = MINBAND*ceil((uvar.freq + (REAL8)uvar.blocksize/(REAL8)uvar.tsft + uvar.freqband + WINGS_FACTOR*(uvar.freq + uvar.freqband)*wings)/MINBAND);
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
  if (XLALReadSFTs(&sftvec,uvar.sftbasename,fmin_read,fband_read,uvar.gpsstart,uvar.gpsend,uvar.tsft)) {
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
  /* NORMALISE THE SFTS */
  /**********************************************************************************/

  if (uvar.blocksize>0) {
    /* compute the background noise using the sfts - this routine uses the running median at the edges to normalise the wings */
    if (XLALNormalizeSFTVect(sftvec,uvar.blocksize,0)) {
      LogPrintf(LOG_CRITICAL,"%s : XLALNormaliseSFTVect() failed with error = %d\n",__func__,xlalErrno);
      return 1;
    }
  }
  else {
    REAL8Vector *means = XLALCreateREAL8Vector(sftvec->length);
    if (uvar.blocksize==0) {
      /* compute the background noise using the sfts - this routine simply divides by the median */
      if (XLALNormalizeSFTVectMedian(sftvec,means,1)) {
        LogPrintf(LOG_CRITICAL,"%s : XLALNormaliseSFTVectMedian() failed with error = %d\n",__func__,xlalErrno);
        return 1;
      }
    }
    else {
      /* compute the background noise using the sfts - this routine simply divides by the mean */
      if (XLALNormalizeSFTVectMean(sftvec,means,1)) {
        LogPrintf(LOG_CRITICAL,"%s : XLALNormaliseSFTVectMean() failed with error = %d\n",__func__,xlalErrno);
        return 1;
      }
    }
    XLALDestroyREAL8Vector(means);
  }
  LogPrintf(LOG_DEBUG,"%s : normalised the SFTs\n",__func__);

  /* for (i=0;i<sftvec->length;i++) {
    for (j=0;j<sftvec->data[i].data->length;j++) {
      if (isnan(crealf(sftvec->data[i].data->data[j]))||isinf(crealf(sftvec->data[i].data->data[j]))||isnan(cimagf(sftvec->data[i].data->data[j]))||isinf(cimagf(sftvec->data[i].data->data[j]))) {
        fprintf(stdout,"SFT %d : %f %e %e\n",i,sftvec->data[i].f0 + j*sftvec->data[i].deltaF,crealf(sftvec->data[i].data->data[j]),cimagf(sftvec->data[i].data->data[j]));
      }
    }
  } */

  /**********************************************************************************/
  /* DEFINE THE BINARY PARAMETER SPACE */
  /**********************************************************************************/

  /* define the binary parameter space */
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

  /* if (XLALOpenSemiCoherentResultsFile(&cfp,newnewtemp,&pspace,clargs,&uvar,1)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALOpenCoherentResultsFile() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : opened coherent results file.\n",__func__); */

   /**********************************************************************************/
  /* COMPUTE THE STATISTICS ON THE COARSE GRID */
  /**********************************************************************************/

  /* compute the demodulated power on the frequency derivitive grid */
  if (XLALCOMPLEX8TimeSeriesArrayToDemodPowerVector(&dmpower,dstimevec,freqgridparams,NULL)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCOMPLEX8TimeSeriesArrayToDemodPowerVector() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  /* fclose(cfp); */
  LogPrintf(LOG_DEBUG,"%s : computed the demodulated power\n",__func__);

  /**********************************************************************************/
  /* OPEN RESULTS FILE */
  /**********************************************************************************/

  if (XLALOpenSemiCoherentResultsFile(&sfp,newnewtemp,&pspace,clargs,&uvar)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALOutputBayesResults() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : output results to file.\n",__func__);

   /**********************************************************************************/
  /* COMPUTE THE STATISTICS ON THE FINE GRID */
  /**********************************************************************************/

  /* compute the semi-coherent detection statistic on the fine grid */
  if (XLALComputeSemiCoherentStat(sfp,dmpower,&pspace,freqgridparams,bingridparams,uvar.ntoplist)) {
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
      LogPrintf(LOG_CRITICAL,"%s : unable to move final results directory %s -> %s.  Exiting.\n",__func__,newnewtemp,newoutputdir);
      return 1;
    }
  }
  LogPrintf(LOG_DEBUG,"%s : moved the results from the temp space.\n",__func__);

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
  uvar->ntoplist = 10;
  uvar->coverage = -1;
  uvar->blocksize = 100;
  uvar->tsft = 256;
  uvar->seed = 1;
  uvar->tempdir = NULL;

  /* initialise all parameter space ranges to zero */
  uvar->freqband = 0;
  uvar->minorbperiod = 0.0;
  uvar->maxorbperiod = 0.0;
  uvar->minasini = 0.0;
  uvar->maxasini = 0.0;
  uvar->tasc = -1.0;
  uvar->deltaorbphase = 2.0*LAL_PI;

  /* ---------- register all user-variables ---------- */
  XLALRegisterUvarMember(sftbasename, 	        STRING, 'i', REQUIRED, "The basename of the input SFT files");
  XLALRegisterUvarMember(outputdir, 	        STRING, 'o', REQUIRED, "The output directory name");
  XLALRegisterUvarMember(comment, 	        STRING, 'C', REQUIRED, "An analysis descriptor string");
  XLALRegisterUvarMember(tempdir,              STRING, 'z', OPTIONAL, "A temporary directory");
  XLALRegisterUvarMember(freq,                   REAL8, 'f', REQUIRED, "The starting frequency (Hz)");
  XLALRegisterUvarMember(freqband,   	        REAL8, 'b', OPTIONAL, "The frequency band (Hz)");
  XLALRegisterUvarMember(minorbperiod,           REAL8, 'p', REQUIRED, "The minimum orbital period value (sec)");
  XLALRegisterUvarMember(maxorbperiod,   	REAL8, 'P', OPTIONAL, "The maximum orbital period value (sec)");
  XLALRegisterUvarMember(minasini,               REAL8, 'a', REQUIRED, "The minimum orbital semi-major axis (sec)");
  XLALRegisterUvarMember(maxasini,       	REAL8, 'A', OPTIONAL, "The maximum orbital semi-major axis (sec)");
  XLALRegisterUvarMember(tasc,                   REAL8, 't', REQUIRED, "The best guess orbital time of ascension (rads)");
  XLALRegisterUvarMember(deltaorbphase,      	REAL8, 'T', OPTIONAL, "The orbital phase uncertainty (cycles)");
  XLALRegisterUvarMember(mismatch,        	REAL8, 'm', OPTIONAL, "The grid mismatch (0->1)");
  XLALRegisterUvarMember(coverage,        	REAL8, 'c', OPTIONAL, "The random template coverage (0->1)");
  XLALRegisterUvarMember(blocksize,        	INT4, 'r', OPTIONAL, "The running median block size");
  XLALRegisterUvarMember(tsft,                    INT4, 'S', OPTIONAL, "The length of the input SFTs in seconds");
  XLALRegisterUvarMember(ntoplist,                INT4, 'x', OPTIONAL, "output the top N results");
  XLALRegisterUvarMember(seed,                    INT4, 'X', OPTIONAL, "The random number seed (0 = clock)");
  XLALRegisterUvarMember(gpsstart,                INT4, 's', OPTIONAL, "The minimum start time (GPS sec)");
  XLALRegisterUvarMember(gpsend,          	INT4, 'e', OPTIONAL, "The maximum end time (GPS sec)");
  XLALRegisterUvarMember(version,                BOOLEAN, 'V', SPECIAL,  "Output code version");

  /* do ALL cmdline and cfgfile handling */
  BOOLEAN should_exit = 0;
  if (XLALUserVarReadAllInput(&should_exit, argc, argv)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALUserVarReadAllInput failed with error = %d\n",__func__,xlalErrno);
    return XLAL_EFAULT;
  }
  if (should_exit) exit(1);

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
				GridParametersVector *fgrid,		/**< UNDOCUMENTED */
				GridParameters *bingrid,                 /**< [in] the grid parameters */
				INT4 ntoplist		/**< UNDOCUMENTED */
				)
{

  toplist TL;					      /* the results toplist */
  Template *bintemp = NULL;                           /* the binary parameter space template */
  Template fdots;                                     /* the freq derivitive template for each segment */
  UINT4 i,j;                                          /* counters */
  UINT4 percent = 0;                                  /* counter for status update */
  REAL8 mean = 0.0;
  REAL8 var = 0.0;

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

  /* allocate memory for the results toplist */
  TL.n = ntoplist;
  if ((TL.data = (REAL8 *)XLALCalloc(TL.n,sizeof(REAL8))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if ((TL.params = (REAL8 **)XLALCalloc(TL.n,sizeof(REAL8 *))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  for (i=0;i<(UINT4)TL.n;i++) {
    TL.params[i] = NULL;
    if ((TL.params[i] = (REAL8 *)XLALCalloc(bingrid->ndim,sizeof(REAL8))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_ENOMEM);
    }
  }
  TL.idx = 0;

  /* allocate memory for the fdots */
  if ((fdots.x = XLALCalloc(fgrid->segment[0]->ndim,sizeof(REAL8))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : XLALCalloc() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  fdots.ndim = fgrid->segment[0]->ndim;

  /* define the chi-squared threshold */
  /* REAL8 thr = gsl_cdf_chisq_Qinv(frac,2*power->length);
  LogPrintf(LOG_DEBUG,"%s : computed the threshold as %f\n",__func__,thr); */

  int (*getnext)(Template **temp,GridParameters *gridparams, ParameterSpace *space,void *);
  INT4 newmax = bingrid->max;
  ParameterSpace *temppspace = NULL;
  if (bingrid->Nr>0) {
    getnext = &XLALGetNextRandomBinaryTemplate;
    newmax = bingrid->Nr;
    temppspace = pspace;
  }
  else getnext = &XLALGetNextTemplate;

  /* single loop over binary templates */
  while (getnext(&bintemp,bingrid,temppspace,r)) {

    REAL8 logLratiosum = 0.0;                       /* initialise likelihood ratio */

    /** loop over segments **********************************************************************************/
    for (i=0;i<power->length;i++) {

      REAL4DemodulatedPower *currentpower = power->segment[i];
      GridParameters *fdotgrid = fgrid->segment[i];
      REAL8 tmid = XLALGPSGetREAL8(&(power->segment[i]->epoch)) + 0.5*pspace->tseg;
      /* REAL8 norm = (REAL8)background->data[i]; */
      INT4 idx = 0;

      /* compute instantaneous frequency derivitives corresponding to the current template for this segment */
      XLALComputeBinaryFreqDerivitives(&fdots,bintemp,tmid);

      /* find indices corresponding to the spin derivitive values for the segment power */
      for (j=0;j<fdots.ndim;j++) {
	UINT4 tempidx = 0.5 + (fdots.x[j] - fdotgrid->grid[j].min)*fdotgrid->grid[j].oneoverdelta;
	idx += tempidx*fdotgrid->prod[j];
      }

      /* define the power at this location in this segment */
      logLratiosum += currentpower->data->data[idx]; /* /norm; */

    } /* end loop over segments */
    /*************************************************************************************/

    /* output semi-coherent statistic for this template if it exceeds the threshold*/
    logLratiosum *= 2.0;     /* make it a true chi-squared variable */
    mean += logLratiosum;
    var += logLratiosum*logLratiosum;
    XLALtoplist(logLratiosum,bintemp,&TL);

    /* if (logLratiosum>thr) {
      for (j=0;j<bingrid->ndim;j++) fprintf(fp,"%6.12f\t",bintemp->x[j]);
      fprintf(fp,"%6.12e\n",logLratiosum);
    } */

    /* output status to screen */
    if (floor(100.0*(REAL8)bintemp->currentidx/(REAL8)newmax) > (REAL8)percent) {
      percent = (UINT4)floor(100*(REAL8)bintemp->currentidx/(REAL8)newmax);
      LogPrintf(LOG_DEBUG,"%s : completed %d%% (%d/%d)\n",__func__,percent,bintemp->currentidx,newmax);
    }

  } /* end loop over templates */
  /*************************************************************************************/

  /* compute mean and variance of results */
  mean = mean/(REAL8)newmax;
  var = var/(REAL8)(newmax-1);
  var = (var - mean*mean);

  /* modify toplist to have correct mean and variance */
  for (i=0;i<(UINT4)TL.n;i++) TL.data[i] = (TL.data[i] - mean)*sqrt(4.0*power->length)/sqrt(var) + 2.0*power->length;

  /* output toplist to file */
  for (i=0;i<(UINT4)TL.n;i++) {
    for (j=0;j<bingrid->ndim;j++) fprintf(fp,"%6.12f\t",TL.params[i][j]);
    fprintf(fp,"%6.12e\n",TL.data[i]);
  }

  /* free template memory */
  XLALFree(fdots.x);
  XLALFree(TL.data);
  for (i=0;i<(UINT4)TL.n;i++) XLALFree(TL.params[i]);
  XLALFree(TL.params);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}

int XLALtoplist(REAL8 x,		/**< [in] the data to add to the toplist */
		Template *params,       /**< [in] the parameters for this result */
		toplist *TL	        /**< [in/out] the toplist */
		)
{

  INT4 i;

  /* check if toplist is allocated */
  if (TL == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Unable to allocate memory for the toplist.  Exiting.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (TL->idx>=TL->n) {
    LogPrintf(LOG_CRITICAL,"%s: toplist index exceeds size of toplist.  Exiting.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* check if the input data is larger than the lowest currently in the list */
  if (x<TL->data[TL->n-1]) return 0;

  /* place the input data into the current location */
  TL->data[TL->idx] = x;
  for (i=0;i<(INT4)params->ndim;i++) TL->params[TL->idx][i] = params->x[i];
  if (TL->idx<TL->n-1) TL->idx++;

  /* find current lowest */
  REAL8 lowest = TL->data[TL->n-1];
  INT4 li = TL->n-1;
  for (i=0;i<TL->n;i++) {
    if (TL->data[i]<lowest) {
      lowest = TL->data[i];
      li = i;
    }
  }

  /* switch current lowest to the end location */
  for (i=0;i<(INT4)params->ndim;i++) {
    REAL8 oldparam = TL->params[TL->n-1][i];
    TL->params[TL->n-1][i] = TL->params[li][i];
    TL->params[li][i] = oldparam;
  }
  REAL8 old = TL->data[TL->n-1];
  TL->data[TL->n-1] = lowest;
  TL->data[li] = old;

  return 0;

}


/** Output the results to file
 *
 * We choose to output all results from a specific analysis to a single file
 *
 */
int XLALOpenSemiCoherentResultsFile(FILE **fp,                  /**< [in] filepointer to output file */
				    CHAR *outputdir,            /**< [in] the output directory name */
				    ParameterSpace *pspace,     /**< [in] the parameter space */
				    CHAR *clargs,               /**< [in] the command line args */
				    UserInput_t *uvar		/**< UNDOCUMENTED */
				    )
{
  CHAR outputfile[LONGSTRINGLENGTH];    /* the output filename */
  time_t curtime = time(NULL);          /* get the current time */
  CHAR *time_string = NULL;             /* stores the current time */
  CHAR *version_string = NULL;          /* pointer to a string containing the git version information */

  /* validate input */
  if (outputdir == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output directory string == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  if (pspace == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, ParameterSpace structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* define the output filename */
  /* the format we adopt is the following SemiCoherentResults-<SOURCE>-<START>_<END>-<MIN_FREQ_INT>_<MIN_FREQ_mHZ>_ <MAX_FREQ_INT>_<MAX_FREQ_mHZ>.txt */
  {
    UINT4 min_freq_int = floor(pspace->space->data[0].min);
    UINT4 max_freq_int = floor(pspace->space->data[0].max);
    UINT4 min_freq_mhz = (UINT4)floor(0.5 + (pspace->space->data[0].min - (REAL8)min_freq_int)*1e3);
    UINT4 max_freq_mhz = (UINT4)floor(0.5 + (pspace->space->data[0].max - (REAL8)max_freq_int)*1e3);
    UINT4 end = (UINT4)ceil(XLALGPSGetREAL8(&(pspace->epoch)) + pspace->span);
    /* if (coherent) snprintf(outputfile,LONGSTRINGLENGTH,"%s/CoherentResults-%s-%d_%d-%04d_%03d_%04d_%03d.txt",
			   outputdir,(CHAR*)uvar->comment,pspace->epoch.gpsSeconds,end,min_freq_int,min_freq_mhz,max_freq_int,max_freq_mhz);
    else */
    snprintf(outputfile,LONGSTRINGLENGTH,"%s/SemiCoherentResults-%s-%d_%d-%04d_%03d_%04d_%03d.txt",
		  outputdir,(CHAR*)uvar->comment,pspace->epoch.gpsSeconds,end,min_freq_int,min_freq_mhz,max_freq_int,max_freq_mhz);
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

  /* free memory */
  XLALFree(time_string);
  XLALFree(version_string);

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
