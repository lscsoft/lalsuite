/*  Copyright (C) 2013 Chris Messenger
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
 * This code is designed to search for intermittent signals in a time series
 *
 */

/***********************************************************************************************/
/* includes */
#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <lal/TimeSeries.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lal/SFTutils.h>
#include <lal/SFTfileIO.h>
#include <lal/ComplexFFT.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lalapps.h>
#include <lal/BandPassTimeSeries.h>

#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/LALConstants.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/RealFFT.h>

#include <lal/TimeSeries.h>
#include <lal/GeneratePulsarSignal.h>
#include <lal/TransientCW_utils.h>
#include <lal/LALString.h>
#include <lal/StringVector.h>

#include <gsl/gsl_rng.h>           /* for random number generation */ 
#include <gsl/gsl_randist.h>       /* for random number generation */ 

#include "SemiCoherent.h"

/** A structure that stores user input variables 
 */
typedef struct { 
  BOOLEAN help;		            /**< trigger output of help string */
  CHAR *outLabel;                   /**< 'misc' entry in SFT-filenames or 'description' entry of frame filenames */
  CHAR *outputdir;                  /**< the output directory */
  CHAR *cachefile;                  /**< the name of the input cache file */
  REAL8 fsamp;                      /**< the sampling frequency of the data */
  INT4 Tmin;                        /**< the minimum length of the SFTs */
  INT4 Tmax;                        /**< the maximum length of the SFTs */
  REAL8 freq;                       /**< the starting frequency */
  REAL8 freqband;                   /**< the band width */
  REAL8 highpassf;                  /**< the high pass filter frequency */
  REAL8 minorbperiod;               /**< the minimum orbital period value */
  REAL8 maxorbperiod;               /**< the maximum orbital period value*/
  REAL8 minasini;                   /**< the minimum orbital semi-major axis value */
  REAL8 maxasini;                   /**< the maximum orbital semi-major axis value */
  REAL8 tasc;                       /**< the best guess orbital time of ascension */
  REAL8 deltaorbphase;              /**< the orbital phase uncertainty (cycles) */
  REAL8 mismatch;                   /**< the grid mismatch */      
  REAL8 thresh;                     /**< the number of results to output */
  CHAR *tempdir;                    /**< a temporary directory for keeping the results */
  BOOLEAN version;	            /**< output version-info */
} UserInput_t;

/***********************************************************************************************/
/* define functions */
int main(int argc,char *argv[]);
int XLALReadUserVars(int argc,char *argv[],UserInput_t *uvar);
int compare_function(const void *a,const void *b);
int XLALInitgslrand(gsl_rng **gslrnd,INT8 seed);
int XLALCopySFT (SFTtype *dest, const SFTtype *src);
int XLALAppendSFT2Vector (SFTVector *vect,const SFTtype *sft);
int XLALDefineBinaryParameterSpace(REAL8Space**,LIGOTimeGPS,REAL8,void* );
int XLALOpenIntermittentResultsFile(FILE **,CHAR *,CHAR *,UserInput_t *uvar,LIGOTimeGPS *start,LIGOTimeGPS *end);	  

/***********************************************************************************************/
/* global variables */
extern int vrbflg;	 	/**< defined in lalapps.c */

/***********************************************************************************************/
/* empty initializers */
UserInput_t empty_UserInput;
static const LALUnit empty_LALUnit;
ParameterSpace empty_ParameterSpace;

/** The main function of Intermittent.c
 *
 */
int main( int argc, char *argv[] )  {

  UserInput_t uvar = empty_UserInput;           /* user input variables */
  INT4 i,j,k;                                     /* counter */
  REAL8 fmin_read,fmax_read,fband_read;         /* the range of frequencies to be read from SFTs */
  CHAR newnewtemp[LONGSTRINGLENGTH];
  FILE *ifp = NULL;
  CHAR *clargs = NULL;                          /* store the command line args */

  vrbflg = 1;	                        /* verbose error-messages */

  /**********************************************************************************/
  /* register and read all user-variables */
  if (XLALReadUserVars(argc,argv,&uvar)) {
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

  /**********************************************************************************/
  /* read in the cache file */
  FILE *cachefp = NULL;
  if ((cachefp = fopen(uvar.cachefile,"r")) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to open binary input file %s\n",__func__,uvar.cachefile);
    return 1;
  }
  i = 0;
  while (fscanf(cachefp,"%*s %*d %*d")!=EOF) i++;
  INT4 Nfiles = i;
  fclose(cachefp);
  LogPrintf(LOG_DEBUG,"%s : counted %d files listed in the cache file.\n",__func__,Nfiles);

  /* allocate memory */
  char **filenames = LALCalloc(Nfiles,sizeof(char*));
  LIGOTimeGPSVector fileStart;
  fileStart.data = LALCalloc(Nfiles,sizeof(LIGOTimeGPS));
  for (i=0;i<Nfiles;i++) filenames[i] = LALCalloc(512,sizeof(char));
  
  if ((cachefp = fopen(uvar.cachefile,"r")) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to open binary input file %s\n",__func__,uvar.cachefile);
    return 1;
  }

  for (i=0;i<Nfiles;i++) {
    fscanf(cachefp,"%s %d %d",filenames[i],&(fileStart.data[i].gpsSeconds),&(fileStart.data[i].gpsNanoSeconds));
  }
  fclose(cachefp);

  /* initialise results vector */
  SFTVector *SFTvect = LALCalloc(1,sizeof(SFTVector));

  /**********************************************************************************/
  /* generate comment string */
  char *VCSInfoString = XLALGetVersionString(0);
  XLAL_CHECK ( VCSInfoString != NULL, XLAL_EFUNC, "XLALGetVersionString(0) failed.\n" );
  CHAR *logstr;
  size_t len;
  XLAL_CHECK ( (logstr = XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE )) != NULL, XLAL_EFUNC );
  char *comment = XLALCalloc ( 1, len = strlen ( logstr ) + strlen(VCSInfoString) + 512 );
  XLAL_CHECK ( comment != NULL, XLAL_ENOMEM, "XLALCalloc(1,%d) failed.\n", len );
  sprintf ( comment, "Generated by:\n%s\n%s\n", logstr, VCSInfoString );

  /* setup the gridding over coherent times */
  INT4 dummyT = uvar.Tmin;
  INT4 NT = 0;
  while (dummyT<uvar.Tmax) {
    dummyT = (INT4)((REAL8)dummyT*(1.0 + uvar.mismatch));
    NT++;
  }
  LogPrintf(LOG_DEBUG,"%s : Going to do %d different coherent lengths.\n",__func__,NT);
 
  /**********************************************************************************/
  /* OPEN INTERMEDIATE RESULTS FILE */
  /**********************************************************************************/
  
  if (XLALOpenIntermittentResultsFile(&ifp,newnewtemp,clargs,&uvar,&(fileStart.data[0]),&(fileStart.data[Nfiles-1]))) {
    LogPrintf(LOG_CRITICAL,"%s : XLALOpenCoherentResultsFile() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : opened coherent results file.\n",__func__);

  /* loop over the different coherent lengths */
  INT4 currentT = uvar.Tmin;
  for (i=0;i<NT;i++) {

    /* define the number of different start times */
    INT4 Ns = ceil(1.0/uvar.mismatch);
    REAL8 Tstep = currentT/(REAL8)Ns;
    
    /* loop over different start times */
    REAL8 currentoffset = 0.0;
    for (k=0;k<Ns;k++) {
      
      SFTVector *sftvec = NULL;                     /* stores the input SFTs */

      /**********************************************************************************/
      /* loop over the input files */
      for (j=0;j<Nfiles;j++) {
		
	REAL8Vector *background = NULL;                  /* estimates of the background for each SFT */
	ParameterSpace pspace = empty_ParameterSpace;    /* the search parameter space */
	COMPLEX8TimeSeriesArray *dstimevec = NULL;       /* contains the downsampled inverse FFT'd SFTs */
	REAL4DemodulatedPowerVector *dmpower = NULL;    /* contains the demodulated power for all SFTs */
	GridParametersVector *freqgridparams = NULL;  /* the coherent grid on the frequency derivitive parameter space */
	GridParameters *bingridparams = NULL;

	/* offset start time */
	LIGOTimeGPS *currentStart = XLALGPSAdd(&(fileStart.data[j]),currentoffset);
	
	/**********************************************************************************/
	/* GENERATE SFTS FROM BINARY INPUT FILE */
	/**********************************************************************************/
	
	/* converts a binary input file into sfts */
	if (XLALBinaryToSFTVector(sftvec,filenames[j],currentStart,currentT,fmin_read,fmax_read,uvar.fsamp,uvar.highpassf)) {
	  LogPrintf(LOG_CRITICAL,"%s : XLALEstimateBackgroundFlux() failed with error = %d\n",__func__,xlalErrno);
	  return 1;
	}
	LogPrintf(LOG_DEBUG,"%s : estimated the background noise from the SFTs\n",__func__);
	exit(0);
/* 	/\* define SFT length and the start and span of the observations plus the definitive segment time *\/ */
/* 	pspace.tseg = 1.0/sftvec->data[0].deltaF; */
/* 	memcpy(&(pspace.epoch),&(sftvec->data[0].epoch),sizeof(LIGOTimeGPS)); */
/* 	pspace.span = XLALGPSDiff(&(sftvec->data[sftvec->length-1].epoch),&(sftvec->data[0].epoch)) + pspace.tseg; */
/* 	LogPrintf(LOG_DEBUG,"%s : SFT length = %f seconds\n",__func__,pspace.tseg); */
/* 	LogPrintf(LOG_DEBUG,"%s : entire dataset starts at GPS time %d contains %d SFTS and spans %.0f seconds\n",__func__,pspace.epoch.gpsSeconds,sftvec->length,pspace.span); */
	
/* 	/\**********************************************************************************\/ */
/* 	/\* DEFINE THE BINARY PARAMETER SPACE *\/ */
/* 	/\**********************************************************************************\/ */
	
/* 	/\*define the binary parameter space *\/ */
/* 	if (XLALDefineBinaryParameterSpace(&(pspace.space),pspace.epoch,pspace.span,&uvar)) { */
/* 	  LogPrintf(LOG_CRITICAL,"%s : XLALDefineBinaryParameterSpace() failed with error = %d\n",__func__,xlalErrno); */
/* 	  return 1; */
/* 	} */
/* 	LogPrintf(LOG_DEBUG,"%s : defined binary parameter prior space\n",__func__); */
	
/* 	/\**********************************************************************************\/ */
/* 	/\* ESTIMATE BACKGROUND NOISE FROM SFTS *\/ */
/* 	/\**********************************************************************************\/ */
	
/* 	/\* compute the background noise using the sfts *\/ */
/* 	if (XLALEstimateBackgroundFlux(&background,sftvec)) { */
/* 	  LogPrintf(LOG_CRITICAL,"%s : XLALEstimateBackgroundFlux() failed with error = %d\n",__func__,xlalErrno); */
/* 	  return 1; */
/* 	} */
/* 	LogPrintf(LOG_DEBUG,"%s : estimated the background noise from the SFTs\n",__func__);  */
	
/* 	/\**********************************************************************************\/ */
/* 	/\* COMPUTE THE COARSE GRID ON FREQUENCY DERIVITIVES *\/ */
/* 	/\**********************************************************************************\/ */
	
/* 	/\* compute the grid parameters for all SFTs *\/ */
/* 	if (XLALComputeFreqGridParamsVector(&freqgridparams,pspace.space,sftvec,uvar.mismatch)) { */
/* 	  LogPrintf(LOG_CRITICAL,"%s : XLALComputeFreqGridParams() failed with error = %d\n",__func__,xlalErrno); */
/* 	  return 1; */
/* 	} */

/* 	/\**********************************************************************************\/ */
/* 	/\* CONVERT ALL SFTS TO DOWNSAMPLED TIMESERIES *\/ */
/* 	/\**********************************************************************************\/ */
	
/* 	/\* convert sfts to downsample dtimeseries *\/ */
/* 	if (XLALSFTVectorToCOMPLEX8TimeSeriesArray(&dstimevec,sftvec)) { */
/* 	  LogPrintf(LOG_CRITICAL,"%s : XLALSFTVectorToCOMPLEX8TimeSeriesArray() failed with error = %d\n",__func__,xlalErrno); */
/* 	  return 1; */
/* 	} */
/* 	LogPrintf(LOG_DEBUG,"%s : converted SFTs to downsampled timeseries\n",__func__); */
	
/* 	/\**********************************************************************************\/ */
/* 	/\* COMPUTE THE STATISTICS ON THE COARSE GRID *\/ */
/* 	/\**********************************************************************************\/ */
	
/* 	/\* compute the demodulated power on the frequency derivitive grid *\/ */
/* 	if (XLALCOMPLEX8TimeSeriesArrayToDemodPowerVector(&dmpower,dstimevec,freqgridparams,background,ifp)) { */
/* 	  LogPrintf(LOG_CRITICAL,"%s : XLALCOMPLEX8TimeSeriesArrayToDemodPowerVector() failed with error = %d\n",__func__,xlalErrno); */
/* 	  return 1; */
/* 	} */
/* 	fclose(cfp); */
/* 	LogPrintf(LOG_DEBUG,"%s : computed the demodulated power\n",__func__); */

/* 	/\**********************************************************************************\/ */
/* 	/\* FREE MEMORY *\/ */
/* 	/\**********************************************************************************\/ */

/* 	/\* free memory inside the loop *\/ */
/* 	XLALDestroyREAL4TimeSeries (Tseries);	 */
/* 	XLALDestroySFTVector ( SFTvect ); */

/* 	/\* clean up the parameter space *\/ */
/* 	if (XLALFreeParameterSpace(&pspace)) { */
/* 	  LogPrintf(LOG_CRITICAL,"%s : XLALFreeParameterSpace() failed with error = %d\n",__func__,xlalErrno); */
/* 	  return 1; */
/* 	} */
/* 	LogPrintf(LOG_DEBUG,"%s : freed the parameter space\n",__func__); */
	
/* 	/\* clean up the demodulated power *\/ */
/* 	if (XLALFreeREAL4DemodulatedPowerVector(dmpower)) { */
/* 	  LogPrintf(LOG_CRITICAL,"%s : XLALFreeREAL4DemodulatedPowerVector() failed with error = %d\n",__func__,xlalErrno); */
/* 	  return 1; */
/* 	} */
/* 	LogPrintf(LOG_DEBUG,"%s : freed the demodulated power\n",__func__); */
	
/* 	/\* free un-needed downsampled timeseries *\/ */
/* 	for (i=0;i<dstimevec->length;i++) { */
/* 	  XLALDestroyCOMPLEX8TimeSeries(dstimevec->data[i]); */
/* 	} */
/* 	XLALFree(dstimevec->data); */
/* 	XLALFree(dstimevec); */
/* 	LogPrintf(LOG_DEBUG,"%s : freed the downsampled timeseries memory\n",__func__); */
	
/* 	/\* free frequency grid - the contents of each segment have been moved to the power structure and are freed later *\/ */
/* 	XLALFree(freqgridparams->segment); */
/* 	XLALFree(freqgridparams); */
	
/* 	/\* free the background estimate *\/ */
/* 	XLALDestroyREAL8Vector(background); */
/* 	LogPrintf(LOG_DEBUG,"%s : freed the background noise estimate\n",__func__); */
	
      } /* end loop over start time offset */

    } /* end loop over observation time length */
    
  }  /* end loop over input files */

  /* move the temporary directory to the final location */
  if (uvar.tempdir) {
    CHAR newoutputdir[LONGSTRINGLENGTH];
    sprintf(newoutputdir,"%s/%.3f-%.3f",uvar.outputdir,uvar.freq,uvar.freq+uvar.freqband);
    if (rename(newnewtemp,newoutputdir)) {
      LogPrintf(LOG_CRITICAL,"%s : unable to move final results directory %s -> %s.  Exiting.\n",newnewtemp,newoutputdir);
      return 1;
    }
  }

  /**********************************************************************************/
  /* FREE MEMORY */
  /**********************************************************************************/
  fclose(ifp);
  XLALFree ( logstr );
  XLALFree ( comment );
  XLALFree(filenames);
  XLALFree(fileStart.data);

  LALCheckMemoryLeaks();

  return 0;

}

/*******************************************************************************/

/** Read in input user arguments
 *
 */
int XLALReadUserVars(int argc,            /**< [in] the command line argument counter */ 
		     char *argv[],        /**< [in] the command line arguments */
		     UserInput_t *uvar    /**< [out] the user input structure */
		     )
{


  /* initialise user variables */
  uvar->outLabel = NULL; 
  uvar->outputdir = NULL;
  uvar->cachefile = NULL;
  uvar->tempdir = NULL;
  uvar->fsamp = 2048;
  uvar->Tmin = 16;
  uvar->Tmax = 256;
  uvar->freq = 550;
  uvar->freqband = 1;
  uvar->highpassf = 40;
  uvar->mismatch = 0.1;
  uvar->thresh = 20;

  /* ---------- register all user-variables ---------- */
  XLALregBOOLUserStruct(help, 		        'h', UVAR_HELP,     "Print this message");
  XLALregSTRINGUserStruct(outLabel, 	        'n', UVAR_REQUIRED, "'misc' entry in SFT-filenames or 'description' entry of frame filenames"); 
  XLALregSTRINGUserStruct(outputdir, 	        'o', UVAR_REQUIRED, "The output directory name"); 
  XLALregSTRINGUserStruct(cachefile, 	        'i', UVAR_REQUIRED, "The input binary file name"); 
  XLALregSTRINGUserStruct(tempdir, 	        'x', UVAR_OPTIONAL, "The temporary directory"); 
  XLALregREALUserStruct(freq,                   'f', UVAR_REQUIRED, "The starting frequency (Hz)");
  XLALregREALUserStruct(freqband,   	        'b', UVAR_REQUIRED, "The frequency band (Hz)");
  XLALregINTUserStruct(Tmin,                    't', UVAR_OPTIONAL, "The min length of segemnts (sec)");
  XLALregINTUserStruct(Tmax,                    'T', UVAR_OPTIONAL, "The max length of segments (sec)");
  XLALregREALUserStruct(fsamp,           	's', UVAR_REQUIRED, "The sampling frequency");
  XLALregREALUserStruct(highpassf,           	'H', UVAR_OPTIONAL, "The high pass filter frequency");
  XLALregREALUserStruct(minorbperiod,          	'p', UVAR_REQUIRED, "The minimum orbital period (sec)");
  XLALregREALUserStruct(maxorbperiod,          	'P', UVAR_REQUIRED, "The maximum orbital period (sec)");
  XLALregREALUserStruct(minasini,          	'a', UVAR_REQUIRED, "The minimum projected orbital radius (sec)");
  XLALregREALUserStruct(maxasini,          	'A', UVAR_REQUIRED, "The maximum projected orbital radius (sec)");
  XLALregREALUserStruct(tasc,                   'c', UVAR_REQUIRED, "The best guess orbital time of ascension (rads)");
  XLALregREALUserStruct(deltaorbphase,      	'C', UVAR_OPTIONAL, "The orbital phase uncertainty (cycles)");
  XLALregREALUserStruct(mismatch,        	'm', UVAR_OPTIONAL, "The grid mismatch (0->1)");
  XLALregREALUserStruct(thresh,         	'z', UVAR_OPTIONAL, "The threshold on Leahy power");

  /* do ALL cmdline and cfgfile handling */
  if (XLALUserVarReadAllInput(argc, argv)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALUserVarReadAllInput() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* if help was requested, we're done here */
  if (uvar->help) exit(0);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}


/* /\*******************************************************************************\/ */
/* /\** Computes the binary parameter space boundaries given the user input args */
/*  * */
/*  * For each search dimension we define the min, max, mid, and span of that dimension */
/*  * plus we give each dimension a name. */
/*  * */
/*  *\/ */
/* int XLALDefineBinaryParameterSpace(REAL8Space **space,                 /\**< [out] the parameter space  *\/ */
/* 				   LIGOTimeGPS epoch,                  /\**< [in] the observation start epoch *\/ */
/* 				   REAL8 span,                         /\**< [in] the observation span *\/ */
/* 				   UserInput_t *uvar                          /\**< [in] the user input variables *\/ */
/* 				   ) */
/* { */
/*   REAL8 midpoint;              /\* the midpoint of the observation *\/ */
/*   REAL8 mintasc;               /\* the minimum time of ascension *\/ */
/*   REAL8 maxtasc;               /\* the maximum time of ascension *\/ */

/*   /\* validate input variables *\/ */
/*   if ((*space) != NULL) { */
/*     LogPrintf(LOG_CRITICAL,"%s : Invalid input, input REAL8Space boundary structure != NULL.\n",__func__); */
/*     XLAL_ERROR(XLAL_EINVAL); */
/*   } */
/*   if (uvar == NULL) { */
/*     LogPrintf(LOG_CRITICAL,"%s : Invalid input, input UserInput_t structure = NULL.\n",__func__); */
/*     XLAL_ERROR(XLAL_EINVAL); */
/*   } */
/*   if (epoch.gpsSeconds < 0) { */
/*     LogPrintf(LOG_CRITICAL,"%s : Invalid input, observation epoch < 0.\n",__func__); */
/*     XLAL_ERROR(XLAL_EINVAL); */
/*   } */
/*   if (span < 0) { */
/*     LogPrintf(LOG_CRITICAL,"%s : Invalid input, observation span < 0.\n",__func__); */
/*     XLAL_ERROR(XLAL_EINVAL); */
/*   } */
  
/*   /\* allocate memory for the parameter space *\/ */
/*   if ( ((*space) = XLALCalloc(1,sizeof(REAL8Space))) == NULL) { */
/*     LogPrintf(LOG_CRITICAL,"%s: XLALCalloc() failed with error = %d\n",__func__,xlalErrno); */
/*     XLAL_ERROR(XLAL_ENOMEM); */
/*   } */
/*   (*space)->ndim = NBINMAX; */
/*   if ( ((*space)->data = XLALCalloc((*space)->ndim,sizeof(REAL8Dimension))) == NULL) { */
/*     LogPrintf(LOG_CRITICAL,"%s: XLALCalloc() failed with error = %d\n",__func__,xlalErrno); */
/*     XLAL_ERROR(XLAL_ENOMEM); */
/*   } */

/*   /\* define observaton midpoint *\/ */
/*   midpoint = XLALGPSGetREAL8(&epoch) + 0.5*span; */

/*   /\* find the closest instance of ascension to the midpoint *\/ */
/*   /\* the midpoint is defined in the detector frame and the time of *\/ */
/*   /\* ascension is in the SSB frame but we only need to be roughly *\/ */
/*   /\* correct in the number of orbits we shift by *\/ */
/*   if (uvar->tasc>0) {     */
/*     REAL8 meanperiod = 0.5*(uvar->maxorbperiod + uvar->minorbperiod); */
/*     INT4 n = (INT4)floor(0.5 + (midpoint - uvar->tasc)/meanperiod); */
/*     REAL8 newtasc = uvar->tasc + n*meanperiod; */
/*     mintasc = newtasc - 0.5*uvar->maxorbperiod*uvar->deltaorbphase; */
/*     maxtasc = newtasc + 0.5*uvar->maxorbperiod*uvar->deltaorbphase; */
/*     LogPrintf(LOG_DEBUG,"%s : shifted tasc by %d orbits\n",__func__,n); */
/*   } */
/*   else {   /\* we have no orbital info so we search a full period at the midpoint *\/ */
/*     mintasc = midpoint - 0.5*uvar->maxorbperiod; */
/*     maxtasc = midpoint + 0.5*uvar->maxorbperiod; */
/*   } */
  
/*   /\* this represents a hyper-cubic parameter space *\/ */
/*   /\* we make sure that parameter ranges are consistent i.e asini > 0 etc.. *\/ */
/*   /\* frequency *\/ */
/*   snprintf((*space)->data[0].name,LALNameLength,"nu"); */
/*   (*space)->data[0].min = uvar->freq; */
/*   (*space)->data[0].max = uvar->freq + uvar->freqband; */
/*   (*space)->data[0].span = (*space)->data[0].max - (*space)->data[0].min; */

/*   /\* asini *\/ */
/*   snprintf((*space)->data[1].name,LALNameLength,"asini"); */
/*   (*space)->data[1].min = uvar->minasini; */
/*   (*space)->data[1].max = uvar->maxasini; */
/*   (*space)->data[1].span = (*space)->data[1].max - (*space)->data[1].min; */
  
/*   /\* orbphase *\/ */
/*   snprintf((*space)->data[2].name,LALNameLength,"tasc"); */
/*   (*space)->data[2].min = mintasc; */
/*   (*space)->data[2].max = maxtasc; */
/*   (*space)->data[2].span = (*space)->data[2].max - (*space)->data[2].min; */
  
/*   /\* omega *\/ */
/*   snprintf((*space)->data[3].name,LALNameLength,"omega"); */
/*   (*space)->data[3].min = LAL_TWOPI/uvar->maxorbperiod; */
/*   (*space)->data[3].max = LAL_TWOPI/uvar->minorbperiod; */
/*   (*space)->data[3].span = (*space)->data[3].max - (*space)->data[3].min; */

/*   /\* output boundaries to screen *\/ */
/*   LogPrintf(LOG_DEBUG,"%s : using flat priors on the following ranges\n",__func__); */
/*   LogPrintf(LOG_DEBUG,"%s : parameter space, %s = [%e -> %e]\n",__func__,(*space)->data[0].name,(*space)->data[0].min,(*space)->data[0].max); */
/*   LogPrintf(LOG_DEBUG,"%s : parameter space, %s = [%e -> %e]\n",__func__,(*space)->data[1].name,(*space)->data[1].min,(*space)->data[1].max); */
/*   LogPrintf(LOG_DEBUG,"%s : parameter space, %s = [%e -> %e]\n",__func__,(*space)->data[2].name,(*space)->data[2].min,(*space)->data[2].max); */
/*   LogPrintf(LOG_DEBUG,"%s : parameter space, %s = [%e -> %e]\n",__func__,(*space)->data[3].name,(*space)->data[3].min,(*space)->data[3].max); */
 
/*   LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__); */
/*   return XLAL_SUCCESS; */
  
/* } */


/* /\** Output the results to file */
/*  * */
/*  * We choose to output all results from a specific analysis to a single file */
/*  * */
/*  *\/ */
/* int XLALOpenIntermittentResultsFile(FILE **fp,                  /\**< [in] filepointer to output file *\/ */
/* 				    CHAR *outputdir,            /\**< [in] the output directory name *\/ */
/* 				    ParameterSpace *pspace,     /\**< [in] the parameter space *\/ */
/* 				    CHAR *clargs,               /\**< [in] the command line args *\/ */
/* 				    UserInput_t *uvar, */
/* 				    INT4 coherent               /\**< [in] flag for outputting coherent results *\/ */
/* 				    ) */
/* { */
/*   CHAR outputfile[LONGSTRINGLENGTH];    /\* the output filename *\/ */
/*   time_t curtime = time(NULL);          /\* get the current time *\/ */
/*   CHAR *time_string = NULL;             /\* stores the current time *\/ */
/*   CHAR *version_string = NULL;          /\* pointer to a string containing the git version information *\/ */

/*   /\* validate input *\/ */
/*   if (outputdir == NULL) { */
/*     LogPrintf(LOG_CRITICAL,"%s: Invalid input, output directory string == NULL.\n",__func__); */
/*     XLAL_ERROR(XLAL_EINVAL); */
/*   } */

/*   if (pspace == NULL) { */
/*     LogPrintf(LOG_CRITICAL,"%s: Invalid input, ParameterSpace structure == NULL.\n",__func__); */
/*     XLAL_ERROR(XLAL_EINVAL); */
/*   } */
  
/*   /\* define the output filename *\/ */
/*   /\* the format we adopt is the following SemiCoherentResults-<SOURCE>-<START>_<END>-<MIN_FREQ_INT>_<MIN_FREQ_mHZ>_ <MAX_FREQ_INT>_<MAX_FREQ_mHZ>.txt *\/ */
/*   { */
/*     UINT4 min_freq_int = floor(pspace->space->data[0].min); */
/*     UINT4 max_freq_int = floor(pspace->space->data[0].max); */
/*     UINT4 min_freq_mhz = (UINT4)floor(0.5 + (pspace->space->data[0].min - (REAL8)min_freq_int)*1e3); */
/*     UINT4 max_freq_mhz = (UINT4)floor(0.5 + (pspace->space->data[0].max - (REAL8)max_freq_int)*1e3); */
/*     UINT4 end = (UINT4)ceil(XLALGPSGetREAL8(&(pspace->epoch)) + pspace->span); */
/*     if (coherent) snprintf(outputfile,LONGSTRINGLENGTH,"%s/CoherentResults-%s-%d_%d-%04d_%03d_%04d_%03d.txt", */
/* 			   outputdir,(CHAR*)uvar->comment,pspace->epoch.gpsSeconds,end,min_freq_int,min_freq_mhz,max_freq_int,max_freq_mhz); */
/*     else snprintf(outputfile,LONGSTRINGLENGTH,"%s/SemiCoherentResults-%s-%d_%d-%04d_%03d_%04d_%03d.txt", */
/* 		  outputdir,(CHAR*)uvar->comment,pspace->epoch.gpsSeconds,end,min_freq_int,min_freq_mhz,max_freq_int,max_freq_mhz); */
/*   } */
/*   LogPrintf(LOG_DEBUG,"%s : output %s\n",__func__,outputfile); */

/*   /\* open the output file *\/ */
/*   if (((*fp) = fopen(outputfile,"w")) == NULL) { */
/*     LogPrintf(LOG_CRITICAL,"%s: Error, failed to open file %s for writing.  Exiting.\n",__func__,outputfile); */
/*     XLAL_ERROR(XLAL_EINVAL); */
/*   } */
 
/*   /\* Convert time to local time representation *\/ */
/*   { */
/*     struct tm *loctime = localtime(&curtime); */
/*     CHAR *temp_time = asctime(loctime); */
/*     UINT4 n = strlen(temp_time); */
/*     time_string = XLALCalloc(n,sizeof(CHAR)); */
/*     snprintf(time_string,n-1,"%s",temp_time); */
/*   } */
 
/*   /\* get GIT version information *\/ */
/*   { */
/*     CHAR *temp_version = XLALGetVersionString(0); */
/*     UINT4 n = strlen(temp_version); */
/*     version_string = XLALCalloc(n,sizeof(CHAR)); */
/*     snprintf(version_string,n-1,"%s",temp_version); */
/*     XLALFree(temp_version); */
/*   } */
  
/*   /\* output header information *\/ */
/*   fprintf((*fp),"%s \n",version_string); */
/*   fprintf((*fp),"%%%% command line args\t\t= %s\n",clargs); */
/*   fprintf((*fp),"%%%% filename\t\t\t\t= %s\n",outputfile); */
/*   fprintf((*fp),"%%%% date\t\t\t\t\t= %s\n",time_string); */
/*   fprintf((*fp),"%%%% start time (GPS sec)\t\t= %d\n",pspace->epoch.gpsSeconds); */
/*   fprintf((*fp),"%%%% observation span (sec)\t= %d\n",(UINT4)pspace->span); */
/*   fprintf((*fp),"%%%% coherent time (sec)\t\t= %d\n",(UINT4)pspace->tseg); */
/*  /\*  fprintf(fp,"%%%% number of segments\t\t= %d\n",Bayes->nsegments); *\/ */
/* /\*   fprintf(fp,"%%%% number of dimensions\t= %d\n",Bayes->gridparams->ndim); *\/ */
/* /\*   if (pspace->ampspace) fprintf(fp,"%%%% amplitude dimension\t\t\t= 1\n"); *\/ */
/* /\*   else fprintf(fp,"%%%% amplitude dimension\t\t\t= 0\n"); *\/ */
/* /\*   fprintf(fp,"%%%% mismatch\t\t\t\t= %6.12f\n",Bayes->gridparams->mismatch); *\/ */
/*   fprintf((*fp),"%%%%\n"); */

/*   /\* free memory *\/ */
/*   XLALFree(time_string); */
/*   XLALFree(version_string); */

/*   LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__); */
/*   return XLAL_SUCCESS; */

/* } */


/** Output the results to file
 *
 * We choose to output all results from a specific analysis to a single file
 *
 */
int XLALOpenIntermittentResultsFile(FILE **fp,                  /**< [in] filepointer to output file */
				    CHAR *outputdir,            /**< [in] the output directory name */
				    CHAR *clargs,               /**< [in] the command line args */
				    UserInput_t *uvar,
				    LIGOTimeGPS *start,
				    LIGOTimeGPS *end
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

  /* define the output filename */
  /* the format we adopt is the following SemiCoherentResults-<SOURCE>-<START>_<END>-<MIN_FREQ_INT>_<MIN_FREQ_mHZ>_ <MAX_FREQ_INT>_<MAX_FREQ_mHZ>.txt */
  {
    UINT4 min_freq_int = floor(uvar->freq);
    UINT4 max_freq_int = floor(uvar->freq+uvar->freqband);
    UINT4 min_freq_mhz = (UINT4)floor(0.5 + (uvar->freq - (REAL8)min_freq_int)*1e3);
    UINT4 max_freq_mhz = (UINT4)floor(0.5 + (uvar->freq + uvar->freqband - (REAL8)max_freq_int)*1e3);
    snprintf(outputfile,LONGSTRINGLENGTH,"%s/IntermittentResults-%s-%d_%d-%04d_%03d_%04d_%03d.txt",
			   outputdir,(CHAR*)uvar->outLabel,start->gpsSeconds,end->gpsSeconds,min_freq_int,min_freq_mhz,max_freq_int,max_freq_mhz);
    
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
  fprintf((*fp),"%%%% start time (GPS sec)\t\t= %d\n",start->gpsSeconds);
/*   fprintf((*fp),"%%%% observation span (sec)\t= %d\n",(UINT4)pspace->span); */
/*   fprintf((*fp),"%%%% coherent time (sec)\t\t= %d\n",(UINT4)pspace->tseg); */
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
