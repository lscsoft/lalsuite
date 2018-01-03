/* Copyright (C) 2013 Chris Messenger
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
 * \ingroup lalapps_pulsar_Xray
 * \file
 * \brief
 * This code is designed to search for intermittent signals in a time series
 *
 */

/***********************************************************************************************/
/* includes */
#include "config.h"
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
  CHAR *outLabel;                   /**< 'misc' entry in SFT-filenames or 'description' entry of frame filenames */
  CHAR *outputdir;                  /**< the output directory */
  CHAR *cachefile;                  /**< the name of the input cache file */
  CHAR *binfile;                    /**< the input binary file */
  REAL8 tsamp;                      /**< the sampling time of the data */
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
  INT4 blocksize;                   /**< the running median blocksize */
  REAL8 thresh;                     /**< the number of results to output */
  CHAR *tempdir;                    /**< a temporary directory for keeping the results */
  BOOLEAN verbose;	            /**< flag for status outputs */
} UserInput_t;

/***********************************************************************************************/
/* define functions */
int main(int argc,char *argv[]);
int XLALReadUserVars(int argc,char *argv[],UserInput_t *uvar);
int compare_function(const void *a,const void *b);
int XLALInitgslrand(gsl_rng **gslrnd,INT8 seed);
int XLALCopySFT (SFTtype *dest, const SFTtype *src);
int XLALAppendSFT2Vector (SFTVector *vect,const SFTtype *sft);
int XLALDefineBinaryParameterSpace(REAL8Space**,LIGOTimeGPS,REAL8,UserInput_t*);
int XLALOpenIntermittentResultsFile(FILE **,CHAR *,CHAR *,CHAR *,UserInput_t *uvar,LIGOTimeGPS *start);

/***********************************************************************************************/
/* global variables */
extern int vrbflg;	 	/**< defined in lalapps.c */

/***********************************************************************************************/
/* empty initializers */
UserInput_t empty_UserInput;
ParameterSpace empty_ParameterSpace;

/** The main function of Intermittent.c
 *
 */
int main( int argc, char *argv[] )  {

  UserInput_t uvar = empty_UserInput;           /* user input variables */
  INT4 i,k,m;                                     /* counter */
  CHAR newtemp[LONGSTRINGLENGTH];
  FILE *ifp = NULL;
  CHAR *clargs = NULL;                          /* store the command line args */

  /**********************************************************************************/
  /* register and read all user-variables */
  if (XLALReadUserVars(argc,argv,&uvar)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadUserVars() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  if (uvar.verbose) fprintf(stdout,"%s : read in uservars\n",__func__);

  /**********************************************************************************/
  /* make temporary directory */
  if (uvar.tempdir) {

    /* initialise the random number generator  - use the clock */
    gsl_rng * q;
    if (XLALInitgslrand(&q,0)) {
      LogPrintf(LOG_CRITICAL,"%s: XLALinitgslrand() failed with error = %d\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_EFAULT);
    }

    INT4 id = (INT4)(1e9*gsl_rng_uniform(q));
    sprintf(newtemp,"%s/%09d",uvar.tempdir,id);
    fprintf(stdout,"temp dir = %s\n",newtemp);
    if (mkdir(newtemp,0755)) {
      if (uvar.verbose) fprintf(stdout,"%s : Unable to make temporary directory %s.  Might be a problem.\n",__func__,newtemp);
      return 1;
    }
  }

  /**********************************************************************************/
  /* read in the cache file and find the correct entry for the binary file */
  FILE *cachefp = NULL;
  if ((cachefp = fopen(uvar.cachefile,"r")) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to open binary input file %s\n",__func__,uvar.cachefile);
    return 1;
  }
  i = 0;
  INT4 idx = -1;
  CHAR filename[LONGSTRINGLENGTH];
  CHAR dummy[LONGSTRINGLENGTH];
  LIGOTimeGPS fileStart;
  INT4 dummysec = 0;
  INT4 dummynan = 0;
  while (fscanf(cachefp,"%s %d %d",dummy,&dummysec,&dummynan)!=EOF) {
    if (strstr(dummy,uvar.binfile)) {
      idx = i;
      strcpy(filename,dummy);
      fileStart.gpsSeconds = dummysec;
      fileStart.gpsNanoSeconds = (INT4)1e9*(floor(1e-9*dummynan/uvar.tsamp + 0.5)*uvar.tsamp);    /* round to make sure we get samples at GPS seconds */
    }
  i++;
  }
  if (idx < 0) {
    LogPrintf(LOG_CRITICAL,"%s : failed to find binary input file %s in cache file %s\n",__func__,filename,uvar.cachefile);
    return 1;
  }
  fclose(cachefp);
  if (uvar.verbose) fprintf(stdout,"%s : found the requested binary file entry in the cache file.\n",__func__);

  /***********************************************************************************/
  /* setup the fixed binaryToSFT parameters */
  BinaryToSFTparams par;
  par.tsamp = uvar.tsamp;
  par.highpassf = uvar.highpassf;
  par.amp_inj = 0.0;
  par.f_inj = 0.0;
  par.asini_inj = 0.0;
  XLALGPSSetREAL8(&(par.tasc_inj),0);
  par.tref = fileStart;
  par.P_inj = 0;
  par.phi_inj = 0;
  par.r = NULL;

  /* setup the gridding over coherent times - which we make sure are integers */
  INT4 dummyT = uvar.Tmin;
  INT4 NT = 0;
  while (dummyT<uvar.Tmax) {
    dummyT = (INT4)((REAL8)dummyT*(1.0 + uvar.mismatch));
    NT++;
  }
  if (uvar.verbose) fprintf(stdout,"%s : Going to do %d different coherent lengths.\n",__func__,NT);

  /**********************************************************************************/
  /* OPEN INTERMEDIATE RESULTS FILE */
  /**********************************************************************************/
  CHAR intname[LONGSTRINGLENGTH];
  if (XLALOpenIntermittentResultsFile(&ifp,intname,newtemp,clargs,&uvar,&fileStart)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALOpenCoherentResultsFile() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  if (uvar.verbose) fprintf(stdout,"%s : opened coherent results file %s.\n",__func__,intname);

  /**********************************************************************************/
  /* loop over the different coherent lengths */
  INT4 currentT = uvar.Tmin;
  for (i=0;i<NT;i++) {

    if (uvar.verbose) fprintf(stdout,"%s : working on segment length %d/%d using a segment time of %d sec\n",__func__,i,NT,currentT);

    /* define SFT frequency bounds */
    REAL8 wings = LAL_TWOPI*uvar.maxasini/uvar.minorbperiod;
    REAL8 fmin_read = MINBAND*floor((uvar.freq - WINGS_FACTOR*uvar.freq*wings - (REAL8)uvar.blocksize/(REAL8)uvar.Tmin)/MINBAND);
    REAL8 fmax_read = MINBAND*ceil((uvar.freq + (REAL8)uvar.blocksize/(REAL8)uvar.Tmin + uvar.freqband + WINGS_FACTOR*(uvar.freq + uvar.freqband)*wings)/MINBAND);
    REAL8 fband_read = fmax_read - fmin_read;
    if (uvar.verbose) fprintf(stdout,"%s : reading in SFT frequency band [%f -> %f]\n",__func__,fmin_read,fmax_read);

    /* define the number of different start times */
    INT4 Ns = ceil(1.0/uvar.mismatch);
    INT4 Tstep = (INT4)floor(currentT/(REAL8)Ns);
    if (Tstep == 0) Ns = 1;
    else Ns = (INT4)ceil((REAL8)currentT/(REAL8)Tstep);
    if (uvar.verbose) fprintf(stdout,"%s : number of start times is %d and time step is %d sec\n",__func__,Ns,Tstep);
    par.tsft = currentT;
    par.freq = fmin_read;
    par.freqband = fband_read;

    /* loop over different start times - make sure they are integer GPS time for simplicity */
    LIGOTimeGPS currentStart;
    currentStart.gpsSeconds = (INT4)ceil((REAL8)fileStart.gpsSeconds + 1e-9*(REAL8)fileStart.gpsNanoSeconds);
    currentStart.gpsNanoSeconds = 0;
    for (k=0;k<Ns;k++) {

      if (uvar.verbose) fprintf(stdout,"%s : working on offset start %d/%d using a start time of %d %d sec\n",__func__,k,Ns,currentStart.gpsSeconds,currentStart.gpsNanoSeconds);
      memcpy(&(par.tstart),&currentStart,sizeof(LIGOTimeGPS));

      ParameterSpace pspace = empty_ParameterSpace;   /* the search parameter space */
      COMPLEX8TimeSeriesArray *dstimevec = NULL;      /* contains the downsampled inverse FFT'd SFTs */
      REAL4DemodulatedPowerVector *dmpower = NULL;   /* contains the demodulated power for all SFTs */
      GridParametersVector *freqgridparams = NULL; /* the coherent grid on the frequency derivitive parameter space */
      SFTVector *sftvec = NULL;
      INT8Vector *np = NULL;
      REAL8Vector *R = NULL;

      /**********************************************************************************/
      /* GENERATE SFTS FROM BINARY INPUT FILE */
      /**********************************************************************************/

      /* converts a binary input file into sfts */
      if (XLALBinaryToSFTVector(&sftvec,filename,&fileStart,&par,&np,&R)) {
        LogPrintf(LOG_CRITICAL,"%s : failed to convert binary input file %s to sfts\n",__func__,uvar.binfile);
        return 1;
      }
      XLALDestroyINT8Vector(np);
      XLALDestroyREAL8Vector(R);
      if (uvar.verbose) fprintf(stdout,"%s : generated SFTs for file %s\n",__func__,filename);

      /* if we have any SFTs */
      if (sftvec->length>0) {

        /* define SFT length and the start and span of the observations plus the definitive segment time */
        pspace.tseg = 1.0/sftvec->data[0].deltaF;
        memcpy(&(pspace.epoch),&(sftvec->data[0].epoch),sizeof(LIGOTimeGPS));
        pspace.span = XLALGPSDiff(&(sftvec->data[sftvec->length-1].epoch),&(sftvec->data[0].epoch)) + pspace.tseg;
        if (uvar.verbose) {
	  fprintf(stdout,"%s : SFT length = %f seconds\n",__func__,pspace.tseg);
 	  fprintf(stdout,"%s : entire dataset starts at GPS time %d contains %d SFTS and spans %.0f seconds\n",__func__,pspace.epoch.gpsSeconds,sftvec->length,pspace.span);
        }

        /**********************************************************************************/
        /* NORMALISE THE SFTS */
        /**********************************************************************************/

        /* compute the background noise using the sfts - this routine uses the running median at the edges to normalise the wings */
        if (XLALNormalizeSFTVect(sftvec,uvar.blocksize, 0)) {
    	  LogPrintf(LOG_CRITICAL,"%s : XLALNormaliseSFTVect() failed with error = %d\n",__func__,xlalErrno);
          return 1;
        }
        if (uvar.verbose) fprintf(stdout,"%s : normalised the SFTs\n",__func__);

        /**********************************************************************************/
        /* DEFINE THE BINARY PARAMETER SPACE */
        /**********************************************************************************/

        /* define the binary parameter space */
        if (XLALDefineBinaryParameterSpace(&(pspace.space),pspace.epoch,pspace.span,&uvar)) {
          LogPrintf(LOG_CRITICAL,"%s : XLALDefineBinaryParameterSpace() failed with error = %d\n",__func__,xlalErrno);
          return 1;
        }
        if (uvar.verbose) fprintf(stdout,"%s : defined binary parameter prior space\n",__func__);

        /**********************************************************************************/
        /* COMPUTE THE COARSE GRID ON FREQUENCY DERIVITIVES */
        /**********************************************************************************/

        /* compute the grid parameters for all SFTs */
        if (XLALComputeFreqGridParamsVector(&freqgridparams,pspace.space,sftvec,uvar.mismatch)) {
          LogPrintf(LOG_CRITICAL,"%s : XLALComputeFreqGridParams() failed with error = %d\n",__func__,xlalErrno);
          return 1;
        }
        if (uvar.verbose) fprintf(stdout,"%s : computed the grid parameters for the sfts\n",__func__);

        /**********************************************************************************/
        /* CONVERT ALL SFTS TO DOWNSAMPLED TIMESERIES */
        /**********************************************************************************/

        if (XLALSFTVectorToCOMPLEX8TimeSeriesArray(&dstimevec,sftvec)) {
          LogPrintf(LOG_CRITICAL,"%s : XLALSFTVectorToCOMPLEX8TimeSeriesArray() failed with error = %d\n",__func__,xlalErrno);
          return 1;
        }
        if (uvar.verbose) fprintf(stdout,"%s : converted SFTs to downsampled timeseries\n",__func__);

        /**********************************************************************************/
        /* COMPUTE THE STATISTICS ON THE COARSE GRID */
        /**********************************************************************************/

        /* compute the demodulated power on the frequency derivitive grid */
        if (XLALCOMPLEX8TimeSeriesArrayToDemodPowerVector(&dmpower,dstimevec,freqgridparams,ifp)) {
          LogPrintf(LOG_CRITICAL,"%s : XLALCOMPLEX8TimeSeriesArrayToDemodPowerVector() failed with error = %d\n",__func__,xlalErrno);
          return 1;
        }
        if (uvar.verbose) fprintf(stdout,"%s : computed the demodulated power\n",__func__);

        /**********************************************************************************/
        /* FREE MEMORY */
        /**********************************************************************************/

        /* free memory inside the loop */
        XLALFreeParameterSpace(&pspace);
        XLALFreeREAL4DemodulatedPowerVector(dmpower);
        for (m=0;m<(INT4)dstimevec->length;m++) {
 	  XLALDestroyCOMPLEX8TimeSeries(dstimevec->data[m]);
        }
        XLALFree(dstimevec->data);
        XLALFree(dstimevec);
        for (m=0;m<(INT4)freqgridparams->length;m++) {
          XLALFree(freqgridparams->segment[m]->grid);
          XLALFree(freqgridparams->segment[m]->prod);
          XLALFree(freqgridparams->segment[m]);
        }
        XLALFree(freqgridparams->segment);
        XLALFree(freqgridparams);
        if (uvar.verbose) fprintf(stdout,"%s : freed memory\n",__func__);
        XLALDestroySFTVector(sftvec);

      } /* end if statement on whether we have SFTs */

      /* update the start time of the segment */
      XLALGPSAdd(&currentStart,(REAL8)Tstep);

    } /* end loop over start times */

    /* update segment length */
    currentT = (INT4)((REAL8)currentT*(1.0 + uvar.mismatch));

  }  /* end loop over segment lengths */

  /* move the temporary directory to the final location */
  if (uvar.tempdir) {
    CHAR newoutputfile[LONGSTRINGLENGTH];
    snprintf(newoutputfile,LONGSTRINGLENGTH,"%s/IntermittentResults-%s-%d.txt",uvar.outputdir,uvar.outLabel,fileStart.gpsSeconds);
    if (rename(intname,newoutputfile)) {
      LogPrintf(LOG_CRITICAL,"%s : unable to move final results file %s -> %s.  Exiting.\n",__func__,intname,newoutputfile);
      return 1;
    }
  }

  /**********************************************************************************/
  /* FREE MEMORY */
  /**********************************************************************************/
  fclose(ifp);

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
  uvar->binfile = NULL;
  uvar->tempdir = NULL;
  uvar->tsamp = 2048;
  uvar->Tmin = 16;
  uvar->Tmax = 256;
  uvar->freq = 550;
  uvar->freqband = 1;
  uvar->highpassf = 40;
  uvar->mismatch = 0.1;
  uvar->thresh = 20;
  uvar->blocksize = 100;
  uvar->verbose = 0;

  /* ---------- register all user-variables ---------- */
  XLALRegisterUvarMember(outLabel, 	        STRING, 'n', REQUIRED, "'misc' entry in SFT-filenames or 'description' entry of frame filenames");
  XLALRegisterUvarMember(outputdir, 	        STRING, 'o', REQUIRED, "The output directory name");
  XLALRegisterUvarMember(binfile, 	        STRING, 'i', REQUIRED, "The input binary file name");
  XLALRegisterUvarMember(cachefile,            STRING, 'e', REQUIRED, "The input cache file name");
  XLALRegisterUvarMember(tempdir, 	        STRING, 'x', OPTIONAL, "The temporary directory");
  XLALRegisterUvarMember(freq,                   REAL8, 'f', REQUIRED, "The starting frequency (Hz)");
  XLALRegisterUvarMember(freqband,   	        REAL8, 'b', REQUIRED, "The frequency band (Hz)");
  XLALRegisterUvarMember(Tmin,                    INT4, 't', OPTIONAL, "The min length of segemnts (sec)");
  XLALRegisterUvarMember(Tmax,                    INT4, 'T', OPTIONAL, "The max length of segments (sec)");
  XLALRegisterUvarMember(tsamp,           	REAL8, 's', REQUIRED, "The sampling time (sec)");
  XLALRegisterUvarMember(highpassf,           	REAL8, 'H', OPTIONAL, "The high pass filter frequency");
  XLALRegisterUvarMember(minorbperiod,          	REAL8, 'p', REQUIRED, "The minimum orbital period (sec)");
  XLALRegisterUvarMember(maxorbperiod,          	REAL8, 'P', REQUIRED, "The maximum orbital period (sec)");
  XLALRegisterUvarMember(minasini,          	REAL8, 'a', REQUIRED, "The minimum projected orbital radius (sec)");
  XLALRegisterUvarMember(maxasini,          	REAL8, 'A', REQUIRED, "The maximum projected orbital radius (sec)");
  XLALRegisterUvarMember(tasc,                   REAL8, 'c', REQUIRED, "The best guess orbital time of ascension (rads)");
  XLALRegisterUvarMember(deltaorbphase,      	REAL8, 'C', OPTIONAL, "The orbital phase uncertainty (cycles)");
  XLALRegisterUvarMember(mismatch,        	REAL8, 'm', OPTIONAL, "The grid mismatch (0->1)");
  XLALRegisterUvarMember(thresh,         	REAL8, 'z', OPTIONAL, "The threshold on Leahy power");
  XLALRegisterUvarMember(blocksize,               INT4, 'B', OPTIONAL, "The running median block size");
  XLALRegisterUvarMember(verbose,                   BOOLEAN, 'v', OPTIONAL, "Output status to standard out");

  /* do ALL cmdline and cfgfile handling */
  BOOLEAN should_exit = 0;
  if (XLALUserVarReadAllInput(&should_exit, argc, argv)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALUserVarReadAllInput failed with error = %d\n",__func__,xlalErrno);
    return XLAL_EFAULT;
  }
  if (should_exit) exit(1);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}


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
				    CHAR *outputfile,           /**< [out] the name of the output file */
				    CHAR *outputdir,            /**< [in] the output directory name */
				    CHAR *clargs,               /**< [in] the command line args */
				    UserInput_t *uvar,		/**< UNDOCUMENTED */
				    LIGOTimeGPS *start		/**< UNDOCUMENTED */
				    )
{
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
    snprintf(outputfile,LONGSTRINGLENGTH,"%s/IntermittentResults-%s-%d-%04d_%03d_%04d_%03d.txt",
			   outputdir,(CHAR*)uvar->outLabel,start->gpsSeconds,min_freq_int,min_freq_mhz,max_freq_int,max_freq_mhz);

  }
  LogPrintf(LOG_DEBUG,"%s : output %s\n",__func__, outputfile);

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
