/*
*  Copyright (C) 2010 -- 2014 Evan Goetz
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
 * \file
 * \ingroup lalapps_pulsar_TwoSpect
 * \author Evan Goetz
 */

#include <sys/stat.h>

#include <lal/UserInput.h>
#include <lal/LALString.h>
#include <lal/Window.h>
#include <lal/DopplerScan.h>

#include <gsl/gsl_math.h>

#include <lalapps.h>

#include "SFTfunctions.h"
#include "IHS.h"
#include "candidates.h"
#include "antenna.h"
#include "templates.h"
#include "TwoSpect.h"
#include "statistics.h"
#include "upperlimits.h"
#include "vectormath.h"
#include "falsealarm.h"


//Global variables
FILE *LOG = NULL;

/**
 * Main program
 */
int main(int argc, char *argv[])
{
   INT4 ii, jj;               //counter variables
   LALStatus XLAL_INIT_DECL(status);
   time_t programstarttime, programendtime;
   struct tm *ptm;

   time(&programstarttime);
   ptm = localtime(&programstarttime);

   //Turn off gsl error handler
   gsl_set_error_handler_off();

   //Set up shop and read in the user input values
   UserInput_t XLAL_INIT_DECL(uvar);
   XLAL_CHECK ( readTwoSpectInputParams(&uvar, argc, argv) == XLAL_SUCCESS, XLAL_EFUNC );

   //Create directory
   INT4 dirstatus = mkdir(uvar.outdirectory, 0777);
   XLAL_CHECK( dirstatus == 0 || (dirstatus == -1 && errno == EEXIST), XLAL_EIO, "Couldn't create directory %s\n", uvar.outdirectory ) ;

   //Filenames for logfile and ULfile
   CHARVector *LOGFILENAME = NULL, *ULFILENAME = NULL;
   XLAL_CHECK( (LOGFILENAME = XLALCreateCHARVector(strlen(uvar.outdirectory)+strlen(uvar.outfilename)+3)) != NULL, XLAL_EFUNC);
   sprintf(LOGFILENAME->data, "%s/%s", uvar.outdirectory, uvar.outfilename);
   XLAL_CHECK( (ULFILENAME = XLALCreateCHARVector(strlen(uvar.outdirectory)+strlen(uvar.ULfilename)+3)) != NULL, XLAL_EFUNC );
   sprintf(ULFILENAME->data, "%s/%s", uvar.outdirectory, uvar.ULfilename);

   //Save args_info
   CHARVector *configfilecopy = NULL;
   XLAL_CHECK( (configfilecopy = XLALCreateCHARVector(strlen(uvar.outdirectory)+strlen(uvar.configCopy)+3)) != NULL, XLAL_EFUNC );
   sprintf(configfilecopy->data, "%s/%s", uvar.outdirectory, uvar.configCopy);
   FILE *INPUTVALS = NULL;
   XLAL_CHECK( (INPUTVALS = fopen(configfilecopy->data, "w")) != NULL, XLAL_EIO, "Failed to fopen %s for writing input parameter values\n", configfilecopy->data);
   CHAR *cfgFileStringCopy = NULL;
   XLAL_CHECK( (cfgFileStringCopy = XLALUserVarGetLog(UVAR_LOGFMT_CFGFILE)) != NULL, XLAL_EFUNC );
   fprintf(INPUTVALS, "%s", cfgFileStringCopy);
   fclose(INPUTVALS);
   XLALFree(cfgFileStringCopy);
   XLALDestroyCHARVector(configfilecopy);

   //Open log file
   XLAL_CHECK( (LOG = fopen(LOGFILENAME->data,"w")) != NULL, XLAL_EIO, "Failed to fopen %s for writing log file\n", LOGFILENAME->data );

   //print VCS info
   CHAR *VCSInfoString;
   XLAL_CHECK( (VCSInfoString = XLALGetVersionString(0)) != NULL, XLAL_EFUNC );
   fprintf(LOG, "%s\n", VCSInfoString);
   fprintf(stderr, "%s\n", VCSInfoString);
   XLALFree(VCSInfoString);

   //print start time
   fprintf(stderr, "Program executed on %s", asctime(ptm));
   fprintf(LOG, "Program executed on %s", asctime(ptm));

   //Set vectormath
   //if (uvar.vectorMath==1) XLAL_CHECK( XLALVectorDeviceSet(VECTORDEVICE_SSE) == XLAL_SUCCESS, XLAL_EFUNC );
   //else if (uvar.vectorMath==2) XLAL_CHECK( XLALVectorDeviceSet(VECTORDEVICE_AVX) == XLAL_SUCCESS, XLAL_EFUNC );
   //else XLAL_CHECK( XLALVectorDeviceSet(VECTORDEVICE_FPU) == XLAL_SUCCESS, XLAL_EFUNC );

   //Print the output directory
   fprintf(stderr, "Output directory: %s\n", uvar.outdirectory);
   fprintf(LOG, "Output directory: %s\n", uvar.outdirectory);

   //Set up the MultiLALDetector structure
   for (ii=0; ii<(INT4)uvar.IFO->length; ii++) {
      if (ii==0) {
         fprintf(LOG,"IFO = %s", uvar.IFO->data[ii]);
         fprintf(stderr,"IFO = %s", uvar.IFO->data[ii]);
      } else {
         fprintf(LOG,"%s", uvar.IFO->data[ii]);
         fprintf(stderr,"%s", uvar.IFO->data[ii]);
      }
      if (ii < (INT4)uvar.IFO->length-1) {
         fprintf(LOG, ",");
         fprintf(stderr, ",");
      } else {
         fprintf(LOG, "\n");
         fprintf(stderr, "\n");
      }
   }
   MultiLALDetector *detectors = NULL;
   XLAL_CHECK( (detectors = XLALMalloc(sizeof(MultiLALDetector))) != NULL, XLAL_ENOMEM );
   //XLAL_CHECK( XLALParseMultiLALDetector(detectors, uvar.IFO) == XLAL_SUCCESS, XLAL_EFUNC );
   detectors->length = uvar.IFO->length;
   for (ii=0; ii<(INT4)detectors->length; ii++) {
      if (strcmp("H1", uvar.IFO->data[ii])==0) {
         detectors->sites[ii] = lalCachedDetectors[LAL_LHO_4K_DETECTOR]; //H1
      } else if (strcmp("L1",uvar.IFO->data[ii])==0) {
         detectors->sites[ii] = lalCachedDetectors[LAL_LLO_4K_DETECTOR]; //L1
      } else if (strcmp("V1", uvar.IFO->data[ii])==0) {
         detectors->sites[ii] = lalCachedDetectors[LAL_VIRGO_DETECTOR];  //V1
      } else if (strcmp("H2", uvar.IFO->data[ii])==0) {
         detectors->sites[ii] = lalCachedDetectors[LAL_LHO_2K_DETECTOR]; //H2
      } else if (strcmp("H2r", uvar.IFO->data[ii])==0) {
         LALDetector H2 = lalCachedDetectors[LAL_LHO_2K_DETECTOR]; //H2 rotated
         H2.frDetector.xArmAzimuthRadians -= 0.25*LAL_PI;
         H2.frDetector.yArmAzimuthRadians -= 0.25*LAL_PI;
         memset(&(H2.frDetector.name), 0, sizeof(CHAR)*LALNameLength);
         snprintf(H2.frDetector.name, LALNameLength, "%s", "LHO_2k_rotatedPiOver4");
         XLAL_CHECK( (XLALCreateDetector(&(detectors->sites[ii]), &(H2.frDetector), LALDETECTORTYPE_IFODIFF)) != NULL, XLAL_EFUNC );
      } else {
         XLAL_ERROR(XLAL_EINVAL, "Not using valid interferometer! Expected 'H1', 'H2', 'H2r' (rotated H2), 'L1', or 'V1' not %s.\n", uvar.IFO->data[ii]);
      }
   }

   //Allocate GSL RNG
   gsl_rng *rng = NULL;
   XLAL_CHECK( (rng = gsl_rng_alloc(gsl_rng_mt19937)) != NULL, XLAL_EFUNC );
   //Random seed value settings for random number generator
   //If the chooseSeed option was given, then:
   //seed = (IFO multiplier)*fabs(round(fmin + fspan + Pmin + Pmax + dfmin + dfmax + alpha + delta))
   if (uvar.chooseSeed) {
      UINT8 IFOmultiplier = 0;
      if (strcmp(uvar.IFO->data[0], "H1")==0) IFOmultiplier = 1;            //H1 gets multiplier 1
      else if (strcmp(uvar.IFO->data[0], "L1")==0) IFOmultiplier = 2;       //L1 gets multiplier 2
      else IFOmultiplier = 3;                                                  //V1 gets multiplier 3
      uvar.randSeed = IFOmultiplier*(UINT8)fabs(round(uvar.fmin + uvar.fspan + uvar.Pmin + uvar.Pmax + uvar.dfmin + uvar.dfmax));
   }
   gsl_rng_set(rng, uvar.randSeed);     //Set the random number generator with the given seed

   //Print to stderr the parameters of the search
   fprintf(stderr, "Tobs = %f sec\n", uvar.Tobs);
   fprintf(stderr, "Tsft = %f sec\n", uvar.Tsft);
   fprintf(stderr, "SFToverlap = %f sec\n", uvar.SFToverlap);
   fprintf(stderr, "fmin = %f Hz\n", uvar.fmin);
   fprintf(stderr, "fspan = %f Hz\n", uvar.fspan);
   fprintf(stderr, "Pmin = %f s\n", uvar.Pmin);
   fprintf(stderr, "Pmax = %f s\n", uvar.Pmax);
   fprintf(stderr, "dfmin = %f Hz\n", uvar.dfmin);
   fprintf(stderr, "dfmax = %f Hz\n", uvar.dfmax);
   fprintf(stderr, "Running median blocksize = %d\n", uvar.blksize);
   fprintf(stderr, "FFT plan flag = %d\n", uvar.FFTplanFlag);
   fprintf(stderr, "RNG seed: %d\n", uvar.randSeed);

   //Initialize ephemeris data structure
   EphemerisData *edat = NULL;
   XLAL_CHECK( (edat = XLALInitBarycenter(uvar.ephemEarth, uvar.ephemSun)) != NULL, XLAL_EFUNC );

   //Maximum detector velocity in units of c from start of observation time - Tsft to end of observation + Tsft
   REAL4 Vmax = 0.0;
   for (ii=0; ii<(INT4)detectors->length; ii++) {
      REAL4 detectorVmax = CompDetectorVmax(uvar.t0-uvar.Tsft, uvar.Tsft, uvar.SFToverlap, uvar.Tobs+2.0*uvar.Tsft, detectors->sites[ii], edat);
      XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC, "CompDetectorVmax() failed\n" );
      if (detectorVmax > Vmax) Vmax = detectorVmax;
   }

   //Assume maximum possible bin shift
   INT4 maxbinshift = (INT4)round(Vmax * (uvar.fmin+uvar.fspan) * uvar.Tsft) + 1;

   //Parameters for the sky-grid from a point/polygon or a sky-grid file
   DopplerSkyScanInit XLAL_INIT_DECL(scanInit);
   DopplerSkyScanState XLAL_INIT_DECL(scan);
   PulsarDopplerParams dopplerpos;
   if (XLALUserVarWasSet(&uvar.skyRegion)) {
      scanInit.gridType = GRID_ISOTROPIC;     //Default value for an approximate-isotropic grid
      scanInit.skyRegionString = uvar.skyRegion;      //"allsky" = Default value for all-sky search
      scanInit.numSkyPartitions = 1;   //Default value so sky is not broken into chunks
      scanInit.Freq = uvar.fmin+0.5*uvar.fspan;  //Midpoint of the frequency band
      scanInit.dAlpha = 0.5/((uvar.fmin+uvar.fspan) * uvar.Tsft * Vmax);
      scanInit.dDelta = scanInit.dAlpha;
      fprintf(LOG, "Sky region = %s\n", uvar.skyRegion);
      fprintf(stderr, "Sky region = %s\n", uvar.skyRegion);
   } else {
      scanInit.gridType = GRID_FILE_SKYGRID;
      scanInit.skyGridFile = uvar.skyRegionFile;
      scanInit.numSkyPartitions = 1;   //Default value so sky is not broken into chunks
      scanInit.Freq = uvar.fmin+0.5*uvar.fspan;  //Midpoint of the frequency band
      fprintf(LOG, "Sky file = %s\n", uvar.skyRegionFile);
      fprintf(stderr, "Sky file = %s\n", uvar.skyRegionFile);
   }

   //Initialize the sky-grid
   InitDopplerSkyScan(&status, &scan, &scanInit);
   XLAL_CHECK( status.statusCode == 0, XLAL_EFUNC );

   //Start at first location
   XLAL_CHECK( XLALNextDopplerSkyPos(&dopplerpos, &scan) == XLAL_SUCCESS, XLAL_EFUNC );

   //Basic units
   REAL4 tempfspan = uvar.fspan + 2.0*uvar.dfmax + (uvar.blksize-1 + 12)/uvar.Tsft;     //= fspan+2*dfmax+extrabins + running median blocksize-1 (Hz)
   INT4 tempnumfbins = (INT4)round(tempfspan*uvar.Tsft)+1;                        //= number of bins in tempfspan
   fprintf(LOG, "FAR for templates = %g\n", uvar.tmplfar);
   fprintf(stderr, "FAR for templates = %g\n", uvar.tmplfar);
   
   //Determine band size to get the SFT data (remember to get extra bins because of the running median and the bin shifts due to detector velocity) with nudge of 0.1/Tsft for rounding issues
   REAL8 minfbin = round(uvar.fmin*uvar.Tsft - uvar.dfmax*uvar.Tsft - 0.5*(uvar.blksize-1) - (REAL8)(maxbinshift) - 6.0)/uvar.Tsft + 0.1/uvar.Tsft;
   REAL8 maxfbin = round((uvar.fmin + uvar.fspan)*uvar.Tsft + uvar.dfmax*uvar.Tsft + 0.5*(uvar.blksize-1) + (REAL8)(maxbinshift) + 6.0)/uvar.Tsft - 0.1/uvar.Tsft;

   //Allocate memory for ffdata structure
   ffdataStruct *ffdata = NULL;
   XLAL_CHECK( (ffdata = new_ffdata(&uvar)) != NULL, XLAL_EFUNC );

   //Maximum number of IHS values to sum = twice the maximum modulation depth
   //Minimum number of IHS values to sum = twice the minimum modulation depth
   INT4 maxrows = (INT4)round(2.0*uvar.dfmax*uvar.Tsft)+1;
   fprintf(LOG, "Maximum row width to be searched = %d\n", maxrows);
   fprintf(stderr, "Maximum row width to be searched = %d\n", maxrows);

   //TF normalization
   ffdata->tfnormalization = 2.0/uvar.Tsft/(uvar.avesqrtSh*uvar.avesqrtSh);

   //Read in the T-F data from SFTs or generate SFTs
   REAL4Vector *tfdata = NULL;
   MultiSFTVector *multiSFTvector = NULL;
   if (!XLALUserVarWasSet(&uvar.injectionSources) && XLALUserVarWasSet(&uvar.inputSFTs) && detectors->length<2 && !uvar.gaussNoiseWithSFTgaps && !XLALUserVarWasSet(&uvar.timestampsFile)) XLAL_CHECK( (tfdata = readInSFTs(&uvar, ffdata->tfnormalization, minfbin, maxfbin)) != NULL, XLAL_EFUNC );
   else if (!XLALUserVarWasSet(&uvar.injectionSources) && XLALUserVarWasSet(&uvar.inputSFTs) && !uvar.gaussNoiseWithSFTgaps && !XLALUserVarWasSet(&uvar.timestampsFile)) XLAL_CHECK( (multiSFTvector = getMultiSFTVector(&uvar, minfbin, maxfbin)) != NULL, XLAL_EFUNC );
   else {
      XLAL_CHECK( (multiSFTvector = generateSFTdata(&uvar, detectors, edat, maxbinshift, rng)) != NULL, XLAL_EFUNC );
      //Convert SFTs to powers
      if (detectors->length==1) {
         XLAL_CHECK( (tfdata = convertSFTdataToPowers(multiSFTvector->data[0], &uvar, ffdata->tfnormalization)) != NULL, XLAL_EFUNC );
         XLALDestroyMultiSFTVector(multiSFTvector);
      }
   } // end load data or generate data

   //Print SFT times, if requested by user
   if (uvar.printSFTtimes) {
      CHARVector *outputfile = NULL;
      XLAL_CHECK( (outputfile = XLALCreateCHARVector(strlen(uvar.outdirectory)+25)) != NULL, XLAL_EFUNC );
      for (ii=0; ii<(INT4)detectors->length; ii++) {
         sprintf(outputfile->data, "%s/%s-%s", uvar.outdirectory, uvar.IFO->data[ii], "inputSFTtimes.dat");
         FILE *INSFTTIMES = NULL;
         XLAL_CHECK( (INSFTTIMES = fopen(outputfile->data, "w")) != NULL, XLAL_EIO, "Couldn't fopen %s for writing input SFT start times", outputfile->data );
         if (detectors->length==1) {
            INT4 sftlength = tfdata->length/ffdata->numffts;
            for (jj=0; jj<ffdata->numffts; jj++) if (tfdata->data[jj*sftlength]!=0.0) fprintf(INSFTTIMES, "%9d 0\n", (INT4)round(uvar.t0+jj*(uvar.Tsft-uvar.SFToverlap)));
         } else {
            for (jj=0; jj<(INT4)(multiSFTvector->data[ii]->length); jj++) fprintf(INSFTTIMES, "%9d 0\n", (INT4)round(XLALGPSGetREAL8(&(multiSFTvector->data[ii]->data[jj].epoch))));
         }
         fclose(INSFTTIMES);
      }
      XLALDestroyCHARVector(outputfile);
   }

   //Line detection only if there is valid noise to be looking at
   INT4Vector *lines = NULL;
   INT4 heavilyContaminatedBand = 0;
   if (XLALUserVarWasSet(&uvar.lineDetection) && !uvar.signalOnly && detectors->length==1) {
      lines = detectLines_simple(tfdata, ffdata, &uvar);
      XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
      if (lines!=NULL) {
         fprintf(LOG, "WARNING: %d line(s) found.\n", lines->length);
         fprintf(stderr, "WARNING: %d line(s) found.\n", lines->length);
         if ((REAL4)lines->length/(ffdata->numfbins + 2*maxbinshift) >= 0.1) {
            heavilyContaminatedBand = 1;
            fprintf(LOG, "WARNING: Band is heavily contaminated by artifacts.\n");
            fprintf(stderr, "WARNING: Band is heavily contaminated by artifacts.\n");
         }
      }
   }

   //If the band is heavily contaminated by lines, don't do any follow up.
   if (heavilyContaminatedBand) uvar.IHSonly = 1;

   //Calculate the running mean values of the SFTs (output here is smaller than initialTFdata). Here,
   //numfbins needs to be the bins you expect to come out of the running means -- the band you are going
   //to search!
   REAL4Vector *background = NULL;
   if (detectors->length==1) {
      XLAL_CHECK( (background = XLALCreateREAL4Vector(ffdata->numffts*(ffdata->numfbins+2*maxbinshift))) != NULL, XLAL_EFUNC );
      if (!uvar.signalOnly) XLAL_CHECK( tfRngMeans(background, tfdata, ffdata->numffts, ffdata->numfbins+2*maxbinshift, uvar.blksize) == XLAL_SUCCESS, XLAL_EFUNC );
      else memset(background->data, 0, sizeof(REAL4)*background->length);
   }

   //TEST: Try cleaning lines
   /* XLAL_CHECK( cleanLines(tfdata, background, lines, &uvar) == XLAL_SUCCESS, XLAL_EFUNC ); */
   /* if (lines!=NULL) { */
   /*    if (!uvar.signalOnly) XLAL_CHECK( tfRngMeans(background, tfdata, ffdata->numffts, ffdata->numfbins + 2*maxbinshift, uvar.blksize) == XLAL_SUCCESS, XLAL_EFUNC ); */
   /*    else memset(background->data, 0, sizeof(REAL4)*background->length); */
   /* } */
   /* XLALDestroyINT4Vector(lines); */
   /* lines = NULL; */

   //Existing SFTs listed in this vector
   INT4Vector *sftexist = NULL;
   INT4Vector *indexValuesOfExistingSFTs = NULL;
   REAL4 frac_tobs_complete = 1.0;
   if (!uvar.signalOnly && detectors->length==1) {
      XLAL_CHECK( (sftexist = existingSFTs(tfdata, (UINT4)ffdata->numffts)) != NULL, XLAL_EFUNC );
      INT4 totalincludedsftnumber = 0;
      for (ii=0; ii<(INT4)sftexist->length; ii++) if (sftexist->data[ii]==1) totalincludedsftnumber++;
      frac_tobs_complete = (REAL4)totalincludedsftnumber/(REAL4)sftexist->length;
      fprintf(LOG, "Duty factor of usable SFTs = %f\n", frac_tobs_complete);
      fprintf(stderr, "Duty factor of usable SFTs = %f\n", frac_tobs_complete);
      if (frac_tobs_complete<0.1) {
         fprintf(stderr, "%s: The useable SFTs cover less than 10 percent of the total observation time\n", __func__);
         fprintf(LOG, "%s: The useable SFTs cover less than 10 percent of the total observation time\n", __func__);
         return 0;
      }

      //Index values of existing SFTs
      XLAL_CHECK( (indexValuesOfExistingSFTs = XLALCreateINT4Vector(totalincludedsftnumber)) != NULL, XLAL_EFUNC );
      jj = 0;
      for (ii=0; ii<(INT4)sftexist->length; ii++) {
         if (sftexist->data[ii] == 1) {
            indexValuesOfExistingSFTs->data[jj] = ii;
            jj++;
         }
      }
   }

   //I wrote this to compensate for a bad input of the expected noise floor
   REAL8 backgroundmeannormfactor = 1.0;
   if (!uvar.signalOnly && detectors->length==1) {
      REAL4 harmonicMean = 0.0;
      XLAL_CHECK( calcHarmonicMean(&harmonicMean, background, ffdata->numfbins + 2*maxbinshift, ffdata->numffts) == XLAL_SUCCESS, XLAL_EFUNC );
      backgroundmeannormfactor = 1.0/harmonicMean;
      ffdata->tfnormalization *= backgroundmeannormfactor;
   }

   //Need to reduce the original TF data to remove the excess bins used for running median calculation. Normalize the TF as the same as the background was normalized
   REAL4Vector *usableTFdata = NULL;
   if (detectors->length==1) {
      XLAL_CHECK( (usableTFdata = XLALCreateREAL4Vector(background->length)) != NULL, XLAL_EFUNC );
      for (ii=0; ii<ffdata->numffts; ii++) memcpy(&(usableTFdata->data[ii*(ffdata->numfbins+2*maxbinshift)]), &(tfdata->data[ii*(tempnumfbins+2*maxbinshift) + (INT4)round(0.5*(uvar.blksize-1))]), sizeof(REAL4)*(ffdata->numfbins+2*maxbinshift));
      for (ii=0; ii<(INT4)usableTFdata->length; ii++) {
         if (usableTFdata->data[ii]!=0.0) {
            usableTFdata->data[ii] *= backgroundmeannormfactor;
            background->data[ii] *= backgroundmeannormfactor;
         }
      }
      //At this point the TF plane and the running median calculation are the same size=numffts*(numfbins + 2*maxbinshift)
      //We can delete the originally loaded SFTs since we have the usableTFdata saved
      XLALDestroyREAL4Vector(tfdata);
   }

   //Print out data product if requested
   if (uvar.printData && detectors->length==1) {
      CHARVector *outputfile = NULL;
      XLAL_CHECK( (outputfile = XLALCreateCHARVector(strlen(uvar.outdirectory)+15)) != NULL, XLAL_EFUNC );
      sprintf(outputfile->data, "%s/%s", uvar.outdirectory, "tfdata.dat");
      XLAL_CHECK( printREAL4Vector2File(usableTFdata, outputfile->data) == XLAL_SUCCESS, XLAL_EFUNC );
      XLALDestroyCHARVector(outputfile);
   }

   //Do mean subtraction of TFdata here--modifies the usableTFdata vector!!!
   if (detectors->length==1) XLAL_CHECK( tfMeanSubtract(usableTFdata, background, ffdata->numffts, ffdata->numfbins+2*maxbinshift) == XLAL_SUCCESS, XLAL_EFUNC );

   //Initialize reused values
   ihsMaximaStruct *ihsmaxima = NULL;
   ihsfarStruct *ihsfarstruct = NULL;
   REAL4VectorAligned *aveNoise = NULL;
   INT4Vector *binshifts = NULL;
   LIGOTimeGPS tStart;
   MultiLIGOTimeGPSVector *multiTimestamps = NULL;
   MultiDetectorStateSeries *multiStateSeries = NULL;
   TwoSpectTemplateVector *templateVec = NULL;
   candidateVector *gaussCandidates1 = NULL, *gaussCandidates2 = NULL, *gaussCandidates3 = NULL, *gaussCandidates4 = NULL, *exactCandidates1 = NULL, *exactCandidates2 = NULL, *ihsCandidates = NULL;
   UpperLimitVector *upperlimits = NULL;
   REAL4FFTPlan *secondFFTplan = NULL;
   XLAL_CHECK( (ihsmaxima = new_ihsMaxima(ffdata->numfbins, maxrows)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihsfarstruct = new_ihsfarStruct(maxrows, &uvar)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (aveNoise = XLALCreateREAL4VectorAligned(ffdata->numfprbins, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (binshifts = XLALCreateINT4Vector(ffdata->numffts)) != NULL, XLAL_EFUNC );
   XLALGPSSetREAL8 ( &tStart, uvar.t0 );
   XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC, "XLALGPSSetREAL8 failed\n" );
   XLAL_CHECK( (multiTimestamps = XLALMakeMultiTimestamps(tStart, uvar.Tobs, uvar.Tsft, uvar.SFToverlap, detectors->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (multiStateSeries = XLALGetMultiDetectorStates(multiTimestamps, detectors, edat, uvar.SFToverlap)) != NULL, XLAL_EFUNC );
   if (XLALUserVarWasSet(&(uvar.templatebankfile))) {
      XLAL_CHECK( (templateVec = readTwoSpectTemplateVector(uvar.templatebankfile)) != NULL, XLAL_EFUNC );
      XLAL_CHECK( templateVec->Tsft==uvar.Tsft && templateVec->SFToverlap==uvar.SFToverlap && templateVec->Tobs==uvar.Tobs, XLAL_EINVAL, "Template bank %s doesn't match input parameters\n", uvar.templatebankfile );
   }
   XLAL_CHECK( (gaussCandidates1 = new_candidateVector(100)) != NULL, XLAL_EFUNC, "new_CandidateVector(%d) failed\n", 100 );
   XLAL_CHECK( (gaussCandidates2 = new_candidateVector(100)) != NULL, XLAL_EFUNC, "new_CandidateVector(%d) failed\n", 100 );
   XLAL_CHECK( (gaussCandidates3 = new_candidateVector(100)) != NULL, XLAL_EFUNC, "new_CandidateVector(%d) failed\n", 100 );
   XLAL_CHECK( (gaussCandidates4 = new_candidateVector(100)) != NULL, XLAL_EFUNC, "new_CandidateVector(%d) failed\n", 100 );
   XLAL_CHECK( (exactCandidates1 = new_candidateVector(100)) != NULL, XLAL_EFUNC, "new_CandidateVector(%d) failed\n", 100 );
   XLAL_CHECK( (exactCandidates2 = new_candidateVector(100)) != NULL, XLAL_EFUNC, "new_CandidateVector(%d) failed\n", 100 );
   XLAL_CHECK( (ihsCandidates = new_candidateVector(100)) != NULL, XLAL_EFUNC, "new_CandidateVector(%d) failed\n", 100 );
   XLAL_CHECK( (upperlimits = new_UpperLimitVector(1)) != NULL, XLAL_EFUNC, "new_UpperLimitVector(%d) failed.\n", 1 );
   XLAL_CHECK( (secondFFTplan = XLALCreateForwardREAL4FFTPlan(ffdata->numffts, uvar.FFTplanFlag)) != NULL, XLAL_EFUNC );
   LIGOTimeGPS refTime = multiTimestamps->data[0]->data[0];

   //Initialize to -1.0 for far just at the start
   ihsfarstruct->ihsfar->data[0] = -1.0;
   REAL4 antweightsrms = 0.0;

   //Davies' algorithm uses an internal error code
   //Later upgrade: remove this because we do our own error handling
   INT4 proberrcode = 0;

   //Set TF normalization
   if (detectors->length==1) ffdata->tfnormalization *= 0.5*uvar.Tsft;

   //Antenna normalization (determined from injections on H1 at ra=0, dec=0, with circular polarization)
   //When doing linear polarizations, the IHS factor needs to be 25.2*1.082 and this antenna weights
   //function needs to be set to use linear polarization.
   SkyPosition skypos0 = {0.0, 0.0, COORDINATESYSTEM_EQUATORIAL};
   REAL4VectorAligned *antweightsforihs2h0 = NULL;
   XLAL_CHECK( (antweightsforihs2h0 = XLALCreateREAL4VectorAligned(ffdata->numffts, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( CompAntennaPatternWeights((REAL4Vector*)antweightsforihs2h0, skypos0, uvar.t0, uvar.Tsft, uvar.SFToverlap, uvar.Tobs, 0, 0.0, lalCachedDetectors[LAL_LHO_4K_DETECTOR]) == XLAL_SUCCESS, XLAL_EFUNC );

   //Set skycounter to -1 at the start
   INT4 skycounter = -1;

   //Print message that we start the analysis
   fprintf(LOG, "Starting TwoSpect analysis...\n");
   fprintf(stderr, "Starting TwoSpect analysis...\n");

   //Search over the sky region (outer loop of single TwoSpect program instance)
   while (scan.state != STATE_FINISHED) {
      fprintf(LOG, "Sky location: RA = %g, DEC = %g\n", dopplerpos.Alpha, dopplerpos.Delta);
      fprintf(stderr, "Sky location: RA = %g, DEC = %g\n", dopplerpos.Alpha, dopplerpos.Delta);
      skycounter++;

      SkyPosition skypos;
      skypos.longitude = dopplerpos.Alpha;
      skypos.latitude = dopplerpos.Delta;
      skypos.system = COORDINATESYSTEM_EQUATORIAL;

      MultiSSBtimes *multissb = NULL;
      XLAL_CHECK( (multissb = XLALGetMultiSSBtimes(multiStateSeries, skypos, refTime, SSBPREC_RELATIVISTICOPT)) != NULL, XLAL_EFUNC );

      //Compute antenna pattern weights. If antennaOff input flag is given, then set all values equal to 1.0
      REAL4VectorAligned *antweights = NULL;
      XLAL_CHECK( (antweights = XLALCreateREAL4VectorAligned(ffdata->numffts, 32)) != NULL, XLAL_EFUNC );
      memset(antweights->data, 0, antweights->length*sizeof(REAL4));
      if (uvar.antennaOff) for (ii=0; ii<(INT4)antweights->length; ii++) antweights->data[ii] = 1.0;
      else {
         if (XLALUserVarWasSet(&uvar.linPolAngle)) XLAL_CHECK( CompAntennaPatternWeights((REAL4Vector*)antweights, skypos, uvar.t0, uvar.Tsft, uvar.SFToverlap, uvar.Tobs, 1, uvar.linPolAngle, detectors->sites[0]) == XLAL_SUCCESS, XLAL_EFUNC );
         else XLAL_CHECK( CompAntennaPatternWeights((REAL4Vector*)antweights, skypos, uvar.t0, uvar.Tsft, uvar.SFToverlap, uvar.Tobs, 0, 0.0, detectors->sites[0]) == XLAL_SUCCESS, XLAL_EFUNC );
      }

      //Compute the bin shifts for each SFT
      XLAL_CHECK( CompBinShifts2(binshifts, multissb->data[0], uvar.fmin + 0.5*uvar.fspan, uvar.Tsft, uvar.dopplerMultiplier) == XLAL_SUCCESS, XLAL_EFUNC );

      //Coherent combination
      if (detectors->length>=2) {
          //TF normalization
         ffdata->tfnormalization = 2.0/uvar.Tsft/(uvar.avesqrtSh*uvar.avesqrtSh);

         //Get multiAMcoefficients
         MultiAMCoeffs *multiAMcoefficients = NULL;
         XLAL_CHECK( (multiAMcoefficients = XLALComputeMultiAMCoeffs(multiStateSeries, NULL, skypos)) != NULL, XLAL_EFUNC );

         //Coherenly combine the SFT data
         if (XLALUserVarWasSet(&uvar.assumeNScosi) && XLALUserVarWasSet(&uvar.assumeNSpsi)) XLAL_CHECK( (tfdata = coherentlyAddSFTs(multiSFTvector, multissb, multiAMcoefficients, uvar.cosiSignCoherent, &uvar.assumeNScosi, &uvar.assumeNSpsi, &uvar)) != NULL, XLAL_EFUNC );
         else if (XLALUserVarWasSet(&uvar.assumeNScosi) && !XLALUserVarWasSet(&uvar.assumeNSpsi)) XLAL_CHECK( (tfdata = coherentlyAddSFTs(multiSFTvector, multissb, multiAMcoefficients, uvar.cosiSignCoherent, &uvar.assumeNScosi, NULL, &uvar)) != NULL, XLAL_EFUNC );
         else if (!XLALUserVarWasSet(&uvar.assumeNScosi) && XLALUserVarWasSet(&uvar.assumeNSpsi)) XLAL_CHECK( (tfdata = coherentlyAddSFTs(multiSFTvector, multissb, multiAMcoefficients, uvar.cosiSignCoherent, NULL, &uvar.assumeNSpsi, &uvar)) != NULL, XLAL_EFUNC );
         else XLAL_CHECK( (tfdata = coherentlyAddSFTs(multiSFTvector, multissb, multiAMcoefficients, uvar.cosiSignCoherent, NULL, NULL, &uvar)) != NULL, XLAL_EFUNC );

         XLALDestroyMultiAMCoeffs(multiAMcoefficients);

         //Continue like there is just one detector

         //Check for lines
         if (XLALUserVarWasSet(&uvar.lineDetection) && !uvar.signalOnly) {
            lines = detectLines_simple(tfdata, ffdata, &uvar);
            XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
            if (lines!=NULL) {
               fprintf(LOG, "WARNING: %d line(s) found.\n", lines->length);
               fprintf(stderr, "WARNING: %d line(s) found.\n", lines->length);
               if ((REAL4)lines->length/(ffdata->numfbins + 2*maxbinshift) >= 0.1) {
                  heavilyContaminatedBand = 1;
                  fprintf(LOG, "WARNING: Band is heavily contaminated by artifacts.\n");
                  fprintf(stderr, "WARNING: Band is heavily contaminated by artifacts.\n");
               }
            }
         }
         if (heavilyContaminatedBand) uvar.IHSonly = 1;

         //Compute background
         XLAL_CHECK( (background = XLALCreateREAL4Vector(ffdata->numffts*(ffdata->numfbins+2*maxbinshift))) != NULL, XLAL_EFUNC );
         if (!uvar.signalOnly) XLAL_CHECK( tfRngMeans(background, tfdata, ffdata->numffts, ffdata->numfbins+2*maxbinshift, uvar.blksize) == XLAL_SUCCESS, XLAL_EFUNC );
         else memset(background->data, 0, sizeof(REAL4)*background->length);

         if (!uvar.signalOnly) {
            XLAL_CHECK( (sftexist = existingSFTs(tfdata, (UINT4)ffdata->numffts)) != NULL, XLAL_EFUNC );
            INT4 totalincludedsftnumber = 0;
            for (ii=0; ii<(INT4)sftexist->length; ii++) if (sftexist->data[ii]==1) totalincludedsftnumber++;
            frac_tobs_complete = (REAL4)totalincludedsftnumber/(REAL4)sftexist->length;
            fprintf(LOG, "Duty factor of usable SFTs = %f\n", frac_tobs_complete);
            fprintf(stderr, "Duty factor of usable SFTs = %f\n", frac_tobs_complete);
            if (frac_tobs_complete<0.1) {
               fprintf(stderr, "%s: The useable SFTs cover less than 10 percent of the total observation time\n", __func__);
               fprintf(LOG, "%s: The useable SFTs cover less than 10 percent of the total observation time\n", __func__);
               return 0;
            }

            //Index values of existing SFTs
            XLAL_CHECK( (indexValuesOfExistingSFTs = XLALCreateINT4Vector(totalincludedsftnumber)) != NULL, XLAL_EFUNC );
            jj = 0;
            for (ii=0; ii<(INT4)sftexist->length; ii++) {
               if (sftexist->data[ii] == 1) {
                  indexValuesOfExistingSFTs->data[jj] = ii;
                  jj++;
               }
            }

            REAL4 harmonicMean = 0.0;
            XLAL_CHECK( calcHarmonicMean(&harmonicMean, background, ffdata->numfbins + 2*maxbinshift, ffdata->numffts) == XLAL_SUCCESS, XLAL_EFUNC );
            backgroundmeannormfactor = 1.0/harmonicMean;
            ffdata->tfnormalization *= backgroundmeannormfactor;
         }

         XLAL_CHECK( (usableTFdata = XLALCreateREAL4Vector(background->length)) != NULL, XLAL_EFUNC );
         for (ii=0; ii<ffdata->numffts; ii++) memcpy(&(usableTFdata->data[ii*(ffdata->numfbins+2*maxbinshift)]), &(tfdata->data[ii*(tempnumfbins+2*maxbinshift) + (INT4)round(0.5*(uvar.blksize-1))]), sizeof(REAL4)*(ffdata->numfbins+2*maxbinshift));
         for (ii=0; ii<(INT4)usableTFdata->length; ii++) {
            if (usableTFdata->data[ii]!=0.0) {
               usableTFdata->data[ii] *= backgroundmeannormfactor;
               background->data[ii] *= backgroundmeannormfactor;
            }
         }
         XLALDestroyREAL4Vector(tfdata);

         //Print out data product if requested
         if (uvar.printData) {
            CHARVector *outputfile = NULL;
            XLAL_CHECK( (outputfile = XLALCreateCHARVector(strlen(uvar.outdirectory)+15)) != NULL, XLAL_EFUNC );
            sprintf(outputfile->data, "%s/%s", uvar.outdirectory, "tfdata.dat");
            XLAL_CHECK( printREAL4Vector2File(usableTFdata, outputfile->data) == XLAL_SUCCESS, XLAL_EFUNC );
            XLALDestroyCHARVector(outputfile);
         }

         //Mean subtraction
         XLAL_CHECK( tfMeanSubtract(usableTFdata, background, ffdata->numffts, ffdata->numfbins+2*maxbinshift) == XLAL_SUCCESS, XLAL_EFUNC );

         ffdata->tfnormalization *= 0.5*uvar.Tsft;
      }
      /////

      XLALDestroyMultiSSBtimes(multissb);

      //Track identified lines
      REAL4 fbin0 = (REAL4)(round(uvar.fmin*uvar.Tsft - uvar.dfmax*uvar.Tsft - 6.0 - 0.5*(uvar.blksize-1) - (REAL8)(maxbinshift))/uvar.Tsft);
      REAL4 df = 1.0/uvar.Tsft;
      REAL4VectorSequence *trackedlines = NULL;
      if (lines!=NULL) {
         XLAL_CHECK( (trackedlines = trackLines(lines, binshifts, fbin0, df)) != NULL, XLAL_EFUNC );
         XLALDestroyINT4Vector(lines);
      }

      //Calculate antenna RMS value
      REAL4 currentAntWeightsRMS = 0.0;
      XLAL_CHECK( calcRms(&currentAntWeightsRMS, (REAL4Vector*)antweights) == XLAL_SUCCESS, XLAL_EFUNC );

      //Slide SFTs here -- need to slide the data and the estimated background
      REAL4Vector *TFdata_slided = NULL, *background_slided = NULL;
      XLAL_CHECK( (TFdata_slided = XLALCreateREAL4Vector(ffdata->numffts*ffdata->numfbins)) != NULL, XLAL_EFUNC );
      XLAL_CHECK( (background_slided = XLALCreateREAL4Vector(TFdata_slided->length)) != NULL, XLAL_EFUNC );
      XLAL_CHECK( slideTFdata(TFdata_slided, &uvar, usableTFdata, binshifts) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( slideTFdata(background_slided, &uvar, background, binshifts) == XLAL_SUCCESS, XLAL_EFUNC );

      if (detectors->length>1) {
         XLALDestroyREAL4Vector(usableTFdata);
         XLALDestroyREAL4Vector(background);
      }

      //Print out data product if requested
      if (uvar.printData) {
         CHARVector *outputfile = NULL;
         XLAL_CHECK( (outputfile = XLALCreateCHARVector(strlen(uvar.outdirectory)+20)) != NULL, XLAL_EFUNC );
         sprintf(outputfile->data, "%s/%s", uvar.outdirectory, "tfbackground.dat");
         XLAL_CHECK( printREAL4Vector2File(background_slided, outputfile->data) == XLAL_SUCCESS, XLAL_EFUNC );
         memset(outputfile->data, 0, outputfile->length*sizeof(CHAR));
         sprintf(outputfile->data, "%s/%s", uvar.outdirectory, "tfdata_slided.dat");
         XLAL_CHECK( printREAL4Vector2File(TFdata_slided, outputfile->data) == XLAL_SUCCESS, XLAL_EFUNC );
         XLALDestroyCHARVector(outputfile);
      }

      //Check the RMS of the antenna weights; if a deviation of more than 1%, then reset the IHS FAR
      if (antweightsrms == 0.0) antweightsrms = currentAntWeightsRMS;
      if ( fabs(currentAntWeightsRMS-antweightsrms)/antweightsrms >= 0.01 ) {
         ihsfarstruct->ihsfar->data[0] = -1.0;
         antweightsrms = currentAntWeightsRMS;
      }

      //Antenna normalization for different sky locations
      REAL8 skypointffnormalization = 1.0;
      XLAL_CHECK( ffPlaneNoise(aveNoise, &uvar, sftexist, background_slided, (REAL4Vector*)antweightsforihs2h0, secondFFTplan, &(skypointffnormalization), rng) == XLAL_SUCCESS, XLAL_EFUNC );

      //Average noise floor of FF plane for each 1st FFT frequency bin
      ffdata->ffnormalization = 1.0;
      XLAL_CHECK( ffPlaneNoise(aveNoise, &uvar, sftexist, background_slided, (REAL4Vector*)antweights, secondFFTplan, &(ffdata->ffnormalization), rng) == XLAL_SUCCESS, XLAL_EFUNC );
      if (uvar.printData) {
         CHARVector *outputfile = NULL;
         XLAL_CHECK( (outputfile = XLALCreateCHARVector(strlen(uvar.outdirectory)+20)) != NULL, XLAL_EFUNC );
         sprintf(outputfile->data, "%s/%s", uvar.outdirectory, "ffbackground.dat");
         XLAL_CHECK( printREAL4Vector2File((REAL4Vector*)aveNoise, outputfile->data) == XLAL_SUCCESS, XLAL_EFUNC );
         XLALDestroyCHARVector(outputfile);
      }

      //Compute the weighted TF data
      REAL4Vector *TFdata_weighted = NULL;
      XLAL_CHECK( (TFdata_weighted = XLALCreateREAL4Vector(ffdata->numffts*ffdata->numfbins)) != NULL, XLAL_EFUNC );
      XLAL_CHECK( tfWeight(TFdata_weighted, TFdata_slided, background_slided, antweights, indexValuesOfExistingSFTs, &uvar) == XLAL_SUCCESS, XLAL_EFUNC );
      XLALDestroyREAL4Vector(TFdata_slided);
      XLALDestroyREAL4VectorAligned(antweights);

      //Print out data product if requested
      if (uvar.printData) {
         CHARVector *outputfile = NULL;
         XLAL_CHECK( (outputfile = XLALCreateCHARVector(strlen(uvar.outdirectory)+20)) != NULL, XLAL_EFUNC );
         sprintf(outputfile->data, "%s/%s", uvar.outdirectory, "procTFdata.dat");
         XLAL_CHECK( printREAL4Vector2File(TFdata_weighted, outputfile->data) == XLAL_SUCCESS, XLAL_EFUNC );
         XLALDestroyCHARVector(outputfile);
      }

      //Calculation of average TF noise per frequency bin ratio to total mean
      //this block of code does not avoid lines when computing the average F-bin ratio. Upper limits remain virtually unchanged
      //when comaring runs that have line finding enabled or disabled
      REAL4Vector *aveTFnoisePerFbinRatio = NULL, *TSofPowers = NULL;
      XLAL_CHECK( (aveTFnoisePerFbinRatio = XLALCreateREAL4Vector(ffdata->numfbins)) != NULL, XLAL_EFUNC );
      XLAL_CHECK( (TSofPowers = XLALCreateREAL4Vector(ffdata->numffts)) != NULL, XLAL_EFUNC );
      memset(TSofPowers->data, 0, sizeof(REAL4)*TSofPowers->length);
      for (ii=0; ii<ffdata->numfbins; ii++) {
         REAL4 totalweightval = 0.0;
         for (jj=0; jj<ffdata->numffts; jj++) {
            if (background_slided->data[jj*ffdata->numfbins + ii]!=0.0) {
               TSofPowers->data[jj] = 1.0/(background_slided->data[jj*ffdata->numfbins + ii]);
               totalweightval += 1.0/(background_slided->data[jj*ffdata->numfbins + ii]*background_slided->data[jj*ffdata->numfbins + ii]);
            }
         }
         aveTFnoisePerFbinRatio->data[ii] = 0.0;
         for (jj=0; jj<ffdata->numffts; jj++) aveTFnoisePerFbinRatio->data[ii] += TSofPowers->data[jj];
         aveTFnoisePerFbinRatio->data[ii] = (aveTFnoisePerFbinRatio->data[ii]/totalweightval);
      }
      REAL4 aveTFaveinv = 1.0/calcMean(aveTFnoisePerFbinRatio);
      for (ii=0; ii<ffdata->numfbins; ii++) {
         aveTFnoisePerFbinRatio->data[ii] *= aveTFaveinv;
         //fprintf(stderr, "%f\n", aveTFnoisePerFbinRatio->data[ii]);
      }
      XLALDestroyREAL4Vector(TSofPowers);
      XLALDestroyREAL4Vector(background_slided);

      //Do the second FFT
      XLAL_CHECK( makeSecondFFT(ffdata, TFdata_weighted, secondFFTplan) == XLAL_SUCCESS, XLAL_EFUNC );
      //Normalize according to LAL PSD spec (also done in ffPlaneNoise() so this doesn't change anything)
      //There is a secret divide by numffts in the weighting of the TF data (sumofweights), so we don't need to do it here;
      //the numffts divisor gets squared when taking the PSD, so it is not applied here
      for (ii=0; ii<(INT4)ffdata->ffdata->length; ii++) ffdata->ffdata->data[ii] *= uvar.Tobs;

      REAL4 secFFTmean = 0.0, secFFTsigma = 0.0;
      secFFTmean = calcMean(ffdata->ffdata);
      XLAL_CHECK( calcStddev(&secFFTsigma, ffdata->ffdata) == XLAL_SUCCESS, XLAL_EFUNC );

      XLALDestroyREAL4Vector(TFdata_weighted);
      TFdata_weighted = NULL;

      //Print out data product if requested
      if (uvar.printData) {
         CHARVector *outputfile = NULL;
         XLAL_CHECK( (outputfile = XLALCreateCHARVector(strlen(uvar.outdirectory)+10)) != NULL, XLAL_EFUNC );
         sprintf(outputfile->data, "%s/%s", uvar.outdirectory, "ffdata.dat");
         XLAL_CHECK( printREAL4Vector2File(ffdata->ffdata, outputfile->data) == XLAL_SUCCESS, XLAL_EFUNC );
         XLALDestroyCHARVector(outputfile);
      }

      fprintf(LOG, "2nd FFT ave = %g, 2nd FFT stddev = %g, expected ave = 1.0\n", secFFTmean, secFFTsigma);
      fprintf(stderr, "2nd FFT ave = %g, 2nd FFT stddev = %g, expected ave = 1.0\n", secFFTmean, secFFTsigma);

      //Exit with failure if there are no SFTs (probably this doesn't get hit)
      XLAL_CHECK( secFFTmean != 0.0, XLAL_FAILURE, "Average second FFT power is 0.0. Perhaps no SFTs are remaining? Program exiting with failure.\n" );

      //If the user wants to test a single, exact template, then we do that here
      if (uvar.templateTest) {
         if (uvar.printData) fprintf(stderr, "numfbins=%d, maxbinshift=%d, numffts=%d, numfprbins=%d\n", ffdata->numfbins, maxbinshift, ffdata->numffts, ffdata->numfprbins);
         fprintf(stderr, "Testing template f=%f, P=%f, Df=%f\n", uvar.templateTestF, uvar.templateTestP, uvar.templateTestDf);
         fprintf(LOG, "Testing template f=%f, P=%f, Df=%f\n", uvar.templateTestF, uvar.templateTestP, uvar.templateTestDf);

         //Load template quantities into a test candidate
         loadCandidateData(&(exactCandidates1->data[0]), uvar.templateTestF, uvar.templateTestP, uvar.templateTestDf, dopplerpos.Alpha, dopplerpos.Delta, 0.0, 0.0, 0.0, 0, 0.0, -1);

         //Resize the output candidate vector if necessary
         if (exactCandidates2->numofcandidates == exactCandidates2->length-1) XLAL_CHECK( (exactCandidates2 = resize_candidateVector(exactCandidates2, 2*exactCandidates2->length)) != NULL, XLAL_EFUNC );

         //Analyze the template stored in the test candidate
         XLAL_CHECK( analyzeOneTemplate(&(exactCandidates2->data[exactCandidates2->numofcandidates]), &(exactCandidates1->data[0]), ffdata, (REAL4Vector*)aveNoise, aveTFnoisePerFbinRatio, &uvar, sftexist, secondFFTplan, rng, !uvar.gaussTemplatesOnly) == XLAL_SUCCESS, XLAL_EFUNC );
         (exactCandidates2->numofcandidates)++;

         //Rescale the h0 output from the normaliztions and amount of observation time present
         exactCandidates2->data[exactCandidates2->numofcandidates-1].h0 /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25);

      } else if (uvar.bruteForceTemplateTest) {
         candidate cand;
         loadCandidateData(&cand, uvar.templateTestF, uvar.templateTestP, uvar.templateTestDf, dopplerpos.Alpha, dopplerpos.Delta, 0.0, 0.0, 0.0, 0, 0.0, -1);
         TwoSpectParamSpaceSearchVals paramspace = {uvar.templateTestF-2.0/uvar.Tsft, uvar.templateTestF+2.0/uvar.Tsft, 21, 5, 5, 0.5, uvar.templateTestDf-2.0/uvar.Tsft,
                                                    uvar.templateTestDf+2.0/uvar.Tsft, 21};
         XLAL_CHECK( bruteForceTemplateTest(&(exactCandidates2), cand, &paramspace, &uvar, ffdata->ffdata, sftexist, (REAL4Vector*)aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, rng, 1) == XLAL_SUCCESS, XLAL_EFUNC );
         for (ii=0; ii<(INT4)exactCandidates2->numofcandidates; ii++) exactCandidates2->data[ii].h0 /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25);
      }

      //If the user wants to do a template search, that is done here
      if (uvar.templateSearch) {
         //fprintf(stderr, "Calling templateSearch\n (in development, last edited 2014-06-09)\n");
         XLAL_CHECK( templateSearch_scox1Style(&exactCandidates2, uvar.fmin, uvar.fspan, uvar.templateSearchP, uvar.templateSearchAsini, uvar.templateSearchAsiniSigma, skypos, &uvar, ffdata->ffdata, sftexist, (REAL4Vector*)aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, rng, 1) == XLAL_SUCCESS, XLAL_EFUNC );
         //fprintf(stderr, "Done calling templateSearch\n");
         for (ii=0; ii<(INT4)exactCandidates2->numofcandidates; ii++) exactCandidates2->data[ii].h0 /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25);
      }

      //Template bank analysis
      if (XLALUserVarWasSet(&(uvar.templatebankfile))) {
         INT4 numofnewcandidates = 0;
         XLAL_CHECK( (numofnewcandidates = testTwoSpectTemplateVector(&exactCandidates2, templateVec, ffdata, (REAL4Vector*)aveNoise, aveTFnoisePerFbinRatio, skypos, &uvar, rng)) >= 0, XLAL_EFUNC );
         for (ii=0; ii<numofnewcandidates; ii++) {
            /* TwoSpectParamSpaceSearchVals paramspace = {exactCandidates2->data[exactCandidates2->numofcandidates-ii-1].fsig-0.5/uvar.Tsft,
                                                       exactCandidates2->data[exactCandidates2->numofcandidates-ii-1].fsig+0.5/uvar.Tsft, 9, 4, 4, 0.4,
                                                       exactCandidates2->data[exactCandidates2->numofcandidates-ii-1].moddepth-0.5/uvar.Tsft,
                                                       exactCandidates2->data[exactCandidates2->numofcandidates-ii-1].moddepth+0.5/uvar.Tsft, 9};
            XLAL_CHECK( bruteForceTemplateSearch(&(exactCandidates2->data[exactCandidates2->numofcandidates-ii-1]),
                                                 exactCandidates2->data[exactCandidates2->numofcandidates-ii-1], &paramspace, &uvar, ffdata->ffdata, sftexist,
                                                 (REAL4Vector*)aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, rng, 1) == XLAL_SUCCESS, XLAL_EFUNC ); */
            exactCandidates2->data[exactCandidates2->numofcandidates-ii-1].h0 /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25); //Scaling here
         }
      }

      if (uvar.signalOnly) return 0;

      //Start of the IHS step!
      //Find the FAR of IHS sum -- only if the templateTest has not been given
      candidateVector *ihsCandidates_reduced = NULL;
      if (!XLALUserVarWasSet(&uvar.templatebankfile) && !uvar.templateTest && !uvar.templateSearch && !uvar.bruteForceTemplateTest) {
         //If the false alarm thresholds need to be computed
         if (ihsfarstruct->ihsfar->data[0]<0.0) XLAL_CHECK( genIhsFar(ihsfarstruct, &uvar, maxrows, (REAL4Vector*)aveNoise, rng) == XLAL_SUCCESS, XLAL_EFUNC );

         //Run the IHS algorithm on the data
         XLAL_CHECK( runIHS(ihsmaxima, ffdata, ihsfarstruct, &uvar, maxrows, (REAL4Vector*)aveNoise, aveTFnoisePerFbinRatio) == XLAL_SUCCESS, XLAL_EFUNC );

         //Find any IHS candidates
         XLAL_CHECK( findIHScandidates(&ihsCandidates, ihsfarstruct, &uvar, ffdata, ihsmaxima, aveTFnoisePerFbinRatio, trackedlines) == XLAL_SUCCESS, XLAL_EFUNC );

         //If requested, keep only the most significant IHS candidates
         if (XLALUserVarWasSet(&uvar.keepOnlyTopNumIHS) && (INT4)ihsCandidates->numofcandidates>uvar.keepOnlyTopNumIHS) {
            XLAL_CHECK( (ihsCandidates_reduced = keepMostSignificantCandidates(ihsCandidates, &uvar)) != NULL, XLAL_EFUNC );
            //Put ihsCandidates_reduced back into a reset ihsCandidates
            ihsCandidates->numofcandidates = 0;
            for (ii=0; ii<(INT4)ihsCandidates_reduced->numofcandidates; ii++) {
               loadCandidateData(&(ihsCandidates->data[ii]), ihsCandidates_reduced->data[ii].fsig, ihsCandidates_reduced->data[ii].period, ihsCandidates_reduced->data[ii].moddepth, dopplerpos.Alpha, dopplerpos.Delta, ihsCandidates_reduced->data[ii].stat, ihsCandidates_reduced->data[ii].h0, ihsCandidates_reduced->data[ii].prob, 0, ihsCandidates_reduced->data[ii].normalization, ihsCandidates_reduced->data[ii].templateVectorIndex);
               (ihsCandidates->numofcandidates)++;
            }
            free_candidateVector(ihsCandidates_reduced);
         }
      }
      //End of the IHS step

      //Start of the Gaussian template search!
      //First check to see if the IHSonly or templateTest or templateSearch was given
      if (uvar.IHSonly && !uvar.templateTest && !uvar.templateSearch && !XLALUserVarWasSet(&uvar.templatebankfile)) {
         //Check the length of the exactCandidates2 vector is large enough and resize if necessary
         if (exactCandidates2->length < exactCandidates2->numofcandidates+ihsCandidates->numofcandidates) XLAL_CHECK( (exactCandidates2 = resize_candidateVector(exactCandidates2, exactCandidates2->numofcandidates+ihsCandidates->numofcandidates)) != NULL, XLAL_EFUNC );

         //Use the typical list
         INT4 numofcandidatesalready = exactCandidates2->numofcandidates;
         for (ii=0; ii<(INT4)ihsCandidates->numofcandidates; ii++) {
            loadCandidateData(&(exactCandidates2->data[ii+numofcandidatesalready]), ihsCandidates->data[ii].fsig, ihsCandidates->data[ii].period, ihsCandidates->data[ii].moddepth, dopplerpos.Alpha, dopplerpos.Delta, ihsCandidates->data[ii].stat, ihsCandidates->data[ii].h0, ihsCandidates->data[ii].prob, 0, ihsCandidates->data[ii].normalization, ihsCandidates->data[ii].templateVectorIndex);
            exactCandidates2->data[ii+numofcandidatesalready].h0 /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25); //Scaling here
            (exactCandidates2->numofcandidates)++;
         }

      } else if (!uvar.templateTest && !uvar.templateSearch && !XLALUserVarWasSet(&uvar.templatebankfile)) {

         //Test the IHS candidates against Gaussian templates in this function
         XLAL_CHECK( testIHScandidates(&gaussCandidates1, ihsCandidates, ffdata, (REAL4Vector*)aveNoise, aveTFnoisePerFbinRatio, skypos, &uvar, rng) == XLAL_SUCCESS, XLAL_EFUNC );

         for (ii=0; ii<(INT4)gaussCandidates1->numofcandidates; ii++) fprintf(stderr, "Candidate %d: f0=%g, P=%g, df=%g\n", ii, gaussCandidates1->data[ii].fsig, gaussCandidates1->data[ii].period, gaussCandidates1->data[ii].moddepth);
      } /* if IHSonly is not given && templateTest not given and templateSearch not given */
      //End of the Gaussian template search

      //Reset IHS candidates, but keep length the same (doesn't reset actual values in the vector)
      ihsCandidates->numofcandidates = 0;

      //Search the candidates further if the number of candidates passing the first Gaussian template test is greater than 0
      if (gaussCandidates1->numofcandidates>0) {
         //Start clustering! Note that the clustering algorithm takes care of the period range of parameter space
         XLAL_CHECK( clusterCandidates(&gaussCandidates2, gaussCandidates1, ffdata, &uvar, (REAL4Vector*)aveNoise, aveTFnoisePerFbinRatio, sftexist, rng, 0) == XLAL_SUCCESS, XLAL_EFUNC );

         for (ii=0; ii<(INT4)gaussCandidates2->numofcandidates; ii++) fprintf(stderr, "Candidate %d: f0=%g, P=%g, df=%g\n", ii, gaussCandidates2->data[ii].fsig, gaussCandidates2->data[ii].period, gaussCandidates2->data[ii].moddepth);

         //Reset first set of Gaussian template candidates, but keep length the same (doesn't reset actual values in the vector)
         gaussCandidates1->numofcandidates = 0;

         //Start detailed Gaussian template search!
         for (ii=0; ii<(INT4)gaussCandidates2->numofcandidates; ii++) {

            if (gaussCandidates3->numofcandidates == gaussCandidates3->length-1) XLAL_CHECK( (gaussCandidates3 = resize_candidateVector(gaussCandidates3, 2*gaussCandidates3->length)) != NULL, XLAL_EFUNC );

            TwoSpectParamSpaceSearchVals paramspace = {gaussCandidates2->data[ii].fsig-1.0/uvar.Tsft, gaussCandidates2->data[ii].fsig+1.0/uvar.Tsft,5, 2, 2, 1.0,
                                                       gaussCandidates2->data[ii].moddepth-1.0/uvar.Tsft, gaussCandidates2->data[ii].moddepth+1.0/uvar.Tsft, 5};
            XLAL_CHECK( bruteForceTemplateSearch(&(gaussCandidates3->data[gaussCandidates3->numofcandidates]), gaussCandidates2->data[ii], &paramspace, &uvar, ffdata->ffdata, 
                                                 sftexist, (REAL4Vector*)aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, rng, 0) == XLAL_SUCCESS, XLAL_EFUNC );
            gaussCandidates3->numofcandidates++;

         } /* for ii < numofcandidates */

         for (ii=0; ii<(INT4)gaussCandidates3->numofcandidates; ii++) fprintf(stderr,"Candidate %d: f0=%g, P=%g, df=%g\n", ii, gaussCandidates3->data[ii].fsig, gaussCandidates3->data[ii].period, gaussCandidates3->data[ii].moddepth);

         //Reset 2nd round of Gaussian template candidates, but keep length the same (doesn't reset actual values in the vector)
         gaussCandidates2->numofcandidates = 0;

         //Start clustering!
         XLAL_CHECK( clusterCandidates(&gaussCandidates4, gaussCandidates3, ffdata, &uvar, (REAL4Vector*)aveNoise, aveTFnoisePerFbinRatio, sftexist, rng, 0) == XLAL_SUCCESS, XLAL_EFUNC );

         for (ii=0; ii<(INT4)gaussCandidates4->numofcandidates; ii++) fprintf(stderr, "Candidate %d: f0=%g, P=%g, df=%g\n", ii, gaussCandidates4->data[ii].fsig, gaussCandidates4->data[ii].period, gaussCandidates4->data[ii].moddepth);

         //Reset 3rd set of Gaussian template candidates
         gaussCandidates3->numofcandidates = 0;

         //Initial check using "exact" template
         for (ii=0; ii<(INT4)gaussCandidates4->numofcandidates; ii++) {
            //Allocate the template memory
            TwoSpectTemplate *template = NULL;
            XLAL_CHECK( (template = new_TwoSpectTemplate(uvar.maxTemplateLength)) != NULL, XLAL_EFUNC );

            //Produce either the Gaussian template or exact template
            if (!uvar.gaussTemplatesOnly) XLAL_CHECK( makeTemplate(template, gaussCandidates4->data[ii], &uvar, sftexist, secondFFTplan) == XLAL_SUCCESS, XLAL_EFUNC );
            else XLAL_CHECK( makeTemplateGaussians(template, gaussCandidates4->data[ii], &uvar) == XLAL_SUCCESS, XLAL_EFUNC );

            //Produce the FAR threshold, if requested
            farStruct *farval = NULL;
            if (uvar.calcRthreshold) {
               XLAL_CHECK( (farval = new_farStruct()) != NULL, XLAL_EFUNC );
               XLAL_CHECK( numericFAR(farval, template, uvar.tmplfar, (REAL4Vector*)aveNoise, aveTFnoisePerFbinRatio, &uvar, rng, uvar.BrentsMethod) == XLAL_SUCCESS, XLAL_EFUNC );
            }

            //Calculate R, false alarm probability and h0 estimate
            REAL8 R = calculateR(ffdata->ffdata, template, (REAL4Vector*)aveNoise, aveTFnoisePerFbinRatio);
            XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
            REAL8 h0 = 0.0;
            REAL8 prob = probR(template, (REAL4Vector*)aveNoise, aveTFnoisePerFbinRatio, R, &uvar, rng, &proberrcode);
            XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
            if ( R > 0.0 ) h0 = 2.7426*pow(R/(uvar.Tsft*uvar.Tobs),0.25);

            //Check that we are above threshold before storing the candidate
            if ((!uvar.calcRthreshold && prob<log10(uvar.tmplfar)) || (uvar.calcRthreshold && R>farval->far)) {
               if (exactCandidates1->numofcandidates == exactCandidates1->length-1) XLAL_CHECK( (exactCandidates1 = resize_candidateVector(exactCandidates1, 2*exactCandidates1->length)) != NULL, XLAL_EFUNC );
               loadCandidateData(&exactCandidates1->data[exactCandidates1->numofcandidates], gaussCandidates4->data[ii].fsig, gaussCandidates4->data[ii].period, gaussCandidates4->data[ii].moddepth, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, R, h0, prob, proberrcode, gaussCandidates4->data[ii].normalization, gaussCandidates4->data[ii].templateVectorIndex);
               exactCandidates1->numofcandidates++;
            }

            free_TwoSpectTemplate(template);
            template = NULL;
            if (uvar.calcRthreshold) {
               free_farStruct(farval);
               farval = NULL;
            }
         } /* for ii < numofcandidates */
         fprintf(LOG, "Number of candidates confirmed with exact templates = %d\n", exactCandidates1->numofcandidates);
         fprintf(stderr, "Number of candidates confirmed with exact templates = %d\n", exactCandidates1->numofcandidates);
         for (ii=0; ii<(INT4)exactCandidates1->numofcandidates; ii++) fprintf(stderr, "Candidate %d: f0=%g, P=%g, df=%g\n", ii, exactCandidates1->data[ii].fsig, exactCandidates1->data[ii].period, exactCandidates1->data[ii].moddepth);
         //Done with initial check

         //Reset 4th set of Gaussian template candidates, but keep length the same (doesn't reset actual values in the vector)
         gaussCandidates4->numofcandidates = 0;

         //Start detailed "exact" template search!
         for (ii=0; ii<(INT4)exactCandidates1->numofcandidates; ii++) {

            if (exactCandidates2->numofcandidates == exactCandidates2->length-1) XLAL_CHECK( (exactCandidates2 = resize_candidateVector(exactCandidates2, 2*exactCandidates2->length)) != NULL, XLAL_EFUNC );

            TwoSpectParamSpaceSearchVals paramspace = {exactCandidates1->data[ii].fsig-0.5/uvar.Tsft, exactCandidates1->data[ii].fsig+0.5/uvar.Tsft, 3, 1, 1, 1.0,
                                                       exactCandidates1->data[ii].moddepth-0.5/uvar.Tsft, exactCandidates1->data[ii].moddepth+0.5/uvar.Tsft, 3};
            if (!uvar.gaussTemplatesOnly) XLAL_CHECK( bruteForceTemplateSearch(&(exactCandidates2->data[exactCandidates2->numofcandidates]), exactCandidates1->data[ii], &paramspace,&uvar, ffdata->ffdata, sftexist, (REAL4Vector*)aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, rng, 1) == XLAL_SUCCESS, XLAL_EFUNC );
            else XLAL_CHECK( bruteForceTemplateSearch(&(exactCandidates2->data[exactCandidates2->numofcandidates]), exactCandidates1->data[ii], &paramspace, &uvar, ffdata->ffdata, sftexist, (REAL4Vector*)aveNoise, aveTFnoisePerFbinRatio, secondFFTplan,rng, 0) == XLAL_SUCCESS, XLAL_EFUNC );
            exactCandidates2->data[exactCandidates2->numofcandidates].h0 /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25);  //Scaling here
            exactCandidates2->numofcandidates++;

            fprintf(stderr,"Candidate %d: f0=%g, P=%g, df=%g\n", ii, exactCandidates2->data[ii].fsig, exactCandidates2->data[ii].period, exactCandidates2->data[ii].moddepth);
         } /* for ii < numofcandidates */
         //End of detailed search

         //Reset first round of exact template candidates, but keep length the same (doesn't reset actual values in the vector)
         exactCandidates1->numofcandidates = 0;

         fprintf(LOG,"Exact step is done with the total number of candidates = %d\n", exactCandidates2->numofcandidates);
         fprintf(stderr,"Exact step is done with the total number of candidates = %d\n", exactCandidates2->numofcandidates);

      } /* if gaussCandidates1->numofcandidates > 0 */

      //Determine upper limits, if the ULoff has not been set
      if (!uvar.ULoff && !uvar.templateTest && !uvar.templateSearch && !XLALUserVarWasSet(&uvar.templatebankfile)) {
         upperlimits->data[upperlimits->length-1].alpha = (REAL4)dopplerpos.Alpha;
         upperlimits->data[upperlimits->length-1].delta = (REAL4)dopplerpos.Delta;
         upperlimits->data[upperlimits->length-1].normalization = ffdata->tfnormalization;

         XLAL_CHECK( skypoint95UL(&(upperlimits->data[upperlimits->length-1]), &uvar, ffdata, ihsmaxima, ihsfarstruct, aveTFnoisePerFbinRatio) == XLAL_SUCCESS, XLAL_EFUNC );

         for (ii=0; ii<(INT4)upperlimits->data[upperlimits->length-1].ULval->length; ii++) upperlimits->data[upperlimits->length-1].ULval->data[ii] /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25);  //Compensation for different duty cycle and antenna pattern weights

         XLAL_CHECK( (upperlimits = resize_UpperLimitVector(upperlimits, upperlimits->length+1)) != NULL, XLAL_EFUNC );
      } /* if producing UL */

      //Destroy stuff
      XLALDestroyREAL4Vector(aveTFnoisePerFbinRatio);
      XLALDestroyREAL4VectorSequence(trackedlines);
      if (detectors->length>1 && !uvar.signalOnly) {
         XLALDestroyINT4Vector(sftexist);
         XLALDestroyINT4Vector(indexValuesOfExistingSFTs);
      }

      //Iterate to next sky location
      XLAL_CHECK( XLALNextDopplerSkyPos(&dopplerpos, &scan) == XLAL_SUCCESS, XLAL_EFUNC );

   } /* while sky scan is not finished */

   if (exactCandidates2->numofcandidates!=0) {
      fprintf(LOG, "\n**Report of candidates:**\n");
      fprintf(stderr, "\n**Report of candidates:**\n");

      for (ii=0; ii<(INT4)exactCandidates2->numofcandidates; ii++) {
         fprintf(LOG, "fsig = %.6f, period = %.6f, df = %.7f, RA = %.4f, DEC = %.4f, R = %.4f, h0 = %g, Prob = %.4f, TF norm = %g\n", exactCandidates2->data[ii].fsig, exactCandidates2->data[ii].period, exactCandidates2->data[ii].moddepth, exactCandidates2->data[ii].ra, exactCandidates2->data[ii].dec, exactCandidates2->data[ii].stat, exactCandidates2->data[ii].h0, exactCandidates2->data[ii].prob, ffdata->tfnormalization);
         fprintf(stderr, "fsig = %.6f, period = %.6f, df = %.7f, RA = %.4f, DEC = %.4f, R = %.4f, h0 = %g, Prob = %.4f, TF norm = %g\n", exactCandidates2->data[ii].fsig, exactCandidates2->data[ii].period, exactCandidates2->data[ii].moddepth, exactCandidates2->data[ii].ra, exactCandidates2->data[ii].dec, exactCandidates2->data[ii].stat, exactCandidates2->data[ii].h0, exactCandidates2->data[ii].prob, ffdata->tfnormalization);
      } /* for ii < exactCandidates2->numofcandidates */
   } /* if exactCandidates2->numofcandidates != 0 */

   //Output upper limits to a file, if ULoff is not given
   if (!uvar.ULoff) {
      for (ii=0; ii<(INT4)upperlimits->length-1; ii++) {
         XLAL_CHECK( outputUpperLimitToFile(ULFILENAME->data, upperlimits->data[ii], uvar.allULvalsPerSkyLoc) == XLAL_SUCCESS, XLAL_EFUNC );
      }
   }

   //Destroy varaibles
   if (detectors->length>1) {
      XLALDestroyMultiSFTVector(multiSFTvector);
   } else {
      XLALDestroyREAL4Vector(background);
      XLALDestroyREAL4Vector(usableTFdata);
      if (!uvar.signalOnly) {
         XLALDestroyINT4Vector(sftexist);
         XLALDestroyINT4Vector(indexValuesOfExistingSFTs);
      }
   }
   if (XLALUserVarWasSet(&(uvar.templatebankfile))) free_TwoSpectTemplateVector(templateVec);
   XLALFree(detectors);
   gsl_rng_free(rng);
   XLALDestroyREAL4VectorAligned(antweightsforihs2h0);
   XLALDestroyREAL4VectorAligned(aveNoise);
   XLALDestroyINT4Vector(binshifts);
   XLALDestroyMultiTimestamps(multiTimestamps);
   XLALDestroyMultiDetectorStateSeries(multiStateSeries);
   free_ffdata(ffdata);
   free_ihsfarStruct(ihsfarstruct);
   free_ihsMaxima(ihsmaxima);
   XLALDestroyREAL4FFTPlan(secondFFTplan);
   XLALDestroyEphemerisData(edat);
   XLALDestroyUserVars();
   XLALDestroyCHARVector(ULFILENAME);
   XLALDestroyCHARVector(LOGFILENAME);
   free_UpperLimitVector(upperlimits);
   free_candidateVector(ihsCandidates);
   free_candidateVector(gaussCandidates1);
   free_candidateVector(gaussCandidates2);
   free_candidateVector(gaussCandidates3);
   free_candidateVector(gaussCandidates4);
   free_candidateVector(exactCandidates1);
   free_candidateVector(exactCandidates2);
   FreeDopplerSkyScan(&status, &scan);
   if (lines!=NULL) XLALDestroyINT4Vector(lines);

   //print end time
   time(&programendtime);
   ptm = localtime(&programendtime);
   fprintf(stderr, "Program finished on %s", asctime(ptm));
   fprintf(LOG, "Program finished on %s", asctime(ptm));

   fclose(LOG);

   //Check for leaks
   LALCheckMemoryLeaks();

   //The end!
   return 0;

} /* main() */


/**
 * Create a new frequency-frequency data structure for the TwoSpect analysis
 * \param [in] params Pointer to the UserInput_t
 * \return Pointer to new ffdataStruct
 */
ffdataStruct * new_ffdata(UserInput_t *params)
{

   XLAL_CHECK_NULL( params != NULL, XLAL_EINVAL );

   ffdataStruct *ffdata = NULL;
   XLAL_CHECK_NULL( (ffdata = XLALMalloc(sizeof(*ffdata))) != NULL, XLAL_ENOMEM );

   ffdata->numfbins = (INT4)(round(params->fspan*params->Tsft + 2.0*params->dfmax*params->Tsft)+12+1);
   ffdata->numffts = (INT4)floor(params->Tobs/(params->Tsft-params->SFToverlap)-1);
   ffdata->numfprbins = (INT4)floorf(ffdata->numffts*0.5) + 1;
   ffdata->tfnormalization = ffdata->ffnormalization = 0.0;

   XLAL_CHECK_NULL ( (ffdata->ffdata = XLALCreateREAL4Vector(ffdata->numfbins * ffdata->numfprbins)) != NULL, XLAL_EFUNC );

   return ffdata;

} /* new_ffdata() */


/**
 * Free the frequency-frequency data structure
 * \param [in] data Pointer to the ffdataStruct
 */
void free_ffdata(ffdataStruct *data)
{

   XLALDestroyREAL4Vector(data->ffdata);
   XLALFree((ffdataStruct*)data);

} /* free_ffdata() */


/**
 * Line detection algorithm
 *
 * 1. Calculate RMS for each fbin as function of time
 * 2. Calculate running median of RMS values
 * 3. Divide RMS values by running median. Lines stick out by being >>1
 * 4. Add up number of lines above threshold and report
 *
 * \param [in] TFdata Pointer to REAL4Vector of SFT powers
 * \param [in] ffdata Pointer to ffdataStruct
 * \param [in] params Pointer to UserInput_t
 * \return Pointer to INT4Vector of bin number of the lines
 */
INT4Vector * detectLines_simple(REAL4Vector *TFdata, ffdataStruct *ffdata, UserInput_t *params)
{

   XLAL_CHECK_NULL( TFdata != NULL && ffdata != NULL && params != NULL, XLAL_EINVAL );

   fprintf(stderr, "Running line detection algorithm... ");

   LALStatus XLAL_INIT_DECL(status);
   INT4 blksize = 11;
   INT4 numlines = 0;
   INT4Vector *lines = NULL;
   //INT4 totalnumfbins = ffdata->numfbins+(params->blksize-1)+2*params->maxbinshift;
   INT4 numffts = (INT4)floor(params->Tobs/(params->Tsft-params->SFToverlap)-1);    //Number of FFTs
   INT4 totalnumfbins = (INT4)TFdata->length/numffts;

   //Blocksize of running median
   LALRunningMedianPar block = {blksize};

   //Compute weights
   REAL4Vector *sftdata = NULL, *weights = NULL;
   XLAL_CHECK_NULL( (sftdata = XLALCreateREAL4Vector(totalnumfbins)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (weights = XLALCreateREAL4Vector(ffdata->numffts)) != NULL, XLAL_EFUNC );
   memset(weights->data, 0, ffdata->numffts*sizeof(REAL4));
   REAL4 sumweights = 0.0;
   for (INT4 ii=0; ii<ffdata->numffts; ii++) {
      if (TFdata->data[ii*totalnumfbins]!=0.0) {
         memcpy(sftdata->data, &(TFdata->data[ii*totalnumfbins]), totalnumfbins*sizeof(REAL4));
         REAL4 stddev = 0.0;
         XLAL_CHECK_NULL( calcStddev(&stddev, sftdata) == XLAL_SUCCESS, XLAL_EFUNC );
         weights->data[ii] = 1.0/(stddev*stddev);
         sumweights += weights->data[ii];
      }
   }
   REAL4 invsumweights = 1.0/sumweights;
   XLALDestroyREAL4Vector(sftdata);

   //Compute RMS for each frequency bin as a function of time
   REAL4Vector *testRMSvals = NULL, *testRngMedian = NULL, *testTSofPowers = NULL;
   XLAL_CHECK_NULL( (testRMSvals = XLALCreateREAL4Vector(totalnumfbins)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (testRngMedian = XLALCreateREAL4Vector(totalnumfbins-(blksize-1))) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (testTSofPowers = XLALCreateREAL4Vector(ffdata->numffts)) != NULL, XLAL_EFUNC );

   for (UINT4 ii=0; ii<testRMSvals->length; ii++) {
      for (INT4 jj=0; jj<ffdata->numffts; jj++) testTSofPowers->data[jj] = TFdata->data[jj*testRMSvals->length + ii]*weights->data[jj]*invsumweights;
      XLAL_CHECK_NULL( calcRms(&(testRMSvals->data[ii]), testTSofPowers) == XLAL_SUCCESS, XLAL_EFUNC ); //This approaches calcMean(TSofPowers) for stationary noise
   }

   //Running median of RMS values
   LALSRunningMedian2(&status, testRngMedian, testRMSvals, block);
   XLAL_CHECK_NULL( status.statusCode == 0, XLAL_EFUNC );

   //Determine which bins are above the threshold and store the bin number of the line
   for (UINT4 ii=0; ii<testRngMedian->length; ii++) {
      REAL4 normrmsval = testRMSvals->data[ii+(blksize-1)/2]/testRngMedian->data[ii];
      if ( (INT4)(ii+(blksize-1)/2) > ((params->blksize-1)/2) && normrmsval > params->lineDetection) {
         XLAL_CHECK_NULL( (lines = XLALResizeINT4Vector(lines, numlines+1)) != NULL, XLAL_EFUNC );
         lines->data[numlines] = ii+(blksize-1)/2;
         numlines++;
      }
   }

   //Destroy stuff
   XLALDestroyREAL4Vector(testTSofPowers);
   XLALDestroyREAL4Vector(testRngMedian);
   XLALDestroyREAL4Vector(testRMSvals);
   XLALDestroyREAL4Vector(weights);

   fprintf(stderr, "done\n");

   return lines;

}

/**
 * Track the lines for the sky position
 * \param [in] lines     Pointer to INT4Vector with list of lines
 * \param [in] binshifts Pointer to INT4Vector of SFT bin shifts
 * \param [in] minfbin   Frequency value of the lowest frequency bin
 * \param [in] df        Spacing of frequency bins (typically 1/Tsft)
 * \return Pointer to REAL4VectorSequence containing lower/mid/upper frequency for each line
 */
REAL4VectorSequence * trackLines(INT4Vector *lines, INT4Vector *binshifts, REAL4 minfbin, REAL4 df)
{

   XLAL_CHECK_NULL( lines != NULL && binshifts != NULL, XLAL_EINVAL );

   fprintf(stderr, "Tracking lines... ");

   REAL4VectorSequence *output = NULL;
   XLAL_CHECK_NULL( (output = XLALCreateREAL4VectorSequence(lines->length, 3)) != NULL, XLAL_EFUNC );

   //REAL4 df = 1.0/params->Tsft;
   //REAL4 minfbin = (REAL4)(round(params->fmin*params->Tsft - params->dfmax*params->Tsft - 6.0 - 0.5*(params->blksize-1) - (REAL8)(maxbinshift))/params->Tsft);

   INT4 maxshiftindex = 0, minshiftindex = 0;
   min_max_index_INT4Vector(binshifts, &minshiftindex, &maxshiftindex);
   INT4 maxshift = binshifts->data[maxshiftindex], minshift = binshifts->data[minshiftindex];

   for (UINT4 ii=0; ii<lines->length; ii++) {
      output->data[ii*3] = lines->data[ii]*df + minfbin;
      //output->data[ii*3 + 1] = (lines->data[ii] + minshift)*df + minfbin;
      //output->data[ii*3 + 2] = (lines->data[ii] + maxshift)*df + minfbin;
      output->data[ii*3 + 1] = (lines->data[ii] + (minshift-1))*df + minfbin;  //Add one extra bin for buffer
      output->data[ii*3 + 2] = (lines->data[ii] + (maxshift+1))*df + minfbin;  //Add one extra bin for buffer
   }

   fprintf(stderr, "done\n");

   return output;

}


/**
 * Test algorithm to clean lines. NOT FULLY TESTED!
 * \param [in,out] TFdata     Pointer to time-frequency data in a REAL4Vector to be cleaned
 * \param [in]     background Pointer to REAL4Vector of running mean data
 * \param [in]     lines      Pointer to INT4Vector of lines
 * \param [in]     params     Pointer to UserInput_t
 * \param [in]     rng        Pointer to gsl_rng
 * \return Status value
 */
INT4 cleanLines(REAL4Vector *TFdata, REAL4Vector *background, INT4Vector *lines, UserInput_t *params, gsl_rng *rng)
{

   if (lines==NULL) return XLAL_SUCCESS;

   INT4 numffts = (INT4)floor(params->Tobs/(params->Tsft-params->SFToverlap)-1);    //Number of FFTs
   INT4 numfbins = (INT4)TFdata->length/numffts;     //Number of frequency bins

   REAL8 prevnoiseval = 0.0;
   for (INT4 ii=0; ii<numffts; ii++) {
      if (TFdata->data[ii*numfbins]!=0.0) {
         for (UINT4 jj=0; jj<lines->length; jj++) {
            if (lines->data[jj]-(params->blksize-1)/2>=0) {
               REAL8 noiseval = expRandNum(background->data[ii*(numfbins-(params->blksize-1)) + lines->data[jj]], rng);
               XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC, "ii=%d, jj=%d, background=%g", ii, jj, background->data[ii*(numfbins-(params->blksize-1)) + lines->data[jj]] );
               if (ii>0 && TFdata->data[(ii-1)*numfbins]!=0.0) {
                  noiseval *= (1.0-0.167*0.167);
                  noiseval += 0.167*prevnoiseval;
               }
               TFdata->data[ii*numfbins + lines->data[jj]] = noiseval;
               prevnoiseval = noiseval;
            }
         } //end loop over vector of lines
      } //if SFT exists
   } //end loop over SFTs

   return XLAL_SUCCESS;

}


/**
 * Compute the second Fourier transform for TwoSpect
 * \param [out] output Pointer to the ffdataStruct to the containers for the second FFT
 * \param [in]  tfdata Pointer REAL4Vector of mean subtracted and weighted data
 * \param [in]  plan   Pointer to REAL4FFTPlan
 * \return Status value
 */
INT4 makeSecondFFT(ffdataStruct *output, REAL4Vector *tfdata, REAL4FFTPlan *plan)
{

   XLAL_CHECK( output != NULL && tfdata != NULL && plan != NULL, XLAL_EINVAL );

   fprintf(stderr, "Computing second FFT over SFTs... ");

   REAL8 winFactor = 8.0/3.0;

   //Do the second FFT
   REAL4Vector *x = NULL, *psd = NULL;
   REAL4Window *win = NULL;
   XLAL_CHECK( (x = XLALCreateREAL4Vector(output->numffts)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (psd = XLALCreateREAL4Vector((UINT4)floor(x->length*0.5)+1)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (win = XLALCreateHannREAL4Window(x->length)) != NULL, XLAL_EFUNC );

   for (INT4 ii=0; ii<output->numfbins; ii++) {

      //Next, loop over times and pick the right frequency bin for each FFT and window
      //for (UINT4 jj=0; jj<x->length; jj++) x->data[jj] = tfdata->data[ii + jj*output->numfbins] * win->data->data[jj];
      fastSSVectorMultiply_with_stride_and_offset(x, tfdata, win->data, output->numfbins, 1, ii, 0);

      //Make the FFT
      XLAL_CHECK( XLALREAL4PowerSpectrum(psd, x, plan) == XLAL_SUCCESS, XLAL_EFUNC );

      //Fix beginning and end values if even, otherwise just the beginning if odd
      if (GSL_IS_EVEN(x->length)==1) {
         psd->data[0] *= 2.0;
         psd->data[psd->length-1] *= 2.0;
      } else {
         psd->data[0] *= 2.0;
      }

      //Scale the data points by 1/N and window factor and (1/fs)
      //Order of vector is by second frequency then first frequency
      //It is possible that when dealing with very loud signals, lines, injections, etc. (e.g., far above the background)
      //then the output power here can be "rounded" because of the cast to nearby integer values.
      //For high (but not too high) power values, this may not be noticed because the cast can round to nearby decimal values.
      for (UINT4 jj=0; jj<psd->length; jj++) output->ffdata->data[psd->length*ii + jj] = (REAL4)(psd->data[jj]*winFactor*output->ffnormalization);

   } /* for ii < numfbins */

   //Destroy stuff
   XLALDestroyREAL4Vector(x);
   XLALDestroyREAL4Vector(psd);
   XLALDestroyREAL4Window(win);

   fprintf(stderr, "done\n");

   return XLAL_SUCCESS;

} /* makeSecondFFT() */


/**
 * Determine the average of the noise power in each frequency bin across the band
 * \param [in] backgrnd Pointer to REAL4Vector of the running mean values
 * \param [in] numfbins Number of frequency bins in the SFTs
 * \param [in] numffts  Number of SFTs in the observation time
 * \param [in] binmin   Minimum SFT bin to look at with this algorithm
 * \param [in] binmax   Maximum SFT bin to look at with this algorithm
 * \return The average of the noise power across the frequency band
 */
REAL4 avgTFdataBand(REAL4Vector *backgrnd, INT4 numfbins, INT4 numffts, INT4 binmin, INT4 binmax)
{

   XLAL_CHECK_REAL4( backgrnd != NULL && numfbins > 0 && numffts > 0 && binmin >= 0 && binmax > binmin, XLAL_EINVAL );

   REAL4Vector *aveNoiseInTime = NULL, *rngMeansOverBand = NULL;
   XLAL_CHECK_REAL4( (aveNoiseInTime = XLALCreateREAL4Vector(numffts)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_REAL4( (rngMeansOverBand = XLALCreateREAL4Vector((UINT4)(binmax-binmin))) != NULL, XLAL_EFUNC );

   for (UINT4 ii=0; ii<aveNoiseInTime->length; ii++) {
      memcpy(rngMeansOverBand->data, &(backgrnd->data[ii*numfbins + binmin]), sizeof(*rngMeansOverBand->data)*rngMeansOverBand->length);
      aveNoiseInTime->data[ii] = calcMean(rngMeansOverBand);
   } /* for ii < aveNoiseInTime->length */

   REAL4 avgTFdata = calcMean(aveNoiseInTime);

   //Destroy stuff
   XLALDestroyREAL4Vector(aveNoiseInTime);
   XLALDestroyREAL4Vector(rngMeansOverBand);

   return avgTFdata;

} /* avgTFdataBand() */


/**
 * Determine the rms of the noise power in each frequency bin across the band
 * \param [in] backgrnd Pointer to REAL4Vector of the running mean values
 * \param [in] numfbins Number of frequency bins in the SFTs
 * \param [in] numffts  Number of SFTs in the observation time
 * \param [in] binmin   Minimum SFT bin to look at with this algorithm
 * \param [in] binmax   Maximum SFT bin to look at with this algorithm
 * \return The RMS of the noise power across the frequency band
 */
REAL4 rmsTFdataBand(REAL4Vector *backgrnd, INT4 numfbins, INT4 numffts, INT4 binmin, INT4 binmax)
{

   XLAL_CHECK_REAL4( backgrnd != NULL && numfbins > 0 && numffts > 0 && binmin >= 0 && binmax > binmin, XLAL_EINVAL );

   REAL4Vector *aveNoiseInTime = NULL, *rngMeansOverBand = NULL;
   XLAL_CHECK_REAL4( (aveNoiseInTime = XLALCreateREAL4Vector(numffts)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_REAL4( (rngMeansOverBand = XLALCreateREAL4Vector((UINT4)(binmax-binmin))) != NULL, XLAL_EFUNC );

   for (UINT4 ii=0; ii<aveNoiseInTime->length; ii++) {
      memcpy(rngMeansOverBand->data, &(backgrnd->data[ii*numfbins + binmin]), sizeof(*rngMeansOverBand->data)*rngMeansOverBand->length);
      aveNoiseInTime->data[ii] = calcMean(rngMeansOverBand);
   } /* for ii < aveNoiseInTime->length */

   REAL4 rmsTFdata = 0.0;
   XLAL_CHECK_REAL4( calcRms(&rmsTFdata, aveNoiseInTime) == XLAL_SUCCESS, XLAL_EFUNC );

   //Destroy stuff
   XLALDestroyREAL4Vector(aveNoiseInTime);
   XLALDestroyREAL4Vector(rngMeansOverBand);

   return rmsTFdata;

} /* rmsTFdataBand() */


/**
 * Measure of the average noise power in each 2nd FFT frequency bin
 * \param [out]    aveNoise      Pointer to REAL4VectorAligned of the expected 2nd FFT powers
 * \param [in]     params        Pointer to UserInput_t
 * \param [in]     sftexist      Pointer to INT4Vector of SFTs existing or not
 * \param [in]     backgrnd      Pointer to REAL4Vector of running means
 * \param [in]     antweights    Pointer to REAL4Vector of antenna pattern weights
 * \param [in]     plan          Pointer to REAL4FFTPlan
 * \param [in,out] normalization Pointer to REAL8 value of the normalization for the 2nd FFT
 * \param [in]     rng           Pointer to gsl_rng
 * \return Status value
 */
INT4 ffPlaneNoise(REAL4VectorAligned *aveNoise, UserInput_t *params, INT4Vector *sftexist, REAL4Vector *backgrnd, REAL4Vector *antweights, REAL4FFTPlan *plan, REAL8 *normalization, gsl_rng *rng)
{

   XLAL_CHECK( aveNoise != NULL && params != NULL && sftexist != NULL && backgrnd != NULL && antweights != NULL && plan != NULL && normalization != NULL, XLAL_EINVAL );

   fprintf(stderr, "Computing noise background estimate... ");

   REAL8 invsumofweights = 0.0, sumofweights = 0.0;

   UINT4 numffts = antweights->length;
   UINT4 numfbins = backgrnd->length/numffts;
   UINT4 numfprbins = aveNoise->length;

   //Set up for making the PSD
   memset(aveNoise->data, 0, sizeof(REAL4)*aveNoise->length);

   //If the user has said there is signal only and no noise in the SFTs, then the noise background of the FF plane will be filled with zeros
   if (!params->signalOnly) {
      //Window and psd allocation
      REAL4Window *win = NULL;
      REAL4VectorAligned *psd = NULL;
      XLAL_CHECK( (win = XLALCreateHannREAL4Window(numffts)) != NULL, XLAL_EFUNC );  //Window function
      XLAL_CHECK( (psd = XLALCreateREAL4VectorAligned(numfprbins, 32)) != NULL, XLAL_EFUNC );   //Current PSD calculation
      REAL4 winFactor = 8.0/3.0;
      REAL8 dutyfactor = 0.0, dutyfactorincrement = 1.0/(REAL8)numffts;

      //Average each SFT across the frequency band, also compute normalization factor
      REAL4Vector *rngMeansOverBand = NULL;
      REAL4VectorAligned *aveNoiseInTime = NULL;
      XLAL_CHECK( (aveNoiseInTime = XLALCreateREAL4VectorAligned(numffts, 32)) != NULL, XLAL_EFUNC );
      XLAL_CHECK( (rngMeansOverBand = XLALCreateREAL4Vector(numfbins)) != NULL, XLAL_EFUNC );

      memset(aveNoiseInTime->data, 0, sizeof(REAL4)*numffts);
      for (UINT4 ii=0; ii<aveNoiseInTime->length; ii++) {
         if (sftexist->data[ii]!=0) {
            memcpy(rngMeansOverBand->data, &(backgrnd->data[ii*numfbins]), sizeof(*rngMeansOverBand->data)*rngMeansOverBand->length);
            XLAL_CHECK( calcMedian(&(aveNoiseInTime->data[ii]), rngMeansOverBand) == XLAL_SUCCESS, XLAL_EFUNC );

            if (!params->noiseWeightOff) sumofweights += (antweights->data[ii]*antweights->data[ii])/(aveNoiseInTime->data[ii]*aveNoiseInTime->data[ii]);
            else sumofweights += (antweights->data[ii]*antweights->data[ii]);

            dutyfactor += dutyfactorincrement;
         }
      } /* for ii < aveNoiseInTime->length */
      invsumofweights = 1.0/sumofweights;

      //Load time series of powers, normalize, mean subtract and Hann window
      REAL4VectorAligned *x = NULL, *multiplicativeFactor = NULL;
      XLAL_CHECK( (x = XLALCreateREAL4VectorAligned(aveNoiseInTime->length, 32)) != NULL, XLAL_EFUNC );
      XLAL_CHECK( (multiplicativeFactor = XLALCreateREAL4VectorAligned(aveNoiseInTime->length, 32)) != NULL, XLAL_EFUNC );

      memset(multiplicativeFactor->data, 0, sizeof(REAL4)*multiplicativeFactor->length);
      REAL4 psdfactor = winFactor*params->Tobs;  //Only multiply by Tobs instead of Tobs/numffts^2 because of the weighting normalization

      for (UINT4 ii=0; ii<x->length; ii++) {
         if (aveNoiseInTime->data[ii] != 0.0) {
            if (!params->noiseWeightOff) multiplicativeFactor->data[ii] = win->data->data[ii]*antweights->data[ii]/aveNoiseInTime->data[ii]*invsumofweights;
            else multiplicativeFactor->data[ii] = win->data->data[ii]*antweights->data[ii]*invsumofweights;
         }
      }

      //Create exponentially distributed values from which to resample for measuring the background distribution
      //Every simulation of background noise samples a subset of these 100*numffts values
      REAL4VectorAligned *randVals = NULL;
      XLAL_CHECK( (randVals = XLALCreateREAL4VectorAligned(100*numffts, 32)) != NULL, XLAL_EFUNC );
      for (UINT4 ii=0; ii<randVals->length; ii++) randVals->data[ii] = expRandNum(1.0, rng);

      //New version. Computes expected background with correlation estimate from Hann windowed and overlapped (must be 50% or 0%) SFTs
      REAL8 correlationfactor = 0.0;
      for (INT4 ii=0; ii<(INT4)floor(win->data->length*(params->SFToverlap/params->Tsft)-1); ii++) correlationfactor += win->data->data[ii]*win->data->data[ii + (INT4)((1.0-(params->SFToverlap/params->Tsft))*win->data->length)];
      correlationfactor /= win->sumofsquares;
      REAL8 corrfactorsquared = correlationfactor*correlationfactor;
      //REAL8 prevnoiseval = 0.0, noiseval = 0.0;
      for (UINT4 ii=0; ii<1000; ii++) {
         memset(x->data, 0, sizeof(REAL4)*x->length);
         /* for (UINT4 jj=0; jj<x->length; jj++) {
            if (sftexist->data[jj] != 0) {
               //To create the correlations
               noiseval = expRandNum(aveNoiseInTime->data[jj], rng);
               XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
               if (jj>0 && sftexist->data[jj-1]!=0) {
                  noiseval *= (1.0-corrfactorsquared);
                  noiseval += corrfactorsquared*prevnoiseval;
               }
               x->data[jj] = (REAL4)(noiseval/aveNoiseInTime->data[jj]-1.0);
               prevnoiseval = noiseval;
            }
         } // for jj < x->length */

         XLAL_CHECK( sampleREAL4Vector((REAL4Vector*)x, (REAL4Vector*)randVals, rng) == XLAL_SUCCESS, XLAL_EFUNC );
         if (params->vectorMath==1) XLAL_CHECK( sseSSVectorMultiply(x, x, aveNoiseInTime) == XLAL_SUCCESS, XLAL_EFUNC );
         else if (params->vectorMath==2) XLAL_CHECK( avxSSVectorMultiply(x, x, aveNoiseInTime) == XLAL_SUCCESS, XLAL_EFUNC );
         else for (UINT4 jj=0; jj<x->length; jj++) x->data[jj] *= aveNoiseInTime->data[jj];
         for (UINT4 jj=0; jj<x->length; jj++) {
            if (sftexist->data[jj] != 0) {
               if (jj>0 && sftexist->data[jj-1]!=0) {
                  x->data[jj] *= (1.0 - corrfactorsquared);
                  x->data[jj] += corrfactorsquared*x->data[jj-1];
               }
            }
         }
         for (UINT4 jj=0; jj<x->length; jj++) if (sftexist->data[jj] != 0) x->data[jj] = x->data[jj]/aveNoiseInTime->data[jj] - 1.0;

         //Reverse the correlations, comment this out below!
         /* memset(x->data, 0, sizeof(REAL4)*x->length);
         for (jj=(INT4)x->length-1; jj>=0; jj--) {
            if (sftexist->data[jj] != 0) {
               noiseval = expRandNum(aveNoiseInTime->data[jj], rng);
               if (XLAL_IS_REAL8_FAIL_NAN(noiseval)) {
                  fprintf(stderr, "%s: expRandNum() failed.\n", __func__);
                  XLAL_ERROR_VOID(XLAL_EFUNC);
               }
               if (jj<(INT4)x->length-1 && sftexist->data[jj+1]!=0) {
                  noiseval *= (1.0-corrfactorsquared);
                  noiseval += corrfactorsquared*prevnoiseval;
               }
               x->data[jj] = (REAL4)(noiseval/aveNoiseInTime->data[jj]-1.0);
               prevnoiseval = noiseval;
            }
            } */

         //Window and rescale because of antenna and noise weights
         if (params->vectorMath==1) XLAL_CHECK( sseSSVectorMultiply(x, x, multiplicativeFactor) == XLAL_SUCCESS, XLAL_EFUNC );
         else if (params->vectorMath==2) XLAL_CHECK( avxSSVectorMultiply(x, x, multiplicativeFactor) == XLAL_SUCCESS, XLAL_EFUNC );
         else for (UINT4 jj=0; jj<x->length; jj++) x->data[jj] *= multiplicativeFactor->data[jj];

         //Do the FFT
         XLAL_CHECK( XLALREAL4PowerSpectrum((REAL4Vector*)psd, (REAL4Vector*)x, plan) == XLAL_SUCCESS, XLAL_EFUNC );

         //Sum into the bins
         if (params->vectorMath==1) XLAL_CHECK( sseSSVectorSum(aveNoise, aveNoise, psd) == XLAL_SUCCESS, XLAL_EFUNC );
         else if (params->vectorMath==2) XLAL_CHECK( avxSSVectorSum(aveNoise, aveNoise, psd) == XLAL_SUCCESS, XLAL_EFUNC );
         else for (UINT4 jj=0; jj<aveNoise->length; jj++) aveNoise->data[jj] += psd->data[jj];
      } /* for ii < 1000 */

      //Average and rescale
      //REAL4 averageRescaleFactor = 2.5e-4*psdfactor*(1.0+2.0*corrfactorsquared);
      REAL4 averageRescaleFactor = 1.0e-3*psdfactor; //*(1.0+2.0*corrfactorsquared);
      if (params->vectorMath==1) XLAL_CHECK( sseScaleREAL4Vector(aveNoise, aveNoise, averageRescaleFactor) == XLAL_SUCCESS, XLAL_EFUNC );
      else if (params->vectorMath==2) XLAL_CHECK( avxScaleREAL4Vector(aveNoise, aveNoise, averageRescaleFactor) == XLAL_SUCCESS, XLAL_EFUNC );
      else for (UINT4 ii=0; ii<aveNoise->length; ii++) aveNoise->data[ii] *= averageRescaleFactor;

      //Fix 0th and end bins (0 only for odd x->length, 0 and end for even x->length)
      if (GSL_IS_EVEN(x->length)==1) {
         aveNoise->data[0] *= 2.0;
         aveNoise->data[aveNoise->length-1] *= 2.0;
      } else {
         aveNoise->data[0] *= 2.0;
      }

      //Compute normalization
      *(normalization) = 1.0/(calcMean((REAL4Vector*)aveNoise));
      for (UINT4 ii=0; ii<aveNoise->length; ii++) aveNoise->data[ii] *= *(normalization);

      //Extra factor for normalization to 1.0 (empirically determined) (Old method)
      //*normalization /= 1.08;

      //Extra factor for normalization to 1.0
      //duty factor because there are zero-valued SFTs (fourth root term)
      //square root term because the uncertainty on the background improves with the size of the running median block
      //correlation term because neighboring sfts are correlated
      *normalization /= (1.0+pow(dutyfactor,0.25)/sqrt((REAL8)params->blksize))*sqrt(1.0+2.0*corrfactorsquared);

      //Destroy stuff
      XLALDestroyREAL4VectorAligned(x);
      XLALDestroyREAL4VectorAligned(psd);
      XLALDestroyREAL4Window(win);
      XLALDestroyREAL4VectorAligned(aveNoiseInTime);
      XLALDestroyREAL4Vector(rngMeansOverBand);
      XLALDestroyREAL4VectorAligned(multiplicativeFactor);
      XLALDestroyREAL4VectorAligned(randVals);
   } else {
     *(normalization) = 1.0;
   }

   fprintf(stderr, "done\n");

   return XLAL_SUCCESS;

} /* ffPlaneNoise() */


INT4 readTwoSpectInputParams(UserInput_t *uvar, int argc, char *argv[])
{

   XLAL_CHECK ( uvar != NULL, XLAL_EINVAL, "Invalid NULL input 'uvar'\n");
   XLAL_CHECK ( argv != NULL, XLAL_EINVAL, "Invalid NULL input 'argv'\n");

   uvar->ephemEarth = XLALStringDuplicate("earth00-19-DE405.dat.gz");
   uvar->ephemSun = XLALStringDuplicate("sun00-19-DE405.dat.gz");
   uvar->blksize = 101;
   uvar->outfilename = XLALStringDuplicate("logfile.txt");
   uvar->configCopy = XLALStringDuplicate("input_args.conf");
   uvar->ULfilename = XLALStringDuplicate("uls.dat");
   uvar->harmonicNumToSearch = 1;
   uvar->periodHarmToCheck = 5;
   uvar->periodFracToCheck = 3;
   uvar->ihsfactor = 5;
   uvar->minTemplateLength = 1;
   uvar->maxTemplateLength = 500;
   uvar->FFTplanFlag = 1;
   uvar->vectorMath = 0;
   uvar->injRandSeed = 0;
   uvar->ULsolver = 0;
   uvar->dopplerMultiplier = 1.0;
   uvar->ihsfar = 1.0;
   uvar->tmplfar = 1.0;
   uvar->ihsfomfar = 1.0;
   uvar->ihsfom = 0.0;
   uvar->keepOnlyTopNumIHS = -1;
   uvar->lineDetection = -1.0;
   uvar->cosiSignCoherent = 0;

   XLALregBOOLUserStruct(help,                         'h', UVAR_HELP    ,  "Print this help/usage message");
   XLALregSTRINGUserStruct(outdirectory,                0 , UVAR_REQUIRED,  "Output directory");
   XLALregLISTUserStruct(IFO,                           0 , UVAR_REQUIRED,  "CSV list of detectors, eg. \"H1,H2,L1,G1, ...\" ");
   XLALregREALUserStruct(Tobs,                          0 , UVAR_REQUIRED,  "Total observation time (in seconds)");
   XLALregREALUserStruct(Tsft,                          0 , UVAR_REQUIRED,  "SFT coherence time (in seconds)");
   XLALregREALUserStruct(SFToverlap,                    0 , UVAR_REQUIRED,  "SFT overlap (in seconds), usually Tsft/2");
   XLALregREALUserStruct(t0,                            0 , UVAR_REQUIRED,  "Start time of the search (in GPS seconds)");
   XLALregREALUserStruct(fmin,                          0 , UVAR_REQUIRED,  "Minimum frequency of band (Hz)");
   XLALregREALUserStruct(fspan,                         0 , UVAR_REQUIRED,  "Frequency span of band (Hz)");
   XLALregREALUserStruct(avesqrtSh,                     0 , UVAR_REQUIRED,  "Expected average of square root of Sh");
   XLALregREALUserStruct(Pmin,                          0 , UVAR_REQUIRED,  "Minimum period to be searched (in seconds)");
   XLALregREALUserStruct(Pmax,                          0 , UVAR_REQUIRED,  "Maximum period to be searched (in seconds)");
   XLALregREALUserStruct(dfmin,                         0 , UVAR_REQUIRED,  "Minimum modulation depth to search (Hz)");
   XLALregREALUserStruct(dfmax,                         0 , UVAR_REQUIRED,  "Maximum modulation depth to search (Hz)");
   XLALregINTUserStruct(blksize,                        0 , UVAR_OPTIONAL,  "Blocksize for running median to determine expected noise of input SFTs");
   XLALregSTRINGUserStruct(outfilename,                 0 , UVAR_OPTIONAL,  "Output file name");
   XLALregSTRINGUserStruct(configCopy,                  0 , UVAR_OPTIONAL,  "Copy of the input values");
   XLALregSTRINGUserStruct(ULfilename,                  0 , UVAR_OPTIONAL,  "Upper limit file name");
   XLALregSTRINGUserStruct(inputSFTs,                   0 , UVAR_OPTIONAL,  "Path and filename of SFTs, conflicts with timestampsFile and segmentFile");
   XLALregSTRINGUserStruct(ephemEarth,                  0 , UVAR_OPTIONAL,  "Earth ephemeris file to use");
   XLALregSTRINGUserStruct(ephemSun,                    0 , UVAR_OPTIONAL,  "Sun ephemeris file to use");
   XLALregSTRINGUserStruct(skyRegion,                   0 , UVAR_OPTIONAL,  "Region of the sky to search (e.g. (ra1,dec1),(ra2,dec2),(ra3,dec3)...) or allsky");
   XLALregSTRINGUserStruct(skyRegionFile,               0 , UVAR_OPTIONAL,  "File with the grid points");
   XLALregREALUserStruct(linPolAngle,                   0 , UVAR_OPTIONAL,  "Polarization angle to search using linear polarization of Fplus (when unspecified default is circular polarization)");
   XLALregSTRINGUserStruct(templatebankfile,            0 , UVAR_OPTIONAL,  "File containing the template data (generated by lalapps_TwoSpectTemplateBank)");
   XLALregINTUserStruct(harmonicNumToSearch,            0 , UVAR_OPTIONAL,  "Number of harmonics of the Pmin to Pmax range to search");
   XLALregINTUserStruct(periodHarmToCheck,              0 , UVAR_OPTIONAL,  "Number of harmonics/sub-harmonics of the IHS candidates to test");
   XLALregINTUserStruct(periodFracToCheck,              0 , UVAR_OPTIONAL,  "Number of fractional periods to check in the sense of [(1...N)+1]/[(1...N)+2]");
   XLALregBOOLUserStruct(templateSearch,                0 , UVAR_OPTIONAL,  "Flag for doing a pure template-based search on search region specified by (sky,f,fspan,P, Asini +- 3 AsiniSigma)");
   XLALregREALUserStruct(templateSearchP,               0 , UVAR_OPTIONAL,  "The template search period; templateSearch flag is required");
   XLALregREALUserStruct(templateSearchAsini,           0 , UVAR_OPTIONAL,  "The template search Asini; templateSearch flag is required");
   XLALregREALUserStruct(templateSearchAsiniSigma,      0 , UVAR_OPTIONAL,  "The template search uncertainty in Asini; templateSearch flag is required");
   XLALregREALUserStruct(assumeNScosi,                  0 , UVAR_OPTIONAL,  "Assume cosi orientation of the source (only used when specifying more than 1 detector)");
   XLALregREALUserStruct(assumeNSpsi,                   0 , UVAR_OPTIONAL,  "Assume psi polarization angle of the source (only used when specifying more than 1 detector)");
   XLALregINTUserStruct(cosiSignCoherent,               0 , UVAR_OPTIONAL,  "For coherent analysis assume [-1,1] values (0), [0,1] values (1), or [-1,0] values (-1) for cosi (Note: unused when assumeNScosi is specified)");
   XLALregINTUserStruct(ihsfactor,                      0 , UVAR_OPTIONAL,  "Number of harmonics to sum in IHS algorithm");
   XLALregREALUserStruct(ihsfar,                        0 , UVAR_OPTIONAL,  "IHS FAR threshold");
   XLALregREALUserStruct(ihsfom,                        0 , UVAR_OPTIONAL,  "IHS FOM = 12*(L_IHS_loc - U_IHS_loc)^2");
   XLALregREALUserStruct(ihsfomfar,                     0 , UVAR_OPTIONAL,  "IHS FOM FAR threshold");
   XLALregREALUserStruct(tmplfar,                       0 , UVAR_OPTIONAL,  "Template FAR threshold");
   XLALregINTUserStruct(keepOnlyTopNumIHS,              0 , UVAR_OPTIONAL,  "Keep the top <number> of IHS candidates based on significance");
   XLALregINTUserStruct(minTemplateLength,              0 , UVAR_OPTIONAL,  "Minimum number of pixels to use in the template");
   XLALregINTUserStruct(maxTemplateLength,              0 , UVAR_OPTIONAL,  "Maximum number of pixels to use in the template");
   XLALregREALUserStruct(ULfmin,                        0 , UVAR_OPTIONAL,  "Minimum signal frequency considered for the upper limit value (Hz)");
   XLALregREALUserStruct(ULfspan,                       0 , UVAR_OPTIONAL,  "Span of signal frequencies considered for the upper limit value (Hz)");
   XLALregREALUserStruct(ULminimumDeltaf,               0 , UVAR_OPTIONAL,  "Minimum modulation depth counted in the upper limit value (Hz)");
   XLALregREALUserStruct(ULmaximumDeltaf,               0 , UVAR_OPTIONAL,  "Maximum modulation depth counted in the upper limit value (Hz)");
   XLALregBOOLUserStruct(allULvalsPerSkyLoc,            0 , UVAR_OPTIONAL,  "Print all UL values in the band specified by ULminimumDeltaf and ULmaximumDeltaf (default prints only the maximum UL)");
   XLALregBOOLUserStruct(markBadSFTs,                   0 , UVAR_OPTIONAL,  "Mark bad SFTs");
   XLALregREALUserStruct(lineDetection,                 0 , UVAR_OPTIONAL,  "Detect stationary lines above threshold, and, if any present, set upper limit only, no template follow-up");
   XLALregINTUserStruct(FFTplanFlag,                    0 , UVAR_OPTIONAL,  "0=Estimate, 1=Measure, 2=Patient, 3=Exhaustive");
   XLALregBOOLUserStruct(fastchisqinv,                  0 , UVAR_OPTIONAL,  "Use a faster central chi-sq inversion function (roughly float precision instead of double)");
   XLALregINTUserStruct(vectorMath,                     0 , UVAR_OPTIONAL,  "Vector math functions: 0=None, 1=SSE, 2=AVX/SSE (Note that user needs to have compiled for SSE or AVX/SSE or program fails)");
   XLALregBOOLUserStruct(followUpOutsideULrange,        0 , UVAR_OPTIONAL,  "Follow up outliers outside the range of the UL values");
   XLALregLISTUserStruct(timestampsFile,                0 , UVAR_OPTIONAL,  "CSV list of files with timestamps, file-format: lines of <GPSsec> <GPSnsec>, conflicts with inputSFTs and segmentFile");
   XLALregLISTUserStruct(segmentFile,                   0 , UVAR_OPTIONAL,  "CSV list of files with segments, file-format: lines with <GPSstart> <GPSend>, conflicts with inputSFTs and timestampsFile");
   XLALregBOOLUserStruct(gaussNoiseWithSFTgaps,         0 , UVAR_OPTIONAL,  "Gaussian noise using avesqrtSh with same gaps as inputSFTs, conflicts with timstampsFile and segmentFile");
   XLALregLISTUserStruct(injectionSources,              0 , UVAR_OPTIONAL,  "CSV file list containing sources to inject or '{Alpha=0;Delta=0;...}'");
   XLALregREALUserStruct(injFmin,                       0 , UVAR_OPTIONAL,  "Minimum frequency of band to create in TwoSpect");
   XLALregREALUserStruct(injBand,                       0 , UVAR_OPTIONAL,  "Width of band to create in TwoSpect");
   XLALregINTUserStruct(injRandSeed,                    0 , UVAR_OPTIONAL,  "Random seed value for reproducable noise, conflicts with inputSFTs");
   XLALregBOOLUserStruct(weightedIHS,                   0 , UVAR_DEVELOPER, "Use the noise-weighted IHS scheme");
   XLALregBOOLUserStruct(signalOnly,                    0 , UVAR_DEVELOPER, "SFTs contain only signal, no noise");
   XLALregBOOLUserStruct(templateTest,                  0 , UVAR_DEVELOPER, "Test the doubly-Fourier transformed data against a single, exact template");
   XLALregREALUserStruct(templateTestF,                 0 , UVAR_DEVELOPER, "The template test frequency; templateTest flag is required");
   XLALregREALUserStruct(templateTestP,                 0 , UVAR_DEVELOPER, "The template test period; templateTest flag is required");
   XLALregREALUserStruct(templateTestDf,                0 , UVAR_DEVELOPER, "The template test modulation depth; templateTest flag is required");
   XLALregBOOLUserStruct(bruteForceTemplateTest,        0 , UVAR_DEVELOPER, "Test a number of different templates using templateTest parameters");
   XLALregINTUserStruct(ULsolver,                       0 , UVAR_DEVELOPER, "Solver function for the upper limit calculation: 0=gsl_ncx2cdf_float_withouttinyprob_solver, 1=gsl_ncx2cdf_withouttinyprob_solver, 2=gsl_ncx2cdf_float_solver, 3=gsl_ncx2cdf_solver, 4=ncx2cdf_float_withouttinyprob_withmatlabchi2cdf_solver, 5=ncx2cdf_withouttinyprob_withmatlabchi2cdf_solver");
   XLALregREALUserStruct(dopplerMultiplier,             0 , UVAR_DEVELOPER, "Multiplier for the Doppler velocity");
   XLALregBOOLUserStruct(IHSonly,                       0 , UVAR_DEVELOPER, "IHS stage only is run. Output statistic is the IHS statistic");
   XLALregBOOLUserStruct(noNotchHarmonics,              0 , UVAR_DEVELOPER, "Do not notch the daily/sidereal harmonics in the IHS step");
   XLALregBOOLUserStruct(calcRthreshold,                0 , UVAR_DEVELOPER, "Calculate the threshold value for R given the template false alarm rate");
   XLALregBOOLUserStruct(antennaOff,                    0 , UVAR_DEVELOPER, "Antenna pattern weights are /NOT/ used if this flag is used");
   XLALregBOOLUserStruct(noiseWeightOff,                0 , UVAR_DEVELOPER, "Turn off noise weighting if this flag is used");
   XLALregBOOLUserStruct(gaussTemplatesOnly,            0 , UVAR_DEVELOPER, "Gaussian templates only throughout the pipeline if this flag is used");
   XLALregBOOLUserStruct(ULoff,                         0 , UVAR_DEVELOPER, "Turn off upper limits computation");
   XLALregBOOLUserStruct(BrentsMethod,                  0 , UVAR_DEVELOPER, "Use Brent's method for root finding of the R threshold value");
   XLALregBOOLUserStruct(printSFTtimes,                 0 , UVAR_DEVELOPER, "Output a list <GPS sec> <GPS nanosec> of SFT start times of SFTs");
   XLALregBOOLUserStruct(printData,                     0 , UVAR_DEVELOPER, "Print to ASCII files the data values");
   XLALregSTRINGUserStruct(saveRvalues,                 0 , UVAR_DEVELOPER, "Print all R values from template bank file");
   XLALregSTRINGUserStruct(printSignalData,             0 , UVAR_DEVELOPER, "Print f0 and h0 per SFT of the signal, used only with injectionSources");
   XLALregSTRINGUserStruct(printMarginalizedSignalData, 0 , UVAR_DEVELOPER, "Print f0 and h0 per SFT of the signal, used only with injectionSources");
   XLALregINTUserStruct(randSeed,                       0 , UVAR_DEVELOPER, "Random seed value");
   XLALregBOOLUserStruct(chooseSeed,                    0 , UVAR_DEVELOPER, "The random seed value is chosen based on the input search parameters");

   //Read all the input from config file and command line (command line has priority)
   //Also checks required variables unless help is requested
   XLAL_CHECK( XLALUserVarReadAllInput(argc, argv) == XLAL_SUCCESS, XLAL_EFUNC );

   //Help and exit
   if (uvar->help) exit(0);

   //Check analysis parameters
   if (ceil(uvar->t0/uvar->Tsft)*uvar->Tsft - uvar->t0 != 0.0) {
      REAL8 oldstart = uvar->t0;
      uvar->t0 = ceil(uvar->t0/uvar->Tsft)*uvar->Tsft;
      fprintf(stderr, "WARNING! Adjusting start time from %f to %f\n", oldstart, uvar->t0);
   }
   XLAL_CHECK( uvar->Pmax >= uvar->Pmin, XLAL_EINVAL, "Pmax is smaller than Pmin\n" );
   XLAL_CHECK( uvar->dfmax >= uvar->dfmin, XLAL_EINVAL, "dfmax is smaller than dfmin\n" );
   if (uvar->Pmax > 0.2*uvar->Tobs) {
      uvar->Pmax = 0.2*uvar->Tobs;
      fprintf(stderr,"WARNING! Adjusting input maximum period to 1/5 the observation time!\n");
   }
   if (uvar->Pmin < 4.0*uvar->Tsft) {
      uvar->Pmin = 4.0*uvar->Tsft;
      fprintf(stderr,"WARNING! Adjusting input minimum period to 4 coherence times!\n");
   }
   if (uvar->dfmax > maxModDepth(uvar->Pmax, uvar->Tsft)) {
      uvar->dfmax = floor(maxModDepth(uvar->Pmax, uvar->Tsft)*(uvar->Tsft))/(uvar->Tsft);
      fprintf(stderr,"WARNING! Adjusting input maximum modulation depth due to maximum modulation depth allowed\n");
   }
   if (uvar->dfmin < 0.5/uvar->Tsft) {
      uvar->dfmin = 0.5/uvar->Tsft;
      fprintf(stderr,"WARNING! Adjusting input minimum modulation depth to 1/2 a frequency bin!\n");
   }
   uvar->dfmin = 0.5*round(2.0*uvar->dfmin*uvar->Tsft)/uvar->Tsft;
   uvar->dfmax = 0.5*round(2.0*uvar->dfmax*uvar->Tsft)/uvar->Tsft;

   //Check blocksize settings for running mean
   if (uvar->blksize % 2 != 1) uvar->blksize += 1;

   //Check timestampsFile/inputSFTs/segmentFile conflicts
   if (XLALUserVarWasSet(&uvar->timestampsFile) && XLALUserVarWasSet(&uvar->inputSFTs)) {
      XLAL_ERROR(XLAL_FAILURE, "Cannot specify timestampsFile and inputSFTs\n");
   } else if (XLALUserVarWasSet(&uvar->segmentFile) && XLALUserVarWasSet(&uvar->inputSFTs)) {
      XLAL_ERROR(XLAL_FAILURE, "Cannot specify segmentFile and inputSFTs\n");
   } else if (XLALUserVarWasSet(&uvar->segmentFile) && XLALUserVarWasSet(&uvar->timestampsFile)) {
      XLAL_ERROR(XLAL_FAILURE, "Cannot specify segmentFile and timestampsFile\n");
   }

   //Check skyRegion and skyRegionFile options
   if ((XLALUserVarWasSet(&uvar->skyRegion) && XLALUserVarWasSet(&uvar->skyRegionFile)) || (!XLALUserVarWasSet(&uvar->skyRegion) && !XLALUserVarWasSet(&uvar->skyRegionFile))) XLAL_ERROR(XLAL_EINVAL, "Specify only one of skyRegion or skyRegionFile\n");

   //Check required options if specifying templateSearch
   if (uvar->templateSearch || XLALUserVarWasSet(&uvar->templateSearchP) || XLALUserVarWasSet(&uvar->templateSearchAsini) || XLALUserVarWasSet(&uvar->templateSearchAsiniSigma)) {
      if (!(uvar->templateSearch && XLALUserVarWasSet(&uvar->templateSearchP) && XLALUserVarWasSet(&uvar->templateSearchAsini) && XLALUserVarWasSet(&uvar->templateSearchAsiniSigma))) XLAL_ERROR(XLAL_FAILURE, "Must specify all or none of templateSearch, templateSearchP, templateSearchAsini, and templateSearchAsiniSigma\n");
   }

   //Error if dfmax is smaller than templateTestDf or if dfmax is smaller than the templateSearch or bruteForceTemplateTest largest modulation depth
   if (uvar->templateTest) XLAL_CHECK( uvar->dfmax >= uvar->templateTestDf, XLAL_EINVAL, "templateTestDf is larger than dfmax\n" );
   if (uvar->templateSearch) XLAL_CHECK( uvar->dfmax >= LAL_TWOPI*(uvar->fmin+uvar->fspan)*(uvar->templateSearchAsini+3.0*uvar->templateSearchAsiniSigma)/uvar->templateSearchP, XLAL_EINVAL, "templateSearch parameters would make the largest modulation depth larger than dfmax\n");
   if (uvar->bruteForceTemplateTest) XLAL_CHECK( uvar->dfmax >= uvar->templateTestDf+2.0/uvar->Tsft, XLAL_EINVAL, "templateTestDf+2/Tsft is larger than dfmax\n" );

   //Check IHS FOM threshold conflicts
   if (XLALUserVarWasSet(&uvar->ihsfomfar) && XLALUserVarWasSet(&uvar->ihsfom)) XLAL_ERROR(XLAL_FAILURE, "Only choose only one of ihsfomfar or ihsfom\n");

   //Check ranges of settings for FAR values
   if (uvar->ihsfar<0.0 || uvar->ihsfar>1.0) XLAL_ERROR(XLAL_EINVAL, "ihsfar must lie between 0 and 1 (inclusive)\n");
   if (uvar->ihsfomfar<0.0 || uvar->ihsfomfar>1.0) XLAL_ERROR(XLAL_EINVAL, "ihsfomfar must lie between 0 and 1 (inclusive)\n");
   if (uvar->tmplfar<0.0 || uvar->tmplfar>1.0) XLAL_ERROR(XLAL_EINVAL, "tmplfar must lie between 0 and 1 (inclusive)\n");

   //params->log10templatefar = log10(params->templatefar);

   //Upper limit settings take span of search values unless specified
   if (!XLALUserVarWasSet(&uvar->ULfmin)) uvar->ULfmin = uvar->fmin;            //Upper limit minimum frequency (Hz)
   if (!XLALUserVarWasSet(&uvar->ULfspan)) uvar->ULfspan = uvar->fspan;         //Upper limit frequency span (Hz)
   if (!XLALUserVarWasSet(&uvar->ULminimumDeltaf)) uvar->ULminimumDeltaf = uvar->dfmin; //Upper limit minimum modulation depth (Hz)
   if (!XLALUserVarWasSet(&uvar->ULmaximumDeltaf)) uvar->ULmaximumDeltaf = uvar->dfmax; //Upper limit maximum modulation depth (Hz)

   //Check SSE/AVX settings
   if (uvar->vectorMath>2 || uvar->vectorMath<0) XLAL_ERROR(XLAL_FAILURE, "Must specify vectorMath to be 0, 1, or 2");

   //Developer options
   if (uvar->templateTest && uvar->bruteForceTemplateTest) XLAL_ERROR(XLAL_FAILURE, "Specify one of templateTest or bruteForceTemplateTest\n");
   if ((uvar->templateTest || uvar->bruteForceTemplateTest) || XLALUserVarWasSet(&uvar->templateTestF) || XLALUserVarWasSet(&uvar->templateTestP) || XLALUserVarWasSet(&uvar->templateTestDf)) {
      if (!((uvar->templateTest || uvar->bruteForceTemplateTest) && XLALUserVarWasSet(&uvar->templateTestF) && XLALUserVarWasSet(&uvar->templateTestP) && XLALUserVarWasSet(&uvar->templateTestDf))) {
         XLAL_ERROR(XLAL_FAILURE, "Must specify all or none of templateTest/bruteForceTemplateTest, templateTestF, templateTestP, and templateTestDf\n");
      }
   }

   //SFT input conflicts with injRandSeed option
   if (XLALUserVarWasSet(&uvar->injRandSeed) && XLALUserVarWasSet(&uvar->inputSFTs) && !uvar->gaussNoiseWithSFTgaps) XLAL_ERROR(XLAL_EINVAL, "When specifying injRandSeed, inputSFTs cannot be used unless specifying gaussNoiseWithSFTgaps\n");

   return XLAL_SUCCESS;

}


/**
 * Print REAL4Vector values to an ASCII file
 * \param [in] vector   Pointer to the REAL4Vector to be printed
 * \param [in] filename CHAR array of the file path and name
 * \return Status value
 */
INT4 printREAL4Vector2File(REAL4Vector *vector, CHAR *filename)
{

   FILE *OUTPUT = NULL;
   XLAL_CHECK( (OUTPUT = fopen(filename, "w")) != NULL, XLAL_EIO, "Couldn't fopen %s for writing.\n", filename );
   for (UINT4 ii=0; ii<vector->length; ii++) fprintf(OUTPUT, "%g\n", vector->data[ii]);
   fclose(OUTPUT);

   return XLAL_SUCCESS;

}
