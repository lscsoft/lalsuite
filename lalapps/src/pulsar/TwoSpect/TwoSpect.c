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
 * \ingroup pulsarApps
 * \author Evan Goetz
 */

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <lal/UserInput.h>
#include <lal/LALString.h>
#include <lal/Window.h>
#include <lal/DopplerScan.h>
#include <lal/VectorOps.h>
#include <lal/CWMakeFakeData.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>

#include <lalapps.h>

#include "IHS.h"
#include "candidates.h"
#include "antenna.h"
#include "templates.h"
#include "TwoSpect.h"
#include "statistics.h"
#include "upperlimits.h"
#include "vectormath.h"


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
   inputParamsStruct *inputParams = NULL;
   XLAL_CHECK( (inputParams = new_inputParams()) != NULL, XLAL_EFUNC );
   XLAL_CHECK ( readTwoSpectInputParams(inputParams, &uvar, argc, argv) == XLAL_SUCCESS, XLAL_EFUNC );

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

   //Print the output directory
   fprintf(stderr, "Output directory: %s\n", uvar.outdirectory);
   fprintf(LOG, "Output directory: %s\n", uvar.outdirectory);

   for (ii=0; ii<inputParams->numofIFOs; ii++) {
      fprintf(LOG,"IFO = %s\n", inputParams->detectors->sites[ii].frDetector.prefix);
      fprintf(stderr,"IFO = %s\n", inputParams->detectors->sites[ii].frDetector.prefix);
   }

   //Print to stderr the parameters of the search
   fprintf(stderr,"Tobs = %f sec\n",inputParams->Tobs);
   fprintf(stderr,"Tcoh = %f sec\n",inputParams->Tcoh);
   fprintf(stderr,"SFToverlap = %f sec\n",inputParams->SFToverlap);
   fprintf(stderr,"fmin = %f Hz\n",inputParams->fmin);
   fprintf(stderr,"fspan = %f Hz\n",inputParams->fspan);
   fprintf(stderr,"Pmin = %f s\n",inputParams->Pmin);
   fprintf(stderr,"Pmax = %f s\n",inputParams->Pmax);
   fprintf(stderr,"dfmin = %f Hz\n",inputParams->dfmin);
   fprintf(stderr,"dfmax = %f Hz\n",inputParams->dfmax);
   fprintf(stderr,"Running median blocksize = %d\n",inputParams->blksize);
   fprintf(stderr,"FFT plan flag = %d\n", inputParams->FFTplanFlag);
   fprintf(stderr,"RNG seed: %d\n", inputParams->randSeed);

   //Initialize ephemeris data structure
   EphemerisData *edat = NULL;
   XLAL_CHECK( (edat = XLALInitBarycenter(uvar.ephemEarth, uvar.ephemSun)) != NULL, XLAL_EFUNC );

   //Maximum detector velocity in units of c from start of observation time - Tcoh to end of observation + Tcoh
   REAL4 Vmax = 0.0;
   for (ii=0; ii<inputParams->numofIFOs; ii++) {
      REAL4 detectorVmax = CompDetectorVmax(inputParams->searchstarttime-inputParams->Tcoh, inputParams->Tcoh, inputParams->SFToverlap, inputParams->Tobs+2.0*inputParams->Tcoh, inputParams->detectors->sites[ii], edat);
      XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC, "CompDetectorVmax() failed\n" );
      if (detectorVmax > Vmax) Vmax = detectorVmax;
   }

   //Assume maximum bin shift possible
   inputParams->maxbinshift = (INT4)round(Vmax * (inputParams->fmin+inputParams->fspan) * inputParams->Tcoh) + 1;

   //Parameters for the sky-grid from a point/polygon or a sky-grid file
   DopplerSkyScanInit XLAL_INIT_DECL(scanInit);
   DopplerSkyScanState XLAL_INIT_DECL(scan);
   PulsarDopplerParams dopplerpos;
   if (XLALUserVarWasSet(&uvar.skyRegion)) {
      scanInit.gridType = GRID_ISOTROPIC;     //Default value for an approximate-isotropic grid
      scanInit.skyRegionString = uvar.skyRegion;      //"allsky" = Default value for all-sky search
      scanInit.numSkyPartitions = 1;   //Default value so sky is not broken into chunks
      scanInit.Freq = inputParams->fmin+0.5*inputParams->fspan;  //Midpoint of the frequency band
      scanInit.dAlpha = 0.5/((inputParams->fmin+inputParams->fspan) * inputParams->Tcoh * Vmax);
      scanInit.dDelta = scanInit.dAlpha;
      fprintf(LOG, "Sky region = %s\n", uvar.skyRegion);
      fprintf(stderr, "Sky region = %s\n", uvar.skyRegion);
   } else {
      scanInit.gridType = GRID_FILE_SKYGRID;
      scanInit.skyGridFile = uvar.skyRegionFile;
      scanInit.numSkyPartitions = 1;   //Default value so sky is not broken into chunks
      scanInit.Freq = inputParams->fmin+0.5*inputParams->fspan;  //Midpoint of the frequency band
      fprintf(LOG, "Sky file = %s\n", uvar.skyRegionFile);
      fprintf(stderr, "Sky file = %s\n", uvar.skyRegionFile);
   }

   //Initialize the sky-grid
   InitDopplerSkyScan(&status, &scan, &scanInit);
   XLAL_CHECK( status.statusCode == 0, XLAL_EFUNC );

   //Start at first location
   XLAL_CHECK( XLALNextDopplerSkyPos(&dopplerpos, &scan) == XLAL_SUCCESS, XLAL_EFUNC );

   //Basic units
   REAL4 tempfspan = inputParams->fspan + 2.0*inputParams->dfmax + (inputParams->blksize-1 + 12)/inputParams->Tcoh;     //= fspan+2*dfmax+extrabins + running median blocksize-1 (Hz)
   INT4 tempnumfbins = (INT4)round(tempfspan*inputParams->Tcoh)+1;                        //= number of bins in tempfspan
   fprintf(LOG, "FAR for templates = %g\n", inputParams->templatefar);
   fprintf(stderr, "FAR for templates = %g\n", inputParams->templatefar);

   //Allocate memory for ffdata structure
   ffdataStruct *ffdata = NULL;
   XLAL_CHECK( (ffdata = new_ffdata(inputParams)) != NULL, XLAL_EFUNC );

   //Maximum number of IHS values to sum = twice the maximum modulation depth
   //Minimum number of IHS values to sum = twice the minimum modulation depth
   INT4 maxrows = (INT4)round(2.0*inputParams->dfmax*inputParams->Tcoh)+1;
   fprintf(LOG, "Maximum row width to be searched = %d\n", maxrows);
   fprintf(stderr, "Maximum row width to be searched = %d\n", maxrows);

   //TF normalization
   ffdata->tfnormalization = 2.0/inputParams->Tcoh/(inputParams->avesqrtSh*inputParams->avesqrtSh);

   //Read in the T-F data from SFTs
   REAL4Vector *tfdata = NULL;
   MultiSFTVector *multiSFTvector = NULL;
   if (!XLALUserVarWasSet(&uvar.injectionSources) && XLALUserVarWasSet(&uvar.inputSFTs) && inputParams->numofIFOs<2 && !uvar.gaussNoiseWithSFTgaps && !XLALUserVarWasSet(&uvar.timestampsFile)) {
      XLAL_CHECK( (tfdata = readInSFTs(inputParams, ffdata->tfnormalization)) != NULL, XLAL_EFUNC );
   } else if (!XLALUserVarWasSet(&uvar.injectionSources) && XLALUserVarWasSet(&uvar.inputSFTs) && !uvar.gaussNoiseWithSFTgaps && !XLALUserVarWasSet(&uvar.timestampsFile)) {
      XLAL_CHECK( (multiSFTvector = getMultiSFTVector(inputParams)) != NULL, XLAL_EFUNC );
   } else {
      //Start by getting timestamps or creating them
      MultiLIGOTimeGPSVector *multiTimestamps = NULL;
      if (XLALUserVarWasSet(&uvar.inputSFTs)) {
         SFTCatalog *catalog = NULL;
         XLAL_CHECK( (catalog = findSFTdata(inputParams)) != NULL, XLAL_EFUNC );
         MultiSFTVector *tmpsftvector = NULL;
         XLAL_CHECK( (tmpsftvector = extractSFTband(inputParams, catalog)) != NULL, XLAL_EFUNC );
         XLAL_CHECK( (multiTimestamps = getMultiTimeStampsFromSFTs(tmpsftvector, inputParams)) != NULL, XLAL_EFUNC );
         XLALDestroyMultiSFTVector(tmpsftvector);
         XLALDestroySFTCatalog(catalog);
      } else if (XLALUserVarWasSet(&uvar.timestampsFile)) {
         XLAL_CHECK( (multiTimestamps = XLALReadMultiTimestampsFiles(uvar.timestampsFile)) != NULL, XLAL_EFUNC );
         for (ii=0; ii<(INT4)multiTimestamps->length; ii++) multiTimestamps->data[ii]->deltaT = inputParams->Tcoh;
      }
      else if (XLALUserVarWasSet(&uvar.segmentFile)) XLAL_CHECK( (multiTimestamps = getMultiTimeStampsFromSegmentsFile(uvar.segmentFile, inputParams)) != NULL, XLAL_EFUNC );
      else {
         LIGOTimeGPS tStart;
         XLALGPSSetREAL8 ( &tStart, inputParams->searchstarttime );
         XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC, "XLALGPSSetREAL8 failed\n" );
         XLAL_CHECK( (multiTimestamps = XLALMakeMultiTimestamps(tStart, inputParams->Tobs, inputParams->Tcoh, inputParams->SFToverlap, inputParams->numofIFOs)) != NULL, XLAL_EFUNC );
      }

      //TwoSpect to analyze:
      REAL8 TwoSpectFmin = round(inputParams->fmin*inputParams->Tcoh - inputParams->dfmax*inputParams->Tcoh - 0.5*(inputParams->blksize-1) - (REAL8)(inputParams->maxbinshift) - 6.0)/inputParams->Tcoh;
      REAL8 TwoSpectBand = round(inputParams->fspan*inputParams->Tcoh + 2.0*inputParams->dfmax*inputParams->Tcoh + (inputParams->blksize-1) + (REAL8)(2.0*inputParams->maxbinshift) + 12.0)/inputParams->Tcoh;

      //Setup the MFD data parameters
      CWMFDataParams DataParams;
      if (XLALUserVarWasSet(&uvar.injFmin) && XLALUserVarWasSet(&uvar.injBand) && uvar.injFmin<=TwoSpectFmin && uvar.injFmin+uvar.injBand>=TwoSpectFmin+TwoSpectBand) {
         DataParams.fMin = uvar.injFmin;
         DataParams.Band = uvar.injBand;
      } else {
         DataParams.fMin = TwoSpectFmin;
         DataParams.Band = TwoSpectBand;
      }
      DataParams.multiIFO.length = inputParams->numofIFOs;
      for (ii=0; ii<inputParams->numofIFOs; ii++) DataParams.multiIFO.sites[ii] = inputParams->detectors->sites[ii];
      DataParams.multiNoiseFloor.length = inputParams->numofIFOs;
      for (ii=0; ii<inputParams->numofIFOs; ii++) DataParams.multiNoiseFloor.sqrtSn[ii] = 0.0;
      DataParams.multiTimestamps = *multiTimestamps;
      DataParams.randSeed = uvar.injRandSeed;
      DataParams.SFTWindowType = "Hann";
      DataParams.SFTWindowBeta = 0;

      MultiSFTVector *signalSFTs = NULL;
      PulsarParamsVector *injectionSources = NULL;
      //If injection sources then read them and make signal sfts
      if (XLALUserVarWasSet(&uvar.injectionSources)) {
         XLAL_CHECK( (injectionSources =  XLALPulsarParamsFromUserInput(uvar.injectionSources)) != NULL, XLAL_EFUNC );

         fprintf(stderr, "Injecting %d signals... ", (INT4)injectionSources->length);

         if (!inputParams->signalOnly) XLAL_CHECK( XLALCWMakeFakeMultiData(&signalSFTs, NULL, injectionSources, &DataParams, edat) == XLAL_SUCCESS, XLAL_EFUNC );
         else XLAL_CHECK( XLALCWMakeFakeMultiData(&multiSFTvector, NULL, injectionSources, &DataParams, edat) == XLAL_SUCCESS, XLAL_EFUNC );

         if (DataParams.fMin != TwoSpectFmin) XLAL_CHECK( XLALMultiSFTVectorResizeBand(signalSFTs, TwoSpectFmin, TwoSpectBand) == XLAL_SUCCESS, XLAL_EFUNC );

         fprintf(stderr, "done\n");
      } // if there are injections

      //If not signal only, create sfts that include noise or extract a band from real data
      if (!inputParams->signalOnly) {
         if (uvar.gaussNoiseWithSFTgaps || XLALUserVarWasSet(&uvar.timestampsFile) || XLALUserVarWasSet(&uvar.segmentFile) || !XLALUserVarWasSet(&uvar.inputSFTs)) {
            for (ii=0; ii<inputParams->numofIFOs; ii++) DataParams.multiNoiseFloor.sqrtSn[ii] = inputParams->avesqrtSh;
            XLAL_CHECK( XLALCWMakeFakeMultiData(&multiSFTvector, NULL, NULL, &DataParams, edat) == XLAL_SUCCESS, XLAL_EFUNC );
            if (DataParams.fMin != TwoSpectFmin) XLAL_CHECK( XLALMultiSFTVectorResizeBand(multiSFTvector, TwoSpectFmin, TwoSpectBand) == XLAL_SUCCESS, XLAL_EFUNC );
         } else {
            SFTCatalog *catalog = NULL;
            XLAL_CHECK( (catalog = findSFTdata(inputParams)) != NULL, XLAL_EFUNC );
            XLAL_CHECK( (multiSFTvector = extractSFTband(inputParams, catalog)) != NULL, XLAL_EFUNC );
            XLALDestroySFTCatalog(catalog);
         }
      } // if not signal only SFTs

      //Add the SFT vectors together
      if (XLALUserVarWasSet(&uvar.injectionSources) && !inputParams->signalOnly) {
         XLAL_CHECK( XLALMultiSFTVectorAdd(multiSFTvector, signalSFTs) == XLAL_SUCCESS, XLAL_EFUNC );
         XLALDestroyMultiSFTVector(signalSFTs);
      }

      //Convert SFTs to powers
      if (inputParams->numofIFOs==1) XLAL_CHECK( (tfdata = convertSFTdataToPowers(multiSFTvector->data[0], inputParams, ffdata->tfnormalization)) != NULL, XLAL_EFUNC );

      //If printing the data outputs, then do that here
      if ((XLALUserVarWasSet(&uvar.printSignalData) || XLALUserVarWasSet(&uvar.printMarginalizedSignalData)) && XLALUserVarWasSet(&uvar.injectionSources)) {
         for (ii=0; ii<inputParams->numofIFOs; ii++) DataParams.multiNoiseFloor.sqrtSn[ii] = 0.0;
         PulsarParamsVector *oneSignal = NULL;
         XLAL_CHECK( (oneSignal = XLALCreatePulsarParamsVector(1)) != NULL, XLAL_EFUNC );

         FILE *SIGNALOUT = NULL, *MARGINALIZEDSIGNALOUT = NULL;
         if (XLALUserVarWasSet(&uvar.printSignalData)) XLAL_CHECK( (SIGNALOUT = fopen(uvar.printSignalData, "w")) != NULL, XLAL_EIO, "Failed to open %s for writing\n", uvar.printSignalData );
         if (XLALUserVarWasSet(&uvar.printMarginalizedSignalData)) XLAL_CHECK( (MARGINALIZEDSIGNALOUT = fopen(uvar.printMarginalizedSignalData, "w")) != NULL, XLAL_EIO, "Failed to open %s for writing\n", uvar.printMarginalizedSignalData );

         for (ii=0; ii<(INT4)injectionSources->length; ii++) {
            memcpy(oneSignal->data, &(injectionSources->data[ii]), sizeof(injectionSources->data[0]));
            if (XLALUserVarWasSet(&uvar.printSignalData)) {
               MultiSFTVector *oneSignalSFTs = NULL;
               XLAL_CHECK( XLALCWMakeFakeMultiData(&oneSignalSFTs, NULL, oneSignal, &DataParams, edat) == XLAL_SUCCESS, XLAL_EFUNC );
               if (DataParams.fMin != TwoSpectFmin) XLAL_CHECK( XLALMultiSFTVectorResizeBand(oneSignalSFTs, TwoSpectFmin, TwoSpectBand) == XLAL_SUCCESS, XLAL_EFUNC );

               REAL8Vector *aveSFTsPower = NULL;
               XLAL_CHECK( (aveSFTsPower = XLALCreateREAL8Vector(multiSFTvector->data[0]->data->data->length)) != NULL, XLAL_EFUNC );
               memset(aveSFTsPower->data, 0, sizeof(REAL8)*aveSFTsPower->length);

               for (jj=0; jj<(INT4)oneSignalSFTs->data[0]->length; jj++) {
                  SFTtype *sft = &(oneSignalSFTs->data[0]->data[jj]);
                  for (INT4 kk=0; kk<(INT4)aveSFTsPower->length; kk++) {
                     REAL8 powerval = 2.0*(creal(sft->data->data[kk])*creal(sft->data->data[kk]) + cimag(sft->data->data[kk])*cimag(sft->data->data[kk]))/inputParams->Tcoh;
                     aveSFTsPower->data[kk] += powerval;
                  }
               }
               for (jj=0; jj<(INT4)aveSFTsPower->length; jj++) fprintf(SIGNALOUT,"%.9g %.9g\n", DataParams.fMin+jj/inputParams->Tcoh, aveSFTsPower->data[jj]/multiSFTvector->data[0]->length);
               XLALDestroyMultiSFTVector(oneSignalSFTs);
               XLALDestroyREAL8Vector(aveSFTsPower);
            }
            if (XLALUserVarWasSet(&uvar.printMarginalizedSignalData)) {
               REAL8Vector *marginalizedSignalData = NULL;
               XLAL_CHECK ( (marginalizedSignalData = XLALCreateREAL8Vector(multiSFTvector->data[0]->data->data->length)) != NULL, XLAL_EFUNC );
               memset(marginalizedSignalData->data, 0, sizeof(REAL8)*marginalizedSignalData->length);
               for (jj=0; jj<300; jj++) {
                  oneSignal->data[0].Amp.cosi = 2.0*gsl_rng_uniform(inputParams->rng) - 1.0;
                  oneSignal->data[0].Amp.psi = LAL_TWOPI*gsl_rng_uniform(inputParams->rng);
                  oneSignal->data[0].Amp.phi0 = LAL_TWOPI*gsl_rng_uniform(inputParams->rng);
                  oneSignal->data[0].Doppler.argp = LAL_TWOPI*gsl_rng_uniform(inputParams->rng);
                  MultiSFTVector *oneSignalSFTs = NULL;
                  XLAL_CHECK( XLALCWMakeFakeMultiData(&oneSignalSFTs, NULL, oneSignal, &DataParams, edat) == XLAL_SUCCESS, XLAL_EFUNC );
                  if (DataParams.fMin != TwoSpectFmin) XLAL_CHECK( XLALMultiSFTVectorResizeBand(oneSignalSFTs, TwoSpectFmin, TwoSpectBand) == XLAL_SUCCESS, XLAL_EFUNC );

                  for (INT4 kk=0; kk<(INT4)oneSignalSFTs->data[0]->length; kk++) {
                     SFTtype *sft = &(oneSignalSFTs->data[0]->data[kk]);
                     for (INT4 ll=0; ll<(INT4)marginalizedSignalData->length; ll++) marginalizedSignalData->data[ll] += (2.0*(creal(sft->data->data[ll])*creal(sft->data->data[ll]) + cimag(sft->data->data[ll])*cimag(sft->data->data[ll]))/inputParams->Tcoh);
                  }
                  XLALDestroyMultiSFTVector(oneSignalSFTs);
               } //Loop over trials
               for (jj=0; jj<(INT4)marginalizedSignalData->length; jj++) {
                  marginalizedSignalData->data[jj] /= 300.0*multiSFTvector->data[0]->length;
                  fprintf(MARGINALIZEDSIGNALOUT,"%.9g %.9g\n", DataParams.fMin+jj/inputParams->Tcoh, marginalizedSignalData->data[jj]);
               }
               XLALDestroyREAL8Vector(marginalizedSignalData);
            } //If printing marginalized data
         } //loop over the number of injected sources
         memset(oneSignal->data, 0, sizeof(injectionSources->data[0]));
         XLALDestroyPulsarParamsVector(oneSignal);
         if (XLALUserVarWasSet(&uvar.printSignalData)) fclose(SIGNALOUT);
         if (XLALUserVarWasSet(&uvar.printMarginalizedSignalData)) fclose(MARGINALIZEDSIGNALOUT);
      } //end printing data

      if (XLALUserVarWasSet(&uvar.injectionSources)) XLALDestroyPulsarParamsVector(injectionSources);
      XLALDestroyMultiTimestamps(multiTimestamps);
      if (inputParams->numofIFOs==1) XLALDestroyMultiSFTVector(multiSFTvector);
   } // end load data or generate data

   //Print SFT times, if requested by user
   if (uvar.printSFTtimes) {
      CHARVector *outputfile = NULL;
      XLAL_CHECK( (outputfile = XLALCreateCHARVector(strlen(uvar.outdirectory)+25)) != NULL, XLAL_EFUNC );
      for (ii=0; ii<inputParams->numofIFOs; ii++) {
         sprintf(outputfile->data, "%s/%s-%s", uvar.outdirectory, uvar.IFO->data[ii], "inputSFTtimes.dat");
         FILE *INSFTTIMES = NULL;
         XLAL_CHECK( (INSFTTIMES = fopen(outputfile->data, "w")) != NULL, XLAL_EIO, "Couldn't fopen %s for writing input SFT start times", outputfile->data );
         if (inputParams->numofIFOs==1) {
            INT4 sftlength = tfdata->length/ffdata->numffts;
            for (jj=0; jj<ffdata->numffts; jj++) if (tfdata->data[jj*sftlength]!=0.0) fprintf(INSFTTIMES, "%9d 0\n", (INT4)round(inputParams->searchstarttime+jj*(inputParams->Tcoh-inputParams->SFToverlap)));
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
   if (XLALUserVarWasSet(&uvar.lineDetection) && !inputParams->signalOnly && inputParams->numofIFOs==1) {
      lines = detectLines_simple(tfdata, ffdata, inputParams);
      XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
      if (lines!=NULL) {
         fprintf(LOG, "WARNING: %d line(s) found.\n", lines->length);
         fprintf(stderr, "WARNING: %d line(s) found.\n", lines->length);
         if ((REAL4)lines->length/(ffdata->numfbins + 2*inputParams->maxbinshift) >= 0.1) {
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
   if (inputParams->numofIFOs==1) {
      XLAL_CHECK( (background = XLALCreateREAL4Vector(ffdata->numffts*(ffdata->numfbins+2*inputParams->maxbinshift))) != NULL, XLAL_EFUNC );
      if (!inputParams->signalOnly) XLAL_CHECK( tfRngMeans(background, tfdata, ffdata->numffts, ffdata->numfbins+2*inputParams->maxbinshift, inputParams->blksize) == XLAL_SUCCESS, XLAL_EFUNC );
      else memset(background->data, 0, sizeof(REAL4)*background->length);
   }

   //TEST: Try cleaning lines
   /* XLAL_CHECK( cleanLines(tfdata, background, lines, inputParams) == XLAL_SUCCESS, XLAL_EFUNC ); */
   /* if (lines!=NULL) { */
   /*    if (!inputParams->signalOnly) XLAL_CHECK( tfRngMeans(background, tfdata, ffdata->numffts, ffdata->numfbins + 2*inputParams->maxbinshift, inputParams->blksize) == XLAL_SUCCESS, XLAL_EFUNC ); */
   /*    else memset(background->data, 0, sizeof(REAL4)*background->length); */
   /* } */
   /* XLALDestroyINT4Vector(lines); */
   /* lines = NULL; */

   //Existing SFTs listed in this vector
   INT4Vector *sftexist = NULL;
   INT4Vector *indexValuesOfExistingSFTs = NULL;
   REAL4 frac_tobs_complete = 1.0;
   if (!inputParams->signalOnly && inputParams->numofIFOs==1) {
      XLAL_CHECK( (sftexist = existingSFTs(tfdata, inputParams, ffdata->numfbins, ffdata->numffts)) != NULL, XLAL_EFUNC );
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
   if (!inputParams->signalOnly && inputParams->numofIFOs==1) {
      REAL4 harmonicMean = 0.0;
      XLAL_CHECK( calcHarmonicMean(&harmonicMean, background, ffdata->numfbins + 2*inputParams->maxbinshift, ffdata->numffts) == XLAL_SUCCESS, XLAL_EFUNC );
      backgroundmeannormfactor = 1.0/harmonicMean;
      ffdata->tfnormalization *= backgroundmeannormfactor;
   }

   //Need to reduce the original TF data to remove the excess bins used for running median calculation. Normalize the TF as the same as the background was normalized
   REAL4Vector *usableTFdata = NULL;
   if (inputParams->numofIFOs==1) {
      XLAL_CHECK( (usableTFdata = XLALCreateREAL4Vector(background->length)) != NULL, XLAL_EFUNC );
      for (ii=0; ii<ffdata->numffts; ii++) memcpy(&(usableTFdata->data[ii*(ffdata->numfbins+2*inputParams->maxbinshift)]), &(tfdata->data[ii*(tempnumfbins+2*inputParams->maxbinshift) + (INT4)round(0.5*(inputParams->blksize-1))]), sizeof(REAL4)*(ffdata->numfbins+2*inputParams->maxbinshift));
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
   if (uvar.printData && inputParams->numofIFOs==1) {
      CHARVector *outputfile = NULL;
      XLAL_CHECK( (outputfile = XLALCreateCHARVector(strlen(uvar.outdirectory)+15)) != NULL, XLAL_EFUNC );
      sprintf(outputfile->data, "%s/%s", uvar.outdirectory, "tfdata.dat");
      XLAL_CHECK( printREAL4Vector2File(usableTFdata, outputfile->data) == XLAL_SUCCESS, XLAL_EFUNC );
      XLALDestroyCHARVector(outputfile);
   }

   //Do mean subtraction of TFdata here--modifies the usableTFdata vector!!!
   if (inputParams->numofIFOs==1) XLAL_CHECK( tfMeanSubtract(usableTFdata, background, ffdata->numffts, ffdata->numfbins+2*inputParams->maxbinshift) == XLAL_SUCCESS, XLAL_EFUNC );

   //Initialize reused values
   ihsMaximaStruct *ihsmaxima = NULL;
   ihsfarStruct *ihsfarstruct = NULL;
   REAL4Vector *aveNoise = NULL;
   INT4Vector *binshifts = NULL;
   LIGOTimeGPS tStart;
   MultiLIGOTimeGPSVector *multiTimestamps = NULL;
   MultiDetectorStateSeries *multiStateSeries = NULL;
   candidateVector *gaussCandidates1 = NULL, *gaussCandidates2 = NULL, *gaussCandidates3 = NULL, *gaussCandidates4 = NULL, *exactCandidates1 = NULL, *exactCandidates2 = NULL, *ihsCandidates = NULL;
   UpperLimitVector *upperlimits = NULL;
   REAL4FFTPlan *secondFFTplan = NULL;
   XLAL_CHECK( (ihsmaxima = new_ihsMaxima(ffdata->numfbins, maxrows)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihsfarstruct = new_ihsfarStruct(maxrows, inputParams)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (aveNoise = XLALCreateREAL4Vector(ffdata->numfprbins)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (binshifts = XLALCreateINT4Vector(ffdata->numffts)) != NULL, XLAL_EFUNC );
   XLALGPSSetREAL8 ( &tStart, inputParams->searchstarttime );
   XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC, "XLALGPSSetREAL8 failed\n" );
   XLAL_CHECK( (multiTimestamps = XLALMakeMultiTimestamps(tStart, inputParams->Tobs, inputParams->Tcoh, inputParams->SFToverlap, inputParams->numofIFOs)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (multiStateSeries = XLALGetMultiDetectorStates(multiTimestamps, inputParams->detectors, edat, inputParams->SFToverlap)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (gaussCandidates1 = new_candidateVector(100)) != NULL, XLAL_EFUNC, "new_CandidateVector(%d) failed\n", 100 );
   XLAL_CHECK( (gaussCandidates2 = new_candidateVector(100)) != NULL, XLAL_EFUNC, "new_CandidateVector(%d) failed\n", 100 );
   XLAL_CHECK( (gaussCandidates3 = new_candidateVector(100)) != NULL, XLAL_EFUNC, "new_CandidateVector(%d) failed\n", 100 );
   XLAL_CHECK( (gaussCandidates4 = new_candidateVector(100)) != NULL, XLAL_EFUNC, "new_CandidateVector(%d) failed\n", 100 );
   XLAL_CHECK( (exactCandidates1 = new_candidateVector(100)) != NULL, XLAL_EFUNC, "new_CandidateVector(%d) failed\n", 100 );
   XLAL_CHECK( (exactCandidates2 = new_candidateVector(100)) != NULL, XLAL_EFUNC, "new_CandidateVector(%d) failed\n", 100 );
   XLAL_CHECK( (ihsCandidates = new_candidateVector(100)) != NULL, XLAL_EFUNC, "new_CandidateVector(%d) failed\n", 100 );
   XLAL_CHECK( (upperlimits = new_UpperLimitVector(1)) != NULL, XLAL_EFUNC, "new_UpperLimitVector(%d) failed.\n", 1 );
   XLAL_CHECK( (secondFFTplan = XLALCreateForwardREAL4FFTPlan(ffdata->numffts, inputParams->FFTplanFlag)) != NULL, XLAL_EFUNC );
   LIGOTimeGPS refTime = multiTimestamps->data[0]->data[0];

   //Initialize to -1.0 for far just at the start
   ihsfarstruct->ihsfar->data[0] = -1.0;
   REAL4 antweightsrms = 0.0;

   //Davies' algorithm uses an internal error code
   //Later upgrade: remove this because we do our own error handling
   INT4 proberrcode = 0;

   //Set TF normalization
   if (inputParams->numofIFOs==1) ffdata->tfnormalization *= 0.5*inputParams->Tcoh;

   //Antenna normalization (determined from injections on H1 at ra=0, dec=0, with circular polarization)
   //When doing linear polarizations, the IHS factor needs to be 25.2*1.082 and this antenna weights
   //function needs to be set to use linear polarization.
   REAL4Vector *antweightsforihs2h0 = NULL;
   XLAL_CHECK( (antweightsforihs2h0 = XLALCreateREAL4Vector(ffdata->numffts)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( CompAntennaPatternWeights(antweightsforihs2h0, 0.0, 0.0, inputParams->searchstarttime, inputParams->Tcoh, inputParams->SFToverlap, inputParams->Tobs, 0, 0.0, lalCachedDetectors[LAL_LHO_4K_DETECTOR]) == XLAL_SUCCESS, XLAL_EFUNC );

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
      REAL4Vector *antweights = NULL;
      XLAL_CHECK( (antweights = XLALCreateREAL4Vector(ffdata->numffts)) != NULL, XLAL_EFUNC );
      memset(antweights->data, 0, antweights->length*sizeof(REAL4));
      if (inputParams->antennaOff) for (ii=0; ii<(INT4)antweights->length; ii++) antweights->data[ii] = 1.0;
      else {
         if (XLALUserVarWasSet(&uvar.linPolAngle)) XLAL_CHECK( CompAntennaPatternWeights(antweights, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, inputParams->searchstarttime, inputParams->Tcoh, inputParams->SFToverlap, inputParams->Tobs, 1, uvar.linPolAngle, inputParams->detectors->sites[0]) == XLAL_SUCCESS, XLAL_EFUNC );
         else XLAL_CHECK( CompAntennaPatternWeights(antweights, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, inputParams->searchstarttime, inputParams->Tcoh, inputParams->SFToverlap, inputParams->Tobs, 0, 0.0, inputParams->detectors->sites[0]) == XLAL_SUCCESS, XLAL_EFUNC );
      }

      //Compute the bin shifts for each SFT
      XLAL_CHECK( CompBinShifts2(binshifts, multissb->data[0], inputParams->fmin + 0.5*inputParams->fspan, inputParams->Tcoh, inputParams->dopplerMultiplier) == XLAL_SUCCESS, XLAL_EFUNC );

      //Coherent combination
      if (inputParams->numofIFOs>=2) {
         //Get multiSSB times and multiAMcoefficients
         MultiAMCoeffs *multiAMcoefficients = NULL;
         XLAL_CHECK( (multiAMcoefficients = XLALComputeMultiAMCoeffs(multiStateSeries, NULL, skypos)) != NULL, XLAL_EFUNC );

         //Coherenly combine the SFT data
         if (XLALUserVarWasSet(&uvar.assumeNScosi) && XLALUserVarWasSet(&uvar.assumeNSpsi)) XLAL_CHECK( (tfdata = coherentlyAddSFTs(multiSFTvector, multissb, multiAMcoefficients, &uvar.assumeNScosi, &uvar.assumeNSpsi, inputParams)) != NULL, XLAL_EFUNC );
         else if (XLALUserVarWasSet(&uvar.assumeNScosi) && !XLALUserVarWasSet(&uvar.assumeNSpsi)) XLAL_CHECK( (tfdata = coherentlyAddSFTs(multiSFTvector, multissb, multiAMcoefficients, &uvar.assumeNScosi, NULL, inputParams)) != NULL, XLAL_EFUNC );
         else if (!XLALUserVarWasSet(&uvar.assumeNScosi) && XLALUserVarWasSet(&uvar.assumeNSpsi)) XLAL_CHECK( (tfdata = coherentlyAddSFTs(multiSFTvector, multissb, multiAMcoefficients, NULL, &uvar.assumeNSpsi, inputParams)) != NULL, XLAL_EFUNC );
         else XLAL_CHECK( (tfdata = coherentlyAddSFTs(multiSFTvector, multissb, multiAMcoefficients, NULL, NULL, inputParams)) != NULL, XLAL_EFUNC );

         //Continue like there is just one detector

         //Check for lines
         if (XLALUserVarWasSet(&uvar.lineDetection) && !inputParams->signalOnly) {
            lines = detectLines_simple(tfdata, ffdata, inputParams);
            XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
            if (lines!=NULL) {
               fprintf(LOG, "WARNING: %d line(s) found.\n", lines->length);
               fprintf(stderr, "WARNING: %d line(s) found.\n", lines->length);
               if ((REAL4)lines->length/(ffdata->numfbins + 2*inputParams->maxbinshift) >= 0.1) {
                  heavilyContaminatedBand = 1;
                  fprintf(LOG, "WARNING: Band is heavily contaminated by artifacts.\n");
                  fprintf(stderr, "WARNING: Band is heavily contaminated by artifacts.\n");
               }
            }
         }
         if (heavilyContaminatedBand) uvar.IHSonly = 1;

         //Compute background
         XLAL_CHECK( (background = XLALCreateREAL4Vector(ffdata->numffts*(ffdata->numfbins+2*inputParams->maxbinshift))) != NULL, XLAL_EFUNC );
         if (!inputParams->signalOnly) XLAL_CHECK( tfRngMeans(background, tfdata, ffdata->numffts, ffdata->numfbins+2*inputParams->maxbinshift, inputParams->blksize) == XLAL_SUCCESS, XLAL_EFUNC );
         else memset(background->data, 0, sizeof(REAL4)*background->length);

         if (!inputParams->signalOnly) {
            XLAL_CHECK( (sftexist = existingSFTs(tfdata, inputParams, ffdata->numfbins, ffdata->numffts)) != NULL, XLAL_EFUNC );
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
            XLAL_CHECK( calcHarmonicMean(&harmonicMean, background, ffdata->numfbins + 2*inputParams->maxbinshift, ffdata->numffts) == XLAL_SUCCESS, XLAL_EFUNC );
            backgroundmeannormfactor = 1.0/harmonicMean;
            ffdata->tfnormalization *= backgroundmeannormfactor;
         }

         XLAL_CHECK( (usableTFdata = XLALCreateREAL4Vector(background->length)) != NULL, XLAL_EFUNC );
         for (ii=0; ii<ffdata->numffts; ii++) memcpy(&(usableTFdata->data[ii*(ffdata->numfbins+2*inputParams->maxbinshift)]), &(tfdata->data[ii*(tempnumfbins+2*inputParams->maxbinshift) + (INT4)round(0.5*(inputParams->blksize-1))]), sizeof(REAL4)*(ffdata->numfbins+2*inputParams->maxbinshift));
         for (ii=0; ii<(INT4)usableTFdata->length; ii++) {
            if (usableTFdata->data[ii]!=0.0) {
               usableTFdata->data[ii] *= backgroundmeannormfactor;
               background->data[ii] *= backgroundmeannormfactor;
            }
         }
         XLALDestroyREAL4Vector(tfdata);

         XLAL_CHECK( tfMeanSubtract(usableTFdata, background, ffdata->numffts, ffdata->numfbins+2*inputParams->maxbinshift) == XLAL_SUCCESS, XLAL_EFUNC );

         ffdata->tfnormalization *= 0.5*inputParams->Tcoh;
      }
      /////

      //Track identified lines
      REAL4VectorSequence *trackedlines = NULL;
      if (lines!=NULL) XLAL_CHECK( (trackedlines = trackLines(lines, binshifts, inputParams)) != NULL, XLAL_EFUNC );

      //Calculate antenna RMS value
      REAL4 currentAntWeightsRMS = 0.0;
      XLAL_CHECK( calcRms(&currentAntWeightsRMS, antweights) == XLAL_SUCCESS, XLAL_EFUNC );

      //Slide SFTs here -- need to slide the data and the estimated background
      REAL4Vector *TFdata_slided = NULL, *background_slided = NULL;
      XLAL_CHECK( (TFdata_slided = XLALCreateREAL4Vector(ffdata->numffts*ffdata->numfbins)) != NULL, XLAL_EFUNC );
      XLAL_CHECK( (background_slided = XLALCreateREAL4Vector(TFdata_slided->length)) != NULL, XLAL_EFUNC );
      XLAL_CHECK( slideTFdata(TFdata_slided, inputParams, usableTFdata, binshifts) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( slideTFdata(background_slided, inputParams, background, binshifts) == XLAL_SUCCESS, XLAL_EFUNC );

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
      XLAL_CHECK( ffPlaneNoise(aveNoise, inputParams, sftexist, background_slided, antweightsforihs2h0, secondFFTplan, &(skypointffnormalization)) == XLAL_SUCCESS, XLAL_EFUNC );

      //Average noise floor of FF plane for each 1st FFT frequency bin
      ffdata->ffnormalization = 1.0;
      XLAL_CHECK( ffPlaneNoise(aveNoise, inputParams, sftexist, background_slided, antweights, secondFFTplan, &(ffdata->ffnormalization)) == XLAL_SUCCESS, XLAL_EFUNC );
      if (uvar.printData) {
         CHAR w[1000];
         snprintf(w, 1000, "%s/%s", uvar.outdirectory, "ffbackground.dat");
         XLAL_CHECK( printREAL4Vector2File(aveNoise, w) == XLAL_SUCCESS, XLAL_EFUNC );
      }

      //Compute the weighted TF data
      REAL4Vector *TFdata_weighted = NULL;
      XLAL_CHECK( (TFdata_weighted = XLALCreateREAL4Vector(ffdata->numffts*ffdata->numfbins)) != NULL, XLAL_EFUNC );
      XLAL_CHECK( tfWeight(TFdata_weighted, TFdata_slided, background_slided, antweights, indexValuesOfExistingSFTs, inputParams) == XLAL_SUCCESS, XLAL_EFUNC );
      XLALDestroyREAL4Vector(TFdata_slided);
      XLALDestroyREAL4Vector(antweights);

      //Print out data product if requested
      if (uvar.printData) {
         CHAR w[1000];
         snprintf(w, 1000, "%s/%s", uvar.outdirectory, "procTFdata.dat");
         XLAL_CHECK( printREAL4Vector2File(TFdata_weighted, w) == XLAL_SUCCESS, XLAL_EFUNC );
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
      for (ii=0; ii<(INT4)ffdata->ffdata->length; ii++) ffdata->ffdata->data[ii] *= inputParams->Tobs;

      REAL4 secFFTmean = 0.0, secFFTsigma = 0.0;
      secFFTmean = calcMean(ffdata->ffdata);
      XLAL_CHECK( calcStddev(&secFFTsigma, ffdata->ffdata) == XLAL_SUCCESS, XLAL_EFUNC );

      XLALDestroyREAL4Vector(TFdata_weighted);
      TFdata_weighted = NULL;

      //Print out data product if requested
      if (uvar.printData) {
         CHAR w[1000];
         snprintf(w, 1000, "%s/%s", uvar.outdirectory, "ffdata.dat");
         XLAL_CHECK( printREAL4Vector2File(ffdata->ffdata, w) == XLAL_SUCCESS, XLAL_EFUNC );
      }

      fprintf(LOG, "2nd FFT ave = %g, 2nd FFT stddev = %g, expected ave = 1.0\n", secFFTmean, secFFTsigma);
      fprintf(stderr, "2nd FFT ave = %g, 2nd FFT stddev = %g, expected ave = 1.0\n", secFFTmean, secFFTsigma);

      //Exit with failure if there are no SFTs (probably this doesn't get hit)
      XLAL_CHECK( secFFTmean != 0.0, XLAL_FAILURE, "Average second FFT power is 0.0. Perhaps no SFTs are remaining? Program exiting with failure.\n" );

      //If the user wants to test a single, exact template, then we do that here
      if (uvar.templateTest) {
         if (uvar.printData) fprintf(stderr, "numfbins=%d, maxbinshift=%d, numffts=%d, numfprbins=%d\n", ffdata->numfbins, inputParams->maxbinshift, ffdata->numffts, ffdata->numfprbins);
         fprintf(stderr, "Testing template f=%f, P=%f, Df=%f\n", uvar.templateTestF, uvar.templateTestP, uvar.templateTestDf);
         fprintf(LOG, "Testing template f=%f, P=%f, Df=%f\n", uvar.templateTestF, uvar.templateTestP, uvar.templateTestDf);

         //Load template quantities into a test candidate
         loadCandidateData(&(exactCandidates1->data[0]), uvar.templateTestF, uvar.templateTestP, uvar.templateTestDf, dopplerpos.Alpha, dopplerpos.Delta, 0.0, 0.0, 0.0, 0, 0.0);

         //Resize the output candidate vector if necessary
         if (exactCandidates2->numofcandidates == exactCandidates2->length-1) XLAL_CHECK( (exactCandidates2 = resize_candidateVector(exactCandidates2, 2*exactCandidates2->length)) != NULL, XLAL_EFUNC );

         //Analyze the template stored in the test candidate
         XLAL_CHECK( analyzeOneTemplate(&(exactCandidates2->data[exactCandidates2->numofcandidates]), &(exactCandidates1->data[0]), ffdata, aveNoise, aveTFnoisePerFbinRatio, inputParams, sftexist, secondFFTplan) == XLAL_SUCCESS, XLAL_EFUNC );
         (exactCandidates2->numofcandidates)++;

         //Rescale the h0 output from the normaliztions and amount of observation time present
         exactCandidates2->data[exactCandidates2->numofcandidates-1].h0 /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25);

      } else if (uvar.bruteForceTemplateTest) {
         candidate cand;
         loadCandidateData(&cand, uvar.templateTestF, uvar.templateTestP, uvar.templateTestDf, dopplerpos.Alpha, dopplerpos.Delta, 0.0, 0.0, 0.0, 0, 0.0);
         XLAL_CHECK( bruteForceTemplateTest(&(exactCandidates2), cand, uvar.templateTestF-2.0/inputParams->Tcoh, uvar.templateTestF+2.0/inputParams->Tcoh, 21, 5, 5, 0.5, uvar.templateTestDf-2.0/inputParams->Tcoh, uvar.templateTestDf+2.0/inputParams->Tcoh, 21, inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 1) == XLAL_SUCCESS, XLAL_EFUNC );
         for (ii=0; ii<(INT4)exactCandidates2->numofcandidates; ii++) exactCandidates2->data[ii].h0 /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25);
      }

      //If the user wants to do a template search, that is done here
      if (uvar.templateSearch) {

         //fprintf(stderr, "Calling templateSearch\n (in development, last edited 2014-06-09)\n");
         XLAL_CHECK( templateSearch_scox1Style(&exactCandidates2, inputParams->fmin, inputParams->fspan, uvar.templateSearchP, uvar.templateSearchAsini, uvar.templateSearchAsiniSigma, dopplerpos.Alpha, dopplerpos.Delta, inputParams, ffdata->ffdata, sftexist, aveNoise,  aveTFnoisePerFbinRatio,  secondFFTplan, 1) == XLAL_SUCCESS, XLAL_EFUNC );
         //fprintf(stderr, "Done calling templateSearch\n");


         for (ii=0; ii<(INT4)exactCandidates2->numofcandidates; ii++) exactCandidates2->data[ii].h0 /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25);
      }

      if (inputParams->signalOnly) return 0;

////////Start of the IHS step!
      candidateVector *ihsCandidates_reduced = NULL;
      //Find the FAR of IHS sum -- only if the templateTest has not been given
      if (!uvar.templateTest && !uvar.templateSearch && !uvar.bruteForceTemplateTest) {
         //If the false alarm thresholds need to be computed
         if (ihsfarstruct->ihsfar->data[0]<0.0) XLAL_CHECK( genIhsFar(ihsfarstruct, inputParams, maxrows, aveNoise) == XLAL_SUCCESS, XLAL_EFUNC );

         //Run the IHS algorithm on the data
         XLAL_CHECK( runIHS(ihsmaxima, ffdata, ihsfarstruct, inputParams, maxrows, aveNoise, aveTFnoisePerFbinRatio) == XLAL_SUCCESS, XLAL_EFUNC );

         //Find any IHS candidates
         XLAL_CHECK( findIHScandidates(&ihsCandidates, ihsfarstruct, inputParams, ffdata, ihsmaxima, aveTFnoisePerFbinRatio, trackedlines) == XLAL_SUCCESS, XLAL_EFUNC );

         //If requested, keep only the most significant IHS candidates
         if (XLALUserVarWasSet(&uvar.keepOnlyTopNumIHS) && (INT4)ihsCandidates->numofcandidates>uvar.keepOnlyTopNumIHS) {
            XLAL_CHECK( (ihsCandidates_reduced = keepMostSignificantCandidates(ihsCandidates, inputParams)) != NULL, XLAL_EFUNC );

            //Put ihsCandidates_reduced back into a reset ihsCandidates
            ihsCandidates->numofcandidates = 0;
            for (ii=0; ii<(INT4)ihsCandidates_reduced->numofcandidates; ii++) {
               loadCandidateData(&(ihsCandidates->data[ii]), ihsCandidates_reduced->data[ii].fsig, ihsCandidates_reduced->data[ii].period, ihsCandidates_reduced->data[ii].moddepth, dopplerpos.Alpha, dopplerpos.Delta, ihsCandidates_reduced->data[ii].stat, ihsCandidates_reduced->data[ii].h0, ihsCandidates_reduced->data[ii].prob, 0, ihsCandidates_reduced->data[ii].normalization);
               (ihsCandidates->numofcandidates)++;
            }

            free_candidateVector(ihsCandidates_reduced);
         }
      }
////////End of the IHS step

////////Start of the Gaussian template search!
      //First check to see if the IHSonly or templateTest or templateSearch was given
      if (uvar.IHSonly && !uvar.templateTest && !uvar.templateSearch) {
         //Check the length of the exactCandidates2 vector is large enough and resize if necessary
         if (exactCandidates2->length < exactCandidates2->numofcandidates+ihsCandidates->numofcandidates) XLAL_CHECK( (exactCandidates2 = resize_candidateVector(exactCandidates2, exactCandidates2->numofcandidates+ihsCandidates->numofcandidates)) != NULL, XLAL_EFUNC );

         //Use the typical list
         INT4 numofcandidatesalready = exactCandidates2->numofcandidates;
         for (ii=0; ii<(INT4)ihsCandidates->numofcandidates; ii++) {
            loadCandidateData(&(exactCandidates2->data[ii+numofcandidatesalready]), ihsCandidates->data[ii].fsig, ihsCandidates->data[ii].period, ihsCandidates->data[ii].moddepth, dopplerpos.Alpha, dopplerpos.Delta, ihsCandidates->data[ii].stat, ihsCandidates->data[ii].h0, ihsCandidates->data[ii].prob, 0, ihsCandidates->data[ii].normalization);
            exactCandidates2->data[ii+numofcandidatesalready].h0 /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25); //Scaling here
            (exactCandidates2->numofcandidates)++;
         }

      } else if (!uvar.templateTest && !uvar.templateSearch) {

         //Test the IHS candidates against Gaussian templates in this function
         XLAL_CHECK( testIHScandidates(&gaussCandidates1, ihsCandidates, ffdata, aveNoise, aveTFnoisePerFbinRatio, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, inputParams) == XLAL_SUCCESS, XLAL_EFUNC );

         for (ii=0; ii<(INT4)gaussCandidates1->numofcandidates; ii++) fprintf(stderr, "Candidate %d: f0=%g, P=%g, df=%g\n", ii, gaussCandidates1->data[ii].fsig, gaussCandidates1->data[ii].period, gaussCandidates1->data[ii].moddepth);
      } /* if IHSonly is not given && templateTest not given and templateSearch not given */
////////End of the Gaussian template search

      //Reset IHS candidates, but keep length the same (doesn't reset actual values in the vector)
      ihsCandidates->numofcandidates = 0;

      //Search the candidates further if the number of candidates passing the first Gaussian template test is greater than 0
      if (gaussCandidates1->numofcandidates>0) {
////////Start clustering! Note that the clustering algorithm takes care of the period range of parameter space
         XLAL_CHECK( clusterCandidates(&gaussCandidates2, gaussCandidates1, ffdata, inputParams, aveNoise, aveTFnoisePerFbinRatio, sftexist, 0) == XLAL_SUCCESS, XLAL_EFUNC );

         for (ii=0; ii<(INT4)gaussCandidates2->numofcandidates; ii++) fprintf(stderr, "Candidate %d: f0=%g, P=%g, df=%g\n", ii, gaussCandidates2->data[ii].fsig, gaussCandidates2->data[ii].period, gaussCandidates2->data[ii].moddepth);
////////End clustering

         //Reset first set of Gaussian template candidates, but keep length the same (doesn't reset actual values in the vector)
         gaussCandidates1->numofcandidates = 0;

////////Start detailed Gaussian template search!
         for (ii=0; ii<(INT4)gaussCandidates2->numofcandidates; ii++) {

            if (gaussCandidates3->numofcandidates == gaussCandidates3->length-1) XLAL_CHECK( (gaussCandidates3 = resize_candidateVector(gaussCandidates3, 2*gaussCandidates3->length)) != NULL, XLAL_EFUNC );
            //bruteForceTemplateSearch(&(gaussCandidates3->data[gaussCandidates3->numofcandidates]), gaussCandidates2->data[ii], gaussCandidates2->data[ii].fsig-2.5/inputParams->Tcoh, gaussCandidates2->data[ii].fsig+2.5/inputParams->Tcoh, 11, 5, gaussCandidates2->data[ii].moddepth-2.5/inputParams->Tcoh, gaussCandidates2->data[ii].moddepth+2.5/inputParams->Tcoh, 11, inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 0);
            //bruteForceTemplateSearch(&(gaussCandidates3->data[gaussCandidates3->numofcandidates]), gaussCandidates2->data[ii], gaussCandidates2->data[ii].fsig-1.5/inputParams->Tcoh, gaussCandidates2->data[ii].fsig+1.5/inputParams->Tcoh, 7, 5, gaussCandidates2->data[ii].moddepth-1.5/inputParams->Tcoh, gaussCandidates2->data[ii].moddepth+1.5/inputParams->Tcoh, 7, inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 0);
            XLAL_CHECK( bruteForceTemplateSearch(&(gaussCandidates3->data[gaussCandidates3->numofcandidates]), //Output candidate
                                                 gaussCandidates2->data[ii],                             //Candidate
                                                 gaussCandidates2->data[ii].fsig-1.0/inputParams->Tcoh,  //Minimum frequency
                                                 gaussCandidates2->data[ii].fsig+1.0/inputParams->Tcoh,  //Maximum frequency
                                                 5,                                                      //Number of frequencies to search in range
                                                 2,                                                      //Number of longer periods to search
                                                 2,                                                      //Number of shorter periods to search
                                                 1.0,                                                    //Period spacing factor
                                                 gaussCandidates2->data[ii].moddepth-1.0/inputParams->Tcoh, //Minimum modulation depth
                                                 gaussCandidates2->data[ii].moddepth+1.0/inputParams->Tcoh, //Maximum modulation depth
                                                 5,                                                      //Number of modulation depths to search in range
                                                 inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 0) == XLAL_SUCCESS, XLAL_EFUNC );
            gaussCandidates3->numofcandidates++;

         } /* for ii < numofcandidates */

         for (ii=0; ii<(INT4)gaussCandidates3->numofcandidates; ii++) fprintf(stderr,"Candidate %d: f0=%g, P=%g, df=%g\n", ii, gaussCandidates3->data[ii].fsig, gaussCandidates3->data[ii].period, gaussCandidates3->data[ii].moddepth);
////////End detailed Gaussian template search

         //Reset 2nd round of Gaussian template candidates, but keep length the same (doesn't reset actual values in the vector)
         gaussCandidates2->numofcandidates = 0;

////////Start clustering!
         XLAL_CHECK( clusterCandidates(&gaussCandidates4, gaussCandidates3, ffdata, inputParams, aveNoise, aveTFnoisePerFbinRatio, sftexist, 0) == XLAL_SUCCESS, XLAL_EFUNC );

         for (ii=0; ii<(INT4)gaussCandidates4->numofcandidates; ii++) fprintf(stderr, "Candidate %d: f0=%g, P=%g, df=%g\n", ii, gaussCandidates4->data[ii].fsig, gaussCandidates4->data[ii].period, gaussCandidates4->data[ii].moddepth);
////////End clustering

         //Reset 3rd set of Gaussian template candidates
         gaussCandidates3->numofcandidates = 0;

////////Initial check using "exact" template
         for (ii=0; ii<(INT4)gaussCandidates4->numofcandidates; ii++) {
            //Allocate the template memory
            templateStruct *template = NULL;
            XLAL_CHECK( (template = new_templateStruct(inputParams->maxtemplatelength)) != NULL, XLAL_EFUNC );

            //Produce either the Gaussian template or exact template
            if (!uvar.gaussTemplatesOnly) XLAL_CHECK( makeTemplate(template, gaussCandidates4->data[ii], inputParams, sftexist, secondFFTplan) == XLAL_SUCCESS, XLAL_EFUNC );
            else XLAL_CHECK( makeTemplateGaussians(template, gaussCandidates4->data[ii], inputParams, ffdata->numfbins, ffdata->numfprbins) == XLAL_SUCCESS, XLAL_EFUNC );

            //Produce the FAR threshold, if requested
            farStruct *farval = NULL;
            if (inputParams->calcRthreshold) {
               XLAL_CHECK( (farval = new_farStruct()) != NULL, XLAL_EFUNC );
               XLAL_CHECK( numericFAR(farval, template, inputParams->templatefar, aveNoise, aveTFnoisePerFbinRatio, inputParams, inputParams->rootFindingMethod) == XLAL_SUCCESS, XLAL_EFUNC );
            }

            //Calculate R, false alarm probability and h0 estimate
            REAL8 R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
            XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
            REAL8 h0 = 0.0;
            REAL8 prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, inputParams, &proberrcode);
            XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
            if ( R > 0.0 ) h0 = 2.7426*pow(R/(inputParams->Tcoh*inputParams->Tobs),0.25);

            //Check that we are above threshold before storing the candidate
            if ((!inputParams->calcRthreshold && prob<inputParams->log10templatefar) || (inputParams->calcRthreshold && R>farval->far)) {
               if (exactCandidates1->numofcandidates == exactCandidates1->length-1) XLAL_CHECK( (exactCandidates1 = resize_candidateVector(exactCandidates1, 2*exactCandidates1->length)) != NULL, XLAL_EFUNC );
               loadCandidateData(&exactCandidates1->data[exactCandidates1->numofcandidates], gaussCandidates4->data[ii].fsig, gaussCandidates4->data[ii].period, gaussCandidates4->data[ii].moddepth, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, R, h0, prob, proberrcode, gaussCandidates4->data[ii].normalization);
               exactCandidates1->numofcandidates++;
            }

            free_templateStruct(template);
            template = NULL;
            if (inputParams->calcRthreshold) {
               free_farStruct(farval);
               farval = NULL;
            }
         } /* for ii < numofcandidates */
         fprintf(LOG, "Number of candidates confirmed with exact templates = %d\n", exactCandidates1->numofcandidates);
         fprintf(stderr, "Number of candidates confirmed with exact templates = %d\n", exactCandidates1->numofcandidates);
         for (ii=0; ii<(INT4)exactCandidates1->numofcandidates; ii++) fprintf(stderr, "Candidate %d: f0=%g, P=%g, df=%g\n", ii, exactCandidates1->data[ii].fsig, exactCandidates1->data[ii].period, exactCandidates1->data[ii].moddepth);
////////Done with initial check

         //Reset 4th set of Gaussian template candidates, but keep length the same (doesn't reset actual values in the vector)
         gaussCandidates4->numofcandidates = 0;

////////Start detailed "exact" template search!
         for (ii=0; ii<(INT4)exactCandidates1->numofcandidates; ii++) {

            if (exactCandidates2->numofcandidates == exactCandidates2->length-1) XLAL_CHECK( (exactCandidates2 = resize_candidateVector(exactCandidates2, 2*exactCandidates2->length)) != NULL, XLAL_EFUNC );

            if (!uvar.gaussTemplatesOnly) {
               //bruteForceTemplateSearch(&(exactCandidates2->data[exactCandidates2->numofcandidates]), exactCandidates1->data[ii], exactCandidates1->data[ii].fsig-1.0/inputParams->Tcoh, exactCandidates1->data[ii].fsig+1.0/inputParams->Tcoh, 5, 5, exactCandidates1->data[ii].moddepth-1.0/inputParams->Tcoh, exactCandidates1->data[ii].moddepth+1.0/inputParams->Tcoh, 5, inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 1);
               //bruteForceTemplateSearch(&(exactCandidates2->data[exactCandidates2->numofcandidates]), exactCandidates1->data[ii], exactCandidates1->data[ii].fsig-1.0/inputParams->Tcoh, exactCandidates1->data[ii].fsig+1.0/inputParams->Tcoh, 5, 3, exactCandidates1->data[ii].moddepth-1.0/inputParams->Tcoh, exactCandidates1->data[ii].moddepth+1.0/inputParams->Tcoh, 5, inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 1);
               XLAL_CHECK( bruteForceTemplateSearch(&(exactCandidates2->data[exactCandidates2->numofcandidates]),
                                                    exactCandidates1->data[ii],
                                                    exactCandidates1->data[ii].fsig-0.5/inputParams->Tcoh,
                                                    exactCandidates1->data[ii].fsig+0.5/inputParams->Tcoh,
                                                    3,
                                                    1,
                                                    1,
                                                    1.0,
                                                    exactCandidates1->data[ii].moddepth-0.5/inputParams->Tcoh,
                                                    exactCandidates1->data[ii].moddepth+0.5/inputParams->Tcoh,
                                                    3,
                                                    inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 1) == XLAL_SUCCESS, XLAL_EFUNC );
            } else {
               //bruteForceTemplateSearch(&(exactCandidates2->data[exactCandidates2->numofcandidates]), exactCandidates1->data[ii], exactCandidates1->data[ii].fsig-1.0/inputParams->Tcoh, exactCandidates1->data[ii].fsig+1.0/inputParams->Tcoh, 5, 5, exactCandidates1->data[ii].moddepth-1.0/inputParams->Tcoh, exactCandidates1->data[ii].moddepth+1.0/inputParams->Tcoh, 5, inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 0);
               //bruteForceTemplateSearch(&(exactCandidates2->data[exactCandidates2->numofcandidates]), exactCandidates1->data[ii], exactCandidates1->data[ii].fsig-1.0/inputParams->Tcoh, exactCandidates1->data[ii].fsig+1.0/inputParams->Tcoh, 5, 3, exactCandidates1->data[ii].moddepth-1.0/inputParams->Tcoh, exactCandidates1->data[ii].moddepth+1.0/inputParams->Tcoh, 5, inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 0);
               XLAL_CHECK( bruteForceTemplateSearch(&(exactCandidates2->data[exactCandidates2->numofcandidates]),
                                                    exactCandidates1->data[ii],
                                                    exactCandidates1->data[ii].fsig-0.5/inputParams->Tcoh,
                                                    exactCandidates1->data[ii].fsig+0.5/inputParams->Tcoh,
                                                    3,
                                                    1,
                                                    1,
                                                    1.0,
                                                    exactCandidates1->data[ii].moddepth-0.5/inputParams->Tcoh,
                                                    exactCandidates1->data[ii].moddepth+0.5/inputParams->Tcoh,
                                                    3,
                                                    inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 0) == XLAL_SUCCESS, XLAL_EFUNC );
            }
            exactCandidates2->data[exactCandidates2->numofcandidates].h0 /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25);  //Scaling here
            exactCandidates2->numofcandidates++;

            fprintf(stderr,"Candidate %d: f0=%g, P=%g, df=%g\n", ii, exactCandidates2->data[ii].fsig, exactCandidates2->data[ii].period, exactCandidates2->data[ii].moddepth);
         } /* for ii < numofcandidates */
////////End of detailed search

         //Reset first round of exact template candidates, but keep length the same (doesn't reset actual values in the vector)
         exactCandidates1->numofcandidates = 0;

         fprintf(LOG,"Exact step is done with the total number of candidates = %d\n", exactCandidates2->numofcandidates);
         fprintf(stderr,"Exact step is done with the total number of candidates = %d\n", exactCandidates2->numofcandidates);

      } /* if gaussCandidates1->numofcandidates > 0 */

      //Determine upper limits, if the ULoff has not been set
      if (!uvar.ULoff && !uvar.templateTest && !uvar.templateSearch) {
         upperlimits->data[upperlimits->length-1].alpha = (REAL4)dopplerpos.Alpha;
         upperlimits->data[upperlimits->length-1].delta = (REAL4)dopplerpos.Delta;
         upperlimits->data[upperlimits->length-1].normalization = ffdata->tfnormalization;

         XLAL_CHECK( skypoint95UL(&(upperlimits->data[upperlimits->length-1]), inputParams, ffdata, ihsmaxima, ihsfarstruct, aveTFnoisePerFbinRatio) == XLAL_SUCCESS, XLAL_EFUNC );

         for (ii=0; ii<(INT4)upperlimits->data[upperlimits->length-1].ULval->length; ii++) upperlimits->data[upperlimits->length-1].ULval->data[ii] /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25);  //Compensation for different duty cycle and antenna pattern weights

         XLAL_CHECK( (upperlimits = resize_UpperLimitVector(upperlimits, upperlimits->length+1)) != NULL, XLAL_EFUNC );
      } /* if producing UL */

      //Destroy stuff
      XLALDestroyREAL4Vector(aveTFnoisePerFbinRatio);
      XLALDestroyREAL4VectorSequence(trackedlines);

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
         XLAL_CHECK( outputUpperLimitToFile(ULFILENAME->data, upperlimits->data[ii], inputParams->printAllULvalues) == XLAL_SUCCESS, XLAL_EFUNC );
      }
   }

   //Destroy varaibles
   XLALDestroyREAL4Vector(antweightsforihs2h0);
   XLALDestroyREAL4Vector(background);
   XLALDestroyREAL4Vector(usableTFdata);
   XLALDestroyREAL4Vector(aveNoise);
   XLALDestroyINT4Vector(binshifts);
   XLALDestroyMultiTimestamps(multiTimestamps);
   XLALDestroyMultiDetectorStateSeries(multiStateSeries);
   XLALDestroyINT4Vector(sftexist);
   XLALDestroyINT4Vector(indexValuesOfExistingSFTs);
   free_ffdata(ffdata);
   free_ihsfarStruct(ihsfarstruct);
   free_inputParams(inputParams);
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

   if (inputParams->numofIFOs>1) XLALDestroyMultiSFTVector(multiSFTvector);
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
 * Create a new inputParamsStruct
 * \param [in] numofIFOs Number of interferometers
 * \return Pointer to new inputParamsStruct
 */
inputParamsStruct * new_inputParams(void)
{
   inputParamsStruct *input = NULL;
   XLAL_CHECK_NULL( (input = XLALMalloc(sizeof(*input))) != NULL, XLAL_ENOMEM );
   XLAL_CHECK_NULL( (input->detectors = XLALMalloc(sizeof(MultiLALDetector))) != NULL, XLAL_ENOMEM );
   XLAL_CHECK_NULL( (input->rng = gsl_rng_alloc(gsl_rng_mt19937)) != NULL, XLAL_EFUNC );
   return input;
} /* new_inputParams() */


/**
 * Free the inputParamsStruct
 * \param [in] params Pointer to the inputParamsStruct
 */
void free_inputParams(inputParamsStruct *params)
{
   gsl_rng_free(params->rng);
   if (params->inputSFTs) XLALFree(params->inputSFTs);
   XLALFree(params->detectors);
   XLALFree((inputParamsStruct*)params);
} /* free_inputParams() */



/**
 * Create a new frequency-frequency data structure for the TwoSpect analysis
 * \param [in] params Pointer to the inputParamsStruct
 * \return Pointer to new ffdataStruct
 */
ffdataStruct * new_ffdata(inputParamsStruct *params)
{

   XLAL_CHECK_NULL( params != NULL, XLAL_EINVAL );

   ffdataStruct *ffdata = NULL;
   XLAL_CHECK_NULL( (ffdata = XLALMalloc(sizeof(*ffdata))) != NULL, XLAL_ENOMEM );

   ffdata->numfbins = (INT4)(round(params->fspan*params->Tcoh + 2.0*params->dfmax*params->Tcoh)+12+1);
   ffdata->numffts = (INT4)floor(params->Tobs/(params->Tcoh-params->SFToverlap)-1);
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
 * Find the SFT data specified by user input
 * \param [in] params Poineter to the inputParamsStruct
 * \return Pointer to SFTCatalog
 */
SFTCatalog * findSFTdata(inputParamsStruct *params)
{

   XLAL_CHECK_NULL( params != NULL, XLAL_EINVAL );

   fprintf(LOG, "Finding SFTs... ");
   fprintf(stderr, "Finding SFTs... ");

   //Set the start and end times in the LIGO GPS format
   LIGOTimeGPS start = LIGOTIMEGPSZERO, end = LIGOTIMEGPSZERO;
   XLALGPSSetREAL8(&start, params->searchstarttime);
   XLAL_CHECK_NULL( xlalErrno == 0, XLAL_EFUNC );
   XLALGPSSetREAL8(&end, params->searchstarttime+params->Tobs);
   XLAL_CHECK_NULL( xlalErrno == 0, XLAL_EFUNC );

   //Setup the constraints
   SFTConstraints XLAL_INIT_DECL(constraints);
   if (params->numofIFOs == 1) constraints.detector = params->detectors->sites[0].frDetector.prefix;
   constraints.minStartTime = &start;
   constraints.maxStartTime = &end;

   //Find SFT files
   SFTCatalog *catalog = NULL;
   XLAL_CHECK_NULL( (catalog = XLALSFTdataFind(params->inputSFTs, &constraints)) != NULL, XLAL_EFUNC );

   fprintf(LOG, "done\n");
   fprintf(stderr, "done\n");

   return catalog;

}


/**
 * Extract the SFT coefficients from the band of interest
 * \param [in] params Pointer to inputParamsStruct
 * \param [in] catalog Pointer to the SFTCatalog
 * \return Pointer to MultiSFTVector
 */
MultiSFTVector * extractSFTband(inputParamsStruct *params, SFTCatalog *catalog)
{

   XLAL_CHECK_NULL( params != NULL && catalog != NULL, XLAL_EINVAL );

   fprintf(LOG, "Extracting band from SFTs... ");
   fprintf(stderr, "Extracting band from SFTs... ");

   //Determine band size (remember to get extra bins because of the running median and the bin shifts due to detector velocity)
   REAL8 minfbin = round(params->fmin*params->Tcoh - params->dfmax*params->Tcoh - 0.5*(params->blksize-1) - (REAL8)(params->maxbinshift) - 6.0)/params->Tcoh;
   REAL8 maxfbin = round((params->fmin + params->fspan)*params->Tcoh + params->dfmax*params->Tcoh + 0.5*(params->blksize-1) + (REAL8)(params->maxbinshift) + 6.0)/params->Tcoh;

   //Now extract the data
   MultiSFTVector *sftvector = NULL;
   XLAL_CHECK_NULL( (sftvector = XLALLoadMultiSFTs(catalog, minfbin+0.1/params->Tcoh, maxfbin-0.1/params->Tcoh)) != NULL, XLAL_EFUNC );

   XLAL_CHECK_NULL( params->numofIFOs == (INT4)sftvector->length, XLAL_FAILURE );

   fprintf(LOG, "done\n");
   fprintf(stderr, "done\n");

   return sftvector;

}


/**
 * Get a MultiSFTVector from the user-input values
 * \param [in] params Pointer to inputParamsStruct
 * \return Pointer to MultiSFTVector
 */
MultiSFTVector * getMultiSFTVector(inputParamsStruct *params)
{
   //Get catalog of SFTs
   SFTCatalog *catalog = NULL;
   XLAL_CHECK_NULL( (catalog = findSFTdata(params)) != NULL, XLAL_EFUNC );

   MultiSFTVector *sftvector = NULL;
   if (params->markBadSFTs && !params->signalOnly) {
      //Extract band to get a MultiSFTVector
      MultiSFTVector *tmpsftvector = NULL;
      XLAL_CHECK_NULL( (tmpsftvector = extractSFTband(params, catalog)) != NULL, XLAL_EFUNC );

      //Get the timestamps of the SFTs applying the KS/Kuipers tests if desired
      MultiLIGOTimeGPSVector *multiTimestamps = NULL;
      XLAL_CHECK_NULL( (multiTimestamps = getMultiTimeStampsFromSFTs(tmpsftvector, params)) != NULL, XLAL_EFUNC );

      //Get the SFT subset
      XLAL_CHECK_NULL( (sftvector = XLALExtractMultiSFTVectorWithMultiTimestamps(tmpsftvector, multiTimestamps)) != NULL, XLAL_EFUNC );

      XLALDestroyMultiSFTVector(tmpsftvector);
      XLALDestroyMultiTimestamps(multiTimestamps);
   } else {
      XLAL_CHECK_NULL( (sftvector = extractSFTband(params, catalog)) != NULL, XLAL_EFUNC );
   }

   XLALDestroySFTCatalog(catalog);

   return sftvector;
}


/**
 * Compute the Dirichlet kernel for large N values
 * \param [in] delta The delta value as the arguement
 * \return Complex valued Dirichlet kernel
 */
COMPLEX16 DirichletKernelLargeN(REAL8 delta)
{
   if (fabs(delta)<1.0e-6) return crect(1.0, 0.0);

   COMPLEX16 val = crect(0.0, LAL_TWOPI*delta);
   return (cexp(val)-1.0)/val;
}


/**
 * Compute the Dirichlet kernel for large N values and the Hann window
 * \param [in] delta The delta value as the arguement
 * \return Complex valued Dirichlet kernel
 */
COMPLEX16 DirichletKernelLargeNHann(REAL8 delta)
{
   return DirichletKernelLargeN(delta) - 0.5*DirichletKernelLargeN(delta+1.0) - 0.5*DirichletKernelLargeN(delta-1.0);
}


/**
 * Phase shift an SFT
 * \param [in,out] sft   Pointer to an SFT
 * \param [in]     shift Amount of the phase shift
 * \return Status value
 */
INT4 PhaseShiftSFT(SFTtype *sft, REAL8 shift)
{
   XLAL_CHECK ( (sft != NULL) && (sft->data != NULL), XLAL_EINVAL );

   if ( shift == 0 ) return XLAL_SUCCESS;

   for ( UINT4 k=0; k < sft->data->length; k++ ) {
      REAL8 fk = sft->f0 + k * sft->deltaF; /* frequency of k-th bin */
      REAL8 shiftCyles = shift * fk;

      COMPLEX8 fact = cexpf(crectf(0.0, (REAL4)(LAL_TWOPI*shiftCyles)));
      sft->data->data[k] *= fact;
   } /* for k < numBins */

   return XLAL_SUCCESS;
}


/**
 * Add SFTs together from a MultiSFTVector
 * \param [in] multiSFTVector      Pointer to a MultiSFTVector containing the SFT data
 * \param [in] multissb            Pointer to a MultiSSBtimes structure
 * \param [in] multiAMcoefficients Pointer to a MultiAMCoeffs structure
 * \param [in] assumeNScosi        Pointer to the assumed cosi value, or NULL if no value assumed
 * \param [in] assumeNSpsi         Pointer to the assumed psi value, or NULL if no value assumed
 * \param [in] params              Pointer to inputParamsStruct
 * \return REAL4Vector of powers of the coherently combined SFTs
 */
REAL4Vector * coherentlyAddSFTs(MultiSFTVector *multiSFTvector, MultiSSBtimes *multissb, MultiAMCoeffs *multiAMcoefficients, REAL8 *assumeNScosi, REAL8 *assumeNSpsi, inputParamsStruct *params) {

   //Sanity check the values
   XLAL_CHECK_NULL( multiSFTvector != NULL, XLAL_EINVAL );
   XLAL_CHECK_NULL( multissb != NULL, XLAL_EINVAL );
   XLAL_CHECK_NULL( multiAMcoefficients != NULL, XLAL_EINVAL );
   XLAL_CHECK_NULL( params != NULL, XLAL_EINVAL );

   fprintf(stderr, "Coherently adding SFT data... ");

   //Determine length of the largest SFTvector and create an SFTvector large enough to hold this
   UINT4 maxSFTs = multiSFTvector->data[0]->length;
   for (INT4 ii=1; ii<params->numofIFOs; ii++) {
      if ( multiSFTvector->data[ii]->length > maxSFTs ) {
         maxSFTs = multiSFTvector->data[ii]->length;
      }
   }
   SFTVector *combinedSFTs = NULL;
   XLAL_CHECK_NULL( (combinedSFTs = XLALCreateSFTVector(maxSFTs, 0)) != NULL, XLAL_EFUNC );

   //Create an INT4Vector to determine at which vector we are using in the multiSFTvector
   INT4Vector *whichSFTinMultiSFTvector = NULL, *whichIFOsToBeUsed = NULL;
   XLAL_CHECK_NULL( (whichSFTinMultiSFTvector = XLALCreateINT4Vector(params->numofIFOs)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (whichIFOsToBeUsed = XLALCreateINT4Vector(params->numofIFOs)) != NULL, XLAL_EFUNC );
   memset(whichSFTinMultiSFTvector->data, 0, whichSFTinMultiSFTvector->length * sizeof(INT4));  //All index values are set to zero at the beginning

   //Loop over the combinedSFTs vector to fill it with single or coherently combined SFTs
   for (INT4 ii=0; ii<(INT4)combinedSFTs->length; ii++) {
      //Loop over the interferometers, determining which ones to use
      memset(whichIFOsToBeUsed->data, 0, whichIFOsToBeUsed->length * sizeof(INT4));  //All index values are set to zero each time
      LIGOTimeGPS smallestGPS = multiSFTvector->data[0]->data[whichSFTinMultiSFTvector->data[0]].epoch;
      INT4 ifoWithSmallestGPS = 0;
      for (INT4 jj=1; jj<params->numofIFOs; jj++) {
         INT4 compareT = XLALGPSCmp(&(multiSFTvector->data[ifoWithSmallestGPS]->data[whichSFTinMultiSFTvector->data[ifoWithSmallestGPS]].epoch), &(multiSFTvector->data[jj]->data[whichSFTinMultiSFTvector->data[jj]].epoch));
         if (compareT>0) {
            smallestGPS = multiSFTvector->data[jj]->data[whichSFTinMultiSFTvector->data[jj]].epoch;
            ifoWithSmallestGPS = jj;
         }
      }
      for (INT4 jj=0; jj<params->numofIFOs; jj++) {
         INT4 compareT = XLALGPSCmp(&(smallestGPS), &(multiSFTvector->data[jj]->data[whichSFTinMultiSFTvector->data[jj]].epoch));
         if (compareT==0) whichIFOsToBeUsed->data[jj] = 1;
      }

      INT4 createSFT = 1;
      for (INT4 jj=0; jj<params->numofIFOs; jj++) {
         if (jj==0 && whichIFOsToBeUsed->data[jj]==1) {
            //Copy the data from the multiSFTvector into the combinedSFTs vector
            XLAL_CHECK_NULL( XLALCopySFT(&(combinedSFTs->data[ii]), &(multiSFTvector->data[jj]->data[whichSFTinMultiSFTvector->data[jj]])) == XLAL_SUCCESS, XLAL_EFUNC );
            createSFT = 0;
            whichSFTinMultiSFTvector->data[jj]++;
         } else if (jj>0 && whichIFOsToBeUsed->data[jj]==1) {
            //Create a copy of the SFT to be shifted since we will manipulate the SFT coefficients
            SFTtype *sftcopy = NULL;
            XLAL_CHECK_NULL( (sftcopy = XLALCreateSFT(0)) != NULL, XLAL_EFUNC );
            XLAL_CHECK_NULL( XLALCopySFT(sftcopy, &(multiSFTvector->data[jj]->data[whichSFTinMultiSFTvector->data[jj]])) == XLAL_SUCCESS, XLAL_EFUNC );
            REAL8 sftstart = XLALGPSGetREAL8(&(sftcopy->epoch));
            XLAL_CHECK_NULL( xlalErrno == 0, XLAL_EFUNC );
            INT4 fftnum = (INT4)((sftstart - params->searchstarttime)/params->SFToverlap);

            //First shift based on wavefront time of arrival
            REAL8 timediff0 = multissb->data[0]->DeltaT->data[fftnum] - 0.5*params->Tcoh*multissb->data[0]->Tdot->data[fftnum];
            REAL8 timediff = multissb->data[jj]->DeltaT->data[fftnum] - 0.5*params->Tcoh*multissb->data[jj]->Tdot->data[fftnum];
            //REAL8 tau = timediff0 - timediff;
            //XLAL_CHECK_NULL( PhaseShiftSFT(sftcopy, tau) == XLAL_SUCCESS, XLAL_EFUNC );
            REAL8 tau = timediff - timediff0;
            XLAL_CHECK_NULL( PhaseShiftSFT(sftcopy, -tau) == XLAL_SUCCESS, XLAL_EFUNC );

            //Average of detector-signal phase relation, unless assumed NS orientation
            //Need to compute arg[(C_1^1+iC_2^1)/(C_1^0+iC_2^0)] to determine an average value of exp(i*arg). Assume magnitude of ratio = 1.
            //This is simply arg[0.25*Fplus1*Fplus0*(1.0+cosi*cosi)*(1.0+cosi*cosi) + Fcross1*Fcross0*cosi*cosi +
            // i*(Fcross1*cosi*0.5*Fplus0*(1.0+cosi*cosi) - 0.5*Fplus1*(1.0+cosi*cosi)*Fcross0*cosi)]
            REAL4 detPhaseArg = 0.0, detPhaseMag = 1.0;
            COMPLEX16 detPhase = crect(0.0, 0.0);
            if (assumeNScosi==NULL && assumeNSpsi==NULL) {
               for (INT4 kk=0; kk<20; kk++) {
                  REAL4 psi = 0.05*kk*LAL_PI;
                  REAL4 Fplus0 = multiAMcoefficients->data[0]->a->data[ii]*cos(2.0*psi) + multiAMcoefficients->data[0]->b->data[ii]*sin(2.0*psi);
                  REAL4 Fcross0 = multiAMcoefficients->data[0]->b->data[ii]*cos(2.0*psi) - multiAMcoefficients->data[0]->a->data[ii]*sin(2.0*psi);
                  REAL4 Fplus1 = multiAMcoefficients->data[jj]->a->data[ii]*cos(2.0*psi) + multiAMcoefficients->data[jj]->b->data[ii]*sin(2.0*psi);
                  REAL4 Fcross1 = multiAMcoefficients->data[jj]->b->data[ii]*cos(2.0*psi) - multiAMcoefficients->data[jj]->a->data[ii]*sin(2.0*psi);
                  for (INT4 ll=0; ll<21; ll++) {
                     REAL4 cosi = 1.0 - 0.05*ll*2.0;
                     REAL4 realterm = 0.25*Fplus1*Fplus0*(1.0+cosi*cosi)*(1.0+cosi*cosi) + Fcross1*Fcross0*cosi*cosi;
                     REAL4 imagterm = Fcross1*cosi*0.5*Fplus0*(1.0+cosi*cosi) - 0.5*Fplus1*(1.0+cosi*cosi)*Fcross0*cosi;
                     COMPLEX16 complexval = crect(realterm, imagterm);
                     detPhase += cexp(crect(0.0, carg(complexval)));
                  }
               }
               detPhase /= 420.0;
               detPhaseArg = (REAL4)gsl_sf_angle_restrict_pos(carg(detPhase));
            } else if (assumeNScosi!=NULL && assumeNSpsi==NULL) {
               for (INT4 kk=0; kk<20; kk++) {
                  REAL4 psi = 0.05*kk*LAL_PI;
                  REAL4 Fplus0 = multiAMcoefficients->data[0]->a->data[ii]*cos(2.0*psi) + multiAMcoefficients->data[0]->b->data[ii]*sin(2.0*psi);
                  REAL4 Fcross0 = multiAMcoefficients->data[0]->b->data[ii]*cos(2.0*psi) - multiAMcoefficients->data[0]->a->data[ii]*sin(2.0*psi);
                  REAL4 Fplus1 = multiAMcoefficients->data[jj]->a->data[ii]*cos(2.0*psi) + multiAMcoefficients->data[jj]->b->data[ii]*sin(2.0*psi);
                  REAL4 Fcross1 = multiAMcoefficients->data[jj]->b->data[ii]*cos(2.0*psi) - multiAMcoefficients->data[jj]->a->data[ii]*sin(2.0*psi);
                  REAL4 cosi = *assumeNScosi;
                  REAL4 realterm = 0.25*Fplus1*Fplus0*(1.0+cosi*cosi)*(1.0+cosi*cosi) + Fcross1*Fcross0*cosi*cosi;
                  REAL4 imagterm = Fcross1*cosi*0.5*Fplus0*(1.0+cosi*cosi) - 0.5*Fplus1*(1.0+cosi*cosi)*Fcross0*cosi;
                  COMPLEX16 complexval = crect(realterm, imagterm);
                  detPhase += cexp(crect(0.0, carg(complexval)));
               }
               detPhase *= 0.05;
               detPhaseArg = (REAL4)gsl_sf_angle_restrict_pos(carg(detPhase));
            } else if (assumeNScosi==NULL && assumeNSpsi!=NULL) {
               REAL4 psi = *assumeNSpsi;
               REAL4 Fplus0 = multiAMcoefficients->data[0]->a->data[ii]*cos(2.0*psi) + multiAMcoefficients->data[0]->b->data[ii]*sin(2.0*psi);
               REAL4 Fcross0 = multiAMcoefficients->data[0]->b->data[ii]*cos(2.0*psi) - multiAMcoefficients->data[0]->a->data[ii]*sin(2.0*psi);
               REAL4 Fplus1 = multiAMcoefficients->data[jj]->a->data[ii]*cos(2.0*psi) + multiAMcoefficients->data[jj]->b->data[ii]*sin(2.0*psi);
               REAL4 Fcross1 = multiAMcoefficients->data[jj]->b->data[ii]*cos(2.0*psi) - multiAMcoefficients->data[jj]->a->data[ii]*sin(2.0*psi);
               for (INT4 kk=0; kk<21; kk++) {
                  REAL4 cosi = 1.0 - 0.05*kk*2.0;
                  REAL4 realterm = 0.25*Fplus1*Fplus0*(1.0+cosi*cosi)*(1.0+cosi*cosi) + Fcross1*Fcross0*cosi*cosi;
                  REAL4 imagterm = Fcross1*cosi*0.5*Fplus0*(1.0+cosi*cosi) - 0.5*Fplus1*(1.0+cosi*cosi)*Fcross0*cosi;
                  COMPLEX16 complexval = crect(realterm, imagterm);
                  detPhase += cexp(crect(0.0, carg(complexval)));
               }
               detPhase /= 21.0;
               detPhaseArg = (REAL4)gsl_sf_angle_restrict_pos(carg(detPhase));
            } else {
               REAL4 psi = *assumeNSpsi;
               REAL4 Fplus0 = multiAMcoefficients->data[0]->a->data[ii]*cos(2.0*psi) + multiAMcoefficients->data[0]->b->data[ii]*sin(2.0*psi);
               REAL4 Fcross0 = multiAMcoefficients->data[0]->b->data[ii]*cos(2.0*psi) - multiAMcoefficients->data[0]->a->data[ii]*sin(2.0*psi);
               REAL4 Fplus1 = multiAMcoefficients->data[jj]->a->data[ii]*cos(2.0*psi) + multiAMcoefficients->data[jj]->b->data[ii]*sin(2.0*psi);
               REAL4 Fcross1 = multiAMcoefficients->data[jj]->b->data[ii]*cos(2.0*psi) - multiAMcoefficients->data[jj]->a->data[ii]*sin(2.0*psi);
               REAL4 cosi = *assumeNScosi;
               COMPLEX16 RatioTerm = crect(0.5*Fplus1*(1.0+cosi*cosi), Fcross1*cosi)/crect(0.5*Fplus0*(1.0+cosi*cosi), Fcross0*cosi);
               detPhaseArg = (REAL4)gsl_sf_angle_restrict_pos(carg(RatioTerm));
               detPhaseMag = (REAL4)cabs(RatioTerm);
            }

            //Now for the frequency difference
            for (INT4 kk=0; kk<(INT4)sftcopy->data->length; kk++) {
               REAL8 delta0 = (sftcopy->f0+sftcopy->deltaF*kk)*params->Tcoh*(multissb->data[0]->Tdot->data[fftnum]-1.0);
               REAL8 delta = (sftcopy->f0+sftcopy->deltaF*kk)*params->Tcoh*(multissb->data[jj]->Tdot->data[fftnum]-1.0);
               COMPLEX16 DirichletRatio = conj(DirichletKernelLargeNHann(delta)/DirichletKernelLargeNHann(delta0));
               REAL8 detArgVal = gsl_sf_angle_restrict_pos(carg(DirichletRatio));
               if (llabs((INT8)delta0-(INT8)delta)>=1) {
                  detArgVal += LAL_PI;
               }
               COMPLEX8 phi = crectf(0.0, (REAL4)(detArgVal+detPhaseArg));
               sftcopy->data->data[kk] *= detPhaseMag*cexpf(phi);
            }

            if (createSFT==1) {
               XLAL_CHECK_NULL( XLALCopySFT(&(combinedSFTs->data[ii]), sftcopy) == XLAL_SUCCESS, XLAL_EFUNC );
               createSFT = 0;
            }
            else XLAL_CHECK_NULL( XLALSFTAdd(&(combinedSFTs->data[ii]), sftcopy) == XLAL_SUCCESS, XLAL_EFUNC );
            XLALDestroySFT(sftcopy);
            whichSFTinMultiSFTvector->data[jj]++;
         }
      }
   }
   XLALDestroyMultiSSBtimes(multissb);
   XLALDestroyMultiAMCoeffs(multiAMcoefficients);
   XLALDestroyINT4Vector(whichSFTinMultiSFTvector);
   XLALDestroyINT4Vector(whichIFOsToBeUsed);

   fprintf(stderr, "done\n");

   REAL4Vector *tfdata = NULL;
   XLAL_CHECK_NULL( (tfdata = convertSFTdataToPowers(combinedSFTs, params, 2.0/params->Tcoh/(params->avesqrtSh*params->avesqrtSh))) != NULL, XLAL_EFUNC );
   XLALDestroySFTVector(combinedSFTs);

   return tfdata;
}


/**
 * Convert the MultiSFTVector sfts into powers
 * \param [in] sfts Pointer to the MultiSFTVector of SFTs
 * \param [in] params Pointer to the inputParamsStruct
 * \param [in] normalization Normalization value to prevent underflow
 * \return Pointer to REAL4Vector containing powers
 */
REAL4Vector * convertSFTdataToPowers(SFTVector *sfts, inputParamsStruct *params, REAL8 normalization)
{

   XLAL_CHECK_NULL( sfts != NULL && params != NULL, XLAL_EINVAL );

   fprintf(LOG, "Converting band to powers... ");
   fprintf(stderr, "Converting band to powers... ");

   INT4 ii, jj;

   //Determine band size (remember to get extra bins because of the running median and the bin shifts due to detector velocity)
   REAL8 minfbin = round(params->fmin*params->Tcoh - params->dfmax*params->Tcoh - 0.5*(params->blksize-1) - (REAL8)(params->maxbinshift) - 6.0)/params->Tcoh;
   REAL8 maxfbin = round((params->fmin + params->fspan)*params->Tcoh + params->dfmax*params->Tcoh + 0.5*(params->blksize-1) + (REAL8)(params->maxbinshift) + 6.0)/params->Tcoh;

   //Now put the power data into the TF plane, looping through each SFT
   //If an SFT doesn't exit, fill the TF pixels of the SFT with zeros
   INT4 numffts = (INT4)floor(params->Tobs/(params->Tcoh-params->SFToverlap)-1);
   INT4 sftlength;
   if (sfts->data[0].data->length == 0) sftlength = (INT4)round(maxfbin*params->Tcoh - minfbin*params->Tcoh + 1);
   else {
      sftlength = sfts->data[0].data->length;
      //Check the length is what we expect
      XLAL_CHECK_NULL( sftlength==(INT4)round(maxfbin*params->Tcoh - minfbin*params->Tcoh + 1), XLAL_EFPINEXCT, "sftlength (%d) is not matching expected length (%d)\n", sftlength, (INT4)round(maxfbin*params->Tcoh - minfbin*params->Tcoh + 1) );
   }
   INT4 nonexistantsft = 0;
   REAL4Vector *tfdata = NULL;
   XLAL_CHECK_NULL( (tfdata = XLALCreateREAL4Vector(numffts*sftlength)) != NULL, XLAL_EFUNC );

   //Load the data into the output vector, roughly normalizing as we go along from the input value
   REAL8 sqrtnorm = sqrt(normalization);
   for (ii=0; ii<numffts; ii++) {
      if (ii-nonexistantsft < (INT4)sfts->length) {
         SFTtype *sft = &(sfts->data[ii - nonexistantsft]);
         if (sft->epoch.gpsSeconds == (INT4)round(ii*(params->Tcoh-params->SFToverlap)+params->searchstarttime)) {
            for (jj=0; jj<sftlength; jj++) {
               COMPLEX8 sftcoeff = sft->data->data[jj];
               tfdata->data[ii*sftlength + jj] = (REAL4)((sqrtnorm*crealf(sftcoeff))*(sqrtnorm*crealf(sftcoeff)) + (sqrtnorm*cimagf(sftcoeff))*(sqrtnorm*cimagf(sftcoeff)));  //power, normalized
            } /* for jj < sftLength */
         } else {
            memset(&(tfdata->data[ii*sftlength]), 0, sizeof(REAL4)*sftlength);
            nonexistantsft++;    //increment the nonexistantsft counter
         }
      } else {
         memset(&(tfdata->data[ii*sftlength]), 0, sizeof(REAL4)*sftlength);
         nonexistantsft++;    //increment the nonexistantsft counter
      }
   } /* for ii < numffts */

   fprintf(LOG, "done\n");
   fprintf(stderr, "done\n");

   fprintf(LOG, "Duty factor = %f\n", 1.0-(REAL4)nonexistantsft/(REAL4)numffts);
   fprintf(stderr, "Duty factor = %f\n", 1.0-(REAL4)nonexistantsft/(REAL4)numffts);

   REAL4 meanTFdata = calcMean(tfdata);
   REAL4 stddev = 0.0;
   XLAL_CHECK_NULL( calcStddev(&stddev, tfdata) == XLAL_SUCCESS, XLAL_EFUNC );
   fprintf(LOG, "TF before weighting, mean subtraction: mean = %g, std. dev. = %g\n", meanTFdata, stddev);
   fprintf(stderr, "TF before weighting, mean subtraction: mean = %g, std. dev. = %g\n", meanTFdata, stddev);

   return tfdata;

}


/**
 * Read in the data SFTs in one function
 * \param [in] params Pointer to inputParamsStruct
 * \param [in] normalization Normalization value determined from expected noise background
 * \return REAL4Vector of SFT powers
 */
REAL4Vector * readInSFTs(inputParamsStruct *params, REAL8 normalization)
{

   XLAL_CHECK_NULL( params != NULL, XLAL_EINVAL );

   MultiSFTVector *sftvector = NULL;
   XLAL_CHECK_NULL( (sftvector = getMultiSFTVector(params)) != NULL, XLAL_EFUNC );

   REAL4Vector *tfdata = NULL;
   XLAL_CHECK_NULL( (tfdata = convertSFTdataToPowers(sftvector->data[0], params, normalization)) != NULL, XLAL_EFUNC );

   //Destroy stuff
   XLALDestroyMultiSFTVector(sftvector);

   return tfdata;

} /* readInSFTs() */


/**
 * Create a list of timestamps from an SFTCatalog
 * \param [in] catalog Pointer to an SFTCatalog
 * \return Pointer to a list of GPS timestamps in a MultiLIGOTimeGPSVector
 */
MultiLIGOTimeGPSVector * getMultiTimeStampsFromSFTCatalog(SFTCatalog *catalog)
{

   XLAL_CHECK_NULL( catalog != NULL, XLAL_EINVAL );

   //Get the MultiSFTCatalogView
   MultiSFTCatalogView *catalogView = NULL;
   XLAL_CHECK_NULL( (catalogView = XLALGetMultiSFTCatalogView(catalog)) != NULL, XLAL_EFUNC );

   MultiLIGOTimeGPSVector *multiTimestamps = NULL;
   XLAL_CHECK_NULL( (multiTimestamps = XLALTimestampsFromMultiSFTCatalogView(catalogView)) != NULL, XLAL_EFUNC );

   XLALDestroyMultiSFTCatalogView(catalogView);

   return multiTimestamps;

}


/**
 * Create a list of timestamps from SFTs that might be a subset from those in an SFTCatalog, applying KS/Kuipers test if desired
 * \param [in] multiSFTVector Pointer to a MultiSFTVector
 * \param [in] params         Pointer to inputParamsStruct
 * \return Pointer to a list of GPS timestamps in a MultiLIGOTimeGPSVector
 */
MultiLIGOTimeGPSVector * getMultiTimeStampsFromSFTs(MultiSFTVector *multiSFTvector, inputParamsStruct *params)
{

   XLAL_CHECK_NULL( params != NULL, XLAL_EINVAL );

   MultiLIGOTimeGPSVector *multiTimestamps = NULL;
   XLAL_CHECK_NULL( (multiTimestamps = XLALCalloc(1, sizeof(*multiTimestamps))) != NULL, XLAL_ENOMEM );
   XLAL_CHECK_NULL( (multiTimestamps->data = XLALCalloc(multiSFTvector->length, sizeof(*multiTimestamps->data))) != NULL, XLAL_ENOMEM );
   multiTimestamps->length = multiSFTvector->length;

   if (params->markBadSFTs && !params->signalOnly) {
      REAL8 tfnormval = 2.0/(params->Tcoh*(params->avesqrtSh*params->avesqrtSh));

      for (INT4 ii=0; ii<(INT4)multiSFTvector->length; ii++) {
         REAL4Vector *tfdata = NULL;
         XLAL_CHECK_NULL( (tfdata = convertSFTdataToPowers(multiSFTvector->data[ii], params, tfnormval)) != NULL, XLAL_EFUNC );

         INT4Vector *removeTheseSFTs = NULL;
         XLAL_CHECK_NULL( (removeTheseSFTs = markBadSFTs(tfdata, params)) != NULL, XLAL_EFUNC );

         INT4 numberofsfts = 0, sftlength = (INT4)multiSFTvector->data[ii]->data[0].data->length;
         for (INT4 jj=0; jj<(INT4)removeTheseSFTs->length; jj++) if (removeTheseSFTs->data[jj]==0 && tfdata->data[jj*sftlength]!=0.0) numberofsfts++;

         XLAL_CHECK_NULL( (multiTimestamps->data[ii] = XLALCreateTimestampVector(numberofsfts)) != NULL, XLAL_EFUNC );

         INT4 kk = 0;
         for (INT4 jj=0; jj<(INT4)removeTheseSFTs->length; jj++) {
            if (removeTheseSFTs->data[jj]==0 && tfdata->data[jj*sftlength]!=0.0) {
               XLALGPSSetREAL8(&(multiTimestamps->data[ii]->data[kk]), params->searchstarttime+jj*(params->Tcoh-params->SFToverlap));
               XLAL_CHECK_NULL( xlalErrno == 0, XLAL_EFUNC );
               kk++;
            }
         }

         multiTimestamps->data[ii]->deltaT = params->Tcoh;

         XLALDestroyINT4Vector(removeTheseSFTs);
         XLALDestroyREAL4Vector(tfdata);
      }
      params->markBadSFTs = 0;
   } else {
      XLAL_CHECK_NULL( (multiTimestamps = XLALExtractMultiTimestampsFromSFTs(multiSFTvector)) != NULL, XLAL_EFUNC );
   }

   return multiTimestamps;

}


/**
 * Create a list of timestamps from a segment list
 * \param [in] file String for the filename
 * \param [in] params Pointer to inputParamsStruct
 * \return Pointer to a list of GPS timestamps in a MultiLIGOTimeGPSVector
 */
MultiLIGOTimeGPSVector * getMultiTimeStampsFromSegmentsFile(LALStringVector *filenames, inputParamsStruct *params)
{

   XLAL_CHECK_NULL( filenames != NULL && params != NULL, XLAL_EINVAL );

   MultiLIGOTimeGPSVector *multiTimestamps = NULL;
   XLAL_CHECK_NULL( (multiTimestamps = XLALCalloc(1, sizeof(MultiLIGOTimeGPSVector))) != NULL, XLAL_ENOMEM );
   multiTimestamps->length = params->numofIFOs;
   XLAL_CHECK_NULL( (multiTimestamps->data = XLALCalloc(multiTimestamps->length, sizeof(multiTimestamps->data[0]))) != NULL, XLAL_ENOMEM );

   for (INT4 ii=0; ii<params->numofIFOs; ii++) {
      LIGOTimeGPSVector *timestamps = NULL;
      XLAL_CHECK_NULL( (timestamps = XLALTimestampsFromSegmentFile(filenames->data[ii], params->Tcoh, params->SFToverlap, 0, 1)) != NULL, XLAL_EFUNC );
      multiTimestamps->data[ii] = timestamps;
   }

   return multiTimestamps;

}



/**
 * Slide the time-frequency data to account for detector motion
 * \param [out] output    Pointer to REAL4Vector of SFT powers that have been corrected
 * \param [in]  params    Pointer to inputParamsStruct
 * \param [in]  tfdata    Pointer to REAL4Vector of SFT powers
 * \param [in]  binshifts Pointer to INT4Vector of bin shift values
 * \return Status value
 */
INT4 slideTFdata(REAL4Vector *output, inputParamsStruct *params, REAL4Vector *tfdata, INT4Vector *binshifts)
{

   XLAL_CHECK( output != NULL && params != NULL && tfdata != NULL && binshifts != NULL, XLAL_EINVAL );

   INT4 ii;
   INT4 numffts = (INT4)floor(params->Tobs/(params->Tcoh-params->SFToverlap)-1);
   INT4 numfbins = (INT4)(round(params->fspan*params->Tcoh+2.0*params->dfmax*params->Tcoh)+12+1);

   for (ii=0; ii<numffts; ii++) {
      XLAL_CHECK( binshifts->data[ii]<=params->maxbinshift, XLAL_EFAILED, "SFT slide value %d is greater than maximum value predicted (%d)", binshifts->data[ii], params->maxbinshift );
      //for (jj=0; jj<numfbins; jj++) output->data[ii*numfbins + jj] = tfdata->data[ii*(numfbins+2*params->maxbinshift) + jj + params->maxbinshift + binshifts->data[ii]];
      memcpy(&(output->data[ii*numfbins]), &(tfdata->data[ii*(numfbins+2*params->maxbinshift) + params->maxbinshift + binshifts->data[ii]]), sizeof(REAL4)*numfbins);
   }

   return XLAL_SUCCESS;

} /* slideTFdata() */




/**
 * Determine the running mean of each SFT
 * \param [out] output   Pointer to REAL4Vector of running mean values of each SFT
 * \param [in]  tfdata   Pointer to REAL4Vector of SFT powers
 * \param [in]  numffts  Number of SFTs in the observation time
 * \param [in]  numfbins Number of frequency bins
 * \param [in]  blksize  Number of bins in the running median
 * \return Status value
 */
INT4 tfRngMeans(REAL4Vector *output, REAL4Vector *tfdata, INT4 numffts, INT4 numfbins, INT4 blksize)
{

   XLAL_CHECK( output != NULL && tfdata != NULL && numffts > 0  && numfbins > 0 && blksize > 0, XLAL_EINVAL );

   fprintf(LOG, "Assessing background... ");
   fprintf(stderr, "Assessing background... ");

   LALStatus XLAL_INIT_DECL(status);
   REAL8 bias;
   INT4 ii, jj;
   INT4 totalfbins = numfbins + blksize - 1;

   //Blocksize of running median
   LALRunningMedianPar block = {blksize};

   //Running median bias calculation
   if (blksize<1000) {
      bias = XLALRngMedBias(blksize);
      XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
   } else  bias = LAL_LN2;
   //REAL8 invbias = 1.0/(bias*1.0099993480677538);  //StackSlide normalization for 101 bins
   REAL8 invbias = 1.0/bias;

   //Allocate for a single SFT data and the medians out of each SFT
   REAL4Vector *inpsd = NULL, *mediansout = NULL;
   XLAL_CHECK( (inpsd = XLALCreateREAL4Vector(totalfbins)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (mediansout = XLALCreateREAL4Vector(numfbins)) != NULL, XLAL_EFUNC );

   //FILE *SFTMEANVALS = fopen("./outputtemp/sftmeanvals.dat","w");

   //Now do the running median
   for (ii=0; ii<numffts; ii++) {
      //If the SFT values were not zero, then compute the running median
      if (tfdata->data[ii*totalfbins]!=0.0) {
         //Determine running median value, convert to mean value
         memcpy(inpsd->data, &(tfdata->data[ii*inpsd->length]), sizeof(REAL4)*inpsd->length);

         //fprintf(SFTMEANVALS, "%g\n", calcMean(inpsd));

         //calculate running median
         LALSRunningMedian2(&status, mediansout, inpsd, block);
         XLAL_CHECK( status.statusCode == 0, XLAL_EFUNC );

         //Now make the output medians into means by multiplying by 1/bias
         for (jj=0; jj<(INT4)mediansout->length; jj++) output->data[ii*numfbins + jj] = (REAL4)(mediansout->data[jj]*invbias);
      } else {
         //Otherwise, set means to zero
         for (jj=0; jj<(INT4)mediansout->length; jj++) output->data[ii*numfbins + jj] = 0.0;
         //fprintf(SFTMEANVALS, "%g\n", 0.0);
      }
   } /* for ii < numffts */

   //fclose(SFTMEANVALS);

   fprintf(LOG, "done\n");
   fprintf(stderr, "done\n");

   fprintf(stderr,"Mean of running means = %g\n", calcMean(output));

   //Destroy stuff
   XLALDestroyREAL4Vector(inpsd);
   XLALDestroyREAL4Vector(mediansout);

   return XLAL_SUCCESS;

} /* tfRngMeans() */


/* Critical values of KS test (from Bickel and Doksum). Does not apply directly (mean determined from distribution)
alpha=0.01
n       10      20      30      40      50      60      80      n>80
        .489    .352    .290    .252    .226    .207    .179    1.628/(sqrt(n)+0.12+0.11/sqrt(n))

alpha=0.05
n       10      20      30      40      50      60      80      n>80
        .409    .294    .242    .210    .188    .172    .150    1.358/(sqrt(n)+0.12+0.11/sqrt(n))

alpha=0.1 (E.G derived using root finding)
n       10      20      30      40      50      60      80      n>80
        .369    .265    .218    .189    .170    .155    .135    1.224/(sqrt(n)+0.12+0.11/sqrt(n))

alpha=0.2 (E.G derived using root finding)
n       10      20      30      40      50      60      80      n>80
                                                                1.073/(sqrt(n)+0.12+0.11/sqrt(n))

alpha=0.5 (E.G derived using root finding)
n       10      20      30      40      50      60      80      n>80
        .249    .179    .147    .128    .115    .105    .091    0.828/(sqrt(n)+0.12+0.11/sqrt(n))

alpha=0.9 (E.G derived using root finding)
n                                                               n>80
                                                                0.571/(sqrt(n)+0.12+0.11/sqrt(n))

Critical values of Kuiper's test using root finding by E.G.
alpha=0.05
n                                                               n>80
                                                                1.747/(sqrt(n)+0.155+0.24/sqrt(n))

alpha=0.1
n                                                               n>80
                                                                1.620/(sqrt(n)+0.155+0.24/sqrt(n))

alpha=0.2
n                                                               n>80
                                                                1.473/(sqrt(n)+0.155+0.24/sqrt(n))
*/
/**
 * Mark the non-Gaussian SFTs using K-S and Kuiper's tests
 * \param [in] tfdata Pointer to REAL4Vector of SFT powers
 * \param [in] params Pointer to inputParamsStruct
 * \return Pointer to an INT4Vector with marked SFTs to be removed with 1 and keep SFT with 0
 */
INT4Vector * markBadSFTs(REAL4Vector *tfdata, inputParamsStruct *params)
{

   XLAL_CHECK_NULL( tfdata != NULL && params != NULL, XLAL_EINVAL );

   fprintf(stderr, "Marking bad SFTs... ");

   INT4 ii;

   INT4 numffts = (INT4)floor(params->Tobs/(params->Tcoh-params->SFToverlap)-1);    //Number of FFTs
   INT4 numfbins = (INT4)(round(params->fspan*params->Tcoh+2.0*params->dfmax*params->Tcoh)+12+1)+2*params->maxbinshift+params->blksize-1;     //Number of frequency bins

   //Allocate output data vector and a single SFT data vector
   INT4Vector *output = NULL;
   XLAL_CHECK_NULL( (output = XLALCreateINT4Vector(numffts)) != NULL, XLAL_EFUNC );
   memset(output->data, 0, sizeof(INT4)*output->length);
   REAL4Vector *tempvect = NULL;
   XLAL_CHECK_NULL( (tempvect = XLALCreateREAL4Vector(numfbins)) != NULL, XLAL_EFUNC );

   //Do the KS and Kuiper test on each SFT
   REAL8 ksthreshold = 1.358/(sqrt(numfbins)+0.12+0.11/sqrt(numfbins));
   //REAL8 ksthreshold = 1.224/(sqrt(numfbins)+0.12+0.11/sqrt(numfbins));  //This is a tighter restriction
   //REAL8 ksthreshold = 1.073/(sqrt(numfbins)+0.12+0.11/sqrt(numfbins));  //This is an even tighter restriction
   REAL8 kuiperthreshold = 1.747/(sqrt(numfbins)+0.155+0.24/sqrt(numfbins));
   //REAL8 kuiperthreshold = 1.620/(sqrt(numfbins)+0.155+0.24/sqrt(numfbins));  //This is a tighter restriction
   //REAL8 kuiperthreshold = 1.473/(sqrt(numfbins)+0.155+0.24/sqrt(numfbins));  //This is an even tighter restriction
   INT4 badsfts = 0, totalsfts = 0;
   REAL8 kstest = 0.0, kuipertest = 0.0;
   for (ii=0; ii<numffts; ii++) {
      if (tfdata->data[ii*numfbins]!=0.0) {
         totalsfts++;
         memcpy(tempvect->data, &(tfdata->data[ii*numfbins]), sizeof(REAL4)*tempvect->length);
         XLAL_CHECK_NULL( ks_test_exp(&kstest, tempvect) == XLAL_SUCCESS, XLAL_EFUNC );
         XLAL_CHECK_NULL( kuipers_test_exp(&kuipertest, tempvect) == XLAL_SUCCESS, XLAL_EFUNC );
         if (kstest>ksthreshold || kuipertest>kuiperthreshold) {
            output->data[ii] = 1;
            badsfts++;
         }
      }
   }

   fprintf(stderr, "done.\n");

   fprintf(stderr, "Fraction excluded in K-S and Kuiper's tests = %f\n", (REAL4)badsfts/(REAL4)totalsfts);

   //Destroy stuff
   XLALDestroyREAL4Vector(tempvect);

   return output;

}


/**
 * Remove the marked SFTs as bad by setting values to 0
 * \param [in,out] tfdata  Pointer to REAL4Vector of SFT powers
 * \param [in]     badsfts Poienter to INT4Vector of bad SFTs
 */
void removeBadSFTs(REAL4Vector *tfdata, INT4Vector *badsfts)
{

   XLAL_CHECK_VOID( tfdata != NULL && badsfts != NULL, XLAL_EINVAL );

   fprintf(stderr, "Removing bad SFTs... ");

   INT4 numfbins_tfdata = tfdata->length/badsfts->length;

   for (INT4 ii=0; ii<(INT4)badsfts->length; ii++) if (badsfts->data[ii]==1) memset(&(tfdata->data[ii*numfbins_tfdata]), 0, sizeof(REAL4)*numfbins_tfdata);

   fprintf(stderr, "done.\n");

}


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
 * \param [in] params Pointer to inputParamsStruct
 * \return Pointer to INT4Vector of bin number of the lines
 */
INT4Vector * detectLines_simple(REAL4Vector *TFdata, ffdataStruct *ffdata, inputParamsStruct *params)
{

   XLAL_CHECK_NULL( TFdata != NULL && ffdata != NULL && params != NULL, XLAL_EINVAL );

   LALStatus XLAL_INIT_DECL(status);
   INT4 blksize = 11, ii, jj;
   INT4 numlines = 0;
   INT4Vector *lines = NULL;
   INT4 totalnumfbins = ffdata->numfbins+(params->blksize-1)+2*params->maxbinshift;

   //Blocksize of running median
   LALRunningMedianPar block = {blksize};

   //Compute weights
   REAL4Vector *sftdata = NULL, *weights = NULL;
   XLAL_CHECK_NULL( (sftdata = XLALCreateREAL4Vector(totalnumfbins)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (weights = XLALCreateREAL4Vector(ffdata->numffts)) != NULL, XLAL_EFUNC );
   memset(weights->data, 0, ffdata->numffts*sizeof(REAL4));
   REAL4 sumweights = 0.0;
   for (ii=0; ii<ffdata->numffts; ii++) {
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

   for (ii=0; ii<(INT4)testRMSvals->length; ii++) {
      for (jj=0; jj<ffdata->numffts; jj++) testTSofPowers->data[jj] = TFdata->data[jj*testRMSvals->length + ii]*weights->data[jj]*invsumweights;
      XLAL_CHECK_NULL( calcRms(&(testRMSvals->data[ii]), testTSofPowers) == XLAL_SUCCESS, XLAL_EFUNC ); //This approaches calcMean(TSofPowers) for stationary noise
   }

   //Running median of RMS values
   LALSRunningMedian2(&status, testRngMedian, testRMSvals, block);
   XLAL_CHECK_NULL( status.statusCode == 0, XLAL_EFUNC );

   //Determine which bins are above the threshold and store the bin number of the line
   for (ii=0; ii<(INT4)testRngMedian->length; ii++) {
      REAL4 normrmsval = testRMSvals->data[ii+(blksize-1)/2]/testRngMedian->data[ii];
      if ( (ii+(blksize-1)/2) > ((params->blksize-1)/2) && normrmsval > params->lineDetection) {
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

   return lines;

}

/**
 * Track the lines for the sky position
 * \param [in] lines     Pointer to INT4Vector with list of lines
 * \param [in] binshifts Pointer to INT4Vector of SFT bin shifts
 * \param [in] params    Pointer to inputParamsStruct
 * \return Pointer to REAL4VectorSequence containing lower/mid/upper frequency for each line
 */
REAL4VectorSequence * trackLines(INT4Vector *lines, INT4Vector *binshifts, inputParamsStruct *params)
{

   XLAL_CHECK_NULL( lines != NULL && binshifts != NULL && params != NULL, XLAL_EINVAL );

   REAL4VectorSequence *output = NULL;
   XLAL_CHECK_NULL( (output = XLALCreateREAL4VectorSequence(lines->length, 3)) != NULL, XLAL_EFUNC );

   REAL4 df = 1.0/params->Tcoh;
   REAL4 minfbin = (REAL4)(round(params->fmin*params->Tcoh - params->dfmax*params->Tcoh - 6.0 - 0.5*(params->blksize-1) - (REAL8)(params->maxbinshift))/params->Tcoh);

   INT4 maxshiftindex = 0, minshiftindex = 0;
   min_max_index_INT4Vector(binshifts, &minshiftindex, &maxshiftindex);
   INT4 maxshift = binshifts->data[maxshiftindex], minshift = binshifts->data[minshiftindex];

   for (INT4 ii=0; ii<(INT4)lines->length; ii++) {
      output->data[ii*3] = lines->data[ii]*df + minfbin;
      //output->data[ii*3 + 1] = (lines->data[ii] + minshift)*df + minfbin;
      //output->data[ii*3 + 2] = (lines->data[ii] + maxshift)*df + minfbin;
      output->data[ii*3 + 1] = (lines->data[ii] + (minshift-1))*df + minfbin;  //Add one extra bin for buffer
      output->data[ii*3 + 2] = (lines->data[ii] + (maxshift+1))*df + minfbin;  //Add one extra bin for buffer
   }

   return output;

}


/**
 * Test algorithm to clean lines. NOT FULLY TESTED!
 * \param [in,out] TFdata     Pointer to time-frequency data in a REAL4Vector to be cleaned
 * \param [in]     background Pointer to REAL4Vector of running mean data
 * \param [in]     lines      Pointer to INT4Vector of lines
 * \param [in]     params     Pointer to inputParamsStruct
 * \return Status value
 */
INT4 cleanLines(REAL4Vector *TFdata, REAL4Vector *background, INT4Vector *lines, inputParamsStruct *params)
{

   if (lines==NULL) return XLAL_SUCCESS;

   INT4 numffts = (INT4)floor(params->Tobs/(params->Tcoh-params->SFToverlap)-1);    //Number of FFTs
   INT4 numfbins = (INT4)TFdata->length/numffts;     //Number of frequency bins

   REAL8 prevnoiseval = 0.0;
   for (INT4 ii=0; ii<(INT4)numffts; ii++) {
      if (TFdata->data[ii*numfbins]!=0.0) {
         for (INT4 jj=0; jj<(INT4)lines->length; jj++) {
            if (lines->data[jj]-(params->blksize-1)/2>=0) {
               REAL8 noiseval = expRandNum(background->data[ii*(numfbins-(params->blksize-1)) + lines->data[jj]], params->rng);
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
 * Determine if the SFTs are existing based on whether the first data point of power in each SFT is 0
 * \param [in] tfdata   Pointer to REAL4Vector of time-frequency data
 * \param [in] params   Pointer to inputParamsStruct
 * \param [in] numfbins Number of frequency bins in each SFT
 * \param [in] numffts  Number of SFTs from the observation time
 * \return Pointer to INT4Vector containing 0 for non-present SFT or 1 for present SFT
 */
INT4Vector * existingSFTs(REAL4Vector *tfdata, inputParamsStruct *params, INT4 numfbins, INT4 numffts)
{

   XLAL_CHECK_NULL( tfdata != NULL && params != NULL && numfbins > 0 && numffts > 0, XLAL_EINVAL );

   INT4 ii;

   INT4Vector *sftexist = NULL;
   XLAL_CHECK_NULL( (sftexist = XLALCreateINT4Vector(numffts)) != NULL, XLAL_EFUNC );

   for (ii=0; ii<numffts; ii++) {
      if (tfdata->data[ii*(numfbins+2*params->maxbinshift+params->blksize-1)] == 0.0) sftexist->data[ii] = 0;
      else {
         sftexist->data[ii] = 1;
      }
   }

   return sftexist;

} /* existingSFTs() */


/**
 * Subtract running mean values from the SFT data, modifying input time-frequency data
 * \param [in,out] tfdata   Pointer to REAL4Vector time-frequency data (modified by this function!)
 * \param [in]     rngMeans Pointer to REAL4Vector of running mean values
 * \param [in]     numffts  Number of SFTs from observation time
 * \param [in]     numfbins Number of SFT frequency bins
 * \return Status value
 */
INT4 tfMeanSubtract(REAL4Vector *tfdata, REAL4Vector *rngMeans, INT4 numffts, INT4 numfbins)
{

   XLAL_CHECK( tfdata != NULL && rngMeans != NULL && numffts > 0 && numfbins > 0, XLAL_EINVAL );

   for (INT4 ii=0; ii<numffts; ii++) if (rngMeans->data[ii*numfbins]!=0.0) for (INT4 jj=0; jj<numfbins; jj++) tfdata->data[ii*numfbins+jj] -= rngMeans->data[ii*numfbins+jj];

   return XLAL_SUCCESS;

} /* tfMeanSubtract() */


/**
 * Weight the SFTs based on antenna pattern and noise variance (Equation 11, assuming the input time-frequency data is already mean subtracted)
 * \param [out] output                    Pointer to REAL4Vector of mean subtracted, noise and antenna pattern weighted SFTs
 * \param [in]  tfdata                    Pointer to REAL4Vector of mean subtracted SFTs
 * \param [in]  rngMeans                  Pointer to REAL4Vector of running mean values
 * \param [in]  antPatternWeights         Pointer to REAL4Vector of antenna pattern weights
 * \param [in]  indexValuesOfExistingSFTs Pointer to INT4Vector of the index values of the existing SFTs
 * \param [in]  params                    Pointer to inputParamsStruct
 * \return Status value
 */
INT4 tfWeight(REAL4Vector *output, REAL4Vector *tfdata, REAL4Vector *rngMeans, REAL4Vector *antPatternWeights, INT4Vector *indexValuesOfExistingSFTs, inputParamsStruct *params)
{

   XLAL_CHECK( output != NULL && tfdata != NULL && rngMeans != NULL && antPatternWeights != NULL && indexValuesOfExistingSFTs != NULL && params != NULL, XLAL_EINVAL );

   INT4 ii, jj;

   INT4 numffts = (INT4)antPatternWeights->length;
   INT4 numfbins = (INT4)(tfdata->length)/numffts;

   //Initially set output to zero
   memset(output->data, 0, sizeof(REAL4)*output->length);

   REAL4Vector *antweightssq = NULL, *rngMeanssq = NULL;
   XLAL_CHECK( (antweightssq = XLALCreateREAL4Vector(numffts)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (rngMeanssq = XLALCreateREAL4Vector(numffts)) != NULL, XLAL_EFUNC );

   //If user specified that SFTs contain signal only, then the values of the backgrnd vector will be zeros.
   //We must set them equal to 1.0 here, then return them to zeros at the end
   if (params->signalOnly) {
      for (ii=0; ii<(INT4)indexValuesOfExistingSFTs->length; ii++) {
         for (jj=0; jj<numfbins; jj++) rngMeans->data[numfbins*indexValuesOfExistingSFTs->data[ii] + jj] = params->avesqrtSh;
      }
   }

   //User specifies whether to use SSE to do the multiplication or not
   if (params->useSSE) XLAL_CHECK( sseSSVectorMultiply(antweightssq, antPatternWeights, antPatternWeights) == XLAL_SUCCESS, XLAL_EFUNC );
   else if (params->useAVX) XLAL_CHECK( avxSSVectorMultiply(antweightssq, antPatternWeights, antPatternWeights) == XLAL_SUCCESS, XLAL_EFUNC );
   else {
      //for (ii=0; ii<numffts; ii++) antweightssq->data[ii] = antPatternWeights->data[ii]*antPatternWeights->data[ii];
      antweightssq = XLALSSVectorMultiply(antweightssq, antPatternWeights, antPatternWeights);
      XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
   }

   //Loop through the SFT frequency bins and weight the data
   for (ii=0; ii<numfbins; ii++) {

      fastSSVectorMultiply_with_stride_and_offset(rngMeanssq, rngMeans, rngMeans, numfbins, numfbins, ii, ii);

      //If noiseWeightOff is given, then set all the noise weights to be 1.0
      if (params->noiseWeightOff) for (jj=0; jj<(INT4)rngMeanssq->length; jj++) if (rngMeanssq->data[jj]!=0.0) rngMeanssq->data[jj] = 1.0;

      //Get sum of antenna pattern weight/variances for each frequency bin as a function of time (only for existant SFTs)
      REAL8 sumofweights = determineSumOfWeights(antweightssq, rngMeanssq);
      REAL8 invsumofweights = 1.0/sumofweights;

      //Now do noise weighting, antenna pattern weighting
      for (jj=0; jj<(INT4)indexValuesOfExistingSFTs->length; jj++) {
         output->data[indexValuesOfExistingSFTs->data[jj]*numfbins+ii] = (REAL4)(invsumofweights*antPatternWeights->data[indexValuesOfExistingSFTs->data[jj]]*tfdata->data[indexValuesOfExistingSFTs->data[jj]*numfbins+ii]/rngMeanssq->data[indexValuesOfExistingSFTs->data[jj]]);
      } /* for jj < indexValuesOfExisitingSFTs->length */
   } /* for ii < numfbins */

   //Remember to reset the backgrnd vector to zero
   if (params->signalOnly) memset(rngMeans->data, 0, sizeof(REAL4)*rngMeans->length);

   //Destroy stuff
   XLALDestroyREAL4Vector(antweightssq);
   XLALDestroyREAL4Vector(rngMeanssq);

   //fprintf(stderr,"TF after weighting, mean subtraction = %g\n",calcMean(output));

   return XLAL_SUCCESS;

} /* tfWeight() */


/**
 * Determine the sum of the weights
 * \param [in] antweightssq Antenna pattern weights squared
 * \param [in] rngMeanssq   Running mean values squared
 * \return Sum of the weights
 */
REAL8 determineSumOfWeights(REAL4Vector *antweightssq, REAL4Vector *rngMeanssq)
{

   XLAL_CHECK_REAL8( antweightssq != NULL && rngMeanssq != NULL, XLAL_EINVAL );

   INT4 ii;
   REAL8 sumofweights = 0.0;
   for (ii=0; ii<(INT4)antweightssq->length; ii++) if (rngMeanssq->data[ii] != 0.0) sumofweights += antweightssq->data[ii]/rngMeanssq->data[ii];

   return sumofweights;

} /* determineSumOfWeights */


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

   INT4 ii, jj;
   REAL8 winFactor = 8.0/3.0;

   //Do the second FFT
   REAL4Vector *x = NULL, *psd = NULL;
   REAL4Window *win = NULL;
   XLAL_CHECK( (x = XLALCreateREAL4Vector(output->numffts)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (psd = XLALCreateREAL4Vector((UINT4)floor(x->length*0.5)+1)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (win = XLALCreateHannREAL4Window(x->length)) != NULL, XLAL_EFUNC );

   for (ii=0; ii<output->numfbins; ii++) {

      //Next, loop over times and pick the right frequency bin for each FFT and window
      //for (jj=0; jj<(INT4)x->length; jj++) x->data[jj] = (tfdata->data[ii + jj*numfbins]*win->data->data[jj]);
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
      for (jj=0; jj<(INT4)psd->length; jj++) output->ffdata->data[psd->length*ii + jj] = (REAL4)(psd->data[jj]*winFactor*output->ffnormalization);

   } /* for ii < numfbins */

   //Destroy stuff
   XLALDestroyREAL4Vector(x);
   XLALDestroyREAL4Vector(psd);
   XLALDestroyREAL4Window(win);

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

   INT4 ii;
   REAL4Vector *aveNoiseInTime = NULL, *rngMeansOverBand = NULL;
   XLAL_CHECK_REAL4( (aveNoiseInTime = XLALCreateREAL4Vector(numffts)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_REAL4( (rngMeansOverBand = XLALCreateREAL4Vector((UINT4)(binmax-binmin))) != NULL, XLAL_EFUNC );

   for (ii=0; ii<(INT4)aveNoiseInTime->length; ii++) {
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

   INT4 ii;
   REAL4Vector *aveNoiseInTime = NULL, *rngMeansOverBand = NULL;
   XLAL_CHECK_REAL4( (aveNoiseInTime = XLALCreateREAL4Vector(numffts)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_REAL4( (rngMeansOverBand = XLALCreateREAL4Vector((UINT4)(binmax-binmin))) != NULL, XLAL_EFUNC );

   for (ii=0; ii<(INT4)aveNoiseInTime->length; ii++) {
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
 * \param [out]    aveNoise      Pointer to REAL4Vector of the expected 2nd FFT powers
 * \param [in]     params        Pointer to inputParamsStruct
 * \param [in]     sftexist      Pointer to INT4Vector of SFTs existing or not
 * \param [in]     backgrnd      Pointer to REAL4Vector of running means
 * \param [in]     antweights    Pointer to REAL4Vector of antenna pattern weights
 * \param [in]     plan          Pointer to REAL4FFTPlan
 * \param [in,out] normalization Pointer to REAL8 value of the normalization for the 2nd FFT
 * \return Status value
 */
INT4 ffPlaneNoise(REAL4Vector *aveNoise, inputParamsStruct *params, INT4Vector *sftexist, REAL4Vector *backgrnd, REAL4Vector *antweights, REAL4FFTPlan *plan, REAL8 *normalization)
{

   XLAL_CHECK( aveNoise != NULL && params != NULL && sftexist != NULL && backgrnd != NULL && antweights != NULL && plan != NULL && normalization != NULL, XLAL_EINVAL );

   INT4 ii, jj, numfbins, numffts, numfprbins;
   REAL8 invsumofweights = 0.0;
   REAL8 sumofweights = 0.0;

   numfbins = (INT4)(round(params->fspan*params->Tcoh+2.0*params->dfmax*params->Tcoh)+12+1);    //Number of frequency bins
   numffts = (INT4)antweights->length;                      //Number of FFTs
   numfprbins = (INT4)floor(numffts*0.5)+1;                 //number of 2nd fft frequency bins

   //Set up for making the PSD
   memset(aveNoise->data, 0, sizeof(REAL4)*aveNoise->length);

   //If the user has said there is signal only and no noise in the SFTs, then the noise background of the FF plane will be filled with zeros
   if (!params->signalOnly) {
      //Window and psd allocation
      REAL4Window *win = NULL;
      REAL4Vector *psd = NULL;
      XLAL_CHECK( (win = XLALCreateHannREAL4Window(numffts)) != NULL, XLAL_EFUNC );  //Window function
      XLAL_CHECK( (psd = XLALCreateREAL4Vector(numfprbins)) != NULL, XLAL_EFUNC );   //Current PSD calculation
      REAL4 winFactor = 8.0/3.0;
      REAL8 dutyfactor = 0.0, dutyfactorincrement = 1.0/(REAL8)numffts;

      //Average each SFT across the frequency band, also compute normalization factor
      REAL4Vector *aveNoiseInTime = NULL, *rngMeansOverBand = NULL;
      XLAL_CHECK( (aveNoiseInTime = XLALCreateREAL4Vector(numffts)) != NULL, XLAL_EFUNC );
      XLAL_CHECK( (rngMeansOverBand = XLALCreateREAL4Vector(numfbins)) != NULL, XLAL_EFUNC );

      memset(aveNoiseInTime->data, 0, sizeof(REAL4)*numffts);
      for (ii=0; ii<(INT4)aveNoiseInTime->length; ii++) {
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
      REAL4Vector *x = NULL, *multiplicativeFactor = NULL;
      XLAL_CHECK( (x = XLALCreateREAL4Vector(aveNoiseInTime->length)) != NULL, XLAL_EFUNC );
      XLAL_CHECK( (multiplicativeFactor = XLALCreateREAL4Vector(aveNoiseInTime->length)) != NULL, XLAL_EFUNC );

      memset(multiplicativeFactor->data, 0, sizeof(REAL4)*multiplicativeFactor->length);
      REAL4 psdfactor = winFactor*params->Tobs;  //Only multiply by Tobs instead of Tobs/numffts^2 because of the weighting normalization

      for (ii=0; ii<(INT4)x->length; ii++) {
         if (aveNoiseInTime->data[ii] != 0.0) {
            if (!params->noiseWeightOff) multiplicativeFactor->data[ii] = win->data->data[ii]*antweights->data[ii]/aveNoiseInTime->data[ii]*invsumofweights;
            else multiplicativeFactor->data[ii] = win->data->data[ii]*antweights->data[ii]*invsumofweights;
         }
      }

      //New version. Computes expected background with correlation estimate from Hann windowed and overlapped (must be 50% or 0%) SFTs
      REAL8 correlationfactor = 0.0;
      for (ii=0; ii<(INT4)floor(win->data->length*(params->SFToverlap/params->Tcoh)-1); ii++) correlationfactor += win->data->data[ii]*win->data->data[ii + (INT4)((1.0-(params->SFToverlap/params->Tcoh))*win->data->length)];
      correlationfactor /= win->sumofsquares;
      REAL8 corrfactorsquared = correlationfactor*correlationfactor;
      REAL8 prevnoiseval = 0.0, noiseval = 0.0;
      for (ii=0; ii<4000; ii++) {
         memset(x->data, 0, sizeof(REAL4)*x->length);
         for (jj=0; jj<(INT4)x->length; jj++) {
            if (sftexist->data[jj] != 0) {
               //To create the correlations
               noiseval = expRandNum(aveNoiseInTime->data[jj], params->rng);
               XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
               if (jj>0 && sftexist->data[jj-1]!=0) {
                  noiseval *= (1.0-corrfactorsquared);
                  noiseval += corrfactorsquared*prevnoiseval;
               }
               x->data[jj] = (REAL4)(noiseval/aveNoiseInTime->data[jj]-1.0);
               prevnoiseval = noiseval;
            }
         } /* for jj < x->length */

         //Reverse the correlations, comment this out below!
         /* memset(x->data, 0, sizeof(REAL4)*x->length);
         for (jj=(INT4)x->length-1; jj>=0; jj--) {
            if (sftexist->data[jj] != 0) {
               noiseval = expRandNum(aveNoiseInTime->data[jj], input->rng);
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
         if (params->useSSE) XLAL_CHECK( sseSSVectorMultiply(x, x, multiplicativeFactor) == XLAL_SUCCESS, XLAL_EFUNC );
         else if (params->useAVX) XLAL_CHECK( avxSSVectorMultiply(x, x, multiplicativeFactor) == XLAL_SUCCESS, XLAL_EFUNC );
         else for (jj=0; jj<(INT4)x->length; jj++) x->data[jj] *= multiplicativeFactor->data[jj];

         //Do the FFT
         XLAL_CHECK( XLALREAL4PowerSpectrum(psd, x, plan) == XLAL_SUCCESS, XLAL_EFUNC );

         //Sum into the bins
         if (params->useSSE) XLAL_CHECK( sseSSVectorSum(aveNoise, aveNoise, psd) == XLAL_SUCCESS, XLAL_EFUNC );
         else if (params->useAVX) XLAL_CHECK( avxSSVectorSum(aveNoise, aveNoise, psd) == XLAL_SUCCESS, XLAL_EFUNC );
         else for (jj=0; jj<(INT4)aveNoise->length; jj++) aveNoise->data[jj] += psd->data[jj];
      } /* for ii < 4000 */

      //Average and rescale
      //REAL4 averageRescaleFactor = 2.5e-4*psdfactor*(1.0+2.0*corrfactorsquared);
      REAL4 averageRescaleFactor = 2.5e-4*psdfactor; //*(1.0+2.0*corrfactorsquared);
      if (params->useSSE) XLAL_CHECK( sseScaleREAL4Vector(aveNoise, aveNoise, averageRescaleFactor) == XLAL_SUCCESS, XLAL_EFUNC );
      else if (params->useAVX) XLAL_CHECK( avxScaleREAL4Vector(aveNoise, aveNoise, averageRescaleFactor) == XLAL_SUCCESS, XLAL_EFUNC );
      else for (ii=0; ii<(INT4)aveNoise->length; ii++) aveNoise->data[ii] *= averageRescaleFactor;

      //Fix 0th and end bins (0 only for odd x->length, 0 and end for even x->length)
      if (GSL_IS_EVEN(x->length)==1) {
         aveNoise->data[0] *= 2.0;
         aveNoise->data[aveNoise->length-1] *= 2.0;
      } else {
         aveNoise->data[0] *= 2.0;
      }

      //Compute normalization
      *(normalization) = 1.0/(calcMean(aveNoise));
      for (ii=0; ii<(INT4)aveNoise->length; ii++) aveNoise->data[ii] *= *(normalization);

      //Extra factor for normalization to 1.0 (empirically determined) (Old method)
      //*normalization /= 1.08;

      //Extra factor for normalization to 1.0
      //duty factor because there are zero-valued SFTs (fourth root term)
      //square root term because the uncertainty on the background improves with the size of the running median block
      //correlation term because neighboring sfts are correlated
      *normalization /= (1.0+pow(dutyfactor,0.25)/sqrt((REAL8)params->blksize))*sqrt(1.0+2.0*corrfactorsquared);

      //Destroy stuff
      XLALDestroyREAL4Vector(x);
      XLALDestroyREAL4Vector(psd);
      XLALDestroyREAL4Window(win);
      XLALDestroyREAL4Vector(aveNoiseInTime);
      XLALDestroyREAL4Vector(rngMeansOverBand);
      XLALDestroyREAL4Vector(multiplicativeFactor);
   } else {
     *(normalization) = 1.0;
   }

   return XLAL_SUCCESS;

} /* ffPlaneNoise() */


INT4 readTwoSpectInputParams(inputParamsStruct *params, UserInput_t *uvar, int argc, char *argv[])
{
   XLAL_CHECK ( uvar != NULL, XLAL_EINVAL, "Invalid NULL input 'uvar'\n");
   XLAL_CHECK ( argv != NULL, XLAL_EINVAL, "Invalid NULL input 'argv'\n");

   INT4 ii;

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
   uvar->injRandSeed = 0;
   uvar->ULsolver = 0;
   uvar->dopplerMultiplier = 1.0;
   uvar->ihsfar = 1.0;
   uvar->tmplfar = 1.0;
   uvar->ihsfomfar = 1.0;
   uvar->ihsfom = 0.0;
   uvar->keepOnlyTopNumIHS = -1;
   uvar->lineDetection = -1.0;

   XLALregBOOLUserStruct(help,                         'h', UVAR_HELP    ,  "Print this help/usage message");
   XLALregSTRINGUserStruct(outdirectory,                0 , UVAR_REQUIRED,  "Output directory");
   XLALregLISTUserStruct(IFO,                           0 , UVAR_REQUIRED,  "CSV list of detectors, eg. \"H1,H2,L1,G1, ...\" ");
   XLALregREALUserStruct(Tobs,                          0 , UVAR_REQUIRED,  "Total observation time (in seconds)");
   XLALregREALUserStruct(Tcoh,                          0 , UVAR_REQUIRED,  "SFT coherence time (in seconds)");
   XLALregREALUserStruct(SFToverlap,                    0 , UVAR_REQUIRED,  "SFT overlap (in seconds), usually Tcoh/2");
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
   XLALregINTUserStruct(harmonicNumToSearch,            0 , UVAR_OPTIONAL,  "Number of harmonics of the Pmin to Pmax range to search");
   XLALregINTUserStruct(periodHarmToCheck,              0 , UVAR_OPTIONAL,  "Number of harmonics/sub-harmonics of the IHS candidates to test");
   XLALregINTUserStruct(periodFracToCheck,              0 , UVAR_OPTIONAL,  "Number of fractional periods to check in the sense of [(1...N)+1]/[(1...N)+2]");
   XLALregBOOLUserStruct(templateSearch,                0 , UVAR_OPTIONAL,  "Flag for doing a pure template-based search on search region specified by (sky,f,fspan,P, Asini +- 3 AsiniSigma)");
   XLALregREALUserStruct(templateSearchP,               0 , UVAR_OPTIONAL,  "The template search period; templateSearch flag is required");
   XLALregREALUserStruct(templateSearchAsini,           0 , UVAR_OPTIONAL,  "The template search Asini; templateSearch flag is required");
   XLALregREALUserStruct(templateSearchAsiniSigma,      0 , UVAR_OPTIONAL,  "The template search uncertainty in Asini; templateSearch flag is required");
   XLALregREALUserStruct(assumeNScosi,                  0 , UVAR_OPTIONAL,  "Assume cosi orientation of the source (only used when specifying more than 1 detector)");
   XLALregREALUserStruct(assumeNSpsi,                   0 , UVAR_OPTIONAL,  "Assume psi polarization angle of the source (only used when specifying more than 1 detector)");
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
   XLALregBOOLUserStruct(useSSE,                        0 , UVAR_OPTIONAL,  "Use SSE functions (caution: user needs to have compiled for SSE or program fails)");
   XLALregBOOLUserStruct(useAVX,                        0 , UVAR_OPTIONAL,  "Use AVX functions (caution: user needs to have compiled for AVX or program fails)");
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
   XLALregBOOLUserStruct(BrentsMethod,                  0 , UVAR_DEVELOPER, "Use Brent's method for root finding");
   XLALregBOOLUserStruct(printSFTtimes,                 0 , UVAR_DEVELOPER, "Output a list <GPS sec> <GPS nanosec> of SFT start times of SFTs");
   XLALregBOOLUserStruct(printData,                     0 , UVAR_DEVELOPER, "Print to ASCII files the data values");
   XLALregSTRINGUserStruct(printSignalData,             0 , UVAR_DEVELOPER, "Print f0 and h0 per SFT of the signal, used only with injectionSources");
   XLALregSTRINGUserStruct(printMarginalizedSignalData, 0 , UVAR_DEVELOPER, "Print f0 and h0 per SFT of the signal, used only with injectionSources");
   XLALregINTUserStruct(randSeed,                       0 , UVAR_DEVELOPER, "Random seed value");
   XLALregBOOLUserStruct(chooseSeed,                    0 , UVAR_DEVELOPER, "The random seed value is chosen based on the input search parameters");

   //Read all the input from config file and command line (command line has priority)
   //Also checks required variables unless help is requested
   XLAL_CHECK( XLALUserVarReadAllInput(argc, argv) == XLAL_SUCCESS, XLAL_EFUNC );

   //Help and exit
   if (uvar->help) exit(0);

   //Set observation parameters
   params->Tcoh = uvar->Tcoh;
   params->Tobs = uvar->Tobs;
   params->SFToverlap = uvar->SFToverlap;
   params->searchstarttime = uvar->t0;
   params->fmin = uvar->fmin;
   params->fspan = uvar->fspan;

   //Set detector parameters
   params->numofIFOs = uvar->IFO->length;
   for (ii=0; ii<params->numofIFOs; ii++) {
      if (strcmp("H1", uvar->IFO->data[ii])==0) {
         params->detectors->sites[ii] = lalCachedDetectors[LAL_LHO_4K_DETECTOR]; //H1
      } else if (strcmp("L1",uvar->IFO->data[ii])==0) {
         params->detectors->sites[ii] = lalCachedDetectors[LAL_LLO_4K_DETECTOR]; //L1
      } else if (strcmp("V1", uvar->IFO->data[ii])==0) {
         params->detectors->sites[ii] = lalCachedDetectors[LAL_VIRGO_DETECTOR];  //V1
      } else if (strcmp("H2", uvar->IFO->data[ii])==0) {
         params->detectors->sites[ii] = lalCachedDetectors[LAL_LHO_2K_DETECTOR]; //H2
      } else if (strcmp("H2r", uvar->IFO->data[ii])==0) {
         LALDetector H2 = lalCachedDetectors[LAL_LHO_2K_DETECTOR]; //H2 rotated
         H2.frDetector.xArmAzimuthRadians -= 0.25*LAL_PI;
         H2.frDetector.yArmAzimuthRadians -= 0.25*LAL_PI;
         memset(&(H2.frDetector.name), 0, sizeof(CHAR)*LALNameLength);
         snprintf(H2.frDetector.name, LALNameLength, "%s", "LHO_2k_rotatedPiOver4");
         XLAL_CHECK( (XLALCreateDetector(&(params->detectors->sites[ii]), &(H2.frDetector), LALDETECTORTYPE_IFODIFF)) != NULL, XLAL_EFUNC );
      } else {
         XLAL_ERROR(XLAL_EINVAL, "Not using valid interferometer! Expected 'H1', 'H2', 'H2r' (rotated H2), 'L1', or 'V1' not %s.\n", uvar->IFO->data[ii]);
      }
   }
   params->detectors->length = params->numofIFOs;
   params->avesqrtSh = uvar->avesqrtSh;

   //Set analysis parameters
   params->Pmax = uvar->Pmax;
   params->Pmin = uvar->Pmin;
   params->dfmax = uvar->dfmax;
   params->dfmin = uvar->dfmin;
   XLAL_CHECK( params->Pmax >= params->Pmin, XLAL_EINVAL, "Pmax is smaller than Pmin\n" );
   XLAL_CHECK( params->dfmax >= params->dfmin, XLAL_EINVAL, "dfmax is smaller than dfmin\n" );
   if (params->Pmax > 0.2*params->Tobs) {
      params->Pmax = 0.2*params->Tobs;
      fprintf(stderr,"WARNING! Adjusting input maximum period to 1/5 the observation time!\n");
   }
   if (params->Pmin < 4.0*params->Tcoh) {
      params->Pmin = 4.0*params->Tcoh;
      fprintf(stderr,"WARNING! Adjusting input minimum period to 4 coherence times!\n");
   }
   if (params->dfmax > maxModDepth(params->Pmax, params->Tcoh)) {
      params->dfmax = floor(maxModDepth(params->Pmax, params->Tcoh)*(params->Tcoh))/(params->Tcoh);
      fprintf(stderr,"WARNING! Adjusting input maximum modulation depth due to maximum modulation depth allowed\n");
   }
   if (params->dfmin < 0.5/params->Tcoh) {
      params->dfmin = 0.5/params->Tcoh;
      fprintf(stderr,"WARNING! Adjusting input minimum modulation depth to 1/2 a frequency bin!\n");
   }
   params->dfmin = 0.5*round(2.0*params->dfmin*params->Tcoh)/params->Tcoh;
   params->dfmax = 0.5*round(2.0*params->dfmax*params->Tcoh)/params->Tcoh;

   //Blocksize settings for running mean
   params->blksize = uvar->blksize;
   if (params->blksize % 2 != 1) params->blksize += 1;

   //SFT input file/directory or we only have a segmentFile/timestampsFile leaving params->inputSFTs==NULL
   if (XLALUserVarWasSet(&uvar->timestampsFile) && XLALUserVarWasSet(&uvar->inputSFTs)) {
      XLAL_ERROR(XLAL_FAILURE, "Cannot specify timestampsFile and inputSFTs\n");
   } else if (XLALUserVarWasSet(&uvar->segmentFile) && XLALUserVarWasSet(&uvar->inputSFTs)) {
      XLAL_ERROR(XLAL_FAILURE, "Cannot specify segmentFile and inputSFTs\n");
   } else if (XLALUserVarWasSet(&uvar->segmentFile) && XLALUserVarWasSet(&uvar->timestampsFile)) {
      XLAL_ERROR(XLAL_FAILURE, "Cannot specify segmentFile and timestampsFile\n");
   } else if (XLALUserVarWasSet(&uvar->inputSFTs)) {
      XLAL_CHECK( (params->inputSFTs = XLALStringDuplicate(uvar->inputSFTs)) != NULL, XLAL_EFUNC );
   } else {
      params->inputSFTs = NULL;
   }

   //Check skyRegion and skyRegionFile options
   if ((XLALUserVarWasSet(&uvar->skyRegion) && XLALUserVarWasSet(&uvar->skyRegionFile)) || (!XLALUserVarWasSet(&uvar->skyRegion) && !XLALUserVarWasSet(&uvar->skyRegionFile))) XLAL_ERROR(XLAL_EINVAL, "Specify only one of skyRegion or skyRegionFile\n");

   //IHS range
   params->harmonicNumToSearch = uvar->harmonicNumToSearch;
   params->periodHarmToCheck = uvar->periodHarmToCheck;
   params->periodFracToCheck = uvar->periodFracToCheck;

   //Check required options if specifying templateSearch
   if (uvar->templateSearch || XLALUserVarWasSet(&uvar->templateSearchP) || XLALUserVarWasSet(&uvar->templateSearchAsini) || XLALUserVarWasSet(&uvar->templateSearchAsiniSigma)) {
      if (!(uvar->templateSearch && XLALUserVarWasSet(&uvar->templateSearchP) && XLALUserVarWasSet(&uvar->templateSearchAsini) && XLALUserVarWasSet(&uvar->templateSearchAsiniSigma))) XLAL_ERROR(XLAL_FAILURE, "Must specify all or none of templateSearch, templateSearchP, templateSearchAsini, and templateSearchAsiniSigma\n");
   }

   //Error if dfmax is smaller than templateTestDf or if dfmax is smaller than the templateSearch or bruteForceTemplateTest largest modulation depth
   if (uvar->templateTest) XLAL_CHECK( params->dfmax >= uvar->templateTestDf, XLAL_EINVAL, "templateTestDf is larger than dfmax\n" );
   if (uvar->templateSearch) XLAL_CHECK( params->dfmax >= LAL_TWOPI*(params->fmin+params->fspan)*(uvar->templateSearchAsini+3.0*uvar->templateSearchAsiniSigma)/uvar->templateSearchP, XLAL_EINVAL, "templateSearch parameters would make the largest modulation depth larger than dfmax\n");
   if (uvar->bruteForceTemplateTest) XLAL_CHECK( params->dfmax >= uvar->templateTestDf+2.0/params->Tcoh, XLAL_EINVAL, "templateTestDf+2/Tcoh is larger than dfmax\n" );

   //IHS folding factor
   params->ihsfactor = uvar->ihsfactor;

   //Statistic threshold settings
   params->ihsfar = uvar->ihsfar;
   params->templatefar = uvar->tmplfar;

   //IHS FOM threshold settings
   if (XLALUserVarWasSet(&uvar->ihsfomfar) && XLALUserVarWasSet(&uvar->ihsfom)) XLAL_ERROR(XLAL_FAILURE, "Only choose only one of ihsfomfar or ihsfom\n");
   params->ihsfomfar = uvar->ihsfomfar;
   params->ihsfom = uvar->ihsfom;

   //Check ranges of settings for FAR values
   if (params->ihsfar<0.0 || params->ihsfar>1.0) XLAL_ERROR(XLAL_EINVAL, "ihsfar must lie between 0 and 1 (inclusive)\n");
   if (params->ihsfomfar<0.0 || params->ihsfomfar>1.0) XLAL_ERROR(XLAL_EINVAL, "ihsfomfar must lie between 0 and 1 (inclusive)\n");
   if (params->templatefar<0.0 || params->templatefar>1.0) XLAL_ERROR(XLAL_EINVAL, "tmplfar must lie between 0 and 1 (inclusive)\n");

   params->log10templatefar = log10(params->templatefar);

   //Keep only top number of IHS candidates
   params->keepOnlyTopNumIHS = uvar->keepOnlyTopNumIHS;

   //Template length settings
   params->mintemplatelength = uvar->minTemplateLength;
   params->maxtemplatelength = uvar->maxTemplateLength;

   //Upper limit settings take span of search values unless specified
   if (XLALUserVarWasSet(&uvar->ULfmin)) params->ULfmin = uvar->ULfmin;            //Upper limit minimum frequency (Hz)
   else params->ULfmin = params->fmin;
   if (XLALUserVarWasSet(&uvar->ULfspan)) params->ULfspan = uvar->ULfspan;         //Upper limit frequency span (Hz)
   else params->ULfspan = params->fspan;
   if (XLALUserVarWasSet(&uvar->ULminimumDeltaf)) params->ULmindf = uvar->ULminimumDeltaf;     //Upper limit minimum modulation depth (Hz)
   else params->ULmindf = params->dfmin;
   if (XLALUserVarWasSet(&uvar->ULmaximumDeltaf)) params->ULmaxdf = uvar->ULmaximumDeltaf;     //Upper limit maximum modulation depth (Hz)
   else params->ULmaxdf = params->dfmax;

   //Print all UL values at each sky location
   if (uvar->allULvalsPerSkyLoc) params->printAllULvalues = 1;
   else params->printAllULvalues = 0;

   //markBadSFTs and lineDetection
   if (uvar->markBadSFTs) params->markBadSFTs = 1;
   else params->markBadSFTs = 0;
   params->lineDetection = uvar->lineDetection;

   //FFT plan flag
   params->FFTplanFlag = uvar->FFTplanFlag;

   //fastchisqinv flag
   if (uvar->fastchisqinv) params->fastchisqinv = 1;
   else params->fastchisqinv = 0;

   //SSE/AVX settings
   if (uvar->useSSE) params->useSSE = 1;
   else params->useSSE = 0;
   if (uvar->useAVX) params->useAVX = 1;
   else params->useAVX = 0;
   if (params->useSSE && params->useAVX) XLAL_ERROR(XLAL_FAILURE, "Can only choose one of useSSE or useAVX, but not both (note that useAVX implies useSSE)");

   //Follow up outside UL range
   if (uvar->followUpOutsideULrange) params->followUpOutsideULrange = 1;
   else params->followUpOutsideULrange = 0;

   //Developer options
   if (uvar->weightedIHS) params->weightedIHS = 1;
   else params->weightedIHS = 0;
   if (uvar->signalOnly) params->signalOnly = 1;
   else params->signalOnly = 0;
   if (uvar->templateTest && uvar->bruteForceTemplateTest) XLAL_ERROR(XLAL_FAILURE, "Specify one of templateTest or bruteForceTemplateTest\n");
   if ((uvar->templateTest || uvar->bruteForceTemplateTest) || XLALUserVarWasSet(&uvar->templateTestF) || XLALUserVarWasSet(&uvar->templateTestP) || XLALUserVarWasSet(&uvar->templateTestDf)) {
      if (!((uvar->templateTest || uvar->bruteForceTemplateTest) && XLALUserVarWasSet(&uvar->templateTestF) && XLALUserVarWasSet(&uvar->templateTestP) && XLALUserVarWasSet(&uvar->templateTestDf))) {
         XLAL_ERROR(XLAL_FAILURE, "Must specify all or none of templateTest/bruteForceTemplateTest, templateTestF, templateTestP, and templateTestDf\n");
      }
   }
   params->ULsolver = uvar->ULsolver;
   params->dopplerMultiplier = uvar->dopplerMultiplier;
   if (uvar->noNotchHarmonics) params->noNotchHarmonics = 1;
   else params->noNotchHarmonics = 0;
   if (uvar->calcRthreshold) params->calcRthreshold = 1;
   else params->calcRthreshold = 0;
   if (uvar->antennaOff) params->antennaOff = 1;
   else params->antennaOff = 0;
   if (uvar->noiseWeightOff) params->noiseWeightOff = 1;
   else params->noiseWeightOff = 0;
   if (uvar->BrentsMethod) params->rootFindingMethod = 1;
   else params->rootFindingMethod = 0;
   params->randSeed = uvar->randSeed;

   //Random seed value settings for random number generator
   //If the chooseSeed option was given, then:
   //seed = (IFO multiplier)*fabs(round(fmin + fspan + Pmin + Pmax + dfmin + dfmax + alpha + delta))
   if (uvar->chooseSeed) {
      UINT8 IFOmultiplier = 0;
      if (strcmp(params->detectors->sites[0].frDetector.prefix, "H1")==0) IFOmultiplier = 1;            //H1 gets multiplier 1
      else if (strcmp(params->detectors->sites[0].frDetector.prefix, "L1")==0) IFOmultiplier = 2;       //L1 gets multiplier 2
      else IFOmultiplier = 3;                                                                           //V1 gets multiplier 3
      params->randSeed = IFOmultiplier*(UINT8)fabs(round(params->fmin + params->fspan + params->Pmin + params->Pmax + params->dfmin + params->dfmax));
   }
   gsl_rng_set(params->rng, params->randSeed);     //Set the random number generator with the given seed

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
   for (INT4 ii=0; ii<(INT4)vector->length; ii++) fprintf(OUTPUT, "%g\n", vector->data[ii]);
   fclose(OUTPUT);

   return XLAL_SUCCESS;

}
