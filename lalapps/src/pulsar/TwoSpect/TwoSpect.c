/*
*  Copyright (C) 2010, 2011, 2012 Evan Goetz
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

#include <lal/Window.h>
#include <lal/LALMalloc.h>
#include <lal/SFTutils.h>
#include <lal/SFTfileIO.h>
#include <lal/DopplerScan.h>
#include <lal/VectorOps.h>

#include <gsl/gsl_math.h>

#include "IHS.h"
#include "candidates.h"
#include "antenna.h"
#include "templates.h"
#include "TwoSpect.h"
#include "statistics.h"
#include "upperlimits.h"
#include "vectormath.h"


//Global variables
FILE *LOG = NULL, *ULFILE = NULL, *NORMRMSOUT = NULL;
CHAR *earth_ephemeris = NULL, *sun_ephemeris = NULL, *sft_dir = NULL;
static const LALStatus empty_status;

//Main program
int main(int argc, char *argv[])
{
   
   INT4 ii, jj;               //counter variables
   LALStatus status = empty_status;  //LALStatus structure
   char s[1000], t[1000], u[1000];   //Path and file name to LOG, ULFILE, and NORMRMSOUT
   time_t programstarttime, programendtime;
   struct tm *ptm;
   
   time(&programstarttime);
   ptm = localtime(&programstarttime);
   
   //Turn off gsl error handler
   gsl_set_error_handler_off();
   
   //Initiate command line interpreter and config file loader
   struct gengetopt_args_info args_info;
   struct cmdline_parser_params *configparams;
   configparams = cmdline_parser_params_create();  //initialize parameters structure
   configparams->check_required = 0;  //don't check for required values at the step
   if ( cmdline_parser_ext(argc, argv, &args_info, configparams) ) {
       fprintf(stderr, "%s: cmdline_parser() failed.\n", __func__);
       XLAL_ERROR(XLAL_FAILURE);
   }
   configparams->initialize = 0;  //don't reinitialize the parameters structure
   if ( args_info.config_given && cmdline_parser_config_file(args_info.config_arg, &args_info, configparams) ) {
      fprintf(stderr, "%s: cmdline_parser_config_file() failed.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   configparams->override = 1;  //override values in the configuration file
   configparams->check_required = 1;  //check for required values now
   if ( cmdline_parser_ext(argc, argv, &args_info, configparams) ) {
      fprintf(stderr, "%s: cmdline_parser_ext() failed.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   
   //Set lalDebugLevel to user input or 0 if no input
   lalDebugLevel = args_info.laldebug_arg;
   
   //Create directory
   mkdir(args_info.outdirectory_arg, 0777);
   snprintf(s, 1000, "%s/%s", args_info.outdirectory_arg, args_info.outfilename_arg);
   snprintf(t, 1000, "%s/%s", args_info.outdirectory_arg, args_info.ULfilename_arg);
   
   //Save args_info
   char v[1000];
   snprintf(v, 1000, "%s/%s", args_info.outdirectory_arg, args_info.configCopy_arg);
   FILE *INPUTVALS = fopen(v, "w");
   if (INPUTVALS==NULL) {
      fprintf(stderr, "%s: Could not save input parameter values.\n", __func__);
      XLAL_ERROR(XLAL_EINVAL);
   }
   if (cmdline_parser_dump(INPUTVALS, &args_info)) {
      fprintf(stderr, "%s: cmdline_parser_dump() failed.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   fclose(INPUTVALS);
   
   //Open log file
   LOG = fopen(s,"w");
   if (LOG==NULL) {
      fprintf(stderr, "%s: Log file could not be opened.\n", __func__);
      XLAL_ERROR(XLAL_EINVAL);
   }
   
   //print start time
   fprintf(stderr, "Program %s %s executed on %s", CMDLINE_PARSER_PACKAGE_NAME, CMDLINE_PARSER_VERSION, asctime(ptm));
   fprintf(LOG, "Program %s %s executed on %s", CMDLINE_PARSER_PACKAGE_NAME, CMDLINE_PARSER_VERSION, asctime(ptm));
   
   //Print out the inputs and outputs
   if (args_info.config_given) {
      fprintf(stderr, "Input parameters file: %s\n", args_info.config_arg);
      fprintf(LOG, "Input parameters file: %s\n", args_info.config_arg);
   }
   if (args_info.sftDir_given) {
      fprintf(stderr, "Input SFTs: %s/%s\n", args_info.sftDir_arg, "*.sft");
      fprintf(LOG, "Input SFTs: %s/%s\n", args_info.sftDir_arg, "*.sft");
   } else if (args_info.sftFile_given) {
      fprintf(stderr, "Input SFT file: %s\n", args_info.sftFile_arg);
      fprintf(LOG, "Input SFT file: %s\n", args_info.sftFile_arg);
   }
   fprintf(stderr, "Output directory: %s\n", args_info.outdirectory_arg);
   fprintf(LOG, "Output directory: %s\n", args_info.outdirectory_arg);
   
   //Allocate input parameters structure memory
   inputParamsStruct *inputParams = new_inputParams(args_info.IFO_given);
   if (inputParams==NULL) {
      fprintf(stderr, "%s: new_inputParams() failed.\n", __func__);
      XLAL_ERROR(XLAL_EFUNC);
   }
   
   //Read TwoSpect input parameters
   if ( (readTwoSpectInputParams(inputParams, args_info)) != 0 ) {
      fprintf(stderr, "%s: readTwoSepctInputParams() failed.\n", __func__);
      XLAL_ERROR(XLAL_EFUNC);
   }
   
   //Initialize ephemeris data structure
   EphemerisData *edat = XLALInitBarycenter(earth_ephemeris, sun_ephemeris);
   if (edat==NULL) {
      fprintf(stderr, "%s: XLALInitBarycenter() failed.\n", __func__);
      XLAL_ERROR(XLAL_EFUNC);
   }
   
   //Maximum orbital earth speed in units of c from start of S6 TwoSpect data for 104 weeks total time
   REAL4 detectorVmax = CompDetectorVmax(931081500.0+inputParams->SFToverlap, inputParams->Tcoh, inputParams->SFToverlap, 62899200.0-inputParams->SFToverlap, inputParams->det[0], edat);
   if (xlalErrno!=0) {
      fprintf(stderr, "%s: CompDetectorVmax() failed.\n", __func__);
      XLAL_ERROR(XLAL_EFUNC);
   }
   
   //Parameters for the sky-grid from a point/polygon or a sky-grid file
   if ((args_info.skyRegion_given && args_info.skyRegionFile_given) || (!args_info.skyRegion_given && !args_info.skyRegionFile_given)) {
      fprintf(stderr, "%s: You must choose either the the sky region (point or polygon) *or* a file.\n", __func__);
      XLAL_ERROR(XLAL_EINVAL);
   }
   CHAR *sky = NULL;
   if (args_info.skyRegion_given) {
      sky = XLALCalloc(strlen(args_info.skyRegion_arg)+1, sizeof(*sky));
      if (sky==NULL) {
         fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", __func__, sizeof(*sky));
         XLAL_ERROR(XLAL_ENOMEM);
      }
      sprintf(sky, "%s", args_info.skyRegion_arg);
      fprintf(LOG, "Sky region = %s\n", sky);
      fprintf(stderr, "Sky region = %s\n", sky);
   } else {
      sky = XLALCalloc(strlen(args_info.skyRegionFile_arg)+1, sizeof(*sky));
      if (sky==NULL) {
         fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", __func__, sizeof(*sky));
         XLAL_ERROR(XLAL_ENOMEM);
      }
      sprintf(sky, "%s", args_info.skyRegionFile_arg);
      fprintf(LOG, "Sky file = %s\n", sky);
      fprintf(stderr, "Sky file = %s\n", sky);
   }
   DopplerSkyScanInit scanInit = empty_DopplerSkyScanInit;
   DopplerSkyScanState scan = empty_DopplerSkyScanState;
   PulsarDopplerParams dopplerpos;
   if (args_info.skyRegion_given) {
      scanInit.gridType = GRID_ISOTROPIC;     //Default value for an approximate-isotropic grid
      scanInit.skyRegionString = sky;      //"allsky" = Default value for all-sky search
      scanInit.numSkyPartitions = 1;   //Default value so sky is not broken into chunks
      scanInit.Freq = args_info.fmin_arg+0.5*args_info.fspan_arg;  //Mid-point of the frequency band
      scanInit.dAlpha = 0.5/((inputParams->fmin+0.5*inputParams->fspan) * inputParams->Tcoh * detectorVmax);
      scanInit.dDelta = scanInit.dAlpha;
   } else {
      scanInit.gridType = GRID_FILE_SKYGRID;
      scanInit.skyGridFile = sky;
      scanInit.numSkyPartitions = 1;   //Default value so sky is not broken into chunks
      scanInit.Freq = args_info.fmin_arg+0.5*args_info.fspan_arg;  //Mid-point of the frequency band
   }
   
   //Initialize the sky-grid
   InitDopplerSkyScan(&status, &scan, &scanInit);
   if (status.statusCode!=0) {
      fprintf(stderr, "%s: InitDopplerSkyScan() failed.\n", __func__);
      XLAL_ERROR(XLAL_EFUNC);
   }
   
   //Start at first location
   if ((XLALNextDopplerSkyPos(&dopplerpos, &scan))!=0) {
      fprintf(stderr, "%s: XLALNextDopplerSkyPos() failed.\n", __func__);
      XLAL_ERROR(XLAL_EFUNC);
   }
   
   
   //Random seed value settings for random number generator
   //If the chooseSeed option was given, then:
   //seed = (IFO multiplier)*fabs(round(fmin + fspan + Pmin + Pmax + dfmin + dfmax + alpha + delta))
   UINT8 IFOmultiplier = 0;
   if (strcmp(inputParams->det[0].frDetector.prefix, "H1")==0) IFOmultiplier = 1;            //H1 gets multiplier 1
   else if (strcmp(inputParams->det[0].frDetector.prefix, "L1")==0) IFOmultiplier = 2;       //L1 gets multiplier 2
   else IFOmultiplier = 3;                                                                   //V1 gets multiplier 3
   if (args_info.randSeed_given && args_info.chooseSeed_given) inputParams->randSeed = args_info.randSeed_arg;
   else if (!args_info.randSeed_given && args_info.chooseSeed_given) inputParams->randSeed = IFOmultiplier*(UINT8)fabs(round(inputParams->fmin + inputParams->fspan + inputParams->Pmin + inputParams->Pmax + inputParams->dfmin + inputParams->dfmax + dopplerpos.Alpha + dopplerpos.Delta));
   else if (args_info.randSeed_given && !args_info.chooseSeed_given) inputParams->randSeed = args_info.randSeed_arg;
   else inputParams->randSeed = 0;
   gsl_rng_set(inputParams->rng, inputParams->randSeed);     //Set the random number generator with the given seed
   
   
   //Basic units
   REAL4 tempfspan = inputParams->fspan + 2.0*inputParams->dfmax + (inputParams->blksize-1 + 12)/inputParams->Tcoh;     //= fspan+2*dfmax+extrabins + running median blocksize-1 (Hz)
   INT4 tempnumfbins = (INT4)round(tempfspan*inputParams->Tcoh)+1;                        //= number of bins in tempfspan
   REAL8 templatefarthresh = args_info.tmplfar_arg;
   fprintf(LOG, "FAR for templates = %g\n", templatefarthresh);
   fprintf(stderr, "FAR for templates = %g\n", templatefarthresh);
   
   //Allocate memory for ffdata structure
   ffdataStruct *ffdata = new_ffdata(inputParams);
   if (ffdata==NULL) {
      fprintf(stderr, "%s: new_ffdata() failed.\n", __func__);
      XLAL_ERROR(XLAL_EFUNC);
   }
   
   //Allocate lists of candidates with initially 100 available slots (will check and rescale later, if necessary)
   //Also allocate for an upperLimitVector of length 1
   candidateVector *gaussCandidates1 = new_candidateVector(100);
   candidateVector *gaussCandidates2 = new_candidateVector(100);
   candidateVector *gaussCandidates3 = new_candidateVector(100);
   candidateVector *gaussCandidates4 = new_candidateVector(100);
   candidateVector *exactCandidates1 = new_candidateVector(100);   
   candidateVector *exactCandidates2 = new_candidateVector(100);
   candidateVector *ihsCandidates = new_candidateVector(100);
   UpperLimitVector *upperlimits = new_UpperLimitVector(1);
   if (gaussCandidates1==NULL) {
      fprintf(stderr, "%s: new_CandidateVector(%d) failed.\n", __func__, 100);
      XLAL_ERROR(XLAL_EFUNC);
   } else if (gaussCandidates2==NULL) {
      fprintf(stderr, "%s: new_CandidateVector(%d) failed.\n", __func__, 100);
      XLAL_ERROR(XLAL_EFUNC);
   } else if (gaussCandidates3==NULL) {
      fprintf(stderr, "%s: new_CandidateVector(%d) failed.\n", __func__, 100);
      XLAL_ERROR(XLAL_EFUNC);
   } else if (gaussCandidates4==NULL) {
      fprintf(stderr, "%s: new_CandidateVector(%d) failed.\n", __func__, 100);
      XLAL_ERROR(XLAL_EFUNC);
   } else if (exactCandidates1==NULL) {
      fprintf(stderr, "%s: new_CandidateVector(%d) failed.\n", __func__, 100);
      XLAL_ERROR(XLAL_EFUNC);
   } else if (exactCandidates2==NULL) {
      fprintf(stderr, "%s: new_CandidateVector(%d) failed.\n", __func__, 100);
      XLAL_ERROR(XLAL_EFUNC);
   } else if (ihsCandidates==NULL) {
      fprintf(stderr, "%s: new_CandidateVector(%d) failed.\n", __func__, 100);
      XLAL_ERROR(XLAL_EFUNC);
   } else if (upperlimits==NULL) {
      fprintf(stderr, "%s: new_UpperLimitVector(%d) failed.\n", __func__, 1);
      XLAL_ERROR(XLAL_EFUNC);
   }
   
   //Second fft plan, only need to make this once for all the exact templates
   REAL4FFTPlan *secondFFTplan = XLALCreateForwardREAL4FFTPlan(ffdata->numffts, inputParams->FFTplanFlag);
   if (secondFFTplan==NULL) {
      fprintf(stderr, "%s: XLALCreateForwardREAL4FFTPlan(%d,%d) failed.\n", __func__, ffdata->numffts, inputParams->FFTplanFlag);
      XLAL_ERROR(XLAL_EFUNC);
   }
   
   //Maximum number of IHS values to sum = twice the maximum modulation depth
   //Minimum number of IHS values to sum = twice the minimum modulation depth
   INT4 maxrows = (INT4)round(2.0*inputParams->dfmax*inputParams->Tcoh)+1;
   
   //Assume maximum bin shift possible
   inputParams->maxbinshift = (INT4)round(detectorVmax * (inputParams->fmin+0.5*inputParams->fspan) * inputParams->Tcoh)+1;
   
   //Read in the T-F data from SFTs
   fprintf(LOG, "Loading in SFTs... ");
   fprintf(stderr, "Loading in SFTs... ");
   ffdata->tfnormalization = 2.0/inputParams->Tcoh/(args_info.avesqrtSh_arg*args_info.avesqrtSh_arg);
   REAL4Vector *tfdata = readInSFTs(inputParams, &(ffdata->tfnormalization));
   if (tfdata==NULL) {
      fprintf(stderr, "\n%s: readInSFTs() failed.\n", __func__);
      XLAL_ERROR(XLAL_EFUNC);
   }
   fprintf(LOG, "done\n");
   fprintf(stderr, "done\n");
   //comment this
   /* FILE *rawtfdata = fopen("./output/rawtfdata.dat","w");
   for (ii=0; ii<(INT4)tfdata->length; ii++) fprintf(rawtfdata, "%f\n", tfdata->data[ii]);
   fclose(rawtfdata); */
   
   //Print SFT times, if requested by user
   if (args_info.printSFTtimes_given) {
      char w[1000];
      snprintf(w, 1000, "%s/%s", args_info.outdirectory_arg, "inputSFTtimes.dat");
      FILE *INSFTTIMES = fopen(w, "w");
      INT4 sftlength = tfdata->length/ffdata->numffts;
      for (ii=0; ii<ffdata->numffts; ii++) {
         if (tfdata->data[ii*sftlength]!=0.0) fprintf(INSFTTIMES, "%9d 0\n", (INT4)round(inputParams->searchstarttime+ii*(inputParams->Tcoh-inputParams->SFToverlap)));
      }
      fclose(INSFTTIMES);
   }
   
   //Removing bad SFTs using K-S test and Kuiper's test
   if (inputParams->markBadSFTs!=0) {
      fprintf(stderr, "Marking and removing bad SFTs... ");
      INT4Vector *removeTheseSFTs = markBadSFTs(tfdata, inputParams);
      if (removeTheseSFTs==NULL) {
         fprintf(stderr, "%s: markBadSFTs() failed.\n", __func__);
         XLAL_ERROR(XLAL_EFUNC);
      }
      removeBadSFTs(tfdata, removeTheseSFTs);
      fprintf(stderr, "done.\n");
      XLALDestroyINT4Vector(removeTheseSFTs);
   }
   /* FILE *rawtfdata = fopen("./output/rawtfdata.dat","w");
   for (ii=0; ii<(INT4)tfdata->length; ii++) fprintf(rawtfdata, "%f\n", tfdata->data[ii]);
   fclose(rawtfdata); */
   if (args_info.printUsedSFTtimes_given) {
      char w[1000];
      snprintf(w, 1000, "%s/%s", args_info.outdirectory_arg, "usedSFTtimes.dat");
      FILE *USEDSFTTIMES = fopen(w, "w");
      INT4 sftlength = tfdata->length/ffdata->numffts;
      for (ii=0; ii<ffdata->numffts; ii++) {
         if (tfdata->data[ii*sftlength]!=0.0) fprintf(USEDSFTTIMES, "%9d 0\n", (INT4)round(inputParams->searchstarttime+ii*(inputParams->Tcoh-inputParams->SFToverlap)));
      }
      fclose(USEDSFTTIMES);
   }
   
   
   //Calculate the running mean values of the SFTs (output here is smaller than initialTFdata). Here,
   //numfbins needs to be the bins you expect to come out of the running means -- the band you are going
   //to search!
   fprintf(LOG, "Assessing background... ");
   fprintf(stderr, "Assessing background... ");
   REAL4Vector *background = XLALCreateREAL4Vector(ffdata->numffts*(ffdata->numfbins + 2*inputParams->maxbinshift));
   if (background==NULL) {
      fprintf(stderr, "\n%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ffdata->numffts*(ffdata->numfbins + 2*inputParams->maxbinshift));
      XLAL_ERROR(XLAL_EFUNC);
   }
   tfRngMeans(background, tfdata, ffdata->numffts, ffdata->numfbins + 2*inputParams->maxbinshift, inputParams->blksize);
   if (xlalErrno!=0) {
      fprintf(stderr, "\n%s: tfRngMeans() failed.\n", __func__);
      XLAL_ERROR(XLAL_EFUNC);
   }
   
   //Existing SFTs listed in this vector
   INT4Vector *sftexist = existingSFTs(tfdata, inputParams, ffdata->numfbins, ffdata->numffts);
   if (sftexist==NULL) {
      fprintf(stderr, "\n%s: existingSFTs() failed.\n", __func__);
      XLAL_ERROR(XLAL_EFUNC);
   }
   INT4 totalincludedsftnumber = 0.0;
   for (ii=0; ii<(INT4)sftexist->length; ii++) if (sftexist->data[ii]==1) totalincludedsftnumber++;
   REAL4 frac_tobs_complete = (REAL4)totalincludedsftnumber/(REAL4)sftexist->length;
   if (frac_tobs_complete<0.1) {
      fprintf(stderr, "%s: The useable SFTs cover less than 10 percent of the total observation time\n", __func__);
      fprintf(LOG, "%s: The useable SFTs cover less than 10 percent of the total observation time\n", __func__);
      return 0;
   }
   
   //Index values of existing SFTs
   INT4Vector *indexValuesOfExistingSFTs = XLALCreateINT4Vector(totalincludedsftnumber);
   if (indexValuesOfExistingSFTs==NULL) {
      fprintf(stderr, "\n%s: XLALCreateINT4Vector(%d) failed.\n", __func__, totalincludedsftnumber);
      XLAL_ERROR(XLAL_EFUNC);
   }
   jj = 0;
   for (ii=0; ii<(INT4)sftexist->length; ii++) {
      if (sftexist->data[ii] == 1) {
         indexValuesOfExistingSFTs->data[jj] = ii;
         jj++;
      }
   }
   
   //I wrote this to compensate for a bad input of the expected noise floor
   REAL8 backgroundmeannormfactor = 0.0;
   /* INT8 avefact = 0;
   for (ii=0; ii<(INT4)background->length; ii++) {
      if (background->data[ii]!=0.0) {
         backgroundmeannormfactor += background->data[ii];
         avefact++;
      }
   }
   backgroundmeannormfactor = (REAL8)avefact/backgroundmeannormfactor; */
   backgroundmeannormfactor = 1.0/calcMedian_ignoreZeros(background);
   ffdata->tfnormalization *= backgroundmeannormfactor;
   fprintf(LOG, "done\n");
   fprintf(stderr, "done\n");
   
   //If necessary, open the NORMRMSOUT file
   if (args_info.normRMSoutput_given) {
      snprintf(u, 1000, "%s/%s", args_info.outdirectory_arg, args_info.normRMSoutput_arg);
      NORMRMSOUT = fopen(u,"w");
      if (NORMRMSOUT==NULL) {
         fprintf(stderr, "%s: normalized RMS data file could not be opened for writing.\n", __func__);
         XLAL_ERROR(XLAL_EINVAL);
      }
   }
   
   //Line detection
   INT4Vector *lines = NULL;
   INT4 heavilyContaminatedBand = 0;
   if (args_info.lineDetection_given) {
      lines = detectLines_simple(tfdata, ffdata, inputParams);
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
   
   //Close the NORMRMSOUT file, if necessary
   if (args_info.normRMSoutput_given) fclose(NORMRMSOUT);
   
   //If the band is heavily contaminated by lines, don't do any follow up.
   if (heavilyContaminatedBand) args_info.IHSonly_given = 1;
   
   //Need to reduce the original TF data to remove the excess bins used for running median calculation. Normalize the TF as the same as the background was normalized
   REAL4Vector *usableTFdata = XLALCreateREAL4Vector(background->length);
   if (usableTFdata==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, background->length);
      XLAL_ERROR(XLAL_EFUNC);
   }
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
   
   
   //Do mean subtraction of TFdata here--modifies the usableTFdata vector!!!
   tfMeanSubtract(usableTFdata, background, ffdata->numffts, ffdata->numfbins+2*inputParams->maxbinshift);
   //fprintf(stderr,"%g\n",calcMean(usableTFdata));
   
   fprintf(LOG, "Maximum row width to be searched = %d\n", maxrows);
   fprintf(stderr, "Maximum row width to be searched = %d\n", maxrows);
   
   //Initialize reused values
   ihsMaximaStruct *ihsmaxima = new_ihsMaxima(ffdata->numfbins, maxrows);
   ihsfarStruct *ihsfarstruct = new_ihsfarStruct(maxrows, inputParams);
   REAL4Vector *detectorVelocities = XLALCreateREAL4Vector(ffdata->numffts);
   INT4Vector *binshifts = XLALCreateINT4Vector(ffdata->numffts);
   REAL4Vector *aveNoise = XLALCreateREAL4Vector(ffdata->numfprbins);
   if (ihsmaxima==NULL) {
      fprintf(stderr, "%s: new_ihsMaxima(%d,%d) failed.\n", __func__, ffdata->numfbins, maxrows);
      XLAL_ERROR(XLAL_EFUNC);
   } else if (ihsfarstruct==NULL) {
      fprintf(stderr, "%s: new_ihsfarStruct(%d) failed.\n", __func__, maxrows);
      XLAL_ERROR(XLAL_EFUNC);
   } else if (detectorVelocities==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ffdata->numffts);
      XLAL_ERROR(XLAL_EFUNC);
   } else if (binshifts==NULL) {
      fprintf(stderr, "%s: XLALCreateINT4Vector(%d) failed.\n", __func__, ffdata->numffts);
      XLAL_ERROR(XLAL_EFUNC);
   } else if (aveNoise==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL8Vector(%d) failed.\n", __func__, ffdata->numfprbins);
      XLAL_ERROR(XLAL_EFUNC);
   } 
   
   //Initialize to zero for far just at the start
   ihsfarstruct->ihsfar->data[0] = 0.0;
   REAL4 antweightsrms = 0.0;
   
   //Davies' algorithm uses an internal error code
   //Later upgrade: remove this because we do our own error handling
   INT4 proberrcode = 0;
   ffdata->tfnormalization *= 0.5*inputParams->Tcoh;
   
   //Print message that we start the analysis
   fprintf(LOG, "Starting TwoSpect analysis...\n");
   fprintf(stderr, "Starting TwoSpect analysis...\n");
   
   
   //Antenna normalization (determined from injections on H1 at ra=0, dec=0, with circular polarization)
   //When doing linear polarizations, the IHS factor needs to be 25.2*1.082 and this antenna weights
   //function needs to be set to use linear polarization.
   REAL4Vector *antweightsforihs2h0 = XLALCreateREAL4Vector(ffdata->numffts);
   CompAntennaPatternWeights(antweightsforihs2h0, 0.0, 0.0, inputParams->searchstarttime, inputParams->Tcoh, inputParams->SFToverlap, inputParams->Tobs, 0, 0.0, lalCachedDetectors[LAL_LHO_4K_DETECTOR]);
   if (xlalErrno!=0) {
      fprintf(stderr, "%s: CompAntennaPatternWeights() failed.\n", __func__);
      XLAL_ERROR(XLAL_EFUNC);
   }
   
   
   //Search over the sky region (outer loop of single TwoSpect program instance)
   while (scan.state != STATE_FINISHED) {
      fprintf(LOG, "Sky location: RA = %g, DEC = %g\n", dopplerpos.Alpha, dopplerpos.Delta);
      fprintf(stderr, "Sky location: RA = %g, DEC = %g\n", dopplerpos.Alpha, dopplerpos.Delta);
      
      //Determine detector velocity w.r.t. a sky location for each SFT
      CompAntennaVelocity(detectorVelocities, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, inputParams->searchstarttime, inputParams->Tcoh, inputParams->SFToverlap, inputParams->Tobs, inputParams->det[0], edat);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: CompAntennaVelocity() failed.\n", __func__);
         XLAL_ERROR(XLAL_EFUNC);
      }
      
      //Compute the bin shifts for each SFT
      CompBinShifts(binshifts, inputParams->fmin+0.5*inputParams->fspan, detectorVelocities, inputParams->Tcoh, inputParams->dopplerMultiplier);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: CompBinShifts() failed.\n", __func__);
         XLAL_ERROR(XLAL_EFUNC);
      }
      
      //Track identified lines
      REAL4VectorSequence *trackedlines = NULL;
      if (lines!=NULL) {
         trackedlines = trackLines(lines, binshifts, inputParams);
      }
      
      //Compute antenna pattern weights. If antennaOff input flag is given, then set all values equal to 1.0
      REAL4Vector *antweights = XLALCreateREAL4Vector(ffdata->numffts);
      if (antweights==NULL) {
         fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ffdata->numffts);
         XLAL_ERROR(XLAL_EFUNC);
      }
      if (args_info.antennaOff_given) {
         for (ii=0; ii<(INT4)antweights->length; ii++) antweights->data[ii] = 1.0;
      } else {
         if (args_info.linPolAngle_given) CompAntennaPatternWeights(antweights, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, inputParams->searchstarttime, inputParams->Tcoh, inputParams->SFToverlap, inputParams->Tobs, 1, args_info.linPolAngle_arg, inputParams->det[0]);
         else CompAntennaPatternWeights(antweights, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, inputParams->searchstarttime, inputParams->Tcoh, inputParams->SFToverlap, inputParams->Tobs, 0, 0.0, inputParams->det[0]);
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: CompAntennaPatternWeights() failed.\n", __func__);
            XLAL_ERROR(XLAL_EFUNC);
         }
      }
      
      //Calculate antenna RMS value
      REAL4 currentAntWeightsRMS = calcRms(antweights);
      if (XLAL_IS_REAL4_FAIL_NAN(currentAntWeightsRMS)) {
         fprintf(stderr, "%s: calcRms() failed.\n", __func__);
         XLAL_ERROR(XLAL_EFUNC);
      }
      
      //Slide SFTs here -- need to slide the data and the estimated background
      REAL4Vector *TFdata_slided = XLALCreateREAL4Vector(ffdata->numffts*ffdata->numfbins);
      if (TFdata_slided==NULL) {
         fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ffdata->numffts*ffdata->numfbins);
         XLAL_ERROR(XLAL_EFUNC);
      }
      REAL4Vector *background_slided = XLALCreateREAL4Vector(TFdata_slided->length);
      if (background_slided==NULL) {
         fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, TFdata_slided->length);
         XLAL_ERROR(XLAL_EFUNC);
      }
      slideTFdata(TFdata_slided, inputParams, usableTFdata, binshifts);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: slideTFdata() failed.\n", __func__);
         XLAL_ERROR(XLAL_EFUNC);
      }
      slideTFdata(background_slided, inputParams, background, binshifts);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: slideTFdata() failed.\n", __func__);
         XLAL_ERROR(XLAL_EFUNC);
      }
      //fprintf(stderr, "Mean of TFdata_slided %g, background_slided %g\n", calcMean(TFdata_slided), calcMean(background_slided));
      /* FILE *TFBACKGROUND = fopen("./output/tfbackground.dat","w");
      for (ii=0; ii<(INT4)background_slided->length; ii++) fprintf(TFBACKGROUND, "%.6g\n", background_slided->data[ii]);
      fclose(TFBACKGROUND); */
      /* FILE *TFSLIDED = fopen("./output/tfslided.dat","w");
      for (ii=0; ii<(INT4)TFdata_slided->length; ii++) fprintf(TFSLIDED, "%.6g\n", TFdata_slided->data[ii]);
      fclose(TFSLIDED); */
      
      //Check the RMS of the antenna weights; if a deviation of more than 1%, then reset the IHS FAR
      if (antweightsrms == 0.0) antweightsrms = currentAntWeightsRMS;
      if ( fabs(currentAntWeightsRMS-antweightsrms)/antweightsrms >= 0.01 ) {
         ihsfarstruct->ihsfar->data[0] = 0.0;
         antweightsrms = currentAntWeightsRMS;
      }
      
      //Antenna normalization for different sky locations
      REAL8 skypointffnormalization = 1.0;
      ffPlaneNoise(aveNoise, inputParams, sftexist, background_slided, antweightsforihs2h0, secondFFTplan, &(skypointffnormalization));
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: ffPlaneNoise() failed.\n", __func__);
         XLAL_ERROR(XLAL_EFUNC);
      }
      
      //Average noise floor of FF plane for each 1st FFT frequency bin
      ffdata->ffnormalization = 1.0;
      ffPlaneNoise(aveNoise, inputParams, sftexist, background_slided, antweights, secondFFTplan, &(ffdata->ffnormalization));
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: ffPlaneNoise() failed.\n", __func__);
         XLAL_ERROR(XLAL_EFUNC);
      }
      
      //Compute the weighted TF data
      REAL4Vector *TFdata_weighted = XLALCreateREAL4Vector(ffdata->numffts*ffdata->numfbins);
      if (TFdata_weighted==NULL) {
         fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ffdata->numffts*ffdata->numfbins);
         XLAL_ERROR(XLAL_EFUNC);
      }
      tfWeight(TFdata_weighted, TFdata_slided, background_slided, antweights, indexValuesOfExistingSFTs, inputParams);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: tfWeight() failed.\n", __func__);
         XLAL_ERROR(XLAL_EFUNC);
      }
      XLALDestroyREAL4Vector(TFdata_slided);
      //XLALDestroyREAL4Vector(background_slided);
      XLALDestroyREAL4Vector(antweights);
      /* FILE *TFDATA = fopen("./output/tfdata.dat","w");
      for (jj=0; jj<(INT4)TFdata_weighted->length; jj++) fprintf(TFDATA,"%.6f\n",TFdata_weighted->data[jj]);
      fclose(TFDATA); */
      
      //Calculation of average TF noise per frequency bin ratio to total mean
      //this block of code does not avoid lines when computing the average F-bin ratio. Upper limits remain virtually unchanged
      //when comaring runs that have line finding enabled or disabled
      REAL4Vector *aveTFnoisePerFbinRatio = XLALCreateREAL4Vector(ffdata->numfbins);
      REAL4Vector *TSofPowers = XLALCreateREAL4Vector(ffdata->numffts);
      if (aveTFnoisePerFbinRatio==NULL) {
         fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ffdata->numfbins);
         XLAL_ERROR(XLAL_EFUNC);
      } else if (TSofPowers==NULL) {
         fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ffdata->numffts);
         XLAL_ERROR(XLAL_EFUNC);
      }
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
      makeSecondFFT(ffdata, TFdata_weighted, secondFFTplan);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: makeSecondFFT() failed.\n", __func__);
         XLAL_ERROR(XLAL_EFUNC);
      }
      
      REAL4 secFFTmean = calcMean(ffdata->ffdata);
      REAL4 secFFTsigma = calcStddev(ffdata->ffdata);
      
      XLALDestroyREAL4Vector(TFdata_weighted);
      TFdata_weighted = NULL;
      
      fprintf(stderr, "2nd FFT ave = %g, 2nd FFT stddev = %g, expected ave = %g\n", secFFTmean, secFFTsigma, 1.0);
      //comment this out
      /* FILE *FFDATA = fopen("./output/ffdata.dat","w");
      for (jj=0; jj<(INT4)ffdata->ffdata->length; jj++) fprintf(FFDATA,"%g\n",ffdata->ffdata->data[jj]);
      fclose(FFDATA); */
      
      //Exit with failure if there are no SFTs (probably this doesn't get hit)
      if (secFFTmean==0.0) {
         fprintf(stderr, "%s: Average second FFT power is 0.0. Perhaps no SFTs are remaining? Program exiting with failure.\n", __func__);
         fprintf(LOG, "%s: Average second FFT power is 0.0. Perhaps no SFTs are remaining? Program exiting with failure.\n", __func__);
         XLAL_ERROR(XLAL_FAILURE);
      }
      
      
////////Start of the IHS step!
      //Find the FAR of IHS sum
      if (ihsfarstruct->ihsfar->data[0]==0.0) {
         fprintf(stderr, "Determining IHS FAR values... ");
         fprintf(LOG, "Determining IHS FAR values... ");
         genIhsFar(ihsfarstruct, inputParams, maxrows, aveNoise);
         if (xlalErrno!=0) {
            fprintf(stderr,"\n%s: genIhsFar() failed.\n", __func__);
            XLAL_ERROR(XLAL_EFUNC);
         }
         fprintf(stderr, "done.\n");
         fprintf(LOG, "done.\n");
      }
      
      //Run the IHS algorithm on the data
      runIHS(ihsmaxima, ffdata, ihsfarstruct, inputParams, maxrows, aveTFnoisePerFbinRatio);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: runIHS() failed.\n", __func__);
         XLAL_ERROR(XLAL_EFUNC);
      }
      
      //Find any IHS candidates
      findIHScandidates(ihsCandidates, ihsfarstruct, inputParams, ffdata, ihsmaxima, aveTFnoisePerFbinRatio, trackedlines);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: findIHScandidates() failed.\n", __func__);
         XLAL_ERROR(XLAL_EFUNC);
      }
      fprintf(LOG, "Candidates found in IHS step = %d\n", ihsCandidates->numofcandidates);
      fprintf(stderr, "Candidates found in IHS step = %d\n", ihsCandidates->numofcandidates);
      //for (ii=0; ii<(INT4)ihsCandidates->numofcandidates; ii++) fprintf(stderr, "%d %g %g %g %g\n", ii, ihsCandidates->data[ii].fsig, ihsCandidates->data[ii].period, ihsCandidates->data[ii].moddepth, ihsCandidates->data[ii].prob);  //comment this
      
      //If requested, keep only the most significant IHS candidates
      candidateVector *ihsCandidates_reduced = NULL;
      if (args_info.keepOnlyTopNumIHS_given && (INT4)ihsCandidates->numofcandidates>args_info.keepOnlyTopNumIHS_arg) {
         fprintf(stderr, "Reducing total number of IHS candidates %d to user input %d\n", ihsCandidates->numofcandidates, args_info.keepOnlyTopNumIHS_arg);
         fprintf(LOG, "Reducing total number of IHS candidates %d to user input %d\n", ihsCandidates->numofcandidates, args_info.keepOnlyTopNumIHS_arg);
         ihsCandidates_reduced = keepMostSignificantCandidates(ihsCandidates, inputParams);
         if (ihsCandidates_reduced==NULL) {
            fprintf(stderr,"%s: keepMostSignificantCandidates() failed.\n", __func__);
            XLAL_ERROR(XLAL_EFUNC);
         }
         //for (ii=0; ii<(INT4)ihsCandidates_reduced->numofcandidates; ii++) fprintf(stderr, "%d %g %g %g %g\n", ii, ihsCandidates_reduced->data[ii].fsig, ihsCandidates_reduced->data[ii].period, ihsCandidates_reduced->data[ii].moddepth, ihsCandidates_reduced->data[ii].prob);  //comment this
      }
      
////////End of the IHS step
      
////////Start of the Gaussian template search!
      if (args_info.IHSonly_given) {
         
         if (args_info.keepOnlyTopNumIHS_given && (INT4)ihsCandidates->numofcandidates>args_info.keepOnlyTopNumIHS_arg) {
            if (exactCandidates2->length < exactCandidates2->numofcandidates+ihsCandidates_reduced->numofcandidates) {
               exactCandidates2 = resize_candidateVector(exactCandidates2, exactCandidates2->numofcandidates+ihsCandidates_reduced->numofcandidates);
               if (exactCandidates2->data==NULL) {
                  fprintf(stderr,"%s: resize_candidateVector(%d) failed.\n", __func__, ihsCandidates_reduced->numofcandidates);
                  XLAL_ERROR(XLAL_EFUNC);
               }
            }
            
            INT4 numofcandidatesalready = exactCandidates2->numofcandidates;
            for (ii=0; ii<(INT4)ihsCandidates_reduced->numofcandidates; ii++) {
               loadCandidateData(&(exactCandidates2->data[ii+numofcandidatesalready]), ihsCandidates_reduced->data[ii].fsig, ihsCandidates_reduced->data[ii].period, ihsCandidates_reduced->data[ii].moddepth, dopplerpos.Alpha, dopplerpos.Delta, ihsCandidates_reduced->data[ii].stat, ihsCandidates_reduced->data[ii].h0, ihsCandidates_reduced->data[ii].prob, 0, ihsCandidates_reduced->data[ii].normalization);
               exactCandidates2->data[ii+numofcandidatesalready].h0 /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25); //Scaling here
               (exactCandidates2->numofcandidates)++;
            }
            
         } else {
            if (exactCandidates2->length < exactCandidates2->numofcandidates+ihsCandidates->numofcandidates) {
               exactCandidates2 = resize_candidateVector(exactCandidates2, exactCandidates2->numofcandidates+ihsCandidates->numofcandidates);
               if (exactCandidates2->data==NULL) {
                  fprintf(stderr,"%s: resize_candidateVector(%d) failed.\n", __func__, ihsCandidates->numofcandidates);
                  XLAL_ERROR(XLAL_EFUNC);
               }
            }
            
            INT4 numofcandidatesalready = exactCandidates2->numofcandidates;
            for (ii=0; ii<(INT4)ihsCandidates->numofcandidates; ii++) {
               loadCandidateData(&(exactCandidates2->data[ii+numofcandidatesalready]), ihsCandidates->data[ii].fsig, ihsCandidates->data[ii].period, ihsCandidates->data[ii].moddepth, dopplerpos.Alpha, dopplerpos.Delta, ihsCandidates->data[ii].stat, ihsCandidates->data[ii].h0, ihsCandidates->data[ii].prob, 0, ihsCandidates->data[ii].normalization);
               exactCandidates2->data[ii+numofcandidatesalready].h0 /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25); //Scaling here
               (exactCandidates2->numofcandidates)++;
            }
            
         }
         
      } else if ((!args_info.simpleBandRejection_given || (args_info.simpleBandRejection_given && secFFTsigma<inputParams->simpleSigmaExclusion))) {
         
         if (args_info.keepOnlyTopNumIHS_given && (INT4)ihsCandidates->numofcandidates>args_info.keepOnlyTopNumIHS_arg) {
            //Test the IHS candidates against Gaussian templates in this function
            if ( testIHScandidates(gaussCandidates1, ihsCandidates_reduced, ffdata, aveNoise, aveTFnoisePerFbinRatio, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, inputParams) != 0 ) {
               fprintf(stderr, "%s: testIHScandidates() failed.\n", __func__);
               XLAL_ERROR(XLAL_EFUNC);
            }
         } else {
            //Test the IHS candidates against Gaussian templates in this function
            if ( testIHScandidates(gaussCandidates1, ihsCandidates, ffdata, aveNoise, aveTFnoisePerFbinRatio, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, inputParams) != 0 ) {
               fprintf(stderr, "%s: testIHScandidates() failed.\n", __func__);
               XLAL_ERROR(XLAL_EFUNC);
            }
         }
         
         fprintf(LOG,"Initial stage done with candidates = %d\n",gaussCandidates1->numofcandidates);
         fprintf(stderr,"Initial stage done with candidates = %d\n",gaussCandidates1->numofcandidates);
         
         for (ii=0; ii<(INT4)gaussCandidates1->numofcandidates; ii++) fprintf(stderr, "Candidate %d: f0=%g, P=%g, df=%g\n", ii, gaussCandidates1->data[ii].fsig, gaussCandidates1->data[ii].period, gaussCandidates1->data[ii].moddepth);
      } /* if IHSonly is not given */
      
      if (args_info.keepOnlyTopNumIHS_given && (INT4)ihsCandidates->numofcandidates>args_info.keepOnlyTopNumIHS_arg) free_candidateVector(ihsCandidates_reduced);
      
////////End of the Gaussian template search

      //Reset IHS candidates, but keep length the same (doesn't reset actual values in the vector)
      ihsCandidates->numofcandidates = 0;
      
      //Search the IHS templates further if user has not specified IHSonly flag
      //if (!args_info.IHSonly_given && ( !args_info.simpleBandRejection_given || ( args_info.simpleBandRejection_given && secFFTsigma<inputParams->simpleSigmaExclusion ) ) && ( !args_info.lineDetection_given || (args_info.lineDetection_given && lines==NULL) )) {
      if (gaussCandidates1->numofcandidates>0) {
////////Start clustering! Note that the clustering algorithm takes care of the period range of parameter space
         clusterCandidates(gaussCandidates2, gaussCandidates1, ffdata, inputParams, aveNoise, aveTFnoisePerFbinRatio, sftexist, 0);
         if (xlalErrno!=0) {
            fprintf(stderr,"%s: clusterCandidates() failed.\n", __func__);
            XLAL_ERROR(XLAL_EFUNC);
         }
         fprintf(LOG, "Clustering done with candidates = %d\n", gaussCandidates2->numofcandidates);
         fprintf(stderr, "Clustering done with candidates = %d\n", gaussCandidates2->numofcandidates);
         
         for (ii=0; ii<(INT4)gaussCandidates2->numofcandidates; ii++) fprintf(stderr, "Candidate %d: f0=%g, P=%g, df=%g\n", ii, gaussCandidates2->data[ii].fsig, gaussCandidates2->data[ii].period, gaussCandidates2->data[ii].moddepth);
////////End clustering
         
         //Reset first set of Gaussian template candidates, but keep length the same (doesn't reset actual values in the vector)
         gaussCandidates1->numofcandidates = 0;
         
////////Start detailed Gaussian template search!
         //REAL4 tcohfactor = 1.49e-3*inputParams->Tcoh + 1.76;
         for (ii=0; ii<(INT4)gaussCandidates2->numofcandidates; ii++) {
            
            if (gaussCandidates3->numofcandidates == gaussCandidates3->length-1) {
               gaussCandidates3 = resize_candidateVector(gaussCandidates3, 2*gaussCandidates3->length);
               if (gaussCandidates3->data==NULL) {
                  fprintf(stderr,"%s: resize_candidateVector(%d) failed.\n", __func__, 2*gaussCandidates3->length);
                  XLAL_ERROR(XLAL_EFUNC);
               }
            }
            //bruteForceTemplateSearch(&(gaussCandidates3->data[gaussCandidates3->numofcandidates]), gaussCandidates2->data[ii], gaussCandidates2->data[ii].fsig-2.5/inputParams->Tcoh, gaussCandidates2->data[ii].fsig+2.5/inputParams->Tcoh, 11, 5, gaussCandidates2->data[ii].moddepth-2.5/inputParams->Tcoh, gaussCandidates2->data[ii].moddepth+2.5/inputParams->Tcoh, 11, inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 0);
            //bruteForceTemplateSearch(&(gaussCandidates3->data[gaussCandidates3->numofcandidates]), gaussCandidates2->data[ii], gaussCandidates2->data[ii].fsig-1.5/inputParams->Tcoh, gaussCandidates2->data[ii].fsig+1.5/inputParams->Tcoh, 7, 5, gaussCandidates2->data[ii].moddepth-1.5/inputParams->Tcoh, gaussCandidates2->data[ii].moddepth+1.5/inputParams->Tcoh, 7, inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 0);
            bruteForceTemplateSearch(&(gaussCandidates3->data[gaussCandidates3->numofcandidates]), //Output candidate
                                     gaussCandidates2->data[ii],                             //Candidate
                                     gaussCandidates2->data[ii].fsig-1.0/inputParams->Tcoh,  //Minimum frequency
                                     gaussCandidates2->data[ii].fsig+1.0/inputParams->Tcoh,  //Maximum frequency
                                     5,                                                      //Number of frequencies to search in range
                                     5,                                                      //Number of periods to search
                                     gaussCandidates2->data[ii].moddepth-1.0/inputParams->Tcoh, //Minimum modulation depth
                                     gaussCandidates2->data[ii].moddepth+1.0/inputParams->Tcoh, //Maximum modulation depth
                                     5,                                                      //Number of modulation depths to search in range
                                     inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 0);
            if (xlalErrno!=0) {
               fprintf(stderr, "%s: bruteForceTemplateSearch() failed.\n", __func__);
               XLAL_ERROR(XLAL_EFUNC);
            }
            gaussCandidates3->numofcandidates++;
            
         } /* for ii < numofcandidates */
          
         for (ii=0; ii<(INT4)gaussCandidates3->numofcandidates; ii++) fprintf(stderr,"Candidate %d: f0=%g, P=%g, df=%g\n", ii, gaussCandidates3->data[ii].fsig, gaussCandidates3->data[ii].period, gaussCandidates3->data[ii].moddepth);
////////End detailed Gaussian template search

         //Reset 2nd round of Gaussian template candidates, but keep length the same (doesn't reset actual values in the vector)
         gaussCandidates2->numofcandidates = 0;

////////Start clustering!
         clusterCandidates(gaussCandidates4, gaussCandidates3, ffdata, inputParams, aveNoise, aveTFnoisePerFbinRatio, sftexist, 0);
         if (xlalErrno!=0) {
            fprintf(stderr,"%s: clusterCandidates() failed.\n", __func__);
            XLAL_ERROR(XLAL_EFUNC);
         }
         fprintf(LOG, "Clustering done with candidates = %d\n", gaussCandidates4->numofcandidates);
         fprintf(stderr, "Clustering done with candidates = %d\n", gaussCandidates4->numofcandidates);
         
         for (ii=0; ii<(INT4)gaussCandidates4->numofcandidates; ii++) fprintf(stderr, "Candidate %d: f0=%g, P=%g, df=%g\n", ii, gaussCandidates4->data[ii].fsig, gaussCandidates4->data[ii].period, gaussCandidates4->data[ii].moddepth);
////////End clustering
         
         //Reset 3rd set of Gaussian template candidates
         gaussCandidates3->numofcandidates = 0;

////////Initial check using "exact" template
         for (ii=0; ii<(INT4)gaussCandidates4->numofcandidates; ii++) {
            
            templateStruct *template = new_templateStruct(inputParams->maxtemplatelength);
            if (template==NULL) {
               fprintf(stderr,"%s: new_templateStruct(%d) failed.\n", __func__, inputParams->maxtemplatelength);
               XLAL_ERROR(XLAL_EFUNC); 
            }
            
            if (!args_info.gaussTemplatesOnly_given) {
               makeTemplate(template, gaussCandidates4->data[ii], inputParams, sftexist, secondFFTplan);
               if (xlalErrno!=0) {
                  fprintf(stderr,"%s: makeTemplate() failed.\n", __func__);
                  XLAL_ERROR(XLAL_EFUNC);
               }
            } else {
               makeTemplateGaussians(template, gaussCandidates4->data[ii], inputParams, ffdata->numfbins, ffdata->numfprbins);
               if (xlalErrno!=0) {
                  fprintf(stderr,"%s: makeTemplateGaussians() failed.\n", __func__);
                  XLAL_ERROR(XLAL_EFUNC);
               }
            }
            
            farStruct *farval = new_farStruct();
            if (farval==NULL) {
               fprintf(stderr,"%s: new_farStruct() failed.\n", __func__);
               XLAL_ERROR(XLAL_EFUNC); 
            }
            
            if (inputParams->calcRthreshold) {
               numericFAR(farval, template, templatefarthresh, aveNoise, aveTFnoisePerFbinRatio, inputParams, inputParams->rootFindingMethod);
               if (xlalErrno!=0) {
                  fprintf(stderr,"%s: numericFAR() failed.\n", __func__);
                  XLAL_ERROR(XLAL_EFUNC);
               }
            }
            
            REAL8 R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
            if (XLAL_IS_REAL8_FAIL_NAN(R)) {
               fprintf(stderr,"%s: calculateR() failed.\n", __func__);
               XLAL_ERROR(XLAL_EFUNC);
            }
            REAL8 prob = 0.0;
            if ( R > 0.0 ) {
               prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, inputParams, &proberrcode);
               if (XLAL_IS_REAL8_FAIL_NAN(prob)) {
                  fprintf(stderr,"%s: probR() failed.\n", __func__);
                  XLAL_ERROR(XLAL_EFUNC);
               }
            }
            REAL8 h0 = 2.7426*pow(R/(inputParams->Tcoh*inputParams->Tobs),0.25);
            if ((!inputParams->calcRthreshold && prob<log10(templatefarthresh)) || (inputParams->calcRthreshold && R>farval->far)) {
               if (exactCandidates1->numofcandidates == exactCandidates1->length-1) {
                  exactCandidates1 = resize_candidateVector(exactCandidates1, 2*exactCandidates1->length);
                  if (exactCandidates1->data==NULL) {
                     fprintf(stderr,"%s: resize_candidateVector(%d) failed.\n", __func__, 2*exactCandidates1->length);
                     XLAL_ERROR(XLAL_EFUNC);
                  }
               }
               loadCandidateData(&exactCandidates1->data[exactCandidates1->numofcandidates], gaussCandidates4->data[ii].fsig, gaussCandidates4->data[ii].period, gaussCandidates4->data[ii].moddepth, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, R, h0, prob, proberrcode, gaussCandidates4->data[ii].normalization);
               exactCandidates1->numofcandidates++;
            }
            
            free_templateStruct(template);
            template = NULL;
            free_farStruct(farval);
            farval = NULL;
         } /* for ii < numofcandidates */
         fprintf(LOG, "Number of candidates confirmed with exact templates = %d\n", exactCandidates1->numofcandidates);
         fprintf(stderr, "Number of candidates confirmed with exact templates = %d\n", exactCandidates1->numofcandidates);
         for (ii=0; ii<(INT4)exactCandidates1->numofcandidates; ii++) fprintf(stderr, "Candidate %d: f0=%g, P=%g, df=%g\n", ii, exactCandidates1->data[ii].fsig, exactCandidates1->data[ii].period, exactCandidates1->data[ii].moddepth);
////////Done with initial check
         
         //Reset 4th set of Gaussian template candidates, but keep length the same (doesn't reset actual values in the vector)
         gaussCandidates4->numofcandidates = 0;

////////Start detailed "exact" template search!
         for (ii=0; ii<(INT4)exactCandidates1->numofcandidates; ii++) {
            
            if (exactCandidates2->numofcandidates == exactCandidates2->length-1) {
               exactCandidates2 = resize_candidateVector(exactCandidates2, 2*exactCandidates2->length);
               if (exactCandidates2->data==NULL) {
                  fprintf(stderr,"%s: resize_candidateVector(%d) failed.\n", __func__, 2*exactCandidates2->length);
                  XLAL_ERROR(XLAL_EFUNC);
               }
            }
            
            if (!args_info.gaussTemplatesOnly_given) {
               //bruteForceTemplateSearch(&(exactCandidates2->data[exactCandidates2->numofcandidates]), exactCandidates1->data[ii], exactCandidates1->data[ii].fsig-1.0/inputParams->Tcoh, exactCandidates1->data[ii].fsig+1.0/inputParams->Tcoh, 5, 5, exactCandidates1->data[ii].moddepth-1.0/inputParams->Tcoh, exactCandidates1->data[ii].moddepth+1.0/inputParams->Tcoh, 5, inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 1);
               //bruteForceTemplateSearch(&(exactCandidates2->data[exactCandidates2->numofcandidates]), exactCandidates1->data[ii], exactCandidates1->data[ii].fsig-1.0/inputParams->Tcoh, exactCandidates1->data[ii].fsig+1.0/inputParams->Tcoh, 5, 3, exactCandidates1->data[ii].moddepth-1.0/inputParams->Tcoh, exactCandidates1->data[ii].moddepth+1.0/inputParams->Tcoh, 5, inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 1);
               bruteForceTemplateSearch(&(exactCandidates2->data[exactCandidates2->numofcandidates]),
                                        exactCandidates1->data[ii],
                                        exactCandidates1->data[ii].fsig-0.5/inputParams->Tcoh,
                                        exactCandidates1->data[ii].fsig+0.5/inputParams->Tcoh,
                                        3,
                                        3,
                                        exactCandidates1->data[ii].moddepth-0.5/inputParams->Tcoh,
                                        exactCandidates1->data[ii].moddepth+0.5/inputParams->Tcoh,
                                        3,
                                        inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 1);
               if (xlalErrno!=0) {
                  fprintf(stderr, "%s: bruteForceTemplateSearch() failed.\n", __func__);
                  XLAL_ERROR(XLAL_EFUNC);
               }
            } else {
               //bruteForceTemplateSearch(&(exactCandidates2->data[exactCandidates2->numofcandidates]), exactCandidates1->data[ii], exactCandidates1->data[ii].fsig-1.0/inputParams->Tcoh, exactCandidates1->data[ii].fsig+1.0/inputParams->Tcoh, 5, 5, exactCandidates1->data[ii].moddepth-1.0/inputParams->Tcoh, exactCandidates1->data[ii].moddepth+1.0/inputParams->Tcoh, 5, inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 0);
               //bruteForceTemplateSearch(&(exactCandidates2->data[exactCandidates2->numofcandidates]), exactCandidates1->data[ii], exactCandidates1->data[ii].fsig-1.0/inputParams->Tcoh, exactCandidates1->data[ii].fsig+1.0/inputParams->Tcoh, 5, 3, exactCandidates1->data[ii].moddepth-1.0/inputParams->Tcoh, exactCandidates1->data[ii].moddepth+1.0/inputParams->Tcoh, 5, inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 0);
               bruteForceTemplateSearch(&(exactCandidates2->data[exactCandidates2->numofcandidates]),
                                        exactCandidates1->data[ii],
                                        exactCandidates1->data[ii].fsig-0.5/inputParams->Tcoh,
                                        exactCandidates1->data[ii].fsig+0.5/inputParams->Tcoh,
                                        3,
                                        3,
                                        exactCandidates1->data[ii].moddepth-0.5/inputParams->Tcoh,
                                        exactCandidates1->data[ii].moddepth+0.5/inputParams->Tcoh,
                                        3,
                                        inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 0);
               if (xlalErrno!=0) {
                  fprintf(stderr, "%s: bruteForceTemplateSearch() failed.\n", __func__);
                  XLAL_ERROR(XLAL_EFUNC);
               }
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
         
      } /* if IHSonly is not given */
      
      //Determine upper limits, if the ULoff has not been set
      if (!args_info.ULoff_given) {
         upperlimits->data[upperlimits->length-1].alpha = (REAL4)dopplerpos.Alpha;
         upperlimits->data[upperlimits->length-1].delta = (REAL4)dopplerpos.Delta;
         upperlimits->data[upperlimits->length-1].normalization = ffdata->tfnormalization;
         skypoint95UL(&(upperlimits->data[upperlimits->length-1]), inputParams, ffdata, ihsmaxima, ihsfarstruct, aveTFnoisePerFbinRatio);
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: skypoint95UL() failed.\n", __func__);
            XLAL_ERROR(XLAL_EFUNC);
         }
         for (ii=0; ii<(INT4)upperlimits->data[upperlimits->length-1].ULval->length; ii++) upperlimits->data[upperlimits->length-1].ULval->data[ii] /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25);  //Compensation for different duty cycle and antenna pattern weights
         upperlimits = resize_UpperLimitVector(upperlimits, upperlimits->length+1);
         if (upperlimits->data==NULL) {
            fprintf(stderr,"%s: resize_UpperLimitVector(%d) failed.\n", __func__, upperlimits->length+1);
            XLAL_ERROR(XLAL_EFUNC);
         }
      } /* if producing UL */
      
      //Destroy stuff
      XLALDestroyREAL4Vector(aveTFnoisePerFbinRatio);
      XLALDestroyREAL4VectorSequence(trackedlines);
      
      //Iterate to next sky location
      if ((XLALNextDopplerSkyPos(&dopplerpos, &scan))!=0) {
         fprintf(stderr,"%s: XLALNextDopplerSkyPos() failed.\n", __func__);
         XLAL_ERROR(XLAL_EFUNC);
      }
      
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
   if (!args_info.ULoff_given) {
      ULFILE = fopen(t,"w");
      if (ULFILE==NULL) {
         fprintf(stderr, "%s: UL file could not be opened.\n", __func__);
         XLAL_ERROR(XLAL_EINVAL);
      }
      for (ii=0; ii<(INT4)upperlimits->length-1; ii++) outputUpperLimitToFile(ULFILE, upperlimits->data[ii], inputParams->printAllULvalues);
      fclose(ULFILE);
   }
   
   //Destroy varaibles
   XLALDestroyREAL4Vector(antweightsforihs2h0);
   XLALDestroyREAL4Vector(background);
   XLALDestroyREAL4Vector(usableTFdata);
   XLALDestroyREAL4Vector(detectorVelocities);
   XLALDestroyREAL4Vector(aveNoise);
   XLALDestroyINT4Vector(binshifts);
   XLALDestroyINT4Vector(sftexist);
   XLALDestroyINT4Vector(indexValuesOfExistingSFTs);
   free_ffdata(ffdata);
   free_ihsfarStruct(ihsfarstruct);
   free_inputParams(inputParams);
   free_ihsMaxima(ihsmaxima);
   XLALDestroyREAL4FFTPlan(secondFFTplan);
   XLALFree((CHAR*)sft_dir);
   XLALFree((CHAR*)earth_ephemeris);
   XLALFree((CHAR*)sun_ephemeris);
   XLALFree((CHAR*)sky);
   XLALDestroyEphemerisData(edat);
   cmdline_parser_free(&args_info);
   XLALFree(configparams);
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



//////////////////////////////////////////////////////////////
// Create new inputParamsStruct  -- done
inputParamsStruct * new_inputParams(INT4 numofIFOs)
{
   
   inputParamsStruct *input = XLALMalloc(sizeof(*input));
   if (input==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%zu) failed.", __func__, sizeof(*input));
      XLAL_ERROR_NULL(XLAL_ENOMEM);
   }
   input->det = XLALMalloc(numofIFOs*sizeof(LALDetector));
   if (input->det==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%zu) failed.", __func__, numofIFOs*sizeof(LALDetector));
      XLAL_ERROR_NULL(XLAL_ENOMEM);
   }
   input->rng = gsl_rng_alloc(gsl_rng_mt19937);
   if (input->rng==NULL) {
      fprintf(stderr,"%s: gsl_rng_alloc() failed.", __func__);
      XLAL_ERROR_NULL(XLAL_ENOMEM);
   }
   
   return input;

} /* new_inputParams() */


//////////////////////////////////////////////////////////////
// Destroy inputParamsStruct  -- done
void free_inputParams(inputParamsStruct *input)
{
   
   XLALFree((CHAR*)input->sftType);
   XLALFree((LALDetector*)input->det);
   gsl_rng_free(input->rng);
   XLALFree((inputParamsStruct*)input);

} /* free_inputParams() */



//////////////////////////////////////////////////////////////
// Allocate ffdataStruct vectors  -- done
ffdataStruct * new_ffdata(inputParamsStruct *input)
{
   
   ffdataStruct *ffdata = XLALMalloc(sizeof(*ffdata));
   if (ffdata==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%zu) failed.\n", __func__, sizeof(*ffdata));
      XLAL_ERROR_NULL(XLAL_ENOMEM);
   }
   
   ffdata->numfbins = (INT4)(round(input->fspan*input->Tcoh + 2.0*input->dfmax*input->Tcoh)+12+1);
   ffdata->numffts = (INT4)floor(input->Tobs/(input->Tcoh-input->SFToverlap)-1);
   ffdata->numfprbins = (INT4)floorf(ffdata->numffts*0.5) + 1;
   
   ffdata->ffdata = XLALCreateREAL4Vector(ffdata->numfbins * ffdata->numfprbins);
   if (ffdata->ffdata==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ffdata->numfbins * ffdata->numfprbins);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   
   ffdata->tfnormalization = ffdata->ffnormalization = 0.0;
   
   return ffdata;
   
} /* new_ffdata() */


//////////////////////////////////////////////////////////////
// Destroy ffdataStruct and vectors  -- done
void free_ffdata(ffdataStruct *data)
{

   XLALDestroyREAL4Vector(data->ffdata);
   XLALFree((ffdataStruct*)data);

} /* free_ffdata() */


//////////////////////////////////////////////////////////////
// Read in SFT data to produce a TF vector
REAL4Vector * readInSFTs(inputParamsStruct *input, REAL8 *normalization)
{
   
   INT4 ii, jj;
   LALStatus status = empty_status;
   SFTCatalog *catalog = NULL;
   
   //Set the start and end times in the LIGO GPS format
   LIGOTimeGPS start = LIGOTIMEGPSZERO, end = LIGOTIMEGPSZERO;
   XLALGPSSetREAL8(&start, input->searchstarttime);
   if (xlalErrno != 0) {
      fprintf(stderr, "%s: XLALGPSSetREAL8() failed on start time = %.9f.\n", __func__, input->searchstarttime);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   XLALGPSSetREAL8(&end, input->searchstarttime+input->Tobs);
   if (xlalErrno != 0) {
      fprintf(stderr, "%s: XLALGPSSetREAL8() failed on end time = %.9f.\n", __func__, input->searchstarttime+input->Tobs);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   
   //Setup the constraints
   SFTConstraints constraints;
   constraints.detector = NULL;
   constraints.startTime = constraints.endTime = NULL;
   constraints.timestamps = NULL;
   constraints.detector = input->det[0].frDetector.prefix;
   constraints.startTime = &start;
   constraints.endTime = &end;
   
   //Find SFT files
   LALSFTdataFind(&status, &catalog, sft_dir, &constraints);
   if (status.statusCode != 0) {
      fprintf(stderr,"%s: LALSFTdataFind() failed with code = %d.\n", __func__, status.statusCode);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   
   //Determine band size (remember to get extra bins because of the running median and the bin shifts due to detector velocity)
   REAL8 minfbin = round(input->fmin*input->Tcoh - input->dfmax*input->Tcoh - 0.5*(input->blksize-1) - (REAL8)(input->maxbinshift) - 6.0)/input->Tcoh;
   REAL8 maxfbin = round((input->fmin + input->fspan)*input->Tcoh + input->dfmax*input->Tcoh + 0.5*(input->blksize-1) + (REAL8)(input->maxbinshift) + 6.0)/input->Tcoh;
   
   //Now extract the data
   SFTVector *sfts = XLALLoadSFTs(catalog, minfbin+0.1/1800.0, maxfbin-0.1/1800.0);
   if (sfts == NULL) {
      fprintf(stderr,"%s: XLALLoadSFTs() failed to load SFTs with given input parameters.\n", __func__);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   
   //Now put the power data into the TF plane, looping through each SFT
   //If an SFT doesn't exit, fill the TF pixels of the SFT with zeros
   INT4 numffts = (INT4)floor(input->Tobs/(input->Tcoh-input->SFToverlap)-1);
   INT4 sftlength;
   if (sfts->length == 0) sftlength = (INT4)round(maxfbin*input->Tcoh - minfbin*input->Tcoh + 1);
   else {
      sftlength = sfts->data->data->length;
      //Check the length is what we expect
      if (sftlength!=(INT4)round(maxfbin*input->Tcoh - minfbin*input->Tcoh + 1)) {
         fprintf(stderr, "%s: sftlength (%d) is not matching expected length (%d).\n", __func__, sftlength, (INT4)round(maxfbin*input->Tcoh - minfbin*input->Tcoh + 1));
         XLAL_ERROR_NULL(XLAL_EFPINEXCT);
      }
   }
   INT4 nonexistantsft = 0;
   REAL4Vector *tfdata = XLALCreateREAL4Vector(numffts*sftlength);
   if (tfdata==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numffts*sftlength);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   
   //Load the data into the output vector, roughly normalizing as we go along from the input value
   REAL8 sqrtnorm = sqrt(*(normalization));
   for (ii=0; ii<numffts; ii++) {
      
      SFTDescriptor *sftdescription = &(catalog->data[ii - nonexistantsft]);
      if (sftdescription->header.epoch.gpsSeconds == (INT4)round(ii*(input->Tcoh-input->SFToverlap)+input->searchstarttime)) {
         
         SFTtype *sft = &(sfts->data[ii - nonexistantsft]);
         for (jj=0; jj<sftlength; jj++) {
            COMPLEX8 sftcoeff = sft->data->data[jj];
            tfdata->data[ii*sftlength + jj] = (REAL4)((sqrtnorm*crealf(sftcoeff))*(sqrtnorm*crealf(sftcoeff)) + (sqrtnorm*cimagf(sftcoeff))*(sqrtnorm*cimagf(sftcoeff)));  //power, normalized
         }
         
      } else {
         for (jj=0; jj<sftlength; jj++) tfdata->data[ii*sftlength + jj] = 0.0;   //Set values to be zero
         nonexistantsft++;    //increment the nonexistantsft counter
      }
      
   } /* for ii < numffts */
   
   //Vladimir's code uses a different SFT normalization factor than MFD
   if (strcmp(input->sftType, "vladimir") == 0) {
      REAL4 vladimirfactor = (REAL4)(0.25*(8.0/3.0));
      for (ii=0; ii<(INT4)tfdata->length; ii++) tfdata->data[ii] *= vladimirfactor;
   }
   
   //Destroy stuff
   XLALDestroySFTCatalog(catalog);
   XLALDestroySFTVector(sfts);
   
   fprintf(stderr,"TF before weighting, mean subtraction: mean = %g, std. dev. = %g\n", calcMean(tfdata), calcStddev(tfdata));
   
   return tfdata;

} /* readInSFTs() */

//Same as readInSFTs(), but now from multiple interferometers.
//!!! Not used, untested !!!
REAL4VectorSequence * readInMultiSFTs(inputParamsStruct *input, REAL8 *normalization)
{
   
   INT4 ii, jj, kk;
   LALStatus status = empty_status;
   SFTCatalog *catalog = NULL;
   
   LIGOTimeGPS start = LIGOTIMEGPSZERO, end = LIGOTIMEGPSZERO;
   XLALGPSSetREAL8(&start, input->searchstarttime);
   if (xlalErrno != 0) {
      fprintf(stderr, "%s: XLALGPSSetREAL8() failed on start time = %.9f.\n", __func__, input->searchstarttime);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   XLALGPSSetREAL8(&end, input->searchstarttime+input->Tobs);
   if (xlalErrno != 0) {
      fprintf(stderr, "%s: XLALGPSSetREAL8() failed on end time = %.9f.\n", __func__, input->searchstarttime+input->Tobs);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   
   SFTConstraints constraints;
   constraints.detector = NULL;
   constraints.startTime = constraints.endTime = NULL;
   constraints.timestamps = NULL;
   constraints.startTime = &start;
   constraints.endTime = &end;
   
   //Find SFT files
   LALSFTdataFind(&status, &catalog, sft_dir, &constraints);
   if (status.statusCode != 0) {
      fprintf(stderr,"%s: LALSFTdataFind() failed with code = %d.\n", __func__, status.statusCode);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   
   //Determine band size (remember to get extra bins because of the running median and the bin shifts due to detector velocity)
   REAL8 minfbin = round(input->fmin*input->Tcoh - 0.5*(input->blksize-1) - (input->maxbinshift))/input->Tcoh;
   REAL8 maxfbin = round((input->fmin + input->fspan)*input->Tcoh + 0.5*(input->blksize-1) + (input->maxbinshift))/input->Tcoh;
   
   //Now extract the data
   MultiSFTVector *sfts = XLALLoadMultiSFTs(catalog, minfbin, maxfbin);
   if (sfts == NULL) {
      fprintf(stderr,"%s: XLALLoadSFTs() failed to load SFTs with given input parameters.\n", __func__);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   
   //Now put the power data into the TF plane, looping through each SFT
   //If an SFT doesn't exit, fill the TF pixels of the SFT with zeros
   INT4 numffts = (INT4)floor(input->Tobs/(input->Tcoh-input->SFToverlap)-1);
   INT4 sftlength;
   if (sfts->length == 0) sftlength = (INT4)(maxfbin*input->Tcoh - minfbin*input->Tcoh + 1);
   else sftlength = sfts->data[0]->data->data->length;
   
   REAL4VectorSequence *multiTFdata = XLALCreateREAL4VectorSequence(input->numofIFOs, (numffts*sftlength));
   if (multiTFdata==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4VectorSequence(%d,%d) failed.\n", __func__, input->numofIFOs, numffts*sftlength);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   
   REAL8 sqrtnorm = sqrt(*(normalization));
   INT4 nonexistantsft = 0;
   INT4Vector *IFOspecificNonexistantsft = XLALCreateINT4Vector(input->numofIFOs);
   if (IFOspecificNonexistantsft==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, input->numofIFOs);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   for (ii=0; ii<input->numofIFOs; ii++) IFOspecificNonexistantsft->data[ii] = 0;
   
   for (ii=0; ii<numffts; ii++) {
      for (jj=0; jj<input->numofIFOs; jj++) {
         SFTDescriptor *sftdescription = &(catalog->data[ii*input->numofIFOs+jj-nonexistantsft]);
         if (sftdescription->header.epoch.gpsSeconds == (INT4)round(ii*(input->Tcoh-input->SFToverlap)+input->searchstarttime)) {
            SFTtype *sft = &(sfts->data[jj]->data[ii-IFOspecificNonexistantsft->data[jj]]);
            for (kk=0; kk<sftlength; kk++) {
               COMPLEX8 sftcoeff = sft->data->data[kk];
               multiTFdata->data[jj*multiTFdata->vectorLength + ii*sftlength + kk] = (REAL4)((sqrtnorm*crealf(sftcoeff))*(sqrtnorm*crealf(sftcoeff)) + (sqrtnorm*cimagf(sftcoeff))*(sqrtnorm*cimagf(sftcoeff)));  //power, normalized
            }
         } else {
            for (kk=0; kk<sftlength; kk++) multiTFdata->data[jj*multiTFdata->vectorLength + ii*sftlength + kk] = 0.0;   //Set values to be zero
            nonexistantsft++;    //increment the nonexistantsft counter
            IFOspecificNonexistantsft->data[jj]++;
         }
      }
   }
   
   //Vladimir's code uses a different SFT normalization factor than MFD
   if (strcmp(input->sftType, "vladimir") == 0) {
      REAL4 vladimirfactor = (REAL4)(0.25*(8.0/3.0));
      for (ii=0; ii<(INT4)(multiTFdata->length*multiTFdata->vectorLength); ii++) multiTFdata->data[ii] *= vladimirfactor;
   }
   
   XLALDestroySFTCatalog(catalog);
   LALDestroyMultiSFTVector(&status, &sfts);
   XLALDestroyINT4Vector(IFOspecificNonexistantsft);
   
   return multiTFdata;
   
} /* readInMultiSFTs() */



//////////////////////////////////////////////////////////////
// Slide SFT TF data
void slideTFdata(REAL4Vector *output, inputParamsStruct *input, REAL4Vector *tfdata, INT4Vector *binshifts)
{
   
   INT4 ii, jj;
   INT4 numffts = (INT4)floor(input->Tobs/(input->Tcoh-input->SFToverlap)-1);
   INT4 numfbins = (INT4)(round(input->fspan*input->Tcoh+2.0*input->dfmax*input->Tcoh)+12+1);
   
   for (ii=0; ii<numffts; ii++) {
      if (binshifts->data[ii]>input->maxbinshift) {
         fprintf(stderr, "%s: SFT slide value %d is greater than maximum value predicted (%d)", __func__, binshifts->data[ii], input->maxbinshift);
         XLAL_ERROR_VOID(XLAL_EFAILED);
      }
      for (jj=0; jj<numfbins; jj++) output->data[ii*numfbins + jj] = tfdata->data[ii*(numfbins+2*input->maxbinshift) + jj + input->maxbinshift + binshifts->data[ii]];
   }
   
} /* slideTFdata() */




//////////////////////////////////////////////////////////////
// Determine the TF running mean of each SFT  -- 
// numffts = number of ffts
// numfbins = number of fbins in the search + 2*maximum bin shift
// blksize = running median blocksize
void tfRngMeans(REAL4Vector *output, REAL4Vector *tfdata, INT4 numffts, INT4 numfbins, INT4 blksize)
{
   
   LALStatus status = empty_status;
   REAL8 bias;
   INT4 ii, jj;
   INT4 totalfbins = numfbins + blksize - 1;
   
   //Blocksize of running median
   LALRunningMedianPar block = {blksize};
   
   //Running median bias calculation
   if (blksize<1000) {
      bias = XLALRngMedBias(blksize);
      if (xlalErrno!=0) {
	 fprintf(stderr,"%s: XLALRngMedBias(%d) failed.\n", __func__, blksize);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
   } else  bias = LAL_LN2;
   REAL8 invbias = 1.0/(bias*1.0099993480677538);  //StackSlide normalization for 101 bins
   
   //Allocate for a single SFT data and the medians out of each SFT
   REAL4Vector *inpsd = XLALCreateREAL4Vector(totalfbins);
   REAL4Vector *mediansout = XLALCreateREAL4Vector(numfbins);
   if (inpsd==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, totalfbins);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (mediansout==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numfbins);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }

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
         if (status.statusCode != 0) {
            fprintf(stderr,"%s: LALSRunningMedian2() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         
         //Now make the output medians into means by multiplying by 1/bias
         for (jj=0; jj<(INT4)mediansout->length; jj++) output->data[ii*numfbins + jj] = (REAL4)(mediansout->data[jj]*invbias);
      } else {
         //Otherwise, set means to zero
         for (jj=0; jj<(INT4)mediansout->length; jj++) output->data[ii*numfbins + jj] = 0.0;
         //fprintf(SFTMEANVALS, "%g\n", 0.0);
      }
   } /* for ii < numffts */
   
   //fclose(SFTMEANVALS);
   
   fprintf(stderr,"Mean of running means = %g\n",calcMean(output));
   
   //Destroy stuff
   XLALDestroyREAL4Vector(inpsd);
   XLALDestroyREAL4Vector(mediansout);

} /* tfRngMeans() */
void multiTFRngMeans(REAL4VectorSequence *output, REAL4VectorSequence *multiTFdata, INT4 numffts, INT4 numfbins, INT4 blksize)
{
   
   LALStatus status = empty_status;
   REAL8 bias;
   INT4 ii, jj, kk;
   INT4 totalfbins = numfbins + blksize - 1;
   
   //Blocksize of running median
   LALRunningMedianPar block = {blksize};
   
   //Running median bias calculation
   if (blksize<1000) {
      LALRngMedBias(&status, &bias, blksize);
      if (status.statusCode != 0) {
         fprintf(stderr,"%s: LALRngMedBias() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
   } else {
      bias = LAL_LN2;
   }
   REAL8 invbias = 1.0/(bias*1.0099993480677538);  //StackSlide normalization for 101 bins
   
   REAL4Vector *inpsd = XLALCreateREAL4Vector(totalfbins);
   REAL4Vector *mediansout = XLALCreateREAL4Vector(numfbins);
   if (inpsd==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, totalfbins);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (mediansout==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numfbins);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)multiTFdata->length; ii++) {
      for (jj=0; jj<numffts; jj++) {
         if (multiTFdata->data[ii*multiTFdata->vectorLength + jj*totalfbins]!=0.0) {
            memcpy(inpsd->data, &(multiTFdata->data[ii*multiTFdata->vectorLength + jj*totalfbins]), sizeof(REAL4)*inpsd->length);
            LALSRunningMedian2(&status, mediansout, inpsd, block);
            if (status.statusCode != 0) {
               fprintf(stderr,"%s: LALSRunningMedian2() failed.\n", __func__);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
            for (kk=0; kk<(INT4)mediansout->length; kk++) output->data[ii*output->vectorLength + jj*mediansout->length + kk] = (REAL4)(mediansout->data[kk]*invbias);
         } else {
            for (kk=0; kk<(INT4)mediansout->length; kk++) output->data[ii*output->vectorLength + jj*mediansout->length + kk] = 0.0;
         }
      }
   }
   
   XLALDestroyREAL4Vector(inpsd);
   XLALDestroyREAL4Vector(mediansout);
   
} /* multiTFRngMeans() */
REAL4Vector * combineMultiTFrngMeans(REAL4VectorSequence *input, INT4 numffts, INT4 numfbins)
{
   
   INT4 ii, jj, kk;
   
   REAL4Vector *output = XLALCreateREAL4Vector(input->vectorLength);
   if (output==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, input->vectorLength);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)input->vectorLength; ii++) output->data[ii] = 0.0;
   
   REAL4Vector *singlepsd = XLALCreateREAL4Vector(numfbins);
   if (singlepsd==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numfbins);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   for (ii=0; ii<numffts; ii++) {
      REAL4 totalweightval = 0.0;
      for (jj=0; jj<(INT4)input->length; jj++) {
         if (input->data[jj*input->vectorLength+ii*numfbins]!=0.0) {
            memcpy(singlepsd->data, &(input->data[jj*input->vectorLength+ii*numfbins]), sizeof(REAL4)*numfbins);
            REAL4 meanval = calcMean(singlepsd);
            REAL4 weightval = 1.0/(meanval*meanval);
            totalweightval += weightval;
            for (kk=0; kk<numfbins; kk++) output->data[ii*numfbins + kk] += singlepsd->data[kk]*weightval;
         }
      }
      REAL4 invtotalweightval = 1.0/totalweightval;
      for (jj=0; jj<numfbins; jj++) output->data[ii*numfbins + jj] *= invtotalweightval;
   }
   
   XLALDestroyREAL4Vector(singlepsd);
   
   return output;
   
} /* combineMultiTFrngMeans() */


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
INT4Vector * markBadSFTs(REAL4Vector *tfdata, inputParamsStruct *params)
{
   
   INT4 ii;
   
   INT4 numffts = (INT4)floor(params->Tobs/(params->Tcoh-params->SFToverlap)-1);    //Number of FFTs
   INT4 numfbins = (INT4)(round(params->fspan*params->Tcoh+2.0*params->dfmax*params->Tcoh)+12+1)+2*params->maxbinshift+params->blksize-1;     //Number of frequency bins
   
   //Allocate output data vector and a single SFT data vector
   INT4Vector *output = XLALCreateINT4Vector(numffts);
   if (output==NULL) {
      fprintf(stderr, "%s: XLALCreateINT4Vector(%d) failed.\n", __func__, numffts);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   memset(output->data, 0, sizeof(INT4)*output->length);
   REAL4Vector *tempvect = XLALCreateREAL4Vector(numfbins);
   if (tempvect==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numfbins);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   
   //Do the KS and Kuiper test on each SFT
   REAL8 ksthreshold = 1.358/(sqrt(numfbins)+0.12+0.11/sqrt(numfbins));
   //REAL8 ksthreshold = 1.224/(sqrt(numfbins)+0.12+0.11/sqrt(numfbins));  //This is a tighter restriction
   //REAL8 ksthreshold = 1.073/(sqrt(numfbins)+0.12+0.11/sqrt(numfbins));  //This is an even tighter restriction
   REAL8 kuiperthreshold = 1.747/(sqrt(numfbins)+0.155+0.24/sqrt(numfbins));
   //REAL8 kuiperthreshold = 1.620/(sqrt(numfbins)+0.155+0.24/sqrt(numfbins));  //This is a tighter restriction
   //REAL8 kuiperthreshold = 1.473/(sqrt(numfbins)+0.155+0.24/sqrt(numfbins));  //This is an even tighter restriction
   for (ii=0; ii<numffts; ii++) {
      if (tfdata->data[ii*numfbins]!=0.0) {
         memcpy(tempvect->data, &(tfdata->data[ii*numfbins]), sizeof(REAL4)*tempvect->length);
         REAL8 kstest = ks_test_exp(tempvect);
         if (XLAL_IS_REAL8_FAIL_NAN(kstest)) {
            fprintf(stderr,"%s: ks_test_exp() failed.\n", __func__);
            XLAL_ERROR_NULL(XLAL_EFUNC);
         }
         
         REAL8 kuipertest = kuipers_test_exp(tempvect);
         if (XLAL_IS_REAL8_FAIL_NAN(kuipertest)) {
            fprintf(stderr,"%s: kuipers_test_exp() failed.\n", __func__);
            XLAL_ERROR_NULL(XLAL_EFUNC);
         }
         
         if (kstest>ksthreshold || kuipertest>kuiperthreshold) output->data[ii] = 1;
      }
   }
   
   //Destroy stuff
   XLALDestroyREAL4Vector(tempvect);
   
   return output;
   
}
INT4VectorSequence * markBadMultiSFTs(REAL4VectorSequence *multiTFdata, inputParamsStruct *params)
{
   
   INT4 ii, jj;
   
   INT4 numffts = (INT4)floor(params->Tobs/(params->Tcoh-params->SFToverlap)-1);    //Number of FFTs
   INT4 numfbins = (INT4)(round(params->fspan*params->Tcoh)+1)+2*params->maxbinshift+params->blksize-1;     //Number of frequency bins
   
   INT4VectorSequence *output = XLALCreateINT4VectorSequence(params->numofIFOs, numffts);
   if (output==NULL) {
      fprintf(stderr, "%s: XLALCreateINT4Vector(%d) failed.\n", __func__, numffts);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)(output->length*output->vectorLength); ii++) output->data[ii] = 0;
   REAL4Vector *tempvect = XLALCreateREAL4Vector(numfbins);
   if (tempvect==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numfbins);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   
   REAL8 ksthreshold = 1.358/(sqrt(numfbins)+0.12+0.11/sqrt(numfbins));
   for (ii=0; ii<params->numofIFOs; ii++) {
      for (jj=0; jj<numffts; jj++) {
         if (multiTFdata->data[ii*multiTFdata->vectorLength+jj*numfbins]!=0.0) {
            memcpy(tempvect->data, &(multiTFdata->data[ii*multiTFdata->vectorLength+jj*numfbins]), sizeof(REAL4)*tempvect->length);
            REAL8 kstest = ks_test_exp(tempvect);
            if (XLAL_IS_REAL8_FAIL_NAN(kstest)) {
               fprintf(stderr,"%s: ks_test_exp() failed.\n", __func__);
               XLAL_ERROR_NULL(XLAL_EFUNC);
            }
            
            if (kstest>ksthreshold) output->data[ii*numffts+jj] = 1;
         }
      }
   }
   
   XLALDestroyREAL4Vector(tempvect);
   
   return output;
   
}

//Those SFTs which have been marked as bad, set them to 0
//This modifies the tfdata vector!!!
void removeBadSFTs(REAL4Vector *tfdata, INT4Vector *badsfts)
{
   
   INT4 ii;
   
   INT4 numfbins_tfdata = tfdata->length/badsfts->length;
   
   for (ii=0; ii<(INT4)badsfts->length; ii++) if (badsfts->data[ii]==1) memset(&(tfdata->data[ii*numfbins_tfdata]), 0, sizeof(REAL4)*numfbins_tfdata);
   
}

//Remove the SFTs (set to zero) the SFTs which fail the K-S and Kuiper's tests
void removeBadMultiSFTs(REAL4VectorSequence *multiTFdata, INT4VectorSequence *badsfts)
{
   
   INT4 ii, jj, kk;
   
   INT4 numfbins_tfdata = multiTFdata->vectorLength/badsfts->vectorLength;
   
   for (ii=0; ii<(INT4)badsfts->length; ii++) {
      for (jj=0; jj<(INT4)badsfts->vectorLength; jj++) {
         if (badsfts->data[ii*badsfts->vectorLength+jj]==1) {
            for (kk=0; kk<numfbins_tfdata; kk++) multiTFdata->data[ii*multiTFdata->vectorLength+jj*numfbins_tfdata+kk] = 0.0;
         }
      }
   }
   
}


//Running median based line detection
// 1. Calculate RMS for each fbin as function of time
// 2. Calculate running median of RMS values
// 3. Divide RMS values by running median. Lines stick out by being >>1
// 4. Add up number of lines above threshold and report
INT4Vector * detectLines_simple(REAL4Vector *TFdata, ffdataStruct *ffdata, inputParamsStruct *params)
{
   
   LALStatus status = empty_status;
   
   INT4 blksize = 11, ii, jj;
   
   INT4 numlines = 0;
   INT4Vector *lines = NULL;
   INT4 totalnumfbins = ffdata->numfbins+(params->blksize-1)+2*params->maxbinshift;
   
   //Blocksize of running median
   LALRunningMedianPar block = {blksize};
   
   //Compute weights
   REAL4Vector *sftdata = XLALCreateREAL4Vector(totalnumfbins);
   REAL4Vector *weights = XLALCreateREAL4Vector(ffdata->numffts);
   if (sftdata==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, totalnumfbins);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   } else if (weights==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ffdata->numffts);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   memset(weights->data, 0, ffdata->numffts*sizeof(REAL4));
   REAL4 sumweights = 0.0;
   for (ii=0; ii<ffdata->numffts; ii++) {
      if (TFdata->data[ii*totalnumfbins]!=0.0) {
         memcpy(sftdata->data, &(TFdata->data[ii*totalnumfbins]), totalnumfbins*sizeof(REAL4));
         REAL4 stddev = calcStddev(sftdata);
         weights->data[ii] = 1.0/(stddev*stddev);
         sumweights += weights->data[ii];
      }
   }
   REAL4 invsumweights = 1.0/sumweights;
   XLALDestroyREAL4Vector(sftdata);
   
   //Compute RMS for each frequency bin as a function of time
   REAL4Vector *testRMSvals = XLALCreateREAL4Vector(totalnumfbins);
   REAL4Vector *testRngMedian = XLALCreateREAL4Vector(totalnumfbins-(blksize-1));
   REAL4Vector *testTSofPowers = XLALCreateREAL4Vector(ffdata->numffts);
   if (testRMSvals==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, totalnumfbins);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   } else if (testRngMedian==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, totalnumfbins-(blksize-1));
      XLAL_ERROR_NULL(XLAL_EFUNC);
   } else if (testTSofPowers==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ffdata->numffts);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)testRMSvals->length; ii++) {
      for (jj=0; jj<ffdata->numffts; jj++) testTSofPowers->data[jj] = TFdata->data[jj*testRMSvals->length + ii]*weights->data[jj]*invsumweights;
      testRMSvals->data[ii] = calcRms(testTSofPowers); //This approaches calcMean(TSofPowers) for stationary noise
   }
   
   //Running median of RMS values
   LALSRunningMedian2(&status, testRngMedian, testRMSvals, block);
   if (status.statusCode != 0) {
      fprintf(stderr,"%s: LALSRunningMedian2() failed.\n", __func__);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   
   //Determine the mid frequency, low frequency and high frequency of the line when considering SFT shifts
   REAL4 f0 = (REAL4)(round(params->fmin*params->Tcoh - params->dfmax*params->Tcoh - 6.0 - 0.5*(params->blksize-1) - (REAL8)(params->maxbinshift) + 0.5*(blksize-1))/params->Tcoh);
   REAL4 df = 1.0/params->Tcoh;
   for (ii=0; ii<(INT4)testRngMedian->length; ii++) {
      REAL4 normrmsval = testRMSvals->data[ii+(blksize-1)/2]/testRngMedian->data[ii];
      if ( (ii+(blksize-1)/2) > ((params->blksize-1)/2) && normrmsval > params->lineDetection) {
         lines = XLALResizeINT4Vector(lines, numlines+1);
         if (lines==NULL) {
            fprintf(stderr,"%s: XLALResizeINT4Vector(lines,%d) failed.\n", __func__, numlines+1);
            XLAL_ERROR_NULL(XLAL_EFUNC);
         }
         lines->data[numlines] = ii+(blksize-1)/2;
         numlines++;
      }
      
      if (NORMRMSOUT!=NULL) fprintf(NORMRMSOUT, "%f %.6f\n", f0+ii*df, normrmsval);
   }
   
   //Destroy stuff
   XLALDestroyREAL4Vector(testTSofPowers);
   XLALDestroyREAL4Vector(testRngMedian);
   XLALDestroyREAL4Vector(testRMSvals);
   XLALDestroyREAL4Vector(weights);
   
   return lines;
   
}

//Track the lines as the SFTs are shifted
REAL4VectorSequence * trackLines(INT4Vector *lines, INT4Vector *binshifts, inputParamsStruct *params)
{
   
   REAL4VectorSequence *output = XLALCreateREAL4VectorSequence(lines->length, 3);
   if (output==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4VectorSequence(%d,3) failed.\n", __func__, lines->length);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   
   REAL4 df = 1.0/params->Tcoh;
   REAL4 minfbin = (REAL4)(round(params->fmin*params->Tcoh - params->dfmax*params->Tcoh - 6.0 - 0.5*(params->blksize-1) - (REAL8)(params->maxbinshift))/params->Tcoh);
   
   INT4 maxshiftindex = 0, minshiftindex = 0;
   min_max_index_INT4Vector(binshifts, &minshiftindex, &maxshiftindex);
   INT4 maxshift = binshifts->data[maxshiftindex], minshift = binshifts->data[minshiftindex];
   
   INT4 ii;
   for (ii=0; ii<(INT4)lines->length; ii++) {
      output->data[ii*3] = lines->data[ii]*df + minfbin;
      //output->data[ii*3 + 1] = (lines->data[ii] + minshift)*df + minfbin;
      //output->data[ii*3 + 2] = (lines->data[ii] + maxshift)*df + minfbin;
      output->data[ii*3 + 1] = (lines->data[ii] + (minshift-1))*df + minfbin;  //Add one extra bin for buffer
      output->data[ii*3 + 2] = (lines->data[ii] + (maxshift+1))*df + minfbin;  //Add one extra bin for buffer
   }
   
   return output;
   
}

//Output a vector of zeros and ones; zero if SFT is missing, 1 if SFT is present
INT4Vector * existingSFTs(REAL4Vector *tfdata, inputParamsStruct *params, INT4 numfbins, INT4 numffts)
{
   
   INT4 ii;
   
   INT4Vector *sftexist = XLALCreateINT4Vector(numffts);
   if (sftexist==NULL) {
      fprintf(stderr, "\n%s: XLALCreateINT4Vector(%d) failed.\n", __func__, numffts);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   
   for (ii=0; ii<numffts; ii++) {
      if (tfdata->data[ii*(numfbins+2*params->maxbinshift+params->blksize-1)] == 0.0) sftexist->data[ii] = 0;
      else {
         sftexist->data[ii] = 1;
      }
   }
   
   return sftexist;
   
} /* existingSFTs() */

//Untested (should not be used) function to list existing SFTs from multiple detectors
INT4VectorSequence * existingMultiSFTs(REAL4VectorSequence *tfdata, inputParamsStruct *params, INT4 numfbins, INT4 numffts)
{
   
   INT4 ii, jj;
   
   INT4VectorSequence *sftexist = XLALCreateINT4VectorSequence(params->numofIFOs, numffts);
   if (sftexist==NULL) {
      fprintf(stderr, "\n%s: XLALCreateINT4VectorSequence(%d,%d) failed.\n", __func__, params->numofIFOs, numffts);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   for (ii=0; ii<params->numofIFOs; ii++) {
      for (jj=0; jj<numffts; jj++) {
         if (tfdata->data[ii*tfdata->vectorLength + jj*(numfbins+2*params->maxbinshift+params->blksize-1)]==0.0) sftexist->data[ii*sftexist->vectorLength + jj] = 0;
         else sftexist->data[ii*sftexist->vectorLength + jj] = 1;
      }
   }
   
   return sftexist;
   
} /* existingMultiSFTs() */

//Untested (should not be used) function to combine the existing SFTs from multiple detectors
INT4Vector * combineExistingMultiSFTs(INT4VectorSequence *input)
{
   
   INT4 ii, jj;
   INT4Vector *sftexist = XLALCreateINT4Vector(input->vectorLength);
   if (sftexist==NULL) {
      fprintf(stderr, "\n%s: XLALCreateINT4Vector(%d) failed.\n", __func__, input->vectorLength);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)sftexist->length; ii++) sftexist->data[ii] = 0;
   for (ii=0; ii<(INT4)input->vectorLength; ii++) {
      jj = 0;
      while (jj<(INT4)input->length && sftexist->data[ii]==0) {
         if (input->data[jj*input->vectorLength + ii]==1) sftexist->data[ii] = 1;
         jj++;
      }
   }
   
   return sftexist;
   
} /* combineExistingMultiSFTs() */


/* Modifies input vector!!! */
//Subtract the mean of the SFTs from the SFT data
void tfMeanSubtract(REAL4Vector *tfdata, REAL4Vector *rngMeans, INT4 numffts, INT4 numfbins)
{
   
   INT4 ii, jj;
   
   for (ii=0; ii<numffts; ii++) if (rngMeans->data[ii*numfbins]!=0.0) for (jj=0; jj<numfbins; jj++) tfdata->data[ii*numfbins+jj] -= rngMeans->data[ii*numfbins+jj];
   
} /* tfMeanSubtract() */


//Weight the SFTs from the antenna pattern weights and the variance of the SFT powers
void tfWeight(REAL4Vector *output, REAL4Vector *tfdata, REAL4Vector *rngMeans, REAL4Vector *antPatternWeights, INT4Vector *indexValuesOfExistingSFTs, inputParamsStruct *input)
{
   
   INT4 ii, jj;
   
   INT4 numffts = (INT4)antPatternWeights->length;
   INT4 numfbins = (INT4)(tfdata->length)/numffts;
   
   REAL4Vector *antweightssq = XLALCreateREAL4Vector(numffts);
   REAL4Vector *rngMeanssq = XLALCreateREAL4Vector(numffts);
   if (antweightssq==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numffts);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (rngMeanssq==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numffts);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   
   //memset(output->data, 0, sizeof(REAL4)*output->length);
   //memset(rngMeanssq->data, 0, sizeof(REAL4)*rngMeanssq->length);
   
   //User specifies whether to use SSE to do the multiplication or not
   if (input->useSSE) {
      antweightssq = sseSSVectorMultiply(antweightssq, antPatternWeights, antPatternWeights);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: sseSSVectorMultiply() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
   } else {
      //for (ii=0; ii<numffts; ii++) antweightssq->data[ii] = antPatternWeights->data[ii]*antPatternWeights->data[ii];
      antweightssq = XLALSSVectorMultiply(antweightssq, antPatternWeights, antPatternWeights);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: XLALSSVectorMultiply() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
   }
   
   //Loop through the SFT frequency bins and weight the data
   for (ii=0; ii<numfbins; ii++) {
      
      rngMeanssq = fastSSVectorMultiply_with_stride_and_offset(rngMeanssq, rngMeans, rngMeans, numfbins, numfbins, ii, ii);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: SSVectorMutiply_with_stride_and_offset() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      
      //If noiseWeightOff is given, then set all the noise weights to be 1.0
      if (input->noiseWeightOff!=0) for (jj=0; jj<(INT4)rngMeanssq->length; jj++) if (rngMeanssq->data[jj]!=0.0) rngMeanssq->data[jj] = 1.0;
      
      //Get sum of antenna pattern weight/variances for each frequency bin as a function of time (only for existant SFTs)
      //REAL8 sumofweights = 0.0;
      //for (jj=0; jj<numffts; jj++) if (rngMeanssq->data[jj] != 0.0) sumofweights += antweightssq->data[jj]/rngMeanssq->data[jj];
      REAL8 sumofweights = determineSumOfWeights(antweightssq, rngMeanssq);
      REAL8 invsumofweights = 1.0/sumofweights;
      
      //Now do noise weighting, antenna pattern weighting
      /* for (jj=0; jj<numffts; jj++) {
         if (rngMeanssq->data[jj] != 0.0) output->data[jj*numfbins+ii] = (REAL4)(invsumofweights*antPatternWeights->data[jj]*tfdata->data[jj*numfbins+ii]/rngMeanssq->data[jj]);
      } */ /* for jj < numffts */
      for (jj=0; jj<(INT4)indexValuesOfExistingSFTs->length; jj++) {
         output->data[indexValuesOfExistingSFTs->data[jj]*numfbins+ii] = (REAL4)(invsumofweights*antPatternWeights->data[indexValuesOfExistingSFTs->data[jj]]*tfdata->data[indexValuesOfExistingSFTs->data[jj]*numfbins+ii]/rngMeanssq->data[indexValuesOfExistingSFTs->data[jj]]);
      }
   } /* for ii < numfbins */
   
   //Destroy stuff
   XLALDestroyREAL4Vector(antweightssq);
   XLALDestroyREAL4Vector(rngMeanssq);
   
   //fprintf(stderr,"TF after weighting, mean subtraction = %g\n",calcMean(output));   
   
} /* tfWeight() */



//////////////////////////////////////////////////////////////
// Do the weighting by noise variance (from tfRngMeans), mean subtraction, and antenna pattern weights  -- done
void tfWeightMeanSubtract(REAL4Vector *output, REAL4Vector *tfdata, REAL4Vector *rngMeans, REAL4Vector *antPatternWeights, inputParamsStruct *input)
{
   
   INT4 ii, jj;
    
   INT4 numffts = (INT4)floor(input->Tobs/(input->Tcoh-input->SFToverlap)-1);    //Number of FFTs
   INT4 numfbins = (INT4)(round(input->fspan*input->Tcoh*2.0*input->dfmax*input->Tcoh)+12+1);                    //Number of frequency bins
   
   REAL4Vector *antweightssq = XLALCreateREAL4Vector(numffts);
   REAL4Vector *rngMeanssq = XLALCreateREAL4Vector(numffts);
   if (antweightssq==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numffts);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (rngMeanssq==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numffts);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   
   for (ii=0; ii<(INT4)rngMeanssq->length; ii++) rngMeanssq->data[ii] = 0.0;
   
   //for (ii=0; ii<numffts; ii++) antweightssq->data[ii] = antPatternWeights->data[ii]*antPatternWeights->data[ii];
   //antweightssq = XLALSSVectorMultiply(antweightssq, antPatternWeights, antPatternWeights);
   antweightssq = fastSSVectorMultiply_with_stride_and_offset(antweightssq, antPatternWeights, antPatternWeights, 1, 1, 0, 0);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: SSVectorMutiply_with_stride_and_offset() failed.\n", __func__);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   
   for (ii=0; ii<numfbins; ii++) {
      
      rngMeanssq = fastSSVectorMultiply_with_stride_and_offset(rngMeanssq, rngMeans, rngMeans, numfbins, numfbins, ii, ii);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: SSVectorMutiply_with_stride_and_offset() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      
      //If noiseWeightOff is given, then set all the noise weights to be 1.0
      if (input->noiseWeightOff!=0) {
         for (jj=0; jj<(INT4)rngMeanssq->length; jj++) {
            if (rngMeanssq->data[jj]!=0.0) {
               rngMeanssq->data[jj] = 1.0;
            }
         }
      }
      
      //Get sum of antenna pattern weight/variances for each frequency bin as a function of time (only for existant SFTs)
      REAL8 sumofweights = 0.0;
      for (jj=0; jj<numffts; jj++) if (rngMeanssq->data[jj] != 0.0) sumofweights += antweightssq->data[jj]/rngMeanssq->data[jj];
      REAL8 invsumofweights = 1.0/sumofweights;
      
      //Now do mean subtraction, noise weighting, antenna pattern weighting
      for (jj=0; jj<numffts; jj++) {
         if (rngMeanssq->data[jj] != 0.0) {
            output->data[jj*numfbins+ii] = (REAL4)(invsumofweights*antPatternWeights->data[jj]*(tfdata->data[jj*numfbins+ii] - rngMeans->data[jj*numfbins+ii])/rngMeanssq->data[jj]);
         } else {
            output->data[jj*numfbins+ii] = 0.0;
         }
      } /* for jj < numffts */
   } /* for ii < numfbins */
   
   //Destroy stuff
   XLALDestroyREAL4Vector(antweightssq);
   XLALDestroyREAL4Vector(rngMeanssq);
   
   //fprintf(stderr,"TF after weighting, mean subtraction = %g\n",calcMean(output));

} /* tfWeightMeanSubtract() */


//Determine the sum of the weights
REAL8 determineSumOfWeights(REAL4Vector *antweightssq, REAL4Vector *rngMeanssq)
{
   
   INT4 ii;
   REAL8 sumofweights = 0.0;
   for (ii=0; ii<(INT4)antweightssq->length; ii++) if (rngMeanssq->data[ii] != 0.0) sumofweights += antweightssq->data[ii]/rngMeanssq->data[ii];
   
   return sumofweights;
   
} /* determineSumOfWeights */


//////////////////////////////////////////////////////////////
// Make the second FFT powers
void makeSecondFFT(ffdataStruct *output, REAL4Vector *tfdata, REAL4FFTPlan *plan)
{
   
   INT4 ii, jj;
   REAL8 winFactor = 8.0/3.0;
   REAL8 psdfactor = winFactor;
   
   //Do the second FFT
   REAL4Vector *x = XLALCreateREAL4Vector(output->numffts);
   if (x==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, output->numffts);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   REAL4Window *win = XLALCreateHannREAL4Window(x->length);
   REAL4Vector *psd = XLALCreateREAL4Vector((UINT4)floor(x->length*0.5)+1);
   if (win==NULL) {
      fprintf(stderr,"%s: XLALCreateHannREAL4Window(%d) failed.\n", __func__, x->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (psd==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, (UINT4)floor(x->length*0.5)+1);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   
   for (ii=0; ii<output->numfbins; ii++) {
      
      //Next, loop over times and pick the right frequency bin for each FFT and window
      //for (jj=0; jj<(INT4)x->length; jj++) x->data[jj] = (tfdata->data[ii + jj*numfbins]*win->data->data[jj]);
      x = fastSSVectorMultiply_with_stride_and_offset(x, tfdata, win->data, output->numfbins, 1, ii, 0);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: SSVectorMutiply_with_stride_and_offset() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      
      //Make the FFT
      if ( (XLALREAL4PowerSpectrum(psd, x, plan)) != 0) {
         fprintf(stderr,"%s: XLALREAL4PowerSpectrum() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      
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
      for (jj=0; jj<(INT4)psd->length; jj++) output->ffdata->data[psd->length*ii + jj] = (REAL4)(psd->data[jj]*psdfactor*output->ffnormalization);
      
   } /* for ii < numfbins */
   
   //Destroy stuff
   XLALDestroyREAL4Vector(x);
   XLALDestroyREAL4Vector(psd);
   XLALDestroyREAL4Window(win);
   
   
} /* makeSecondFFT() */



//////////////////////////////////////////////////////////////
// Determine the average of the noise power in each frequency bin across the band
REAL4 avgTFdataBand(REAL4Vector *backgrnd, INT4 numfbins, INT4 numffts, INT4 binmin, INT4 binmax)
{
   
   INT4 ii;
   REAL4Vector *aveNoiseInTime = XLALCreateREAL4Vector(numffts);
   REAL4Vector *rngMeansOverBand = XLALCreateREAL4Vector((UINT4)(binmax-binmin));
   if (aveNoiseInTime==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numffts);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   } else if (rngMeansOverBand==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, (UINT4)(binmax-binmin));
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)aveNoiseInTime->length; ii++) {
      memcpy(rngMeansOverBand->data, &(backgrnd->data[ii*numfbins + binmin]), sizeof(*rngMeansOverBand->data)*rngMeansOverBand->length);
      aveNoiseInTime->data[ii] = calcMean(rngMeansOverBand);
      if (XLAL_IS_REAL4_FAIL_NAN(aveNoiseInTime->data[ii])) {
         fprintf(stderr,"%s: calcMean() failed.\n", __func__);
         XLAL_ERROR_REAL4(XLAL_EFUNC);
      }
   } /* for ii < aveNoiseInTime->length */
   
   REAL4 avgTFdata = calcMean(aveNoiseInTime);
   if (XLAL_IS_REAL4_FAIL_NAN(avgTFdata)) {
      fprintf(stderr,"%s: calcMean() failed.\n", __func__);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }
   
   //Destroy stuff
   XLALDestroyREAL4Vector(aveNoiseInTime);
   XLALDestroyREAL4Vector(rngMeansOverBand);
   
   return avgTFdata;
   
} /* avgTFdataBand() */


//////////////////////////////////////////////////////////////
// Determine the rms of the noise power in each frequency bin across the band
REAL4 rmsTFdataBand(REAL4Vector *backgrnd, INT4 numfbins, INT4 numffts, INT4 binmin, INT4 binmax)
{
   
   INT4 ii;
   REAL4Vector *aveNoiseInTime = XLALCreateREAL4Vector(numffts);
   REAL4Vector *rngMeansOverBand = XLALCreateREAL4Vector((UINT4)(binmax-binmin));
   if (aveNoiseInTime==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numffts);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   } else if (rngMeansOverBand==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, (UINT4)(binmax-binmin));
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)aveNoiseInTime->length; ii++) {
      memcpy(rngMeansOverBand->data, &(backgrnd->data[ii*numfbins + binmin]), sizeof(*rngMeansOverBand->data)*rngMeansOverBand->length);
      aveNoiseInTime->data[ii] = calcMean(rngMeansOverBand);
      if (XLAL_IS_REAL4_FAIL_NAN(aveNoiseInTime->data[ii])) {
         fprintf(stderr,"%s: calcMean() failed.\n", __func__);
         XLAL_ERROR_REAL4(XLAL_EFUNC);
      }
   } /* for ii < aveNoiseInTime->length */
   
   REAL4 rmsTFdata = calcRms(aveNoiseInTime);
   if (XLAL_IS_REAL4_FAIL_NAN(rmsTFdata)) {
      fprintf(stderr,"%s: calcRms() failed.\n", __func__);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }
   
   //Destroy stuff
   XLALDestroyREAL4Vector(aveNoiseInTime);
   XLALDestroyREAL4Vector(rngMeansOverBand);
   
   return rmsTFdata;
   
} /* rmsTFdataBand() */


//////////////////////////////////////////////////////////////
// Measure of the average noise power in each 2st FFT frequency bin
void ffPlaneNoise(REAL4Vector *aveNoise, inputParamsStruct *input, INT4Vector *sftexist, REAL4Vector *backgrnd, REAL4Vector *antweights, REAL4FFTPlan *plan, REAL8 *normalization)
{
   
   INT4 ii, jj, numfbins, numffts, numfprbins;
   REAL8 invsumofweights = 0.0;
   REAL8 sumofweights = 0.0;
   
   numfbins = (INT4)(round(input->fspan*input->Tcoh+2.0*input->dfmax*input->Tcoh)+12+1);    //Number of frequency bins
   numffts = (INT4)antweights->length;                      //Number of FFTs
   numfprbins = (INT4)floor(numffts*0.5)+1;                 //number of 2nd fft frequency bins
   
   //Set up for making the PSD
   memset(aveNoise->data, 0, sizeof(REAL4)*aveNoise->length);
   
   REAL4Window *win = XLALCreateHannREAL4Window(numffts);  //Window function
   REAL4Vector *psd = XLALCreateREAL4Vector(numfprbins);   //Current PSD calculation
   if (win==NULL) {
      fprintf(stderr,"%s: XLALCreateHannREAL4Window(%d) failed.\n", __func__, numffts);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (psd==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numfprbins);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   REAL4 winFactor = 8.0/3.0;

   //Average each SFT across the frequency band, also compute normalization factor
   REAL4Vector *aveNoiseInTime = XLALCreateREAL4Vector(numffts);
   REAL4Vector *rngMeansOverBand = XLALCreateREAL4Vector(numfbins);
   if (aveNoiseInTime==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numffts);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (rngMeansOverBand==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numfbins);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   memset(aveNoiseInTime->data, 0, sizeof(REAL4)*numffts);
   for (ii=0; ii<(INT4)aveNoiseInTime->length; ii++) {
      if (sftexist->data[ii]!=0) {
         memcpy(rngMeansOverBand->data, &(backgrnd->data[ii*numfbins]), sizeof(*rngMeansOverBand->data)*rngMeansOverBand->length);
         //aveNoiseInTime->data[ii] = calcMean(rngMeansOverBand);  //comment out?
         aveNoiseInTime->data[ii] = calcMedian(rngMeansOverBand);
         //aveNoiseInTime->data[ii] = (REAL4)(calcRms(rngMeansOverBand));  //For exp dist and large blksize this approaches the mean
         if (XLAL_IS_REAL4_FAIL_NAN(aveNoiseInTime->data[ii])) {
            fprintf(stderr,"%s: calcMean() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         
         if (input->noiseWeightOff==0) sumofweights += (antweights->data[ii]*antweights->data[ii])/(aveNoiseInTime->data[ii]*aveNoiseInTime->data[ii]);
         else sumofweights += (antweights->data[ii]*antweights->data[ii]);
      }
   } /* for ii < aveNoiseInTime->length */
   invsumofweights = 1.0/sumofweights;
   
   //Load time series of powers, normalize, mean subtract and Hann window
   REAL4Vector *x = XLALCreateREAL4Vector(aveNoiseInTime->length);
   REAL4Vector *multiplicativeFactor = XLALCreateREAL4Vector(aveNoiseInTime->length);
   if (x==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, aveNoiseInTime->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (multiplicativeFactor==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, aveNoiseInTime->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   memset(multiplicativeFactor->data, 0, sizeof(REAL4)*multiplicativeFactor->length);
   REAL4 psdfactor = winFactor;
   
   for (ii=0; ii<(INT4)x->length; ii++) {
      if (aveNoiseInTime->data[ii] != 0.0) {
         if (input->noiseWeightOff==0) multiplicativeFactor->data[ii] = win->data->data[ii]*antweights->data[ii]/aveNoiseInTime->data[ii]*invsumofweights;
         else multiplicativeFactor->data[ii] = win->data->data[ii]*antweights->data[ii]*invsumofweights;
      }
   }
   
   //Previous version with no correlation between SFTs when there is overlap
   /* for (ii=0; ii<100000; ii++) {
      for (jj=0; jj<(INT4)x->length; jj++) {
         if (aveNoiseInTime->data[jj] != 0.0) x->data[jj] = multiplicativeFactor->data[jj]*(expRandNum(aveNoiseInTime->data[jj], rng)/aveNoiseInTime->data[jj]-1.0);
         else x->data[jj] = 0.0;
      }
      
      //Do the FFT
      INT4 check = XLALREAL8PowerSpectrum(psd,x,plan);
      if (check != 0) {
         printf("Something wrong with second PSD in background estimate...\n");
         fprintf(LOG,"Something wrong with second PSD in background estimate...\n");
      }
      
      //Rescale and sum into the bins
      for (jj=0; jj<(INT4)aveNoise->length; jj++) {
         aveNoise->data[jj] += psd->data[jj];
      }
   } */
   
   //New version. Computes expected background with correlation estimate from Hann windowed and overlapped (must be 50% or 0%) SFTs
   REAL8 correlationfactor = 0.0;
   for (ii=0; ii<(INT4)floor(win->data->length*(input->SFToverlap/input->Tcoh)-1); ii++) correlationfactor += win->data->data[ii]*win->data->data[ii + (INT4)((1.0-(input->SFToverlap/input->Tcoh))*win->data->length)];
   correlationfactor /= win->sumofsquares;
   REAL8 corrfactorsquared = correlationfactor*correlationfactor;
   REAL8 prevnoiseval = 0.0;
   REAL8 noiseval = 0.0;
   for (ii=0; ii<4000; ii++) {
      memset(x->data, 0, sizeof(REAL4)*x->length);
      for (jj=0; jj<(INT4)x->length; jj++) {
         if (sftexist->data[jj] != 0) {
            //To create the correlations
            noiseval = expRandNum(aveNoiseInTime->data[jj], input->rng);
            if (XLAL_IS_REAL8_FAIL_NAN(noiseval)) {
               fprintf(stderr, "%s: expRandNum() failed.\n", __func__);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
            if (jj>0 && sftexist->data[jj-1]!=0) {
               noiseval *= (1.0-corrfactorsquared);
               noiseval += corrfactorsquared*prevnoiseval;
            }
            x->data[jj] = (REAL4)(noiseval/aveNoiseInTime->data[jj]-1.0);
            prevnoiseval = noiseval;
         }
      } /* for jj < x->length */
      
      //Window and rescale because of antenna and noise weights
      if (input->useSSE) {
         sseSSVectorMultiply(x, x, multiplicativeFactor);
         if (xlalErrno!=0) {
            fprintf(stderr,"%s: sseSSVectorMultiply() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
      }
      else for (jj=0; jj<(INT4)x->length; jj++) x->data[jj] *= multiplicativeFactor->data[jj];
      
      //Do the FFT
      if ( (XLALREAL4PowerSpectrum(psd, x, plan)) != 0) {
         fprintf(stderr,"%s: XLALREAL4PowerSpectrum() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      
      //Sum into the bins
      if (input->useSSE) {
         sseSSVectorSum(aveNoise, aveNoise, psd);
         if (xlalErrno!=0) {
            fprintf(stderr,"%s: sseSSVectorSum() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
      }
      else for (jj=0; jj<(INT4)aveNoise->length; jj++) aveNoise->data[jj] += psd->data[jj];
   } /* for ii < 4000 */
   
   //Average and rescale
   REAL4 averageRescaleFactor = 2.5e-4*psdfactor*(1.0+2.0*corrfactorsquared);
   if (input->useSSE) {
      sseScaleREAL4Vector(aveNoise, aveNoise, averageRescaleFactor);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: sseScaleREAL4Vector() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
   }
   else for (ii=0; ii<(INT4)aveNoise->length; ii++) aveNoise->data[ii] *= averageRescaleFactor;
   
   //Fix 0th and end bins (0 only for odd x->length, 0 and end for even x->length)
   if (GSL_IS_EVEN(x->length)==1) {
      aveNoise->data[0] *= 2.0;
      aveNoise->data[aveNoise->length-1] *= 2.0;
   } else {
      aveNoise->data[0] *= 2.0;
   }
   
   //comment this out
   //FILE *BACKGRND = fopen("./output/background.dat","w");
   
   //Compute normalization
   *(normalization) = 1.0/(calcMean(aveNoise));
   for (ii=0; ii<(INT4)aveNoise->length; ii++) {
      aveNoise->data[ii] *= *(normalization);
      //fprintf(BACKGRND,"%.6f\n",aveNoise->data[ii]);
   }
   
   //Extra factor for normalization to 1.0 (empirically determined)
   // *(normalization) /= 1.0245545525190294;
   //*normalization /= 1.040916688722758;
   *normalization /= 1.08;
   
   //fclose(BACKGRND);
   
   //Destroy stuff
   XLALDestroyREAL4Vector(x);
   XLALDestroyREAL4Vector(psd);
   XLALDestroyREAL4Window(win);
   XLALDestroyREAL4Vector(aveNoiseInTime);
   XLALDestroyREAL4Vector(rngMeansOverBand);
   XLALDestroyREAL4Vector(multiplicativeFactor);

} /* ffPlaneNoise() */



//For testing purposes only!!!! Don't use
REAL4Vector * simpleTFdata(REAL8 fsig, REAL8 period, REAL8 moddepth, REAL8 Tcoh, REAL8 Tobs, REAL8 SFToverlap, REAL8 fminimum, REAL8 fmaximum, REAL8 sqrtSh)
{
   
   INT4 numfbins = (INT4)(round((fmaximum-fminimum)*Tcoh)+1);   //Number of frequency bins
   INT4 numffts = (INT4)floor(Tobs/(Tcoh-SFToverlap)-1); //Number of FFTs
   
   REAL4Vector *output = XLALCreateREAL4Vector(numfbins*numffts);
   
   //Initialize the random number generator
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
   if (rng==NULL) {
      fprintf(stderr,"%s: gsl_rng_alloc() failed.\n", __func__);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   gsl_rng_set(rng, 0);
   
   INT4 ii, jj;
   REAL8 correlationfactor = 0.167, corrfactorsquared = correlationfactor*correlationfactor;
   for (ii=0; ii<numffts; ii++) {
      for (jj=0; jj<numfbins; jj++) {
         if (ii==0) {
            output->data[jj] = expRandNum(sqrtSh, rng);
         } else {
            output->data[ii*numfbins + jj] = corrfactorsquared*output->data[(ii-1)*numfbins + jj] + (1.0-corrfactorsquared)*expRandNum(sqrtSh, rng);
         }

      }
   }
   
   for (ii=0; ii<numffts; ii++) {
      REAL8 fbin = fsig + moddepth*sin(LAL_TWOPI*((ii+1)*SFToverlap)/period) - fminimum;
      for (jj=0; jj<numfbins; jj++) output->data[ii*numfbins + jj] += 0.03*(2.0/3.0)*Tcoh*sqsincxoverxsqminusone(fbin*Tcoh-(REAL8)jj);
   }
   
   //Destroy stuff
   gsl_rng_free(rng);
   
   return output;
   
}



//Convert the gengetopt_args_info struct into something less complicated and select appropriate IFO(s)
INT4 readTwoSpectInputParams(inputParamsStruct *params, struct gengetopt_args_info args_info)
{
   
   INT4 ii;
   
   //Defaults given or option passed
   params->blksize = args_info.blksize_arg;                                      //Block size of SFT running median (bins)
   params->dopplerMultiplier = args_info.dopplerMultiplier_arg;                  //Velocity of Earth multiplier (default = 1.0)
   params->mintemplatelength = args_info.minTemplateLength_arg;                  //Minimum number of template weights (pixels)
   params->maxtemplatelength = args_info.maxTemplateLength_arg;                  //Maximum number of template weights (pixels)
   params->ihsfactor = args_info.ihsfactor_arg;                                  //IHS folding factor (default = 5)
   params->rootFindingMethod = args_info.BrentsMethod_given;                     //Use Brent's method (default = 0)
   params->antennaOff = args_info.antennaOff_given;                              //Antenna pattern off (default = 0)
   params->noiseWeightOff = args_info.noiseWeightOff_given;                      //Noise weighting off (default = 0)
   params->calcRthreshold = args_info.calcRthreshold_given;                      //Directly calculate the R threshold value (defualt = 0)
   params->markBadSFTs = args_info.markBadSFTs_given;                            //Mark bad SFTs (default = 0)
   params->FFTplanFlag = args_info.FFTplanFlag_arg;                              //FFTW plan flag
   params->printAllULvalues = args_info.allULvalsPerSkyLoc_given;                //Output all UL values at each sky location (default = 0)
   params->fastchisqinv = args_info.fastchisqinv_given;                          //Use faster chi-sq inversion (default = 0)
   params->useSSE = args_info.useSSE_given;                                      //Use SSE optimized functions (dafualt = 0)
   params->followUpOutsideULrange = args_info.followUpOutsideULrange_given;      //Follow up outliers outside of UL range (default = 0)
   params->validateSSE = args_info.validateSSE_given;                            //Validate SSE functions (default = 0)
   params->noNotchHarmonics = args_info.noNotchHarmonics_given;                  //Do not notch the daily/sidereal harmonics (default = 0)
   params->harmonicNumToSearch = args_info.harmonicNumToSearch_arg;              //Search the number of harmonics specified by the Pmin-->Pmax range (default = 1 meaning search only the range of Pmin-->Pmax)
   params->ULsolver = args_info.ULsolver_arg;                                    //Solver function for UL calculation (default = 0)
   
   //Non-default arguments
   if (args_info.Tcoh_given) params->Tcoh = args_info.Tcoh_arg;                  //SFT coherence time (s)
   else {
      fprintf(stderr, "%s: a SFT coherence time must be specified.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   if (args_info.SFToverlap_given) params->SFToverlap = args_info.SFToverlap_arg; //SFT overlap (s)
   else {
      fprintf(stderr, "%s: the SFT overlap time must be specified.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   if (args_info.Tobs_given) params->Tobs = args_info.Tobs_arg;                  //Total observation time (s)
   else {
      fprintf(stderr, "%s: an observation time must be specified.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   if (args_info.fmin_given) params->fmin = args_info.fmin_arg;                  //Minimum frequency to search (Hz)
   else {
      fprintf(stderr, "%s: a minimum frequency to search must be specified.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   if (args_info.fspan_given) params->fspan = args_info.fspan_arg;               //Maximum frequency to search (Hz)
   else {
      fprintf(stderr, "%s: a frequency span must be specified.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   if (args_info.t0_given) params->searchstarttime = args_info.t0_arg;           //GPS start time of the search (s)
   else {
      fprintf(stderr, "%s: a search start time must be specified.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   if (args_info.Pmin_given) params->Pmin = args_info.Pmin_arg;                  //Minimum period to search (s)
   else {
      fprintf(stderr, "%s: a minimum period to search must be specified.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   if (args_info.Pmax_given) params->Pmax = args_info.Pmax_arg;                  //Maximum period to search (s)
   else {
      fprintf(stderr, "%s: a maximum period to search must be specified.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   if (args_info.dfmin_given) params->dfmin = args_info.dfmin_arg;               //Minimum modulation depth to search (Hz)
   else {
      fprintf(stderr, "%s: a minimum modulation depth to search must be specified.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   if (args_info.dfmax_given) params->dfmax = args_info.dfmax_arg;               //Maximum modulation depth to search (Hz)
   else {
      fprintf(stderr, "%s: a maximum modulation depth to search must be specified.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   if (args_info.ihsfar_given) params->ihsfar = args_info.ihsfar_arg;            //Incoherent harmonic sum false alarm rate
   else {
      fprintf(stderr, "%s: the IHS FAR must be specified.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   if (args_info.tmplfar_given) params->templatefar = args_info.tmplfar_arg;     //Template false alarm rate
   else {
      fprintf(stderr, "%s: the template FAR must be specified.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   
   if (args_info.keepOnlyTopNumIHS_given) params->keepOnlyTopNumIHS = args_info.keepOnlyTopNumIHS_arg;         //Keep only top X IHS candidates
   else params->keepOnlyTopNumIHS = -1;
   if (args_info.simpleBandRejection_given) params->simpleSigmaExclusion = args_info.simpleBandRejection_arg;  //Simple band rejection (default off)
   if (args_info.lineDetection_given) params->lineDetection = args_info.lineDetection_arg;                     //Line detection
   
   
   params->log10templatefar = log10(params->templatefar);                        //log_10(template FAR)
   
   //Settings for IHS FOM
   //Exit with error if neither is chosen
   if (args_info.ihsfomfar_given) params->ihsfomfar = args_info.ihsfomfar_arg;   //IHS figure of merit false alarm rate
   else params->ihsfomfar = 0.0;
   if (args_info.ihsfom_given) params->ihsfom = args_info.ihsfom_arg;            //IHS figure of merit threshold value
   else params->ihsfom = 0.0;
   if ((params->ihsfom!=0.0 && params->ihsfomfar!=0.0) || (params->ihsfom==0.0 && params->ihsfomfar==0.0)) {
      fprintf(stderr, "%s: You must choose either the IHS FOM FAR argument or the IHS FOM argument.\n", __func__);
      XLAL_ERROR(XLAL_EINVAL);
   }
   
   //Blocksize should be an odd number
   if (params->blksize % 2 != 1) params->blksize += 1;
   
   // Warnings when using hidden flags
   if (args_info.antennaOff_given) {
      fprintf(LOG,"WARNING: Antenna pattern weights are all being set to 1.0\n");
      fprintf(stderr,"WARNING: Antenna pattern weights are all being set to 1.0\n");
   }
   if (args_info.noiseWeightOff_given) {
      fprintf(LOG,"WARNING: Noise weights are all being set to 1.0\n");
      fprintf(stderr,"WARNING: Noise weights are all being set to 1.0\n");
   }
   if (args_info.gaussTemplatesOnly_given) {
      fprintf(LOG,"WARNING: Only Gaussian templates will be used\n");
      fprintf(stderr,"WARNING: Only Gaussian templates will be used\n");
   }
   if (args_info.IHSonly_given) {
      fprintf(LOG,"WARNING: Only IHS stage is being used\n");
      fprintf(stderr,"WARNING: Only IHS stage is being used\n");
   }
   if (args_info.markBadSFTs_given) {
      fprintf(LOG,"WARNING: Marking bad SFTs\n");
      fprintf(stderr,"WARNING: Marking bad SFTs\n");
   }
   
   //Adjust parameter space search values, if necessary
   if (params->Pmax < params->Pmin) {
      REAL4 tempP = params->Pmax;
      params->Pmax = params->Pmin;
      params->Pmin = tempP;
      fprintf(LOG,"WARNING! Maximum period is smaller than minimum period... switching the two\n");
      fprintf(stderr,"WARNING! Maximum period is smaller than minimum period... switching the two\n");
   }
   if (params->dfmax < params->dfmin) {
      REAL4 tempdf = params->dfmax;
      params->dfmax = params->dfmin;
      params->dfmin = tempdf;
      fprintf(LOG,"WARNING! Maximum modulation depth is smaller than minimum modulation depth... switching the two\n");
      fprintf(stderr,"WARNING! Maximum modulation depth is smaller than minimum modulation depth... switching the two\n");
   }
   if (params->Pmax > 0.2*(params->Tobs)) {
      params->Pmax = 0.2*(params->Tobs);
      fprintf(LOG,"WARNING! Adjusting input maximum period to 1/5 the observation time!\n");
      fprintf(stderr,"WARNING! Adjusting input maximum period to 1/5 the observation time!\n");
   }
   if (params->Pmin < 2.0*3600) {
      params->Pmin = 2.0*3600;
      fprintf(LOG,"WARNING! Adjusting input minimum period to 2 hours!\n");
      fprintf(stderr,"WARNING! Adjusting input minimum period to 2 hours!\n");
   }
   if (params->dfmax > maxModDepth(params->Pmax, params->Tcoh)) {
      params->dfmax = floor(maxModDepth(params->Pmax, params->Tcoh)*(params->Tcoh))/(params->Tcoh);
      fprintf(LOG,"WARNING! Adjusting input maximum modulation depth due to maximum modulation depth allowed\n");
      fprintf(stderr,"WARNING! Adjusting input maximum modulation depth due to maximum modulation depth allowed\n");
   }
   if (params->dfmin < 0.5/(params->Tcoh)) {
      params->dfmin = 0.5/(params->Tcoh);
      fprintf(LOG,"WARNING! Adjusting input minimum modulation depth to 1/2 a frequency bin!\n");
      fprintf(stderr,"WARNING! Adjusting input minimum modulation depth to 1/2 a frequency bin!\n");
   }
   
   //Adjustments for improper modulation depth inputs
   params->dfmin = 0.5*round(2.0*params->dfmin*params->Tcoh)/params->Tcoh;
   params->dfmax = 0.5*round(2.0*params->dfmax*params->Tcoh)/params->Tcoh;
   
   //Upper limit settings
   if (args_info.ULfmin_given) params->ULfmin = args_info.ULfmin_arg;            //Upper limit minimum frequency (Hz)
   else params->ULfmin = params->fmin;
   if (args_info.ULfspan_given) params->ULfspan = args_info.ULfspan_arg;         //Upper limit maximum frequency (Hz)
   else params->ULfspan = params->fspan;
   if (args_info.ULminimumDeltaf_given) params->ULmindf = args_info.ULminimumDeltaf_arg;     //Upper limit minimum modulation depth (Hz)
   else params->ULmindf = params->dfmin;
   if (args_info.ULmaximumDeltaf_given) params->ULmaxdf = args_info.ULmaximumDeltaf_arg;     //Upper limit maximum modulation depth (Hz)
   else params->ULmaxdf = params->dfmax;
   
   //Print to stderr the parameters of the search
   fprintf(stderr,"Tobs = %f sec\n",params->Tobs);
   fprintf(stderr,"Tcoh = %f sec\n",params->Tcoh);
   fprintf(stderr,"SFToverlap = %f sec\n",params->SFToverlap);
   fprintf(stderr,"fmin = %f Hz\n",params->fmin);
   fprintf(stderr,"fspan = %f Hz\n",params->fspan);
   fprintf(stderr,"Pmin = %f s\n",params->Pmin);
   fprintf(stderr,"Pmax = %f s\n",params->Pmax);
   fprintf(stderr,"dfmin = %f Hz\n",params->dfmin);
   fprintf(stderr,"dfmax = %f Hz\n",params->dfmax);
   fprintf(stderr,"Running median blocksize = %d\n",params->blksize);
   fprintf(stderr,"FFT plan flag = %d\n", params->FFTplanFlag);
   if (args_info.ihsfomfar_given) fprintf(stderr,"IHS FOM FAR = %f\n", params->ihsfomfar);
   else fprintf(stderr,"IHS FOM = %f\n", params->ihsfom);

   //Root finding method
   if (args_info.BrentsMethod_given == 0) {
      fprintf(LOG,"Using Newton's method for root finding.\n");
      fprintf(stderr,"Using Newton's method for root finding.\n");
   } else {
      fprintf(LOG,"Using Brent's method for root finding.\n");
      fprintf(stderr,"Using Brent's method for root finding.\n");
   }
   
   //SFT type standard or vladimir (Vladimir's SFT generation program has a different normalization factor than standard v2)
   params->sftType = XLALCalloc(strlen(args_info.sftType_arg)+1, sizeof(*(params->sftType)));
   if (params->sftType==NULL) {
      fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", __func__, sizeof(*(params->sftType)));
      XLAL_ERROR(XLAL_ENOMEM);
   }
   sprintf(params->sftType, "%s", args_info.sftType_arg);
   if (strcmp(params->sftType, "standard")==0) {
      fprintf(LOG,"sftType = %s\n", params->sftType);
      fprintf(stderr,"sftType = %s\n", params->sftType);
   } else if (strcmp(params->sftType, "vladimir")==0) {
      fprintf(LOG,"sftType = %s\n", params->sftType);
      fprintf(stderr,"sftType = %s\n", params->sftType);
   } else {
      fprintf(stderr, "%s: Not using valid type of SFT! Expected 'standard' or 'vladimir' not %s.\n", __func__, params->sftType);
      XLAL_ERROR(XLAL_EINVAL);
   }
   
   //Interferometer from the IFO given parameter. Could be more than 1, but this functionality is not fully implemented, so we exit with
   //an error if this is done
   params->numofIFOs = args_info.IFO_given;
   if (params->numofIFOs>1) {
      fprintf(stderr, "%s: Only one IFO is allowed at the present time.\n", __func__);
      XLAL_ERROR(XLAL_EINVAL);
   }
   if (params->numofIFOs==0) {
      fprintf(stderr, "%s: You must specify an IFO.\n", __func__);
      XLAL_ERROR(XLAL_EINVAL);
   }
   CHAR *IFO = NULL;
   for (ii=0; ii<params->numofIFOs; ii++) {
      IFO = XLALCalloc(strlen(args_info.IFO_arg[ii])+1, sizeof(*IFO));
      if (IFO==NULL) {
         fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", __func__, sizeof(*IFO));
         XLAL_ERROR(XLAL_ENOMEM);
      }
      sprintf(IFO, "%s", args_info.IFO_arg[ii]);
      if (strcmp("L1", IFO)==0) {
         fprintf(LOG,"IFO = %s\n", IFO);
         fprintf(stderr,"IFO = %s\n", IFO);
         params->det[ii] = lalCachedDetectors[LAL_LLO_4K_DETECTOR]; //L1
      } else if (strcmp("H1", IFO)==0) {
         fprintf(LOG,"IFO = %s\n", IFO);
         fprintf(stderr,"IFO = %s\n", IFO);
         params->det[ii] = lalCachedDetectors[LAL_LHO_4K_DETECTOR]; //H1
      } else if (strcmp("V1", IFO)==0) {
         fprintf(LOG,"IFO = %s\n", IFO);
         fprintf(stderr,"IFO = %s\n", IFO);
         params->det[ii] = lalCachedDetectors[LAL_VIRGO_DETECTOR]; //V1
      } else {
         fprintf(stderr, "%s: Not using valid interferometer! Expected 'H1', 'L1', or 'V1' not %s.\n", __func__, IFO);
         XLAL_ERROR(XLAL_EINVAL);
      }
      XLALFree((CHAR*)IFO);
   }
   
   
   //Read in file names for ephemeris files
   if (!args_info.ephemDir_given) {
      fprintf(stderr, "%s: An ephemeris directory path must be specified.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   if (!args_info.ephemYear_given) {
      fprintf(stderr, "%s: An ephemeris year/type suffix must be specified.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   earth_ephemeris = XLALCalloc(strlen(args_info.ephemDir_arg)+25, sizeof(*earth_ephemeris));
   sun_ephemeris = XLALCalloc(strlen(args_info.ephemDir_arg)+25, sizeof(*sun_ephemeris));
   if (earth_ephemeris==NULL) {
      fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", __func__, sizeof(*earth_ephemeris));
      XLAL_ERROR(XLAL_ENOMEM);
   } else if (sun_ephemeris==NULL) {
      fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", __func__, sizeof(*sun_ephemeris));
      XLAL_ERROR(XLAL_ENOMEM);
   }
   sprintf(earth_ephemeris, "%s/earth%s.dat", args_info.ephemDir_arg, args_info.ephemYear_arg);
   sprintf(sun_ephemeris, "%s/sun%s.dat", args_info.ephemDir_arg, args_info.ephemYear_arg);

   //SFT input
   if (args_info.sftDir_given && !args_info.sftFile_given) {
      sft_dir = XLALCalloc(strlen(args_info.sftDir_arg)+20, sizeof(*sft_dir));
      if (sft_dir==NULL) {
	 fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", __func__, sizeof(*sft_dir));
	 XLAL_ERROR(XLAL_ENOMEM);
      }
      sprintf(sft_dir, "%s/*.sft", args_info.sftDir_arg);
   } else if (!args_info.sftDir_given && args_info.sftFile_given) {
      sft_dir = XLALCalloc(strlen(args_info.sftFile_arg)+2, sizeof(*sft_dir));
      if (sft_dir==NULL) {
	 fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", __func__, sizeof(*sft_dir));
	 XLAL_ERROR(XLAL_ENOMEM);
      }
      sprintf(sft_dir, "%s", args_info.sftFile_arg);
   } else if ((args_info.sftDir_given && args_info.sftFile_given) || !(args_info.sftDir_given && args_info.sftFile_given)) {
      fprintf(stderr, "%s: One of either sftDir or sftFile must be given but not both or neither.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   
   return 0;
   
} /* readTwoSepctInputParams() */



