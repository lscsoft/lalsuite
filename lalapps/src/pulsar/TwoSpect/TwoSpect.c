/*
*  Copyright (C) 2010, 2011 Evan Goetz
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


//Global variables
FILE *LOG = NULL, *ULFILE = NULL;
CHAR *earth_ephemeris = NULL, *sun_ephemeris = NULL, *sft_dir = NULL;


//Main program
int main(int argc, char *argv[])
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj;//, kk, ll;       //counter variables
   LALStatus status;          //LALStatus structure
   status.statusPtr = NULL;   //Set statuspointer to NULL
   char s[20000], t[20000];     //Path and file name to LOG and ULFILE
   time_t programstarttime, programendtime;
   struct tm *ptm;
   
   time(&programstarttime);
   ptm = localtime(&programstarttime);
   
   //Turn off gsl error handler
   gsl_set_error_handler_off();
   
   //Initiate command line interpreter and config file loader
   struct gengetopt_args_info args_info;
   struct cmdline_parser_params *configparams;
   configparams = cmdline_parser_params_create();
   configparams->initialize = 0;
   configparams->override = 1;
   if ( cmdline_parser(argc, argv, &args_info) ) {
      fprintf(stderr, "%s: cmdline_parser() failed.\n", fn);
      XLAL_ERROR(fn, XLAL_FAILURE);
   }
   if ( args_info.config_given ) {
      if ( cmdline_parser_config_file(args_info.config_arg, &args_info, configparams) ) {
         fprintf(stderr, "%s: cmdline_parser_config_file() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_FAILURE);
      }
   }
   
   //Set lalDebugLevel to user input or 0 if no input
   lalDebugLevel = args_info.verbosity_arg;
   
   //Create directory
   mkdir(args_info.outdirectory_arg, 0777);
   snprintf(s, 20000, "%s/%s", args_info.outdirectory_arg, args_info.outfilename_arg);
   snprintf(t, 20000, "%s/%s", args_info.outdirectory_arg, args_info.ULfilename_arg);
   
   //Open log file
   LOG = fopen(s,"w");
   if (LOG==NULL) {
      fprintf(stderr, "%s: Log file could not be opened.\n", fn);
      XLAL_ERROR(fn, XLAL_EINVAL);
   }
   
   //print start time
   fprintf(stderr, "Program executed on %s\n", asctime(ptm));
   fprintf(LOG, "Program executed on %s\n", asctime(ptm));
   
   //Print out the inputs and outputs
   fprintf(stderr, "Input parameters file: %s\n", args_info.config_arg);
   fprintf(LOG, "Input parameters file: %s\n", args_info.config_arg);
   fprintf(stderr, "Input SFTs: %s/%s\n", args_info.sftDir_arg, "*.sft");
   fprintf(LOG, "Input SFTs: %s/%s\n", args_info.sftDir_arg, "*.sft");
   fprintf(stderr, "Output directory: %s\n", args_info.outdirectory_arg);
   fprintf(LOG, "Output directory: %s\n", args_info.outdirectory_arg);
   
   //Allocate input parameters structure memory
   inputParamsStruct *inputParams = new_inputParams();
   if (inputParams==NULL) {
      fprintf(stderr, "%s: new_inputParams() failed.\n", fn);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   
   //Read TwoSpect input parameters
   if ( (readTwoSpectInputParams(inputParams, args_info)) != 0 ) {
      fprintf(stderr, "%s: readTwoSepctInputParams() failed.\n", fn);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   
   //Initialize ephemeris data structure
   EphemerisData *edat = XLALInitBarycenter(earth_ephemeris, sun_ephemeris);
   if (edat==NULL) {
      fprintf(stderr, "%s: XLALInitBarycenter() failed.\n", fn);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   
   //Maximum orbital earth speed in units of c from start of S6 TwoSpect data for 104 weeks total time
   REAL4 detectorVmax = CompDetectorVmax(931081500.0+inputParams->SFToverlap, inputParams->Tcoh, inputParams->SFToverlap, 62899200.0-inputParams->SFToverlap, inputParams->det, edat);
   if (xlalErrno!=0) {
      fprintf(stderr, "%s: CompDetectorVmax() failed.\n", fn);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   
   //Parameters for the sky-grid from a point/polygon or a sky-grid file
   if ((args_info.skyRegion_given && args_info.skyRegionFile_given) || (!args_info.skyRegion_given && !args_info.skyRegionFile_given)) {
      fprintf(stderr, "%s: You must choose either the the sky region (point or polygon) *or* a file.\n", fn);
      XLAL_ERROR(fn, XLAL_EINVAL);
   }
   CHAR *sky = NULL;
   if (args_info.skyRegion_given) {
      sky = XLALCalloc(strlen(args_info.skyRegion_arg)+1, sizeof(*sky));
      if (sky==NULL) {
         fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", fn, sizeof(*sky));
         XLAL_ERROR(fn, XLAL_ENOMEM);
      }
      sprintf(sky, "%s", args_info.skyRegion_arg);
      fprintf(LOG, "Sky region = %s\n", sky);
      fprintf(stderr, "Sky region = %s\n", sky);
   } else {
      sky = XLALCalloc(strlen(args_info.skyRegionFile_arg)+1, sizeof(*sky));
      if (sky==NULL) {
         fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", fn, sizeof(*sky));
         XLAL_ERROR(fn, XLAL_ENOMEM);
      }
      sprintf(sky, "%s", args_info.skyRegionFile_arg);
      fprintf(LOG, "Sky file = %s\n", sky);
      fprintf(stderr, "Sky file = %s\n", sky);
   }
   DopplerSkyScanInit scanInit = empty_DopplerSkyScanInit;
   DopplerSkyScanState scan = empty_DopplerSkyScanState;
   PulsarDopplerParams dopplerpos;
   if (args_info.skyRegion_given) {
      scanInit.gridType = 1;     //Default value for an approximate-isotropic grid
      scanInit.skyRegionString = sky;      //"allsky" = Default value for all-sky search
      scanInit.numSkyPartitions = 1;   //Default value so sky is not broken into chunks
      scanInit.Freq = args_info.fmin_arg+0.5*args_info.fspan_arg;  //Mid-point of the frequency band
      scanInit.dAlpha = 0.5/((inputParams->fmin+0.5*inputParams->fspan) * inputParams->Tcoh * detectorVmax);
      scanInit.dDelta = scanInit.dAlpha;
   } else {
      scanInit.gridType = 3;
      scanInit.skyGridFile = sky;
      scanInit.numSkyPartitions = 1;   //Default value so sky is not broken into chunks
      scanInit.Freq = args_info.fmin_arg+0.5*args_info.fspan_arg;  //Mid-point of the frequency band
   }
   
   //Initialize the sky-grid
   InitDopplerSkyScan(&status, &scan, &scanInit);
   if (status.statusCode!=0) {
      fprintf(stderr, "%s: InitDopplerSkyScan() failed.\n", fn);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   
   //Start at first location
   if ((XLALNextDopplerSkyPos(&dopplerpos, &scan))!=0) {
      fprintf(stderr, "%s: XLALNextDopplerSkyPos() failed.\n", fn);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   
   //Basic units
   REAL4 tempfspan = inputParams->fspan + (inputParams->blksize-1)/inputParams->Tcoh;
   INT4 tempnumfbins = (INT4)round(tempfspan*inputParams->Tcoh)+1;
   REAL8 templatefarthresh = args_info.tmplfar_arg;
   fprintf(LOG, "FAR for templates = %g\n", templatefarthresh);
   fprintf(stderr, "FAR for templates = %g\n", templatefarthresh);
   
   //Allocate memory for ffdata structure
   ffdataStruct *ffdata = new_ffdata(inputParams);
   if (ffdata==NULL) {
      fprintf(stderr, "%s: new_ffdata() failed.\n", fn);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   
   //Allocate lists of candidates with initially 100 available slots (will check and rescale later, if necessary)
   candidateVector *gaussCandidates1 = new_candidateVector(100);
   candidateVector *gaussCandidates2 = new_candidateVector(100);
   candidateVector *gaussCandidates3 = new_candidateVector(100);
   candidateVector *gaussCandidates4 = new_candidateVector(100);
   candidateVector *exactCandidates1 = new_candidateVector(100);   
   candidateVector *exactCandidates2 = new_candidateVector(100);
   candidateVector *ihsCandidates = new_candidateVector(100);
   UpperLimitVector *upperlimits = new_UpperLimitVector(1);
   if (gaussCandidates1==NULL) {
      fprintf(stderr, "%s: new_CandidateVector(%d) failed.\n", fn, 100);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (gaussCandidates2==NULL) {
      fprintf(stderr, "%s: new_CandidateVector(%d) failed.\n", fn, 100);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (gaussCandidates3==NULL) {
      fprintf(stderr, "%s: new_CandidateVector(%d) failed.\n", fn, 100);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (gaussCandidates4==NULL) {
      fprintf(stderr, "%s: new_CandidateVector(%d) failed.\n", fn, 100);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (exactCandidates1==NULL) {
      fprintf(stderr, "%s: new_CandidateVector(%d) failed.\n", fn, 100);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (exactCandidates2==NULL) {
      fprintf(stderr, "%s: new_CandidateVector(%d) failed.\n", fn, 100);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (ihsCandidates==NULL) {
      fprintf(stderr, "%s: new_CandidateVector(%d) failed.\n", fn, 100);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (upperlimits==NULL) {
      fprintf(stderr, "%s: new_UpperLimitVector(%d) failed.\n", fn, 1);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   
   //Second fft plan, only need to make this once for all the exact templates
   REAL4FFTPlan *secondFFTplan = XLALCreateForwardREAL4FFTPlan(ffdata->numffts, inputParams->FFTplanFlag);
   if (secondFFTplan==NULL) {
      fprintf(stderr, "%s: XLALCreateForwardREAL4FFTPlan(%d,%d) failed.\n", fn, ffdata->numffts, inputParams->FFTplanFlag);
      XLAL_ERROR(fn, XLAL_EFUNC);
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
   /* XLALDestroyREAL4Vector(tfdata);
   tfdata = NULL;
   tfdata = simpleTFdata(100.0, 513864.0, 20.0*3.667e-3, inputParams->Tcoh, inputParams->Tobs, inputParams->SFToverlap, inputParams->fmin-(inputParams->maxbinshift+(inputParams->blksize-1)/2)/inputParams->Tcoh, inputParams->fmin+inputParams->fspan+(inputParams->maxbinshift+(inputParams->blksize-1)/2)/inputParams->Tcoh, 1.0); */
   if (tfdata==NULL) {
      fprintf(stderr, "\n%s: readInSFTs() failed.\n", fn);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   fprintf(LOG, "done\n");
   fprintf(stderr, "done\n");
   /* FILE *rawtfdata = fopen("./output/rawtfdata.dat","w");
   for (ii=0; ii<(INT4)tfdata->length; ii++) fprintf(rawtfdata, "%f\n", tfdata->data[ii]);
   fclose(rawtfdata); */
   
   //Removing bad SFTs
   if (inputParams->markBadSFTs!=0) {
      fprintf(stderr, "Marking and removing bad SFTs... ");
      INT4Vector *removeTheseSFTs = markBadSFTs(tfdata, inputParams);
      if (removeTheseSFTs==NULL) {
         fprintf(stderr, "%s: markBadSFTs() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      removeBadSFTs(tfdata, removeTheseSFTs);
      fprintf(stderr, "done.\n");
      XLALDestroyINT4Vector(removeTheseSFTs);
   }
   /* FILE *rawtfdata = fopen("./output/rawtfdata.dat","w");
   for (ii=0; ii<(INT4)tfdata->length; ii++) fprintf(rawtfdata, "%f\n", tfdata->data[ii]);
   fclose(rawtfdata); */
   
   //Existing SFTs listed in this vector
   INT4Vector *sftexist = XLALCreateINT4Vector(ffdata->numffts);
   if (sftexist==NULL) {
      fprintf(stderr, "\n%s: XLALCreateINT4Vector(%d) failed.\n", fn, ffdata->numffts);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   REAL4 totalincludedsftnumber = 0.0;
   //FILE *timestamps = fopen("./output/timestamps.dat","w");
   for (ii=0; ii<ffdata->numffts; ii++) {
      if (tfdata->data[ii*(ffdata->numfbins+2*inputParams->maxbinshift+inputParams->blksize-1)] == 0.0) sftexist->data[ii] = 0;
      else {
         sftexist->data[ii] = 1;
         totalincludedsftnumber += 1.0;
         //fprintf(timestamps, "%d %d\n", (INT4)round(inputParams->searchstarttime+ii*(inputParams->Tcoh-inputParams->SFToverlap)), 0);
      }
   }
   //fclose(timestamps);
   REAL4 frac_tobs_complete = totalincludedsftnumber/sftexist->length;
   
   //Calculate the running mean values of the SFTs (output here is smaller than initialTFdata). Here,
   //numfbins needs to be the bins you expect to come out of the running means -- the band you are going
   //to search!
   fprintf(LOG, "Assessing background... ");
   fprintf(stderr, "Assessing background... ");
   REAL4Vector *background = XLALCreateREAL4Vector(ffdata->numffts*(ffdata->numfbins + 2*inputParams->maxbinshift));
   if (background==NULL) {
      fprintf(stderr, "\n%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ffdata->numffts*(ffdata->numfbins + 2*inputParams->maxbinshift));
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   tfRngMeans(background, tfdata, ffdata->numffts, ffdata->numfbins + 2*inputParams->maxbinshift, inputParams->blksize);
   if (xlalErrno!=0) {
      fprintf(stderr, "\n%s: tfRngMeans() failed.\n", fn);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   
   //I wrote this to compensate for a bad input of the expected noise floor
   REAL8 backgroundmeannormfactor = 0.0;
   INT8 avefact = 0;
   for (ii=0; ii<(INT4)background->length; ii++) {
      if (background->data[ii]!=0.0) {
         backgroundmeannormfactor += background->data[ii];
         avefact++;
      }
   }
   backgroundmeannormfactor = (REAL8)avefact/backgroundmeannormfactor;
   ffdata->tfnormalization *= backgroundmeannormfactor;
   fprintf(LOG, "done\n");
   fprintf(stderr, "done\n");
   
   //Need to reduce the original TF data to remove the excess bins used for running median calculation. Normalize the TF as the same as the background was normalized
   REAL4Vector *usableTFdata = XLALCreateREAL4Vector(background->length);
   if (usableTFdata==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, background->length);
      XLAL_ERROR(fn, XLAL_EFUNC);
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
   ihsfarStruct *ihsfarstruct = new_ihsfarStruct(maxrows);
   REAL4Vector *detectorVelocities = XLALCreateREAL4Vector(ffdata->numffts);
   INT4Vector *binshifts = XLALCreateINT4Vector(ffdata->numffts);
   REAL4Vector *aveNoise = XLALCreateREAL4Vector(ffdata->numfprbins);
   if (ihsmaxima==NULL) {
      fprintf(stderr, "%s: new_ihsMaxima(%d,%d) failed.\n", fn, ffdata->numfbins, maxrows);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (ihsfarstruct==NULL) {
      fprintf(stderr, "%s: new_ihsfarStruct(%d) failed.\n", fn, maxrows);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (detectorVelocities==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ffdata->numffts);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (binshifts==NULL) {
      fprintf(stderr, "%s: XLALCreateINT4Vector(%d) failed.\n", fn, ffdata->numffts);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (aveNoise==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL8Vector(%d) failed.\n", fn, ffdata->numfprbins);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } 
   
   //Initialize to zero for far just at the start
   ihsfarstruct->ihsfar->data[0] = 0.0;
   REAL4 antweightsrms = 0.0;
   
   INT4 proberrcode = 0;
   ffdata->tfnormalization *= 0.5*inputParams->Tcoh;
   
   //Print message that we start the analysis
   fprintf(LOG, "Starting TwoSpect analysis...\n");
   fprintf(stderr, "Starting TwoSpect analysis...\n");
   
   
   REAL4Vector *antweightsforihs2h0 = XLALCreateREAL4Vector(ffdata->numffts);
   if (args_info.antennaOff_given) for (ii=0; ii<(INT4)antweightsforihs2h0->length; ii++) antweightsforihs2h0->data[ii] = 1.0;
   else {
      CompAntennaPatternWeights(antweightsforihs2h0, 0.0, 0.0, inputParams->searchstarttime, inputParams->Tcoh, inputParams->SFToverlap, inputParams->Tobs, inputParams->det);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: CompAntennaPatternWeights() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
   }

   
   
   //Search over the sky region
   while (scan.state != STATE_FINISHED) {
      fprintf(LOG, "Sky location: RA = %g, DEC = %g\n", dopplerpos.Alpha, dopplerpos.Delta);
      fprintf(stderr, "Sky location: RA = %g, DEC = %g\n", dopplerpos.Alpha, dopplerpos.Delta);
      
      //Determine detector velocity w.r.t. a sky location for each SFT
      CompAntennaVelocity(detectorVelocities, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, inputParams->searchstarttime, inputParams->Tcoh, inputParams->SFToverlap, inputParams->Tobs, inputParams->det, edat);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: CompAntennaVelocity() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      
      //Compute the bin shifts for each SFT
      CompBinShifts(binshifts, inputParams->fmin+0.5*inputParams->fspan, detectorVelocities, inputParams->Tcoh, inputParams->dopplerMultiplier);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: CompBinShifts() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      
      //Compute antenna pattern weights. If antennaOff input flag is given, then set all values equal to 1.0
      REAL4Vector *antweights = XLALCreateREAL4Vector((UINT4)ffdata->numffts);
      if (antweights==NULL) {
         fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ffdata->numffts);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      if (args_info.antennaOff_given) {
         for (ii=0; ii<(INT4)antweights->length; ii++) antweights->data[ii] = 1.0;
      } else {
         CompAntennaPatternWeights(antweights, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, inputParams->searchstarttime, inputParams->Tcoh, inputParams->SFToverlap, inputParams->Tobs, inputParams->det);
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: CompAntennaPatternWeights() failed.\n", fn);
            XLAL_ERROR(fn, XLAL_EFUNC);
         }
      }
      
      //Calculate antenna RMS value
      REAL4 currentAntWeightsRMS = calcRms(antweights);
      if (XLAL_IS_REAL4_FAIL_NAN(currentAntWeightsRMS)) {
         fprintf(stderr, "%s: calcRms() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      
      //Slide SFTs here -- need to slide the data and the estimated background
      REAL4Vector *TFdata_slided = XLALCreateREAL4Vector((UINT4)(ffdata->numffts*ffdata->numfbins));
      REAL4Vector *background_slided = XLALCreateREAL4Vector(TFdata_slided->length);
      if (TFdata_slided==NULL) {
         fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (UINT4)(ffdata->numffts*ffdata->numfbins));
         XLAL_ERROR(fn, XLAL_EFUNC);
      } else if (background_slided==NULL) {
         fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, TFdata_slided->length);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      slideTFdata(TFdata_slided, inputParams, usableTFdata, binshifts);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: slideTFdata() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      slideTFdata(background_slided, inputParams, background, binshifts);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: slideTFdata() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      //fprintf(stderr, "Mean of TFdata_slided %g, background_slided %g\n", calcMean(TFdata_slided), calcMean(background_slided));
      /* FILE *TFBACKGROUND = fopen("./output/tfbackground.dat","w");
      for (ii=0; ii<(INT4)background_slided->length; ii++) fprintf(TFBACKGROUND, "%.6g\n", background_slided->data[ii]);
      fclose(TFBACKGROUND); */
      
      //Check the RMS of the antenna weights, if bigger than standard deviation then reset the IHS FAR and the average noise background of the 2nd FFT
      if (antweightsrms == 0.0) antweightsrms = currentAntWeightsRMS;
      if ( fabs(currentAntWeightsRMS-antweightsrms)/antweightsrms >= 0.01 ) {
         ihsfarstruct->ihsfar->data[0] = 0.0;
         antweightsrms = currentAntWeightsRMS;
      }
      
      //TODO: Test normalization
      REAL8 skypointffnormalization = 1.0;
      ffPlaneNoise(aveNoise, inputParams, background_slided, antweightsforihs2h0, secondFFTplan, &(skypointffnormalization));
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: ffPlaneNoise() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      
      //Average noise floor of FF plane for each 1st FFT frequency bin
      ffdata->ffnormalization = 1.0;
      ffPlaneNoise(aveNoise, inputParams, background_slided, antweights, secondFFTplan, &(ffdata->ffnormalization));
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: ffPlaneNoise() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      
      //Compute the weighted TF data
      REAL4Vector *TFdata_weighted = XLALCreateREAL4Vector((UINT4)(ffdata->numffts*ffdata->numfbins));
      if (TFdata_weighted==NULL) {
         fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (UINT4)(ffdata->numffts*ffdata->numfbins));
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      tfWeight(TFdata_weighted, TFdata_slided, background_slided, antweights, inputParams);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: tfWeight() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      XLALDestroyREAL4Vector(TFdata_slided);
      XLALDestroyREAL4Vector(background_slided);
      XLALDestroyREAL4Vector(antweights);
      /* FILE *TFDATA = fopen("./output/tfdata.dat","w");
      for (jj=0; jj<(INT4)TFdata_weighted->length; jj++) fprintf(TFDATA,"%.6f\n",TFdata_weighted->data[jj]);
      fclose(TFDATA); */
      
      //Calculation of average TF noise per frequency bin ratio to total mean
      REAL4Vector *aveTFnoisePerFbinRatio = XLALCreateREAL4Vector(ffdata->numfbins);
      REAL4Vector *TSofPowers = XLALCreateREAL4Vector((UINT4)ffdata->numffts);
      if (aveTFnoisePerFbinRatio==NULL) {
         fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ffdata->numfbins);
         XLAL_ERROR(fn, XLAL_EFUNC);
      } else if (TSofPowers==NULL) {
         fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ffdata->numffts);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      for (ii=0; ii<ffdata->numfbins; ii++) {
         for (jj=0; jj<ffdata->numffts; jj++) TSofPowers->data[jj] = TFdata_weighted->data[jj*ffdata->numfbins + ii];
         aveTFnoisePerFbinRatio->data[ii] = calcRms(TSofPowers); //This approaches calcMean(TSofPowers) for stationary noise
      }
      REAL4 aveTFaveinv = 1.0/calcMean(aveTFnoisePerFbinRatio);
      for (ii=0; ii<ffdata->numfbins; ii++) aveTFnoisePerFbinRatio->data[ii] *= aveTFaveinv;
      XLALDestroyREAL4Vector(TSofPowers);
      
      
      //Do the second FFT
      makeSecondFFT(ffdata, TFdata_weighted, secondFFTplan);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: makeSecondFFT() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      
      REAL4 secFFTmean = calcMean(ffdata->ffdata);
      
      XLALDestroyREAL4Vector(TFdata_weighted);
      fprintf(stderr, "2nd FFT ave = %g, 2nd FFT stddev = %g, expected ave = %g\n", secFFTmean, calcStddev(ffdata->ffdata), calcMean(aveNoise)*calcMean(aveTFnoisePerFbinRatio));
      /* FILE *FFDATA = fopen(./output/ffdata.dat,"w");
      for (jj=0; jj<(INT4)ffdata->ffdata->length; jj++) fprintf(FFDATA,"%g\n",ffdata->ffdata->data[jj]);
      fclose(FFDATA); */
      
      if (secFFTmean==0.0) {
         fprintf(stderr, "Apparently, no SFTs were read in (Average power value is zero). Program exiting with failure.\n");
         XLAL_ERROR(fn, XLAL_FAILURE);
      }
      
      
////////Start of the IHS step!
      //Find the FAR of IHS sum
      if (ihsfarstruct->ihsfar->data[0]==0.0) {
         fprintf(stderr, "Determining IHS FAR values... ");
         fprintf(LOG, "Determining IHS FAR values... ");
         genIhsFar(ihsfarstruct, inputParams, maxrows, aveNoise);
         if (xlalErrno!=0) {
            fprintf(stderr,"\n%s: genIhsFar() failed.\n", fn);
            XLAL_ERROR(fn, XLAL_EFUNC);
         }
         fprintf(stderr, "done.\n");
         fprintf(LOG, "done.\n");
      }
      
      //Run the IHS algorithm on the data
      runIHS(ihsmaxima, ffdata, inputParams, maxrows, aveTFnoisePerFbinRatio);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: runIHS() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      
      //Find any IHS candidates
      findIHScandidates(ihsCandidates, ihsfarstruct, inputParams, ffdata, ihsmaxima, aveNoise, aveTFnoisePerFbinRatio);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: findIHScandidates() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      fprintf(LOG, "Candidates found in IHS step = %d\n", ihsCandidates->numofcandidates);
      fprintf(stderr, "Candidates found in IHS step = %d\n", ihsCandidates->numofcandidates);
      for (ii=0; ii<(INT4)ihsCandidates->numofcandidates; ii++) fprintf(stderr, "%d %g %g %g %g\n", ii, ihsCandidates->data[ii].fsig, ihsCandidates->data[ii].period, ihsCandidates->data[ii].moddepth, ihsCandidates->data[ii].h0);
////////End of the IHS step
      
////////Start of the Gaussian template search!
      if (!args_info.IHSonly_given) {
         
         //Test the IHS candidates against Gaussian templates in this function
         if ( testIHScandidates(gaussCandidates1, ihsCandidates, ffdata, aveNoise, aveTFnoisePerFbinRatio, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, inputParams) != 0 ) {
            fprintf(stderr, "%s: testIHScandidates() failed.\n", fn);
            XLAL_ERROR(fn, XLAL_EFUNC);
         }
         
         fprintf(LOG,"Initial stage done with candidates = %d\n",gaussCandidates1->numofcandidates);
         fprintf(stderr,"Initial stage done with candidates = %d\n",gaussCandidates1->numofcandidates);
         
         for (ii=0; ii<(INT4)gaussCandidates1->numofcandidates; ii++) fprintf(stderr, "Candidate %d: f0=%g, P=%g, df=%g\n", ii, gaussCandidates1->data[ii].fsig, gaussCandidates1->data[ii].period, gaussCandidates1->data[ii].moddepth);
      } /* if IHSonly is not given */
      else {
         if (exactCandidates2->length < exactCandidates2->numofcandidates+ihsCandidates->numofcandidates) {
            exactCandidates2 = resize_candidateVector(exactCandidates2, exactCandidates2->numofcandidates+ihsCandidates->numofcandidates);
            if (exactCandidates2->data==NULL) {
               fprintf(stderr,"%s: resize_candidateVector(%d) failed.\n", fn, ihsCandidates->numofcandidates);
               XLAL_ERROR(fn, XLAL_EFUNC);
            }
         }
         INT4 numofcandidatesalready = exactCandidates2->numofcandidates;
         for (ii=0; ii<(INT4)ihsCandidates->numofcandidates; ii++) {
            loadCandidateData(&(exactCandidates2->data[ii+numofcandidatesalready]), ihsCandidates->data[ii].fsig, ihsCandidates->data[ii].period, ihsCandidates->data[ii].moddepth, ihsCandidates->data[ii].ra, ihsCandidates->data[ii].dec, ihsCandidates->data[ii].stat, ihsCandidates->data[ii].h0, 0.0, 0, ihsCandidates->data[ii].normalization);
            exactCandidates2->data[ii+numofcandidatesalready].h0 /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25); //TODO: test moving this here
            (exactCandidates2->numofcandidates)++;
         }
      } /* if IHSonly is given */

////////End of the Gaussian template search

      //Reset IHS candidates, but keep length the same (doesn't reset actual values in the vector)
      ihsCandidates->numofcandidates = 0;
      
      //Search the IHS templates further if user has not specified IHSonly flag
      if (!args_info.IHSonly_given) {
////////Start clustering! Note that the clustering algorithm takes care of the period range of parameter space
         clusterCandidates(gaussCandidates2, gaussCandidates1, ffdata, inputParams, aveNoise, aveTFnoisePerFbinRatio, sftexist, 0);
         if (xlalErrno!=0) {
            fprintf(stderr,"%s: clusterCandidates() failed.\n", fn);
            XLAL_ERROR(fn, XLAL_EFUNC);
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
                  fprintf(stderr,"%s: resize_candidateVector(%d) failed.\n", fn, 2*gaussCandidates3->length);
                  XLAL_ERROR(fn, XLAL_EFUNC);
               }
            }
            //efficientTemplateSearch(&(gaussCandidates3->data[gaussCandidates3->numofcandidates]), gaussCandidates2->data[ii], gaussCandidates2->data[ii].fsig-2.5/inputParams->Tcoh, gaussCandidates2->data[ii].fsig+2.5/inputParams->Tcoh, 0.125/inputParams->Tcoh, 5, gaussCandidates2->data[ii].moddepth-2.5/inputParams->Tcoh, gaussCandidates2->data[ii].moddepth+2.5/inputParams->Tcoh, 0.125/inputParams->Tcoh, inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 0);
            bruteForceTemplateSearch(&(gaussCandidates3->data[gaussCandidates3->numofcandidates]), gaussCandidates2->data[ii], gaussCandidates2->data[ii].fsig-2.5/inputParams->Tcoh, gaussCandidates2->data[ii].fsig+2.5/inputParams->Tcoh, 11, 5, gaussCandidates2->data[ii].moddepth-2.5/inputParams->Tcoh, gaussCandidates2->data[ii].moddepth+2.5/inputParams->Tcoh, 11, inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 0);
            gaussCandidates3->numofcandidates++;
            
         } /* for ii < numofcandidates */
          
         for (ii=0; ii<(INT4)gaussCandidates3->numofcandidates; ii++) fprintf(stderr,"Candidate %d: f0=%g, P=%g, df=%g\n", ii, gaussCandidates3->data[ii].fsig, gaussCandidates3->data[ii].period, gaussCandidates3->data[ii].moddepth);
////////End detailed Gaussian template search

         //Reset 2nd round of Gaussian template candidates, but keep length the same (doesn't reset actual values in the vector)
         gaussCandidates2->numofcandidates = 0;

////////Start clustering!
         clusterCandidates(gaussCandidates4, gaussCandidates3, ffdata, inputParams, aveNoise, aveTFnoisePerFbinRatio, sftexist, 0);
         if (xlalErrno!=0) {
            fprintf(stderr,"%s: clusterCandidates() failed.\n", fn);
            XLAL_ERROR(fn, XLAL_EFUNC);
         }
         fprintf(LOG, "Clustering done with candidates = %d\n", gaussCandidates4->numofcandidates);
         fprintf(stderr, "Clustering done with candidates = %d\n", gaussCandidates4->numofcandidates);
         
         for (ii=0; ii<(INT4)gaussCandidates4->numofcandidates; ii++) fprintf(stderr, "Candidate %d: f0=%g, P=%g, df=%g\n", ii, gaussCandidates4->data[ii].fsig, gaussCandidates4->data[ii].period, gaussCandidates4->data[ii].moddepth);
////////End clustering
         
         //Reset 3rd set of Gaussian template candidates
         gaussCandidates3->numofcandidates = 0;

////////Initial check using "exact" template
         for (ii=0; ii<(INT4)gaussCandidates4->numofcandidates; ii++) {
            
            templateStruct *template = new_templateStruct(inputParams->templatelength);
            if (template==NULL) {
               fprintf(stderr,"%s: new_templateStruct(%d) failed.\n", fn, inputParams->templatelength);
               XLAL_ERROR(fn, XLAL_EFUNC); 
            }
            
            if (!args_info.gaussTemplatesOnly_given) {
               makeTemplate(template, gaussCandidates4->data[ii], inputParams, sftexist, secondFFTplan);
               if (xlalErrno!=0) {
                  fprintf(stderr,"%s: makeTemplate() failed.\n", fn);
                  XLAL_ERROR(fn, XLAL_EFUNC);
               }
            } else {
               makeTemplateGaussians(template, gaussCandidates4->data[ii], inputParams);
               if (xlalErrno!=0) {
                  fprintf(stderr,"%s: makeTemplateGaussians() failed.\n", fn);
                  XLAL_ERROR(fn, XLAL_EFUNC);
               }
            }
            
            farStruct *farval = new_farStruct();
            if (farval==NULL) {
               fprintf(stderr,"%s: new_farStruct() failed.\n", fn);
               XLAL_ERROR(fn, XLAL_EFUNC); 
            }
            
            if (inputParams->calcRthreshold) {
               numericFAR(farval, template, templatefarthresh, aveNoise, aveTFnoisePerFbinRatio, inputParams->rootFindingMethod);
               if (xlalErrno!=0) {
                  fprintf(stderr,"%s: numericFAR() failed.\n", fn);
                  XLAL_ERROR(fn, XLAL_EFUNC);
               }
            }
            
            REAL8 R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
            if (XLAL_IS_REAL8_FAIL_NAN(R)) {
               fprintf(stderr,"%s: calculateR() failed.\n", fn);
               XLAL_ERROR(fn, XLAL_EFUNC);
            }
            REAL8 prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, &proberrcode);
            if (XLAL_IS_REAL8_FAIL_NAN(prob)) {
               fprintf(stderr,"%s: probR() failed.\n", fn);
               XLAL_ERROR(fn, XLAL_EFUNC);
            }
            REAL8 h0 = 2.7426*pow(R/(inputParams->Tcoh*inputParams->Tobs),0.25);
            if ((!inputParams->calcRthreshold && prob<log10(templatefarthresh)) || (inputParams->calcRthreshold && R>farval->far)) {
               if (exactCandidates1->numofcandidates == exactCandidates1->length-1) {
                  exactCandidates1 = resize_candidateVector(exactCandidates1, 2*exactCandidates1->length);
                  if (exactCandidates1->data==NULL) {
                     fprintf(stderr,"%s: resize_candidateVector(%d) failed.\n", fn, 2*exactCandidates1->length);
                     XLAL_ERROR(fn, XLAL_EFUNC);
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
                  fprintf(stderr,"%s: resize_candidateVector(%d) failed.\n", fn, 2*exactCandidates2->length);
                  XLAL_ERROR(fn, XLAL_EFUNC);
               }
            }
            
            if (!args_info.gaussTemplatesOnly_given) {
               bruteForceTemplateSearch(&(exactCandidates2->data[exactCandidates2->numofcandidates]), exactCandidates1->data[ii], exactCandidates1->data[ii].fsig-1.0/inputParams->Tcoh, exactCandidates1->data[ii].fsig+1.0/inputParams->Tcoh, 5, 5, exactCandidates1->data[ii].moddepth-1.0/inputParams->Tcoh, exactCandidates1->data[ii].moddepth+1.0/inputParams->Tcoh, 5, inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 1);
            } else {
               bruteForceTemplateSearch(&(exactCandidates2->data[exactCandidates2->numofcandidates]), exactCandidates1->data[ii], exactCandidates1->data[ii].fsig-1.0/inputParams->Tcoh, exactCandidates1->data[ii].fsig+1.0/inputParams->Tcoh, 5, 5, exactCandidates1->data[ii].moddepth-1.0/inputParams->Tcoh, exactCandidates1->data[ii].moddepth+1.0/inputParams->Tcoh, 5, inputParams, ffdata->ffdata, sftexist, aveNoise, aveTFnoisePerFbinRatio, secondFFTplan, 0);
            }
            exactCandidates2->data[exactCandidates2->numofcandidates].h0 /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25);  //TODO: test moving this here
            exactCandidates2->numofcandidates++;
            
            fprintf(stderr,"Candidate %d: f0=%g, P=%g, df=%g\n", ii, exactCandidates2->data[ii].fsig, exactCandidates2->data[ii].period, exactCandidates2->data[ii].moddepth);
         } /* for ii < numofcandidates */
////////End of detailed search
         
         //Reset first round of exact template candidates, but keep length the same (doesn't reset actual values in the vector)
         exactCandidates1->numofcandidates = 0;
         
         fprintf(LOG,"Exact step is done with the total number of candidates = %d\n", exactCandidates2->numofcandidates);
         fprintf(stderr,"Exact step is done with the total number of candidates = %d\n", exactCandidates2->numofcandidates);
         
      } /* if IHSonly is not given */
      
      //Determine upper limits
      upperlimits->data[upperlimits->length-1].alpha = (REAL4)dopplerpos.Alpha;
      upperlimits->data[upperlimits->length-1].delta = (REAL4)dopplerpos.Delta;
      upperlimits->data[upperlimits->length-1].normalization = ffdata->tfnormalization;
      //skypoint95UL(&(upperlimits->data[upperlimits->length-1]), ihsfarstruct, inputParams, ffdata, ihsmaxima, aveNoise, aveTFnoisePerFbinRatio);
      skypoint95UL(&(upperlimits->data[upperlimits->length-1]), inputParams, ffdata, ihsmaxima, aveNoise, aveTFnoisePerFbinRatio);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: skypoint95UL() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      //for (ii=0; ii<(INT4)upperlimits->data[upperlimits->length-1].ULval->length; ii++) upperlimits->data[upperlimits->length-1].ULval->data[ii] /= sqrt(ffdata->tfnormalization);
      //for (ii=0; ii<(INT4)upperlimits->data[upperlimits->length-1].ULval->length; ii++) upperlimits->data[upperlimits->length-1].ULval->data[ii] /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete,0.25);   //TODO: verify the sft-loss parameter correction
      for (ii=0; ii<(INT4)upperlimits->data[upperlimits->length-1].ULval->length; ii++) upperlimits->data[upperlimits->length-1].ULval->data[ii] /= sqrt(ffdata->tfnormalization)*pow(frac_tobs_complete*ffdata->ffnormalization/skypointffnormalization,0.25);   //TODO: verify the sft-loss parameter correction and sky position normalization difference
      upperlimits = resize_UpperLimitVector(upperlimits, upperlimits->length+1);
      if (upperlimits->data==NULL) {
         fprintf(stderr,"%s: resize_UpperLimitVector(%d) failed.\n", fn, upperlimits->length+1);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      
      //Destroy stuff
      XLALDestroyREAL4Vector(aveTFnoisePerFbinRatio);
      
      //Iterate to next sky location
      if ((XLALNextDopplerSkyPos(&dopplerpos, &scan))!=0) {
         fprintf(stderr,"%s: XLALNextDopplerSkyPos() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
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
   
   ULFILE = fopen(t,"w");
   if (ULFILE==NULL) {
      fprintf(stderr, "%s: UL file could not be opened.\n", fn);
      XLAL_ERROR(fn, XLAL_EINVAL);
   }
   for (ii=0; ii<(INT4)upperlimits->length-1; ii++) outputUpperLimitToFile(ULFILE, upperlimits->data[ii], inputParams->ULmindf, inputParams->ULmaxdf, inputParams->printAllULvalues);
   fclose(ULFILE);
   
   //Destroy varaibles
   XLALDestroyREAL4Vector(antweightsforihs2h0);
   XLALDestroyREAL4Vector(background);
   XLALDestroyREAL4Vector(usableTFdata);
   XLALDestroyREAL4Vector(detectorVelocities);
   XLALDestroyREAL4Vector(aveNoise);
   XLALDestroyINT4Vector(binshifts);
   XLALDestroyINT4Vector(sftexist);
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
   
   //print end time
   time(&programendtime);
   ptm = localtime(&programendtime);
   fprintf(stderr, "Program finished on %s\n", asctime(ptm));
   fprintf(LOG, "Program finished on %s\n", asctime(ptm));
   
   fclose(LOG);
   
   LALCheckMemoryLeaks();
   
   return 0;

} /* main() */



//////////////////////////////////////////////////////////////
// Create new inputParamsStruct  -- done
inputParamsStruct * new_inputParams(void)
{
   
   const CHAR *fn = __func__;
   
   inputParamsStruct *input = XLALMalloc(sizeof(*input));
   if (input==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%zu) failed.", fn, sizeof(*input));
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
   }
   
   return input;

} /* new_inputParams() */


//////////////////////////////////////////////////////////////
// Destroy inputParamsStruct  -- done
void free_inputParams(inputParamsStruct *input)
{
   
   XLALFree((CHAR*)input->sftType);
   XLALFree((inputParamsStruct*)input);

} /* free_inputParams() */



//////////////////////////////////////////////////////////////
// Allocate ffdataStruct vectors  -- done
ffdataStruct * new_ffdata(inputParamsStruct *input)
{
   
   const CHAR *fn = __func__;
   
   ffdataStruct *ffdata = XLALMalloc(sizeof(*ffdata));
   if (ffdata==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%zu) failed.\n", fn, sizeof(*ffdata));
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
   }
   
   ffdata->numfbins = (INT4)(round(input->fspan*input->Tcoh)+1);
   ffdata->numffts = (INT4)floor(input->Tobs/(input->Tcoh-input->SFToverlap)-1);
   ffdata->numfprbins = (INT4)floorf(ffdata->numffts*0.5) + 1;
   
   ffdata->ffdata = XLALCreateREAL4Vector(ffdata->numfbins * ffdata->numfprbins);
   if (ffdata->ffdata==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ffdata->numfbins * ffdata->numfprbins);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
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
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj;
   LALStatus status;
   status.statusPtr = NULL;
   SFTCatalog *catalog = NULL;
   SFTVector *sfts = NULL;
   
   LIGOTimeGPS start = LIGOTIMEGPSZERO, end = LIGOTIMEGPSZERO;
   XLALGPSSetREAL8(&start, input->searchstarttime);
   if (xlalErrno != 0) {
      fprintf(stderr, "%s: XLALGPSSetREAL8() failed on start time = %.9f.\n", fn, input->searchstarttime);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   }
   XLALGPSSetREAL8(&end, input->searchstarttime+input->Tobs);
   if (xlalErrno != 0) {
      fprintf(stderr, "%s: XLALGPSSetREAL8() failed on end time = %.9f.\n", fn, input->searchstarttime+input->Tobs);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   }
   
   SFTConstraints constraints;
   constraints.detector = NULL;
   constraints.startTime = constraints.endTime = NULL;
   constraints.timestamps = NULL;
   constraints.detector = input->det.frDetector.prefix;
   constraints.startTime = &start;
   constraints.endTime = &end;
   
   //Find SFT files
   LALSFTdataFind(&status, &catalog, sft_dir, &constraints);
   if (status.statusCode != 0) {
      fprintf(stderr,"%s: LALSFTdataFind() failed with code = %d.\n", fn, status.statusCode);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   }
   
   //Determine band size (remember to get extra bins because of the running median and the bin shifts due to detector velocity)
   REAL8 minfbin = round(input->fmin*input->Tcoh - 0.5*(input->blksize-1) - (input->maxbinshift))/input->Tcoh;
   REAL8 maxfbin = round((input->fmin + input->fspan)*input->Tcoh + 0.5*(input->blksize-1) + (input->maxbinshift))/input->Tcoh;
   
   //Now extract the data
   LALLoadSFTs(&status, &sfts, catalog, minfbin, maxfbin);
   if (status.statusCode != 0) {
      fprintf(stderr,"%s: LALLoadSFTs() failed with code = %d.\n", fn, status.statusCode);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if (sfts == NULL) {
      fprintf(stderr,"%s: LALLoadSFTs() failed to load SFTs with given input parameters.\n", fn);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   }
   
   //Now put the power data into the TF plane, looping through each SFT
   //If an SFT doesn't exit, fill the TF pixels of the SFT with zeros
   //INT4 numffts = (INT4)floor(2*(input->Tobs/input->Tcoh)-1);
   INT4 numffts = (INT4)floor(input->Tobs/(input->Tcoh-input->SFToverlap)-1);
   INT4 sftlength;
   if (sfts->length == 0) sftlength = (INT4)(maxfbin*input->Tcoh - minfbin*input->Tcoh + 1);
   else sftlength = sfts->data->data->length;
   INT4 nonexistantsft = 0;
   REAL4Vector *tfdata = XLALCreateREAL4Vector((UINT4)(numffts*sftlength));
   if (tfdata==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numffts*sftlength);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   }
   
   REAL8 sqrtnorm = sqrt(*(normalization));
   for (ii=0; ii<numffts; ii++) {
      
      SFTDescriptor *sftdescription = &(catalog->data[ii - nonexistantsft]);
      if (sftdescription->header.epoch.gpsSeconds == (INT4)round(ii*(input->Tcoh-input->SFToverlap)+input->searchstarttime)) {
         SFTtype *sft = &(sfts->data[ii - nonexistantsft]);
         for (jj=0; jj<sftlength; jj++) {
            COMPLEX8 sftcoeff = sft->data->data[jj];
            tfdata->data[ii*sftlength + jj] = (REAL4)((sqrtnorm*sftcoeff.re)*(sqrtnorm*sftcoeff.re) + (sqrtnorm*sftcoeff.im)*(sqrtnorm*sftcoeff.im));  //power, normalized
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
   
   LALDestroySFTCatalog(&status, &catalog);
   XLALDestroySFTVector(sfts);
   
   fprintf(stderr,"TF before weighting, mean subtraction = %g\n",calcMean(tfdata));
   
   return tfdata;

} /* readInSFTs() */




//////////////////////////////////////////////////////////////
// Slide SFT TF data  -- 
void slideTFdata(REAL4Vector *output, inputParamsStruct *input, REAL4Vector *tfdata, INT4Vector *binshifts)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj;
   INT4 numffts = (INT4)floor(input->Tobs/(input->Tcoh-input->SFToverlap)-1);
   INT4 numfbins = (INT4)(round(input->fspan*input->Tcoh)+1);
   
   for (ii=0; ii<numffts; ii++) {
      if (binshifts->data[ii]>input->maxbinshift) {
         fprintf(stderr, "%s: SFT slide value %d is greater than maximum value predicted (%d)", fn, binshifts->data[ii], input->maxbinshift);
         XLAL_ERROR_VOID(fn, XLAL_EFAILED);
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
   
   const CHAR *fn = __func__;
   
   LALStatus status;
   status.statusPtr = NULL;
   REAL8 bias;
   INT4 ii, jj;
   INT4 totalfbins = numfbins + blksize - 1;
   
   //Blocksize of running median
   LALRunningMedianPar block = {blksize};
   
   //Running median bias calculation
   if (blksize<1000) {
      LALRngMedBias(&status, &bias, blksize);
      if (status.statusCode != 0) {
         fprintf(stderr,"%s: LALRngMedBias() failed.\n", fn);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
   } else {
      bias = LAL_LN2;
   }
   REAL8 invbias = 1.0/(bias*1.0099993480677538);  //TODO: StackSlide normalization for 101 bins
   
   REAL4Vector *inpsd = XLALCreateREAL4Vector(totalfbins);
   REAL4Vector *mediansout = XLALCreateREAL4Vector(numfbins);
   if (inpsd==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, totalfbins);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (mediansout==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numfbins);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   for (ii=0; ii<numffts; ii++) {
      //If the SFT values were not zero, then compute the running median
      if (tfdata->data[ii*totalfbins]!=0.0) {
         //Determine running median value, convert to mean value
         memcpy(inpsd->data, &(tfdata->data[ii*inpsd->length]), sizeof(REAL4)*inpsd->length);
         
         //calculate running median
         LALSRunningMedian2(&status, mediansout, inpsd, block);
         if (status.statusCode != 0) {
            fprintf(stderr,"%s: LALSRunningMedian2() failed.\n", fn);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         }
         
         //Now make the output medians into means by multiplying by 1/bias
         for (jj=0; jj<(INT4)mediansout->length; jj++) output->data[ii*numfbins + jj] = (REAL4)(mediansout->data[jj]*invbias);
      } else {
         //Otherwise, set means to zero
         for (jj=0; jj<(INT4)mediansout->length; jj++) output->data[ii*numfbins + jj] = 0.0;
      }
   } /* for ii < numffts */
   
   fprintf(stderr,"Mean of running means = %g\n",calcMean(output));
   
   XLALDestroyREAL4Vector(inpsd);
   XLALDestroyREAL4Vector(mediansout);

} /* tfRngMeans() */



/* Critical values of KS test (from Bickel and Doksum). Does not apply directly (mean determined from distribution)
alpha=0.01
n       10      20      30      40      50      60      80      n>80
        .489    .352    .290    .252    .226    .207    .179    1.628/(sqrt(n)+0.12+0.11/sqrt(n))
        
alpha=0.05
n       10      20      30      40      50      60      80      n>80
        .409    .294    .242    .210    .188    .172    .150    1.358/(sqrt(n)+0.12+0.11/sqrt(n))
*/
INT4Vector * markBadSFTs(REAL4Vector *tfdata, inputParamsStruct *params)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii;
   
   INT4 numffts = (INT4)floor(params->Tobs/(params->Tcoh-params->SFToverlap)-1);    //Number of FFTs
   INT4 numfbins = (INT4)(round(params->fspan*params->Tcoh)+1)+2*params->maxbinshift+params->blksize-1;     //Number of frequency bins
   
   INT4Vector *output = XLALCreateINT4Vector(numffts);
   if (output==NULL) {
      fprintf(stderr, "%s: XLALCreateINT4Vector(%d) failed.\n", fn, numffts);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   }
   for (ii=0; ii<numffts; ii++) output->data[ii] = 0;
   REAL4Vector *tempvect = XLALCreateREAL4Vector(numfbins);
   if (tempvect==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numfbins);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   }
   
   REAL8 ksthreshold = 1.358/(sqrt(numfbins)+0.12+0.11/sqrt(numfbins));
   for (ii=0; ii<numffts; ii++) {
      if (tfdata->data[ii*numfbins]!=0.0) {
         memcpy(tempvect->data, &(tfdata->data[ii*numfbins]), sizeof(REAL4)*tempvect->length);
         REAL8 kstest = ks_test_exp(tempvect);
         if (XLAL_IS_REAL8_FAIL_NAN(kstest)) {
            fprintf(stderr,"%s: ks_test_exp() failed.\n", fn);
            XLAL_ERROR_NULL(fn, XLAL_EFUNC);
         }
         
         if (kstest>ksthreshold) output->data[ii] = 1;
      }
   }
   
   XLALDestroyREAL4Vector(tempvect);
   
   return output;
   
}



void removeBadSFTs(REAL4Vector *tfdata, INT4Vector *badsfts)
{
   
   INT4 ii, jj;
   
   INT4 numfbins_tfdata = tfdata->length/badsfts->length;
   
   for (ii=0; ii<(INT4)badsfts->length; ii++) {
      if (badsfts->data[ii]==1) {
         for (jj=0; jj<numfbins_tfdata; jj++) tfdata->data[ii*numfbins_tfdata + jj] = 0.0;
      }
   }
   
}

/* Modifies input vector!!! */
void tfMeanSubtract(REAL4Vector *tfdata, REAL4Vector *rngMeans, INT4 numffts, INT4 numfbins)
{
   
   INT4 ii, jj;
   
   for (ii=0; ii<numffts; ii++) if (rngMeans->data[ii*numfbins]!=0.0) for (jj=0; jj<numfbins; jj++) tfdata->data[ii*numfbins+jj] -= rngMeans->data[ii*numfbins+jj];
   
} /* tfMeanSubtract() */

void tfWeight(REAL4Vector *output, REAL4Vector *tfdata, REAL4Vector *rngMeans, REAL4Vector *antPatternWeights, inputParamsStruct *input)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj;
   
   INT4 numffts = (INT4)floor(input->Tobs/(input->Tcoh-input->SFToverlap)-1);    //Number of FFTs
   INT4 numfbins = (INT4)(round(input->fspan*input->Tcoh)+1);     //Number of frequency bins
   
   REAL4Vector *antweightssq = XLALCreateREAL4Vector(numffts);
   REAL4Vector *rngMeanssq = XLALCreateREAL4Vector(numffts);
   if (antweightssq==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numffts);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (rngMeanssq==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numffts);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   for (ii=0; ii<(INT4)rngMeanssq->length; ii++) rngMeanssq->data[ii] = 0.0;
   
   //for (ii=0; ii<numffts; ii++) antweightssq->data[ii] = antPatternWeights->data[ii]*antPatternWeights->data[ii];
   antweightssq = XLALSSVectorMultiply(antweightssq, antPatternWeights, antPatternWeights);
   //antweightssq = SSVectorMultiply_with_stride_and_offset(antweightssq, antPatternWeights, antPatternWeights, 1, 1, 0, 0);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: SSVectorMutiply_with_stride_and_offset() failed.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   for (ii=0; ii<numfbins; ii++) {
      
      rngMeanssq = SSVectorMultiply_with_stride_and_offset(rngMeanssq, rngMeans, rngMeans, numfbins, numfbins, ii, ii);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: SSVectorMutiply_with_stride_and_offset() failed.\n", fn);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      
      //If noiseWeightOff is given, then set all the noise weights to be 1.0
      if (input->noiseWeightOff!=0) for (jj=0; jj<(INT4)rngMeanssq->length; jj++) if (rngMeanssq->data[jj]!=0.0) rngMeanssq->data[jj] = 1.0;
      
      //Get sum of antenna pattern weight/variances for each frequency bin as a function of time (only for existant SFTs)
      REAL8 sumofweights = 0.0;
      for (jj=0; jj<numffts; jj++) if (rngMeanssq->data[jj] != 0.0) sumofweights += antweightssq->data[jj]/rngMeanssq->data[jj];
      REAL8 invsumofweights = 1.0/sumofweights;
      
      //Now do noise weighting, antenna pattern weighting
      for (jj=0; jj<numffts; jj++) {
         if (rngMeanssq->data[jj] != 0.0) output->data[jj*numfbins+ii] = (REAL4)(invsumofweights*antPatternWeights->data[jj]*tfdata->data[jj*numfbins+ii]/rngMeanssq->data[jj]);
         else output->data[jj*numfbins+ii] = 0.0;
      } /* for jj < numffts */
   } /* for ii < numfbins */
   
   
   XLALDestroyREAL4Vector(antweightssq);
   XLALDestroyREAL4Vector(rngMeanssq);
   
   //fprintf(stderr,"TF after weighting, mean subtraction = %g\n",calcMean(output));   
   
} /* tfWeight() */



//////////////////////////////////////////////////////////////
// Do the weighting by noise variance (from tfRngMeans), mean subtraction, and antenna pattern weights  -- done
void tfWeightMeanSubtract(REAL4Vector *output, REAL4Vector *tfdata, REAL4Vector *rngMeans, REAL4Vector *antPatternWeights, inputParamsStruct *input)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj;
    
   INT4 numffts = (INT4)floor(input->Tobs/(input->Tcoh-input->SFToverlap)-1);    //Number of FFTs
   INT4 numfbins = (INT4)(round(input->fspan*input->Tcoh)+1);     //Number of frequency bins
   
   REAL4Vector *antweightssq = XLALCreateREAL4Vector(numffts);
   REAL4Vector *rngMeanssq = XLALCreateREAL4Vector(numffts);
   if (antweightssq==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numffts);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (rngMeanssq==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numffts);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   for (ii=0; ii<(INT4)rngMeanssq->length; ii++) rngMeanssq->data[ii] = 0.0;
   
   //for (ii=0; ii<numffts; ii++) antweightssq->data[ii] = antPatternWeights->data[ii]*antPatternWeights->data[ii];
   //antweightssq = XLALSSVectorMultiply(antweightssq, antPatternWeights, antPatternWeights);
   antweightssq = SSVectorMultiply_with_stride_and_offset(antweightssq, antPatternWeights, antPatternWeights, 1, 1, 0, 0);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: SSVectorMutiply_with_stride_and_offset() failed.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   for (ii=0; ii<numfbins; ii++) {
      
      rngMeanssq = SSVectorMultiply_with_stride_and_offset(rngMeanssq, rngMeans, rngMeans, numfbins, numfbins, ii, ii);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: SSVectorMutiply_with_stride_and_offset() failed.\n", fn);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
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
   
   
   XLALDestroyREAL4Vector(antweightssq);
   XLALDestroyREAL4Vector(rngMeanssq);
   
   //fprintf(stderr,"TF after weighting, mean subtraction = %g\n",calcMean(output));

} /* tfWeightMeanSubtract() */


//////////////////////////////////////////////////////////////
// Make the second FFT powers
void makeSecondFFT(ffdataStruct *output, REAL4Vector *tfdata, REAL4FFTPlan *plan)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj;
   
   //Do the second FFT
   REAL4Vector *x = XLALCreateREAL4Vector(output->numffts);
   REAL4Window *win = XLALCreateHannREAL4Window(x->length);
   REAL4Vector *psd = XLALCreateREAL4Vector((UINT4)floor(x->length*0.5)+1);
   if (x==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, output->numffts);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (win==NULL) {
      fprintf(stderr,"%s: XLALCreateHannREAL4Window(%d) failed.\n", fn, x->length);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (psd==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (UINT4)floor(x->length*0.5)+1);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   REAL8 winFactor = 8.0/3.0;
   REAL8 psdfactor = winFactor;
   //First loop over frequencies
   for (ii=0; ii<output->numfbins; ii++) {
   
      //Next, loop over times and pick the right frequency bin for each FFT and window
      //for (jj=0; jj<(INT4)x->length; jj++) x->data[jj] = (tfdata->data[ii + jj*numfbins]*win->data->data[jj]);
      x = SSVectorMultiply_with_stride_and_offset(x, tfdata, win->data, output->numfbins, 1, ii, 0);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: SSVectorMutiply_with_stride_and_offset() failed.\n", fn);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      
      //Make the FFT
      if ( (XLALREAL4PowerSpectrum(psd, x, plan)) != 0) {
         fprintf(stderr,"%s: XLALREAL4PowerSpectrum() failed.\n", fn);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
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
      for (jj=0; jj<(INT4)psd->length; jj++) output->ffdata->data[psd->length*ii + jj] = (REAL4)(psd->data[jj]*psdfactor*output->ffnormalization);
      
   } /* for ii < numfbins */
   
   XLALDestroyREAL4Vector(x);
   XLALDestroyREAL4Vector(psd);
   XLALDestroyREAL4Window(win);
   
} /* makeSecondFFT() */



//////////////////////////////////////////////////////////////
// Determine the average of the noise power in each frequency bin across the band
REAL4 avgTFdataBand(REAL4Vector *backgrnd, INT4 numfbins, INT4 numffts, INT4 binmin, INT4 binmax)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii;
   REAL4Vector *aveNoiseInTime = XLALCreateREAL4Vector((UINT4)numffts);
   REAL4Vector *rngMeansOverBand = XLALCreateREAL4Vector((UINT4)(binmax-binmin));
   if (aveNoiseInTime==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numffts);
      XLAL_ERROR_REAL4(fn, XLAL_EFUNC);
   } else if (rngMeansOverBand==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (UINT4)(binmax-binmin));
      XLAL_ERROR_REAL4(fn, XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)aveNoiseInTime->length; ii++) {
      memcpy(rngMeansOverBand->data, &(backgrnd->data[ii*numfbins + binmin]), sizeof(*rngMeansOverBand->data)*rngMeansOverBand->length);
      aveNoiseInTime->data[ii] = calcMean(rngMeansOverBand);
      if (XLAL_IS_REAL4_FAIL_NAN(aveNoiseInTime->data[ii])) {
         fprintf(stderr,"%s: calcMean() failed.\n", fn);
         XLAL_ERROR_REAL4(fn, XLAL_EFUNC);
      }
   } /* for ii < aveNoiseInTime->length */
   
   REAL4 avgTFdata = calcMean(aveNoiseInTime);
   if (XLAL_IS_REAL4_FAIL_NAN(avgTFdata)) {
      fprintf(stderr,"%s: calcMean() failed.\n", fn);
      XLAL_ERROR_REAL4(fn, XLAL_EFUNC);
   }
   
   XLALDestroyREAL4Vector(aveNoiseInTime);
   XLALDestroyREAL4Vector(rngMeansOverBand);
   
   return avgTFdata;
   
} /* avgTFdataBand() */


//////////////////////////////////////////////////////////////
// Determine the rms of the noise power in each frequency bin across the band
REAL4 rmsTFdataBand(REAL4Vector *backgrnd, INT4 numfbins, INT4 numffts, INT4 binmin, INT4 binmax)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii;
   REAL4Vector *aveNoiseInTime = XLALCreateREAL4Vector((UINT4)numffts);
   REAL4Vector *rngMeansOverBand = XLALCreateREAL4Vector((UINT4)(binmax-binmin));
   if (aveNoiseInTime==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numffts);
      XLAL_ERROR_REAL4(fn, XLAL_EFUNC);
   } else if (rngMeansOverBand==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (UINT4)(binmax-binmin));
      XLAL_ERROR_REAL4(fn, XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)aveNoiseInTime->length; ii++) {
      memcpy(rngMeansOverBand->data, &(backgrnd->data[ii*numfbins + binmin]), sizeof(*rngMeansOverBand->data)*rngMeansOverBand->length);
      aveNoiseInTime->data[ii] = calcMean(rngMeansOverBand);
      if (XLAL_IS_REAL4_FAIL_NAN(aveNoiseInTime->data[ii])) {
         fprintf(stderr,"%s: calcMean() failed.\n", fn);
         XLAL_ERROR_REAL4(fn, XLAL_EFUNC);
      }
   } /* for ii < aveNoiseInTime->length */
   
   REAL4 rmsTFdata = calcRms(aveNoiseInTime);
   if (XLAL_IS_REAL4_FAIL_NAN(rmsTFdata)) {
      fprintf(stderr,"%s: calcRms() failed.\n", fn);
      XLAL_ERROR_REAL4(fn, XLAL_EFUNC);
   }
   
   XLALDestroyREAL4Vector(aveNoiseInTime);
   XLALDestroyREAL4Vector(rngMeansOverBand);
   
   return rmsTFdata;
   
} /* rmsTFdataBand() */


//////////////////////////////////////////////////////////////
// Measure of the average noise power in each 2st FFT frequency bin  -- 
void ffPlaneNoise(REAL4Vector *aveNoise, inputParamsStruct *input, REAL4Vector *backgrnd, REAL4Vector *antweights, REAL4FFTPlan *plan, REAL8 *normalization)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj, numfbins, numffts, numfprbins;
   REAL8 invsumofweights = 0.0;
   REAL8 sumofweights = 0.0;
   
   numfbins = (INT4)(round(input->fspan*input->Tcoh)+1);   //Number of frequency bins
   numffts = (INT4)floor(input->Tobs/(input->Tcoh-input->SFToverlap)-1); //Number of FFTs
   numfprbins = (INT4)floor(numffts*0.5)+1;     //number of 2nd fft frequency bins
   
   //Initialize the random number generator
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
   if (rng==NULL) {
      fprintf(stderr,"%s: gsl_rng_alloc() failed.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   /* srand(time(NULL));
   UINT8 randseed = rand();
   gsl_rng_set(rng, randseed); */
   gsl_rng_set(rng, 0);
   
   //Set up for making the PSD
   for (ii=0; ii<(INT4)aveNoise->length; ii++) aveNoise->data[ii] = 0.0;
   REAL4Window *win = XLALCreateHannREAL4Window((UINT4)numffts);
   REAL4Vector *psd = XLALCreateREAL4Vector((UINT4)numfprbins);   //Current PSD calculation
   if (win==NULL) {
      fprintf(stderr,"%s: XLALCreateHannREAL4Window(%d) failed.\n", fn, numffts);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (psd==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (UINT4)numfprbins);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   REAL4 winFactor = 8.0/3.0;

   //Average each SFT across the frequency band, also compute normalization factor
   REAL4Vector *aveNoiseInTime = XLALCreateREAL4Vector((UINT4)numffts);
   REAL4Vector *rngMeansOverBand = XLALCreateREAL4Vector((UINT4)numfbins);
   if (aveNoiseInTime==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numffts);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (rngMeansOverBand==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numfbins);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)aveNoiseInTime->length; ii++) {
      
      if (backgrnd->data[ii*numfbins]!=0.0) {
         memcpy(rngMeansOverBand->data, &(backgrnd->data[ii*numfbins]), sizeof(*rngMeansOverBand->data)*rngMeansOverBand->length);
         //aveNoiseInTime->data[ii] = calcMean(rngMeansOverBand);
         aveNoiseInTime->data[ii] = (REAL4)(calcRms(rngMeansOverBand));  //For exp dist and large blksize this approaches the mean
         if (XLAL_IS_REAL4_FAIL_NAN(aveNoiseInTime->data[ii])) {
            fprintf(stderr,"%s: calcMean() failed.\n", fn);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         }
         
         if (input->noiseWeightOff==0) sumofweights += (antweights->data[ii]*antweights->data[ii])/(aveNoiseInTime->data[ii]*aveNoiseInTime->data[ii]);
         else sumofweights += (antweights->data[ii]*antweights->data[ii]);
      } else {
         aveNoiseInTime->data[ii] = 0.0;
      }
      
   } /* for ii < aveNoiseInTime->length */
   invsumofweights = 1.0/sumofweights;
   
   //Load time series of powers, normalize, mean subtract and Hann window
   REAL4Vector *x = XLALCreateREAL4Vector(aveNoiseInTime->length);
   REAL8Vector *multiplicativeFactor = XLALCreateREAL8Vector(x->length);
   REAL4 psdfactor = winFactor;
   
   for (ii=0; ii<(INT4)x->length; ii++) {
      if (aveNoiseInTime->data[ii] != 0.0) {
         if (input->noiseWeightOff==0) multiplicativeFactor->data[ii] = win->data->data[ii]*antweights->data[ii]/aveNoiseInTime->data[ii]*invsumofweights;
         else multiplicativeFactor->data[ii] = win->data->data[ii]*antweights->data[ii]*invsumofweights;
      } else {
         multiplicativeFactor->data[ii] = 0.0;
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
   //FILE *BACKGRNDX = fopen("./output/backgroundx.dat","w");
   for (ii=0; ii<400; ii++) {
      for (jj=0; jj<(INT4)x->length; jj++) {
         if (aveNoiseInTime->data[jj] != 0.0) {
            noiseval = expRandNum(aveNoiseInTime->data[jj], rng);
            if (jj==0 || (jj>0 && aveNoiseInTime->data[jj-1] == 0.0)) {
               x->data[jj] = (REAL4)(multiplicativeFactor->data[jj]*(noiseval/aveNoiseInTime->data[jj]-1.0));
               //fprintf(BACKGRNDX, "%f\n", (noiseval/aveNoiseInTime->data[jj]-1.0)*antweights->data[jj]/aveNoiseInTime->data[jj]*invsumofweights);
            } else {
               REAL8 newnoiseval = (1.0-corrfactorsquared)*noiseval + corrfactorsquared*prevnoiseval;
               REAL8 newavenoise = (1.0-corrfactorsquared)*aveNoiseInTime->data[jj] + corrfactorsquared*aveNoiseInTime->data[jj-1];
               x->data[jj] = (REAL4)(multiplicativeFactor->data[jj]*(newnoiseval/newavenoise-1.0));
               //fprintf(BACKGRNDX, "%f\n", (newnoiseval/newavenoise-1.0)*antweights->data[jj]/aveNoiseInTime->data[jj]*invsumofweights);
            }
            prevnoiseval = noiseval;
         } else {
            x->data[jj] = 0.0;
            //fprintf(BACKGRNDX, "%f\n", 0.0);
         }
      } /* for jj < x->length */
      
      //Do the FFT
      if ( (XLALREAL4PowerSpectrum(psd, x, plan)) != 0) {
         fprintf(stderr,"%s: XLALREAL4PowerSpectrum() failed.\n", fn);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      
      //Rescale and sum into the bins
      for (jj=0; jj<(INT4)aveNoise->length; jj++) {
         aveNoise->data[jj] += psd->data[jj];
      }
   } /* for ii < 200 */
   //fclose(BACKGRNDX);
   
   //Average
   for (ii=0; ii<(INT4)aveNoise->length; ii++) aveNoise->data[ii] *= (REAL4)(2.5e-3*psdfactor*(1.0+2.0*corrfactorsquared));
   
   //Fix 0th and end bins (0 only for odd x->length, 0 and end for even x->length)
   if (GSL_IS_EVEN(x->length)==1) {
      aveNoise->data[0] *= 2.0;
      aveNoise->data[aveNoise->length-1] *= 2.0;
   } else {
      aveNoise->data[0] *= 2.0;
   }
   
   //FILE *BACKGRND = fopen("./output/background.dat","w");
   
   //Compute normalization
   *(normalization) = 1.0/(calcMean(aveNoise));
   for (ii=0; ii<(INT4)aveNoise->length; ii++) {
      aveNoise->data[ii] *= *(normalization);
      //fprintf(BACKGRND,"%.6f\n",aveNoise->data[ii]);
   }
   
   //TODO: remove this extra factor
   // *(normalization) /= 1.0245545525190294;
   *normalization /= 1.040916688722758;
   
   //fclose(BACKGRND);

   XLALDestroyREAL4Vector(x);
   XLALDestroyREAL4Vector(psd);
   XLALDestroyREAL4Window(win);
   XLALDestroyREAL4Vector(aveNoiseInTime);
   XLALDestroyREAL4Vector(rngMeansOverBand);
   XLALDestroyREAL8Vector(multiplicativeFactor);
   gsl_rng_free(rng);

} /* ffPlaneNoise() */


//For testing purposes only!!!!
REAL4Vector * simpleTFdata(REAL8 fsig, REAL8 period, REAL8 moddepth, REAL8 Tcoh, REAL8 Tobs, REAL8 SFToverlap, REAL8 fminimum, REAL8 fmaximum, REAL8 sqrtSh)
{
   
   const CHAR *fn = __func__;
   
   INT4 numfbins = (INT4)(round((fmaximum-fminimum)*Tcoh)+1);   //Number of frequency bins
   INT4 numffts = (INT4)floor(Tobs/(Tcoh-SFToverlap)-1); //Number of FFTs
   
   REAL4Vector *output = XLALCreateREAL4Vector(numfbins*numffts);
   
   //Initialize the random number generator
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
   if (rng==NULL) {
      fprintf(stderr,"%s: gsl_rng_alloc() failed.\n", fn);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
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
   
   gsl_rng_free(rng);
   
   return output;
   
}




INT4 readTwoSpectInputParams(inputParamsStruct *params, struct gengetopt_args_info args_info)
{
   
   const CHAR *fn = __func__;
   
   //Defaults given or option passed
   params->Tcoh = args_info.Tcoh_arg;
   params->SFToverlap = args_info.SFToverlap_arg;
   params->Pmin = args_info.Pmin_arg;
   params->blksize = args_info.blksize_arg;
   params->dopplerMultiplier = args_info.dopplerMultiplier_arg;
   params->templatelength = args_info.templateLength_arg;
   params->ihsfar = args_info.ihsfar_arg;
   params->templatefar = args_info.tmplfar_arg;
   params->ULmindf = args_info.ULminimumDeltaf_arg;
   params->ULmaxdf = args_info.ULmaximumDeltaf_arg;
   params->ihsfactor = args_info.ihsfactor_arg;
   params->rootFindingMethod = args_info.BrentsMethod_given;
   params->antennaOff = args_info.antennaOff_given;
   params->noiseWeightOff = args_info.noiseWeightOff_given;
   params->calcRthreshold = args_info.calcRthreshold_given;
   params->markBadSFTs = args_info.markBadSFTs_given;
   params->FFTplanFlag = args_info.FFTplanFlag_arg;
   params->printAllULvalues = args_info.allULvalsPerSkyLoc_given;
   params->fastchisqinv = args_info.fastchisqinv_given;
   
   //Non-default arguments
   if (args_info.Tobs_given) params->Tobs = args_info.Tobs_arg;
   else params->Tobs = 10*168*3600;
   if (args_info.fmin_given) params->fmin = args_info.fmin_arg;
   else params->fmin = 99.9;
   if (args_info.fspan_given) params->fspan = args_info.fspan_arg;
   else params->fspan = 0.2;
   if (args_info.t0_given) params->searchstarttime = args_info.t0_arg;
   else params->searchstarttime = 900000000.0;
   if (args_info.Pmax_given) params->Pmax = args_info.Pmax_arg;
   else params->Pmax = 0.2*(params->Tobs);
   if (args_info.dfmin_given) params->dfmin = args_info.dfmin_arg;
   else params->dfmin = 0.5/(params->Tcoh);
   if (args_info.dfmax_given) params->dfmax = args_info.dfmax_arg;
   else params->dfmax = maxModDepth(params->Pmax, params->Tcoh);
   if (args_info.ULfmin_given) params->ULfmin = args_info.ULfmin_arg;
   else params->ULfmin = params->fmin;
   if (args_info.ULfspan_given) params->ULfspan = args_info.ULfspan_arg;
   else params->ULfspan = params->fspan;
   
   //Settings for IHS FOM
   if (args_info.ihsfomfar_given) params->ihsfomfar = args_info.ihsfomfar_arg;
   else params->ihsfomfar = 0.0;
   if (args_info.ihsfom_given) params->ihsfom = args_info.ihsfom_arg;
   else params->ihsfom = 0.0;
   if ((params->ihsfom!=0.0 && params->ihsfomfar!=0.0) || (params->ihsfom==0.0 && params->ihsfomfar==0.0)) {
      fprintf(stderr, "%s: You must choose either the IHS FOM FAR argument or the IHS FOM argument.\n", fn);
      XLAL_ERROR(fn, XLAL_EINVAL);
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
   if (2.0*(params->dfmax)+6.0/(params->Tcoh) > params->fspan) {
      params->dfmax = floor(0.5*(params->fspan)*(params->Tcoh))/(params->Tcoh) - 6.0/(params->Tcoh);
      fprintf(LOG,"WARNING! Adjusting input maximum modulation depth due to frequency span of band\n");
      fprintf(stderr,"WARNING! Adjusting input maximum modulation depth due to frequency span of band\n");
   }
   if (params->dfmin < 0.5/(params->Tcoh)) {
      params->dfmin = 0.5/(params->Tcoh);
      fprintf(LOG,"WARNING! Adjusting input minimum modulation depth to 1/2 a frequency bin!\n");
      fprintf(stderr,"WARNING! Adjusting input minimum modulation depth to 1/2 a frequency bin!\n");
   }
   
   //Print to log file and stderr the parameters of the search
   fprintf(LOG,"Tobs = %f sec\n",params->Tobs);
   fprintf(LOG,"Tcoh = %f sec\n",params->Tcoh);
   fprintf(LOG,"SFToverlap = %f sec\n",params->SFToverlap);
   fprintf(LOG,"fmin = %f Hz\n",params->fmin);
   fprintf(LOG,"fspan = %f Hz\n",params->fspan);
   fprintf(LOG,"Pmin = %f s\n",params->Pmin);
   fprintf(LOG,"Pmax = %f s\n",params->Pmax);
   fprintf(LOG,"dfmin = %f Hz\n",params->dfmin);
   fprintf(LOG,"dfmax = %f Hz\n",params->dfmax);
   fprintf(LOG,"Running median blocksize = %d\n",params->blksize);
   fprintf(LOG,"FFT plan flag = %d\n", params->FFTplanFlag);
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
   if (args_info.ihsfomfar_given) {
      fprintf(LOG,"IHS FOM FAR = %f\n", params->ihsfomfar);
      fprintf(stderr,"IHS FOM FAR = %f\n", params->ihsfomfar);
   } else {
      fprintf(LOG,"IHS FOM = %f\n", params->ihsfom);
      fprintf(stderr,"IHS FOM = %f\n", params->ihsfom);
   }

   
   //Root finding method
   if (args_info.BrentsMethod_given == 0) {
      fprintf(LOG,"Using Newton's method for root finding.\n");
      fprintf(stderr,"Using Newton's method for root finding.\n");
   } else {
      fprintf(LOG,"Using Brent's method for root finding.\n");
      fprintf(stderr,"Using Brent's method for root finding.\n");
   }
   
   //SFT type MFD or vladimir
   params->sftType = XLALCalloc(strlen(args_info.sftType_arg)+1, sizeof(*(params->sftType)));
   if (params->sftType==NULL) {
      fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", fn, sizeof(*(params->sftType)));
      XLAL_ERROR(fn, XLAL_ENOMEM);
   }
   sprintf(params->sftType, "%s", args_info.sftType_arg);
   if (strcmp(params->sftType, "MFD")==0) {
      fprintf(LOG,"sftType = %s\n", params->sftType);
      fprintf(stderr,"sftType = %s\n", params->sftType);
   } else if (strcmp(params->sftType, "vladimir")==0) {
      fprintf(LOG,"sftType = %s\n", params->sftType);
      fprintf(stderr,"sftType = %s\n", params->sftType);
   } else {
      fprintf(stderr, "%s: Not using valid type of SFT! Expected 'MFD' or 'vladimir' not %s.\n", fn, params->sftType);
      XLAL_ERROR(fn, XLAL_EINVAL);
   }
   
   //Interferometer
   CHAR *IFO = XLALCalloc(strlen(args_info.IFO_arg)+1, sizeof(*IFO));
   if (IFO==NULL) {
      fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", fn, sizeof(*IFO));
      XLAL_ERROR(fn, XLAL_ENOMEM);
   }
   sprintf(IFO, "%s", args_info.IFO_arg);
   if (strcmp("L1", IFO)==0) {
      fprintf(LOG,"IFO = %s\n", IFO);
      fprintf(stderr,"IFO = %s\n", IFO);
      params->det = lalCachedDetectors[LAL_LLO_4K_DETECTOR]; //L1
   } else if (strcmp("H1", IFO)==0) {
      fprintf(LOG,"IFO = %s\n", IFO);
      fprintf(stderr,"IFO = %s\n", IFO);
      params->det = lalCachedDetectors[LAL_LHO_4K_DETECTOR]; //H1
   } else if (strcmp("V1", IFO)==0) {
      fprintf(LOG,"IFO = %s\n", IFO);
      fprintf(stderr,"IFO = %s\n", IFO);
      params->det = lalCachedDetectors[LAL_VIRGO_DETECTOR]; //V1
   } else {
      fprintf(stderr, "%s: Not using valid interferometer! Expected 'H1', 'L1', or 'V1' not %s.\n", fn, IFO);
      XLAL_ERROR(fn, XLAL_EINVAL);
   }
   XLALFree((CHAR*)IFO);
   
   //Allocate memory for files and directory
   earth_ephemeris = XLALCalloc(strlen(args_info.ephemDir_arg)+25, sizeof(*earth_ephemeris));
   sun_ephemeris = XLALCalloc(strlen(args_info.ephemDir_arg)+25, sizeof(*sun_ephemeris));
   sft_dir = XLALCalloc(strlen(args_info.sftDir_arg)+20, sizeof(*sft_dir));
   if (earth_ephemeris==NULL) {
      fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", fn, sizeof(*earth_ephemeris));
      XLAL_ERROR(fn, XLAL_ENOMEM);
   } else if (sun_ephemeris==NULL) {
      fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", fn, sizeof(*sun_ephemeris));
      XLAL_ERROR(fn, XLAL_ENOMEM);
   } else if (sft_dir==NULL) {
      fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", fn, sizeof(*sft_dir));
      XLAL_ERROR(fn, XLAL_ENOMEM);
   }
   sprintf(earth_ephemeris, "%s/earth%s.dat", args_info.ephemDir_arg, args_info.ephemYear_arg);
   sprintf(sun_ephemeris, "%s/sun%s.dat", args_info.ephemDir_arg, args_info.ephemYear_arg);
   sprintf(sft_dir, "%s/*.sft", args_info.sftDir_arg);
   
   return 0;
   
} /* readTwoSepctInputParams() */



REAL4Vector * SSVectorMultiply_with_stride_and_offset(REAL4Vector *output, REAL4Vector *input1, REAL4Vector *input2, INT4 stride1, INT4 stride2, INT4 offset1, INT4 offset2)
{
   
   REAL4 *a, *b, *c;
   INT4   n;
   
   a = input1->data + offset1;
   b = input2->data + offset2;
   c = output->data;
   n = output->length;
   
   while (n-- > 0) {
      *c = (*a)*(*b);
      a = a + stride1;
      b = b + stride2;
      c++;
   }
   
   return output;
   
} /* SSVectorMultiply_with_stride_and_offset() */

REAL8Vector * DDVectorMultiply_with_stride_and_offset(REAL8Vector *output, REAL8Vector *input1, REAL8Vector *input2, INT4 stride1, INT4 stride2, INT4 offset1, INT4 offset2)
{
   
   REAL8 *a;
   REAL8 *b;
   REAL8 *c;
   INT4   n;
   
   a = input1->data + offset1;
   b = input2->data + offset2;
   c = output->data;
   n = output->length;
   
   while (n-- > 0) {
      *c = (*a)*(*b);
      a = a + stride1;
      b = b + stride2;
      c++;
   }
   
   return output;
   
} /* DDVectorMultiply_with_stride_and_offset() */

REAL8Vector * SDVectorMultiply_with_stride_and_offset(REAL8Vector *output, REAL4Vector *input1, REAL8Vector *input2, INT4 stride1, INT4 stride2, INT4 offset1, INT4 offset2)
{
   
   REAL4 *a;
   REAL8 *b;
   REAL8 *c;
   INT4   n;
   
   a = input1->data + offset1;
   b = input2->data + offset2;
   c = output->data;
   n = output->length;
   
   while (n-- > 0) {
      *c = (*a)*(*b);
      a = a + stride1;
      b = b + stride2;
      c++;
   }
   
   return output;
   
} /* SDVectorMultiply_with_stride_and_offset() */

