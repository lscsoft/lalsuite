/*
*  Copyright (C) 2010 Evan Goetz
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

#include <stdio.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>


#include <lal/LALConstants.h>
#include <lal/Window.h>
#include <lal/LALMalloc.h>
#include <lal/SFTutils.h>
#include <lal/SFTfileIO.h>
#include <lal/DopplerScan.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>

#include "cmdline.h"
#include "IHS.h"
#include "candidates.h"
#include "antenna.h"
#include "templates.h"
#include "TwoSpect.h"


//Global variables
inputParamsStruct *inputParams;

FILE *LOG = NULL, *TFDATA = NULL, *FFDATA = NULL;

CHAR *earth_ephemeris = NULL, *sun_ephemeris = NULL, *sft_dir = NULL;



//Main program
int main(int argc, char *argv[])
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj, kk, ll;
   REAL4 ihsfarthresh, templatefarthresh, inputExpectedSqrtSh;
   LALStatus status;
   status.statusPtr = NULL;
   char s[20000], t[20000], u[20000];
   
   //Turn off gsl error handler
   gsl_set_error_handler_off();
   
   struct gengetopt_args_info args_info;
   struct cmdline_parser_params *configparams;
   configparams = cmdline_parser_params_create();
   if ( cmdline_parser(argc, argv, &args_info) ) exit(-1);
   if ( args_info.config_given ) {
      if ( cmdline_parser_config_file(args_info.config_arg, &args_info, configparams) ) exit(-1);
   }
   
   //Set lalDebugLevel to user input or 0 if no input
   lalDebugLevel = args_info.verbosity_arg;
   
   //Create directory
   if (args_info.outdirectory_given) {
      mkdir(args_info.outdirectory_arg, 0777);
      snprintf(s, 20000, "%s/logfile.txt", args_info.outdirectory_arg);
      snprintf(t, 20000, "%s/tfdata.dat", args_info.outdirectory_arg);
      snprintf(u, 20000, "%s/ffdata.dat", args_info.outdirectory_arg);
   } else {
      mkdir("output",0777);
      snprintf(s, 20000, "%s/logfile.txt", "output");
      snprintf(t, 20000, "%s/tfdata.dat", "output");
      snprintf(u, 20000, "%s/ffdata.dat", "output");
   }
   
   //Open log file
   LOG = fopen(s,"w");
   if (LOG==NULL) {
      fprintf(stderr,"Log file could not be opened.\n");
      exit(-1);
   }
   
   //Allocate input parameters structure memory
   inputParams = new_inputParams();
   if (inputParams==NULL) {
      fprintf(stderr,"%s: new_inputParams() failed.\n", fn);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   
    //Defaults given or option passed
   inputParams->Tcoh = args_info.Tcoh_arg;
   inputParams->SFToverlap = args_info.SFToverlap_arg;
   inputParams->Pmin = args_info.Pmin_arg;
   ihsfarthresh = args_info.ihsfar_arg;
   templatefarthresh = args_info.tmplfar_arg;
   inputExpectedSqrtSh =  args_info.avesqrtSh_arg;
   inputParams->blksize = args_info.blksize_arg;
   inputParams->dopplerMultiplier = args_info.dopplerMultiplier_arg;
   inputParams->templatelength = args_info.templateLength_arg;
   earth_ephemeris = XLALCalloc(strlen(args_info.ephemDir_arg)+20, sizeof(*earth_ephemeris));
   sun_ephemeris = XLALCalloc(strlen(args_info.ephemDir_arg)+20, sizeof(*sun_ephemeris));
   sft_dir = XLALCalloc(strlen(args_info.sftDir_arg)+20, sizeof(*sft_dir));
   if (earth_ephemeris==NULL) {
      fprintf(stderr,"%s: XLALCalloc(%lu) failed.\n", fn, sizeof(*earth_ephemeris));
      XLAL_ERROR(fn, XLAL_ENOMEM);
   } else if (sun_ephemeris==NULL) {
      fprintf(stderr,"%s: XLALCalloc(%lu) failed.\n", fn, sizeof(*sun_ephemeris));
      XLAL_ERROR(fn, XLAL_ENOMEM);
   } else if (sft_dir==NULL) {
      fprintf(stderr,"%s: XLALCalloc(%lu) failed.\n", fn, sizeof(*sft_dir));
      XLAL_ERROR(fn, XLAL_ENOMEM);
   }
   
   sprintf(earth_ephemeris,"%s/earth05-09.dat",args_info.ephemDir_arg);
   sprintf(sun_ephemeris,"%s/sun05-09.dat",args_info.ephemDir_arg);
   sprintf(sft_dir,"%s/*.sft",args_info.sftDir_arg);
   
   //Non-default arguments
   if (args_info.Tobs_given) inputParams->Tobs = args_info.Tobs_arg;
   else inputParams->Tobs = 10*168*3600;
   if (args_info.fmin_given) inputParams->fmin = args_info.fmin_arg;
   else inputParams->fmin = 99.9;
   if (args_info.fspan_given) inputParams->fspan = args_info.fspan_arg;
   else inputParams->fspan = 0.2;
   if (args_info.t0_given) inputParams->searchstarttime = args_info.t0_arg;
   else inputParams->searchstarttime = 900000000.0;
   if (args_info.Pmax_given) inputParams->Pmax = args_info.Pmax_arg;
   else inputParams->Pmax = 0.2*inputParams->Tobs;
   if (args_info.dfmin_given) inputParams->dfmin = args_info.dfmin_arg;
   else inputParams->dfmin = 0.5/inputParams->Tcoh;
   if (args_info.dfmax_given) inputParams->dfmax = args_info.dfmax_arg;
   else inputParams->dfmax = maxModDepth(inputParams->Pmax,inputParams->Tcoh);
   
   //Blocksize should be an odd number
   if (inputParams->blksize % 2 != 1) inputParams->blksize += 1;
   
   // Warnings when using hidden flags
   if (args_info.antennaOff_given) {
      fprintf(LOG,"WARNING: Antenna pattern weights are all being set to 1.0\n");
      fprintf(stderr,"WARNING: Antenna pattern weights are all being set to 1.0\n");
   }
   if (args_info.gaussTemplatesOnly_given) {
      fprintf(LOG,"WARNING: Only Gaussian templates will be used\n");
      fprintf(stderr,"WARNING: Only Gaussian templates will be used\n");
   }
   if (args_info.IHSonly_given) {
      fprintf(LOG,"WARNING: Only IHS stage is being used\n");
      fprintf(stderr,"WARNING: Only IHS stage is being used\n");
   }
   
   //Adjust parameter space search values, if necessary
   if (inputParams->Pmax < inputParams->Pmin) {
      REAL4 tempP = inputParams->Pmax;
      inputParams->Pmax = inputParams->Pmin;
      inputParams->Pmin = tempP;
      fprintf(LOG,"WARNING! Maximum period is smaller than minimum period... switching the two\n");
      fprintf(stderr,"WARNING! Maximum period is smaller than minimum period... switching the two\n");
   }
   if (inputParams->dfmax < inputParams->dfmin) {
      REAL4 tempdf = inputParams->dfmax;
      inputParams->dfmax = inputParams->dfmin;
      inputParams->dfmin = tempdf;
      fprintf(LOG,"WARNING! Maximum modulation depth is smaller than minimum modulation depth... switching the two\n");
      fprintf(stderr,"WARNING! Maximum modulation depth is smaller than minimum modulation depth... switching the two\n");
   }
   if (inputParams->Pmax > 0.2*inputParams->Tobs) {
      inputParams->Pmax = 0.2*inputParams->Tobs;
      fprintf(LOG,"WARNING! Adjusting input maximum period to 1/5 the observation time!\n");
      fprintf(stderr,"WARNING! Adjusting input maximum period to 1/5 the observation time!\n");
   }
   if (inputParams->Pmin < 2.0*3600) {
      inputParams->Pmin = 2.0*3600;
      fprintf(LOG,"WARNING! Adjusting input minimum period to 2 hours!\n");
      fprintf(stderr,"WARNING! Adjusting input minimum period to 2 hours!\n");
   }
   if (inputParams->dfmax > maxModDepth(inputParams->Pmax,inputParams->Tcoh)) {
      inputParams->dfmax = 0.5*floor(2.0*maxModDepth(inputParams->Pmax,inputParams->Tcoh)*inputParams->Tcoh)/inputParams->Tcoh;
      fprintf(LOG,"WARNING! Adjusting input maximum modulation depth due to maximum modulation depth allowed\n");
      fprintf(stderr,"WARNING! Adjusting input maximum modulation depth due to maximum modulation depth allowed\n");
   }
   if (2.0*inputParams->dfmax+6.0/inputParams->Tcoh > inputParams->fspan) {
      inputParams->dfmax = floor(0.5*(inputParams->fspan - 6.0/inputParams->Tcoh)*inputParams->Tcoh)/inputParams->Tcoh;
      fprintf(LOG,"WARNING! Adjusting input maximum modulation depth due to frequency span of band\n");
      fprintf(stderr,"WARNING! Adjusting input maximum modulation depth due to frequency span of band\n");
   }
   if (inputParams->dfmin < 0.5/inputParams->Tcoh) {
      inputParams->dfmin = 0.5/inputParams->Tcoh;
      fprintf(LOG,"WARNING! Adjusting input minimum modulation depth to 1/2 a frequency bin!\n");
      fprintf(stderr,"WARNING! Adjusting input minimum modulation depth to 1/2 a frequency bin!\n");
   }
   
   //Parameters for the sky-grid
   CHAR *sky = XLALCalloc(strlen(args_info.skyRegion_arg)+1, sizeof(*sky));
   if (sky==NULL) {
      fprintf(stderr,"%s: XLALCalloc(%lu) failed.\n", fn, sizeof(*sky));
      XLAL_ERROR(fn, XLAL_ENOMEM);
   }
   CHAR *IFO = XLALCalloc(strlen(args_info.IFO_arg)+1, sizeof(*IFO));
   if (IFO==NULL) {
      fprintf(stderr,"%s: XLALCalloc(%lu) failed.\n", fn, sizeof(*IFO));
      XLAL_ERROR(fn, XLAL_ENOMEM);
   }
   sprintf(sky,"%s",args_info.skyRegion_arg);
   sprintf(IFO,"%s",args_info.IFO_arg);
   DopplerSkyScanInit scanInit = empty_DopplerSkyScanInit;
   DopplerSkyScanState scan = empty_DopplerSkyScanState;
   PulsarDopplerParams dopplerpos;
   scanInit.gridType = 1;     //Default value for an approximate-isotropic grid
   scanInit.skyRegionString = sky;      //"allsky" = Default value for all-sky search
   scanInit.numSkyPartitions = 1;   //Default value so sky is not broken into chunks
   scanInit.Freq = args_info.fmin_arg;  
   
   //Detector velocity during SFTs
   LALDetector det;
   if (strcmp("L1", IFO)==0) {
      det = lalCachedDetectors[LALDetectorIndexLLODIFF]; //L1
   } else {
      det = lalCachedDetectors[LALDetectorIndexLHODIFF]; //H1
   }
   inputParams->det = &det;
   //EphemerisData *edat = new_Ephemeris(earth_ephemeris, sun_ephemeris);
   EphemerisData *edat = XLALInitBarycenter(earth_ephemeris, sun_ephemeris);
   if (edat==NULL) {
      fprintf(stderr,"%s: XLALInitBarycenter() failed.\n", fn);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   REAL4 detectorVmax = 9.93e-5;   //Average orbital earth speed in units of c
   
   //Initialize the sky-grid
   scanInit.dAlpha = 0.5/(inputParams->fmin * inputParams->Tcoh * detectorVmax);
   scanInit.dDelta = scanInit.dAlpha;
   InitDopplerSkyScan(&status, &scan, &scanInit);
   if (status.statusCode!=0) {
      fprintf(stderr,"%s: InitDopplerSkyScan() failed.\n", fn);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   
   //Start at first location
   if ((XLALNextDopplerSkyPos(&dopplerpos, &scan))!=0) {
      fprintf(stderr,"%s: XLALNextDopplerSkyPos() failed.\n", fn);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   
   //Print to log file and stderr the parameters of the search
   fprintf(LOG,"Tobs = %f sec\n",inputParams->Tobs);
   fprintf(LOG,"Tcoh = %f sec\n",inputParams->Tcoh);
   fprintf(LOG,"SFToverlap = %f sec\n",inputParams->SFToverlap);
   fprintf(LOG,"fmin = %f Hz\n",inputParams->fmin);
   fprintf(LOG,"fspan = %f Hz\n",inputParams->fspan);
   fprintf(LOG,"Pmin = %f s\n",inputParams->Pmin);
   fprintf(LOG,"Pmax = %f s\n",inputParams->Pmax);
   fprintf(LOG,"dfmin = %f Hz\n",inputParams->dfmin);
   fprintf(LOG,"dfmax = %f Hz\n",inputParams->dfmax);
   fprintf(LOG,"Sky region = %s\n",sky);
   fprintf(LOG,"Running median blocksize = %d\n",inputParams->blksize);
   fprintf(LOG,"FAR for IHS = %f, for templates = %f\n",ihsfarthresh,templatefarthresh);
   fprintf(stderr,"Tobs = %f sec\n",inputParams->Tobs);
   fprintf(stderr,"Tcoh = %f sec\n",inputParams->Tcoh);
   fprintf(stderr,"SFToverlap = %f sec\n",inputParams->SFToverlap);
   fprintf(stderr,"fmin = %f Hz\n",inputParams->fmin);
   fprintf(stderr,"fspan = %f Hz\n",inputParams->fspan);
   fprintf(stderr,"Pmin = %f s\n",inputParams->Pmin);
   fprintf(stderr,"Pmax = %f s\n",inputParams->Pmax);
   fprintf(stderr,"dfmin = %f Hz\n",inputParams->dfmin);
   fprintf(stderr,"dfmax = %f Hz\n",inputParams->dfmax);
   fprintf(stderr,"Sky region = %s\n",sky);
   fprintf(stderr,"Running median blocksize = %d\n",inputParams->blksize);
   fprintf(stderr,"FAR for IHS = %f, for templates = %f\n",ihsfarthresh,templatefarthresh);
   
   //Basic units
   INT4 numffts = (INT4)floor(inputParams->Tobs/(inputParams->Tcoh-inputParams->SFToverlap)-1);
   INT4 numfbins = (INT4)round(inputParams->fspan*inputParams->Tcoh)+1;
   INT4 numfprbins = (INT4)floor(numffts*0.5) + 1;
   REAL4 tempfspan = inputParams->fspan + (inputParams->blksize-1)/inputParams->Tcoh;
   INT4 tempnumfbins = (INT4)round(tempfspan*inputParams->Tcoh)+1;
   
   //Allocate memory for ffdata structure
   ffdataStruct *ffdata = new_ffdata(inputParams);
   if (ffdata==NULL) {
      fprintf(stderr,"%s: new_ffdata() failed.\n", fn);
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
   if (gaussCandidates1==NULL) {
      fprintf(stderr,"%s: new_CandidateVector(%d) failed.\n", fn, 100);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (gaussCandidates2==NULL) {
      fprintf(stderr,"%s: new_CandidateVector(%d) failed.\n", fn, 100);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (gaussCandidates3==NULL) {
      fprintf(stderr,"%s: new_CandidateVector(%d) failed.\n", fn, 100);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (gaussCandidates4==NULL) {
      fprintf(stderr,"%s: new_CandidateVector(%d) failed.\n", fn, 100);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (exactCandidates1==NULL) {
      fprintf(stderr,"%s: new_CandidateVector(%d) failed.\n", fn, 100);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (exactCandidates2==NULL) {
      fprintf(stderr,"%s: new_CandidateVector(%d) failed.\n", fn, 100);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (ihsCandidates==NULL) {
      fprintf(stderr,"%s: new_CandidateVector(%d) failed.\n", fn, 100);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   
   //Second fft plan, only need to make this once for all the exact templates
   REAL4FFTPlan *secondFFTplan = XLALCreateForwardREAL4FFTPlan((UINT4)numffts, 0);
   if (secondFFTplan==NULL) {
      fprintf(stderr,"%s: XLALCreateForwardREAL4FFTPlan(%d,%d) failed.\n", fn, numffts, 0);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   
   //Maximum number of IHS values to sum = twice the maximum modulation depth
   //Minimum number of IHS values to sum = twice the minimum modulation depth
   INT4 maxcols = (INT4)floor(2.0*inputParams->dfmax*inputParams->Tcoh)+1;
   
   //Assume maximum bin shift possible
   inputParams->maxbinshift = (INT4)round(detectorVmax * inputParams->fmin * inputParams->Tcoh)+1; //TODO: better way to do this?
   
   //Read in the T-F data from SFTs
   fprintf(LOG,"Loading in SFTs... ");
   fprintf(stderr,"Loading in SFTs... ");
   ffdata->tfnormalization = 2.0/inputParams->Tcoh/(inputExpectedSqrtSh*inputExpectedSqrtSh);
   REAL4Vector *tfdata = readInSFTs(inputParams, &(ffdata->tfnormalization));
   if (tfdata==NULL) {
      fprintf(stderr,"%s: readInSFTs() failed.\n", fn);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   fprintf(LOG,"done\n");
   fprintf(stderr,"done\n");
   
   //Calculate the running mean values of the SFTs (output here is smaller than initialTFdata). Here,
   //numfbins needs to be the bins you expect to come out of the running means -- the band you are going
   //to search!
   fprintf(LOG,"Assessing background... ");
   fprintf(stderr,"Assessing background... ");
   REAL4Vector *background = XLALCreateREAL4Vector((UINT4)(numffts*(numfbins + 2*inputParams->maxbinshift)));
   if (background==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (UINT4)(numffts*(numfbins + 2*inputParams->maxbinshift)));
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   tfRngMeans(background, tfdata, numffts, numfbins + 2*inputParams->maxbinshift, inputParams->blksize);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: tfRngMeans() failed.\n", fn);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   
   //I wrote this to compensate for a bad input of the expected noise floor and for non-present SFTs
   REAL8 backgroundmeannormfactor = 0.0;
   INT4 avefact = 0;
   for (ii=0; ii<(INT4)background->length; ii++) {
      if (background->data[ii]!=0.0) {
         backgroundmeannormfactor += background->data[ii];
         avefact++;
      }
   }
   backgroundmeannormfactor = avefact/backgroundmeannormfactor;
   for (ii=0; ii<(INT4)background->length; ii++) background->data[ii] *= backgroundmeannormfactor;
   ffdata->tfnormalization *= backgroundmeannormfactor;
   fprintf(LOG,"done\n");
   fprintf(stderr,"done\n");
   
   //Need to reduce the original TF data to remove the excess bins used for running median calculation. Normalize the TF as the same as the background was normalized
   REAL4Vector *usableTFdata = XLALCreateREAL4Vector(background->length);
   if (usableTFdata==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, background->length);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   for (ii=0; ii<numffts; ii++) {
      for (jj=0; jj<numfbins+2*inputParams->maxbinshift; jj++) usableTFdata->data[ii*(numfbins+2*inputParams->maxbinshift) + jj] = backgroundmeannormfactor * tfdata->data[ii*(tempnumfbins+2*inputParams->maxbinshift) + jj + (INT4)(0.5*(inputParams->blksize-1))];
   }
   //At this point the TF plane and the running median calculation are the same size=numffts*(numfbins + 2*maxbinshift)
   //We can delete the originally loaded SFTs since we have the usableTFdata saved
   XLALDestroyREAL4Vector(tfdata);
   
   //IHS FOM (allows a relative offset of +/- 1 bin between maximum values
   REAL4 ihsfomfar = 6.0;
   fprintf(LOG,"IHS FOM FAR = %f\n",ihsfomfar);
   fprintf(stderr,"IHS FOM FAR = %f\n",ihsfomfar);
   
   fprintf(LOG,"Maximum column width to be searched = %d\n",maxcols);
   fprintf(stderr,"Maximum column width to be searched = %d\n",maxcols);
   
   //Set number of candidates to be zero
   //numofcandidates = numofcandidates2 = numofcandidatesadded = 0;
   
   //Initialize reused values
   ihsMaximaStruct *ihsmaxima = new_ihsMaxima(numfbins, maxcols);
   ihsfarStruct *ihsfarstruct = new_ihsfarStruct(maxcols);
   REAL4Vector *detectorVelocities = XLALCreateREAL4Vector((UINT4)numffts);
   INT4Vector *binshifts = XLALCreateINT4Vector((UINT4)numffts);
   REAL4Vector *aveNoise = XLALCreateREAL4Vector((UINT4)numfprbins);
   if (ihsmaxima==NULL) {
      fprintf(stderr,"%s: new_ihsMaxima(%d,%d) failed.\n", fn, numfbins, maxcols);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (ihsfarstruct==NULL) {
      fprintf(stderr,"%s: new_ihsfarStruct(%d) failed.\n", fn, maxcols);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (detectorVelocities==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numffts);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (binshifts==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, numffts);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } else if (aveNoise==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numfprbins);
      XLAL_ERROR(fn, XLAL_EFUNC);
   } 
   
   //Initialize to zero for far just at the start
   ihsfarstruct->ihsfar->data[0] = 0.0;
   REAL4 antweightsrms = 0.0;
   
   //R probability calculator errorcode
   INT4 proberrcode = 0;
   
   //Print message that we start the analysis
   fprintf(LOG,"Starting TwoSpect analysis...\n");
   fprintf(stderr,"Starting TwoSpect analysis...\n");
   
   //Search over the sky region
   while (scan.state != STATE_FINISHED) {
      fprintf(LOG,"Sky location: RA = %g, DEC = %g\n", dopplerpos.Alpha, dopplerpos.Delta);
      fprintf(stderr,"Sky location: RA = %g, DEC = %g\n", dopplerpos.Alpha, dopplerpos.Delta);
      
      //Determine detector velocity w.r.t. a sky location for each SFT
      CompAntennaVelocity(detectorVelocities, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, inputParams->searchstarttime, inputParams->Tcoh, inputParams->SFToverlap, inputParams->Tobs, det, edat);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: CompAntennaVelocity() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      
      //Compute the bin shifts for each SFT
      CompBinShifts(binshifts, inputParams->fmin+inputParams->fspan*0.5, detectorVelocities, inputParams->Tcoh, inputParams->dopplerMultiplier);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: CompBinShifts() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      
      //Compute antenna pattern weights. If antennaOff input flag is given, then set all values equal to 1.0
      REAL4Vector *antweights = XLALCreateREAL4Vector((UINT4)numffts);
      if (antweights==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numffts);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      CompAntennaPatternWeights(antweights, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, inputParams->searchstarttime, inputParams->Tcoh, inputParams->SFToverlap, inputParams->Tobs, det);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: CompAntennaPatternWeights() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      //If antenna weights are set to be "off" then set all weights equal to 1.0
      if (args_info.antennaOff_given) {
         for (ii=0; ii<(INT4)antweights->length; ii++) antweights->data[ii] = 1.0;
      }
      
      //Calculate antenna RMS value
      REAL4 currentAntWeightsRMS = calcRms(antweights);
      if (XLAL_IS_REAL4_FAIL_NAN(currentAntWeightsRMS)) {
         fprintf(stderr,"%s: calcRms() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      
      //Slide SFTs here -- need to slide the data and the estimated background
      REAL4Vector *TFdata_slided = XLALCreateREAL4Vector((UINT4)(numffts*numfbins));
      REAL4Vector *background_slided = XLALCreateREAL4Vector(TFdata_slided->length);
      if (TFdata_slided==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (UINT4)(numffts*numfbins));
         XLAL_ERROR(fn, XLAL_EFUNC);
      } else if (background_slided==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, TFdata_slided->length);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      slideTFdata(TFdata_slided, inputParams, usableTFdata, binshifts);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: slideTFdata() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      slideTFdata(background_slided, inputParams, background, binshifts);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: slideTFdata() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      fprintf(stderr,"Mean of TFdata_slided %g, background_slided %g\n",calcMean(TFdata_slided),calcMean(background_slided));
      
      //Check the RMS of the antenna weights, if bigger than standard deviation then reset the IHS FAR and the average noise background of the 2nd FFT
      if (antweightsrms == 0.0) {
         antweightsrms = currentAntWeightsRMS;
      }
      if ( fabs(currentAntWeightsRMS-antweightsrms)/antweightsrms >= 0.03 ) {
         ihsfarstruct->ihsfar->data[0] = 0.0;
         antweightsrms = currentAntWeightsRMS;
      }
      
      //Average noise floor of FF plane for each 1st FFT frequency bin
      ffdata->ffnormalization = 1.0;
      ffPlaneNoise(aveNoise, inputParams, background_slided, antweights, &(ffdata->ffnormalization));
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: ffPlaneNoise() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      
      //Calculation of average TF noise per frequency bin ratio to total mean
      REAL4 aveTFinv = 1.0/avgTFdataBand(background_slided, numfbins, numffts, 0, numfbins);
      REAL4 rmsTFinv = 1.0/avgTFdataBand(background_slided, numfbins, numffts, 0, numfbins);
      REAL4Vector *aveTFnoisePerFbinRatio = XLALCreateREAL4Vector((UINT4)numfbins);
      REAL4Vector *rmsTFnoisePerFbinRatio = XLALCreateREAL4Vector((UINT4)numfbins);
      REAL4Vector *TSofPowers = XLALCreateREAL4Vector((UINT4)numffts);
      if (aveTFnoisePerFbinRatio==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numfbins);
         XLAL_ERROR(fn, XLAL_EFUNC);
      } else if (rmsTFnoisePerFbinRatio==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numfbins);
         XLAL_ERROR(fn, XLAL_EFUNC);
      } else if (TSofPowers==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numffts);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      for (ii=0; ii<numfbins; ii++) {
         for (jj=0; jj<numffts; jj++) TSofPowers->data[jj] = background_slided->data[jj*numfbins + ii];
         aveTFnoisePerFbinRatio->data[ii] = calcMean(TSofPowers)*aveTFinv;
         rmsTFnoisePerFbinRatio->data[ii] = calcRms(TSofPowers)*rmsTFinv;
      }
      XLALDestroyREAL4Vector(TSofPowers);
      
      //Compute the weighted TF data
      REAL4Vector *TFdata_weighted = XLALCreateREAL4Vector((UINT4)(numffts*numfbins));
      if (TFdata_weighted==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (UINT4)(numffts*numfbins));
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      tfWeightMeanSubtract(TFdata_weighted, TFdata_slided, background_slided, antweights, inputParams);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: tfWeightMeanSubtract() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      XLALDestroyREAL4Vector(TFdata_slided);
      XLALDestroyREAL4Vector(background_slided);
      XLALDestroyREAL4Vector(antweights);
      /* TFDATA = fopen(t,"w");
      for (jj=0; jj<(INT4)TFdata_weighted->length; jj++) fprintf(TFDATA,"%g\n",TFdata_weighted->data[jj]);
      fclose(TFDATA); */
      
      //Do the second FFT
      makeSecondFFT(ffdata->ffdata, ffdata->ffnormalization, TFdata_weighted, inputParams, secondFFTplan);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: makeSecondFFT() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      REAL4 secFFTmean = calcMean(ffdata->ffdata);
      
      /* REAL8 secfftnormfactor = 1.0/calcMean(ffdata->ffdata);
      ffdata->ffnormalization *= secfftnormfactor;
      for (ii=0; ii<(INT4)ffdata->ffdata->length; ii++) ffdata->ffdata->data[ii] *= secfftnormfactor; */
      
      XLALDestroyREAL4Vector(TFdata_weighted);
      fprintf(stderr,"2nd FFT ave = %g, expected ave = %g\n",secFFTmean,calcMean(aveNoise)*calcMean(aveTFnoisePerFbinRatio));
      /* FFDATA = fopen(u,"w");
      for (jj=0; jj<(INT4)ffdata->ffdata->length; jj++) fprintf(FFDATA,"%g\n",ffdata->ffdata->data[jj]);
      fclose(FFDATA); */
      
      if (secFFTmean==0.0) {
         fprintf(stderr,"Apparently, no SFTs were read in (Average power value is zero). Program exiting with failure.\n");
         XLAL_ERROR(fn, XLAL_FAILURE);
      }
      
      
////////Start of the IHS step!
      //Find the FAR of IHS sum
      if (ihsfarstruct->ihsfar->data[0]==0.0) {
         fprintf(stderr,"Determining IHS FAR values...\n");
         fprintf(LOG,"Determining IHS FAR values...\n");
         genIhsFar(ihsfarstruct, maxcols, ihsfarthresh, aveNoise, inputParams->Tobs);
         if (xlalErrno!=0) {
            fprintf(stderr,"%s: genIhsFar() failed.\n", fn);
            XLAL_ERROR(fn, XLAL_EFUNC);
         }
      }
      
      //Run the IHS algorithm on the data
      runIHS(ihsmaxima, ffdata, inputParams, maxcols);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: runIHS() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      
      //Find any IHS candidates
      findIHScandidates(ihsCandidates, ihsfarstruct, inputParams, ffdata, ihsmaxima, aveTFnoisePerFbinRatio, rmsTFnoisePerFbinRatio);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: findIHScandidates() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      fprintf(LOG,"Candidates found in IHS step = %d\n",ihsCandidates->numofcandidates);
      fprintf(stderr,"Candidates found in IHS step = %d\n",ihsCandidates->numofcandidates);
////////End of the IHS step
      
////////Start of the Gaussian template search!
      if (!args_info.IHSonly_given) {
         farStruct *farval = NULL;
         for (ii=0; ii<(INT4)ihsCandidates->numofcandidates; ii++) {
            
            //Now assess the IHS candidate if the signal is away from the band edges, the modulation depth is greater or equal to minimum specified and less than or equal to the maximum specified, and if the period/modulation depth combo is within allowable limits for a template to be made. We will cut the period space in the next step.
            if (ihsCandidates->data[ii].fsig-ihsCandidates->data[ii].moddepth-6.0/inputParams->Tcoh > inputParams->fmin && ihsCandidates->data[ii].fsig+ihsCandidates->data[ii].moddepth+6.0/inputParams->Tcoh < inputParams->fmin+inputParams->fspan && ihsCandidates->data[ii].moddepth >= inputParams->dfmin && ihsCandidates->data[ii].moddepth <= inputParams->dfmax && ihsCandidates->data[ii].moddepth < maxModDepth(ihsCandidates->data[ii].period,inputParams->Tcoh) && ihsCandidates->data[ii].period >= 2.0*3600.0) {
               
               //Allocate memory for template
               templateStruct *template = new_templateStruct(inputParams->templatelength);
               if (template==NULL) {
                  fprintf(stderr,"%s: new_templateStruct(%d) failed.\n", fn, inputParams->templatelength);
                  XLAL_ERROR(fn, XLAL_EFUNC); 
               }
               
               //Make a Gaussian train template
               makeTemplateGaussians(template, ihsCandidates->data[ii], inputParams);
               if (xlalErrno!=0) {
                  fprintf(stderr,"%s: makeTemplateGaussians() failed.\n", fn);
                  XLAL_ERROR(fn, XLAL_EFUNC);
               }
               
               //Estimate the FAR for these bin weights
               farval = new_farStruct();
               if (farval==NULL) {
                  fprintf(stderr,"%s: new_farStruct() failed.\n", fn);
                  XLAL_ERROR(fn, XLAL_EFUNC); 
               }
               numericFAR(farval, template, templatefarthresh, aveNoise, aveTFnoisePerFbinRatio);
               if (xlalErrno!=0) {
                  fprintf(stderr,"%s: numericFAR() failed.\n", fn);
                  XLAL_ERROR(fn, XLAL_EFUNC);
               }
               
               //Caclulate R, probability noise caused the candidate, and estimate of h0
               REAL8 R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
               REAL8 prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, &proberrcode);
               if (XLAL_IS_REAL8_FAIL_NAN(R)) {
                  fprintf(stderr,"%s: calculateR() failed.\n", fn);
                  XLAL_ERROR(fn, XLAL_EFUNC);
               } else if (XLAL_IS_REAL8_FAIL_NAN(prob)) {
                  fprintf(stderr,"%s: probR() failed.\n", fn);
                  XLAL_ERROR(fn, XLAL_EFUNC);
               }
               REAL8 h0 = 2.9569*pow(R/(inputParams->Tcoh*inputParams->Tobs),0.25);
               
               /* if (prob!=0.0) {
                  FILE *TEMPLATEOUT = fopen("./templatevalues.dat","w");
                  INT4 mm;
                  for (mm=0; mm<(INT4)template->templatedata->length; mm++) fprintf(TEMPLATEOUT,"%g %g %g %d %d %d\n",template->templatedata->data[mm],ffdata->ffdata->data[template->pixellocations->data[mm]],aveNoise->data[template->secondfftfrequencies->data[mm]]*aveTFnoisePerFbinRatio->data[template->firstfftfrequenciesofpixels->data[mm]],template->pixellocations->data[mm],template->firstfftfrequenciesofpixels->data[mm],template->secondfftfrequencies->data[mm]);
                  fclose(TEMPLATEOUT);
               } */
               
               //Destroy unneeded things
               free_templateStruct(template);
               template = NULL;
               
               /* Log the candidate if R exceeds the FAR or check other possibilities of different 
               periods. Use same farval.far because the FAR is ~independent on period. */
               REAL8 bestPeriod = 0.0;
               REAL8 besth0 = 0.0;
               REAL8 bestR = 0.0;
               REAL8 bestProb = 1.0;
               INT4 bestproberrcode = 0;
               REAL8 initialFAR = farval->far;
               REAL8 Rfirst = R;
               if (R > farval->far) {
                  bestR = R;
                  besth0 = h0;
                  bestProb = prob;
                  bestPeriod = ihsCandidates->data[ii].period;
                  if (gaussCandidates1->numofcandidates == gaussCandidates1->length-1) gaussCandidates1 = resize_candidateVector(gaussCandidates1, 2*gaussCandidates1->length);
                  loadCandidateData(&gaussCandidates1->data[gaussCandidates1->numofcandidates], ihsCandidates->data[ii].fsig, ihsCandidates->data[ii].period, ihsCandidates->data[ii].moddepth, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, R, besth0, bestProb, proberrcode, ihsCandidates->data[ii].normalization);
                  gaussCandidates1->numofcandidates++;
               }
               
               //Try shifting period by fractions, if necessary
               if (bestProb == 1.0) {
                  for (jj=0; jj<3; jj++) {
                     REAL8 periodfact = (jj+1.0)/(jj+2.0);
                     if ( periodfact*ihsCandidates->data[ii].period > minPeriod(ihsCandidates->data[ii].moddepth, inputParams->Tcoh) && periodfact*ihsCandidates->data[ii].period>2.0*3600.0) {
                        ihsCandidates->data[ii].period *= periodfact;
                        template = new_templateStruct(inputParams->templatelength);
                        if (template==NULL) {
                           fprintf(stderr,"%s: new_templateStruct(%d) failed.\n", fn, inputParams->templatelength);
                           XLAL_ERROR(fn, XLAL_EFUNC); 
                        }
                        makeTemplateGaussians(template, ihsCandidates->data[ii], inputParams);
                        if (xlalErrno!=0) {
                           fprintf(stderr,"%s: makeTemplateGaussians() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        }
                        R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
                        prob = (probR(template, aveNoise, aveTFnoisePerFbinRatio, R, &proberrcode));
                        if (XLAL_IS_REAL8_FAIL_NAN(R)) {
                           fprintf(stderr,"%s: calculateR() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        } else if (XLAL_IS_REAL8_FAIL_NAN(prob)) {
                           fprintf(stderr,"%s: probR() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        }
                        h0 = 2.9569*pow(R/(inputParams->Tcoh*inputParams->Tobs),0.25);
                        if (R>farval->far && prob < bestProb) {
                           bestPeriod = ihsCandidates->data[ii].period;
                           besth0 = h0;
                           bestR = R;
                           bestProb = prob;
                           bestproberrcode = proberrcode;
                        }
                        free_templateStruct(template);
                        template = NULL;
                        ihsCandidates->data[ii].period /= periodfact;
                     } /* shift shorter period */
                     periodfact = 1.0/periodfact;
                     if ( periodfact*ihsCandidates->data[ii].period < 0.2*inputParams->Tobs ) {
                        ihsCandidates->data[ii].period *= periodfact;
                        template = new_templateStruct(inputParams->templatelength);
                        if (template==NULL) {
                           fprintf(stderr,"%s: new_templateStruct(%d) failed.\n", fn, inputParams->templatelength);
                           XLAL_ERROR(fn, XLAL_EFUNC); 
                        }
                        makeTemplateGaussians(template, ihsCandidates->data[ii], inputParams);
                        if (xlalErrno!=0) {
                           fprintf(stderr,"%s: makeTemplateGaussians() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        }
                        R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
                        prob = (probR(template, aveNoise, aveTFnoisePerFbinRatio, R, &proberrcode));
                        if (XLAL_IS_REAL8_FAIL_NAN(R)) {
                           fprintf(stderr,"%s: calculateR() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        } else if (XLAL_IS_REAL8_FAIL_NAN(prob)) {
                           fprintf(stderr,"%s: probR() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        }
                        h0 = 2.9569*pow(R/(inputParams->Tcoh*inputParams->Tobs),0.25);
                        if (R>farval->far && prob < bestProb) {
                           bestPeriod = ihsCandidates->data[ii].period;
                           besth0 = h0;
                           bestR = R;
                           bestProb = prob;
                           bestproberrcode = proberrcode;
                        }
                        free_templateStruct(template);
                        template = NULL;
                        ihsCandidates->data[ii].period /= periodfact;
                     } /* shift longer period */
                  } /* for jj < 3 (fractions of period) */
               } /* if a best probability was not found */
               
               //If a potentially interesting candidate has been found (exceeding the FAR threshold) then try some new periods
               if (bestR != 0.0) {
                  //Shift period by harmonics
                  for (jj=2; jj<6; jj++) {
                     if (ihsCandidates->data[ii].period/jj > minPeriod(ihsCandidates->data[ii].moddepth, inputParams->Tcoh) && ihsCandidates->data[ii].period/jj >= 2.0*3600.0) {
                        ihsCandidates->data[ii].period /= (REAL8)jj;
                        template = new_templateStruct(inputParams->templatelength);
                        if (template==NULL) {
                           fprintf(stderr,"%s: new_templateStruct(%d) failed.\n", fn, inputParams->templatelength);
                           XLAL_ERROR(fn, XLAL_EFUNC); 
                        }
                        makeTemplateGaussians(template, ihsCandidates->data[ii], inputParams);
                        if (xlalErrno!=0) {
                           fprintf(stderr,"%s: makeTemplateGaussians() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        }
                        numericFAR(farval, template, templatefarthresh, aveNoise, aveTFnoisePerFbinRatio);
                        R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
                        prob = (probR(template, aveNoise, aveTFnoisePerFbinRatio, R, &proberrcode));
                        if (XLAL_IS_REAL8_FAIL_NAN(R)) {
                           fprintf(stderr,"%s: calculateR() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        } else if (XLAL_IS_REAL8_FAIL_NAN(prob)) {
                           fprintf(stderr,"%s: probR() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        }
                        h0 = 2.9569*pow(R/(inputParams->Tcoh*inputParams->Tobs),0.25);
                        if (R>farval->far && prob < bestProb) {
                           bestPeriod = ihsCandidates->data[ii].period;
                           besth0 = h0;
                           bestR = R;
                           bestProb = prob;
                           bestproberrcode = proberrcode;
                        }
                        free_templateStruct(template);
                        template = NULL;
                        ihsCandidates->data[ii].period *= (REAL8)jj;
                     } /* shorter period harmonics */
                     if (ihsCandidates->data[ii].period*jj < 0.2*inputParams->Tobs) {
                        ihsCandidates->data[ii].period *= (REAL8)jj;
                        template = new_templateStruct(inputParams->templatelength);
                        if (template==NULL) {
                           fprintf(stderr,"%s: new_templateStruct(%d) failed.\n", fn, inputParams->templatelength);
                           XLAL_ERROR(fn, XLAL_EFUNC); 
                        }
                        makeTemplateGaussians(template, ihsCandidates->data[ii], inputParams);
                        if (xlalErrno!=0) {
                           fprintf(stderr,"%s: makeTemplateGaussians() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        }
                        numericFAR(farval, template, templatefarthresh, aveNoise, aveTFnoisePerFbinRatio);
                        R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
                        prob = (probR(template, aveNoise, aveTFnoisePerFbinRatio, R, &proberrcode));
                        if (XLAL_IS_REAL8_FAIL_NAN(R)) {
                           fprintf(stderr,"%s: calculateR() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        } else if (XLAL_IS_REAL8_FAIL_NAN(prob)) {
                           fprintf(stderr,"%s: probR() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        }
                        h0 = 2.9569*pow(R/(inputParams->Tcoh*inputParams->Tobs),0.25);
                        if (R>farval->far && prob < bestProb) {
                           bestPeriod = ihsCandidates->data[ii].period;
                           besth0 = h0;
                           bestR = R;
                           bestProb = prob;
                           bestproberrcode = proberrcode;
                        }
                        free_templateStruct(template);
                        template = NULL;
                        ihsCandidates->data[ii].period /= (REAL8)jj;
                     } /* longer period harmonics */
                  } /* shift by harmonics */
                  
                  //Shift period by fractions
                  for (jj=0; jj<10; jj++) {
                     REAL8 periodfact = (jj+1.0)/(jj+2.0);
                     if ( periodfact*ihsCandidates->data[ii].period > minPeriod(ihsCandidates->data[ii].moddepth, inputParams->Tcoh) && periodfact*ihsCandidates->data[ii].period>2.0*3600.0) {
                        ihsCandidates->data[ii].period *= periodfact;
                        template = new_templateStruct(inputParams->templatelength);
                        if (template==NULL) {
                           fprintf(stderr,"%s: new_templateStruct(%d) failed.\n", fn, inputParams->templatelength);
                           XLAL_ERROR(fn, XLAL_EFUNC); 
                        }
                        makeTemplateGaussians(template, ihsCandidates->data[ii], inputParams);
                        if (xlalErrno!=0) {
                           fprintf(stderr,"%s: makeTemplateGaussians() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        }
                        R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
                        prob = (probR(template, aveNoise, aveTFnoisePerFbinRatio, R, &proberrcode));
                        if (XLAL_IS_REAL8_FAIL_NAN(R)) {
                           fprintf(stderr,"%s: calculateR() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        } else if (XLAL_IS_REAL8_FAIL_NAN(prob)) {
                           fprintf(stderr,"%s: probR() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        }
                        h0 = 2.9569*pow(R/(inputParams->Tcoh*inputParams->Tobs),0.25);
                        if (R>farval->far && prob < bestProb) {
                           bestPeriod = ihsCandidates->data[ii].period;
                           besth0 = h0;
                           bestR = R;
                           bestProb = prob;
                           bestproberrcode = proberrcode;
                        }
                        free_templateStruct(template);
                        template = NULL;
                        ihsCandidates->data[ii].period /= periodfact;
                     } /* shift to shorter period */
                     periodfact = 1.0/periodfact;
                     if ( periodfact*ihsCandidates->data[ii].period < 0.2*inputParams->Tobs ) {
                        ihsCandidates->data[ii].period *= periodfact;
                        template = new_templateStruct(inputParams->templatelength);
                        if (template==NULL) {
                           fprintf(stderr,"%s: new_templateStruct(%d) failed.\n", fn, inputParams->templatelength);
                           XLAL_ERROR(fn, XLAL_EFUNC); 
                        }
                        makeTemplateGaussians(template, ihsCandidates->data[ii], inputParams);
                        if (xlalErrno!=0) {
                           fprintf(stderr,"%s: makeTemplateGaussians() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        }
                        R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
                        prob = (probR(template, aveNoise, aveTFnoisePerFbinRatio, R, &proberrcode));
                        if (XLAL_IS_REAL8_FAIL_NAN(R)) {
                           fprintf(stderr,"%s: calculateR() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        } else if (XLAL_IS_REAL8_FAIL_NAN(prob)) {
                           fprintf(stderr,"%s: probR() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        }
                        h0 = 2.9569*pow(R/(inputParams->Tcoh*inputParams->Tobs),0.25);
                        if (R>farval->far && prob < bestProb) {
                           bestPeriod = ihsCandidates->data[ii].period;
                           besth0 = h0;
                           bestR = R;
                           bestProb = prob;
                           bestproberrcode = proberrcode;
                        }
                        free_templateStruct(template);
                        template = NULL;
                        ihsCandidates->data[ii].period /= periodfact;
                     } /* shift to longer period */
                  } /* for jj < 10 (period fractions) */
               } /* if bestR != 0.0 */
               
               free_farStruct(farval);
               farval = NULL;
               
               if (bestR != 0.0) {
                  //If a better period was found, then make sure to save it
                  if (bestPeriod != 0.0) ihsCandidates->data[ii].period = bestPeriod;
                  
                  if (Rfirst > initialFAR && bestProb < gaussCandidates1->data[gaussCandidates1->numofcandidates-1].prob) {
                     //free_candidate(gaussCandidates1[numofcandidates2-1]);
                     //gaussCandidates1[numofcandidates2-1] = NULL;
                     //gaussCandidates1[numofcandidates2-1] = new_candidate();
                     //loadCandidateData(gaussCandidates1[numofcandidates2-1], ihsCandidates[ii]->fsig, ihsCandidates[ii]->period, ihsCandidates[ii]->moddepth, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, bestR, besth0, bestProb, bestproberrcode, ihsCandidates[ii]->normalization);
                     loadCandidateData(&gaussCandidates1->data[gaussCandidates1->numofcandidates-1], ihsCandidates->data[ii].fsig, ihsCandidates->data[ii].period, ihsCandidates->data[ii].moddepth, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, bestR, besth0, bestProb, bestproberrcode, ihsCandidates->data[ii].normalization);
                  } else if (Rfirst <= initialFAR && bestProb != 1.0) {
                     //gaussCandidates1[numofcandidates2] = new_candidate();
                     if (gaussCandidates1->numofcandidates == gaussCandidates1->length-1) gaussCandidates1 = resize_candidateVector(gaussCandidates1, 2*gaussCandidates1->length);
                     loadCandidateData(&gaussCandidates1->data[gaussCandidates1->numofcandidates], ihsCandidates->data[ii].fsig, ihsCandidates->data[ii].period, ihsCandidates->data[ii].moddepth, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, bestR, besth0, bestProb, bestproberrcode, ihsCandidates->data[ii].normalization);
                     gaussCandidates1->numofcandidates++;
                  }
               } /* if bestR != 0.0, add candidate or replace if something better is found */
            } /* if within boundaries */
         } /* for ii < numofcandidates */
         fprintf(LOG,"Initial stage done with candidates = %d\n",gaussCandidates1->numofcandidates);
         fprintf(stderr,"Initial stage done with candidates = %d\n",gaussCandidates1->numofcandidates);
         
         for (ii=0; ii<(INT4)gaussCandidates1->numofcandidates; ii++) fprintf(stderr, "Candidate %d: f0=%g, P=%g, df=%g\n", ii, gaussCandidates1->data[ii].fsig, gaussCandidates1->data[ii].period, gaussCandidates1->data[ii].moddepth);
      } /* if IHSonly is not given */
      else {
         if (exactCandidates2->length < ihsCandidates->numofcandidates) exactCandidates2 = resize_candidateVector(exactCandidates2, ihsCandidates->numofcandidates);
         for (ii=0; ii<(INT4)ihsCandidates->numofcandidates; ii++) {
            loadCandidateData(&exactCandidates2->data[ii], ihsCandidates->data[ii].fsig, ihsCandidates->data[ii].period, ihsCandidates->data[ii].moddepth, ihsCandidates->data[ii].ra, ihsCandidates->data[ii].dec, ihsCandidates->data[ii].stat, 0.0, 0.0, 0, ihsCandidates->data[ii].normalization);
         }
      } /* if IHSonly is given */

////////End of the Gaussian template search

      //Reset IHS candidates, but keep length the same (doesn't reset actual values in the vector)
      ihsCandidates->numofcandidates = 0;
      
      //Search the IHS templates further if user has not specified IHSonly flag
      if (!args_info.IHSonly_given) {
////////Start clustering! Note that the clustering algorithm takes care of the period range of parameter space
         clusterCandidates(gaussCandidates2, gaussCandidates1, ffdata, inputParams, aveNoise, aveTFnoisePerFbinRatio, 0);
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
         REAL4 tcohfactor = 1.49e-3*inputParams->Tcoh + 1.76;
         for (ii=0; ii<(INT4)gaussCandidates2->numofcandidates; ii++) {
            
            REAL8Vector *trialf, *trialb, *trialp;
            REAL8 minf, maxf, minb, maxb;
            UINT4 numf, numb, nump;
            
            //Set up parameters of modulation depth search
            minb = gaussCandidates2->data[ii].moddepth-1.0/inputParams->Tcoh;
            maxb = gaussCandidates2->data[ii].moddepth+1.0/inputParams->Tcoh;
            if (minb<0.25/inputParams->Tcoh) minb = 0.5/inputParams->Tcoh;
            numb = (UINT4)round(2*(maxb-minb)*inputParams->Tcoh)+1;
            trialb = XLALCreateREAL8Vector(numb);
            if (trialb==NULL) {
               fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", fn, numb);
               XLAL_ERROR(fn, XLAL_EFUNC);
            }
            for (jj=0; jj<(INT4)numb; jj++) trialb->data[jj] = minb + 0.5*jj/inputParams->Tcoh;
            //trialb = XLALCreateREAL4Vector(1);
            //trialb->data[0] = gaussCandidates2[ii]->moddepth;
            
            //Set up parameters of signal frequency search
            minf = gaussCandidates2->data[ii].fsig-1.0/inputParams->Tcoh;
            maxf = gaussCandidates2->data[ii].fsig+1.0/inputParams->Tcoh;
            if (minf<inputParams->fmin) minf = inputParams->fmin;
            if (maxf>inputParams->fmin+inputParams->fspan) maxf = inputParams->fmin+inputParams->fspan;
            numf = (UINT4)round(2*(maxf-minf)*inputParams->Tcoh)+1;
            trialf = XLALCreateREAL8Vector(numf);
            if (trialf==NULL) {
               fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", fn, numf);
               XLAL_ERROR(fn, XLAL_EFUNC);
            }
            for (jj=0; jj<(INT4)numf; jj++) trialf->data[jj] = minf + 0.5*jj/inputParams->Tcoh;
            /* minf = gaussCandidates2[ii]->fsig-8.0/inputParams->Tcoh;
            maxf = gaussCandidates2[ii]->fsig+8.0/inputParams->Tcoh;
            if (minf<inputParams->fmin) minf = inputParams->fmin;
            if (maxf>inputParams->fmin+inputParams->fspan) maxf = inputParams->fmin+inputParams->fspan;
            numf = (UINT4)roundf(4*(maxf-minf)*inputParams->Tcoh)+1;
            trialf = XLALCreateREAL4Vector(numf);
            for (jj=0; jj<(INT4)numf; jj++) trialf->data[jj] = minf + 0.25*jj/inputParams->Tcoh; */
            
            //Search over 9 different periods
            nump = 9;
            trialp = XLALCreateREAL8Vector(nump);
            if (trialp==NULL) {
               fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", fn, nump);
               XLAL_ERROR(fn, XLAL_EFUNC);
            }
            
            //FILE *Rtemplatevals = fopen("./Rtemplatevals.dat","w");
            
            //Now search over the parameter space. Frequency, then modulation depth, then period
            //Initialze best values as the initial point we are searching around
            INT4 bestproberrcode;
            REAL8 bestf, bestp, bestdf, bestR, besth0, bestProb;
            bestf = gaussCandidates2->data[ii].fsig;
            bestp = gaussCandidates2->data[ii].period;
            bestdf = gaussCandidates2->data[ii].moddepth;
            bestR = gaussCandidates2->data[ii].stat;
            besth0 = gaussCandidates2->data[ii].h0;
            bestProb = gaussCandidates2->data[ii].prob;
            bestproberrcode = gaussCandidates2->data[ii].proberrcode;
            for (jj=0; jj<(INT4)trialf->length; jj++) {
               for (kk=0; kk<(INT4)trialb->length; kk++) {
                  //Start with period of the first guess, then determine nearest neighbor from the
                  //modulation depth amplitude to find the other period guesses. These parameters 
                  //are determined from simulation to scale the N.N. distance w.r.t. mod. depth with
                  //20% mismatch parameter
                  INT4 midposition = (INT4)((nump-1)*0.5);
                  trialp->data[midposition] = gaussCandidates2->data[ii].period;
                  for (ll=0; ll<midposition; ll++) {
                     REAL8 nnp = trialp->data[midposition+ll]*trialp->data[midposition+ll]*(1+trialp->data[midposition+ll]/tcohfactor/inputParams->Tobs)/tcohfactor/inputParams->Tobs*sqrt(3.6e-3/trialb->data[kk]);
                     trialp->data[midposition+(ll+1)] = trialp->data[midposition+ll] + nnp;
                     nnp = trialp->data[midposition-ll]*trialp->data[midposition-ll]*(1+trialp->data[midposition-ll]/tcohfactor/inputParams->Tobs)/tcohfactor/inputParams->Tobs*sqrt(3.6e-3/trialb->data[kk]);
                     trialp->data[midposition-(ll+1)] = trialp->data[midposition-ll] - nnp;
                  }
                  
                  //Take the mean period and compute a template/FAR pair.
                  REAL8 tempP = calcMeanD(trialp);
                  candidate cand;
                  loadCandidateData(&cand, trialf->data[jj], tempP, trialb->data[kk], (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, 0, 0, 0.0, 0, 0.0);
                  templateStruct *template = new_templateStruct(inputParams->templatelength);
                  if (template==NULL) {
                     fprintf(stderr,"%s: new_templateStruct(%d) failed.\n", fn, inputParams->templatelength);
                     XLAL_ERROR(fn, XLAL_EFUNC); 
                  }
                  makeTemplateGaussians(template, cand, inputParams);
                  if (xlalErrno!=0) {
                     fprintf(stderr,"%s: makeTemplateGaussians() failed.\n", fn);
                     XLAL_ERROR(fn, XLAL_EFUNC);
                  }
                  farStruct *farval = new_farStruct();
                  if (farval==NULL) {
                     fprintf(stderr,"%s: new_farStruct() failed.\n", fn);
                     XLAL_ERROR(fn, XLAL_EFUNC); 
                  }
                  numericFAR(farval, template, templatefarthresh, aveNoise, aveTFnoisePerFbinRatio);
                  if (xlalErrno!=0) {
                     fprintf(stderr,"%s: numericFAR() failed.\n", fn);
                     XLAL_ERROR(fn, XLAL_EFUNC);
                  }
                  free_templateStruct(template);
                  template = NULL;
                  for (ll=0; ll<(INT4)trialp->length; ll++) {
                     if ( trialf->data[jj]-trialb->data[kk]-6/inputParams->Tcoh > inputParams->fmin && trialf->data[jj]+trialb->data[kk]+6/inputParams->Tcoh < inputParams->fmin+inputParams->fspan && trialb->data[kk] < maxModDepth(trialp->data[ll], inputParams->Tcoh) && trialp->data[ll] > minPeriod(trialb->data[kk], inputParams->Tcoh) && inputParams->Tobs/trialp->data[ll] > 5.0 && trialp->data[ll] >= 2.0*3600.0) {
                        loadCandidateData(&cand, trialf->data[jj], trialp->data[ll], trialb->data[kk], (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, 0, 0, 0.0, 0, 0.0);
                        template = new_templateStruct(inputParams->templatelength);
                        if (template==NULL) {
                           fprintf(stderr,"%s: new_templateStruct(%d) failed.\n", fn, inputParams->templatelength);
                           XLAL_ERROR(fn, XLAL_EFUNC); 
                        }
                        makeTemplateGaussians(template, cand, inputParams);
                        if (xlalErrno!=0) {
                           fprintf(stderr,"%s: makeTemplateGaussians() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        }
                        REAL8 R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
                        REAL8 prob = (probR(template, aveNoise, aveTFnoisePerFbinRatio, R, &proberrcode));
                        if (XLAL_IS_REAL8_FAIL_NAN(R)) {
                           fprintf(stderr,"%s: calculateR() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        } else if (XLAL_IS_REAL8_FAIL_NAN(prob)) {
                           fprintf(stderr,"%s: probR() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        }
                        REAL8 h0 = 2.9569*pow(R/(inputParams->Tcoh*inputParams->Tobs),0.25);
                        //fprintf(Rtemplatevals,"%.9g %.9g %.9g %.9g %.9g %.9g\n",trialf->data[jj], trialp->data[ll], trialb->data[kk], R, h0, prob);
                        //if (ll==1) fprintf(stderr,"%f %g %g\n",trialf->data[jj],R,snr);
                        if (R > farval->far && prob < bestProb) {
                           bestf = trialf->data[jj];
                           bestp = trialp->data[ll];
                           bestdf = trialb->data[kk];
                           bestR = R;
                           besth0 = h0;
                           bestProb = prob;
                           bestproberrcode = proberrcode;
                        }
                        free_templateStruct(template);
                        template = NULL;
                     } /* if within boundaries */
                  } /* for ll < trialp */
                  
                  free_farStruct(farval);
                  farval = NULL;
               } /* for kk < trialb */
            } /* for jj < trialf */
            
            //fclose(Rtemplatevals);
            
            if (bestf!=0.0) {
               if (gaussCandidates3->numofcandidates == gaussCandidates3->length-1) gaussCandidates3 = resize_candidateVector(gaussCandidates3, 2*gaussCandidates3->length);
               loadCandidateData(&gaussCandidates3->data[gaussCandidates3->numofcandidates], bestf, bestp, bestdf, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, bestR, besth0, bestProb, bestproberrcode, gaussCandidates2->data[0].normalization);
               gaussCandidates3->numofcandidates++;
            } else {
               fprintf(stderr,"WTF?!\n");
            }
            
            XLALDestroyREAL8Vector(trialf);
            XLALDestroyREAL8Vector(trialb);
            XLALDestroyREAL8Vector(trialp);
            trialf = NULL;
            trialb = NULL;
            trialp = NULL;
         } /* for ii < numofcandidates */
          
         for (ii=0; ii<(INT4)gaussCandidates3->numofcandidates; ii++) fprintf(stderr,"Candidate %d: f0=%g, P=%g, df=%g\n", ii, gaussCandidates3->data[ii].fsig, gaussCandidates3->data[ii].period, gaussCandidates3->data[ii].moddepth);
////////End detailed Gaussian template search

         //Reset 2nd round of Gaussian template candidates, but keep length the same (doesn't reset actual values in the vector)
         gaussCandidates2->numofcandidates = 0;

////////Start clustering!
         clusterCandidates(gaussCandidates4, gaussCandidates3, ffdata, inputParams, aveNoise, aveTFnoisePerFbinRatio, 0);
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
               makeTemplate(template, gaussCandidates4->data[ii], inputParams, secondFFTplan);
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
            
            numericFAR(farval, template, templatefarthresh, aveNoise, aveTFnoisePerFbinRatio);
            if (xlalErrno!=0) {
               fprintf(stderr,"%s: numericFAR() failed.\n", fn);
               XLAL_ERROR(fn, XLAL_EFUNC);
            }
            
            REAL8 R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
            if (XLAL_IS_REAL8_FAIL_NAN(R)) {
               fprintf(stderr,"%s: calculateR() failed.\n", fn);
               XLAL_ERROR(fn, XLAL_EFUNC);
            }
            REAL8 h0 = 2.9569*pow(R/(inputParams->Tcoh*inputParams->Tobs),0.25);
            if (R > farval->far) {
               if (exactCandidates1->numofcandidates == exactCandidates1->length-1) exactCandidates1 = resize_candidateVector(exactCandidates1, 2*exactCandidates1->length);
               loadCandidateData(&exactCandidates1->data[exactCandidates1->numofcandidates], gaussCandidates4->data[ii].fsig, gaussCandidates4->data[ii].period, gaussCandidates4->data[ii].moddepth, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, R, h0, 0.0, 0, gaussCandidates4->data[ii].normalization);
               exactCandidates1->numofcandidates++;
            }
            
            free_templateStruct(template);
            template = NULL;
            free_farStruct(farval);
            farval = NULL;
         } /* for ii < numofcandidates */
         fprintf(LOG, "Number of candidates confirmed with exact templates = %d\n", exactCandidates1->numofcandidates);
         fprintf(stderr, "Number of candidates confirmed with exact templates = %d\n", exactCandidates1->numofcandidates);
////////Done with initial check
         
         //Reset 4th set of Gaussian template candidates, but keep length the same (doesn't reset actual values in the vector)
         gaussCandidates4->numofcandidates = 0;

////////Start detailed "exact" template search!
         for (ii=0; ii<(INT4)exactCandidates1->numofcandidates; ii++) {
         
            REAL8Vector *trialf, *trialb, *trialp;
            REAL8 minf, maxf, minb, maxb;
            UINT4 numf, numb, nump;
            
            minb = exactCandidates1->data[ii].moddepth-1.0/inputParams->Tcoh;
            maxb = exactCandidates1->data[ii].moddepth+1.0/inputParams->Tcoh;
            if (minb<0.25/inputParams->Tcoh) minb = 0.5/inputParams->Tcoh;
            numb = (UINT4)round(2*(maxb-minb)*inputParams->Tcoh)+1;
            trialb = XLALCreateREAL8Vector(numb);
            if (trialb==NULL) {
               fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", fn, numb);
               XLAL_ERROR(fn, XLAL_EFUNC);
            }
            for (jj=0; jj<(INT4)numb; jj++) trialb->data[jj] = minb + 0.5*jj/inputParams->Tcoh;

            minf = exactCandidates1->data[ii].fsig-1.0/inputParams->Tcoh;
            maxf = exactCandidates1->data[ii].fsig+1.0/inputParams->Tcoh;
            if (minf<inputParams->fmin) minf = inputParams->fmin;
            if (maxf>inputParams->fmin+inputParams->fspan) maxf = inputParams->fmin+inputParams->fspan;
            numf = (UINT4)round(2*(maxf-minf)*inputParams->Tcoh)+1;
            trialf = XLALCreateREAL8Vector(numf);
            if (trialb==NULL) {
               fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", fn, numf);
               XLAL_ERROR(fn, XLAL_EFUNC);
            }
            for (jj=0; jj<(INT4)numf; jj++) trialf->data[jj] = minf + 0.5*jj/inputParams->Tcoh;
            
            //This time only use 5 period templates
            nump = 5;
            trialp = XLALCreateREAL8Vector(nump);
            if (trialb==NULL) {
               fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", fn, nump);
               XLAL_ERROR(fn, XLAL_EFUNC);
            }
            
            //Same as before
            INT4 bestproberrcode;
            REAL8 bestf, bestp, bestdf, bestR, besth0, bestProb;
            bestf = exactCandidates1->data[ii].fsig;
            bestp = exactCandidates1->data[ii].period;
            bestdf = exactCandidates1->data[ii].moddepth;
            bestR = exactCandidates1->data[ii].stat;
            besth0 = exactCandidates1->data[ii].h0;
            bestProb = exactCandidates1->data[ii].prob;
            bestproberrcode = exactCandidates1->data[ii].proberrcode;
            for (jj=0; jj<(INT4)trialf->length; jj++) {
               for (kk=0; kk<(INT4)trialb->length; kk++) {
                  INT4 midposition = (INT4)((nump-1)*0.5);
                  trialp->data[midposition] = exactCandidates1->data[ii].period;
                  for (ll=0; ll<midposition; ll++) {
                     REAL8 nnp = trialp->data[midposition+ll]*trialp->data[midposition+ll]*(1+trialp->data[midposition+ll]/tcohfactor/inputParams->Tobs)/tcohfactor/inputParams->Tobs*sqrt(3.6e-3/trialb->data[kk]);
                     trialp->data[midposition+(ll+1)] = trialp->data[midposition+ll] + nnp;
                     nnp = trialp->data[midposition-ll]*trialp->data[midposition-ll]*(1+trialp->data[midposition-ll]/tcohfactor/inputParams->Tobs)/tcohfactor/inputParams->Tobs*sqrt(3.6e-3/trialb->data[kk]);
                     trialp->data[midposition-(ll+1)] = trialp->data[midposition-ll] - nnp;
                  }
                  
                  REAL8 tempP = calcMeanD(trialp);
                  candidate cand;
                  loadCandidateData(&cand, trialf->data[jj], tempP, trialb->data[kk], (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, 0, 0, 0.0, 0, 0.0);
                  templateStruct *template = new_templateStruct(inputParams->templatelength);
                  if (template==NULL) {
                     fprintf(stderr,"%s: new_templateStruct(%d) failed.\n", fn, inputParams->templatelength);
                     XLAL_ERROR(fn, XLAL_EFUNC); 
                  }
                  if (!args_info.gaussTemplatesOnly_given) {
                     makeTemplate(template, cand, inputParams, secondFFTplan);
                     if (xlalErrno!=0) {
                        fprintf(stderr,"%s: makeTemplate() failed.\n", fn);
                        XLAL_ERROR(fn, XLAL_EFUNC);
                     }
                  } else {
                     makeTemplateGaussians(template, cand, inputParams);
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
                  numericFAR(farval, template, templatefarthresh, aveNoise, aveTFnoisePerFbinRatio);
                  if (xlalErrno!=0) {
                     fprintf(stderr,"%s: numericFAR() failed.\n", fn);
                     XLAL_ERROR(fn, XLAL_EFUNC);
                  }
                  free_templateStruct(template);
                  template = NULL;
                  for (ll=0; ll<(INT4)trialp->length; ll++) {
                     if ( trialf->data[jj]-trialb->data[kk]-6/inputParams->Tcoh > inputParams->fmin && trialf->data[jj]+trialb->data[kk]+6/inputParams->Tcoh < inputParams->fmin+inputParams->fspan && trialb->data[kk]<maxModDepth(trialp->data[ll], inputParams->Tcoh) && trialp->data[ll] > minPeriod(trialb->data[kk], inputParams->Tcoh) && inputParams->Tobs/trialp->data[ll]>5 && trialp->data[ll] >= 2.0*3600.0) {
                        
                        loadCandidateData(&cand, trialf->data[jj], trialp->data[ll], trialb->data[kk], (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, 0, 0, 0.0, 0, 0.0);
                        template = new_templateStruct(inputParams->templatelength);
                        if (template==NULL) {
                           fprintf(stderr,"%s: new_templateStruct(%d) failed.\n", fn, inputParams->templatelength);
                           XLAL_ERROR(fn, XLAL_EFUNC); 
                        }
                        if (!args_info.gaussTemplatesOnly_given) {
                           makeTemplate(template, cand, inputParams, secondFFTplan);
                           if (xlalErrno!=0) {
                              fprintf(stderr,"%s: makeTemplate() failed.\n", fn);
                              XLAL_ERROR(fn, XLAL_EFUNC);
                           }
                        } else {
                           makeTemplateGaussians(template, cand, inputParams);
                           if (xlalErrno!=0) {
                              fprintf(stderr,"%s: makeTemplateGaussians() failed.\n", fn);
                              XLAL_ERROR(fn, XLAL_EFUNC);
                           }
                        }
                        
                        REAL8 R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
                        REAL8 prob = (probR(template, aveNoise, aveTFnoisePerFbinRatio, R, &proberrcode));
                        if (XLAL_IS_REAL8_FAIL_NAN(R)) {
                           fprintf(stderr,"%s: calculateR() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        } else if (XLAL_IS_REAL8_FAIL_NAN(prob)) {
                           fprintf(stderr,"%s: probR() failed.\n", fn);
                           XLAL_ERROR(fn, XLAL_EFUNC);
                        }
                        REAL8 h0 = 2.9569*pow(R/(inputParams->Tcoh*inputParams->Tobs),0.25);
                        
                        if (R > farval->far && prob < bestProb) {
                           bestf = trialf->data[jj];
                           bestp = trialp->data[ll];
                           bestdf = trialb->data[kk];
                           bestR = R;
                           besth0 = h0;
                           bestProb = prob;
                           bestproberrcode = proberrcode;
                           
                           /* FILE *TEMPLATEOUT = fopen("./templatevalues.dat","w");
                           INT4 mm;
                           for (mm=0; mm<(INT4)template->templatedata->length; mm++) fprintf(TEMPLATEOUT,"%g %g %g %d %d %d\n",template->templatedata->data[mm],ffdata->ffdata->data[template->pixellocations->data[mm]],aveNoise->data[template->secondfftfrequencies->data[mm]]*aveTFnoisePerFbinRatio->data[template->firstfftfrequenciesofpixels->data[mm]],template->pixellocations->data[mm],template->firstfftfrequenciesofpixels->data[mm],template->secondfftfrequencies->data[mm]);
                           fclose(TEMPLATEOUT); */
                        }
                        free_templateStruct(template);
                        template = NULL;
                     } /* if within boundaries */
                  } /* for ll < trialp */
                  
                  free_farStruct(farval);
                  farval = NULL;
               } /* for kk < trialb */
            } /* for jj < trialf */
            
            //Load candidate
            if (exactCandidates2->numofcandidates == exactCandidates2->length-1) exactCandidates2 = resize_candidateVector(exactCandidates2, 2*exactCandidates2->length);
            loadCandidateData(&exactCandidates2->data[exactCandidates2->numofcandidates], bestf, bestp, bestdf, (REAL4)dopplerpos.Alpha, (REAL4)dopplerpos.Delta, bestR, besth0, bestProb, bestproberrcode, exactCandidates1->data[0].normalization);
            
            exactCandidates2->numofcandidates++;
            
            //Destroy parameter space values
            XLALDestroyREAL8Vector(trialf);
            XLALDestroyREAL8Vector(trialb);
            XLALDestroyREAL8Vector(trialp);
            trialf = NULL;
            trialb = NULL;
            trialp = NULL;
         } /* for ii < numofcandidates */
         fprintf(LOG,"Exact step is done with the total number of candidates = %d\n", exactCandidates2->numofcandidates);
         fprintf(stderr,"Exact step is done with the total number of candidates = %d\n", exactCandidates2->numofcandidates);
////////End of detailed search
         
         //Reset first round of exact template candidates, but keep length the same (doesn't reset actual values in the vector)
         exactCandidates1->numofcandidates = 0;
         
      }
      
      //Destroy stuff
      //XLALDestroyREAL8Vector(aveNoise);
      XLALDestroyREAL4Vector(aveTFnoisePerFbinRatio);
      XLALDestroyREAL4Vector(rmsTFnoisePerFbinRatio);
      
      //Iterate to next sky location
      //XLALNextDopplerSkyPos(&dopplerpos, &scan);
      if ((XLALNextDopplerSkyPos(&dopplerpos, &scan))!=0) {
         fprintf(stderr,"%s: XLALNextDopplerSkyPos() failed.\n", fn);
         XLAL_ERROR(fn, XLAL_EFUNC);
      }
      
   } /* while sky scan is not finished */
   
   if (exactCandidates2->numofcandidates!=0) {
      fprintf(LOG, "\n**Report of candidates:**\n");
      fprintf(stderr, "\n**Report of candidates:**\n");
      
      for (ii=0; ii<(INT4)exactCandidates2->numofcandidates; ii++) {
         fprintf(LOG, "fsig = %g, period = %g, df = %g, RA = %g, DEC = %g, R = %g, h0 = %g, Prob = %g, Error code = %d, FF normalization = %g\n", exactCandidates2->data[ii].fsig, exactCandidates2->data[ii].period, exactCandidates2->data[ii].moddepth, exactCandidates2->data[ii].ra, exactCandidates2->data[ii].dec, exactCandidates2->data[ii].stat, exactCandidates2->data[ii].h0, exactCandidates2->data[ii].prob, exactCandidates2->data[ii].proberrcode, exactCandidates2->data[ii].normalization);
         fprintf(stderr, "fsig = %g, period = %g, df = %g, RA = %g, DEC = %g, R = %g, h0 = %g, Prob = %g, Error code = %d, FF normalization = %g\n", exactCandidates2->data[ii].fsig, exactCandidates2->data[ii].period, exactCandidates2->data[ii].moddepth, exactCandidates2->data[ii].ra, exactCandidates2->data[ii].dec, exactCandidates2->data[ii].stat, exactCandidates2->data[ii].h0, exactCandidates2->data[ii].prob, exactCandidates2->data[ii].proberrcode, exactCandidates2->data[ii].normalization);
      }
   }
   
   //Destroy varaibles
   XLALDestroyREAL4Vector(background);
   XLALDestroyREAL4Vector(usableTFdata);
   XLALDestroyREAL4Vector(detectorVelocities);
   XLALDestroyREAL4Vector(aveNoise);
   XLALDestroyINT4Vector(binshifts);
   free_ffdata(ffdata);
   free_ihsfarStruct(ihsfarstruct);
   free_inputParams(inputParams);
   free_ihsMaxima(ihsmaxima);
   XLALDestroyREAL4FFTPlan(secondFFTplan);
   XLALFree((CHAR*)sft_dir);
   XLALDestroyEphemerisData(edat);
   cmdline_parser_free(&args_info);
   XLALFree(configparams);
   free_candidateVector(ihsCandidates);
   free_candidateVector(gaussCandidates1);
   free_candidateVector(gaussCandidates2);
   free_candidateVector(gaussCandidates3);
   free_candidateVector(gaussCandidates4);
   free_candidateVector(exactCandidates1);
   free_candidateVector(exactCandidates2);
   
   fclose(LOG);
   
   return 0;

}
/*** End of main() ***/



//////////////////////////////////////////////////////////////
// Create new inputParamsStruct  -- done
inputParamsStruct * new_inputParams(void)
{
   
   const CHAR *fn = __func__;
   
   inputParamsStruct *input = XLALMalloc(sizeof(*input));
   if (input==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%lu) failed.", fn, sizeof(*input));
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
   }
   
   return input;

}


//////////////////////////////////////////////////////////////
// Destroy inputParamsStruct  -- done
void free_inputParams(inputParamsStruct *input)
{

   XLALFree((inputParamsStruct*)input);

}



//////////////////////////////////////////////////////////////
// Create a exponentially distributed noise value  -- done
REAL8 expRandNum(REAL8 mu, gsl_rng *ptrToGenerator)
{
   
   const CHAR *fn = __func__;
   
   if (mu<=0.0) {
      fprintf(stderr,"%s: expRandNum(%f, %p) failed.\n", fn, mu, ptrToGenerator);
      XLAL_ERROR_REAL8(fn, XLAL_EINVAL);
   } else if (ptrToGenerator==NULL) {
      fprintf(stderr,"%s: expRandNum(%f, %p) failed.\n", fn, mu, ptrToGenerator);
      XLAL_ERROR_REAL8(fn, XLAL_EFAULT);
   }
   
   REAL8 noise = gsl_ran_exponential(ptrToGenerator, mu);
   
   return noise;

}



//////////////////////////////////////////////////////////////
// Allocate ffdataStruct vectors  -- done
ffdataStruct * new_ffdata(inputParamsStruct *input)
{
   
   const CHAR *fn = __func__;
   
   ffdataStruct *ffdata = XLALMalloc(sizeof(*ffdata));
   if (ffdata==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%lu) failed.\n", fn, sizeof(*ffdata));
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
   }
   
   UINT4 numfbins = (UINT4)(round(input->fspan*input->Tcoh)+1);
   UINT4 numffts = (UINT4)floor(input->Tobs/(input->Tcoh-input->SFToverlap)-1);
   UINT4 numfprbins = (UINT4)floorf(numffts*0.5) + 1;
   
   ffdata->ffdata = XLALCreateREAL4Vector(numfbins * numfprbins);
   if (ffdata->ffdata==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numfbins * numfprbins);
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
   }
   
   ffdata->tfnormalization = ffdata->ffnormalization = 0.0;
   
   return ffdata;
   
}


//////////////////////////////////////////////////////////////
// Destroy ffdataStruct and vectors  -- done
void free_ffdata(ffdataStruct *data)
{

   XLALDestroyREAL4Vector(data->ffdata);
   XLALFree((ffdataStruct*)data);

} 


//////////////////////////////////////////////////////////////
// Read in SFT data to produce a TF vector
REAL4Vector * readInSFTs(inputParamsStruct *input, REAL8 *normalization)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj;
   LALStatus status;
   status.statusPtr = NULL;
   SFTCatalog *catalog = NULL;
   SFTConstraints *constraints = NULL;
   SFTVector *sfts = NULL;
   
   //Find SFT files
   LALSFTdataFind(&status, &catalog, sft_dir, constraints);
   if (status.statusCode != 0) {
      fprintf(stderr,"%s: LAALSFTdataFind() failed with code = %d.\n", fn, status.statusCode);
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
   INT4 sftlength = sfts->data->data->length;
   INT4 nonexistantsft = 0;
   REAL4Vector *tfdata = XLALCreateREAL4Vector((UINT4)(numffts*sftlength));
   if (tfdata==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numffts*sftlength);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   }
   REAL8 sqrtnorm = sqrt(*(normalization));
   for (ii=0; ii<numffts; ii++) {
      
      SFTDescriptor *sftdescription = &(catalog->data[ii - nonexistantsft]);
      if (sftdescription->header.epoch.gpsSeconds == (INT4)(ii*(input->Tcoh-input->SFToverlap)+input->searchstarttime)) {
         SFTtype *sft = &(sfts->data[ii - nonexistantsft]);
         for (jj=0; jj<sftlength; jj++) {
            COMPLEX8 sftcoeff = sft->data->data[jj];
            tfdata->data[ii*sftlength + jj] = (REAL4)((sqrtnorm*sftcoeff.re)*(sqrtnorm*sftcoeff.re) + (sqrtnorm*sftcoeff.im)*(sqrtnorm*sftcoeff.im));  //power, normalized
         }
      } else {
         for (jj=0; jj<sftlength; jj++) tfdata->data[ii*sftlength + jj] = 0.0;   //Set values to be zero
         nonexistantsft++;    //increment the nonexistantsft counter
      }
      
   }
   
   LALDestroySFTCatalog(&status, &catalog);
   XLALDestroySFTVector(sfts);
   
   fprintf(stderr,"TF before weighting, mean subtraction = %g\n",calcMean(tfdata));
   
   return tfdata;

}




//////////////////////////////////////////////////////////////
// Slide SFT TF data  -- 
void slideTFdata(REAL4Vector *output, inputParamsStruct *input, REAL4Vector *tfdata, INT4Vector *binshifts)
{
   
   INT4 ii, jj;
   INT4 numffts = (INT4)floor(input->Tobs/(input->Tcoh-input->SFToverlap)-1);
   INT4 numfbins = (INT4)(round(input->fspan*input->Tcoh)+1);
   
   for (ii=0; ii<numffts; ii++) {
      for (jj=0; jj<numfbins; jj++) output->data[ii*numfbins + jj] = tfdata->data[ii*(numfbins+2*input->maxbinshift) + jj + input->maxbinshift + binshifts->data[ii]];
   }
   
}



//////////////////////////////////////////////////////////////
// Determine the TF running mean of each SFT  -- 
void tfRngMeans(REAL4Vector *output, REAL4Vector *tfdata, INT4 numffts, INT4 numfbins, INT4 blksize)
{
   
   const CHAR *fn = __func__;
   
   LALStatus status;
   status.statusPtr = NULL;
   REAL8 bias;
   INT4 ii, jj;
   
   //Blocksize of running median. This needs to be >=501 bins because of accuracy in
   //determination of background is essential for 2nd PSD computation.
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
   REAL8 invbias = 1.0/bias;
   
   REAL4Vector *inpsd = XLALCreateREAL4Vector((UINT4)(numfbins+blksize-1));
   REAL4Vector *mediansout = XLALCreateREAL4Vector((UINT4)numfbins);
   if (inpsd==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (UINT4)(numfbins+blksize-1));
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (mediansout==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numfbins);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   for (ii=0; ii<numffts; ii++) {
      //If the SFT values were not zero, then compute the running median
      if (tfdata->data[ii*(numfbins+blksize-1)]!=0.0) {
         //Determine running median value, convert to mean value
         //for (jj=0; jj<(INT4)inpsd->length; jj++) inpsd->data[jj] = tfdata->data[ii*(numfbins+blksize-1) + jj];
         memcpy(inpsd->data, &(tfdata->data[ii*(numfbins+blksize-1)]), sizeof(REAL4)*inpsd->length);
         
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
   }
   
   fprintf(stderr,"Mean of running means = %g\n",calcMean(output));
   
   XLALDestroyREAL4Vector(inpsd);
   XLALDestroyREAL4Vector(mediansout);

}



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
      
      //Get sum of antenna pattern weight/variances for each frequency bin as a function of time (only for existant SFTs
      REAL8 sumofweights = 0.0;
      for (jj=0; jj<numffts; jj++) {
         if (rngMeanssq->data[jj] != 0.0) sumofweights += antweightssq->data[jj]/rngMeanssq->data[jj];
      }
      REAL8 invsumofweights = 1.0/sumofweights;
      
      //Now do mean subtraction, noise weighting, antenna pattern weighting
      for (jj=0; jj<numffts; jj++) {
         if (rngMeans->data[jj*numfbins+ii] != 0.0) {
            output->data[jj*numfbins+ii] = (REAL4)(invsumofweights*antPatternWeights->data[jj]*(tfdata->data[jj*numfbins+ii] - rngMeans->data[jj*numfbins+ii])/(rngMeanssq->data[jj]));
         } else {
            output->data[jj*numfbins+ii] = 0.0;
         }
      }
   }
   
   XLALDestroyREAL4Vector(antweightssq);
   XLALDestroyREAL4Vector(rngMeanssq);
   
   fprintf(stderr,"TF after weighting, mean subtraction = %g\n",calcMean(output));

}


//////////////////////////////////////////////////////////////
// Make the second FFT powers
void makeSecondFFT(REAL4Vector *output, REAL8 normalization, REAL4Vector *tfdata, inputParamsStruct *input, REAL4FFTPlan *plan)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj;
   
   INT4 numffts = (INT4)floor(input->Tobs/(input->Tcoh-input->SFToverlap)-1);    //Number of FFTs
   INT4 numfbins = (INT4)(round(input->fspan*input->Tcoh)+1);    //Number of frequency bins
   
   //Do the second FFT
   REAL4Vector *x = XLALCreateREAL4Vector((UINT4)numffts);
   REAL4Window *win = XLALCreateHannREAL4Window(x->length);
   REAL4Vector *psd = XLALCreateREAL4Vector((UINT4)floor(x->length*0.5)+1);
   if (x==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numffts);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (win==NULL) {
      fprintf(stderr,"%s: XLALCreateHannREAL4Window(%d) failed.\n", fn, x->length);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (psd==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (UINT4)floor(x->length*0.5)+1);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   REAL4 winFactor = 8.0/3.0;
   //REAL8 psdfactor = winFactor*0.5*params->Tcoh;
   //REAL8 psdfactor = (REAL8)(winFactor*0.5*params->Tobs*params->Tcoh);
   //REAL8 psdfactor = winFactor/(x->length*x->length);
   REAL4 psdfactor = winFactor;
   //First loop over frequencies
   for (ii=0; ii<numfbins; ii++) {
   
      //Next, loop over times and pick the right frequency bin for each FFT and window
      //for (jj=0; jj<(INT4)x->length; jj++) x->data[jj] = (REAL4)(tfdata->data[ii + jj*numfbins]*win->data->data[jj]);
      x = SSVectorMultiply_with_stride_and_offset(x, tfdata, win->data, numfbins, 1, ii, 0);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: SSVectorMutiply_with_stride_and_offset() failed.\n", fn);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      
      //Make the FFT
      //INT4 check = XLALREAL4PowerSpectrum(psd,x,plan);
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
      for (jj=0; jj<(INT4)psd->length; jj++) output->data[psd->length*ii + jj] = (REAL4)(psd->data[jj]*psdfactor*normalization);
      
   }
   
   XLALDestroyREAL4Vector(x);
   XLALDestroyREAL4Vector(psd);
   XLALDestroyREAL4Window(win);
   
}



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
   }
   
   REAL4 avgTFdata = calcMean(aveNoiseInTime);
   if (XLAL_IS_REAL4_FAIL_NAN(avgTFdata)) {
      fprintf(stderr,"%s: calcMean() failed.\n", fn);
      XLAL_ERROR_REAL4(fn, XLAL_EFUNC);
   }
   
   XLALDestroyREAL4Vector(aveNoiseInTime);
   XLALDestroyREAL4Vector(rngMeansOverBand);
   
   return avgTFdata;
   
}


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
   }
   
   REAL4 rmsTFdata = calcRms(aveNoiseInTime);
   if (XLAL_IS_REAL4_FAIL_NAN(rmsTFdata)) {
      fprintf(stderr,"%s: calcRms() failed.\n", fn);
      XLAL_ERROR_REAL4(fn, XLAL_EFUNC);
   }
   
   XLALDestroyREAL4Vector(aveNoiseInTime);
   XLALDestroyREAL4Vector(rngMeansOverBand);
   
   return rmsTFdata;
   
}


//////////////////////////////////////////////////////////////
// Measure of the average noise power in each 2st FFT frequency bin  -- 
void ffPlaneNoise(REAL4Vector *aveNoise, inputParamsStruct *input, REAL4Vector *backgrnd, REAL4Vector *antweights, REAL8 *normalization)
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
   //srand(time(NULL));
   //UINT8 randseed = rand();
   //gsl_rng_set(rng, randseed);
   gsl_rng_set(rng, 0);
   
   //Set up for making the PSD
   for (ii=0; ii<(INT4)aveNoise->length; ii++) aveNoise->data[ii] = 0.0;
   REAL4Window *win = XLALCreateHannREAL4Window((UINT4)numffts);
   REAL4FFTPlan *plan = XLALCreateForwardREAL4FFTPlan((UINT4)numffts, 0 );
   REAL4Vector *psd = XLALCreateREAL4Vector((UINT4)numfprbins);   //Current PSD calculation
   if (plan==NULL) {
      fprintf(stderr,"%s: XLALCreateForwardREAL4FFTPlan(%d, 0) failed.\n", fn, numffts);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (win==NULL) {
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
         //for (jj=0; jj<(INT4)rngMeansOverBand->length; jj++) rngMeansOverBand->data[jj] = backgrnd->data[ii*numfbins + jj];
         memcpy(rngMeansOverBand->data, &(backgrnd->data[ii*numfbins]), sizeof(*rngMeansOverBand->data)*rngMeansOverBand->length);
         aveNoiseInTime->data[ii] = calcMean(rngMeansOverBand);
         if (XLAL_IS_REAL4_FAIL_NAN(aveNoiseInTime->data[ii])) {
            fprintf(stderr,"%s: calcMean() failed.\n", fn);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         }
         
         sumofweights += (antweights->data[ii]*antweights->data[ii])/(aveNoiseInTime->data[ii]*aveNoiseInTime->data[ii]);
      } else {
         aveNoiseInTime->data[ii] = 0.0;
      }
      
   }
   invsumofweights = 1.0/sumofweights;
   
   //FILE *BACKGRND = fopen("./background.dat","w");
   
   //Load time series of powers, normalize, mean subtract and Hann window
   REAL4Vector *x = XLALCreateREAL4Vector(aveNoiseInTime->length);
   REAL8Vector *multiplicativeFactor = XLALCreateREAL8Vector(x->length);
   //REAL8 psdfactor = winFactor*0.5*param->Tcoh;
   //REAL8 psdfactor = (REAL8)(winFactor*0.5*param->Tobs*param->Tcoh);
   //REAL8 psdfactor = winFactor/(x->length*x->length);
   REAL4 psdfactor = winFactor;
   
   for (ii=0; ii<(INT4)x->length; ii++) {
      if (aveNoiseInTime->data[ii] != 0.0) multiplicativeFactor->data[ii] = win->data->data[ii]*antweights->data[ii]/aveNoiseInTime->data[ii]*invsumofweights;
      else multiplicativeFactor->data[ii] = 0.0;
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
   for (ii=0; ii<200; ii++) {
      for (jj=0; jj<(INT4)x->length; jj++) {
         if (aveNoiseInTime->data[jj] != 0.0) {
            noiseval = expRandNum(aveNoiseInTime->data[jj], rng);
            if (jj==0 || (jj>0 && aveNoiseInTime->data[jj-1] == 0.0)) {
               x->data[jj] = (REAL4)(multiplicativeFactor->data[jj]*(noiseval/aveNoiseInTime->data[jj]-1.0));
            } else {
               REAL8 newnoiseval = (1.0-corrfactorsquared)*noiseval + corrfactorsquared*prevnoiseval;
               REAL8 newavenoise = (1.0-corrfactorsquared)*aveNoiseInTime->data[jj] + corrfactorsquared*aveNoiseInTime->data[jj-1];
               x->data[jj] = (REAL4)(multiplicativeFactor->data[jj]*(newnoiseval/newavenoise-1.0));
            }
            prevnoiseval = noiseval;
         } else {
            x->data[jj] = 0.0;
         }
      }
      
      //Do the FFT
      if ( (XLALREAL4PowerSpectrum(psd, x, plan)) != 0) {
         fprintf(stderr,"%s: XLALREAL4PowerSpectrum() failed.\n", fn);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      
      //Rescale and sum into the bins
      for (jj=0; jj<(INT4)aveNoise->length; jj++) {
         aveNoise->data[jj] += psd->data[jj];
      }
   }
   
   //Average
   for (ii=0; ii<(INT4)aveNoise->length; ii++) aveNoise->data[ii] *= (REAL4)(5.0e-3*psdfactor*(1.0+2.0*corrfactorsquared));
   
   //Fix 0th and end bins (0 only for odd x->length, 0 and end for even x->length)
   if (GSL_IS_EVEN(x->length)==1) {
      aveNoise->data[0] *= 2.0;
      aveNoise->data[aveNoise->length-1] *= 2.0;
   } else {
      aveNoise->data[0] *= 2.0;
   }
   
   //Compute normalization
   *(normalization) = 1.0/calcMean(aveNoise);
   for (ii=0; ii<(INT4)aveNoise->length; ii++) {
      aveNoise->data[ii] *= *(normalization);
      //fprintf(BACKGRND,"%g\n",aveNoise->data[ii]);
   }
   
   //fclose(BACKGRND);

   XLALDestroyREAL4Vector(x);
   XLALDestroyREAL4Vector(psd);
   XLALDestroyREAL4Window(win);
   XLALDestroyREAL4FFTPlan(plan);
   XLALDestroyREAL4Vector(aveNoiseInTime);
   XLALDestroyREAL4Vector(rngMeansOverBand);
   XLALDestroyREAL8Vector(multiplicativeFactor);
   gsl_rng_free(rng);

}



//////////////////////////////////////////////////////////////
// Calculate the R statistic
REAL8 calculateR(REAL4Vector *ffdata, templateStruct *templatestruct, REAL4Vector *noise, REAL4Vector *fbinaveratios)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii;
   
   REAL8 sumofsqweights = 0.0;
   for (ii=0; ii<(INT4)templatestruct->templatedata->length; ii++) if (templatestruct->templatedata->data[ii]!=0.0) sumofsqweights += (templatestruct->templatedata->data[ii]*templatestruct->templatedata->data[ii]);
   if (sumofsqweights==0.0) {
      fprintf(stderr,"%s: Sum of weights squared = 0.\n", fn);
      XLAL_ERROR_REAL8(fn, XLAL_EFPDIV0);
   }
   REAL8 sumofsqweightsinv = 1.0/sumofsqweights;
   
   REAL8 R = 0.0;
   for (ii=0; ii<(INT4)templatestruct->templatedata->length; ii++) {
      if (templatestruct->templatedata->data[ii]!=0.0) {
         R += (ffdata->data[ templatestruct->pixellocations->data[ii] ] - noise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ])*templatestruct->templatedata->data[ii]*sumofsqweightsinv;
         //fprintf(stderr,"%g %g %g %g\n", templatestruct->templatedata->data[ii]*sumofsqweightsinv, ffdata->data[ templatestruct->pixellocations->data[ii] ], noise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ], R);
      }
   }
   
   return R;
   
}



//////////////////////////////////////////////////////////////
// Calculates maximum modulation depth
REAL8 maxModDepth(REAL8 period, REAL8 cohtime)
{

   REAL8 maxB = 0.5*period/cohtime/cohtime;
   
   return maxB;

}


//////////////////////////////////////////////////////////////
// Calculates minimum period allowable for modulation depth and Tcoh
REAL8 minPeriod(REAL8 moddepth, REAL8 cohtime)
{

   REAL8 minP = 2.0*moddepth*cohtime*cohtime;
   
   return minP;

}



//////////////////////////////////////////////////////////////
// Compute the mean value of a vector of values
REAL4 calcMean(REAL4Vector *vector)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii;
   
   double *gslarray = XLALMalloc(sizeof(double)*vector->length);
   if (gslarray==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%d) failed.\n", fn, vector->length);
      XLAL_ERROR_REAL4(fn, XLAL_ENOMEM);
   }
   for (ii=0; ii<(INT4)vector->length; ii++) gslarray[ii] = (double)vector->data[ii];
   REAL4 meanval = (REAL4)gsl_stats_mean(gslarray, 1, vector->length);
   
   XLALFree((double*)gslarray);
   
   return meanval;

}


//////////////////////////////////////////////////////////////
// Compute the standard deviation of a vector of values
REAL4 calcStddev(REAL4Vector *vector)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii;
   
   double *gslarray = XLALMalloc(sizeof(double)*vector->length);
   if (gslarray==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%d) failed.\n", fn, vector->length);
      XLAL_ERROR_REAL4(fn, XLAL_ENOMEM);
   }
   for (ii=0; ii<(INT4)vector->length; ii++) gslarray[ii] = (double)vector->data[ii];
   REAL4 stddev = (REAL4)gsl_stats_sd(gslarray, 1, vector->length);
   
   XLALFree((double*)gslarray);   
   
   return stddev;

}



//////////////////////////////////////////////////////////////
// Compute the RMS of a vector of values
REAL4 calcRms(REAL4Vector *vector)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii;
   REAL8Vector *sqvector = XLALCreateREAL8Vector(vector->length);
   if (sqvector==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", fn, vector->length);
      XLAL_ERROR_REAL4(fn, XLAL_ENOMEM);
   }
   for (ii=0; ii<(INT4)vector->length; ii++) sqvector->data[ii] = (REAL8)(vector->data[ii]*vector->data[ii]);
   REAL4 rms = (REAL4)sqrt(calcMeanD(sqvector));
   
   /* double *gslarray = (double*)XLALMalloc(sizeof(double)*vector->length);
   for (ii=0; ii<(INT4)vector->length; ii++) gslarray[ii] = (double)vector->data[ii];
   REAL4 rms = (REAL4)sqrt(gsl_stats_tss_m(gslarray, 1, vector->length, 0.0)/vector->length); */
   
   XLALDestroyREAL8Vector(sqvector);
   //XLALFree((double*)gslarray);
   
   return rms;

}



//////////////////////////////////////////////////////////////
// Compute the mean value of a vector of REAL8 values
REAL8 calcMeanD(REAL8Vector *vector)
{
   
   REAL8 meanval = gsl_stats_mean((double*)vector->data, 1, vector->length);
   
   return meanval;
   
}


//////////////////////////////////////////////////////////////
// Compute the standard deviation of a vector of REAL8 values
REAL8 calcStddevD(REAL8Vector *vector)
{
     
   REAL8 stddev = gsl_stats_sd((double*)vector->data, 1, vector->length);
   
   return stddev;
   
}



REAL4Vector * SSVectorMultiply_with_stride_and_offset(REAL4Vector *output, REAL4Vector *input1, REAL4Vector *input2, INT4 stride1, INT4 stride2, INT4 offset1, INT4 offset2)
{
   
   const CHAR *fn = __func__;
   
   REAL4 *a;
   REAL4 *b;
   REAL4 *c;
   INT4   n;
   
   if ( ! output || ! input1 || !input2 || ! output->data || ! input1->data || ! input2->data ) {
      XLAL_ERROR_NULL( fn, XLAL_EFAULT );
   } else if ( ! output->length ) {
      XLAL_ERROR_NULL( fn, XLAL_EBADLEN );
   }
   
   a = input1->data + offset1;
   b = input2->data + offset2;
   c = output->data;
   n = output->length;
   
   while (n-- > 0) {
      *c++ = (*a)*(*b);
      a = a + stride1;
      b = b + stride2;
   }
   
   return output;
   
}



