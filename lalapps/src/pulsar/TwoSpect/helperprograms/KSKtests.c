/*
*  Copyright (C) 2013 Evan Goetz
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
#include <gsl/gsl_cdf.h>

#include "../antenna.h"
#include "../TwoSpectTypes.h"
#include "../statistics.h"
#include "../cmdline.h"

//Prototypes
inputParamsStruct * new_inputParams(INT4 numofIFOs);
void free_inputParams(inputParamsStruct *input);
REAL4Vector * readInSFTs(inputParamsStruct *input, REAL8 *normalization);
INT4 readTwoSpectInputParams(inputParamsStruct *params, struct gengetopt_args_info args_info);
INT4 qsort_REAL4_compar(const void *a, const void *b);

//Global variables
CHAR *earth_ephemeris = NULL, *sun_ephemeris = NULL, *sft_dir_file = NULL;
static const LALStatus empty_status;

//Main program
int main(int argc, char *argv[])
{
   
  INT4 ii, jj;               //counter variables
   
   //Turn off gsl error handler
   gsl_set_error_handler_off();
   
   //Initiate command line interpreter and config file loader
   struct gengetopt_args_info args_info;
   struct cmdline_parser_params *configparams;
   configparams = cmdline_parser_params_create();  //initialize parameters structure
   configparams->check_required = 0;  //don't check for required values at the step
   if ( cmdline_parser_ext(argc, argv, &args_info, configparams) ) {
       fprintf(stderr, "%s: cmdline_parser_ext() failed.\n", __func__);
       XLAL_ERROR(XLAL_FAILURE);
   }
   configparams->initialize = 0;  //don't reinitialize the parameters structure
   if ( args_info.config_given && cmdline_parser_config_file(args_info.config_arg, &args_info, configparams) ) {
      fprintf(stderr, "%s: cmdline_parser_config_file() failed.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   //Check required
   if ( cmdline_parser_required(&args_info, argv[0]) ) {
      fprintf(stderr, "%s: cmdline_parser_required() failed.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   
   //Set lalDebugLevel to user input or 0 if no input
   lalDebugLevel = args_info.laldebug_arg;
   
   //Allocate input parameters structure memory
   inputParamsStruct *inputParams = new_inputParams(args_info.IFO_given);
   if (inputParams==NULL) {
      fprintf(stderr, "%s: new_inputParams() failed.\n", __func__);
      XLAL_ERROR(XLAL_EFUNC);
   }
   
   //Read TwoSpect input parameters
   if ( (readTwoSpectInputParams(inputParams, args_info)) != 0 ) {
      fprintf(stderr, "%s: readTwoSpectInputParams() failed.\n", __func__);
      XLAL_ERROR(XLAL_EFUNC);
   }
   
   //Initialize ephemeris data structure
   EphemerisData *edat = XLALInitBarycenter(earth_ephemeris, sun_ephemeris);
   if (edat==NULL) {
      fprintf(stderr, "%s: XLALInitBarycenter() failed.\n", __func__);
      XLAL_ERROR(XLAL_EFUNC);
   }
   
   //Maximum detector velocity in units of c from start of observation time - Tcoh to end of observation + Tcoh
   REAL4 detectorVmax = CompDetectorVmax(inputParams->searchstarttime-inputParams->Tcoh, inputParams->Tcoh, inputParams->SFToverlap, inputParams->Tobs+2.0*inputParams->Tcoh, inputParams->det[0], edat);
   if (xlalErrno!=0) {
      fprintf(stderr, "%s: CompDetectorVmax() failed.\n", __func__);
      XLAL_ERROR(XLAL_EFUNC);
   }
   
   //Assume maximum bin shift possible
   inputParams->maxbinshift = (INT4)round(detectorVmax * (inputParams->fmin+0.5*inputParams->fspan) * inputParams->Tcoh)+1;

   //Read in the T-F data from SFTs
   fprintf(stderr, "Loading in SFTs... ");
   REAL8 tfnormalization = 2.0/inputParams->Tcoh/(args_info.avesqrtSh_arg*args_info.avesqrtSh_arg);
   REAL4Vector *tfdata = readInSFTs(inputParams, &(tfnormalization));
   if (tfdata==NULL) {
      fprintf(stderr, "\n%s: readInSFTs() failed.\n", __func__);
      XLAL_ERROR(XLAL_EFUNC);
   }
   fprintf(stderr, "done\n");
   
   //Removing bad SFTs using K-S test and Kuiper's test
   if (inputParams->markBadSFTs!=0 && inputParams->signalOnly==0) {
      fprintf(stderr, "Marking and removing bad SFTs... ");
      INT4 numffts = (INT4)floor(inputParams->Tobs/(inputParams->Tcoh-inputParams->SFToverlap)-1);    //Number of FFTs
      INT4 numfbins = (INT4)(round(inputParams->fspan*inputParams->Tcoh+2.0*inputParams->dfmax*inputParams->Tcoh)+12+1)+2*inputParams->maxbinshift+inputParams->blksize-1;     //Number of frequency bins
      REAL4Vector *tempvect = XLALCreateREAL4Vector(numfbins);
      if (tempvect==NULL) {
         fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numfbins);
         XLAL_ERROR(XLAL_EFUNC);
      }
      REAL8 ksthreshold = 1.358/(sqrt(numfbins)+0.12+0.11/sqrt(numfbins));
      REAL8 kuiperthreshold = 1.747/(sqrt(numfbins)+0.155+0.24/sqrt(numfbins));
      INT4 badsfts = 0, totalsfts = 0;
      FILE *OUTPUT = fopen("./output/kskoutput.dat","w");
      for (ii=0; ii<numffts; ii++) {
         if (tfdata->data[ii*numfbins]!=0.0) {
            totalsfts++;
            memcpy(tempvect->data, &(tfdata->data[ii*numfbins]), sizeof(REAL4)*tempvect->length);
            qsort(tempvect->data, tempvect->length, sizeof(REAL4), qsort_REAL4_compar);
            REAL4 vector_median = 0.0;
            if (tempvect->length % 2 != 1) vector_median = 0.5*(tempvect->data[(INT4)(0.5*tempvect->length)-1] + tempvect->data[(INT4)(0.5*tempvect->length)]);
            else vector_median = tempvect->data[(INT4)(0.5*tempvect->length)];
            REAL4 vector_mean = (REAL4)(vector_median/LAL_LN2);

            REAL8 ksvalue = 0.0, testval1, testval2, testval;
            REAL8 oneoverlength = 1.0/tempvect->length;
            for (jj=0; jj<(INT4)tempvect->length; jj++) {
               testval1 = fabs((1.0+jj)*oneoverlength - gsl_cdf_exponential_P(tempvect->data[jj], vector_mean));
               testval2 = fabs(jj*oneoverlength - gsl_cdf_exponential_P(tempvect->data[jj], vector_mean));
               testval = fmax(testval1, testval2);
               if (testval>ksvalue) ksvalue = testval;
            }

            REAL8 loval = 0.0, hival = 0.0;
            for (jj=0; jj<(INT4)tempvect->length; jj++) {
              testval1 = (1.0+jj)*oneoverlength - gsl_cdf_exponential_P(tempvect->data[jj], vector_mean);
              testval2 = jj*oneoverlength - gsl_cdf_exponential_P(tempvect->data[jj], vector_mean);
              if (hival<testval1) hival = testval1;
              if (loval<testval2) loval = testval2;
            }
            REAL8 kuiperval1 = hival + loval;

            loval = -1.0, hival = -1.0;
            for (jj=0; jj<(INT4)tempvect->length; jj++) {
              testval1 = (1.0+jj)*oneoverlength - gsl_cdf_exponential_P(tempvect->data[jj], vector_mean);
              testval2 = jj*oneoverlength - gsl_cdf_exponential_P(tempvect->data[jj], vector_mean);
              if (hival<testval1) hival = testval1;
              if (hival<testval2) hival = testval2;
              if (loval<-testval1) loval = -testval1;
              if (loval<-testval2) loval = -testval2;
            }
            REAL8 kuiperval = hival + loval;

            fprintf(OUTPUT, "%g %g %g\n", ksvalue, kuiperval1, kuiperval);

            if (ksvalue>ksthreshold || kuiperval>kuiperthreshold) {
               badsfts++;
            }
         }
      }
      fclose(OUTPUT);
      fprintf(stderr, "Fraction excluded in K-S and Kuiper's tests = %f\n", (REAL4)badsfts/(REAL4)totalsfts);
      XLALDestroyREAL4Vector(tempvect);
   }

   XLALDestroyREAL4Vector(tfdata);
   XLALDestroyEphemerisData(edat);
   cmdline_parser_free(&args_info);
   free_inputParams(inputParams);

   return 0;

}

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

void free_inputParams(inputParamsStruct *input)
{
   
   XLALFree((CHAR*)input->sftType);
   XLALFree((LALDetector*)input->det);
   gsl_rng_free(input->rng);
   XLALFree((inputParamsStruct*)input);

} /* free_inputParams() */

REAL4Vector * readInSFTs(inputParamsStruct *input, REAL8 *normalization)
{
   
   INT4 ii, jj;
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
   catalog = XLALSFTdataFind(sft_dir_file, &constraints);
   if (catalog==NULL) {
      fprintf(stderr,"%s: XLALSFTdataFind() failed.\n", __func__);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   
   //Determine band size (remember to get extra bins because of the running median and the bin shifts due to detector velocity)
   REAL8 minfbin = round(input->fmin*input->Tcoh - input->dfmax*input->Tcoh - 0.5*(input->blksize-1) - (REAL8)(input->maxbinshift) - 6.0)/input->Tcoh;
   REAL8 maxfbin = round((input->fmin + input->fspan)*input->Tcoh + input->dfmax*input->Tcoh + 0.5*(input->blksize-1) + (REAL8)(input->maxbinshift) + 6.0)/input->Tcoh;
   
   //Now extract the data
   SFTVector *sfts = XLALLoadSFTs(catalog, minfbin+0.1/input->Tcoh, maxfbin-0.1/input->Tcoh);
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
   
   fprintf(stderr, "Duty factor = %f\n", 1.0-(REAL4)nonexistantsft/(REAL4)numffts);
   
   return tfdata;

} /* readInSFTs() */

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
   params->signalOnly = args_info.signalOnly_given;                              //SFTs contain only signal, no noise (default = 0)

   //Parameters without default arguments but are required
   //gengetopt has already checked to see they are present and would have failed
   params->Tobs = args_info.Tobs_arg;
   params->Tcoh = args_info.Tcoh_arg;
   params->SFToverlap = args_info.SFToverlap_arg;
   params->searchstarttime = args_info.t0_arg;
   params->fmin = args_info.fmin_arg;
   params->fspan = args_info.fspan_arg;
   params->Pmin = args_info.Pmin_arg;
   params->Pmax = args_info.Pmax_arg;
   params->dfmin = args_info.dfmin_arg;
   params->dfmax = args_info.dfmax_arg;
   
   //Non-default arguments (but nonetheless required in certain circumstances)
   if (args_info.ihsfar_given) params->ihsfar = args_info.ihsfar_arg;            //Incoherent harmonic sum false alarm rate
   else if (!args_info.ihsfar_given && !args_info.templateTest_given) {
      fprintf(stderr, "%s: the IHS FAR must be specified.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   if (args_info.tmplfar_given) params->templatefar = args_info.tmplfar_arg;     //Template false alarm rate
   else if (!args_info.tmplfar_given && args_info.templateTest_given) {
      fprintf(stderr, "%s: the template FAR must be specified.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }

   //Arguments that might not be given, but we should just set them to an invalid value if they are not given
   if (args_info.keepOnlyTopNumIHS_given) params->keepOnlyTopNumIHS = args_info.keepOnlyTopNumIHS_arg;         //Keep only top X IHS candidates
   else params->keepOnlyTopNumIHS = -1;
   if (args_info.lineDetection_given) params->lineDetection = args_info.lineDetection_arg;                     //Line detection
   else params->lineDetection = -1.0;

   //Settings for IHS FOM
   //Exit with error if both or neither is chosen unless we are only doing the template test
   //When the parameters are not given, set them to zero
   if (!args_info.templateTest_given && ((args_info.ihsfomfar_given && args_info.ihsfom_given) || (!args_info.ihsfomfar_given && !args_info.ihsfom_given))) {
      fprintf(stderr, "%s: You must choose only one of the IHS FOM FAR argument or the IHS FOM argument.\n", __func__);
      XLAL_ERROR(XLAL_EINVAL);
   }
   if (args_info.ihsfomfar_given) params->ihsfomfar = args_info.ihsfomfar_arg;   //IHS figure of merit false alarm rate
   else params->ihsfomfar = 0.0;
   if (args_info.ihsfom_given) params->ihsfom = args_info.ihsfom_arg;            //IHS figure of merit threshold value
   else params->ihsfom = 0.0;
   
   //log10(template FAR)
   params->log10templatefar = log10(params->templatefar);
   
   //Blocksize should be an odd number
   if (params->blksize % 2 != 1) params->blksize += 1;
   
   // Warnings when using hidden flags
   if (args_info.signalOnly_given) {
      fprintf(stderr,"WARNING: --signalOnly argument has been specified\n");
   }
   if (args_info.templateTest_given) {
      fprintf(stderr,"WARNING: --templateTest argument has been specified\n");
   }
   if (args_info.ULsolver_arg!=0) {
      fprintf(stderr,"WARNING: --ULsolver = %d instead of the default value of 0\n", args_info.ULsolver_arg);
   }
   if (args_info.dopplerMultiplier_given) {
      fprintf(stderr,"WARNING: --dopplerMultiplier = %g instead of the default value of 1.0\n", args_info.dopplerMultiplier_arg);
   }
   if (args_info.IHSonly_given) {
      fprintf(stderr,"WARNING: Only IHS stage is being used\n");
   }
   if (args_info.noNotchHarmonics_given) {
      fprintf(stderr,"WARNING: Daily and sidereal period harmonics are not being notched out\n");
   }
   if (args_info.calcRthreshold_given) {
      fprintf(stderr,"WARNING: R threshold values for templates is being calculated with Monte Carlo simulations\n");
   }
   if (args_info.BrentsMethod_given) {
      fprintf(stderr,"WARNING: Using Brent's method for root finding instead of Newton's method.\n");
   }
   if (args_info.antennaOff_given) {
      fprintf(stderr,"WARNING: Antenna pattern weights are all being set to 1.0\n");
   }
   if (args_info.noiseWeightOff_given) {
      fprintf(stderr,"WARNING: Noise weights are all being set to 1.0\n");
   }
   if (args_info.gaussTemplatesOnly_given) {
      fprintf(stderr,"WARNING: Only Gaussian templates will be used\n");
   }
   if (args_info.validateSSE_given) {
      fprintf(stderr,"WARNING: SSE computations will be validated\n");
   }
   if (args_info.ULoff_given) {
      fprintf(stderr,"WARNING: --ULoff has been specifed; no upper limits will be produced\n");
   }
   if (args_info.printSFTtimes_given) {
      fprintf(stderr,"WARNING: input SFT start times are being saved\n");
   }
   if (args_info.printUsedSFTtimes_given) {
      fprintf(stderr,"WARNING: used SFT start times are being saved\n");
   }
   if (args_info.randSeed_given) {
      fprintf(stderr,"NOTE: random seed value %d is being used\n", args_info.randSeed_arg);
   }
   if (args_info.chooseSeed_given) {
      fprintf(stderr,"NOTE: random seed valueis being chosen based on the input search parameters\n");
   }

   //Extra warning that bad SFTs are being marked and removed
   if (args_info.markBadSFTs_given) {
      fprintf(stderr,"WARNING: Marking bad SFTs\n");
   }
   
   //Adjust parameter space search values, if necessary
   if (params->Pmax < params->Pmin) {
      REAL4 tempP = params->Pmax;
      params->Pmax = params->Pmin;
      params->Pmin = tempP;
      fprintf(stderr,"WARNING! Maximum period is smaller than minimum period... switching the two\n");
   }
   if (params->dfmax < params->dfmin) {
      REAL4 tempdf = params->dfmax;
      params->dfmax = params->dfmin;
      params->dfmin = tempdf;
      fprintf(stderr,"WARNING! Maximum modulation depth is smaller than minimum modulation depth... switching the two\n");
   }
   if (params->Pmax > 0.2*(params->Tobs)) {
      params->Pmax = 0.2*(params->Tobs);
      fprintf(stderr,"WARNING! Adjusting input maximum period to 1/5 the observation time!\n");
   }
   if (params->Pmin < 2.0*3600) {
      params->Pmin = 2.0*3600;
      fprintf(stderr,"WARNING! Adjusting input minimum period to 2 hours!\n");
   }
   if (params->dfmax > (0.5*params->Pmax/params->Tcoh/params->Tcoh)) {
      params->dfmax = floor((0.5*params->Pmax/params->Tcoh/params->Tcoh)*(params->Tcoh))/(params->Tcoh);
      fprintf(stderr,"WARNING! Adjusting input maximum modulation depth due to maximum modulation depth allowed\n");
   }
   if (params->dfmin < 0.5/(params->Tcoh)) {
      params->dfmin = 0.5/(params->Tcoh);
      fprintf(stderr,"WARNING! Adjusting input minimum modulation depth to 1/2 a frequency bin!\n");
   }
   
   //Adjustments for improper modulation depth inputs
   params->dfmin = 0.5*round(2.0*params->dfmin*params->Tcoh)/params->Tcoh;
   params->dfmax = 0.5*round(2.0*params->dfmax*params->Tcoh)/params->Tcoh;
   
   //Upper limit settings take span of search values unless specified
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
   
   //SFT type standard or vladimir (Vladimir's SFT generation program has a different normalization factor than standard v2)
   params->sftType = XLALCalloc(strlen(args_info.sftType_arg)+1, sizeof(*(params->sftType)));
   if (params->sftType==NULL) {
      fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", __func__, sizeof(*(params->sftType)));
      XLAL_ERROR(XLAL_ENOMEM);
   }
   sprintf(params->sftType, "%s", args_info.sftType_arg);
   if (strcmp(params->sftType, "standard")==0) {
      fprintf(stderr,"sftType = %s\n", params->sftType);
   } else if (strcmp(params->sftType, "vladimir")==0) {
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
         fprintf(stderr,"IFO = %s\n", IFO);
         params->det[ii] = lalCachedDetectors[LAL_LLO_4K_DETECTOR]; //L1
      } else if (strcmp("H1", IFO)==0) {
         fprintf(stderr,"IFO = %s\n", IFO);
         params->det[ii] = lalCachedDetectors[LAL_LHO_4K_DETECTOR]; //H1
      } else if (strcmp("V1", IFO)==0) {
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
   earth_ephemeris = XLALCalloc(strlen(args_info.ephemDir_arg)+strlen(args_info.ephemYear_arg)+12, sizeof(*earth_ephemeris));
   sun_ephemeris = XLALCalloc(strlen(args_info.ephemDir_arg)+strlen(args_info.ephemYear_arg)+12, sizeof(*sun_ephemeris));
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
   if ((args_info.sftDir_given && args_info.sftFile_given) || (!args_info.sftDir_given && !args_info.sftFile_given)) {
      fprintf(stderr, "%s: One of either sftDir or sftFile must be given but not both or neither.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   } else if (args_info.sftDir_given && !args_info.sftFile_given) {
      sft_dir_file = XLALCalloc(strlen(args_info.sftDir_arg)+20, sizeof(*sft_dir_file));
      if (sft_dir_file==NULL) {
	 fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", __func__, sizeof(*sft_dir_file));
	 XLAL_ERROR(XLAL_ENOMEM);
      }
      sprintf(sft_dir_file, "%s/*.sft", args_info.sftDir_arg);
   } else if (!args_info.sftDir_given && args_info.sftFile_given) {
      sft_dir_file = XLALCalloc(strlen(args_info.sftFile_arg)+2, sizeof(*sft_dir_file));
      if (sft_dir_file==NULL) {
	 fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", __func__, sizeof(*sft_dir_file));
	 XLAL_ERROR(XLAL_ENOMEM);
      }
      sprintf(sft_dir_file, "%s", args_info.sftFile_arg);
   }
   
   return 0;
   
} /* readTwoSepctInputParams() */

INT4 qsort_REAL4_compar(const void *a, const void *b)
{
   const REAL4 *y = a;
   const REAL4 *z = b;

   if ( *y < *z ) return -1;
   if ( *y > *z ) return 1;
   return 0;
} /* qsort_REAL4_compar() */

