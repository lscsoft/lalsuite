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
#include "../TwoSpect.h"
#include "../statistics.h"


//Global variables
CHAR *earth_ephemeris = NULL, *sun_ephemeris = NULL, *sft_dir_file = NULL;
static const LALStatus empty_status;

//Main program
int main(int argc, char *argv[])
{
   
   INT4 ii;               //counter variables
   
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
   
   //Allocate memory for ffdata structure
   ffdataStruct *ffdata = new_ffdata(inputParams);
   if (ffdata==NULL) {
      fprintf(stderr, "%s: new_ffdata() failed.\n", __func__);
      XLAL_ERROR(XLAL_EFUNC);
   }
   
   //Assume maximum bin shift possible
   inputParams->maxbinshift = (INT4)round(detectorVmax * (inputParams->fmin+0.5*inputParams->fspan) * inputParams->Tcoh)+1;

   //Read in the T-F data from SFTs
   fprintf(stderr, "Loading in SFTs... ");
   ffdata->tfnormalization = 2.0/inputParams->Tcoh/(args_info.avesqrtSh_arg*args_info.avesqrtSh_arg);
   REAL4Vector *tfdata = readInSFTs(inputParams, &(ffdata->tfnormalization));
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
      INT4Vector *removeTheseSFTs = XLALCreateINT4Vector(numffts);
      if (removeTheseSFTs==NULL) {
         fprintf(stderr, "%s: XLALCreateINT4Vector(%d) failed.\n", __func__, numffts);
         XLAL_ERROR(XLAL_EFUNC);
      }
      memset(removeTheseSFTs->data, 0, sizeof(INT4)*removeTheseSFTs->length);
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
            REAL8 kstest = ks_test_exp(tempvect);
            if (XLAL_IS_REAL8_FAIL_NAN(kstest)) {
               fprintf(stderr,"%s: ks_test_exp() failed.\n", __func__);
               XLAL_ERROR(XLAL_EFUNC);
            }

            qsort(tempvect->data, tempvect->length, sizeof(REAL4), qsort_REAL4_compar);
            REAL4 vector_median = 0.0;
            if (tempvect->length % 2 != 1) vector_median = 0.5*(tempvect->data[(INT4)(0.5*tempvect->length)-1] + tempvect->data[(INT4)(0.5*tempvect->length)]);
            else vector_median = tempvect->data[(INT4)(0.5*tempvect->length)];
            REAL4 vector_mean = (REAL4)(vector_median/LAL_LN2);
            REAL8 loval = 0.0, hival = 0.0;
            REAL8 oneoverlength = 1.0/tempvect->length;
            for (ii=0; ii<(INT4)tempvect->length; ii++) {
              REAL8 testval1 = (1.0+ii)*oneoverlength - gsl_cdf_exponential_P(tempvect->data[ii], vector_mean);
              REAL8 testval2 = ii*oneoverlength - gsl_cdf_exponential_P(tempvect->data[ii], vector_mean);
              if (hival<testval1) hival = testval1;
              if (loval<testval2) loval = testval2;
            }
            REAL8 kuiperval1 = hival + loval;
            loval = -1.0, hival = -1.0;
            for (ii=0; ii<(INT4)tempvect->length; ii++) {
              REAL8 testval1 = (1.0+ii)*oneoverlength - gsl_cdf_exponential_P(tempvect->data[ii], vector_mean);
              REAL8 testval2 = ii*oneoverlength - gsl_cdf_exponential_P(tempvect->data[ii], vector_mean);
              if (hival<testval1) hival = testval1;
              if (hival<testval2) hival = testval2;
              if (loval<-testval1) loval = -testval1;
              if (loval<-testval2) loval = -testval2;
            }
            REAL8 kuiperval = hival + loval;

            fprintf(OUTPUT, "%g %g %g\n", kstest, kuiperval1, kuiperval);

            if (kstest>ksthreshold || kuiperval>kuiperthreshold) {
               removeTheseSFTs->data[ii] = 1;
               badsfts++;
            }
         }
      }
      fclose(OUTPUT);
      fprintf(stderr, "Fraction excluded in K-S and Kuiper's tests = %f\n", (REAL4)badsfts/(REAL4)totalsfts);
      XLALDestroyREAL4Vector(tempvect);
      removeBadSFTs(tfdata, removeTheseSFTs);
      fprintf(stderr, "done.\n");
      XLALDestroyINT4Vector(removeTheseSFTs);
   }

   XLALDestroyREAL4Vector(tfdata);
   free_ffdata(ffdata);
   XLALDestroyEphemerisData(edat);
   cmdline_parser_free(&args_info);
   free_inputParams(inputParams);

   return 0;

}
