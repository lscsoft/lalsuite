/*
*  Copyright (C) 2013, 2014 Evan Goetz
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
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_roots.h>

#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>

#include "cmdline_compareCandidates.h"

REAL8 f_diff(double indexval, void *params);
INT4 readInCoincidentOutliers(double **output, int **output_job, int *output_numOutliers, const char *infilename);
INT4 readInIFOoutliersAndSort(double **output, int **output_job, int *output_numCoincident, const char *infilename);

struct solver_params {
   double fvalue;
   double *data_array;
};
double f_diff(double indexval, void *params) {
   struct solver_params *p = (struct solver_params *)params;
   return p->fvalue - p->data_array[(int)round(indexval)*9];
}
INT4 readInCoincidentOutliers(double **output, int **output_job, int *output_numCoincident, const char *infilename) {
   FILE *CANDS = NULL;
   XLAL_CHECK( (CANDS = fopen(infilename,"r")) != NULL, XLAL_EIO, "Can't fopen %s", infilename );
   //Determines number of candidates in the file
   INT4 ch, count = 0;
   do {
      ch = fgetc(CANDS);
      if (ch == '\n') count++;
   } while (ch != EOF);
   double *allcands = NULL;
   XLAL_CHECK( (allcands = (double*)XLALMalloc(sizeof(double)*count*18)) != NULL, XLAL_ENOMEM );
   int *allcands_job = NULL;
   XLAL_CHECK( (allcands_job = (int*)XLALMalloc(sizeof(int)*count*2)) != NULL, XLAL_ENOMEM );
   rewind(CANDS);
   //Put the data into the array
   for (INT4 ii=0; ii<count; ii++) {
      fscanf(CANDS, "%la %la %la %la %la %la %la %la %la %d %la %la %la %la %la %la %la %la %la %d", &(allcands[ii*18]), &(allcands[ii*18+1]), &(allcands[ii*18+2]), &(allcands[ii*18+3]), &(allcands[ii*18+4]), &(allcands[ii*18+5]), &(allcands[ii*18+6]), &(allcands[ii*18+7]), &(allcands[ii*18+8]), &(allcands_job[2*ii]), &(allcands[ii*18+9]), &(allcands[ii*18+10]), &(allcands[ii*18+11]), &(allcands[ii*18+12]), &(allcands[ii*18+13]), &(allcands[ii*18+14]), &(allcands[ii*18+15]), &(allcands[ii*18+16]), &(allcands[ii*18+17]), &(allcands_job[ii*2+1]));
   }
   fclose(CANDS);
   (*output) = allcands;
   (*output_job) = allcands_job;
   (*output_numCoincident) = count;
   return XLAL_SUCCESS;
}
INT4 readInIFOoutliersAndSort(double **output, int **output_job, int *output_numOutliers, const char *infilename) {
   FILE *IFOCANDS = NULL;
   XLAL_CHECK( (IFOCANDS = fopen(infilename,"r")) != NULL, XLAL_EIO, "Can't fopen %s", infilename );
   //Determines number of candidates in the file
   INT4 ch, ifocount = 0;
   do {
      ch = fgetc(IFOCANDS);
      if (ch == '\n') ifocount++;
   } while (ch != EOF);
   double *allIFOcands = NULL;
   XLAL_CHECK( (allIFOcands = (double*)XLALMalloc(sizeof(double)*ifocount*9)) != NULL, XLAL_ENOMEM );
   int *allIFOcands_job = NULL;
   XLAL_CHECK( (allIFOcands_job = (int*)XLALMalloc(sizeof(int)*ifocount)) != NULL, XLAL_ENOMEM );
   rewind(IFOCANDS);
   //Put data in array
   for (INT4 ii=0; ii<ifocount; ii++) {
      fscanf(IFOCANDS, "%la %la %la %la %la %la %la %la %la %d", &(allIFOcands[ii*9]), &(allIFOcands[ii*9+1]), &(allIFOcands[ii*9+2]), &(allIFOcands[ii*9+3]), &(allIFOcands[ii*9+4]), &(allIFOcands[ii*9+5]), &(allIFOcands[ii*9+6]), &(allIFOcands[ii*9+7]), &(allIFOcands[ii*9+8]), &(allIFOcands_job[ii]));
   }
   fclose(IFOCANDS);
   //Sort the array based on the frequency
   size_t *sorted_index = NULL;
   XLAL_CHECK( (sorted_index = (size_t*)XLALMalloc(sizeof(size_t)*ifocount)) != NULL, XLAL_ENOMEM );
   double *allIFOcands_sorted = NULL;
   XLAL_CHECK( (allIFOcands_sorted = (double*)XLALMalloc(sizeof(double)*ifocount*9)) != NULL, XLAL_ENOMEM );
   int *allIFOcands_job_sorted = NULL;
   XLAL_CHECK( (allIFOcands_job_sorted = (int*)XLALMalloc(sizeof(int)*ifocount)) != NULL, XLAL_ENOMEM );
   gsl_sort_index(sorted_index, allIFOcands, 9, ifocount);
   for (INT4 ii=0; ii<ifocount; ii++) {
      memcpy(&(allIFOcands_sorted[ii*9]), &(allIFOcands[sorted_index[ii]*9]), sizeof(double)*9);
      allIFOcands_job_sorted[ii] = allIFOcands_job[sorted_index[ii]];
   }
   XLALFree(allIFOcands);
   XLALFree(sorted_index);
   XLALFree(allIFOcands_job);
   (*output) = allIFOcands_sorted;
   (*output_job) = allIFOcands_job_sorted;
   (*output_numOutliers) = ifocount;
   return XLAL_SUCCESS;
}

INT4 main(int argc, char *argv[]) {

   //Turn off gsl error handler
   gsl_set_error_handler_off();

   struct gengetopt_args_info args_info;
   struct cmdline_parser_params *configparams;
   configparams = cmdline_parser_params_create();
   configparams->check_required = 0;  //don't check for required values at the step
   cmdline_parser_ext(argc, argv, &args_info, configparams);
   configparams->initialize = 0;  //don't reinitialize the parameters structure
   if (args_info.config_given) cmdline_parser_config_file(args_info.config_arg, &args_info, configparams);
   cmdline_parser_required(&args_info, argv[0]);

   double *allIFO1cands_sorted = NULL, *allIFO2cands_sorted = NULL;
   int *allIFO1cands_job_sorted = NULL, *allIFO2cands_job_sorted = NULL;
   int ifo1count = 0, ifo2count = 0;
   XLAL_CHECK( readInIFOoutliersAndSort(&allIFO1cands_sorted, &allIFO1cands_job_sorted, &ifo1count, args_info.infile1_arg) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( readInIFOoutliersAndSort(&allIFO2cands_sorted, &allIFO2cands_job_sorted, &ifo2count, args_info.infile2_arg) == XLAL_SUCCESS, XLAL_EFUNC );

   //Open a file to save the output data
   FILE *CANDS = NULL;
   XLAL_CHECK( (CANDS = fopen(args_info.outfile1_arg,"w")) != NULL, XLAL_EIO, "Can't fopen %s", args_info.outfile1_arg );

   //Setup and allocate the solver
   int status = GSL_CONTINUE;
   int max_iter = 100;
   double x_lo = 0.0, x_hi = (double)(ifo2count-1);
   double foundIndex = -1.0;
   gsl_function F;
   F.function = &f_diff;
   const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
   gsl_root_fsolver *s = NULL;
   XLAL_CHECK( (s = gsl_root_fsolver_alloc(T)) != NULL, XLAL_EFUNC );

   double tobs = args_info.Tobs_arg;
   double fdiff_allowed = args_info.fdiff_allowed_arg;
   double dfdiff_allowed = args_info.dfdiff_allowed_arg;
   double skydiff_allowed = args_info.skydiff_allowed_arg*200.0;
   double periodTcohFactor = 2.7*(args_info.Tcoh_arg/1800.0) + 1.8;
   int numpassingf = 0, numpassingdf = 0, numpassingP = 0, numpassingskyloc = 0;
   INT4 ii, jj;
   for (ii=0; ii<ifo1count; ii++) {

      //Check that the frequency of the H1 candidate we look at is within the frequency span of the L1 candidates
      if ( allIFO1cands_sorted[ii*9] >= allIFO2cands_sorted[0] && allIFO1cands_sorted[ii*9] <= allIFO2cands_sorted[(ifo2count-1)*9] ) {
         if (foundIndex < 0.0 || status != GSL_SUCCESS || (status==GSL_SUCCESS && fabs(allIFO1cands_sorted[ii*9]-allIFO2cands_sorted[(int)round(gsl_root_fsolver_root(s))*9])>fdiff_allowed)) {
            //Do a root finding search for the closest L1 candidate in frequency to the H1 candidate frequency
            int iter = 0;
            struct solver_params params = {allIFO1cands_sorted[ii*9], allIFO2cands_sorted};
            F.params = &params;
            XLAL_CHECK( gsl_root_fsolver_set(s, &F, x_lo, x_hi) == GSL_SUCCESS, XLAL_EFUNC );
            do {
               iter++;
               status = gsl_root_fsolver_iterate(s);
               XLAL_CHECK( status == GSL_SUCCESS, XLAL_EFUNC );
               foundIndex = gsl_root_fsolver_root(s);
               status = gsl_root_test_residual(f_diff(foundIndex, &params), 1.05*fdiff_allowed);
               XLAL_CHECK( status == GSL_SUCCESS || status == GSL_CONTINUE, XLAL_EFUNC );
            } while (status == GSL_CONTINUE && iter < max_iter);
         }

         //If the search was successful, then we step through the L1 candidates to find matching candidates
         if (status == GSL_SUCCESS) {
            jj = (int)round(foundIndex);   //start at the index of the L1 candidate found in the root finding

            //Step backwards in L1 candidates until we are definitely below the H1 candidate in frequency (or at the start of the L1 list)
            while (jj>0 && (allIFO1cands_sorted[ii*9]-allIFO2cands_sorted[jj*9])<1.05*fdiff_allowed) jj--;

            //Set the foundIndex value to the start of where we should search from
            foundIndex = jj;

            //Starting from the L1 candidate below the H1 candidate frequency
            for ( ; jj<ifo2count; jj++) {
               //Check that if the frequency of L1 candidate is above the H1 value by greater than the allowed value, break the loop
               if (allIFO2cands_sorted[jj*9]-allIFO1cands_sorted[ii*9] > 1.05*fdiff_allowed) break;

               //If the H1 and L1 frequency values are near enough, proceed with checking more parameter values
               if (fabs(allIFO1cands_sorted[ii*9]-allIFO2cands_sorted[jj*9])<=fdiff_allowed) {
                  numpassingf++;
                  //Check the modulation depth
                  if (fabs(allIFO1cands_sorted[ii*9+2]-allIFO2cands_sorted[jj*9+2])<=dfdiff_allowed) {
                     numpassingdf++;
                     //Check the period and harmonic values
                     double Pdiff_allowed = allIFO1cands_sorted[ii*9+1]*allIFO1cands_sorted[ii*9+1]*sqrt(3.6e-3/allIFO1cands_sorted[ii*9+2])/(periodTcohFactor*tobs);
                     double Pdiff_allowed_2 = allIFO2cands_sorted[jj*9+1]*allIFO2cands_sorted[jj*9+1]*sqrt(3.6e-3/allIFO2cands_sorted[jj*9+2])/(periodTcohFactor*tobs);
                     int foundmatch = 0, passedtestP = 0;
                     for (int kk=1; kk<=7; kk++) {
                        double P1factor = 0.0;
                        if (kk==1) P1factor = 1.0;
                        else if (kk<5) P1factor = 1.0/kk;
                        else P1factor = (double)(kk-3);

                        for (int ll=1; ll<=7; ll++) {
                           double P2factor = 0.0;
                           if (ll==1) P2factor = 1.0;
                           else if (ll<5) P2factor = 1.0/ll;
                           else P2factor = (double)(ll-3);

                           if (fabs(P1factor*allIFO1cands_sorted[ii*9+1]-P2factor*allIFO2cands_sorted[jj*9+1])<=Pdiff_allowed*(P1factor*P1factor) || fabs(P1factor*allIFO1cands_sorted[ii*9+1]-P2factor*allIFO2cands_sorted[jj*9+1])<=Pdiff_allowed_2*(P2factor*P2factor)) {
                              if (!passedtestP) {
                                 numpassingP++;
                                 passedtestP = 1;
                              }
                              //Check the sky location
                              double absd1mPo2 = fabs(allIFO1cands_sorted[ii*9+4]-M_PI_2);
                              double absd2mPo2 = fabs(allIFO2cands_sorted[jj*9+4]-M_PI_2);
                              double dist = acos(sin(absd1mPo2)*sin(absd2mPo2)*cos(allIFO1cands_sorted[ii*9+3]-allIFO2cands_sorted[jj*9+3])+cos(absd1mPo2)*cos(absd2mPo2));

                              if (dist<=2.0*skydiff_allowed/(allIFO1cands_sorted[ii*9]+allIFO2cands_sorted[jj*9])) {
                                 foundmatch = 1;
                                 numpassingskyloc++;
                                 fprintf(CANDS, "%f %f %f %f %f %f %g %f %g %d %f %f %f %f %f %f %g %f %g %d\n",  allIFO1cands_sorted[ii*9], allIFO1cands_sorted[ii*9+1], allIFO1cands_sorted[ii*9+2], allIFO1cands_sorted[ii*9+3], allIFO1cands_sorted[ii*9+4], allIFO1cands_sorted[ii*9+5], allIFO1cands_sorted[ii*9+6], allIFO1cands_sorted[ii*9+7], allIFO1cands_sorted[ii*9+8], allIFO1cands_job_sorted[ii], allIFO2cands_sorted[jj*9], allIFO2cands_sorted[jj*9+1], allIFO2cands_sorted[jj*9+2], allIFO2cands_sorted[jj*9+3], allIFO2cands_sorted[jj*9+4], allIFO2cands_sorted[jj*9+5], allIFO2cands_sorted[jj*9+6], allIFO2cands_sorted[jj*9+7], allIFO2cands_sorted[jj*9+8], allIFO2cands_job_sorted[jj]);
                              } //end sky check
                           } //end period check
                           if (foundmatch) break;
                        } //end test different P2factors
                        if (foundmatch) break;
                     } //end test different P1factors
                  } //end modulation depth check
               } //end frequency check
            } //end test against L1 values
         } //end successful search
      } else if ((allIFO2cands_sorted[0]-allIFO1cands_sorted[ii*9]) <= fdiff_allowed && allIFO1cands_sorted[ii*9] <= allIFO2cands_sorted[0]) {

         for (jj=0; jj<ifo2count; jj++) {
            if (allIFO2cands_sorted[jj*9]-allIFO1cands_sorted[ii*9] > 1.05*fdiff_allowed) break;

            if (fabs(allIFO1cands_sorted[ii*9]-allIFO2cands_sorted[jj*9])<=fdiff_allowed) {
               numpassingf++;
               if (fabs(allIFO1cands_sorted[ii*9+2]-allIFO2cands_sorted[jj*9+2])<=dfdiff_allowed) {
                  numpassingdf++;
                  double Pdiff_allowed = allIFO1cands_sorted[ii*9+1]*allIFO1cands_sorted[ii*9+1]*sqrt(3.6e-3/allIFO1cands_sorted[ii*9+2])/(periodTcohFactor*tobs);
                  double Pdiff_allowed_2 = allIFO2cands_sorted[jj*9+1]*allIFO2cands_sorted[jj*9+1]*sqrt(3.6e-3/allIFO2cands_sorted[jj*9+2])/(periodTcohFactor*tobs);
                  int foundmatch = 0, passedtestP = 0;
                  for (int kk=1; kk<=7; kk++) {
                     double P1factor = 0.0;
                     if (kk==1) P1factor = 1.0;
                     else if (kk<5) P1factor = 1.0/kk;
                     else P1factor = (double)(kk-3);

                     for (int ll=1; ll<=7; ll++) {
                        double P2factor = 0.0;
                        if (ll==1) P2factor = 1.0;
                        else if (ll<5) P2factor = 1.0/ll;
                        else P2factor = (double)(ll-3);

                        if (fabs(P1factor*allIFO1cands_sorted[ii*9+1]-P2factor*allIFO2cands_sorted[jj*9+1])<=Pdiff_allowed*(P1factor*P1factor) || fabs(P1factor*allIFO1cands_sorted[ii*9+1]-P2factor*allIFO2cands_sorted[jj*9+1])<=Pdiff_allowed_2*(P2factor*P2factor)) {
                           if (!passedtestP) {
                              numpassingP++;
                              passedtestP = 1;
                           }
                           //Check the sky location
                           double absd1mPo2 = fabs(allIFO1cands_sorted[ii*9+4]-M_PI_2);
                           double absd2mPo2 = fabs(allIFO2cands_sorted[jj*9+4]-M_PI_2);
                           double dist = acos(sin(absd1mPo2)*sin(absd2mPo2)*cos(allIFO1cands_sorted[ii*9+3]-allIFO2cands_sorted[jj*9+3])+cos(absd1mPo2)*cos(absd2mPo2));

                           if (dist<=2.0*skydiff_allowed/(allIFO1cands_sorted[ii*9]+allIFO2cands_sorted[jj*9])) {
                              foundmatch = 1;
                              numpassingskyloc++;
                              fprintf(CANDS, "%f %f %f %f %f %f %g %f %g %d %f %f %f %f %f %f %g %f %g %d\n",  allIFO1cands_sorted[ii*9], allIFO1cands_sorted[ii*9+1], allIFO1cands_sorted[ii*9+2], allIFO1cands_sorted[ii*9+3], allIFO1cands_sorted[ii*9+4], allIFO1cands_sorted[ii*9+5], allIFO1cands_sorted[ii*9+6], allIFO1cands_sorted[ii*9+7], allIFO1cands_sorted[ii*9+8], allIFO1cands_job_sorted[ii], allIFO2cands_sorted[jj*9], allIFO2cands_sorted[jj*9+1], allIFO2cands_sorted[jj*9+2], allIFO2cands_sorted[jj*9+3], allIFO2cands_sorted[jj*9+4], allIFO2cands_sorted[jj*9+5], allIFO2cands_sorted[jj*9+6], allIFO2cands_sorted[jj*9+7], allIFO2cands_sorted[jj*9+8], allIFO2cands_job_sorted[jj]);
                           } //end sky check
                        } //end period check
                        if (foundmatch) break;
                     } //end test different P2factors
                     if (foundmatch) break;
                  } //end test different P1factors
               } //end modulation depth check
            } //end frequency check
         } //end test against L1 values
      } else if ((allIFO1cands_sorted[ii*9]-allIFO2cands_sorted[(ifo2count-1)*9]) <= fdiff_allowed && allIFO1cands_sorted[ii*9] >= allIFO2cands_sorted[(ifo2count-1)*9]) {
         jj = ifo2count-1;
         while (jj>0 && (allIFO1cands_sorted[ii*9]-allIFO2cands_sorted[jj*9])<1.05*fdiff_allowed) jj--;

         for ( ; jj<ifo2count; jj++) {
            if (fabs(allIFO1cands_sorted[ii*9]-allIFO2cands_sorted[jj*9])<=fdiff_allowed) {
               numpassingf++;
               if (fabs(allIFO1cands_sorted[ii*9+2]-allIFO2cands_sorted[jj*9+2])<=dfdiff_allowed) {
                  numpassingdf++;
                  double Pdiff_allowed = allIFO1cands_sorted[ii*9+1]*allIFO1cands_sorted[ii*9+1]*sqrt(3.6e-3/allIFO1cands_sorted[ii*9+2])/(periodTcohFactor*tobs);
                  double Pdiff_allowed_2 = allIFO2cands_sorted[jj*9+1]*allIFO2cands_sorted[jj*9+1]*sqrt(3.6e-3/allIFO2cands_sorted[jj*9+2])/(periodTcohFactor*tobs);
                  int foundmatch = 0, passedtestP = 0;
                  for (int kk=1; kk<=7; kk++) {
                     double P1factor = 0.0;
                     if (kk==1) P1factor = 1.0;
                     else if (kk<5) P1factor = 1.0/kk;
                     else P1factor = (double)(kk-3);

                     for (int ll=1; ll<=7; ll++) {
                        double P2factor = 0.0;
                        if (ll==1) P2factor = 1.0;
                        else if (ll<5) P2factor = 1.0/ll;
                        else P2factor = (double)(ll-3);

                        if (fabs(P1factor*allIFO1cands_sorted[ii*9+1]-P2factor*allIFO2cands_sorted[jj*9+1])<=Pdiff_allowed*(P1factor*P1factor) || fabs(P1factor*allIFO1cands_sorted[ii*9+1]-P2factor*allIFO2cands_sorted[jj*9+1])<=Pdiff_allowed_2*(P2factor*P2factor)) {
                           if (!passedtestP) {
                              numpassingP++;
                              passedtestP = 1;
                           }
                           //Check the sky location
                           double absd1mPo2 = fabs(allIFO1cands_sorted[ii*9+4]-M_PI_2);
                           double absd2mPo2 = fabs(allIFO2cands_sorted[jj*9+4]-M_PI_2);
                           double dist = acos(sin(absd1mPo2)*sin(absd2mPo2)*cos(allIFO1cands_sorted[ii*9+3]-allIFO2cands_sorted[jj*9+3])+cos(absd1mPo2)*cos(absd2mPo2));

                           if (dist<=2.0*skydiff_allowed/(allIFO1cands_sorted[ii*9]+allIFO2cands_sorted[jj*9])) {
                              foundmatch = 1;
                              numpassingskyloc++;
                              fprintf(CANDS, "%f %f %f %f %f %f %g %f %g %d %f %f %f %f %f %f %g %f %g %d\n",  allIFO1cands_sorted[ii*9], allIFO1cands_sorted[ii*9+1], allIFO1cands_sorted[ii*9+2], allIFO1cands_sorted[ii*9+3], allIFO1cands_sorted[ii*9+4], allIFO1cands_sorted[ii*9+5], allIFO1cands_sorted[ii*9+6], allIFO1cands_sorted[ii*9+7], allIFO1cands_sorted[ii*9+8], allIFO1cands_job_sorted[ii], allIFO2cands_sorted[jj*9], allIFO2cands_sorted[jj*9+1], allIFO2cands_sorted[jj*9+2], allIFO2cands_sorted[jj*9+3], allIFO2cands_sorted[jj*9+4], allIFO2cands_sorted[jj*9+5], allIFO2cands_sorted[jj*9+6], allIFO2cands_sorted[jj*9+7], allIFO2cands_sorted[jj*9+8], allIFO2cands_job_sorted[jj]);
                           } //end sky check
                        } //end period check
                        if (foundmatch) break;
                     } //end test different P2factors
                     if (foundmatch) break;
                  } //end test different P1factors
               } //end modulation depth check
            } //end frequency check
         } //end check against L1 values
      } //end if H1 candidate is barely outside L1 frequency range
   } //end loop over H1 values

   //Close combined list file
   fclose(CANDS);
   CANDS = NULL;

   gsl_root_fsolver_free(s);
   XLALFree(allIFO1cands_sorted);
   XLALFree(allIFO2cands_sorted);
   XLALFree(allIFO1cands_job_sorted);
   XLALFree(allIFO2cands_job_sorted);

   fprintf(stderr, "Passed f = %d, passed df = %d, passed P = %d, passed sky loc %d\n", numpassingf, numpassingdf, numpassingP, numpassingskyloc);


   /// PART TWO: ///

   double *allcands = NULL;
   int *allcands_job = NULL;
   int count = 0;
   XLAL_CHECK( readInCoincidentOutliers(&allcands, &allcands_job, &count, args_info.outfile1_arg) == XLAL_SUCCESS, XLAL_EFUNC );

   //Allocate usedvalue array
   INT4Vector *usedvalue = NULL;
   XLAL_CHECK( (usedvalue = XLALCreateINT4Vector(count)) != NULL, XLAL_EFUNC );
   memset(usedvalue->data, 0, sizeof(INT4)*count);

   //Open a file to save the output data
   FILE *NEWCANDS = NULL;
   XLAL_CHECK( (NEWCANDS = fopen(args_info.outfile2_arg,"w")) != NULL, XLAL_EIO, "Couldn't fopen %s\n", args_info.outfile2_arg );

   for (ii=0; ii<count; ii++) {
      if (!usedvalue->data[ii]) {
         int bestcand = ii;
         double bestcandprob = allcands[ii*18+16];
         for (jj=0; jj<count; jj++) {
            if (usedvalue->data[jj] || jj==ii) continue;
            if (allcands[jj*18] == allcands[bestcand*18] && allcands[jj*18+1] == allcands[bestcand*18+1] && allcands[jj*18+2] == allcands[bestcand*18+2] && allcands[jj*18+3] == allcands[bestcand*18+3] && allcands[jj*18+4] == allcands[bestcand*18+4]) {
               if (allcands[jj*18+16]<bestcandprob) {
                  usedvalue->data[bestcand] = 1;
                  bestcandprob = allcands[jj*18+16];
                  bestcand = jj;
               } else {
                  usedvalue->data[jj] = 1;
               }
            }
         }
         fprintf(NEWCANDS, "%f %f %f %f %f %f %g %f %g %d %f %f %f %f %f %f %g %f %g %d\n", allcands[bestcand*18], allcands[bestcand*18+1], allcands[bestcand*18+2], allcands[bestcand*18+3], allcands[bestcand*18+4], allcands[bestcand*18+5], allcands[bestcand*18+6], allcands[bestcand*18+7], allcands[bestcand*18+8], allcands_job[bestcand*2], allcands[bestcand*18+9], allcands[bestcand*18+10], allcands[bestcand*18+11], allcands[bestcand*18+12], allcands[bestcand*18+13], allcands[bestcand*18+14], allcands[bestcand*18+15], allcands[bestcand*18+16], allcands[bestcand*18+17], allcands_job[bestcand*2+1]);
         usedvalue->data[bestcand] = 1;
      }
   }

   fclose(NEWCANDS);
   NEWCANDS = NULL;

   XLALFree(allcands);
   XLALFree(allcands_job);
   XLALDestroyINT4Vector(usedvalue);

   /// PART THREE: ///

   XLAL_CHECK( readInCoincidentOutliers(&allcands, &allcands_job, &count, args_info.outfile2_arg) == XLAL_SUCCESS, XLAL_EFUNC );

   XLAL_CHECK( (usedvalue = XLALCreateINT4Vector(count)) != NULL, XLAL_EFUNC );
   memset(usedvalue->data, 0, sizeof(INT4)*count);

   //Open a file to save the output data
   XLAL_CHECK( (NEWCANDS = fopen(args_info.finalOutfile_arg,"w")) != NULL, XLAL_EIO, "Couldn't fopen %s\n", args_info.finalOutfile_arg );

   for (ii=0; ii<count; ii++) {
      if (!usedvalue->data[ii]) {
         int bestcand = ii;
         double bestcandprob = allcands[ii*18+7];
         for (jj=0; jj<count; jj++) {
            if (usedvalue->data[jj] || jj==ii) continue;
            if (allcands[jj*18+9] == allcands[bestcand*18+9] && allcands[jj*18+10] == allcands[bestcand*18+10] && allcands[jj*18+11] == allcands[bestcand*18+11] && allcands[jj*18+12] == allcands[bestcand*18+12] && allcands[jj*18+13] == allcands[bestcand*18+13]) {
               if (allcands[jj*18+7]<bestcandprob) {
                  usedvalue->data[bestcand] = 1;
                  bestcandprob = allcands[jj*18+7];
                  bestcand = jj;
               } else {
                  usedvalue->data[jj] = 1;
               }
            }
         }
         fprintf(NEWCANDS, "%f %f %f %f %f %f %g %f %g %d %f %f %f %f %f %f %g %f %g %d\n", allcands[bestcand*18], allcands[bestcand*18+1], allcands[bestcand*18+2], allcands[bestcand*18+3], allcands[bestcand*18+4], allcands[bestcand*18+5], allcands[bestcand*18+6], allcands[bestcand*18+7], allcands[bestcand*18+8], allcands_job[bestcand*2], allcands[bestcand*18+9], allcands[bestcand*18+10], allcands[bestcand*18+11], allcands[bestcand*18+12], allcands[bestcand*18+13], allcands[bestcand*18+14], allcands[bestcand*18+15], allcands[bestcand*18+16], allcands[bestcand*18+17], allcands_job[bestcand*2+1]);
         usedvalue->data[bestcand] = 1;
      }
   }

   fclose(NEWCANDS);

   XLALFree(allcands);
   XLALFree(allcands_job);
   XLALDestroyINT4Vector(usedvalue);


   //PART FOUR: Swap input values
   XLAL_CHECK( readInCoincidentOutliers(&allcands, &allcands_job, &count, args_info.outfile1_arg) == XLAL_SUCCESS, XLAL_EFUNC );

   XLAL_CHECK( (usedvalue = XLALCreateINT4Vector(count)) != NULL, XLAL_EFUNC );
   memset(usedvalue->data, 0, sizeof(INT4)*count);

   //Open a file to save the output data
   XLAL_CHECK( (NEWCANDS = fopen(args_info.outfile3_arg,"w")) != NULL, XLAL_EIO, "Couldn't fopen %s\n", args_info.outfile3_arg );

   for (ii=0; ii<count; ii++) {
      if (!usedvalue->data[ii]) {
         int bestcand = ii;
         double bestcandprob = allcands[ii*18+7];
         for (jj=0; jj<count; jj++) {
            if (usedvalue->data[jj] || jj==ii) continue;
            if (allcands[jj*18+9] == allcands[bestcand*18+9] && allcands[jj*18+10] == allcands[bestcand*18+10] && allcands[jj*18+11] == allcands[bestcand*18+11] && allcands[jj*18+12] == allcands[bestcand*18+12] && allcands[jj*18+13] == allcands[bestcand*18+13]) {
               if (allcands[jj*18+7]<bestcandprob) {
                  usedvalue->data[bestcand] = 1;
                  bestcandprob = allcands[jj*18+7];
                  bestcand = jj;
               } else {
                   usedvalue->data[jj] = 1;
               }
            }
         }
         fprintf(NEWCANDS, "%f %f %f %f %f %f %g %f %g %d %f %f %f %f %f %f %g %f %g %d\n", allcands[bestcand*18], allcands[bestcand*18+1], allcands[bestcand*18+2], allcands[bestcand*18+3], allcands[bestcand*18+4], allcands[bestcand*18+5], allcands[bestcand*18+6], allcands[bestcand*18+7], allcands[bestcand*18+8], allcands_job[bestcand*2], allcands[bestcand*18+9], allcands[bestcand*18+10], allcands[bestcand*18+11], allcands[bestcand*18+12], allcands[bestcand*18+13], allcands[bestcand*18+14], allcands[bestcand*18+15], allcands[bestcand*18+16], allcands[bestcand*18+17], allcands_job[bestcand*2+1]);
          usedvalue->data[bestcand] = 1;
      }
   }

   fclose(NEWCANDS);
   NEWCANDS = NULL;

   XLALFree(allcands);
   XLALFree(allcands_job);
   XLALDestroyINT4Vector(usedvalue);


   XLAL_CHECK( readInCoincidentOutliers(&allcands, &allcands_job, &count, args_info.outfile3_arg) == XLAL_SUCCESS, XLAL_EFUNC );

   XLAL_CHECK( (usedvalue = XLALCreateINT4Vector(count)) != NULL, XLAL_EFUNC );
   memset(usedvalue->data, 0, sizeof(INT4)*count);

   //Open a file to save the output data
   XLAL_CHECK( (NEWCANDS = fopen(args_info.finalOutfile_arg,"a")) != NULL, XLAL_EIO, "Couldn't fopen %s\n", args_info.finalOutfile_arg );

   for (ii=0; ii<count; ii++) {
      if (!usedvalue->data[ii]) {
         int bestcand = ii;
         double bestcandprob = allcands[ii*18+16];
         for (jj=0; jj<count; jj++) {
            if (usedvalue->data[jj] || jj==ii) continue;
            if (allcands[jj*18] == allcands[bestcand*18] && allcands[jj*18+1] == allcands[bestcand*18+1] && allcands[jj*18+2] == allcands[bestcand*18+2] && allcands[jj*18+3] == allcands[bestcand*18+3] && allcands[jj*18+4] == allcands[bestcand*18+4]) {
               if (allcands[jj*18+16]<bestcandprob) {
                  usedvalue->data[bestcand] = 1;
                  bestcandprob = allcands[jj*18+16];
                  bestcand = jj;
               } else {
                  usedvalue->data[jj] = 1;
               }
            }
         }
         fprintf(NEWCANDS, "%f %f %f %f %f %f %g %f %g %d %f %f %f %f %f %f %g %f %g %d\n", allcands[bestcand*18], allcands[bestcand*18+1], allcands[bestcand*18+2], allcands[bestcand*18+3], allcands[bestcand*18+4], allcands[bestcand*18+5], allcands[bestcand*18+6], allcands[bestcand*18+7], allcands[bestcand*18+8], allcands_job[bestcand*2], allcands[bestcand*18+9], allcands[bestcand*18+10], allcands[bestcand*18+11], allcands[bestcand*18+12], allcands[bestcand*18+13], allcands[bestcand*18+14], allcands[bestcand*18+15], allcands[bestcand*18+16], allcands[bestcand*18+17], allcands_job[bestcand*2+1]);
         usedvalue->data[bestcand] = 1;
      }
   }

   fclose(NEWCANDS);
   NEWCANDS = NULL;

   XLALFree(allcands);
   XLALFree(allcands_job);
   XLALDestroyINT4Vector(usedvalue);

   return 0;
}
