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
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_roots.h>

#include <lal/LALStdlib.h>

REAL8 f_diff(double indexval, void *params);

struct solver_params {
   double fvalue;
   double *data_array;
};
double f_diff(double indexval, void *params) {
   struct solver_params *p = (struct solver_params *)params;
   return p->fvalue - p->data_array[(int)round(indexval)*9];
}
 
INT4 main(void) {

   //Turn off gsl error handler
   gsl_set_error_handler_off();

   FILE *H1CANDS, *L1CANDS;
   const char *infile1 = "/Users/evgoet/Documents/MATLAB/pulsar/S6/500-520HzH1Candidates.dat";
   const char *infile2 = "/Users/evgoet/Documents/MATLAB/pulsar/S6/500-520HzL1Candidates.dat";
   const char *outfile1 = "/Users/evgoet/Documents/MATLAB/pulsar/S6/500-520HzCandidates_output1.dat";
   const char *outfile2 = "/Users/evgoet/Documents/MATLAB/pulsar/S6/500-520HzCandidates_output2.dat";
   const char *outfile3 = "/Users/evgoet/Documents/MATLAB/pulsar/S6/500-520HzCandidates_final_H1L1.dat";
   const char *outfile4 = "/Users/evgoet/Documents/MATLAB/pulsar/S6/500-520HzCandidates_output3.dat";

   XLAL_CHECK( (H1CANDS = fopen(infile1,"r")) != NULL, XLAL_EIO, "Can't fopen %s", infile1 );

   //Determines number of candidates in the file
   int ch, h1count = 0, l1count = 0;
   do {
      ch = fgetc(H1CANDS);
      if (ch == '\n') h1count++;
   } while (ch != EOF);

   XLAL_CHECK( (L1CANDS = fopen(infile2,"r")) != NULL, XLAL_EIO, "Can't fopen %s", infile2 );

   //Determines number of candidates in the file
   do {
      ch = fgetc(L1CANDS);
      if (ch == '\n') l1count++;
   } while (ch != EOF);

   double *allh1cands = NULL, *alll1cands = NULL;
   int *allh1cands_job = NULL, *alll1cands_job = NULL;
   XLAL_CHECK( (allh1cands = (double*)XLALMalloc(sizeof(double)*h1count*9)) != NULL, XLAL_ENOMEM );
   XLAL_CHECK( (alll1cands = (double*)XLALMalloc(sizeof(double)*l1count*9)) != NULL, XLAL_ENOMEM );
   XLAL_CHECK( (allh1cands_job = (int*)XLALMalloc(sizeof(int)*h1count)) != NULL, XLAL_ENOMEM );
   XLAL_CHECK( (alll1cands_job = (int*)XLALMalloc(sizeof(int)*l1count)) != NULL, XLAL_ENOMEM );

   //Reset the pointer in the streams
   rewind(H1CANDS);
   rewind(L1CANDS);

   //Put the data into the array
   int ii, jj;
   for (ii=0; ii<h1count; ii++) {
      fscanf(H1CANDS, "%la %la %la %la %la %la %la %la %la %d", &(allh1cands[ii*9]), &(allh1cands[ii*9+1]), &(allh1cands[ii*9+2]), &(allh1cands[ii*9+3]), &(allh1cands[ii*9+4]), &(allh1cands[ii*9+5]), &(allh1cands[ii*9+6]), &(allh1cands[ii*9+7]), &(allh1cands[ii*9+8]), &(allh1cands_job[ii]));
   }

   //Sort the array based on the frequency
   size_t *sorted_index = NULL;
   XLAL_CHECK( (sorted_index = (size_t*)XLALMalloc(sizeof(size_t)*h1count)) != NULL, XLAL_ENOMEM );
   double *allh1cands_sorted = NULL;
   XLAL_CHECK( (allh1cands_sorted = (double*)XLALMalloc(sizeof(double)*h1count*9)) != NULL, XLAL_ENOMEM );
   int *allh1cands_job_sorted = NULL;
   XLAL_CHECK( (allh1cands_job_sorted = (int*)XLALMalloc(sizeof(int)*h1count)) != NULL, XLAL_ENOMEM );
   gsl_sort_index(sorted_index, allh1cands, 9, h1count);
   for (ii=0; ii<h1count; ii++) {
      memcpy(&(allh1cands_sorted[ii*9]), &(allh1cands[sorted_index[ii]*9]), sizeof(double)*9);
      allh1cands_job_sorted[ii] = allh1cands_job[sorted_index[ii]];
   }
   XLALFree(allh1cands);
   XLALFree(sorted_index);
   XLALFree(allh1cands_job);

   //Put the data into the array
   for (ii=0; ii<l1count; ii++) {
      fscanf(L1CANDS, "%la %la %la %la %la %la %la %la %la %d", &(alll1cands[ii*9]), &(alll1cands[ii*9+1]), &(alll1cands[ii*9+2]), &(alll1cands[ii*9+3]), &(alll1cands[ii*9+4]), &(alll1cands[ii*9+5]), &(alll1cands[ii*9+6]), &(alll1cands[ii*9+7]), &(alll1cands[ii*9+8]), &(alll1cands_job[ii]));
   }

   //TODO: Remove this!
   // for (ii=0; ii<l1count; ii++) alll1cands[ii*9] += 1.0;

   //Sort the array based on the frequency
   XLAL_CHECK( (sorted_index = (size_t*)XLALMalloc(sizeof(size_t)*l1count)) != NULL, XLAL_ENOMEM );
   double *alll1cands_sorted = NULL;
   XLAL_CHECK( (alll1cands_sorted = (double*)XLALMalloc(sizeof(double)*l1count*9)) != NULL, XLAL_ENOMEM );
   int *alll1cands_job_sorted = NULL;
   XLAL_CHECK( (alll1cands_job_sorted = (int*)XLALMalloc(sizeof(int)*l1count)) != NULL, XLAL_ENOMEM );
   gsl_sort_index(sorted_index, alll1cands, 9, l1count);
   for (ii=0; ii<l1count; ii++) {
      memcpy(&(alll1cands_sorted[ii*9]), &(alll1cands[sorted_index[ii]*9]), sizeof(double)*9);
      alll1cands_job_sorted[ii] = alll1cands_job[sorted_index[ii]];
   }
   XLALFree(alll1cands);
   XLALFree(sorted_index);
   XLALFree(alll1cands_job);

   //Close the streams
   fclose(H1CANDS);
   fclose(L1CANDS);

   //Open a file to save the output data
   FILE *CANDS = NULL;
   XLAL_CHECK( (CANDS = fopen(outfile1,"w")) != NULL, XLAL_EIO, "Can't fopen %s", outfile1 );

   //Setup and allocate the solver
   int status = GSL_CONTINUE;
   int max_iter = 100;
   double x_lo = 0.0, x_hi = (double)(l1count-1);
   double foundIndex = -1.0;
   gsl_function F;
   F.function = &f_diff;
   const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
   gsl_root_fsolver *s = NULL;
   XLAL_CHECK( (s = gsl_root_fsolver_alloc(T)) != NULL, XLAL_EFUNC );

   double tobs = 40551300.0;
   double fdiff_allowed = 1.0/1800.0;
   double dfdiff_allowed = fdiff_allowed;
   double skydiff_allowed = 0.04*200.0;
   int numpassingf = 0, numpassingdf = 0, numpassingP = 0, numpassingskyloc = 0;
   for (ii=0; ii<h1count; ii++) {

      //Check that the frequency of the H1 candidate we look at is within the frequency span of the L1 candidates
      if ( allh1cands_sorted[ii*9] >= alll1cands_sorted[0] && allh1cands_sorted[ii*9] <= alll1cands_sorted[(l1count-1)*9] ) {
         if (foundIndex < 0.0 || status != GSL_SUCCESS || (status==GSL_SUCCESS && fabs(allh1cands_sorted[ii*9]-alll1cands_sorted[(int)round(gsl_root_fsolver_root(s))*9])>fdiff_allowed)) {
            //Do a root finding search for the closest L1 candidate in frequency to the H1 candidate frequency
            int iter = 0;
            struct solver_params params = {allh1cands_sorted[ii*9], alll1cands_sorted};
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
            while (jj>0 && (allh1cands_sorted[ii*9]-alll1cands_sorted[jj*9])<1.05*fdiff_allowed) jj--;

            //Set the foundIndex value to the start of where we should search from
            foundIndex = jj;

            //Starting from the L1 candidate below the H1 candidate frequency
            for ( ; jj<l1count; jj++) {
               //Check that if the frequency of L1 candidate is above the H1 value by greater than the allowed value, break the loop
               if (alll1cands_sorted[jj*9]-allh1cands_sorted[ii*9] > 1.05*fdiff_allowed) break;

               //If the H1 and L1 frequency values are near enough, proceed with checking more parameter values
               if (fabs(allh1cands_sorted[ii*9]-alll1cands_sorted[jj*9])<=fdiff_allowed) {
                  numpassingf++;
                  //Check the modulation depth
                  if (fabs(allh1cands_sorted[ii*9+2]-alll1cands_sorted[jj*9+2])<=dfdiff_allowed) {
                     numpassingdf++;
                     //Check the period and harmonic values
                     double Pdiff_allowed = allh1cands_sorted[ii*9+1]*allh1cands_sorted[ii*9+1]*sqrt(3.6e-3/allh1cands_sorted[ii*9+2])/(4.5*tobs);
                     double Pdiff_allowed_2 = alll1cands_sorted[jj*9+1]*alll1cands_sorted[jj*9+1]*sqrt(3.6e-3/alll1cands_sorted[jj*9+2])/(4.5*tobs);
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

                           if (fabs(P1factor*allh1cands_sorted[ii*9+1]-P2factor*alll1cands_sorted[jj*9+1])<=Pdiff_allowed*(P1factor*P1factor) || fabs(P1factor*allh1cands_sorted[ii*9+1]-P2factor*alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2*(P2factor*P2factor)) {
                              if (!passedtestP) {
                                 numpassingP++;
                                 passedtestP = 1;
                              }
                              //Check the sky location
                              double absd1mPo2 = fabs(allh1cands_sorted[ii*9+4]-M_PI_2);
                              double absd2mPo2 = fabs(alll1cands_sorted[jj*9+4]-M_PI_2);
                              double dist = acos(sin(absd1mPo2)*sin(absd2mPo2)*cos(allh1cands_sorted[ii*9+3]-alll1cands_sorted[jj*9+3])+cos(absd1mPo2)*cos(absd2mPo2));

                              if (dist<=2.0*skydiff_allowed/(allh1cands_sorted[ii*9]+alll1cands_sorted[jj*9])) {
                                 foundmatch = 1;
                                 numpassingskyloc++;
                                 fprintf(CANDS, "%f %f %f %f %f %f %g %f %g %d %f %f %f %f %f %f %g %f %g %d\n",  (float)allh1cands_sorted[ii*9], (float)allh1cands_sorted[ii*9+1], (float)allh1cands_sorted[ii*9+2], (float)allh1cands_sorted[ii*9+3], (float)allh1cands_sorted[ii*9+4], (float)allh1cands_sorted[ii*9+5], allh1cands_sorted[ii*9+6], (float)allh1cands_sorted[ii*9+7], allh1cands_sorted[ii*9+8], allh1cands_job_sorted[ii], (float)alll1cands_sorted[jj*9], (float)alll1cands_sorted[jj*9+1], (float)alll1cands_sorted[jj*9+2], (float)alll1cands_sorted[jj*9+3], (float)alll1cands_sorted[jj*9+4], (float)alll1cands_sorted[jj*9+5], alll1cands_sorted[jj*9+6], (float)alll1cands_sorted[jj*9+7], alll1cands_sorted[jj*9+8], alll1cands_job_sorted[jj]);
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
      } else if ((alll1cands_sorted[0]-allh1cands_sorted[ii*9]) <= fdiff_allowed && allh1cands_sorted[ii*9] <= alll1cands_sorted[0]) {

         for (jj=0; jj<l1count; jj++) {
            if (alll1cands_sorted[jj*9]-allh1cands_sorted[ii*9] > 1.05*fdiff_allowed) break;

            if (fabs(allh1cands_sorted[ii*9]-alll1cands_sorted[jj*9])<=fdiff_allowed) {
               numpassingf++;
               if (fabs(allh1cands_sorted[ii*9+2]-alll1cands_sorted[jj*9+2])<=dfdiff_allowed) {
                  numpassingdf++;
                  double Pdiff_allowed = allh1cands_sorted[ii*9+1]*allh1cands_sorted[ii*9+1]*sqrt(3.6e-3/allh1cands_sorted[ii*9+2])/(4.5*tobs);
                  double Pdiff_allowed_2 = alll1cands_sorted[jj*9+1]*alll1cands_sorted[jj*9+1]*sqrt(3.6e-3/alll1cands_sorted[jj*9+2])/(4.5*tobs);
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

                        if (fabs(P1factor*allh1cands_sorted[ii*9+1]-P2factor*alll1cands_sorted[jj*9+1])<=Pdiff_allowed*(P1factor*P1factor) || fabs(P1factor*allh1cands_sorted[ii*9+1]-P2factor*alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2*(P2factor*P2factor)) {
                           if (!passedtestP) {
                              numpassingP++;
                              passedtestP = 1;
                           }
                           //Check the sky location
                           double absd1mPo2 = fabs(allh1cands_sorted[ii*9+4]-M_PI_2);
                           double absd2mPo2 = fabs(alll1cands_sorted[jj*9+4]-M_PI_2);
                           double dist = acos(sin(absd1mPo2)*sin(absd2mPo2)*cos(allh1cands_sorted[ii*9+3]-alll1cands_sorted[jj*9+3])+cos(absd1mPo2)*cos(absd2mPo2));

                           if (dist<=2.0*skydiff_allowed/(allh1cands_sorted[ii*9]+alll1cands_sorted[jj*9])) {
                              foundmatch = 1;
                              numpassingskyloc++;
                              fprintf(CANDS, "%f %f %f %f %f %f %g %f %g %d %f %f %f %f %f %f %g %f %g %d\n",  (float)allh1cands_sorted[ii*9], (float)allh1cands_sorted[ii*9+1], (float)allh1cands_sorted[ii*9+2], (float)allh1cands_sorted[ii*9+3], (float)allh1cands_sorted[ii*9+4], (float)allh1cands_sorted[ii*9+5], allh1cands_sorted[ii*9+6], (float)allh1cands_sorted[ii*9+7], allh1cands_sorted[ii*9+8], allh1cands_job_sorted[ii], (float)alll1cands_sorted[jj*9], (float)alll1cands_sorted[jj*9+1], (float)alll1cands_sorted[jj*9+2], (float)alll1cands_sorted[jj*9+3], (float)alll1cands_sorted[jj*9+4], (float)alll1cands_sorted[jj*9+5], alll1cands_sorted[jj*9+6], (float)alll1cands_sorted[jj*9+7], alll1cands_sorted[jj*9+8], alll1cands_job_sorted[jj]);
                           } //end sky check
                        } //end period check
                        if (foundmatch) break;
                     } //end test different P2factors
                     if (foundmatch) break;
                  } //end test different P1factors
               } //end modulation depth check
            } //end frequency check
         } //end test against L1 values
      } else if ((allh1cands_sorted[ii*9]-alll1cands_sorted[(l1count-1)*9]) <= fdiff_allowed && allh1cands_sorted[ii*9] >= alll1cands_sorted[(l1count-1)*9]) {
         jj = l1count-1;
         while (jj>0 && (allh1cands_sorted[ii*9]-alll1cands_sorted[jj*9])<1.05*fdiff_allowed) jj--;

         for ( ; jj<l1count; jj++) {
            if (fabs(allh1cands_sorted[ii*9]-alll1cands_sorted[jj*9])<=fdiff_allowed) {
               numpassingf++;
               if (fabs(allh1cands_sorted[ii*9+2]-alll1cands_sorted[jj*9+2])<=dfdiff_allowed) {
                  numpassingdf++;
                  double Pdiff_allowed = allh1cands_sorted[ii*9+1]*allh1cands_sorted[ii*9+1]*sqrt(3.6e-3/allh1cands_sorted[ii*9+2])/(4.5*tobs);
                  double Pdiff_allowed_2 = alll1cands_sorted[jj*9+1]*alll1cands_sorted[jj*9+1]*sqrt(3.6e-3/alll1cands_sorted[jj*9+2])/(4.5*tobs);
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

                        if (fabs(P1factor*allh1cands_sorted[ii*9+1]-P2factor*alll1cands_sorted[jj*9+1])<=Pdiff_allowed*(P1factor*P1factor) || fabs(P1factor*allh1cands_sorted[ii*9+1]-P2factor*alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2*(P2factor*P2factor)) {
                           if (!passedtestP) {
                              numpassingP++;
                              passedtestP = 1;
                           }
                           //Check the sky location
                           double absd1mPo2 = fabs(allh1cands_sorted[ii*9+4]-M_PI_2);
                           double absd2mPo2 = fabs(alll1cands_sorted[jj*9+4]-M_PI_2);
                           double dist = acos(sin(absd1mPo2)*sin(absd2mPo2)*cos(allh1cands_sorted[ii*9+3]-alll1cands_sorted[jj*9+3])+cos(absd1mPo2)*cos(absd2mPo2));

                           if (dist<=2.0*skydiff_allowed/(allh1cands_sorted[ii*9]+alll1cands_sorted[jj*9])) {
                              foundmatch = 1;
                              numpassingskyloc++;
                              fprintf(CANDS, "%f %f %f %f %f %f %g %f %g %d %f %f %f %f %f %f %g %f %g %d\n",  (float)allh1cands_sorted[ii*9], (float)allh1cands_sorted[ii*9+1], (float)allh1cands_sorted[ii*9+2], (float)allh1cands_sorted[ii*9+3], (float)allh1cands_sorted[ii*9+4], (float)allh1cands_sorted[ii*9+5], allh1cands_sorted[ii*9+6], (float)allh1cands_sorted[ii*9+7], allh1cands_sorted[ii*9+8], allh1cands_job_sorted[ii], (float)alll1cands_sorted[jj*9], (float)alll1cands_sorted[jj*9+1], (float)alll1cands_sorted[jj*9+2], (float)alll1cands_sorted[jj*9+3], (float)alll1cands_sorted[jj*9+4], (float)alll1cands_sorted[jj*9+5], alll1cands_sorted[jj*9+6], (float)alll1cands_sorted[jj*9+7], alll1cands_sorted[jj*9+8], alll1cands_job_sorted[jj]);
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
   XLALFree(allh1cands_sorted);
   XLALFree(alll1cands_sorted);
   XLALFree(allh1cands_job_sorted);
   XLALFree(alll1cands_job_sorted);

   fprintf(stderr, "Passed f = %d, passed df = %d, passed P = %d, passed sky loc %d\n", numpassingf, numpassingdf, numpassingP, numpassingskyloc);


   /// PART TWO: ///

   //open list for reading
   XLAL_CHECK( (CANDS = fopen(outfile1,"r")) != NULL, XLAL_EIO, "Couldn't fopen %s\n", outfile1 );

   //Determines number of candidates in the file
   int count = 0;
   do {
      ch = fgetc(CANDS);
      if (ch == '\n') count++;
   } while (ch != EOF);

   double *allcands = NULL;
   XLAL_CHECK( (allcands = (double*)XLALMalloc(sizeof(double)*count*18)) != NULL, XLAL_ENOMEM );
   int *allcands_job = NULL;
   XLAL_CHECK( (allcands_job = (int*)XLALMalloc(sizeof(int)*count*2)) != NULL, XLAL_ENOMEM );

   //Reset the pointer in the stream
   rewind(CANDS);

   //Put the data into the array
   for (ii=0; ii<count; ii++) {
      fscanf(CANDS, "%la %la %la %la %la %la %la %la %la %d %la %la %la %la %la %la %la %la %la %d", &(allcands[ii*18]), &(allcands[ii*18+1]), &(allcands[ii*18+2]), &(allcands[ii*18+3]), &(allcands[ii*18+4]), &(allcands[ii*18+5]), &(allcands[ii*18+6]), &(allcands[ii*18+7]), &(allcands[ii*18+8]), &(allcands_job[2*ii]), &(allcands[ii*18+9]), &(allcands[ii*18+10]), &(allcands[ii*18+11]), &(allcands[ii*18+12]), &(allcands[ii*18+13]), &(allcands[ii*18+14]), &(allcands[ii*18+15]), &(allcands[ii*18+16]), &(allcands[ii*18+17]), &(allcands_job[ii*2+1]));
   }

   //Close the stream
   fclose(CANDS);
   CANDS = NULL;

   //Open a file to save the output data
   FILE *NEWCANDS = NULL;
   XLAL_CHECK( (NEWCANDS = fopen(outfile2,"w")) != NULL, XLAL_EIO, "Couldn't fopen %s\n", outfile2 );

   for (ii=0; ii<count; ii++) {
      if (allcands[ii*18]!=0.0) {
         int bestcand = ii;
         double bestcandprob = allcands[ii*18+16];
         for (jj=0; jj<count; jj++) {
            if (jj==ii) continue;
            if (allcands[jj*18]!=0.0 && allcands[jj*18] == allcands[ii*18] && allcands[jj*18+1] == allcands[ii*18+1] && allcands[jj*18+2] == allcands[ii*18+2] && allcands[jj*18+3] == allcands[ii*18+3] && allcands[jj*18+4] == allcands[ii*18+4]) {
               if (allcands[jj*18+16]<bestcandprob) {
                  bestcandprob = allcands[jj*18+16];
                  allcands[bestcand*18] = 0.0;
                  bestcand = jj;
               } else {
                  allcands[jj*18] = 0.0;
               }
            }
         }
         fprintf(NEWCANDS, "%f %f %f %f %f %f %g %f %g %d %f %f %f %f %f %f %g %f %g %d\n",  (float)allcands[bestcand*18], (float)allcands[bestcand*18+1], (float)allcands[bestcand*18+2], (float)allcands[bestcand*18+3], (float)allcands[bestcand*18+4], (float)allcands[bestcand*18+5], allcands[bestcand*18+6], (float)allcands[bestcand*18+7], allcands[bestcand*18+8], allcands_job[bestcand*2], (float)allcands[bestcand*18+9], (float)allcands[bestcand*18+10], (float)allcands[bestcand*18+11], (float)allcands[bestcand*18+12], (float)allcands[bestcand*18+13], (float)allcands[bestcand*18+14], allcands[bestcand*18+15], (float)allcands[bestcand*18+16], allcands[bestcand*18+17], allcands_job[bestcand*2+1]);
         allcands[bestcand*18] = 0.0;
      }
   }

   fclose(NEWCANDS);
   NEWCANDS = NULL;

   XLALFree(allcands);
   XLALFree(allcands_job);

   /// PART THREE: ///

   //open list for reading
   XLAL_CHECK( (CANDS = fopen(outfile2,"r")) != NULL, XLAL_EIO, "Couldn't fopen %s\n", outfile2 );

   //Determines number of candidates in the file
   count = 0;
   do {
      ch = fgetc(CANDS);
      if (ch == '\n') count++;
   } while (ch != EOF);

   XLAL_CHECK( (allcands = (double*)XLALMalloc(sizeof(double)*count*18)) != NULL, XLAL_ENOMEM );
   XLAL_CHECK( (allcands_job = (int*)XLALMalloc(sizeof(int)*count*2)) != NULL, XLAL_ENOMEM );

   //Reset the pointer in the stream
   rewind(CANDS);

   //Put the data into the array
   for (ii=0; ii<count; ii++) {
      fscanf(CANDS, "%la %la %la %la %la %la %la %la %la %d %la %la %la %la %la %la %la %la %la %d", &(allcands[ii*18]), &(allcands[ii*18+1]), &(allcands[ii*18+2]), &(allcands[ii*18+3]), &(allcands[ii*18+4]), &(allcands[ii*18+5]), &(allcands[ii*18+6]), &(allcands[ii*18+7]), &(allcands[ii*18+8]), &(allcands_job[2*ii]), &(allcands[ii*18+9]), &(allcands[ii*18+10]), &(allcands[ii*18+11]), &(allcands[ii*18+12]), &(allcands[ii*18+13]), &(allcands[ii*18+14]), &(allcands[ii*18+15]), &(allcands[ii*18+16]), &(allcands[ii*18+17]), &(allcands_job[ii*2+1]));
   }

   //Close the stream
   fclose(CANDS);
   CANDS = NULL;

   //Open a file to save the output data
   XLAL_CHECK( (NEWCANDS = fopen(outfile3,"w")) != NULL, XLAL_EIO, "Couldn't fopen %s\n", outfile3 );

   for (ii=0; ii<count; ii++) {
      if (allcands[ii*18]!=0.0) {
         int bestcand = ii;
         double bestcandprob = allcands[ii*18+7];
         for (jj=0; jj<count; jj++) {
            if (jj==ii) continue;
            if (allcands[jj*18]!=0.0 && allcands[jj*18+9] == allcands[ii*18+9] && allcands[jj*18+10] == allcands[ii*18+10] && allcands[jj*18+11] == allcands[ii*18+11] && allcands[jj*18+12] == allcands[ii*18+12] && allcands[jj*18+13] == allcands[ii*18+13]) {
               if (allcands[jj*18+7]<bestcandprob) {
                  bestcandprob = allcands[jj*18+7];
                  allcands[bestcand*18] = 0.0;
                  bestcand = jj;
               } else {
                  allcands[jj*18] = 0.0;
               }
            }
         }
         fprintf(NEWCANDS, "%f %f %f %f %f %f %g %f %g %d %f %f %f %f %f %f %g %f %g %d\n",  (float)allcands[bestcand*18], (float)allcands[bestcand*18+1], (float)allcands[bestcand*18+2], (float)allcands[bestcand*18+3], (float)allcands[bestcand*18+4], (float)allcands[bestcand*18+5], allcands[bestcand*18+6], (float)allcands[bestcand*18+7], allcands[bestcand*18+8], allcands_job[bestcand*2], (float)allcands[bestcand*18+9], (float)allcands[bestcand*18+10], (float)allcands[bestcand*18+11], (float)allcands[bestcand*18+12], (float)allcands[bestcand*18+13], (float)allcands[bestcand*18+14], allcands[bestcand*18+15], (float)allcands[bestcand*18+16], allcands[bestcand*18+17], allcands_job[bestcand*2+1]);
         allcands[bestcand*18] = 0.0;
      }
   }

   fclose(NEWCANDS);

   XLALFree(allcands);
   XLALFree(allcands_job);


   //PART FOUR: Swap input values
   XLAL_CHECK( (CANDS = fopen(outfile1,"r")) != NULL, XLAL_EIO, "Couldn't fopen %s\n", outfile1 );

   count = 0;
   do {
      ch = fgetc(CANDS);
      if (ch == '\n') count++;
   } while (ch != EOF);

   allcands = NULL;
   XLAL_CHECK( (allcands = (double*)XLALMalloc(sizeof(double)*count*18)) != NULL, XLAL_ENOMEM );
   allcands_job = NULL;
   XLAL_CHECK( (allcands_job = (int*)XLALMalloc(sizeof(int)*count*2)) != NULL, XLAL_ENOMEM );

   //Reset the pointer in the stream
   rewind(CANDS);

   //Put the data into the array
   for (ii=0; ii<count; ii++) {
      fscanf(CANDS, "%la %la %la %la %la %la %la %la %la %d %la %la %la %la %la %la %la %la %la %d", &(allcands[ii*18]), &(allcands[ii*18+1]), &(allcands[ii*18+2]), &(allcands[ii*18+3]), &(allcands[ii*18+4]), &(allcands[ii*18+5]), &(allcands[ii*18+6]), &(allcands[ii*18+7]), &(allcands[ii*18+8]), &(allcands_job[2*ii]), &(allcands[ii*18+9]), &(allcands[ii*18+10]), &(allcands[ii*18+11]), &(allcands[ii*18+12]), &(allcands[ii*18+13]), &(allcands[ii*18+14]), &(allcands[ii*18+15]), &(allcands[ii*18+16]), &(allcands[ii*18+17]), &(allcands_job[ii*2+1]));
   }

   //Close the stream
   fclose(CANDS);
   CANDS = NULL;

   //Open a file to save the output data
   XLAL_CHECK( (NEWCANDS = fopen(outfile4,"w")) != NULL, XLAL_EIO, "Couldn't fopen %s\n", outfile4 );

   for (ii=0; ii<count; ii++) {
      if (allcands[ii*18]!=0.0) {
         int bestcand = ii;
         double bestcandprob = allcands[ii*18+7];
         for (jj=0; jj<count; jj++) {
            if (jj==ii) continue;
            if (allcands[jj*18]!=0.0 && allcands[jj*18+9] == allcands[ii*18+9] && allcands[jj*18+10] == allcands[ii*18+10] && allcands[jj*18+11] == allcands[ii*18+11] && allcands[jj*18+12] == allcands[ii*18+12] && allcands[jj*18+13] == allcands[ii*18+13]) {
               if (allcands[jj*18+7]<bestcandprob) {
                  bestcandprob = allcands[jj*18+7];
                  allcands[bestcand*18] = 0.0;
                  bestcand = jj;
               } else {
                  allcands[jj*18] = 0.0;
               }
            }
         }
         fprintf(NEWCANDS, "%f %f %f %f %f %f %g %f %g %d %f %f %f %f %f %f %g %f %g %d\n",  (float)allcands[bestcand*18], (float)allcands[bestcand*18+1], (float)allcands[bestcand*18+2], (float)allcands[bestcand*18+3], (float)allcands[bestcand*18+4], (float)allcands[bestcand*18+5], allcands[bestcand*18+6], (float)allcands[bestcand*18+7], allcands[bestcand*18+8], allcands_job[bestcand*2], (float)allcands[bestcand*18+9], (float)allcands[bestcand*18+10], (float)allcands[bestcand*18+11], (float)allcands[bestcand*18+12], (float)allcands[bestcand*18+13], (float)allcands[bestcand*18+14], allcands[bestcand*18+15], (float)allcands[bestcand*18+16], allcands[bestcand*18+17], allcands_job[bestcand*2+1]);
         allcands[bestcand*18] = 0.0;
      }
   }

   fclose(NEWCANDS);
   NEWCANDS = NULL;

   XLALFree(allcands);
   XLALFree(allcands_job);

   XLAL_CHECK( (CANDS = fopen(outfile4,"r")) != NULL, XLAL_EIO, "Couldn't fopen %s\n", outfile4 );

   //Determines number of candidates in the file
   count = 0;
   do {
      ch = fgetc(CANDS);
      if (ch == '\n') count++;
   } while (ch != EOF);

   XLAL_CHECK( (allcands = (double*)XLALMalloc(sizeof(double)*count*18)) != NULL, XLAL_ENOMEM );
   XLAL_CHECK( (allcands_job = (int*)XLALMalloc(sizeof(int)*count*2)) != NULL, XLAL_ENOMEM );

   //Reset the pointer in the stream
   rewind(CANDS);

   //Put the data into the array
   for (ii=0; ii<count; ii++) {
      fscanf(CANDS, "%la %la %la %la %la %la %la %la %la %d %la %la %la %la %la %la %la %la %la %d", &(allcands[ii*18]), &(allcands[ii*18+1]), &(allcands[ii*18+2]), &(allcands[ii*18+3]), &(allcands[ii*18+4]), &(allcands[ii*18+5]), &(allcands[ii*18+6]), &(allcands[ii*18+7]), &(allcands[ii*18+8]), &(allcands_job[2*ii]), &(allcands[ii*18+9]), &(allcands[ii*18+10]), &(allcands[ii*18+11]), &(allcands[ii*18+12]), &(allcands[ii*18+13]), &(allcands[ii*18+14]), &(allcands[ii*18+15]), &(allcands[ii*18+16]), &(allcands[ii*18+17]), &(allcands_job[ii*2+1]));
   }

   //Close the stream
   fclose(CANDS);
   CANDS = NULL;

   //Open a file to save the output data
   XLAL_CHECK( (NEWCANDS = fopen(outfile3,"a")) != NULL, XLAL_EIO, "Couldn't fopen %s\n", outfile3 );

   for (ii=0; ii<count; ii++) {
      if (allcands[ii*18]!=0.0) {
         int bestcand = ii;
         double bestcandprob = allcands[ii*18+16];
         for (jj=0; jj<count; jj++) {
            if (jj==ii) continue;
            if (allcands[jj*18]!=0.0 && allcands[jj*18] == allcands[ii*18] && allcands[jj*18+1] == allcands[ii*18+1] && allcands[jj*18+2] == allcands[ii*18+2] && allcands[jj*18+3] == allcands[ii*18+3] && allcands[jj*18+4] == allcands[ii*18+4]) {
               if (allcands[jj*18+16]<bestcandprob) {
                  bestcandprob = allcands[jj*18+16];
                  allcands[bestcand*18] = 0.0;
                  bestcand = jj;
               } else {
                  allcands[jj*18] = 0.0;
               }
            }
         }
         fprintf(NEWCANDS, "%f %f %f %f %f %f %g %f %g %d %f %f %f %f %f %f %g %f %g %d\n",  (float)allcands[bestcand*18], (float)allcands[bestcand*18+1], (float)allcands[bestcand*18+2], (float)allcands[bestcand*18+3], (float)allcands[bestcand*18+4], (float)allcands[bestcand*18+5], allcands[bestcand*18+6], (float)allcands[bestcand*18+7], allcands[bestcand*18+8], allcands_job[bestcand*2], (float)allcands[bestcand*18+9], (float)allcands[bestcand*18+10], (float)allcands[bestcand*18+11], (float)allcands[bestcand*18+12], (float)allcands[bestcand*18+13], (float)allcands[bestcand*18+14], allcands[bestcand*18+15], (float)allcands[bestcand*18+16], allcands[bestcand*18+17], allcands_job[bestcand*2+1]);
         allcands[bestcand*18] = 0.0;
      }
   }

   fclose(NEWCANDS);
   NEWCANDS = NULL;

   XLALFree(allcands);
   XLALFree(allcands_job);

   return 0;
}
