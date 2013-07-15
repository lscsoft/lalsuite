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
#include <string.h>
#include <math.h>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_roots.h>

struct solver_params {
   double fvalue;
   double *data_array;
};
double fdiff(double index, void *params) {
   struct solver_params *p = (struct solver_params *)params;
   return p->fvalue - p->data_array[(int)round(index)*9];
}
 
int main(void) {
   FILE *H1CANDS, *L1CANDS;
   H1CANDS = fopen("/Users/evgoet/Documents/MATLAB/pulsar/S6/50-350HzH1candidates.dat","r");
   if (H1CANDS == NULL) {
      fprintf(stderr, "%s: %s does not exist\n", __func__, "/Users/evgoet/Documents/MATLAB/pulsar/S6/50-350HzH1candidates.dat");
      exit(1);
   }

   //Determines number of candidates in the file
   int ch, h1count = 0, l1count = 0;
   do {
      ch = fgetc(H1CANDS);
      if (ch == '\n') h1count++;
   } while (ch != EOF);

   L1CANDS = fopen("/Users/evgoet/Documents/MATLAB/pulsar/S6/50-350HzL1candidates.dat","r");
   if (L1CANDS == NULL) {
      fprintf(stderr, "%s: %s does not exist\n", __func__, "/Users/evgoet/Documents/MATLAB/pulsar/S6/50-350HzL1candidates.dat");
      exit(1);
   }

   //Determines number of candidates in the file
   do {
      ch = fgetc(L1CANDS);
      if (ch == '\n') l1count++;
   } while (ch != EOF);

   double *allh1cands = (double*)malloc(sizeof(double)*h1count*9);
   double *alll1cands = (double*)malloc(sizeof(double)*l1count*9);
   int *allh1cands_job = (int*)malloc(sizeof(int)*h1count);
   int *alll1cands_job = (int*)malloc(sizeof(int)*l1count);

   //Reset the pointer in the streams
   rewind(H1CANDS);
   rewind(L1CANDS);

   //Put the data into the array
   int ii, jj;
   for (ii=0; ii<h1count; ii++) {
      fscanf(H1CANDS, "%la %la %la %la %la %la %la %la %la %d", &(allh1cands[ii*9]), &(allh1cands[ii*9+1]), &(allh1cands[ii*9+2]), &(allh1cands[ii*9+3]), &(allh1cands[ii*9+4]), &(allh1cands[ii*9+5]), &(allh1cands[ii*9+6]), &(allh1cands[ii*9+7]), &(allh1cands[ii*9+8]), &(allh1cands_job[ii]));
   }

   //Sort the array based on the frequency
   size_t *sorted_index = (size_t*)malloc(sizeof(size_t)*h1count);
   double *allh1cands_sorted = (double*)malloc(sizeof(double)*h1count*9);
   int *allh1cands_job_sorted = (int*)malloc(sizeof(int)*h1count);
   gsl_sort_index(sorted_index, allh1cands, 9, h1count);
   for (ii=0; ii<h1count; ii++) {
      memcpy(&(allh1cands_sorted[ii*9]), &(allh1cands[sorted_index[ii]*9]), sizeof(double)*9);
      allh1cands_job_sorted[ii] = allh1cands_job[sorted_index[ii]];
   }
   free(allh1cands);
   free(sorted_index);
   free(allh1cands_job);

   //Put the data into the array
   for (ii=0; ii<l1count; ii++) {
      fscanf(L1CANDS, "%la %la %la %la %la %la %la %la %la %d", &(alll1cands[ii*9]), &(alll1cands[ii*9+1]), &(alll1cands[ii*9+2]), &(alll1cands[ii*9+3]), &(alll1cands[ii*9+4]), &(alll1cands[ii*9+5]), &(alll1cands[ii*9+6]), &(alll1cands[ii*9+7]), &(alll1cands[ii*9+8]), &(alll1cands_job[ii]));
   }

   //Sort the array based on the frequency
   sorted_index = (size_t*)malloc(sizeof(size_t)*l1count);
   double *alll1cands_sorted = (double*)malloc(sizeof(double)*l1count*9);
   int *alll1cands_job_sorted = (int*)malloc(sizeof(int)*l1count);
   gsl_sort_index(sorted_index, alll1cands, 9, l1count);
   for (ii=0; ii<l1count; ii++) {
      memcpy(&(alll1cands_sorted[ii*9]), &(alll1cands[sorted_index[ii]*9]), sizeof(double)*9);
      alll1cands_job_sorted[ii] = alll1cands_job[sorted_index[ii]];
   }
   free(alll1cands);
   free(sorted_index);
   free(alll1cands_job);

   //Close the streams
   fclose(H1CANDS);
   fclose(L1CANDS);

   //Open a file to save the output data
   FILE *CANDS = fopen("/Users/evgoet/Documents/MATLAB/pulsar/S6/50-350HzCandidates.dat","w");

   //Setup and allocate the solver
   int status;
   int max_iter = 100;
   double x_lo = 0.0, x_hi = (double)(l1count-1);
   double foundIndex = 0.0;
   gsl_function F;
   F.function = &fdiff;
   const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
   gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

   double tobs = 40551300.0;
   double fdiff_allowed = 1.0/1800.0;
   double dfdiff_allowed = fdiff_allowed;
   double skydiff_allowed = 0.04*200.0;
   for (ii=0; ii<h1count; ii++) {

      //Check that the frequency of the H1 candidate we look at is within the frequency span of the L1 candidates
      if ( (allh1cands_sorted[ii*9] - alll1cands_sorted[0]) >= 0.0 && (allh1cands_sorted[ii*9] - alll1cands_sorted[(l1count-1)*9]) <= 0.0) {
         //Do a root finding search for the closest L1 candidate in frequency to the H1 candidate frequency
         int iter = 0;
         struct solver_params params = {allh1cands_sorted[ii*9], alll1cands_sorted};
         F.params = &params;
         gsl_root_fsolver_set (s, &F, x_lo, x_hi);
         do {
            iter++;
            status = gsl_root_fsolver_iterate(s);
            foundIndex = gsl_root_fsolver_root(s);
            status = gsl_root_test_residual(fdiff(foundIndex, &params), 1.05*fdiff_allowed);
         } while (status == GSL_CONTINUE && iter < max_iter);

         //If the search was successful, then we step through the L1 candidates to find matching candidates
         if (status == GSL_SUCCESS && iter < max_iter) {
            jj = (int)round(foundIndex);   //start at the index of the L1 candidate found in the root finding

            int bestmatch = -1;
            double bestmatchprob = 0.0;

            //Step backwards in L1 candidates until we are definitely below the H1 candidate in frequency (or at the start of the L1 list)
            while (jj>0 && (allh1cands_sorted[ii*9]-alll1cands_sorted[jj*9])<1.05*fdiff_allowed) jj--;

            //Starting from the L1 candidate below the H1 candidate frequency
            for (/* jj value */; jj<l1count; jj++) {
               //Check that if the frequency of L1 candidate is above the H1 value by greater than the allowed value, break the loop
               if (allh1cands_sorted[ii*9]-alll1cands_sorted[jj*9]<-1.05*fdiff_allowed) break;

               //If the H1 and L1 frequency values are near enough, proceed with checking more parameter values
               if (fabs(allh1cands_sorted[ii*9]-alll1cands_sorted[jj*9])<=fdiff_allowed) {
                  //Check the modulation depth
                  if (fabs(allh1cands_sorted[ii*9+2]-alll1cands_sorted[jj*9+2])<=dfdiff_allowed) {
                     //Check the period and harmonic values
                     double Pdiff_allowed = 1.5*allh1cands_sorted[ii*9+1]*allh1cands_sorted[ii*9+1]*sqrt(3.6e-3/allh1cands_sorted[ii*9+2])/(4.5*tobs);
                     double Pdiff_allowed_2 = 1.5*alll1cands_sorted[jj*9+1]*alll1cands_sorted[jj*9+1]*sqrt(3.6e-3/alll1cands_sorted[jj*9+2])/(4.5*tobs);
                     if (fabs(allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed || fabs(allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2 || fabs(allh1cands_sorted[ii*9+1]-2.0*alll1cands_sorted[jj*9+1])<=Pdiff_allowed || fabs(2.0*allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2 || fabs(allh1cands_sorted[ii*9+1]-3.0*alll1cands_sorted[jj*9+1])<=Pdiff_allowed || fabs(3.0*allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2 || fabs(allh1cands_sorted[ii*9+1]-4.0*alll1cands_sorted[jj*9+1])<=Pdiff_allowed || fabs(4.0*allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2 || fabs(allh1cands_sorted[ii*9+1]-0.5*alll1cands_sorted[jj*9+1])<=Pdiff_allowed || fabs(0.5*allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2 || fabs(allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1]/3.0)<=Pdiff_allowed || fabs(allh1cands_sorted[ii*9+1]/3.0-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2 || fabs(allh1cands_sorted[ii*9+1]-0.25*alll1cands_sorted[jj*9+1])<=Pdiff_allowed || fabs(0.25*allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2) {
                        //Check the sky location
                        double absd1mPo2 = fabs(allh1cands_sorted[ii*9+4]-M_PI_2);
                        double absd2mPo2 = fabs(alll1cands_sorted[jj*9+4]-M_PI_2);
                        double dist = acos(sin(absd1mPo2)*sin(absd2mPo2)*cos(allh1cands_sorted[ii*9+3]-alll1cands_sorted[jj*9+3])+cos(absd1mPo2)*cos(absd2mPo2));

                        if (dist<=2.0*skydiff_allowed/(allh1cands_sorted[ii*9]+alll1cands_sorted[jj*9])) {
                           if (bestmatchprob>alll1cands_sorted[jj*9+7]) {
                              bestmatchprob = alll1cands_sorted[jj*9+7];
                              bestmatch = jj;
                           }
                        } //end sky check
                     } //end period check
                  } //end modulation depth check
               } //end frequency check
            } //end test against L1 values
            if (bestmatch>=0) {
               fprintf(CANDS, "%f %f %f %f %f %f %g %f %g %d %f %f %f %f %f %f %g %f %g %d\n",  (float)allh1cands_sorted[ii*9], (float)allh1cands_sorted[ii*9+1], (float)allh1cands_sorted[ii*9+2], (float)allh1cands_sorted[ii*9+3], (float)allh1cands_sorted[ii*9+4], (float)allh1cands_sorted[ii*9+5], allh1cands_sorted[ii*9+6], (float)allh1cands_sorted[ii*9+7], allh1cands_sorted[ii*9+8], allh1cands_job_sorted[ii], (float)alll1cands_sorted[bestmatch*9], (float)alll1cands_sorted[bestmatch*9+1], (float)alll1cands_sorted[bestmatch*9+2], (float)alll1cands_sorted[bestmatch*9+3], (float)alll1cands_sorted[bestmatch*9+4], (float)alll1cands_sorted[bestmatch*9+5], alll1cands_sorted[bestmatch*9+6], (float)alll1cands_sorted[bestmatch*9+7], alll1cands_sorted[bestmatch*9+8], alll1cands_job_sorted[bestmatch]);
            }
         } //end successful search
      } else if ((alll1cands_sorted[0]-allh1cands_sorted[ii*9])<=fdiff_allowed && (alll1cands_sorted[0]-allh1cands_sorted[ii*9])>=0.0) {
         int bestmatch = -1;
         double bestmatchprob = 0.0;

         for (jj=0; jj<l1count; jj++) {
            if (allh1cands_sorted[ii*9]-alll1cands_sorted[jj*9]<-1.05*fdiff_allowed) break;

            if (fabs(allh1cands_sorted[ii*9]-alll1cands_sorted[jj*9])<=fdiff_allowed) {
               if (fabs(allh1cands_sorted[ii*9+2]-alll1cands_sorted[jj*9+2])<=dfdiff_allowed) {
                  double Pdiff_allowed = allh1cands_sorted[ii*9+1]*allh1cands_sorted[ii*9+1]*sqrt(3.6e-3/allh1cands_sorted[ii*9+2])/(4.5*tobs);
                  double Pdiff_allowed_2 = alll1cands_sorted[jj*9+1]*alll1cands_sorted[jj*9+1]*sqrt(3.6e-3/alll1cands_sorted[jj*9+2])/(4.5*tobs);
                  if (fabs(allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed || fabs(allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2 || fabs(allh1cands_sorted[ii*9+1]-2.0*alll1cands_sorted[jj*9+1])<=Pdiff_allowed || fabs(2.0*allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2 || fabs(allh1cands_sorted[ii*9+1]-3.0*alll1cands_sorted[jj*9+1])<=Pdiff_allowed || fabs(3.0*allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2 || fabs(allh1cands_sorted[ii*9+1]-4.0*alll1cands_sorted[jj*9+1])<=Pdiff_allowed || fabs(4.0*allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2 || fabs(allh1cands_sorted[ii*9+1]-0.5*alll1cands_sorted[jj*9+1])<=Pdiff_allowed || fabs(0.5*allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2 || fabs(allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1]/3.0)<=Pdiff_allowed || fabs(allh1cands_sorted[ii*9+1]/3.0-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2 || fabs(allh1cands_sorted[ii*9+1]-0.25*alll1cands_sorted[jj*9+1])<=Pdiff_allowed || fabs(0.25*allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2) {
                     double absd1mPo2 = fabs(allh1cands_sorted[ii*9+4]-M_PI_2);
                     double absd2mPo2 = fabs(alll1cands_sorted[jj*9+4]-M_PI_2);
                     double dist = acos(sin(absd1mPo2)*sin(absd2mPo2)*cos(allh1cands_sorted[ii*9+3]-alll1cands_sorted[jj*9+3])+cos(absd1mPo2)*cos(absd2mPo2));

                     if (dist<=2.0*skydiff_allowed/(allh1cands_sorted[ii*9]+alll1cands_sorted[jj*9])) {
                        if (bestmatchprob>alll1cands_sorted[jj*9+7]) {
                           bestmatchprob = alll1cands_sorted[jj*9+7];
                           bestmatch = jj;
                        }
                     } //end sky check
                  } //end period check
               } //end modulation depth check
            } //end frequency check
         } //end test against L1 values
         if (bestmatch>=0) {
            fprintf(CANDS, "%f %f %f %f %f %f %g %f %g %d %f %f %f %f %f %f %g %f %g %d\n",  (float)allh1cands_sorted[ii*9], (float)allh1cands_sorted[ii*9+1], (float)allh1cands_sorted[ii*9+2], (float)allh1cands_sorted[ii*9+3], (float)allh1cands_sorted[ii*9+4], (float)allh1cands_sorted[ii*9+5], allh1cands_sorted[ii*9+6], (float)allh1cands_sorted[ii*9+7], allh1cands_sorted[ii*9+8], allh1cands_job_sorted[ii], (float)alll1cands_sorted[bestmatch*9], (float)alll1cands_sorted[bestmatch*9+1], (float)alll1cands_sorted[bestmatch*9+2], (float)alll1cands_sorted[bestmatch*9+3], (float)alll1cands_sorted[bestmatch*9+4], (float)alll1cands_sorted[bestmatch*9+5], alll1cands_sorted[bestmatch*9+6], (float)alll1cands_sorted[bestmatch*9+7], alll1cands_sorted[bestmatch*9+8], alll1cands_job_sorted[bestmatch]);
         }
      } else if ((allh1cands_sorted[ii*9]-alll1cands_sorted[(l1count-1)*9])<=fdiff_allowed && (allh1cands_sorted[ii*9]-alll1cands_sorted[(l1count-1)*9])>=0.0) {
         int bestmatch = -1;
         double bestmatchprob = 0.0;

         jj = l1count-1;
         while (l1count>0 && (allh1cands_sorted[ii*9]-alll1cands_sorted[jj*9])<1.05*fdiff_allowed) jj--;

         for (/* jj value */; jj<l1count; jj++) {
            if (allh1cands_sorted[ii*9]-alll1cands_sorted[jj*9]<-1.05*fdiff_allowed) break;

            if (fabs(allh1cands_sorted[ii*9]-alll1cands_sorted[jj*9])<=fdiff_allowed) {
               if (fabs(allh1cands_sorted[ii*9+2]-alll1cands_sorted[jj*9+2])<=dfdiff_allowed) {
                  double Pdiff_allowed = allh1cands_sorted[ii*9+1]*allh1cands_sorted[ii*9+1]*sqrt(3.6e-3/allh1cands_sorted[ii*9+2])/(4.5*tobs);
                  double Pdiff_allowed_2 = alll1cands_sorted[jj*9+1]*alll1cands_sorted[jj*9+1]*sqrt(3.6e-3/alll1cands_sorted[jj*9+2])/(4.5*tobs);
                  if (fabs(allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed || fabs(allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2 || fabs(allh1cands_sorted[ii*9+1]-2.0*alll1cands_sorted[jj*9+1])<=Pdiff_allowed || fabs(2.0*allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2 || fabs(allh1cands_sorted[ii*9+1]-3.0*alll1cands_sorted[jj*9+1])<=Pdiff_allowed || fabs(3.0*allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2 || fabs(allh1cands_sorted[ii*9+1]-4.0*alll1cands_sorted[jj*9+1])<=Pdiff_allowed || fabs(4.0*allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2 || fabs(allh1cands_sorted[ii*9+1]-0.5*alll1cands_sorted[jj*9+1])<=Pdiff_allowed || fabs(0.5*allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2 || fabs(allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1]/3.0)<=Pdiff_allowed || fabs(allh1cands_sorted[ii*9+1]/3.0-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2 || fabs(allh1cands_sorted[ii*9+1]-0.25*alll1cands_sorted[jj*9+1])<=Pdiff_allowed || fabs(0.25*allh1cands_sorted[ii*9+1]-alll1cands_sorted[jj*9+1])<=Pdiff_allowed_2) {
                     double absd1mPo2 = fabs(allh1cands_sorted[ii*9+4]-M_PI_2);
                     double absd2mPo2 = fabs(alll1cands_sorted[jj*9+4]-M_PI_2);
                     double dist = acos(sin(absd1mPo2)*sin(absd2mPo2)*cos(allh1cands_sorted[ii*9+3]-alll1cands_sorted[jj*9+3])+cos(absd1mPo2)*cos(absd2mPo2));

                     if (dist<=2.0*skydiff_allowed/(allh1cands_sorted[ii*9]+alll1cands_sorted[jj*9])) {
                        if (bestmatchprob>alll1cands_sorted[jj*9+7]) {
                           bestmatchprob = alll1cands_sorted[jj*9+7];
                           bestmatch = jj;
                        }
                     } //end sky check
                  } //end period check
               } //end modulation depth check
            } //end frequency check
         } //end check against L1 values
         if (bestmatch>=0) {
            fprintf(CANDS, "%f %f %f %f %f %f %g %f %g %d %f %f %f %f %f %f %g %f %g %d\n",  (float)allh1cands_sorted[ii*9], (float)allh1cands_sorted[ii*9+1], (float)allh1cands_sorted[ii*9+2], (float)allh1cands_sorted[ii*9+3], (float)allh1cands_sorted[ii*9+4], (float)allh1cands_sorted[ii*9+5], allh1cands_sorted[ii*9+6], (float)allh1cands_sorted[ii*9+7], allh1cands_sorted[ii*9+8], allh1cands_job_sorted[ii], (float)alll1cands_sorted[bestmatch*9], (float)alll1cands_sorted[bestmatch*9+1], (float)alll1cands_sorted[bestmatch*9+2], (float)alll1cands_sorted[bestmatch*9+3], (float)alll1cands_sorted[bestmatch*9+4], (float)alll1cands_sorted[bestmatch*9+5], alll1cands_sorted[bestmatch*9+6], (float)alll1cands_sorted[bestmatch*9+7], alll1cands_sorted[bestmatch*9+8], alll1cands_job_sorted[bestmatch]);
         }
      } //end if H1 candidate is barely outside L1 frequency range
   } //end loop over H1 values

   fclose(CANDS);
   gsl_root_fsolver_free(s);
   free(allh1cands_sorted);
   free(alll1cands_sorted);
   free(allh1cands_job_sorted);
   free(alll1cands_job_sorted);

   return 0;
}
