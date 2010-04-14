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

#include <math.h>

#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_roots.h>

#include <lal/LALConstants.h>
#include <lal/LALMalloc.h>
#include <lal/Window.h>

#include "templates.h"

//////////////////////////////////////////////////////////////
// Allocate memory for farStruct struct  -- done
farStruct * new_farStruct(void)
{

   farStruct *farstruct = (farStruct*)XLALMalloc(sizeof(farStruct));
   farstruct->topRvalues = NULL;

   return farstruct;

}

//////////////////////////////////////////////////////////////
// Destroy farStruct struct  -- done
void free_farStruct(farStruct *farstruct)
{
   
   XLALDestroyREAL8Vector(farstruct->topRvalues);
   farstruct->topRvalues = NULL;
   
   XLALFree((farStruct*)farstruct);

}


//////////////////////////////////////////////////////////////
// Estimate the FAR of the R statistic from the weights
//void estimateFAR(farStruct *out, REAL4Vector *weights, topbinsStruct *topbinsstruct, REAL4 thresh, REAL4Vector *ffplanenoise)
void estimateFAR(farStruct *out, templateStruct *templatestruct, INT4 trials, REAL8 thresh, REAL8Vector *ffplanenoise, REAL8Vector *fbinaveratios)
{
   
   INT4 ii, jj;
   
   REAL8 sumofsqweights = 0.0;
   for (ii=0; ii<(INT4)templatestruct->templatedata->length; ii++) sumofsqweights += (templatestruct->templatedata->data[ii]*templatestruct->templatedata->data[ii]);
   REAL8 sumofsqweightsinv = 1.0/sumofsqweights;
   
   //INT4 trials = (INT4)roundf(100000*0.01/thresh);    //Number of trials to determine FAR value
   REAL8Vector *Rs = XLALCreateREAL8Vector((UINT4)trials);
   
   //RandomParams *param = XLALCreateRandomParams(0);
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
   srand(time(NULL));
   UINT8 randseed = rand();
   gsl_rng_set(rng, randseed);
   
   for (ii=0; ii<trials; ii++) {
      //Create noise value and R value
      REAL8 R = 0.0;
      for (jj=0; jj<(INT4)templatestruct->templatedata->length; jj++) {
         REAL8 noise = expRandNum(ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ], rng);
         R += (noise - ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ])*templatestruct->templatedata->data[jj];
      }
      Rs->data[ii] = R*sumofsqweightsinv;
   }
   REAL8 mean = calcMean(Rs);
   REAL8 sigma = calcStddev(Rs);
   
   //Do an insertion sort. At best this is O(thresh*trials), at worst this is O(thresh*trials*trials).
   if (out->topRvalues == NULL) out->topRvalues = XLALCreateREAL8Vector((UINT4)roundf(thresh*trials)+1);
   out->topRvalues->data[0] = Rs->data[0];
   for (ii=1; ii<(INT4)out->topRvalues->length; ii++) {
      INT4 insertionpoint = ii;
      while (insertionpoint > 0 && Rs->data[ii] > out->topRvalues->data[insertionpoint-1]) insertionpoint--;
      
      for (jj=out->topRvalues->length-1; jj>insertionpoint; jj--) out->topRvalues->data[jj] = out->topRvalues->data[jj-1];
      out->topRvalues->data[insertionpoint] = Rs->data[ii];
   }
   for (ii=out->topRvalues->length; ii<trials; ii++) {
      if (Rs->data[ii] > out->topRvalues->data[out->topRvalues->length - 1]) {
         INT4 insertionpoint = out->topRvalues->length - 1;
         while (insertionpoint > 0 && Rs->data[ii] > out->topRvalues->data[insertionpoint-1]) insertionpoint--;
         
         for (jj=out->topRvalues->length-1; jj>insertionpoint; jj--) out->topRvalues->data[jj] = out->topRvalues->data[jj-1];
         out->topRvalues->data[insertionpoint] = Rs->data[ii];
      }
   }
   
   out->far = out->topRvalues->data[out->topRvalues->length - 1];
   out->distMean = mean;
   out->distSigma = sigma;
   
   //Destroy
   XLALDestroyREAL8Vector(Rs);
   gsl_rng_free(rng);

}



void numericFAR(farStruct *out, templateStruct *templatestruct, REAL8 thresh, REAL8Vector *ffplanenoise)
{
   
   const gsl_root_fdfsolver_type *T = gsl_root_fdfsolver_steffenson;
   gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc(T);
   gsl_function_fdf FDF;
   struct gsl_probR_pars params = {templatestruct, ffplanenoise, thresh};
   
   FDF.f = &gsl_probR;
   FDF.df = &gsl_dprobRdR;
   FDF.fdf = &gsl_probRtimesDprobRdR;
   FDF.params = &params;
   
   REAL8 rootguess = 0.0;
   REAL8 sumwsq = 0.0;
   INT4 ii;
   for (ii=0; ii<(INT4)templatestruct->templatedata->length; ii++) sumwsq += templatestruct->templatedata->data[ii]*templatestruct->templatedata->data[ii];
   for (ii=0; ii<(INT4)templatestruct->templatedata->length; ii++) rootguess += templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]/sumwsq;
   
   gsl_root_fdfsolver_set (s, &FDF, rootguess);
   
   ii = 0;
   INT4 max_iter = 100;
   INT4 status;
   REAL8 root;
   do {
      ii++;
      status = gsl_root_fdfsolver_iterate(s);
      root = rootguess;
      rootguess = gsl_root_fdfsolver_root(s);
      status = gsl_root_test_delta (rootguess, root, 0, thresh*0.1);
   } while (status == GSL_CONTINUE && ii<max_iter);
   
   out->far = rootguess;
   out->distMean = 0.0;
   out->distSigma = 1.0; //TODO: Get the real value of sigma
   
   gsl_root_fdfsolver_free(s);
   
}
REAL8 gsl_probR(REAL8 R, void *param)
{
   
   struct gsl_probR_pars *pars = (struct gsl_probR_pars*)param;
   INT4 ii, jj;
   REAL8 prob = 0.0;
   REAL8 sumwsq = 0.0;
   for (ii=0; ii<(INT4)pars->templatestruct->templatedata->length; ii++) sumwsq += pars->templatestruct->templatedata->data[ii]*pars->templatestruct->templatedata->data[ii];
   
   REAL8 fact1, fact2, sumval, prodval, shiftamt;
   for (ii=0; ii<(INT4)pars->templatestruct->templatedata->length; ii++) {
      fact1 = fact2 = sumval = shiftamt = 0.0;
      prodval = 0.0;
      
      INT4 prodfactorposneg = 1;
      if (ii==0) {
         for (jj=1; jj<(INT4)pars->templatestruct->templatedata->length; jj++) {
            if ((pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ])>0.0) {
               prodval += log(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ]);
               //prodfactorposneg *= 1;
            } else if ((pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ])<0.0) {
               prodval += log(-(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ]));
               prodfactorposneg *= -1;
            }
         }
      } else if (ii==(INT4)pars->templatestruct->templatedata->length-1) {
         for (jj=0; jj<(INT4)pars->templatestruct->templatedata->length-1; jj++) {
            if ((pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ])>0.0) {
               prodval += log(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ]);
               //prodfactorposneg *= 1;
            } else if ((pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ])<0.0) {
               prodval += log(-(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ]));
               prodfactorposneg *= -1;
            }
         }
      } else {
         for (jj=0; jj<ii; jj++) {
            if ((pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ])>0.0) {
               prodval += log(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ]);
               //prodfactorposneg *= 1;
            } else if ((pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ])<0.0) {
               prodval += log(-(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ]));
               prodfactorposneg *= -1;
            }
         }
         for (jj=ii+1; jj<(INT4)pars->templatestruct->templatedata->length; jj++) {
            if ((pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ])>0.0) {
               prodval += log(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ]);
               //prodfactorposneg *= 1;
            } else if ((pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ])<0.0) {
               prodval += log(-(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ]));
               prodfactorposneg *= -1;
            }
         }
      }
      
      prodval = prodfactorposneg*exp(prodval);
      
      if (ii==0) {
         for (jj=1; jj<(INT4)pars->templatestruct->templatedata->length; jj++) sumval += pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ];
      } else if (ii==(INT4)pars->templatestruct->templatedata->length-1) {
         for (jj=0; jj<(INT4)pars->templatestruct->templatedata->length-1; jj++) sumval += pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ];
      } else {
         for (jj=0; jj<ii; jj++) sumval += pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ];
         for (jj=ii+1; jj<(INT4)pars->templatestruct->templatedata->length; jj++) sumval += pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ];
      }
      
      sumval += shiftamt;
      
      fact1 = (-R*sumwsq-sumval)/(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ]);
      fact2 = ((REAL8)pars->templatestruct->templatedata->length-1.0)*log(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ]);
      
      prob += exp(fact1+fact2)/prodval;
   }
   prob /= LAL_E;
   
   return prob-pars->threshold;
   
}
REAL8 gsl_dprobRdR(REAL8 R, void *param)
{
   
   struct gsl_probR_pars *pars = (struct gsl_probR_pars*)param;
   INT4 ii;
   REAL8 sumwsq = 0.0;
   for (ii=0; ii<(INT4)pars->templatestruct->templatedata->length; ii++) sumwsq += pars->templatestruct->templatedata->data[ii]*pars->templatestruct->templatedata->data[ii];
   
   REAL8 dprobRdR = -(gsl_probR(R, pars)+pars->threshold)*sumwsq;
   
   return dprobRdR;
   
}
void gsl_probRtimesDprobRdR(REAL8 R, void *param, REAL8 *probR, REAL8 *dprobRdR)
{
   
   struct gsl_probR_pars *pars = (struct gsl_probR_pars*)param;
   INT4 ii, jj;
   REAL8 prob = 0.0;
   REAL8 sumwsq = 0.0;
   for (ii=0; ii<(INT4)pars->templatestruct->templatedata->length; ii++) sumwsq += pars->templatestruct->templatedata->data[ii]*pars->templatestruct->templatedata->data[ii];
   
   REAL8 fact1, fact2, sumval, prodval, shiftamt;
   for (ii=0; ii<(INT4)pars->templatestruct->templatedata->length; ii++) {
      fact1 = fact2 = sumval = shiftamt = 0.0;
      prodval = 0.0;
      
      INT4 prodfactorposneg = 1;
      if (ii==0) {
         for (jj=1; jj<(INT4)pars->templatestruct->templatedata->length; jj++) {
            if ((pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ])>0.0) {
               prodval += log(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ]);
               //prodfactorposneg *= 1;
            } else if ((pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ])<0.0) {
               prodval += log(-(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ]));
               prodfactorposneg *= -1;
            }
         }
      } else if (ii==(INT4)pars->templatestruct->templatedata->length-1) {
         for (jj=0; jj<(INT4)pars->templatestruct->templatedata->length-1; jj++) {
            if ((pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ])>0.0) {
               prodval += log(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ]);
               //prodfactorposneg *= 1;
            } else if ((pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ])<0.0) {
               prodval += log(-(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ]));
               prodfactorposneg *= -1;
            }
         }
      } else {
         for (jj=0; jj<ii; jj++) {
            if ((pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ])>0.0) {
               prodval += log(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ]);
               //prodfactorposneg *= 1;
            } else if ((pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ])<0.0) {
               prodval += log(-(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ]));
               prodfactorposneg *= -1;
            }
         }
         for (jj=ii+1; jj<(INT4)pars->templatestruct->templatedata->length; jj++) {
            if ((pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ])>0.0) {
               prodval += log(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ]);
               //prodfactorposneg *= 1;
            } else if ((pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ])<0.0) {
               prodval += log(-(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ] - pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ]));
               prodfactorposneg *= -1;
            }
         }
      }
      
      prodval = prodfactorposneg*exp(prodval);
      
      if (ii==0) {
         for (jj=1; jj<(INT4)pars->templatestruct->templatedata->length; jj++) sumval += pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ];
      } else if (ii==(INT4)pars->templatestruct->templatedata->length-1) {
         for (jj=0; jj<(INT4)pars->templatestruct->templatedata->length-1; jj++) sumval += pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ];
      } else {
         for (jj=0; jj<ii; jj++) sumval += pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ];
         for (jj=ii+1; jj<(INT4)pars->templatestruct->templatedata->length; jj++) sumval += pars->templatestruct->templatedata->data[jj]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[jj] ];
      }
      
      sumval += shiftamt;
      
      fact1 = (-R*sumwsq-sumval)/(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ]);
      fact2 = ((REAL8)pars->templatestruct->templatedata->length-1.0)*log(pars->templatestruct->templatedata->data[ii]*pars->ffplanenoise->data[ pars->templatestruct->secondfftfrequencies->data[ii] ]);
      
      prob += exp(fact1+fact2)/prodval;
   }
   prob /= LAL_E;
   
   *probR = prob-pars->threshold;
   *dprobRdR = -(prob+pars->threshold)*sumwsq;
   
}


//////////////////////////////////////////////////////////////
// Analytically calculate the probability of a true signal
REAL8 probR(templateStruct *templatestruct, REAL8Vector *ffplanenoise, REAL8Vector *fbinaveratios, REAL8 R)
{
   
   INT4 ii, jj;
   REAL8 prob = 0.0;
   REAL8 sumwsq = 0.0;
   for (ii=0; ii<(INT4)templatestruct->templatedata->length; ii++) sumwsq += templatestruct->templatedata->data[ii]*templatestruct->templatedata->data[ii];
   
   REAL8 fact1, fact2, sumval, prodval, shiftamt;
   for (ii=0; ii<(INT4)templatestruct->templatedata->length; ii++) {
      fact1 = fact2 = sumval = shiftamt = 0.0;
      prodval = 0.0;
      
      INT4 prodfactorposneg = 1;
      if (ii==0) {
         for (jj=1; jj<(INT4)templatestruct->templatedata->length; jj++) {
            if ((templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ] - templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ])>0.0) {
               prodval += log(templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ] - templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ]);
               //prodfactorposneg *= 1;
            } else if ((templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ] - templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ])<0.0) {
               prodval += log(-(templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ] - templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ]));
               prodfactorposneg *= -1;
            }
         }
      } else if (ii==(INT4)templatestruct->templatedata->length-1) {
         for (jj=0; jj<(INT4)templatestruct->templatedata->length-1; jj++) {
            if ((templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ] - templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ])>0.0) {
               prodval += log(templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ] - templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ]);
               //prodfactorposneg *= 1;
            } else if ((templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ] - templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ])<0.0) {
               prodval += log(-(templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ] - templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ]));
               prodfactorposneg *= -1;
            }
         }
      } else {
         for (jj=0; jj<ii; jj++) {
            if ((templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ] - templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ])>0.0) {
               prodval += log(templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ] - templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ]);
               //prodfactorposneg *= 1;
            } else if ((templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ] - templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ])<0.0) {
               prodval += log(-(templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ] - templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ]));
               prodfactorposneg *= -1;
            }
         }
         for (jj=ii+1; jj<(INT4)templatestruct->templatedata->length; jj++) {
            if ((templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ] - templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ])>0.0) {
               prodval += log(templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ] - templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ]);
               //prodfactorposneg *= 1;
            } else if ((templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ] - templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ])<0.0) {
               prodval += log(-(templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ] - templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ]));
               prodfactorposneg *= -1;
            }
         }
      }
      
      prodval = prodfactorposneg*exp(prodval);
      
      if (ii==0) {
         for (jj=1; jj<(INT4)templatestruct->templatedata->length; jj++) sumval += templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ];
      } else if (ii==(INT4)templatestruct->templatedata->length-1) {
         for (jj=0; jj<(INT4)templatestruct->templatedata->length-1; jj++) sumval += templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ];
      } else {
         for (jj=0; jj<ii; jj++) sumval += templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ];
         for (jj=ii+1; jj<(INT4)templatestruct->templatedata->length; jj++) sumval += templatestruct->templatedata->data[jj]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ];
      }
      
      sumval += shiftamt;
      
      fact1 = (-R*sumwsq-sumval)/(templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ]);
      fact2 = ((REAL8)templatestruct->templatedata->length-1.0)*log(templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ]);
      
      prob += exp(fact1+fact2)/prodval;
   }
   prob /= LAL_E;
   
   prob = log10(prob);
   
   return prob;
   
}


templateStruct * new_templateStruct(INT4 length)
{
   
   INT4 ii;
   
   templateStruct *templatestruct = XLALMalloc(sizeof(templateStruct));
   
   templatestruct->templatedata = XLALCreateREAL8Vector((UINT4)length);
   templatestruct->pixellocations = XLALCreateINT4Vector((UINT4)length);
   templatestruct->firstfftfrequenciesofpixels = XLALCreateINT4Vector((UINT4)length);
   templatestruct->secondfftfrequencies = XLALCreateINT4Vector((UINT4)length);
   for (ii=0; ii<length; ii++) {
      templatestruct->templatedata->data[ii] = 0.0;
      templatestruct->pixellocations->data[ii] = 0;
      templatestruct->firstfftfrequenciesofpixels->data[ii] = 0;
      templatestruct->secondfftfrequencies->data[ii] = 0;
   }
   
   return templatestruct;
   
}


void free_templateStruct(templateStruct *nameoftemplate)
{
   
   XLALDestroyREAL8Vector(nameoftemplate->templatedata);
   XLALDestroyINT4Vector(nameoftemplate->pixellocations);
   XLALDestroyINT4Vector(nameoftemplate->firstfftfrequenciesofpixels);
   XLALDestroyINT4Vector(nameoftemplate->secondfftfrequencies);
   
   XLALFree((templateStruct*)nameoftemplate);
   
}




//////////////////////////////////////////////////////////////
// Make an estimated template based on FFT of train of Gaussians
//void makeTemplateGaussians(ffdataStruct *out, candidate *in)
void makeTemplateGaussians(templateStruct *out, candidate *in, inputParamsStruct *params)
{

   INT4 ii, jj, kk, numfbins, numffts, N;
   
   numfbins = (INT4)(round(params->fspan*params->Tcoh)+1);   //Number of frequency bins
   numffts = (INT4)floor(2*(params->Tobs/params->Tcoh)-1);     //Number of FFTs
   N = (INT4)floor(params->Tobs/in->period);     //Number of Gaussians
   
   REAL8 periodf = 1.0/in->period;
   
   //Set up frequencies and determine separation in time of peaks for each frequency
   REAL8Vector *phi_actual = XLALCreateREAL8Vector((UINT4)numfbins);
   for (ii=0; ii<(INT4)phi_actual->length; ii++) {
      //out->f->data[ii] = in->fmin + ii/in->Tcoh;
      if ( fabs(params->fmin + ii/params->Tcoh - in->fsig)/in->moddepth <= 1.0 ) {
         phi_actual->data[ii] = 0.5*in->period - asin(fabs(params->fmin + ii/params->Tcoh - in->fsig)/
            in->moddepth)*LAL_1_PI*in->period;
      } else {
         phi_actual->data[ii] = 0.0;
      }
   }
   
   //Create second FFT frequencies
   REAL8Vector *fpr = XLALCreateREAL8Vector((UINT4)floor(numffts*0.5)+1);
   for (ii=0; ii<(INT4)fpr->length; ii++) fpr->data[ii] = (REAL8)ii/params->Tobs;
   
   //Scale used for "spillover" into bins outside of phi_actual
   REAL8 k = in->moddepth*params->Tcoh;    //amplitude of modulation in units of bins
   REAL8Vector *scale = XLALCreateREAL8Vector((UINT4)numfbins);      //the scaling factor
   INT4 m0 = (INT4)round(in->fsig*params->Tcoh) - (INT4)round(params->fmin*params->Tcoh);   //central frequency bin
   INT4 mextent = (INT4)floor(in->moddepth*params->Tcoh);   //Bins filled by modulation
   REAL8 overage = (k-(REAL8)mextent)-1.0;
   INT4 fnumstart = -1;
   INT4 fnumend = -1;
   for (ii=0; ii<(INT4)scale->length; ii++) {
      if (mextent != 0) {
         if (ii < m0-mextent-2 || ii > m0+mextent+2) {
            scale->data[ii] = 0.0;
         } else if (ii == m0-mextent-2 || ii == m0+mextent+2) {
            scale->data[ii] = sincxoverxsqminusone(overage-1)*sincxoverxsqminusone(overage-1);
         } else if (ii == m0-mextent-1 || ii == m0+mextent+1) {
            scale->data[ii] = sincxoverxsqminusone(overage)*sincxoverxsqminusone(overage);
         } else {
            scale->data[ii] = 1.0;
         }
      } else {
         if (ii < m0-2 || ii > m0+2) {
            scale->data[ii] = 0.0;
         } else if (ii == m0-2 || ii == m0+2) {
            scale->data[ii] = sincxoverxsqminusone(overage-1)*sincxoverxsqminusone(overage-1);
         } else if (ii == m0-1 || ii == m0+1) {
            scale->data[ii] = sincxoverxsqminusone(overage)*sincxoverxsqminusone(overage);
         } else {
            scale->data[ii] = 1.0;
         }
      }
   }
   for (ii=0; ii<(INT4)scale->length; ii++) {
      if (scale->data[ii] != 0.0 && fnumstart == -1) fnumstart = ii;
      if (scale->data[ii] == 0.0 && fnumstart != -1 && fnumend==-1) fnumend = ii-1;
   }
   if (fnumend==-1) {
      exit(-1);
   }
   
   //Make sigmas for each frequency
   REAL8Vector *sigmas = XLALCreateREAL8Vector((UINT4)(fnumend-fnumstart+1));
   REAL8Vector *wvals = XLALCreateREAL8Vector((UINT4)floor(2.0*in->period/params->Tcoh));
   REAL8Vector *allsigmas = XLALCreateREAL8Vector(wvals->length * sigmas->length);
   for (ii=0; ii<(INT4)wvals->length; ii++) {         //t = (ii+1)*in->Tcoh*0.5
      REAL8 sigbin = (in->moddepth*cos(LAL_TWOPI*periodf*((ii+1)*params->Tcoh*0.5))+in->fsig)*params->Tcoh;
      REAL8 sigbinvelocity = fabs(-in->moddepth*sin(LAL_TWOPI*periodf*((ii+1)*params->Tcoh*0.5))*params->Tcoh*0.5*params->Tcoh*LAL_TWOPI*periodf);
      REAL8 sigma = 0.5 * params->Tcoh * ((383.85*LAL_1_PI)*(0.5*6.1e-3) / ((sigbinvelocity+0.1769)*(sigbinvelocity+0.1769)+(0.5*6.1e-3)*(0.5*6.1e-3)) + 0.3736);   //Derived fit from simulation
      for (jj=0; jj<(INT4)sigmas->length; jj++) {
         allsigmas->data[ii*sigmas->length + jj] = sincxoverxsqminusone(sigbin-round(params->fmin*params->Tcoh+jj+fnumstart))*sincxoverxsqminusone(sigbin-round(params->fmin*params->Tcoh+jj+fnumstart))*sigma;
      }
   }
   for (ii=0; ii<(INT4)sigmas->length; ii++) {
      for (jj=0; jj<(INT4)wvals->length; jj++) wvals->data[jj] = allsigmas->data[ii + jj*sigmas->length]*allsigmas->data[ii + jj*sigmas->length];
      sigmas->data[ii] = sqrt(calcMean(wvals));
   }
   
   
   //Create template
   REAL8 sum = 0.0;
   //REAL8Vector *fulltemplate = XLALCreateREAL8Vector(sigmas->length*fpr->length);
   REAL8 dataval;
   for (ii=0; ii<(INT4)sigmas->length; ii++) {
      REAL8 s = sigmas->data[ii];
      REAL8 scale1 = 1.0/(1.0+exp(-phi_actual->data[ii+fnumstart]*phi_actual->data[ii+fnumstart]*0.5/(s*s)));
      for (jj=0; jj<(INT4)fpr->length; jj++) {
         
         if (jj==0 || jj==1) {
            dataval = 0.0;
         } else if (fabs(cos(in->period*LAL_TWOPI*fpr->data[jj])-1.0)<1e-5) {
            dataval = scale->data[ii+fnumstart] * scale1 * 2.0 * LAL_TWOPI * s * s * exp(-s * s * LAL_TWOPI * LAL_TWOPI * fpr->data[jj] * fpr->data[jj]) * (cos(phi_actual->data[ii+fnumstart] * LAL_TWOPI * fpr->data[jj]) + 1.0) * N * N;
         } else {
            dataval = scale->data[ii+fnumstart] * scale1 * 2.0 * LAL_TWOPI * s * s * exp(-s * s * LAL_TWOPI * LAL_TWOPI * fpr->data[jj] * fpr->data[jj]) * (cos(N * in->period * LAL_TWOPI * fpr->data[jj]) - 1.0) * (cos(phi_actual->data[ii+fnumstart] * LAL_TWOPI * fpr->data[jj]) + 1.0) / (cos(in->period * LAL_TWOPI * fpr->data[jj]) - 1.0);
         }
         
         //Set any bin below 1e-12 to 0.0 and the DC bins (jj=0 and jj=1) to 0.0
         //if (fulltemplate->data[ii*fpr->length + jj] <= 1e-12 || jj==0 || jj==1) fulltemplate->data[ii*fpr->length + jj] = 0.0;
         if (dataval <= 1e-12 || jj==0 || jj==1) dataval = 0.0;
         
         //Sum up the weights in total
         //sum += fulltemplate->data[ii*fpr->length + jj];
         sum += dataval;
         
         //Compare with weakest top bins and if larger, launch a search to find insertion spot
         if (dataval > out->templatedata->data[out->templatedata->length-1]) {
            INT4 insertionpoint = (INT4)out->templatedata->length-1;
            while (insertionpoint > 0 && dataval > out->templatedata->data[insertionpoint-1]) insertionpoint--;
            
            for (kk=out->templatedata->length-1; kk>insertionpoint; kk--) {
               out->templatedata->data[kk] = out->templatedata->data[kk-1];
               out->pixellocations->data[kk] = out->pixellocations->data[kk-1];
               out->firstfftfrequenciesofpixels->data[kk] = out->firstfftfrequenciesofpixels->data[kk-1];
               out->secondfftfrequencies->data[kk] = out->secondfftfrequencies->data[kk-1];
            }
            out->templatedata->data[insertionpoint] = dataval;
            out->pixellocations->data[insertionpoint] = (ii+fnumstart)*fpr->length + jj;
            out->firstfftfrequenciesofpixels->data[insertionpoint] = ii+fnumstart;
            out->secondfftfrequencies->data[insertionpoint] = jj;
         }
      }
   }
   
   //Normalize
   for (ii=0; ii<(INT4)out->templatedata->length; ii++) out->templatedata->data[ii] /= sum;
   
   //Destroy variables
   XLALDestroyREAL8Vector(phi_actual);
   XLALDestroyREAL8Vector(scale);
   XLALDestroyREAL8Vector(sigmas);
   XLALDestroyREAL8Vector(allsigmas);
   XLALDestroyREAL8Vector(wvals);
   XLALDestroyREAL8Vector(fpr);

}


//////////////////////////////////////////////////////////////
// Make an template based on FFT of sinc squared functions  -- done
//void makeTemplate(ffdataStruct *out, candidate *in, REAL4FFTPlan *plan)
void makeTemplate(templateStruct *out, candidate *in, inputParamsStruct *params, REAL8FFTPlan *plan)
{
   
   INT4 ii, jj, kk, numfbins, numffts;
   
   numfbins = (INT4)(round(params->fspan*params->Tcoh)+1);   //Number of frequency bins
   numffts = (INT4)floor(2*(params->Tobs/params->Tcoh)-1);     //Number of FFTs
   
   REAL8Vector *psd1 = XLALCreateREAL8Vector((UINT4)(numfbins*numffts));
   INT4Vector *freqbins = XLALCreateINT4Vector((UINT4)numfbins);
   
   REAL8 periodf = 1.0/in->period;
   REAL8 B = in->moddepth*params->Tcoh;
   
   //Bin numbers of the frequencies
   for (ii=0; ii<numfbins; ii++) freqbins->data[ii] = (INT4)roundf(params->fmin*params->Tcoh) + ii;
   
   //Determine the signal modulation in bins with time at center of coherence time and create
   //Hann windowed PSDs
   for (ii=0; ii<numffts; ii++) {
      REAL8 t = 0.5*params->Tcoh*ii;
      REAL8 n0 = B*sin(LAL_TWOPI*periodf*t) + in->fsig*params->Tcoh;
      for (jj=0; jj<numfbins; jj++) {
         //Create windowed PSD values
         if ( fabs(n0-freqbins->data[jj]) <= 5.0 ) psd1->data[ii*numfbins + jj] = 2.0/3.0*params->Tcoh*sincxoverxsqminusone(n0-freqbins->data[jj])*sincxoverxsqminusone(n0-freqbins->data[jj]);
         else psd1->data[ii*numfbins + jj] = 0.0;
      }
   }
   
   
   //Do the second FFT
   REAL8Vector *x = XLALCreateREAL8Vector((UINT4)numffts);
   REAL8Window *win = XLALCreateHannREAL8Window(x->length);
   REAL8 winFactor = 8.0/3.0;
   REAL8Vector *psd = XLALCreateREAL8Vector((UINT4)floor(x->length*0.5)+1);
   REAL8 sum = 0.0;
   INT4 doSecondFFT;
   //First loop over frequencies
   for (ii=0; ii<numfbins; ii++) {
      //Set doSecondFFT check flag to 0. Value becomes 1 if at least one element in frequency row is non-zero
      doSecondFFT = 0;
   
      //Next, loop over times
      for (jj=0; jj<(INT4)x->length; jj++) {
         //Pick the right frequency bin for each FFT
         x->data[jj] = psd1->data[ii+jj*numfbins];
         
         //Check, do we need to do the second FFT...?
         if (doSecondFFT==0 && x->data[jj]>0.0) doSecondFFT = 1;
         
         //window
         x->data[jj] *= win->data->data[jj];
      }
      
      //Make the FFT
      INT4 check = 0;
      if (doSecondFFT==1) check = XLALREAL8PowerSpectrum(psd,x,plan);
      if (check != 0) printf("Something wrong with second PSD...\n");
      
      //Scale the data points by 1/N and window factor and (1/fs)
      //Order of vector is by second frequency then first frequency
      //Ignore the DC and 1st frequency bins
      if (doSecondFFT==1) {
         for (jj=2; jj<(INT4)psd->length; jj++) {
            
            REAL8 correctedValue = psd->data[jj]*winFactor/x->length*0.5*params->Tcoh;
            
            sum += correctedValue;
            
            //If value is largest than smallest logged bin, then launch a simple search to find the place to insert it
            if (correctedValue > out->templatedata->data[out->templatedata->length-1]) {
               INT4 insertionpoint = (INT4)out->templatedata->length-1;
               while (insertionpoint > 0 && correctedValue > out->templatedata->data[insertionpoint-1]) insertionpoint--;
               
               for (kk=out->templatedata->length-1; kk>insertionpoint; kk--) {
                  out->templatedata->data[kk] = out->templatedata->data[kk-1];
                  out->pixellocations->data[kk] = out->pixellocations->data[kk-1];
                  out->firstfftfrequenciesofpixels->data[kk] = out->firstfftfrequenciesofpixels->data[kk-1];
                  out->secondfftfrequencies->data[kk] = out->secondfftfrequencies->data[kk-1];
               }
               out->templatedata->data[insertionpoint] = correctedValue;
               out->pixellocations->data[insertionpoint] = ii*psd->length + jj;
               out->firstfftfrequenciesofpixels->data[insertionpoint] = ii;
               out->secondfftfrequencies->data[insertionpoint] = jj;
            }
         }
      }
      
   }
   
   //Normalize
   for (ii=0; ii<(INT4)out->templatedata->length; ii++) out->templatedata->data[ii] /= sum;
   
   XLALDestroyREAL8Vector(psd1);
   XLALDestroyINT4Vector(freqbins);
   XLALDestroyREAL8Vector(x);
   XLALDestroyREAL8Window(win);
   XLALDestroyREAL8Vector(psd);
   
}




//////////////////////////////////////////////////////////////
// Calculates y = sin(pi*x)/(pi*x)/(x^2-1)
REAL8 sincxoverxsqminusone(REAL8 x)
{
   
   REAL8 val;
   
   if (x==1.0 || x==-1.0) val = -0.5;
   else val = gsl_sf_sinc(x)/(x*x-1);
   
   return val;
   
}



