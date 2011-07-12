/*
 *  Copyright (C) 2011 Evan Goetz
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

#include <lal/LALMalloc.h>

#include <gsl/gsl_roots.h>

#include "upperlimits.h"
#include "statistics.h"
#include "IHS.h"

UpperLimitVector * new_UpperLimitVector(UINT4 length)
{
   
   const CHAR *fn = __func__;
   
   UpperLimitVector *vector = XLALMalloc(sizeof(*vector));
   if (vector==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%zu) failed.\n", fn, sizeof(*vector));
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
   }
   
   vector->length = length;
   if (length==0) vector->data = NULL;
   else {
      vector->data = XLALMalloc( length*sizeof(*vector->data) );
      if (vector->data==NULL) {
         XLALFree((UpperLimitVector*)vector);
         fprintf(stderr,"%s: XLALMalloc(%zu) failed.\n", fn, length*sizeof(*vector->data));
         XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
      }
   }
   
   return vector;
   
} /* new_UpperLimitVector() */



UpperLimitVector * resize_UpperLimitVector(UpperLimitVector *vector, UINT4 length)
{
   
   const CHAR *fn = __func__;
   
   if (vector==NULL) return new_UpperLimitVector(length);
   if (length==0) {
      free_UpperLimitVector(vector);
      return NULL;
   }
   
   vector->data = XLALRealloc(vector->data, length*sizeof(*vector->data));
   if (vector->data==NULL) {
      vector->length = 0;
      fprintf(stderr,"%s: XLALRealloc(%zu) failed.\n", fn, length*sizeof(*vector->data));
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
   }
   vector->length = length;
   
   return vector;
   
} /* resize_UpperLimitVector() */



void free_UpperLimitVector(UpperLimitVector *vector)
{
   
   const CHAR *fn = __func__;
   
   if (vector==NULL) return;
   if ((!vector->length || !vector->data) && (vector->length || vector->data)) XLAL_ERROR_VOID(fn, XLAL_EINVAL);
   if (vector->data) XLALFree((UpperLimit*)vector->data);
   vector->data = NULL;
   XLALFree((UpperLimitVector*)vector);
   return;
   
} /* free_UpperLimitVector() */



//void skypoint95UL(UpperLimit *ul, ihsfarStruct *ihsfarstruct, inputParamsStruct *params, ffdataStruct *ffdata, ihsMaximaStruct *ihsmaxima, REAL4Vector *aveNoise, REAL4Vector *fbinavgs)
void skypoint95UL(UpperLimit *ul, inputParamsStruct *params, ffdataStruct *ffdata, ihsMaximaStruct *ihsmaxima, REAL4Vector *aveNoise, REAL4Vector *fbinavgs)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj, kk;
   
   INT4 minrows = (INT4)round(2.0*params->dfmin*params->Tcoh)+1;
   
   REAL4Vector *twiceAveNoise = XLALCreateREAL4Vector(aveNoise->length);
   if (twiceAveNoise==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   memcpy(twiceAveNoise->data, aveNoise->data, sizeof(REAL4)*aveNoise->length);
   for (ii=0; ii<(INT4)aveNoise->length; ii++) twiceAveNoise->data[ii] *= 2.0;
   
   //Initialize solver
   const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
   gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);
   gsl_function F;
   struct ncx2cdf_solver_params pars;
   
   INT4 totaliterations = 0;
   REAL8 highesth0 = 0.0, fsig = 0.0, period = 0.0, moddepth = 0.0;
   REAL8 dailyharmonic = params->Tobs/(24.0*3600);
   REAL8 dailyharmonic2 = dailyharmonic*2.0, dailyharmonic3 = dailyharmonic*3.0, dailyharmonic4 = dailyharmonic*4.0;
   for (ii=minrows; ii<=ihsmaxima->rows; ii++) {
      REAL8 loudestoutlier = 0.0, loudestoutlierminusnoise = 0.0, loudestoutliernoise = 0.0;
      INT4 jjbinofloudestoutlier = 0, locationofloudestoutlier = 0;
      for (jj=0; jj<ffdata->numfbins-(ii-1); jj++) {
         INT4 locationinmaximavector = (ii-2)*ffdata->numfbins - ((ii-1)*(ii-1)-(ii-1))/2 + jj;
         INT4 location = ihsmaxima->locations->data[locationinmaximavector];
         
         REAL8 noise = aveNoise->data[location] + aveNoise->data[location*2] + aveNoise->data[location*3] + aveNoise->data[location*4] + aveNoise->data[location*5];
         for (kk=1; kk<=5; kk++) if (fabs(dailyharmonic-kk*location)<=1.0 || fabs(dailyharmonic2-kk*location)<=1.0 || fabs(dailyharmonic3-kk*location)<=1.0 || fabs(dailyharmonic4-kk*location)<=1.0) noise -= aveNoise->data[location*kk];
         REAL8 totalnoise = 0.0;
         for (kk=0; kk<ii; kk++) totalnoise += noise*fbinavgs->data[jj+kk];
         REAL8 ihsminusnoise = ihsmaxima->maxima->data[locationinmaximavector] - totalnoise;
         
         if (ihsminusnoise>loudestoutlierminusnoise) {
            loudestoutlier = ihsmaxima->maxima->data[locationinmaximavector];
            loudestoutliernoise = totalnoise;
            loudestoutlierminusnoise = ihsminusnoise;
            locationofloudestoutlier = ihsmaxima->locations->data[locationinmaximavector];
            jjbinofloudestoutlier = jj;
         }
      }
      
      /* REAL4Vector *tempfbinavgs = XLALCreateREAL4Vector(ii);
      if (tempfbinavgs==NULL) {
         fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      memcpy(tempfbinavgs->data, &(fbinavgs->data[jjbinofloudestoutlier]), sizeof(REAL4)*ii);
      REAL4 avenoiseinrange = calcMean(tempfbinavgs);
      XLALDestroyREAL4Vector(tempfbinavgs); */
      
      //REAL8 initialguess = ncx2inv(0.95, 2.0*avenoiseinrange*ihsfarstruct->ihsdistMean->data[ii-2], 2.0*(loudestoutlier-avenoiseinrange*ihsfarstruct->ihsdistMean->data[ii-2]));
      REAL8 initialguess = ncx2inv(0.95, 2.0*loudestoutliernoise, 2.0*loudestoutlierminusnoise);
      if (XLAL_IS_REAL8_FAIL_NAN(initialguess)) {
         //fprintf(stderr, "%s: ncx2inv(%f,%f,%f) failed.\n",fn, 0.95, 2.0*avenoiseinrange*ihsfarstruct->ihsdistMean->data[ii-2], 2.0*(loudestoutlier-avenoiseinrange*ihsfarstruct->ihsdistMean->data[ii-2]));
         fprintf(stderr, "%s: ncx2inv(%f,%f,%f) failed.\n",fn, 0.95, 2.0*loudestoutliernoise, 2.0*loudestoutlierminusnoise);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      REAL8 lo = 0.125*initialguess, hi = 8.0*initialguess;
      pars.val = 2.0*loudestoutlier;
      //pars.dof = 2.0*avenoiseinrange*ihsfarstruct->ihsdistMean->data[ii-2];
      pars.dof = 2.0*loudestoutliernoise;
      pars.ULpercent = 0.95;
      F.function = &gsl_ncx2cdf_float_solver;
      F.params = &pars;
      if (gsl_root_fsolver_set(s, &F, lo, hi) != 0) {
         fprintf(stderr,"%s: gsl_root_fsolver_set() failed.\n", fn);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      
      INT4 status = GSL_CONTINUE;
      INT4 max_iter = 100;
      REAL8 root = 0.0;
      jj = 0;
      while (status==GSL_CONTINUE && jj<max_iter) {
         jj++;
         status = gsl_root_fsolver_iterate(s);
         if (status!=GSL_CONTINUE && status!=GSL_SUCCESS) {
            fprintf(stderr,"%s: gsl_root_fsolver_iterate() failed with code %d.\n", fn, status);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         }
         root = gsl_root_fsolver_root(s);
         //fprintf(stderr, "root = %.8f\n", root);
         lo = gsl_root_fsolver_x_lower(s);
         hi = gsl_root_fsolver_x_upper(s);
         status = gsl_root_test_interval(lo, hi, 0.0, 0.001);
         if (status!=GSL_CONTINUE && status!=GSL_SUCCESS) {
            fprintf(stderr,"%s: gsl_root_test_interval() failed with code %d.\n", fn, status);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         }
      }
      if (status != GSL_SUCCESS) {
         fprintf(stderr, "%s: Root finding iteration (%d/%d) failed with code %d. Current root = %f\n", fn, jj, max_iter, status, root);
         XLAL_ERROR_VOID(fn, XLAL_FAILURE);
      } else if (jj==max_iter) {
         fprintf(stderr, "%s: Root finding failed to converge after %d iterations", fn, jj);
         XLAL_ERROR_VOID(fn, XLAL_FAILURE);
      }
      
      totaliterations += jj;
      
      //REAL8 h0 = ihs2h0(0.5*root+loudestoutlier, locationofloudestoutlier, jjbinofloudestoutlier, ii, params, aveNoise, fbinavgs);
      REAL8 h0 = ihs2h0(root+pars.dof, params);
      if (XLAL_IS_REAL8_FAIL_NAN(h0)) {
         fprintf(stderr, "%s: ihs2h0() failed.\n", fn);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      if (h0>highesth0) {
         highesth0 = h0;
         fsig = params->fmin + (0.5*(ii-1) + jjbinofloudestoutlier)/params->Tcoh;
         period = params->Tobs/locationofloudestoutlier;
         moddepth = 0.5*(ii-1)/params->Tcoh;
      }
      //fprintf(stderr, "%d done\n", ii);
   } /* for ii=minrows --> maximum rows */
   
   ul->ULval = highesth0;
   ul->fsig = fsig;
   ul->period = period;
   ul->moddepth = moddepth;
   ul->iterations2reachUL = totaliterations;
   
   gsl_root_fsolver_free(s);
   XLALDestroyREAL4Vector(twiceAveNoise);
   
}
REAL8 gsl_ncx2cdf_solver(REAL8 x, void *p)
{
   
   struct ncx2cdf_solver_params *params = (struct ncx2cdf_solver_params*)p;
   return ncx2cdf(params->val, params->dof, x) - (1.0-params->ULpercent);
   
}
REAL8 gsl_ncx2cdf_float_solver(REAL8 x, void *p)
{
   
   struct ncx2cdf_solver_params *params = (struct ncx2cdf_solver_params*)p;
   return (REAL8)(ncx2cdf_float((REAL4)params->val, (REAL4)params->dof, (REAL4)x) - (1.0-params->ULpercent));
   
}



void outputUpperLimitsToFile(FILE *outputfile, UpperLimitVector *ulvector)
{
   
   INT4 ii;
   for (ii=0; ii<(INT4)ulvector->length; ii++) {
      fprintf(outputfile,"%.5f %.5f %.6f", ulvector->data[ii].alpha, ulvector->data[ii].delta, ulvector->data[ii].ULval);
   }
   
}



