/*
*  Copyright (C) 2010 Michele Vallisneri, Will Farr, Evan Ochsner, 2014 A. Klein
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
#ifndef _LALSIMIMRSPINALIGNEDEOBGSLOPTIMIZEDINTERPOLATION_C
#define _LALSIMIMRSPINALIGNEDEOBGSLOPTIMIZEDINTERPOLATION_C

#include <stdio.h>
#include <math.h>
#include "stddef.h"

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSpinEOB.h"

/* static size_t optimized_gsl_interp_bsearch(const double x_array[], double x, size_t index_lo, size_t index_hi); */
/* static inline size_t optimized_gsl_interp_accel_find(gsl_interp_accel * a, const double xa[], size_t len, double x); */
/* static inline void optimized_coeff_calc (const double c_array[], double dy, double dx, size_t index, double * b, double * c, double * d); */
/* static int optimized_cspline_eval (const void * vstate, const double x_array[], const double y_array[], size_t size, double x, gsl_interp_accel * a, double *y,unsigned int *index_old, double *x_lo_old,double *y_lo_old,double *b_i_old,double *c_i_old,double *d_i_old); */
/* static int optimized_gsl_spline_eval_e(const gsl_spline * spline, double interptime, gsl_interp_accel * accel, double * output,unsigned int *index_old, double *x_lo_old,double *y_lo_old,double *b_i_old,double *c_i_old,double *d_i_old){ */

typedef struct
{
  double * c;
  double * g;
  double * diag;
  double * offdiag;
} cspline_state_t;

static size_t optimized_gsl_interp_bsearch(const double x_array[], double x, size_t index_lo, size_t index_hi)
{
  size_t ilo = index_lo;
  size_t ihi = index_hi;
  while(ihi > ilo + 1) {
    size_t i = (ihi + ilo)/2;
    if(x_array[i] > x)
      ihi = i;
    else
      ilo = i;
  }

  return ilo;
}

static inline size_t optimized_gsl_interp_accel_find(gsl_interp_accel * a, const double xa[], size_t len, double x)
{
  size_t x_index = a->cache;

  if(x < xa[x_index]) { //OPTIMIZED: This shouldn't occur often, if at all
    a->miss_count++;
    a->cache = optimized_gsl_interp_bsearch(xa, x, 0, x_index);
  }
  else
    {
      if(x >= xa[x_index + 1])
        {
	  if( x>=xa[x_index + 2]) //OPTIMIZED: Check the next one. It'll probably be there.
            {
	      a->miss_count++;
	      a->cache = optimized_gsl_interp_bsearch(xa, x, x_index, len-1);
            }
	  else
            {
	      a->hit_count++;
	      a->cache = x_index + 1;
            }
        }
      else
        {
	  a->hit_count++;
        }
    }
  return a->cache;
}


static inline void optimized_coeff_calc (const double c_array[], double dy, double dx, size_t index, double * b, double * c, double * d)
{
  const double c_i = c_array[index];
  const double c_ip1 = c_array[index + 1];
  *b = (dy / dx) - dx * (c_ip1 + 2.0 * c_i) / 3.0; //OPTIMIZED: It seems c_array is an array of derivatives
  *c = c_i;
  *d = (c_ip1 - c_i) / (3.0 * dx);
}

static int optimized_cspline_eval (const void * vstate, const double x_array[], const double y_array[], size_t size, double x, gsl_interp_accel * a, double *y,unsigned int *index_old, double *x_lo_old,double *y_lo_old,double *b_i_old,double *c_i_old,double *d_i_old)
{
  const cspline_state_t *state = (const cspline_state_t *) vstate;

  double x_lo, x_hi;
  double dx;
  size_t index;

  index = optimized_gsl_interp_accel_find (a, x_array, size, x);
  if(index==*index_old && (*index_old)>0) {
    double delx = x - *x_lo_old;
    *y = *y_lo_old + delx * (*b_i_old + delx * (*c_i_old + delx * *d_i_old));
  } else {
    /* evaluate */
    x_hi = x_array[index + 1];
    x_lo = x_array[index];
    dx = x_hi - x_lo;
    const double y_lo = y_array[index];
    const double y_hi = y_array[index + 1];
    const double dy = y_hi - y_lo;
    double delx = x - x_lo;
    double b_i, c_i, d_i;
    optimized_coeff_calc(state->c, dy, dx, index,  &b_i, &c_i, &d_i);//OPTIMIZED: It seems state->c is an array of derivatives at the indexed locations
    *y = y_lo + delx * (b_i + delx * (c_i + delx * d_i));

    *b_i_old = b_i;
    *c_i_old = c_i;
    *d_i_old = d_i;
    *x_lo_old = x_lo;
    *y_lo_old = y_lo;
    *index_old=index;
  }

  return GSL_SUCCESS;
}

static int optimized_gsl_spline_eval_e(const gsl_spline * spline, double interptime, gsl_interp_accel * accel, double * output,unsigned int *index_old, double *x_lo_old,double *y_lo_old,double *b_i_old,double *c_i_old,double *d_i_old){
  return optimized_cspline_eval(spline->interp->state, spline->x, spline->y, spline->interp->size, interptime, accel, output,index_old,x_lo_old,y_lo_old,b_i_old,c_i_old,d_i_old);
}


static int SEOBNRv2OptimizedInterpolatorIncludeAmpPhase(REAL8Array *yin, REAL8 tinit, REAL8 deltat, UINT4 num_input_times, REAL8Array ** yout)
{
    int errnum = 0;

    /* needed for the integration */
    size_t dim=4;

    /* needed for the final interpolation */
    gsl_spline *interp = NULL;
    gsl_interp_accel *accel = NULL;
    int outputlen = 0;
    REAL8Array *output = NULL;
    REAL8 *times, *vector;      /* aliases */

    /* note: for speed, this replaces the single CALLGSL wrapper applied before each GSL call */
    interp = gsl_spline_alloc(gsl_interp_cspline, num_input_times);
    accel = gsl_interp_accel_alloc();

    //printf("heya?\n");

    outputlen = (int)(yin->data[num_input_times-1] / deltat) + 1;
    output = XLALCreateREAL8ArrayL(2, dim + 3, outputlen);/* Original (dim+1), Optimized (dim+3), since we're adding amp & phase */

    if (!interp || !accel || !output) {
      errnum = XLAL_ENOMEM;   /* ouch again, ran out of memory */
      if (output)
        XLALDestroyREAL8Array(output);
      outputlen = 0;
      goto bail_out;
    }

    /* make an array of times */
    times = output->data;
    for (int j = 0; j < outputlen; j++)
      times[j] = tinit + deltat * j;

    /* interpolate! */
    for (unsigned int i = 1; i <= dim+2; i++) { /* Original (dim), Optimized (dim+2), since we're also interpolating amp & phase */
     //gsl_spline_init(interp, &yin->data[0], &yin->data[num_input_times * i], num_input_times + 1);
     gsl_spline_init(interp, &yin->data[0], &yin->data[num_input_times * i], num_input_times);

      vector = output->data + outputlen * i;
      unsigned int index_old=0;
      double x_lo_old=0,y_lo_old=0,b_i_old=0,c_i_old=0,d_i_old=0;
      for (int j = 0; j < outputlen; j++) {
        optimized_gsl_spline_eval_e(interp,times[j],accel, &(vector[j]),&index_old,&x_lo_old,&y_lo_old,&b_i_old,&c_i_old,&d_i_old);
      }
    }

    /* deallocate stuff and return */
  bail_out:

    if (interp)
        XLAL_CALLGSL(gsl_spline_free(interp));
    if (accel)
        XLAL_CALLGSL(gsl_interp_accel_free(accel));

    if (errnum)
        XLAL_ERROR(errnum);

    *yout = output;
    return outputlen;
}

#endif /* _LALSIMIMRSPINALIGNEDEOBGSLOPTIMIZEDINTERPOLATION_C */
