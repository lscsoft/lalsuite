#ifndef _LALSIMIMRSPINPRECEOBGSLOPTIMIZEDINTERPOLATION_C
#define _LALSIMIMRSPINPRECEOBGSLOPTIMIZEDINTERPOLATION_C

#include "LALSimIMRSpinAlignedEOBOptimizedInterpolatorGeneral.c"

/*------------------------------------------------------------------------------------------
 *
 *    Definitions of functions (only one in this file, so no prototypes needed).
 *
 *------------------------------------------------------------------------------------------
 */
static int SEOBNRv3OptimizedInterpolatorGeneral(
                REAL8 * yin, /**<< Data to be interpolated; time first */
                REAL8 tinit, /**<< time at which to begin interpolating */
                REAL8 deltat, /**<< Spacing between interpolated times */
                UINT4 num_input_times, /**<< The number of input times */
                REAL8Array ** yout, /**<< Interpolation output */
                size_t dim /**<< Number of quantities interpolated (e.g. if yin = {t,x,y,z} then dim 3) */
                )
{
    int errnum = 0;

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

    outputlen = (int)(yin[num_input_times-1] / deltat) + 1;

    output = XLALCreateREAL8ArrayL(2, dim+1, outputlen);/* Original (dim+1)*/

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
    for (unsigned int i = 1; i <= dim; i++) { /* Original (dim)  */
     //gsl_spline_init(interp, &yin->data[0], &yin->data[num_input_times * i], num_input_times + 1);
     gsl_spline_init(interp, yin, &(yin[num_input_times * i]), num_input_times);

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

#endif // _LALSIMIMRSPINPRECEOBGSLOPTIMIZEDINTERPOLATION_C
