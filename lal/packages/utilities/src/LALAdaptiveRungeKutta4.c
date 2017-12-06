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

#include <lal/LALAdaptiveRungeKutta4.h>

ark4GSLIntegrator *XLALAdaptiveRungeKutta4Init( int dim,
                             int (* dydt) (double t, const double y[], double dydt[], void * params),  /* These are XLAL functions! */
                             int (* stop) (double t, const double y[], double dydt[], void * params),
                             double eps_abs, double eps_rel
                  				 )
{
  ark4GSLIntegrator *integrator;

	/* allocate our custom integrator structure */
  if (!(integrator = (ark4GSLIntegrator *) LALCalloc(1, sizeof(ark4GSLIntegrator))))
	{
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	}

  /* allocate the GSL ODE components */
  XLAL_CALLGSL( integrator->step    = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45,dim) );
  XLAL_CALLGSL( integrator->control = gsl_odeiv_control_y_new(eps_abs,eps_rel) );
  XLAL_CALLGSL( integrator->evolve  = gsl_odeiv_evolve_alloc(dim) );

  /* allocate the GSL system (functions, etc.) */
	integrator->sys = (gsl_odeiv_system *)LALCalloc(1, sizeof(gsl_odeiv_system));

  /* if something failed to be allocated, bail out */
  if ( !(integrator->step) || !(integrator->control) || !(integrator->evolve) || !(integrator->sys) )
  {
    XLALAdaptiveRungeKutta4Free(integrator);
    XLAL_ERROR_NULL(XLAL_ENOMEM);
  }

  integrator->dydt = dydt;
  integrator->stop = stop;

  integrator->sys->function  = dydt;
  integrator->sys->jacobian  = NULL;
  integrator->sys->dimension = dim;
  integrator->sys->params    = NULL;

  integrator->retries = 6;
	integrator->stopontestonly = 0;

  return integrator;
}


void XLALAdaptiveRungeKutta4Free( ark4GSLIntegrator *integrator )
{
  if (!integrator) return;

  if (integrator->evolve)  XLAL_CALLGSL( gsl_odeiv_evolve_free(integrator->evolve) );
  if (integrator->control) XLAL_CALLGSL( gsl_odeiv_control_free(integrator->control) );
  if (integrator->step)    XLAL_CALLGSL( gsl_odeiv_step_free(integrator->step) );

  LALFree( integrator->sys );
  LALFree( integrator );

  return;
}

/* Local function to store interpolated step in output array */
static int storeStateInOutput(REAL8Array **output, REAL8 t, REAL8 *y, size_t dim, int *outputlen, int count) {
  REAL8Array *out = *output;
  int len = *outputlen;
  size_t i;

  if (count > len) {
    /* Resize array! Have to make a new array, and copy over. */
    REAL8Array *new = XLALCreateREAL8ArrayL(2, dim+1, 2*len);

    if (!new) {
      return XLAL_ENOMEM;
    }

    for (i = 0; i < dim+1; i++) {
      memcpy(&(new->data[2*i*len]), &(out->data[i*len]), len*sizeof(REAL8));
    }

    XLALDestroyREAL8Array(out);
    out = new;
    len *= 2;
  }

  /* Store the current step. */
  out->data[count-1] = t;
  for (i = 1; i < dim+1; i++) {
    out->data[i*len + count - 1] = y[i-1];
  }

  *output = out;
  *outputlen = len;
  return GSL_SUCCESS;
}

/* Local function to shrink output array to proper size before it's returned */
static int shrinkOutput(REAL8Array **output, int *outputlen, int count, size_t dim) {
  REAL8Array *out = *output;
  int len = *outputlen;

  REAL8Array *new = XLALCreateREAL8ArrayL(2, dim+1, count);

  if (!new) {
    return XLAL_ENOMEM;
  }

  size_t i;
  for (i = 0; i < dim+1; i++) {
    memcpy(&(new->data[i*count]), &(out->data[i*len]), count*sizeof(REAL8));
  }

  *output = new;
  *outputlen = count;

  XLALDestroyREAL8Array(out);

  return GSL_SUCCESS;
}

/* Copied from GSL rkf45.c */
typedef struct
{
  double *k1;
  double *k2;
  double *k3;
  double *k4;
  double *k5;
  double *k6;
  double *y0;
  double *ytmp;
}
rkf45_state_t;

/**
 * Fourth-order Runge-Kutta ODE integrator using Runge-Kutta-Fehlberg (RKF45)
 * steps with adaptive step size control.  Intended for use in various
 * waveform generation routines such as SpinTaylorT4 and various EOB models.
 *
 * The method is described in
 *
 * Abramowitz & Stegun, Handbook of Mathematical Functions, Tenth Printing,
 * National Bureau of Standards, Washington, DC, 1972
 * (available online at http://people.math.sfu.ca/~cbm/aands/ )
 *
 * This function also includes "on-the-fly" interpolation of the
 * differential equations at regular intervals in-between integration
 * steps. This "on-the-fly" interpolation method is derived and
 * described in the Mathematica notebook "RKF_with_interpolation.nb";
 * see
 * https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/InspiralPipelineDevelopment/120312111836InspiralPipelineDevelopmentImproved%20Adaptive%20Runge-Kutta%20integrator
 *
 * This method is functionally equivalent to XLALAdaptiveRungeKutta4,
 * but is nearly always faster due to the improved interpolation.
 */
int XLALAdaptiveRungeKutta4Hermite( ark4GSLIntegrator *integrator,	/**< struct holding dydt, stopping test, stepper, etc. */
                                    void *params,			/**< params struct used to compute dydt and stopping test */
                                    REAL8 *yinit,			/**< pass in initial values of all variables - overwritten to final values */
                                    REAL8 tinit,			/**< integration start time */
                                    REAL8 tend_in,			/**< maximum integration time */
                                    REAL8 deltat,			/**< step size for evenly sampled output */
                                    REAL8Array **yout			/**< array holding the evenly sampled output */
                                    )
{
  int status;
  size_t dim, retries, i;
  int outputlen = 0, count = 0;

  REAL8Array *output = NULL;

  REAL8 t, tintp, h;

  REAL8 *ytemp = NULL;

  REAL8 tend = tend_in;

  /* If want to stop only on test, then tend = +/-infinity; otherwise
     tend_in */
  if (integrator->stopontestonly) {
    if (tend < tinit) tend = -1.0/0.0;
    else tend = 1.0/0.0;
  }

  dim = integrator->sys->dimension;

  outputlen = ((int)(tend_in - tinit)/deltat);
  //if (outputlen < 0) outputlen = -outputlen;
  if (outputlen < 0) {
    XLALPrintError("XLAL Error - %s: (tend_in - tinit) and deltat must have the same sign\ntend_in: %f, tinit: %f, deltat: %f\n", __func__, tend_in, tinit, deltat);
    status = XLAL_EINVAL;
    goto bail_out;
  }
  outputlen += 2;
  
  output = XLALCreateREAL8ArrayL(2, (dim+1), outputlen);

  if (!output) {
    status = XLAL_ENOMEM;
    goto bail_out;
  }

  ytemp = XLALCalloc(dim, sizeof(REAL8));

  if (!ytemp) {
    status = XLAL_ENOMEM;
    goto bail_out;
  }

  /* Setup. */
  integrator->sys->params = params;
  integrator->returncode = 0;
  retries = integrator->retries;
  t = tinit;
  tintp = tinit;
  h = deltat;

  /* Copy over first step. */
  output->data[0] = tinit;
  for (i = 1; i <= dim; i++) output->data[i*outputlen] = yinit[i-1];
  count = 1;

  XLAL_BEGINGSL;

  /* We are starting a fresh integration; clear GSL step and evolve
     objects. */
  gsl_odeiv_step_reset(integrator->step);
  gsl_odeiv_evolve_reset(integrator->evolve);

  /* Enter evolution loop.  NOTE: we *always* take at least one
     step. */
  while (1) {
    REAL8 told = t;

    status = gsl_odeiv_evolve_apply(integrator->evolve, integrator->control, integrator->step, integrator->sys, &t, tend, &h, yinit);

    /* Check for failure, retry if haven't retried too many times
       already. */
    if (status != GSL_SUCCESS) {
      if (retries--) {
        /* Retries to spare; reduce h, try again.*/
        h /= 10.0;
        continue;
      } else {
        /* Out of retries, bail with status code. */
        integrator->returncode = status;
        break;
      }
    } else {
      /* Successful step, reset retry counter. */
      retries = integrator->retries;
    }

    /* Now interpolate until we would go past the current integrator time, t.
     * Note we square to get an absolute value, because we may be 
     * integrating t in the positive or negative direction */
    while ( (tintp + deltat)*(tintp + deltat) < t*t) {
      tintp += deltat;

      /* tintp = told + (t-told)*theta, 0 <= theta <= 1.  We have to
         compute h = (t-told) because the integrator returns a
         suggested next h, not the actual stepsize taken. */
      REAL8 hUsed = t - told;
      REAL8 theta = (tintp - told)/hUsed;

      /* These are the interpolating coefficients for y(t + h*theta) =
         ynew + i1*h*k1 + i5*h*k5 + i6*h*k6 + O(h^4). */
      REAL8 i0 = 1.0 + theta*theta*(3.0-4.0*theta);
      REAL8 i1 = -theta*(theta-1.0);
      REAL8 i6 = -4.0*theta*theta*(theta-1.0);
      REAL8 iend = theta*theta*(4.0*theta - 3.0);

      /* Grab the k's from the integrator state. */
      rkf45_state_t *rkfState = integrator->step->state;
      REAL8 *k1 = rkfState->k1;
      REAL8 *k6 = rkfState->k6;
      REAL8 *y0 = rkfState->y0;

      for (i = 0; i < dim; i++) {
        ytemp[i] = i0*y0[i] + iend*yinit[i] + hUsed*i1*k1[i] + hUsed*i6*k6[i];
      }

      /* Store the interpolated value in the output array. */
      count++;
      if ((status = storeStateInOutput(&output, tintp, ytemp, dim, &outputlen, count)) == XLAL_ENOMEM) goto bail_out;
    }

    /* Now that we have recorded the last interpolated step that we
       could, check for termination criteria. */
    if (!integrator->stopontestonly && t >= tend) break;

    /* If there is a stopping function in integrator, call it with the
       last value of y and dydt from the integrator. */
    if (integrator->stop) {
      if ((status = integrator->stop(t, yinit, integrator->evolve->dydt_out, params)) != GSL_SUCCESS) {
        integrator->returncode = status;
        break;
      }
    }
  }

  XLAL_ENDGSL;

  /* Now that the interpolation is done, shrink the output array down
     to exactly count samples. */
  shrinkOutput(&output, &outputlen, count, dim);

  /* Store the final *interpolated* sample in yinit. */
  for (i = 0; i < dim; i++) {
    yinit[i] = output->data[(i+2)*outputlen - 1];
  }


 bail_out:
  /* If we have an error, then we should free allocated memory, and
     then return.  Currently, the only errors we encounter in this
     function are XLAL_ENOMEM, so we just test for that. */
  XLALFree(ytemp);

  if (status == XLAL_ENOMEM || status == XLAL_EINVAL) {
    if (output) XLALDestroyREAL8Array(output);
    *yout = NULL;
    XLAL_ERROR(status);
  }

  *yout = output;
  return outputlen;
}

int XLALAdaptiveRungeKutta4( ark4GSLIntegrator *integrator,
                         void *params,
                         REAL8 *yinit,
                         REAL8 tinit, REAL8 tend, REAL8 deltat,
                         REAL8Array **yout
                         )
{
  int status; /* used throughout */

	/* needed for the integration */
  size_t dim, bufferlength, cnt, retries;
	REAL8 t, tnew, h0;
  REAL8Array *buffers = NULL;
  REAL8 *temp = NULL, *y, *y0, *dydt_in, *dydt_in0, *dydt_out, *yerr; /* aliases */

	/* needed for the final interpolation */
  gsl_spline *interp = NULL;
  gsl_interp_accel *accel = NULL;
  int outputlen = 0;
	REAL8Array *output = NULL;
  REAL8 *times, *vector;	/* aliases */

  /* allocate the buffers!
	   note: REAL8Array has a field dimLength (UINT4Vector) with dimensions, and a field data that points to a single memory block;
	         dimLength itself has fields length and data */
  dim = integrator->sys->dimension;
  bufferlength = (int)((tend - tinit)/deltat) + 2;    		/* allow for the initial value and possibly a final semi-step */
  buffers = XLALCreateREAL8ArrayL(2,dim+1,bufferlength);	/* 2-dimensional array, (dim+1) x bufferlength */
  temp = LALCalloc(6*dim,sizeof(REAL8));

  if (!buffers || !temp) {
    status = XLAL_ENOMEM;
    goto bail_out;
  }

	y = temp; y0 = temp + dim; dydt_in = temp + 2*dim; dydt_in0 = temp + 3*dim; dydt_out = temp + 4*dim; yerr = temp + 5*dim; /* aliases */

  /* set up to get started */
  integrator->sys->params = params;

	integrator->returncode = 0;

  cnt = 0;
  retries = integrator->retries;

	t  = tinit; h0 = deltat;
  memcpy(y,yinit,dim*sizeof(REAL8));

  /* store the first data point */
  buffers->data[0] = t;
  for(unsigned int i=1;i<=dim;i++) buffers->data[i*bufferlength] = y[i-1];

  /* compute derivatives at the initial time (dydt_in), bail out if impossible */
  if ((status = integrator->dydt(t,y,dydt_in,params)) != GSL_SUCCESS) {
		integrator->returncode = status;
		goto bail_out;
	}

	/* note: for speed, this replaces the single CALLGSL wrapper applied before each GSL call */
	XLAL_BEGINGSL;

  while(1) {

     if (!integrator->stopontestonly && t >= tend) {
                        break;
     }

		if (integrator->stop) {
      if ((status = integrator->stop(t,y,dydt_in,params)) != GSL_SUCCESS) {
				integrator->returncode = status;
				break;
			}
		}

		/* ready to try stepping! */
    try_step:

		/* if we would be stepping beyond the final time, stop there instead... */
    if (!integrator->stopontestonly && t + h0 > tend) h0 = tend - t;

    memcpy(y0      ,y,      dim*sizeof(REAL8));   /* save y to y0, dydt_in to dydt_in0 */
		memcpy(dydt_in0,dydt_in,dim*sizeof(REAL8));

    /* call the GSL stepper function */
		status = gsl_odeiv_step_apply(integrator->step,t,h0,y,yerr,dydt_in,dydt_out,integrator->sys);
		/* note: If the user-supplied functions defined in the system dydt return a status other than GSL_SUCCESS,
		   the step will be aborted. In this case, the elements of y will be restored to their pre-step values,
		   and the error code from the user-supplied function will be returned. */

		/* did the stepper report a derivative-evaluation error? */
    if (status != GSL_SUCCESS) {
      if (retries--) {
        h0 = h0 / 10.0;               /* if we have singularity retries left, reduce the timestep and try again */
        goto try_step;
      } else {
				integrator->returncode = status;
				break;                        /* otherwise exit the loop */
      }
    } else {
      retries = integrator->retries;  /* we stepped successfully, reset the singularity retries */
    }

    tnew = t + h0;

    /* call the GSL error-checking function */
    status = gsl_odeiv_control_hadjust(integrator->control,integrator->step,y,yerr,dydt_out,&h0);

		/* did the error-checker reduce the stepsize?
		   note: other possible return codes are GSL_ODEIV_HADJ_INC if it was increased,
																						 GSL_ODEIV_HADJ_NIL if it was unchanged */
    if (status == GSL_ODEIV_HADJ_DEC) {
      memcpy(y,      y0,      dim*sizeof(REAL8));	/* if so, undo the step, and try again */
			memcpy(dydt_in,dydt_in0,dim*sizeof(REAL8));
			goto try_step;
    }

		/* update the current time and input derivatives */
    t = tnew;
    memcpy(dydt_in,dydt_out,dim*sizeof(REAL8));
    cnt++;

    /* check if interpolation buffers need to be extended */
    if (cnt >= bufferlength) {
			REAL8Array *rebuffers;

				/* sadly, we cannot use ResizeREAL8Array, because it would only work if we extended the first array dimension,
					 so we need to copy everything over and switch the buffers. Oh well. */
      if (!(rebuffers = XLALCreateREAL8ArrayL(2,dim+1,2*bufferlength))) {
        status = XLAL_ENOMEM;	/* ouch, that hurt */
        goto bail_out;
      } else {
				for(unsigned int i=0;i<=dim;i++) memcpy(&rebuffers->data[i*2*bufferlength],&buffers->data[i*bufferlength],cnt*sizeof(REAL8));
				XLALDestroyREAL8Array(buffers); buffers = rebuffers;
				bufferlength *= 2;
			}
    }

    /* copy time and state into interpolation buffers */
    buffers->data[cnt] = t;
    for(unsigned int i=1;i<=dim;i++) buffers->data[i*bufferlength + cnt] = y[i-1]; /* y does not have time */
  }
	XLAL_ENDGSL;

	/* copy the final state into yinit */

	memcpy(yinit,y,dim*sizeof(REAL8));

  /* if we have completed at least one step, allocate the GSL interpolation object and the output array */
  if (cnt == 0) goto bail_out;

  XLAL_CALLGSL( interp = gsl_spline_alloc(gsl_interp_cspline,cnt+1) );
  XLAL_CALLGSL( accel  = gsl_interp_accel_alloc() );

  outputlen = (int)(t / deltat) + 1;
  output = XLALCreateREAL8ArrayL(2,dim+1,outputlen);

  if (!interp || !accel || !output) {
    status = XLAL_ENOMEM;	/* ouch again, ran out of memory */
	  if (output) XLALDestroyREAL8Array(output);
		outputlen = 0;
    goto bail_out;
  }

  /* make an array of times */
  times = output->data;
  for(int j=0;j<outputlen;j++) times[j] = tinit + deltat * j;

  /* interpolate! */
	XLAL_BEGINGSL;
  for(unsigned int i=1;i<=dim;i++) {
    gsl_spline_init(interp,&buffers->data[0],&buffers->data[bufferlength*i],cnt+1);

    vector = output->data + outputlen*i;
    for(int j=0;j<outputlen;j++) {
      gsl_spline_eval_e(interp,times[j],accel, &(vector[j]));
    }
  }
	XLAL_ENDGSL;

  /* deallocate stuff and return */
  bail_out:

  if (buffers) XLALDestroyREAL8Array(buffers);   /* let's be careful, although all these checks may not be needed */
  if (temp)    LALFree(temp);

  if (interp)  XLAL_CALLGSL( gsl_spline_free(interp) );
  if (accel)   XLAL_CALLGSL( gsl_interp_accel_free(accel) );

  if (status == XLAL_ENOMEM) XLAL_ERROR(XLAL_ENOMEM);	/* TO DO: will this return? */

	*yout = output;
  return outputlen; /* TO DO: check XLAL error reporting conventions */
}





/**
 * Fourth-order Runge-Kutta ODE integrator using Runge-Kutta-Fehlberg (RKF45)
 * steps with adaptive step size control.  Intended for use in Fourier domain
 * waveform generation routines based on SpinTaylorTN models.
 *
 * The method is described in
 *
 * Abramowitz & Stegun, Handbook of Mathematical Functions, Tenth Printing,
 * National Bureau of Standards, Washington, DC, 1972
 * (available online at http://people.math.sfu.ca/~cbm/aands/ )
 *
 * This method is equivalent to XLALAdaptiveRungeKutta4 and
 * XLALAdaptiveRungeKutta4Hermite, but does not includes any interpolation.
 *
 */
int XLALAdaptiveRungeKutta4IrregularIntervals( ark4GSLIntegrator *integrator,      /**< struct holding dydt, stopping test, stepper, etc. */
                                    void *params,                       /**< params struct used to compute dydt and stopping test */
                                    REAL8 *yinit,                       /**< pass in initial values of all variables - overwritten to final values */
                                    REAL8 tinit,                        /**< integration start time */
                                    REAL8 tend_in,                      /**< maximum integration time */
                                    REAL8Array **yout                   /**< array holding the unevenly sampled output */
                                    )
{
  UINT4 MaxRK4Steps = 1000000;
  int status; /* used throughout */
  unsigned int i;

  // if we integrate forward in time, we start at i=0, and if we integrate backwards, we start at i=MaxRK4Steps-1
  unsigned int iStart, iCurrent;
  int di;

  REAL8 tend = tend_in;

        /* needed for the integration */
  size_t dim, bufferlength, cnt, retries;
        REAL8 t, tnew, h0;
  REAL8Array *buffers = NULL;
  REAL8 *temp = NULL, *y, *y0, *dydt_in, *dydt_in0, *dydt_out, *yerr; /* aliases */

  int outputlen = 0;
  REAL8Array *output = NULL;

  /* allocate the buffers!
           note: REAL8Array has a field dimLength (UINT4Vector) with dimensions, and a field data that points to a single memory block;
                 dimLength itself has fields length and data */
  dim = integrator->sys->dimension;
  bufferlength = MaxRK4Steps;
  buffers = XLALCreateREAL8ArrayL(2, dim+2, bufferlength);        /* 2-dimensional array, (dim+2) x bufferlength */
  temp = LALCalloc(6*dim,sizeof(REAL8));

  if (!buffers || !temp) {
    status = XLAL_ENOMEM;
    goto bail_out;
  }

        y = temp; y0 = temp + dim; dydt_in = temp + 2*dim; dydt_in0 = temp + 3*dim; dydt_out = temp + 4*dim; yerr = temp + 5*dim; /* aliases */

  /* set up to get started */
  integrator->sys->params = params;

        integrator->returncode = 0;

  cnt = 0;
  retries = integrator->retries;

  t = tinit;
  if(tend > tinit)
  {
    h0 = 1.;
  }
  else
  {
    h0 = -1.;
  }
  memcpy(y,yinit,dim*sizeof(REAL8));


  // current index starts at iStart, and increases in steps of di. We've reached the memory limit if cnt == bufferlength.
  if(h0 > 0.)
  {
    iStart = 0;
    di = 1;
  }
  else
  {
    iStart = MaxRK4Steps-1;
    di = -1;
  }

  iCurrent = iStart;


  /* store the first data point */
  buffers->data[iCurrent] = t;
  for(i=1;i<=dim;i++) buffers->data[i*bufferlength+iCurrent] = y[i-1];

  /* compute derivatives at the initial time (dydt_in), bail out if impossible */
  if ((status = integrator->dydt(t,y,dydt_in,params)) != GSL_SUCCESS) {
                integrator->returncode = status;
                goto bail_out;
        }

  buffers->data[i*bufferlength + iCurrent] = dydt_in[1]; /* add domega/dt. here i=dim+1 */

        /* note: for speed, this replaces the single CALLGSL wrapper applied before each GSL call */
        XLAL_BEGINGSL;

  while(1) {

     if (!integrator->stopontestonly && t >= tend) {
                        break;
     }

                if (integrator->stop) {
      if ((status = integrator->stop(t,y,dydt_in,params)) != GSL_SUCCESS) {
                                integrator->returncode = status;
                                break;
                        }
                }

                /* ready to try stepping! */
    try_step:

                /* if we would be stepping beyond the final time, stop there instead... */
    if (!integrator->stopontestonly && t + h0 > tend) h0 = tend - t;

    memcpy(y0      ,y,      dim*sizeof(REAL8));   /* save y to y0, dydt_in to dydt_in0 */
                memcpy(dydt_in0,dydt_in,dim*sizeof(REAL8));

    /* call the GSL stepper function */
                status = gsl_odeiv_step_apply(integrator->step,t,h0,y,yerr,dydt_in,dydt_out,integrator->sys);
                /* note: If the user-supplied functions defined in the system dydt return a status other than GSL_SUCCESS,
                   the step will be aborted. In this case, the elements of y will be restored to their pre-step values,
                   and the error code from the user-supplied function will be returned. */

                /* did the stepper report a derivative-evaluation error? */
    if (status != GSL_SUCCESS) {
      if (retries--) {
        h0 = h0 / 10.0;               /* if we have singularity retries left, reduce the timestep and try again */
        goto try_step;
      } else {
                                integrator->returncode = status;
                                break;                        /* otherwise exit the loop */
      }
    } else {
      retries = integrator->retries;  /* we stepped successfully, reset the singularity retries */
    }

    tnew = t + h0;

    /* call the GSL error-checking function */
    status = gsl_odeiv_control_hadjust(integrator->control,integrator->step,y,yerr,dydt_out,&h0);

                /* did the error-checker reduce the stepsize?
                   note: other possible return codes are GSL_ODEIV_HADJ_INC if it was increased, GSL_ODEIV_HADJ_NIL if it was unchanged */
    if (status == GSL_ODEIV_HADJ_DEC) {
      memcpy(y,      y0,      dim*sizeof(REAL8));       /* if so, undo the step, and try again */
                        memcpy(dydt_in,dydt_in0,dim*sizeof(REAL8));
                        goto try_step;
    }

                /* update the current time and input derivatives */
    t = tnew;
    memcpy(dydt_in,dydt_out,dim*sizeof(REAL8));
    cnt++;
    iCurrent += di;

    /* check if interpolation buffers need to be extended */
    if (cnt >= bufferlength) {
                        REAL8Array *rebuffers;

                                /* sadly, we cannot use ResizeREAL8Array, because it would only work if we extended the first array dimension,
                                         so we need to copy everything over and switch the buffers. Oh well. */
      if (!(rebuffers = XLALCreateREAL8ArrayL(2,dim+1,2*bufferlength)))
      {
        status = XLAL_ENOMEM;   /* ouch, that hurt */
        goto bail_out;
      }
      else
      {
        for(i=0; i<=dim; i++) memcpy(&rebuffers->data[i*2*bufferlength],&buffers->data[i*bufferlength],cnt*sizeof(REAL8));
        XLALDestroyREAL8Array(buffers); buffers = rebuffers;
        bufferlength *= 2;
        if(di < 0)
        {
          iCurrent += bufferlength;
        }
      }
    }

    /* copy time and state into interpolation buffers */
    buffers->data[iCurrent] = t;
    for(i=1; i<=dim; i++) buffers->data[i*bufferlength + iCurrent] = y[i-1]; /* y does not have time */

    buffers->data[i*bufferlength + iCurrent] = dydt_in[1]; /* add domega/dt. here i=dim+1 */
  }

        XLAL_ENDGSL;

        /* copy the final state into yinit */

        memcpy(yinit,y,dim*sizeof(REAL8));

  /* if we have completed at least one step, allocate the output array */
  if (cnt == 0) goto bail_out;


  outputlen = cnt+1;
  output = XLALCreateREAL8ArrayL(2,dim+2,outputlen);

  // depending on the direction of integration, we copy starting from 0 or from the current index
  if(di > 0)
  {
    iStart = 0;
  }
  else
  {
    iStart = iCurrent;
  }


  if (!output) {
    status = XLAL_ENOMEM;       /* ouch again, ran out of memory */
          if (output) XLALDestroyREAL8Array(output);
                outputlen = 0;
    goto bail_out;
  }

  for(i = 0; i <= dim+1; i++)
  {
    memcpy(&(output->data[i*outputlen]), &(buffers->data[i*bufferlength + iStart]), outputlen*sizeof(REAL8));
  }

  /* deallocate stuff and return */
  bail_out:

  if (buffers) XLALDestroyREAL8Array(buffers);   /* let's be careful, although all these checks may not be needed */
  if (temp)    LALFree(temp);

  if (status == XLAL_ENOMEM) XLAL_ERROR(XLAL_ENOMEM);   /* TO DO: will this return? */

        *yout = output;
  return outputlen; /* TO DO: check XLAL error reporting conventions */
}
