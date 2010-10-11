/*
*  Copyright (C) 2010 Michele Vallisneri
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

/*  <lalVerbatim file="LALAdaptiveRungeKutta4CV">
Author: Vallisneri, M.
$Id:$
</lalVerbatim>  */

/* TO DO: XLAL or LAL in these docfile declarations? */

/*  <lalLaTeX>

\subsection{Module \texttt{LALAdaptiveRungeKutta4.c}}
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALAdaptiveRungeKutta4CP}
\idx{XLALAdaptiveRungeKutta4()}

\begin{itemize}
\item {\tt integrator} Integration structure (quasi-class). Created using {\tt XLALAdaptiveRungeKutta4Init()}.
...
\end{itemize}

\subsubsection*{Description}
The code \texttt{LALAdaptiveRungeKutta4.c} evolves a system of $n$ coupled first--order differential equations.
Internally, it uses GSL routines to perform adaptive-step evolution, and then interpolates the resulting
trajectories to a fixed step size.

Prior to evolving a system using {\tt XLALAdaptiveRungeKutta4()}, it is necessary to create an integrator structure using
{\tt XLALAdaptiveRungeKutta4Init()}. Once you are done with the integrator, free it with {\tt XLALAdaptiveRungeKutta4Free()}.
\subsubsection*{Algorithm}
TBF.

\subsubsection*{Uses}
For updated SpinTaylor waveforms.

\subsubsection*{Notes}
None so far...

\vfill{\footnotesize\input{LALAdaptiveRungeKutta4CV}}

</lalLaTeX>  */

#include <LALAdaptiveRungeKutta4.h>

NRCSID (LALADAPTIVERUNGEKUTTA4C, "$Id$");

/* <lalVerbatim file="XLALAdaptiveRungeKutta4CP"> */
ark4GSLIntegrator *
XLALAdaptiveRungeKutta4Init( int dim,
									 					 int (* dydt) (double t, const double y[], double dydt[], void * params),  /* These are XLAL functions! */
									 					 int (* stop) (double t, const double y[], double dydt[], void * params),
									 					 double eps_abs, double eps_rel
                  				 )
{ /* </lalVerbatim> */
  static const char *func = "XLALAdaptiveRungeKutta4Init"; /* TO DO: is this correct XLAL etiquette? */
  
  ark4GSLIntegrator *integrator;
  
	/* allocate our custom integrator structure */
  if (!(integrator = (ark4GSLIntegrator *) LALCalloc(1, sizeof(ark4GSLIntegrator))))
	{
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
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
    XLAL_ERROR_NULL(func, XLAL_ENOMEM);
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

/* <lalVerbatim file="XLALAdaptiveRungeKutta4CP"> */
void
XLALAdaptiveRungeKutta4Free( ark4GSLIntegrator *integrator )
{	/* </lalVerbatim>  */
  if (!integrator) return;

  if (integrator->evolve)  XLAL_CALLGSL( gsl_odeiv_evolve_free(integrator->evolve) );
  if (integrator->control) XLAL_CALLGSL( gsl_odeiv_control_free(integrator->control) );
  if (integrator->step)    XLAL_CALLGSL( gsl_odeiv_step_free(integrator->step) );

  LALFree( integrator->sys );
  LALFree( integrator );

  return;
}

/* <lalVerbatim file="XLALAdaptiveRungeKutta4CP"> */
unsigned int
XLALAdaptiveRungeKutta4( ark4GSLIntegrator *integrator,
						 						 void *params,
						 						 REAL8 *yinit,
						 						 REAL8 tinit, REAL8 tend, REAL8 deltat,
						 						 REAL8Array **yout
					   					 )
{ /* </lalVerbatim> */
  static const char *func = "XLALAdaptiveRungeKutta4"; /* TO DO: is this correct XLAL etiquette? */
  
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
		if (integrator->stop) {
      if ((status = integrator->stop(t,y,dydt_in,params)) != GSL_SUCCESS) {
				integrator->returncode = status;
				break;
			}
    } else if (!integrator->stopontestonly && t >= tend) {
			break;
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
				for(unsigned int i=0;i<=dim;i++) memcpy(&rebuffers->data[i*2*bufferlength],&buffers->data[i*bufferlength],(cnt-1)*sizeof(REAL8));
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

  XLAL_CALLGSL( interp = gsl_spline_alloc(gsl_interp_cspline,cnt) );
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
    gsl_spline_init(interp,&buffers->data[0],&buffers->data[bufferlength*i],cnt);

    vector = output->data + outputlen*i;
    for(int j=0;j<outputlen;j++) {
      vector[j] = gsl_spline_eval(interp,times[j],accel);
    }
  }
	XLAL_ENDGSL;

  /* deallocate stuff and return */
  bail_out:
  
  if (buffers) XLALDestroyREAL8Array(buffers);   /* let's be careful, although all these checks may not be needed */
  if (temp)    LALFree(temp);

  if (interp)  XLAL_CALLGSL( gsl_spline_free(interp) );
  if (accel)   XLAL_CALLGSL( gsl_interp_accel_free(accel) );
  
  if (status == XLAL_ENOMEM) XLAL_ERROR(func, XLAL_ENOMEM);	/* TO DO: will this return? */

	*yout = output;
  return outputlen; /* TO DO: check XLAL error reporting conventions */
}

