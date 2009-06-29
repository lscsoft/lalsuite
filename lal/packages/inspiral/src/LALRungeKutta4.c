/*
*  Copyright (C) 2007 Duncan Brown, Gareth Jones, Jolien Creighton, B.S. Sathyaprakash, Craig Robinson
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

/*  <lalVerbatim file="LALRungeKutta4CV">
Author: Robinson, C. A.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALRungeKutta4.c}}
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALRungeKutta4CP}
\idx{LALRungeKutta4()}

\begin{itemize}
\item {\tt n:} The number of coupled equations being integrated.
\item {\tt yout:} The output values for the system after the time-step.
\item {\tt input:} The input for the system
\item {\tt integrator} Required for the GSL integratior. Created using {\tt XLALRungeKutta4Init()}.
\item {\tt params} Parameters to be passed to the derivative function
\end{itemize}

\subsubsection*{Description}
The code \texttt{LALRungeKutta4.c} solves a system of $n$ coupled first--order differential equations.
Internally, it uses the gsl routines for performing adaptive step evolution of the system, but to the outside
user, it returns results for a fixed step size.

Prior to evolving a system using {\tt LALRungeKutta4()}, it is necessary to create the GSL integrator using
{\tt XLALRungeKutta4Init()}. Once the evolution of the system has finished, this integrator should then
be freed using {\tt XLALRungeKutta4Free()}.
\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALRungeKutta4CV}}

</lalLaTeX>  */


#include <lal/LALInspiral.h>

/* macro to "use" unused function parameters */
#define UNUSED(expr) do { (void)(expr); } while (0)

struct RungeGSLParams {
  rk4In *input;
  void  *params;
};

static int derivativeGSLWrapper(
                                REAL8 t,
                                const REAL8 y[],
                                REAL8 dydx[],
                                void *params);



NRCSID (LALRUNGEKUTTA4C, "$Id$");

/* Function for allocating memory and setting up the GSL integrator */
/* <lalVerbatim file="LALRungeKutta4CP"> */
rk4GSLIntegrator * XLALRungeKutta4Init( INT4 n,
                                        rk4In *input
                                      )
{ /* </lalVerbatim>  */

  static const char *func = "XLALRungeKutta4Init";
  rk4GSLIntegrator  *integrator = NULL;

  /* Check we have an input */
  if (!input)
    XLAL_ERROR_NULL(func, XLAL_EFAULT);

  /* Allocate memory for the integrator structure */
  if (!(integrator = (rk4GSLIntegrator *) LALCalloc(1, sizeof(rk4GSLIntegrator))))
  {
    XLAL_ERROR_NULL(func, XLAL_ENOMEM);
  }

  integrator->input = input;

  /* Set the algorithm to 4th-order Runge-Kutta */
  integrator->type = gsl_odeiv_step_rkf45;

  /* Allocate memory for data values */
  if (!(integrator->y = (REAL8 *) LALMalloc(n * sizeof(REAL8))))
  {
    LALFree(integrator);
    XLAL_ERROR_NULL(func, XLAL_ENOMEM);
  }

  /* Initialise GSL integrator */
  XLAL_CALLGSL( integrator->step    = gsl_odeiv_step_alloc(integrator->type, n) );
  XLAL_CALLGSL( integrator->control = gsl_odeiv_control_standard_new(1.0e-2, 1.0e-2, 1.0, 1.0) );
  XLAL_CALLGSL( integrator->evolve  = gsl_odeiv_evolve_alloc(n) );

  /* Check the integrator is allocated correctly */
  if (!(integrator->step) || !(integrator->control) || !(integrator->evolve))
  {
    XLALRungeKutta4Free( integrator );
    XLAL_ERROR_NULL(func, XLAL_ENOMEM);
  }

  return integrator;
}


/*  <lalVerbatim file="LALRungeKutta4CP"> */
void
LALRungeKutta4(
   LALStatus        *status,
   REAL8Vector      *yout,
   rk4GSLIntegrator *integrator,
   void             *params
   )
{ /* </lalVerbatim>  */

   INT4 i;
   REAL8 t = 0.0;
   struct RungeGSLParams gslParams;
   rk4In *input = NULL;
   REAL8 h;
   gsl_odeiv_system sys;

   INITSTATUS(status, "LALRungeKutta4", LALRUNGEKUTTA4C);
   ATTATCHSTATUSPTR(status);

   ASSERT (yout, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (yout->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (integrator, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (integrator->input, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  /* Initialise GSL integrator */

  input = integrator->input;
  h     = input->h;

  gslParams.input = input;
  gslParams.params = params;

  sys.function = derivativeGSLWrapper;
  sys.jacobian =  NULL;
  sys.dimension = input->n;
  sys.params    = &gslParams;


  memcpy( integrator->y, input->y->data, input->n * sizeof(REAL8));

   /* Evolve the system */
  while (t < input->h)
  {
    REAL8 tOld = t;
    CALLGSL( gsl_odeiv_evolve_apply(integrator->evolve, integrator->control,
                 integrator->step, &sys,
				&t, input->h, &h, integrator->y), status );
    /*printf("h = %e, t = %e\n", h, t);*/
    BEGINFAIL(status)
    {
        ABORT(status, LALINSPIRALH_ESTOPPED, LALINSPIRALH_MSGESTOPPED);
    }
    ENDFAIL(status);

    /* In case integration becomes degenerate */
    if (t == tOld)
    {
         for (i=0; i<input->n; i++)
           yout->data[i] = 0.0;

         ABORT(status, LALINSPIRALH_ESTOPPED, LALINSPIRALH_MSGESTOPPED);
    }
  }

  memcpy( yout->data, integrator->y, input->n * sizeof(REAL8));

  DETATCHSTATUSPTR(status);
  RETURN (status);
}


/* Function for freeing up memory for the GSL integrator */
/*  <lalVerbatim file="LALRungeKutta4CP"> */
void XLALRungeKutta4Free( rk4GSLIntegrator *integrator )
{ /* </lalVerbatim> */

  static const char *func = "XLALRungeKutta4Free";

  if (!integrator) XLAL_ERROR_VOID(func, XLAL_EFAULT);

  /* Free the GSL integrator controls etc */
  if (integrator->evolve)  XLAL_CALLGSL( gsl_odeiv_evolve_free(integrator->evolve) );
  if (integrator->control) XLAL_CALLGSL( gsl_odeiv_control_free(integrator->control) );
  if (integrator->step)    XLAL_CALLGSL( gsl_odeiv_step_free(integrator->step) );
  LALFree( integrator->y );
  LALFree( integrator );
  return;
}


/* A simple wrapper function to allow GSL to use the LAL
   derivative functions */
static int derivativeGSLWrapper(
				REAL8 t,
				const REAL8 y[],
				REAL8 dydx[],
				void *params)
{
  struct RungeGSLParams *in = (struct RungeGSLParams *)params;
  REAL8Vector dyVect;
  REAL8Vector *yVect = in->input->yt;

  /* t is unused in this function */
  UNUSED(t);

  memcpy(yVect->data, y, in->input->n * sizeof(REAL8));

  dyVect.length = in->input->n;
  dyVect.data = dydx;
  in->input->function(yVect, &dyVect, in->params);
  return GSL_SUCCESS;
}
