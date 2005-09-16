/*  <lalVerbatim file="LALRungeKutta4AdaptCV">
Author: Robinson, C. A.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALRungeKutta4Adapt.c}}

The code \texttt{LALRungeKutta4Adapt.c} solves a system of $n$ coupled first--order differential equations.
Internally, it uses the gsl routines for performing adaptive step evolution of the system, but to the outside
user, it returns results for a fixed step size.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALRungeKutta4AdaptCP}
\idx{LALRungeKutta4Adapt()}

\subsubsection*{Description}


\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALRungeKutta4AdaptCV}}

</lalLaTeX>  */


#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>

struct RungeGSLParams {
  rk4In *input;
  void  *params;
};

static int derivativeGSLWrapper(
                                REAL8 t,
                                const REAL8 y[],
                                REAL8 dydx[],
                                void *params);


NRCSID (LALRUNGEKUTTA4ADAPTC, "$Id$");

/*  <lalVerbatim file="LALRungeKutta4AdaptCP"> */
void 
LALRungeKutta4(
   LALStatus   *status,
   REAL8Vector *yout,
   rk4In       *input,
   void        *params
   )
{ /* </lalVerbatim>  */

   INT4 i;
   REAL8 t = 0.0;
   struct RungeGSLParams gslParams;
   REAL8 y[input->n];
   REAL8 h = input->h;
   const gsl_odeiv_step_type * type = gsl_odeiv_step_rk4;
   gsl_odeiv_step *step;
   gsl_odeiv_control *control;
   gsl_odeiv_evolve *evolve;
   gsl_odeiv_system sys; 
  
   INITSTATUS(status, "LALRungeKutta4Adapt", LALRUNGEKUTTA4ADAPTC);
   ATTATCHSTATUSPTR(status);

   ASSERT (yout, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (yout->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (input, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  /* Initialise GSL integrator */
  step    = gsl_odeiv_step_alloc(type, input->n);
  control = gsl_odeiv_control_standard_new(0.0, 1.0e-12, 0.5, 0.5);
  evolve  = gsl_odeiv_evolve_alloc(input->n);



  gslParams.input = input;
  gslParams.params = params;

  sys.function = derivativeGSLWrapper;
  sys.jacobian =  NULL;
  sys.dimension = input->n; 
  sys.params    = &gslParams;


  for (i = 0; i < input->n; i++)
    y[i] = input->y->data[i];

   /* Evolve the system */
  while (t < input->h)
  {
    int gslStatus = gsl_odeiv_evolve_apply(evolve, control, step, &sys,
				&t, input->h, &h, y);
     /*printf("h = %e, t = %e\n", h, t);*/
    if (gslStatus != GSL_SUCCESS)
    {  
  	gsl_odeiv_evolve_free(evolve);
   	gsl_odeiv_control_free(control);
   	gsl_odeiv_step_free(step);
        ABORT(status, LALINSPIRALH_ESTOPPED, LALINSPIRALH_MSGESTOPPED);
     }  
  }

  for (i=0; i<input->n; i++)
     yout->data[i] = y[i];

  gsl_odeiv_evolve_free(evolve);
  gsl_odeiv_control_free(control);
  gsl_odeiv_step_free(step);
   
  DETATCHSTATUSPTR(status);
  RETURN (status);
}


/* A simple wrapper function to allow GSL to use the LAL
   derivative functions */
int derivativeGSLWrapper(
				REAL8 t, 
				const REAL8 y[], 
				REAL8 dydx[],
				void *params)
{
  struct RungeGSLParams *in = (struct RungeGSLParams *)params;
  REAL8Vector yVect;
  REAL8Vector dyVect;

  yVect.length = dyVect.length = in->input->n;
  yVect.data = (REAL8 *)LALMalloc(in->input->n * sizeof(REAL8));
  memcpy(yVect.data, y, in->input->n * sizeof(REAL8));
  dyVect.data = dydx;
  in->input->function(&yVect, &dyVect, in->params);
  LALFree(yVect.data);
  return GSL_SUCCESS;
}
