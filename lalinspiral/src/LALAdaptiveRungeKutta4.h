#ifndef _LALADAPTIVERUNGEKUTTA4_H
#define _LALADAPTIVERUNGEKUTTA4_H

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>

#include <lal/LALGSL.h>
#include <lal/SeqFactories.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

#define XLAL_BEGINGSL \
        { \
          gsl_error_handler_t *saveGSLErrorHandler_; \
          XLALGSL_PTHREAD_MUTEX_LOCK; \
          saveGSLErrorHandler_ = gsl_set_error_handler_off();

#define XLAL_ENDGSL \
          gsl_set_error_handler( saveGSLErrorHandler_ ); \
          XLALGSL_PTHREAD_MUTEX_UNLOCK; \
        }

/* <lalVerbatim file="LALInspiralRungeKuttaH">  */
typedef struct
tagark4GSLIntegrator
{
  gsl_odeiv_step    *step;
  gsl_odeiv_control *control;
  gsl_odeiv_evolve  *evolve;

  gsl_odeiv_system  *sys;   

  int (* dydt) (double t, const double y[], double dydt[], void * params);
  int (* stop) (double t, const double y[], double dydt[], void * params);
    
  int retries;		/* retries with smaller step when derivatives encounter singularity */
	int stopontestonly;	/* stop only on test, use tend to size buffers only */

	int returncode;
} ark4GSLIntegrator;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{ark4GSLIntegrator}}
</lalLaTeX>  */

ark4GSLIntegrator *
XLALAdaptiveRungeKutta4Init( int dim,
									 					 int (* dydt) (double t, const double y[], double dydt[], void * params),
									 					 int (* stop) (double t, const double y[], double dydt[], void * params),
									 					 double eps_abs, double eps_rel
													 );

void
XLALAdaptiveRungeKutta4Free( ark4GSLIntegrator *integrator );

unsigned int
XLALAdaptiveRungeKutta4( ark4GSLIntegrator *integrator,
						 						 void *params,
						 						 REAL8 *yinit,
						 						 REAL8 tinit, REAL8 tend, REAL8 deltat,
						 						 REAL8Array **yout
											 );

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _LALADAPTIVERUNGEKUTTA4_H */
