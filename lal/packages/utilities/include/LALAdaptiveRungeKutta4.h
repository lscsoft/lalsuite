#ifndef _LALADAPTIVERUNGEKUTTA4_H
#define _LALADAPTIVERUNGEKUTTA4_H

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>

#include <lal/LALGSL.h>
#include <lal/SeqFactories.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/* remove SWIG interface directives */
#if !defined(SWIG) && !defined(SWIGLAL_STRUCT)
#define SWIGLAL_STRUCT(...)
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


typedef struct
tagark4GSLIntegrator
{
  SWIGLAL_STRUCT(ark4GSLIntegrator);
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


ark4GSLIntegrator *XLALAdaptiveRungeKutta4Init( int dim,
                             int (* dydt) (double t, const double y[], double dydt[], void * params),
                             int (* stop) (double t, const double y[], double dydt[], void * params),
                             double eps_abs, double eps_rel
                             );

void XLALAdaptiveRungeKutta4Free( ark4GSLIntegrator *integrator );

int XLALAdaptiveRungeKutta4( ark4GSLIntegrator *integrator,
                         void *params,
                         REAL8 *yinit,
                         REAL8 tinit, REAL8 tend, REAL8 deltat,
                         REAL8Array **yout
                         );
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
int XLALAdaptiveRungeKutta4Hermite( ark4GSLIntegrator *integrator, /**< struct holding dydt, stopping test, stepper, etc. */
                                    void *params, /**< params struct used to compute dydt and stopping test */
                                    REAL8 *yinit, /**< pass in initial values of all variables - overwritten to final values */
                                    REAL8 tinit, /**< integration start time */
                                    REAL8 tend_in, /**< maximum integration time */
                                    REAL8 deltat, /**< step size for evenly sampled output */
                                    REAL8Array **yout /**< array holding the evenly sampled output */
                                    );

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALADAPTIVERUNGEKUTTA4_H */
