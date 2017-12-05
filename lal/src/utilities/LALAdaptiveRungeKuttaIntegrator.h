#ifndef _LALADAPTIVERUNGEKUTTAINTEGRATOR_H
#define _LALADAPTIVERUNGEKUTTAINTEGRATOR_H

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>

#include <lal/LALGSL.h>
#include <lal/SeqFactories.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \defgroup LALAdaptiveRungeKuttaIntegrator_h Header LALAdaptiveRungeKuttaIntegrator.h
 * \ingroup lal_utilities
 * \author Vallisneri, M.
 * \brief Adaptive Runge-Kutta4
 *
 * <ul>
 * <li> \c integrator Integration structure (quasi-class). Created using <tt>XLALAdaptiveRungeKuttaIntegratorInit()</tt>.
 * ...</li>
 * </ul>
 *
 * ### Description ###
 *
 * The code \ref LALAdaptiveRungeKuttaIntegrator.c evolves a system of \f$n\f$ coupled first--order differential equations.
 * Internally, it uses GSL routines to perform adaptive-step evolution, and then interpolates the resulting
 * trajectories to a fixed step size.
 *
 * Prior to evolving a system using <tt>XLALAdaptiveRungeKutta4()</tt>, it is necessary to create an integrator structure using
 * <tt>XLALAdaptiveRungeKuttaIntegratorInit()</tt>. Once you are done with the integrator, free it with <tt>XLALAdaptiveRungeKuttaIntegratorFree()</tt>.
 *
 * ### Algorithm ###
 *
 * TBF.
 *
 * ### Uses ###
 *
 * For updated SpinTaylor waveforms.
 *
 * ### Notes ###
 *
 * None so far...
 *
 */
/*@{*/

typedef struct tagLALAdaptiveRungeKuttaIntegrator
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
} LALAdaptiveRungeKuttaIntegrator;

LALAdaptiveRungeKuttaIntegrator *XLALAdaptiveRungeKutta4Init( int dim,
                             int (* dydt) (double t, const double y[], double dydt[], void * params),
                             int (* stop) (double t, const double y[], double dydt[], void * params),
                             double eps_abs, double eps_rel
                             );

/* OPTIMIZED */
/**
 * Eighth-order Runge-Kutta ODE integrator using Runge-Kutta-Fehlberg steps
 * with adaptive step size control.  Intended for use in time domain
 * waveform generation routines based on SEOBNRv2,3,4 models.
 */

LALAdaptiveRungeKuttaIntegrator *XLALAdaptiveRungeKutta4InitEighthOrderInstead( int dim,
                             int (* dydt) (double t, const double y[], double dydt[], void * params),
                             int (* stop) (double t, const double y[], double dydt[], void * params),
                             double eps_abs, double eps_rel
                             );
/* END OPTIMIZED */

void XLALAdaptiveRungeKuttaFree( LALAdaptiveRungeKuttaIntegrator *integrator );

int XLALAdaptiveRungeKutta4( LALAdaptiveRungeKuttaIntegrator *integrator,
                         void *params,
                         REAL8 *yinit,
                         REAL8 tinit, REAL8 tend, REAL8 deltat,
                         REAL8Array **yout
                         );
/* OPTIMIZED */
/**
 * Fourth-order Runge-Kutta ODE integrator using Runge-Kutta-Fehlberg steps
 * with adaptive step size control.  Intended for use in time domain
 * waveform generation routines based on SEOBNRv2,3,4 models.  This method
 * does not includes any interpolation.
 */
int XLALAdaptiveRungeKutta4NoInterpolate(LALAdaptiveRungeKuttaIntegrator * integrator,
         void * params, REAL8 * yinit, REAL8 tinit, REAL8 tend, REAL8 deltat_or_h0,
					 REAL8Array ** t_and_yout,INT4 EOBversion);
/* END OPTIMIZED */

int XLALAdaptiveRungeKutta4Hermite( LALAdaptiveRungeKuttaIntegrator *integrator,
                                    void *params,
                                    REAL8 *yinit,
                                    REAL8 tinit,
                                    REAL8 tend_in,
                                    REAL8 deltat,
                                    REAL8Array **yout
                                    );

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
 * Memory is allocated in steps of LAL_MAX_RK4_STEPS.
 */
int XLALAdaptiveRungeKutta4IrregularIntervals( LALAdaptiveRungeKuttaIntegrator *integrator,      /**< struct holding dydt, stopping test, stepper, etc. */
                                    void *params,                       /**< params struct used to compute dydt and stopping test */
                                    REAL8 *yinit,                       /**< pass in initial values of all variables - overwritten to final values */
                                    REAL8 tinit,                        /**< integration start time */
                                    REAL8 tend_in,                      /**< maximum integration time */
                                    REAL8Array **yout                   /**< array holding the unevenly sampled output */
                                    );

/*@}*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALADAPTIVERUNGEKUTTAINTEGRATOR_H */
