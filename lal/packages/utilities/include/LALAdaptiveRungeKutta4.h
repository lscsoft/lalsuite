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

/**
 * \addtogroup LALAdaptiveRungeKutta4_h
 * \author Vallisneri, M.
 * \brief Adaptive Runge-Kutta4
 *
 * <ul>
 * <li> \c integrator Integration structure (quasi-class). Created using <tt>XLALAdaptiveRungeKutta4Init()</tt>.
 * ...</li>
 * </ul>
 *
 * ### Description ###
 *
 * The code \ref LALAdaptiveRungeKutta4.c evolves a system of \f$n\f$ coupled first--order differential equations.
 * Internally, it uses GSL routines to perform adaptive-step evolution, and then interpolates the resulting
 * trajectories to a fixed step size.
 *
 * Prior to evolving a system using <tt>XLALAdaptiveRungeKutta4()</tt>, it is necessary to create an integrator structure using
 * <tt>XLALAdaptiveRungeKutta4Init()</tt>. Once you are done with the integrator, free it with <tt>XLALAdaptiveRungeKutta4Free()</tt>.
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
int XLALAdaptiveRungeKutta4Hermite( ark4GSLIntegrator *integrator,
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
int XLALAdaptiveRungeKutta4IrregularIntervals( ark4GSLIntegrator *integrator,      /**< struct holding dydt, stopping test, stepper, etc. */
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

#endif /* _LALADAPTIVERUNGEKUTTA4_H */
