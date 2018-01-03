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
 * \defgroup LALAdaptiveRungeKutta4_h Header LALAdaptiveRungeKutta4.h
 * \ingroup lal_utilities
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

typedef struct tagLALAdaptiveRungeKutta4Integrator
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
} LALAdaptiveRungeKutta4Integrator;

LALAdaptiveRungeKutta4Integrator *XLALAdaptiveRungeKutta4Init( int dim,
                             int (* dydt) (double t, const double y[], double dydt[], void * params),
                             int (* stop) (double t, const double y[], double dydt[], void * params),
                             double eps_abs, double eps_rel
                             );

/* OPTIMIZED */
LALAdaptiveRungeKutta4Integrator *XLALAdaptiveRungeKutta4InitEighthOrderInstead( int dim,
                             int (* dydt) (double t, const double y[], double dydt[], void * params),
                             int (* stop) (double t, const double y[], double dydt[], void * params),
                             double eps_abs, double eps_rel
                             );
/* END OPTIMIZED */

void XLALAdaptiveRungeKutta4Free( LALAdaptiveRungeKutta4Integrator *integrator );

int XLALAdaptiveRungeKutta4( LALAdaptiveRungeKutta4Integrator *integrator,
                         void *params,
                         REAL8 *yinit,
                         REAL8 tinit, REAL8 tend, REAL8 deltat,
                         REAL8Array **yout
                         );
/* OPTIMIZED */
int XLALAdaptiveRungeKutta4NoInterpolate(LALAdaptiveRungeKutta4Integrator * integrator,
         void * params, REAL8 * yinit, REAL8 tinit, REAL8 tend, REAL8 deltat_or_h0,
         REAL8Array ** t_and_yout);
/* END OPTIMIZED */

int XLALAdaptiveRungeKutta4Hermite( LALAdaptiveRungeKutta4Integrator *integrator,
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
int XLALAdaptiveRungeKutta4IrregularIntervals( LALAdaptiveRungeKutta4Integrator *integrator,      /**< struct holding dydt, stopping test, stepper, etc. */
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
