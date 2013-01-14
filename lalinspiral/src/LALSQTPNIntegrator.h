#ifndef LALSQTPNINTEGRATOR_H
#define LALSQTPNINTEGRATOR_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>

#include <lal/LALInspiral.h>
#include <lal/LALGSL.h>
#include <lal/LALStatusMacros.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup LALSQTPNIntegrator_h Header LALSQTPNIntegrator.h
 * \ingroup pkg_CBC_NEW
 *
 * \brief Contains the function declarations and structures needed by the integration method.
 * \author László Veréb
 * \date 2010.05.21.
 */
/*@{*/

/** The structure contains the integration method and its settings.
 */
typedef struct tagLALSQTPNIntegratorSystem{
	const gsl_odeiv_step_type* type;
	gsl_odeiv_step* step;
	gsl_odeiv_control* control;
	gsl_odeiv_evolve* evolve;
	gsl_odeiv_system system;
} LALSQTPNIntegratorSystem;

/**		The function initialize the integration method.
 * @param[out]	integrator	: the structure containing the integration method
 * @param[in]		num		: the number of the dynamic variables
 * @param[in]		params	: the parameters used in the derivative function
 * @param[in]	derivator	: pointer to the derivative function
 */
int XLALSQTPNIntegratorInit(LALSQTPNIntegratorSystem *integrator, INT2 num, void *params,
		int(*derivator)(REAL8, const REAL8[], REAL8[], void *));

/**		The function deallocates the memory allocated for the integrator
 * function.
 * @param[in]	integrator	: the structure containing the integration method
 */
void XLALSQTPNIntegratorFree(LALSQTPNIntegratorSystem *integrator);

/**		The function evolves the system with the given time-step.
 * @param[in,out]	values	: as input parameters the system's actual position,
 * as ouput the system's next position.
 * @param[in]	integrator	: the integration method
 * @param[in]		step	: the step size
 */
int XLALSQTPNIntegratorFunc(REAL8 values[], LALSQTPNIntegratorSystem *integrator, REAL8 step);

  /*@}*//* end: LALSQTPNIntegrator_h */


#ifdef __cplusplus
}
#endif

#endif /* LALSQTPNINTEGRATOR_H */
