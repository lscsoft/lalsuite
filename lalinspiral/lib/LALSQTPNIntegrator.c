/**
 * @file LALSQTPNIntegrator.c
 * Contains the function definitions needed by the integration method.
 * @author László Veréb
 * @date 2010.05.21.
 */

#include <lal/LALSQTPNIntegrator.h>
#include <lal/LALSQTPNWaveform.h>

int XLALSQTPNIntegratorInit(LALSQTPNIntegratorSystem *integrator, INT2 num, 
		void *params, int(*derivator)(REAL8, const REAL8[], REAL8[], 
		void *)) {
	
	// Check for input errors
	if (num <= 0) {
		XLAL_ERROR(XLAL_EBADLEN);
	}
	if (!params) {
		XLAL_ERROR(XLAL_EFAULT);
	}
	if (!derivator) {
		XLAL_ERROR(XLAL_EFAULT);
	}

	// Initialise GSL integrator
	integrator->type = gsl_odeiv_step_rkf45;
	integrator->system.jacobian = NULL;
	integrator->system.dimension = num;
	integrator->system.params = params;
	integrator->system.function = derivator;
	XLAL_CALLGSL(integrator->step = gsl_odeiv_step_alloc(integrator->type, 
		num));
	XLAL_CALLGSL(integrator->control = gsl_odeiv_control_standard_new(
		1.0e-2, 1.0e-2, 1., 1.));
	XLAL_CALLGSL(integrator->evolve = gsl_odeiv_evolve_alloc(num));

	// Check if the integrator is correctly allocated
	if (!(integrator->step) || !(integrator->control) 
			|| !(integrator->evolve)) {
		XLALSQTPNIntegratorFree(integrator);
		XLAL_ERROR(XLAL_ENOMEM);
	}
	return XLAL_SUCCESS;
}

void XLALSQTPNIntegratorFree(LALSQTPNIntegratorSystem *integrator) {
	if (integrator->evolve) {
		XLAL_CALLGSL(gsl_odeiv_evolve_free(integrator->evolve));
	}
	if (integrator->control) {
		XLAL_CALLGSL(gsl_odeiv_control_free(integrator->control));
	}
	if (integrator->step) {
		XLAL_CALLGSL(gsl_odeiv_step_free(integrator->step));
	}
}

int XLALSQTPNIntegratorFunc(REAL8 values[], 
		LALSQTPNIntegratorSystem *integrator, REAL8 step) {
	REAL8 time = 0., time_Old, step_X = step;
	while (time < step) {
		time_Old = time;
		XLAL_CALLGSL(gsl_odeiv_evolve_apply(integrator->evolve,
				integrator->control, integrator->step,
				&(integrator->system), &time, step, &step_X, 
				values));
		if (time == time_Old) {
			memset(values, 0, 
				integrator->system.dimension * sizeof(REAL8));
			XLAL_ERROR(XLAL_EFUNC);
		}
	}
	return XLAL_SUCCESS;
}
