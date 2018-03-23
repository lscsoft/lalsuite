#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include "cubic_interp.h"

typedef struct {
		bicubic_interp *region0;
		cubic_interp *region1;
		cubic_interp *region2;
		double xmax, ymax, vmax, r1, r2;
		int k;
} log_radial_integrator;

log_radial_integrator *log_radial_integrator_init(double r1, double r2, int k, int cosmology, double pmax, size_t size);
void log_radial_integrator_free(log_radial_integrator *integrator);
double log_radial_integrator_eval(const log_radial_integrator *integrator, double p, double b, double log_p, double log_b);

