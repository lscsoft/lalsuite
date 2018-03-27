#ifndef DISTANCE_INTEGRATOR_H
#define DISTANCE_INTEGRATOR_H

/* exclude from SWIG interface and C++ code */
#if !defined(SWIG) && !defined(__cplusplus)

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <lal/cubic_interp.h>

typedef struct taglog_radial_integrator {
		bicubic_interp *region0;
		cubic_interp *region1;
		cubic_interp *region2;
		double xmax, ymax, vmax, r1, r2;
		int k;
} log_radial_integrator;

typedef struct tagradial_integrand_params {
    double scale;
    double p;
    double b;
    int k, cosmology;
} radial_integrand_params;

double log_dVC_dVL(double DL);
void dVC_dVL_init(void);

double log_radial_integrand(double r, void *params);
double log_radial_integral(double r1, double r2, double p, double b, int k, int cosmology);

/**
 * Distance integrator for marginalisation. Assumes a besselI0-type marginalised phase likelihood.
 * @param r1 Minimum distance (Mpc)
 * @param r2 Maximum distance (Mpc)
 * @param k  Exponent of distance prior \f$ p(r) propto r^k \f$
 * @param cosmology 0: Euclidean, 1: use co-moving volume prior
 * @param pmax: The maximum optimal SNR to allow
 * @param size: Size of lookup table
 */
log_radial_integrator *log_radial_integrator_init(double r1, double r2, int k, int cosmology, double pmax, size_t size);

/**
 * Free an integrator
 */
void log_radial_integrator_free(log_radial_integrator *integrator);

/**
 * Evaluate the log distance integrator for given SNRs.
 * With a template at reference distance (1Mpc), compute the marginal likelihood
 * over distance. Uses the two SNRs \f$ p=sqrt(<h|h>) \f$ and \f$ b=<d|h> \f$.
 * @param integrator a log_radial_integrator
 * @param p The optimal SNR \f$ p = sqrt(<h|h>) \f$
 * @param b match between template and data \f$ b = <h|d> \f$
 * @param log_p log(p)
 * @param log_b log(b)
 */
double log_radial_integrator_eval(const log_radial_integrator *integrator, double p, double b, double log_p, double log_b);

#endif /* !defined(SWIG) && !defined(__cplusplus) */

#endif /* DISTANCE_INTEGRATOR_H */
