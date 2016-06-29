/*
 * Copyright (C) 2015-2016  Leo Singer
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 */


#include "bayestar_distance.h"

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_erf.h>

#include <chealpix.h>


static void integrals(
    double z,
    double *x2, double *x3, double *x4,
    double *dx2, double *dx3, double *dx4)
{
    const double H = gsl_sf_hazard(- z);
    const double Hp = - H * (z + H);
    const double z2 = gsl_pow_2(z);
    *x2 = z2 + 1 + z * H;
    *x3 = z * (z2 + 3) + (z2 + 2) * H;
    *x4 = z2 * (z2 + 6) + 3 + z * (z2 + 5) * H;
    *dx2 = 2 * z + H + z * Hp;
    *dx3 = 3 * (z2 + 1) + 2 * z * H + (z2 + 2) * Hp;
    *dx4 = 4 * z * (z2 + 3) + (3 * z2 + 5) * H + z * (z2 + 5) * Hp;
}


static void fdf(double z, void *params, double *fval, double *dfval)
{
    const double mean_std = *(double *)params;
    const double target = 1 / gsl_pow_2(mean_std) + 1;
    double x2, x3, x4, dx2, dx3, dx4;
    integrals(z, &x2, &x3, &x4, &dx2, &dx3, &dx4);
    *fval = target * gsl_pow_2(x3) - x4 * x2;
    *dfval = target * 2 * x3 * dx3 - x4 * dx2 - dx4 * x2;
}


static double f(double z, void *params)
{
    double fval, dfval;
    fdf(z, params, &fval, &dfval);
    return fval;
}


static double df(double z, void *params)
{
    double fval, dfval;
    fdf(z, params, &fval, &dfval);
    return dfval;
}


static int solve_z(double mean_std, double *result)
{
    /* Set up variables for tracking progress toward the solution. */
    static const int max_iter = 50;
    int iter = 0;
    double z = mean_std;
    int status;

    /* Set up solver (on stack). */
    const gsl_root_fdfsolver_type *algo = gsl_root_fdfsolver_steffenson;
    char state[algo->size];
    gsl_root_fdfsolver solver = {algo, NULL, 0, state};
    gsl_function_fdf fun = {f, df, fdf, &mean_std};
    gsl_root_fdfsolver_set(&solver, &fun, z);

    do
    {
        const double zold = z;
        status = gsl_root_fdfsolver_iterate(&solver);
        z = gsl_root_fdfsolver_root(&solver);
        status = gsl_root_test_delta (z, zold, 0, GSL_SQRT_DBL_EPSILON);
        iter++;
    } while (status == GSL_CONTINUE && iter < max_iter);

    *result = z;
    return status;
}


int bayestar_distance_moments_to_parameters(
    double mean, double std, double *mu, double *sigma, double *norm)
{
    /* Set up function to solve. */
    double mean_std = mean / std;
    /* Minimum value of (mean/std) for a quadratically weighted
     * normal distribution. The limit of (mean/std) as (mu/sigma) goes to -inf
     * is sqrt(3). We limit (mean/std) to a little bit more than sqrt(3),
     * because as (mu/sigma) becomes more and more negative the normalization
     * has to get very large.
     */
    static const double min_mean_std = M_SQRT3 + 0.1;
    int status;

    if (gsl_finite(mean_std) && mean_std >= min_mean_std)
    {
        double z, x2, x3, x4, dx2, dx3, dx4;
        status = solve_z(mean_std, &z);
        integrals(z, &x2, &x3, &x4, &dx2, &dx3, &dx4);
        *sigma = mean * x2 / x3;
        *mu = *sigma * z;
        *norm = 1 / (gsl_pow_2(*sigma) * x2 * gsl_sf_erf_Q(-z));
    } else {
        status = GSL_SUCCESS;
        *mu = INFINITY;
        *sigma = 1;
        *norm = 0;
    }

    return status;
}


void bayestar_distance_parameters_to_moments(
    double mu, double sigma, double *mean, double *std, double *norm)
{
    if (gsl_finite(mu / sigma))
    {
        const double z = mu / sigma;
        double x2, x3, x4, dx2, dx3, dx4;

        integrals(z, &x2, &x3, &x4, &dx2, &dx3, &dx4);

        *mean = sigma * x3 / x2;
        *std = *mean * sqrt(x4 * x2 / gsl_pow_2(x3) - 1);
        *norm = 1 / (gsl_pow_2(sigma) * x2 * gsl_sf_erf_Q(-z));
    } else {
        *mean = INFINITY;
        *std = 1;
        *norm = 0;
    }
}


double bayestar_volume_render_kernel(
    double x, double y, double max_distance, int axis0, int axis1,
    const double *R, long nside, int nest,
    const double *prob, const double *mu,
    const double *sigma, const double *norm)
{
    /* Determine which axis to integrate over
     * (the one that isn't in the args) */
    int axis2;
    int axes[] = {0, 0, 0};
    axes[axis0] = 1;
    axes[axis1] = 1;
    for (axis2 = 0; axes[axis2]; axis2++)
        ; /* loop body intentionally no-op */

    /* Construct grid in theta, the elevation angle from the
     * spatial origin to the plane of the screen. */

    /* Transverse distance from origin to point on screen */
    double a = sqrt(gsl_pow_2(x) + gsl_pow_2(y));

    /* Maximum value of theta (at edge of screen-aligned cube) */
    double theta_max = atan2(max_distance, a);
    double dtheta = 0.5 * M_PI / nside / 4;

    /* Construct regular grid from -theta_max to +theta_max */
    double ret = 0;
    for (double theta = -theta_max; theta <= theta_max; theta += dtheta)
    {
        /* Differential z = a tan(theta),
         * dz = dz/dtheta dtheta = a tan'(theta) dtheta = a sec^2(theta) dtheta,
         * and dtheta = const */
        double dz_dtheta = a / gsl_pow_2(cos(theta));
        double z = a * tan(theta);
        double xyz[3];
        xyz[axis0] = x;
        xyz[axis1] = y;
        xyz[axis2] = z;

       /* Transform from screen-aligned cube to celestial coordinates before
        * looking up pixel indices. */
        double vec[3];
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, R, 3, xyz, 1, 0, vec, 1);
        long ipix;
        if (nest)
            vec2pix_nest(nside, vec, &ipix);
        else
            vec2pix_ring(nside, vec, &ipix);
        double r = sqrt(gsl_pow_2(x) + gsl_pow_2(y) + gsl_pow_2(z));

        ret += prob[ipix] / sigma[ipix] * exp(-0.5 * gsl_pow_2((r - mu[ipix]) / sigma[ipix])) * norm[ipix] * dz_dtheta;
    }

    return ret / (M_SQRT2 * M_SQRTPI);
}


double bayestar_marginal_distance_distribution(
    double r, long npix,
    const double *prob, const double *mu,
    const double *sigma, const double *norm)
{
    double sum = 0;

    for (long i = 0; i < npix; i ++)
    {
        sum += prob[i] * norm[i] / sigma[i] * exp(-0.5 * gsl_pow_2((r - mu[i]) / sigma[i])) * gsl_pow_2(r);
    }

    return sum / (M_SQRT2 * M_SQRTPI);
}
