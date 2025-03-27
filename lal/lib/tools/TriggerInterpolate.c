/*
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
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA
 *
 * Copyright (C) 2012-2020 Leo Singer
 */

#include <complex.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_version.h>

#include <lal/TriggerInterpolate.h>

#define STRINGIZE(s) #s


/*
 * General functions
 */


/* Unwrap phase angles. All elements of arg must be between -M_PI and M_PI. */
static void unwrap(double *arg, size_t n)
{
    double prev = arg[0];
    long winds = 0;
    for (size_t i = 1; i < n; i ++)
    {
        double next = arg[i], delta = next - prev;

        if (delta > M_PI)
            winds -= 1;
        else if (delta < -M_PI)
            winds += 1;
        arg[i] += winds * 2 * M_PI;

        prev = next;
    }
}


static void cabs_carg_unwrapped(double *amp, double *arg, const double complex *y, size_t n)
{
    for (size_t i = 0; i < n; i ++)
    {
        amp[i] = cabs(y[i]);
        arg[i] = carg(y[i]);
    }

    unwrap(arg, n);
}


/* Strip leading zero coefficients of a polynomial.
 * Pass the number of coefficients as n. The new number of coefficients is returned. */
static size_t poly_strip(const double *a, size_t n)
{
    for (; n > 0 && a[n - 1] == 0; n--)
        ; /* loop body intentionally a no-op */
    return n;
}


/* Polynomial multiply-accumulate.
 * The resulting polynomial has length (na + nb - 1). */
static void poly_mac(double *y, const double *a, size_t na, const double *b, size_t nb)
{
    /* Compute convolution. */
    for (size_t ia = 0; ia < na; ia++)
        for (size_t ib = 0; ib < nb; ib++)
            y[ia + ib] += a[ia] * b[ib];
}


/* Compute derivative of a polynomial.
 * The resulting polynomial has length (n - 1). */
static void poly_der(double *y, const double *a, size_t n)
{
    for (size_t i = 0; i < n - 1; i ++)
        y[i] = (i + 1) * a[i + 1];
}


/* Construct Catmull-Rom inteprolating polynomial. */
static void poly_interp(double *a, const double y[4])
{
    a[3] = 0.5 * (y[3] - y[0]) + 1.5 * (y[1] - y[2]);
    a[2] = y[0] - 2.5 * y[1] + 2 * y[2] - 0.5 * y[3];
    a[1] = 0.5 * (y[2] - y[0]);
    a[0] = y[1];
}


/* Find all real roots of a polynomial.
 * For quartic and higher order polynomials, use gsl_poly_complex_solve,
 * but allocate the workspace on the stack.
 *
 * Since the workspace for a polynomial of n coefficients is an (n-1)x(n-1)
 * matrix, to avoid exhausting stack memory, this function should only be
 * used for low-order (n<100) polynomials.
 */
static int poly_real_roots(size_t *nroots, double *r, const double *a, size_t n) {
    n = poly_strip(a, n);

    if (n <= 3) {

       /* Handle polynomials of quadratic or lower order.
        * gsl_poly_solve_quadratic itself handles n=0, n=1, and n=2
        * as special cases. */
        *nroots = gsl_poly_solve_quadratic(
            a[2], a[1], a[0],
            &r[0], &r[1]);

        return GSL_SUCCESS;

    } else if (n <= 4) {

        /* Handle cubic polynomials. */
        *nroots = gsl_poly_solve_cubic(
            a[2] / a[3], a[1] / a[3], a[0] / a[3],
            &r[0], &r[1], &r[2]);

        return GSL_SUCCESS;

    } else {

        /* Handle polynomials of quartic or higher order. */
        double complex complex_roots[n - 1];
        double workspace_matrix[n - 1][n - 1];
        gsl_poly_complex_workspace workspace = {n - 1, *workspace_matrix};

        int result = gsl_poly_complex_solve(a, n, &workspace,
            (gsl_complex_packed_ptr) complex_roots);

        if (result == GSL_SUCCESS)
        {
            *nroots = 0;
            for (size_t i = 0; i + 1 < n; i ++)
                if (cimag(complex_roots[i]) == 0)
                    r[(*nroots)++] = creal(complex_roots[i]);
        }
        return result;

    }
}


/* Find least-squares fit of a polynomial to data.
 * For quadratic and higher order polynomials, workspaces are allocated
 * on the stack, so this function should not be used for large problems
 * (e.g., polynomial order times number of data points much greater than 10k).
 */
static int poly_fit(double *a, size_t na, const double *x, const double *y, size_t n)
{
    if (na > n)
    {
        GSL_ERROR("Requested polynomial with more coefficients than data points", GSL_EINVAL);
    } else if (na == 0) {
        GSL_ERROR("Requested a polynomial with zero coefficients", GSL_EINVAL);
    } else if (na == 1) {
        *a = gsl_stats_mean(y, 1, n);
        return GSL_SUCCESS;
    } else if (na == 2) {
        double unused;
        return gsl_fit_linear(x, 1, y, 1, n, &a[0], &a[1], &unused, &unused, &unused, &unused);
    } else {
        /* Set up solver using stack memory. */
        double X[n][na], A[n][na], Q[na][na], QSI[na][na], C[na][na], S[na], t[n], xt[na], D[na], unused;
        gsl_matrix_view Xview = gsl_matrix_view_array(*X, n, na);
        gsl_matrix_view Aview = gsl_matrix_view_array(*A, n, na);
        gsl_matrix_view Qview = gsl_matrix_view_array(*Q, na, na);
        gsl_matrix_view QSIview = gsl_matrix_view_array(*QSI, na, na);
        gsl_matrix_view Cview = gsl_matrix_view_array(*C, na, na);
        gsl_vector_view Sview = gsl_vector_view_array(S, na);
        gsl_vector_view tview = gsl_vector_view_array(t, n);
        gsl_vector_view xtview = gsl_vector_view_array(xt, na);
        gsl_vector_view Dview = gsl_vector_view_array(D, na);
        gsl_vector_const_view yview = gsl_vector_const_view_array(y, n);
        gsl_vector_view cview = gsl_vector_view_array(a, na);

        gsl_multifit_linear_workspace workspace = {
            .A = &Aview.matrix, .Q = &Qview.matrix, .QSI = &QSIview.matrix, .S
            = &Sview.vector, .t = &tview.vector, .xt = &xtview.vector, .D =
            &Dview.vector};

        #if GSL_MAJOR_VERSION < 2
            workspace.n = n;
            workspace.p = na;
        #else
            workspace.nmax = n;
            workspace.pmax = na;
        #endif

        for (size_t i = 0; i < n; i ++)
            for (size_t ia = 0; ia < na; ia ++)
                X[i][ia] = gsl_pow_int(x[i], ia);

        return gsl_multifit_linear(&Xview.matrix, &yview.vector,
            &cview.vector, &Cview.matrix, &unused, &workspace);
    }
}


/*
 * CubicSpline
 */


#define CUBIC_SPLINE_WINDOW_SIZE 2


int XLALTriggerInterpolateCubicSpline(
    double *t,
    double complex *y,
    const double complex *data,
    unsigned int window)
{
    if (window < CUBIC_SPLINE_WINDOW_SIZE)
        GSL_ERROR(
            "Window size must be >= " STRINGIZE(CUBIC_SPLINE_WINDOW_SIZE),
            GSL_EINVAL);

    int result;

    /* Find which among samples -1, 0, +1 has the maximum amplitude. */
    double complex y_max = data[-1];
    double amp_max = cabs(y_max);
    double t_max = -1;
    for (int i = 0; i < 2; i ++)
    {
        double amp = cabs(data[i]);
        if (amp > amp_max)
        {
            t_max = i;
            amp_max = amp;
            y_max = data[i];
        }
    }

    /* Fit interpolant to samples -2 through 1 and samples -1 through 2. */
    for (int i = -1; i < 1; i ++)
    {
        double polys[2][4], deriv[6] = {0}, roots[5];
        size_t nroots;

        /* Take derivative of absolute value of interpolant.
         * If the interpolant is z(t) = a(t) + i * b(t),
         * then d(|z(t)|^2)/dt = 2 a(t) a'(t) + 2 b(t) b'(t).
         * The factor of 2 is dropped because we are only interested in finding
         * the zeros. */
        for (int j = 0; j < 2; j ++)
        {
            double data_part[4], deriv_part[3];
            for (int k = 0; k < 4; k ++)
                data_part[k] = ((const double *) data)[(k + i - 1) * 2 + j];
            poly_interp(polys[j], data_part);
            poly_der(deriv_part, polys[j], 4);
            poly_mac(deriv, polys[j], 4, deriv_part, 3);
       }

       result = poly_real_roots(&nroots, roots, deriv, 6);
       if (result != GSL_SUCCESS)
           return result;

       for (size_t j = 0; j < nroots; j ++)
       {
           if (roots[j] > 0 && roots[j] < 1)
           {
               double complex y_new = crect(
                   gsl_poly_eval(polys[0], 4, roots[j]),
                   gsl_poly_eval(polys[1], 4, roots[j]));
               double amp = cabs(y_new);
               if (amp > amp_max)
               {
                   t_max = roots[j] + i;
                   y_max = y_new;
                   amp_max = amp;
               }
           }
       }
    }

    *t = t_max;
    *y = y_max;
    return result;
}


/*
 * CubicSplineAmpPhase
 */


int XLALTriggerInterpolateCubicSplineAmpPhase(
    double *t,
    double complex *y,
    const double complex *data,
    unsigned int window)
{
    if (window < CUBIC_SPLINE_WINDOW_SIZE)
        GSL_ERROR(
            "Window size must be >= " STRINGIZE(CUBIC_SPLINE_WINDOW_SIZE),
            GSL_EINVAL);

    double amps[2 * CUBIC_SPLINE_WINDOW_SIZE + 1];
    double args[2 * CUBIC_SPLINE_WINDOW_SIZE + 1];
    double t_max, amp_max, arg_max;
    cabs_carg_unwrapped(
        amps, args,
        data - CUBIC_SPLINE_WINDOW_SIZE,
        2 * CUBIC_SPLINE_WINDOW_SIZE + 1);

    /* Find which among samples -1, 0, +1 has the maximum amplitude. */
    t_max = -1;
    amp_max = amps[1];
    arg_max = args[1];
    for (int i = 0; i < 2; i ++)
    {
        if (amps[i + 2] > amp_max)
        {
            t_max = i;
            amp_max = amps[i + 2];
            arg_max = args[i + 2];
        }
    }

    /* Fit interpolant to samples -2 through 1 and samples -1 through 2. */
    for (int i = -1; i < 1; i ++)
    {
        double amp_poly[4], arg_poly[4], der[3], roots[2];
        size_t nroots;
        int result;

        poly_interp(amp_poly, &amps[i + 1]);
        poly_interp(arg_poly, &args[i + 1]);
        poly_der(der, amp_poly, 4);

        result = poly_real_roots(&nroots, roots, der, 3);
        if (result != GSL_SUCCESS)
            return result;

        for (size_t j = 0; j < nroots; j ++)
        {
            double amp, arg;
            if (roots[j] < 0 || roots[j] > 1)
                continue;

            amp = gsl_poly_eval(amp_poly, 4, roots[j]);
            arg = gsl_poly_eval(arg_poly, 4, roots[j]);
            if (amp > amp_max)
            {
                amp_max = amp;
                arg_max = arg;
                t_max = roots[j] + i;
            }
        }
    }

    *t = t_max;
    *y = cpolar(amp_max, arg_max);
    return GSL_SUCCESS;
}


/*
 * Lanczos
 */


/* Data structure providing arguments for minimizer cost function. */
typedef struct {
    const double complex *data;
    unsigned int window;
} lanczos_params;


/* The Lanczos kernel. */
static double lanczos(double t, double a)
{
    if (t < -a || t > a)
        return 0;

    return gsl_sf_sinc(t) * gsl_sf_sinc(t / a);
}


/* The Lanczos reconstruction filter interpolant. */
static double complex lanczos_interpolant(double t, const lanczos_params *params)
{
    double complex ret = 0;
    for (int i = -(int)params->window; i <= (int)params->window; i ++)
        ret += lanczos(t - i, params->window) * params->data[i];
    return ret;
}


/* The cost function to minimize. */
static double lanczos_cost(double t, void *params)
{
    return -cabs(lanczos_interpolant(t, params));
}


int XLALTriggerInterpolateLanczos(
    double *tmax,
    double complex *ymax,
    const double complex *data,
    unsigned int window)
{
    if (window < 1)
        GSL_ERROR("Window size must be >= 1", GSL_EINVAL);

    /* Set up minimizer state on stack. */
    const gsl_min_fminimizer_type *fminimizer_type = gsl_min_fminimizer_brent;
    char fminimizer_state[fminimizer_type->size];
    gsl_min_fminimizer fminimizer = {.type = fminimizer_type, .state = fminimizer_state};

    lanczos_params params = {data, window};
    gsl_function func = {lanczos_cost, &params};
    int result;

    result = gsl_min_fminimizer_set_with_values(&fminimizer, &func,
        0, -cabs(data[0]), -1, -cabs(data[-1]), 1, -cabs(data[1]));
    if (result != GSL_SUCCESS)
        GSL_ERROR("failed to initialize minimizer", result);

    do {
        result = gsl_min_fminimizer_iterate(&fminimizer);
        if (result != GSL_SUCCESS)
            GSL_ERROR("failed to perform minimizer iteration", result);

        double t1 = gsl_min_fminimizer_x_lower(&fminimizer);
        double t2 = gsl_min_fminimizer_x_upper(&fminimizer);
        result = gsl_min_test_interval(t1, t2, 1e-6, 0);
    } while (result == GSL_CONTINUE);

    if (result != GSL_SUCCESS)
        GSL_ERROR("failed to perform minimizer convergence test", result);

    *tmax = gsl_min_fminimizer_x_minimum(&fminimizer);
    *ymax = lanczos_interpolant(*tmax, &params);
    return result;
}


/*
 * Nearest neighbor
 */


int XLALTriggerInterpolateNearestNeighbor(
    double *tmax,
    double complex *ymax,
    const double complex *data,
    __attribute__ ((unused)) unsigned int window)
{
    *tmax = 0;
    *ymax = *data;
    return GSL_SUCCESS;
}


/*
 * Quadratic fit
 */


int XLALTriggerInterpolateQuadraticFit(
    double *tmax,
    double complex *ymax,
    const double complex *data,
    unsigned int window)
{
    if (window < 1)
        GSL_ERROR("Window size must be >= 1", GSL_EINVAL);

    int result, i, nsamples = 2 * window + 1;
    double t, t_new, amp, arg, x[nsamples], amps[nsamples], args[nsamples], amp_poly[3], arg_poly[2];

    cabs_carg_unwrapped(amps, args, data - (int)window, nsamples);

    i = (amps[window - 1] > amps[window + 1]) ? -1 : 1;
    t = i;
    amp = amps[window + i];
    arg = args[window + i];

    for (i = 0; i < nsamples; i ++)
        x[i] = i - (int)window;

    result = poly_fit(amp_poly, 3, x, amps, nsamples);
    if (result != GSL_SUCCESS)
        return result;

    result = poly_fit(arg_poly, 2, x, args, nsamples);
    if (result != GSL_SUCCESS)
        return result;

    if (amp_poly[2] == 0)
        t_new = GSL_SIGN(amp_poly[1]);
    else
        t_new = -0.5 * amp_poly[1] / amp_poly[2];

    if (t_new > -1 && t_new < 1)
    {
        double amp_new = gsl_poly_eval(amp_poly, 3, t_new);
        double arg_new = gsl_poly_eval(arg_poly, 2, t_new);
        if (amp_new > amp)
        {
            t = t_new;
            amp = amp_new;
            arg = arg_new;
        }
    }

    *tmax = t;
    *ymax = cpolar(amp, arg);
    return result;
}
