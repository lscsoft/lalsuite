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
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 *
 * Copyright (C) 2012 Leo Singer
 */

#include <complex.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_trig.h>

#include <lal/TriggerInterpolation.h>


/*
 * Helpers for declaring apply functions for other data types
 */


typedef int (*XLALCOMPLEX16ApplyFunc)(void *, double *, COMPLEX16 *, const COMPLEX16 *);


static int XLALCOMPLEX8ApplyTriggerInterpolant(
    void *interp,
    XLALCOMPLEX16ApplyFunc applyfunc,
    int window,
    double *tmax,
    COMPLEX8 *ymax,
    const COMPLEX8 *y)
{
    int i;
    int ret;
    double tmax_full;
    COMPLEX16 ymax_full;
    COMPLEX16 data_full[2 * window + 1];
    COMPLEX16 *y_full = &data_full[window];

    for (i = -window; i <= window; i ++)
        y_full[i] = y[i];

    ret = applyfunc(interp, &tmax_full, &ymax_full, y_full);
    if (ret == GSL_SUCCESS)
    {
        *tmax = tmax_full;
        *ymax = ymax_full;
    }
    return ret;
}


static int XLALREAL8ApplyTriggerInterpolant(
    void *interp,
    XLALCOMPLEX16ApplyFunc applyfunc,
    int window,
    double *tmax,
    REAL8 *ymax,
    const REAL8 *y)
{
    int i;
    int ret;
    double tmax_full;
    COMPLEX16 ymax_full;
    COMPLEX16 data_full[2 * window + 1];
    COMPLEX16 *y_full = &data_full[window];

    for (i = -window; i <= window; i ++)
        y_full[i] = y[i];

    ret = applyfunc(interp, &tmax_full, &ymax_full, y_full);
    if (ret == GSL_SUCCESS)
    {
        *tmax = tmax_full;
        *ymax = creal(ymax_full);
    }
    return ret;
}


static int XLALREAL4ApplyTriggerInterpolant(
    void *interp,
    XLALCOMPLEX16ApplyFunc applyfunc,
    int window,
    double *tmax,
    REAL4 *ymax,
    const REAL4 *y)
{
    int i;
    int ret;
    double tmax_full;
    COMPLEX16 ymax_full;
    COMPLEX16 data_full[2 * window + 1];
    COMPLEX16 *y_full = &data_full[window];

    for (i = -window; i <= window; i ++)
        y_full[i] = y[i];

    ret = applyfunc(interp, &tmax_full, &ymax_full, y_full);
    if (ret == GSL_SUCCESS)
    {
        *tmax = tmax_full;
        *ymax = creal(ymax_full);
    }
    return ret;
}


/*
 * General functions
 */


/* Return |z|^2 for a complex number z. */
static double cabs2(COMPLEX16 z)
{
    return gsl_pow_2(creal(z)) + gsl_pow_2(cimag(z));
}


/*
 * CubicSpline
 */


/* Provide declaration of opaque data structure to hold Lanzos interpolant state. */
struct tagCubicSplineTriggerInterpolant {
    gsl_poly_complex_workspace *workspace;
};


/* Data structure providing arguments for minimizer cost function. */
typedef struct {
    const COMPLEX16 *data;
    unsigned int window;
} CubicSplineTriggerInterpolantParams;


/* Strip leading zero coefficients of a polynomial.
 * Pass the number of coefficients as n. The new number of coefficients is returned. */
static size_t poly_strip(const double *a, size_t n)
{
    for (; n > 0 && a[n - 1] == 0; n--)
        ; /* loop body intentionally a no-op */
    return n;
}


/* Polynomial multiply-accumulate. */
static void poly_mac(double *y, const double *a, size_t na, const double *b, size_t nb)
{
    size_t ia, ib;

    /* Compute convolution. */
    for (ia = 0; ia < na; ia++)
        for (ib = 0; ib < nb; ib++)
            y[ia + ib] += a[ia] * b[ib];
}


/* Compute derivative of a polynomial. */
static void poly_der(double *a, size_t n)
{
    size_t i;
    for (i = 0; i < n; i ++)
        a[i] *= i;
}


/* Construct inteprolating polynomial. */
static void poly_interp(double *a, const double y_0, const double y_1, const double y_2, const double y_3)
{
    a[3] = 0.5 * (y_3 - y_0) + 1.5 * (y_1 - y_2);
    a[2] = y_0 - 2.5 * y_1 + 2 * y_2 - 0.5 * y_3;
    a[1] = 0.5 * (y_2 - y_0);
    a[0] = y_1;
}


/**
 * Treat \c are and \c aim as the real and imaginary parts of a polynomial with
 * \n complex coefficients. Find all local extrema of the absolute value of the
 * polynomial.
 */
static int interp_find_roots(size_t *nroots, COMPLEX16 *roots, const double *are, const double *aim, size_t n, gsl_poly_complex_workspace *workspace)
{
    int ret = GSL_EFAILED;
    double b[2 * n - 2];

    /* Compute the coefficients of the polynomial
     * b = D[|a|^2 + |b|^2, t]. */
    {
        double ad[n];
        memset(b, 0, sizeof(double) * (2 * n - 2));

        memcpy(ad, are, sizeof(double) * n);
        poly_der(ad, n);
        poly_mac(b, &ad[1], n - 1, are, n);

        memcpy(ad, aim, sizeof(double) * n);
        poly_der(ad, n);
        poly_mac(b, &ad[1], n - 1, aim, n);
    }

    n = poly_strip(b, 2 * n - 2);

    if (n == 0) {
        *nroots = 0;
        return GSL_SUCCESS;
    }

    ret = gsl_poly_complex_solve(b, n, workspace, (double *) roots);
    if (ret == GSL_SUCCESS)
        *nroots = n - 1;

    return ret;
}


/**
 * Perform cubic spline interpolation of a trigger. Since cubic splines
 * inherently take into account an even number of samples, we are going to
 * have to call this function twice: once for the first four of the five samples
 * surrounding the trigger and once for the last four of the five samples
 * surrounding the trigger.
 */
static int cubic_interp_1(gsl_poly_complex_workspace *workspace, double *t, COMPLEX16 *val, const COMPLEX16 *y)
{
    double argmax = NAN, new_argmax;
    COMPLEX16 maxval, new_maxval;
    double max_abs2, new_max_abs2;

    size_t n = 4;
    double are[n], aim[n];
    COMPLEX16 roots[2 * n - 3];

    size_t nroots, iroot;
    int result;

    /* Compute coefficients of interpolating polynomials for real and imaginary
     * parts of data. */
    poly_interp(are, creal(y[0]), creal(y[1]), creal(y[2]), creal(y[3]));
    poly_interp(aim, cimag(y[0]), cimag(y[1]), cimag(y[2]), cimag(y[3]));

    /* Find local maxima of (|a|^2 + |b|^2). */
    result = interp_find_roots(&nroots, roots, are, aim, n, workspace);
    if (result != GSL_SUCCESS)
        return result;

    /* Determine which of the endpoints is greater. */
    argmax = 0;
    maxval = y[1];
    max_abs2 = cabs2(maxval);

    new_argmax = 1;
    new_maxval = y[2];
    new_max_abs2 = cabs2(new_maxval);
    if (new_max_abs2 > max_abs2) {
        argmax = new_argmax;
        maxval = new_maxval;
        max_abs2 = new_max_abs2;
    }

    /* See if there is a local extremum that is greater than the endpoints. */
    for (iroot = 0; iroot < nroots; iroot++) {
        /* Skip this root if it is complex. */
        if (cimag(roots[iroot]) != 0)
            continue;

        new_argmax = creal(roots[iroot]);

        /* Skip this root if it is outside of (0, 1) */
        if (new_argmax <= 0 || new_argmax >= 1)
            continue;

        new_maxval = gsl_poly_eval(are, n, new_argmax) + gsl_poly_eval(aim, n, new_argmax) * I;
        new_max_abs2 = cabs2(new_maxval);

        if (new_max_abs2 > max_abs2) {
            argmax = new_argmax;
            maxval = new_maxval;
        }
    }

    *t = argmax;
    *val = maxval;
    return GSL_SUCCESS;
}


CubicSplineTriggerInterpolant *XLALCreateCubicSplineTriggerInterpolant(unsigned int window)
{
    CubicSplineTriggerInterpolant *interp = calloc(1, sizeof(CubicSplineTriggerInterpolant));

    if (window < 2)
        goto fail;

    if (!interp)
        goto fail;

    interp->workspace = gsl_poly_complex_workspace_alloc(6);
    if (!interp->workspace)
        goto fail;

    return interp;
fail:
    XLALDestroyCubicSplineTriggerInterpolant(interp);
    return NULL;
}


void XLALDestroyCubicSplineTriggerInterpolant(CubicSplineTriggerInterpolant *interp)
{
    if (interp)
    {
        gsl_poly_complex_workspace_free(interp->workspace);
        interp->workspace = NULL;
    }
    free(interp);
}


int XLALCOMPLEX16ApplyCubicSplineTriggerInterpolant(
    CubicSplineTriggerInterpolant *interp,
    double *t,
    COMPLEX16 *y,
    const COMPLEX16 *data)
{
    COMPLEX16 max1, max2;
    double max1_abs1, max2_abs2;
    double argmax1, argmax2;
    int result;

    result = cubic_interp_1(interp->workspace, &argmax1, &max1, &data[-2]);
    if (result != GSL_SUCCESS)
        return result;
    result = cubic_interp_1(interp->workspace, &argmax2, &max2, &data[-1]);
    if (result != GSL_SUCCESS)
        return result;
    max1_abs1 = cabs2(max1);
    max2_abs2 = cabs2(max2);

    if (max1_abs1 > max2_abs2) {
        *t = argmax1 - 1;
        *y = max1;
    } else {
        *t = argmax2;
        *y = max2;
    }

    return GSL_SUCCESS;
}


int XLALCOMPLEX8ApplyCubicSplineTriggerInterpolant(
    CubicSplineTriggerInterpolant *interp,
    double *tmax,
    COMPLEX8 *ymax,
    const COMPLEX8 *y)
{
    return XLALCOMPLEX8ApplyTriggerInterpolant(interp,
        (XLALCOMPLEX16ApplyFunc) XLALCOMPLEX16ApplyCubicSplineTriggerInterpolant,
        2, tmax, ymax, y);
}


int XLALREAL8ApplyCubicSplineTriggerInterpolant(
    CubicSplineTriggerInterpolant *interp,
    double *tmax,
    REAL8 *ymax,
    const REAL8 *y)
{
    return XLALREAL8ApplyTriggerInterpolant(interp,
        (XLALCOMPLEX16ApplyFunc) XLALCOMPLEX16ApplyCubicSplineTriggerInterpolant,
        2, tmax, ymax, y);
}


int XLALREAL4ApplyCubicSplineTriggerInterpolant(
    CubicSplineTriggerInterpolant *interp,
    double *tmax,
    REAL4 *ymax,
    const REAL4 *y)
{
    return XLALREAL4ApplyTriggerInterpolant(interp,
        (XLALCOMPLEX16ApplyFunc) XLALCOMPLEX16ApplyCubicSplineTriggerInterpolant,
        2, tmax, ymax, y);
}


/*
 * Lanczos
 */


/* Provide declaration of opaque data structure to hold Lanzos interpolant state. */
struct tagLanczosTriggerInterpolant {
    gsl_min_fminimizer *fminimizer;
    unsigned int window;
};


/* Data structure providing arguments for minimizer cost function. */
typedef struct {
    const COMPLEX16 *data;
    unsigned int window;
} LanczosTriggerInterpolantParams;


/* The Lanczos kernel. */
static double lanczos(double t, double a)
{
    return gsl_sf_sinc(t) * gsl_sf_sinc(t / a);
}


/* The Lanczos reconstruction filter interpolant. */
static COMPLEX16 lanczos_interpolant(double t, const LanczosTriggerInterpolantParams *params)
{
    COMPLEX16 ret;
    int i;

    for (ret = 0, i = -(int)params->window; i <= (int)params->window; i ++)
        ret += lanczos(t - i, params->window) * params->data[i];

    return ret;
}


/* The cost function to minimize. */
static double lanczos_cost(double t, void *params)
{
    return -cabs2(lanczos_interpolant(t, params));
}


LanczosTriggerInterpolant *XLALCreateLanczosTriggerInterpolant(unsigned int window)
{
    LanczosTriggerInterpolant *interp = calloc(1, sizeof(LanczosTriggerInterpolant));

    if (!interp)
        goto fail;

    interp->fminimizer = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
    if (!interp->fminimizer)
        goto fail;

    interp->window = window;

    return interp;
fail:
    XLALDestroyLanczosTriggerInterpolant(interp);
    return NULL;
}


void XLALDestroyLanczosTriggerInterpolant(LanczosTriggerInterpolant *interp)
{
    if (interp)
    {
        gsl_min_fminimizer_free(interp->fminimizer);
        interp->fminimizer = NULL;
    }
    free(interp);
}


int XLALCOMPLEX16ApplyLanczosTriggerInterpolant(
    LanczosTriggerInterpolant *interp,
    double *t,
    COMPLEX16 *y,
    const COMPLEX16 *data)
{
    static const double epsabs = 1e-5;

    LanczosTriggerInterpolantParams params = {data, interp->window};
    gsl_function func = {lanczos_cost, &params};
    double t1, t2;
    int result;

    result = gsl_min_fminimizer_set_with_values(interp->fminimizer, &func,
        0, -cabs2(data[0]), -1, -cabs2(data[-1]), 1, -cabs2(data[1]));
    if (result != GSL_SUCCESS)
        GSL_ERROR("failed to initialize minimizer", result);

    do {
        result = gsl_min_fminimizer_iterate(interp->fminimizer);
        if (result != GSL_SUCCESS)
            GSL_ERROR("failed to perform minimizer iteration", result);

        t1 = gsl_min_fminimizer_x_lower(interp->fminimizer);
        t2 = gsl_min_fminimizer_x_upper(interp->fminimizer);
    } while (t2 - t1 > epsabs);

    *t = gsl_min_fminimizer_x_minimum(interp->fminimizer);
    *y = lanczos_interpolant(*t, &params);

    return GSL_SUCCESS;
}


int XLALCOMPLEX8ApplyLanczosTriggerInterpolant(
    LanczosTriggerInterpolant *interp,
    double *tmax,
    COMPLEX8 *ymax,
    const COMPLEX8 *y)
{
    return XLALCOMPLEX8ApplyTriggerInterpolant(interp,
        (XLALCOMPLEX16ApplyFunc) XLALCOMPLEX16ApplyLanczosTriggerInterpolant,
        interp->window, tmax, ymax, y);
}


int XLALREAL8ApplyLanczosTriggerInterpolant(
    LanczosTriggerInterpolant *interp,
    double *tmax,
    REAL8 *ymax,
    const REAL8 *y)
{
    return XLALREAL8ApplyTriggerInterpolant(interp,
        (XLALCOMPLEX16ApplyFunc) XLALCOMPLEX16ApplyLanczosTriggerInterpolant,
        interp->window, tmax, ymax, y);
}


int XLALREAL4ApplyLanczosTriggerInterpolant(
    LanczosTriggerInterpolant *interp,
    double *tmax,
    REAL4 *ymax,
    const REAL4 *y)
{
    return XLALREAL4ApplyTriggerInterpolant(interp,
        (XLALCOMPLEX16ApplyFunc) XLALCOMPLEX16ApplyLanczosTriggerInterpolant,
        interp->window, tmax, ymax, y);
}


/*
 * Nearest neighbor
 */


struct tagNearestNeighborTriggerInterpolant {
    unsigned int window;
};


NearestNeighborTriggerInterpolant *XLALCreateNearestNeighborTriggerInterpolant(unsigned int window)
{
    NearestNeighborTriggerInterpolant *interp = calloc(1, sizeof(NearestNeighborTriggerInterpolant));

    if (!interp)
        goto fail;

    if (window != 0)
        goto fail;

    interp->window = window;

    return interp;
fail:
    XLALDestroyNearestNeighborTriggerInterpolant(interp);
    return NULL;
}


void XLALDestroyNearestNeighborTriggerInterpolant(NearestNeighborTriggerInterpolant *interp)
{
    free(interp);
}


int XLALCOMPLEX16ApplyNearestNeighborTriggerInterpolant(
    __attribute__((unused)) NearestNeighborTriggerInterpolant *interp,
    double *t,
    COMPLEX16 *y,
    const COMPLEX16 *data)
{
    *t = 0;
    *y = data[0];
    return GSL_SUCCESS;
}


int XLALCOMPLEX8ApplyNearestNeighborTriggerInterpolant(
    NearestNeighborTriggerInterpolant *interp,
    double *tmax,
    COMPLEX8 *ymax,
    const COMPLEX8 *y)
{
    return XLALCOMPLEX8ApplyTriggerInterpolant(interp,
        (XLALCOMPLEX16ApplyFunc) XLALCOMPLEX16ApplyNearestNeighborTriggerInterpolant,
        interp->window, tmax, ymax, y);
}


int XLALREAL8ApplyNearestNeighborTriggerInterpolant(
    NearestNeighborTriggerInterpolant *interp,
    double *tmax,
    REAL8 *ymax,
    const REAL8 *y)
{
    return XLALREAL8ApplyTriggerInterpolant(interp,
        (XLALCOMPLEX16ApplyFunc) XLALCOMPLEX16ApplyNearestNeighborTriggerInterpolant,
        interp->window, tmax, ymax, y);
}


int XLALREAL4ApplyNearestNeighborTriggerInterpolant(
    NearestNeighborTriggerInterpolant *interp,
    double *tmax,
    REAL4 *ymax,
    const REAL4 *y)
{
    return XLALREAL4ApplyTriggerInterpolant(interp,
        (XLALCOMPLEX16ApplyFunc) XLALCOMPLEX16ApplyNearestNeighborTriggerInterpolant,
        interp->window, tmax, ymax, y);
}


/*
 * Quadratic fit
 */


struct tagQuadraticFitTriggerInterpolant {
    gsl_multifit_linear_workspace *workspace;
    gsl_matrix *X;
    gsl_matrix *cov;
    gsl_vector *y;
    gsl_vector *c;
    unsigned int window;
};


QuadraticFitTriggerInterpolant *XLALCreateQuadraticFitTriggerInterpolant(unsigned int window)
{
    QuadraticFitTriggerInterpolant *interp = calloc(1, sizeof(QuadraticFitTriggerInterpolant));
    int i;

    if (!interp)
        goto fail;

    if (window < 1)
        goto fail;

    interp->window = window;

    interp->workspace = gsl_multifit_linear_alloc(2 * window + 1, 3);
    if (!interp->workspace)
        goto fail;

    interp->X = gsl_matrix_alloc(2 * window + 1, 3);
    if (!interp->X)
        goto fail;

    for (i = 0; i < 2 * (int)window + 1; i ++)
    {
        double x = i - window;
        gsl_matrix_set(interp->X, i, 0, 1);
        gsl_matrix_set(interp->X, i, 1, x);
        gsl_matrix_set(interp->X, i, 2, gsl_pow_2(x));
    }

    interp->cov = gsl_matrix_alloc(3, 3);
    if (!interp->cov)
        goto fail;

    interp->y = gsl_vector_alloc(2 * window + 1);
    if (!interp->y)
        goto fail;

    interp->c = gsl_vector_alloc(3);
    if (!interp->c)
        goto fail;

    interp->window = window;

    return interp;
fail:
    XLALDestroyQuadraticFitTriggerInterpolant(interp);
    return NULL;
}


void XLALDestroyQuadraticFitTriggerInterpolant(QuadraticFitTriggerInterpolant *interp)
{
    if (interp)
    {
        gsl_multifit_linear_free(interp->workspace);
        interp->workspace = NULL;
        gsl_matrix_free(interp->X);
        interp->X = NULL;
        gsl_matrix_free(interp->cov);
        interp->cov = NULL;
        gsl_vector_free(interp->y);
        interp->y = NULL;
        gsl_vector_free(interp->c);
        interp->c = NULL;
    }
    free(interp);
}


int XLALCOMPLEX16ApplyQuadraticFitTriggerInterpolant(
    QuadraticFitTriggerInterpolant *interp,
    double *t,
    COMPLEX16 *y,
    const COMPLEX16 *data)
{
    int i;
    int result;
    double chisq;
    double a, b, tmax;

    for (i = -(int)interp->window; i <= (int)interp->window; i ++)
        gsl_vector_set(interp->y, i + interp->window, cabs(y[i]));

    result = gsl_multifit_linear(interp->X, interp->y, interp->c, interp->cov, &chisq, interp->workspace);
    if (result != GSL_SUCCESS)
        GSL_ERROR("regression failed", result);

    a = gsl_vector_get(interp->c, 2);
    b = gsl_vector_get(interp->c, 1);
    tmax = -0.5 * b / a;

    if ((a > 0 || (a == 0 && b > 0)) && tmax > -1 && tmax < 1)
        *t = tmax;
    else
        *t = 0;

    *y = data[0];

    return GSL_SUCCESS;
}


int XLALCOMPLEX8ApplyQuadraticFitTriggerInterpolant(
    QuadraticFitTriggerInterpolant *interp,
    double *tmax,
    COMPLEX8 *ymax,
    const COMPLEX8 *y)
{
    return XLALCOMPLEX8ApplyTriggerInterpolant(interp,
        (XLALCOMPLEX16ApplyFunc) XLALCOMPLEX16ApplyQuadraticFitTriggerInterpolant,
        interp->window, tmax, ymax, y);
}


int XLALREAL8ApplyQuadraticFitTriggerInterpolant(
    QuadraticFitTriggerInterpolant *interp,
    double *tmax,
    REAL8 *ymax,
    const REAL8 *y)
{
    return XLALREAL8ApplyTriggerInterpolant(interp,
        (XLALCOMPLEX16ApplyFunc) XLALCOMPLEX16ApplyQuadraticFitTriggerInterpolant,
        interp->window, tmax, ymax, y);
}


int XLALREAL4ApplyQuadraticFitTriggerInterpolant(
    QuadraticFitTriggerInterpolant *interp,
    double *tmax,
    REAL4 *ymax,
    const REAL4 *y)
{
    return XLALREAL4ApplyTriggerInterpolant(interp,
        (XLALCOMPLEX16ApplyFunc) XLALCOMPLEX16ApplyQuadraticFitTriggerInterpolant,
        interp->window, tmax, ymax, y);
}
