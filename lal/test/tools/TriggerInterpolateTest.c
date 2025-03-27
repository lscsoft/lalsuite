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
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_test.h>

#include <lal/TriggerInterpolate.h>


static double cpart(double complex z, int i)
{
    return ((double *)&z)[i];
}


static double to_phase_angle(double arg)
{
    arg = fmod(arg, 2 * M_PI);
    if (arg < -M_PI)
        arg += 2 * M_PI;
    else if (arg > M_PI)
        arg -= 2 * M_PI;
    return arg;
}


// Derivative of normalized sinc function.
static double dsinc(double t)
{
    return (cos(M_PI * t) - gsl_sf_sinc(t)) / t;
}


// Second derivative of normalized sinc function.
static double d2sinc(double t)
{
    return -(sin(M_PI * t) * M_PI + 2 * dsinc(t)) / t;
}


// Value and derivative of Lanczos kernel.
static void fdf_lanczos(double fdf[3], double t, double a)
{
    if (t < -a || t > a)
    {
        fdf[0] = fdf[1] = fdf[2] = 0;
    } else {
        double ta = t / a;
        double sinct = gsl_sf_sinc(t), sincta = gsl_sf_sinc(ta);
        double dsinct = dsinc(t), dsincta = dsinc(ta) / a;
        double d2sinct = d2sinc(t), d2sincta = d2sinc(ta) / gsl_pow_2(a);
        fdf[0] = sinct * sincta;
        fdf[1] = sinct * dsincta + sincta * dsinct;
        fdf[2] = 2 * dsinct * dsincta + sinct * d2sincta + sincta * d2sinct;
    }
}


static void print_gsl_error(
    const char *reason,
    const char *file,
    int line,
    __attribute__ ((unused)) int gsl_errno)
{
    gsl_stream_printf("ERROR", file, line, reason);
}


#define TEST_WINDOW_SIZE_TOO_SMALL(FUNC, MIN_WINDOW) \
{ \
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off(); \
    for (unsigned int window = 0; window < MIN_WINDOW; window ++) \
    { \
        double complex y[2 * window + 1]; \
        result = FUNC(&tmax, &ymax, &y[window], window); \
        gsl_test(result != GSL_EINVAL, "%s window=%u returns GSL_EINVAL", #FUNC, window); \
    } \
    gsl_set_error_handler(old_handler); \
}


// Draw a random quadratic polynomial for the amplitude:
//   - It has an absolute maximum at t_max, which is in (-1, +1).
//   - The value at the maximum is amp_max.
//   - It is greater than zero over the range (-window, window).
static void ran_amp_poly_quadratic(
    gsl_rng *rng,
    double poly[3],
    double *t_max,
    double *amp_max,
    unsigned int window
) {
    double t = gsl_ran_flat(rng, -1, 1);
    double a = -gsl_ran_gamma(rng, 2, 1);
    double b = -2 * a * t;
    double c = gsl_ran_gamma(rng, 2, 1) - a * window * (window + 2 * fabs(t));
    poly[0] = c;
    poly[1] = b;
    poly[2] = a;
    *t_max = t;
    *amp_max = c - 0.25 * gsl_pow_2(b) / a;
}


// Draw a random linear polynomial for the phase:
//   - The absolute value of the slope is no greater than pi.
static void ran_arg_poly_linear(gsl_rng *rng, double poly[2])
{
    poly[0] = gsl_ran_flat(rng, -M_PI, M_PI);
    poly[1] = gsl_ran_flat(rng, -M_PI, M_PI);
}


int main(__attribute__ ((unused)) int argc, __attribute__ ((unused)) char **argv)
{
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    int result;
    double complex ymax;
    double tmax;

    gsl_set_error_handler(print_gsl_error);

    TEST_WINDOW_SIZE_TOO_SMALL(XLALTriggerInterpolateCubicSpline, 2)
    TEST_WINDOW_SIZE_TOO_SMALL(XLALTriggerInterpolateCubicSplineAmpPhase, 2)
    TEST_WINDOW_SIZE_TOO_SMALL(XLALTriggerInterpolateQuadraticFit, 1)
    TEST_WINDOW_SIZE_TOO_SMALL(XLALTriggerInterpolateLanczos, 1)

    /* Test Lanczos interpolant. */
    for (unsigned int window = 1; window < 50; window ++)
    {
        unsigned int len = 2 * window + 1;
        double complex y[len], (*data) = &y[window], central_datum;

        // Random values for all samples
        for (unsigned int i = 0; i < len; i ++)
            y[i] = crect(gsl_ran_gaussian(rng, 1), gsl_ran_gaussian(rng, 1));

        // Normalized value of central datum
        central_datum = data[0];
        central_datum /= cabs(central_datum);

        // Pick a random maximum time and phase.
        double tmax_expected = gsl_ran_flat(rng, -1, 1);

        // Solve for the absolute value of the central point that produces
        // a stationary point at tmax_expected.
        double complex pdp[3] = {0}, qdq[3];
        double a = 0, b = 0, c = 0, a1 = 0, b1 = 0, c1 = 0, roots[2], fdf[3];
        for (unsigned int i = 1; i <= window; i ++)
        {
            for (int sign = -1; sign <= 1; sign += 2)
            {
                double complex datum = data[sign * (int)i];
                fdf_lanczos(fdf, tmax_expected - sign * (int)i, window);
                for (int deriv = 0; deriv < 3; deriv++)
                    pdp[deriv] += datum * fdf[deriv];
            }
        }
        fdf_lanczos(fdf, tmax_expected, window);
        for (int deriv = 0; deriv < 3; deriv ++)
            qdq[deriv] = central_datum * fdf[deriv];
        for (int part = 0; part < 2; part ++)
        {
            // Quadratic equation for stationary point
            a += cpart(qdq[0], part) * cpart(qdq[1], part);
            b += cpart(pdp[0], part) * cpart(qdq[1], part) + cpart(pdp[1], part) * cpart(qdq[0], part);
            c += cpart(pdp[0], part) * cpart(pdp[1], part);

            // Quadratic equation for second derivative at stationary point
            a1 += gsl_pow_2(cpart(qdq[1], part)) + cpart(qdq[0], part) * cpart(qdq[2], part);
            b1 += 2 * cpart(pdp[1], part) * cpart(qdq[1], part) + cpart(pdp[0], part) * cpart(qdq[2], part) + cpart(pdp[2], part) * cpart(qdq[0], part);
            c1 += gsl_pow_2(cpart(pdp[1], part)) + cpart(pdp[0], part) * cpart(pdp[2], part);
        }
        int nroots = gsl_poly_solve_quadratic(a, b, c, &roots[0], &roots[1]);
        for (int iroot = 0; iroot < nroots; iroot ++)
        {
            double root = roots[iroot];
            double complex ymax_expected = pdp[0] + root * qdq[0];
            double abs_expected = cabs(ymax_expected);
            data[0] = root * central_datum;

            double second_deriv = a1 * gsl_pow_2(root) + b1 * root + c1;

            // Skip if this is a saddle point or local minimum.
            if (second_deriv >= 0)
                continue;

            // Skip if the local maximum is not greater than the sample points.
            if (abs_expected <= GSL_MAX_DBL(cabs(data[0]), GSL_MAX_DBL(cabs(data[-1]), cabs(data[1]))))
                continue;

            // Skip if the new sample point is not greater than the
            // two enclosing points.
            if (fabs(root) <= GSL_MAX_DBL(cabs(data[-1]), cabs(data[1])))
                continue;

            result = XLALTriggerInterpolateLanczos(&tmax, &ymax, data, window);
            gsl_test(result, "XLALTriggerInterpolateLanczos random signal window=%u return value", window);
            gsl_test_abs(tmax, tmax_expected, 1e-6, "XLALTriggerInterpolateLanczos random signal window=%u tmax", window);
            gsl_test_abs(cabs(ymax), cabs(ymax_expected), 1e-6, "XLALTriggerInterpolateLanczos random signal window=%u abs", window);
            gsl_test_abs(carg(ymax), carg(ymax_expected), 1e-6, "XLALTriggerInterpolateLanczos random signal window=%u arg", window);
        }
    }

    /* Test on a signal whose absolute value is itself quadratic, and that has a maximum in [-1, 1]. */
    for (unsigned int window = 1; window < 100; window ++)
    {
        unsigned int len = 2 * window + 1;
        double amp_poly[3], arg_poly[2], tmax_expected, amp_expected;
        ran_amp_poly_quadratic(rng, amp_poly, &tmax_expected, &amp_expected, window);
        ran_arg_poly_linear(rng, arg_poly);
        double arg_expected = to_phase_angle(gsl_poly_eval(arg_poly, 2, tmax_expected));
        double complex y[len];
        for (unsigned int i = 0; i < len; i ++)
        {
            double x = (int)i - (int)window;
            double amp = gsl_poly_eval(amp_poly, 3, x);
            double arg = gsl_poly_eval(arg_poly, 2, x);
            y[i] = cpolar(amp, arg);
        }

        result = XLALTriggerInterpolateQuadraticFit(&tmax, &ymax, &y[window], window);
        gsl_test(result, "XLALTriggerInterpolateQuadraticFit window=%u return value", window);
        gsl_test_abs(tmax, tmax_expected, 1e3 * GSL_DBL_EPSILON, "XLALTriggerInterpolateQuadraticFit quadratic signal window=%u tmax", window);
        gsl_test_abs(cabs(ymax), amp_expected, 1e5 * GSL_DBL_EPSILON, "XLALTriggerInterpolateQuadraticFit quadratic signal window=%u abs", window);
        gsl_test_abs(carg(ymax), arg_expected, 1e5 * GSL_DBL_EPSILON, "XLALTriggerInterpolateQuadraticFit quadratic signal window=%u arg", window);

        if (window >= 2)
        {
            result = XLALTriggerInterpolateCubicSplineAmpPhase(&tmax, &ymax, &y[window], window);
            gsl_test(result, "XLALTriggerInterpolateCubicSplineAmpPhase quadratic signal window=%u return value", window);
            gsl_test_abs(tmax, tmax_expected, 1e5 * GSL_DBL_EPSILON, "XLALTriggerInterpolateCubicSplineAmpPhase quadratic signal window=%u tmax", window);
            gsl_test_abs(cabs(ymax), amp_expected, 1e5 * GSL_DBL_EPSILON, "XLALTriggerInterpolateCubicSplineAmpPhase quadratic signal window=%u abs", window);
            gsl_test_abs(carg(ymax), arg_expected, 1e5 * GSL_DBL_EPSILON, "XLALTriggerInterpolateCubicSplineAmpPhase quadratic signal window=%u arg", window);
        }
    }

    /* Test on a signal whose real part is quadratic. */
    for (int run = 0; run < 100; run ++)
    {
        unsigned int window = 2, len = 2 * window + 1;
        double amp_poly[3], tmax_expected, re_expected, im_expected = 0;
        ran_amp_poly_quadratic(rng, amp_poly, &tmax_expected, &re_expected, window);
        double complex y[len];
        for (unsigned int i = 0; i < len; i ++)
        {
            double x = (int)i - (int)window;
            y[i] = gsl_poly_eval(amp_poly, 3, x);
        }

        result = XLALTriggerInterpolateCubicSpline(&tmax, &ymax, &y[window], window);
        gsl_test(result, "XLALTriggerInterpolateCubicSpline quadratic real part window=%u return value", window);
        gsl_test_abs(tmax, tmax_expected, 1e6 * GSL_DBL_EPSILON, "XLALTriggerInterpolateCubicSpline quadratic real part window=%u tmax", window);
        gsl_test_abs(creal(ymax), re_expected, 1e6 * GSL_DBL_EPSILON, "XLALTriggerInterpolateCubicSpline quadratic real part window=%u real", window);
        gsl_test_abs(cimag(ymax), im_expected, 0, "XLALTriggerInterpolateCubicSpline quadratic real part window=%u imag", window);
    }

    /* Test on a signal whose absolute value is linear. */
    for (unsigned int window = 1; window < 100; window ++)
    {
        unsigned int len = 2 * window + 1;
        int sign = 2 * gsl_ran_bernoulli(rng, 0.5) - 1;
        double complex y[len];
        double tmax_expected = sign;
        double real_expected = window + 1;
        double imag_expected = 0;
        for (unsigned int i = 0; i < len; i ++)
        {
            double x = (int)i - (int)window;
            y[i] = crect(sign * x + window, 0);
        }

        result = XLALTriggerInterpolateQuadraticFit(&tmax, &ymax, &y[window], window);
        gsl_test(result, "XLALTriggerInterpolateQuadraticFit linear signal window=%u return value", window);
        gsl_test_abs(tmax, tmax_expected, 0, "XLALTriggerInterpolateQuadraticFit linear signalwindow=%u tmax", window);
        gsl_test_abs(creal(ymax), real_expected, 0, "XLALTriggerInterpolateQuadraticFit linear signal window=%u real", window);
        gsl_test_abs(cimag(ymax), imag_expected, 0, "XLALTriggerInterpolateQuadraticFit linear signal window=%u imag", window);

        if (window >= 2)
        {
            result = XLALTriggerInterpolateCubicSpline(&tmax, &ymax, &y[window], window);
            gsl_test(result, "XLALTriggerInterpolateCubicSpline linear signal window=%u return value", window);
            gsl_test_abs(tmax, tmax_expected, 0, "XLALTriggerInterpolateCubicSpline linear signal window=%u tmax", window);
            gsl_test_abs(creal(ymax), real_expected, 0, "XLALTriggerInterpolateCubicSpline linear signal window=%u imag", window);
            gsl_test_abs(cimag(ymax), imag_expected, 0, "XLALTriggerInterpolateCubicSpline linear signal window=%u imag", window);

            result = XLALTriggerInterpolateCubicSplineAmpPhase(&tmax, &ymax, &y[window], window);
            gsl_test(result, "XLALTriggerInterpolateCubicSplineAmpPhase linear signal window=%u return value", window);
            gsl_test_abs(tmax, tmax_expected, 0, "XLALTriggerInterpolateCubicSplineAmpPhase linear signal window=%u tmax", window);
            gsl_test_abs(creal(ymax), real_expected, 0, "XLALTriggerInterpolateCubicSplineAmpPhase linear signal window=%u imag", window);
            gsl_test_abs(cimag(ymax), imag_expected, 0, "XLALTriggerInterpolateCubicSplineAmpPhase linear signal window=%u imag", window);
        }
    }

    /* NearestNeighbor */
    for (int i = 0; i < 100; i ++)
    {
        double complex y = crect(gsl_ran_gaussian(rng, 1), gsl_ran_gaussian(rng, 1));
        result = XLALTriggerInterpolateNearestNeighbor(&tmax, &ymax, &y, 0);
        gsl_test(result, "XLALTriggerInterpolateNearestNeighbor return value");
        gsl_test_abs(tmax, 0, 0, "XLALTriggerInterpolateNearestNeighbor tmax");
        gsl_test_abs(creal(ymax), creal(y), 0, "XLALTriggerInterpolateNearestNeighbor real part");
        gsl_test_abs(cimag(ymax), cimag(y), 0, "XLALTriggerInterpolateNearestNeighbor imag part");
    }

    return gsl_test_summary();
}
