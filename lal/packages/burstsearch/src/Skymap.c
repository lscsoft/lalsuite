/*
 * $Id$
 *
 * Copyright (C) 2008  Antony Searle
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <lal/LALConstants.h>
#include <lal/LALMalloc.h>
#include <lal/XLALError.h>
#include <lal/DetResponse.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Skymap.h>

// Convenience functions for tiny stack vectors and matrices

// Dot product of 3-vectors

static double dot3(double a[3], double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Inverse of a 2x2 matrix

static void inv22(double a[2][2], double b[2][2])
{
    double c = b[0][0] * b[1][1] - b[0][1] * b[1][0];
    a[0][0] =  b[1][1] / c;
    a[0][1] = -b[0][1] / c;
    a[1][0] = -b[1][0] / c;
    a[1][1] =  b[0][0] / c;
}

// Determinant of a 2x2 matrix

static double det22(double a[2][2])
{
    return a[0][0] * a[1][1] - a[0][1] * a[1][0];
}

// Coordinate transformation theta,phi -> x,y,z

void XLALSkymapCartesianFromSpherical(double a[3], double b[2])
{
    a[0] = sin(b[0]) * cos(b[1]);
    a[1] = sin(b[0]) * sin(b[1]);
    a[2] = cos(b[0]);
}

// Coordinate transformation x,y,z -> theta,phi

void XLALSkymapSphericalFromCartesian(double a[2], double b[3])
{
    a[0] = acos(b[2]);
    a[1] = atan2(b[1], b[0]);
}

// Time of arrival at a detector, relative to center of Earth, for signals
// from a given (unit) direction

static double site_time(LALDetector* site, double direction[3])
{
    return -dot3(site->location, direction) / LAL_C_SI;
}

// Plus and cross polarization response of a detector for signals from a
// given direction

static void site_response(double f[2], LALDetector* site, double direction[3])
{
    double thetaphi[2];
    XLALSkymapSphericalFromCartesian(thetaphi, direction);
    XLALComputeDetAMResponse(&f[0], &f[1], site->response, thetaphi[1], LAL_PI_2 - thetaphi[0], 0, 0);
}

// XLALSkymapLogSumExp(a, b) computes log(exp(a) + exp(b)) but will not
// overflow for a or b > ~300
//
// For a > b, we use the identity
//
// log(exp(a) + exp(b)
//     = log(exp(a) * (1 + exp(b) / exp(a)))
//     = a + log(1 + exp(b - a))
//     = a + log1p(exp(b - a))
//
// where b - a < 0 and exp(b - a) cannot overflow (though it may
// underflow).
//
// And for a < b
//
// log(exp(a) + exp(b))
//     = log(exp(b) * (exp(a) / exp(b) + 1))
//     = b + log(exp(a - b) + 1)
//     = b + log1p(exp(a - b))
//
// If neither a < b nor a > b, we either have equality or both values
// are (the same) plus or minus infinity.  Forming (a - b) in the case of
// infinities results in a NaN, so we must use a third expression
//
// log(exp(a) + exp(b))
//     = log(exp(a) + exp(a))
//     = log(2 * exp(a))
//     = log(2) + a

double XLALSkymapLogSumExp(double a, double b)
{
    return (a < b) ?
        (b + log1p(exp(a - b))) :
        ((b < a) ? (a + log1p(exp(b - a))) : (a + log(2)));
}

// Find the maximum of a sequence  of doubles and return a pointer to it

static double* findmax(double* begin, double* end)
{
    double* p;
    double* m;
    m = begin;
    for (p = begin; p != end; ++p)
    {
        if (*m < *p)
        {
            m = p;
        }
    }
    return m;
}

// Find log sum_i exp(a[i])) given max_i a[i], using the same technique as
// XLALSkymapLogSumExp to accumulate the result without overflowing for
// a[i] > ~300

static double logtotalexpwithmax(double* begin, double* end, double m)
{
    double t;
    double* p;
    t = 0;
    for (p = begin; p != end; ++p)
    {
        t += exp(*p - m);
    }
    return m + log(t);
}

// Find log sum_i exp(a[i]) using findmax and logtotalexpwithmax

double XLALSkymapLogTotalExp(double* begin, double* end)
{
    return logtotalexpwithmax(begin, end, *findmax(begin, end));
}

// To cubic interpolate
//     x(0 <= t <= 1)
// from
//     x[-1], x[0], x[1], x[2]
// we compute
//     w_t[-1], w_t[0], w_t[1], w_t[2]
// such that
//     x(t) = sum_i x[i] w_t[i]

double XLALSkymapInterpolate(double t, double* x)
{
    int whole = floor(t);
    t -= whole;

    double h[4], w[4];

    // Hermite basis functions at t

    h[0] = (1. + 2. * t) * (1. - t) * (1. - t);
    h[1] = t * (1. - t) * (1. - t);
    h[2] = t * t * (3. - 2. * t);
    h[3] = t * t * (t - 1.);

    // Weights

    w[0] = -0.5 * h[1];
    w[1] = h[0] - 0.5 * h[3];
    w[2] = h[2] + 0.5 * h[1];
    w[3] = 0.5 * h[3];

    double y = 0;
    int i;

    for (i = 0; i != 4; ++i)
    {
        y += x[whole + i - 1] * w[i];
    }

    return y;


}

// Construct an XLALSkymap2PlanType in the given memory from
//     sample frequency of matched filter timeseries
//     number of detectors
//     list of detector LAL ID numbers

void XLALSkymapPlanConstruct(int sampleFrequency, int n, int* detectors, XLALSkymapPlanType* plan)
{
    int i;

    plan->sampleFrequency = sampleFrequency;
    plan->n = n;

    for (i = 0; i != plan->n; ++i)
    {
        plan->site[i] = lalCachedDetectors[detectors[i]];
    }

}

// Construct an XLALSkymap2DirectionProperties object in the given memory
// from
//     a plan
//     a theta, phi tuple

void XLALSkymapDirectionPropertiesConstruct(
        XLALSkymapPlanType* plan,
        double* directions,
        XLALSkymapDirectionPropertiesType* properties
        )
{
    double x[3];
    int j;

    // Convert theta, phi direction to cartesian
    XLALSkymapCartesianFromSpherical(x, directions);

    for (j = 0; j != plan->n; ++j)
    {
        // Delay (in seconds)
        properties->delay[j] = site_time(plan->site + j, x);
        // Plus and cross polarization responses
        site_response(properties->f[j], plan->site + j, x);
    }
}

// Construct a XLALSkymap2Kernel object in the given memory from
//     a plan
//     direction properties (time delays, interpolation weights and
//         antenna patterns)
//     the noise-weighted inner product of the template with itself

void XLALSkymapKernelConstruct(
        XLALSkymapPlanType* plan,
        XLALSkymapDirectionPropertiesType* properties,
        double* wSw,
        XLALSkymapKernelType* kernel
        )
{

    int i, j, k, l;

    //
    // W = diag(w.S_j^{-1}.w)
    //
    // F (F^T W F + I) F^T
    //

    double a[2][2]; // F^T W F
    double b[2][2]; // inv(A)

    // Compute the kernel

    // F^T W F + I

    for (i = 0; i != 2; ++i)
    {
        for (j = 0; j != 2; ++j)
        {
            a[i][j] = ((i == j) ? 1.0 : 0.0);
            for (k = 0; k != plan->n; ++k)
            {
                a[i][j] += properties->f[k][i] * wSw[k] * properties->f[k][j];
            }
        }
    }

    // (F^T W F + I)^{-1}

    inv22(b, a);

    // F (F^T W F + I)^{-1} F^T

    for (i = 0; i != plan->n; ++i)
    {
        for (j = 0; j != plan->n; ++j)
        {
            kernel->k[i][j] = 0.0;
            for (k = 0; k != 2; ++k)
            {
                for (l = 0; l != 2; ++l)
                {
                    kernel->k[i][j] += properties->f[i][k] * b[k][l] * properties->f[j][l];
                }
            }
        }
    }

    kernel->logNormalization = 0.5 * log(det22(b));

}

// Compute the marginalization integral over a_plus and a_cross for the
// system described by
//     a plan
//     a direction's properties
//     a kernel
//     a matched filter time series for each detector
//     a signal arrival time

void XLALSkymapApply(
        XLALSkymapPlanType* plan,
        XLALSkymapDirectionPropertiesType* properties,
        XLALSkymapKernelType* kernel,
        double** xSw,
        double tau,
        double* posterior
        )
{
    double a;
    int i, j;

    double x[XLALSKYMAP_N];

    // Interpolate the matched filter values

    for (i = 0; i != plan->n; ++i)
    {
        x[i] = XLALSkymapInterpolate((tau + properties->delay[i]) * plan->sampleFrequency, xSw[i]);
    }

    // This implementation does not exploit the symmetry of the expression

    // Compute x^T.K.x

    a = 0;

    for (i = 0; i != plan->n; ++i)
    {
        for (j = 0; j != plan->n; ++j)
        {
            a += x[i] * kernel->k[i][j] * x[j];
        }
    }

    // Scale and apply the normalization

    *posterior = 0.5 * a + kernel->logNormalization;

}


