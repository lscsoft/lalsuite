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

#define pi LAL_PI

/* functions to handle vectors and matrices with dimensions that are
   (a) small and
   (b) known at compile time */

#define max(A,B) (((A) > (B)) ? (A) : (B))

static double sq(double a)
{
    return a * a;
}

static void set2(double a[2], double v0, double v1)
{
    a[0] = v0;
    a[1] = v1;
}

static void set3(double a[3], double v0, double v1, double v2)
{
    a[0] = v0;
    a[1] = v1;
    a[2] = v2;
}

static void add3(double a[3], double b[3], double c[3])
{
    a[0] = b[0] + c[0];
    a[1] = b[1] + c[1];
    a[2] = b[2] + c[2];
}

static void sub3(double a[3], double b[3], double c[3])
{
    a[0] = b[0] - c[0];
    a[1] = b[1] - c[1];
    a[2] = b[2] - c[2];
}

/* unused
static void mul3(double a[3], double b[3], double c)
{
    a[0] = b[0] * c;
    a[1] = b[1] * c;
    a[2] = b[2] * c;
}
*/

static void div2(double a[2], double b[2], double c)
{
    a[0] = b[0] / c;
    a[1] = b[1] / c;
}

static void div3(double a[3], double b[3], double c)
{
    a[0] = b[0] / c;
    a[1] = b[1] / c;
    a[2] = b[2] / c;
}

static double dot3(double a[3], double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static double sqn3(double a[3])
{
    return dot3(a, a);
}

static double nrm3(double a[3])
{
    return sqrt(sqn3(a));
}

static void normalize3(double a[3])
{
    div3(a, a, nrm3(a));
}

static void cross(double a[3], double b[3], double c[3])
{
    a[0] = b[1] * c[2] - b[2] * c[1];
    a[1] = b[2] * c[0] - b[0] * c[2];
    a[2] = b[0] * c[1] - b[1] * c[0];
}

static void set32(double a[3][2],
           double v00, double v01,
           double v10, double v11,
           double v20, double v21)
{
    set2(a[0], v00, v01);
    set2(a[1], v10, v11);
    set2(a[2], v20, v21);
}

static void div22(double a[2][2], double b[2][2], double c)
{
    div2(a[0], b[0], c);
    div2(a[1], b[1], c);
}

static void transpose32(double a[2][3], double b[3][2])
{
    a[0][0] = b[0][0];
    a[0][1] = b[1][0];
    a[0][2] = b[2][0];
    a[1][0] = b[0][1];
    a[1][1] = b[1][1];
    a[1][2] = b[2][1];
}

static void set22(double a[2][2], double a00, double a01, double a10, double a11)
{
    a[0][0] = a00;
    a[0][1] = a01;
    a[1][0] = a10;
    a[1][1] = a11;
}

static void inv22(double a[2][2], double b[2][2])
{
    a[0][0] = b[1][1];
    a[0][1] = -b[0][1];
    a[1][0] = -b[1][0];
    a[1][1] = b[0][0];
    div22(a, a, b[0][0] * b[1][1] - b[0][1] * b[1][0]);
}

static void mul322(double a[3][2], double b[3][2], double c[2][2])
{
    int i, j, k;
    for (i = 0; i != 3; ++i)
        for (k = 0; k != 2; ++k)
        {
            a[i][k] = 0;
            for (j = 0; j != 2; ++j)
                a[i][k] += b[i][j] * c[j][k];
        }
}

static void mul323(double a[3][3], double b[3][2], double c[2][3])
{
    int i, j, k;
    for (i = 0; i != 3; ++i)
        for (k = 0; k != 3; ++k)
        {
            a[i][k] = 0;
            for (j = 0; j != 2; ++j)
                a[i][k] += b[i][j] * c[j][k];
        }
}

static void mul222(double a[2][2], double b[2][2], double c[2][2])
{
    int i, j, k;
    for (i = 0; i != 2; ++i)
        for (k = 0; k != 2; ++k)
        {
            a[i][k] = 0;
            for (j = 0; j != 2; ++j)
                a[i][k] += b[i][j] * c[j][k];
        }
}

static void mul122(double a[1][2], double b[1][2], double c[2][2])
{
    int i, j, k;
    for (i = 0; i != 1; ++i)
        for (k = 0; k != 2; ++k)
        {
            a[i][k] = 0;
            for (j = 0; j != 2; ++j)
                a[i][k] += b[i][j] * c[j][k];
        }
}

static void mul121(double a[1][1], double b[1][2], double c[2][1])
{
    int i, j, k;
    for (i = 0; i != 1; ++i)
        for (k = 0; k != 1; ++k)
        {
            a[i][k] = 0;
            for (j = 0; j != 2; ++j)
                a[i][k] += b[i][j] * c[j][k];
        }
}

static double det22(double a[2][2])
{
    return a[0][0] * a[1][1] - a[0][1] * a[1][0];
}

/* log functions
 *
 * These functions assist in working with values that might overflow if we
 * ever explicitly constructed their exponentials.
 */

/*
 * FIXME: when the C99 transition is complete, remove log1p declaration
 */

/*
 * declare the C99 function log1p to suppress warning
 *
 * log1p(x) = log(1 + x), but is more accurate for |x| << 1
 */
double log1p(double x);

/*
 * XLALSkymapLogSumExp(a, b) computes log(exp(a) + exp(b)) but will not
 * overflow for a or b > ~300
 *
 * For a > b, we use the identity
 *
 * log(exp(a) + exp(b)
 *     = log(exp(a) * (1 + exp(b) / exp(a)))
 *     = a + log(1 + exp(b - a))
 *     = a + log1p(exp(b - a))
 *
 * where b - a < 0 and exp(b - a) cannot overflow (though it may
 * underflow).
 *
 * And for a < b
 *
 * log(exp(a) + exp(b))
 *     = log(exp(b) * (exp(a) / exp(b) + 1))
 *     = b + log(exp(a - b) + 1)
 *     = b + log1p(exp(a - b))
 *
 * If neither a < b nor a > b, we either have equality or both values
 * are (the same) plus or minus infinity.  Forming (a - b) in the case of
 * infinities results in a NaN, so we must use a third expression
 *
 * log(exp(a) + exp(b))
 *     = log(exp(a) + exp(a))
 *     = log(2 * exp(a))
 *     = log(2) + a
 */

double XLALSkymapLogSumExp(double a, double b)
{
    return
        (a < b) ? (b + log1p(exp(a - b))) : (
	    (b < a) ? (a + log1p(exp(b - a))) : (a + log(2))
	    );
}

double XLALSkymapLogDifferenceExp(double a, double b)
{
    return a + log(1.0 - exp(b - a));
}

static double* findmax(double* begin, double* end)
{
    double* p;
    double* m;
    m = begin;
    for (p = begin; p != end; ++p)
	if (*m < *p)
	    m = p;
    return m;
}

static double logtotalexpwithmax(double* begin, double* end, double m)
{
    double t;
    double* p;
    t = 0;
    for (p = begin; p != end; ++p)
	t += exp(*p - m);
    return m + log(t);
}

double XLALSkymapLogTotalExp(double* begin, double* end)
{
    return logtotalexpwithmax(begin, end, *findmax(begin, end));
}

/* coordinate transformations */

void XLALSkymapCartesianFromSpherical(double a[3], double theta, double phi)
{
    set3(a, sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

void XLALSkymapSphericalFromCartesian(double a[2], double b[3])
{
    set2(a, acos(b[2]), atan2(b[1], b[0]));
}

/* geometrical properties of an interferometer site */

static double site_time(LALDetector* site, double direction[3])
{
    return -dot3(site->location, direction) / LAL_C_SI;
}

static void site_response(double f[2], LALDetector* site, double direction[3])
{
    double thetaphi[2];
    XLALSkymapSphericalFromCartesian(thetaphi, direction);
    XLALComputeDetAMResponse(&f[0], &f[1], site->response, thetaphi[1], LAL_PI_2 - thetaphi[0], 0, 0);
}

static void construct_hlv(LALDetector site[3])
{
    site[0] = lalCachedDetectors[LAL_LHO_4K_DETECTOR];
    site[1] = lalCachedDetectors[LAL_LLO_4K_DETECTOR];
    site[2] = lalCachedDetectors[LAL_VIRGO_DETECTOR];
}

/* properties of each pixel in the tiling */

static void construct_pixel(XLALSkymapPixelType* pixel)
{
    pixel->area = 0;
    set3(pixel->direction, 0, 0, 0);
}

/* properties of a network and a sample rate */

/*
 * These functions convert between integer time delays that index a
 * conceptual three-dimensional array, and a linear index corresponding
 * to how the array is actually packed in memory
 *
 * For a simple L x M x N array, the indices (i, j, k) map to a linear
   index
 *
 *     p = (M * N) * i + N * j + k
 *
 * To store pixels, we have an array with
 *
 *     L = 2                 (hemispheres)
 *     M = 2 * plan->hv + 1  (delays from -plan->hv to +plan->hv)
 *     N = 2 * plan->hl + 1  (delays from -plan->hl to +plan->hl)
 *
 * and indices that are offset delays (so that they are >= 0)
 *
 *     i = hemisphere
 *     j = hv + plan->hv
 *     k = hl + plan->hl
 *
 */

static int index_from_delays(XLALSkymapPlanType* plan, int hl, int hv, int hemisphere)
{
    return
	(hl + plan->hl) +
	(hv + plan->hv) * (plan->hl * 2 + 1) +
	(hemisphere   ) * (plan->hl * 2 + 1) * (plan->hv * 2 + 1);
}

static void delays_from_index(XLALSkymapPlanType *plan, int lal_index, int delays[3])
{
    delays[2] = lal_index / ((plan->hl * 2 + 1) * (plan->hv * 2 + 1));
    lal_index %= ((plan->hl * 2 + 1) * (plan->hv * 2 + 1));
    delays[1] = (lal_index / (plan->hl * 2 + 1)) - plan->hv;
    delays[0] = (lal_index % (plan->hl * 2 + 1)) - plan->hl;
}

/*
static int index_from_times(XLALSkymapPlanType* plan, int times[3], int hemisphere)
{
    return index_from_delays(plan, times[1] - times[0], times[2] - times[0], hemisphere);
}
*/

void XLALSkymapDelaysFromDirection(XLALSkymapPlanType* plan, int delays[3], double direction[3])
{
    /* compute the delays rounded to the nearest sample */
    delays[0] = (int) floor((
        site_time(&plan->site[1], direction) -
        site_time(&plan->site[0], direction)) *
        plan->sampleFrequency + 0.5);
    delays[1] = (int) floor((
        site_time(&plan->site[2], direction) -
        site_time(&plan->site[0], direction)) *
        plan->sampleFrequency + 0.5);
    delays[2] = dot3(plan->siteNormal, direction) < 0 ? 0 : 1;
}

int XLALSkymapIndexFromDirection(XLALSkymapPlanType* plan, double direction[3])
{
    int i[3];
    XLALSkymapDelaysFromDirection(plan, i, direction);
    return index_from_delays(plan, i[0], i[1], i[2]);
}

XLALSkymapPlanType* XLALSkymapConstructPlanMN(int sampleFrequency, int m, int n)
{
    static const char func[] = "XLALSkymapConstructPlanMN";
    XLALSkymapPlanType* plan;

    if (sampleFrequency <= 0)
    {
      XLALPrintError("%s(): invalid sample frequency %d\n", func, sampleFrequency);
      XLAL_ERROR_NULL(func, XLAL_EINVAL);
    }

    if(m <= 0 || n <= 0)
    {
        XLALPrintError("%s(): invalid raster size %d x %d\n", func, m, n);
        XLAL_ERROR_NULL(func, XLAL_EINVAL);
    }

    plan = (XLALSkymapPlanType*) XLALMalloc(sizeof(*plan));

    /* possible runtime error */
    if (!plan)
    {
      XLAL_ERROR_NULL(func, XLAL_EFUNC);
    }

    /* define the sample frequency */
    plan->sampleFrequency = sampleFrequency;
    plan->m = m;
    plan->n = n;

    /* set up the hanford-livingston-virgo network */
    construct_hlv(plan->site);

    {   /* compute the baselines */
        double v_hl[3], v_hv[3], v_lv[3];
        sub3(v_hl, plan->site[1].location, plan->site[0].location);
        sub3(v_hv, plan->site[2].location, plan->site[0].location);
        sub3(v_lv, plan->site[2].location, plan->site[1].location);
        div3(v_hl, v_hl, LAL_C_SI);
        div3(v_hv, v_hv, LAL_C_SI);
        div3(v_lv, v_lv, LAL_C_SI);
        cross(plan->siteNormal, v_hl, v_hv);
        /* compute the maximum delay rounded to nearest sample */
        plan->hl = (int) floor(nrm3(v_hl) * plan->sampleFrequency + 0.5);
        plan->hv = (int) floor(nrm3(v_hv) * plan->sampleFrequency + 0.5);
        plan->lv = (int) floor(nrm3(v_lv) * plan->sampleFrequency + 0.5);
    }

    plan->pixelCount = (plan->hl * 2 + 1) * (plan->hv * 2 + 1) * 2;
    plan->pixel = (XLALSkymapPixelType*) XLALMalloc(sizeof(*plan->pixel) * plan->pixelCount);

    /* possible runtime error */
    if (!plan->pixel)
    {
        XLALFree(plan);
        XLAL_ERROR_NULL(func, XLAL_EFUNC);
    }

    /* compute the pixel properties */
    {
        int i, j;
        double area = 0;

        /* initialize the pixels */
        for (i = 0; i != plan->pixelCount; ++i)
        {
            construct_pixel(&plan->pixel[i]);
        }

        /* scan over the sky */
        for (i = 0; i != m; ++i)
        {
            for (j = 0; j != n; ++j)
            {
                double direction[3];
                int k;
                /* compute theta and phi */
                double theta = acos(1. - (i + 0.5) * (2. / plan->m));
                double phi   = (j + 0.5) * (2. * pi / plan->n);
                /* compute direction */
                XLALSkymapCartesianFromSpherical(direction, theta, phi);
                /* determine the corresponding pixel */
                k = XLALSkymapIndexFromDirection(plan, direction);
                /* increase the total area */
                area += 1.0;
                /* increase the pixel area */
                plan->pixel[k].area += 1.0;
                /* accumulate the direction */
                add3(plan->pixel[k].direction, plan->pixel[k].direction, direction);
            }
        }

        /* compute the relative area, average direction, and antenna pattern for each physical pixel */
        for (i = 0; i != plan->pixelCount; ++i)
        {   /* is the direction physical? */
            if (plan->pixel[i].area > 0)
            {   /* compute the average direction */
                normalize3(plan->pixel[i].direction);
                /* normalize the area */
                plan->pixel[i].area /= area;
                for (j = 0; j != 3; ++j)
                {
                    site_response(plan->pixel[i].f[j], &plan->site[j], plan->pixel[i].direction);
                }
            }
            else
            {
                /* zero out the unused memory */
                set32(plan->pixel[i].f, 0, 0, 0, 0, 0, 0);
            }
        }
    }

    return plan;
}

void XLALSkymapDestroyPlan(XLALSkymapPlanType* plan)
{
    if (plan)
        XLALFree(plan->pixel);
    XLALFree(plan);
}

static void compute_kernel(XLALSkymapPlanType* plan, int lal_index, double sigma, double w[3], double kernel[3][3], double* log_normalization)
{
    double f[2][3];
    double wfsfw[2][2] = {{0, 0}, {0, 0}};
    double wfsfwa[2][2];
    double iwfsfwa[2][2];
    double fiwfsfwa[3][2];
    double wfsfwiwfsfwa[2][2];
    int i, j;

    transpose32(f, plan->pixel[lal_index].f);
    for (i = 0; i != 3; ++i)
    {
        wfsfw[0][0] += sq(f[0][i] * w[i]);
        wfsfw[0][1] += f[0][i] * f[1][i] * sq(w[i]);
        wfsfw[1][1] += sq(f[1][i] * w[i]);
    }
    wfsfw[1][0] = wfsfw[0][1];

    set22(wfsfwa, wfsfw[0][0] + 1/sq(sigma), wfsfw[0][1]    ,
    wfsfw[1][0]    , wfsfw[1][1] + 1/sq(sigma));
    inv22(iwfsfwa, wfsfwa);
    mul322(fiwfsfwa, plan->pixel[lal_index].f, iwfsfwa);
    mul323(kernel, fiwfsfwa, f);
    for (i = 0; i != 3; ++i)
    {
        for (j = 0; j != 3; ++j)
        {
            kernel[i][j] *= w[i] * w[j];
        }
    }

    mul222(wfsfwiwfsfwa, wfsfw, iwfsfwa);
    wfsfwiwfsfwa[0][0] -= 1;
    wfsfwiwfsfwa[1][1] -= 1;
    *log_normalization = log(det22(wfsfwiwfsfwa)) * 0.5;
}

static void compute_kernel2(XLALSkymapPlanType* plan, int lal_index, double sigma, double w[3], double kernel[2][2], double* log_normalization, double** x)
{
    double f[2][2];
    double wfsfw[2][2] = {{0, 0}, {0, 0}};
    double wfsfwa[2][2];
    double iwfsfwa[2][2];
    double fiwfsfwa[2][2];
    double wfsfwiwfsfwa[2][2];
    int i, j;
    double w2[2];

    /* transpose32(f, plan->pixel[lal_index].f); */
    if (x[0])
    {
        /* we have hanford data */
        f[0][0] = plan->pixel[lal_index].f[0][0];
        f[1][0] = plan->pixel[lal_index].f[0][1];
        w2[0] = w[0];
        if (x[1])
        {
            /* we have livingston data */
            f[0][1] = plan->pixel[lal_index].f[1][0];
            f[1][1] = plan->pixel[lal_index].f[1][1];
            w2[1] = w[1];
        }
        else
        {
            /* we must have virgo data */
            f[0][1] = plan->pixel[lal_index].f[2][0];
            f[1][1] = plan->pixel[lal_index].f[2][1];
            w2[1] = w[2];
        }
    }
    else
    {
        /* we must have livingston and virgo data */
        f[0][0] = plan->pixel[lal_index].f[1][0];
        f[1][0] = plan->pixel[lal_index].f[1][1];
        w2[0] = w[1];
        f[0][1] = plan->pixel[lal_index].f[2][0];
        f[1][1] = plan->pixel[lal_index].f[2][1];
        w2[1] = w[2];
    }

    for (i = 0; i != 2; ++i)
    {
        wfsfw[0][0] += sq(f[0][i] * w2[i]);
        wfsfw[0][1] += f[0][i] * f[1][i] * sq(w2[i]);
        wfsfw[1][1] += sq(f[1][i] * w2[i]);
    }
    wfsfw[1][0] = wfsfw[0][1];

    set22(wfsfwa, wfsfw[0][0] + 1/sq(sigma), wfsfw[0][1]    ,
    wfsfw[1][0]    , wfsfw[1][1] + 1/sq(sigma));
    inv22(iwfsfwa, wfsfwa);
    mul222(fiwfsfwa, plan->pixel[lal_index].f, iwfsfwa);
    mul222(kernel, fiwfsfwa, f);
    for (i = 0; i != 2; ++i)
    {
        for (j = 0; j != 2; ++j)
        {
            kernel[i][j] *= w2[i] * w2[j];
        }
    }

    mul222(wfsfwiwfsfwa, wfsfw, iwfsfwa);
    wfsfwiwfsfwa[0][0] -= 1;
    wfsfwiwfsfwa[1][1] -= 1;
    *log_normalization = log(det22(wfsfwiwfsfwa)) * 0.5;
}

static void compute_kernel1(XLALSkymapPlanType* plan, int lal_index, double sigma, double w[3], double kernel[1][1], double* log_normalization, double** x)
{
    double f[2][1];
    double wfsfw[2][2] = {{0, 0}, {0, 0}};
    double wfsfwa[2][2];
    double iwfsfwa[2][2];
    double fiwfsfwa[1][2];
    double wfsfwiwfsfwa[2][2];
    int i;

    if (x[0])
    {
        i = 0;
    }
    else
    {
        if (x[1])
        {
            i = 1;
        }
        else
        {
            i = 2;
        }
    }

    /* transpose32(f, plan->pixel[lal_index].f); */
    f[0][0] = plan->pixel[lal_index].f[i][0];
    f[1][0] = plan->pixel[lal_index].f[i][1];

    wfsfw[0][0] += sq(f[0][0] * w[i]);
    wfsfw[0][1] += f[0][0] * f[1][0] * sq(w[i]);
    wfsfw[1][1] += sq(f[1][0] * w[i]);

    wfsfw[1][0] = wfsfw[0][1];

    set22(wfsfwa, wfsfw[0][0] + 1/sq(sigma), wfsfw[0][1]    ,
    wfsfw[1][0]    , wfsfw[1][1] + 1/sq(sigma));
    inv22(iwfsfwa, wfsfwa);
    mul122(fiwfsfwa, plan->pixel[lal_index].f, iwfsfwa);
    mul121(kernel, fiwfsfwa, f);
    kernel[0][0] *= w[i] * w[i];

    mul222(wfsfwiwfsfwa, wfsfw, iwfsfwa);
    wfsfwiwfsfwa[0][0] -= 1;
    wfsfwiwfsfwa[1][1] -= 1;
    *log_normalization = log(det22(wfsfwiwfsfwa)) * 0.5;
}

int XLALSkymapGlitchHypothesis(XLALSkymapPlanType* plan, double *p, double sigma, double w[3], int begin[3], int end[3], double** x)
{
    /* static const char func[] = "XLALSkymapGlitchHypothesis"; */
    double* buffer;
    int i;

    /* allocate working memory */
    buffer = (double*) XLALMalloc(sizeof(double)* max(max(end[0], end[1]), end[2]));

    /* analyze each detector individually */
    for (i = 0; i != 3; ++i)
    {
        /* was data supplied for this detector? */
        if (x[i])
        {
            int t;
            double k;
            double *m;
            /* compute the kernel */
            k = 0.5 / (sq(w[i]) + 1.0 / sq(sigma));
            /* loop over all times */
            for (t = begin[i]; t != end[i]; ++t)
            {
                buffer[t] = k * (sq(x[i][t]) + sq(x[i + 3][t]));
            }
            /* find the maximum to prevent over or underflow when accumulating exponentials */
            m = findmax(buffer + begin[i], buffer + end[i]);
            /* marginalize over time */
            p[i] = logtotalexpwithmax(buffer + begin[i], buffer + end[i], *m);
            /* apply normalization */
            p[i] -= log(sq(sigma) * sq(w[i]) + 1.0);
        }
        else
        {
            /* if there is no data we draw no conclusion */
            p[i] = 0.0;
        }
    }
    XLALFree(buffer);
    return 0;
}

/* a new interface is required that instead produces a per-tile product */

/* return not a skymap but a grid of odds ratios for each direction individually */
/* to make a skymap requires correcting for unequal areas and unequal time ranges */

int XLALSkymapSignalHypothesis(XLALSkymapPlanType* plan, double* p, double sigma, double w[3], int begin[3], int end[3], double** x, int *counts, int *modes)
{
    int delay_limits[6];
    delay_limits[0] = -plan->hl;
    delay_limits[1] =  plan->hl;
    delay_limits[2] = -plan->hv;
    delay_limits[3] =  plan->hv;
    delay_limits[4] = -plan->lv - 1;
    delay_limits[5] =  plan->lv + 1;
    return XLALSkymapSignalHypothesisWithLimits(plan, p, sigma, w, begin, end, x, counts, modes, delay_limits);
}

int XLALSkymapSignalHypothesisWithLimits(XLALSkymapPlanType* plan, double* p, double sigma, double w[3], int begin[3], int end[3], double** x, int *counts, int *modes, int delay_limits[6])
{
    /* static const char func[] = "XLALSkymapSignalHypothesis"; */

    double* buffer;
    int hl;
    double total_normalization;

    /* double* times[3]; Arrival time posterior for each detector */

    total_normalization = 0;

    buffer = 0;
    if (x[0])
    {
        buffer = (double*) XLALMalloc(sizeof(double) * (end[0] - begin[0]));
    }
    else
    {
        if (x[1])
        {
            buffer = (double*) XLALMalloc(sizeof(double) * (end[1] - begin[1]));
        }
        else
        {
            if (x[2])
            {
                buffer = (double*) XLALMalloc(sizeof(double) * (end[2] - begin[2]));
            }
            else
            {
                return 1;
            }
        }
    }

    /* Initialize the arrival time posterior for each detector
    {
        int d;
        for (d = 0; d != 3; ++d)
        {
            if (x[d])
            {
                int i;
                times[d] = (double*) XLALMalloc(sizeof(double) * end[d]);
                for (i = 0; i != end[d]; ++i)
                {
                    times[d][i] = log(0);
                }
            }
        }
    }
     */

    /* loop over hanford-livingston delay */
    for (hl = -plan->hl; hl <= plan->hl; ++hl)
    {
        /* loop over hanford-virgo delay */
        int hv;
        for (hv = -plan->hv; hv <= plan->hv; ++hv)
        {
            /* loop over hemisphere */
            int hemisphere;
            for (hemisphere = 0; hemisphere != 2; ++hemisphere)
            {
                /* compute the index into the buffers from delays */
                int lal_index = index_from_delays(plan, hl, hv, hemisphere);
                /* (log) zero the output buffers */
                p[lal_index] = log(0);
                counts[lal_index] = 0;
                modes[lal_index] = 0;
                /* test if the delays are physical and relevant */
                if (
                    (plan->pixel[lal_index].area > 0)
                    &&
                    (hl >= delay_limits[0])
                    &&
                    (hl <= delay_limits[1])
                    &&
                    (hv >= delay_limits[2])
                    &&
                    (hv <= delay_limits[3])
                    &&
                    ((hv - hl) >= delay_limits[4])
                    &&
                    ((hv - hl) <= delay_limits[5])
                    )
                {
                    /* compute the begin and end times */
                    int b;
                    int e;

                    if (x[0])
                    { /* if we have hanford data */
                        b = begin[0];
                        e = end[0];
                        if (x[1])
                        { /* if we have livingston data */
                            if (b + hl < begin[1])
                            {   /* livingston constrains the begin time */
                                b = begin[1] - hl;
                            }
                            if (e + hl > end[1])
                            {   /* livingston constrains the end time */
                                e = end[1] - hl;
                            }
                        }
                        if (x[2])
                        { /* if we have virgo data */
                            if (b + hv < begin[2])
                            {   /* virgo constrains the begin time */
                                b = begin[2] - hv;
                            }
                            if (e + hv > end[2])
                            {   /* virgo constrains the end time */
                                e = end[2] - hv;
                            }
                        }
                    }
                    else
                    {
                        if (x[1])
                        { /* if we have livingston data */
                            b = begin[1] - hl;
                            e = end[1] -hl;
                            if (x[2])
                            { /* if we have virgo data */
                                if (b + hv < begin[2])
                                {   /* virgo constrains the begin time */
                                    b = begin[2] - hv;
                                }
                                if (e + hv > end[2])
                                {   /* virgo constrains the end time */
                                    e = end[2] - hv;
                                }
                            }
                        }
                        else
                        { /* if we have virgo data */
                            if (x[2])
                            {
                                b = begin[2] - hv;
                                e = end[2] - hv;
                            }
                            else
                            { /* we have no data */
                                b = 0;
                                e = 0;
                            }
                        }
                    }

                    /* test if there is enough data to analyze */
                    if (b < e)
                    {
                        double log_normalization;
                        double *hr, *lr, *vr, *hi, *li, *vi;
                        double* stop;
                        double* q;
                        double* m;

                        /* create offset pointers to simplify inner loop */
                        hr = x[0] + b     ; /* hanford real */
                        lr = x[1] + b + hl; /* livingston real */
                        vr = x[2] + b + hv; /* virgo real */
                        hi = x[3] + b     ; /* hanford imag */
                        li = x[4] + b + hl; /* livingston imag */
                        vi = x[5] + b + hv; /* virgo imag */

                        stop = x[0] + e; /* hanford real stop */

                        q = buffer;

                        if (x[0])
                        {
                            if (x[1])
                            {
                                if (x[2])
                                {
                                    /* hanford-livingston-virgo analysis */
                                    double kernel[3][3];
                                    /* compute the kernel */
                                    compute_kernel(plan, lal_index, sigma, w, kernel, &log_normalization);
                                    /* loop over arrival times */
                                    for (; hr != stop; ++hr, ++lr, ++vr, ++hi, ++li, ++vi, ++q)
                                    {
                                        *q = 0.5 * (
                                                kernel[0][0] * (sq(*hr) + sq(*hi)) +
                                                kernel[0][1] * (*hr * *lr + *hi * *li) * 2 +
                                                kernel[0][2] * (*hr * *vr + *hi * *vi) * 2 +
                                                kernel[1][1] * (sq(*lr) + sq(*li)) +
                                                kernel[1][2] * (*lr * *vr + *li * *vi) * 2 +
                                                kernel[2][2] * (sq(*vr) + sq(*vi))
                                                );
                                    } /* end loop over arrival times */
                                }
                                else
                                {
                                    /* hanford-livingston analysis */
                                    double kernel[2][2];
                                    /* compute the kernel */
                                    compute_kernel2(plan, lal_index, sigma, w, kernel, &log_normalization, x);
                                    /* loop over arrival times */
                                    for (; hr != stop; ++hr, ++lr, ++vr, ++hi, ++li, ++vi, ++q)
                                    {
                                        *q = 0.5 * (
                                                kernel[0][0] * (sq(*hr) + sq(*hi)) +
                                                kernel[0][1] * (*hr * *lr + *hi * *li) * 2 +
                                                kernel[1][1] * (sq(*lr) + sq(*li))
                                                );
                                    } /* end loop over arrival times */
                                }
                            }
                            else
                            {
                                if (x[2])
                                {
                                    /* hanford-virgo analysis */
                                    double kernel[2][2];
                                    /* compute the kernel */
                                    compute_kernel2(plan, lal_index, sigma, w, kernel, &log_normalization, x);
                                    /* loop over arrival times */
                                    for (; hr != stop; ++hr, ++lr, ++vr, ++hi, ++li, ++vi, ++q)
                                    {
                                        *q = 0.5 * (
                                                kernel[0][0] * (sq(*hr) + sq(*hi)) +
                                                kernel[0][1] * (*hr * *vr + *hi * *vi) * 2 +
                                                kernel[1][1] * (sq(*vr) + sq(*vi))
                                                );
                                    } /* end loop over arrival times */
                                }
                                else
                                {
                                    /* hanford analysis */
                                    double kernel[1][1];
                                    /* compute the kernel */
                                    compute_kernel1(plan, lal_index, sigma, w, kernel, &log_normalization, x);
                                    /* loop over arrival times */
                                    for (; hr != stop; ++hr, ++lr, ++vr, ++hi, ++li, ++vi, ++q)
                                    {
                                        *q = 0.5 * (
                                                kernel[0][0] * (sq(*hr) + sq(*hi))
                                                );
                                    } /* end loop over arrival times */
                                }
                            }
                        }
                        else
                        {
                            if (x[1])
                            {
                                if (x[2])
                                {
                                    /* livingston-virgo analysis */
                                    double kernel[2][2];
                                    /* compute the kernel */
                                    compute_kernel2(plan, lal_index, sigma, w, kernel, &log_normalization, x);
                                    /* loop over arrival times */
                                    for (; hr != stop; ++hr, ++lr, ++vr, ++hi, ++li, ++vi, ++q)
                                    {
                                        *q = 0.5 * (
                                                kernel[0][0] * (sq(*lr) + sq(*li)) +
                                                kernel[0][1] * (*lr * *vr + *li * *vi) * 2 +
                                                kernel[1][1] * (sq(*vr) + sq(*vi))
                                                );
                                    } /* end loop over arrival times */
                                }
                                else
                                {
                                    /* livingston analysis */
                                    double kernel[1][1];
                                    /* compute the kernel */
                                    compute_kernel1(plan, lal_index, sigma, w, kernel, &log_normalization, x);
                                    /* loop over arrival times */
                                    for (; hr != stop; ++hr, ++lr, ++vr, ++hi, ++li, ++vi, ++q)
                                    {
                                        *q = 0.5 * (
                                                kernel[0][0] * (sq(*lr) + sq(*li))
                                                );
                                    } /* end loop over arrival times */
                                }
                            }
                            else
                            {
                                if (x[2])
                                {
                                    /* virgo analysis */
                                    double kernel[1][1];
                                    /* compute the kernel */
                                    compute_kernel1(plan, lal_index, sigma, w, kernel, &log_normalization, x);
                                    /* loop over arrival times */
                                    for (; hr != stop; ++hr, ++lr, ++vr, ++hi, ++li, ++vi, ++q)
                                    {
                                        *q = 0.5 * (
                                                kernel[0][0] * (sq(*vr) + sq(*vi))
                                                );
                                    } /* end loop over arrival times */
                                }
                                else
                                {
                                    /* no data at all */
                                    return 1;
                                }
                            }
                        }

                        /* compute the plausibility for the direction */
                        m = findmax(buffer, buffer + e - b);
                        p[lal_index] =
                                logtotalexpwithmax(buffer, buffer + e - b, *m) +
                                log_normalization * 2 - log(e - b);
                        counts[lal_index] = e - b;
                        modes[lal_index] = b + (m - buffer);

                        /* Compute the arrival time posterior for each detector
                        {
                            int i;
                            for (i = b; i != e; ++i)
                            {
                                double np;
                                np = buffer[i - b] + log_normalization * 2 + log(plan->pixel[lal_index].area);
                                times[0][i] = XLALSkymapLogSumExp(times[0][i], np);
                                times[1][i + hl] = XLALSkymapLogSumExp(times[1][i + hl], np);
                                times[2][i + hv] = XLALSkymapLogSumExp(times[2][i + hv], np);
                            }
                        }
                        */
                    } /* end test if there is enough data to analyze */
                } /* end test if the delays are physical */
            } /* end loop over hemisphere */
        } /* end loop over hanford-virgo delay */
    } /* end loop over hanford-livingston delay */

    /* release working memory */
    XLALFree(buffer);
    return 0;
}

int XLALSkymapEllipticalHypothesis(XLALSkymapPlanType* plan, double* p, double sigma, double w[3], int begin[3], int end[3], double** x, int* bests)
{
    /* indicate that a detector has no data by x[i] == 0 */

    /* static const char func[] = "XLALSkymapEllipticalHypothesis"; */

    int *counts;
    int *modes;
    double c;
    int i;

    counts = (int *) XLALMalloc(sizeof(int) * plan->pixelCount);
    modes = (int *) XLALMalloc(sizeof(int) * plan->pixelCount);

    XLALSkymapSignalHypothesis(plan, p, sigma, w, begin, end, x, counts, modes);

    /* the prior of each pixel is the product of its area and the length
     * arrival times it represents */

    /* compute the normalization factor for the prior */
    c = 0;
    for (i = 0; i != plan->pixelCount; ++i)
    {
        if (plan->pixel[i].area > 0)
        {
            c += plan->pixel[i].area * counts[i];
        }
    }
    /* apply the prior */
    for (i = 0; i != plan->pixelCount; ++i)
    {
        if (plan->pixel[i].area > 0)
        {
            p[i] += log(plan->pixel[i].area * counts[i] / c);
        }
    }
    /* now the sum over the skymap is the posterior odds of a signal from
     * any direction */

    if (bests) {
        /* find the most likely arrival times */
        double *a;
        int j;
        int delays[3];
        a = findmax(p, p + plan->pixelCount);
        j = a - p;
        printf("i = %d\n", j);
        delays_from_index(plan, j, delays);
        printf("delays[] = { %d, %d, %d}\n", delays[0], delays[1], delays[2]);
        printf("index_from_delays(delays) = %d\n", index_from_delays(plan, delays[0], delays[1], delays[2]));
        bests[0] = modes[j];
        bests[1] = bests[0] + delays[0];
        bests[2] = bests[0] + delays[1];
        bests[3] = delays[2];
    }

    XLALFree(modes);
    XLALFree(counts);

    return 0;
}


void XLALSkymapSum(XLALSkymapPlanType* plan, double* a, const double* b, const double* c)
{
    int i;
    for (i = 0; i != plan->pixelCount; ++i)
    {   /* check to see if the pixel is valid */
        if (plan->pixel[i].area > 0)
        {   /* sum the log-represented values */
            a[i] = XLALSkymapLogSumExp(b[i], c[i]);
        }
    }
}

void XLALSkymapModeThetaPhi(XLALSkymapPlanType* plan, double* p, double thetaphi[2])
{
    double mode = 0;
    int mode_i = -1;
    int i;
    /* for each pixel */
    for (i = 0; i < plan->pixelCount; ++i)
    {   /* only consider valid pixels */
        if (plan->pixel[i].area > 0)
        {   /* if the greatest or first valid pixel */
            if ((mode_i == -1) || (p[i] - log(plan->pixel[i].area)) > mode)
            {
                mode   = p[i] - log(plan->pixel[i].area);
                mode_i = i;
            }
        }
    }
    XLALSkymapSphericalFromCartesian(thetaphi, plan->pixel[mode_i].direction);
}

int XLALSkymapRender(double* q, XLALSkymapPlanType* plan, double* p)
{
    static const char func[] = "XLALSkymapRender";
    int i, j;

    /* scan over the sky */
    for (i = 0; i != plan->m; ++i)
    {
        for (j = 0; j != plan->n; ++j)
        {
            double direction[3];
            int k;
            /* compute theta and phi */
            double theta = acos(1. - (i + 0.5) * (2. / plan->m));
            double phi   = (j + 0.5) * (2. * pi / plan->n);
            /* compute direction */
            XLALSkymapCartesianFromSpherical(direction, theta, phi);
            /* determine the corresponding pixel */
            k = XLALSkymapIndexFromDirection(plan, direction);

            if (plan->pixel[k].area > 0)
            {
                if (p[k] < -log(0))
                {
                    q[i + j * plan->m] = p[k];
                    /* apply area corrections */
                    q[i + j * plan->m] -= log(plan->pixel[k].area);
                }
                else
                {
                    XLALPrintError("%s(): attempted to render from a pixel with value +inf or nan to (%i, %j)\n", func, i, j);
                    XLAL_ERROR(func, XLAL_EINVAL);
                }
            }
            else
            {
                XLALPrintError("%s(): attempted to render from a pixel with zero area to (%i, %j)\n", func, i, j);
                XLAL_ERROR(func, XLAL_EINVAL);
            }
        }
    }

    return 0;
}

