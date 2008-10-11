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

static void mul3(double a[3], double b[3], double c)
{
    a[0] = b[0] * c;
    a[1] = b[1] * c;
    a[2] = b[2] * c;
}

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
    a[0][1] = -b[1][0];
    a[1][0] = -b[0][1];
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

static double det22(double a[2][2])
{
    return a[0][0] * a[1][1] - a[0][1] * a[1][0];
}

/* coordinate transformations */

static void cartesian_from_spherical(double a[3], double theta, double phi)
{
    set3(a, sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

static void spherical_from_cartesian(double a[2], double b[3])
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
    spherical_from_cartesian(thetaphi, direction);
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

static int index_from_times(XLALSkymapPlanType* plan, int times[3], int h)
{
    int i[3];
    i[0] = times[1] - times[0];
    i[1] = times[2] - times[0];
    i[2] = h;
    return 
        (i[0] + plan->hl) + 
        (i[1] + plan->hv) * (plan->hl * 2 + 1) +
        (i[2]           ) * (plan->hl * 2 + 1) * (plan->hv * 2 + 1);
}

static int index_from_direction(XLALSkymapPlanType* plan, double direction[3])
{
    int i[3];
    /* compute the delays rounded to the nearest sample */
    i[0] = (int) floor((
        site_time(&plan->site[1], direction) - 
        site_time(&plan->site[0], direction)) * 
        plan->sampleFrequency + 0.5);
    i[1] = (int) floor((
        site_time(&plan->site[2], direction) - 
        site_time(&plan->site[0], direction)) * 
        plan->sampleFrequency + 0.5);
    i[2] = dot3(plan->siteNormal, direction) < 0 ? 0 : 1;
    return 
        (i[0] + plan->hl) + 
        (i[1] + plan->hv) * (plan->hl * 2 + 1) +
        (i[2]               ) * (plan->hl * 2 + 1) * (plan->hv * 2 + 1);
}

XLALSkymapPlanType* XLALSkymapConstructPlan(int sampleFrequency)
{
    static const char func[] = "XLALSkymapConstructPlan";
    XLALSkymapPlanType* plan;

    if (sampleFrequency <= 0)
    {
      XLALPrintError("%s(): invalid sample frequency %d\n", func, sampleFrequency);
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

    /* set up the hanford-livingston-virgo network */
    construct_hlv(plan->site);
    
    {   /* compute the baselines */
        double v_hl[3], v_hv[3];
        sub3(v_hl, plan->site[1].location, plan->site[0].location);
        sub3(v_hv, plan->site[2].location, plan->site[0].location);
        div3(v_hl, v_hl, LAL_C_SI);
        div3(v_hv, v_hv, LAL_C_SI);
        cross(plan->siteNormal, v_hl, v_hv);
        /* compute the maximum delay rounded to nearest sample */ 
        plan->hl = (int) floor(nrm3(v_hl) * plan->sampleFrequency + 0.5);
        plan->hv = (int) floor(nrm3(v_hv) * plan->sampleFrequency + 0.5);
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
        int m = 1024;
        int n = 2048;
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
                double theta = (i + 0.5) * (pi / m);
                double phi   = (j + 0.5) * (2 * pi / n);
                /* compute direction */
                cartesian_from_spherical(direction, theta, phi);
                /* determine the corresponding pixel */
                k = index_from_direction(plan, direction);
                /* increase the total area */
                area += sin(theta);
                /* increase the pixel area */
                plan->pixel[k].area += sin(theta);
                /* accumulate the direction */
                mul3(direction, direction, sin(theta));
                add3(plan->pixel[k].direction, plan->pixel[k].direction, direction);
            }
        }

        /* compute the relative area, average direction, and antenna pattern for each physical pixel */
        for (i = 0; i != plan->pixelCount; ++i)
        {   /* is the direction physical? */
            if (plan->pixel[i].area > 0)
            {   /* compute the average direction */
                div3(plan->pixel[i].direction, plan->pixel[i].direction, plan->pixel[i].area);
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

static void compute_kernel(XLALSkymapPlanType* plan, int index, double sigma, double w[3], double kernel[3][3], double* log_normalization)
{
    double f[2][3];
    double wfsfw[2][2] = {{0, 0}, {0, 0}};
    double wfsfwa[2][2];
    double iwfsfwa[2][2];
    double fiwfsfwa[3][2];
    double wfsfwiwfsfwa[2][2];
    int i, j;
    
    transpose32(f, plan->pixel[index].f);
    for (i = 0; i != 3; ++i)
    {
        wfsfw[0][0] += sq(plan->pixel[index].f[i][0] * w[i]);
        wfsfw[0][1] += plan->pixel[index].f[i][0] * plan->pixel[index].f[i][1] * sq(w[i]);
        wfsfw[1][1] += sq(plan->pixel[index].f[i][1] * w[i]);
    }
    wfsfw[1][0] = wfsfw[0][1];
    
    set22(wfsfwa, wfsfw[0][0] + 1/sq(sigma), wfsfw[0][1]    ,
    wfsfw[1][0]    , wfsfw[1][1] + 1/sq(sigma));
    inv22(iwfsfwa, wfsfwa);
    mul322(fiwfsfwa, plan->pixel[index].f, iwfsfwa);
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
    *log_normalization = log(0.5 * det22(wfsfwiwfsfwa));
}

int XLALSkymapAnalyzeElliptical(double* p, XLALSkymapPlanType* plan, double sigma, double w[3], int n, double** x)
{
    static const char func[] = "XLALSkymapAnalyzeElliptical";
    int hl, hv, hemisphere;
    double* buffer;

    if (n <= 0)
        XLAL_ERROR(func, XLAL_EINVAL);

    /* initialize storage */
    buffer = (double*) XLALMalloc(sizeof(*buffer) * n);    
    
    if (!buffer) /* handle out of memory */
    {            /* no resources to free */
        XLAL_ERROR(func, XLAL_EFUNC);
    }    
    
    for (hl = -plan->hl; hl <= plan->hl; ++hl)
    {
        for (hv = -plan->hv; hv <= plan->hv; ++hv)
        {
            int t[3];
            t[0] = 0;
            t[1] = hl;
            t[2] = hv;
            for (hemisphere = 0; hemisphere != 2; ++hemisphere)
            {
                int index = index_from_times(plan, t, hemisphere);
                if (plan->pixel[index].area > 0)
                {
                    int start = max(0, max(-hl, -hv));
                    int stop  = n - max(0, max(hl, hv));
                    if (stop > start)
                    {
                        /* compute the necessary terms */
                        double log_normalization;
                        double kernel[3][3];
                        
                        int tau;
                        double buffer_max = 0;

                        compute_kernel(plan, index, sigma, w, kernel, &log_normalization);
                        
                        for (tau = start; tau < stop; ++tau)
                        {
                            /* compute x . K . x */
                            buffer[tau] = 0.5 * (
                                kernel[0][0] * (sq(x[0][tau]) + sq(x[3][tau])) + 
                                kernel[0][1] * (x[0][tau] * x[1][tau + hl] + x[3][tau] * x[4][tau + hl]) * 2 +
                                kernel[0][2] * (x[0][tau] * x[2][tau + hv] + x[3][tau] * x[5][tau + hv]) * 2 +
                                kernel[1][1] * (sq(x[1][tau + hl]) + sq(x[4][tau + hl])) + 
                                kernel[1][2] * (x[1][tau + hl] * x[2][tau + hv] + x[4][tau + hl] * x[5][tau + hv]) * 2 +
                                kernel[2][2] * (sq(x[2][tau + hv]) + sq(x[5][tau + hv]))
                                );
                            /* track the maximum */
                            buffer_max = (tau == start) ? buffer[tau] : max(buffer[tau], buffer_max);
                        }
                        /* guarding against overflows, compute the sum */
                        for (tau = start; tau < stop; ++tau)
                        {
                            p[index] += exp(buffer[tau] - buffer_max);
                        }
                        /* return to log space */
                        p[index] = log(p[index]) + buffer_max;
                        /* normalize */
                        p[index] += 2 * log_normalization;
                        /* apply directional prior due to pixel geometry */
                        p[index] += log(plan->pixel[index].area);
                    }
                }
            }
        }
    }
    
    XLALFree(buffer);
    
    return 0;
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
    spherical_from_cartesian(thetaphi, plan->pixel[mode_i].direction);
}

int XLALSkymapRenderEquirectangular(int m, int n, double* q, XLALSkymapPlanType* plan, double* p)
{
    static const char func[] = "XLALSkymapRenderEquirectangular";
    int i, j;

    if(m <= 0 || n <= 0)
        XLAL_ERROR(func, XLAL_EINVAL);

    /* scan over the sky */
    for (i = 0; i != m; ++i)
    {
        for (j = 0; j != n; ++j)
        {
            double direction[3];
            int k;
            /* compute theta and phi */
            double theta = (i + 0.5) * (pi / m);
            double phi   = (j + 0.5) * (2 * pi / n);
            /* compute direction */
            cartesian_from_spherical(direction, theta, phi);
            /* determine the corresponding pixel */
            k = index_from_direction(plan, direction);

            q[i + j * m] = p[k];
            /* apply area corrections */
            q[i + j * m] += log(sin(theta) / plan->pixel[k].area);
        }
    }

    return 0;
}

int XLALSkymapRenderEqualArea(int m, int n, double* q, XLALSkymapPlanType* plan, double* p)
{
    static const char func[] = "XLALSkymapRenderEqualArea";
    int i, j;

    if(m <= 0 || n <= 0)
        XLAL_ERROR(func, XLAL_EINVAL);

    /* scan over the sky */
    for (i = 0; i != m; ++i)
    {
        for (j = 0; j != n; ++j)
        {
            double direction[3];
            int k;
            /* compute theta and phi */
            double theta = acos(1. - (i + 0.5) * (2. / m));
            double phi   = (j + 0.5) * (2. * pi / n);
            /* compute direction */
            cartesian_from_spherical(direction, theta, phi);
            /* determine the corresponding pixel */
            k = index_from_direction(plan, direction);

            if (plan->pixel[k].area > 0)
            {
                q[i + j * m] = p[k];
                /* apply area corrections */
                q[i + j * m] -= log(plan->pixel[k].area);
            }
            else
            {
                q[i + j * m] = log(0);
            }
        }
    }

    return 0;
}


int XLALSkymapRenderMollweide(int m, int n, double* q, XLALSkymapPlanType* plan, double* p)
{
    static const char func[] = "XLALSkymapRenderMollweide";
    int i, j;

    if(m <= 0 || n <= 0)
        XLAL_ERROR(func, XLAL_EINVAL);

    /* scan over the sky */
    for (i = 0; i != m; ++i)
    {
        for (j = 0; j != n; ++j)
        {
            double theta = asin(((i + 0.5) / m) * 2 - 1);
            double lambda = pi*(((j + 0.5) / n) * 2 - 1)/cos(theta);
            if ((lambda >= -pi) && (lambda <= pi))
            {
                double phi = asin(((2*theta)+sin(2*theta))/pi);
                double direction[3];
                /* compute direction */
                cartesian_from_spherical(direction, pi/2 - phi, lambda);
                /* determine the corresponding pixel */
                q[i + j * m] = p[index_from_direction(plan, direction)];
                /* undo area prior since this is an equiarea projection */
                q[i + j * m] -= log(plan->pixel[index_from_direction(plan, direction)].area);
            }
        }
    }

    return 0;
}

