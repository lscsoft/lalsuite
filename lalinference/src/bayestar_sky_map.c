/*                                           >y#
                                            ~'#o+
                                           '~~~md~
                '|+>#!~'':::::....        .~'~'cY#
            .+oy+>|##!~~~''':::......     ~:'':md! .
          #rcmory+>|#~''':::'::...::.::. :..'''Yr:...
        'coRRaamuyb>|!~'''::::':...........  .+n|.::..
       !maMMNMRYmuybb|!~'''':.........::::::: ro'..::..
      .cODDMYouuurub!':::...........:::~'.. |o>::...:..
      >BDNCYYmroyb>|#~:::::::::::::~':.:: :ob::::::::..
      uOCCNAa#'''''||':::.                :oy':::::::::.
    :rRDn!  :~::'y+::':  ... ...:::.     :ob':::::::::::.
   yMYy:   :>yooCY'.':.   .:'':......    ~u+~::::::::::::.
  >>:'. .~>yBDMo!.'': . .:'':.   .      >u|!:::::::::::::.
    ':'~|mYu#:'~'''. :.~:':...         yy>|~:::::::::::::..
    :!ydu>|!rDu::'. +'#~::!#'.~:     |r++>#':::::::::::::..
    mn>>>>>YNo:'': !# >'::::...  ..:cyb++>!:::::::::..:::...
    :ouooyodu:'': .!:.!:::.       yobbbb+>~::::::::....:....
     'cacumo~''' .'~ :~'.::.    :aybbbbbb>':::'~''::::....
      .mamd>'''. :~' :':'.:.   om>bbbyyyb>'.#b>|#~~~'':..
      .yYYo''': .:~' .'::'   .ny>+++byyoao!b+|||#!~~~''''''::.
      .#RUb:''. .:'' .:':   |a#|>>>>yBMdb #yb++b|':::::''':'::::::.
      .'CO!'''  .:'' .'    uu~##|+mMYy>+:|yyo+:::'::.         .::::::
      .:RB~''' ..::'.':   o>~!#uOOu>bby'|yB>.'::  '~!!!!!~':. ..  .::::
       :Rm''': ..:~:!:  'c~~+YNnbyyybb~'mr.':  !+yoy+>||!~'::.       :::.
      ..Oo''': .'' ~:  !+|BDCryuuuuub|#B!::  !rnYaocob|#!~'':.  ..    .::.
      . nB''': :  .'  |dNNduroomnddnuun::.  ydNAMMOary+>#~:.:::...      .:
       .uC~'''    :. yNRmmmadYUROMMBmm.:   bnNDDDMRBoy>|#~':....:.      .:
                 :' ymrmnYUROMAAAAMYn::. .!oYNDDMYmub|!~'::....:..     :
                 !'#booBRMMANDDDNNMO!:. !~#ooRNNAMMOOmuy+#!':::.......    :.
                .!'!#>ynCMNDDDDDNMRu.. '|:!raRMNAMOOdooy+|!~:::........   .:
                 : .'rdbcRMNNNNAMRB!:  |!:~bycmdYYBaoryy+|!~':::.::::.  ..
                 ..~|RMADnnONAMMRdy:. .>#::yyoroccruuybb>#!~'::':...::.
                  :'oMOMOYNMnybyuo!.  :>#::b+youuoyyy+>>|!~':.    :::::
                  ''YMCOYYNMOCCCRdoy##~~~: !b>bb+>>>||#~:..:::     ::::.
                  .:OMRCoRNAMOCROYYUdoy|>~:.~!!~!~~':...:'::::.   :::::.
                  ''oNOYyMNAMMMRYnory+|!!!:.....     ::.  :'::::::::::::
                 .:..uNabOAMMCOdcyb+|!~':::.          !!'.. :~:::::'''':.
                  .   +Y>nOORYauyy>!!'':....           !#~..  .~:''''''':.

****************  ____  _____  ______________________    ____     **************
***************  / __ )/   \ \/ / ____/ ___/_  __/   |  / __ \   ***************
**************  / __  / /| |\  / __/  \__ \ / / / /| | / /_/ /  ****************
*************  / /_/ / ___ |/ / /___ ___/ // / / ___ |/ _, _/  *****************
************  /_____/_/  |_/_/_____//____//_/ /_/  |_/_/ |_|  ******************
*/


/*
 * Copyright (C) 2013  Leo Singer
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

#include "bayestar_sky_map.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <lal/DetResponse.h>
#include <lal/InspiralInjectionParams.h>
#include <lal/LALSimulation.h>

#include <chealpix.h>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_test.h>

#include "logaddexp.h"

#ifndef _OPENMP
#define omp ignore
#endif


/* Storage for old GSL error handler. */
static gsl_error_handler_t *old_handler;


/* Custom, reentrant GSL error handler that simply prints the error message. */
static void
my_gsl_error (const char *reason, const char *file, int line, int gsl_errno)
{
    (void)gsl_errno;
    fprintf(stderr, "gsl: %s:%d: %s: %s\n", file, line, "ERROR", reason);
}


/* Custom error handler that ignores underflow errors. */
static void
ignore_underflow (const char *reason, const char *file, int line, int gsl_errno)
{
    if (gsl_errno != GSL_EUNDRFLW)
        old_handler(reason, file, line, gsl_errno);
}


/* Compute |z|^2. Hopefully a little faster than gsl_pow_2(cabs(z)), because no
 * square roots are necessary. */
static double cabs2(double complex z) {
    return gsl_pow_2(creal(z)) + gsl_pow_2(cimag(z));
}


/*
 * Catmull-Rom cubic spline interpolant of x(t) for regularly gridded
 * samples x_i(t_i), assuming:
 *
 *     t_0 = -1, x_0 = x[0],
 *     t_1 = 0,  x_1 = x[1],
 *     t_2 = 1,  x_2 = x[2],
 *     t_3 = 2,  x_3 = x[3].
 */
static double complex complex_catrom(
    double complex x0,
    double complex x1,
    double complex x2,
    double complex x3,
    double t
) {
    return x1
        + t*(-0.5*x0 + 0.5*x2
        + t*(x0 - 2.5*x1 + 2.*x2 - 0.5*x3
        + t*(-0.5*x0 + 1.5*x1 - 1.5*x2 + 0.5*x3)));
}


/*
 * Catmull-Rom cubic spline interpolant of x(t) for regularly gridded
 * samples x_i(t_i), assuming:
 *
 *     t_0 = -1, x_0 = x[0],
 *     t_1 = 0,  x_1 = x[1],
 *     t_2 = 1,  x_2 = x[2],
 *     t_3 = 2,  x_3 = x[3].
 *
 * I am careful to handle certain situations where some of
 * the samples are infinite or not-a-number:
 *
 *     * If t <= 0, then return x1.
 *     * If t >= 1, then return x2.
 *     * Otherwise, if either x0 or x3 are non-real, then fall back to
 *       linear interpolation between x1 and x2.
 *     * Otherwise, if x1 and/or x2 are infinite and both have the same
 *       then return infinity of that sign.
 *     * If x1 and x2 are infinities of different signs or if either are
 *       NaN, then return NaN.
 *     * Otherwise, if x0, x1, x2, and x3 are all finite, then return
 *       the standard Catmull-Rom formula.
 *
 * ***IMPORTANT NOTE*** I use the ISO C99 function isfinite() instead of
 * isinf(), the non-standard finite(), or gsl_finite(), because isfinite()
 * is reliably faster than the alternatives on Mac OS and Scientific Linux.
 *
 */
static double real_catrom(
    double x0,
    double x1,
    double x2,
    double x3,
    double t
) {
    double result;

    if (t <= 0)
        result = x1;
    else if (t >= 1)
        result = x2;
    else if (!(isfinite(x0 + x3)))
        result = x1 + (1 - t) * x2;
    else if (isfinite(x1) && isfinite(x2))
        result = x1
            + t*(-0.5*x0 + 0.5*x2
            + t*(x0 - 2.5*x1 + 2.*x2 - 0.5*x3
            + t*(-0.5*x0 + 1.5*x1 - 1.5*x2 + 0.5*x3)));
    else
        result = x1 + x2;

    return result;
}


/* Evaluate a complex time series using cubic spline interpolation, assuming
 * that the vector x gives the samples of the time series at times
 * 0, 1, ..., nsamples-1, and that the time series is given by the complex
 * conjugate at negative times. */
static double complex eval_acor(
    const double complex *x,
    size_t nsamples,
    double t
) {
    size_t i;
    double f;
    double complex y;

    /* Break |t| into integer and fractional parts. */
    {
        double dbl_i;
        f = modf(fabs(t), &dbl_i);
        i = dbl_i;
    }

    if (i == 0)
        y = complex_catrom(conj(x[1]), x[0], x[1], x[2], f);
    else if (i < nsamples - 2)
        y = complex_catrom(x[i-1], x[i], x[i+1], x[i+2], f);
    else
        y = 0;

    if (t < 0)
        y = conj(y);

    return y;
}


typedef struct {
    size_t size;
    double dx, xmin;
    double y[];
} cubic_interp;


static double cubic_interp_eval(const cubic_interp *interp, double x)
{
    double i;
    const double u = modf((x - interp->xmin) / interp->dx, &i);

    #define CLAMP(m) ((size_t) ((m) <= 0 ? 0 : ((m) >= interp->size - 1 ? interp->size - 1 : m)))
    const double y0 = interp->y[CLAMP(i - 1)];
    const double y1 = interp->y[CLAMP(i)];
    const double y2 = interp->y[CLAMP(i + 1)];
    const double y3 = interp->y[CLAMP(i + 2)];
    #undef CLAMP

    return real_catrom(y0, y1, y2, y3, u);
}


typedef struct {
    size_t xsize, ysize;
    double xmin, ymin, dx, dy;
    double z[];
} bicubic_interp;


static double bicubic_interp_eval(const bicubic_interp *interp, double x, double y)
{
    double i, j;
    const double u = modf((x - interp->xmin) / interp->dx, &i);
    const double v = modf((y - interp->ymin) / interp->dy, &j);
    double z[4];

    #define CLAMPX(m) ((size_t) ((m) <= 0 ? 0 : ((m) >= interp->xsize - 1 ? interp->xsize - 1 : m)))
    #define CLAMPY(m) ((size_t) ((m) <= 0 ? 0 : ((m) >= interp->ysize - 1 ? interp->ysize - 1 : m)))
    for (unsigned int k = 0; k < 4; k ++)
    {
        const size_t ii = CLAMPX(i + k - 1);
        const double z0 = interp->z[ii * interp->xsize + CLAMPY(j - 1)];
        const double z1 = interp->z[ii * interp->xsize + CLAMPY(j)];
        const double z2 = interp->z[ii * interp->xsize + CLAMPY(j + 1)];
        const double z3 = interp->z[ii * interp->xsize + CLAMPY(j + 2)];
        z[k] = real_catrom(z0, z1, z2, z3, v);
    }
    #undef CLAMPX
    #undef CLAMPY

    return real_catrom(z[0], z[1], z[2], z[3], u);
}


struct log_radial_integrator_t {
    bicubic_interp *region0;
    cubic_interp *region1;
    cubic_interp *region2;
    double xmax, ymax, vmax, r1, r2;
    int k;
};


typedef struct {
    double scale;
    double p;
    double b;
    int k;
} radial_integrand_params;


static double radial_integrand(double r, void *params)
{
    const radial_integrand_params *integrand_params = params;
    const double scale = integrand_params->scale;
    const double p = integrand_params->p;
    const double b = integrand_params->b;
    const int k = integrand_params->k;
    return gsl_sf_exp_mult(
        scale - gsl_pow_2(p / r - 0.5 * b / p),
        gsl_sf_bessel_I0_scaled(b / r) * gsl_pow_int(r, k));
}


static double log_radial_integrand(double r, void *params)
{
    const radial_integrand_params *integrand_params = params;
    const double scale = integrand_params->scale;
    const double p = integrand_params->p;
    const double b = integrand_params->b;
    const int k = integrand_params->k;
    return log(gsl_sf_bessel_I0_scaled(b / r) * gsl_pow_int(r, k))
        + scale - gsl_pow_2(p / r - 0.5 * b / p);
}


static double log_radial_integral(double r1, double r2, double p, double b, int k)
{
    radial_integrand_params params = {0, p, b, k};
    double breakpoints[5];
    unsigned char nbreakpoints = 0;
    double result, abserr, log_offset = -INFINITY;
    int ret;

    if (b != 0) {
        /* Calculate the approximate distance at which the integrand attains a
         * maximum (middle) and a fraction eta of the maximum (left and right).
         * This neglects the scaled Bessel function factors and the power-law
         * distance prior. It assumes that the likelihood is approximately of
         * the form
         *
         *    -p^2/r^2 + B/r.
         *
         * Then the middle breakpoint occurs at 1/r = -B/2A, and the left and
         * right breakpoints occur when
         *
         *   A/r^2 + B/r = log(eta) - B^2/4A.
         */

        static const double eta = 0.01;
        const double middle = 2 * gsl_pow_2(p) / b;
        const double left = 1 / (1 / middle + sqrt(-log(eta)) / p);
        const double right = 1 / (1 / middle - sqrt(-log(eta)) / p);

        /* Use whichever of the middle, left, and right points lie within the
         * integration limits as initial subdivisions for the adaptive
         * integrator. */

        breakpoints[nbreakpoints++] = r1;
        if(left > breakpoints[nbreakpoints-1] && left < r2)
            breakpoints[nbreakpoints++] = left;
        if(middle > breakpoints[nbreakpoints-1] && middle < r2)
            breakpoints[nbreakpoints++] = middle;
        if(right > breakpoints[nbreakpoints-1] && right < r2)
            breakpoints[nbreakpoints++] = right;
        breakpoints[nbreakpoints++] = r2;
    } else {
        /* Inner breakpoints are undefined because b = 0. */
        breakpoints[nbreakpoints++] = r1;
        breakpoints[nbreakpoints++] = r2;
    }

    /* Re-scale the integrand so that the maximum value at any of the
     * breakpoints is 1. Note that the initial value of the constant term
     * is overwritten. */

    for (unsigned char i = 0; i < nbreakpoints; i++)
    {
        double new_log_offset = log_radial_integrand(breakpoints[i], &params);
        if (new_log_offset > log_offset)
            log_offset = new_log_offset;
    }

    /* If the largest value of the log integrand was -INFINITY, then the
     * integrand is 0 everywhere. Set log_offset to 0, because subtracting
     * -INFINITY would make the integrand infinite. */
    if (log_offset == -INFINITY)
        log_offset = 0;

    params.scale = -log_offset;

    {
        /* Maximum number of subdivisions for adaptive integration. */
        static const size_t n = 64;

        /* Allocate workspace on stack. Hopefully, a little bit faster than
         * using the heap in multi-threaded code. */

        double alist[n];
        double blist[n];
        double rlist[n];
        double elist[n];
        size_t order[n];
        size_t level[n];
        gsl_integration_workspace workspace = {
            .alist = alist,
            .blist = blist,
            .rlist = rlist,
            .elist = elist,
            .order = order,
            .level = level,
            .limit = n
        };

        /* Set up integrand data structure. */
        const gsl_function func = {radial_integrand, &params};

        /* Perform adaptive Gaussian quadrature. */
        ret = gsl_integration_qagp(&func, breakpoints, nbreakpoints,
            DBL_MIN, 1e-8, n, &workspace, &result, &abserr);

        /* FIXME: do we care to keep the error estimate around? */
    }

    /* FIXME: do something with ret */
    (void)ret;

    /* Done! */
    return log_offset + log(result);
}


static const size_t default_log_radial_integrator_size = 400;


log_radial_integrator *log_radial_integrator_init(double r1, double r2, int k, double pmax, size_t size)
{
    if (size <= 1)
        GSL_ERROR_NULL("size must be > 1", GSL_EINVAL);

    const double alpha = 4;
    const double p0 = 0.5 * (k >= 0 ? r2 : r1);
    const double xmax = log(pmax);
    const double x0 = GSL_MIN_DBL(log(p0), xmax);
    const double xmin = x0 - (1 + M_SQRT2) * alpha;
    const double ymax = x0 + alpha;
    const double ymin = 2 * x0 - M_SQRT2 * alpha - xmax;
    const size_t len = size * size;
    const double d = (xmax - xmin) / (size - 1); /* dx = dy = du */
    const double umin = - (1 + M_SQRT1_2) * alpha;
    const double vmax = x0 - M_SQRT1_2 * alpha;
    /* const double umax = xmax - vmax; */ /* unused */

    log_radial_integrator *integrator = malloc(sizeof(log_radial_integrator));
    if (!integrator)
        GSL_ERROR_NULL("not enough memory to allocate integrator", GSL_ENOMEM);
    integrator->region0 = malloc(sizeof(bicubic_interp) + len * sizeof(double));
    if (!integrator->region0)
    {
        free(integrator);
        GSL_ERROR_NULL("not enough memory to allocate integrator", GSL_ENOMEM);
    }
    integrator->region1 = malloc(sizeof(cubic_interp) + size * sizeof(double));
    if (!integrator->region1)
    {
        free(integrator->region0);
        free(integrator);
        GSL_ERROR_NULL("not enough memory to allocate integrator", GSL_ENOMEM);
    }
    integrator->region2 = malloc(sizeof(cubic_interp) + size * sizeof(double));
    if (!integrator->region2)
    {
        free(integrator->region0);
        free(integrator->region1);
        free(integrator);
        GSL_ERROR_NULL("not enough memory to allocate integrator", GSL_ENOMEM);
    }

    integrator->region0->xsize = integrator->region0->ysize = size;
    integrator->region0->xmin = xmin;
    integrator->region0->ymin = ymin;
    integrator->region0->dx = integrator->region0->dy = d;

    old_handler = gsl_set_error_handler(ignore_underflow);

    #pragma omp parallel for
    for (size_t i = 0; i < len; i ++)
    {
        const size_t ix = i / size;
        const size_t iy = i % size;
        const double x = xmin + ix * d;
        const double y = ymin + iy * d;
        const double p = exp(x);
        const double r0 = exp(y);
        const double b = 2 * gsl_pow_2(p) / r0;
        /* Note: using this where p > r0; could reduce evaluations by half */
        integrator->region0->z[i] = log_radial_integral(r1, r2, p, b, k);
    }

    gsl_set_error_handler(old_handler);

    integrator->region1->size = size;
    integrator->region1->xmin = xmin;
    integrator->region1->dx = d;

    for (size_t i = 0; i < size; i ++)
    {
        integrator->region1->y[i] = integrator->region0->z[i * size + (size - 1)];
    }

    integrator->region2->size = size;
    integrator->region2->xmin = umin;
    integrator->region2->dx = d;

    for (size_t i = 0; i < size; i ++)
    {
        integrator->region2->y[i] = integrator->region0->z[i * size + (size - 1 - i)];
    }

    integrator->xmax = xmax;
    integrator->ymax = ymax;
    integrator->vmax = vmax;
    integrator->r1 = r1;
    integrator->r2 = r2;
    integrator->k = k;
    return integrator;
}


void log_radial_integrator_free(log_radial_integrator *integrator)
{
    if (integrator)
    {
        free(integrator->region0);
        integrator->region0 = NULL;
        free(integrator->region1);
        integrator->region1 = NULL;
        free(integrator->region2);
        integrator->region2 = NULL;
        free(integrator);
    }
}


double log_radial_integrator_eval(const log_radial_integrator *integrator, double p, double b)
{
    const double r0 = 2 * gsl_pow_2(p) / b;
    const double x = log(p);
    const double y = log(r0);
    double result;
    assert(x <= integrator->xmax);

    if (p == 0) {
        /* note: p2 == 0 implies b == 0 */
        assert(b == 0);
        int k1 = integrator->k + 1;

        if (k1 == 0)
        {
            result = log(log(integrator->r2 / integrator->r1));
        } else {
            result = log((gsl_pow_int(integrator->r2, k1) - gsl_pow_int(integrator->r1, k1)) / k1);
        }
    } else {
        if (y >= integrator->ymax) {
            result = cubic_interp_eval(integrator->region1, x);
        } else {
            const double v = 0.5 * (x + y);
            if (v <= integrator->vmax)
            {
                const double u = 0.5 * (x - y);
                result = cubic_interp_eval(integrator->region2, u);
            } else {
                result = bicubic_interp_eval(integrator->region0, x, y);
            }
        }
        result += gsl_pow_2(0.5 * b / p);
    }

    return result;
}


/* Find error in time of arrival. */
static void toa_errors(
    double *dt,
    double theta,
    double phi,
    double gmst,
    int nifos,
    const double **locs, /* Input: detector position. */
    const double *toas /* Input: time of arrival. */
) {
    /* Convert to Cartesian coordinates. */
    double n[3];
    ang2vec(theta, phi - gmst, n);

    for (int i = 0; i < nifos; i ++)
    {
        double dot = 0;
        for (int j = 0; j < 3; j ++)
        {
            dot += locs[i][j] * n[j];
        }
        dt[i] = toas[i] + dot / LAL_C_SI;
    }
}


/* Compute antenna factors from the detector response tensor and source
 * sky location, and return as a complex number F_plus + i F_cross. */
static double complex complex_antenna_factor(
    const float response[3][3],
    double ra,
    double dec,
    double gmst
) {
    double complex F;
    XLALComputeDetAMResponse(
        (double *)&F,     /* Type-punned real part */
        1 + (double *)&F, /* Type-punned imag part */
        response, ra, dec, 0, gmst);
    return F;
}


/* Expression for complex amplitude on arrival (without 1/distance factor) */
static double complex signal_amplitude_model(
    double complex F,               /* Antenna factor */
    double complex exp_i_twopsi,    /* e^(i*2*psi), for polarization angle psi */
    double u,                       /* cos(inclination) */
    double u2                       /* cos^2(inclination */
) {
    const double complex tmp = F * conj(exp_i_twopsi);
    return 0.5 * (1 + u2) * creal(tmp) + I * u * cimag(tmp);
}


static const unsigned int nglfixed = 10;
static const unsigned int ntwopsi = 10;


static double complex exp_i(double phi) {
    return cos(phi) + I * sin(phi);
}


/* Data structure to store a pixel in an adaptively refined sky map. */
typedef struct {
    double log4p;           /* Logarithm base 4 of probability */
    unsigned char order;    /* HEALPix resolution order */
    unsigned long ipix;     /* HEALPix nested pixel index */
} adaptive_sky_map_pixel;


/* Data structure to store an adaptively refined sky map. */
typedef struct {
    unsigned char max_order;
    size_t len;
    adaptive_sky_map_pixel pixels[];
} adaptive_sky_map;


/* Compare two pixels by contained probability. */
static int adaptive_sky_map_pixel_compare(const void *a, const void *b)
{
    int ret;
    const adaptive_sky_map_pixel *apix = a;
    const adaptive_sky_map_pixel *bpix = b;

    const double delta_log4p = apix->log4p - bpix->log4p;
    const char delta_order = apix->order - bpix->order;

    if (delta_log4p < delta_order)
        ret = -1;
    else if (delta_log4p > delta_order)
        ret = 1;
    else
        ret = 0;

    return ret;
}


/* Sort pixels by ascending contained probability. */
static void adaptive_sky_map_sort(adaptive_sky_map *map)
{
    qsort(map->pixels, map->len, sizeof(map->pixels[0]),
        adaptive_sky_map_pixel_compare);
}


static void *realloc_or_free(void *ptr, size_t size)
{
    void *new_ptr = realloc(ptr, size);
    if (!new_ptr)
    {
        free(ptr);
        GSL_ERROR_NULL("not enough memory to resize array", GSL_ENOMEM);
    }
    return new_ptr;
}


/* Subdivide the final last_n pixels of an adaptively refined sky map. */
static adaptive_sky_map *adaptive_sky_map_refine(
    adaptive_sky_map *map,
    size_t last_n
) {
    assert(last_n <= map->len);

    /* New length: adding 4*last_n new pixels, removing last_n old pixels. */
    const size_t new_len = map->len + 3 * last_n;
    const size_t new_size = sizeof(adaptive_sky_map)
        + new_len * sizeof(adaptive_sky_map_pixel);

    map = realloc_or_free(map, new_size);
    if (map)
    {
        for (size_t i = 0; i < last_n; i ++)
        {
            adaptive_sky_map_pixel *const old_pixel
                = &map->pixels[map->len - i - 1];
            const unsigned char order = 1 + old_pixel->order;
            if (order > map->max_order)
                map->max_order = order;
            const unsigned long ipix = 4 * old_pixel->ipix;
            for (unsigned char j = 0; j < 4; j ++)
            {
                adaptive_sky_map_pixel *const new_pixel
                    = &map->pixels[new_len - (4 * i + j) - 1];
                new_pixel->log4p = old_pixel->log4p;
                new_pixel->order = order;
                new_pixel->ipix = j + ipix;
            }
        }
        map->len = new_len;
    }
    return map;
}


static adaptive_sky_map *adaptive_sky_map_alloc(unsigned char order)
{
    unsigned long nside = (unsigned long)1 << order;
    unsigned long npix = nside2npix(nside);
    size_t size = sizeof(adaptive_sky_map)
        + npix * sizeof(adaptive_sky_map_pixel);

    adaptive_sky_map *map = malloc(size);
    if (!map)
        GSL_ERROR_NULL("not enough memory to allocate sky map", GSL_ENOMEM);

    map->len = npix;
    map->max_order = order;
    for (unsigned long ipix = 0; ipix < npix; ipix ++)
    {
        map->pixels[ipix].log4p = GSL_NAN;
        map->pixels[ipix].order = order;
        map->pixels[ipix].ipix = ipix;
    }
    return map;
}


static const double M_LN4 = 2 * M_LN2;
static const double M_1_LN4 = 0.5 / M_LN2;


static double *adaptive_sky_map_rasterize(adaptive_sky_map *map, long *out_npix)
{
    const unsigned char order = map->max_order;
    const unsigned long nside = (unsigned long)1 << order;
    const long npix = nside2npix(nside);

    double *P = malloc(npix * sizeof(double));
    if (!P)
        GSL_ERROR_NULL("not enough memory to allocate image", GSL_ENOMEM);

    double norm = 0;
    const double max_log4p = map->pixels[map->len - 1].log4p;

    /* Rescale so that log(max) = 0, and convert from log base 4 to
     * natural logarithm. */
    for (ssize_t i = (ssize_t)map->len - 1; i >= 0; i --)
    {
        map->pixels[i].log4p -= max_log4p;
        map->pixels[i].log4p *= M_LN4;
    }

    for (ssize_t i = (ssize_t)map->len - 1; i >= 0; i --)
    {
        const unsigned long reps = (unsigned long)1 << 2 * (order - map->pixels[i].order);
        const double dP = gsl_sf_exp_mult(map->pixels[i].log4p, reps);
        if (dP <= 0)
            break; /* We have reached underflow. */
        norm += dP;
    }
    norm = 1 / norm;
    for (ssize_t i = (ssize_t)map->len - 1; i >= 0; i --)
    {
        const double value = gsl_sf_exp_mult(map->pixels[i].log4p, norm);
        const unsigned long base_ipix = map->pixels[i].ipix;
        const unsigned long reps = (unsigned long)1 << 2 * (order - map->pixels[i].order);
        const unsigned long start = base_ipix * reps;
        const unsigned long stop = (base_ipix + 1) * reps;
        for (unsigned long ipix_nest = start; ipix_nest < stop; ipix_nest ++)
            P[ipix_nest] = value;
    }
    *out_npix = npix;

    return P;
}


double *bayestar_sky_map_toa_phoa_snr(
    long *inout_npix,
    /* Prior */
    double min_distance,            /* Minimum distance */
    double max_distance,            /* Maximum distance */
    int prior_distance_power,       /* Power of distance in prior */
    /* Detector network */
    double gmst,                    /* GMST (rad) */
    unsigned int nifos,             /* Number of detectors */
    unsigned long nsamples,         /* Length of autocorrelation sequence */
    double sample_rate,             /* Sample rate in seconds */
    const double complex **acors,   /* Autocorrelation sequences */
    const float (**responses)[3],   /* Detector responses */
    const double **locations,       /* Barycentered Cartesian geographic detector positions (m) */
    const double *horizons,         /* SNR=1 horizon distances for each detector */
    /* Observations */
    const double *toas,             /* Arrival time differences relative to network barycenter (s) */
    const double *phoas,            /* Phases on arrival */
    const double *snrs              /* SNRs */
) {
    double complex exp_i_phoas[nifos];
    for (unsigned int iifo = 0; iifo < nifos; iifo ++)
        exp_i_phoas[iifo] = exp_i(phoas[iifo]);

    log_radial_integrator *integrator;
    {
        double pmax = 0;
        for (unsigned int iifo = 0; iifo < nifos; iifo ++)
        {
            pmax += gsl_pow_2(horizons[iifo] / max_distance);
        }
        pmax = sqrt(0.5 * pmax);
        integrator = log_radial_integrator_init(
            min_distance, max_distance, prior_distance_power, pmax,
            default_log_radial_integrator_size);
        if (!integrator)
            return NULL;
    }

    static const unsigned char order0 = 4;
    unsigned char level = order0;
    adaptive_sky_map *map = adaptive_sky_map_alloc(order0);
    if (!map)
    {
        free(integrator);
        return NULL;
    }
    const unsigned long npix0 = map->len;

    /* Look up Gauss-Legendre quadrature rule for integral over cos(i). */
    gsl_integration_glfixed_table *gltable
        = gsl_integration_glfixed_table_alloc(nglfixed);

    /* Don't bother checking the return value. GSL has static, precomputed
     * values for certain orders, and for the order I have picked it will
     * return a pointer to one of these. See:
     *
     * http://git.savannah.gnu.org/cgit/gsl.git/tree/integration/glfixed.c
     */
    assert(gltable);
    assert(gltable->precomputed); /* We don't have to free it. */

    while (1)
    {
        /* Use our own error handler while in parallel section to avoid
         * concurrent calls to the GSL error handler, which if provided by the
         * user may not be threadsafe. */
        old_handler = gsl_set_error_handler(my_gsl_error);

        #pragma omp parallel for
        for (unsigned long i = 0; i < npix0; i ++)
        {
            adaptive_sky_map_pixel *const pixel = &map->pixels[map->len - npix0 + i];
            double complex F[nifos];
            double dt[nifos];
            double accum = -INFINITY;

            {
                double theta, phi;
                const unsigned long nside = (unsigned long)1 << pixel->order;
                pix2ang_nest(nside, pixel->ipix, &theta, &phi);

                /* Look up antenna factors */
                for (unsigned int iifo = 0; iifo < nifos; iifo++)
                    F[iifo] = complex_antenna_factor(
                        responses[iifo], phi, M_PI_2-theta, gmst) * horizons[iifo];

                toa_errors(dt, theta, phi, gmst, nifos, locations, toas);
            }

            /* Integrate over 2*psi */
            for (unsigned int itwopsi = 0; itwopsi < ntwopsi; itwopsi++)
            {
                const double twopsi = (2 * M_PI / ntwopsi) * itwopsi;
                const double complex exp_i_twopsi = exp_i(twopsi);
                double accum1 = -INFINITY;

                /* Integrate over u from -1 to 1. */
                for (unsigned int iu = 0; iu < nglfixed; iu++)
                {
                    double u, weight;
                    double accum2 = -INFINITY;
                    {
                        /* Look up Gauss-Legendre abscissa and weight. */
                        int ret = gsl_integration_glfixed_point(
                            -1, 1, iu, &u, &weight, gltable);

                        /* Don't bother checking return value; the only
                         * possible failure is in index bounds checking. */
                        assert(ret == GSL_SUCCESS);
                    }
                    const double u2 = gsl_pow_2(u);
                    double complex z_times_r[nifos];
                    double p2 = 0;

                    for (unsigned int iifo = 0; iifo < nifos; iifo ++)
                    {
                        p2 += cabs2(
                            z_times_r[iifo] = signal_amplitude_model(
                                F[iifo], exp_i_twopsi, u, u2));
                    }
                    p2 *= 0.5;
                    double p = sqrt(p2);

                    for (long isample = 1 - (long)nsamples;
                        isample < (long)nsamples; isample++)
                    {
                        double b;
                        {
                            double complex I0arg_complex_times_r = 0;
                            for (unsigned int iifo = 0; iifo < nifos; iifo ++)
                            {
                                I0arg_complex_times_r += snrs[iifo]
                                    * conj(z_times_r[iifo]) * exp_i_phoas[iifo]
                                    * eval_acor(acors[iifo], nsamples,
                                        dt[iifo] * sample_rate + isample);
                            }
                            b = cabs(I0arg_complex_times_r);
                        }

                        double result = log_radial_integrator_eval(
                            integrator, p, b);
                        accum2 = logaddexp(accum2, result);
                    }

                    accum1 = logaddexp(accum1, accum2 + log(weight));
                }

                accum = logaddexp(accum, accum1);
            }

            /* Record logarithm base 4 of posterior. */
            pixel->log4p = M_1_LN4 * accum;
        }

        /* Restore old error handler. */
        gsl_set_error_handler(old_handler);

        /* Sort pixels by ascending posterior probability. */
        adaptive_sky_map_sort(map);

        /* If we have reached order=11 (nside=2048), stop. */
        if (level++ >= 11)
            break;

        /* Adaptively refine the pixels that contain the most probability. */
        map = adaptive_sky_map_refine(map, npix0 / 4);
        if (!map)
            return NULL;
    }

    log_radial_integrator_free(integrator);

    /* Set error handler to ignore underflow errors, but invoke the user's
     * error handler for all other errors. */
    old_handler = gsl_set_error_handler(ignore_underflow);

    /* Flatten sky map to an image. */
    double *P = adaptive_sky_map_rasterize(map, inout_npix);
    free(map);

    /* Restore old error handler. */
    gsl_set_error_handler(old_handler);

    /* Done! */
    return P;
}


double bayestar_log_likelihood_toa_snr(
    /* Parameters */
    double ra,                      /* Right ascension (rad) */
    double sin_dec,                 /* Sin(declination) */
    double distance,                /* Distance */
    double u,                       /* Cos(inclination) */
    double twopsi,                  /* Twice polarization angle (rad) */
    double t,                       /* Barycentered arrival time (s) */
    /* Detector network */
    double gmst,                    /* GMST (rad) */
    unsigned int nifos,             /* Number of detectors */
    unsigned long nsamples,         /* Length of autocorrelation sequence */
    double sample_rate,             /* Sample rate in seconds */
    const double complex **acors,   /* Autocorrelation sequences */
    const float (**responses)[3],   /* Detector responses */
    const double **locations,       /* Barycentered Cartesian geographic detector positions (m) */
    const double *horizons,         /* SNR=1 horizon distances for each detector */
    /* Observations */
    const double *toas,             /* Arrival time differences relative to network barycenter (s) */
    const double *snrs              /* SNRs */
) {
    const double dec = asin(sin_dec);
    const double u2 = gsl_pow_2(u);
    const double complex exp_i_twopsi = exp_i(twopsi);
    const double one_by_r = 1 / distance;

    /* Compute time of arrival errors */
    double dt[nifos];
    toa_errors(dt, M_PI_2 - dec, ra, gmst, nifos, locations, toas);
    for (unsigned int iifo = 0; iifo < nifos; iifo++)
        dt[iifo] -= t;

    double A = 0, B = 0, product = 1;

    /* Loop over detectors */
    for (unsigned int iifo = 0; iifo < nifos; iifo++)
    {
        const double complex F = complex_antenna_factor(
            responses[iifo], ra, dec, gmst) * horizons[iifo];
        const double complex z_times_r =
            signal_amplitude_model(F, exp_i_twopsi, u, u2);
        const double rho2_times_r2 = cabs2(z_times_r);
        const double acor2 =
            cabs2(eval_acor(acors[iifo], nsamples, dt[iifo] * sample_rate));
        const double i0arg_times_r = snrs[iifo] * sqrt(rho2_times_r2 * acor2);

        A += rho2_times_r2;
        B += i0arg_times_r;
        product *= gsl_sf_bessel_I0_scaled(i0arg_times_r * one_by_r);
    }
    A *= -0.5;

    return (A * one_by_r + B) * one_by_r + log(product);
}


double bayestar_log_likelihood_toa_phoa_snr(
    /* Parameters */
    double ra,                      /* Right ascension (rad) */
    double sin_dec,                 /* Sin(declination) */
    double distance,                /* Distance */
    double u,                       /* Cos(inclination) */
    double twopsi,                  /* Twice polarization angle (rad) */
    double t,                       /* Barycentered arrival time (s) */
    /* Detector network */
    double gmst,                    /* GMST (rad) */
    unsigned int nifos,             /* Number of detectors */
    unsigned long nsamples,         /* Length of autocorrelation sequence */
    double sample_rate,             /* Sample rate in seconds */
    const double complex **acors,   /* Autocorrelation sequences */
    const float (**responses)[3],   /* Detector responses */
    const double **locations,       /* Barycentered Cartesian geographic detector positions (m) */
    const double *horizons,         /* SNR=1 horizon distances for each detector */
    /* Observations */
    const double *toas,             /* Arrival time differences relative to network barycenter (s) */
    const double *phoas,            /* Phases on arrival */
    const double *snrs              /* SNRs */
) {
    const double dec = asin(sin_dec);
    const double u2 = gsl_pow_2(u);
    const double complex exp_i_twopsi = exp_i(twopsi);
    const double one_by_r = 1 / distance;

    /* Compute time of arrival errors */
    double dt[nifos];
    toa_errors(dt, M_PI_2 - dec, ra, gmst, nifos, locations, toas);
    for (unsigned int iifo = 0; iifo < nifos; iifo++)
        dt[iifo] -= t;

    double complex i0arg_complex_times_r = 0;
    double A = 0;

    /* Loop over detectors */
    for (unsigned int iifo = 0; iifo < nifos; iifo++)
    {
        const double complex F = complex_antenna_factor(
            responses[iifo], ra, dec, gmst) * horizons[iifo];

        const double complex z_times_r =
             signal_amplitude_model(F, exp_i_twopsi, u, u2);
        const double complex zhat = snrs[iifo] * exp_i(phoas[iifo]);

        i0arg_complex_times_r += zhat * conj(z_times_r)
            * eval_acor(acors[iifo], nsamples, dt[iifo] * sample_rate);
        A += cabs2(z_times_r);
    }
    A *= -0.5;

    const double i0arg_times_r = cabs(i0arg_complex_times_r);

    return (A * one_by_r + i0arg_times_r) * one_by_r
        + log(gsl_sf_bessel_I0_scaled(i0arg_times_r * one_by_r));
}


/*
 * Unit tests
 */


static void test_cabs2(double complex z)
{
    double result = cabs2(z);
    double expected = gsl_pow_2(cabs(z));
    gsl_test_abs(result, expected, 2 * GSL_DBL_EPSILON,
        "testing cabs2(%g + %g j)", creal(z), cimag(z));
}


static void test_complex_catrom(void)
{
    for (double t = 0; t <= 1; t += 0.01)
    {
        double complex result = complex_catrom(0, 0, 0, 0, t);
        double complex expected = 0;
        gsl_test_abs(creal(result), creal(expected), 0,
            "testing complex Catmull-rom interpolant for zero input");
        gsl_test_abs(cimag(result), cimag(expected), 0,
            "testing complex Catmull-rom interpolant for zero input");
    }

    for (double t = 0; t <= 1; t += 0.01)
    {
        double complex result = complex_catrom(1, 1, 1, 1, t);
        double complex expected = 1;
        gsl_test_abs(creal(result), creal(expected), 0,
            "testing complex Catmull-rom interpolant for unit input");
        gsl_test_abs(cimag(result), cimag(expected), 0,
            "testing complex Catmull-rom interpolant for unit input");
    }

    for (double t = 0; t <= 1; t += 0.01)
    {
        double complex result = complex_catrom(1.0j, 1.0j, 1.0j, 1.0j, t);
        double complex expected = 1.0j;
        gsl_test_abs(creal(result), creal(expected), 0,
            "testing complex Catmull-rom interpolant for unit imaginary input");
        gsl_test_abs(cimag(result), cimag(expected), 1,
            "testing complex Catmull-rom interpolant for unit imaginary input");
    }

    for (double t = 0; t <= 1; t += 0.01)
    {
        double complex result = complex_catrom(1.0+1.0j, 1.0+1.0j, 1.0+1.0j, 1.0+1.0j, t);
        double complex expected = 1.0+1.0j;
        gsl_test_abs(creal(result), creal(expected), 0,
            "testing complex Catmull-rom interpolant for unit real + imaginary input");
        gsl_test_abs(cimag(result), cimag(expected), 0,
            "testing complex Catmull-rom interpolant for unit real + imaginary input");
    }

    for (double t = 0; t <= 1; t += 0.01)
    {
        double complex result = complex_catrom(1, 0, 1, 4, t);
        double complex expected = gsl_pow_2(t);
        gsl_test_abs(creal(result), creal(expected), 0,
            "testing complex Catmull-rom interpolant for quadratic real input");
        gsl_test_abs(cimag(result), cimag(expected), 0,
            "testing complex Catmull-rom interpolant for quadratic real input");
    }
}


static void test_real_catrom(void)
{
    for (double t = 0; t <= 1; t += 0.01)
    {
        double result = real_catrom(0, 0, 0, 0, t);
        double expected = 0;
        gsl_test_abs(result, expected, 0,
            "testing Catmull-rom interpolant for zero input");
    }

    for (double t = 0; t <= 1; t += 0.01)
    {
        double result = real_catrom(1, 1, 1, 1, t);
        double expected = 1;
        gsl_test_abs(result, expected, 0,
            "testing Catmull-rom interpolant for unit input");
    }

    for (double t = 0; t <= 1; t += 0.01)
    {
        double result = real_catrom(1, 0, 1, 4, t);
        double expected = gsl_pow_2(t);
        gsl_test_abs(result, expected, 0,
            "testing Catmull-rom interpolant for quadratic input");
    }

    for (double t = 0; t <= 1; t += 0.01)
    {
        double result = real_catrom(
            GSL_POSINF, GSL_POSINF, GSL_POSINF, GSL_POSINF, t);
        double expected = GSL_POSINF;
        gsl_test_abs(result, expected, 0,
            "testing Catmull-rom interpolant for +inf input");
    }

    for (double t = 0; t <= 1; t += 0.01)
    {
        double result = real_catrom(
            0, GSL_POSINF, GSL_POSINF, GSL_POSINF, t);
        double expected = GSL_POSINF;
        gsl_test_abs(result, expected, 0,
            "testing Catmull-rom interpolant for +inf input");
    }

    for (double t = 0; t <= 1; t += 0.01)
    {
        double result = real_catrom(
            GSL_POSINF, GSL_POSINF, GSL_POSINF, 0, t);
        double expected = GSL_POSINF;
        gsl_test_abs(result, expected, 0,
            "testing Catmull-rom interpolant for +inf input");
    }

    for (double t = 0; t <= 1; t += 0.01)
    {
        double result = real_catrom(
            0, GSL_POSINF, GSL_POSINF, 0, t);
        double expected = GSL_POSINF;
        gsl_test_abs(result, expected, 0,
            "testing Catmull-rom interpolant for +inf input");
    }

    for (double t = 0.01; t <= 1; t += 0.01)
    {
        double result = real_catrom(
            0, 0, GSL_POSINF, 0, t);
        double expected = GSL_POSINF;
        gsl_test_abs(result, expected, 0,
            "testing Catmull-rom interpolant for +inf input");
    }

    {
        double result = real_catrom(
            0, GSL_NEGINF, GSL_POSINF, 0, 1);
        double expected = GSL_POSINF;
        gsl_test_abs(result, expected, 0,
            "testing Catmull-rom interpolant for +inf input");
    }

    {
        double result = real_catrom(
            0, GSL_POSINF, GSL_NEGINF, 0, 0);
        double expected = GSL_POSINF;
        gsl_test_abs(result, expected, 0,
            "testing Catmull-rom interpolant for +inf input");
    }

    for (double t = 0; t <= 1; t += 0.01)
    {
        double result = real_catrom(
            0, GSL_NEGINF, GSL_NEGINF, GSL_NEGINF, t);
        double expected = GSL_NEGINF;
        gsl_test_abs(result, expected, 0,
            "testing Catmull-rom interpolant for -inf input");
    }

    for (double t = 0; t <= 1; t += 0.01)
    {
        double result = real_catrom(
            GSL_NEGINF, GSL_NEGINF, GSL_NEGINF, 0, t);
        double expected = GSL_NEGINF;
        gsl_test_abs(result, expected, 0,
            "testing Catmull-rom interpolant for -inf input");
    }

    for (double t = 0; t <= 1; t += 0.01)
    {
        double result = real_catrom(
            0, GSL_NEGINF, GSL_NEGINF, 0, t);
        double expected = GSL_NEGINF;
        gsl_test_abs(result, expected, 0,
            "testing Catmull-rom interpolant for -inf input");
    }

    for (double t = 0.01; t <= 1; t += 0.01)
    {
        double result = real_catrom(
            0, 0, GSL_NEGINF, 0, t);
        double expected = GSL_NEGINF;
        gsl_test_abs(result, expected, 0,
            "testing Catmull-rom interpolant for -inf input");
    }

    {
        double result = real_catrom(
            0, GSL_NEGINF, GSL_POSINF, 0, 0);
        double expected = GSL_NEGINF;
        gsl_test_abs(result, expected, 0,
            "testing Catmull-rom interpolant for -inf input");
    }

    {
        double result = real_catrom(
            0, GSL_POSINF, GSL_NEGINF, 0, 1);
        double expected = GSL_NEGINF;
        gsl_test_abs(result, expected, 0,
            "testing Catmull-rom interpolant for -inf input");
    }

    for (double t = 0.01; t < 1; t += 0.01)
    {
        double result = real_catrom(
            0, GSL_NEGINF, GSL_POSINF, 0, t);
        double expected = GSL_NAN;
        gsl_test_abs(result, expected, 0,
            "testing Catmull-rom interpolant for indeterminate input");
    }
}


static void test_eval_acor(void)
{
    size_t nsamples = 64;
    double complex x[nsamples];

    /* Populate data with samples of x(t) = t^2 + t j */
    for (size_t i = 0; i < nsamples; i ++)
        x[i] = gsl_pow_2(i) + i * 1.0j;

    for (double t = -nsamples; t <= nsamples; t += 0.1)
    {
        double result = eval_acor(x, nsamples, t);
        double expected = (fabs(t) < nsamples) ? (gsl_pow_2(t) + t*1.0j) : 0;
        gsl_test_abs(creal(result), creal(expected), 0,
            "testing real part of eval_acor(%g) for x(t) = t^2 + t j", t);
        gsl_test_abs(cimag(result), cimag(expected), 0,
            "testing imaginary part of eval_acor(%g) for x(t) = t^2 + t j", t);
    }
}


static void test_signal_amplitude_model(
    double ra,
    double dec,
    double inclination,
    double polarization,
    unsigned long epoch_nanos,
    const char *instrument
) {
    const LALDetector *detector = XLALDetectorPrefixToLALDetector(
        instrument);
    LIGOTimeGPS epoch;
    double gmst;
    XLALINT8NSToGPS(&epoch, epoch_nanos);
    gmst = XLALGreenwichMeanSiderealTime(&epoch);

    double complex result;
    {
        double complex F;
        const double complex exp_i_twopsi = exp_i(2 * polarization);
        const double u = cos(inclination);
        const double u2 = gsl_pow_2(u);
        XLALComputeDetAMResponse(
            (double *)&F,     /* Type-punned real part */
            1 + (double *)&F, /* Type-punned imag part */
            detector->response, ra, dec, 0, gmst);

        result = signal_amplitude_model(F, exp_i_twopsi, u, u2);
    }

    float abs_expected;
    {
        SimInspiralTable sim_inspiral;
        memset(&sim_inspiral, 0, sizeof(sim_inspiral));
        sim_inspiral.distance = 1;
        sim_inspiral.longitude = ra;
        sim_inspiral.latitude = dec;
        sim_inspiral.inclination = inclination;
        sim_inspiral.polarization = polarization;
        sim_inspiral.geocent_end_time = epoch;
        sim_inspiral.end_time_gmst = gmst;

        XLALInspiralSiteTimeAndDist(
            &sim_inspiral, detector, &epoch, &abs_expected);
        abs_expected = 1 / abs_expected;
    }

    double complex expected;
    {
        double Fplus, Fcross;
        const double u = cos(inclination);
        const double u2 = gsl_pow_2(u);
        XLALComputeDetAMResponse(
            &Fplus, &Fcross, detector->response, ra, dec, polarization, gmst);
        expected = 0.5 * (1 + u2) * Fplus + I * u * Fcross;
    }

    /* Test to nearly float (32-bit) precision because
     * XLALInspiralSiteTimeAndDist returns result as float. */
    gsl_test_abs(cabs(result), abs_expected, 1.5 * GSL_FLT_EPSILON,
        "testing abs(signal_amplitude_model(ra=%0.03f, dec=%0.03f, "
        "inc=%0.03f, pol=%0.03f, epoch_nanos=%lu, detector=%s))",
        ra, dec, inclination, polarization, epoch_nanos, instrument);

    gsl_test_abs(creal(result), creal(expected), 2.5 * GSL_DBL_EPSILON,
        "testing re(signal_amplitude_model(ra=%0.03f, dec=%0.03f, "
        "inc=%0.03f, pol=%0.03f, epoch_nanos=%lu, detector=%s))",
        ra, dec, inclination, polarization, epoch_nanos, instrument);

    gsl_test_abs(cimag(result), cimag(expected), 2.5 * GSL_DBL_EPSILON,
        "testing im(signal_amplitude_model(ra=%0.03f, dec=%0.03f, "
        "inc=%0.03f, pol=%0.03f, epoch_nanos=%lu, detector=%s))",
        ra, dec, inclination, polarization, epoch_nanos, instrument);
}


static void test_log_radial_integral(
    double expected, double tol, double r1, double r2, double p2, double b, int k)
{
    const double p = sqrt(p2);
    log_radial_integrator *integrator = log_radial_integrator_init(
        r1, r2, k, p + 0.5, default_log_radial_integrator_size);

    gsl_test(!integrator, "testing that integrator object is non-NULL");
    if (integrator)
    {
        const double result = log_radial_integrator_eval(integrator, p, b);

        gsl_test_rel(
            result, expected, tol,
            "testing toa_phoa_snr_log_radial_integral("
            "r1=%g, r2=%g, p2=%g, b=%g, k=%d)", r1, r2, p2, b, k);
    }

    log_radial_integrator_free(integrator);
}


int bayestar_test(void)
{
    for (double re = -1; re < 1; re += 0.1)
        for (double im = -1; im < 1; im += 0.1)
            test_cabs2(re + im * 1.0j);

    test_complex_catrom();
    test_real_catrom();
    test_eval_acor();

    for (double ra = -M_PI; ra <= M_PI; ra += 0.1 * M_PI)
        for (double dec = -M_PI_2; dec <= M_PI_2; dec += 0.1 * M_PI)
            for (double inc = -M_PI; inc <= M_PI; inc += 0.1 * M_PI)
                for (double pol = -M_PI; pol <= M_PI; pol += 0.1 * M_PI)
                    for (unsigned long epoch = 1000000000000000000;
                         epoch <= 1000086400000000000;
                         epoch += 3600000000000)
                        test_signal_amplitude_model(ra, dec, inc, pol, epoch, "H1");

    /* Tests of radial integrand with p2=0, b=0. */
    test_log_radial_integral(0, 0, 0, 1, 0, 0, 0);
    test_log_radial_integral(0, 0, exp(1), exp(2), 0, 0, -1);
    test_log_radial_integral(log(63), 0, 3, 6, 0, 0, 2);
    /* Te integrand with p2>0, b=0 (from Mathematica). */
    test_log_radial_integral(-0.480238, 1e-3, 1, 2, 1, 0, 0);
    test_log_radial_integral(0.432919, 1e-3, 1, 2, 1, 0, 2);
    test_log_radial_integral(-2.76076, 1e-3, 0, 1, 1, 0, 2);
    test_log_radial_integral(61.07118, 1e-3, 0, 1e9, 1, 0, 2);
    test_log_radial_integral(-112.23053, 5e-2, 0, 0.1, 1, 0, 2);
    /* Note: this test underflows, so we test that the log is -inf. */
    /* test_log_radial_integral(-1.00004e6, 1e-8, 0, 1e-3, 1, 0, 2); */
    test_log_radial_integral(-INFINITY, 1e-3, 0, 1e-3, 1, 0, 2);

    /* Tests of radial integrand with p2>0, b>0 with ML peak outside
     * of integration limits (true values from Mathematica NIntegrate). */
    test_log_radial_integral(2.94548, 1e-4, 0, 4, 1, 1, 2);
    test_log_radial_integral(2.94545, 1e-4, 0.5, 4, 1, 1, 2);
    test_log_radial_integral(2.94085, 1e-4, 1, 4, 1, 1, 2);
    /* Tests of radial integrand with p2>0, b>0 with ML peak outside
     * of integration limits (true values from Mathematica NIntegrate). */
    test_log_radial_integral(-2.43264, 1e-5, 0, 1, 1, 1, 2);
    test_log_radial_integral(-2.43808, 1e-5, 0.5, 1, 1, 1, 2);
    test_log_radial_integral(-0.707038, 1e-5, 1, 1.5, 1, 1, 2);

    {
        const double r1 = 0.0, r2 = 0.25, pmax = 1.0;
        const int k = 2;
        const double tol = 1e-5;
        log_radial_integrator *integrator = log_radial_integrator_init(
            r1, r2, k, pmax, default_log_radial_integrator_size);
        old_handler = gsl_set_error_handler(ignore_underflow);
        for (double p = 0.01; p <= pmax; p += 0.01)
        {
            for (double b = 0.0; b <= 2 * pmax; b += 0.01)
            {
                const double r0 = 2 * gsl_pow_2(p) / b;
                const double x = log(p);
                const double y = log(r0);
                const double expected = exp(log_radial_integral(r1, r2, p, b, k));
                const double result = exp(log_radial_integrator_eval(integrator, p, b) - gsl_pow_2(0.5 * b / p));
                gsl_test_abs(
                    result, expected, tol, "testing log_radial_integrator_eval("
                    "r1=%g, r2=%g, p=%g, b=%g, k=%d, x=%g, y=%g)", r1, r2, p, b, k, x, y);
            }
        }
        gsl_set_error_handler(old_handler);
        log_radial_integrator_free(integrator);
    }

    return gsl_test_summary();
}
