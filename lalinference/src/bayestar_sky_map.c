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
#include <gsl/gsl_sort_vector_double.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_version.h> /* FIXME: remove, depend on GSL >= 1.15 */

#include "logaddexp.h"


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


/* For integers x and y, compute floor(x/y) using integer arithmetic. */
static int divfloor(int x, int y)
{
    if (y < 0)
    {
        x = -x;
        y = -y;
    }

    if (x >= 0)
        return x / y;
    else
        return (x - y + 1) / y;
}


/* Raise a real number x to a half-integer power twon/2.
 * This method uses repeated multiplication and, if necessary, one square root,
 * which is hopefully more efficient than pow(). */
static double pow_halfint(double x, int twon)
{
    double result = gsl_pow_int(x, divfloor(twon, 2));

    if (GSL_IS_ODD(twon))
        result *= sqrt(x);

    return result;
}


/* Exponential integral for half-integers n/2, scaled by e^x.
 * Avoids underflow or overflow by using asymptotic approximations. */
static double halfinteger_En_scaled(int twon, double x)
{
    if (x > 0.8 * GSL_LOG_DBL_MAX) {
        return 1 / x;
    } else if (x < 2 * GSL_SQRT_DBL_EPSILON) {
        if (twon < 2) {
            return gsl_sf_gamma(1 - 0.5 * twon) * pow_halfint(x, twon - 2);
        } else if (twon == 2) {
            return -log(M_SQRTPI * x);
        } else /* n > 2 */ {
            return 1 / (0.5 * twon - 1);
        }
    } else {
        if (GSL_IS_EVEN(twon) && twon >= 0)
        {
            return gsl_sf_expint_En_scaled(twon / 2, x);
        } else {
            return gsl_sf_exp_mult(x, gsl_sf_gamma_inc(1 - 0.5 * twon, x)
                * pow_halfint(x, twon - 2));
        }
    }
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
static double complex catrom(
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
        y = catrom(conj(x[1]), x[0], x[1], x[2], f);
    else if (i < nsamples - 2)
        y = catrom(x[i-1], x[i], x[i+1], x[i+2], f);
    else
        y = 0;

    if (t < 0)
        y = conj(y);

    return y;
}


/* Return indices of sorted pixels from greatest to smallest. */
static gsl_permutation *get_pixel_ranks(long npix, double *P)
{
    gsl_permutation *pix_perm = gsl_permutation_alloc(npix);
    if (pix_perm)
    {
        gsl_vector_view P_vector = gsl_vector_view_array(P, npix);
        gsl_sort_vector_index(pix_perm, &P_vector.vector);
        gsl_permutation_reverse(pix_perm);
    }
    return pix_perm;
}


/* Exponentiate and normalize a log probability sky map. */
static void exp_normalize(long npix, double *P, gsl_permutation *pix_perm)
{
    long i;
    double accum, max_log_p;

    /* Find the value of the greatest log probability. */
    max_log_p = P[gsl_permutation_get(pix_perm, 0)];

    /* Subtract it off. */
    for (i = 0; i < npix; i ++)
        P[i] -= max_log_p;

    /* Exponentiate to convert from log probability to probability. */
    for (i = 0; i < npix; i ++)
        P[i] = exp(P[i]);

    /* Sum entire sky map to find normalization. */
    for (accum = 0, i = 0; i < npix; i ++)
        accum += P[gsl_permutation_get(pix_perm, i)];

    /* Normalize. */
    for (i = 0; i < npix; i ++)
        P[i] /= accum;
}


static long indexof_confidence_level(
    long npix, double *P, double level, gsl_permutation *pix_perm)
{
    double accum;
    long maxpix;

    for (accum = 0, maxpix = 0; maxpix < npix && accum <= level; maxpix ++)
        accum += P[gsl_permutation_get(pix_perm, maxpix)];

    return maxpix;
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


static int bayestar_sky_map_toa_not_normalized_log(
    long npix,                      /* Number of HEALPix pixels */
    double *P,                      /* Preallocated posterior pixel map */
    /* Detector network */
    double gmst,                    /* GMST (rad) */
    unsigned int nifos,             /* Number of detectors */
    unsigned long nsamples,         /* Length of autocorrelation sequence */
    double sample_rate,             /* Sample rate in seconds */
    const double complex **acors,   /* Autocorrelation sequences */
    const double **locations,       /* Barycentered Cartesian geographic detector positions (m) */
    /* Observations */
    const double *toas,             /* Arrival time differences relative to network barycenter (s) */
    const double *snrs              /* SNRs */
) {
    long nside;

    /* Determine the lateral HEALPix resolution. */
    nside = npix2nside(npix);
    if (nside < 0)
        GSL_ERROR("output is not a valid HEALPix array", GSL_EINVAL);

    /* Loop over pixels. */
    #pragma omp parallel for
    for (long i = 0; i < npix; i ++)
    {
        double log_posterior = -INFINITY;
        double dt[nifos];

        /* Compute time of arrival errors */
        {
            double theta, phi;
            pix2ang_ring(nside, i, &theta, &phi);
            toa_errors(dt, theta, phi, gmst, nifos, locations, toas);
        }

        /* Loop over arrival time */
        for (long isample = 1 - (long)nsamples; isample < (long)nsamples; isample++)
        {
            double product = 1;
            double scale = 0;

            /* Loop over detectors */
            for (unsigned int iifo = 0; iifo < nifos; iifo++)
            {
                const double I0arg = gsl_pow_2(snrs[iifo]) * cabs(eval_acor(
                    acors[iifo], nsamples, dt[iifo] * sample_rate + isample));
                product *= gsl_sf_bessel_I0_scaled(I0arg);
                scale += I0arg;
            }

            log_posterior = logaddexp(log_posterior, scale + log(product));
        }

        P[i] = log_posterior;
    }

    /* Done! */
    return GSL_SUCCESS;
}


static const double autoresolution_confidence_level = 0.9999;
static const long autoresolution_count_pix = 12288;
static const unsigned int nglfixed = 10;
static const unsigned int ntwopsi = 10;


#if GSL_MINOR_VERSION >= 15 /* FIXME: add dependency on GSL >= 1.15 */
static double toa_phoa_snr_log_radial_integral(
    double r1, double r2, double p2, double p, double b, int k,
    const gsl_integration_glfixed_table *gltable)
{
    double result;

    if (p2 == 0) {
        /* note: p2 == 0 implies b == 0 */
        assert(b == 0);
        int k1 = k + 1;

        if (k == -1)
        {
            result = log(log(r2 / r1));
        } else {
            result = log((gsl_pow_int(r2, k1) - gsl_pow_int(r1, k1)) / k1);
        }
    } else if (b == 0) {
        int k1 = k + 1;
        int k3 = k + 3;
        double arg1 = p2 / gsl_pow_2(r1);
        double arg2 = p2 / gsl_pow_2(r2);
        double E1 = halfinteger_En_scaled(k3, arg1);
        double E2 = halfinteger_En_scaled(k3, arg2);
        result = log1p(
            -gsl_sf_exp_mult(-(arg1 - arg2), gsl_pow_int(r1/r2, k1) * E1/E2))
            + k1 * log(r2) - arg2 + log(E2) - M_LN2;
    } else {
        double one_by_ml_r = 0.5 * b / p2;
        double ml_r = 1 / one_by_ml_r;
        if (ml_r > r1 &&
            ml_r < r2)
        {
            result = 0;
            double Gmin = M_SQRTPI * gsl_cdf_ugaussian_Q(
                M_SQRT2 * (p / r1 - 0.5 * b / p)) / p;
            double Gmax = M_SQRTPI * gsl_cdf_ugaussian_Q(
                M_SQRT2 * (p / r2 - 0.5 * b / p)) / p;
            for (unsigned int iG = 0; iG < gltable->n; iG ++)
            {
                double G, Gweight;
                {
                    /* Look up Gauss-Legendre abscissa and weight. */
                    int ret = gsl_integration_glfixed_point(
                        Gmin, Gmax, iG, &G, &Gweight, gltable);

                    /* Don't bother checking return value; the only
                     * possible failure is in index bounds checking. */
                    assert(ret == GSL_SUCCESS);
                }
                double r = 1 / (one_by_ml_r + M_SQRT1_2 / p
                    * gsl_cdf_ugaussian_Qinv(p * G / M_SQRTPI));
                result += Gweight * gsl_sf_bessel_I0_scaled(b / r)
                    * gsl_pow_int(r * one_by_ml_r, k) * gsl_pow_2(r);
            }
            result = log(result) + 0.5 * b * one_by_ml_r + k * log(ml_r);
        } else {
            result = 0;
            double rstar = k < 0 ? r1 : r2;
            double one_by_rstar = 1 / rstar;
            double log_scale =
                (b - p2 * one_by_rstar) * one_by_rstar +
                log(gsl_sf_bessel_I0_scaled(b * one_by_rstar)) +
                k * log(rstar);

            for (unsigned int ir = 0; ir < gltable->n; ir ++)
            {
                double r, rweight;
                {
                    /* Look up Gauss-Legendre abscissa and weight. */
                    int ret = gsl_integration_glfixed_point(
                        r1, r2, ir, &r,
                        &rweight, gltable);

                    /* Don't bother checking return value; the only
                     * possible failure is in index bounds checking. */
                    assert(ret == GSL_SUCCESS);
                }
                double one_by_r = 1 / r;
                result += rweight
                    * exp((b - p2 * one_by_r) * one_by_r - log_scale)
                    * gsl_sf_bessel_I0_scaled(b * one_by_r) * gsl_pow_int(r, k);
            }
            result = log(result) + log_scale;
        }
    }

    return result;
}
#endif /* GSL_MINOR_VERSION >= 1.15 */


static double *bayestar_sky_map_toa_adapt_resolution(
    gsl_permutation **pix_perm,
    long *maxpix,
    long *npix,
    /* Detector network */
    double gmst,                    /* GMST (rad) */
    unsigned int nifos,             /* Number of detectors */
    unsigned long nsamples,         /* Length of autocorrelation sequence */
    double sample_rate,             /* Sample rate in seconds */
    const double complex **acors,   /* Autocorrelation sequences */
    const double **locations,       /* Barycentered Cartesian geographic detector positions (m) */
    /* Observations */
    const double *toas,             /* Arrival time differences relative to network barycenter (s) */
    const double *snrs              /* SNRs */
) {
    int ret;
    double *P = NULL;
    long my_npix = *npix;
    long my_maxpix = *npix;
    gsl_permutation *my_pix_perm = NULL;

    if (my_npix == -1)
    {
        my_npix = autoresolution_count_pix / 4;
        do {
            my_npix *= 4;

            free(P);
            gsl_permutation_free(my_pix_perm);

            P = malloc(my_npix * sizeof(double));
            if (!P)
                GSL_ERROR_NULL("failed to allocate output array", GSL_ENOMEM);
            ret = bayestar_sky_map_toa_not_normalized_log(my_npix, P, gmst,
                nifos, nsamples, sample_rate, acors, locations, toas, snrs);
            if (ret != GSL_SUCCESS)
            {
                free(P);
                P = NULL;
                goto fail;
            }
            my_pix_perm = get_pixel_ranks(my_npix, P);
            if (!my_pix_perm)
            {
                free(P);
                P = NULL;
                goto fail;
            }
            exp_normalize(my_npix, P, my_pix_perm);

            my_maxpix = indexof_confidence_level(my_npix, P,
                autoresolution_confidence_level, my_pix_perm);
        } while (my_maxpix < autoresolution_count_pix);
    } else {
        P = malloc(my_npix * sizeof(double));
        if (!P)
            GSL_ERROR_NULL("failed to allocate output array", GSL_ENOMEM);
        ret = bayestar_sky_map_toa_not_normalized_log(my_npix, P, gmst,
            nifos, nsamples, sample_rate, acors, locations, toas, snrs);
        if (ret != GSL_SUCCESS)
        {
            free(P);
            P = NULL;
            goto fail;
        }
        my_pix_perm = get_pixel_ranks(my_npix, P);
        if (!my_pix_perm)
        {
            free(P);
            P = NULL;
            goto fail;
        }
        exp_normalize(my_npix, P, my_pix_perm);

        my_maxpix = indexof_confidence_level(my_npix, P,
            autoresolution_confidence_level, my_pix_perm);
    }

    *npix = my_npix;
    *pix_perm = my_pix_perm;
    *maxpix = my_maxpix;
fail:
    return P;
}


double *bayestar_sky_map_toa(
    long *npix,
    /* Detector network */
    double gmst,                    /* GMST (rad) */
    unsigned int nifos,             /* Number of detectors */
    unsigned long nsamples,         /* Length of autocorrelation sequence */
    double sample_rate,             /* Sample rate in seconds */
    const double complex **acors,   /* Autocorrelation sequences */
    const double **locations,       /* Barycentered Cartesian geographic detector positions (m) */
    /* Observations */
    const double *toas,             /* Arrival time differences relative to network barycenter (s) */
    const double *snrs              /* SNRs */
) {
    long maxpix;
    gsl_permutation *pix_perm = NULL;
    double *ret = bayestar_sky_map_toa_adapt_resolution(&pix_perm, &maxpix,
        npix, gmst, nifos, nsamples, sample_rate, acors, locations, toas, snrs);
    gsl_permutation_free(pix_perm);
    return ret;
}


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
    for (ssize_t i = (ssize_t)map->len - 1; i >= 0; i --)
    {
        map->pixels[i].log4p -= max_log4p;
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
        {
            long ipix_ring;
            nest2ring(nside, ipix_nest, &ipix_ring);
            P[ipix_ring] = value;
        }
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
#if GSL_MINOR_VERSION >= 15 /* FIXME: add dependency on GSL >= 1.15 */
    double complex exp_i_phoas[nifos];
    for (unsigned int iifo = 0; iifo < nifos; iifo ++)
        exp_i_phoas[iifo] = exp_i(phoas[iifo]);

    static const unsigned char order0 = 4;
    unsigned char level = order0;
    adaptive_sky_map *map = adaptive_sky_map_alloc(order0);
    if (!map)
        return NULL;
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

                        double result = toa_phoa_snr_log_radial_integral(
                            min_distance, max_distance, p2, p, b,
                            prior_distance_power, gltable);
                        accum2 = logaddexp(accum2, result);
                    }

                    accum1 = logaddexp(accum1, accum2 + log(weight));
                }

                accum = logaddexp(accum, accum1);
            }

            /* Record logarithm base 4 of posterior. */
            pixel->log4p = accum - 2 * M_LN2;
        }

        /* Restore old error handler. */
        gsl_set_error_handler(old_handler);

        /* Sort pixels by descending log posterior. */
        adaptive_sky_map_sort(map);

        /* If we have reached order=11 (nside=2048), stop. */
        if (level++ >= 11)
            break;

        /* Adaptively refine the pixels that contain the most probability. */
        map = adaptive_sky_map_refine(map, npix0 / 4);
        if (!map)
            return NULL;
    }

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
#else /* GSL_MINOR_VERSION < 15 */
    /* Unused arguments */
    (void)inout_npix,
    (void)min_distance;
    (void)max_distance;
    (void)prior_distance_power;
    (void)gmst;
    (void)nifos;
    (void)nsamples;
    (void)sample_rate;
    (void)acors;
    (void)responses;
    (void)locations;
    (void)horizons;
    (void)toas;
    (void)phoas;
    (void)snrs;
    /* Unused functions */
    (void)my_gsl_error;
    (void)ignore_underflow;
    (void)adaptive_sky_map_sort;
    (void)adaptive_sky_map_refine;
    (void)adaptive_sky_map_alloc;
    (void)adaptive_sky_map_rasterize;
    GSL_ERROR_NULL("requires GSL >= 1.15", GSL_EFAILED);
#endif /* GSL_MINOR_VERSION >= 15 */
}


double bayestar_log_likelihood_toa(
    /* Parameters */
    double ra,                      /* Right ascension (rad) */
    double sin_dec,                 /* Sin(declination) */
    double t,                       /* Barycentered arrival time (s) */
    /* Detector network */
    double gmst,                    /* GMST (rad) */
    unsigned int nifos,             /* Number of detectors */
    unsigned long nsamples,         /* Length of autocorrelation sequence */
    double sample_rate,             /* Sample rate in seconds */
    const double complex **acors,   /* Autocorrelation sequences */
    const double **locations,       /* Barycentered Cartesian geographic detector positions (m) */
    /* Observations */
    const double *toas,             /* Arrival time differences relative to network barycenter (s) */
    const double *snrs              /* SNRs */
) {
    const double dec = asin(sin_dec);

    /* Compute time of arrival errors */
    double dt[nifos];
    toa_errors(dt, M_PI_2 - dec, ra, gmst, nifos, locations, toas);
    for (unsigned int iifo = 0; iifo < nifos; iifo++)
        dt[iifo] -= t;

    double product = 1;
    double scale = 0;

    /* Loop over detectors */
    for (unsigned int iifo = 0; iifo < nifos; iifo++)
    {
        const double I0arg = gsl_pow_2(snrs[iifo]) * cabs(eval_acor(
            acors[iifo], nsamples, dt[iifo] * sample_rate));
        product *= gsl_sf_bessel_I0_scaled(I0arg);
        scale += I0arg;
    }

    return log(product) + scale;
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


static void test_divfloor(int x, int y)
{
    int result = divfloor(x, y);
    int expected = floor((double)x / y);
    gsl_test_int(result, expected, "testing divfloor(%d, %d)", x, y);
}


static void test_pow_halfint(double x, int twon)
{
    double result = pow_halfint(x, twon);
    double expected = pow(x, 0.5 * twon);
    gsl_test_abs(result, expected, 1e-8, "testing pow_halfint(%g, %d)", x, twon);
}


static void test_halfinteger_En_scaled(int twon, double x)
{
    double expected, result = halfinteger_En_scaled(twon, x);

    if (GSL_IS_EVEN(twon) && twon >= 0)
    {
        expected = gsl_sf_expint_En_scaled(twon / 2, x);
    } else {
        expected = gsl_sf_exp_mult(x, pow_halfint(x, twon - 2)
            * gsl_sf_gamma_inc(1 - 0.5 * twon, x));
    }

    gsl_test_abs(result, expected, 1e-8, "testing E_{%d/2}(%g)", twon, x);
}


static void test_cabs2(double complex z)
{
    double result = cabs2(z);
    double expected = gsl_pow_2(cabs(z));
    gsl_test_abs(result, expected, 2 * GSL_DBL_EPSILON,
        "testing cabs2(%g + %g j)", creal(z), cimag(z));
}


static void test_catrom(void)
{
    for (double t = 0; t <= 1; t += 0.01)
    {
        double complex result = catrom(0, 0, 0, 0, t);
        double complex expected = 0;
        gsl_test_abs(creal(result), creal(expected), 0,
            "testing Catmull-rom interpolant for zero input");
        gsl_test_abs(cimag(result), cimag(expected), 0,
            "testing Catmull-rom interpolant for zero input");
    }

    for (double t = 0; t <= 1; t += 0.01)
    {
        double complex result = catrom(1, 1, 1, 1, t);
        double complex expected = 1;
        gsl_test_abs(creal(result), creal(expected), 0,
            "testing Catmull-rom interpolant for unit input");
        gsl_test_abs(cimag(result), cimag(expected), 0,
            "testing Catmull-rom interpolant for unit input");
    }

    for (double t = 0; t <= 1; t += 0.01)
    {
        double complex result = catrom(1.0j, 1.0j, 1.0j, 1.0j, t);
        double complex expected = 1.0j;
        gsl_test_abs(creal(result), creal(expected), 0,
            "testing Catmull-rom interpolant for unit imaginary input");
        gsl_test_abs(cimag(result), cimag(expected), 1,
            "testing Catmull-rom interpolant for unit imaginary input");
    }

    for (double t = 0; t <= 1; t += 0.01)
    {
        double complex result = catrom(1.0+1.0j, 1.0+1.0j, 1.0+1.0j, 1.0+1.0j, t);
        double complex expected = 1.0+1.0j;
        gsl_test_abs(creal(result), creal(expected), 0,
            "testing Catmull-rom interpolant for unit real + imaginary input");
        gsl_test_abs(cimag(result), cimag(expected), 0,
            "testing Catmull-rom interpolant for unit real + imaginary input");
    }

    for (double t = 0; t <= 1; t += 0.01)
    {
        double complex result = catrom(1, 0, 1, 4, t);
        double complex expected = gsl_pow_2(t);
        gsl_test_abs(creal(result), creal(expected), 0,
            "testing Catmull-rom interpolant for quadratic real input");
        gsl_test_abs(cimag(result), cimag(expected), 0,
            "testing Catmull-rom interpolant for quadratic real input");
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


static void test_toa_phoa_snr_log_radial_integral(
    double expected, double tol, double r1, double r2, double p2, double b, int k)
{
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

    double result = toa_phoa_snr_log_radial_integral(
        r1, r2, p2, sqrt(p2), b, k, gltable);

    gsl_test_rel(
        result, expected, tol, "testing toa_phoa_snr_log_radial_integral("
        "r1=%g, r2=%g, p2=%g, b=%g, k=%d)", r1, r2, p2, b, k);
}


int bayestar_test(void)
{
    for (int x = -512; x < 512; x ++)
        for (int y = -512; y < 512; y ++)
            if (y != 0)
                test_divfloor(x, y);

    for (int twon = -5; twon <= 5; twon ++)
        for (double log10x = -2; log10x <= 2; log10x += 0.1)
            test_pow_halfint(pow(10, log10x), twon);

    for (int twon = -5; twon <= 5; twon ++)
        for (double log10x = -2; log10x <= 2; log10x += 0.1)
            test_halfinteger_En_scaled(twon, pow(10, log10x));

    for (double re = -1; re < 1; re += 0.1)
        for (double im = -1; im < 1; im += 0.1)
            test_cabs2(re + im * 1.0j);

    test_catrom();
    test_eval_acor();

    for (double ra = -M_PI; ra <= M_PI; ra += 0.1 * M_PI)
        for (double dec = -M_PI_2; dec <= M_PI_2; dec += 0.1 * M_PI)
            for (double inc = -M_PI; inc <= M_PI; inc += 0.1 * M_PI)
                for (double pol = -M_PI; pol <= M_PI; pol += 0.1 * M_PI)
                    for (unsigned long epoch = 1000000000000000000; epoch <= 1000086400000000000; epoch += 3600000000000)
                        test_signal_amplitude_model(ra, dec, inc, pol, epoch, "H1");

    /* Tests of radial integrand with p2=0, b=0. */
    test_toa_phoa_snr_log_radial_integral(0, 0, 0, 1, 0, 0, 0);
    test_toa_phoa_snr_log_radial_integral(0, 0, exp(1), exp(2), 0, 0, -1);
    test_toa_phoa_snr_log_radial_integral(log(63), 0, 3, 6, 0, 0, 2);
    /* Tests of radial integrand with p2>0, b=0 (from Mathematica). */
    test_toa_phoa_snr_log_radial_integral(-0.480238, 1e-5, 1, 2, 1, 0, 0);
    test_toa_phoa_snr_log_radial_integral(0.432919, 1e-5, 1, 2, 1, 0, 2);
    test_toa_phoa_snr_log_radial_integral(-2.76076, 1e-5, 0, 1, 1, 0, 2);
    test_toa_phoa_snr_log_radial_integral(61.07118, 1e-5, 0, 1e9, 1, 0, 2);
    test_toa_phoa_snr_log_radial_integral(-112.23053, 1e-5, 0, 0.1, 1, 0, 2);
    test_toa_phoa_snr_log_radial_integral(-1.00004e6, 1e-5, 0, 1e-3, 1, 0, 2);
    /* Tests of radial integrand with p2>0, b>0 with ML peak outside
     * of integration limits (true values from Mathematica NIntegrate). */
    test_toa_phoa_snr_log_radial_integral(2.94548, 1e-4, 0, 4, 1, 1, 2);
    test_toa_phoa_snr_log_radial_integral(2.94545, 1e-4, 0.5, 4, 1, 1, 2);
    test_toa_phoa_snr_log_radial_integral(2.94085, 1e-4, 1, 4, 1, 1, 2);
    /* Tests of radial integrand with p2>0, b>0 with ML peak outside
     * of integration limits (true values from Mathematica NIntegrate). */
    test_toa_phoa_snr_log_radial_integral(-2.43264, 1e-5, 0, 1, 1, 1, 2);
    test_toa_phoa_snr_log_radial_integral(-2.43808, 1e-5, 0.5, 1, 1, 1, 2);
    test_toa_phoa_snr_log_radial_integral(-0.707038, 1e-5, 1, 1.5, 1, 1, 2);

    return gsl_test_summary();
}
