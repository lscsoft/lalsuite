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

#include <float.h>
#include <math.h>
#include <string.h>

#include <lal/DetResponse.h>

#include <chealpix.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sort_vector_double.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_vector.h>

#include "logaddexp.h"


/* Custom, reentrant GSL error handler that simply prints the error message. */
static void
my_gsl_error (const char *reason, const char *file, int line, int gsl_errno)
{
    (void)gsl_errno;
    fprintf(stderr, "gsl: %s:%d: %s: %s\n", file, line, "ERROR", reason);
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


static double complex eval_acor(
    const double complex *x,
    size_t nsamples,
    double t
) {
    /* Break |t| into integer and fractional parts. */
    size_t i;
    double f;
    double complex y;

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


/* Expression for complex amplitude on arrival (without 1/distance factor) */
static double complex signal_amplitude_model(
    double complex F,               /* Antenna factor */
    double complex exp_i_twopsi,    /* e^(i*2*psi), for polarization angle psi */
    double u,                       /* cos(inclination) */
    double u2                       /* cos^2(inclination */
) {
    const double complex tmp = F * exp_i_twopsi;
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
static const long autoresolution_count_pix_toa_snr = 3072;
static const long autoresolution_count_pix_toa_phoa_snr = 12288;
static const unsigned int nu = 16;
static const unsigned int ntwopsi = 16;


static double *bayestar_sky_map_toa_adapt_resolution(
    gsl_permutation **pix_perm,
    long *maxpix,
    long *npix,
    long autoresolution_count_pix,
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
        npix, autoresolution_count_pix_toa_snr, gmst, nifos, nsamples,
        sample_rate, acors, locations, toas, snrs);
    gsl_permutation_free(pix_perm);
    return ret;
}


static double complex exp_i(double phi) {
    return cos(phi) + I * sin(phi);
}


static double cabs2(double complex z) {
    return gsl_pow_2(creal(z)) + gsl_pow_2(cimag(z));
}

typedef struct {
    int prior_distance_power;
    double A;
    double B;
    double C;
} base_radial_integrand_params;


typedef struct {
    base_radial_integrand_params base;
    unsigned int n;
    double *I0args_times_r;
} toa_snr_radial_integrand_params;


typedef struct {
    base_radial_integrand_params base;
    double I0arg_times_r;
} toa_phoa_snr_radial_integrand_params;


static double base_log_radial_integrand(
    double r,
    const base_radial_integrand_params *p
) {
    const double one_by_r = 1 / r;

    return (p->A * one_by_r + p->B) * one_by_r + p->C
        + p->prior_distance_power * log(r);
}


static double toa_snr_log_radial_integrand(double r, void *params)
{
    const toa_snr_radial_integrand_params *p = params;
    double product = 1;
    double one_by_r = 1 / r;

    for (unsigned int i = 0; i < p->n; i ++)
        product *= gsl_sf_bessel_I0_scaled(p->I0args_times_r[i] * one_by_r);

    return base_log_radial_integrand(r, &p->base) + log(product);
}


static double toa_phoa_snr_log_radial_integrand(double r, void *params)
{
    const toa_phoa_snr_radial_integrand_params *p = params;
    return base_log_radial_integrand(r, &p->base)
        + log(gsl_sf_bessel_I0_scaled(p->I0arg_times_r / r));
}


static double base_radial_integrand(
    double r,
    const base_radial_integrand_params *p
) {
    const double one_by_r = 1 / r;
    return gsl_sf_exp_mult((p->A * one_by_r + p->B) * one_by_r + p->C,
        gsl_pow_int(r, p->prior_distance_power));
}


static double toa_snr_radial_integrand(double r, void *params)
{
    const toa_snr_radial_integrand_params *p = params;
    double product = 1;
    double one_by_r = 1 / r;

    for (unsigned int i = 0; i < p->n; i ++)
        product *= gsl_sf_bessel_I0_scaled(p->I0args_times_r[i] * one_by_r);

    return base_radial_integrand(r, &p->base) * product;
}


static double toa_phoa_snr_radial_integrand(double r, void *params)
{
    const toa_phoa_snr_radial_integrand_params *p = params;
    return base_radial_integrand(r, &p->base)
        * gsl_sf_bessel_I0_scaled(p->I0arg_times_r / r);
}


static int log_radial_integral(
    double *result,
    base_radial_integrand_params *p,
    double (*radial_integrand) (double, void *),
    double (*log_radial_integrand) (double, void *),
    double min_distance,
    double max_distance
) {
    double breakpoints[5];
    unsigned char nbreakpoints = 0;
    double integral, log_offset = -INFINITY;
    int ret;

    {
        /* Calculate the approximate distance at which the integrand attains a
         * maximum (middle) and a fraction eta of the maximum (left and right).
         * This neglects the scaled Bessel function factors and the power-law
         * distance prior. It assumes that the likelihood is approximately of
         * the form
         *
         *    A/r^2 + B/r.
         *
         * Then the middle breakpoint occurs at 1/r = -B/2A, and the left and
         * right breakpoints occur when
         *
         *   A/r^2 + B/r = log(eta) - B^2/4A.
         */

        static const double eta = 0.01;
        const double middle = -2 * p->A / p->B;
        const double left = 1 / (1 / middle + sqrt(log(eta) / p->A));
        const double right = 1 / (1 / middle - sqrt(log(eta) / p->A));

        /* Use whichever of the middle, left, and right points lie within the
         * integration limits as initial subdivisions for the adaptive
         * integrator. */

        breakpoints[nbreakpoints++] = min_distance;
        if(left > breakpoints[nbreakpoints-1] && left < max_distance)
            breakpoints[nbreakpoints++] = left;
        if(middle > breakpoints[nbreakpoints-1] && middle < max_distance)
            breakpoints[nbreakpoints++] = middle;
        if(right > breakpoints[nbreakpoints-1] && right < max_distance)
            breakpoints[nbreakpoints++] = right;
        breakpoints[nbreakpoints++] = max_distance;
    }

    /* Re-scale the integrand so that the maximum value at any of the
     * breakpoints is 1. Note that the initial value of the constant term
     * is overwritten. */

    for (unsigned char i = 0; i < nbreakpoints; i++)
    {
        double new_log_offset = log_radial_integrand(breakpoints[i], p);
        if (new_log_offset > log_offset)
            log_offset = new_log_offset;
    }

    /* If the largets value of the log integrand was -INFINITY, then the
     * integrand is 0 everywhere. Set log_offset to 0, because subtracting
     * -INFINITY would make the integrand infinite. */
    if (log_offset == -INFINITY)
        log_offset = 0;

    p->C -= log_offset;

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
        const gsl_function func = {radial_integrand, p};

        /* FIXME: do we care to keep the error estimate around? */
        double abserr;

        /* Perform adaptive Gaussian quadrature. */
        ret = gsl_integration_qagp(&func, breakpoints, nbreakpoints,
            DBL_MIN, 0.05, n, &workspace, &integral, &abserr);
    }

    /* If integrator succceeded, then restore scale factor and write output. */
    if (ret == GSL_SUCCESS)
        *result = log(integral) + log_offset;

    /* Done! */
    return ret;
}


double *bayestar_sky_map_toa_snr(
    long *npix,
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
    const double *snrs              /* SNRs */
) {
    long nside;
    long maxpix;
    double *P;
    gsl_permutation *pix_perm;

    /* Hold GSL return values for any thread that fails. */
    int gsl_errno = GSL_SUCCESS;

    /* Storage for old GSL error handler. */
    gsl_error_handler_t *old_handler;

    /* Evaluate posterior term only first. */
    P = bayestar_sky_map_toa_adapt_resolution(&pix_perm, &maxpix, npix,
        autoresolution_count_pix_toa_snr, gmst, nifos, nsamples, sample_rate,
        acors, locations, toas, snrs);
    if (!P)
        return NULL;

    /* Determine the lateral HEALPix resolution. */
    nside = npix2nside(*npix);

    /* Zero pixels that didn't meet the TDOA cut. */
    for (long i = maxpix; i < *npix; i ++)
    {
        long ipix = gsl_permutation_get(pix_perm, i);
        P[ipix] = -INFINITY;
    }

    /* Use our own error handler while in parallel section to avoid concurrent
     * calls to the GSL error handler, which if provided by the user may not
     * be threadsafe. */
    old_handler = gsl_set_error_handler(my_gsl_error);

    /* Compute posterior factor for amplitude consistency. */
    #pragma omp parallel for firstprivate(gsl_errno) lastprivate(gsl_errno)
    for (long i = 0; i < maxpix; i ++)
    {
        /* Cancel further computation if a GSL error condition has occurred.
         *
         * Note: if one thread sets gsl_errno, not necessarily all thread will
         * get the updated value. That's OK, because most failure modes will
         * cause GSL error conditions on all threads. If we cared to have any
         * failure on any thread terminate all of the other threads as quickly
         * as possible, then we would want to insert the following pragma here:
         *
         *     #pragma omp flush(gsl_errno)
         *
         * and likewise before any point where we set gsl_errno.
         */

        if (gsl_errno != GSL_SUCCESS)
            goto skip;

        {
            long ipix = gsl_permutation_get(pix_perm, i);
            double complex F[nifos];
            double dt[nifos];
            double accum = -INFINITY;

            {
                double theta, phi;
                pix2ang_ring(nside, ipix, &theta, &phi);

                /* Look up antenna factors */
                for (unsigned int iifo = 0; iifo < nifos; iifo++)
                {
                    XLALComputeDetAMResponse(
                        (double *)&F[iifo],     /* Type-punned real part */
                        1 + (double *)&F[iifo], /* Type-punned imag part */
                        responses[iifo], phi, M_PI_2 - theta, 0, gmst);
                    F[iifo] *= horizons[iifo];
                }

                toa_errors(dt, theta, phi, gmst, nifos, locations, toas);
            }

            /* Integrate over 2*psi */
            for (unsigned int itwopsi = 0; itwopsi < ntwopsi; itwopsi++)
            {
                const double twopsi = (2 * M_PI / ntwopsi) * itwopsi;
                const double complex exp_i_twopsi = exp_i(twopsi);

                /* Integrate over u; since integrand only depends on u^2 we only
                 * have to go from u=0 to u=1. We want to include u=1, so the upper
                 * limit has to be <= */
                for (unsigned int iu = 0; iu <= nu; iu++)
                {
                    const double u = (double)iu / nu;
                    const double u2 = gsl_pow_2(u);
                    double rho2_times_r2[nifos];
                    double A = 0;

                    for (unsigned int iifo = 0; iifo < nifos; iifo ++)
                    {
                        A += rho2_times_r2[iifo] = cabs2(signal_amplitude_model(
                            F[iifo], exp_i_twopsi, u, u2));
                    }
                    A *= -0.5;

                    for (long isample = 1 - (long)nsamples;
                        isample < (long)nsamples; isample++)
                    {
                        double B = 0;
                        double I0args_times_r[nifos];

                        for (unsigned int iifo = 0; iifo < nifos; iifo ++)
                        {
                            const double acor2 = cabs2(eval_acor(acors[iifo],
                                nsamples, dt[iifo] * sample_rate + isample));
                            const double I0arg_times_r = snrs[iifo]
                                * sqrt(acor2 * rho2_times_r2[iifo]);
                            B += I0args_times_r[iifo] = I0arg_times_r;
                        }

                        double result;
                        toa_snr_radial_integrand_params params = {
                            .base = {
                                .A = A,
                                .B = B,
                                .C = 0,
                                .prior_distance_power = prior_distance_power
                            },
                            .I0args_times_r = I0args_times_r,
                            .n = nifos
                        };
                        int ret = log_radial_integral(
                            &result, &params.base,
                            toa_snr_radial_integrand,
                            toa_snr_log_radial_integrand,
                            min_distance, max_distance);

                        /* If the integrator failed, then record the GSL error
                         * value for later reporting when we leave the parallel
                         * section. Then, break out of the loop. */
                        if (ret != GSL_SUCCESS)
                        {
                            gsl_errno = ret;
                            goto skip;
                        }

                        accum = logaddexp(accum, result);
                    }
                }
            }

            /* Accumulate (log) posterior terms for SNR and TDOA. */
            P[ipix] += accum;
        }

        skip: /* this statement intentionally left blank */;
    }

    /* Restore old error handler. */
    gsl_set_error_handler(old_handler);

    /* Free permutation. */
    gsl_permutation_free(pix_perm);

    /* Check if there was an error in any thread evaluating any pixel. If there
     * was, raise the error and return. */
    if (gsl_errno != GSL_SUCCESS)
    {
        free(P);
        GSL_ERROR_NULL(gsl_strerror(gsl_errno), gsl_errno);
    }

    /* Exponentiate and normalize posterior. */
    pix_perm = get_pixel_ranks(*npix, P);
    if (!pix_perm)
    {
        free(P);
        return NULL;
    }
    exp_normalize(*npix, P, pix_perm);
    gsl_permutation_free(pix_perm);

    return P;
}


double *bayestar_sky_map_toa_phoa_snr(
    long *npix,
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
    long nside;
    long maxpix;
    double *P;
    gsl_permutation *pix_perm;

    /* Hold GSL return values for any thread that fails. */
    int gsl_errno = GSL_SUCCESS;

    /* Storage for old GSL error handler. */
    gsl_error_handler_t *old_handler;

    double complex exp_i_phoas[nifos];
    for (unsigned int iifo = 0; iifo < nifos; iifo ++)
        exp_i_phoas[iifo] = exp_i(phoas[iifo]);

    /* Evaluate posterior term only first. */
    P = bayestar_sky_map_toa_adapt_resolution(&pix_perm, &maxpix, npix,
        autoresolution_count_pix_toa_phoa_snr, gmst, nifos, nsamples,
        sample_rate, acors, locations, toas, snrs);
    if (!P)
        return NULL;

    /* Determine the lateral HEALPix resolution. */
    nside = npix2nside(*npix);

    /* Zero pixels that didn't meet the TDOA cut. */
    for (long i = maxpix; i < *npix; i ++)
    {
        long ipix = gsl_permutation_get(pix_perm, i);
        P[ipix] = -INFINITY;
    }

    /* Use our own error handler while in parallel section to avoid concurrent
     * calls to the GSL error handler, which if provided by the user may not
     * be threadsafe. */
    old_handler = gsl_set_error_handler(my_gsl_error);

    /* Compute posterior factor for amplitude consistency. */
    #pragma omp parallel for firstprivate(gsl_errno) lastprivate(gsl_errno)
    for (long i = 0; i < maxpix; i ++)
    {
        /* Cancel further computation if a GSL error condition has occurred.
         *
         * Note: if one thread sets gsl_errno, not necessarily all thread will
         * get the updated value. That's OK, because most failure modes will
         * cause GSL error conditions on all threads. If we cared to have any
         * failure on any thread terminate all of the other threads as quickly
         * as possible, then we would want to insert the following pragma here:
         *
         *     #pragma omp flush(gsl_errno)
         *
         * and likewise before any point where we set gsl_errno.
         */

        if (gsl_errno != GSL_SUCCESS)
            goto skip;

        {
            long ipix = gsl_permutation_get(pix_perm, i);
            double complex F[nifos];
            double dt[nifos];
            double accum = -INFINITY;

            {
                double theta, phi;
                pix2ang_ring(nside, ipix, &theta, &phi);

                /* Look up antenna factors */
                for (unsigned int iifo = 0; iifo < nifos; iifo++)
                {
                    XLALComputeDetAMResponse(
                        (double *)&F[iifo],     /* Type-punned real part */
                        1 + (double *)&F[iifo], /* Type-punned imag part */
                        responses[iifo], phi, M_PI_2 - theta, 0, gmst);
                    F[iifo] *= horizons[iifo];
                }

                toa_errors(dt, theta, phi, gmst, nifos, locations, toas);
            }

            /* Integrate over 2*psi */
            for (unsigned int itwopsi = 0; itwopsi < ntwopsi; itwopsi++)
            {
                const double twopsi = (2 * M_PI / ntwopsi) * itwopsi;
                const double complex exp_i_twopsi = exp_i(twopsi);

                /* Integrate over u from u=-1 to u=1. */
                for (int iu = -(int)nu; iu <= (int)nu; iu++)
                {
                    const double u = (double)iu / nu;
                    const double u2 = gsl_pow_2(u);
                    double complex z_times_r[nifos];
                    double A = 0;

                    for (unsigned int iifo = 0; iifo < nifos; iifo ++)
                    {
                        A += cabs2(
                            z_times_r[iifo] = signal_amplitude_model(
                                F[iifo], exp_i_twopsi, u, u2));
                    }
                    A *= -0.5;

                    for (long isample = 1 - (long)nsamples;
                        isample < (long)nsamples; isample++)
                    {
                        double complex I0arg_complex_times_r = 0;

                        for (unsigned int iifo = 0; iifo < nifos; iifo ++)
                        {
                            I0arg_complex_times_r += snrs[iifo]
                                * z_times_r[iifo] * exp_i_phoas[iifo]
                                * eval_acor(acors[iifo], nsamples,
                                    dt[iifo] * sample_rate + isample);
                        }

                        double result;
                        double I0arg_times_r = cabs(I0arg_complex_times_r);
                        toa_phoa_snr_radial_integrand_params params = {
                            .base = {
                                .A = A,
                                .B = I0arg_times_r,
                                .C = 0,
                                .prior_distance_power = prior_distance_power
                            },
                            .I0arg_times_r = I0arg_times_r,
                        };
                        int ret = log_radial_integral(
                            &result, &params.base,
                            toa_phoa_snr_radial_integrand,
                            toa_phoa_snr_log_radial_integrand,
                            min_distance, max_distance);

                        /* If the integrator failed, then record the GSL error
                         * value for later reporting when we leave the parallel
                         * section. Then, break out of the loop. */
                        if (ret != GSL_SUCCESS)
                        {
                            gsl_errno = ret;
                            goto skip;
                        }

                        accum = logaddexp(accum, result);
                    }
                }
            }

            /* Accumulate (log) posterior terms for SNR and TDOA. */
            P[ipix] += accum;
        }

        skip: /* this statement intentionally left blank */;
    }

    /* Restore old error handler. */
    gsl_set_error_handler(old_handler);

    /* Free permutation. */
    gsl_permutation_free(pix_perm);

    /* Check if there was an error in any thread evaluating any pixel. If there
     * was, raise the error and return. */
    if (gsl_errno != GSL_SUCCESS)
    {
        free(P);
        GSL_ERROR_NULL(gsl_strerror(gsl_errno), gsl_errno);
    }

    /* Exponentiate and normalize posterior. */
    pix_perm = get_pixel_ranks(*npix, P);
    if (!pix_perm)
    {
        free(P);
        return NULL;
    }
    exp_normalize(*npix, P, pix_perm);
    gsl_permutation_free(pix_perm);

    return P;
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
        double complex F;
        XLALComputeDetAMResponse(
            (double *)&F,     /* Type-punned real part */
            1 + (double *)&F, /* Type-punned imag part */
            responses[iifo], ra, dec, 0, gmst);
        F *= horizons[iifo];

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
        double complex F;
        XLALComputeDetAMResponse(
            (double *)&F,     /* Type-punned real part */
            1 + (double *)&F, /* Type-punned imag part */
            responses[iifo], ra, dec, 0, gmst);
        F *= horizons[iifo];

        const double complex z_times_r =
             signal_amplitude_model(F, exp_i_twopsi, u, u2);
        const double complex zhat = snrs[iifo] * exp_i(phoas[iifo]);

        i0arg_complex_times_r += zhat * z_times_r
            * eval_acor(acors[iifo], nsamples, dt[iifo] * sample_rate);
        A += cabs2(z_times_r);
    }
    A *= -0.5;

    const double i0arg_times_r = cabs(i0arg_complex_times_r);

    return (A * one_by_r + i0arg_times_r) * one_by_r
        + log(gsl_sf_bessel_I0_scaled(i0arg_times_r * one_by_r));
}
