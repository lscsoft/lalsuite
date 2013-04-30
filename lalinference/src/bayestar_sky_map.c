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
#include <lal/LALSimulation.h>

#include <chealpix.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sort_vector_double.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_vector.h>

#include "logaddexp.h"


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


static long indexof_confidence_level(long npix, double *P, double level, gsl_permutation *pix_perm)
{
    double accum;
    long maxpix;

    for (accum = 0, maxpix = 0; maxpix < npix && accum <= level; maxpix ++)
        accum += P[gsl_permutation_get(pix_perm, maxpix)];

    return maxpix;
}


/* Perform sky localization based on TDOAs alone. Returns log probability; not normalized. */
static int bayestar_sky_map_tdoa_not_normalized_log(
    long npix, /* Input: number of HEALPix pixels. */
    double *P, /* Output: pre-allocated array of length npix to store posterior map. */
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    const double **locs, /* Input: array of detector positions. */
    const double *toas, /* Input: array of times of arrival. */
    const double *w_toas /* Input: sum-of-squares weights, (1/TOA variance)^2. */
) {
    long nside;
    long i;

    /* Determine the lateral HEALPix resolution. */
    nside = npix2nside(npix);
    if (nside < 0)
        GSL_ERROR("output is not a valid HEALPix array", GSL_EINVAL);

    /* Loop over pixels. */
    for (i = 0; i < npix; i ++)
    {
        int j;

        /* Determine polar coordinates of this pixel. */
        double theta, phi;
        pix2ang_ring(nside, i, &theta, &phi);

        /* Convert from equatorial to geographic coordinates. */
        phi -= gmst;

        /* Convert to Cartesian coordinates. */
        double n[3];
        ang2vec(theta, phi, n);

        /* Loop over detectors. */
        double dt[nifos];
        for (j = 0; j < nifos; j ++)
            dt[j] = toas[j] + cblas_ddot(3, n, 1, locs[j], 1) / LAL_C_SI;

        /* Evaluate the (un-normalized) Gaussian log likelihood. */
        P[i] = -0.5 * gsl_stats_wtss(w_toas, 1, dt, 1, nifos);
    }

    /* Done! */
    return GSL_SUCCESS;
}


static const double autoresolution_confidence_level = 0.9999;
static const long autoresolution_count_pix = 3072;


/* Perform sky localization based on TDOAs alone. */
static double *bayestar_sky_map_tdoa_adapt_resolution(
    gsl_permutation **pix_perm,
    long *maxpix,
    long *npix, /* In/out: number of HEALPix pixels. */
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    const double **locs, /* Input: array of detector positions. */
    const double *toas, /* Input: array of times of arrival. */
    const double *w_toas /* Input: sum-of-squares weights, (1/TOA variance)^2. */
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
            ret = bayestar_sky_map_tdoa_not_normalized_log(my_npix, P, gmst, nifos, locs, toas, w_toas);
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

            my_maxpix = indexof_confidence_level(my_npix, P, autoresolution_confidence_level, my_pix_perm);
        } while (my_maxpix < autoresolution_count_pix);
    } else {
        P = malloc(my_npix * sizeof(double));
        if (!P)
            GSL_ERROR_NULL("failed to allocate output array", GSL_ENOMEM);
        ret = bayestar_sky_map_tdoa_not_normalized_log(my_npix, P, gmst, nifos, locs, toas, w_toas);
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
    }

    *npix = my_npix;
    *pix_perm = my_pix_perm;
    *maxpix = my_maxpix;
fail:
    return P;
}


/* Perform sky localization based on TDOAs alone. */
double *bayestar_sky_map_tdoa(
    long *npix, /* In/out: number of HEALPix pixels. */
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    const double **locs, /* Input: array of detector positions. */
    const double *toas, /* Input: array of times of arrival. */
    const double *w_toas /* Input: sum-of-squares weights, (1/TOA variance)^2. */
) {
    long maxpix;
    gsl_permutation *pix_perm = NULL;
    double *ret = bayestar_sky_map_tdoa_adapt_resolution(&pix_perm, &maxpix, npix, gmst, nifos, locs, toas, w_toas);
    gsl_permutation_free(pix_perm);
    return ret;
}


typedef struct {
    double A;
    double B;
    double log_offset;
} inner_integrand_params;


/* Radial integrand for uniform-in-log-distance prior. */
static double radial_integrand_uniform_in_log_distance(double r, void *params)
{
    const inner_integrand_params *integrand_params = (const inner_integrand_params *) params;

    const double onebyr = 1 / r;
    const double onebyr2 = gsl_pow_2(onebyr);
    return exp(integrand_params->A * onebyr2 + integrand_params->B * onebyr - integrand_params->log_offset) * onebyr;
}


/* Radial integrand for uniform-in-volume prior. */
static double radial_integrand_uniform_in_volume(double r, void *params)
{
    const inner_integrand_params *integrand_params = (const inner_integrand_params *) params;

    const double onebyr = 1 / r;
    const double onebyr2 = gsl_pow_2(onebyr);
    return exp(integrand_params->A * onebyr2 + integrand_params->B * onebyr - integrand_params->log_offset) * gsl_pow_2(r);
}


double *bayestar_sky_map_tdoa_snr(
    long *npix, /* Input: number of HEALPix pixels. */
    double gmst, /* Greenwich mean sidereal time in radians. */
    int nifos, /* Input: number of detectors. */
    /* FIXME: make const; change XLALComputeDetAMResponse prototype */
    /* const */ float **responses, /* Pointers to detector responses. */
    const double **locations, /* Pointers to locations of detectors in Cartesian geographic coordinates. */
    const double *toas, /* Input: array of times of arrival with arbitrary relative offset. (Make toas[0] == 0.) */
    const double *snrs, /* Input: array of SNRs. */
    const double *w_toas, /* Input: sum-of-squares weights, (1/TOA variance)^2. */
    const double *horizons, /* Distances at which a source would produce an SNR of 1 in each detector. */
    double min_distance,
    double max_distance,
    bayestar_prior_t prior)
{
    long nside;
    long maxpix;
    long i;
    double d1[nifos];
    double *P;
    gsl_permutation *pix_perm;

    /* Function pointer to hold radial integrand. */
    double (* radial_integrand) (double x, void *params);

    /* Will point to memory for storing GSL return values for each thread. */
    int *gsl_errnos;

    /* Storage for old GSL error handler. */
    gsl_error_handler_t *old_handler;

    /* Maximum number of subdivisions for adaptive integration. */
    static const size_t subdivision_limit = 64;

    /* Subdivide radial integral where likelihood is this fraction of the maximum,
     * will be used in solving the quadratic to find the breakpoints */
    static const double eta = 0.01;

    /* Use this many integration steps in 2*psi  */
    static const int ntwopsi = 16;

    /* Number of integration steps in cos(inclination) */
    static const int nu = 16;

    /* Choose radial integrand function based on selected prior. */
    switch (prior)
    {
        case BAYESTAR_PRIOR_UNIFORM_IN_LOG_DISTANCE:
            radial_integrand = radial_integrand_uniform_in_log_distance;
            break;
        case BAYESTAR_PRIOR_UNIFORM_IN_VOLUME:
            radial_integrand = radial_integrand_uniform_in_volume;
            break;
        default:
            GSL_ERROR_NULL("unrecognized choice of prior", GSL_EINVAL);
            break;
    }

    /* Rescale distances so that furthest horizon distance is 1. */
    {
        double d1max;
        memcpy(d1, horizons, sizeof(d1));
        for (d1max = d1[0], i = 1; i < nifos; i ++)
            if (d1[i] > d1max)
                d1max = d1[i];
        for (i = 0; i < nifos; i ++)
            d1[i] /= d1max;
        min_distance /= d1max;
        max_distance /= d1max;
    }

    /* Evaluate posterior term only first. */
    P = bayestar_sky_map_tdoa_adapt_resolution(&pix_perm, &maxpix, npix, gmst, nifos, locations, toas, w_toas);
    if (!P)
        return NULL;

    /* Determine the lateral HEALPix resolution. */
    nside = npix2nside(*npix);

    /* Zero pixels that didn't meet the TDOA cut. */
    for (i = 0; i < maxpix; i ++)
    {
        long ipix = gsl_permutation_get(pix_perm, i);
        P[ipix] = log(P[ipix]);
    }
    for (; i < *npix; i ++)
    {
        long ipix = gsl_permutation_get(pix_perm, i);
        P[ipix] = -INFINITY;
    }

    /* Allocate space to store per-pixel, per-thread error value. */
    gsl_errnos = calloc(maxpix, sizeof(int));
    if (!gsl_errnos)
    {
        free(P);
        gsl_permutation_free(pix_perm);
        GSL_ERROR_NULL("failed to allocate space for pixel error status", GSL_ENOMEM);
    }

    /* Turn off error handler while in parallel section to avoid concurrent
     * calls to the GSL error handler, which if provided by the user may not
     * be threadsafe. */
    old_handler = gsl_set_error_handler_off();

    /* Compute posterior factor for amplitude consistency. */
    #pragma omp parallel for
    for (i = 0; i < maxpix; i ++)
    {
        long ipix = gsl_permutation_get(pix_perm, i);
        double F[nifos][2];
        double theta, phi;
        int itwopsi, iu, iifo;
        double accum = -INFINITY;

        /* Prepare workspace for adaptive integrator. */
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(subdivision_limit);

        /* If the workspace could not be allocated, then record the GSL
         * error value for later reporting when we leave the parallel
         * section. Then, skip to the next loop iteration. */
        if (!workspace)
        {
            gsl_errnos[i] = GSL_ENOMEM;
            continue;
        }

        /* Look up polar coordinates of this pixel */
        pix2ang_ring(nside, ipix, &theta, &phi);

        /* Look up antenna factors */
        for (iifo = 0; iifo < nifos; iifo ++)
        {
            XLALComputeDetAMResponse(&F[iifo][0], &F[iifo][1], (float (*)[3])responses[iifo], phi, M_PI_2 - theta, 0, gmst);
            F[iifo][0] *= d1[iifo];
            F[iifo][1] *= d1[iifo];
        }

        /* Integrate over 2*psi */
        for (itwopsi = 0; itwopsi < ntwopsi; itwopsi++)
        {
            const double twopsi = (2 * M_PI / ntwopsi) * itwopsi;
            const double costwopsi = cos(twopsi);
            const double sintwopsi = sin(twopsi);

            /* Integrate over u; since integrand only depends on u^2 we only
             * have to go from u=0 to u=1. We want to include u=1, so the upper
             * limit has to be <= */
            for (iu = 0; iu <= nu; iu++)
            {
                const double u = (double)iu / nu;
                const double u2 = gsl_pow_2(u);
                const double u4 = gsl_pow_2(u2);

                double A = 0, B = 0;
                double breakpoints[5];
                int num_breakpoints = 0;

                /* The log-likelihood is quadratic in the estimated and true
                 * values of the SNR, and in 1/r. It is of the form A/r^2 + B/r,
                 * where A depends only on the true values of the SNR and is
                 * strictly negative and B depends on both the true values and
                 * the estimates and is strictly positive.
                 *
                 * The middle breakpoint is at the maximum of the log-likelihood,
                 * occurring at 1/r=-B/2A. The lower and upper breakpoints occur
                 * when the likelihood becomes eta times its maximum value. This
                 * occurs when
                 *
                 *   A/r^2 + B/r = log(eta) - B^2/4A.
                 *
                 */

                /* Loop over detectors */
                for (iifo = 0; iifo < nifos; iifo++)
                {
                    const double Fp = F[iifo][0]; /* `plus' antenna factor times r */
                    const double Fx = F[iifo][1]; /* `cross' antenna factor times r */
                    const double FpFp = gsl_pow_2(Fp);
                    const double FxFx = gsl_pow_2(Fx);
                    const double FpFx = Fp * Fx;
                    const double rhotimesr2 = 0.125 * ((FpFp + FxFx) * (1 + 6*u2 + u4) + gsl_pow_2(1 - u2) * ((FpFp - FxFx) * costwopsi + 2 * FpFx * sintwopsi));
                    const double rhotimesr = sqrt(rhotimesr2);

                    /* FIXME: due to roundoff, rhotimesr2 can be very small and
                     * negative rather than simply zero. If this happens, don't
                     accumulate the log-likelihood terms for this detector. */
                    if (rhotimesr2 > 0)
                    {
                        A += rhotimesr2;
                        B += rhotimesr * snrs[iifo];
                    }
                }
                A *= -0.5;

                {
                    const double middle_breakpoint = -2 * A / B;
                    const double lower_breakpoint = 1 / (1 / middle_breakpoint + sqrt(log(eta) / A));
                    const double upper_breakpoint = 1 / (1 / middle_breakpoint - sqrt(log(eta) / A));
                    breakpoints[num_breakpoints++] = min_distance;
                    if(lower_breakpoint > breakpoints[num_breakpoints-1] && lower_breakpoint < max_distance)
                        breakpoints[num_breakpoints++] = lower_breakpoint;
                    if(middle_breakpoint > breakpoints[num_breakpoints-1] && middle_breakpoint < max_distance)
                        breakpoints[num_breakpoints++] = middle_breakpoint;
                    if(upper_breakpoint > breakpoints[num_breakpoints-1] && upper_breakpoint < max_distance)
                        breakpoints[num_breakpoints++] = upper_breakpoint;
                    breakpoints[num_breakpoints++] = max_distance;
                }

                {
                    /* Perform adaptive integration. Stop when a relative
                     * accuracy of 0.05 has been reached. */
                    inner_integrand_params integrand_params = {A, B, -0.25 * gsl_pow_2(B) / A};
                    const gsl_function func = {radial_integrand, &integrand_params};
                    double result, abserr;
                    int ret = gsl_integration_qagp(&func, &breakpoints[0], num_breakpoints, DBL_MIN, 0.05, subdivision_limit, workspace, &result, &abserr);

                    /* If the integrator failed, then record the GSL error
                     * value for later reporting when we leave the parallel
                     * section. Then, break out of the inner loop. */
                    if (ret != GSL_SUCCESS)
                    {
                        gsl_errnos[i] = ret;
                        break;
                    }

                    /* Take the logarithm and put the log-normalization back in. */
                    result = log(result) + integrand_params.log_offset;

                    /* Accumulate result. */
                    accum = logaddexp(accum, result);
                }
            }
        }
        /* Discard workspace for adaptive integrator. */
        gsl_integration_workspace_free(workspace);

        /* Accumulate (log) posterior terms for SNR and TDOA. */
        P[ipix] += accum;
    }

    /* Restore old error handler. */
    gsl_set_error_handler(old_handler);

    /* Free permutation. */
    gsl_permutation_free(pix_perm);

    /* Check if there was an error in any thread evaluating any pixel. If there
     * was, raise the error and return. */
    for (i = 0; i < maxpix; i ++)
    {
        int gsl_errno = gsl_errnos[i];
        if (gsl_errno != GSL_SUCCESS)
        {
            free(gsl_errnos);
            free(P);
            GSL_ERROR_NULL(gsl_strerror(gsl_errno), gsl_errno);
        }
    }

    /* Discard array of GSL error values. */
    free(gsl_errnos);

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
