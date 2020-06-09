/*  Copyright (C) 2019 Vijay Varma
 *  Evaluates NRSur7dq4Remnant model for remnant BH mass, spin and recoil kick
 *  for generically precessing BBH.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

/**
 * \author Vijay Varma
 *
 * \file
 *
 * \brief NRSur7dq4Remnant model for remnant BH mass, spin and recoil kick for
 * generically precessing BBH.
 *
 * The binary data file is available at:
 * https://dcc.ligo.org/LIGO-T1900393/public.
 * Get the lalsuite-extra repo or put the data into a location in your
 * LAL_DATA_PATH.
 *
 * **Paper**: https://arxiv.org/abs/1905.09300,
 * https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.1.033015
 *
 * **Parameter ranges**:
 *
 *   q = [1, 6.01]
 *
 *   \f$\chi_{1}, \chi_{2}\f$ = [-1, 1]
 *
 * **Training parameter ranges**:
 *
 *   q = [1, 4.01]
 *
 *   \f$\chi_{1}, \chi_{2}\f$ = [-0.81, 0.81]
 *
 *   But extrapolates reasonably to the above mass ratios and spins. However,
 *   if a guarantee of accuracy is required, this model should be used within
 *   the training parameter range.
 */

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


#include <libgen.h>

#include <assert.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_complex_math.h>

#include <lal/SeqFactories.h>
#include <lal/FileIO.h>
#include <lal/H5FileIO.h>
#include <lal/XLALError.h>
#include <lal/LALSimIMR.h>

#include "LALSimNRSurRemnantUtils.h"

#include "LALSimNRSur7dq4Remnant.h"


#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
#endif

#ifdef LAL_PTHREAD_LOCK
static pthread_once_t NRSur7dq4Remnant_is_initialized = PTHREAD_ONCE_INIT;
#endif


/**
 * Global surrogate data. This data will be loaded at most once. Any
 * executable which calls NRSur7dq4Remnant_Init_LALDATA directly or by calling
 * the XLALNRSur7dq4Remnant function will have a memory leak according to
 * valgrind, because we never free this memory. This memory, however, gets freed
 * after the executable terminates.
 */
static PrecessingRemnantFitData __lalsim_NRSur7dq4Remnant_data;

//*************************************************************************/
//************************* function definitions **************************/
//*************************************************************************/

/**
 * Helper function to check if the NRSur7dq4Remnant model has been initialized.
 */
static bool NRSur7dq4Remnant_IsSetup(void) {
  if(__lalsim_NRSur7dq4Remnant_data.setup)
    return true;
  else
    return false;
}

/**
 * Surrogate initialization.
 *
 * This needs to be called once, before __lalsim_NRSur7dq4Remnant_data is used.
 * It finds the H5 data file with the NRSur7dq4Remnant data and loads the
 * surrogate.  Can be called multiple times, will just return true if surrogate
 * is already loaded.
 */
static void NRSur7dq4Remnant_Init_LALDATA(void) {

    if (NRSur7dq4Remnant_IsSetup()) return;

    LALH5File *file = NULL;
    NRSurRemnant_LoadH5File(&file, NRSur7dq4Remnant_DATAFILE);

    int ret = PrecessingNRSurRemnant_Init(&__lalsim_NRSur7dq4Remnant_data,
            file);

    XLALH5FileClose(file);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR_VOID(XLAL_FAILURE, "Failure loading data from %s\n",
                NRSur7dq4Remnant_DATAFILE);
    }
}

/**
 * Map from mass ratio and spins to surrogate fit parameters.
 *
 * The fit parameters are \f$[log_e(q), \chi_{1x}, \chi_{1y}, \hat{\chi},
 *                              \chi_{2x}, \chi_{2y} \chi_a]\f$.
 * \f$\hat{\chi}\f$ is defined in Eq.(8) of arxiv:1905.09300.
 * \f$\chi_a = (\chi_{1z} - \chi_{2z})/2 \f$.
 *
 * The spins must be specified in the coorbital frame at t=-100M from
 * the total amplitude peak, as described in Sec.IV.C of arxiv:1905.09300.
 */
static int NRSur7dq4Remnant_fitParams(
    gsl_vector* fit_params, /**< Output: mapped fit parameters. */
    const REAL8 q,          /**< Mass ratio m1 / m2 >= 1. */
    const REAL8 chiAx,      /**< Dimless x-spin of heavier BH. */
    const REAL8 chiAy,      /**< Dimless y-spin of heavier BH. */
    const REAL8 chiAz,      /**< Dimless z-spin of heavier BH. */
    const REAL8 chiBx,      /**< Dimless x-spin of lighter BH. */
    const REAL8 chiBy,      /**< Dimless y-spin of lighter BH. */
    const REAL8 chiBz,      /**< Dimless z-spin of lighter BH. */
    LALDict* LALparams      /**< Dict with extra parameters */
) {

    // By default we do not allow unlimited_extrapolation
    UINT4 unlim_extrap = 0;
    if (LALparams != NULL &&
            XLALDictContains(LALparams, "unlimited_extrapolation")) {
        // Unless the user asks for it
        unlim_extrap
            = XLALDictLookupUINT4Value(LALparams, "unlimited_extrapolation");
    }

    //// Sanity check parameter ranges
    REAL8 chiAmag = sqrt(chiAx*chiAx + chiAy*chiAy + chiAz*chiAz);
    REAL8 chiBmag = sqrt(chiBx*chiBx + chiBy*chiBy + chiBz*chiBz);
    // These are the limits of the training parameter space, so print a
    // warning when extrapolated beyond.
    REAL8 q_max_soft = 4.01;
    REAL8 chi_max_soft = 0.81;
    // These are the limits beyond which extrapolation is expected to be
    // very bad, so raise an error.
    REAL8 q_max_hard = 6.01;

    if (q < 1) {
        XLAL_ERROR(XLAL_FAILURE, "Invalid mass ratio q = %0.4f < 1\n", q);
    }
    if ((q > q_max_hard) && (unlim_extrap == 0)) {
        XLAL_ERROR(XLAL_FAILURE,
            "Too much extrapolation in mass ratio; q = %0.4f > %0.4f\n",
            q, q_max_hard);
    }
    if (q > q_max_soft) {
        print_warning(
            "Extrapolating outside training range q = %0.4f > %0.4f\n",
            q, q_max_soft);
    }

    if (chiAmag > 1) {
        XLAL_ERROR(XLAL_FAILURE,
            "Invalid spin magnitude |chiA| = %0.4f > 1\n", chiAmag);
    }
    if (chiBmag > 1) {
        XLAL_ERROR(XLAL_FAILURE,
            "Invalid spin magnitude |chiB| = %0.4f > 1\n", chiBmag);
    }
    if (chiAmag > chi_max_soft) {
        print_warning(
            "Extrapolating outside training range |chiA| = %0.4f > %0.4f\n",
            chiAmag, chi_max_soft);
    }
    if (chiBmag > chi_max_soft) {
        print_warning(
            "Extrapolating outside training range |chiB| = %0.4f > %0.4f\n",
            chiBmag, chi_max_soft);
    }


    const REAL8 eta = q/(1.+q)/(1.+q);
    const REAL8 chi_wtAvg = (q*chiAz+chiBz)/(1.+q);
    const REAL8 chi_hat
        = (chi_wtAvg - 38.*eta/113.*(chiAz + chiBz))/(1. - 76.*eta/113.);
    const REAL8 chi_a = (chiAz - chiBz)/2.;

    XLAL_CHECK((fit_params != NULL) && (fit_params->size == 7), XLAL_EDIMS,
        "Size of fit_params should be 7, not %zu.\n", fit_params->size);

    gsl_vector_set(fit_params, 0, log(q));
    gsl_vector_set(fit_params, 1, chiAx);
    gsl_vector_set(fit_params, 2, chiAy);
    gsl_vector_set(fit_params, 3, chi_hat);
    gsl_vector_set(fit_params, 4, chiBx);
    gsl_vector_set(fit_params, 5, chiBy);
    gsl_vector_set(fit_params, 6, chi_a);

    return XLAL_SUCCESS;
}

//************************ main functions ***************************/

/**
 * Returns the remnant BH's mass, spin, or kick according to NRSur7dq4Remnant
 * model.
 *
 * Reference: arxiv:1905.09300.
 *
 * The input spins must be dimensionless, and given in the coorbital frame at
 * t=-100M from the total amplitude peak, as described in the paper. The
 * remnant mass is returned in units of total mass, the remnant spin is
 * dimensionless, and the remnant kick is in units of c. The remnant spin and
 * kick vectors are also returned in the same frame.
 *
 * remnant_property can be one of "mf", "chif" or "vf", respectively, for
 * the remnant mass, spin or kick.
 *
 * Usage examples:
 * mf = lalsim.NRSur7dq4Remnant(q, s1x, s1y, s1z, s2x, s2y, s2z, "mf").data[0]
 * chif = lalsim.NRSur7dq4Remnant(q, s1x, s1y, s1z, s2x, s2y, s2z, "chif").data
 */
int XLALNRSur7dq4Remnant(
    gsl_vector **result,        /**<Output: The requested remnant property. */
    REAL8 q,                    /**< Mass ratio of Bh1/Bh2. q>=1. */
    REAL8 s1x,                  /**< S1x in coorbital frame at t=-100M */
    REAL8 s1y,                  /**< S1y in coorbital frame at t=-100M */
    REAL8 s1z,                  /**< S1z in coorbital frame at t=-100M */
    REAL8 s2x,                  /**< S2x in coorbital frame at t=-100M */
    REAL8 s2y,                  /**< S2y in coorbital frame at t=-100M */
    REAL8 s2z,                  /**< S2z in coorbital frame at t=-100M */
    char *remnant_property,     /**< One of "mf", "chif" or "vf" */
    LALDict* LALparams          /**< Dict with extra parameters */
) {

#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&NRSur7dq4Remnant_is_initialized,
                      NRSur7dq4Remnant_Init_LALDATA);
#else
    NRSur7dq4Remnant_Init_LALDATA();
#endif

    // Loaded surrogate data
    const PrecessingRemnantFitData *sur_data = &__lalsim_NRSur7dq4Remnant_data;
    if (!sur_data->setup) {
        XLAL_ERROR(XLAL_EFAILED, "Error loading surrogate data.\n");
    }

    // assign size to dummy_worker
    gsl_vector *dummy_worker = gsl_vector_alloc(sur_data->params_dim);

    // Compute fit_params (initialize with dimension of the surrogate)
    gsl_vector *fit_params = gsl_vector_alloc(sur_data->params_dim);
    int ret = NRSur7dq4Remnant_fitParams(fit_params, q, s1x, s1y, s1z,
            s2x, s2y, s2z, LALparams);
    if(ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to evaluate fit_params.");
    }

    REAL8 tmp;
    // For final mass, we have a single scalar fit
    if (strcmp(remnant_property, "mf") == 0) {

        ScalarFitData *fit_data = sur_data->mf_data;

        *result = gsl_vector_alloc(1);
        tmp = NRHybSur_eval_fit(fit_data, fit_params, sur_data->x_train,
                dummy_worker);
        gsl_vector_set(*result, 0, tmp);

        return XLAL_SUCCESS;

    // For final spin and kick, we have a vector fit
    } else {

        VectorFitData *vec_data;
        if (strcmp(remnant_property, "chif") == 0) {
            vec_data = sur_data->chif_data;
        } else if (strcmp(remnant_property, "vf") == 0) {
            vec_data = sur_data->vf_data;
        } else {
            XLAL_ERROR(XLAL_EINVAL, "Invalid remnant_property, should be one "
                    "of 'mf', 'chif' or 'vf'");
        }

        *result = gsl_vector_alloc(vec_data->vec_dim);
        for (UINT4 i=0; i<vec_data->vec_dim; i++) {
            tmp = NRHybSur_eval_fit(vec_data->fit_data[i], fit_params,
                    sur_data->x_train, dummy_worker);
            gsl_vector_set(*result, i, tmp);
        }

        return XLAL_SUCCESS;
    }
}
