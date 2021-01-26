/*  Copyright (C) 2019 Vijay Varma
 *  Evaluates NRSur3dq8Remnant model for remnant BH mass, spin and recoil kick
 *  for aligned-spin BBH.
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

/**
 * \author Vijay Varma
 *
 * \file
 *
 * \brief NRSur3dq8Remnant model for remnant BH mass, spin and recoil kick for
 * aligned-spin BBH.
 *
 * The binary data file is available at:
 * https://dcc.ligo.org/LIGO-T1900034/public.
 * Get the lalsuite-extra repo or put the data into a location in your
 * LAL_DATA_PATH.
 *
 * **Paper**: https://arxiv.org/abs/1809.09125,
 *    https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.011101.
 *    The model is referred to as surfinBH3dq8 in the paper.
 *
 * **Parameter ranges**:
 *
 *   q = [1, 9.1]
 *   \f$\chi_{1z}, \chi_{2z}\f$ = [-0.91, 0.91]
 *
 *   OR
 *
 *   q = [1, 10.1]
 *   \f$\chi_{1z}, \chi_{2z}\f$ = [-0.81, 0.81]
 *
 * **Training parameter ranges**:
 *
 *   q = [1, 8]
 *   \f$\chi_{1z}, \chi_{2z}\f$ = [-0.81, 0.81]
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

#include "LALSimNRSur3dq8Remnant.h"



#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
#endif

#ifdef LAL_PTHREAD_LOCK
static pthread_once_t NRSur3dq8Remnant_is_initialized = PTHREAD_ONCE_INIT;
#endif


/**
 * Global surrogate data. This data will be loaded at most once. Any
 * executable which calls NRSur3dq8Remnant_Init_LALDATA directly or by calling
 * the XLALNRSur3dq8Remnant function will have a memory leak according to
 * valgrind, because we never free this memory. This memory, however, gets freed
 * after the executable terminates.
 */
static AlignedSpinRemnantFitData __lalsim_NRSur3dq8Remnant_data;

//*************************************************************************/
//************************* function definitions **************************/
//*************************************************************************/

/**
 * Helper function to check if the NRSur3dq8Remnant model has been initialized.
 */
static bool NRSur3dq8Remnant_IsSetup(void) {
  if(__lalsim_NRSur3dq8Remnant_data.setup)
    return true;
  else
    return false;
}

/**
 * Surrogate initialization.
 *
 * This needs to be called once, before __lalsim_NRSur3dq8Remnant_data is used.
 * It finds the H5 data file with the NRSur3dq8Remnant data and loads the
 * surrogate.  Can be called multiple times, will just return true if surrogate
 * is already loaded.
 */
static void NRSur3dq8Remnant_Init_LALDATA(void) {

    if (NRSur3dq8Remnant_IsSetup()) return;

    LALH5File *file = NULL;
    NRSurRemnant_LoadH5File(&file, NRSur3dq8Remnant_DATAFILE);

    int ret = AlignedSpinNRSurRemnant_Init(&__lalsim_NRSur3dq8Remnant_data,
            file);

    XLALH5FileClose(file);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR_VOID(XLAL_FAILURE, "Failure loading data from %s\n",
                NRSur3dq8Remnant_DATAFILE);
    }
}

/**
 * Map from mass ratio and spins to surrogate fit parameters.
 *
 * The fit parameters are \f$[log_e(q), \hat{\chi}, \chi_a]\f$.
 * \f$\hat{\chi}\f$ is defined in Eq.(S5) of arxiv:1809.09125.
 * \f$\chi_a = (\chi_{1z} - \chi_{2z})/2 \f$.
 */
static int NRSur3dq8Remnant_fitParams(
    gsl_vector* fit_params, /**< Output: mapped fit parameters. */
    const REAL8 q,          /**< Mass ratio m1 / m2 >= 1. */
    const REAL8 chiAz,      /**< Dimless z-spin of heavier BH. */
    const REAL8 chiBz,      /**< Dimless z-spin of lighter BH. */
    LALDict* LALparams      /**< Dict with extra parameters */
) {

    //// Sanity check parameter ranges
    const char *param_validity = "This model is valid for q <= 9.1 & "
        "|chiAz,chiBz| <= 0.91, or q <= 10.1 & |chiAz,chiBz| <= 0.81";

    // By default we do not allow unlimited_extrapolation
    UINT4 unlim_extrap = 0;
    if (LALparams != NULL &&
            XLALDictContains(LALparams, "unlimited_extrapolation")) {
        // Unless the user asks for it
        unlim_extrap
            = XLALDictLookupUINT4Value(LALparams, "unlimited_extrapolation");
    }

    if (q < 1) {
        XLAL_ERROR(XLAL_FAILURE, "Invalid mass ratio q = %0.4f < 1\n", q);
    }
    if (fabs(chiAz) > 1) {
        XLAL_ERROR(XLAL_FAILURE,
            "Invalid spin magnitude |chiA| = %0.4f > 1\n", fabs(chiAz));
    }
    if (fabs(chiBz) > 1) {
        XLAL_ERROR(XLAL_FAILURE,
            "Invalid spin magnitude |chiB| = %0.4f > 1\n", fabs(chiBz));
    }
    if ((q > 10.1) && (unlim_extrap == 0)) {
        XLAL_ERROR(XLAL_EINVAL,
            "Too much extrapolation in mass ratio; q=%0.4f > 10.1\n%s\n", q,
            param_validity);
    }
    if (q > 8) {
        print_warning(
            "Extrapolating outside training range q=%0.4f > 8\n", q);
    }
    if (fabs(chiAz) > 0.91 && (unlim_extrap == 0)) {
        XLAL_ERROR(XLAL_EINVAL,
            "Too much extrapolation; |chiAz|=%0.4f > 0.91\n%s\n", fabs(chiAz),
            param_validity);
    }
    if (fabs(chiBz) > 0.91 && (unlim_extrap == 0)) {
        XLAL_ERROR(XLAL_EINVAL,
            "Too much extrapolation; |chiBz|=%0.4f > 0.91\n%s\n", fabs(chiBz),
            param_validity);
    }
    if (fabs(chiAz) > 0.81) {
        if ((q > 9.1) && (unlim_extrap == 0)) {
            XLAL_ERROR(XLAL_EINVAL,
                "Too much extrapolation; q=%0.4f > 9.1 & |chiAz|=%.04f"
                " >0.81\n%s\n", q, fabs(chiAz), param_validity);
        }
        print_warning(
            "Extrapolating outside training range |chiAz|=%0.4f > 0.81\n",
            fabs(chiAz));
    }
    if (fabs(chiBz) > 0.81) {
        if ((q > 9.1) && (unlim_extrap == 0)) {
            XLAL_ERROR(XLAL_EINVAL,
                "Too much extrapolation; q=%0.4f > 9.1 & |chiBz|=%.04f"
                " >0.81\n%s\n", q, fabs(chiBz), param_validity);
        }
        print_warning(
            "Extrapolating outside training range |chiBz|=%0.4f > 0.81\n",
            fabs(chiBz));
    }

    const REAL8 eta = q/(1.+q)/(1.+q);
    const REAL8 chi_wtAvg = (q*chiAz+chiBz)/(1.+q);
    const REAL8 chi_hat
        = (chi_wtAvg - 38.*eta/113.*(chiAz + chiBz))/(1. - 76.*eta/113.);
    const REAL8 chi_a = (chiAz - chiBz)/2.;

    XLAL_CHECK((fit_params != NULL) && (fit_params->size == 3), XLAL_EDIMS,
        "Size of fit_params should be 3, not %zu.\n", fit_params->size);

    gsl_vector_set(fit_params, 0, log(q));
    gsl_vector_set(fit_params, 1, chi_hat);
    gsl_vector_set(fit_params, 2, chi_a);

    return XLAL_SUCCESS;
}

//************************ main functions ***************************/

/**
 * Returns the remnant BH's mass, spin, or kick according to NRSur3dq8Remnant
 * model.
 *
 * Reference: arxiv:1809.09125. This model is referred to as surfinBH3dq8 in
 * the paper.
 *
 * The remnant mass is returned in units of total mass, the remnant spin is
 * dimensionless, and the recoil kick is in units of c.
 *
 * remnant_property can be one of "mf", "chifz", "vfx" or "vfy", respectively,
 * for the remnant mass, the z-component of remnant spin, and the x and y
 * components of the kick. The kick components are given in the coorbital frame
 * at t=-100M from the total amplitude peak, as described in the paper. All
 * other components of the spin and kick are zero due to the symmetries of
 * aligned-spin systems and are not modelled by the surrogate.
 *
 * Usage examples:
 * mf = lalsim.NRSur3dq8Remnant(q, s1z, s2z, "mf")
 * chifz = lalsim.NRSur3dq8Remnant(q, s1z, s2z, "chifz")
 */
int XLALNRSur3dq8Remnant(
    REAL8 *result,              /**<Output: The requested remnant property. */
    REAL8 q,                    /**< Mass ratio of Bh1/Bh2. q>=1. */
    REAL8 s1z,                  /**< S1z z-spin of Bh1 */
    REAL8 s2z,                  /**< S2z z-spin of Bh2 */
    char *remnant_property,     /**< One of "mf", "chifz", "vfx" or "vfy" */
    LALDict* LALparams          /**< Dict with extra parameters */
) {

#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&NRSur3dq8Remnant_is_initialized,
                      NRSur3dq8Remnant_Init_LALDATA);
#else
    NRSur3dq8Remnant_Init_LALDATA();
#endif

    // Loaded surrogate data
    const AlignedSpinRemnantFitData *sur_data = &__lalsim_NRSur3dq8Remnant_data;
    if (!sur_data->setup) {
        XLAL_ERROR(XLAL_EFAILED, "Error loading surrogate data.\n");
    }

    // assign size to dummy_worker
    gsl_vector *dummy_worker = gsl_vector_alloc(sur_data->params_dim);

    // Compute fit_params (initialize with dimension of the surrogate)
    gsl_vector *fit_params = gsl_vector_alloc(sur_data->params_dim);
    int ret = NRSur3dq8Remnant_fitParams(fit_params, q, s1z, s2z, LALparams);
    if(ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to evaluate fit_params.");
    }

    // All quantities are treated as scalars
    ScalarFitData *fit_data = NULL;
    if (strcmp(remnant_property, "mf") == 0) {
        // final mass
        fit_data = sur_data->mf_data;
    } else if (strcmp(remnant_property, "chifz") == 0) {
        // z-component of final spin (other components are zero)
        fit_data = sur_data->chifz_data;
    } else if (strcmp(remnant_property, "vfx") == 0) {
        // x-component of recoil kick
        fit_data = sur_data->vfx_data;
    } else if (strcmp(remnant_property, "vfy") == 0) {
        // y-component of recoil kick, z-component is zero
        fit_data = sur_data->vfy_data;
    } else {
        XLAL_ERROR(XLAL_EINVAL, "Invalid remnant_property, should be one "
                "of 'mf', 'chifz', 'vfx' or 'vfy'");
    }

    *result = NRHybSur_eval_fit(fit_data, fit_params, sur_data->x_train,
            dummy_worker);
    return XLAL_SUCCESS;

}
