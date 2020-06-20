/*
 * Copyright (C) 2018 Vijay Varma
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
 * \brief C code for NRHybSur3dq8 waveform model, an NR-hybrid surrogate model
 * for aligned-spin BBH.
 *
 * The binary data file is available at https://dcc.ligo.org/LIGO-T1900034.
 * Get the lalsuite-extra repo or put the data into a location in your
 * LAL_DATA_PATH.
 *
 * **Paper**: https://arxiv.org/abs/1812.07865
 *
 * **Parameter ranges**:
 *
 *   q = [1, 10.1] and \f$\chi_{1z}, \chi_{2z}\f$ = [-0.81, 0.81]
 *   or
 *   q = [1, 9.1] and \f$\chi_{1z}, \chi_{2z}\f$ = [-0.91, 0.91]
 *
 *   modes: \f$ \ell \leq 4, m \geq 0 \f$, and (5,5), but not (4,1) or (4,0).
 *   m<0 modes are determined from the m \f$\geq0\f$ modes.
 *
 *   \f$M \geq 2.25 M_{\odot} \f$, for fstart=20Hz, for all modes.
 *
 * **Training parameter ranges**:
 *
 *   q = [1, 8]
 *
 *   \f$\chi_{1z}, \chi_{2z}\f$ = [-0.8, 0.8]
 *
 *   But extrapolates reasonably to the above mass ratios and spins.
 */


#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include "LALSimIMRNRHybSur3dq8.h"


#include <libgen.h>

#include <lal/FileIO.h>
#include <lal/SphericalHarmonics.h>
#include <lal/H5FileIO.h>


#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
#endif

#ifdef LAL_PTHREAD_LOCK
static pthread_once_t NRHybSur3dq8_is_initialized = PTHREAD_ONCE_INIT;
#endif

/**
 * Global surrogate data.
 * This data will be loaded at most once. Any executable which calls
 * NRHybSur3dq8_Init_LALDATA directly or by calling any XLAL NRHybSur3dq8
 * function will have a memory leak according to valgrind, because we never
 * free this memory.
 */
static NRHybSurData __lalsim_NRHybSur3dq8_data;

//*************************************************************************/
//************************* function definitions **************************/
//*************************************************************************/

/**
 * Helper function to check if the NRHybSur3dq8 model has been initialized.
 */
static bool NRHybSur3dq8_IsSetup(void) {
  if(__lalsim_NRHybSur3dq8_data.setup)
    return true;
  else
    return false;
}


/**
 * Surrogate initialization.
 *
 * This needs to be called once, before __lalsim_NRHybSur3dq8_data is used.  It
 * finds the H5 data file with the NRHybSur3dq8 data and loads the surrogate.
 * Can be called multiple times, will just return true if surrogate is already
 * loaded.
 */
static void NRHybSur3dq8_Init_LALDATA(void) {

    if (NRHybSur3dq8_IsSetup()) return;

    char *path = XLALFileResolvePathLong(NRHybSur3dq8_DATAFILE, PKG_DATA_DIR);
    if (path==NULL) {
        XLAL_ERROR_VOID(XLAL_ENOENT,
            "Unable to find data file %s in $LAL_DATA_PATH\n",
            NRHybSur3dq8_DATAFILE);
    }


    char *dir = dirname(path);
    const UINT4 size = strlen(dir) + strlen(NRHybSur3dq8_DATAFILE) + 2;
    char *file_path = XLALMalloc(size);
    snprintf(file_path, size, "%s/%s", dir, NRHybSur3dq8_DATAFILE);

    LALH5File *file = XLALH5FileOpen(file_path, "r");
    if (file==NULL) {
        XLAL_ERROR_VOID(XLAL_EIO,
            "Unable to load data file %s in $LAL_DATA_PATH."
            " File may be corrupted.\n", NRHybSur3dq8_DATAFILE);
    }

    int ret = NRHybSur_Init(&__lalsim_NRHybSur3dq8_data, file);
    XLALH5FileClose(file);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR_VOID(XLAL_FAILURE, "Failure loading data from %s\n",
                file_path);
    }

    XLALFree(path);
    XLALFree(file_path);
}


/**
 * Map from mass ratio and spins to surrogate fit parameters.
 *
 * The fit parameters are \f$[log_e(q), \hat{\chi}, \chi_a]\f$.
 * \f$\hat{\chi}\f$ is defined in Eq.(46) of arxiv:1812.07865.
 * \f$\chi_a = (\chi_{1z} - \chi_{2z})/2 \f$.
 * These are described in Sec.VI.C.3 of arxiv:1812.07865.
 */
int NRHybSur3dq8_fitParams(
    gsl_vector* fit_params,  /**< Output: mapped fit parameters. */
    const REAL8 q,          /**< Mass ratio m1 / m2 >= 1. */
    const REAL8 chi1z,      /**< Dimless spin of heavier BH. */
    const REAL8 chi2z       /**< Dimless spin of lighter BH. */
) {

    const REAL8 eta = q/(1.+q)/(1.+q);
    const REAL8 chi_wtAvg = (q*chi1z+chi2z)/(1.+q);
    const REAL8 chiHat
        = (chi_wtAvg - 38.*eta/113.*(chi1z + chi2z))/(1. - 76.*eta/113.);
    const REAL8 chi_a = (chi1z - chi2z)/2.;

    XLAL_CHECK((fit_params != NULL) && (fit_params->size == 3), XLAL_EDIMS,
        "NRHybSur3dq8_fitParams(): size of fit_params should be 3, not %zu.\n",
        fit_params->size);

    gsl_vector_set(fit_params, 0, log(q));
    gsl_vector_set(fit_params, 1, chiHat);
    gsl_vector_set(fit_params, 2, chi_a);

    return XLAL_SUCCESS;
}


/**
 * This is the core function of the NRHybSur3dq8 model.
 *
 * It evaluates all required waveform modes. For each mode, it evaluates each
 * waveform data piece. The different data pieces are described in Sec.VI.B of
 * arxiv:1812.07865.
 * For the (2,2) mode the data pieces are the amplitude and phase. Note that we
 * model the phase residual but add back the 0PN term at evaluation time.
 * For all other modes the data pieces are the real and imaginary parts of the
 * strain in the coorbital frame.
 *
 * The reference point is the time at which the (2,2) mode frequency equals
 * fRef. If fRef=0, sets fRef=fMin. We set the orbital phase to 0 at the
 * reference point. The orbital phase is obtained as \f$\phi_{22}/2\f$, so this
 * leaves a pi ambiguity. But the surrogate data is already aligned such that
 * the heavier BH is on the +ve x-axis at t=-1000M. See Sec.VI.A.4 of
 * arxiv:1812.07865.  This resolves the pi ambiguity. This means that after the
 * realignment, the orbital phase at reference frequency fRef is 0, or Bh1 is
 * on the +ve x-axis. Note that this is alignment is done using only waveform
 * quantities, so this doesn't necessarily correspond to the (gauge dependent)
 * NR Bh trajectories. The modes are returned in this frame, which agrees with
 * the LAL convention. When combining the modes to get the polarizations, the
 * Ylms should be evaluated at (inclination, pi/2 - phiRef), following the LAL
 * convention.
 *
 * Only uses data at (2,2) mode frequencies >= fMin. This determines the start
 * time. The start time, along with the step size deltaT, is used to determine
 * the output_times. Uses cubic spline interpolation to interpolate from the
 * surrogate's time array to output_times.
 *
 * NOTE: If mass ratio q<1, the labels of the two BHs are swapped internally.
 *
 * **Returned values**:
 *
 * phi_22: contains the phase of the (2,2) mode. This is always evaluated as
 *      this is required for other modes as well to transform from coorbital
 *      frame to inertial frame.
 *
 * evaluated_mode_dps: Contains all other data pieces. This is a list of size
 *      num_modes_incl, the number of modes to include. For each mode this
 *      contains the amplitude, and real and imaginary parts of the coorbital
 *      frame strain. For the (2,2) mode only the amplitude is required. For
 *      all other modes only the coorbital frame strain is required. So, we
 *      evaluate only the required data pieces of each mode. The amplitude and
 *      coorbital frame strain is in units of rh/M and needs to be rescaled to
 *      physical units.
 *
 * epoch: Initial time value, w.r.t. the peak (t=0 at the peak) of the total
 * waveform amplitude, as defined in Eq.38 of arxiv:1812.07865.
 */
int NRHybSur3dq8_core(
    gsl_vector **phi_22,  /**< Output: phase of (2,2) mode. */
    EvaluatedDataPieces **evaluated_mode_dps, /**< Output: All other data
                            pieces. */
    LIGOTimeGPS *epoch,     /**< Output: Initial time value, where t=0 is at
                                the peak of the total waveform amplitude. */
    const REAL8 deltaT,    /**< Sampling interval (s). */
    const REAL8 fMin,      /**< Start GW frequency (Hz). */
    const REAL8 fRef,      /**< Reference GW frequency (Hz). */
    REAL8 q,               /**< Mass ratio m1/m2. */
    const REAL8 Mtot_sec,  /**< Total mass in geometric units (s). */
    REAL8 chi1z,           /**< Dimensionless spin of Bh1. */
    REAL8 chi2z,           /**< Dimensionless spin of Bh2. */
    LALValue* ModeArray,   /**< Container for (ell, m) modes to generate. */
    LALDict* LALparams     /**< Dict with extra parameters */
) {

    const NRHybSurData *NR_hybsur_data = &__lalsim_NRHybSur3dq8_data;

    REAL8 init_orbphase = 0;    // In the LAL convention the larger BH should
                                // be on the +ve x-axis at fRef, this is done
                                // by setting orbphase = 0 at fRef.
    // For the surrogate Bh1 is defined to be the one with the larger mass,
    // Bh2 with the smaller mass. Therefore, if q<1, we swap the labels of the
    // two Bhs
    if (q < 1) {
        q = 1./q;
        REAL8 tmp = chi1z;
        chi1z = chi2z;
        chi2z = tmp;

        // The above swap generates the requested waveform but rotated by
        // pi about the z-axis, we undo this by adding pi to init_orbphase.
        init_orbphase += LAL_PI;
    }

    const char *param_validity = "This model is valid for q <= 9.1 & "
        "|chi1z,chi2z| <= 0.91, or q <= 10.1 & |chi1z,chi2z| <= 0.81";

    // By default we do not allow unlimited_extrapolation
    UINT4 unlim_extrap = 0;
    if (LALparams != NULL &&
            XLALDictContains(LALparams, "unlimited_extrapolation")) {
        // Unless the user asks for it
        unlim_extrap
            = XLALDictLookupUINT4Value(LALparams, "unlimited_extrapolation");
    }

    // Sanity checks and warnings
    if ((q > 10.1) && (unlim_extrap == 0)) {
        XLAL_ERROR(XLAL_EDOM,
            "Too much extrapolation in mass ratio; q=%0.4f > 10.1\n%s\n", q,
            param_validity);
    }
    if (q > 8.1) {
        XLAL_PRINT_WARNING(
            "Extrapolating outside training range q=%0.4f > 8.1\n", q);
    }
    if ((fabs(chi1z) > 0.91) && (unlim_extrap == 0)) {
        XLAL_ERROR(XLAL_EDOM,
            "Too much extrapolation; |chi1z|=%0.4f > 0.91\n%s\n", fabs(chi1z),
            param_validity);
    }
    if ((fabs(chi2z) > 0.91) && (unlim_extrap == 0)) {
        XLAL_ERROR(XLAL_EDOM,
            "Too much extrapolation; |chi2z|=%0.4f > 0.91\n%s\n", fabs(chi2z),
            param_validity);
    }
    if (fabs(chi1z) > 0.81) {
        if ((q > 9.1) && (unlim_extrap == 0)) {
            XLAL_ERROR(XLAL_EDOM,
                "Too much extrapolation; q=%0.4f > 9.1 & |chi1z|=%.04f"
                " >0.81\n%s\n", q, fabs(chi1z), param_validity);
        }
        XLAL_PRINT_WARNING(
            "Extrapolating outside training range |chi1z|=%0.4f > 0.81\n",
            fabs(chi1z));
    }
    if (fabs(chi2z) > 0.81) {
        if ((q > 9.1) && (unlim_extrap == 0)) {
            XLAL_ERROR(XLAL_EDOM,
                "Too much extrapolation; q=%0.4f > 9.1 & |chi2z|=%.04f"
                " >0.81\n%s\n", q, fabs(chi2z), param_validity);
        }
        XLAL_PRINT_WARNING(
            "Extrapolating outside training range |chi2z|=%0.4f > 0.81\n",
            fabs(chi2z));
    }

    // Get dimensionless start and reference frequencies. Note that cases where
    // fRef<fMin, deltaT is too small, and fRef or fMin are too
    // small/large (i.e. outside of the surrogate range) are handled in
    // NRHybSur_eval_phase_22().

    // dimensionless start angular frequency of (2,2) mode in rad/M
    const REAL8 omegaM_22_min = 2*LAL_PI*fMin*Mtot_sec;

    // dimensionless reference angular frequency of (2,2) mode in rad/M
    REAL8 omegaM_22_ref;
    if (fRef == 0) {
        // If fRef is 0, set it to fMin
        omegaM_22_ref = omegaM_22_min;
    }
    else {
        omegaM_22_ref = 2*LAL_PI*fRef*Mtot_sec;
    }

    // dimensionless time step size
    const REAL8 deltaTOverM = deltaT/Mtot_sec;

    // Compute fit_params (initialize with dimension of the surrogate)
    gsl_vector *fit_params = gsl_vector_alloc(NR_hybsur_data->params_dim);
    int ret = NRHybSur3dq8_fitParams(fit_params, q, chi1z, chi2z);
    if(ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to evaluate fit_params.");
    }

    // Training data parameters, required for GPR fits
    const gsl_matrix *x_train = NR_hybsur_data->x_train;

    // sparse time array of the surrogate
    const gsl_vector *domain = NR_hybsur_data->domain;

    // assign size to dummy_worker and dummy_dp
    gsl_vector *dummy_worker = gsl_vector_alloc(NR_hybsur_data->params_dim);
    gsl_vector *dummy_dp = gsl_vector_alloc(domain->size);

    gsl_vector *output_times = NULL;

    // symmetric mass ratio
    const REAL8 eta = q/(1.+q)/(1.+q);

    // Evaluate phase of (2,2) mode
    // The phase of the (2,2,) mode is always evaluated as this is
    // needed for the transformation from coorbital to inertial frame for all
    // modes
    ret = NRHybSur_eval_phase_22(phi_22, &output_times, eta, fit_params,
        omegaM_22_min, deltaTOverM, init_orbphase, omegaM_22_ref, dummy_dp,
        x_train, dummy_worker, NR_hybsur_data);
    if(ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC,
            "Failed to evaluate phi_22 data piece");
    }

    // set epoch to initial time. Note that t=0 is at the peak of the total
    // waveform amplitude, as defined in Eq.38 of arxiv:1812.07865.
    const REAL8 t0 = gsl_vector_get(output_times, 0);
    XLALGPSAdd(epoch, Mtot_sec * t0);

    // Evaluate other data pieces for required modes
    const gsl_matrix_long *mode_list = NR_hybsur_data->mode_list;
    const UINT4 num_modes_modeled = NR_hybsur_data->num_modes_modeled;
    UINT4 incl_mode_idx = 0;   // This tracks the output modes
    for (UINT4 mode_idx = 0; mode_idx < num_modes_modeled; mode_idx++){

        const UINT4 ell = gsl_matrix_long_get(mode_list, mode_idx, 0);
        const UINT4 m = gsl_matrix_long_get(mode_list, mode_idx, 1);

        // Evaluate a mode only if it is required
        if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, m) == 1) {

            const ModeDataPieces *data_pieces
                = NR_hybsur_data->mode_data_pieces[mode_idx];

            if((ell != data_pieces->ell) || (m != data_pieces->m)){
                XLAL_ERROR(XLAL_EDATA, "Modes do not agree");
            }

            evaluated_mode_dps[incl_mode_idx]
                = (EvaluatedDataPieces *)
                XLALMalloc(sizeof(EvaluatedDataPieces));

            ret = NRHybSur_eval_mode_data_pieces(
                &evaluated_mode_dps[incl_mode_idx], ell, m,
                data_pieces, output_times, fit_params, dummy_dp,
                x_train, dummy_worker, NR_hybsur_data);
            if(ret != XLAL_SUCCESS) {
                XLAL_ERROR(XLAL_EFUNC,
                    "Failed to evaluate (%u, %u) mode", ell, m);
            }

            incl_mode_idx += 1;
        }
    }

    gsl_vector_free(fit_params);
    gsl_vector_free(dummy_dp);
    gsl_vector_free(dummy_worker);
    gsl_vector_free(output_times);

    return XLAL_SUCCESS;
}






//************************ main functions ***************************/

/**
 * Reference: arxiv:1812.07865.
 *
 * Returns the plus and cross polarizations of NRHybSur3dq8 waveform model.
 *
 * The reference point is the time at which the (2,2) mode frequency equals
 * fRef. If fRef=0, sets fRef=fMin. We set the orbital phase to 0 at the
 * reference point. The orbital phase is obtained as \f$\phi_{22}/2\f$, so this
 * leaves a pi ambiguity. But the surrogate data is already aligned such that
 * the heavier BH is on the +ve x-axis at t=-1000M. See Sec.VI.A.4 of
 * arxiv:1812.07865.  This resolves the pi ambiguity. This means that after the
 * realignment, the orbital phase at reference frequency fRef is 0, or Bh1 is
 * on the +ve x-axis. Note that this is alignment is done using only waveform
 * quantities, so this doesn't necessarily correspond to the (gauge dependent)
 * NR Bh trajectories. The polarizations of the waveform are returned in the sky
 * of this frame at (inclination, pi/2 - phiRef). This agrees with the LAL
 * convention.
 *
 * Only uses data at (2,2) mode frequencies >= fMin. This determines the start
 * time. The start time, along with the step size deltaT, is used to determine
 * the output_times. Uses cubic spline interpolation to interpolate from the
 * surrogate's time array to output_times.
 *
 * By default, uses all available modes of the surrogate, that is \f$ \ell \leq
 * 4, m \geq 0 \f$, and (5,5), but not (4,1) or (4,0).  For m>0 modes, the
 * contribution from the m<0 counterparts is automatically added.
 *
 * If a custom ModeArray is given, only those modes are used. Note that it only
 * looks for m>=0 modes in ModeArray, and will ignore m<0 modes even if present.
 * The m<0 modes automatically get added.
 *
 * This surrogate model is trained on the following range.
 *   q <= 8, |chi_1|, |chi_2| <= 0.8
 *   If you want a guarantee of accuracy you should stick to the above ranges.
 *
 *   We allow extrapolation to the following ranges, but with a warning:
 *   q = [1, 10.1] and \f$\chi_{1z}, \chi_{2z}\f$ = [-0.81, 0.81]
 *   or
 *   q = [1, 9.1] and \f$\chi_{1z}, \chi_{2z}\f$ = [-0.91, 0.91]
 *   We expect the model to be reasonable when extrapolated to these ranges.
 *
 *   Beyond the above ranges, we raise an error. If you want to extrapolate
 *   beyond these limits you can specify unlimited_extrapolation = 1 in your
 *   dictParams as follows:
 *       # USE AT YOUR OWN RISK!!
 *       lal.DictInsertUINT4Value(dictParams, "unlimited_extrapolation", 1)
 */
INT4 XLALSimIMRNRHybSur3dq8Polarizations(
    REAL8TimeSeries **hplus,        /**<Output: \f$h_+\f$ polarization. */
    REAL8TimeSeries **hcross,       /**<Output: \f$h_{\times}\f$ polarization.*/
    REAL8 phiRef,                   /**< azimuthal angle for Ylms */
    REAL8 inclination,              /**< Inclination angle. */
    REAL8 deltaT,                   /**< Sampling interval (s). */
    REAL8 m1,                       /**< Mass of Bh1 (kg). */
    REAL8 m2,                       /**< Mass of Bh2 (kg). */
    REAL8 distance,                 /**< Distance of source (m). */
    REAL8 fMin,                     /**< Start GW frequency (Hz). */
    REAL8 fRef,                     /**< Reference GW frequency (Hz). */
    REAL8 chi1z,                    /**< Dimensionless spin of Bh1. */
    REAL8 chi2z,                    /**< Dimensionless spin of Bh2. */
    LALDict* LALparams              /**< Dict with extra parameters */
)
{
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&NRHybSur3dq8_is_initialized, NRHybSur3dq8_Init_LALDATA);
#else
    NRHybSur3dq8_Init_LALDATA();
#endif

    // Loaded surrogate data
    const NRHybSurData *NR_hybsur_data = &__lalsim_NRHybSur3dq8_data;

    if (NR_hybsur_data->setup != 1){
        XLAL_ERROR(XLAL_FAILURE, "Surrogate data is not loaded.");
    }

    // If ModeArray is not specified, use all available modes
    LALValue* ModeArray
        = XLALSimInspiralWaveformParamsLookupModeArray(LALparams);
    if (ModeArray == NULL) {
        ModeArray = XLALSimInspiralCreateModeArray();
        NRHybSur_set_default_modes(ModeArray, NR_hybsur_data);
    }

    // Make sure we didn't request any unavailable modes, and get number of
    // modes to include
    UINT4 num_modes_incl, max_ell;
    int ret = NRHybSur_check_mode_array(&num_modes_incl, &max_ell, ModeArray,
            NR_hybsur_data);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Inappropriate ModeArray.");
    }

    // Check Nyquist frequency is larger than the ringdow frequency
    // of the (max_ell,max_ell,0) QNM overtone, where max_ell is the
    // highest ell index included. Raises warning only, not error.
    ret = NRHybSur_sanity_check_sample_rate(deltaT, m1, m2, chi1z, chi2z,
        max_ell);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "check_sample_rate failed.");
    }

    // Total mass in seconds
    const REAL8 Mtot_sec = (m1 + m2)/LAL_MSUN_SI * LAL_MTSUN_SI;

    // mass ratio
    const REAL8 q = m1/m2;

    // rescale factor for the amplitude based on the distance and total mass
    const REAL8 a0 = Mtot_sec * LAL_C_SI / distance;

    // Initialize vector to hold phase of the (2,2,) mode
    gsl_vector *phi_22 = NULL;

    // Initialize data pieces for amplitude of (2,2) mode, and coorbital
    // frame strain for other modes.
    EvaluatedDataPieces **evaluated_mode_dps
        = XLALMalloc(num_modes_incl * sizeof(*evaluated_mode_dps));

    LIGOTimeGPS epoch = LIGOTIMEGPSZERO;

    // Evaluate the surrogate
    ret = NRHybSur3dq8_core(&phi_22, evaluated_mode_dps,
        &epoch, deltaT, fMin, fRef, q, Mtot_sec, chi1z, chi2z,
        ModeArray, LALparams);
    if(ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Surrogate evaluation failed.");
    }

    const UINT4 num_times = phi_22->size;

    // initialize the hp and hc vectors to 0
    REAL8TimeSeries *hPlusTS =
      XLALCreateREAL8TimeSeries("H_PLUS", &epoch, 0.0, deltaT, &lalStrainUnit,
            num_times);
    REAL8TimeSeries *hCrossTS =
      XLALCreateREAL8TimeSeries("H_CROSS", &epoch, 0.0, deltaT, &lalStrainUnit,
            num_times);
    memset(hPlusTS->data->data, 0, hPlusTS->data->length*sizeof(REAL8));
    memset(hCrossTS->data->data, 0, hCrossTS->data->length*sizeof(REAL8));

    // Sum up all required modes
    COMPLEX16 swsh_coef;
    COMPLEX16 swsh_coef_negM = 0;   // some compilers complain if this is
                                    // not initialized here
    COMPLEX16 hlm;
    COMPLEX16 hcomplex;
    int pre_factor;
    UINT4 ell, m;
    for (UINT4 mode_idx = 0; mode_idx < num_modes_incl; mode_idx++){

        EvaluatedDataPieces *this_mode_eval_dp = evaluated_mode_dps[mode_idx];
        ell = this_mode_eval_dp->ell;
        m = this_mode_eval_dp->m;

        /* See Harald Pfeiffer, T18002260-v1 for frame diagram. The surrogate
         * frame has its z direction aligned with the orbital angular momentum,
         * which agrees with the lowercase source frame z in the diagram. The
         * surrogate x direction, however, points along the line of ascending
         * nodes (the funny omega with circles on the ends). The uppercase Z,
         * which is the direction along which we want to evaluate the waveform,
         * is always in the surrogate frame (y, z) plane. Z is rotated from z
         * towards the +ive y surrogate axis, so we should always evaluate at
         * (inclination, pi/2-phiRef). */
        swsh_coef = XLALSpinWeightedSphericalHarmonic(inclination,
                0.5 * LAL_PI - phiRef, -2, ell, m);

        if (m > 0) {
            // Note: Each time a new mode is added to the strain, the value of
            // swsh_coef_negM persists. When m>0 the old value is overwritten,
            // while if m=0 it has a meaningless value from the previous
            // iteration. This is not a bug, but can be confusing, so just
            // adding this note for anyone that is debugging this code.
            swsh_coef_negM = XLALSpinWeightedSphericalHarmonic(inclination,
                    0.5 * LAL_PI - phiRef, -2, ell, -m);
        }

        for (UINT4 j=0; j < num_times; j++) {

            // The (2,2) phase is required for all modes; for the other modes,
            // we use it to transform from the coorbital frame to the inertial
            // frame
            const REAL8 tmp_phi_22 = gsl_vector_get(phi_22, j);

            if (ell == 2 && m ==2){
                const REAL8 tmp_Amp_22 = gsl_vector_get(
                    this_mode_eval_dp->ampl_eval, j);
                // (2,2) mode strain in the inertial frame
                hlm = tmp_Amp_22 * (cos(tmp_phi_22) - I*sin(tmp_phi_22));

            } else {

                // set the real and imaginary parts (in coorbital frame) to
                // zero by default
                REAL8 coorb_hlm_re = 0.0;
                REAL8 coorb_hlm_im = 0.0;

                // for m=0, l=odd, the real part is zero. For all other
                // cases get the real part
                if (m != 0 || ell % 2 == 0) {
                    coorb_hlm_re = gsl_vector_get(
                        this_mode_eval_dp->coorb_re_eval, j);
                }

                // for m=0, l=even, the imaginary part is zero. For all other
                // cases get the imaginary part
                if (m != 0 || ell % 2 == 1) {
                    coorb_hlm_im = gsl_vector_get(
                        this_mode_eval_dp->coorb_im_eval, j);
                }

                // strain in the inertial frame for all non (2,2) modes.
                // See Eq.(39) of arxiv:1812.07865.for definition of coorbital
                // frame strain
                hlm = (coorb_hlm_re + I * coorb_hlm_im)
                        * (cos(m*tmp_phi_22/2.) - I*sin(m*tmp_phi_22/2.));
            }

            hcomplex = hlm * swsh_coef;
            if (m > 0) {
                // For m>0 modes, also add contribution from m<0 modes as
                // (-1)**l * conj(hlm) * Y_{ell, -m}
                if (ell % 2 == 0){
                    pre_factor = 1;
                }
                else {
                    pre_factor = -1;
                }
                hcomplex += pre_factor * conj(hlm) * swsh_coef_negM;
            }

            // hcomplex = h_+ - i*h_x
            hPlusTS->data->data[j] += creal(hcomplex) * a0;
            hCrossTS->data->data[j] -= cimag(hcomplex) * a0;
        }
    }

    // Point the output pointers to the relevant time series and return
    (*hplus) = hPlusTS;
    (*hcross) = hCrossTS;

    // Cleanup
    NRHybSur_DestroyEvaluatedDataPieces(phi_22, evaluated_mode_dps,
        num_modes_incl);
    if(ModeArray) {
        XLALDestroyValue(ModeArray);
    }


    return XLAL_SUCCESS;
}



/**
 * Reference: arxiv:1812.07865.
 *
 * Returns the spin-weighted spherical harmonic modes of NRHybSur3dq8 waveform
 * model.
 *
 * The reference point is the time at which the (2,2) mode frequency equals
 * fRef. If fRef=0, sets fRef=fMin. We set the orbital phase to 0 at the
 * reference point. The orbital phase is obtained as \f$\phi_{22}/2\f$, so this
 * leaves a pi ambiguity. But the surrogate data is already aligned such that
 * the heavier BH is on the +ve x-axis at t=-1000M. See Sec.VI.A.4 of
 * arxiv:1812.07865.  This resolves the pi ambiguity. This means that after the
 * realignment, the orbital phase at reference frequency fRef is 0, or Bh1 is
 * on the +ve x-axis. Note that this is alignment is done using only waveform
 * quantities, so this doesn't necessarily correspond to the (gauge dependent)
 * NR Bh trajectories. The modes are returned in this frame, which agrees with
 * the LAL convention. When combining the modes to get the polarizations, the
 * Ylms should be evaluated at (inclination, pi/2 - phiRef), following the LAL
 * convention.
 *
 * Only uses data at (2,2) mode frequencies >= fMin. This determines the start
 * time. The start time, along with the step size deltaT, is used to determine
 * the output_times. Uses cubic spline interpolation to interpolate from the
 * surrogate's time array to output_times.
 *
 * By default, returns all available modes of the surrogate, that is \f$ \ell
 * \leq 4, m \geq 0 \f$, and (5,5), but not (4,1) or (4,0).
 *
 * If a custom ModeArray is given, only those modes are used. Note that it only
 * looks for m>=0 modes in ModeArray, and will ignore m<0 modes even if present.
 *
 * This surrogate model is trained on the following range.
 *   q <= 8, |chi_1|, |chi_2| <= 0.8
 *   If you want a guarantee of accuracy you should stick to the above ranges.
 *
 *   We allow extrapolation to the following ranges, but with a warning:
 *   q = [1, 10.1] and \f$\chi_{1z}, \chi_{2z}\f$ = [-0.81, 0.81]
 *   or
 *   q = [1, 9.1] and \f$\chi_{1z}, \chi_{2z}\f$ = [-0.91, 0.91]
 *   We expect the model to be reasonable when extrapolated to these ranges.
 *
 *   Beyond the above ranges, we raise an error. If you want to extrapolate
 *   beyond these limits you can specify unlimited_extrapolation = 1 in your
 *   dictParams as follows:
 *       # USE AT YOUR OWN RISK!!
 *       lal.DictInsertUINT4Value(dictParams, "unlimited_extrapolation", 1)
 */
SphHarmTimeSeries *XLALSimIMRNRHybSur3dq8Modes(
    REAL8 deltaT,                   /**< Sampling interval (s). */
    REAL8 m1,                       /**< Mass of Bh1 (kg). */
    REAL8 m2,                       /**< Mass of Bh2 (kg). */
    REAL8 chi1z,                    /**< Dimensionless spin of Bh1. */
    REAL8 chi2z,                    /**< Dimensionless spin of Bh2. */
    REAL8 fMin,                     /**< Start GW frequency (Hz). */
    REAL8 fRef,                     /**< Reference GW frequency (Hz). */
    REAL8 distance,                 /**< Distance of source (m). */
    LALDict* LALparams              /**< Dict with extra parameters */
) {
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&NRHybSur3dq8_is_initialized, NRHybSur3dq8_Init_LALDATA);
#else
    NRHybSur3dq8_Init_LALDATA();
#endif

    // Loaded surrogate data
    const NRHybSurData *NR_hybsur_data = &__lalsim_NRHybSur3dq8_data;

    if (NR_hybsur_data->setup != 1){
        XLAL_ERROR_NULL(XLAL_FAILURE, "Surrogate data is not loaded.");
    }

    // If ModeArray is not specified, use all available modes
    LALValue* ModeArray
        = XLALSimInspiralWaveformParamsLookupModeArray(LALparams);
    if (ModeArray == NULL) {
        ModeArray = XLALSimInspiralCreateModeArray();
        NRHybSur_set_default_modes(ModeArray, NR_hybsur_data);
    }

    // Make sure we didn't request any unavailable modes, and get number of
    // modes to include
    UINT4 num_modes_incl, max_ell;
    int ret = NRHybSur_check_mode_array(&num_modes_incl, &max_ell, ModeArray,
            NR_hybsur_data);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR_NULL(XLAL_EFUNC, "Inappropriate ModeArray.");
    }

    // Check Nyquist frequency is larger than the ringdow frequency
    // of the (max_ell,max_ell,0) QNM overtone, where max_ell is the
    // highest ell index included. Raises warning only, not error.
    ret = NRHybSur_sanity_check_sample_rate(deltaT, m1, m2, chi1z, chi2z,
        max_ell);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR_NULL(XLAL_EFUNC, "check_sample_rate failed.");
    }

    // Total mass in seconds
    const REAL8 Mtot_sec = (m1 + m2)/LAL_MSUN_SI * LAL_MTSUN_SI;

    // mass ratio
    const REAL8 q = m1/m2;

    // rescale factor for the amplitude based on the distance and total mass
    const REAL8 a0 = Mtot_sec * LAL_C_SI / distance;

    // Initialize vector to hold phase of the (2,2,) mode
    gsl_vector *phi_22 = NULL;

    // Initialize data pieces for amplitude of (2,2) mode, and coorbital
    // frame strain for other modes.
    EvaluatedDataPieces **evaluated_mode_dps
        = XLALMalloc(num_modes_incl * sizeof(*evaluated_mode_dps));

    LIGOTimeGPS epoch = LIGOTIMEGPSZERO;

    // Evaluate the surrogate
    ret = NRHybSur3dq8_core(&phi_22, evaluated_mode_dps,
        &epoch, deltaT, fMin, fRef, q, Mtot_sec, chi1z, chi2z,
        ModeArray, LALparams);
    if(ret != XLAL_SUCCESS) {
        XLAL_ERROR_NULL(XLAL_EFUNC, "Surrogate evaluation failed.");
    }

    const UINT4 num_times = phi_22->size;

    // Evaluate all required modes
    SphHarmTimeSeries *hlms = NULL;
    COMPLEX16TimeSeries *tmp_mode;
    UINT4 ell, m;
    char mode_name[32];
    for (UINT4 mode_idx = 0; mode_idx < num_modes_incl; mode_idx++){

        EvaluatedDataPieces *this_mode_eval_dp = evaluated_mode_dps[mode_idx];
        ell = this_mode_eval_dp->ell;
        m = this_mode_eval_dp->m;

        snprintf(mode_name, sizeof(mode_name), "(%u, %u) mode", ell, m);
        tmp_mode = XLALCreateCOMPLEX16TimeSeries(mode_name, &epoch, 0.0,
                deltaT, &lalStrainUnit, num_times);

        for (UINT4 j=0; j < num_times; j++) {

            const REAL8 tmp_phi_22 = gsl_vector_get(phi_22, j);

            if (ell == 2 && m ==2){
                const REAL8 tmp_Amp_22 = gsl_vector_get(
                    this_mode_eval_dp->ampl_eval, j);

                // (2,2) mode strain in the inertial frame
                tmp_mode->data->data[j]
                    = a0* tmp_Amp_22 * (cos(tmp_phi_22) - I*sin(tmp_phi_22));

            } else {

                // set the real and imaginary parts (in coorbital frame) to
                // zero by default
                REAL8 coorb_hlm_re = 0.0;
                REAL8 coorb_hlm_im = 0.0;

                // for m=0, l=odd, the real part is zero. For all other
                // cases get the real part
                if (m != 0 || ell % 2 == 0) {
                    coorb_hlm_re = gsl_vector_get(
                        this_mode_eval_dp->coorb_re_eval, j);
                }

                // for m=0, l=even, the imaginary part is zero. For all other
                // cases get the imaginary part
                if (m != 0 || ell % 2 == 1) {
                    coorb_hlm_im = gsl_vector_get(
                        this_mode_eval_dp->coorb_im_eval, j);
                }

                // strain in the inertial frame for all non (2,2) modes.
                // See Eq.(39) of arxiv:1812.07865.for definition of coorbital
                // frame strain
                tmp_mode->data->data[j]
                    = a0 * (coorb_hlm_re + I * coorb_hlm_im)
                        * (cos(m*tmp_phi_22/2.) - I*sin(m*tmp_phi_22/2.));
            }
        }

        hlms = XLALSphHarmTimeSeriesAddMode(hlms, tmp_mode, ell, m);
        XLALDestroyCOMPLEX16TimeSeries(tmp_mode);
    }

    // Cleanup
    NRHybSur_DestroyEvaluatedDataPieces(phi_22, evaluated_mode_dps,
        num_modes_incl);
    if(ModeArray) {
        XLALDestroyValue(ModeArray);
    }

    return hlms;
}
