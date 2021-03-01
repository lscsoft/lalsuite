/*
 * Copyright (C) 2019 Jonathan Thompson, Edward Fauchon-Jones, Sebastian Khan
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

#include "LALSimIMRPhenomNSBH.h"

/* This is ugly, but allows us to reuse internal IMRPhenomC and IMRPhenomD functions without making those functions XLAL */
#include "LALSimIMRPhenomC_internals.c"
#include "LALSimIMRPhenomD_internals.c"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_poly.h>

#ifndef _OPENMP
#define omp ignore
#endif


int IMRPhenomNSBH_Core(
    COMPLEX16FrequencySeries **htilde, /**< Output: Frequency-domain waveform h+ */
    REAL8 phiRef,                      /**< Phase at reference time */
    REAL8 fRef,                        /**< Reference frequency (Hz); 0 defaults to fLow */
    REAL8 distance,                    /**< Distance of source (m) */
    REAL8 mBH_SI,                      /**< Mass of BH (kg) */
    REAL8 mNS_SI,                      /**< Mass of neutron star 2 (kg) */
    REAL8 chi_BH,                      /**< Dimensionless aligned component spin of Black Hole */
    REAL8 chi_NS,                      /**< Dimensionless aligned component spin of NS */
    LALDict *extraParams,              /**< extra params */
    const REAL8Sequence *freqs_in,     /**< Frequency points at which to evaluate the waveform (Hz) */
    REAL8 deltaF)                      /**< Sampling frequency (Hz) */
{

    XLAL_CHECK(mBH_SI >= mNS_SI, XLAL_EDOM, "BH mass must be larger or equal to NS mass\n");
    XLAL_CHECK(mBH_SI <= 100.0*mNS_SI, XLAL_EDOM, "Mass ratio must be less than or equal to 100");
    if (chi_NS != 0.0)
        XLALPrintWarning("IMRPhenomNSBH is not calibrated for a non-zero NS spin\n");

    /*sanity checks*/
    /* Check output arrays */
    if (!htilde)
        XLAL_ERROR(XLAL_EFAULT, "htilde must not be a NULL pointer\n");
    if (*htilde)
    {
        XLALPrintError("(*htilde) is supposed to be NULL, but got %p", (*htilde));
        XLAL_ERROR(XLAL_EFAULT);
    }

    /* Find frequency bounds */
    if (!freqs_in)
        XLAL_ERROR(XLAL_EFAULT, "freqs_in must not be a NULL pointer\n");
    double f_min = freqs_in->data[0];
    double f_max = freqs_in->data[freqs_in->length - 1];
    XLAL_CHECK(f_min > 0, XLAL_EDOM, "Minimum frequency must be positive.\n");
    XLAL_CHECK(f_max >= 0, XLAL_EDOM, "Maximum frequency must be non-negative.\n");
    if (f_max > 0)
        XLAL_CHECK(f_max >= f_min, XLAL_EDOM, "Maximum frequency must not be less than minimum frequency.\n");
    if (fRef == 0.0)
        fRef = f_min;

    IMRPhenomDPhaseCoefficients *pPhi = NULL;
    PNPhasingSeries *pn = NULL;
    REAL8Sequence *freqs = NULL;
    LALDict *extraParams_in = extraParams;

    int retcode = XLALSimInspiralSetQuadMonParamsFromLambdas(extraParams);
    XLAL_CHECK(retcode == XLAL_SUCCESS, XLAL_EFUNC, "Failed to set quadparams from Universal relation.\n");
    REAL8 lambda_NS = XLALSimInspiralWaveformParamsLookupTidalLambda2(extraParams);
    XLAL_CHECK(lambda_NS <= 5000, XLAL_EDOM, "lambda2 must be less than or equal to 5000");

    /* External units: SI; internal units: solar masses */
    const REAL8 mBH = mBH_SI / LAL_MSUN_SI;
    const REAL8 mNS = mNS_SI / LAL_MSUN_SI;
    const REAL8 M = mBH + mNS;
    const REAL8 M_sec = M * LAL_MTSUN_SI; /* Total mass in seconds */
    REAL8 eta = mBH * mNS / (M * M);        /* Symmetric mass-ratio */
    if (eta > 0.25)
        PhenomInternal_nudge(&eta, 0.25, 1e-6);
    if (eta > 0.25 || eta < 0.0)
        XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");
    XLAL_CHECK(mNS <= 3.0, XLAL_EDOM, "NS mass must be less than or equal to 3 solar masses");
    if(mNS < 1)
        XLALPrintWarning("IMRPhenomNSBH is not calibrated for a NS mass less than 1\n");

    const REAL8 chi_eff = XLALSimIMRPhenomBComputeChi(mBH, mNS, chi_BH, chi_NS);
    LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0, 0}

    UINT4 offset;
    size_t n_full;

    // compute phenomC attributes
    BBHPhenomCParams *params = ComputeIMRPhenomCParams(mBH, mNS, chi_eff, extraParams);
    if (!params)
        XLAL_ERROR(XLAL_EFAULT, "PhenomC properties was returned as a NULL pointer");
    if (params->g1 < 0.0) {
      XLALPrintWarning("IMRPhenomNSBH phenomenological coefficient gamma_1 = %f. gamma_1 has been increased to 0.0 to remove unphysical zeros in the amplitude\n", params->g1);
      params->g1 = 0.0;
    }
    if (params->del1 < 0.0) {
      XLALPrintWarning("IMRPhenomNSBH phenomenological coefficient delta_1 = %f. delta_1 has been increased to 0.0 to remove unphysical zeros in the amplitude\n", params->del1);
      params->del1 = 0.0;
    }
    if (params->del2 < 1.0E-4) {
      XLALPrintWarning("IMRPhenomNSBH phenomenological coefficient delta_2 = %f. delta_2 has been increased to 1e-4 to remove unphysical zeros in the amplitude\n", params->del2);
      params->del2 = 1.0E-4;
    }

    // compute phenomNSBH terms; this only takes the BH spin, not chi_eff
    BBHPhenomNSBHParams *NSBH_params = ComputeIMRPhenomNSBHParams(mBH, mNS, chi_BH, lambda_NS, params);
    if (!NSBH_params)
        XLAL_ERROR(XLAL_EFAULT, "PhenomNSBH properties was returned as a NULL pointer");

    // Prepare frequeny series
    if (f_max == 0) {
      f_max = 0.2 / M_sec;
    }

    if (deltaF > 0){ // freqs contains uniform frequency grid with spacing deltaF; we start at frequency 0

        n_full = PhenomInternal_NextPow2(f_max / deltaF) + 1;

        XLAL_CHECK(XLALGPSAdd(&ligotimegps_zero, -1.0 / deltaF), XLAL_EFUNC, "Failed to shift coalescence time to t=0, tried to apply shift of -1.0/deltaF with deltaF=%g.", deltaF);
        *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &ligotimegps_zero, 0.0, deltaF, &lalStrainUnit, n_full);
        XLAL_CHECK(*htilde, XLAL_ENOMEM, "Failed to allocate waveform COMPLEX16FrequencySeries of length %zu for f_max=%f, deltaF=%g.", n_full, f_max, deltaF);

        // Recreate freqs using only the lower and upper bounds
        UINT4 iStart = (UINT4)(f_min / deltaF);
        UINT4 iStop = (UINT4)(f_max / deltaF);
        freqs = XLALCreateREAL8Sequence(iStop - iStart);
        if (!freqs)
            XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
        for (UINT4 i = iStart; i < iStop; i++)
            freqs->data[i - iStart] = i * deltaF;

        offset = iStart;

    } else { // freqs contains frequencies with non-uniform spacing; we start at lowest given frequency
        n_full = freqs_in->length;

        *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &ligotimegps_zero, f_min, deltaF, &lalStrainUnit, n_full);
        XLAL_CHECK(*htilde, XLAL_ENOMEM, "Failed to allocated waveform COMPLEX16FrequencySeries of length %zu from sequence.", n_full);

        freqs = XLALCreateREAL8Sequence(freqs_in->length);
        if (!freqs)
            XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
        for (UINT4 i = 0; i < n_full; i++)
            freqs->data[i] = freqs_in->data[i]; // just copy input
        offset = 0;
    }

    //n_full is now the length for the FrequencySeries
    memset((*htilde)->data->data, 0, n_full * sizeof(COMPLEX16));
    XLALUnitMultiply(&((*htilde)->sampleUnits), &((*htilde)->sampleUnits), &lalSecondUnit);

    // compute phenomD phase
    int errcode = init_useful_powers(&powers_of_pi, LAL_PI);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_useful_powers() failed.");

    // IMRPhenomD assumes that m1 >= m2.
    pPhi = XLALMalloc(sizeof(IMRPhenomDPhaseCoefficients));

    // final spin
    REAL8 finspin = NSBH_params->chif;

    ComputeIMRPhenomDPhaseCoefficients(pPhi, eta, chi_BH, chi_NS, finspin, extraParams);
    if (extraParams == NULL)
    {
        extraParams = XLALCreateDict();
    }
    XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams, LAL_SIM_INSPIRAL_SPIN_ORDER_35PN);
    XLALSimInspiralTaylorF2AlignedPhasing(&pn, mBH, mNS, chi_BH, chi_NS, extraParams);

    // Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
    // (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
    // was not available when PhenomD was tuned.
    pn->v[6] -= (Subtract3PNSS(mBH, mNS, M, eta, chi_BH, chi_NS) * pn->v[0]);

    PhiInsPrefactors phi_prefactors;
    errcode = init_phi_ins_prefactors(&phi_prefactors, pPhi, pn);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_phi_ins_prefactors failed");

    ComputeIMRPhenDPhaseConnectionCoefficients(pPhi, pn, &phi_prefactors, 1.0, 1.0);

    // incorporating fRef
    const REAL8 MfRef = M_sec * fRef;
    UsefulPowers powers_of_fRef;
    int status = init_useful_powers(&powers_of_fRef, MfRef);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers failed for MfRef");
    const REAL8 phi0 = IMRPhenDPhase(MfRef, pPhi, pn, &powers_of_fRef, &phi_prefactors, 1.0, 1.0);

    // factor of 2 b/c phi0 is orbital phase
    const REAL8 phi_precalc = 2. * phiRef + phi0;

    // use fmax frequency to compute t0
    REAL8 t0 = DPhiMRD(f_max * M_sec, pPhi, 1.0, 1.0);

    // get tidal phase
    REAL8Sequence *phi_tidal = XLALCreateREAL8Sequence(n_full);
    // amp_tidal and planck_taper are not used but this function calculates it anyway.
    REAL8Sequence *amp_tidal = XLALCreateREAL8Sequence(n_full);
    REAL8Sequence *planck_taper = XLALCreateREAL8Sequence(n_full);
    status = XLALSimNRTunedTidesFDTidalPhaseFrequencySeries(
        phi_tidal, amp_tidal, planck_taper, freqs,
        mBH_SI, mNS_SI, 0.0, lambda_NS, NRTidalv2_V);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "XLALSimNRTunedTidesFDTidalPhaseFrequencySeries Failed.");

    const REAL8 amp0 = XLALSimPhenomUtilsFDamp0(M, distance);
    const REAL8 ylm_fac = 2. * sqrt(5. / (64. * LAL_PI));
    // loop over frequencies
    int status_in_for = XLAL_SUCCESS;
    for (UINT4 i = 0; i < freqs->length; i++)
    { // loop over frequency points in sequence
        double fHz = freqs->data[i];
        double Mf = M_sec * fHz;
        int j = i + offset; // shift index for frequency series if needed

        UsefulPowers powers_of_f;
        status_in_for = init_useful_powers(&powers_of_f, Mf);
        if (XLAL_SUCCESS != status_in_for)
        {
            XLALPrintError("init_useful_powers failed for Mf, status_in_for=%d", status_in_for);
            status = status_in_for;
        }
        else{

            REAL8 amp = PhenomNSBHAmplitudeOneFrequency(fHz, params, NSBH_params);

            REAL8 phi = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_f, &phi_prefactors, 1.0, 1.0);
            phi += phi_tidal->data[i] - t0 * (Mf - MfRef) - phi_precalc;
            ((*htilde)->data->data)[j] = amp0 * amp * cexp(-I * phi) * ylm_fac;
        }
    }

    // clean up and return

    /* If extraParams was allocated in this function and not passed in
   * we need to free it to prevent a leak */
    if (extraParams && !extraParams_in)
    {
        XLALDestroyDict(extraParams);
    }
    else
    {
        XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams, LAL_SIM_INSPIRAL_SPIN_ORDER_ALL);
    }

    if (params)
        XLALFree(params);
    if (NSBH_params)
        XLALFree(NSBH_params);
    if (pPhi)
        XLALFree(pPhi);
    if (pn)
        XLALFree(pn);

    if (freqs)
        XLALDestroyREAL8Sequence(freqs);

    if (phi_tidal)
        XLALDestroyREAL8Sequence(phi_tidal);
    if (amp_tidal)
        XLALDestroyREAL8Sequence(amp_tidal);
    if (planck_taper)
        XLALDestroyREAL8Sequence(planck_taper);

    return XLAL_SUCCESS;
}


/**
 * Computes PhenomNSBH Amplitude at a single frequency
 */
static REAL8 PhenomNSBHAmplitudeOneFrequency(
    const REAL8 f,                         /**< frequency Hz */
    const BBHPhenomCParams *params,        /**< pointer to Object storing coefficients and constants: PhenomC */
    const BBHPhenomNSBHParams *params_NSBH /**< pointer to Object storing coefficients and constants: PhenomC_NSBH */
)
{

    REAL8 v = cbrt(params->piM * f);
    REAL8 Mf = params->m_sec * f;

    REAL8 v2 = v * v;
    REAL8 v3 = v * v2;
    REAL8 v4 = v2 * v2;
    REAL8 v5 = v3 * v2;
    REAL8 v6 = v3 * v3;
    REAL8 v7 = v4 * v3;
    REAL8 v10 = v5 * v5;

    /* Get the amplitude */
    REAL8 xdot = 1. + params->xdota2 * v2 + params->xdota3 * v3 + params->xdota4 * v4 +
                 params->xdota5 * v5 + (params->xdota6 + params->xdota6log * log(v2)) * v6 +
                 params->xdota7 * v7;
    xdot *= (params->xdotaN * v10);

    if (xdot < 0.0 && f < params->f1)
    {
        XLALPrintError("omegaDot < 0, while frequency is below SPA-PM matching freq.");
        XLAL_ERROR(XLAL_EDOM);
    }

    REAL8 aPN = 0.0;
    REAL8 aPM = 0.0;

    /* Following Emma's code, take only the absolute value of omegaDot, when
   * computing the amplitude */
    REAL8 omgdot = 1.5 * v * xdot;
    REAL8 ampfac = sqrt(fabs(LAL_PI / omgdot));

    /* Get the real and imaginary part of the PM amplitude */
    REAL8 AmpPNre = ampfac * params->AN * v2 * (1. + params->A2 * v2 + params->A3 * v3 + params->A4 * v4 + params->A5 * v5 + (params->A6 + params->A6log * log(v2)) * v6);
    REAL8 AmpPNim = ampfac * params->AN * v2 * (params->A5imag * v5 + params->A6imag * v6);

    /* Following Emma's code, we take the absolute part of the complex SPA
   * amplitude, and keep that as the amplitude */
    aPN = sqrt(AmpPNre * AmpPNre + AmpPNim * AmpPNim);

    aPM = (params_NSBH->gamma_correction * params->g1 * pow(Mf, (5. / 6.))); //gamma_correction is 1 in BBH case

    /* The Ring-down aamplitude */
    REAL8 Mfrd = 0.0;
    REAL8 sig2 = 0.0;

    Mfrd = params_NSBH->f_RD;
    sig2 = params_NSBH->sigma * params_NSBH->sigma;

    REAL8 L = sig2 / ((Mf - Mfrd) * (Mf - Mfrd) + sig2 / 4.);

    REAL8 aRD = params_NSBH->epsilon_tide * params->del1 * L * pow(Mf, (-7. / 6.)); //epsilon_tide is 1 in BBH case

    /* final amplitude contributions */
    REAL8 wMinusf0_PN = 0.0;
    REAL8 wMinusf0_PM = 0.0;
    REAL8 wPlusf0 = 0.0;

    wMinusf0_PN = wMinus(f, params_NSBH->epsilon_ins * params_NSBH->f0_tilde_PN, params->d0 + params_NSBH->sigma_tide, params);
    wMinusf0_PM = wMinus(f, params_NSBH->f0_tilde_PM, params->d0 + params_NSBH->sigma_tide, params);
    wPlusf0 = wPlus(f, params_NSBH->f0_tilde_RD, params->d0 + params_NSBH->sigma_tide, params);

    REAL8 amp = -(aPN * wMinusf0_PN + aPM * wMinusf0_PM + aRD * wPlusf0);

    return amp;
}

/**
 * PhenomC parameters for NSBH amplitude model
 */
static BBHPhenomNSBHParams *ComputeIMRPhenomNSBHParams(
    const REAL8 m1,                     /**< Mass of companion 1 (solar masses) */
    const REAL8 m2,                     /**< Mass of companion 2 (solar masses) */
    const REAL8 chi,                    /**< Dimensionless spin of black hole */
    const REAL8 lambda,                 /**< Dimensionless tidal deformability of NS */
    const BBHPhenomCParams *params      /**< pointer to Object storing coefficients and constants: PhenomC */
)
{

    BBHPhenomNSBHParams *p = (BBHPhenomNSBHParams *)XLALMalloc(sizeof(BBHPhenomNSBHParams));
    if (!p)
        XLAL_ERROR_NULL(XLAL_EFUNC);
    memset(p, 0, sizeof(BBHPhenomNSBHParams));

    /**
     * Compute NSBH parameters
     */

    /**
     * mixed_merger
     */

    const REAL8 q = m1 / m2;
    const REAL8 MBH = m1;
    const REAL8 MNS = m2;
    p->m_sec = (MBH+MNS) * LAL_MTSUN_SI;
    p->lambda = lambda;

    // Get NS compactness and baryonic mass
    p->C = XLALSimNSBH_compactness_from_lambda(lambda);

    // Get Torus mass/NS baryonic mass
    p->Mtorus = XLALSimNSBH_torus_mass_fit(q, chi, p->C);

    // Get remnant spin for assumed aligned spin system
    p->chif = XLALBHNS_spin_aligned(MBH, MNS, chi, lambda);

    // Get remnant mass scaled to a total (initial) mass of 1
    const REAL8 Mtot = MBH + MNS;
    p->final_mass = XLALBHNS_mass_aligned(MBH, MNS, chi, lambda) / Mtot;

    const COMPLEX16 omega_tilde = XLALSimIMRPhenomNSBH_omega_tilde(p->chif);

    // Prepare remnant dependant quantities
    p->f_RD = creal(omega_tilde)/2.0/LAL_PI/p->final_mass;
    const REAL8 mu = q * p->C;
    const REAL8 xi = XLALSimNSBH_xi_tide(q, chi, mu);

    // In MBH units, 1-2C follows arXiv:1207.6304v1 on the torus mass
    const REAL8 rtide = xi * (1.0 - 2.0 * p->C) / mu;

    p->q_factor = creal(omega_tilde)/cimag(omega_tilde)/2.0;
    p->f_tide = fabs( XLALSimNSBH_fGWinKerr(rtide, 1.0, chi) * (1.0 + 1.0 / q) );

    /**
   * amp
   */

    const REAL8 f_RD_tilde = 0.99 * 0.98 * p->f_RD;

    // NSBH systems seem to require a correction to gamma fitted to BBH cases
    const REAL8 lambda2 = lambda*lambda;
    if (lambda > 1.0)
        p->gamma_correction = 1.25;
    else
        p->gamma_correction = 1.0 + 0.5*lambda - 0.25*lambda2;

    // Ringdown damping w/ fudge factor
    if (lambda > 1.0) {
        p->delta_2_prime = XLALSimIMRPhenomNSBH_delta2_prime(p->f_tide, f_RD_tilde);
    } else {
        const REAL8 c_2 = params->del2 - 0.81248;
        p->delta_2_prime = params->del2 - 2*c_2*lambda + c_2*lambda2;
    }
    p->sigma = p->delta_2_prime * p->f_RD / p->q_factor;

    // Determine the type of merger we see
    if (p->f_tide < p->f_RD)
    {
        // mildly disruptive or totally disruptive, f_tide < f_RD
        p->epsilon_tide = 0.0;

        p->epsilon_ins = XLALSimIMRPhenomNSBH_epsilon_ins_with_torus_mass(p->Mtorus, p->C, q, chi);
        if (p->epsilon_ins > 1.)
            p->epsilon_ins = 1.;

        p->sigma_tide = XLALSimIMRPhenomNSBH_sigma_tide_with_torus_mass(p->Mtorus, p->C, q, chi);
        if (p->Mtorus > 0.0)
        {
            // totally disruptive merger, Mtorus > 0
            p->merger_type = DISRUPTIVE;
            p->f0_tilde_PN = p->f_tide / p->m_sec;
            p->f0_tilde_PM = p->f_tide / p->m_sec;
            p->f0_tilde_RD = 0.0;
        }
        else
        {
            // mildly disruptive with no remnant, Mtorus == 0
            p->merger_type = MILDLY_DISRUPTIVE_NO_TORUS_REMNANT;

            const REAL8 f1 = (1.0 - 1.0 / q) * f_RD_tilde + p->epsilon_ins * p->f_tide / q;
            const REAL8 f2 = (1.0 - 1.0 / q) * f_RD_tilde + p->f_tide / q;

            p->f0_tilde_PN = f1 / p->m_sec;
            p->f0_tilde_PM = f2 / p->m_sec;
            p->f0_tilde_RD = 0.0;

            // for this case, sigma_tide is computed by averaging the
            // disruptive and non-disruptive values
            const REAL8 sigma_tide_ND = XLALSimIMRPhenomNSBH_sigma_tide_ND(
                XLALSimIMRPhenomNSBH_x_ND_prime(p->f_tide, f_RD_tilde, p->C, chi));
            p->sigma_tide = (p->sigma_tide + sigma_tide_ND) / 2.0;
        }
    }
    else
    {
        // mildly disruptive or non-disruptive, f_tide >= f_RD
        if (lambda > 1.0) {
            p->f0_tilde_PN = f_RD_tilde / p->m_sec;
            p->f0_tilde_PM = f_RD_tilde / p->m_sec;
            p->f0_tilde_RD = f_RD_tilde / p->m_sec;
        } else {
            const REAL8 f0 = (1.0 - 0.02*lambda + 0.01*lambda2)*0.98*p->f_RD;
            p->f0_tilde_PN = f0 / p->m_sec;
            p->f0_tilde_PM = f0 / p->m_sec;
            p->f0_tilde_RD = f0 / p->m_sec;
        }

        p->epsilon_tide = XLALSimIMRPhenomNSBH_epsilon_tide_ND(
            XLALSimIMRPhenomNSBH_x_ND(p->f_tide, f_RD_tilde, p->C, chi));
        p->sigma_tide = XLALSimIMRPhenomNSBH_sigma_tide_ND(
            XLALSimIMRPhenomNSBH_x_ND_prime(p->f_tide, f_RD_tilde, p->C, chi));

        if (p->Mtorus > 0.0)
        {
            // mildly disruptive with remnant mass, f_tide > f_RD AND Mtorus > 0
            p->merger_type = MILDLY_DISRUPTIVE_TORUS_REMNANT;
            p->epsilon_ins = XLALSimIMRPhenomNSBH_epsilon_ins_with_torus_mass(p->Mtorus, p->C, q, chi);
        }
        else
        {
            // non-disruptive, f_tide > f_RD AND Mtorus == 0
            p->merger_type = NON_DISRUPTIVE;
            p->epsilon_ins = 1.0;
        }
    }

    return p;
}

/**
 * @addtogroup LALSimIMRPhenom_c
 *
 * @{
 *
 * @name Routines for IMR Phenomenological Model "NSBH"
 *
 * @{
 *
 * @author Jonathan Thompson, Edward Fauchon-Jones, Sebastian Khan
 *
 * @brief C code for <code>IMRPhenomNSBH</code> phenomenological waveform model.
 *
 * This is a single-spin, non-precessing frequency domain model. This model is
 * based on the amplitude model described by \cite Pannarale:2015jka and the
 * <code>IMRPhenomD</code> based NRTidal phase model described by
 * \cite Dietrich:2019kaq. Please see
 * <a href="https://dcc.ligo.org/LIGO-T1900729">LIGO-T1900729</a> for a
 * technical description of the implemented model.
 *
 * @note The model can be evaluated within the following parameter space
 * boundary outside of which an <code>XLAL_EDOM</code> error will be thrown
 * * \f$ m_{\mathrm{NS}} \leq 3 M_{\odot} \f$
 * * \f$ 1 \leq q \leq 100 \f$
 * * \f$ 0 \leq \Lambda_2 \leq 5000 \f$
 *
 * @note The model will throw a warning if it is evaluated inside the above
 * parameter space boundary but violates any of the following conditions
 * * \f$ \chi_{\mathrm{NS}} = 0 \f$
 * * \f$ m_{\mathrm{NS}} \geq 1 M_{\odot} \f$
 * * \f$ \delta_1 \geq 0 \f$
 * * \f$ \delta_2 \geq 10^{-4} \f$
 * * \f$ \gamma_1 \geq 0 \f$
 *
 * @note If any of the conditions on the phenomenological coefficient \f$
 * \delta_1, \delta_2, \gamma_1 \f$ are violated then they are increased to the
 * values \f$ 0, 10^{-4}, 0 \f$ respectively to remove unphysical zeros in the
 * amplitude.
 *
 * @note The models amplitude was calibrated to mass-ratios [1:2,1:3,1:4,1:5].
 * * Along the mass-ratio 1:2 line it was calibrated to BH spins [-0.5, 0.75].
 * * Along the mass-ratio 1:3 line it was calibrated to BH spins [-0.5, 0.75].
 * * Along the mass-ratio 1:4 line it was calibrated to BH spins [0, 0.75].
 * * Along the mass-ratio 1:5 line it was calibrated to BH spins [0, 0.75].
 *
 * @note Please see \cite Yamamoto:2008js, \cite Lackey:2013axa and
 * \cite Pannarale:2015jka for full details of the NR data used to calibrate
 * the amplitude for this model.
 *
 * @note The models phase uses the phase of <code>IMRPhenomD_NRTidalv2</code>. For
 * full details of the NR data used to calibrate the phase for this model please
 * see \cite Husa:2015iqa, \cite Khan:2015jqa and \cite Dietrich:2019kaq
 *
 * @attention The model is usable outside this parameter range,
 * and tests have shown that the model produces sensible results. However the
 * current set of numerical relativity simulations for NSBH systems is limited.
 * In particular they do not cover the mass ratio ranges and spin ranges of
 * numerical relativity simulations that are available for BBH systems. As such
 * you should carefully consider applications of this model for use case when
 * evaluated outside the suggested parameter space. For more information please
 * see the review wiki which can be found at
 * https://git.ligo.org/waveforms/reviews/nsbh-models/wikis/home.
 *
 */

/**
 * Compute waveform in LAL format at specified frequencies for the IMRPhenomNSBH model.
 *
 * XLALSimIMRPhenomNSBH() returns the plus and cross polarizations as a complex
 * frequency series with equal spacing deltaF and contains zeros from zero frequency
 * to the starting frequency and zeros beyond the cutoff frequency in the ringdown.
 *
 * In contrast, XLALSimIMRPhenomNSBHFrequencySequence() returns a
 * complex frequency series with entries exactly at the frequencies specified in
 * the sequence freqs (which can be unequally spaced). No zeros are added.
 *
 * This function is designed as an entry point for reduced order quadratures.
 */
int XLALSimIMRPhenomNSBHFrequencySequence(
    COMPLEX16FrequencySeries **htilde, /**< Output: Frequency-domain waveform h+ */
    const REAL8Sequence *freqs,        /**< Frequency points at which to evaluate the waveform (Hz) */
    REAL8 phiRef,                      /**< Phase at reference time */
    REAL8 fRef,                        /**< Reference frequency (Hz); 0 defaults to fLow */
    REAL8 distance,                    /**< Distance of source (m) */
    REAL8 mBH_SI,                      /**< Mass of BH (kg) */
    REAL8 mNS_SI,                      /**< Mass of neutron star 2 (kg) */
    REAL8 chi_BH,                      /**< Dimensionless aligned component spin of Black Hole */
    REAL8 chi_NS,                      /**< Dimensionless aligned component spin of NS */
    LALDict *extraParams               /**< linked list containing the extra testing GR parameters and tidal parameters */
)
{
    if (!freqs)
        XLAL_ERROR(XLAL_EFAULT, "freqs must not be a NULL pointer\n");

    // Call the internal core function with deltaF = 0 to indicate that freqs is non-uniformly
    // spaced and we want the strain only at these frequencies
    int retcode = IMRPhenomNSBH_Core(htilde,
                                          phiRef, fRef, distance, mBH_SI, mNS_SI, chi_BH, chi_NS, extraParams, freqs, 0);

    return (retcode);
}

/**
 * Driver routine to compute the single-spin, non-precessing,
 * neutron-star-black-hole, inspiral-merger-ringdown phenomenological waveform
 * IMRPhenomNSBH in the frequency domain in LAL format.
 *
 * All input parameters should be in SI units. Angles should be in radians.
 *
 * Returns the plus and cross polarizations as a complex frequency series with
 * equal spacing deltaF and contains zeros from zero frequency to the starting
 * frequency fLow and zeros beyond the cutoff frequency in the ringdown.
 */
int XLALSimIMRPhenomNSBH(
    COMPLEX16FrequencySeries **htilde, /**< Output: Frequency-domain waveform h+ */
    REAL8 phiRef,                      /**< Phase at reference time */
    REAL8 deltaF,                      /**< Sampling frequency (Hz) */
    REAL8 fLow,                        /**< Starting GW frequency (Hz) */
    REAL8 fHigh,                       /**< End frequency; 0 defaults to Mf=0.2 */
    REAL8 fRef,                        /**< Reference frequency (Hz); 0 defaults to fLow */
    REAL8 distance,                    /**< Distance of source (m) */
    REAL8 mBH_SI,                       /**< Mass of BH (kg) */
    REAL8 mNS_SI,                       /**< Mass of neutron star 2 (kg) */
    REAL8 chi_BH,                        /**< Dimensionless aligned component spin of Black Hole */
    REAL8 chi_NS,                        /**< Dimensionless aligned component spin of NS */
    LALDict *extraParams               /**< linked list containing the extra testing GR parameters and tidal parameters */
)
{
    // Use fLow, fHigh, deltaF to compute freqs sequence
    // Instead of building a full sequence we only transfer the boundaries and let
    // the internal core function do the rest (and properly take care of corner cases).
    REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
    freqs->data[0] = fLow;
    freqs->data[1] = fHigh;

    int retcode = IMRPhenomNSBH_Core(htilde,
                                          phiRef, fRef, distance, mBH_SI, mNS_SI, chi_BH, chi_NS, extraParams, freqs, deltaF);

    XLALDestroyREAL8Sequence(freqs);

    return (retcode);
}


/***********************************************/
/* Utility functions for PhenomNSBH parameters */
/***********************************************/

/**
 * @addtogroup LALSimIMRPhenomNSBHUtility XLALSimIMRPhenomNSBHUtility
 * @brief C code for utility routines for <code>IMRPhenomNSBH</code>
 * phenomenological waveform model.
 *
 * @{
 *
 * @name Utility routines for IMR Phenomenological Model "NSBH"
 *
 * @{
 *
 * @author Jonathan Thompson, Edward Fauchon-Jones, Sebastian Khan
 */

/**
 * Convenience function for expression appearing in disruptive merger
 *
 * Combination of parameters that helps understand and handle disruptive cases
 * which produce tori. This is Eq. (23) in arXiv:1509.00512 on pg. (7).
 */
double XLALSimIMRPhenomNSBH_x_D(
    const REAL8 Mtorus,  /**< Baryonic mass of the torus remnant of a BH-NS merger in units of the NS baryonic mass. */
    const REAL8 C,       /**< Neutron star compactness */
    const REAL8 q,       /**< Mass ratio of the NSBH system M_BH/M_NS > 1 */
    const REAL8 chi      /**< Dimensionless spin parameter -1 <= chi <= 1 */
) {
  REAL8 nu = XLALSimIMRPhenomNSBH_eta_from_q(q);
  return Mtorus + 0.424912*C + 0.363604*sqrt(nu) - 0.0605591*chi;
}

/**
 * Correction to the inspiral transition frequency with spin contributions
 *
 * Correction factor that multiplies the inspiral transition frequency when
 * the mixed phenom GW is generated. See https://arxiv.org/abs/1509.00512 Eq.
 * (23) (p7) for it's definition.
 */
double XLALSimIMRPhenomNSBH_epsilon_ins_with_torus_mass(
    const REAL8 Mtorus,  /**< Baryonic mass of the torus remnant of a BH-NS merger in units of the NS baryonic mass. */
    const REAL8 C,       /**< Neutron star compactness */
    const REAL8 q,       /**< Mass ratio of the NSBH system M_BH/M_NS > 1 */
    const REAL8 chi      /**< Dimensionless spin parameter -1 <= chi <= 1 */
) {
  return 1.29971 - 1.61724 * XLALSimIMRPhenomNSBH_x_D(Mtorus, C, q, chi);
}

/**
 * Convinience function for expression appearing in disruptive merger
 *
 * Combination of parameters that helps understand and handle disruptive cases
 * which produce tori.
 *
 * See Eq. (25) of arXiv:1509:00512
 */
double XLALSimIMRPhenomNSBH_x_D_prime(
    const REAL8 Mtorus,  /**< Baryonic mass of the torus remnant of a BH-NS merger in units of the NS baryonic mass. */
    const REAL8 C,       /**< Neutron star compactness */
    const REAL8 q,       /**< Mass ratio of the NSBH system M_BH/M_NS > 1 */
    const REAL8 chi      /**< Dimensionless spin parameter -1 <= chi <= 1 */
) {
  REAL8 nu = XLALSimIMRPhenomNSBH_eta_from_q(q);
  return Mtorus - 0.132754*C + 0.576669*sqrt(nu) - 0.0603749*chi - 0.0601185*pow(chi, 2.0) - 0.0729134*pow(chi, 3.0);
}

/**
 * Correction to ringdown Lorentzian width for disruptive mergers
 *
 * Correction to the Lorentzian width parameter in the ringdown ansatz used when
 * the mixed phenom GW is generated for disruptive mergers that do produce tori.
 * See https://arxiv.org/abs/1509.00512 Eq. (24) (p7) for it's definition.
 */
double XLALSimIMRPhenomNSBH_sigma_tide_with_torus_mass(
    const REAL8 Mtorus,  /**< Baryonic mass of the torus remnant of a BH-NS merger in units of the NS baryonic mass. */
    const REAL8 C,       /**< Neutron star compactness */
    const REAL8 q,       /**< Mass ratio of the NSBH system M_BH/M_NS > 1 */
    const REAL8 chi      /**< Dimensionless spin parameter -1 <= chi <= 1 */
) {
  return 0.137722 - 0.293237*XLALSimIMRPhenomNSBH_x_D_prime(Mtorus, C, q, chi);
}

/**
 * PhenomC parameter delta_1 NSBH correction factor
 *
 * Factor that mutltiplies the (PhenomC) parameter delta_1 when the mixed phenom
 * GW is generated. See https://arxiv.org/abs/1509.00512 Eq. (16) (p5) for its
 * definition.
 */
double XLALSimIMRPhenomNSBH_epsilon_tide_ND(
    const REAL8 x_ND /**< Dimensionless fit parameters x_ND. See https://arxiv.org/abs/1506.00512 Eq. (17) (p6) for it's original definition (x_ND). */
) {
  return 2.0*XLALSimIMRPhenomNSBH_window_plus(x_ND, -0.0796251, 0.0801192);
}

/**
 * Correction to ringdown Lorentzian width for nondisruptive mergers
 *
 * Correction to the Lorentzian width parameter in the ringdown ansatz used when
 * the mixed phenom GW is generated for nondisruptive mergers that do not
 * produce tori. See https://arxiv.org/abs/1509.00512 Eq. (18) (p6) for its
 * definition.
 */
double XLALSimIMRPhenomNSBH_sigma_tide_ND(
    const REAL8 x_ND_prime /**< Dimensionless fit parameter. See https://arxiv.org/abs/1509.00512 Eq. (19) (p6) for it's original definition (x_ND'). */
) {
  return 2.0*XLALSimIMRPhenomNSBH_window_minus(x_ND_prime, -0.206465, 0.226844);
}

/**
 * Convinience function for expression appearing in disruptive merger
 *
 * Combination of parameters that helps understand and handle non-disruptive
 * cases. See Eq. (17) of arXiv:1509:00512
 */
double XLALSimIMRPhenomNSBH_x_ND(
    const REAL8 f_tide,        /**< Frequency at which the tidal interactions occur */
    const REAL8 f_RD_tilde,    /**< Scaled ringdown frequency */
    const REAL8 C,             /**< Neutron star compactness */
    const REAL8 chi            /**< Dimensionless spin parameter -1 <= chi <= 1 */
) {
  return pow((f_tide - f_RD_tilde)/f_RD_tilde, 2.0) - 0.571505*C - 0.00508451*chi;
}

/**
 * Convinience function for expression appearing in disruptive merger
 *
 * Combination of parameters that helps understand and handle non-disruptive
 * cases. See Eq. (19) of arXiv:1509:00512
 */
double XLALSimIMRPhenomNSBH_x_ND_prime(
    const REAL8 f_tide,        /**< Frequency at which the tidal interactions occur */
    const REAL8 f_RD_tilde,    /**< Scaled ringdown frequency */
    const REAL8 C,             /**< Neutron star compactness */
    const REAL8 chi            /**< Dimensionless spin parameter -1 <= chi <= 1 */
) {
  return pow((f_tide - f_RD_tilde)/f_RD_tilde, 2.0) - 0.657424*C - 0.0259977*chi;
}

/**
 * Fitted coefficient for PhenomC Lorentzian
 *
 * See Eq. (20) of arXiv:1509:00512
 */
double XLALSimIMRPhenomNSBH_delta2_prime(
    const REAL8 f_tide,        /**< Frequency at which the tidal interactions occur */
    const REAL8 f_RD_tilde     /**< Scaled ringdown frequency */
) {
  return 1.62496*XLALSimIMRPhenomNSBH_window_plus((f_tide - f_RD_tilde)/f_RD_tilde, 0.0188092, 0.338737);
}

/**
 * Hyperbolic tangent sigmoid function
 *
 * This function approaches 0.5 as f approaches infinity. It approaches 0 as
 * f approaches -infinity. It's value will be 0.25 at f0. As d approaches 0
 * it approaches a step function.
 */
double XLALSimIMRPhenomNSBH_window_plus(
  const REAL8 f,   /**< Value at which to evaluate sigmoid function */
  const REAL8 f0,  /**< Center of sigmoid function */
  const REAL8 d    /**< Width of sigmoid function */
) {
  return 0.25*(1 + tanh(4.0*(f-f0)/d));
}

/**
 * Hyperbolic tangent sigmoid function
 *
 * This function approaches 0.5 as f approaches -infinity. It approaches 0 as
 * f approaches infinity. It's value will be 0.25 at f0. As d approaches 0
 * it approaches a step function.
 */
double XLALSimIMRPhenomNSBH_window_minus(
  const REAL8 f,   /**< Value at which to evaluate sigmoid function */
  const REAL8 f0,  /**< Center of sigmoid function */
  const REAL8 d    /**< Width of sigmoid function */
) {
  return 0.25*(1 - tanh(4.0*(f-f0)/d));
}

/**
 * Convenience function to calculate symmetric mass ratio from q
 */
double XLALSimIMRPhenomNSBH_eta_from_q(
    const REAL8 q   /**< Mass ratio */
) {
  return q/pow(1.0+q, 2.0);
}

/**
 * NS baryonic mass as a function of NS gravitational mass
 *
 * NS baryonic mass as a function of NS gravitational mass. Both are in solar
 * mass units. See Eq. (22) in https://arxiv.org/abs/1601.06083 for the
 * definition of this expression.
 */
double XLALSimIMRPhenomNSBH_baryonic_mass_from_C(
  const REAL8 C,   /**< Compactness of neutron star. */
  const REAL8 Mg   /**< Neutron star gravitational mass. */
) {
  REAL8 d1 = 6.19E-1;
  REAL8 d2 = 1.359E-1;
  return Mg + Mg*(d1*C + d2*pow(C,2.0));
}

/**
 * 220 quasi-normal mode dimensionless frequency
 *
 * See Eq. (22) (p8) in https://arxiv.org/abs/1810.03550 for the definition of
 * this quantity.
 */
COMPLEX16 XLALSimIMRPhenomNSBH_omega_tilde(
  const REAL8 a  /**< Dimensionless final spin parameter of the BH remnant */
) {
  REAL8 kappa = pow(log(2-a)/log(3), 0.5);
  COMPLEX16 omega_tilde = (1.0 + kappa*(
    1.5578*cexp(I*2.9031)
    + 1.9510*cexp(I*5.9210)*kappa
    + 2.0997*cexp(I*2.7606)*pow(kappa,2)
    + 1.4109*cexp(I*5.9143)*pow(kappa,3)
    + 0.4106*cexp(I*2.7952)*pow(kappa,4)));
  return omega_tilde;
}

/* Write a wrapper function to return the NSBH parameters */
int XLALSimIMRPhenomNSBHProperties(
    REAL8 *f_RD,                        /**< Output: NSBH ringdown frequency [Hz] */
    REAL8 *f_tide,                      /**< Output: NSBH tidal disruption frequency [Hz] */
    REAL8 *torus_mass,                  /**< Output: Torus remnant mass (kg) */
    REAL8 *compactness,                 /**< Output: Compactness of neutron star */
    REAL8 *final_mass,                  /**< Output: final mass after merger (kg) */
    REAL8 *chif,                        /**< Output: final dimensionless spin */
    REAL8 mBH_SI,                       /**< Mass of BH (kg) */
    REAL8 mNS_SI,                       /**< Mass of neutron star 2 (kg) */
    REAL8 chi_BH,                        /**< Dimensionless aligned component spin of Black Hole */
    REAL8 lambda_NS                     /**< Dimensionless tidal deformability of NS */
)
{

  const double mBH = mBH_SI / LAL_MSUN_SI;
  const double mNS = mNS_SI / LAL_MSUN_SI;
  const double M_sec = (mBH + mNS) * LAL_MTSUN_SI;

  const REAL8 chi_eff = XLALSimIMRPhenomBComputeChi(mBH, mNS, chi_BH, 0.0);
  BBHPhenomCParams *params = ComputeIMRPhenomCParams(mBH, mNS, chi_eff, NULL);
  if (!params)
      XLAL_ERROR(XLAL_EFUNC);
  if (params->g1 < 0.0)
      params->g1 = 0.0;
  if (params->del1 < 0.0)
      params->del1 = 0.0;
  if (params->del2 < 0.0)
      params->del2 = 0.0;
  BBHPhenomNSBHParams *NSBH_params = ComputeIMRPhenomNSBHParams(mBH, mNS, chi_BH, lambda_NS, params);
  if (!NSBH_params)
      XLAL_ERROR(XLAL_EFAULT, "PhenomNSBH properties was returned as a NULL pointer");

  const double Mb = XLALSimIMRPhenomNSBH_baryonic_mass_from_C(NSBH_params->C, mNS);

  *f_RD = NSBH_params->f_RD / M_sec;
  *f_tide = NSBH_params->f_tide / M_sec;
  *torus_mass = NSBH_params->Mtorus * Mb * LAL_MSUN_SI;
  *compactness = NSBH_params->C;
  *final_mass = NSBH_params->final_mass * LAL_MSUN_SI * (mBH + mNS);
  *chif = NSBH_params->chif;

  if (params)
      XLALFree(params);
  if (NSBH_params)
      XLALFree(NSBH_params);

  return XLAL_SUCCESS;
}

/** @} */

/** @} */

/** @} */

/** @} */
