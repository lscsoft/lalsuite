/*
 *  Copyright (C) 2019 HyunWon Lee, JeongCho Kim, Chunglee Kim, Marc Favata, K.G. Arun
 *  Assembled from code found in:
 *    - LALInspiralStationaryPhaseApproximation2.c
 *    - LALInspiralChooseModel.c
 *    - LALInspiralSetup.c
 *    - LALSimInspiralTaylorF2ReducedSpin.c
 *    - LALSimInspiralTaylorF2.c
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
/*
#include <stdlib.h>
#include <math.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/LALConstants.h>
#include <lal/Sequence.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiralEOS.h>
#include <lal/LALSimInspiral.h>
#include <lal/Units.h>
#include <lal/XLALError.h>
#include <lal/AVFactories.h>
#include "LALSimInspiralPNCoefficients.c"
*/

#ifndef _OPENMP
#define omp ignore
#endif

/**
 * @addtogroup LALSimInspiralTaylorF2Ecc_c
 * @brief Routines for generating eccentric TaylorF2 waveforms
 * @{
 *
 * @review TaylorF2Ecc is reviewed, with review statement and git has details available at https://git.ligo.org/waveforms/reviews/taylorf2ecc/wikis/home.
 *
 * @name Routines for TaylorF2Ecc Waveforms
 * @sa
 * Section IIIF of Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 *
 * Section IV of Marc, et al paper Phys. Rev. D 93, 124061 (2016), arXiv:1605.00304.
 * review page is https://git.ligo.org/waveforms/reviews/taylorf2ecc/wikis/Eccentric-phase-PN-coefficient-form.
 *
 * @{
 */

/**
 * \author Jeongcho Kim, Chunglee Kim, Hyung Won Lee, Marc Favata, K.G. Arun
 * \file
 *
 * \brief Module to compute the eccentric TaylorF2 inspiral waveform for small eccentricity.
 * Code is based on Section IV of Marc, et al paper Phys. Rev. D 93, 124061 (2016), arXiv:1605.00304.
 * Code review page is https://git.ligo.org/waveforms/reviews/taylorf2ecc/wikis/Eccentric-phase-PN-coefficient-form.
 */

int XLALSimInspiralTaylorF2CoreEcc(
        COMPLEX16FrequencySeries **htilde_out, /**< FD waveform */
	const REAL8Sequence *freqs,            /**< frequency points at which to evaluate the waveform (Hz) */
        const REAL8 phi_ref,                   /**< reference orbital phase (rad) */
        const REAL8 m1_SI,                     /**< mass of companion 1 (kg) */
        const REAL8 m2_SI,                     /**< mass of companion 2 (kg) */
        const REAL8 f_ref,                     /**< Reference GW frequency (Hz) - if 0 reference point is coalescence */
	const REAL8 shft,		       /**< time shift to be applied to frequency-domain phase (sec)*/
        const REAL8 r,                         /**< distance of source (m) */
        const REAL8 eccentricity,              /**< eccentricity effect control < 0 : no eccentricity effect */
        LALDict *p,                       /**< Linked list containing the extra parameters >**/
        PNPhasingSeries *pfaP /**< Phasing coefficients >**/
        )
{

    if (!htilde_out) XLAL_ERROR(XLAL_EFAULT);
    if (!freqs) XLAL_ERROR(XLAL_EFAULT);
    /* external: SI; internal: solar masses */
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 eta = m1 * m2 / (m * m);
    const REAL8 piM = LAL_PI * m_sec;
    REAL8 amp0;
    size_t i;
    COMPLEX16 *data = NULL;
    LIGOTimeGPS tC = {0, 0};
    INT4 iStart = 0;

    COMPLEX16FrequencySeries *htilde = NULL;

    if (*htilde_out) { //case when htilde_out has been allocated in XLALSimInspiralTaylorF2
	    htilde = *htilde_out;
	    iStart = htilde->data->length - freqs->length; //index shift to fill pre-allocated data
	    if(iStart < 0) XLAL_ERROR(XLAL_EFAULT);
    }
    else { //otherwise allocate memory here
	    htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &tC, freqs->data[0], 0., &lalStrainUnit, freqs->length);
	    if (!htilde) XLAL_ERROR(XLAL_EFUNC);
	    XLALUnitMultiply(&htilde->sampleUnits, &htilde->sampleUnits, &lalSecondUnit);
    }

    /* phasing coefficients */
    PNPhasingSeries pfa = *pfaP;

    REAL8 pfaN = 0.; REAL8 pfa1 = 0.;
    REAL8 pfa2 = 0.; REAL8 pfa3 = 0.; REAL8 pfa4 = 0.;
    REAL8 pfa5 = 0.; REAL8 pfl5 = 0.;
    REAL8 pfa6 = 0.; REAL8 pfl6 = 0.;
    REAL8 pfa7 = 0.;

    INT4 phaseO=XLALSimInspiralWaveformParamsLookupPNPhaseOrder(p);
    switch (phaseO)
    {
        case -1:
        case 7:
            pfa7 = pfa.v[7];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 6:
            pfa6 = pfa.v[6];
            pfl6 = pfa.vlogv[6];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 5:
            pfa5 = pfa.v[5];
            pfl5 = pfa.vlogv[5];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 4:
            pfa4 = pfa.v[4];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 3:
            pfa3 = pfa.v[3];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 2:
            pfa2 = pfa.v[2];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 1:
            pfa1 = pfa.v[1];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 0:
            pfaN = pfa.v[0];
            break;
        default:
            XLAL_ERROR(XLAL_ETYPE, "Invalid phase PN order %d", phaseO);
    }

    /* Validate expansion order arguments.
     * This must be done here instead of in the OpenMP parallel loop
     * because when OpenMP parallelization is turned on, early exits
     * from loops (via return or break statements) are not permitted.
     */

    /* Validate amplitude PN order. */
    INT4 amplitudeO=XLALSimInspiralWaveformParamsLookupPNAmplitudeOrder(p);
    switch (amplitudeO)
    {
        case -1:
        case 7:
        case 6:
        case 5:
        case 4:
        case 3:
        case 2:
        case 0:
            break;
        default:
            XLAL_ERROR(XLAL_ETYPE, "Invalid amplitude PN order %d", amplitudeO);
    }

    /* Generate tidal terms separately.
     * Enums specifying tidal order are in LALSimInspiralWaveformFlags.h
     */
    REAL8 pft10 = 0.;
    REAL8 pft12 = 0.;
    REAL8 pft13 = 0.;
    REAL8 pft14 = 0.;
    REAL8 pft15 = 0.;
    switch( XLALSimInspiralWaveformParamsLookupPNTidalOrder(p) )
    {
	case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_75PN:
            pft15 = pfa.v[15];
#if __GNUC__ >= 7
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_7PN:
            pft14 = pfa.v[14];
#if __GNUC__ >= 7
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_65PN:
            pft13 = pfa.v[13];
#if __GNUC__ >= 7
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
	    pft12 = pfa.v[12];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
            pft10 = pfa.v[10];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLAL_ERROR(XLAL_EINVAL, "Invalid tidal PN order %d", XLALSimInspiralWaveformParamsLookupPNTidalOrder(p));
    }

    /* The flux and energy coefficients below are used to compute SPA amplitude corrections */

    /* flux coefficients */
    const REAL8 FTaN = XLALSimInspiralPNFlux_0PNCoeff(eta);
    const REAL8 FTa2 = XLALSimInspiralPNFlux_2PNCoeff(eta);
    const REAL8 FTa3 = XLALSimInspiralPNFlux_3PNCoeff(eta);
    const REAL8 FTa4 = XLALSimInspiralPNFlux_4PNCoeff(eta);
    const REAL8 FTa5 = XLALSimInspiralPNFlux_5PNCoeff(eta);
    const REAL8 FTl6 = XLALSimInspiralPNFlux_6PNLogCoeff(eta);
    const REAL8 FTa6 = XLALSimInspiralPNFlux_6PNCoeff(eta);
    const REAL8 FTa7 = XLALSimInspiralPNFlux_7PNCoeff(eta);

    /* energy coefficients */
    const REAL8 dETaN = 2. * XLALSimInspiralPNEnergy_0PNCoeff(eta);
    const REAL8 dETa1 = 2. * XLALSimInspiralPNEnergy_2PNCoeff(eta);
    const REAL8 dETa2 = 3. * XLALSimInspiralPNEnergy_4PNCoeff(eta);
    const REAL8 dETa3 = 4. * XLALSimInspiralPNEnergy_6PNCoeff(eta);


    /* Perform some initial checks */
    if (m1_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (f_ref < 0) XLAL_ERROR(XLAL_EDOM);
    if (r <= 0) XLAL_ERROR(XLAL_EDOM);

    /* extrinsic parameters */
    amp0 = -4. * m1 * m2 / r * LAL_MRSUN_SI * LAL_MTSUN_SI * sqrt(LAL_PI/12.L);

    data = htilde->data->data;

    REAL8 v_ecc_ref = 0.0;
    REAL8 f_ecc = XLALSimInspiralWaveformParamsLookupEccentricityFreq(p);
    if( eccentricity > 0) {
        v_ecc_ref = cbrt(piM*f_ecc);
    }

    /* Compute the SPA phase at the reference point
     * N.B. f_ref == 0 means we define the reference time/phase at "coalescence"
     * when the frequency approaches infinity. In that case,
     * the integrals Eq. 3.15 of arXiv:0907.0700 vanish when evaluated at
     * f_ref == infinity. If f_ref is finite, we must compute the SPA phase
     * evaluated at f_ref, store it as ref_phasing and subtract it off.
     */
    REAL8 ref_phasing = 0.;
    INT4 ecc_order = XLALSimInspiralWaveformParamsLookupPNEccentricityOrder(p);
    if( f_ref != 0. ) {
        const REAL8 vref = cbrt(piM*f_ref);
        const REAL8 logvref = log(vref);
        const REAL8 v2ref = vref * vref;
        const REAL8 v3ref = vref * v2ref;
        const REAL8 v4ref = vref * v3ref;
        const REAL8 v5ref = vref * v4ref;
        const REAL8 v6ref = vref * v5ref;
        const REAL8 v7ref = vref * v6ref;
        const REAL8 v8ref = vref * v7ref;
        const REAL8 v9ref = vref * v8ref;
        const REAL8 v10ref = vref * v9ref;
        const REAL8 v12ref = v2ref * v10ref;
        const REAL8 v13ref = vref * v12ref;
        const REAL8 v14ref = vref * v13ref;
        const REAL8 v15ref = vref * v14ref;
        ref_phasing += pfa7 * v7ref;
        ref_phasing += (pfa6 + pfl6 * logvref) * v6ref;
        ref_phasing += (pfa5 + pfl5 * logvref) * v5ref;
        ref_phasing += pfa4 * v4ref;
        ref_phasing += pfa3 * v3ref;
        ref_phasing += pfa2 * v2ref;
        ref_phasing += pfa1 * vref;
        ref_phasing += pfaN;

        /* Tidal terms in reference phasing */
        ref_phasing += pft15 * v15ref;
        ref_phasing += pft14 * v14ref;
        ref_phasing += pft13 * v13ref;
        ref_phasing += pft12 * v12ref;
        ref_phasing += pft10 * v10ref;

        /* Eccentricity terms in phasing */
        if( eccentricity > 0 ) {
          ref_phasing += eccentricityPhasing_F2(vref, v_ecc_ref, eccentricity, eta, ecc_order);
        }

        ref_phasing /= v5ref;
    } /* End of if(f_ref != 0) block */

    #pragma omp parallel for
    for (i = 0; i < freqs->length; i++) {
        const REAL8 f = freqs->data[i];
        const REAL8 v = cbrt(piM*f);
        const REAL8 logv = log(v);
        const REAL8 v2 = v * v;
        const REAL8 v3 = v * v2;
        const REAL8 v4 = v * v3;
        const REAL8 v5 = v * v4;
        const REAL8 v6 = v * v5;
        const REAL8 v7 = v * v6;
        const REAL8 v8 = v * v7;
        const REAL8 v9 = v * v8;
        const REAL8 v10 = v * v9;
        const REAL8 v12 = v2 * v10;
        const REAL8 v13 = v * v12;
        const REAL8 v14 = v * v13;
        const REAL8 v15 = v * v14;
        REAL8 phasing = 0.;
        REAL8 dEnergy = 0.;
        REAL8 flux = 0.;
        REAL8 amp;

        phasing += pfa7 * v7;
        phasing += (pfa6 + pfl6 * logv) * v6;
        phasing += (pfa5 + pfl5 * logv) * v5;
        phasing += pfa4 * v4;
        phasing += pfa3 * v3;
        phasing += pfa2 * v2;
        phasing += pfa1 * v;
        phasing += pfaN;

        /* Tidal terms in phasing */
        phasing += pft15 * v15;
        phasing += pft14 * v14;
        phasing += pft13 * v13;
        phasing += pft12 * v12;
        phasing += pft10 * v10;

        /* Eccentricity terms in phasing */
        if( eccentricity > 0 ) {
          phasing += eccentricityPhasing_F2(v, v_ecc_ref, eccentricity, eta, ecc_order);
        }
        phasing /= v5;

    /* WARNING! Amplitude orders beyond 0 have NOT been reviewed!
     * Use at your own risk. The default is to turn them off.
     * These do not currently include spin corrections.
     * Note that these are not higher PN corrections to the amplitude.
     * They are the corrections to the leading-order amplitude arising
     * from the stationary phase approximation. See for instance
     * Eq 6.9 of arXiv:0810.5336
     */
	switch (amplitudeO)
        {
            case 7:
                flux += FTa7 * v7;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            case 6:
                flux += (FTa6 + FTl6*logv) * v6;
                dEnergy += dETa3 * v6;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            case 5:
                flux += FTa5 * v5;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            case 4:
                flux += FTa4 * v4;
                dEnergy += dETa2 * v4;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            case 3:
                flux += FTa3 * v3;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            case 2:
                flux += FTa2 * v2;
                dEnergy += dETa1 * v2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            case -1: /* Default to no SPA amplitude corrections */
            case 0:
                flux += 1.;
                dEnergy += 1.;
        }

        flux *= FTaN * v10;
        dEnergy *= dETaN * v;
        // Note the factor of 2 b/c phi_ref is orbital phase
        phasing += shft * f - 2.*phi_ref - ref_phasing;
        amp = amp0 * sqrt(-dEnergy/flux) * v;
        data[i+iStart] = amp * cos(phasing - LAL_PI_4)
                - amp * sin(phasing - LAL_PI_4) * 1.0j;
    }

    *htilde_out = htilde;
    return XLAL_SUCCESS;
}

/**
 * Computes the stationary phase approximation to the Fourier transform of
 * a chirp waveform with eccentric correction. The amplitude is given by expanding \f$1/\sqrt{\dot{F}}\f$.
 * If the PN order is set to -1, then the highest implemented order is used.
 *
 * @note f_ref is the GW frequency at which phi_ref is defined. The most common
 * choice in the literature is to choose the reference point as "coalescence",
 * when the frequency becomes infinite. This is the behavior of the code when
 * f_ref==0. If f_ref > 0, phi_ref sets the orbital phase at that GW frequency.
 *
 * See arXiv:0810.5336 and arXiv:astro-ph/0504538 for spin corrections
 * to the phasing.
 * See arXiv:1303.7412 for spin-orbit phasing corrections at 3 and 3.5PN order
 * See Phys. Rev. Lett. 112, 101101(2014) for eccentric phasing corrections upto 3PN order
 *
 * The spin and tidal order enums are defined in LALSimInspiralWaveformFlags.h
 */
int XLALSimInspiralTaylorF2Ecc(
        COMPLEX16FrequencySeries **htilde_out, /**< FD waveform */
        const REAL8 phi_ref,                   /**< reference orbital phase (rad) */
        const REAL8 deltaF,                    /**< frequency resolution */
        const REAL8 m1_SI,                     /**< mass of companion 1 (kg) */
        const REAL8 m2_SI,                     /**< mass of companion 2 (kg) */
        const REAL8 S1z,                       /**<  z component of the spin of companion 1 */
        const REAL8 S2z,                       /**<  z component of the spin of companion 2  */
        const REAL8 fStart,                    /**< start GW frequency (Hz) */
        const REAL8 fEnd,                      /**< highest GW frequency (Hz) of waveform generation - if 0, end at Schwarzschild ISCO */
        const REAL8 f_ref,                     /**< Reference GW frequency (Hz) - if 0 reference point is coalescence */
        const REAL8 r,                         /**< distance of source (m) */
        const REAL8 eccentricity,                       /**< eccentricity effect control < 0 : no eccentricity effect */
        LALDict *p                             /**< Linked list containing the extra parameters >**/
        )
{
    /* external: SI; internal: solar masses */
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    // const REAL8 eta = m1 * m2 / (m * m);
    const REAL8 piM = LAL_PI * m_sec;
    const REAL8 vISCO = 1. / sqrt(6.);
    const REAL8 fISCO = vISCO * vISCO * vISCO / piM;
    //const REAL8 m1OverM = m1 / m;
    // const REAL8 m2OverM = m2 / m;
    REAL8 shft, f_max;
    size_t i, n;
    INT4 iStart;
    REAL8Sequence *freqs = NULL;
    LIGOTimeGPS tC = {0, 0};
    int ret;
    int retcode;
    REAL8 fCONT;
    INT4 tideO = XLALSimInspiralWaveformParamsLookupPNTidalOrder(p);
    REAL8 lambda1 = XLALSimInspiralWaveformParamsLookupTidalLambda1(p);
    REAL8 lambda2 = XLALSimInspiralWaveformParamsLookupTidalLambda2(p);
    retcode = XLALSimInspiralSetQuadMonParamsFromLambdas(p);
    XLAL_CHECK(retcode == XLAL_SUCCESS, XLAL_EFUNC, "Failed to set quadparams from Universal relation.\n");

    COMPLEX16FrequencySeries *htilde = NULL;

    /* Perform some initial checks */
    if (!htilde_out) XLAL_ERROR(XLAL_EFAULT);
    if (*htilde_out) XLAL_ERROR(XLAL_EFAULT);
    if (m1_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
    if (f_ref < 0) XLAL_ERROR(XLAL_EDOM);
    if (r <= 0) XLAL_ERROR(XLAL_EDOM);
    if (eccentricity < 0.0 || eccentricity >= 1.0) XLAL_ERROR(XLAL_EDOM);

    /* allocate htilde */
    if ( (fEnd == 0.) && ( tideO == 0)) // End at ISCO
        f_max = fISCO;
    else if ( (fEnd == 0.) && ( tideO != 0 )) { // End at the minimum of the contact and ISCO frequencies only when tides are enabled
        fCONT = XLALSimInspiralContactFrequency(m1, lambda1, m2, lambda2); /* Contact frequency of two compact objects */
        f_max = (fCONT > fISCO) ? fISCO : fCONT;
    }
    else // End at user-specified freq.
        f_max = fEnd;
    if (f_max <= fStart) XLAL_ERROR(XLAL_EDOM);

    n = (size_t) (f_max / deltaF + 1);
    XLALGPSAdd(&tC, -1 / deltaF);  /* coalesce at t=0 */
    htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, n);
    if (!htilde) XLAL_ERROR(XLAL_EFUNC);
    memset(htilde->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitMultiply(&htilde->sampleUnits, &htilde->sampleUnits, &lalSecondUnit);

    /* Fill with non-zero vals from fStart to f_max */
    iStart = (INT4) ceil(fStart / deltaF);

    /* Sequence of frequencies where waveform model is to be evaluated */
    freqs = XLALCreateREAL8Sequence(n - iStart);

    /* extrinsic parameters */
    shft = LAL_TWOPI * (tC.gpsSeconds + 1e-9 * tC.gpsNanoSeconds);

    #pragma omp parallel for
    for (i = iStart; i < n; i++) {
        freqs->data[i-iStart] = i * deltaF;
    }

    /* phasing coefficients */
    PNPhasingSeries pfa;
    XLALSimInspiralPNPhasing_F2(&pfa, m1, m2, S1z, S2z, S1z*S1z, S2z*S2z, S1z*S2z, p);
    ret = XLALSimInspiralTaylorF2CoreEcc(&htilde, freqs, phi_ref, m1_SI, m2_SI,
                                      f_ref, shft, r, eccentricity, p, &pfa);

    XLALDestroyREAL8Sequence(freqs);

    *htilde_out = htilde;

    return ret;
}

/** @} */
/** @} */
