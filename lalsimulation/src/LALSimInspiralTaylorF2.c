/*
 *  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
 *  Copyright (C) 2012 Leo Singer, Evan Ochsner, Les Wade, Alex Nitz
 *  Assembled from code found in:
 *    - LALInspiralStationaryPhaseApproximation2.c
 *    - LALInspiralChooseModel.c
 *    - LALInspiralSetup.c
 *    - LALSimInspiralTaylorF2ReducedSpin.c
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

#include <stdlib.h>
#include <math.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <lal/Units.h>
#include <lal/XLALError.h>
#include "LALSimInspiralPNCoefficients.c"

int XLALSimInspiralTaylorF2Phasing(
        PNPhasingSeries **pn,
        const REAL8 m1,
        const REAL8 m2,
        const REAL8 chi1L,
        const REAL8 chi2L,
        const LALSimInspiralSpinOrder spinO
	)
{
    PNPhasingSeries *pfa;

    if (!pn) XLAL_ERROR(XLAL_EFAULT);
    if (*pn) XLAL_ERROR(XLAL_EFAULT);


    pfa = (PNPhasingSeries *) LALMalloc(sizeof(PNPhasingSeries));

    XLALSimInspiralPNPhasing_F2WithSO(pfa, m1, m2, chi1L, chi2L, spinO);

    *pn = pfa;

    return XLAL_SUCCESS;
}


/**
 * Computes the stationary phase approximation to the Fourier transform of
 * a chirp waveform with phase given by \eqref{eq_InspiralFourierPhase_f2}
 * and amplitude given by expanding \f$1/\sqrt{\dot{F}}\f$. If the PN order is
 * set to -1, then the highest implemented order is used.
 *
 * See arXiv:0810.5336 and arXiv:astro-ph/0504538 for spin corrections
 * to the phasing.
 * See arXiv:1303.7412 for spin-orbit phasing corrections at 3 and 3.5PN order
 */
int XLALSimInspiralTaylorF2(
        COMPLEX16FrequencySeries **htilde_out, /**< FD waveform */
        const REAL8 phi_ref,                   /**< reference orbital phase (rad) */
        const REAL8 deltaF,                    /**< frequency resolution */
        const REAL8 m1_SI,                     /**< mass of companion 1 (kg) */
        const REAL8 m2_SI,                     /**< mass of companion 2 (kg) */
        const REAL8 S1z,                       /**<  z component of the spin of companion 1 */
        const REAL8 S2z,                       /**<  z component of the spin of companion 2  */
        const REAL8 fStart,                    /**< start GW frequency (Hz) */
        const REAL8 f_ref,                     /**< Reference GW frequency (Hz) - if 0 reference point is coalescence */
        const REAL8 fEnd,                      /**< highest GW frequency (Hz) of waveform generation - if 0, end at Schwarzschild ISCO */
        const REAL8 r,                         /**< distance of source (m) */
        const REAL8 lambda1,                   /**< (tidal deformation of body 1)/(mass of body 1)^5 */
        const REAL8 lambda2,                   /**< (tidal deformation of body 2)/(mass of body 2)^5 */
        const LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
        const LALSimInspiralTidalOrder tideO,  /**< flag to control tidal effects */
        const INT4 phaseO,                     /**< twice PN phase order */
        const INT4 amplitudeO                  /**< twice PN amplitude order */
        )
{
    const REAL8 lambda = -1987./3080.;
    const REAL8 theta = -11831./9240.;

    /* external: SI; internal: solar masses */
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 eta = m1 * m2 / (m * m);
    const REAL8 piM = LAL_PI * m_sec;
    const REAL8 vISCO = 1. / sqrt(6.);
    const REAL8 fISCO = vISCO * vISCO * vISCO / piM;
    const REAL8 v0 = 1.0; /* v0=c */
    const REAL8 chi1 = m1 / m;
    const REAL8 chi2 = m2 / m;
    const REAL8 lam1 = lambda1;
    const REAL8 lam2 = lambda2;
    REAL8 shft, amp0, f_max;
    size_t i, n, iStart;
    COMPLEX16 *data = NULL;
    LIGOTimeGPS tC = {0, 0};

    /* phasing coefficients */
    const REAL8 pfaN = 3.L/(128.L * eta);
    const REAL8 pfa2 = 5.L*(743.L/84.L + 11.L * eta)/9.L;
    const REAL8 pfa3 = -16.L*LAL_PI;
    const REAL8 pfa4 = 5.L*(3058.673L/7.056L + 5429.L/7.L * eta
                     + 617.L * eta*eta)/72.L;
    const REAL8 pfa5 = 5.L/9.L * (7729.L/84.L - 13.L * eta) * LAL_PI;
    const REAL8 pfl5 = 5.L/3.L * (7729.L/84.L - 13.L * eta) * LAL_PI;
    const REAL8 pfa6 = (11583.231236531L/4.694215680L
                     - 640.L/3.L * LAL_PI * LAL_PI - 6848.L/21.L*LAL_GAMMA)
                     + eta * (-15335.597827L/3.048192L
                     + 2255./12. * LAL_PI * LAL_PI
                     - 1760./3.*theta +12320./9.*lambda)
                     + eta*eta * 76055.L/1728.L - eta*eta*eta * 127825.L/1296.L;
    const REAL8 pfl6 = -6848.L/21.L;
    const REAL8 pfa7 = LAL_PI * 5.L/756.L * ( 15419335.L/336.L
                     + 75703.L/2.L * eta - 14809.L * eta*eta);

    /* Spin coefficients */
    REAL8 pn_beta = 0;
    REAL8 pn_sigma = 0;
    REAL8 pn_gamma = 0;
    REAL8 psiSO3 = 0., psiSO35 = 0.; // 3PN and 3.5PN spin-orbit phasing terms

    REAL8 d = (m1 - m2) / (m1 + m2);
    REAL8 xs = .5 * (S1z + S2z);
    REAL8 xa = .5 * (S1z - S2z);

    REAL8 qm_def1 = 1; /* The QM deformability parameters */
    REAL8 qm_def2 = 1; /* This is 1 for black holes and larger for neutron stars */

    /* Validate expansion order arguments.
     * This must be done here instead of in the OpenMP parallel loop
     * because when OpenMP parallelization is turned on, early exits
     * from loops (via return or break statements) are not permitted.
     */

    /* Validate phase PN order. */
    switch (phaseO)
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
            XLAL_ERROR(XLAL_ETYPE, "Invalid phase PN order %s", phaseO);
    }

    /* Validate amplitude PN order. */
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
            XLAL_ERROR(XLAL_ETYPE, "Invalid amplitude PN order %s", amplitudeO);
    }

    switch( spinO )
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
            psiSO35 = (chi1 * (-8980424995./762048. + 6586595.*eta/756.
                    - 305.*eta*eta/36.) + d * (170978035./48384.
                    - 2876425.*eta/672. - 4735.*eta*eta/144.) ) * chi1 * S1z
                    + (chi2 * (-8980424995./762048. + 6586595.*eta/756.
                    - 305.*eta*eta/36.) - d * (170978035./48384.
                    - 2876425.*eta/672. - 4735.*eta*eta/144.) ) * chi2 * S2z;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
            psiSO3 = LAL_PI * ( (260.*chi1 + 1490./3.) * chi1 * S1z
                    + (260.*chi2 + 1490./3.) * chi2 * S2z);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            /* Compute 2.5PN SO correction */
            // See Eq. (6.25) in arXiv:0810.5336
            pn_gamma = (732985.L/2268.L - 24260.L/81.L * eta - 340.L/9.L * eta * eta ) * xs;
            pn_gamma += (732985.L/2268.L +140.L/9.0L * eta) * xa * d;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            /* Compute 2.0PN SS, QM, and self-spin */
            // See Eq. (6.24) in arXiv:0810.5336
            // 9b,c,d in arXiv:astro-ph/0504538
            pn_sigma = eta * (721.L/48.L *S1z*S2z-247.L/48.L*S1z*S2z);
            pn_sigma += (720*qm_def1 - 1)/96.0 * (chi1*chi1*S1z*S1z);
            pn_sigma += (720*qm_def2 - 1)/96.0 * (chi2*chi2*S2z*S2z);
            pn_sigma -= (240*qm_def1 - 7)/96.0 * (chi1*chi1*S1z*S1z);
            pn_sigma -= (240*qm_def2 - 7)/96.0 * (chi2*chi2*S2z*S2z);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            /* Compute 1.5PN SO correction */
            // Eq. (6.23) in arXiv:0810.5336
            pn_beta = (113.L/12.L- 19.L/3.L * eta) * xs + 113.L/12.L * d * xa;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLAL_ERROR(XLAL_EINVAL, "Invalid spin PN order %s", spinO);
    }

    /* Tidal coefficients for phasing, fluz, and energy */
    REAL8 pft10 = 0.;
    REAL8 pft12 = 0.;
    switch( tideO )
    {
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
            pft12 = - 5.L * lam2 * chi2*chi2*chi2*chi2 * (3179.L - 919.L*chi2
                    - 2286.L*chi2*chi2 + 260.L*chi2*chi2*chi2)/28.L
                    - 5.L * lam1 * chi1*chi1*chi1*chi1 * (3179.L - 919.L*chi1
                    - 2286.L*chi1*chi1 + 260.L*chi1*chi1*chi1)/28.L;
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
            pft10 = - 24.L * lam2 * chi2*chi2*chi2*chi2 * (1.L + 11.L*chi1)
                    - 24.L * lam1 * chi1*chi1*chi1*chi1 * (1.L + 11.L*chi2);
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLAL_ERROR(XLAL_EINVAL, "Invalid tidal PN order %s", tideO);
    }

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


    COMPLEX16FrequencySeries *htilde;

    /* Perform some initial checks */
    if (!htilde_out) XLAL_ERROR(XLAL_EFAULT);
    if (*htilde_out) XLAL_ERROR(XLAL_EFAULT);
    if (m1_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
    if (r <= 0) XLAL_ERROR(XLAL_EDOM);

    /* allocate htilde */
    if ( fEnd == 0. ) // End at ISCO
        f_max = fISCO;
    else // End at user-specified freq.
        f_max = fEnd;
    n = (size_t) (f_max / deltaF + 1);
    XLALGPSAdd(&tC, -1 / deltaF);  /* coalesce at t=0 */
    htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, n);
    if (!htilde) XLAL_ERROR(XLAL_EFUNC);
    memset(htilde->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitDivide(&htilde->sampleUnits, &htilde->sampleUnits, &lalSecondUnit);

    /* extrinsic parameters */
    amp0 = -4. * m1 * m2 / r * LAL_MRSUN_SI * LAL_MTSUN_SI * sqrt(LAL_PI/12.L);
    shft = LAL_TWOPI * (tC.gpsSeconds + 1e-9 * tC.gpsNanoSeconds);

    const REAL8 log4=log(4.0);
    const REAL8 logv0=log(v0);
    
    /* Fill with non-zero vals from fStart to f_max */
    iStart = (size_t) ceil(fStart / deltaF);
    data = htilde->data->data;

    /* Compute the SPA phase at the reference point */
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
    REAL8 ref_phasing = 0.;
    switch (phaseO)
    {
        case -1:
        case 7:
            ref_phasing += pfa7 * v7ref;
        case 6:
            ref_phasing += (pfa6 + pfl6 * (log4+logvref)) * v6ref;
        case 5:
            ref_phasing += (pfa5 + pfl5 * (logvref-logv0)) * v5ref;
        case 4:
            ref_phasing += pfa4 * v4ref;
        case 3:
            ref_phasing += pfa3 * v3ref;
        case 2:
            ref_phasing += pfa2 * v2ref;
        case 0:
            ref_phasing += 1.;
    }
    switch( spinO )
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
            ref_phasing += psiSO35 * v7ref;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
            ref_phasing += psiSO3 * v6ref;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            ref_phasing += -pn_gamma * (1 + 3*(logvref-logv0)) * v5ref;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            ref_phasing += -10.L*pn_sigma * v4ref;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            ref_phasing += 4.L*pn_beta * v3ref;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            ;
    }
    switch( tideO )
    {
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
            ref_phasing += pft12 * v12ref;
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
            ref_phasing += pft10 * v10ref;
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            ;
    }
    ref_phasing *= pfaN / v5ref;

    #pragma omp parallel for
    for (i = iStart; i < n; i++) {
        const REAL8 f = i * deltaF;
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
        REAL8 phasing = 0.;
        REAL8 dEnergy = 0.;
        REAL8 flux = 0.;
        REAL8 amp;

        switch (phaseO)
        {
            case -1:
            case 7:
                phasing += pfa7 * v7;
            case 6:
                phasing += (pfa6 + pfl6 * (log4+logv)) * v6;
            case 5:
                phasing += (pfa5 + pfl5 * (logv-logv0)) * v5;
            case 4:
                phasing += pfa4 * v4;
            case 3:
                phasing += pfa3 * v3;
            case 2:
                phasing += pfa2 * v2;
            case 0:
                phasing += 1.;
        }
        switch (amplitudeO)
        {
            case -1:
            case 7:
                flux += FTa7 * v7;
            case 6:
                flux += (FTa6 + FTl6*2.0*(log4+logv)) * v6;
                dEnergy += dETa3 * v6;
            case 5:
                flux += FTa5 * v5;
            case 4:
                flux += FTa4 * v4;
                dEnergy += dETa2 * v4;
            case 3:
                flux += FTa3 * v3;
            case 2:
                flux += FTa2 * v2;
                dEnergy += dETa1 * v2;
            case 0:
                flux += 1.;
                dEnergy += 1.;
        }

        switch( spinO )
        {
            case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
            case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
                phasing += psiSO35 * v7;
            case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
                phasing += psiSO3 * v6;
            case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
                phasing += -pn_gamma * (1 + 3*(logv-logv0)) * v5;
            case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
                phasing += -10.L*pn_sigma * v4;
            case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
                phasing += 4.L*pn_beta * v3;
            case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
            case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
            case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
                ;
        }

        switch( tideO )
        {
            case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
            case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
                phasing += pft12 * v12;
            case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
                phasing += pft10 * v10;
            case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
                ;
        }

        phasing *= pfaN / v5;
        flux *= FTaN * v10;
        dEnergy *= dETaN * v;
        // Note the factor of 2 b/c phi_ref is orbital phase
        phasing += shft * f - 2.*phi_ref - ref_phasing;
        amp = amp0 * sqrt(-dEnergy/flux) * v;
        data[i] = amp * cos(phasing - LAL_PI_4)
                - amp * sin(phasing - LAL_PI_4) * 1.0j;
    }

    *htilde_out = htilde;
    return XLAL_SUCCESS;
}
