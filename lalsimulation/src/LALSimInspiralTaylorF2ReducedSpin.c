/*
*  Copyright (C) 2011 P. Ajith, N. Fotopoulos
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

#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <lal/Units.h>
#include <lal/XLALError.h>

#define Pi_p2 9.8696044010893586188344909998761511
#define Pi_p2by3 2.1450293971110256000774441009412356
#define log4 1.3862943611198906188344642429163531

static size_t NextPow2(const size_t n) {
  return 1 << (size_t) ceil(log2(n));
}

/**
 * Compute the dimensionless, aligned-spin parameter chi as used in the
 * TaylorF2RedSpin waveform. This is different from chi in IMRPhenomB!
 * Reference: http://arxiv.org/abs/1107.1267, paragraph 3.
 */
REAL8 XLALSimInspiralTaylorF2ReducedSpinComputeChi(
    const REAL8 m1,                          /**< mass of companion 1 */
    const REAL8 m2,                          /**< mass of companion 2 */
    const REAL8 s1z,                         /**< dimensionless spin of companion 1 */
    const REAL8 s2z                          /**< dimensionless spin of companion 2 */
) {
    const REAL8 m = m1 + m2;
    const REAL8 eta = m1 * m2 / (m * m);
    const REAL8 delta = (m1 - m2) / m;
    const REAL8 chi_s = (s1z + s2z) / 2.;
    const REAL8 chi_a = (s1z - s2z) / 2.;

    /* sanity checks */
    if (m1 <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2 <= 0) XLAL_ERROR(XLAL_EDOM);
    if (fabs(s1z) > 1) XLAL_ERROR(XLAL_EDOM);
    if (fabs(s2z) > 1) XLAL_ERROR(XLAL_EDOM);

    return chi_s * (1. - 76. * eta / 113.) + delta * chi_a;
}

/**
 * Driver routine to compute a non-precessing post-Newtonian inspiral waveform
 * in the frequency domain, described in http://arxiv.org/abs/1107.1267.
 *
 * The chi parameter should be determined from
 * XLALSimInspiralTaylorF2ReducedSpinComputeChi.
 *
 * A note from Evan Ochsner on differences with respect to TaylorF2:
 *
 * The amplitude-corrected SPA/F2 waveforms are derived and explicitly given in
 * <http://arxiv.org/abs/gr-qc/0607092> Sec. II and Appendix A (non-spinning)
 * and <http://arxiv.org/abs/0810.5336> Sec. VI and Appendix D (spin-aligned).
 *
 * The difference between F2 and F2ReducedSpin is that F2ReducedSpin always
 * keeps only the leading-order TD amplitude multiplying the 2nd harmonic (
 * A_(2,0)(t) in Eq. 2.3 of the first paper OR alpha/beta_2^(0)(t) in Eq. 6.7
 * of the second paper) but expands out the \f$1/\sqrt{\dot{F}}\f$ ( Eq. 5.3 OR Eq.
 * 6.10-6.11 resp.) to whichever order is given as 'ampO' in the code.
 *
 * On the other hand, the F2 model in the papers above will PN expand BOTH the
 * TD amplitude and the factor \f$1/\sqrt{\dot{F}}\f$, take their product, and keep
 * all terms up to the desired amplitude order, as in Eq. 6.13-6.14 of the
 * second paper.
 *
 * In particular, the F2ReducedSpin will always have only the 2nd harmonic, but
 * F2 will have multiple harmonics starting at ampO = 0.5PN. Even if you were
 * to compare just the 2nd harmonic, you would have a difference starting at
 * 1PN ampO, because the F2 has a 1PN TD amp. correction to the 2nd harmonic
 * (alpha/beta_2^(2)(t)) which will not be accounted for by the F2ReducedSpin.
 * So, the two should agree when ampO=0, but will be different in any other
 * case.
 */
int XLALSimInspiralTaylorF2ReducedSpin(
    COMPLEX16FrequencySeries **htilde,   /**< FD waveform */
    const REAL8 phic,                /**< orbital coalescence phase (rad) */
    const REAL8 deltaF,              /**< frequency resolution (Hz) */
    const REAL8 m1_SI,               /**< mass of companion 1 (kg) */
    const REAL8 m2_SI,               /**< mass of companion 2 (kg) */
    const REAL8 chi,                 /**< dimensionless aligned-spin param */
    const REAL8 fStart,              /**< start GW frequency (Hz) */
    const REAL8 r,                   /**< distance of source (m) */
    const INT4 phaseO,               /**< twice PN phase order */
    const INT4 ampO                  /**< twice PN amplitude order */
    ) {
    /* external: SI; internal: solar masses */
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 eta = m1 * m2 / (m * m);
    const REAL8 etap2 = eta * eta;
    const REAL8 etap3 = etap2 * eta;
    const REAL8 piM = LAL_PI * m_sec;
    const REAL8 mSevenBySix = -7./6.;
    const REAL8 vISCO = 1. / sqrt(6.);
    const REAL8 fISCO = vISCO * vISCO * vISCO / piM;
    REAL8 v0 = cbrt(piM * fStart);
    REAL8 logv0 = log(v0);
    REAL8 shft, amp0, f_max;
    REAL8 psiNewt, psi2, psi3, psi4, psi5, psi6, psi6L, psi7, psi3S, psi4S, psi5S;
    REAL8 alpha2, alpha3, alpha4, alpha5, alpha6, alpha6L, alpha7, alpha3S, alpha4S, alpha5S;
    REAL8 eta_fac = -113. + 76. * eta;
    size_t i, n, iStart, iISCO;
    COMPLEX16 *data = NULL;
    LIGOTimeGPS tStart = {0, 0};

    /* check inputs for sanity */
    if (*htilde) XLAL_ERROR(XLAL_EFAULT);
    if (deltaF <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m1_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (fabs(chi) > 1) XLAL_ERROR(XLAL_EDOM);
    if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
    if (r <= 0) XLAL_ERROR(XLAL_EDOM);
    if (ampO > 7) XLAL_ERROR(XLAL_EDOM); /* only implemented to pN 3.5 */
    if (phaseO > 7) XLAL_ERROR(XLAL_EDOM); /* only implemented to pN 3.5 */

    /* allocate htilde */
    f_max = NextPow2(fISCO);
    n = f_max / deltaF + 1;
    XLALGPSAdd(&tStart, -1 / deltaF);  /* coalesce at t=0 */
    *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &tStart, 0.0, deltaF, &lalStrainUnit, n);
    if (!(*htilde)) XLAL_ERROR(XLAL_EFUNC);
    memset((*htilde)->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitDivide(&((*htilde)->sampleUnits), &((*htilde)->sampleUnits), &lalSecondUnit);

    /* extrinsic parameters */
    amp0 = -pow(m_sec, 5./6.) * sqrt(5.*eta / 24.) / (Pi_p2by3 * r / LAL_C_SI);
    shft = LAL_TWOPI * (tStart.gpsSeconds + 1e-9 * tStart.gpsNanoSeconds);

    /* spin terms in the amplitude and phase (in terms of the reduced
     * spin parameter */
    psi3S = 113.*chi/3.;
    psi4S = 63845.*(-81. + 4.*eta)*chi*chi/(8. * eta_fac * eta_fac);
    psi5S = -565.*(-146597. + 135856.*eta + 17136.*etap2)*chi/(2268.*eta_fac);

    alpha3S = (113.*chi)/24.;
    alpha4S = (12769.*chi*chi*(-81. + 4.*eta))/(32. * eta_fac * eta_fac);
    alpha5S = (-113.*chi*(502429. - 591368.*eta + 1680*etap2))/(16128.*eta_fac);

    /* coefficients of the phase at PN orders from 0 to 3.5PN */
    psiNewt = 3./(128.*eta);
    psi2 = 3715./756. + 55.*eta/9.;
    psi3 = psi3S - 16.*LAL_PI;
    psi4 = 15293365./508032. + 27145.*eta/504. + 3085.*eta*eta/72. + psi4S;
    psi5 = (38645.*LAL_PI/756. - 65.*LAL_PI*eta/9. + psi5S);
    psi6 = 11583231236531./4694215680. - (640.*Pi_p2)/3. - (6848.*LAL_GAMMA)/21.
             + (-5162.983708047263 + 2255.*Pi_p2/12.)*eta
             + (76055.*etap2)/1728. - (127825.*etap3)/1296.;
    psi6L = -6848./21.;
    psi7 = (77096675.*LAL_PI)/254016. + (378515.*LAL_PI*eta)/1512.
             - (74045.*LAL_PI*eta*eta)/756.;

    /* amplitude coefficients */
    alpha2 = 1.1056547619047619 + (11*eta)/8.;
    alpha3 = -LAL_TWOPI + alpha3S;
    alpha4 = 0.8939214212884228 + (18913*eta)/16128. + (1379*etap2)/1152. + alpha4S;
    alpha5 = (-4757*LAL_PI)/1344. + (57*eta*LAL_PI)/16. + alpha5S;
    alpha6 = -58.601030974347324 + (3526813753*eta)/2.7869184e7 -
                (1041557*etap2)/258048. + (67999*etap3)/82944. +
                (10*Pi_p2)/3. - (451*eta*Pi_p2)/96.;
    alpha6L = 856/105.;
    alpha7 = (-5111593*LAL_PI)/2.709504e6 - (72221*eta*LAL_PI)/24192. -
                (1349*etap2*LAL_PI)/24192.;

    /* select the terms according to the PN order chosen */
    switch (ampO) {
        case 0:
        case 1:
            alpha2 = 0.;
        case 2:
            alpha3 = 0.;
        case 3:
            alpha4 = 0.;
        case 4:
            alpha5 = 0.;
        case 5:
            alpha6 = 0.;
            alpha6L = 0.;
        case 6:
            alpha7 = 0.;
        default:
            break;
    }

    switch (phaseO) {
        case 0:
        case 1:
            psi2 = 0.;
        case 2:
            psi3 = 0.;
        case 3:
            psi4 = 0.;
        case 4:
            psi5 = 0.;
        case 5:
            psi6 = 0.;
            psi6L = 0.;
        case 6:
            psi7 = 0.;
        default:
            break;
    }

    iStart = (size_t) ceil(fStart / deltaF);
    iISCO = (size_t) (fISCO / deltaF);
    iISCO = (iISCO < n) ? iISCO : n;  /* overflow protection; should we warn? */
    data = (*htilde)->data->data;
    for (i = iStart; i < iISCO; i++) {
        /* fourier frequency corresponding to this bin */
        const REAL8 f = i * deltaF;
        const REAL8 v3 = piM*f;

        /* PN expansion parameter */
        REAL8 v, v2, v4, v5, v6, v7, logv, Psi, amp;
        v = cbrt(v3);
        v2 = v*v; v4 = v3*v; v5 = v4*v; v6 = v3*v3; v7 = v6*v;
        logv = log(v);

        /* compute the phase and amplitude */
        Psi = psiNewt / v5 * (1.
            + psi2 * v2 + psi3 * v3 + psi4 * v4
            + psi5 * v5 * (1. + 3. * (logv - logv0))
            + (psi6 + psi6L * (log4 + logv)) * v6 + psi7 * v7);

        amp = amp0 * pow(f, mSevenBySix) * (1.
            + alpha2 * v2 + alpha3 * v3 + alpha4 * v4 + alpha5 * v5
            + (alpha6 + alpha6L * (LAL_GAMMA + log4 + logv)) * v6
            + alpha7 * v7);

        data[i] = amp * cos(Psi + shft * f - 2.*phic - LAL_PI_4)
            - I * (amp * sin(Psi + shft * f - 2.*phic - LAL_PI_4));
    }

    return XLAL_SUCCESS;
}

/**
Compute the chirp time of the "reduced-spin" templates
*/
REAL8 XLALSimInspiralTaylorF2ReducedSpinChirpTime(
    const REAL8 fStart,  /**< start GW frequency (Hz) */
    const REAL8 m1_SI,   /**< mass of companion 1 (kg) */
    const REAL8 m2_SI,   /**< mass of companion 2 (kg) */
    const REAL8 chi,     /**< dimensionless aligned-spin param */
    const INT4 O         /**< twice PN phase order */
    ) {
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1+m2;
    const REAL8 eta = m1*m2/(m*m);
    const REAL8 eta2 = eta*eta;
    const REAL8 chi2 = chi*chi;
    const REAL8 sigma0 = (-12769*(-81. + 4.*eta))/(16.*(-113. + 76.*eta)*(-113. + 76.*eta));
    const REAL8 gamma0 = (565*(-146597. + 135856.*eta + 17136.*eta2))/(2268.*(-113. + 76.*eta));

    REAL8 v = cbrt(LAL_PI * m * LAL_MTSUN_SI * fStart);
    REAL8 vk = v;  /* v^k */
    REAL8 tau = 1.;  /* chirp time */
    REAL8 tk[8];  /* chirp time coefficients up to 3.5 PN */
    size_t k;
    UINT4 order = (UINT4) O;

    if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m1_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (fabs(chi) > 1) XLAL_ERROR(XLAL_EDOM);
    if (O > 7) XLAL_ERROR(XLAL_EDOM); /* only implemented to pN 3.5 */
    if (O == -1) order = 7; /* to be changed to use #define MAX_PHASE_ORDER in LALSimInspiral.h */

    /* chirp time coefficients up to 3.5PN  */
    tk[0] = (5.*m*LAL_MTSUN_SI)/(256.*pow(v,8)*eta);
    tk[1] = 0.;
    tk[2] = 2.9484126984126986 + (11*eta)/3.;
    tk[3] = (-32*LAL_PI)/5. + (226.*chi)/15.;
    tk[4] = 6.020630590199042 - 2*sigma0*chi2 + (5429*eta)/504. + (617*eta2)/72.;
    tk[5] = (3*gamma0*chi)/5. - (7729*LAL_PI)/252. + (13*LAL_PI*eta)/3.;
    tk[6] = -428.291776175525 + (128*Pi_p2)/3. + (6848*LAL_GAMMA)/105. + (3147553127*eta)/3.048192e6 -
       (451*Pi_p2*eta)/12. - (15211*eta2)/1728. + (25565*eta2*eta)/1296. + (6848*log(4*v))/105.;
    tk[7] = (-15419335*LAL_PI)/127008. - (75703*LAL_PI*eta)/756. + (14809*LAL_PI*eta2)/378.;

    /* compute chirp time */
    for (k = 2; k <= order; k++) {
        vk *= v;
        tau += tk[k] * vk;
    }
    tau *= tk[0];

    return tau;
}
