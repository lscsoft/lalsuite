/*
 *  Copyright (C) 2014 Andrew Lundgren
 *  Based on code in LALSimInspiralTaylorF2.c
 *  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
 *  Copyright (C) 2012 Leo Singer, Evan Ochsner, Les Wade, Alex Nitz
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
#include <stdio.h>
#include <stdbool.h>

typedef struct tagLALSimInspiralSF2Orientation
{
    REAL8 thetaJ, psiJ;
    REAL8 chi, kappa, alpha0;
    REAL8 v_ref;
} LALSimInspiralSF2Orientation;

typedef struct tagLALSimInspiralSF2Coeffs
{
    REAL8 m1, m2, mtot, eta;
    REAL8 kappa, kappa_perp, gamma0;
    REAL8 prec_fac;
    REAL8 aclog1, aclog2;
    REAL8 ac[6];
    REAL8 zc[7];
} LALSimInspiralSF2Coeffs;

// Prototypes

REAL8 safe_atan2(REAL8 val1, REAL8 val2);

void XLALSimInspiralSF2CalculateOrientation(
    LALSimInspiralSF2Orientation *orientation,
    REAL8 m1, REAL8 m2, REAL8 f_ref,
    REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz,
    REAL8 s1x, REAL8 s1y, REAL8 s1z);

void XLALSimInspiralSF2CalculateCoeffs(
    LALSimInspiralSF2Coeffs *coeffs,
    REAL8 m1, REAL8 m2, REAL8 chi, REAL8 kappa);

REAL8 XLALSimInspiralSF2Alpha(REAL8 v, LALSimInspiralSF2Coeffs coeffs);
REAL8 XLALSimInspiralSF2Beta(REAL8 v, LALSimInspiralSF2Coeffs coeffs);
REAL8 XLALSimInspiralSF2Zeta(REAL8 v, LALSimInspiralSF2Coeffs coeffs);

COMPLEX16 XLALSimInspiralSF2Polarization(
    REAL8 thetaJ, REAL8 psiJ, int mm);

void XLALSimInspiralSF2Emission(
    REAL8 *emission, REAL8 v, LALSimInspiralSF2Coeffs coeffs);

REAL8 safe_atan2(REAL8 val1, REAL8 val2)
{
    if (val1 == 0. && val2 == 0.)
    { return 0.; }
    else
    { return atan2(val1, val2); }
}

void XLALSimInspiralSF2CalculateOrientation(
    LALSimInspiralSF2Orientation *orientation,
    REAL8 m1, REAL8 m2, REAL8 v_ref,
    REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz,
    REAL8 s1x, REAL8 s1y, REAL8 s1z )
{
    orientation->v_ref = v_ref;
    const REAL8 chi1 = sqrt(s1x*s1x+s1y*s1y+s1z*s1z);
    orientation->chi = chi1;
    orientation->kappa = (chi1 > 0.) ? (lnhatx*s1x+lnhaty*s1y+lnhatz*s1z)/chi1 : 1.;
    const REAL8 Jx0 = m1*m2*lnhatx/v_ref + m1*m1*s1x;
    const REAL8 Jy0 = m1*m2*lnhaty/v_ref + m1*m1*s1y;
    const REAL8 Jz0 = m1*m2*lnhatz/v_ref + m1*m1*s1z;
    const REAL8 thetaJ = acos(Jz0 / sqrt(Jx0*Jx0+Jy0*Jy0+Jz0*Jz0));
    orientation->thetaJ = thetaJ;
    const REAL8 psiJ = safe_atan2(Jy0, -Jx0);
    orientation->psiJ = psiJ;

    /* Rotate Lnhat back to frame where J is along z, to figure out initial alpha */
    const REAL8 rotLx = lnhatx*cos(thetaJ)*cos(psiJ) - lnhaty*cos(thetaJ)*sin(psiJ) + lnhatz*sin(thetaJ);
    const REAL8 rotLy = lnhatx*sin(psiJ) + lnhaty*cos(psiJ);
    orientation->alpha0 = safe_atan2(rotLy, rotLx);
}

void XLALSimInspiralSF2CalculateCoeffs(
    LALSimInspiralSF2Coeffs *coeffs,
    REAL8 m1, REAL8 m2, REAL8 chi, REAL8 kappa )
{
    const REAL8 quadparam = 1.;
    coeffs->m1 = m1;
    coeffs->m2 = m2;
    const REAL8 mtot = m1+m2;
    coeffs->mtot = mtot;
    const REAL8 eta = m1*m2/(mtot*mtot);
    coeffs->eta = eta;
    const REAL8 gamma0 = m1*chi/m2;
    coeffs->kappa = kappa;
    coeffs->kappa_perp = sqrt(1.-kappa*kappa);
    coeffs->gamma0 = gamma0;

    const REAL8 pn_beta = (113.*m1/(12.*mtot) - 19.*eta/6.)*chi*kappa;
    const REAL8 pn_sigma = (  (5.*quadparam*(3.*kappa*kappa-1.)/2.)
                           + (7. - kappa*kappa)/96.  )
                       * (m1*m1*chi*chi/mtot/mtot);
    const REAL8 pn_gamma = (5.*(146597. + 7056.*eta)*m1/(2268.*mtot)
                        - 10.*eta*(1276. + 153.*eta)/81.)*chi*kappa;

    coeffs->prec_fac = 5.*(4. + 3.*m2/m1)/64.;
    const REAL8 dtdv2 = 743./336. + 11.*eta/4.;
    const REAL8 dtdv3 = -4.*LAL_PI + pn_beta;
    const REAL8 dtdv4 = 3058673./1016064. + 5429.*eta/1008. + 617.*eta*eta/144. - pn_sigma;
    const REAL8 dtdv5 = (-7729./672.+13.*eta/8.)*LAL_PI + 9.*pn_gamma/40.;

    coeffs->aclog1 = kappa*(1.-kappa*kappa)*gamma0*gamma0*gamma0/2.-dtdv2*kappa*gamma0-dtdv3;
    coeffs->aclog2 = dtdv2*gamma0+dtdv3*kappa+(1.-kappa*kappa)*(dtdv4-dtdv5*kappa/gamma0)/gamma0/2.;
    coeffs->ac[0] = -1./3.;
    coeffs->ac[1] = -gamma0*kappa/6.;
    coeffs->ac[2] = gamma0*gamma0*(-1./3.+kappa*kappa/2.) - dtdv2;
    coeffs->ac[3] = dtdv3+dtdv4*kappa/gamma0/2.+dtdv5*(1./3.-kappa*kappa/2.)/gamma0/gamma0;
    coeffs->ac[4] = dtdv4/2.+dtdv5*kappa/gamma0/6.;
    coeffs->ac[5] = dtdv5/3.;

    coeffs->zc[0] = -1./3.;
    coeffs->zc[1] = -gamma0*kappa/2.;
    coeffs->zc[2] = -dtdv2;
    coeffs->zc[3] = dtdv2*gamma0*kappa+dtdv3;
    coeffs->zc[4] = dtdv3*gamma0*kappa+dtdv4;
    coeffs->zc[5] = (dtdv4*gamma0*kappa+dtdv5)/2.;
    coeffs->zc[6] = dtdv5*gamma0*kappa/3.;
}

REAL8 XLALSimInspiralSF2Alpha(
    REAL8 v, LALSimInspiralSF2Coeffs coeffs)
{
    const REAL8 gam = coeffs.gamma0*v;
    const REAL8 kappa = coeffs.kappa;

    const REAL8 sqrtfac = sqrt(1. + 2.*kappa*gam + gam*gam);
    const REAL8 logfac1 = log((1. + kappa*gam + sqrtfac)/v);
    const REAL8 logfac2 = log(kappa + gam + sqrtfac);

    const REAL8 aclog1 = coeffs.aclog1;
    const REAL8 aclog2 = coeffs.aclog2;
    const REAL8 *ac = coeffs.ac;

    return coeffs.prec_fac * (aclog1*logfac1 + aclog2*logfac2 \
                    + (((ac[0]/v+ac[1])/v+ac[2])/v + ac[3] \
                        + (ac[4]+ac[5]*v)*v)*sqrtfac );
}

REAL8 XLALSimInspiralSF2Zeta(
    REAL8 v, LALSimInspiralSF2Coeffs coeffs)
{
    const REAL8 *zc = coeffs.zc;

    return coeffs.prec_fac*(((zc[0]/v+zc[1])/v+zc[2])/v + zc[3]*log(v) + \
                        (zc[4]+(zc[5]+zc[6]*v)*v)*v);
}

REAL8 XLALSimInspiralSF2Beta(
    REAL8 v, LALSimInspiralSF2Coeffs coeffs)
{
    const REAL8 kappa = coeffs.kappa;
    const REAL8 kappa_perp = coeffs.kappa_perp;
    const REAL8 gamma0 = coeffs.gamma0;

    REAL8 beta = safe_atan2(gamma0*v*kappa_perp, 1. + kappa*gamma0*v);

    return beta;
}

COMPLEX16 XLALSimInspiralSF2Polarization(
    REAL8 thetaJ, REAL8 psiJ, int mm)
{
    COMPLEX16 plus_fac, cross_fac;
    switch (mm)
    {
        case 2:
            plus_fac = (1.+cos(thetaJ)*cos(thetaJ))/2.;
            cross_fac = -1.j*cos(thetaJ);
            break;
        case 1:
            plus_fac = sin(2.*thetaJ);
            cross_fac = -2.j*sin(thetaJ);
            break;
        case 0:
            plus_fac = 3.*sin(thetaJ)*sin(thetaJ);
            cross_fac = 0.;
            break;
        case -1:
           plus_fac = -sin(2.*thetaJ);
            cross_fac = -2.j*sin(thetaJ);
            break;
        case -2:
            plus_fac = (1.+cos(thetaJ)*cos(thetaJ))/2.;
            cross_fac = 1.j*cos(thetaJ);
            break;
        default:
            plus_fac = 0.;
            cross_fac = 0.;
    }
    
    return plus_fac*cos(2.*psiJ) + cross_fac*sin(2.*psiJ);
}

void XLALSimInspiralSF2Emission(
    REAL8 *emission, REAL8 v, LALSimInspiralSF2Coeffs coeffs)
{
    const REAL8 gam = coeffs.gamma0*v;
    const REAL8 kappa = coeffs.kappa;
    const REAL8 kappa_perp = coeffs.kappa_perp;

    const REAL8 sqrtfac = sqrt(1. + 2.*kappa*gam + gam*gam);
    const REAL8 cosbeta = (1. + kappa*gam)/sqrtfac;
    const REAL8 sinbeta = (kappa_perp*gam)/sqrtfac;

    emission[0] = (1.+cosbeta)*(1.+cosbeta)/4.;
    emission[1] = (1.+cosbeta)*sinbeta/4.;
    emission[2] = sinbeta*sinbeta/4.;
    emission[3] = (1.-cosbeta)*sinbeta/4.;
    emission[4] = (1.-cosbeta)*(1.-cosbeta)/4.;

    return; 
}


/** FIXME Is this needed?
 * Find the least nonnegative integer power of 2 that is
 * greater than or equal to n.  Inspired by similar routine
 * in gstlal.
static size_t CeilPow2(double n) {
    double signif;
    int exponent;
    signif = frexp(n, &exponent);
    if (signif < 0)
        return 1;
    if (signif == 0.5)
        exponent -= 1;
    return ((size_t) 1) << exponent;
} */ 

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
int XLALSimInspiralSpinTaylorF2(
        COMPLEX16FrequencySeries **hplus_out,  /**< FD hplus waveform */
        COMPLEX16FrequencySeries **hcross_out, /**< FD hcross waveform */
        const REAL8 phi_ref,                   /**< reference orbital phase (rad) */
        const REAL8 deltaF,                    /**< frequency resolution */
        const REAL8 m1_SI,                     /**< mass of companion 1 (kg) */
        const REAL8 m2_SI,                     /**< mass of companion 2 (kg) */
        const REAL8 s1x,                             /**< initial value of S1x */
        const REAL8 s1y,                             /**< initial value of S1y */
        const REAL8 s1z,                             /**< initial value of S1z */
        const REAL8 lnhatx,                          /**< initial value of LNhatx */
        const REAL8 lnhaty,                          /**< initial value of LNhaty */
        const REAL8 lnhatz,                          /**< initial value of LNhatz */
        const REAL8 fStart,                    /**< start GW frequency (Hz) */
        const REAL8 fEnd,                      /**< highest GW frequency (Hz) of waveform generation - if 0, end at Schwarzschild ISCO */
        const REAL8 f_ref,                     /**< Reference GW frequency (Hz) - if 0 reference point is coalescence */
        const REAL8 r,                         /**< distance of source (m) */
        LALSimInspiralTestGRParam *moreParams, /**< Linked list of extra. Pass in NULL (or None in python) for standard waveform. Set "sideband",m to get a single sideband (m=-2..2) */
        const LALSimInspiralSpinOrder spinO,   /**< twice PN order of spin effects */
        const INT4 phaseO,                     /**< twice PN phase order */
        const INT4 amplitudeO                  /**< twice PN amplitude order */
        )
{
    /* external: SI; internal: solar masses */
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 eta = m1 * m2 / (m * m);
    const REAL8 piM = LAL_PI * m_sec;
    const REAL8 vISCO = 1. / sqrt(6.);
    const REAL8 fISCO = vISCO * vISCO * vISCO / piM;
    REAL8 shft, amp0, f_max;
    size_t i, n, iStart;
    COMPLEX16 *data_plus = NULL;
    COMPLEX16 *data_cross = NULL;
    LIGOTimeGPS tC = {0, 0};

    /* If f_ref = 0, use f_ref = f_low for everything except the phase offset */
    const REAL8 v_ref = f_ref > 0. ? cbrt(piM*f_ref) : cbrt(piM*fStart);

    REAL8 alpha, alpha_ref;
    COMPLEX16 prec_plus, prec_cross, phasing_fac;
    bool enable_precession = true; /* Handle the non-spinning case separately */
    int mm;

    LALSimInspiralSF2Orientation orientation;
    XLALSimInspiralSF2CalculateOrientation(&orientation, m1, m2, v_ref, lnhatx, lnhaty, lnhatz, s1x, s1y, s1z);

    LALSimInspiralSF2Coeffs coeffs;
    XLALSimInspiralSF2CalculateCoeffs(&coeffs, m1, m2, orientation.chi, orientation.kappa);
    enable_precession = orientation.chi != 0. && orientation.kappa != 1.;

    alpha_ref = enable_precession ? XLALSimInspiralSF2Alpha(v_ref, coeffs) - orientation.alpha0 : 0.;

    COMPLEX16 SBplus[5]; /* complex sideband factors for plus pol, mm=2 is first entry */
    COMPLEX16 SBcross[5]; /* complex sideband factors for cross pol, mm=2 is first entry */
    REAL8 emission[5]; /* emission factor for each sideband */
    if ( !XLALSimInspiralTestGRParamExists(moreParams, "sideband") )
    {
        for(mm = -2; mm <= 2; mm++)
        {
            SBplus[2-mm] = XLALSimInspiralSF2Polarization(orientation.thetaJ, orientation.psiJ, mm);
            SBcross[2-mm] = XLALSimInspiralSF2Polarization(orientation.thetaJ, orientation.psiJ+LAL_PI/4., mm);
        }
    }
    else
    {
        memset(SBplus, 0, 5 * sizeof(COMPLEX16));
        memset(SBcross, 0, 5 * sizeof(COMPLEX16));
        mm = (int) XLALSimInspiralGetTestGRParam(moreParams, "sideband");
        SBplus[2-mm] = XLALSimInspiralSF2Polarization(orientation.thetaJ, orientation.psiJ, mm);
        SBcross[2-mm] = XLALSimInspiralSF2Polarization(orientation.thetaJ, orientation.psiJ+LAL_PI/4., mm);
    }

    const REAL8 chi1L = orientation.chi*orientation.kappa;
    const REAL8 chi1sq = orientation.chi*orientation.chi;
    /* FIXME: Cannot yet set QM constant in ChooseFDWaveform interface */
    const REAL8 quadparam1 = 1.;
    /* phasing coefficients */
    PNPhasingSeries pfa;
    XLALSimInspiralPNPhasing_F2(&pfa, m1, m2, chi1L, 0., chi1sq, 0., 0., quadparam1, 0., spinO);

    REAL8 pfaN = 0.;
    REAL8 pfa2 = 0.; REAL8 pfa3 = 0.; REAL8 pfa4 = 0.;
    REAL8 pfa5 = 0.; REAL8 pfl5 = 0.;
    REAL8 pfa6 = 0.; REAL8 pfl6 = 0.;
    REAL8 pfa7 = 0.; REAL8 pfa8 = 0.;

    switch (phaseO)
    {
        case -1:
        case 7:
            pfa7 = pfa.v[7];
        case 6:
            pfa6 = pfa.v[6];
            pfl6 = pfa.vlogv[6];
        case 5:
            pfa5 = pfa.v[5];
            pfl5 = pfa.vlogv[5];
        case 4:
            pfa4 = pfa.v[4];
        case 3:
            pfa3 = pfa.v[3];
        case 2:
            pfa2 = pfa.v[2];
        case 0:
            pfaN = pfa.v[0];
            break;
        default:
            XLAL_ERROR(XLAL_ETYPE, "Invalid phase PN order %s", phaseO);
    }

    /* Add the zeta factor to the phasing, since it looks like a pN series.
     * This is the third Euler angle after alpha and beta.
     */
    if (enable_precession)
    {
        pfa2 += 2.*coeffs.prec_fac*coeffs.zc[0];
        pfa3 += 2.*coeffs.prec_fac*coeffs.zc[1];
        pfa4 += 2.*coeffs.prec_fac*coeffs.zc[2];
        pfl5 += 2.*coeffs.prec_fac*coeffs.zc[3];
        pfa6 += 2.*coeffs.prec_fac*coeffs.zc[4];
        pfa7 += 2.*coeffs.prec_fac*coeffs.zc[5];
        pfa8 += 2.*coeffs.prec_fac*coeffs.zc[6];
    }

    /* Validate expansion order arguments.
     * This must be done here instead of in the OpenMP parallel loop
     * because when OpenMP parallelization is turned on, early exits
     * from loops (via return or break statements) are not permitted.
     */

    /* Validate amplitude PN order. */
    if (amplitudeO != 0) { XLAL_ERROR(XLAL_ETYPE, "Invalid amplitude PN order %s", amplitudeO); }

    /* energy coefficients - not currently used, but could for MECO
    const REAL8 dETaN = 2. * XLALSimInspiralPNEnergy_0PNCoeff(eta);
    const REAL8 dETa1 = 2. * XLALSimInspiralPNEnergy_2PNCoeff(eta);
    const REAL8 dETa2 = 3. * XLALSimInspiralPNEnergy_4PNCoeff(eta);
    const REAL8 dETa3 = 4. * XLALSimInspiralPNEnergy_6PNCoeff(eta);
    */

    COMPLEX16FrequencySeries *hplus;
    COMPLEX16FrequencySeries *hcross;

    /* Perform some initial checks */
    if (!hplus_out) XLAL_ERROR(XLAL_EFAULT);
    if (!hcross_out) XLAL_ERROR(XLAL_EFAULT);
    if (*hplus_out) XLAL_ERROR(XLAL_EFAULT);
    if (*hcross_out) XLAL_ERROR(XLAL_EFAULT);
    if (m1_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
    if (f_ref < 0) XLAL_ERROR(XLAL_EDOM);
    if (r <= 0) XLAL_ERROR(XLAL_EDOM);

    /* allocate htilde */
    if ( fEnd == 0. ) // End at ISCO
        f_max = fISCO;
    else // End at user-specified freq.
        f_max = fEnd;
    n = (size_t) (f_max / deltaF + 1);
    XLALGPSAdd(&tC, -1 / deltaF);  /* coalesce at t=0 */
    /* Allocate hplus and hcross */
    hplus = XLALCreateCOMPLEX16FrequencySeries("hplus: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, n);
    if (!hplus) XLAL_ERROR(XLAL_EFUNC);
    memset(hplus->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitDivide(&hplus->sampleUnits, &hplus->sampleUnits, &lalSecondUnit);
    hcross = XLALCreateCOMPLEX16FrequencySeries("hcross: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, n);
    if (!hcross) XLAL_ERROR(XLAL_EFUNC);
    memset(hcross->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitDivide(&hcross->sampleUnits, &hcross->sampleUnits, &lalSecondUnit);

    /* extrinsic parameters */
    amp0 = -4. * m1 * m2 / r * LAL_MRSUN_SI * LAL_MTSUN_SI * sqrt(LAL_PI/12.L);
    amp0 *= sqrt(-2. * XLALSimInspiralPNEnergy_0PNCoeff(eta)/XLALSimInspiralPNFlux_0PNCoeff(eta));
    shft = LAL_TWOPI * (tC.gpsSeconds + 1e-9 * tC.gpsNanoSeconds);

    /* Fill with non-zero vals from fStart to f_max */
    iStart = (size_t) ceil(fStart / deltaF);
    data_plus = hplus->data->data;
    data_cross = hcross->data->data;

    /* Compute the SPA phase at the reference point */
    REAL8 ref_phasing = 0.;
    if (f_ref > 0.)
    {
        const REAL8 logvref = log(v_ref);
        const REAL8 v2ref = v_ref * v_ref;
        const REAL8 v3ref = v_ref * v2ref;
        const REAL8 v4ref = v_ref * v3ref;
        const REAL8 v5ref = v_ref * v4ref;
        ref_phasing = (pfaN + pfa2 * v2ref + pfa3 * v3ref + pfa4 * v4ref) / v5ref + (pfa5 + pfl5 * logvref) + (pfa6 + pfl6 * logvref) * v_ref + pfa7 * v2ref + pfa8 * v3ref;
    } /* end of if (f_ref > 0.) */

    #pragma omp parallel for
    for (i = iStart; i < n; i++) {
        const REAL8 f = i * deltaF;
        const REAL8 v = cbrt(piM*f);
        const REAL8 logv = log(v);
        const REAL8 v2 = v * v;
        const REAL8 v3 = v * v2;
        const REAL8 v4 = v * v3;
        const REAL8 v5 = v * v4;
        REAL8 phasing = (pfaN + pfa2 * v2 + pfa3 * v3 + pfa4 * v4) / v5 + (pfa5 + pfl5 * logv) + (pfa6 + pfl6 * logv) * v + pfa7 * v2 + pfa8 * v3;
        COMPLEX16 amp = amp0 / (v3 * sqrt(v));

        alpha = enable_precession ? XLALSimInspiralSF2Alpha(v, coeffs) - alpha_ref : 0.;

        COMPLEX16 u = cos(alpha) + 1.0j*sin(alpha);
        XLALSimInspiralSF2Emission(emission, v, coeffs);
        prec_plus = SBplus[0]*emission[0]*u*u
                    + SBplus[1]*emission[1]*u
                    + SBplus[2]*emission[2]
                    + SBplus[3]*emission[3]/u
                    + SBplus[4]*emission[4]/u/u;
        prec_cross= SBcross[0]*emission[0]*u*u
                    + SBcross[1]*emission[1]*u
                    + SBcross[2]*emission[2]
                    + SBcross[3]*emission[3]/u
                    + SBcross[4]*emission[4]/u/u;
        // Note the factor of 2 b/c phi_ref is orbital phase
        phasing += shft * f - 2.*phi_ref - ref_phasing - LAL_PI_4;
        phasing_fac = cos(phasing) - 1.0j*sin(phasing);
        data_plus[i] = amp * prec_plus * phasing_fac;
        data_cross[i] = amp * prec_cross * phasing_fac;
    }

    *hplus_out = hplus;
    *hcross_out = hcross;
    return XLAL_SUCCESS;
}

