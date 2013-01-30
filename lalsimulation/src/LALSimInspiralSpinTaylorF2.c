/*
 *  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
 *  Copyright (C) 2012 Leo Singer
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
#include <stdio.h>


typedef struct tagLALSimInspiralSF2Orientation
{
    REAL8 thetaJ, psiJ;
    REAL8 chi, kappa, alpha0;
    REAL8 v_ref;
} LALSimInspiralSF2Orientation;

typedef struct tagLALSimInspiralSF2Coeffs
{
    REAL8 m1, m2, mtot, eta;
    REAL8 kappa, gamma0;
    REAL8 prec_fac;
    REAL8 pn_beta, pn_sigma, pn_gamma;
    REAL8 dtdv2, dtdv3, dtdv4, dtdv5;
} LALSimInspiralSF2Coeffs;

// Prototypes 
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

COMPLEX16 XLALSimInspiralSF2Emission(
    REAL8 beta, int mm);


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
    const REAL8 psiJ = atan2(Jy0, -Jx0); /* FIXME: check that Jy0 and Jx0 are not both 0 */
    orientation->psiJ = psiJ;

    /* Rotate Lnhat back to frame where J is along z, to figure out initial alpha */
    const REAL8 rotLx = lnhatx*cos(thetaJ)*cos(psiJ) - lnhaty*cos(thetaJ)*sin(psiJ) + lnhatz*sin(thetaJ);
    const REAL8 rotLy = lnhatx*sin(psiJ) + lnhaty*cos(psiJ);
    orientation->alpha0 = atan2(rotLy, rotLx); /* FIXME: check that rotLy and rotLx are not both 0 */
}

void XLALSimInspiralSF2CalculateCoeffs(
    LALSimInspiralSF2Coeffs *coeffs,
    REAL8 m1, REAL8 m2, REAL8 chi, REAL8 kappa )
{
    coeffs->m1 = m1;
    coeffs->m2 = m2;
    const REAL8 mtot = m1+m2;
    coeffs->mtot = mtot;
    const REAL8 eta = m1*m2/(mtot*mtot);
    coeffs->eta = eta;
    coeffs->kappa = kappa;
    coeffs->gamma0 = m1*chi/m2;

    coeffs->pn_beta = (113.*m1/(12.*mtot) - 19.*eta/6.)*chi*kappa;
    coeffs->pn_sigma = 0.; // Add these in later
    coeffs->pn_gamma = 0.;    

    coeffs->prec_fac = 5.*(4. + 3.*m2/m1)/64.;
    coeffs->dtdv2 = 743./336. + 11.*eta/4.;
    coeffs->dtdv3 = -4.*LAL_PI + coeffs->pn_beta;
    coeffs->dtdv4 = 3058673./1016064. + 5429.*eta/1008. + 617.*eta*eta/144.;
    coeffs->dtdv5 = (-7729./672.+13.*eta/8.)*LAL_PI;
}

REAL8 XLALSimInspiralSF2Alpha(
    REAL8 v, LALSimInspiralSF2Coeffs coeffs)
{
    const REAL8 gamma0 = coeffs.gamma0;
    const REAL8 gam = gamma0*v;
    const REAL8 kappa = coeffs.kappa;

    const REAL8 sqrtfac = sqrt(1. + 2.*kappa*gam + gam*gam);
    const REAL8 logv = log(v);
    const REAL8 logfac1 = log(1. + kappa*gam + sqrtfac);
    const REAL8 logfac2 = log(kappa + gam + sqrtfac);

    const REAL8 prec_fac = coeffs.prec_fac;
    const REAL8 dtdv2 = coeffs.dtdv2;
    const REAL8 dtdv3 = coeffs.dtdv3;
    const REAL8 dtdv4 = coeffs.dtdv4;
    const REAL8 dtdv5 = coeffs.dtdv5;

    return prec_fac*(logfac2*(dtdv2*gamma0*pow(1,2) + dtdv3*kappa*pow(1,3) - (dtdv5*kappa*pow(1,5)*pow(gamma0,-2))/2. + (dtdv4*pow(1,4)*pow(gamma0,-1))/2. - (dtdv4*pow(1,4)*pow(gamma0,-1)*pow(kappa,2))/2. + (dtdv5*pow(1,5)*pow(gamma0,-2)*pow(kappa,3))/2.) + logfac1*(-(dtdv2*gamma0*kappa*pow(1,2)) - dtdv3*pow(1,3) + (kappa*pow(gamma0,3))/2. - (pow(gamma0,3)*pow(kappa,3))/2.) + logv*(dtdv2*gamma0*kappa*pow(1,2) + dtdv3*pow(1,3) - (kappa*pow(gamma0,3))/2. + (pow(gamma0,3)*pow(kappa,3))/2.) + sqrtfac*(dtdv3*pow(1,3) + (dtdv4*v*pow(1,4))/2. + (dtdv5*pow(1,5)*pow(gamma0,-2))/3. + (dtdv4*kappa*pow(1,4)*pow(gamma0,-1))/2. + (dtdv5*kappa*v*pow(1,5)*pow(gamma0,-1))/6. - (dtdv5*pow(1,5)*pow(gamma0,-2)*pow(kappa,2))/2. - pow(v,-3)/3. - (gamma0*kappa*pow(v,-2))/6. - dtdv2*pow(1,2)*pow(v,-1) - (pow(gamma0,2)*pow(v,-1))/3. + (pow(gamma0,2)*pow(kappa,2)*pow(v,-1))/2. + (dtdv5*pow(1,5)*pow(v,2))/3.));
}

REAL8 XLALSimInspiralSF2Zeta(
    REAL8 v, LALSimInspiralSF2Coeffs coeffs)
{
    const REAL8 gamma0 = coeffs.gamma0;
    const REAL8 gam = gamma0*v;
    const REAL8 kappa = coeffs.kappa;

    const REAL8 sqrtfac = sqrt(1. + 2.*kappa*gam + gam*gam);
    const REAL8 logv = log(v);
    const REAL8 logfac1 = log(1. + kappa*gam + sqrtfac);
    const REAL8 logfac2 = log(kappa + gam + sqrtfac);

    const REAL8 prec_fac = coeffs.prec_fac;
    const REAL8 dtdv2 = coeffs.dtdv2;
    const REAL8 dtdv3 = coeffs.dtdv3;
    const REAL8 dtdv4 = coeffs.dtdv4;
    const REAL8 dtdv5 = coeffs.dtdv5;

    return prec_fac*(dtdv3*gamma0*kappa*v*pow(1,3) + dtdv4*v*pow(1,4) + logfac2*(-(dtdv2*gamma0*pow(1,2)) - dtdv3*kappa*pow(1,3) + (dtdv5*kappa*pow(1,5)*pow(gamma0,-2))/2. - (dtdv4*pow(1,4)*pow(gamma0,-1))/2. + (dtdv4*pow(1,4)*pow(gamma0,-1)*pow(kappa,2))/2. - (dtdv5*pow(1,5)*pow(gamma0,-2)*pow(kappa,3))/2.) + logv*((kappa*pow(gamma0,3))/2. - (pow(gamma0,3)*pow(kappa,3))/2.) + logfac1*(dtdv2*gamma0*kappa*pow(1,2) + dtdv3*pow(1,3) - (kappa*pow(gamma0,3))/2. + (pow(gamma0,3)*pow(kappa,3))/2.) - pow(v,-3)/3. - (gamma0*kappa*pow(v,-2))/2. - dtdv2*pow(1,2)*pow(v,-1) + (dtdv4*gamma0*kappa*pow(1,4)*pow(v,2))/2. + (dtdv5*pow(1,5)*pow(v,2))/2. + sqrtfac*(-(dtdv3*pow(1,3)) - (dtdv4*v*pow(1,4))/2. - (dtdv5*pow(1,5)*pow(gamma0,-2))/3. - (dtdv4*kappa*pow(1,4)*pow(gamma0,-1))/2. - (dtdv5*kappa*v*pow(1,5)*pow(gamma0,-1))/6. + (dtdv5*pow(1,5)*pow(gamma0,-2)*pow(kappa,2))/2. + pow(v,-3)/3. + (gamma0*kappa*pow(v,-2))/6. + dtdv2*pow(1,2)*pow(v,-1) + (pow(gamma0,2)*pow(v,-1))/3. - (pow(gamma0,2)*pow(kappa,2)*pow(v,-1))/2. - (dtdv5*pow(1,5)*pow(v,2))/3.) + (dtdv5*gamma0*kappa*pow(1,5)*pow(v,3))/3.);
}

REAL8 XLALSimInspiralSF2Beta(
    REAL8 v, LALSimInspiralSF2Coeffs coeffs)
{
    const REAL8 kappa = coeffs.kappa;
    const REAL8 gamma0 = coeffs.gamma0;

    return acos((1. + kappa*gamma0*v)/sqrt(1. + 2.*kappa*gamma0*v + gamma0*gamma0*v*v));
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

COMPLEX16 XLALSimInspiralSF2Emission(
    REAL8 beta, int mm)
{
    return pow(cos(beta/2.), 2+mm) * pow(sin(beta/2.), 2-mm);
}


/**
 * Find the least nonnegative integer power of 2 that is
 * greater than or equal to n.  Inspired by similar routine
 * in gstlal.
 */
static size_t CeilPow2(double n) {
    double signif;
    int exponent;
    signif = frexp(n, &exponent);
    if (signif < 0)
        return 1;
    if (signif == 0.5)
        exponent -= 1;
    return ((size_t) 1) << exponent;
}

/* Calculate the spin corrections for TaylorF2 
        reference -> <http://arxiv.org/pdf/0810.5336v3.pdf>
*/

typedef struct {
    REAL8 beta;
    REAL8 sigma;
    REAL8 gamma;
}sf2_spin_corr;

sf2_spin_corr sf2_spin_corrections(
        const REAL8 m1,               /**< mass of companion 1 (solar masses) */
        const REAL8 m2,               /**< mass of companion 2 (solar masses) */
		const REAL8 S1z,               /**< z component of the spin of companion 1 */
		const REAL8 S2z ) ;            /**< z component of the spin of companion 2  */


sf2_spin_corr sf2_spin_corrections(
        const REAL8 m1,               /**< mass of companion 1 (solar masses) */
        const REAL8 m2,               /**< mass of companion 2 (solar masses) */
		const REAL8 S1z,               /**< z component of the spin of companion 1 */
		const REAL8 S2z )             /**< z component of the spin of companion 2  */
{
    sf2_spin_corr spin_corrections;

    REAL8 M = m1 + m2;
    REAL8 v = m1 * m2 / (M * M);
    REAL8 d = (m1 - m2) / (m1 + m2);
    REAL8 xs = .5 * (S1z + S2z);
    REAL8 xa = .5 * (S1z - S2z);

    REAL8 sf2_beta = (113.L/12.L- 19.L/3.L * v) * xs + 113.L/12.L * d * xa;
    
    REAL8 sf2_sigma = v * (721.L/48.L * (xs*xs - xa*xa)-247.L/48.L*(xs*xs - xa*xa));
    sf2_sigma += (1-2*v)* (719/96.0 * (xs*xs + xa*xa) - 233.L/96.L * (xs*xs + xa*xa));
    sf2_sigma += d * (719/48.0 * xs*xa - 233.L/48.L * xs * xa);
    
    REAL8 sf2_gamma = (732985.L/2268.L - 24260.L/81.L * v - 340.L/9.L * v * v ) * xs;
    sf2_gamma += (732985.L/2268.L +140.L/9.0L * v) * xa * d;

    spin_corrections.beta = sf2_beta;
    spin_corrections.sigma = sf2_sigma;
    spin_corrections.gamma = sf2_gamma;
    return spin_corrections;
}
/**
 * Computes the stationary phase approximation to the Fourier transform of
 * a chirp waveform with phase given by Eq.\eqref{eq_InspiralFourierPhase_f2}
 * and amplitude given by expanding \f$1/\sqrt{\dot{F}}\f$. If the PN order is
 * set to -1, then the highest implemented order is used.
 * \author B.S. Sathyaprakash
 */
int XLALSimInspiralSpinTaylorF2(
	COMPLEX16FrequencySeries **htilde_out, /**< frequency-domain waveform */
	REAL8 psi,                      /**< desired polarization */
	REAL8 phic,                     /**< coalescence GW phase */
	REAL8 deltaF,                   /**< sampling frequency (Hz) */
	REAL8 m1_SI,                    /**< mass of companion 1 (kg) */
	REAL8 m2_SI,                    /**< mass of companion 2 (kg) */
	REAL8 fStart,                   /**< start GW frequency (Hz) */
	REAL8 r,                        /**< distance of source (m) */
	REAL8 s1x,                      /**< initial value of S1x */
	REAL8 s1y,                      /**< initial value of S1y */
	REAL8 s1z,                      /**< initial value of S1z */
	REAL8 lnhatx,                   /**< initial value of LNhatx */
	REAL8 lnhaty,                   /**< initial value of LNhaty */
	REAL8 lnhatz,                   /**< initial value of LNhatz */
	int phaseO,                     /**< twice PN phase order */
	int amplitudeO                  /**< twice PN amplitude order */
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
    const REAL8 v0 = cbrt(piM * fStart);
    REAL8 shft, phi0, amp0, f_max;
    size_t i, n, iStart, iISCO;
    int mm;
    COMPLEX16 *data = NULL;
    LIGOTimeGPS tC = {0, 0};

    REAL8 alpha, alpha_ref, zeta, zeta_ref, beta;
    COMPLEX16 prec_fac;

    LALSimInspiralSF2Orientation orientation;
    XLALSimInspiralSF2CalculateOrientation(&orientation, m1, m2, v0, lnhatx, lnhaty, lnhatz, s1x, s1y, s1z);

    orientation.psiJ = orientation.psiJ + psi;

    /* fprintf(stdout, "thetaJ = %f\n", orientation.thetaJ * 180./LAL_PI);
    fprintf(stdout, "psiJ = %f\n", orientation.psiJ * 180./LAL_PI);
    fprintf(stdout, "alpha0 = %f\n", orientation.alpha0 * 180./LAL_PI);
    fprintf(stdout, "chi = %f\n", orientation.chi);
    fprintf(stdout, "kappa = %f\n", orientation.kappa);
    fprintf(stdout, "v0 = %f\n", v0); */


    LALSimInspiralSF2Coeffs coeffs;
    XLALSimInspiralSF2CalculateCoeffs(&coeffs, m1, m2, orientation.chi, orientation.kappa);

    alpha_ref = XLALSimInspiralSF2Alpha(v0, coeffs) - orientation.alpha0;
    zeta_ref = XLALSimInspiralSF2Zeta(v0, coeffs);
 
    COMPLEX16 SBfac[5]; /* complex sideband factors, mm=2 is first entry */
    for(mm = -2; mm <= 2; mm++)
    {
        SBfac[2-mm] = XLALSimInspiralSF2Polarization(orientation.thetaJ, orientation.psiJ, mm);
        // fprintf(stdout, "m = %i, %f, %f\n", mm, creal(SBfac[2-mm]), cimag(SBfac[2-mm]));
    }

    const REAL8 pn_beta = coeffs.pn_beta;
    const REAL8 pn_sigma = 0.;
    const REAL8 pn_gamma = 0.; /* FIXME: will add higher PN later */

    /* phasing coefficients */
    const REAL8 pfaN = 3.L/(128.L * eta);
    const REAL8 pfa2 = 5.L*(743.L/84.L + 11.L * eta)/9.L;
    const REAL8 pfa3 = -16.L*LAL_PI + 4.L*pn_beta;
    const REAL8 pfa4 = 5.L*(3058.673L/7.056L + 5429.L/7.L * eta
                     + 617.L * eta*eta)/72.L - 10.L*pn_sigma;
    const REAL8 pfa5 = 5.L/9.L * (7729.L/84.L - 13.L * eta) * LAL_PI - pn_gamma;
    const REAL8 pfl5 = 5.L/3.L * (7729.L/84.L - 13.L * eta) * LAL_PI - pn_gamma *3.0;
    const REAL8 pfa6 = (11583.231236531L/4.694215680L - 640.L/3.L * LAL_PI * LAL_PI - 6848.L/21.L*LAL_GAMMA)
                     + eta * (-15335.597827L/3.048192L + 2255./12. * LAL_PI * LAL_PI - 1760./3.*theta +12320./9.*lambda)
                     + eta*eta * 76055.L/1728.L
                     - eta*eta*eta*  127825.L/1296.L ;
    const REAL8 pfl6 = -6848.L/21.L;
    const REAL8 pfa7 = LAL_PI * 5.L/756.L * ( 15419335.L/336.L + 75703.L/2.L * eta - 14809.L * eta*eta);

    /* flux coefficients */
    const REAL8 FTaN = XLALSimInspiralTaylorT1Flux_0PNCoeff(eta);
    const REAL8 FTa2 = XLALSimInspiralTaylorT1Flux_2PNCoeff(eta);
    const REAL8 FTa3 = XLALSimInspiralTaylorT1Flux_3PNCoeff(eta);
    const REAL8 FTa4 = XLALSimInspiralTaylorT1Flux_4PNCoeff(eta);
    const REAL8 FTa5 = XLALSimInspiralTaylorT1Flux_5PNCoeff(eta);
    const REAL8 FTl6 = XLALSimInspiralTaylorT1Flux_6PNLogCoeff(eta);
    const REAL8 FTa6 = XLALSimInspiralTaylorT1Flux_6PNCoeff(eta);
    const REAL8 FTa7 = XLALSimInspiralTaylorT1Flux_7PNCoeff(eta);

    /* energy coefficients */
    const REAL8 dETaN = 2. * XLALSimInspiralPNEnergy_0PNCoeff(eta);
    const REAL8 dETa1 = 2. * XLALSimInspiralPNEnergy_2PNCoeff(eta);
    const REAL8 dETa2 = 3. * XLALSimInspiralPNEnergy_4PNCoeff(eta);
    const REAL8 dETa3 = 4. * XLALSimInspiralPNEnergy_6PNCoeff(eta);

    COMPLEX16FrequencySeries *htilde;

    /* Perform some initial checks */
    if (!htilde_out) XLAL_ERROR(XLAL_EFAULT);
    if (*htilde_out) XLAL_ERROR(XLAL_EFAULT);
    if (phic < 0) XLAL_ERROR(XLAL_EDOM);
    if (m1_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
    if (r <= 0) XLAL_ERROR(XLAL_EDOM);

    /* allocate htilde */
    f_max = CeilPow2(fISCO);
    n = f_max / deltaF + 1;
    XLALGPSAdd(&tC, -1 / deltaF);  /* coalesce at t=0 */
    htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, n);
    if (!htilde) XLAL_ERROR(XLAL_EFUNC);
    memset(htilde->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitDivide(&htilde->sampleUnits, &htilde->sampleUnits, &lalSecondUnit);

    /* extrinsic parameters */
    phi0 = phic;
    amp0 = 4. * m1 * m2 / r * LAL_MRSUN_SI * LAL_MTSUN_SI * sqrt(LAL_PI/12.L); /* Why was there a factor of deltaF in the lalinspiral version? */
    shft = -LAL_TWOPI * (tC.gpsSeconds + 1e-9 * tC.gpsNanoSeconds);

    iStart = (size_t) ceil(fStart / deltaF);
    iISCO = (size_t) (fISCO / deltaF);
    iISCO = (iISCO < n) ? iISCO : n;  /* overflow protection; should we warn? */
    data = htilde->data->data;
    for (i = iStart; i < iISCO; i++) {
        const REAL8 f = i * deltaF;
        const REAL8 v = cbrt(piM*f);
        const REAL8 v2 = v * v;
        const REAL8 v3 = v * v2;
        const REAL8 v4 = v * v3;
        const REAL8 v5 = v * v4;
        const REAL8 v6 = v * v5;
        const REAL8 v7 = v * v6;
        const REAL8 v8 = v * v7;
        const REAL8 v9 = v * v8;
        const REAL8 v10 = v * v9;
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
                phasing += (pfa6 + pfl6 * log(4.*v) ) * v6;
            case 5:
                phasing += (pfa5 + pfl5 * log(v/v0)) * v5;
            case 4:
                phasing += pfa4 * v4;
            case 3:
                phasing += pfa3 * v3;
            case 2:
                phasing += pfa2 * v2;
            case 0:
                phasing += 1.;
                break;
            default:
                XLALDestroyCOMPLEX16FrequencySeries(htilde);
                XLAL_ERROR(XLAL_ETYPE);
        }
        switch (amplitudeO)
        {
            case -1:
            case 7:
                flux += FTa7 * v7;
            case 6:
                flux += (FTa6 + FTl6*log(16.*v2)) * v6;
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
                break;
            default:
                XLALDestroyCOMPLEX16FrequencySeries(htilde);
                XLAL_ERROR(XLAL_ETYPE);
        }
        phasing *= pfaN / v5;
        flux *= FTaN * v10;
        dEnergy *= dETaN * v;

        alpha = XLALSimInspiralSF2Alpha(v, coeffs) - alpha_ref;
        beta = XLALSimInspiralSF2Beta(v, coeffs);
        zeta = XLALSimInspiralSF2Zeta(v, coeffs) - zeta_ref;

        prec_fac = 0.j;
        for(mm = -2; mm <= 2; mm++)
        {
            prec_fac += 
                XLALSimInspiralSF2Emission(beta, mm)
                * SBfac[2-mm]
                * ( cos( (mm-2.) * alpha) + sin( (mm-2.) * alpha)*1.0j );
        }

        phasing += shft * f + phi0;
        phasing += 2.*zeta;
        amp = amp0 * sqrt(-dEnergy/flux) * v;
        data[i] = ((COMPLEX16)amp)*prec_fac*(cos(phasing + LAL_PI_4) - sin(phasing + LAL_PI_4) * 1.0j);
    }

    *htilde_out = htilde;
    return XLAL_SUCCESS;
}
