/*
*  Copyright (C) 2011 P. Ajith
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

#include <lal/AVFactories.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <lal/Units.h>
#include <lal/XLALError.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

/* Computed in Mathematica with N[\[Pi]^2, 35] */
#define Pi_p2 9.8696044010893586188344909998761511
#define Pi_p3 31.006276680299820175476315067101395
#define Pi_p4 97.409091034002437236440332688705111
#define Pi_p5 306.01968478528145326274131004343561
#define Pi_p4by3 4.6011511144704899609836928301589280
#define Pi_p5by3 6.7388085956981412852178979102191822
#define Pi_p10by3 45.411541289455155012188835024560824
#define Pi_p7by3 14.454942539276981063196177096150014
#define Pi_p1by3 1.4645918875615232630201425272637904
#define Pi_p2by3 2.1450293971110256000774441009412356
#define Pi_11by3 66.509374974201175524807943232081902
#define tenbyPi_p2by3 2.1638812222639821850478371484217674
#define twoPi_p2by3 3.4050219214767547368607032470366207

/**
 * Function to compute the metric elements using waveform derivatives
 */
static REAL8 MetricCoeffs(REAL8Vector *A, REAL8Vector *dPsii, REAL8Vector *dPsij,
        REAL8Vector *dAi, REAL8Vector*dAj, REAL8Vector *Sh, REAL8 hSqr, REAL8 df) {
    size_t k = A->length;
    REAL8 gij = 0.;
    for (;k--;) {
        gij += (A->data[k]*A->data[k]*dPsii->data[k]*dPsij->data[k]
                + dAi->data[k]*dAj->data[k])/Sh->data[k];
    }
    return 4.*df*gij/(2.*hSqr);
}

/**
 * Frequency domain amplitude of the TaylorF2 Reduced Spin waveforms
 */
static REAL8 XLALSimInspiralTaylorF2RedSpinAofF(
    const REAL8 mc_msun,   /**< chirp mass (solar masses) */
    const REAL8 eta,   /**< symmetric mass ratio  */
    const REAL8 chi,   /**< reduced-spin parameter */
    const REAL8 f   /**< Fourier frequency (Hz) */
) {
    const REAL8 mc = mc_msun * LAL_MTSUN_SI;
    const REAL8 etap2 = eta*eta;
    const REAL8 etap3 = etap2*eta;
    const REAL8 chip2 = chi*chi;
    const REAL8 mcp2  = mc*mc;

    REAL8 eta_three_fifths = pow(eta, 0.6);
    REAL8 f_mc_eta = (f*mc)/eta_three_fifths;
    REAL8 f_mc_eta_one_third = cbrt(f_mc_eta);

    return (sqrt(0.8333333333333334)*pow(mc/eta_three_fifths,0.8333333333333334)*
        sqrt(eta)*(1 + (Pi_p4by3*
         f_mc_eta_one_third * f_mc_eta_one_third *
         (743 + 924*eta))/672. -
        (Pi_p10by3 *
         f_mc_eta * f_mc_eta * f_mc_eta_one_third*
         (5111593 + 8088752*eta + 151088*etap2))/2.709504e6 +
        (f*mc*LAL_PI*(-2*LAL_PI + (113*chi)/24.))/eta_three_fifths +
        (Pi_p5by3 *
         f_mc_eta * f_mc_eta / f_mc_eta_one_third *
         (12*LAL_PI*(-4757 + 4788*eta) -
           (113*(502429 - 591368*eta + 1680*etap2)*chi)/
            (-113 + 76*eta)))/16128. +
        (Pi_p5 * Pi_p1by3 *
         f_mc_eta * f_mc_eta_one_third *
         (7266251 + 9532152*eta + 9730224*etap2 +
           (3243530304*(-81 + 4*eta)*chip2)/
            pow(113 - 76*eta,2)))/8.128512e6 +
        (f*f*mcp2*Pi_p2*
         (-58.601030974347324 + (10*Pi_p2)/3. +
           (3526813753*eta)/2.7869184e7 -
           (451*Pi_p2*eta)/96. -
           (1041557*etap2)/258048. +
           (67999*etap3)/82944. +
           (856*(3*LAL_GAMMA + log((64*f*mc*LAL_PI)/eta_three_fifths)))/315.))/(eta_three_fifths * eta_three_fifths)))/
        (2.*pow(f,1.1666666666666667)*Pi_p4by3);
}

/**
 * Derivative of the amplitude with respect to \f$\chi\f$ (reduced-spin parameter)
 */
static REAL8 XLALSimInspiralTaylorF2RedSpinDerivAChi(
        const REAL8 mc_msun,   /**< chirp mass (solar masses) */
        const REAL8 eta,  /**< symmetric mass ratio  */
        const REAL8 chi,  /**< reduced-spin parameter */
        const REAL8 f     /**< Fourier frequency (Hz) */
) {
    const REAL8 mc = mc_msun * LAL_MTSUN_SI;
    const REAL8 etap2 = eta*eta;
    REAL8 etap35 = pow(eta, 0.6);

    return (113*sqrt(0.8333333333333334)*Pi_p1by3*
        pow(mc/etap35,0.8333333333333334)*sqrt(eta)*
        ((672*f*mc)/etap35 -
        (Pi_p1by3*Pi_p1by3*
         pow((f*mc)/etap35,1.6666666666666667)*
         (502429 - 591368*eta + 1680*etap2))/
        (-113 + 76*eta) + (113904*Pi_p1by3*
         pow((f*mc)/etap35,1.3333333333333333)*
         (-81 + 4*eta)*chi)/pow(113 - 76*eta,2)))/
    (32256.*pow(f,1.1666666666666667));
}

/**
 * Derivative of the amplitude with respect to \f$\eta\f$ (symm. mass ratio)
 */
static REAL8 XLALSimInspiralTaylorF2RedSpinDerivAEta(
    const REAL8 mc_msun,   /**< chirp mass (M_sun) */
    const REAL8 eta,  /**< symmetric mass ratio  */
    const REAL8 chi,  /**< reduced-spin parameter */
    const REAL8 f     /**< Fourier frequency (Hz) */
) {
    const REAL8 mc = mc_msun * LAL_MTSUN_SI;
    const REAL8 etap2 = eta*eta;
    const REAL8 etap3 = etap2*eta;
    const REAL8 etap4 = etap3*eta;
    const REAL8 etap5 = etap4*eta;
    const REAL8 chip2 = chi*chi;
    const REAL8 mcp2  = mc*mc;
    const REAL8 eta_fac = (-113 + 76 * eta) * (-113 + 76 * eta) * (-113 + 76 * eta);
    REAL8 eta_three_fifths = pow(eta, 0.6);
    REAL8 f_mc_eta_one_third = cbrt((f*mc)/eta_three_fifths);

    return (pow(mc/eta_three_fifths,0.8333333333333334)*
        (2235340800*f_mc_eta_one_third * f_mc_eta_one_third*
        eta_three_fifths*eta_three_fifths*eta_fac*(-743 + 1386*eta) +
        f*f*mcp2*Pi_p4by3*
        eta_fac *
        (257959397806029 - 36738273116160*LAL_GAMMA -
         95047630643350*eta - 12126223216800*etap2 +
         5541700903200*etap3 +
         7823692800*Pi_p2*(-1920 + 451*eta) -
         1940400*Pi_p4by3*
          f_mc_eta_one_third*
          (-5111593 - 2311072*eta + 64752*etap2)) +
        369600*f*mc*Pi_p1by3*eta_three_fifths*
        (12192768*LAL_PI*eta_fac +
         35962920*Pi_p5by3*
          f_mc_eta_one_third * f_mc_eta_one_third *
          eta_fac -
         28703808*eta_fac*chi -
         71190*Pi_p2by3*
          f_mc_eta_one_third * f_mc_eta_one_third *
          (-6415515901 + 12944580756*eta -
            10861276272*etap2 + 3401313728*etap3)*
          chi + Pi_p1by3*
          f_mc_eta_one_third *
          (34636018030144*etap3 -
            27532505500416*etap4 +
            6407002215936*etap5 -
            1442897*(-7266251 + 20575296*chip2) -
            21696*etap2*(-4887203 + 102257316*chip2) +
            76614*eta*(-320998087 + 907387488*chip2))) -
        12246091038720*f*f*mcp2*Pi_p4by3*
        eta_fac*log((64*f*mc*LAL_PI)/eta_three_fifths)))/
        (1.5021490176e12*sqrt(30)*pow(f,1.1666666666666667)*pow(eta,1.7)*
        eta_fac);
}

/**
 * Derivative of the amplitude with respect to the chirp mass
 */
static REAL8 XLALSimInspiralTaylorF2RedSpinDerivAMChirp(
    const REAL8 mc_msun,   /**< chirp mass (M_sun) */
    const REAL8 eta,  /**< symmetric mass ratio  */
    const REAL8 chi,  /**< reduced-spin parameter */
    const REAL8 f     /**< Fourier frequency (Hz) */
) {
    const REAL8 mc = mc_msun * LAL_MTSUN_SI;
    const REAL8 etap2 = eta*eta;
    const REAL8 etap3 = etap2*eta;
    const REAL8 etap4 = etap3*eta;
    const REAL8 chip2 = chi*chi;
    const REAL8 mcp2  = mc*mc;
    REAL8 eta_three_fifths = pow(eta, 0.6);

    return (sqrt(0.8333333333333334)*(5 -
      (17*f*f*mcp2*Pi_p2*Pi_p2*(-320 + 451*eta))/
       (96.*eta_three_fifths*eta_three_fifths) +
      (3*Pi_p2by3*
         pow((f*mc)/eta_three_fifths,0.6666666666666666)*
         (743 + 924*eta))/224. +
      (5*Pi_p3/Pi_p1by3*
         pow((f*mc)/eta_three_fifths,1.6666666666666667)*
         (-4757 + 4788*eta))/448. -
      (19*pow(LAL_PI,3.3333333333333335)*
         pow((f*mc)/eta_three_fifths,2.3333333333333335)*
         (5111593 + 8088752*eta + 151088*etap2))/2.709504e6 +
      (f*mc*Pi_p2*(-33047278387200*eta_three_fifths +
           f*mc*(-1471974996766431 + 208183547658240*LAL_GAMMA +
              3231619441873900*eta - 103072897342800*etap2 +
              20935314523200*etap3)))/
       (1.5021490176e12*eta_three_fifths*eta_three_fifths) +
      (1243*f*mc*LAL_PI*chi)/(24.*eta_three_fifths) -
      (565*Pi_p5by3*
         pow((f*mc)/eta_three_fifths,1.6666666666666667)*
         (502429 - 591368*eta + 1680*etap2)*chi)/
       (5376.*(-113 + 76*eta)) +
      (13*Pi_p4by3*
         pow((f*mc)/eta_three_fifths,1.3333333333333333)*
         (2490853280*etap2 - 112068617472*etap3 +
           56201773824*etap4 +
           1808*eta*(-1708561 + 7175952*chip2) -
           12769*(-7266251 + 20575296*chip2)))/
       (8.128512e6*pow(113 - 76*eta,2)) +
      (14552*f*f*mcp2*Pi_p2*
         log((64*f*mc*LAL_PI)/eta_three_fifths))/(315.*eta_three_fifths*eta_three_fifths)))/
    (12.*pow(f,1.1666666666666667)*Pi_p2by3*
    pow(mc/eta_three_fifths,0.16666666666666666)*pow(eta,0.1));
}

/**
 * Derivative of the phasae with respect to \f$\chi\f$ (reduced spin parameter)
 */
static REAL8 XLALSimInspiralTaylorF2RedSpinDerivPsiChi(
    const REAL8 mc_msun,   /**< chirp mass (M_sun) */
    const REAL8 eta,  /**< symmetric mass ratio  */
    const REAL8 chi,  /**< reduced-spin parameter */
    const REAL8 f     /**< Fourier frequency (Hz) */
) {
    const REAL8 mc = mc_msun * LAL_MTSUN_SI;
    const REAL8 etap2 = eta*eta;
    REAL8 eta_three_fifths = pow(eta, 0.6);

    return (113*((756*f*mc)/eta_three_fifths +
      (320355*Pi_p1by3*
         pow((f*mc)/eta_three_fifths,1.3333333333333333)*
         (-81 + 4*eta)*chi)/pow(113 - 76*eta,2) -
      (5*Pi_p2by3*
         pow((f*mc)/eta_three_fifths,1.6666666666666667)*
         (-146597 + 135856*eta + 17136*etap2)*
         (1 + log(f / 40.)))/
       (-113 + 76*eta)))/
    (96768.*Pi_p2by3*
    pow((f*mc)/eta_three_fifths,1.6666666666666667)*eta);
}

/**
 * Derivative of the phasae with respect to \f$\eta\f$ (symmetric mass ratio)
 */
static REAL8 XLALSimInspiralTaylorF2RedSpinDerivPsiEta(
    const REAL8 mc_msun,   /**< chirp mass (M_sun) */
    const REAL8 eta,  /**< symmetric mass ratio  */
    const REAL8 chi,  /**< reduced-spin parameter */
    const REAL8 f     /**< Fourier frequency (Hz) */
) {
    const REAL8 mc = mc_msun * LAL_MTSUN_SI;
    const REAL8 etap2 = eta*eta;
    const REAL8 etap3 = etap2*eta;
    const REAL8 etap4 = etap3*eta;
    const REAL8 etap5 = etap4*eta;
    const REAL8 chip2 = chi*chi;
    const REAL8 mcp2  = mc*mc;
    REAL8 eta_three_fifths = pow(eta, 0.6);
    const REAL8 eta_fac = -113 + 76 * eta;
    const REAL8 fmc_fac = cbrt(f * mc / eta_three_fifths);

    return (-77616000*f*mc*LAL_PI*eta_fac*
     (23187*LAL_PI*eta_fac*eta_fac -
       113*(16565461 - 22282744*eta + 12261424*etap2)*chi)*
     log(fmc_fac/
       (2.*cbrt(5)*
         cbrt(mc/eta_three_fifths))) +
    (31046400*fmc_fac*fmc_fac*
        eta_three_fifths * eta_three_fifths *eta_fac*eta_fac*eta_fac*(-743 + 1386*eta) -
       f*f*mcp2*LAL_PI * Pi_p1by3 *
        eta_fac * eta_fac * eta_fac *
        (33984313019673 - 4592284139520*LAL_GAMMA -
          12118079538950*eta - 413215941600*etap2 +
          2083465692000*etap3 +
          977961600*Pi_p2*(-3072 + 451*eta) +
          323400*LAL_PI*Pi_p1by3*
           fmc_fac*
           (15419335 + 3633744*eta + 2132496*etap2)) -
       18480*f*mc*Pi_p1by3*eta_three_fifths*
        (-6096384*LAL_PI*eta_fac * eta_fac * eta_fac +
          32461800*Pi_p2/Pi_p1by3*
           fmc_fac*fmc_fac*
           eta_fac * eta_fac * eta_fac +
          14351904*eta_fac * eta_fac * eta_fac*chi -
          158200*LAL_PI/Pi_p1by3*
           fmc_fac*fmc_fac*
           (-1871897093 + 3776925108*eta -
             3079029456*etap2 + 931868224*etap3)*
           chi - 5*Pi_p1by3*
           fmc_fac*
           (14990425815136*etap3 -
             12186233587584*etap4 +
             2866657264128*etap5 -
             1442897*(-3058673 + 5143824*chip2) +
             612912*eta*(-17749451 + 28355859*chip2) -
             2712*etap2*(-202619251 + 204514632*chip2)))
        + 1530761379840*f*f*mcp2*LAL_PI*Pi_p1by3*
        eta_fac * eta_fac * eta_fac*log((64*f*mc*LAL_PI)/eta_three_fifths))/
     (fmc_fac * fmc_fac *eta_three_fifths))/
    (5.007163392e11*f*mc*LAL_PI*etap2*eta_fac * eta_fac * eta_fac);
}

/**
 * Derivative of the phase with respect to the chirp mass
 */
static REAL8 XLALSimInspiralTaylorF2RedSpinDerivPsiMChirp(
    const REAL8 mc_msun,   /**< chirp mass (M_sun) */
    const REAL8 eta,  /**< symmetric mass ratio  */
    const REAL8 chi,  /**< reduced-spin parameter */
    const REAL8 f     /**< Fourier frequency (Hz) */
) {
    const REAL8 mc = mc_msun * LAL_MTSUN_SI;
    const REAL8 etap2 = eta*eta;
    const REAL8 etap3 = etap2*eta;
    const REAL8 etap4 = etap3*eta;
    const REAL8 chip2 = chi*chi;
    const REAL8 mcp2  = mc*mc;
    const REAL8 mcp3  = mcp2*mc;
    REAL8 eta_three_fifths = pow(eta, 0.6);

    return (cbrt((f*mc)/eta_three_fifths)*
    (-93139200*pow(113 - 76*eta,2)*eta_three_fifths * eta_three_fifths *
       (252 + Pi_p2by3*
          pow((f*mc)/eta_three_fifths,0.6666666666666666)*
          (743 + 924*eta)) +
      f*f*mcp2*LAL_PI*LAL_PI*pow(113 - 76*eta,2)*
       (10052469856691 - 1530761379840*LAL_GAMMA -
         24236159077900*eta + 206607970800*etap2 -
         462992376000*etap3 +
         1955923200*Pi_p2*(-512 + 451*eta) -
         184800*Pi_p4by3*
          cbrt((f*mc)/eta_three_fifths)*
          (-15419335 - 12718104*eta + 4975824*etap2)) +
      9240*f*mc*LAL_PI*eta_three_fifths*
       (16257024*LAL_PI*pow(113 - 76*eta,2) -
         38271744*pow(113 - 76*eta,2)*chi -
         5*cbrt(LAL_PI)*
          cbrt((f*mc)/eta_three_fifths)*
          (-20737091296*etap2 - 43167841920*etap3 +
            25146116352*etap4 +
            904*eta*(19183315 + 3587976*chip2) -
            12769*(-3058673 + 5143824*chip2))) -
      510253793280*f*f*mcp2*LAL_PI*LAL_PI*
       pow(113 - 76*eta,2)*log((64*f*mc*LAL_PI)/eta_three_fifths)))/
    (6.0085960704e11*f*f*mcp3*pow(LAL_PI,1.6666666666666667)*
    pow(113 - 76*eta,2)*eta);
}


/**
 * Compute the template-space metric of "reduced-spin" PN templates in
 * Mchirp-eta-chi parameter space.
 */
int XLALSimInspiralTaylorF2RedSpinMetricMChirpEtaChi(
    REAL8 *gamma00,  /**< template metric coeff. 00 in mChirp-eta-chi */
    REAL8 *gamma01,  /**< template metric coeff. 01/10 in mChirp-eta-chi */
    REAL8 *gamma02,  /**< template metric coeff. 02/20 in mChirp-eta-chi */
    REAL8 *gamma11,  /**< template metric coeff. 11 in mChirp-eta-chi */
    REAL8 *gamma12,  /**< template metric coeff. 12/21 in mChirp-eta-chi */
    REAL8 *gamma22,  /**< template metric coeff. 22 in mChirp-eta-chi */
    const REAL8 mc,     /**< chirp mass (in solar mass) */
    const REAL8 eta,    /**< symmetric mass ratio */
    const REAL8 chi,    /**< reduced-spin parameter */
    const REAL8 fLow,   /**< low-frequency cutoff (Hz) */
    const REAL8FrequencySeries *Sh  /**< PSD in strain per root Hertz */
) {
    REAL8Vector *A=NULL, *dAMc=NULL, *dAEta=NULL;
    REAL8Vector *dAChi=NULL, *dAT0=NULL, *dAPhi0=NULL, *dPsiMc=NULL;
    REAL8Vector *dPsiEta=NULL, *dPsiChi=NULL, *dPsiT0=NULL, *dPsiPhi0=NULL;

    const REAL8 df = Sh->deltaF;
    REAL8 hSqr = 0.;
    int s = 0;

    /* compute the Schwarzschild ISCO frequency */
    REAL8 fCut = 1. / (LAL_PI*mc*pow(eta, -0.6) * sqrt(6 * 6 * 6) * LAL_MTSUN_SI);

    /* create a view of the PSD between fLow and fCut */
    size_t nBins = (fCut - fLow) / df;
    size_t k = nBins;
    REAL8Vector Shdata = {nBins, Sh->data->data + (size_t) (fLow / df)}; /* copy the Vector, including its pointer to the actual data */
    /* drop low-frequency samples */
    Shdata.length = nBins;  /* drop high-frequency samples */

    /* allocate memory for various vectors */
    A = XLALCreateREAL8Vector(nBins);
    dAMc = XLALCreateREAL8Vector(nBins);
    dAEta = XLALCreateREAL8Vector(nBins);
    dAChi = XLALCreateREAL8Vector(nBins);
    dAT0 = XLALCreateREAL8Vector(nBins);
    dAPhi0 = XLALCreateREAL8Vector(nBins);
    dPsiMc = XLALCreateREAL8Vector(nBins);
    dPsiEta = XLALCreateREAL8Vector(nBins);
    dPsiChi = XLALCreateREAL8Vector(nBins);
    dPsiT0 = XLALCreateREAL8Vector(nBins);
    dPsiPhi0 = XLALCreateREAL8Vector(nBins);

    /* create a frequency vector from fLow to fCut with frequency resolution df */
    for (;k--;) {
        const REAL8 f = fLow + k * df;

        /* compute the amplitude of the frequency-domain waveform */
        A->data[k] = XLALSimInspiralTaylorF2RedSpinAofF(mc, eta, chi, f);

        /* compute the square of the template norm */
        hSqr += A->data[k] * A->data[k] / Shdata.data[k];

        /* compute the waveform deratives with respect to the parameters */
        dAMc->data[k] = XLALSimInspiralTaylorF2RedSpinDerivAMChirp(mc, eta, chi, f);
        dAEta->data[k] = XLALSimInspiralTaylorF2RedSpinDerivAEta(mc, eta, chi, f);
        dAChi->data[k] = XLALSimInspiralTaylorF2RedSpinDerivAChi(mc, eta, chi, f);
        dAT0->data[k] = 0.;
        dAPhi0->data[k] = 0.;
        dPsiMc->data[k] = XLALSimInspiralTaylorF2RedSpinDerivPsiMChirp(mc, eta, chi, f);
        dPsiEta->data[k] = XLALSimInspiralTaylorF2RedSpinDerivPsiEta(mc, eta, chi, f);
        dPsiChi->data[k] = XLALSimInspiralTaylorF2RedSpinDerivPsiChi(mc, eta, chi, f);
        dPsiT0->data[k] = LAL_TWOPI * f;
        dPsiPhi0->data[k] = 1.;
    }
    hSqr *= 4 * df;

    /* allocate memory, and initialize the Fisher matrix */
    gsl_matrix * g = gsl_matrix_calloc (5, 5);

    /* compute the components of the Fisher matrix in coordinates mc, eta, chi, t0, phi0 */
    gsl_matrix_set (g, 0,0, MetricCoeffs(A, dPsiMc, dPsiMc, dAMc, dAMc, &Shdata, hSqr, df));
    gsl_matrix_set (g, 0,1, MetricCoeffs(A, dPsiMc, dPsiEta, dAMc, dAEta, &Shdata, hSqr, df));
    gsl_matrix_set (g, 0,2, MetricCoeffs(A, dPsiMc, dPsiChi, dAMc, dAChi, &Shdata, hSqr, df));
    gsl_matrix_set (g, 0,3, MetricCoeffs(A, dPsiMc, dPsiT0, dAMc, dAT0, &Shdata, hSqr, df));
    gsl_matrix_set (g, 0,4, MetricCoeffs(A, dPsiMc, dPsiPhi0, dAMc, dAPhi0, &Shdata, hSqr, df));

    gsl_matrix_set (g, 1,0, gsl_matrix_get(g, 0,1));
    gsl_matrix_set (g, 1,1, MetricCoeffs(A, dPsiEta, dPsiEta, dAEta, dAEta, &Shdata, hSqr, df));
    gsl_matrix_set (g, 1,2, MetricCoeffs(A, dPsiEta, dPsiChi, dAEta, dAChi, &Shdata, hSqr, df));
    gsl_matrix_set (g, 1,3, MetricCoeffs(A, dPsiEta, dPsiT0, dAEta, dAT0, &Shdata, hSqr, df));
    gsl_matrix_set (g, 1,4, MetricCoeffs(A, dPsiEta, dPsiPhi0, dAEta, dAPhi0, &Shdata, hSqr, df));

    gsl_matrix_set (g, 2,0, gsl_matrix_get(g, 0,2));
    gsl_matrix_set (g, 2,1, gsl_matrix_get(g, 1,2));
    gsl_matrix_set (g, 2,2, MetricCoeffs(A, dPsiChi, dPsiChi, dAChi, dAChi, &Shdata, hSqr, df));
    gsl_matrix_set (g, 2,3, MetricCoeffs(A, dPsiChi, dPsiT0, dAChi, dAT0, &Shdata, hSqr, df));
    gsl_matrix_set (g, 2,4, MetricCoeffs(A, dPsiChi, dPsiPhi0, dAChi, dAPhi0, &Shdata, hSqr, df));

    gsl_matrix_set (g, 3,0, gsl_matrix_get(g, 0,3));
    gsl_matrix_set (g, 3,1, gsl_matrix_get(g, 1,3));
    gsl_matrix_set (g, 3,2, gsl_matrix_get(g, 2,3));
    gsl_matrix_set (g, 3,3, MetricCoeffs(A, dPsiT0, dPsiT0, dAT0, dAT0, &Shdata, hSqr, df));
    gsl_matrix_set (g, 3,4, MetricCoeffs(A, dPsiT0, dPsiPhi0, dAT0, dAPhi0, &Shdata, hSqr, df));

    gsl_matrix_set (g, 4,0, gsl_matrix_get(g, 0,4));
    gsl_matrix_set (g, 4,1, gsl_matrix_get(g, 1,4));
    gsl_matrix_set (g, 4,2, gsl_matrix_get(g, 2,4));
    gsl_matrix_set (g, 4,3, gsl_matrix_get(g, 3,4));
    gsl_matrix_set (g, 4,4, MetricCoeffs(A, dPsiPhi0, dPsiPhi0, dAPhi0, dAPhi0, &Shdata, hSqr, df));

    /* free the memory */
    XLALDestroyREAL8Vector(A);
    XLALDestroyREAL8Vector(dAMc);
    XLALDestroyREAL8Vector(dAEta);
    XLALDestroyREAL8Vector(dAChi);
    XLALDestroyREAL8Vector(dAT0);
    XLALDestroyREAL8Vector(dAPhi0);
    XLALDestroyREAL8Vector(dPsiMc);
    XLALDestroyREAL8Vector(dPsiEta);
    XLALDestroyREAL8Vector(dPsiChi);
    XLALDestroyREAL8Vector(dPsiT0);
    XLALDestroyREAL8Vector(dPsiPhi0);

    {
    /* Form submatrices g1, g2, g3, g4, defined as:
     *              g = [ g1 g2
     *                    g4 g3 ]                           */
    gsl_matrix_view g1v = gsl_matrix_submatrix (g, 0, 0, 3, 3);
    gsl_matrix_view g2v = gsl_matrix_submatrix (g, 0, 3, 3, 2);
    gsl_matrix_view g3v = gsl_matrix_submatrix (g, 3, 3, 2, 2);
    gsl_matrix_view g4v = gsl_matrix_submatrix (g, 3, 0, 2, 3);

    gsl_matrix * g1 = gsl_matrix_calloc (3, 3);
    gsl_matrix * g2 = gsl_matrix_calloc (3, 2);
    gsl_matrix * g3 = gsl_matrix_calloc (2, 2);
    gsl_matrix * g4 = gsl_matrix_calloc (2, 3);
    gsl_matrix * g3invg4 = gsl_matrix_calloc (2, 3);
    gsl_matrix * g2g3invg4 = gsl_matrix_calloc (3, 3);

    /* Project out the t0 and phi0 dimensions: gamma =  g1 - g2 g3^{-1} g4 */
    gsl_matrix * g3inv = gsl_matrix_calloc (2, 2);
    gsl_permutation * p = gsl_permutation_calloc (2);

    gsl_matrix_memcpy (g1, &g1v.matrix);
    gsl_matrix_memcpy (g2, &g2v.matrix);
    gsl_matrix_memcpy (g3, &g3v.matrix);
    gsl_matrix_memcpy (g4, &g4v.matrix);
    gsl_matrix_free (g);

    gsl_linalg_LU_decomp (g3, p, &s);
    gsl_linalg_LU_invert (g3, p, g3inv);
    gsl_permutation_free (p);
    gsl_matrix_free (g3);

    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, g3inv, g4,  0.0, g3invg4);
    gsl_matrix_free (g4);
    gsl_matrix_free (g3inv);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, g2, g3invg4,  0.0, g2g3invg4);
    gsl_matrix_free (g2);
    gsl_matrix_free (g3invg4);

    gsl_matrix_sub (g1, g2g3invg4);
    gsl_matrix_free (g2g3invg4);

    *gamma00 = gsl_matrix_get(g1, 0, 0);
    *gamma01 = gsl_matrix_get(g1, 0, 1);
    *gamma02 = gsl_matrix_get(g1, 0, 2);
    *gamma11 = gsl_matrix_get(g1, 1, 1);
    *gamma12 = gsl_matrix_get(g1, 1, 2);
    *gamma22 = gsl_matrix_get(g1, 2, 2);
    gsl_matrix_free (g1);
    }

    return XLAL_SUCCESS;
}


/**
 * Compute the template-space metric of "reduced-spin" PN templates in
 * theta0, theta3, theta3s parameter space.
 */
int XLALSimInspiralTaylorF2RedSpinComputeNoiseMoments(
    REAL8Vector *momI_0,    /**< noise moments: \f$momI_0(f) = \int_{f0}^f (f'/f0)^{(0-17)/3} df'\f$ */
    REAL8Vector *momI_2,    /**< noise moments: \f$momI_2(f) = \int_{f0}^f (f'/f0)^{(2-17)/3} df'\f$ */
    REAL8Vector *momI_3,    /**< noise moments: \f$momI_3(f) = \int_{f0}^f (f'/f0)^{(3-17)/3} df'\f$ */
    REAL8Vector *momI_4,    /**< noise moments: \f$momI_4(f) = \int_{f0}^f (f'/f0)^{(4-17)/3} df'\f$ */
    REAL8Vector *momI_5,    /**< noise moments: \f$momI_5(f) = \int_{f0}^f (f'/f0)^{(5-17)/3} df'\f$ */
    REAL8Vector *momI_6,    /**< noise moments: \f$momI_6(f) = \int_{f0}^f (f'/f0)^{(6-17)/3} df'\f$ */
    REAL8Vector *momI_7,    /**< noise moments: \f$momI_7(f) = \int_{f0}^f (f'/f0)^{(7-17)/3} df'\f$ */
    REAL8Vector *momI_8,    /**< noise moments: \f$momI_8(f) = \int_{f0}^f (f'/f0)^{(8-17)/3} df'\f$ */
    REAL8Vector *momI_9,    /**< noise moments: \f$momI_9(f) = \int_{f0}^f (f'/f0)^{(9-17)/3} df'\f$ */
    REAL8Vector *momI_10,   /**< noise moments: \f$momI_10(f) = \int_{f0}^f (f'/f0)^{(10-17)/3} df'\f$ */
    REAL8Vector *momI_11,   /**< noise moments: \f$momI_11(f) = \int_{f0}^f (f'/f0)^{(11-17)/3} df'\f$ */
    REAL8Vector *momI_12,   /**< noise moments: \f$momI_12(f) = \int_{f0}^f (f'/f0)^{(12-17)/3} df'\f$ */
    REAL8Vector *momI_13,   /**< noise moments: \f$momI_13(f) = \int_{f0}^f (f'/f0)^{(13-17)/3} df'\f$ */
    REAL8Vector *momI_14,   /**< noise moments: \f$momI_14(f) = \int_{f0}^f (f'/f0)^{(14-17)/3} df'\f$ */
    REAL8Vector *momI_15,   /**< noise moments: \f$momI_15(f) = \int_{f0}^f (f'/f0)^{(15-17)/3} df'\f$ */
    REAL8Vector *momI_16,   /**< noise moments: \f$momI_16(f) = \int_{f0}^f (f'/f0)^{(16-17)/3} df'\f$ */
    REAL8Vector *momJ_5,    /**< noise moments: \f$momJ_5(f) = \int_{f0}^f (f'/f0)^{(5-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_6,    /**< noise moments: \f$momJ_6(f) = \int_{f0}^f (f'/f0)^{(6-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_7,    /**< noise moments: \f$momJ_7(f) = \int_{f0}^f (f'/f0)^{(7-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_8,    /**< noise moments: \f$momJ_8(f) = \int_{f0}^f (f'/f0)^{(8-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_9,    /**< noise moments: \f$momJ_9(f) = \int_{f0}^f (f'/f0)^{(9-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_10,   /**< noise moments: \f$momJ_10(f) = \int_{f0}^f (f'/f0)^{(10-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_11,   /**< noise moments: \f$momJ_11(f) = \int_{f0}^f (f'/f0)^{(11-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_12,   /**< noise moments: \f$momJ_12(f) = \int_{f0}^f (f'/f0)^{(12-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_13,   /**< noise moments: \f$momJ_13(f) = \int_{f0}^f (f'/f0)^{(13-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_14,   /**< noise moments: \f$momJ_14(f) = \int_{f0}^f (f'/f0)^{(14-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momK_10,   /**< noise moments: \f$momK_10(f) = \int_{f0}^f (f'/f0)^{(10-17)/3} log(f'/f0) log(f'/f0) df'\f$ */
    REAL8Vector *momK_11,   /**< noise moments: \f$momK_11(f) = \int_{f0}^f (f'/f0)^{(11-17)/3} log(f'/f0) log(f'/f0) df'\f$ */
    REAL8Vector *momK_12,   /**< noise moments: \f$momK_12(f) = \int_{f0}^f (f'/f0)^{(12-17)/3} log(f'/f0) log(f'/f0) df'\f$ */
    REAL8Vector *Sh,         /**< one sided PSD of the detector noise: Sh(f) for f = [fLow, fNyq] */
    REAL8 fLow,             /**< low frequency cutoff (Hz) */
    REAL8 df)               /**< frequency resolution of the psd vector (Hz) */
{
    size_t i;
    const size_t ilow = (size_t) (fLow / df);  /* PSD starts at 0 Hz; others start at fLow */

    /* some minimum error checking */
    if ( !momI_0 || !momI_2 || !momI_3 || !momI_4 || !momI_5 || !momI_6 || !momI_7 || !momI_8 || !momI_9 || !momI_10 || !momI_11 || !momI_12 || !momI_13 || !momI_14 || !momI_15 || !momI_16 || !momJ_5 || !momJ_6 || !momJ_7 || !momJ_8 || !momJ_9 || !momJ_10 || !momJ_11 || !momJ_12 || !momJ_13 || !momJ_14 || !momK_10 || !momK_11 || !momK_12) {
        XLALPrintError("Moments not initialized");
        XLAL_ERROR(XLAL_EFAULT);
    }
    if (momI_0->length > (Sh->length - ilow)) {
        XLALPrintError("Sh not long enough to fill moment vectors");
        XLAL_ERROR(XLAL_EDOM);
    }

    momI_0->data[0] = 0.;
    momI_2->data[0] = 0.;
    momI_3->data[0] = 0.;
    momI_4->data[0] = 0.;
    momI_5->data[0] = 0.;
    momI_6->data[0] = 0.;
    momI_7->data[0] = 0.;
    momI_8->data[0] = 0.;
    momI_9->data[0] = 0.;
    momI_10->data[0] = 0.;
    momI_11->data[0] = 0.;
    momI_12->data[0] = 0.;
    momI_13->data[0] = 0.;
    momI_14->data[0] = 0.;
    momI_15->data[0] = 0.;
    momI_16->data[0] = 0.;

    momJ_5->data[0] = 0.;
    momJ_6->data[0] = 0.;
    momJ_7->data[0] = 0.;
    momJ_8->data[0] = 0.;
    momJ_9->data[0] = 0.;
    momJ_10->data[0] = 0.;
    momJ_11->data[0] = 0.;
    momJ_12->data[0] = 0.;
    momJ_13->data[0] = 0.;
    momJ_14->data[0] = 0.;

    momK_10->data[0] = 0.;
    momK_11->data[0] = 0.;
    momK_12->data[0] = 0.;

    const REAL8 fLowmSevenBythree = pow(fLow, -7./3.);
    const REAL8 fLowFac = fLowmSevenBythree * df;
    for (i=1; i < momI_0->length; i++) {
         const REAL8 psdfac = Sh->data[i + ilow];
         if (psdfac) {
             const REAL8 fbyfLow = (i*df+fLow)/fLow;
             const REAL8 logfbyfLow = log(fbyfLow);
             const REAL8 logfbyfLowSq = logfbyfLow * logfbyfLow;
             const REAL8 fbyfLowp2 = fbyfLow * fbyfLow;
             const REAL8 fbyfLowp3 = fbyfLowp2 * fbyfLow;
             const REAL8 fbyfLowp4 = fbyfLowp3 * fbyfLow;
             const REAL8 fbyfLowp5 = fbyfLowp4 * fbyfLow;
             const REAL8 cbrtfbyfLow = cbrt(fbyfLow);

             momI_0->data[i]  = momI_0->data[i-1]  + fLowFac/(fbyfLowp5*cbrtfbyfLow*cbrtfbyfLow)/psdfac; /* f^{-17/3} */
             momI_2->data[i]  = momI_2->data[i-1]  + fLowFac/fbyfLowp5/psdfac; /* f^{-15/3} */
             momI_3->data[i]  = momI_3->data[i-1]  + fLowFac/fbyfLowp5*cbrtfbyfLow/psdfac; /* f^{-14/3} */
             momI_4->data[i]  = momI_4->data[i-1]  + fLowFac/fbyfLowp4/cbrtfbyfLow/psdfac; /* f^{-13/3} */
             momI_5->data[i]  = momI_5->data[i-1]  + fLowFac/fbyfLowp4/psdfac; /* f^{-12/3} */
             momI_6->data[i]  = momI_6->data[i-1]  + fLowFac/fbyfLowp4*cbrtfbyfLow/psdfac; /* f^{-11/3} */
             momI_7->data[i]  = momI_7->data[i-1]  + fLowFac/fbyfLowp3/cbrtfbyfLow/psdfac; /* f^{-10/3} */
             momI_8->data[i]  = momI_8->data[i-1]  + fLowFac/fbyfLowp3/psdfac; /* f^{-9/3} */
             momI_9->data[i]  = momI_9->data[i-1]  + fLowFac/fbyfLowp3*cbrtfbyfLow/psdfac; /* f^{-8/3} */
             momI_10->data[i] = momI_10->data[i-1] + fLowFac/fbyfLowp2/cbrtfbyfLow/psdfac; /* f^{-7/3} */
             momI_11->data[i] = momI_11->data[i-1] + fLowFac/fbyfLowp2/psdfac; /* f^{-6/3} */
             momI_12->data[i] = momI_12->data[i-1] + fLowFac/fbyfLowp2*cbrtfbyfLow/psdfac; /* f^{-5/3} */
             momI_13->data[i] = momI_13->data[i-1] + fLowFac/fbyfLow/cbrtfbyfLow/psdfac; /* f^{-4/3} */
             momI_14->data[i] = momI_14->data[i-1] + fLowFac/fbyfLow/psdfac; /* f^{-3/3} */
             momI_15->data[i] = momI_15->data[i-1] + fLowFac/fbyfLow*cbrtfbyfLow/psdfac; /* f^{-2/3} */
             momI_16->data[i] = momI_16->data[i-1] + fLowFac/cbrtfbyfLow/psdfac; /* f^{-1/3} */

             momJ_5->data[i]  = momJ_5->data[i-1]  + fLowFac*logfbyfLow/fbyfLowp4/psdfac; /* f^{-12/3} */
             momJ_6->data[i]  = momJ_6->data[i-1]  + fLowFac*logfbyfLow/fbyfLowp4*cbrtfbyfLow/psdfac; /* f^{-11/3} */
             momJ_7->data[i]  = momJ_7->data[i-1]  + fLowFac*logfbyfLow/fbyfLowp3/cbrtfbyfLow/psdfac; /* f^{-10/3} */
             momJ_8->data[i]  = momJ_8->data[i-1]  + fLowFac*logfbyfLow/fbyfLowp3/psdfac; /* f^{-9/3} */
             momJ_9->data[i]  = momJ_9->data[i-1]  + fLowFac*logfbyfLow/fbyfLowp3*cbrtfbyfLow/psdfac; /* f^{-8/3} */
             momJ_10->data[i] = momJ_10->data[i-1] + fLowFac*logfbyfLow/fbyfLowp2/cbrtfbyfLow/psdfac; /* f^{-7/3} */
             momJ_11->data[i] = momJ_11->data[i-1] + fLowFac*logfbyfLow/fbyfLowp2/psdfac; /* f^{-6/3} */
             momJ_12->data[i] = momJ_12->data[i-1] + fLowFac*logfbyfLow/fbyfLowp2*cbrtfbyfLow/psdfac; /* f^{-5/3} */
             momJ_13->data[i] = momJ_13->data[i-1] + fLowFac*logfbyfLow/fbyfLow/cbrtfbyfLow/psdfac; /* f^{-4/3} */
             momJ_14->data[i] = momJ_14->data[i-1] + fLowFac*logfbyfLow/fbyfLow/psdfac; /* f^{-3/3} */

             momK_10->data[i] = momK_10->data[i-1] + fLowFac*logfbyfLowSq/fbyfLowp2/cbrtfbyfLow/psdfac; /* f^{-7/3} */
             momK_11->data[i] = momK_11->data[i-1] + fLowFac*logfbyfLowSq/fbyfLowp2/psdfac; /* f^{-6/3} */
             momK_12->data[i] = momK_12->data[i-1] + fLowFac*logfbyfLowSq/fbyfLowp2*cbrtfbyfLow/psdfac; /* f^{-5/3} */
         }
    }

    return XLAL_SUCCESS;
}


/**
 * Compute the Fisher information matrix of "reduced-spin" PN templates in
 * theta0, theta3, theta3s, t0, phi0 parameter space, for an SNR=1/sqrt(2) signal.
 */
gsl_matrix *XLALSimInspiralTaylorF2RedSpinFisherMatrixChirpTimes(
    const REAL8 theta0,     /**< dimensionless parameter related to the chirp time by theta0 = 2 pi fLow tau0 */
    const REAL8 theta3,     /**< dimensionless parameter related to the chirp time by theta3 = -2 pi fLow tau3 */
    const REAL8 theta3s,    /**< dimensionless parameter related to the chirp time by theta3s = 2 pi fLow tau3s */
    const REAL8 fLow,       /**< low-frequency cutoff (Hz) */
    const REAL8 df,         /**< frequency resolution of the noise moment vectors (Hz) */
    REAL8Vector *momI_0,     /**< noise moments: \f$momI_0(f) = \int_{f0}^f (f'/f0)^{(0-17)/3} df'\f$ */
    REAL8Vector *momI_2,     /**< noise moments: \f$momI_2(f) = \int_{f0}^f (f'/f0)^{(2-17)/3} df'\f$ */
    REAL8Vector *momI_3,     /**< noise moments: \f$momI_3(f) = \int_{f0}^f (f'/f0)^{(3-17)/3} df'\f$ */
    REAL8Vector *momI_4,     /**< noise moments: \f$momI_4(f) = \int_{f0}^f (f'/f0)^{(4-17)/3} df'\f$ */
    REAL8Vector *momI_5,     /**< noise moments: \f$momI_5(f) = \int_{f0}^f (f'/f0)^{(5-17)/3} df'\f$ */
    REAL8Vector *momI_6,     /**< noise moments: \f$momI_6(f) = \int_{f0}^f (f'/f0)^{(6-17)/3} df'\f$ */
    REAL8Vector *momI_7,     /**< noise moments: \f$momI_7(f) = \int_{f0}^f (f'/f0)^{(7-17)/3} df'\f$ */
    REAL8Vector *momI_8,     /**< noise moments: \f$momI_8(f) = \int_{f0}^f (f'/f0)^{(8-17)/3} df'\f$ */
    REAL8Vector *momI_9,     /**< noise moments: \f$momI_9(f) = \int_{f0}^f (f'/f0)^{(9-17)/3} df'\f$ */
    REAL8Vector *momI_10,    /**< noise moments: \f$momI_10(f) = \int_{f0}^f (f'/f0)^{(10-17)/3} df'\f$ */
    REAL8Vector *momI_11,    /**< noise moments: \f$momI_11(f) = \int_{f0}^f (f'/f0)^{(11-17)/3} df'\f$ */
    REAL8Vector *momI_12,    /**< noise moments: \f$momI_12(f) = \int_{f0}^f (f'/f0)^{(12-17)/3} df'\f$ */
    REAL8Vector *momI_13,    /**< noise moments: \f$momI_13(f) = \int_{f0}^f (f'/f0)^{(13-17)/3} df'\f$ */
    REAL8Vector *momI_14,    /**< noise moments: \f$momI_14(f) = \int_{f0}^f (f'/f0)^{(14-17)/3} df'\f$ */
    REAL8Vector *momI_15,    /**< noise moments: \f$momI_15(f) = \int_{f0}^f (f'/f0)^{(15-17)/3} df'\f$ */
    REAL8Vector *momI_16,    /**< noise moments: \f$momI_16(f) = \int_{f0}^f (f'/f0)^{(16-17)/3} df'\f$ */
    REAL8Vector *momJ_5,     /**< noise moments: \f$momJ_5(f) = \int_{f0}^f (f'/f0)^{(5-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_6,     /**< noise moments: \f$momJ_6(f) = \int_{f0}^f (f'/f0)^{(6-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_7,     /**< noise moments: \f$momJ_7(f) = \int_{f0}^f (f'/f0)^{(7-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_8,     /**< noise moments: \f$momJ_8(f) = \int_{f0}^f (f'/f0)^{(8-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_9,     /**< noise moments: \f$momJ_9(f) = \int_{f0}^f (f'/f0)^{(9-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_10,    /**< noise moments: \f$momJ_10(f) = \int_{f0}^f (f'/f0)^{(10-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_11,    /**< noise moments: \f$momJ_11(f) = \int_{f0}^f (f'/f0)^{(11-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_12,    /**< noise moments: \f$momJ_12(f) = \int_{f0}^f (f'/f0)^{(12-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_13,    /**< noise moments: \f$momJ_13(f) = \int_{f0}^f (f'/f0)^{(13-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_14,    /**< noise moments: \f$momJ_14(f) = \int_{f0}^f (f'/f0)^{(14-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momK_10,    /**< noise moments: \f$momK_14(f) = \int_{f0}^f (f'/f0)^{(14-17)/3} log(f'/f0) log(f'/f0) df'\f$ */
    REAL8Vector *momK_11,    /**< noise moments: \f$momK_15(f) = \int_{f0}^f (f'/f0)^{(15-17)/3} log(f'/f0) log(f'/f0) df'\f$ */
    REAL8Vector *momK_12     /**< noise moments: \f$momK_16(f) = \int_{f0}^f (f'/f0)^{(16-17)/3} log(f'/f0) log(f'/f0) df'\f$ */
) {

    if (theta0 <= 0) XLAL_ERROR_NULL(XLAL_EDOM);
    if (theta3 <= 0) XLAL_ERROR_NULL(XLAL_EDOM);
    if (fLow <= 0) XLAL_ERROR_NULL(XLAL_EDOM);
    if (df <= 0) XLAL_ERROR_NULL(XLAL_EDOM);

    REAL8 theta0_p2 = theta0 * theta0;
    REAL8 theta3_p2 = theta3 * theta3;
    REAL8 theta3s_p2 = theta3s * theta3s;
    REAL8 theta3_p3 = theta3_p2 * theta3;
    REAL8 theta3_p4 = theta3_p3 * theta3;
    REAL8 theta3_p5 = theta3_p4 * theta3;
    REAL8 theta0_p1by3 = cbrt(theta0);
    REAL8 theta3_p1by3 = cbrt(theta3);
    REAL8 theta0_p2by3 = theta0 / theta0_p1by3;
    REAL8 theta3_p2by3 = theta3 / theta3_p1by3;
    REAL8 theta3_p5by3 = theta3_p2 / theta3_p1by3;
    REAL8 theta0_p4by3 = theta0 * theta0_p1by3;
    REAL8 theta3_p4by3 = theta3 * theta3_p1by3;
    REAL8 theta0_p5by3 = theta0_p2 / theta0_p1by3;
    REAL8 theta3_p8by3 = theta3_p3 / theta3_p1by3;
    REAL8 theta0_p7by3 = theta0_p2 * theta0_p1by3;
    REAL8 theta3_p10by3 = theta3_p3 * theta3_p1by3;
    REAL8 theta3bytheta0_p2by3 = cbrt(theta3_p2/theta0_p2);

    /* harder to name combinations that come up a lot */
    REAL8 theta03_fac = -152*2.15443469003188*Pi_p5by3*theta0_p2by3 + 565*theta3_p5by3;
    REAL8 theta03_fac_p2 = theta03_fac * theta03_fac;
    REAL8 theta03_fac_p3 = theta03_fac_p2 * theta03_fac;

    /* compute the frequency bin corresponding to the cutoff frequency
     * (Schwarzschild test particle ISCO). Note that the first frequency
     * bin correspond to a frequency of fLow  */
    REAL8 fISCO = (8*sqrt(0.6666666666666666)*fLow*LAL_PI*theta0)/(15.*theta3);
    if (fISCO <= fLow) {
        XLALPrintError("fISCO <= fLow");
        XLAL_ERROR_NULL(XLAL_EDOM);
    }
    size_t iCut = (fISCO - fLow) / df;
    if (iCut > momI_10->length - 1) {
        XLALPrintWarning("moments do not cover fISCO (%g Hz); truncating to (%g Hz)", fISCO, (momI_10->length - 1) * df + fLow);
        iCut = momI_10->length - 1;
    }

    /* template norm */
    REAL8 hSqr = 2.*momI_10->data[iCut];

    /*******************************************************************/
    /* derivatives of the phasing coefficients w.r.t. parameters       */
    /*******************************************************************/
    REAL8 dPhi0byDtheta0 = 0.6;

    REAL8 dPhi2byDtheta0 = (11088*LAL_PI + (743*tenbyPi_p2by3*theta3_p5by3)/theta0_p2by3)/(12096.*theta3);

    REAL8 dPhi4byDtheta0 = (5429*2.92401773821287*cbrt(LAL_PI/(2.*theta0_p2*theta3)))/16128. - (15293365*1.7099759466767)/(3.2514048e7*1.5874010519682*Pi_p4by3/theta3bytheta0_p2by3/theta3bytheta0_p2by3) + (617*Pi_p2)/(384.*theta3_p2) - (13596750*4.64158883361278*Pi_p7by3*theta3_p8by3*theta3s_p2)/(theta0_p2by3*theta03_fac_p3) + (91125*1.7099759466767*Pi_p2by3*theta3_p8by3*theta3s_p2)/(2.*1.5874010519682*theta0_p4by3*theta03_fac_p2) + (2052000*Pi_p4*theta3*theta3s_p2)/theta03_fac_p3;

    REAL8 dPhi5byDtheta0 = (5*1.7099759466767*theta3_p2by3*(-873377*2.15443469003188*theta3 + 2345552*2.15443469003188*theta3s + (344829849600*Pi_p10by3*theta0_p4by3*theta3s)/theta03_fac_p2))/(1.0934784e7*twoPi_p2by3*theta0_p5by3);

    REAL8 dPhi6byDtheta0 = (-888945361920*Pi_p5*theta0_p2 + 33057275328*4.64158883361278*Pi_p10by3*theta0_p4by3*theta3_p5by3 + 9694463631160*2.15443469003188*Pi_p5by3*theta0_p2by3*theta3_p10by3 - 352848545280*2.15443469003188*Pi_11by3*theta0_p2by3*theta3_p10by3 + 3*(-11072977443251 + 1530761379840*LAL_GAMMA)*theta3_p5 + 3004298035200*Pi_p2*theta3_p5 + 1530761379840*theta3_p5*log((10*theta3)/(LAL_PI*theta0)))/(9.61375371264e11*Pi_p2*theta0_p2*theta3_p3);

    REAL8 dPhi7byDtheta0 = (-5*theta3_p2by3*(12718104*4.64158883361278*Pi_p5by3*theta0_p2by3 + 77096675*2.15443469003188*theta3_p5by3))/(2.60112384e8*Pi_p4by3*theta0_p7by3);

    REAL8 dPhi2byDtheta3 = (-5544*LAL_PI*theta0 + 743*tenbyPi_p2by3*theta0_p1by3*theta3_p5by3)/(6048.*theta3_p2);

    REAL8 dPhi3byDtheta3 = -1.5;

    REAL8 dPhi4byDtheta3 = ((-52242624*Pi_p2*theta0_p4by3)/theta3_p3 - (2736216*4.64158883361278*Pi_p1by3*theta0_p2by3)/theta3_p4by3 + (15293365*2.15443469003188*theta3_p1by3)/Pi_p4by3 + (740710656000*2.15443469003188*Pi_p2by3*theta3_p5by3*(608*2.15443469003188*Pi_p5by3*theta0_p2by3 + 565*theta3_p5by3)*theta3s_p2)/theta03_fac_p3 - (7315660800*4.64158883361278*Pi_p7by3*theta0_p2by3*(456*2.15443469003188*Pi_p5by3*theta0_p2by3 + 3955*theta3_p5by3)*theta3s_p2)/theta03_fac_p3)/(1.6257024e7*theta0_p1by3);

    REAL8 dPhi5byDtheta3 = (5*(640*2.15443469003188*Pi_p10by3*theta0_p4by3*theta3_p5by3*(13950845*theta3 - 36141056*theta3s) - 4000*Pi_p5by3*theta0_p2by3*theta3_p10by3*(16594163*theta3 - 12773768*theta3s) + 2825*4.64158883361278*theta3_p5*(4366885*theta3 - 4691104*theta3s) + 1000194048*4.64158883361278*Pi_p5*theta0_p2*theta3s))/(387072.*Pi_p2by3*theta0_p2by3*pow(152*2.15443469003188*Pi_p5by3*theta0_p2by3*theta3 - 565*theta3_p8by3,2));

    REAL8 dPhi6byDtheta3 = (1333418042880*Pi_p5*theta0_p2 - 66114550656*4.64158883361278*Pi_p10by3*theta0_p4by3*theta3_p5by3 - 4847231815580*2.15443469003188*Pi_p5by3*theta0_p2by3*theta3_p10by3 + 176424272640*2.15443469003188*Pi_11by3*theta0_p2by3*theta3_p10by3 + 3*(11328104339891 - 1530761379840*LAL_GAMMA)*theta3_p5 - 3004298035200*Pi_p2*theta3_p5 - 1530761379840*theta3_p5*log((10*theta3)/(LAL_PI*theta0)))/(4.80687685632e11*Pi_p2*theta0*theta3_p4);

    REAL8 dPhi7byDtheta3 = (5*(17059968*Pi_p10by3 + (7267488*4.64158883361278*Pi_p5by3*theta3_p5by3)/theta0_p2by3 + (77096675*2.15443469003188*theta3_p10by3)/theta0_p4by3))/(1.48635648e8*Pi_p4by3*theta3_p2);

    REAL8 dPhi3byDtheta3s = 1.5;

    REAL8 dPhi4byDtheta3s = (675*theta3*(8*4.64158883361278*Pi_p7by3*theta0_p2by3 - 405*2.15443469003188*Pi_p2by3*theta3_p5by3)*theta3s)/(2.*theta0_p1by3*theta03_fac_p2);

    REAL8 dPhi5byDtheta3s = (5*(-2785343*tenbyPi_p2by3*theta3bytheta0_p2by3 + (816*LAL_PI*(-2373 + (3301450*theta3_p5by3)/theta03_fac))/theta3))/1.7313408e7;

    REAL8 dPhi8byDt0 = 2.*fLow*LAL_PI;

    REAL8 dPhi5byDphi0 = 1.;

    REAL8 dPhiL5byDtheta0 = (5*1.7099759466767*theta3_p2by3*(-873377*2.15443469003188*theta3 + 2345552*2.15443469003188*theta3s + (344829849600*Pi_p10by3*theta0_p4by3*theta3s)/theta03_fac_p2))/(1.0934784e7*twoPi_p2by3*theta0_p5by3);

    REAL8 dPhiL6byDtheta0 = (535*theta3_p2)/(336.*Pi_p2*theta0_p2);

    REAL8 dPhiL5byDtheta3 = (5*(640*2.15443469003188*Pi_p10by3*theta0_p4by3*theta3_p5by3*(13950845*theta3 - 36141056*theta3s) - 4000*Pi_p5by3*theta0_p2by3*theta3_p10by3*(16594163*theta3 - 12773768*theta3s) + 2825*4.64158883361278*theta3_p5*(4366885*theta3 - 4691104*theta3s) + 1000194048*4.64158883361278*Pi_p5*theta0_p2*theta3s))/(387072.*Pi_p2by3*theta0_p2by3*pow(152*2.15443469003188*Pi_p5by3*theta0_p2by3*theta3 - 565*theta3_p8by3,2));

    REAL8 dPhiL6byDtheta3 = (-535*theta3)/(168.*Pi_p2*theta0);

    REAL8 dPhiL5byDtheta3s = (5*(-2785343*tenbyPi_p2by3*theta3bytheta0_p2by3 + (816*LAL_PI*(-2373 + (3301450*theta3_p5by3)/theta03_fac))/theta3))/1.7313408e7;

    /*******************************************************************/
    /*              (theta0, theta0) component of the metric           */
    /*******************************************************************/

    REAL8 g_theta0theta0 = 	dPhi0byDtheta0 * dPhi0byDtheta0 * momI_0->data[iCut] +
        dPhi0byDtheta0 * dPhi2byDtheta0 * momI_2->data[iCut] +
        dPhi0byDtheta0 * dPhi4byDtheta0 * momI_4->data[iCut] +
        dPhi0byDtheta0 * dPhi5byDtheta0 * momI_5->data[iCut] +
            dPhi0byDtheta0 * dPhiL5byDtheta0 * momJ_5->data[iCut] +
        dPhi0byDtheta0 * dPhi6byDtheta0 * momI_6->data[iCut] +
            dPhi0byDtheta0 * dPhiL6byDtheta0 * momJ_6->data[iCut] +
        dPhi0byDtheta0 * dPhi7byDtheta0 * momI_7->data[iCut] +
        dPhi2byDtheta0 * dPhi0byDtheta0 * momI_2->data[iCut] +
        dPhi2byDtheta0 * dPhi2byDtheta0 * momI_4->data[iCut] +
        dPhi2byDtheta0 * dPhi4byDtheta0 * momI_6->data[iCut] +
        dPhi2byDtheta0 * dPhi5byDtheta0 * momI_7->data[iCut] +
            dPhi2byDtheta0 * dPhiL5byDtheta0 * momJ_7->data[iCut] +
        dPhi2byDtheta0 * dPhi6byDtheta0 * momI_8->data[iCut] +
            dPhi2byDtheta0 * dPhiL6byDtheta0 * momJ_8->data[iCut] +
        dPhi2byDtheta0 * dPhi7byDtheta0 * momI_9->data[iCut] +
        dPhi4byDtheta0 * dPhi0byDtheta0 * momI_4->data[iCut] +
        dPhi4byDtheta0 * dPhi2byDtheta0 * momI_6->data[iCut] +
        dPhi4byDtheta0 * dPhi4byDtheta0 * momI_8->data[iCut] +
        dPhi4byDtheta0 * dPhi5byDtheta0 * momI_9->data[iCut] +
            dPhi4byDtheta0 * dPhiL5byDtheta0 * momJ_9->data[iCut] +
        dPhi4byDtheta0 * dPhi6byDtheta0 * momI_10->data[iCut] +
            dPhi4byDtheta0 * dPhiL6byDtheta0 * momJ_10->data[iCut] +
        dPhi4byDtheta0 * dPhi7byDtheta0 * momI_11->data[iCut] +
        dPhi5byDtheta0 * dPhi0byDtheta0 * momI_5->data[iCut] +
            dPhiL5byDtheta0 * dPhi0byDtheta0 * momJ_5->data[iCut] +
        dPhi5byDtheta0 * dPhi2byDtheta0 * momI_7->data[iCut] +
            dPhiL5byDtheta0 * dPhi2byDtheta0 * momJ_7->data[iCut] +
        dPhi5byDtheta0 * dPhi4byDtheta0 * momI_9->data[iCut] +
            dPhiL5byDtheta0 * dPhi4byDtheta0 * momJ_9->data[iCut] +
        dPhi5byDtheta0 * dPhi5byDtheta0 * momI_10->data[iCut] +
            dPhi5byDtheta0 * dPhiL5byDtheta0 * momJ_10->data[iCut] +
            dPhiL5byDtheta0 * dPhi5byDtheta0 * momJ_10->data[iCut] +
            dPhiL5byDtheta0 * dPhiL5byDtheta0 * momK_10->data[iCut] +
        dPhi5byDtheta0 * dPhi6byDtheta0 * momI_11->data[iCut] +
            dPhi5byDtheta0 * dPhiL6byDtheta0 * momJ_11->data[iCut] +
            dPhiL5byDtheta0 * dPhi6byDtheta0 * momJ_11->data[iCut] +
            dPhiL5byDtheta0 * dPhiL6byDtheta0 * momK_11->data[iCut] +
        dPhi5byDtheta0 * dPhi7byDtheta0 * momI_12->data[iCut] +
            dPhiL5byDtheta0 * dPhi7byDtheta0 * momJ_12->data[iCut] +
        dPhi6byDtheta0 * dPhi0byDtheta0 * momI_6->data[iCut] +
            dPhiL6byDtheta0 * dPhi0byDtheta0 * momJ_6->data[iCut] +
        dPhi6byDtheta0 * dPhi2byDtheta0 * momI_8->data[iCut] +
            dPhiL6byDtheta0 * dPhi2byDtheta0 * momJ_8->data[iCut] +
        dPhi6byDtheta0 * dPhi4byDtheta0 * momI_10->data[iCut] +
            dPhiL6byDtheta0 * dPhi4byDtheta0 * momJ_10->data[iCut] +
        dPhi6byDtheta0 * dPhi5byDtheta0 * momI_11->data[iCut] +
            dPhi6byDtheta0 * dPhiL5byDtheta0 * momJ_11->data[iCut] +
            dPhiL6byDtheta0 * dPhi5byDtheta0 * momJ_11->data[iCut] +
            dPhiL6byDtheta0 * dPhiL5byDtheta0 * momK_11->data[iCut] +
        dPhi6byDtheta0 * dPhi6byDtheta0 * momI_12->data[iCut] +
            dPhi6byDtheta0 * dPhiL6byDtheta0 * momJ_12->data[iCut] +
            dPhiL6byDtheta0 * dPhi6byDtheta0 * momJ_12->data[iCut] +
            dPhiL6byDtheta0 * dPhiL6byDtheta0 * momK_12->data[iCut] +
        dPhi6byDtheta0 * dPhi7byDtheta0 * momI_13->data[iCut] +
            dPhiL6byDtheta0 * dPhi7byDtheta0 * momJ_13->data[iCut] +
        dPhi7byDtheta0 * dPhi0byDtheta0 * momI_7->data[iCut] +
        dPhi7byDtheta0 * dPhi2byDtheta0 * momI_9->data[iCut] +
        dPhi7byDtheta0 * dPhi4byDtheta0 * momI_11->data[iCut] +
        dPhi7byDtheta0 * dPhi5byDtheta0 * momI_12->data[iCut] +
            dPhi7byDtheta0 * dPhiL5byDtheta0 * momJ_12->data[iCut] +
        dPhi7byDtheta0 * dPhi6byDtheta0 * momI_13->data[iCut] +
            dPhi7byDtheta0 * dPhiL6byDtheta0 * momJ_13->data[iCut] +
        dPhi7byDtheta0 * dPhi7byDtheta0 * momI_14->data[iCut] +
            0. ;

    /*******************************************************************/
    /*              (theta0, theta3) component of the metric           */
    /*******************************************************************/

    REAL8 g_theta0theta3 = 	dPhi0byDtheta0 * dPhi2byDtheta3 * momI_2->data[iCut] +
        dPhi0byDtheta0 * dPhi3byDtheta3 * momI_3->data[iCut] +
        dPhi0byDtheta0 * dPhi4byDtheta3 * momI_4->data[iCut] +
        dPhi0byDtheta0 * dPhi5byDtheta3 * momI_5->data[iCut] +
            dPhi0byDtheta0 * dPhiL5byDtheta3 * momJ_5->data[iCut] +
        dPhi0byDtheta0 * dPhi6byDtheta3 * momI_6->data[iCut] +
            dPhi0byDtheta0 * dPhiL6byDtheta3 * momJ_6->data[iCut] +
        dPhi0byDtheta0 * dPhi7byDtheta3 * momI_7->data[iCut] +
        dPhi2byDtheta0 * dPhi2byDtheta3 * momI_4->data[iCut] +
        dPhi2byDtheta0 * dPhi3byDtheta3 * momI_5->data[iCut] +
        dPhi2byDtheta0 * dPhi4byDtheta3 * momI_6->data[iCut] +
        dPhi2byDtheta0 * dPhi5byDtheta3 * momI_7->data[iCut] +
            dPhi2byDtheta0 * dPhiL5byDtheta3 * momJ_7->data[iCut] +
        dPhi2byDtheta0 * dPhi6byDtheta3 * momI_8->data[iCut] +
            dPhi2byDtheta0 * dPhiL6byDtheta3 * momJ_8->data[iCut] +
        dPhi2byDtheta0 * dPhi7byDtheta3 * momI_9->data[iCut] +
        dPhi4byDtheta0 * dPhi2byDtheta3 * momI_6->data[iCut] +
        dPhi4byDtheta0 * dPhi3byDtheta3 * momI_7->data[iCut] +
        dPhi4byDtheta0 * dPhi4byDtheta3 * momI_8->data[iCut] +
        dPhi4byDtheta0 * dPhi5byDtheta3 * momI_9->data[iCut] +
            dPhi4byDtheta0 * dPhiL5byDtheta3 * momJ_9->data[iCut] +
        dPhi4byDtheta0 * dPhi6byDtheta3 * momI_10->data[iCut] +
            dPhi4byDtheta0 * dPhiL6byDtheta3 * momJ_10->data[iCut] +
        dPhi4byDtheta0 * dPhi7byDtheta3 * momI_11->data[iCut] +
        dPhi5byDtheta0 * dPhi2byDtheta3 * momI_7->data[iCut] +
            dPhiL5byDtheta0 * dPhi2byDtheta3 * momJ_7->data[iCut] +
        dPhi5byDtheta0 * dPhi3byDtheta3 * momI_8->data[iCut] +
            dPhiL5byDtheta0 * dPhi3byDtheta3 * momJ_8->data[iCut] +
        dPhi5byDtheta0 * dPhi4byDtheta3 * momI_9->data[iCut] +
            dPhiL5byDtheta0 * dPhi4byDtheta3 * momJ_9->data[iCut] +
        dPhi5byDtheta0 * dPhi5byDtheta3 * momI_10->data[iCut] +
            dPhi5byDtheta0 * dPhiL5byDtheta3 * momJ_10->data[iCut] +
            dPhiL5byDtheta0 * dPhi5byDtheta3 * momJ_10->data[iCut] +
            dPhiL5byDtheta0 * dPhiL5byDtheta3 * momK_10->data[iCut] +
        dPhi5byDtheta0 * dPhi6byDtheta3 * momI_11->data[iCut] +
            dPhi5byDtheta0 * dPhiL6byDtheta3 * momJ_11->data[iCut] +
            dPhiL5byDtheta0 * dPhi6byDtheta3 * momJ_11->data[iCut] +
            dPhiL5byDtheta0 * dPhiL6byDtheta3 * momK_11->data[iCut] +
        dPhi5byDtheta0 * dPhi7byDtheta3 * momI_12->data[iCut] +
            dPhiL5byDtheta0 * dPhi7byDtheta3 * momJ_12->data[iCut] +
        dPhi6byDtheta0 * dPhi2byDtheta3 * momI_8->data[iCut] +
            dPhiL6byDtheta0 * dPhi2byDtheta3 * momJ_8->data[iCut] +
        dPhi6byDtheta0 * dPhi3byDtheta3 * momI_9->data[iCut] +
            dPhiL6byDtheta0 * dPhi3byDtheta3 * momJ_9->data[iCut] +
        dPhi6byDtheta0 * dPhi4byDtheta3 * momI_10->data[iCut] +
            dPhiL6byDtheta0 * dPhi4byDtheta3 * momJ_10->data[iCut] +
        dPhi6byDtheta0 * dPhi5byDtheta3 * momI_11->data[iCut] +
            dPhi6byDtheta0 * dPhiL5byDtheta3 * momJ_11->data[iCut] +
            dPhiL6byDtheta0 * dPhi5byDtheta3 * momJ_11->data[iCut] +
            dPhiL6byDtheta0 * dPhiL5byDtheta3 * momK_11->data[iCut] +
        dPhi6byDtheta0 * dPhi6byDtheta3 * momI_12->data[iCut] +
            dPhi6byDtheta0 * dPhiL6byDtheta3 * momJ_12->data[iCut] +
            dPhiL6byDtheta0 * dPhi6byDtheta3 * momJ_12->data[iCut] +
            dPhiL6byDtheta0 * dPhiL6byDtheta3 * momK_12->data[iCut] +
        dPhi6byDtheta0 * dPhi7byDtheta3 * momI_13->data[iCut] +
            dPhiL6byDtheta0 * dPhi7byDtheta3 * momJ_13->data[iCut] +
        dPhi7byDtheta0 * dPhi2byDtheta3 * momI_9->data[iCut] +
        dPhi7byDtheta0 * dPhi3byDtheta3 * momI_10->data[iCut] +
        dPhi7byDtheta0 * dPhi4byDtheta3 * momI_11->data[iCut] +
        dPhi7byDtheta0 * dPhi5byDtheta3 * momI_12->data[iCut] +
            dPhi7byDtheta0 * dPhiL5byDtheta3 * momJ_12->data[iCut] +
        dPhi7byDtheta0 * dPhi6byDtheta3 * momI_13->data[iCut] +
            dPhi7byDtheta0 * dPhiL6byDtheta3 * momJ_13->data[iCut] +
        dPhi7byDtheta0 * dPhi7byDtheta3 * momI_14->data[iCut] +
            0. ;

    /*******************************************************************/
    /*              (theta0, theta3s) component of the metric           */
    /*******************************************************************/

    REAL8 g_theta0theta3s = 	dPhi0byDtheta0 * dPhi3byDtheta3s * momI_3->data[iCut] +
        dPhi0byDtheta0 * dPhi4byDtheta3s * momI_4->data[iCut] +
        dPhi0byDtheta0 * dPhi5byDtheta3s * momI_5->data[iCut] +
            dPhi0byDtheta0 * dPhiL5byDtheta3s * momJ_5->data[iCut] +
        dPhi2byDtheta0 * dPhi3byDtheta3s * momI_5->data[iCut] +
        dPhi2byDtheta0 * dPhi4byDtheta3s * momI_6->data[iCut] +
        dPhi2byDtheta0 * dPhi5byDtheta3s * momI_7->data[iCut] +
            dPhi2byDtheta0 * dPhiL5byDtheta3s * momJ_7->data[iCut] +
        dPhi4byDtheta0 * dPhi3byDtheta3s * momI_7->data[iCut] +
        dPhi4byDtheta0 * dPhi4byDtheta3s * momI_8->data[iCut] +
        dPhi4byDtheta0 * dPhi5byDtheta3s * momI_9->data[iCut] +
            dPhi4byDtheta0 * dPhiL5byDtheta3s * momJ_9->data[iCut] +
        dPhi5byDtheta0 * dPhi3byDtheta3s * momI_8->data[iCut] +
            dPhiL5byDtheta0 * dPhi3byDtheta3s * momJ_8->data[iCut] +
        dPhi5byDtheta0 * dPhi4byDtheta3s * momI_9->data[iCut] +
            dPhiL5byDtheta0 * dPhi4byDtheta3s * momJ_9->data[iCut] +
        dPhi5byDtheta0 * dPhi5byDtheta3s * momI_10->data[iCut] +
            dPhi5byDtheta0 * dPhiL5byDtheta3s * momJ_10->data[iCut] +
            dPhiL5byDtheta0 * dPhi5byDtheta3s * momJ_10->data[iCut] +
            dPhiL5byDtheta0 * dPhiL5byDtheta3s * momK_10->data[iCut] +
        dPhi6byDtheta0 * dPhi3byDtheta3s * momI_9->data[iCut] +
            dPhiL6byDtheta0 * dPhi3byDtheta3s * momJ_9->data[iCut] +
        dPhi6byDtheta0 * dPhi4byDtheta3s * momI_10->data[iCut] +
            dPhiL6byDtheta0 * dPhi4byDtheta3s * momJ_10->data[iCut] +
        dPhi6byDtheta0 * dPhi5byDtheta3s * momI_11->data[iCut] +
            dPhi6byDtheta0 * dPhiL5byDtheta3s * momJ_11->data[iCut] +
            dPhiL6byDtheta0 * dPhi5byDtheta3s * momJ_11->data[iCut] +
            dPhiL6byDtheta0 * dPhiL5byDtheta3s * momK_11->data[iCut] +
        dPhi7byDtheta0 * dPhi3byDtheta3s * momI_10->data[iCut] +
        dPhi7byDtheta0 * dPhi4byDtheta3s * momI_11->data[iCut] +
        dPhi7byDtheta0 * dPhi5byDtheta3s * momI_12->data[iCut] +
            dPhi7byDtheta0 * dPhiL5byDtheta3s * momJ_12->data[iCut] +
            0. ;

    /*******************************************************************/
    /*              (theta0, t0) component of the metric           */
    /*******************************************************************/

    REAL8 g_theta0t0 = 	dPhi0byDtheta0 * dPhi8byDt0 * momI_8->data[iCut] +
        dPhi2byDtheta0 * dPhi8byDt0 * momI_10->data[iCut] +
        dPhi4byDtheta0 * dPhi8byDt0 * momI_12->data[iCut] +
        dPhi5byDtheta0 * dPhi8byDt0 * momI_13->data[iCut] +
            dPhiL5byDtheta0 * dPhi8byDt0 * momJ_13->data[iCut] +
        dPhi6byDtheta0 * dPhi8byDt0 * momI_14->data[iCut] +
            dPhiL6byDtheta0 * dPhi8byDt0 * momJ_14->data[iCut] +
        dPhi7byDtheta0 * dPhi8byDt0 * momI_15->data[iCut] +
            0. ;

    /*******************************************************************/
    /*              (theta0, phi0) component of the metric           */
    /*******************************************************************/

    REAL8 g_theta0phi0 = 	dPhi0byDtheta0 * dPhi5byDphi0 * momI_5->data[iCut] +
        dPhi2byDtheta0 * dPhi5byDphi0 * momI_7->data[iCut] +
        dPhi4byDtheta0 * dPhi5byDphi0 * momI_9->data[iCut] +
        dPhi5byDtheta0 * dPhi5byDphi0 * momI_10->data[iCut] +
            dPhiL5byDtheta0 * dPhi5byDphi0 * momJ_10->data[iCut] +
        dPhi6byDtheta0 * dPhi5byDphi0 * momI_11->data[iCut] +
            dPhiL6byDtheta0 * dPhi5byDphi0 * momJ_11->data[iCut] +
        dPhi7byDtheta0 * dPhi5byDphi0 * momI_12->data[iCut] +
            0. ;

    /*******************************************************************/
    /*              (theta3, theta3) component of the metric           */
    /*******************************************************************/

    REAL8 g_theta3theta3 = 	dPhi2byDtheta3 * dPhi2byDtheta3 * momI_4->data[iCut] +
        dPhi2byDtheta3 * dPhi3byDtheta3 * momI_5->data[iCut] +
        dPhi2byDtheta3 * dPhi4byDtheta3 * momI_6->data[iCut] +
        dPhi2byDtheta3 * dPhi5byDtheta3 * momI_7->data[iCut] +
            dPhi2byDtheta3 * dPhiL5byDtheta3 * momJ_7->data[iCut] +
        dPhi2byDtheta3 * dPhi6byDtheta3 * momI_8->data[iCut] +
            dPhi2byDtheta3 * dPhiL6byDtheta3 * momJ_8->data[iCut] +
        dPhi2byDtheta3 * dPhi7byDtheta3 * momI_9->data[iCut] +
        dPhi3byDtheta3 * dPhi2byDtheta3 * momI_5->data[iCut] +
        dPhi3byDtheta3 * dPhi3byDtheta3 * momI_6->data[iCut] +
        dPhi3byDtheta3 * dPhi4byDtheta3 * momI_7->data[iCut] +
        dPhi3byDtheta3 * dPhi5byDtheta3 * momI_8->data[iCut] +
            dPhi3byDtheta3 * dPhiL5byDtheta3 * momJ_8->data[iCut] +
        dPhi3byDtheta3 * dPhi6byDtheta3 * momI_9->data[iCut] +
            dPhi3byDtheta3 * dPhiL6byDtheta3 * momJ_9->data[iCut] +
        dPhi3byDtheta3 * dPhi7byDtheta3 * momI_10->data[iCut] +
        dPhi4byDtheta3 * dPhi2byDtheta3 * momI_6->data[iCut] +
        dPhi4byDtheta3 * dPhi3byDtheta3 * momI_7->data[iCut] +
        dPhi4byDtheta3 * dPhi4byDtheta3 * momI_8->data[iCut] +
        dPhi4byDtheta3 * dPhi5byDtheta3 * momI_9->data[iCut] +
            dPhi4byDtheta3 * dPhiL5byDtheta3 * momJ_9->data[iCut] +
        dPhi4byDtheta3 * dPhi6byDtheta3 * momI_10->data[iCut] +
            dPhi4byDtheta3 * dPhiL6byDtheta3 * momJ_10->data[iCut] +
        dPhi4byDtheta3 * dPhi7byDtheta3 * momI_11->data[iCut] +
        dPhi5byDtheta3 * dPhi2byDtheta3 * momI_7->data[iCut] +
            dPhiL5byDtheta3 * dPhi2byDtheta3 * momJ_7->data[iCut] +
        dPhi5byDtheta3 * dPhi3byDtheta3 * momI_8->data[iCut] +
            dPhiL5byDtheta3 * dPhi3byDtheta3 * momJ_8->data[iCut] +
        dPhi5byDtheta3 * dPhi4byDtheta3 * momI_9->data[iCut] +
            dPhiL5byDtheta3 * dPhi4byDtheta3 * momJ_9->data[iCut] +
        dPhi5byDtheta3 * dPhi5byDtheta3 * momI_10->data[iCut] +
            dPhi5byDtheta3 * dPhiL5byDtheta3 * momJ_10->data[iCut] +
            dPhiL5byDtheta3 * dPhi5byDtheta3 * momJ_10->data[iCut] +
            dPhiL5byDtheta3 * dPhiL5byDtheta3 * momK_10->data[iCut] +
        dPhi5byDtheta3 * dPhi6byDtheta3 * momI_11->data[iCut] +
            dPhi5byDtheta3 * dPhiL6byDtheta3 * momJ_11->data[iCut] +
            dPhiL5byDtheta3 * dPhi6byDtheta3 * momJ_11->data[iCut] +
            dPhiL5byDtheta3 * dPhiL6byDtheta3 * momK_11->data[iCut] +
        dPhi5byDtheta3 * dPhi7byDtheta3 * momI_12->data[iCut] +
            dPhiL5byDtheta3 * dPhi7byDtheta3 * momJ_12->data[iCut] +
        dPhi6byDtheta3 * dPhi2byDtheta3 * momI_8->data[iCut] +
            dPhiL6byDtheta3 * dPhi2byDtheta3 * momJ_8->data[iCut] +
        dPhi6byDtheta3 * dPhi3byDtheta3 * momI_9->data[iCut] +
            dPhiL6byDtheta3 * dPhi3byDtheta3 * momJ_9->data[iCut] +
        dPhi6byDtheta3 * dPhi4byDtheta3 * momI_10->data[iCut] +
            dPhiL6byDtheta3 * dPhi4byDtheta3 * momJ_10->data[iCut] +
        dPhi6byDtheta3 * dPhi5byDtheta3 * momI_11->data[iCut] +
            dPhi6byDtheta3 * dPhiL5byDtheta3 * momJ_11->data[iCut] +
            dPhiL6byDtheta3 * dPhi5byDtheta3 * momJ_11->data[iCut] +
            dPhiL6byDtheta3 * dPhiL5byDtheta3 * momK_11->data[iCut] +
        dPhi6byDtheta3 * dPhi6byDtheta3 * momI_12->data[iCut] +
            dPhi6byDtheta3 * dPhiL6byDtheta3 * momJ_12->data[iCut] +
            dPhiL6byDtheta3 * dPhi6byDtheta3 * momJ_12->data[iCut] +
            dPhiL6byDtheta3 * dPhiL6byDtheta3 * momK_12->data[iCut] +
        dPhi6byDtheta3 * dPhi7byDtheta3 * momI_13->data[iCut] +
            dPhiL6byDtheta3 * dPhi7byDtheta3 * momJ_13->data[iCut] +
        dPhi7byDtheta3 * dPhi2byDtheta3 * momI_9->data[iCut] +
        dPhi7byDtheta3 * dPhi3byDtheta3 * momI_10->data[iCut] +
        dPhi7byDtheta3 * dPhi4byDtheta3 * momI_11->data[iCut] +
        dPhi7byDtheta3 * dPhi5byDtheta3 * momI_12->data[iCut] +
            dPhi7byDtheta3 * dPhiL5byDtheta3 * momJ_12->data[iCut] +
        dPhi7byDtheta3 * dPhi6byDtheta3 * momI_13->data[iCut] +
            dPhi7byDtheta3 * dPhiL6byDtheta3 * momJ_13->data[iCut] +
        dPhi7byDtheta3 * dPhi7byDtheta3 * momI_14->data[iCut] +
            0. ;

    /*******************************************************************/
    /*              (theta3, theta3s) component of the metric           */
    /*******************************************************************/

    REAL8 g_theta3theta3s = 	dPhi2byDtheta3 * dPhi3byDtheta3s * momI_5->data[iCut] +
        dPhi2byDtheta3 * dPhi4byDtheta3s * momI_6->data[iCut] +
        dPhi2byDtheta3 * dPhi5byDtheta3s * momI_7->data[iCut] +
            dPhi2byDtheta3 * dPhiL5byDtheta3s * momJ_7->data[iCut] +
        dPhi3byDtheta3 * dPhi3byDtheta3s * momI_6->data[iCut] +
        dPhi3byDtheta3 * dPhi4byDtheta3s * momI_7->data[iCut] +
        dPhi3byDtheta3 * dPhi5byDtheta3s * momI_8->data[iCut] +
            dPhi3byDtheta3 * dPhiL5byDtheta3s * momJ_8->data[iCut] +
        dPhi4byDtheta3 * dPhi3byDtheta3s * momI_7->data[iCut] +
        dPhi4byDtheta3 * dPhi4byDtheta3s * momI_8->data[iCut] +
        dPhi4byDtheta3 * dPhi5byDtheta3s * momI_9->data[iCut] +
            dPhi4byDtheta3 * dPhiL5byDtheta3s * momJ_9->data[iCut] +
        dPhi5byDtheta3 * dPhi3byDtheta3s * momI_8->data[iCut] +
            dPhiL5byDtheta3 * dPhi3byDtheta3s * momJ_8->data[iCut] +
        dPhi5byDtheta3 * dPhi4byDtheta3s * momI_9->data[iCut] +
            dPhiL5byDtheta3 * dPhi4byDtheta3s * momJ_9->data[iCut] +
        dPhi5byDtheta3 * dPhi5byDtheta3s * momI_10->data[iCut] +
            dPhi5byDtheta3 * dPhiL5byDtheta3s * momJ_10->data[iCut] +
            dPhiL5byDtheta3 * dPhi5byDtheta3s * momJ_10->data[iCut] +
            dPhiL5byDtheta3 * dPhiL5byDtheta3s * momK_10->data[iCut] +
        dPhi6byDtheta3 * dPhi3byDtheta3s * momI_9->data[iCut] +
            dPhiL6byDtheta3 * dPhi3byDtheta3s * momJ_9->data[iCut] +
        dPhi6byDtheta3 * dPhi4byDtheta3s * momI_10->data[iCut] +
            dPhiL6byDtheta3 * dPhi4byDtheta3s * momJ_10->data[iCut] +
        dPhi6byDtheta3 * dPhi5byDtheta3s * momI_11->data[iCut] +
            dPhi6byDtheta3 * dPhiL5byDtheta3s * momJ_11->data[iCut] +
            dPhiL6byDtheta3 * dPhi5byDtheta3s * momJ_11->data[iCut] +
            dPhiL6byDtheta3 * dPhiL5byDtheta3s * momK_11->data[iCut] +
        dPhi7byDtheta3 * dPhi3byDtheta3s * momI_10->data[iCut] +
        dPhi7byDtheta3 * dPhi4byDtheta3s * momI_11->data[iCut] +
        dPhi7byDtheta3 * dPhi5byDtheta3s * momI_12->data[iCut] +
            dPhi7byDtheta3 * dPhiL5byDtheta3s * momJ_12->data[iCut] +
            0. ;

    /*******************************************************************/
    /*              (theta3, t0) component of the metric           */
    /*******************************************************************/

    REAL8 g_theta3t0 = 	dPhi2byDtheta3 * dPhi8byDt0 * momI_10->data[iCut] +
        dPhi3byDtheta3 * dPhi8byDt0 * momI_11->data[iCut] +
        dPhi4byDtheta3 * dPhi8byDt0 * momI_12->data[iCut] +
        dPhi5byDtheta3 * dPhi8byDt0 * momI_13->data[iCut] +
            dPhiL5byDtheta3 * dPhi8byDt0 * momJ_13->data[iCut] +
        dPhi6byDtheta3 * dPhi8byDt0 * momI_14->data[iCut] +
            dPhiL6byDtheta3 * dPhi8byDt0 * momJ_14->data[iCut] +
        dPhi7byDtheta3 * dPhi8byDt0 * momI_15->data[iCut] +
            0. ;

    /*******************************************************************/
    /*              (theta3, phi0) component of the metric           */
    /*******************************************************************/

    REAL8 g_theta3phi0 = 	dPhi2byDtheta3 * dPhi5byDphi0 * momI_7->data[iCut] +
        dPhi3byDtheta3 * dPhi5byDphi0 * momI_8->data[iCut] +
        dPhi4byDtheta3 * dPhi5byDphi0 * momI_9->data[iCut] +
        dPhi5byDtheta3 * dPhi5byDphi0 * momI_10->data[iCut] +
            dPhiL5byDtheta3 * dPhi5byDphi0 * momJ_10->data[iCut] +
        dPhi6byDtheta3 * dPhi5byDphi0 * momI_11->data[iCut] +
            dPhiL6byDtheta3 * dPhi5byDphi0 * momJ_11->data[iCut] +
        dPhi7byDtheta3 * dPhi5byDphi0 * momI_12->data[iCut] +
            0. ;

    /*******************************************************************/
    /*              (theta3s, theta3s) component of the metric           */
    /*******************************************************************/

    REAL8 g_theta3stheta3s = 	dPhi3byDtheta3s * dPhi3byDtheta3s * momI_6->data[iCut] +
        dPhi3byDtheta3s * dPhi4byDtheta3s * momI_7->data[iCut] +
        dPhi3byDtheta3s * dPhi5byDtheta3s * momI_8->data[iCut] +
            dPhi3byDtheta3s * dPhiL5byDtheta3s * momJ_8->data[iCut] +
        dPhi4byDtheta3s * dPhi3byDtheta3s * momI_7->data[iCut] +
        dPhi4byDtheta3s * dPhi4byDtheta3s * momI_8->data[iCut] +
        dPhi4byDtheta3s * dPhi5byDtheta3s * momI_9->data[iCut] +
            dPhi4byDtheta3s * dPhiL5byDtheta3s * momJ_9->data[iCut] +
        dPhi5byDtheta3s * dPhi3byDtheta3s * momI_8->data[iCut] +
            dPhiL5byDtheta3s * dPhi3byDtheta3s * momJ_8->data[iCut] +
        dPhi5byDtheta3s * dPhi4byDtheta3s * momI_9->data[iCut] +
            dPhiL5byDtheta3s * dPhi4byDtheta3s * momJ_9->data[iCut] +
        dPhi5byDtheta3s * dPhi5byDtheta3s * momI_10->data[iCut] +
            dPhi5byDtheta3s * dPhiL5byDtheta3s * momJ_10->data[iCut] +
            dPhiL5byDtheta3s * dPhi5byDtheta3s * momJ_10->data[iCut] +
            dPhiL5byDtheta3s * dPhiL5byDtheta3s * momK_10->data[iCut] +
            0. ;

    /*******************************************************************/
    /*              (theta3s, t0) component of the metric           */
    /*******************************************************************/

    REAL8 g_theta3st0 = 	dPhi3byDtheta3s * dPhi8byDt0 * momI_11->data[iCut] +
        dPhi4byDtheta3s * dPhi8byDt0 * momI_12->data[iCut] +
        dPhi5byDtheta3s * dPhi8byDt0 * momI_13->data[iCut] +
            dPhiL5byDtheta3s * dPhi8byDt0 * momJ_13->data[iCut] +
            0. ;

    /*******************************************************************/
    /*              (theta3s, phi0) component of the metric           */
    /*******************************************************************/

    REAL8 g_theta3sphi0 = 	dPhi3byDtheta3s * dPhi5byDphi0 * momI_8->data[iCut] +
        dPhi4byDtheta3s * dPhi5byDphi0 * momI_9->data[iCut] +
        dPhi5byDtheta3s * dPhi5byDphi0 * momI_10->data[iCut] +
            dPhiL5byDtheta3s * dPhi5byDphi0 * momJ_10->data[iCut] +
            0. ;

    /*******************************************************************/
    /*              (t0, t0) component of the metric           */
    /*******************************************************************/

    REAL8 g_t0t0 = 	dPhi8byDt0 * dPhi8byDt0 * momI_16->data[iCut] +
            0. ;

    /*******************************************************************/
    /*              (t0, phi0) component of the metric           */
    /*******************************************************************/

    REAL8 g_t0phi0 = 	dPhi8byDt0 * dPhi5byDphi0 * momI_13->data[iCut] +
            0. ;

    /*******************************************************************/
    /*              (phi0, phi0) component of the metric           */
    /*******************************************************************/

    REAL8 g_phi0phi0 = 	dPhi5byDphi0 * dPhi5byDphi0 * momI_10->data[iCut] +
            0. ;

    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    /* allocate memory, and initialize the Fisher matrix */
    gsl_matrix * g = gsl_matrix_calloc (5, 5);

    /* compute the components of the Fisher matrix in coordinates theta0, theta3, theta3s, t0, phi0 */
    gsl_matrix_set (g, 0,0, g_theta0theta0/hSqr);
    gsl_matrix_set (g, 0,1, g_theta0theta3/hSqr);
    gsl_matrix_set (g, 0,2, g_theta0theta3s/hSqr);
    gsl_matrix_set (g, 0,3, g_theta0t0/hSqr);
    gsl_matrix_set (g, 0,4, g_theta0phi0/hSqr);

    gsl_matrix_set (g, 1,0, gsl_matrix_get(g, 0,1));
    gsl_matrix_set (g, 1,1, g_theta3theta3/hSqr);
    gsl_matrix_set (g, 1,2, g_theta3theta3s/hSqr);
    gsl_matrix_set (g, 1,3, g_theta3t0/hSqr);
    gsl_matrix_set (g, 1,4, g_theta3phi0/hSqr);

    gsl_matrix_set (g, 2,0, gsl_matrix_get(g, 0,2));
    gsl_matrix_set (g, 2,1, gsl_matrix_get(g, 1,2));
    gsl_matrix_set (g, 2,2, g_theta3stheta3s/hSqr);
    gsl_matrix_set (g, 2,3, g_theta3st0/hSqr);
    gsl_matrix_set (g, 2,4, g_theta3sphi0/hSqr);

    gsl_matrix_set (g, 3,0, gsl_matrix_get(g, 0,3));
    gsl_matrix_set (g, 3,1, gsl_matrix_get(g, 1,3));
    gsl_matrix_set (g, 3,2, gsl_matrix_get(g, 2,3));
    gsl_matrix_set (g, 3,3, g_t0t0/hSqr);
    gsl_matrix_set (g, 3,4, g_t0phi0/hSqr);

    gsl_matrix_set (g, 4,0, gsl_matrix_get(g, 0,4));
    gsl_matrix_set (g, 4,1, gsl_matrix_get(g, 1,4));
    gsl_matrix_set (g, 4,2, gsl_matrix_get(g, 2,4));
    gsl_matrix_set (g, 4,3, gsl_matrix_get(g, 3,4));
    gsl_matrix_set (g, 4,4, g_phi0phi0/hSqr);

    return g;
}


/**
 * Compute the template-space metric of "reduced-spin" PN templates in
 * theta0, theta3, theta3s parameter space.
 */
int XLALSimInspiralTaylorF2RedSpinMetricChirpTimes(
    REAL8 *gamma00,         /**< template metric coeff. 00 in theta0-theta3-theta3s*/
    REAL8 *gamma01,         /**< template metric coeff. 01/10 in theta0-theta3-theta3s */
    REAL8 *gamma02,         /**< template metric coeff. 02/20 in theta0-theta3-theta3s */
    REAL8 *gamma11,         /**< template metric coeff. 11 in theta0-theta3-theta3s */
    REAL8 *gamma12,         /**< template metric coeff. 12/21 in theta0-theta3-theta3s */
    REAL8 *gamma22,         /**< template metric coeff. 22 in theta0-theta3-theta3s */
    const REAL8 theta0,     /**< dimensionless parameter related to the chirp time by theta0 = 2 pi fLow tau0 */
    const REAL8 theta3,     /**< dimensionless parameter related to the chirp time by theta3 = -2 pi fLow tau3 */
    const REAL8 theta3s,    /**< dimensionless parameter related to the chirp time by theta3s = 2 pi fLow tau3s */
    const REAL8 fLow,       /**< low-frequency cutoff (Hz) */
    const REAL8 df,         /**< frequency resolution of the noise moment vectors (Hz) */
    REAL8Vector *momI_0,     /**< noise moments: \f$momI_0(f) = \int_{f0}^f (f'/f0)^{(0-17)/3} df'\f$ */
    REAL8Vector *momI_2,     /**< noise moments: \f$momI_2(f) = \int_{f0}^f (f'/f0)^{(2-17)/3} df'\f$ */
    REAL8Vector *momI_3,     /**< noise moments: \f$momI_3(f) = \int_{f0}^f (f'/f0)^{(3-17)/3} df'\f$ */
    REAL8Vector *momI_4,     /**< noise moments: \f$momI_4(f) = \int_{f0}^f (f'/f0)^{(4-17)/3} df'\f$ */
    REAL8Vector *momI_5,     /**< noise moments: \f$momI_5(f) = \int_{f0}^f (f'/f0)^{(5-17)/3} df'\f$ */
    REAL8Vector *momI_6,     /**< noise moments: \f$momI_6(f) = \int_{f0}^f (f'/f0)^{(6-17)/3} df'\f$ */
    REAL8Vector *momI_7,     /**< noise moments: \f$momI_7(f) = \int_{f0}^f (f'/f0)^{(7-17)/3} df'\f$ */
    REAL8Vector *momI_8,     /**< noise moments: \f$momI_8(f) = \int_{f0}^f (f'/f0)^{(8-17)/3} df'\f$ */
    REAL8Vector *momI_9,     /**< noise moments: \f$momI_9(f) = \int_{f0}^f (f'/f0)^{(9-17)/3} df'\f$ */
    REAL8Vector *momI_10,    /**< noise moments: \f$momI_10(f) = \int_{f0}^f (f'/f0)^{(10-17)/3} df'\f$ */
    REAL8Vector *momI_11,    /**< noise moments: \f$momI_11(f) = \int_{f0}^f (f'/f0)^{(11-17)/3} df'\f$ */
    REAL8Vector *momI_12,    /**< noise moments: \f$momI_12(f) = \int_{f0}^f (f'/f0)^{(12-17)/3} df'\f$ */
    REAL8Vector *momI_13,    /**< noise moments: \f$momI_13(f) = \int_{f0}^f (f'/f0)^{(13-17)/3} df'\f$ */
    REAL8Vector *momI_14,    /**< noise moments: \f$momI_14(f) = \int_{f0}^f (f'/f0)^{(14-17)/3} df'\f$ */
    REAL8Vector *momI_15,    /**< noise moments: \f$momI_15(f) = \int_{f0}^f (f'/f0)^{(15-17)/3} df'\f$ */
    REAL8Vector *momI_16,    /**< noise moments: \f$momI_16(f) = \int_{f0}^f (f'/f0)^{(16-17)/3} df'\f$ */
    REAL8Vector *momJ_5,     /**< noise moments: \f$momJ_5(f) = \int_{f0}^f (f'/f0)^{(5-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_6,     /**< noise moments: \f$momJ_6(f) = \int_{f0}^f (f'/f0)^{(6-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_7,     /**< noise moments: \f$momJ_7(f) = \int_{f0}^f (f'/f0)^{(7-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_8,     /**< noise moments: \f$momJ_8(f) = \int_{f0}^f (f'/f0)^{(8-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_9,     /**< noise moments: \f$momJ_9(f) = \int_{f0}^f (f'/f0)^{(9-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_10,    /**< noise moments: \f$momJ_10(f) = \int_{f0}^f (f'/f0)^{(10-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_11,    /**< noise moments: \f$momJ_11(f) = \int_{f0}^f (f'/f0)^{(11-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_12,    /**< noise moments: \f$momJ_12(f) = \int_{f0}^f (f'/f0)^{(12-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_13,    /**< noise moments: \f$momJ_13(f) = \int_{f0}^f (f'/f0)^{(13-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momJ_14,    /**< noise moments: \f$momJ_14(f) = \int_{f0}^f (f'/f0)^{(14-17)/3} log(f'/f0) df'\f$ */
    REAL8Vector *momK_10,    /**< noise moments: \f$momK_14(f) = \int_{f0}^f (f'/f0)^{(14-17)/3} log(f'/f0) log(f'/f0) df'\f$ */
    REAL8Vector *momK_11,    /**< noise moments: \f$momK_15(f) = \int_{f0}^f (f'/f0)^{(15-17)/3} log(f'/f0) log(f'/f0) df'\f$ */
    REAL8Vector *momK_12     /**< noise moments: \f$momK_16(f) = \int_{f0}^f (f'/f0)^{(16-17)/3} log(f'/f0) log(f'/f0) df'\f$ */
) {
    gsl_matrix * g = XLALSimInspiralTaylorF2RedSpinFisherMatrixChirpTimes(
        theta0, theta3, theta3s, fLow, df, momI_0, momI_2, momI_3, momI_4,
        momI_5, momI_6, momI_7, momI_8, momI_9, momI_10, momI_11, momI_12,
        momI_13, momI_14, momI_15, momI_16, momJ_5, momJ_6, momJ_7, momJ_8,
        momJ_9, momJ_10, momJ_11, momJ_12, momJ_13, momJ_14, momK_10, momK_11,
        momK_12);

    int s = 0;

    if (!g)
        XLAL_ERROR(XLAL_FAILURE);

    /* Form submatrices g1, g2, g3, g4, defined as:
     *              g = [ g1 g2
     *                    g4 g3 ]                           */
    gsl_matrix_view g1v = gsl_matrix_submatrix (g, 0, 0, 3, 3);
    gsl_matrix_view g2v = gsl_matrix_submatrix (g, 0, 3, 3, 2);
    gsl_matrix_view g3v = gsl_matrix_submatrix (g, 3, 3, 2, 2);
    gsl_matrix_view g4v = gsl_matrix_submatrix (g, 3, 0, 2, 3);

    gsl_matrix * g1 = gsl_matrix_calloc (3, 3);
    gsl_matrix * g2 = gsl_matrix_calloc (3, 2);
    gsl_matrix * g3 = gsl_matrix_calloc (2, 2);
    gsl_matrix * g4 = gsl_matrix_calloc (2, 3);
    gsl_matrix * g3invg4 = gsl_matrix_calloc (2, 3);
    gsl_matrix * g2g3invg4 = gsl_matrix_calloc (3, 3);

    /* Project out the t0 and phi0 dimensions: gamma =  g1 - g2 g3^{-1} g4 */
    gsl_matrix * g3inv = gsl_matrix_calloc (2, 2);
    gsl_permutation * p = gsl_permutation_calloc (2);

    gsl_matrix_memcpy (g1, &g1v.matrix);
    gsl_matrix_memcpy (g2, &g2v.matrix);
    gsl_matrix_memcpy (g3, &g3v.matrix);
    gsl_matrix_memcpy (g4, &g4v.matrix);
    gsl_matrix_free (g);

    gsl_linalg_LU_decomp (g3, p, &s);
    gsl_linalg_LU_invert (g3, p, g3inv);
    gsl_permutation_free (p);
    gsl_matrix_free (g3);

    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, g3inv, g4,  0.0, g3invg4);
    gsl_matrix_free (g4);
    gsl_matrix_free (g3inv);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, g2, g3invg4,  0.0, g2g3invg4);
    gsl_matrix_free (g2);
    gsl_matrix_free (g3invg4);

    gsl_matrix_sub (g1, g2g3invg4);
    gsl_matrix_free (g2g3invg4);

    *gamma00 = gsl_matrix_get(g1, 0, 0);
    *gamma01 = gsl_matrix_get(g1, 0, 1);
    *gamma02 = gsl_matrix_get(g1, 0, 2);
    *gamma11 = gsl_matrix_get(g1, 1, 1);
    *gamma12 = gsl_matrix_get(g1, 1, 2);
    *gamma22 = gsl_matrix_get(g1, 2, 2);
    gsl_matrix_free (g1);

    return XLAL_SUCCESS;
}


/* compute theta0, theta3, theta3s from mc, eta, chi */
void XLALSimInspiralTaylorF2RedSpinChirpTimesFromMchirpEtaChi(
    double *theta0, /**< dimensionless parameter related to the chirp time by theta0 = 2 pi fLow tau0 */
    double *theta3, /**< dimensionless parameter related to the chirp time by theta3 = -2 pi fLow tau3 */
    double *theta3s,/**< dimensionless parameter related to the chirp time by theta3s = 2 pi fLow tau3s */
    double mc,      /**< chirp mass (M_sun) */
    double eta,     /**< symmetric mass ratio  */
    double chi,     /**< reduced-spin parameter */
    double fLow)    /**< low-frequency cutoff (Hz) */
{
    *theta0 = 5./(128.*pow(LAL_PI*mc*LAL_MTSUN_SI*fLow,1.6666666666666667));
    *theta3 = cbrt(LAL_PI)/(4.*pow(mc*LAL_MTSUN_SI*fLow,0.6666666666666666)*pow(eta,0.6));
    *theta3s = (113.*chi)/(192.*pow(LAL_PI*mc*LAL_MTSUN_SI*fLow,0.6666666666666666)*pow(eta,0.6));
}


/* compute mc, eta, chi from theta0, theta3, theta3s */
void XLALSimInspiralTaylorF2RedSpinMchirpEtaChiFromChirpTimes(
    double *mc,     /**< chirp mass (M_sun) */
    double *eta,    /**< symmetric mass ratio  */
    double *chi,    /**< reduced-spin parameter */
    double theta0,  /**< dimensionless parameter related to the chirp time by theta0 = 2 pi fLow tau0 */
    double theta3,  /**< dimensionless parameter related to the chirp time by theta3 = -2 pi fLow tau3 */
    double theta3s, /**< dimensionless parameter related to the chirp time by theta3s = 2 pi fLow tau3s */
    double fLow)    /**< low-frequency cutoff (Hz) */
{
    *mc = pow(2.*theta0*theta0*theta0/125., -1./5.)/(16.*LAL_PI*fLow)/LAL_MTSUN_SI;
    *eta = cbrt(16.*pow(LAL_PI,5.)*theta0*theta0/ (25.*pow(theta3,5.)));
    *chi = 48.*LAL_PI*theta3s/(113.*theta3);
}
