/*
*  Copyright (C) 2011 Drew Keppel, 2012 Riccardo Sturani
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

#include <lal/LALConstants.h>
#include <lal/LALAtomicDatatypes.h>

#include <math.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * Computes the PN Coefficients for using in the PN energy equation.
 *
 * Terms given in equation 3.1 of: Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 * For the spin terms a good reference are (3.15) and (3.16) of 1303.7412
 *
 * In the latest version coefficients of the terms n.S and L.S are reported
 * "Averaged" spin coefficients refer to the ones obtained by orbital averaging,
 * i.e. by using
 * n_i n_j = 1/2 (\f$\delta_{ij} - \hat LN_i \hat LN_j\f$)
 * However such orbital averaging at 2PN would introduce corrections
 * at 3PN, as LNh is not constant.
 */

static REAL8 UNUSED
XLALSimInspiralPNEnergy_0PNCoeff(
	REAL8 eta)
{
	return -eta / 2.0;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_2PNCoeff(
	REAL8 eta)
{
	return -(0.75 + eta/12.0);
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_4PNCoeff(
	REAL8 eta)
{
	return -(27.0/8.0 - 19.0/8.0 * eta + 1./24.0 * eta*eta);
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_6PNCoeff(
	REAL8 eta)
{
	return -(67.5/6.4 - (344.45/5.76 - 20.5/9.6 * LAL_PI*LAL_PI) * eta + 15.5/9.6 * eta*eta + 3.5/518.4 * eta*eta*eta);
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_8PNCoeff(
        REAL8 eta)
{
        return (-39.69/1.28 + (-123.671/5.76 + 9.037/1.536 *LAL_PI*LAL_PI+ 1792./15.*log(2)+89.6/1.5*LAL_GAMMA)* eta + (-498.449/3.456 +31.57/5.76*LAL_PI*LAL_PI)*eta*eta + 3.01/17.28 *eta*eta*eta + 7.7/3110.4*eta*eta*eta*eta);
        /*see arXiv:1305.4884, or eq.(26) of arXiv:1309.3474
          note that in getting a4 from PRD 62, 084011 (2000),
          the first reference is using the fact that \omega_{static} = 0
          (see arXiv:gr-qc/0105038) */
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_8PNLogCoeff(
        REAL8 eta)
{
	return 896./15.*eta;
        /* arXiv:1305.4884 has a misprint, it should have a factor of nu
           See for instance arXiv:1002.0726
           Also note that this term is usually given as 448*log(x)/15
           since x=v^2 the log(v) term is twice this */
}

/*  Eq. (4.6) of arXiv:1212.5520
 */
static REAL8 UNUSED
XLALSimInspiralPNEnergy_3PNSOCoeff(
	REAL8 mByM)
{
	return 2. / 3. + 2. / mByM;
}

/*  Eq. (6) of arXiv:astro-ph/0504538v2
 */
static REAL8 UNUSED
XLALSimInspiralPNEnergy_4PNS1S2CoeffAvg(
	REAL8 eta)
{
	return 1./eta;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_4PNS1S2Coeff(
	REAL8 eta)
{
	return -2./eta;
}

/*  Eq. (6) of arXiv:astro-ph/0504538v2
 */
static REAL8 UNUSED
XLALSimInspiralPNEnergy_4PNS1OS2OCoeffAvg(
	REAL8 eta)
{
	return -3./eta;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_4PNS1nS2nCoeff(
	REAL8 UNUSED eta)
{
	return 6./eta;
}

/*  Eq. (6) of arXiv:astro-ph/0504538v2
 */
static REAL8 UNUSED
XLALSimInspiralPNEnergy_4PNQMS1S1CoeffAvg(
	REAL8 mByM)
{
	return .5/mByM/mByM;
}

/*  Eq. (6) of arXiv:astro-ph/0504538v2
 */
static REAL8 UNUSED
XLALSimInspiralPNEnergy_4PNQMS1OS1OCoeffAvg(
	REAL8 mByM)
{
	return -1.5/mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_4PNQMS1S1Coeff(
	REAL8 mByM)
{
	return -1./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_4PNQMS1nS1nCoeff(
	REAL8 mByM)
{
	return 3./mByM/mByM;
}

/*  Eq. 4.6 of arXiv:1212.5520
 */
static REAL8 UNUSED
XLALSimInspiralPNEnergy_5PNSOCoeff(
	REAL8 mByM)
{
        return 5./3. + 3./mByM + 29.*mByM/9. + mByM*mByM/9.;
}

/*  From (3.30) of arXiv:1501.01529
 */
static REAL8 UNUSED
XLALSimInspiralPNEnergy_6PNS1S2Coeff(
	REAL8 eta)
{
	return -7./eta -1./3.;
}

/*  From (3.30) of arXiv:1501.01529
 */
static REAL8 UNUSED
XLALSimInspiralPNEnergy_6PNS1OS2OCoeff(
	REAL8 eta)
{
	return 16./3./eta - 2./9.;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_6PNS1nS2nCoeff(
	REAL8 eta)
{
        return 13./eta - 3.;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_6PNS1vS2vCoeff(
	REAL8 eta)
{
        return 5./eta;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_6PNS1S2CoeffAvg(
	REAL8 eta)
{
	return 2./eta -11./6.;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_6PNS1OS2OCoeffAvg(
	REAL8 eta)
{
  return -11./3./eta + 2.3/1.8;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_6PNS1S1Coeff(
	REAL8 mByM)
{
        return 2./(mByM*mByM) - 1./mByM -1.;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_6PNS1OS1OCoeff(
	REAL8 mByM)
{
        return 3./(mByM*mByM) -2./(3.*mByM) -1./9.;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_6PNS1nS1nCoeff(
	REAL8 mByM)
{
         return -8./(mByM*mByM)+11./3./mByM+1.;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_6PNS1vS1vCoeff(
	REAL8 mByM)
{
        return 2./(mByM*mByM)-2./mByM;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_6PNS1S1CoeffAvg(
	REAL8 mByM)
{
        return -1./(mByM*mByM) - 1./6./mByM -0.5;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_6PNS1OS1OCoeffAvg(
	REAL8 mByM)
{
        return 6./(mByM*mByM) -1.5/mByM -1.1/1.8;
}

/*  From (3.30) of arXiv:1501.01529
 */
static REAL8 UNUSED
XLALSimInspiralPNEnergy_6PNQMS1S1Coeff(
	REAL8 mByM)
{
	return -2.5/mByM/mByM - 2.5/mByM - 5./6.;
}

/*  From (3.30) of arXiv:1501.01529
 */
static REAL8 UNUSED
XLALSimInspiralPNEnergy_6PNQMS1nS1nCoeff(
	REAL8 mByM)
{
	return 6.5/mByM/mByM + 8./5./mByM + 2.5;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_6PNQMS1vS1vCoeff(
	REAL8 mByM)
{
	return 1./mByM/mByM -1./mByM;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_6PNQMS1S1CoeffAvg(
	REAL8 mByM)
{
	return 1.25/mByM/mByM + 1.25/mByM + 5./12.;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_6PNQMS1OS1OCoeffAvg(
	REAL8 mByM)
{
	return -3.75/mByM/mByM - 3.75/mByM - 1.25;
}

/*  Eq. (4.6) of arXiv:1212.5520
 *  Symbol definitions right above eq. (3.1)
 */
static REAL8 UNUSED
XLALSimInspiralPNEnergy_7PNSOCoeff(
	REAL8 mByM)
{
	return -75./4. + 27./(4.*mByM) + 53.*mByM/2. + 67*mByM*mByM/6. + 17.*mByM*mByM*mByM/12. - mByM*mByM*mByM*mByM/12.;
}

/*
 * Tidal correction coefficients to Energy
 */

static REAL8 UNUSED
XLALSimInspiralPNEnergy_10PNTidalCoeff(
	REAL8 mByM)
{
       return -9.0 * mByM*mByM*mByM*mByM*(1.-mByM);
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_12PNTidalCoeff(
	REAL8 mByM)
{
  return (-33./2. + 11./2.*mByM - 11./2.*mByM*mByM + 33./2.*mByM*mByM*mByM)*mByM*mByM*mByM*mByM;
}

/**
 * Computes the flux PN Coefficients.
 *
 * Terms given in equation 3.2 of: Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 * For terms involving spins see eq.(3.13) of arXiv:1303.7412
 */

static REAL8 UNUSED
XLALSimInspiralPNFlux_0PNCoeff(
	REAL8 eta)
{
	return 32.0 * eta*eta / 5.0;
}


static REAL8 UNUSED
XLALSimInspiralPNFlux_2PNCoeff(
	REAL8 eta)
{
	return -(12.47/3.36 + 3.5/1.2 * eta);
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_3PNCoeff(
	REAL8 UNUSED eta)
{
	return 4.0 * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_3PNSOCoeff(
	REAL8 mByM)
{
	return -3./2. - 5./4./mByM;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNCoeff(
	REAL8 eta)
{
	return -(44.711/9.072 - 92.71/5.04 * eta - 6.5/1.8 * eta*eta);
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNS1S2Coeff(
    REAL8 eta)
{
    return 31./8./eta;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNS1nS2nCoeff(
    REAL8 eta)
{
    return -15./eta;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNS1vS2vCoeff(
    REAL8 eta)
{
    return 71./24./eta;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNS1S1Coeff(
    REAL8 mByM)
{
    return 1./(16.*mByM*mByM);
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNS1vS1vCoeff(
    REAL8 mByM)
{
    return 1./48./(mByM*mByM);
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNS1S2CoeffAvg(
    REAL8 eta)
{
    return -103./48./eta;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNS1OS2OCoeffAvg(
    REAL8 eta)
{
    return 289./48./eta;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNS1S1CoeffAvg(
    REAL8 mByM)
{
    return 7./96./(mByM*mByM);
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNS1OS1OCoeffAvg(
    REAL8 mByM)
{
    return -1./96./(mByM*mByM);
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNQMS1S1Coeff(
    REAL8 mByM)
{
    return 2./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNQMS1nS1nCoeff(
    REAL8 mByM)
{
    return -7.5/mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNQMS1vS1vCoeff(
    REAL8 mByM)
{
    return 1.5/mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNQMS1S1CoeffAvg(
    REAL8 mByM)
{
    return -1./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNQMS1OS1OCoeffAvg(
    REAL8 mByM)
{
    return 3./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_5PNCoeff(
	REAL8 eta)
{
	return -(81.91/6.72 + 58.3/2.4 * eta) * LAL_PI;
}

/* Eq. (4.9) of arXiv:1307.6793
 */
static REAL8 UNUSED
XLALSimInspiralPNFlux_5PNSOCoeff(
	REAL8 mByM)
{
	return 63./8. - 13./(16.*mByM) - (73.*mByM)/36. - (157.*mByM*mByM)/18.;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNCoeff(
	REAL8 eta)
{
        return (664.3739519/6.9854400 + 16.0/3.0 * LAL_PI*LAL_PI - 17.12/1.05 * LAL_GAMMA - 17.12/1.05*log(4.) + (4.1/4.8 * LAL_PI*LAL_PI - 134.543/7.776) * eta - 94.403/3.024 * eta*eta - 7.75/3.24 * eta*eta*eta);
}

/* Note that this coefficient multiplies log(v)*/
static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNLogCoeff(
	REAL8 UNUSED eta)
{
	return -17.12/1.05;
}

/* Eq. (4.9) of arXiv:1307.6793
 * (symbol definitions around eq. 3.1)
 */

static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNSOCoeff(
	REAL8 mByM)
{
	return LAL_PI*( -17./3. - 31./(6.*mByM) );
}

/* From (4.12) of 1501.01529
 */
static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNS1S2Coeff(
        REAL8 eta)
{
	return -2.9/16.8/eta - 25.9/1.2;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNS1OS2OCoeff(
        REAL8 eta)
{
        return -49./6./eta + 44./9.;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNS1nS2nCoeff(
        REAL8 eta)
{
        return 349.9/4.2/eta +  117.59/1.68;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNS1vS2vCoeff(
        REAL8 eta)
{
        return -38.9/1.2/eta - 20.27/5.04;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNS1S2CoeffAvg(
        REAL8 eta)
{
        return 212.3/8.4/eta + 82.1/7.2;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNS1OS2OCoeffAvg(
        REAL8 eta)
{
        return -564.7/16.8/eta - 202.3/7.2;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNS1S1Coeff(
        REAL8 mByM)
{
        return -21./(8.*mByM*mByM) + 21.5/2.4/mByM - 1./24.;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNS1OS1OCoeff(
        REAL8 mByM)
{
        return -.5/(mByM*mByM) - 43./(6.*mByM) + 22./9.;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNS1nS1nCoeff(
        REAL8 mByM)
{
        return 23.53/(1.12*mByM*mByM) - 6.47/1.68/mByM + 2.27/3.36;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNS1vS1vCoeff(
        REAL8 mByM)
{
        return 8.81/(1.12*mByM*mByM) - 36.67/2.52/mByM + 6.1/100.8;
}


static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNS1S1CoeffAvg(
        REAL8 mByM)
{
        return 18.9/(1.6*mByM*mByM) - 3.5/14.4/mByM +4.7/14.4;
}


static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNS1OS1OCoeffAvg(
        REAL8 mByM)
{
        return -23.9/(1.6*mByM*mByM) + 2.93/1.44/mByM + 2.99/1.44;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNQMS1S1Coeff(
        REAL8 mByM)
{
        return  -27.9/(5.6*mByM*mByM) - 45./(8.*mByM) + 43./4.;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNQMS1nS1nCoeff(
        REAL8 mByM)
{
        return 154.1/(4.2*mByM*mByM) + 15.17/(1.68*mByM) - 96.1/2.8;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNQMS1vS1vCoeff(
    REAL8 mByM)
{
    return  -36.53/(1.68*mByM*mByM) + 65.9/(8.4*mByM) + 2.9/1.4;
}


static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNQMS1S1CoeffAvg(
        REAL8 mByM)
{
        return  27.9/(11.2*mByM*mByM) + 4.5/(1.6*mByM) - 43./8.;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNQMS1OS1OCoeffAvg(
        REAL8 mByM)
{
        return -8.37/(1.12*mByM*mByM) - 13.5/(1.6*mByM) + 12.9/8.;
}

/*
 */
static REAL8 UNUSED
XLALSimInspiralPNFlux_7PNCoeff(
	REAL8 eta)
{
	return -(162.85/5.04 - 214.745/1.728 * eta - 193.385/3.024 * eta*eta) * LAL_PI;
}

/* Eq. (4.9) of arXiv:1307.6793
 */
static REAL8 UNUSED
XLALSimInspiralPNFlux_7PNSOCoeff(
	REAL8 mByM)
{
        return (380.647/13.608) + 95.35/(3.36*mByM) - 401.15*mByM/7.56 + 3742.*mByM*mByM/63. - 35.*mByM*mByM*mByM/108. - 1117.*mByM*mByM*mByM*mByM/54.;
}

/* Eq. (4.9) of arXiv:1307.6793
 */
static REAL8 UNUSED
XLALSimInspiralPNFlux_8PNSOCoeff(
	REAL8 mByM)
{
        return LAL_PI * (125.47/2.52 - 71.63/(6.72*mByM) -3.137*mByM/2.016 - 212.41*mByM*mByM/3.36);
}

/*
 * Tidal correction coefficients to Flux
 */

static REAL8 UNUSED
XLALSimInspiralPNFlux_10PNTidalCoeff(
	REAL8 mByM)
{
        return 6. *(3. - 2.*mByM) * mByM*mByM*mByM*mByM;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_12PNTidalCoeff(
	REAL8 mByM)
{
        return (-176./7. - 1803./28.*mByM + 643./4.*mByM*mByM -155./2.*mByM*mByM*mByM) * mByM*mByM*mByM*mByM;
}

/* Non-spin phasing terms - see arXiv:0907.0700, Eq. 3.18 */
static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_2PNCoeff(
        REAL8 eta
    )
{
        return 5.*(74.3/8.4 + 11.*eta)/9.;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_3PNCoeff(
        REAL8 UNUSED eta)
{
        return -16.*LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_4PNCoeff(
        REAL8 eta
    )
{
        return 5.*(3058.673/7.056 + 5429./7.*eta+617.*eta*eta)/72.;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_5PNCoeff(
        REAL8 eta
    )
{
        return 5./9.*(772.9/8.4-13.*eta)*LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_5PNLogCoeff(
        REAL8 eta
    )
{
        return 5./3.*(772.9/8.4-13.*eta)*LAL_PI;
}

static REAL8 UNUSED XLALSimInspiralTaylorF2Phasing_6PNLogCoeff(
        REAL8 UNUSED eta
    )
{
  return -684.8/2.1;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_6PNCoeff(
        REAL8 eta
    )
{
  return 11583.231236531/4.694215680 - 640./3.*LAL_PI*LAL_PI - 684.8/2.1*LAL_GAMMA + eta*(-15737.765635/3.048192 + 225.5/1.2*LAL_PI*LAL_PI) + eta*eta*76.055/1.728 - eta*eta*eta*127.825/1.296 + XLALSimInspiralTaylorF2Phasing_6PNLogCoeff(eta)*log(4.);
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_7PNCoeff(
        REAL8 eta
    )
{
        return LAL_PI*(770.96675/2.54016 + 378.515/1.512*eta - 740.45/7.56*eta*eta);
}

/* Spin-orbit terms - can be derived from arXiv:1303.7412, Eq. 3.15-16 */

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_3PNSOCoeff(
        REAL8 mByM
    )
{
        return  mByM*(25.+38./3.*mByM);
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_5PNSOCoeff(
        REAL8 mByM
    )
{
  return -mByM*(1391.5/8.4-mByM*(1.-mByM)*10./3.+ mByM*(1276./8.1+mByM*(1.-mByM)*170./9.));

}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_6PNSOCoeff(
        REAL8 mByM)
{
  return LAL_PI*mByM*(1490./3. + mByM*260.);
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_7PNSOCoeff(
        REAL8 mByM
    )
{
  REAL8 eta=mByM*(1.-mByM);
  return mByM*(-17097.8035/4.8384+eta*28764.25/6.72+eta*eta*47.35/1.44 + mByM*(-7189.233785/1.524096+eta*458.555/3.024-eta*eta*534.5/7.2));
}

/*
 * Spin-squared corrections to TF2 phasing
 * Compute 2.0PN SS, QM, and self-spin
 * See Eq. (6.24) in arXiv:0810.5336
 * 9b,c,d in arXiv:astro-ph/0504538
 * Note that these terms are understood to multiply
 * dimensionless spin magnitudes \chi_i=S_i/m_i^2
 * differently from the convention adopted for SpinTaylorTX
 * whose spinning coefficients multiply \chi_LAL=S_i/M^2
 * where M=m_1+m_2.
 * See also https://dcc.ligo.org/T1800298
 */

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_4PNS1S2Coeff(
        REAL8 eta
    )
{
  return  247./4.8*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_4PNS1S2OCoeff(
        REAL8 eta
    )
{
  return  -721./4.8*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_4PNQM2SOCoeff(
        REAL8 mByM
    )
{
  return -720./9.6*mByM*mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_4PNSelf2SOCoeff(
        REAL8 mByM
    )
{
  return 1./9.6*mByM*mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_4PNQM2SCoeff(
        REAL8 mByM
    )
{
  return 240./9.6*mByM*mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_4PNSelf2SCoeff(
        REAL8 mByM
    )
{
  return -7./9.6*mByM*mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_6PNS1S2OCoeff(
        REAL8 eta
    )
{
  return (326.75/1.12L + 557.5/1.8*eta)*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_6PNSelf2SCoeff(
        REAL8 mByM
    )
{
  return (-4108.25/6.72-108.5/1.2*mByM+125.5/3.6*mByM*mByM)*mByM*mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_6PNQM2SCoeff(
        REAL8 mByM
    )
{
  return (4703.5/8.4+2935./6.*mByM-120.*mByM*mByM)*mByM*mByM;
}

/*
 * Tidal corrections to F2 phasing
 * See arXiv:1101.1673
 */

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_10PNTidalCoeff(
        REAL8 mByM /**< ratio of object mass to total mass */
    )
{
  return (-288. + 264.*mByM)*mByM*mByM*mByM*mByM;

}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_12PNTidalCoeff(
        REAL8 mByM /**< ratio of object mass to total mass */
    )
{
  return (-15895./28. + 4595./28.*mByM + 5715./14.*mByM*mByM - 325./7.*mByM*mByM*mByM)*mByM*mByM*mByM*mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_13PNTidalCoeff(
      REAL8 mByM /**< ratio of object mass to total mass */
    )
/*  literature: Agathos et al (arxiv 1503.0545) eq (5)
 * the coefficient mByM4 conversion & transformation (6.5PN, 7PN, 7.5PN):
 * mByM=mA/M: mA= mass star A, M is total mass (mA+mB)
 * Lambda (unitless) = lambda(m) / mA^5
 * to call the function:
 * Lambda * XLALSimInspiralTaylorF2Phasing_13PNTidalCoeff
 * lambda(m)*mByM^4/mA^5= lambda(m)*(mA/M)^4/(mA)^5= lambda/(M^4*mA)
 * =lambda/(mByM*M^5) eq (5)
 */
{
  return mByM*mByM*mByM*mByM * 24.L*(12.L - 11.L*mByM)*LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_14PNTidalCoeff(
      REAL8 mByM /**< ratio of object mass to total mass */
    )
/* literature: Agathos et al (arxiv 1503.0545) eq (5)
 * caveat: these are incomplete terms
 * conversion see XLALSimInspiralTaylorF2Phasing_13PNTidalCoeff above
 * --> completed by the terms given in equation (4) of :
 * Tatsuya Narikawa, Nami Uchikata, Takahiro Tanaka,
 * "Gravitational-wave constraints on the GWTC-2 events by measuring
 * the tidal deformability and the spin-induced quadrupole moment",
 * Phys. Rev. D 104, 084056 (2021), arXiv:2106.09193
 */
{
  REAL8 mByM3 = mByM*mByM*mByM;
  REAL8 mByM4 = mByM3 * mByM;
  return - mByM4 * 5.L*(193986935.L/571536.L - 14415613.L/381024.L*mByM - 57859.L/378.L*mByM*mByM - 209495.L/1512.L*mByM3 + 965.L/54.L*mByM4 - 4.L*mByM4*mByM);
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_15PNTidalCoeff(
      REAL8 mByM /**< ratio of object mass to total mass */
    )
/* literature: Agathos et al (arxiv 1503.0545) eq (5)
 * conversion see XLALSimInspiralTaylorF2Phasing_13PNTidalCoeff above
 * --> corrected by the terms given in equation (4) of :
 * Tatsuya Narikawa, Nami Uchikata, Takahiro Tanaka,
 * "Gravitational-wave constraints on the GWTC-2 events by measuring
 * the tidal deformability and the spin-induced quadrupole moment",
 * Phys. Rev. D 104, 084056 (2021), arXiv:2106.09193
 */
{
  REAL8 mByM2 = mByM*mByM;
  REAL8 mByM3 = mByM2*mByM;
  REAL8 mByM4 = mByM3*mByM;
  return mByM4 * 1.L/28.L*LAL_PI*(27719.L - 22415.L*mByM + 7598.L*mByM2 - 10520.L*mByM3) ;
}

/* The phasing function for TaylorF2 frequency-domain waveform.
 * This function is tested in ../test/PNCoefficients.c for consistency
 * with the energy and flux in this file.
 */
static void UNUSED
XLALSimInspiralPNPhasing_F2(
	PNPhasingSeries *pfa, /**< \todo UNDOCUMENTED */
	const REAL8 m1, /**< Mass of body 1, in Msol */
	const REAL8 m2, /**< Mass of body 2, in Msol */
	const REAL8 chi1L, /**< Component of dimensionless spin 1 along Lhat */
	const REAL8 chi2L, /**< Component of dimensionless spin 2 along Lhat */
	const REAL8 chi1sq,/**< Magnitude of dimensionless spin 1 */
	const REAL8 chi2sq, /**< Magnitude of dimensionless spin 2 */
	const REAL8 chi1dotchi2, /**< Dot product of dimensionles spin 1 and spin 2 */
	LALDict *p /**< LAL dictionary containing accessory parameters */
	)
{
    const REAL8 mtot = m1 + m2;
    const REAL8 eta = m1*m2/mtot/mtot;
    const REAL8 m1M = m1/mtot;
    const REAL8 m2M = m2/mtot;

    const REAL8 pfaN = 3.L/(128.L * eta);

    memset(pfa, 0, sizeof(PNPhasingSeries));

    pfa->v[0] = 1.L;
    pfa->v[1] = 0.L;
    pfa->v[2] = XLALSimInspiralTaylorF2Phasing_2PNCoeff(eta);
    pfa->v[3] = XLALSimInspiralTaylorF2Phasing_3PNCoeff(eta);
    pfa->v[4] = XLALSimInspiralTaylorF2Phasing_4PNCoeff(eta);
    pfa->v[5] = XLALSimInspiralTaylorF2Phasing_5PNCoeff(eta);
    pfa->vlogv[5] = XLALSimInspiralTaylorF2Phasing_5PNLogCoeff(eta);
    pfa->v[6] =  XLALSimInspiralTaylorF2Phasing_6PNCoeff(eta);
    pfa->vlogv[6] = XLALSimInspiralTaylorF2Phasing_6PNLogCoeff(eta);
    pfa->v[7] = XLALSimInspiralTaylorF2Phasing_7PNCoeff(eta);

    /* modify the PN coefficients if a non null LALSimInspiralTestGRParam structure is passed */
    /* BEWARE: this is for the non-spinning case only!*/
    pfa->v[0]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi0(p));
    pfa->v[1] = XLALSimInspiralWaveformParamsLookupNonGRDChi1(p);
    pfa->v[2]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi2(p));
    pfa->v[3]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi3(p));
    pfa->v[4]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi4(p));
    pfa->v[5]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi5(p));
    pfa->vlogv[5]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi5L(p));
    pfa->v[6]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi6(p));
    pfa->vlogv[6]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi6L(p));
    pfa->v[7]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi7(p));

    const REAL8 qm_def1=1.+XLALSimInspiralWaveformParamsLookupdQuadMon1(p);
    const REAL8 qm_def2=1.+XLALSimInspiralWaveformParamsLookupdQuadMon2(p);

    switch( XLALSimInspiralWaveformParamsLookupPNSpinOrder(p) )
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
	    pfa->v[7] += XLALSimInspiralTaylorF2Phasing_7PNSOCoeff(m1M)*chi1L+XLALSimInspiralTaylorF2Phasing_7PNSOCoeff(m2M)*chi2L;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
            pfa->v[6] += XLALSimInspiralTaylorF2Phasing_6PNSOCoeff(m1M)*chi1L
	      + XLALSimInspiralTaylorF2Phasing_6PNSOCoeff(m2M)*chi2L
	      + XLALSimInspiralTaylorF2Phasing_6PNS1S2OCoeff(eta)*chi1L*chi2L
	      + (XLALSimInspiralTaylorF2Phasing_6PNQM2SCoeff(m1M)*qm_def1+XLALSimInspiralTaylorF2Phasing_6PNSelf2SCoeff(m1M))*chi1sq
	      + (XLALSimInspiralTaylorF2Phasing_6PNQM2SCoeff(m2M)*qm_def2+XLALSimInspiralTaylorF2Phasing_6PNSelf2SCoeff(m2M))*chi2sq;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            pfa->v[5] += XLALSimInspiralTaylorF2Phasing_5PNSOCoeff(m1M)*chi1L
	      + XLALSimInspiralTaylorF2Phasing_5PNSOCoeff(m2M)*chi2L;
            pfa->vlogv[5] += 3.*(XLALSimInspiralTaylorF2Phasing_5PNSOCoeff(m1M)*chi1L
				 + XLALSimInspiralTaylorF2Phasing_5PNSOCoeff(m2M)*chi2L);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
	    /* 2PN SS, QM, and self-spin */
            pfa->v[4] += XLALSimInspiralTaylorF2Phasing_4PNS1S2Coeff(eta)*chi1dotchi2+XLALSimInspiralTaylorF2Phasing_4PNS1S2OCoeff(eta)*chi1L*chi2L
	      + (XLALSimInspiralTaylorF2Phasing_4PNQM2SOCoeff(m1M)*qm_def1+XLALSimInspiralTaylorF2Phasing_4PNSelf2SOCoeff(m1M))*chi1L*chi1L
	      + (XLALSimInspiralTaylorF2Phasing_4PNQM2SOCoeff(m2M)*qm_def2+XLALSimInspiralTaylorF2Phasing_4PNSelf2SOCoeff(m2M))*chi2L*chi2L
	      + (XLALSimInspiralTaylorF2Phasing_4PNQM2SCoeff(m1M)*qm_def1+XLALSimInspiralTaylorF2Phasing_4PNSelf2SCoeff(m1M))*chi1sq
	      + (XLALSimInspiralTaylorF2Phasing_4PNQM2SCoeff(m2M)*qm_def2+XLALSimInspiralTaylorF2Phasing_4PNSelf2SCoeff(m2M))*chi2sq;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            pfa->v[3] += XLALSimInspiralTaylorF2Phasing_3PNSOCoeff(m1M)*chi1L+XLALSimInspiralTaylorF2Phasing_3PNSOCoeff(m2M)*chi2L;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %i\n",
			   __func__, XLALSimInspiralWaveformParamsLookupPNSpinOrder(p) );
            XLAL_ERROR_VOID(XLAL_EINVAL);
            break;
    }

    REAL8 lambda1=XLALSimInspiralWaveformParamsLookupTidalLambda1(p);
    REAL8 lambda2=XLALSimInspiralWaveformParamsLookupTidalLambda2(p);
    switch( XLALSimInspiralWaveformParamsLookupPNTidalOrder(p) )
    {
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_75PN:
            pfa->v[15] = (lambda1*XLALSimInspiralTaylorF2Phasing_15PNTidalCoeff(m1M) + lambda2*XLALSimInspiralTaylorF2Phasing_15PNTidalCoeff(m2M));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_DEFAULT:
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_7PN:
            pfa->v[14] = (lambda1*XLALSimInspiralTaylorF2Phasing_14PNTidalCoeff(m1M) + lambda2*XLALSimInspiralTaylorF2Phasing_14PNTidalCoeff(m2M));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_65PN:
            pfa->v[13] = (lambda1*XLALSimInspiralTaylorF2Phasing_13PNTidalCoeff(m1M) + lambda2*XLALSimInspiralTaylorF2Phasing_13PNTidalCoeff(m2M));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
            pfa->v[12] = (lambda1*XLALSimInspiralTaylorF2Phasing_12PNTidalCoeff(m1M) + lambda2*XLALSimInspiralTaylorF2Phasing_12PNTidalCoeff(m2M) );
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
            pfa->v[10] = ( lambda1*XLALSimInspiralTaylorF2Phasing_10PNTidalCoeff(m1M) + lambda2*XLALSimInspiralTaylorF2Phasing_10PNTidalCoeff(m2M) );
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid tidal PN order %i\n",
                           __func__, XLALSimInspiralWaveformParamsLookupPNTidalOrder(p) );
            XLAL_ERROR_VOID(XLAL_EINVAL);
    }


    /* At the very end, multiply everything in the series by pfaN */
    for(int ii = 0; ii <= PN_PHASING_SERIES_MAX_ORDER; ii++)
    {
        pfa->v[ii] *= pfaN;
        pfa->vlogv[ii] *= pfaN;
        pfa->vlogvsq[ii] *= pfaN;
    }
}

/**
 * Computes the PN Coefficients for using in the TaylorT2 phasing equation.
 *
 * Terms given in equation 3.8a of: Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 */

static REAL8 UNUSED
XLALSimInspiralTaylorT2Phasing_0PNCoeff(
	REAL8 eta)
{
	return -1./(32.*eta);
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Phasing_2PNCoeff(
	REAL8 eta)
{
	return 3.715/1.008 + 5.5/1.2 * eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Phasing_3PNCoeff(
	REAL8 UNUSED eta)
{
	return -10. * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Phasing_4PNCoeff(
	REAL8 eta)
{
	return 15.293365/1.016064 + 27.145/1.008 * eta + 30.85/1.44 * eta*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Phasing_5PNCoeff(
	REAL8 eta)
{
	return (386.45/6.72 - 65./8. * eta) * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Phasing_6PNCoeff(
	REAL8 eta)
{
	return 1234.8611926451/1.8776862720 - 160./3. * LAL_PI*LAL_PI - 171.2/2.1 * LAL_GAMMA
		+ (225.5/4.8 * LAL_PI*LAL_PI - 1573.7765635/1.2192768) * eta
		+ 76.055/6.912 * eta*eta - 127.825/5.184 * eta*eta*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Phasing_6PNLogCoeff(
	REAL8 UNUSED eta)
{
	return -85.6/2.1;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Phasing_7PNCoeff(
	REAL8 eta)
{
	return (77.096675/2.032128 + 37.8515/1.2096 * eta - 74.045/6.048 * eta*eta) * LAL_PI;
}

/*
 * TaylorT2 derivatives dt/dv
 */

/* The expression for dt/dv has an extra factor of M not implemented here */
static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_0PNCoeff(
    REAL8 eta)
{
    return 5./(32.*eta);
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_2PNCoeff(
    REAL8 eta)
{
    return 743./336. + 11.*eta/4.;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_3PNCoeff(
    REAL8 UNUSED eta)
{
    return -4.*LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNCoeff(
    REAL8 eta)
{
    return 3058673./1016064. + 5429.*eta/1008. + 617.*eta*eta/144.;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_5PNCoeff(
    REAL8 eta)
{
    return (-7729./672.+13.*eta/8.)*LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNCoeff(
    REAL8 eta)
{
    return -10817850546611./93884313600. + 32.*LAL_PI*LAL_PI/3.
            + 1712.*LAL_GAMMA/105.
            + (3147553127./12192768. - 451.*LAL_PI*LAL_PI/48.)*eta
            - 15211.*eta*eta/6912. + 25565.*eta*eta*eta/5184.
            + 856.*log(16.)/105.;
}

/* The convention here is that this is the coefficient in front of v^6 log(v)
 * in the dt/dv expansion, NOT the one proportional to v^6 log(16 v^2).
 * Hence the term above containing log(16).
 */
static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNLogCoeff(
    REAL8 UNUSED eta)
{
    return 1712./105.;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_7PNCoeff(
    REAL8 eta)
{
    return LAL_PI*(-15419335./1016064. -75703.*eta/6048. + 14809.*eta*eta/3024);
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_3PNSOCoeff(
    REAL8 mByM)
{
    return 19./6. + 25./mByM/4.;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_5PNSOCoeff(
    REAL8 mByM)
{
    return -17.*mByM*mByM/4. + 5.*mByM + 1249./36. + 8349./mByM/224.;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNSOCoeff(
    REAL8 mByM)
{
    return LAL_PI*( -13. - 149./mByM/6.);
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_7PNSOCoeff(
    REAL8 mByM)
{
    const REAL8 mByMsq = mByM*mByM;
    return 1069.*mByMsq*mByMsq/288. - 1741.*mByMsq*mByM/192. + 176383.*mByMsq/12096. + 707767.*mByM/3456. + 133100377./6096384. + 34195607./mByM/193536.;
}

/* At 2 PN there are several spin^2 terms; see arXiv:astro-ph/0504538
 * The dt/dv spin^2 term at 2 PN is just -sigma (Eq. 9b-9d)
 * The terms 4PNSS and 4PNSSL are spin1-spin1 terms.
 */
static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNS1S2Coeff(
    REAL8 eta)
{
    return -79./8./eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNS1vS2vCoeff(
    REAL8 eta)
{
    return -71./24./eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNS1nS2nCoeff(
    REAL8 eta)
{
    return 33./eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNS1S2CoeffAvg(
    REAL8 eta)
{
    return 247./48./eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNS1OS2OCoeffAvg(
    REAL8 eta)
{
    return -721./48./eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNS1S1Coeff(
    REAL8 mByM)
{
    return -1./16./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNS1vS1vCoeff(
        REAL8 mByM)
{
    return -1./48./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNS1S1CoeffAvg(
    REAL8 mByM)
{
    return - 7./96./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNS1OS1OCoeffAvg(
        REAL8 mByM)
{
    return 1./96./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNQMS1S1Coeff(
        REAL8 mByM)
{
        return -5./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNQMS1nS1nCoeff(
	REAL8 mByM)
{
	return 16.5/mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNQMS1vS1vCoeff(
	REAL8 mByM)
{
	return -1.5/mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNQMS1S1CoeffAvg(
        REAL8 mByM)
{
        return 2.5/mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNQMS1OS1OCoeffAvg(
	REAL8 mByM)
{
	return -7.5/mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNS1S2Coeff(
    REAL8 eta)
{
    return -98.173/(1.344*eta) - 46.1/2.4;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNS1OS2OCoeff(
    REAL8 eta)
{
    return 140.3/(2.4*eta) + 6.7/1.8;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNS1nS2nCoeff(
    REAL8 eta)
{
    return 373./(3.*eta) + 93.25/1.68;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNS1vS2vCoeff(
    REAL8 eta)
{
    return 140.699/(4.032*eta) - 32.11/2.52;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNS1S1Coeff(
    REAL8 mByM)
{
    return 27.565/(2.688*mByM*mByM) - 21.3/(1.6*mByM) - 17.3/4.8;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNS1OS1OCoeff(
    REAL8 mByM)
{
    return 32.5/(1.6*mByM*mByM) + 107./6./mByM + 67./36.;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNS1nS1nCoeff(
    REAL8 mByM)
{
    return -59.37/(1.12*mByM*mByM) + 103.7/5.6/mByM + 11.17/3.36;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNS1vS1vCoeff(
    REAL8 mByM)
{
    return  8.5/(806.4*mByM*mByM) + 6.485/1.008/mByM +2.9/50.4;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNQMS1S1Coeff(
    REAL8 mByM)
{
    return  -94.07/(3.36*mByM*mByM) -58.7/2.4/mByM + 6.;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNQMS1nS1nCoeff(
    REAL8 mByM)
{
    return  563.9/(8.4*mByM*mByM) + 157.45/1.68/mByM - 171./7;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNQMS1vS1vCoeff(
    REAL8 mByM)
{
    return  56.65/(3.36*mByM*mByM) - 170.9/8.4/mByM + 45./7.;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNS1S2CoeffAvg(
    REAL8 eta)
{
    return 52.973/(8.064*eta) + 3.13/1.44;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNS1OS2OCoeffAvg(
    REAL8 eta)
{
    return -170.603/(8.064*eta) - 25.43/1.44;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNS1S1CoeffAvg(
    REAL8 mByM)
{
    return -37.427/(2.304*mByM*mByM) - 2.41/(2.88*mByM) - 5.51/2.88;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNS1OS1OCoeffAvg(
    REAL8 mByM)
{
    return 75.4979/(1.6128*mByM*mByM) + 15.43/2.88/mByM + 4.9/28.8;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNQMS1S1CoeffAvg(
    REAL8 mByM)
{
    return  94.07/(6.72*mByM*mByM) + 58.7/4.8/mByM -3.;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_6PNQMS1OS1OCoeffAvg(
    REAL8 mByM)
{
    return  -94.07/(2.24*mByM*mByM) -58.7/1.6/mByM + 9.;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_10PNTidalCoeff(
        REAL8 mByM)
{
        return 6.*mByM*mByM*mByM*mByM * (-12.+11.*mByM);
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_12PNTidalCoeff(
        REAL8 mByM)
{
        return mByM*mByM*mByM*mByM * (-3179/8. + 919/8.*mByM + 1143/4.*mByM*mByM - 65./2.*mByM*mByM*mByM);
}

/*
 * Tidal correction coefficients to Phasing
 */

static REAL8 UNUSED
XLALSimInspiralTaylorT2Phasing_10PNTidalCoeff(
	REAL8 chi,
	REAL8 lambda)
{
	return lambda * (-66.*chi + 72.) * chi*chi*chi*chi;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Phasing_12PNTidalCoeff(
	REAL8 eta,
	REAL8 chi,
	REAL8 lambda)
{
	return lambda * chi*chi*chi*chi * ( -1497.5*chi/5.6 - 225.5*eta*chi/1.4 + 1589.5/5.6
		+ 259.5*eta/1.4 + 398.5*chi*chi/2.8 - 965.*chi*chi*chi/7.);
}

/**
 * Computes the PN Coefficients for using in the TaylorT2 timing equation.
 *
 * Terms given in equation 3.8b of: Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 */

static REAL8 UNUSED
XLALSimInspiralTaylorT2Timing_0PNCoeff(
	REAL8 totalmass,
	REAL8 eta)
{
	totalmass *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert totalmass from kilograms to seconds */
	return -5.*totalmass/(256.*eta);
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Timing_2PNCoeff(
	REAL8 eta)
{
	return 7.43/2.52 + 11./3. * eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Timing_3PNCoeff(
	REAL8 UNUSED eta)
{
	return -32./5. * LAL_PI;;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Timing_4PNCoeff(
	REAL8 eta)
{
	return 30.58673/5.08032 + 54.29/5.04*eta + 61.7/7.2*eta*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Timing_5PNCoeff(
	REAL8 eta)
{
	return -(77.29/2.52 -13./3.*eta) * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Timing_6PNCoeff(
	REAL8 eta)
{
	return -1005.2469856691/2.3471078400 + 128./3. * LAL_PI*LAL_PI + 68.48/1.05 * LAL_GAMMA
		+ (3147.553127/3.048192 - 45.1/1.2 * LAL_PI*LAL_PI) * eta
		- 15.211/1.728 * eta*eta + 25.565/1.296 * eta*eta*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Timing_6PNLogCoeff(
	REAL8 UNUSED eta)
{
	return 34.24/1.05;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Timing_7PNCoeff(
	REAL8 eta)
{
	return (-154.19335/1.27008 - 757.03/7.56 * eta + 148.09/3.78 * eta*eta) * LAL_PI;
}

/*
 * Tidal correction coefficients to Timing
 */

static REAL8 UNUSED
XLALSimInspiralTaylorT2Timing_10PNTidalCoeff(
	REAL8 chi,
	REAL8 lambda)
{
	return lambda * (-264.*chi + 288.) * chi*chi*chi*chi;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2Timing_12PNTidalCoeff(
	REAL8 eta,
	REAL8 chi,
	REAL8 lambda)
{
	return lambda * chi*chi*chi*chi * (-2995.*chi/4. - 451.*eta*chi + 3179./4. + 519.*eta
		+ (797.*chi*chi)/2. - 386.*chi*chi*chi);
}



/**
 * Computes the PN Coefficients for using in the TaylorT3 phasing equation.
 *
 * Terms given in equation 3.10a of: Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 */


static REAL8 UNUSED
XLALSimInspiralTaylorT3Phasing_0PNCoeff(
	REAL8 eta)
{
	return -1./eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Phasing_2PNCoeff(
	REAL8 eta)
{
	return 3.715/8.064 + 5.5/9.6 * eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Phasing_3PNCoeff(
	REAL8 UNUSED eta)
{
	return -3./4. * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Phasing_4PNCoeff(
	REAL8 eta)
{
	return 9.275495/14.450688 + 2.84875/2.58048 * eta + 1.855/2.048 * eta*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Phasing_5PNCoeff(
	REAL8 eta)
{
	return (3.8645/2.1504 - 6.5/25.6 * eta) * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Phasing_6PNCoeff(
	REAL8 eta)
{
	return 83.1032450749357/5.7682522275840 - 5.3/4.0 * LAL_PI*LAL_PI - 10.7/5.6 * LAL_GAMMA
		+ (-126.510089885/4.161798144 + 2.255/2.048 * LAL_PI*LAL_PI) * eta
		+ 1.54565/18.35008 * eta*eta - 1.179625/1.769472 * eta*eta*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Phasing_6PNLogCoeff(
	REAL8 UNUSED eta)
{
	return -10.7/5.6;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Phasing_7PNCoeff(
	REAL8 eta)
{
	return (1.88516689/1.73408256 + 4.88825/5.16096 * eta - 1.41769/5.16096 * eta*eta) * LAL_PI;
}

/*
 * Tidal correction coefficients to Phasing
 */

static REAL8 UNUSED
XLALSimInspiralTaylorT3Phasing_10PNTidalCoeff(
	REAL8 chi,
	REAL8 lambda)
{
	return lambda * (-3.3*chi/51.2 + 9./128.) * chi*chi*chi*chi;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Phasing_12PNTidalCoeff(
	REAL8 eta,
	REAL8 chi,
	REAL8 lambda)
{
	return lambda * chi*chi*chi*chi * (-1.30715*chi/13.76256
		- 8.745*eta*chi/114.688 + 2.3325/22.9376 + 4.905*eta/57.344
		+ 3.985*chi*chi/114.688 - 9.65*chi*chi*chi/286.72);
}



/**
 * Computes the PN Coefficients for using in the TaylorT3 frequency equation.
 *
 * Terms given in equation 3.10b of: Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 */

static REAL8 UNUSED
XLALSimInspiralTaylorT3Frequency_0PNCoeff(
	REAL8 totalmass)
{
	totalmass *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert totalmass from kilograms to seconds */
	return 1. / (8. * LAL_PI * totalmass);
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Frequency_2PNCoeff(
	REAL8 eta)
{
	return 7.43/26.88 + 1.1/3.2 * eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Frequency_3PNCoeff(
	REAL8 UNUSED eta)
{
	return -3./10. * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Frequency_4PNCoeff(
	REAL8 eta)
{
	return 1.855099/14.450688 + 5.6975/25.8048 * eta + 3.71/20.48 * eta*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Frequency_5PNCoeff(
	REAL8 eta)
{
	return (-7.729/21.504 + 1.3/25.6 * eta) * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Frequency_6PNCoeff(
	REAL8 eta)
{
	return -7.20817631400877/2.88412611379200 + 5.3/20.0 * LAL_PI*LAL_PI + 1.07/2.80 * LAL_GAMMA
		+ (25.302017977/4.161798144 - 4.51/20.48 * LAL_PI*LAL_PI) * eta
		- 3.0913/183.5008 * eta*eta + 2.35925/17.69472 * eta*eta*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Frequency_6PNLogCoeff(
	REAL8 UNUSED eta)
{
	return 1.07/2.80;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Frequency_7PNCoeff(
	REAL8 eta)
{
	return (-1.88516689/4.33520640 - 9.7765/25.8048 * eta + 1.41769/12.90240 * eta*eta) * LAL_PI;
}

/*
 * Tidal correction coefficients to Frequency
 */

static REAL8 UNUSED
XLALSimInspiralTaylorT3Frequency_10PNTidalCoeff(
	REAL8 chi,
	REAL8 lambda)
{
	return lambda * (-9.9*chi/102.4 + 2.7/25.6) * chi*chi*chi*chi;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT3Frequency_12PNTidalCoeff(
	REAL8 eta,
	REAL8 chi,
	REAL8 lambda)
{
	return lambda * chi*chi*chi*chi * (-8.579*chi/65.536 - 1.947*eta*chi/16.384
		+ 1.8453/13.1072 + 4.329*eta/32.768 + 2.391*chi*chi/65.536
		- 5.79*chi*chi*chi/163.84);
}

/**
 * Computes the PN Coefficients for using in the TaylorT4 frequency equation.
 *
 * Spin-less terms given in equation 3.6 of: Alessandra Buonanno, Bala R Iyer,
 * Evan Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 * ( with (domega/dt)/omega = 3 (dv/dt)/v ),
 * spin-orbit terms in eq. 4.10a of S.Marsat, A.Bohe, L.Blanchet, A.Buonanno
 * Class.Quant.Grav. 31 (2014) 025023, arXiv:1307.6793,
 * spin-spin terms in eqs. 5.15-17 of E.Racine, A.Buonanno, L.Kidder PRD80(2009)
 * 044010, arXiv:0812.4413
 *** UNREVIEWED ****
 */

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_0PNCoeff(
	REAL8 eta)
{
	return 96. / 5. * eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_2PNCoeff(
	REAL8 eta)
{
	return ( -(1.0/336.0) * (743.0 + 924.0*eta) );
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_3PNCoeff(
	REAL8 UNUSED eta)
{
	return 4.0 * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_3PNSOCoeff(
	REAL8 mByM)
{
	return - 19./6. - 25./4./mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNCoeff(
	REAL8 eta)
{
	return (34103. + 122949.*eta + 59472.*eta*eta)/18144.;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNS1S2Coeff(
	REAL8 eta)
{
       return 79. / (8. * eta);
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNS1nS2nCoeff(
	REAL8 eta)
{
	return -33. / eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNS1vS2vCoeff(
	REAL8 eta)
{
	return 7.1 / 2.4 / eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNS1S2CoeffAvg(
	REAL8 eta)
{
	return - 247. / 48. / eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNS1OS2OCoeffAvg(
	REAL8 eta)
{
	return 721. / 48. / eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNS1S1Coeff(
        REAL8 mByM)
{
        return 1./(16.*mByM*mByM);
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNS1vS1vCoeff(
        REAL8 mByM)
{
        return 1./48./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNS1S1CoeffAvg(
        REAL8 mByM)
{
        return 7./96./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNS1OS1OCoeffAvg(
        REAL8 mByM)
{
        return -1./96./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNQMS1S1Coeff(
        REAL8 mByM)
{
        return 5./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNQMS1nS1nCoeff(
        REAL8 mByM)
{
        return -16.5/mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNQMS1vS1vCoeff(
        REAL8 mByM)
{
        return 1.5/mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNQMS1S1CoeffAvg(
        REAL8 mByM)
{
        return -2.5/mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNQMS1OS1OCoeffAvg(
        REAL8 mByM)
{
        return 7.5/mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_5PNCoeff(
	REAL8 eta)
{
	return ( -(1.0/672.0) * LAL_PI * (4159.0 + 15876.0*eta) );
	/* coefficient 15876 corrected (from 14532) according
	   to 2005 erratas for L. Blanchet, Phys. Rev. D 54, 1417 (1996)
	   (see Phys. Rev. D 71 129904 (E) (2005)) and L. Blanchet,
	   B. R. Iyer, and B. Joguet, Phys. Rev. D 65, 064005 (2002)
	   (see Phys. Rev. D 71 129903 (E) (2005)).
	   See errata for Arun et al., Phys. Rev. D 71, 084008
	   (2005) (see  Phys. Rev. D 72 069903 (E) (2005))
	   for corrected coefficients
	*/
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_5PNSOCoeff(
	REAL8 mByM)
{
	return  -809./(84.*mByM) + 13.795/1.008 - 527.*mByM/24. - 79.*mByM*mByM/6.;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNCoeff(
	REAL8 eta)
{
	return ( 16447.322263/139.7088 - 1712./105.
		* LAL_GAMMA - 561.98689/2.17728 * eta + LAL_PI * LAL_PI
		* (16./3. + 451./48. * eta) + 541./896. * eta * eta
		- 5605./2592. * eta * eta * eta - 856./105. * log(16.) );
}

/* The coefficient of the log is normalized for the argument of the log to be v=(M omega)^(1/3) */
static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNLogCoeff(
	REAL8 UNUSED eta)
{
	return -(1712.0/105.0);
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNSOCoeff(
	REAL8 mByM)
{
	return LAL_PI * ( -37./3. - 151./6./mByM );
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNS1S2Coeff(
	REAL8 eta)
{
	return 98.69/(3.36*eta) - 168.5/4.8;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNS1OS2OCoeff(
	REAL8 eta)
{
	return 237./(4.*eta) + 49./3.;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNS1nS2nCoeff(
	REAL8 eta)
{
	return 36.31/(1.68*eta) + 211.67/1.68;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNS1vS2vCoeff(
	REAL8 eta)
{
	return -230.3/(4.8*eta) - 3.557/1.008;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNS1S2CoeffAvg(
	REAL8 eta)
{
	return 108.79/(6.72*eta) + 75.25/2.88;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNS1OS2OCoeffAvg(
	REAL8 eta)
{
        return 162.25/(2.24*eta) - 129.31/2.88;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNS1S1Coeff(
	REAL8 mByM)
{
        return -33.7/(3.2*mByM*mByM) + 41.5/(3.2*mByM) + 37.9/9.6;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNS1OS1OCoeff(
	REAL8 mByM)
{
	return 75./(4.*mByM*mByM) + 87./(4.*mByM) + 49./6.;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNS1nS1nCoeff(
	REAL8 mByM)
{
	return 59.37/(1.12*mByM*mByM) - 103.7/(5.6*mByM) - 11.17/3.36;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNS1vS1vCoeff(
	REAL8 mByM)
{
	return -2.3/(22.4*mByM*mByM) -13.201/(2.016*mByM) + 1.15/20.16;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNS1S1CoeffAvg(
	REAL8 mByM)
{
        return 101.9/(6.4*mByM*mByM) + 2.51/(5.76*mByM) + 13.33/5.76;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNS1OS1OCoeffAvg(
	REAL8 mByM)
{
	return -49.3/(6.4*mByM*mByM) + 197.47/(5.76*mByM) + 56.45/5.76;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNQMS1S1Coeff(
	REAL8 mByM)
{
	return 6.59/(1.12*mByM*mByM) - 7.3/(2.4*mByM) + 21.5;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNQMS1nS1nCoeff(
	REAL8 mByM)
{
	return 19.63/(3.36*mByM*mByM) - 4.99/(1.68*mByM) - 185.7/2.8;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNQMS1vS1vCoeff(
	REAL8 mByM)
{
	return -39.47/(1.68*mByM*mByM) + 25.4/(2.1*mByM) + 5.1/2.8;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNQMS1S1CoeffAvg(
	REAL8 mByM)
{
	return -6.59/(2.24*mByM*mByM) + 7.3/(4.8*mByM) - 43./4.;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_6PNQMS1OS1OCoeffAvg(
	REAL8 mByM)
{
	return 19.77/(2.24*mByM*mByM) - 7.3/(1.6*mByM) + 129./4.;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_7PNCoeff(
	REAL8 eta)
{
	return (LAL_PI/12096.0) * (-13245.0 + 717350.0*eta + 731960.0*eta*eta);
	/* coefficients 717350 and 731960 corrected (from 661775 and 599156) according
	   to 2005 erratas for L. Blanchet, Phys. Rev. D 54, 1417 (1996)
	   (see Phys. Rev. D 71 129904 (E) (2005)) and L. Blanchet,
	   B. R. Iyer, and B. Joguet, Phys. Rev. D 65, 064005 (2002)
	   (see Phys. Rev. D 71 129903 (E) (2005)).
	   See errata for Arun et al., Phys. Rev. D 71, 084008
	   (2005) (see  Phys. Rev. D 72 069903 (E) (2005))
	   for corrected coefficients
	*/
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_7PNSOCoeff(
	REAL8 mByM)
{
       return  -1195.759 / 18.144 / mByM + 2694.373 / 18.144 - 432.2 / 2.1 * mByM + 1425.7 / 86.4 *mByM*mByM - 351.05 / 8.64 *mByM*mByM*mByM - 108.19 / 4.32 *mByM*mByM*mByM*mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_8PNSOCoeff(
	REAL8 mByM)
{
        return LAL_PI*(266.519/2.016 - 166.5/(2.8*mByM) - 828.43*mByM/6.72 -343.03*mByM*mByM/3.36);
        // see eq.(4.10a) of arXiv:1307.6793
}

/*
 * Tidal correction coefficients to DOmega/dt
 */

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_10PNTidalCoeff(
        REAL8 mByM)
{
        return 6.*mByM*mByM*mByM*mByM * (12.-11.*mByM);
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_12PNTidalCoeff(
        REAL8 mByM)
{
        return mByM*mByM*mByM*mByM * (4421./56. - 12263./56.*mByM + 1893./4.*mByM*mByM - 661./2.*mByM*mByM*mByM);
}

/*
 * For L see eq. 2.9 of arxiv:gr-qc/9506022
 */

static REAL8 UNUSED
XLALSimInspiralLN(REAL8 M,
		  REAL8 eta,
		  REAL8 v)
{
        return M*M*eta/v;
}
/* eq. 4.7 of http://arxiv.org/pdf/1212.5520.pdf */
static REAL8 UNUSED
XLALSimInspiralL_2PN(
        REAL8 eta)
{
        return 1.5 + eta/6.;
}

/* Orbital averaged from eq.(2.9c) of
 * \cite Kidder:1995zr, see also eq.(4.7) of \cite Bohe:2012mr.
 * Explicit formula can be found in eq.(9) of https://dcc.ligo.org/T1500554/public.
 */

static REAL8 UNUSED
XLALSimInspiralL_3PNSicoeff(
        REAL8 mByM)
{
        return -(1.+1./mByM);
}

static REAL8 UNUSED
XLALSimInspiralL_3PNSincoeff(
        REAL8 mByM)
{
        return 0.5*(1.+3./mByM);
}

static REAL8 UNUSED
XLALSimInspiralL_3PNSiLcoeff(
        REAL8 mByM)
{
        return 0.5*(1./3.-3./mByM);
}

static REAL8 UNUSED
XLALSimInspiralL_3PNSicoeffAvg(
        REAL8 mByM)
{
        return -0.75-0.25/mByM;
}

static REAL8 UNUSED
XLALSimInspiralL_3PNSiLcoeffAvg(
        REAL8 mByM)
{
        return -(1./3.+9./mByM)/4.;
}

static REAL8 UNUSED
XLALSimInspiralL_5PNSicoeff(
        REAL8 mByM)
{
  REAL8 eta=mByM*(1.-mByM);
  return -0.5*(5.+1./mByM)+eta/3.*(1.+4./mByM);
}

static REAL8 UNUSED
XLALSimInspiralL_5PNSincoeff(
        REAL8 mByM)
{
  REAL8 eta=mByM*(1.-mByM);
  return 3./8.*(3.+5./mByM)-7.*eta/4.*(1./6.+1./mByM);
}

static REAL8 UNUSED
XLALSimInspiralL_5PNSiLcoeff(
        REAL8 mByM)
{
  REAL8 eta=mByM*(1.-mByM);
  return -1./8.*(15.+17./mByM)+eta/4.*(-1.7/1.8+19./3./mByM);
}

/* eq. 4.7 of http://arxiv.org/pdf/1212.5520.pdf */
static REAL8 UNUSED
XLALSimInspiralL_4PN(
        REAL8 eta)
{
        return 27./8. - 19./8.*eta + eta*eta/24.;
}

static REAL8 UNUSED
XLALSimInspiralL_6PN(
        REAL8 eta)
{
  return 13.5/1.6 + (-68.89/1.44 + 4.1/2.4 * LAL_PI*LAL_PI)*eta + 3.1/2.4*eta*eta + 7./1296.*eta*eta*eta;
}

/*
 * dLh
 *
 * \f$d \hat{L_N}/d \hat{t} = M * d\hat{L_N} / dt = \Omega_L x \hat{L_N}\f$
 * This is Eq. (10) of gr-qc/0405090 ( times M b/c we use \f$\hat{t}\f$)
 */

static REAL8 UNUSED
XLALSimInspiralLDot_3PNSOCoeff(
        REAL8 mByM)
{
        return 0.5+1.5/mByM;
}

/* Using spin-self^2 derivatives at v^6 from
 * eq. A.2 of Blanchet et al. 1501.01529
 * and relating to LNh derivative through (4.25)
 * of arXiv:0812.4413.
 */
static REAL8 UNUSED
XLALSimInspiralLDot_4PNS1S2CoeffAvg(
        REAL8 eta)
{
        return -1.5/eta;
}

static REAL8 UNUSED
XLALSimInspiralLDot_4PNQMSSCoeff(
        REAL8 mByM)
{
        return -1.5/(mByM*mByM);
}

/* Using spin derivatives at v^7 from
 * eq. 7.8 of Blanchet et al. gr-qc/0605140
 * and relating to LNh derivative through (4.25)
 * of arXiv:0812.4413.
 */
static REAL8 UNUSED
XLALSimInspiralLDot_5PNSOCoeff(
        REAL8 mByM)
{
  return ( 9./8./mByM + 5./8 + 29./24.*mByM +mByM*mByM/24.);
}

// See (3.4) of arXiv:1212.5520
static REAL8 UNUSED
XLALSimInspiralLDot_7PNSOCoeff(
        REAL8 mByM)
{
        return -7.5/1.6 + 2.7/(1.6*mByM) + 53.*mByM/8. + 6.7*mByM*mByM/2.4 + 1.7*mByM*mByM*mByM/4.8 - mByM*mByM*mByM*mByM/48.;
}

/*
 * dS1
 * d S_1 / d \hat{t} = M * d S_1 / dt = \Omega_{S1,S2,LN,v} x S_1
 * However, that paper uses spin variables which are M^2 times our spins
 */

/* dS1, 1.5PN: eq. (8) of gr-qc/0405090.
 */
static REAL8 UNUSED
XLALSimInspiralSpinDot_3PNCoeff(
	REAL8 mByM)
{
	return 3./2. -mByM - mByM*mByM/2.;
}

/* S1S2 contribution
 * see. eq. A.2 of arXiv:1501.01529
 */
static const REAL8 UNUSED
XLALSimInspiralSpinDot_4PNS2Coeff=-1.;

static const REAL8 UNUSED
XLALSimInspiralSpinDot_4PNS2nCoeff=3.;

static const REAL8 UNUSED
XLALSimInspiralSpinDot_4PNS2CoeffAvg=0.5;

static const REAL8 UNUSED
XLALSimInspiralSpinDot_4PNS2OCoeffAvg=-1.5;

/* S1S1 contribution
 * again eq. A.2 of arXiv:1501.01529
 */
static REAL8 UNUSED
XLALSimInspiralSpinDot_4PNQMSOCoeffAvg(
	REAL8 mByM)
{
	return 1.5 * (1. - 1./mByM);
}

static REAL8 UNUSED
XLALSimInspiralSpinDot_4PNQMSnCoeff(
	REAL8 mByM)
{
	return 3. * (1./mByM - 1.);
}

/* dS1, 2.5PN
 * eq. 7.8 of Blanchet et al. gr-qc/0605140
 */
static REAL8 UNUSED
XLALSimInspiralSpinDot_5PNCoeff(
	REAL8 mByM)
{
	return 9./8. - mByM/2. + 7.*mByM*mByM/12. - 7.*mByM*mByM*mByM/6. - mByM*mByM*mByM*mByM/24.;
}

/* S1S2 contribution
 * again eq. A.2 of arXiv:1501.01529
 */

static REAL8 UNUSED
XLALSimInspiralSpinDot_6PNS2Coeff(
	REAL8 mByM)
{
        return -1.5 -mByM;
}

static REAL8 UNUSED
XLALSimInspiralSpinDot_6PNS1nCoeff(
	REAL8 mByM)
{
        return 3.5-3./mByM-.5*mByM*mByM;
}

static REAL8 UNUSED
XLALSimInspiralSpinDot_6PNS2nCoeff(
	REAL8 mByM)
{
        return 1.5 +2.*mByM+mByM*mByM;
}

static REAL8 UNUSED
XLALSimInspiralSpinDot_6PNS1vCoeff(
	REAL8 mByM)
{
        return 3. -1.5*mByM-1.5/mByM;
}

static REAL8 UNUSED
XLALSimInspiralSpinDot_6PNS2vCoeff(
	REAL8 mByM)
{
        return 1.5 +mByM;
}

static REAL8 UNUSED
XLALSimInspiralSpinDot_6PNS2CoeffAvg(
	REAL8 mByM)
{
  return XLALSimInspiralSpinDot_6PNS2Coeff(mByM)+0.5*(XLALSimInspiralSpinDot_6PNS2nCoeff(mByM)+XLALSimInspiralSpinDot_6PNS2vCoeff(mByM));
}

static REAL8 UNUSED
XLALSimInspiralSpinDot_6PNS1OCoeffAvg(
	REAL8 mByM)
{
  return -0.5*(XLALSimInspiralSpinDot_6PNS1nCoeff(mByM)+XLALSimInspiralSpinDot_6PNS1vCoeff(mByM));
}

static REAL8 UNUSED
XLALSimInspiralSpinDot_6PNS2OCoeffAvg(
	REAL8 mByM)
{
  return -0.5*(XLALSimInspiralSpinDot_6PNS2nCoeff(mByM)+XLALSimInspiralSpinDot_6PNS2vCoeff(mByM));
}

/* S1S1 contribution
 * again eq. A.2 of arXiv:1501.01529
 */

static REAL8 UNUSED
XLALSimInspiralSpinDot_6PNQMSnCoeff(
	REAL8 mByM)
{
        return 3. * (.5/mByM + 1. - mByM - .5*mByM*mByM);
}

static REAL8 UNUSED
XLALSimInspiralSpinDot_6PNQMSvCoeff(
	REAL8 mByM)
{
        return 3. * (1./mByM -1.);
}

static REAL8 UNUSED
XLALSimInspiralSpinDot_6PNQMSOCoeffAvg(
	REAL8 mByM)
{
  return -0.5*(XLALSimInspiralSpinDot_6PNQMSnCoeff(mByM)+XLALSimInspiralSpinDot_6PNQMSvCoeff(mByM));
}

/* dS1, 3.5PN
 * eq. 3.4 of Bohe' et al. arXiv:1212.5520
 */

static REAL8 UNUSED
XLALSimInspiralSpinDot_7PNCoeff(
	REAL8 mByM)
{
  return (mByM*mByM*mByM*mByM*mByM*mByM/48. - 3./8.*mByM*mByM*mByM*mByM*mByM - 3.9/1.6*mByM*mByM*mByM*mByM - 23./6.*mByM*mByM*mByM +18.1/1.6*mByM*mByM -51./8.*mByM + 2.7/1.6);
}

/**
 * Computes the PN Coefficients for using in the TaylorEt v(zeta) equation,
 * which is the square root of the x(zeta) equation.
 *
 * Terms given in equation 3.11 of: Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 */

static REAL8 UNUSED
XLALSimInspiralTaylorEtVOfZeta_2PNCoeff(
	REAL8 eta)
{
	return (3.0/4.0 + 1.0/12.0 * eta);
}

static REAL8 UNUSED
XLALSimInspiralTaylorEtVOfZeta_4PNCoeff(
	REAL8 eta)
{
	return (9.0/2.0 - 17.0/8.0 * eta + 1.0/18.0 * eta*eta);
}

static REAL8 UNUSED
XLALSimInspiralTaylorEtVOfZeta_6PNCoeff(
	REAL8 eta)
{
	return (40.5/1.6 + (20.5/9.6 * LAL_PI*LAL_PI - 479.5/7.2) * eta
		+ 5.5/6.4 * eta*eta + 3.5/129.6 * eta*eta*eta);
}


/**
 * Computes the PN Coefficients for using in the TaylorEt dPhase/dt equation.
 *
 * Terms given in equation 3.13a of: Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 */

static REAL8 UNUSED
XLALSimInspiralTaylorEtPhasing_0PNCoeff(
	REAL8 m)
{
	return 1.0/m;
}

static REAL8 UNUSED
XLALSimInspiralTaylorEtPhasing_2PNCoeff(
	REAL8 eta)
{
	return (9.0/8.0 + 1.0/8.0 * eta);
}

static REAL8 UNUSED
XLALSimInspiralTaylorEtPhasing_4PNCoeff(
	REAL8 eta)
{
	return (8.91/1.28 - 20.1/6.4 * eta + 1.1/12.8 * eta*eta);
}

static REAL8 UNUSED
XLALSimInspiralTaylorEtPhasing_6PNCoeff(
	REAL8 eta)
{
	return (41.445/1.024 - (309.715/3.072 - 20.5/6.4 * LAL_PI*LAL_PI) * eta
		+ 1.215/1.024 * eta*eta + 4.5/102.4 * eta*eta*eta);
}


/**
 * Computes the PN Coefficients for using in the TaylorEt dZeta/dt equation.
 *
 * Terms given in equation 3.13b of: Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 */

static REAL8 UNUSED
XLALSimInspiralTaylorEtZeta_0PNCoeff(
	REAL8 m,
	REAL8 eta)
{
	return 64.0 * eta / (5.0 * m);
}

static REAL8 UNUSED
XLALSimInspiralTaylorEtZeta_2PNCoeff(
	REAL8 eta)
{
	return (1.3/33.6 - 5.0/2.0 * eta);
}

static REAL8 UNUSED
XLALSimInspiralTaylorEtZeta_3PNCoeff(
	REAL8 UNUSED eta)
{
	return 4.0 * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorEtZeta_4PNCoeff(
	REAL8 eta)
{
	return (11.7857/1.8144 - 12.017/2.016 * eta + 5.0/2.0 * eta*eta);
}

static REAL8 UNUSED
XLALSimInspiralTaylorEtZeta_5PNCoeff(
	REAL8 eta)
{
	return (49.13/6.72 - 177.0/8.0 * eta) * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralTaylorEtZeta_6PNCoeff(
	REAL8 eta)
{
	return (379.99588601/2.79417600 + 16.0/3.0 * LAL_PI*LAL_PI - 17.12/1.05 * LAL_GAMMA
		+ (36.9/3.2 * LAL_PI*LAL_PI - 2486.1497/7.2576) * eta
		+ 48.8849/1.6128 * eta*eta - 8.5/6.4 * eta*eta*eta);
}

static REAL8 UNUSED
XLALSimInspiralTaylorEtZeta_6PNLogCoeff(
	REAL8 UNUSED eta)
{
	return -8.56/1.05;
}

static REAL8 UNUSED
XLALSimInspiralTaylorEtZeta_7PNCoeff(
	REAL8 eta)
{
	return (129.817/2.304 - 320.7739/4.8384 * eta + 61.3373/1.2096 * eta*eta) * LAL_PI;
}


/**
 * Computes the PN Coefficients for using in the TaylorF2Ecc equation.
 *
 * 3-dimensional REAL8 array eccPNCoeffs[ORDER][v_power][v0_power] are calculated,
 * where ORDER is relative PN order, v_power is power of v, and v0_power is power of v0.
 * Note that ORDER = v_power + v0_power.
 *
 * Terms given in equation 6.26 of: Blake Moore, Marc Favata,
 * K.G.Arun, and Chandra Kant Mishra, "Gravitational-wave phasing
 * for low-eccentricity inspiralling compact binaries to 3PN order",
 * Phys. Rev. D 93, 124061 (2016), arXiv:1605.00304
 */

static INT4 UNUSED
eccentricityPNCoeffs_F2(REAL8 eta, REAL8 eccPNCoeffs[LAL_MAX_ECC_PN_ORDER+1][LAL_MAX_ECC_PN_ORDER+1][LAL_MAX_ECC_PN_ORDER+1])
{
  INT4 ret = 0;
  memset(eccPNCoeffs, 0x00, (LAL_MAX_ECC_PN_ORDER+1)*(LAL_MAX_ECC_PN_ORDER+1)*(LAL_MAX_ECC_PN_ORDER+1)*sizeof(REAL8));
  eccPNCoeffs[0][0][0] = 1.0; // lowest order constant term

  eccPNCoeffs[2][2][0] = 29.9076223/8.1976608 + 18.766963/2.927736*eta; //v^2 term
  eccPNCoeffs[2][1][1] = 0.0; //v*v0 term
  eccPNCoeffs[2][0][2] = 2.833/1.008 - 19.7/3.6*eta; //v0^2 term

  eccPNCoeffs[3][3][0] = -28.19123/2.82600*LAL_PI; //v^3 term
  eccPNCoeffs[3][0][3] = 37.7/7.2*LAL_PI; //v0^3 term

  eccPNCoeffs[4][4][0] = 16.237683263/3.330429696 + 241.33060753/9.71375328*eta+156.2608261/6.9383952*eta*eta; //v^4 term
  eccPNCoeffs[4][2][2] = 84.7282939759/8.2632420864-7.18901219/3.68894736*eta-36.97091711/1.05398496*eta*eta; //v^2*v0^2 term
  eccPNCoeffs[4][0][4] = -1.193251/3.048192 - 66.317/9.072*eta +18.155/1.296*eta*eta;  //v0^4 term

  eccPNCoeffs[5][5][0] = -28.31492681/1.18395270*LAL_PI - 115.52066831/2.70617760*LAL_PI*eta; //v^5 term
  eccPNCoeffs[5][3][2] = -79.86575459/2.84860800*LAL_PI + 55.5367231/1.0173600*LAL_PI*eta; //v^3*v0^2 term
  eccPNCoeffs[5][2][3] = 112.751736071/5.902315776*LAL_PI + 70.75145051/2.10796992*LAL_PI*eta; //v^2*v0^3 term
  eccPNCoeffs[5][0][5] = 76.4881/9.0720*LAL_PI - 94.9457/2.2680*LAL_PI*eta;  //v0^5 term

  eccPNCoeffs[6][6][0] = -436.03153867072577087/1.32658535116800000 + 53.6803271/1.9782000*LAL_GAMMA + 157.22503703/3.25555200*LAL_PI*LAL_PI
                         +(2991.72861614477/6.89135247360 - 15.075413/1.446912*LAL_PI*LAL_PI)*eta
                         +345.5209264991/4.1019955200*eta*eta + 506.12671711/8.78999040*eta*eta*eta
                         + 384.3505163/5.9346000*log(2.0) - 112.1397129/1.7584000*log(3.0); //v^6 term except log(16*v^2) term
  eccPNCoeffs[6][4][2] = 46.001356684079/3.357073133568 + 253.471410141755/5.874877983744*eta
                         - 169.3852244423/2.3313007872*eta*eta - 307.833827417/2.497822272*eta*eta*eta; //v^4*v0^2 term
  eccPNCoeffs[6][3][3] = -106.2809371/2.0347200*LAL_PI*LAL_PI; //v^3*v0^3 term
  eccPNCoeffs[6][2][4] = -3.56873002170973/2.49880440692736 - 260.399751935005/8.924301453312*eta
                         + 15.0484695827/3.5413894656*eta*eta + 340.714213265/3.794345856*eta*eta*eta; //v^2*v0^4 term
  eccPNCoeffs[6][0][6] = 265.31900578691/1.68991764480 - 33.17/1.26*LAL_GAMMA + 12.2833/1.0368*LAL_PI*LAL_PI
                         + (91.55185261/5.48674560 - 3.977/1.152*LAL_PI*LAL_PI)*eta - 5.732473/1.306368*eta*eta
                         - 30.90307/1.39968*eta*eta*eta + 87.419/1.890*log(2.0) - 260.01/5.60*log(3.0);  //v0^6 term except log(16*v0^2) term
  //printPNCoeffs_F2(eccPNCoeffs);
  return ret;
}

/**
 * Compute eccentric phase correction term using eccPNCeoffs[k][i][j]
 *
 */

static REAL8 UNUSED
eccentricityPhasing_F2(REAL8 v, REAL8 v0, REAL8 ecc, REAL8 eta, INT4 ecc_order)
{
  static REAL8 v0_power[LAL_MAX_ECC_PN_ORDER+1];
  /* following code is not efficient in memory usage, need to be improved later */
  static REAL8 eccPNCoeffs[LAL_MAX_ECC_PN_ORDER+1][LAL_MAX_ECC_PN_ORDER+1][LAL_MAX_ECC_PN_ORDER+1]; // we want to calculate just one time
  REAL8 v_power[LAL_MAX_ECC_PN_ORDER+1];
  REAL8 phasing = 0.0;
  REAL8 global_factor;
  v0_power[0] = 1.0;
  for(int i=1; i<=LAL_MAX_ECC_PN_ORDER; i++)
  {
    v0_power[i] = v0_power[i-1]*v0;
  }
  eccentricityPNCoeffs_F2(eta, eccPNCoeffs);
  //printPNCoeffs_F2(eccPNCoeffs);
  v_power[0] = 1.0;
  for(int i=1; i<=LAL_MAX_ECC_PN_ORDER; i++)
  {
    v_power[i] = v_power[i-1]*v;
  }

  global_factor = -2.355/1.462*ecc*ecc*pow(v0/v, 19.0/3.0);
  global_factor *= (3.0/128.0/eta);  // overall factor except v^-5 in phase term, this is Newtonian phase term
  if(ecc_order == -1) {
    ecc_order = LAL_MAX_ECC_PN_ORDER;
  }
  if(ecc_order > LAL_MAX_ECC_PN_ORDER) {
    return XLAL_REAL8_FAIL_NAN;
  }

  REAL8 phaseOrder = 0;
  for(int i=0; i<=ecc_order; i++)
  {
    phaseOrder = 0;
    INT4 k = 0;
    for(int j=i; j>=0; j--)
    {
      k = i - j;
      if( j==6 )
      {
        phaseOrder += (eccPNCoeffs[i][j][k]+53.6803271/3.9564000*log(16.0*v_power[2]))*v_power[j]*v0_power[k];
        //phasing += (eccPNCoeffs[i][j][k]+53.6803271/3.9564000*log(16.0*v_power[2]))*v_power[j]*v0_power[k];
      }
      else if( k == 6 )
      {
        phaseOrder += (eccPNCoeffs[i][j][k] - 33.17/2.52*log(16.0*v0_power[2]))*v_power[j]*v0_power[k];
        //phasing += (eccPNCoeffs[i][j][k] - 33.17/2.52*log(16.0*v0_power[2]))*v_power[j]*v0_power[k];
      }
      else
      {
        phaseOrder += eccPNCoeffs[i][j][k]*v_power[j]*v0_power[k];
        //phasing += eccPNCoeffs[i][j][k]*v_power[j]*v0_power[k];
      }
    }
      phasing += phaseOrder;
      //ecc_phase_order[i] = phaseOrder*global_factor;
  }
  //fprintf(stdout, "======== DEBUG for eccentricity ================\n");
  //fprintf(stdout, "eccentricityPhasing_F2 phasing = %g, global_factor = %g, ecc_order = %d, ecc = %g\n", phasing, global_factor, ecc_order, ecc);
  return phasing*global_factor;
}
