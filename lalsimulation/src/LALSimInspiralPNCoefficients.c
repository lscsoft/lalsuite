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
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
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

static REAL8 UNUSED
XLALSimInspiralPNEnergy_3PNSOCoeff(
	REAL8 mByM)
{
	return 2. / 3. + 2. / mByM;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_4PNS1S2Coeff(
	REAL8 eta)
{
	return 1./eta;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_4PNS1S2OCoeff(
	REAL8 eta)
{
	return -3./eta;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_4PNQM2SCoeff(
	REAL8 mByM)
{
	return (1./mByM/mByM) / 2.;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_4PNQM2SOCoeff(
	REAL8 mByM)
{
	return -3. * (1./mByM/mByM) / 2.;
}

static REAL8 UNUSED
XLALSimInspiralPNEnergy_5PNSOCoeff(
	REAL8 mByM)
{
	return 5./3. + 3./mByM + 29.*mByM/9. + mByM*mByM/9.;
}

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
    return -103./48./eta;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNS1S2OCoeff(
    REAL8 eta)
{
    return 289./48./eta;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNQM2SCoeff(
    REAL8 mByM)
{
    return -1./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNQM2SOCoeff(
    REAL8 mByM)
{
    return 3./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNSelf2SCoeff(
    REAL8 mByM)
{
    return 7./96./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_4PNSelf2SOCoeff(
    REAL8 mByM)
{
    return -1./96./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_5PNCoeff(
	REAL8 eta)
{
	return -(81.91/6.72 + 58.3/2.4 * eta) * LAL_PI;
}

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

static REAL8 UNUSED
XLALSimInspiralPNFlux_6PNSOCoeff(
	REAL8 mByM)
{
	return LAL_PI*( -17./3. - 31./(6.*mByM) );
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_7PNCoeff(
	REAL8 eta)
{
	return -(162.85/5.04 - 214.745/1.728 * eta - 193.385/3.024 * eta*eta) * LAL_PI;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_7PNSOCoeff(
	REAL8 mByM)
{
        return (380.647/13.608) + 9535./(336.*mByM) - 40115.*mByM/756. + 3742.*mByM*mByM/63. - 35.*mByM*mByM*mByM/108. - 1117.*mByM*mByM*mByM*mByM/54.;
}

static REAL8 UNUSED
XLALSimInspiralPNFlux_8PNSOCoeff(
	REAL8 mByM)
{
        return LAL_PI * (125.47/2.52 - 71.63/(6.72*mByM) -3.137*mByM/2.016 - 212.41*mByM*mByM/3.36);
       // see eq. (4.9) of arXiv:1307.6793
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
	const REAL8 chi1sq, /**< Magnitude of dimensionless spin 1 */
	const REAL8 chi2sq, /**< Magnitude of dimensionless spin 2 */
	const REAL8 chi1dotchi2, /**< Dot product of dimensionles spin 1 and spin 2 */
	const REAL8 qm_def1, /**< Quadrupole deformation parameter of body 1 (dimensionless) */
	const REAL8 qm_def2, /**< Quadrupole deformation parameter of body 2 (dimensionless) */
	const LALSimInspiralSpinOrder spinO /**< Enums specifying spin order are in LALSimInspiralWaveformFlags.h */
	)
{
    const REAL8 mtot = m1 + m2;
    const REAL8 d = (m1 - m2) / (m1 + m2);
    const REAL8 eta = m1*m2/mtot/mtot;
    const REAL8 m1M = m1/mtot;
    const REAL8 m2M = m2/mtot;
    /* Use the spin-orbit variables from arXiv:1303.7412, Eq. 3.9
     * We write dSigmaL for their (\delta m/m) * \Sigma_\ell
     * There's a division by mtotal^2 in both the energy and flux terms
     * We just absorb the division by mtotal^2 into SL and dSigmaL
     */
    const REAL8 SL = m1M*m1M*chi1L + m2M*m2M*chi2L;
    const REAL8 dSigmaL = d*(m2M*chi2L - m1M*chi1L);

    const REAL8 pfaN = 3.L/(128.L * eta);

    memset(pfa, 0, sizeof(PNPhasingSeries));

    /* Non-spin phasing terms - see arXiv:0907.0700, Eq. 3.18 */
    pfa->v[0] = 1.L;
    pfa->v[2] = 5.L*(743.L/84.L + 11.L * eta)/9.L; 
    pfa->v[3] = -16.L*LAL_PI;
    pfa->v[4] = 5.L*(3058.673L/7.056L + 5429.L/7.L * eta
                     + 617.L * eta*eta)/72.L;
    pfa->v[5] = 5.L/9.L * (7729.L/84.L - 13.L * eta) * LAL_PI;
    pfa->vlogv[5] = 5.L/3.L * (7729.L/84.L - 13.L * eta) * LAL_PI;
    pfa->v[6] = (11583.231236531L/4.694215680L
                     - 640.L/3.L * LAL_PI * LAL_PI - 6848.L/21.L*LAL_GAMMA)
                 + eta * (-15737.765635L/3.048192L
                     + 2255./12. * LAL_PI * LAL_PI)
                 + eta*eta * 76055.L/1728.L
                 - eta*eta*eta * 127825.L/1296.L;
    pfa->v[6] += (-6848.L/21.L)*log(4.);
    pfa->vlogv[6] = -6848.L/21.L;
    pfa->v[7] = LAL_PI * ( 77096675.L/254016.L
                     + 378515.L/1512.L * eta - 74045.L/756.L * eta*eta);


    /* Compute 2.0PN SS, QM, and self-spin */
    // See Eq. (6.24) in arXiv:0810.5336
    // 9b,c,d in arXiv:astro-ph/0504538
    REAL8 pn_sigma = eta * (721.L/48.L*chi1L*chi2L - 247.L/48.L*chi1dotchi2);
    pn_sigma += (720.L*qm_def1 - 1.L)/96.0L * m1M * m1M * chi1L * chi1L;
    pn_sigma += (720.L*qm_def2 - 1.L)/96.0L * m2M * m2M * chi2L * chi2L;
    pn_sigma -= (240.L*qm_def1 - 7.L)/96.0L * m1M * m1M * chi1sq;
    pn_sigma -= (240.L*qm_def2 - 7.L)/96.0L * m2M * m2M * chi2sq;

    /* Spin-orbit terms - can be derived from arXiv:1303.7412, Eq. 3.15-16 */
    const REAL8 pn_gamma = (554345.L/1134.L + 110.L*eta/9.L)*SL + (13915.L/84.L - 10.L*eta/3.)*dSigmaL;
    switch( spinO )
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
            pfa->v[7] += (-8980424995.L/762048.L + 6586595.L*eta/756.L - 305.L*eta*eta/36.L)*SL - (170978035.L/48384.L - 2876425.L*eta/672.L - 4735.L*eta*eta/144.L) * dSigmaL;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
            pfa->v[6] += LAL_PI * (3760.L*SL + 1490.L*dSigmaL)/3.L;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            pfa->v[5] += -1.L * pn_gamma;
            pfa->vlogv[5] += -3.L * pn_gamma;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            pfa->v[4] += -10.L * pn_sigma;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            pfa->v[3] += 188.L*SL/3.L + 25.L*dSigmaL;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %i\n",
                    __func__, spinO );
            XLAL_ERROR_VOID(XLAL_EINVAL);
            break;
    }

    /* At the very end, multiply everything in the series by pfaN */
    for(int ii = 0; ii <= PN_PHASING_SERIES_MAX_ORDER; ii++)
    {
        pfa->v[ii] *= pfaN;
        pfa->vlogv[ii] *= pfaN;
        pfa->vlogvsq[ii] *= pfaN;
    }
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
 * The terms 4PNSS and 4PNSSL are spin1-spin1 terms, and also depend
 * on qm_def, the dimensionless quadrupole-monopole deformation
 * parameter (1 for BH, larger for NS
 */
static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNS1S2Coeff(
    REAL8 eta)
{
    return 247./48./eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNS1S2LCoeff(
    REAL8 eta)
{
    return -721./48./eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNSSCoeff(
    REAL8 mByM, REAL8 qm_def)
{
    return (5.*qm_def/2. - 7./96.)/mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT2dtdv_4PNSSLCoeff(
    REAL8 mByM, REAL8 qm_def)
{
    return (-15.*qm_def/2. + 1./96.)/mByM/mByM;
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
	return - 247. / 48. / eta;
}

static REAL8 UNUSED XLALSimInspiralTaylorT4wdot_4PNS1S2OCoeff(
	REAL8 eta)
{
	return 721. / 48. / eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNSelfSSCoeff(
	REAL8 mByM)
{
	return 7./96./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNSelfSSOCoeff(
	REAL8 mByM)
{
	return -1./96./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNQMCoeff(
	REAL8 mByM)
{
	return -5./2./mByM/mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorT4wdot_4PNQMSOCoeff(
	REAL8 mByM)
{
	return 15./2./mByM/mByM;
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

static REAL8 UNUSED
XLALSimInspiralSpinDot_3PNCoeff(
	REAL8 mByM)
{
	return 3./2. -mByM - mByM*mByM/2.;
}

static const REAL8 UNUSED
XLALSimInspiralSpinDot_4PNCoeffS1S2=0.5;

static const REAL8 UNUSED
XLALSimInspiralSpinDot_4PNCoeffLS1LS2=-1.5;

static REAL8 UNUSED
XLALSimInspiralSpinDot_4PNCoeffLSOSself(
	REAL8 mByM)
{
	return 1.5 * (1./mByM - 1.);
}

static REAL8 UNUSED
XLALSimInspiralSpinDot_5PNCoeff(
	REAL8 mByM)
{
	return 9./8. - mByM/2. + 7.*mByM*mByM/12. - 7.*mByM*mByM*mByM/6. - mByM*mByM*mByM*mByM/24.;
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
