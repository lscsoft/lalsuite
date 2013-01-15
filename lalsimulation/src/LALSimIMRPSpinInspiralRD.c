/*
 * Copyright (C) 2011 Riccardo Sturani, John Veitch, Drew Keppel
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
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/LALConstants.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/SphericalHarmonics.h>
#include <math.h>
#include <lal/LALAdaptiveRungeKutta4.h>
#include <lal/LALSimInspiral.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define minIntLen        8

/* use error codes above 1024 to avoid conflicts with GSL */
#define LALPSIRDPN_TEST_ENERGY		1025
#define LALPSIRDPN_TEST_OMEGADOT	1026
#define LALPSIRDPN_TEST_OMEGANAN	1028
#define LALPSIRDPN_TEST_OMEGAMATCH      1029
#define LALPSIRDPN_TEST_OMEGANONPOS     1031

typedef struct LALPSpinInspiralRDstructparams {
  REAL8 dt;
  REAL8 eta;                  ///< symmetric mass ratio
  REAL8 dm;                   ///< \f$m_1-m_2\f$
  REAL8 m1m2;                 ///< \f$m_1/m_2\f$
  REAL8 m2m1;                 ///< \f$m_2/m_1\f$
  REAL8 m1m;
  REAL8 m2m;
  REAL8 m1msq;
  REAL8 m2msq;
  REAL8 m;
  REAL8 wdotorb[8];           ///< Coefficients of the analytic PN expansion of \f$ \dot\omega_orb\f$
  REAL8 wdotorblog;           ///< Log coefficient of the PN expansion of of \f$\dot\omega_orb\f$
  REAL8 wdotspin15S1LNh;
  REAL8 wdotspin15S2LNh;
  REAL8 wdotspin20S1S2;
  REAL8 wdotspin20S1S1;       ///< Coeff. of the \f$s_1s_1\f$ cntrb. to \f$\dot\omega\f$
  REAL8 wdotspin20S1S2LNh;
  REAL8 wdotspin25S1LNh;
  REAL8 wdotspin25S2LNh;      ///< Coeff. of the \f$s_2\cdot \hat L_N\f$ cntrb. to \f$\dot\omega\f$
  REAL8 wdotspin30S1LNh;
  REAL8 wdotspin30S2LNh;
  REAL8 S1dot15;
  REAL8 S2dot15;
  REAL8 Sdot20;
  REAL8 Sdot20S;
  REAL8 S1dot25;
  REAL8 S2dot25;
  REAL8 LNhdot15;
  REAL8 epnorb[4];           ///< Coefficients of the PN expansion of the energy
  REAL8 epnspin15S1dotLNh;   ///< Coeff. of the \f$S_1\cdot L\f$ term in energy
  REAL8 epnspin15S2dotLNh;   ///< Coeff. of the \f$S_2\cdot L\f$ term in energy
  REAL8 epnspin20S1S2;       ///< Coeff. of the \f$S_1\cdot S_2\f$ term in energy
  REAL8 epnspin20S1S2dotLNh; ///< Coeff. of the \f$S_{1,2}\cdot L\f$ term in energy
  REAL8 epnspin20S1S1;       ///< Coeff. of the \f$S_1\cdot S_1\f$ term in energy
  REAL8 epnspin20S1S1dotLNh;
  REAL8 epnspin20S2S2;       ///< Coeff. of the \f$S_2\cdot S_2\f$ term in energy
  REAL8 epnspin20S2S2dotLNh;
  REAL8 epnspin25S1dotLNh;
  REAL8 epnspin25S2dotLNh;
  REAL8 lengths;
  REAL8 omOffset;
  REAL8 polarization;
  int length;
  UINT4 inspiralOnly;
} LALPSpinInspiralRDparams;

static REAL8 OmMatch(REAL8 LNhS1, REAL8 LNhS2, REAL8 S1S1, REAL8 S1S2, REAL8 S2S2) {

  const REAL8 omM       = 0.0555;
  const REAL8 omMsz12   =    9.97e-4;
  const REAL8 omMs1d2   =  -2.032e-3;
  const REAL8 omMssq    =   5.629e-3;
  const REAL8 omMsz1d2  =   8.646e-3;
  const REAL8 omMszsq   =  -5.909e-3;
  const REAL8 omM3s1d2  =   1.801e-3;
  const REAL8 omM3ssq   = -1.4059e-2;
  const REAL8 omM3sz1d2 =  1.5483e-2;
  const REAL8 omM3szsq  =   8.922e-3;

  return omM + /*6.05e-3 * sqrtOneMinus4Eta +*/
    omMsz12   * (LNhS1 + LNhS2) +
    omMs1d2   * (S1S2) +
    omMssq    * (S1S1 + S2S2) +
    omMsz1d2  * (LNhS1 * LNhS2) +
    omMszsq   * (LNhS1 * LNhS1 + LNhS2 * LNhS2) +
    omM3s1d2  * (LNhS1 + LNhS2) * (S1S2) +
    omM3ssq   * (LNhS1 + LNhS2) * (S1S1+S2S2) +
    omM3sz1d2 * (LNhS1 + LNhS2) * (LNhS1*LNhS2) +
    omM3szsq  * (LNhS1 + LNhS2) * (LNhS1*LNhS1+LNhS2*LNhS2);
}

static REAL8 fracRD(REAL8 LNhS1, REAL8 LNhS2, REAL8 S1S1, REAL8 S1S2, REAL8 S2S2) {

  const double frac0      = 0.840;
  const double fracsz12   = -2.145e-2;
  const double fracs1d2   = -4.421e-2;
  const double fracssq    = -2.643e-2;
  const double fracsz1d2  = -5.876e-2;
  const double fracszsq   = -2.215e-2;

  return frac0 +
    fracsz12   * (LNhS1 + LNhS2) +
    fracs1d2   * (S1S2) +
    fracssq    * (S1S1 + S2S2) +
    fracsz1d2  * (LNhS1 * LNhS2) +
    fracszsq   * (LNhS1 * LNhS1 + LNhS2 * LNhS2);
}


typedef struct LALPSpinInspiralPhenParsStruct {
  REAL8 endtime;
  REAL8 Psi;
  REAL8 alpha;
  REAL8 ci;
  REAL8 omega;
  REAL8 domega;
  REAL8 ddomega;
  REAL8 diota;
  REAL8 ddiota;
  REAL8 dalpha;
  REAL8 ddalpha;
  REAL8 energy;
  REAL8 LNhS1;
  REAL8 LNhS2;
  REAL8 S1S1;
  REAL8 S1S2;
  REAL8 S2S2;
  INT4  countback;
  INT4  intreturn;
} LALPSpinInspiralPhenPars;

typedef struct LALSpinInspiralAngleStruct {
  REAL8 ci2;
  REAL8 si2;
  REAL8 ci;
  REAL8 si;
  REAL8 c2i;
  REAL8 s2i;
  REAL8 c2i2;
  REAL8 s2i2;
  REAL8 c3i2;
  REAL8 s3i2;
  REAL8 c4i2;
  REAL8 s4i2;
  REAL8 c5i2;
  REAL8 s5i2;
  REAL8 c6i2;
  REAL8 s6i2;
  REAL8 c8i2;
  REAL8 s8i2;
  REAL8 cdi;
  REAL8 sdi;
} LALSpinInspiralAngle;


/******************** PN term calculators *******************************/

/* value of thetahat set according to
Blanchet et. al, Phys. Rev. Lett. 93, 091101 (2004) */
#define thetahat 1039.0/4620.0


static int XLALSimInspiralSpinTaylorCoeffs(REAL8 ST[9], /**< Output: spin coefficients array */
		      REAL8 eta /**< Symmetric mass ratio */ )
{
  if(!ST) XLAL_ERROR(XLAL_EFAULT);
  
   ST[0/*LAL_PNORDER_NEWTONIAN*/] = 1.0;
   ST[1/*LAL_PNORDER_HALF*/] = 0.0;
   ST[2/*LAL_PNORDER_ONE*/] = ( -(1.0/336.0) * (743.0 + 924.0*eta) );
   ST[3/*LAL_PNORDER_ONE_POINT_FIVE*/] = ( 4.0 * LAL_PI );
   ST[4/*LAL_PNORDER_TWO*/] =  ( (34103.0 + 122949.0*eta + 59472.0*eta*eta)/18144.0 );

   ST[5/*LAL_PNORDER_TWO_POINT_FIVE*/] = ( -(1.0/672.0) * LAL_PI * (4159.0 + 15876.0*eta) );
   /* coefficient 15876 corrected (from 14532) according
      to 2005 erratas for L. Blanchet, Phys. Rev. D 54, 1417 (1996)
      (see Phys. Rev. D 71 129904 (E) (2005)) and L. Blanchet,
      B. R. Iyer, and B. Joguet, Phys. Rev. D 65, 064005 (2002)
      (see Phys. Rev. D 71 129903 (E) (2005)).
      See errata for Arun et al., Phys. Rev. D 71, 084008
      (2005) (see  Phys. Rev. D 72 069903 (E) (2005))
      for corrected coefficients
   */

   /* both ak->ST[6] and [7] are stored for the threePN contribution */

   ST[6/*LAL_PNORDER_THREE*/] = ( (16447322263.0/139708800.0)
		- (1712.0/105.0)* LAL_GAMMA
		- (273811877.0/1088640.0)*eta - (88.0/3.0)*thetahat*eta
		+ (541.0/896.0)*eta*eta - (5605.0/2592.0)*eta*eta*eta
		+ (1.0/48.0) * LAL_PI*LAL_PI * (256.0 + 451.0*eta)
		- (856.0/105.0)*log(16.0) );
   ST[7/*LAL_PNORDER_THREE+1*/] = ( -(1712.0/315.0) );     /* extra 3PN component */
   /* sT[8] is the LAL_PNORDER_THREE_POINT_FIVE contribution */
   ST[8] = (LAL_PI/12096.0) * (-13245.0 + 717350.0*eta + 731960.0*eta*eta);
   /* coefficients 717350 and 731960 corrected (from 661775 and 599156) according
      to 2005 erratas for L. Blanchet, Phys. Rev. D 54, 1417 (1996)
      (see Phys. Rev. D 71 129904 (E) (2005)) and L. Blanchet,
      B. R. Iyer, and B. Joguet, Phys. Rev. D 65, 064005 (2002)
      (see Phys. Rev. D 71 129903 (E) (2005)).
      See errata for Arun et al., Phys. Rev. D 71, 084008
      (2005) (see  Phys. Rev. D 72 069903 (E) (2005))
      for corrected coefficients
   */
   return XLAL_SUCCESS;
}

#define lambda -11831./9240.

REAL8 ieta=1.0; /** Comparable mass limit */

REAL8 ETa3(REAL8 eta);
REAL8 ETa3(REAL8 eta)
{
 return (-675./64. + (209323./4032. - 205.*LAL_PI*LAL_PI/96.
            - 110./9. * lambda)*ieta*eta
            - 155./96. * ieta*eta*eta - 35./5184. * ieta*eta*eta*eta);
}
REAL8 ETa2(REAL8 eta);
REAL8 ETa2(REAL8 eta)
{
     return (-(27. - 19*ieta*eta + ieta*eta*eta/3.)/8.);
}
REAL8 ETa1(REAL8 eta);
REAL8 ETa1(REAL8 eta)
{
  return((9. + ieta*eta)/12.);
}
REAL8 ETaN(REAL8 eta);
REAL8 ETaN(REAL8 eta)
{
  return ( -eta/2.);
}

/***********************************************************************/



/**
 * Convenience function to set up LALPSpinInspiralRDparams struct
 */

static int XLALPSpinInspiralRDparamsSetup(
    LALPSpinInspiralRDparams *mparams,  /** Output: RDparams structure */
    UINT4 inspiralOnly,                 /** Only generate inspiral */
    REAL8 deltaT,                       /** sampling interval */
    REAL8 fLow,                         /** Starting frequency */
    REAL8 m1,                           /** Mass 1 */
    REAL8 m2,                           /** Mass 2 */
    LALSimInspiralSpinOrder spinO,      /** twice PN order of spin effects */
    UINT4 order                         /** twice PN Order in Phase */
    )
{
  REAL8 totalMass = m1+m2;
  REAL8 eta = m1*m2/(totalMass * totalMass);
  REAL8 chirpMass = pow(m1*m2,0.6)/pow(totalMass,0.2);
  REAL8 ST[9]; /* SpinTaylor terms */
  
  XLALSimInspiralSpinTaylorCoeffs(ST,eta);
  
  mparams->inspiralOnly = inspiralOnly;
  mparams->dt           = deltaT;
  mparams->lengths      = (5.0 / 256.0) / LAL_PI * pow(LAL_PI * chirpMass * LAL_MTSUN_SI * fLow,-5.0 / 3.0) / fLow;
  mparams->omOffset     = 0.006;

  /* setup coefficients for PN equations */
  mparams->m     = totalMass;
  mparams->m2m1  = m2 / m1;
  mparams->m1m2  = m1 / m2;
  mparams->m1m   = m1 / totalMass;
  mparams->m2m   = m2 / totalMass;
  mparams->m1msq = mparams->m1m * mparams->m1m;
  mparams->m2msq = mparams->m2m * mparams->m2m;
  mparams->dm    = (m1 - m2) / totalMass;
  mparams->eta = eta;

  switch (order) {

    case -1: // Use the highest PN order available. Move if higher terms added.
    case 7:
      mparams->wdotorb[7] = ST[8];

    case 6:
      mparams->epnorb[3] = ETa3(eta);
      mparams->wdotorb[6] = ST[6];
      mparams->wdotorblog = ST[7];
      mparams->wdotspin30S1LNh = -LAL_PI/3. * ( 188. - 151./2./mparams->m1m);
      mparams->wdotspin30S2LNh = -LAL_PI/3. * ( 188. + 151./2./mparams->m2m);

    case 5:
      mparams->wdotorb[5] = ST[5];
      mparams->epnspin25S1dotLNh = 8. - 31. / 9. * mparams->eta + (3. - 10. / 3. * mparams->eta) * mparams->m2m1;
      mparams->epnspin25S2dotLNh = 8. - 31. / 9. * mparams->eta + (3. - 10. / 3. * mparams->eta) * mparams->m1m2;
      mparams->wdotspin25S1LNh = -31319. / 1008. + 1159. / 24. * mparams->eta + (-809. / 84. + 281. / 8. * mparams->eta) * mparams->m2m1;
      mparams->wdotspin25S2LNh = -31319. / 1008. + 1159. / 24. * mparams->eta + (-809. / 84. + 281. / 8. * mparams->eta) * mparams->m1m2;
      mparams->S1dot25 = 0.5625 + 1.25 * mparams->eta - mparams->eta * mparams->eta / 24. + mparams->dm * (-0.5625 + 0.625 * mparams->eta);
      mparams->S2dot25 = 0.5625 + 1.25 * mparams->eta - mparams->eta * mparams->eta / 24. - mparams->dm * (-0.5625 + 0.625 * mparams->eta);

    case 4:
      mparams->epnorb[2] = ETa2(eta);
      mparams->wdotorb[4] = ST[4];
      mparams->wdotspin20S1S2 = -(1.0 / 48.0) / mparams->eta;
      mparams->wdotspin20S1S1 = 1. / 96.;
      mparams->Sdot20  = 0.5;
      mparams->Sdot20S = 0.5;
      mparams->epnspin20S1S2 = 1. / mparams->eta;
      mparams->epnspin20S1S2dotLNh = -3. / mparams->eta;
      mparams->epnspin20S1S1 = (1. + mparams->m2m1) * (1. + mparams->m2m1) / 2.;
      mparams->epnspin20S2S2 = (1. + mparams->m1m2) * (1. + mparams->m1m2) / 2.;
      mparams->epnspin20S1S1dotLNh = -3. * (1. + mparams->m2m1) * (1. + mparams->m2m1) / 2.;
      mparams->epnspin20S2S2dotLNh = -3. * (1. + mparams->m1m2) * (1. + mparams->m1m2) / 2.;

    case 3:
      mparams->wdotorb[3] = ST[3];
      mparams->epnspin15S1dotLNh = 8. / 3. + 2. * mparams->m2m1;
      mparams->epnspin15S2dotLNh = 8. / 3. + 2. * mparams->m1m2;
      mparams->wdotspin15S1LNh = -(113.0 + 75.0 * mparams->m2m1) / 12.0;
      mparams->wdotspin15S2LNh = -(113.0 + 75.0 * mparams->m1m2) / 12.0;
      mparams->LNhdot15 = 0.5;
      mparams->S1dot15 = (4.0 + 3.0 * mparams->m2m1) / 2.0 * mparams->eta;
      mparams->S2dot15 = (4.0 + 3.0 * mparams->m1m2) / 2.0 * mparams->eta;

    case 2:
      mparams->epnorb[1] = ETa1(eta);
      mparams->wdotorb[2] = ST[2];

    case 1:
      mparams->wdotorb[1] = ST[1];

    case 0:
      mparams->epnorb[0] = ETaN(eta);
      mparams->wdotorb[0] = ST[0];
      break;

    case 8:
      XLALPrintError("*** LALPhenSpinInspiralRD ERROR: PhenSpin approximant not available at pseudo4PN order\n");
			XLAL_ERROR(XLAL_EDOM);
      break;

    case 9:
      XLALPrintError("*** LALPhenSpinInspiralRD ERROR: NUM_ORDER not a valid PN order\n");
			XLAL_ERROR(XLAL_EDOM);
			break;

    default:
      XLALPrintError("*** LALPhenSpinInspiralRD ERROR: Impossible to create waveform with %d order\n",order);
			XLAL_ERROR(XLAL_EFAILED);
      break;
  }

  switch (spinO) {

    case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
    case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
    case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
      /*This kills all spin effects in the phase. Still there are spin effects
	in the waveform due to orbital plane precession*/      
      mparams->epnspin15S1dotLNh = 0.;
      mparams->epnspin15S2dotLNh = 0.;
      mparams->wdotspin15S1LNh   = 0.;
      mparams->wdotspin15S2LNh   = 0.;
      mparams->S1dot15           = 0.;
      mparams->S2dot15           = 0.;

    case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
      /* This keeps only the leading spin-orbit interactions*/
      mparams->wdotspin20S1S2      = 0.;
      mparams->epnspin20S1S2       = 0.;
      mparams->epnspin20S1S2dotLNh = 0.;

      mparams->wdotspin20S1S1 = 0.;
      mparams->epnspin20S1S1 = 0.;
      mparams->epnspin20S2S2 = 0.;
      mparams->Sdot20S = 0.;
      mparams->epnspin20S1S1 = 0.;
      mparams->epnspin20S2S2 = 0.;
      mparams->epnspin20S1S1dotLNh = 0.;
      mparams->epnspin20S2S2dotLNh = 0.;

    case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
      /* This kills all spin interaction intervening at 2.5PN order or higher*/
      mparams->epnspin25S1dotLNh   = 0.;
      mparams->epnspin25S2dotLNh   = 0.;
      mparams->wdotspin25S1LNh     = 0.;
      mparams->wdotspin25S2LNh     = 0.;
      mparams->S1dot25             = 0.;
      mparams->S2dot25             = 0.;

    case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
      /* This kills all spin interaction intervening at 3PN order or higher*/
      mparams->wdotspin30S1LNh     = 0.;
      mparams->wdotspin30S2LNh     = 0.;

    case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
    case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
    default:
      break;
  }
  return XLAL_SUCCESS;
}


static int XLALSpinInspiralDerivatives(
  UNUSED double t,
  const double values[],
  double dvalues[],
  void *mparams)
{
  REAL8 omega;                // time-derivative of the orbital phase
  REAL8 LNhx, LNhy, LNhz;     // orbital angular momentum unit vector
  REAL8 S1x, S1y, S1z;        // dimension-less spin variable S/M^2
  REAL8 S2x, S2y, S2z;
  REAL8 LNhS1, LNhS2;         // scalar products
  REAL8 domega;               // derivative of omega
  REAL8 dLNhx, dLNhy, dLNhz;  // derivatives of \f$\hat L_N\f$ components
  REAL8 dS1x, dS1y, dS1z;     // derivative of \f$S_i\f$
  REAL8 dS2x, dS2y, dS2z;
  REAL8 energy,energyold;

  /* auxiliary variables*/
  REAL8 S1S2, S1S1, S2S2;     // Scalar products
  REAL8 alphadotcosi;         // alpha is the right ascension of L, i(iota) the angle between L and J
  REAL8 v, v2, v3, v4, v5, v6, v7;
  REAL8 tmpx, tmpy, tmpz, cross1x, cross1y, cross1z, cross2x, cross2y, cross2z, LNhxy;

  LALPSpinInspiralRDparams *params = (LALPSpinInspiralRDparams *) mparams;

  /* --- computation start here --- */
  omega = values[1];

  LNhx = values[2];
  LNhy = values[3];
  LNhz = values[4];

  S1x = values[5];
  S1y = values[6];
  S1z = values[7];

  S2x = values[8];
  S2y = values[9];
  S2z = values[10];

  energyold = values[11];

  v = cbrt(omega);
  v2 = v * v;
  v3 = omega;
  v4 = omega * v;
  v5 = omega * v2;
  v6 = omega * omega;
  v7 = v6 * v;

  // Omega derivative without spin effects up to 3.5 PN
  // this formula does not include the 1.5PN shift mentioned in arXiv:0810.5336, five lines below (3.11)
  domega = params->wdotorb[0]
          + v * (params->wdotorb[1]
                 + v * (params->wdotorb[2]
                        + v * (params->wdotorb[3]
                               + v * (params->wdotorb[4]
                                      + v * (params->wdotorb[5]
                                             + v * (params->wdotorb[6] + params->wdotorblog * log(omega)
                                                    + v * params->wdotorb[7]))))));

  energy = (1. + v2 * (params->epnorb[1] +
                       v2 * (params->epnorb[2] +
                             v2 * params->epnorb[3])));

  // Adding spin effects
  // L dot S1,2
  LNhS1 = (LNhx * S1x + LNhy * S1y + LNhz * S1z);
  LNhS2 = (LNhx * S2x + LNhy * S2y + LNhz * S2z);

  // wdotspin15SiLNh = -1/12 (...)
  domega += v3 * (params->wdotspin15S1LNh * LNhS1 + params->wdotspin15S2LNh * LNhS2); // see e.g. Buonanno et al. gr-qc/0211087

  energy += v3 * (params->epnspin15S1dotLNh * LNhS1 + params->epnspin15S2dotLNh * LNhS2);  // see e.g. Blanchet et al. gr-qc/0605140

  // wdotspin20S1S1 = -1/48 eta
  S1S1 = (S1x * S1x + S1y * S1y + S1z * S1z);
  S2S2 = (S2x * S2x + S2y * S2y + S2z * S2z);
  S1S2 = (S1x * S2x + S1y * S2y + S1z * S2z);
  domega += params->wdotspin20S1S2 * v4 * (247.0 * S1S2 - 721.0 * LNhS1 * LNhS2);	// see e.g. Buonanno et al. arXiv:0810.5336
  domega += params->wdotspin20S1S1 * v4 * (719. * (LNhS1 * LNhS1 + LNhS2 * LNhS2) - 233. * (S1S1 + S2S2));   // see Racine et al. arXiv:0812.4413

  energy += v4 * (params->epnspin20S1S2 * S1S2 + params->epnspin20S1S2dotLNh * LNhS1 * LNhS2);    // see e.g. Buonanno et al. as above
  energy += v4 * (params->epnspin20S1S1 * S1S1 + params->epnspin20S2S2 * S2S2 + params->epnspin20S1S1dotLNh * LNhS1 * LNhS1 + params->epnspin20S2S2 * LNhS2 * LNhS2);	// see Racine et al. as above

  // wdotspin25SiLNh = see below
  domega += v5 * (params->wdotspin25S1LNh * LNhS1 + params->wdotspin25S2LNh * LNhS2);	//see (8.3) of Blanchet et al.
  energy += v5 * (params->epnspin25S1dotLNh * LNhS1 + params->epnspin25S2dotLNh * LNhS2);    //see (7.9) of Blanchet et al.

  domega += v6 * (params->wdotspin30S1LNh * LNhS1 + params->wdotspin30S2LNh * LNhS2); // see (6.5) of arXiv:1104.5659

  // Setting the right pre-factor
  domega *= 96. / 5. * params->eta * v5 * omega* omega;

  energy *= params->epnorb[0] * v2;

  /*Derivative of the angular momentum and spins */

  /* dS1, 1.5PN */
  /* S1dot15= (4+3m2/m1)/2 * eta */
  cross1x = (LNhy * S1z - LNhz * S1y);
  cross1y = (LNhz * S1x - LNhx * S1z);
  cross1z = (LNhx * S1y - LNhy * S1x);

  dS1x = params->S1dot15 * v5 * cross1x;
  dS1y = params->S1dot15 * v5 * cross1y;
  dS1z = params->S1dot15 * v5 * cross1z;

  /* dS1, 2PN */
  /* Sdot20= 0.5 */
  tmpx = S1z * S2y - S1y * S2z;
  tmpy = S1x * S2z - S1z * S2x;
  tmpz = S1y * S2x - S1x * S2y;

  // S1S2 contribution see. eq. 2.23 of arXiv:0812.4413
  dS1x += params->Sdot20 * v6 * (tmpx - 3. * LNhS2 * cross1x);
  dS1y += params->Sdot20 * v6 * (tmpy - 3. * LNhS2 * cross1y);
  dS1z += params->Sdot20 * v6 * (tmpz - 3. * LNhS2 * cross1z);
  // S1S1 contribution
  dS1x -= 3. * params->Sdot20S * v6 * LNhS1 * cross1x * (1. + params->m2m1) * params->m2m;
  dS1y -= 3. * params->Sdot20S * v6 * LNhS1 * cross1y * (1. + params->m2m1) * params->m2m;
  dS1z -= 3. * params->Sdot20S * v6 * LNhS1 * cross1z * (1. + params->m2m1) * params->m2m;

  // dS1, 2.5PN, eq. 7.8 of Blanchet et al. gr-qc/0605140
  // S1dot25= 9/8-eta/2.+eta+mparams->eta*29./24.+mparams->m1m2*(-9./8.+5./4.*mparams->eta)
  dS1x += params->S1dot25 * v7 * cross1x;
  dS1y += params->S1dot25 * v7 * cross1y;
  dS1z += params->S1dot25 * v7 * cross1z;

  /* dS2, 1.5PN */
  cross2x = (LNhy * S2z - LNhz * S2y);
  cross2y = (LNhz * S2x - LNhx * S2z);
  cross2z = (LNhx * S2y - LNhy * S2x);

  dS2x = params->S2dot15 * v5 * cross2x;
  dS2y = params->S2dot15 * v5 * cross2y;
  dS2z = params->S2dot15 * v5 * cross2z;

  /* dS2, 2PN */
  dS2x += params->Sdot20 * v6 * (-tmpx - 3.0 * LNhS1 * cross2x);
  dS2y += params->Sdot20 * v6 * (-tmpy - 3.0 * LNhS1 * cross2y);
  dS2z += params->Sdot20 * v6 * (-tmpz - 3.0 * LNhS1 * cross2z);
  // S2S2 contribution
  dS2x -= 3. * params->Sdot20S * v6 * LNhS2 * cross2x * (1. + params->m1m2) * params->m1m;
  dS2y -= 3. * params->Sdot20S * v6 * LNhS2 * cross2y * (1. + params->m1m2) * params->m1m;
  dS2z -= 3. * params->Sdot20S * v6 * LNhS2 * cross2z * (1. + params->m1m2) * params->m1m;

  // dS2, 2.5PN, eq. 7.8 of Blanchet et al. gr-qc/0605140
  dS2x += params->S2dot25 * v7 * cross2x;
  dS2y += params->S2dot25 * v7 * cross2y;
  dS2z += params->S2dot25 * v7 * cross2z;

  dLNhx = -(dS1x + dS2x) * v / params->eta;
  dLNhy = -(dS1y + dS2y) * v / params->eta;
  dLNhz = -(dS1z + dS2z) * v / params->eta;

  /* dphi */
  LNhxy = LNhx * LNhx + LNhy * LNhy;

  if (LNhxy > 0.0)
    alphadotcosi = LNhz * (LNhx * dLNhy - LNhy * dLNhx) / LNhxy;
  else
  {
    //XLALPrintWarning("*** LALPSpinInspiralRD WARNING ***: alphadot set to 0 by hand LNh:(%12.4e %12.4e %12.4e)\n",LNhx,LNhy,LNhz);
    alphadotcosi = 0.;
  }

  /* dvalues->data[0] is the phase derivative */
  /* omega is the derivative of the orbital phase omega \neq dvalues->data[0] */
  dvalues[0] = omega - alphadotcosi;
  dvalues[1] = domega;

  dvalues[2] = dLNhx;
  dvalues[3] = dLNhy;
  dvalues[4] = dLNhz;

  dvalues[5] = dS1x;
  dvalues[6] = dS1y;
  dvalues[7] = dS1z;

  dvalues[8] = dS2x;
  dvalues[9] = dS2y;
  dvalues[10] = dS2z;

  dvalues[11] = (energy-energyold)/params->dt*params->m*LAL_MTSUN_SI;

  return GSL_SUCCESS;
} /* end of XLALSpinInspiralDerivatives */

static int XLALGenerateWaveDerivative (
    REAL8Vector *dwave,
    REAL8Vector *wave,
    REAL8 dt
    )
{
  /* XLAL error handling */
  int errcode = XLAL_SUCCESS;

  /* For checking GSL return codes */
  INT4 gslStatus;

  UINT4 j;
  double *x, *y;
  double dy;
  gsl_interp_accel *acc;
  gsl_spline *spline;

  if (wave->length!=dwave->length)
    XLAL_ERROR( XLAL_EFUNC );

  /* Getting interpolation and derivatives of the waveform using gsl spline routine */
  /* Initialize arrays and supporting variables for gsl */

  x = (double *) LALMalloc(wave->length * sizeof(double));
  y = (double *) LALMalloc(wave->length * sizeof(double));

  if ( !x || !y )
  {
    if ( x ) LALFree (x);
    if ( y ) LALFree (y);
    XLAL_ERROR( XLAL_ENOMEM );
  }

  for (j = 0; j < wave->length; ++j)
  {
		x[j] = j;
		y[j] = wave->data[j];
  }

  XLAL_CALLGSL( acc = (gsl_interp_accel*) gsl_interp_accel_alloc() );
  XLAL_CALLGSL( spline = (gsl_spline*) gsl_spline_alloc(gsl_interp_cspline, wave->length) );
  if ( !acc || !spline )
  {
    if ( acc )    gsl_interp_accel_free(acc);
    if ( spline ) gsl_spline_free(spline);
    LALFree( x );
    LALFree( y );
    XLAL_ERROR( XLAL_ENOMEM );
  }

  /* Gall gsl spline interpolation */
  XLAL_CALLGSL( gslStatus = gsl_spline_init(spline, x, y, wave->length) );
  if ( gslStatus != GSL_SUCCESS )
  {
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    LALFree( x );
    LALFree( y );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Getting first and second order time derivatives from gsl interpolations */
  for (j = 0; j < wave->length; ++j)
  {
    XLAL_CALLGSL(gslStatus = gsl_spline_eval_deriv_e( spline, j, acc, &dy ) );
    if (gslStatus != GSL_SUCCESS )
    {
      gsl_spline_free(spline);
      gsl_interp_accel_free(acc);
      LALFree( x );
      LALFree( y );
      XLAL_ERROR( XLAL_EFUNC );
    }
    dwave->data[j]  = (REAL8)(dy / dt);
  }

  /* Free gsl variables */
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  LALFree(x);
  LALFree(y);
	
  return errcode;
}

static int XLALSpinInspiralTest(UNUSED double t, const double values[], double dvalues[], void *mparams) {
	
  LALPSpinInspiralRDparams *params = (LALPSpinInspiralRDparams *) mparams;
	
  REAL8 omega;
  REAL8 energy;
  REAL8 denergy;
	
  omega   = values[1];
  energy  = values[11];
  denergy = dvalues[11];
	
  REAL8 LNhS1=(values[2]*values[5]+values[3]*values[6]+values[4]*values[7])/params->m1msq;
  REAL8 LNhS2=(values[2]*values[8]+values[3]*values[9]+values[4]*values[10])/params->m2msq;
  REAL8 S1sq=(values[5]*values[5]+values[6]*values[6]+values[7]*values[7])/params->m1msq/params->m1msq;
  REAL8 S2sq=(values[8]*values[8]+values[9]*values[9]+values[10]*values[10])/params->m2msq/params->m2msq;
  REAL8 S1S2=(values[5]*values[8]+values[6]*values[9]+values[7]*values[10])/params->m1msq/params->m2msq;
	
  REAL8 omegaMatch=OmMatch(LNhS1,LNhS2,S1sq,S1S2,S2sq)+params->omOffset;
	
  if ( (energy > 0.0) || (( denergy > - 0.01*energy/params->dt*params->m*LAL_MTSUN_SI )&&(energy<0.) ) ) {
    /*energy increase*/
    /*XLALPrintWarning("** LALPSpinInspiralRD WARNING **: Energy increases dE %12.6e E %12.6e  m/M:(%12.4e, %12.4e)  om %12.6e vs. omM %12.6e\n",denergy, 0.01*energy, params->m1m, params->m2m, omega, omegaMatch);*/
    return LALPSIRDPN_TEST_ENERGY;
  }
  else if (omega < 0.0) {
    fprintf(stderr,"** LALPSpinInspiralRD WARNING **: Omega has become -ve, this should lead to nan's \n");
    return LALPSIRDPN_TEST_OMEGANONPOS;
  }
  else if (dvalues[1] < 0.0) {
    /* omegadot < 0 */
    return LALPSIRDPN_TEST_OMEGADOT;
  }
  else if (isnan(omega)) {
    /* omega is nan */
    return LALPSIRDPN_TEST_OMEGANAN;
  } 
  else if ((params->inspiralOnly!=1)&&(omega>omegaMatch)) {
    return LALPSIRDPN_TEST_OMEGAMATCH;
  }
  else
    return GSL_SUCCESS;
}

static int XLALSpinInspiralFillH2Modes(
    REAL8Vector* h2P2,
    REAL8Vector* h2M2,
    REAL8Vector* h2P1,
    REAL8Vector* h2M1,
    REAL8Vector* h20,
    UINT4 j,
    REAL4 amp,
    REAL4 v,
    REAL4 eta,
    REAL4 dm,
    REAL8 Psi,
    REAL8 alpha,
    LALSpinInspiralAngle *an
    )
{	
  REAL8 amp20 = amp * sqrt(1.5);
  REAL8 v2    = v*v;
  REAL8 omega = v*v*v;
  const REAL8 omegaC = 0.05;
  const REAL8 Afac   = .5;
  REAL8 damp  = omega > omegaC ? 1. : exp(-Afac*(1.-omegaC/omega)*(1.-omegaC/omega));
	
  h2P2->data[2 * j] = amp * ( ( 1. - damp * v2 / 42. * (107. - 55. * eta) )*
														 ( cos(2. * (Psi + alpha)) * an->c4i2 + cos(2. * (Psi - alpha)) * an->s4i2) +
														 v * dm / 3. * an->si *
														 ( cos(Psi - 2. * alpha) * an->s2i2 + cos(Psi + 2. * alpha) * an->c2i2 ) );
	
  h2M2->data[2 * j] = amp * ( ( 1. - damp * v2 / 42. * (107. - 55. * eta) )*
														 ( cos(2. * (Psi + alpha)) * an->c4i2 + cos(2. * (Psi - alpha)) * an->s4i2) -
														 v * dm / 3. * an->si *
														 ( cos(Psi - 2. * alpha) * an->s2i2 + cos(Psi + 2. * alpha) * an->c2i2) );
	
  h2P2->data[2 * j + 1] = amp * ( (1. - damp * v2 / 42. * (107. - 55. * eta) )*
                                 ( -sin(2. * (Psi + alpha)) * an->c4i2 + sin(2. * (Psi - alpha)) * an->s4i2) +
																 v * dm / 3. * an->si *
																 ( sin(Psi - 2. * alpha) * an->s2i2 - sin(Psi + 2. * alpha) * an->c2i2) );
	
  h2M2->data[2 * j + 1] = amp * ( (1. - damp * v2 / 42. * (107. - 55. * eta) )*
																 ( sin(2. * (Psi + alpha)) * an->c4i2 - sin(2. * (Psi - alpha)) * an->s4i2) +
																 v * dm / 3. * an->si *
																 ( sin(Psi - 2. * alpha) * an->s2i2 - sin(Psi + 2. * alpha) * an->c2i2) );
	
  h2P1->data[2 * j] = amp * (an->si * (1. - damp * v2 / 42. * (107. - 55. * eta) ) *
														 ( -cos(2. * Psi - alpha) * an->s2i2 + cos(2. * Psi + alpha) * an->c2i2) +
														 v * dm / 3. *
														 (-cos(Psi + alpha) * (an->ci + an->cdi)/2. -
															cos(Psi - alpha) * an->s2i2 * (1. + 2. * an->ci) ) );
	
  h2M1->data[2 * j] = amp * (an->si * (1. - damp * v2 / 42. * (107. - 55. * eta)) *
														 ( cos(2. * Psi - alpha) * an->s2i2 - cos(2. * Psi + alpha) * an->c2i2) +
														 v * dm / 3. *
														 (-cos(Psi + alpha) * (an->ci + an->cdi)/2. -
															cos(Psi - alpha) * an->s2i2 * (1. + 2. * an->ci) ) );
	
  h2P1->data[2 * j + 1] = amp * (an->si * (1. - damp * v2 / 42. * (107. - 55. * eta) ) *
																 ( -sin(2. * Psi - alpha ) * an->s2i2 - sin(2. * Psi + alpha) * an->c2i2) +
																 v * dm / 3. *
																 (sin(Psi + alpha) * (an->ci + an->cdi)/2. -
																	sin(Psi - alpha) * an->s2i2 * (1. + 2. * an->ci) ) );
	
  h2M1->data[2 * j + 1] = amp * (an->si * (1. - damp * v2 / 42. * (107. - 55. * eta)) *
																 ( -sin(2. * Psi - alpha) * an->s2i2 - sin(2. * Psi + alpha) * an->c2i2) -
																 v * dm / 3. *
																 (sin(Psi + alpha) * (an->ci + an->cdi) / 2. -
																	sin(Psi - alpha) * an->s2i2 * (1. + 2. * an->ci) ) );
	
  h20->data[2 * j] = amp20 * ( an->s2i * (1.- damp * v2/42. * (107. - 55.*eta) ) * cos(2. * Psi) );
	
  h20->data[2 * j + 1] = amp20 * ( v * dm / 3. * an->sdi * sin(Psi) );
	
  return XLAL_SUCCESS;
	
}

static int XLALSpinInspiralFillH3Modes(
    REAL8Vector* h3P3,
    REAL8Vector* h3M3,
    REAL8Vector* h3P2,
    REAL8Vector* h3M2,
    REAL8Vector* h3P1,
    REAL8Vector* h3M1,
    REAL8Vector* h30,
    UINT4 j,
    REAL4 amp,
    REAL4 v,
    REAL4 eta,
    REAL4 dm,
    REAL8 Psi,
    REAL8 alpha,
    LALSpinInspiralAngle *an
    )
{	
  REAL8 amp32 = amp * sqrt(1.5);
  REAL8 amp31 = amp * sqrt(0.15);
  REAL8 amp30 = amp / sqrt(5)/2.;
  REAL8 v2    = v*v;
	
  h3P3->data[2 * j] = amp * ( v * dm *
														 (-9. * cos(3. * (Psi - alpha)) * an->s6i2 -
															cos(  Psi - 3. * alpha) * an->s4i2 * an->c2i2 +
															cos(  Psi + 3. * alpha) * an->s2i2 * an->c4i2 +
															9. * cos(3. * (Psi + alpha)) * an->c6i2) +
														 v2 * 4. * an->si * (1. - 3. * eta) *
														 (-cos(2. * Psi - 3. * alpha) * an->s4i2 +
															cos(2. * Psi + 3. * alpha) * an->c4i2 ) );
	
  h3M3->data[2 * j] = amp * (-v * dm *
														 (-9. * cos(3. * (Psi - alpha)) * an->s6i2 -
															cos(  Psi - 3. * alpha) * an->s4i2 * an->c2i2 +
															cos(  Psi + 3. * alpha) * an->s2i2 * an->c4i2 +
															9. * cos(3. * (Psi + alpha)) * an->c6i2) +
														 v2 * 4. * an->si * (1. - 3. * eta) * 
														 (-cos(2. * Psi - 3. * alpha) * an->s4i2 +
															cos(2. * Psi + 3. * alpha) * an->c4i2 ) );
	
  h3P3->data[2 * j + 1] = amp * ( v * dm *
                                 (-9. * sin(3. * (Psi - alpha)) * an->s6i2 -
                                  sin(  Psi - 3. * alpha) * an->s4i2 * an->c2i2 -
                                  sin(  Psi + 3. * alpha) * an->s2i2 * an->c4i2 -
																	9. * sin(3. * (Psi + alpha)) * an->c6i2) +
																 v2 * 4. * an->si * (1. - 3. * eta) * 
																 (-sin(2. * Psi - 3. * alpha) * an->s4i2
																	-sin(2. * Psi + 3. * alpha) * an->c4i2 ) );
	
  h3M3->data[2 * j + 1] = amp * ( v * dm *
																 (-9. * sin(3. * (Psi - alpha)) * an->s6i2 -
																	sin(  Psi - 3. * alpha) * an->s4i2 * an->c2i2 -
																	sin(  Psi + 3. * alpha) * an->s2i2 * an->c4i2 -
																	9. * sin(3. * (Psi + alpha)) * an->c6i2) -
																 v2 * 4. * an->si * (1. - 3. * eta) *
																 (   - sin(2. * Psi - 3. * alpha) * an->s4i2
																	- sin(2. * Psi + 3. * alpha) * an->c4i2 ) );
	
  h3P2->data[2 * j] = amp32 * ( v * dm / 3. *
															 ( 27. * cos(3. * Psi - 2. * alpha) * an->si*an->s4i2 +
																27. * cos(3. * Psi + 2. * alpha) * an->si*an->c4i2 +
																cos( Psi + 2. * alpha) * an->c3i2 * (5.*an->si2-3.*an->si*an->ci2-3.*an->ci*an->si2) /2. +
																cos( Psi - 2. * alpha) * an->s3i2 * (5.*an->ci2+3.*an->ci*an->ci2-3.*an->si*an->si2) /2.) +
															 v2*(1./3.-eta) *
															 ( - 8.*an->c4i2*(3.*an->ci-2.)*cos(2.*(Psi+alpha)) +
																8.*an->s4i2*(3.*an->ci+2.)*cos(2.*(Psi-alpha)) ) );
	
  h3M2->data[2 * j] = amp32 * ( v * dm / 3. *
															 ( 27. * cos(3. * Psi - 2. * alpha) * an->si*an->s4i2 +
																27. * cos(3. * Psi + 2. * alpha) * an->si*an->c4i2 +
																cos( Psi + 2. * alpha) * an->c3i2 * (5.*an->si2-3.*an->si*an->ci2-3.*an->ci*an->si2) /2. +
																cos( Psi - 2. * alpha) * an->s3i2 * (5.*an->ci2+3.*an->ci*an->ci2-3.*an->si*an->si2) /2.) -
															 v2*(1./3.-eta) *
															 ( 8.*an->c4i2*(3.*an->ci-2.)*cos(2.*(Psi+alpha)) -
																8.*an->s4i2*(3.*an->ci+2.)*cos(2.*(Psi-alpha)) ) );
	
  h3P2->data[2 * j + 1 ] = amp32 * ( v * dm / 3. *
																		( 27. * sin(3. * Psi - 2. * alpha) * an->si*an->s4i2 -
																		 27. * cos(3. * Psi + 2. * alpha) * an->si*an->c4i2 -
																		 sin( Psi + 2. * alpha) * an->c3i2 * (5.*an->si2-3.*an->si*an->ci2-3.*an->ci*an->si2) /2. +
																		 sin( Psi - 2. * alpha) * an->s3i2 * (5.*an->ci2+3.*an->ci*an->ci2-3.*an->si*an->si2) /2.) +
																		v2*(1./3.-eta) *
																		( 8.*an->c4i2*(3.*an->ci-2.)*sin(2.*(Psi+alpha)) +
																		 8.*an->s4i2*(3.*an->ci+2.)*sin(2.*(Psi-alpha)) ) );
	
  h3M2->data[2 * j + 1 ] = amp32 * ( -v * dm / 3. *
																		( 27. * sin(3. * Psi - 2. * alpha) * an->si*an->s4i2 -
																		 27. * cos(3. * Psi + 2. * alpha) * an->si*an->c4i2 -
																		 sin( Psi + 2. * alpha) * an->c3i2 * (5.*an->si2-3.*an->si*an->ci2-3.*an->ci*an->si2) /2.+
																		 sin( Psi - 2. * alpha) * an->s3i2 * (5.*an->ci2+3.*an->ci*an->ci2-3.*an->si*an->si2) /2.)+
																		v2*(1./3.-eta) *
																		( 8.*an->c4i2*(3.*an->ci-2.)*sin(2.*(Psi+alpha)) +
																		 8.*an->s4i2*(3.*an->ci+2.)*sin(2.*(Psi-alpha)) ) );
	
  h3P1->data[2 * j] = amp31 * ( v * dm / 6. *
															 ( -135. * cos(3.*Psi - alpha) * an->s2i*an->s2i2 +
																135. * cos(3.*Psi + alpha)  * an->s2i*an->c2i2 +
																cos(Psi+alpha) * an->c2i2*(15.*an->cdi-20.*an->ci+13.)/2.-
																cos(Psi-alpha) * an->s2i2*(15.*an->cdi+20.*an->ci+13.)/2. )
															 -v2 * (1./3.-eta)*
															 ( 20.*an->c3i2*cos(2.*Psi+alpha)*(3.*(an->si2*an->ci+an->ci2*an->si)-5.*an->si2) +
																20.*an->s3i2*cos(2.*Psi-alpha)*(3.*(an->ci2*an->ci-an->si2*an->si)+5.*an->ci2) ) );
	
  h3M1->data[2 * j] = amp31 * (-v * dm / 6. *
															 ( -135. * cos(3.*Psi - alpha) * an->s2i*an->s2i2 +
																135. * cos(3.*Psi + alpha) * an->s2i*an->c2i2 +
																cos(Psi+alpha) * an->c2i2*(15.*an->cdi-20.*an->ci+13.)/2.-
																cos(Psi-alpha) * an->s2i2*(15.*an->cdi+20.*an->ci+13.)/2. )
															 -v2 * (1./3.-eta)*
															 ( 20.*an->c3i2*cos(2.*Psi+alpha)*(3.*(an->si2*an->ci+an->ci2*an->si)-5.*an->si2) +
																20.*an->s3i2*cos(2.*Psi-alpha)*(3.*(an->ci2*an->ci-an->si2*an->si)+5.*an->ci2) ) );
	
  h3P1->data[2 * j + 1] = amp31 * ( v * dm / 6. *
																	 ( -135. * sin(3.*Psi - alpha) * an->s2i*an->s2i2 -
																		135.* sin(3.*Psi + alpha) * an->s2i*an->c2i2 -
																		sin(Psi+alpha) * an->c2i2*(15.*an->cdi-20.*an->ci+13.)/2.-
																		sin(Psi-alpha) * an->s2i2*(15.*an->cdi+20.*an->ci+13.)/2. )
																	 +v2 * (1./3.-eta)*
																	 ( 20.*an->c3i2*sin(2.*Psi+alpha)*(3.*(an->si2*an->ci+an->ci2*an->si)-5.*an->si2)
																		-20.*an->s3i2*sin(2.*Psi-alpha)*(3.*(an->ci2*an->ci-an->si2*an->si)+5.*an->ci2) ) );
	
  h3M1->data[2 * j + 1] = amp31 * ( v * dm / 6. *
																	 ( -135. * sin(3.*Psi - alpha) *an->s2i*an->s2i2 -
																		135. * sin(3.*Psi + alpha) *an->s2i*an->c2i2 -
																		sin(Psi+alpha) * an->c2i2*(15.*an->cdi-20.*an->ci+13.)/2.-
																		sin(Psi-alpha) * an->s2i2*(15.*an->cdi+20.*an->ci+13.)/2. )
																	 -v2 * (1./3.-eta)*
																	 ( 20.*an->c3i2*sin(2.*Psi+alpha)*(3.*(an->si2*an->ci+an->ci2*an->si)-5.*an->si2)
																		-20.*an->s3i2*sin(2.*Psi-alpha)*(3.*(an->ci2*an->ci-an->si2*an->si)+5.*an->ci2) ) );
	
  h30->data[2 * j] = 0.;
	
  h30->data[2 * j + 1] =  amp30 * ( v * dm *
																	 ( cos(Psi) * an->si*(cos(2.*Psi)*(45.*an->s2i)-(25.*an->cdi-21.) ) )
																	 +v2 * (1.-3.*eta) *
																	 (80. * an->s2i*an->c2i*sin(2.*Psi) ) );
	
	return XLAL_SUCCESS;
	
}

static int XLALSpinInspiralFillH4Modes(
    REAL8Vector* h4P4,
    REAL8Vector* h4M4,
    REAL8Vector* h4P3,
    REAL8Vector* h4M3,
    REAL8Vector* h4P2,
    REAL8Vector* h4M2,
    REAL8Vector* h4P1,
    REAL8Vector* h4M1,
    REAL8Vector* h40,
    INT4  j,
    REAL8 amp44,
    UNUSED REAL8 v,
    REAL8 eta,
    UNUSED REAL8 dm,
    REAL8 Psi,
    REAL8 alpha,
    LALSpinInspiralAngle *an
    )
{		
  REAL8 amp43 = - amp44 * sqrt(2.);
  REAL8 amp42 = amp44 * sqrt(7.)/2.;
  REAL8 amp41 = amp44 * sqrt(3.5)/4.;
  REAL8 amp40 = amp44 * sqrt(17.5)/16.;
	
  h4P4->data[2 * j] = amp44 * (1. - 3.*eta) *
	( 4.* an->s8i2 * cos(4.*(Psi-alpha)) + cos(2.*Psi-4.*alpha) *an->s6i2*an->c2i2
	 + an->s2i2*an->c6i2* cos(2.*Psi+4.*alpha) + 4.*an->c8i2* cos(4.*(Psi+alpha)) );
	
  h4M4->data[2 * j] = amp44 * (1. - 3.*eta) *
	( 4.* an->s8i2 * cos(4.*(Psi-alpha)) + cos(2.*Psi-4.*alpha) *an->s6i2*an->c2i2
	 + an->s2i2*an->c6i2* cos(2.*Psi+4.*alpha) + 4.*an->c8i2* cos(4.*(Psi+alpha)) );
	
  h4P4->data[2 * j + 1] = amp44 * (1. - 3.*eta) *
	( 4.* an->s8i2 * sin(4.*(Psi-alpha)) + sin(2.*Psi-4.*alpha) *an->s6i2*an->c2i2
	 - an->s2i2*an->c6i2* sin(2.*Psi+4.*alpha) - 4.*an->c8i2* sin(4.*(Psi+alpha)) );
	
  h4M4->data[2 * j + 1] = - amp44 * (1. - 3.*eta) *
	( 4.* an->s8i2 * sin(4.*(Psi-alpha)) + sin(2*Psi-4.*alpha) *an->s6i2*an->c2i2
	 - an->s2i2*an->c6i2* sin(2.*Psi+4.*alpha) - 4.*an->c8i2* sin(4.*(Psi+alpha)) );
	
  h4P3->data[2 * j] = amp43 * (1. - 3.*eta) * an->si *
	( 4.*an->s6i2* cos(4.*Psi-3.*alpha) - 4.*an->c6i2* cos(4.*Psi+3.*alpha) -
	 an->s4i2*(an->ci+0.5)/2. * cos(2.*Psi-3.*alpha) - an->c4i2*(an->ci-0.5) * cos(2.*Psi+3.*alpha) );
	
  h4M3->data[2 * j] = - amp43 * (1. - 3.*eta) * an->si *
	( 4.*an->s6i2* cos(4.*Psi-3.*alpha) - 4.*an->c6i2* cos(4.*Psi+3.*alpha) -
	 an->s4i2*(an->ci+0.5)/2. * cos(2.*Psi-3.*alpha) - an->c4i2*(an->ci-0.5) * cos(2.*Psi+3.*alpha) );
	
  h4P3->data[2 * j + 1] = amp43 * (1. - 3.*eta) * an->si *
	( 4.*an->s6i2* sin(4.*Psi-3.*alpha) + 4.*an->c6i2* sin(4.*Psi+3.*alpha) -
	 an->s4i2*(an->ci+0.5)/2. * sin(2.*Psi-3.*alpha) + an->c4i2*(an->ci-0.5) * sin(2.*Psi+3.*alpha) );
	
  h4M3->data[2 * j + 1] = amp43 * (1. - 3.*eta) * an->si *
	( 4.*an->s6i2* sin(4.*Psi-3.*alpha) + 4.*an->c6i2* sin(4.*Psi+3.*alpha) -
	 an->s4i2*(an->ci+0.5)/2. * sin(2.*Psi-3.*alpha) + an->c4i2*(an->ci-0.5) * sin(2.*Psi+3.*alpha) );
	
  h4P2->data[2 * j] = amp42 * (1. - 3.*eta) *
	( 16.*an->s6i2*an->c2i2 * cos(4.*Psi-2.*alpha) + 16.*an->c6i2*an->s2i2 * cos(4.*Psi+2.*alpha)
	 - an->c4i2 * cos(2.*(Psi+alpha))*(an->cdi-2.*an->ci+9./7.)/2. - an->s4i2 * cos(2.*(Psi-alpha))*(an->cdi+2.*an->ci+9./7.)/2. );
	
  h4M2->data[2 * j] = amp42 * (1. - 3.*eta) *
	( 16.*an->s6i2*an->c2i2 * cos(4.*Psi-2.*alpha) + 16.*an->c6i2*an->s2i2 * cos(4.*Psi+2.*alpha)
	 - an->c4i2 * cos(2.*(Psi+alpha))*(an->cdi-2.*an->ci+9./7.)/2. - an->s4i2 * cos(2.*(Psi-alpha))*(an->cdi+2.*an->ci+9./7.)/2. );
	
  h4P2->data[2 * j + 1] = amp42 * (1. - 3.*eta) *
	( 16.*an->s6i2*an->c2i2 * sin(4.*Psi-2.*alpha) - 16.*an->c6i2*an->s2i2 * sin(4.*Psi+2.*alpha)
	 + an->c4i2 * sin(2.*(Psi+alpha))*(an->cdi-2.*an->ci+9./7.)/2. - an->s4i2 * sin(2.*(Psi-alpha))*(an->cdi+2.*an->ci+9./7.)/2. );
	
  h4M2->data[2 * j + 1] = -amp42 * (1. - 3.*eta) *
	( 16.*an->s6i2*an->c2i2 * sin(4.*Psi-2.*alpha) - 16.*an->c6i2*an->s2i2 * sin(4.*Psi+2.*alpha)
	 + an->c4i2 * sin(2.*(Psi+alpha))*(an->cdi-2.*an->ci+9./7.)/2. - an->s4i2 * sin(2.*(Psi-alpha))*(an->cdi+2.*an->ci+9./7.)/2. );
	
  h4P1->data[2 * j] = amp41 * (1. - 3.*eta) *
	( -64.*an->s5i2*an->c3i2 * cos(4.*Psi-alpha) +64.*an->s3i2*an->c5i2 * cos(4.*Psi+alpha) -
	 an->s3i2*cos(2.*Psi-alpha)*((an->cdi*an->ci2-an->sdi*an->si2)+2.*(an->ci2*an->ci-an->si2*an->si)+19./7.*an->ci2) +
	 an->c3i2*cos(2.*Psi+alpha)*((an->cdi*an->si2+an->sdi*an->ci2)-2.*(an->si*an->ci2+an->ci*an->si2)+19./7.*an->ci2) );
	
  h4M1->data[2 * j] = -amp41 * (1. - 3.*eta) *
	( -64*an->s5i2*an->c3i2 * cos(4.*Psi-alpha) +64.*an->s3i2*an->c5i2 * cos(4.*Psi+alpha) -
	 an->s3i2*cos(2.*Psi-alpha)*((an->cdi*an->ci2-an->sdi*an->si2)+2.*(an->ci2*an->ci-an->si2*an->si)+19./7.*an->ci2) +
	 an->c3i2*cos(2.*Psi+alpha)*((an->cdi*an->si2+an->sdi*an->ci2)-2.*(an->si*an->ci2+an->ci*an->si2)+19./7.*an->ci2) );
	
  h4P1->data[2 * j + 1] = amp41 * (1. - 3.*eta) *
	( -64.*an->s5i2*an->c3i2 * sin(4.*Psi-alpha) - 64.*an->s3i2*an->c5i2 * sin(4.*Psi+alpha) -
	 an->s3i2*sin(2.*Psi-alpha)*((an->cdi*an->ci2-an->sdi*an->si2)+2.*(an->ci2*an->ci-an->si2*an->si)+19./7.*an->ci2) -
	 an->c3i2*sin(2.*Psi+alpha)*((an->cdi*an->si2+an->sdi*an->ci2)-2.*(an->si*an->ci2+an->ci*an->si2)+19./7.*an->ci2) );
	
  h4M1->data[2 * j + 1] = amp41 * (1. - 3.*eta) *
	( -64.*an->s5i2*an->c3i2 * sin(4.*Psi-alpha) - 64.*an->s3i2*an->c5i2 * sin(4.*Psi+alpha) -
	 an->s3i2*sin(2.*Psi-alpha)*((an->cdi*an->ci2-an->sdi*an->si2)+2.*(an->ci2*an->ci-an->si2*an->si)+19./7.*an->ci2) -
	 an->c3i2*sin(2.*Psi+alpha)*((an->cdi*an->si2+an->sdi*an->ci2)-2.*(an->si*an->ci2+an->ci*an->si2)+19./7.*an->ci2) );
	
  h40->data[2 * j] = amp40 * (1.-3.*eta) * an->s2i * (8.*an->s2i*cos(4.*Psi) +
																										 cos(2.*Psi)*(an->cdi+5./7.) );
  h40->data[2 * j +1] = 0.;
	
  return XLAL_SUCCESS;
}


static int XLALSpinInspiralAdaptiveEngine(
    const UINT4 neqs, 
    const REAL8 yinit[],  
    REAL8 amp22ini, 
    LALPSpinInspiralRDparams *mparams,
    REAL8Vector* h2P2,
    REAL8Vector* h2M2,
    REAL8Vector* h2P1,
    REAL8Vector* h2M1,
    REAL8Vector* h20,
    REAL8Vector* h3P3,
    REAL8Vector* h3M3,
    REAL8Vector* h3P2,
    REAL8Vector* h3M2,
    REAL8Vector* h3P1,
    REAL8Vector* h3M1,
    REAL8Vector* h30,
    REAL8Vector* h4P4,
    REAL8Vector* h4M4,
    REAL8Vector* h4P3,
    REAL8Vector* h4M3,
    REAL8Vector* h4P2,
    REAL8Vector* h4M2,
    REAL8Vector* h4P1,
    REAL8Vector* h4M1,
    REAL8Vector* h40,
    REAL8Vector* freq,
    REAL8Vector* phase,
    LALPSpinInspiralPhenPars *phenPars 
    )
{
	
  INT4 j;
  INT4 k;
  INT4 kMatch=0;
  INT4 jMatch=0;
  INT4 Npoints=10;
  INT4 intlen;
  INT4 intreturn;
	
  LALSpinInspiralAngle trigAngle;

  REAL8Array *yout;
  ark4GSLIntegrator *integrator;
	
  REAL8 Psi;
  REAL8 alpha=0.;
  REAL8 alphaold;
  REAL8 v,v2;
  REAL8 dt;
  REAL8 Mass;
  REAL8 amp22;
  REAL8 amp33;
  REAL8 amp44;
	
  REAL8 LNhxy;
  REAL8 LNhS1;
  REAL8 LNhS2;
  REAL8 S1S1;
  REAL8 S1S2;
  REAL8 S2S2;
  REAL8 omegaMatch;
  REAL8 c1,c2;

  REAL8 *yin = (REAL8 *) LALMalloc(sizeof(REAL8) * neqs);
	
  /* allocate the integrator */
  integrator = XLALAdaptiveRungeKutta4Init(neqs,XLALSpinInspiralDerivatives,XLALSpinInspiralTest,1.0e-6,1.0e-6);
  if (!integrator) {
    fprintf(stderr,"**** LALPSpinInspiralRD ERROR ****: Cannot allocate adaptive integrator.\n");
    if (XLALClearErrno() == XLAL_ENOMEM)
      XLAL_ERROR( XLAL_ENOMEM );
    else
      XLAL_ERROR( XLAL_EDOM );
  }
	
  /* stop the integration only when the test is true */
  integrator->stopontestonly = 1;
	
  /* run the integration; note: time is measured in units of total mass */
	
  Mass = mparams->m * LAL_MTSUN_SI;
  dt   = mparams->dt;
	
  for (UINT4 jeq=0; jeq<neqs; jeq++) yin[jeq]=yinit[jeq];

  REAL8 S1x0=yinit[5];
  REAL8 S1y0=yinit[6];
  REAL8 S1z0=yinit[7];
  REAL8 S2x0=yinit[8];
  REAL8 S2y0=yinit[9];
  REAL8 S2z0=yinit[10];
	
  intlen = XLALAdaptiveRungeKutta4Hermite(integrator,(void *)mparams,yin,0.0,mparams->lengths/Mass,dt/Mass,&yout);

  intreturn = integrator->returncode;
  XLALAdaptiveRungeKutta4Free(integrator);
  if (intlen == XLAL_FAILURE)
  {
    XLALPrintError("Error in Adaptive Integrator\n");
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* End integration*/
	
  /* Start of the integration checks*/
  if (!intlen) {
    phenPars->intreturn=intreturn;
    if (XLALClearErrno() == XLAL_ENOMEM) {
      XLAL_ERROR(  XLAL_ENOMEM);
    } else {
      fprintf(stderr,"**** LALPSpinInspiralRD ERROR ****: integration failed with errorcode %d, integration length %d\n",intreturn,intlen);
      XLAL_ERROR( XLAL_EFAILED);
    }
  }
	
  /* if we have enough space, compute the waveform components; otherwise abort */
  if ( intlen >= mparams->length ) {
    fprintf(stderr,"**** LALPSpinInspiralRD ERROR ****: no space to write in waveforms: %d vs. %d\n",intlen,mparams->length);
    XLALPrintError("**** LALPSpinInspiralRD ERROR ****: error generating second derivatives\n");
    XLALPrintError("                     m:           : %12.5f  %12.5f\n",mparams->m1m*mparams->m,mparams->m2m*mparams->m);
    XLALPrintError("              S1:                 : %12.5f  %12.5f  %12.5f\n",S1x0,S1y0,S1z0);
    XLALPrintError("              S2:                 : %12.5f  %12.5f  %12.5f\n",S2x0,S2y0,S2z0);
    XLAL_ERROR(XLAL_ESIZE);
  }
	
  if ( intlen < minIntLen ) {
    fprintf(stderr,"**** LALPSpinInspiralRD ERROR ****: incorrect integration with length %d\n",intlen);
    XLAL_ERROR(XLAL_ESIZE);
  }
  /* End of integration checks*/
	
  REAL8 *Phi    = &yout->data[1*intlen];
  REAL8 *omega  = &yout->data[2*intlen];
  REAL8 *LNhx   = &yout->data[3*intlen];
  REAL8 *LNhy   = &yout->data[4*intlen];
  REAL8 *LNhz   = &yout->data[5*intlen];
  REAL8 *S1x    = &yout->data[6*intlen];
  REAL8 *S1y    = &yout->data[7*intlen];
  REAL8 *S1z    = &yout->data[8*intlen];
  REAL8 *S2x    = &yout->data[9*intlen];
  REAL8 *S2y    = &yout->data[10*intlen];
  REAL8 *S2z    = &yout->data[11*intlen];
  REAL8 *energy = &yout->data[12*intlen];
	
  mparams->polarization=2.*atan2(LNhy[1],LNhx[1])-atan2(LNhy[2],LNhx[2]);
	
  if (mparams->inspiralOnly!=1) {

    INT4 errcode=XLAL_SUCCESS;

    j=intlen;
	
    do {
      j--;
      LNhS1=(LNhx[j]*S1x[j]+LNhy[j]*S1y[j]+LNhz[j]*S1z[j])/mparams->m1msq;
      LNhS2=(LNhx[j]*S2x[j]+LNhy[j]*S2y[j]+LNhz[j]*S2z[j])/mparams->m2msq;
      S1S1=(S1x[j]*S1x[j]+S1y[j]*S1y[j]+S1z[j]*S1z[j])/mparams->m1msq/mparams->m1msq;
      S1S2=(S1x[j]*S2x[j]+S1y[j]*S2y[j]+S1z[j]*S2z[j])/mparams->m1msq/mparams->m2msq;
      S2S2=(S2x[j]*S2x[j]+S2y[j]*S2y[j]+S2z[j]*S2z[j])/mparams->m2msq/mparams->m2msq;
      omegaMatch=OmMatch(LNhS1,LNhS2,S1S1,S1S2,S2S2);
      if (omegaMatch>omega[j]) {
				if (omega[j-1]<omega[j]) jMatch=j;
				// The numerical integrator sometimes stops and stores twice the last
				// omega value, this 'if' instruction avoids keeping two identical 
				// values of omega at the end of the integration.
      }
    } while ((j>0)&&(jMatch==0));
		
    if (omegaMatch<omega[jMatch]) {
      fprintf(stderr,"*** LALPSpinInspiralRD ERROR ***: Impossible to attach phenom. part\n");
      XLAL_ERROR(XLAL_EFAILED);
    }

    // Data structure are copied into Npoints-long
    // REAL8Array for interpolation and derivative computation
    if ( ((INT4)Npoints) > intlen) Npoints = intlen;

    if ( (omega[jMatch+1]>omega[jMatch]) && ((jMatch+1)<intlen) )
      kMatch=Npoints-2;
    else
      kMatch=Npoints-1;
    //We keep until the point where omega > omegaMatch for better derivative
    //computation, but do the matching at the last point at which 
    // omega < omegaMatch
		
    REAL8Vector *omega_s   = XLALCreateREAL8Vector(Npoints);
    REAL8Vector *LNhx_s    = XLALCreateREAL8Vector(Npoints);
    REAL8Vector *LNhy_s    = XLALCreateREAL8Vector(Npoints);
    REAL8Vector *LNhz_s    = XLALCreateREAL8Vector(Npoints);
    REAL8Vector *alpha_s   = XLALCreateREAL8Vector(Npoints);
		
    REAL8Vector *domega    = XLALCreateREAL8Vector(Npoints);
    REAL8Vector *dLNhx     = XLALCreateREAL8Vector(Npoints);
    REAL8Vector *dLNhy     = XLALCreateREAL8Vector(Npoints);
    REAL8Vector *dLNhz     = XLALCreateREAL8Vector(Npoints);
    REAL8Vector *diota     = XLALCreateREAL8Vector(Npoints);
    REAL8Vector *dalpha    = XLALCreateREAL8Vector(Npoints);
		
    REAL8Vector *ddomega   = XLALCreateREAL8Vector(Npoints);
    REAL8Vector *ddiota    = XLALCreateREAL8Vector(Npoints);
    REAL8Vector *ddalpha   = XLALCreateREAL8Vector(Npoints);
		
    for (k=0;k<Npoints;k++) {
      j=k+jMatch-kMatch;
      omega_s->data[k]  = omega[j];
      LNhx_s->data[k]   = LNhx[j];
      LNhy_s->data[k]   = LNhy[j];
      LNhz_s->data[k]   = LNhz[j];
    }

    errcode  = XLALGenerateWaveDerivative(domega,omega_s,dt);
    errcode += XLALGenerateWaveDerivative(dLNhx,LNhx_s,dt);
    errcode += XLALGenerateWaveDerivative(dLNhy,LNhy_s,dt);
    errcode += XLALGenerateWaveDerivative(dLNhz,LNhz_s,dt);
    if (errcode != 0) {
      fprintf(stderr,"**** LALPSpinInspiralRD ERROR ****: error generating first derivatives: #points %d\n",Npoints);
      XLAL_ERROR(XLAL_EFAILED);
    }

    for (k=0;k<Npoints;k++) {
      LNhxy = sqrt(LNhx_s->data[k] * LNhx_s->data[k] + LNhy_s->data[k] * LNhy_s->data[k]);
      if (LNhxy > 0.) {
				diota->data[k]  = -dLNhz->data[k] / LNhxy;
				dalpha->data[k] = (LNhx_s->data[k] * dLNhy->data[k] - LNhy_s->data[k] * dLNhx->data[k]) / LNhxy;
      } else {
				diota->data[k]  = 0.;
				dalpha->data[k] = 0.;
      }
    }

    errcode  = XLALGenerateWaveDerivative(ddiota,diota,dt);
    errcode += XLALGenerateWaveDerivative(ddalpha,dalpha,dt);
    errcode += XLALGenerateWaveDerivative(ddomega,domega,dt);
    if (errcode != 0) {
      XLALPrintError("**** LALPSpinInspiralRD ERROR ****: error generating second derivatives\n");
      XLALPrintError("                     m:           : %12.5f  %12.5f\n",mparams->m1m*mparams->m,mparams->m2m*mparams->m);
      XLALPrintError("              S1:                 : %12.5f  %12.5f  %12.5f\n",S1x0,S1y0,S1z0);
      XLALPrintError("              S2:                 : %12.5f  %12.5f  %12.5f\n",S2x0,S2y0,S2z0);
      XLALPrintError("     omM %12.5f   om[%d] %12.5f\n",omegaMatch,jMatch,omega);
      XLAL_ERROR(XLAL_EFAILED);
    }
		
    if (ddomega->data[kMatch]<0.) {
      fprintf(stdout,"*** LALPSpinInspiralRD WARNING: the attach of the phenom. phase has been shifted back: m1 %12.6f  m2 %12.6f\n",mparams->m1m*mparams->m,mparams->m2m*mparams->m);
      fprintf(stdout,"  Integration returned %d\n   1025: Energy increases\n   1026: Omegadot -ve\n   1028: Omega NAN\n   1029: Omega > Omegamatch\n   1031: Omega -ve\n",intreturn); 
      while ((kMatch>0)&&(ddomega->data[kMatch]<0.)) {
				kMatch--;
				jMatch--;
      } 
    }

    phenPars->intreturn = intreturn;
    phenPars->energy    = energy[jMatch];
    phenPars->omega     = omega_s->data[kMatch];
    phenPars->domega    = domega->data[kMatch];
    phenPars->ddomega   = ddomega->data[kMatch];
    phenPars->diota     = diota->data[kMatch];
    phenPars->ddiota    = ddiota->data[kMatch];
    phenPars->dalpha    = dalpha->data[kMatch];
    phenPars->ddalpha   = ddalpha->data[kMatch];
    phenPars->countback = jMatch;
    phenPars->Psi       = Phi[jMatch];
    phenPars->endtime   = ((REAL8) jMatch)*dt;
    phenPars->ci        = LNhz[jMatch];
    phenPars->LNhS1     = LNhS1;
    phenPars->LNhS2     = LNhS2;
    phenPars->S1S2      = S1S2;
    phenPars->S1S1      = S1S1;
    phenPars->S2S2      = S2S2;
		
    XLALDestroyREAL8Vector(omega_s);
    XLALDestroyREAL8Vector(LNhx_s);
    XLALDestroyREAL8Vector(LNhy_s);
    XLALDestroyREAL8Vector(LNhz_s);
    XLALDestroyREAL8Vector(alpha_s);
    XLALDestroyREAL8Vector(dLNhx);
    XLALDestroyREAL8Vector(dLNhy);
    XLALDestroyREAL8Vector(dLNhz);
    XLALDestroyREAL8Vector(diota);
    XLALDestroyREAL8Vector(dalpha);
    XLALDestroyREAL8Vector(domega);
    XLALDestroyREAL8Vector(ddomega);
    XLALDestroyREAL8Vector(ddiota);
    XLALDestroyREAL8Vector(ddalpha);
  }
  else {
    jMatch=intlen-1;
    phenPars->intreturn = intreturn;
    phenPars->energy    = 0.;
    phenPars->omega     = 0.;
    phenPars->domega    = 0.;
    phenPars->ddomega   = 0.;
    phenPars->diota     = 0.;
    phenPars->ddiota    = 0.;
    phenPars->dalpha    = 0.;
    phenPars->ddalpha   = 0.;
    phenPars->countback = intlen-1;
    phenPars->Psi       = 0.;
    phenPars->endtime   = 0.;
    phenPars->ci        = 0.;
    phenPars->LNhS1     = 0.;
    phenPars->LNhS2     = 0.;
    phenPars->S1S2      = 0.;
    phenPars->S1S1      = 0.;
    phenPars->S2S2      = 0.;
  }
	
  /* Now fill the Hlm waveform structures*/
	
  alphaold=alpha;
  //REAL8 alphaoold = 0.;
  /*The follwoing if structure is necessary in the case L is initially parallel to 
    N so that alpha is undefined at the beginning but different from zero at the first 
    step (if the spins are not aligned with L). 
    Such a discontinuity of alpha would induce 
    a discontinuity of the waveform between its initial value and its value after the
    first integration step. This does not happen during the integration as in that 
    case alpha can be safely set to the previous value, just before L becomes parallel 
    to N. In the case L stays all the time parallel to N than alpha can be 
    safely set to zero, as it is.*/
  if ((LNhy[0]*LNhy[0]+LNhx[0]*LNhx[0])>0.) 
    alpha=atan2(LNhy[0],LNhx[0]);
  else {
    if ((S1x[0]*S1x[0]+S1y[0]*S1y[0]+S2x[0]*S2x[0]+S2y[0]*S2y[0])>0.) {
      c1=0.75+mparams->eta/2-0.75*mparams->dm;
      c2=0.75+mparams->eta/2+0.75*mparams->dm;
      alpha=atan2(-c1*S1x[0]-c2*S2x[0],c1*S1y[0]+c2*S2y[0]);
    }
    else
      alpha=0.;  
  }

  for (j=0;j<=jMatch;j++) {
		
    freq->data[j]=omega[j];
    v=cbrt(omega[j]);
    v2=v*v;
		
    // amp22= -2.0 * params->mu * LAL_MRSUN_SI/(params->distance) * sqrt( 16.*LAL_PI/5.)*v2;
    // amp20 = amp22*sqrt(3/2)
    // Y22 \pm Y2-2= sqrt(5/PI)    ((1+cos^2 t)/4, (cos t)/2)
    // Y21 \pm Y2-1= sqrt(5/PI)    ((sin t)/2, (sin 2t)/4)
    // Y20         = sqrt(15/2 PI) (sin^2 t)/4
		
    amp22 = amp22ini * v2;
    amp33 = -amp22 / 4. * sqrt(5./42.); 
    amp44 = amp22 * sqrt(5./7.) * 2./9.* v2;
		
    Psi = phase->data[j] = Phi[j];// - 2. * omega[j] * log(omega[j]);
	
    trigAngle.ci   = (LNhz[j]);
    trigAngle.cdi  = 2. * trigAngle.ci * trigAngle.ci - 1.;
    trigAngle.c2i  = trigAngle.ci * trigAngle.ci;
    trigAngle.s2i  = 1. - trigAngle.ci * trigAngle.ci;
    trigAngle.si   = sqrt(trigAngle.s2i);
    trigAngle.sdi  = 2. * trigAngle.ci * trigAngle.si;
    trigAngle.c2i2 = (1. + trigAngle.ci) / 2.;
    trigAngle.s2i2 = (1. - trigAngle.ci) / 2.;
    trigAngle.ci2  = sqrt(trigAngle.c2i2);
    trigAngle.si2  = sqrt(trigAngle.s2i2);
    trigAngle.c3i2 = trigAngle.c2i2 * trigAngle.ci2;
    trigAngle.s3i2 = trigAngle.s2i2 * trigAngle.si2;
    trigAngle.c4i2 = trigAngle.c2i2 * trigAngle.c2i2;
    trigAngle.s4i2 = trigAngle.s2i2 * trigAngle.s2i2;
    trigAngle.c5i2 = trigAngle.c4i2 * trigAngle.ci2;
    trigAngle.s5i2 = trigAngle.s4i2 * trigAngle.si2;
    trigAngle.c6i2 = trigAngle.c4i2 * trigAngle.c2i2;
    trigAngle.s6i2 = trigAngle.s4i2 * trigAngle.s2i2;
    trigAngle.c8i2 = trigAngle.c4i2 * trigAngle.c4i2;
    trigAngle.s8i2 = trigAngle.s4i2 * trigAngle.s4i2;
		
    //alphaoold = alphaold;
    alphaold  = alpha;
    if ((LNhy[j]*LNhy[j]+LNhx[j]*LNhx[j])>0.) {
      alpha = atan2(LNhy[j], LNhx[j]);
    }
    else alpha = alphaold;
		
    XLALSpinInspiralFillH2Modes(h2P2,h2M2,h2P1,h2M1,h20,j,amp22,v,mparams->eta,mparams->dm,Psi,alpha,&trigAngle);
		
    /*if (j>2) {
		 if ((alphaold*alphaoold)<0.) {
		 if ( fabs(cos(2.*(Phi[j-1]+alphaold))-cos(2.*(Phi[j-2]+alphaoold)))>0.2) {
		 fprintf(stdout,"*** LALPSpinInspiralRD WARNING ***: Possible problem with coordinate singularity:\n Step %d  LNhy: %12.6e LNhx: %12.6e  Psi+alpha: %12.6e alpha %12.6e\n Step %d  LNhy: %12.6e  LNhx: %12.6e  Psi+alpha: %12.6e  alpha %12.6e\n Step %d  LNhy: %12.6e  LNhx: %12.6e  Psi+alpha: %12.6e  alpha %12.6e\n",j,LNhy[j],LNhx[j],(Phi[j]+alpha)/LAL_PI,alpha/LAL_PI,j-1,LNhy[j-1],LNhx[j-1],(Phi[j-1]+alphaold)/LAL_PI,alphaold/LAL_PI,j-2,LNhy[j-2],LNhx[j-2],(Phi[j-2]+alphaoold)/LAL_PI,alphaoold/LAL_PI);
		 fprintf(stdout,"            m: (%12.6e,%12.6e)\n", mparams->m1m*mparams->m, mparams->m2m*mparams->m);
		 fprintf(stdout,"            S1: (%9.6f,%9.6f,%9.6f)\n",yinit[5]/mparams->m1msq,yinit[6]/mparams->m1msq,yinit[7]/mparams->m1msq);
		 fprintf(stdout,"            S2: (%9.6f,%9.6f,%9.6f)\n",yinit[8]/mparams->m2msq,yinit[9]/mparams->m2msq,yinit[10]/mparams->m2msq);
		 }
		 }
		 }*/
		
    XLALSpinInspiralFillH3Modes(h3P3,h3M3,h3P2,h3M2,h3P1,h3M1,h30,j,amp33,v,mparams->eta,mparams->dm,Psi,alpha,&trigAngle);
		
    XLALSpinInspiralFillH4Modes(h4P4,h4M4,h4P3,h4M3,h4P2,h4M2,h4P1,h4M1,h40,j,amp44,v,mparams->eta,mparams->dm,Psi,alpha,&trigAngle);
		
  }
	
  if (yin) LALFree(yin);
  if (yout) XLALDestroyREAL8Array(yout);
	
  phenPars->alpha=alpha;
	
  return XLAL_SUCCESS;
	
} /* End of the inspiral part created via the adaptive integration method */

int XLALSimIMRPSpinFinalMassSpin(
    REAL8 *finalMass,
    REAL8 *finalSpin,
    REAL8 m1,
    REAL8 m2,
    REAL8 s1x,
    REAL8 s1y,
    REAL8 s1z,
    REAL8 s2x,
    REAL8 s2y,
    REAL8 s2z,
    REAL8 energy,
    REAL8 LNhvecx,
    REAL8 LNhvecy,
    REAL8 LNhvecz
    )
{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;
  REAL8 qq,ll,eta;
	
  /* See eq.(6) in arXiv:0904.2577 */
  REAL8 ma1,ma2,a12,a12l;
  REAL8 cosa1=0.;
  REAL8 cosa2=0.;
  REAL8 cosa12=0.;
	
  REAL8 t0=-2.9;
  REAL8 t3=2.6;
  REAL8 s4=-0.123;
  REAL8 s5=0.45;
  REAL8 t2=16.*(0.6865-t3/64.-sqrt(3.)/2.);
	
  /* get a local copy of the intrinstic parameters */
  qq = m2/m1;
  eta = m1*m2/((m1+m2)*(m1+m2));
  /* done */
  ma1 = sqrt( s1x*s1x + s1y*s1y + s1z*s1z );
  ma2 = sqrt( s2x*s2x + s2y*s2y + s2z*s2z );
	
  if (ma1>0.) cosa1 = (s1x*LNhvecx + s1y*LNhvecy + s1z*LNhvecz)/ma1;
  else cosa1=0.;
  if (ma2>0.) cosa2 = (s2x*LNhvecx + s2y*LNhvecy + s2z*LNhvecz)/ma2;
  else cosa2=0.;
  if ((ma1>0.)&&(ma2>0.)) {
    cosa12  = (s1x*s2x + s1y*s2y + s1z*s2z)/ma1/ma2;
  }
  else cosa12=0.;
	
  a12  = ma1*ma1 + ma2*ma2*qq*qq*qq*qq + 2.*ma1*ma2*qq*qq*cosa12 ;
  a12l = ma1*cosa1 + ma2*cosa2*qq*qq ;
  ll = 2.*sqrt(3.)+ t2*eta + t3*eta*eta + s4*a12/(1.+qq*qq)/(1.+qq*qq) + (s5*eta+t0+2.)/(1.+qq*qq)*a12l;
	
  /* Estimate final mass by adding the negative binding energy to the rest mass*/
  *finalMass = 1. + energy;

  /* Estimate final spin */
  *finalSpin = sqrt( a12 + 2.*ll*qq*a12l + ll*ll*qq*qq)/(1.+qq)/(1.+qq);

  /* Check value of finalMass */
  if (*finalMass < 0.) {
    fprintf(stderr,"*** LALPSpinInspiralRingdownWave ERROR: Estimated final mass <0 : %12.6f\n ",*finalMass);
    fprintf(stderr,"***                                    Final mass set to initial mass\n");
    XLAL_ERROR( XLAL_ERANGE);
    *finalMass = 1.;
  }
	
  /* Check value of finalSpin */
  if ((*finalSpin > 1.)||(*finalSpin < 0.)) {
    if ((*finalSpin>=1.)&&(*finalSpin<1.01)) {
			fprintf(stderr,"*** LALPSpinInspiralRingdownWave WARNING: Estimated final Spin slightly >1 : %11.3e\n ",*finalSpin);
			fprintf(stderr,"      (m1=%8.3f  m2=%8.3f s1=(%8.3f,%8.3f,%8.3f) s2=(%8.3f,%8.3f,%8.3f) ) final spin set to 1 and code goes on\n",m1,m2,s1x,s1y,s1z,s2x,s2y,s2z);
			*finalSpin = .99999;
		}
    else {
      fprintf(stderr,"*** LALPSpinInspiralRingdownWave ERROR: Unphysical estimation of final Spin : %11.3e\n ",*finalSpin);
			fprintf(stderr,"      (m1=%8.3f  m2=%8.3f s1=(%8.3f,%8.3f,%8.3f) s2=(%8.3f,%8.3f,%8.3f) )\n",m1,m2,s1x,s1y,s1z,s2x,s2y,s2z); 
			fprintf(stderr,"***                                    Code aborts\n");
      *finalSpin = 0.;
      XLAL_ERROR( XLAL_ERANGE);
    }
  }
	
  /*For reference these are the formula used in the EOBNR construction*/
  //*finalMass = 1. - 0.057191 * eta - 0.498 * eta*eta;
  //*finalSpin = 3.464102 * eta - 2.9 * eta*eta;
	
  return errcode;
}

static int XLALPSpinGenerateQNMFreq(
    COMPLEX8Vector *modefreqs,
    UINT4 l,
    INT4 m,
    UINT4 nmodes,
    REAL8 finalMass,
    REAL8 finalSpin,
    REAL8 totalMass
    )
{
	
	
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;
  UINT4 i;
  /* Fitting coefficients for QNM frequencies from PRD73, 064030, gr-qc/0512160, tables VIII and IX */
  REAL4 BCW22re[3][3]  = { {1.5251, -1.1568,  0.1292}, {1.3673, -1.0260,  0.1628}, { 1.3223, -1.0257,  0.1860} };
  REAL4 BCW22im[3][3]  = { {0.7000,  1.4187, -0.4990}, {0.1000,  0.5436, -0.4731}, {-0.1000,  0.4206, -0.4256} };
	
  /*REAL4 BCW2m2re[3][3] = { {0.2938,  0.0782,  1.3546}, {0.2528,  0.0921,  1.3344}, { 0.1873,  0.1117,  1.3322} };
	 REAL4 BCW2m2im[3][3] = { {1.6700,  0.4192,  1.4700}, {0.4550,  0.1729,  1.3617}, { 0.1850,  0.1266,  1.3661} };*/
	
  REAL4 BCW21re[3][3]  = { {0.60000, -0.2339, 0.4175}, {0.5800, -0.2416, 0.4708}, { 0.5660, -0.2740, 0.4960} };
  REAL4 BCW21im[3][3]  = { {-0.30000, 2.3561, -0.2277}, {-0.3300, 0.9501, -0.2072}, { -0.1000, 0.4173, -0.2774} };
	
  /*REAL4 BCW2m1re[3][3] = { {0.3441, 0.0293, 2.0010}, {0.3165, 0.0301, 2.3415}, {0.2696, 0.0315, 2.7755} };
	 REAL4 BCW2m1im[3][3] = { {2.0000, 0.1078, 5.0069}, {0.6100, 0.0276, 13.1683}, {0.2900, 0.0276, 6.4715} };*/
	
  REAL4 BCW20re[3][3]  = { {0.4437, -0.0739,  0.3350}, {0.4185, -0.0768,  0.4355}, { 0.3734, -0.0794,  0.6306} };
  REAL4 BCW20im[3][3]  = { {4.0000,  -1.9550, 0.1420}, {1.2500,  -0.6359, 0.1614}, {0.5600,  -0.2589, -0.3034} };
	
  REAL4 BCW33re[3][3]  = { {1.8596, -1.3043, 0.1818}, {1.8566, -1.2818, 0.1934}, {1.8004, -1.2558, 0.2133} };
  REAL4 BCW33im[3][3]  = { {0.9000, 2.3430, -0.4810}, {0.2274, 0.8173, -0.4731}, {0.0400, 0.5445, -0.4539} };
	
  /*REAL4 BCW3m3re[3][3] = { {0.4673, 0.1296, 1.3255}, {0.4413, 0.1387, 1.3178}, {0.3933, 0.1555, 1.3037} };
	 REAL4 BCW3m3im[3][3] = { {2.5500, 0.6576, 1.3378}, {0.7900, 0.2381, 1.3706}, {0.4070, 0.1637, 1.3819} };*/
	
  REAL4 BCW32re[3][3]  = { {1.1481, -0.5552, 0.3002}, {1.1226, -0.5471, 0.3264}, {1.0989, -0.5550, 0.3569} };
  REAL4 BCW32im[3][3]  = { {0.8313, 2.3773, -0.3655}, {0.2300, 0.8025, -0.3684}, {0.1000, 0.4804, -0.3784}};
	
  /*REAL4 BCW3m2re[3][3] = { {0.5158, 0.8195, 1.408}, {0.4413, 0.1378, 1.3178}, {0.4567, 0.09300, 1.4469} };
	 REAL4 BCW3m2im[3][3] = { {2.9000, 0.3365, 2.3050}, {0.9000, 0.1295, 1.6142}, {0.4900, 0.0848, 1.9737} };*/
	
  REAL4 BCW31re[3][3]  = { {0.8345, -0.2405, 0.4095}, {0.8105, -0.2342, 0.4660}, {0.7684, -0.2252, 0.5805} };
  REAL4 BCW31im[3][3]  = { {23.8450, -20.724, 0.03837}, {8.8530, -7.8506, 0.03418}, {2.1800, -1.6273, 0.1163} };
	
  /*REAL4 BCW3m1re[3][3] = { {0.5751, 0.02508, 3.1360}, {0.5584, 0.02514, 3.4154}, {0.5271, 0.02561, 3.8011} };
	 REAL4 BCW3m1im[3][3] = { {3.0464, 0.1162, -0.2812}, {1.2000, -0.1928, 0.1037}, {1.0000, -0.4424, 0.02467} };*/
	
  REAL4 BCW30re[3][3]  = { {0.6873, -0.09282, 0.3479}, {0.6687, -0.09155, 0.4021}, {0.6343, -0.08915, 0.5117} };
  REAL4 BCW30im[3][3]  = { {6.7841, -3.6112, 0.09480}, {2.0075, -0.9930, 0.1197}, {0.9000, -0.3409, 0.2679} };
	
  REAL4 BCW44re[3][3]  = { {2.3, -1.5056, 0.2244}, {2.3, -1.5173, 0.2271}, {2.3, -1.5397, 0.2321} };
  REAL4 BCW44im[3][3]  = { {1.1929, 3.1191, -0.4825}, {0.3, 1.1034, -0.4703}, {0.11, 0.6997, -0.4607} };
	
  /*REAL4 BCW4m4re[3][3]  = { {0.6256, 0.18, 1.3218}, {0.6061, 0.1869, 1.3168}, {0.5686, 0.2003, 1.3068} };
	 REAL4 BCW4m4im[3][3]  = { {3.4, 0.8696, 1.4074}, {1.08, 0.3095, 1.3279}, {0.5980, 0.2015, 1.3765} };*/
	
  REAL4 BCW43re[3][3] = { {1.6869, -0.8862, 0.2822}, {1.6722, -0.8843, 0.2923}, {1.6526, -0.8888, 0.3081} };
  REAL4 BCW43im[3][3] = { {1.4812, 2.8096, -0.4271}, {0.4451, 0.9569, -0.425}, {0.22, 0.5904, -0.4236} };
	
  /*REAL4 BCW4m3re[3][3] = { {0.6728, 0.1338, 1.3413}, {0.6562, 0.1377, 1.3456}, {0.6244, 0.1454, 1.3513} };
	 REAL4 BCW4m3im[3][3] = { {3.7, 0.5829, 1.6681}, {1.18, 0.2111, 1.4129}, {0.66, 0.1385, 1.3742} };*/
	
  REAL4 BCW42re[3][3]  = { {1.2702, -0.4685, 0.3835}, {1.2462, -0.4580, 0.4139}, {1.2025, -0.4401, 0.4769} };
  REAL4 BCW42im[3][3]  = { {-3.6, 7.7749, -0.1491}, {-1.5, 2.8601, -0.1392}, {-1.5, 2.2784, -0.1124}};
	
  /*REAL4 BCW4m2re[3][3] = { {0.7294, 0.07842, 1.5646}, {0.7154, 0.07979, 1.5852}, {0.6885, 0.08259, 1.6136} };
	 REAL4 BCW4m2im[3][3] = { {4., 0.2777, 2.0647}, {1.32, 0.08694, 4.3255}, {0.75, 0.05803, 3.7971} };*/
	
  REAL4 BCW41re[3][3]  = { {1.0507, -0.2478, 0.4348}, {1.0337, -0.2439, 0.4695}, {1.0019, -0.2374, 0.5397} };
  REAL4 BCW41im[3][3]  = { {14., -9.8240, 0.09047}, {4.2, -2.8399, 0.1081}, {2.2, -1.4195, 0.1372} };
	
  /*REAL4 BCW4m1re[3][3] = { {0.7908, 0.02024, 5.4628}, {0.7785, 0.02005, 5.8547}, {0.7549, 0.01985, 6.5272} };
	 REAL4 BCW4m1im[3][3] = { {4.6, -0.4038, 0.4629}, {1.6, -0.2323, 0.2306}, {1.6, -0.8136, 0.03163} };*/
	
  REAL4 BCW40re[3][3]  = { {0.9175, -0.1144, 0.3511}, {0.9028, -0.1127, 0.3843}, {0.8751, -0.1096, 0.4516} };
  REAL4 BCW40im[3][3]  = { {7.0, -2.7934, 0.1708}, {2.2, -0.8308, 0.2023}, {1.2, -0.4159, 0.2687} };
	
  /* QNM frequencies from the fitting given in PRD73, 064030 */
	
  if ((l==2)&&(abs(m)==2)) {
    for (i = 0; i < nmodes; i++)
		{
			modefreqs->data[i] = BCW22re[i][0] + BCW22re[i][1] * pow(1.- finalSpin, BCW22re[i][2]);
			modefreqs->data[i] += I * creal(modefreqs->data[i]) / 2.
			/ (BCW22im[i][0] + BCW22im[i][1] * pow(1.- finalSpin, BCW22im[i][2]));
			modefreqs->data[i] *= 1./ finalMass / (totalMass * LAL_MTSUN_SI);
		}
  }
  else {
    if ((l==2)&&(m==0)) {
      for (i = 0; i < nmodes; i++)
			{
				modefreqs->data[i] = BCW20re[i][0] + BCW20re[i][1] * pow(1.- finalSpin, BCW20re[i][2]);
				modefreqs->data[i] += I * creal(modefreqs->data[i]) / 2.
				/ (BCW20im[i][0] + BCW20im[i][1] * pow(1.- finalSpin, BCW20im[i][2]));
				modefreqs->data[i] /= finalMass * totalMass * LAL_MTSUN_SI;
			}
    }
    else {
      if ((l==2)&&(abs(m)==1)) {
				for (i = 0; i < nmodes; i++) {
					modefreqs->data[i] = BCW21re[i][0] + BCW21re[i][1] * pow(1.- finalSpin, BCW21re[i][2]);
					modefreqs->data[i] += I * creal(modefreqs->data[i]) / 2.
					/ (BCW21im[i][0] + BCW21im[i][1] * pow(1.- finalSpin, BCW21im[i][2]));
					modefreqs->data[i] /= finalMass * totalMass * LAL_MTSUN_SI;
				}
      }
      else {
				if ((l==3)&&(abs(m)==3)) {
					for (i = 0; i < nmodes; i++) {
						modefreqs->data[i] = BCW33re[i][0] + BCW33re[i][1] * pow(1.- finalSpin, BCW33re[i][2]);
						modefreqs->data[i] += I * creal(modefreqs->data[i]) / 2.
						/ (BCW33im[i][0] + BCW33im[i][1] * pow(1.- finalSpin, BCW33im[i][2]));
						modefreqs->data[i] /= finalMass * totalMass * LAL_MTSUN_SI;
					}
				}
				else
					if ((l==3)&&(abs(m)==2)) {
						for (i = 0; i < nmodes; i++) {
							modefreqs->data[i] = BCW32re[i][0] + BCW32re[i][1] * pow(1.- finalSpin, BCW32re[i][2]);
							modefreqs->data[i] += I * creal(modefreqs->data[i]) / 2.
							/ (BCW32im[i][0] + BCW32im[i][1] * pow(1.- finalSpin, BCW32im[i][2]));
							modefreqs->data[i] /= finalMass * totalMass * LAL_MTSUN_SI;
						}
					}
					else {
						if ((l==3)&&(abs(m)==1)) {
							for (i = 0; i < nmodes; i++) {
								modefreqs->data[i] = BCW31re[i][0] + BCW31re[i][1] * pow(1.- finalSpin, BCW31re[i][2]);
								modefreqs->data[i] += I * creal(modefreqs->data[i]) / 2.
								/ (BCW31im[i][0] + BCW31im[i][1] * pow(1.- finalSpin, BCW31im[i][2]));
								modefreqs->data[i] /= finalMass * totalMass * LAL_MTSUN_SI;
							}
						}
						else {
							if ((l==3)&&(m==0)) {
								for (i = 0; i < nmodes; i++) {
									modefreqs->data[i] = BCW30re[i][0] + BCW30re[i][1] * pow(1.- finalSpin, BCW30re[i][2]);
									modefreqs->data[i] += I * creal(modefreqs->data[i]) / 2.
									/ (BCW30im[i][0] + BCW30im[i][1] * pow(1.- finalSpin, BCW30im[i][2]));
									modefreqs->data[i] /= finalMass * totalMass * LAL_MTSUN_SI;
								}
							}
							else {
								if ((l==4)&&(abs(m)==4)) {
									for (i = 0; i < nmodes; i++) {
										modefreqs->data[i] = BCW44re[i][0] + BCW44re[i][1] * pow(1.- finalSpin, BCW44re[i][2]);
										modefreqs->data[i] += I * creal(modefreqs->data[i]) / 2.
										/ (BCW44im[i][0] + BCW44im[i][1] * pow(1.- finalSpin, BCW44im[i][2]));
										modefreqs->data[i] /= finalMass * totalMass * LAL_MTSUN_SI;
									}
								}
								else {
									if ((l==4)&&(abs(m)==3)) {
										for (i = 0; i < nmodes; i++) {
											modefreqs->data[i] = BCW43re[i][0] + BCW43re[i][1] * pow(1.- finalSpin, BCW43re[i][2]);
											modefreqs->data[i] += I * creal(modefreqs->data[i]) / 2.
											/ (BCW43im[i][0] + BCW43im[i][1] * pow(1.- finalSpin, BCW43im[i][2]));
											modefreqs->data[i] /= finalMass * totalMass * LAL_MTSUN_SI;
										}
									}
									else {
										if ((l==4)&&(abs(m)==2)) {
											for (i = 0; i < nmodes; i++) {
												modefreqs->data[i] = BCW42re[i][0] + BCW42re[i][1] * pow(1.- finalSpin, BCW42re[i][2]);
												modefreqs->data[i] += I * creal(modefreqs->data[i]) / 2.
												/ (BCW42im[i][0] + BCW42im[i][1] * pow(1.- finalSpin, BCW42im[i][2]));
												modefreqs->data[i] /= finalMass * totalMass * LAL_MTSUN_SI;
											}
										}
										else {
											if ((l==4)&&(abs(m)==1)) {
												for (i = 0; i < nmodes; i++) {
													modefreqs->data[i] = BCW41re[i][0] + BCW41re[i][1] * pow(1.- finalSpin, BCW41re[i][2]);
													modefreqs->data[i] += I * creal(modefreqs->data[i]) / 2.
													/ (BCW41im[i][0] + BCW41im[i][1] * pow(1.- finalSpin, BCW41im[i][2]));
													modefreqs->data[i] /= finalMass * totalMass * LAL_MTSUN_SI;
												}
											}
											else {
												if ((l==4)&&(m==0)) {
													for (i = 0; i < nmodes; i++) {
														modefreqs->data[i] = BCW40re[i][0] + BCW40re[i][1] * pow(1.- finalSpin, BCW40re[i][2]);
														modefreqs->data[i] += I * creal(modefreqs->data[i]) / 2.
														/ (BCW40im[i][0] + BCW40im[i][1] * pow(1.- finalSpin, BCW40im[i][2]));
														modefreqs->data[i] /= finalMass * totalMass * LAL_MTSUN_SI;
													}
												}
												else {
													fprintf(stderr,"*** LALPSpinInspiralRingdownWave ERROR: Ringdown modes for l=%d m=%d not available\n",l,m);
													XLAL_ERROR( XLAL_EDOM );
												}
											}
										}
									}
								}
							}
						}
					}
      }
    }
  }
	
  return errcode;
}

static int XLALPSpinInspiralRingdownWave (
    REAL8Vector *rdwave,
    REAL8Vector *matchinspwave,
    COMPLEX8Vector *modefreqs,
    UINT4 nmodes,
    REAL8 tsampling
    )
{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;
	
  /* Needed to check GSL return codes */
  INT4 gslStatus;
	
  UINT4 i, j;
	
  /* Sampling rate from input */
  REAL8 dt;
  gsl_matrix *coef;
  gsl_vector *hderivs;
  gsl_vector *x;
  gsl_permutation *p;
  REAL8Vector *modeamps;
	
  int s;
  REAL8 tj;
	
  dt = 1.0 / tsampling;
	
  if ( modefreqs->length != nmodes )
	{
		XLAL_ERROR( XLAL_EBADLEN );
	}
	
  /* Solving the linear system for QNMs amplitude coefficients using gsl routine */
  /* Initialize matrices and supporting variables */
  XLAL_CALLGSL( coef = (gsl_matrix *) gsl_matrix_alloc(2 * nmodes, 2 * nmodes) );
  XLAL_CALLGSL( hderivs = (gsl_vector *) gsl_vector_alloc(2 * nmodes) );
  XLAL_CALLGSL( x = (gsl_vector *) gsl_vector_alloc(2 * nmodes) );
  XLAL_CALLGSL( p = (gsl_permutation *) gsl_permutation_alloc(2 * nmodes) );
	
  /* Check all matrices and variables were allocated */
  if ( !coef || !hderivs || !x || !p )
	{
		if (coef)    gsl_matrix_free(coef);
		if (hderivs) gsl_vector_free(hderivs);
		if (x)       gsl_vector_free(x);
		if (p)       gsl_permutation_free(p);
		XLAL_ERROR( XLAL_ENOMEM );
	}
	
  /* Define the linear system Ax=y */
  /* Matrix A (2*nmodes by 2*nmodes) has block symmetry. Define half of A here as "coef" */
  /* Define y here as "hderivs" */
	
  j=0;
  while (j<nmodes) {
    if (j==0) {
      for (i = 0; i < nmodes; i++) {
				gsl_matrix_set(coef, 2*j, i, 1.);
				gsl_matrix_set(coef, 2*j, i+nmodes, 0.);
				gsl_matrix_set(coef, 2*j+1, i, -cimag(modefreqs->data[i]));
				gsl_matrix_set(coef, 2*j+1, i+nmodes, creal(modefreqs->data[i]));
      }
    }
    else {
      if (j==1) {
				for (i = 0; i < nmodes; i++) {
					gsl_matrix_set(coef, 2*j, i, cimag(modefreqs->data[i])*cimag(modefreqs->data[i])-creal(modefreqs->data[i])*creal(modefreqs->data[i]));
					gsl_matrix_set(coef, 2*j, i+nmodes, -2.*cimag(modefreqs->data[i])*creal(modefreqs->data[i]));
					gsl_matrix_set(coef, 2*j+1, i, -cimag(modefreqs->data[i])*cimag(modefreqs->data[i])*cimag(modefreqs->data[i])+3.*cimag(modefreqs->data[i])*creal(modefreqs->data[i])*creal(modefreqs->data[i]));
					gsl_matrix_set(coef, 2*j+1, i+nmodes, -creal(modefreqs->data[i])*creal(modefreqs->data[i])*creal(modefreqs->data[i])+3.*creal(modefreqs->data[i])*cimag(modefreqs->data[i])*cimag(modefreqs->data[i]));
				}
      }
      else {
				if (j==2) {
					for (i = 0; i < nmodes; i++) {
						gsl_matrix_set(coef, 2*j, i, pow(cimag(modefreqs->data[i]),4.)+pow(creal(modefreqs->data[i]),4.)-6.*pow(creal(modefreqs->data[i])*cimag(modefreqs->data[i]),2.));
						gsl_matrix_set(coef, 2*j, i+nmodes, -4.*pow(cimag(modefreqs->data[i]),3.)*creal(modefreqs->data[i])+4.*pow(creal(modefreqs->data[i]),3.)*cimag(modefreqs->data[i]));
						gsl_matrix_set(coef, 2*j+1, i, -pow(cimag(modefreqs->data[i]),5.)+10.*pow(cimag(modefreqs->data[i]),3.)*pow(creal(modefreqs->data[i]),2.)-5.*cimag(modefreqs->data[i])*pow(creal(modefreqs->data[i]),4.));
						gsl_matrix_set(coef, 2*j+1, i+nmodes, 5.*pow(cimag(modefreqs->data[i]),4.)*creal(modefreqs->data[i])-10.*pow(cimag(modefreqs->data[i]),2.)*pow(creal(modefreqs->data[i]),3.)+pow(creal(modefreqs->data[i]),5.));
					}
				}
				else {
					fprintf(stderr,"*** LALPSpinInspiralRingDown ERROR ***: nmode must be <=2, %d selected\n",nmodes);
					gsl_matrix_free(coef);
					gsl_vector_free(hderivs);
					gsl_vector_free(x);
					gsl_permutation_free(p);
					XLAL_ERROR( XLAL_EDOM );
				}
      }
    }
    gsl_vector_set(hderivs, 2*j, matchinspwave->data[2*j]);
    gsl_vector_set(hderivs, 2*j+1, matchinspwave->data[2*j+1]);
    j++;
  }
	
  /* Call gsl LU decomposition to solve the linear system */
  XLAL_CALLGSL( gslStatus = gsl_linalg_LU_decomp(coef, p, &s) );
  if ( gslStatus == GSL_SUCCESS )
  {
    XLAL_CALLGSL( gslStatus = gsl_linalg_LU_solve(coef, p, hderivs, x) );
  }
	
  if ( gslStatus != GSL_SUCCESS )
  {
    gsl_matrix_free(coef);
    gsl_vector_free(hderivs);
    gsl_vector_free(x);
    gsl_permutation_free(p);
    XLAL_ERROR( XLAL_EFUNC );
  }
	
  /* Putting solution to an XLAL vector */
  modeamps = XLALCreateREAL8Vector(2*nmodes);
	
  if ( !modeamps )
  {
    gsl_matrix_free(coef);
    gsl_vector_free(hderivs);
    gsl_vector_free(x);
    gsl_permutation_free(p);
    XLAL_ERROR( XLAL_ENOMEM );
  }
	
  for (i = 0; i < 2*nmodes; i++) {
    modeamps->data[i] = gsl_vector_get(x, i);
  }
	
  /* Free all gsl linear algebra objects */
  gsl_matrix_free(coef);
  gsl_vector_free(hderivs);
  gsl_vector_free(x);
  gsl_permutation_free(p);
	
  /* Build ring-down waveforms */
  UINT4 Nrdwave=rdwave->length;
  for (j = 0; j < Nrdwave; j++) {
    tj = j * dt;
    rdwave->data[j] = 0.;
    for (i = 0; i < nmodes; i++) {
      rdwave->data[j] += exp(- tj * cimag(modefreqs->data[i]))
			* ( modeamps->data[i] * cos(tj * creal(modefreqs->data[i]))
				 +   modeamps->data[i + nmodes] * sin(tj * creal(modefreqs->data[i])) );
    }
  }
	
  XLALDestroyREAL8Vector(modeamps);
  return errcode;
} /*End of XLALPSpinInspiralRingdownWave */

static int XLALPSpinInspiralAttachRingdownWave (
    REAL8Vector *sigl,
    UINT4 *attpos,
    UINT4 nmodes,
    UINT4 l,
    INT4 m,
    REAL8 finalMass,
    REAL8 finalSpin,
    REAL8 tSampling,
    REAL8 totalMass
    )
{
  const UINT4 Npatch=40;
  const UINT4 offsetAttch = 2;
	
  COMPLEX8Vector *modefreqs;
  UINT4 Nrdwave;
	
  UINT4 i=0;
  UINT4 j=0;
  UINT4 k=0;
  UINT4 atpos;
  INT4 errcode;
	
  REAL8Vector	*rdwave;
  REAL8Vector	*inspwave,*dinspwave;
  REAL8Vector	*matchinspwave;
  REAL8 dt;
	
  dt = 1./tSampling;
  atpos=(*attpos);
	
  /* Create memory for the QNM frequencies */
  modefreqs = XLALCreateCOMPLEX8Vector( nmodes );
  if ( !modefreqs )
  {
    XLAL_ERROR( XLAL_ENOMEM );
  }
  errcode = XLALPSpinGenerateQNMFreq( modefreqs, l, m, nmodes, finalMass, finalSpin, totalMass);
  if ( errcode != XLAL_SUCCESS )
  {
    XLALDestroyCOMPLEX8Vector( modefreqs );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Ringdown signal length: 10 times the decay time of the n=0 mode */
  Nrdwave = (INT4) (10. / cimag(modefreqs->data[0]) / dt);
  /* Patch length, centered around the matching point "attpos" */

  (*attpos)+=Nrdwave;
	
  /* Check the value of attpos, to prevent memory access problems later */
  if ( atpos < Npatch || atpos + Npatch >= sigl->length )
  {
  XLALPrintError( "Value of attpos inconsistent with given value of Npatch: atpos=%d  Npatch=%d, sign->length=%d, totalMass=11.5f\n",atpos,Npatch,sigl->length,totalMass);
  XLALDestroyCOMPLEX8Vector( modefreqs );
  XLAL_ERROR( XLAL_EFAILED );
  }
	
  /* Create memory for the ring-down and full waveforms, derivatives of inspirals 
  and waveforms and its derivative values at the attach point */
	
  rdwave = XLALCreateREAL8Vector( Nrdwave );
  inspwave = XLALCreateREAL8Vector( Npatch );
  dinspwave = XLALCreateREAL8Vector( Npatch );
  matchinspwave = XLALCreateREAL8Vector( 2*nmodes );
	
  /* Check memory was allocated */
  if ( !rdwave || !inspwave || !dinspwave || !matchinspwave )
  {
    XLALDestroyCOMPLEX8Vector( modefreqs );
    if (rdwave) XLALDestroyREAL8Vector( rdwave );
    if (inspwave) XLALDestroyREAL8Vector( inspwave );
    if (dinspwave) XLALDestroyREAL8Vector( dinspwave );
    if (matchinspwave) XLALDestroyREAL8Vector( matchinspwave );
    XLAL_ERROR( XLAL_ENOMEM );
  }
	
  /* Generate derivatives of the last part of inspiral waves */
  /* Take the last part of sigl1 */
	
  for (i=0; i<2; i++) {
    /* i=0(1) for real(imaginary) part */
    for (j = 0; j < Npatch; j++) {
      inspwave->data[j]    = sigl->data[2*(atpos - Npatch + j)+i];	
    }
		
    for (k=0;k<2*nmodes;k++) {
      matchinspwave->data[k] = inspwave->data[Npatch-1-offsetAttch];
      if ((k+1)<2*nmodes) {
        errcode = XLALGenerateWaveDerivative( dinspwave, inspwave, dt);
        if ( (errcode != XLAL_SUCCESS) ) {
          XLALDestroyCOMPLEX8Vector( modefreqs );
          XLALDestroyREAL8Vector( rdwave );
          XLALDestroyREAL8Vector( inspwave );
          XLALDestroyREAL8Vector( dinspwave );
          XLALDestroyREAL8Vector( matchinspwave );
          XLAL_ERROR( XLAL_EFUNC );
        }
        for (j=0; j<Npatch; j++) {
          inspwave->data[j]=dinspwave->data[j];
        }
      }
    }
		
    errcode = XLALPSpinInspiralRingdownWave (rdwave,matchinspwave,modefreqs,nmodes,tSampling);

    if ( errcode != XLAL_SUCCESS ) {
      XLALDestroyCOMPLEX8Vector( modefreqs );
      XLALDestroyREAL8Vector( rdwave );
      XLALDestroyREAL8Vector( inspwave );
      XLALDestroyREAL8Vector( dinspwave );
      XLALDestroyREAL8Vector( matchinspwave );
      XLAL_ERROR( XLAL_EFUNC );
    }
    /* Generate full waveforms, by stitching inspiral and ring-down waveforms */
		
    for (j = 0; j < Nrdwave; j++) {
      sigl->data[2*j + 2*(atpos - 1 - offsetAttch) + i ] = rdwave->data[j];
    }	
  }
	
  /* Free memory */
  XLALDestroyCOMPLEX8Vector( modefreqs );
  XLALDestroyREAL8Vector( rdwave );
  XLALDestroyREAL8Vector( inspwave );
  XLALDestroyREAL8Vector( dinspwave );
  XLALDestroyREAL8Vector( matchinspwave );
	
  return errcode;
}


/**
 * Driver routine to compute the PhenSpin Inspiral Ringdown waveform.
 *
 * All units are SI units.
 */
int XLALSimIMRPSpinInspiralRDGenerator(
    REAL8TimeSeries **hplus,	/**< +-polarization waveform */
    REAL8TimeSeries **hcross,	/**< x-polarization waveform */
    REAL8 phi_start,            /**< start phase */
    REAL8 deltaT,               /**< sampling interval */
    REAL8 m1,                   /**< mass of companion 1 */
    REAL8 m2,                   /**< mass of companion 2 */
    REAL8 f_min,                /**< start frequency */
    REAL8 r,                    /**< distance of source */
    REAL8 iota,                 /**< inclination of source (rad) */
    REAL8 s1x,                  /**< x-component of dimensionless spin for object 1 */
    REAL8 s1y,                  /**< y-component of dimensionless spin for object 1 */
    REAL8 s1z,                  /**< z-component of dimensionless spin for object 1 */
    REAL8 s2x,                  /**< x-component of dimensionless spin for object 2 */
    REAL8 s2y,                  /**< y-component of dimensionless spin for object 2 */
    REAL8 s2z,                  /**< z-component of dimensionless spin for object 2 */
    int phaseO,                 /**< twice post-Newtonian phase order */
    LALSimInspiralFrameAxis axisChoice, /**< Choice of axis for input spin params */
    int inspiralOnly            /**< 0 generate RD, 1 generate inspiralOnly*/
    )
{

  LIGOTimeGPS t_start = LIGOTIMEGPSZERO;
  UINT4 length;                // signal vector length
  REAL8 tn = XLALSimInspiralTaylorLength(deltaT, m1, m2, f_min, phaseO);
  REAL8 x = 1.1 * (tn + 1. ) / deltaT;
  length = ceil(log10(x)/log10(2.));
  length = pow(2, length);

  if (*hplus == NULL)
    *hplus = XLALCreateREAL8TimeSeries("H+", &t_start, 0.0, deltaT, &lalDimensionlessUnit, length);

  if (*hcross == NULL)
    *hcross = XLALCreateREAL8TimeSeries("Hx", &t_start, 0.0, deltaT, &lalDimensionlessUnit, length);

  if(*hplus == NULL || *hcross == NULL)
    XLAL_ERROR(XLAL_ENOMEM);

  /* declare code parameters and variables */
  const INT4 neqs = 11+1;      // number of dynamical variables plus the energy function
  UINT4 count,apcount;         // integration steps performed
  UINT4 i, j, k, l;            // counters

  REAL8 v = 0.;
  REAL8 v2 = 0.;
  REAL8 v2old;
  REAL8 tim;                   // time (units of total mass)
  REAL8 unitHz;
  REAL8 initomega,initphi;
  REAL8 inc;
  REAL8 LNhmag,initJmag;
  REAL8 initS1[3],initS2[3],initLNh[3],initJ[3];
  REAL8 iS1[3],iS2[3];
  REAL8 phiJ,thetaJ;
  REAL8 ry[3][3],rz[3][3];
  REAL8 dt;
  REAL8 tSampling = 1./deltaT;
  INT4  intreturn;
  REAL8 yinit[neqs];

  REAL8Vector* h2P2;
  REAL8Vector* h2M2;
  REAL8Vector* h2P1;
  REAL8Vector* h2M1;
  REAL8Vector* h20;
  REAL8Vector* h3P3;
  REAL8Vector* h3M3;
  REAL8Vector* h3P2;
  REAL8Vector* h3M2;
  REAL8Vector* h3P1;
  REAL8Vector* h3M1;
  REAL8Vector* h30;
  REAL8Vector* h4P4;
  REAL8Vector* h4M4;
  REAL8Vector* h4P3;
  REAL8Vector* h4M3;
  REAL8Vector* h4P2;
  REAL8Vector* h4M2;
  REAL8Vector* h4P1;
  REAL8Vector* h4M1;
  REAL8Vector* h40;

  REAL8Vector* sigp;
  REAL8Vector* sigc;
  REAL8Vector* fap;
  REAL8Vector* phap;

  LALPSpinInspiralPhenPars phenPars;
  LALSpinInspiralAngle trigAngle;
  LALPSpinInspiralRDparams mparams;

  REAL8 Psi=0.;
  REAL8 amp22ini,amp22,amp33,amp44;

  REAL8 alpha=0.;
  REAL8 t0,tAs;
  REAL8 om0,om1,om;
  REAL8 Psi0,alpha0;
  REAL8 dalpha0,dalpha1;
  //REAL8 omold,iota0,diota0,diota1;
  REAL8 LNhS1,LNhS2;
  REAL8 S1S1,S1S2,S2S2;

  COMPLEX8Vector *modefreqs;
  COMPLEX16 MultSphHarmP;       // Generic spin-weighted spherical harmonics
  COMPLEX16 MultSphHarmM;       // Generic spin-weighted spherical harmonics
  REAL8 x0, x1, x2, x3;

  /* The number of Ring Down modes is hard-coded here */
  const UINT4 nmodes=2;
  /* Nmodes should be restricted to either 1 or 2*/

  UINT4 errcode;

  REAL8 finalMass,finalSpin;
  REAL8 energy=0.;
  REAL8 omegaMatch;
  REAL8 frOmRD,omegaRD;
  REAL8 mass1 = m1 / LAL_MSUN_SI; /* mass of object 1 in solar masses */
  REAL8 mass2 = m2 / LAL_MSUN_SI; /* mass of object 2 in solar masses */
  REAL8 totalMass = mass1 + mass2;
  REAL8 eta = mass1 * mass2 / (totalMass * totalMass);
  REAL8 mass = (m1 + m2) * LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m from kilograms to seconds */

  REAL8 mu = eta * totalMass; /* Reduced mass in solar masses */
  unitHz = mass * (REAL8) LAL_PI;

  /*
	 if ((signalvec2)||(hh))
    params->nStartPad = 0;*/    /* must be zero for templates and injection */
  /* -- length in seconds from Newtonian formula; */

  dt = deltaT;

  /* setup coefficients for PN equations */
  if(XLALPSpinInspiralRDparamsSetup(&mparams, /** Output: RDparams structure */
			     inspiralOnly, 	/** Do not Only generate inspiral */
			     deltaT, 		/** sampling interval */
			     f_min,		/** Starting frequency */
			     mass1,		/** Mass 1 */
			     mass2,		/** Mass 2 */
			     LAL_SIM_INSPIRAL_SPIN_ORDER_ALL,	/** Spin order */
			     phaseO		/** PN Order in Phase */ ))
    XLAL_ERROR(XLAL_EFUNC);

  /* Check that initial frequency is smaller than omegamatch ~ xxyy for m=100 Msun */
  initphi   = phi_start;
  initomega = f_min*unitHz;

  /* Check that initial frequency is smaller than omegamatch ~ xxyy for m=100 Msun */

  LNhS1 = s1z;
  LNhS2 = s2z;
  S1S1 = s1x*s1x + s1y*s1y + s1z*s1z;
  S1S2 = s1x*s2x + s1y*s2y + s1z*s2z;
  S2S2 = s2x*s2x + s2y*s2y + s2z*s2z;


  omegaMatch = OmMatch(LNhS1,LNhS2,S1S1,S1S2,S2S2);

  if ( initomega > omegaMatch ) {
    /*if ((params->spin1[0]==params->spin1[1])&&(params->spin1[1]==params->spin2[0])&&(params->spin2[0]==params->spin2[1])&&(params->spin2[1]==0.)) {
      //Beware, this correspond to a shift of the initial phase!
      initomega = 0.95*omegaMatch;
      fprintf(stdout,"*** LALPSpinInspiralRD WARNING ***: Initial frequency reset from %12.6e to %12.6e Hz, m:(%12.4e,%12.4e)\n",params->fLower,initomega/unitHz,params->mass1,params->mass2);
      }*/
    /*else {*/
    XLALPrintError("**** LALPSpinInspiralRD ERROR ****: the product of initial frequency times masses is too high: %11.5e for omM ~ %11.5e\n",f_min*mass*LAL_PI,omegaMatch);
    XLALPrintError("                                    please consider decreasing initial freq %8.3f Hz or m:(%8.3f, %8.3f) Msun\n",f_min,mass1,mass2);
    XLAL_ERROR(XLAL_EFAILED);
    /*}*/
  }

  /* Here we use the following convention:
     the coordinates of the spin vectors spin1,2 and the inclination 
     variable refers to different physical parameters according to the value of 
     axisChoice:

     * LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L: inclination denotes the angle 
                 between the view direction N and the initial L 
                 (initial L//z, N in the x-z plane) and the spin 
		 coordinates are given with respect to initial L.
     * LAL_SIM_INSPIRAL_FRAME_AXIS_TOTAL_J:   inclination denotes the angle 
                 between the view direction and J (J is constant during the 
                 evolution, J//z, both N and initial L are in the x-z plane) 
                 and the spin coordinates are given wrt initial ** L **.
     * LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW:     inclination denotes the angle 
                 between the initial L and N (N//z, initial L in the x-z plane)
                 and the spin coordinates are given with respect to N.

     In order to reproduce the results of the SpinTaylor code View must be chosen.
     The spin magnitude are normalized to the individual mass^2, i.e.
     they are dimension-less.
     The modulus of the initial angular momentum is fixed by m1,m2 and
     initial frequency.
     The polarization angle is not used here, it enters the pattern
     functions along with the angles marking the sky position of the
     source. However a misalignment of the line of nodes wrt to the standard
     convention can take place when L is initially aligned with N, this is 
     re-absorbed by a redefinition of the psi angle.*/

  // Physical magnitude of the orbital angular momentum
  LNhmag = eta * totalMass * totalMass / cbrt(initomega);

  // Physical values of the spins
  initS1[0] = s1x * mass1 * mass1;
  initS1[1] = s1y * mass1 * mass1;
  initS1[2] = s1z * mass1 * mass1;
  initS2[0] = s2x * mass2 * mass2;
  initS2[1] = s2y * mass2 * mass2;
  initS2[2] = s2z * mass2 * mass2;

  switch (axisChoice) {

  case LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L:
    //printf("*** OrbitalL ***\n");
    initLNh[0] = 0.;
    initLNh[1] = 0.;
    initLNh[2] = 1.;
    inc = iota;
    break;

  case LAL_SIM_INSPIRAL_FRAME_AXIS_TOTAL_J:
    //printf("*** TotalJ ***\n");
    for (j=0;j<3;j++) {
      iS1[j] = initS1[j];
      iS2[j] = initS2[j];
      initJ[j] = iS1[j] + iS2[j];
      initS1[j] = initS2[j]=0.;
      initLNh[j] = 0.;
    }
    initJ[2] += LNhmag;
    initJmag = sqrt(initJ[0] * initJ[0] + initJ[1] * initJ[1] + initJ[2] * initJ[2]);
    if (initJ[1])
      phiJ = atan2(initJ[1], initJ[0]);
    else
      phiJ = 0.;
    thetaJ = acos(initJ[2]/initJmag);
    rz[0][0] = -cos(phiJ);
    rz[0][1] = -sin(phiJ);
    rz[0][2] = 0.;
    rz[1][0] = sin(phiJ);
    rz[1][1] = -cos(phiJ);
    rz[1][2] = 0.;
    rz[2][0] = 0.;
    rz[2][1] = 0.;
    rz[2][2] = 1.;
    ry[0][0] = cos(thetaJ);
    ry[0][1] = 0;
    ry[0][2] = sin(thetaJ);
    ry[1][0] = 0.;
    ry[1][1] = 1.;
    ry[1][2] = 0.;
    ry[2][0] = -sin(thetaJ);
    ry[2][1] = 0.;
    ry[2][2] = cos(thetaJ);
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
	initLNh[j] += ry[j][k] * rz[k][2];
	for (l = 0; l < 3; l++) {
           initS1[j] += ry[j][k] * rz[k][l] * iS1[l];
           initS2[j] += ry[j][k] * rz[k][l] * iS2[l];
	}
      }
    }
    inc = iota;
    break;

  case LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW:
  default:
    //printf("*** View ***\n");
    initLNh[0] = sin(iota);
    initLNh[1] = 0.;
    initLNh[2] = cos(iota);
    inc = 0.;
    break;
  }

  /*All the PN formulas used in the differential equation integration 
    assume that the spin variables are the physical ones divided by
    totalmasss^2, here we introduce the correct normalization, changing the
    input one, where spin components were normalized on individual mass. */
  for (j = 0; j < 3; j++) {
    initS1[j] /= totalMass * totalMass;
    initS2[j] /= totalMass * totalMass;
  }

  mparams.length = length;

  /* Allocate memory for temporary arrays */

  h2P2 = XLALCreateREAL8Vector(length * 2);
  h2M2 = XLALCreateREAL8Vector(length * 2);
  h2P1 = XLALCreateREAL8Vector(length * 2);
  h2M1 = XLALCreateREAL8Vector(length * 2);
  h20  = XLALCreateREAL8Vector(length * 2);
  h3P3 = XLALCreateREAL8Vector(length * 2);
  h3M3 = XLALCreateREAL8Vector(length * 2);
  h3P2 = XLALCreateREAL8Vector(length * 2);
  h3M2 = XLALCreateREAL8Vector(length * 2);
  h3P1 = XLALCreateREAL8Vector(length * 2);
  h3M1 = XLALCreateREAL8Vector(length * 2);
  h30  = XLALCreateREAL8Vector(length * 2);
  h4P4 = XLALCreateREAL8Vector(length * 2);
  h4M4 = XLALCreateREAL8Vector(length * 2);
  h4P3 = XLALCreateREAL8Vector(length * 2);
  h4M3 = XLALCreateREAL8Vector(length * 2);
  h4P2 = XLALCreateREAL8Vector(length * 2);
  h4M2 = XLALCreateREAL8Vector(length * 2);
  h4P1 = XLALCreateREAL8Vector(length * 2);
  h4M1 = XLALCreateREAL8Vector(length * 2);
  h40  = XLALCreateREAL8Vector(length * 2);
  sigp = XLALCreateREAL8Vector(length);
  sigc = XLALCreateREAL8Vector(length);
  fap  = XLALCreateREAL8Vector(length);
  phap = XLALCreateREAL8Vector(length);

  if (!(h2P2 && h2M2 && h2P1 && h2M1 && h20 && sigp && sigc && fap && phap && h3P3 && h3M3 && h3P2 && h3M2 && h3P1 && h3M1 && h30 && h4P4 && h4M4 && h4P3 && h4M3 && h4P2 && h4M2 && h4P1 && h4M1 && h40) ) {
    if (h2P2)
      XLALDestroyREAL8Vector(h2P2);
    if (h2M2)
      XLALDestroyREAL8Vector(h2M2);
    if (h2P1)
      XLALDestroyREAL8Vector(h2P1);
    if (h2M2)
      XLALDestroyREAL8Vector(h2M1);
    if (h20)
      XLALDestroyREAL8Vector(h20);
    if (h3P3)
      XLALDestroyREAL8Vector(h3P3);
    if (h3M3)
      XLALDestroyREAL8Vector(h3M3);
    if (h3P2)
      XLALDestroyREAL8Vector(h3P2);
    if (h3M2)
      XLALDestroyREAL8Vector(h3M2);
    if (h3P1)
      XLALDestroyREAL8Vector(h3P1);
    if (h3M1)
      XLALDestroyREAL8Vector(h3M1);
    if (h30)
      XLALDestroyREAL8Vector(h30);
    if (h4P4)
      XLALDestroyREAL8Vector(h4P4);
    if (h4M4)
      XLALDestroyREAL8Vector(h4M4);
    if (h4P3)
      XLALDestroyREAL8Vector(h4P3);
    if (h4M3)
      XLALDestroyREAL8Vector(h4M3);
    if (h4P2)
      XLALDestroyREAL8Vector(h4P2);
    if (h4M2)
      XLALDestroyREAL8Vector(h4M2);
    if (h4P1)
      XLALDestroyREAL8Vector(h4P1);
    if (h4M1)
      XLALDestroyREAL8Vector(h4M1);
    if (h40)
      XLALDestroyREAL8Vector(h40);
    if (sigp)
      XLALDestroyREAL8Vector(sigp);
    if (sigc)
      XLALDestroyREAL8Vector(sigc);
    if (fap)
      XLALDestroyREAL8Vector(fap);
    if (phap)
      XLALDestroyREAL8Vector(phap);
    XLAL_ERROR(XLAL_ENOMEM);
  }

  memset(h2P2->data, 0, h2P2->length * sizeof(REAL8));
  memset(h2M2->data, 0, h2M2->length * sizeof(REAL8));
  memset(h2P1->data, 0, h2P1->length * sizeof(REAL8));
  memset(h2M1->data, 0, h2P1->length * sizeof(REAL8));
  memset(h20->data,  0, h20->length  * sizeof(REAL8));
  memset(h3P3->data, 0, h3P3->length * sizeof(REAL8));
  memset(h3M3->data, 0, h3M3->length * sizeof(REAL8));
  memset(h3P2->data, 0, h3P2->length * sizeof(REAL8));
  memset(h3M2->data, 0, h3M2->length * sizeof(REAL8));
  memset(h3P1->data, 0, h3P1->length * sizeof(REAL8));
  memset(h3M1->data, 0, h3M1->length * sizeof(REAL8));
  memset(h30->data,  0, h30->length  * sizeof(REAL8));
  memset(h4P4->data, 0, h3P3->length * sizeof(REAL8));
  memset(h4M4->data, 0, h3M3->length * sizeof(REAL8));
  memset(h4P3->data, 0, h3P3->length * sizeof(REAL8));
  memset(h4M3->data, 0, h3M3->length * sizeof(REAL8));
  memset(h4P2->data, 0, h3P2->length * sizeof(REAL8));
  memset(h4M2->data, 0, h3M2->length * sizeof(REAL8));
  memset(h4P1->data, 0, h3P1->length * sizeof(REAL8));
  memset(h4M1->data, 0, h3M1->length * sizeof(REAL8));
  memset(h40->data,  0, h30->length  * sizeof(REAL8));
  memset(sigp->data, 0, sigp->length * sizeof(REAL8));
  memset(sigc->data, 0, sigc->length * sizeof(REAL8));
  memset(fap->data,  0, fap->length  * sizeof(REAL8));
  memset(phap->data, 0, phap->length * sizeof(REAL8));

  /* Here there used to be a check that OmegaRD is smaller than Nyquist, it
     has been taken out */

    amp22ini = -2.0 * mu * LAL_MRSUN_SI / r * sqrt(16. * LAL_PI / 5.);
 
  /* initialize the coordinates */
  yinit[0] = initphi;     /* phi */
  yinit[1] = initomega;   /* omega (really pi M f) */
  yinit[2] = initLNh[0];   /* LNh(x,y,z) */
  yinit[3] = initLNh[1];
  yinit[4] = initLNh[2];

  yinit[5] = initS1[0];   /* Spin1(x,y,z) */
  yinit[6] = initS1[1];
  yinit[7] = initS1[2];

  yinit[8] = initS2[0];   /* Spin2(x,y,z) */
  yinit[9] = initS2[1];
  yinit[10]= initS2[2];

  yinit[11]= 0.;

  phenPars.intreturn = 0;
  phenPars.energy    = 0.;
  phenPars.omega     = 0.;
  phenPars.domega    = 0.;
  phenPars.ddomega   = 0.;
  phenPars.diota     = 0.;
  phenPars.ddiota    = 0.;
  phenPars.dalpha    = 0.;
  phenPars.ddalpha   = 0.;
  phenPars.countback = 0;
  phenPars.endtime   = 0.;
  phenPars.Psi       = 0.;
  phenPars.alpha     = 0.;
  phenPars.ci        = 0.;
  phenPars.LNhS1     = 0.;
  phenPars.LNhS2     = 0.;
  phenPars.S1S1      = 0.;
  phenPars.S1S2      = 0.;
  phenPars.S2S2      = 0.;

  errcode = XLALSpinInspiralAdaptiveEngine(neqs,yinit,amp22ini,&mparams,h2P2,h2M2,h2P1,h2M1,h20,h3P3,h3M3,h3P2,h3M2,h3P1,h3M1,h30,h4P4,h4M4,h4P3,h4M3,h4P2,h4M2,h4P1,h4M1,h40,fap,phap,&phenPars);
  if(errcode) XLAL_ERROR(XLAL_EFUNC);
  intreturn=phenPars.intreturn;

  count= phenPars.countback;

  /* report on abnormal termination:
     Termination is fine if omegamatch is passed or if energy starts 
     increasing  */

  if ( (intreturn != LALPSIRDPN_TEST_OMEGAMATCH) && (intreturn != LALPSIRDPN_TEST_ENERGY) )
    {
      XLALPrintWarning("** LALPSpinInspiralRD WARNING **: integration terminated with code %d.\n",intreturn);
      fprintf(stderr,"  1025: Energy increases\n  1026: Omegadot -ve\n  1028: Omega NAN\n  1029: Omega > Omegamatch\n  1031: Omega -ve\n");
      XLALPrintWarning("  Waveform parameters were m1 = %14.6e, m2 = %14.6e, inc = %10.6f,\n", mass1, mass2, iota);
      XLALPrintWarning("                           S1 = (%10.6f,%10.6f,%10.6f)\n", s1x, s1y, s1z);
      XLALPrintWarning("                           S2 = (%10.6f,%10.6f,%10.6f)\n", s2x, s2y, s2z);
    }

  if (intreturn==LALPSIRDPN_TEST_OMEGAMATCH) {

    tim = t0 = phenPars.endtime;
    tAs = t0 + 2. * phenPars.domega / phenPars.ddomega;
    om1 = phenPars.domega * tAs * (1. - t0 / tAs) * (1. - t0 / tAs);
    om0 = phenPars.omega - om1 / (1. - t0 / tAs);
    om  = phenPars.omega;

    //diota1 = phenPars.ddiota * tAs * (1. - t0 / tAs) * (1. - t0 / tAs);
    //diota0 = phenPars.diota - diota1 / (1. - t0 / tAs);

    dalpha1 = phenPars.ddalpha * tAs * (1. - t0 / tAs) * (1. - t0 / tAs);
    dalpha0 = phenPars.dalpha - dalpha1 / (1. - t0 / tAs);

    //printf("time %12.6e  count %d\n",tim,phenPars.count);

    if ((tAs < t0) || (om1 < 0.)) {
      XLALPrintError("**** LALPSpinInspiralRD ERROR ****: Could not attach phen part for:\n");
      XLALPrintError(" tAs %12.6e  t0 %12.6e  om1 %12.6e\n",tAs,t0,om1);
      XLALPrintError("   m1 = %14.6e, m2 = %14.6e, inc = %10.6f,\n", mass1, mass2, iota);
      XLALPrintError("   S1 = (%10.6f,%10.6f,%10.6f)\n", s1x, s1y, s1z);
      XLALPrintError("   S2 = (%10.6f,%10.6f,%10.6f)\n", s2x, s2y, s2z);
      XLALDestroyREAL8Vector(h2P2);
      XLALDestroyREAL8Vector(h2M2);
      XLALDestroyREAL8Vector(h2P1);
      XLALDestroyREAL8Vector(h2M1);
      XLALDestroyREAL8Vector(h20);
      XLALDestroyREAL8Vector(h3P3);
      XLALDestroyREAL8Vector(h3M3);
      XLALDestroyREAL8Vector(h3P2);
      XLALDestroyREAL8Vector(h3M2);
      XLALDestroyREAL8Vector(h3P1);
      XLALDestroyREAL8Vector(h3M1);
      XLALDestroyREAL8Vector(h30);
      XLALDestroyREAL8Vector(h4P4);
      XLALDestroyREAL8Vector(h4M4);
      XLALDestroyREAL8Vector(h4P3);
      XLALDestroyREAL8Vector(h4M3);
      XLALDestroyREAL8Vector(h4P2);
      XLALDestroyREAL8Vector(h4M2);
      XLALDestroyREAL8Vector(h4P1);
      XLALDestroyREAL8Vector(h4M1);
      XLALDestroyREAL8Vector(h40);
      XLALDestroyREAL8Vector(fap);
      XLALDestroyREAL8Vector(phap);
      XLALDestroyREAL8Vector(sigp);
      XLALDestroyREAL8Vector(sigc);
      XLAL_ERROR(XLAL_EFAILED);
    }
    else {
      trigAngle.ci = phenPars.ci;
      Psi    = phenPars.Psi;// - 2. * om * log(om);
      Psi0   = Psi + tAs * (om1/mass -dalpha1*trigAngle.ci) * log(1. - t0 / tAs);
      alpha0 = phenPars.alpha + tAs * dalpha1 * log(1. - t0 / tAs);
      //iota0  = acos(phenPars.ci) + diota1 * tAs * log(1. - t0 / tAs);
      energy = phenPars.energy;
      count = phenPars.countback;

      /* Get QNM frequencies */
      errcode = XLALSimIMRPSpinFinalMassSpin(&finalMass,&finalSpin,m1,m2,s1x,s1y,s1z,s2x,s2y,s2z,energy,initLNh[0],initLNh[1],initLNh[2]);
      modefreqs=XLALCreateCOMPLEX8Vector(nmodes);
      errcode+=XLALPSpinGenerateQNMFreq(modefreqs, 2, 2, nmodes, finalMass, finalSpin, totalMass);
      if (errcode != XLAL_SUCCESS) {
        XLALPrintError("**** LALPhenSpinInspiralRD ERROR ****: impossible to generate RingDown frequency\n");
        XLALPrintError( "   m  (%11.4e  %11.4e)  f0 %11.4e\n",mass1, mass2, f_min);
        XLALPrintError( "   S1 (%8.4f  %8.4f  %8.4f)\n", initS1[0],initS1[1], initS1[2]);
        XLALPrintError( "   S2 (%8.4f  %8.4f  %8.4f)\n", initS2[0],initS2[1], initS2[2]);
        XLALDestroyREAL8Vector(h2P2);
        XLALDestroyREAL8Vector(h2M2);
        XLALDestroyREAL8Vector(h2P1);
        XLALDestroyREAL8Vector(h2M1);
        XLALDestroyREAL8Vector(h20);
        XLALDestroyREAL8Vector(h3P3);
        XLALDestroyREAL8Vector(h3M3);
        XLALDestroyREAL8Vector(h3P2);
        XLALDestroyREAL8Vector(h3M2);
        XLALDestroyREAL8Vector(h3P1);
        XLALDestroyREAL8Vector(h3M1);
        XLALDestroyREAL8Vector(h30);
        XLALDestroyREAL8Vector(h4P4);
        XLALDestroyREAL8Vector(h4M4);
        XLALDestroyREAL8Vector(h4P3);
        XLALDestroyREAL8Vector(h4M3);
        XLALDestroyREAL8Vector(h4P2);
        XLALDestroyREAL8Vector(h4M2);
        XLALDestroyREAL8Vector(h4P1);
        XLALDestroyREAL8Vector(h4M1);
        XLALDestroyREAL8Vector(h40);
        XLALDestroyREAL8Vector(fap);
        XLALDestroyREAL8Vector(phap);
        XLALDestroyREAL8Vector(sigp);
        XLALDestroyREAL8Vector(sigc);
        XLALDestroyCOMPLEX8Vector(modefreqs);
        XLAL_ERROR(XLAL_EFAILED);
      }

      omegaRD = creal(modefreqs->data[0]) * unitHz / LAL_PI / 2.;
      frOmRD = fracRD(phenPars.LNhS1,phenPars.LNhS2,phenPars.S1S1,phenPars.S1S2,phenPars.S2S2)*omegaRD;

      v     = cbrt(om);
      v2    = v*v;
      amp22 = amp22ini*v2;

      do {
        count++;
        if (count >= length) {
          XLALPrintError("**** LALPhenSpinInspiralRD ERROR ****: phen. part exceeds array length");
          XLALPrintError( "   m  (%11.4e  %11.4e)  f0 %11.4e\n",mass1, mass2, f_min);
          XLALPrintError( "   S1 (%8.4f  %8.4f  %8.4f)\n", initS1[0],initS1[1], initS1[2]);
          XLALPrintError( "   S2 (%8.4f  %8.4f  %8.4f)\n", initS2[0],initS2[1], initS2[2]);
          XLALDestroyREAL8Vector(h2P2);
          XLALDestroyREAL8Vector(h2M2);
          XLALDestroyREAL8Vector(h2P1);
          XLALDestroyREAL8Vector(h2M1);
          XLALDestroyREAL8Vector(h20);
          XLALDestroyREAL8Vector(h3P3);
          XLALDestroyREAL8Vector(h3M3);
          XLALDestroyREAL8Vector(h3P2);
          XLALDestroyREAL8Vector(h3M2);
          XLALDestroyREAL8Vector(h3P1);
          XLALDestroyREAL8Vector(h3M1);
          XLALDestroyREAL8Vector(h30);
          XLALDestroyREAL8Vector(h4P4);
          XLALDestroyREAL8Vector(h4M4);
          XLALDestroyREAL8Vector(h4P3);
          XLALDestroyREAL8Vector(h4M3);
          XLALDestroyREAL8Vector(h4P2);
          XLALDestroyREAL8Vector(h4M2);
          XLALDestroyREAL8Vector(h4P1);
          XLALDestroyREAL8Vector(h4M1);
          XLALDestroyREAL8Vector(h40);
          XLALDestroyREAL8Vector(fap);
          XLALDestroyREAL8Vector(phap);
          XLALDestroyREAL8Vector(sigp);
          XLALDestroyREAL8Vector(sigc);
          XLALDestroyCOMPLEX8Vector(modefreqs);
          XLAL_ERROR(XLAL_ENOMEM);
        }

        tim += dt;
        v2old = v2;
        //omold = om;
        om = om1 / (1. - tim / tAs) + om0;
        fap->data[count] = om;
        Psi  = Psi0 + (- tAs * (om1/mass-dalpha1*trigAngle.ci) * log(1. - tim / tAs) + (om0/mass-dalpha0*trigAngle.ci) * (tim - t0) );
        //trigAngle.ci = cos(diota0 * (tim - t0) - diota1 * tAs * log(1. - tim / tAs) + iota0);
        alpha = alpha0 + ( dalpha0 * (tim - t0) - dalpha1 * tAs * log(1. - tim / tAs) );

        v  = cbrt(om);
        v2 = v*v;
        amp22 *= v2 / v2old;

        amp33 = -amp22 / 4. * sqrt(5. / 42.);
        amp44 = amp22 * sqrt(5./7.) * 2./9.   * v2;

        trigAngle.cdi  = 2. * trigAngle.ci * trigAngle.ci - 1.;
        trigAngle.c2i  = trigAngle.ci * trigAngle.ci;
        trigAngle.s2i  = 1. - trigAngle.ci * trigAngle.ci;
        trigAngle.si   = sqrt(trigAngle.s2i);
        trigAngle.sdi  = 2. * trigAngle.ci * trigAngle.si;
        trigAngle.c2i2 = (1. + trigAngle.ci) / 2.;
        trigAngle.s2i2 = (1. - trigAngle.ci) / 2.;
        trigAngle.ci2  = sqrt(trigAngle.c2i2);
        trigAngle.si2  = sqrt(trigAngle.s2i2);
        trigAngle.c3i2 = trigAngle.c2i2 * trigAngle.ci2;
        trigAngle.s3i2 = trigAngle.s2i2 * trigAngle.si2;
        trigAngle.c4i2 = trigAngle.c2i2 * trigAngle.c2i2;
        trigAngle.s4i2 = trigAngle.s2i2 * trigAngle.s2i2;
        trigAngle.c5i2 = trigAngle.c4i2 * trigAngle.ci2;
        trigAngle.s5i2 = trigAngle.s4i2 * trigAngle.si2;
        trigAngle.c6i2 = trigAngle.c4i2 * trigAngle.c2i2;
        trigAngle.s6i2 = trigAngle.s4i2 * trigAngle.s2i2;
        trigAngle.c8i2 = trigAngle.c4i2 * trigAngle.c4i2;
        trigAngle.s8i2 = trigAngle.s4i2 * trigAngle.s4i2;

        errcode=XLALSpinInspiralFillH2Modes(h2P2,h2M2,h2P1,h2M1,h20,count,amp22,v,mparams.eta,mparams.dm,Psi,alpha,&trigAngle);

        errcode += XLALSpinInspiralFillH3Modes(h3P3,h3M3,h3P2,h3M2,h3P1,h3M1,h30,count,amp33,v,mparams.eta,mparams.dm,Psi,alpha,&trigAngle);

        errcode += XLALSpinInspiralFillH4Modes(h4P4,h4M4,h4P3,h4M3,h4P2,h4M2,h4P1,h4M1,h40,count,amp44,v,mparams.eta,mparams.dm,Psi,alpha,&trigAngle);

        fap->data[count] = om;
        phap->data[count] = Psi;

      } while ( (om < frOmRD) && (tim < tAs) );

      XLALDestroyCOMPLEX8Vector(modefreqs);

      /*--------------------------------------------------------------
       * Attach the ringdown waveform to the end of inspiral
       -------------------------------------------------------------*/

      //printf("time %12.6e  count %d\n",tim,count);

      apcount  = count;
      errcode  = XLALPSpinInspiralAttachRingdownWave(h2P2, &apcount, nmodes, 2, 2, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h2P2->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h2M2, &apcount, nmodes, 2, -2, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h2M2->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h2P1, &apcount, nmodes, 2, 1, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h2P1->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h2M1, &apcount, nmodes, 2, -1, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h2M1->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h20, &apcount, nmodes, 2, 0, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h20->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h3P3, &apcount, nmodes, 3, 3, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h3P3->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h3M3, &apcount, nmodes, 3, -3, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h3M3->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h3P2, &apcount, nmodes, 3, 2, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h3P2->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h3M2, &apcount, nmodes, 3, -2, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h3P2->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h3P1, &apcount, nmodes, 3, 1, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h3P1->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h3M1, &apcount, nmodes, 3, -1, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h3M1->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h30, &apcount, nmodes, 3, 0, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h30->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4P4, &apcount, nmodes, 4, 4, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h4P4->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4M4, &apcount, nmodes, 4, -4, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h4M4->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4P3, &apcount, nmodes, 4, 3, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h4P3->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4M3, &apcount, nmodes, 4, -3, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h4M3->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4P2, &apcount, nmodes, 4, 2, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h4P4->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4M2, &apcount, nmodes, 4, -2, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h4M4->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4P1, &apcount, nmodes, 4, 1, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h4P3->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4M1, &apcount, nmodes, 4, -1, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h4M3->data[i] = 0.;

      apcount  = count;
      errcode += XLALPSpinInspiralAttachRingdownWave(h40 , &apcount, nmodes, 4, 0, finalMass, finalSpin, tSampling, totalMass);
      for (i = 2 * apcount; i < 2 * length; i++) h40->data[i] = 0.;

      if (errcode != XLAL_SUCCESS) {
	fprintf(stderr,"**** LALPSpinInspiralRD ERROR ****: impossible to create RingDownWave\n");
	XLALDestroyREAL8Vector(h2P2);
	XLALDestroyREAL8Vector(h2M2);
	XLALDestroyREAL8Vector(h2P1);
	XLALDestroyREAL8Vector(h2M1);
	XLALDestroyREAL8Vector(h20);
	XLALDestroyREAL8Vector(h3P3);
	XLALDestroyREAL8Vector(h3M3);
	XLALDestroyREAL8Vector(h3P2);
	XLALDestroyREAL8Vector(h3M2);
	XLALDestroyREAL8Vector(h3P1);
	XLALDestroyREAL8Vector(h3M1);
	XLALDestroyREAL8Vector(h30);
	XLALDestroyREAL8Vector(h4P4);
	XLALDestroyREAL8Vector(h4M4);
	XLALDestroyREAL8Vector(h4P3);
	XLALDestroyREAL8Vector(h4M3);
	XLALDestroyREAL8Vector(h4P2);
	XLALDestroyREAL8Vector(h4M2);
	XLALDestroyREAL8Vector(h4P1);
	XLALDestroyREAL8Vector(h4M1);
	XLALDestroyREAL8Vector(h40);
	XLALDestroyREAL8Vector(fap);
	XLALDestroyREAL8Vector(phap);
	XLAL_ERROR(XLAL_EFAILED);
      }
    }

  } /*End of if intreturn==LALPSIRDPN_TEST_OMEGAMATCH*/

  /*-------------------------------------------------------------------
   * Compute the spherical harmonics required for constructing (h+,hx).
   -------------------------------------------------------------------*/

  /* The angles theta for the spherical harmonics has been set according to 
     the input inclination parameter and the axisChoice */

  for (i = 0; i < length; i++) {
    fap->data[i] /= unitHz;
    sigp->data[i] = 0.;
    sigc->data[i] = 0.;
  }

  MultSphHarmP=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 2, 2);
  MultSphHarmM=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 2, -2);
  for (i = 0; i < length; i++) {
    x0 = h2P2->data[2 * i];
    x1 = h2P2->data[2 * i + 1];
    x2 = h2M2->data[2 * i];
    x3 = h2M2->data[2 * i + 1];
    sigp->data[i] +=   x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP) + x2 * creal(MultSphHarmM) - x3 * cimag(MultSphHarmM);
    sigc->data[i] += - x0 * cimag(MultSphHarmP) - x1 * creal(MultSphHarmP) - x2 * cimag(MultSphHarmM) - x3 * creal(MultSphHarmM);
  }

  MultSphHarmP=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 2, 1);
  MultSphHarmM=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 2, -1);
  for (i = 0; i < length; i++) {
    x0 = h2P1->data[2 * i];
    x1 = h2P1->data[2 * i + 1];
    x2 = h2M1->data[2 * i];
    x3 = h2M1->data[2 * i + 1];
    sigp->data[i] +=   x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP) + x2 * creal(MultSphHarmM) - x3 * cimag(MultSphHarmM);
    sigc->data[i] += - x0 * cimag(MultSphHarmP) - x1 * creal(MultSphHarmP) - x2 * cimag(MultSphHarmM) - x3 * creal(MultSphHarmM);
  }

  MultSphHarmP=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 2, 0);
  for (i = 0; i < length; i++) {
    x0 = h20->data[2 * i];
    x1 = h20->data[2 * i + 1];
    sigp->data[i] += x1 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP);
    sigc->data[i] -= x1 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP);
  }

  MultSphHarmP=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 3, 3);
  MultSphHarmM=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 3, -3);
  for (i = 0; i < length; i++) {
    x0 = h3P3->data[2 * i];
    x1 = h3P3->data[2 * i + 1];
    x2 = h3M3->data[2 * i];
    x3 = h3M3->data[2 * i + 1];
    sigp->data[i] += x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP) + x2 * creal(MultSphHarmM) - x3 * cimag(MultSphHarmM);
    sigc->data[i] -= x0 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP) + x2 * cimag(MultSphHarmM) + x3 * creal(MultSphHarmM);
  }

  MultSphHarmP=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 3, 2);
  MultSphHarmM=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 3, -2);
  for (i = 0; i < length; i++) {
    x0 = h3P2->data[2 * i];
    x1 = h3P2->data[2 * i + 1];
    x2 = h3M2->data[2 * i];
    x3 = h3M2->data[2 * i + 1];
    sigp->data[i] += x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP) + x2 * creal(MultSphHarmM) - x3 * cimag(MultSphHarmM);
    sigc->data[i] -= x0 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP) + x2 * cimag(MultSphHarmM) + x3 * creal(MultSphHarmM);
  }

  MultSphHarmP=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 3, 1);
  MultSphHarmM=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 3, -1);
  for (i = 0; i < length; i++) {
    x0 = h3P1->data[2 * i];
    x1 = h3P1->data[2 * i + 1];
    x2 = h3M1->data[2 * i];
    x3 = h3M1->data[2 * i + 1];
    sigp->data[i] += x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP) + x2 * creal(MultSphHarmM) - x3 * cimag(MultSphHarmM);
    sigc->data[i] -= x0 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP) + x2 * cimag(MultSphHarmM) + x3 * creal(MultSphHarmM);
  }

  MultSphHarmP=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 3, 0);
  for (i = 0; i < length; i++) {
    x0 = h30->data[2 * i];
    x1 = h30->data[2 * i + 1];    
    sigp->data[i] += x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP);
    sigc->data[i] -= x0 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP);
  }

  MultSphHarmP=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 4, 4);
  MultSphHarmM=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 4, -4);
  for (i = 0; i < length; i++) {
    x0 = h4P4->data[2 * i];
    x1 = h4P4->data[2 * i + 1];
    x2 = h4P4->data[2 * i];
    x3 = h4M4->data[2 * i + 1];
    sigp->data[i] += x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP) + x2 * creal(MultSphHarmM) - x3 * cimag(MultSphHarmM);
    sigc->data[i] -= x0 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP) + x2 * cimag(MultSphHarmM) + x3 * creal(MultSphHarmM);
  }

  MultSphHarmP=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 4, 3);
  MultSphHarmM=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 4, -3);
  for (i = 0; i < length; i++) {
    x0 = h4P3->data[2 * i];
    x1 = h4P3->data[2 * i + 1];
    x2 = h4M3->data[2 * i];
    x3 = h4M3->data[2 * i + 1];
    sigp->data[i] += x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP) + x2 * creal(MultSphHarmM) - x3 * cimag(MultSphHarmM);
    sigc->data[i] -= x0 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP) + x2 * cimag(MultSphHarmM) + x3 * creal(MultSphHarmM);
  }

  MultSphHarmP=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 4, 2);
  MultSphHarmM=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 4, -2);
  for (i = 0; i < length; i++) {
    x0 = h4P2->data[2 * i];
    x1 = h4P2->data[2 * i + 1];
    x2 = h4M2->data[2 * i];
    x3 = h4M2->data[2 * i + 1];
    sigp->data[i] += x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP) + x2 * creal(MultSphHarmM) - x3 * cimag(MultSphHarmM);
    sigc->data[i] -= x0 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP) + x2 * cimag(MultSphHarmM) + x3 * creal(MultSphHarmM);
  }

  MultSphHarmP=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 4, 1);
  MultSphHarmM=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 4, -1);
  for (i = 0; i < length; i++) {
    x0 = h4P1->data[2 * i];
    x1 = h4P1->data[2 * i + 1];
    x2 = h4M1->data[2 * i];
    x3 = h4M1->data[2 * i + 1];
    sigp->data[i] += x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP) + x2 * creal(MultSphHarmM) - x3 * cimag(MultSphHarmM);
    sigc->data[i] -= x0 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP) + x2 * cimag(MultSphHarmM) + x3 * creal(MultSphHarmM);
  }

  MultSphHarmP=XLALSpinWeightedSphericalHarmonic(inc, 0., -2, 4, 0);
  for (i = 0; i < length; i++) {
    x0 = h40->data[2 * i];
    x1 = h40->data[2 * i + 1];
    sigp->data[i] += x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP);
    sigc->data[i] -= x0 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP);
  }

  /*------------------------------------------------------
   * If required by the user copy other data sets to the
   * relevant arrays
   ------------------------------------------------------*/

  memcpy((*hplus)->data->data, sigp->data, length * (sizeof(REAL8)));
  memcpy((*hcross)->data->data, sigc->data, length * (sizeof(REAL8)));

  /* Clean up */
  XLALDestroyREAL8Vector(h2P2);
  XLALDestroyREAL8Vector(h2M2);
  XLALDestroyREAL8Vector(h2P1);
  XLALDestroyREAL8Vector(h2M1);
  XLALDestroyREAL8Vector(h20);
  XLALDestroyREAL8Vector(h3P3);
  XLALDestroyREAL8Vector(h3M3);
  XLALDestroyREAL8Vector(h3P2);
  XLALDestroyREAL8Vector(h3M2);
  XLALDestroyREAL8Vector(h3P1);
  XLALDestroyREAL8Vector(h3M1);
  XLALDestroyREAL8Vector(h30);
  XLALDestroyREAL8Vector(h4P4);
  XLALDestroyREAL8Vector(h4M4);
  XLALDestroyREAL8Vector(h4P3);
  XLALDestroyREAL8Vector(h4M3);
  XLALDestroyREAL8Vector(h4P2);
  XLALDestroyREAL8Vector(h4M2);
  XLALDestroyREAL8Vector(h4P1);
  XLALDestroyREAL8Vector(h4M1);
  XLALDestroyREAL8Vector(h40);
  XLALDestroyREAL8Vector(fap);
  XLALDestroyREAL8Vector(phap);
  XLALDestroyREAL8Vector(sigp);
  XLALDestroyREAL8Vector(sigc);

  /* Careful here: count reports the bin number at the attachment of the
     RD, useful for instance to determine tC, it is NOT the actual number
     of waveform points*/
  return count;
  /*End */
}
