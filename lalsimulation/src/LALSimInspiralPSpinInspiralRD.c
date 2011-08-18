/*
 * Copyright (C) 2011 Riccardo Sturani, John Veitch
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
  REAL8 wdotorb[8];           ///< Coefficients of the analytic PN expansion of \f$ \dot\omega_orb\f $
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
  REAL8 OmCutoff;
  REAL8 lengths;
  REAL8 omOffset;
  REAL8 polarization;
  UINT4 length;
  UINT4 inspiralOnly;
} LALPSpinInspiralRDparams;

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

int XLALSimInspiralSpinTaylorCoeffs(REAL8 ST[9], /** Output: spin coefficients array */
		      eta, /** Symmetric mass ratio */ )
{
  if(!ST) XLAL_ERROR(__func__,XLAL_EFAULT);
  
   ST[LAL_PNORDER_NEWTONIAN] = 1.0;
   ST[LAL_PNORDER_HALF] = 0.0;
   ST[LAL_PNORDER_ONE] = ( -(1.0/336.0) * (743.0 + 924.0*eta) );
   ST[LAL_PNORDER_ONE_POINT_FIVE] = ( 4.0 * LAL_PI );
   ST[LAL_PNORDER_TWO] =  ( (34103.0 + 122949.0*eta + 59472.0*eta*eta)/18144.0 );

   ST[LAL_PNORDER_TWO_POINT_FIVE] = ( -(1.0/672.0) * LAL_PI * (4159.0 + 15876.0*eta) );
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

   ST[LAL_PNORDER_THREE] = ( (16447322263.0/139708800.0)
		- (1712.0/105.0)* LAL_GAMMA
		- (273811877.0/1088640.0)*eta - (88.0/3.0)*thetahat*eta
		+ (541.0/896.0)*eta*eta - (5605.0/2592.0)*eta*eta*eta
		+ (1.0/48.0) * LAL_PI*LAL_PI * (256.0 + 451.0*eta)
		- (856.0/105.0)*log(16.0) );
   ST[LAL_PNORDER_THREE+1] = ( -(1712.0/315.0) );     /* extra 3PN component */
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

REAL8 ETa3(REAL8 eta)
{
 return (-675./64. + (209323./4032. - 205.*LAL_PI*LAL_PI/96.
            - 110./9. * lambda)*ieta*eta
            - 155./96. * ieta*eta*eta - 35./5184. * ieta*eta*eta*eta);
}

REAL8 ETa2(REAL8 eta)
{
     return (-(27. - 19*ieta*eta + ieta*eta*eta/3.)/8.);
}

REAL8 ETa1(REAL8 eta)
{
  return((9. + ieta*eta)/12.);
}

REAL8 ETaN(REAL8 eta)
{
  return ( -eta/2.);
}

/***********************************************************************/



/**
 * Convenience function to set up LALPSpinInspiralRDparams struct
 */

int XLALPSpinInspiralRDparamsSetup(LALPSpinInspiralRDparams *mparams, /** Output: RDparams structure */
			     UINT4 inspiralOnly, 	/** Only generate inspiral */
			     REAL8 deltaT, 		/** sampling interval */
			     REAL8 fLow,		/** Starting frequency */
			     REAL8 fCutoff,		/** CHECKME: Cutoff frequency? */
			     REAL8 m1,			/** Mass 1 */
			     REAL8 m2,			/** Mass 2 */
			     LALSpinInteraction spinInteraction,	/** Spin interaction */
			     LALPNOrder order		/** PN Order in Phase */ )
{
  REAL8 totalMass = m1+m2;
  REAL8 tSampling = 1.0/deltaT;
  REAL8 eta = m1*m2/(totalMass * totalMass);
  REAL8 chirpMass = pow(m1*m2,0.6)/pow(totalMass,0.2);
  REAL8 ST[9]; /* SpinTaylor terms */
  
  XLALSimInspiralSpinTaylorCoeffs(ST,eta);
  
  mparams->inspiralOnly = inspiralOnly;
  mparams->dt           = deltaT;
  mparams->OmCutoff     = fCutoff*totalMass * LAL_MTSUN_SI * (REAL8) LAL_PI;
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

    case LAL_PNORDER_THREE_POINT_FIVE:
      mparams->wdotorb[7] = ST[8];

    case LAL_PNORDER_THREE:
      mparams->epnorb[3] = ETa3(eta);
      mparams->wdotorb[6] = ST[6];
      mparams->wdotorblog = ST[7];
      mparams->wdotspin30S1LNh = -LAL_PI/3. * ( 188. - 151./2./mparams->m1m);
      mparams->wdotspin30S2LNh = -LAL_PI/3. * ( 188. + 151./2./mparams->m2m);

    case LAL_PNORDER_TWO_POINT_FIVE:
      mparams->wdotorb[5] = paramsInit->ak.ST[5];
      mparams->epnspin25S1dotLNh = 8. - 31. / 9. * mparams->eta + (3. - 10. / 3. * mparams->eta) * mparams->m2m1;
      mparams->epnspin25S2dotLNh = 8. - 31. / 9. * mparams->eta + (3. - 10. / 3. * mparams->eta) * mparams->m1m2;
      mparams->wdotspin25S1LNh = -31319. / 1008. + 1159. / 24. * mparams->eta + (-809. / 84. + 281. / 8. * mparams->eta) * mparams->m2m1;
      mparams->wdotspin25S2LNh = -31319. / 1008. + 1159. / 24. * mparams->eta + (-809. / 84. + 281. / 8. * mparams->eta) * mparams->m1m2;
      mparams->S1dot25 = 0.5625 + 1.25 * mparams->eta - mparams->eta * mparams->eta / 24. + mparams->dm * (-0.5625 + 0.625 * mparams->eta);
      mparams->S2dot25 = 0.5625 + 1.25 * mparams->eta - mparams->eta * mparams->eta / 24. - mparams->dm * (-0.5625 + 0.625 * mparams->eta);

    case LAL_PNORDER_TWO:
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

    case LAL_PNORDER_ONE_POINT_FIVE:
      mparams->wdotorb[3] = ST[3];
      mparams->epnspin15S1dotLNh = 8. / 3. + 2. * mparams->m2m1;
      mparams->epnspin15S2dotLNh = 8. / 3. + 2. * mparams->m1m2;
      mparams->wdotspin15S1LNh = -(113.0 + 75.0 * mparams->m2m1) / 12.0;
      mparams->wdotspin15S2LNh = -(113.0 + 75.0 * mparams->m1m2) / 12.0;
      mparams->LNhdot15 = 0.5;
      mparams->S1dot15 = (4.0 + 3.0 * mparams->m2m1) / 2.0 * mparams->eta;
      mparams->S2dot15 = (4.0 + 3.0 * mparams->m1m2) / 2.0 * mparams->eta;

    case LAL_PNORDER_ONE:
      mparams->epnorb[1] = ETa1(eta);
      mparams->wdotorb[2] = ST[2];

    case LAL_PNORDER_HALF:
      mparams->wdotorb[1] = ST[1];

    case LAL_PNORDER_NEWTONIAN:
      mparams->epnorb[0] = ETaN(eta);
      mparams->wdotorb[0] = ST[0];
      break;

    case LAL_PNORDER_PSEUDO_FOUR:
      XLALPrintError("*** LALPhenSpinInspiralRD ERROR: PhenSpin approximant not available at pseudo4PN order\n");
			XLAL_ERROR(__func__,XLAL_EDOM);
      break;

    case LAL_PNORDER_NUM_ORDER:
      XLALPrintError("*** LALPhenSpinInspiralRD ERROR: NUM_ORDER not a valid PN order\n");
			XLAL_ERROR(__func__,XLAL_EDOM);
			break;

    default:
      XLALPrintError("*** LALPhenSpinInspiralRD ERROR: Impossible to create waveform with %d order\n",params->order);
			XLAL_ERROR(__func__,XLAL_EFAILED);
      break;
  }

  switch (params->spinInteraction) {

  case LAL_NOInter:
    mparams->wdotspin30S1LNh   = 0.;
    mparams->wdotspin30S2LNh   = 0.;
    mparams->epnspin25S1dotLNh = 0.;
    mparams->epnspin25S2dotLNh = 0.;
    mparams->wdotspin25S1LNh   = 0.;
    mparams->wdotspin25S2LNh   = 0.;
    mparams->S1dot25           = 0.;
    mparams->S2dot25           = 0.;
    mparams->epnspin15S1dotLNh = 0.;
    mparams->epnspin15S2dotLNh = 0.;
    mparams->wdotspin15S1LNh   = 0.;
    mparams->wdotspin15S2LNh   = 0.;
    mparams->S1dot15           = 0.;
    mparams->S2dot15           = 0.;

  case LAL_SOInter:
    mparams->wdotspin20S1S2      = 0.;
    mparams->epnspin20S1S2       = 0.;
    mparams->epnspin20S1S2dotLNh = 0.;

  case LAL_SSInter:
    mparams->wdotspin20S1S1 = 0.;
    mparams->epnspin20S1S1 = 0.;
    mparams->epnspin20S2S2 = 0.;
    mparams->Sdot20S = 0.;
    mparams->epnspin20S1S1 = 0.;
    mparams->epnspin20S2S2 = 0.;
    mparams->epnspin20S1S1dotLNh = 0.;
    mparams->epnspin20S2S2dotLNh = 0.;
    break;

  case LAL_SSselfInter:
    break;
  case LAL_QMInter:
    break;
  case LAL_AllInter:
    break;
  default:
    break;
  }
  return(XLAL_SUCCESS);
}


INT4 XLALSpinInspiralDerivatives(double t, const double values[], double dvalues[], void *mparams) {

    REAL8 omega;                // time-derivative of the orbital phase
    REAL8 LNhx, LNhy, LNhz;     // orbital angolar momentum unit vector
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

    UNUSED(t);

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
	//fprintf(stderr,"*** LALPSpinInspiralRD WARNING ***: alphadot set to 0 by hand LNh:(%12.4e %12.4e %12.4e)\n",LNhx,LNhy,LNhz);
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

    REAL8 omega;                // time-derivative of the orbital phase
    REAL8 LNhx, LNhy, LNhz;     // orbital angolar momentum unit vector
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

    UNUSED(t);

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
	//fprintf(stderr,"*** LALPSpinInspiralRD WARNING ***: alphadot set to 0 by hand LNh:(%12.4e %12.4e %12.4e)\n",LNhx,LNhy,LNhz);
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


void LALSpinInspiralDerivatives(REAL8Vector * values, REAL8Vector * dvalues, void *mparams)
{
  XLALSpinInspiralDerivatives(0.,values->data,dvalues->data,mparams);
}				/* end of LALSpinInspiralDerivatives */


/**
 * Driver routine to compute the PhenSpin Inspiral Ringdown waveform.
 *
 */
int XLALSimInspiralPSpinInspiralRDGenerator(
		REAL8TimeSeries **hplus,  /**< +-polarization waveform */
	       	REAL8TimeSeries **hcross, /**< x-polarization waveform */
	       	LIGOTimeGPS *tc,          /**< coalescence time */
	       	REAL8 phi0,               /**< start phase */
	       	REAL8 x0,                 /**< tail-term gauge choice thing (if you don't know, just set it to zero) */
	       	REAL8 deltaT,             /**< sampling interval */
	       	REAL8 m1,                 /**< mass of companion 1 */
	       	REAL8 m2,                 /**< mass of companion 2 */
	       	REAL8 f_min,              /**< start frequency */
	       	REAL8 r,                  /**< distance of source */
	       	REAL8 iota,                  /**< inclination of source (rad) */
		REAL8 spin1[3],		  /**< Spin vector on mass1 */
		REAL8 spin2[3],		  /**< Spin vector on mass2 */
	       	int amplitudeO,           /**< twice post-Newtonian amplitude order */
	       	int phaseO                /**< twice post-Newtonian phase order */
		);

INT4 XLALPSpinInspiralRDEngine(
			REAL8Vector * hh,
			REAL8Vector * ff,
			REAL8Vector * phi,
			UINT4 *countback,
			InspiralTemplate *params,
			InspiralInit     *paramsInit)
{

  REAL8Vector *signalvec1,*signalvec2;
  if(hplus) if(*hplus) signalvec1=*hplus->data;
  if(hcross) if(*hcross) signalvec2=*hcross->data;
  REAL8Vector *hh, *ff, *phi;
  UINT4 countback;
  
  
  /* declare code parameters and variables */
  const INT4 neqs = 11+1;      // number of dynamical variables plus the energy function
  UINT4 count,apcount;         // integration steps performed
  UINT4 length;                // signal vector length
  UINT4 i, j, k, l;            // counters

  REAL8 v = 0.;
  REAL8 v2 = 0.;
  REAL8 v2old;
  REAL8 mass;                  // Total mass in SI units
  REAL8 tim;                   // time (units of total mass)
  REAL8 unitHz;
  REAL8 initomega,initphi;
  REAL8 inc;
  REAL8 LNhmag,initJmag;
  REAL8 LNhxy;
  REAL8 initS1[3],initS2[3],initLNh[3],initJ[3];
  REAL8 iS1[3],iS2[3];
  REAL8 phiJ,thetaJ;
  REAL8 ry[3][3],rz[3][3];
  REAL8 dt;

  INT4  intreturn;
  REAL8 yinit[neqs];

  LALPSpinInspiralRDparams mparams;

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
  REAL8Vector* hap;
  REAL8Vector* phap;

  LALPSpinInspiralPhenPars phenPars;

  REAL8 Psi=0.;
  REAL8 amp22ini,amp22,amp33,amp44;

  LALSpinInspiralAngle trigAngle;

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
  REAL8 fCutoff = 0.0; /* Set only for inspiral only waveforms */
  REAL8 totalMass = m1+m2;

  if(!params) XLAL_ERROR(__func__, XLAL_EFAULT);

  mass = (m1+m2) * LAL_MTSUN_SI;
  unitHz = (m1+m2) * LAL_MTSUN_SI * (REAL8) LAL_PI;

  if ((signalvec2)||(hh))
    params->nStartPad = 0;    /* must be zero for templates and injection */
  /* -- length in seconds from Newtonian formula; */

  dt = deltaT;

  /* setup coefficients for PN equations */
  if(XLALPSpinInspiralRDparamsSetup(&mparams, /** Output: RDparams structure */
			     inspiralOnly, 	/** Only generate inspiral */
			     deltaT, 		/** sampling interval */
			     f_min,		/** Starting frequency */
			     fCutoff,		/** CHECKME: Cutoff frequency? */
			     m1,			/** Mass 1 */
			     m2,			/** Mass 2 */
			     LAL_AllInter,	/** Spin interaction */
			     phaseO		/** PN Order in Phase */ ))
		XLAL_ERROR(__func__,XLAL_EFUNC);

  /* Check that initial frequency is smaller than omegamatch ~ xxyy for m=100 Msun */
  initphi   = phi0;
  initomega = f_min*unitHz;

  /* Check that initial frequency is smaller than omegamatch ~ xxyy for m=100 Msun */

  LNhS1=spin1[2];
  LNhS2=spin2[2];
  S1S1=spin1[0]*spin1[0]+spin1[1]*spin1[1]+spin1[2]*spin1[2];
  S1S2=spin1[0]*spin2[0]+spin1[1]*spin2[1]+spin1[2]*spin2[2];
  S2S2=spin2[0]*spin2[0]+spin2[1]*spin2[1]+spin2[2]*spin2[2];

  omegaMatch = OmMatch(LNhS1,LNhS2,S1S1,S1S2,S2S2);

  if ( initomega > omegaMatch ) {
    /*if ((params->spin1[0]==params->spin1[1])&&(params->spin1[1]==params->spin2[0])&&(params->spin2[0]==params->spin2[1])&&(params->spin2[1]==0.)) {
      //Beware, this correspond to a shift of the initial phase!
      initomega = 0.95*omegaMatch;
      fprintf(stdout,"*** LALPSpinInspiralRD WARNING ***: Initial frequency reset from %12.6e to %12.6e Hz, m:(%12.4e,%12.4e)\n",params->fLower,initomega/unitHz,params->mass1,params->mass2);
      }*/
    /*else {*/
    XLALPrintError("**** LALPSpinInspiralRD ERROR ****: Initial frequency too high: %11.5e for omM ~ %11.5e and m:(%8.3f, %8.3f)\n",f_min,omegaMatch/unitHz,m1,m2);
    XLAL_ERROR(__func__,XLAL_EFAILED);
    /*}*/
  }

  /* Here we use the following convention:
     the coordinates of the spin vectors params->spin1,2 and the params->inclination 
     variable refers to different physical parameters according to the value of 
     params->axisChoice:
     * OrbitalL: params->inclination denotes the angle between the view direction
                 N and the initial L (initial L//z, N in the x-z plane) and the spin 
		 coordinates are given with respect to initial L.
     * View:     params->inclination denotes the angle between the initial L and N 
                 (N//z, initial L in the x-z plane) and the spin coordinates 
		 are given with respect to N.
     * TotalJ:   params->inclination denotes the angle between the view directoin 
                 and J (J is constant during the evolution, J//z, both N and initial 
		 L are in the x-z plane) and the spin coordinates are given wrt 
		 initial ** L **.
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
  for (i = 0; i < 3; i++) {
    initS1[i] = spin1[i] * m1 * m1;
    initS2[i] = spin2[i] * m2 * m2;
  }

  switch (params->axisChoice) {

  case OrbitalL:
    //printf("*** OrbitalL ***\n");
    initLNh[0] = 0.;
    initLNh[1] = 0.;
    initLNh[2] = 1.;
    inc = iota;
    break;

  case View:
    //printf("*** View ***\n");
    initLNh[0] = sin(iota);
    initLNh[1] = 0.;
    initLNh[2] = cos(iota);
    inc = 0.;
    break;

  default:
    //case TotalJ:
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
  }

  if (initS1[0]*initS1[0]+initS1[1]*initS1[1]+initS2[0]*initS2[0]+initS2[1]*initS2[1]) LNhxy=sqrt(initLNh[0]*initLNh[0]+initLNh[1]*initLNh[1]);
  else LNhxy=1.;

  /*All the PN formulas used in the differential equation integration 
    assume that the spin variables are the physical ones divided by
    totalmasss^2, here we introduce the correct normalization, changing the
    input one, where spin components were normalized on individual mass. */
  for (j = 0; j < 3; j++) {
    initS1[j] /= totalMass * totalMass;
    initS2[j] /= totalMass * totalMass;
  }

  if (signalvec1) {
    length = signalvec1->length;
  } else {
    if (ff)
      length = ff->length;
    else
      length = 0;
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
  hap  = XLALCreateREAL8Vector(length * 2);
  fap  = XLALCreateREAL8Vector(length);
  phap = XLALCreateREAL8Vector(length);

  if (!h2P2 || !h2M2 || !h2P1 || !h2M1 || !h20 || !sigp || !sigc || !fap || !phap || !hap || !h3P3 || !h3M3 || !h3P2 || !h3M2 || !h3P1 || !h3M1 || !h30 || !h4P4 || !h4M4 || !h4P3 || !h4M3 || !h4P2 || !h4M2 || !h4P1 || !h4M1 || !h40 ) {
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
    if (hap)
      XLALDestroyREAL8Vector(hap);
    if (phap)
      XLALDestroyREAL8Vector(phap);
    XLAL_ERROR(__func__,XLAL_ENOMEM);
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
  memset(hap->data,  0, hap->length  * sizeof(REAL8));
  memset(fap->data,  0, fap->length  * sizeof(REAL8));
  memset(phap->data, 0, phap->length * sizeof(REAL8));

  /* Here there used to be a check that OmegaRD is smaller than Nyquist, it
     has been taken out */

    amp22ini = -2.0 * params->mu * LAL_MRSUN_SI / r * sqrt(16. * LAL_PI / 5.);
 
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

  if (params->fixedStep == 1) {
    //printf("Fixed step integration\n");
   if( XLALSpinInspiralEngine(neqs,yinit,amp22ini,&mparams,h2P2,h2M2,h2P1,h2M1,h20,h3P3,h3M3,h3P2,h3M2,h3P1,h3M1,h30,h4P4,h4M4,h4P3,h4M3,h4P2,h4M2,h4P1,h4M1,h40,fap,phap,&phenPars))
		 XLAL_ERROR(__func__,XLAL_EFUNC);
  }
	else {
    //printf("Adaptive integration\n");
    errcode = XLALSpinInspiralAdaptiveEngine(neqs,yinit,amp22ini,&mparams,h2P2,h2M2,h2P1,h2M1,h20,h3P3,h3M3,h3P2,h3M2,h3P1,h3M1,h30,h4P4,h4M4,h4P3,h4M3,h4P2,h4M2,h4P1,h4M1,h40,fap,phap,&phenPars);
		if(errcode) XLAL_ERROR(__func__,XLAL_EFUNC);
  }
  intreturn=phenPars.intreturn;

  /* report on abnormal termination:
     Termination is fine if omegamatch is passed or if energy starts 
     increasing  */

  if ( (intreturn!=LALPSIRDPN_TEST_OMEGACUT) && (intreturn != LALPSIRDPN_TEST_OMEGAMATCH) && (intreturn != LALPSIRDPN_TEST_ENERGY) )
    {
      fprintf(stderr,"** LALPSpinInspiralRD WARNING **: integration terminated with code %d.\n",intreturn);
      fprintf(stderr,"  1025: Energy increases\n  1026: Omegadot -ve\n  1028: Omega NAN\n  1029: Omega > Omegamatch\n  1031: Omega -ve\n  1032: Omega > OmegaCut %12.6e\n",mparams.OmCutoff);
      fprintf(stderr,"  Waveform parameters were m1 = %14.6e, m2 = %14.6e, inc = %10.6f,\n", params->mass1, params->mass2, params->inclination);
      fprintf(stderr,"                           S1 = (%10.6f,%10.6f,%10.6f)\n", params->spin1[0], params->spin1[1], params->spin1[2]);
      fprintf(stderr,"                           S2 = (%10.6f,%10.6f,%10.6f)\n", params->spin2[0], params->spin2[1], params->spin2[2]);
    }

  if ((params->inspiralOnly != 1)&&(intreturn==LALPSIRDPN_TEST_OMEGAMATCH)) {

    tim = t0 = phenPars.endtime;
    tAs = t0 + 2. * phenPars.domega / phenPars.ddomega;
    om1 = phenPars.domega * tAs * (1. - t0 / tAs) * (1. - t0 / tAs);
    om0 = phenPars.omega - om1 / (1. - t0 / tAs);
    om  = phenPars.omega;

    //diota1 = phenPars.ddiota * tAs * (1. - t0 / tAs) * (1. - t0 / tAs);
    //diota0 = phenPars.diota - diota1 / (1. - t0 / tAs);

    dalpha1 = phenPars.ddalpha * tAs * (1. - t0 / tAs) * (1. - t0 / tAs);
    dalpha0 = phenPars.dalpha - dalpha1 / (1. - t0 / tAs);

    //printf("time %12.6e  count %d\n",tim,phenPars.countback);

    if ((tAs < t0) || (om1 < 0.)) {
      XLALPrintError("**** LALPSpinInspiralRD ERROR ****: Could not attach phen part for:\n");
      XLALPrintError(" tAs %12.6e  t0 %12.6e  om1 %12.6e\n",tAs,t0,om1);
      XLALPrintError("   m1 = %14.6e, m2 = %14.6e, inc = %10.6f,\n", params->mass1, params->mass2, params->inclination);
      XLALPrintError("   S1 = (%10.6f,%10.6f,%10.6f)\n", params->spin1[0], params->spin1[1], params->spin1[2]);
      XLALPrintError("   S2 = (%10.6f,%10.6f,%10.6f)\n", params->spin2[0], params->spin2[1], params->spin2[2]);
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
      XLALDestroyREAL8Vector(hap);
      XLALDestroyREAL8Vector(fap);
      XLALDestroyREAL8Vector(phap);
      XLALDestroyREAL8Vector(sigp);
      XLALDestroyREAL8Vector(sigc);
      XLAL_ERROR(__func__, XLAL_EFAILED);
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
      errcode = XLALPSpinFinalMassSpin(&finalMass, &finalSpin, params, energy, initLNh);
      modefreqs=XLALCreateCOMPLEX8Vector(nmodes);
      errcode+=XLALPSpinGenerateQNMFreq(modefreqs, params, 2, 2, nmodes, finalMass, finalSpin);
      if (errcode != XLAL_SUCCESS) {
	XLALPrintError("**** LALPhenSpinInspiralRD ERROR ****: impossible to generate RingDown frequency\n");
	XLALPrintError( "   m  (%11.4e  %11.4e)  f0 %11.4e\n",params->mass1, params->mass2, params->fLower);
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
	XLALDestroyREAL8Vector(hap);
	XLALDestroyREAL8Vector(sigp);
	XLALDestroyREAL8Vector(sigc);
	XLALDestroyCOMPLEX8Vector(modefreqs);
	XLAL_ERROR(__func__,XLAL_EFAILED);
      }

      omegaRD = modefreqs->data[0].re * unitHz / LAL_PI / 2.;
      frOmRD = fracRD(phenPars.LNhS1,phenPars.LNhS2,phenPars.S1S1,phenPars.S1S2,phenPars.S2S2)*omegaRD;

      v     = cbrt(om);
      v2    = v*v;
      amp22 = amp22ini*v2;

      do {

	count++;
	if (count >= length) {
          XLALPrintError("**** LALPhenSpinInspiralRD ERROR ****: phen. part exceeds array length");
          XLALPrintError( "   m  (%11.4e  %11.4e)  f0 %11.4e\n",params->mass1, params->mass2, params->fLower);
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
          XLALDestroyREAL8Vector(hap);
          XLALDestroyREAL8Vector(fap);
          XLALDestroyREAL8Vector(phap);
          XLALDestroyREAL8Vector(sigp);
          XLALDestroyREAL8Vector(sigc);
          XLALDestroyCOMPLEX8Vector(modefreqs);
          XLAL_ERROR(__func__,XLAL_ENOMEM);
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

	errcode=XLALSpinInspiralFillH2Modes(h2P2,h2M2,h2P1,h2M1,h20,count,amp22,v,mparams.eta,mparams.dm,Psi,alpha,trigAngle);

	errcode += XLALSpinInspiralFillH3Modes(h3P3,h3M3,h3P2,h3M2,h3P1,h3M1,h30,count,amp33,v,mparams.eta,mparams.dm,Psi,alpha,trigAngle);

	errcode += XLALSpinInspiralFillH4Modes(h4P4,h4M4,h4P3,h4M3,h4P2,h4M2,h4P1,h4M1,h40,count,amp44,v,mparams.eta,mparams.dm,Psi,alpha,trigAngle);

	fap->data[count] = om;
	phap->data[count] = Psi;

      } while ( (om < frOmRD) && (tim < tAs) );

      XLALDestroyCOMPLEX8Vector(modefreqs);
      *countback=count;

      /*--------------------------------------------------------------
       * Attach the ringdown waveform to the end of inspiral
       -------------------------------------------------------------*/

      //printf("time %12.6e  count %d\n",tim,count);

      apcount  = *countback;
      errcode  = XLALPSpinInspiralAttachRingdownWave(h2P2, params, &apcount, nmodes, 2, 2, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h2P2->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h2M2, params, &apcount, nmodes, 2, -2, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h2M2->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h2P1, params, &apcount, nmodes, 2, 1, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h2P1->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h2M1, params, &apcount, nmodes, 2, -1, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h2M1->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h20, params, &apcount, nmodes, 2, 0, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h20->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h3P3, params, &apcount, nmodes, 3, 3, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h3P3->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h3M3, params, &apcount, nmodes, 3, -3, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h3M3->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h3P2, params, &apcount, nmodes, 3, 2, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h3P2->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h3M2, params, &apcount, nmodes, 3, -2, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h3P2->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h3P1, params, &apcount, nmodes, 3, 1, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h3P1->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h3M1, params, &apcount, nmodes, 3, -1, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h3M1->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h30, params, &apcount, nmodes, 3, 0, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h30->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4P4, params, &apcount, nmodes, 4, 4, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h4P4->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4M4, params, &apcount, nmodes, 4, -4, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h4M4->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4P3, params, &apcount, nmodes, 4, 3, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h4P3->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4M3, params, &apcount, nmodes, 4, -3, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h4M3->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4P2, params, &apcount, nmodes, 4, 2, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h4P4->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4M2, params, &apcount, nmodes, 4, -2, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h4M4->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4P1, params, &apcount, nmodes, 4, 1, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h4P3->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4M1, params, &apcount, nmodes, 4, -1, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h4M3->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = *countback;
      errcode += XLALPSpinInspiralAttachRingdownWave(h40, params, &apcount, nmodes, 4, 0, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h40->data[i] = 0.;
      if (apcount > count) count = apcount;

      *countback=count;

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
	XLALDestroyREAL8Vector(hap);
	XLALDestroyREAL8Vector(fap);
	XLALDestroyREAL8Vector(phap);
	XLAL_ERROR(__func__,XLAL_EFAILED);
      }
    }

  } /*End of if not inspiralonly and test_omegamatch*/

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

  errcode  = XLALSphHarm(&MultSphHarmP, 2, 2, inc, 0.);
  errcode += XLALSphHarm(&MultSphHarmM, 2, -2, inc, 0.);
  if (errcode != XLAL_SUCCESS) {
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
    XLALDestroyREAL8Vector(hap);
    XLALDestroyREAL8Vector(fap);
    XLALDestroyREAL8Vector(phap);
    XLALPrintError("**** LALPSpinInspiralRD ERROR ****: impossible to create Y22 or Y2-2\n");
    XLAL_ERROR(__func__,XLAL_EFAILED);
  }
  for (i = 0; i < length; i++) {
    x0 = h2P2->data[2 * i];
    x1 = h2P2->data[2 * i + 1];
    x2 = h2M2->data[2 * i];
    x3 = h2M2->data[2 * i + 1];
    sigp->data[i] +=   x0 * MultSphHarmP.re - x1 * MultSphHarmP.im + x2 * MultSphHarmM.re - x3 * MultSphHarmM.im;
    sigc->data[i] += - x0 * MultSphHarmP.im - x1 * MultSphHarmP.re - x2 * MultSphHarmM.im - x3 * MultSphHarmM.re;
  }

  errcode  = XLALSphHarm(&MultSphHarmP, 2, 1, inc, 0.);
  errcode += XLALSphHarm(&MultSphHarmM, 2, -1, inc, 0.);
  if (errcode != XLAL_SUCCESS){
    XLALDestroyREAL8Vector(h2P1);
    XLALDestroyREAL8Vector(h2M1);
    XLALPrintWarning("** LALPSpinInspiralRD WARNING **: impossible to create Y21\n");
  } else {
    for (i = 0; i < length; i++) {
      x0 = h2P1->data[2 * i];
      x1 = h2P1->data[2 * i + 1];
      x2 = h2M1->data[2 * i];
      x3 = h2M1->data[2 * i + 1];
      sigp->data[i] +=   x0 * MultSphHarmP.re - x1 * MultSphHarmP.im + x2 * MultSphHarmM.re - x3 * MultSphHarmM.im;
      sigc->data[i] += - x0 * MultSphHarmP.im - x1 * MultSphHarmP.re - x2 * MultSphHarmM.im - x3 * MultSphHarmM.re;
    }
  }

  errcode = XLALSphHarm(&MultSphHarmP, 2, 0, inc, 0.);
  if (errcode != XLAL_SUCCESS) {
    XLALDestroyREAL8Vector(h20);
    XLALPrintWarning("** LALPSpinInspiralRD WARNING **: impossible to create Y20\n");
  } else {
    for (i = 0; i < length; i++) {
      x0 = h20->data[2 * i];
      x1 = h20->data[2 * i + 1];
      sigp->data[i] += x1 * MultSphHarmP.re - x1 * MultSphHarmP.im;
      sigc->data[i] -= x1 * MultSphHarmP.im + x1 * MultSphHarmP.re;
    }
  }

  errcode  = XLALSphHarm(&MultSphHarmP, 3, 3, inc, 0.);
  errcode += XLALSphHarm(&MultSphHarmM, 3, -3, inc, 0.);

  if (errcode != XLAL_SUCCESS) {
    XLALDestroyREAL8Vector(h3P3);
    XLALDestroyREAL8Vector(h3M3);
    XLALPrintWarning("** LALPSpinInspiralRD WARNING **: impossible to create Y33,Y3-3\n");
  } else {
    for (i = 0; i < length; i++) {
      x0 = h3P3->data[2 * i];
      x1 = h3P3->data[2 * i + 1];
      x2 = h3M3->data[2 * i];
      x3 = h3M3->data[2 * i + 1];
      sigp->data[i] += x0 * MultSphHarmP.re - x1 * MultSphHarmP.im + x2 * MultSphHarmM.re - x3 * MultSphHarmM.im;
      sigc->data[i] -= x0 * MultSphHarmP.im + x1 * MultSphHarmP.re + x2 * MultSphHarmM.im + x3 * MultSphHarmM.re;
    }
  }

  errcode  = XLALSphHarm(&MultSphHarmP, 3, 2, inc, 0.);
  errcode += XLALSphHarm(&MultSphHarmM, 3, -2, inc, 0.);
  if (errcode != XLAL_SUCCESS) {
    XLALDestroyREAL8Vector(h3P2);
    XLALDestroyREAL8Vector(h3M2);
    XLALPrintWarning("** LALPSpinInspiralRD WARNING **: impossible to create Y32,Y3-2\n");
  } else {
    for (i = 0; i < length; i++) {
      x0 = h3P2->data[2 * i];
      x1 = h3P2->data[2 * i + 1];
      x2 = h3M2->data[2 * i];
      x3 = h3M2->data[2 * i + 1];
      sigp->data[i] += x0 * MultSphHarmP.re - x1 * MultSphHarmP.im + x2 * MultSphHarmM.re - x3 * MultSphHarmM.im;
      sigc->data[i] -= x0 * MultSphHarmP.im + x1 * MultSphHarmP.re + x2 * MultSphHarmM.im + x3 * MultSphHarmM.re;
    }
  }

  errcode  = XLALSphHarm(&MultSphHarmP, 3, 1, inc, 0.);
  errcode += XLALSphHarm(&MultSphHarmM, 3, -1, inc, 0.);
  if (errcode != XLAL_SUCCESS) {
    XLALDestroyREAL8Vector(h3P1);
    XLALDestroyREAL8Vector(h3M1);
    XLALPrintWarning("** LALPSpinInspiralRD WARNING **: impossible to create Y31,Y3-1\n");
  } else {
    for (i = 0; i < length; i++) {
      x0 = h3P1->data[2 * i];
      x1 = h3P1->data[2 * i + 1];
      x2 = h3M1->data[2 * i];
      x3 = h3M1->data[2 * i + 1];
      sigp->data[i] += x0 * MultSphHarmP.re - x1 * MultSphHarmM.im + x2 * MultSphHarmM.re - x3 * MultSphHarmM.im;
      sigc->data[i] -= x0 * MultSphHarmP.im + x1 * MultSphHarmP.re + x2 * MultSphHarmM.im + x3 * MultSphHarmM.re;
    }
  }

  errcode  = XLALSphHarm(&MultSphHarmP, 3, 0, inc, 0.);
  if (errcode != XLAL_SUCCESS) {
    XLALDestroyREAL8Vector(h30);
    XLALPrintWarning("** LALPSpinInspiralRD WARNING **: impossible to create Y30\n");
  } else {
    for (i = 0; i < length; i++) {
      x0 = h30->data[2 * i];
      x1 = h30->data[2 * i + 1];    
      sigp->data[i] += x0 * MultSphHarmP.re - x1 * MultSphHarmM.im;
      sigc->data[i] -= x0 * MultSphHarmP.im + x1 * MultSphHarmP.re;
    }
  }

  errcode  = XLALSphHarm(&MultSphHarmP, 4, 4, inc, 0.);
  errcode += XLALSphHarm(&MultSphHarmM, 4, -4, inc, 0.);
  if (errcode != XLAL_SUCCESS) {
    XLALDestroyREAL8Vector(h4P4);
    XLALDestroyREAL8Vector(h4M4);
    XLALPrintWarning("** LALPSpinInspiralRD WARNING **: impossible to create Y44,Y4-4\n");
  } else {
    for (i = 0; i < length; i++) {
      x0 = h4P4->data[2 * i];
      x1 = h4P4->data[2 * i + 1];
      x2 = h4P4->data[2 * i];
      x3 = h4M4->data[2 * i + 1];
      sigp->data[i] += x0 * MultSphHarmP.re - x1 * MultSphHarmP.im + x2 * MultSphHarmM.re - x3 * MultSphHarmM.im;
      sigc->data[i] -= x0 * MultSphHarmP.im + x1 * MultSphHarmP.re + x2 * MultSphHarmM.im + x3 * MultSphHarmM.re;
    }
  }

  errcode  = XLALSphHarm(&MultSphHarmP, 4, 3, inc, 0.);
  errcode += XLALSphHarm(&MultSphHarmM, 4, -3, inc, 0.);
  if (errcode != XLAL_SUCCESS) {
    XLALDestroyREAL8Vector(h4P3);
    XLALDestroyREAL8Vector(h4M3);
    XLALPrintWarning("** LALPSpinInspiralRD WARNING **: impossible to create Y43,Y4-3\n");
  } else {
    for (i = 0; i < length; i++) {
      x0 = h4P3->data[2 * i];
      x1 = h4P3->data[2 * i + 1];
      x2 = h4M3->data[2 * i];
      x3 = h4M3->data[2 * i + 1];
      sigp->data[i] += x0 * MultSphHarmP.re - x1 * MultSphHarmP.im + x2 * MultSphHarmM.re - x3 * MultSphHarmM.im;
      sigc->data[i] -= x0 * MultSphHarmP.im + x1 * MultSphHarmP.re + x2 * MultSphHarmM.im + x3 * MultSphHarmM.re;
    }
  }

  errcode  = XLALSphHarm(&MultSphHarmP, 4, 2, inc, 0.);
  errcode += XLALSphHarm(&MultSphHarmM, 4, -2, inc, 0.);
  if (errcode != XLAL_SUCCESS) {
    XLALDestroyREAL8Vector(h4P2);
    XLALDestroyREAL8Vector(h4M2);
    XLALPrintWarning("** LALPSpinInspiralRD WARNING **: impossible to create Y42,Y4-2\n");
  } else {
    for (i = 0; i < length; i++) {
      x0 = h4P2->data[2 * i];
      x1 = h4P2->data[2 * i + 1];
      x2 = h4M2->data[2 * i];
      x3 = h4M2->data[2 * i + 1];
      sigp->data[i] += x0 * MultSphHarmP.re - x1 * MultSphHarmP.im + x2 * MultSphHarmM.re - x3 * MultSphHarmM.im;
      sigc->data[i] -= x0 * MultSphHarmP.im + x1 * MultSphHarmP.re + x2 * MultSphHarmM.im + x3 * MultSphHarmM.re;
    }
  }

  errcode  = XLALSphHarm(&MultSphHarmP, 4, 1, inc, 0.);
  errcode += XLALSphHarm(&MultSphHarmM, 4, -1, inc, 0.);
  if (errcode != XLAL_SUCCESS) {
    XLALDestroyREAL8Vector(h4P1);
    XLALDestroyREAL8Vector(h4M1);
    XLALPrintWarning("** LALPSpinInspiralRD WARNING **: impossible to create Y41,Y4-1\n");
  } else {
    for (i = 0; i < length; i++) {
      x0 = h4P1->data[2 * i];
      x1 = h4P1->data[2 * i + 1];
      x2 = h4M1->data[2 * i];
      x3 = h4M1->data[2 * i + 1];
      sigp->data[i] += x0 * MultSphHarmP.re - x1 * MultSphHarmM.im + x2 * MultSphHarmM.re - x3 * MultSphHarmM.im;
      sigc->data[i] -= x0 * MultSphHarmP.im + x1 * MultSphHarmP.re + x2 * MultSphHarmM.im + x3 * MultSphHarmM.re;
    }
  }

  errcode  = XLALSphHarm(&MultSphHarmP, 4, 0, inc, 0.);
  if (errcode != XLAL_SUCCESS) {
    XLALDestroyREAL8Vector(h40);
    XLALPrintWarning("** LALPSpinInspiralRD WARNING **: impossible to create Y40\n");
  } else {
    for (i = 0; i < length; i++) {
      x0 = h40->data[2 * i];
      x1 = h40->data[2 * i + 1];
      sigp->data[i] += x0 * MultSphHarmP.re - x1 * MultSphHarmM.im;
      sigc->data[i] -= x0 * MultSphHarmP.im + x1 * MultSphHarmP.re;
    }
  }

  params->fFinal = params->tSampling / 2.;
  if (LNhxy) params->polarisationAngle = 0.;
  else       params->polarisationAngle = mparams.polarization;

  /*------------------------------------------------------
   * If required by the user copy other data sets to the
   * relevant arrays
   ------------------------------------------------------*/

  if (hh) {
    for (i = 0; i < length; i++) {
      j = 2 * i;
      k = 2 * i + 1;
      hap->data[j] = sigp->data[i];
      hap->data[k] = sigc->data[i];
    }
  }

  if (signalvec1)
    memcpy(signalvec1->data, sigp->data, length * (sizeof(REAL8)));
  if (signalvec2)
    memcpy(signalvec2->data, sigc->data, length * (sizeof(REAL8)));
  if (hh)
    memcpy(hh->data,         hap->data, 2 * length * (sizeof(REAL8)));
  if (ff)
    memcpy(ff->data,         fap->data, length * (sizeof(REAL8)));
  if (phi)
    memcpy(phi->data,        phap->data, length * (sizeof(REAL8)));

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
  XLALDestroyREAL8Vector(hap);
  XLALDestroyREAL8Vector(sigp);
  XLALDestroyREAL8Vector(sigc);
  return(countback);
  /*End */

}

