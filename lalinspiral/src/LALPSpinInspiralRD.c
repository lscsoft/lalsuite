
/*
*  Copyright (C) 2011 Riccardo Sturani
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

/**
\file
\ingroup LALInspiral_h

 * \brief Module to generate generic spinning binaries waveforms complete with ring-down
 *
 * \heading{Prototypes}
 *
 * <tt>LALPSpinInspiralRD()</tt>
 * <dl>
 * <dt>status:</dt><dd> Input/Output
 * <dt>signalvec:</dt><dd> Output containing the inspiral waveform.</dd>
 * <dt>params:</dt><dd> Input containing binary chirp parameters.</dd>
 * </dl>
 *
 *
 * <tt>LALPSpinInspiralRDTemplates()</tt>
 * <dl>
 * <dt>status:</dt><dd>Input/Output
 * <dt>signalvec1:</dt><dd>Output containing the \f$+\f$ inspiral waveform.</dd>
 * <dt>signalvec2:</dt><dd>Output containing the \f$\times\f$ inspiral waveform.</dd>
 * <dt>params:</dt><dd>Input containing binary chirp parameters.</dd>
 * </dl>
 *
 *
 * <tt>LALPSpinInspiralRDInjection()</tt>
 * <dl>
 * <dt>status:</dt><dd> Input/Output
 * <dt>signalvec:</dt><dd>Output containing the inspiral waveform.</dd>
 * <dt>params:</dt><dd>Input containing binary chirp parameters.</dd>
 * </dl>
 *
 * <tt>LALPSpinInspiralRDFreqDom()</tt>
 * <dl>
 * <dt>status:</dt><dd> Input/Output
 * <dt>signalvec:</dt><dd>Output containing the inspiral waveform in frequency domain.</dd>
 * <dt>params:</dt><dd>Input containing binary chirp parameters.</dd>
 * </dl>
 *
 * \heading{Description}
 * This codes provide complete waveforms for generically spinning binary systems.
 * In order to construct the waveforms three phases are joined together:
 * an initial inspiral phase, a phenomenological phase encompassing the description
 * of the merger and the ring-down of the final black hole.
 * During the inspiral phase the system is evolved according to the standard
 * PN formulas, valid up to 3.5PN for the orbital motion,
 * to 3PN level for spin-orbital momentum and to 2PN for spin-spin contributions
 * to the orbital phase.
 * Then a phenomenological phase is added during which the frequency of the
 * waveform has a pole-like behaviour. The stitching is performed in order to
 * ensure continuity of the phase, the frequency and its first and second
 * derivatives. Finally a ring-down phase is attached.
 *
 * \heading{Algorithm}
 *
 * \heading{Uses}
 * \code
 * LALPSpinInspiralRDderivatives()
 * LALInspiralSetup()
 * LALInspiralChooseModel()
 * LALRungeKutta4()
 * OmMatch()
 * fracRD()
 * XLALPSpinInspiralRDSetParams()
 * XLALSpinInspiralTest()
 * XLALSpinInspiralDerivatives()
 * LALSpinInspiralDerivatives()
 * XLALSpinInspiralFillH2Modes()
 * XLALSpinInspiralFillH3Modes()
 * XLALSpinInspiralFillH4Modes()
 * \endcode
 *
 * \heading{Notes}
 *
*/

/** \defgroup psird Complete phenomenological spin-inspiral waveforms
 * \ingroup ch_inspiral
 *
 * This code provides complete waveforms for generically spinning binary
 * systems.
 *
 */

#include <lal/LALPSpinInspiralRD.h>
#include <lal/LALAdaptiveRungeKutta4.h>

#include <lal/Units.h>
#include <lal/LALInspiral.h>
#include <lal/SeqFactories.h>
#include <lal/NRWaveInject.h>
#include <lal/RealFFT.h>

#define sqrtOnePointFive 1.22474
#define sqrtPoint15      0.387298
#define sqrtFiveOver2    1.1183
#define minIntLen        8

#define UNUSED(expr) do { (void)(expr); } while (0)

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
  REAL8 wdotorb[8];           ///< Coefficients of the analytic PN expansion of \f$\dot\omega_{orb}\f$
  REAL8 wdotorblog;           ///< Log coefficient of the PN expansion of of \f$\dot\omega_{orb}\f$
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

static int XLALPSpinInspiralRDSetParams(LALPSpinInspiralRDparams *mparams,InspiralTemplate *params, InspiralInit *paramsInit) {

  mparams->inspiralOnly = params->inspiralOnly;
  mparams->dt           = 1./params->tSampling;
  mparams->OmCutoff     = params->fCutoff*params->totalMass * LAL_MTSUN_SI * (REAL8) LAL_PI;
  mparams->lengths      = (5.0 / 256.0) / LAL_PI * pow(LAL_PI * params->chirpMass * LAL_MTSUN_SI * params->fLower,-5.0 / 3.0) / params->fLower;
  mparams->omOffset     = 0.006;

  /* setup coefficients for PN equations */
  mparams->m     = params->totalMass;
  mparams->m2m1  = params->mass2 / params->mass1;
  mparams->m1m2  = params->mass1 / params->mass2;
  mparams->m1m   = params->mass1 / params->totalMass;
  mparams->m2m   = params->mass2 / params->totalMass;
  mparams->m1msq = mparams->m1m * mparams->m1m;
  mparams->m2msq = mparams->m2m * mparams->m2m;
  mparams->dm    = (params->mass1 - params->mass2) / params->totalMass;
  
  /* params->eta might have been set up before but just for safety, we
   * recompute it here below.*/
  params->eta = (params->mass1 * params->mass2) / (params->mass1 + params->mass2) / (params->mass1 + params->mass2);
  mparams->eta = params->eta;

  switch (params->order) {

    case LAL_PNORDER_THREE_POINT_FIVE:
      mparams->wdotorb[7] = paramsInit->ak.ST[8];

    case LAL_PNORDER_THREE:
      mparams->epnorb[3] = paramsInit->ak.ETa3;
      mparams->wdotorb[6] = paramsInit->ak.ST[6];
      mparams->wdotorblog = paramsInit->ak.ST[7];
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
      mparams->epnorb[2] = paramsInit->ak.ETa2;
      mparams->wdotorb[4] = paramsInit->ak.ST[4];
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
      mparams->wdotorb[3] = paramsInit->ak.ST[3];
      mparams->epnspin15S1dotLNh = 8. / 3. + 2. * mparams->m2m1;
      mparams->epnspin15S2dotLNh = 8. / 3. + 2. * mparams->m1m2;
      mparams->wdotspin15S1LNh = -(113.0 + 75.0 * mparams->m2m1) / 12.0;
      mparams->wdotspin15S2LNh = -(113.0 + 75.0 * mparams->m1m2) / 12.0;
      mparams->LNhdot15 = 0.5;
      mparams->S1dot15 = (4.0 + 3.0 * mparams->m2m1) / 2.0 * mparams->eta;
      mparams->S2dot15 = (4.0 + 3.0 * mparams->m1m2) / 2.0 * mparams->eta;

    case LAL_PNORDER_ONE:
      mparams->epnorb[1] = paramsInit->ak.ETa1;
      mparams->wdotorb[2] = paramsInit->ak.ST[2];

    case LAL_PNORDER_HALF:
      mparams->wdotorb[1] = paramsInit->ak.ST[1];

    case LAL_PNORDER_NEWTONIAN:
      mparams->epnorb[0] = paramsInit->ak.ETaN;
      mparams->wdotorb[0] = paramsInit->ak.ST[0];
      break;

    case LAL_PNORDER_PSEUDO_FOUR:
      XLALPrintError("*** LALPhenSpinInspiralRD ERROR: PhenSpin approximant not available at pseudo4PN order\n");
      break;

    case LAL_PNORDER_NUM_ORDER:
      XLALPrintError("*** LALPhenSpinInspiralRD ERROR: NUM_ORDER not a valid PN order\n");

    default:
      XLALPrintError("*** LALPhenSpinInspiralRD ERROR: Impossible to create waveform with %d order\n",params->order);
      break;
  }

  switch (params->interaction) {

    case LAL_SIM_INSPIRAL_INTERACTION_NONE:
      /*This kills all spin effects in the phase. Still there are spin effects
	in the waveform due to orbital plane precession*/
      mparams->epnspin15S1dotLNh = 0.;
      mparams->epnspin15S2dotLNh = 0.;
      mparams->wdotspin15S1LNh   = 0.;
      mparams->wdotspin15S2LNh   = 0.;
      mparams->S1dot15           = 0.;
      mparams->S2dot15           = 0.;

    case LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_15PN:  
      /* This keeps only the leading spin-orbit interactions*/
      mparams->wdotspin20S1S2      = 0.;
      mparams->epnspin20S1S2       = 0.;
      mparams->epnspin20S1S2dotLNh = 0.;

    case LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_2PN:
      /* This keeps S1-S2 interactions but kill spin self-interactions*/
      mparams->wdotspin20S1S1 = 0.;
      mparams->epnspin20S1S1 = 0.;
      mparams->epnspin20S2S2 = 0.;
      mparams->Sdot20S = 0.;
      mparams->epnspin20S1S1 = 0.;
      mparams->epnspin20S2S2 = 0.;
      mparams->epnspin20S1S1dotLNh = 0.;
      mparams->epnspin20S2S2dotLNh = 0.;

    case LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_SELF_2PN: 
      /* This kills all spin interaction intervening at 2.5PN order or higher*/
      mparams->epnspin25S1dotLNh   = 0.;
      mparams->epnspin25S2dotLNh   = 0.;
      mparams->wdotspin25S1LNh     = 0.;
      mparams->wdotspin25S2LNh     = 0.;
      mparams->S1dot25             = 0.;
      mparams->S2dot25             = 0.;

    case LAL_SIM_INSPIRAL_INTERACTION_QUAD_MONO_2PN:

    case LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_25PN:
      mparams->wdotspin30S1LNh     = 0.;
      mparams->wdotspin30S2LNh     = 0.;

    case LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_3PN:

    case LAL_SIM_INSPIRAL_INTERACTION_ALL_SPIN:

    case LAL_SIM_INSPIRAL_INTERACTION_TIDAL_5PN:

    case LAL_SIM_INSPIRAL_INTERACTION_TIDAL_6PN:

    case LAL_SIM_INSPIRAL_INTERACTION_ALL: 

    default:
      break;
  }

  return XLAL_SUCCESS;
}

/*
 *
 * function to set derivatives: values and mparams input, dvalues output
 *
 */

static int XLALSpinInspiralTest(double t, const double values[], double dvalues[], void *mparams) {

  LALPSpinInspiralRDparams *params = (LALPSpinInspiralRDparams *) mparams;
  REAL8 omega;
  REAL8 energy;
  REAL8 denergy;
  INT4 returnint;

  UNUSED(t);
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
    returnint= LALPSIRDPN_TEST_ENERGY;
  }
  else if (omega < 0.0) {
    XLALPrintWarning("** LALPSpinInspiralRD WARNING **: Omega has become -ve, this should lead to nan's \n");
    returnint= LALPSIRDPN_TEST_OMEGANONPOS;
  }
  else if (dvalues[1] < 0.0) {
    /* omegadot < 0 */
    returnint= LALPSIRDPN_TEST_OMEGADOT;
  }
  else if (isnan(omega)) {
    /* omega is nan */
    returnint= LALPSIRDPN_TEST_OMEGANAN;
  } 
  else if ((params->inspiralOnly==1)&&(omega>params->OmCutoff)) {
    returnint= LALPSIRDPN_TEST_OMEGACUT;
  }
  else if ((params->inspiralOnly!=1)&&(omega>omegaMatch)) {
    returnint= LALPSIRDPN_TEST_OMEGAMATCH;
  }
  else
    returnint= GSL_SUCCESS;

  return returnint;

}

/**
 * \ingroup psird
 * \brief Module to compute detivative of dynamical variables
 */

static int XLALSpinInspiralDerivatives(double t, const double values[], double dvalues[], void *mparams) {

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

    LALPSpinInspiralRDparams *params= (LALPSpinInspiralRDparams *) mparams;
    
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


void LALSpinInspiralDerivatives(REAL8Vector * values, REAL8Vector * dvalues, void *mparams)
{
  XLALSpinInspiralDerivatives(0.,values->data,dvalues->data,mparams);
}				/* end of LALSpinInspiralDerivatives */

/**
 * \ingroup psird
 * \brief Main module to produce waveforms
 */

static int XLALPSpinInspiralRDEngine(
			REAL8Vector * signalvec1,
			REAL8Vector * signalvec2,
			REAL8Vector * hh,
			REAL8Vector * ff,
			REAL8Vector * phi,
			InspiralTemplate *params,
			InspiralInit     *paramsInit);


void LALPSpinInspiralRD(LALStatus * status, REAL4Vector * signalvec, InspiralTemplate * params)
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  if (XLALPSpinInspiralRD(signalvec, params))
    ABORTXLAL(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}


int XLALPSpinInspiralRD(REAL4Vector * signalvec, InspiralTemplate * params)
{
  INT4 count=0;
  InspiralInit paramsInit;

  if(!signalvec) XLAL_ERROR(XLAL_EFAULT);
  if(!signalvec->data) XLAL_ERROR(XLAL_EFAULT);
  if(!params) XLAL_ERROR(XLAL_EFAULT);
  if(params->nStartPad<0) XLAL_ERROR(XLAL_EBADLEN);
  if(params->nEndPad<0) XLAL_ERROR(XLAL_EBADLEN);
  if(params->fLower <=0) XLAL_ERROR(XLAL_EFREQ);
  if(params->tSampling<=0) XLAL_ERROR(XLAL_ETIME);
  if(params->totalMass<=0) XLAL_ERROR(XLAL_EDOM);

  if (XLALInspiralSetup(&(paramsInit.ak), params))
    XLAL_ERROR(XLAL_EFUNC);
  if (XLALInspiralChooseModel(&(paramsInit.func),&(paramsInit.ak), params))
    XLAL_ERROR(XLAL_EFUNC);

  REAL8Vector *s=XLALCreateREAL8Vector(signalvec->length);
  memset(s->data, 0, s->length * sizeof(REAL8));

  /* Call the engine function */
  count = XLALPSpinInspiralRDEngine(s, NULL, NULL, NULL, NULL, params, &paramsInit);
  if (count == XLAL_FAILURE)
    XLAL_ERROR(XLAL_EFUNC);

  UINT4 i;
  for (i=0;i<s->length;i++) {
    signalvec->data[i]=(REAL4) s->data[i];
  }

  return XLAL_SUCCESS;
}

/**
 * \ingroup psird
 * \brief Module to produce waveform templates
 */

void LALPSpinInspiralRDTemplates(LALStatus * status,
         REAL4Vector * signalvec1,
         REAL4Vector * signalvec2,
         InspiralTemplate * params)
{
    INITSTATUS(status);
    ATTATCHSTATUSPTR(status);

    if (XLALPSpinInspiralRDTemplates(signalvec1, signalvec2, params))
      ABORTXLAL(status);

    DETATCHSTATUSPTR(status);
    RETURN(status);
}

int XLALPSpinInspiralRDTemplates(
         REAL4Vector * signalvec1,
         REAL4Vector * signalvec2,
         InspiralTemplate * params)
{
    INT4 count=0;
    InspiralInit paramsInit;

    if(!signalvec1) XLAL_ERROR(XLAL_EFAULT);
    if(!signalvec1->data) XLAL_ERROR(XLAL_EFAULT);
    if(!signalvec2) XLAL_ERROR(XLAL_EFAULT);
    if(!signalvec2->data) XLAL_ERROR(XLAL_EFAULT);
    if(!params) XLAL_ERROR(XLAL_EFAULT);
    if(params->nStartPad<0) XLAL_ERROR(XLAL_EBADLEN);
    if(params->nEndPad<0) XLAL_ERROR(XLAL_EBADLEN);
    if(params->fLower <=0) XLAL_ERROR(XLAL_EFREQ);
    if(params->tSampling<=0) XLAL_ERROR(XLAL_ETIME);
    if(params->totalMass<=0) XLAL_ERROR(XLAL_EDOM);

    if (XLALInspiralSetup(&(paramsInit.ak), params))
      XLAL_ERROR(XLAL_EFUNC);
    if (XLALInspiralChooseModel(&(paramsInit.func), &(paramsInit.ak), params))
      XLAL_ERROR(XLAL_EFUNC);

    REAL8Vector* s1=XLALCreateREAL8Vector(signalvec1->length);
    REAL8Vector* s2=XLALCreateREAL8Vector(signalvec2->length);

    memset(s1->data, 0, signalvec1->length * sizeof(REAL8));
    memset(s2->data, 0, signalvec2->length * sizeof(REAL8));

    /* Call the engine function*/
    count = XLALPSpinInspiralRDEngine(s1, s2, NULL, NULL, NULL, params, &paramsInit);
    if (count == XLAL_FAILURE)
      XLAL_ERROR(XLAL_EFUNC);

    UINT4 i;
    for (i=0;i<s1->length;i++) {
      signalvec1->data[i]=(REAL4) s1->data[i];
    }
    for (i=0;i<s2->length;i++) {
      signalvec2->data[i]=(REAL4) s2->data[i];
    }

    XLALDestroyREAL8Vector(s1);
    XLALDestroyREAL8Vector(s2);

    return XLAL_SUCCESS;
}

/**
 * \ingroup psird
 * \brief Module to produce injection waveforms
 */

void LALPSpinInspiralRDForInjection(LALStatus        * status,
            CoherentGW       * waveform,
            InspiralTemplate * params,
            PPNParamStruc    * ppnParams)
{
    INITSTATUS(status);
    ATTATCHSTATUSPTR(status);

    if (XLALPSpinInspiralRDForInjection(waveform, params, ppnParams))
      ABORTXLAL(status);

    DETATCHSTATUSPTR(status);
    RETURN(status);
}

int XLALPSpinInspiralRDForInjection(
            CoherentGW       * waveform,
            InspiralTemplate * params,
            PPNParamStruc    * ppnParams)
{
    REAL8Vector *h = NULL;	/* pointers to generated amplitude  data */
    REAL8Vector *f = NULL;	/* pointers to generated  frequency data */
    REAL8Vector *phi = NULL;	/* pointer to generated phase data */

    InspiralInit paramsInit;
    INT4 nbins,count;
    UINT4 i;

    /* Make sure parameter and waveform structures exist. */
    if(!params) XLAL_ERROR(XLAL_EFAULT);
    if(!waveform) XLAL_ERROR(XLAL_EFAULT);
    /* Make sure waveform fields don't exist. */
    if(waveform->a) XLAL_ERROR(XLAL_EFAULT);
    if(waveform->f) XLAL_ERROR(XLAL_EFAULT);
    if(waveform->phi) XLAL_ERROR(XLAL_EFAULT);
    if(waveform->h) XLAL_ERROR(XLAL_EFAULT);

    /* Compute some parameters */
    XLALInspiralInit(params, &paramsInit);
    if (paramsInit.nbins == 0) {
      return XLAL_SUCCESS;
    }

    nbins = paramsInit.nbins;

    /* Now we can allocate memory and vector for coherentGW structure */
    f = XLALCreateREAL8Vector(nbins);
    h = XLALCreateREAL8Vector(2 * nbins);
    phi = XLALCreateREAL8Vector(nbins);

    /* By default the waveform is empty */
    memset(f->data,   0,     nbins * sizeof(REAL8));
    memset(h->data,   0, 2 * nbins * sizeof(REAL8));
    memset(phi->data, 0,     nbins * sizeof(REAL8));

    /* Call the engine function */
    count = XLALPSpinInspiralRDEngine(NULL, NULL, h, f, phi, params, &paramsInit);
    if (count == XLAL_FAILURE)
      XLAL_ERROR(XLAL_EFUNC);

    /* Check an empty waveform hasn't been returned */
    for (i = 0; i < phi->length; i++) {
      if (phi->data[i] != 0.0)
	break;
      if (i == phi->length - 1) {
	XLALDestroyREAL8Vector(f);
	XLALDestroyREAL8Vector(h);
	XLALDestroyREAL8Vector(phi);
      }
    }

    /* Allocate the waveform structures. */
    if ((waveform->h = (REAL4TimeVectorSeries *)
         LALMalloc(sizeof(REAL4TimeVectorSeries))) == NULL) {
      XLAL_ERROR(XLAL_ENOMEM);
    }
    memset(waveform->h, 0, sizeof(REAL4TimeVectorSeries));

    if ((waveform->f = (REAL4TimeSeries *)
         LALMalloc(sizeof(REAL4TimeSeries))) == NULL) {
      LALFree(waveform->h);
      waveform->a = NULL;
      XLAL_ERROR(XLAL_ENOMEM);
    }
    memset(waveform->f, 0, sizeof(REAL4TimeSeries));

    if ((waveform->phi = (REAL8TimeSeries *)
         LALMalloc(sizeof(REAL8TimeSeries))) == NULL) {
      LALFree(waveform->h);
      waveform->h = NULL;
      LALFree(waveform->f);
      waveform->f = NULL;
      XLAL_ERROR(XLAL_ENOMEM);
    }
    memset(waveform->phi, 0, sizeof(REAL8TimeSeries));

    waveform->h->data = XLALCreateREAL4VectorSequence(count, 2);
    waveform->f->data = XLALCreateREAL4Vector(count);
    waveform->phi->data = XLALCreateREAL8Vector(count);

    for (i=0; i<(UINT4)count; i++) {
      waveform->f->data->data[i]     =(REAL4) f->data[i];
      waveform->h->data->data[2*i]   =(REAL4) h->data[2*i];
      waveform->h->data->data[2*i+1] =(REAL4) h->data[2*i+1];
      waveform->phi->data->data[i]   = phi->data[i];
    }

    waveform->h->deltaT = waveform->f->deltaT = waveform->phi->deltaT = 1. / params->tSampling;

    waveform->h->sampleUnits = lalStrainUnit;
    waveform->f->sampleUnits = lalHertzUnit;
    waveform->phi->sampleUnits = lalDimensionlessUnit;

    waveform->position = ppnParams->position;
    waveform->psi = ppnParams->psi;

    snprintf(waveform->h->name, LALNameLength,  "PSpinInspiralRD amplitudes");
    snprintf(waveform->f->name, LALNameLength,  "PSpinInspiralRD frequency");
    snprintf(waveform->phi->name, LALNameLength,"PSpinInspiralRD main phase");

    /* --- fill some output --- */
    ppnParams->tc = params->tC;
    ppnParams->length = count;
    ppnParams->dfdt = ((REAL4) (waveform->f->data->data[count - 1]- waveform->f->data->data[count - 2])) * ppnParams->deltaT;
    ppnParams->fStop = params->fFinal;
    ppnParams->termCode = GENERATEPPNINSPIRALH_EFSTOP;
    ppnParams->termDescription = GENERATEPPNINSPIRALH_MSGEFSTOP;

    ppnParams->fStart = ppnParams->fStartIn;

    /* --- free memory --- */

    XLALDestroyREAL8Vector(f);
    XLALDestroyREAL8Vector(h);
    XLALDestroyREAL8Vector(phi);

    return XLAL_SUCCESS;
} /* End LALPSpinInspiralRDForInjection */

void LALPSpinInspiralRDFreqDom(LALStatus * status,
			       REAL4Vector * signalvec,
			       InspiralTemplate * params)
{
    INITSTATUS(status);
    ATTATCHSTATUSPTR(status);

    if (XLALPSpinInspiralRDFreqDom(signalvec, params))
      ABORTXLAL(status);

    DETATCHSTATUSPTR(status);
    RETURN(status);
}

int XLALPSpinInspiralRDFreqDom(
			       REAL4Vector * signalvec,
			       InspiralTemplate * params)
{
    REAL8Vector *tsignalvec = NULL;
    REAL4Vector *fsignalvec = NULL;
    REAL4FFTPlan *forwPlan = NULL;

    InspiralInit paramsInit;

    INT4 count;
    UINT4 nbins,idx;
    UINT4 length = signalvec->length;

    if(!signalvec) XLAL_ERROR(XLAL_EFAULT);
    if(!signalvec->data) XLAL_ERROR(XLAL_EFAULT);
    if(!params) XLAL_ERROR(XLAL_EFAULT);
    if(params->nStartPad<0) XLAL_ERROR(XLAL_EBADLEN);
    if(params->nEndPad<0) XLAL_ERROR(XLAL_EBADLEN);
    if(params->fLower <=0) XLAL_ERROR(XLAL_EFREQ);
    if(params->tSampling<=0) XLAL_ERROR(XLAL_ETIME);
    if(params->totalMass<=0) XLAL_ERROR(XLAL_EDOM);

    XLALInspiralInit(params, &paramsInit);

    if (paramsInit.nbins == 0) {
      return XLAL_SUCCESS;
    }
    nbins = paramsInit.nbins;
    if (nbins < length) nbins = length;
    tsignalvec = XLALCreateREAL8Vector(nbins);
    fsignalvec = XLALCreateREAL4Vector(nbins);

    memset(signalvec->data, 0, length * sizeof(REAL4));
    memset(tsignalvec->data, 0, nbins * sizeof(REAL8));
    memset(fsignalvec->data, 0, nbins * sizeof(REAL4));

    /* Call the engine function */
    count = XLALPSpinInspiralRDEngine(tsignalvec, NULL, NULL, NULL, NULL, params, &paramsInit);
    if (count == XLAL_FAILURE)
      XLAL_ERROR(XLAL_EFUNC);

    REAL4Vector *tsigR4=XLALCreateREAL4Vector(nbins);
    for (idx=0;idx<nbins;idx++) tsigR4->data[idx]=(REAL4) tsignalvec->data[idx];
    XLALDestroyREAL8Vector(tsignalvec);
    XLALSimInspiralREAL4WaveTaper(tsigR4, (LALSimInspiralApplyTaper) 3);

    forwPlan = XLALCreateForwardREAL4FFTPlan(nbins, 0);
    if (forwPlan == NULL) {
      XLALDestroyREAL4Vector(fsignalvec);
      XLALDestroyREAL8Vector(tsignalvec);
      XLALDestroyREAL4Vector(tsigR4);
      XLALDestroyREAL4FFTPlan(forwPlan);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    XLALREAL4VectorFFT(fsignalvec, tsigR4, forwPlan);
    XLALDestroyREAL4Vector(tsigR4);

    for (idx = 0; idx < nbins; idx++) fsignalvec->data[idx] /= params->tSampling;

    if (nbins > length) {
      /*do interpolation*/
      REAL8Vector *fsigRe  = XLALCreateREAL8Vector(nbins/2);
      REAL8Vector *fsigIm  = XLALCreateREAL8Vector(nbins/2);
      REAL8Vector *freq    = XLALCreateREAL8Vector(length/2);
      REAL8Vector *freqSup = XLALCreateREAL8Vector(nbins/2);

      REAL8 dF = 1./params->tSampling/((REAL8) length);
      REAL8 dFsup = 1./params->tSampling/((REAL8) nbins);

      for (idx=0; idx< length/2; idx++) freq->data[idx]=((REAL4) idx) *dF;

      for (idx = 0; idx < nbins/2; idx++) {
	fsigRe->data[idx]=(REAL8)fsignalvec->data[idx];
	fsigIm->data[idx]=(REAL8)fsignalvec->data[nbins-idx];
	freqSup->data[idx]=((REAL8) idx) *dFsup;
      }

      gsl_interp_accel* acc;
      gsl_spline* spline_real;
      gsl_spline* spline_imag;
      XLAL_CALLGSL( spline_imag  = (gsl_spline*) gsl_spline_alloc(gsl_interp_cspline, nbins/2));
      XLAL_CALLGSL( spline_real  = (gsl_spline*) gsl_spline_alloc(gsl_interp_cspline, nbins/2));
      XLAL_CALLGSL( acc          = (gsl_interp_accel*) gsl_interp_accel_alloc());
      XLAL_CALLGSL( gsl_spline_init(spline_real, freqSup->data, fsigRe->data, nbins/2));
      XLAL_CALLGSL( gsl_spline_init(spline_imag, freqSup->data, fsigIm->data, nbins/2));

      for (idx=0;idx<length/2;idx++){
        signalvec->data[idx]        = gsl_spline_eval ( spline_real , freq->data[idx] , acc);
        signalvec->data[length-idx] = gsl_spline_eval ( spline_imag , freq->data[idx] , acc);
      }
      signalvec->data[0] = 0.;
      signalvec->data[length / 2] = 0.;

      XLALDestroyREAL8Vector(fsigRe);
      XLALDestroyREAL8Vector(fsigIm);
      XLALDestroyREAL8Vector(freq);
      XLALDestroyREAL8Vector(freqSup);
      gsl_spline_free(spline_real);
      gsl_spline_free(spline_imag);
      gsl_interp_accel_free(acc);

    }
    else {
      for (idx=0; idx<length; idx++)
	signalvec->data[idx]=fsignalvec->data[idx];
    }

    XLALDestroyREAL4Vector(fsignalvec);
    XLALDestroyREAL4FFTPlan(forwPlan);

    return XLAL_SUCCESS;
}

/**
 * \ingroup psird
 * \brief Module actually computing PSIRD waveforms
 */

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
				){

  REAL8 amp20 = amp * sqrtOnePointFive;
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

  return 0;

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
				){

  REAL8 amp32 = amp * sqrtOnePointFive;
  REAL8 amp31 = amp * sqrtPoint15;
  REAL8 amp30 = amp / sqrtFiveOver2;
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

    return 0;

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
				REAL8 v,
				REAL8 eta,
				REAL8 dm,
				REAL8 Psi,
				REAL8 alpha,
				LALSpinInspiralAngle *an
				){

  UNUSED(v);
  UNUSED(dm);

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

  return 0;
}

static int XLALSpinInspiralEngine(
				UINT4 neqs, 
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
  INT4 intreturn;
  UINT4 count=0;
  UINT4 write=0;
  UINT4 modcount=0;
  UINT4 j;

  REAL8 dt,tm,timewrite;
  REAL8 Mass,dm;
  REAL8 v,v2;
  REAL8 Phi=0.;
  REAL8 omega;

  REAL8 LNhx,LNhy,LNhz;
  REAL8 S1x,S1y,S1z;
  REAL8 S2x,S2y,S2z;
  REAL8 LNhxy;
  REAL8 energy      = 0.;
  REAL8 energywrite = 0.;

  REAL8 LNhS1=0.;
  REAL8 LNhS2=0.;
  REAL8 S1S1=0.;
  REAL8 S1S2=0.;
  REAL8 S2S2=0.;
  REAL8 dLNhx,dLNhy,dLNhz;
  REAL8 LNhS1w = 0.;
  REAL8 LNhS2w = 0.;
  REAL8 S1S1w  = 0.;
  REAL8 S2S2w  = 0.;
  REAL8 S1S2w  = 0.;

  //REAL8 Phiold;
  REAL8 alpha,alphaold;

  REAL8 Phiwrite     = 0.;
  REAL8 alphawrite   = 0.;
  REAL8 omegawrite   = 0.;

  REAL8 amp22,amp33,amp44;
  REAL8 unitHz;
  LALSpinInspiralAngle trigAngle;
  UINT4 subsampling=1;

  UINT4 Npoints = 20;
  INT4 errcode;

  rk4In in4;      // used to setup the Runge-Kutta integration
  rk4GSLIntegrator *integrator;

  REAL8Vector dummy, values, dvalues, newvalues, yt, dym, dyt;
  // support variables

  dummy.length = neqs * 6;

  values.length = dvalues.length = newvalues.length = yt.length = dym.length = dyt.length = neqs;
  
  if (!(dummy.data = (REAL8 *) LALMalloc(sizeof(REAL8) * neqs * 6))) {
    XLAL_ERROR(XLAL_ENOMEM);
  }

  dt=mparams->dt;
  Mass=mparams->m * LAL_MTSUN_SI;
  dm=mparams->dm;
  unitHz= Mass * (REAL8) LAL_PI;

  values.data    = &dummy.data[0];
  dvalues.data   = &dummy.data[neqs];
  newvalues.data = &dummy.data[2 * neqs];
  yt.data        = &dummy.data[3 * neqs];
  dym.data       = &dummy.data[4 * neqs];
  dyt.data       = &dummy.data[5 * neqs];

  REAL8 S1x0=values.data[5];
  REAL8 S1y0=values.data[6];
  REAL8 S1z0=values.data[7];
  REAL8 S2x0=values.data[8];
  REAL8 S2y0=values.data[9];
  REAL8 S2z0=values.data[10];

  /* Variables initializations */
  for (j=0;j<neqs;j++) values.data[j]=yinit[j];

  omega  = values.data[1];
  //Phiold = Phi;
  Phi    = values.data[0];
  v      = cbrt(omega);
  v2     = v*v;
  alpha  = atan2(values.data[3],values.data[2]);
  trigAngle.ci = LNhz = values.data[4];

  LNhS1 = (values.data[2]*values.data[5]+values.data[3]*values.data[6]+values.data[4]*values.data[7])/mparams->m1msq;
  LNhS2 = (values.data[2]*values.data[8]+values.data[3]*values.data[9]+values.data[4]*values.data[10])/mparams->m2msq;
  S1S1  = (values.data[5]*values.data[5]+values.data[6]*values.data[6]+values.data[7]*values.data[7])/mparams->m1msq/mparams->m1msq;
  S2S2  = (values.data[8]*values.data[8]+values.data[9]*values.data[9]+values.data[10]*values.data[10])/mparams->m2msq/mparams->m2msq;
  S1S2= (values.data[5]*values.data[8]+values.data[6]*values.data[9]+values.data[7]*values.data[10])/mparams->m1msq/mparams->m2msq;;

  while ( (OmMatch(LNhS1,LNhS1,S1S1,S1S2,S2S2) * 16. / unitHz) > (REAL4) (subsampling) / dt ) {
    subsampling *= 2;
    dt /= (REAL8) (subsampling);
  }

  in4.function = LALSpinInspiralDerivatives;
  in4.y        = &values;
  in4.dydx     = &dvalues;
  in4.h        = dt / Mass;
  in4.n        = neqs;
  in4.yt       = &yt;
  in4.dym      = &dym;
  in4.dyt      = &dyt;

  /* Start of the calculation of the inspiral part via the fixed step integration method */

  /* Initialize GSL integrator */
  if (!(integrator = XLALRungeKutta4Init(neqs, &in4))) {
    XLALFree(dummy.data);
    INT4 errNum = XLALClearErrno();
    if (errNum == XLAL_ENOMEM)
      XLAL_ERROR(XLAL_ENOMEM);
    else
    XLAL_ERROR(XLAL_EDOM);
  }

  count = write= 0;
  tm = timewrite = 0.;

  LALSpinInspiralDerivatives(&values, &dvalues, (void *) &mparams);

  omega  = values.data[1];
  //Phiold = Phi;
  Phi    = values.data[0];
  v      = cbrt(omega);
  v2     = v*v;
  alpha  = atan2(values.data[3],values.data[2]);

  REAL8Vector *domega  = XLALCreateREAL8Vector(Npoints);
  REAL8Vector *diota   = XLALCreateREAL8Vector(Npoints);
  REAL8Vector *dalpha  = XLALCreateREAL8Vector(Npoints);
  REAL8Vector *ddomega = XLALCreateREAL8Vector(Npoints);
  REAL8Vector *ddiota  = XLALCreateREAL8Vector(Npoints);
  REAL8Vector *ddalpha = XLALCreateREAL8Vector(Npoints);

  do {

    if (count%subsampling==0) {

      modcount=0;

      if (write >= mparams->length) {
	XLALRungeKutta4Free(integrator);
	XLALFree(dummy.data);
	XLAL_ERROR(XLAL_EDOM);
      }

      amp22 = amp22ini * v2;
      amp33 = -amp22 / 4. * sqrt(5. / 42.);
      amp44 = amp22 * v2 * 2.*sqrt(5./7.)/9.;

      Phiwrite   = Phi;
      omegawrite = omega;
      alphawrite = alpha;
      energywrite  = energy;
      LNhS1w = LNhS1;
      LNhS2w = LNhS2;
      S1S1w  = S1S1;
      S1S2w  = S1S2;
      S2S2w  = S2S2;

      trigAngle.ci   = LNhz;
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

      XLALSpinInspiralFillH2Modes(h2P2,h2M2,h2P1,h2M1,h20,write,amp22,v,mparams->eta,dm,Phiwrite,alphawrite,&trigAngle);

      XLALSpinInspiralFillH3Modes(h3P3,h3M3,h3P2,h3M2,h3P1,h3M1,h30,write,amp33,v,mparams->eta,dm,Phiwrite,alphawrite,&trigAngle);

      XLALSpinInspiralFillH4Modes(h4P4,h4M4,h4P3,h4M3,h4P2,h4M2,h4P1,h4M1,h40,write,amp44,v,mparams->eta,dm,Phiwrite,alphawrite,&trigAngle);

      freq->data[write]=omega;
      phase->data[write]=Phi;

      write++;

      timewrite+=mparams->dt;

    }

    in4.x = tm / mparams->m;

    if (XLALRungeKutta4(&newvalues, integrator,(void *) mparams) == XLAL_FAILURE)
      XLAL_ERROR(XLAL_EFUNC);
    /* updating values of dynamical variables */

    Phi = values.data[0] = newvalues.data[0];
    omega = values.data[1] = newvalues.data[1];

    LNhx = values.data[2] = newvalues.data[2];
    LNhy = values.data[3] = newvalues.data[3];
    LNhz = values.data[4] = newvalues.data[4];

    S1x = values.data[5] = newvalues.data[5];
    S1y = values.data[6] = newvalues.data[6];
    S1z = values.data[7] = newvalues.data[7];

    S2x = values.data[8] = newvalues.data[8];
    S2y = values.data[9] = newvalues.data[9];
    S2z = values.data[10] = newvalues.data[10];

    energy = values.data[11] = newvalues.data[11];

    alphaold = alpha;
    LNhxy = sqrt(LNhx * LNhx + LNhy * LNhy);
    if (LNhxy>0.)
      alpha = atan2(LNhy, LNhx);
    else
      alpha = alphaold;

    /*if (count>1) {
      if ((alpha*alphaold)<0.) {
	if (fabs(cos(2.*(Phi+alpha))-cos(2.*(Phiold+alphaold)))>0.2) {
	fprintf(stdout,"*** LALPSpinInspiralRD WARNING ***: Possible problem with coordinate singularity:\n Step %d  LNhy: %12.6e LNhx: %12.6e  Psi+alpha: %12.6e\n Step %d      Psiold+alphaold %12.6e\n",write,LNhy,LNhx,(Phi+alpha)/LAL_PI,write-1,(Phiold+alphaold)/LAL_PI);
	fprintf(stdout,"            m: (%12.6e,%12.6e)\n", mparams->m1m*mparams->m, mparams->m2m*mparams->m);
	fprintf(stdout,"            S1: (%9.6f,%9.6f,%9.6f)\n",yinit[5]/mparams->m1msq,yinit[6]/mparams->m1msq,yinit[7]/mparams->m1msq);
	fprintf(stdout,"            S2: (%9.6f,%9.6f,%9.6f)\n",yinit[8]/mparams->m2msq,yinit[9]/mparams->m2msq,yinit[10]/mparams->m2msq);
	}
      }
      }*/

    LNhS1 = (S1x * LNhx + S1y * LNhy + S1z * LNhz) / mparams->m1msq;
    LNhS2 = (S2x * LNhx + S2y * LNhy + S2z * LNhz) / mparams->m2msq;
    S1S1  = (S1x * S1x + S1y * S1y + S1z * S1z)/mparams->m1msq/mparams->m1msq;
    S2S2  = (S2x * S2x + S2y * S2y + S2z * S2z)/mparams->m2msq/mparams->m2msq;
    S1S2  = (S1x * S2x + S1y * S2y + S1z * S2z)/mparams->m1msq/mparams->m2msq;

    LALSpinInspiralDerivatives(&values, &dvalues, (void *) mparams);

    dLNhx = dvalues.data[2];
    dLNhy = dvalues.data[3];
    dLNhz = dvalues.data[4];
    if (LNhxy > 0.) {
      diota->data[Npoints-1]  = -dLNhz / LNhxy;
      dalpha->data[Npoints-1] = (LNhx * dLNhy - LNhy * dLNhx) / LNhxy;
    } else {
      diota->data[Npoints-1]  = 0.;
      dalpha->data[Npoints-1] = 0.;
    }

    v  = cbrt(omega);
    v2 = v*v;

    for (j=0;j<Npoints-1;j++) {
      domega->data[j] = domega->data[j+1];
      diota->data[j]  = diota->data[j+1];
      dalpha->data[j] = dalpha->data[j+1];
    }
    domega->data[Npoints-1] = dvalues.data[1];

    tm += dt;

    count++;
    modcount++;

    intreturn=XLALSpinInspiralTest(0.,values.data,dvalues.data,mparams);

  } while (intreturn==GSL_SUCCESS);

  XLALRungeKutta4Free(integrator);
  XLALFree(dummy.data);

  if (count<Npoints) {
    XLALPrintError("*** LALPSpinInspiralRD WARNING: inspiral integration vey short: %12.f sec\n",tm);
    XLAL_ERROR(XLAL_EFAILED);
  }

  errcode = XLALGenerateWaveDerivative(ddomega,domega,dt);
  errcode += XLALGenerateWaveDerivative(ddalpha,dalpha,dt);
  errcode += XLALGenerateWaveDerivative(ddiota,diota,dt);

  if (errcode != 0) {
    XLALPrintError("**** LALPSpinInspiralRD ERROR ****: error generating derivatives\n");
    XLALPrintError("                     m:           : %12.5f  %12.5f\n",mparams->m1m*mparams->m,mparams->m2m*mparams->m);
    XLALPrintError("              S1:                 : %12.5f  %12.5f  %12.5f\n",S1x0,S1y0,S1z0);
    XLALPrintError("              S2:                 : %12.5f  %12.5f  %12.5f\n",S2x0,S2y0,S2z0);
    XLAL_ERROR(XLAL_EFUNC);
  }

  phenPars->endtime   = timewrite-mparams->dt;
  phenPars->intreturn = intreturn;
  phenPars->Psi       = Phiwrite;
  phenPars->alpha     = alphawrite;
  phenPars->energy    = energywrite;
  phenPars->omega     = omegawrite;
  phenPars->domega    = domega->data[Npoints-1-modcount]/Mass;
  phenPars->ddomega   = ddomega->data[Npoints-1-modcount]/Mass;
  phenPars->diota     = diota->data[Npoints-1-modcount];
  phenPars->ddiota    = ddiota->data[Npoints-1-modcount];
  phenPars->dalpha    = dalpha->data[Npoints-1-modcount] ;
  phenPars->ddalpha   = ddalpha->data[Npoints-1-modcount];
  phenPars->ci        = trigAngle.ci;
  phenPars->countback = write-1;
  phenPars->LNhS1     = LNhS1w;
  phenPars->LNhS2     = LNhS2w;
  phenPars->S1S2      = S1S2w;
  phenPars->S1S1      = S1S1w;
  phenPars->S2S2      = S2S2w;

  XLALDestroyREAL8Vector(domega);
  XLALDestroyREAL8Vector(diota);
  XLALDestroyREAL8Vector(dalpha);
  XLALDestroyREAL8Vector(ddomega);
  XLALDestroyREAL8Vector(ddiota);
  XLALDestroyREAL8Vector(ddalpha);

  return XLAL_SUCCESS;
} /* End of XLALSpinInspiralEngine*/

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

  UINT4 j;
  UINT4 k;
  UINT4 kMatch=0;
  UINT4 jMatch=0;
  UINT4 Npoints=10;
  UINT4 intlen;
  UINT4 intreturn;

  LALSpinInspiralAngle trigAngle;

  REAL8Array *yout=NULL;
  //yout=malloc(sizeof(REAL8Array));
  ark4GSLIntegrator *integrator=NULL;
  //integrator=malloc(sizeof(ark4GSLIntegrator));
 //memset(&integrator,0,sizeof(ark4GSLIntegrator)+1);

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

  INT4 errcode;

  REAL8 *yin = XLALMalloc(sizeof(REAL8) * neqs);
  if (!yin) XLAL_ERROR(XLAL_ENOMEM);

  /* allocate the integrator */
  integrator = XLALAdaptiveRungeKutta4Init(neqs,XLALSpinInspiralDerivatives,XLALSpinInspiralTest,1.0e-6,1.0e-6);
  if (!integrator) {
    XLALPrintError("**** LALPSpinInspiralRD ERROR ****: Cannot allocate adaptive integrator.\n");
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

  for (j=0; j<neqs; j++) yin[j]=yinit[j];

  REAL8 S1x0=yinit[5];
  REAL8 S1y0=yinit[6];
  REAL8 S1z0=yinit[7];
  REAL8 S2x0=yinit[8];
  REAL8 S2y0=yinit[9];
  REAL8 S2z0=yinit[10];

  intlen = XLALAdaptiveRungeKutta4Hermite(integrator,(void *)mparams,yin,0.0,mparams->lengths/Mass,dt/Mass,&yout);

  intreturn = integrator->returncode;
  XLALAdaptiveRungeKutta4Free(integrator);

  /* End integration*/

  /* Start of the integration checks*/
  if (!intlen) {
    if (XLALClearErrno() == XLAL_ENOMEM) {
      XLAL_ERROR(  XLAL_ENOMEM);
    } else {
      XLALPrintError("**** LALPSpinInspiralRD ERROR ****: integration failed with errorcode %d, integration length %d\n",intreturn,intlen);
      XLAL_ERROR( XLAL_EFAILED);
    }
  }

  /* if we have enough space, compute the waveform components; otherwise abort */
  if ( intlen >= mparams->length ) {
    XLALPrintError("**** LALPSpinInspiralRD ERROR ****: no space to write in waveforms: %d vs. %d\n",intlen,mparams->length);
    XLAL_ERROR(XLAL_ESIZE);
  }

  if ( intlen < minIntLen ) {
    XLALPrintError("**** LALPSpinInspiralRD ERROR ****: incorrect integration with length %d\n",intlen);
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

  if (mparams->inspiralOnly!=1) {

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
      XLALPrintError("*** LALPSpinInspiralRD ERROR ***: Impossible to attach phenom. part\n");
      XLAL_ERROR(XLAL_EFAILED);
    }

    // Data structure are copied into Npoints-long
    // REAL8Array for interpolation and derivative computation
    if (Npoints > intlen) Npoints = intlen;

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
    if (errcode != XLAL_SUCCESS) {
      XLALPrintError("**** LALPSpinInspiralRD ERROR ****: error generating first derivatives: #points %d\n",Npoints);      
      XLALPrintError("                     m:           : %12.5f  %12.5f\n",mparams->m1m*mparams->m,mparams->m2m*mparams->m);
      XLALPrintError("              S1:                 : %12.5f  %12.5f  %12.5f\n",S1x0,S1y0,S1z0);
      XLALPrintError("              S2:                 : %12.5f  %12.5f  %12.5f\n",S2x0,S2y0,S2z0);
      XLALPrintError("     omM %12.5f   om[%d] %12.5f\n",omegaMatch,jMatch,omega);
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
    if (errcode != XLAL_SUCCESS) {
      XLALPrintError("**** LALPSpinInspiralRD ERROR ****: error generating second derivatives\n");
      XLALPrintError("                     m:           : %12.5f  %12.5f\n",mparams->m1m*mparams->m,mparams->m2m*mparams->m);
      XLALPrintError("              S1:                 : %12.5f  %12.5f  %12.5f\n",S1x0,S1y0,S1z0);
      XLALPrintError("              S2:                 : %12.5f  %12.5f  %12.5f\n",S2x0,S2y0,S2z0);
      XLALPrintError("     omM %12.5f   om[%d] %12.5f\n",omegaMatch,jMatch,omega);
      XLAL_ERROR(XLAL_EFAILED);
    }

    if (ddomega->data[kMatch]<0.) {
      XLALPrintWarning("*** LALPSpinInspiralRD WARNING: the attach of the phenom. phase has been shifted back: m1 %12.6f  m2 %12.6f\n",mparams->m1m*mparams->m,mparams->m2m*mparams->m);
      XLALPrintWarning("  Integration returned %d\n   1025: Energy increases\n   1026: Omegadot -ve\n   1028: Omega NAN\n   1029: Omega > Omegamatch\n   1031: Omega -ve\n   1032: Omega > OmegaCut %12.6e\n",intreturn,mparams->OmCutoff); 
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

  //REAL8 alphaoold = 0.;
  alphaold=alpha;
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

    errcode  = XLALSpinInspiralFillH2Modes(h2P2,h2M2,h2P1,h2M1,h20,j,amp22,v,mparams->eta,mparams->dm,Psi,alpha,&trigAngle);

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

    errcode += XLALSpinInspiralFillH3Modes(h3P3,h3M3,h3P2,h3M2,h3P1,h3M1,h30,j,amp33,v,mparams->eta,mparams->dm,Psi,alpha,&trigAngle);

    errcode += XLALSpinInspiralFillH4Modes(h4P4,h4M4,h4P3,h4M3,h4P2,h4M2,h4P1,h4M1,h40,j,amp44,v,mparams->eta,mparams->dm,Psi,alpha,&trigAngle);

    if (errcode != XLAL_SUCCESS)
      XLAL_ERROR(XLAL_EFUNC);
  }

  phenPars->alpha=alpha;

  if (yin)  XLALFree(yin);
  if (yout) XLALDestroyREAL8Array(yout);  

  return XLAL_SUCCESS;

} /* End of the inspiral part created via the adaptive integration method */


static int XLALPSpinInspiralRDEngine(
			REAL8Vector * signalvec1,
			REAL8Vector * signalvec2,
			REAL8Vector * hh,
			REAL8Vector * ff,
			REAL8Vector * phi,
			InspiralTemplate *params,
			InspiralInit     *paramsInit)
{

  /* declare code parameters and variables */
  const INT4 neqs = 11+1;      // number of dynamical variables plus the energy function
  UINT4 origcount,apcount;     // integration steps performed
  UINT4 count = 0;             // integration steps performed
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
  REAL8 initS1[3],initS2[3],initLNh[3],initJ[3];
  REAL8 iS1[3],iS2[3];
  REAL8 phiJ,thetaJ;
  REAL8 ry[3][3],rz[3][3];
  REAL8 dt;

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
  REAL8Vector* hap;
  REAL8Vector* phap;

  LALPSpinInspiralRDparams mparams;
  LALPSpinInspiralPhenPars phenPars;
  LALSpinInspiralAngle trigAngle;

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


  if(!params)
    XLAL_ERROR(XLAL_EFAULT);

  if ((params->fCutoff<=0.)&&(params->inspiralOnly==1)) {
    XLALPrintError("*** LALPSIRD ERROR ***: fCutoff %12.6e, with inspiral flag on it is mandatory to specify a positive cutoff frequency\n",params->fCutoff);
    XLAL_ERROR(XLAL_EDOM);
  }

  mass = params->totalMass * LAL_MTSUN_SI;
  unitHz = params->totalMass * LAL_MTSUN_SI * (REAL8) LAL_PI;

  if ((signalvec2)||(hh))
    params->nStartPad = 0;    /* must be zero for templates and injection */
  /* -- length in seconds from Newtonian formula; */

  dt = 1. / params->tSampling;

  /* setup coefficients for PN equations */
  XLALPSpinInspiralRDSetParams(&mparams,params,paramsInit);

  /* Check that initial frequency is smaller than omegamatch ~ xxyy for m=100 Msun */
  initphi   = params->startPhase/2.;
  initomega = params->fLower*unitHz;

  /* Check that initial frequency is smaller than omegamatch ~ xxyy for m=100 Msun */

  LNhS1=params->spin1[2];
  LNhS2=params->spin2[2];
  S1S1=params->spin1[0]*params->spin1[0]+params->spin1[1]*params->spin1[1]+params->spin1[2]*params->spin1[2];
  S1S2=params->spin1[0]*params->spin2[0]+params->spin1[1]*params->spin2[1]+params->spin1[2]*params->spin2[2];
  S2S2=params->spin2[0]*params->spin2[0]+params->spin2[1]*params->spin2[1]+params->spin2[2]*params->spin2[2];

  omegaMatch = OmMatch(LNhS1,LNhS2,S1S1,S1S2,S2S2);

  if ( initomega > omegaMatch ) {
    /*if ((params->spin1[0]==params->spin1[1])&&(params->spin1[1]==params->spin2[0])&&(params->spin2[0]==params->spin2[1])&&(params->spin2[1]==0.)) {
      //Beware, this correspond to a shift of the initial phase!
      initomega = 0.95*omegaMatch;
      fprintf(stdout,"*** LALPSpinInspiralRD WARNING ***: Initial frequency reset from %12.6e to %12.6e Hz, m:(%12.4e,%12.4e)\n",params->fLower,initomega/unitHz,params->mass1,params->mass2);
      }*/
    /*else {*/
    XLALPrintError("**** LALPSpinInspiralRD ERROR ****: Initial frequency too high: %11.5e for omM ~ %11.5e and m:(%8.3f, %8.3f)\n",params->fLower,omegaMatch/unitHz,params->mass1,params->mass2);
    XLAL_ERROR(XLAL_EFAILED);
    /*}*/
  }

  /* Here we use the following convention:
     the coordinates of the spin vectors params->spin1,2 and the params->inclination 
     variable refers to different physical parameters according to the value of 
     params->axisChoice:

     * OrbitalL: params->inclination denotes the angle between the view direction
                 N and the initial L (initial L//z, N in the x-z plane) and the spin 
		 coordinates are given with respect to initial L.
     * TotalJ:   params->inclination denotes the angle between the view directoin 
                 and J (J is constant during the evolution, J//z, both N and initial 
		 L are in the x-z plane) and the spin coordinates are given wrt 
		 initial ** L **.

     * View:     params->inclination denotes the angle between the initial L and N 
                 (N//z, initial L in the x-z plane) and the spin coordinates 
		 are given with respect to N.

     In order to reproduce the results of the SpinTaylor code View must be chosen.
     The spin magnitude are normalized to the individual mass^2, i.e.
     they are dimension-less.
     The modulus of the initial angular momentum is fixed by m1,m2 and
     initial frequency.
     The polarization angle is not used here, it enters the pattern
     functions along with the angles marking the sky position of the
     source. */

  // Physical magnitude of the orbital angular momentum
  LNhmag = params->eta * params->totalMass * params->totalMass / cbrt(initomega);

  // Physical values of the spins
  for (i = 0; i < 3; i++) {
    initS1[i] = params->spin1[i] * params->mass1 * params->mass1;
    initS2[i] = params->spin2[i] * params->mass2 * params->mass2;
  }

  switch (params->axisChoice) {

  case LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L:
    //printf("*** OrbitalL ***\n");
    initLNh[0] = 0.;
    initLNh[1] = 0.;
    initLNh[2] = 1.;
    inc = params->inclination;
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
    inc = params->inclination;
    break;

  case LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW:
  default:
    //printf("*** View ***\n");
    initLNh[0] = sin(params->inclination);
    initLNh[1] = 0.;
    initLNh[2] = cos(params->inclination);
    inc = 0.;
    break;
  }

  /*All the PN formulas used in the differential equation integration 
    assume that the spin variables are the physical ones divided by
    totalmasss^2, here we introduce the correct normalization, changing the
    input one, where spin components were normalized on individual mass. */
  for (j = 0; j < 3; j++) {
    initS1[j] /= params->totalMass * params->totalMass;
    initS2[j] /= params->totalMass * params->totalMass;
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
  memset(hap->data,  0, hap->length  * sizeof(REAL8));
  memset(fap->data,  0, fap->length  * sizeof(REAL8));
  memset(phap->data, 0, phap->length * sizeof(REAL8));

  /* Here there used to be a check that OmegaRD is smaller than Nyquist, it
     has been taken out */

  params->ampOrder = (LALPNOrder) 1;
  if (params->distance > 0.)
    amp22ini = -2.0 * params->mu * LAL_MRSUN_SI / params->distance * sqrt(16. * LAL_PI / 5.);
  else
    amp22ini = 2. * sqrt(LAL_PI / 5.0) * params->signalAmplitude;

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
    if (XLALSpinInspiralEngine(neqs,yinit,amp22ini,&mparams,h2P2,h2M2,h2P1,h2M1,h20,h3P3,h3M3,h3P2,h3M2,h3P1,h3M1,h30,h4P4,h4M4,h4P3,h4M3,h4P2,h4M2,h4P1,h4M1,h40,fap,phap,&phenPars) == XLAL_FAILURE) {
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
      XLALDestroyREAL8Vector(sigp);
      XLALDestroyREAL8Vector(sigc);
      XLALDestroyREAL8Vector(fap);
      XLALDestroyREAL8Vector(hap);
      XLALDestroyREAL8Vector(phap);
      XLAL_ERROR(XLAL_EFUNC);
    }
  }  else {
    if (XLALSpinInspiralAdaptiveEngine(neqs,yinit,amp22ini,&mparams,h2P2,h2M2,h2P1,h2M1,h20,h3P3,h3M3,h3P2,h3M2,h3P1,h3M1,h30,h4P4,h4M4,h4P3,h4M3,h4P2,h4M2,h4P1,h4M1,h40,fap,phap,&phenPars) == XLAL_FAILURE) {
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
      XLALDestroyREAL8Vector(sigp);
      XLALDestroyREAL8Vector(sigc);
      XLALDestroyREAL8Vector(fap);
      XLALDestroyREAL8Vector(hap);
      XLALDestroyREAL8Vector(phap);
      XLAL_ERROR(XLAL_EFUNC);      
    }
  }
  intreturn=phenPars.intreturn;
  /* report on abnormal termination:
     Termination is fine if omegamatch is passed or if energy starts 
     increasing  */

  if ( (intreturn!=LALPSIRDPN_TEST_OMEGACUT) && (intreturn != LALPSIRDPN_TEST_OMEGAMATCH) && (intreturn != LALPSIRDPN_TEST_ENERGY) )
    {
      XLALPrintWarning("** LALPSpinInspiralRD WARNING **: integration terminated with code %d.\n",intreturn);
      XLALPrintWarning("  1025: Energy increases\n  1026: Omegadot -ve\n  1028: Omega NAN\n  1029: Omega > Omegamatch\n  1031: Omega -ve\n  1032: Omega > OmegaCut %12.6e\n",mparams.OmCutoff);
      XLALPrintWarning("  Waveform parameters were m1 = %14.6e, m2 = %14.6e, inc = %10.6f,\n", params->mass1, params->mass2, params->inclination);
      XLALPrintWarning("                           S1 = (%10.6f,%10.6f,%10.6f)\n", params->spin1[0], params->spin1[1], params->spin1[2]);
      XLALPrintWarning("                           S2 = (%10.6f,%10.6f,%10.6f)\n", params->spin2[0], params->spin2[1], params->spin2[2]);
    }

  count = phenPars.countback;
  params->tC = ((REAL8) count) * dt;

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
      XLAL_ERROR(XLAL_EFAILED);
    }
    else {
      trigAngle.ci   = phenPars.ci;
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

      Psi    = phenPars.Psi;// - 2. * om * log(om);
      Psi0   = Psi + tAs * (om1/mass -dalpha1*trigAngle.ci) * log(1. - t0 / tAs);
      alpha0 = phenPars.alpha + tAs * dalpha1 * log(1. - t0 / tAs);
      //iota0  = acos(phenPars.ci) + diota1 * tAs * log(1. - t0 / tAs);
      energy = phenPars.energy;

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
	XLAL_ERROR(XLAL_EFAILED);
      }

      omegaRD = crealf(modefreqs->data[0]) * unitHz / LAL_PI / 2.;
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

	errcode=XLALSpinInspiralFillH2Modes(h2P2,h2M2,h2P1,h2M1,h20,count,amp22,v,mparams.eta,mparams.dm,Psi,alpha,&trigAngle);

	errcode += XLALSpinInspiralFillH3Modes(h3P3,h3M3,h3P2,h3M2,h3P1,h3M1,h30,count,amp33,v,mparams.eta,mparams.dm,Psi,alpha,&trigAngle);

	errcode += XLALSpinInspiralFillH4Modes(h4P4,h4M4,h4P3,h4M3,h4P2,h4M2,h4P1,h4M1,h40,count,amp44,v,mparams.eta,mparams.dm,Psi,alpha,&trigAngle);

        if (errcode != XLAL_SUCCESS)
          XLAL_ERROR(XLAL_EFUNC);

	fap->data[count] = om;
	phap->data[count] = Psi;

      } while ( (om < frOmRD) && (tim < tAs) );

      XLALDestroyCOMPLEX8Vector(modefreqs);
      origcount=count;
      params->tC = ((REAL8) count) * dt;

      /*--------------------------------------------------------------
       * Attach the ringdown waveform to the end of inspiral
       -------------------------------------------------------------*/

      //printf("time %12.6e  count %d\n",tim,count);

      apcount  = origcount;
      errcode  = XLALPSpinInspiralAttachRingdownWave(h2P2, params, &apcount, nmodes, 2, 2, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h2P2->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h2M2, params, &apcount, nmodes, 2, -2, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h2M2->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h2P1, params, &apcount, nmodes, 2, 1, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h2P1->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h2M1, params, &apcount, nmodes, 2, -1, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h2M1->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h20, params, &apcount, nmodes, 2, 0, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h20->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h3P3, params, &apcount, nmodes, 3, 3, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h3P3->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h3M3, params, &apcount, nmodes, 3, -3, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h3M3->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h3P2, params, &apcount, nmodes, 3, 2, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h3P2->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h3M2, params, &apcount, nmodes, 3, -2, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h3P2->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h3P1, params, &apcount, nmodes, 3, 1, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h3P1->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h3M1, params, &apcount, nmodes, 3, -1, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h3M1->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h30, params, &apcount, nmodes, 3, 0, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h30->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4P4, params, &apcount, nmodes, 4, 4, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h4P4->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4M4, params, &apcount, nmodes, 4, -4, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h4M4->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4P3, params, &apcount, nmodes, 4, 3, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h4P3->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4M3, params, &apcount, nmodes, 4, -3, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h4M3->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4P2, params, &apcount, nmodes, 4, 2, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h4P4->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4M2, params, &apcount, nmodes, 4, -2, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h4M4->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4P1, params, &apcount, nmodes, 4, 1, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h4P3->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h4M1, params, &apcount, nmodes, 4, -1, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h4M3->data[i] = 0.;
      if (apcount > count) count = apcount;

      apcount  = origcount;
      errcode += XLALPSpinInspiralAttachRingdownWave(h40, params, &apcount, nmodes, 4, 0, finalMass, finalSpin);
      for (i = 2 * apcount; i < 2 * length; i++) h40->data[i] = 0.;
      if (apcount > count) count = apcount;

      if (errcode != XLAL_SUCCESS) {
	XLALPrintError("**** LALPSpinInspiralRD ERROR ****: impossible to create RingDownWave\n");
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
	XLAL_ERROR(XLAL_EFAILED);
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
    XLAL_ERROR(XLAL_EFAILED);
  }
  for (i = 0; i < length; i++) {
    x0 = h2P2->data[2 * i];
    x1 = h2P2->data[2 * i + 1];
    x2 = h2M2->data[2 * i];
    x3 = h2M2->data[2 * i + 1];
    sigp->data[i] +=   x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP) + x2 * creal(MultSphHarmM) - x3 * cimag(MultSphHarmM);
    sigc->data[i] += - x0 * cimag(MultSphHarmP) - x1 * creal(MultSphHarmP) - x2 * cimag(MultSphHarmM) - x3 * creal(MultSphHarmM);
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
      sigp->data[i] +=   x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP) + x2 * creal(MultSphHarmM) - x3 * cimag(MultSphHarmM);
      sigc->data[i] += - x0 * cimag(MultSphHarmP) - x1 * creal(MultSphHarmP) - x2 * cimag(MultSphHarmM) - x3 * creal(MultSphHarmM);
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
      sigp->data[i] += x1 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP);
      sigc->data[i] -= x1 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP);
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
      sigp->data[i] += x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP) + x2 * creal(MultSphHarmM) - x3 * cimag(MultSphHarmM);
      sigc->data[i] -= x0 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP) + x2 * cimag(MultSphHarmM) + x3 * creal(MultSphHarmM);
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
      sigp->data[i] += x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP) + x2 * creal(MultSphHarmM) - x3 * cimag(MultSphHarmM);
      sigc->data[i] -= x0 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP) + x2 * cimag(MultSphHarmM) + x3 * creal(MultSphHarmM);
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
      sigp->data[i] += x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP) + x2 * creal(MultSphHarmM) - x3 * cimag(MultSphHarmM);
      sigc->data[i] -= x0 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP) + x2 * cimag(MultSphHarmM) + x3 * creal(MultSphHarmM);
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
      sigp->data[i] += x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmM);
      sigc->data[i] -= x0 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP);
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
      sigp->data[i] += x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP) + x2 * creal(MultSphHarmM) - x3 * cimag(MultSphHarmM);
      sigc->data[i] -= x0 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP) + x2 * cimag(MultSphHarmM) + x3 * creal(MultSphHarmM);
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
      sigp->data[i] += x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP) + x2 * creal(MultSphHarmM) - x3 * cimag(MultSphHarmM);
      sigc->data[i] -= x0 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP) + x2 * cimag(MultSphHarmM) + x3 * creal(MultSphHarmM);
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
      sigp->data[i] += x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP) + x2 * creal(MultSphHarmM) - x3 * cimag(MultSphHarmM);
      sigc->data[i] -= x0 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP) + x2 * cimag(MultSphHarmM) + x3 * creal(MultSphHarmM);
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
      sigp->data[i] += x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP) + x2 * creal(MultSphHarmM) - x3 * cimag(MultSphHarmM);
      sigc->data[i] -= x0 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP) + x2 * cimag(MultSphHarmM) + x3 * creal(MultSphHarmM);
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
      sigp->data[i] += x0 * creal(MultSphHarmP) - x1 * cimag(MultSphHarmP);
      sigc->data[i] -= x0 * cimag(MultSphHarmP) + x1 * creal(MultSphHarmP);
    }
  }

  params->fFinal = params->tSampling / 2.;

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

  return count;

  /*End */
}
