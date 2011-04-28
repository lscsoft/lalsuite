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

/****  <lalVerbatim file="LALPSpinInspiralRDCV">
 * $Id$
 **** </lalVerbatim> */

/****  <lalLaTeX>
 *
 * \subsection{Module \texttt{LALPSpinInspiralRD.c},
 * \texttt{LALPSpinInspiralTemplates} and \texttt{LALPSpinInspiralForInjection}}
 * \label{ss:LALPSpinInspiralRD.c}
 *
 * Module to generate generic spinning binaries waveforms complete with ring-down
 *
 * \subsubsection*{Prototypes}
 * \vspace{0.1in}
 * \input{LALPSpinInspiralRDCP}
 * \idx{\verb&LALPSpinInspiralRD()&}
 * \begin{description}
 * \item {\tt signalvec:} Output containing the inspiral waveform.
 * \item {\tt params:} Input containing binary chirp parameters.
 * \end{description}
 *
 * \input{LALPSpinInspiralRDTemplatesCP}
 * \idx{\verb&LALPSpinInspiralRDTemplates()&}
 * \begin{description}
 * \item {\tt signalvec1:} Output containing the $+$ inspiral waveform.
 * \item {\tt signalvec2:} Output containing the $\times$ inspiral waveform.
 * \item {\tt params:} Input containing binary chirp parameters.
 * \end{description}
 *
 * \input{LALPSpinInspiralRDInjectionCP}
 * \idx{\verb&LALPSpinInspiralRDInjection()&}
 * \begin{description}
 * \item {\tt signalvec:} Output containing the inspiral waveform.
 * \item {\tt params:} Input containing binary chirp parameters.
 * \end{description}
 *
 * \subsubsection*{Description}
 * This codes provide complete waveforms for generically spinning binary systems.
 * In order to construct the waveforms three phases are joined together:
 * an initial inspiral phase, a phenomenological phase encompassing the description
 * of the merger and the ring-down of the final black hole.
 * During the inspiral phase the system is evolved according to the standard
 * PN formulas, valid up to 3.5PN for the orbital motion,
 * to 2.5PN level for spin-orbital momentum and to 2PN for spin-spin contributions
 * to the orbital phase.
 * Then a phenomenological phase is added during which the frequency of the
 * waveform has a pole-like behaviour. The stitching is performed in order to
 * ensure continuity of the phase, the frequency and its first and second
 * derivatives. Finally a ring-down phase is attached.
 *
 * \subsubsection*{Algorithm}
 *
 * \subsubsection*{Uses}
 * \begin{verbatim}
 * LALPSpinInspiralRDderivatives
 * LALInspiralSetup
 * LALInspiralChooseModel
 * LALRungeKutta4
 * \end{verbatim}
 *
 * \subsubsection*{Notes}
 *
 * \vfill{\footnotesize\input{LALPSpinInspiralRDCV}}
 *
 **** </lalLaTeX>  */

/** \defgroup psird Complete phenomenological spin-inspiral waveforms
 *
 * This code provides complete waveforms for generically spinning binary
 * systems.
 *
 */

#include <LALPSpinInspiralRD.h>
#include <LALAdaptiveRungeKutta4.h>

#include <lal/Units.h>
#include <lal/LALInspiral.h>
#include <lal/SeqFactories.h>
#include <lal/NRWaveInject.h>
#include <lal/RealFFT.h>

NRCSID(LALPSPININSPIRALRDC, "$Id$");

#define omM0 0.065

typedef struct LALPSpinInspiralRDstructparams {
    REAL8 eta;			///< symmetric mass ratio
    REAL8 dm;			///< \f$m_1-m_2\f$
    REAL8 m1m2;			///< \f$m_1/m_2\f$
    REAL8 m2m1;			///< \f$m_2/m_1\f$
    REAL8 m2m;
    REAL8 m1m;
    REAL8 wdotorb[8];		///< Coefficients of the analytic PN expansion of \f$\dot\omega_{orb}\f$
    REAL8 wdotorblog;		///< Log coefficient of the PN expansion of of \f$\dot\omega_{orb}\f$
    REAL8 wdotspin15S1Lh;
    REAL8 wdotspin15S2Lh;
    REAL8 wdotspin20S1S2;
    REAL8 wdotspin20S1S1;	///< Coeff. of the \f$s_1s_1\f$ cntrb. to \f$\dot\omega\f$
    REAL8 wdotspin20S1S2Lh;
    REAL8 wdotspin25S1Lh;
    REAL8 wdotspin25S2Lh;	///< Coeff. of the \f$s_2\cdot \hat L_N\f$ cntrb. to \f$\dot\omega\f$
    REAL8 S1dot15;
    REAL8 S2dot15;
    REAL8 Sdot20;
    REAL8 S1dot25;
    REAL8 S2dot25;
    REAL8 Lhdot15;
    REAL8 Lhdot20;
    REAL8 epnorb[4];		///< Coefficients of the PN expansion of the energy
    REAL8 epnspin15S1dotLh;	///< Coeff. of the \f$S_1\cdot L\f$ term in energy
    REAL8 epnspin15S2dotLh;	///< Coeff. of the \f$S_2\cdot L\f$ term in energy
    REAL8 epnspin20S1S2;	///< Coeff. of the \f$S_1\cdot S_2\f$ term in energy
    REAL8 epnspin20S1S2dotLh;	///< Coeff. of the \f$S_{1,2}\cdot L\f$ term in energy
    REAL8 epnspin20S1S1;	///< Coeff. of the \f$S_1\cdot S_1\f$ term in energy
    REAL8 epnspin20S1S1dotLh;
    REAL8 epnspin20S2S2;	///< Coeff. of the \f$S_2\cdot S_2\f$ term in energy
    REAL8 epnspin20S2S2dotLh;
    REAL8 epnspin25S1dotLh;
    REAL8 epnspin25S2dotLh;

} LALPSpinInspiralRDparams;

static void XLALPSpinInspiralRDAdaptiveSetParams(LALPSpinInspiralRDparams *mparams,InspiralTemplate *params, InspiralInit *paramsInit) {

  UINT4 j;

  /* setup coefficients for PN equations */
  mparams->m2m1 = params->mass2 / params->mass1;
  mparams->m1m2 = params->mass1 / params->mass2;
  mparams->m1m = params->mass1 / params->totalMass;
  mparams->m2m = params->mass2 / params->totalMass;
  mparams->dm = (params->mass1 - params->mass2) / params->totalMass;
  
  /* params->eta might have been set up before but just for safety, we
   * recompute it here below.*/
  params->eta = (params->mass1 * params->mass2) / (params->mass1 + params->mass2) / (params->mass1 + params->mass2);
  mparams->eta = params->eta;
  
  for (j = LAL_PNORDER_NEWTONIAN; j <= params->order; j++) {
    mparams->wdotorb[j] = paramsInit->ak.ST[j];
  }
  mparams->wdotorblog = 0.;
  for (j = params->order + 1; j < 8; j++) {
    mparams->wdotorb[j] = 0.;
  }
  if ((params->order) >= 6) {
    mparams->wdotorblog = paramsInit->ak.ST[7];
    if ((params->order) == 7)
      mparams->wdotorb[7] = paramsInit->ak.ST[8];
  }
  
  mparams->epnorb[0] = paramsInit->ak.ETaN;
  
  switch (params->order) {
  case LAL_PNORDER_NEWTONIAN:
  case LAL_PNORDER_HALF:
    break;
  case LAL_PNORDER_ONE:
    mparams->epnorb[1] = paramsInit->ak.ETa1;
    break;
  case LAL_PNORDER_ONE_POINT_FIVE:
    
    mparams->epnorb[1] = paramsInit->ak.ETa1;
    mparams->epnspin15S1dotLh = 8. / 3. + 2. * mparams->m2m1;
    mparams->epnspin15S2dotLh = 8. / 3. + 2. * mparams->m1m2;
    
    mparams->wdotspin15S1Lh = -(113.0 + 75.0 * mparams->m2m1) / 12.0;
    mparams->wdotspin15S2Lh = -(113.0 + 75.0 * mparams->m1m2) / 12.0;
    
    mparams->Lhdot15 = 0.5;
    
    mparams->S1dot15 = (4.0 + 3.0 * mparams->m2m1) / 2.0 * mparams->eta;
    mparams->S2dot15 = (4.0 + 3.0 * mparams->m1m2) / 2.0 * mparams->eta;

    mparams->Lhdot15 = 0.5;
    mparams->S1dot15 = (4.0 + 3.0 * mparams->m2m1) / 2.0 * mparams->eta;
    mparams->S2dot15 = (4.0 + 3.0 * mparams->m1m2) / 2.0 * mparams->eta;
    
    break;
  case LAL_PNORDER_TWO:
    
    mparams->epnorb[1] = paramsInit->ak.ETa1;
    mparams->epnspin15S1dotLh = 8. / 3. + 2. * mparams->m2m1;
    mparams->epnspin15S2dotLh = 8. / 3. + 2. * mparams->m1m2;
    mparams->epnorb[2] = paramsInit->ak.ETa2;
    
    mparams->wdotspin15S1Lh = -(113.0 + 75.0 * mparams->m2m1) / 12.0;
    mparams->wdotspin15S2Lh = -(113.0 + 75.0 * mparams->m1m2) / 12.0;
    mparams->wdotspin20S1S2 = -(1.0 / 48.0) / mparams->eta;
    mparams->wdotspin20S1S1 = 1. / 96.;
    
    mparams->Lhdot15 = 0.5;
    mparams->Lhdot20 = -1.5 / mparams->eta;
    
    mparams->S1dot15 = (4.0 + 3.0 * mparams->m2m1) / 2.0 * mparams->eta;
    mparams->S2dot15 = (4.0 + 3.0 * mparams->m1m2) / 2.0 * mparams->eta;
    mparams->Sdot20 = 0.5;
    break;
    
  case LAL_PNORDER_TWO_POINT_FIVE:
    
    mparams->epnorb[1] = paramsInit->ak.ETa1;
    mparams->epnspin15S1dotLh = 8. / 3. + 2. * mparams->m2m1;
    mparams->epnspin15S2dotLh = 8. / 3. + 2. * mparams->m1m2;
    mparams->epnorb[2] = paramsInit->ak.ETa2;
    mparams->epnspin25S1dotLh = 8. - 31. / 9. * mparams->eta + (3. - 10. / 3. * mparams->eta) * mparams->m2m1;
    mparams->epnspin25S2dotLh = 8. - 31. / 9. * mparams->eta + (3. - 10. / 3. * mparams->eta) * mparams->m1m2;

    mparams->wdotspin15S1Lh = -(113.0 + 75.0 * mparams->m2m1) / 12.0;
    mparams->wdotspin15S2Lh = -(113.0 + 75.0 * mparams->m1m2) / 12.0;
    mparams->wdotspin20S1S2 = -(1.0 / 48.0) / mparams->eta;
    mparams->wdotspin20S1S1 = 1. / 96.;
    mparams->wdotspin25S1Lh = -31319. / 1008. + 1159. / 24. * mparams->eta + (-809. / 84. + 281. / 8. * mparams->eta) * mparams->m2m1;
    mparams->wdotspin25S2Lh = -31319. / 1008. + 1159. / 24. * mparams->eta + (-809. / 84. + 281. / 8. * mparams->eta) * mparams->m1m2;

    mparams->Lhdot15 = 0.5;
    mparams->Lhdot20 = -1.5 / mparams->eta;

    mparams->S1dot15 = (4.0 + 3.0 * mparams->m2m1) / 2.0 * mparams->eta;
    mparams->S2dot15 = (4.0 + 3.0 * mparams->m1m2) / 2.0 * mparams->eta;
    mparams->Sdot20 = 0.5;
    mparams->S1dot25 = 0.5625 + 1.25 * mparams->eta - mparams->eta * mparams->eta / 24. + mparams->dm * (-0.5625 + 0.625 * mparams->eta);
    mparams->S2dot25 = 0.5625 + 1.25 * mparams->eta - mparams->eta * mparams->eta / 24. - mparams->dm * (-0.5625 + 0.625 * mparams->eta);
	break;

    case LAL_PNORDER_THREE:
    case LAL_PNORDER_THREE_POINT_FIVE:
      mparams->epnorb[1] = paramsInit->ak.ETa1;
      mparams->epnspin15S1dotLh = 8. / 3. + 2. * mparams->m2m1;
      mparams->epnspin15S2dotLh = 8. / 3. + 2. * mparams->m1m2;
      mparams->epnorb[2] = paramsInit->ak.ETa2;
      mparams->epnspin20S1S2 = 1. / mparams->eta;
      mparams->epnspin20S1S2dotLh = -3. / mparams->eta;
      mparams->epnspin20S1S1 = (1. + mparams->m2m1) * (1. + mparams->m2m1) / 2.;
      mparams->epnspin20S2S2 = (1. + mparams->m1m2) * (1. + mparams->m1m2) / 2.;
      mparams->epnspin20S1S1dotLh = -3. * (1. + mparams->m2m1) * (1. + mparams->m2m1) / 2.;
      mparams->epnspin20S2S2dotLh = -3. * (1. + mparams->m1m2) * (1. + mparams->m1m2) / 2.;
      mparams->epnspin25S1dotLh = 8. - 31. / 9. * mparams->eta + (3. - 10. / 3. * mparams->eta) *  mparams->m2m1;
      mparams->epnspin25S2dotLh = 8. - 31. / 9. * mparams->eta + (3. - 10. / 3. * mparams->eta) *  mparams->m1m2;
      mparams->epnorb[3] = paramsInit->ak.ETa3;
      
      mparams->wdotspin15S1Lh = -(113.0 + 75.0 * mparams->m2m1) / 12.0;
      mparams->wdotspin15S2Lh = -(113.0 + 75.0 * mparams->m1m2) / 12.0;
      mparams->wdotspin20S1S2 = -(1.0 / 48.0) / mparams->eta;
      mparams->wdotspin20S1S1 = 1. / 96.;
      mparams->wdotspin25S1Lh = -31319. / 1008. + 1159. / 24. * mparams->eta + (-809. / 84. + 281. / 8. * mparams->eta) * mparams->m2m1;
	mparams->wdotspin25S2Lh = -31319. / 1008. + 1159. / 24. * mparams->eta + (-809. / 84. + 281. / 8. * mparams->eta) * mparams->m1m2;
	mparams->S1dot15 = (4.0 + 3.0 * mparams->m2m1) / 2.0 * mparams->eta;
	mparams->S2dot15 = (4.0 + 3.0 * mparams->m1m2) / 2.0 * mparams->eta;
	mparams->Sdot20 = 0.5;
	mparams->S1dot25 = 0.5625 + 1.25 * mparams->eta - mparams->eta * mparams->eta / 24. + mparams->dm * (-0.5625 + 0.625 * mparams->eta);
	mparams->S2dot25 = 0.5625 + 1.25 * mparams->eta - mparams->eta * mparams->eta / 24. - mparams->dm * (-0.5625 + 0.625 * mparams->eta);
	break;
    case LAL_PNORDER_PSEUDO_FOUR:
	fprintf(stderr,
		"*** LALPhenSpinInspiralRD ERROR: PhenSpin approximant not available at pseudo4PN order\n");
	break;
    case LAL_PNORDER_NUM_ORDER:
	fprintf(stderr,
		"*** LALPhenSpinInspiralRD ERROR: NUM_ORDER not a valid PN order\n");
    }
}

#define UNUSED(expr) do { (void)(expr); } while (0)

/*
 *
 * function to set derivatives: values and mparams input, dvalues output
 *
 */

static int XLALPSpinInspiralRDAdaptiveTest(double t,const double values[],double dvalues[],void *mparams) {
  REAL8 omega, v, v2, v3, v4, v5, energy, energyold;
  REAL8 Lhx, Lhy, Lhz, S1x, S1y, S1z, S2x, S2y, S2z;
  REAL8 LhS1, LhS2,S1S2,S1S1,S2S2;
  LALPSpinInspiralRDparams *params = (LALPSpinInspiralRDparams*)mparams;

  UNUSED(t);

  omega = values[1];	
  v = pow(omega,oneby3);
  v2 = v*v;
  v3 = omega;
  v4 = omega*v;
  v5 = omega*v2;

  // energy=-eta/2 * v^2 * [1-(9+\eta)/12 v^2 +...] up to 3PN without spin effects
  energy = (1. + v2 * (params->epnorb[1] +
		       v2 * (params->epnorb[2] +
			     v2 * params->epnorb[3])));
	
  Lhx = values[2]; 
  Lhy = values[3]; 
  Lhz = values[4];
  S1x  = values[5]; 
  S1y  = values[6]; 
  S1z  = values[7];
  S2x  = values[8]; 
  S2y  = values[9];
  S2z  = values[10];
  energyold = values[11];
		
  LhS1 = (Lhx*S1x + Lhy*S1y + Lhz*S1z);
  LhS2 = (Lhx*S2x + Lhy*S2y + Lhz*S2z);		
			
  S1S2 = (S1x*S2x + S1y*S2y + S1z*S2z);
  S1S1 = (S1x*S1x + S1y*S1y + S1z*S1z);
  S2S2 = (S2x*S2x + S2y*S2y + S2z*S2z);

  energy += v3 * (params->epnspin15S1dotLh * LhS1 + params->epnspin15S2dotLh * LhS2);	// see e.g. Blanchet et al. gr-qc/0605140		

  energy += v4 * (params->epnspin20S1S2 * S1S2 + params->epnspin20S1S2dotLh * LhS1 * LhS2);	// see e.g. Buonanno et al. as above
  energy += v4 * (params->epnspin20S1S1 * S1S1 + params->epnspin20S2S2 * S2S2 + params->epnspin20S1S1dotLh * LhS1 * LhS1 + params->epnspin20S2S2 * LhS2 * LhS2);	// see Racine et al. as above

  energy += v5 * (params->epnspin25S1dotLh * LhS1 + params->epnspin25S2dotLh * LhS2);	//see (7.9) of Blanchet et al.

  if (energy<0.) printf("En %12.6e  om=%12.6e\n",energy,omega);
  energy *= params->epnorb[0] * v2;	

  if ( (energy > 0.0) || ( energy > 0.99*energyold) ) { 
    /* energy increase*/
    printf("En increases %12.6e  %12.6e  %11.4e\n",energy,energyold,params->epnorb[0]);
    return LALPSIRDPN_TEST_ENERGY;
  } 
  else if (dvalues[1] < 0.0) {
    /* omegadot < 0 */
    return LALPSIRDPN_TEST_OMEGADOT;
  } 
  else if (isnan(omega)) {
    /* omega is nan */
    return LALPSIRDPN_TEST_OMEGANAN;
  } else 
    if (omega<omM0) {
      return GSL_SUCCESS;
    }
    else { 
      printf("omM %12.6e  omega %12.6e\n",omM0,omega);
      return LALPSIRDPN_TEST_OMEGAMATCH;
    }
}

/**
 * \ingroup psird
 * \brief Module to compute detivative of dynamical variables
 */

static int XLALPSpinInspiralRDAdaptiveDerivatives(double t, const double values[], double dvalues[], void *mparams)
{

    REAL8 Phi;			// half of the main GW phase, this is \f$Phi\f$ of eq.3.11 of arXiv:0810.5336
    REAL8 omega;		// time-derivative of the orbital phase
    REAL8 Lhx, Lhy, Lhz;	// orbital angolar momentum unit vector
    REAL8 S1x, S1y, S1z;	// dimension-less spin variable S/M^2
    REAL8 S2x, S2y, S2z;
    REAL8 LhS1, LhS2;		// scalar products
    REAL8 domega;		// derivative of omega
    REAL8 dLhx, dLhy, dLhz;	// derivatives of \f$\hat L_N\f$ components
    REAL8 dS1x, dS1y, dS1z;	// derivative of \f$S_i\f$
    REAL8 dS2x, dS2y, dS2z;

    /* auxiliary variables*/
    REAL8 S1S2, S1S1, S2S2;	// Scalar products
    REAL8 alphadotcosi;		// alpha is the right ascension of L, i(iota) the angle between L and J
    REAL8 v, v2, v3, v4, v5, v6, v7;
    REAL8 tmpx, tmpy, tmpz, cross1x, cross1y, cross1z, cross2x, cross2y, cross2z, Lhxy;

    LALPSpinInspiralRDparams *params = (LALPSpinInspiralRDparams *) mparams;

    UNUSED(t);

    /* --- computation start here --- */
    Phi = values[0];
    omega = values[1];

    if (omega < 0.0) {
      fprintf(stderr,"WARNING: Omega has become -ve, this should lead to nan's \n");
      return LALPSIRD_DERIVATIVE_OMEGANONPOS;
    }

    Lhx = values[2];
    Lhy = values[3];
    Lhz = values[4];

    S1x = values[5];
    S1y = values[6];
    S1z = values[7];

    S2x = values[8];
    S2y = values[9];
    S2z = values[10];

    v = pow(omega, oneby3);
    v2 = v * v;
    v3 = omega;
    v4 = omega * v;
    v5 = omega * v2;
    v6 = omega * omega;
    v7 = omega * v;

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


    // Adding spin effects
    // L dot S1,2
    LhS1 = (Lhx * S1x + Lhy * S1y + Lhz * S1z);
    LhS2 = (Lhx * S2x + Lhy * S2y + Lhz * S2z);

    // wdotspin15SiLh = -1/12 (...)
    domega += v3 * (params->wdotspin15S1Lh * LhS1 + params->wdotspin15S2Lh * LhS2);	// see e.g. Buonanno et al. gr-qc/0211087

    // wdotspin20S1S1 = -1/48 eta
    S1S1 = (S1x * S1x + S1y * S1y + S1z * S1z);
    S2S2 = (S2x * S2x + S2y * S2y + S2z * S2z);
    S1S2 = (S1x * S2x + S1y * S2y + S1z * S2z);
    domega += params->wdotspin20S1S2 * v4 * (247.0 * S1S2 - 721.0 * LhS1 * LhS2);	// see e.g. Buonanno et al. arXiv:0810.5336
    domega += params->wdotspin20S1S1 * v4 * (719. * (LhS1 * LhS1 + LhS2 * LhS2) - 233. * (S1S1 + S2S2));	// see Racine et al. arXiv:0812.4413

    // wdotspin25SiLh = see below
    domega += v5 * (params->wdotspin25S1Lh * LhS1 + params->wdotspin25S2Lh * LhS2);	//see (8.3) of Blanchet et al.

    // Setting the right pre-factor

    domega *= 96. / 5. * params->eta * v5 * omega* omega;

    /*Derivative of the angular momentum and spins */

    /* dS1, 1.5PN */
    /* S1dot15= (4+3m2/m1)/2 * eta */

    cross1x = (Lhy * S1z - Lhz * S1y);
    cross1y = (Lhz * S1x - Lhx * S1z);
    cross1z = (Lhx * S1y - Lhy * S1x);

    cross2x = (Lhy * S2z - Lhz * S2y);
    cross2y = (Lhz * S2x - Lhx * S2z);
    cross2z = (Lhx * S2y - Lhy * S2x);

    dS1x = params->S1dot15 * v5 * cross1x;
    dS1y = params->S1dot15 * v5 * cross1y;
    dS1z = params->S1dot15 * v5 * cross1z;

    /* dS1, 2PN */
    /* Sdot20= 0.5 */
    tmpx = S1z * S2y - S1y * S2z;
    tmpy = S1x * S2z - S1z * S2x;
    tmpz = S1y * S2x - S1x * S2y;

    // S1S2 contribution
    dS1x += params->Sdot20 * v6 * (tmpx - 3. * LhS2 * cross1x);
    dS1y += params->Sdot20 * v6 * (tmpy - 3. * LhS2 * cross1y);
    dS1z += params->Sdot20 * v6 * (tmpz - 3. * LhS2 * cross1z);
    // S1S1 contribution
    dS1x -= 3. * params->Sdot20 * v6 * LhS1 * cross1x * (1. + params->m2m1) * params->m2m;
    dS1y -= 3. * params->Sdot20 * v6 * LhS1 * cross1y * (1. + params->m2m1) * params->m2m;
    dS1z -= 3. * params->Sdot20 * v6 * LhS1 * cross1z * (1. + params->m2m1) * params->m2m;

    // dS1, 2.5PN, eq. 7.8 of Blanchet et al. gr-qc/0605140
    // S1dot25= 9/8-eta/2.+eta+mparams->eta*29./24.+mparams->m1m2*(-9./8.+5./4.*mparams->eta)
    dS1x += params->S1dot25 * v7 * cross1x;
    dS1y += params->S1dot25 * v7 * cross1y;
    dS1z += params->S1dot25 * v7 * cross1z;

    /* dS2, 1.5PN */
    dS2x = params->S2dot15 * v5 * cross2x;
    dS2y = params->S2dot15 * v5 * cross2y;
    dS2z = params->S2dot15 * v5 * cross2z;

    /* dS2, 2PN */
    dS2x += params->Sdot20 * v6 * (-tmpx - 3.0 * LhS1 * cross2x);
    dS2y += params->Sdot20 * v6 * (-tmpy - 3.0 * LhS1 * cross2y);
    dS2z += params->Sdot20 * v6 * (-tmpz - 3.0 * LhS1 * cross2z);
    // S2S2 contribution
    dS2x -= 3. * params->Sdot20 * v6 * LhS2 * cross2x * (1. + params->m1m2) * params->m1m;
    dS2y -= 3. * params->Sdot20 * v6 * LhS2 * cross2y * (1. + params->m1m2) * params->m1m;
    dS2z -= 3. * params->Sdot20 * v6 * LhS2 * cross2z * (1. + params->m1m2) * params->m1m;

    // dS2, 2.5PN, eq. 7.8 of Blanchet et al. gr-qc/0605140
    dS2x += params->S2dot25 * v7 * cross2x;
    dS2y += params->S2dot25 * v7 * cross2y;
    dS2z += params->S2dot25 * v7 * cross2z;

    dLhx = -(dS1x + dS2x) * v / params->eta;
    dLhy = -(dS1y + dS2y) * v / params->eta;
    dLhz = -(dS1z + dS2z) * v / params->eta;

    /* dphi */
    Lhxy = Lhx * Lhx + Lhy * Lhy;

    if (Lhxy > 0.0)
	alphadotcosi = Lhz * (Lhx * dLhy - Lhy * dLhx) / Lhxy;
    else
	alphadotcosi = 0.;

    /* dvalues->data[0] is the phase derivative */
    /* omega is the derivative of the orbital phase omega \neq dvalues->data[0] */

    dvalues[0] = omega - alphadotcosi;
    dvalues[1] = domega;

    dvalues[2] = dLhx;
    dvalues[3] = dLhy;
    dvalues[4] = dLhz;

    dvalues[5] = dS1x;
    dvalues[6] = dS1y;
    dvalues[7] = dS1z;

    dvalues[8] = dS2x;
    dvalues[9] = dS2y;
    dvalues[10] = dS2z;

    return GSL_SUCCESS;

} /* end of LALPSpinInspiralRDderivatives */


/**
 * \ingroup psird
 * \brief Main module to produce waveforms 
 */

void LALPSpinInspiralRD(LALStatus * status, REAL4Vector * signalvec, InspiralTemplate * params)
{

  UINT4 count=0;
  InspiralInit paramsInit;
  
  INITSTATUS(status, "LALPSpinInspiralRD", LALPSPININSPIRALRDC);
  ATTATCHSTATUSPTR(status);
  
  ASSERT(signalvec,       status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(signalvec->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params,          status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params->nStartPad >= 0,  status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->nEndPad >= 0,    status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->fLower > 0,      status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->tSampling > 0,   status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->totalMass > 0.,  status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  LALInspiralSetup(status->statusPtr, &(paramsInit.ak), params);
  CHECKSTATUSPTR(status);
  LALInspiralChooseModel(status->statusPtr, &(paramsInit.func),&(paramsInit.ak), params);
  CHECKSTATUSPTR(status);
  
  REAL8Vector *s=XLALCreateREAL8Vector(signalvec->length);
  memset(s->data, 0, s->length * sizeof(REAL8));

  /* Call the engine function */
  XLALPSpinInspiralRDEngine(s, NULL, NULL, NULL, NULL, &count, params, &paramsInit);

  UINT4 i;
  for (i=0;i<s->length;i++) {
    signalvec->data[i]=(REAL4) s->data[i];
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);

}


NRCSID(LALPSPININSPIRALRDTEMPLATESC, "$Id$");

/**
 * \ingroup psird
 * \brief Module to produce waveform templates
 */

void LALPSpinInspiralRDTemplates(LALStatus * status,
				 REAL4Vector * signalvec1,
				 REAL4Vector * signalvec2,
				 InspiralTemplate * params)
{
    UINT4 count=0;
    InspiralInit paramsInit;

    INITSTATUS(status, "LALPSpinInspiralRDTemplates",LALPSPININSPIRALRDTEMPLATESC);
    ATTATCHSTATUSPTR(status);

    ASSERT(signalvec1,       status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
    ASSERT(signalvec1->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
    ASSERT(signalvec2,       status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
    ASSERT(signalvec2->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
    ASSERT(params,           status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
    ASSERT(params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
    ASSERT(params->nEndPad >= 0,   status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
    ASSERT(params->fLower > 0,     status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
    ASSERT(params->tSampling > 0,  status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
    ASSERT(params->totalMass > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

    LALInspiralSetup(status->statusPtr, &(paramsInit.ak), params);
    CHECKSTATUSPTR(status);
    LALInspiralChooseModel(status->statusPtr, &(paramsInit.func), &(paramsInit.ak), params);
    CHECKSTATUSPTR(status);

    REAL8Vector* s1=XLALCreateREAL8Vector(signalvec1->length);
    REAL8Vector* s2=XLALCreateREAL8Vector(signalvec2->length);

    memset(s1->data, 0, signalvec1->length * sizeof(REAL8));
    memset(s2->data, 0, signalvec2->length * sizeof(REAL8));

    XLALPSpinInspiralRDEngine(s1, s2, NULL, NULL, NULL, &count, params, &paramsInit);

    UINT4 i;
    for (i=0;i<s1->length;i++) {
      signalvec1->data[i]=(REAL4) s1->data[i];
    }
    for (i=0;i<s2->length;i++) {
      signalvec2->data[i]=(REAL4) s2->data[i];
    }

    XLALDestroyREAL8Vector(s1);
    XLALDestroyREAL8Vector(s2);

    DETATCHSTATUSPTR(status);
    RETURN(status);
}


NRCSID(LALPSPININSPIRALRDINJECTIONC, "$Id$");

/**
 * \ingroup psird
 * \brief Module to produce injection waveforms
 */

void LALPSpinInspiralRDForInjection(LALStatus        * status,
				    CoherentGW       * waveform,
				    InspiralTemplate * params,
				    PPNParamStruc    * ppnParams)
{

    REAL8Vector *h = NULL;	/* pointers to generated amplitude  data */
    REAL8Vector *f = NULL;	/* pointers to generated  frequency data */
    REAL8Vector *phi = NULL;	/* pointer to generated phase data */

    InspiralInit paramsInit;
    UINT4 nbins,count,i;


    INITSTATUS(status, "LALPSpinInspiralRDInjection", LALPSPININSPIRALRDINJECTIONC);
    ATTATCHSTATUSPTR(status);

    /* Make sure parameter and waveform structures exist. */
    ASSERT(params,   status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
    ASSERT(waveform, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
    /* Make sure waveform fields don't exist. */
    ASSERT(!(waveform->a),   status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
    ASSERT(!(waveform->f),   status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
    ASSERT(!(waveform->phi), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
    ASSERT(!(waveform->h),   status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

    /* Compute some parameters */
    LALInspiralInit(status->statusPtr, params, &paramsInit);
    CHECKSTATUSPTR(status);
    if (paramsInit.nbins == 0) {
	DETATCHSTATUSPTR(status);
	RETURN(status);
    }

    nbins = 2 * paramsInit.nbins;

    /* Now we can allocate memory and vector for coherentGW structure */
    f = XLALCreateREAL8Vector(nbins);
    h = XLALCreateREAL8Vector(2 * nbins);
    phi = XLALCreateREAL8Vector(nbins);

    /* By default the waveform is empty */
    memset(f->data,   0,     nbins * sizeof(REAL8));
    memset(h->data,   0, 2 * nbins * sizeof(REAL8));
    memset(phi->data, 0,     nbins * sizeof(REAL8));

    /* Call the engine function */
    XLALPSpinInspiralRDEngine(NULL, NULL, h, f, phi, &count, params, &paramsInit);

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
	ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
    }
    memset(waveform->h, 0, sizeof(REAL4TimeVectorSeries));

    if ((waveform->f = (REAL4TimeSeries *)
	 LALMalloc(sizeof(REAL4TimeSeries))) == NULL) {
	LALFree(waveform->h);
	waveform->a = NULL;
	ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
    }
    memset(waveform->f, 0, sizeof(REAL4TimeSeries));

    if ((waveform->phi = (REAL8TimeSeries *)
	 LALMalloc(sizeof(REAL8TimeSeries))) == NULL) {
	LALFree(waveform->h);
	waveform->h = NULL;
	LALFree(waveform->f);
	waveform->f = NULL;
	ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
    }
    memset(waveform->phi, 0, sizeof(REAL8TimeSeries));

    waveform->h->data = XLALCreateREAL4VectorSequence(count, 2);
    waveform->f->data = XLALCreateREAL4Vector(count);
    waveform->phi->data = XLALCreateREAL8Vector(count);

    for (i=0; i<count; i++) {
      waveform->f->data->data[i]=(REAL4) f->data[i];
      waveform->h->data->data[i]=(REAL4) h->data[i];
      waveform->h->data->data[i+nbins]=(REAL4) h->data[i+nbins];
    }

    memcpy(waveform->phi->data->data, phi->data, count * (sizeof(REAL8)));

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
    ppnParams->tc = (double) (count - 1) / params->tSampling;
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

    DETATCHSTATUSPTR(status);
    RETURN(status);

}				/* End LALPSpinInspiralRDForInjection */

void LALPSpinInspiralRDFreqDom(LALStatus * status,
			       REAL4Vector * signalvec,
			       InspiralTemplate * params)
{

    REAL8Vector *tsignalvec = NULL;
    REAL4Vector *fsignalvec = NULL;
    REAL4FFTPlan *forwPlan = NULL;

    InspiralInit paramsInit;

    UINT4 nbins,count,idx;
    UINT4 length = signalvec->length;

    INITSTATUS(status, "LALPSpinInspiralRDFReqDom", LALPSPININSPIRALRDC);
    ATTATCHSTATUSPTR(status);

    ASSERT(signalvec,       status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
    ASSERT(signalvec->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
    ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
    ASSERT(params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
    ASSERT(params->nEndPad >= 0,   status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
    ASSERT(params->fLower > 0,     status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
    ASSERT(params->tSampling > 0,  status,  LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
    ASSERT(params->totalMass > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

    LALInspiralInit(status->statusPtr, params, &paramsInit);
    CHECKSTATUSPTR(status);

    if (paramsInit.nbins == 0) {
	DETATCHSTATUSPTR(status);
	RETURN(status);
    }
    nbins = paramsInit.nbins;
    if (nbins < signalvec->length) nbins = signalvec->length;
    tsignalvec = XLALCreateREAL8Vector(nbins);
    fsignalvec = XLALCreateREAL4Vector(nbins);

    memset(signalvec->data, 0, length * sizeof(REAL4));
    memset(tsignalvec->data, 0, nbins * sizeof(REAL8));
    memset(fsignalvec->data, 0, nbins * sizeof(REAL8));

    /* Call the engine function */

    XLALPSpinInspiralRDEngine(tsignalvec, NULL, NULL, NULL, NULL, &count, params, &paramsInit);
    REAL4Vector *tsigR4=XLALCreateREAL4Vector(nbins);
    for (idx=0;idx<nbins;idx++) tsigR4->data[idx]=(REAL4) tsignalvec->data[idx];
    XLALDestroyREAL8Vector(tsignalvec);
    XLALInspiralWaveTaper(tsigR4, 3);

    forwPlan = XLALCreateForwardREAL4FFTPlan(nbins, 0);
    if (forwPlan == NULL) {
      XLALDestroyREAL4Vector(fsignalvec);
      XLALDestroyREAL4Vector(tsigR4);
      ABORTXLAL(status);
    }
    XLALREAL4VectorFFT(fsignalvec, tsigR4, forwPlan);
    XLALDestroyREAL4Vector(tsigR4);

    for (idx = 0; idx < nbins/2; idx++) fsignalvec->data[idx]/=params->tSampling/2.;

    if (nbins>signalvec->length) {
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

      gsl_interp_accel *acc    = (gsl_interp_accel*) gsl_interp_accel_alloc();
      gsl_spline* spline_real = (gsl_spline*) gsl_spline_alloc(gsl_interp_cspline, nbins/2);
      gsl_spline* spline_imag = (gsl_spline*) gsl_spline_alloc(gsl_interp_cspline, nbins/2);
      gsl_spline_init(spline_real, freqSup->data, fsigRe->data, nbins/2);
      gsl_spline_init(spline_imag, freqSup->data, fsigIm->data, nbins/2);

      for (idx=0;idx<length/2;idx++){
	    signalvec->data[idx] = gsl_spline_eval ( spline_real , freq->data[idx] , acc);
	    signalvec->data[length-idx] = gsl_spline_eval ( spline_imag , freq->data[idx] , acc);

      }
      signalvec->data[0] = 0.;
      signalvec->data[nbins / 2] = 0.;

      XLALDestroyREAL8Vector(fsigRe);
      XLALDestroyREAL8Vector(fsigIm);
      XLALDestroyREAL8Vector(freq);
      XLALDestroyREAL8Vector(freqSup);

    }

    XLALDestroyREAL4Vector(fsignalvec);
    XLALDestroyREAL4FFTPlan(forwPlan);

    DETATCHSTATUSPTR(status);
    RETURN(status);
}


/*
 *
 * Main function
 *
 */

NRCSID(LALPSPININSPIRALRDENGINEC, "$Id$");

/**
 * \ingroup psird
 * \brief Module actually computing PSIRD waveforms
 */


int XLALPSpinInspiralRDEngine(
			      REAL8Vector * signalvec1,
			      REAL8Vector * signalvec2,
			      REAL8Vector * hh,
			      REAL8Vector * ff,
			      REAL8Vector * phi,
			      UINT4 * count,
			      InspiralTemplate * params,
			      InspiralInit     * paramsInit)
{

  /* declare code parameters and variables */
  const INT4 nn = 11;	          // number of dynamical variables
  UINT4 apcount;                  // integration steps performed
  UINT4 length;		          // signal vector length
  UINT4 i, j, k, l;		  // counters
  
  expnCoeffs ak;		  // Coefficients in a generic PN expansion (E, flux...)
  
  REAL8 lengths;		  // Wf length in seconds
  REAL8 v = 0.;
  REAL8 v2 = 0.;
  REAL8 v2old;
  REAL8 mass;			  // Total mass in SI units
  REAL8 tim;			  // time (units of total mass)
  REAL8 unitHz;
  REAL8 dt;
  REAL8 initomega,initphi;
  REAL8 inc;
  REAL8 Lhmag,initJmag;
  REAL8 initS1[3],initS2[3],initLh[3],initJ[3];
  REAL8 iS1[3],iS2[3];
  REAL8 phiJ,thetaJ;
  REAL8 ry[3][3],rz[3][3];

  REAL8Array *yout;               // Auxiliary variable for integration
  ark4GSLIntegrator *integrator;
  UINT4 intlen;
  UINT4 kend=0;
  UINT4 jend=0;
  INT4  intreturn;
  REAL8 yinit[nn];

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

  REAL8 Psi;
  REAL8 amp22ini,amp22,amp20,amp33,amp32,amp31,amp30,amp44,amp43,amp42,amp41,amp40;
  REAL8 ci,si,c2i,s2i,cdi,sdi;
  REAL8 ci2,si2,c2i2,s2i2,c3i2,s3i2,c4i2,s4i2,c5i2,s5i2,c6i2,s6i2,c8i2,s8i2;

  REAL8 Lhxy;
  REAL8 alpha=0.;
  REAL8 alphaold;

  const UINT4 Npoints = 50;

  REAL8 t0,tAs;
  REAL8 om0,om1,omold,om;
  REAL8 Psi0,alpha0,iota0;
  REAL8 diota1,diota0,dalpha0,dalpha1;

  COMPLEX8Vector *modefreqs;
  COMPLEX16 MultSphHarmP;	// Generic spin-weighted spherical harmonics
  COMPLEX16 MultSphHarmM;	// Generic spin-weighted spherical harmonics
  REAL8 x0, x1, x2, x3;
  UINT4 nmodes,errcode;

  REAL8 finalMass,finalSpin;
  REAL8 energy=0.;
  REAL8 fracRD,omegaRD;

  static const char *func = "XLALPSpinInspiralEDEngine";

  ak = paramsInit->ak;

  mass = params->totalMass * LAL_MTSUN_SI;
  unitHz = params->totalMass * LAL_MTSUN_SI * (REAL8) LAL_PI;
  dt = 1.0 / params->tSampling;   /*    tSampling is in Hz, so dt is in seconds */

  if ((signalvec2)||(hh))
    params->nStartPad = 0;	/* must be zero for templates and injection */
  /* -- length in seconds from Newtonian formula; */
  lengths = (5.0 / 256.0) / LAL_PI * pow(LAL_PI * params->chirpMass * LAL_MTSUN_SI * params->fLower,-5.0 / 3.0) / params->fLower;

  /* setup coefficients for PN equations */
  XLALPSpinInspiralRDAdaptiveSetParams(&mparams,params,paramsInit);

  /* Check that initial frequency is smaller than omegamatch ~ xxyy for m=100 Msun */
  initphi=params->startPhase;
  initomega=params->fLower*unitHz;
  if ( initomega > (0.85 * omM0) ) {
    fprintf(stderr,"*** LALPSpinInspiralRD ERROR: Initial frequency too high: %11.5e for omM ~ %11.5e and m:(%8.3f, %8.3f)\n",params->fLower,omM0/unitHz,params->mass1,params->mass2);
    XLAL_ERROR( func , XLAL_EDOM );
  }
  v=pow(initomega,oneby3);
  v2 = v*v;

  /* Here we use the following convention:
     the coordinates of the spin vectors params->spin1,2 are assumed to
     be fed to this engine macro either in the frame set by the inital orbital 
     angular momentum (params->directionChoice= OrbitalL or TotalJ) or in the frame
     set by the viewing direction (params->directionChoice=View).
     The gw spherical modes are computed in whichever of the three frames
     specified by the params->directionChoice variable.
     The spin magnitude are normalized to the individual mass^2, i.e.
     they are dimension-less.
     The modulus of the initial angular momentum is fixed by m1,m2 and
     initial frequency, the inclination is the
     angle between the view direction and the orbital angular momentum.
     The polarization angle is not used here, it enters the pattern
     functions along with the angles marking the sky position of the
     source. */

  Lhmag = params->eta * params->totalMass * params->totalMass / v;

  for (i = 0; i < 3; i++) {
    initS1[i] = params->spin1[i] * params->mass1 * params->mass1;
    initS2[i] = params->spin2[i] * params->mass2 * params->mass2;
  }

  switch (params->axisChoice) {

  case OrbitalL:
    initLh[0] = 0.;
    initLh[1] = 0.;
    initLh[2] = 1.;
    inc = params->inclination;
    break;

  case View:
    initLh[0] = sin(params->inclination);
    initLh[1] = 0.;
    initLh[2] = cos(params->inclination);
    inc = 0.;
    break;

  default:
    //case TotalJ:
    for (j = 0; j < 3; j++) {
      iS1[j] = initS1[j];
      iS2[j] = initS2[j];
      initJ[j] = iS1[j] + iS2[j];
      initS1[j] = initS2[j]=0.;
      initLh[j] = 0.;
    }
    initJ[2] += Lhmag;
    initJmag = sqrt(initJ[0] * initJ[0] + initJ[1] * initJ[1] + initJ[2] * initJ[2]);
    if (initJ[0] == 0.)
      phiJ = 0.;
    else
      phiJ = atan2(initJ[1], initJ[0]);
    thetaJ = acos(initJ[2]/initJmag);
    rz[0][0] = cos(phiJ);
    rz[0][1] = sin(phiJ);
    rz[0][2] = 0.;
    rz[1][0] = -sin(phiJ);
    rz[1][1] = cos(phiJ);
    rz[1][2] = 0.;
    rz[2][0] = 0.;
    rz[2][1] = 0.;
    rz[2][2] = 1.;
    ry[0][0] = cos(thetaJ);
    ry[0][1] = 0;
    ry[0][2] = -sin(thetaJ);
    ry[1][0] = 0.;
    ry[1][1] = 1.;
    ry[1][2] = 0.;
    ry[2][0] = sin(thetaJ);
    ry[2][1] = 0.;
    ry[2][2] = cos(thetaJ);
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
	initLh[j] += ry[j][k] * rz[k][2];
	for (l = 0; l < 3; l++) {
	  initS1[j] += ry[j][k] * rz[k][l] * iS1[l];
	  initS2[j] += ry[j][k] * rz[k][l] * iS2[l];
	}
      }
    }
    inc = params->inclination;
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
  
  /* initialize the coordinates */
  yinit[0] = initphi;     /* phi */
  yinit[1] = initomega;   /* omega (really pi M f) */
  yinit[2] = initLh[0];   /* Lh(x,y,z) */
  yinit[3] = initLh[1];
  yinit[4] = initLh[2];
	
  yinit[5] = initS1[0];   /* Spin1(x,y,z) */
  yinit[6] = initS1[1];
  yinit[7] = initS1[2];
 
  yinit[8] = initS2[0];   /* Spin2(x,y,z) */
  yinit[9] = initS2[1]; 
  yinit[10]= initS2[2]; 

  yinit[11] = initS2[0];   /* dLh(x,y,z) */
  yinit[12] = initS2[1]; 
  yinit[13] = initS2[2]; 


  /* allocate the integrator */
  integrator = XLALAdaptiveRungeKutta4Init(14,XLALPSpinInspiralRDAdaptiveDerivatives,XLALPSpinInspiralRDAdaptiveTest,1.0e-6,1.0e-6);							
  if (!integrator) {
    fprintf(stderr,"LALPSpinInspiralRD: Cannot allocate integrator.\n");
    if (XLALClearErrno() == XLAL_ENOMEM)
      XLAL_ERROR( func, XLAL_ENOMEM );
    else
      XLAL_ERROR( func, XLAL_EDOM );
  }
  
  /* stop the integration only when the test is true */
  integrator->stopontestonly = 1;
		
  /* run the integration; note: time is measured in units of total mass */
  intlen = XLALAdaptiveRungeKutta4(integrator,(void *)&mparams,yinit,0.0,lengths/mass,dt/mass,&yout);
  
  intreturn = integrator->returncode;
  XLALAdaptiveRungeKutta4Free(integrator);
	
  if (!intlen) {
    if (XLALClearErrno() == XLAL_ENOMEM) {
      XLAL_ERROR( func,  XLAL_ENOMEM);
    } else {
      fprintf(stderr,"LALPSpinInspiralRD: integration failed with errorcode %d.\n",intreturn);
      XLAL_ERROR( func, XLAL_EFAILED);
    }
  }
  
  /* report on abnormal termination */ 
  if (intreturn != LALPSIRDPN_TEST_OMEGAMATCH)// && intreturn != LALPSIRDPN_TEST_ENERGY && intreturn != LALPSIRDPN_TEST_OMEGADOT) 
    {
      fprintf(stderr,"LALPSpinInspiralRD WARNING: integration terminated with code %d.\n",intreturn);
      fprintf(stderr," Waveform parameters were m1 = %14.6e, m2 = %14.6e, inc = %10.6f,\n", params->mass1, params->mass2, params->inclination);
      fprintf(stderr,"                          S1 = (%10.6f,%10.6f,%10.6f)\n", params->spin1[0], params->spin1[1], params->spin1[2]);
      fprintf(stderr,"                          S2 = (%10.6f,%10.6f,%10.6f)\n", params->spin2[0], params->spin2[1], params->spin2[2]);
    }
  
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
    XLAL_ERROR( func, XLAL_ENOMEM );
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
  
  /* The number of Ring Down modes is hard-coded here, it cannot exceed 3 */
  nmodes = 2;
  
  /* For RD, check that the 220 QNM freq. is less than the Nyquist freq. */
  /* Get QNM frequencies */
  modefreqs = XLALCreateCOMPLEX8Vector(nmodes);

  /* Call XLALFinalMassSpin() to get mass and spin of the final black hole */
  errcode = XLALPSpinFinalMassSpin(&finalMass, &finalSpin, params, energy, initLh);
  if (errcode != XLAL_SUCCESS) {
    XLALDestroyCOMPLEX8Vector(modefreqs);
    XLAL_ERROR(func,XLAL_EFAILED);
  }
  
  errcode = XLALPSpinGenerateQNMFreq(modefreqs, params, 2, 2, nmodes, finalMass, finalSpin);
  if (errcode != XLAL_SUCCESS) {
    XLALDestroyCOMPLEX8Vector(modefreqs);
    XLAL_ERROR(func,XLAL_EFAILED);
  }
  
  omegaRD = modefreqs->data[0].re * unitHz / LAL_PI / 2.;
  /* If Nyquist freq. <  220 QNM freq., one could print a warning message */
  /* Note that we cancelled a factor of 2 occuring on both sides */
  /* if (params->tSampling < modefreqs->data[0].re / LAL_PI)
     fprintf(stdout,"LALPhenSpin WARNING : Estimated ringdown freq larger than Nyquist freq.\n"); */
  
  params->ampOrder = 1;
  if (params->distance > 0.) 
    amp22ini = -2.0 * params->mu * LAL_MRSUN_SI / params->distance * sqrt(16. * LAL_PI / 5.);
  else
    amp22ini = 2. * sqrt(LAL_PI / 5.0) * params->signalAmplitude;

  /* if we have enough space, compute the waveform components; otherwise abort */
  if ((signalvec1 && intlen >= signalvec1->length) || (ff && intlen >= ff->length)) {
    if (signalvec1) {
      fprintf(stderr,"LALPSpinInspiralRD: no space to write in signalvec1: %d vs. %d\n",intlen,signalvec1->length);
    } else if (ff) {
      fprintf(stderr,"LALPSpinInspiralRD: no space to write in f: %d vs. %d\n",intlen,ff->length);
    } else {
      fprintf(stderr,"LALPSpinInspiralRD: no space to write anywhere!\n");
    }
    XLAL_ERROR(func, XLAL_ESIZE);
  }
		
    
  REAL8 *Phi    = &yout->data[1*intlen];
  REAL8 *omega  = &yout->data[2*intlen];
  REAL8 *Lhx    = &yout->data[3*intlen];
  REAL8 *Lhy    = &yout->data[4*intlen];
  REAL8 *Lhz    = &yout->data[5*intlen];
  REAL8 *S1x    = &yout->data[6*intlen];
  REAL8 *S1y    = &yout->data[7*intlen];
  REAL8 *S1z    = &yout->data[8*intlen];
  REAL8 *S2x    = &yout->data[9*intlen];
  REAL8 *S2y    = &yout->data[10*intlen];
  REAL8 *S2z    = &yout->data[11*intlen];


  for (j=0;j<intlen;j++) {
    
    v=pow(omega[j],oneby3);
    v2=v*v;

    // amp22= -2.0 * params->mu * LAL_MRSUN_SI/(params->distance) * sqrt( 16.*LAL_PI/5.)*v2;
    // amp20 = amp22*sqrt(3/2)
    // Y22 \pm Y2-2= sqrt(5/PI)    ((1+cos^2 t)/4, (cos t)/2)
    // Y21 \pm Y2-1= sqrt(5/PI)    ((sin t)/2, (sin 2t)/4)
    // Y20         = sqrt(15/2 PI) (sin^2 t)/4
    
    amp22 = amp22ini * v2;
    amp20 = amp22 * sqrt(1.5);
    amp33 = -amp22 / 4. * sqrt(5. / 42.);
    amp32 = amp33 * sqrt(1.5);
    amp31 = amp33 * sqrt(0.6)/2.;
    amp30 = amp33 / sqrt(5.) /2.;
    amp44 = amp22 * sqrt(5./7.) * 2./9.* v2;
    amp43 = - amp22 * sqrt(10./7.) / 36. * v2; 
    amp42 = - amp22 * sqrt(5.) / 2. / 63. * v2;
    amp41 = - amp22 * sqrt(2.5) / 126. * v2; 
    amp40 = - amp22 * sqrt(12.5) / 252. * v2;

    Psi = Phi[j] - 2. * omega[j] * log(omega[j]);
    
    ci = (Lhz[j]);
    cdi = 2.*ci*ci-1.;
    c2i = ci * ci;
    s2i = 1. - ci * ci;
    si = sqrt(s2i);
    sdi = 2.*ci*si;
    c2i2 = (1. + ci) / 2.;
    s2i2 = (1. - ci) / 2.;
    ci2 = sqrt(c2i2);
    si2 = sqrt(s2i2);
    c3i2=c2i2*ci2;
    s3i2=s2i2*si2;
    c4i2 = c2i2 * c2i2;
    s4i2 = s2i2 * s2i2;
    c5i2 = c4i2 * ci2;
    s5i2 = s4i2 * si2;
    c6i2 = c4i2 * c2i2;
    s6i2 = s4i2 * s2i2;
    c8i2 = c4i2 * c4i2;
    s8i2 = s4i2 * s4i2;
    
    alphaold = alpha;
    if ((Lhy[j]*Lhy[j]+Lhx[j]*Lhx[j])>0.)
      alpha = atan2(Lhy[j], Lhx[j]);
    else 
      alpha = alphaold;
    
    // In the even(odd) number entries of hlm the real(imaginary) value of 
    // hlm is stored as taken from AppendixB of arXiv:0810.5336. 
    // hl-m=(-1)^(l+m) hlm*(Psi+pi)
    // The components h20, h32, h31 and h30 are actually not in the paper,
    // but I(Riccardo) got them from Evan Ochsner (will they be ever 
    // published?).
    // After combination with the spherical harmonics Yxy, the polarizations 
    // h+ and hx are obtained,
    // according to 
    // h+ - ihx = Sum_{l,m} h_{lm}Y_{lm}. 
    // Note that by adding all the l=2 modes, eqs. (31-37) of 
    // arXiv:gr-qc/0211087 are obtained,  **BUT** for an overall minus sign 
    // in both h+ and hx.
    
    h2P2->data[2 * j] = amp22 * ( 1. /(1. + v2 / 42. * (107. - 55. * mparams.eta) )* 
				  ( cos(2. * (Psi + alpha)) * c4i2 + cos(2. * (Psi - alpha)) * s4i2) +
				  v * mparams.dm / 3. * si *
				  ( cos(Psi - 2. * alpha) * s2i2 + cos(Psi + 2. * alpha) * c2i2 ) );
    
    h2M2->data[2 * j] = amp22 * ( 1. /(1. + v2 / 42. * (107. - 55. * mparams.eta) )* 
				  ( cos(2. * (Psi + alpha)) * c4i2 + cos(2. * (Psi - alpha)) * s4i2) -
				  v * mparams.dm / 3. * si * 
				  ( cos(Psi - 2. * alpha) * s2i2 + cos(Psi + 2. * alpha) * c2i2) );
    
    h2P2->data[2 * j + 1] = amp22 * (1. /(1. + v2 / 42. * (107. - 55. * mparams.eta) )* 
				     ( -sin(2. * (Psi + alpha)) * c4i2 + sin(2. * (Psi - alpha)) * s4i2) +
				     v * mparams.dm / 3. * si *
				     ( sin(Psi - 2. * alpha) * s2i2 - sin(Psi + 2. * alpha) * c2i2) );
    
    h2M2->data[2 * j + 1] = amp22 * (1. /(1. + v2 / 42. * (107. - 55. * mparams.eta) )* 
				     ( sin(2. * (Psi + alpha)) * c4i2 - sin(2. * (Psi - alpha)) * s4i2) +
				     v * mparams.dm / 3. * si *
				     ( sin(Psi - 2. * alpha) * s2i2 - sin(Psi + 2. * alpha) * c2i2) );
    
    h2P1->data[2 * j] = amp22 * (si /(1. + v2 / 42. * (107. - 55. * mparams.eta) ) *
				 ( -cos(2. * Psi - alpha) * s2i2 + cos(2. * Psi + alpha) * c2i2) +
				 v * mparams.dm / 3. *
				 (-cos(Psi + alpha) * (ci + cdi)/2. - 
				  cos(Psi - alpha) * s2i2 * (1. + 2. * ci) ) );
    
    h2M1->data[2 * j] = amp22 * (si / (1. + v2 / 42. * (107. - 55. * mparams.eta)) *
				 ( cos(2. * Psi - alpha) * s2i2 - cos(2. * Psi + alpha) * c2i2) +
				 v * mparams.dm / 3. * 
				 (-cos(Psi + alpha) * (ci + cdi)/2. -
				  cos(Psi - alpha) * s2i2 * (1. + 2. * ci) ) );
    
    h2P1->data[2 * j + 1] = amp22 * (si /(1. + v2 / 42. * (107. - 55. * mparams.eta) ) * 
				     ( -sin(2. * Psi - alpha ) * s2i2 - sin(2. * Psi + alpha) * c2i2) +
				     v * mparams.dm / 3. * 
				     (sin(Psi + alpha) * (ci + cdi)/2. - 
				      sin(Psi - alpha) * s2i2 * (1. + 2. * ci) ) );
    
    
    h2M1->data[2 * j + 1] = amp22 * (si / (1. + v2 / 42. * (107. - 55. * mparams.eta)) * 
				     ( -sin(2. * Psi - alpha) * s2i2 - sin(2. * Psi + alpha) * c2i2) -
				     v * mparams.dm / 3. * 
				     (sin(Psi + alpha) * (ci + cdi) / 2. - 
				      sin(Psi - alpha) * s2i2 * (1. + 2. * ci) ) );
    
    h20->data[2 * j] = amp20 * ( s2i / (1.+ v2/42. * (107. - 55.*mparams.eta) ) * cos(2. * Psi) );
    
    h20->data[2 * j + 1] = amp20 * ( v * mparams.dm / 3. * sdi * sin(Psi) );
    
    h3P3->data[2 * j] = amp33 * ( v * mparams.dm *
				  (-9. * cos(3. * (Psi - alpha)) * s6i2 -
				   cos(  Psi - 3. * alpha) * s4i2 * c2i2 +
				   cos(  Psi + 3. * alpha) * s2i2 * c4i2 +
				   9. * cos(3. * (Psi + alpha)) * c6i2) +
				  v2 * 4. * si * (1. - 3. * mparams.eta) * 
				  (    - cos(2. * Psi - 3. * alpha) * s4i2 +
				       cos(2. * Psi + 3. * alpha) * c4i2 ) );
    
    h3M3->data[2 * j] = amp33 * (-v * mparams.dm *
				 (-9. * cos(3. * (Psi - alpha)) * s6i2 -
				  cos(  Psi - 3. * alpha) * s4i2 * c2i2 +
				  cos(  Psi + 3. * alpha) * s2i2 * c4i2 +
				  9. * cos(3. * (Psi + alpha)) * c6i2) +
				 v2 * 4. * si * (1. - 3. * mparams.eta) * 
				 (    - cos(2. * Psi - 3. * alpha) * s4i2 +
				      cos(2. * Psi + 3. * alpha) * c4i2 ) );
    
    h3P3->data[2 * j + 1] = amp33 * ( v * mparams.dm *
				      (-9. * sin(3. * (Psi - alpha)) * s6i2 -
				       sin(  Psi - 3. * alpha) * s4i2 * c2i2 -
				       sin(  Psi + 3. * alpha) * s2i2 * c4i2 -
				       9. * sin(3. * (Psi + alpha)) * c6i2) +
				      v2 * 4. * si * (1. - 3. * mparams.eta) * 
				      (  - sin(2. * Psi - 3. * alpha) * s4i2
					 - sin(2. * Psi + 3. * alpha) * c4i2 ) );
    
    h3M3->data[2 * j + 1] = amp33 * ( v * mparams.dm *
				      (-9. * sin(3. * (Psi - alpha)) * s6i2 -
				       sin(  Psi - 3. * alpha) * s4i2 * c2i2 -
				       sin(  Psi + 3. * alpha) * s2i2 * c4i2 -
				       9. * sin(3. * (Psi + alpha)) * c6i2) -
				      v2 * 4. * si * (1. - 3. * mparams.eta) * 
				      (   - sin(2. * Psi - 3. * alpha) * s5i2 * ci2
					  - sin(2. * Psi + 3. * alpha) * c5i2 * si2 ) );

    h3P2->data[2 * j] = amp32 * ( v * mparams.dm / 3. * 
				  ( 27. * cos(3. * Psi - 2. * alpha) * si*s4i2 + 
				    27. * cos(3. * Psi + 2. * alpha) * si*c4i2 +
				    cos( Psi + 2. * alpha) * c3i2 * (5.*si2-3.*si*ci2-3.*ci*si2) /2. +
				    cos( Psi - 2. * alpha) * s3i2 * (5.*ci2+3.*ci*ci2-3.*si*si2) /2.) +
				  v2*(1./3.-mparams.eta) * 
				  ( - 8.*c4i2*(3.*ci-2.)*cos(2.*(Psi+alpha)) +
				    8.*s4i2*(3.*ci+2.)*cos(2.*(Psi-alpha)) ) );

    h3M2->data[2 * j] = amp32 * ( v * mparams.dm / 3. *
				  ( 27. * cos(3. * Psi - 2. * alpha) * si*s4i2 +
				    27. * cos(3. * Psi + 2. * alpha) * si*c4i2 +
				    cos( Psi + 2. * alpha) * c3i2 * (5.*si2-3.*si*ci2-3.*ci*si2) /2. +
				    cos( Psi - 2. * alpha) * s3i2 * (5.*ci2+3.*ci*ci2-3.*si*si2) /2.) +
				  v2*(1./3.-mparams.eta) * 
				  ( 8.*c4i2*(3.*ci-2.)*cos(2.*(Psi+alpha)) -
				    8.*s4i2*(3.*ci+2.)*cos(2.*(Psi-alpha)) ) );

    h3P2->data[2 * j + 1 ] = amp33 * ( v * mparams.dm / 3. * 
				       ( 27. * sin(3. * Psi - 2. * alpha) * si*s4i2 - 
					 27. * cos(3. * Psi + 2. * alpha) * si*c4i2 -
					 sin( Psi + 2. * alpha) * c3i2 * (5.*si2-3.*si*ci2-3.*ci*si2) /2. +
					 sin( Psi - 2. * alpha) * s3i2 * (5.*ci2+3.*ci*ci2-3.*si*si2) /2.) +
				       v2*(1./3.-mparams.eta) *
				       ( 8.*c4i2*(3.*ci-2.)*sin(2.*(Psi+alpha)) +
					 8.*s4i2*(3.*ci+2.)*sin(2.*(Psi-alpha)) ) );

    
    h3P2->data[2 * j + 1 ] = amp32 * ( -v * mparams.dm / 3. * 
				       ( 27. * sin(3. * Psi - 2. * alpha) * si*s4i2 - 
					 27. * cos(3. * Psi + 2. * alpha) * si*c4i2 -
					 sin( Psi + 2. * alpha) * c3i2 * (5.*si2-3.*si*ci2-3.*ci*si2) /2.+
					 sin( Psi - 2. * alpha) * s3i2 * (5.*ci2+3.*ci*ci2-3.*si*si2) /2.)+
				       v2*(1./3.-mparams.eta) *
				       ( 8.*c4i2*(3.*ci-2.)*sin(2.*(Psi+alpha)) +
					 8.*s4i2*(3.*ci+2.)*sin(2.*(Psi-alpha)) ) );

    h3P1->data[2 * j] = amp31 * ( v * mparams.dm / 6. * 
				  ( -135. * cos(3.*Psi - alpha) * s2i*s2i2 +
				    135. * cos(3.*Psi + alpha)  * s2i*c2i2 +
				    cos(Psi+alpha) * c2i2*(15.*cdi-20.*ci+13.)/2.-
				    cos(Psi-alpha) * s2i2*(15.*cdi+20.*ci+13.)/2. )
				  -v2 * (1./3.-mparams.eta)* 
				  ( 20.*c3i2*cos(2.*Psi+alpha)*(3.*(si2*ci+ci2*si)-5.*si2) +
				    20.*s3i2*cos(2.*Psi-alpha)*(3.*(ci2*ci-si2*si)+5.*ci2) ) );

    h3M1->data[2 * j] = amp31 * (-v * mparams.dm / 6. *
				 ( -135. * cos(3.*Psi - alpha) * s2i*s2i2 +
				   135. * cos(3.*Psi + alpha) * s2i*c2i2 +
				   cos(Psi+alpha) * c2i2*(15.*cdi-20.*ci+13.)/2.-
				   cos(Psi-alpha) * s2i2*(15.*cdi+20.*ci+13.)/2. )
				 -v2 * (1./3.-mparams.eta)*
				 ( 20.*c3i2*cos(2.*Psi+alpha)*(3.*(si2*ci+ci2*si)-5.*si2) +
				   20.*s3i2*cos(2.*Psi-alpha)*(3.*(ci2*ci-si2*si)+5.*ci2) ) );

    h3P1->data[2 * j + 1] = amp31 * ( v * mparams.dm / 6. *
				      ( -135. * sin(3.*Psi - alpha) * s2i*s2i2 -
					135.* sin(3.*Psi + alpha) * s2i*c2i2 -
					sin(Psi+alpha) * c2i2*(15.*cdi-20.*ci+13.)/2.-
					sin(Psi-alpha) * s2i2*(15.*cdi+20.*ci+13.)/2. )
				      +v2 * (1./3.-mparams.eta)*
				      ( 20.*c3i2*sin(2.*Psi+alpha)*(3.*(si2*ci+ci2*si)-5.*si2)
					-20.*s3i2*sin(2.*Psi-alpha)*(3.*(ci2*ci-si2*si)+5.*ci2) ) );

    h3M1->data[2 * j + 1] = amp31 * ( v * mparams.dm / 6. *
				      ( -135. * sin(3.*Psi - alpha) *s2i*s2i2 -
					135. * sin(3.*Psi + alpha) *s2i*c2i2 -
					sin(Psi+alpha) * c2i2*(15.*cdi-20.*ci+13.)/2.-
					sin(Psi-alpha) * s2i2*(15.*cdi+20.*ci+13.)/2. )
				      -v2 * (1./3.-mparams.eta)* 
				      ( 20.*c3i2*sin(2.*Psi+alpha)*(3.*(si2*ci+ci2*si)-5.*si2)
					-20.*s3i2*sin(2.*Psi-alpha)*(3.*(ci2*ci-si2*si)+5.*ci2) ) );

    h30->data[2 * j] = 0.;
    
    h30->data[2 * j + 1] = amp30 * ( v * mparams.dm *
				     ( cos(Psi) * si*(cos(2.*Psi)*(45.*s2i)-(25.*cdi-21.) ) )
				     +v2 * (1.-3.*mparams.eta) *
				     (80. * s2i*c2i*sin(2.*Psi) ) );

    h4P4->data[2 * j] = amp44 * (1. - 3.*mparams.eta) * 
                                  ( 4.* s8i2 * cos(4.*(Psi-alpha)) + cos(2.*Psi-4.*alpha) *s6i2*c2i2 
				    + s2i2*c6i2* cos(2.*Psi+4.*alpha) + 4.*c8i2* cos(4.*(Psi+alpha)) );

    h4M4->data[2 * j] = amp44 * (1. - 3.*mparams.eta) * 
                                  ( 4.* s8i2 * cos(4.*(Psi-alpha)) + cos(2.*Psi-4.*alpha) *s6i2*c2i2 
				    + s2i2*c6i2* cos(2.*Psi+4.*alpha) + 4.*c8i2* cos(4.*(Psi+alpha)) );

    h4P4->data[2 * j + 1] = amp44 * (1. - 3.*mparams.eta) * 
                                  ( 4.* s8i2 * sin(4.*(Psi-alpha)) + sin(2.*Psi-4.*alpha) *s6i2*c2i2 
				    - s2i2*c6i2* sin(2.*Psi+4.*alpha) - 4.*c8i2* sin(4.*(Psi+alpha)) );

    h4M4->data[2 * j + 1] = - amp44 * (1. - 3.*mparams.eta) * 
                                  ( 4.* s8i2 * sin(4.*(Psi-alpha)) + sin(2*Psi-4.*alpha) *s6i2*c2i2 
				    - s2i2*c6i2* sin(2.*Psi+4.*alpha) - 4.*c8i2* sin(4.*(Psi+alpha)) );

    h4P3->data[2 * j] = amp43 * (1. - 3.*mparams.eta) * 
      ( 4.*s6i2* cos(4.*Psi-3.*alpha) - 4.*c6i2* cos(4.*Psi+3.*alpha) - 
	s4i2*(ci+0.5)/2. * cos(2.*Psi-3.*alpha) - c4i2*(ci-0.5) * cos(2.*Psi+3.*alpha) );

    h4M3->data[2 * j] = - amp43 * (1. - 3.*mparams.eta) * 
      ( 4.*s6i2* cos(4.*Psi-3.*alpha) - 4.*c6i2* cos(4.*Psi+3.*alpha) - 
	s4i2*(ci+0.5)/2. * cos(2.*Psi-3.*alpha) - c4i2*(ci-0.5) * cos(2.*Psi+3.*alpha) );

    h4P3->data[2 * j + 1] = amp43 * (1. - 3.*mparams.eta) * 
      ( 4.*s6i2* sin(4.*Psi-3.*alpha) + 4.*c6i2* sin(4.*Psi+3.*alpha) - 
	s4i2*(ci+0.5)/2. * sin(2.*Psi-3.*alpha) + c4i2*(ci-0.5) * sin(2.*Psi+3.*alpha) );

    h4M3->data[2 * j + 1] = amp43 * (1. - 3.*mparams.eta) * 
      ( 4.*s6i2* sin(4.*Psi-3.*alpha) + 4.*c6i2* sin(4.*Psi+3.*alpha) - 
	s4i2*(ci+0.5)/2. * sin(2.*Psi-3.*alpha) + c4i2*(ci-0.5) * sin(2.*Psi+3.*alpha) );

    h4P2->data[2 * j] = amp42 * (1. - 3.*mparams.eta) *
      ( 224.*s6i2*c2i2 * cos(4.*Psi-2.*alpha) - 224.*c6i2*s2i2 * cos(4.*Psi+2.*alpha) 
	- c4i2 * cos(2.*(Psi+alpha))*(7.*cdi-14.*ci+9) - s4i2 * cos(2.*(Psi-alpha))*(7.*cdi+14*ci+9) );

    h4M2->data[2 * j] = amp42 * (1. - 3.*mparams.eta) *
      ( 224.*s6i2*c2i2 * cos(4.*Psi-2.*alpha) - 224.*c6i2*s2i2 * cos(4.*Psi+2.*alpha) 
	- c4i2 * cos(2.*(Psi+alpha))*(7.*cdi-14.*ci+9) - s4i2 * cos(2.*(Psi-alpha))*(7.*cdi+14*ci+9) );

    h4P2->data[2 * j + 1] = amp42 * (1. - 3.*mparams.eta) *
      ( 224.*s6i2*c2i2 * sin(4.*Psi-2.*alpha) + 224.*c6i2*s2i2 * sin(4.*Psi+2.*alpha) 
	+ c4i2 * sin(2.*(Psi+alpha))*(7.*cdi-14.*ci+9) - s4i2 * sin(2.*(Psi-alpha))*(7.*cdi+14*ci+9) );

    h4M2->data[2 * j + 1] = -amp42 * (1. - 3.*mparams.eta) *
      ( 224.*s6i2*c2i2 * sin(4.*Psi-2.*alpha) + 224.*c6i2*s2i2 * sin(4.*Psi+2.*alpha) 
	+ c4i2 * sin(2.*(Psi+alpha))*(7.*cdi-14.*ci+9) - s4i2 * sin(2.*(Psi-alpha))*(7.*cdi+14*ci+9) );
    
    h4P1->data[2 * j] = amp41 * (1. - 3.*mparams.eta) * 
      ( 448.*s5i2*c3i2 * cos(4.*Psi-alpha) - 448.*s3i2*c5i2 * cos(4.*Psi+alpha) +
	s3i2*cos(2.*Psi-alpha)*(7.*(cdi*ci2-sdi*si2)+14.*(ci2*ci-si2*si)+19.*ci2) -
	c3i2*cos(2.*Psi+alpha)*(7.*(cdi*si2+sdi*ci2)-14.*(si*ci2+ci*si2)+19.*ci2) );

    h4M1->data[2 * j] = -amp41 * (1. - 3.*mparams.eta) * 
      ( 448.*s5i2*c3i2 * cos(4.*Psi-alpha) - 448.*s3i2*c5i2 * cos(4.*Psi+alpha) +
	s3i2*cos(2.*Psi-alpha)*(7.*(cdi*ci2-sdi*si2)+14.*(ci2*ci-si2*si)+19.*ci2) -
	c3i2*cos(2.*Psi+alpha)*(7.*(cdi*si2+sdi*ci2)-14.*(si*ci2+ci*si2)+19.*ci2) );

    h4P1->data[2 * j + 1] = amp41 * (1. - 3.*mparams.eta) * 
      ( 448.*s5i2*c3i2 * sin(4.*Psi-alpha) + 448.*s3i2*c5i2 * sin(4.*Psi+alpha) +
	s3i2*sin(2.*Psi-alpha)*(7.*(cdi*ci2-sdi*si2)+14.*(ci2*ci-si2*si)+19.*ci2) +
	c3i2*sin(2.*Psi+alpha)*(7.*(cdi*si2+sdi*ci2)-14.*(si*ci2+ci*si2)+19.*ci2) );

    h4M1->data[2 * j + 1] = amp41 * (1. - 3.*mparams.eta) * 
      ( 448.*s5i2*c3i2 * sin(4.*Psi-alpha) + 448.*s3i2*c5i2 * sin(4.*Psi+alpha) +
	s3i2*sin(2.*Psi-alpha)*(7.*(cdi*ci2-sdi*si2)+14.*(ci2*ci-si2*si)+19.*ci2) +
	c3i2*sin(2.*Psi+alpha)*(7.*(cdi*si2+sdi*ci2)-14.*(si*ci2+ci*si2)+19.*ci2) );

    h40->data[2 * j] = amp40 * (1.-3.*mparams.eta) * s2i * (-56.*cos(4.*Psi) - 
							       cos(2.*Psi)*(7.*cdi+5.) );
    h40->data[2 * j +1] = 0.;
    
    fap->data[j] = omega[j];
    phap->data[j] = Phi[j];
 
  }

  REAL8 sqrtOneMinus4Eta;
  REAL8 omegamatch;

  const REAL8 omM     = 0.055;
  const REAL8 omMz1p2 = -5.02e-3;
  const REAL8 omM12   = -4.29e-4;
  const REAL8 omMsq   = 5.78e-3;
  const REAL8 omMz12  = 2.66e-3;
  const REAL8 omMzsq  = -9.27e-3;

  if (params->eta<0.25) sqrtOneMinus4Eta=sqrt(1.-4.*params->eta);
  else                  sqrtOneMinus4Eta=0.;

  REAL8Vector *omega_s = XLALCreateREAL8Vector(Npoints);
  REAL8Vector *Lhx_s   = XLALCreateREAL8Vector(Npoints);
  REAL8Vector *Lhy_s   = XLALCreateREAL8Vector(Npoints);
  REAL8Vector *Lhz_s   = XLALCreateREAL8Vector(Npoints);

  REAL8Vector *domega_s = XLALCreateREAL8Vector(Npoints);
  REAL8Vector *dLhx_s   = XLALCreateREAL8Vector(Npoints);
  REAL8Vector *dLhy_s   = XLALCreateREAL8Vector(Npoints);
  REAL8Vector *dLhz_s   = XLALCreateREAL8Vector(Npoints);
  REAL8Vector *diota    = XLALCreateREAL8Vector(Npoints);
  REAL8Vector *dalpha   = XLALCreateREAL8Vector(Npoints);

  REAL8Vector *ddomega_s = XLALCreateREAL8Vector(Npoints);
  REAL8Vector *ddiota    = XLALCreateREAL8Vector(Npoints);
  REAL8Vector *ddalpha   = XLALCreateREAL8Vector(Npoints);

  REAL8 LhS1,LhS2,S1S2,S1S1,S2S2;

  for (j=0; j<Npoints; j++) {
    k=j+intlen-Npoints;
    omega_s->data[j] = omega[k];
    Lhx_s->data[j]   = Lhx[k];
    Lhy_s->data[j]   = Lhy[k];
    Lhz_s->data[j]   = Lhz[k];

    LhS1=(Lhx[k]*S1x[k]+Lhy[k]*S1y[k]+Lhz[k]*S1z[k])/mparams.m1m/mparams.m1m;
    LhS2=(Lhx[k]*S2x[k]+Lhy[k]*S2y[k]+Lhz[k]*S2z[k])/mparams.m2m/mparams.m2m;;
    S1S2=(S1x[k]*S2x[k]+S1y[k]*S2y[k]+S1z[k]*S2z[k])/pow(mparams.m1m*mparams.m2m,2.);
    S1S1=(S1x[k]*S1x[k]+S1y[k]*S1y[k]+S1z[k]*S1z[k])/pow(mparams.m1m,4.);
    S2S2=(S2x[k]*S2x[k]+S2y[k]*S2y[k]+S2z[k]*S2z[k])/pow(mparams.m2m,4.);

    omegamatch= omM + 6.05e-3 * sqrtOneMinus4Eta + 
      omMz1p2 * (LhS1 + LhS2) + omM12 * (S1S2 - LhS1 * LhS2) + 
      omMsq   * (S1S1 + S2S2 - LhS1 * LhS1 - LhS2 * LhS2) +
      omMz12  * (LhS1 * LhS2) + omMzsq * (LhS1 * LhS1 + LhS2 * LhS2);

    if ((omegamatch<omega_s->data[j])&&(jend==0)) {
      jend=j;
      kend=k;
    }
  }

  printf("omM %12.6e  %12.6e  %d\n",omegamatch,omega_s->data[jend],jend);

  if (jend==0) fprintf(stdout,"***** LALPSIRD ERROR *****: jend not set\n");

  errcode  = XLALGenerateWaveDerivative(domega_s,omega_s,dt);
  errcode += XLALGenerateWaveDerivative(dLhx_s,Lhx_s,dt);
  errcode += XLALGenerateWaveDerivative(dLhy_s,Lhy_s,dt);
  errcode += XLALGenerateWaveDerivative(dLhz_s,Lhz_s,dt);
  if (errcode != 0) {
    fprintf(stderr,"*** LALPSpinInspiralRD: ERROR generating derivative ***\n");
    XLAL_ERROR(func,XLAL_EFAILED);
  }

  for (j=0; j<Npoints; j++) {
    Lhxy = sqrt(Lhx_s->data[j] * Lhx_s->data[j] + Lhy_s->data[j] * Lhy_s->data[j]);
    if (Lhxy > 0.) {
      diota->data[j]  = -dLhz_s->data[j] / Lhxy;
      dalpha->data[j] = (Lhx_s->data[j] * dLhy_s->data[j] - Lhy_s->data[j] * dLhx_s->data[j]) / Lhxy;
    } else {
      diota->data[j]  = 0.;
      dalpha->data[j] = 0.;
    }
  }

  errcode  = XLALGenerateWaveDerivative(ddiota,diota,dt);
  errcode += XLALGenerateWaveDerivative(ddalpha,dalpha,dt);
  errcode += XLALGenerateWaveDerivative(ddomega_s,domega_s,dt);
  if (errcode != 0) {
    fprintf(stderr,"*** LALPSpinInspiralRD: ERROR generating derivative ***\n");
    XLAL_ERROR(func,XLAL_EFAILED);
  }
  
  tim = t0 = yout->data[kend] * mass;
  tAs = t0 + 2. * domega_s->data[jend] / ddomega_s->data[jend];
  om1 = domega_s->data[jend] * tAs * (1. - t0 / tAs) * (1. - t0 / tAs);
  om0 = omega[kend] - om1 / (1. - t0 / tAs);
  om  = omega[kend];
  
  printf("t0 %12.6e tAs %12.6e dt %12.6e  m %12.6e  iota %12.6e  %12.6e\n",t0,tAs,dt,mass,ddiota->data[jend],diota->data[jend]);

  diota1 = ddiota->data[jend] * tAs * (1. - t0 / tAs) * (1. - t0 / tAs);
  diota0 = diota->data[jend] - diota1 / (1. - t0 / tAs);
  
  dalpha1 = ddalpha->data[jend] * tAs * (1. - t0 / tAs) * (1. - t0 / tAs);
  dalpha0 = dalpha->data[jend] - dalpha1 / (1. - t0 / tAs);

  if ((tAs < t0) || (om1 < 0.)) {
    fprintf(stderr,"** LALPSpinInspiralRD ERROR **: Could not attach phen part\n");
    XLAL_ERROR(func, XLAL_EFAILED);
  }
  
  printf("***** kend %d  jend %d  %d %d\n",kend,jend,intlen,Npoints);

  Psi0 = Psi + tAs * om1 * log(1. - t0 / tAs);
  alpha0 = alpha + tAs * dalpha1 * log(1. - t0 / tAs);
  iota0 = acos(Lhz[intlen-1]) + diota1 * tAs * log(1. - t0 / tAs);
  
  /* Get QNM frequencies */
  modefreqs = XLALCreateCOMPLEX8Vector(nmodes);

  errcode = XLALPSpinFinalMassSpin(&finalMass, &finalSpin, params, energy, initLh);
  if (errcode != XLAL_SUCCESS) {
    XLAL_ERROR(func,XLAL_EFAILED);
  }

  errcode =XLALPSpinGenerateQNMFreq(modefreqs, params, 2, 2, nmodes, finalMass, finalSpin);
  if (errcode != XLAL_SUCCESS) {
    XLALDestroyCOMPLEX8Vector(modefreqs);
    XLAL_ERROR(func,XLAL_EFAILED);
  }

  omegaRD = modefreqs->data[0].re * unitHz / LAL_PI / 2.;
  XLALDestroyCOMPLEX8Vector(modefreqs);

  const double frac0 = 0.57;
  const double frac1p2 = 1.42e-2;
  const double frac12 = -3.71e-3;
  const double fracsq = -1.201e-2;
  const double fracz12 = -2.447e-2;
  const double fraczsq = -1.930e-2;

  fracRD = frac0 + frac1p2 * (LhS1 + LhS2) + frac12 * (S1S2 - LhS1 * LhS1) +
    fracsq * (S1S1 + S2S2 - LhS1 * LhS1 - LhS2 * LhS2) +
    fracz12 * (LhS1 * LhS2) + fraczsq * (LhS1 * LhS1 + LhS2 * LhS2);

  j = kend+1;

  do {

    if (j >= length) {
      fprintf(stderr,"** LALPhenSpinInspiralRD ERROR**: phen part exceeds array length");
      fprintf(stderr, "** m (%11.4e  %11.4e)  f0 %11.4e\n",params->mass1, params->mass2, params->fLower);
      fprintf(stderr, "** S1 (%8.4f  %8.4f  %8.4f)\n", initS1[0],initS1[1], initS1[2]);
      fprintf(stderr, "** S2 (%8.4f  %8.4f  %8.4f)\n", initS2[0],initS2[1], initS2[2]);
      XLAL_ERROR(func,XLAL_ENOMEM);
    }

    tim += dt;
    omold = om;
    om = om1 / (1. - tim / tAs) + om0;
    if (j%10==0) printf("** P ** om[%d] %12.6e\n",j,om);
    fap->data[j] = om;
    Psi  = Psi0 - tAs * om1 * log(1. - tim / tAs) + om0 * (tim - t0);
    ci = cos(diota0 * (tim - t0) - diota1 * tAs * log(1. - tim / tAs) + iota0);
    alpha = alpha0 + dalpha0 * (tim - t0) - dalpha1 * tAs / mass * log(1. - tim / tAs);
    v2old = v2;
    v  = pow(om, oneby3);
    v2 = v*v;
    amp22 *= v2 / v2old;

    amp20 = amp22 * sqrt(1.5);
    amp33 = -amp22 / 4. * sqrt(5. / 42.);
    amp32 = amp33 * sqrt(1.5);
    amp31 = amp33 * sqrt(0.6)/2.;
    amp30 = amp33 / sqrt(5.) /2.;
    amp44 = amp22 * sqrt(5./7.) * 2./9.* v2;
    amp43 = - amp22 * sqrt(10./7.) / 36. * v2; 
    amp42 = - amp22 * sqrt(5.) / 2. / 63. * v2;
    amp41 = - amp22 * sqrt(2.5) / 126. * v2; 
    amp40 = - amp22 * sqrt(12.5) / 252. * v2;
    
    ci = (Lhz[j]);
    cdi = 2.*ci*ci-1.;
    c2i = ci * ci;
    s2i = 1. - ci * ci;
    si = sqrt(s2i);
    sdi = 2.*ci*si;
    c2i2 = (1. + ci) / 2.;
    s2i2 = (1. - ci) / 2.;
    ci2 = sqrt(c2i2);
    si2 = sqrt(s2i2);
    c3i2=c2i2*ci2;
    s3i2=s2i2*si2;
    c4i2 = c2i2 * c2i2;
    s4i2 = s2i2 * s2i2;
    c5i2 = c4i2 * ci2;
    s5i2 = s4i2 * si2;
    c6i2 = c4i2 * c2i2;
    s6i2 = s4i2 * s2i2;
    c8i2 = c4i2 * c4i2;
    s8i2 = s4i2 * s4i2;
    
    h2P2->data[2 * j] = amp22 * ( 1. /(1. + v2 / 42. * (107. - 55. * mparams.eta) )* 
				  ( cos(2. * (Psi + alpha)) * c4i2 + cos(2. * (Psi - alpha)) * s4i2) +
				  v * mparams.dm / 3. * si *
				  ( cos(Psi - 2. * alpha) * s2i2 + cos(Psi + 2. * alpha) * c2i2 ) );
    
    h2M2->data[2 * j] = amp22 * ( 1. /(1. + v2 / 42. * (107. - 55. * mparams.eta) )* 
				  ( cos(2. * (Psi + alpha)) * c4i2 + cos(2. * (Psi - alpha)) * s4i2) -
				  v * mparams.dm / 3. * si * 
				  ( cos(Psi - 2. * alpha) * s2i2 + cos(Psi + 2. * alpha) * c2i2) );
    
    h2P2->data[2 * j + 1] = amp22 * (1. /(1. + v2 / 42. * (107. - 55. * mparams.eta) )* 
				     ( -sin(2. * (Psi + alpha)) * c4i2 + sin(2. * (Psi - alpha)) * s4i2) +
				     v * mparams.dm / 3. * si *
				     ( sin(Psi - 2. * alpha) * s2i2 - sin(Psi + 2. * alpha) * c2i2) );
    
    h2M2->data[2 * j + 1] = amp22 * (1. /(1. + v2 / 42. * (107. - 55. * mparams.eta) )* 
				     ( sin(2. * (Psi + alpha)) * c4i2 - sin(2. * (Psi - alpha)) * s4i2) +
				     v * mparams.dm / 3. * si *
				     ( sin(Psi - 2. * alpha) * s2i2 - sin(Psi + 2. * alpha) * c2i2) );
    
    h2P1->data[2 * j] = amp22 * (si /(1. + v2 / 42. * (107. - 55. * mparams.eta) ) *
				 ( -cos(2. * Psi - alpha) * s2i2 + cos(2. * Psi + alpha) * c2i2) +
				 v * mparams.dm / 3. *
				 (-cos(Psi + alpha) * (ci + cdi)/2. - 
				  cos(Psi - alpha) * s2i2 * (1. + 2. * ci) ) );
    
    h2M1->data[2 * j] = amp22 * (si / (1. + v2 / 42. * (107. - 55. * mparams.eta)) *
				 ( cos(2. * Psi - alpha) * s2i2 - cos(2. * Psi + alpha) * c2i2) +
				 v * mparams.dm / 3. * 
				 (-cos(Psi + alpha) * (ci + cdi)/2. -
				  cos(Psi - alpha) * s2i2 * (1. + 2. * ci) ) );
    
    h2P1->data[2 * j + 1] = amp22 * (si /(1. + v2 / 42. * (107. - 55. * mparams.eta) ) * 
				     ( -sin(2. * Psi - alpha ) * s2i2 - sin(2. * Psi + alpha) * c2i2) +
				     v * mparams.dm / 3. * 
				     (sin(Psi + alpha) * (ci + cdi)/2. - 
				      sin(Psi - alpha) * s2i2 * (1. + 2. * ci) ) );
    
    
    h2M1->data[2 * j + 1] = amp22 * (si / (1. + v2 / 42. * (107. - 55. * mparams.eta)) * 
				     ( -sin(2. * Psi - alpha) * s2i2 - sin(2. * Psi + alpha) * c2i2) -
				     v * mparams.dm / 3. * 
				     (sin(Psi + alpha) * (ci + cdi) / 2. - 
				      sin(Psi - alpha) * s2i2 * (1. + 2. * ci) ) );
    
    h20->data[2 * j] = amp20 * ( s2i / (1.+ v2/42. * (107. - 55.*mparams.eta) ) * cos(2. * Psi) );
    
    h20->data[2 * j + 1] = amp20 * ( v * mparams.dm * 2. / 3. * si*ci * sin(Psi) );
    
    h3P3->data[2 * j] = amp33 * ( v * mparams.dm *
				  (-9. * cos(3. * (Psi - alpha)) * s6i2 -
				   cos(  Psi - 3. * alpha) * s4i2 * c2i2 +
				   cos(  Psi + 3. * alpha) * s2i2 * c4i2 +
				   9. * cos(3. * (Psi + alpha)) * c6i2) +
				  v2 * 4. * si * (1. - 3. * mparams.eta) * 
				  (    - cos(2. * Psi - 3. * alpha) * s4i2 +
				       cos(2. * Psi + 3. * alpha) * c4i2 ) );
    
    h3M3->data[2 * j] = amp33 * (-v * mparams.dm *
				 (-9. * cos(3. * (Psi - alpha)) * s6i2 -
				  cos(  Psi - 3. * alpha) * s4i2 * c2i2 +
				  cos(  Psi + 3. * alpha) * s2i2 * c4i2 +
				  9. * cos(3. * (Psi + alpha)) * c6i2) +
				 v2 * 4. * si * (1. - 3. * mparams.eta) * 
				 (    - cos(2. * Psi - 3. * alpha) * s4i2 +
				      cos(2. * Psi + 3. * alpha) * c4i2 ) );
    
    h3P3->data[2 * j + 1] = amp33 * ( v * mparams.dm *
				      (-9. * sin(3. * (Psi - alpha)) * s6i2 -
				       sin(  Psi - 3. * alpha) * s4i2 * c2i2 -
				       sin(  Psi + 3. * alpha) * s2i2 * c4i2 -
				       9. * sin(3. * (Psi + alpha)) * c6i2) +
				      v2 * 4. * si * (1. - 3. * mparams.eta) * 
				      (  - sin(2. * Psi - 3. * alpha) * s4i2
					 - sin(2. * Psi + 3. * alpha) * c4i2 ) );
    
    h3M3->data[2 * j + 1] = amp33 * ( v * mparams.dm *
				      (-9. * sin(3. * (Psi - alpha)) * s6i2 -
				       sin(  Psi - 3. * alpha) * s4i2 * c2i2 -
				       sin(  Psi + 3. * alpha) * s2i2 * c4i2 -
				       9. * sin(3. * (Psi + alpha)) * c6i2) -
				      v2 * 4. * si * (1. - 3. * mparams.eta) * 
				      (   - sin(2. * Psi - 3. * alpha) * s5i2 * ci2
					  - sin(2. * Psi + 3. * alpha) * c5i2 * si2 ) );

    h3P2->data[2 * j] = amp32 * ( v * mparams.dm / 3. * 
				  ( 27. * cos(3. * Psi - 2. * alpha) * si*s4i2 + 
				    27. * cos(3. * Psi + 2. * alpha) * si*c4i2 +
				    cos( Psi + 2. * alpha) * c3i2 * (5.*si2-3.*si*ci2-3.*ci*si2) /2. +
				    cos( Psi - 2. * alpha) * s3i2 * (5.*ci2+3.*ci*ci2-3.*si*si2) /2.) +
				  v2*(1./3.-mparams.eta) * 
				  ( - 8.*c4i2*(3.*ci-2.)*cos(2.*(Psi+alpha)) +
				    8.*s4i2*(3.*ci+2.)*cos(2.*(Psi-alpha)) ) );

    h3M2->data[2 * j] = amp32 * ( v * mparams.dm / 3. *
				  ( 27. * cos(3. * Psi - 2. * alpha) * si*s4i2 +
				    27. * cos(3. * Psi + 2. * alpha) * si*c4i2 +
				    cos( Psi + 2. * alpha) * c3i2 * (5.*si2-3.*si*ci2-3.*ci*si2) /2. +
				    cos( Psi - 2. * alpha) * s3i2 * (5.*ci2+3.*ci*ci2-3.*si*si2) /2.) +
				  v2*(1./3.-mparams.eta) * 
				  ( 8.*c4i2*(3.*ci-2.)*cos(2.*(Psi+alpha)) -
				    8.*s4i2*(3.*ci+2.)*cos(2.*(Psi-alpha)) ) );

    h3P2->data[2 * j + 1 ] = amp33 * ( v * mparams.dm / 3. * 
				       ( 27. * sin(3. * Psi - 2. * alpha) * si*s4i2 - 
					 27. * cos(3. * Psi + 2. * alpha) * si*c4i2 -
					 sin( Psi + 2. * alpha) * c3i2 * (5.*si2-3.*si*ci2-3.*ci*si2) /2. +
					 sin( Psi - 2. * alpha) * s3i2 * (5.*ci2+3.*ci*ci2-3.*si*si2) /2.) +
				       v2*(1./3.-mparams.eta) *
				       ( 8.*c4i2*(3.*ci-2.)*sin(2.*(Psi+alpha)) +
					 8.*s4i2*(3.*ci+2.)*sin(2.*(Psi-alpha)) ) );

    
    h3P2->data[2 * j + 1 ] = amp32 * ( -v * mparams.dm / 3. * 
				       ( 27. * sin(3. * Psi - 2. * alpha) * si*s4i2 - 
					 27. * cos(3. * Psi + 2. * alpha) * si*c4i2 -
					 sin( Psi + 2. * alpha) * c3i2 * (5.*si2-3.*si*ci2-3.*ci*si2) /2.+
					 sin( Psi - 2. * alpha) * s3i2 * (5.*ci2+3.*ci*ci2-3.*si*si2) /2.)+
				       v2*(1./3.-mparams.eta) *
				       ( 8.*c4i2*(3.*ci-2.)*sin(2.*(Psi+alpha)) +
					 8.*s4i2*(3.*ci+2.)*sin(2.*(Psi-alpha)) ) );

    h3P1->data[2 * j] = amp31 * ( v * mparams.dm / 6. * 
				  ( -135. * cos(3.*Psi - alpha) * s2i*s2i2 +
				    135. * cos(3.*Psi + alpha)  * s2i*c2i2 +
				    cos(Psi+alpha) * c2i2*(15.*cdi-20.*ci+13.)/2.-
				    cos(Psi-alpha) * s2i2*(15.*cdi+20.*ci+13.)/2. )
				  -v2 * (1./3.-mparams.eta)* 
				  ( 20.*c3i2*cos(2.*Psi+alpha)*(3.*(si2*ci+ci2*si)-5.*si2) +
				    20.*s3i2*cos(2.*Psi-alpha)*(3.*(ci2*ci-si2*si)+5.*ci2) ) );

    h3M1->data[2 * j] = amp31 * (-v * mparams.dm / 6. *
				 ( -135. * cos(3.*Psi - alpha) * s2i*s2i2 +
				   135. * cos(3.*Psi + alpha) * s2i*c2i2 +
				   cos(Psi+alpha) * c2i2*(15.*cdi-20.*ci+13.)/2.-
				   cos(Psi-alpha) * s2i2*(15.*cdi+20.*ci+13.)/2. )
				 -v2 * (1./3.-mparams.eta)*
				 ( 20.*c3i2*cos(2.*Psi+alpha)*(3.*(si2*ci+ci2*si)-5.*si2) +
				   20.*s3i2*cos(2.*Psi-alpha)*(3.*(ci2*ci-si2*si)+5.*ci2) ) );

    h3P1->data[2 * j + 1] = amp31 * ( v * mparams.dm / 6. *
				      ( -135. * sin(3.*Psi - alpha) * s2i*s2i2 -
					135.* sin(3.*Psi + alpha) * s2i*c2i2 -
					sin(Psi+alpha) * c2i2*(15.*cdi-20.*ci+13.)/2.-
					sin(Psi-alpha) * s2i2*(15.*cdi+20.*ci+13.)/2. )
				      +v2 * (1./3.-mparams.eta)*
				      ( 20.*c3i2*sin(2.*Psi+alpha)*(3.*(si2*ci+ci2*si)-5.*si2)
					-20.*s3i2*sin(2.*Psi-alpha)*(3.*(ci2*ci-si2*si)+5.*ci2) ) );

    h3M1->data[2 * j + 1] = amp31 * ( v * mparams.dm / 6. *
				      ( -135. * sin(3.*Psi - alpha) *s2i*s2i2 -
					135. * sin(3.*Psi + alpha) *s2i*c2i2 -
					sin(Psi+alpha) * c2i2*(15.*cdi-20.*ci+13.)/2.-
					sin(Psi-alpha) * s2i2*(15.*cdi+20.*ci+13.)/2. )
				      -v2 * (1./3.-mparams.eta)* 
				      ( 20.*c3i2*sin(2.*Psi+alpha)*(3.*(si2*ci+ci2*si)-5.*si2)
					-20.*s3i2*sin(2.*Psi-alpha)*(3.*(ci2*ci-si2*si)+5.*ci2) ) );

    h30->data[2 * j] = 0.;
    
    h30->data[2 * j + 1] = amp30 * ( v * mparams.dm *
				     ( cos(Psi) * si*(cos(2.*Psi)*(45.*s2i)-(25.*cdi-21.) ) )
				     +v2 * (1.-3.*mparams.eta) *
				     (80. * s2i*c2i*sin(2.*Psi) ) );

    h4P4->data[2 * j] = amp44 * (1. - 3.*mparams.eta) * 
                                  ( 4.* s8i2 * cos(4.*(Psi-alpha)) + cos(2.*Psi-4.*alpha) *s6i2*c2i2 
				    + s2i2*c6i2* cos(2.*Psi+4.*alpha) + 4.*c8i2* cos(4.*(Psi+alpha)) );

    h4M4->data[2 * j] = amp44 * (1. - 3.*mparams.eta) * 
                                  ( 4.* s8i2 * cos(4.*(Psi-alpha)) + cos(2.*Psi-4.*alpha) *s6i2*c2i2 
				    + s2i2*c6i2* cos(2.*Psi+4.*alpha) + 4.*c8i2* cos(4.*(Psi+alpha)) );

    h4P4->data[2 * j + 1] = amp44 * (1. - 3.*mparams.eta) * 
                                  ( 4.* s8i2 * sin(4.*(Psi-alpha)) + sin(2.*Psi-4.*alpha) *s6i2*c2i2 
				    - s2i2*c6i2* sin(2.*Psi+4.*alpha) - 4.*c8i2* sin(4.*(Psi+alpha)) );

    h4M4->data[2 * j + 1] = - amp44 * (1. - 3.*mparams.eta) * 
                                  ( 4.* s8i2 * sin(4.*(Psi-alpha)) + sin(2*Psi-4.*alpha) *s6i2*c2i2 
				    - s2i2*c6i2* sin(2.*Psi+4.*alpha) - 4.*c8i2* sin(4.*(Psi+alpha)) );

    h4P3->data[2 * j] = amp43 * (1. - 3.*mparams.eta) * 
      ( 4.*s6i2* cos(4.*Psi-3.*alpha) - 4.*c6i2* cos(4.*Psi+3.*alpha) - 
	s4i2*(ci+0.5)/2. * cos(2.*Psi-3.*alpha) - c4i2*(ci-0.5) * cos(2.*Psi+3.*alpha) );

    h4M3->data[2 * j] = - amp43 * (1. - 3.*mparams.eta) * 
      ( 4.*s6i2* cos(4.*Psi-3.*alpha) - 4.*c6i2* cos(4.*Psi+3.*alpha) - 
	s4i2*(ci+0.5)/2. * cos(2.*Psi-3.*alpha) - c4i2*(ci-0.5) * cos(2.*Psi+3.*alpha) );

    h4P3->data[2 * j + 1] = amp43 * (1. - 3.*mparams.eta) * 
      ( 4.*s6i2* sin(4.*Psi-3.*alpha) + 4.*c6i2* sin(4.*Psi+3.*alpha) - 
	s4i2*(ci+0.5)/2. * sin(2.*Psi-3.*alpha) + c4i2*(ci-0.5) * sin(2.*Psi+3.*alpha) );

    h4M3->data[2 * j + 1] = amp43 * (1. - 3.*mparams.eta) * 
      ( 4.*s6i2* sin(4.*Psi-3.*alpha) + 4.*c6i2* sin(4.*Psi+3.*alpha) - 
	s4i2*(ci+0.5)/2. * sin(2.*Psi-3.*alpha) + c4i2*(ci-0.5) * sin(2.*Psi+3.*alpha) );

    h4P2->data[2 * j] = amp42 * (1. - 3.*mparams.eta) *
      ( 224.*s6i2*c2i2 * cos(4.*Psi-2.*alpha) - 224.*c6i2*s2i2 * cos(4.*Psi+2.*alpha) 
	- c4i2 * cos(2.*(Psi+alpha))*(7.*cdi-14.*ci+9) - s4i2 * cos(2.*(Psi-alpha))*(7.*cdi+14*ci+9) );

    h4M2->data[2 * j] = amp42 * (1. - 3.*mparams.eta) *
      ( 224.*s6i2*c2i2 * cos(4.*Psi-2.*alpha) - 224.*c6i2*s2i2 * cos(4.*Psi+2.*alpha) 
	- c4i2 * cos(2.*(Psi+alpha))*(7.*cdi-14.*ci+9) - s4i2 * cos(2.*(Psi-alpha))*(7.*cdi+14*ci+9) );

    h4P2->data[2 * j + 1] = amp42 * (1. - 3.*mparams.eta) *
      ( 224.*s6i2*c2i2 * sin(4.*Psi-2.*alpha) + 224.*c6i2*s2i2 * sin(4.*Psi+2.*alpha) 
	+ c4i2 * sin(2.*(Psi+alpha))*(7.*cdi-14.*ci+9) - s4i2 * sin(2.*(Psi-alpha))*(7.*cdi+14*ci+9) );

    h4M2->data[2 * j + 1] = -amp42 * (1. - 3.*mparams.eta) *
      ( 224.*s6i2*c2i2 * sin(4.*Psi-2.*alpha) + 224.*c6i2*s2i2 * sin(4.*Psi+2.*alpha) 
	+ c4i2 * sin(2.*(Psi+alpha))*(7.*cdi-14.*ci+9) - s4i2 * sin(2.*(Psi-alpha))*(7.*cdi+14*ci+9) );
    
    h4P1->data[2 * j] = amp41 * (1. - 3.*mparams.eta) * 
      ( 448.*s5i2*c3i2 * cos(4.*Psi-alpha) - 448.*s3i2*c5i2 * cos(4.*Psi+alpha) +
	s3i2*cos(2.*Psi-alpha)*(7.*(cdi*ci2-sdi*si2)+14.*(ci2*ci-si2*si)+19.*ci2) -
	c3i2*cos(2.*Psi+alpha)*(7.*(cdi*si2+sdi*ci2)-14.*(si*ci2+ci*si2)+19.*ci2) );

    h4M1->data[2 * j] = -amp41 * (1. - 3.*mparams.eta) * 
      ( 448.*s5i2*c3i2 * cos(4.*Psi-alpha) - 448.*s3i2*c5i2 * cos(4.*Psi+alpha) +
	s3i2*cos(2.*Psi-alpha)*(7.*(cdi*ci2-sdi*si2)+14.*(ci2*ci-si2*si)+19.*ci2) -
	c3i2*cos(2.*Psi+alpha)*(7.*(cdi*si2+sdi*ci2)-14.*(si*ci2+ci*si2)+19.*ci2) );

    h4P1->data[2 * j + 1] = amp41 * (1. - 3.*mparams.eta) * 
      ( 448.*s5i2*c3i2 * sin(4.*Psi-alpha) + 448.*s3i2*c5i2 * sin(4.*Psi+alpha) +
	s3i2*sin(2.*Psi-alpha)*(7.*(cdi*ci2-sdi*si2)+14.*(ci2*ci-si2*si)+19.*ci2) +
	c3i2*sin(2.*Psi+alpha)*(7.*(cdi*si2+sdi*ci2)-14.*(si*ci2+ci*si2)+19.*ci2) );

    h4M1->data[2 * j + 1] = amp41 * (1. - 3.*mparams.eta) * 
      ( 448.*s5i2*c3i2 * sin(4.*Psi-alpha) + 448.*s3i2*c5i2 * sin(4.*Psi+alpha) +
	s3i2*sin(2.*Psi-alpha)*(7.*(cdi*ci2-sdi*si2)+14.*(ci2*ci-si2*si)+19.*ci2) +
	c3i2*sin(2.*Psi+alpha)*(7.*(cdi*si2+sdi*ci2)-14.*(si*ci2+ci*si2)+19.*ci2) );

    h40->data[2 * j] = amp40 * (1.-3.*mparams.eta) * s2i * (-56.*cos(4.*Psi) - 
							       cos(2.*Psi)*(7.*cdi+5.) );
    h40->data[2 * j +1] = 0.;    
    
    fap->data[j] = om;
    phap->data[j] = Psi;
    
    j++;
    
  } while (om < (fracRD * omegaRD));

  XLALDestroyREAL8Vector(omega_s);
  XLALDestroyREAL8Vector(Lhx_s);
  XLALDestroyREAL8Vector(Lhy_s);
  XLALDestroyREAL8Vector(Lhz_s);

  XLALDestroyREAL8Vector(domega_s);
  XLALDestroyREAL8Vector(dLhx_s);
  XLALDestroyREAL8Vector(dLhy_s);
  XLALDestroyREAL8Vector(dLhz_s);
  XLALDestroyREAL8Vector(diota);
  XLALDestroyREAL8Vector(dalpha);

  XLALDestroyREAL8Vector(ddomega_s);
  XLALDestroyREAL8Vector(ddiota);
  XLALDestroyREAL8Vector(ddalpha);
  
  *count=j;
  
  /*Now ringdown is inserted */
  
  /*--------------------------------------------------------------
   * Attach the ringdown waveform to the end of inspiral
   -------------------------------------------------------------*/
  
  apcount  = j;
  errcode  = XLALPSpinInspiralAttachRingdownWave(h2P2, params, &apcount, nmodes, 2, 2, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h2P2->data[i] = 0.;
  if (apcount > *count) *count = apcount;

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h2M2, params, &apcount, nmodes, 2, -2, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h2M2->data[i] = 0.;
  if (apcount > *count) *count = apcount;
  
  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h2P1, params, &apcount, nmodes, 2, 1, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h2P1->data[i] = 0.;  

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h2M1, params, &apcount, nmodes, 2, -1, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h2M1->data[i] = 0.;  

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h20, params, &apcount, nmodes, 2, 0, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h20->data[i] = 0.;  

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h3P3, params, &apcount, nmodes, 3, 3, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h3P3->data[i] = 0.;  

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h3M3, params, &apcount, nmodes, 3, -3, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h3M3->data[i] = 0.;  

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h3P2, params, &apcount, nmodes, 3, 2, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h3P2->data[i] = 0.;  

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h3M2, params, &apcount, nmodes, 3, -2, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h3P2->data[i] = 0.;  

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h3P1, params, &apcount, nmodes, 3, 1, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h3P1->data[i] = 0.;  

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h3M1, params, &apcount, nmodes, 3, -1, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h3M1->data[i] = 0.;  

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h30, params, &apcount, nmodes, 3, 0, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h30->data[i] = 0.;  

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h4P4, params, &apcount, nmodes, 4, 4, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h4P4->data[i] = 0.;  

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h4M4, params, &apcount, nmodes, 4, -4, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h4M4->data[i] = 0.;  

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h4P3, params, &apcount, nmodes, 4, 3, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h4P3->data[i] = 0.;  

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h4M3, params, &apcount, nmodes, 4, -3, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h4M3->data[i] = 0.;  

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h4P2, params, &apcount, nmodes, 4, 2, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h4P4->data[i] = 0.;  

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h4M2, params, &apcount, nmodes, 4, -2, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h4M4->data[i] = 0.;  

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h4P1, params, &apcount, nmodes, 4, 1, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h4P3->data[i] = 0.;  

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h4M1, params, &apcount, nmodes, 4, -1, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h4M3->data[i] = 0.;  

  apcount  = j;
  errcode += XLALPSpinInspiralAttachRingdownWave(h40, params, &apcount, nmodes, 4, 0, finalMass, finalSpin);
  for (i = 2 * apcount; i < 2 * length; i++) h40->data[i] = 0.;  

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
    fprintf(stderr,"** LALPSpinInspiralRD ERROR **: impossible to create RingDownWave\n");
    XLAL_ERROR(func,XLAL_EFAILED);
  }

  /*-------------------------------------------------------------------
   * Compute the spherical harmonics required for constructing (h+,hx).
   -------------------------------------------------------------------*/
  
  /* The angles theta for the spherical harmonics has been set according to 
     the input inclination parameter and the axisChoice */
  
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
    fprintf(stderr,"** LALPSpinInspiralRD ERROR **: impossible to create Y22 or Y2-2\n");
    XLAL_ERROR(func,XLAL_EFAILED);
  }
  for (i = 0; i < length; i++) {
    fap->data[i] /= unitHz;
    x0 = h2P2->data[2 * i];
    x1 = h2P2->data[2 * i + 1];
    x2 = h2M2->data[2 * i];
    x3 = h2M2->data[2 * i + 1];
    sigp->data[i] =   x0 * MultSphHarmP.re - x1 * MultSphHarmP.im + x2 * MultSphHarmM.re - x3 * MultSphHarmM.im;
    sigc->data[i] = - x0 * MultSphHarmP.im - x1 * MultSphHarmP.re - x2 * MultSphHarmM.im - x3 * MultSphHarmM.re;
  }
  
  errcode  = XLALSphHarm(&MultSphHarmP, 2, 1, inc, 0.);
  errcode += XLALSphHarm(&MultSphHarmM, 2, -1, inc, 0.);
  if (errcode != XLAL_SUCCESS){
    XLALDestroyREAL8Vector(h2P1);
    XLALDestroyREAL8Vector(h2M1);
    fprintf(stdout, "** LALPSpinInspiralRD WARNING **: impossible to create Y21\n");
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
    fprintf(stdout, "** LALPSpinInspiralRD WARNING **: impossible to create Y20\n");
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
    fprintf(stdout, "** LALPSpinInspiralRD WARNING **: impossible to create Y33,Y3-3\n");
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
    fprintf(stdout, "** LALPSpinInspiralRD WARNING **: impossible to create Y32,Y3-2\n");
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
    fprintf(stdout, "** LALPSpinInspiralRD WARNING **: impossible to create Y31,Y3-1\n");
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
    fprintf(stdout, "** LALPSpinInspiralRD WARNING **: impossible to create Y30\n");
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
    fprintf(stdout, "** LALPSpinInspiralRD WARNING **: impossible to create Y44,Y4-3\n");
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
    fprintf(stdout, "** LALPSpinInspiralRD WARNING **: impossible to create Y43,Y4-3\n");
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
    fprintf(stdout, "** LALPSpinInspiralRD WARNING **: impossible to create Y32,Y3-2\n");
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
    fprintf(stdout, "** LALPSpinInspiralRD WARNING **: impossible to create Y31,Y3-1\n");
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
    fprintf(stdout, "** LALPSpinInspiralRD WARNING **: impossible to create Y30\n");
  } else {
    for (i = 0; i < length; i++) {
      x0 = h40->data[2 * i];
      x1 = h40->data[2 * i + 1];    
      sigp->data[i] += x0 * MultSphHarmP.re - x1 * MultSphHarmM.im;
      sigc->data[i] -= x0 * MultSphHarmP.im + x1 * MultSphHarmP.re;
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

  return errcode;

  /*End */

}
