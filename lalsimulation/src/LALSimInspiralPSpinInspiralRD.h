/*
 *  LALSimInspiralPSpinInspiralRD.h
 *  
 *
 *  Created by John Veitch on 23/08/2011.
 *  (C) John Veitch 2011, all rights reserved.
 *
 */

/**
 * Convenience function to set up LALPSpinInspiralRDparams struct
 */

#include "LALSimInspiral.h"

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

typedef enum {
  TotalJ,
  View,
  OrbitalL,
} InputAxis;


int XLALPSpinInspiralRDparamsSetup(LALPSpinInspiralRDparams *mparams, /** Output: RDparams structure */
																	 UINT4 inspiralOnly, 	/** Only generate inspiral */
																	 REAL8 deltaT, 		/** sampling interval */
																	 REAL8 fLow,		/** Starting frequency */
																	 REAL8 fCutoff,		/** CHECKME: Cutoff frequency? */
																	 REAL8 m1,			/** Mass 1 */
																	 REAL8 m2,			/** Mass 2 */
																	 LALSimSpinInteraction spinInteraction,	/** Spin interaction */
																	 UINT4 order		/** twice PN Order in Phase */ 
																	 );

int XLALSimInspiralPSpinInspiralRDGenerator(
																						REAL8TimeSeries *hplus,  /**< +-polarization waveform */
																						REAL8TimeSeries *hcross, /**< x-polarization waveform */
																						REAL8 phi0,               /**< start phase */
																						REAL8 deltaT,             /**< sampling interval */
																						REAL8 m1,                 /**< mass of companion 1 */
																						REAL8 m2,                 /**< mass of companion 2 */
																						REAL8 f_min,              /**< start frequency */
																						REAL8 r,                  /**< distance of source */
																						REAL8 iota,               /**< inclination of source (rad) */
																						REAL8 spin1[3],						/**< Spin vector on mass1 */
																						REAL8 spin2[3],						/**< Spin vector on mass2 */
																						int phaseO,               /**< twice post-Newtonian phase order */
																						InputAxis axisChoice			/**< Choice of axis for input spin params */
																						);
