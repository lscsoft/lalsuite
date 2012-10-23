/*
 * Copyright (C) 2011 Drew Keppel, J. Creighton, S. Fairhurst, B. Krishnan, L. Santamaria, Stas Babak, David Churches, B.S. Sathyaprakash, Craig Robinson , Thomas Cokelaer
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

#include <math.h>

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALSimInspiral.h>
#define LAL_USE_COMPLEX_SHORT_MACROS
#include <lal/FindRoot.h>
#include <lal/LALComplex.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#include "LALSimInspiraldEnergyFlux.c"
#include "LALSimInspiralPNCoefficients.c"

#include "check_series_macros.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

typedef struct
tagexpnCoeffsTaylorT3 {
   int ieta;
   /* Taylor expansion coefficents in phi(t)*/
   REAL8 ptaN, pta2, pta3, pta4, pta5, pta6, pta7, ptl6, pta10, pta12;
   /* Taylor expansion coefficents in f(t)*/
   REAL8 ftaN, fta2, fta3, fta4, fta5, fta6, fta7, ftl6, fta10, fta12;

   /* sampling interval*/
   REAL8 samplinginterval;
   /* symmetric mass ratio, total mass, component masses*/
   REAL8 eta, totalmass, m1, m2, chi1, chi2;

   /* initial and final values of frequency, time, velocity; lso
    values of velocity and frequency; final phase.*/
   REAL8 f0, fn, t0, theta_lso, v0, vn, vf, vlso, flso, phiC;

   /* last stable orbit and pole defined by various Taylor and P-approximants*/
   REAL8 vlsoT0, vlsoT2, vlsoT4, vlsoT6;
}  expnCoeffsTaylorT3;

typedef REAL8 (SimInspiralPhasing3)(
   REAL8 td,
   expnCoeffsTaylorT3 *ak);

typedef REAL8 (SimInspiralFrequency3)(
   REAL8 td,
   expnCoeffsTaylorT3 *ak);

typedef struct
tagexpnFuncTaylorT3
{
   SimInspiralPhasing3 *phasing3;
   SimInspiralFrequency3 *frequency3;
} expnFuncTaylorT3;

typedef struct
{
	REAL8 (*func)(REAL8 tC, expnCoeffsTaylorT3 *ak);
	expnCoeffsTaylorT3 ak;
}
FreqInFromChirptime;

static REAL8 XLALInspiralFrequency3Wrapper(REAL8 tC, void *pars)
{
  FreqInFromChirptime *in;
  REAL8 freq, f;

  in = (FreqInFromChirptime *) pars;
  freq = in->func(tC, &(in->ak));
  if (XLAL_IS_REAL8_FAIL_NAN(freq))
    XLAL_ERROR_REAL8(XLAL_EFUNC);
  f = freq - in->ak.f0;

  /*
  fprintf(stderr, "Here freq=%e f=%e tc=%e f0=%e\n", freq, *f, tC, in->ak.f0);
   */

  return f;
}

static REAL8
XLALSimInspiralFrequency3_0PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta3;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta3 = theta*theta*theta;

  frequency = theta3*ak->ftaN;

  return frequency;
}


static REAL8
XLALSimInspiralFrequency3_2PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3, theta10, theta12;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta10 = theta2*theta2*theta2*theta2*theta2;
  theta12 = theta10*theta2;

  frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta10*theta10
             + ak->fta12*theta12);

  return frequency;
}


static REAL8
XLALSimInspiralFrequency3_3PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3, theta10, theta12;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta10 = theta3*theta3*theta3*theta;
  theta12 = theta10*theta2;

  frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3
             + ak->fta10*theta10
             + ak->fta12*theta12);

  return frequency;
}


static REAL8
XLALSimInspiralFrequency3_4PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3,theta4, theta10, theta12;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta10 = theta4*theta4*theta2;
  theta12 = theta10*theta2;

  frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3
             + ak->fta4*theta4
             + ak->fta10*theta10
             + ak->fta12*theta12);

  return frequency;
}


static REAL8
XLALSimInspiralFrequency3_5PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5, theta10, theta12;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta10 = theta5*theta5;
  theta12 = theta10*theta2;

  frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3
             + ak->fta4*theta4
             + ak->fta5*theta5
             + ak->fta10*theta10
             + ak->fta12*theta12);

  return frequency;
}


static REAL8
XLALSimInspiralFrequency3_6PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5,theta6, theta10, theta12;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta6 = theta5*theta;
  theta10 = theta6*theta4;
  theta12 = theta10*theta2;

  frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3
             + ak->fta4*theta4
             + ak->fta5*theta5
             + (ak->fta6 + ak->ftl6*log(2.*theta))*theta6
             + ak->fta10*theta10
             + ak->fta12*theta12);

  return frequency;
}


static REAL8
XLALSimInspiralFrequency3_7PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5,theta6,theta7, theta10, theta12;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta6 = theta5*theta;
  theta7 = theta6*theta;
  theta10 = theta7*theta3;
  theta12 = theta10*theta2;

  frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3
             + ak->fta4*theta4
             + ak->fta5*theta5
             + (ak->fta6 + ak->ftl6*log(2.*theta))*theta6
             + ak->fta7*theta7
             + ak->fta10*theta10
             + ak->fta12*theta12);

  return frequency;
}


static REAL8
XLALSimInspiralPhasing3_0PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta5;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta5 = pow(td,-0.625);
  phase = (ak->ptaN/theta5);

  return phase;
}


static REAL8
XLALSimInspiralPhasing3_2PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta5, theta10, theta12;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta5 = theta2*theta2*theta;
  theta10 = theta5*theta5;
  theta12 = theta10*theta2;

  phase = (ak->ptaN/theta5) * (1.
         + ak->pta2*theta2
         + ak->pta10*theta10
         + ak->pta12*theta12);

  return phase;
}


static REAL8
XLALSimInspiralPhasing3_3PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3,theta5, theta10, theta12;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta5 = theta2*theta3;
  theta10 = theta5*theta5;
  theta12 = theta10*theta2;

  phase = (ak->ptaN/theta5) * (1.
         + ak->pta2*theta2
         + ak->pta3*theta3
         + ak->pta10*theta10
         + ak->pta12*theta12);

  return phase;
}


static REAL8
XLALSimInspiralPhasing3_4PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5, theta10, theta12;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta10 = theta5*theta5;
  theta12 = theta10*theta2;

  phase = (ak->ptaN/theta5) * (1.
         + ak->pta2*theta2
         + ak->pta3*theta3
         + ak->pta4*theta4
         + ak->pta10*theta10
         + ak->pta12*theta12);

  return phase;
}


static REAL8
XLALSimInspiralPhasing3_5PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5, theta10, theta12;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta10 = theta5*theta5;
  theta12 = theta10*theta2;

  phase = (ak->ptaN/theta5) * (1.
         + ak->pta2*theta2
         + ak->pta3*theta3
         + ak->pta4*theta4
         + ak->pta5 * log(theta/ak->theta_lso) * theta5
         + ak->pta10*theta10
         + ak->pta12*theta12);

  return phase;
}


static REAL8
XLALSimInspiralPhasing3_6PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5,theta6, theta10, theta12;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta6 = theta5*theta;
  theta10 = theta6*theta4;
  theta12 = theta10*theta2;

  phase = (ak->ptaN/theta5) * (1.
         + ak->pta2*theta2
         + ak->pta3*theta3
         + ak->pta4*theta4
         + ak->pta5*log(theta/ak->theta_lso)*theta5
         +(ak->ptl6*log(2.*theta) + ak->pta6)*theta6
         + ak->pta10*theta10
         + ak->pta12*theta12);

  return phase;
}


static REAL8
XLALSimInspiralPhasing3_7PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5,theta6,theta7, theta10, theta12;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta6 = theta5*theta;
  theta7 = theta6*theta;
  theta10 = theta7*theta3;
  theta12 = theta10*theta2;

  phase = (ak->ptaN/theta5) * (1.
         + ak->pta2*theta2
         + ak->pta3*theta3
         + ak->pta4*theta4
         + ak->pta5*log(theta/ak->theta_lso)*theta5
         +(ak->ptl6*log(2.*theta) + ak->pta6)*theta6
         + ak->pta7*theta7
         + ak->pta10*theta10
         + ak->pta12*theta12);

  return phase;
}


/**
 * Returns the sum of chirp times to a specified order.
 *
 * Computes the sum of the chirp times to a specified order. Inputs given in SI
 * units.
 */
static REAL8 XLALSimInspiralChirpLength(
		REAL8 m1,					/**< mass of companion 1 */
		REAL8 m2,					/**< mass of companion 2 */
		REAL8 lambda1,					/**< (tidal deformability of body 1)/(mass of body 1)^5 */
		REAL8 lambda2,					/**< (tidal deformability of body 2)/(mass of body 2)^5 */
		LALSimInspiralInteraction interactionFlags, 	/**< flag to control spin and tidal effects */
		REAL8 f_min,					/**< start frequency */
		int O						/**< twice post-Newtonian order */
		)
{
	REAL8 tN, t2, t3 ,t4, t5, t6, t6l, t7, tC, t10, t12;
	REAL8 v, v2, v3, v4, v5, v6, v7, v8, v10, v12;
	REAL8 piFl = LAL_PI * f_min;
	REAL8 m = m1 + m2;
	REAL8 chi1 = m1/(m1+m2), chi2 = m2/(m1+m2);
	REAL8 eta = m1 * m2 / m / m;
	m *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m from kilograms to seconds */

	/* Should use the coefficients from LALInspiraSetup.c to avoid errors.
	 */
	v = cbrt(piFl * m);
	v2 = v*v;
	v3 = v*v2;
	v4 = v*v3;
	v5 = v*v4;
	v6 = v*v5;
	v7 = v*v6;
	v8 = v*v7;
	v10 = v2*v8;
	v12 = v2*v10;

	tN = XLALSimInspiralTaylorT2Timing_0PNCoeff(m1+m2, eta);
	t2 = XLALSimInspiralTaylorT2Timing_2PNCoeff(eta);
	t3 = XLALSimInspiralTaylorT2Timing_3PNCoeff(eta);
	t4 = XLALSimInspiralTaylorT2Timing_4PNCoeff(eta);
	t5 = XLALSimInspiralTaylorT2Timing_5PNCoeff(eta);
	t6 = XLALSimInspiralTaylorT2Timing_6PNCoeff(eta);
	t7 = XLALSimInspiralTaylorT2Timing_7PNCoeff(eta);
	t6l = XLALSimInspiralTaylorT2Timing_6PNLogCoeff(eta);

	/* Tidal co-efficients for t(v). */
	t10 = 0.;
	t12 = 0.;
	if( interactionFlags >= LAL_SIM_INSPIRAL_INTERACTION_TIDAL_5PN)	
	{
		t10 = XLALSimInspiralTaylorT2Timing_10PNTidalCoeff(chi1,lambda1)
		    + XLALSimInspiralTaylorT2Timing_10PNTidalCoeff(chi2,lambda2);
	}
	if( interactionFlags >= LAL_SIM_INSPIRAL_INTERACTION_TIDAL_6PN )
	{
		t12 = XLALSimInspiralTaylorT2Timing_12PNTidalCoeff(eta,chi1,lambda1)
		    + XLALSimInspiralTaylorT2Timing_12PNTidalCoeff(eta,chi2,lambda2);
	}

	switch (O) {
		case 0:
		case 1:
			t2 = 0.;
		case 2:
			t3 = 0.;
		case 3:
			t4 = 0.;
		case 4:
			t5 = 0.;
		case 5:
			t6 = 0.;
			t6l = 0.;
		case 6:
			t7 = 0.;
		case 7:
        case -1: // Use the max PN order, move if higher orders implemented
			break;
		case 8:
			XLALPrintError("XLAL Error - %s: Not supported for requested PN order\n", __func__);
			XLAL_ERROR_REAL8(XLAL_EINVAL);
			break;
		default:
			XLALPrintError("XLAL Error - %s: Unknown PN order in switch\n", __func__);
			XLAL_ERROR_REAL8(XLAL_EINVAL);
	}

	tC = -tN / v8 * (1.
		+ t2 * v2
		+ t3 * v3
		+ t4 * v4
		+ t5 * v5
		+ (t6 + t6l*log(16*v2)) * v6
		+ t7 * v7
		+ t10 * v10
		+ t12 * v12);

	return tC;
}


/**
 * Set up the expnCoeffsTaylorT3 and expnFuncTaylorT3 structures for
 * generating a TaylorT3 waveform.
 *
 * Inputs given in SI units.
 */
static int XLALSimInspiralTaylorT3Setup(
		expnCoeffsTaylorT3 *ak,				/**< coefficients for TaylorT3 evolution [modified] */
		expnFuncTaylorT3 *f,				/**< functions for TaylorT3 evolution [modified] */
	       	REAL8 deltaT,					/**< sampling interval */
		REAL8 m1,					/**< mass of companion 1 */
		REAL8 m2,					/**< mass of companion 2 */
		REAL8 lambda1,					/**< (tidal deformability of body 1)/(mass of body 1)^5 */
		REAL8 lambda2,					/**< (tidal deformability of body 2)/(mass of body 2)^5 */
		LALSimInspiralInteraction interactionFlags, 	/**< flag to control spin and tidal effects */
		REAL8 f_min,					/**< start frequency */
		int O						/**< twice post-Newtonian order */
		)
{
  REAL8 eta, tn, chi1, chi2;

  ak->t0 = 0;
  ak->m1 = m1;
  ak->m2 = m2;
  ak->totalmass = ak->m1 + ak->m2;
  eta = ak->eta = m1 * m2 / (ak->totalmass * ak->totalmass);
  chi1 = ak->chi1 = m1/ak->totalmass;
  chi2 = ak->chi2 = m2/ak->totalmass;
  ak->totalmass *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert totalmass from kilograms to seconds */

  ak->f0 = f_min;
  ak->samplinginterval = deltaT;
  ak->fn = 1. / (2. * ak->samplinginterval);
  ak->v0 = cbrt(LAL_PI * ak->totalmass * f_min);

  ak->ptaN = XLALSimInspiralTaylorT3Phasing_0PNCoeff(eta);
  ak->pta2 = XLALSimInspiralTaylorT3Phasing_2PNCoeff(eta);
  ak->pta3 = XLALSimInspiralTaylorT3Phasing_3PNCoeff(eta);
  ak->pta4 = XLALSimInspiralTaylorT3Phasing_4PNCoeff(eta);
  ak->pta5 = XLALSimInspiralTaylorT3Phasing_5PNCoeff(eta);
  ak->pta6 = XLALSimInspiralTaylorT3Phasing_6PNCoeff(eta);
  ak->pta7 = XLALSimInspiralTaylorT3Phasing_7PNCoeff(eta);
  ak->ptl6 = XLALSimInspiralTaylorT3Phasing_6PNLogCoeff(eta);

  ak->ftaN = XLALSimInspiralTaylorT3Frequency_0PNCoeff(m1+m2);
  ak->fta2 = XLALSimInspiralTaylorT3Frequency_2PNCoeff(eta);
  ak->fta3 = XLALSimInspiralTaylorT3Frequency_3PNCoeff(eta);
  ak->fta4 = XLALSimInspiralTaylorT3Frequency_4PNCoeff(eta);
  ak->fta5 = XLALSimInspiralTaylorT3Frequency_5PNCoeff(eta);
  ak->fta6 = XLALSimInspiralTaylorT3Frequency_6PNCoeff(eta);
  ak->fta7 = XLALSimInspiralTaylorT3Frequency_7PNCoeff(eta);
  ak->ftl6 = XLALSimInspiralTaylorT3Frequency_6PNLogCoeff(eta);

  /* Tidal co-efficients for phasing and frequency */
  ak->pta10 = 0.;
  ak->pta12 = 0.;
  ak->fta10 = 0.;
  ak->fta12 = 0.;
  if( interactionFlags >= LAL_SIM_INSPIRAL_INTERACTION_TIDAL_5PN)
  {
    ak->pta10 = XLALSimInspiralTaylorT3Phasing_10PNTidalCoeff(chi1,lambda1)
              + XLALSimInspiralTaylorT3Phasing_10PNTidalCoeff(chi2,lambda2);
    ak->fta10 = XLALSimInspiralTaylorT3Frequency_10PNTidalCoeff(chi1,lambda1)
              + XLALSimInspiralTaylorT3Frequency_10PNTidalCoeff(chi2,lambda2);
  }
  if( interactionFlags >= LAL_SIM_INSPIRAL_INTERACTION_TIDAL_6PN )
  {
    ak->pta12 = XLALSimInspiralTaylorT3Phasing_12PNTidalCoeff(eta,chi1,lambda1)
              + XLALSimInspiralTaylorT3Phasing_12PNTidalCoeff(eta,chi2,lambda2);
    ak->fta12 = XLALSimInspiralTaylorT3Frequency_12PNTidalCoeff(eta,chi1,lambda1)
              + XLALSimInspiralTaylorT3Frequency_12PNTidalCoeff(eta,chi2,lambda2);
  }

  switch (O)
  {
     case 0:
           f->phasing3 = &XLALSimInspiralPhasing3_0PN;
           f->frequency3 = &XLALSimInspiralFrequency3_0PN;
           break;
     case 1:
           XLALPrintError("XLAL Error - %s: PN approximant not supported for requested PN order\n", __func__);
           XLAL_ERROR(XLAL_EINVAL);
           break;
     case 2:
           f->phasing3 = &XLALSimInspiralPhasing3_2PN;
           f->frequency3 = &XLALSimInspiralFrequency3_2PN;
           break;
     case 3:
           f->phasing3 = &XLALSimInspiralPhasing3_3PN;
           f->frequency3 = &XLALSimInspiralFrequency3_3PN;
           break;
     case 4:
           f->phasing3 = &XLALSimInspiralPhasing3_4PN;
           f->frequency3 = &XLALSimInspiralFrequency3_4PN;
           break;
     case 5:
           f->phasing3 = &XLALSimInspiralPhasing3_5PN;
           f->frequency3 = &XLALSimInspiralFrequency3_5PN;
           break;
     case 6:
           f->phasing3 = &XLALSimInspiralPhasing3_6PN;
           f->frequency3 = &XLALSimInspiralFrequency3_6PN;
           break;
     case 7:
     case -1: // Use the highest PN order available, move if higher terms added
           f->phasing3 = &XLALSimInspiralPhasing3_7PN;
           f->frequency3 = &XLALSimInspiralFrequency3_7PN;
           break;
     case 8:
           XLALPrintError("XLAL Error - %s: PN approximant not supported for requested PN order\n", __func__);
           XLAL_ERROR(XLAL_EINVAL);
           break;
     default:
        XLALPrintError("XLAL Error - %s: Unknown PN order in switch\n", __func__);
        XLAL_ERROR(XLAL_EINVAL);
  }

  tn = XLALSimInspiralTaylorLength(deltaT, m1, m2, f_min, O);
  ak->theta_lso = pow(tn, -0.125);

  return XLAL_SUCCESS;
}


/**
 * Computes a post-Newtonian orbit using the Taylor T3 method.
 */
int XLALSimInspiralTaylorT3PNEvolveOrbit(
		REAL8TimeSeries **V,                        /**< post-Newtonian parameter [returned] */
		REAL8TimeSeries **phi,                      /**< orbital phase [returned] */
		REAL8 phiRef,                               /**< reference orbital phase (rad) */
		REAL8 deltaT,                               /**< sampling interval (s) */
		REAL8 m1,                                   /**< mass of companion 1 (kg) */
		REAL8 m2,                                   /**< mass of companion 2 (kg) */
		REAL8 f_min,                                /**< starting GW frequency (Hz) */
		REAL8 fRef,                                 /**< reference GW frequency (Hz) */
		REAL8 lambda1,                              /**< (tidal deformability of body 1)/(mass of body 1)^5 */
		REAL8 lambda2,                              /**< (tidal deformability of body 2)/(mass of body 2)^5 */
		LALSimInspiralInteraction interactionFlags, /**< flag to control spin and tidal effects */
		int O                                       /**< twice post-Newtonian order */
		)
{

	const UINT4 blocklen = 1024;
	const REAL8 visco = sqrt(1.0/6.0);
	REAL8 m = m1 + m2, VRef = 0.;
	REAL8 nu = m1 * m2 / m / m;
	m *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m from kilograms to seconds */
	REAL8 tmptC, tC, c1, xmin, xmax, xacc, v, phase, fOld, t, td, temp, tempMin = 0, tempMax = 0;
	REAL8 (*freqfunction)(REAL8, void *);
	UINT4 j, len, idxRef = 0;
	LIGOTimeGPS t_end = LIGOTIMEGPSZERO;
	REAL8 f;
	void *pars;

	expnFuncTaylorT3 expnfunc;
	expnCoeffsTaylorT3 ak;
	FreqInFromChirptime timeIn;

	/* allocate memory */

	*V = XLALCreateREAL8TimeSeries("ORBITAL_FREQUENCY_PARAMETER", &t_end, 0.0, deltaT, &lalDimensionlessUnit,
		blocklen);
	*phi = XLALCreateREAL8TimeSeries("ORBITAL_PHASE", &t_end, 0.0, deltaT, &lalDimensionlessUnit, blocklen);
	if (!V || !phi)
		XLAL_ERROR(XLAL_EFUNC);


	/* initialize expnCoeffsTaylorT3 and expnFuncTaylorT3 structures */
	if (XLALSimInspiralTaylorT3Setup(&ak, &expnfunc, deltaT, m1, m2, lambda1,
		lambda2, interactionFlags, f_min, O))
		XLAL_ERROR(XLAL_EFUNC);

	tC = XLALSimInspiralChirpLength(m1, m2, lambda1, lambda2, interactionFlags, f_min, O);
	c1 = nu/(5.*m);

	/*
	 * In Jan 2003 we realized that the tC determined as a sum of chirp
	 * times is not quite the tC that should enter the definition of Theta
	 * in the expression for the frequency as a function of time (see DIS3,
	 * 2000). This is because chirp times are obtained by inverting t(f).
	 * Rather tC should be obtained by solving the equation f0 - f(tC) = 0.
	 * This is what is implemented below.
	 */

	timeIn.func = expnfunc.frequency3;
	timeIn.ak = ak;
	freqfunction = &XLALInspiralFrequency3Wrapper;
	xmin = c1*tC/2.;
	xmax = c1*tC*2.;
	xacc = 1.e-6;
	pars = (void*) &timeIn;
	/* tc is the instant of coalescence */

	/* we add 5 so that if tC is small then xmax
	 * is always greater than a given value (here 5)*/
	xmax = c1*tC*3 + 5.;

	/* for x in [xmin, xmax], we search the value which gives the max
	 * frequency.  and keep the corresponding rootIn.xmin. */

	for (tmptC = c1*tC/1000.; tmptC < xmax; tmptC+=c1*tC/1000.){
		temp = XLALInspiralFrequency3Wrapper(tmptC , pars);
		if (XLAL_IS_REAL8_FAIL_NAN(temp))
			XLAL_ERROR(XLAL_EFUNC);
		if (temp > tempMax) {
			xmin = tmptC;
			tempMax = temp;
		}
		if (temp < tempMin) {
			tempMin = temp;
		}
	}

	/* if we have found a value positive then everything should be fine in
	 * the BissectionFindRoot function */
	if (tempMax > 0  &&  tempMin < 0){
		tC = XLALDBisectionFindRoot (freqfunction, xmin, xmax, xacc, pars);
		if (XLAL_IS_REAL8_FAIL_NAN(tC))
			XLAL_ERROR(XLAL_EFUNC);
	}
	else{
		XLALPrintError("Can't find good bracket for BisectionFindRoot");
		XLAL_ERROR(XLAL_EMAXITER);
	}

	tC /= c1;

	/* start waveform generation */

	t = 0.;
	td = c1 * (tC - t);
	phase = expnfunc.phasing3(td, &ak);
	if (XLAL_IS_REAL8_FAIL_NAN(phase))
		XLAL_ERROR(XLAL_EFUNC);
	f = expnfunc.frequency3(td, &ak);
	if (XLAL_IS_REAL8_FAIL_NAN(f))
		XLAL_ERROR(XLAL_EFUNC);

	v = cbrt(f * LAL_PI * m);
	(*V)->data->data[0] = v;
	(*phi)->data->data[0] = phase;

	j = 0;
	while (1) {

		/* make one step */

		j++;
		fOld = f;
		t = j * deltaT;
		td = c1 * (tC - t);
		phase = expnfunc.phasing3(td, &ak);
		if (XLAL_IS_REAL8_FAIL_NAN(phase))
			XLAL_ERROR(XLAL_EFUNC);
		f = expnfunc.frequency3(td, &ak);
		if (XLAL_IS_REAL8_FAIL_NAN(f))
			XLAL_ERROR(XLAL_EFUNC);
		v = cbrt(f * LAL_PI * m);

		/* check termination conditions */

		if (v >= visco) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated at ISCO\n", __func__);
			break;
		}
		if (f <= fOld) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated when frequency stalled\n", __func__);
			break;
		}
	
		/* save current values in vectors but first make sure we don't write past end of vectors */

		if ( j >= (*V)->data->length ) {
			if ( ! XLALResizeREAL8TimeSeries(*V, 0, (*V)->data->length + blocklen) )
				XLAL_ERROR(XLAL_EFUNC);
			if ( ! XLALResizeREAL8TimeSeries(*phi, 0, (*phi)->data->length + blocklen) )
				XLAL_ERROR(XLAL_EFUNC);
		}
		(*V)->data->data[j] = v;
		(*phi)->data->data[j] = phase;
	}

	/* make the correct length */

	if ( ! XLALResizeREAL8TimeSeries(*V, 0, j) )
		XLAL_ERROR(XLAL_EFUNC);
	if ( ! XLALResizeREAL8TimeSeries(*phi, 0, j) )
		XLAL_ERROR(XLAL_EFUNC);

	/* adjust to correct time */

	XLALGPSAdd(&(*phi)->epoch, -1.0*j*deltaT);
	XLALGPSAdd(&(*V)->epoch, -1.0*j*deltaT);

	/* Do a constant phase shift to get desired value of phiRef */
	len = (*phi)->data->length;
	/* For fRef==0, phiRef is phase of last sample */
	if( fRef == 0. )
		phiRef -= (*phi)->data->data[len-1];
	/* For fRef==fmin, phiRef is phase of first sample */
	else if( fRef == f_min )
		phiRef -= (*phi)->data->data[0];
	/* phiRef is phase when f==fRef */
	else
	{
		VRef = pow(LAL_PI * LAL_G_SI*(m1+m2) * fRef, 1./3.) / LAL_C_SI;
		j = 0;
		do {
			idxRef = j;
			j++;
		} while ((*V)->data->data[j] <= VRef);
		phiRef -= (*phi)->data->data[idxRef];
	}
	for (j = 0; j < len; ++j)
		(*phi)->data->data[j] += phiRef;

	return (int)(*V)->data->length;
}


/**
 * Driver routine to compute the post-Newtonian inspiral waveform.
 *
 * This routine allows the user to specify different pN orders
 * for phasing calcuation vs. amplitude calculations.
 */
int XLALSimInspiralTaylorT3PNGenerator(
		REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
		REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
		REAL8 phiRef,                               /**< reference orbital phase (rad) */
		REAL8 v0,                                   /**< tail-term gauge choice (default = 1) */
		REAL8 deltaT,                               /**< sampling interval (s) */
		REAL8 m1,                                   /**< mass of companion 1 (kg) */
		REAL8 m2,                                   /**< mass of companion 2 (kg) */
		REAL8 f_min,                                /**< starting GW frequency (Hz) */
		REAL8 fRef,                                 /**< reference GW frequency (Hz) */
		REAL8 r,                                    /**< distance of source (m) */
		REAL8 i,                                    /**< inclination of source (rad) */
		REAL8 lambda1,                              /**< (tidal deformability of body 1)/(mass of body 1)^5 */
		REAL8 lambda2,                              /**< (tidal deformability of body 2)/(mass of body 2)^5 */
		LALSimInspiralInteraction interactionFlags, /**< flag to control spin and tidal effects */
		int amplitudeO,                             /**< twice post-Newtonian amplitude order */
		int phaseO                                  /**< twice post-Newtonian phase order */
		)
{
	REAL8TimeSeries *V;
	REAL8TimeSeries *phi;
	int status;
	int n;
	n = XLALSimInspiralTaylorT3PNEvolveOrbit(&V, &phi, phiRef, deltaT,
			m1, m2, f_min, fRef, lambda1, lambda2,
			interactionFlags, phaseO);
	if ( n < 0 )
		XLAL_ERROR(XLAL_EFUNC);
	status = XLALSimInspiralPNPolarizationWaveforms(hplus, hcross, V, phi,
			v0, m1, m2, r, i, amplitudeO);
	XLALDestroyREAL8TimeSeries(phi);
	XLALDestroyREAL8TimeSeries(V);
	if ( status < 0 )
		XLAL_ERROR(XLAL_EFUNC);
	return n;
}

/**
 * Driver routine to compute the -2 spin-weighted spherical harmonic mode
 * using TaylorT3 phasing.
 */
COMPLEX16TimeSeries *XLALSimInspiralTaylorT3PNModes(
		REAL8 phiRef,                               /**< reference orbital phase (rad) */
		REAL8 v0,                                   /**< tail-term gauge choice (default = 1) */
		REAL8 deltaT,                               /**< sampling interval (s) */
		REAL8 m1,                                   /**< mass of companion 1 (kg) */
		REAL8 m2,                                   /**< mass of companion 2 (kg) */
		REAL8 f_min,                                /**< starting GW frequency (Hz) */
		REAL8 fRef,                                 /**< reference GW frequency (Hz) */
		REAL8 r,                                    /**< distance of source (m) */
		REAL8 lambda1,                              /**< (tidal deformability of body 1)/(mass of body 1)^5 */
		REAL8 lambda2,                              /**< (tidal deformability of body 2)/(mass of body 2)^5 */
		LALSimInspiralInteraction interactionFlags, /**< flag to control spin and tidal effects */
		int amplitudeO,                             /**< twice post-Newtonian amplitude order */
		int phaseO,                                 /**< twice post-Newtonian phase order */
		int l,                                      /**< l index of mode */
		int m                                       /**< m index of mode */
		)
{
	COMPLEX16TimeSeries *hlm;
	/* The Schwarzschild ISCO frequency - for sanity checking fRef */
	REAL8 fISCO = pow(LAL_C_SI,3) / (pow(6.,3./2.)*LAL_PI*(m1+m2)*LAL_G_SI);

	/* Sanity check fRef value */
	if( fRef < 0. )
	{
		XLALPrintError("XLAL Error - %s: fRef = %f must be >= 0\n", 
				__func__, fRef);
		XLAL_ERROR_NULL(XLAL_EINVAL);
	}
	if( fRef != 0. && fRef < f_min )
	{
		XLALPrintError("XLAL Error - %s: fRef = %f must be > fStart = %f\n", 
				__func__, fRef, f_min);
		XLAL_ERROR_NULL(XLAL_EINVAL);
	}
	if( fRef >= fISCO )
	{
		XLALPrintError("XLAL Error - %s: fRef = %f must be < Schwar. ISCO=%f\n",
				__func__, fRef, fISCO);
		XLAL_ERROR_NULL(XLAL_EINVAL);
	}

	REAL8TimeSeries *V;
	REAL8TimeSeries *phi;
	int n;
	n = XLALSimInspiralTaylorT3PNEvolveOrbit(&V, &phi, phiRef, deltaT,
			m1, m2, f_min, fRef, lambda1, lambda2,
			interactionFlags, phaseO);
	if ( n < 0 )
		XLAL_ERROR_NULL(XLAL_EFUNC);
	hlm = XLALCreateSimInspiralPNModeCOMPLEX16TimeSeries(V, phi,
			v0, m1, m2, r, amplitudeO, l, m);
	XLALDestroyREAL8TimeSeries(phi);
	XLALDestroyREAL8TimeSeries(V);
	return hlm;
}

/**
 * Driver routine to compute the post-Newtonian inspiral waveform.
 *
 * This routine uses the same pN order for phasing and amplitude
 * (unless the order is -1 in which case the highest available
 * order is used for both of these -- which might not be the same).
 *
 * Constant log term in amplitude set to 1.  This is a gauge choice.
 */
int XLALSimInspiralTaylorT3PN(
		REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
		REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
		REAL8 phiRef,                               /**< reference orbital phase (rad) */
		REAL8 deltaT,                               /**< sampling interval (s) */
		REAL8 m1,                                   /**< mass of companion 1 (kg) */
		REAL8 m2,                                   /**< mass of companion 2 (kg) */
		REAL8 f_min,                                /**< starting GW frequency (Hz) */
		REAL8 fRef,                                 /**< reference GW frequency (Hz) */
		REAL8 r,                                    /**< distance of source (m) */
		REAL8 i,                                    /**< inclination of source (rad) */
		REAL8 lambda1,                              /**< (tidal deformability of body 1)/(mass of body 1)^5 */
		REAL8 lambda2,                              /**< (tidal deformability of body 2)/(mass of body 2)^5 */
		LALSimInspiralInteraction interactionFlags, /**< flag to control spin and tidal effects */
		int O                                       /**< twice post-Newtonian order */
		)
{
	/* set v0 to default value 1 */
	return XLALSimInspiralTaylorT3PNGenerator(hplus, hcross, phiRef, 1.0,
			deltaT, m1, m2, f_min, fRef, r, i, lambda1, lambda2,
			interactionFlags, O, O);
}


/**
 * Driver routine to compute the restricted post-Newtonian inspiral waveform.
 *
 * This routine computes the phasing to the specified order, but
 * only computes the amplitudes to the Newtonian (quadrupole) order.
 *
 * Constant log term in amplitude set to 1.  This is a gauge choice.
 */
int XLALSimInspiralTaylorT3PNRestricted(
		REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
		REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
		REAL8 phiRef,                               /**< reference orbital phase (rad) */
		REAL8 deltaT,                               /**< sampling interval (s) */
		REAL8 m1,                                   /**< mass of companion 1 (kg) */
		REAL8 m2,                                   /**< mass of companion 2 (kg) */
		REAL8 f_min,                                /**< starting GW frequency (Hz) */
		REAL8 fRef,                                 /**< reference GW frequency (Hz) */
		REAL8 r,                                    /**< distance of source (m)*/
		REAL8 i,                                    /**< inclination of source (rad) */
		REAL8 lambda1,                              /**< (tidal deformability of body 1)/(mass of body 1)^5 */
		REAL8 lambda2,                              /**< (tidal deformability of body 2)/(mass of body 2)^5 */
		LALSimInspiralInteraction interactionFlags, /**< flag to control spin and tidal effects */
		int O                                       /**< twice post-Newtonian phase order */
		)
{
	/* use Newtonian order for amplitude */
	/* set v0 to default value 1 */
	return XLALSimInspiralTaylorT3PNGenerator(hplus, hcross, phiRef, 1.0,
			deltaT, m1, m2, f_min, fRef, r, i, lambda1, lambda2,
			interactionFlags, 0, O);
}
