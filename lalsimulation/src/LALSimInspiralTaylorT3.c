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

#include <lal/LALSimInspiral.h>
#define LAL_USE_COMPLEX_SHORT_MACROS
#include <lal/FindRoot.h>
#include <lal/LALComplex.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/LALSimInspiraldEnergyFlux.h>

#include "check_series_macros.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

NRCSID(LALSIMINSPIRALTAYLORT3C, "$Id$");

typedef struct
tagexpnCoeffsTaylorT3 {
   int ieta;
   /* Taylor expansion coefficents in phi(t)*/
   REAL8 ptaN, pta2, pta3, pta4, pta5, pta6, pta7, ptl6;
   /* Taylor expansion coefficents in f(t)*/
   REAL8 ftaN, fta2, fta3, fta4, fta5, fta6, fta7, ftl6;

   /* sampling interval*/
   REAL8 samplinginterval;
   /* symmetric mass ratio, total mass, component masses*/
   REAL8 eta, totalmass, m1, m2;
   /* unknown 3PN parameters, euler constant*/
   REAL8 lambda, theta, EulerC;

   /* initial and final values of frequency, time, velocity; lso
    values of velocity and frequency; final phase.*/
   REAL8 f0, fn, t0, tn, v0, vn, vf, vlso, flso, phiC;

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
  static const char *func = "XLALInspiralFrequency3Wrapper";

  FreqInFromChirptime *in;
  REAL8 freq, f;

  in = (FreqInFromChirptime *) pars;
  freq = in->func(tC, &(in->ak));
  if (XLAL_IS_REAL8_FAIL_NAN(freq))
    XLAL_ERROR_REAL8(func, XLAL_EFUNC);
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
    XLAL_ERROR_REAL8(__func__, XLAL_EFAULT);

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
  REAL8 theta,theta2,theta3;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(__func__, XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;

  frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2);

  return frequency;
}


static REAL8
XLALSimInspiralFrequency3_3PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(__func__, XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;

  frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3);

  return frequency;
}


static REAL8
XLALSimInspiralFrequency3_4PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3,theta4;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(__func__, XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;

  frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3
             + ak->fta4*theta4);

  return frequency;
}


static REAL8
XLALSimInspiralFrequency3_5PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(__func__, XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;

  frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3
             + ak->fta4*theta4
             + ak->fta5*theta5);

  return frequency;
}


static REAL8
XLALSimInspiralFrequency3_6PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5,theta6;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(__func__, XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta6 = theta5*theta;

  frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3
             + ak->fta4*theta4
             + ak->fta5*theta5
             + (ak->fta6 + ak->ftl6*log(td))*theta6);

  return frequency;
}


static REAL8
XLALSimInspiralFrequency3_7PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5,theta6,theta7;
  REAL8 frequency;

  if (ak == NULL)
    XLAL_ERROR_REAL8(__func__, XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta6 = theta5*theta;
  theta7 = theta6*theta;

  frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3
             + ak->fta4*theta4
             + ak->fta5*theta5
             + (ak->fta6 + ak->ftl6*log(td))*theta6
             + ak->fta7*theta7);

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
    XLAL_ERROR_REAL8(__func__, XLAL_EFAULT);

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
  REAL8 theta,theta2,theta5;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(__func__, XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta5 = theta2*theta2*theta;

  phase = (ak->ptaN/theta5) * (1.
         + ak->pta2*theta2);

  return phase;
}


static REAL8
XLALSimInspiralPhasing3_3PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3,theta5;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(__func__, XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta5 = theta2*theta3;

  phase = (ak->ptaN/theta5) * (1.
         + ak->pta2*theta2
         + ak->pta3*theta3);

  return phase;
}


static REAL8
XLALSimInspiralPhasing3_4PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(__func__, XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;

  phase = (ak->ptaN/theta5) * (1.
         + ak->pta2*theta2
         + ak->pta3*theta3
         + ak->pta4*theta4);

  return phase;
}


static REAL8
XLALSimInspiralPhasing3_5PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(__func__, XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;

  phase = (ak->ptaN/theta5) * (1.
         + ak->pta2*theta2
         + ak->pta3*theta3
         + ak->pta4*theta4
         + ak->pta5 * log(td/ak->tn) * theta5);

  return phase;
}


static REAL8
XLALSimInspiralPhasing3_6PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5,theta6;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(__func__, XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta6 = theta5*theta;

  phase = (ak->ptaN/theta5) * (1.
         + ak->pta2*theta2
         + ak->pta3*theta3
         + ak->pta4*theta4
         + ak->pta5*log(td/ak->tn)*theta5
         +(ak->ptl6*log(td/256.) + ak->pta6)*theta6);

  return phase;
}


static REAL8
XLALSimInspiralPhasing3_7PN (
   REAL8       td,
   expnCoeffsTaylorT3 *ak
   )
{
  REAL8 theta,theta2,theta3,theta4,theta5,theta6,theta7;
  REAL8 phase;

  if (ak == NULL)
    XLAL_ERROR_REAL8(__func__, XLAL_EFAULT);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta6 = theta5*theta;
  theta7 = theta6*theta;

  phase = (ak->ptaN/theta5) * (1.
         + ak->pta2*theta2
         + ak->pta3*theta3
         + ak->pta4*theta4
         + ak->pta5*log(td/ak->tn)*theta5
         +(ak->ptl6*log(td/256.) + ak->pta6)*theta6
         + ak->pta7*theta7);

  return phase;
}


/**
 * Returns the sum of chirp times to a specified order.
 *
 * Computes the sum of the chirp times to a specified order. Inputs given in SI
 * units.
 */
static REAL8 XLALSimInspiralChirpLength(
		REAL8 m1,		/**< mass of companion 1 */
		REAL8 m2,		/**< mass of companion 2 */
		REAL8 f_min,		/**< start frequency */
		int O			/**< twice post-Newtonian order */
		)
{
	REAL8 v, tN, t0, t2, t3 ,t4, t5, t6, t7, tC;
	REAL8 piFl = LAL_PI * f_min;
	REAL8 m = m1 + m2;
	REAL8 eta = m1 * m2 / m / m;
	m *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m from kilograms to seconds */

	REAL8 oneby3 = 1.L/3.L;
	REAL8 twoby3 = 2.L/3.L;
	REAL8 fourby3 = 4.L/3.L;
	REAL8 fiveby3 = 5.L/3.L;
	REAL8 eightby3 = 8.L/3.L;

	/* Should use the coefficients from LALInspiraSetup.c to avoid errors.
	 */
	v = pow(piFl * m, 1.L/3.L);
	tN = 5.L/256.L / eta * m / pow(v,8.L);

	t0 = 5.0L/(256.0L*eta*pow(m,fiveby3)*pow(piFl,eightby3));
	t2 = (3715.0L + (4620.0L*eta))/(64512.0*eta*m*pow(piFl,2.0));
	t3 = LAL_PI/(8.0*eta*pow(m,twoby3)*pow(piFl,fiveby3));
	t4 = (5.0/(128.0*eta*pow(m,oneby3)*pow(piFl,fourby3)))
		* (3058673./1016064. + 5429.*eta/1008. +617.*eta*eta/144.);
	t5 = -5.*(7729./252. - 13./3.*eta)/(256.*eta*f_min);
	/* This is a ddraft. t6 and t7 need to be checked propely. For now set to zero for all orders below. */
	t6 = -10052469856691./23471078400. + 128./3.*LAL_PI*LAL_PI
		+(15335597827.L/15240960.L-451.L/12.L*LAL_PI*LAL_PI+352./3.*-11831.L/9240.L-2464.L/9.L*-1987.L/3080.L)*eta
		+6848.L/105.L* LAL_GAMMA
		-15211.L/1728.L*eta*eta+25565.L/1296.L*eta*eta*eta;
	t6 = tN * (t6  + 6848.L/105.L*log(4.*v)) * pow(v,6);
	t7 = (-15419335.L/127008.L-75703.L/756.L*eta+14809.L/378.L*eta*eta) * LAL_PI * tN * pow(v,7);

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
		case 6:
		case 7:
			t6 = 0.;
			t7 = 0.;
			break;
		case 8:
			XLALPrintError("XLAL Error - %s: Not supported for requested PN order\n", __func__);
			XLAL_ERROR_REAL8(__func__, XLAL_EINVAL);
			break;
		default:
			XLALPrintError("XLAL Error - %s: Unknown PN order in switch\n", __func__);
			XLAL_ERROR_REAL8(__func__, XLAL_EINVAL);
	}

	tC = t0 + t2 - t3 + t4 - t5 + t6 - t7;

	return tC;
}


/**
 * Set up the expnCoeffsTaylorT3 and expnFuncTaylorT3 structures for
 * generating a TaylorT3 waveform.
 *
 * Inputs given in SI units.
 */
static int XLALSimInspiralTaylorT3Setup(
		expnCoeffsTaylorT3 *ak,	/**< coefficients for TaylorT3 evolution [modified] */
		expnFuncTaylorT3 *f,	/**< functions for TaylorT3 evolution [modified] */
	       	REAL8 deltaT,		/**< sampling interval */
		REAL8 m1,		/**< mass of companion 1 */
		REAL8 m2,		/**< mass of companion 2 */
		REAL8 f_min,		/**< start frequency */
		int O			/**< twice post-Newtonian order */
		)
{
  expnCoeffsdEnergyFlux akEF;
  REAL8 eta, lso, vn, tofv,vlso;//vpole - set but not used
  REAL8 oneby6 = 1./6.;
  TofVIn in1;
  void *in2;


  //vpole = 0.0; - set but not used
  ak->EulerC = LAL_GAMMA;
  ak->lambda = -1987./3080.;
  ak->theta = -11831./9240.;
  ak->t0 = 0;
  ak->m1 = m1;
  ak->m2 = m2;
  ak->totalmass = ak->m1 + ak->m2;
  eta = ak->eta = m1 * m2 / (ak->totalmass * ak->totalmass);
  ak->totalmass = ak->totalmass * LAL_MTSUN_SI;
  ak->f0 = f_min;
  ak->samplinginterval = deltaT;
  ak->fn = 1. / (2. * ak->samplinginterval);
  ak->v0 = cbrt(LAL_PI * ak->totalmass * f_min);

  ak->ptaN = -2./eta;
  ak->pta2 = 3715./8064. + 55.*eta/96.;
  ak->pta3 = -0.75*LAL_PI;
  ak->pta4 = 9.275495/14.450688 + 2.84875*eta/2.58048
           + 1855.*eta*eta/2048.;
  ak->pta5 = -(3.8645/17.2032 - 65./2048.*eta) * LAL_PI;
  ak->pta6 =  (83.1032450749357/5.7682522275840 - 53./40.*LAL_PI*LAL_PI
           - 107./56. * ak->EulerC + (-123.292747421/4.161798144
           + 2.255/2.048 *LAL_PI*LAL_PI + 385./48. * ak->lambda
           - 55./16.* ak->theta) * eta + 1.54565/18.35008 * eta*eta
           - 1.179625/1.769472 * eta*eta*eta);

  ak->pta7 =  (1.88516689/1.73408256 + 488825./516096. * eta
           - 141769./516096. * eta*eta) * LAL_PI;
  ak->ptl6 = 107./448.;

  ak->ftaN = 1./(8.*LAL_PI*ak->totalmass);
  ak->fta2 = 743./2688.+(11.*eta)/32.;
  ak->fta3 = -0.3*LAL_PI;
  ak->fta4 = 1.855099/14.450688 +  5.6975*eta/25.8048
           + 3.71*eta*eta/20.48;
  ak->fta5 = -(7.729/21.504 - 13./256.*eta) * LAL_PI;
  ak->fta6 = (-7.20817631400877/2.88412611379200 + (53./200.) * LAL_PI*LAL_PI
           + 1.07/2.80 * ak->EulerC + 1.07/2.80 * log(2.)
           + (1.23292747421/.20808990720 - 4.51/20.48 * LAL_PI*LAL_PI
           - 77./48. * ak->lambda + 11./16. * ak->theta) * eta
           - 3.0913/183.5008 * eta*eta
           + 2.35925/17.69472 * eta*eta*eta);

  ak->fta7 = (-1.88516689/4.33520640 - 97765./258048. * eta
           + 141769./1290240. * eta*eta)*LAL_PI;
  ak->ftl6 = -107./2240.;


/* Taylor coefficients of E(x) */
  akEF.ETaN = -eta/2.;
  akEF.ETa1 = -(9. + eta)/12.;
  akEF.ETa2 = -(27. - 19*eta + eta*eta/3.)/8.;
  akEF.ETa3 = -675./64. + (209323./4032. - 205.*LAL_PI*LAL_PI/96.
           - 110./9. * ak->lambda)*eta
           - 155./96. * eta*eta - 35./5184. * eta*eta*eta;

/* Taylor coefficients of dE(v)/dv. (NOTE v and NOT x) */
  akEF.dETaN = -eta;
  akEF.dETa1 = 2.*akEF.ETa1;
  akEF.dETa2 = 3.*akEF.ETa2;
  akEF.dETa3 = 4.*akEF.ETa3;

/* Taylor coefficients of flux. */
  akEF.FTaN = 32.*eta*eta/5.;
  akEF.FTa1 = 0.;
  akEF.FTa2 = -1247./336.-35.*eta/12.;
  akEF.FTa3 = 4.*LAL_PI;
  akEF.FTa4 = -44711./9072.+9271.*eta/504.+65.*eta*eta/18.;
  akEF.FTa5 = -(8191./672.+583./24.*eta)*LAL_PI;
  akEF.FTl6 = -1712./105.;
  akEF.FTa6 = 6643739519./69854400. + 16.*LAL_PI*LAL_PI/3. + akEF.FTl6 * log (4.L)
           - 1712./105.*ak->EulerC+ (-134543./7776. + 41.*LAL_PI*LAL_PI/48.) * eta
           - 94403./3024. * eta*eta - 775./324. * eta*eta*eta;
  akEF.FTa7 = LAL_PI * (-16285./504. + 214745./1728. * eta
              + 193385./3024.* eta*eta);
  akEF.FTa8 = - 117.5043907226773;
  akEF.FTl8 =   52.74308390022676;

  lso = sqrt(oneby6);
  vlso = 0; //- set but not used

/* Location of the 0PN and 1PN T- and P-approximant last stable orbit: */
  akEF.vlsoT0 = lso;
  akEF.vlsoP0 = lso;
  akEF.vlsoP2 = lso;
/*
  vlsoT2 =  6./(9.+eta);
  This correct value makes vlso too large for vlsoT2 hence use 1/sqrt(6)
*/
  akEF.vlsoT2 = lso;

  switch (O)
  {
     case 0:
           vn = akEF.vlso = vlso = akEF.vlsoT0;
           in1.dEnergy = dEt0;
           in1.flux = Ft0;
           f->phasing3 = &XLALSimInspiralPhasing3_0PN;
           f->frequency3 = &XLALSimInspiralFrequency3_0PN;
           break;
     case 1:
       XLALPrintError("XLAL Error - %s: PN approximant not supported for requested PN order\n", __func__);
       XLAL_ERROR(__func__, XLAL_EINVAL);
       break;
     case 2:
           vn = akEF.vlso = vlso = akEF.vlsoT2;
           in1.dEnergy = dEt2;
           in1.flux = Ft2;
           f->phasing3 = &XLALSimInspiralPhasing3_2PN;
           f->frequency3 = &XLALSimInspiralFrequency3_2PN;
           break;
     case 3:
           vn = akEF.vlso = vlso = akEF.vlsoT2;
           in1.dEnergy = dEt2;
           in1.flux = Ft3;
           f->phasing3 = &XLALSimInspiralPhasing3_3PN;
           f->frequency3 = &XLALSimInspiralFrequency3_3PN;
           break;
     case 4:
/*
   The value vlsoT4 is too large and doesn't work sometimes;
   so we use vlsoT2.
*/
           vn = akEF.vlso = vlso = akEF.vlsoT2;
           in1.dEnergy = dEt4;
           in1.flux = Ft4;
           f->phasing3 = &XLALSimInspiralPhasing3_4PN;
           f->frequency3 = &XLALSimInspiralFrequency3_4PN;
           break;
     case 5:
/*
   The value vlsoT4 is too large and doesn't work with 2.5 PN
   Taylor approximant; so we use vlsoT2.
*/
           vn = akEF.vlso = vlso = akEF.vlsoT2;
           in1.dEnergy = dEt4;
           in1.flux = Ft5;
           f->phasing3 = &XLALSimInspiralPhasing3_5PN;
           f->frequency3 = &XLALSimInspiralFrequency3_5PN;
           break;
     case 6:
/*
   vlsoT6 is as yet undetermined and vlsoT4 is too large in
   certain cases (TaylorT2 crashes for (1.4,10)); using vlsoT2;
*/
           vn = akEF.vlso = vlso = akEF.vlsoT2;
           in1.dEnergy = dEt6;
           in1.flux = Ft6;
           f->phasing3 = &XLALSimInspiralPhasing3_6PN;
           f->frequency3 = &XLALSimInspiralFrequency3_6PN;
           break;
     case 7:
           vn = akEF.vlso = vlso = akEF.vlsoT2;
           in1.dEnergy = dEt6;
           in1.flux = Ft7;
           f->phasing3 = &XLALSimInspiralPhasing3_7PN;
           f->frequency3 = &XLALSimInspiralFrequency3_7PN;
           break;
     case 8:
           XLALPrintError("XLAL Error - %s: PN approximant not supported for requested PN order\n", __func__);
           XLAL_ERROR(__func__, XLAL_EINVAL);
           break;
     default:
        XLALPrintError("XLAL Error - %s: Unknown PN order in switch\n", __func__);
        XLAL_ERROR(__func__, XLAL_EINVAL);
  }

  if (ak->fn)
  {
    vn = cbrt(LAL_PI * ak->totalmass * ak->fn);
    vn = (vn < vlso) ? vn :  vlso;
  }

  in1.t=0.0;
  in1.v0=ak->v0;
  in1.t0=ak->t0;
  in1.vlso=akEF.vlso;
  in1.totalmass = ak->totalmass;
  in1.coeffs = &akEF;

  in2 = (void *) &in1;

  tofv = XLALSimInspiralTofV(vn, in2);
  if (XLAL_IS_REAL8_FAIL_NAN(tofv))
    XLAL_ERROR(__func__, XLAL_EFUNC);

  ak->tn = -tofv - ak->samplinginterval;

  return 0;
}


/**
 * Computes a post-Newtonian orbit using the Taylor T3 method.
 */
int XLALSimInspiralTaylorT3PNEvolveOrbit(
		REAL8TimeSeries **x,   /**< post-Newtonian parameter [returned] */
	       	REAL8TimeSeries **phi, /**< orbital phase [returned] */
	       	LIGOTimeGPS *tc,       /**< coalescence time */
	       	REAL8 phi0,            /**< initial phase */
	       	REAL8 deltaT,          /**< sampling interval */
		REAL8 m1,              /**< mass of companion 1 */
		REAL8 m2,              /**< mass of companion 2 */
		REAL8 f_min,           /**< start frequency */
		int O                  /**< twice post-Newtonian order */
		)
{
	const UINT4 blocklen = 1024;
	const REAL8 xisco = 1.0 / 6.0;
	REAL8 m = m1 + m2;
	REAL8 nu = m1 * m2 / m / m;
	m *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m from kilograms to seconds */
	REAL8 mass1 = m1 / LAL_MSUN_SI; /* convert m1 from kilograms to solar masses */
	REAL8 mass2 = m2 / LAL_MSUN_SI; /* convert m2 from kilograms to solar masses */
	REAL8 tmptC, tC, c1, xmin, xmax, xacc, v, v2, phase, fOld, t, td, temp, tempMin = 0, tempMax = 0;
	REAL8 (*freqfunction)(REAL8, void *);
	UINT4 j;
	REAL8 f;
	void *pars;

	expnFuncTaylorT3 expnfunc;
	expnCoeffsTaylorT3 ak;
	FreqInFromChirptime timeIn;

	/* allocate memory */

	*x = XLALCreateREAL8TimeSeries("ORBITAL_FREQUENCY_PARAMETER", tc, 0.0, deltaT, &lalDimensionlessUnit,
		blocklen);
	*phi = XLALCreateREAL8TimeSeries("ORBITAL_PHASE", tc, 0.0, deltaT, &lalDimensionlessUnit, blocklen);
	if (!x || !phi)
		XLAL_ERROR(__func__, XLAL_EFUNC);


	/* initialize expnCoeffsTaylorT3 and expnFuncTaylorT3 structures */
	if (XLALSimInspiralTaylorT3Setup(&ak, &expnfunc, deltaT, mass1, mass2, f_min, O))
		XLAL_ERROR(__func__, XLAL_EFUNC);

	tC = XLALSimInspiralChirpLength(m1, m2, f_min, O);
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
			XLAL_ERROR(__func__, XLAL_EFUNC);
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
			XLAL_ERROR(__func__, XLAL_EFUNC);
	}
	else{
		XLALPrintError("Can't find good bracket for BisectionFindRoot");
		XLAL_ERROR(__func__, XLAL_EMAXITER);
	}

	tC /= c1;

	/* start waveform generation */

	t = 0.;
	td = c1 * (tC - t);
	phase = expnfunc.phasing3(td, &ak);
	if (XLAL_IS_REAL8_FAIL_NAN(phase))
		XLAL_ERROR(__func__, XLAL_EFUNC);
	f = expnfunc.frequency3(td, &ak);
	if (XLAL_IS_REAL8_FAIL_NAN(f))
		XLAL_ERROR(__func__, XLAL_EFUNC);

	v = cbrt(f * LAL_PI * m);
	v2 = v * v;
	(*x)->data->data[0] = v2;
	(*phi)->data->data[0] = phase / 2.;

	j = 0;
	while (1) {

		/* make one step */

		j++;
		fOld = f;
		t = j * deltaT;
		td = c1 * (tC - t);
		phase = expnfunc.phasing3(td, &ak);
		if (XLAL_IS_REAL8_FAIL_NAN(phase))
			XLAL_ERROR(__func__, XLAL_EFUNC);
		f = expnfunc.frequency3(td, &ak);
		if (XLAL_IS_REAL8_FAIL_NAN(f))
			XLAL_ERROR(__func__, XLAL_EFUNC);
		v = cbrt(f * LAL_PI * m);
		v2 = v * v;

		/* check termination conditions */

		if (t >= tC) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated at coalesence time\n", __func__);
			break;
		}
		if (v2 >= xisco) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated at ISCO\n", __func__);
			break;
		}
		if (f <= fOld) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated when frequency stalled\n", __func__);
			break;
		}
	
		/* save current values in vectors but first make sure we don't write past end of vectors */

		if ( j >= (*x)->data->length ) {
			if ( ! XLALResizeREAL8TimeSeries(*x, 0, (*x)->data->length + blocklen) )
				XLAL_ERROR(__func__, XLAL_EFUNC);
			if ( ! XLALResizeREAL8TimeSeries(*phi, 0, (*phi)->data->length + blocklen) )
				XLAL_ERROR(__func__, XLAL_EFUNC);
		}
		(*x)->data->data[j] = v2;
		(*phi)->data->data[j] = phase / 2.;
	}

	/* make the correct length */

	if ( ! XLALResizeREAL8TimeSeries(*x, 0, j) )
		XLAL_ERROR(__func__, XLAL_EFUNC);
	if ( ! XLALResizeREAL8TimeSeries(*phi, 0, j) )
		XLAL_ERROR(__func__, XLAL_EFUNC);

	/* adjust to correct tc and phic */

	XLALGPSAdd(&(*phi)->epoch, -1.0*j*deltaT);
	XLALGPSAdd(&(*x)->epoch, -1.0*j*deltaT);

	/* adjust so initial phase is phi0 */

	phi0 -= (*phi)->data->data[0];
	for (j = 0; j < (*phi)->data->length; ++j)
		(*phi)->data->data[j] += phi0;

	return (int)(*x)->data->length;
}


/**
 * Driver routine to compute the post-Newtonian inspiral waveform.
 *
 * This routine allows the user to specify different pN orders
 * for phasing calcuation vs. amplitude calculations.
 */
int XLALSimInspiralTaylorT3PNGenerator(
		REAL8TimeSeries **hplus,  /**< +-polarization waveform */
	       	REAL8TimeSeries **hcross, /**< x-polarization waveform */
	       	LIGOTimeGPS *tc,          /**< coalescence time */
	       	REAL8 phic,               /**< coalescence phase */
	       	REAL8 x0,                 /**< tail-term gauge choice thing (if you don't know, just set it to zero) */
	       	REAL8 deltaT,             /**< sampling interval */
	       	REAL8 m1,                 /**< mass of companion 1 */
	       	REAL8 m2,                 /**< mass of companion 2 */
	       	REAL8 f_min,              /**< start frequency */
	       	REAL8 r,                  /**< distance of source */
	       	REAL8 i,                  /**< inclination of source (rad) */
	       	int amplitudeO,           /**< twice post-Newtonian amplitude order */
	       	int phaseO                /**< twice post-Newtonian phase order */
		)
{
	static const char *func = "XLALSimInspiralPNGenerator";
	REAL8TimeSeries *x;
	REAL8TimeSeries *phi;
	int status;
	int n;
	n = XLALSimInspiralTaylorT3PNEvolveOrbit(&x, &phi, tc, phic, deltaT, m1, m2, f_min, phaseO);
	if ( n < 0 )
		XLAL_ERROR(func, XLAL_EFUNC);
	status = XLALSimInspiralPNPolarizationWaveforms(hplus, hcross, x, phi, x0, m1, m2, r, i, amplitudeO);
	XLALDestroyREAL8TimeSeries(phi);
	XLALDestroyREAL8TimeSeries(x);
	if ( status < 0 )
		XLAL_ERROR(func, XLAL_EFUNC);
	return n;
}


/**
 * Driver routine to compute the post-Newtonian inspiral waveform.
 *
 * This routine uses the same pN order for phasing and amplitude
 * (unless the order is -1 in which case the highest available
 * order is used for both of these -- which might not be the same).
 *
 * Log terms in amplitudes are ignored.  This is a gauge choice.
 */
int XLALSimInspiralTaylorT3PN(
		REAL8TimeSeries **hplus,  /**< +-polarization waveform */
	       	REAL8TimeSeries **hcross, /**< x-polarization waveform */
	       	LIGOTimeGPS *tc,          /**< coalescence time */
	       	REAL8 phic,               /**< coalescence phase */
	       	REAL8 deltaT,             /**< sampling interval */
	       	REAL8 m1,                 /**< mass of companion 1 */
	       	REAL8 m2,                 /**< mass of companion 2 */
	       	REAL8 f_min,              /**< start frequency */
	       	REAL8 r,                  /**< distance of source */
	       	REAL8 i,                  /**< inclination of source (rad) */
	       	int O                     /**< twice post-Newtonian order */
		)
{
	/* set x0=0 to ignore log terms */
	return XLALSimInspiralTaylorT3PNGenerator(hplus, hcross, tc, phic, 0.0, deltaT, m1, m2, f_min, r, i, O, O);
}


/**
 * Driver routine to compute the restricted post-Newtonian inspiral waveform.
 *
 * This routine computes the phasing to the specified order, but
 * only computes the amplitudes to the Newtonian (quadrupole) order.
 *
 * Log terms in amplitudes are ignored.  This is a gauge choice.
 */
int XLALSimInspiralTaylorT3PNRestricted(
		REAL8TimeSeries **hplus,  /**< +-polarization waveform */
	       	REAL8TimeSeries **hcross, /**< x-polarization waveform */
	       	LIGOTimeGPS *tc,          /**< coalescence time */
	       	REAL8 phic,               /**< coalescence phase */
	       	REAL8 deltaT,             /**< sampling interval */
	       	REAL8 m1,                 /**< mass of companion 1 */
	       	REAL8 m2,                 /**< mass of companion 2 */
	       	REAL8 f_min,              /**< start frequency */
	       	REAL8 r,                  /**< distance of source */
	       	REAL8 i,                  /**< inclination of source (rad) */
	       	int O                     /**< twice post-Newtonian phase order */
		)
{
	/* use Newtonian order for amplitude */
	/* set x0=0 to ignore log terms */
	return XLALSimInspiralTaylorT3PNGenerator(hplus, hcross, tc, phic, 0.0, deltaT, m1, m2, f_min, r, i, 0, O);
}
