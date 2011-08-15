/*
 * Copyright (C) 2011 Drew Keppel, J. Creighton, S. Fairhurst, B. Krishnan, L. Santamaria, Stas Babak, David Churches, B.S. Sathyaprakash, Craig Robinson , Thomas Cokelaer, Duncan Brown, Gareth Jones, David McKechan, Riccardo Sturani, Laszlo Vereb
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
\author Keppel, D.
\file
\ingroup LALSimInspiralSetup_h

\brief Code to set up the structures and coefficients for generating the different PN waveforms

<tt>REAL8 XLALSimInspiralChirpLength()</tt>
<ul>
<li> \c m1: Mass of companion 1 in SI units. </li>
<li> \c m2: Mass of companion 2 in SI units. </li>
<li> \c f_min: Starting frequency. </li>
<li> \c O: PN order for computing sum of chirp times. </li>
</ul>


\heading{Description}


\heading{Algorithm}
FILL ME

\heading{Uses}
FILL ME

\heading{Notes}
FILL ME



*/

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <LALSimInspiralTaylorT3.h>
#include <LALSimInspiraldEnergyFlux.h>
#include <lal/LALSimInspiralPhasing3.h>
#include <lal/LALSimInspiralFrequency3.h>


/**
 * Returns the sum of chirp times to a specified order.
 *
 * Computes the sum of the chirp times to a specified order. Inputs given in SI
 * units.
 */
REAL8 XLALSimInspiralChirpLength(
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
int XLALSimInspiralTaylorT3Setup(
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
  REAL8 vpole, eta, lso, vlso, vn, tofv;
  REAL8 oneby6 = 1./6.;
  TofVIn in1;
  void *in2;


  vpole = 0.0;
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
  vlso = 0;

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

