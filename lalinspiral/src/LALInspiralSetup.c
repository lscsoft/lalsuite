/*
*  Copyright (C) 2007 David Churches, Drew Keppel, Duncan Brown, Gareth Jones, Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
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
 * \author Sathyaprakash, B. S.
 * \file
 * \ingroup LALInspiral_h
 *
 * \brief Module to generate all the Taylor and Pade coefficients needed in
 * waveform generation.
 *
 * ### Prototypes ###
 *
 * <tt>XLALInspiralSetup()</tt>
 * <ul>
 * <li> \c ak: Output containing PN expansion coefficients of various physical
 * quantities such as energy, flux, frequency, phase and timing.</li>
 * <li> \c params: Input containing binary chirp parameters.</li>
 * </ul>
 *
 * ### Description ###
 *
 * Module to generate all the coefficiants needed in the Taylor and Pade expressions
 * for the energy and flux functions \f$E^{\prime}(v)\f$ and \f$\mathcal{F}(v)\f$.
 * These are used to solve the gravitational wave phasing formula.
 * The coefficients are used by the function \c LALInspiralChooseModel to define
 * the energy and flux functions by accessing the structure \c ak and are tabulated
 * in the two Tables\tableref{table_energy} and\tableref{table_flux}.
 *
 * ### Algorithm ###
 *
 * None.
 *
 * ### Uses ###
 *
 * None.
 *
 * ### Notes ###
 *
 */



#include <math.h>
#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>

/* static void LALPadeCoeffs7(int n, double *padecoeffs, double *taylorcoeffs); */

void
LALInspiralSetup (
   LALStatus        *status,
   expnCoeffs       *ak,
   InspiralTemplate *params
   )
{
   XLALPrintDeprecationWarning("LALInspiralSetup", "XLALInspiralSetup");

   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   if ( XLALInspiralSetup(ak, params) == XLAL_FAILURE )
   {
      ABORTXLAL(status);
   }
   DETATCHSTATUSPTR(status);
   RETURN (status);
}

int
XLALInspiralSetup (
   expnCoeffs       *ak,
   InspiralTemplate *params
   )
{
   INT4 ieta;
   /*INT4  pnorder=7;
    * */
   REAL8 lso, eta, vpole, vlso;
   REAL8 a1, a2, a3, a4, a5, a6, a7, a8;
   REAL8 c1, c2, c3, c4, c5, c6, c7, c8;
   REAL8 a12, a22, a32, a42, a52, a62, a72, a23, a33, a43, a53, a34, a44;
   REAL8 oneby6=1.0/6.0;
   REAL8 beta = 0.L;
   REAL8 sigma = 0.L;
   REAL8 chi1, chi2;

   if (ak == NULL)
      XLAL_ERROR(XLAL_EFAULT);
   if (params == NULL)
      XLAL_ERROR(XLAL_EFAULT);
   if (params->mass1 <= 0)
      XLAL_ERROR(XLAL_EDOM);
   if (params->mass2 <= 0)
      XLAL_ERROR(XLAL_EDOM);
   if (params->fLower <= 0)
      XLAL_ERROR(XLAL_EDOM);
   if (params->fCutoff <= 0)
      XLAL_ERROR(XLAL_EDOM);
   if (params->tSampling <= 0)
      XLAL_ERROR(XLAL_EDOM);
   if (params->tSampling <= 2*params->fCutoff)
      XLAL_ERROR(XLAL_EDOM);

   vpole = 0.0;
   ak->omegaS = params->OmegaS;
   ak->zeta2  = params->Zeta2;
   ak->ieta = params->ieta;
   ak->EulerC = LAL_GAMMA;
   ak->lambda = -1987./3080.;
   ak->theta = -11831./9240.;
   ak->t0 = params->startTime;
   ak->m1 = params->mass1;
   ak->m2 = params->mass2;
   ak->f0 = params->fLower;
   ak->fn = params->fCutoff;
   ak->samplingrate = params->tSampling;
   ak->samplinginterval = 1./ak->samplingrate;

   lso = sqrt(oneby6);
/*
   ieta determines the nature of the waveforms:
   ieta=0 testmass waveforms
   ieta=1 comparable mass waveforms.
*/
   ieta = ak->ieta;

/* Compute the total mass and eta from m1 and m2 */
   ak->totalmass = (ak->m1 + ak->m2);
   ak->eta = (ak->m1*ak->m2) / (ak->totalmass*ak->totalmass);
   eta = ak->eta;
   ak->totalmass = ak->totalmass * LAL_MTSUN_SI;


/* Aligned spin corrections (Poisson and Will PRD 52 848 (1995))
   Use the z components of the spins */
  chi1 = params->spin1[2];
  chi2 = params->spin2[2];
  if (eta <= 0.25L)
  {
    /* Chi1 is spin on larger mass
       m1 = mtot * (1 + sqrt(1 - 4 eta)) / 2
       m2 = mtot * (1 - sqrt(1 - 4 eta)) / 2 */
    beta = ((113.L - 76.L * ieta * eta) * (chi1 + chi2) / 24.L)
         + (113.L * sqrt(1.L - 4.L * ieta * eta) * (chi1 - chi2) / 24.L);
    sigma = 474.L * ieta * eta * chi1 * chi2 / 48.L;
  }


/* Set initial velocity according to initial frequency */

   ak->v0 = pow (LAL_PI * ak->totalmass * ak->f0, (1./3.));

/* Taylor coefficients of E(x) */
   ak->ETaN = -eta/2.;
   ak->ETa1 = -(9. + ieta*eta)/12.;
   ak->ETa2 = -(27. - 19*ieta*eta + ieta*eta*eta/3.)/8.;
   ak->ETa3 = -675./64. + (209323./4032. - 205.*LAL_PI*LAL_PI/96.
            - 110./9. * ak->lambda)*ieta*eta
            - 155./96. * ieta*eta*eta - 35./5184. * ieta*eta*eta*eta;

/* Taylor coefficients of e(x) */
   ak->eTaN = -1.;
   ak->eTa1 = -1.-ieta*eta/3.;
   ak->eTa2 = -3.+35.*ieta*eta/12.;
   ak->eTa3 = -9.+(26189./504. - 205./96.*LAL_PI*LAL_PI -110./9.*ak->lambda)*ieta*eta
              - 103./36.*ieta*eta*eta + ieta*eta*eta*eta/81.;

/* Taylor coefficients of dE(v)/dv. (NOTE v and NOT x) */
   ak->dETaN = -eta;
   ak->dETa1 = 2.*ak->ETa1;
   ak->dETa2 = 3.*ak->ETa2;
   ak->dETa3 = 4.*ak->ETa3;

/* Pade coefficients of e(x)-function. */
   ak->ePaN = -1.;
   c1 = ak->ePa1 = 1.+ieta*eta/3.;
   c2 = ak->ePa2 = -(144. - 81.*ieta*eta + 4.*ieta*eta*eta) / (36.+12*ieta*eta);
   c3 = ak->ePa3 = (ak->eTa1*ak->eTa3 - ak->eTa2*ak->eTa2)
                 / (ak->eTa1*(ak->eTa1*ak->eTa1 - ak->eTa2));

/* Location of the 0PN and 1PN T- and P-approximant last stable orbit: */
   ak->vlsoT0 = lso;
   ak->vlsoP0 = lso;
   ak->vlsoP2 = lso;
/*
   vlsoT2 =  6./(9.+ieta*eta);
   This correct value makes vlso too large for vlsoT2 hence use 1/sqrt(6)
*/
   ak->vlsoT2 = lso;

/* vlsoT4 is also too large for certain T-approximants: most notably
   TaylorT3, TaylorT2(3PN (1.4,10) crashes)
*/
   ak->vlsoT4 = pow(-ak->ETa1 + pow(ak->ETa1*ak->ETa1 - 3*ak->ETa2,0.5)/(3*ak->ETa2), 0.5);
   ak->vlsoP4 = pow((-1.+pow(-ak->ePa1/ak->ePa2,0.5))/(ak->ePa1 + ak->ePa2), 0.5);
   ak->vpoleP4 = pow(4.*(3.+ieta*eta)/(36.-35.*ieta*eta), 0.5);
/* THE VALUE NEEDS TO BE CHANGED AFTER THE CALCULATION OF CORRECT LSOs */
   ak->vlsoT6 = ak->vlsoT4;
/*
   vlsoP6 comes from solving a quadratic equation; the plus-root
   is negative and can't be used; the postive root is:
*/
   ak->vlsoP6 = sqrt(0.5/(c2*c2+c3*c3+c2*c1+2.*c2*c3)
              *(-2.*c3 - 2.*c2 - 2.*sqrt(-c2*c1)));
/* The 3PN pole doesn't exist for eta=1/4. So decided to use 2PN pole always */
   ak->vpoleP6 = ak->vpoleP4;

/* For the EOBPP model, vpole and vlso are tuned to NR */
   ak->vpolePP = 0.85;
   ak->vlsoPP  = 1.0;
/*
   a = c1*c3;
   b = c1+c2+c3;
   vpole = (-b + sqrt(b*b-4.*a))/(2.*a);
   fprintf(stderr, "%e %e %e %e %e vpolePlus=%e\n", c1, c2, c3, b*b, 4.*a, vpole);
   vpole = (-b - sqrt(b*b-4.*a))/(2.*a);
   fprintf(stderr, "%e %e vpoleMinus=%e\n", ak->vlsoP4, ak->vlsoP6, vpole);
   ak->vpoleP6 = vpole;
*/

/* Taylor coefficients of flux. */
   ak->fTaN = ak->fPaN = ak->FTaN = 32.*eta*eta/5.;
   ak->FTa1 = 0.;
   ak->FTa2 = -1247./336.-35.*ieta*eta/12.;
   ak->FTa3 = 4.*LAL_PI;
   ak->FTa4 = -44711./9072.+9271.*ieta*eta/504.+65.*ieta*eta*eta/18.;
   ak->FTa5 = -(8191./672.+583./24.*ieta*eta)*LAL_PI;
   ak->FTl6 = -1712./105.;
   ak->FTa6 = 6643739519./69854400. + 16.*LAL_PI*LAL_PI/3. + ak->FTl6 * log (4.L)
            - 1712./105.*ak->EulerC+ (-134543./7776. + 41.*LAL_PI*LAL_PI/48.) * ieta*eta
            - 94403./3024. * ieta*eta*eta - 775./324. * ieta * eta*eta*eta;
   ak->FTa7 = LAL_PI * (-16285./504. + 214745./1728. * ieta*eta
               + 193385./3024.* ieta*eta*eta);
   ak->FTa8 = - 117.5043907226773;
   ak->FTl8 =   52.74308390022676;

/* Initialize members of the structures which get fed into
   phasing3() and frequency3().
*/

  ak->ptaN = -2./eta;
  ak->pta2 = 3715./8064. + 55.*ieta*eta/96.;
  ak->pta3 = -0.75*LAL_PI;
  ak->pta4 = 9.275495/14.450688 + 2.84875*ieta*eta/2.58048
           + 1855.*ieta*eta*eta/2048.;
  ak->pta5 = -(3.8645/17.2032 - 65./2048.*ieta*eta) * LAL_PI;
  ak->pta6 =  (83.1032450749357/5.7682522275840 - 53./40.*LAL_PI*LAL_PI
           - 107./56. * ak->EulerC + (-123.292747421/4.161798144
           + 2.255/2.048 *LAL_PI*LAL_PI + 385./48. * ak->lambda
           - 55./16.* ak->theta) * ieta * eta + 1.54565/18.35008 * ieta*eta*eta
           - 1.179625/1.769472 * ieta*eta*eta*eta);

  ak->pta7 =  (1.88516689/1.73408256 + 488825./516096. * ieta*eta
           - 141769./516096. * ieta*eta*eta) * LAL_PI;
  ak->ptl6 = 107./448.;

  ak->ftaN = 1./(8.*LAL_PI*ak->totalmass);
  ak->fta2 = 743./2688.+(11.*ieta*eta)/32.;
  ak->fta3 = -0.3*LAL_PI;
  ak->fta4 = 1.855099/14.450688 +  5.6975*ieta*eta/25.8048
           + 3.71*ieta*eta*eta/20.48;
  ak->fta5 = -(7.729/21.504 - 13./256.*ieta*eta) * LAL_PI;
  ak->fta6 = (-7.20817631400877/2.88412611379200 + (53./200.) * LAL_PI*LAL_PI
           + 1.07/2.80 * ak->EulerC + 1.07/2.80 * log(2.)
           + (1.23292747421/.20808990720 - 4.51/20.48 * LAL_PI*LAL_PI
           - 77./48. * ak->lambda + 11./16. * ak->theta) * ieta*eta
           - 3.0913/183.5008 * ieta*eta*eta
           + 2.35925/17.69472 * ieta*eta*eta*eta);

  ak->fta7 = (-1.88516689/4.33520640 - 97765./258048. * ieta*eta
           + 141769./1290240. * ieta*eta*eta)*LAL_PI;
  ak->ftl6 = -107./2240.;

/* Initialize members of the structures which get fed into
   phasing2() and frequency2().
*/

  ak->tvaN = -5.*ak->totalmass/(256.*eta);
  ak->tva2 = (743./252. + 11./3. * ieta * eta);
  ak->tva3 = -32./5. * LAL_PI;
  ak->tva4 = 30.58673/5.08032 + 54.29/5.04*ieta*eta + 61.7/7.2*ieta*eta*eta;
  ak->tva5 = -(77.29/2.52 -13./3.*ieta*eta) * LAL_PI;
  ak->tva6 = - 1005.2469856691/2.3471078400
           + (128./3.- 451./12.*ieta*eta) * LAL_PI*LAL_PI
           + 68.48/1.05*ak->EulerC
           + 1533.5597827/1.5240960 * ieta*eta
           - 15.211/1.728*ieta*eta*eta
           + 25.565/1.296*ieta*eta*eta*eta
           + 352./3.*ieta*eta*ak->theta
           - 2464./9.*ieta*eta*ak->lambda;
  ak->tva7 =  (-154.19335/1.27008  - 75703./756.*ieta*eta
      +  14809./378.*ieta*eta*eta) * LAL_PI;
  ak->tvl6 = 68.48/1.05;

  ak->pvaN = -1./(16.*eta);
  ak->pva2 = (3715./1008. + 55./12.*ieta*eta);
  ak->pva3 = -10.*LAL_PI;
  ak->pva4 = (15293365./1016064. + 27145./1008.*ieta*eta
             + 3085./144.*ieta*eta*eta);
  ak->pva5 = (38645./672. - 65./8. * ieta*eta) * LAL_PI;
  ak->pva6 = 1234.8611926451/1.8776862720 - 160./3. * LAL_PI*LAL_PI
           - 1712./21. * ak->EulerC
           + (-15335597827./12192768. + 2255./48. * LAL_PI*LAL_PI
           + 3080./9.*ak->lambda - 440./3. * ak->theta) * ieta*eta
           + 76055./6912. * ieta*eta*eta - 127825./5184. * ieta*eta*eta*eta;
  ak->pva7 = (77096675./2032128. + 378515./12096.*ieta*eta
           - 74045./6048.*ieta*eta*eta) * LAL_PI;
  ak->pvl6 = - 1712./21.;

/*
   Expansion coefficients for the Fourier domain phase of the usual
   stationary phase approximation
*/

  /* pfak equal 8/3* pvak - 5/3*tvak, which can be used to check that the
    * coefficients that are coded here below are exact e.g.,
   ak->pfa6 = 8./3. * ak->pva6 - 5./3. * ak->tva6;
   */

  /* For TaylorT2 and TaylorF2, we used the following references :
   * PRD72,029901(E)2005 for the 7th order which completes (corrects)
   * PRD66,027502,2002, where there are typos in Eq 1.1, 1.2, 2.9, 2.10, and
   * 2.11.. */

   ak->pfaN = 3.L/(128.L * eta);
   ak->pfa2 = 5.L*(743.L/84.L + 11.L * ieta*eta)/9.L;
   ak->pfa3 = -16.L*LAL_PI + 4.L*beta;
   ak->pfa4 = 5.L*(3058.673L/7.056L + 5429.L/7.L * ieta*eta
		   + 617.L * ieta*eta*eta)/72.L - 10.L*sigma;
   ak->pfa5 = 5.L/9.L * (7729.L/84.L - 13.L * ieta*eta) * LAL_PI;
   ak->pfl5 = 5.L/3.L * (7729.L/84.L - 13.L * ieta*eta) * LAL_PI;

   ak->pfa6 = (11583.231236531L/4.694215680L - 640.L/3.L * LAL_PI * LAL_PI - 6848.L/21.L*ak->EulerC)
     + ieta*eta * (-15335.597827L/3.048192L + 2255./12. * LAL_PI * LAL_PI - 1760./3.*ak->theta +12320./9.*ak->lambda)
    + ieta*eta*eta * 76055.L/1728.L
    - ieta*eta*eta*eta*  127825.L/1296.L ;

   ak->pfl6 = -6848.L/21.L;

   ak->pfa7 = LAL_PI * 5.L/756.L * ( 15419335.L/336.L + 75703.L/2.L * ieta*eta - 14809.L * ieta*eta*eta);

/*
  Taylor coefficients of f(v)=(1-v/vpole)F(v)
*/

   switch (params->order)
   {
      case LAL_PNORDER_NEWTONIAN:
      case LAL_PNORDER_HALF:
      case LAL_PNORDER_ONE:
      case LAL_PNORDER_ONE_POINT_FIVE:
      case LAL_PNORDER_TWO:
      case LAL_PNORDER_TWO_POINT_FIVE:
         vpole = ak->vpoleP4;
         break;
      case LAL_PNORDER_THREE:
      case LAL_PNORDER_THREE_POINT_FIVE:
         vpole = ak->vpoleP6;
         break;
      case LAL_PNORDER_PSEUDO_FOUR:
         if ( params->approximant == EOBNRv2 || params->approximant == EOBNRv2HM )
         {
           vpole = ak->vpolePP;
         }
         else
         {
           vpole = ak->vpoleP6;
         }
         break;
      default:
         XLALPrintError("XLAL Error - %s: Unknown PN order in switch\n", __func__);
         XLAL_ERROR(XLAL_EINVAL);
         break;
   }

   /* We need a different vlso for the PP model */
   if ( params->approximant == EOBNRv2  || params->approximant == EOBNRv2HM )
   {
     vlso = ak->vlsoPP;
   }
   else
   {
     vlso = ak->vlsoP4;
   }

   ak->fTa1 = ak->FTa1 - 1./vpole;
   ak->fTa2 = ak->FTa2 - ak->FTa1/vpole;
   ak->fTa3 = ak->FTa3 - ak->FTa2/vpole;
   ak->fTa4 = ak->FTa4 - ak->FTa3/vpole;
   ak->fTa5 = ak->FTa5 - ak->FTa4/vpole;
   ak->fTa6 = ak->FTa6 - ak->FTa5/vpole + ak->FTl6*log(vlso);
   ak->fTa7 = ak->FTa7 - ( ak->FTa6 + ak->FTl6*log(vlso))/vpole;
   ak->fTa8 = ak->FTa8 - ak->FTa7/vpole + ak->FTl8*log(vlso);
/*
   Pade coefficients of f(v);  assumes that a0=1 => c0=1
*/

   a1 = ak->fTa1;
   a2 = ak->fTa2;
   a3 = ak->fTa3;
   a4 = ak->fTa4;
   a5 = ak->fTa5;
   a6 = ak->fTa6;
   a7 = ak->fTa7;
   a8 = ak->fTa8;

   c1 = -a1;
   c2 = -(c1*c1 - a2)/c1;
   c3 = -(c1*pow(c2+c1,2.) + a3)/(c1*c2);
   c4 = -(c1*pow(c2+c1,3.) + c1*c2*c3*(c3+2*c2+2*c1) - a4)/(c1*c2*c3);
   c5 = -(c1*pow(pow(c1+c2,2.)+c2*c3,2.) + c1*c2*c3*pow(c1+c2+c3+c4,2.) + a5)
      /(c1*c2*c3*c4);
   c6 = -(c1*pow(c1+c2, 5.)
      + c1*c2*c3*(pow(c1+c2+c3,3.) + 3.*pow(c1+c2,3.)
      + 3.*c2*c3*(c1+c2) + c3*c3*(c2-c1))
      + c1*c2*c3*c4*(pow(c1+c2+c3+c4,2.)
      + 2.*pow(c1+c2+c3,2.) + c3*(c4-2.*c1))
      + c1*c2*c3*c4*c5*(c5 + 2.*(c1+c2+c3+c4))
      - a6)/(c1*c2*c3*c4*c5);

   a12 = a1*a1;
   a22 = a2*a2;
   a32 = a3*a3;
   a42 = a4*a4;
   a52 = a5*a5;
   a62 = a6*a6;
   a72 = a7*a7;
   a23 = a22*a2;
   a33 = a32*a3;
   a43 = a42*a4;
   a53 = a52*a5;
   a34 = a33*a3;
   a44 = a43*a4;

   c7 = ((a23 + a32 + a12*a4 - a2*(2*a1*a3 + a4)) *
      (-a44 - a32*a52 + a1*a53 - a22*a62 + a33*a7
      + a22*a5*a7 + a42*(3*a3*a5 + 2*a2*a6 + a1*a7)
      + a3*(2*a2*a5*a6 + a1*a62 - a1*a5*a7)
      - 2*a4*((a32 + a1*a5)*a6 + a2*(a52 + a3*a7))))/
      ((a33 + a1*a42 + a22*a5 - a3*(2*a2*a4 + a1*a5))*
      (a34 + a22*a42 - a43 - 2*a1*a2*a4*a5 + a12*a52
      - a2*a52 - a23*a6 - a12*a4*a6 + a2*a4*a6
      - a32*(3*a2*a4 + 2*a1*a5 + a6)
      + 2*a3*((a22 + a4)*a5 + a1*(a42 + a2*a6))));


   c8 = -(((a33 + a1*a42 + a22*a5 - a3*(2*a2*a4 + a1*a5)) *
      (pow(a4,5) + pow(a5,4) - 2*a33*a5*a6 + 2*a1*a3*a52*a6
      + 2*a3*a5*a62 + a12*a62*a6 - 2*a3*a52*a7
      + 2*a1*a32*a6*a7 - 2*a12*a5*a6*a7 + a32*a72
      + pow(a3,4)*a8 - 2*a1*a32*a5*a8 + a12*a52*a8
      - a32*a6*a8 - a43*(4*a3*a5 + 3*a2*a6 + 2*a1*a7 + a8)
      + a42*(3*a2*a52 + 3*a32*a6 + 4*a1*a5*a6 + a62
      + 4*a2*a3*a7 + 2*a5*a7 + a22*a8 + 2*a1*a3*a8)
      + a22*(a52*a6 - 2*a3*a6*a7 + 2*a3*a5*a8)
      + a22*a2*(a72 - a6*a8) - a4*(3*a52*a6 - 2*a22*a62 + 2*a33*a7
      + 4*a22*a5*a7 + a2*a72 - a2*a6*a8 - 3*a32*(a52 - a2*a8)
      + 2*a3*(a2*a5*a6 + 2*a1*a62 + a6*a7 - a5*a8)
      + 2*a1*(a52*a5 - a2*a6*a7 + a2*a5*a8) + a12*(-a72 + a6*a8))
      + a2*(-a62*a6 + 2*a5*a6*a7 + 2*a1*a5*(-a62 + a5*a7)
      + a32*(a62 + 2*a5*a7) - a52*a8
      - 2*a3*(a52*a5 + a1*(a72 - a6*a8)))))
      / ((pow(a3,4) + a22*a42 - a43 - 2*a1*a2*a4*a5
      + a12*a52 - a2*a52 - a22*a2*a6 - a12*a4*a6 + a2*a4*a6
      - a32*(3*a2*a4 + 2*a1*a5 + a6)
      + 2*a3*((a22 + a4)*a5 + a1*(a42 + a2*a6)))
      * (-pow(a4,4) - a32*a52 + a1*a52*a5 - a22*a62 + a33*a7
      + a22*a5*a7 + a42*(3*a3*a5 + 2*a2*a6 + a1*a7)
      + a3*(2*a2*a5*a6 + a1*a62 - a1*a5*a7)
      - 2*a4*((a32 + a1*a5)*a6 + a2*(a52 + a3*a7)))));

   ak->fPa1 = c1;
   ak->fPa2 = c2;
   ak->fPa3 = c3;
   ak->fPa4 = c4;
   ak->fPa5 = c5;
   ak->fPa6 = c6;
   ak->fPa7 = c7;
   ak->fPa8 = c8;

   /* spinning case */
   ak->thetahat = 1039.0/4620.0;  /* value of thetahat set according to
                         Blanchet et. al, Phys. Rev. Lett. 93, 091101 (2004) */

   ak->ST[LAL_PNORDER_NEWTONIAN] = 1.0;
   ak->ST[LAL_PNORDER_HALF] = 0.0;
   ak->ST[LAL_PNORDER_ONE] = ( -(1.0/336.0) * (743.0 + 924.0*eta) );
   ak->ST[LAL_PNORDER_ONE_POINT_FIVE] = ( 4.0 * LAL_PI );
   ak->ST[LAL_PNORDER_TWO] =  ( (34103.0 + 122949.0*eta + 59472.0*eta*eta)/18144.0 );

   ak->ST[LAL_PNORDER_TWO_POINT_FIVE] = ( -(1.0/672.0) * LAL_PI * (4159.0 + 15876.0*eta) );
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

   ak->ST[LAL_PNORDER_THREE] = ( (16447322263.0/139708800.0)
		- (1712.0/105.0)* ak->EulerC
		- (273811877.0/1088640.0)*eta - (88.0/3.0)*ak->thetahat*eta
		+ (541.0/896.0)*eta*eta - (5605.0/2592.0)*eta*eta*eta
		+ (1.0/48.0) * LAL_PI*LAL_PI * (256.0 + 451.0*eta)
		- (856.0/105.0)*log(16.0) );
   ak->ST[LAL_PNORDER_THREE+1] = ( -(1712.0/315.0) );     /* extra 3PN component */
   /* sT[8] is the LAL_PNORDER_THREE_POINT_FIVE contribution */
   ak->ST[8] = (LAL_PI/12096.0) * (-13245.0 + 717350.0*eta + 731960.0*eta*eta);
   /* coefficients 717350 and 731960 corrected (from 661775 and 599156) according
      to 2005 erratas for L. Blanchet, Phys. Rev. D 54, 1417 (1996)
      (see Phys. Rev. D 71 129904 (E) (2005)) and L. Blanchet,
      B. R. Iyer, and B. Joguet, Phys. Rev. D 65, 064005 (2002)
      (see Phys. Rev. D 71 129903 (E) (2005)).
      See errata for Arun et al., Phys. Rev. D 71, 084008
      (2005) (see  Phys. Rev. D 72 069903 (E) (2005))
      for corrected coefficients
   */

/*
   taylorcoeffs[0] = 1.;
   taylorcoeffs[1] = a1;
   taylorcoeffs[2] = a2;
   taylorcoeffs[3] = a3;
   taylorcoeffs[4] = a4;
   taylorcoeffs[5] = a5;
   taylorcoeffs[6] = a6;
   taylorcoeffs[7] = a7;
   LALPadeCoeffs7(pnorder, padecoeffs, taylorcoeffs);
   fprintf(stderr, "%e %e %e %e %e %e %e\n", a1, a2, a3, a4, a5, a6, a7);
   fprintf(stderr, "%e %e %e %e %e %e %e\n", c1, c2, c3, c4, c5, c6, c7);
   fprintf(stderr, "%e %e %e %e %e %e %e\n", padecoeffs[1], padecoeffs[2],
   padecoeffs[3], padecoeffs[4], padecoeffs[5], padecoeffs[6], padecoeffs[7]);

*/
  return XLAL_SUCCESS;

}

/* pade_coeffs.f -- translated by f2c (version 20000531).  */

/* ---------------------------------------------------------------------- */
/* {\sf pade\_coeffs (int n, real*8 as(0:n), real*8 cs(0:n))} */
/* ---------------------------------------------------------------------- */
/* 	First Version: 4.1.97. */
/* 	Purpose: Given Taylor coefficients as(0)...as(n) (n<=11) to */
/* 	        find (near) diagonal Pade coefficients cc(0)...cc(n) of the */
/*             continued fraction form at various orders. */
/* 	AUTHORS: T.Damour, B.R. Iyer and B.S.Sathyaprakash. */
/* 	Revision History: Last Update:23.6.97. */
/* 	DEPENDENCIES: None. */
/* 	INPUTS: */
/* 	   n: Array size of Taylor coefficients. */
/*           At the moment n can atmost be 11. */
/* 	   as: Array contianining the Taylor coefficients. */
/* 	OUTPUTS: */
/* 	   cs: Array contianining the Pade coefficients. */
/* 	NOTES: Damour, Iyer and Sathyaprakash (1997). */
/* ---------------------------------------------------------------------- */
/* Subroutine */
#if 0
static void LALPadeCoeffs7(int n, double *cs, double *as)
{
    /* System generated locals */
    int i__1;
    double d__1;

    /* Local variables */
    static int i__;
    static double a0, a1, c0, c1, a2, c2, a3, c3, a4, c4, a5,
	    c5, a6, c6, a7;

    static double t3, t4, t8, t9, t10, t11, t12, t13, t14, t15, t21, t18, t19,
	    t22, t23, t24, t25, t26, t27, t28, t29, t30, t32, t33, t34, t37,
	    t38, t40, t41, t42, t43, t44, t45, t46, t47, t48, t50, t51, t52,
	    t53, t54, t55, t56, t57, t58, t59, t60, t61, t62, t63, t64, t65,
	    t66, t67, t68, t69, t70, t71, t72, t73, t76, t77, t78, t80, t81,
	    t82, t83, t84, t85, t86, t87, t88, t89, t90, t91, t92, t93, t95,
	    t96, t97, t98, t99,
            t100, t101, t103, t104, t106, t107, t108,
	    t109, t110, t111, t112, t113, t114, t115, t116, t117, t118, t119,
	    t120, t121, t123, t124, t125, t126, t127, t128, t130, t131, t133,
	    t134, t135, t136, t137, t138, t139, t140, t144, t145, t148, t150,
	    t151, t153, t155, t156, t157, t158, t160, t161, t162, t163, t164,
	    t165, t167, t168, t170, t171, t173, t174, t176, t177, t179, t181,
	    t182, t185, t188, t190, t191, t192, t193, t194, t196, t199, t201,
	    t202, t203, t205, t207, t210, t212, t214, t216, t217, t219, t229,
	    t231, t233, t234, t235, t236, t237, t238, t240, t243, t244, t245,
	    t247, t248, t249, t250, t251, t252, t254, t255, t260, t261, t262,
	    t265, t267, t268, t270, t271, t273, t275, t276, t278, t280, t281,
	    t283, t284, t286, t289;

/* ---------------------------------------------------------------------- */
    a0 = as[0];
    i__1 = n;
    for (i__ = 0; i__ <= i__1; ++i__) {
	cs[i__] = 0.;
    }
    cs[0] = a0;
    if (n == 0 || cs[0] == 0.) {
	return ;
    }
    c0 = cs[0];
    a1 = as[1];
    t3 = 1 / c0;
    cs[1] = -a1 / a0;
    if (n == 1 || cs[1] == 0.) {
	return ;
    }
    c1 = cs[1];
    a2 = as[2];
    t4 = 1 / c1;
    t8 = c0 * c1;
    cs[2] = t3 * t4 * a2 - c1;
    if (n == 2 || cs[2] == 0.) {
	return ;
    }
    c2 = cs[2];
    a3 = as[3];
/* Computing 2nd power */
    d__1 = c2;
    t9 = d__1 * d__1;
    t10 = t8 * t9;
/* Computing 2nd power */
    d__1 = c1;
    t11 = d__1 * d__1;
    t12 = c0 * t11;
    t13 = t12 * c2;
    t14 = t11 * c1;
    t15 = c0 * t14;
    t18 = 1 / c2;
    t19 = t4 * t18;
    cs[3] = -(a3 + t10 + t13 * 2 + t15) * t3 * t19;
    if (n == 3 || cs[3] == 0.) {
	return ;
    }
    c3 = cs[3];
    a4 = as[4];
/* Computing 2nd power */
    d__1 = c3;
    t21 = d__1 * d__1;
    t22 = c2 * t21;
    t23 = t8 * t22;
    t24 = t9 * c3;
    t25 = t8 * t24;
    t26 = t9 * c2;
    t27 = t8 * t26;
    t28 = c2 * c3;
    t29 = t12 * t28;
    t30 = t12 * t9;
    t32 = t14 * c2 * c0;
/* Computing 2nd power */
    d__1 = t11;
    t33 = d__1 * d__1;
    t34 = c0 * t33;
    t37 = 1 / c3;
    t38 = t19 * t37;
    t40 = t11 * t26;
    t41 = t40 * c0;
/* Computing 2nd power */
    d__1 = t9;
    t42 = d__1 * d__1;
    t43 = t8 * t42;
    t44 = t11 * t9;
    t45 = c0 * c3;
    t46 = t44 * t45;
    t47 = t14 * t9;
    t48 = t47 * c0;
    t50 = t33 * c2 * c0;
    t51 = t33 * c1;
    t52 = c0 * t51;
    t53 = t15 * t28;
    cs[4] = -(-a4 + t23 + t25 * 2 + t27 + t29 * 2 + t30 * 3 + t32 * 3 + t34) *
	     t3 * t38;
    if (n == 4 || cs[4] == 0.) {
	return ;
    }
    c4 = cs[4];
    a5 = as[5];
    t54 = t28 * c4;
    t55 = t12 * t54;
    t56 = t12 * t22;
    t57 = c0 * t9;
    t58 = t21 * c1;
    t59 = t57 * t58;
    t60 = c0 * t26;
    t61 = c3 * c1;
    t62 = t60 * t61;
    t63 = t24 * c4;
    t64 = t8 * t63;
/* Computing 2nd power */
    d__1 = c4;
    t65 = d__1 * d__1;
    t66 = t28 * t65;
    t67 = t8 * t66;
    t68 = t22 * c4;
    t69 = t8 * t68;
    t70 = t21 * c3;
    t71 = c2 * t70;
    t72 = t8 * t71;
    t73 = a5 + t41 * 4 + t43 + t46 * 6 + t48 * 6 + t50 * 4 + t52 + t53 * 3 +
	    t55 * 2 + t56 * 2 + t59 * 3 + t62 * 3 + t64 * 2 + t67 + t69 * 2 +
	    t72;
    t76 = t18 * t37;
    t77 = 1 / c4;
    t78 = t76 * t77;
    t80 = t21 * t65;
    t81 = t80 * c2;
    t82 = t8 * t81;
    t83 = t26 * c3;
    t84 = t12 * t83;
    t85 = t9 * t21;
    t86 = t85 * c4;
    t87 = t8 * t86;
    t88 = c0 * t42;
    t89 = t88 * t61;
    t90 = t26 * t21;
    t91 = t8 * t90;
    t92 = c3 * c4;
    cs[5] = -t73 * t3 * t4 * t78;
    if (n == 5 || cs[5] == 0.) {
	return ;
    }
    c5 = cs[5];
    a6 = as[6];
    t93 = t92 * c5;
    t95 = t34 * t28;
    t96 = t65 * c4;
    t97 = t28 * t96;
    t98 = t8 * t97;
    t99 = t8 * c2;
/* Computing 2nd power */
    d__1 = c5;
    t100 = d__1 * d__1;
    t101 = t92 * t100;
    t103 = t21 * c4;
    t104 = t103 * c5;
    t106 = t70 * c4;
    t107 = t106 * c2;
    t108 = t8 * t107;
    t109 = t12 * t85;
    t110 = t15 * t22;
    t111 = t15 * t24;
    t112 = t12 * t66;
    t113 = t82 * -3 + a6 - t84 * 12 - t87 * 6 - t89 * 4 - t91 * 6 - t13 * 2 *
	    t93 - t95 * 4 - t98 - t99 * t101 - t99 * 2 * t104 - t108 * 3 -
	    t109 * 9 - t110 * 3 - t111 * 12 - t112 * 2;
    t114 = t24 * t65;
    t115 = t8 * t114;
    t116 = t42 * c2;
    t117 = t8 * t116;
    t118 = t15 * t54;
    t119 = t33 * t11;
    t120 = c0 * t119;
    t121 = t12 * t63;
    t123 = t51 * c2 * c0;
    t124 = t15 * t26;
/* Computing 2nd power */
    d__1 = t21;
    t125 = d__1 * d__1;
    t126 = c2 * t125;
    t127 = t8 * t126;
    t128 = t34 * t9;
    t130 = c3 * t65;
    t131 = t130 * c5;
    t133 = t12 * t42;
    t134 = t83 * c4;
    t135 = t8 * t134;
    t136 = t12 * t68;
    t137 = t9 * t70;
    t138 = t8 * t137;
    t139 = t12 * t71;
    t140 = t115 * -2 - t117 - t118 * 3 - t120 - t121 * 6 - t123 * 5 - t124 *
	    10 - t127 - t128 * 10 - t10 * 2 * t93 - t99 * 2 * t131 - t133 * 5
	    - t135 * 3 - t136 * 4 - t138 * 4 - t139 * 2;
    t144 = 1 / c5;
    t145 = t77 * t144;
    t148 = t106 * c5;
    t150 = t125 * c4;
    t151 = t150 * c2;
    t153 = t80 * c5;
    t155 = t26 * t70;
    t156 = t8 * t155;
    cs[6] = (t113 + t140) * t3 * t4 * t76 * t145;
    if (n == 6 || cs[6] == 0.) {
	return ;
    }
    c6 = cs[6];
    a7 = as[7];
    t157 = c5 * c6;
    t158 = t130 * t157;
    t160 = t125 * c3;
    t161 = c2 * t160;
    t162 = t8 * t161;
/* Computing 2nd power */
    d__1 = c6;
    t163 = d__1 * d__1;
    t164 = c5 * t163;
    t165 = t92 * t164;
    t167 = t65 * t100;
    t168 = t167 * c3;
    t170 = t96 * c5;
    t171 = t170 * c3;
    t173 = t100 * c6;
    t174 = t92 * t173;
    t176 = t100 * c5;
    t177 = t92 * t176;
    t179 = t92 * t157;
/* Computing 2nd power */
    d__1 = t65;
    t181 = d__1 * d__1;
    t182 = t28 * t181;
    t185 = t24 * t96;
    t188 = t99 * 3 * t148 + t8 * 4 * t151 + t99 * 6 * t153 + t156 * 10 + t99 *
	     2 * t158 + t162 + t99 * t165 + t99 * 3 * t168 + t99 * 3 * t171 +
	    t99 * 2 * t174 + t99 * t177 + t10 * 2 * t179 + t8 * t182 + t10 *
	    2 * t101 + t8 * 2 * t185 + t10 * 4 * t131;
    t190 = t42 * t21;
    t191 = t8 * t190;
    t192 = c0 * t116;
    t193 = t192 * t61;
    t194 = t85 * t65;
    t196 = t137 * c4;
    t199 = t83 * t65;
    t201 = t9 * t125;
    t202 = t8 * t201;
    t203 = t90 * c4;
    t205 = t22 * t96;
    t207 = t71 * t65;
    t210 = t52 * t28;
    t212 = t15 * t71;
    t214 = t10 * 6 * t104 + t191 * 10 + t193 * 5 + t8 * 9 * t194 + t8 * 12 *
	    t196 + t27 * 3 * t93 + t8 * 3 * t199 + t202 * 5 + t8 * 12 * t203
	    + t8 * 4 * t205 + t8 * 6 * t207 + t32 * 3 * t93 + t210 * 5 + t15 *
	     3 * t66 + t212 * 3 + t15 * 6 * t68;
    t216 = t42 * c3;
    t217 = t216 * c4;
    t219 = t12 * t126;
    t229 = t103 * t157;
    t231 = t103 * t100;
    t233 = t34 * t26;
    t234 = t15 * t42;
    t235 = t12 * t116;
    t236 = t8 * 4 * t217 + t219 * 2 + t30 * 6 * t93 + t12 * 6 * t81 + t12 * 6
	    * t114 + t12 * 6 * t107 + t13 * 4 * t104 + t13 * 2 * t179 + t13 *
	    2 * t101 + t13 * 4 * t131 + t12 * 2 * t97 + t99 * 2 * t229 + t99 *
	     2 * t231 + t233 * 20 + t234 * 15 + t235 * 6;
    t237 = t12 * t216;
    t238 = t15 * t85;
    t240 = t12 * t90;
    t243 = t33 * t14;
    t244 = c0 * t243;
    t245 = t34 * t24;
    t247 = t15 * t83;
    t248 = t34 * t22;
    t249 = t12 * t137;
    t250 = t42 * t9;
    t251 = t8 * t250;
    t252 = t52 * t9;
    t254 = t119 * c2 * c0;
    t255 = t237 * 20 + t238 * 18 + t12 * 18 * t86 + t240 * 24 + t15 * 12 *
	    t63 + t12 * 12 * t134 + t244 + t245 * 20 + t34 * 4 * t54 + a7 +
	    t247 * 30 + t248 * 4 + t249 * 12 + t251 + t252 * 15 + t254 * 6;
    t260 = t37 * t77;
    t261 = 1 / c6;
    t262 = t144 * t261;
    t265 = t12 * t250;
    t267 = t243 * c2 * c0;
    t268 = t71 * t96;
    t270 = t181 * c4;
    t271 = t28 * t270;
    t273 = t126 * t65;
    t275 = t160 * c4;
    t276 = t275 * c2;
    t278 = t106 * t157;
    t280 = t21 * t96;
    t281 = t280 * c5;
    t283 = t70 * t65;
    t284 = t283 * c5;
    t286 = t106 * t100;
    t289 = c4 * c5;
    cs[7] = -(t188 + t214 + t236 + t255) * t3 * t19 * t260 * t262;
    return ;
} /* LALPadeCoeffs7 */
   c7 = -(15.*pow(c1,5.)*pow(c2,2.)
      + 6.*pow(c1,2.)*c2*pow(c3,3.)*c4
      + 6.*c1*c2*pow(c3,3.)*pow(c4,2.)
      + a7
      + c1*pow(c2,6.)
      + 2.*c1*c2*pow(c3,2.)*c4*c5*c6
      + 3.*c1*c2*pow(c4,2.)*pow(c5,2.)*c3
      + c1*c2*c3*c4*c5*pow(c6,2.)
      + 6.*c1*c2*pow(c3,2.)*pow(c4,2.)*c5
      + 2.*c1*c2*c3*pow(c4,2.)*c5*c6
      + 3.*c1*c2*pow(c3,3.)*c4*c5
      + 3.*c1*c2*pow(c4,3.)*c5*c3
      + c1*c2*c3*c4*pow(c5,3.)
      + 4.*pow(c2,2.)*c1*c3*pow(c4,2.)*c5
      + 2.*pow(c2,2.)*c1*c3*c4*c5*c6
      + 5.*c1*pow(c2,2.)*pow(c3,4.)
      + 10.*c1*pow(c2,3.)*pow(c3,3.)
      + 10.*c1*pow(c2,4.)*pow(c3,2.)
      + 5.*pow(c2,5.)*c3*c1
      + pow(c1,7.)
      + 4.*pow(c1,4.)*c2*c3*c4
      + 18.*pow(c1,2.)*pow(c2,2.)*pow(c3,2.)*c4
      + 12.*pow(c1,2.)*pow(c2,3.)*c3*c4
      + 12.*pow(c1,3.)*pow(c2,2.)*c3*c4
      + 15.*pow(c1,3.)*pow(c2,4.)
      + 12.*pow(c1,2.)*pow(c2,2.)*pow(c3,3.)
      + 24.*pow(c1,2.)*pow(c2,3.)*pow(c3,2.)
      + 20.*pow(c1,2.)*pow(c2,4.)*c3
      + 18.*pow(c1,3.)*pow(c2,2.)*pow(c3,2.)
      + 30.*pow(c1,3.)*pow(c2,3.)*c3
      + 6.*pow(c1,2.)*pow(c2,5.)
      + 4.*pow(c1,4.)*c2*pow(c3,2.)
      + 20.*pow(c1,4.)*pow(c2,2.)*c3
      + c1*c2*pow(c3,5.)
      + 20.*pow(c1,4.)*pow(c2,3.)
      + 6.*c2*pow(c1,6.)
      + 5.*pow(c1,5.)*c2*c3
      + 3.*pow(c1,3.)*c2*c3*pow(c4,2.)
      + 6.*pow(c1,3.)*c2*pow(c3,2.)*c4
      + 3.*pow(c1,3.)*c2*pow(c3,3.)
      + 3.*pow(c1,3.)*c2*c3*c4*c5
      + 2.*pow(c1,2.)*c2*c3*c4*pow(c5,2.)
      + 6.*pow(c1,2.)*c2*pow(c3,2.)*pow(c4,2.)
      + 2.*pow(c1,2.)*c2*c3*pow(c4,3.)
      + 4.*pow(c1,2.)*c2*pow(c3,2.)*c4*c5
      + 4.*pow(c1,2.)*c2*c3*pow(c4,2.)*c5
      + 6.*c1*pow(c2,2.)*pow(c3,2.)*c4*c5
      + 2.*pow(c1,2.)*c2*c3*c4*c5*c6
      + 6.*pow(c1,2.)*pow(c2,2.)*c3*c4*c5
      + 2.*pow(c1,2.)*c2*pow(c3,4.)
      + 6.*pow(c1,2.)*pow(c2,2.)*c3*pow(c4,2.)
      + 2.*c1*c2*pow(c3,2.)*c4*pow(c5,2.)
      + 3.*c1*pow(c2,3.)*c3*c4*c5
      + 9.*c1*pow(c2,2.)*pow(c3,2.)*pow(c4,2.)
      + 12.*c1*pow(c2,2.)*pow(c3,3.)*c4
      + 3.*c1*pow(c2,3.)*c3*pow(c4,2.)
      + 2.*c1*c2*c3*c4*pow(c5,2.)*c6
      + 12.*c1*pow(c2,3.)*pow(c3,2.)*c4
      + c1*c2*c3*pow(c4,4.)
      + 4.*c1*c2*pow(c3,2.)*pow(c4,3.)
      + 4.*c1*pow(c3,4.)*c4*c2
      + 4.*c1*pow(c2,4.)*c3*c4
      + 2.*pow(c2,2.)*c1*c3*pow(c4,3.)
      + 2.*pow(c2,2.)*c1*c3*c4*pow(c5,2.))/(c1*c2*c3*c4*c5*c6);

   fprintf(stderr, "c7 from old formula=%e\n", c7);
#endif
