/* 
   setup code for giving all the taylor and pade coefficients
   call this first before any template generation
   Created: February 17, 2000.
   Author: B.S.Sathyaprakash, Cardiff University
   Revision History: 
	Feb 22: Added total mass and ieta in expnCoeffs.
      March 19: inspiralsetup should now be called with InspiralTemplate
                as an arguemnt. All params are setup here except the model,
                which still will have to be set using chooseModel.
   Purpose: To set various constants for use by functions 
      generating Pade and Taylor approximants of flux and energy.
   Dependecies: 
   Inputs: inspiral signal parameter array.
   Outputs: 
     Taylor coefficients of E(v).
     Taylor coefficients of dE(v)/dv.
     Taylor coefficients of e(x).
     Pade coefficients of e(x).
     Taylor coefficients of F(v).
     Taylor coefficients of f(v).
     Pade coefficients of f(v).
*/
   
#include <math.h>
#include "Inspiral.h"
#include "LALStdlib.h"

NRCSID (INSPIRALSETUPC, "$Id$");


void LALInspiralSetup (LALStatus *status,
		    expnCoeffs *ak,
		    InspiralTemplate *params) {
   
   int ieta;
   REAL8 lso, eta, vpole;
   REAL8 a1, a2, a3, a4, a5;
   REAL8 c1, c2, c3, c4, c5;
   REAL8 oneby6=1.0/6.0;

   INITSTATUS (status, "LALInspiralSetup", INSPIRALSETUPC);

   ASSERT (ak,  status, INSPIRALSETUP_ENULL, INSPIRALSETUP_MSGENULL);
   ASSERT (params,  status, INSPIRALSETUP_ENULL, INSPIRALSETUP_MSGENULL);


   ak->ieta = params->ieta;
   ak->t0 = params->startTime;
   ak->m1 = params->mass1;
   ak->m2 = params->mass2;
   ak->f0 = params->fLower;
   ak->fn = params->fCutoff;
   ak->samplingrate = params->tSampling;
   ak->samplinginterval = 1./ak->samplingrate;

   lso = pow(oneby6, 0.5);
/* 
   ieta determines the nature of the waveforms: 
   ieta=0 testmass waveforms 
   ieta=1 comparable mass waveforms.
*/
   ieta = ak->ieta;

/* Compute the total mass and eta from m1 and m2 */
   ak->totalmass = (ak->m1 + ak->m2) * LAL_MTSUN_SI;
   ak->eta = ak->m1 * ak->m2 / pow(ak->m1 + ak->m2, 2.);
   eta = ak->eta;

/* Set initial velocity according to initial frequency */

   ak->v0 = pow (LAL_PI * ak->totalmass * ak->f0, 1./3.);

/* Taylor coefficients of E(x) */
   ak->ETaN = -eta/2.;
   ak->ETa1 = -(9.+ieta*eta)/12.;
   ak->ETa2 = -(27.-19*ieta*eta+ieta*eta*eta/3.)/8.;

/* Taylor coefficients of e(x) */
   ak->eTaN = -1.;
   ak->eTa1 = -1.-ieta*eta/3.;
   ak->eTa2 = -3.+35.*ieta*eta/12.;

/* Taylor coefficients of dE(v)/dv. (NOTE v and NOT x) */
   ak->dETaN = -eta;
   ak->dETa1 = -(9. + ieta*eta)/6.;
   ak->dETa2 = -(3./8.) * (27. - 19.*ieta*eta + ieta*eta*eta/3.);

/* Pade coefficients of e(x)-function. */
   ak->ePaN = -1.;
   ak->ePa1 = 1.+ieta*eta/3.;
   ak->ePa2 = -(144. - 81.*ieta*eta + 4.*ieta*eta*eta) / (36.+12*ieta*eta);

/* Location of the 2PN T- and P-approximant last stable orbit and pole: */
   ak->vlsoT0 = lso;
   ak->vlsoP0 = lso;
   ak->vlsoP2 = lso;
/* 
   vlsoT2 =  6./(9.+ieta*eta); 
   This correct value makes vlso too large for vlsoT2 hence use 1/sqrt(6)
*/
   ak->vlsoT2 = pow(1./6., 0.5);   
   ak->vlsoT4 = pow(-ak->ETa1 + pow(ak->ETa1*ak->ETa1 - 3*ak->ETa2,0.5)/(3*ak->ETa2), 0.5);
   ak->vlsoP4 = pow((-1.+pow(-ak->ePa1/ak->ePa2,0.5))/(ak->ePa1 + ak->ePa2), 0.5);
   ak->vpoleP4 = vpole = pow(4.*(3.+ieta*eta)/(36.-35.*ieta*eta), 0.5);

/* Taylor coefficients of flux. */
   ak->fTaN = ak->fPaN = ak->FTaN = 32.*eta*eta/5.;
   ak->FTa1 = 0.;
   ak->FTa2 = -1247./336.-35.*ieta*eta/12.;
   ak->FTa3 = 4.*LAL_PI;
   ak->FTa4 = -44711./9072.+9271.*ieta*eta/504.+65.*ieta*eta*eta/18.;
   ak->FTa5 = -(8191./672.+535./24.*ieta*eta)*LAL_PI;

/* Taylor coefficients of f(v)=(1-v/vpole)F(v) */

   ak->fTa1 = ak->FTa1 - 1./vpole;
   ak->fTa2 = ak->FTa2 - ak->FTa1/vpole;
   ak->fTa3 = ak->FTa3 - ak->FTa2/vpole;
   ak->fTa4 = ak->FTa4 - ak->FTa3/vpole;
   ak->fTa5 = ak->FTa5 - ak->FTa4/vpole;

/* Pade coefficients of f(v);  assumes that a0=1 => c0=1 */

   a1 = ak->fTa1;
   a2 = ak->fTa2;
   a3 = ak->fTa3;
   a4 = ak->fTa4;
   a5 = ak->fTa5;

   c1 = -a1;
   c2 = -(c1*c1 - a2)/c1;
   c3 = -(c1*pow(c2+c1,2.) + a3)/(c1*c2);
   c4 = -(c1*pow(c2+c1,3.) + c1*c2*c3*(c3+2*c2+2*c1) - a4)/(c1*c2*c3);
   c5 = -(c1*pow(pow(c1+c2,2.)+c2*c3,2.) + c1*c2*c3*pow(c1+c2+c3+c4,2.) + a5)/(c1*c2*c3*c4);

   ak->fPa1 = c1;
   ak->fPa2 = c2;
   ak->fPa3 = c3;
   ak->fPa4 = c4;
   ak->fPa5 = c5;

   RETURN (status);

}
