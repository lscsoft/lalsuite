/**** <lalVerbatim file="LALInspiralStationaryPhaseApprox2CV">
 * Author: B.S. Sathyaprakash
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{LALInspiralStationaryPhaseApprox2.c}}
 * This module computes the usual stationary phase approximation to the
 * Fourier transform of a chirp waveform.
 * %% A one-line description of the function(s) defined in this module.
 *
 * \subsubsection*{Prototypes}
 * \input{LALStationaryPhaseApprox2CP}
 * \idx{LALStationaryPhaseApprox2()}
 *
 * \subsubsection*{Description}
 *
 * %% A description of the data analysis task performed by this function;
 * %% this is the main place to document the module.
 *
 * \subsubsection*{Algorithm}
 *
 * %% A description of the method used to perform the calculation.
 *
 * \subsubsection*{Uses}
 *
 * %% List of any external functions called by this function.
 * \begin{verbatim}
 * None
 * \end{verbatim}
 * \subsubsection*{Notes}
 *
 * %% Any relevant notes.
 *
 * \vfill{\footnotesize\input{LALInspiralStationaryPhaseApprox2CV}}
 *
 **** </lalLaTeX> */

#include "LALInspiral.h"

static void LALInspiralTaylorF2Phasing0PN (REAL8 v, REAL8 *phase, expnCoeffs *ak);
static void LALInspiralTaylorF2Phasing2PN (REAL8 v, REAL8 *phase, expnCoeffs *ak);
static void LALInspiralTaylorF2Phasing3PN (REAL8 v, REAL8 *phase, expnCoeffs *ak);
static void LALInspiralTaylorF2Phasing4PN (REAL8 v, REAL8 *phase, expnCoeffs *ak);
static void LALInspiralTaylorF2Phasing5PN (REAL8 v, REAL8 *phase, expnCoeffs *ak);
static void LALInspiralTaylorF2Phasing6PN (REAL8 v, REAL8 *phase, expnCoeffs *ak);
static void LALInspiralTaylorF2Phasing7PN (REAL8 v, REAL8 *phase, expnCoeffs *ak);

NRCSID (LALINSPIRALSTATIONARYPHASEAPPROX2C, "$Id$");

/*  <lalVerbatim file="LALInspiralStationaryPhaseApprox2CP"> */
void 
LALInspiralStationaryPhaseApprox2 
(
 LALStatus *status,
 REAL4Vector *signal,
 InspiralTemplate *params) 
 { /* </lalVerbatim>  */
   REAL8 Oneby3, h1, h2, pimmc, f, v, df, shft, phi, amp0, amp, psif, psi;
   INT4 n, nby2, i, f0, fn;
   expnCoeffs ak;
   expnFunc func;
   void (*LALInspiralTaylorF2Phasing)(REAL8 v, REAL8 *phase, expnCoeffs *ak);

   INITSTATUS (status, "LALInspiralStationaryPhaseApprox2", LALINSPIRALSTATIONARYPHASEAPPROX2C);
   ATTATCHSTATUSPTR(status);
   ASSERT (signal,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal->length>2,  status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);

   switch (params->order) 
   {
	   case newtonian:
	   case oneHalfPN:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing0PN;
		   break;
	   case onePN:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing2PN;
		   break;
	   case onePointFivePN:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing3PN;
		   break;
	   case twoPN:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing4PN;
		   break;
	   case twoPointFivePN:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing5PN;
		   break;
	   case threePN:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing6PN;
		   break;
	   case threePointFivePN:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing7PN;
		   break;
   }
   n = signal->length;
   nby2 = n/2;
   LALInspiralSetup(status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);

   Oneby3 = 1.L/3.L;
   df = params->tSampling/signal->length;
   pimmc = LAL_PI * params->totalMass * LAL_MTSUN_SI;
   f0 = params->fLower;
   fn = (params->fCutoff < ak.fn) ? params->fCutoff : ak.fn; 
   v = (pimmc*f0, Oneby3); 

   /* If we want to pad with zeroes in the beginning then the instant of
    * coalescence will be the chirp time + the duration for which padding
    * is needed. Thus, in the equation below nStartPad occurs with a +ve sign.
    * This code doesn't support non-zero start-time. i.e. params->startTime
    * should be necessarily zero.
    */
   shft = 2.L*LAL_PI * (ak.tn + params->nStartPad/params->tSampling);
   phi =  params->startPhase + LAL_PI/4.L;
   amp0 = params->signalAmplitude * ak.totalmass * pow(LAL_PI/12.L, 0.5L) * df;
/* 
   Compute the standard stationary phase approximation. 
*/
   h1 = signal->data[0] = 0.L;
   h2 = signal->data[nby2] = 0.L;
   for (i=1; i<nby2; i++) {
      f = i * df;
      if (f < f0 || f > fn)
      {
	      /* 
	       * All frequency components below f0 and above fn are set to zero  
	       */
	      signal->data[i] = 0.;
	      signal->data[n-i] = 0.;
      }
      else
      {
	      v = pow(pimmc * f, Oneby3);
	      LALInspiralTaylorF2Phasing(v, &psif, &ak);
	      psi = shft * f + phi + psif;
	      /* 
		 dEnergybyFlux computes 1/(dv/dt) while what we need is 1/(dF/dt):
		 dF/dt=(dF/dv)(dv/dt)=-dEnergybyFlux/(dF/dv)=-dEnergybyFlux 3v^2/(pi m^2)
		 Note that our energy is defined as E=-eta v^2/2 and NOT -eta m v^2/2.
		 This is the reason why there is an extra m in the last equation above
		 amp = amp0 * pow(-dEnergybyFlux(v)/v^2, 0.5) * v^2;
		     = amp0 * pow(-dEnergybyFlux(v), 0.5) * v;
		     
	       */
	      amp = amp0 * pow(-func.dEnergy(v,&ak)/func.flux(v,&ak),0.5L) * v;
	      signal->data[i] = (REAL4) (amp * cos(psi));
	      signal->data[n-i] = (REAL4) (-amp * sin(psi));

      }

      /*
	 printf ("%e %e \n", v, psif);
	 printf ("%e %e %e %e %e\n", f, pow(h1,2.)+pow(h2,2.), h2, psi, psif);
	 printf ("&\n");
       */
   
   }
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

/*  <lalVerbatim file="LALInspiralTaylorF2Phasing0PNCP"> */
static void LALInspiralTaylorF2Phasing0PN (REAL8 v, REAL8 *phase, expnCoeffs *ak) {
/* </lalVerbatim>  */

   REAL8 x;
   x = v*v;
   *phase = ak->pfaN/pow(v,5.);
}

/*  <lalVerbatim file="LALInspiralTaylorF2Phasing2PNCP"> */
static void LALInspiralTaylorF2Phasing2PN (REAL8 v, REAL8 *phase, expnCoeffs *ak) {
/* </lalVerbatim>  */
   REAL8 x;
   x = v*v;
   *phase = ak->pfaN * (1. + ak->pfa2 * x)/pow(v,5.);
}

/*  <lalVerbatim file="LALInspiralTaylorF2Phasing3PNCP"> */
static void LALInspiralTaylorF2Phasing3PN (REAL8 v, REAL8 *phase, expnCoeffs *ak) {
/* </lalVerbatim>  */
   REAL8 x;
   x = v*v;
   *phase = ak->pfaN * (1. + ak->pfa2 * x + ak->pfa3 * v*x)/pow(v,5.);
}

/*  <lalVerbatim file="LALInspiralTaylorF2Phasing4PNCP"> */
static void LALInspiralTaylorF2Phasing4PN (REAL8 v, REAL8 *phase, expnCoeffs *ak) {
/* </lalVerbatim>  */
   REAL8 x;
   x = v*v;
   *phase = ak->pfaN * (1. + ak->pfa2 * x + ak->pfa3 * v*x + ak->pfa4 * x*x)/pow(v,5.);
}

/*  <lalVerbatim file="LALInspiralTaylorF2Phasing5PNCP"> */
static void LALInspiralTaylorF2Phasing5PN (REAL8 v, REAL8 *phase, expnCoeffs *ak) {
/* </lalVerbatim>  */
   REAL8 x, y;
   x = v*v;
   y = x*x;
   *phase = ak->pfaN * (1. + ak->pfa2 * x + ak->pfa3 * v*x + ak->pfa4 * y
         + (ak->pfa5 + ak->pfl5 * log(v/ak->v0)) * y*v)/pow(v,5.);
}

/*  <lalVerbatim file="LALInspiralTaylorF2Phasing6PNCP"> */
static void LALInspiralTaylorF2Phasing6PN (REAL8 v, REAL8 *phase, expnCoeffs *ak) {
/* </lalVerbatim>  */
   REAL8 x, y, z;
   x = v*v;
   y = x*x;
   z = y*x;
   /*
   fprintf(stderr, "This order is currently not implemented\n");
   exit(0);
   */
   *phase = ak->pfaN * (1. + ak->pfa2 * x + ak->pfa3 * v*x + ak->pfa4 * y
         + (ak->pfa5 + ak->pfl5 * log(v/ak->v0)) * y*v + ak->pfa6 * z)/pow(v,5.);
}

/*  <lalVerbatim file="LALInspiralTaylorF2Phasing7PNCP"> */
static void LALInspiralTaylorF2Phasing7PN (REAL8 v, REAL8 *phase, expnCoeffs *ak) {
/* </lalVerbatim>  */
   REAL8 x, y, z;
   x = v*v;
   y = x*x;
   /*
   fprintf(stderr, "This order is currently not implemented\n");
   exit(0);
   */
   *phase = ak->pfaN * (1. + ak->pfa2 * x + ak->pfa3 * v*x + ak->pfa4 * y
         + (ak->pfa5 + ak->pfl5 * log(v/ak->v0)) * y*v + ak->pfa6 * z + ak->pfa7*z*v)/pow(v,5.);
}
