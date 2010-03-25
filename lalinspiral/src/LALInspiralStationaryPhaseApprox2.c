/*
*  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
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

/**** <lalVerbatim file="LALInspiralStationaryPhaseApprox2CV">
 * Author: B.S. Sathyaprakash
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 *
 * \subsection{Module \texttt{LALInspiralStationaryPhaseApprox2.c}}
 * %% A one-line description of the function(s) defined in this module.
 * This module computes the usual stationary phase approximation to the
 * Fourier transform of a chirp waveform
 * given by Eq.~(\ref{eq:InspiralFourierPhase:f2}).
 *
 * \subsubsection*{Prototypes}
 * \input{LALInspiralStationaryPhaseApprox2CP}
 * \idx{LALInspiralStationaryPhaseApprox2()}
 * \begin{itemize}
 * \item {\tt signalvec:} Output containing the inspiral waveform.
 * \item {\tt params:} Input containing binary chirp parameters.
 * \end{itemize}
 *
 * \subsubsection*{Description}
 *
 * %% A description of the data analysis task performed by this function;
 * %% this is the main place to document the module.
 * Computes the Fourier transform of the chirp signal in the stationary
 * phase approximation and returns the real and imagninary parts of the
 * Fourier domain signal in the convention of fftw. For a signal vector
 * of length {\tt n=signalvec->length} ({\tt n} even):
 * \begin{itemize}
 * \item {\tt signalvec->data[0]} is the {\it real} 0th frequency component of the Fourier transform.
 * \item {\tt signalvec->data[n/2]} is the {\it real} Nyquist frequency component of the Fourier transform.
 * \item {\tt signalvec->data[k]} and {\tt signalvec->data[n-k],} for {\tt k=1,\ldots, n/2-1,} are
 * the real and imaginary parts of the Fourier transform at a frequency $k\Delta f=k/T,$ $T$ being
 * the duration of the signal and $\Delta f=1/T$ is the frequency resolution.
 * \end{itemize}
 *
 * \subsubsection*{Algorithm}
 *
 * %% A description of the method used to perform the calculation.
 *
 * The standard SPA is given by Eq.~(\ref{eq:InspiralFourierPhase:f2}).
 * We define a variable function pointer {\tt LALInspiralTaylorF2Phasing} and point
 * it to one of the {\texttt static} functions defined within this function
 * that explicitly calculates the Fourier phase at the PN order chosen by the user.
 * The reference points are chosen so that on inverse Fourier transforming
 * the time-domain waveform will
 * \begin{itemize}
 * \item be padded with zeroes in the first {\tt params->nStartPad} bins,
 * \item begin with a phase shift of {\tt params->nStartPhase} radians,
 * \item have an amplitude of ${\tt n} v^2.$
 * \end{itemize}
 * \subsubsection*{Uses}
 * \begin{verbatim}
   LALInspiralSetup
   LALInspiralChooseModel
   LALInspiralTaylorF2Phasing[0234567]PN
 * \end{verbatim}
 *
 * %% List of any external functions called by this function.
 * \begin{verbatim}
 * None
 * \end{verbatim}
 * \subsubsection*{Notes}
 *
 * %% Any relevant notes.
 *
 * If it is required to compare the output of this module with a time domain
 * signal one should use an inverse Fourier transform routine that packs data
 * in the same way as fftw. Moreover, one should divide the resulting inverse
 * Fourier transform by a factor ${\tt n}/2$ to be consistent with the
 * amplitude used in time-domain signal models.
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
LALInspiralStationaryPhaseApprox2 (
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
{ /* </lalVerbatim>  */
   REAL8 Oneby3, h1, h2, pimmc, f, v, df, shft, phi, amp0, amp, psif, psi;
   INT4 n, nby2, i, f0, fn;
   expnCoeffs ak;
   expnFunc func;
   void (*LALInspiralTaylorF2Phasing)(REAL8, REAL8 *, expnCoeffs *) = NULL;

   INITSTATUS (status, "LALInspiralStationaryPhaseApprox2", LALINSPIRALSTATIONARYPHASEAPPROX2C);
   ATTATCHSTATUSPTR(status);
   ASSERT (signalvec,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signalvec->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signalvec->length>2,  status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);

   switch (params->order)
   {
	   case LAL_PNORDER_NEWTONIAN:
	   case LAL_PNORDER_HALF:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing0PN;
		   break;
	   case LAL_PNORDER_ONE:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing2PN;
		   break;
	   case LAL_PNORDER_ONE_POINT_FIVE:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing3PN;
		   break;
	   case LAL_PNORDER_TWO:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing4PN;
		   break;
	   case LAL_PNORDER_TWO_POINT_FIVE:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing5PN;
		   break;
	   case LAL_PNORDER_THREE:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing6PN;
		   break;
	   case LAL_PNORDER_THREE_POINT_FIVE:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing7PN;
		   break;
           default:
                   ABORT( status, LALINSPIRALH_EORDERMISSING, LALINSPIRALH_MSGEORDERMISSING );
   }
   n = signalvec->length;
   nby2 = n/2;
   memset( &ak, 0, sizeof( ak ) );
   LALInspiralSetup(status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);

   Oneby3 = 1.L/3.L;
   df = params->tSampling/signalvec->length;
   pimmc = LAL_PI * params->totalMass * LAL_MTSUN_SI;
   f0 = params->fLower;
   fn = (params->fCutoff < ak.fn) ? params->fCutoff : ak.fn;
   v = pow(pimmc*f0, Oneby3);

   /* If we want to pad with zeroes in the beginning then the instant of
    * coalescence will be the chirp time + the duration for which padding
    * is needed. Thus, in the equation below nStartPad occurs with a +ve sign.
    * This code doesn't support non-zero start-time. i.e. params->startTime
    * should be necessarily zero.
    */
   shft = 2.L*LAL_PI * (ak.tn + params->nStartPad/params->tSampling + params->startTime);
   phi =  params->startPhase + LAL_PI/4.L;
   amp0 = params->signalAmplitude * ak.totalmass * pow(LAL_PI/12.L, 0.5L) * df;
/*
   Compute the standard stationary phase approximation.
*/
   h1 = signalvec->data[0] = 0.L;
   h2 = signalvec->data[nby2] = 0.L;
   for (i=1; i<nby2; i++) {
      f = i * df;
      if (f < f0 || f > fn)
      {
	      /*
	       * All frequency components below f0 and above fn are set to zero
	       */
	      signalvec->data[i] = 0.;
	      signalvec->data[n-i] = 0.;
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
	      signalvec->data[i] = (REAL4) (amp * cos(psi));
	      signalvec->data[n-i] = (REAL4) (-amp * sin(psi));

      }
      /*
	 printf ("%e %e \n", v, psif);
	 printf ("%e %e %e %e %e\n", f, pow(h1,2.)+pow(h2,2.), h2, psi, psif);
	 printf ("&\n");
       */

   }
   params->fFinal = fn;
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

   *phase = ak->pfaN * (1. + ak->pfa2 * x + ak->pfa3 * v*x + ak->pfa4 * y
         + (ak->pfa5 + ak->pfl5 * log(v/ak->v0)) * y*v
         + (ak->pfa6 + ak->pfl6 * log(4.*v) ) * z)
     /pow(v,5.);
}

/*  <lalVerbatim file="LALInspiralTaylorF2Phasing7PNCP"> */
static void LALInspiralTaylorF2Phasing7PN (REAL8 v, REAL8 *phase, expnCoeffs *ak) {
/* </lalVerbatim>  */
   REAL8 x, y, z;
   x = v*v;
   y = x*x;
   z = y*x;

   *phase = ak->pfaN * (1. + ak->pfa2 * x + ak->pfa3 * v*x + ak->pfa4 * y
         + (ak->pfa5 + ak->pfl5 * log(v/ak->v0)) * y*v
         + (ak->pfa6 + ak->pfl6 * log(4.*v) ) * z
         + ak->pfa7 * z*v)
     /pow(v,5.);
}
