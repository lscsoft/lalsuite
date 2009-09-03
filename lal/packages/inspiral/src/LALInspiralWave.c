/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, David McKechan, B.S. Sathyaprakash, Thomas Cokelaer
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

/*  <lalVerbatim file="LALInspiralWaveCV">
Author: Churches, D. K and Sathyaprakash, B.S. and Cokelaer. T
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWave.c} and \texttt{LALInspiralWaveTemplates.c}}

Interface routine needed to generate all waveforms in the {\tt inspiral} package.

To generate a waveform
a user is noramlly required to (a) choose the binary parameters, starting frequency, number of
bins of leading and trailing zero-padding, etc., in the structure {\tt InspiralTemplate params}
and (b) call the following three functions in the order given:
{\tt LALInspiralParameterCalc,} {\tt LALInspiralWaveLength,} and {\tt LALInspiralWave.}
Either a time- or a frequency-domain signal is returned depending upon the
{\tt approximant} requested (see Notes below).

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWaveCP}
\idx{LALInspiralWave()}
\begin{itemize}
\item {\tt signalvec:} Output containing the inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\vspace{0.1in}
\input{LALInspiralWaveTemplatesCP}
\idx{LALInspiralWaveTemplates()}
\begin{itemize}
\item {\tt signalvec1:} Output containing the 0-phase inspiral waveform.
\item {\tt signalvec2:} Output containing the $\pi/2$-phase inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\subsubsection*{Description}

The code \texttt{LALInspiralWave} is the user interface to the inspiral codes. It takes from the user all
the physical parameters which specify the binary, and calls the relevent wave generation function.
Currently nine different approximants are fully implemented. These are {\tt TaylorT1, TaylorT2,
TaylorT3, TaylorF1, TaylorF2, PadeT1, EOB, BCV, SpinTaylorT3.}
{\tt Taylor} approximants can all be generated at seven different post-Newtonian orders,
from Newtonian to 3.5 PN order, {\tt PadeT1} exists at order 1.5PN and higher,
{\tt EOB} at orders 2 and higher. {\tt SpinTaylorT3} is implemented only at 2PN order
by solving the evolution equations for the spin and orbital angular momenta and a
time-domain phasing formula. Finally, PN order is undefined for {\tt BCV.}
The approximant and the order are set up by the enums \texttt{Approximant} and \texttt{Order,}
respectively.

The waveforms are all terminated one bin before the last stable orbit is reached.
The last stable orbit corresponding to a given {\tt Approximant} and {\tt Order} is
defined as follows: For all {\tt Taylor} approximants at orders 0PN, 1PN and 1.5PN
$v_{\rm lso}^2=1/6,$ and at 2PN, 2.5PN, 3PN and 3.5PN
$v_{\rm lso}^2 = x^{\rm lso}_{T_4},$ where $x^{\rm lso}_{T_4}$ is
defined in Table \ref{table:energy}.  In the case of {\tt Pade} approximant
at 1.5PN order $v_{\rm lso}^2=1/6,$ and at orders 2PN, 2.5PN, 3PN and 3.5PN
$v_{\rm lso}^2 = x^{\rm lso}_{P_4},$ where $x^{\rm lso}_{P_4}$ is
defined in Table \ref{table:energy}. In the case of {\tt EOB} approximant,
defined only at orders greater than 2PN, the plunge waveform is
terminated at the light-ring orbit defined by Equation~(\ref{eq:LightRing}).

In the case of {\tt LALInspiralWaveTemplates} {\tt *signalvec1}
contains the `0-phase' inspiral template and {\tt *signalvec2} contains
a signal that is $\pi/2$ out of phase with respect to {\tt *signalvec1.}
Currently, a template pair is generated only for the following {\tt approximants:}
{\tt TaylorT1, TaylorT2, TaylorT3, PadeT1, EOB.}

See the test codes for examples of how to generate different approximations.

\subsubsection*{Algorithm}
Simple use of {\tt switch} statement to access different PN approximations.

\subsubsection*{Uses}
Depending on the user inputs one of the following functions is called:\\
\texttt{LALInspiralWave1\\}
\texttt{LALInspiralWave2\\}
\texttt{LALInspiralWave3\\}
\texttt{LALInspiralStationaryPhaseApprox1\\}
\texttt{LALInspiralStationaryPhaseApprox2\\}
\texttt{LALEOBWaveform\\}
\texttt{LALBCVWaveform\\}
\texttt{LALInspiralSpinModulatedWave\\}

\subsubsection*{Notes}

\begin{itemize}
    \item A time-domain waveform is returned when the {\tt approximant} is one of
        {\tt TaylorT1, TaylorT2, TaylorT3, PadeT1, EOB, SpinTaylorT3}
    \item A frequency-domain waveform is returned when the {\tt approximant} is one of
        {\tt TaylorF1, TaylorF2, BCV}.
         In these cases the code returns the real and imagninary parts of the
         Fourier domain signal in the convention of fftw. For a signal vector
         of length {\tt n=signalvec->length} ({\tt n} even):
     \begin{itemize}
         \item {\tt signalvec->data[0]} is the {\it real} 0th frequency component of the Fourier transform.
         \item {\tt signalvec->data[n/2]} is the {\it real} Nyquist frequency component of the Fourier transform.
         \item {\tt signalvec->data[k]} and {\tt signalvec->data[n-k],} for {\tt k=1,\ldots, n/2-1,} are
             the real and imaginary parts of the Fourier transform at a frequency $k\Delta f=k/T,$ $T$ being
             the duration of the signal and $\Delta f=1/T$ is the frequency resolution.
     \end{itemize}
\end{itemize}
\vfill{\footnotesize\input{LALInspiralWaveCV}}

</lalLaTeX>  */

#include <lal/LALInspiral.h>
#include <lal/LALNoiseModels.h>
#include <lal/LALStdlib.h>
#include <lal/GeneratePPNInspiral.h>

NRCSID (LALINSPIRALWAVEC, "$Id$");

/*  <lalVerbatim file="LALInspiralWaveCP"> */
void
LALInspiralWave(
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
{ /* </lalVerbatim>  */

   INITSTATUS(status, "LALInspiralWave", LALINSPIRALWAVEC);
   ATTATCHSTATUSPTR(status);

   ASSERT (signalvec,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signalvec->length >= 2, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (signalvec->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);


   ASSERT((INT4)params->approximant >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT((INT4)params->approximant < NumApproximants, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT((int)params->order >= 0, status, LALINSPIRALH_EORDER, LALINSPIRALH_MSGEORDER);
   ASSERT(params->order <= 8, status, LALINSPIRALH_EORDER, LALINSPIRALH_MSGEORDER);

   switch (params->approximant)
   {
      case TaylorT1:
      case PadeT1:
           LALInspiralWave1(status->statusPtr, signalvec, params);
           CHECKSTATUSPTR(status);
	   break;
      case TaylorT2:
           LALInspiralWave2(status->statusPtr, signalvec, params);
           CHECKSTATUSPTR(status);
	   break;
      case TaylorT3:
           LALInspiralWave3(status->statusPtr, signalvec, params);
           CHECKSTATUSPTR(status);
	   break;
      case EOB:
      case EOBNR:
           LALEOBWaveform(status->statusPtr, signalvec, params);
           CHECKSTATUSPTR(status);
	   break;
      case IMRPhenomA:
           LALBBHPhenWaveTimeDom(status->statusPtr, signalvec, params);
           CHECKSTATUSPTR(status);
       break;
      case BCV:
           LALBCVWaveform(status->statusPtr, signalvec, params);
           CHECKSTATUSPTR(status);
	   break;
      case BCVSpin:
           LALBCVSpinWaveform(status->statusPtr, signalvec, params);
           CHECKSTATUSPTR(status);
	   break;
      case TaylorF1:
	   LALInspiralStationaryPhaseApprox1(status->statusPtr, signalvec, params);
           CHECKSTATUSPTR(status);
	   break;
      case TaylorF2:
      case FindChirpSP:
           LALInspiralStationaryPhaseApprox2(status->statusPtr, signalvec, params);
           CHECKSTATUSPTR(status);
	   break;
      case PadeF1:
           ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
	   break;
      case SpinTaylorT3:
           LALInspiralSpinModulatedWave(status->statusPtr, signalvec, params);
           CHECKSTATUSPTR(status);
           break;
      case SpinTaylor:
           /*GenerateTimeDomainWaveformForInjection (status->statusPtr, signalvec, params);
           CHECKSTATUSPTR(status);*/
	   /* this generate the h+ waveform whereas the one above takes into
            * account h+, hx and orientation of the system*/
           LALSTPNWaveform(status->statusPtr, signalvec, params);
           CHECKSTATUSPTR(status);
	   break;
      case AmpCorPPN:
   	   LALInspiralAmplitudeCorrectedWave(status->statusPtr, signalvec, params);
	   CHECKSTATUSPTR(status);
      	   break;
      case Eccentricity:
   	   LALInspiralEccentricity(status->statusPtr, signalvec, params);
	   CHECKSTATUSPTR(status);
      	   break;
      case TaylorEt:
   	   LALTaylorEtWaveform(status->statusPtr, signalvec, params);
	   CHECKSTATUSPTR(status);
      	   break;
      case TaylorT4:
   	   LALTaylorT4Waveform(status->statusPtr, signalvec, params);
	   CHECKSTATUSPTR(status);
      	   break;
      case TaylorN:
   	   LALTaylorNWaveform(status->statusPtr, signalvec, params);
	   CHECKSTATUSPTR(status);
      	   break;
      default:
           ABORT( status, 9999, "Unknown case in switch." );
   }

   DETATCHSTATUSPTR(status);
   RETURN (status);
}


NRCSID (LALINSPIRALWAVETEMPLATESC, "$Id$");

/*  <lalVerbatim file="LALInspiralWaveTemplatesCP"> */
void
LALInspiralWaveTemplates(
   LALStatus        *status,
   REAL4Vector      *signalvec1,
   REAL4Vector      *signalvec2,
   InspiralTemplate *params
   )
{ /* </lalVerbatim>  */

   INITSTATUS(status, "LALInspiralWaveTemplates", LALINSPIRALWAVETEMPLATESC);
   ATTATCHSTATUSPTR(status);

   ASSERT (signalvec1,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signalvec2,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signalvec1->length >= 2, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (signalvec2->length >= 2, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (signalvec1->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signalvec2->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);


   ASSERT((INT4)params->approximant >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT((INT4)params->approximant < NumApproximants, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT((INT4)params->order >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT((INT4)params->order <= 8, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   switch (params->approximant)
   {
      case TaylorT1:
      case PadeT1:
           LALInspiralWave1Templates(status->statusPtr, signalvec1, signalvec2, params);
           CHECKSTATUSPTR(status);
      	   break;
      case TaylorT2:
           LALInspiralWave2Templates(status->statusPtr, signalvec1, signalvec2, params);
           CHECKSTATUSPTR(status);
      break;
      case TaylorT3:
           LALInspiralWave3Templates(status->statusPtr, signalvec1, signalvec2, params);
           CHECKSTATUSPTR(status);
      	   break;
      case EOB:
      case EOBNR:
           LALEOBWaveformTemplates(status->statusPtr, signalvec1, signalvec2, params);
           CHECKSTATUSPTR(status);
           break;
      case IMRPhenomA:
           LALBBHPhenWaveTimeDomTemplates(status->statusPtr, signalvec1, signalvec2, params);
           CHECKSTATUSPTR(status);
           break;
      case TaylorF1:
      case TaylorF2:
      case FindChirpSP:
      case PadeF1:
      case BCV:
      case BCVSpin:
      case SpinTaylor:
           LALSTPNWaveformTemplates(status->statusPtr, signalvec1, signalvec2, params);
           CHECKSTATUSPTR(status);
           break;
      case AmpCorPPN:
      	   LALInspiralAmplitudeCorrectedWaveTemplates(status->statusPtr, signalvec1, signalvec2, params);
	   CHECKSTATUSPTR(status);
	   break;
      case Eccentricity:
           LALInspiralEccentricityTemplates(status->statusPtr, signalvec1, signalvec2, params);
           CHECKSTATUSPTR(status);
      	   break;
      default:
           ABORT( status, 9999, "Unknown case in switch." );

   }

   DETATCHSTATUSPTR(status);
   RETURN (status);
}



NRCSID (LALINSPIRALWAVEFORINJECTIONC, "$Id$");

/*  <lalVerbatim file="LALInspiralWaveForInjectionCP"> */
void
LALInspiralWaveForInjection(
   LALStatus        *status,
   CoherentGW       *waveform,
   InspiralTemplate *inspiralParams,
   PPNParamStruc  *ppnParams)
{ /* </lalVerbatim>  */

   INITSTATUS(status, "LALInspiralWaveForInjection", LALINSPIRALWAVEFORINJECTIONC);
   ATTATCHSTATUSPTR(status);

   ASSERT (inspiralParams,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (ppnParams,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);


   switch (inspiralParams->approximant)
     {
     case TaylorT1:
     case PadeT1:
       LALInspiralWave1ForInjection(status->statusPtr, waveform, inspiralParams, ppnParams);
       CHECKSTATUSPTR(status);
       break;
     case TaylorT2:
       LALInspiralWave2ForInjection(status->statusPtr, waveform, inspiralParams, ppnParams);
       CHECKSTATUSPTR(status);
       break;
     case TaylorT3:
       LALInspiralWave3ForInjection(status->statusPtr, waveform, inspiralParams, ppnParams);
       CHECKSTATUSPTR(status);
       break;
     case EOB:
     case EOBNR:
       LALEOBWaveformForInjection(status->statusPtr, waveform, inspiralParams, ppnParams);
       CHECKSTATUSPTR(status);
	   break;
     case IMRPhenomA:
       LALBBHPhenWaveTimeDomForInjection (status->statusPtr, waveform, inspiralParams, ppnParams);
       CHECKSTATUSPTR(status);
       break;
      case BCV:
	/*       LALBCVWaveformForInjection(status->statusPtr, waveform, inspiralParams, ppnParams);*/
      case BCVSpin:
      case TaylorF1:
      case TaylorF2:
      case FindChirpSP:
      case PadeF1:
           ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
	   break;
      case SpinTaylor:
           LALSTPNWaveformForInjection(status->statusPtr, waveform, inspiralParams, ppnParams);
           CHECKSTATUSPTR(status);
           break;
      case AmpCorPPN:
	   LALInspiralAmplitudeCorrectedWaveForInjection(status->statusPtr, waveform, inspiralParams, ppnParams);
	   CHECKSTATUSPTR(status);
     	   break;
      case TaylorT4:
           LALTaylorT4WaveformForInjection(status->statusPtr, waveform, inspiralParams, ppnParams);
           CHECKSTATUSPTR(status);
           break;
      default:
           ABORT( status, 9999, "Unknown case in switch." );

   }

   DETATCHSTATUSPTR(status);
   RETURN (status);
}

