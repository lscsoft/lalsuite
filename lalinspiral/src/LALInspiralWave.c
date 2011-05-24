/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, David McKechan, B.S. Sathyaprakash, Thomas Cokelaer, Riccardo Sturani, Laszlo Vereb
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
\author Churches, D. K and Sathyaprakash, B.S. and Cokelaer. T
\file
\ingroup LALInspiral_h

\brief Interface routine needed to generate all waveforms in the \ref CBC_inspiral package.

To generate a waveform
a user is noramlly required to (a) choose the binary parameters, starting frequency, number of
bins of leading and trailing zero-padding, etc., in the structure <tt>InspiralTemplate params</tt>
and (b) call the following three functions in the order given:
LALInspiralParameterCalc(), LALInspiralWaveLength() and LALInspiralWave().
Either a time- or a frequency-domain signal is returned depending upon the
\c approximant requested (see Notes below).

\heading{Prototypes}

<tt>LALInspiralWave()</tt>
<ul>
<li> \c signalvec: Output containing the inspiral waveform.</li>
<li> \c params: Input containing binary chirp parameters.</li>
</ul>



<tt>LALInspiralWaveTemplates()</tt>
<ul>
<li> \c signalvec1: Output containing the 0-phase inspiral waveform.</li>
<li> \c signalvec2: Output containing the \f$\pi/2\f$-phase inspiral waveform.</li>
<li> \c params: Input containing binary chirp parameters.</li>
</ul>

\heading{Description}

The code LALInspiralWave() is the user interface to the inspiral codes. It takes from the user all
the physical parameters which specify the binary, and calls the relevent wave generation function.
Currently ten different approximants are fully implemented. These are #TaylorT1, #TaylorT2,
#TaylorT3, #TaylorF1, #TaylorF2, #PadeT1, #EOB, #BCV, #SpinTaylorT3, #PhenSpinTaylorRD.
\c Taylor approximants can all be generated at seven different post-Newtonian orders,
from Newtonian to 3.5 PN order, #PadeT1 exists at order 1.5PN and higher,
#EOB at orders 2 and higher. #SpinTaylorT3 is implemented only at 2PN order
by solving the evolution equations for the spin and orbital angular momenta and a
time-domain phasing formula. Finally, PN order is undefined for #BCV.
The approximant and the order are set up by the enums ::Approximant and ::LALPNOrder,
respectively.

The waveforms are all terminated one bin before the last stable orbit is reached.
The last stable orbit corresponding to a given ::Approximant and ::LALPNOrder is
defined as follows: For all \c Taylor approximants at orders 0PN, 1PN and 1.5PN
\f$v_\textrm{lso}^2=1/6,\f$ and at 2PN, 2.5PN, 3PN and 3.5PN
\f$v_\textrm{lso}^2 = x^\textrm{lso}_{T_4},\f$ where \f$x^\textrm{lso}_{T_4}\f$ is
defined in Table.\tableref{table_energy}.  In the case of \c Pade approximant
at 1.5PN order \f$v_\textrm{lso}^2=1/6,\f$ and at orders 2PN, 2.5PN, 3PN and 3.5PN
\f$v_\textrm{lso}^2 = x^\textrm{lso}_{P_4},\f$ where \f$x^\textrm{lso}_{P_4}\f$ is
defined in Table.\tableref{table_energy}. In the case of #EOB approximant,
defined only at orders greater than 2PN, the plunge waveform is
terminated at the light-ring orbit defined by Equation.\eqref{eq_LightRing}.

In the case of LALInspiralWaveTemplates() <tt>*signalvec1</tt>
contains the `0-phase' inspiral template and <tt>*signalvec2</tt> contains
a signal that is \f$\pi/2\f$ out of phase with respect to <tt>*signalvec1.</tt>
Currently, a template pair is generated only for the following \c approximants:
#TaylorT1, #TaylorT2, #TaylorT3, #PadeT1, #EOB.

See the test codes for examples of how to generate different approximations.

\heading{Algorithm}
Simple use of \c switch statement to access different PN approximations.

\heading{Uses}
Depending on the user inputs one of the following functions is called:
\code
LALInspiralWave1()
LALInspiralWave2()
LALInspiralWave3()
LALInspiralStationaryPhaseApprox1()
LALInspiralStationaryPhaseApprox2()
LALEOBWaveform()
LALBCVWaveform()
LALInspiralSpinModulatedWave()
\endcode

\heading{Notes}

<ul>
    <li> A time-domain waveform is returned when the ::Approximant is one of
    #TaylorT1, #TaylorT2, #TaylorT3, #PadeT1, #EOB, #SpinTaylorT3, #PhenSpinTaylorRD, #SpinQuadTaylor
    </li><li> A frequency-domain waveform is returned when the ::Approximant is one of
         #TaylorF1, #TaylorF2, #BCV.
         In these cases the code returns the real and imagninary parts of the
         Fourier domain signal in the convention of fftw. For a signal vector
         of length <tt>n=signalvec->length</tt> (\c n even):</li>
         <ul>
         <li> <tt>signalvec->data[0]</tt> is the \e real 0th frequency component of the Fourier transform.</li>
         <li> <tt>signalvec->data[n/2]</tt> is the \e real Nyquist frequency component of the Fourier transform.</li>
         <li> <tt>signalvec->data[k]</tt> and <tt>signalvec->data[n-k],</tt> for <tt>k=1,..., n/2-1,</tt> are
             the real and imaginary parts of the Fourier transform at a frequency \f$k\Delta f=k/T,\f$ \f$T\f$ being
             the duration of the signal and \f$\Delta f=1/T\f$ is the frequency resolution.</li>
         </ul>
</ul>

*/

#include <lal/LALInspiral.h>
#include <lal/LALNoiseModels.h>
#include <lal/LALStdlib.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/LALSQTPNWaveformInterface.h>

NRCSID (LALINSPIRALWAVEC, "$Id$");


void
LALInspiralWave(
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
{

   INITSTATUS(status, "LALInspiralWave", LALINSPIRALWAVEC);
   ATTATCHSTATUSPTR(status);

   ASSERT (signalvec,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signalvec->length >= 2, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (signalvec->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);


   ASSERT((INT4)params->approximant >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT((INT4)params->approximant < NumApproximants, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT((UINT4)params->order < (UINT4)LAL_PNORDER_NUM_ORDER,
            status, LALINSPIRALH_EORDER, LALINSPIRALH_MSGEORDER);

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
      case IMRPhenomB:
           LALBBHPhenWaveTimeDom(status->statusPtr, signalvec, params);
           CHECKSTATUSPTR(status);
       break;
      case IMRPhenomFA:
      case IMRPhenomFB:
           LALBBHPhenWaveFreqDom(status->statusPtr, signalvec, params);
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
      case PhenSpinTaylorRD:
           LALPSpinInspiralRD(status->statusPtr, signalvec, params);
           CHECKSTATUSPTR(status);
	   break;
      case PhenSpinTaylorRDF:
           LALPSpinInspiralRDFreqDom(status->statusPtr, signalvec, params);
           CHECKSTATUSPTR(status);
           break;
	  case SpinQuadTaylor:
	   		TRY(LALSQTPNWaveform(status->statusPtr, signalvec, params), status);
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
           ABORT( status, LALINSPIRALH_ESWITCH, LALINSPIRALH_MSGESWITCH );
   }

   DETATCHSTATUSPTR(status);
   RETURN (status);
}


NRCSID (LALINSPIRALWAVETEMPLATESC, "$Id$");


void
LALInspiralWaveTemplates(
   LALStatus        *status,
   REAL4Vector      *signalvec1,
   REAL4Vector      *signalvec2,
   InspiralTemplate *params
   )
{

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
   ASSERT((UINT4)params->order < (UINT4)LAL_PNORDER_NUM_ORDER,
               status, LALINSPIRALH_EORDER, LALINSPIRALH_MSGEORDER);

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
      case IMRPhenomB:
           LALBBHPhenWaveTimeDomTemplates(status->statusPtr, signalvec1, signalvec2, params);
           CHECKSTATUSPTR(status);
           break;
      case IMRPhenomFA:
      case IMRPhenomFB:
           LALBBHPhenWaveFreqDomTemplates(status->statusPtr, signalvec1, signalvec2, params);
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
      case PhenSpinTaylorRD:
	   LALPSpinInspiralRDTemplates(status->statusPtr, signalvec1, signalvec2, params);
           CHECKSTATUSPTR(status);
	   break;
      case SpinQuadTaylor:
           TRY(LALSTPNWaveformTemplates(status->statusPtr, signalvec1, signalvec2, params), status);
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
           ABORT( status, LALINSPIRALH_ESWITCH, LALINSPIRALH_MSGESWITCH );

   }

   DETATCHSTATUSPTR(status);
   RETURN (status);
}



NRCSID (LALINSPIRALWAVEFORINJECTIONC, "$Id$");


void
LALInspiralWaveForInjection(
   LALStatus        *status,
   CoherentGW       *waveform,
   InspiralTemplate *inspiralParams,
   PPNParamStruc  *ppnParams)
{

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
     case IMRPhenomB:
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
      case SpinTaylorFrameless:
           LALSTPNWaveformFramelessForInjection(status->statusPtr, waveform, inspiralParams, ppnParams);
           CHECKSTATUSPTR(status);
     break;
      case SpinTaylor:
           LALSTPNWaveformForInjection(status->statusPtr, waveform, inspiralParams, ppnParams);
           CHECKSTATUSPTR(status);
     break;
      case PhenSpinTaylorRD:
	   LALPSpinInspiralRDForInjection(status->statusPtr, waveform, inspiralParams, ppnParams);
           CHECKSTATUSPTR(status);
	   break;
	  case SpinQuadTaylor:
		   TRY(LALSQTPNWaveformForInjection(status->statusPtr, waveform, inspiralParams, ppnParams), status);
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
           ABORT( status, LALINSPIRALH_ESWITCH, LALINSPIRALH_MSGESWITCH );

   }

   DETATCHSTATUSPTR(status);
   RETURN (status);
}

