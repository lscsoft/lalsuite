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
 * \author Churches, D. K and Sathyaprakash, B.S. and Cokelaer. T
 * \file
 * \ingroup LALInspiral_h
 *
 * \brief Interface routine needed to generate all waveforms in the \ref pkg_inspiral package.
 *
 * To generate a waveform
 * a user is noramlly required to (a) choose the binary parameters, starting frequency, number of
 * bins of leading and trailing zero-padding, etc., in the structure <tt>InspiralTemplate params</tt>
 * and (b) call the following three functions in the order given:
 * LALInspiralParameterCalc(), LALInspiralWaveLength() and LALInspiralWave().
 * Either a time- or a frequency-domain signal is returned depending upon the
 * \c approximant requested (see Notes below).
 *
 * \heading{Prototypes}
 *
 * <tt>LALInspiralWave()</tt>
 * <ul>
 * <li> \c signalvec: Output containing the inspiral waveform.</li>
 * <li> \c params: Input containing binary chirp parameters.</li>
 * </ul>
 *
 * <tt>LALInspiralWaveTemplates()</tt>
 * <ul>
 * <li> \c signalvec1: Output containing the 0-phase inspiral waveform.</li>
 * <li> \c signalvec2: Output containing the \f$\pi/2\f$-phase inspiral waveform.</li>
 * <li> \c params: Input containing binary chirp parameters.</li>
 * </ul>
 *
 * \heading{Description}
 *
 * The code LALInspiralWave() is the user interface to the inspiral codes. It takes from the user all
 * the physical parameters which specify the binary, and calls the relevent wave generation function.
 * Currently ten different approximants are fully implemented. These are #TaylorT1, #TaylorT2,
 * #TaylorT3, #TaylorF1, #TaylorF2, #TaylorF2RedSpin, #PadeT1, #EOB, #BCV, #SpinTaylorT3, #PhenSpinTaylorRD.
 * \c Taylor approximants can all be generated at seven different post-Newtonian orders,
 * from Newtonian to 3.5 PN order, #PadeT1 exists at order 1.5PN and higher,
 * #EOB at orders 2 and higher. #SpinTaylorT3 is implemented only at 2PN order
 * by solving the evolution equations for the spin and orbital angular momenta and a
 * time-domain phasing formula. Finally, PN order is undefined for #BCV.
 * The approximant and the order are set up by the enums ::Approximant and ::LALPNOrder,
 * respectively.
 *
 * The waveforms are all terminated one bin before the last stable orbit is reached.
 * The last stable orbit corresponding to a given ::Approximant and ::LALPNOrder is
 * defined as follows: For all \c Taylor approximants at orders 0PN, 1PN and 1.5PN
 * \f$v_\textrm{lso}^2=1/6,\f$ and at 2PN, 2.5PN, 3PN and 3.5PN
 * \f$v_\textrm{lso}^2 = x^\textrm{lso}_{T_4},\f$ where \f$x^\textrm{lso}_{T_4}\f$ is
 * defined in Table.\tableref{table_energy}.  In the case of \c Pade approximant
 * at 1.5PN order \f$v_\textrm{lso}^2=1/6,\f$ and at orders 2PN, 2.5PN, 3PN and 3.5PN
 * \f$v_\textrm{lso}^2 = x^\textrm{lso}_{P_4},\f$ where \f$x^\textrm{lso}_{P_4}\f$ is
 * defined in Table.\tableref{table_energy}. In the case of #EOB approximant,
 * defined only at orders greater than 2PN, the plunge waveform is
 * terminated at the light-ring orbit defined by Equation.\eqref{eq_LightRing}.
 *
 * In the case of LALInspiralWaveTemplates() <tt>*signalvec1</tt>
 * contains the `0-phase' inspiral template and <tt>*signalvec2</tt> contains
 * a signal that is \f$\pi/2\f$ out of phase with respect to <tt>*signalvec1.</tt>
 * Currently, a template pair is generated only for the following \c approximants:
 * #TaylorT1, #TaylorT2, #TaylorT3, #PadeT1, #EOB.
 *
 * See the test codes for examples of how to generate different approximations.
 *
 * \heading{Algorithm}
 * Simple use of \c switch statement to access different PN approximations.
 *
 * \heading{Uses}
 * Depending on the user inputs one of the following functions is called:
 * \code
 * LALInspiralWave1()
 * LALInspiralWave2()
 * LALInspiralWave3()
 * XLALInspiralStationaryPhaseApprox1()
 * XLALInspiralStationaryPhaseApprox2()
 * LALEOBWaveform()
 * LALBCVWaveform()
 * LALInspiralSpinModulatedWave()
 * \endcode
 *
 * \heading{Notes}
 *
 * <ul>
 * <li> A time-domain waveform is returned when the ::Approximant is one of
 * #TaylorT1, #TaylorT2, #TaylorT3, #PadeT1, #EOB, #SpinTaylorT3, #PhenSpinTaylorRD, #SpinQuadTaylor
 * </li><li> A frequency-domain waveform is returned when the ::Approximant is one of
 * #TaylorF1, #TaylorF2, #TaylorF2RedSpin, #BCV.
 * In these cases the code returns the real and imagninary parts of the
 * Fourier domain signal in the convention of fftw. For a signal vector
 * of length <tt>n=signalvec->length</tt> (\c n even):</li>
 * <ul>
 * <li> <tt>signalvec->data[0]</tt> is the \e real 0th frequency component of the Fourier transform.</li>
 * <li> <tt>signalvec->data[n/2]</tt> is the \e real Nyquist frequency component of the Fourier transform.</li>
 * <li> <tt>signalvec->data[k]</tt> and <tt>signalvec->data[n-k],</tt> for <tt>k=1,..., n/2-1,</tt> are
 * the real and imaginary parts of the Fourier transform at a frequency \f$k\Delta f=k/T,\f$ \f$T\f$ being
 * the duration of the signal and \f$\Delta f=1/T\f$ is the frequency resolution.</li>
 * </ul>
 * </ul>
 *
 */

#include <lal/LALInspiral.h>
#include <lal/LALNoiseModels.h>
#include <lal/LALStdlib.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/LALSQTPNWaveformInterface.h>
#include <lal/TimeSeries.h>

/**
 * Generate the plus and cross polarizations for a waveform
 * form a row of the sim_inspiral table.
 *
 * Parses a row from the sim_inspiral table and passes the appropriate members
 * to XLALSimInspiralChooseWaveform().
 *
 * FIXME: this should eventually be moved to lalsimulation
 * along with the appropriate string parsing functions
 */
int XLALSimInspiralChooseWaveformFromSimInspiral(
    REAL8TimeSeries **hplus,	/**< +-polarization waveform */
    REAL8TimeSeries **hcross,	/**< x-polarization waveform */
    SimInspiralTable *thisRow,	/**< row from the sim_inspiral table containing waveform parameters */
    REAL8 deltaT		/**< time step */
    )
{
   int ret;
   LALPNOrder order;
   Approximant approximant;
   LALSimInspiralApplyTaper taper;

   REAL8 phi0 = thisRow->coa_phase;
   REAL8 m1 = thisRow->mass1 * LAL_MSUN_SI;
   REAL8 m2 = thisRow->mass2 * LAL_MSUN_SI;
   REAL8 S1x = thisRow->spin1x;
   REAL8 S1y = thisRow->spin1y;
   REAL8 S1z = thisRow->spin1z;
   REAL8 S2x = thisRow->spin2x;
   REAL8 S2y = thisRow->spin2y;
   REAL8 S2z = thisRow->spin2z;
   REAL8 f_min = thisRow->f_lower;
   REAL8 f_ref = 0.;
   REAL8 r = thisRow->distance * LAL_PC_SI * 1e6;
   REAL8 i = thisRow->inclination;
   REAL8 lambda1 = 0.; /* FIXME:0 turns these terms off, these should be obtained by some other means */
   REAL8 lambda2 = 0.; /* FIXME:0 turns these terms off, these should be obtained by some other means */
   LALSimInspiralWaveformFlags *waveFlags=XLALSimInspiralCreateWaveformFlags();
   LALSimInspiralTestGRParam *nonGRparams = NULL;
   int amplitudeO = thisRow->amp_order;

   /* get approximant */
   approximant = XLALGetApproximantFromString(thisRow->waveform);
   if ( (int) approximant == XLAL_FAILURE)
      XLAL_ERROR(XLAL_EFUNC);

   /* get phase PN order; this is an enum such that the value is twice the PN order */
   order = XLALGetOrderFromString(thisRow->waveform);
   if ( (int) order == XLAL_FAILURE)
      XLAL_ERROR(XLAL_EFUNC);

   /* get taper option */
   taper = XLALGetTaperFromString(thisRow->taper);
   if ( (int) taper == XLAL_FAILURE)
      XLAL_ERROR(XLAL_EFUNC);

   /* generate +,x waveforms */
   /* special case for NR waveforms */
   switch(approximant)
   {
      case NumRelNinja2:
         if (XLALNRInjectionFromSimInspiral(hplus, hcross, thisRow, deltaT) == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
         break;

      default:
         ret = XLALSimInspiralChooseTDWaveform(hplus, hcross, phi0, deltaT,
               m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, i,
               lambda1, lambda2, waveFlags, nonGRparams, amplitudeO,
               order, approximant);
         XLALSimInspiralDestroyWaveformFlags(waveFlags);
         XLALSimInspiralDestroyTestGRParam(nonGRparams);
         if( ret == XLAL_FAILURE )
            XLAL_ERROR(XLAL_EFUNC);
   }

   /* taper the waveforms */
   if (XLALSimInspiralREAL8WaveTaper((*hplus)->data, taper) == XLAL_FAILURE)
      XLAL_ERROR(XLAL_EFUNC);

   if (XLALSimInspiralREAL8WaveTaper((*hcross)->data, taper) == XLAL_FAILURE)
      XLAL_ERROR(XLAL_EFUNC);

   return XLAL_SUCCESS;
}

/**
 * Generate the plus and cross polarizations for a waveform
 * form a row of the InspiralTemplate structure.
 *
 * Parses the InspiralTemplate stucture and passes the appropriate members
 * to XLALSimInspiralChooseWaveform().
 */
int
XLALSimInspiralChooseWaveformFromInspiralTemplate(
   REAL8TimeSeries **hplus,	/**< +-polarization waveform */
   REAL8TimeSeries **hcross,	/**< x-polarization waveform */
   InspiralTemplate *params	/**< stucture containing waveform parameters */
   )
{
  int ret;
  REAL8 deltaT = 1./params->tSampling;
  REAL8 phi0 = params->startPhase; /* startPhase is used as the peak amplitude phase here */
  REAL8 m1 = params->mass1 * LAL_MSUN_SI;
  REAL8 m2 = params->mass2 * LAL_MSUN_SI;
  REAL8 S1x = params->spin1[0];
  REAL8 S1y = params->spin1[1];
  REAL8 S1z = params->spin1[2];
  REAL8 S2x = params->spin2[0];
  REAL8 S2y = params->spin2[1];
  REAL8 S2z = params->spin2[2];
  REAL8 f_min = params->fLower;
  REAL8 f_ref = 0.;
  /* Value of 'distance' fed to lalsim is conventional to obtain a correct template norm */
  REAL8 r = params->distance;
  REAL8 i = params->inclination;
  REAL8 lambda1 = 0.; /* FIXME:0 turns these terms off, these should be obtained by some other means */
  REAL8 lambda2 = 0.; /* FIXME:0 turns these terms off, these should be obtained by some other means */
  LALSimInspiralWaveformFlags *waveFlags = XLALSimInspiralCreateWaveformFlags();
  LALSimInspiralTestGRParam *nonGRparams = NULL;
  LALPNOrder amplitudeO = params->ampOrder;
  LALPNOrder order = params->order;
  Approximant approximant = params->approximant;

  /* generate +,x waveforms */
  ret = XLALSimInspiralChooseTDWaveform(hplus, hcross, phi0, deltaT, m1, m2,
            S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, i, lambda1, lambda2,
            waveFlags, nonGRparams, amplitudeO, order, approximant);
  XLALSimInspiralDestroyWaveformFlags(waveFlags);
  XLALSimInspiralDestroyTestGRParam(nonGRparams);
  if( ret == XLAL_FAILURE)
    XLAL_ERROR(XLAL_EFUNC);

  return XLAL_SUCCESS;
}


void
LALInspiralWave(
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
{

   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   ASSERT (signalvec,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signalvec->length >= 2, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (signalvec->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);


   ASSERT((INT4)params->approximant >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT((INT4)params->approximant < NumApproximants, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT((UINT4)params->order < (UINT4)LAL_PNORDER_NUM_ORDER,
            status, LALINSPIRALH_EORDER, LALINSPIRALH_MSGEORDER);

   if ( XLALSimInspiralImplementedTDApproximants(params->approximant) )
   {
      REAL8TimeSeries *hplus = NULL;
      REAL8TimeSeries *hcross = NULL;
      unsigned int idx;

      /* generate hplus and hcross */
      if (XLALSimInspiralChooseWaveformFromInspiralTemplate(&hplus, &hcross, params) == XLAL_FAILURE)
         ABORTXLAL(status);

      /* check length of waveform compared to signalvec */
      if (hplus->data->length > signalvec->length)
      {
         XLALDestroyREAL8TimeSeries(hplus);
         XLALDestroyREAL8TimeSeries(hcross);
         XLALPrintError("XLALSimInspiralChooseWaveformFromInspiralTemplate generated a waveform longer than output vector.\n");
         ABORT(status, LALINSPIRALH_EVECTOR, LALINSPIRALH_MSGEVECTOR);
      }

      /* initialize template vector to zero since hplus may be shorter */
      memset(signalvec->data, 0.0, signalvec->length * sizeof(signalvec->data[0]));

      /* convert REAL8 hplus waveform to REAL4 and store it in signalvec->data */
      for(idx = 0; idx < signalvec->length && idx < hplus->data->length ; idx++)
      {
          signalvec->data[idx] = (REAL4) hplus->data->data[idx];
      }

      params->tC = - ( ((REAL8)hplus->epoch.gpsSeconds) + 1.e-9*((REAL8)hplus->epoch.gpsNanoSeconds) );

      /* free hplus and hcross */
      XLALDestroyREAL8TimeSeries(hplus);
      XLALDestroyREAL8TimeSeries(hcross);
   }
   else switch (params->approximant)
   {
      case TaylorT1:
      case PadeT1:
           if (XLALInspiralWave1(signalvec, params) == XLAL_FAILURE) ABORTXLAL(status);
	   break;
      case TaylorT2:
           if (XLALInspiralWave2(signalvec, params) == XLAL_FAILURE) ABORTXLAL(status);
	   break;
      case TaylorT3:
           if (XLALInspiralWave3(signalvec, params) == XLAL_FAILURE) ABORTXLAL(status);
	   break;
      case EOB:
      case EOBNR:
           if (XLALEOBWaveform(signalvec, params) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case EOBNRv2:
      case EOBNRv2HM:
           if (XLALEOBPPWaveform(signalvec, params) == XLAL_FAILURE) ABORTXLAL(status);
	   break;
      case IMRPhenomA:
      case IMRPhenomB:
           if (XLALBBHPhenWaveTimeDom(signalvec, params) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case IMRPhenomFA:
           if (XLALBBHPhenWaveAFreqDom(signalvec, params) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case IMRPhenomFB:
           if (XLALBBHPhenWaveBFreqDom(signalvec, params) == XLAL_FAILURE) ABORTXLAL(status);
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
	   if (XLALInspiralStationaryPhaseApprox1(signalvec, params) == XLAL_FAILURE) ABORTXLAL(status);
	   break;
      case TaylorF2:
      case FindChirpSP:
           if (XLALInspiralStationaryPhaseApprox2(signalvec, params) == XLAL_FAILURE) ABORTXLAL(status);
	   break;
      case TaylorF2RedSpin:
           if (XLALTaylorF2ReducedSpin(signalvec, params) == XLAL_FAILURE) ABORTXLAL(status);
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
      case SpinTaylorFrameless:
           if (XLALSTPNFramelessWaveform(signalvec, params) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case PhenSpinTaylorRD:
           if (XLALPSpinInspiralRD(signalvec, params) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case SpinQuadTaylor:
           if (XLALSQTPNWaveform(signalvec, params) == XLAL_FAILURE) ABORTXLAL(status);
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
   	   if (XLALTaylorEtWaveform(signalvec, params) == XLAL_FAILURE) ABORTXLAL(status);
      	   break;
      case TaylorT4:
   	   LALTaylorT4Waveform(status->statusPtr, signalvec, params);
	   CHECKSTATUSPTR(status);
      	   break;
      case TaylorN:
   	   if (XLALTaylorNWaveform(signalvec, params) == XLAL_FAILURE) ABORTXLAL(status);
      	   break;
      default:
           ABORT( status, LALINSPIRALH_ESWITCH, LALINSPIRALH_MSGESWITCH );
   }

   DETATCHSTATUSPTR(status);
   RETURN (status);
}

void
LALInspiralWaveTemplates(
   LALStatus        *status,
   REAL4Vector      *signalvec1,
   REAL4Vector      *signalvec2,
   InspiralTemplate *params
   )
{

   INITSTATUS(status);
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

   if ( XLALSimInspiralImplementedTDApproximants(params->approximant) )
   {
      REAL8TimeSeries *hplus = NULL;
      REAL8TimeSeries *hcross = NULL;
      unsigned int idx;

      /* generate hplus and hcross */
      if (XLALSimInspiralChooseWaveformFromInspiralTemplate(&hplus, &hcross, params) == XLAL_FAILURE)
         ABORTXLAL(status);

      /* check length of waveform compared to signalvec */
      if (hplus->data->length > signalvec1->length)
      {
         XLALDestroyREAL8TimeSeries(hplus);
         XLALDestroyREAL8TimeSeries(hcross);
         XLALPrintError("XLALSimInspiralChooseWaveformFromInspiralTemplate generated a waveform longer than output vector.\n");
         ABORT(status, LALINSPIRALH_EVECTOR, LALINSPIRALH_MSGEVECTOR);
      }

      /* initialize template vectors to zero since hplus and hcross may be shorter */
      memset(signalvec1->data, 0.0, signalvec1->length * sizeof(signalvec1->data[0]));
      memset(signalvec2->data, 0.0, signalvec2->length * sizeof(signalvec2->data[0]));

      /* convert REAL8 hplus and hcross waveforms to REAL4 and store them in signalvec[1,2]->data */
      for(idx = 0; idx < signalvec1->length; idx++)
      {
          signalvec1->data[idx] = (REAL4) hplus->data->data[idx];
          signalvec2->data[idx] = (REAL4) hcross->data->data[idx];
      }

      params->tC = - ( ((REAL8)hplus->epoch.gpsSeconds) + 1.e-9*((REAL8)hplus->epoch.gpsNanoSeconds) );

      /* free hplus and hcross */
      XLALDestroyREAL8TimeSeries(hplus);
      XLALDestroyREAL8TimeSeries(hcross);
   }
   else switch (params->approximant)
   {
      case TaylorT1:
      case PadeT1:
           if (XLALInspiralWave1Templates(signalvec1, signalvec2, params) == XLAL_FAILURE) ABORTXLAL(status);
      	   break;
      case TaylorT2:
           if (XLALInspiralWave2Templates(signalvec1, signalvec2, params) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case TaylorT3:
           if (XLALInspiralWave3Templates(signalvec1, signalvec2, params) == XLAL_FAILURE) ABORTXLAL(status);
      	   break;
      case TaylorEt:
           if (XLALTaylorEtWaveformTemplates(signalvec1, signalvec2, params) == XLAL_FAILURE) ABORTXLAL(status);
      	   break;
      case EOB:
      case EOBNR:
           if (XLALEOBWaveformTemplates(signalvec1, signalvec2, params) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case EOBNRv2:
      case EOBNRv2HM:
           if (XLALEOBPPWaveformTemplates(signalvec1, signalvec2, params) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case IMRPhenomA:
      case IMRPhenomB:
           if (XLALBBHPhenWaveTimeDomTemplates(signalvec1, signalvec2, params) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case IMRPhenomFA:
           if (XLALBBHPhenWaveAFreqDomTemplates(signalvec1, signalvec2, params) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case IMRPhenomFB:
           if (XLALBBHPhenWaveBFreqDomTemplates(signalvec1, signalvec2, params) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case TaylorF2RedSpin:
           if (XLALTaylorF2ReducedSpinTemplates(signalvec1, signalvec2, params) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case TaylorF1:
      case TaylorF2:
      case FindChirpSP:
      case PadeF1:
      case BCV:
      case BCVSpin:
           ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
	   break;
      case SpinTaylor:
           LALSTPNWaveformTemplates(status->statusPtr, signalvec1, signalvec2, params);
           CHECKSTATUSPTR(status);
           break;
      case SpinTaylorFrameless:
           if (XLALSTPNFramelessWaveformTemplates(signalvec1, signalvec2, params) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case PhenSpinTaylorRD:
	   if (XLALPSpinInspiralRDTemplates(signalvec1, signalvec2, params) == XLAL_FAILURE) ABORTXLAL(status);
	   break;
      case SpinQuadTaylor:
           LALSTPNWaveformTemplates(status->statusPtr, signalvec1, signalvec2, params);
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

void
LALInspiralWaveForInjection(
   LALStatus        *status,
   CoherentGW       *waveform,
   InspiralTemplate *inspiralParams,
   PPNParamStruc    *ppnParams)
{

   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   ASSERT (inspiralParams,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (ppnParams,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);


   if ( XLALSimInspiralImplementedTDApproximants(inspiralParams->approximant) )
   {
      REAL8TimeSeries *hplus = NULL;
      REAL8TimeSeries *hcross = NULL;
      unsigned int idx;

      /* generate hplus and hcross */
      if (XLALSimInspiralChooseWaveformFromInspiralTemplate(&hplus, &hcross, inspiralParams) == XLAL_FAILURE)
        ABORTXLAL(status);

      /* allocate the waveform data stucture h in the CoherentGW structure */
      waveform->h = (REAL4TimeVectorSeries *) LALMalloc( sizeof(REAL4TimeVectorSeries) );
      if ( waveform->h == NULL )
      {
         XLALDestroyREAL8TimeSeries(hplus);
         XLALDestroyREAL8TimeSeries(hcross);
         ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
      }
      memset( waveform->h, 0, sizeof(REAL4TimeVectorSeries) );

      /* populate the waveform data stucture h in the CoherentGW structure */
      waveform->h->data = XLALCreateREAL4VectorSequence(hplus->data->length, 2);
      if (waveform->h == NULL)
      {
         CHAR warnMsg[1024];
         XLALDestroyREAL8TimeSeries(hplus);
         XLALDestroyREAL8TimeSeries(hcross);
         LALFree(waveform->h);
         snprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
             "Memory allocation error when allocating CoherentGW REAL4VectorSequence.\n");
         LALInfo( status, warnMsg );
         ABORT( status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM );
      }
      waveform->h->f0 = hplus->f0;
      waveform->h->deltaT = hplus->deltaT;
      waveform->h->epoch = hplus->epoch;
      waveform->h->sampleUnits = hplus->sampleUnits;
      waveform->position = ppnParams->position;
      waveform->psi = ppnParams->psi;

      /* convert REAL8 hplus and hcross waveforms to REAL4 and store them in waveform->h->data */
      for(idx = 0; idx < waveform->h->data->length; idx++)
      {
         waveform->h->data->data[idx * 2] = (REAL4) hplus->data->data[idx];
         waveform->h->data->data[idx * 2 + 1] = (REAL4) hcross->data->data[idx];
      }

      /* free hplus and hcross */
      XLALDestroyREAL8TimeSeries(hplus);
      XLALDestroyREAL8TimeSeries(hcross);
   }
   else switch (inspiralParams->approximant)
   {
      case TaylorT1:
      case PadeT1:
           if (XLALInspiralWave1ForInjection(waveform, inspiralParams, ppnParams) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case TaylorT2:
           if (XLALInspiralWave2ForInjection(waveform, inspiralParams, ppnParams) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case TaylorT3:
           if (XLALInspiralWave3ForInjection(waveform, inspiralParams, ppnParams) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case EOB:
      case EOBNR:
           if (XLALEOBWaveformForInjection(waveform, inspiralParams, ppnParams) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case EOBNRv2:
      case EOBNRv2HM:
           if (XLALEOBPPWaveformForInjection(waveform, inspiralParams, ppnParams) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case IMRPhenomA:
      case IMRPhenomB:
           if (XLALBBHPhenWaveTimeDomForInjection(waveform, inspiralParams, ppnParams) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case BCV:
           /*LALBCVWaveformForInjection(status->statusPtr, waveform, inspiralParams, ppnParams);*/
      case BCVSpin:
      case TaylorF1:
      case TaylorF2:
      case TaylorF2RedSpin:
      case FindChirpSP:
      case PadeF1:
           ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
           break;
      case SpinTaylorFrameless:
           if (XLALSTPNFramelessWaveformForInjection(waveform, inspiralParams, ppnParams) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case SpinTaylor:
           LALSTPNWaveformForInjection(status->statusPtr, waveform, inspiralParams, ppnParams);
           CHECKSTATUSPTR(status);
           break;
      case PhenSpinTaylorRD:
           if (XLALPSpinInspiralRDForInjection(waveform, inspiralParams, ppnParams) == XLAL_FAILURE) ABORTXLAL(status);
           break;
      case SpinQuadTaylor:
           if (XLALSQTPNWaveformForInjection(waveform, inspiralParams, ppnParams) == XLAL_FAILURE) ABORTXLAL(status);
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

