/*
*  Copyright (C) 2007 Teviet Creighton
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
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/Units.h>
#include <lal/Inject.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/SimulateInspiral.h>

/* The lower cutoff frequency, if not specified as an input paramster,
   is defined as the point where the sensitivity is reduced by a
   factor of SIMULATEINSPIRALC_CUTOFF: */
#define SIMULATEINSPIRALC_CUTOFF (0.000001)

/**
 * \author Creighton, T. D.
 *
 * \brief Injects inspiral waveforms into detector output.
 *
 * ### Description ###
 *
 * <tt>LALSimulateInspiral()</tt>:
 * This function generates a binary inspiral signal using the parameters
 * in <tt>*params</tt>, simulates an instrument's response to that signal
 * using the instrument transfer function <tt>*transfer</tt>, and injects
 * the resulting waveform into detector output stored in <tt>*output</tt>.
 *
 * The <tt>*output</tt> time series should have all of its fields set to
 * their desired values, and should have a data sequence already
 * allocated; the function <tt>LALSimulateInspiral()</tt> simply adds the
 * inspiral waveform on top of the existing data. The \c epoch and
 * \c deltaT fields must be set, as they are used to determine the
 * sample rate and time positioning of the injected signal.  The
 * \c sampleUnits field must be set to \c lalADCCountUnit for
 * consistency.
 *
 * The <tt>*transfer</tt> frequency series should define the complex
 * frequency response function \f$T(f)\f$ taking the differential strain
 * signal \f$\tilde{h}(f)\f$ to detector response
 * \f$\tilde{o}(f)=T(f)\tilde{h}(f)\f$, and should have units of ADC counts
 * per strain.  It is treated as zero outside its frequency domain, and
 * is linearly interpolated between its frequency samples.
 *
 * The <tt>*params</tt> structure represents the parameters of an inspiral
 * signal to be injected (if <tt>params->next</tt>=\c NULL), or the
 * head of a linked list of parameter structures for multiple injections.
 * For each structure, if the \c signalAmplitude field is \f$\geq0\f$,
 * the injected waveform will be scaled to give it the correct
 * characteristic detection amplitude, and the \c effDist field is
 * set appropriately.  If \c signalAmplitude\f$<0\f$ and
 * \c effDist\f$>0\f$, the waveform is injected with that effective
 * distance, and the \c signalAmplitude field is set appropriately.
 * If \c signalAmplitude\f$<0\f$ and \c effDist\f$\leq0\f$, an error is
 * returned (that and all subsequent injections are skipped).
 *
 * An error is also returned (and no injections performed) if any of the
 * fields of <tt>*output</tt> and <tt>*transfer</tt> are not set to usable
 * values, including such things as wrong units or bad sampling
 * intervals.
 *
 * ### Usage ###
 *
 * One of the most useful applications of this routine is to generate
 * simulated noise containing a signal.  The following code snippet
 * generates white Gaussian noise with rms amplitude \c SIGMA, and
 * injects a signal with intrinsic signal-to-noise ratio
 * \f$\sqrt{(h|h)}=\f$\c SNR into it, coalescing at a time \c DT
 * seconds from the start of the time series, with a wave phase
 * \c PHI at coalescence.  The <tt>REAL4TimeSeries output</tt> and
 * <tt>COMPLEX8FrequencySeries transfer</tt> structures are assumed to be
 * defined and allocated outside of this block.
 *
 * \code
 * {
 * UINT4 i;
 * SimulateInspiralParamStruc inspParams;
 * RandomParams *randParams = NULL;
 *
 * // Generate white Gaussian noise.
 * LALCreateRandomParams( status->statusPtr, &randParams, 0 );
 * LALNormalDeviates( status->statusPtr, output.data, randParams );
 * for ( i = 0; i < output.data->length; i++ )
 * output.data->data[i] *= SIGMA;
 * LALDestroyRandomParams( status->statusPtr, &randParams );
 *
 * // Inject signal.
 * inspParams.timeC = output.epoch;
 * inspParams.timeC.gpsSeconds += DT;
 * inspParams.phiC = PHI; inspParams.mass1 = M1; inspParams.mass2 = M2;
 * inspParams.signalAmplitude = SNR*SIGMA;
 * inspParams.next = NULL;
 * LALSimulateInspiral( status->statusPtr, &output, &transfer, &inspParams );
 * }
 * \endcode
 *
 * ### Algorithm ###
 *
 * The default mode of operation, when one specifies the desired
 * amplitude, is as follows:
 *
 * First, <tt>LALGeneratePPNInspiral()</tt> is called to generate the
 * signal, placing the source at a distance of 1Mpc with optimal
 * orientation.  For lack of anything better, the amplitude and phase
 * functions are sampled at the full sampling interval as the output data
 * stream.  This function call also returns the time at coalescence.
 *
 * Second, the waveform produced by the signal in the detector output
 * stream is calculated.  The basic algorithm is the same as that in
 * <tt>LALSimulateCoherentGW()</tt>, but we can simplify it significantly
 * because we ignore polarization responses, time delays, and
 * interpolation between time samples.  Thus we have only to compute the
 * effect of the frequency transfer function \f${\cal T}(f)\f$.  As stated in
 * the \ref SimulateCoherentGW.h header, for quasiperiodic waveforms
 * \f$h(t)=\mathrm{Re}[{\cal H}(t)e^{i\phi(t)}]\f$ we can approximate the
 * instrument response (in the absence of noise) as:
 * \f[
 * o(t) \approx \mathrm{Re}[{\cal T}\{f(t)\}{\cal H}(t)e^{i\phi(t)}] \;.
 * \f]
 * In our case we are only sensitive to a single polarization (let's say
 * \f$h_+\f$), so we take \f${\cal H}(t)=A_+(t)\f$, where the phase of \f$\cal H\f$
 * is absorbed into the coalescence phase of \f$\phi\f$.  Then we can write
 * the instrument response as:
 * \f[
 * o(t) \approx A_+(t)[ T_\mathrm{re}\{f(t)\}\cos\phi(t)
 * - T_\mathrm{im}\{f(t)\}\sin\phi(t) ] \;.
 * \f]
 * This calculation can be done in place and stored in one of the arrays
 * for the amplitude, phase, or frequency functions, since they are
 * already sampled at the correct rate and have the correct length.
 *
 * Third, the characteristic detection amplitude is computed, and the
 * whole waveform is scaled so that it has the correct value.
 * Simultaneously, the effective distance is set to 1Mpc/(the scale factor).
 * The epoch is also adjusted to give the waveform the correct
 * coalescence time.
 *
 * Finally, LALSSInjectTimeSeries() is called to inject the
 * waveform into the output time series.  The whole procedure is repeated
 * for any other nodes in the linked list of parameters.
 *
 * If any parameter structure specifies the effective distance in place
 * of the characteristic detection amplitude, then the signal is injected
 * with that effective distance and is not rescaled.  The characteristic
 * detection amplitude field is set to the measured value.
 *
 */
void
LALSimulateInspiral( LALStatus                  *stat,
		     REAL4TimeSeries            *output,
		     COMPLEX8FrequencySeries    *transfer,
		     SimulateInspiralParamStruc *params )
{
  CHAR name[LALNameLength]; /* name of output time series */
  UINT4 i;                  /* an index */
  COMPLEX8 *tData;          /* pointer to transfer function data */
  PPNParamStruc ppnParams;  /* the parameters of the inspiral */
  CoherentGW signalvec;     /* the signal generated */

  /* Frequency interpolation constants.  By linear interpolation, the
     transfer function at any frequency f is given by:

                            f - f_k                  f_{k+1} - f
     T(f) = T(f_{k+1}) * -------------  +  T(f_k) * -------------
                         f_{k+1} - f_k              f_{k+1} - f_k

          = x*T_{k+1}  +  ( 1 - x )*T_k

     where:

     y = ( f - f0 )/deltaF = fOffset + dfInv*f
     k = (int)( y )
     x = y - k

  */
  REAL8 dfInv, fOffset;

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Make sure argument structures and their fields exist. */
  ASSERT( output, stat, SIMULATEINSPIRALH_ENUL,
	  SIMULATEINSPIRALH_MSGENUL );
  ASSERT( output->data, stat, SIMULATEINSPIRALH_ENUL,
	  SIMULATEINSPIRALH_MSGENUL );
  ASSERT( output->data->data, stat, SIMULATEINSPIRALH_ENUL,
	  SIMULATEINSPIRALH_MSGENUL );
  ASSERT( transfer, stat, SIMULATEINSPIRALH_ENUL,
	  SIMULATEINSPIRALH_MSGENUL );
  ASSERT( transfer->data, stat, SIMULATEINSPIRALH_ENUL,
	  SIMULATEINSPIRALH_MSGENUL );
  ASSERT( transfer->data->data, stat, SIMULATEINSPIRALH_ENUL,
	  SIMULATEINSPIRALH_MSGENUL );

  /* Set up interpolation constants for the transfer function. */
  ASSERT( transfer->deltaF != 0.0, stat, SIMULATEINSPIRALH_EDF,
	  SIMULATEINSPIRALH_MSGEDF );
  dfInv = 1.0/transfer->deltaF;
  fOffset = -dfInv*transfer->f0;

  /* If the low-frequency cutoff is not specified, find the point
     where the amplitude of the response drops to
     SIMULATEINSPIRALC_CUTOFF times its maximum. */
  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  if ( params->fStart > 0.0 )
    ppnParams.fStartIn = params->fStart;
  else {
    REAL4 xferMax = 0.0; /* maximum amplitude of transfer function */
    tData = transfer->data->data;
    i = transfer->data->length;
    while ( i-- ) {
      REAL4 xfer = fabs( crealf(*tData) ) + fabs( cimagf(*tData) );
      if ( xfer > xferMax )
	xferMax = xfer;
      tData++;
    }
    xferMax *= SIMULATEINSPIRALC_CUTOFF;
    tData = transfer->data->data;
    i = 0;
    while ( fabs( crealf(*tData) ) + fabs( cimagf(*tData) ) < xferMax ) {
      tData++;
      i++;
    }
    ppnParams.fStartIn = transfer->f0 + i*transfer->deltaF;
  }

  /* Set up other parameters that won't change between injections. */
  memset( &signalvec, 0, sizeof(CoherentGW) );
  ppnParams.deltaT = output->deltaT;
  tData = transfer->data->data;
  strncpy( name, output->name, LALNameLength );

  /* Inject for each node in the list of parameters. */
  while ( params ) {
    BOOLEAN fFlag = 0; /* whether we stepped outside transfer domain */
    INT8 tc;           /* coalescence time, in GPS nanoseconds */
    REAL4 amp2 = 0.0;  /* characteristic detection amplitude squared */
    REAL4 *aData;      /* pointer to amplitude data */
    REAL4 *fData;      /* pointer to frequency data */
    REAL8 *phiData;    /* pointer to phase data */

    /* First, generate the inspiral amplitude and phase as functions
       of time. */
    ppnParams.mTot = params->mass1 + params->mass2;
    ppnParams.eta = params->mass1*params->mass2
      /ppnParams.mTot/ppnParams.mTot;
    ppnParams.phi = params->phiC;
    ppnParams.d = 1000000.0*LAL_PC_SI;
    if ( params->signalAmplitude < 0.0 ) {
      if ( params->effDist <= 0.0 ) {
	ABORT( stat, SIMULATEINSPIRALH_EBAD, SIMULATEINSPIRALH_MSGEBAD );
      }
      ppnParams.d *= params->effDist;
    }
    TRY( LALGeneratePPNInspiral( stat->statusPtr, &signalvec,
				 &ppnParams ), stat );

    /* Next, simulate the instrument response.  To save memory, this
       will be done in place, overwriting the frequency function
       signalvec.f since it has the right type and length.  Also compute
       characteristic detection amplitude. */
    aData = signalvec.a->data->data;
    fData = signalvec.f->data->data;
    phiData = signalvec.phi->data->data;
    for ( i = 0; i < signalvec.f->data->length; i++ ) {
      REAL8 y = fOffset + *fData*dfInv;
      if ( y < 0.0 || y >= transfer->data->length ) {
	*fData = 0.0;
	fFlag = 1;
      } else {
	UINT4 k = (UINT4)( y );
	REAL8 x = y - k;
	REAL4 tRe = x*crealf(tData[k+1]) + ( 1.0 - x )*crealf(tData[k]);
	REAL4 tIm = x*cimagf(tData[k+1]) + ( 1.0 - x )*cimagf(tData[k]);
	*fData = *aData*( tRe*cos( *phiData ) - tIm*sin( *phiData ) );
	amp2 += (*fData)*(*fData);
      }
      aData+=2;
      fData++;
      phiData++;
    }
    signalvec.f->sampleUnits = lalADCCountUnit;

    /* Warn if we ever stepped outside of the frequency domain of the
       transfer function. */
    if ( fFlag )
      LALWarning( stat, "Signal passed outside of the frequency domain"
		  " of the transfer function (transfer function is"
		  " treated as zero outside its specified domain)" );

    /* Rescale either amplitude or distance, and shift time. */
    if ( params->signalAmplitude < 0.0 )
      params->signalAmplitude = sqrt( amp2 );
    else {
      amp2 = params->signalAmplitude / sqrt( amp2 );
      fData = signalvec.f->data->data;
      i = signalvec.f->data->length;
      while ( i-- ) {
	*fData *= amp2;
	fData++;
      }
      if ( amp2 == 0.0 )
	params->effDist = LAL_REAL4_MAX;
      else
	params->effDist = 1.0 / amp2;
    }
    tc = params->timeC.gpsSeconds;
    tc *= 1000000000L;
    tc += params->timeC.gpsNanoSeconds;
    tc -= (INT8)( 1000000000.0L*ppnParams.tc );
    signalvec.f->epoch.gpsSeconds = tc / 1000000000L;
    signalvec.f->epoch.gpsNanoSeconds = tc % 1000000000L;

    /* Inject the waveform into the output data, and reset everything
       for the next injection. */
    TRY( LALSDestroyVectorSequence( stat->statusPtr, &(signalvec.a->data) ),
	 stat );
    TRY( LALDDestroyVector( stat->statusPtr, &(signalvec.phi->data) ),
	 stat );
    LALFree( signalvec.a ); signalvec.a = NULL;
    LALFree( signalvec.phi ); signalvec.phi = NULL;
    LALSSInjectTimeSeries( stat->statusPtr, output, signalvec.f );
    BEGINFAIL( stat ) {
      TRY( LALSDestroyVector( stat->statusPtr, &(signalvec.f->data) ),
	   stat );
      LALFree( signalvec.f ); signalvec.f = NULL;
    } ENDFAIL( stat );
    TRY( LALSDestroyVector( stat->statusPtr, &(signalvec.f->data) ),
	 stat );
    LALFree( signalvec.f ); signalvec.f = NULL;
    strncpy( output->name, name, LALNameLength );
    params = params->next;
  }

  /* End of parameter list; no further cleanup should be necessary. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
