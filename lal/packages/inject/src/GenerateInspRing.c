/** \file GenerateInspRing.c
 *  \ingroup GenerateInspRing
 *  \author S.Fairhurst
 * 
 *  \brief Functions for adding a (realistic?) merger ringdown to the end of
 *  and inspiral waveform
 *
 * $Id$ 
 *
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/Units.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GenerateRing.h>
#include <lal/Ring.h>

/** Takes an inspiral waveform, and a simInspiralTable and generates a ringdown
 * at an appropriate frequency and quality */
CoherentGW *
XLALGenerateInspRing(
    CoherentGW		    *waveform,    /**< the inspiral waveform */
    SimInspiralTable	*inspiralInj  /**< details of the inspiral */
    )
{
  static const char *func = "XLALClusterSnglInspiralTable";

  REAL4 *a, *f, *shift; /* pointers to generated amplitude and frequency data */
  REAL8 *phi;   /* pointer to generated phase data */

  /* length of the waveform */
  INT4 inputLength;
  INT4 outputLength;
  INT4 linearLength;
  INT4 ringLength;

  /* waveform parameters */
  REAL8 phase;
  REAL4 ampPlus;
  REAL4 ampPlusDot;
  REAL4 a2Plus;
  REAL4 ampCross;
  REAL4 ampCrossDot;
  REAL4 a2Cross;
  REAL4 freq;
  REAL4 freqDot;
  REAL4 polarization;
  REAL4 polDot;
  REAL4 dampFac;
  REAL4 phaseFac;
  REAL4 A, lambda;
  INT4 n, N;

  /* ringdown details */
  SimRingdownTable     *ringInj = NULL;
  REAL4                 mTot = 0;
  REAL4                 orbAngMom, totalAngMom;
  REAL4                 Jx, Jy, Jz, tmp_x;


  /*
   *
   * Work out where the inspiral ended
   *
   */

  /* check that the inspiral waveform exists */
  if ( !waveform || !(waveform->a->data) || !(waveform->f->data) ||
      !(waveform->phi->data) )
  {
    XLALPrintError("Invalid waveform passed as input to %s\n", func);
    XLAL_ERROR_NULL(func,XLAL_EIO);
  }

  if ( (waveform->a->data->length < 2) )
  {
    XLALPrintError(
        "Length of waveform input to %s must be at least 2 points\n", func);
    XLAL_ERROR_NULL(func,XLAL_EIO);
  }

  if ( !inspiralInj )
  {
    XLALPrintError("No sim inspiral table passed as input to %s\n", func);
    XLAL_ERROR_NULL(func,XLAL_EIO);
  }

  /* read in the number of points already present */
  inputLength = waveform->a->data->length;

  /* record the final frequency, amplitude and phase */
  ampPlus = waveform->a->data->data[2*inputLength - 2];
  ampPlusDot = ampPlus - waveform->a->data->data[2*inputLength - 4];

  ampCross = waveform->a->data->data[2*inputLength - 1];
  ampCrossDot = ampCross - waveform->a->data->data[2*inputLength - 3];

  freq = waveform->f->data->data[inputLength - 1];
  freqDot = freq - waveform->f->data->data[inputLength - 2];

  phase = waveform->phi->data->data[inputLength - 1];

  polarization = waveform->shift->data->data[inputLength - 1];
  polDot = polarization - waveform->shift->data->data[inputLength - 2];


  /*
   *
   * Work out where we want the ringdown to start
   *
   */

  /* Calculate the Ringdown parameters */
  XLALPrintInfo( "Calculating (approximate) parameters for the ringdown\n" );

  ringInj = (SimRingdownTable *) XLALCalloc( 1, sizeof(SimRingdownTable) );
  if ( ! ringInj )
  {
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );
  }

  memcpy( ringInj->waveform, "ringdown", LIGOMETA_WAVEFORM_MAX * sizeof(CHAR));
  LALSnprintf( ringInj->coordinates, LIGOMETA_COORDINATES_MAX * sizeof(CHAR), 
      "EQUATORIAL");
  ringInj->geocent_start_time = inspiralInj->geocent_end_time;
  ringInj->geocent_start_time.gpsNanoSeconds;

  ringInj->longitude    = inspiralInj->longitude;
  ringInj->latitude     = inspiralInj->latitude;
  ringInj->distance     = inspiralInj->distance;

  XLALPrintInfo(
      "Ringdown longitude = %.2f, latitude = %.2f, distance = %.2f\n",
      ringInj->longitude, ringInj->latitude, ringInj->distance );


  /* estimate the final angular momentum */
  orbAngMom = 4 * inspiralInj->eta * 0.7;
  mTot = inspiralInj->mass1 + inspiralInj->mass2;

  Jx = orbAngMom * sin(inspiralInj->inclination) + 
    (inspiralInj->spin1x * inspiralInj->mass1 * inspiralInj->mass1 + 
     inspiralInj->spin2x * inspiralInj->mass2 * inspiralInj->mass2) / 
    (mTot * mTot) ;
  Jy = (inspiralInj->spin1y * inspiralInj->mass1 * inspiralInj->mass1 + 
      inspiralInj->spin2y * inspiralInj->mass2 * inspiralInj->mass2) / 
    (mTot * mTot) ;
  Jz = orbAngMom * cos(inspiralInj->inclination) + 
    (inspiralInj->spin1z * inspiralInj->mass1 * inspiralInj->mass1 + 
     inspiralInj->spin2z * inspiralInj->mass2 * inspiralInj->mass2) / 
    (mTot * mTot);
  totalAngMom = pow(Jx * Jx + Jy * Jy + Jz * Jz, 0.5);

  XLALPrintInfo( "Approximate orbital ang mom = %.2f\n", orbAngMom);
  XLALPrintInfo( "Approximate total angular momentum of inspiral:"
      "Jx = %.2f, Jy = %.2f, Jz = %.2f\n", Jx, Jy, Jz);
  XLALPrintInfo( "Estimated Final Spin = %.2f\n", totalAngMom );

  /* estimate the final mass */
  XLALPrintInfo( "Total inspiral mass = %.2f\n", mTot );
  mTot *= (1 - 0.01 * (1.0 + 6.0 * totalAngMom * totalAngMom ) );
  ringInj->mass = mTot;
  XLALPrintInfo( "Estimated Final Mass = %.2f\n", ringInj->mass );

  /* populate the ring params */
  ringInj->spin = totalAngMom < 0.99 ? totalAngMom : 0.95;

  ringInj->frequency = pow( LAL_C_SI, 3) / LAL_G_SI / LAL_MSUN_SI / 2.0 / 
    LAL_PI * ( 1.0 - 0.63 * pow( ( 1.0 - ringInj->spin ), 0.3 ) ) / 
    ringInj->mass;
  ringInj->quality = 2.0 * pow( ( 1.0 - ringInj->spin ), -0.45 );

  XLALPrintInfo( "Estimated Frequency = %.2f, Quality = %.2f\n", 
      ringInj->frequency, ringInj->quality );

  /* XXX For now, we don't acutally correctly calculate the ring epsilon, 
   * amplitude or hrss
   ringInj->epsilon = 0.01 * (16 * inspiralInj->eta * inspiralInj->eta); 
   ringInj->amplitude = XLALBlackHoleRingAmplitude( ringInj->frequency,
   ringInj->quality, ringInj->distance, ringInj->epsilon );
   ringInj->hrss = ringInj->amplitude * sqrt( 2 / LAL_PI / ringInj->frequency ) 
   * pow( ( 2.0 * pow( ringInj->quality, 3.0 ) + ringInj->quality ) / 
   ( 1.0 + 4.0 * pow ( ringInj->quality, 2 ) ) , 0.5); XXX */

  /*
   *
   * Work out the length of the full signal
   *
   */

  /* We extend the frequency as follows:
   * 1) Increase the frequency linearly to time T s.t.
   *      frequency = 0.75 * ringdown frequency
   * 2) After that, exponentially approach the ring frequency.
   *    -- The exponential decay is determined by matching the derivative
   */

  /* We extend the amplitude as follows:
   * 1) Use a parabolic fit to the amplitude which fits 3 data points:
   * f(0), f'(0) from inspiral f'(2T) / f(2T) from ringdown
   */

  /* calculate the number of points needed to get to 0.75 * ringFreq */
  linearLength = ceil( (0.75 * ringInj->frequency - freq) / freqDot);

  /* calculate number of additional points necessary to damp by exp(12). */
  phaseFac = LAL_TWOPI * ringInj->frequency * waveform->a->deltaT;
  dampFac = exp( - phaseFac / (2 * ringInj->quality));
  ringLength = ceil((24 * ringInj->quality)/( phaseFac));  

  outputLength = inputLength + 2 * linearLength + ringLength;

  /* extend the data structures */
  waveform->a->data->length = outputLength;
  waveform->a->data->data = (REAL4 *) XLALRealloc( ( waveform->a->data->data ), 
      2*outputLength*sizeof(REAL4) );  
  if ( !waveform->a->data->data )
  {
    XLALFree( waveform->a);
    XLAL_ERROR_NULL( "XLALCalculateNRStrain", XLAL_ENOMEM );
  }

  memset(waveform->a->data->data + inputLength, 0, 
      (outputLength - inputLength) * sizeof(REAL4 *) );
  XLALResizeREAL8TimeSeries( waveform->phi, 0, outputLength);
  if ( !waveform->phi->data )
  {
    XLALFree( waveform->a);
    XLALFree( waveform->phi);
    XLAL_ERROR_NULL( "XLALCalculateNRStrain", XLAL_ENOMEM );
  }

  XLALResizeREAL4TimeSeries( waveform->f, 0, outputLength);
  if ( !waveform->f->data->data )
  {
    XLALFree( waveform->a );
    XLALFree( waveform->phi );
    XLALFree( waveform->f );
    XLAL_ERROR_NULL( "XLALCalculateNRStrain", XLAL_ENOMEM );
  }

  if ( waveform->shift )
  {
    XLALResizeREAL4TimeSeries( waveform->shift, 0, outputLength);
    if ( !waveform->f->data->data )
    {
      XLALFree( waveform->a );
      XLALFree( waveform->phi );
      XLALFree( waveform->f );
      XLALFree( waveform->shift );
    }
  }

  a = &(waveform->a->data->data[2*inputLength]);
  phi = &(waveform->phi->data->data[inputLength]);
  f = &(waveform->f->data->data[inputLength]);

  if ( waveform->shift )
  {
    shift = &(waveform->shift->data->data[inputLength]);
  }

  /* 
   * set frequency, phase and shift (if needed) up to the start of ringdown
   */

  /* run frequency close to ring frequency */
  for ( n = 0; n < linearLength; n++ )
  {
    phase = *(phi++) = phase + LAL_TWOPI * freq * waveform->a->deltaT;
    freq = *(f++) = freq + freqDot;
    if ( shift )
    {
      polarization = *(shift++) = polarization + polDot;
    }
  }

  /* converge to ring frequency and constant polarization */
  /* fit to f = f_ring - A * exp( - lambda * n )
   * then A = freq + fDot
   *      A * lambda = fDot
   */
  A = ringInj->frequency - (freq + freqDot);
  lambda = freqDot / A;

  for ( n = 0; n < linearLength + ringLength; n++ )
  {
    phase = *(phi++) = phase + LAL_TWOPI * freq * waveform->a->deltaT;
    freq = *(f++) = ringInj->frequency - A * exp( - n * lambda );
    if ( shift )
    {
      polDot /= exp( - lambda );
      polarization = *(shift++) = polarization + polDot;
    }
  }

  /*
   * set amplitude up to the start of ringdown
   */

  /* we fit with a quadratic.  The three pieces of data we need to fit are:
   * 1) The initial amplitude 
   * 2) The initial amplitude derivative
   * 3) The final amplitude derivative / amplitude
   *
   * so, if A_i = a_0 + a_1 * i + a_2 * i^2,
   * 
   * a_0 = amp
   * a_1 = ampDot
   * 
   *     a_1 + 2 * N * a_2
   * -------------------------- = (dampFac - 1)
   * a_0 + N * a_1  + N^2 * a_2
   *
   *       (dampFac - 1) (a_0 + N * a_1) - a_1
   * a_2 = ----------------------------
   *       (1 - dampFac) * N^2 + 2 * N 
   */

  N = 2 * linearLength;

  ampPlus += ampPlusDot;
  ampCross += ampCrossDot;

  a2Plus = ((dampFac - 1) * ( ampPlus + N * ampPlusDot ) - ampPlusDot) /
    ( 2*N - N*N*(dampFac - 1) );
  a2Cross = ((dampFac - 1) * ( ampCross + N * ampCrossDot ) - ampCrossDot) /
    ( 2*N - N*N*(dampFac - 1) );


  /* quadratic part */
  for ( n = 0; n < N; n++ )
  {
    *(a++) = ampPlus + ampPlusDot * n + a2Plus * n * n;
    *(a++) = ampCross + ampCrossDot * n + a2Cross * n * n;
  }

  /* ringdown part */
  ampPlus = ampPlus + ampPlusDot * N + a2Plus * N * N;
  ampCross = ampCross + ampCrossDot * N + a2Cross * N * N;

  for ( n = 0; n < ringLength; n++ )
  {
    ampPlus = *(a++) = ampPlus * dampFac;
    ampCross = *(a++) = ampCross * dampFac;
  }

  XLALFree( ringInj );
  return( waveform );
}

