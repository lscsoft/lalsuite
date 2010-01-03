/* *  Copyright (C) 2007 Bernd Machenschalk, Stephen Fairhurst
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
#include <lal/Date.h>
#include <lal/DetResponse.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Units.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GenerateRing.h>
#include <lal/Ring.h>
#include <lal/TimeDelay.h>
#include <lal/TimeSeries.h>
#include <lal/GenerateInspRing.h>

NRCSID (GENERATEINSPRINGC,"$Id$");

static SimRingdownTable*
XLALDeriveRingdownParameters(
    SimRingdownTable    *ringInj
    );

/** Takes an inspiral waveform, and a simInspiralTable and generates a ringdown
 * at an appropriate frequency and quality */
CoherentGW *
XLALGenerateInspRing(
    CoherentGW          *waveform,     /**< the inspiral waveform */
    SimInspiralTable    *inspiralInj,  /**< details of the inspiral */
    SimRingdownTable    *ringInj,      /**< ringdown details to be populated */
    int                  injectSignalType /**< type of injection */
    )
{
  static const char *func = "XLALGenerateInspRing";

  REAL4 *a, *f; /* pointers to generated amplitude and frequency data */
  REAL4 *shift = NULL; /* pointers to generated shift (polarization) data */
  REAL8 *phi;   /* pointer to generated phase data */

  /* length of the waveform */
  INT4 inputLength;
  INT4 outputLength;
  INT4 mergerLength;
  INT4 ringLength;
  INT4 endMerger = 0;
  INT4 condt;

  /* waveform parameters */
  REAL8 phase;
  REAL8 mergerPhase;
  REAL8 phi0;
  REAL4 dt;
  REAL4 ampPlus;
  REAL4 ampPlusDot;
  REAL4 a2Plus;
  REAL4 ampCross;
  REAL4 ampCrossDot;
  REAL4 a2Cross;
  REAL4 freq;
  REAL4 freqDot;
  REAL4 freq0;
  REAL4 freqDot0;
  REAL4 polarization = 0;
  REAL4 polDot = 0;
  REAL4 dampFac;
  REAL4 phaseFac;
  REAL4 A, lambda;
  REAL4 tToRing;
  REAL4 tRing;
  REAL4 amp;
  INT4 n, N;

  /* ringdown details */
  REAL4                 mTot = 0;
  REAL4                 orbAngMom, totalAngMom;
  REAL4                 Jx, Jy, Jz;

  /*REAL4                 splus, scross, cosiotaA, cosiotaB;*/
  REAL4                 cosiota;

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
  dt = waveform->f->deltaT;

  ampPlus = waveform->a->data->data[2*inputLength - 2];
  ampPlusDot = ampPlus - waveform->a->data->data[2*inputLength - 4];

  ampCross = waveform->a->data->data[2*inputLength - 1];
  ampCrossDot = ampCross - waveform->a->data->data[2*inputLength - 3];

  freq0 = waveform->f->data->data[inputLength - 1];
  freqDot0 = freq0 - waveform->f->data->data[inputLength - 2];
  freqDot0 /= dt;

  phase = waveform->phi->data->data[inputLength - 1];


  if ( waveform->shift )
  {
    polarization = waveform->shift->data->data[inputLength - 1];
    polDot = polarization - waveform->shift->data->data[inputLength - 2];
  }

  XLALPrintInfo(
      "Final inspiral frequency = %.2f Hz, fdot = %.2e Hz/s\n"
      "\n", freq0, freqDot0);

  /*
   *
   * Work out where we want the ringdown to start
   *
   */



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
  XLALPrintInfo( "Approximate total angular momentum of inspiral:\n"
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

  XLALPrintInfo(
      "Estimated Ringdown Frequency = %.2f Hz, and Quality = %.2f\n"
      "\n", ringInj->frequency, ringInj->quality );

  ringInj->longitude    = inspiralInj->longitude;
  ringInj->latitude     = inspiralInj->latitude;
  ringInj->distance     = inspiralInj->distance;

  XLALPrintInfo(
      "Ringdown longitude = %.2f, latitude = %.2f, distance = %.2f\n",
      ringInj->longitude, ringInj->latitude, ringInj->distance );


  /*
   *
   * Work out the length of the full signal
   *
   */

  /* We use a very simple model to extend the frequency of the waveform.
   * This is done in two parts:
   *
   * 1) Continue the frequency by assuming that fdot is proportional f^(11/3)
   *
   *    using this, we obtain the frequency and phase at any later time
   *    given f, fdot at t=0 as:
   *
   *    f(t) = f(0) * ( 1 - (8 fdot(0) t / 3 f(0) ) )^(-3/8)
   *
   *    phi(t) = phi0 - (2 pi) * 3/5 * f(0)^2 / fdot(0) *
   *                     (1 - (8 fdot(0) t / 3 f(0) ) )^(5/8)
   *
   *    where phi0 = phi(0) + ( 2 pi ) * ( 3 f(0)^2 / 5 fdot(0) )
   *
   *
   * 2) At the point before the ring frequency would be exceeded,
   *    exponentially approach the ring frequency.
   *    -- The exponential decay is determined by matching the derivative
   */

  /* calculate the number of points to reach 0.9 * ringdown frequency */
  tToRing = 3 * freq0 / ( 8 * freqDot0 ) *
    (1 - pow( freq0 / (0.9 * ringInj->frequency), 8.0/3.0) );
  mergerLength = ceil( tToRing / dt ) - 1;
  phi0 = phase + LAL_TWOPI * 3 * freq0 * freq0 / (5 * freqDot0);

  mergerPhase = phi0 - phase
    - ( LAL_TWOPI) * (3.0 * freq0 * freq0 ) / ( 5.0 * freqDot0 ) *
    pow( 1 - ( 8.0 * freqDot0 * tToRing) / ( 3.0 * freq0 ) , 5.0/8.0 );


  XLALPrintInfo( "Adding by hand a merger and ringdown\n");
  XLALPrintInfo( "Time to reach 0.9 of ringdown frequency = %.4f seconds,\n"
      "%d data points, %.2f radians in GW phase\n"
      "We then add the same length to asymptote to ringdown values\n",
      tToRing, mergerLength, mergerPhase);

  if ( mergerPhase > 200 * LAL_TWOPI )
  {
    XLALPrintError("Failed to add a decent merger and ringdown\n"
        "The merger had a length of %.2f radians in GW phase (only allow 400 pi)\n"
        "Returning null from %s\n",
        mergerPhase, func);
    XLALFree( waveform->a );
    XLALFree( waveform->phi );
    XLALFree( waveform->f );
    XLALFree( waveform->shift );
    XLAL_ERROR_NULL(func,XLAL_EFAILED);
  }


  /* calculate number of additional points necessary to damp by exp(12). */
  phaseFac = LAL_TWOPI * ringInj->frequency * dt;
  dampFac = exp( - phaseFac / (2 * ringInj->quality));

  ringLength = ceil((24 * ringInj->quality)/( phaseFac));
  tRing = ringLength * dt;

  XLALPrintInfo( "Time to damp the ringdown by exp(12) = %.4f seconds,\n"
      "%d data points, %.2f radians in GW phase\n"
      "\n", tRing, ringLength, ringLength * phaseFac);

  /* Total extra length is merger length + ringdown length + some length over
   * which we smooth the merger to the ring */

  outputLength = inputLength + 2 * mergerLength - 1 + ringLength;




  /*
   *
   * extend the data structures
   *
   */

  waveform->a->data->length = outputLength;
  waveform->a->data->data = (REAL4 *) XLALRealloc( ( waveform->a->data->data ),
      2*outputLength*sizeof(REAL4) );
  if ( !waveform->a->data->data )
  {
    XLALFree( waveform->a);
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );
  }

  memset(waveform->a->data->data + 2 * inputLength, 0,
      2 * (outputLength - inputLength) * sizeof(REAL4) );
  XLALResizeREAL8TimeSeries( waveform->phi, 0, outputLength);
  if ( !waveform->phi->data )
  {
    XLALFree( waveform->a);
    XLALFree( waveform->phi);
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );
  }

  XLALResizeREAL4TimeSeries( waveform->f, 0, outputLength);
  if ( !waveform->f->data->data )
  {
    XLALFree( waveform->a );
    XLALFree( waveform->phi );
    XLALFree( waveform->f );
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );
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
      XLAL_ERROR_NULL( func, XLAL_ENOMEM );
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
   *
   * generate frequency, phase and shift (if needed) for merger and ringdown
   *
   */

  freq = freq0;

  /* run frequency close to ring frequency */
  for ( n = 1; n < mergerLength + 1; n++ )
  {
    REAL4 factor = 1 - (8.0 * freqDot0 * n * dt) / (3.0 * freq0 );
    freq = *(f++) = freq0 * pow( factor , -3.0/8.0 );
    phase = *(phi++) = phi0 -
      LAL_TWOPI * 3.0 * freq0 * freq0 / (5.0 * freqDot0) *
      pow( factor , 5.0/8.0);
    if ( shift )
    {
      polarization = *(shift++) = polarization + polDot;
    }
  }

  /* populate ringdown phase */
  ringInj->phase = phase;

  /* converge to ring frequency and constant polarization */
  /* fit to freq = f_ring - A * exp( - lambda * t )
   * then A = f_ring - freq
   *      A * lambda = freqDot
   */
  freqDot = (freq - *(f - 2)) / dt;
  A = ringInj->frequency - freq;
  lambda = freqDot / A;

  XLALPrintInfo(
      "Starting to asymptote to ringdown\n"
      "Final 'merger' frequency = %.2f Hz, fdot = %.2e Hz/s\n",
      freq, freqDot);
  XLALPrintInfo(
      "Frequency evolution fitted to freq = f_ring - A exp(-lambda t)\n"
      "A = %.2f, lambda = %.2e\n"
      "\n", A, lambda);

  condt = 0;
  for ( n = 1; n < mergerLength + ringLength; n++ )
  {
    phase = *(phi++) = phase + LAL_TWOPI * freq * dt;
    freq = *(f++) = ringInj->frequency - A * exp( - n * dt * lambda );
    if ( (freq == ringInj->frequency) && (condt == 0))
    {
      endMerger = n - 1.0;
      condt = 1.0;
    }
    if ( shift )
    {
      polDot *= exp( - lambda * dt);
      polarization = *(shift++) = polarization + polDot;
    }
  }

  /* ringdown start time */
  {
    REAL8 tRing2;
    tRing2 = XLALGPSGetREAL8(&(inspiralInj->geocent_end_time));
    tRing2 += ( mergerLength + endMerger + 1.0) * dt;
    XLALGPSSetREAL8( &(ringInj->geocent_start_time), tRing2 );
  }
  XLALPrintInfo(
      "Ringdown start time = %d.%d\n", ringInj->geocent_start_time.gpsSeconds,
      ringInj->geocent_start_time.gpsNanoSeconds);

  /* ringdown polarization */
  ringInj->polarization = polarization;


  /*
   *
   * set amplitude for merger and ringdown
   *
   */

  /* From end of inspiral to start of ringdown, we fit the amplitudes with
   * a quadratic.  The three pieces of data we need to fit are:
   * 1) The initial amplitude
   * 2) The initial amplitude derivative
   * 3) The ringdown amplitude derivative / amplitude =  1 - dampFac
   *
   * so, if A(n) = a_0 + a_1 * n + a_2 * n^2,
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

  /* the number of points */
  N = 2 * mergerLength;

  a2Plus = ((dampFac - 1) * ( ampPlus + N * ampPlusDot ) - ampPlusDot) /
    ( 2*N - N*N*(dampFac - 1) );
  a2Cross = ((dampFac - 1) * ( ampCross + N * ampCrossDot ) - ampCrossDot) /
    ( 2*N - N*N*(dampFac - 1) );

  XLALPrintInfo( "Fitting amplitude evolution to a quadratic\n"
      "A = a_0 + a_1 * t + a_2 * t^2\n"
      "For plus polarization, a_0 = %.2e, a_1 = %.2e, a_2 = %.2e\n"
      "For cross polarization, a_0 = %.2e, a_1 = %.2e, a_2 = %.2e\n"
      "\n", ampPlus, ampPlusDot / dt, a2Plus / (dt * dt),
      ampCross, ampCrossDot / dt, a2Cross / (dt * dt) );

  /* quadratic part */
  for ( n = 1; n < N; n++ )
  {
    *(a++) = ampPlus + ampPlusDot * n + a2Plus * n * n;
    *(a++) = ampCross + ampCrossDot * n + a2Cross * n * n;
  }

  /* set the final amplitudes */
  ampPlus = *(a - 2);
  ampCross = *(a - 1);

  for ( n = 0; n < ringLength; n++ )
  {
    ampPlus = *(a++) = ampPlus * dampFac;
    ampCross = *(a++) = ampCross * dampFac;
  }

  /* h0 */
  n = 2.0*(inputLength + mergerLength + endMerger );
  ampPlus = waveform->a->data->data[n];
  n = 2.0*(inputLength + mergerLength + endMerger ) + 1.0;
  ampCross = waveform->a->data->data[n];
  /* FIXME This calculation leads to nans and thus failures is there a typo?
  cosiotaA = ampPlus / ampCross *
    (1.0 + sqrt(1-ampCross*ampCross/ampPlus/ampPlus));
  cosiotaB = ampPlus / ampCross *
    (1.0 - sqrt(1-ampCross*ampCross/ampPlus/ampPlus));

  if ( cosiotaA > -1.0 && cosiotaA < 1.0)
    cosiota = cosiotaA;
  else if ( cosiotaB > -1.0 && cosiotaB < 1.0)
    cosiota = cosiotaB;
  else
  {
    XLALPrintError("inclination angle out of range\n");
    XLALFree( waveform->a );
    XLALFree( waveform->phi );
    XLALFree( waveform->f );
    XLALFree( waveform->shift );
    XLAL_ERROR_NULL(func,XLAL_EFAILED);
  }
  */
  /* ringdown inclination and amplitude */
  /*ringInj->inclination = acos( cosiota );*/
  /* Since the inclination angle doens't seem to be reliably calculated by
     the above calculation I am setting it to zero */
  cosiota = 0.0;
  ringInj->inclination = 0.0;
  amp =  ampPlus / ( 1.0 + cosiota * cosiota );
  ringInj->amplitude = sqrt( amp*amp);

  /* zero out inspiral and merger if we only want to inject a ringdown*/
  if ( injectSignalType ==  imr_ring_inject ||
      strstr(inspiralInj->waveform, "KludgeRingOnly") )
  {
    for ( n = 0; n < 2.0*(inputLength + mergerLength + endMerger ); n++ )
    {
      waveform->a->data->data[n]= 0;
    }
  }
  ringInj = XLALDeriveRingdownParameters(ringInj);

  return( waveform );
}

static SimRingdownTable*
XLALDeriveRingdownParameters(
    SimRingdownTable    *ringInj      /**< ringdown details to be populated */
    )
{
  /*static const char *func = "XLALDeriveRingdownParameters";*/
  REAL4 cosiota;
  double splus, scross, fplus, fcross;
  double tDelay, tStart;

  LALDetector lho = lalCachedDetectors[LALDetectorIndexLHODIFF];
  LALDetector llo = lalCachedDetectors[LALDetectorIndexLLODIFF];

  /* Calculate the Ringdown parameters */
  XLALPrintInfo( "Calculating (approximate) parameters for the ringdown\n" );

  /* waveform */
  memcpy( ringInj->waveform, "Ringdown",
      LIGOMETA_WAVEFORM_MAX * sizeof(CHAR));
  snprintf( ringInj->coordinates, LIGOMETA_COORDINATES_MAX * sizeof(CHAR),
      "EQUATORIAL");

  /* calculate hrss */
  ringInj->hrss = ringInj->amplitude * sqrt( 2 / LAL_PI / ringInj->frequency ) *
    pow( ( 2.0 * pow( ringInj->quality, 3.0 ) + ringInj->quality ) /
        ( 1.0 + 4.0 * pow ( ringInj->quality, 2 ) ) , 0.5);

  XLALPrintInfo( "Ringdown h_rss = %.2g\n", ringInj->hrss);

  /* and epsilon */
  ringInj->epsilon = XLALBlackHoleRingEpsilon( ringInj->frequency,
      ringInj->quality, ringInj->distance, ringInj->amplitude );
  XLALPrintInfo( "Ringdown epsilon = %.2g\n", ringInj->epsilon);

  /* calculate start time at sites */
  tStart =  XLALGPSGetREAL8(&(ringInj->geocent_start_time) );
  ringInj->start_time_gmst = XLALGreenwichSiderealTime(
      &(ringInj->geocent_start_time), 0.0 );

  /* lho */
  tDelay = XLALTimeDelayFromEarthCenter( lho.location, ringInj->longitude,
      ringInj->latitude, &(ringInj->geocent_start_time) );
  XLALGPSSetREAL8( &(ringInj->h_start_time), tDelay + tStart );
  /* llo */
  tDelay = XLALTimeDelayFromEarthCenter( llo.location, ringInj->longitude,
      ringInj->latitude, &(ringInj->geocent_start_time) );
  XLALGPSSetREAL8( &(ringInj->l_start_time), tDelay + tStart );

  XLALPrintInfo( "Site start times:\n"
      "LHO: %d.%d GPS seconds\n"
      "LLO: %d.%d GPS seconds\n",
      ringInj->h_start_time.gpsSeconds, ringInj->h_start_time.gpsNanoSeconds,
      ringInj->l_start_time.gpsSeconds, ringInj->l_start_time.gpsNanoSeconds);

  /* compute effective distance and hrss for sites */
  cosiota = cos(ringInj->inclination);
  splus = -( 1.0 + cosiota * cosiota );
  scross = -2.0 * cosiota;
  /* LHO */
  XLALComputeDetAMResponse(&fplus, &fcross, lho.response, ringInj->latitude,
      ringInj->longitude, ringInj->polarization, ringInj->start_time_gmst);
  ringInj->eff_dist_h = 2.0 * ringInj->distance;
  ringInj->eff_dist_h /= sqrt( splus*splus*fplus*fplus +
      scross*scross*fcross*fcross );
  ringInj->hrss_h = ringInj->amplitude * pow ( (
        (2*pow(ringInj->quality,3)+ringInj->quality ) *
        splus*splus*fplus*fplus +
        2*pow(ringInj->quality,2) * splus*scross*fplus*fcross +
        2*pow(ringInj->quality,3) * scross*scross*fcross*fcross )
      /  2.0 / LAL_PI / ringInj->frequency /
      ( 1.0 + 4.0 * pow ( ringInj->quality, 2 ) ) , 0.5 );
  XLALPrintInfo( "LHO response: F+ = %f, Fx = %f\n", fplus, fcross);

  /* compute hrss at LLO */
  /* LLO */
  XLALComputeDetAMResponse(&fplus, &fcross, llo.response, ringInj->longitude,
      ringInj->latitude, ringInj->polarization, ringInj->start_time_gmst);
  ringInj->eff_dist_l = 2.0 * ringInj->distance;
  ringInj->eff_dist_l /= sqrt( splus*splus*fplus*fplus +
      scross*scross*fcross*fcross );
  ringInj->hrss_l = ringInj->amplitude * pow ( (
        (2*pow(ringInj->quality,3)+ringInj->quality ) *
        splus*splus*fplus*fplus +
        2*pow(ringInj->quality,2) * splus*scross*fplus*fcross +
        2*pow(ringInj->quality,3) * scross*scross*fcross*fcross )
      /  2.0 / LAL_PI / ringInj->frequency /
      ( 1.0 + 4.0 * pow ( ringInj->quality, 2 ) ) , 0.5 );

  XLALPrintInfo( "LLO response: F+ = %f, Fx = %f\n", fplus, fcross);
  XLALPrintInfo( "Site effective distances:\n"
      "LHO: %f Mpc\n"
      "LLO: %f Mpc\n",
      ringInj->eff_dist_h, ringInj->eff_dist_l);

  XLALPrintInfo( "Site h_rss:\n"
      "LHO: %g\n"
      "LLO: %g\n",
      ringInj->hrss_l, ringInj->hrss_l);
  return( ringInj );
}

