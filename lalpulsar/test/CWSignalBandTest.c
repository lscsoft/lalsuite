/*
 *
 * Copyright (C) 2018 Reinhard Prix
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

/*********************************************************************************/
/**
 * \author Reinhard Prix
 * \file
 * \brief Test for XLALCWSignal[Covering]Band().
 *
 */
#include <math.h>

#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/PulsarDataTypes.h>
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LogPrintf.h>

int test_SignalCoveringBand(UINT4 Ntrials);

int main(void)
{

  UINT4 seed = 1;
  srand ( seed );
  UINT4 Ntrials = 100;

  // test functions predicting CW signal bands
  // 1) XLALCWSignalCoveringBand(): total covering band for all signals over the sky and for given spin ranges
  // 2a) XLALSingleCWSignalCoveringBand(): actual covering band for a specific signal
  // 2b) XLALPrepareSingleCWSignalCoveringBand(): helper function compute detStates + 'worst' skyposition for widest Doppler band
  LIGOTimeGPS t0 = { 1206527521, 0};
  REAL8 Tspan = 1e6;
  LIGOTimeGPS t1 = t0;
  XLALGPSAdd ( &t1, Tspan );
  LIGOTimeGPS refTime = t0;
  XLALGPSAdd ( &t1, -Tspan );	// reftime outside of search-band for generality

  REAL8 asini = 0;
  REAL8 period = 0;
  REAL8 ecc = 0;
  LIGOTimeGPS tp = {0,0};
  REAL8 argp = 0;

  PulsarSpinRange spinRange;
  PulsarSpins XLAL_INIT_DECL(fkdot);
  fkdot[0] = 1000;
  fkdot[1] = -3e-9;
  fkdot[2] = 0;
  XLAL_CHECK_MAIN ( XLALInitPulsarSpinRangeFromSpins ( &spinRange, &refTime, fkdot, fkdot ) == XLAL_SUCCESS, XLAL_EFUNC );

  FILE *fp = stderr;
  if ( lalDebugLevel & LALINFOBIT ) {
    XLAL_CHECK_MAIN ( (fp = fopen ( "test_SignalCoveringBand.dat", "wb" )) != NULL, XLAL_EFUNC );
  }

  // ----- 1) get covering band for all signals over the sky for given spin-range
  REAL8 minCoverFreq, maxCoverFreq;
  REAL8 tic = XLALGetCPUTime();
  XLAL_CHECK_MAIN ( XLALCWSignalCoveringBand ( &minCoverFreq, &maxCoverFreq, &t0, &t1, &spinRange, asini, period, ecc ) == XLAL_SUCCESS, XLAL_EFUNC );
  REAL8 time_XLALCWSignalCoveringBand = XLALGetCPUTime() - tic;
  REAL8 coveringFreqBand = maxCoverFreq - minCoverFreq;
  fprintf ( fp, "coveringFreqBand = %g;\n", coveringFreqBand );

  // ----- 2a+b) get covering band for specific signals drawn randomly over the sky, check they satisfy 'width<=worst' and lie within total coverage band
  DetectorStateSeries *detStates;
  SkyPosition skypos_maxdoppler;
  REAL8 dT = 1800;
  LALDetector *detector;
  XLAL_CHECK_MAIN ( (detector = XLALGetSiteInfo ( "H1" )) != NULL, XLAL_EFUNC );
  EphemerisData *edat;
  XLAL_CHECK_MAIN ( (edat = XLALInitBarycenter ( TEST_DATA_DIR "earth00-19-DE405.dat.gz", TEST_DATA_DIR "sun00-19-DE405.dat.gz" )) != NULL, XLAL_EFUNC );
  tic = XLALGetCPUTime();
  XLAL_CHECK_MAIN ( (detStates = XLALPrepareCWSignalBand ( &skypos_maxdoppler, t0, Tspan, dT, detector, edat )) != NULL, XLAL_EFUNC );
  REAL8 time_XLALPrepareCWSignalBand = XLALGetCPUTime() - tic;
  XLALFree ( detector );
  REAL8 minFreq, maxFreq;
  PulsarDopplerParams XLAL_INIT_DECL(doppler);
  doppler.refTime = refTime;
  doppler.Alpha = skypos_maxdoppler.longitude;
  doppler.Delta = skypos_maxdoppler.latitude;
  memcpy ( doppler.fkdot, fkdot, sizeof(fkdot) );
  doppler.asini = asini;
  doppler.period = period;
  doppler.ecc = ecc;
  doppler.tp = tp;
  doppler.argp = argp;

  XLAL_CHECK_MAIN ( XLALCWSignalBand ( &minFreq, &maxFreq, detStates, &doppler ) == XLAL_SUCCESS, XLAL_EFUNC );
  REAL8 maxFreqBandWidth = maxFreq - minFreq;

  fprintf ( fp, "Alpha_max = %g; Delta_max = %g; minFreq = %g; maxFreq = %g;\nmaxFreqBandWidth = %g;\n",
            doppler.Alpha, doppler.Delta, minFreq, maxFreq, maxFreqBandWidth );
  XLAL_CHECK_MAIN ( maxFreqBandWidth <= coveringFreqBand, XLAL_EFAILED,
                    "(maxFreqBandWidth = %g) > (coveringFreqBand = %g)\n",
                    maxFreqBandWidth, coveringFreqBand);

  fprintf ( fp, "%%%% Alpha      Delta       minFreq    maxFreq   FreqBandWidth\n" );
  fprintf ( fp, "results = [\n");
  REAL8 time_XLALCWSignalBand = 0;
  for ( UINT4 i = 0; i < Ntrials; i ++ )
    {
      // draw random skyposition isotropically over the sky
      doppler.Alpha = LAL_TWOPI * (1.0 * rand() / ( RAND_MAX + 1.0 ) );  // alpha uniform in [0, 2pi)
      doppler.Delta = LAL_PI_2 - acos ( 1 - 2.0 * rand()/RAND_MAX );	// sin(delta) uniform in [-1,1]
      tic = XLALGetCPUTime();
      XLAL_CHECK_MAIN ( XLALCWSignalBand ( &minFreq, &maxFreq, detStates, &doppler ) == XLAL_SUCCESS, XLAL_EFUNC );
      time_XLALCWSignalBand += (XLALGetCPUTime() - tic);
      REAL8 FreqBandWidth = maxFreq - minFreq;
      XLAL_CHECK_MAIN ( FreqBandWidth <= maxFreqBandWidth, XLAL_EFAILED,
                   "Alpha = %g, Delta = %g: (FreqBandWidth = %g) > (maxFreqBandWidth = %g)\n",
                   doppler.Alpha, doppler.Delta, FreqBandWidth, maxFreqBandWidth );
      fprintf ( fp, "%10g, %10g %10.8g %10.8g %10.2e;\n", doppler.Alpha, doppler.Delta, minFreq, maxFreq, FreqBandWidth );
    }
  fprintf ( fp, "];\n");
  time_XLALCWSignalBand /= Ntrials;

  fprintf ( fp, "time_XLALCWSignalCoveringBand = %.2e;\n", time_XLALCWSignalCoveringBand );
  fprintf ( fp, "time_XLALPrepareCWSignalBand = %.2e;\n", time_XLALPrepareCWSignalBand );
  fprintf ( fp, "time_XLALCWSignalBand = %.2e;\n", time_XLALCWSignalBand );

  XLALDestroyDetectorStateSeries ( detStates );
  XLALDestroyEphemerisData ( edat );

  if ( fp )  {
    fclose ( fp );
  }

  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} /* main() */

