/*
*  Copyright (C) 2014 Reinhard Prix
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

#include <lal/XLALError.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/UniversalDopplerMetric.h>
#include <lal/LogPrintf.h>
#include <lal/CWMakeFakeData.h>
#include <lal/LALConstants.h>
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/ComputeFstat.h>
#include <lal/DetectorStates.h>

// basic consistency checks of the ComputeFstat module: compare F-stat results for all
// *available* Fstat methods against each other

static int compareFstatResults ( const FstatResults *result1, const FstatResults *result2 );

// ---------- main ----------
int
main ( int argc, char *argv[] )
{
  XLAL_CHECK ( argc == 1, XLAL_EINVAL, "No input arguments allowed.\n" );
  XLAL_CHECK ( argv != NULL, XLAL_EINVAL );


  EphemerisData *ephem;
  XLAL_CHECK ( (ephem = XLALInitBarycenter ( TEST_DATA_DIR "earth00-19-DE405.dat.gz", TEST_DATA_DIR "sun00-19-DE405.dat.gz" )) != NULL, XLAL_EFUNC );

  // ----- setup injection and data parameters
  UINT4 numDetectors = 2;

  // use signal-only injections to avoid usage of different noise bins to contaminate error-comparison
  MultiNoiseFloor XLAL_INIT_DECL(injectNoiseFloor);
  injectNoiseFloor.length = numDetectors;
  injectNoiseFloor.sqrtSn[0] = 1;
  injectNoiseFloor.sqrtSn[1] = 3;

  MultiNoiseFloor XLAL_INIT_DECL(assumeNoiseFloor);
  assumeNoiseFloor.length = numDetectors;
  assumeNoiseFloor.sqrtSn[0] = 1;
  assumeNoiseFloor.sqrtSn[1] = 3;


  LALStringVector *detNames = NULL;
  XLAL_CHECK ( (detNames = XLALCreateStringVector ( "H1", "L1", NULL )) != NULL, XLAL_EFUNC );
  MultiLALDetector XLAL_INIT_DECL(detInfo);
  XLAL_CHECK ( XLALParseMultiLALDetector ( &detInfo, detNames ) == XLAL_SUCCESS, XLAL_EFUNC );

  LIGOTimeGPS startTime = {711595934, 0};
  REAL8 Tspan = 100 * 3600;
  LIGOTimeGPS endTime = startTime;
  XLALGPSAdd( &endTime, Tspan );
  REAL8 Tsft = 1800;

  LIGOTimeGPS refTime = { startTime.gpsSeconds + round(0.5*Tspan), 0 };	// reftime in middle of segment

  MultiLIGOTimeGPSVector *multiTimestamps;
  XLAL_CHECK ( (multiTimestamps = XLALMakeMultiTimestamps ( startTime, Tspan, Tsft, 0, numDetectors )) != NULL, XLAL_EFUNC );

  SFTCatalog *catalog;
  XLAL_CHECK ( (catalog = XLALMultiAddToFakeSFTCatalog ( NULL, detNames, multiTimestamps )) != NULL, XLAL_EFUNC );

  // ----- CW sources to injet ----------
  REAL8 Freq = 100.0;

  PulsarParamsVector *sources;
  XLAL_CHECK ( (sources = XLALCreatePulsarParamsVector(1)) != NULL, XLAL_EFUNC );

  sources->data[0].Amp.h0   = 1;
  sources->data[0].Amp.cosi = 0.5;
  sources->data[0].Amp.psi  = 0.1;
  sources->data[0].Amp.phi0 = 1.2;

  REAL8 asini = 2.9;
  REAL8 Period = 10 * 3600;

  PulsarDopplerParams XLAL_INIT_DECL(Doppler);
  Doppler.Alpha = 0.5;
  Doppler.Delta = -0.5;
  Doppler.fkdot[0] = Freq;
  Doppler.fkdot[1] = -1e-9;
  Doppler.refTime = refTime;

  Doppler.asini = asini;
  Doppler.ecc = 0.1;
  Doppler.tp = startTime;
  Doppler.period = Period;
  Doppler.argp = 0.5;

  sources->data[0].Doppler = Doppler;

  REAL8 dFreq = 0.1 / Tspan;		// 10x finer than native FFT resolution
  REAL8 mis = 0.5;
  REAL8 df1dot = sqrt( 720.0 * mis ) / (LAL_PI * Tspan * Tspan);	// metric (f-projected) stepsize for given mismatch mis
  REAL8 dSky = 1e4 / (Freq * Tspan);	// rough estimate of a 'metric' sky step, eg. Eq.(118) in \cite Prix07

  UINT4 numFreqBins = 1000;
  UINT4 numf1dotBins = 2;
  UINT4 numSkyPoints = 2;

  PulsarSpinRange XLAL_INIT_DECL(spinRange);
  spinRange.refTime = refTime;
  memcpy ( &spinRange.fkdot, &sources->data[0].Doppler.fkdot, sizeof(spinRange.fkdot) );
  spinRange.fkdotBand[0] = (numFreqBins - 1)*dFreq - 10*LAL_REAL8_EPS;
  spinRange.fkdotBand[1] = (numf1dotBins - 1)*df1dot - 10*LAL_REAL8_EPS;

  REAL8 minCoverFreq, maxCoverFreq;
  XLAL_CHECK ( XLALCWSignalCoveringBand ( &minCoverFreq, &maxCoverFreq, &startTime, &endTime, &spinRange, asini, Period ) == XLAL_SUCCESS, XLAL_EFUNC );

  // ----- setup extra Fstat method params
  FstatExtraParams XLAL_INIT_DECL(extraParams);
  extraParams.randSeed  = 1;
  extraParams.SSBprec = SSBPREC_RELATIVISTICOPT;
  extraParams.Dterms = 8;	// constant value that works for all Demod methods

  // ----- prepare input data with injection for all available methods
  FstatInput *input[FMETHOD_END];
  FstatResults *results[FMETHOD_END];
  for ( UINT4 iMethod = FMETHOD_START; iMethod < FMETHOD_END; iMethod ++ )
    {
      results[iMethod] = NULL;
      if ( !XLALFstatMethodIsAvailable(iMethod) ) {
        continue;
      }
      input[iMethod] = XLALCreateFstatInput ( catalog, minCoverFreq, maxCoverFreq, sources, &injectNoiseFloor, &assumeNoiseFloor, 50, ephem, iMethod, &extraParams );
      XLAL_CHECK ( input[iMethod] != NULL, XLAL_EFUNC );
    }

  // ----- loop over all templates {sky, f1dot}
  for ( UINT4 iSky = 0; iSky < numSkyPoints; iSky ++ )
    {
      for ( UINT4 if1dot = 0; if1dot < numf1dotBins; if1dot ++ )
        {
          // ----- loop over all available methods and compare Fstat results
          FstatMethodType firstMethod = FMETHOD_START;
          for ( UINT4 iMethod = FMETHOD_START; iMethod < FMETHOD_END; iMethod ++ )
            {
              if ( !XLALFstatMethodIsAvailable(iMethod) ) {
                continue;
              }
              if ( firstMethod == FMETHOD_START ) {	// keep track of first available method found
                firstMethod = iMethod;
              }
              XLAL_CHECK ( XLALComputeFstat ( &results[iMethod], input[iMethod], &Doppler, dFreq, numFreqBins, FSTATQ_2F ) == XLAL_SUCCESS, XLAL_EFUNC );

              // compare to first result
              if ( compareFstatResults ( results[firstMethod], results[iMethod] ) != XLAL_SUCCESS )
                {
                  XLALPrintError ("Comparison between method '%s' and '%s' failed comparison\n",
                                  XLALGetFstatMethodName(firstMethod), XLALGetFstatMethodName(iMethod) );
                  XLAL_ERROR ( XLAL_EFUNC );
                }

            }  // for i < FMETHOD_END

          Doppler.fkdot[1] += df1dot;

        } // for if1dot < numf1dotBins

      Doppler.Alpha += dSky;

    } // for iSky < numSkyPoints

  // free remaining memory
  for ( UINT4 iMethod=FMETHOD_START; iMethod < FMETHOD_END; iMethod ++ )
    {
      if ( !XLALFstatMethodIsAvailable(iMethod) ) {
        continue;
      }
      XLALDestroyFstatInput ( input[iMethod] );
      XLALDestroyFstatResults ( results[iMethod] );
    } // for i < FMETHOD_END

  XLALDestroyPulsarParamsVector ( sources );
  XLALDestroySFTCatalog ( catalog );
  XLALDestroyMultiTimestamps ( multiTimestamps );
  XLALDestroyStringVector ( detNames );
  XLALDestroyEphemerisData ( ephem );
  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} // main()

static REAL8
relError ( REAL8 val1, REAL8 val2 )
{
  return ( (val1-val2) / ( 0.5 * (val1+val2) ) );
} // relError()

static int
compareFstatResults ( const FstatResults *result1, const FstatResults *result2 )
{
  XLAL_CHECK ( (result1 != NULL) && (result2 != NULL), XLAL_EINVAL );

  XLAL_CHECK ( result1->whatWasComputed == result2->whatWasComputed, XLAL_EINVAL );
  XLAL_CHECK ( result1->whatWasComputed == FSTATQ_2F, XLAL_EINVAL, "Currently on 2F comparison supported.\n");
  XLAL_CHECK ( result1->dFreq == result2->dFreq, XLAL_EINVAL );
  XLAL_CHECK ( result1->numFreqBins == result2->numFreqBins, XLAL_EINVAL );
  UINT4 numFreqBins = result1->numFreqBins;
  XLAL_CHECK ( result1->numDetectors == result2->numDetectors, XLAL_EINVAL );

  REAL8 tolMaxRelErr = 1.2;
  REAL8 tolMeanRelErr = 0.1;
  REAL8 tolMeanAbsRelErr = 0.15;
  REAL8 tolRelErrSum = 0.1;

  REAL8 maxRelErr = 0;
  REAL8 meanRelErr = 0;
  REAL8 meanAbsRelErr = 0;
  REAL8 relErrSum = 0;
  REAL8 sum1 = 0, sum2 = 0;
  REAL8 maxTwoF1 = 0, twoF2atMax = 0;
  for ( UINT4 k = 0; k < numFreqBins; k ++ )
    {
      REAL8 twoF1 = result1->twoF[k];
      REAL8 twoF2 = result2->twoF[k];

      // keep track of highest 2F_1 peak and corresponding 2F_2 value
      if ( twoF1 > maxTwoF1 ) {
        maxTwoF1 = twoF1;
        twoF2atMax = twoF2;
      }

      REAL8 relErr = relError ( twoF1, twoF2 );

      meanRelErr += relErr;
      meanAbsRelErr += fabs(relErr);
      maxRelErr = fmax ( maxRelErr, fabs ( relErr ) );
      sum1 += twoF1;
      sum2 += twoF2;

    } // for k < numFreqBins

  meanRelErr = fabs(meanRelErr) / numFreqBins;
  meanAbsRelErr = meanAbsRelErr / numFreqBins;
  relErrSum = fabs ( relError ( sum1, sum2 ) );

  REAL8 relErrAtMax = fabs ( relError ( maxTwoF1, twoF2atMax ) );
  REAL8 tolRelErrAtMax = 0.1;

  XLALPrintWarning ( "maxRelErr = %g (tol=%g), relErrAtMax = %g (tol=%g), meanRelErr = %g (tol=%g), meanAbsRelErr = %g (tol=%g), relErrSum = %g (tol=%g)\n",
                     maxRelErr, tolMaxRelErr, relErrAtMax, tolRelErrAtMax, meanRelErr, tolMeanRelErr, meanAbsRelErr, tolMeanAbsRelErr, relErrSum, tolRelErrSum );

  if ( (maxRelErr > tolMaxRelErr) || (relErrAtMax > tolRelErrAtMax) || (meanRelErr > tolMeanRelErr) || (meanAbsRelErr > tolMeanAbsRelErr) || (relErrSum > tolRelErrSum) )
    {
      XLALPrintError ( "maxRelErr = %g (tol=%g), relErrAtMax = %g (tol=%g), meanRelErr = %g (tol=%g), meanAbsRelErr = %g (tol=%g), relErrSum = %g (tol=%g)\n",
                       maxRelErr, tolMaxRelErr, relErrAtMax, tolRelErrAtMax, meanRelErr, tolMeanRelErr, meanAbsRelErr, tolMeanAbsRelErr, relErrSum, tolRelErrSum );

      FILE *fp = fopen("failed-fstats.dat", "wb");
      fprintf (fp, "%%%% Freq             F1        F2\n");
      for ( UINT4 k1 = 0; k1 < result1->numFreqBins; k1 ++ )
        {
          REAL8 fi = result1->doppler.fkdot[0] + k1 * result1->dFreq;
          fprintf ( fp, "%20.16g %10.5g %10.5g\n", fi, result1->twoF[k1], result2->twoF[k1] );
        }
      fclose(fp);

      XLAL_ERROR ( XLAL_ETOL );
    } // if error > tolerance

  return XLAL_SUCCESS;

} // compareFstatResults()
