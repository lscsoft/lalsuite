/*
 * Copyright (C) 2013 Reinhard Prix
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */


/**
 * \author Reinhard Prix
 * \ingroup CWMakeFakeData_h
 * \brief Functions to generate 'fake' data containing CW signals and/or Gaussian noise.
 * These basically present a high-level wrapper API to the lower-level CW signal-generation
 * functions in lalsuite.
 */

// ---------- includes
#include <math.h>

// GSL includes

// LAL includes
#include <lal/CWMakeFakeData.h>
#include <lal/TimeSeries.h>
#include <lal/GeneratePulsarSignal.h>
#include <lal/TransientCW_utils.h>
#include <lal/LALString.h>
#include <lal/StringVector.h>
#include <lal/Units.h>
#include <lal/ConfigFile.h>
#include <lal/LFTandTSutils.h>
#include <fftw3.h>
#include <lal/FFTWMutex.h>
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/ConfigFile.h>

// ---------- local defines

// ---------- local macro definitions
#define SQ(x) ( (x) * (x) )
// ---------- local type definitions

// ---------- Global variables
const REAL8 eps = 10 * LAL_REAL8_EPS;

const char *const InjectionSourcesHelpString = "Source parameters to inject for simulated signal(s).\n"
"This is a comma-separated list of file patterns for configuration files,\n"
"or else direct configuration strings in the following format:\n"
"  * Enclose with curly braces ('{}').\n"
"  * Give pulsar parameters as key=value pairs with a '=' separator.\n"
"  * Separate each key=value pair with a semicolon (';').\n"
"Available parameters are:\n"
"  * Required parameters: Alpha, Delta, Freq, refTime\n"
"  * Optional parameters:\n"
"    - Injection amplitudes: either (h0, cosi) or (aPlus, aCross), psi, phi0\n"
"    - Higher-order spindowns: f1dot, f2dot, ... f6dot\n"
"    - Binary sources: orbitTp, orbitArgp, orbitasini, orbitEcc, orbitPeriod\n"
"    - Transient injections: transientWindowType, transientStartTime, transientTau\n"
"Examples:\n"
"  * '{Alpha=0; Delta=0; Freq=50; f1dot=1e-11; f2dot=0; refTime=1000000000; h0=1.00000000e-23; cosi=0; psi=0; phi0=0;}'\n"
"  * 'file1.dat,someFiles*.txt,{Alpha=0;Delta=0;Freq=0;refTime=1000000000;},someOtherFiles[0-9].dat'\n\n";

// ---------- local prototypes
static UINT4 gcd (UINT4 numer, UINT4 denom);
int XLALcorrect_phase ( SFTtype *sft, LIGOTimeGPS tHeterodyne );
int XLALCheckConfigFileWasFullyParsed ( const char *fname, const LALParsedDataFile *cfgdata );

// ==================== FUNCTION DEFINITIONS ====================


/**
 * Generate fake 'CW' data, returned either as SFTVector or REAL4Timeseries or both,
 * for given CW-signal ("pulsar") parameters and output parameters (frequency band etc)
 */
int
XLALCWMakeFakeMultiData ( MultiSFTVector **multiSFTs,			///< [out] pointer to optional SFT-vector for output
                          MultiREAL8TimeSeries **multiTseries,		///< [out] pointer to optional timeseries-vector for output
                          const PulsarParamsVector *injectionSources,	///< [in] (optional) array of sources inject
                          const CWMFDataParams *dataParams,		///< [in] parameters specifying the type of data to generate
                          const EphemerisData *edat			///< [in] ephemeris data
                          )
{
  XLAL_CHECK ( (multiSFTs == NULL) || ((*multiSFTs) == NULL ), XLAL_EINVAL );
  XLAL_CHECK ( (multiTseries == NULL) || ((*multiTseries) == NULL ), XLAL_EINVAL );
  XLAL_CHECK ( (multiSFTs != NULL) || (multiTseries != NULL), XLAL_EINVAL );

  XLAL_CHECK ( dataParams != NULL, XLAL_EINVAL );
  XLAL_CHECK ( edat != NULL, XLAL_EINVAL );

  const MultiLIGOTimeGPSVector *multiTimestamps = dataParams->multiTimestamps;

  // check multi-detector input
  XLAL_CHECK ( dataParams->multiIFO.length >= 1, XLAL_EINVAL );
  UINT4 numDet = dataParams->multiIFO.length;
  XLAL_CHECK ( multiTimestamps->length == numDet, XLAL_EINVAL, "Inconsistent number of IFOs: detInfo says '%d', multiTimestamps says '%d'\n", numDet, multiTimestamps->length );
  XLAL_CHECK ( dataParams->multiNoiseFloor.length == numDet, XLAL_EINVAL );

  // check Tsft, consistent over detectors
  REAL8 Tsft = multiTimestamps->data[0]->deltaT;
  XLAL_CHECK ( Tsft > 0, XLAL_EINVAL, "Got invalid Tsft = %g must be > 0\n", Tsft );
  for ( UINT4 X=0; X < numDet; X ++ ) {
    XLAL_CHECK ( multiTimestamps->data[X]->deltaT == Tsft, XLAL_EINVAL, "Inconsistent Tsft, for Tsft[X=0]=%g, while Tsft[X=%d]=%g\n", Tsft, X, multiTimestamps->data[X]->deltaT );
  }

  // ----- prepare output containers, as required
  MultiSFTVector *outMSFTs = NULL;
  MultiREAL8TimeSeries *outMTS = NULL;
  if ( multiSFTs != NULL )
    {
      XLAL_CHECK ( (outMSFTs = XLALCalloc ( 1, sizeof(*outMSFTs) )) != NULL, XLAL_ENOMEM );
      XLAL_CHECK ( (outMSFTs->data = XLALCalloc ( numDet, sizeof(outMSFTs->data[0]))) != NULL, XLAL_ENOMEM );
      outMSFTs->length = numDet;
    } // if multiSFTs
  if ( multiTseries != NULL )
    {
      XLAL_CHECK ( (outMTS = XLALCalloc ( 1, sizeof(*outMTS) )) != NULL, XLAL_ENOMEM );
      XLAL_CHECK ( (outMTS->data = XLALCalloc ( numDet, sizeof(outMTS->data[0]))) != NULL, XLAL_ENOMEM );
      outMTS->length = numDet;
    } // if multiTseries

  for ( UINT4 X=0; X < numDet; X ++ )
    {
      SFTVector **svp = NULL;
      REAL8TimeSeries **tsp = NULL;
      if ( outMSFTs != NULL ) {
        svp = &(outMSFTs->data[X]);
      }
      if ( outMTS != NULL ) {
        tsp = &(outMTS->data[X]);
      }
      XLAL_CHECK ( XLALCWMakeFakeData ( svp, tsp, injectionSources, dataParams, X, edat ) == XLAL_SUCCESS, XLAL_EFUNC );

    } // for X < numDet

  // return multi-Timeseries if requested
  if ( multiTseries ) {
    (*multiTseries) = outMTS;
  }
  // return multi-SFTs if requested
  if ( multiSFTs ) {
    (*multiSFTs) = outMSFTs;
  }

  return XLAL_SUCCESS;

} // XLALCWMakeFakeMultiData()

/**
 * Single-IFO version of XLALCWMakeFakeMultiData(), handling the actual
 * work, but same input API. The 'detectorIndex' has the index of the detector
 * to be used from the multi-IFO arrays.
 */
int
XLALCWMakeFakeData ( SFTVector **SFTvect,
                     REAL8TimeSeries **Tseries,
                     const PulsarParamsVector *injectionSources,
                     const CWMFDataParams *dataParams,
                     UINT4 detectorIndex,	/* index for current detector in dataParams */
                     const EphemerisData *edat
                     )
{
  XLAL_CHECK ( (SFTvect == NULL) || ((*SFTvect) == NULL ), XLAL_EINVAL );
  XLAL_CHECK ( (Tseries == NULL) || ((*Tseries) == NULL ), XLAL_EINVAL );
  XLAL_CHECK ( (SFTvect != NULL) || (Tseries != NULL), XLAL_EINVAL );
  XLAL_CHECK ( edat != NULL, XLAL_EINVAL );

  XLAL_CHECK ( dataParams != NULL, XLAL_EINVAL );
  XLAL_CHECK ( detectorIndex < dataParams->multiIFO.length, XLAL_EINVAL );
  XLAL_CHECK ( detectorIndex < dataParams->multiNoiseFloor.length, XLAL_EINVAL );
  XLAL_CHECK ( detectorIndex < dataParams->multiTimestamps->length, XLAL_EINVAL );
  XLAL_CHECK ( (dataParams->inputMultiTS == NULL) || (detectorIndex < dataParams->inputMultiTS->length), XLAL_EINVAL );

  // initial default values fMin, sampling rate from caller input or timeseries
  REAL8 fMin  = dataParams->fMin;
  REAL8 fBand = dataParams->Band;
  REAL8 fSamp = 2.0 * fBand;
  if ( dataParams->inputMultiTS != NULL )
    {
      const REAL8TimeSeries *ts = dataParams->inputMultiTS->data[detectorIndex];
      XLAL_CHECK ( ts != NULL, XLAL_EINVAL );
      REAL8 dt = ts->deltaT;
      fMin = ts->f0;
      fSamp = 1.0 / dt;
      fBand = 0.5 * fSamp;
      XLAL_CHECK ( ( dataParams->fMin >= fMin ) && ( dataParams->fMin + dataParams->Band <= fMin + fBand ), XLAL_EINVAL, "Requested fMin=%f and fBand=%f are not covered by what the input timeseries can provide (fMin=%f, fBand=%f).", dataParams->fMin, dataParams->Band, fMin, fBand );
    }

  const LIGOTimeGPSVector *timestamps = dataParams->multiTimestamps->data[detectorIndex];
  const LALDetector *site = &dataParams->multiIFO.sites[detectorIndex];
  REAL8 Tsft = timestamps->deltaT;

  // if SFT output requested: need *effective* fMin and Band consistent with SFT bins
  // Note: this band is only used for internal data operations; ultimately SFTs covering
  // the half-open interval dataParams->[fMin,fMin+Band) are returned to the user using
  // XLALExtractStrictBandFromSFTVector()
  if ( SFTvect != NULL )
    {
      UINT4 firstBinEff, numBinsEff;
      XLAL_CHECK ( XLALFindCoveringSFTBins ( &firstBinEff, &numBinsEff, fMin, fBand, Tsft ) == XLAL_SUCCESS, XLAL_EFUNC );

      REAL8 fBand_eff = (numBinsEff - 1.0) / Tsft;
      REAL8 fMin_eff  = firstBinEff / Tsft;
      REAL8 fMax = fMin + dataParams->Band;
      REAL8 fMax_eff = fMin_eff + fBand_eff;
      if ( (fMin_eff != fMin) || (fBand_eff != fBand ) ) {
        XLALPrintWarning("Caller asked for Band [%.16g, %.16g] Hz, effective SFT-Band produced is [%.16g, %.16g] Hz\n",
                         fMin, fMax, fMin_eff, fMax_eff );
        XLAL_CHECK ( dataParams->inputMultiTS == NULL, XLAL_EINVAL, "Cannot expand effective frequency band with input timeseries given. Timeseries seems inconsistent with SFTs\n");
        fMin = fMin_eff;		// (potentially) lower minimal frequency to fit SFT bins
        fBand = fBand_eff;
        fSamp = 2.0 * fBand_eff;	// (potentially) higher sampling rate required to fit SFT bins
      } // if (fMin_eff != fMin) || (fBand_eff != fBand)
    } // if SFT-output requested

  // characterize the output time-series
  UINT4 n0_fSamp = (UINT4) round ( Tsft * fSamp );

  // by construction, fSamp * Tsft = integer, but if there are gaps between SFTs,
  // then we might have to sample at higher rates in order for all SFT-boundaries to
  // fall on exact timesteps of the timeseries.
  // ie we start from fsamp0 = n0_fSamp/Tsft, and then try to find the smallest
  // n1_fSamp >= n0_fSamp, such that for fsamp1 = n1_fSamp/Tsft, for all gaps i: Dt_i * fsamp1 = int
  UINT4 n1_fSamp;
  XLAL_CHECK ( XLALFindSmallestValidSamplingRate ( &n1_fSamp, n0_fSamp, timestamps ) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( (SFTvect != NULL) && (n1_fSamp != n0_fSamp) )
    {
      REAL8 fSamp1 = n1_fSamp / Tsft;	// increased sampling rate to fit all gaps
      XLALPrintWarning ( "GAPS: Initial SFT sampling frequency fSamp0= %d/%.0f = %g had to be increased to fSamp1 = %d/%.0f = %g\n",
                         n0_fSamp, Tsft, fSamp, n1_fSamp, Tsft, fSamp1 );
      XLAL_CHECK ( dataParams->inputMultiTS == NULL, XLAL_EINVAL, "Cannot expand effective frequency band with input timeseries given. Timeseries seems inconsistent with SFT timestamps\n");
      fSamp = fSamp1;
    } // if higher effective sampling rate required

  // ----- start-time and duration -----
  LIGOTimeGPS XLAL_INIT_DECL(firstGPS);
  LIGOTimeGPS XLAL_INIT_DECL(lastGPS);
  REAL8 duration;
  // NOTE: use time-interval of timeseries (if given) otherwise timestamps (must be subset of timeseries)
  if ( dataParams->inputMultiTS != NULL )
    {
      const REAL8TimeSeries *ts = dataParams->inputMultiTS->data[detectorIndex];
      XLAL_CHECK ( ts != NULL, XLAL_EINVAL );
      firstGPS = ts->epoch;
      duration = ts->data->length * ts->deltaT;
      lastGPS = firstGPS;
      XLALGPSAdd ( &lastGPS, duration );
    }
  else // use input timestamps
    {
      firstGPS = timestamps->data[0];
      lastGPS = timestamps->data [ timestamps->length - 1 ];
      XLALGPSAdd( &lastGPS, Tsft );
      duration = XLALGPSDiff ( &lastGPS, &firstGPS );
    }
  REAL8 firstGPS_REAL8 = XLALGPSGetREAL8 ( &firstGPS );
  REAL8 lastGPS_REAL8  = XLALGPSGetREAL8 ( &lastGPS );

  // start with an empty output time-series
  REAL4TimeSeries *Tseries_sum;
  {
    REAL8 numSteps = ceil ( fSamp * duration );
    XLAL_CHECK ( numSteps < (REAL8)LAL_UINT4_MAX, XLAL_EDOM, "Sorry, time-series of %g samples too long to fit into REAL4TimeSeries (maxLen = %g)\n", numSteps, (REAL8)LAL_UINT4_MAX );
    REAL8 dt = 1.0 / fSamp;
    REAL8 fHeterodyne = fMin;	// heterodyne signals at lower end of frequency-band
    CHAR *detPrefix = XLALGetChannelPrefix ( site->frDetector.name );
    XLAL_CHECK ( (Tseries_sum = XLALCreateREAL4TimeSeries ( detPrefix, &firstGPS, fHeterodyne, dt, &lalStrainUnit, (UINT4)numSteps )) != NULL, XLAL_EFUNC );
    memset ( Tseries_sum->data->data, 0, Tseries_sum->data->length * sizeof(Tseries_sum->data->data[0]) );
    XLALFree ( detPrefix );
  } // generate empty timeseries

  // add CW signals, if any
  UINT4 numPulsars = injectionSources ? injectionSources->length : 0;
  for ( UINT4 iInj = 0; iInj < numPulsars; iInj ++ )
    {
      // truncate any transient-CW timeseries to the actual support of the transient signal,
      // in order to make the generation more efficient, these 'partial timeseries'
      // will then be added to the full timeseries
      const PulsarParams *pulsarParams = &( injectionSources->data[iInj] );
      UINT4 t0, t1;
      XLAL_CHECK ( XLALGetTransientWindowTimespan ( &t0, &t1, pulsarParams->Transient ) == XLAL_SUCCESS, XLAL_EFUNC );

      // use latest possible start-time: max(t0,firstGPS), but not later than than lastGPS
      LIGOTimeGPS XLAL_INIT_DECL(signalStartGPS);
      if ( t0 <= firstGPS_REAL8 ) {
        signalStartGPS = firstGPS;
      } else if ( t0 >= lastGPS_REAL8 ) {
        signalStartGPS = lastGPS;
      } else { // firstGPS < t0 < lastGPS:
        // make sure signal start-time is an integer multiple of deltaT from timeseries start
        // to allow safe adding of resulting signal timeseries
        REAL8 offs0_aligned = round ( (t0 - firstGPS_REAL8) * fSamp ) / fSamp;
        signalStartGPS = firstGPS;
        XLALGPSAdd( &signalStartGPS, offs0_aligned );
      }

      // use earliest possible end-time: min(t1,lastGPS), but not earlier than firstGPS
      LIGOTimeGPS XLAL_INIT_DECL(signalEndGPS);
      if ( t1 >= lastGPS_REAL8 ) {
        signalEndGPS = lastGPS;
      } else if ( t1 <= firstGPS_REAL8 ) {
        signalEndGPS = firstGPS;
      } else {
        signalEndGPS.gpsSeconds = t1;
      }

      REAL8 fCoverMin, fCoverMax;
      const PulsarSpins *fkdot = &(pulsarParams->Doppler.fkdot);
      PulsarSpinRange XLAL_INIT_DECL ( spinRange );
      spinRange.refTime = pulsarParams->Doppler.refTime;
      memcpy ( spinRange.fkdot, fkdot, sizeof(spinRange.fkdot) );
      XLAL_INIT_MEM ( spinRange.fkdotBand );

      XLAL_CHECK ( XLALCWSignalCoveringBand ( &fCoverMin, &fCoverMax, &signalStartGPS, &signalEndGPS, &spinRange, pulsarParams->Doppler.asini, pulsarParams->Doppler.period, pulsarParams->Doppler.ecc ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK ( (fCoverMin >= fMin) && (fCoverMax < fMin + fBand), XLAL_EINVAL, "Error: injection signal %d:'%s' needs frequency band [%f,%f]Hz, injecting into [%f,%f]Hz\n",
                   iInj, pulsarParams->name, fCoverMin, fCoverMax, fMin, fMin + fBand );

      REAL8 signalDuration = XLALGPSDiff ( &signalEndGPS, &signalStartGPS );
      XLAL_CHECK ( signalDuration >= 0, XLAL_EFAILED, "Something went wrong, got negative signal duration = %g\n", signalDuration );
      if ( signalDuration > 0 )	// only need to do sth if transient-window had finite overlap with output TS
        {
          REAL4TimeSeries *Tseries_i = NULL;
          XLAL_CHECK ( (Tseries_i = XLALGenerateCWSignalTS ( pulsarParams, site, signalStartGPS, signalDuration, fSamp, fMin, edat, dataParams->sourceDeltaT )) != NULL, XLAL_EFUNC );

          // since XLALAddREAL4TimeSeries() does not enforce strict sample alignment,
          // we do our own safety check here
          REAL8 Delta_epoch = XLALGPSDiff(&Tseries_sum->epoch, &Tseries_i->epoch);
          REAL8 bin_mismatch = fabs(Delta_epoch / Tseries_sum->deltaT);
          REAL8 mismatch = fabs(bin_mismatch - round(bin_mismatch))*Tseries_sum->deltaT;
          XLAL_CHECK ( mismatch <= 1e-9, XLAL_EDATA, "Incompatible start-times when adding signal time series %d of %d, bins misaligned by %g seconds (%g bins).\n", iInj, numPulsars, mismatch, bin_mismatch );

          XLAL_CHECK ( (Tseries_sum = XLALAddREAL4TimeSeries ( Tseries_sum, Tseries_i )) != NULL, XLAL_EFUNC );
          XLALDestroyREAL4TimeSeries ( Tseries_i );
        }
    } // for iInj < numSources

  /* add Gaussian noise if requested */
  REAL8 sqrtSn = dataParams->multiNoiseFloor.sqrtSn[detectorIndex];
  if ( sqrtSn > 0)
    {
      REAL8 noiseSigma = sqrtSn * sqrt ( 0.5 * fSamp );
      INT4 randSeed = (dataParams->randSeed == 0) ? 0 : (dataParams->randSeed + detectorIndex);	// seed=0 means to use /dev/urandom, so don't touch it
      XLAL_CHECK ( XLALAddGaussianNoise ( Tseries_sum, noiseSigma, randSeed ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

  // convert final signal+Gaussian-noise timeseries into REAL8 precision:
  REAL8TimeSeries *outTS;
  XLAL_CHECK ( (outTS = XLALConvertREAL4TimeSeriesToREAL8 ( Tseries_sum )) != NULL, XLAL_EFUNC );
  XLALDestroyREAL4TimeSeries ( Tseries_sum );

  // add input noise time-series here if given
  if ( dataParams->inputMultiTS != NULL ) {
    // since XLALAddREAL8TimeSeries() does not enforce strict sample alignment,
    // we do our own safety check here
    REAL8 Delta_epoch = XLALGPSDiff(&outTS->epoch, &dataParams->inputMultiTS->data[detectorIndex]->epoch);;
    REAL8 bin_mismatch = fabs(Delta_epoch / outTS->deltaT);
    REAL8 mismatch = fabs(bin_mismatch - round(bin_mismatch))*outTS->deltaT;
    XLAL_CHECK ( mismatch <= 1e-9, XLAL_EDATA, "Incompatible start-times when adding input noise time-series, bins misaligned by %g seconds (%g bins).\n", mismatch, bin_mismatch );
    XLAL_CHECK ( (outTS = XLALAddREAL8TimeSeries ( outTS, dataParams->inputMultiTS->data[detectorIndex] )) != NULL, XLAL_EFUNC );
  }

  // turn final timeseries into SFTs, if requested
  if ( SFTvect != NULL )
    {
      // compute SFTs from timeseries
      SFTVector *sftVect;
      XLAL_CHECK ( (sftVect = XLALMakeSFTsFromREAL8TimeSeries ( outTS, timestamps, dataParams->SFTWindowType, dataParams->SFTWindowParam)) != NULL, XLAL_EFUNC );

      // extract requested band
      XLAL_CHECK ( ((*SFTvect) = XLALExtractStrictBandFromSFTVector ( sftVect, dataParams->fMin, dataParams->Band )) != NULL, XLAL_EFUNC );
      XLALDestroySFTVector ( sftVect );
    } // if SFTvect

  // return timeseries if requested
  if ( Tseries != NULL ) {
    (*Tseries) = outTS;
  } else {
    XLALDestroyREAL8TimeSeries ( outTS );
  }

  return XLAL_SUCCESS;

} // XLALCWMakeFakeData()


/**
 * Generate a (heterodyned) REAL4 timeseries of a CW signal for given pulsarParams,
 * site, start-time, duration, and sampling-rate
 *
 * NOTE: this is mostly an API-wrapper to the more 'old-style' function
 * XLALGeneratePulsarSignal() [which will become deprecated in the future],
 * extended for the option to generate transient-CW signals
 */
REAL4TimeSeries *
XLALGenerateCWSignalTS ( const PulsarParams *pulsarParams,	///< input CW pulsar-signal parameters
                         const LALDetector *site,		///< detector
                         LIGOTimeGPS startTime,			///< time-series start-time GPS
                         REAL8 duration,			///< time-series duration to generate
                         REAL8 fSamp,				///< sampling frequency
                         REAL8 fHet,				///< heterodyning frequency
                         const EphemerisData *edat,		///< ephemeris data
                         REAL8 sourceDeltaT                     ///< source-frame sampling period (optional: 0 == previous internal defaults)
                         )
{
  XLAL_CHECK_NULL ( pulsarParams != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( site != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( duration > 0, XLAL_EDOM );
  XLAL_CHECK_NULL ( fSamp > 0, XLAL_EDOM );
  XLAL_CHECK_NULL ( fHet >= 0, XLAL_EDOM );

  // translate amplitude params
  REAL8 aPlus  = pulsarParams->Amp.aPlus;
  REAL8 aCross = pulsarParams->Amp.aCross;
  // translate 'modern' fkdot into 'old-style' spindown-vector
  UINT4 s_max;
  for ( s_max = PULSAR_MAX_SPINS-1; s_max > 0; s_max -- )
    {
      if ( pulsarParams->Doppler.fkdot[s_max] != 0 )
        break;
    } // for s_max = max ... 0
  REAL8Vector *spindown = NULL;
  if ( s_max > 0 )
    {
      XLAL_CHECK_NULL ( (spindown = XLALCreateREAL8Vector ( s_max )) != NULL, XLAL_EFUNC );
      for ( UINT4 s = 0; s < s_max; s ++ ) {
        spindown->data[s] = pulsarParams->Doppler.fkdot[s+1];
      }
    }

  /*----------------------------------------
   * fill old-style PulsarSignalParams struct
   *----------------------------------------*/
  PulsarSignalParams XLAL_INIT_DECL(params);
  params.pulsar.refTime            = pulsarParams->Doppler.refTime;
  params.pulsar.position.system    = COORDINATESYSTEM_EQUATORIAL;
  params.pulsar.position.longitude = pulsarParams->Doppler.Alpha;
  params.pulsar.position.latitude  = pulsarParams->Doppler.Delta;
  params.pulsar.aPlus              = aPlus;
  params.pulsar.aCross             = aCross;
  params.pulsar.phi0               = pulsarParams->Amp.phi0;
  params.pulsar.psi                = pulsarParams->Amp.psi;
  params.pulsar.f0                 = pulsarParams->Doppler.fkdot[0];
  params.pulsar.spindown           = spindown;
  params.orbit.tp                  = pulsarParams->Doppler.tp;
  params.orbit.argp                = pulsarParams->Doppler.argp;
  params.orbit.asini               = pulsarParams->Doppler.asini;
  params.orbit.ecc                 = pulsarParams->Doppler.ecc;
  params.orbit.period              = pulsarParams->Doppler.period;
  params.transfer                  = NULL;
  params.ephemerides               = edat;
  params.fHeterodyne               = fHet;
  params.sourceDeltaT              = sourceDeltaT;

  // detector-specific settings
  params.startTimeGPS              = startTime;
  params.duration                  = ceil ( duration );
  params.samplingRate              = fSamp;
  params.site                      = site;

  /*----------------------------------------
   * generate the signal time-series
   *----------------------------------------*/
  REAL4TimeSeries *Tseries;
  XLAL_CHECK_NULL ( (Tseries = XLALGeneratePulsarSignal ( &params )) != NULL, XLAL_EFUNC );
  // ----- free internal memory
  XLALDestroyREAL8Vector ( spindown );

  // ----- apply transient-CW window
  XLAL_CHECK_NULL ( XLALApplyTransientWindow ( Tseries, pulsarParams->Transient ) == XLAL_SUCCESS, XLAL_EFUNC );

  return Tseries;

} // XLALGenerateCWSignalTS()


///
/// Make SFTs from given REAL8TimeSeries at given timestamps, potentially applying a time-domain window on each timestretch first
///
SFTVector *
XLALMakeSFTsFromREAL8TimeSeries ( const REAL8TimeSeries *timeseries,	//!< input time-series
                                  const LIGOTimeGPSVector *timestamps, 	//!< timestamps to produce SFTs for (can be NULL), if given must all lies within timeseries' time-span
                                  const char *windowType,		//!< optional time-domain window function to apply before FFTing
                                  REAL8 windowParam			//!< window parameter, if any
                                  )
{
  XLAL_CHECK_NULL ( timeseries != NULL, XLAL_EINVAL, "Invalid NULL input 'timeseries'\n");
  XLAL_CHECK_NULL ( timestamps != NULL, XLAL_EINVAL, "Invalid NULL input 'timestamps'\n");

  REAL8 dt = timeseries->deltaT;	// timeseries timestep */
  REAL8 Tsft = timestamps->deltaT;
  REAL8 df = 1.0 / Tsft;		// SFT frequency spacing

  // make sure that number of timesamples/SFT is an integer (up to possible rounding error 'eps')
  REAL8 timestepsSFT0 = Tsft / dt;
  UINT4 timestepsSFT  = lround ( timestepsSFT0 );
  XLAL_CHECK_NULL ( fabs ( timestepsSFT0 - timestepsSFT ) / timestepsSFT0 < eps, XLAL_ETOL,
                    "Inconsistent sampling-step (dt=%g) and Tsft=%g: must be integer multiple Tsft/dt = %g >= %g\n",
                    dt, Tsft, timestepsSFT0, eps );

  // prepare window function if requested
  REAL8Window *window = NULL;
  if ( windowType != NULL ) {
    XLAL_CHECK_NULL ( (window = XLALCreateNamedREAL8Window ( windowType, windowParam, timestepsSFT )) != NULL, XLAL_EFUNC );
  }

  // ---------- Prepare FFT ----------
  REAL8Vector *timeStretchCopy;	// input array of length N
  XLAL_CHECK_NULL ( (timeStretchCopy = XLALCreateREAL8Vector ( timestepsSFT )) != NULL, XLAL_EFUNC, "XLALCreateREAL4Vector(%d) failed.\n", timestepsSFT );
  UINT4 numSFTBins = timestepsSFT / 2 + 1;	// number of positive frequency-bins + 'DC' to be stored in SFT
  fftw_complex *fftOut;	// output array of length N/2 + 1
  XLAL_CHECK_NULL ( (fftOut = fftw_malloc ( numSFTBins * sizeof(fftOut[0]) )) != NULL, XLAL_ENOMEM, "fftw_malloc(%d*sizeof(complex)) failed\n", numSFTBins );
  fftw_plan fftplan;	// FFTW plan
  LAL_FFTW_WISDOM_LOCK;
  XLAL_CHECK_NULL ( (fftplan = fftw_plan_dft_r2c_1d ( timestepsSFT, timeStretchCopy->data, fftOut, FFTW_ESTIMATE)) != NULL, XLAL_EFUNC );	// FIXME: or try FFTW_MEASURE
  LAL_FFTW_WISDOM_UNLOCK;

  LIGOTimeGPS tStart = timeseries->epoch;

  // get last possible start-time for an SFT
  REAL8 duration =  round ( timeseries->data->length * dt ); // rounded to seconds
  LIGOTimeGPS tLast = tStart;
  XLALGPSAdd( &tLast, duration - Tsft );

  // check that all timestamps lie within [tStart, tLast]
  for ( UINT4 i = 0; i < timestamps->length; i ++ )
    {
      char buf1[256], buf2[256];
      XLAL_CHECK_NULL ( XLALGPSDiff ( &tStart, &(timestamps->data[i]) ) <= 0, XLAL_EDOM, "Timestamp i=%d: %s before start-time %s\n",
                        i, XLALGPSToStr ( buf1, &(timestamps->data[i]) ), XLALGPSToStr ( buf2, &tStart ) );
      XLAL_CHECK_NULL ( XLALGPSDiff ( &tLast,   &(timestamps->data[i]) ) >=0, XLAL_EDOM, "Timestamp i=%d: %s after last start-time %s\n",
                        i, XLALGPSToStr ( buf1, &(timestamps->data[i]) ), XLALGPSToStr ( buf2, &tLast ) );
    }

  UINT4 numSFTs = timestamps->length;

  // prepare output SFT-vector
  SFTVector *sftvect;
  XLAL_CHECK_NULL ( (sftvect = XLALCreateSFTVector ( numSFTs, numSFTBins )) != NULL, XLAL_EFUNC,
                    "XLALCreateSFTVector(numSFTs=%d, numBins=%d) failed.\n", numSFTs, numSFTBins );

  // main loop: apply FFT to the requested time-stretches and store in output SFTs
  for ( UINT4 iSFT = 0; iSFT < numSFTs; iSFT++ )
    {
      SFTtype *thisSFT = &(sftvect->data[iSFT]);	// point to current SFT-slot to store output in

      // find the start-bin for this SFT in the time-series
      REAL8 offset = XLALGPSDiff ( &(timestamps->data[iSFT]), &tStart );
      INT4 offsetBins = lround ( offset / dt );

      // copy timeseries-data for that SFT into local buffer
      memcpy ( timeStretchCopy->data, timeseries->data->data + offsetBins, timeStretchCopy->length * sizeof(timeStretchCopy->data[0]) );

      // window the current time series stretch if required
      REAL8 sigma_window = 1;
      if ( window != NULL )
        {
	  sigma_window = sqrt ( window->sumofsquares / window->data->length );
	  for( UINT4 iBin = 0; iBin < timeStretchCopy->length; iBin++ ) {
            timeStretchCopy->data[iBin] *= window->data->data[iBin];
          }
        } // if window

      // FFT this time-stretch
      fftw_execute ( fftplan );

      // fill the header of the i'th output SFT */
      strcpy ( thisSFT->name, timeseries->name );
      thisSFT->epoch = timestamps->data[iSFT];
      thisSFT->f0 = timeseries->f0;			// SFT starts at heterodyning frequency
      thisSFT->deltaF = df;

      // normalize DFT-data to conform to SFT specification ==> multiply DFT by (dt/sigma{window})
      // the SFT normalization in case of windowing follows the conventions detailed in \cite SFT-spec
      REAL8 norm = dt / sigma_window;
      for ( UINT4 k = 0; k < numSFTBins ; k ++ ) {
        thisSFT->data->data[k] = (COMPLEX8) ( norm * fftOut[k] );
      }

      // correct heterodyning-phase, IF NECESSARY: ie if (fHet * tStart) is not an integer, such that phase-corr = multiple of 2pi
      if ( ( (INT4)timeseries->f0 != timeseries->f0  ) || (timeseries->epoch.gpsNanoSeconds != 0) || (thisSFT->epoch.gpsNanoSeconds != 0) ) {
        XLAL_CHECK_NULL ( XLALcorrect_phase ( thisSFT, timeseries->epoch) == XLAL_SUCCESS, XLAL_EFUNC );
      }

    } // for iSFT < numSFTs

  // free memory
  fftw_free ( fftOut );
  LAL_FFTW_WISDOM_LOCK;
  fftw_destroy_plan ( fftplan );
  LAL_FFTW_WISDOM_UNLOCK;
  XLALDestroyREAL8Vector ( timeStretchCopy );
  XLALDestroyREAL8Window ( window );

  return sftvect;

} // XLALMakeSFTsFromREAL8TimeSeries()



/**
 * Find the smallest sampling rate of the form fsamp = n / Tsft, with n>=n0,
 * such that all gap sizes Dt_i between SFTs of the given timestamps are also
 * exactly resolved, ie. that Dt_i * fsamp = integer, for all i
 *
 * The smallest allowed sampling rate is the user-specified fsamp0 = n0 / Tsft,
 * which guarantees by construction that fSamp0 * Tsft = n0 = integer
 * This sampling rate would be valid if there are no gaps between SFTs,
 * so it's only in cases of gaps that are non-integer multiples of Tsft that we'll
 * (potentially) have to increase the sampling rate.
 *
 * NOTE: This approach replaces the old mfdv4 habit of 'nudging' noise SFTs start-times
 * to fall on integer timesteps of the fsamp0 timeseries. The purpose of this function
 * is to avoid this behaviour, by appropriately increasing the sampling rate
 * as required.
 *
 * NOTE2: we only allow integer-second gaps, everything else will be
 * rejected with an error-message.
 *
 * NOTE3: Tsft=timestamps->deltaT must be integer seconds, everything else is rejected
 * with an error as well
 */
int
XLALFindSmallestValidSamplingRate ( UINT4 *n1,				//< [out] minimal valid sampling rate n1/Tsft
                                    UINT4 n0, 				//< [in] minimal sampling rate n0/Tsft
                                    const LIGOTimeGPSVector *timestamps	//< [in] start-timestamps and length Tsft of SFTs
                                    )
{
  XLAL_CHECK ( n1 != NULL, XLAL_EINVAL );
  XLAL_CHECK ( n0 > 0, XLAL_EINVAL );
  XLAL_CHECK ( timestamps && (timestamps->length > 0), XLAL_EINVAL );
  REAL8 TsftREAL = timestamps->deltaT;
  XLAL_CHECK ( TsftREAL == round(TsftREAL), XLAL_EDOM, "Only exact integer-second Tsft allowed, got Tsft = %.16g s\n", TsftREAL );
  UINT4 Tsft = (UINT4)TsftREAL;
  XLAL_CHECK ( Tsft > 0, XLAL_EINVAL );

  // NOTE: all 'valid' sampling rates are of the form  fSamp = n / Tsft, where n >= n0,
  // therefore we label sampling rates by their index 'n'. The purpose is to find the
  // smallest 'valid' n, which is the max of the required 'n' over all gaps
  UINT4 nCur = n0;

  // ----- We just step through the vector and figure out for each gap if we need to
  // decrease the stepsize Tsft/n from the previous value
  UINT4 numSFTs = timestamps->length;

  for ( UINT4 i = 1; i < numSFTs; i ++ )
    {
      LIGOTimeGPS *t1 = &(timestamps->data[i]);
      LIGOTimeGPS *t0 = &(timestamps->data[i-1]);

      INT4 nsdiff = t1->gpsNanoSeconds - t0->gpsNanoSeconds;
      XLAL_CHECK ( nsdiff == 0, XLAL_EDOM, "Only integer-second gaps allowed, found %d ns excess in gap between i=%d and i-1\n", nsdiff, i );

      INT4 gap_i0 = t1->gpsSeconds - t0->gpsSeconds;
      XLAL_CHECK ( gap_i0 > 0, XLAL_EDOM, "Timestamps must be sorted in increasing order, found negative gap %d s between i=%d and i-1\n", gap_i0, i );

      // now reduce gap to remainder wrt Tsft
      INT4 gap_i = gap_i0 % Tsft;

      if ( ((INT8)gap_i * nCur) % Tsft == 0 ) {
        continue;
      }

      // otherwise:
      // solve for required new smaller step-size 'dt = Tsft/nNew' to fit integer cycles
      // both into Tsft (by construction) and into (gap_i mod Tsft), with nNew > nCur
      //
      // gap[i] == (t[i+1] - t[i]) mod Tsft
      // such that 0 <= gap[i] < Tsft
      // we want integers nNew, m , such that nNew > nCur, and m < nNew, with
      // nNew * dtNew = Tsft   AND  m * dtNew = gap[i]
      // ==> m / nNew = gap[i] / Tsft
      // This could be solved easily by rounding nCur to next highest
      // multiple of Tsft: nNew' = ceil(nCur/Tsft) * Tsft > nCur
      // but this can be wasteful if the fraction simplifies: so we compute
      // the greatest common divisor 'g'
      // g = gcd(gap[i], Tsft), and then use
      // nNew = ceil ( nCur  * g  / Tsft ) * Tsft / g > nCur

      UINT4 g = gcd ( gap_i, Tsft );
      REAL8 Tg = TsftREAL / g;
      UINT4 nNew = (UINT4) ceil ( nCur / Tg ) * Tg;

      XLAL_CHECK ( nNew > nCur, XLAL_ETOL, "This seems wrong: nNew = %d !> nCur = %d, but should be greater!\n", nNew, nCur );
      XLALPrintInfo ("Need to increase from fSamp = %d/Tsft = %g to fSamp = %d / Tsft = %g\n", nCur, nCur/TsftREAL, nNew, nNew / TsftREAL );

      nCur = nNew;

    } // for i < numSFTs


  // our final minimal valid sampling rate is therefore n1/Tsft
  (*n1) = nCur;

  return XLAL_SUCCESS;

} // XLALFindSmallestValidSamplingRate()


/* Find greatest common divisor between two numbers,
 * where numer <= denom
 * this is an implementation of the Euclidean Algorithm,
 * taken from John's UnitNormalize.c, and extended to UINT4's
 *
 * For reference, see also
 * https://en.wikipedia.org/wiki/Euclidean_algorithm#Implementations
 */
static UINT4
gcd (UINT4 numer, UINT4 denom)
{
  UINT4 next_numer, next_denom, rmdr;

  next_numer = numer;
  next_denom = denom;
  while ( next_denom != 0 )
    {
      rmdr = next_numer % next_denom;
      next_numer = next_denom;
      next_denom = rmdr;
    }
  return next_numer;
} // gcd

/**
 * Create *zero-initialized* PulsarParamsVector for numPulsars
 */
PulsarParamsVector *
XLALCreatePulsarParamsVector ( UINT4 numPulsars )
{
  PulsarParamsVector *ret;
  XLAL_CHECK_NULL ( ( ret = XLALCalloc ( 1, sizeof(*ret))) != NULL, XLAL_ENOMEM );

  ret->length = numPulsars;
  if ( numPulsars > 0 ) {
    XLAL_CHECK_NULL ( (ret->data = XLALCalloc ( numPulsars, sizeof(ret->data[0]))) != NULL, XLAL_ENOMEM );
  }

  return ret;

} // XLALCreatePulsarParamsVector()

/**
 * Destructor for PulsarParamsVector type
 */
void
XLALDestroyPulsarParamsVector ( PulsarParamsVector *ppvect )
{
  if ( ppvect == NULL ) {
    return;
  }

  if ( ppvect->data != NULL )
    {
      XLALFree ( ppvect->data );
    }

  XLALFree ( ppvect );

  return;
} // XLALDestroyPulsarParamsVector()


/**
 * Function to parse a config-file-type string (or section thereof)
 * into a PulsarParams struct.
 *
 * NOTE: The section-name is optional, and can be given as NULL,
 * in which case the top of the file (ie the "default section")
 * is used.
 *
 * NOTE2: eventually ATNF/TEMPO2-style 'par-file' variables will also
 * be understood by this function, but we start out with a simpler version
 * that just deals with our 'CW-style' input variable for now
 *
 * NOTE3: The config-file must be of a special "SourceParamsIO" form,
 * defining the following required and optional parameters:
 *
 * REQUIRED:
 * Alpha, Delta, Freq, refTime (unless refTimeDef != NULL)
 *
 * OPTIONAL:
 * f1dot, f2dot, f3dot, f4dot, f5dot, f6dot
 * {h0, cosi} or {aPlus, aCross}, psi, phi0
 * transientWindowType, transientStartTime, transientTau
 *
 * Other config-variables found in the file will ... ?? error or accept?
 */
int
XLALReadPulsarParams ( PulsarParams *pulsarParams,	///< [out] pulsar parameters to fill in from config string
                       LALParsedDataFile *cfgdata,      ///< [in] pre-parsed "SourceParamsIO" config-file contents
                       const CHAR *secName,		///< [in] section-name to use from config-file string (can be NULL)
                       const LIGOTimeGPS *refTimeDef	///< [in] default reference time if refTime is not given
                       )
{
  XLAL_CHECK ( pulsarParams != NULL, XLAL_EINVAL );
  XLAL_CHECK ( cfgdata != NULL, XLAL_EINVAL );

  XLAL_INIT_MEM ( (*pulsarParams) );	// wipe input struct clean

  // ---------- PulsarAmplitudeParams ----------
  // ----- h0, cosi
  REAL8 h0 = 0; BOOLEAN have_h0;
  REAL8 cosi = 0; BOOLEAN have_cosi;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &h0, cfgdata, secName, "h0", &have_h0 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &cosi, cfgdata, secName, "cosi", &have_cosi ) == XLAL_SUCCESS, XLAL_EFUNC );
  // ----- ALTERNATIVE: aPlus, aCross
  REAL8 aPlus = 0; BOOLEAN have_aPlus;
  REAL8 aCross = 0; BOOLEAN have_aCross;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &aPlus, cfgdata, secName, "aPlus", &have_aPlus ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &aCross, cfgdata, secName, "aCross", &have_aCross ) == XLAL_SUCCESS, XLAL_EFUNC );

  // if h0 then also need cosi, and vice-versa
  XLAL_CHECK ( (have_h0 && have_cosi) || ( !have_h0 && !have_cosi ), XLAL_EINVAL );
  // if aPlus then also need aCross, and vice-versa
  XLAL_CHECK ( (have_aPlus && have_aCross) || ( !have_aPlus && !have_aCross ), XLAL_EINVAL );
  // {h0,cosi} or {aPlus, aCross} mutually exclusive sets
  XLAL_CHECK ( ! ( have_h0 && have_aPlus ), XLAL_EINVAL );

  if ( have_h0 )	/* translate {h_0, cosi} into A_{+,x}*/
    {
        XLAL_CHECK ( h0 >= 0, XLAL_EDOM );
        XLAL_CHECK ( (cosi >= -1) && (cosi <= 1), XLAL_EDOM );
        aPlus = 0.5 * h0 * (1.0 + SQ(cosi));
        aCross = h0 * cosi;
    } //if {h0, cosi}
  pulsarParams->Amp.aPlus = aPlus;
  pulsarParams->Amp.aCross = aCross;

  // ----- psi
  REAL8 psi = 0; BOOLEAN have_psi;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &psi, cfgdata, secName, "psi", &have_psi ) == XLAL_SUCCESS, XLAL_EFUNC );
  pulsarParams->Amp.psi = psi;

  // ----- phi0
  REAL8 phi0 = 0; BOOLEAN have_phi0;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &phi0, cfgdata, secName, "phi0", &have_phi0 ) == XLAL_SUCCESS, XLAL_EFUNC );
  pulsarParams->Amp.phi0 = phi0;

  // ---------- PulsarDopplerParams ----------

  // ----- refTime
  LIGOTimeGPS refTime_GPS; BOOLEAN have_refTime;
  XLAL_CHECK ( XLALReadConfigEPOCHVariable ( &refTime_GPS, cfgdata, secName, "refTime", &have_refTime ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( have_refTime ) {
    pulsarParams->Doppler.refTime = refTime_GPS;
  } else if ( refTimeDef != NULL ) {
    pulsarParams->Doppler.refTime = *refTimeDef;
  } else {
    XLAL_ERROR ( XLAL_EINVAL, "missing value refTime, and no default is available" );
  }

  // ----- Alpha
  REAL8 Alpha_Rad = 0; BOOLEAN have_Alpha;
  XLAL_CHECK ( XLALReadConfigRAJVariable ( &Alpha_Rad, cfgdata, secName, "Alpha", &have_Alpha ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( have_Alpha, XLAL_EINVAL, "missing required value Alpha" );

  XLAL_CHECK ( (Alpha_Rad >= 0) && (Alpha_Rad < LAL_TWOPI), XLAL_EDOM );
  pulsarParams->Doppler.Alpha = Alpha_Rad;

  // ----- Delta
  REAL8 Delta_Rad = 0; BOOLEAN have_Delta;
  XLAL_CHECK ( XLALReadConfigDECJVariable ( &Delta_Rad, cfgdata, secName, "Delta", &have_Delta ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( have_Delta, XLAL_EINVAL, "missing required value Delta" );

  XLAL_CHECK ( (Delta_Rad >= -LAL_PI_2) && (Delta_Rad <= LAL_PI_2), XLAL_EDOM );
  pulsarParams->Doppler.Delta = Delta_Rad;

  // ----- fkdot
  // Freq
  REAL8 Freq = 0; BOOLEAN have_Freq;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &Freq, cfgdata, secName, "Freq", &have_Freq ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( have_Freq, XLAL_EINVAL, "missing required value Freq" );

  XLAL_CHECK ( Freq > 0, XLAL_EDOM );
  pulsarParams->Doppler.fkdot[0] = Freq;

  // f1dot
  REAL8 f1dot = 0; BOOLEAN have_f1dot;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &f1dot, cfgdata, secName, "f1dot", &have_f1dot ) == XLAL_SUCCESS, XLAL_EFUNC );
  pulsarParams->Doppler.fkdot[1] = f1dot;
  // f2dot
  REAL8 f2dot = 0; BOOLEAN have_f2dot;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &f2dot, cfgdata, secName, "f2dot", &have_f2dot ) == XLAL_SUCCESS, XLAL_EFUNC );
  pulsarParams->Doppler.fkdot[2] = f2dot;
  // f3dot
  REAL8 f3dot = 0; BOOLEAN have_f3dot;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &f3dot, cfgdata, secName, "f3dot", &have_f3dot ) == XLAL_SUCCESS, XLAL_EFUNC );
  pulsarParams->Doppler.fkdot[3] = f3dot;
  // f4dot
  REAL8 f4dot = 0; BOOLEAN have_f4dot;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &f4dot, cfgdata, secName, "f4dot", &have_f4dot ) == XLAL_SUCCESS, XLAL_EFUNC );
  pulsarParams->Doppler.fkdot[4] = f4dot;
  // f5dot
  REAL8 f5dot = 0; BOOLEAN have_f5dot;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &f5dot, cfgdata, secName, "f5dot", &have_f5dot ) == XLAL_SUCCESS, XLAL_EFUNC );
  pulsarParams->Doppler.fkdot[5] = f5dot;
  // f6dot
  REAL8 f6dot = 0; BOOLEAN have_f6dot;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &f6dot, cfgdata, secName, "f6dot", &have_f6dot ) == XLAL_SUCCESS, XLAL_EFUNC );
  pulsarParams->Doppler.fkdot[6] = f6dot;

  // ----- orbit
  LIGOTimeGPS orbitTp; BOOLEAN have_orbitTp;
  XLAL_CHECK ( XLALReadConfigEPOCHVariable ( &orbitTp, cfgdata, secName, "orbitTp", &have_orbitTp ) == XLAL_SUCCESS, XLAL_EFUNC );
  REAL8 orbitArgp = 0; 	BOOLEAN have_orbitArgp;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &orbitArgp, cfgdata, secName, "orbitArgp", &have_orbitArgp ) == XLAL_SUCCESS, XLAL_EFUNC );
  REAL8 orbitasini = 0 /* isolated pulsar */; BOOLEAN have_orbitasini;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &orbitasini, cfgdata, secName, "orbitasini", &have_orbitasini ) == XLAL_SUCCESS, XLAL_EFUNC );
  REAL8 orbitEcc = 0;  	BOOLEAN have_orbitEcc;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &orbitEcc, cfgdata, secName, "orbitEcc", &have_orbitEcc ) == XLAL_SUCCESS, XLAL_EFUNC );
  REAL8 orbitPeriod = 0;BOOLEAN have_orbitPeriod;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &orbitPeriod, cfgdata, secName, "orbitPeriod", &have_orbitPeriod ) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( have_orbitasini || have_orbitEcc || have_orbitPeriod || have_orbitArgp || have_orbitTp )
    {
      XLAL_CHECK ( orbitasini >= 0, XLAL_EDOM );
      XLAL_CHECK ( (orbitasini == 0) || ( have_orbitPeriod && (orbitPeriod > 0) && have_orbitTp ), XLAL_EINVAL, "If orbitasini>0 then we also need 'orbitPeriod>0' and 'orbitTp'\n" );
      XLAL_CHECK ( (orbitEcc >= 0) && (orbitEcc <= 1), XLAL_EDOM );

      /* fill in orbital parameter structure */
      pulsarParams->Doppler.tp 		= orbitTp;
      pulsarParams->Doppler.argp 	= orbitArgp;
      pulsarParams->Doppler.asini 	= orbitasini;
      pulsarParams->Doppler.ecc 	= orbitEcc;
      pulsarParams->Doppler.period 	= orbitPeriod;

    } // if have non-trivial orbit

  // ---------- transientWindow_t ----------

  // ----- type
  char *transientWindowType = NULL; BOOLEAN have_transientWindowType;
  XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &transientWindowType, cfgdata, secName, "transientWindowType", &have_transientWindowType ) == XLAL_SUCCESS, XLAL_EFUNC );
  // ----- t0
  UINT4 transientStartTime = 0; BOOLEAN have_transientStartTime;
  XLAL_CHECK ( XLALReadConfigUINT4Variable ( &transientStartTime, cfgdata, secName, "transientStartTime", &have_transientStartTime ) == XLAL_SUCCESS, XLAL_EFUNC );
  // ----- tau (still keeping deprecated Days variant for backwards compatibility, for now)
  UINT4 transientTau = 0; BOOLEAN have_transientTau;
  XLAL_CHECK ( XLALReadConfigUINT4Variable ( &transientTau, cfgdata, secName, "transientTau", &have_transientTau ) == XLAL_SUCCESS, XLAL_EFUNC );
  REAL8 transientTauDays = 0; BOOLEAN have_transientTauDays;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &transientTauDays, cfgdata, secName, "transientTauDays", &have_transientTauDays ) == XLAL_SUCCESS, XLAL_EFUNC );

  int twtype = TRANSIENT_NONE;
  if ( have_transientWindowType ) {
    XLAL_CHECK ( (twtype = XLALParseTransientWindowName ( transientWindowType )) >= 0, XLAL_EFUNC );
    XLALFree ( transientWindowType );
  }
  pulsarParams->Transient.type = twtype;

  if ( pulsarParams->Transient.type != TRANSIENT_NONE )
    {
      XLAL_CHECK ( have_transientStartTime && (have_transientTau || have_transientTauDays), XLAL_EINVAL, "For transientWindowType!=None, we also need transientStartTime and either transientTau (deprecated) or transientTau.");
      XLAL_CHECK ( !(have_transientTau && have_transientTauDays), XLAL_EINVAL, "Cannot have both transientTau and transientTauDays; the latter is deprecated." );

      pulsarParams->Transient.t0   = transientStartTime;
      if ( have_transientTauDays ) {
        XLAL_CHECK ( transientTauDays > 0, XLAL_EDOM );
        printf("Warning: Option transientTauDays is deprecated, please switch to using transientTau [UINT4, in seconds] instead.\n");
        pulsarParams->Transient.tau  = (UINT4) round ( transientTauDays * 86400 );
      }
      else {
        pulsarParams->Transient.tau  = transientTau;
      }
    } /* if transient window != none */
  else
    {
      XLAL_CHECK ( !(have_transientStartTime || have_transientTau || have_transientTauDays), XLAL_EINVAL, "Cannot use transientStartTime, transientTau or transientTauDays without transientWindowType!" );
    }

  return XLAL_SUCCESS;
} // XLALParsePulsarParams()


/**
 * Parse a given 'CWsources' config file for PulsarParams, return vector
 * of all pulsar definitions found [using sections]
 */
PulsarParamsVector *
XLALPulsarParamsFromFile ( const char *fname, 			///< [in] 'CWsources' config file name
                           const LIGOTimeGPS *refTimeDef	///< [in] default reference time if refTime is not given
                           )
{
  XLAL_CHECK_NULL ( fname != NULL, XLAL_EINVAL );

  LALParsedDataFile *cfgdata = NULL;
  XLAL_CHECK_NULL ( XLALParseDataFile ( &cfgdata, fname ) == XLAL_SUCCESS, XLAL_EFUNC );

  LALStringVector *sections;
  XLAL_CHECK_NULL ( (sections = XLALListConfigFileSections ( cfgdata )) != NULL, XLAL_EFUNC );

  UINT4 numPulsars = sections->length;	// currently only single-section defs supported! FIXME

  PulsarParamsVector *sources;
  XLAL_CHECK_NULL ( (sources = XLALCreatePulsarParamsVector ( numPulsars )) != NULL, XLAL_EFUNC );

  for ( UINT4 i = 0; i < numPulsars; i ++ )
    {
      const char *sec_i = sections->data[i];

      if ( strcmp ( sec_i, "default" ) == 0 ) {	// special handling of 'default' section
        sec_i = NULL;
      }
      XLAL_CHECK_NULL ( XLALReadPulsarParams ( &sources->data[i], cfgdata, sec_i, refTimeDef ) == XLAL_SUCCESS, XLAL_EFUNC );

      // ----- source naming convention: 'filename:section'
      snprintf ( sources->data[i].name, sizeof(sources->data[i].name), "%s:%s", fname, sections->data[i] );
      sources->data[i].name[sizeof(sources->data[i].name)-1] = '\0';

    } // for i < numPulsars

  XLALDestroyStringVector ( sections );

  XLAL_CHECK_NULL ( XLALCheckConfigFileWasFullyParsed( fname, cfgdata ) == XLAL_SUCCESS, XLAL_EINVAL );
  XLALDestroyParsedDataFile ( cfgdata );

  return sources;

} // XLALPulsarParamsFromFile()


// internal helper function: check that config-file was fully parsed (no leftover un-parsed lines),
// throw error if not and list unparsed entries
int
XLALCheckConfigFileWasFullyParsed ( const char *fname, const LALParsedDataFile *cfgdata )
{
  XLAL_CHECK ( cfgdata != NULL, XLAL_EINVAL );

  UINT4Vector *n_unread = XLALConfigFileGetUnreadEntries ( cfgdata );
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC );

  if ( n_unread != NULL )
    {
      XLALPrintError ( "ERROR: Pulsar params config file '%s' contained '%d' unknown entries:\n", fname, n_unread->length );
      for ( UINT4 i = 0; i < n_unread->length; i ++ )
        {
          XLALPrintError ( "%s'%s'", i > 0 ? ", " : "", cfgdata->lines->tokens[ n_unread->data[i] ] );
        }
      XLALPrintError ( "\n" );
      XLAL_ERROR ( XLAL_EINVAL );
    }

  return XLAL_SUCCESS;

} // XLALCheckConfigFileWasFullyParsed()


/**
 * Function to determine the PulsarParamsVector input from a user-input defining CW sources.
 *
 * This option supports a dual-type feature: if any string in the list is of the form '{...}', then
 * it determines the *contents* of a config-file, otherwise the name-pattern of config-files to be parsed by XLALFindFiles(),
 * NOTE: when specifying file-contents, options can be separated by ';' and/or newlines)
 */
PulsarParamsVector *
XLALPulsarParamsFromUserInput ( const LALStringVector *UserInput,	///< [in] user-input CSV list defining 'CW sources'
                                const LIGOTimeGPS *refTimeDef		///< [in] default reference time if refTime is not given
                                )
{
  XLAL_CHECK_NULL ( UserInput, XLAL_EINVAL );
  XLAL_CHECK_NULL ( UserInput->length > 0, XLAL_EINVAL );

  PulsarParamsVector *allSources = NULL;

  for ( UINT4 l = 0; l < UserInput->length; l ++ )
    {
      const char *thisInput = UserInput->data[l];

      if ( thisInput[0] != '{' )	// if it's an actual file-specification
        {
          LALStringVector *file_list;
          XLAL_CHECK_NULL ( ( file_list = XLALFindFiles ( &thisInput[0] )) != NULL, XLAL_EFUNC );
          UINT4 numFiles = file_list->length;
          for ( UINT4 i = 0; i < numFiles; i ++ )
            {
              PulsarParamsVector *sources_i;
              XLAL_CHECK_NULL ( (sources_i = XLALPulsarParamsFromFile ( file_list->data[i], refTimeDef )) != NULL, XLAL_EFUNC );

              XLAL_CHECK_NULL ( (allSources = XLALPulsarParamsVectorAppend ( allSources, sources_i )) != NULL, XLAL_EFUNC );
              XLALDestroyPulsarParamsVector ( sources_i );
            } // for i < numFiles

          XLALDestroyStringVector ( file_list );

        } // if file-pattern given
      else
        {
          UINT4 len = strlen(thisInput);
          XLAL_CHECK_NULL ( (thisInput[0] == '{') && (thisInput[len-1] == '}'), XLAL_EINVAL, "Invalid file-content input:\n%s\n", thisInput );
          char *buf;
          XLAL_CHECK_NULL ( (buf = XLALStringDuplicate ( &thisInput[1] )) != NULL, XLAL_EFUNC );
          len = strlen(buf);
          buf[len-1] = 0;	// remove trailing '}'

          LALParsedDataFile *cfgdata = NULL;
          XLAL_CHECK_NULL ( XLALParseDataFileContent ( &cfgdata, buf ) == XLAL_SUCCESS, XLAL_EFUNC );
          XLALFree ( buf );

          PulsarParamsVector *addSource;
          XLAL_CHECK_NULL ( (addSource = XLALCreatePulsarParamsVector ( 1 )) != NULL, XLAL_EFUNC );

          XLAL_CHECK_NULL ( XLALReadPulsarParams ( &addSource->data[0], cfgdata, NULL, refTimeDef ) == XLAL_SUCCESS, XLAL_EFUNC );
          strncpy ( addSource->data[0].name, "direct-string-input", sizeof(addSource->data[0].name) );

          XLAL_CHECK_NULL ( XLALCheckConfigFileWasFullyParsed( "{command-line}", cfgdata ) == XLAL_SUCCESS, XLAL_EINVAL );
          XLALDestroyParsedDataFile ( cfgdata );

          XLAL_CHECK_NULL ( (allSources = XLALPulsarParamsVectorAppend ( allSources, addSource )) != NULL, XLAL_EFUNC );
          XLALDestroyPulsarParamsVector ( addSource );

        } // if direct config-string given

    } // for l < len(UserInput)

  return allSources;

} // XLALPulsarParamsFromUserInput()

/**
 * Append the given PulsarParamsVector 'add' to the vector 'list' ( which can be NULL), return resulting list
 * with new elements 'add' appended at the end of 'list'.
 */
PulsarParamsVector *
XLALPulsarParamsVectorAppend ( PulsarParamsVector *list, const PulsarParamsVector *add )
{
  XLAL_CHECK_NULL ( add != NULL, XLAL_EINVAL );

  PulsarParamsVector *ret;
  if ( list == NULL )
    {
      XLAL_CHECK_NULL ( (ret = XLALCalloc ( 1, sizeof(*ret))) != NULL, XLAL_ENOMEM );
    }
  else
    {
      ret = list;
    }

  UINT4 oldlen = ret->length;
  UINT4 addlen = add->length;
  UINT4 newlen = oldlen + addlen;
  ret->length = newlen;
  XLAL_CHECK_NULL ( (ret->data = XLALRealloc ( ret->data, newlen * sizeof(ret->data[0]) )) != NULL, XLAL_ENOMEM );
  memcpy ( ret->data + oldlen, add->data, addlen * sizeof(ret->data[0]) );
  // we have to properly copy the 'name' string fields
  for ( UINT4 i = 0; i < addlen; i ++ )
    {
      strncpy ( ret->data[oldlen + i].name, add->data[i].name, sizeof(ret->data[oldlen + i].name) );
    }

  return ret;

} // XLALPulsarParamsVectorAppend()

/**
 * Write a PulsarParamsVector to a FITS file
 */
int
XLALFITSWritePulsarParamsVector ( FITSFile *file, const CHAR *tableName, const PulsarParamsVector *list )
{
  XLAL_CHECK ( file != NULL, XLAL_EFAULT );
  XLAL_CHECK ( tableName != NULL, XLAL_EFAULT );
  XLAL_CHECK ( list != NULL, XLAL_EFAULT );

  // Begin FITS table
  XLAL_CHECK ( XLALFITSTableOpenWrite ( file, tableName, "list of injections" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Describe FITS table
  XLAL_FITS_TABLE_COLUMN_BEGIN( PulsarParams );
  XLAL_CHECK ( XLAL_FITS_TABLE_COLUMN_ADD_ARRAY( file, CHAR, name ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, Amp.psi, "psi [rad]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, Amp.phi0, "phi0 [rad]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, Amp.aPlus, "aPlus" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, Amp.aCross, "aCross" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, GPSTime, Doppler.refTime, "refTime" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, Doppler.Alpha, "Alpha [rad]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, Doppler.Delta, "Delta [rad]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, Doppler.fkdot[0], "Freq [Hz]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  for ( size_t k = 1; k < PULSAR_MAX_SPINS; ++k ) {
    char col_name[64];
    snprintf( col_name, sizeof( col_name ), "f%zudot [Hz/s^%zu]", k, k );
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, Doppler.fkdot[k], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
  }
  XLAL_CHECK ( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, Doppler.asini, "orbitasini [s]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, Doppler.period, "orbitPeriod [s]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, Doppler.ecc, "orbitEcc" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, GPSTime, Doppler.tp, "orbitTp" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, Doppler.argp, "orbitArgp [rad]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, UINT4, Transient.type, "transientWindowType" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, UINT4, Transient.t0, "transientStartTime" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, UINT4, Transient.tau, "transientTau" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write FITS table
  for ( size_t i = 0; i < list->length; ++i ) {
    XLAL_CHECK( XLALFITSTableWriteRow( file, &list->data[i] ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

} // XLALFITSWritePulsarParamsVector()

/**
 * Destructor for a CWMFDataParams type
 *
 * \note: This is mostly useful for the SWIG wrappers
 */
void
XLALDestroyCWMFDataParams ( CWMFDataParams *params )
{
  if ( params ) {
    fflush(stdout);
    XLALDestroyMultiTimestamps ( params->multiTimestamps );
    XLALDestroyMultiREAL8TimeSeries ( params->inputMultiTS );
    XLALFree ( params );
  }
} // XLALDestroyCWMFDataParams()
