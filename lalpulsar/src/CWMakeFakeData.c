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
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
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

// ---------- local defines

// ---------- local macro definitions
#define SQ(x) ( (x) * (x) )
#define INIT_MEM(x) memset(&(x), 0, sizeof((x)))
// ---------- local type definitions

// ---------- empty initializers
const CWMFDataParams empty_CWMFDataParams;

// ---------- Global variables

// ---------- local prototypes
static UINT4 gcd (UINT4 numer, UINT4 denom);

// ==================== FUNCTION DEFINITIONS ====================


/**
 * Generate fake 'CW' data, returned either as SFTVector or REAL4Timeseries or both,
 * for given CW-signal ("pulsar") parameters and output parameters (frequency band etc)
 */
int
XLALCWMakeFakeMultiData ( MultiSFTVector **multiSFTs,			///< [out] pointer to optional SFT-vector for output
                          MultiREAL4TimeSeries **multiTseries,		///< [out] pointer to optional timeseries-vector for output
                          const PulsarParamsVector *injectionSources,	///< [in] array of sources inject
                          const CWMFDataParams *dataParams,		///< [in] parameters specifying the type of data to generate
                          const EphemerisData *edat			///< [in] ephemeris data
                          )
{
  XLAL_CHECK ( (multiSFTs == NULL) || ((*multiSFTs) == NULL ), XLAL_EINVAL );
  XLAL_CHECK ( (multiTseries == NULL) || ((*multiTseries) == NULL ), XLAL_EINVAL );
  XLAL_CHECK ( (multiSFTs != NULL) || (multiTseries != NULL), XLAL_EINVAL );

  XLAL_CHECK ( injectionSources != NULL, XLAL_EINVAL );
  XLAL_CHECK ( dataParams != NULL, XLAL_EINVAL );
  XLAL_CHECK ( edat != NULL, XLAL_EINVAL );

  const MultiLIGOTimeGPSVector *multiTimestamps = &(dataParams->multiTimestamps);

  // check multi-detector input
  XLAL_CHECK ( dataParams->detInfo.length >= 1, XLAL_EINVAL );
  UINT4 numDet = dataParams->detInfo.length;
  XLAL_CHECK ( multiTimestamps->length == numDet, XLAL_EINVAL, "Inconsistent number of IFOs: detInfo says '%d', multiTimestamps says '%d'\n", numDet, multiTimestamps->length );

  // check Tsft, consistent over detectors
  REAL8 Tsft = multiTimestamps->data[0]->deltaT;
  XLAL_CHECK ( Tsft > 0, XLAL_EINVAL, "Got invalid Tsft = %g must be > 0\n", Tsft );
  for ( UINT4 X=0; X < numDet; X ++ ) {
    XLAL_CHECK ( multiTimestamps->data[X]->deltaT == Tsft, XLAL_EINVAL, "Inconsistent Tsft, for Tsft[X=0]=%g, while Tsft[X=%d]=%g\n", Tsft, X, multiTimestamps->data[X]->deltaT );
  }

  // ----- prepare output containers, as required
  MultiSFTVector *outMSFTs = NULL;
  MultiREAL4TimeSeries *outMTS = NULL;
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
      /* detector params */
      CWMFDataParams dataParamsX = (*dataParams); // struct-copy for general settings
      dataParamsX.detInfo.length = 1;
      dataParamsX.detInfo.sites[0] = dataParams->detInfo.sites[X];
      dataParamsX.detInfo.sqrtSn[0] = dataParams->detInfo.sqrtSn[X];
      MultiLIGOTimeGPSVector mTimestamps = empty_MultiLIGOTimeGPSVector;
      mTimestamps.length = 1;
      mTimestamps.data = &(multiTimestamps->data[X]); // such that pointer mTimestamps.data[0] = multiTimestamps->data[X]
      dataParamsX.multiTimestamps = mTimestamps;
      dataParamsX.randSeed = dataParams->randSeed + X;	// increase seed in deterministic way: allows comparison w mfd_v4 !!

      SFTVector **svp = NULL;
      REAL4TimeSeries **tsp = NULL;
      if ( outMSFTs != NULL ) {
        svp = &(outMSFTs->data[X]);
      }
      if ( outMTS != NULL ) {
        tsp = &(outMTS->data[X]);
      }
      XLAL_CHECK ( XLALCWMakeFakeData ( svp, tsp, injectionSources, &dataParamsX, edat ) == XLAL_SUCCESS, XLAL_EFUNC );

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
 * work, but same input API. The input detector-arrays must all contain
 * only a single detector, otherwise an error is returned.
 */
int
XLALCWMakeFakeData ( SFTVector **SFTvect,
                     REAL4TimeSeries **Tseries,
                     const PulsarParamsVector *injectionSources,
                     const CWMFDataParams *dataParams,
                     const EphemerisData *edat
                     )
{
  XLAL_CHECK ( (SFTvect == NULL) || ((*SFTvect) == NULL ), XLAL_EINVAL );
  XLAL_CHECK ( (Tseries == NULL) || ((*Tseries) == NULL ), XLAL_EINVAL );
  XLAL_CHECK ( (SFTvect != NULL) || (Tseries != NULL), XLAL_EINVAL );

  XLAL_CHECK ( injectionSources != NULL, XLAL_EINVAL );
  XLAL_CHECK ( dataParams != NULL, XLAL_EINVAL );
  XLAL_CHECK ( edat != NULL, XLAL_EINVAL );
  XLAL_CHECK ( dataParams->detInfo.length ==1, XLAL_EINVAL );
  XLAL_CHECK ( dataParams->multiTimestamps.length == 1, XLAL_EINVAL );
  XLAL_CHECK ( dataParams->multiTimestamps.data[0] != NULL, XLAL_EINVAL );

  // initial default values fMin, sampling rate from caller input
  REAL8 fMin  = dataParams->fMin;
  REAL8 fBand = dataParams->Band;
  REAL8 fSamp = 2.0 * fBand;

  const LIGOTimeGPSVector *timestamps = dataParams->multiTimestamps.data[0];
  const LALDetector *site = &dataParams->detInfo.sites[0];
  REAL8 Tsft = timestamps->deltaT;

  // if SFT output requested: need *effective* fMin and Band consistent with SFT bins
  if ( SFTvect )
    {
      UINT4 firstBinEff, numBinsEff;
      XLAL_CHECK ( XLALFindCoveringSFTBins ( &firstBinEff, &numBinsEff, dataParams->fMin, dataParams->Band, Tsft ) == XLAL_SUCCESS, XLAL_EFUNC );

      REAL8 fBand_eff = (numBinsEff - 1.0) / Tsft;
      REAL8 fMin_eff  = firstBinEff / Tsft;
      REAL8 fMax = fMin + dataParams->Band;
      REAL8 fMax_eff = fMin_eff + fBand_eff;
      if ( (fMin_eff != fMin) || (fBand_eff != dataParams->Band ) ) {
        XLALPrintWarning("Caller asked for Band [%.16g, %.16g] Hz, effective SFT-Band produced is [%.16g, %.16g] Hz\n",
                         fMin, fMax, fMin_eff, fMax_eff );
      }
      fMin = fMin_eff;		// (potentially) lower minimal frequency to fit SFT bins
      fBand = fBand_eff;
      fSamp = 2.0 * fBand_eff;	// (potentially) higher sampling rate required to fit SFT bins
    } // if SFT-output

  /* characterize the output time-series */
  UINT4 n0_fSamp = (UINT4) round ( Tsft * fSamp );

  // by construction, fSamp * Tsft = integer, but if there are gaps between SFTs,
  // then we might have to sample at higher rates in order for all SFT-boundaries to
  // fall on exact timesteps of the timeseries.
  // ie we start from fsamp0 = n0_fSamp/Tsft, and then try to find the smallest
  // n1_fSamp >= n0_fSamp, such that for fsamp1 = n1_fSamp/Tsft, for all gaps i: Dt_i * fsamp1 = int
  UINT4 n1_fSamp;
  XLAL_CHECK ( XLALFindSmallestValidSamplingRate ( &n1_fSamp, n0_fSamp, timestamps ) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( n1_fSamp != n0_fSamp )
    {
      REAL8 fSamp1 = n1_fSamp / Tsft;	// increased sampling rate to fit all gaps
      XLALPrintWarning ( "GAPS: Initial SFT sampling frequency fSamp0= %d/%.0f = %g had to be increased to fSamp1 = %d/%.0f = %g\n",
                         n0_fSamp, Tsft, fSamp, n1_fSamp, Tsft, fSamp1 );
      fSamp = fSamp1;
    } // if higher effective sampling rate required

   /* ----- start-time and duration ----- */
  LIGOTimeGPS firstGPS = timestamps->data[0];
  REAL8 firstGPS_REAL8 = XLALGPSGetREAL8 ( &firstGPS );
  LIGOTimeGPS lastGPS  = timestamps->data [ timestamps->length - 1 ];
  REAL8 lastGPS_REAL8 = XLALGPSGetREAL8 ( &lastGPS );
  XLALGPSAdd( &lastGPS, Tsft );
  REAL8 duration = XLALGPSDiff ( &lastGPS, &firstGPS );

  REAL4TimeSeries *Tseries_sum = NULL;
  XLAL_CHECK ( (Tseries_sum = XLALGenerateCWSignalTS ( &injectionSources->data[0], site, firstGPS, duration, fSamp, fMin, edat )) != NULL, XLAL_EFUNC );

  UINT4 numPulsars = injectionSources->length;
  for ( UINT4 iInj = 1; iInj < numPulsars; iInj ++ )
    {
      // for all but the first time-series, we truncate any transient-CW timeseries to the actual support
      // of the transient signal, in order to make the generation more efficient, these 'partial timeseries'
      // will then be added to the full timeseries of the 1st signal (with iInj=0)
      const PulsarParams *pulsarParams = &( injectionSources->data[iInj] );
      UINT4 t0, t1;
      XLAL_CHECK ( XLALGetTransientWindowTimespan ( &t0, &t1, pulsarParams->Transient ) == XLAL_SUCCESS, XLAL_EFUNC );

      // use latest possible start-time: max(t0,firstGPS), but not later than than lastGPS
      LIGOTimeGPS signalStartGPS; INIT_MEM ( signalStartGPS );
      if ( t0 <= firstGPS_REAL8 ) {
        signalStartGPS = firstGPS;
      } else if ( t0 >= lastGPS_REAL8 ) {
        signalStartGPS = lastGPS;
      }
      else {
        signalStartGPS.gpsSeconds = t0;
      }

      // use earliest possible end-time: min(t1,lastGPS), but not earlier than firstGPS
      LIGOTimeGPS signalEndGPS; INIT_MEM ( signalEndGPS );
      if ( t1 >= lastGPS_REAL8 ) {
        signalEndGPS = lastGPS;
      } else if ( t1 <= firstGPS_REAL8 ) {
        signalEndGPS = firstGPS;
      } else {
        signalEndGPS.gpsSeconds = t1;
      }
      REAL8 signalDuration = XLALGPSDiff ( &signalEndGPS, &signalStartGPS );
      XLAL_CHECK ( signalDuration >= 0, XLAL_EFAILED, "Something went wrong, got negative signal duration = %g\n", signalDuration );
      if ( signalDuration > 0 )	// only need to do sth if transient-window had finite overlap with output TS
        {
          REAL4TimeSeries *Tseries_i = NULL;
          XLAL_CHECK ( (Tseries_i = XLALGenerateCWSignalTS ( pulsarParams, site, signalStartGPS, signalDuration, fSamp, fMin, edat )) != NULL, XLAL_EFUNC );

          XLAL_CHECK ( (Tseries_sum = XLALAddREAL4TimeSeries ( Tseries_sum, Tseries_i )) != NULL, XLAL_EFUNC );
          XLALDestroyREAL4TimeSeries ( Tseries_i );
        }
    } // for iInj = 1 ... (numPulsars-1)

  /* add Gaussian noise if requested */
  REAL8 sqrtSn = dataParams->detInfo.sqrtSn[0];
  if ( sqrtSn > 0)
    {
      REAL8 noiseSigma = sqrtSn * sqrt ( fBand );
      XLAL_CHECK ( XLALAddGaussianNoise ( Tseries_sum, noiseSigma, dataParams->randSeed ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

  /*---------------------------------------------
   * turn this timeseries into SFTs, if requested
   *---------------------------------------------*/
  if ( SFTvect )
    {
      // Prepare windowing of time series for SFTs
      REAL4Window *window = NULL;
      if ( dataParams->SFTWindowType )
        {
          REAL8 dt = Tseries_sum->deltaT;
          UINT4 numTimesteps = round ( Tsft / dt );	/* number of time-samples in an Tsft (should be exact integer!) */

          XLAL_CHECK ( (window = XLALCreateNamedREAL4Window ( dataParams->SFTWindowType, dataParams->SFTWindowBeta, numTimesteps )) != NULL, XLAL_EFUNC );
        } // if uvar->SFTwindowName

      SFTParams sftParams = empty_SFTParams;
      sftParams.Tsft = Tsft;
      sftParams.timestamps = timestamps;
      sftParams.noiseSFTs = NULL;	// not used here any more!
      sftParams.window = window;

      /* compute SFTs from timeseries */
      SFTVector *sftVect;
      XLAL_CHECK ( (sftVect = XLALSignalToSFTs (Tseries_sum, &sftParams)) != NULL, XLAL_EFUNC );

      // extract effective band from this, if neccessary (ie if faster-sampled output SFTs)
      if ( n1_fSamp != n0_fSamp )
        {
          XLAL_CHECK ( ((*SFTvect) = XLALExtractBandFromSFTVector ( sftVect, fMin, fBand )) != NULL, XLAL_EFUNC );
          XLALDestroySFTVector ( sftVect );
        }
      else
        {
          (*SFTvect) = sftVect;
        }

      XLALDestroyREAL4Window ( window );

    } // if SFTvect

  // return timeseries if requested
  if ( Tseries )
    {
      (*Tseries) = Tseries_sum;
    }
  else
    {
      XLALDestroyREAL4TimeSeries ( Tseries_sum );
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
                         const EphemerisData *edat		///< ephemeris data
                         )
{
  XLAL_CHECK_NULL ( pulsarParams != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( site != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( duration > 0, XLAL_EDOM );
  XLAL_CHECK_NULL ( fSamp > 0, XLAL_EDOM );
  XLAL_CHECK_NULL ( fHet >= 0, XLAL_EDOM );

  // translate amplitude params
  REAL8 h0     = pulsarParams->Amp.h0;
  REAL8 cosi   = pulsarParams->Amp.cosi;
  REAL8 aPlus  = 0.5 * h0 * ( 1.0 + SQ(cosi) );
  REAL8 aCross = h0 * cosi;
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
  PulsarSignalParams params = empty_PulsarSignalParams;
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
  params.orbit                     = pulsarParams->Doppler.orbit;
  params.transfer                  = NULL;
  params.ephemerides               = edat;
  params.fHeterodyne               = fHet;

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
  XLAL_CHECK ( TsftREAL == round(TsftREAL), XLAL_EDOM, "Only exact integer-second Tsft allowed, got Tsft = %g s\n", TsftREAL );
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

      XLALPrintInfo ("gap_i = %d s, remainder wrt Tsft=%d s = %d s: ", gap_i0, Tsft, gap_i );
      if ( (gap_i * nCur) % Tsft == 0 ) {
        XLALPrintInfo ("Fits exactly with fSamp = %d / Tsft = %g\n", nCur, nCur / TsftREAL );
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
      XLALPrintInfo ("Need to increase to fSamp = %d / Tsft = %g\n", nNew, nNew / TsftREAL );

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
 * Create PulsarParamsVector for numPulsars
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

  UINT4 numPulsars = ppvect->length;
  if ( ppvect->data != NULL )
    {
      for ( UINT4 i = 0 ; i < numPulsars; i ++ )
        {
          BinaryOrbitParams *orbit = ppvect->data[i].Doppler.orbit;
          XLALFree ( orbit );
          XLALFree ( ppvect->data[i].name );
        } // for i < numPulsars
      XLALFree ( ppvect->data );
    }

  XLALFree ( ppvect );

  return;
} // XLALDestroyPulsarParamsVector()

/**
 * Destructor for PulsarParams types
 */
void
XLALDestroyPulsarParams ( PulsarParams *params )
{
  if ( params == NULL ) {
    return;
  }
  XLALFree ( params->name );
  XLALFree ( params->Doppler.orbit );
  XLALFree ( params );

  return;

} // XLALDestroyPulsarParams()


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
 * Alpha, Delta, Freq, refTime
 *
 * OPTIONAL:
 * f1dot, f2dot, f3dot, f4dot, f5dot, f6dot
 * {h0, cosi} or {aPlus, aCross}, psi, phi0
 * transientWindowType, transientStartTime, transientTauDays
 *
 * Other config-variables found in the file will ... ?? error or accept?
 */
int
XLALReadPulsarParams ( PulsarParams *pulsarParams,	///< [out] pulsar parameters to fill in from config string
                       const LALParsedDataFile *cfgdata,///< [in] pre-parsed "SourceParamsIO" config-file contents
                       const CHAR *secName		///< [in] section-name to use from config-file string (can be NULL)
                       )
{
  XLAL_CHECK ( pulsarParams != NULL, XLAL_EINVAL );
  XLAL_CHECK ( pulsarParams->Doppler.orbit == NULL, XLAL_EINVAL );
  XLAL_CHECK ( cfgdata != NULL, XLAL_EINVAL );

  INIT_MEM ( (*pulsarParams) );	// wipe input struct clean

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

  if ( have_aPlus )	/* translate A_{+,x} into {h_0, cosi} */
    {
      XLAL_CHECK ( fabs ( aCross ) <= aPlus, XLAL_EDOM, "ERROR: |aCross| (= %g) must be <= aPlus (= %g).\n", fabs(aCross), aPlus );
      REAL8 disc = sqrt ( SQ(aPlus) - SQ(aCross) );
      h0 = aPlus + disc;
      if ( h0 > 0 ) {
        cosi = aCross / h0;	// avoid division by 0!
      }
    } //if {aPlus, aCross}

  XLAL_CHECK ( h0 >= 0, XLAL_EDOM );
  pulsarParams->Amp.h0 	= h0;

  XLAL_CHECK ( (cosi >= -1) && (cosi <= 1), XLAL_EDOM );
  pulsarParams->Amp.cosi= cosi;

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
  REAL8 refTime_GPS = 0; BOOLEAN have_refTime;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &refTime_GPS, cfgdata, secName, "refTime", &have_refTime ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( have_refTime, XLAL_EINVAL );

  XLAL_CHECK ( XLALGPSSetREAL8 ( & pulsarParams->Doppler.refTime, refTime_GPS ) != NULL, XLAL_EFUNC );

  // ----- Alpha
  REAL8 Alpha_Rad = 0; BOOLEAN have_Alpha;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &Alpha_Rad, cfgdata, secName, "Alpha", &have_Alpha ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( have_Alpha, XLAL_EINVAL );

  XLAL_CHECK ( (Alpha_Rad >= 0) && (Alpha_Rad < LAL_TWOPI), XLAL_EDOM );
  pulsarParams->Doppler.Alpha = Alpha_Rad;

  // ----- Delta
  REAL8 Delta_Rad = 0; BOOLEAN have_Delta;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &Delta_Rad, cfgdata, secName, "Delta", &have_Delta ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( have_Delta, XLAL_EINVAL );

  XLAL_CHECK ( (Delta_Rad >= -LAL_PI_2) && (Delta_Rad <= LAL_PI_2), XLAL_EDOM );
  pulsarParams->Doppler.Delta = Delta_Rad;

  // ----- fkdot
  // Freq
  REAL8 Freq = 0; BOOLEAN have_Freq;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &Freq, cfgdata, secName, "Freq", &have_Freq ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( have_Freq, XLAL_EINVAL );

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
  REAL8 orbitTpSSB = 0;	BOOLEAN have_orbitTpSSB;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &orbitTpSSB, cfgdata, secName, "orbitTpSSB", &have_orbitTpSSB ) == XLAL_SUCCESS, XLAL_EFUNC );
  REAL8 orbitArgp = 0; 	BOOLEAN have_orbitArgp;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &orbitArgp, cfgdata, secName, "orbitArgp", &have_orbitArgp ) == XLAL_SUCCESS, XLAL_EFUNC );
  REAL8 orbitasini = 0; BOOLEAN have_orbitasini;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &orbitasini, cfgdata, secName, "orbitasini", &have_orbitasini ) == XLAL_SUCCESS, XLAL_EFUNC );
  REAL8 orbitEcc = 0;  	BOOLEAN have_orbitEcc;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &orbitEcc, cfgdata, secName, "orbitEcc", &have_orbitEcc ) == XLAL_SUCCESS, XLAL_EFUNC );
  REAL8 orbitPeriod = 0;BOOLEAN have_orbitPeriod;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &orbitPeriod, cfgdata, secName, "orbitPeriod", &have_orbitPeriod ) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( have_orbitasini || have_orbitEcc || have_orbitPeriod || have_orbitArgp || have_orbitTpSSB )
    {
      XLAL_CHECK ( orbitasini >= 0, XLAL_EDOM );
      XLAL_CHECK ( (orbitasini == 0) || ( have_orbitEcc && have_orbitPeriod && have_orbitArgp && have_orbitTpSSB ), XLAL_EINVAL );
      XLAL_CHECK ( (orbitEcc >= 0) && (orbitEcc <= 1), XLAL_EDOM );

      BinaryOrbitParams *orbit;
      XLAL_CHECK ( ( orbit = XLALCalloc ( 1, sizeof(BinaryOrbitParams))) != NULL, XLAL_ENOMEM );

      /* fill in orbital parameter structure */
      XLAL_CHECK ( XLALGPSSetREAL8 ( &(orbit->tp), orbitTpSSB ) != NULL, XLAL_EFUNC );
      orbit->argp 	= orbitArgp;
      orbit->asini 	= orbitasini;
      orbit->ecc 	= orbitEcc;
      orbit->period 	= orbitPeriod;

      pulsarParams->Doppler.orbit = orbit;
    } // if have non-trivial orbit

  // ---------- transientWindow_t ----------

  // ----- type
  char *transientWindowType = NULL; BOOLEAN have_transientWindowType;
  XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &transientWindowType, cfgdata, secName, "transientWindowType", &have_transientWindowType ) == XLAL_SUCCESS, XLAL_EFUNC );
  // ----- t0
  REAL8 transientStartTime = 0; BOOLEAN have_transientStartTime;
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &transientStartTime, cfgdata, secName, "transientStartTime", &have_transientStartTime ) == XLAL_SUCCESS, XLAL_EFUNC );
  // ----- tau
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
      XLAL_CHECK ( have_transientStartTime && have_transientTauDays, XLAL_EINVAL );
      XLAL_CHECK ( transientStartTime >= 0, XLAL_EDOM );
      XLAL_CHECK ( transientTauDays > 0, XLAL_EDOM );

      pulsarParams->Transient.t0   = (UINT4) transientStartTime;
      pulsarParams->Transient.tau  = (UINT4) ( transientTauDays * LAL_DAYSID_SI );
    } /* if transient window != none */
  else
    {
      XLAL_CHECK ( !(have_transientStartTime || have_transientTauDays), XLAL_EINVAL );
    }

  return XLAL_SUCCESS;
} // XLALParsePulsarParams()


/**
 * Parse a given 'CWsources' config file for PulsarParams, return vector
 * of all pulsar definitions found [using sections]
 */
PulsarParamsVector *
XLALPulsarParamsFromFile ( const char *fname 		///< [in] 'CWsources' config file name
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
      XLAL_CHECK_NULL ( XLALReadPulsarParams ( &sources->data[i], cfgdata, sec_i ) == XLAL_SUCCESS, XLAL_EFUNC );

      // ----- source naming convention: 'filename:section'
      char *name;
      size_t len = strlen(fname) + strlen(sections->data[i]) + 2;
      XLAL_CHECK_NULL ( (name = XLALCalloc(1, len)) != NULL, XLAL_ENOMEM );
      sprintf ( name, "%s:%s", fname, sections->data[i] );
      sources->data[i].name = name;

    } // for i < numPulsars

  XLALDestroyStringVector ( sections );
  XLALDestroyParsedDataFile ( cfgdata );

  return sources;

} // XLALPulsarParamsFromFile()

/**
 * Function to determine the PulsarParamsVector input from a user-input defining CW sources.
 *
 * This option supports a dual-type feature: if the input starts with a '@', then
 * it determines a list of input config-files to be parsed by XLALFindFiles(),
 * otherwise the input will be interpreted as a config-file string directly (with
 * options separated by ';' and/or newlines)
 */
PulsarParamsVector *
XLALPulsarParamsFromUserInput ( const char *UserInput		///< [in] user-input string defining 'CW sources'
                                )
{
  XLAL_CHECK_NULL ( UserInput, XLAL_EINVAL );

  PulsarParamsVector *sources = NULL;

  if ( UserInput[0] == '@' )
    {
      LALStringVector *file_list;
      XLAL_CHECK_NULL ( ( file_list = XLALFindFiles ( &UserInput[1] )) != NULL, XLAL_EFUNC );
      UINT4 numFiles = file_list->length;
      for ( UINT4 i = 0; i < numFiles; i ++ )
        {
          PulsarParamsVector *sources_i;
          XLAL_CHECK_NULL ( (sources_i = XLALPulsarParamsFromFile ( file_list->data[i] )) != NULL, XLAL_EFUNC );

          if ( sources == NULL )
            {
              sources = sources_i;
            }
          else
            {
              UINT4 oldlen = sources->length;
              UINT4 addlen = sources_i->length;
              UINT4 newlen = oldlen + addlen;
              sources->length = newlen;
              XLAL_CHECK_NULL ( (sources->data = XLALRealloc ( sources->data, newlen * sizeof(sources->data[0]) )) != NULL, XLAL_ENOMEM );
              memcpy ( sources->data + oldlen, sources_i->data, addlen * sizeof(sources->data[0]) );
              XLALFree ( sources_i->data );
              XLALFree ( sources_i );
            }
        } // for i < numFiles

      XLALDestroyStringVector ( file_list );

    } // if file-name spec given
  else
    {
      LALParsedDataFile *cfgdata = NULL;
      XLAL_CHECK_NULL ( XLALParseDataFileContent ( &cfgdata, UserInput ) == XLAL_SUCCESS, XLAL_EFUNC );

      XLAL_CHECK_NULL ( (sources = XLALCreatePulsarParamsVector ( 1 )) != NULL, XLAL_EFUNC );

      XLAL_CHECK_NULL ( XLALReadPulsarParams ( &sources->data[0], cfgdata, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_NULL ( (sources->data[0].name = XLALStringDuplicate ( "from-commandline" )) != NULL, XLAL_EFUNC );
      XLALDestroyParsedDataFile ( cfgdata );

    } // if direct config-string given

  return sources;

} // XLALPulsarParamsFromUserInput()
