//
// Copyright (C) 2012, 2013, 2014 Karl Wette
// Copyright (C) 2009 Chris Messenger, Reinhard Prix, Pinkesh Patel, Xavier Siemens, Holger Pletsch
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
//

// This file implements the F-statistic resampling algorithm. It is not compiled directly, but
// included from ComputeFstat.c

#include <lal/LogPrintf.h>

// ========== Resamp internals ==========

// ----- local constants
static LALUnit emptyLALUnit;

// ----- local types ----------
struct tagFstatInput_Resamp {
  MultiCOMPLEX8TimeSeries  *multiTimeSeries_DET;	// input SFTs converted into a heterodyned timeseries

  UINT4 prev_numSamples_SRC;				// keep track of previous SRC-frame samples (ie dFreqOut)
  UINT4 prev_numFreqBinsOut;				// keep track of previous number of output frequency bins
  // ----- workspace ----------
  MultiCOMPLEX8TimeSeries *ws_multiTimeSeries_SRC;	// keep allocated memory for modifying this SRC-frame TS

  MultiCOMPLEX8TimeSeries *ws_multiFa_SRC;
  MultiCOMPLEX8TimeSeries *ws_multiFb_SRC;

  COMPLEX8Vector *ws_outaX;				// hold results of FTT
  COMPLEX8Vector *ws_outbX;
  ComplexFFTPlan *ws_fftplan;
  COMPLEX8 *ws_Fa_k;
  COMPLEX8 *ws_Fb_k;

  // ----- buffering -----
  PulsarDopplerParams prev_doppler;			// buffering: previous phase-evolution ("doppler") parameters

  MultiAMCoeffs *prev_multiAMcoef;			// buffering: previous AM-coeffs, unique to skypos
  MultiSSBtimes *prev_multiSSBsky;			// buffering: previous sky-only multiSSB times, depends on skypos and reftime

  MultiLIGOTimeGPSVector *prev_multiTimestamps_SRC;	// buffering: SFT timestamps translated into SRC frame
  MultiCOMPLEX8TimeSeries *prev_multiTimeSeries_SRC;	// buffering: multi-detector SRC-frame timeseries
};


// ----- local prototypes ----------
static int
XLALAntennaWeightCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *ax, COMPLEX8TimeSeries *bx,
                                      const COMPLEX8TimeSeries *timeseries,
                                      const AMCoeffs *AMcoef,
                                      const LIGOTimeGPSVector *TS );

static int
XLALAntennaWeightMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *m_ax, MultiCOMPLEX8TimeSeries *m_bx,
                                           const MultiCOMPLEX8TimeSeries *multiTimeSeries,
                                           const MultiAMCoeffs *multiAMcoef,
                                           const MultiLIGOTimeGPSVector *multiTS );

static int
XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries **mTimeSeries_SRC,
                                                 MultiLIGOTimeGPSVector **mTimestamps_SRC,
                                                 const MultiCOMPLEX8TimeSeries *mTimeSeries_DET,
                                                 const MultiLIGOTimeGPSVector *mTimestamps_DET,
                                                 const MultiSSBtimes *mSRC_timing,
                                                 const REAL8 deltaF
                                                 );

static int
XLALBarycentricResampleCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *TimeSeries_SRC,
                                            LIGOTimeGPSVector *Timestamps_SRC,
                                            const COMPLEX8TimeSeries *TimeSeries_DET,
                                            const LIGOTimeGPSVector *Timestamps_DET,
                                            const SSBtimes *SRC_timing
                                            );

// ==================== function definitions

static void
DestroyFstatInput_Resamp ( FstatInput_Resamp* resamp )
{
  XLALDestroyMultiCOMPLEX8TimeSeries (resamp->multiTimeSeries_DET );

  // ----- free workspace
  XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->ws_multiTimeSeries_SRC );

  XLALDestroyCOMPLEX8Vector ( resamp->ws_outaX );
  XLALDestroyCOMPLEX8Vector ( resamp->ws_outbX );

  XLALDestroyCOMPLEX8FFTPlan ( resamp->ws_fftplan );

  XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->ws_multiFa_SRC );
  XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->ws_multiFb_SRC );

  XLALFree ( resamp->ws_Fa_k );
  XLALFree ( resamp->ws_Fb_k );

  // ----- free buffer
  XLALDestroyMultiAMCoeffs ( resamp->prev_multiAMcoef );
  XLALDestroyMultiSSBtimes ( resamp->prev_multiSSBsky );
  XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->prev_multiTimeSeries_SRC );
  XLALDestroyMultiTimestamps ( resamp->prev_multiTimestamps_SRC );

  XLALFree ( resamp );

  return;
} // DestroyFstatInput_Resamp()

static int
SetupFstatInput_Resamp ( FstatInput_Resamp *resamp,
                         const FstatInput_Common *common,
                         MultiSFTVector *multiSFTs
                         )
{
  // Check input
  XLAL_CHECK(common != NULL, XLAL_EFAULT);
  XLAL_CHECK(resamp != NULL, XLAL_EFAULT);
  XLAL_CHECK(multiSFTs != NULL, XLAL_EFAULT);

  // Convert SFTs into heterodyned complex timeseries [in detector frame]
  /* generate multiple coincident timeseries - one for each detector spanning start -> end */
  /* we need each timeseries to span the exact same amount of time and to start at the same time */
  /* because for the multi-detector Fstat we need frequency bins to be coincident */
  /* The memory allocated here is freed when the buffer is cleared in the calling program */
  /* generate complex heterodyned timeseries from the input SFTs */
  XLAL_CHECK ( (resamp->multiTimeSeries_DET = XLALMultiSFTVectorToCOMPLEX8TimeSeries ( multiSFTs )) != NULL, XLAL_EFUNC );

  XLALDestroyMultiSFTVector ( multiSFTs );	// don't need them SFTs any more ...

  return XLAL_SUCCESS;

} // SetupFstatInput_Resamp()


static int
GetFstatExtraBins_Resamp ( FstatInput_Resamp* resamp )
{
  XLAL_CHECK(resamp != NULL, XLAL_EFAULT);
  return 8;	// use 8 extra bins to give better agreement with LALDemod(w Dterms=8) near the boundaries
} // GetFstatExtraBins_Resamp()


static int
ComputeFstat_Resamp ( FstatResults* Fstats,
                      const FstatInput_Common *common,
                      FstatInput_Resamp* resamp
                      )
{
  // Check input
  XLAL_CHECK ( Fstats != NULL, XLAL_EFAULT );
  XLAL_CHECK ( common != NULL, XLAL_EFAULT );
  XLAL_CHECK ( resamp != NULL, XLAL_EFAULT );

  const FstatQuantities whatToCompute = Fstats->whatWasComputed;
  XLAL_CHECK ( !(whatToCompute & FSTATQ_ATOMS_PER_DET), XLAL_EINVAL, "Resampling does not currently support atoms per detector" );

  // ----- handy shortcuts ----------
  PulsarDopplerParams thisPoint = Fstats->doppler;
  const MultiCOMPLEX8TimeSeries *multiTimeSeries_DET = resamp->multiTimeSeries_DET;
  REAL8 dt_DET = multiTimeSeries_DET->data[0]->deltaT;
  UINT4 numSamples_DET = multiTimeSeries_DET->data[0]->data->length;
  REAL8 Tspan_DET = numSamples_DET * dt_DET;

  MultiAMCoeffs *multiAMcoef;
  MultiLIGOTimeGPSVector *multiTimestamps_SRC = NULL;
  MultiCOMPLEX8TimeSeries *multiTimeSeries_SRC = NULL;

  /* determine resampled timeseries parameters */
  REAL8 dFreqOut = ( Fstats->dFreq > 0 ) ? Fstats->dFreq : 1.0 / Tspan_DET;

  REAL8 Tspan_SRC = 1.0 / dFreqOut;                                       /* the effective observation time based on the requested frequency resolution (for zero padding) */
  UINT4 numSamples_SRC = (UINT4) ceil ( Tspan_SRC / dt_DET );      /* we use ceil() so that we artificially widen the band rather than reduce it */

  // ============================== BEGIN: handle buffering =============================
  // ----- is it the same skyposition and reference time and frequency-resolution as last call ? -----
  if ( (resamp->prev_doppler.Alpha == thisPoint.Alpha) &&
       (resamp->prev_doppler.Delta == thisPoint.Delta) &&
       (XLALGPSDiff ( &resamp->prev_doppler.refTime, &thisPoint.refTime ) == 0 ) &&
       (resamp->prev_numSamples_SRC == numSamples_SRC)
       )
    {
      multiAMcoef = resamp->prev_multiAMcoef;
      MultiSSBtimes *multiSSBsky = resamp->prev_multiSSBsky;
      // ----- is it the same binary-orbital parameters as last call ? -----
      if ( (resamp->prev_doppler.asini == thisPoint.asini) &&
           (resamp->prev_doppler.period == thisPoint.period) &&
           (resamp->prev_doppler.ecc == thisPoint.ecc) &&
           (XLALGPSCmp( &resamp->prev_doppler.tp, &thisPoint.tp )==0 ) &&
           (resamp->prev_doppler.argp == thisPoint.argp)
           )
        { // ----- no changes in sky + binary ==> reuse SRC-frame timeseries and SFT timestamps
          multiTimestamps_SRC = resamp->prev_multiTimestamps_SRC;
          multiTimeSeries_SRC = resamp->prev_multiTimeSeries_SRC;
        }
      else
        {  // ----- same skypos but changes in binary-orbital parameters: recompute just those

          if ( thisPoint.asini > 0 )
            {
              // add binary time corrections to the SSB time delays and SSB time derivitive
              MultiSSBtimes *multiBinary = NULL;
              XLAL_CHECK ( XLALAddMultiBinaryTimes ( &multiBinary, multiSSBsky, &thisPoint ) == XLAL_SUCCESS, XLAL_EFUNC );
              XLAL_CHECK ( XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( &multiTimeSeries_SRC, &multiTimestamps_SRC, multiTimeSeries_DET, common->timestamps, multiBinary, dFreqOut)
                           == XLAL_SUCCESS, XLAL_EFUNC );
              XLALDestroyMultiSSBtimes ( multiBinary );
            } // if asini > 0
          else
            {
              XLAL_CHECK ( XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( &multiTimeSeries_SRC, &multiTimestamps_SRC, multiTimeSeries_DET, common->timestamps, multiSSBsky, dFreqOut)
                           == XLAL_SUCCESS, XLAL_EFUNC );
            } // if asini==0

          // ----- store new weighted SRC timeseries in buffer ----------
          resamp->prev_doppler = thisPoint;

          XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->prev_multiTimeSeries_SRC );
          resamp->prev_multiTimeSeries_SRC = multiTimeSeries_SRC;

          XLALDestroyMultiTimestamps ( resamp->prev_multiTimestamps_SRC );
          resamp->prev_multiTimestamps_SRC = multiTimestamps_SRC;
        } // end: if changed binary parameters

    } // end: if identical sky-position and reftime
  else
    { // ----- changed sky-position: compute SSB + AMcoef for this skyposition
      SkyPosition skypos;
      skypos.system = COORDINATESYSTEM_EQUATORIAL;
      skypos.longitude = thisPoint.Alpha;
      skypos.latitude  = thisPoint.Delta;
      const MultiDetectorStateSeries *multiDetStates = common->detectorStates;
      const MultiNoiseWeights *multiWeights = common->noiseWeights;

      MultiSSBtimes *multiSSBsky;
      XLAL_CHECK ( (multiSSBsky = XLALGetMultiSSBtimes ( multiDetStates, skypos, thisPoint.refTime, common->SSBprec )) != NULL, XLAL_EFUNC );

      if ( thisPoint.asini > 0 )
        { // add binary time corrections to the SSB time delays and SSB time derivitive
          MultiSSBtimes *multiBinary = NULL;
          XLAL_CHECK ( XLALAddMultiBinaryTimes ( &multiBinary, multiSSBsky, &thisPoint ) == XLAL_SUCCESS, XLAL_EFUNC );
          XLAL_CHECK ( XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( &multiTimeSeries_SRC, &multiTimestamps_SRC, multiTimeSeries_DET, common->timestamps, multiBinary, dFreqOut)
                       == XLAL_SUCCESS, XLAL_EFUNC );
          XLALDestroyMultiSSBtimes ( multiBinary );
        } // if asini > 0
      else
        {
          XLAL_CHECK ( XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( &multiTimeSeries_SRC, &multiTimestamps_SRC, multiTimeSeries_DET, common->timestamps, multiSSBsky, dFreqOut)
                       == XLAL_SUCCESS, XLAL_EFUNC );
        } // if asini==0

      XLAL_CHECK ( (multiAMcoef = XLALComputeMultiAMCoeffs ( multiDetStates, multiWeights, skypos )) != NULL, XLAL_EFUNC );

      // ----- store everything in buffer ----------
      resamp->prev_doppler = thisPoint;

      XLALDestroyMultiAMCoeffs ( resamp->prev_multiAMcoef );
      resamp->prev_multiAMcoef = multiAMcoef;

      XLALDestroyMultiSSBtimes ( resamp->prev_multiSSBsky );
      resamp->prev_multiSSBsky = multiSSBsky;

      XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->prev_multiTimeSeries_SRC );
      resamp->prev_multiTimeSeries_SRC = multiTimeSeries_SRC;

      XLALDestroyMultiTimestamps ( resamp->prev_multiTimestamps_SRC );
      resamp->prev_multiTimestamps_SRC = multiTimestamps_SRC;

    } // end: if could not reuse any previously buffered quantites

  UINT4 numFreqBinsOut = Fstats->numFreqBins;
  REAL8 dt_SRC = multiTimeSeries_SRC->data[0]->deltaT;

  // ============================== check workspace is properly allocated and initialized ===========
  if ( (resamp->ws_multiTimeSeries_SRC == NULL) || (resamp->prev_numSamples_SRC != numSamples_SRC) )
    {
      XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->ws_multiTimeSeries_SRC );
      XLAL_CHECK ( (resamp->ws_multiTimeSeries_SRC = XLALDuplicateMultiCOMPLEX8TimeSeries ( multiTimeSeries_SRC )) != NULL, XLAL_EFUNC );
      XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->ws_multiFa_SRC );
      XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->ws_multiFb_SRC );
      XLAL_CHECK ( (resamp->ws_multiFa_SRC = XLALDuplicateMultiCOMPLEX8TimeSeries ( multiTimeSeries_SRC )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( (resamp->ws_multiFb_SRC = XLALDuplicateMultiCOMPLEX8TimeSeries ( multiTimeSeries_SRC )) != NULL, XLAL_EFUNC );

      XLALDestroyCOMPLEX8Vector ( resamp->ws_outaX );
      XLALDestroyCOMPLEX8Vector ( resamp->ws_outbX );
      XLAL_CHECK ( (resamp->ws_outaX = XLALCreateCOMPLEX8Vector ( numSamples_SRC )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( (resamp->ws_outbX = XLALCreateCOMPLEX8Vector ( numSamples_SRC )) != NULL, XLAL_EFUNC );

      XLALDestroyCOMPLEX8FFTPlan ( resamp->ws_fftplan );
      XLAL_CHECK ( (resamp->ws_fftplan = XLALCreateCOMPLEX8FFTPlan ( numSamples_SRC, 1, 0) ) != NULL, XLAL_EFUNC );

      XLALFree ( resamp->ws_Fa_k );
      XLALFree ( resamp->ws_Fb_k );
      XLAL_CHECK ( (resamp->ws_Fa_k = XLALCalloc ( numFreqBinsOut, sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );
      XLAL_CHECK ( (resamp->ws_Fb_k = XLALCalloc ( numFreqBinsOut, sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );

      resamp->prev_numSamples_SRC = numSamples_SRC;
    } // if number of SRC samples has changed
  else
    {
      XLAL_CHECK ( XLALCopyMultiCOMPLEX8TimeSeries ( resamp->ws_multiTimeSeries_SRC, multiTimeSeries_SRC ) == XLAL_SUCCESS, XLAL_EFUNC );

      memset ( resamp->ws_Fa_k, 0, numFreqBinsOut * sizeof(COMPLEX8) );
      memset ( resamp->ws_Fb_k, 0, numFreqBinsOut * sizeof(COMPLEX8) );
    }
  // update number of prev samples only at the end!
  resamp->prev_numSamples_SRC = numSamples_SRC;
  // ====================================================================================================

  /* compute the fractional bin offset between the user requested initial frequency */
  /* and the closest output frequency bin */

  REAL8 diff = multiTimeSeries_DET->data[0]->f0 - thisPoint.fkdot[0]; /* the difference between the new timeseries heterodyne frequency and the user requested lowest frequency */

  // use given frequency resolution or exactly 'diff' if dFreq=0 // FIXME: temporary fix until we properly figure out 1-bin resampling efficiently
  INT4  diff_bins = (INT4)lround( diff / dFreqOut );           /* the rounded number of output frequency bins difference */
  REAL8 shift = diff - dFreqOut * diff_bins;                       /* the fractional bin frequency offset */

  /* store AM coefficient integrals in local variables */
  REAL4 Ad = multiAMcoef->Mmunu.Ad;
  REAL4 Bd = multiAMcoef->Mmunu.Bd;
  REAL4 Cd = multiAMcoef->Mmunu.Cd;
  REAL4 Ed = multiAMcoef->Mmunu.Ed;
  REAL4 Dd = multiAMcoef->Mmunu.Dd;
  REAL4 Dd_inv = 1.0f / Dd;

  /* shift the timeseries by a fraction of a frequency bin so that user requested frequency is exactly resolved */
  XLAL_CHECK ( XLALFrequencyShiftMultiCOMPLEX8TimeSeries ( resamp->ws_multiTimeSeries_SRC, shift ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* apply spin derivitive correction to resampled timeseries */
  /* this function only applies a correction if there are any non-zero spin derivitives */
  XLAL_CHECK ( XLALSpinDownCorrectionMultiTS ( resamp->ws_multiTimeSeries_SRC, &thisPoint ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALAntennaWeightMultiCOMPLEX8TimeSeries ( resamp->ws_multiFa_SRC, resamp->ws_multiFb_SRC, resamp->ws_multiTimeSeries_SRC, multiAMcoef, multiTimestamps_SRC ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* we now compute the FFTs of the resampled functions Fa and Fb for each detector */
  /* and combine them into the multi-detector F-statistic */

  /* define new initial frequency of the frequency domain representations of Fa and Fb */
  /* before the shift the zero bin was the heterodyne frequency */
  /* now we've shifted it by N - NhalfPosDC(N) bins */
  REAL8 f0_shifted = resamp->ws_multiTimeSeries_SRC->data[0]->f0 - NhalfNeg(numSamples_SRC) * dFreqOut;
  /* define number of bins offset from the internal start frequency bin to the user requested bin */
  UINT4 offset_bins = (UINT4) lround ( ( thisPoint.fkdot[0] - f0_shifted ) / dFreqOut );

  UINT4 numDetectors = resamp->multiTimeSeries_DET->length;
  /* loop over detectors */
  for ( UINT4 X=0; X < numDetectors; X++ )
    {
      const COMPLEX8Vector *ina = resamp->ws_multiFa_SRC->data[X]->data; /* we point the input to the current detector Fa timeseries */
      const COMPLEX8Vector *inb = resamp->ws_multiFb_SRC->data[X]->data; /* we point the input to the current detector Fb timeseries */

      /* Fourier transform the resampled Fa(t) and Fb(t) */
      XLAL_CHECK ( XLALCOMPLEX8VectorFFT ( resamp->ws_outaX, ina, resamp->ws_fftplan ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK ( XLALCOMPLEX8VectorFFT ( resamp->ws_outbX, inb, resamp->ws_fftplan ) == XLAL_SUCCESS, XLAL_EFUNC );

      /* the complex FFT output is shifted such that the heterodyne frequency is at DC */
      /* we need to shift the negative frequencies to before the positive ones */
      XLAL_CHECK ( XLALReorderFFTWtoSFT ( resamp->ws_outaX ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK ( XLALReorderFFTWtoSFT ( resamp->ws_outbX ) == XLAL_SUCCESS, XLAL_EFUNC );

      REAL4 AdX = multiAMcoef->data[X]->A;
      REAL4 BdX = multiAMcoef->data[X]->B;
      REAL4 CdX = multiAMcoef->data[X]->C;
      REAL4 EdX = 0; // FIXME
      REAL4 DdX_inv = 1.0 / multiAMcoef->data[X]->D;

      /* compute final Fa,Fb and Fstats (per-detector and combined) */
      for ( UINT4 k = 0; k < numFreqBinsOut; k++ )
        {
          UINT4 idy = k + offset_bins;
          COMPLEX8 FaX_k = dt_SRC * resamp->ws_outaX->data[idy];
          COMPLEX8 FbX_k = dt_SRC * resamp->ws_outbX->data[idy];

          resamp->ws_Fa_k[k] += FaX_k;
          resamp->ws_Fb_k[k] += FbX_k;

          if ( whatToCompute & FSTATQ_FAFB_PER_DET )
            {
              Fstats->FaPerDet[X][k] = FaX_k;
              Fstats->FbPerDet[X][k] = FbX_k;
            }

          if ( whatToCompute & FSTATQ_2F_PER_DET )
            {
              Fstats->twoFPerDet[X][k] = XLALComputeFstatFromFaFb ( FaX_k, FbX_k, AdX, BdX, CdX, EdX, DdX_inv );
            }
        } // for k < numFreqBinsOut

    } // for X < numDetectors

  if ( whatToCompute & FSTATQ_FAFB )
    {
      for ( UINT4 k=0; k < numFreqBinsOut; k ++ )
        {
          Fstats->Fa[k] = resamp->ws_Fa_k[k];
          Fstats->Fb[k] = resamp->ws_Fb_k[k];
        } // for k < numFreqBinsOut
    } // if FSTATQ_FAFB

  if ( whatToCompute & FSTATQ_2F )
    {
      for ( UINT4 k=0; k < numFreqBinsOut; k++ )
        {
          Fstats->twoF[k] = XLALComputeFstatFromFaFb ( resamp->ws_Fa_k[k], resamp->ws_Fb_k[k], Ad, Bd, Cd, Ed, Dd_inv );
        } // for k < numFreqBinsOut
    } // if FSTATQ_2F

  // Return F-atoms per detector
  if (whatToCompute & FSTATQ_ATOMS_PER_DET) {
    XLAL_ERROR(XLAL_EFAILED, "NOT implemented!");
  }

  Fstats->Mmunu = multiAMcoef->Mmunu;

  return XLAL_SUCCESS;

} // ComputeFstat_Resamp()


/**
 * Computed the weighted timeseries Fa(t) = x(t).a(t) and Fb(t) = x(t).b(t) for a multi-detector timeseries
 */
int
XLALAntennaWeightMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *multi_ax,                      /**< [out] the timeseries weighted by a(t) */
                                           MultiCOMPLEX8TimeSeries *multi_bx,                      /**< [out] the timeseries weighted by b(t) */
                                           const MultiCOMPLEX8TimeSeries *multiTimeSeries,         /**< [in] the input multi-detector timeseries */
                                           const MultiAMCoeffs *multiAMcoef,                       /**< [in] the multi-detector AM coefficients */
                                           const MultiLIGOTimeGPSVector *multiTimestamps           /**< [in] the multi SFT timestamps */
                                           )
{
  // input sanity checks
  XLAL_CHECK ( (multi_ax != NULL) && (multi_bx != NULL), XLAL_EINVAL );
  XLAL_CHECK ( multiTimeSeries != NULL, XLAL_EINVAL );
  XLAL_CHECK ( multiAMcoef != NULL, XLAL_EINVAL );
  XLAL_CHECK ( multiTimestamps != NULL, XLAL_EINVAL );
  UINT4 numDetectors = multiTimeSeries->length;
  XLAL_CHECK ( (numDetectors > 0) && (multiAMcoef->length == numDetectors) && (multiTimestamps->length == numDetectors), XLAL_EINVAL );
  XLAL_CHECK ( (multi_ax->length == numDetectors) && (multi_bx->length == numDetectors), XLAL_EINVAL );

  /* loop over detectors */
  for ( UINT4 X=0; X < numDetectors; X++)
    {
      /* point to current detector params */
      COMPLEX8TimeSeries *timeseries = multiTimeSeries->data[X];
      AMCoeffs *AMcoef = multiAMcoef->data[X];
      LIGOTimeGPSVector *timestamps = multiTimestamps->data[X];
      COMPLEX8TimeSeries *ax = multi_ax->data[X];
      COMPLEX8TimeSeries *bx = multi_bx->data[X];

      XLAL_CHECK ( XLALAntennaWeightCOMPLEX8TimeSeries ( ax, bx, timeseries, AMcoef, timestamps ) == XLAL_SUCCESS, XLAL_EFUNC );
    } // for X < numDetectors

  /* success */
  return XLAL_SUCCESS;

} // XLALAntennaWeightMultiCOMPLEX8TimeSeries()


/**
 * Computed the weighted timeseries Fa(t) = x(t).a(t) and Fb(t) = x(t).b(t) for a multi-detector timeseries
 */
int
XLALAntennaWeightCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *ax,                       /**< [out] the timeseries weighted by a(t) */
                                      COMPLEX8TimeSeries *bx,                       /**< [out] the timeseries weighted by b(t) */
                                      const COMPLEX8TimeSeries *timeseries,         /**< [in] the input timeseries */
                                      const AMCoeffs *AMcoef,                       /**< [in] the AM coefficients */
                                      const LIGOTimeGPSVector *timestamps           /**< [in] SFT timestamps */
                                      )
{
  // check input sanity
  XLAL_CHECK ( (ax != NULL) && ( bx != NULL ), XLAL_EINVAL );
  XLAL_CHECK ( (timeseries != NULL) && (timeseries->data != NULL) && ( timeseries->data->length > 0 ), XLAL_EINVAL );
  XLAL_CHECK ( (AMcoef != NULL) && (AMcoef->a != NULL) && (AMcoef->b != NULL), XLAL_EINVAL );
  XLAL_CHECK ( timestamps != NULL, XLAL_EINVAL );
  UINT4 numSamples = timeseries->data->length;
  XLAL_CHECK ( (ax->data->length == numSamples) && (bx->data->length == numSamples), XLAL_EINVAL );
  UINT4 numSFTs = timestamps->length;
  XLAL_CHECK ( (AMcoef->a->length == numSFTs) && (AMcoef->b->length == numSFTs), XLAL_EINVAL );

  /* local copies */
  REAL8 start = XLALGPSGetREAL8(&timeseries->epoch);
  REAL8 deltaT = timeseries->deltaT;

  REAL8 Tsft = timestamps->deltaT;
  UINT4 nbins = lround ( Tsft / deltaT );

  // set output TS to 0, to avoid garbage data in case of gaps
  memset ( ax->data->data, 0, numSamples * sizeof(ax->data->data[0]) );
  memset ( bx->data->data, 0, numSamples * sizeof(bx->data->data[0]) );

  /* loop over SFT timestamps */
  for ( UINT4 j=0; j < numSFTs; j++ )
    {
      REAL8 t = XLALGPSGetREAL8 ( &(timestamps->data[j]) );                     /* the GPS time at the start of the SFT */
      UINT4 start_index = (UINT4) lround ( (t - start) / deltaT );      /* index of timesample corresponding to the start of the SFT */
      REAL4 a = AMcoef->a->data[j];                              /* value of the antenna pattern a(t) at the MID-POINT of the SFT */
      REAL4 b = AMcoef->b->data[j];                              /* value of the antenna pattern b(t) at the MID-POINT of the SFT */

      /* loop over samples from this SFT */
      for ( UINT4 k=0; k < nbins; k++ )
        {
          UINT4 time_index = start_index + k;
          if ( time_index >= numSamples ) {
            break;
          }
          /* weight the complex timeseries by the antenna patterns */
          ax->data->data[time_index] = a * timeseries->data->data[time_index];
          bx->data->data[time_index] = b * timeseries->data->data[time_index];
        } // for k < nbins

    } // for j < numSFTs

  /* success */
  return XLAL_SUCCESS;

} // XLALAntennaWeightCOMPLEX8TimeSeries()


/**
 * Performs barycentric resampling on a multi-detector timeseries
 */
int
XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries **mTimeSeries_SRC,		///< [out] resampled timeseries in the source (SRC) frame
                                                 MultiLIGOTimeGPSVector **mTimestamps_SRC,		///< [out] multi SFT timestamps in the SRC frame
                                                 const MultiCOMPLEX8TimeSeries *mTimeSeries_DET,	///< [in] detector frame (DET) timeseries
                                                 const MultiLIGOTimeGPSVector *mTimestamps_DET,		///< [in] multi SFT timestamps in the DET frame
                                                 const MultiSSBtimes *mSRC_timing,			///< [in] multi-detector SRC timing data (time offsets+derivatives)
                                                 const REAL8 deltaF					///< [in] user defined output frequency resolution
                                                 )
{
  // check input sanity
  XLAL_CHECK ( (mTimeSeries_SRC != NULL) && ( (*mTimeSeries_SRC) == NULL ), XLAL_EINVAL );
  XLAL_CHECK ( (mTimestamps_SRC != NULL) && ((*mTimestamps_SRC) == NULL), XLAL_EINVAL );

  XLAL_CHECK ( mTimeSeries_DET != NULL, XLAL_EINVAL );
  XLAL_CHECK ( mTimestamps_DET != NULL, XLAL_EINVAL );
  XLAL_CHECK ( mSRC_timing != NULL, XLAL_EINVAL );

  XLAL_CHECK ( deltaF > 0, XLAL_EINVAL );
  UINT4 numDetectors = mTimeSeries_DET->length;
  XLAL_CHECK ( (numDetectors >0) && (mSRC_timing->length == numDetectors) && (mTimestamps_DET->length == numDetectors), XLAL_EINVAL );

  /* define the length of an SFT (assuming 1/T resolution) */
  REAL8 Tsft = mTimestamps_DET->data[0]->deltaT;

  /* find earliest and latest SRC time */
  LIGOTimeGPS earliestSRC, latestSRC;
  XLAL_CHECK ( XLALEarliestMultiSSBtime ( &earliestSRC, mSRC_timing, Tsft ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALLatestMultiSSBtime ( &latestSRC, mSRC_timing, Tsft ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* determine resampled timeseries parameters */
  REAL8 Teff = 1.0 / deltaF;                                       /* the effective observation time based on the requested frequency resolution (for zero padding) */
  REAL8 fHet = mTimeSeries_DET->data[0]->f0;                               /* the input timeseries heterodyne frequency */

  /* redefine sample rate and compute number of samples in the new timeseries */
  REAL8 deltaT = mTimeSeries_DET->data[0]->deltaT;                         /* the sample rate of the downsampled detector frame timeseries */
  UINT4 numSamplesOut = (UINT4) ceil ( Teff / deltaT );      /* we use ceil() so that we artificially widen the band rather than reduce it */
  REAL8 deltaTEff = Teff / numSamplesOut;

  // allocate memory for the output resampled timeseries
  XLAL_CHECK ( ((*mTimeSeries_SRC) = XLALMalloc ( sizeof(*(*mTimeSeries_SRC)) )) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ((*mTimeSeries_SRC)->data = XLALMalloc ( numDetectors * sizeof((*mTimeSeries_SRC)->data[0]) )) != NULL, XLAL_ENOMEM );
  (*mTimeSeries_SRC)->length = numDetectors;
  XLAL_CHECK ( ((*mTimestamps_SRC) = XLALCalloc ( 1, sizeof( *(*mTimestamps_SRC) ))) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( ((*mTimestamps_SRC)->data = XLALCalloc ( numDetectors, sizeof( (*mTimestamps_SRC)->data[0] ))) != NULL, XLAL_EFUNC );
  (*mTimestamps_SRC)->length = numDetectors;

  for ( UINT4 X=0; X < numDetectors; X++)
    {
      // shorthand pointers
      SSBtimes *SRCtimingX = mSRC_timing->data[X];
      COMPLEX8TimeSeries *TimeSeries_DETX = mTimeSeries_DET->data[X];
      LIGOTimeGPSVector *Timestamps_DETX = mTimestamps_DET->data[X];

      // create empty timeseries structures for the resampled timeseries
      XLAL_CHECK ( ((*mTimeSeries_SRC)->data[X] = XLALCreateCOMPLEX8TimeSeries ( TimeSeries_DETX->name, &earliestSRC, fHet, deltaTEff, &emptyLALUnit, numSamplesOut )) != NULL, XLAL_EFUNC );
      memset ( (*mTimeSeries_SRC)->data[X]->data->data, 0, numSamplesOut * sizeof(COMPLEX8)) ; 	// set all time-samples to zero (in case there are gaps)
      XLAL_CHECK ( ((*mTimestamps_SRC)->data[X] = XLALCreateTimestampVector ( Timestamps_DETX->length )) != NULL, XLAL_EFUNC );
      (*mTimestamps_SRC)->data[X]->deltaT = Timestamps_DETX->deltaT;

      // perform resampling on current detector timeseries */
      XLAL_CHECK ( XLALBarycentricResampleCOMPLEX8TimeSeries ( (*mTimeSeries_SRC)->data[X], (*mTimestamps_SRC)->data[X], TimeSeries_DETX, Timestamps_DETX, SRCtimingX ) == XLAL_SUCCESS, XLAL_EFUNC );
    } // for X < numDetectors

  /* success */
  return XLAL_SUCCESS;

} // XLALBarycentricResampleMultiCOMPLEX8TimeSeries()


/**
 * Performs barycentric resampling of a timeseries in the detector frame
 * into a uniformly-sampled timeseries at the source frame.
 *
 * We expect that the output timeseries has already been allocated correctly,
 * *and* carry the correct start-time epoch for the output! (FIXME!)
 *
 */
int
XLALBarycentricResampleCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *TimeSeries_SRC,		///< [out] resampled timeseries in the source (SRC) frame x(t(t_SRC)), must be alloced correctly already!
                                            LIGOTimeGPSVector *Timestamps_SRC,		///< [out] SFT timestamps translated into the source frame
                                            const COMPLEX8TimeSeries *TimeSeries_DET,	///< [in] the input detector-frame timeseries x(t)
                                            const LIGOTimeGPSVector *Timestamps_DET,	///< [in] the SFT timestamps in the detector frame
                                            const SSBtimes *SRC_timing			///< [in] the source-frame time-shifts and time-derivatives at the SFT midpoints
                                            )
{
  // check input sanity
  XLAL_CHECK ( (TimeSeries_DET != NULL) && (TimeSeries_DET->data != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (TimeSeries_SRC != NULL) && (TimeSeries_SRC->data != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (Timestamps_DET != NULL) && (Timestamps_DET->data != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (Timestamps_SRC != NULL) && (Timestamps_SRC->data != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (SRC_timing != NULL) && (SRC_timing->DeltaT != NULL) && (SRC_timing->Tdot != NULL), XLAL_EINVAL );

  UINT4 numSamples_DET = TimeSeries_DET->data->length;
  XLAL_CHECK ( numSamples_DET > 0, XLAL_EINVAL );

  UINT4 numSFTs = Timestamps_DET->length;
  XLAL_CHECK ( (numSFTs > 0) && (SRC_timing->DeltaT->length == numSFTs) && (SRC_timing->Tdot->length == numSFTs), XLAL_EINVAL );
  XLAL_CHECK ( Timestamps_SRC->length == numSFTs, XLAL_EINVAL );

  /* define some useful shorthands */
  REAL8 Tsft    = Timestamps_DET->deltaT;
  REAL8 refTime = XLALGPSGetREAL8 ( &SRC_timing->refTime );
  REAL8 fHet    = TimeSeries_DET->f0;

  REAL8 start_DET   = XLALGPSGetREAL8 ( &TimeSeries_DET->epoch );
  REAL8 deltaT_DET  = TimeSeries_DET->deltaT;
  REAL8 end_DET     = start_DET + (numSamples_DET - 1) * deltaT_DET;	// time of *last sample* in detector-frame timeseries

  REAL8 start_SRC   = XLALGPSGetREAL8 ( &(TimeSeries_SRC->epoch) );
  REAL8 deltaT_SRC  = TimeSeries_SRC->deltaT;

  UINT4 numSamples_SRC = TimeSeries_SRC->data->length;

  /* allocate memory for the uniformly sampled detector time samples (Fa and Fb real and imaginary) */
  REAL8Vector* ts[2]; // store real and imaginary parts of input timeseries as separate real vectors
  XLAL_CHECK ( (ts[0] = XLALCreateREAL8Vector ( numSamples_DET )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (ts[1] = XLALCreateREAL8Vector ( numSamples_DET )) != NULL, XLAL_EFUNC );

  /* allocate memory for the *uniform* detector time vector required for interpolation */
  REAL8Vector *t_DET; // a vector of *uniform* time values in the detector frame (for interpolation)
  XLAL_CHECK ( (t_DET = XLALCreateREAL8Vector ( numSamples_DET )) != NULL, XLAL_EFUNC );

  /* place the timeseries into REAL8Vectors for gsl to be able to interpolate them */
  for ( UINT4 j=0; j < numSamples_DET; j++ )
    {
      t_DET->data[j] = start_DET + j * deltaT_DET;
      ts[0]->data[j] = crealf ( TimeSeries_DET->data->data[j] );
      ts[1]->data[j] = cimagf ( TimeSeries_DET->data->data[j] );
    } // for j < numSamples_DET

  /* initialise the gsl spline interpolation for each of the 2 timeseries */
  gsl_spline* spline_ts[2]; XLAL_INIT_MEM(spline_ts);
  XLAL_CHECK ( XLALGSLInitInterpolateREAL8Vector ( &(spline_ts[0]), t_DET, ts[0] ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALGSLInitInterpolateREAL8Vector ( &(spline_ts[1]), t_DET, ts[1] ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* loop over SFT timestamps to compute the detector frame time samples corresponding to uniformly sampled SRC time samples */
  for ( UINT4 j=0; j < numSFTs; j++ )
    {
      /* define some useful shorthands */
      REAL8 Tdot         = SRC_timing->Tdot->data[j];                                         /* the instantaneous time derivitive dt_SRC/dt_DET at the MID-POINT of the SFT */
      REAL8 SFTmid_SRC   = refTime + SRC_timing->DeltaT->data[j];                             /* MID-POINT time of the SFT at the SRC */
      REAL8 SFTstart_SRC = SFTmid_SRC - 0.5*Tsft*Tdot;                                 /* START time of the SFT at the SRC */
      XLALGPSSetREAL8 ( &(Timestamps_SRC->data[j]), SFTstart_SRC );
      REAL8 SFTend_SRC   = SFTmid_SRC + 0.5*Tsft*Tdot;                                 /* END time of the SFT at the SRC */
      REAL8 SFTstart_DET = XLALGPSGetREAL8 ( &(Timestamps_DET->data[j]) );                   /* START time of the SFT at the detector */
      REAL8 SFTmid_DET   = SFTstart_DET + 0.5*Tsft;                                    /* MID-POINT time of the SFT at the detector */

      /* define some indices */
      UINT4 SFTidx_start_SRC  = lround ( (SFTstart_SRC - start_SRC) / deltaT_SRC );       /* the index of the resampled timeseries corresponding to the start of the SFT */
      UINT4 SFTidx_end_SRC    = lround ( (SFTend_SRC - start_SRC) / deltaT_SRC );         /* the index of the resampled timeseries corresponding to the end of the SFT */
      UINT4 SFTnumSamples_SRC = SFTidx_end_SRC - SFTidx_start_SRC + 1;                          /* the number of samples in the SRC-frame for this SFT */

      /* allocate memory for the *non-uniform* detector time samples for this SFT */
      /* have to allocate it inside the loop because it may have different lengths for each SFT */
      REAL8Vector *detectortimes; // a vector of *non-uniform* time values in the detector frame (used for interpolation) */
      XLAL_CHECK ( (detectortimes = XLALCreateREAL8Vector ( SFTnumSamples_SRC )) != NULL, XLAL_EFUNC );

      /* for each time sample in the SRC frame for this SFT we estimate the detector time. */
      /* We use a linear approximation expanding around the midpoint of an SFT where  */
      /* t_DET = SFTmid_DET + (t_SRC - SFTmid_SRC)*dt_DET/dt_SRC */
      for ( UINT4 k=0; k < SFTnumSamples_SRC; k++ )
        {
          REAL8 t_SRC = start_SRC + ( k + SFTidx_start_SRC ) * deltaT_SRC;                 /* the SRC time of the current resampled time sample */
          detectortimes->data[k] = SFTmid_DET + ( t_SRC - SFTmid_SRC ) / Tdot;          /* the approximated DET time of the current resampled time sample */

          /*
           * NOTE: we need to be careful that none of the times falls outside
           * of the range of detector timesamples, in order to avoid problems in the interpolation
           * therefore we truncate the detector-times to fully fall within the detector timeseries span
           */
          if ( detectortimes->data[k] > end_DET )
            {
              detectortimes->data[k] = end_DET;
              XLALPrintWarning ("%s: time-sample jSFT=%d, kSample=%d at t=%f to interpolate is *after* detector-timeseries, nudged back to end (end=%f)\n",
                                __func__, j, k, detectortimes->data[k], end_DET );
            }
          if ( detectortimes->data[k] < start_DET )
            {
              detectortimes->data[k] = start_DET;
              XLALPrintWarning ("%s: time-sample jSFT=%d, kSample=%d at t=%f to interpolate is *before* detector-timeseries, nudged to beginning (start=%f)\n",
                                __func__, j, k, detectortimes->data[k], start_DET );
            }
        } /* for k < SFTnumSamples_SRC */

      // interpolate on the non-uniformly sampled detector time vector for this SFT for re and im parts input timeseries
      REAL8Vector *out_ts[2]; XLAL_INIT_MEM(out_ts);
      XLAL_CHECK ( XLALGSLInterpolateREAL8Vector ( &(out_ts[0]), detectortimes, spline_ts[0] ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK ( XLALGSLInterpolateREAL8Vector ( &(out_ts[1]), detectortimes, spline_ts[1] ) == XLAL_SUCCESS, XLAL_EFUNC );

      /* place these interpolated timeseries into the output */
      /* and apply correction due to non-zero heterodyne frequency of input */
      for ( UINT4 k=0; k < SFTnumSamples_SRC; k++ )
        {
          UINT4 idx = k + SFTidx_start_SRC;                                                                     /* the full resampled timeseries index */
          if ( idx >= numSamples_SRC ) {	// temporary FIX to avoid writing outside of memory bounds (FIXME!)
            break;
          }
          REAL8 tDiff = start_SRC + idx * deltaT_SRC - detectortimes->data[k];                              /* the difference between t_SRC and t_DET */
          REAL8 cycles = fmod ( fHet * tDiff, 1 );                                                          /* the accumulated heterodyne cycles */

          /* use a look-up-table for speed to compute real and imaginary phase */
          REAL4 cosphase, sinphase;                                                                         /* the real and imaginary parts of the phase correction */
          XLAL_CHECK( XLALSinCos2PiLUT ( &sinphase, &cosphase, -cycles ) == XLAL_SUCCESS, XLAL_EFUNC );

          TimeSeries_SRC->data->data[idx] = crectf( out_ts[0]->data[k]*cosphase - out_ts[1]->data[k]*sinphase, out_ts[1]->data[k]*cosphase + out_ts[0]->data[k]*sinphase );
        } // for k < SFTnumSamples_SRC

      /* free memory used for this SFT */
      XLALDestroyREAL8Vector ( out_ts[0] );
      XLALDestroyREAL8Vector ( out_ts[1] );
      XLALDestroyREAL8Vector ( detectortimes );

    } // for j < numSFTs

  /* free memory */
  XLALDestroyREAL8Vector ( ts[0] );
  XLALDestroyREAL8Vector ( ts[1] );
  gsl_spline_free ( spline_ts[0] );
  gsl_spline_free ( spline_ts[1] );
  XLALDestroyREAL8Vector ( t_DET );

  /* success */
  return XLAL_SUCCESS;

} // XLALBarycentricResampleCOMPLEX8TimeSeries()
