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

// internal Resampling-specific parameters
struct tagFstatInput_Resamp {
  MultiCOMPLEX8TimeSeries  *multiTimeSeries_DET;	// input SFTs converted into a heterodyned timeseries

  // ----- buffering -----
  PulsarDopplerParams prev_doppler;			// buffering: previous phase-evolution ("doppler") parameters

  MultiAMCoeffs *prev_multiAMcoef;			// buffering: previous AM-coeffs, unique to skypos
  MultiSSBtimes *prev_multiSSBsky;			// buffering: previous sky-only multiSSB times, depends on skypos and reftime

  MultiCOMPLEX8TimeSeries *prev_multiFa_SRC; 		// buffering: final multi-detector SRC timeseries weighted by a(t)
  MultiCOMPLEX8TimeSeries *prev_multiFb_SRC; 		// buffering: final multi-detector SRC timeseries weighted by b(t)
};

static void
DestroyFstatInput_Resamp ( FstatInput_Resamp* resamp )
{
  XLALDestroyMultiCOMPLEX8TimeSeries (resamp->multiTimeSeries_DET );

  // ----- free buffer
  XLALDestroyMultiAMCoeffs ( resamp->prev_multiAMcoef );
  XLALDestroyMultiSSBtimes ( resamp->prev_multiSSBsky );
  XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->prev_multiFa_SRC );
  XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->prev_multiFb_SRC );

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

#if 0
  /* compute the correct SFT start-times:
   * the function XLALMultiSFTVectorToCOMPLEX8TimeSeries may have (internally) shifted SFT start times
   * to coincide with time-samples of the time series,
   * and since these times will be used later on for the resampling we also need to generate
   * correspondingly 'nudged' SFT timestamps.
   */
  // correct SFT start-times to integer bins of the timeseries:
  REAL8 t0 = XLALGPSGetREAL8 ( &(resamp->multiTimeSeries_DET->data[0]->epoch) );
  REAL8 dt = resamp->multiTimeSeries_DET->data[0]->deltaT;
  for ( UINT4 X = 0; X < common->timestamps->length; X ++ )
      {
        LIGOTimeGPSVector *ts = common->timestamps->data[X];
        for ( UINT4 j = 0; j < ts->length; j ++ )
          {
            REAL8 tj = XLALGPSGetREAL8 ( &(ts->data[j]) );
            REAL8 Dtj = tj - t0;
            REAL8 DtjNudged = round ( Dtj / dt ) * dt;
            REAL8 tjNudged = t0 + DtjNudged;
            if ( fabs ( Dtj - DtjNudged ) > 1e-9 ) {
              XLALPrintInfo ("Nudged Dtj = %.16g / %.16g = %.16g to DtjNudged = %.16g / %.16g = %.16g: Delta = %.16g\n",
                             Dtj, dt, Dtj / dt, DtjNudged, dt, DtjNudged / dt, Dtj - DtjNudged );
            }
            XLALGPSSetREAL8 ( &(ts->data[j]), tjNudged );
          } // for j < numTimeSamples
      } // for X < numDetectors
#endif

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
  UINT4 numTimeSamplesIn = multiTimeSeries_DET->data[0]->data->length;
  REAL8 dtIn = multiTimeSeries_DET->data[0]->deltaT;
  REAL8 TspanIn = numTimeSamplesIn * dtIn;
  REAL8 dFreqOut = (Fstats->dFreq > 0) ? Fstats->dFreq : 1.0 / TspanIn;

  MultiAMCoeffs *multiAMcoef;
  MultiCOMPLEX8TimeSeries *multiFa_SRC = NULL;
  MultiCOMPLEX8TimeSeries *multiFb_SRC = NULL;
  // ============================== BEGIN: handle buffering =============================
  // ----- is it the same skyposition and reference time as last call ? -----
  if ( (resamp->prev_doppler.Alpha == thisPoint.Alpha) &&
       (resamp->prev_doppler.Delta == thisPoint.Delta) &&
       (XLALGPSDiff ( &resamp->prev_doppler.refTime, &thisPoint.refTime ) == 0 )
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
        { // ----- no changes in sky + binary ==> reuse everything
          multiFa_SRC = resamp->prev_multiFa_SRC;
          multiFb_SRC = resamp->prev_multiFb_SRC;
        }
      else
        {  // ----- same skypos but changes in binary-orbital parameters: recompute just those
          MultiCOMPLEX8TimeSeries *multiTimeSeries_SRC = NULL;
          MultiLIGOTimeGPSVector *multiTimestamps_SRC = NULL;
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

          XLAL_CHECK ( XLALAntennaWeightMultiCOMPLEX8TimeSeries ( &multiFa_SRC, &multiFb_SRC, multiTimeSeries_SRC, multiAMcoef, multiTimestamps_SRC ) == XLAL_SUCCESS, XLAL_EFUNC );

          XLALDestroyMultiCOMPLEX8TimeSeries ( multiTimeSeries_SRC );
          XLALDestroyMultiTimestamps ( multiTimestamps_SRC );

          // ----- store new weighted SRC timeseries in buffer ----------
          resamp->prev_doppler = thisPoint;
          XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->prev_multiFa_SRC );
          XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->prev_multiFb_SRC );
          resamp->prev_multiFa_SRC = multiFa_SRC;
          resamp->prev_multiFb_SRC = multiFb_SRC;
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

      MultiCOMPLEX8TimeSeries *multiTimeSeries_SRC = NULL;
      MultiLIGOTimeGPSVector *multiTimestamps_SRC = NULL;

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

      // antenna-weighting
      XLAL_CHECK ( (multiAMcoef = XLALComputeMultiAMCoeffs ( multiDetStates, multiWeights, skypos )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( XLALAntennaWeightMultiCOMPLEX8TimeSeries ( &multiFa_SRC, &multiFb_SRC, multiTimeSeries_SRC, multiAMcoef, multiTimestamps_SRC ) == XLAL_SUCCESS, XLAL_EFUNC );

      XLALDestroyMultiCOMPLEX8TimeSeries ( multiTimeSeries_SRC );
      XLALDestroyMultiTimestamps ( multiTimestamps_SRC );

      // ----- store everything in buffer ----------
      XLALDestroyMultiAMCoeffs ( resamp->prev_multiAMcoef );
      XLALDestroyMultiSSBtimes ( resamp->prev_multiSSBsky );
      XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->prev_multiFa_SRC );
      XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->prev_multiFb_SRC );

      resamp->prev_doppler = thisPoint;
      resamp->prev_multiAMcoef = multiAMcoef;
      resamp->prev_multiSSBsky = multiSSBsky;
      resamp->prev_multiFa_SRC = multiFa_SRC;
      resamp->prev_multiFb_SRC = multiFb_SRC;

    } // end: if could not reuse any previously buffered quantites
  // ============================== END: handle buffering =============================


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

  // *copy* complete resampled multi-complex8 timeseries so we can apply spindown-corrections to it
  MultiCOMPLEX8TimeSeries *multiFa_spin, *multiFb_spin;
  XLAL_CHECK ( (multiFa_spin = XLALDuplicateMultiCOMPLEX8TimeSeries ( multiFa_SRC )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (multiFb_spin = XLALDuplicateMultiCOMPLEX8TimeSeries ( multiFb_SRC )) != NULL, XLAL_EFUNC );

  /* shift the timeseries by a fraction of a frequency bin so that user requested frequency is exactly resolved */
  if (shift != 0.0)
    {
      XLAL_CHECK ( XLALFrequencyShiftMultiCOMPLEX8TimeSeries ( &multiFa_spin, shift ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK ( XLALFrequencyShiftMultiCOMPLEX8TimeSeries ( &multiFb_spin, shift ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

  /* apply spin derivitive correction to resampled timeseries */
  /* this function only applies a correction if there are any non-zero spin derivitives */
  XLAL_CHECK ( XLALSpinDownCorrectionMultiFaFb ( &multiFa_spin, &multiFb_spin, &thisPoint ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* we now compute the FFTs of the resampled functions Fa and Fb for each detector */
  /* and combine them into the multi-detector F-statistic */

  /* we use the first detector Fa time series to obtain the number of time samples and the sampling time */
  /* these should be the same for all Fa and Fb timeseries */
  UINT4 numSamples = multiFa_spin->data[0]->data->length;
  REAL8 dt = multiFa_spin->data[0]->deltaT;

  /* allocate memory for individual-detector FFT outputs */
  COMPLEX8Vector *outaX, *outbX;
  XLAL_CHECK ( (outaX = XLALCreateCOMPLEX8Vector(numSamples)) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (outbX = XLALCreateCOMPLEX8Vector(numSamples)) != NULL, XLAL_EFUNC );

  /* make forwards FFT plan - this will be re-used for each detector */
  ComplexFFTPlan *pfwd;
  XLAL_CHECK ( (pfwd = XLALCreateCOMPLEX8FFTPlan ( numSamples, 1, 0) ) != NULL, XLAL_EFUNC );

  UINT4 numFreqBins = Fstats->numFreqBins;

  /* define new initial frequency of the frequency domain representations of Fa and Fb */
  /* before the shift the zero bin was the heterodyne frequency */
  /* now we've shifted it by N - NhalfPosDC(N) bins */
  REAL8 f0_shifted = multiFa_spin->data[0]->f0 - NhalfNeg(numSamples) * dFreqOut;
  /* define number of bins offset from the internal start frequency bin to the user requested bin */
  UINT4 offset_bins = (UINT4) lround ( ( thisPoint.fkdot[0] - f0_shifted ) / dFreqOut );

  COMPLEX8 *Fa_k, *Fb_k;
  XLAL_CHECK ( (Fa_k = XLALCalloc ( numFreqBins, sizeof(*Fa_k))) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( (Fb_k = XLALCalloc ( numFreqBins, sizeof(*Fa_k))) != NULL, XLAL_ENOMEM );

  UINT4 numDetectors = resamp->multiTimeSeries_DET->length;
  /* loop over detectors */
  for ( UINT4 X=0; X < numDetectors; X++ )
    {
      COMPLEX8Vector *ina = multiFa_spin->data[X]->data; /* we point the input to the current detector Fa timeseries */
      COMPLEX8Vector *inb = multiFb_spin->data[X]->data; /* we point the input to the current detector Fb timeseries */

      /* Fourier transform the resampled Fa(t) and Fb(t) */
      XLAL_CHECK ( XLALCOMPLEX8VectorFFT ( outaX, ina, pfwd ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK ( XLALCOMPLEX8VectorFFT ( outbX, inb, pfwd ) == XLAL_SUCCESS, XLAL_EFUNC );

      /* the complex FFT output is shifted such that the heterodyne frequency is at DC */
      /* we need to shift the negative frequencies to before the positive ones */
      XLAL_CHECK ( XLALFFTShiftCOMPLEX8Vector ( &outaX ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK ( XLALFFTShiftCOMPLEX8Vector ( &outbX ) == XLAL_SUCCESS, XLAL_EFUNC );

      REAL4 AdX = multiAMcoef->data[X]->A;
      REAL4 BdX = multiAMcoef->data[X]->B;
      REAL4 CdX = multiAMcoef->data[X]->C;
      REAL4 EdX = 0; // FIXME
      REAL4 DdX_inv = 1.0 / multiAMcoef->data[X]->D;

      /* compute final Fa,Fb and Fstats (per-detector and combined) */
      for ( UINT4 k = 0; k < numFreqBins; k++ )
        {
          UINT4 idy = k + offset_bins;
          COMPLEX8 FaX_k = dt * outaX->data[idy];
          COMPLEX8 FbX_k = dt * outbX->data[idy];

          Fa_k[k] += FaX_k;
          Fb_k[k] += FbX_k;

          if ( whatToCompute & FSTATQ_FAFB_PER_DET )
            {
              Fstats->FaPerDet[X][k] = FaX_k;
              Fstats->FbPerDet[X][k] = FbX_k;
            }

          if ( whatToCompute & FSTATQ_2F_PER_DET )
            {
              Fstats->twoFPerDet[X][k] = XLALComputeFstatFromFaFb ( FaX_k, FbX_k, AdX, BdX, CdX, EdX, DdX_inv );
            }
        } // for k < numFreqBins

    } // for X < numDetectors

  if ( whatToCompute & FSTATQ_FAFB )
    {
      for ( UINT4 k=0; k < numFreqBins; k ++ )
        {
          Fstats->Fa[k] = Fa_k[k];
          Fstats->Fb[k] = Fb_k[k];
        } // for k < numFreqBins
    } // if FSTATQ_FAFB

  if ( whatToCompute & FSTATQ_2F )
    {
      for ( UINT4 k=0; k < numFreqBins; k++ )
        {
          Fstats->twoF[k] = XLALComputeFstatFromFaFb ( Fa_k[k], Fb_k[k], Ad, Bd, Cd, Ed, Dd_inv );
        } // for k < numFreqBins
    } // if FSTATQ_2F

  // free memory not stored in the buffer
  XLALFree ( Fa_k );
  XLALFree ( Fb_k );
  XLALDestroyCOMPLEX8Vector ( outaX );
  XLALDestroyCOMPLEX8Vector ( outbX );
  XLALDestroyCOMPLEX8FFTPlan ( pfwd );

  XLALDestroyMultiCOMPLEX8TimeSeries ( multiFa_spin );
  XLALDestroyMultiCOMPLEX8TimeSeries ( multiFb_spin );

  // Return F-atoms per detector
  if (whatToCompute & FSTATQ_ATOMS_PER_DET) {
    XLAL_ERROR(XLAL_EFAILED, "NOT implemented!");
  }

  Fstats->Mmunu = multiAMcoef->Mmunu;

  return XLAL_SUCCESS;

} // ComputeFstat_Resamp()
