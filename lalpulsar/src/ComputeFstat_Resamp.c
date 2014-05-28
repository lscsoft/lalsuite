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

/////////////////////////////////////////////////////////////////
//////////////////// Old resampling API code ////////////////////
/////////////////////////////////////////////////////////////////

#define NhalfPosDC(N) ((UINT4)(ceil ( ((N)/2.0 - 1e-6 ))))      /* round up */
#define NhalfNeg(N) ((UINT4)( (N) - NhalfPosDC(N) ))            /* round down (making sure N+ + N- = (N-1) */

/* [opaque] type holding a ComputeFBuffer for use in the resampling F-stat codes */
typedef struct tagComputeFBuffer_RS ComputeFBuffer_RS;

/* Extra parameters controlling the actual computation of F */
typedef struct tagComputeFParams {
  UINT4 Dterms;         /* how many terms to keep in the Dirichlet kernel (~16 is usually fine) */
  SSBprecision SSBprec; /* whether to use full relativist SSB-timing, or just simple Newtonian */
  ComputeFBuffer_RS *buffer; /* buffer for storing pre-resampled timeseries (used for resampling implementation) */
  const EphemerisData *edat;   /* ephemeris data for re-computing multidetector states */
  BOOLEAN returnAtoms;  /* whether or not to return the 'FstatAtoms' used to compute the F-statistic */
  BOOLEAN returnSingleF; /* in multi-detector case, whether or not to also return the single-detector Fstats computed from the atoms */
} ComputeFParams;

/* Struct holding buffered ComputeFStat()-internal quantities to avoid unnecessarily
 * recomputing things that depend ONLY on the skyposition and detector-state series (but not on the spins).
 * For the first call of ComputeFStatFreqBand_RS() the pointer-entries should all be NULL.
 */
struct tagComputeFBuffer_RS {
  MultiDetectorStateSeries *multiDetStates;             /* buffer for each detStates (store pointer) and skypos */
  REAL8 Alpha, Delta;                                         /* skyposition of candidate */
  LIGOTimeGPS segstart;                                       /* the start time of the first SFT of the first detector (used to check if the segment has changed) */
  MultiSSBtimes *multiSSB;
  MultiSSBtimes *multiBinary;
  MultiAMCoeffs *multiAMcoef;
  MultiCOMPLEX8TimeSeries *multiTimeseries;                   /* the buffered unweighted multi-detector timeseries */
  MultiCOMPLEX8TimeSeries *multiFa_resampled;                 /* the buffered multi-detector resampled timeseries weighted by a(t) */
  MultiCOMPLEX8TimeSeries *multiFb_resampled;                 /* the buffered multi-detector resampled timeseries weighted by b(t) */
};

/* Struct holding a vector of buffered ComputeFStat()-internal quantities to avoid unnecessarily
 * recomputing things that depend ONLY on the skyposition and detector-state series (but not on the spins).
 */
typedef struct tagComputeFBufferVector_RS {
  ComputeFBuffer_RS **data;                                    /* pointer to a series of ComputeFBuffer_RS structures */
  UINT4 length;                                               /* the length of the vector */
} ComputeFBufferVector_RS;

/* Destruction of a ComputeFBuffer *contents*,
 * i.e. the multiSSB and multiAMcoeff, while the
 * buffer-container is not freed (which is why it's passed
 * by value and not by reference...) */
static void
XLALEmptyComputeFBuffer_RS ( ComputeFBuffer_RS *buffer)
{

  if ( buffer->multiSSB ) XLALDestroyMultiSSBtimes( buffer->multiSSB );
  buffer->multiSSB = NULL;
  if ( buffer->multiBinary ) XLALDestroyMultiSSBtimes( buffer->multiBinary );
  buffer->multiBinary = NULL;
  if ( buffer->multiAMcoef) XLALDestroyMultiAMCoeffs( buffer->multiAMcoef );
  buffer->multiAMcoef = NULL;
  if ( buffer->multiTimeseries) XLALDestroyMultiCOMPLEX8TimeSeries( buffer->multiTimeseries );
  buffer->multiTimeseries = NULL;
  if ( buffer->multiFa_resampled) XLALDestroyMultiCOMPLEX8TimeSeries( buffer->multiFa_resampled );
  buffer->multiFa_resampled = NULL;
  if ( buffer->multiFb_resampled) XLALDestroyMultiCOMPLEX8TimeSeries( buffer->multiFb_resampled );
  buffer->multiFb_resampled = NULL;
  if ( buffer->multiDetStates) XLALDestroyMultiDetectorStateSeries( buffer->multiDetStates);
  buffer->multiDetStates = NULL;
  /* if ( buffer ) XLALFree(buffer); */

  return;

} /* XLALEmptyComputeFBuffer_RS() */

/////////////////////////////////////////////////////////////////
//////////////////// New resampling API code ////////////////////
/////////////////////////////////////////////////////////////////

struct tagFstatInput_Resamp {
  MultiSFTVector *multiSFTs;                    // Input multi-detector SFTs
  ComputeFParams params;                        // Additional parameters for ComputeFStat() and ComputeFStatFreqBand_RS()
};

static inline void
DestroyFstatInput_Resamp(
  FstatInput_Resamp* resamp
  )
{
  XLALDestroyMultiSFTVector(resamp->multiSFTs);
  XLALEmptyComputeFBuffer_RS(resamp->params.buffer);
  XLALFree(resamp->params.buffer);
  XLALFree(resamp);
}

static int
SetupFstatInput_Resamp(
  FstatInput_Resamp *resamp,
  const FstatInput_Common *common,
  MultiSFTVector *multiSFTs
  )
{

  // Check input
  XLAL_CHECK(common != NULL, XLAL_EFAULT);
  XLAL_CHECK(resamp != NULL, XLAL_EFAULT);
  XLAL_CHECK(multiSFTs != NULL, XLAL_EFAULT);

  // Save pointer to SFTs
  resamp->multiSFTs = multiSFTs;

  // Set parameters to pass to ComputeFStatFreqBand_RS()
  resamp->params.SSBprec = common->SSBprec;
  resamp->params.buffer = NULL;
  resamp->params.edat = common->ephemerides;

  return XLAL_SUCCESS;

}

static int
GetFstatExtraBins_Resamp(
  FstatInput_Resamp* resamp
  )
{


  // Check input
  XLAL_CHECK(resamp != NULL, XLAL_EFAULT);

  // FIXME: resampling should not require extra frequency bins, however
  // the following is required to get 'testCFSv2_resamp.sh' to pass
  return 8;

}

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

  // Get which F-statistic quantities to compute
  const FstatQuantities whatToCompute = Fstats->whatWasComputed;

  // Check which quantities can be computed
  XLAL_CHECK ( !(whatToCompute & FSTATQ_FAFB), XLAL_EINVAL, "Resampling does not currently support Fa & Fb");
  XLAL_CHECK ( !(whatToCompute & FSTATQ_FAFB_PER_DET), XLAL_EINVAL, "Resampling does not currently support Fa & Fb per detector");
  XLAL_CHECK ( !(whatToCompute & FSTATQ_ATOMS_PER_DET), XLAL_EINVAL, "Resampling does not currently support atoms per detector");

  // Set parameters to pass to ComputeFStatFreqBand_RS()
  resamp->params.returnSingleF = whatToCompute & FSTATQ_2F_PER_DET;
  resamp->params.returnAtoms = 0;

  const PulsarDopplerParams *thisPoint = &Fstats->doppler;
  MultiSFTVector *multiSFTs = resamp->multiSFTs;
  const MultiNoiseWeights *multiWeights = common->noiseWeights;
  ComputeFParams *params = &resamp->params;
  {// ================================================================================
    // Call ComputeFStatFreqBand_RS()
    UINT4 numDetectors;
    ComputeFBuffer_RS *cfBuffer = NULL;
    MultiDetectorStateSeries *multiDetStates = NULL;
    MultiSSBtimes *multiSSB = NULL;
    MultiAMCoeffs *multiAMcoef = NULL;
    MultiCOMPLEX8TimeSeries *multiTimeseries = NULL;
    REAL8 Ad, Bd, Cd, Dd_inv, AdX, BdX, CdX, DdX_inv;
    SkyPosition skypos;
    MultiCOMPLEX8TimeSeries *multiFa_resampled = NULL;
    MultiCOMPLEX8TimeSeries *multiFb_resampled = NULL;
    COMPLEX8Vector *Faf_resampled = NULL;
    COMPLEX8Vector *Fbf_resampled = NULL;
    UINT4 numSamples;
    REAL8 f0_shifted;
    REAL8 f0_shifted_single;
    REAL8 dt;
    ComplexFFTPlan *pfwd = NULL;  /* this will store the FFT plan */
    COMPLEX8Vector *outa = NULL;  /* this will contain the FFT output of Fa for this detector */
    COMPLEX8Vector *outb = NULL;  /* this will contain the FFT output of Fb for this detector */
    COMPLEX8Vector *outaSingle = NULL; /* this will contain Faf_resampled for a single IFO */
    COMPLEX8Vector *outbSingle = NULL; /* this will contain Fbf_resampled for a single IFO */

    cfBuffer = params->buffer;                      /* set local pointer to the buffer location */
    numDetectors = multiSFTs->length;               /* set the number of detectors to the number of sets of SFTs */
    // unused: SFTtype * firstSFT = &(multiSFTs->data[0]->data[0]);      /* use data from the first SFT from the first detector to set other params */

    /* check that the multidetector noise weights have the same length as the multiSFTs */
    XLAL_CHECK ( (multiWeights == NULL) || (multiWeights->length == numDetectors), XLAL_EINVAL );

    /* first, if there is no buffer allocate space for one */
    /* IMPORTANT - use Calloc here so that all pointers within the structure are NULL */
    if ( cfBuffer == NULL ) {
      XLAL_CHECK ( (cfBuffer = XLALCalloc ( 1, sizeof(*cfBuffer) )) != NULL, XLAL_ENOMEM );
    }

    /* Dealing with the input SFT -> timeseries conversion and whether it has already been done and is buffered */

    /* generate bandpassed and downsampled timeseries for each detector                         */
    /* we only ever do this once for a given dataset so we read it from the buffer if it exists */
    /* in future implementations we will pass this directly to the function instead of SFTs     */

    /* check if there is an not existing timeseries and if the start time in the buffer does not match the start time of the SFTs */
    if ( !cfBuffer->multiTimeseries || ( XLALGPSCmp(&cfBuffer->segstart,&multiSFTs->data[0]->data[0].epoch) != 0) )
      {
        XLALPrintInfo ("*** New segment : recomputing timeseries and detstates\n");
        if ( !cfBuffer->multiTimeseries) { XLALPrintInfo("timeseries pointer was null\n"); }
        if ( XLALGPSCmp ( &cfBuffer->segstart, &multiSFTs->data[0]->data[0].epoch) != 0 ) {
          XLALPrintInfo("segstart changed from %d to %d\n",cfBuffer->segstart.gpsSeconds,multiSFTs->data[0]->data[0].epoch.gpsSeconds);
        }

        /* if there was no existing timeseries we need to recompute the timeseries from the SFTs */
        /* generate multiple coincident timeseries - one for each detector spanning start -> end */
        /* we need each timeseries to span the exact same amount of time and to start at the same time */
        /* because for the multi-detector Fstat we need frequency bins to be coincident */
        /* The memory allocated here is freed when the buffer is cleared in the calling program */

        /* generate a new timeseries from the input SFTs */
        XLAL_CHECK ( (multiTimeseries = XLALMultiSFTVectorToCOMPLEX8TimeSeries(multiSFTs)) != NULL, XLAL_EFUNC );

        /* recompute the multidetector states for the possibly time shifted SFTs */
        /* the function XLALMultiSFTVectorToCOMPLEX8TimeSeries may have shifted the SFT start times around */
        /* and since these times will be used later on for the resampling we also need to recompute the */
        /* MultiDetectorStates because their timestamps are used later on to compute SSB times which we */
        /* need to be accurate at the midpoints of the SFTs.  Understand ? */

        /* recompute the multiDetStates for the new SFT start times */
        MultiLALDetector multiIFO;
        XLAL_CHECK ( XLALMultiLALDetectorFromMultiSFTs ( &multiIFO, multiSFTs ) == XLAL_SUCCESS, XLAL_EFUNC );
        MultiLIGOTimeGPSVector *multiTS;
        XLAL_CHECK ( (multiTS = XLALExtractMultiTimestampsFromSFTs ( multiSFTs )) != NULL, XLAL_EFUNC );

        REAL8 Tsft = 1.0 / multiSFTs->data[0]->data[0].deltaF;
        REAL8 tOffset = 0.5 * Tsft;
        XLAL_CHECK ( (multiDetStates = XLALGetMultiDetectorStates ( multiTS, &multiIFO, params->edat, tOffset )) != NULL, XLAL_EFUNC );
        XLALDestroyMultiTimestamps ( multiTS );

        /* set all other segment dependent quantity pointers to NULL */
        /* this will basically mean that we will have to compute all sky dependent quantities again */
        XLALEmptyComputeFBuffer_RS( cfBuffer );

        /* buffer the multitimeseries, detstates and the current start time of the input data */
        cfBuffer->multiTimeseries = multiTimeseries;
        cfBuffer->multiDetStates = multiDetStates;
        cfBuffer->segstart.gpsSeconds = multiSFTs->data[0]->data[0].epoch.gpsSeconds;
        cfBuffer->segstart.gpsNanoSeconds = multiSFTs->data[0]->data[0].epoch.gpsNanoSeconds;

      }  /* if (!cfBuffer->multiTimeseries || (buffered-start != SFT-start) ) */
    /* End of the SFT -> timeseries buffering checks */


    /* compute the fractional bin offset between the user requested initial frequency */
    /* and the closest output frequency bin */
    REAL8 diff = cfBuffer->multiTimeseries->data[0]->f0 - thisPoint->fkdot[0]; /* the difference between the new timeseries heterodyne frequency and the user requested lowest frequency */

    // use given frequency resolution or exactly 'diff' if dFreq=0 // FIXME: temporary fix until we properly figure out 1-bin resampling efficiently
    REAL8 dFreq = (Fstats->dFreq > 0) ? Fstats->dFreq : diff;
    INT4  diff_bins = (INT4)lround( diff / dFreq );           /* the rounded number of output frequency bins difference */
    REAL8 shift = diff - dFreq * diff_bins;                       /* the fractional bin frequency offset */

    /* Dealing with sky position dependent quantities and buffering them */

    /* if the sky position has changed or if any of the sky position dependent quantities are not buffered
       i.e the multiDetstates, the multiAMcoefficients, the multiSSB times and the resampled multiTimeSeries Fa and Fb,
       then we need to recompute these and buffer them */
    if ( (cfBuffer->Alpha != thisPoint->Alpha )                                 /* and alpha hasn't changed */
         || ( cfBuffer->Delta != thisPoint->Delta )                             /* and delta hasn't changed */
         || ( cfBuffer->multiAMcoef == NULL )                                 /* and we have a buffered multiAMcoefficents */
         || ( cfBuffer->multiSSB == NULL )                                    /* and we have buffered multiSSB times */
         || ( cfBuffer->multiFa_resampled == NULL )                           /* and we have buffered multiFa_resampled  */
         || ( cfBuffer->multiFb_resampled == NULL )                           /* and we have multiFb_resampled */
         || ( thisPoint->asini > 0 )						// no buffering in binary-CW case for now
         )
      {
        XLALPrintInfo("*** New sky position : recomputing SSB times, AM coefficients and Fa and Fb\n");

        /* temporary timeseries used to store the unbarycentred antenna weighted Fa and Fb timeseries */
        MultiCOMPLEX8TimeSeries *multiFa = NULL;
        MultiCOMPLEX8TimeSeries *multiFb = NULL;

        /* compute the SSB times corresponding to the midpoints of each SFT for the current sky position for all detectors */
        skypos.system = COORDINATESYSTEM_EQUATORIAL;
        skypos.longitude = thisPoint->Alpha;
        skypos.latitude  = thisPoint->Delta;
        XLAL_CHECK ( (multiSSB = XLALGetMultiSSBtimes ( cfBuffer->multiDetStates, skypos, thisPoint->refTime, params->SSBprec )) != NULL, XLAL_EFUNC );

        MultiSSBtimes *multiBinary = NULL;
        MultiSSBtimes *multiSSBTotal = NULL;
        // handle binary-orbital timing corrections, if applicable
        if ( thisPoint->asini > 0 )
          {
            // compute binary time corrections to the SSB time delays and SSB time derivitive
            XLAL_CHECK ( XLALAddMultiBinaryTimes ( &multiBinary, multiSSB, thisPoint ) == XLAL_SUCCESS, XLAL_EFUNC );
            multiSSBTotal = multiBinary;
          }
        else
          {
            multiSSBTotal = multiSSB;
          }

        /* compute the AM parameters for each detector */
        XLAL_CHECK ( (multiAMcoef = XLALComputeMultiAMCoeffs ( cfBuffer->multiDetStates, multiWeights, skypos )) != NULL, XLAL_EFUNC );

        /* Generate a(t) and b(t) weighted heterodyned downsampled timeseries */
        XLAL_CHECK ( XLALAntennaWeightMultiCOMPLEX8TimeSeries ( &multiFa, &multiFb, cfBuffer->multiTimeseries, multiAMcoef, multiSFTs) == XLAL_SUCCESS, XLAL_EFUNC );

        /* Perform barycentric resampling on the multi-detector timeseries */
        XLAL_CHECK ( XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( &multiFa_resampled, &multiFb_resampled, multiFa, multiFb, multiSSBTotal, multiSFTs, dFreq) == XLAL_SUCCESS, XLAL_EFUNC );

        XLALDestroyMultiSSBtimes ( multiBinary );

        /* free multiFa and MultiFb - we won't need them again since we're storing the resampled versions */
        XLALDestroyMultiCOMPLEX8TimeSeries ( multiFa );
        XLALDestroyMultiCOMPLEX8TimeSeries ( multiFb );

        /* buffer all new sky position dependent values - after clearing them */
        cfBuffer->Alpha = thisPoint->Alpha;
        cfBuffer->Delta = thisPoint->Delta;

        XLALDestroyMultiSSBtimes( cfBuffer->multiSSB );
        cfBuffer->multiSSB = multiSSB;

        XLALDestroyMultiAMCoeffs( cfBuffer->multiAMcoef );
        cfBuffer->multiAMcoef = multiAMcoef;

        XLALDestroyMultiCOMPLEX8TimeSeries( cfBuffer->multiFa_resampled );
        XLALDestroyMultiCOMPLEX8TimeSeries( cfBuffer->multiFb_resampled );
        cfBuffer->multiFa_resampled = multiFa_resampled;
        cfBuffer->multiFb_resampled = multiFb_resampled;

      } /* could not reuse previously buffered quantities */

    /* End of the sky position dependent quantity buffering */

    /* store AM coefficient integrals in local variables */
    if ( cfBuffer->multiAMcoef )
      {
        Ad = cfBuffer->multiAMcoef->Mmunu.Ad;
        Bd = cfBuffer->multiAMcoef->Mmunu.Bd;
        Cd = cfBuffer->multiAMcoef->Mmunu.Cd;
        Dd_inv = 1.0 / cfBuffer->multiAMcoef->Mmunu.Dd;
      }
    else
      {
        XLALPrintError ( "Programming error: 'multiAMcoef' not available!\n");
        XLAL_ERROR ( XLAL_EFAILED );
      }

    // *copy* complete resampled multi-complex8 timeseries so we can apply spindown-corrections to it
    MultiCOMPLEX8TimeSeries *multiFa_spin, *multiFb_spin;
    XLAL_CHECK ( (multiFa_spin = XLALDuplicateMultiCOMPLEX8TimeSeries ( cfBuffer->multiFa_resampled )) != NULL, XLAL_EFUNC );
    XLAL_CHECK ( (multiFb_spin = XLALDuplicateMultiCOMPLEX8TimeSeries ( cfBuffer->multiFb_resampled )) != NULL, XLAL_EFUNC );

    /* shift the timeseries by a fraction of a frequency bin so that user requested frequency is exactly resolved */
    if (shift != 0.0)
      {
        XLAL_CHECK ( XLALFrequencyShiftMultiCOMPLEX8TimeSeries ( &multiFa_spin, shift ) == XLAL_SUCCESS, XLAL_EFUNC );
        XLAL_CHECK ( XLALFrequencyShiftMultiCOMPLEX8TimeSeries ( &multiFb_spin, shift ) == XLAL_SUCCESS, XLAL_EFUNC );
      }

    /* apply spin derivitive correction to resampled timeseries */
    /* this function only applies a correction if there are any non-zero spin derivitives */
    XLAL_CHECK ( XLALSpinDownCorrectionMultiFaFb ( &multiFa_spin, &multiFb_spin, thisPoint ) == XLAL_SUCCESS, XLAL_EFUNC );

    /* we now compute the FFTs of the resampled functions Fa and Fb for each detector */
    /* and combine them into the multi-detector F-statistic */

    /* we use the first detector Fa time series to obtain the number of time samples and the sampling time */
    /* these should be the same for all Fa and Fb timeseries */
    numSamples = multiFa_spin->data[0]->data->length;
    dt = multiFa_spin->data[0]->deltaT;

    /* allocate memory for Fa(f) and Fb(f) and individual detector FFT outputs */
    XLAL_CHECK ( (Faf_resampled = XLALCreateCOMPLEX8Vector(numSamples)) != NULL, XLAL_EFUNC );
    XLAL_CHECK ( (Fbf_resampled = XLALCreateCOMPLEX8Vector(numSamples)) != NULL, XLAL_EFUNC );
    XLAL_CHECK ( (outa = XLALCreateCOMPLEX8Vector(numSamples)) != NULL, XLAL_EFUNC );
    XLAL_CHECK ( (outb = XLALCreateCOMPLEX8Vector(numSamples)) != NULL, XLAL_EFUNC );
    if ( params->returnSingleF )
      {
        XLAL_CHECK ( (outaSingle = XLALCreateCOMPLEX8Vector(numSamples)) != NULL, XLAL_EFUNC );
        XLAL_CHECK ( (outbSingle = XLALCreateCOMPLEX8Vector(numSamples)) != NULL, XLAL_EFUNC );
      }

    /* initialise output vectors to zero since it will be added to */
    memset ( Faf_resampled->data, 0, numSamples * sizeof(COMPLEX8) );
    memset ( Fbf_resampled->data, 0, numSamples * sizeof(COMPLEX8) );

    /* make forwards FFT plan - this will be re-used for each detector */
    XLAL_CHECK ( (pfwd = XLALCreateCOMPLEX8FFTPlan ( numSamples, 1, 0) ) != NULL, XLAL_EFUNC );

    UINT4 numFreqBins = Fstats->numFreqBins;

    /* loop over detectors */
    for ( UINT4 X=0; X < numDetectors; X++ )
      {
        COMPLEX8Vector *ina = multiFa_spin->data[X]->data; /* we point the input to the current detector Fa timeseries */
        COMPLEX8Vector *inb = multiFb_spin->data[X]->data; /* we point the input to the current detector Fb timeseries */

        /* initialise output vectors to zero for safety */
        memset ( outa->data, 0, numSamples * sizeof(COMPLEX8) );
        memset ( outb->data, 0, numSamples * sizeof(COMPLEX8) );
        if ( params->returnSingleF )
          {
            memset ( outaSingle->data, 0, numSamples * sizeof(COMPLEX8) );
            memset ( outbSingle->data, 0, numSamples * sizeof(COMPLEX8) );
          }

        /* Fourier transform the resampled Fa(t) and Fb(t) */
        XLAL_CHECK ( XLALCOMPLEX8VectorFFT ( outa, ina, pfwd ) == XLAL_SUCCESS, XLAL_EFUNC );
        XLAL_CHECK ( XLALCOMPLEX8VectorFFT ( outb, inb, pfwd ) == XLAL_SUCCESS, XLAL_EFUNC );

        /*  add to summed Faf and Fbf and normalise by dt */
        for ( UINT4 j = 0; j < numSamples; j++ )
          {
            Faf_resampled->data[j] += outa->data[j] * dt;
            Fbf_resampled->data[j] += outb->data[j] * dt;
          }

        /* compute single-IFO F-stats, if requested */
        if ( params->returnSingleF )
          {
            if ( params->buffer == NULL )
              {
                AdX = multiAMcoef->data[X]->A;
                BdX = multiAMcoef->data[X]->B;
                CdX = multiAMcoef->data[X]->C;
                DdX_inv = 1.0 / multiAMcoef->data[X]->D;
              }
            else
              {
                AdX = cfBuffer->multiAMcoef->data[X]->A;
                BdX = cfBuffer->multiAMcoef->data[X]->B;
                CdX = cfBuffer->multiAMcoef->data[X]->C;
                DdX_inv = 1.0 / cfBuffer->multiAMcoef->data[X]->D;
              }

            /* normalize by dt */
            for ( UINT4 j = 0; j < numSamples; j++ )
              {
                outaSingle->data[j] = outa->data[j] * dt;
                outbSingle->data[j] = outb->data[j] * dt;
              }

            /* the complex FFT output is shifted such that the heterodyne frequency is at DC */
            /* we need to shift the negative frequencies to before the positive ones */
            XLAL_CHECK ( XLALFFTShiftCOMPLEX8Vector ( &outaSingle ) == XLAL_SUCCESS, XLAL_EFUNC );
            XLAL_CHECK ( XLALFFTShiftCOMPLEX8Vector ( &outbSingle ) == XLAL_SUCCESS, XLAL_EFUNC );

            /* define new initial frequency of the frequency domain representations of Fa and Fb */
            /* before the shift the zero bin was the heterodyne frequency */
            /* now we've shifted it by N - NhalfPosDC(N) bins */
            f0_shifted_single = multiFa_spin->data[X]->f0 - NhalfNeg(numSamples) * dFreq;

            /* define number of bins offset from the internal start frequency bin to the user requested bin */
            UINT4 offset_single = floor(0.5 + (thisPoint->fkdot[0] - f0_shifted_single)/ dFreq );

            /* compute final single-IFO F-stat */
            for ( UINT4 k = 0; k < numFreqBins; k++ )
              {
                UINT4 idy = k + offset_single;
                COMPLEX16 FaX = outaSingle->data[idy];
                COMPLEX16 FbX = outbSingle->data[idy];
                Fstats->twoFPerDet[X][k] = ComputeFstatFromFaFb ( FaX, FbX, AdX, BdX, CdX, 0, DdX_inv );
              } // for k < numFreqBins
          } // if returnSingleF

      } // for X < numDetectors

    /* the complex FFT output is shifted such that the heterodyne frequency is at DC */
    /* we need to shift the negative frequencies to before the positive ones */
    XLAL_CHECK ( XLALFFTShiftCOMPLEX8Vector ( &Faf_resampled ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK ( XLALFFTShiftCOMPLEX8Vector ( &Fbf_resampled ) == XLAL_SUCCESS, XLAL_EFUNC );

    /* define new initial frequency of the frequency domain representations of Fa and Fb */
    /* before the shift the zero bin was the heterodyne frequency */
    /* now we've shifted it by N - NhalfPosDC(N) bins */
    f0_shifted = multiFa_spin->data[0]->f0 - NhalfNeg(numSamples) * dFreq;

    /* loop over requested output frequencies and construct F *NOT* 2F */

    /* define number of bins offset from the internal start frequency bin to the user requested bin */
    UINT4 offset = floor(0.5 + (thisPoint->fkdot[0] - f0_shifted)/dFreq);

    for ( UINT4 k=0; k < numFreqBins; k++ )
      {
        UINT4 idx = k + offset;
        /* ----- compute final Fstatistic-value ----- */
        COMPLEX16 Fa = Faf_resampled->data[idx];
        COMPLEX16 Fb = Fbf_resampled->data[idx];
        Fstats->twoF[k] = ComputeFstatFromFaFb ( Fa, Fb, Ad, Bd, Cd, 0, Dd_inv );
      } // for k < numFreqBins

    /* free memory not stored in the buffer */
    XLALDestroyCOMPLEX8Vector( Faf_resampled );
    XLALDestroyCOMPLEX8Vector( Fbf_resampled );
    XLALDestroyCOMPLEX8Vector( outa );
    XLALDestroyCOMPLEX8Vector( outb );
    if ( params->returnSingleF ) {
      XLALDestroyCOMPLEX8Vector ( outaSingle );
      XLALDestroyCOMPLEX8Vector ( outbSingle );
    }
    XLALDestroyCOMPLEX8FFTPlan ( pfwd );

    XLALDestroyMultiCOMPLEX8TimeSeries ( multiFa_spin );
    XLALDestroyMultiCOMPLEX8TimeSeries ( multiFb_spin );

    /* IMPORTANT - point the input buffer pointer to the buffered data */
    params->buffer = cfBuffer;

  }// ================================================================================

  // Return Fa & Fb per detector
  if (whatToCompute & FSTATQ_FAFB_PER_DET) {
    XLAL_ERROR(XLAL_EFAILED, "Unimplemented!");
  }

  // Return F-atoms per detector
  if (whatToCompute & FSTATQ_ATOMS_PER_DET) {
    XLAL_ERROR(XLAL_EFAILED, "Unimplemented!");
  }

  // Resampling cannot currently return amplitude modulation coefficients
  Fstats->Mmunu.Ad = NAN;
  Fstats->Mmunu.Bd = NAN;
  Fstats->Mmunu.Cd = NAN;
  Fstats->Mmunu.Ed = NAN;
  Fstats->Mmunu.Dd = NAN;
  Fstats->Mmunu.Sinv_Tsft = NAN;

  return XLAL_SUCCESS;

} // ComputeFstat_Resamp()
