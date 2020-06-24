/*
 *  Copyright (C) 2012, 2013 John Whelan, Shane Larson and Badri Krishnan
 *  Copyright (C) 2013, 2014 Badri Krishnan, John Whelan, Yuanhao Zhang
 *  Copyright (C) 2016, 2017 Grant David Meadors
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
#ifndef _PULSARCROSSCORRV2_H
#define _PULSARCROSSCORRV2_H

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup PulsarCrossCorr_v2_h Header PulsarCrossCorr_v2.h
 * \ingroup lalpulsar_crosscorr
 * \author John Whelan, Yuanhao Zhang, Shane Larson, Badri Krishnan
 * \date 2012, 2013, 2014
 * \brief Header-file for XLAL routines for v2 CW cross-correlation searches
 *
 */
/** @{ */

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if HAVE_GLOB_H
#include <glob.h>
#endif
#include <time.h>
#include <errno.h>
#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/DetectorSite.h>
#include <lal/LALDatatypes.h>
#include <lal/LALHough.h>
#include <lal/RngMedBias.h>
#include <lal/LALRunningMedian.h>
#include <lal/Velocity.h>
#include <lal/Statistics.h>
#include <lal/ComputeFstat.h>
#include <lal/LALConstants.h>
#include <lal/SFTfileIO.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/LALInitBarycenter.h>
#include <lal/SFTClean.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_trig.h>
#include <lal/FrequencySeries.h>
#include <lal/Sequence.h>
#include <lal/SinCosLUT.h>
#include <lal/LogPrintf.h>
#include <lal/UniversalDopplerMetric.h>
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/FFTWMutex.h>
#include <fftw3.h>
#include "ComputeFstat.h"
#include "GeneratePulsarSignal.h"

/* ******************************************************************
 *  Structure, enum, union, etc., typdefs.
 */

/** Index to refer to an SFT given a set of SFTs from several different detectors */
  typedef struct tagSFTIndex {
    UINT4 detInd; /**< index of detector in list */
    UINT4 sftInd; /**< index of SFT in list for this detector */
  } SFTIndex;

/** List of SFT indices */
  typedef struct tagSFTIndexList {
    UINT4    length; /**< number of SFTs */
    SFTIndex *data; /**< array of SFT indices */
  } SFTIndexList;

/** Index to refer to a pair of SFTs */
  typedef struct tagSFTPairIndex {
#if 0
    SFTIndex sftInd1; /**< index of 1st SFT in pair */
    SFTIndex sftInd2; /**< index of 2nd SFT in pair */
#endif
    UINT4 sftNum[2]; /**< ordinal numbers of first and second SFTs */
  } SFTPairIndex;

/** List of SFT pair indices */
  typedef struct tagSFTPairIndexList {
    UINT4    length; /**< number of SFT Pairs */
    SFTPairIndex *data; /**< array of SFT Pair indices */
  } SFTPairIndexList;


/** Resampling Counter of matching SFTs for a given detector Y_K_X matching SFT K_X */
  typedef struct tagSFTCount {
    UINT4 detInd; /**< original vector index of detector Y */
    UINT4 sftCount; /**< number of matching SFTs */
  } SFTCount;

/* Special multi-types for resampling */
/* These are MULTI-structs */


/** INNER List of SFT indices */
  typedef struct tagResampSFTMultiCountList {
    UINT4       length; /**< number of detectors Y */
    UINT4       sftInd; /**< original vector index of sft K */
    SFTCount  * data;   /**<  arrays of count of SFTs L_Y, at given detector Y, that match each given SFT K_X */
  } ResampSFTMultiCountList;

/** MIDDLE List of SFT indices */
  typedef struct tagResampSFTMultiCountListDet {
    UINT4                      length; /**< number of SFTs K_X */
    UINT4                      detInd; /**< original vector index of detector X */
    ResampSFTMultiCountList  * data;   /**< array of SFT L_Y indices for given SFT K_X at detector X */
  } ResampSFTMultiCountListDet;

/** OUTER List of SFT indices */
  typedef struct tagMultiResampSFTMultiCountList {
    UINT4                         length; /**< number of detectors X */
    ResampSFTMultiCountListDet  * data;   /**< list per-detectors X */
  } MultiResampSFTMultiCountList;


/** Resampling: Index to refer to fields of an SFT given a specific index L_Y_K_X */
  typedef struct tagResampSFTIndex {
    UINT4 detInd;  /**< index of detector in list */
    UINT4 sftInd;  /**< index of SFT in list for this detector L */
    UINT4 flatInd; /**< index in the flat SFT list */
    UINT4 pairInd; /**< index of SFT among overall pairs */
    REAL4 sciFlag; /**< science-or-not value: 1 is science mode, 0 is in a gap */
  } ResampSFTIndex;

/** Resampling: List of SFT indices L for a given detector Y_K_X: indexing method is nominally original vectors but may be affected by gaps */
  typedef struct tagResampSFTIndexList {
    UINT4             length; /**< number of SFTs */
    UINT4             detInd; /**< original vector index of detector Y */
    ResampSFTIndex  * data; /**< array of SFT indices */
  } ResampSFTIndexList;

/** Resampling: Multi List of indices of SFT L_Y, for a given sft K_X  */
  typedef struct tagResampSFTMultiIndexList {
    UINT4                 length;  /**< number of detectors Y */
    UINT4                 sftInd;  /**< original vector index of sft K */
    UINT4                 flatInd; /**< index in the flat SFT list */
    REAL4                 sciFlag; /**< science-or-not value: 1 is science mode, 0 is in a gap */
    ResampSFTIndexList  * data;    /**< array, per detector, of lists of SFT indices */
  } ResampSFTMultiIndexList;

/** Resampling Multi List of SFT pair indices (L_Y_K), for a given detector X */
  typedef struct tagResampSFTPairMultiIndexList {
    UINT4                      length; /**< number of K SFTs at detector X */
    UINT4                      detInd; /**< original vector index of detector X */
    ResampSFTMultiIndexList  * data;  /**< array of SFT L_Y indices for given SFT K_X at detector X */
  } ResampSFTPairMultiIndexList;


/** Resampling Multi List of all paired SFTs L_Y_K_X, top level (multiple pairs, multiple detectors) */
  typedef struct tagMultiResampSFTPairIndexMultiList {
    UINT4                          allPairCount;     /**< count of all pairs */
    UINT4                          oldPairCount;     /**< count of sft pairs, old-style */
    UINT4                          sftTotalCount;    /**< count of all sfts */
    REAL8                          maxLag;           /**< Maximum allowed lag time */
    BOOLEAN                        inclAutoCorr;     /**< Do we include auto-correlation? */
    BOOLEAN                        inclSameDetector; /**< Do we include same detectors? */
    REAL8                          Tsft;             /**< Duration of one SFT */
    REAL8                          Tshort;           /**< Duration of resamp Tshort */
    SFTIndexList                 * indexList;        /**< Make an overall flat index list */
    SFTPairIndexList             * pairIndexList;    /**< Make a classic pair index list */
    UINT4                          length;           /**< number of detectors X */
    ResampSFTPairMultiIndexList  * data;             /**< list per-detector X */
  } MultiResampSFTPairMultiIndexList;

// ----- local types ----------

typedef struct tagCrossCorrTimings_t
{
  REAL8 Total;          //!< total time spent in XLALComputeFstatResamp()
  REAL8 Bary;           //!< time spent in barycentric resampling
  REAL8 Spin;           //!< time spent in spindown+frequency correction
  REAL8 FFT;            //!< time spent in FFT
  REAL8 Norm;           //!< time spent normalizing the final Fa,Fb
  REAL8 Fab2F;          //!< time to compute Fstat from {Fa,Fb}
  REAL8 Mem;            //!< time to realloc and memset-0 arrays
  REAL8 SumFabX;        //!< time to sum_X Fab^X
  REAL8 F1Buf;          //!< Resampling timing 'constant': Fstat time per template per detector for a 'buffered' case (same skypos, same numFreqBins)
  REAL8 F1NoBuf;        //!< Resampling timing 'constant': Fstat time per template per detector for an 'unbuffered' usage (different skypos and numFreqBins)
} CrossCorrTimings_t;

typedef struct tagResampCrossCorrTimingInfo
{ // NOTE: all times refer to a single-detector timing case
  BOOLEAN collectTiming;        //!< turn on/off the collection of F-stat-method-specific timing-data (stored in workspace)

  UINT4 numFreqBins;            //!< number of frequency bins to compute F-stat for
  UINT4 numSFTs;                //!< total number of SFTs used
  UINT4 numDetectors;           //!< number of detectors
  UINT4 numSamplesFFT0;         //!< 'original' length of barycentered timeseries to be FFT'ed
  UINT4 numSamplesFFT;          //!< actual length of FFT, potentially rounded up to power-of-2 for efficiency
  CrossCorrTimings_t tau;
} ResampCrossCorrTimingInfo;

// ----- workspace ----------

typedef struct tagResampCrossCorrWorkspace
{
  // intermediate quantities to interpolate and operate on SRC-frame timeseries
  COMPLEX8Vector *TStmp1_SRC;   //!< can hold a single-detector SRC-frame spindown-corrected timeseries [without zero-padding]
  COMPLEX8Vector *TStmp2_SRC;   //!< can hold a single-detector SRC-frame spindown-corrected timeseries [without zero-padding]
  REAL8Vector *SRCtimes_DET;    //!< holds uniformly-spaced SRC-frame timesteps translated into detector frame [for interpolation]

  // input padded timeseries ts(t) and output Fab(f) of length 'numSamplesFFT' and corresponding fftw plan
  UINT4 numSamplesFFT;          //!< allocated number of zero-padded SRC-frame time samples (related to dFreq)
  UINT4 decimateFFT;            //!< output every n-th frequency bin, with n>1 iff (dFreq > 1/Tspan), and was internally decreased by n
  fftwf_plan fftplan;           //!< buffer FFT plan for given numSamplesOut length
  COMPLEX8 *TS_FFT;             //!< zero-padded, spindown-corr SRC-frame TS
  COMPLEX8 *FabX_Raw;           //!< raw full-band FFT result Fa,Fb

  // arrays of size numFreqBinsOut over frequency bins f_k:
  UINT4 numFreqBinsOut;         //!< number of output frequency bins {f_k}
  COMPLEX8 *FaX_k;              //!< properly normalized F_a^X(f_k) over output bins
  COMPLEX8 *FbX_k;              //!< properly normalized F_b^X(f_k) over output bins
  COMPLEX8 *Fa_k;               //!< properly normalized F_a(f_k) over output bins
  COMPLEX8 *Fb_k;               //!< properly normalized F_b(f_k) over output bins
  UINT4 numFreqBinsAlloc;       //!< internal: keep track of allocated length of frequency-arrays

  ResampCrossCorrTimingInfo *timingInfo; //!< pointer to storage for collecting timing data (which lives in ResampMethodData)
} ResampCrossCorrWorkspace;

/* end Resampling multi-types */

/** A collection of UINT4Vectors -- one for each IFO  */
  /* Probably belongs in SFTUtils.h */
typedef struct tagMultiUINT4Vector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(MultiUINT4Vector, UINT4Vector*, data, UINT4, length));
#endif /* SWIG */
  UINT4        length;  /**< number of ifos */
  UINT4Vector  **data; 	/**< unit4vector for each ifo */
} MultiUINT4Vector;

  /*
 *  Functions Declarations (i.e., prototypes).
 */

int XLALGetDopplerShiftedFrequencyInfo
(
   REAL8Vector            *shiftedFreqs,
   UINT4Vector              *lowestBins,
   COMPLEX8Vector      *expSignalPhases,
   REAL8VectorSequence        *sincList,
   UINT4                        numBins,
   PulsarDopplerParams            *dopp,
   SFTIndexList                   *sfts,
   MultiSFTVector            *inputSFTs,
   MultiSSBtimes            *multiTimes,
   MultiUINT4Vector            *badBins,
   REAL8                           Tsft
   )
  ;

int XLALCreateSFTIndexListFromMultiSFTVect
(
   SFTIndexList        **indexList,
   MultiSFTVector            *sfts
 )
  ;

int XLALCreateSFTPairIndexList
(
   SFTPairIndexList  **pairIndexList,
   SFTIndexList           *indexList,
   MultiSFTVector              *sfts,
   REAL8                      maxLag,
   BOOLEAN              inclAutoCorr
   )
  ;

int XLALCreateSFTPairIndexListResamp
(
   MultiResampSFTPairMultiIndexList  ** resampPairIndexList,
   SFTPairIndexList                  ** pairIndexList,
   SFTIndexList                       * indexList,
   MultiSFTVector                     * sfts,
   REAL8                                maxLag,
   BOOLEAN                              inclAutoCorr,
   BOOLEAN                              inclSameDetector,
   REAL8                                Tsft,
   REAL8                                Tshort
   )
  ;

int XLALCreateSFTPairIndexListShortResamp
(
   MultiResampSFTPairMultiIndexList  **        resampPairIndexList,
   const REAL8                                 maxLag,
   const BOOLEAN                               inclAutoCorr,
   const BOOLEAN                               inclSameDetector,
   const REAL8                                 Tsft,
   const MultiLIGOTimeGPSVector      *_LAL_RESTRICT_ multiTimes
   )
  ;

int XLALTestResampPairIndexList
  (
    MultiResampSFTPairMultiIndexList  * resampMultiPairIndexList
  );

int XLALCalculateCrossCorrGammas
  (
   REAL8Vector      ** Gamma_ave,
   REAL8Vector      ** Gamma_circ,
   SFTPairIndexList  * pairIndexList,
   SFTIndexList      * indexList,
   MultiAMCoeffs     * multiCoeffs
  )
 ;

int XLALCalculateCrossCorrGammasResamp
  (
   REAL8Vector                      ** Gamma_ave,
   REAL8Vector                      ** Gamma_circ,
   MultiResampSFTPairMultiIndexList  * resampMultiPairIndexList,
   MultiAMCoeffs                     * multiCoeffs
  )
 ;

int XLALCalculateCrossCorrGammasResampShort
  (
   REAL8Vector                      ** Gamma_ave,
   REAL8Vector                      ** Gamma_circ,
   MultiResampSFTPairMultiIndexList  * resampMultiPairIndexList,
   MultiAMCoeffs                     * multiCoeffs
  )
 ;

int XLALCalculatePulsarCrossCorrStatistic
  (
   REAL8                         *ccStat,
   REAL8                      *evSquared,
   REAL8Vector                *curlyGAmp,
   COMPLEX8Vector       *expSignalPhases,
   UINT4Vector               *lowestBins,
   REAL8VectorSequence         *sincList,
   SFTPairIndexList            *sftPairs,
   SFTIndexList              *sftIndices,
   MultiSFTVector             *inputSFTs,
   MultiNoiseWeights       *multiWeights,
   UINT4                         numBins
   )
  ;

int XLALCalculatePulsarCrossCorrStatisticResamp
  (
   REAL8Vector                             *_LAL_RESTRICT_ ccStatVector,
   REAL8Vector                             *_LAL_RESTRICT_ evSquaredVector,
   REAL8Vector                             *_LAL_RESTRICT_ numeEquivAve,
   REAL8Vector                             *_LAL_RESTRICT_ numeEquivCirc,
   const REAL8Vector                       *_LAL_RESTRICT_ resampCurlyGAmp,
   const MultiResampSFTPairMultiIndexList  *_LAL_RESTRICT_ resampPairs,
   const MultiNoiseWeights                 *_LAL_RESTRICT_ multiWeights,
   const PulsarDopplerParams               *_LAL_RESTRICT_ binaryTemplateSpacings,
   const PulsarDopplerParams               *_LAL_RESTRICT_ dopplerpos,
   const MultiCOMPLEX8TimeSeries           *_LAL_RESTRICT_ multiTimeSeries_SRC_a,
   const MultiCOMPLEX8TimeSeries           *_LAL_RESTRICT_ multiTimeSeries_SRC_b,
   ResampCrossCorrWorkspace                *_LAL_RESTRICT_ ws,
   COMPLEX8                                *_LAL_RESTRICT_ ws1KFaX_k,
   COMPLEX8                                *_LAL_RESTRICT_ ws1KFbX_k,
   COMPLEX8                                *_LAL_RESTRICT_ ws2LFaX_k,
   COMPLEX8                                *_LAL_RESTRICT_ ws2LFbX_k
   )
  ;

int XLALCalculateCrossCorrPhaseDerivatives
  (
   REAL8VectorSequence           ** phaseDerivs,
   const PulsarDopplerParams      * dopplerPoint,
   const EphemerisData            * edat,
   SFTIndexList                   * indexList,
   MultiSSBtimes                  * multiTimes,
   const DopplerCoordinateSystem  * coordSys
   )
  ;

int XLALCalculateCrossCorrPhaseDerivativesShort
  (
   REAL8VectorSequence              ** resampPhaseDerivs,
   const PulsarDopplerParams         * dopplerPoint,
   const EphemerisData               * edat,
   SFTIndexList                      * indexList,
   MultiResampSFTPairMultiIndexList  * resampMultiPairIndexList,
   MultiSSBtimes                     * multiTimes,
   const DopplerCoordinateSystem     * coordSys
   )
  ;

int XLALCalculateCrossCorrPhaseMetric
  (
   gsl_matrix                        **g_ij,
   gsl_vector                       **eps_i,
   REAL8                        *sumGammaSq,
   const REAL8VectorSequence   *phaseDerivs,
   const SFTPairIndexList    *pairIndexList,
   const REAL8Vector             *Gamma_ave,
   const REAL8Vector            *Gamma_circ,
   const DopplerCoordinateSystem  *coordSys
   );

int XLALCalculateCrossCorrPhaseMetricShort
  (
   gsl_matrix                            ** g_ij,
   gsl_vector                            ** eps_i,
   REAL8                                  * sumGammaSq,
   const REAL8VectorSequence              * phaseDerivs,
   const MultiResampSFTPairMultiIndexList * resampMultiPairs,
   const REAL8Vector                      * Gamma_ave,
   const REAL8Vector                      * Gamma_circ,
   const DopplerCoordinateSystem          * coordSys
   );

int XLALCalculateLMXBCrossCorrDiagMetric
  (
   REAL8                      *hSens,
   REAL8                       *g_ff,
   REAL8                       *g_aa,
   REAL8                       *g_TT,
   REAL8                       *g_pp,
   REAL8             *weightedMuTAve,
   PulsarDopplerParams DopplerParams,
   REAL8Vector              *G_alpha,
   SFTPairIndexList   *pairIndexList,
   SFTIndexList           *indexList,
   MultiSFTVector              *sfts,
   MultiNoiseWeights   *multiWeights
   )
  ;

int XLALCalculateLMXBCrossCorrDiagMetricShort
  (
   REAL8                                   *         hSens,
   REAL8                                   *         g_ff,
   REAL8                                   *         g_aa,
   REAL8                                   *         g_TT,
   REAL8                                   *         g_pp,
   const PulsarDopplerParams                         DopplerParams,
   const REAL8Vector                       *_LAL_RESTRICT_ G_alpha,
   const MultiResampSFTPairMultiIndexList  *_LAL_RESTRICT_ resampMultiPairIndexList,
   const MultiLIGOTimeGPSVector            *_LAL_RESTRICT_ timestamps,
   const MultiNoiseWeights                 *_LAL_RESTRICT_ multiWeights
   )
  ;

/* (possible future function) */
//int
//XLALBesselCrossCorrOrbitalSpaceStep( 
//     COMPLEX8             * xTildeOut,
//     COMPLEX8             * xTildeIn,
//     PulsarDopplerParams  * dopplerpos,
//     PulsarDopplerParams  * binaryTemplateSpacings,
//     INT4                 aStep,
//     INT4                 pStep,
//     INT4                 tStep,
//     INT4                 nFFT,
//     UINT4                nTermsMax
//);

LIGOTimeGPSVector 
*XLALModifyTimestampsFromSFTsShort ( 
    REAL8TimeSeries          **        sciFlag,
    const LIGOTimeGPSVector  *_LAL_RESTRICT_ Times,
    const REAL8                        tShort,
    const UINT4                        numShortPerDet
);

LIGOTimeGPSVector
*XLALExtractTimestampsFromSFTsShort (
    REAL8TimeSeries ** sciFlag,
    const SFTVector  * sfts,
    REAL8              tShort,
    UINT4              numShortPerDet
);

MultiLIGOTimeGPSVector 
*XLALModifyMultiTimestampsFromSFTs ( 
    MultiREAL8TimeSeries          **        scienceFlagVect,
    const MultiLIGOTimeGPSVector  *_LAL_RESTRICT_ multiTimes,
    const REAL8                             tShort,
    const UINT4                             numShortPerDet
);

MultiLIGOTimeGPSVector 
*XLALExtractMultiTimestampsFromSFTsShort ( 
    MultiREAL8TimeSeries  ** scienceFlagVect,
    const MultiSFTVector  * multiSFTs,
    REAL8                   tShort,
    UINT4                   numShortPerDet
);

int 
XLALFillDetectorTensorShort (
    DetectorState      * detState, 
    const LALDetector  * detector 
);

DetectorStateSeries* 
XLALGetDetectorStatesShort ( 
    const LIGOTimeGPSVector  * timestamps, 
    const LALDetector        * detector, 
    const EphemerisData      * edat, 
    REAL8                      tOffset,
    REAL8                      tShort,
    UINT4                      numShortPerDet
);

MultiDetectorStateSeries* 
XLALGetMultiDetectorStatesShort ( 
    const MultiLIGOTimeGPSVector  * multiTS, 
    const MultiLALDetector        * multiIFO, 
    const EphemerisData           * edat, 
    REAL8                           tOffset,
    REAL8                           tShort,
    UINT4                           numShortPerDet
);


int
XLALModifyAMCoeffsWeights (
    REAL8Vector                   **         resampMultiWeightsX,
    const MultiNoiseWeights       *_LAL_RESTRICT_  multiWeights,
    const REAL8                              tShort,
    const REAL8                              tSFTOld,
    const UINT4                              numShortPerDet,
    const MultiLIGOTimeGPSVector  *_LAL_RESTRICT_  multiTimes,
    const UINT4                              maxNumStepsOldIfGapless,
    const UINT4                              X
);

int
XLALModifyMultiAMCoeffsWeights(
     MultiNoiseWeights             **        multiWeights,
     const REAL8                             tShort,
     const REAL8                             tSFTOld,
     const UINT4                             numShortPerDet,
     const MultiLIGOTimeGPSVector  *_LAL_RESTRICT_ multiTimes
);


int 
XLALWeightMultiAMCoeffsShort (  
     MultiAMCoeffs            * multiAMcoef, 
     const MultiNoiseWeights  * multiWeights, 
     REAL8                      tShort,
     REAL8                      tSFTOld,
     UINT4                      numShortPerDet,
     MultiLIGOTimeGPSVector   * multiTimes
);

MultiAMCoeffs 
*XLALComputeMultiAMCoeffsShort ( 
     const MultiDetectorStateSeries  * multiDetStates, 
     const MultiNoiseWeights         * multiWeights, 
     SkyPosition                       skypos,
     REAL8                             tShort,
     REAL8                             tSFTOld,
     UINT4                             numShortPerDet,
     MultiLIGOTimeGPSVector          * multiTimes
);

AMCoeffs 
*XLALComputeAMCoeffsShort ( 
     const DetectorStateSeries  * DetectorStates, 
     SkyPosition                  skypos 
);


UINT4
XLALCrossCorrNumShortPerDetector( 
    const REAL8 resampTshort,
    const INT4  startTime,
    const INT4  endTime
);

REAL8TimeSeries
*XLALCrossCorrGapFinderResamp(
    LIGOTimeGPSVector        *_LAL_RESTRICT_ timestamps,
    const LIGOTimeGPSVector  *_LAL_RESTRICT_ Times
);

REAL8TimeSeries
*XLALCrossCorrGapFinderResampAlt(
    LIGOTimeGPSVector *_LAL_RESTRICT_ timestamps,
    const SFTVector   *_LAL_RESTRICT_ sfts
);

int
XLALEquipCrossCorrPairsWithScienceFlags(
    MultiResampSFTPairMultiIndexList  * resampMultiPairs,
     MultiREAL8TimeSeries             * scienceFlagVect
);

///* (possible future function) */
//MultiSFTVector
//*XLALModifyCrossCorrTimestampsIntoSFTVector(
//    const MultiLIGOTimeGPSVector *multiTimes
//);

int
XLALCreateCrossCorrWorkspace( 
    ResampCrossCorrWorkspace   **        ws,
    COMPLEX8                   **        ws1KFaX_kOut,
    COMPLEX8                   **        ws1KFbX_kOut,
    COMPLEX8                   **        ws2LFaX_kOut,
    COMPLEX8                   **        ws2LFbX_kOut,
    MultiCOMPLEX8TimeSeries    **        multiTimeSeries_SRC_a,
    MultiCOMPLEX8TimeSeries    **        multiTimeSeries_SRC_b,
    const PulsarDopplerParams            binaryTemplateSpacings,
    const FstatInput           *_LAL_RESTRICT_ resampFstatInput, 
    const UINT4                          numFreqBins,
    const REAL8                          tCoh,
    const BOOLEAN                        treatWarningsAsErrors
);

/** @} */

void XLALDestroySFTIndexList ( SFTIndexList *sftIndices );

void XLALDestroySFTPairIndexList ( SFTPairIndexList *sftPairs );

void XLALDestroyResampSFTIndexList( ResampSFTIndexList *sftResampList );

void XLALDestroyResampSFTMultiIndexList( ResampSFTMultiIndexList *sftResampMultiList );

void XLALDestroyResampSFTPairMultiIndexList( ResampSFTPairMultiIndexList *sftResampPairMultiList  );

void XLALDestroyMultiResampSFTPairMultiIndexList ( MultiResampSFTPairMultiIndexList *sftMultiPairsResamp );

void XLALDestroyMultiMatchList( MultiResampSFTMultiCountList * localMultiListOfLmatchingGivenMultiK );

void XLALDestroyResampCrossCorrWorkspace ( void *workspace );

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif


#endif     /* Close double-include protection _PULSARCROSSCORR_H */
