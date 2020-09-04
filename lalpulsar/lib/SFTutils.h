/*
 *  Copyright (C) 2014 Karl Wette
 *  Copyright (C) 2009 Chris Messenger
 *  Copyright (C) 2005 Reinhard Prix
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
#ifndef _SFTUTILS_H  /* Double-include protection. */
#define _SFTUTILS_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup SFTutils_h Header SFTutils.h
 * \ingroup lalpulsar_sft
 * \author Reinhard Prix, Badri Krishnan, Badri Krishnan, Iraj Gholami,
 * Reinhard Prix, Alicia Sintes, Karl Wette, David Keitel
 * \date 2005-2020
 * \brief Utility functions for handling of SFTs and associated structures
 *
 * This module contains various helper functions to create, handle, combine,
 * and destroy SFTs (Short Fourier Transforms) and related data structures,
 * including SFTtype, SFTVector, SFTCatalog
 * (and their multi-detector generalizations)
 * as well as tools for dealing with timestamps, segments
 * and ASD/PSD (Amplitude/Power Spectral Density) estimates.
 *
 */
/** @{ */

/*---------- INCLUDES ----------*/
#include <stdarg.h>

#include <lal/LALDatatypes.h>
#include <lal/DetectorSite.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/RngMedBias.h>
#include <lal/LALRunningMedian.h>
#include <lal/Segments.h>
#include <lal/SFTfileIO.h>
#include <lal/StringVector.h>
#include <lal/UserInputParse.h>

/*---------- DEFINES ----------*/

/*---------- exported types ----------*/

/** A vector of REAL8FrequencySeries */
typedef struct tagREAL8FrequencySeriesVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(REAL8FrequencySeriesVector, REAL8FrequencySeries, data, UINT4, length));
#endif /* SWIG */
  UINT4                  length;
  REAL8FrequencySeries   *data;
} REAL8FrequencySeriesVector;


/** Special type for holding a PSD vector (over several SFTs) */
typedef REAL8FrequencySeriesVector PSDVector;

/** A collection of PSD vectors -- one for each IFO in a multi-IFO search */
typedef struct tagMultiPSDVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(MultiPSDVector, PSDVector*, data, UINT4, length));
#endif /* SWIG */
  UINT4      length;  	/**< number of ifos */
  PSDVector  **data; 	/**< sftvector for each ifo */
} MultiPSDVector;

/** One noise-weight (number) per SFT (therefore indexed over IFOs and SFTs */
typedef struct tagMultiNoiseWeights {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(MultiNoiseWeights, REAL8Vector*, data, UINT4, length));
#endif /* SWIG */
  UINT4 length;		/**< number of detectors */
  REAL8Vector **data;	/**< weights-vector for each detector */
  REAL8 Sinv_Tsft;	/**< normalization factor used: \f$\mathcal{S}^{-1}\,T_\mathrm{SFT}\f$ (using single-sided PSD!) */
  BOOLEAN isNotNormalized;  /**< if true: weights are saved unnormalized (divide by Sinv_Tsft to get normalized version). */
} MultiNoiseWeights;

/** common types of mathematical operations over an array */
typedef enum tagMathOpType {
  MATH_OP_ARITHMETIC_SUM = 0,   /**< \f$\sum_k x_k\f$      */
  MATH_OP_ARITHMETIC_MEAN,      /**< \f$\sum_k x_k / N\f$ */
  MATH_OP_ARITHMETIC_MEDIAN,    /**< \f$x_1 \leq \dots \leq  x_{N/2} \leq \dots \leq x_n\f$ */
  MATH_OP_HARMONIC_SUM,         /**< \f$1 / \sum_k (1/x_k)\f$ */
  MATH_OP_HARMONIC_MEAN,        /**< \f$N / \sum_k (1/x_k)\f$ */
  MATH_OP_POWERMINUS2_SUM,      /**< \f$1 / \sqrt{ \sum_k (1/x_k^2) }\f$ */
  MATH_OP_POWERMINUS2_MEAN,     /**< \f$1 / \sqrt{ \sum_k (1/x_k^2) / N }\f$ */
  MATH_OP_MINIMUM,              /**< \f$\min_k(x_k)\f$ */
  MATH_OP_MAXIMUM,              /**< \f$\max_k(x_k)\f$ */
  MATH_OP_LAST
} MathOpType;

extern const UserChoices MathOpTypeChoices;

/*---------- Global variables ----------*/

/*---------- exported prototypes [API] ----------*/
/* ----------------------------------------------------------------------
 *  some prototypes for general functions handling these data-types
 *----------------------------------------------------------------------*/
SFTtype* XLALCreateSFT ( UINT4 numBins );
SFTVector* XLALCreateSFTVector (UINT4 numSFTs, UINT4 numBins );
MultiSFTVector *XLALCreateMultiSFTVector ( UINT4 length, UINT4Vector *numsft );

int XLALAppendSFT2Vector (SFTVector *vect, const SFTtype *sft );

void XLALDestroySFT (SFTtype *sft);
void XLALDestroySFTVector (SFTVector *vect);
void XLALDestroyMultiSFTVector ( MultiSFTVector *multvect );

SFTVector *XLALDuplicateSFTVector ( const SFTVector *sftsIn );

int XLALReorderMultiSFTVector( MultiSFTVector *multiSFTs, const LALStringVector *IFOs);

COMPLEX8Vector *XLALrefineCOMPLEX8Vector (const COMPLEX8Vector *in, UINT4 refineby, UINT4 Dterms);

int XLALExtractBandFromSFT ( SFTtype **outSFT, const SFTtype *inSFT, REAL8 fMin, REAL8 Band );
SFTVector *XLALExtractBandFromSFTVector ( const SFTVector *inSFTs, REAL8 fMin, REAL8 Band );
MultiSFTVector *XLALExtractBandFromMultiSFTVector ( const MultiSFTVector *inSFTs, REAL8 fMin, REAL8 Band );
int XLALFindCoveringSFTBins ( UINT4 *firstBin, UINT4 *numBins, REAL8 fMinIn, REAL8 BandIn, REAL8 Tsft );

LIGOTimeGPSVector *XLALCreateTimestampVector (UINT4 len);
LIGOTimeGPSVector *XLALResizeTimestampVector ( LIGOTimeGPSVector *vector, UINT4 length );
LIGOTimeGPSVector *XLALMakeTimestamps ( LIGOTimeGPS tStart, REAL8 Tspan, REAL8 Tsft, REAL8 Toverlap );
MultiLIGOTimeGPSVector *XLALMakeMultiTimestamps ( LIGOTimeGPS tStart, REAL8 Tspan, REAL8 Tsft, REAL8 Toverlap, UINT4 numDet );

LIGOTimeGPSVector *XLALExtractTimestampsFromSFTs ( const SFTVector *sfts );
MultiLIGOTimeGPSVector *XLALExtractMultiTimestampsFromSFTs ( const MultiSFTVector *multiSFTs );

LIGOTimeGPSVector *XLALTimestampsFromSFTCatalog ( const SFTCatalog *catalog );
MultiLIGOTimeGPSVector *XLALTimestampsFromMultiSFTCatalogView ( const MultiSFTCatalogView *multiView );

LIGOTimeGPSVector *XLALTimestampsFromSegmentFile( const char *filename, REAL8 Tsft, REAL8 Toverlap, BOOLEAN adjustSegExtraTime, BOOLEAN synchronize );

void XLALDestroyTimestampVector (LIGOTimeGPSVector *vect);
void XLALDestroyMultiTimestamps ( MultiLIGOTimeGPSVector *multiTS );

int XLALFindCWDetector ( CHAR** prefix, INT4 *lalCachedIndex, const CHAR *name, const BOOLEAN exactMatch );
BOOLEAN XLALIsValidCWDetector ( const CHAR *name );
CHAR *XLALGetChannelPrefix ( const CHAR *name );
LALDetector *XLALGetSiteInfo ( const CHAR *name );

LALSegList *XLALReadSegmentsFromFile ( const char *fname );

// adding SFTs
int XLALMultiSFTVectorAdd ( MultiSFTVector *a, const MultiSFTVector *b );
int XLALSFTVectorAdd ( SFTVector *a, const SFTVector *b );
int XLALSFTAdd ( SFTtype *a, const SFTtype *b );

int XLALEarliestMultiSFTsample ( LIGOTimeGPS *out, const MultiSFTVector *multisfts );
int XLALLatestMultiSFTsample ( LIGOTimeGPS *out, const MultiSFTVector *multisfts );

// resizing SFTs
int XLALMultiSFTVectorResizeBand ( MultiSFTVector *multiSFTs, REAL8 f0, REAL8 Band );
int XLALSFTVectorResizeBand ( SFTVector *SFTs, REAL8 f0, REAL8 Band );
int XLALSFTResizeBand ( SFTtype *SFT, REAL8 f0, REAL8 Band );

// destructors
void XLALDestroyPSDVector ( PSDVector *vect );
void XLALDestroyMultiPSDVector ( MultiPSDVector *multvect );

MultiNoiseWeights *XLALComputeMultiNoiseWeights ( const MultiPSDVector *rngmed, UINT4 blocksRngMed, UINT4 excludePercentile);

void XLALDestroyMultiNoiseWeights ( MultiNoiseWeights *weights );

SFTCatalog *XLALAddToFakeSFTCatalog( SFTCatalog *catalog, const CHAR *detector, const LIGOTimeGPSVector *timestamps );
SFTCatalog *XLALMultiAddToFakeSFTCatalog( SFTCatalog *catalog, const LALStringVector *detectors, const MultiLIGOTimeGPSVector *timestamps );
int XLALCopySFT ( SFTtype *dest, const SFTtype *src );

int XLALSFTCatalogTimeslice( SFTCatalog *slice, const SFTCatalog *catalog, const LIGOTimeGPS *minStartGPS, const LIGOTimeGPS *maxStartGPS );
int XLALFindTimesliceBounds ( UINT4 *iStart, UINT4 *iEnd, const LIGOTimeGPSVector *timestamps, const LIGOTimeGPS *minStartGPS, const LIGOTimeGPS *maxStartGPS );

SFTVector *XLALExtractSFTVectorWithTimestamps ( const SFTVector *sfts, const LIGOTimeGPSVector *timestamps );
MultiSFTVector *XLALExtractMultiSFTVectorWithMultiTimestamps ( const MultiSFTVector *multiSFTs, const MultiLIGOTimeGPSVector *multiTimestamps );

int XLALValidateSFTFile ( const char *fname );

// compute and work with PSDs
int XLALComputePSDandNormSFTPower ( REAL8Vector **finalPSD,
                                    MultiPSDVector **multiPSDVector,
                                    REAL8Vector **normSFT,
                                    MultiSFTVector *inputSFTs,
                                    const BOOLEAN returnMultiPSDVector,
                                    const BOOLEAN returnNormSFT,
                                    const UINT4 blocksRngMed,
                                    const MathOpType PSDmthopSFTs,
                                    const MathOpType PSDmthopIFOs,
                                    const MathOpType nSFTmthopSFTs,
                                    const MathOpType nSFTmthopIFOs,
                                    const BOOLEAN normalizeByTotalNumSFTs,
                                    const REAL8 FreqMin,
                                    const REAL8 FreqBand
                              );
int XLALDumpMultiPSDVector ( const CHAR *outbname, const MultiPSDVector *multiPSDVect );
int XLALCropMultiPSDandSFTVectors ( MultiPSDVector *multiPSDVect, MultiSFTVector *multiSFTVect, UINT4 firstBin, UINT4 lastBin );
REAL8FrequencySeries *XLALComputeSegmentDataQ ( const MultiPSDVector *multiPSDVect, LALSeg segment );
REAL8 XLALMathOpOverArray(const REAL8* data, const size_t length, const MathOpType optype);
REAL8 XLALGetMathOpNormalizationFactorFromTotalNumberOfSFTs (  const UINT4 totalNumSFTs, const MathOpType optypeSFTs );
int XLALWritePSDtoFilePointer ( FILE *fpOut,
                                REAL8Vector *PSDVect,
                                REAL8Vector *normSFTVect,
                                BOOLEAN outputNormSFT,
                                BOOLEAN outFreqBinEnd,
                                INT4 PSDmthopBins,
                                INT4 nSFTmthopBins,
                                INT4 binSize,
                                INT4 binStep,
                                REAL8 Freq0,
                                REAL8 dFreq
                              );

/** @} */

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
