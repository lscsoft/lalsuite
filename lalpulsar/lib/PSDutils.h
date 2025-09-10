/*
 * Copyright (C) 2013, 2020 David Keitel
 * Copyright (C) 2010, 2013, 2014, 2016, 2017, 2022 Karl Wette
 * Copyright (C) 2010 Chris Messenger
 * Copyright (C) 2007, 2013 John T. Whelan
 * Copyright (C) 2006 Badri Krishnan
 * Copyright (C) 2004--2006, 2008--2013, 2016 Reinhard Prix
 * Copyright (C) 2004, 2005 Bernd Machenschalk
 * Copyright (C) 2004, 2005 Alicia Sintes
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA
 */

#ifndef _PSDUTILS_H     /* Double-include protection. */
#define _PSDUTILS_H

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

/* includes */
#include <stdlib.h>

#include <lal/LALStdlib.h>
#include <lal/SFTfileIO.h>

/**
 * \defgroup PSDutils_h Header PSDutils.h
 * \ingroup lalpulsar_sft
 *
 * \brief Module for computing PSD (Power Spectral Density) estimates and handling related structures
 *
 * Overview:
 * \ref PSD-type-cdtor-func "create/destroy functions",
 * \ref PSD-type-mod-func "modify functions",
 * \ref PSD-type-prop-func "property functions",
 * \ref PSD-type-gen-func "generation functions",
 * \ref PSD-file-write-func "file writing functions".
 */

/** @{ */

/*---------- exported types ----------*/

/** Special type for holding a PSD vector (over several SFTs) */
typedef REAL8FrequencySeriesVector PSDVector;

/** A collection of PSD vectors -- one for each IFO in a multi-IFO search */
typedef struct tagMultiPSDVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL( ARRAY_1D( MultiPSDVector, PSDVector *, data, UINT4, length ) );
#endif /* SWIG */
  UINT4      length;    /**< number of ifos */
  PSDVector  **data;    /**< sftvector for each ifo */
} MultiPSDVector;

/** One noise-weight (number) per SFT (therefore indexed over IFOs and SFTs */
typedef struct tagMultiNoiseWeights {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL( ARRAY_1D( MultiNoiseWeights, REAL8Vector *, data, UINT4, length ) );
#endif /* SWIG */
  UINT4 length;         /**< number of detectors */
  REAL8Vector **data;   /**< weights-vector for each detector */
  REAL8 Sinv_Tsft;      /**< normalization factor used: \f$ \mathcal{S}^{-1}\,T_\mathrm{SFT} \f$ (using single-sided PSD!) */
  BOOLEAN isNotNormalized;  /**< if true: weights are saved unnormalized (divide by Sinv_Tsft to get normalized version). */
} MultiNoiseWeights;

/** common types of mathematical operations over an array */
typedef enum tagMathOpType {
  MATH_OP_ARITHMETIC_SUM = 0,   /**< \f$ \sum_k x_k \f$ */
  MATH_OP_ARITHMETIC_MEAN,      /**< \f$ \sum_k x_k / N \f$ */
  MATH_OP_ARITHMETIC_MEDIAN,    /**< \f$ x_1 \leq \dots \leq  x_{N/2} \leq \dots \leq x_n \f$ */
  MATH_OP_HARMONIC_SUM,         /**< \f$ 1 / \sum_k (1/x_k) \f$ */
  MATH_OP_HARMONIC_MEAN,        /**< \f$ N / \sum_k (1/x_k) \f$ */
  MATH_OP_POWERMINUS2_SUM,      /**< \f$ 1 / \sqrt{ \sum_k (1/x_k^2) } \f$ */
  MATH_OP_POWERMINUS2_MEAN,     /**< \f$ 1 / \sqrt{ \sum_k (1/x_k^2) / N } \f$ */
  MATH_OP_MINIMUM,              /**< \f$ \min_k(x_k) \f$ */
  MATH_OP_MAXIMUM,              /**< \f$ \max_k(x_k) \f$ */
  MATH_OP_LAST
} MathOpType;

/*---------- global variables ----------*/

extern const UserChoices MathOpTypeChoices;

/*---------- exported prototypes [API] ----------*/

/**
 * \name PSD type create/destroy functions
 * \anchor PSD-type-cdtor-func
 */
/** @{ */

void XLALDestroyPSDVector( PSDVector *vect );
void XLALDestroyMultiPSDVector( MultiPSDVector *multvect );

MultiNoiseWeights *XLALCreateMultiNoiseWeights( const UINT4 length );
MultiNoiseWeights *XLALCopyMultiNoiseWeights( const MultiNoiseWeights *multiWeights );
void XLALDestroyMultiNoiseWeights( MultiNoiseWeights *weights );

/** @} */

/**
 * \name PSD type modify functions
 * \anchor PSD-type-mod-func
 */
/** @{ */

int XLALCropMultiPSDandSFTVectors( MultiPSDVector *multiPSDVect, MultiSFTVector *multiSFTVect, UINT4 firstBin, UINT4 lastBin );

/** @} */

/**
 * \name PSD type property functions
 * \anchor PSD-type-prop-func
 */
/** @{ */

MultiNoiseWeights *XLALComputeMultiNoiseWeights( const MultiPSDVector *rngmed, UINT4 blocksRngMed, UINT4 excludePercentile );

/** @} */

/**
 * \name PSD generation functions
 * \anchor PSD-type-gen-func
 */
/** @{ */

REAL8FrequencySeries *XLALComputeSegmentDataQ( const MultiPSDVector *multiPSDVect, LALSeg segment );
REAL8 XLALMathOpOverArray( const REAL8 *data, const size_t length, const MathOpType optype );
REAL8 XLALMathOpOverREAL8Vector( const REAL8Vector *data, const MathOpType optype );
REAL8 XLALGetMathOpNormalizationFactorFromTotalNumberOfSFTs( const UINT4 totalNumSFTs, const MathOpType optypeSFTs );

int XLALComputePSDandNormSFTPower( REAL8Vector **finalPSD, MultiPSDVector **multiPSDVector, REAL8Vector **normSFT, MultiSFTVector *inputSFTs, const BOOLEAN returnMultiPSDVector, const BOOLEAN returnNormSFT,
                                   const UINT4 blocksRngMed, const MathOpType PSDmthopSFTs, const MathOpType PSDmthopIFOs, const MathOpType nSFTmthopSFTs, const MathOpType nSFTmthopIFOs,
                                   const BOOLEAN normalizeByTotalNumSFTs, const REAL8 FreqMin, const REAL8 FreqBand, const BOOLEAN normalizeSFTsInPlace );
int XLALComputePSDfromSFTs( REAL8Vector **finalPSD, MultiSFTVector *inputSFTs, const UINT4 blocksRngMed, const MathOpType PSDmthopSFTs, const MathOpType PSDmthopIFOs,
                            const BOOLEAN normalizeByTotalNumSFTs, const REAL8 FreqMin, const REAL8 FreqBand );

/** @} */

/**
 * \name PSD file writing functions
 * \anchor PSD-file-write-func
 */
/** @{ */

int XLALDumpMultiPSDVector( const CHAR *outbname, const MultiPSDVector *multiPSDVect );

int XLALWritePSDtoFilePointer( FILE *fpOut, REAL8Vector *PSDVect, REAL8Vector *normSFTVect, BOOLEAN outputNormSFT, BOOLEAN outFreqBinEnd, INT4 PSDmthopBins, INT4 nSFTmthopBins,
                               INT4 binSize, INT4 binStep, REAL8 Freq0, REAL8 dFreq );

/** @} */

/** @} */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection */
