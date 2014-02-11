/*
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
 * \ingroup pkg_SFTIO
 * \author Reinhard Prix, Badri Krishnan
 * \date 2005
 * \brief Utility functions for handling of SFTtype and SFTVectors
 *
 * The helper functions XLALCreateSFT(), XLALDestroySFT(), XLALCreateSFTVector()
 * and XLALDestroySFTVector() respectively allocate and free SFT-structs and SFT-vectors.
 * Similarly, XLALCreateTimestampVector() and XLALDestroyTimestampVector() allocate and free
 * a bunch of GPS-timestamps.
 *
 */
/*@{*/

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
  UINT4 length;		/**< number of ifos */
  REAL8Vector **data;	/**< weights-vector for each SFTs */
  REAL8 Sinv_Tsft;	/**< normalization factor used: \f$\mathcal{S}^{-1}\,T_\mathrm{SFT}\f$ (using single-sided PSD!) */
} MultiNoiseWeights;

/*---------- Global variables ----------*/
/* empty init-structs for the types defined in here */
extern const PSDVector empty_PSDVector;
extern const MultiPSDVector empty_MultiPSDVector;
extern const MultiNoiseWeights empty_MultiNoiseWeights;


// ---------- obsolete LAL-API was moved into external file
#ifndef SWIG /* exclude from SWIG interface */
#include "SFTutils-LAL.h"
#endif
// ------------------------------

/*---------- exported prototypes [API] ----------*/
/* ----------------------------------------------------------------------
 *  some prototypes for general functions handling these data-types
 *----------------------------------------------------------------------*/
SFTtype* XLALCreateSFT ( UINT4 numBins );
SFTVector* XLALCreateSFTVector (UINT4 numSFTs, UINT4 numBins );

void XLALDestroySFT (SFTtype *sft);
void XLALDestroySFTVector (SFTVector *vect);

COMPLEX8Vector *XLALrefineCOMPLEX8Vector (const COMPLEX8Vector *in, UINT4 refineby, UINT4 Dterms);
SFTVector *XLALExtractBandFromSFTVector ( const SFTVector *inSFTs, REAL8 fMin, REAL8 Band );
int XLALFindCoveringSFTBins ( UINT4 *firstBin, UINT4 *numBins, REAL8 fMinIn, REAL8 BandIn, REAL8 Tsft );

LIGOTimeGPSVector *XLALCreateTimestampVector (UINT4 len);
LIGOTimeGPSVector *XLALMakeTimestamps ( LIGOTimeGPS tStart, REAL8 Tspan, REAL8 Tsft, REAL8 Toverlap );
MultiLIGOTimeGPSVector *XLALMakeMultiTimestamps ( LIGOTimeGPS tStart, REAL8 Tspan, REAL8 Tsft, REAL8 Toverlap, UINT4 numDet );

LIGOTimeGPSVector *XLALExtractTimestampsFromSFTs ( const SFTVector *sfts );
MultiLIGOTimeGPSVector *XLALExtractMultiTimestampsFromSFTs ( const MultiSFTVector *multiSFTs );

LIGOTimeGPSVector *XLALTimestampsFromSFTCatalog ( const SFTCatalog *catalog );
MultiLIGOTimeGPSVector *XLALTimestampsFromMultiSFTCatalogView ( const MultiSFTCatalogView *multiView );

LIGOTimeGPSVector *XLALTimestampsFromSegmentFile( const char *filename, REAL8 Tsft, REAL8 Toverlap, INT4 adjustSegExtraTime, INT4 synchronize );

void XLALDestroyTimestampVector (LIGOTimeGPSVector *vect);
void XLALDestroyMultiTimestamps ( MultiLIGOTimeGPSVector *multiTS );

CHAR *XLALGetChannelPrefix ( const CHAR *name );
LALDetector *XLALGetSiteInfo ( const CHAR *name );

LALSegList *XLALReadSegmentsFromFile ( const char *fname );

// adding SFTs
int XLALMultiSFTVectorAdd ( MultiSFTVector *a, const MultiSFTVector *b );
int XLALSFTVectorAdd ( SFTVector *a, const SFTVector *b );
int XLALSFTAdd ( SFTtype *a, const SFTtype *b );

int XLALEarliestMultiSFTsample ( LIGOTimeGPS *out, MultiSFTVector *multisfts );
int XLALLatestMultiSFTsample ( LIGOTimeGPS *out, MultiSFTVector *multisfts );

// destructors
void XLALDestroyPSDVector ( PSDVector *vect );
void XLALDestroyMultiSFTVector ( MultiSFTVector *multvect );
void XLALDestroyMultiPSDVector ( MultiPSDVector *multvect );

MultiNoiseWeights *XLALComputeMultiNoiseWeights ( const MultiPSDVector *rngmed, UINT4 blocksRngMed, UINT4 excludePercentile);

void XLALDestroyMultiNoiseWeights ( MultiNoiseWeights *weights );

/*@}*/

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
