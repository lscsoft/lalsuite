/*
 * Copyright (C) 2005 Reinhard Prix
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

/*---------- DEFINES ----------*/

/*---------- exported types ----------*/

/** A vector of COMPLEX8FrequencySeries */
typedef struct tagCOMPLEX8FrequencySeriesVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(COMPLEX8FrequencySeries, data, UINT4, length));
#endif /* SWIG */
  UINT4 			length;		/**< number of SFTs */
  COMPLEX8FrequencySeries 	*data;		/**< array of SFTs */
} COMPLEX8FrequencySeriesVector;

/** A vector of REAL8FrequencySeries */
typedef struct tagREAL8FrequencySeriesVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(REAL8FrequencySeries, data, UINT4, length));
#endif /* SWIG */
  UINT4                  length;
  REAL8FrequencySeries   *data;
} REAL8FrequencySeriesVector;

/** A vector of REAL4FrequencySeries */
typedef struct tagREAL4FrequencySeriesVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(REAL4FrequencySeries, data, UINT4, length));
#endif /* SWIG */
  UINT4                  length;
  REAL4FrequencySeries   *data;
} REAL4FrequencySeriesVector;


/** A so-called 'SFT' (short-Fourier-transform) will be stored in a COMPLEX8FrequencySeries */
typedef COMPLEX8FrequencySeries 	SFTtype;


/** The corresponding vector-type to hold a vector of 'SFTs' */
typedef COMPLEX8FrequencySeriesVector 	SFTVector;

/** Special type for holding a PSD vector (over several SFTs) */
typedef REAL8FrequencySeriesVector PSDVector;

/** A collection of SFT vectors -- one for each IFO in a multi-IFO search */
typedef struct tagMultiSFTVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(SFTVector*, data, UINT4, length));
#endif /* SWIG */
  UINT4      length;  	/**< number of ifos */
  SFTVector  **data; 	/**< sftvector for each ifo */
} MultiSFTVector;


/** A collection of PSD vectors -- one for each IFO in a multi-IFO search */
typedef struct tagMultiPSDVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(PSDVector*, data, UINT4, length));
#endif /* SWIG */
  UINT4      length;  	/**< number of ifos */
  PSDVector  **data; 	/**< sftvector for each ifo */
} MultiPSDVector;

/** One noise-weight (number) per SFT (therefore indexed over IFOs and SFTs */
typedef struct tagMultiNoiseWeights {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(REAL8Vector*, data, UINT4, length));
#endif /* SWIG */
  UINT4 length;		/**< number of ifos */
  REAL8Vector **data;	/**< weights-vector for each SFTs */
  REAL8 Sinv_Tsft;	/**< normalization factor used: \f$\mathcal{S}^{-1}\,T_\mathrm{SFT}\f$ (using single-sided PSD!) */
} MultiNoiseWeights;

/** A collection of (multi-IFO) time-series */
typedef struct tagMultiREAL4TimeSeries {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(REAL4TimeSeries*, data, UINT4, length));
#endif /* SWIG */
  UINT4 length;			/**< number of ifos */
  REAL4TimeSeries **data;	/**< vector of REAL4 timeseries */
} MultiREAL4TimeSeries;

/** A vector of 'timestamps' of type LIGOTimeGPS */
typedef struct tagLIGOTimeGPSVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(LIGOTimeGPS, data, UINT4, length));
#endif /* SWIG */
  UINT4 	length;		/**< number of timestamps */
  LIGOTimeGPS 	*data;		/**< array of timestamps */
  REAL8		deltaT;		/**< 'length' of each timestamp (e.g. typically Tsft) */
} LIGOTimeGPSVector;

/** A vector of 'timestamps' of type LIGOTimeGPS */
typedef struct tagMultiLIGOTimeGPSVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(LIGOTimeGPSVector*, data, UINT4, length));
#endif /* SWIG */
  UINT4 	        length;	   /**< number of timestamps vectors or ifos */
  LIGOTimeGPSVector 	**data;    /**< timestamps vector for each ifo */
} MultiLIGOTimeGPSVector;

/*---------- Global variables ----------*/
/* empty init-structs for the types defined in here */
extern const SFTtype empty_SFTtype;
extern const SFTVector empty_SFTVector;
extern const PSDVector empty_PSDVector;
extern const MultiSFTVector empty_MultiSFTVector;
extern const MultiPSDVector empty_MultiPSDVector;
extern const MultiNoiseWeights empty_MultiNoiseWeights;
extern const MultiREAL4TimeSeries empty_MultiREAL4TimeSeries;
extern const LIGOTimeGPSVector empty_LIGOTimeGPSVector;
extern const MultiLIGOTimeGPSVector empty_MultiLIGOTimeGPSVector;

// ---------- obsolete LAL-API was moved into external file
#include "SFTutils-LAL.h"
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

SFTVector* XLALExtractBandfromSFTs ( const SFTVector *sfts, REAL8 fMin, REAL8 fMax );

LIGOTimeGPSVector *XLALCreateTimestampVector (UINT4 len);
LIGOTimeGPSVector *XLALMakeTimestamps ( LIGOTimeGPS tStart, REAL8 duration, REAL8 tStep );
LIGOTimeGPSVector *XLALExtractTimestampsFromSFTs ( const SFTVector *sfts );
MultiLIGOTimeGPSVector *XLALExtractMultiTimestampsFromSFTs ( const MultiSFTVector *multiSFTs );

void XLALDestroyTimestampVector (LIGOTimeGPSVector *vect);
void XLALDestroyMultiTimestamps ( MultiLIGOTimeGPSVector *multiTS );

CHAR *XLALGetChannelPrefix ( const CHAR *name );
LALDetector *XLALGetSiteInfo ( const CHAR *name );

LALSegList *XLALReadSegmentsFromFile ( const char *fname );
void XLALDestroyPSDVector ( PSDVector *vect );
void XLALDestroyMultiSFTVector ( MultiSFTVector *multvect );
void XLALDestroyMultiPSDVector ( MultiPSDVector *multvect );
void XLALDestroyMultiNoiseWeights ( MultiNoiseWeights *weights );

/*@}*/

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
