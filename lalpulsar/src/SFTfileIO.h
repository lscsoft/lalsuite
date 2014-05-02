/*
 * Copyright (C) 2010 Karl Wette
 * Copyright (C) 2004, 2005 R. Prix, B. Machenschalk, A.M. Sintes
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
#ifndef _SFTFILEIO_H  	/* Double-include protection. */
#define _SFTFILEIO_H

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

/* includes */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>

/**
 * \defgroup SFTfileIO_h Header SFTfileIO.h
 * \ingroup pkg_SFTIO
 * \author R. Prix, B. Machenschalk, A.M. Sintes, B. Krishnan
 * \brief Module for reading/writing/manipulating SFTs (Short Fourier transforms)
 *
 * This implements the SFTv2 standard defined in LIGO-T040164-01-Z
 * A previous non-LAL implementation of this standard is found in the "SFT reference library"
 * gravity.phys.uwm.edu:2402/usr/local/cvs/lscsoft sftlib, Copyright (C) 2004 Bruce Allen
 *
 * <p> <h3> Overview:</h3>
 * - SFT-reading: LALSFTdataFind(), LALLoadSFTs(), LALLoadMultiSFTs()
 * - SFT-writing: LALWriteSFT2file(), LALWrite_v2SFT_to_v1file()
 * - SFT-checking: LALCheckSFTs(): complete check of SFT-validity including CRC64 checksum
 * - free SFT-catalog: LALDestroySFTCatalog()
 * - general manipulation of SFTVectors:
 * - LALDestroySFTVector(): free up a complete SFT-vector
 * - LALDestroyMultiSFTVector(): free a multi-IFO vector of SFT-vectors
 * - LALConcatSFTVectors(): concatenate two ::SFTVector's
 * - LALAppendSFT2Vector(): append a single SFT (::SFTtype) to an ::SFTVector
 *
 * <p>
 * <h2>Usage: Reading of SFT-files</h2>
 *
 * The basic operation of <b>reading SFTs</b> from files proceeds in two simple steps:
 *
 * -# LALSFTdataFind(): get an '::SFTCatalog' of SFTs matching certain requirements (free with LALDestroySFTCatalog())
 * -# LALLoadSFTs(): load a frequency-band into a single-IFO SFTVector defined by the catalogue, OR <br>
 * LALLoadMultiSFTs(): load a frequency-band into a multi-IFO vector of SFTVectors defined by the catalogue
 *
 * <b>Note 1:</b> currently supported SFT file-formats are (merged or single) SFT-v1 and SFT-v2 files.
 * This might be extended in the future to support further file-formats (frames?).
 * None of the following API depends on the details of the underlying file-format. This will ensure that
 * codes using the following functions will NOT have to be changed irrespective of SFT file-format used.
 *
 * <b>Note 2:</b> irrespective of the underlying SFT file-format, the returned SFTs (::SFTVector) will
 * <em>ALWAYS</em> be normalized according the the LAL-specification for frequency-series
 * (<tt>LIGO-T010095-00</tt>), that is the pure DFT of the time-series \f$x_j\f$ is <em>multiplied</em>
 * by the time-step \f$\Delta t\f$:
 * \f[
 * \mathrm{data}[k] = X^\mathrm{d}_k = \Delta t \,\sum_{j=0}^{N-1} x_j \,e^{-i2\pi \,k \,j / N}
 * \f]
 *
 * <h4>Details to 1: find matching SFTs and get the SFTCatalog:</h4>
 *
 * \code
 * LALSFTdataFind(LALStatus *, SFTCatalog **catalog, const CHAR *file_pattern, SFTConstraints *constraints);
 * \endcode
 *
 * This function returns an SFTCatalog of matching SFTs for a given file-pattern
 * (e.g. "SFT.*", "SFT.000", "/some/path/some_files_[0-9]?.sft", etc ) and additional, optional SFTConstraints.
 *
 * The optional constraints are:
 * - detector-prefix (e.g. "H1", "H2", "L1", "G1", "V1", etc..)  [\em required for v1-SFTs!]
 * - GPS start-time + end-time
 * - a list of GPS-timestamps
 *
 * <b>Note 1:</b> Any constraint can be specified as \c NULL, all given constraints will be
 * combined by logical \c AND.
 *
 * <b>Note 2:</b> if a timestamps-list is given, *ALL* timestamps within
 * <tt>[minStartTime, maxStartTime)</tt> MUST be found!]
 *
 * <b>Note 3:</b> LALSFTdataFind() will refuse to return any SFTs without their detector-name
 * properly set. This applies only to v1-SFTs, for which you have to use constraints->detector,
 * so that the detector-name gets properly set.
 *
 * <b>Note 4:</b> One special constraint->detector is "??", which acts as if
 * constraints->detector==NULL, except that it allows v1-SFTs to be returned with
 * detector-name set to "??"'.
 *
 * The returned SFTCatalog is a vector of 'SFTDescriptor's describing one SFT, with the fields
 * - \c locator:  an opaque data-type describing where to read this SFT from.
 * - \c header:	the SFts header
 * - \c comment: the comment-string found in the SFT, if any
 * - \c numBins: the number of frequency-bins in the SFT
 * - \c version: version-number of SFT file-format
 * - \c crc64: the crc64 checksum reported by this SFT
 *
 * One can use the following catalog-handling API functions:
 * - LALDestroySFTCatalog(): free up a complete SFT-catalog
 * - LALSFTtimestampsFromCatalog(): extract the list of SFT timestamps found in the ::SFTCatalog
 * - LALDestroyTimestampVector(): free up a timestamps-vector (::LIGOTimeGPSVector)
 * - XLALshowSFTLocator(): [*debugging only*] show a static string describing the 'locator'
 *
 * <b>NOTE:</b> The SFTs in the returned catalogue are \em guaranteed to
 * - be sorted in order of increasing GPS-epoch
 * - contain a valid detector-name, except if constraints->detector=="??"
 *
 * <h4>Details to 2: load frequency-band from SFTs described in an SFTCatalog</h4>
 *
 * \code
 * LALLoadSFTs ( LALStatus *, SFTVector **sfts, const SFTCatalog *catalog, REAL8 fMin, REAL8 fMax);
 * \endcode
 *
 * This function takes an ::SFTCatalog and reads the smallest frequency-band containing <tt>[fMin, fMax]</tt>
 * from the SFTs, returning the resulting ::SFTVector. Note that this function will return an error if the
 * SFTCatalog contains SFTs from different detectors, for which LALLoadMultiSFTs() must be used.
 *
 * The frequency-bounds are optional and \c -1 can be used to specify an 'open bound', i.e.<br>
 * <tt>[-1, fMax]</tt>: read from first frequency-bin in the SFT up to \c fMax.<br>
 * <tt>[fMin, -1]</tt>: read from \c fMin up to last frequency-bin in the SFTS<br>
 * <tt>[-1, -1]</tt>: read ALL frequency-bins from SFT.
 *
 * \code
 * LALLoadMultiSFTs ( LALStatus *, MultiSFTVector **sfts, const SFTCatalog *catalog, REAL8 fMin, REAL8 fMax);
 * \endcode
 *
 * This function is similar to the above, except that it accepts an ::SFTCatalog with different detectors,
 * and returns corresponding multi-IFO vector of SFTVectors.
 *
 * <p><h2>Usage: Writing of SFT-files</h2>
 *
 * For <b>writing SFTs</b> there are two functions, depending on the desired output-format (v1 or v2  SFTs):
 * - LALWriteSFT2file(): write a single SFT (::SFTtype) into an SFT-file following the specification v2
 * (<tt>LIGO-T040164-01-Z</tt>).
 *
 * - LALWrite_v2SFT_to_v1file(): write a single ::SFTtype into an SFT-v1 file. Note: this is provided
 * for backwards-compatibility, and assumes the input SFT-data in memory to be correctly normalized
 * according to the v2-specification (i.e. data = dt x DFT).
 *
 * Note: in addition to these two function which take properly normalized SFTs as input, there is a DEPRECATED
 * legacy-function, LALWriteSFTfile(), which writes an v1-SFT file, but *without* changing the data-normalization,
 * i.e. this will only be correct for v1-normalized data (i.e. data = DFT)
 *
 */

/*@{*/

// ---------- exported types ----------

/** A vector of COMPLEX8FrequencySeries */
typedef struct tagCOMPLEX8FrequencySeriesVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(COMPLEX8FrequencySeriesVector, COMPLEX8FrequencySeries, data, UINT4, length));
#endif /* SWIG */
  UINT4 			length;		/**< number of SFTs */
  COMPLEX8FrequencySeries 	*data;		/**< array of SFTs */
} COMPLEX8FrequencySeriesVector;


/** A vector of REAL4FrequencySeries */
typedef struct tagREAL4FrequencySeriesVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(REAL4FrequencySeriesVector, REAL4FrequencySeries, data, UINT4, length));
#endif /* SWIG */
  UINT4                  length;
  REAL4FrequencySeries   *data;
} REAL4FrequencySeriesVector;

/** A collection of (multi-IFO) time-series */
typedef struct tagMultiREAL4TimeSeries {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(MultiREAL4TimeSeries, REAL4TimeSeries*, data, UINT4, length));
#endif /* SWIG */
  UINT4 length;			/**< number of ifos */
  REAL4TimeSeries **data;	/**< vector of REAL4 timeseries */
} MultiREAL4TimeSeries;


/** A vector of 'timestamps' of type LIGOTimeGPS */
typedef struct tagLIGOTimeGPSVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(LIGOTimeGPSVector, LIGOTimeGPS, data, UINT4, length));
#endif /* SWIG */
  UINT4 	length;		/**< number of timestamps */
  LIGOTimeGPS 	*data;		/**< array of timestamps */
  REAL8		deltaT;		/**< 'length' of each timestamp (e.g. typically Tsft) */
} LIGOTimeGPSVector;

/** A vector of 'timestamps' of type LIGOTimeGPS */
typedef struct tagMultiLIGOTimeGPSVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(MultiLIGOTimeGPSVector, LIGOTimeGPSVector*, data, UINT4, length));
#endif /* SWIG */
  UINT4 	        length;	   /**< number of timestamps vectors or ifos */
  LIGOTimeGPSVector 	**data;    /**< timestamps vector for each ifo */
} MultiLIGOTimeGPSVector;


/** A so-called 'SFT' (short-Fourier-transform) will be stored in a COMPLEX8FrequencySeries */
typedef COMPLEX8FrequencySeries 	SFTtype;


/** The corresponding vector-type to hold a vector of 'SFTs' */
typedef COMPLEX8FrequencySeriesVector 	SFTVector;

/** A collection of SFT vectors -- one for each IFO in a multi-IFO search */
typedef struct tagMultiSFTVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(MultiSFTVector, SFTVector*, data, UINT4, length));
#endif /* SWIG */
  UINT4      length;  	/**< number of ifos */
  SFTVector  **data; 	/**< sftvector for each ifo */
} MultiSFTVector;


/**
 * 'Constraints' for SFT-matching: which detector, within which time-stretch and which
 * timestamps exactly should be loaded ?
 * Any of the entries is optional, and they will be combined by logical AND.
 * Note however, that *ALL* timestamps within [minStartTime, maxStartTime) MUST be found if specified.
 */
typedef struct tagSFTConstraints
{
  CHAR *detector;			/**< 2-char channel-prefix describing the detector (eg 'H1', 'H2', 'L1', 'G1' etc) */
  LIGOTimeGPS *minStartTime;		/**< only include SFTs whose epoch is >= minStartTime */
  LIGOTimeGPS *maxStartTime;		/**< only include SFTs whose epoch is <  maxStartTime */
  LIGOTimeGPSVector *timestamps;	/**< list of timestamps  */
} SFTConstraints;


/**
 * A 'descriptor' of an SFT: basically containing the header-info plus an opaque description
 * of where exactly to load this SFT from.
 */
typedef struct tagSFTDescriptor
{
  struct tagSFTLocator *locator; 	/**< *internal* description of where to find this SFT [opaque!] */
  SFTtype header;			/**< SFT-header info */
  CHAR *comment;			/**< comment-entry in SFT-header (v2 only) */
  UINT4 numBins;			/**< number of frequency-bins in this SFT */
  UINT4 version;			/**< SFT-specification version */
  UINT8 crc64;				/**< crc64 checksum */
} SFTDescriptor;


/** An "SFT-catalogue": a vector of SFTdescriptors, as returned by LALSFTdataFind() */
typedef struct tagSFTCatalog
{
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(SFTCatalog, SFTDescriptor, data, UINT4, length));
#endif /* SWIG */
  UINT4 length;			/**< number of SFTs in catalog */
  SFTDescriptor *data;		/**< array of data-entries describing matched SFTs */
} SFTCatalog;

/**
 * A multi-SFT-catalogue "view": a multi-IFO vector of SFT-catalogs
 *
 * Note: this is only a multi-IFO "view" of an existing SFTCatalog,
 * various allocated memory of the original catalog is only
 * pointed to, not duplicated!
 * This means one must not free the original catalog
 * while this multi-view is still in use!
 */
typedef struct tagMultiSFTCatalogView
{
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(MultiSFTCatalogView, SFTCatalog, data, UINT4, length));
#endif /* SWIG */
  UINT4 length;			/**< number of detectors */
  SFTCatalog *data;		/**< array of SFT-catalog pointers */
} MultiSFTCatalogView;


/*---------- Global variables ----------*/
/* empty init-structs for the types defined in here */
extern const SFTConstraints empty_SFTConstraints;
extern const SFTCatalog empty_SFTCatalog;
extern const SFTtype empty_SFTtype;
extern const SFTVector empty_SFTVector;
extern const MultiSFTVector empty_MultiSFTVector;
extern const MultiREAL4TimeSeries empty_MultiREAL4TimeSeries;
extern const LIGOTimeGPSVector empty_LIGOTimeGPSVector;
extern const MultiLIGOTimeGPSVector empty_MultiLIGOTimeGPSVector;

/*
 * Functions Declarations (i.e., prototypes).
 */

int XLALCWGPSinRange( const LIGOTimeGPS gps, const LIGOTimeGPS* minGPS, const LIGOTimeGPS* maxGPS );

LALStringVector *XLALFindFiles (const CHAR *globstring);

SFTCatalog *XLALSFTdataFind ( const CHAR *file_pattern, const SFTConstraints *constraints );

int XLALWriteSFTVector2Dir  ( const SFTVector *sftVect, const CHAR *dirname, const CHAR *SFTcomment, const CHAR *Misc );
int XLALWriteSFTVector2File ( const SFTVector *sftVect, const CHAR *dirname, const CHAR *SFTcomment, const CHAR *Misc );
int XLALWriteSFTVector2NamedFile ( const SFTVector *sftVect, const CHAR *filename, const CHAR *SFTcomment );
int XLALWriteSFT2fp   ( const SFTtype *sft, FILE *fp, const CHAR *SFTcomment );
int XLALWriteSFT2file ( const SFTtype *sft, const CHAR *fname, const CHAR *SFTcomment );

LIGOTimeGPSVector *XLALReadTimestampsFile ( const CHAR *fname );
MultiLIGOTimeGPSVector *XLALReadMultiTimestampsFiles ( const LALStringVector *fnames );

SFTVector* XLALLoadSFTs (const SFTCatalog *catalog, REAL8 fMin, REAL8 fMax);

MultiSFTVector* XLALLoadMultiSFTs (const SFTCatalog *catalog, REAL8 fMin, REAL8 fMax);
MultiSFTVector *XLALLoadMultiSFTsFromView ( const MultiSFTCatalogView *multiCatalogView, REAL8 fMin, REAL8 fMax );

void XLALDestroySFTCatalog ( SFTCatalog *catalog );
INT4 XLALCountIFOsInCatalog( const SFTCatalog *catalog);
const CHAR * XLALshowSFTLocator ( const struct tagSFTLocator *locator );

void XLALDestroyMultiSFTCatalogView ( MultiSFTCatalogView *multiView );
MultiSFTCatalogView *XLALGetMultiSFTCatalogView ( const SFTCatalog *catalog );

char *XLALGetOfficialName4SFT ( const SFTtype *sft, const char *Misc );
char *XLALGetOfficialName4MergedSFTs ( const SFTVector *sfts, const char *Misc );
char *XLALOfficialSFTFilename ( char site, char channel, UINT4 numSFTs, UINT4 Tsft, UINT4 GPS_start, UINT4 Tspan, const char *Misc );
int XLALCheckValidDescriptionField ( const char *desc );

/*@}*/

// ---------- obsolete LAL-API was moved into external file
#include <lal/SFTfileIO-LAL.h>
// ------------------------------

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _SFTBIN_H */
