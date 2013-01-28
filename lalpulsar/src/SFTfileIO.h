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
#include <lal/PulsarDataTypes.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

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
 *
 *
 <p> <h3> Overview:</h3>
 - SFT-reading: LALSFTdataFind(), LALLoadSFTs(), LALLoadMultiSFTs()
 - SFT-writing: LALWriteSFT2file(), LALWrite_v2SFT_to_v1file()
 - SFT-checking: LALCheckSFTs(): complete check of SFT-validity including CRC64 checksum
 - free SFT-catalog: LALDestroySFTCatalog()
 - general manipulation of SFTVectors:
 	- LALDestroySFTVector(): free up a complete SFT-vector
	- LALDestroyMultiSFTVector(): free a multi-IFO vector of SFT-vectors
	- LALConcatSFTVectors(): concatenate two ::SFTVector's
	- LALAppendSFT2Vector(): append a single SFT (::SFTtype) to an ::SFTVector

 <p>
 <h2>Usage: Reading of SFT-files</h2>

 The basic operation of <b>reading SFTs</b> from files proceeds in two simple steps:

 	-# LALSFTdataFind(): get an '::SFTCatalog' of SFTs matching certain requirements (free with LALDestroySFTCatalog())
	-# LALLoadSFTs(): load a frequency-band into a single-IFO SFTVector defined by the catalogue, OR <br>
	LALLoadMultiSFTs(): load a frequency-band into a multi-IFO vector of SFTVectors defined by the catalogue

 <b>Note 1:</b> currently supported SFT file-formats are (merged or single) SFT-v1 and SFT-v2 files.
 This might be extended in the future to support further file-formats (frames?).
 None of the following API depends on the details of the underlying file-format. This will ensure that
 codes using the following functions will NOT have to be changed irrespective of SFT file-format used.

 <b>Note 2:</b> irrespective of the underlying SFT file-format, the returned SFTs (::SFTVector) will
 <em>ALWAYS</em> be normalized according the the LAL-specification for frequency-series
 (<tt>LIGO-T010095-00</tt>), that is the pure DFT of the time-series \f$x_j\f$ is <em>multiplied</em>
 by the time-step \f$\Delta t\f$:
 \f[
 \mathrm{data}[k] = X^\mathrm{d}_k = \Delta t \,\sum_{j=0}^{N-1} x_j \,e^{-i2\pi \,k \,j / N}
 \f]

 <h4>Details to 1: find matching SFTs and get the SFTCatalog:</h4>

 \code
 LALSFTdataFind(LALStatus *, SFTCatalog **catalog, const CHAR *file_pattern, SFTConstraints *constraints);
 \endcode

 This function returns an SFTCatalog of matching SFTs for a given file-pattern
 (e.g. "SFT.*", "SFT.000", "/some/path/some_files_[0-9]?.sft", etc ) and additional, optional SFTConstraints.

 The optional constraints are:
 - detector-prefix (e.g. "H1", "H2", "L1", "G1", "V1", etc..)  [\em required for v1-SFTs!]
 - GPS start-time + end-time
 - a list of GPS-timestamps

 <b>Note 1:</b> Any constraint can be specified as \c NULL, all given constraints will be
 combined by logical \c AND.

 <b>Note 2:</b> if a timestamps-list is given, *ALL* timestamps within
 <tt>[startTime, endTime]</tt> MUST be found!]

 <b>Note 3:</b> LALSFTdataFind() will refuse to return any SFTs without their detector-name
 properly set. This applies only to v1-SFTs, for which you have to use constraints->detector,
 so that the detector-name gets properly set.

 <b>Note 4:</b> One special constraint->detector is "??", which acts as if
 constraints->detector==NULL, except that it allows v1-SFTs to be returned with
 detector-name set to "??"'.


 The returned SFTCatalog is a vector of 'SFTDescriptor's describing one SFT, with the fields
 - \c locator:  an opaque data-type describing where to read this SFT from.
 - \c header:	the SFts header
 - \c comment: the comment-string found in the SFT, if any
 - \c numBins: the number of frequency-bins in the SFT
 - \c version: version-number of SFT file-format
 - \c crc64: the crc64 checksum reported by this SFT


 One can use the following catalog-handling API functions:
 - LALDestroySFTCatalog(): free up a complete SFT-catalog
 - LALSFTtimestampsFromCatalog(): extract the list of SFT timestamps found in the ::SFTCatalog
 - LALDestroyTimestampVector(): free up a timestamps-vector (::LIGOTimeGPSVector)
 - XLALshowSFTLocator(): [*debugging only*] show a static string describing the 'locator'


 <b>NOTE:</b> The SFTs in the returned catalogue are \em guaranteed to
 - be sorted in order of increasing GPS-epoch
 - contain a valid detector-name, except if constraints->detector=="??"


 <h4>Details to 2: load frequency-band from SFTs described in an SFTCatalog</h4>

 \code
 LALLoadSFTs ( LALStatus *, SFTVector **sfts, const SFTCatalog *catalog, REAL8 fMin, REAL8 fMax);
 \endcode

 This function takes an ::SFTCatalog and reads the smallest frequency-band containing <tt>[fMin, fMax]</tt>
 from the SFTs, returning the resulting ::SFTVector. Note that this function will return an error if the
 SFTCatalog contains SFTs from different detectors, for which LALLoadMultiSFTs() must be used.

 The frequency-bounds are optional and \c -1 can be used to specify an 'open bound', i.e.<br>
 <tt>[-1, fMax]</tt>: read from first frequency-bin in the SFT up to \c fMax.<br>
 <tt>[fMin, -1]</tt>: read from \c fMin up to last frequency-bin in the SFTS<br>
 <tt>[-1, -1]</tt>: read ALL frequency-bins from SFT.

 \code
 LALLoadMultiSFTs ( LALStatus *, MultiSFTVector **sfts, const SFTCatalog *catalog, REAL8 fMin, REAL8 fMax);
 \endcode

 This function is similar to the above, except that it accepts an ::SFTCatalog with different detectors,
 and returns corresponding multi-IFO vector of SFTVectors.

 <p><h2>Usage: Writing of SFT-files</h2>

 For <b>writing SFTs</b> there are two functions, depending on the desired output-format (v1 or v2  SFTs):
 	- LALWriteSFT2file(): write a single SFT (::SFTtype) into an SFT-file following the specification v2
	(<tt>LIGO-T040164-01-Z</tt>).

	- LALWrite_v2SFT_to_v1file(): write a single ::SFTtype into an SFT-v1 file. Note: this is provided
	for backwards-compatibility, and assumes the input SFT-data in memory to be correctly normalized
	according to the v2-specification (i.e. data = dt x DFT).


Note: in addition to these two function which take properly normalized SFTs as input, there is a DEPRECATED
legacy-function, LALWriteSFTfile(), which writes an v1-SFT file, but *without* changing the data-normalization,
i.e. this will only be correct for v1-normalized data (i.e. data = DFT)

*/

/*@{*/

/** \name Error codes */
/*@{*/
#define SFTFILEIO_ENULL 	1
#define SFTFILEIO_EFILE 	2
#define SFTFILEIO_EHEADER 	3
#define SFTFILEIO_EVERSION 	4
#define SFTFILEIO_EVAL 		5

#define SFTFILEIO_EDIFFTSFT	6
#define SFTFILEIO_EDIFFDET	7
#define SFTFILEIO_EDETECTOR	8

#define SFTFILEIO_ENONULL 	12
#define SFTFILEIO_EFREQBAND 	13
#define SFTFILEIO_EMEM 		14
#define SFTFILEIO_EGLOB 	15
#define SFTFILEIO_EDIFFLENGTH 	17
#define SFTFILEIO_ESFTFORMAT	18
#define SFTFILEIO_ESFTWRITE	19
#define SFTFILEIO_ECONSTRAINTS  20
#define SFTFILEIO_EMERGEDSFT  	21
#define SFTFILEIO_ECRC64  	22

#define SFTFILEIO_MSGENULL 	"Null pointer"
#define SFTFILEIO_MSGEFILE 	"Error in file-IO"
#define SFTFILEIO_MSGEHEADER 	"Incorrect header in file"
#define SFTFILEIO_MSGEVERSION 	"This SFT-version is not currently supported"
#define SFTFILEIO_MSGEVAL  	"Invalid value"

#define SFTFILEIO_MSGEDIFFTSFT	"Inconsistent values of Tsft for matched SFTs"
#define SFTFILEIO_MSGEDIFFDET	"Inconsistent detector-values for matched SFTs"
#define SFTFILEIO_MSGEDETECTOR	"Illegal detector name"

#define SFTFILEIO_MSGENONULL  	"Output pointer not NULL"
#define SFTFILEIO_MSGEFREQBAND 	"Required frequency-band is not in SFT"
#define SFTFILEIO_MSGEMEM 	"Out of memory"
#define SFTFILEIO_MSGEGLOB 	"Failed to get filelist from directory/pattern"
#define SFTFILEIO_MSGEDIFFLENGTH "Sorry, can only read SFTs of identical length (currently)"
#define SFTFILEIO_MSGESFTFORMAT	 "Illegal SFT-format"
#define SFTFILEIO_MSGESFTWRITE	 "Failed to write SFT to file"
#define SFTFILEIO_MSGECONSTRAINTS "Could not satisfy the requested SFT-query constraints"
#define SFTFILEIO_MSGEMERGEDSFT   "Inconsistent blocks in merged SFT"
#define SFTFILEIO_MSGECRC64	"Invalid CRC64 checksum in SFT"
/*@}*/

/** 'Constraints' for SFT-matching: which detector, within which time-stretch and which
 * timestamps exactly should be loaded ?
 * Any of the entries is optional, and they will be combined by logical AND.
 * Note however, that *ALL* timestamps within [startTime, endTime] MUST be found if specified.
 */
typedef struct tagSFTConstraints
{
  CHAR *detector;			/**< 2-char channel-prefix describing the detector (eg 'H1', 'H2', 'L1', 'G1' etc) */
  LIGOTimeGPS *startTime;		/**< only include SFTs starting >= startTime */
  LIGOTimeGPS *endTime;			/**< only include SFTs starting <= endTime */
  LIGOTimeGPSVector *timestamps;	/**< list of timestamps  */
} SFTConstraints;


/** A 'descriptor' of an SFT: basically containing the header-info plus an opaque description
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
  UINT4 length;			/**< number of SFTs in catalog */
  SFTDescriptor *data;		/**< array of data-entries describing matched SFTs */
} SFTCatalog;

/*---------- Global variables ----------*/
/* empty init-structs for the types defined in here */
extern const SFTConstraints empty_SFTConstraints;
extern const SFTCatalog empty_SFTCatalog;

/*
 * Functions Declarations (i.e., prototypes).
 */

/*================================================================================
 * NEW API: allowing for SFT-v2
 *================================================================================*/
SFTCatalog *XLALSFTdataFind ( const CHAR *file_pattern, const SFTConstraints *constraints );

int  XLALWriteSFT2fp (const SFTtype *sft, FILE *fp, const CHAR *SFTcomment );
int  XLALWriteSFT2file (const SFTtype *sft, const CHAR *fname, const CHAR *SFTcomment );
int  XLALWriteSFTVector2Dir (const SFTVector *sftVect, const CHAR *basename, const CHAR *SFTcomment, const CHAR *description);
int  XLALWriteSFTVector2File(const SFTVector *sftVect, const CHAR *filename, const CHAR *SFTcomment);
LIGOTimeGPSVector *XLALReadTimestampsFile ( const CHAR *fname );

SFTVector* XLALLoadSFTs (const SFTCatalog *catalog, REAL8 fMin, REAL8 fMax);
MultiSFTVector* XLALLoadMultiSFTs (const SFTCatalog *catalog, REAL8 fMin, REAL8 fMax);
void XLALDestroySFTCatalog ( SFTCatalog *catalog );
INT4 XLALCountIFOsInCatalog( const SFTCatalog *catalog);
const CHAR * XLALshowSFTLocator ( const struct tagSFTLocator *locator );

/*================================================================================
 * DEPRECATED LAL-API [use XLAL-API whenever available]
 *================================================================================*/
void LALReadTimestampsFile (LALStatus* , LIGOTimeGPSVector **timestamps, const CHAR *fname); /* use XLALReadTimestampsFile() instead! */
void LALSFTdataFind (LALStatus *, SFTCatalog **catalog, const CHAR *file_pattern, SFTConstraints *constraints);
void LALLoadSFTs ( LALStatus *, SFTVector **sfts, const SFTCatalog *catalog, REAL8 fMin, REAL8 fMax);
void LALLoadMultiSFTs ( LALStatus *status, MultiSFTVector **sfts, const SFTCatalog *catalog, REAL8 fMin, REAL8 fMax);
void LALWriteSFT2file (LALStatus *, const SFTtype *sft, const CHAR *fname, const CHAR *SFTcomment );
void LALWriteSFTVector2Dir (LALStatus *, const SFTVector *sftVect, const CHAR *basename, const CHAR *SFTcomment, const CHAR *description);

void LALWrite_v2SFT_to_v1file (LALStatus *, const SFTtype *sft, const CHAR *fname);
void LALCheckSFTs ( LALStatus *, INT4 *check_result, const CHAR *file_pattern, SFTConstraints *constraints );
void LALCheckSFTCatalog ( LALStatus *status, INT4 *check_result, SFTCatalog *catalog );

void LALDestroySFTCatalog ( LALStatus *status, SFTCatalog **catalog );
void LALSFTtimestampsFromCatalog (LALStatus *, LIGOTimeGPSVector **timestamps, const SFTCatalog *catalog );

/*================================================================================
 * OBSOLETE v1-only API [DEPRECATED!]
 *================================================================================*/

/** [DEPRECATED] This structure contains the header-info contained in an SFT-file of specification
 * version v1.0.
 */
typedef struct tagSFTHeader {
  REAL8  version;		/**< SFT version-number (currently only 1.0 allowed )*/
  INT4   gpsSeconds;		/**< gps start-time (seconds)*/
  INT4   gpsNanoSeconds;	/**< gps start-time (nanoseconds) */
  REAL8  timeBase;		/**< length of data-stretch in seconds */
  INT4   fminBinIndex;		/**< first frequency-index contained in SFT */
  INT4   length;  		/**< number of frequency bins */
} SFTHeader;



void LALReadSFTheader (LALStatus *, SFTHeader *header, const CHAR *fname);
void LALReadSFTdata (LALStatus *, SFTtype *sft, const CHAR *fname, INT4 fminBinIndex);
void LALWriteSFTfile (LALStatus *, const SFTtype *sft, const CHAR *outfname);
void LALReadSFTfile (LALStatus *, SFTtype **sft, REAL8 fMin, REAL8 fMax, const CHAR *fname);

void LALReadSFTfiles (LALStatus *,
		      SFTVector **sftvect,
		      REAL8 fMin,
		      REAL8 fMax,
		      UINT4 wingBins,
		      const CHAR *fpattern);

void
LALGetSFTheaders (LALStatus *,
		  SFTVector **headers,
		  const CHAR *fpattern,
		  const LIGOTimeGPS *startTime,
		  const LIGOTimeGPS *endTime);


/*@}*/

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _SFTBIN_H */
