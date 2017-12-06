#ifdef _SFTFILEIO_H /* Only include via SFTfileIO.h */
#ifndef SWIG /* exclude from SWIG interface */
#ifndef _SFTFILEIOLAL_H
#define _SFTFILEIOLAL_H
/** \cond DONT_DOXYGEN */

/** \name Error codes */
/*@{*/
#define SFTFILEIO_ENULL         1
#define SFTFILEIO_EFILE         2
#define SFTFILEIO_EHEADER       3
#define SFTFILEIO_EVERSION      4
#define SFTFILEIO_EVAL          5

#define SFTFILEIO_EDIFFTSFT	6
#define SFTFILEIO_EDIFFDET	7
#define SFTFILEIO_EDETECTOR	8

#define SFTFILEIO_ENONULL       12
#define SFTFILEIO_EFREQBAND     13
#define SFTFILEIO_EMEM          14
#define SFTFILEIO_EGLOB         15
#define SFTFILEIO_EDIFFLENGTH   17
#define SFTFILEIO_ESFTFORMAT	18
#define SFTFILEIO_ESFTWRITE	19
#define SFTFILEIO_ECONSTRAINTS  20
#define SFTFILEIO_EMERGEDSFT    21
#define SFTFILEIO_ECRC64        22

#define SFTFILEIO_MSGENULL      "Null pointer"
#define SFTFILEIO_MSGEFILE      "Error in file-IO"
#define SFTFILEIO_MSGEHEADER    "Incorrect header in file"
#define SFTFILEIO_MSGEVERSION   "This SFT-version is not currently supported"
#define SFTFILEIO_MSGEVAL       "Invalid value"

#define SFTFILEIO_MSGEDIFFTSFT	"Inconsistent values of Tsft for matched SFTs"
#define SFTFILEIO_MSGEDIFFDET	"Inconsistent detector-values for matched SFTs"
#define SFTFILEIO_MSGEDETECTOR	"Illegal detector name"

#define SFTFILEIO_MSGENONULL    "Output pointer not NULL"
#define SFTFILEIO_MSGEFREQBAND  "Required frequency-band is not in SFT"
#define SFTFILEIO_MSGEMEM       "Out of memory"
#define SFTFILEIO_MSGEGLOB      "Failed to get filelist from directory/pattern"
#define SFTFILEIO_MSGEDIFFLENGTH "Sorry, can only read SFTs of identical length (currently)"
#define SFTFILEIO_MSGESFTFORMAT  "Illegal SFT-format"
#define SFTFILEIO_MSGESFTWRITE   "Failed to write SFT to file"
#define SFTFILEIO_MSGECONSTRAINTS "Could not satisfy the requested SFT-query constraints"
#define SFTFILEIO_MSGEMERGEDSFT   "Inconsistent blocks in merged SFT"
#define SFTFILEIO_MSGECRC64	"Invalid CRC64 checksum in SFT"
/*@}*/

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

/**
 * [DEPRECATED] This structure contains the header-info contained in an SFT-file of specification
 * version v1.0.
 */
typedef struct tagSFTHeader {
  REAL8  version;		/**< SFT version-number (currently only 1.0 allowed )*/
  INT4   gpsSeconds;		/**< gps start-time (seconds)*/
  INT4   gpsNanoSeconds;	/**< gps start-time (nanoseconds) */
  REAL8  timeBase;		/**< length of data-stretch in seconds */
  INT4   fminBinIndex;		/**< first frequency-index contained in SFT */
  INT4   length;                /**< number of frequency bins */
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

/** \endcond */
#endif /* _SFTFILEIOLAL_H */
#endif /* SWIG */
#endif /* _SFTFILEIO_H */
