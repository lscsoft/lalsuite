/*
 * Copyright (C) 2010, 2013, 2014, 2016, 2019--2022 Karl Wette
 * Copyright (C) 2010 Chris Messenger
 * Copyright (C) 2009, 2011 Adam Mercer
 * Copyright (C) 2006 Badri Krishnan
 * Copyright (C) 2004--2007, 2010, 2012, 2017 Bernd Machenschalk
 * Copyright (C) 2004--2008, 2010--2016 Reinhard Prix
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

/*---------- includes ----------*/

#include "SFTinternal.h"
#include "SFTReferenceLibrary.h"

/*---------- constants ----------*/

/** blocksize used in SFT-reading for the CRC-checksum computation (has to be multiple of 8 !!) */
#define BLOCKSIZE 8192 * 8

/*---------- internal types ----------*/

typedef struct
{
  REAL8 version;
  INT4 gps_sec;
  INT4 gps_nsec;
  REAL8 tbase;
  INT4 first_frequency_index;
  INT4 nsamples;
  UINT8 crc64;
  CHAR detector[2];
  UINT2 windowspec;
  INT4 comment_length;
} _SFT_header_t;

/** segments read so far from one SFT */
typedef struct {
  UINT4 first;                     /**< first bin in this segment */
  UINT4 last;                      /**< last bin in this segment */
  LIGOTimeGPS epoch;               /**< timestamp of this SFT */
  struct tagSFTLocator *lastfrom;  /**< last bin read from this locator */
} SFTReadSegment;

/*---------- internal prototypes ----------*/

static int read_header_from_fp ( FILE *fp, SFTtype *header, UINT4 *nsamples, UINT8 *header_crc64, UINT8 *ref_crc64, UINT2 *SFTwindowspec, CHAR **SFTcomment, BOOLEAN swapEndian);

/*========== function definitions ==========*/

/// \addtogroup SFTfileIO_h
/// @{

/**
 * Load the given frequency-band <tt>[fMin, fMax)</tt> (half-open) from the SFT-files listed in the
 * SFT-'catalogue' ( returned by XLALSFTdataFind() ).
 *
 * Note: \a fMin (or \a fMax) is allowed to be set to \c -1, which means to read in all
 * Frequency-bins from the lowest (or up to the highest) found in all SFT-files of the catalog.
 *
 * Note 2: The returned frequency-interval is guaranteed to contain <tt>[fMin, fMax)</tt>,
 * but is allowed to be larger, as it must be an interval of discrete frequency-bins as found
 * in the SFT-file.
 *
 * Note 3: This function has the capability to read sequences of SFT segments and
 * putting them together to single SFTs while reading.
 *
 * Note 4: The 'fudge region' allowing for numerical noise is fudge= 10*LAL_REAL8_EPS ~2e-15
 * relative deviation: ie if the SFT contains a bin at 'fi', then we consider for example
 * "fMin == fi" if  fabs(fi - fMin)/fi < fudge.
 */
SFTVector*
XLALLoadSFTs (const SFTCatalog *catalog,   /**< The 'catalogue' of SFTs to load */
	      REAL8 fMin,		   /**< minumum requested frequency (-1 = read from lowest) */
	      REAL8 fMax		   /**< maximum requested frequency (-1 = read up to highest) */
	      )
{
  UINT4 catPos;                    /**< current file in catalog */
  UINT4 firstbin, lastbin;         /**< the first and last bin we want to read */
  UINT4 minbin, maxbin;            /**< min and max bin of all SFTs in the catalog */
  UINT4 nSFTs = 1;                 /**< number of SFTs, i.e. different GPS timestamps */
  REAL8 deltaF;                    /**< frequency spacing of SFT */
  SFTCatalog locatalog;            /**< local copy of the catalog to be sorted by 'locator' */
  SFTVector* sftVector = NULL;     /**< the vector of SFTs to be returned */
  SFTReadSegment*segments = NULL;  /**< array of segments already read of an SFT */
  char empty = '\0';               /**< empty string */
  char* fname = &empty;            /**< name of currently open file, initially "" */
  FILE* fp = NULL;                 /**< open file */
  SFTtype* thisSFT = NULL;         /**< SFT to read from file */

  /* error handler: free memory and return with error */
#define XLALLOADSFTSERROR(eno)	{		\
    if(fp)					\
      fclose(fp);				\
    if(segments) 				\
      XLALFree(segments);			\
    if(locatalog.data)				\
      XLALFree(locatalog.data);			\
    if(thisSFT)					\
      XLALDestroySFT(thisSFT);			\
    if(sftVector)				\
      XLALDestroySFTVector(sftVector);		\
    XLAL_ERROR_NULL(eno);	                \
  }

  /* initialize locatalog.data so it doesn't get free()d on early error */
  locatalog.data = NULL;

  /* check function parameters */
  if(!catalog)
    XLALLOADSFTSERROR(XLAL_EINVAL);

  /* determine number of SFTs, i.e. number of different GPS timestamps.
     The catalog should be sorted by GPS time, so just count changes.
     Record the 'index' of GPS time in the 'isft' field of the locator,
     so that we know later in which SFT to put this segment

     while at it, record max and min bin of all SFTs in the catalog */

  LIGOTimeGPS epoch = catalog->data[0].header.epoch;
  catalog->data[0].locator->isft = nSFTs - 1;
  deltaF = catalog->data[0].header.deltaF; /* Hz/bin */
  minbin = firstbin = lround ( catalog->data[0].header.f0 / deltaF );
  maxbin = lastbin = firstbin + catalog->data[0].numBins - 1;
  for(catPos = 1; catPos < catalog->length; catPos++) {
    firstbin = lround ( catalog->data[catPos].header.f0 / deltaF );
    lastbin = firstbin + catalog->data[catPos].numBins - 1;
    if (firstbin < minbin)
      minbin = firstbin;
    if (lastbin > maxbin)
      maxbin = lastbin;
    if(!GPSEQUAL(epoch, catalog->data[catPos].header.epoch)) {
      epoch = catalog->data[catPos].header.epoch;
      nSFTs++;
    }
    catalog->data[catPos].locator->isft = nSFTs - 1;
  }
  XLALPrintInfo("%s: fMin: %f, fMax: %f, deltaF: %f, minbin: %u, maxbin: %u\n", __func__, fMin, fMax, deltaF, minbin, maxbin);

  /* calculate first and last frequency bin to read */
  if (fMin < 0)
    firstbin = minbin;
  else
    firstbin = XLALRoundFrequencyDownToSFTBin ( fMin, deltaF );
  if (fMax < 0)
    lastbin = maxbin;
  else {
    lastbin = XLALRoundFrequencyUpToSFTBin ( fMax, deltaF ) - 1;
    if((lastbin == 0) && (fMax != 0)) {
      XLALPrintError("ERROR: last bin to read is 0 (fMax: %f, deltaF: %f)\n", fMax, deltaF);
      XLALLOADSFTSERROR(XLAL_EINVAL);
    }
  }
  XLALPrintInfo ( "%s: Reading from first bin: %u, last bin: %u\n", __func__, firstbin, lastbin);

  /* allocate the SFT vector that will be returned */
  if (!(sftVector = XLALCreateSFTVector (nSFTs, lastbin + 1 - firstbin))) {
    XLALPrintError("ERROR: Couldn't create sftVector\n");
    XLALLOADSFTSERROR(XLAL_ENOMEM);
  }

  /* allocate an additional single SFT where SFTs are read in */
  if(!(thisSFT = XLALCreateSFT (lastbin + 1 - firstbin))) {
    XLALPrintError("ERROR: Couldn't create thisSFT\n");
    XLALLOADSFTSERROR(XLAL_ENOMEM);
  }

  /* make a copy of the catalog that gets sorted by locator.
     Eases maintaing a correctly (epoch-)sorted catalog, particulary in case of errors
     note: only the pointers to the SFTdescriptors are copied & sorted, not the descriptors */
  locatalog.length = catalog->length;
  {
    UINT4 size = catalog->length * sizeof(catalog->data[0]);
    if(!(locatalog.data = XLALMalloc(size))) {
      XLALPrintError("ERROR: Couldn't allocate locatalog.data\n");
      XLALLOADSFTSERROR(XLAL_ENOMEM);
    }
    memcpy(locatalog.data, catalog->data, size);
  }

  /* sort catalog by f0, locator */
  qsort( (void*)locatalog.data, locatalog.length, sizeof( locatalog.data[0] ), compareSFTloc );

  /* allocate segment vector, one element per final SFT */
  if(!(segments = XLALCalloc(nSFTs, sizeof(SFTReadSegment)))) {
    XLALPrintError("ERROR: Couldn't allocate locatalog.data\n");
    XLALLOADSFTSERROR(XLAL_ENOMEM);
  }

  /* loop over all files (actually locators) in the catalog */
  for(catPos = 0; catPos < catalog->length; catPos++) {

    struct tagSFTLocator*locator = locatalog.data[catPos].locator;
    UINT4 isft = locator->isft;;
    UINT4 firstBinRead;
    UINT4 lastBinRead;

    if (locatalog.data[catPos].header.data) {
      /* the SFT data has already been read into the catalog,
	 copy the relevant part to thisSFT */

      volatile REAL8 tmp = locatalog.data[catPos].header.f0 / deltaF;
      UINT4 firstSFTbin = lround ( tmp );
      UINT4 lastSFTbin = firstSFTbin + locatalog.data[catPos].numBins - 1;
      UINT4 firstBin2read = firstbin;
      UINT4 lastBin2read = lastbin;
      UINT4 numBins2read, offsetBins;

      /* limit the interval to be read to what's actually in the SFT */
      if ( firstBin2read < firstSFTbin )
	firstBin2read = firstSFTbin;
      if ( lastBin2read > lastSFTbin )
	lastBin2read = lastSFTbin;

      /* check that requested interval is found in SFT */
      if ( firstBin2read <= lastBin2read ) {

	 /* keep a copy of the data pointer */
	COMPLEX8Sequence*data = thisSFT->data;

	firstBinRead = firstBin2read;
	lastBinRead = lastBin2read;
	offsetBins = firstBin2read - firstSFTbin;
	numBins2read = lastBin2read - firstBin2read + 1;

	/* copy the header */
	*thisSFT = locatalog.data[catPos].header;
	/* restore data pointer */
	thisSFT->data = data;
	/* copy the data */
	memcpy(thisSFT->data->data,
	       locatalog.data[catPos].header.data + offsetBins,
	       numBins2read * sizeof(COMPLEX8));

	/* update the start-frequency entry in the SFT-header to the new value */
	thisSFT->f0 = 1.0 * firstBin2read * thisSFT->deltaF;

      } else {
	/* no data was needed from this SFT (segment) */
	firstBinRead = 0;
	lastBinRead = 0;
      }

    } else {
      /* SFT data had not yet been read - read it */

      /* open and close a file only when necessary, i.e. reading a different file */
      if(strcmp(fname, locator->fname)) {
	if(fp) {
	  fclose(fp);
	  fp = NULL;
	}
	fname = locator->fname;
	fp = fopen(fname,"rb");
	XLALPrintInfo("%s: Opening file '%s'\n", __func__, fname);
	if(!fp) {
	  XLALPrintError("ERROR: Couldn't open file '%s'\n", fname);
	  XLALLOADSFTSERROR(XLAL_EIO);
	}
      }

      /* seek to the position of the SFT in the file (if necessary) */
      if ( locator->offset )
	if ( fseek( fp, locator->offset, SEEK_SET ) == -1 ) {
	  XLALPrintError("ERROR: Couldn't seek to position %ld in file '%s'\n",
			 locator->offset, fname);
	  XLALLOADSFTSERROR(XLAL_EIO);
	}

      /* read SFT data */
      lastBinRead = read_sft_bins_from_fp ( thisSFT, &firstBinRead, firstbin, lastbin, fp );
      XLALPrintInfo ("%s: Read data from %s:%lu: %u - %u\n", __func__, locator->fname, locator->offset, firstBinRead, lastBinRead);
    }
    /* SFT data has been read from file or taken from catalog */

    if(lastBinRead) {
	/* data was actually read */

	if(segments[isft].last == 0) {

	  /* no data was read for this SFT yet: must be first segment */
	  if(firstBinRead != firstbin) {
	    XLALPrintError("ERROR: data gap or overlap at first bin of SFT#%u (GPS %lf)"
			   " expected bin %u, bin %u read from file '%s'\n",
			   isft, GPS2REAL8(thisSFT->epoch),
			   firstbin, firstBinRead, fname);
	    XLALLOADSFTSERROR(XLAL_EIO);
	  }
	  segments[isft].first = firstBinRead;
	  segments[isft].epoch = thisSFT->epoch;

	/* if not first segment, segment must fit at the end of previous data */
	} else if(firstBinRead != segments[isft].last + 1) {
	  XLALPrintError("ERROR: data gap or overlap in SFT#%u (GPS %lf)"
			 " between bin %u read from file '%s' and bin %u read from file '%s'\n",
			 isft, GPS2REAL8(thisSFT->epoch),
			 segments[isft].last, segments[isft].lastfrom->fname,
			 firstBinRead, fname);
	  XLALLOADSFTSERROR(XLAL_EIO);
	}

	/* consistency checks */
	if(deltaF != thisSFT->deltaF) {
	  XLALPrintError("ERROR: deltaF mismatch (%f/%f) in SFT read from file '%s'\n",
			 thisSFT->deltaF, deltaF, fname);
	  XLALLOADSFTSERROR(XLAL_EIO);
	}
	if(!GPSEQUAL(segments[isft].epoch, thisSFT->epoch)) {
	  XLALPrintError("ERROR: GPS epoch mismatch (%f/%f) in SFT read from file '%s'\n",
			 GPS2REAL8(segments[isft].epoch), GPS2REAL8(thisSFT->epoch), fname);
	  XLALLOADSFTSERROR(XLAL_EIO);
	}

	/* data is ok, add to SFT */
	segments[isft].last               = lastBinRead;
	segments[isft].lastfrom           = locator;
        memcpy( sftVector->data[isft].name, locatalog.data[catPos].header.name, sizeof(sftVector->data[isft].name));
	sftVector->data[isft].sampleUnits = locatalog.data[catPos].header.sampleUnits;
	memcpy(sftVector->data[isft].data->data + (firstBinRead - firstbin),
	       thisSFT->data->data,
	       (lastBinRead - firstBinRead + 1) * sizeof(COMPLEX8));

      } else if(!firstBinRead) {
	/* no needed data had been in this segment */
        XLALPrintInfo ( "%s: No data read from %s:%lu\n", __func__, locator->fname, locator->offset);

	/* set epoch if not yet set, if already set, check it */
	if(GPSZERO(segments[isft].epoch))
	  segments[isft].epoch = thisSFT->epoch;
	else if (!GPSEQUAL(segments[isft].epoch, thisSFT->epoch)) {
	  XLALPrintError("ERROR: GPS epoch mismatch (%f/%f) in SFT read from file '%s'\n",
			 GPS2REAL8(segments[isft].epoch), GPS2REAL8(thisSFT->epoch), fname);
	  XLALLOADSFTSERROR(XLAL_EIO);
	}

      } else {
	/* failed to read data */

	XLALPrintError("ERROR: Error (%u) reading SFT from file '%s'\n", firstBinRead, fname);
	XLALLOADSFTSERROR(XLAL_EIO);
      }
  }

  /* close the last file */
  if(fp) {
    fclose(fp);
    fp = NULL;
  }

  /* check that all SFTs are complete */
  for(UINT4 isft = 0; isft < nSFTs; isft++) {
    if(segments[isft].last == lastbin) {
      sftVector->data[isft].f0 = 1.0 * firstbin * deltaF;
      sftVector->data[isft].epoch = segments[isft].epoch;
      sftVector->data[isft].deltaF = deltaF;
    } else {
      if (segments[isft].last)
	XLALPrintError("ERROR: data missing at end of SFT#%u (GPS %lf)"
		       " expected bin %u, bin %u read from file '%s'\n",
		       isft, GPS2REAL8(segments[isft].epoch),
		       lastbin, segments[isft].last,
		       segments[isft].lastfrom->fname);
      else
	XLALPrintError("ERROR: no data could be read for SFT#%u (GPS %lf)\n",
		       isft, GPS2REAL8(segments[isft].epoch));
      XLALLOADSFTSERROR(XLAL_EIO);
    }
  }

  /* cleanup  */
  XLALFree(segments);
  XLALFree(locatalog.data);
  XLALDestroySFT(thisSFT);

  return(sftVector);

} /* XLALLoadSFTs() */


/**
 * Function to load a catalog of SFTs from possibly different detectors.
 * This is similar to XLALLoadSFTs except that the input SFT catalog is
 * allowed to contain multiple ifos. The output is the structure
 * MultiSFTVector which is a vector of (pointers to) SFTVectors, one for
 * each ifo found in the catalog. As in XLALLoadSFTs, fMin and fMax can be
 * set to -1 to get the full SFT from the lowest to the highest frequency
 * bin found in the SFT.
 *
 * output SFTvectors are sorted alphabetically by detector-name
 *
 * NOTE: this is basically a backwards-compatible API wrapper to
 * XLALLoadMultiSFTsFromView(), which takes a MultiSFTCatalogView as input instead.
 *
 */
MultiSFTVector *
XLALLoadMultiSFTs (const SFTCatalog *inputCatalog,   /**< The 'catalogue' of SFTs to load */
		   REAL8 fMin,		             /**< minumum requested frequency (-1 = read from lowest) */
		   REAL8 fMax		             /**< maximum requested frequency (-1 = read up to highest) */
		   )

{
  XLAL_CHECK_NULL ( (inputCatalog != NULL) && (inputCatalog->length != 0), XLAL_EINVAL );

  MultiSFTCatalogView *multiCatalogView;
  // get the (alphabetically-sorted!) multiSFTCatalogView
  XLAL_CHECK_NULL ( (multiCatalogView = XLALGetMultiSFTCatalogView ( inputCatalog )) != NULL, XLAL_EFUNC );

  MultiSFTVector *multiSFTs;
  XLAL_CHECK_NULL ( ( multiSFTs = XLALLoadMultiSFTsFromView ( multiCatalogView, fMin, fMax )) != NULL, XLAL_EFUNC );

  /* free memory and exit */
  XLALDestroyMultiSFTCatalogView ( multiCatalogView );

  return multiSFTs;

} /* XLALLoadMultiSFTs() */


/**
 * This function loads a MultiSFTVector from a given input MultiSFTCatalogView,
 * otherwise the documentation of XLALLoadMultiSFTs() applies.
 *
 * Note: this is basically the core-function of XLALLoadMultiSFTs() doing the
 * actual work.
 *
 * Note2: we keep the IFO sort-order of the input multiCatalogView
 */
MultiSFTVector *
XLALLoadMultiSFTsFromView ( const MultiSFTCatalogView *multiCatalogView,/**< The multi-SFT catalogue view of SFTs to load */
                            REAL8 fMin,		             		/**< minumum requested frequency (-1 = read from lowest) */
                            REAL8 fMax		             		/**< maximum requested frequency (-1 = read up to highest) */
                            )
{
  XLAL_CHECK_NULL ( multiCatalogView != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( multiCatalogView->length != 0, XLAL_EINVAL );

  UINT4 numIFOs = multiCatalogView->length;

  /* create multi sft vector */
  MultiSFTVector *multiSFTs;
  XLAL_CHECK_NULL ( (multiSFTs = XLALCalloc(1, sizeof(*multiSFTs))) != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL ( (multiSFTs->data = XLALCalloc ( numIFOs, sizeof(*multiSFTs->data))) != NULL, XLAL_ENOMEM );
  multiSFTs->length = numIFOs;

  for ( UINT4 X = 0; X < numIFOs; X++ )
    {
      if( ( multiSFTs->data[X] = XLALLoadSFTs ( &(multiCatalogView->data[X]), fMin, fMax ) ) == NULL )
        {
          /* free sft vectors created previously in loop */
          XLALDestroyMultiSFTVector ( multiSFTs );
          XLAL_ERROR_NULL ( XLAL_EFUNC, "Failed to XLALLoadSFTs() for IFO X = %d\n", X );
        } // if XLALLoadSFTs() failed

    } // for X < numIFOs

  // return final multi-SFT vector
  return multiSFTs;

} // XLALLoadMultiSFTsFromView()


/**
 * Write the given SFTtype to a FILE pointer.
 * Add the comment to SFT if SFTcomment != NULL.
 *
 * NOTE: the comment written into the SFT-file contains the 'sft->name' field concatenated with
 * the user-specified 'SFTcomment'
 */
int
XLALWriteSFT2FilePointer(
  const SFTtype *sft,	        /**< SFT to write to disk */
  FILE *fp,		        /**< pointer to open file */
  const CHAR* SFTwindowtype,    /**< window applied to SFT, if any */
  const REAL8 SFTwindowparam,   /**< parameter of window */
  const CHAR *SFTcomment        /**< optional comment */
  )
{
  UINT4 comment_len = 0;
  CHAR *_SFTcomment;
  UINT4 pad_len = 0;
  CHAR pad[] = {0, 0, 0, 0, 0, 0, 0};	/* for comment-padding */
  _SFT_header_t rawheader;

  /* check input consistency */
  if (!sft || !sft->data || sft->deltaF <= 0 || sft->f0 < 0 || sft->data->length ==0 )
    XLAL_ERROR ( XLAL_EINVAL );
  if (!( (sft->epoch.gpsSeconds >= 0) && (sft->epoch.gpsNanoSeconds >= 0) ))
    XLAL_ERROR ( XLAL_EINVAL );
  if (!( sft->epoch.gpsNanoSeconds < 1000000000 ))
    XLAL_ERROR ( XLAL_EINVAL );
  if ( !XLALIsValidCWDetector(sft->name) ) {
    XLALPrintError ("\nInvalid detector prefix '%c%c'\n\n", sft->name[0], sft->name[1] );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  if ( !fp )
    XLAL_ERROR ( XLAL_EINVAL );

  XLAL_CHECK ( SFTwindowtype != NULL, XLAL_EFAULT );

  /* concat sft->name + SFTcomment for SFT-file comment-field */
  comment_len = strlen(sft->name) + 1;
  if ( SFTcomment )
    comment_len += strlen(SFTcomment) + 1;	/* separate by "\n" */

  if ( (_SFTcomment = XLALCalloc( comment_len, sizeof(CHAR) )) == NULL ) {
    XLAL_ERROR( XLAL_ENOMEM );
  }
  strcpy ( _SFTcomment, sft->name );
  if ( SFTcomment ) {
    strcat ( _SFTcomment, "\n" );
    strcat ( _SFTcomment, SFTcomment );
  }

  /* comment length including null terminator to string must be an
   * integer multiple of eight bytes.
   */
  pad_len = (8 - (comment_len % 8)) % 8;

  /* ----- fill out header */
  rawheader.version        		= MAX_SFT_VERSION;
  rawheader.gps_sec        		= sft->epoch.gpsSeconds;
  rawheader.gps_nsec       		= sft->epoch.gpsNanoSeconds;
  rawheader.tbase          		= TSFTfromDFreq ( sft->deltaF );
  rawheader.first_frequency_index 	= lround ( sft->f0 / sft->deltaF );
  rawheader.nsamples       		= sft->data->length;
  rawheader.crc64          		= 0;	/* set to 0 for crc-calculation */
  rawheader.detector[0]    		= sft->name[0];
  rawheader.detector[1]    		= sft->name[1];
  rawheader.comment_length 		= comment_len + pad_len;

  XLAL_CHECK ( build_sft_windowspec( &rawheader.windowspec, NULL, SFTwindowtype, SFTwindowparam ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* ----- compute CRC */
  rawheader.crc64 = crc64((const unsigned char*)&rawheader, sizeof(rawheader), ~(0ULL));

  rawheader.crc64 = crc64((const unsigned char*)_SFTcomment, comment_len, rawheader.crc64);
  rawheader.crc64 = crc64((const unsigned char*)pad, pad_len, rawheader.crc64);

  rawheader.crc64 = crc64((const unsigned char*) sft->data->data, sft->data->length * sizeof( *sft->data->data ), rawheader.crc64);

  /* ----- write the header to file */
  if (1 != fwrite( &rawheader, sizeof(rawheader), 1, fp) ) {
    XLAL_ERROR ( XLAL_EIO );
  }

  /* ----- write the comment to file */
  if ( comment_len != fwrite( _SFTcomment, 1, comment_len, fp) ) {
    XLAL_ERROR ( XLAL_EIO );
  }
  if (pad_len != fwrite( pad, 1, pad_len, fp) ) {
    XLAL_ERROR ( XLAL_EIO );
  }

  XLALFree ( _SFTcomment );

  /* write the data to the file.  Data must be packed REAL,IMAG,REAL,IMAG,... */
  if ( sft->data->length != fwrite( sft->data->data, sizeof(*sft->data->data), sft->data->length, fp) ) {
    XLAL_ERROR ( XLAL_EIO );
  }

  return XLAL_SUCCESS;

} /* XLALWriteSFT2FilePointer() */


/**
 * Write the given SFTtype to a SFT file with the supplied filename.
 * Add the comment to SFT if SFTcomment != NULL.
 *
 * NOTE: the comment written into the SFT-file contains the 'sft->name' field concatenated with
 * the user-specified 'SFTcomment'
 */
int
XLALWriteSFT2NamedFile(
  const SFTtype *sft,           /**< SFT to write to disk */
  const CHAR *SFTfilename,      /**< SFT filename */
  const CHAR* SFTwindowtype,    /**< window applied to SFT, if any */
  const REAL8 SFTwindowparam,   /**< parameter of window */
  const CHAR *SFTcomment        /**< optional comment */
  )
{
  FILE  *fp = NULL;

  /*   Make sure the arguments are not NULL */
  if (!( sft ))
    XLAL_ERROR ( XLAL_EINVAL );
  if (!( sft->data ))
    XLAL_ERROR ( XLAL_EINVAL );
  if (!( SFTfilename ))
    XLAL_ERROR ( XLAL_EINVAL );
 
  if ( !XLALIsValidCWDetector(sft->name) ) {
    XLALPrintError ("\nInvalid detector prefix '%c%c'\n\n", sft->name[0], sft->name[1] );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  XLAL_CHECK ( SFTwindowtype != NULL, XLAL_EFAULT );

  /* open SFT-file for writing */
  if ( (fp = fopen ( SFTfilename, "wb" )) == NULL )
    {
      XLALPrintError ("\nFailed to open file '%s' for writing: %s\n\n", SFTfilename, strerror(errno));
      XLAL_ERROR ( XLAL_EIO );
    }

  /* write SFT to file */
  if ( XLALWriteSFT2FilePointer (sft, fp, SFTwindowtype, SFTwindowparam, SFTcomment) != XLAL_SUCCESS ) {
    XLAL_ERROR ( XLAL_EIO );
  }

  fclose(fp);

  return XLAL_SUCCESS;

} /* XLALWriteSFT2NamedFile() */


/**
 * Write the given SFTtype to a SFT file with a standard (\cite SFT-spec) filename.
 * Add the comment to SFT if SFTcomment != NULL.
 *
 * NOTE: the comment written into the SFT-file contains the 'sft->name' field concatenated with
 * the user-specified 'SFTcomment'
 *
 * NOTE: The SFT filename spec is updated to reflect the filename of the written SFT, by
 * setting the \c numSFTs, \c detector, \c SFTtimebase, \c gpsStart, and \c SFTspan fields.
 * if needed, the SFT filename can be reconstructed with XLALBuildSFTFilenameFromSpec()
 */
int
XLALWriteSFT2StandardFile(
  const SFTtype *sft,           /**< SFT to write to disk */
  SFTFilenameSpec *SFTfnspec,   /**< SFT filename specification used to construct filename */
  const CHAR *SFTcomment        /**< optional comment */
  )
{

  /*   Make sure the arguments are not NULL */
  if (!( sft ))
    XLAL_ERROR ( XLAL_EINVAL );
  if (!( sft->data ))
    XLAL_ERROR ( XLAL_EINVAL );
 
  if ( !XLALIsValidCWDetector(sft->name) ) {
    XLALPrintError ("\nInvalid detector prefix '%c%c'\n\n", sft->name[0], sft->name[1] );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  XLAL_CHECK ( SFTfnspec != NULL, XLAL_EFAULT );

  const UINT4 Tsft = (UINT4) round ( 1.0 / sft->deltaF );

  SFTfnspec->numSFTs         = 1;
  SFTfnspec->detector[0]     = sft->name[0];
  SFTfnspec->detector[1]     = sft->name[1];
  SFTfnspec->detector[2]     = 0;
  SFTfnspec->SFTtimebase     = Tsft;
  SFTfnspec->gpsStart        = sft->epoch.gpsSeconds;

  /* calculate sft 'duration' -- may be different from timebase if nanosecond of sft-epoch is non-zero */
  SFTfnspec->SFTspan = Tsft;
  if ( sft->epoch.gpsNanoSeconds > 0) {
    SFTfnspec->SFTspan += 1;
  }

  // build SFT filename
  char *SFTfilename = XLALBuildSFTFilenameFromSpec ( SFTfnspec );
  XLAL_CHECK ( SFTfilename != NULL, XLAL_EFUNC );

  // write SFT
  XLAL_CHECK ( XLALWriteSFT2NamedFile ( sft, SFTfilename, SFTfnspec->window_type, SFTfnspec->window_param, SFTcomment ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLALFree( SFTfilename );
  
  return XLAL_SUCCESS;

} /* XLALWriteSFT2StandardFile() */


/**
 * Write the given SFTVector to a single merged SFT file with the supplied filename.
 * Add the comment to SFT if SFTcomment != NULL.
 */
int
XLALWriteSFTVector2NamedFile(
  const SFTVector *sftVect,     /**< SFT vector to write to disk */
  const CHAR *SFTfilename,      /**< SFT filename */
  const CHAR* SFTwindowtype,    /**< window applied to SFT, if any */
  const REAL8 SFTwindowparam,   /**< parameter of window */
  const CHAR *SFTcomment        /**< optional comment */
  )
{
  XLAL_CHECK ( sftVect != NULL, XLAL_EINVAL );
  XLAL_CHECK ( sftVect->data != NULL, XLAL_EINVAL );
  XLAL_CHECK ( sftVect->length > 0, XLAL_EINVAL );
  XLAL_CHECK ( SFTfilename != NULL, XLAL_EINVAL );

  const UINT4 numSFTs = sftVect->length;

  /* open SFT-file for writing */
  FILE *fp;
  XLAL_CHECK ( (fp = fopen ( SFTfilename, "wb" )) != NULL, XLAL_EIO, "Failed to open '%s' for writing: %s\n\n", SFTfilename, strerror(errno));

  for ( UINT4 k = 0; k < numSFTs; k++ )
    {
      SFTtype *sft = &( sftVect->data[k] );

      /* write the k^th sft */
      XLAL_CHECK ( XLALWriteSFT2FilePointer ( sft, fp, SFTwindowtype, SFTwindowparam, SFTcomment ) == XLAL_SUCCESS, XLAL_EFUNC );

    } // for k < numSFTs

  fclose(fp);

  return XLAL_SUCCESS;

} /* XLALWriteSFTVector2NamedFile() */


/**
 * Write the given SFTVector to SFT file(s) with a standard (\cite SFT-spec) filename(s).
 * Add the comment to SFT if SFTcomment != NULL.
 *
 * NOTE: The SFT filename spec is updated to reflect the filename of the written SFT, by
 * setting the \c numSFTs, \c detector, \c SFTtimebase, \c gpsStart, and \c SFTspan fields.
 * if needed, the SFT filename can be reconstructed with XLALBuildSFTFilenameFromSpec()
 */
int
XLALWriteSFTVector2StandardFile(
  const SFTVector *sftVect,     /**< SFT vector to write to disk */
  SFTFilenameSpec *SFTfnspec,   /**< SFT filename specification used to construct filename(s) */
  const CHAR *SFTcomment,       /**< optional comment */
  const BOOLEAN merged          /**< If true, write a single merged SFT file; otherwise, write individual SFT files */
  )
{
  XLAL_CHECK ( sftVect != NULL, XLAL_EINVAL );
  XLAL_CHECK ( sftVect->data != NULL, XLAL_EINVAL );
  XLAL_CHECK ( sftVect->length > 0, XLAL_EINVAL );

  XLAL_CHECK ( SFTfnspec != NULL, XLAL_EFAULT );

  const UINT4 numSFTs           = sftVect->length;
  const SFTtype *sftStart       = &(sftVect->data[0]);
  const SFTtype *sftEnd         = &(sftVect->data[numSFTs-1]);
  const LIGOTimeGPS *epochStart = &(sftStart->epoch);
  const LIGOTimeGPS *epochEnd   = &(sftEnd->epoch);

  const UINT4 Tsft = (UINT4) round ( 1.0 / sftStart->deltaF );

  if (merged) {   // write a single merged SFT file

    SFTfnspec->numSFTs         = numSFTs;
    SFTfnspec->detector[0]     = sftStart->name[0];
    SFTfnspec->detector[1]     = sftStart->name[1];
    SFTfnspec->detector[2]     = 0;
    SFTfnspec->SFTtimebase     = Tsft;
    SFTfnspec->gpsStart        = sftStart->epoch.gpsSeconds;

    /* calculate time interval covered -- may be different from timebase if nanosecond of sft-epochs are non-zero */
    SFTfnspec->SFTspan = epochEnd->gpsSeconds - epochStart->gpsSeconds + Tsft;
    if ( epochStart->gpsNanoSeconds > 0) {
      SFTfnspec->SFTspan += 1;
    }
    if ( epochEnd->gpsNanoSeconds > 0) {
      SFTfnspec->SFTspan += 1;
    }

    // build SFT filename
    char *SFTfilename = XLALBuildSFTFilenameFromSpec ( SFTfnspec );
    XLAL_CHECK ( SFTfilename != NULL, XLAL_EFUNC );

    // write SFT
    XLAL_CHECK ( XLALWriteSFTVector2NamedFile ( sftVect, SFTfilename, SFTfnspec->window_type, SFTfnspec->window_param, SFTcomment ) == XLAL_SUCCESS, XLAL_EFUNC );

    XLALFree ( SFTfilename );

  } else {   // write individual SFT files

    // local copy of 'SFTfnspec' so that we can return name of first written SFT
    SFTFilenameSpec spec = *SFTfnspec;

    for ( UINT4 k = 0; k < numSFTs; k++ ) {

      SFTtype *sft = &(sftVect->data[k]);

      // return spec of first written SFT in 'SFTfnspec'; otherwise use local copy
      SFTFilenameSpec *p_spec = ( k == 0 ) ? SFTfnspec : &spec;

      XLAL_CHECK ( XLALWriteSFT2StandardFile( sft, p_spec, SFTcomment ) == XLAL_SUCCESS, XLAL_EFUNC );

    } // for k < numSFTs

  }

  return XLAL_SUCCESS;

} /* XLALWriteSFTVector2StandardFile() */


/**
 * Verify that the contents of a SFT file are valid.
 *
 * This is just an XLAL wrapper to the SFTReferenceLibrary function ValidateSFTFile().
 *
 * \return: XLAL_SUCCESS if no validation errors encountered.
 */
int
XLALCheckSFTFileIsValid ( const char *fname )
{
    int errcode = ValidateSFTFile(fname);
    XLAL_CHECK ( errcode==0, XLAL_EFUNC, "SFT validation error on file '%s': code %d (%s)", fname, errcode, SFTErrorMessage( errcode ) );
    return XLAL_SUCCESS;
} /* XLALCheckSFTFileIsValid */


/* a little endian-swapper needed for SFT reading/writing */
void
endian_swap(CHAR * pdata, size_t dsize, size_t nelements)
{
  UINT4 i, j, indx;
  CHAR tempbyte;

  if (dsize <= 1) return;

  for (i=0; i<nelements; i++)
    {
      indx = dsize;
      for (j=0; j<dsize/2; j++)
	{
	  tempbyte = pdata[j];
	  indx = indx - 1;
	  pdata[j] = pdata[indx];
	  pdata[indx] = tempbyte;
	}

      pdata = pdata + dsize;
    }

  return;

} /* endian swap */


/**
 * Open an "SFT" defined by the SFT-locator, return a FILE-pointer to the beginning of this SFT.
 * \note The returned filepointer could point to an SFT-block within a merged SFT-file,
 * so you should not assume that SEEK_SET takes you to the beginning of this block!
 * (instead you'd have to save the current position returned by this function, which \
 * points to the beginning of the block!)
 *
 * NOTE: Ideally this should be the *ONLY* function using the internal structure of the opaque
 * SFTLocator type
 *
 */
FILE *
fopen_SFTLocator ( const struct tagSFTLocator *locator )
{
  FILE *fp = NULL;
  CHAR *fname;

  if ( !locator )
    return NULL;

  fname = locator->fname;
  if ( (fp = fopen( fname, "rb" )) == NULL )
    {
      XLALPrintError ("\nFailed to open SFT '%s' for reading: %s\n\n", fname, strerror(errno) );
      return NULL;
    }

  if ( fseek( fp, locator->offset, SEEK_SET ) == -1 )
    {
      XLALPrintError ("\nFailed to set fp-offset to '%ld': %s\n\n", locator->offset, strerror(errno) );
      fclose(fp);
      return NULL;
    }

  return fp;

} /* fopen_SFTLocator() */


/**
 * Read valid SFT version-number at position fp, and determine if we need to
 * endian-swap the data.
 * Restores filepointer to original position before returning.
 *
 * RETURN: 0 = OK, -1 = ERROR
 */
int
read_SFTversion_from_fp ( UINT4 *version, BOOLEAN *need_swap, FILE *fp )
{
  long save_filepos;
  REAL8 ver;

  /* store fileposition for restoring in case of failure */
  if ( ( save_filepos = ftell(fp) ) == -1 )
    {
      XLALPrintError ("\nftell() failed: %s\n\n", strerror(errno) );
      return -1;
    }

  /* read version-number */
  if  ( 1 != fread ( &ver, sizeof( ver ), 1, fp) )
    {
      if (lalDebugLevel) XLALPrintError ("\nCould not read version-number from file\n\n");
      goto failed;
    }


  /* figure out endian-ness and check version-range */
  for ( *version = MAX_SFT_VERSION; *version >= MIN_SFT_VERSION; --(*version) )
    {
      REAL8 vertest = *version;
      if ( ! memcmp( &ver, &vertest, sizeof( ver ) ) ) {
	*need_swap = FALSE;
	break;
      }
      endian_swap( (char*)(&vertest), sizeof( vertest ), 1 );
      if ( ! memcmp( &ver, &vertest, sizeof( ver ) ) ) {
	*need_swap = TRUE;
	break;
      }
    }
  if ( *version < MIN_SFT_VERSION ) {
    if ( lalDebugLevel ) {
      unsigned char *v = (unsigned char*)(&ver);
      XLALPrintError( "\nERROR: illegal SFT-version (%X %X %X %X %X %X %X %X) not within [%.0f, %.0f]\n",
		     v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],
		     (float)MIN_SFT_VERSION, (float)MAX_SFT_VERSION );
    }
    goto failed;
  }

  /* restore initial filepointer position */
  if ( fseek ( fp, save_filepos, SEEK_SET ) == -1 )
    {
      XLALPrintError ("\nfseek() failed to return to intial fileposition: %s\n\n", strerror(errno) );
      goto failed;
    }

  return 0;

 failed:
  fseek ( fp, save_filepos, SEEK_SET );
  return -1;

} /* read_SFTversion_from_fp() */


/* Try to read an SFT-header (of ANY VALID SFT-VERSION) at the given FILE-pointer fp,
 * and return the SFT-header, SFT-version-number and number of frequency-samples in the SFT.
 *
 * Sets the filepointer fp at the end of the header if successful, leaves it at
 * initial position if not.
 *
 * RETURN 0 = OK, -1 on ERROR
 *
 * We do basic checking of compliance with the SFT-spec (<tt>LIGO-T04164-01-Z</tt>)
 * as far as a single header is concerned.
 *
 * NOTE: fatal errors will produce a XLALPrintError() error-message, but
 * non-fatal 'SFT-format'-errors will only output error-messages if lalDebugLevel > 0.
 * --> this function can therefore be used to check if a given file actually contains SFTs
 *
 *
 */
int
read_sft_header_from_fp (FILE *fp, SFTtype *header, UINT4 *version, UINT8 *crc64, UINT2 *SFTwindowspec, BOOLEAN *swapEndian, CHAR **SFTcomment, UINT4 *numBins )
{
  SFTtype XLAL_INIT_DECL(head);
  UINT4 nsamples;
  CHAR *comm = NULL;
  UINT8 ref_crc = 0;
  UINT8 header_crc;
  UINT2 windowspec;

  UINT4 ver;
  BOOLEAN need_swap;
  long save_filepos;

  if ( !header || !version || !numBins || !fp  )
    {
      XLALPrintError ("\nERROR read_sft_header_from_fp(): called with NULL input\n\n");
      return -1;
    }
  if ( SFTcomment && ((*SFTcomment) != NULL) )
    {
      XLALPrintError ("\nERROR: Comment-string passed to read_sft_header_from_fp() is not empty!\n\n");
      return -1;
    }

  /* store fileposition for restoring in case of failure */
  if ( ( save_filepos = ftell(fp) ) == -1 )
    {
      XLALPrintError ("\nftell() failed: %s\n\n", strerror(errno) );
      return -1;
    }

  if ( read_SFTversion_from_fp ( &ver, &need_swap, fp ) != 0 )
    return -1;

  /* read this SFT-header */
  XLAL_INIT_MEM(head);

  if ( MIN_SFT_VERSION <= ver && ver <= MAX_SFT_VERSION )
    {
      if ( read_header_from_fp ( fp, &head, &nsamples, &header_crc, &ref_crc, &windowspec, &comm, need_swap ) != 0 )
        goto failed;
    }
  else
    {
      XLALPrintError ("\nUnsupported SFT-version %d.\n\n", ver);
      goto failed;
    } /* switch(ver) */


  /* ----- some general SFT-header consistency-checks */
  if ( (head.epoch.gpsSeconds < 0 ) || ( head.epoch.gpsNanoSeconds < 0 ) || ( head.epoch.gpsNanoSeconds >= 1000000000)  )
    {
      XLALPrintError ("\nInvalid GPS-epoch in SFT : [%d, %d]!\n\n",
					  head.epoch.gpsSeconds, head.epoch.gpsNanoSeconds );
      goto failed;
    }

  if ( head.deltaF <= 0 )
    {
      XLALPrintError ("\nNegative frequency-spacing in SFT!\n\n");
      goto failed;
    }

  if ( head.f0 < 0 )
    {
      XLALPrintError ("\nNegative start-frequency in SFT!\n\n");
      goto failed;
    }

  /* ok */
  (*header) = head;
  (*version) = ver;

  if ( SFTcomment )	  /* return of comment is optional */
    (*SFTcomment) = comm;
  else
    if ( comm ) LALFree(comm);

  (*swapEndian) = need_swap;
  (*crc64) = ref_crc;
  (*numBins) = nsamples;

  if ( SFTwindowspec ) {     /* return of windowspec is optional */
    (*SFTwindowspec) = windowspec;
  }

  return 0;

  /* ---------- */
 failed:
  /* restore filepointer initial position  */
  if ( fseek ( fp, save_filepos, SEEK_SET ) == -1 )
    XLALPrintError ("\nfseek() failed to return to intial fileposition: %s\n\n", strerror(errno) );

  /* free comment  if we allocated one */
  if ( comm )
    LALFree ( comm );

  return -1;

} /* read_sft_header_from_fp() */


/* ----- SFT header-reading function:
 *
 * return general SFTtype header, place filepointer at the end of the header if it succeeds,
 * set fp to initial position if it fails.
 * RETURN: 0 = OK, -1 = ERROR
 *
 * NOTE: fatal errors will produce a XLALPrintError() error-message, but
 * non-fatal 'SFT-format'-errors will only output error-messages if lalDebugLevel > 0.
 *
 */
int
read_header_from_fp ( FILE *fp, SFTtype *header, UINT4 *nsamples, UINT8 *header_crc64, UINT8 *ref_crc64, UINT2 *SFTwindowspec, CHAR **SFTcomment, BOOLEAN swapEndian)
{
  _SFT_header_t rawheader;
  long save_filepos;
  CHAR *comm = NULL;
  UINT8 crc;


  /* check input-consistency */
  if ( !fp || !header || !nsamples || !SFTcomment )
    {
      XLALPrintError ( "\nERROR read_header_from_fp(): called with NULL input!\n\n");
      return -1;
    }
  if ( SFTcomment && (*SFTcomment != NULL) )
    {
      XLALPrintError ("\nERROR: Comment-string passed to read_header_from_fp() is not NULL!\n\n");
      return -1;
    }

  /* store fileposition for restoring in case of failure */
  if ( ( save_filepos = ftell(fp) ) == -1 )
    {
      XLALPrintError ("\nERROR: ftell() failed: %s\n\n", strerror(errno) );
      return -1;
    }

  /* read the whole header */
  if (fread( &rawheader, sizeof(rawheader), 1, fp) != 1)
    {
      if (lalDebugLevel) XLALPrintError ("\nCould not read header. %s\n\n", strerror(errno) );
      goto failed;
    }

  /* ----- compute CRC for the header:
   * NOTE: the CRC checksum is computed on the *bytes*, not the numbers,
   * so this must be computed before any endian-swapping.
   */
  {
    UINT8 save_crc = rawheader.crc64;
    rawheader.crc64 = 0;

    crc = crc64((const unsigned char*)&rawheader, sizeof(rawheader), ~(0ULL));

    rawheader.crc64 = save_crc;
    /* NOTE: we're not done with crc yet, because we also need to
     * include the comment's CRC , see below
     */
  }/* compute crc64 checksum */

  /* ----- swap endian-ness if required ----- */
  if (swapEndian)
    {
      endian_swap((CHAR*)(&rawheader.version), 			sizeof(rawheader.version) 		, 1);
      endian_swap((CHAR*)(&rawheader.gps_sec), 			sizeof(rawheader.gps_sec) 		, 1);
      endian_swap((CHAR*)(&rawheader.gps_nsec), 		sizeof(rawheader.gps_nsec) 		, 1);
      endian_swap((CHAR*)(&rawheader.tbase), 			sizeof(rawheader.tbase) 		, 1);
      endian_swap((CHAR*)(&rawheader.first_frequency_index), 	sizeof(rawheader.first_frequency_index) , 1);
      endian_swap((CHAR*)(&rawheader.nsamples), 		sizeof(rawheader.nsamples) 		, 1);
      endian_swap((CHAR*)(&rawheader.crc64),			sizeof(rawheader.crc64)			, 1);
      endian_swap((CHAR*)(&rawheader.windowspec),		sizeof(rawheader.windowspec)		, 1);
      endian_swap((CHAR*)(&rawheader.comment_length),		sizeof(rawheader.comment_length)	, 1);
      /* ----- */

    } /* if endian_swap */

  /* double-check version-number */
  if ( !( MIN_SFT_VERSION <= rawheader.version && rawheader.version <= MAX_SFT_VERSION ) )
    {
      XLALPrintError ("\nWrong SFT-version %g in read_header_from_fp()\n\n", rawheader.version );
      goto failed;
    }

  if ( rawheader.nsamples <= 0 )
    {
      XLALPrintError ("\nNon-positive number of samples in SFT!\n\n");
      goto failed;
    }

  if ( rawheader.comment_length < 0 )
    {
      XLALPrintError ("\nNegative comment-length in SFT!\n\n");
      goto failed;
    }

  if ( rawheader.comment_length % 8 != 0 )
    {
      XLALPrintError ("\nComment-length must be multiple of 8 bytes!\n\n");
      goto failed;
    }

  const CHAR detector[3] = { rawheader.detector[0], rawheader.detector[1], 0 };
  if ( ! XLALIsValidCWDetector ( detector ) )
    {
      XLALPrintError ("\nIllegal detector-name in SFT: '%s'\n\n", detector );
      goto failed;
    }

  if ( rawheader.version == 2 )
    {
      rawheader.windowspec = 0;   /* window of version 2 SFTs is unknown */
    }

  /* ----- Now read comment (if any) ----- */
  comm = NULL;
  if ( rawheader.comment_length )
    {
      CHAR *ptr;
      if ( (comm = LALCalloc(1, rawheader.comment_length )) == NULL )
	{
	  XLALPrintError ("\nFATAL: out of memory ...\n\n");
	  goto failed;
	}
      if ( (size_t)rawheader.comment_length != fread ( comm, 1, rawheader.comment_length, fp ) )
	{
	  XLALPrintError ("\nCould not read %d-bytes comment\n\n", rawheader.comment_length);
	  goto failed;
	}

      /* check that comment is 0-terminated */
      if ( comm[ rawheader.comment_length - 1] != 0 )
	{
	  XLALPrintError ("\nComment is not properly 0-terminated!\n\n");
	  goto failed;
	}

      /* check that no NON-NULL bytes after first NULL in comment (->spec) */
      ptr = strchr ( comm, 0 );	/* guaranteed to find sth, after previous check */
      while ( ptr < (comm + rawheader.comment_length - 1) )
	if ( *ptr++ != 0 )
	  {
	    XLALPrintError ("\nNon-NULL bytes found after comment-end!\n\n");
	    goto failed;
	  }

      /* comment length including null terminator to string must be an
       * integer multiple of eight bytes. comment==NULL means 'no
       * comment'
       */
      if (SFTcomment)
	{
	  CHAR pad[] = {0, 0, 0, 0, 0, 0, 0};	/* for comment-padding */
	  UINT4 comment_len = strlen(comm) + 1;
	  UINT4 pad_len = (8 - (comment_len % 8)) % 8;

	  crc = crc64((const unsigned char*)comm, comment_len, crc);
	  crc = crc64((const unsigned char*)pad, pad_len, crc);
	}

    } /* if comment_length > 0 */

  /*  ok: */
  memset ( header, 0, sizeof( *header ) );

  strncpy( header->name, detector, sizeof( header->name ) );

  header->epoch.gpsSeconds 	= rawheader.gps_sec;
  header->epoch.gpsNanoSeconds 	= rawheader.gps_nsec;

  header->f0 			= rawheader.first_frequency_index / rawheader.tbase;
  header->deltaF 		= 1.0 / rawheader.tbase;

  (*nsamples) = rawheader.nsamples;
  (*ref_crc64) = rawheader.crc64;
  (*SFTcomment) = comm;
  (*header_crc64) = crc;
  (*SFTwindowspec) = rawheader.windowspec;

  return 0;

 failed:
  /* restore filepointer initial position  */
  if ( fseek ( fp, save_filepos, SEEK_SET ) == -1 )
    XLALPrintError ("\nfseek() failed to return to intial fileposition: %s\n\n", strerror(errno) );

  /* free comment  if we allocated one */
  if ( comm )
    LALFree (comm);

  return -1;

} /* read_header_from_fp() */


/*
   This function reads an SFT (segment) from an open file pointer into a buffer.
   firstBin2read specifies the first bin to read from the SFT, lastBin2read is the last bin.
   If the SFT contains fewer bins than specified, all bins from the SFT are read.
   The function returns the last bin actually read, firstBinRead
   is set to the first bin actually read. In case of an error, 0 is returned
   and firstBinRead is set to a code further decribing the error condition.
*/
UINT4
read_sft_bins_from_fp ( SFTtype *ret, UINT4 *firstBinRead, UINT4 firstBin2read, UINT4 lastBin2read , FILE *fp )
{
  UINT4 version;
  UINT8 crc64;
  BOOLEAN swapEndian;
  UINT4 numBins2read;
  UINT4 firstSFTbin, lastSFTbin, numSFTbins;
  INT4 offsetBins;
  long offsetBytes;
  volatile REAL8 tmp;	/* intermediate results: try to force IEEE-arithmetic */

  if (!firstBinRead)
    {
      XLALPrintError ( "read_sft_bins_from_fp(): got passed NULL *firstBinRead\n" );
      return((UINT4)-1);
    }

  *firstBinRead = 0;

  if ((ret == NULL) ||
      (ret->data == NULL) ||
      (ret->data->data == NULL))
    {
      XLALPrintError ( "read_sft_bins_from_fp(): got passed NULL SFT*\n" );
      *firstBinRead = 1;
      return(0);
    }

  if (!fp)
    {
      XLALPrintError ( "read_sft_bins_from_fp(): got passed NULL FILE*\n" );
      *firstBinRead = 1;
      return(0);
    }

  if ( firstBin2read > lastBin2read )
    {
      XLALPrintError ("read_sft_bins_from_fp(): Empty frequency-interval requested [%d, %d] bins\n",
		      firstBin2read, lastBin2read );
      *firstBinRead = 1;
      return(0);
    }

  {
    COMPLEX8Sequence*data = ret->data;
    if ( read_sft_header_from_fp (fp, ret, &version, &crc64, NULL, &swapEndian, NULL, &numSFTbins ) != 0 )
      {
	XLALPrintError ("read_sft_bins_from_fp(): Failed to read SFT-header!\n");
	*firstBinRead = 2;
	return(0);
      }
    ret->data = data;
  }

  tmp = ret->f0 / ret->deltaF;
  firstSFTbin = lround ( tmp );
  lastSFTbin = firstSFTbin + numSFTbins - 1;

  /* limit the interval to be read to what's actually in the SFT */
  if ( firstBin2read < firstSFTbin )
    firstBin2read = firstSFTbin;
  if ( lastBin2read > lastSFTbin )
    lastBin2read = lastSFTbin;

  /* check that requested interval is found in SFT */
  /* return 0 (no bins read) if this isn't the case */
  if ( firstBin2read > lastBin2read ) {
    *firstBinRead = 0;
    return(0);
  }

  *firstBinRead = firstBin2read;

  offsetBins = firstBin2read - firstSFTbin;
  offsetBytes = offsetBins * 2 * sizeof( REAL4 );
  numBins2read = lastBin2read - firstBin2read + 1;

  if ( ret->data->length < numBins2read )
    {
      XLALPrintError ("read_sft_bins_from_fp(): passed SFT has not enough bins (%u/%u)\n",
		      ret->data->length, numBins2read );
      *firstBinRead = 1;
      return(0);
    }

  /* seek to the desired bins */
  if ( fseek ( fp, offsetBytes, SEEK_CUR ) != 0 )
    {
      XLALPrintError ( "read_sft_bins_from_fp(): Failed to fseek() to first frequency-bin %d: %s\n",
		       firstBin2read, strerror(errno) );
      *firstBinRead = 3;
      return(0);
    }

  /* actually read the data */
  if ( numBins2read != fread ( ret->data->data, sizeof( COMPLEX8 ), numBins2read, fp ) )
    {
      XLALPrintError ("read_sft_bins_from_fp(): Failed to read %d bins from SFT!\n", numBins2read );
      *firstBinRead = 4;
      return(0);
    }

  /* update the start-frequency entry in the SFT-header to the new value */
  ret->f0 = 1.0 * firstBin2read * ret->deltaF;

  /* take care of normalization and endian-swapping */
  if ( swapEndian )
    {
      UINT4 i;

      for ( i=0; i < numBins2read; i ++ )
	{
	  REAL4 re = crealf(ret->data->data[i]);
	  REAL4 im = cimagf(ret->data->data[i]);

	  if ( swapEndian )
	    {
	      endian_swap( (CHAR *) &re, sizeof ( re ), 1 );
	      endian_swap( (CHAR *) &im, sizeof ( im ), 1 );
	    }

          ret->data->data[i] = crectf( re, im );
	} /* for i < numBins2read */
    } /* swapEndian */

  /* return last bin read */
  return(lastBin2read);

} /* read_sft_bins_from_fp() */


/**
 * Check the SFT-block starting at fp for valid crc64 checksum.
 * Restores filepointer before leaving.
 */
BOOLEAN
has_valid_crc64 ( FILE *fp )
{
  long save_filepos;
  UINT8 computed_crc, ref_crc;
  SFTtype header;
  UINT4 numBins;
  UINT2 windowspec;
  CHAR *SFTcomment = NULL;
  UINT4 data_len;
  char block[BLOCKSIZE];
  UINT4 version;
  BOOLEAN need_swap;

  /* input consistency */
  if ( !fp )
    {
      XLALPrintError ("\nhas_valid_crc64() was called with NULL filepointer!\n\n");
      return FALSE;
    }

  /* store fileposition for restoring in case of failure */
  if ( ( save_filepos = ftell(fp) ) == -1 )
    {
      XLALPrintError ("\nERROR: ftell() failed: %s\n\n", strerror(errno) );
      return -1;
    }

  if ( read_SFTversion_from_fp ( &version, &need_swap, fp ) != 0 )
    return -1;

  if ( !( MIN_SFT_VERSION <= version && version <= MAX_SFT_VERSION ) )
    {
      XLALPrintError ("\nhas_valid_crc64() was called on an invalid version(=%u) SFT.\n\n", version);
      return -1;
    }

  /* ----- compute CRC ----- */
  /* read the header, unswapped, only to obtain it's crc64 checksum */
  if ( read_header_from_fp ( fp, &header, &numBins, &computed_crc, &ref_crc, &windowspec, &SFTcomment, need_swap ) != 0 )
    {
      if ( SFTcomment ) LALFree ( SFTcomment );
      return FALSE;
    }
  if ( SFTcomment ) LALFree ( SFTcomment );

  /* read data in blocks of BLOCKSIZE, computing CRC */
  data_len = numBins * 8 ;	/* COMPLEX8 data */
  while ( data_len > 0 )
    {
      /* read either BLOCKSIZE or amount remaining */
      int toread = (BLOCKSIZE < data_len) ? BLOCKSIZE : data_len;
      if (toread != (int)fread ( block, 1, toread, fp) )
	{
	  XLALPrintError ("\nFailed to read all frequency-bins from SFT.\n\n");
	  return FALSE;
	}
      data_len -= toread;

      /* compute CRC64: don't endian-swap for that! */
      computed_crc = crc64((const unsigned char*)block, toread, computed_crc );

    } /* while data */

  /* check that checksum is consistent */
  return ( computed_crc == ref_crc );

} /* has_valid_crc64 */

/// @}
