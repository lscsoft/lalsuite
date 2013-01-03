/*
 * Copyright (C) 2010 Karl Wette
 * Copyright (C) 2004, 2005 R. Prix, B. Machenschalk, A.M. Sintes
 *
 * crc64() taken from SFTReferenceLibrary.c Copyright (C) 2004 Bruce Allen
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

/**
 * \addtogroup SFTfileIO_h
 * \author R. Prix, B. Machenschalk, A.M. Sintes, B. Krishnan
 *
 * \brief IO library for reading/writing "Short Fourier transform" (SFT) data files.
 *
 * This implements the SFTv2 standard defined in LIGO-T040164-01-Z
 * A previous non-LAL implementation of this standard is found in the "SFT reference library"
 * gravity.phys.uwm.edu:2402/usr/local/cvs/lscsoft sftlib, Copyright (C) 2004 Bruce Allen
 *
 * The function calc_crc64() here is based on crc64() in SFTReferenceLibrary.c.
 *
 *
 */

/*---------- INCLUDES ----------*/
#include <sys/types.h>
#include <errno.h>
#include <string.h>
#include <strings.h>

#ifndef _MSC_VER
#include <dirent.h>
#else
#include <io.h>
#endif

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/SFTfileIO.h>
#include <lal/StringVector.h>
#include <lal/Sequence.h>
#include <lal/ConfigFile.h>
#include <lal/LogPrintf.h>

/*---------- DEFINES ----------*/

#define MIN_SFT_VERSION 1
#define MAX_SFT_VERSION 2

#define TRUE    1
#define FALSE   0

/** blocksize used in SFT-reading for the CRC-checksum computation (has to be multiple of 8 !!) */
#define BLOCKSIZE 8192 * 8

/** size of blocks allocated for SFT data. For Einstein\@home SFTs this should be set to 8000 (externally) */
#ifndef SFTFILEIO_REALLOC_BLOCKSIZE
#define SFTFILEIO_REALLOC_BLOCKSIZE 100
#endif

/*----- Macros ----- */

#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

#define GPSEQUAL(gps1,gps2) (((gps1).gpsSeconds == (gps2).gpsSeconds) && ((gps1).gpsNanoSeconds == (gps2).gpsNanoSeconds))

#define GPSZERO(gps) (((gps).gpsSeconds == 0) && ((gps).gpsNanoSeconds == 0))

/* rounding for positive numbers! */
#define MYROUND(x) ( floor( (x) + 0.5 ) )

/*---------- internal types ----------*/

/* NOTE: the locator is implemented as an OPAQUE type in order to enforce encapsulation
 * of the actual physical storage of SFTs and to ease future extensions of the interface.
 * DO NOT TRY TO USE THIS TYPE OUTSIDE OF THIS MODULE!!
 */
struct tagSFTLocator
{
  CHAR *fname;		/* name of file containing this SFT */
  long offset;		/* SFT-offset with respect to a merged-SFT */
  UINT4 isft;           /* index of SFT this locator belongs to, used only in XLALLoadSFTs() */
};

typedef struct
{
  REAL8 version;
  INT4 gps_sec;
  INT4 gps_nsec;
  REAL8 tbase;
  INT4 first_frequency_index;
  INT4 nsamples;
} _SFT_header_v1_t;

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
  CHAR padding[2];
  INT4 comment_length;
} _SFT_header_v2_t;


/*---------- Global variables ----------*/
/* empty struct initializers */
static LALStatus empty_status;
const SFTConstraints empty_SFTConstraints;
const SFTCatalog empty_SFTCatalog;

/*---------- internal prototypes ----------*/
static LALStringVector *find_files (const CHAR *fpattern);

static void endian_swap(CHAR * pdata, size_t dsize, size_t nelements);
static int amatch(char *str, char *p);	/* glob pattern-matcher (public domain)*/
static BOOLEAN is_pattern(const char*c); /* filename string is a glob-style pattern */

static BOOLEAN is_valid_detector (const char *channel);
static BOOLEAN consistent_mSFT_header ( SFTtype header1, UINT4 version1, UINT4 nsamples1, SFTtype header2, UINT4 version2, UINT4 nsamples2 );
static BOOLEAN timestamp_in_list( LIGOTimeGPS timestamp, LIGOTimeGPSVector *list );
static long get_file_len ( FILE *fp );

static FILE * fopen_SFTLocator ( const struct tagSFTLocator *locator );
static BOOLEAN has_valid_v2_crc64 (FILE *fp );

static void lal_read_sft_bins_from_fp ( LALStatus *status, SFTtype **sft, UINT4 *binsread, UINT4 firstBin2read, UINT4 lastBin2read , FILE *fp );
static UINT4 read_sft_bins_from_fp ( SFTtype *ret, UINT4 *firstBinRead, UINT4 firstBin2read, UINT4 lastBin2read , FILE *fp );

static int read_sft_header_from_fp (FILE *fp, SFTtype  *header, UINT4 *version, UINT8 *crc64, BOOLEAN *swapEndian, CHAR **SFTcomment, UINT4 *numBins );
static int read_v2_header_from_fp ( FILE *fp, SFTtype *header, UINT4 *nsamples, UINT8 *header_crc64, UINT8 *ref_crc64, CHAR **SFTcomment, BOOLEAN swapEndian);
static int read_v1_header_from_fp ( FILE *fp, SFTtype *header, UINT4 *nsamples, BOOLEAN swapEndian);
static int compareSFTdesc(const void *ptr1, const void *ptr2);
static int compareSFTloc(const void *ptr1, const void *ptr2);
static int compareDetName(const void *ptr1, const void *ptr2);
static UINT8 calc_crc64(const CHAR *data, UINT4 length, UINT8 crc);
int read_SFTversion_from_fp ( UINT4 *version, BOOLEAN *need_swap, FILE *fp );

/*==================== FUNCTION DEFINITIONS ====================*/

/** Find the list of SFTs matching the \a file_pattern and satisfying the given \a constraints,
 * return an \c SFTCatalog of the matching SFTs.
 *
 * The optional \a constraints that can be specified are (type SFTConstraints)
 * - 'detector':      	which detector
 * - 'time-span':    	GPS start- and end-times
 * - 'timestamps':    	list of GPS start-times
 *
 *
 * ==> The returned SFTCatalog can be used directly as input to LALLoadSFTs()
 * to load a single-IFO SFTVector, or LALLoadMultiSFTs() to load a
 * multi-IFO vector of SFTVectors
 *
 * Except for the 'file_pattern' input, all the other constraints are optional
 * and can be passed as NULL (either globally constraings==NULL, or individually).
 *
 * Note that the constraints are combined by 'AND' and the resulting full constraint
 * MUST be satisfied (in particular: if 'timestamps' is given, all timestamps within
 * [startTime, endTime] MUST be found!.
 *
 * The returned SFTs in the catalogue are sorted by increasing GPS-epochs !
 *
 */
void
LALSFTdataFind (LALStatus *status,			/**< pointer to LALStatus structure */
		SFTCatalog **catalog,		/**< [out] SFT-catalogue of matching SFTs */
		const CHAR *file_pattern,	/**< which SFT-files */
		SFTConstraints *constraints	/**< additional constraints for SFT-selection */
		)
{
  LALStringVector *fnames;
  UINT4 i, numFiles, numSFTs;
  SFTCatalog *ret = NULL;
  SFTtype first_header = empty_SFTtype;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* ----- check input */
  ASSERT ( catalog, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( (*catalog) == NULL, status, SFTFILEIO_ENONULL, SFTFILEIO_MSGENONULL );

  ASSERT ( file_pattern, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );

  if ( constraints && constraints->detector )
    if ( strncmp(constraints->detector, "??", 2) && !is_valid_detector(constraints->detector) )
      {
	XLALPrintError( "\nInvalid detector-constraint '%s'\n\n", constraints->detector );
	ABORT ( status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL );
      }

  /* prepare return-catalog */
  if ( (ret = LALCalloc ( 1, sizeof ( SFTCatalog ))) == NULL ) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }

  /* find matching filenames */
  if ( (fnames = find_files (file_pattern)) == NULL) {
    XLALPrintError ("\nFailed to get filelist for pattern '%s'.\n\n", file_pattern);
    ABORT (status, SFTFILEIO_EGLOB, SFTFILEIO_MSGEGLOB);
  }
  numFiles = fnames->length;

  numSFTs = 0;
  /* ----- main loop: parse all matching files */
  for ( i = 0; i < numFiles; i ++ )
    {
      CHAR *fname = fnames->data[i];

      FILE *fp;
      long file_len;

      /* merged SFTs need to satisfy stronger consistency-constraints (-> see spec) */
      BOOLEAN mfirst_block = TRUE;
      UINT4   mprev_version = 0;
      SFTtype mprev_header;
      REAL8   mprev_nsamples = 0;


      if ( ( fp = LALFopen( fname, "rb" ) ) == NULL )
	{
	  XLALPrintError ( "\nFailed to open matched file '%s'\n\n", fname );
	  XLALDestroyStringVector ( fnames );
	  LALDestroySFTCatalog ( status->statusPtr, &ret );
	  ABORT( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
	}
      if ( (file_len = get_file_len(fp)) == 0 )
	{
	  XLALPrintError ( "\nGot file-len == 0 for '%s'\n\n", fname );
	  XLALDestroyStringVector ( fnames );
	  LALDestroySFTCatalog ( status->statusPtr, &ret );
	  fclose(fp);
	  ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
	}
      /* go through SFT-blocks in fp */
      while ( ftell(fp) < file_len )
	{
	  SFTtype this_header;
	  UINT4 this_version;
	  UINT4 this_nsamples;
	  UINT8 this_crc;
	  CHAR *this_comment = NULL;
	  BOOLEAN endian;
	  BOOLEAN want_this_block = FALSE;

	  long this_filepos;
	  if ( (this_filepos = ftell(fp)) == -1 )
	    {
	      XLALPrintError ("\nftell() failed for '%s'\n\n", fname );
	      XLALDestroyStringVector ( fnames );
	      LALDestroySFTCatalog ( status->statusPtr, &ret );
	      fclose (fp);
	      ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
	    }

	  if ( read_sft_header_from_fp (fp, &this_header, &this_version, &this_crc,
					&endian, &this_comment, &this_nsamples ) != 0 )
	    {
	      XLALPrintError ("\nERROR:File-block '%s:%ld' is not a valid SFT!\n\n", fname, ftell(fp));
	      XLALDestroyStringVector ( fnames );
	      if ( this_comment ) LALFree ( this_comment );
	      LALDestroySFTCatalog ( status->statusPtr, &ret );
	      fclose(fp);
	      ABORT ( status, SFTFILEIO_EHEADER, SFTFILEIO_MSGEHEADER );
	    }

	  /* if merged-SFT: check consistency constraints */
	  if ( !mfirst_block )
	    {
	      if ( ! consistent_mSFT_header ( mprev_header, mprev_version, mprev_nsamples,
					      this_header, this_version, this_nsamples ) )
		{
		  XLALPrintError ("merged SFT-file '%s' contains inconsistent SFT-blocks!\n\n", fname);
		  if ( this_comment ) LALFree ( this_comment );
		  XLALDestroyStringVector ( fnames );
		  LALDestroySFTCatalog ( status->statusPtr, &ret );
		  fclose(fp);
		  ABORT ( status, SFTFILEIO_EMERGEDSFT, SFTFILEIO_MSGEMERGEDSFT );
		}
	    } /* if !mfirst_block */

	  mprev_header = this_header;
	  mprev_version = this_version;
	  mprev_nsamples = this_nsamples;

	  want_this_block = TRUE;	/* default */
	  /* but does this SFT-block satisfy the user-constraints ? */
	  if ( constraints )
	    {
	      if ( constraints->detector && strncmp(constraints->detector, "??", 2) )
		{
		  /* v1-SFTs have '??' as detector-name */
		  if ( ! strncmp (this_header.name, "??", 2 ) )
		    strncpy ( this_header.name, constraints->detector, 2 );	/* SET to constraint! */
		  else if ( strncmp( constraints->detector, this_header.name, 2) )
		    want_this_block = FALSE;
		}

	      if ( constraints->startTime && ( GPS2REAL8(this_header.epoch) < GPS2REAL8( *constraints->startTime)))
		want_this_block = FALSE;

	      if ( constraints->endTime && ( GPS2REAL8(this_header.epoch) > GPS2REAL8( *constraints->endTime ) ) )
		want_this_block = FALSE;

	      if ( constraints->timestamps && !timestamp_in_list(this_header.epoch, constraints->timestamps) )
		want_this_block = FALSE;

	    } /* if constraints */

	  if ( want_this_block )
	    {
	      SFTDescriptor *desc;

	      numSFTs ++;

	      /* do we need to alloc more memory for the SFTs? */
	      if (  numSFTs > ret->length )
		{
		  UINT4 j;
		  /* we realloc SFT-memory blockwise in order to
		   * improve speed in debug-mode (using LALMalloc/LALFree)
		   */
		  ret->data = LALRealloc ( ret->data, (ret->length + SFTFILEIO_REALLOC_BLOCKSIZE) * sizeof( *(ret->data) ) );
		  if ( ret->data == NULL )
		    {
		      XLALPrintError("Memeory reallocation for SFTs failed: nSFT:%d, length:%d, add:%d\n",
				     numSFTs, ret->length, SFTFILEIO_REALLOC_BLOCKSIZE);
		      XLALDestroyStringVector ( fnames );
		      LALDestroySFTCatalog ( status->statusPtr, &ret );
		      if ( this_comment ) LALFree ( this_comment );
		      fclose(fp);
		      ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);
		    }
		  /* properly initialize data-fields pointers to NULL to avoid SegV when Freeing */
		  for ( j=0; j < SFTFILEIO_REALLOC_BLOCKSIZE; j ++ )
		    memset ( &(ret->data[ret->length + j]), 0, sizeof( ret->data[0] ) );

		  ret->length += SFTFILEIO_REALLOC_BLOCKSIZE;
		}

	      desc = &(ret->data[numSFTs - 1]);

	      desc->locator = LALCalloc ( 1, sizeof ( *(desc->locator) ) );
	      if ( desc->locator )
		desc->locator->fname = LALCalloc( 1, strlen(fname) + 1 );
	      if ( (desc->locator == NULL) || (desc->locator->fname == NULL ) )
		{
		  XLALDestroyStringVector ( fnames );
		  LALDestroySFTCatalog ( status->statusPtr, &ret );
		  fclose(fp);
		  ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);
		}
	      strcpy ( desc->locator->fname, fname );
	      desc->locator->offset = this_filepos;

	      desc->header  = this_header;
	      desc->comment = this_comment;
	      desc->numBins = this_nsamples;
	      desc->version = this_version;
	      desc->crc64   = this_crc;

	    } /* if want_this_block */
	  else
	    {
	      if ( this_comment ) LALFree ( this_comment );
	    }

	  /* seek to end of SFT data-entries in file  */
	  if ( fseek ( fp, this_nsamples * 8 , SEEK_CUR ) == -1 )
	    {
	      XLALPrintError ( "\nFailed to skip DATA field for SFT '%s': %s\n", fname, strerror(errno) );
	      if ( this_comment ) LALFree ( this_comment );
	      XLALDestroyStringVector ( fnames );
	      LALDestroySFTCatalog ( status->statusPtr, &ret );
	      fclose(fp);
	      ABORT ( status, SFTFILEIO_ESFTFORMAT, SFTFILEIO_MSGESFTFORMAT);
	    }

	  mfirst_block = FALSE;

	} /* while !feof */

      fclose(fp);

    } /* for i < numFiles */

  /* free matched filenames */
  XLALDestroyStringVector ( fnames );

  /* now realloc SFT-vector (alloc'ed blockwise) to its *actual size* */
  if ( (ret->data = LALRealloc ( ret->data, numSFTs * sizeof( *(ret->data) ) )) == NULL )
    {
      LALDestroySFTCatalog ( status->statusPtr, &ret );
      ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);
    }
  ret->length = numSFTs;

  /* ----- final consistency-checks: ----- */

  /* did we find exactly the timestamps that lie within [startTime, endTime]? */
  if ( constraints && constraints->timestamps )
    {
      LIGOTimeGPSVector *ts = constraints->timestamps;
      UINT4 numRequested = 0;
      REAL8 t0, t1;
      if ( constraints->startTime )
	t0 = GPS2REAL8 ( (*constraints->startTime) );
      else
	t0 = -1;
      if ( constraints->endTime )
	t1 = GPS2REAL8 ( (*constraints->endTime) );
      else
	t1 = LAL_REAL4_MAX;	/* large enough */

      for (i=0; i < ts->length; i ++ )
	{
	  REAL8 ti = GPS2REAL8(ts->data[i]);
	  if ( (t0 <= ti) && ( ti <= t1 ) )
	    numRequested ++;
	}
      if ( numRequested != ret->length )
	{
	  XLALPrintError ("\nERROR: found %s SFTs (%d) than given timestamps within [%f, %f] (%d)\n\n",
			  (ret->length < numRequested )?"fewer":"more",
			  ret->length, t0, t1, numRequested );
	  ABORT ( status, SFTFILEIO_ECONSTRAINTS, SFTFILEIO_MSGECONSTRAINTS );
	}
    } /* if constraints->timestamps */

  /* have all matched SFTs identical dFreq values ? */
  for ( i = 0; i < ret->length; i ++ )
    {
      SFTtype this_header = ret->data[i].header;

      if ( i == 0 )
	first_header = this_header;

      /* dont give out v1-SFTs without detector-entry, except if constraint->detector="??" ! */
      if ( !constraints || !constraints->detector || strncmp(constraints->detector, "??", 2) )
	{
	  if ( !strncmp ( this_header.name, "??", 2 ) )
	    {
	      LALDestroySFTCatalog ( status->statusPtr, &ret );
	      XLALPrintError ("\nERROR: '%s' matched v1-SFTs but no detector-constraint given!\n\n",
			      file_pattern);
	      ABORT ( status, SFTFILEIO_EDETECTOR, SFTFILEIO_MSGEDETECTOR );
	    }
	} /* if detector-constraint was not '??' */

      if ( this_header.deltaF != first_header.deltaF )
	{
	  LALDestroySFTCatalog ( status->statusPtr, &ret );
	  XLALPrintError("\nERROR: file-pattern '%s' matched SFTs with inconsistent deltaF: %.18g != %.18g!\n\n",
                         file_pattern, this_header.deltaF, first_header.deltaF );
	  ABORT ( status, SFTFILEIO_EDIFFTSFT, SFTFILEIO_MSGEDIFFTSFT );
	}

    } /* for i < numSFTs */


  /* sort catalog in order of increasing GPS-time */
  qsort( (void*)ret->data, ret->length, sizeof( ret->data[0] ), compareSFTdesc );


  /* return result catalog (=sft-vect and locator-vect) */
  (*catalog) = ret;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALSFTdataFind() */


/** Extract a timstamps-vector from the given SFTCatalog.
 *
 * \note A list of *unique* timestamps is returned, i.e. only a single copy of a timestamp
 * is returned, even if there are multiple occurrances of this timestamp in the catalog,
 * e.g. for multiple IFOs or multiple frequency-bands...
 *
 */
void
LALSFTtimestampsFromCatalog (LALStatus *status,			/**< pointer to LALStatus structure */
			     LIGOTimeGPSVector **timestamps,	/**< [out] extracted timestamps */
			     const SFTCatalog *catalog )	/**< input SFT-catalogue */
{
  UINT4 numSFTs, numTS;
  UINT4 i;
  LIGOTimeGPSVector *ret = NULL;
  REAL8 Tsft;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( timestamps, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( catalog, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( catalog->length > 0, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( *timestamps == NULL, status, SFTFILEIO_ENONULL, SFTFILEIO_MSGENONULL );

  numSFTs = catalog->length;

  Tsft = 1.0 / catalog->data[0].header.deltaF;

  /* large enough for all timestamps, but we might have fewer than that */
  TRY ( LALCreateTimestampVector ( status->statusPtr, &ret, numSFTs ), status );

  numTS = 0;
  for ( i=0; i < numSFTs; i ++ )
    {
      LIGOTimeGPS *thisTs = &(catalog->data[i].header.epoch);
      if ( i && (thisTs->gpsSeconds == ret->data[i-1].gpsSeconds) && (thisTs->gpsNanoSeconds == ret->data[i-1].gpsNanoSeconds) )
	continue;	/* skip double timestamp */

      ret->data[numTS++] = catalog->data[i].header.epoch;	/* new timestamp */

    } /* for i < numSFTs */

  /* realloc data-segment in timestamps-vector to the actual length */
  if ( (ret->data = LALRealloc ( ret->data, numTS * sizeof ( *ret->data ) )) == NULL ) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }
  ret->length = numTS;
  ret->deltaT = Tsft;

  /* done: return Ts-vector */
  (*timestamps) = ret;

  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* LALTimestampsFromSFTCatalog() */



/*
   This function reads an SFT (segment) from an open file pointer into a buffer.
   firstBin2read specifies the first bin to read from the SFT, lastBin2read is the last bin.
   If the SFT contains fewer bins than specified, all bins from the SFT are read.
   The function returns the last bin actually read, firstBinRead
   is set to the first bin actually read. In case of an error, 0 is returned
   and firstBinRead is set to a code further decribing the error condition.
*/
static UINT4
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
    if ( read_sft_header_from_fp (fp, ret, &version, &crc64, &swapEndian, NULL, &numSFTbins ) != 0 )
      {
	XLALPrintError ("read_sft_bins_from_fp(): Failed to read SFT-header!\n");
	*firstBinRead = 2;
	return(0);
      }
    ret->data = data;
  }

  tmp = ret->f0 / ret->deltaF;
  firstSFTbin = MYROUND ( tmp );
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
  if ( version == 1 || swapEndian )
    {
      UINT4 i;
      REAL8 band = 1.0 * numSFTbins * ret->deltaF;/* need the TOTAL frequency-band in the SFT-file! */
      REAL8 fsamp = 2.0 * band;
      REAL8 dt = 1.0 / fsamp;

      for ( i=0; i < numBins2read; i ++ )
	{
	  REAL4 *rep, *imp;

	  rep = &(ret->data->data[i].re);
	  imp = &(ret->data->data[i].im);

	  if ( swapEndian )
	    {
	      endian_swap( (CHAR *) rep, sizeof ( *rep ), 1 );
	      endian_swap( (CHAR *) imp, sizeof ( *imp ), 1 );
	    }

	  /* if the SFT-file was in v1-Format: need to renormalize the data now by 'Delta t'
	   * in order to follow the correct SFT-normalization
	   * (see LIGO-T040164-01-Z, and LIGO-T010095-00)
	   */
	  if ( version == 1 )
	    {
	      (*rep) *= dt;
	      (*imp) *= dt;
	    }
	} /* for i < numBins2read */
    } /* if SFT-v1 */

  /* return last bin read */
  return(lastBin2read);

} /* read_sft_bins_from_fp() */


/** segments read so far from one SFT */
typedef struct {
  UINT4 first;                     /**< first bin in this segment */
  UINT4 last;                      /**< last bin in this segment */
  LIGOTimeGPS epoch;               /**< timestamp of this SFT */
  struct tagSFTLocator *lastfrom;  /**< last bin read from this locator */
} SFTReadSegment;


/** Load the given frequency-band <tt>[fMin, fMax]</tt> (inclusively) from the SFT-files listed in the
 * SFT-'catalogue' ( returned by LALSFTdataFind() ).
 *
 * Note: \a fMin (or \a fMax) is allowed to be set to \c -1, which means to read in all
 * Frequency-bins from the lowest (or up to the highest) found in all SFT-files of the catalog.
 *
 * Note 2: The returned frequency-interval is guaranteed to contain <tt>[fMin, fMax]</tt>,
 * but is allowed to be larger, as it must be an interval of discrete frequency-bins as found
 * in the SFT-file.
 *
 * Note 3: This function has the capability to read sequences of (v2-)SFT segments and
 * putting them together to single SFTs while reading.
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
  minbin = firstbin = MYROUND( catalog->data[0].header.f0 / deltaF );
  maxbin = lastbin = firstbin + catalog->data[0].numBins - 1;
  for(catPos = 1; catPos < catalog->length; catPos++) {
    firstbin = MYROUND( catalog->data[catPos].header.f0 / deltaF );
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
  LogPrintf(LOG_DETAIL, "XLALLoadSFTs(): fMin: %f, fMax: %f, deltaF: %f, minbin: %u, maxbin: %u\n",
	    fMin, fMax, deltaF, minbin, maxbin);

  /* calculate first and last frequency bin to read */
  if (fMin < 0)
    firstbin = minbin;
  else
    firstbin = floor (fMin / deltaF);
  if (fMax < 0)
    lastbin = maxbin;
  else {
    lastbin = ceil (fMax / deltaF);
    if((lastbin == 0) && (fMax != 0)) {
      XLALPrintError("ERROR: last bin to read is 0 (fMax: %f, deltaF: %f)\n", fMax, deltaF);
      XLALLOADSFTSERROR(XLAL_EINVAL);
    }
  }
  LogPrintf(LOG_DETAIL, "XLALLoadSFTs(): Reading from first bin: %u, last bin: %u\n", firstbin, lastbin);

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
      UINT4 firstSFTbin = MYROUND ( tmp );
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
	LogPrintf(LOG_DETAIL, "XLALLoadSFTs(): Opening file '%s'\n", fname);
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
      LogPrintf(LOG_DETAIL, "XLALLoadSFTs(): Read data from %s:%lu: %u - %u\n",
		locator->fname, locator->offset, firstBinRead, lastBinRead);
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
	if(locatalog.data[catPos].header.name && *(locatalog.data[catPos].header.name))
	  strcpy(sftVector->data[isft].name, locatalog.data[catPos].header.name);
	sftVector->data[isft].sampleUnits = locatalog.data[catPos].header.sampleUnits;
	memcpy(sftVector->data[isft].data->data + (firstBinRead - firstbin),
	       thisSFT->data->data,
	       (lastBinRead - firstBinRead + 1) * sizeof(COMPLEX8));

      } else if(!firstBinRead) {
	/* no needed data had been in this segment */
	LogPrintf(LOG_DETAIL, "XLALLoadSFTs(): No data read from %s:%lu\n", locator->fname, locator->offset);

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



/** Function to load a catalog of SFTs from possibly different detectors.
    This is similar to LALLoadSFTs except that the input SFT catalog is
    allowed to contain multiple ifos. The output is the structure
    MultiSFTVector which is a vector of (pointers to) SFTVectors, one for
    each ifo found in the catalog. As in LALLoadSFTs, fMin and fMax can be
    set to -1 to get the full SFT from the lowest to the highest frequency
    bin found in the SFT.
    *
    * output SFTvectors are sorted alphabetically by detector-name
    *
 */
MultiSFTVector *
XLALLoadMultiSFTs (const SFTCatalog *inputCatalog,   /**< The 'catalogue' of SFTs to load */
		   REAL8 fMin,		             /**< minumum requested frequency (-1 = read from lowest) */
		   REAL8 fMax		             /**< maximum requested frequency (-1 = read up to highest) */
		   )

{
  UINT4 k, j, i, length;
  UINT4 numifo=0; /* number of ifos */
  UINT4 numifoMax, numifoMaxNew; /* for memory allocation purposes */
  CHAR  *name=NULL;
  CHAR  **ifolist=NULL; /* list of ifo names */
  UINT4  *numsfts=NULL; /* number of sfts for each ifo */
  UINT4 **sftLocationInCatalog=NULL; /* location of sfts in catalog for each ifo */
  SFTCatalog **catalog=NULL;
  MultiSFTVector *multSFTVec=NULL;

  if(!inputCatalog || !(inputCatalog->length)) {
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  length = inputCatalog->length;
  if ( (name = (CHAR *)XLALCalloc(3, sizeof(CHAR))) == NULL ) {
    XLAL_ERROR_NULL(XLAL_ENOMEM);
  }

  /* the number of ifos can be at most equal to length */
  /* each ifo name is 2 characters + \0 */

  numifoMax = 3; /* should be sufficient -- realloc used later in case required */

  if ( (ifolist = (CHAR **)XLALCalloc( 1, numifoMax * sizeof(CHAR *))) == NULL) {
    XLAL_ERROR_NULL(XLAL_ENOMEM);
  }
  if ( (sftLocationInCatalog = (UINT4 **)XLALCalloc( 1, numifoMax * sizeof(UINT4 *))) == NULL) {
    XLAL_ERROR_NULL(XLAL_ENOMEM);
  }
  if ( (numsfts = (UINT4 *)XLALCalloc( 1, numifoMax * sizeof(UINT4))) == NULL) {
    XLAL_ERROR_NULL(XLAL_ENOMEM);
  }

  for ( k = 0; k < numifoMax; k++) {
    if ( (ifolist[k] = (CHAR *)XLALCalloc( 1, 3*sizeof(CHAR))) == NULL) {
      XLAL_ERROR_NULL(XLAL_ENOMEM);
    }
    if ( (sftLocationInCatalog[k] = (UINT4 *)XLALCalloc( 1, length*sizeof(UINT4))) == NULL) {
      XLAL_ERROR_NULL(XLAL_ENOMEM);
    }
  }

  /* loop over sfts in catalog and look at ifo names and
     find number of different ifos and number of sfts for each ifo
     Also find location of sft in catalog  */
  for ( k = 0; k < length; k++)
    {
      strncpy( name, inputCatalog->data[k].header.name, 3 );

      /* go through list of ifos till a match is found or list is exhausted */
      for ( j = 0; ( j < numifo ) && strncmp( name, ifolist[j], 3); j++ )
	;

      if ( j < numifo )
	{
	  /* match found with jth existing ifo */
	  sftLocationInCatalog[j][ numsfts[j] ] = k;
	  numsfts[j]++;
	}
      else
	{
	  /* add ifo to list of ifos */

	  /* first check if number of ifos is larger than numifomax */
	  /* and realloc if necessary */
	  if ( numifo >= numifoMax )
	    {
	      numifoMaxNew = numifoMax + 3;
	      if ( (ifolist = (CHAR **)XLALRealloc( ifolist, numifoMaxNew * sizeof(CHAR *))) == NULL) {
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	      }
	      if ( (sftLocationInCatalog = (UINT4 **)XLALRealloc( sftLocationInCatalog, numifoMaxNew * sizeof(UINT4 *))) == NULL) {
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	      }
	      if ( (numsfts = (UINT4 *)XLALRealloc( numsfts, numifoMaxNew * sizeof(UINT4))) == NULL) {
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	      }

	      for ( i = numifoMax; i < numifoMaxNew; i++) {
		if ( (ifolist[i] = (CHAR *)XLALCalloc( 1, 3*sizeof(CHAR))) == NULL) {
		  XLAL_ERROR_NULL(XLAL_ENOMEM);
		}
		if ( (sftLocationInCatalog[i] = (UINT4 *)XLALCalloc( 1, length*sizeof(UINT4))) == NULL) {
		  XLAL_ERROR_NULL(XLAL_ENOMEM);
		}
	      } /* loop from numifoMax to numifoMaxNew */

	      /* reset numifoMax */
	      numifoMax = numifoMaxNew;
	    } /* if ( numifo >= numifoMax) -- end of realloc */

	  strncpy( ifolist[numifo], name, 3);
	  sftLocationInCatalog[j][0] = k;
	  numsfts[numifo] = 1;
	  numifo++;

	} /* else part of if ( j < numifo ) */
    } /*  for ( k = 0; k < length; k++) */

  /* now we can create the catalogs */
  if ( (catalog = (SFTCatalog **)XLALCalloc( numifo, sizeof(SFTCatalog *))) == NULL) {
    XLAL_ERROR_NULL(XLAL_ENOMEM);
  }

  for ( j = 0; j < numifo; j++)
    {
      if ( (catalog[j] = (SFTCatalog *)XLALCalloc(1, sizeof(SFTCatalog))) == NULL) {
	XLAL_ERROR_NULL(XLAL_ENOMEM);
      }
      catalog[j]->length = numsfts[j];
      if ( (catalog[j]->data = (SFTDescriptor *)XLALCalloc( numsfts[j], sizeof(SFTDescriptor))) == NULL) {
	XLAL_ERROR_NULL(XLAL_ENOMEM);
      }

      for ( k = 0; k < numsfts[j]; k++)
	{
	  UINT4 location = sftLocationInCatalog[j][k];
	  catalog[j]->data[k] = inputCatalog->data[location];
	}
    }

  /* create multi sft vector */
  if ( (multSFTVec = (MultiSFTVector *)XLALCalloc(1, sizeof(MultiSFTVector))) == NULL){
    XLAL_ERROR_NULL(XLAL_ENOMEM);
  }
  multSFTVec->length = numifo;

  if ( (multSFTVec->data = (SFTVector **)XLALCalloc(numifo, sizeof(SFTVector *))) == NULL) {
    XLAL_ERROR_NULL(XLAL_ENOMEM);
  }
  for ( j = 0; j < numifo; j++) {
    if( ! ( multSFTVec->data[j] = XLALLoadSFTs ( catalog[j], fMin, fMax ) ) )
    {
      /* free sft vectors created previously in loop */
      for ( i = 0; (INT4)i < (INT4)j-1; i++)
	XLALDestroySFTVector(multSFTVec->data[i]);
      XLALFree(multSFTVec->data);
      XLALFree(multSFTVec);

      /* also free catalog and other memory allocated earlier */
      for ( i = 0; i < numifo; i++) {
	XLALFree(catalog[i]->data);
	XLALFree(catalog[i]);
      }
      XLALFree( catalog);

      for ( i = 0; i < numifoMax; i++) {
	XLALFree(ifolist[i]);
	XLALFree(sftLocationInCatalog[i]);
      }
      XLALFree(ifolist);
      XLALFree(sftLocationInCatalog);

      XLALFree(numsfts);
      XLALFree(name);

      XLAL_ERROR_NULL ( XLAL_EFUNC );

    }
  }

  /* sort final multi-SFT vector by detector-name */
  qsort ( multSFTVec->data, multSFTVec->length, sizeof( multSFTVec->data[0] ), compareDetName );

  /* free memory and exit */
  for ( j = 0; j < numifo; j++) {
    XLALFree(catalog[j]->data);
    XLALFree(catalog[j]);
  }
  XLALFree(catalog);

  for ( k = 0; k < numifoMax; k++) {
    XLALFree(ifolist[k]);
    XLALFree(sftLocationInCatalog[k]);
  }
  XLALFree(ifolist);
  XLALFree(sftLocationInCatalog);

  XLALFree(numsfts);
  XLALFree(name);

  return(multSFTVec);

} /* XLALLoadMultiSFTs() */



/** Load the given frequency-band <tt>[fMin, fMax]</tt> (inclusively) from the SFT-files listed in the
 * SFT-'catalogue' ( returned by LALSFTdataFind() ).
 *
 * Note: \a fMin (or \a fMax) is allowed to be set to \c -1, which means to read in all
 * Frequency-bins from the lowest (or up to the highest) found in the SFT-file.
 *
 * Note 2: The returned frequency-interval is guaranteed to contain <tt>[fMin, fMax]</tt>,
 * but is allowed to be larger, as it must be an interval of discrete frequency-bins as found
 * in the SFT-file.
 *
 * Note 3: This function has the capability to read sequences of (v2-)SFT segments and
 * putting them together to single SFTs while reading.
*/
void
LALLoadSFTs ( LALStatus *status,	/**< pointer to LALStatus structure */
	      SFTVector **outsfts,	   /**< [out] vector of read-in SFTs */
	      const SFTCatalog *catalog,   /**< The 'catalogue' of SFTs to load */
	      REAL8 fMin,		   /**< minumum requested frequency (-1 = read from lowest) */
	      REAL8 fMax		   /**< maximum requested frequency (-1 = read up to highest) */
	      )
{
  UINT4 catFile = 0;           /* current file in catalog */
  LIGOTimeGPS epoch;           /* current timestamp */
  UINT4 firstbin, lastbin;     /* the first and last bin we want to read */
  UINT4 nextbin;               /* the Bin we expect to read next */
  UINT4 binsread;              /* number of bins actually read from an sft */
  REAL8 deltaF;
  SFTtype *onesft = NULL;      /* a single SFT that was read */
  COMPLEX8Vector *sftbins;     /* bins of the SFT that is constructed */
  SFTVector *sfts;             /* the SFTVector to return */
  FILE *fp;                    /* filepointer to read an SFT from */
  UINT4 i;                     /* loop counter */
  UINT4 firstInSFT, lastInSFT=0; /* first and last bin in current SFT */

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( outsfts, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( *outsfts == NULL, status, SFTFILEIO_ENONULL, SFTFILEIO_MSGENONULL );
  ASSERT ( catalog, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( fMin <= fMax, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL );

  sfts = (SFTVector*)LALMalloc(sizeof(SFTVector));
  if (!sfts) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }
  sfts->length = 0;
  sfts->data = NULL;

#ifdef SFTFILEIO_DEBUG
  fprintf(stderr, ": catalog has %u files\n", catalog->length);
  for(i=0; i < catalog->length; i++)
    fprintf(stderr, ": %s\n", XLALshowSFTLocator ( catalog->data[i].locator ) );
#endif

  /* while there are files in catalog */
  while (catFile < catalog->length)
    {
      /* calculate first and last frequency bin to read */
      /* the patch for fMin/fMax == -1 should work as it did before
	 with 'single' SFT files. However with segmented SFT files,
	 the function reads only the frequency range of the
	 FIRST SEGMENT
      */
      deltaF = catalog->data[catFile].header.deltaF; /* Hz/bin */
      if (fMin < 0)
	firstbin = MYROUND( catalog->data[catFile].header.f0 / deltaF );
      else
	firstbin = floor(fMin / deltaF);
      if (fMax < 0)
	lastbin = firstbin + catalog->data[catFile].numBins - 1;
      else
	lastbin = ceil (fMax / deltaF);
      nextbin = firstbin;

      /* get first timestamp */
      epoch = catalog->data[catFile].header.epoch;

#ifdef SFTFILEIO_DEBUG
      fprintf(stderr,"+ firstbin: %u, lastbin: %u\n", firstbin, lastbin);
#endif

      /* allocate space for the frequency bins */
      sftbins = XLALCreateCOMPLEX8Vector(lastbin-firstbin+1);
      if (!sftbins) {
	LALDestroySFTVector (status->statusPtr, &sfts);
	ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
      }
      /* add a new SFT to the output SFTVector */
      sfts->data = (SFTtype*)LALRealloc(sfts->data, (sfts->length+1)*sizeof(SFTtype));
      if (!sfts->data) {
	LALFree(sfts);
	ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
      }
      /* initialize metadata, copying the header information from the first SFT for defaults */
      sfts->data[sfts->length] = catalog->data[catFile].header;;
      sfts->data[sfts->length].f0 = firstbin * deltaF;
      /* attach the bin space to the sft vector */
      sfts->data[sfts->length].data = sftbins;
      /* vector lenth has increased */
      sfts->length++;

      /* while there are files with this timestamp */
      while ((catFile < catalog->length) &&
	     (GPSEQUAL(epoch,catalog->data[catFile].header.epoch)))
	{
	  /* deltaF consistency check */
	  if ( deltaF != catalog->data[catFile].header.deltaF ) {
	    XLALPrintError ( "Frequency spacing dosn't match in SFT '%s (%f %f)\n",
			    XLALshowSFTLocator ( catalog->data[catFile].locator ),
			    deltaF, catalog->data[catFile].header.deltaF );
	    LALDestroySFTVector (status->statusPtr, &sfts);
	    ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
	  }

	  firstInSFT = MYROUND( catalog->data[catFile].header.f0 / deltaF );
	  lastInSFT  = firstInSFT + catalog->data[catFile].numBins - 1;

#ifdef SFTFILEIO_DEBUG
	  fprintf(stderr, "- nextbin: %u, firstInSFT: %u, lastInSFT: %u, %s\n",
		  nextbin, firstInSFT, lastInSFT,
		  XLALshowSFTLocator ( catalog->data[catFile].locator ) );
#endif

	  /* if fmax is in this SFT, this is either a single SFT or the last in a sequence */
	  if (( lastbin >= firstInSFT) &&
	      ( lastbin <= lastInSFT)) {

	    /* issue an error if the file neither contains firstbin nor starts with nextbin */
	    if (( firstbin < firstInSFT ) &&
		( nextbin != firstInSFT )) {
	      XLALPrintError ( "Starting frequency %f not contained in SFT '%s'\n"
			      "   (or sequence broken at this last file)\n",
			      fMin, XLALshowSFTLocator ( catalog->data[catFile].locator ) );
	      LALDestroySFTVector (status->statusPtr, &sfts);
	      ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
	    }

	    /* read from nextbin to lastbin */
	    if ( (fp = fopen_SFTLocator ( catalog->data[catFile].locator )) == NULL ) {
	      XLALPrintError ( "Failed to open locator '%s'\n",
			      XLALshowSFTLocator ( catalog->data[catFile].locator ) );
	      LALDestroySFTVector (status->statusPtr, &sfts);
	      ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
	    }
	    lal_read_sft_bins_from_fp (status->statusPtr, &onesft, &binsread, nextbin, lastbin, fp);
	    fclose(fp);
	    if ( status->statusPtr->statusCode ) {
	      XLALPrintError ( "Failed to read from locator '%s'\n",
			      XLALshowSFTLocator ( catalog->data[catFile].locator ) );
	      LALDestroySFTVector (status->statusPtr, &sfts);
	      ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
	    }
	    /* insert read bins into the sft vector to retun */
	    for(i=0; i < binsread; i++)
	      sftbins->data[nextbin - firstbin + i] = onesft->data->data[i];
	    LALDestroySFTtype(status->statusPtr,&onesft);

	    /* skip remaining catalog files with same timestamp (must have higher frequency) */
	    while ((catFile < catalog->length) &&
		   (GPSEQUAL(epoch,catalog->data[catFile].header.epoch)))
	      catFile++;



	  /* if fmin is in this SFT, this SFT starts a sequence (single SFT was previous case) */
	  } else if (( firstbin >= firstInSFT ) &&
		     ( firstbin <= lastInSFT  )) {
	    /* read from firstbin to end */
	    if ( (fp = fopen_SFTLocator ( catalog->data[catFile].locator )) == NULL ) {
	      XLALPrintError ( "Failed to open locator '%s'\n",
			      XLALshowSFTLocator ( catalog->data[catFile].locator ) );
	      LALDestroySFTVector (status->statusPtr, &sfts);
	      ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
	    }
	    lal_read_sft_bins_from_fp (status->statusPtr, &onesft, &binsread, firstbin, lastbin, fp);
	    fclose(fp);
	    if ( status->statusPtr->statusCode ) {
	      XLALPrintError ( "Failed to read from locator '%s'\n",
			      XLALshowSFTLocator ( catalog->data[catFile].locator ) );
	      LALDestroySFTVector (status->statusPtr, &sfts);
	      ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
	    }
	    /* insert read bins into the sft vector to retun */
	    for(i=0; i < binsread; i++)
	      sftbins->data[i] = onesft->data->data[i];
	    LALDestroySFTtype(status->statusPtr,&onesft);

	    /* read from nextbin on in next file */
	    nextbin += binsread;
	    catFile++;



          /* if this SFT starts with nextbin, it is the next SFT we expect to read */
	  } else if ( nextbin == firstInSFT) {
	    /* read whole file */
	    if ( (fp = fopen_SFTLocator ( catalog->data[catFile].locator )) == NULL ) {
	      XLALPrintError ( "Failed to open locator '%s'\n",
			      XLALshowSFTLocator ( catalog->data[catFile].locator ) );
	      LALDestroySFTVector (status->statusPtr, &sfts);
	      ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
	    }
	    lal_read_sft_bins_from_fp (status->statusPtr, &onesft, &binsread, nextbin, lastbin, fp);

	    fclose(fp);
	    if ( status->statusPtr->statusCode ) {
	      XLALPrintError ( "Failed to read from locator '%s'\n",
			      XLALshowSFTLocator ( catalog->data[catFile].locator ) );
	      LALDestroySFTVector (status->statusPtr, &sfts);
	      ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
	    }
	    /* insert read bins into the sft vector to retun */
	    for(i=0; i < binsread; i++)
	      sftbins->data[nextbin - firstbin + i] = onesft->data->data[i];
	    LALDestroySFTtype(status->statusPtr,&onesft);

	    /* read from nextbin on in next file */
	    nextbin += binsread;
	    catFile++;



	  /* if the frequency range of this file is completely lower than fmin, skip this file */
	  } else if ( firstbin > lastInSFT ) {
	    /* skip this file */
	    catFile++;



	  /* if none of the above applies, something must be wrong with this sequence */
	  } else {
	    XLALPrintError ( "Error in SFT sequence at locator '%s'\n"
			    "  (expected bin:%u, SFT has bins %u to %u)\n",
			    XLALshowSFTLocator ( catalog->data[catFile].locator ),
			    nextbin, firstInSFT, lastInSFT);
	    ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
	  }


	} /* while there are files with this timestamp */


      /* last check: did we find fMax at all? */
      if ( lastbin > lastInSFT) {
	XLALPrintError ( "Ending frequency %f not contained in SFT (sequence)\n"
			"  (expected bin: %u, last bin in SFT sequence: %u)\n",
			fMax, lastbin, lastInSFT);
	ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
      }


    } /* while files in catalog */

  /* assign output value */
  *outsfts = sfts;
  /* cleanup & return */
  DETATCHSTATUSPTR (status);
  RETURN(status);
} /* LALLoadSFTs() */


/** Function to load a catalog of SFTs from possibly different detectors.
    This is similar to LALLoadSFTs except that the input SFT catalog is
    allowed to contain multiple ifos.  The output is the structure
    MultiSFTVector which is a vector of (pointers to) SFTVectors, one for
    each ifo found in the catalog.   As in LALLoadSFTs, fMin and fMax can be
    set to -1 to get the full SFT from the lowest to the highest frequency
    bin found in the SFT.
    *
    * output SFTvectors are sorted alphabetically by detector-name
    *
 */
void LALLoadMultiSFTs ( LALStatus *status,		/**< pointer to LALStatus structure */
			MultiSFTVector **out,             /**< [out] vector of read-in SFTs -- one sft vector for each ifo found in catalog*/
			const SFTCatalog *inputCatalog,   /**< The 'catalogue' of SFTs to load */
			REAL8 fMin,		          /**< minumum requested frequency (-1 = read from lowest) */
			REAL8 fMax		          /**< maximum requested frequency (-1 = read up to highest) */
			)

{
  UINT4 k, j, i, length;
  UINT4 numifo=0; /* number of ifos */
  UINT4 numifoMax, numifoMaxNew; /* for memory allocation purposes */
  CHAR  *name=NULL;
  CHAR  **ifolist=NULL; /* list of ifo names */
  UINT4  *numsfts=NULL; /* number of sfts for each ifo */
  UINT4 **sftLocationInCatalog=NULL; /* location of sfts in catalog for each ifo */
  SFTCatalog **catalog=NULL;
  MultiSFTVector *multSFTVec=NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( out, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( *out == NULL, status, SFTFILEIO_ENONULL, SFTFILEIO_MSGENONULL );
  ASSERT ( inputCatalog, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( inputCatalog->length, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL );

  length = inputCatalog->length;
  if ( (name = (CHAR *)LALCalloc(3, sizeof(CHAR))) == NULL ) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }

  /* the number of ifos can be at most equal to length */
  /* each ifo name is 2 characters + \0 */

  numifoMax = 3; /* should be sufficient -- realloc used later in case required */

  if ( (ifolist = (CHAR **)LALCalloc( 1, numifoMax * sizeof(CHAR *))) == NULL) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }
  if ( (sftLocationInCatalog = (UINT4 **)LALCalloc( 1, numifoMax * sizeof(UINT4 *))) == NULL) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }
  if ( (numsfts = (UINT4 *)LALCalloc( 1, numifoMax * sizeof(UINT4))) == NULL) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }

  for ( k = 0; k < numifoMax; k++) {
    if ( (ifolist[k] = (CHAR *)LALCalloc( 1, 3*sizeof(CHAR))) == NULL) {
      ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
    }
    if ( (sftLocationInCatalog[k] = (UINT4 *)LALCalloc( 1, length*sizeof(UINT4))) == NULL) {
      ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
    }
  }

  /* loop over sfts in catalog and look at ifo names and
     find number of different ifos and number of sfts for each ifo
     Also find location of sft in catalog  */
  for ( k = 0; k < length; k++)
    {
      strncpy( name, inputCatalog->data[k].header.name, 3 );

      /* go through list of ifos till a match is found or list is exhausted */
      for ( j = 0; ( j < numifo ) && strncmp( name, ifolist[j], 3); j++ )
	;

      if ( j < numifo )
	{
	  /* match found with jth existing ifo */
	  sftLocationInCatalog[j][ numsfts[j] ] = k;
	  numsfts[j]++;
	}
      else
	{
	  /* add ifo to list of ifos */

	  /* first check if number of ifos is larger than numifomax */
	  /* and realloc if necessary */
	  if ( numifo >= numifoMax )
	    {
	      numifoMaxNew = numifoMax + 3;
	      if ( (ifolist = (CHAR **)LALRealloc( ifolist, numifoMaxNew * sizeof(CHAR *))) == NULL) {
		ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
	      }
	      if ( (sftLocationInCatalog = (UINT4 **)LALRealloc( sftLocationInCatalog, numifoMaxNew * sizeof(UINT4 *))) == NULL) {
		ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
	      }
	      if ( (numsfts = (UINT4 *)LALRealloc( numsfts, numifoMaxNew * sizeof(UINT4))) == NULL) {
		ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
	      }

	      for ( i = numifoMax; i < numifoMaxNew; i++) {
		if ( (ifolist[i] = (CHAR *)LALCalloc( 1, 3*sizeof(CHAR))) == NULL) {
		  ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
		}
		if ( (sftLocationInCatalog[i] = (UINT4 *)LALCalloc( 1, length*sizeof(UINT4))) == NULL) {
		  ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
		}
	      } /* loop from numifoMax to numifoMaxNew */

	      /* reset numifoMax */
	      numifoMax = numifoMaxNew;
	    } /* if ( numifo >= numifoMax) -- end of realloc */

	  strncpy( ifolist[numifo], name, 3);
	  sftLocationInCatalog[j][0] = k;
	  numsfts[numifo] = 1;
	  numifo++;

	} /* else part of if ( j < numifo ) */
    } /*  for ( k = 0; k < length; k++) */

  /* now we can create the catalogs */
  if ( (catalog = (SFTCatalog **)LALCalloc( numifo, sizeof(SFTCatalog *))) == NULL) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }

  for ( j = 0; j < numifo; j++)
    {
      if ( (catalog[j] = (SFTCatalog *)LALCalloc(1, sizeof(SFTCatalog))) == NULL) {
	ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
      }
      catalog[j]->length = numsfts[j];
      if ( (catalog[j]->data = (SFTDescriptor *)LALCalloc( numsfts[j], sizeof(SFTDescriptor))) == NULL) {
	ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
      }

      for ( k = 0; k < numsfts[j]; k++)
	{
	  UINT4 location = sftLocationInCatalog[j][k];
	  catalog[j]->data[k] = inputCatalog->data[location];
	}
    }

  /* create multi sft vector */
  if ( (multSFTVec = (MultiSFTVector *)LALCalloc(1, sizeof(MultiSFTVector))) == NULL){
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }
  multSFTVec->length = numifo;

  if ( (multSFTVec->data = (SFTVector **)LALCalloc(numifo, sizeof(SFTVector *))) == NULL) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }
  for ( j = 0; j < numifo; j++) {
#ifdef USEXLALLOADSFTS
    if( ! ( multSFTVec->data[j] = XLALLoadSFTs ( catalog[j], fMin, fMax ) ) )
#else
    LALLoadSFTs ( status->statusPtr, multSFTVec->data + j, catalog[j], fMin, fMax );
    BEGINFAIL ( status )
#endif
    {
      /* free sft vectors created previously in loop */
      for ( i = 0; (INT4)i < (INT4)j-1; i++)
	LALDestroySFTVector ( status->statusPtr, multSFTVec->data + i);
      LALFree(multSFTVec->data);
      LALFree(multSFTVec);

      /* also free catalog and other memory allocated earlier */
      for ( i = 0; i < numifo; i++) {
	LALFree(catalog[i]->data);
	LALFree(catalog[i]);
      }
      LALFree( catalog);

      for ( i = 0; i < numifoMax; i++) {
	LALFree(ifolist[i]);
	LALFree(sftLocationInCatalog[i]);
      }
      LALFree(ifolist);
      LALFree(sftLocationInCatalog);

      LALFree(numsfts);
      LALFree(name);

#ifdef USEXLALLOADSFTS
      ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
#endif
    }
#ifndef USEXLALLOADSFTS
    ENDFAIL ( status );
#endif
  }

  /* sort final multi-SFT vector by detector-name */
  qsort ( multSFTVec->data, multSFTVec->length, sizeof( multSFTVec->data[0] ), compareDetName );

  /* free memory and exit */
  for ( j = 0; j < numifo; j++) {
    LALFree(catalog[j]->data);
    LALFree(catalog[j]);
  }
  LALFree( catalog);

  for ( k = 0; k < numifoMax; k++) {
    LALFree(ifolist[k]);
    LALFree(sftLocationInCatalog[k]);
  }
  LALFree(ifolist);
  LALFree(sftLocationInCatalog);

  LALFree(numsfts);
  LALFree(name);

  *out = multSFTVec;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALLoadMultiSFTs() */


/** Function to check validity of SFTs listed in catalog.
 * This function simply reads in those SFTs and checks their CRC64 checksum, which
 * is the only check that has not yet been done by the operations up to this point.
 *
 * Returns the LAL-return code of a failure (or 0 on success) in 'check_result'.
 *
 * \note: because this function has to read the complete SFT-data into memory for the
 *  whole set of matching SFTs, it is potentially slow and memory-intensive.
 *
 * This function will NOT fail if one of the SFT-operations fails, instead it returns
 * the status-code of this failure in 'check_result'. The function will fail, however,
 * on invalid input.
 *
 */
void
LALCheckSFTs ( LALStatus *status,			/**< pointer to LALStatus structure */
	       INT4 *check_result, 	     /**< LAL-status of SFT-operations */
	       const CHAR *file_pattern,     /**< where to find the SFTs: normally a path+file-pattern */
	       SFTConstraints *constraints   /**< additional constraints for SFT-selection */
	       )
{
  LALStatus sft_status = empty_status;
  SFTCatalog *catalog = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( check_result, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( file_pattern, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );

  if ( constraints && constraints->detector && ! is_valid_detector(constraints->detector) )
    {
      XLALPrintError( "\nInvalid detector-constraint '%s'\n\n", constraints->detector );
      ABORT ( status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL );
    }

  /* Step 1: find the catalog of matching SFTs */
  LALSFTdataFind ( &sft_status, &catalog, file_pattern, constraints );
  if ( ((*check_result) = sft_status.statusCode) == 0 ) {

    /* Step 2: step through SFTs and check CRC64 */
    if ( catalog ) {
      TRY ( LALCheckSFTCatalog ( status->statusPtr, check_result, catalog ), status );
    }
  }

  if ( catalog ) {
    TRY ( LALDestroySFTCatalog ( status->statusPtr, &catalog ), status );
  }

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* LALCheckSFTs() */


/* checks the SFTs in a given SFTcatalog */
void
LALCheckSFTCatalog ( LALStatus *status,			/**< pointer to LALStatus structure */
		     INT4 *check_result,  /**< LAL-status of SFT-operations */
		     SFTCatalog *catalog  /**< catalog of SFTs to check */
		     )
{
  UINT4 i;

  INITSTATUS(status);

  ASSERT ( check_result, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( catalog,      status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );

  (*check_result) = 0;

  /* step through SFTs and check CRC64 */
  for ( i=0; i < catalog->length; i ++ )
    {
      FILE *fp;

      switch ( catalog->data[i].version  )
	{
	case 1:	/* version 1 had no CRC  */
	  continue;
	case 2:
	  if ( (fp = fopen_SFTLocator ( catalog->data[i].locator )) == NULL )
	    {
	      XLALPrintError ( "Failed to open locator '%s'\n",
			      XLALshowSFTLocator ( catalog->data[i].locator ) );
	      (*check_result) = SFTFILEIO_EFILE;
	      goto sft_failed;
	    }
	  if ( ! has_valid_v2_crc64 ( fp ) != 0 )
	    {
	      XLALPrintError ( "CRC64 checksum failure for SFT '%s'\n",
			      XLALshowSFTLocator ( catalog->data[i].locator ) );
	      (*check_result) = SFTFILEIO_ECRC64;
	      fclose(fp);
	      goto sft_failed;
	    }
	  fclose(fp);
	  break;

	default:
	  XLALPrintError ( "Illegal SFT-version encountered : %d\n", catalog->data[i].version );
	  (*check_result) = SFTFILEIO_EVERSION;
	  goto sft_failed;
	  break;
	} /* switch (version ) */

    } /* for i < numSFTs */

 sft_failed:

  RETURN ( status );

} /* LALCheckSFTCatalog() */

/** Load timestamps file into LIGOTimeGPSVector struct, allocated here.
 *
 * The timestamps file is of the format: <repeated lines of the form "seconds nano-seconds">
 * allowing for '%#' as comments, which are ignored.
 *
 */
LIGOTimeGPSVector *
XLALReadTimestampsFile ( const CHAR *fname )
{
  LIGOTimeGPSVector *timestamps = NULL;

  /** check input consistency */
  if ( !fname ) {
    XLALPrintError ( "%s: NULL input 'fname'", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /* read and parse timestamps-list file contents*/
  LALParsedDataFile *flines = NULL;
  if ( XLALParseDataFile ( &flines, fname ) != XLAL_SUCCESS )
    XLAL_ERROR_NULL ( XLAL_EFUNC );

  UINT4 numTS = flines->lines->nTokens;
  /* allocate and initialized segment list */
  if ( ( timestamps = XLALCreateTimestampVector ( numTS )) == NULL )
    XLAL_ERROR_NULL ( XLAL_EFUNC );

  UINT4 iTS;
  for ( iTS = 0; iTS < numTS; iTS ++ )
    {
      INT4 secs, ns;
      if ( sscanf ( flines->lines->tokens[iTS], "%d %d\n", &secs, &ns ) != 2 ) {
        XLALPrintError ("%s: failed to parse data-line %d: '%s' in timestamps-file '%s' (needs to be in format 'sec ns')\n", __func__, iTS + 1, flines->lines->tokens[iTS], fname );
        XLALDestroyTimestampVector ( timestamps );
        XLALDestroyParsedDataFile ( flines );
        XLAL_ERROR_NULL ( XLAL_ESYS );
      }
      if ( ( secs < 0 ) || ( ns < 0 ) ) {
	XLALPrintError ("%s: timestamps-file contains negative time-entry in line %d : s = %d, ns = %d\n", __func__, iTS, secs, ns );
	XLAL_ERROR_NULL ( XLAL_EDOM );
      }
      if ( ns > 999999999 ) {
        XLALPrintError ("%s: timestamps-file contains nano-seconds entry >= 1 billion ns in line %d: s = %d, ns = %d\n", __func__, iTS, secs, ns );
	XLAL_ERROR_NULL ( XLAL_EDOM );
      }

      timestamps->data[iTS].gpsSeconds = secs;
      timestamps->data[iTS].gpsNanoSeconds = ns;

    } /* for iTS < numTS */

  /* free parsed segment file contents */
  XLALDestroyParsedDataFile ( flines );

  return timestamps;

} /* XLALReadTimestampsFile() */



/**  [DEPRECATED] Read timestamps file and returns timestamps vector (alloc'ed in here!).
 *
 * \deprecated Use XLALReadTimestampsFile() instead.
 */
void
LALReadTimestampsFile (LALStatus* status, LIGOTimeGPSVector **timestamps, const CHAR *fname)
{
  FILE *fp;
  LIGOTimeGPSVector *ts = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (fname, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT (timestamps, status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT (*timestamps == NULL, status, SFTFILEIO_ENONULL,  SFTFILEIO_MSGENONULL);

  XLALPrintDeprecationWarning("LALReadTimestampsFile()", "XLALReadTimestampsFile()");

  if ( (fp = LALFopen( fname, "r")) == NULL) {
    XLALPrintError("\nUnable to open timestampsname file %s\n\n", fname);
    ABORT (status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE);
  }

  /* initialize empty timestamps-vector*/
  if ( (ts = LALCalloc(1, sizeof(LIGOTimeGPSVector))) == NULL) {
    ABORT (status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);
  }

  while(1)
    {
      INT4 secs, ns;
      if (fscanf ( fp, "%d  %d\n", &secs, &ns ) != 2)
	break;

      if ( ( secs < 0 ) || ( ns < 0 ) ) {
	XLALPrintError ("\nERROR: timestamps-file contained negative time-entry in line %d \n\n",
		       ts->length);
	ABORT ( status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
      }

      /* make space for the new entry */
      ts->length ++;
      if ( (ts->data = LALRealloc(ts->data, ts->length * sizeof(ts->data[0])) ) == NULL) {
	ABORT (status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);
      }

      ts->data[ts->length - 1].gpsSeconds = (UINT4) secs;
      ts->data[ts->length - 1].gpsNanoSeconds = (UINT4) ns;

    } /* while entries found */
  fclose(fp);

  /* hand over timestamps vector */
  (*timestamps) = ts;

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALReadTimestampsFile() */


/** Write the given *v2-normalized* (i.e. dt x DFT) SFTtype to a FILE pointer.
 *  Add the comment to SFT if SFTcomment != NULL.
 *
 * NOTE: Currently this only supports writing v2-SFTs.
 * If you need to write a v1-SFT, you should use LALWrite_v2SFT_to_v1file()
 *
 * NOTE2: the comment written into the SFT-file contains the 'sft->name' field concatenated with
 * the user-specified 'SFTcomment'
 *
 */
int
XLALWriteSFT2fp ( const SFTtype *sft,	/**< SFT to write to disk */
                  FILE *fp,		/**< pointer to open file */
                  const CHAR *SFTcomment)/**< optional comment (for v2 only) */
{
  UINT4 comment_len = 0;
  CHAR *_SFTcomment;
  UINT4 pad_len = 0;
  CHAR pad[] = {0, 0, 0, 0, 0, 0, 0};	/* for comment-padding */
  _SFT_header_v2_t rawheader;

  /* check input consistency */
  if (!sft || !sft->data || sft->deltaF <= 0 || sft->f0 < 0 || sft->data->length ==0 )
    XLAL_ERROR ( XLAL_EINVAL );
  if (!( (sft->epoch.gpsSeconds >= 0) && (sft->epoch.gpsNanoSeconds >= 0) ))
    XLAL_ERROR ( XLAL_EINVAL );
  if (!( sft->epoch.gpsNanoSeconds < 1000000000 ))
    XLAL_ERROR ( XLAL_EINVAL );
  if ( !is_valid_detector(sft->name) ) {
    XLALPrintError ("\nInvalid detector prefix '%c%c'\n\n", sft->name[0], sft->name[1] );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  if ( !fp )
    XLAL_ERROR ( XLAL_EINVAL );


  /* concat sft->name + SFTcomment for SFT-file comment-field */
  comment_len = strlen(sft->name) + 1;
  if ( SFTcomment )
    comment_len += strlen(SFTcomment) + 2;	/* separate by "; " */

  if ( (_SFTcomment = XLALCalloc( comment_len, sizeof(CHAR) )) == NULL ) {
    XLAL_ERROR( XLAL_ENOMEM );
  }
  strcpy ( _SFTcomment, sft->name );
  if ( SFTcomment ) {
    strcat ( _SFTcomment, "; " );
    strcat ( _SFTcomment, SFTcomment );
  }

  /* comment length including null terminator to string must be an
   * integer multiple of eight bytes.
   */
  pad_len = (8 - (comment_len % 8)) % 8;

  /* ----- fill out header */
  rawheader.version        		= 2;
  rawheader.gps_sec        		= sft->epoch.gpsSeconds;
  rawheader.gps_nsec       		= sft->epoch.gpsNanoSeconds;
  rawheader.tbase          		= 1.0 / sft->deltaF;
  rawheader.first_frequency_index 	= MYROUND( sft->f0 / sft->deltaF );
  rawheader.nsamples       		= sft->data->length;
  rawheader.crc64          		= 0;	/* set to 0 for crc-calculation */
  rawheader.detector[0]    		= sft->name[0];
  rawheader.detector[1]    		= sft->name[1];
  rawheader.padding[0]     		= 0;
  rawheader.padding[1]     		= 0;
  rawheader.comment_length 		= comment_len + pad_len;

  /* ----- compute CRC */
  rawheader.crc64 = calc_crc64((const CHAR*)&rawheader, sizeof(rawheader), ~(0ULL));

  rawheader.crc64 = calc_crc64((const CHAR*)_SFTcomment, comment_len, rawheader.crc64);
  rawheader.crc64 = calc_crc64((const CHAR*)pad, pad_len, rawheader.crc64);

  rawheader.crc64 = calc_crc64((const CHAR*) sft->data->data, sft->data->length * sizeof( *sft->data->data ), rawheader.crc64);

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

} /* XLALWriteSFT2fp() */

/** Write the given *v2-normalized* (i.e. dt x DFT) SFTtype to a v2-SFT file.
 *  Add the comment to SFT if SFTcomment != NULL.
 *
 * NOTE: Currently this only supports writing v2-SFTs.
 * If you need to write a v1-SFT, you should use LALWrite_v2SFT_to_v1file()
 *
 * NOTE2: the comment written into the SFT-file contains the 'sft->name' field concatenated with
 * the user-specified 'SFTcomment'
 *
 */
int
XLALWriteSFT2file(
		  const SFTtype *sft,		/**< SFT to write to disk */
		  const CHAR *fname,		/**< filename */
		  const CHAR *SFTcomment)	/**< optional comment (for v2 only) */
{
  FILE  *fp = NULL;

  /*   Make sure the arguments are not NULL */
  if (!( sft ))
    XLAL_ERROR ( XLAL_EINVAL );
  if (!( sft->data ))
    XLAL_ERROR ( XLAL_EINVAL );
  if (!( fname ))
    XLAL_ERROR ( XLAL_EINVAL );
 
  if ( !is_valid_detector(sft->name) ) {
    XLALPrintError ("\nInvalid detector prefix '%c%c'\n\n", sft->name[0], sft->name[1] );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  /* open SFT-file for writing */
  if ( (fp = LALFopen ( fname, "wb" )) == NULL )
    {
      XLALPrintError ("\nFailed to open file '%s' for writing: %s\n\n", fname, strerror(errno));
      XLAL_ERROR ( XLAL_EIO );
    }

  /* write SFT to file */
  if ( XLALWriteSFT2fp (sft, fp, SFTcomment) != XLAL_SUCCESS ) {
    XLAL_ERROR ( XLAL_EIO );
  }

  fclose(fp);

  return XLAL_SUCCESS;

} /* XLALWriteSFT2file() */

void
LALWriteSFT2file (LALStatus *status,			/**< pointer to LALStatus structure */
		  const SFTtype *sft,		/**< SFT to write to disk */
		  const CHAR *fname,		/**< filename */
		  const CHAR *SFTcomment)	/**< optional comment (for v2 only) */
{
  XLALPrintDeprecationWarning("LALWriteSFT2file", "XLALWriteSFT2file");
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  if ( XLALWriteSFT2file( sft, fname, SFTcomment ) != XLAL_SUCCESS ) {
    ABORT ( status, LAL_EXLAL, LAL_MSGEXLAL );
  }
  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* WriteSFTtoFile() */



/** Write the given *v2-normalized* (i.e. dt x DFT) SFTVector to a directory.
 *  Add the comment to SFT if SFTcomment != NULL.
 *
 * NOTE: Currently this only supports writing v2-SFTs.
 * If you need to write a v1-SFT, you should use LALWriteSFTfile()
 *
 * Output SFTs have naming convention following LIGO-T040164-01
 */
int
XLALWriteSFTVector2Dir(
		       const SFTVector *sftVect,	/**< SFT vector to write to disk */
		       const CHAR *dirname,		/**< base filename (including directory path)*/
		       const CHAR *SFTcomment,		/**< optional comment (for v2 only) */
		       const CHAR *description)         /**< optional sft description to go in the filename */
{
  UINT4 length, k;
  CHAR *filename = NULL;
  CHAR filenumber[16];
  SFTtype *sft;
  UINT4 timeBase, duration;
  UINT4 filenamelen;
  LIGOTimeGPS time0;

  if (! (sftVect) ) XLAL_ERROR ( XLAL_EINVAL );
  if (! (sftVect->data) ) XLAL_ERROR ( XLAL_EINVAL );
  if (! (sftVect->length > 0) ) XLAL_ERROR ( XLAL_EINVAL );
  if (! (dirname) ) XLAL_ERROR ( XLAL_EINVAL );

  length = sftVect->length;

  filenamelen = 128 + strlen(dirname);
  if ( description )
    filenamelen += strlen ( description );

  if ( (filename = (CHAR *)XLALCalloc(1, filenamelen )) == NULL) {
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  /* will not be same as actual sft timebase if it is not
     an integer number of seconds */
  timeBase = floor(1.0/sftVect->data[0].deltaF + 0.5);

  for ( k = 0; k < length; k++) {

    sft = sftVect->data + k;
    if ( sft == NULL ) {
      XLAL_ERROR ( XLAL_EFAULT );
    }

    if ( sft->name == NULL ) {
      XLAL_ERROR ( XLAL_EFAULT );
    }


    time0 = sft->epoch;

    /* calculate sft 'duration' -- may be different from timebase if nanosecond
       of sft-epoch is non-zero */
    duration = timeBase;
    if ( time0.gpsNanoSeconds > 0) {
      duration += 1;
    }

    /* create the k^th filename following naming convention
       -- try to simplify this*/
    strcpy( filename, dirname);
    strcat( filename, "/");
    strncat( filename, sft->name, 1);
    strcat( filename, "-1_"); /* single (not merged) sft */
    strncat( filename, sft->name, 2); /* full detector name */
    strcat( filename, "_");
    sprintf( filenumber, "%d", timeBase); /* sft timebase */
    strcat( filename, filenumber);
    strcat( filename, "SFT");
    if ( description ) {
      strcat( filename, "_");
      strcat( filename, description);
    }
    strcat( filename, "-");
    sprintf( filenumber, "%09d", sft->epoch.gpsSeconds);
    strncat( filename, filenumber, 9);
    strcat( filename, "-");
    sprintf( filenumber, "%d", duration);
    strcat( filename, filenumber);
    strcat( filename, ".sft");

    /* write the k^th sft */
    if ( XLALWriteSFT2file( sft, filename, SFTcomment ) != XLAL_SUCCESS ) {
      XLAL_ERROR ( xlalErrno );
    }
  }

  XLALFree(filename);

  return XLAL_SUCCESS;

} /* XLALWriteSFTVector2Dir() */

void
LALWriteSFTVector2Dir (LALStatus *status,			/**< pointer to LALStatus structure */
		       const SFTVector *sftVect,	/**< SFT vector to write to disk */
		       const CHAR *dirname,		/**< base filename (including directory path)*/
		       const CHAR *SFTcomment,		/**< optional comment (for v2 only) */
		       const CHAR *description)         /**< optional sft description to go in the filename */
{
  XLALPrintDeprecationWarning("LALWriteSFTVector2Dir", "XLALWriteSFTVector2Dir");
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  if ( XLALWriteSFTVector2Dir( sftVect, dirname, SFTcomment, description ) != XLAL_SUCCESS ) {
    ABORT ( status, LAL_EXLAL, LAL_MSGEXLAL );
  }
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



/** Write the given *v2-normalized* (i.e. dt x DFT) SFTVector to a single concatenated SFT file.
 *  Add the comment to SFT if SFTcomment != NULL.
 *
 * NOTE: Currently this only supports writing v2-SFTs.
 * If you need to write a v1-SFT, you should use LALWriteSFTfile()
 */
int
XLALWriteSFTVector2File(
		       const SFTVector *sftVect,	/**< SFT vector to write to disk */
		       const CHAR *filename,		/**< filename of concatenated SFT */
		       const CHAR *SFTcomment)		/**< optional comment (for v2 only) */
{
  UINT4 length, k;
  FILE *fp = NULL;
  SFTtype *sft;

  if (! (sftVect) ) XLAL_ERROR ( XLAL_EINVAL );
  if (! (sftVect->data) ) XLAL_ERROR ( XLAL_EINVAL );
  if (! (sftVect->length > 0) ) XLAL_ERROR ( XLAL_EINVAL );
  if (! (filename) ) XLAL_ERROR ( XLAL_EINVAL );

  length = sftVect->length;

  /* open SFT-file for writing */
  if ( (fp = LALFopen ( filename, "wb" )) == NULL )
    {
      XLALPrintError ("\nFailed to open file '%s' for writing: %s\n\n", filename, strerror(errno));
      XLAL_ERROR ( XLAL_EIO );
    }

  for ( k = 0; k < length; k++) {

    sft = sftVect->data + k;
    if ( sft == NULL ) {
      XLAL_ERROR ( XLAL_EFAULT );
    }

    if ( sft->name == NULL ) {
      XLAL_ERROR ( XLAL_EFAULT );
    }

    /* write the k^th sft */
    if ( XLALWriteSFT2fp ( sft, fp, SFTcomment ) != XLAL_SUCCESS ) {
      XLAL_ERROR ( xlalErrno );
    }
  }

  fclose(fp);

  return XLAL_SUCCESS;

} /* XLALWriteSFTVector2File() */



/** For backwards-compatibility: write a *v2-normalized* (ie dt x DFT) SFTtype
 * to a v1-SFT file.
 *
 * NOTE: the only difference to WriteSFTfile() is that the data-normalization
 *      is changed back to v1-type 'DFT', by dividing the dt corresponding to the
 *      frequency-band contained in the SFTtype.
 */
void
LALWrite_v2SFT_to_v1file (LALStatus *status,			/**< pointer to LALStatus structure */
			  const SFTtype *sft,		/**< SFT to write to disk */
			  const CHAR *fname)		/**< filename */
{
  UINT4 i, numBins;
  REAL8 Band, dt;
  SFTtype v1SFT;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL and perform basic checks*/
  ASSERT (sft,   status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);
  ASSERT (sft->data,  status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT (sft->deltaF > 0, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT (sft->f0 >= 0, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT ( (sft->epoch.gpsSeconds >= 0) && (sft->epoch.gpsNanoSeconds >= 0), status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT ( sft->epoch.gpsNanoSeconds < 1000000000, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT (sft->data->length > 0, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);

  ASSERT (fname, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);

  if ( !is_valid_detector(sft->name) ) {
    ABORT ( status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL );
  }

  numBins = sft->data->length;
  Band = sft->deltaF * numBins ;
  dt = 1.0 / (2.0 * Band);

  v1SFT.data = NULL;
  TRY ( LALCopySFT (status->statusPtr, &v1SFT, sft ), status );

  for ( i=0; i < numBins; i ++ )
    {
      v1SFT.data->data[i].re = (REAL4) ( (REAL8)v1SFT.data->data[i].re / dt );
      v1SFT.data->data[i].im = (REAL4) ( (REAL8)v1SFT.data->data[i].im / dt );
    }

  TRY ( LALWriteSFTfile (status->statusPtr, &v1SFT, fname ), status );

  XLALDestroyCOMPLEX8Vector ( v1SFT.data );

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* LALWrite_v2SFT_to_v1file() */




/** [OBSOLETE] Write a *v1-normalized* (i.e. raw DFT) SFTtype to a SFT-v1 file.
 *
 * \note:only SFT-spec v1.0 is supported, and the SFTtype must follow the
 * *obsolete* v1-normalization. => Use LALWriteSFT2file() to write v2 SFTs !
 *
 */
void
LALWriteSFTfile (LALStatus  *status,			/**< pointer to LALStatus structure */
		 const SFTtype *sft,		/**< SFT to write to disk */
		 const CHAR *outfname)		/**< filename */
{
  FILE  *fp = NULL;
  COMPLEX8  *inData;
  INT4  i;
  UINT4 datalen;
  REAL4  *rawdata;
  CHAR *rawheader, *ptr;
  SFTHeader header;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL and perform basic checks*/
  ASSERT (sft,   status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);
  ASSERT (sft->data,  status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT (sft->deltaF > 0, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT (outfname, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);

  /* fill in the header information */
  header.version = 1.0;
  header.gpsSeconds = sft->epoch.gpsSeconds;
  header.gpsNanoSeconds = sft->epoch.gpsNanoSeconds;
  header.timeBase = 1.0 / sft->deltaF;
  header.fminBinIndex = (INT4) floor (sft->f0 / sft->deltaF + 0.5);	/* round to closest int! */
  header.length = sft->data->length;

  /* build raw header for writing to disk */
  rawheader = LALCalloc (1, sizeof(_SFT_header_v1_t) );
  if (rawheader == NULL) {
    ABORT (status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);
  }
  ptr = rawheader;
  memcpy( ptr, &header.version, sizeof(REAL8) );
  ptr += sizeof (REAL8);
  memcpy( ptr, &header.gpsSeconds, sizeof(INT4) );
  ptr += sizeof (INT4);
  memcpy( ptr, &header.gpsNanoSeconds, sizeof(INT4) );
  ptr += sizeof (INT4);
  memcpy( ptr, &header.timeBase, sizeof(REAL8) );
  ptr += sizeof (REAL8);
  memcpy( ptr, &header.fminBinIndex, sizeof(INT4) );
  ptr += sizeof (INT4);
  memcpy( ptr, &header.length, sizeof(INT4) );

  /* write data into a contiguous REAL4-array */
  datalen = 2 * header.length * sizeof(REAL4);	/* amount of bytes for SFT-data */

  rawdata = LALCalloc (1, datalen);
  if (rawdata == NULL) {
    LALFree (rawheader);
    ABORT (status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);
  }

  inData = sft->data->data;
  for ( i = 0; i < header.length; i++)
    {
      rawdata[2 * i]     = inData[i].re;
      rawdata[2 * i + 1] = inData[i].im;
    } /* for i < length */


  /* open the file for writing */
  fp = LALFopen(outfname, "wb");
  if (fp == NULL) {
    LALFree (rawheader);
    LALFree (rawdata);
    XLALPrintError ("\nFailed to open file '%s' for writing!\n\n", outfname );
    ABORT (status, SFTFILEIO_EFILE,  SFTFILEIO_MSGEFILE);
  }

  /* write the header*/
  if( fwrite( rawheader, sizeof(_SFT_header_v1_t), 1, fp) != 1) {
    LALFree (rawheader);
    LALFree (rawdata);
    fclose (fp);
    ABORT (status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE);
  }

  /* write the data */
  if (fwrite( rawdata, datalen, 1, fp) != 1) {
    LALFree (rawheader);
    LALFree (rawdata);
    fclose (fp);
    ABORT (status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE);
  }

  /* done */
  fclose(fp);
  LALFree (rawheader);
  LALFree (rawdata);


  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* WriteSFTtoFile() */



/** Free an 'SFT-catalogue' */
void
LALDestroySFTCatalog ( LALStatus *status,			/**< pointer to LALStatus structure */
		       SFTCatalog **catalog )	/**< the 'catalogue' to free */
{
  INITSTATUS(status);

  ASSERT ( catalog, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );

  if ( *catalog )
    {
      if ( (*catalog) -> data )
	{
	  UINT4 i;
	  for ( i=0; i < (*catalog)->length; i ++ )
	    {
	      SFTDescriptor *ptr = &( (*catalog)->data[i] );
	      if ( ptr->locator )
		{
		  if ( ptr->locator->fname )
		    LALFree ( ptr->locator->fname );
		  LALFree ( ptr->locator );
		}
	      if ( ptr->comment )
		LALFree ( ptr->comment );

	      /* this should not happen, but just in case: free data-entry in SFT-header */
	      if ( ptr->header.data )
		XLALDestroyCOMPLEX8Sequence (ptr->header.data);
	    } /* for i < length */

	  LALFree ( (*catalog)->data );

	} /* if *catalog->data */

      LALFree ( *catalog );

    } /* if *catalog */

  (*catalog) = NULL;

  RETURN ( status );

} /* LALDestroySFTCatalog() */



/** Free an 'SFT-catalogue' */
void
XLALDestroySFTCatalog ( SFTCatalog *catalog  /**< the 'catalogue' to free */ )
{
  if ( catalog ) {

    if ( catalog -> data )
      {
        UINT4 i;
        for ( i=0; i < catalog->length; i ++ )
          {
            SFTDescriptor *ptr = &( catalog->data[i] );
            if ( ptr->locator )
              {
                if ( ptr->locator->fname )
                  XLALFree ( ptr->locator->fname );
                XLALFree ( ptr->locator );
              }
            if ( ptr->comment )
              XLALFree ( ptr->comment );

            /* this should not happen, but just in case: free data-entry in SFT-header */
            if ( ptr->header.data )
              XLALDestroyCOMPLEX8Sequence (ptr->header.data);
          } /* for i < length */

        catalog->length = 0;

        XLALFree ( catalog->data );

      } /* if catalog->data */

    XLALFree ( catalog );

  } /* if catalog */

} /* XLALDestroySFTCatalog() */


/** Mostly for *debugging* purposes: provide a user-API to allow inspecting the SFT-locator
 * [which is an OPAQUE entry in the SFTCatalog!]
 *
 * NOTE: this returns a STATIC string, so don't try to FREE it, and make a copy if you need
 * to keep it beyond one call of this function!
 *
 */
const CHAR *
XLALshowSFTLocator ( const struct tagSFTLocator *locator )
{
  static CHAR ret[512];

  if ( !locator )
    return NULL;

  snprintf ( ret, sizeof(ret), "%s : %ld", locator->fname, locator->offset );
  ret[ sizeof(ret) - 1 ] = 0;

  return ret;

} /* XLALshowSFTLocator() */


INT4 XLALCountIFOsInCatalog( const SFTCatalog *catalog)
{

  UINT4 k, j, numifo=0, length;
  CHAR  *name=NULL;
  CHAR  **ifolist=NULL; /* list of ifo names */

  length = catalog->length;

  name = (CHAR *)LALCalloc(3, sizeof(CHAR));

  ifolist = (CHAR **)LALCalloc( length, sizeof(CHAR *));
  for ( k = 0; k < length; k++)
    ifolist[k] = (CHAR *)LALCalloc( 3, sizeof(CHAR));

  /* go through catalog and look at each ifo name */
  for ( k = 0; k < length; k++)
    {
      strncpy( name, catalog->data[k].header.name, 3 );

      /* go through list of ifos till a match is found or list is exhausted */
      for ( j = 0; ( j < numifo ) && strncmp( name, ifolist[j], 3); j++ )
	;

      if ( j >= numifo )
	{
	  /* add ifo to list of ifos */
	  strncpy( ifolist[numifo], name, 3);
	  numifo++;
	}

    }

  LALFree(name);
  for ( j = 0; j < catalog->length; j++)
    LALFree(ifolist[j]);
  LALFree(ifolist);

  return numifo;

}

/*================================================================================
 * OBSOLETE and deprecated SFT-v1 API :
 *================================================================================*/

/** [DEPRECATED] Function to read and return a list of SFT-headers for given
 * filepattern and start/end times.
 *
 * The \a startTime and \a endTime entries are allowed to be NULL,
 * in which case they are ignored.
 *
 * \note We return the headers as an SFTVector, but with empty data-fields.
 */
void
LALGetSFTheaders (LALStatus *status,			/**< pointer to LALStatus structure */
		  SFTVector **headers,		/**< [out] Vector of SFT-headers */
		  const CHAR *fpattern,		/**< path/filepattern */
		  const LIGOTimeGPS *startTime,	/**< include only SFTs after this time (can be NULL) */
		  const LIGOTimeGPS *endTime)	/**< include only SFTs before this (can be NULL)*/
{
  UINT4 i, numSFTs;
  SFTVector *out = NULL;
  LALStringVector *fnames;
  REAL8 t0, t1;		/* start- and end-times as reals */
  UINT4 numHeaders;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* check input */
  ASSERT ( headers, status,  SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT ( *headers == NULL, status,  SFTFILEIO_ENONULL,  SFTFILEIO_MSGENONULL);
  ASSERT ( fpattern, status,  SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);

  if ( startTime )
    t0 = GPS2REAL8( *startTime );
  else
    t0 = 0;

  if ( endTime )
    t1 = GPS2REAL8 ( *endTime );
  else
    t1 = 0;

  /* get filelist of files matching fpattern */
  if ( (fnames = find_files (fpattern)) == NULL) {
    ABORT (status, SFTFILEIO_EGLOB, SFTFILEIO_MSGEGLOB);
  }

  /* prepare output-vector */
  if ( (out = LALCalloc (1, sizeof(*out) )) == NULL ) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }

  numSFTs = fnames->length;

  /* main loop: load all SFT-headers and put them into an SFTvector */
  numHeaders = 0;
  for (i=0; i < numSFTs; i++)
    {
      SFTHeader header;
      REAL8 epoch;
      SFTtype *sft;

      TRY ( LALReadSFTheader (status->statusPtr, &header, fnames->data[i]), status );
      epoch = 1.0 * header.gpsSeconds + 1.0e-9 * header.gpsNanoSeconds;

      if ( t0 && ( epoch < t0 ) )	/* epoch earlier than start-time ==> skip */
	continue;
      if ( t1 && ( epoch > t1 ) )	/* epoch later than end-time ==> skip */
	continue;

      /* found suitable SFT ==> add to SFTVector of headers */
      numHeaders ++;
      if ( ( out->data = LALRealloc(out->data, numHeaders * sizeof( *out->data) )) == NULL )
	{
	  LALFree (out);
	  ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
	}

      sft = &(out->data[numHeaders - 1]);

      /* fill in header-data */
      strncpy (sft->name, fnames->data[i], LALNameLength);
      sft->name[LALNameLength - 1 ] = '\0';	/* make sure it's 0-terminated */
      sft->deltaF  			= 1.0 / header.timeBase;
      sft->f0      			= header.fminBinIndex / header.timeBase;
      sft->epoch.gpsSeconds     	= header.gpsSeconds;
      sft->epoch.gpsNanoSeconds 	= header.gpsNanoSeconds;

      sft->data = NULL;	/* no SFT-data proper */

    } /* for i < numSFTs */

  out->length = numHeaders;

  XLALDestroyStringVector (fnames);

  *headers = out;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALGetSFTheaders() */



/** [DEPRECATED] Basic SFT reading-function.
 *  Given a filename \a fname and frequency-limits [\a fMin, \a fMax],
 *  returns an SFTtype \a sft containing the SFT-data.
 *
 * \note 1) the actual returned frequency-band is
 * <tt>[floor(Tsft * fMin), ceil(Tsft * fMax)] / Tsft</tt>,
 * i.e. we need to round to an integer frequency-bin within the SFT, but
 * the requested frequency-band is guaranteed to be contained in the output
 * (if present in the SFT-file), but can be slightly larger.
 *
 * 2) The special input <tt>fMin=fMax=0</tt> means to read and
 * return the <em>whole</em> frequency-band contained in the SFT-file.
 *
 * 3) Currently only SFTv1 are supported!!
 *
 */
void
LALReadSFTfile (LALStatus *status,			/**< pointer to LALStatus structure */
		SFTtype **sft, 		/**< [out] output SFT */
		REAL8 fMin, 		/**< lower frequency-limit */
		REAL8 fMax,		/**< upper frequency-limit */
		const CHAR *fname)	/**< path+filename */
{
  SFTHeader  header;		/* SFT file-header version1 */
  UINT4 readlen;
  INT4 fminBinIndex, fmaxBinIndex;
  SFTtype *outputSFT = NULL;
  UINT4 i;
  REAL4 renorm;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (sft, status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT (*sft == NULL, status, SFTFILEIO_ENONULL, SFTFILEIO_MSGENONULL);
  ASSERT (fname,  status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT (fMin <= fMax, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);

  /* read the header */
  TRY ( LALReadSFTheader (status->statusPtr, &header, fname), status);

  /* ----- figure out which data we want to read ----- */

  /* special case: fMin==fMax==0 means "read all" */
  if ( (fMin == 0) && (fMax == 0) )
    {
      fminBinIndex = header.fminBinIndex;
      fmaxBinIndex = fminBinIndex + header.length - 1;
    }
  else
    {
      /* find the right frequency-bin and number of bins
       * The rounding here is chosen such that the required
       * frequency-interval is _guaranteed_ to lie within the
       * returned range  */
      fminBinIndex = (INT4) floor (fMin * header.timeBase);  /* round this down */
      fmaxBinIndex = (INT4) ceil  (fMax * header.timeBase);  /* round up */
    }

  readlen = (UINT4)(fmaxBinIndex - fminBinIndex) + 1;	/* number of bins to read */

  /* allocate the final SFT to be returned */
  TRY ( LALCreateSFTtype (status->statusPtr, &outputSFT, readlen), status);


  /* and read it, using the lower-level function: */
  LALReadSFTdata (status->statusPtr, outputSFT, fname, fminBinIndex);
  BEGINFAIL (status) {
    LALDestroySFTtype (status->statusPtr, &outputSFT);
  } ENDFAIL (status);


  /*
   * NOTE: the following renormalization is necessary for v1-SFTs
   * as the data are not multiplied by dt, therefore extracting a
   * sub-band B' of the total band B requires mulitplication of
   * the data by B'/B.
   * SFTv2 will store FFT-data multiplied by dt and then this will
   * not be necessary any more.
   *
   */
  renorm = 1.0 * readlen / header.length;

  /* let's re-normalize and fill data into output-vector */
  if (renorm != 1)
    for (i=0; i < readlen; i++)
      {
	outputSFT->data->data[i].re *= renorm;
	outputSFT->data->data[i].im *= renorm;
      }

  /* that's it: return */
  *sft = outputSFT;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALReadSFTfile() */



/** [DEPRECATED] Higher-level SFT-reading function to read a whole vector of SFT files
 * and return an SFTvector \a sftvect. The handling of
 * [fMin, fMax] is identical to LALReadSFTfile().
 *
 * \note 1) the file-pattern \a fpattern can use a wide range of glob-patterns.
 *
 * 2) currently the SFTs matching the pattern are required to have the same
 * number of frequency bins, otherwise an error will be returned.
 * (This might be relaxed in the future).
 *
 * 3) This function does <em>not</em> use <tt>glob()</tt> and should therefore
 * be safe even under condor.
 *
 */
void
LALReadSFTfiles (LALStatus *status,			/**< pointer to LALStatus structure */
		 SFTVector **sftvect,	/**< [out] output SFT vector */
		 REAL8 fMin,	       	/**< lower frequency-limit */
		 REAL8 fMax,		/**< upper frequency-limit */
		 UINT4 wingBins,	/**< number of frequency-bins to be added left and right. */
		 const CHAR *fpattern)	/**< path/filepattern */
{

  UINT4 i, numSFTs;
  SFTVector *out = NULL;
  SFTtype *oneSFT = NULL;
  SFTHeader header;
  REAL8 dFreq; 		/* frequency spacing in SFT */
  REAL8 fWing;		/* frequency band to be added as "wings" to the 'physical band' */
  LALStringVector *fnames;
  UINT4 firstlen = 0;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (sftvect, status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT (*sftvect == NULL, status, SFTFILEIO_ENONULL, SFTFILEIO_MSGENONULL);
  ASSERT (fpattern,  status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT (fMin <= fMax, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);

  /* make filelist
   * NOTE: we don't use glob() as it was reported to fail under condor */
  if ( (fnames = find_files (fpattern)) == NULL) {
    ABORT (status, SFTFILEIO_EGLOB, SFTFILEIO_MSGEGLOB);
  }

  /* allocate head of output SFT-vector */
  if ( (out = LALCalloc ( 1, sizeof(SFTVector) )) == NULL ) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }

  /* read header of first sft to determine Tsft, and therefore dfreq */
  LALReadSFTheader(status->statusPtr, &header, fnames->data[0]);
  BEGINFAIL(status) {
    XLALDestroyStringVector (fnames);
  } ENDFAIL(status);

  numSFTs = fnames->length;
  dFreq = 1.0 / header.timeBase;
  fWing = wingBins * dFreq;

  /* main loop: load all SFTs and put them into the SFTvector */
  for (i=0; i < numSFTs; i++)
    {
      LALReadSFTfile (status->statusPtr, &oneSFT, fMin-fWing, fMax+fWing, fnames->data[i]);
      BEGINFAIL (status) {
	XLALDestroyStringVector (fnames);
	LALDestroySFTtype (status->statusPtr, &oneSFT);
	if (out) LALDestroySFTVector (status->statusPtr, &out);
      } ENDFAIL (status);

      if ( !firstlen )
	firstlen = oneSFT->data->length;
      /* make sure all SFTs have same length */
      if ( oneSFT->data->length != firstlen )
	{
	  XLALDestroyStringVector (fnames);
	  LALDestroySFTtype (status->statusPtr, &oneSFT);
	  LALDestroySFTVector (status->statusPtr, &out);
	  ABORT (status, SFTFILEIO_EDIFFLENGTH, SFTFILEIO_MSGEDIFFLENGTH);
	} /* if length(thisSFT) != common length */

      LALAppendSFT2Vector ( status->statusPtr, out, oneSFT );
      BEGINFAIL(status) {
	XLALDestroyStringVector (fnames);
	LALDestroySFTtype (status->statusPtr, &oneSFT);
	LALDestroySFTVector (status->statusPtr, &out);
      } ENDFAIL(status);

      LALDestroySFTtype (status->statusPtr, &oneSFT);
      oneSFT = NULL;	/* important for next call of LALReadSFTfile()! */

    } /* for i < numSFTs */

  XLALDestroyStringVector (fnames);

  *sftvect = out;

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALReadSFTfiles () */



/*================================================================================
 * LOW-level internal SFT-handling functions, should *NOT* be used outside this file!
 *================================================================================*/

/** Open an "SFT" defined by the SFT-locator, return a FILE-pointer to the beginning of this SFT.
 *
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
  if ( (fp = LALFopen( fname, "rb" )) == NULL )
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


/** Check the v2 SFT-block starting at fp for valid crc64 checksum.
 *  Restores filepointer before leaving.
 */
BOOLEAN
has_valid_v2_crc64 ( FILE *fp )
{
  long save_filepos;
  UINT8 computed_crc, ref_crc;
  SFTtype header;
  UINT4 numBins;
  CHAR *SFTcomment = NULL;
  UINT4 data_len;
  char block[BLOCKSIZE];
  UINT4 version;
  BOOLEAN need_swap;

  /* input consistency */
  if ( !fp )
    {
      XLALPrintError ("\nhas_valid_v2_crc64() was called with NULL filepointer!\n\n");
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

  if ( version != 2 )
    {
      XLALPrintError ("\nhas_valid_v2_crc64() was called on non-v2 SFT.\n\n");
      return -1;
    }

  /* ----- compute CRC ----- */
  /* read the header, unswapped, only to obtain it's crc64 checksum */
  if ( read_v2_header_from_fp ( fp, &header, &numBins, &computed_crc, &ref_crc, &SFTcomment, need_swap ) != 0 )
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
      computed_crc = calc_crc64( (const CHAR*)block, toread, computed_crc );

    } /* while data */

  /* check that checksum is consistent */
  return ( computed_crc == ref_crc );

} /* has_valid_v2_crc64 */


/** [DEPRECATED]: Low-level function to read only the SFT-header of a given file.
 *
 * NOTE: don't use! This function is obsolete and SFT-v1-specific and only kept for
 * backwards-compatibility with Hough-codes.
 */
void
LALReadSFTheader (LALStatus  *status,			/**< pointer to LALStatus structure */
		  SFTHeader   *header,	/**< [out] returned header */
		  const CHAR  *fname)	/**< path+filename */
{
  FILE *fp = NULL;
  SFTHeader  header1;
  CHAR *rawheader = NULL;
  CHAR *ptr = NULL;
  CHAR inVersion[8];
  REAL8 version;
  BOOLEAN swapEndian = 0;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (header, status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT (fname,  status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);

  /* opening the SFT binary file */
  fp = LALOpenDataFile( fname );
  if (fp == NULL) {
    ABORT (status, SFTFILEIO_EFILE,  SFTFILEIO_MSGEFILE);
  }

  /* read version-number */
  if  (fread (inVersion, sizeof(inVersion), 1, fp) != 1) {
    fclose (fp);
    if (lalDebugLevel) XLALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* try version 1.0 */
  version = 1.0;
  /* try swapping the version if it is not equal */
  if(memcmp(inVersion,&version,sizeof(version))){
    endian_swap (inVersion, sizeof(inVersion),1);
    swapEndian = 1;
    /* fail if still not true */
    if(memcmp(inVersion,&version,sizeof(version))){
      fclose (fp);
      if (lalDebugLevel) XLALPrintError ("\nOnly v1-SFTs supported at the moment!: %s\n\n", fname);
      ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
    }
  }

  /* read the whole header */
  rawheader = LALCalloc (1, sizeof(_SFT_header_v1_t) );
  if (rawheader == NULL) {
    fclose (fp);
    ABORT (status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);
  }

  rewind (fp);	/* go back to start */
  if (fread( rawheader, sizeof(_SFT_header_v1_t), 1, fp) != 1) {
    fclose (fp);
    LALFree (rawheader);
    if (lalDebugLevel) XLALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  fclose(fp);

  /* now fill-in the header-struct with the appropriate fields */
  /* NOTE: we have to do it this way, because the struct can have
   * padding in memory, so the fields are not guaranteed to lie 'close'
   * Endian-swapping ist done here if necessary
   */
  ptr = rawheader;
  if (swapEndian) endian_swap((CHAR*)ptr,sizeof(REAL8),1);
  memcpy( &header1.version, ptr, sizeof(REAL8) );
  ptr += sizeof(REAL8);
  if (swapEndian) endian_swap((CHAR*)ptr,sizeof(INT4),1);
  memcpy( &header1.gpsSeconds, ptr, sizeof(INT4) );
  ptr += sizeof(INT4);
  if (swapEndian) endian_swap((CHAR*)ptr,sizeof(INT4),1);
  memcpy( &header1.gpsNanoSeconds, ptr, sizeof(INT4) );
  ptr += sizeof(INT4);
  if (swapEndian) endian_swap((CHAR*)ptr,sizeof(REAL8),1);
  memcpy( &header1.timeBase, ptr, sizeof(REAL8) );
  ptr += sizeof(REAL8);
  if (swapEndian) endian_swap((CHAR*)ptr,sizeof(INT4),1);
  memcpy( &header1.fminBinIndex, ptr, sizeof(INT4) );
  ptr += sizeof(INT4);
  if (swapEndian) endian_swap((CHAR*)ptr,sizeof(INT4),1);
  memcpy( &header1.length, ptr, sizeof(INT4) );

  LALFree (rawheader);

  /* ----- do some consistency-checks on the header-fields: ----- */

  /* gps_sec and gps_nsec >= 0 */
  if ( (header1.gpsSeconds < 0) || (header1.gpsNanoSeconds <0) ) {
    if (lalDebugLevel) XLALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* tbase > 0 */
  if ( header1.timeBase <= 0 ) {
    if (lalDebugLevel) XLALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* fminindex >= 0 */
  if (header1.fminBinIndex < 0) {
    if (lalDebugLevel) XLALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* nsamples >= 0 */
  if (header1.length < 0) {
    if (lalDebugLevel) XLALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* ok, the SFT-header seems consistent, so let's return it */
  *header = header1;

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALReadSFTheader() */



/** [DEPRECATED] This is a function for low-level SFT data-reading:
 * the SFT-data is read starting from fminBinIndex and filled
 * into the pre-allocate vector sft of length N
 *
 * NOTE: !! NO re-normalization is done here!! this remains up
 *       to the caller of this function!!
 *
 */
void
LALReadSFTdata(LALStatus *status,			/**< pointer to LALStatus structure */
	       SFTtype    *sft,    /**< [out] output-SFT: assuming memory is allocated  */
	       const CHAR *fname,  /**< path+filename */
	       INT4 fminBinIndex)  /**< minimun frequency-index to read */
{
  FILE        *fp = NULL;
  SFTHeader  header;
  UINT4 offset, readlen;
  REAL4 *rawdata = NULL;
  CHAR inVersion[8];
  REAL8 version;
  BOOLEAN swapEndian = 0;
  UINT4 i;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (sft,   status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);
  ASSERT (sft->data, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);
  ASSERT (fname, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);

  /* Read header */
  TRY ( LALReadSFTheader (status->statusPtr, &header, fname), status);

  /* check that the required frequency-interval is part of the SFT */
  readlen = sft->data->length;
  if ( (fminBinIndex < header.fminBinIndex)
       || (fminBinIndex + (INT4)readlen > header.fminBinIndex + header.length) ) {
    ABORT (status, SFTFILEIO_EFREQBAND, SFTFILEIO_MSGEFREQBAND);
  }

  /* how many frequency-bins to skip */
  offset = fminBinIndex - header.fminBinIndex;

  /* open file for reading */
  if ( (fp = LALOpenDataFile( fname )) == NULL) {
    ABORT (status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE);
  }

  /* read version-number */
  if  (fread (inVersion, sizeof(inVersion), 1, fp) != 1) {
    fclose (fp);
    if (lalDebugLevel) XLALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* set invalid version */
  version = -1.0;

  if (version < 0) {
    /* try version 1.0 */
    version = 1.0;
    /* try swapping the version if it is not equal */
    if(memcmp(inVersion,&version,sizeof(version))){
      endian_swap (inVersion, sizeof(inVersion),1);
      swapEndian = 1;
      /* set invalid version if still not true */
      if(memcmp(inVersion,&version,sizeof(version))){
	version = -1;
	endian_swap (inVersion, sizeof(inVersion),1);
	swapEndian = 0;
      }
    }
  }

  /* fail if the version is invalid */
  if (version < 0) {
    fclose (fp);
    if (lalDebugLevel) XLALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* check compatibility of version with this function */
  if (version != 1) {
    fclose (fp);
    ABORT (status, SFTFILEIO_EVERSION, SFTFILEIO_MSGEVERSION);
  }

  /* skip SFT-header in file */
  rewind (fp);
  if (fseek(fp, sizeof(_SFT_header_v1_t), SEEK_SET) != 0) {
    fclose (fp);
    ABORT (status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE);
  }

  /* skip offset data points to the correct frequency-bin */
  if (fseek(fp, offset * 2 * sizeof(REAL4), SEEK_CUR) != 0) {
    fclose (fp);
    ABORT (status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE);
  }

  /* ----- prepare memory for data-reading ----- */
  rawdata = LALCalloc (1, 2 * readlen *sizeof(REAL4) );
  if (rawdata == NULL) {
    fclose (fp);
    ABORT (status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);
  }

  /* we don't rely on memory-packing, so we read into a REAL4 array first */
  if (fread( rawdata, 2 * readlen * sizeof(REAL4), 1, fp) != 1) {
    fclose(fp);
    LALFree (rawdata);
    ABORT (status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE);
  }
  fclose(fp);

  /* now fill data into output-vector */
  for (i=0; i < readlen; i++)
    {
      if (swapEndian)
	endian_swap((CHAR*)&rawdata[2*i], sizeof(REAL4), 2);
      sft->data->data[i].re = rawdata[2 * i];
      sft->data->data[i].im = rawdata[2 * i + 1];
    }

  LALFree (rawdata);

  /* now fill in the header-info */
  strncpy (sft->name, fname, LALNameLength);
  sft->name[LALNameLength - 1 ] = '\0';	/* make sure it's 0-terminated */
  sft->deltaF  			= 1.0 / header.timeBase;
  sft->f0      			= fminBinIndex / header.timeBase;
  sft->epoch.gpsSeconds     	= header.gpsSeconds;
  sft->epoch.gpsNanoSeconds 	= header.gpsNanoSeconds;


  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALReadSFTdata() */



/***********************************************************************
 * internal helper functions
 ***********************************************************************/

static BOOLEAN
timestamp_in_list( LIGOTimeGPS timestamp, LIGOTimeGPSVector *list )
{
  UINT4 i;
  LIGOTimeGPS *el;

  if ( !list )
    return FALSE;

  el = &(list->data[0]);
  for ( i=0; i < list->length; i++, el++)
    {
      if ( (timestamp.gpsSeconds == el->gpsSeconds) && ( timestamp.gpsNanoSeconds == el->gpsNanoSeconds ) )
	return TRUE;
    } /* for i < length */

  return FALSE;

} /* timestamp_in_list() */


/* check consistency constraints for SFT-blocks within a merged SFT-file,
 * see SFT-v2 spec */
static BOOLEAN
consistent_mSFT_header ( SFTtype header1, UINT4 version1, UINT4 nsamples1, SFTtype header2, UINT4 version2, UINT4 nsamples2 )
{
  /* 1) identical detector */
  if ( (header1.name[0] != header2.name[0]) || (header1.name[1] != header2.name[1]) )
    {
      XLALPrintError ("\nInvalid merged SFT: non-identical detectors\n\n");
      return FALSE;
    }

  /* 2) identical version-number */
  if ( version1 != version2 )
    {
      XLALPrintError ("\nInvalid merged SFT: non-identical version-numbers\n\n");
      return FALSE;
    }

  /* 3) increasing GPS-times */
  if ( GPS2REAL8 ( header1.epoch ) >= GPS2REAL8 ( header2.epoch ) )
    {
      XLALPrintError ("\nInvalid merged SFT: non-increasing GPS epochs \n\n" );
      return FALSE;
    }

  /* 4) identical tbase */
  if ( header1.deltaF != header2.deltaF )
    {
      XLALPrintError ("\nInvalid merged SFT: non-identical time baselines\n\n");
      return FALSE;
    }

  /* 5) identical start-frequency */
  if ( header1.f0 != header2.f0 )
    {
      XLALPrintError ("\nInvalid merged SFT: non-identical start-frequencies\n\n");
      return FALSE;
    }

  /* 6) identical number of frequency-bins */
  if ( nsamples1 != nsamples2 )
    {
      XLALPrintError ("\nInvalid merged SFT: non-identical number of frequency-bins\n\n" );
      return FALSE;
    }

  return TRUE;

} /* consistent_mSFT_header() */


/** Read bins from an SFT, leave filepointer at the end of the read SFT if successful,
 * leave fp at initial position if failure
 *
 * This uses read_sft_header_from_fp() for reading the header, and then reads the data.
 *
 * firstBin2read is the first bin to read from the SFT the fp points to,
 * lastBin2read is the last bin.
 * An error is issued if firstBin2read is not contained in the sft,
 * but NO ERROR is issued if lastBin2read is not contained.
 *
 * binread is set to the number of bins actually read (from firtsBin to either
 * lastBin or the end of the SFT
 *
 * NOTE: we DO NOT check the crc64 checksum in here, a separate "check-SFT" function
 * should be used for that.
 *
 * NOTE2:  The returned SFT is always normalized correctly according to SFT-v2 spec
 * (see LIGO-T040164-01-Z, and LIGO-T010095-00), irrespective of the input-files.
*/
static void
lal_read_sft_bins_from_fp ( LALStatus *status, SFTtype **sft, UINT4 *binsread, UINT4 firstBin2read, UINT4 lastBin2read , FILE *fp )
{
  SFTtype *ret = NULL;
  UINT4 version;
  UINT8 crc64;
  BOOLEAN swapEndian;
  UINT4 numBins2read;
  UINT4 firstSFTbin, lastSFTbin, numSFTbins;
  INT4 offsetBins;
  long offsetBytes;
  volatile REAL8 tmp;	/* intermediate results: try to force IEEE-arithmetic */

  INITSTATUS(status);
  ATTATCHSTATUSPTR ( status );

  ASSERT ( sft, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( *sft == NULL, status, SFTFILEIO_ENONULL, SFTFILEIO_MSGENONULL );

  TRY ( LALCreateSFTtype ( status->statusPtr, &ret, 0 ), status );

  if ( read_sft_header_from_fp (fp, ret, &version, &crc64, &swapEndian, NULL, &numSFTbins ) != 0 )
    {
      XLALPrintError ("\nFailed to read SFT-header!\n\n");
      LALDestroySFTtype ( status->statusPtr, &ret );
      ABORT ( status, SFTFILEIO_EHEADER, SFTFILEIO_MSGEHEADER );
    }

  tmp = ret->f0 / ret->deltaF;
  firstSFTbin = MYROUND ( tmp );
  lastSFTbin = firstSFTbin + numSFTbins - 1;

  if ( firstBin2read > lastBin2read )
    {
      if ( lalDebugLevel )
	XLALPrintError ("\nEmpty frequency-interval requested [%d, %d] bins\n\n",
		       firstBin2read, lastBin2read );
      ABORT ( status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL );
    }

  /* check that requested interval is found in SFT */
  if ( firstBin2read < firstSFTbin )
    {
      if ( lalDebugLevel )
	XLALPrintError ( "\nRequested bin is not contained in SFT! (%d < %d)\n\n",
			firstBin2read, firstSFTbin );
      LALDestroySFTtype ( status->statusPtr, &ret );
      ABORT ( status, SFTFILEIO_EFREQBAND, SFTFILEIO_MSGEFREQBAND );
    }
  if ( lastBin2read > lastSFTbin )
    lastBin2read = lastSFTbin;
  *binsread = lastBin2read - firstBin2read + 1;

  offsetBins = firstBin2read - firstSFTbin;
  offsetBytes = offsetBins * 2 * sizeof( REAL4 );
  numBins2read = lastBin2read - firstBin2read + 1;

  if ( fseek ( fp, offsetBytes, SEEK_CUR ) != 0 )
    {
      if ( lalDebugLevel )
	XLALPrintError ( "\nFailed to fseek() to first frequency-bin %d: %s\n\n",
			firstBin2read, strerror(errno) );
      LALDestroySFTtype ( status->statusPtr, &ret );
      ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
    }

  if ( (ret->data = XLALCreateCOMPLEX8Vector ( numBins2read )) == NULL ) {
    LALDestroySFTtype ( status->statusPtr, &ret );
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }

  if ( numBins2read != fread ( ret->data->data, 2*sizeof( REAL4 ), numBins2read, fp ) )
    {
      if (lalDebugLevel) XLALPrintError ("\nFailed to read %d bins from SFT!\n\n", numBins2read );
      LALDestroySFTtype ( status->statusPtr, &ret );
      ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
    }

  /* update the start-frequency entry in the SFT-header to the new value */
  ret->f0 = 1.0 * firstBin2read * ret->deltaF;

  /* take care of normalization and endian-swapping */
  if ( version == 1 || swapEndian )
    {
      UINT4 i;
      REAL8 band = 1.0 * numSFTbins * ret->deltaF;/* need the TOTAL frequency-band in the SFT-file! */
      REAL8 fsamp = 2.0 * band;
      REAL8 dt = 1.0 / fsamp;

      for ( i=0; i < numBins2read; i ++ )
	{
	  REAL4 *rep, *imp;

	  rep = &(ret->data->data[i].re);
	  imp = &(ret->data->data[i].im);

	  if ( swapEndian )
	    {
	      endian_swap( (CHAR *) rep, sizeof ( *rep ), 1 );
	      endian_swap( (CHAR *) imp, sizeof ( *imp ), 1 );
	    }

	  /* if the SFT-file was in v1-Format: need to renormalize the data now by 'Delta t'
	   * in order to follow the correct SFT-normalization
	   * (see LIGO-T040164-01-Z, and LIGO-T010095-00)
	   */
	  if ( version == 1 )
	    {
	      (*rep) *= dt;
	      (*imp) *= dt;
	    }
	} /* for i < numBins2read */
    } /* if SFT-v1 */

  /* return resulting SFT */
  (*sft) = ret;

  DETATCHSTATUSPTR ( status );
  RETURN (status);
} /* lal_read_sft_bins_from_fp() */



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
static int
read_sft_header_from_fp (FILE *fp, SFTtype *header, UINT4 *version, UINT8 *crc64, BOOLEAN *swapEndian, CHAR **SFTcomment, UINT4 *numBins )
{
  SFTtype head = empty_SFTtype;
  UINT4 nsamples;
  CHAR *comm = NULL;
  UINT8 ref_crc = 0;
  UINT8 header_crc;

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


  /* read this SFT-header with version-specific code */
  head = empty_SFTtype;

  switch( ver )
    {
    case 1:
      if ( read_v1_header_from_fp ( fp, &head, &nsamples, need_swap ) != 0 )
	goto failed;
      break;

    case 2:
      if ( read_v2_header_from_fp ( fp, &head, &nsamples, &header_crc, &ref_crc, &comm, need_swap ) != 0 )
	goto failed;
      break;

    default:
      XLALPrintError ("\nUnsupported SFT-version %d.\n\n", ver);
      goto failed;
      break;
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



/* ----- SFT v2 -specific header-reading function:
 *
 * return general SFTtype header, place filepointer at the end of the header if it succeeds,
 * set fp to initial position if it fails.
 * RETURN: 0 = OK, -1 = ERROR
 *
 * NOTE: fatal errors will produce a XLALPrintError() error-message, but
 * non-fatal 'SFT-format'-errors will only output error-messages if lalDebugLevel > 0.
 *
 */
static int
read_v2_header_from_fp ( FILE *fp, SFTtype *header, UINT4 *nsamples, UINT8 *header_crc64, UINT8 *ref_crc64, CHAR **SFTcomment, BOOLEAN swapEndian)
{
  _SFT_header_v2_t rawheader;
  long save_filepos;
  CHAR *comm = NULL;
  UINT8 crc;


  /* check input-consistency */
  if ( !fp || !header || !nsamples || !SFTcomment )
    {
      XLALPrintError ( "\nERROR read_v2_header_from_fp(): called with NULL input!\n\n");
      return -1;
    }
  if ( SFTcomment && (*SFTcomment != NULL) )
    {
      XLALPrintError ("\nERROR: Comment-string passed to read_v2_header_from_fp() is not NULL!\n\n");
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
      if (lalDebugLevel) XLALPrintError ("\nCould not read v2-header. %s\n\n", strerror(errno) );
      goto failed;
    }

  /* ----- compute CRC for the header:
   * NOTE: the CRC checksum is computed on the *bytes*, not the numbers,
   * so this must be computed before any endian-swapping.
   */
  {
    UINT8 save_crc = rawheader.crc64;
    rawheader.crc64 = 0;

    crc = calc_crc64((const CHAR*)&rawheader, sizeof(rawheader), ~(0ULL));

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

      /* v2-specific */
      endian_swap((CHAR*)(&rawheader.crc64),			sizeof(rawheader.crc64)			, 1);
      endian_swap((CHAR*)(&rawheader.comment_length),		sizeof(rawheader.comment_length)	, 1);
      /* ----- */

    } /* if endian_swap */

  /* double-check version-number */
  if ( rawheader.version != 2 )
    {
      XLALPrintError ("\nWrong SFT-version %d in read_v2_header_from_fp()\n\n", rawheader.version );
      goto failed;
    }

  if ( rawheader.nsamples <= 0 )
    {
      XLALPrintError ("\nNon-positive number of samples in SFT!\n\n");
      goto failed;
    }

  /* ----- v2-specific consistency-checks ----- */
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

  if ( ! is_valid_detector ( rawheader.detector ) )
    {
      XLALPrintError ("\nIllegal detector-name in SFT: '%c%c'\n\n",
					  rawheader.detector[0], rawheader.detector[1] );
      goto failed;
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
	  XLALPrintError ("\nCould not read %d-bytes comment\n\n");
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

	  crc = calc_crc64((const CHAR*)comm, comment_len, crc);
	  crc = calc_crc64((const CHAR*)pad, pad_len, crc);
	}

    } /* if comment_length > 0 */

  /*  ok: */
  memset ( header, 0, sizeof( *header ) );

  header->name[0]		= rawheader.detector[0];
  header->name[1]		= rawheader.detector[1];
  header->name[2]		= 0;

  header->epoch.gpsSeconds 	= rawheader.gps_sec;
  header->epoch.gpsNanoSeconds 	= rawheader.gps_nsec;

  header->f0 			= rawheader.first_frequency_index / rawheader.tbase;
  header->deltaF 		= 1.0 / rawheader.tbase;

  (*nsamples) = rawheader.nsamples;
  (*ref_crc64) = rawheader.crc64;
  (*SFTcomment) = comm;
  (*header_crc64) = crc;


  return 0;

 failed:
  /* restore filepointer initial position  */
  if ( fseek ( fp, save_filepos, SEEK_SET ) == -1 )
    XLALPrintError ("\nfseek() failed to return to intial fileposition: %s\n\n", strerror(errno) );

  /* free comment  if we allocated one */
  if ( comm )
    LALFree (comm);

  return -1;

} /* read_v2_header_from_fp() */


/* ----- SFT v1 -specific header-reading function:
 *
 * return general SFTtype header, place filepointer at the end of the header if it succeeds,
 * set fp to initial position if it fails.
 * RETURN: 0 = OK, -1 = ERROR
 *
 * NOTE: fatal errors will produce a XLALPrintError() error-message, but
 * non-fatal 'SFT-format'-errors will only output error-messages if lalDebugLevel > 0.
 *
 */
static int
read_v1_header_from_fp ( FILE *fp, SFTtype *header, UINT4 *nsamples, BOOLEAN swapEndian)
{
  _SFT_header_v1_t rawheader;
  long save_filepos;

  if ( !fp || !header || !nsamples)
    {
      XLALPrintError ( "\nERROR read_v1_header_from_fp(): called with NULL input!\n\n");
      return -1;
    }

  /* store fileposition for restoring in case of failure */
  if ( ( save_filepos = ftell(fp) ) == -1 )
    {
      XLALPrintError ("\nftell() failed: %s\n\n", strerror(errno) );
      return -1;
    }

  /* read the whole header */
  if (fread( &rawheader, sizeof(rawheader), 1, fp) != 1)
    {
      if (lalDebugLevel) XLALPrintError ("\nCould not read v1-header. %s\n\n", strerror(errno) );
      goto failed;
    }

  if (swapEndian)
    {
      endian_swap((CHAR*)(&rawheader.version), 			sizeof(rawheader.version) 		, 1);
      endian_swap((CHAR*)(&rawheader.gps_sec), 			sizeof(rawheader.gps_sec) 		, 1);
      endian_swap((CHAR*)(&rawheader.gps_nsec), 		sizeof(rawheader.gps_nsec) 		, 1);
      endian_swap((CHAR*)(&rawheader.tbase), 			sizeof(rawheader.tbase) 		, 1);
      endian_swap((CHAR*)(&rawheader.first_frequency_index), 	sizeof(rawheader.first_frequency_index) , 1);
      endian_swap((CHAR*)(&rawheader.nsamples), 		sizeof(rawheader.nsamples) 		, 1);
    }

  /* double-check version-number */
  if ( rawheader.version != 1 )
    {
      XLALPrintError ("\nWrong SFT-version %d in read_v1_header_from_fp()\n\n", rawheader.version );
      goto failed;
    }

  if ( rawheader.nsamples <= 0 )
    {
      XLALPrintError ("\nNon-positive number of samples in SFT!\n\n");
      goto failed;
    }


  /* ok: */
  memset ( header, 0, sizeof( *header ) );

  /* NOTE: v1-SFTs don't contain a detector-name, in which case we set it to '??' */
  strcpy ( header->name, "??" );

  header->epoch.gpsSeconds 	= rawheader.gps_sec;
  header->epoch.gpsNanoSeconds 	= rawheader.gps_nsec;
  header->deltaF 		= 1.0 / rawheader.tbase;
  header->f0 			= rawheader.first_frequency_index / rawheader.tbase;

  (*nsamples) = rawheader.nsamples;

  return 0;

 failed:
  /* restore filepointer initial position  */
  if ( fseek ( fp, save_filepos, SEEK_SET ) == -1 )
    XLALPrintError ("\nfseek() failed to return to intial fileposition: %s\n\n", strerror(errno) );

  return -1;

} /* read_v1_header_from_fp() */

/* check that channel-prefix defines a 'known' detector.  The list of
 * known detectors implemented here for now follows the list in
 * Appendix D of LIGO-T970130-F-E:
 *
 * returns TRUE if valid, FALSE otherwise */
static BOOLEAN
is_valid_detector (const char *channel)
{
  int i;
  const char *knownDetectors[] =
    {
      "A1",       /* ALLEGRO */
      "B1",       /* NIOBE */
      "E1",       /* EXPLORER */
      "G1",       /* GEO_600 */
      "H1",       /* LHO_4k */
      "H2",       /* LHO_2k */
      "K1",       /* ACIGA */
      "L1",       /* LLO_4k */
      "N1",       /* Nautilus */
      "O1",       /* AURIGA */
      "P1",       /* CIT_40 */
      "T1",       /* TAMA_300 */
      "V1",       /* Virgo_CITF */
      "V2",       /* Virgo (3km) */
      "Z1",	  /* LISA effective IFO 1 */
      "Z2",	  /* LISA effective IFO 2 */
      "Z3",	  /* LISA effective IFO 3 */
      "Z4",	  /* LISA effective IFO 2 minus 3 */
      "Z5",	  /* LISA effective IFO 3 minus 1 */
      "Z6",	  /* LISA effective IFO 1 minus 2 */
      "Z7",	  /* LISA pseudo TDI A */
      "Z8",	  /* LISA pseudo TDI E */
      "Z9",	  /* LISA pseudo TDI T */
      "X1",       /* RXTE PCA */
      "X2",       /* RXTE ASM */
      NULL
    };

  if ( !channel )
    return FALSE;

  if ( strlen(channel) < 2 )
    return FALSE;

  for ( i = 0; knownDetectors[i]; i ++ )
    {
      if ( ( knownDetectors[i][0] == channel[0] ) && ( knownDetectors[i][1] == channel[1] )  )
	return TRUE;
    }

  return FALSE;

} /* is_valid_detector() */


/* a little endian-swapper needed for SFT reading/writing */
static void
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

/*----------------------------------------------------------------------
 * glob() has been reported to fail under condor, so we use our own
 * function to get a filelist from a directory, using a glob-like pattern.
 * also can read a list of file names from a "list file".
 *
 * looks pretty ugly with all the #ifdefs for the Microsoft C compiler
 *
 * NOTE: the list of filenames is returned SORTED ALPHABETICALLY !
 *
 *----------------------------------------------------------------------*/
static LALStringVector *
find_files (const CHAR *globdir)
{
#ifndef _MSC_VER
  DIR *dir;
  struct dirent *entry;
#else
  intptr_t dir;
  struct _finddata_t entry;
  CHAR* ptr3;
#endif
  CHAR *dname;
  const CHAR *ptr1, *ptr2;
  CHAR *fpattern;
  size_t dirlen;
  CHAR **filelist = NULL;
  UINT4 numFiles = 0, newNumFiles = 0;
  LALStringVector *ret = NULL;
  UINT4 j;
  UINT4 namelen;
  CHAR *thisFname = NULL;

#define FILE_SEPARATOR ';'
  if ( (ptr2 = strchr (globdir, FILE_SEPARATOR)) )
    { /* globdir is multi-pattern ("pattern1;pattern2;pattern3") */
      /* call find_files() with every pattern found in globdir */

      ptr1 = (const CHAR*)globdir;
      while ( (ptr2 = strchr (ptr1, FILE_SEPARATOR)) )
	{
	  /* ptr1 points to the beginning of a pattern, ptr2 to the end */

	  /* copy the current name to thisFname */
	  namelen = ptr2 - ptr1;
	  if ((thisFname = LALRealloc(thisFname, (namelen+1)*sizeof(CHAR))) == NULL) {
	    for (j=0; j < numFiles; j++)
	      LALFree (filelist[j]);
	    if(filelist)
	      LALFree (filelist);
	    return(NULL);
	  }
	  strncpy(thisFname,ptr1,namelen);
	  thisFname[namelen] = '\0';

	  /* call find_files(thisFname) */
	  ret = find_files(thisFname);

	  /* append the output (if any) to the existing filelist */
	  if (ret) {
	    newNumFiles = numFiles + ret->length;

	    if ((filelist = LALRealloc (filelist, (newNumFiles) * sizeof(CHAR*))) == NULL) {
	      XLALDestroyStringVector(ret);
	      LALFree(thisFname);
	      return (NULL);
	    }

	    for(j=0; j < ret->length; j++)
	      filelist[numFiles+j] = ret->data[j];
	    LALFree(ret->data);
	    LALFree(ret);
	    numFiles = newNumFiles;
	  } else {
	    for (j=0; j < numFiles; j++)
	      LALFree (filelist[j]);
	    if(filelist)
	      LALFree (filelist);
	    LALFree(thisFname);
	    return(NULL);
	  }

	  /* skip the separator */
	  ptr1 = ptr2 + 1;
	} /* while */

      LALFree(thisFname);

      ret = find_files(ptr1);
      if (ret) {
	newNumFiles = numFiles + ret->length;

	if ((filelist = LALRealloc (filelist, (newNumFiles) * sizeof(CHAR*))) == NULL) {
	  XLALDestroyStringVector(ret);
	  return (NULL);
	}

	for(j=0; j < ret->length; j++)
	  filelist[numFiles+j] = ret->data[j];
	LALFree(ret->data);
	LALFree(ret);
	numFiles = newNumFiles;
      }

    } /* if multi-pattern */

  /* read list of file names from a "list file" */
#define LIST_PREFIX "list:"
  else if (strncmp(globdir, LIST_PREFIX, strlen(LIST_PREFIX)) == 0) {
    LALParsedDataFile *list = NULL;
    CHAR* listfname = NULL;

    /* create list file name
       prefix with "./" if not an absolute file name (see LALOpenDataFile()) */
    if ((listfname = LALCalloc(1, strlen(globdir) + 3)) == NULL) {
      return NULL;
    }
    ptr1 = globdir + strlen(LIST_PREFIX);
    if (*ptr1 == '/')
      *listfname = '\0';
    else
      strcpy(listfname, "./");
    strcat(listfname, ptr1);
#undef LIST_PREFIX

    /* read list of file names from file */
    if (XLALParseDataFile(&list, listfname) != XLAL_SUCCESS) {
      XLALPrintError("\n%s: Could not parse list file '%s'\n", __func__, listfname);
      return NULL;
    }

    /* allocate "filelist" */
    numFiles = list->lines->nTokens;
    if (numFiles == 0) {
      XLALPrintWarning("\n%s: List file '%s' contains no file names\n", __func__, listfname);
      LALFree(listfname);
      XLALDestroyParsedDataFile(list);
      return NULL;
    }
    if ((filelist = LALRealloc (filelist, numFiles * sizeof(CHAR*))) == NULL) {
      LALFree(listfname);
      XLALDestroyParsedDataFile(list);
      return NULL;
    }

    /* copy file names from "list" to "filelist" */
    for (j = 0; j < numFiles; ++j) {
      ptr1 = list->lines->tokens[j];

      /* these prefixes are added to file names by e.g. ligo_data_find */
#define FILE_PREFIX "file://localhost/"
      if (strncmp(ptr1, FILE_PREFIX, strlen(FILE_PREFIX)) == 0) {
	ptr1 += strlen(FILE_PREFIX) - 1;
      }
#undef FILE_PREFIX
      else
#define FILE_PREFIX "file:///"
      if (strncmp(ptr1, FILE_PREFIX, strlen(FILE_PREFIX)) == 0) {
	ptr1 += strlen(FILE_PREFIX) - 1;
      }
#undef FILE_PREFIX

      /* allocate "filelist", and cleanup if it fails  */
      if ((filelist[j] = LALCalloc(1, strlen(ptr1) + 1)) == NULL) {
	while (j-- > 0)
	  LALFree(filelist[j]);
	LALFree(filelist);
	LALFree(listfname);
	XLALDestroyParsedDataFile(list);
	return NULL;
      }

      /* copy string */
      strcpy(filelist[j], ptr1);

    }

    /* cleanup */
    LALFree(listfname);
    XLALDestroyParsedDataFile(list);

  } /* if list file */

  else if (is_pattern(globdir))

    { /* globdir is a single glob-style pattern */

      /* First we separate the globdir into directory-path and file-pattern */

#ifndef _WIN32
#define DIR_SEPARATOR '/'
#else
#define DIR_SEPARATOR '\\'
#endif

      /* any path specified or not ? */
      ptr1 = strrchr (globdir, DIR_SEPARATOR);
      if (ptr1)
	{ /* yes, copy directory-path */
	  dirlen = (size_t)(ptr1 - globdir) + 1;
	  if ( (dname = LALCalloc (1, dirlen)) == NULL)
	    return (NULL);
	  strncpy (dname, globdir, dirlen);
	  dname[dirlen-1] = '\0';

	  ptr1 ++; /* skip dir-separator */
	  /* copy the rest as a glob-pattern for matching */
	  if ( (fpattern = LALCalloc (1, strlen(ptr1) + 1)) == NULL )
	    {
	      LALFree (dname);
	      return (NULL);
	    }
	  strcpy (fpattern, ptr1);

	} /* if ptr1 */
      else /* no pathname given, assume "." */
	{
	  if ( (dname = LALCalloc(1, 2)) == NULL)
	    return (NULL);
	  strcpy (dname, ".");

	  if ( (fpattern = LALCalloc(1, strlen(globdir)+1)) == NULL)
	    {
	      LALFree (dname);
	      return (NULL);
	    }
	  strcpy (fpattern, globdir);	/* just file-pattern given */
	} /* if !ptr */


#ifndef _MSC_VER
      /* now go through the file-list in this directory */
      if ( (dir = opendir(dname)) == NULL) {
	XLALPrintError ("Can't open data-directory `%s`\n", dname);
	LALFree (dname);
	return (NULL);
      }
#else
      if ((ptr3 = (CHAR*)LALMalloc(strlen(dname)+3)) == NULL)
	return(NULL);
      sprintf(ptr3,"%s\\*",dname);
      dir = _findfirst(ptr3,&entry);
      LALFree(ptr3);
      if (dir == -1) {
	XLALPrintError ("Can't find file for pattern `%s`\n", ptr3);
	LALFree (dname);
	return (NULL);
      }
#endif

#ifndef _MSC_VER
      while ( (entry = readdir (dir)) != NULL )
#else
      do
#endif
	{
#ifndef _MSC_VER
	  thisFname = entry->d_name;
#else
	  thisFname = entry.name;
#endif

	  /* now check if glob-pattern fpattern matches the current filename */
	  if ( amatch(thisFname, fpattern)
	       /* and check if we didnt' match some obvious garbage like "." or ".." : */
	       && strcmp( thisFname, ".") && strcmp( thisFname, "..") )
	    {

	      numFiles ++;
	      if ( (filelist = LALRealloc (filelist, numFiles * sizeof(CHAR*))) == NULL) {
		LALFree (dname);
		LALFree (fpattern);
		return (NULL);
	      }

	      namelen = strlen(thisFname) + strlen(dname) + 2 ;

	      if ( (filelist[ numFiles - 1 ] = LALCalloc (1, namelen)) == NULL) {
		for (j=0; j < numFiles; j++)
		  LALFree (filelist[j]);
		LALFree (filelist);
		LALFree (dname);
		LALFree (fpattern);
		return (NULL);
	      }

	      sprintf(filelist[numFiles-1], "%s%c%s", dname, DIR_SEPARATOR, thisFname);

	    } /* if filename matched pattern */

	} /* while more directory entries */
#ifdef _MSC_VER
      while ( _findnext (dir,&entry) == 0 );
#endif

#ifndef _MSC_VER
      closedir (dir);
#else
      _findclose(dir);
#endif

      LALFree (dname);
      LALFree (fpattern);

    } /* if is_pattern */

  else

    { /* globdir is a single simple filename */
      /* add it to the list of filenames as it is */

      numFiles++;
      if ( (filelist = LALRealloc (filelist, numFiles * sizeof(CHAR*))) == NULL) {
	return (NULL);
      }
      namelen = strlen(globdir) + 1;
      if ( (filelist[ numFiles - 1 ] = LALCalloc (1, namelen)) == NULL) {
	LALFree (filelist);
	return (NULL);
      }
      strcpy(filelist[numFiles-1], globdir );
    }

  /* ok, did we find anything? */
  if (numFiles == 0)
    return (NULL);

  /* make a LALStringVector from the list of filenames */
  if ( (ret = LALCalloc (1, sizeof (LALStringVector) )) == NULL)
    {
      for (j=0; j<numFiles; j++)
	LALFree (filelist[j]);
      LALFree (filelist);
      return (NULL);
    }
  ret->length = numFiles;
  ret->data = filelist;

  /* sort this alphabetically (in-place) */
  if(numFiles>1)
    XLALSortStringVector (ret);

  return (ret);
} /* find_files() */

/* portable file-len function */
static long get_file_len ( FILE *fp )
{
  long save_fp;
  long len;

  if ( (save_fp = ftell(fp)) == -1 )
    return 0;

  if ( fseek ( fp, 0, SEEK_END ) == -1 )
    return 0;

  len = ftell(fp);

  if ( fseek ( fp, save_fp, SEEK_SET ) == -1 )
    return 0;

  return len;

} /* get_file_len() */





/* ----- the following function crc64() was taken from SFTReferenceLibrary.c
 * and adapted to LAL .
 *
 *  The quantity below is: D800000000000000 (base-16) =
 *  1101100000000000000000000000000000000000000000000000000000000000
 *  (base-2).  The primitive polynomial is x^64 + x^4 + x^3 + x + 1.
 */
#define POLY64 0xd800000000000000ULL
#define TABLELEN 256

/* The crc64 checksum of M bytes of data at address data is returned
 * by crc64(data, M, ~(0ULL)). Call the function multiple times to
 * compute the checksum of data made in contiguous chunks, setting
 * final argument to the previously accumulated checksum value. */
static UINT8
calc_crc64(const CHAR *data, UINT4 length, UINT8 crc)
{
  UINT8 CRCTable[TABLELEN];
  UINT4 i;

  /* is there is no data, simply return previous checksum value */
  if (!length || !data )
    return crc;

  /* initialize the CRC table for fast computation.  We could keep
     this table in memory to make the computation faster, but that is
     not re-entrant for multi-threaded code.
  */
  for (i = 0; i < TABLELEN; i++) {
    UINT4 j;
    UINT8 part = i;
    for (j = 0; j < 8; j++) {
      if (part & 1)
        part = (part >> 1) ^ POLY64;
      else
        part >>= 1;
    }
    CRCTable[i] = part;
  }

  /* compute the CRC-64 code */
  for (i=0; i<length; i++) {
    UINT8 temp1 = crc >> 8;
    UINT8 temp2 = CRCTable[(crc ^ (UINT8) data[i]) & 0xff];
    crc = temp1 ^ temp2;
  }

  return crc;

} /* calc_crc64() */

/* compare two SFT-descriptors by their GPS-epoch, then starting frequency */
static int
compareSFTdesc(const void *ptr1, const void *ptr2)
{
  const SFTDescriptor *desc1 = ptr1;
  const SFTDescriptor *desc2 = ptr2;

  if      ( GPS2REAL8( desc1->header.epoch ) < GPS2REAL8 ( desc2->header.epoch ) )
    return -1;
  else if ( GPS2REAL8( desc1->header.epoch ) > GPS2REAL8 ( desc2->header.epoch ) )
    return 1;
  else if ( desc1->header.f0 < desc2->header.f0 )
    return -1;
  else if ( desc1->header.f0 > desc2->header.f0 )
    return 1;
  else
    return 0;
} /* compareSFTdesc() */


/* compare two SFT-descriptors by their locator (f0, file, position) */
static int
compareSFTloc(const void *ptr1, const void *ptr2)
{
  const SFTDescriptor *desc1 = ptr1;
  const SFTDescriptor *desc2 = ptr2;
  int s;
  if ( desc1->header.f0 < desc2->header.f0 )
    return -1;
  else if ( desc1->header.f0 > desc2->header.f0 )
    return 1;
  s = strcmp(desc1->locator->fname, desc2->locator->fname);
  if(!s) {
    if (desc1->locator->offset < desc2->locator->offset)
      return(-1);
    else if (desc1->locator->offset > desc2->locator->offset)
      return(1);
    else
      return(0);
  }
  return(s);
} /* compareSFTloc() */


/* compare two SFT-vectors by detector name in alphabetic order */
static int
compareDetName(const void *ptr1, const void *ptr2)
{
  SFTVector const* const* sftvect1 = (SFTVector const* const*)ptr1;
  SFTVector const* const* sftvect2 = (SFTVector const* const*)ptr2;

  if ( (*sftvect1)->data[0].name[0] < (*sftvect2)->data[0].name[0] )
    return -1;
  else if ( (*sftvect1)->data[0].name[0] > (*sftvect2)->data[0].name[0] )
    return 1;
  else if ( (*sftvect1)->data[0].name[1] < (*sftvect2)->data[0].name[1] )
    return -1;
  else if ( (*sftvect1)->data[0].name[1] > (*sftvect2)->data[0].name[1] )
    return 1;
  else
    return 0;

} /* compareDetName() */



/** Read valid SFT version-number at position fp, and determine if we need to
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

 /* filename string is a glob-style pattern, i.e. it contains '*' or '?' or '[' */
static BOOLEAN is_pattern(const char*c) {
  while((*c != '\0') && (*c != '*') && (*c != '?') && (*c != '['))
    c++;
  return(*c != '\0');
}


/*======================================================================*/
/*
 * robust glob pattern matcher
 * ozan s. yigit/dec 1994
 * public domain
 *
 * glob patterns:
 *	*	matches zero or more characters
 *	?	matches any single character
 *	[set]	matches any character in the set
 *	[^set]	matches any character NOT in the set
 *		where a set is a group of characters or ranges. a range
 *		is written as two characters seperated with a hyphen: a-z denotes
 *		all characters between a to z inclusive.
 *	[-set]	set matches a literal hypen and any character in the set
 *	[]set]	matches a literal close bracket and any character in the set
 *
 *	char	matches itself except where char is '*' or '?' or '['
 *	\char	matches char, including any pattern character
 *
 * examples:
 *	a*c		ac abc abbc ...
 *	a?c		acc abc aXc ...
 *	a[a-z]c		aac abc acc ...
 *	a[-a-z]c	a-c aac abc ...
 *
 */

#ifndef NEGATE
#define NEGATE	'^'			/* std cset negation char */
#endif

static int amatch(char *str, char *p);

static int
amatch(char *str, char *p)
{
	int negate;
	int match;
	int c;

	while (*p) {
		if (!*str && *p != '*')
			return FALSE;

		switch (c = *p++) {

		case '*':
			while (*p == '*')
				p++;

			if (!*p)
				return TRUE;

			if (*p != '?' && *p != '[' && *p != '\\')
				while (*str && *p != *str)
					str++;

			while (*str) {
				if (amatch(str, p))
					return TRUE;
				str++;
			}
			return FALSE;

		case '?':
			if (*str)
				break;
			return FALSE;
/*
 * set specification is inclusive, that is [a-z] is a, z and
 * everything in between. this means [z-a] may be interpreted
 * as a set that contains z, a and nothing in between.
 */
		case '[':
			if (*p != NEGATE)
				negate = FALSE;
			else {
				negate = TRUE;
				p++;
			}

			match = FALSE;

			while (!match && (c = *p++)) {
				if (!*p)
					return FALSE;
				if (*p == '-') {	/* c-c */
					if (!*++p)
						return FALSE;
					if (*p != ']') {
						if (*str == c || *str == *p ||
						    (*str > c && *str < *p))
							match = TRUE;
					}
					else {		/* c-] */
						if (*str >= c)
							match = TRUE;
						break;
					}
				}
				else {			/* cc or c] */
					if (c == *str)
						match = TRUE;
					if (*p != ']') {
						if (*p == *str)
							match = TRUE;
					}
					else
						break;
				}
			}

			if (negate == match)
				return FALSE;
/*
 * if there is a match, skip past the cset and continue on
 */
			while (*p && *p != ']')
				p++;
			if (!*p++)	/* oops! */
				return FALSE;
			break;

		case '\\':
			if (*p)
				c = *p++;
		default:
			if (c != *str)
				return FALSE;
			break;

		}
		str++;
	}

	return !*str;
}
