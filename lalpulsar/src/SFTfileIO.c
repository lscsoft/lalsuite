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
 */

/*---------- INCLUDES ----------*/
#include <sys/types.h>
#include <errno.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>

#ifndef _MSC_VER
#include <dirent.h>
#else
#include <io.h>
#endif

#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/SFTfileIO.h>
#include <lal/StringVector.h>
#include <lal/Sequence.h>
#include <lal/ConfigFile.h>
#include <lal/LogPrintf.h>
#include <lal/SFTutils.h>

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

/** segments read so far from one SFT */
typedef struct {
  UINT4 first;                     /**< first bin in this segment */
  UINT4 last;                      /**< last bin in this segment */
  LIGOTimeGPS epoch;               /**< timestamp of this SFT */
  struct tagSFTLocator *lastfrom;  /**< last bin read from this locator */
} SFTReadSegment;

/*---------- Global variables ----------*/
/* empty struct initializers */
static LALStatus empty_status;
const SFTConstraints empty_SFTConstraints;
const SFTCatalog empty_SFTCatalog;
const SFTtype empty_SFTtype;
const SFTVector empty_SFTVector;
const MultiSFTVector empty_MultiSFTVector;
const MultiREAL4TimeSeries empty_MultiREAL4TimeSeries;
const LIGOTimeGPSVector empty_LIGOTimeGPSVector;
const MultiLIGOTimeGPSVector empty_MultiLIGOTimeGPSVector;


static REAL8 fudge_up   = 1 + 10 * LAL_REAL8_EPS;	// about ~1 + 2e-15
static REAL8 fudge_down = 1 - 10 * LAL_REAL8_EPS;	// about ~1 - 2e-15


/*---------- internal prototypes ----------*/
static void endian_swap(CHAR * pdata, size_t dsize, size_t nelements);
static int amatch(char *str, char *p);	/* glob pattern-matcher (public domain)*/
static BOOLEAN is_pattern(const char*c); /* filename string is a glob-style pattern */

static BOOLEAN is_valid_detector (const char *channel);
static BOOLEAN consistent_mSFT_header ( SFTtype header1, UINT4 version1, UINT4 nsamples1, SFTtype header2, UINT4 version2, UINT4 nsamples2 );
static BOOLEAN timestamp_in_list( LIGOTimeGPS timestamp, LIGOTimeGPSVector *list );
static long get_file_len ( FILE *fp );

static FILE * fopen_SFTLocator ( const struct tagSFTLocator *locator );

static UINT4 read_sft_bins_from_fp ( SFTtype *ret, UINT4 *firstBinRead, UINT4 firstBin2read, UINT4 lastBin2read , FILE *fp );
static int read_sft_header_from_fp (FILE *fp, SFTtype  *header, UINT4 *version, UINT8 *crc64, BOOLEAN *swapEndian, CHAR **SFTcomment, UINT4 *numBins );
static int read_v2_header_from_fp ( FILE *fp, SFTtype *header, UINT4 *nsamples, UINT8 *header_crc64, UINT8 *ref_crc64, CHAR **SFTcomment, BOOLEAN swapEndian);
static int read_v1_header_from_fp ( FILE *fp, SFTtype *header, UINT4 *nsamples, BOOLEAN swapEndian);

static int compareSFTdesc(const void *ptr1, const void *ptr2);
static int compareSFTloc(const void *ptr1, const void *ptr2);
static int compareDetNameCatalogs ( const void *ptr1, const void *ptr2 );

static UINT8 calc_crc64(const CHAR *data, UINT4 length, UINT8 crc);
static int read_SFTversion_from_fp ( UINT4 *version, BOOLEAN *need_swap, FILE *fp );

/*==================== FUNCTION DEFINITIONS ====================*/

// ---------- obsolete LAL-API was moved into external file
#include "SFTfileIO-LAL.c"
// ------------------------------

/**
 * Find the list of SFTs matching the \a file_pattern and satisfying the given \a constraints,
 * return an \c SFTCatalog of the matching SFTs.
 *
 * The optional \a constraints that can be specified are (type SFTConstraints)
 * - 'detector':      	which detector
 * - 'time-span':    	GPS start- and end-times
 * - 'timestamps':    	list of GPS start-times
 *
 * ==> The returned SFTCatalog can be used directly as input to XLALLoadSFTs()
 * to load a single-IFO SFTVector, or XLALLoadMultiSFTs() to load a
 * multi-IFO vector of SFTVectors
 *
 * Except for the 'file_pattern' input, all the other constraints are optional
 * and can be passed as NULL (either globally constraings==NULL, or individually).
 *
 * Note that the constraints are combined by 'AND' and the resulting full constraint
 * MUST be satisfied (in particular: if 'timestamps' is given, all timestamps within
 * [startTime, endTime) MUST be found!.
 *
 * The returned SFTs in the catalogue are sorted by increasing GPS-epochs !
 *
 */
SFTCatalog *
XLALSFTdataFind ( const CHAR *file_pattern,		/**< which SFT-files */
                  const SFTConstraints *constraints	/**< additional constraints for SFT-selection */
                  )
{
  /* ----- check input */
  XLAL_CHECK_NULL ( file_pattern != NULL, XLAL_EINVAL );

  if ( constraints && constraints->detector )
    {
      if ( (strncmp(constraints->detector, "??", 2) != 0) && !is_valid_detector ( constraints->detector ) )
        {
          XLAL_ERROR_NULL ( XLAL_EDOM, "Invalid detector-constraint '%s'\n\n", constraints->detector );
        }
    }

  /* prepare return-catalog */
  SFTCatalog *ret;
  XLAL_CHECK_NULL ( (ret = LALCalloc ( 1, sizeof (*ret) )) != NULL, XLAL_ENOMEM );

  /* find matching filenames */
  LALStringVector *fnames;
  XLAL_CHECK_NULL ( (fnames = XLALFindFiles (file_pattern)) != NULL, XLAL_EFUNC, "Failed to find filelist matching pattern '%s'.\n\n", file_pattern );
  UINT4 numFiles = fnames->length;

  UINT4 numSFTs = 0;
  /* ----- main loop: parse all matching files */
  for ( UINT4 i = 0; i < numFiles; i ++ )
    {
      const CHAR *fname = fnames->data[i];

      /* merged SFTs need to satisfy stronger consistency-constraints (-> see spec) */
      BOOLEAN mfirst_block = TRUE;
      UINT4   mprev_version = 0;
      SFTtype mprev_header;
      REAL8   mprev_nsamples = 0;

      FILE *fp;
      if ( ( fp = LALFopen( fname, "rb" ) ) == NULL )
	{
          XLALPrintError ("ERROR: Failed to open matched file '%s'\n\n", fname );
	  XLALDestroyStringVector ( fnames );
	  XLALDestroySFTCatalog ( ret );
	  XLAL_ERROR_NULL ( XLAL_EIO );
	}

      long file_len;
      if ( (file_len = get_file_len(fp)) == 0 )
	{
          XLALPrintError ("ERROR: got file-len == 0 for '%s'\n\n", fname );
	  XLALDestroyStringVector ( fnames );
	  XLALDestroySFTCatalog ( ret );
	  fclose(fp);
	  XLAL_ERROR_NULL ( XLAL_EIO );
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
              XLALPrintError ("ERROR: ftell() failed for '%s'\n\n", fname );
	      XLALDestroyStringVector ( fnames );
	      XLALDestroySFTCatalog ( ret );
	      fclose (fp);
	      XLAL_ERROR_NULL ( XLAL_EIO );
	    }

	  if ( read_sft_header_from_fp (fp, &this_header, &this_version, &this_crc, &endian, &this_comment, &this_nsamples ) != 0 )
	    {
              XLALPrintError ("ERROR: File-block '%s:%ld' is not a valid SFT!\n\n", fname, ftell(fp));
	      XLALDestroyStringVector ( fnames );
	      XLALFree ( this_comment );
	      XLALDestroySFTCatalog ( ret );
	      fclose(fp);
	      XLAL_ERROR_NULL ( XLAL_EDATA );
	    }

	  /* if merged-SFT: check consistency constraints */
	  if ( !mfirst_block )
	    {
	      if ( ! consistent_mSFT_header ( mprev_header, mprev_version, mprev_nsamples, this_header, this_version, this_nsamples ) )
		{
                  XLALPrintError ( "ERROR: merged SFT-file '%s' contains inconsistent SFT-blocks!\n\n", fname);
		  XLALFree ( this_comment );
		  XLALDestroyStringVector ( fnames );
		  XLALDestroySFTCatalog ( ret );
		  fclose(fp);
		  XLAL_ERROR_NULL ( XLAL_EDATA );
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
		  if ( ! strncmp (this_header.name, "??", 2 ) ) {
		    strncpy ( this_header.name, constraints->detector, 2 );	/* SET to constraint! */
                  }
		  else if ( strncmp( constraints->detector, this_header.name, 2) ) {
		    want_this_block = FALSE;
                  }
		}

	      if ( constraints->startTime && ( GPS2REAL8(this_header.epoch) < GPS2REAL8( *constraints->startTime))) {
		want_this_block = FALSE;
              }

	      if ( constraints->endTime && ( GPS2REAL8(this_header.epoch) >= GPS2REAL8( *constraints->endTime ) ) ) {
		want_this_block = FALSE;
              }

	      if ( constraints->timestamps && !timestamp_in_list(this_header.epoch, constraints->timestamps) ) {
		want_this_block = FALSE;
              }

	    } /* if constraints */

	  if ( want_this_block )
	    {
	      numSFTs ++;

	      /* do we need to alloc more memory for the SFTs? */
	      if (  numSFTs > ret->length )
		{
		  /* we realloc SFT-memory blockwise in order to
		   * improve speed in debug-mode (using LALMalloc/LALFree)
		   */
                  int len = (ret->length + SFTFILEIO_REALLOC_BLOCKSIZE) * sizeof( *(ret->data) );
                  if ( (ret->data = LALRealloc ( ret->data, len )) == NULL )
		    {
                      XLALPrintError ("ERROR: SFT memory reallocation failed: nSFT:%d, len = %d\n", numSFTs, len );
		      XLALDestroyStringVector ( fnames );
		      XLALDestroySFTCatalog ( ret );
		      XLALFree ( this_comment );
		      fclose(fp);
		      XLAL_ERROR_NULL ( XLAL_ENOMEM );
		    }

		  /* properly initialize data-fields pointers to NULL to avoid SegV when Freeing */
		  for ( UINT4 j=0; j < SFTFILEIO_REALLOC_BLOCKSIZE; j ++ ) {
		    memset ( &(ret->data[ret->length + j]), 0, sizeof( ret->data[0] ) );
                  }

		  ret->length += SFTFILEIO_REALLOC_BLOCKSIZE;
		} // if numSFTs > ret->length

	      SFTDescriptor *desc = &(ret->data[numSFTs - 1]);

	      desc->locator = XLALCalloc ( 1, sizeof ( *(desc->locator) ) );
	      if ( desc->locator ) {
		desc->locator->fname = XLALCalloc( 1, strlen(fname) + 1 );
              }
	      if ( (desc->locator == NULL) || (desc->locator->fname == NULL ) )
		{
                  XLALPrintError ("ERROR: XLALCalloc() failed\n" );
		  XLALDestroyStringVector ( fnames );
		  XLALDestroySFTCatalog ( ret );
		  fclose(fp);
		  XLAL_ERROR_NULL ( XLAL_ENOMEM );
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
	      XLALFree ( this_comment );
	    }

	  /* seek to end of SFT data-entries in file  */
	  if ( fseek ( fp, this_nsamples * 8 , SEEK_CUR ) == -1 )
	    {
              XLALPrintError ("ERROR: Failed to skip DATA field for SFT '%s': %s\n", fname, strerror(errno) );
	      XLALFree ( this_comment );
	      XLALDestroyStringVector ( fnames );
	      XLALDestroySFTCatalog ( ret );
	      fclose(fp);
	      XLAL_ERROR_NULL ( XLAL_EIO );
	    }

	  mfirst_block = FALSE;

	} /* while !feof */

      fclose(fp);

    } /* for i < numFiles */

  /* free matched filenames */
  XLALDestroyStringVector ( fnames );

  /* now realloc SFT-vector (alloc'ed blockwise) to its *actual size* */
  int len;
  if ( (ret->data = XLALRealloc ( ret->data, len = numSFTs * sizeof( *(ret->data) ))) == NULL )
    {
      XLALDestroySFTCatalog ( ret );
      XLAL_ERROR_NULL ( XLAL_ENOMEM, "XLALRecalloc ( %d ) failed.\n", len );
    }
  ret->length = numSFTs;

  /* ----- final consistency-checks: ----- */

  /* did we find all timestamps that lie within [startTime, endTime)? */
  if ( constraints && constraints->timestamps )
    {
      LIGOTimeGPSVector *ts = constraints->timestamps;
      REAL8 t0, t1;
      if ( constraints->startTime ) {
	t0 = GPS2REAL8 ( (*constraints->startTime) );
      }
      else {
	t0 = 0;
      }
      if ( constraints->endTime ) {
	t1 = GPS2REAL8 ( (*constraints->endTime) );
      }
      else {
	t1 = LAL_REAL4_MAX;	/* large enough */
      }

      for ( UINT4 i = 0; i < ts->length; i ++ )
	{
          const LIGOTimeGPS *ts_i = &(ts->data[i]);
	  REAL8 ti = GPS2REAL8((*ts_i));
	  if ( (t0 <= ti) && ( ti < t1 ) )
            {
              UINT4 j;
              for ( j = 0; j < ret->length; j ++ )
                {
                  const LIGOTimeGPS *sft_i = &(ret->data[j].header.epoch);
                  if ( (ts_i->gpsSeconds == sft_i->gpsSeconds) && ( ts_i->gpsNanoSeconds == sft_i->gpsNanoSeconds ) ) {
                    break;
                  }
                }
              XLAL_CHECK_NULL ( j < ret->length, XLAL_EFAILED, "Timestamp %d : [%d, %d] did not find a matching SFT\n\n", (i+1), ts_i->gpsSeconds, ts_i->gpsNanoSeconds );
            } // if timestamp ti within startTime/endTime constraint
	} // for i < ts->length

    } /* if constraints->timestamps */

  SFTtype first_header = empty_SFTtype;
  /* have all matched SFTs identical dFreq values ? */
  for ( UINT4 i = 0; i < ret->length; i ++ )
    {
      SFTtype this_header = ret->data[i].header;

      if ( i == 0 ) {
	first_header = this_header;
      }

      /* dont give out v1-SFTs without detector-entry, except if constraint->detector="??" ! */
      if ( !constraints || !constraints->detector || strncmp(constraints->detector, "??", 2) )
	{
	  if ( !strncmp ( this_header.name, "??", 2 ) )
	    {
	      XLALDestroySFTCatalog ( ret );
	      XLAL_ERROR_NULL ( XLAL_EINVAL, "Pattern '%s' matched v1-SFTs but no detector-constraint given!\n\n", file_pattern);
	    }
	} /* if detector-constraint was not '??' */

      if ( this_header.deltaF != first_header.deltaF )
	{
	  XLALDestroySFTCatalog ( ret );
	  XLAL_ERROR_NULL ( XLAL_EDATA, "Pattern '%s' matched SFTs with inconsistent deltaF: %.18g != %.18g!\n\n",
                            file_pattern, this_header.deltaF, first_header.deltaF );
	}

    } /* for i < numSFTs */


  /* sort catalog in order of increasing GPS-time */
  qsort( (void*)ret->data, ret->length, sizeof( ret->data[0] ), compareSFTdesc );


  /* return result catalog (=sft-vect and locator-vect) */
  return ret;

} /* XLALSFTdataFind() */


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
	  REAL4 re = crealf(ret->data->data[i]);
	  REAL4 im = cimagf(ret->data->data[i]);

	  if ( swapEndian )
	    {
	      endian_swap( (CHAR *) &re, sizeof ( re ), 1 );
	      endian_swap( (CHAR *) &im, sizeof ( im ), 1 );
	    }

	  /* if the SFT-file was in v1-Format: need to renormalize the data now by 'Delta t'
	   * in order to follow the correct SFT-normalization
	   * (see LIGO-T040164-01-Z, and LIGO-T010095-00)
	   */
	  if ( version == 1 )
	    {
	      re *= dt;
	      im *= dt;
	    }

          ret->data->data[i] = crectf( re, im );
	} /* for i < numBins2read */
    } /* if SFT-v1 */

  /* return last bin read */
  return(lastBin2read);

} /* read_sft_bins_from_fp() */


/**
 * Load the given frequency-band <tt>[fMin, fMax]</tt> (inclusively) from the SFT-files listed in the
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
  XLALPrintInfo("%s: fMin: %f, fMax: %f, deltaF: %f, minbin: %u, maxbin: %u\n", __func__, fMin, fMax, deltaF, minbin, maxbin);

  /* calculate first and last frequency bin to read */
  if (fMin < 0)
    firstbin = minbin;
  else
    firstbin = (UINT4) floor (fMin / deltaF * fudge_up);	// round *down*, but allow for 10*eps 'fudge'
  if (fMax < 0)
    lastbin = maxbin;
  else {
    lastbin = (UINT4) ceil (fMax / deltaF * fudge_down);	// round *up*, but allow for 10*eps fudge
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
	if(locatalog.data[catPos].header.name && *(locatalog.data[catPos].header.name))
	  strcpy(sftVector->data[isft].name, locatalog.data[catPos].header.name);
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
 * This is similar to LALLoadSFTs except that the input SFT catalog is
 * allowed to contain multiple ifos. The output is the structure
 * MultiSFTVector which is a vector of (pointers to) SFTVectors, one for
 * each ifo found in the catalog. As in LALLoadSFTs, fMin and fMax can be
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
 * Load timestamps file 'fname' into LIGOTimeGPSVector struct, allocated here.
 *
 * The timestamps file is of the format: <repeated lines of the form "seconds nano-seconds">
 * allowing for '%#' as comments, which are ignored.
 *
 */
LIGOTimeGPSVector *
XLALReadTimestampsFile ( const CHAR *fname )
{
  /** check input consistency */
  XLAL_CHECK_NULL ( fname != NULL, XLAL_EINVAL );

  /* read and parse timestamps-list file contents*/
  LALParsedDataFile *flines = NULL;
  XLAL_CHECK_NULL ( XLALParseDataFile ( &flines, fname ) == XLAL_SUCCESS, XLAL_EFUNC );

  UINT4 numTS = flines->lines->nTokens;

  /* allocate and initialized segment list */
  LIGOTimeGPSVector *timestamps = NULL;
  XLAL_CHECK_NULL ( ( timestamps = XLALCreateTimestampVector ( numTS )) != NULL, XLAL_EFUNC );

  for ( UINT4 iTS = 0; iTS < numTS; iTS ++ )
    {
      INT4 secs, ns;
      char junk[11] = "";
      if ( sscanf ( flines->lines->tokens[iTS], "%d %d%10s\n", &secs, &ns, junk ) != 2 ) {
        if ( junk[0] != 0 ) { XLALPrintError ("Found junk '%s' after timestamp in line %d/%d\n", junk, iTS+1, numTS ); }
        XLALPrintError ( "Failed to parse data-line %d: '%s' in timestamps-file '%s' (needs to be in format 'sec ns')\n",
                         iTS + 1, flines->lines->tokens[iTS], fname );
        XLALDestroyTimestampVector ( timestamps );
        XLALDestroyParsedDataFile ( flines );
        XLAL_ERROR_NULL ( XLAL_ESYS );
      }

      XLAL_CHECK_NULL ( ( secs >= 0 ) && ( ns >= 0 ), XLAL_EDOM,
                        "Timestamps-file '%s' contains negative time-entry in line %d : s = %d, ns = %d\n", fname, iTS, secs, ns );

      XLAL_CHECK_NULL ( ns <= 999999999, XLAL_EDOM,
                        "Timestamps-file '%s' contains nano-seconds entry >= 1 billion ns in line %d: s = %d, ns = %d\n", fname, iTS, secs, ns );

      timestamps->data[iTS].gpsSeconds     = secs;
      timestamps->data[iTS].gpsNanoSeconds = ns;

    } /* for iTS < numTS */

  /* free parsed segment file contents */
  XLALDestroyParsedDataFile ( flines );

  return timestamps;

} /* XLALReadTimestampsFile() */


/**
 * Load several timestamps files, return a MultiLIGOTimeGPSVector struct, allocated here.
 *
 * The timestamps files are of the format: <repeated lines of the form "seconds nano-seconds">
 * allowing for '%#' as comments, which are ignored.
 *
 */
MultiLIGOTimeGPSVector *
XLALReadMultiTimestampsFiles ( const LALStringVector *fnames )
{
  XLAL_CHECK_NULL ( fnames != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( fnames->data != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( fnames->length > 0, XLAL_EDOM );

  UINT4 numDet = fnames->length;

  // ----- prepare output container
  MultiLIGOTimeGPSVector *multiTS;
  XLAL_CHECK_NULL ( ( multiTS = XLALCalloc ( 1, sizeof(*multiTS) )) != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL ( ( multiTS->data = XLALCalloc ( numDet, sizeof(multiTS->data[0]))) != NULL, XLAL_ENOMEM );
  multiTS->length = numDet;

  for ( UINT4 X=0; X < numDet; X ++ )
    {
      XLAL_CHECK_NULL ( fnames->data[X] != NULL, XLAL_EINVAL );
      XLAL_CHECK_NULL ( ( multiTS->data[X] = XLALReadTimestampsFile ( fnames->data[X] )) != NULL, XLAL_EFUNC );
    } // for X < numDet

  return multiTS;

} // XLALReadMultiTimestampsFiles()


/**
 * Write the given *v2-normalized* (i.e. dt x DFT) SFTtype to a FILE pointer.
 * Add the comment to SFT if SFTcomment != NULL.
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


/**
 * Write the given *v2-normalized* (i.e. dt x DFT) SFTtype to a v2-SFT file.
 * Add the comment to SFT if SFTcomment != NULL.
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


/**
 * Write the given *v2-normalized* (i.e. dt x DFT) SFTVector to a directory.
 * Add the comment to SFT if SFTcomment != NULL.
 *
 * NOTE: Currently this only supports writing v2-SFTs.
 * If you need to write a v1-SFT, you should use LALWriteSFTfile()
 *
 * Output SFTs have naming convention following LIGO-T040164-01
 */
int
XLALWriteSFTVector2Dir ( const SFTVector *sftVect,	/**< SFT vector to write to disk */
                         const CHAR *dirname,		/**< base filename (including directory path)*/
                         const CHAR *SFTcomment,	/**< optional comment */
                         const CHAR *Misc         	/**< optional 'Misc' field in SFT description (can be NULL) */
                         )
{
  XLAL_CHECK ( sftVect != NULL, XLAL_EINVAL );
  XLAL_CHECK ( sftVect->data != NULL, XLAL_EINVAL );
  XLAL_CHECK ( sftVect->length > 0, XLAL_EINVAL );
  XLAL_CHECK ( dirname != NULL, XLAL_EINVAL );

  UINT4 numSFTs = sftVect->length;

  for ( UINT4 k = 0; k < numSFTs; k++ )
    {
      SFTtype *sft = &(sftVect->data[k]);

      CHAR *filename;
      XLAL_CHECK ( (filename = XLALGetOfficialName4SFT ( sft, Misc )) != NULL, XLAL_EFUNC );

      CHAR *path;
      int len = strlen ( dirname ) + 1 + strlen ( filename ) + 1;
      XLAL_CHECK ( (path = XLALCalloc ( 1, len )) != NULL, XLAL_ENOMEM );
      sprintf ( path, "%s/%s", dirname, filename );
      XLAL_CHECK ( XLALWriteSFT2file( sft, path, SFTcomment ) == XLAL_SUCCESS, XLAL_EFUNC );

      XLALFree ( path );
      XLALFree ( filename );

    } // for k < numSFTs

  return XLAL_SUCCESS;

} /* XLALWriteSFTVector2Dir() */


/**
 * Write the given *v2-normalized* (i.e. dt x DFT) SFTVector to a single concatenated SFT file.
 * Add the comment to SFT if SFTcomment != NULL.
 *
 * NOTE: user specifies output directory, but the output SFT-filename follows the SFT-v2 naming convention,
 * see XLALOfficialSFTFilename() for details.
 */
int
XLALWriteSFTVector2File ( const SFTVector *sftVect,	//!< SFT vector to write to disk */
                          const CHAR *dirname,		//!< base filename (including directory path)*/
                          const CHAR *SFTcomment,	//!< optional comment (can be NULL) */
                          const CHAR *Misc         	//!< optional 'Misc' field in SFT description (can be NULL) */
                          )
{
  XLAL_CHECK ( sftVect != NULL, XLAL_EINVAL );
  XLAL_CHECK ( dirname != NULL, XLAL_EINVAL );

  char *filename;
  XLAL_CHECK ( (filename = XLALGetOfficialName4MergedSFTs ( sftVect, Misc )) != NULL, XLAL_EFUNC );

  CHAR *path;
  int len = strlen ( dirname ) + 1 + strlen ( filename ) + 1;
  XLAL_CHECK ( (path = XLALCalloc ( 1, len )) != NULL, XLAL_ENOMEM );
  sprintf ( path, "%s/%s", dirname, filename );

  XLAL_CHECK ( XLALWriteSFTVector2NamedFile( sftVect, path, SFTcomment ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLALFree ( path );
  XLALFree ( filename );

  return XLAL_SUCCESS;

} // XLALWriteSFTVector2File()


/**
 * Write the given *v2-normalized* (i.e. dt x DFT) SFTVector to a single concatenated SFT file.
 * Add the comment to SFT if SFTcomment != NULL.
 *
 * Allows specifying a filename for the output merged-SFT file.
 */
int
XLALWriteSFTVector2NamedFile ( const SFTVector *sftVect,	/**< SFT vector to write to disk */
                               const CHAR *filename,		/**< complete path+filename for concatenated SFT */
                               const CHAR *SFTcomment 		/**< optional comment */
                               )
{
  XLAL_CHECK ( sftVect != NULL, XLAL_EINVAL );
  XLAL_CHECK ( sftVect->data != NULL, XLAL_EINVAL );
  XLAL_CHECK ( sftVect->length > 0, XLAL_EINVAL );
  XLAL_CHECK ( filename != NULL, XLAL_EINVAL );

  UINT4 numSFTs = sftVect->length;

  /* open SFT-file for writing */
  FILE *fp;
  XLAL_CHECK ( (fp = LALFopen ( filename, "wb" )) != NULL, XLAL_EIO, "Failed to open '%s' for writing: %s\n\n", filename, strerror(errno));

  for ( UINT4 k = 0; k < numSFTs; k++ )
    {
      SFTtype *sft = &( sftVect->data[k] );

      /* write the k^th sft */
      XLAL_CHECK ( XLALWriteSFT2fp ( sft, fp, SFTcomment ) == XLAL_SUCCESS, XLAL_EFUNC );

    } // for k < numSFTs

  fclose(fp);

  return XLAL_SUCCESS;

} /* XLALWriteSFTVector2NamedFile() */


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


/**
 * Mostly for *debugging* purposes: provide a user-API to allow inspecting the SFT-locator
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

} // XLALCountIFOsInCatalog()


/**
 * Return a MultiSFTCatalogView generated from an input SFTCatalog.
 *
 * The input catalog can describe SFTs from several IFOs in one vector,
 * while the returned multi-Catalog view contains an array of single-IFO SFTCatalogs.
 *
 * NOTE: remember that this is only a multi-IFO "view" of the existing SFTCatalog,
 * various allocated memory of the original catalog is only pointed to, not duplicated!
 * This means one must not free the original catalog while this multi-view is still in use!
 *
 * NOTE2: the returned multi-IFO catalog is sorted alphabetically by detector-name
 *
 */
MultiSFTCatalogView *
XLALGetMultiSFTCatalogView ( const SFTCatalog *catalog )
{
  XLAL_CHECK_NULL ( catalog != NULL, XLAL_EINVAL );

  UINT4 numSFTsTotal = catalog->length;

  /* the number of ifos can be at most equal to numSFTsTotal */
  /* each ifo name is 2 characters + \0 */
  UINT4 numIFOsMax = 3; /* should be sufficient -- realloc used later in case required */
  UINT4 numIFOsMaxNew; /* for memory allocation purposes */

  CHAR  **ifolist;	/* list of ifo names */
  XLAL_CHECK_NULL ( (ifolist = XLALCalloc( 1, numIFOsMax * sizeof(*ifolist))) != NULL, XLAL_ENOMEM );

  UINT4 **sftLocationInCatalog;	/* location of sfts in catalog for each ifo */
  XLAL_CHECK_NULL ( (sftLocationInCatalog = XLALCalloc( 1, numIFOsMax * sizeof(*sftLocationInCatalog)) ) != NULL, XLAL_ENOMEM );

  UINT4  *numSFTsPerIFO;	/* number of sfts for each ifo 'X' */
  XLAL_CHECK_NULL ( (numSFTsPerIFO = XLALCalloc( 1, numIFOsMax * sizeof(*numSFTsPerIFO))) != NULL, XLAL_ENOMEM );

  for ( UINT4 X = 0; X < numIFOsMax; X++ )
    {
      XLAL_CHECK_NULL ( (ifolist[X] = XLALCalloc( 1, 3 * sizeof(*ifolist[X]) )) != NULL, XLAL_ENOMEM );
      XLAL_CHECK_NULL ( (sftLocationInCatalog[X] = XLALCalloc( 1, numSFTsTotal * sizeof(*sftLocationInCatalog[X]))) != NULL, XLAL_ENOMEM );
    } // for k < numIFOsMax

  UINT4 numIFOs = 0; /* number of ifos found so far */

  /* loop over sfts in catalog and look at ifo names and
   * find number of different ifos and number of sfts for each ifo
   * Also find location of sft in catalog */
  for ( UINT4 k = 0; k < numSFTsTotal; k++)
    {
      CHAR  name[3];
      strncpy( name, catalog->data[k].header.name, 3 );

      UINT4 X;
      /* go through list of ifos till a match is found or list is exhausted */
      for ( X = 0; ( X < numIFOs ) && strncmp( name, ifolist[X], 3); X++ )
	;

      if ( X < numIFOs )
	{
	  /* match found with X-th existing ifo */
	  sftLocationInCatalog[X][ numSFTsPerIFO[X] ] = k;
	  numSFTsPerIFO[X] ++;
	}
      else
	{
	  /* add ifo to list of ifos */

	  /* first check if number of ifos is larger than numIFOsmax */
	  /* and realloc if necessary */
	  if ( numIFOs >= numIFOsMax )
	    {
	      numIFOsMaxNew = numIFOsMax + 3;
	      XLAL_CHECK_NULL ( (ifolist = XLALRealloc( ifolist, numIFOsMaxNew * sizeof(*ifolist))) != NULL, XLAL_ENOMEM );

	      XLAL_CHECK_NULL ( (sftLocationInCatalog = XLALRealloc( sftLocationInCatalog, numIFOsMaxNew * sizeof(*sftLocationInCatalog))) != NULL, XLAL_ENOMEM );

	      XLAL_CHECK_NULL ( (numSFTsPerIFO = XLALRealloc( numSFTsPerIFO, numIFOsMaxNew * sizeof(*numSFTsPerIFO))) != NULL, XLAL_ENOMEM );

	      for ( UINT4 Y = numIFOsMax; Y < numIFOsMaxNew; Y++ )
                {
                  XLAL_CHECK_NULL ( (ifolist[Y] = XLALCalloc( 1,  3 * sizeof(*ifolist[Y]))) != NULL, XLAL_ENOMEM );
                  XLAL_CHECK_NULL ( (sftLocationInCatalog[Y] = XLALCalloc( 1, numSFTsTotal * sizeof(*sftLocationInCatalog[Y]))) != NULL, XLAL_ENOMEM );
                } // endfor X=numIFOsMax < numIFOsMaxNew

	      numIFOsMax = numIFOsMaxNew; // reset numIFOsMax

	    } /* endif ( numIFOs >= numIFOsMax) -- end of realloc */

	  strncpy( ifolist[numIFOs], name, 3);
	  sftLocationInCatalog[X][0] = k;
	  numSFTsPerIFO[numIFOs] = 1;
	  numIFOs ++;

	} /* else part of if ( X < numIFOs ) */

    } /*  for ( k = 0; k < numSFTsTotal; k++) */

  /* now we can create the return multi-SFT catalog view */
  MultiSFTCatalogView *ret;

  XLAL_CHECK_NULL ( (ret = XLALCalloc( 1, sizeof(*ret))) != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL ( (ret->data = XLALCalloc( numIFOs, sizeof(*ret->data))) != NULL, XLAL_ENOMEM );
  ret->length = numIFOs;

  for ( UINT4 X = 0; X < numIFOs; X++ )
    {
      ret->data[X].length = numSFTsPerIFO[X];

      XLAL_CHECK_NULL ( (ret->data[X].data = XLALCalloc( numSFTsPerIFO[X], sizeof(*(ret->data[X].data)) )) != NULL, XLAL_ENOMEM );

      for ( UINT4 k = 0; k < numSFTsPerIFO[X]; k++ )
	{
	  UINT4 location = sftLocationInCatalog[X][k];
	  ret->data[X].data[k] = catalog->data[location];	// struct copy, but keep all original pointers in struct!
	} // for k < numSFTsPerIFO[X]

    } // for X < numIFOs

  // free all temporary internal memory
  for ( UINT4 X = 0; X < numIFOsMax; X ++)
    {
      XLALFree ( ifolist[X] );
      XLALFree ( sftLocationInCatalog[X] );
    }
  XLALFree ( ifolist );
  XLALFree ( sftLocationInCatalog );
  XLALFree ( numSFTsPerIFO );

  // sort final multi-catalog view alphabetically by detector name
  qsort ( ret->data, ret->length, sizeof( ret->data[0] ), compareDetNameCatalogs );

  return ret;

} // XLALGetMultiSFTCatalogView()


/**
 * Destroys a MultiSFTCatalogView, without freeing the original
 * catalog that the 'view' was referring to, which
 * must be destroyed separately using XLALDestroySFTCatalog().
 */
void
XLALDestroyMultiSFTCatalogView ( MultiSFTCatalogView *multiView )
{
  if ( !multiView ) {
    return;
  }

  for ( UINT4 X = 0; X < multiView->length; X ++ )
    {
      XLALFree ( multiView->data[X].data );
    }

  XLALFree ( multiView->data );
  XLALFree ( multiView );

  return;

} // XLALDestroyMultiSFTCatalog()


/**
 * Return the 'official' file name for a given SFT, folllowing the SFT-v2 naming convention
 * LIGO-T040164-01 https://dcc.ligo.org/cgi-bin/DocDB/ShowDocument?docid=27385,
 * see also XLALOfficialSFTFilename() for details.
 */
char *
XLALGetOfficialName4SFT ( const SFTtype *sft,	//!< [in] input SFT to generate name for
                          const char *Misc	//!< [in] optional 'Misc' entry in the SFT 'D' field (can be NULL)
                          )
{
  XLAL_CHECK_NULL ( sft != NULL, XLAL_EINVAL );


  UINT4 Tsft = (UINT4) round ( 1.0 / sft->deltaF );

  /* calculate sft 'duration' -- may be different from timebase if nanosecond of sft-epoch is non-zero */
  UINT4 Tspan = Tsft;
  if ( sft->epoch.gpsNanoSeconds > 0) {
    Tspan += 1;
  }

  char *filename;
  XLAL_CHECK_NULL ( (filename = XLALOfficialSFTFilename ( sft->name[0], sft->name[1], 1, Tsft, sft->epoch.gpsSeconds, Tspan, Misc )) != NULL, XLAL_EFUNC );

  return filename;

} // XLALGetOfficialName4SFT()


/**
 * Return the 'official' file name for a given SFT-vector written into a single "merged SFT-file",
 * folllowing the SFT-v2 naming convention
 * LIGO-T040164-01 https://dcc.ligo.org/cgi-bin/DocDB/ShowDocument?docid=27385,
 * see also XLALOfficialSFTFilename() for details.
 */
char *
XLALGetOfficialName4MergedSFTs ( const SFTVector *sfts,	//!< [in] input SFT vector to generate name for
                                 const char *Misc	//!< [in] optional 'Misc' entry in the SFT 'D' field (can be NULL)
                                 )
{
  XLAL_CHECK_NULL ( sfts != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( sfts->length > 0, XLAL_EINVAL );

  UINT4 numSFTs = sfts->length;
  SFTtype *sftStart       = &(sfts->data[0]);
  SFTtype *sftEnd         = &(sfts->data[numSFTs-1]);
  LIGOTimeGPS *epochStart = &(sftStart->epoch);
  LIGOTimeGPS *epochEnd   = &(sftEnd->epoch);

  const char *name = sftStart->name;
  UINT4 Tsft = (UINT4) round ( 1.0 / sftStart->deltaF );

  /* calculate time interval covered -- may be different from timebase if nanosecond of sft-epochs are non-zero */
  UINT4 Tspan = epochEnd->gpsSeconds - epochStart->gpsSeconds + Tsft;
  if ( epochStart->gpsNanoSeconds > 0) {
    Tspan += 1;
  }
  if ( epochEnd->gpsNanoSeconds > 0) {
    Tspan += 1;
  }

  char *filename;
  XLAL_CHECK_NULL ( (filename = XLALOfficialSFTFilename ( name[0], name[1], numSFTs, Tsft, epochStart->gpsSeconds, Tspan, Misc )) != NULL, XLAL_EFUNC );

  return filename;

} // XLALGetOfficialName4MergedSFTs()


/**
 * Return the 'official' file name for a given SFT, folllowing the SFT-v2 naming convention
 * LIGO-T040164-01 https://dcc.ligo.org/cgi-bin/DocDB/ShowDocument?docid=27385, namely
 *
 * name = S-D-G-T.sft
 * where
 * S = Source: upper-case single letter site designation 'G', 'H', 'L', 'V', ...
 * D = description: a free-form string of alphanumerics and {_, +, #}
 * G = GPS start time of first SFT in seconds (9- or 10-digit number)
 * T = total time interval covered by the data in this file
 *
 * furthermore, the v2-spec uses the following convention for the description field 'D':
 * D = numSFTs_IFO_SFTtype[_Misc]
 * where
 * numSFTs : number of SFTs in the file
 * IFO     : 2-character detector name, eg 'G1', 'H1', 'H2', 'L1', 'V1', ...
 * SFTtype : SFT-timebase, in the form '[T]SFT', where [T] is the SFT-duration in seconds, eg "1800SFT"
 * Misc    : optional string providing additional information
 */
char *
XLALOfficialSFTFilename ( char site,		//!< site-character 'G', 'H', 'L', ...
                          char channel,	//!< channel character '1', '2', ...
                          UINT4 numSFTs,	//!< number of SFTs in SFT-file
                          UINT4 Tsft,		//!< time-baseline in (integer) seconds
                          UINT4 GPS_start,	//!< GPS seconds of first SFT start time
                          UINT4 Tspan,		//!< total time-spanned by all SFTs in seconds
                          const char *Misc	//!< [in] optional 'Misc' entry in the SFT 'D' field (can be NULL)
                          )
{
  if ( Misc != NULL ) {
    XLAL_CHECK_NULL ( XLALCheckValidDescriptionField ( Misc ) == XLAL_SUCCESS, XLAL_EINVAL );
  }

  // ----- S
  char S[2] = { site, 0 };

  // ----- D
  char D[512];
  char IFO[2] = { site, channel };
  size_t written = snprintf ( D, sizeof(D), "%d_%c%c_%dSFT%s%s", numSFTs, IFO[0], IFO[1], Tsft, Misc ? "_" : "", Misc ? Misc : "" );
  XLAL_CHECK_NULL ( written < sizeof(D), XLAL_EINVAL, "Description field length of %d exceeds buffer length of %d characters\n", written, sizeof(D)-1 );

  // ----- G
  char G[11];
  written = snprintf ( G, sizeof(G), "%09d", GPS_start );
  XLAL_CHECK_NULL ( written < sizeof(G), XLAL_EINVAL, "GPS seconds %d exceed buffer length of %d characters\n", GPS_start, sizeof(G)-1 );

  // ----- T
  char T[10];
  written = snprintf ( T, sizeof(T), "%d", Tspan );
  XLAL_CHECK_NULL ( written < sizeof(T), XLAL_EINVAL, "Tspan=%d s exceed buffer length of %d characters\n", Tspan, sizeof(T)-1 );

  // S-D-G-T.sft
  size_t len = strlen(S) + 1 + strlen(D) + 1 + strlen(G) + 1 + strlen(T) + 4 + 1;
  char *filename;
  XLAL_CHECK_NULL ( (filename = XLALCalloc ( 1, len )) != NULL, XLAL_ENOMEM );

  written = snprintf ( filename, len, "%s-%s-%s-%s.sft", S, D, G, T );
  XLAL_CHECK_NULL ( written < len, XLAL_EFAILED, "Miscounted string-length, expected %d characters but got %d\n", len - 1, written );

  return filename;

} // XLALGetOfficialName4SFT()


/**
 * Check whether given string qualifies as a valid 'description' field of a FRAME (or SFT)
 * filename, according to  LIGO-T010150-00-E "Naming Convention for Frame Files which are to be Processed by LDAS",
 * LIGO-T040164-01 at https://dcc.ligo.org/LIGO-T040164-x0/public
 *
 */
int
XLALCheckValidDescriptionField ( const char *desc )
{
  XLAL_CHECK ( desc != NULL, XLAL_EINVAL );

  size_t len = strlen ( desc );

  if ( len == 1 && isupper(desc[0]) ) {
    XLAL_ERROR ( XLAL_EINVAL, "Single uppercase description reserved for class-1 raw frames!\n" );
  }

  for ( UINT4 i=0; i < len; i ++ )
    {
      int c = desc[i];
      if ( !isalnum(c) && (c!='_') && (c!='+') && (c!='#') ) {	// all the valid characters allowed
        XLAL_ERROR ( XLAL_EINVAL, "Invalid chacter '%c' found, only alphanumeric and ['_', '+', '#'] are allowed\n", c );
      }
    } // for i < len

  return XLAL_SUCCESS;

} // XLALCheckValidDescriptionField()


/*================================================================================
 * LOW-level internal SFT-handling functions, should *NOT* be used outside this file!
 *================================================================================*/


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
static FILE *
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
LALStringVector *
XLALFindFiles (const CHAR *globstring)
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

  XLAL_CHECK_NULL ( globstring != NULL, XLAL_EINVAL );

#define FILE_SEPARATOR ';'
  if ( (ptr2 = strchr (globstring, FILE_SEPARATOR)) )
    { /* globstring is multi-pattern ("pattern1;pattern2;pattern3") */
      /* call XLALFindFiles() with every pattern found in globstring */

      ptr1 = (const CHAR*)globstring;
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
	    XLAL_ERROR_NULL ( XLAL_ENOMEM );
	  }
	  strncpy(thisFname,ptr1,namelen);
	  thisFname[namelen] = '\0';

	  /* call XLALFindFiles(thisFname) */
	  ret = XLALFindFiles(thisFname);

	  /* append the output (if any) to the existing filelist */
	  if (ret) {
	    newNumFiles = numFiles + ret->length;

	    if ((filelist = LALRealloc (filelist, (newNumFiles) * sizeof(CHAR*))) == NULL) {
	      XLALDestroyStringVector(ret);
	      LALFree(thisFname);
	      XLAL_ERROR_NULL ( XLAL_ENOMEM );
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
	    XLAL_ERROR_NULL ( XLAL_EFUNC);
	  }

	  /* skip the separator */
	  ptr1 = ptr2 + 1;
	} /* while */

      LALFree(thisFname);

      ret = XLALFindFiles(ptr1);
      if (ret) {
	newNumFiles = numFiles + ret->length;

	if ((filelist = LALRealloc (filelist, (newNumFiles) * sizeof(CHAR*))) == NULL) {
	  XLALDestroyStringVector(ret);
	  XLAL_ERROR_NULL ( XLAL_ENOMEM );
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
  else if (strncmp(globstring, LIST_PREFIX, strlen(LIST_PREFIX)) == 0) {
    LALParsedDataFile *list = NULL;
    CHAR* listfname = NULL;

    /* create list file name
       prefix with "./" if not an absolute file name (see LALOpenDataFile()) */
    if ((listfname = LALCalloc(1, strlen(globstring) + 3)) == NULL) {
      XLAL_ERROR_NULL ( XLAL_ENOMEM ) ;
    }
    ptr1 = globstring + strlen(LIST_PREFIX);
    if (*ptr1 == '/')
      *listfname = '\0';
    else
      strcpy(listfname, "./");
    strcat(listfname, ptr1);
#undef LIST_PREFIX

    /* read list of file names from file */
    if (XLALParseDataFile(&list, listfname) != XLAL_SUCCESS) {
      XLAL_ERROR_NULL ( XLAL_EFUNC, "Could not parse list file '%s'\n",listfname );
    }

    /* allocate "filelist" */
    numFiles = list->lines->nTokens;
    if (numFiles == 0) {
      XLALPrintWarning("\n%s: List file '%s' contains no file names\n", __func__, listfname);
      LALFree(listfname);
      XLALDestroyParsedDataFile(list);
      XLAL_ERROR_NULL ( XLAL_EINVAL );
    }
    if ((filelist = LALRealloc (filelist, numFiles * sizeof(CHAR*))) == NULL) {
      LALFree(listfname);
      XLALDestroyParsedDataFile(list);
      XLAL_ERROR_NULL ( XLAL_ENOMEM );
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
	XLAL_ERROR_NULL ( XLAL_ENOMEM );
      }

      /* copy string */
      strcpy(filelist[j], ptr1);

    }

    /* cleanup */
    LALFree(listfname);
    XLALDestroyParsedDataFile(list);

  } /* if list file */

  else if (is_pattern(globstring))

    { /* globstring is a single glob-style pattern */

      /* First we separate the globstring into directory-path and file-pattern */

#ifndef _WIN32
#define DIR_SEPARATOR '/'
#else
#define DIR_SEPARATOR '\\'
#endif

      /* any path specified or not ? */
      ptr1 = strrchr (globstring, DIR_SEPARATOR);
      if (ptr1)
	{ /* yes, copy directory-path */
	  dirlen = (size_t)(ptr1 - globstring) + 1;
	  if ( (dname = LALCalloc (1, dirlen)) == NULL)
	    XLAL_ERROR_NULL ( XLAL_ENOMEM );
	  strncpy (dname, globstring, dirlen);
	  dname[dirlen-1] = '\0';

	  ptr1 ++; /* skip dir-separator */
	  /* copy the rest as a glob-pattern for matching */
	  if ( (fpattern = LALCalloc (1, strlen(ptr1) + 1)) == NULL )
	    {
	      LALFree (dname);
	      XLAL_ERROR_NULL ( XLAL_ENOMEM );
	    }
	  strcpy (fpattern, ptr1);

	} /* if ptr1 */
      else /* no pathname given, assume "." */
	{
	  if ( (dname = LALCalloc(1, 2)) == NULL)
            XLAL_ERROR_NULL ( XLAL_ENOMEM );
	  strcpy (dname, ".");

	  if ( (fpattern = LALCalloc(1, strlen(globstring)+1)) == NULL)
	    {
	      LALFree (dname);
              XLAL_ERROR_NULL ( XLAL_ENOMEM );
	    }
	  strcpy (fpattern, globstring);	/* just file-pattern given */
	} /* if !ptr */


#ifndef _MSC_VER
      /* now go through the file-list in this directory */
      if ( (dir = opendir(dname)) == NULL) {
	XLALPrintError ("Can't open data-directory `%s`\n", dname);
	LALFree (dname);
        XLAL_ERROR_NULL ( XLAL_EIO );
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
        XLAL_ERROR_NULL ( XLAL_EIO );
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
                XLAL_ERROR_NULL ( XLAL_ENOMEM );
	      }

	      namelen = strlen(thisFname) + strlen(dname) + 2 ;

	      if ( (filelist[ numFiles - 1 ] = LALCalloc (1, namelen)) == NULL) {
		for (j=0; j < numFiles; j++)
		  LALFree (filelist[j]);
		LALFree (filelist);
		LALFree (dname);
		LALFree (fpattern);
                XLAL_ERROR_NULL ( XLAL_ENOMEM );
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

    { /* globstring is a single simple filename */
      /* add it to the list of filenames as it is */

      numFiles++;
      if ( (filelist = LALRealloc (filelist, numFiles * sizeof(CHAR*))) == NULL) {
        XLAL_ERROR_NULL ( XLAL_ENOMEM );
      }
      namelen = strlen(globstring) + 1;
      if ( (filelist[ numFiles - 1 ] = LALCalloc (1, namelen)) == NULL) {
	LALFree (filelist);
        XLAL_ERROR_NULL ( XLAL_ENOMEM );
      }
      strcpy(filelist[numFiles-1], globstring );
    }

  /* ok, did we find anything? */
  if (numFiles == 0)
    XLAL_ERROR_NULL ( XLAL_EINVAL );


  /* make a LALStringVector from the list of filenames */
  if ( (ret = LALCalloc (1, sizeof (LALStringVector) )) == NULL)
    {
      for (j=0; j<numFiles; j++)
	LALFree (filelist[j]);
      LALFree (filelist);
      XLAL_ERROR_NULL ( XLAL_ENOMEM );
    }
  ret->length = numFiles;
  ret->data = filelist;

  /* sort this alphabetically (in-place) */
  if(numFiles>1)
    XLALSortStringVector (ret);

  return (ret);

} /* XLALFindFiles() */


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


/* compare two SFT-catalog by detector name in alphabetic order */
static int
compareDetNameCatalogs ( const void *ptr1, const void *ptr2 )
{
  SFTCatalog const* cat1 = (SFTCatalog const*)ptr1;
  SFTCatalog const* cat2 = (SFTCatalog const*)ptr2;
  const char *name1 = cat1->data[0].header.name;
  const char *name2 = cat2->data[0].header.name;

  if ( name1[0] < name2[0] )
    return -1;
  else if ( name1[0] > name2[0] )
    return 1;
  else if ( name1[1] < name2[1] )
    return -1;
  else if ( name1[1] > name2[1] )
    return 1;
  else
    return 0;

} /* compareDetNameCatalogs() */


/**
 * Read valid SFT version-number at position fp, and determine if we need to
 * endian-swap the data.
 * Restores filepointer to original position before returning.
 *
 * RETURN: 0 = OK, -1 = ERROR
 */
static int
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
