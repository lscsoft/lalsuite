/*
 * Copyright (C) 2010, 2014, 2016, 2019, 2020, 2022 Karl Wette
 * Copyright (C) 2010 Chris Messenger
 * Copyright (C) 2009 Adam Mercer
 * Copyright (C) 2004--2008, 2013, 2014, 2017 Reinhard Prix
 * Copyright (C) 2004, 2005, 2009, 2010, 2017 Bernd Machenschalk
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

#include <sys/types.h>
#include <sys/stat.h>

#include <lal/Units.h>
#include <lal/Sequence.h>

#include "SFTinternal.h"

/*---------- internal prototypes ----------*/

static long get_file_len ( FILE *fp );

static BOOLEAN consistent_mSFT_header ( SFTtype header1, UINT4 version1, UINT4 nsamples1, SFTtype header2, UINT4 version2, UINT4 nsamples2 );
static BOOLEAN timestamp_in_list( LIGOTimeGPS timestamp, LIGOTimeGPSVector *list );

/*========== function definitions ==========*/

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
 * [minStartTime, maxStartTime) MUST be found!.
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
      if ( !XLALIsValidCWDetector ( constraints->detector ) )
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
      SFTtype XLAL_INIT_DECL( mprev_header );
      REAL8   mprev_nsamples = 0;

      FILE *fp;
      if ( ( fp = fopen( fname, "rb" ) ) == NULL )
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
	      if ( constraints->detector )
		{
                  if ( strncmp( constraints->detector, this_header.name, 2) ) {
		    want_this_block = FALSE;
                  }
		}

	      if ( XLALCWGPSinRange(this_header.epoch, constraints->minStartTime, constraints->maxStartTime) != 0 ) {
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

	  mfirst_block = FALSE;

	  /* skip seeking if we know we would reach the end */
	  if ( ftell ( fp ) + (long)this_nsamples * 8 >= file_len )
	    break;

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

  /* did we find all timestamps that lie within [minStartTime, maxStartTime)? */
  if ( constraints && constraints->timestamps )
    {
      LIGOTimeGPSVector *ts = constraints->timestamps;

      for ( UINT4 i = 0; i < ts->length; i ++ )
	{
          const LIGOTimeGPS *ts_i = &(ts->data[i]);
	  if ( XLALCWGPSinRange(*ts_i, constraints->minStartTime, constraints->maxStartTime) == 0 )
            {
              UINT4 j;
              for ( j = 0; j < ret->length; j ++ )
                {
                  const LIGOTimeGPS *sft_i = &(ret->data[j].header.epoch);
                  if ( (ts_i->gpsSeconds == sft_i->gpsSeconds) && ( ts_i->gpsNanoSeconds == sft_i->gpsNanoSeconds ) ) {
                    break;
                  }
                }
              XLAL_CHECK_NULL ( j < ret->length, XLAL_EFAILED,
                                "Timestamp %d : [%d, %d] did not find a matching SFT\n\n", (i+1), ts_i->gpsSeconds, ts_i->gpsNanoSeconds );
            }
	} // for i < ts->length

    } /* if constraints->timestamps */

  SFTtype XLAL_INIT_DECL(first_header);
  /* have all matched SFTs identical dFreq values ? */
  for ( UINT4 i = 0; i < ret->length; i ++ )
    {
      SFTtype this_header = ret->data[i].header;

      if ( i == 0 ) {
	first_header = this_header;
      }

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
      memcpy( name, catalog->data[k].header.name, 3*sizeof(CHAR) );

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
 * This function reads in the SFTs in the catalog and validates their CRC64 checksums.
 * The result of the validation is returned in '*crc_check'. The function itself returns
 * XLAL_SUCCESS if the operation suceeds (even if the checksums fail to validate),
 * and XLAL_FAILURE otherwise.
 *
 * \note: because this function has to read the complete SFT data into memory it is
 * potentially slow and memory-intensive.
 */
int
XLALCheckCRCSFTCatalog(
  BOOLEAN *crc_check,  /**< set to true if checksum validation passes */
  SFTCatalog *catalog  /**< catalog of SFTs to check */
  )
{

  XLAL_CHECK( crc_check != NULL, XLAL_EINVAL );
  XLAL_CHECK( catalog != NULL, XLAL_EINVAL );

  /* CRC checks are assumed to pass until one fails */
  *crc_check = 1;

  /* step through SFTs and check CRC64 */
  for ( UINT4 i=0; i < catalog->length; i ++ )
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
              return XLAL_FAILURE;
	    }
	  if ( !(has_valid_crc64 ( fp ) != 0) )
	    {
	      XLALPrintError ( "CRC64 checksum failure for SFT '%s'\n",
			      XLALshowSFTLocator ( catalog->data[i].locator ) );
              *crc_check = 0;
	      fclose(fp);
              return XLAL_SUCCESS;
	    }
	  fclose(fp);
	  break;

	default:
	  XLALPrintError ( "Illegal SFT-version encountered : %d\n", catalog->data[i].version );
          return XLAL_FAILURE;
	  break;
	} /* switch (version ) */

    } /* for i < numSFTs */

  return XLAL_SUCCESS;

} /* XLALCheckCRCSFTCatalog() */


/**
 * Return a sorted string vector listing the unique IFOs in the given catalog.
 */
LALStringVector *XLALListIFOsInCatalog( const SFTCatalog *catalog )
{
  XLAL_CHECK_NULL( catalog != NULL, XLAL_EFAULT );
  LALStringVector *ifos = NULL;
  for ( UINT4 k = 0; k < catalog->length; ++k )
    {
      char *name = XLALGetChannelPrefix( catalog->data[k].header.name );
      if ( XLALFindStringInVector( name, ifos ) < 0 )
        {
          ifos = XLALAppendString2Vector( ifos, name );
          XLAL_CHECK_NULL( ifos != NULL, XLAL_EFUNC );
        }
      XLALFree( name );
    }
  XLAL_CHECK_NULL( XLALSortStringVector( ifos ) == XLAL_SUCCESS, XLAL_EFUNC );
  return ifos;
} // XLALListIFOsInCatalog()


/**
 * Count the number of the unique IFOs in the given catalog.
 */
INT4 XLALCountIFOsInCatalog( const SFTCatalog *catalog )
{
  XLAL_CHECK( catalog != NULL, XLAL_EFAULT );
  LALStringVector *ifos = XLALListIFOsInCatalog( catalog );
  XLAL_CHECK( ifos != NULL, XLAL_EFUNC );
  UINT4 nifos = ifos->length;
  XLALDestroyStringVector( ifos );
  return nifos;
} // XLALCountIFOsInCatalog()


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
  XLAL_LAST_ELEM(ret) = 0;

  return ret;

} /* XLALshowSFTLocator() */


///
/// Set a SFT catalog 'slice' to a timeslice of a larger SFT catalog 'catalog', with entries
/// restricted to the interval ['minStartGPS','maxStartGPS') according to XLALCWGPSinRange().
/// The catalog 'slice' just points to existing data in 'catalog', and therefore should not
/// be deallocated.
///
int XLALSFTCatalogTimeslice(
  SFTCatalog *slice,			///< [out] Timeslice of SFT catalog
  const SFTCatalog *catalog,		///< [in] SFT catalog
  const LIGOTimeGPS *minStartGPS,	///< [in] Minimum starting GPS time
  const LIGOTimeGPS *maxStartGPS	///< [in] Maximum starting GPS time
  )
{

  // Check input
  XLAL_CHECK( slice != NULL, XLAL_EFAULT );
  XLAL_CHECK( catalog != NULL, XLAL_EFAULT );
  XLAL_CHECK( minStartGPS != NULL && maxStartGPS != NULL, XLAL_EFAULT );
  XLAL_CHECK( catalog->length > 0, XLAL_EINVAL );
  XLAL_CHECK( XLALGPSCmp( minStartGPS, maxStartGPS ) < 1 , XLAL_EINVAL , "minStartGPS (%"LAL_GPS_FORMAT") is greater than maxStartGPS (%"LAL_GPS_FORMAT")\n",
              LAL_GPS_PRINT(*minStartGPS), LAL_GPS_PRINT(*maxStartGPS) );

  // get a temporary timestamps vector with SFT epochs so we can call XLALFindTimesliceBounds()
  LIGOTimeGPSVector timestamps;
  timestamps.length = catalog->length;
  XLAL_CHECK ( (timestamps.data = XLALCalloc ( timestamps.length, sizeof(timestamps.data[0]) )) != NULL, XLAL_ENOMEM );
  for ( UINT4 i = 0; i < timestamps.length; i ++ ) {
    timestamps.data[i] = catalog->data[i].header.epoch;
  }

  UINT4 iStart, iEnd;
  XLAL_CHECK ( XLALFindTimesliceBounds ( &iStart, &iEnd, &timestamps, minStartGPS, maxStartGPS ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLALFree ( timestamps.data );

  // Initialise timeslice of SFT catalog
  XLAL_INIT_MEM(*slice);

  // If not empty: set timeslice of SFT catalog
  if ( iStart <= iEnd )
    {
      slice->length = iEnd - iStart + 1;
      slice->data = &catalog->data[iStart];
    }

  return XLAL_SUCCESS;

} // XLALSFTCatalogTimeslice()


//
// Like XLALSFTCatalogTimeslice, but actually returning the resulting SFTCatalog
// so that memory can be properly freed whenever using SWIG bindings.
//
SFTCatalog *XLALReturnSFTCatalogTimeslice(
  const SFTCatalog *catalog,		///< [in] SFT catalog
  const LIGOTimeGPS *minStartGPS,	///< [in] Minimum starting GPS time
  const LIGOTimeGPS *maxStartGPS	///< [in] Maximum starting GPS time
  )
{

  // Allocate the output SFTCatalog
  SFTCatalog *slice = NULL;
  XLAL_CHECK_NULL ( (slice = LALCalloc ( 1, sizeof (*slice) )) != NULL, XLAL_ENOMEM );

  XLAL_CHECK_NULL ( XLALSFTCatalogTimeslice(slice, catalog, minStartGPS, maxStartGPS) == XLAL_SUCCESS, XLAL_EFUNC );

  return slice;

} // XLALReturnSFTCatalogTimeslice()


/**
 * Create a 'fake' SFT catalog which contains only detector and timestamp information.
 */
SFTCatalog *
XLALAddToFakeSFTCatalog ( SFTCatalog *catalog,                          /**< [in] SFT catalog; if NULL, a new catalog is created */
                          const CHAR *detector,                         /**< [in] Name of detector to set fake catalog entries to */
                          const LIGOTimeGPSVector *timestamps           /**< [in] Timestamps of each fake catalog entry */
                          )
{

  // Check input
  XLAL_CHECK_NULL( detector != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( timestamps != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( timestamps->length > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL( timestamps->data != NULL, XLAL_EFAULT );

  // Get channel prefix
  CHAR *channel = XLALGetChannelPrefix( detector );
  XLAL_CHECK_NULL( channel != NULL, XLAL_EFUNC );

  // If catalog is NULL, create a new fake catalog
  if (catalog == NULL) {
    catalog = XLALCalloc(1, sizeof(*catalog));
    XLAL_CHECK_NULL( catalog != NULL, XLAL_ENOMEM );
  }

  // Extend catalog to add new timestamps
  const UINT4 new_length = catalog->length + timestamps->length;
  catalog->data = XLALRealloc(catalog->data, new_length * sizeof(catalog->data[0]));
  XLAL_CHECK_NULL( catalog->data != NULL, XLAL_ENOMEM );

  // Fill out new SFT descriptors with channel name and timestamps info
  for (UINT4 i = 0; i < timestamps->length; ++i) {
    SFTDescriptor *desc = &catalog->data[catalog->length + i];
    memset(desc, 0, sizeof(*desc));
    strncpy(desc->header.name, channel, 2);
    desc->header.epoch = timestamps->data[i];
    desc->header.deltaF = 1.0 / timestamps->deltaT;
    desc->header.sampleUnits = lalDimensionlessUnit;
    desc->version = 2;
  }

  // Set new catalog length
  catalog->length = new_length;

  // Sort catalog
  qsort( (void*)catalog->data, catalog->length, sizeof( catalog->data[0] ), compareSFTdesc );

  // Cleanup
  XLALFree(channel);

  return catalog;

} /* XLALAddToFakeSFTCatalog() */


/**
 * Multi-detector and multi-timestamp wrapper of XLALAddToFakeSFTCatalog().
 */
SFTCatalog *
XLALMultiAddToFakeSFTCatalog ( SFTCatalog *catalog,                          /**< [in] SFT catalog; if NULL, a new catalog is created */
                               const LALStringVector *detectors,             /**< [in] Detector names to set fake catalog entries to */
                               const MultiLIGOTimeGPSVector *timestamps      /**< [in] Timestamps for each detector of each fake catalog entry */
                               )
{
  // Check input
  XLAL_CHECK_NULL( detectors != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( detectors->length > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL( detectors->data != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( timestamps != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( timestamps->length == detectors->length, XLAL_EINVAL );
  XLAL_CHECK_NULL( timestamps->data != NULL, XLAL_EFAULT );

  // Loop over detectors, calling XLALAddToFakeSFTCatalog()
  for (UINT4 X = 0; X < detectors->length; ++X) {
    catalog = XLALAddToFakeSFTCatalog( catalog, detectors->data[X], timestamps->data[X] );
    XLAL_CHECK_NULL( catalog != NULL, XLAL_EFUNC );
  }

  return catalog;

} /* XLALMultiAddToFakeSFTCatalog() */


/* portable file-len function */
static long get_file_len ( FILE *fp )
{
#ifdef _WIN32

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

#else

  struct stat st;

  if ( fstat(fileno(fp), &st) )
    return 0;

  return st.st_size;

#endif
} /* get_file_len() */


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
