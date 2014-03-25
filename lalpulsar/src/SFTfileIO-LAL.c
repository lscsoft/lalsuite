/** \cond DONT_DOXYGEN */

static void lal_read_sft_bins_from_fp ( LALStatus *status, SFTtype **sft, UINT4 *binsread, UINT4 firstBin2read, UINT4 lastBin2read , FILE *fp );
static int compareDetName(const void *ptr1, const void *ptr2);
static BOOLEAN has_valid_v2_crc64 (FILE *fp );

/**
 * [DEPRECATED] Read timestamps file and returns timestamps vector (alloc'ed in here!).
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

  XLAL_PRINT_DEPRECATION_WARNING("XLALReadTimestampsFile()");

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


/**
 * \deprecated Use XLALSFTdataFind() instead.
 */
void
LALSFTdataFind ( LALStatus *status,			/**< pointer to LALStatus structure */
                 SFTCatalog **catalog,		/**< [out] SFT-catalogue of matching SFTs */
                 const CHAR *file_pattern,	/**< which SFT-files */
                 SFTConstraints *constraints	/**< additional constraints for SFT-selection */
                 )
{
  INITSTATUS(status);

  ASSERT ( catalog, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( (*catalog) == NULL, status, SFTFILEIO_ENONULL, SFTFILEIO_MSGENONULL );

  SFTCatalog *ret = XLALSFTdataFind ( file_pattern, constraints );
  if ( ret == NULL ) {
    XLALPrintError ("XLALSFTdataFind() failed with xlalErrno = %d\n", xlalErrno );
    ABORTXLAL(status);
  }

  (*catalog) = ret;

  RETURN(status);
} /* LALSFTdataFind() */


/**
 * Load the given frequency-band <tt>[fMin, fMax]</tt> (inclusively) from the SFT-files listed in the
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


/**
 * Function to load a catalog of SFTs from possibly different detectors.
 * This is similar to LALLoadSFTs except that the input SFT catalog is
 * allowed to contain multiple ifos.  The output is the structure
 * MultiSFTVector which is a vector of (pointers to) SFTVectors, one for
 * each ifo found in the catalog.   As in LALLoadSFTs, fMin and fMax can be
 * set to -1 to get the full SFT from the lowest to the highest frequency
 * bin found in the SFT.
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


void
LALWriteSFT2file (LALStatus *status,			/**< pointer to LALStatus structure */
		  const SFTtype *sft,		/**< SFT to write to disk */
		  const CHAR *fname,		/**< filename */
		  const CHAR *SFTcomment)	/**< optional comment (for v2 only) */
{
  XLAL_PRINT_DEPRECATION_WARNING("XLALWriteSFT2file");
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  if ( XLALWriteSFT2file( sft, fname, SFTcomment ) != XLAL_SUCCESS ) {
    ABORTXLAL(status);
  }
  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* WriteSFTtoFile() */


void
LALWriteSFTVector2Dir (LALStatus *status,			/**< pointer to LALStatus structure */
		       const SFTVector *sftVect,	/**< SFT vector to write to disk */
		       const CHAR *dirname,		/**< base filename (including directory path)*/
		       const CHAR *SFTcomment,		/**< optional comment (for v2 only) */
		       const CHAR *description)         /**< optional sft description to go in the filename */
{
  XLAL_PRINT_DEPRECATION_WARNING("XLALWriteSFTVector2Dir");
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  if ( XLALWriteSFTVector2Dir( sftVect, dirname, SFTcomment, description ) != XLAL_SUCCESS ) {
    ABORTXLAL(status);
  }
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/**
 * For backwards-compatibility: write a *v2-normalized* (ie dt x DFT) SFTtype
 * to a v1-SFT file.
 *
 * NOTE: the only difference to WriteSFTfile() is that the data-normalization
 * is changed back to v1-type 'DFT', by dividing the dt corresponding to the
 * frequency-band contained in the SFTtype.
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
      v1SFT.data->data[i] = crectf( (REAL4) ( (REAL8)crealf(v1SFT.data->data[i]) / dt ), (REAL4) ( (REAL8)cimagf(v1SFT.data->data[i]) / dt ) );
    }

  TRY ( LALWriteSFTfile (status->statusPtr, &v1SFT, fname ), status );

  XLALDestroyCOMPLEX8Vector ( v1SFT.data );

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* LALWrite_v2SFT_to_v1file() */


/**
 * Function to check validity of SFTs listed in catalog.
 * This function simply reads in those SFTs and checks their CRC64 checksum, which
 * is the only check that has not yet been done by the operations up to this point.
 *
 * Returns the LAL-return code of a failure (or 0 on success) in 'check_result'.
 *
 * \note: because this function has to read the complete SFT-data into memory for the
 * whole set of matching SFTs, it is potentially slow and memory-intensive.
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



/**
 * Extract a timstamps-vector from the given SFTCatalog.
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


/**
 * [DEPRECATED]: Low-level function to read only the SFT-header of a given file.
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


/**
 * [DEPRECATED] This is a function for low-level SFT data-reading:
 * the SFT-data is read starting from fminBinIndex and filled
 * into the pre-allocate vector sft of length N
 *
 * NOTE: !! NO re-normalization is done here!! this remains up
 * to the caller of this function!!
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
      sft->data->data[i] = crectf( rawdata[2 * i], rawdata[2 * i + 1] );
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


/**
 * [OBSOLETE] Write a *v1-normalized* (i.e. raw DFT) SFTtype to a SFT-v1 file.
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
      rawdata[2 * i]     = crealf(inData[i]);
      rawdata[2 * i + 1] = cimagf(inData[i]);
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


/**
 * [DEPRECATED] Basic SFT reading-function.
 * Given a filename \a fname and frequency-limits [\a fMin, \a fMax],
 * returns an SFTtype \a sft containing the SFT-data.
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
	outputSFT->data->data[i] *= ((REAL4) renorm);
      }

  /* that's it: return */
  *sft = outputSFT;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALReadSFTfile() */


/**
 * [DEPRECATED] Higher-level SFT-reading function to read a whole vector of SFT files
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
  if ( (fnames = XLALFindFiles (fpattern)) == NULL) {
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


/**
 * [DEPRECATED] Function to read and return a list of SFT-headers for given
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
  if ( (fnames = XLALFindFiles (fpattern)) == NULL) {
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

/**
 * Read bins from an SFT, leave filepointer at the end of the read SFT if successful,
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

  /* return resulting SFT */
  (*sft) = ret;

  DETATCHSTATUSPTR ( status );
  RETURN (status);
} /* lal_read_sft_bins_from_fp() */


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


/**
 * Check the v2 SFT-block starting at fp for valid crc64 checksum.
 * Restores filepointer before leaving.
 */
static BOOLEAN
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


/** \endcond */
