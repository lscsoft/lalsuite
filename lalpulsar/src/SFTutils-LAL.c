/** \cond DONT_DOXYGEN */

/**
 * \deprecated Use XLALCreateSFT() instead
 * Allows for numBins == 0.
 */
void
LALCreateSFTtype (LALStatus *status,	/**< pointer to LALStatus structure */
		  SFTtype **output, 	/**< [out] allocated SFT-struct */
		  UINT4 numBins)	/**< number of frequency-bins */
{
  SFTtype *sft = NULL;

  INITSTATUS(status);

  ASSERT (output != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT (*output == NULL, status, SFTUTILS_ENONULL,  SFTUTILS_MSGENONULL);

  if ( (sft = XLALCreateSFT ( numBins )) == NULL )
    {
      XLALPrintError ("XLALCreateSFT() failed with xlalErrno = %d\n", xlalErrno );
      ABORT ( status, SFTUTILS_EFUNC, SFTUTILS_MSGEFUNC );
    }

  (*output) = sft;

  RETURN (status);

} /* LALCreateSFTtype() */

/**
 * \deprecated Use XLALCreateSFTVector() instead
 */
void
LALCreateSFTVector (LALStatus *status,	/**< pointer to LALStatus structure */
		    SFTVector **output, /**< [out] allocated SFT-vector */
		    UINT4 numSFTs, 	/**< number of SFTs */
		    UINT4 numBins)	/**< number of frequency-bins per SFT */
{
  INITSTATUS(status);

  ASSERT (output != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT (*output == NULL, status, SFTUTILS_ENONULL,  SFTUTILS_MSGENONULL);

  SFTVector *vect;
  if ( (vect = XLALCreateSFTVector ( numSFTs, numBins )) == NULL )
    {
      XLALPrintError ("XLALCreateSFTVector() failed with xlalErrno = %d\n", xlalErrno );
      ABORT (status, SFTUTILS_EFUNC, SFTUTILS_MSGEFUNC);
    }

  (*output) = vect;

  RETURN (status);

} /* LALCreateSFTVector() */

/**
 * Create an empty multi-IFO SFT vector for given number of IFOs and number of SFTs per IFO
 */
void LALCreateMultiSFTVector ( LALStatus *status,     /**< pointer to LALStatus structure */
			       MultiSFTVector **out,  /**< [out] multi sft vector created */
			       UINT4 length,          /**< number of sft data points */
			       UINT4Vector *numsft    /**< number of sfts in each sftvect */
			       )
{

  UINT4 k, j, numifo;
  MultiSFTVector *multSFTVec=NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( out, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL );
  ASSERT ( *out == NULL, status, SFTUTILS_ENONULL, SFTUTILS_MSGENONULL );
  ASSERT ( length, status, SFTUTILS_EINPUT, SFTUTILS_MSGEINPUT );
  ASSERT ( numsft, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL );
  ASSERT ( numsft->length > 0, status, SFTUTILS_EINPUT, SFTUTILS_MSGEINPUT );
  ASSERT ( numsft->data, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL );

  if ( (multSFTVec = (MultiSFTVector *)LALCalloc(1, sizeof(MultiSFTVector))) == NULL){
    ABORT ( status, SFTUTILS_EMEM, SFTUTILS_MSGEMEM );
  }

  numifo = numsft->length;
  multSFTVec->length = numifo;

  if ( (multSFTVec->data = (SFTVector **)LALCalloc( 1, numifo*sizeof(SFTVector *))) == NULL) {
    ABORT ( status, SFTUTILS_EMEM, SFTUTILS_MSGEMEM );
  }

  for ( k = 0; k < numifo; k++) {
    LALCreateSFTVector (status->statusPtr, multSFTVec->data + k, numsft->data[k], length);
      BEGINFAIL ( status ) {
	for ( j = 0; j < k-1; j++)
	  LALDestroySFTVector ( status->statusPtr, multSFTVec->data + j );
	LALFree( multSFTVec->data);
	LALFree( multSFTVec);
      } ENDFAIL(status);
  } /* loop over ifos */

  *out = multSFTVec;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALCreateMultiSFTVector() */

/**
 * \deprecated Use XLALDestroySFT() instead.
 */
void
LALDestroySFTtype (LALStatus *status,	/**< pointer to LALStatus structure */
		   SFTtype **sft)	/**< SFT-struct to free */
{

  INITSTATUS(status);

  ASSERT (sft != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  XLALDestroySFT ( (*sft) );

  (*sft) = NULL;

  RETURN (status);

} /* LALDestroySFTtype() */

/**
 * \deprecated Use XLALDestroySFTVector() instead.
 */
void
LALDestroySFTVector (LALStatus *status,	/**< pointer to LALStatus structure */
		     SFTVector **vect)	/**< the SFT-vector to free */
{
  INITSTATUS(status);

  ASSERT (vect != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  XLALDestroySFTVector ( *vect );

  (*vect) = NULL;

  RETURN (status);

} /* LALDestroySFTVector() */

/**
 * \deprecated Use XLALDestroyPSDVector() instead
 */
void
LALDestroyPSDVector (LALStatus *status,	/**< pointer to LALStatus structure */
		     PSDVector **vect)	/**< the SFT-vector to free */
{
  INITSTATUS(status);

  ASSERT (vect != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  XLALDestroyPSDVector ( (*vect) );
  (*vect) = NULL;

  RETURN (status);

} /* LALDestroyPSDVector() */

/**
 * \deprecated Use XLALDestroyMultiSFTVector() instead.
 */
void
LALDestroyMultiSFTVector (LALStatus *status,		/**< pointer to LALStatus structure */
		          MultiSFTVector **multvect)	/**< the SFT-vector to free */
{
  INITSTATUS(status);
  ASSERT (multvect != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  XLALDestroyMultiSFTVector ( (*multvect) );

  (*multvect) = NULL;

  RETURN (status);

} /* LALDestroyMultiSFTVector() */

/**
 * \deprecated Use XLALDestroyMultiPSDVector() instead.
 */
void
LALDestroyMultiPSDVector (LALStatus *status,		/**< pointer to LALStatus structure */
		          MultiPSDVector **multvect)	/**< the SFT-vector to free */
{
  INITSTATUS(status);
  ASSERT (multvect != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  XLALDestroyMultiPSDVector ( (*multvect) );
  (*multvect) = NULL;

  RETURN (status);

} /* LALDestroyMultiPSDVector() */

/**
 * \deprecated Use XLALCopySFT()
 *
 * Copy an entire SFT-type into another.
 * We require the destination-SFT to have a NULL data-entry, as the
 * corresponding data-vector will be allocated here and copied into
 *
 * Note: the source-SFT is allowed to have a NULL data-entry,
 * in which case only the header is copied.
 */
void
LALCopySFT (LALStatus *status,	/**< pointer to LALStatus structure */
	    SFTtype *dest, 	/**< [out] copied SFT (needs to be allocated already) */
	    const SFTtype *src)	/**< input-SFT to be copied */
{
  INITSTATUS(status);

  XLAL_PRINT_DEPRECATION_WARNING("XLALCopySFT");

  int ret = XLALCopySFT ( dest, src );
  if ( ret != XLAL_SUCCESS ) {
    ABORT ( status, SFTUTILS_EFUNC, SFTUTILS_MSGEFUNC );
  }

  RETURN (status);

} /* LALCopySFT() */

/** Append the given SFTtype to the SFT-vector (no SFT-specific checks are done!) */
void
LALAppendSFT2Vector (LALStatus *status,		/**< pointer to LALStatus structure */
		     SFTVector *vect,		/**< destinatino SFTVector to append to */
		     const SFTtype *sft)	/**< the SFT to append */
{
  UINT4 oldlen;
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( sft, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL );
  ASSERT ( vect, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL );

  oldlen = vect->length;

  if ( (vect->data = LALRealloc ( vect->data, (oldlen + 1)*sizeof( *vect->data ) )) == NULL ) {
    ABORT ( status, SFTUTILS_EMEM, SFTUTILS_MSGEMEM );
  }
  memset ( &(vect->data[oldlen]), 0, sizeof( vect->data[0] ) );
  vect->length ++;

  TRY ( LALCopySFT( status->statusPtr, &vect->data[oldlen], sft ), status);

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* LALAppendSFT2Vector() */


/** \deprecated LAL-interface: use XLALCreateTimestampVector() instead */
void
LALCreateTimestampVector (LALStatus *status,		/**< pointer to LALStatus structure */
			  LIGOTimeGPSVector **vect, 	/**< [out] allocated timestamp-vector  */
			  UINT4 length)			/**< number of elements */
{
  INITSTATUS(status);

  ASSERT (vect != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT (*vect == NULL, status, SFTUTILS_ENONULL,  SFTUTILS_MSGENONULL);

  LIGOTimeGPSVector *out = NULL;
  if ( (out = XLALCreateTimestampVector( length )) == NULL ) {
    XLALPrintError ("XLALCreateTimestampVector() failed with xlalErrno = %d\n", xlalErrno );
    ABORT (status,  SFTUTILS_EMEM,  SFTUTILS_MSGEMEM);
  }

  (*vect) = out;

  RETURN (status);

} /* LALCreateTimestampVector() */

/**
 * \deprecated Use XLALDestroyTimestampVector() instead.
 */
void
LALDestroyTimestampVector (LALStatus *status,		/**< pointer to LALStatus structure */
			   LIGOTimeGPSVector **vect)	/**< timestamps-vector to be freed */
{
  INITSTATUS(status);

  ASSERT (vect != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  if ( (*vect) == NULL )
    goto finished;

  XLALDestroyTimestampVector ( (*vect) );

  (*vect) = NULL;

 finished:
  RETURN (status);

} /* LALDestroyTimestampVector() */

/**
 * \deprecated Use XLALMakeTimestamps() instead.
 */
void
LALMakeTimestamps ( LALStatus *status,			/**< pointer to LALStatus structure */
                    LIGOTimeGPSVector **timestamps, 	/**< [out] timestamps-vector */
                    LIGOTimeGPS tStart,			/**< GPS start-time */
                    REAL8 duration, 			/**< duration in seconds */
                    REAL8 tStep				/**< length of one (SFT) timestretch in seconds */
                    )
{
  INITSTATUS(status);

  ASSERT (timestamps != NULL, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);
  ASSERT (*timestamps == NULL,status, SFTUTILS_ENONULL, SFTUTILS_MSGENONULL);

  LIGOTimeGPSVector *ts;
  ts = XLALMakeTimestamps ( tStart, duration, tStep, 0 );
  if ( ts == NULL )
    {
      XLALPrintError ("XLALMakeTimestamps() failed with xlalErrno = %d\n", xlalErrno );
      ABORT ( status,  SFTUTILS_EFUNC,  SFTUTILS_MSGEFUNC );
    }

  (*timestamps) = ts;

  RETURN( status );

} /* LALMakeTimestamps() */

/**
 * \deprecated LAL wrapper to XLALExtractTimestampsFromSFTs()
 */
void
LALGetSFTtimestamps (LALStatus *status,			/**< pointer to LALStatus structure */
		     LIGOTimeGPSVector **timestamps,	/**< [out] extracted timestamps */
		     const SFTVector *sfts )		/**< input SFT-vector  */
{
  LIGOTimeGPSVector *ret = NULL;

  INITSTATUS(status);

  ASSERT ( timestamps, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL );
  ASSERT ( sfts, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL );
  ASSERT ( sfts->length > 0, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL );
  ASSERT ( *timestamps == NULL, status, SFTUTILS_ENONULL, SFTUTILS_MSGENONULL );

  if ( ( ret = XLALExtractTimestampsFromSFTs ( sfts )) == NULL ) {
    XLALPrintError ("%s: call to XLALExtractTimestampsFromSFTs() failed with code %d\n", __func__, xlalErrno );
    ABORT (status, SFTUTILS_EFUNC, SFTUTILS_MSGEFUNC);
  }

  /* done: return Ts-vector */
  (*timestamps) = ret;

  RETURN(status);

} /* LALGetSFTtimestamps() */

/**
 * Computes weight factors arising from SFTs with different noise
 * floors -- it multiplies an existing weight vector
 */
void
LALComputeNoiseWeights  (LALStatus        *status,
                          REAL8Vector      *weightV,
                          const SFTVector  *sftVect,
                          INT4             blkSize,
                          UINT4            excludePercentile)
{

  UINT4 lengthVect, lengthSFT, lengthPSD, halfLengthPSD;
  UINT4 j, excludeIndex;
  SFTtype *sft;
  REAL8FrequencySeries periodo;
  REAL8Sequence mediansV, inputV;
  LALRunningMedianPar rngMedPar;

  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (weightV, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);
  ASSERT (sftVect, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);
  ASSERT (blkSize > 0, status,  SFTUTILS_EINPUT, SFTUTILS_MSGEINPUT);
  ASSERT (weightV->data,status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);
  ASSERT (sftVect->data,status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);
  ASSERT (excludePercentile <= 100, status, SFTUTILS_EINPUT, SFTUTILS_MSGEINPUT);
  /* -------------------------------------------   */

  /* Make sure there is no size mismatch */
  ASSERT (weightV->length == sftVect->length, status, SFTUTILS_EINPUT, SFTUTILS_MSGEINPUT);
  /* -------------------------------------------   */

  /* Make sure there are elements to be computed*/
  ASSERT (sftVect->length, status, SFTUTILS_EINPUT, SFTUTILS_MSGEINPUT);


  /* set various lengths */
  lengthVect = sftVect->length;
  lengthSFT = sftVect->data->data->length;
  ASSERT( lengthSFT > 0, status,  SFTUTILS_EINPUT, SFTUTILS_MSGEINPUT);
  lengthPSD = lengthSFT - blkSize + 1;

  /* make sure blksize is not too big */
  ASSERT(lengthPSD > 0, status, SFTUTILS_EINPUT, SFTUTILS_MSGEINPUT);

  halfLengthPSD = lengthPSD/2; /* integer division */

  /* allocate memory for periodogram */
  periodo.data = NULL;
  periodo.data = (REAL8Sequence *)LALMalloc(sizeof(REAL8Sequence));
  periodo.data->length = lengthSFT;
  periodo.data->data = (REAL8 *)LALMalloc( lengthSFT * sizeof(REAL8));

  /* allocate memory for vector of medians */
  mediansV.length = lengthPSD;
  mediansV.data = (REAL8 *)LALMalloc(lengthPSD * sizeof(REAL8));

  /* rng med block size */
  rngMedPar.blocksize = blkSize;

  /* calculate index in psd medians vector from which to calculate mean */
  excludeIndex =  (excludePercentile * halfLengthPSD) ; /* integer arithmetic */
  excludeIndex /= 100; /* integer arithmetic */

  /* loop over sfts and calculate weights */
  for (j=0; j<lengthVect; j++) {
    REAL8 sumMed = 0.0;
    UINT4 k;

    sft = sftVect->data + j;

    /* calculate the periodogram */
    TRY (LALSFTtoPeriodogram (status->statusPtr, &periodo, sft), status);

    /* calculate the running median */
    inputV.length = lengthSFT;
    inputV.data = periodo.data->data;
    TRY( LALDRunningMedian2(status->statusPtr, &mediansV, &inputV, rngMedPar), status);

    /* now sort the mediansV.data vector and exclude the top and last percentiles */
    gsl_sort(mediansV.data, 1, mediansV.length);

    /* sum median excluding appropriate elements */
    for (k = excludeIndex; k < lengthPSD - excludeIndex; k++) {
      sumMed += mediansV.data[k];
    }

    /* weight is proportional to 1/sumMed */
    weightV->data[j] /= sumMed;

  } /* end of loop over sfts */

  /* remember to normalize weights immediately after leaving this function */

  /* free memory */
  LALFree(mediansV.data);
  LALFree(periodo.data->data);
  LALFree(periodo.data);

  DETATCHSTATUSPTR (status);
   /* normal exit */
  RETURN (status);

} /* LALComputeNoiseWeights() */


/**
 * \deprecated Use XLALComputeMultiNoiseWeights() instead
 */
void LALComputeMultiNoiseWeights  (LALStatus             *status,
				   MultiNoiseWeights     **out,
				   const MultiPSDVector  *rngmed,
				   UINT4                 blocksRngMed,
				   UINT4                 excludePercentile)
{

  MultiNoiseWeights *weights = NULL;

  INITSTATUS(status);

  ASSERT ( rngmed, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);
  ASSERT ( rngmed->data, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);
  ASSERT ( rngmed->length, status, SFTUTILS_EINPUT, SFTUTILS_MSGEINPUT);

  ASSERT ( out, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);
  ASSERT ( *out == NULL, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);

  if ( (weights = XLALComputeMultiNoiseWeights ( rngmed, blocksRngMed, excludePercentile )) == NULL )
    {
      XLALPrintError ("XLALComputeMultiNoiseWeights() failed with xlalErrno = %d\n", xlalErrno );
      ABORT ( status, SFTUTILS_EFUNC, SFTUTILS_MSGEFUNC );
    }

  *out = weights;

  RETURN (status);

} /* LALComputeMultiNoiseWeights() */


/** \deprecated Use XLALDestroyMultiNoiseWeights() instead */
void
LALDestroyMultiNoiseWeights  (LALStatus         *status,
			      MultiNoiseWeights **weights)
{
  INITSTATUS(status);

  ASSERT ( weights != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  XLALDestroyMultiNoiseWeights ( (*weights) );

  (*weights) = NULL;

  RETURN (status);

} /* LALDestroyMultiNoiseWeights() */


/* ==================================================
 * SFT up-sampling routines
 * ==================================================
 */

/**
 * upsample a given multi-SFTvector by the given (integer) factor,
 * _replacing_ the original SFTs
 */
void
upsampleMultiSFTVector (LALStatus *status,		/**< pointer to LALStatus structure */
			  MultiSFTVector *inout,	/**< [in,out]: upsampled multi SFT-vector */
			  UINT4 upsample, 		/**< integer factor to upsample by */
			  UINT4 Dterms			/**< number of terms in Dirichlet kernel [on each side] */
			  )
{
  UINT4 X, numDet;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( inout, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);
  ASSERT ( inout->length, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);

  if ( upsample < 2 ) 	/* nothing to do */
    goto done;

  numDet = inout->length;

  for ( X=0; X < numDet; X ++ )
    {
      SFTVector *thisSFTvect = inout->data[X];
      TRY ( upsampleSFTVector ( status->statusPtr, thisSFTvect, upsample, Dterms ), status );
    } /* for X < numDet */

 done:
  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* upsampleMultiSFTVector() */


void
upsampleSFTVector (LALStatus *status,		/**< pointer to LALStatus structure */
		     SFTVector *inout,		/**< [in,out]: upsampled SFT-vector */
		     UINT4 upsample, 		/**< integer factor to upsample by */
		     UINT4 Dterms		/**< number of terms in Dirichlet kernel [on each side] */
		     )
{
  UINT4 alpha, numSFTs;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( inout, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);
  ASSERT ( inout->length, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);

  numSFTs = inout->length;

  for ( alpha=0; alpha < numSFTs; alpha ++ )
    {
      COMPLEX8Vector *this_data = inout->data[alpha].data;
      COMPLEX8Vector *new_data;
      if ( (new_data = XLALrefineCOMPLEX8Vector ( this_data, upsample, Dterms )) == NULL ) {
	XLALPrintError ("\nSFT oversampling failed ... \n\n");
	ABORT ( status, SFTUTILS_EFUNC,SFTUTILS_MSGEFUNC );
      }

      /* now replace old SFT with new upsampled one */
      XLALDestroyCOMPLEX8Vector ( this_data );
      inout->data[alpha].data = new_data;
      /*
      inout->data[alpha].deltaF /= 1.0 * upsample;
      */

    } /* for alpha < numSFTs */



  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* upsampleSFTVector() */

/** \endcond */
