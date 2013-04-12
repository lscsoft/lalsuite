/** \deprecated Use XLALCreateSFT() instead
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

/** \deprecated Use XLALCreateSFTVector() instead
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

/** Create an empty multi-IFO SFT vector for given number of IFOs and number of SFTs per IFO
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

/** \deprecated Use XLALDestroySFT() instead.
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

/** \deprecated Use XLALDestroySFTVector() instead.
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

/** \deprecated Use XLALDestroyPSDVector() instead
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

/** \deprecated Use XLALDestroyMultiSFTVector() instead.
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

/** \deprecate Use XLALDestroyMultiPSDVector() instead.
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

/** Copy an entire SFT-type into another.
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
  ATTATCHSTATUSPTR ( status );

  ASSERT (dest,  status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT (dest->data == NULL, status, SFTUTILS_ENONULL, SFTUTILS_MSGENONULL );
  ASSERT (src, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  /* copy complete head (including data-pointer, but this will be separately alloc'ed and copied in the next step) */
  memcpy ( dest, src, sizeof(*dest) );

  /* copy data (if there's any )*/
  if ( src->data )
    {
      UINT4 numBins = src->data->length;
      if ( (dest->data = XLALCreateCOMPLEX8Vector ( numBins )) == NULL ) {
	ABORT ( status, SFTUTILS_EMEM, SFTUTILS_MSGEMEM );
      }
      memcpy (dest->data->data, src->data->data, numBins * sizeof (src->data->data[0]));
    }

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALCopySFT() */



/** Subtract two SFT-vectors and put the results in a new one (which it allocates).
 *
 */
void
LALSubtractSFTVectors (LALStatus *status,	/**< pointer to LALStatus structure */
		     SFTVector **outVect,	/**< [out] difference of SFT-vectors */
		     const SFTVector *inVect1,	/**< input-vector 1 */
		     const SFTVector *inVect2 ) /**< input-vector 2 */
{
  UINT4 numSFTs1, numSFTs2;
  UINT4 i, j;
  SFTVector *ret = NULL;
  CHAR name1Trunc[LALNameLength];
  CHAR name2Trunc[LALNameLength];
  CHAR prefix[LALNameLength];
  UINT4 halfNameLength;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (outVect,  status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT ( *outVect == NULL,  status, SFTUTILS_ENONULL,  SFTUTILS_MSGENONULL);
  ASSERT (inVect1 && inVect1->data, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT (inVect2 && inVect2->data, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT ( inVect1->data[0].data && inVect2->data[0].data, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  numSFTs1 = inVect1 -> length;
  numSFTs2 = inVect2 -> length;

  if ( numSFTs1 != numSFTs2 )
    {
      XLALPrintError ("\nERROR: the SFT-vectors must have the same number of SFTs!\n\n");
      ABORT ( status, SFTUTILS_EINPUT,  SFTUTILS_MSGEINPUT);
    }

  TRY ( LALCreateSFTVector ( status->statusPtr, &ret, numSFTs1, inVect1->data[0].data->length ), status );

  halfNameLength = (LALNameLength - strlen("Xn:{}-{}"))/2;

  /* copy the SFTs and subtract their data one-by-one */
  for (i=0; i < numSFTs1; i ++)
    {
      UINT4 numBins1, numBins2;
      LIGOTimeGPS epoch1, epoch2;
      REAL8 Freq1, Freq2, deltaF1, deltaF2;
      numBins1 = inVect1->data[i].data->length;
      numBins2 = inVect2->data[i].data->length;
      epoch1   = inVect1->data[i].epoch;
      epoch2   = inVect2->data[i].epoch;
      Freq1    = inVect1->data[i].f0;
      Freq2    = inVect2->data[i].f0;
      deltaF1  = inVect1->data[i].deltaF;
      deltaF2  = inVect2->data[i].deltaF;

      if ( numBins1 != numBins2 ) {
	XLALPrintError ("\nERROR: the SFTs must have the same number of frequency-bins!\n\n");
	goto failed;
      }
      if ( (epoch1.gpsSeconds != epoch2.gpsSeconds) || ( epoch1.gpsNanoSeconds != epoch2.gpsNanoSeconds ) ) {
	XLALPrintError ("\nERROR: the SFTs must have the same epochs!\n\n");
	goto failed;
      }
      if ( Freq1 != Freq2 ) {
	XLALPrintError ("\nERROR: the SFTs must have the same start frequency!\n\n");
	goto failed;
      }
      if ( deltaF1 != deltaF2 ) {
	XLALPrintError ("\nERROR: the SFTs must have the same frequency-steps!\n\n");
	goto failed;
      }
      /* copy header info */
      ret->data[i].epoch  = epoch1;
      ret->data[i].f0     = Freq1;
      ret->data[i].deltaF = deltaF1;

      for (j=0; j < numBins1; j++)
	{
	  ret->data[i].data->data[j].realf_FIXME = crealf(inVect1->data[i].data->data[j]) - crealf(inVect2->data[i].data->data[j]);
	  ret->data[i].data->data[j].imagf_FIXME = cimagf(inVect1->data[i].data->data[j]) - cimagf(inVect2->data[i].data->data[j]);
	}  /* for j < numBins1 */

      snprintf ( name1Trunc, halfNameLength, "%s", inVect1->data[i].name );
      snprintf ( name2Trunc, halfNameLength, "%s", inVect2->data[i].name );
      snprintf ( prefix, (strlen("Xn:") + 1), "%s", inVect1->data[i].name );
      snprintf ( ret->data[i].name, LALNameLength, "%s{%s}-{%s}", prefix, name1Trunc, name2Trunc );
    } /* for i < numSFTs1 */

  /* success: */
  (*outVect) = ret;
  DETATCHSTATUSPTR (status);
  RETURN (status);

 failed:
  LALDestroySFTVector (  status->statusPtr, &ret );
  ABORT ( status, SFTUTILS_EINPUT,  SFTUTILS_MSGEINPUT);

} /* LALSubtractSFTVectors() */



/** Linearly combine two or more SFT-vectors and put the results in a new one (which it allocates).
 *
 */
void
LALLinearlyCombineSFTVectors
(LALStatus *status,		/**< pointer to LALStatus structure */
 SFTVector **outVect,	          /**< [out] linear combo of SFT-vectors */
 SFTVector **inVects,	  /**< array of SFT-vectors */
 const COMPLEX16Vector *weights,  /**< vector of SFT-weights */
 const CHAR *outName)             /**< name for output vector */
{
  UINT4 numSFTs, numSFTVects;
  UINT4 i, j, k;
  SFTVector *ret = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (outVect,  status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT ( *outVect == NULL,  status, SFTUTILS_ENONULL,  SFTUTILS_MSGENONULL);
  ASSERT (inVects && inVects[0] && inVects[0]->data
	  && inVects[0]->data[0].data,
	  status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);
  ASSERT (weights && weights->data,
	  status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  numSFTVects = weights->length;

  if ( numSFTVects < 1 )
    {
      XLALPrintError ("\nERROR: must be combining at least one SFT Vector!\n\n");
      ABORT ( status, SFTUTILS_EINPUT,  SFTUTILS_MSGEINPUT);
    }

  numSFTs = inVects[0] -> length;

  TRY ( LALCreateSFTVector ( status->statusPtr, &ret, numSFTs, inVects[0]->data[0].data->length ), status );

  /* copy the SFTs from the first vector */
  for (i=0; i < numSFTs; i ++)
    {
      UINT4 numBins1, numBins2;
      LIGOTimeGPS epoch1, epoch2;
      REAL8 Freq1, Freq2, deltaF1, deltaF2;
      numBins1 = inVects[0]->data[i].data->length;
      epoch1   = inVects[0]->data[i].epoch;
      Freq1    = inVects[0]->data[i].f0;
      deltaF1  = inVects[0]->data[i].deltaF;

      /* copy header info */
      ret->data[i].epoch  = epoch1;
      ret->data[i].f0     = Freq1;
      ret->data[i].deltaF = deltaF1;

      for (k=0; k < numBins1; k++)
	{
	  ret->data[i].data->data[k].realf_FIXME
	    = weights->data[0].re * crealf(inVects[0]->data[i].data->data[k])
	    - weights->data[0].im * cimagf(inVects[0]->data[i].data->data[k]);
	  ret->data[i].data->data[k].imagf_FIXME
	    = weights->data[0].re * cimagf(inVects[0]->data[i].data->data[k])
	    + weights->data[0].im * crealf(inVects[0]->data[i].data->data[k]);
	}  /* for k < numBins1 */

      /* add in the other SFTs one-by-one */
      for (j=1; j < numSFTVects; j++)
	{
	  numBins2 = inVects[j]->data[i].data->length;
	  epoch2   = inVects[j]->data[i].epoch;
	  Freq2    = inVects[j]->data[i].f0;
	  deltaF2  = inVects[j]->data[i].deltaF;

	  if ( numBins1 != numBins2 ) {
	    XLALPrintError ("\nERROR: the SFTs must have the same number of frequency-bins!\n\n");
	    goto failed;
	  }
	  if ( (epoch1.gpsSeconds != epoch2.gpsSeconds) || ( epoch1.gpsNanoSeconds != epoch2.gpsNanoSeconds ) ) {
	    XLALPrintError ("\nERROR: the SFTs must have the same epochs!\n\n");
	    goto failed;
	  }
	  if ( Freq1 != Freq2 ) {
	    XLALPrintError ("\nERROR: the SFTs must have the same start frequency!\n\n");
	    goto failed;
	  }
	  if ( deltaF1 != deltaF2 ) {
	    XLALPrintError ("\nERROR: the SFTs must have the same frequency-steps!\n\n");
	    goto failed;
	  }


	  for (k=0; k < numBins1; k++)
	    {
	      ret->data[i].data->data[k].realf_FIXME
		+= weights->data[j].re * crealf(inVects[j]->data[i].data->data[k])
		- weights->data[j].im * cimagf(inVects[j]->data[i].data->data[k]);
	      ret->data[i].data->data[k].imagf_FIXME
		+= weights->data[j].re * cimagf(inVects[j]->data[i].data->data[k])
		+ weights->data[j].im * crealf(inVects[j]->data[i].data->data[k]);
	    }  /* for k < numBins1 */

	} /* for j < numSFTVects */
      memcpy ( ret->data[i].name, outName, LALNameLength*sizeof(CHAR) );
    } /* for i < numSFTs */

  /* success: */
  (*outVect) = ret;
  DETATCHSTATUSPTR (status);
  RETURN (status);

 failed:
  LALDestroySFTVector (  status->statusPtr, &ret );
  ABORT ( status, SFTUTILS_EINPUT,  SFTUTILS_MSGEINPUT);

} /* LALLinearlyCombineSFTVectors() */



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

/** \deprecated Use XLALDestroyTimestampVector() instead.
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

/** \deprecated Use XLALMakeTimestamps() instead.
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
  ts = XLALMakeTimestamps ( tStart, duration, tStep );
  if ( ts == NULL )
    {
      XLALPrintError ("XLALMakeTimestamps() failed with xlalErrno = %d\n", xlalErrno );
      ABORT ( status,  SFTUTILS_EFUNC,  SFTUTILS_MSGEFUNC );
    }

  (*timestamps) = ts;

  RETURN( status );

} /* LALMakeTimestamps() */

/** \deprecated LAL wrapper to XLALExtractTimestampsFromSFTs()
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

/** Computes weight factors arising from SFTs with different noise
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


/** Computes weight factors arising from MultiSFTs with different noise
 * floors
 */
void LALComputeMultiNoiseWeights  (LALStatus             *status,
				   MultiNoiseWeights     **out,
				   const MultiPSDVector  *rngmed,
				   UINT4                 blocksRngMed,
				   UINT4                 excludePercentile)
{
  UINT4 Y, X, alpha, k, numifos, numsfts, lengthsft, numsftsTot;
  MultiNoiseWeights *weights;
  REAL8 Tsft = 1.0 / rngmed->data[0]->data[0].deltaF;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( rngmed, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);
  ASSERT ( rngmed->data, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);
  ASSERT ( rngmed->length, status, SFTUTILS_EINPUT, SFTUTILS_MSGEINPUT);

  ASSERT ( out, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);
  ASSERT ( *out == NULL, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);

  numifos = rngmed->length;

  if ( (weights = (MultiNoiseWeights *)LALCalloc(1, sizeof(MultiNoiseWeights))) == NULL ){
    ABORT (status,  SFTUTILS_EMEM,  SFTUTILS_MSGEMEM);
  }

  weights->length = numifos;
  if ( (weights->data = (REAL8Vector **)LALCalloc( numifos, sizeof(REAL8Vector *))) == NULL) {
    ABORT (status,  SFTUTILS_EMEM,  SFTUTILS_MSGEMEM);
  }

  numsftsTot = 0;
  REAL8 sumWeights = 0;

  for ( X = 0; X < numifos; X++)
    {
      numsfts = rngmed->data[X]->length;
      numsftsTot += numsfts;

      /* create k^th weights vector */
      LALDCreateVector ( status->statusPtr, &(weights->data[X]), numsfts);
      BEGINFAIL( status ) {
	for ( Y = 0; Y < X-1; Y++)
	  LALDDestroyVector (status->statusPtr, &(weights->data[Y]));
	LALFree (weights->data);
	LALFree (weights);
      } ENDFAIL(status);

      /* loop over rngmeds and calculate weights -- one for each sft */
      for ( alpha = 0; alpha < numsfts; alpha++)
	{
	  REAL8FrequencySeries *thisrm;
	  UINT4 halfBlock = blocksRngMed/2;
	  UINT4 excludeIndex, halfLength, length;
          REAL8 wXa;

	  thisrm = &(rngmed->data[X]->data[alpha]);

	  lengthsft = thisrm->data->length;
	  if ( lengthsft < blocksRngMed ) {
	    ABORT ( status, SFTUTILS_EINPUT, SFTUTILS_MSGEINPUT);
	  }

	  length = lengthsft - blocksRngMed + 1;
	  halfLength = length/2;

	  /* calculate index in power medians vector from which to calculate mean */
	  excludeIndex =  excludePercentile * halfLength ; /* integer arithmetic */
	  excludeIndex /= 100; /* integer arithmetic */

	  REAL8 Tsft_avgS2 = 0.0;	// 'S2' refers to double-sided PSD
	  for ( k = halfBlock + excludeIndex; k < lengthsft - halfBlock - excludeIndex; k++)
	    Tsft_avgS2 += thisrm->data->data[k];
	  Tsft_avgS2 /= lengthsft - 2*halfBlock - 2*excludeIndex;

          wXa = 1.0/Tsft_avgS2;	// unnormalized weight
	  weights->data[X]->data[alpha] = wXa;

	  sumWeights += wXa;	// sum the weights to normalize this at the end
	} /* end loop over sfts for each ifo */

    } /* end loop over ifos */

  /* overall noise-normalization factor Sinv = 1/Nsft sum_Xa Sinv_Xa,
   * see Eq.(60) in CFSv2 notes:
   * https://dcc.ligo.org/cgi-bin/private/DocDB/ShowDocument?docid=1665&version=3
   */
  REAL8 TsftS2_inv = sumWeights / numsftsTot;	// this is double-sided PSD 'S2'

  /* make weights of order unity by normalizing with TsftS2_inv, see Eq.(58) in CFSv2 notes (v3) */
  for ( X = 0; X < numifos; X ++) {
    numsfts = weights->data[X]->length;
    for ( alpha = 0; alpha < numsfts; alpha ++)
      weights->data[X]->data[alpha] /= TsftS2_inv;
  }

  weights->Sinv_Tsft = 0.5 * Tsft*Tsft * TsftS2_inv;		/* 'Sinv * Tsft' refers to single-sided PSD!! Eq.(60) in CFSv2 notes (v3)*/

  *out = weights;


  DETATCHSTATUSPTR (status);
   /* normal exit */
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

/** upsample a given multi-SFTvector by the given (integer) factor,
 *  _replacing_ the original SFTs
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
