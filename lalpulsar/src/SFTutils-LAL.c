/** \cond DONT_DOXYGEN */

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
