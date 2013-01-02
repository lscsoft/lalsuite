/*
 * Copyright (C) 2005 Reinhard Prix
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

/*---------- INCLUDES ----------*/
#include <stdarg.h>

#include <gsl/gsl_sort_double.h>

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/FrequencySeries.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/LISAspecifics.h>

#include "SFTutils.h"

/*---------- DEFINES ----------*/
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/
/* empty struct initializers */
const SFTtype empty_SFTtype;
const SFTVector empty_SFTVector;
const PSDVector empty_PSDVector;
const MultiSFTVector empty_MultiSFTVector;
const MultiPSDVector empty_MultiPSDVector;
const MultiNoiseWeights empty_MultiNoiseWeights;
const MultiREAL4TimeSeries empty_MultiREAL4TimeSeries;
const LIGOTimeGPSVector empty_LIGOTimeGPSVector;
const MultiLIGOTimeGPSVector empty_MultiLIGOTimeGPSVector;

/*---------- internal prototypes ----------*/

/*==================== FUNCTION DEFINITIONS ====================*/

/** Extract a frequency band from an SFTVector, returning a new SFTvector
 * Note: fMin < 0 implies to start from lowest frequency bin,
 *       fMax < 0 implies to include up to highest frequency bin
 *
 * if fMin, fMax > 0, the corresponding frequency MUST be contained in the
 * input SFT, otherwise an error is returned.
 *
 * We guarantee that both fMin and fMax will be *contained* in the returned SFT,
 * which means the actual min(f) can be < fMin, and max(f) > fMax is possible.
 *
 */
SFTVector *
XLALExtractBandfromSFTs ( const SFTVector *sfts, REAL8 fMin, REAL8 fMax )
{
  REAL8 dFreq;
  UINT4 iMin, iMax, i0, numSFTs, numBinsIn, numBinsOut, iSFT;
  SFTVector *out;
  REAL8 f0Out;
  COMPLEX8Vector *sav;

  if ( !sfts || !sfts->data || (sfts->length==0) || !sfts->data[0].data ) {
    XLAL_ERROR_NULL( XLAL_EINVAL );
  }

  dFreq = sfts->data[0].deltaF;

  i0 = floor ( sfts->data[0].f0 / dFreq + 0.5 );	/* round to nearest bin */
  numBinsIn = sfts->data[0].data->length;

  if ( fMin < 0 )
    iMin = i0;
  else
    iMin = floor ( fMin / dFreq + 1e-6 );	/* round down */

  if ( fMax < 0 )
    iMax = i0 + numBinsIn - 1;
  else
    iMax = ceil ( fMax / dFreq - 1e-6 );	/* round up */

  if ( iMax < iMin ) {
    XLALPrintError ("Resulting SFT has no bins iMax (%d) < iMin (%d)!\n", iMax, iMin );
    XLAL_ERROR_NULL( XLAL_EINVAL );
  }

  numSFTs = sfts->length;
  f0Out = iMin * dFreq;
  numBinsOut = iMax - iMin + 1;

  if ( (out = XLALCreateSFTVector ( numSFTs, numBinsOut )) == NULL ) {
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  /* now copy heads and all requested bins */
  for ( iSFT = 0; iSFT < numSFTs; iSFT ++ )
    {
      /* first copy complete head (saving data-pointer first, which will be copied in the next step) */
      sav = out->data[iSFT].data;
      memcpy ( &(out->data[iSFT]), &(sfts->data[iSFT]), sizeof(sfts->data[0]) );
      out->data[iSFT].data = sav;

      /* fix new header information */
      out->data[iSFT].f0 = f0Out;

      /* copy data */
      memcpy (out->data[iSFT].data->data, sfts->data[iSFT].data->data + iMin - i0, numBinsOut * sizeof (sfts->data[iSFT].data->data[0] ) );
      out->data[iSFT].data->length = numBinsOut;

    } /* for iSFT < numSFTs */

  return ( out );

} /* XLALExtractBandfromSFTs() */


/** XLAL function to create one SFT-struct.
 *
 * Note: Allows for numBins == 0, in which case only the header is
 * allocated, with a NULL data pointer.
 */
SFTtype *
XLALCreateSFT ( UINT4 numBins )
{
  SFTtype *sft;

  if ( (sft = XLALCalloc (1, sizeof(*sft) )) == NULL )
    XLAL_ERROR_NULL ( XLAL_ENOMEM, "XLALCalloc (1, %d) failed.\n", sizeof(*sft) );

  if ( numBins )
    {
      if ( (sft->data = XLALCreateCOMPLEX8Vector ( numBins )) == NULL )
        {
          XLALFree ( sft );
          XLAL_ERROR_NULL ( XLAL_EFUNC, "XLALCreateCOMPLEX8Vector ( %s ) failed. xlalErrno = %d\n", numBins, xlalErrno );
        }
    }
  else
    sft->data = NULL;	/* no data, just header */

  return sft;

} /* XLALCreateSFT() */

/** Destructor for one SFT */
void
XLALDestroySFT ( SFTtype *sft )
{
  if ( !sft )
    return;

  if ( sft->data )
    XLALDestroyCOMPLEX8Vector ( sft->data );

  XLALFree ( sft );

  return;

} /* XLALDestroySFT() */



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


/** XLAL function to create an SFTVector of \c numSFT SFTs with \c SFTlen frequency-bins
 */
SFTVector *
XLALCreateSFTVector (UINT4 numSFTs, 	/**< number of SFTs */
		     UINT4 numBins	/**< number of frequency-bins per SFT */
		     )
{
  UINT4 iSFT;
  SFTVector *vect;

  if ( (vect = XLALCalloc ( 1, sizeof(*vect) )) == NULL ) {
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  vect->length = numSFTs;
  if ( (vect->data = XLALCalloc (1, numSFTs * sizeof ( *vect->data ) )) == NULL ) {
    XLALFree (vect);
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  for ( iSFT = 0; iSFT < numSFTs; iSFT ++)
    {
      COMPLEX8Vector *data = NULL;

      /* allow SFTs with 0 bins: only header */
      if ( numBins )
	{
	  if ( (data = XLALCreateCOMPLEX8Vector ( numBins )) == NULL )
	    {
	      UINT4 j;
	      for ( j = 0; j < iSFT; j++ )
		XLALDestroyCOMPLEX8Vector ( vect->data[j].data );
	      XLALFree (vect->data);
	      XLALFree (vect);
	      XLAL_ERROR_NULL( XLAL_ENOMEM );
	    }
	}

      vect->data[iSFT].data = data;

    } /* for iSFT < numSFTs */

  return vect;

} /* XLALCreateSFTVector() */



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


/** XLAL interface to destroy an SFTVector
 */
void
XLALDestroySFTVector ( SFTVector *vect )
{
  if ( !vect )
    return;

  for ( UINT4 i=0; i < vect->length; i++ )
    {
      SFTtype *sft = &( vect->data[i] );
      if ( sft->data )
	{
	  if ( sft->data->data )
	    XLALFree ( sft->data->data );
	  XLALFree ( sft->data );
	}
    } // for i < numSFTs

  XLALFree ( vect->data );
  XLALFree ( vect );

  return;

} /* XLALDestroySFTVector() */


/** Destroy a PSD-vector
 */
void
XLALDestroyPSDVector ( PSDVector *vect )	/**< the PSD-vector to free */
{
  if ( vect == NULL )	/* nothing to be done */
    return;

  for ( UINT4 i=0; i < vect->length; i++ )
    {
      REAL8FrequencySeries *psd = &( vect->data[i] );
      if ( psd->data )
	{
	  if ( psd->data->data )
	    XLALFree ( psd->data->data );
	  XLALFree ( psd->data );
	}
    } // for i < numPSDs

  XLALFree ( vect->data );
  XLALFree ( vect );

  return;

} /* XLALDestroyPSDVector() */

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


/** Destroy a multi SFT-vector
 */
void
XLALDestroyMultiSFTVector ( MultiSFTVector *multvect )	/**< the SFT-vector to free */
{
  if ( multvect == NULL )	/* nothing to be done */
    return;

  for ( UINT4 i = 0; i < multvect->length; i++ )
    XLALDestroySFTVector ( multvect->data[i] );

  XLALFree( multvect->data );
  XLALFree( multvect );

  return;

} /* XLALDestroyMultiSFTVector() */


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



/** Destroy a multi PSD-vector
 */
void
XLALDestroyMultiPSDVector ( MultiPSDVector *multvect )	/**< the SFT-vector to free */
{
  if ( multvect == NULL )
    return;

  for ( UINT4 i = 0; i < multvect->length; i++ )
    XLALDestroyPSDVector ( multvect->data[i] );

  XLALFree( multvect->data );
  XLALFree( multvect );

  return;

} /* XLALDestroyMultiPSDVector() */


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
	  ret->data[i].data->data[j].re = inVect1->data[i].data->data[j].re - inVect2->data[i].data->data[j].re;
	  ret->data[i].data->data[j].im = inVect1->data[i].data->data[j].im - inVect2->data[i].data->data[j].im;
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
	  ret->data[i].data->data[k].re
	    = weights->data[0].re * inVects[0]->data[i].data->data[k].re
	    - weights->data[0].im * inVects[0]->data[i].data->data[k].im;
	  ret->data[i].data->data[k].im
	    = weights->data[0].re * inVects[0]->data[i].data->data[k].im
	    + weights->data[0].im * inVects[0]->data[i].data->data[k].re;
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
	      ret->data[i].data->data[k].re
		+= weights->data[j].re * inVects[j]->data[i].data->data[k].re
		- weights->data[j].im * inVects[j]->data[i].data->data[k].im;
	      ret->data[i].data->data[k].im
		+= weights->data[j].re * inVects[j]->data[i].data->data[k].im
		+ weights->data[j].im * inVects[j]->data[i].data->data[k].re;
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



/** Allocate a LIGOTimeGPSVector */
LIGOTimeGPSVector *
XLALCreateTimestampVector ( UINT4 length )
{
  int len;
  LIGOTimeGPSVector *out = XLALCalloc (1, len = sizeof(LIGOTimeGPSVector));
  if (out == NULL)
    XLAL_ERROR_NULL ( XLAL_ENOMEM, "Failed to allocate XLALCalloc(1,%d)\n", len );

  out->length = length;
  out->data = XLALCalloc (1, len = length * sizeof(LIGOTimeGPS));
  if (out->data == NULL) {
    XLALFree (out);
    XLAL_ERROR_NULL ( XLAL_ENOMEM, "Failed to allocate XLALCalloc(1,%d)\n", len );
  }

  return out;

} /* XLALCreateTimestampVector() */


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


/** De-allocate a LIGOTimeGPSVector */
void
XLALDestroyTimestampVector ( LIGOTimeGPSVector *vect)
{
  if ( !vect )
    return;

  XLALFree ( vect->data );
  XLALFree ( vect );

  return;

} /* XLALDestroyTimestampVector() */


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



/** Given a start-time, duration and 'stepsize' tStep, returns a list of timestamps
 * covering this time-stretch.
 *
 * NOTE: boundary-handling: the returned list of timestamps are guaranteed to *cover* the
 * interval [tStart, tStart+duration], assuming a each timestamp covers a length of 'tStep'
 * This implies that the actual timestamps-coverage can extend up to 'tStep' beyond 'tStart+duration'.
 */
LIGOTimeGPSVector *
XLALMakeTimestamps ( LIGOTimeGPS tStart,		/**< GPS start-time */
                     REAL8 duration, 			/**< duration in seconds */
                     REAL8 tStep			/**< length of one (SFT) timestretch in seconds */
                     )
{
  XLAL_CHECK_NULL ( tStep > 0, XLAL_EDOM, "Invalid non-positive input 'tStart = %g'\n", tStep );
  XLAL_CHECK_NULL ( duration > 0, XLAL_EDOM, "Invalid non-positive input 'duration = %g'\n", duration );

  UINT4 numSFTs = ceil( duration / tStep );			/* >= 1 !*/

  LIGOTimeGPSVector *ts;
  XLAL_CHECK_NULL ( (ts = XLALCreateTimestampVector ( numSFTs )) != NULL, XLAL_EFUNC );

  ts->deltaT = tStep;

  LIGOTimeGPS tt = tStart;	/* initialize to start-time */
  for (UINT4 i = 0; i < numSFTs; i++)
    {
      ts->data[i] = tt;
      /* get next time-stamp */
      /* NOTE: we add the interval tStep successively (rounded correctly to ns each time!)
       * instead of using iSFT*Tsft, in order to avoid possible ns-rounding problems
       * with REAL8 intervals, which becomes critial from about 100days on...
       */
      XLAL_CHECK_NULL ( XLALGPSAdd ( &tt, tStep ) != NULL, XLAL_EFUNC );

    } /* for i < numSFTs */

  return ts;

} /* XLALMakeTimestamps() */


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



/** Extract timstamps-vector from the given SFTVector
 */
LIGOTimeGPSVector *
XLALExtractTimestampsFromSFTs ( const SFTVector *sfts )		/**< [in] input SFT-vector  */
{
  /* check input consistency */
  if ( !sfts ) {
    XLALPrintError ("%s: invalid NULL input 'sfts'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  UINT4 numSFTs = sfts->length;
  /* create output vector */
  LIGOTimeGPSVector *ret = NULL;
  if ( ( ret = XLALCreateTimestampVector ( numSFTs )) == NULL ) {
    XLALPrintError ("%s: XLALCreateTimestampVector(%d) failed.\n", __func__, numSFTs );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  REAL8 Tsft = 1.0 / sfts->data[0].deltaF;
  ret->deltaT = Tsft;

  UINT4 i;
  for ( i=0; i < numSFTs; i ++ )
    ret->data[i] = sfts->data[i].epoch;

  /* done: return Ts-vector */
  return ret;

} /* XLALExtractTimestampsFromSFTs() */


/** Given a multi-SFT vector, return a MultiLIGOTimeGPSVector holding the
 * SFT timestamps
 */
MultiLIGOTimeGPSVector *
XLALExtractMultiTimestampsFromSFTs ( const MultiSFTVector *multiSFTs )
{
  /* check input consistency */
  if ( !multiSFTs || multiSFTs->length == 0 ) {
    XLALPrintError ("%s: illegal NULL or empty input 'multiSFTs'.\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }
  UINT4 numIFOs = multiSFTs->length;

  /* create output vector */
  MultiLIGOTimeGPSVector *ret = NULL;
  if ( (ret = XLALCalloc ( 1, sizeof(*ret) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc ( 1, %d ).\n", __func__, sizeof(*ret));
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  if ( (ret->data = XLALCalloc ( numIFOs, sizeof(*ret->data) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc ( %d, %d ).\n", __func__, numIFOs, sizeof(ret->data[0]) );
    XLALFree (ret);
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  ret->length = numIFOs;

  /* now extract timestamps vector from each SFT-vector */
  UINT4 X;
  for ( X=0; X < numIFOs; X ++ )
    {
      if ( (ret->data[X] = XLALExtractTimestampsFromSFTs ( multiSFTs->data[X] )) == NULL ) {
        XLALPrintError ("%s: XLALExtractTimestampsFromSFTs() failed for X=%d\n", __func__, X );
        XLALDestroyMultiTimestamps ( ret );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }

    } /* for X < numIFOs */

  return ret;

} /* XLALExtractMultiTimestampsFromSFTs() */


/** Destroy a MultiLIGOTimeGPSVector timestamps vector
 */
void
XLALDestroyMultiTimestamps ( MultiLIGOTimeGPSVector *multiTS )
{
  UINT4 numIFOs, X;

  if ( !multiTS )
    return;

  numIFOs = multiTS->length;
  for ( X=0; X < numIFOs; X ++ ) {
    XLALDestroyTimestampVector ( multiTS->data[X] );
  }

  XLALFree ( multiTS->data );
  XLALFree ( multiTS );

  return;

} /* XLALDestroyMultiTimestamps() */




/** Extract/construct the unique 2-character "channel prefix" from the given
 * "detector-name", which unfortunately will not always follow any of the
 * official detector-naming conventions given in the Frames-Spec LIGO-T970130-F-E
 * This function therefore sometime has to do some creative guessing:
 *
 * NOTE: in case the channel-number can not be deduced from the name,
 * it is set to '1', and a warning will be printed if lalDebugLevel > 0.
 *
 * NOTE2: the returned string is allocated here!
 *
 * Note3: if more than one valid detector-string is found in the input, an error is returned
 *
 */
CHAR *
XLALGetChannelPrefix ( const CHAR *name )
{
  CHAR *channel = XLALCalloc( 3, sizeof(CHAR) );  /* 2 chars + \0 */

#define CHECK_UNIQUE do { if ( channel[0] != 0 ) XLAL_ERROR_NULL ( XLAL_EINVAL, "More than one matching detector name found in '%s'", name ); } while(0)

  if ( !channel ) {
    XLAL_ERROR_NULL ( XLAL_ENOMEM, "Failed to calloc(3)!\n" );
  }
  if ( !name ) {
    XLALFree ( channel );
    XLAL_ERROR_NULL ( XLAL_EINVAL, "Invalid NULL input 'name'" );
  }

  /* first handle (currently) unambiguous ones */
  if ( strstr( name, "ALLEGRO") || strstr ( name, "A1") ) {
    CHECK_UNIQUE;
    strcpy ( channel, "A1");
  }
  if ( strstr(name, "NIOBE") || strstr( name, "B1") ) {
    CHECK_UNIQUE;
    strcpy ( channel, "B1");
  }
  if ( strstr(name, "EXPLORER") || strstr( name, "E1") ) {
    CHECK_UNIQUE;
    strcpy ( channel, "E1");
  }
  if ( strstr(name, "GEO") || strstr(name, "G1") ) {
    CHECK_UNIQUE;
    strcpy ( channel, "G1" );
  }
  if ( strstr(name, "ACIGA") || strstr (name, "K1") ) {
    CHECK_UNIQUE;
    strcpy ( channel, "K1" );
  }
  if ( strstr(name, "LLO") || strstr(name, "Livingston") || strstr(name, "L1") ) {
    CHECK_UNIQUE;
    strcpy ( channel, "L1" );
  }
  if ( strstr(name, "Nautilus") || strstr(name, "N1") ) {
    CHECK_UNIQUE;
    strcpy ( channel, "N1" );
  }
  if ( strstr(name, "AURIGA") || strstr(name,"O1") ) {
    CHECK_UNIQUE;
    strcpy ( channel, "O1" );
  }
  if ( strstr(name, "CIT_40") || strstr(name, "Caltech-40") || strstr(name, "P1") ) {
    CHECK_UNIQUE;
    strcpy ( channel, "P1" );
  }
  if ( strstr(name, "TAMA") || strstr(name, "T1") ) {
    CHECK_UNIQUE;
    strcpy (channel, "T1" );
  }
  /* currently the only real ambiguity arises with H1 vs H2 */
  if ( strstr(name, "LHO") || strstr(name, "Hanford") || strstr(name, "H1") || strstr(name, "H2") ) {
    if ( strstr(name, "LHO_2k") ||  strstr(name, "H2") )
      {
        CHECK_UNIQUE;
        strcpy ( channel, "H2" );
      }
    if ( strstr(name, "LHO_4k") ||  strstr(name, "H1") )
      {
        CHECK_UNIQUE;
        strcpy ( channel, "H1" );
      }
    /* otherwise: guess */
    if ( channel[0] == 0 )
      {
        strcpy ( channel, "H1" );
        XLALPrintWarning("Detector-name '%s' ambiguous, guessing '%s'\n", name, channel );
      }
  } /* if LHO */
  /* LISA channel names are simply left unchanged */
  if ( strstr(name, "Z1") || strstr(name, "Z2") || strstr(name, "Z3")
       || strstr(name, "Z4") || strstr(name, "Z5") || strstr(name, "Z6")
       || strstr(name, "Z7") || strstr(name, "Z8") || strstr(name, "Z9") )
    {
      CHECK_UNIQUE;
      strncpy ( channel, name, 2);
      channel[2] = 0;
    }
  if ( strstr(name, "Virgo") || strstr(name, "VIRGO") || strstr(name, "V1") || strstr(name, "V2") )
    {
      if ( strstr(name, "Virgo_CITF") || strstr(name, "V2") )
        {
          CHECK_UNIQUE;
          strcpy ( channel, "V2" );
        }
      if ( strstr(name, "Virgo") || strstr(name, "VIRGO") || strstr(name, "V1") )
        {
          CHECK_UNIQUE;
          strcpy ( channel, "V1" );
        }
    } /* if Virgo */

  /* Did we fail to find any matches? */
  if ( channel[0] == 0 )
    XLAL_ERROR_NULL ( XLAL_EINVAL, "Unknown detector-name '%s'", name );
  else
    return channel;

} /* XLALGetChannelPrefix() */


/** Find the site geometry-information 'LALDetector' (mis-nomer!) given a detector-name.
 * The LALDetector struct is allocated here.
 */
LALDetector *
XLALGetSiteInfo ( const CHAR *name )
{
  CHAR *channel;
  LALDetector *site;

  /* first turn the free-form 'detector-name' into a well-defined channel-prefix */
  if ( ( channel = XLALGetChannelPrefix ( name ) ) == NULL ) {
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  if ( ( site = XLALCalloc ( 1, sizeof( *site) )) == NULL ) {
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  switch ( channel[0] )
    {
    case 'T':
      (*site) = lalCachedDetectors[LAL_TAMA_300_DETECTOR];
      break;
    case 'V':
      (*site) = lalCachedDetectors[LAL_VIRGO_DETECTOR];
      break;
    case 'G':
      (*site) = lalCachedDetectors[LAL_GEO_600_DETECTOR];
      break;
    case 'H':
      if ( channel[1] == '1' )
	(*site) = lalCachedDetectors[LAL_LHO_4K_DETECTOR];
      else
	(*site) = lalCachedDetectors[LAL_LHO_2K_DETECTOR];
      break;
    case 'L':
      (*site) = lalCachedDetectors[LAL_LLO_4K_DETECTOR];
      break;
    case 'P':
      (*site) = lalCachedDetectors[LAL_CIT_40_DETECTOR];
      break;
    case 'N':
      (*site) = lalCachedDetectors[LAL_NAUTILUS_DETECTOR];
      break;

    case 'Z':       /* create dummy-sites for LISA  */
      if ( XLALcreateLISA ( site, channel[1] ) != 0 )
	{
	  XLALPrintError("\nFailed to created LISA detector '%d'\n\n", channel[1]);
	  XLALFree ( site );
	  XLALFree ( channel );
	  XLAL_ERROR_NULL ( XLAL_EFUNC );
	}
      break;

    default:
      XLALPrintError ( "\nSorry, I don't have the site-info for '%c%c'\n\n", channel[0], channel[1]);
      XLALFree(site);
      XLALFree(channel);
      XLAL_ERROR_NULL ( XLAL_EINVAL );
      break;
    } /* switch channel[0] */

  XLALFree ( channel );

  return site;

} /* XLALGetSiteInfo() */

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


/** Destroy a MultiNoiseWeights object */
void
XLALDestroyMultiNoiseWeights ( MultiNoiseWeights *weights )
{
  if ( weights == NULL)
    return;

  for ( UINT4 k = 0; k < weights->length; k++ )
    XLALDestroyREAL8Vector ( weights->data[k] );

  XLALFree ( weights->data );
  XLALFree ( weights );

  return;

} /* XLALDestroyMultiNoiseWeights() */

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


#define OOTWOPI		(1.0 / LAL_TWOPI )
/** Interpolate frequency-series to newLen frequency-bins.
 *  This is using DFT-interpolation (derived from zero-padding).
 */
COMPLEX8Vector *
XLALrefineCOMPLEX8Vector (const COMPLEX8Vector *in,
			  UINT4 refineby,
			  UINT4 Dterms)
{
  UINT4 newLen, oldLen, l;
  COMPLEX8Vector *ret = NULL;

  if ( !in )
    return NULL;

  oldLen = in->length;
  newLen = oldLen * refineby;

  /* the following are used to speed things up in the innermost loop */
  if ( (ret = XLALCreateCOMPLEX8Vector ( newLen )) == NULL )
    return NULL;

  for (l=0; l < newLen; l++)
    {

      REAL8 kappa_l_k;
      REAL8 remain, kstarREAL;
      UINT4 kstar, kmin, kmax, k;
      REAL8 sink, coskm1;
      REAL8 Yk_re, Yk_im, Xd_re, Xd_im;

      kstarREAL = 1.0 * l  / refineby;
      kstar = (INT4)( kstarREAL + 0.5);	/* round to closest bin */
      kstar = MIN ( kstar, oldLen - 1 );	/* stay within the old SFT index-bounds */
      remain = kstarREAL - kstar;

      /* boundaries for innermost loop */
      kmin = MAX( 0, (INT4)kstar - (INT4)Dterms );
      kmax = MIN( oldLen, kstar + Dterms );

      Yk_re = Yk_im = 0;
      if ( fabs(remain) > 1e-5 )	/* denominater doens't vanish */
	{
	  /* Optimization: sin(2pi*kappa(l,k)) = sin(2pi*kappa(l,0) and idem for cos */
	  sink = sin ( LAL_TWOPI * remain );
	  coskm1 = cos ( LAL_TWOPI * remain ) - 1.0;

	  /* ---------- innermost loop: k over 2*Dterms around kstar ---------- */
	  for (k = kmin; k < kmax; k++)
	    {
	      REAL8 Plk_re, Plk_im;

	      Xd_re = in->data[k].re;
	      Xd_im = in->data[k].im;

	      kappa_l_k = kstarREAL - k;

	      Plk_re = sink / kappa_l_k;
	      Plk_im = coskm1 / kappa_l_k;

	      Yk_re += Plk_re * Xd_re - Plk_im * Xd_im;
	      Yk_im += Plk_re * Xd_im + Plk_im * Xd_re;

	    } /* hotloop over Dterms */
	}
      else	/* kappa -> 0: Plk = 2pi delta(k, l) */
	{
	  Yk_re = LAL_TWOPI * in->data[kstar].re;
	  Yk_im = LAL_TWOPI * in->data[kstar].im;
	}

      ret->data[l].re = OOTWOPI * Yk_re;
      ret->data[l].im = OOTWOPI * Yk_im;

    }  /* for l < newlen */

  return ret;

} /* XLALrefineCOMPLEX8Vector() */


/**
 * Function to read a segment list from given filename, returns a *sorted* SegmentList
 *
 * The segment-list format parse here is consistent with Xavie's segment lists used previously
 * and follows the format <repeated lines of form "startGPS endGPS duration[h] NumSFTs">,
 * allowed comment-characters are '%' and '#'
 *
 * \note we (ab)use the integer 'id' field in LALSeg to carry the total number of SFTs
 * contained in that segment. This can be used as a consistency check when loading SFTs for these segments.
 *
 */
LALSegList *
XLALReadSegmentsFromFile ( const char *fname	/**< name of file containing segment list */
                           )
{
  LALSegList *segList = NULL;

  /** check input consistency */
  if ( !fname ) {
    XLALPrintError ( "%s: NULL input 'fname'", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /* read and parse segment-list file contents*/
  LALParsedDataFile *flines = NULL;
  if ( XLALParseDataFile ( &flines, fname ) != XLAL_SUCCESS )
    XLAL_ERROR_NULL ( XLAL_EFUNC );

  UINT4 numSegments = flines->lines->nTokens;
  /* allocate and initialized segment list */
  if ( (segList = XLALCalloc ( 1, sizeof(*segList) )) == NULL )
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  if ( XLALSegListInit ( segList ) != XLAL_SUCCESS )
    XLAL_ERROR_NULL ( XLAL_EFUNC );


  UINT4 iSeg;
  for ( iSeg = 0; iSeg < numSegments; iSeg ++ )
    {
      REAL8 t0, t1, TspanHours;
      INT4 NSFT;
      LALSeg thisSeg;
      int ret;
      ret = sscanf ( flines->lines->tokens[iSeg], "%lf %lf %lf %d", &t0, &t1, &TspanHours, &NSFT );
      if ( ret != 4 ) {
        XLALPrintError ("%s: failed to parse data-line %d (%d) in segment-list %s: '%s'\n", __func__, iSeg, ret, fname, flines->lines->tokens[iSeg] );
        XLALSegListClear ( segList );
        XLALFree ( segList );
        XLALDestroyParsedDataFile ( flines );
        XLAL_ERROR_NULL ( XLAL_ESYS );
      }
      /* check internal consistency of these numbers */
      REAL8 hours = 3600.0;
      if ( fabs ( t1 - t0 - TspanHours * hours ) >= 1.0 ) {
        XLALPrintError ("%s: Inconsistent segment list, in line %d: t0 = %f, t1 = %f, Tspan = %f != t1 - t0 (to within 1s)\n", __func__, iSeg, t0, t1, TspanHours );
        XLAL_ERROR_NULL ( XLAL_EDOM );
      }

      LIGOTimeGPS start, end;
      XLALGPSSetREAL8( &start, t0 );
      XLALGPSSetREAL8( &end,   t1 );

      /* we set number of SFTs as 'id' field, as we have no other use for it */
      if ( XLALSegSet ( &thisSeg, &start, &end, NSFT ) != XLAL_SUCCESS )
        XLAL_ERROR_NULL ( XLAL_EFUNC );

      if ( XLALSegListAppend ( segList, &thisSeg ) != XLAL_SUCCESS )
        XLAL_ERROR_NULL ( XLAL_EFUNC );

    } /* for iSeg < numSegments */

  /* sort final segment list in increasing GPS start-times */
  if ( XLALSegListSort( segList ) != XLAL_SUCCESS )
    XLAL_ERROR_NULL ( XLAL_EFUNC );

  /* free parsed segment file contents */
  XLALDestroyParsedDataFile ( flines );

  return segList;

} /* XLALReadSegmentsFromFile() */

/* ============================================================
 * deprecated LAL interface API follow below
 * mostly these are now just LAL-wrappers to the corresponding
 * XLAL-inteface functions
 * ============================================================
 */

