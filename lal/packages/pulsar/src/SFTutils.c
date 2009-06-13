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

/**
 * \author Reinhard Prix, Badri Krishnan, John T. Whelan
 * \date 2005, 2007
 * \file
 * \ingroup SFTfileIO
 * \brief Utility functions for handling of SFTtype and SFTVector's.
 *
 * $Id$
 *
 */

/*---------- INCLUDES ----------*/
#include <stdarg.h>

#include <gsl/gsl_sort_double.h>

#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/FrequencySeries.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/LISAspecifics.h>

#include "SFTutils.h"

NRCSID( SFTUTILSC, "$Id$" );

/*---------- DEFINES ----------*/
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))

#define TRUE (1==1)
#define FALSE (1==0)

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
  const CHAR *fn = "XLALExtractBandfromSFTs";
  REAL8 dFreq;
  UINT4 iMin, iMax, i0, numSFTs, numBinsIn, numBinsOut, iSFT;
  SFTVector *out;
  REAL8 f0Out;
  COMPLEX8Vector *sav;

  if ( !sfts || !sfts->data || (sfts->length==0) || !sfts->data[0].data ) {
    XLAL_ERROR_NULL( fn, XLAL_EINVAL );
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
    XLAL_ERROR_NULL( fn, XLAL_EINVAL );
  }

  numSFTs = sfts->length;
  f0Out = iMin * dFreq;
  numBinsOut = iMax - iMin + 1;

  if ( (out = XLALCreateSFTVector ( numSFTs, numBinsOut )) == NULL ) {
    XLAL_ERROR_NULL( fn, XLAL_ENOMEM );
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



/** Create one SFT-struct. Allows for numBins == 0.
 */
void
LALCreateSFTtype (LALStatus *status,
		  SFTtype **output, 	/**< [out] allocated SFT-struct */
		  UINT4 numBins)	/**< number of frequency-bins */
{
  SFTtype *sft = NULL;

  INITSTATUS( status, "LALCreateSFTtype", SFTUTILSC);
  ATTATCHSTATUSPTR( status );

  ASSERT (output != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT (*output == NULL, status, SFTUTILS_ENONULL,  SFTUTILS_MSGENONULL);

  sft = LALCalloc (1, sizeof(*sft) );
  if (sft == NULL) {
    ABORT (status, SFTUTILS_EMEM, SFTUTILS_MSGEMEM);
  }

  if ( numBins )
    {
      LALCCreateVector (status->statusPtr, &(sft->data), numBins);
      BEGINFAIL (status) {
	LALFree (sft);
      } ENDFAIL (status);
    }
  else
    sft->data = NULL;	/* no data, just header */

  *output = sft;

  DETATCHSTATUSPTR( status );
  RETURN (status);

} /* LALCreateSFTtype() */


/** Create a whole vector of \c numSFT SFTs with \c SFTlen frequency-bins
 */
void
LALCreateSFTVector (LALStatus *status,
		    SFTVector **output, /**< [out] allocated SFT-vector */
		    UINT4 numSFTs, 	/**< number of SFTs */
		    UINT4 numBins)	/**< number of frequency-bins per SFT */
{
  SFTVector *vect;

  INITSTATUS( status, "LALCreateSFTVector", SFTUTILSC);
  ATTATCHSTATUSPTR( status );

  ASSERT (output != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT (*output == NULL, status, SFTUTILS_ENONULL,  SFTUTILS_MSGENONULL);

  if ( (vect = XLALCreateSFTVector ( numSFTs, numBins )) == NULL ) {
    ABORT (status, SFTUTILS_EMEM, SFTUTILS_MSGEMEM);
  }

  *output = vect;

  DETATCHSTATUSPTR( status );
  RETURN (status);

} /* LALCreateSFTVector() */


/** XLAL function to create an SFTVector of \c numSFT SFTs with \c SFTlen frequency-bins
 */
SFTVector *
XLALCreateSFTVector (UINT4 numSFTs, 	/**< number of SFTs */
		     UINT4 numBins	/**< number of frequency-bins per SFT */
		     )
{
  const CHAR *fn = "XLALCreateSFTVector()";
  UINT4 iSFT;
  SFTVector *vect;

  if ( (vect = LALCalloc ( 1, sizeof(*vect) )) == NULL ) {
    XLAL_ERROR_NULL( fn, XLAL_ENOMEM );
  }

  vect->length = numSFTs;
  if ( (vect->data = LALCalloc (1, numSFTs * sizeof ( *vect->data ) )) == NULL ) {
    LALFree (vect);
    XLAL_ERROR_NULL( fn, XLAL_ENOMEM );
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
	      LALFree (vect->data);
	      LALFree (vect);
	      XLAL_ERROR_NULL( fn, XLAL_ENOMEM );
	    }
	}

      vect->data[iSFT].data = data;

    } /* for iSFT < numSFTs */

  return vect;

} /* XLALCreateSFTVector() */



void LALCreateMultiSFTVector ( LALStatus *status,
			       MultiSFTVector **out,  /**< [out] multi sft vector created */
			       UINT4 length,          /**< number of sft data points */
			       UINT4Vector *numsft    /**< number of sfts in each sftvect */
			       )
{

  UINT4 k, j, numifo;
  MultiSFTVector *multSFTVec=NULL;

  INITSTATUS (status, "LALCreateMultiSFTs", SFTUTILSC);
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

} /* LALLoadMultiSFTs() */


/** Destroy an SFT-struct.
 */
void
LALDestroySFTtype (LALStatus *status,
		   SFTtype **sft)	/**< SFT-struct to free */
{

  INITSTATUS( status, "LALDestroySFTtype", SFTUTILSC);
  ATTATCHSTATUSPTR (status);

  ASSERT (sft != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  if (*sft == NULL)
    goto finished;

  if ( (*sft)->data )
    {
      if ( (*sft)->data->data )
	LALFree ( (*sft)->data->data );
      LALFree ( (*sft)->data );
    }

  LALFree ( (*sft) );

  *sft = NULL;

 finished:
  DETATCHSTATUSPTR( status );
  RETURN (status);

} /* LALDestroySFTtype() */


/** Destroy an SFT-vector
 */
void
LALDestroySFTVector (LALStatus *status,
		     SFTVector **vect)	/**< the SFT-vector to free */
{
  INITSTATUS( status, "LALDestroySFTVector", SFTUTILSC);

  ASSERT (vect != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  XLALDestroySFTVector ( *vect );

  *vect = NULL;

  RETURN (status);

} /* LALDestroySFTVector() */


/** XLAL interface to destroy an SFTVector
 */
void
XLALDestroySFTVector (SFTVector *vect)
{
  UINT4 i;
  SFTtype *sft;

  if ( !vect )
    return;

  for (i=0; i < vect->length; i++)
    {
      sft = &( vect->data[i] );
      if ( sft->data )
	{
	  if ( sft->data->data )
	    LALFree ( sft->data->data );
	  LALFree ( sft->data );
	}
    }

  LALFree ( vect->data );
  LALFree ( vect );

  return;

} /* XLALDestroySFTVector() */



/** Destroy a PSD-vector
 */
void
LALDestroyPSDVector (LALStatus *status,
		     PSDVector **vect)	/**< the SFT-vector to free */
{
  UINT4 i;
  REAL8FrequencySeries *psd;

  INITSTATUS( status, "LALDestroyPSDVector", SFTUTILSC);
  ATTATCHSTATUSPTR( status );

  ASSERT (vect != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  if ( *vect == NULL )	/* nothing to be done */
    goto finished;

  for (i=0; i < (*vect)->length; i++)
    {
      psd = &( (*vect)->data[i] );
      if ( psd->data )
	{
	  if ( psd->data->data )
	    LALFree ( psd->data->data );
	  LALFree ( psd->data );
	}
    }

  LALFree ( (*vect)->data );
  LALFree ( *vect );

  *vect = NULL;

 finished:
  DETATCHSTATUSPTR( status );
  RETURN (status);

} /* LALDestroyPSDVector() */


/** Destroy a multi SFT-vector
 */
void
LALDestroyMultiSFTVector (LALStatus *status,
		          MultiSFTVector **multvect)	/**< the SFT-vector to free */
{
  UINT4 i;

  INITSTATUS( status, "LALDestroyMultiSFTVector", SFTUTILSC);
  ATTATCHSTATUSPTR( status );

  ASSERT (multvect != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  if ( *multvect == NULL )	/* nothing to be done */
    goto finished;

  for ( i = 0; i < (*multvect)->length; i++)
      LALDestroySFTVector( status->statusPtr, (*multvect)->data + i);

  LALFree( (*multvect)->data );
  LALFree( *multvect );

  *multvect = NULL;

 finished:
  DETATCHSTATUSPTR( status );
  RETURN (status);

} /* LALDestroyMultiSFTVector() */



/** Destroy a multi PSD-vector
 */
void
LALDestroyMultiPSDVector (LALStatus *status,
		          MultiPSDVector **multvect)	/**< the SFT-vector to free */
{
  UINT4 i;

  INITSTATUS( status, "LALDestroyMultiPSDVector", SFTUTILSC);
  ATTATCHSTATUSPTR( status );

  ASSERT (multvect != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  if ( *multvect == NULL )
    goto finished;

  for ( i = 0; i < (*multvect)->length; i++) {
    LALDestroyPSDVector( status->statusPtr, (*multvect)->data + i);
  }

  LALFree( (*multvect)->data );
  LALFree( *multvect );

  *multvect = NULL;

 finished:
  DETATCHSTATUSPTR( status );
  RETURN (status);

} /* LALDestroySFTVector() */



/** Copy an entire SFT-type into another.
 * We require the destination-SFT to have a NULL data-entry, as the
 * corresponding data-vector will be allocated here and copied into
 *
 * Note: the source-SFT is allowed to have a NULL data-entry,
 * in which case only the header is copied.
 */
void
LALCopySFT (LALStatus *status,
	    SFTtype *dest, 	/**< [out] copied SFT (needs to be allocated already) */
	    const SFTtype *src)	/**< input-SFT to be copied */
{


  INITSTATUS( status, "LALCopySFT", SFTUTILSC);
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
LALSubtractSFTVectors (LALStatus *status,
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

  INITSTATUS( status, "LALSubtractSFTVectors", SFTUTILSC);
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
      LALPrintError ("\nERROR: the SFT-vectors must have the same number of SFTs!\n\n");
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
	LALPrintError ("\nERROR: the SFTs must have the same number of frequency-bins!\n\n");
	goto failed;
      }
      if ( (epoch1.gpsSeconds != epoch2.gpsSeconds) || ( epoch1.gpsNanoSeconds != epoch2.gpsNanoSeconds ) ) {
	LALPrintError ("\nERROR: the SFTs must have the same epochs!\n\n");
	goto failed;
      }
      if ( Freq1 != Freq2 ) {
	LALPrintError ("\nERROR: the SFTs must have the same start frequency!\n\n");
	goto failed;
      }
      if ( deltaF1 != deltaF2 ) {
	LALPrintError ("\nERROR: the SFTs must have the same frequency-steps!\n\n");
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

      LALSnprintf ( name1Trunc, halfNameLength, "%s", inVect1->data[i].name );
      LALSnprintf ( name2Trunc, halfNameLength, "%s", inVect2->data[i].name );
      LALSnprintf ( prefix, (strlen("Xn:") + 1), "%s", inVect1->data[i].name );
      LALSnprintf ( ret->data[i].name, LALNameLength, "%s{%s}-{%s}", prefix, name1Trunc, name2Trunc );
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
(LALStatus *status,
 SFTVector **outVect,	          /**< [out] linear combo of SFT-vectors */
 SFTVector **inVects,	  /**< array of SFT-vectors */
 const COMPLEX16Vector *weights,  /**< vector of SFT-weights */
 const CHAR *outName)             /**< name for output vector */
{
  UINT4 numSFTs, numSFTVects;
  UINT4 i, j, k;
  SFTVector *ret = NULL;

  INITSTATUS( status, "LALLinearlyCombineSFTVectors", SFTUTILSC);
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
      LALPrintError ("\nERROR: must be combining at least one SFT Vector!\n\n");
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
	    LALPrintError ("\nERROR: the SFTs must have the same number of frequency-bins!\n\n");
	    goto failed;
	  }
	  if ( (epoch1.gpsSeconds != epoch2.gpsSeconds) || ( epoch1.gpsNanoSeconds != epoch2.gpsNanoSeconds ) ) {
	    LALPrintError ("\nERROR: the SFTs must have the same epochs!\n\n");
	    goto failed;
	  }
	  if ( Freq1 != Freq2 ) {
	    LALPrintError ("\nERROR: the SFTs must have the same start frequency!\n\n");
	    goto failed;
	  }
	  if ( deltaF1 != deltaF2 ) {
	    LALPrintError ("\nERROR: the SFTs must have the same frequency-steps!\n\n");
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
LALAppendSFT2Vector (LALStatus *status,
		     SFTVector *vect,		/**< destinatino SFTVector to append to */
		     const SFTtype *sft)	/**< the SFT to append */
{
  UINT4 oldlen;
  INITSTATUS( status, "LALAppendSFT2Vector", SFTUTILSC);
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
XLALCreateTimestampVector (UINT4 length)
{
  LIGOTimeGPSVector *out = NULL;

  out = LALCalloc (1, sizeof(LIGOTimeGPSVector));
  if (out == NULL)
    XLAL_ERROR_NULL ( "XLALCreateTimestampVector", XLAL_ENOMEM );

  out->length = length;
  out->data = LALCalloc (1, length * sizeof(LIGOTimeGPS));
  if (out->data == NULL) {
    LALFree (out);
    XLAL_ERROR_NULL ( "XLALCreateTimestampVector", XLAL_ENOMEM );
  }

  return out;

} /* XLALCreateTimestampVector() */



/** LAL-interface: Allocate a LIGOTimeGPSVector */
void
LALCreateTimestampVector (LALStatus *status,
			  LIGOTimeGPSVector **vect, 	/**< [out] allocated timestamp-vector  */
			  UINT4 length)			/**< number of elements */
{
  LIGOTimeGPSVector *out = NULL;

  INITSTATUS( status, "LALCreateTimestampVector", SFTUTILSC);

  ASSERT (vect != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT (*vect == NULL, status, SFTUTILS_ENONULL,  SFTUTILS_MSGENONULL);

  if ( (out = XLALCreateTimestampVector( length )) == NULL ) {
    XLALClearErrno();
    ABORT (status,  SFTUTILS_EMEM,  SFTUTILS_MSGEMEM);
  }

  *vect = out;

  RETURN (status);

} /* LALCreateTimestampVector() */


/** De-allocate a LIGOTimeGPSVector */
void
XLALDestroyTimestampVector ( LIGOTimeGPSVector *vect)
{
  if ( !vect )
    return;

  LALFree ( vect->data );
  LALFree ( vect );

  return;

} /* XLALDestroyTimestampVector() */


/** De-allocate a LIGOTimeGPSVector
 */
void
LALDestroyTimestampVector (LALStatus *status,
			   LIGOTimeGPSVector **vect)	/**< timestamps-vector to be freed */
{
  INITSTATUS( status, "LALDestroyTimestampVector", SFTUTILSC);

  ASSERT (vect != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  if ( *vect == NULL )
    goto finished;

  XLALDestroyTimestampVector ( (*vect) );

  (*vect) = NULL;

 finished:
  RETURN (status);

} /* LALDestroyTimestampVector() */



/** Given a start-time, duration and 'stepsize' tStep, returns a list of timestamps
 * covering this time-stretch.
 */
void
LALMakeTimestamps(LALStatus *status,
		  LIGOTimeGPSVector **timestamps, 	/**< [out] timestamps-vector */
		  LIGOTimeGPS tStart, 			/**< GPS start-time */
		  REAL8 duration, 			/**< duration in seconds */
		  REAL8 tStep)				/**< length of one (SFT) timestretch in seconds */
{
  UINT4 i;
  UINT4 numSFTs;
  LIGOTimeGPS tt;
  LIGOTimeGPSVector *ts = NULL;

  INITSTATUS( status, "LALMakeTimestamps", SFTUTILSC);
  ATTATCHSTATUSPTR (status);

  ASSERT (timestamps != NULL, status, SFTUTILS_ENULL,
	  SFTUTILS_MSGENULL);
  ASSERT (*timestamps == NULL,status, SFTUTILS_ENONULL,
	  SFTUTILS_MSGENONULL);

  numSFTs = ceil( duration / tStep );			/* >= 1 !*/
  if ( (ts = LALCalloc (1, sizeof( *ts )) ) == NULL ) {
    ABORT (status,  SFTUTILS_EMEM,  SFTUTILS_MSGEMEM);
  }

  ts->length = numSFTs;
  if ( (ts->data = LALCalloc (1, numSFTs * sizeof (*ts->data) )) == NULL) {
    ABORT (status,  SFTUTILS_EMEM,  SFTUTILS_MSGEMEM);
  }

  tt = tStart;	/* initialize to start-time */
  for (i = 0; i < numSFTs; i++)
    {
      ts->data[i] = tt;
      /* get next time-stamp */
      /* NOTE: we add the interval tStep successively (rounded correctly to ns each time!)
       * instead of using iSFT*Tsft, in order to avoid possible ns-rounding problems
       * with REAL8 intervals, which becomes critial from about 100days on...
       */
      XLALGPSAdd(&tt, tStep);	/* can't fail */

    } /* for i < numSFTs */

  *timestamps = ts;

  DETATCHSTATUSPTR( status );
  RETURN( status );

} /* LALMakeTimestamps() */


/** Extract timstamps-vector from the given SFTVector.
 */
void
LALGetSFTtimestamps (LALStatus *status,
		     LIGOTimeGPSVector **timestamps,	/**< [out] extracted timestamps */
		     const SFTVector *sfts )		/**< input SFT-vector  */
{
  UINT4 numSFTs;
  UINT4 i;
  LIGOTimeGPSVector *ret = NULL;

  INITSTATUS (status, "LALGetSFTtimestamps", SFTUTILSC );
  ATTATCHSTATUSPTR (status);

  ASSERT ( timestamps, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL );
  ASSERT ( sfts, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL );
  ASSERT ( sfts->length > 0, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL );
  ASSERT ( *timestamps == NULL, status, SFTUTILS_ENONULL, SFTUTILS_MSGENONULL );

  numSFTs = sfts->length;

  TRY ( LALCreateTimestampVector ( status->statusPtr, &ret, numSFTs ), status );

  for ( i=0; i < numSFTs; i ++ )
    ret->data[i] = sfts->data[i].epoch;

  /* done: return Ts-vector */
  (*timestamps) = ret;

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* LALGetSFTtimestamps() */




/** Extract/construct the unique 2-character "channel prefix" from the given
 * "detector-name", which unfortunately will not always follow any of the
 * official detector-naming conventions given in the Frames-Spec LIGO-T970130-F-E
 * This function therefore sometime has to do some creative guessing:
 *
 * NOTE: in case the channel-number can not be deduced from the name,
 * it is set to '1', and a warning will be printed if lalDebugLevel > 0.
 *
 * NOTE2: the returned string is allocated here!
 */
CHAR *
XLALGetChannelPrefix ( const CHAR *name )
{
  CHAR *channel = LALCalloc( 3, sizeof(CHAR) );  /* 2 chars + \0 */

  if ( !channel ) {
    XLAL_ERROR_NULL ( "XLALGetChannelPrefix", XLAL_ENOMEM );
  }
  if ( !name ) {
    LALFree ( channel );
    XLAL_ERROR_NULL ( "XLALGetChannelPrefix", XLAL_EINVAL );
  }

  /* first handle (currently) unambiguous ones */
  if ( strstr( name, "ALLEGRO") || strstr ( name, "A1") )
    strcpy ( channel, "A1");
  else if ( strstr(name, "NIOBE") || strstr( name, "B1") )
    strcpy ( channel, "B1");
  else if ( strstr(name, "EXPLORER") || strstr( name, "E1") )
    strcpy ( channel, "E1");
  else if ( strstr(name, "GEO") || strstr(name, "G1") )
    strcpy ( channel, "G1" );
  else if ( strstr(name, "ACIGA") || strstr (name, "K1") )
    strcpy ( channel, "K1" );
  else if ( strstr(name, "LLO") || strstr(name, "Livingston") || strstr(name, "L1") )
    strcpy ( channel, "L1" );
  else if ( strstr(name, "Nautilus") || strstr(name, "N1") )
    strcpy ( channel, "N1" );
  else if ( strstr(name, "AURIGA") || strstr(name,"O1") )
    strcpy ( channel, "O1" );
  else if ( strstr(name, "CIT_40") || strstr(name, "Caltech-40") || strstr(name, "P1") )
    strcpy ( channel, "P1" );
  else if ( strstr(name, "TAMA") || strstr(name, "T1") )
    strcpy (channel, "T1" );
  /* currently the only real ambiguity arises with H1 vs H2 */
  else if ( strstr(name, "LHO") || strstr(name, "Hanford") || strstr(name, "H1") || strstr(name, "H2") )
    {
      if ( strstr(name, "LHO_2k") ||  strstr(name, "H2") )
	strcpy ( channel, "H2" );
      else if ( strstr(name, "LHO_4k") ||  strstr(name, "H1") )
	strcpy ( channel, "H1" );
      else /* otherwise: guess */
	{
	  strcpy ( channel, "H1" );
	  if ( lalDebugLevel )
	    LALPrintError("WARNING: Detector-name '%s' not unique, guessing '%s'\n", name, channel );
	}
    } /* if LHO */
  /* LISA channel names are simply left unchanged */
  else if ( strstr(name, "Z1") || strstr(name, "Z2") || strstr(name, "Z3")
	    || strstr(name, "Z4") || strstr(name, "Z5") || strstr(name, "Z6")
	    || strstr(name, "Z7") || strstr(name, "Z8") || strstr(name, "Z9") )
    {
      strncpy ( channel, name, 2);
      channel[2] = 0;
    }
  /* try matching VIRGO last, because 'V1','V2' might be used as version-numbers
   * also in some input-strings */
  else if ( strstr(name, "Virgo") || strstr(name, "VIRGO") || strstr(name, "V1") || strstr(name, "V2") )
    {
      if ( strstr(name, "Virgo_CITF") || strstr(name, "V1") )
	strcpy ( channel, "V1" );
      else if ( strstr(name, "Virgo") || strstr(name, "VIRGO") || strstr(name, "V2") )
	strcpy ( channel, "V2" );
    } /* if Virgo */


  if ( channel[0] == 0 )
    {
      if ( lalDebugLevel ) LALPrintError ( "\nERROR: unknown detector-name '%s'\n\n", name );
      XLAL_ERROR_NULL ( "XLALGetChannelPrefix", XLAL_EINVAL );
    }
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
    XLAL_ERROR_NULL ( "XLALGetSiteInfo()", XLAL_EINVAL );
  }

  if ( ( site = LALCalloc ( 1, sizeof( *site) )) == NULL ) {
    XLAL_ERROR_NULL ( "XLALGetSiteInfo()", XLAL_ENOMEM );
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
	  LALPrintError("\nFailed to created LISA detector '%d'\n\n", channel[1]);
	  LALFree ( site );
	  LALFree ( channel );
	  XLAL_ERROR_NULL ( "XLALGetSiteInfo()", XLAL_EFUNC );
	}
      break;

    default:
      LALPrintError ( "\nSorry, I don't have the site-info for '%c%c'\n\n", channel[0], channel[1]);
      LALFree(site);
      LALFree(channel);
      XLAL_ERROR_NULL ( "XLALGetSiteInfo()", XLAL_EINVAL );
      break;
    } /* switch channel[0] */

  LALFree ( channel );

  return site;

} /* XLALGetSiteInfo() */


/** Computes weight factors arising from SFTs with different noise
    floors -- it multiplies an existing weight vector */
void LALComputeNoiseWeights  (LALStatus        *status,
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
  INITSTATUS (status, "LALComputeNoiseWeights", SFTUTILSC);
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
  ASSERT (weightV->length == sftVect->length, status,
	  SFTUTILS_EINPUT, SFTUTILS_MSGEINPUT);
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
 * floors -- it multiplies an existing weight vector
 */
void LALComputeMultiNoiseWeights  (LALStatus             *status,
				   MultiNoiseWeights     **out,
				   const MultiPSDVector  *rngmed,
				   UINT4                 blocksRngMed,
				   UINT4                 excludePercentile)
{
  REAL8 Tsft_Sn=0.0, Tsft_sumSn=0.0;
  UINT4 Y, X, alpha, k, numifos, numsfts, lengthsft, numsftsTot;
  MultiNoiseWeights *weights;
  REAL8 Tsft = 1.0 / rngmed->data[0]->data[0].deltaF;
  REAL8 Tsft_calS;	/* overall noise-normalization */

  INITSTATUS (status, "LALComputeMultiNoiseWeights", SFTUTILSC);
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

	  Tsft_Sn = 0.0;
	  for ( k = halfBlock + excludeIndex; k < lengthsft - halfBlock - excludeIndex; k++)
	    Tsft_Sn += thisrm->data->data[k];
	  Tsft_Sn /= lengthsft - 2*halfBlock - 2*excludeIndex;

	  Tsft_sumSn += Tsft_Sn; /* sumSn is just a normalization factor */

	  weights->data[X]->data[alpha] = 1.0/Tsft_Sn;
	} /* end loop over sfts for each ifo */

    } /* end loop over ifos */

  Tsft_calS = Tsft_sumSn/numsftsTot;	/* overall noise-normalization is the average Tsft * <S_{X,\alpha}> */

  /* make weights of order unity by myltiplying by sumSn/total number of sfts */
  for ( X = 0; X < numifos; X ++) {
    numsfts = weights->data[X]->length;
    for ( alpha = 0; alpha < numsfts; alpha ++)
      weights->data[X]->data[alpha] *= Tsft_calS;
  }

  weights->Sinv_Tsft = 0.5 * Tsft*Tsft / Tsft_calS;		/* 'Sinv * Tsft' normalization factor uses single-sided PSD!! */

  *out = weights;


  DETATCHSTATUSPTR (status);
   /* normal exit */
  RETURN (status);
}


void
LALDestroyMultiNoiseWeights  (LALStatus         *status,
			      MultiNoiseWeights **weights)
{
  UINT4 k;

  INITSTATUS (status, "LALDestroyMultiNoiseWeights", SFTUTILSC);
  ATTATCHSTATUSPTR (status);

  ASSERT ( weights != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  if (*weights == NULL)
    goto finished;

  for ( k = 0; k < (*weights)->length; k++)
    LALDDestroyVector (status->statusPtr, (*weights)->data + k);

  LALFree( (*weights)->data );
  LALFree(*weights);

  *weights = NULL;

 finished:
  DETATCHSTATUSPTR (status);
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
upsampleMultiSFTVector (LALStatus *status,
			  MultiSFTVector *inout,	/**< [in,out]: upsampled multi SFT-vector */
			  UINT4 upsample, 		/**< integer factor to upsample by */
			  UINT4 Dterms			/**< number of terms in Dirichlet kernel [on each side] */
			  )
{
  UINT4 X, numDet;

  INITSTATUS( status, "upsampleMultiSFTVector", SFTUTILSC );
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
upsampleSFTVector (LALStatus *status,
		     SFTVector *inout,		/**< [in,out]: upsampled SFT-vector */
		     UINT4 upsample, 		/**< integer factor to upsample by */
		     UINT4 Dterms		/**< number of terms in Dirichlet kernel [on each side] */
		     )
{
  UINT4 alpha, numSFTs;

  INITSTATUS( status, "upsampleSFTVector", SFTUTILSC );
  ATTATCHSTATUSPTR (status);

  ASSERT ( inout, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);
  ASSERT ( inout->length, status, SFTUTILS_ENULL, SFTUTILS_MSGENULL);

  numSFTs = inout->length;

  for ( alpha=0; alpha < numSFTs; alpha ++ )
    {
      COMPLEX8Vector *this_data = inout->data[alpha].data;
      COMPLEX8Vector *new_data;
      if ( (new_data = XLALrefineCOMPLEX8Vector ( this_data, upsample, Dterms )) == NULL ) {
	LALPrintError ("\nSFT oversampling failed ... \n\n");
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

#define LD_SMALL4       (1.0e-6)		/**< "small" number for REAL4*/
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
