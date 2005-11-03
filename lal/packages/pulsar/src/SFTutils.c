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
 * \author Reinhard Prix
 * \date 2005
 * \file 
 * \brief Utility functions for handling of SFTtype and SFTVector's.
 *
 * $Id$
 *
 */

/*---------- INCLUDES ----------*/
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>


#include "SFTutils.h"

NRCSID( SFTUTILSC, "$Id$" );

/*---------- DEFINES ----------*/
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))

#define TRUE (1==1)
#define FALSE (1==0)

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- empty initializers ---------- */

/*---------- Global variables ----------*/

/*---------- internal prototypes ----------*/

/*==================== FUNCTION DEFINITIONS ====================*/

/** Create one SFT-struct
 */
void
LALCreateSFTtype (LALStatus *status, 
		  SFTtype **output, 	/**< [out] allocated SFT-struct */
		  UINT4 SFTlen)		/**< number of frequency-bins */
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

  LALCCreateVector (status->statusPtr, &(sft->data), SFTlen);
  BEGINFAIL (status) { 
    LALFree (sft);
  } ENDFAIL (status);

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
		    UINT4 SFTlen)	/**< number of frequency-bins per SFT */
{
  UINT4 iSFT, j;
  SFTVector *vect;	/* vector to be returned */
  COMPLEX8Vector *data = NULL;

  INITSTATUS( status, "LALCreateSFTVector", SFTUTILSC);
  ATTATCHSTATUSPTR( status );

  ASSERT (output != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT (*output == NULL, status, SFTUTILS_ENONULL,  SFTUTILS_MSGENONULL);

  vect = LALCalloc (1, sizeof(*vect) );
  if (vect == NULL) {
    ABORT (status, SFTUTILS_EMEM, SFTUTILS_MSGEMEM);
  }

  vect->length = numSFTs;

  vect->data = LALCalloc (1, numSFTs * sizeof ( *vect->data ) );
  if (vect->data == NULL) {
    LALFree (vect);
    ABORT (status, SFTUTILS_EMEM, SFTUTILS_MSGEMEM);
  }

  for (iSFT=0; iSFT < numSFTs; iSFT ++)
    {
      LALCCreateVector (status->statusPtr, &data , SFTlen);
      BEGINFAIL (status) { /* crap, we have to de-allocate as far as we got so far... */
	for (j=0; j<iSFT; j++)
	  LALCDestroyVector (status->statusPtr, (COMPLEX8Vector**)&(vect->data[j].data) );
	LALFree (vect->data);
	LALFree (vect);
      } ENDFAIL (status);

      vect->data[iSFT].data = data;
      data = NULL;

    } /* for iSFT < numSFTs */

  *output = vect;

  DETATCHSTATUSPTR( status );
  RETURN (status);

} /* LALCreateSFTVector() */

/** Destroy an SFT-struct.
 */
void
LALDestroySFTtype (LALStatus *status, 
		   SFTtype **sft)	/**< SFT-struct to free */
{

  INITSTATUS( status, "LALDestroySFTtype", SFTUTILSC);
  ATTATCHSTATUSPTR (status);

  ASSERT (sft != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  /* be flexible: if points to null, nothing to do here.. (like 'free()') */
  if (*sft == NULL) {
    DETATCHSTATUSPTR( status );
    RETURN (status);
  }

  if ( (*sft)->data )
    {
      if ( (*sft)->data->data )
	LALFree ( (*sft)->data->data );
      LALFree ( (*sft)->data );
    }
  
  LALFree ( (*sft) );

  *sft = NULL;

  DETATCHSTATUSPTR( status );
  RETURN (status);

} /* LALDestroySFTtype() */


/** Destroy an SFT-vector
 */
void
LALDestroySFTVector (LALStatus *status, 
		     SFTVector **vect)	/**< the SFT-vector to free */
{
  UINT4 i;
  SFTtype *sft;

  INITSTATUS( status, "LALDestroySFTVector", SFTUTILSC);
  ATTATCHSTATUSPTR( status );

  ASSERT (vect != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT (*vect != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  
  for (i=0; i < (*vect)->length; i++)
    {
      sft = &( (*vect)->data[i] );
      if ( sft->data )
	{
	  if ( sft->data->data )
	    LALFree ( sft->data->data );
	  LALFree ( sft->data );
	}
    }

  LALFree ( (*vect)->data );
  LALFree ( *vect );

  *vect = NULL;

  DETATCHSTATUSPTR( status );
  RETURN (status);

} /* LALDestroySFTVector() */


/** Copy an entire SFT-type into another. 
 * We require the destination to have at least as many frequency-bins
 * as the source, but it can have less..
 */
void
LALCopySFT (LALStatus *status, 
	    SFTtype *dest, 	/**< [out] copied SFT (needs to be allocated already) */
	    const SFTtype *src)	/**< input-SFT to be copied */
{

  INITSTATUS( status, "LALCopySFT", SFTUTILSC);

  ASSERT (dest,  status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT (dest->data,  status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT (src, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT (src->data, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  /* some hard requirements */
  if ( dest->data->length < src->data->length ) {
    LALPrintError ("\nERROR: target-SFT has wrong size !\n\n");
    ABORT (status, SFTUTILS_EINPUT, SFTUTILS_MSGEINPUT);
  }
  
  /* copy head */
  strcpy (dest->name, src->name);
  dest->epoch = src->epoch;
  dest->f0 = src->f0;
  dest->deltaF = src->deltaF;
  dest->sampleUnits = src->sampleUnits;
  /* copy data */
  memcpy (dest->data->data, src->data->data, dest->data->length * sizeof (dest->data->data[0]));
  
  RETURN (status);

} /* LALCopySFT() */


/** Concatenate two SFT-vectors to a new one.
 */
void
LALConcatSFTVectors (LALStatus *status,
		     SFTVector **outVect,
		     const SFTVector *inVect1,
		     const SFTVector *inVect2 )
{
  UINT4 newlen;
  SFTVector *ret = NULL;

  INITSTATUS( status, "LALDestroySFTVector", SFTUTILSC);
  ATTATCHSTATUSPTR (status); 

  ASSERT (outVect,  status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT ( *outVect == NULL,  status, SFTUTILS_ENONULL,  SFTUTILS_MSGENONULL);
  ASSERT (inVect1,  status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT (inVect2,  status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  newlen = inVect1 -> length + inVect2 -> length;
    
  
  DETATCHSTATUSPTR (status); 
  RETURN (status);

} /* LALConcatSFTVectors() */






/** Allocate a LIGOTimeGPSVector
 */
void
LALCreateTimestampVector (LALStatus *status, 
			  LIGOTimeGPSVector **vect, 	/**< [out] allocated timestamp-vector  */
			  UINT4 length)			/**< number of elements */
{
  LIGOTimeGPSVector *out = NULL;

  INITSTATUS( status, "LALCreateTimestampVector", SFTUTILSC);

  ASSERT (vect != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT (*vect == NULL, status, SFTUTILS_ENONULL,  SFTUTILS_MSGENONULL);

  out = LALCalloc (1, sizeof(LIGOTimeGPSVector));
  if (out == NULL) {
    ABORT (status,  SFTUTILS_EMEM,  SFTUTILS_MSGEMEM);
  }
  out->length = length;
  out->data = LALCalloc (1, length * sizeof(LIGOTimeGPS));
  if (out->data == NULL) {
    LALFree (out);
    ABORT (status,  SFTUTILS_EMEM,  SFTUTILS_MSGEMEM);
  }

  *vect = out;

  RETURN (status);
  
} /* LALCreateTimestampVector() */


/** De-allocate a LIGOTimeGPSVector
 */
void
LALDestroyTimestampVector (LALStatus *status, 
			   LIGOTimeGPSVector **vect)	/**< timestamps-vector to be freed */
{
  INITSTATUS( status, "LALDestroyTimestampVector", SFTUTILSC);

  ASSERT (vect != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);
  ASSERT (*vect != NULL, status, SFTUTILS_ENULL,  SFTUTILS_MSGENULL);

  LALFree ( (*vect)->data);
  LALFree ( *vect );
  
  *vect = NULL;

  RETURN (status);
  
} /* LALDestroyTimestampVector() */



/** Given a start-time, duration and Tsft, returns a list of timestamps
 * covering this time-stretch.
 * returns NULL on out-of-memory
 */
void
LALMakeTimestamps(LALStatus *status,
		  LIGOTimeGPSVector **timestamps, 	/**< [out] timestamps-vector */
		  LIGOTimeGPS tStart, 			/**< GPS start-time */
		  REAL8 duration, 			/**< duration in seconds */
		  REAL8 Tsft)				/**< length of one (SFT) timestretch in seconds */
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

  ASSERT ( duration >= Tsft, status, SFTUTILS_EINPUT, SFTUTILS_MSGEINPUT);

  numSFTs = (UINT4)( duration / Tsft );			/* floor */
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
      /* NOTE: we add the interval Tsft successively instead of
       * via iSFT*Tsft, because this way we avoid possible ns-rounding problems
       * with a REAL8 interval, which becomes critial from about 100days on...
       */
      LALAddFloatToGPS( status->statusPtr, &tt, &tt, Tsft);
      BEGINFAIL( status ) {
	LALFree (ts->data);
	LALFree (ts);
      } ENDFAIL(status);

    } /* for i < numSFTs */

  *timestamps = ts;

  DETATCHSTATUSPTR( status );
  RETURN( status );
  
} /* LALMakeTimestamps() */
