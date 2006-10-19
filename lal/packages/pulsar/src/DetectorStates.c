/*
 * Copyright (C) 2005, 2006 Reinhard Prix
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

/** \author R. Prix
 * \file 
 * \brief
 * functions to handle DetectorStatesSeries: positions, velocities, detector-tensors
 * of detector as function of time.
 *
 */

/*---------- INCLUDES ----------*/
#define __USE_ISOC99 1
#include <math.h>

#include "DetectorStates.h"

NRCSID( DETECTORSTATESC, "$Id$");

/*---------- local DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)

/*----- Macros ----- */
/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- empty initializers ---------- */

/*---------- Global variables ----------*/

/*---------- internal prototypes ----------*/
int XLALFillDetectorTensor (DetectorState *detState, const LALDetector *detector );	/* no need to export this ... */

/*==================== FUNCTION DEFINITIONS ====================*/

/** Get the 'detector state' (ie detector-tensor, position, velocity, etc) for the given
 * vector of timestamps, shifted by a common time-shift \a tOffset.
 *
 * This function just calls LALBarycenterEarth() and LALBarycenter() for the
 * given vector of timestamps (shifted by tOffset) and returns the positions, 
 * velocities and LMSTs of the detector, stored in a DetectorStateSeries. 
 * There is also an entry containing the EarthState at each timestamp, which 
 * can be used as input for subsequent calls to LALBarycenter().
 *
 * \a tOffset allows one to easily use the midpoints of SFT-timestamps, for example.
 *
 * \note the DetectorStateSeries is allocated here and should be free'ed with
 * LALDestroyDetectorStateSeries().
 *
 */
void
LALGetDetectorStates (LALStatus *status,
		      DetectorStateSeries **DetectorStates,	/**< [out] series of DetectorStates */
		      const LIGOTimeGPSVector *timestamps,	/**< array of GPS timestamps t_i */
		      const LALDetector *detector,		/**< detector info */
		      const EphemerisData *edat,		/**< ephemeris file data */	
		      REAL8 tOffset
		      )
{
  UINT4 i, j, numSteps;
  DetectorStateSeries *ret = NULL;

  INITSTATUS( status, "LALGetDetectorStates", DETECTORSTATESC );
  ATTATCHSTATUSPTR (status);

  ASSERT ( DetectorStates != NULL, status,  DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);
  ASSERT ( *DetectorStates == NULL, status,  DETECTORSTATES_ENONULL, DETECTORSTATES_MSGENONULL);

  ASSERT ( timestamps, status, DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);
  ASSERT ( detector, status, DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);
  ASSERT ( edat, status, DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);

  numSteps = timestamps->length;

  TRY ( LALCreateDetectorStateSeries (status->statusPtr, &ret, numSteps), status);

  /* enter detector-info into the head of the state-vector */
  ret->detector = (*detector);


  /* now fill all the vector-entries corresponding to different timestamps */
  for ( i=0; i < numSteps; i++ )
    {
      BarycenterInput baryinput;
      EmissionTime emit;
      DetectorState *state = &(ret->data[i]);
      EarthState *earth = &(state->earthState);
      LIGOTimeGPS tgps;

      /* shift timestamp by tOffset */
      TRY ( LALAddFloatToGPS(status->statusPtr, &tgps, &timestamps->data[i], tOffset ), status);

      /*----- first get earth-state */
      LALBarycenterEarth (status->statusPtr, earth, &tgps, edat );
      BEGINFAIL(status){
	TRY ( LALDestroyDetectorStateSeries(status->statusPtr, &ret), status);
      }ENDFAIL(status);
      /*----- then get detector-specific info */
      baryinput.tgps = tgps;			
      baryinput.site = (*detector);
      baryinput.site.location[0] /= LAL_C_SI;
      baryinput.site.location[1] /= LAL_C_SI;
      baryinput.site.location[2] /= LAL_C_SI;
      baryinput.alpha = baryinput.delta = 0;	/* irrelevant */
      baryinput.dInv = 0;

      LALBarycenter (status->statusPtr, &emit, &baryinput, earth);
      BEGINFAIL(status) {
	TRY ( LALDestroyDetectorStateSeries(status->statusPtr, &ret), status);
      }ENDFAIL(status);

      /*----- extract the output-data from this */
      for (j=0; j < 3; j++)	/* copy detector's position and velocity */
	{
	  state->rDetector[j] = emit.rDetector[j];
	  state->vDetector[j] = emit.vDetector[j];
	} /* for j < 3 */

      /* local mean sidereal time = GMST + longitude */
      state->LMST = earth->gmstRad + detector->frDetector.vertexLongitudeRadians;
      state->LMST = fmod (state->LMST, LAL_TWOPI );	/* normalize */

      /* insert timestamp */
      state->tGPS = tgps;

      /* compute the detector-tensor at this time-stamp in SSB-fixed Cartesian coordinates
       * [EQUATORIAL for Earth-based, ECLIPTIC for LISA]  
       */
      if ( XLALFillDetectorTensor ( state, detector ) != 0 ) {
	LALPrintError ( "\nXLALFillDetectorTensor() failed ... errno = %d\n\n", xlalErrno );
	TRY ( LALDestroyDetectorStateSeries(status->statusPtr, &ret), status);
	ABORT ( status, DETECTORSTATES_EXLAL, DETECTORSTATES_MSGEXLAL );
      }

    } /* for i < numSteps */

  /* return result */
  (*DetectorStates) = ret;

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALGetDetectorStates() */





/** Function to compute the LWL detector-tensor for the given \a detector in 
 * SSB-fixed cartesian coordinates at time tgps.
 * The coordinates used are: EQUATORIAL for Earth-based detectors, but ECLIPTIC for LISA.
 * RETURN: 0 = OK, -1 = ERROR
 */
int
XLALFillDetectorTensor (DetectorState *detState,	/**< [out,in]: detector state: fill in detector-tensor */
			const LALDetector *detector	/**< [in]: which detector */
			)
{
  CHAR *channel;

  if ( !detState || !detector ) {
    xlalErrno = XLAL_EINVAL;
    return -1;
  }

  /* get a normalized descriptor for this detector from its name */
  if ( (channel = XLALGetChannelPrefix ( detector->frDetector.name )) == NULL ) {
    LALPrintError ("\nXLALGetChannelPrefix() failed for detector name '%s'\n\n", detector->frDetector.name);
    xlalErrno = XLAL_EINVAL;
    return -1;
  }

  /* we need to distinguish two cases: space-borne (i.e. LISA) and Earth-based detectors */
  if ( channel[0] == 'Z' )	/* LISA */
    { 
      LISAarmT armA, armB;
      /* we distinuish (currently) 3 different TDI observables: X (channel=1), Y (channel=2), Z (channel=3) */
      switch ( channel[1] )
	{
	case 1: 	/* TDI observable 'X' */
	  armA = LISA_ARM3; armB = LISA_ARM2;
	  break;
	case 2:		/* TDI observable 'Y' */
	  armA = LISA_ARM1; armB = LISA_ARM3;
	  break;
	case 3:		/* TDI observable 'Z' */
	  armA = LISA_ARM2; armB = LISA_ARM1;
	  break;
	default:	/* unknown */
	  xlalErrno = XLAL_EINVAL;
	  return -1;
	  break;
	} /* switch channel[1] */

      if ( XLALgetLISAtwoArmIFO ( &(detState->detT), detState->tGPS, armA, armB ) != 0 ) {
	LALPrintError ("\nXLALgetLISAtwoArmIFO() failed !\n\n");
	xlalErrno = XLAL_EINVAL;
	return -1;
      }


      
    } /* if LISA */
  else
    {
      /* ==> FIXME: EARTH-based case goes here <== */
      ;
    } /* if Earth-based */
    

  LALFree ( channel );

  return 0;

} /* XLALFillDetectorTensor() */



/* ===== Multi-IFO versions of some of the above functions ===== */

/** Get the detector-time series for the given MultiSFTVector. 
 * (see LALGetDetectorStates() for more comments).
 *
 * \note The time-series is based on the <em>midpoints</em> of the SFT-timestamps.
 */
void
LALGetMultiDetectorStates( LALStatus *status, 
			   MultiDetectorStateSeries **mdetStates, /**< [out] multi-IFO detector-states */
			   const MultiSFTVector *multiSFTs, 		/**< [in] multi-IFO SFTs */
			   const EphemerisData *edat )		/**< ephemeris files data nix nix*/
{
  UINT4 X, numDetectors;
  MultiDetectorStateSeries *ret = NULL;
  REAL8 t0=LAL_REAL4_MAX, t1=0;

  INITSTATUS (status, "LALGetMultiDetectorStates", DETECTORSTATESC );
  ATTATCHSTATUSPTR (status);

  ASSERT ( mdetStates, status, DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);
  ASSERT ( multiSFTs, status, DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);
  ASSERT ( *mdetStates == NULL, status, DETECTORSTATES_ENONULL, DETECTORSTATES_MSGENONULL);

  numDetectors = multiSFTs->length;

  /* prepare return-structure */
  if ( ( ret = LALCalloc ( 1, sizeof( *ret ) )) == NULL ) {
    ABORT (status, DETECTORSTATES_EMEM, DETECTORSTATES_MSGEMEM);
  }
  if ( ( ret->data = LALCalloc ( numDetectors, sizeof( *(ret->data) ) )) == NULL ) {
    LALFree ( ret );
    ABORT (status, DETECTORSTATES_EMEM, DETECTORSTATES_MSGEMEM);
  }
  ret->length = numDetectors;

  /* loop over detectors */
  for ( X=0; X < numDetectors; X ++ )
    {
      LIGOTimeGPSVector *ts = NULL;
      LALDetector *det = NULL;

      SFTVector *this_sftvect = multiSFTs->data[X];
      REAL8 Tsft = 1.0 / this_sftvect->data[0].deltaF;
      REAL8 tOffs = 0.5 * Tsft ;	/* shift all timestamps by Tsft / 2 */
      REAL8 ti;
      /* get earliest start-time, and latest end-time */
      ti = GPS2REAL8(this_sftvect->data[0].epoch);
      t0 = MYMIN(t0, ti );
      ti = GPS2REAL8(this_sftvect->data[this_sftvect->length-1].epoch) + Tsft;
      t1 = MYMAX(t1, ti );

      /* timestamps from SFTVector  of detector X */
      LALGetSFTtimestamps ( status->statusPtr, &ts, this_sftvect );
      if ( status->statusPtr->statusCode ) 
	{
	  LALPrintError ( "\nCall to LALGetSFTtimestamps() has failed ... \n\n");
	  goto failed;
	}

      /* LALDetector struct for this detector */
      if ( (det = XLALGetSiteInfo ( this_sftvect->data[0].name )) == NULL ) 
	{
	  LALPrintError ("\nCall to XLALGetSiteInfo() has failed ... \n\n");
	  XLALDestroyTimestampVector ( ts );
	  goto failed;
	}
      /* fill in the detector-state series for this detector */
      LALGetDetectorStates (status->statusPtr, &(ret->data[X]), ts, det, edat, tOffs );
      if ( status->statusPtr->statusCode ) 
	{
	  LALPrintError ( "\nCall to LALGetDetectorStates() has failed ... \n\n");
	  XLALDestroyTimestampVector ( ts );
	  LALFree ( det );
	  goto failed;
	}

      /* free temporary mem */
      XLALDestroyTimestampVector ( ts );
      ts = NULL;
      LALFree ( det );

    } /* for X < numDetectors */

  ret->Tspan = t1 - t0;		/* total time spanned by all SFTs */

  goto success;

 failed:
  /* free complete MultiDetectorStateSeries built up so far */
  XLALDestroyMultiDetectorStateSeries ( ret );	/* NOTE: this function is "NULL-robust" */
  ABORT ( status, -1, "LALGetMultiDetectorStates failed" );
  
 success:

  (*mdetStates) = ret;
  
  DETATCHSTATUSPTR (status);
  RETURN ( status );

} /* LALGetMultiDetectorStates() */



/* ===== Object creation/destruction functions ===== */

/** Create a DetectorStateSeries */
void
LALCreateDetectorStateSeries (LALStatus *status, 
			      DetectorStateSeries **vect,	/**< output vector */
			      UINT4 length )			/**< number of entries */
{
  DetectorStateSeries *ret = NULL;

  INITSTATUS (status, "LALCreateDetectorStateSeries", DETECTORSTATESC );

  ASSERT ( vect, status, DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);
  ASSERT ( *vect == NULL, status, DETECTORSTATES_ENONULL, DETECTORSTATES_MSGENONULL);

  if ( (ret = LALCalloc(1, sizeof(DetectorStateSeries) )) == NULL ) {
    ABORT (status, DETECTORSTATES_EMEM, DETECTORSTATES_MSGEMEM);
  }

  if ( (ret->data = LALCalloc (length, sizeof(DetectorState) )) == NULL ) {
    LALFree (ret);
    ABORT (status, DETECTORSTATES_EMEM, DETECTORSTATES_MSGEMEM);
  }

  ret->length = length;

  /* return result */
  (*vect) = ret;

  RETURN (status);

} /* LALCreateDetectorStateSeries() */

/* Get rid of a DetectorStateSeries */
void
XLALDestroyDetectorStateSeries ( DetectorStateSeries *detStates )
{
  if ( !detStates )
    return;

  if ( detStates->data ) LALFree ( detStates->data );
  LALFree ( detStates );

  return;

} /* XLALDestroyDetectorStateSeries() */

/** Destroy a DetectorStateSeries (and set it to NULL) */
void
LALDestroyDetectorStateSeries (LALStatus *status, 
			       DetectorStateSeries **detStates ) /**< pointer to vector to be destroyed */
{
  INITSTATUS (status, "LALDestroyDetectorStateSeries", DETECTORSTATESC );

  ASSERT ( detStates, status, DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);

  XLALDestroyDetectorStateSeries ( (*detStates) );

  (*detStates) = NULL;

  RETURN (status);
} /* LALDestroyDetectorStateSeries() */

/** Helper function to get rid of a multi-IFO DetectorStateSeries 
 * Note, this is "NULL-robust" in the sense that it will not crash 
 * on NULL-entries anywhere in this struct, so it can be used
 * for failure-cleanup even on incomplete structs.
 */
void
XLALDestroyMultiDetectorStateSeries ( MultiDetectorStateSeries *mdetStates )
{
  UINT4 X, numDet;

  if ( !mdetStates )
    return;

  numDet = mdetStates->length;
  if ( mdetStates->data )
    {
      for ( X=0; X < numDet ; X ++ )
	XLALDestroyDetectorStateSeries ( mdetStates->data[X] );

      LALFree ( mdetStates->data );
    }

  LALFree ( mdetStates );

  return;

} /* XLALDestroyMultiDetectorStateSeries() */

