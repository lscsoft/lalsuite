/*
 * Copyright (C) 2006, 2008 John T. Whelan
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

/*---------- INCLUDES ----------*/
#define __USE_ISOC99 1
#include <math.h>

#include <lal/CWFastMath.h>
#include <lal/DetectorStates.h>
#include <lal/LISAspecifics.h>

/*---------- local DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)

/*----- Macros ----- */
#define SQUARE(x) ((x) * (x))
#define INIT_MEM(x) memset(&(x), 0, sizeof((x)))

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- empty initializers ---------- */

/*---------- Global variables ----------*/

/*---------- internal prototypes ----------*/
int XLALFillDetectorTensor (DetectorState *detState, const LALDetector *detector );	/* no need to export this ... */

/*==================== FUNCTION DEFINITIONS ====================*/

/**
 * \deprecated Use XLALGetDetectorStates() instead
 */
void
LALGetDetectorStates (LALStatus *status,			/**< pointer to LALStatus structure */
		      DetectorStateSeries **DetectorStates,	/**< [out] series of DetectorStates */
		      const LIGOTimeGPSVector *timestamps,	/**< array of GPS timestamps t_i */
		      const LALDetector *detector,		/**< detector info */
		      const EphemerisData *edat,		/**< ephemeris file data */
		      REAL8 tOffset				/**< compute detector states at timestamps SHIFTED by tOffset */
		      )
{
  INITSTATUS(status);

  ASSERT ( DetectorStates != NULL, status,  DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);
  ASSERT ( *DetectorStates == NULL, status,  DETECTORSTATES_ENONULL, DETECTORSTATES_MSGENONULL);

  ASSERT ( timestamps, status, DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);
  ASSERT ( detector, status, DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);
  ASSERT ( edat, status, DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);

  /* call XLAL function */
  DetectorStateSeries *ret = NULL;
  if ( ( ret = XLALGetDetectorStates ( timestamps, detector, edat, tOffset )) == NULL ) {
    XLALPrintError ("%s: XLALGetDetectorStates() failed with code=%d\n", __func__, xlalErrno );
    ABORT ( status, DETECTORSTATES_EXLAL, DETECTORSTATES_MSGEXLAL );
  }

  /* return result */
  (*DetectorStates) = ret;

  RETURN (status);

} /* LALGetDetectorStates() */


/**
 * Function to compute the LWL detector-tensor for the given \a detector in
 * SSB-fixed cartesian coordinates at time tgps.
 * The coordinates used are: EQUATORIAL for Earth-based detectors, but ECLIPTIC for LISA.
 * RETURN: 0 = OK, -1 = ERROR
 */
int
XLALFillDetectorTensor (DetectorState *detState,	/**< [out,in]: detector state: fill in detector-tensor */
			const LALDetector *detector	/**< [in]: which detector */
			)
{
  const CHAR *prefix;

  if ( !detState || !detector ) {
    xlalErrno = XLAL_EINVAL;
    return -1;
  }

  prefix = detector->frDetector.prefix;

  /* we need to distinguish two cases: space-borne (i.e. LISA) and Earth-based detectors */
  if ( prefix[0] == 'Z' )	/* LISA */
    {
      if ( XLALprecomputeLISAarms ( detState ) != 0 ) {
	XLALPrintError ("\nXLALprecomputeLISAarms() failed !\n\n");
	xlalErrno = XLAL_EINVAL;
	return -1;
      }

      if ( XLALgetLISADetectorTensorLWL ( &(detState->detT), detState->detArms, prefix[1] ) != 0 ) {
	XLALPrintError ("\nXLALgetLISADetectorTensorLWL() failed !\n\n");
	xlalErrno = XLAL_EINVAL;
	return -1;
      }

    } /* if LISA */
  else
    {
      REAL4 sinG, cosG, sinGcosG, sinGsinG, cosGcosG;
      SymmTensor3 *detT = &(detState->detT);

      XLAL_CHECK( XLALSinCosLUT ( &sinG, &cosG, detState->earthState.gmstRad ) == XLAL_SUCCESS, XLAL_EFUNC );
      sinGsinG = sinG * sinG;
      sinGcosG = sinG * cosG;
      cosGcosG = cosG * cosG;

      /*
      printf("GMST = %fdeg; cosG = %f, sinG= %f\n",
	     LAL_180_PI * atan2(sinG,cosG), cosG, sinG);
      */

      detT->d11 = detector->response[0][0] * cosGcosG
            - 2 * detector->response[0][1] * sinGcosG
                + detector->response[1][1] * sinGsinG;
      detT->d22 = detector->response[0][0] * sinGsinG
            + 2 * detector->response[0][1] * sinGcosG
                + detector->response[1][1] * cosGcosG;
      detT->d12 = (detector->response[0][0] - detector->response[1][1])
                                           * sinGcosG
	        + detector->response[0][1] * (cosGcosG - sinGsinG);
      detT->d13 = detector->response[0][2] * cosG
                - detector->response[1][2] * sinG;
      detT->d23 = detector->response[0][2] * sinG
                + detector->response[1][2] * cosG;
      detT->d33 = detector->response[2][2];

      /*
      printf("d = (%f %f %f\n",detT->d11,detT->d12,detT->d13);
      printf("     %f %f %f\n",detT->d12,detT->d22,detT->d23);
      printf("     %f %f %f)\n",detT->d13,detT->d23,detT->d33);

      printf("d*= (%f %f %f\n",detector->response[0][0],
	     detector->response[0][1],detector->response[0][2]);
      printf("     %f %f %f\n",detector->response[1][0],
	     detector->response[1][1],detector->response[1][2]);
      printf("     %f %f %f)\n",detector->response[2][0],
	     detector->response[2][1],detector->response[2][2]);
      */

    } /* if Earth-based */

  return 0;

} /* XLALFillDetectorTensor() */

/**
 * Compute the "squared-tensor" v x v for given vector v,
 * the result is returned in a "detectorTensor" struct
 */
int
XLALTensorSquareVector3 ( SymmTensor3 *vxv, REAL4 v[3] )
{
  if ( !vxv )
    return -1;

  vxv->d11 = v[0] * v[0];
  vxv->d12 = v[0] * v[1];
  vxv->d13 = v[0] * v[2];

  vxv->d22 = v[1] * v[1];
  vxv->d23 = v[1] * v[2];

  vxv->d33 = v[2] * v[2];

  return 0;

} /* XLALTensorSquareVector3() */

/**
 * Compute the symmetrized tensor product T = v x w + w x v
 */
int
XLALSymmetricTensorProduct3 ( SymmTensor3 *vxw, REAL4 v[3], REAL4 w[3] )
{
  if ( !vxw )
    return -1;

  vxw->d11 = 2.0f * v[0] * w[0];
  vxw->d12 = v[0] * w[1] + w[0] * v[1];
  vxw->d13 = v[0] * w[2] + w[0] * v[2];

  vxw->d22 = 2.0f * v[1] * w[1];
  vxw->d23 = v[1] * w[2] + w[1] * v[2];

  vxw->d33 = 2.0f * v[2] * w[2];

  return 0;

} /* XLALSymmTensorProduct() */

/**
 * Convenience function for adding two SymmTensor3s: aT + bT
 * NOTE: it *is* safe to have sum point to the same tensor-struct as either aT or bT.
 */
int
XLALAddSymmTensor3s ( SymmTensor3 *sum, const SymmTensor3 *aT, const SymmTensor3 *bT )
{
  if ( !sum || !aT || !bT )
    return -1;

  sum->d11 = aT->d11  + bT->d11;
  sum->d12 = aT->d12  + bT->d12;
  sum->d13 = aT->d13  + bT->d13;

  sum->d22 = aT->d22  + bT->d22;
  sum->d23 = aT->d23  + bT->d23;

  sum->d33 = aT->d33  + bT->d33;

  return 0;

} /* XLALAddSymmTensor3s() */

/**
 * Convenience function for subtracting two SymmTensor3s: aT - bT
 * NOTE: it *is* safe to have diff point to the same tensor-struct as either aT or bT.
 */
int
XLALSubtractSymmTensor3s ( SymmTensor3 *diff, const SymmTensor3 *aT, const SymmTensor3 *bT )
{
  if ( !diff || !aT || !bT )
    return -1;

  diff->d11 = aT->d11  - bT->d11;
  diff->d12 = aT->d12  - bT->d12;
  diff->d13 = aT->d13  - bT->d13;

  diff->d22 = aT->d22  - bT->d22;
  diff->d23 = aT->d23  - bT->d23;

  diff->d33 = aT->d33  - bT->d33;

  return 0;

} /* XLALSubtractSymmTensor3s() */

/**
 * Convenience function for multiplying a SymmTensor3 by a scalar factor.
 * NOTE: it *is* safe to have aT and mult point to the same tensor-struct
 */
int
XLALScaleSymmTensor3 ( SymmTensor3 *mult, const SymmTensor3 *aT, REAL4 factor )
{
  if ( !mult || !aT )
    return -1;

  mult->d11 = factor * aT->d11;
  mult->d12 = factor * aT->d12;
  mult->d13 = factor * aT->d13;

  mult->d22 = factor * aT->d22;
  mult->d23 = factor * aT->d23;

  mult->d33 = factor * aT->d33;

  return 0;

} /* XLALScaleSymmTensor3() */

/**
 * Contract two symmetric tensors over both indices T1 : T2
 */
REAL4
XLALContractSymmTensor3s ( const SymmTensor3 *T1, const SymmTensor3 *T2 )
{
  REAL4 ret;

  if ( !T1 || !T2 )
    XLALREAL4FailNaN();

  ret = T1->d11 * T2->d11
    + T1->d22 * T2->d22
    + T1->d33 * T2->d33
    + 2.0f * ( T1->d12 * T2->d12
	       + T1->d13 * T2->d13
	       + T1->d23 * T2->d23 );

  return ret;

} /* XLALContractSymmTensor3s() */


/* ===== Multi-IFO versions of some of the above functions ===== */

/**
 * Deprecated LAL wrapper to XLALGetMultiDetectorStates().
 * Get the detector-time series for the given MultiSFTVector.
 * (see LALGetDetectorStates() for more comments).
 *
 * \note The time-series is based on the <em>midpoints</em> of the SFT-timestamps.
 */
void
LALGetMultiDetectorStates( LALStatus *status,				/**< pointer to LALStatus structure */
			   MultiDetectorStateSeries **mdetStates, 	/**< [out] multi-IFO detector-states */
			   const MultiSFTVector *multiSFTs, 		/**< [in] multi-IFO SFTs */
			   const EphemerisData *edat )			/**< ephemeris data */
{
  INITSTATUS(status);

  ASSERT ( mdetStates, status, DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);
  ASSERT ( multiSFTs, status, DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);
  ASSERT ( *mdetStates == NULL, status, DETECTORSTATES_ENONULL, DETECTORSTATES_MSGENONULL);

  /* NOTE API change: XLAL-function wants detector-array and timestamps-vectors directly,
   * instead of a multi-SFT vector. We therefore need to extract this info from the
   * multi-SFT vector first
   */
  MultiLALDetector multiIFO;
  if ( XLALMultiLALDetectorFromMultiSFTs ( &multiIFO, multiSFTs ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALMultiLALDetectorFromMultiSFTs() failed with code %d\n", __func__, xlalErrno );
    ABORT ( status, DETECTORSTATES_EXLAL, DETECTORSTATES_MSGEXLAL );
  }

  MultiLIGOTimeGPSVector *multiTS;
  if ( ( multiTS = XLALExtractMultiTimestampsFromSFTs ( multiSFTs )) == NULL ) {
    XLALPrintError ("%s: XLALExtractMultiTimestampsFromSFTs() failed with code %d\n", __func__, xlalErrno );
    ABORT ( status, DETECTORSTATES_EXLAL, DETECTORSTATES_MSGEXLAL );
  }

  /* ---------- central wrapper: call XLAL function XLALGetMultiDetectorStates() ---------- */
  MultiDetectorStateSeries *ret = NULL;
  /* the API of this LAL interface specifies a hardcoded shift by Tsft/2 of all timestamps */
  REAL8 Tsft = 1.0 / multiSFTs->data[0]->data[0].deltaF;
  REAL8 tOffset = 0.5 * Tsft;
  if ( ( ret = XLALGetMultiDetectorStates( multiTS, &multiIFO, edat, tOffset )) == NULL ) {
    XLALPrintError ("%s: XLALGetMultiDetectorStates() failed with code %d\n", __func__, xlalErrno );
    XLALDestroyMultiTimestamps ( multiTS );
    ABORT ( status, DETECTORSTATES_EXLAL, DETECTORSTATES_MSGEXLAL );
  }

  /* free temporary mem */
  XLALDestroyMultiTimestamps ( multiTS );

  (*mdetStates) = ret;

  RETURN ( status );

} /* LALGetMultiDetectorStates() */



/* ===== Object creation/destruction functions ===== */

/** Create a DetectorStateSeries */
void
LALCreateDetectorStateSeries (LALStatus *status,		/**< pointer to LALStatus structure */
			      DetectorStateSeries **vect,	/**< output vector */
			      UINT4 length )			/**< number of entries */
{
  DetectorStateSeries *ret = NULL;

  INITSTATUS(status);

  ASSERT ( vect, status, DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);
  ASSERT ( *vect == NULL, status, DETECTORSTATES_ENONULL, DETECTORSTATES_MSGENONULL);

  if ( (ret = XLALCreateDetectorStateSeries ( length )) == NULL ) {
    ABORT ( status, DETECTORSTATES_EXLAL, DETECTORSTATES_MSGEXLAL );
  }

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
LALDestroyDetectorStateSeries (LALStatus *status,		/**< pointer to LALStatus structure */
			       DetectorStateSeries **detStates ) /**< pointer to vector to be destroyed */
{
  INITSTATUS(status);

  ASSERT ( detStates, status, DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);

  XLALDestroyDetectorStateSeries ( (*detStates) );

  (*detStates) = NULL;

  RETURN (status);
} /* LALDestroyDetectorStateSeries() */

/**
 * Helper function to get rid of a multi-IFO DetectorStateSeries
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


/**
 * Helper funxtion to copy velocity, time and position vectors out of the
 * multi-detector state series
 */
void LALGetMultiDetectorVelTimePos(LALStatus                *status,
				   REAL8VectorSequence      **outVel,
				   REAL8VectorSequence      **outPos,
				   LIGOTimeGPSVector        **outTime,
				   MultiDetectorStateSeries *in)
{

  UINT4 numifo, len, numsft, iIFO, iSFT, j;
  REAL8VectorSequence *velV = NULL;
  REAL8VectorSequence *posV = NULL;
  LIGOTimeGPSVector *timeV = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (in, status, DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);
  ASSERT (in->length > 0, status, DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);

  ASSERT (*outVel == NULL, status, DETECTORSTATES_ENONULL, DETECTORSTATES_MSGENONULL);
  ASSERT (*outPos == NULL, status, DETECTORSTATES_ENONULL, DETECTORSTATES_MSGENONULL);
  ASSERT (*outTime == NULL, status, DETECTORSTATES_ENONULL, DETECTORSTATES_MSGENONULL);

  /* number of ifos */
  numifo = in->length;

  len = 0;
  /* calculate total number of data points */
  for (j = 0, iIFO = 0; iIFO < numifo; iIFO++ ) {
    ASSERT (in->data[iIFO], status, DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);
    len += in->data[iIFO]->length;
  }

  /* allocate memory for vectors */
  velV =  XLALCreateREAL8VectorSequence ( len, 3 );
  posV =  XLALCreateREAL8VectorSequence ( len, 3 );
  TRY (LALCreateTimestampVector ( status->statusPtr, &timeV,  len), status);

  /* copy the timestamps, weights, and velocity vector */
  for (j = 0, iIFO = 0; iIFO < numifo; iIFO++ ) {

    ASSERT (in->data[iIFO], status, DETECTORSTATES_ENULL, DETECTORSTATES_MSGENULL);

    numsft = in->data[iIFO]->length;
    for ( iSFT = 0; iSFT < numsft; iSFT++, j++) {

      velV->data[3*j] = in->data[iIFO]->data[iSFT].vDetector[0];
      velV->data[3*j+1] = in->data[iIFO]->data[iSFT].vDetector[1];
      velV->data[3*j+2] = in->data[iIFO]->data[iSFT].vDetector[2];

      posV->data[3*j] = in->data[iIFO]->data[iSFT].rDetector[0];
      posV->data[3*j+1] = in->data[iIFO]->data[iSFT].rDetector[1];
      posV->data[3*j+2] = in->data[iIFO]->data[iSFT].rDetector[2];

      /* mid time of sfts */
      timeV->data[j] = in->data[iIFO]->data[iSFT].tGPS;

    } /* loop over SFTs */

  } /* loop over IFOs */

  *outVel = velV;
  *outPos = posV;
  *outTime = timeV;

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}

/** Create a DetectorStateSeries with length entries */
DetectorStateSeries*
XLALCreateDetectorStateSeries ( UINT4 length )		/**< number of entries */
{
  DetectorStateSeries *ret = NULL;

  if ( (ret = LALCalloc(1, sizeof(DetectorStateSeries) )) == NULL ) {
    XLALPrintError ("%s: failed to LALCalloc(1, %d)\n", __func__, sizeof(DetectorStateSeries) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  if ( (ret->data = LALCalloc (length, sizeof(DetectorState) )) == NULL ) {
    XLALFree (ret);
    XLALPrintError ("%s: failed to LALCalloc(%d, %d)\n", __func__, length, sizeof(DetectorState) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  ret->length = length;

  return ret;

} /* XLALCreateDetectorStateSeries() */


/**
 * Get the 'detector state' (ie detector-tensor, position, velocity, etc) for the given
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
 */
DetectorStateSeries *
XLALGetDetectorStates ( const LIGOTimeGPSVector *timestamps,	/**< array of GPS timestamps t_i */
                        const LALDetector *detector,		/**< detector info */
                        const EphemerisData *edat,		/**< ephemeris file data */
                        REAL8 tOffset				/**< compute detector states at timestamps SHIFTED by tOffset */
                        )
{
  /* check input consistency */
  if ( !timestamps || !detector || !edat ) {
    XLALPrintError ("%s: invalid NULL input, timestamps=%p, detector=%p, edat=%p\n", __func__, timestamps, detector, edat );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /* prepare return vector */
  UINT4 numSteps = timestamps->length;
  DetectorStateSeries *ret = NULL;
  if ( ( ret = XLALCreateDetectorStateSeries ( numSteps )) == NULL ) {
    XLALPrintError ("%s: XLALCreateDetectorStateSeries(%d) failed.\n", __func__, numSteps );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  /* enter detector-info into the head of the state-vector */
  ret->detector = (*detector);

  /* set 'time-span' associated with each timestamp */
  ret->deltaT = timestamps->deltaT;

  /* set SSB coordinate system used: EQUATORIAL for Earth-based, ECLIPTIC for LISA */
  if ( detector->frDetector.prefix[0] == 'Z' )	/* LISA */
    ret->system = COORDINATESYSTEM_ECLIPTIC;
  else	/* Earth-based */
    ret->system = COORDINATESYSTEM_EQUATORIAL;

  /* now fill all the vector-entries corresponding to different timestamps */
  UINT4 i;
  for ( i=0; i < numSteps; i++ )
    {
      BarycenterInput baryinput;
      EmissionTime emit;
      DetectorState *state = &(ret->data[i]);
      EarthState *earth = &(state->earthState);
      LIGOTimeGPS tgps;

      /* shift timestamp by tOffset */
      tgps = timestamps->data[i];
      XLALGPSAdd(&tgps, tOffset);

      /*----- first get earth-state */
      if ( XLALBarycenterEarth ( earth, &tgps, edat ) != XLAL_SUCCESS ) {
        XLALDestroyDetectorStateSeries ( ret );
        XLALPrintError("%s: XLALBarycenterEarth() failed with xlalErrno=%d\n", __func__, xlalErrno );
        XLAL_ERROR_NULL ( XLAL_EFAILED );
      }

      /*----- then get detector-specific info */
      baryinput.tgps = tgps;
      baryinput.site = (*detector);
      baryinput.site.location[0] /= LAL_C_SI;
      baryinput.site.location[1] /= LAL_C_SI;
      baryinput.site.location[2] /= LAL_C_SI;
      baryinput.alpha = baryinput.delta = 0;	/* irrelevant */
      baryinput.dInv = 0;

      if ( XLALBarycenter ( &emit, &baryinput, earth) != XLAL_SUCCESS ) {
	XLALDestroyDetectorStateSeries( ret );
        XLALPrintError("%s: XLALBarycenterEarth() failed with xlalErrno=%d\n", __func__, xlalErrno );
        XLAL_ERROR_NULL ( XLAL_EFAILED );
      }

      /*----- extract the output-data from this */
      UINT4 j;
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
	XLALDestroyDetectorStateSeries(ret);
	XLALPrintError ( "%s: XLALFillDetectorTensor() failed ... errno = %d\n\n", __func__, xlalErrno );
	XLAL_ERROR_NULL ( XLAL_EFUNC );
      }

    } /* for i < numSteps */

  /* return result */
  return ret;

} /* XLALGetDetectorStates() */


/**
 * Get the detector-time series for the given MultiLIGOTimeGPSVector.
 * NOTE: contrary to the deprecated LALGetMultiDetectorStates() interface, this
 * function computes detector-states at the given timestamps shifted by tOffset
 *
 */
MultiDetectorStateSeries *
XLALGetMultiDetectorStates( const MultiLIGOTimeGPSVector *multiTS, /**< [in] multi-IFO timestamps */
                            const MultiLALDetector *multiIFO, 	   /**< [in] multi-IFO array holding detector info */
                            const EphemerisData *edat,		   /**< [in] ephemeris data */
                            REAL8 tOffset			   /**< [in] shift all timestamps by this amount */
                            )
{
  /* check input consistency */
  if ( !multiIFO || !multiTS || !edat ) {
    XLALPrintError ("%s: invalid NULL input (multiIFO=%p, multiTS=%p or edat=%p)\n", __func__, multiIFO, multiTS, edat );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  UINT4 numDetectors;
  numDetectors = multiIFO->length;
  if ( numDetectors != multiTS->length ) {
    XLALPrintError ("%s: inconsistent number of IFOs in 'multiIFO' (%d) and 'multiTS' (%d)\n", __func__, multiIFO->length, multiTS->length );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /* prepare return-structure */
  MultiDetectorStateSeries *ret = NULL;
  if ( ( ret = LALCalloc ( 1, sizeof( *ret ) )) == NULL ) {
    XLALPrintError ("%s: LALCalloc ( 1, sizeof(%d)) failed\n", __func__, sizeof(*ret) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  if ( ( ret->data = LALCalloc ( numDetectors, sizeof( *(ret->data) ) )) == NULL ) {
    XLALFree ( ret );
    XLALPrintError ("%s: LALCalloc ( %d, sizeof(%d)) failed\n", __func__, numDetectors, sizeof(*(ret->data)) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  ret->length = numDetectors;

  REAL8 t0=LAL_REAL4_MAX;
  REAL8 t1=0;
  REAL8 deltaT = multiTS->data[0]->deltaT;
  LIGOTimeGPS startTime = {0, 0};
  /* loop over detectors */
  UINT4 X;
  for ( X=0; X < numDetectors; X ++ )
    {
      LIGOTimeGPSVector *tsX = multiTS->data[X];
      const LALDetector *detX = &(multiIFO->sites[X]);

      if ( !tsX || !detX ) {
        XLALPrintError ("%s: invalid NULL data-vector tsX[%d] = %p, detX[%d] = %p\n", __func__, X, tsX, X, detX );
        XLAL_ERROR_NULL ( XLAL_EINVAL );
      }
      if ( multiTS->data[X]->deltaT != deltaT ) {
        XLALPrintError ("%s: inconsistent time-base multi-timeseries deltaT[%d]=%f  !=  deltaT[0] = %f\n", __func__, X, multiTS->data[X]->deltaT, deltaT );
        XLAL_ERROR_NULL ( XLAL_EINVAL );
      }

      /* fill in the detector-state series for this detector */
      if ( ( ret->data[X] = XLALGetDetectorStates ( tsX, detX, edat, tOffset )) == NULL ) {
        XLALPrintError ("%s: XLALGetDetectorStates() failed.\n", __func__ );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }

      /* keep track of earliest/latest timestamp in order to determine total Tspan */
      UINT4 numTS = tsX->length;
      REAL8 t0_X = XLALGPSGetREAL8( &tsX->data[0] );
      REAL8 t1_X = XLALGPSGetREAL8( &tsX->data[numTS-1] );

      if ( t0_X < t0 ) {
        t0 = t0_X;
        startTime = tsX->data[0];
      }
      if ( t1_X > t1 ) t1 = t1_X;

    } /* for X < numDetectors */

  ret->Tspan = t1 - t0 + deltaT;	/* total time spanned by all SFTs */
  ret->startTime = startTime;		/* earliest start-time of observation */

  return ret;

} /* XLALGetMultiDetectorStates() */

/**
 * Parse string-vectors (typically input by user) of N detector noise-floors \f$\sqrt{S_X}\f$
 * for detectors \f$X=1\ldots N\f$, where here we assume equal number of SFTs per detector
 * such that \f$S^{-1} = \frac{1}{N}\sum_{X=0}^{N-1} S_X^{-1}\f$.
 *
 * \note input of length(sqrtSX)=1 < numDetectors is valid: use that single number for all detectors,
 *  otherwise we enforce length(sqrtSX) == numDetectors.
 *
 * returns result in MultiNoiseFloor struct 'multiNoiseFloor'.
 */
int
XLALParseMultiNoiseFloor ( MultiNoiseFloor *multiNoiseFloor,	/**< [out] parsed multi-IFO noise floor info */
                           const LALStringVector *sqrtSX,	/**< [in] string-list of \f$\sqrt{S_X}\f$ for detectors \f$X\f$ */
                           UINT4 numDetectors			/**< [in] number of detectors. NOTE: length[sqrtSX] must be EITHER =numDetectors OR =1 */
                           )
{
  XLAL_CHECK ( multiNoiseFloor != NULL, XLAL_EINVAL );
  XLAL_CHECK ( sqrtSX != NULL, XLAL_EINVAL );
  XLAL_CHECK ( (numDetectors > 0) && (numDetectors <= PULSAR_MAX_DETECTORS), XLAL_EINVAL );
  UINT4 numSqrtSX = sqrtSX->length;
  XLAL_CHECK ( (numSqrtSX == numDetectors) || (numSqrtSX == 1), XLAL_EINVAL );

  /* initialize empty return struct */
  multiNoiseFloor->length = numDetectors;

  /* parse input strings and fill multiNoiseFloor */
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      UINT4 X0 = X % numSqrtSX;		// always = 0 if (numSqrtSX == 1), otherwise = X if (numSqrtSX==numDetectors)
      const char *sqrtSnStr = sqrtSX->data[X0];
      REAL8 sqrtSn;
      XLAL_CHECK ( sscanf ( sqrtSnStr , "%lf", &sqrtSn ) == 1, XLAL_EINVAL, "Failed to parse '%s' into REAL8\n", sqrtSnStr );
      XLAL_CHECK ( sqrtSn >= 0, XLAL_EDOM );
      multiNoiseFloor->sqrtSn[X] = sqrtSn;
    } /* for X < numDetectors */

  return XLAL_SUCCESS;

} /* XLALParseMultiNoiseFloor() */


/**
 * Parse string-vectors (typically input by user) of N detector-names
 * for detectors \f$X=1\ldots N\f$, returns a MultiLALDetector struct
 *
 */
int
XLALParseMultiLALDetector ( MultiLALDetector *detInfo,        /**< [out] parsed detector-info struct */
                            const LALStringVector *detNames   /**< [in] list of detector names */
                            )
{
  XLAL_CHECK ( detInfo != NULL, XLAL_EINVAL );
  XLAL_CHECK ( detNames != NULL, XLAL_EINVAL );
  UINT4 numDet = detNames->length;
  XLAL_CHECK ( (numDet > 0) && (numDet <= PULSAR_MAX_DETECTORS), XLAL_EINVAL );

  detInfo->length = numDet;

  /* parse input strings and fill detInfo */
  for ( UINT4 X = 0; X < numDet; X ++ )
    {
      LALDetector *ifo;
      /* first parse detector name */
      XLAL_CHECK ( (ifo = XLALGetSiteInfo ( detNames->data[X] ) ) != NULL, XLAL_EINVAL, "Failed to parse detector-name '%s'\n", detNames->data[X] );
      detInfo->sites[X] = (*ifo);	// struct copy
      XLALFree ( ifo );

    } /* for X < numDet */

  return XLAL_SUCCESS;

} /* XLALParseMultiLALDetector() */


/**
 * Extract multi detector-info from a given multi SFTCatalog view
*/
int
XLALMultiLALDetectorFromMultiSFTCatalogView ( MultiLALDetector *multiIFO,		//!< [out] list of detectors found in catalog
                                              const MultiSFTCatalogView *multiView	//!< [in] multi-IFO view of SFT catalog
                                              )
{
  // check input consistency
  XLAL_CHECK ( multiIFO != NULL, XLAL_EINVAL );
  XLAL_CHECK ( multiView != NULL, XLAL_EINVAL );
  UINT4 numIFOs = multiView->length;
  XLAL_CHECK ( (numIFOs > 0) && (numIFOs <= PULSAR_MAX_DETECTORS), XLAL_EINVAL );



  multiIFO->length = numIFOs;
  for ( UINT4 X=0; X < numIFOs; X ++ )
    {
      LALDetector *site;
      XLAL_CHECK ( (site = XLALGetSiteInfo ( multiView->data[X].data->header.name )) != NULL, XLAL_EFUNC );
      multiIFO->sites[X] = (*site);	 // struct-copy
      XLALFree ( site );
    } /* for X < numIFOs */

  return XLAL_SUCCESS;

} /* XLALMultiLALDetectorFromMultiSFTCatalogView() */


/**
 * Extract multi detector-info from a given multi SFTCatalog view
 */
int
XLALMultiLALDetectorFromMultiSFTs ( MultiLALDetector *multiIFO,	//!< [out] list of detectors found in catalog
                                    const MultiSFTVector *multiSFTs	//!< [in] multi-IFO SFT vector
                                    )
{
  // check input consistency
  XLAL_CHECK ( multiIFO != NULL, XLAL_EINVAL );
  XLAL_CHECK ( multiSFTs != NULL, XLAL_EINVAL );
  XLAL_CHECK ( multiSFTs->length > 0, XLAL_EINVAL );

  UINT4 numIFOs = multiSFTs->length;

  multiIFO->length = numIFOs;
  for ( UINT4 X=0; X < numIFOs; X ++ )
    {
      LALDetector *site;
      XLAL_CHECK ( (site = XLALGetSiteInfo ( multiSFTs->data[X]->data[0].name )) != NULL, XLAL_EFUNC );
      multiIFO->sites[X] = (*site);	 // struct-copy
      XLALFree ( site );
    } /* for X < numIFOs */

  return XLAL_SUCCESS;

} /* XLALMultiLALDetectorFromMultiSFTVector() */
