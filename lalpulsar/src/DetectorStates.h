/*
 * Copyright (C) 2006 Reinhard Prix
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
 * \ingroup pulsarCommon
 * \date 2005
 * \file
 * \ingroup pulsar
 * \brief API for the DetectorStates.c functions.
 *
 * $Id$
 *
 */

#ifndef _DETECTORSTATES_H  /* Double-include protection. */
#define _DETECTORSTATES_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

#include <lal/LALRCSID.h>
NRCSID( DETECTORSTATESH, "$Id$" );

/*---------- exported INCLUDES ----------*/
#include <lal/LALComputeAM.h>
#include <lal/PulsarDataTypes.h>
#include <lal/SeqFactories.h>

/*---------- exported DEFINES ----------*/

/*----- Error-codes -----*/
#define DETECTORSTATES_ENULL 		1
#define DETECTORSTATES_ENONULL 		2
#define DETECTORSTATES_EINPUT  		3
#define DETECTORSTATES_EMEM   		4
#define DETECTORSTATES_EXLAL		5
#define DETECTORSTATES_EIEEE		6

#define DETECTORSTATES_MSGENULL		"Arguments contained an unexpected null pointer"
#define DETECTORSTATES_MSGENONULL 	"Output pointer is non-NULL"
#define DETECTORSTATES_MSGEINPUT   	"Invalid input"
#define DETECTORSTATES_MSGEMEM   	"Out of memory. Bad."
#define DETECTORSTATES_MSGEXLAL		"XLAL function call failed"
#define DETECTORSTATES_MSGEIEEE		"Floating point failure"

/*---------- exported types ----------*/

/** A symmetric 3x3 tensor (such as detector-tensors), storing only the upper triangle.
 */
typedef struct
{
  REAL4 d11;   REAL4 d12;   REAL4 d13;
               REAL4 d22;   REAL4 d23;
                            REAL4 d33;
} SymmTensor3;


/** A symmetric 3x3 tensor (such as detector-tensors), storing only the upper triangle, using REAL8 precision
 */
typedef struct
{
  REAL8 d11;   REAL8 d12;   REAL8 d13;
               REAL8 d22;   REAL8 d23;
                            REAL8 d33;
} SymmTensor3d;


/** Struct containing pre-computed quantites describing a
 * single detector arm: unit-vector along detector-arm, arm-length,
 * and arm "basis-tensor" n x n. This is used to speed up the
 * computation of LISA detector tensors in the rigid-adiabatic approximation.
 */
typedef struct
{
  REAL4 n[3];			/**< unit vector pointing along this arm */
  SymmTensor3 basisT;		/**< arm "basis-tensor" (n x n) */
  REAL4 armlength_c;		/**< armlengths in seconds L / c */
} DetectorArm;

typedef DetectorArm Detector3Arms[3];	/**< used to allow functions some type/size checking */

/* ----- Output types for LALGetDetectorStates() */
/** State-info about position, velocity and LMST of a detector together
 * with corresponding EarthState.
 */
typedef struct
{
  LIGOTimeGPS tGPS;		/**< GPS timestamps corresponding to this entry */
  REAL8 rDetector[3];		/**< Cartesian coords of detector position in ICRS J2000. Units=sec */
  REAL8 vDetector[3];		/**< Cart. coords. of detector velocity, in dimensionless units (v/c)*/
  REAL8 LMST;			/**< local mean sidereal time at the detector-location in radians */
  EarthState earthState;	/**< EarthState information */
  Detector3Arms detArms;	/**< include up to three arms to allow describing LISA */
  SymmTensor3 detT;		/**< Detector-tensor components in SSB-fixed, Cartesian coordinates */
} DetectorState;


/** Timeseries of DetectorState's, representing the detector-info at different timestamps.
 * In addition to the standard 'vector'-fields we also store the detector-info in here.
 */
typedef struct
{
  UINT4 length;			/**< total number of entries */
  DetectorState *data;		/**< array of DetectorState entries */
  LALDetector detector;		/**< detector-info corresponding to this timeseries */
  CoordinateSystem system; 	/**< The coordinate system used for detector's position/velocity and detector-tensor */
} DetectorStateSeries;

/** Multi-IFO time-series of DetectorStates */
typedef struct
{
  UINT4 length;			/**< number of detectors */
  DetectorStateSeries **data;	/**< vector of pointers to DetectorStateSeries */
  LIGOTimeGPS startTime;	/**< (earliest) startTime of the observation */
  REAL8 Tspan;			/**< total spanned duration of the observation */
} MultiDetectorStateSeries;

/*---------- exported Global variables ----------*/

/*---------- exported prototypes [API] ----------*/
void
LALGetDetectorStates (LALStatus *,
		      DetectorStateSeries **DetectorStates,
		      const LIGOTimeGPSVector *timestamps,
		      const LALDetector *detector,
		      const EphemerisData *edat,
		      REAL8 tOffset);

void
LALGetMultiDetectorStates( LALStatus *,
			   MultiDetectorStateSeries **mdetStates,
			   const MultiSFTVector *multiSFTs,
			   const EphemerisData *edat );

void LALCreateDetectorStateSeries (LALStatus *, DetectorStateSeries **vect, UINT4 length );

int XLALAddSymmTensor3s ( SymmTensor3 *sum, const SymmTensor3 *aT, const SymmTensor3 *bT );
int XLALSubtractSymmTensor3s ( SymmTensor3 *diff, const SymmTensor3 *aT, const SymmTensor3 *bT );
int XLALScaleSymmTensor3   ( SymmTensor3 *mult, const SymmTensor3 *aT, REAL4 factor );
int XLALTensorSquareVector3 ( SymmTensor3 *vxv, REAL4 v[3] );
int XLALSymmetricTensorProduct3 ( SymmTensor3 *vxw, REAL4 v[3], REAL4 w[3] );
REAL4 XLALContractSymmTensor3s ( const SymmTensor3 *T1, const SymmTensor3 *T2 );


/* destructors */
void XLALDestroyDetectorStateSeries ( DetectorStateSeries *detStates );
void LALDestroyDetectorStateSeries(LALStatus *, DetectorStateSeries **vect );
void XLALDestroyMultiDetectorStateSeries ( MultiDetectorStateSeries *mdetStates );


/* helpers */

void LALGetMultiDetectorVelTimePos(LALStatus                *status,
				   REAL8VectorSequence      **outVel,
				   REAL8VectorSequence      **outPos,
				   LIGOTimeGPSVector        **outTime,
				   MultiDetectorStateSeries *in);

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
