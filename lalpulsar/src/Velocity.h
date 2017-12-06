/*
 *  Copyright (C) 2005 Badri Krishnan, Alicia Sintes
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

/* double inclusion protection */
#ifndef _VELOCITY_H
#define _VELOCITY_H

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \author Krishnan, B., Sintes, A.M.
 * \defgroup Velocity_h Header Velocity.h
 * \ingroup pkg_pulsarHough
 * \brief Computation of instant and averaged velocities for a given detector and the like.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/Velocity.h>
 * \endcode
 *
 * To find the velocity of a given detetector at a given time, or the averaged
 * velocity  of a detector in a certain time interval.
 *
 */
/*@{*/


/* *************
 *    Includes. This header may include others; if so, they go immediately
 *    after include-loop protection. Includes should appear in the following
 *    order:
 *    a. Standard library includes
 *    b. LDAS includes
 *    c. LAL includes
 */
#include<lal/Date.h>
#include<lal/LALDatatypes.h>
#include<lal/ComputeSky.h>
#include<lal/LALInitBarycenter.h>
#include<lal/LALBarycenter.h>
#include<lal/LALStdlib.h>
#include<lal/LALConstants.h>


/* ***************************************
 *   Error codes and messages. This must be auto-extracted for
 *    inclusion in the documentation.
 */
/**\name Error Codes */
/*@{*/
#define VELOCITYH_ENULL 1
#define VELOCITYH_EVAL 2
#define VELOCITYH_MSGENULL "Null Pointer"
#define VELOCITYH_MSGEVAL "Invalid Value"
/*@}*/

/* *****************************************************
 *   Structure, enum, union, etc., typdefs.
 */

/**
 * This structure stores the parameters required by LALBarycenter() to calculate
 * Earth velocity at a given detector location.
 */
typedef struct tagVelocityPar {
  LALDetector    detector; 	/**< the detector */
  EphemerisData  *edat;  	/**< ephemeris data pointer from LALInitBarycenter() */
  LIGOTimeGPS    startTime; 	/**< start of time interval */
  REAL8          tBase; 	/**< duration of interval */
  REAL8          vTol;  	/**< fractional accuracy required for velocity (redundant for average velocity calculation) */
} VelocityPar;

/* ***************************************************
 *  Functions Declarations (i.e., prototypes).
 */
void LALAvgDetectorVel(LALStatus    *status,
		    REAL8        v[3], /* output vector representing average velocity */
		    VelocityPar  *in); /* parameters required to calculate V */

void LALAvgDetectorPos(LALStatus    *status,
		    REAL8        x[3], /* output vector representing average position */
		    VelocityPar  *in); /* parameters required to calculate position */

void LALDetectorVel(LALStatus   *status,
		 REAL8       v[3],  /* output velocity vector */
		 LIGOTimeGPS *time0, /* time at which velocity is calculated */
		 LALDetector  detector, /* detector */
		 EphemerisData *edat);

void LALDetectorPos(LALStatus   *status,
		 REAL8       x[3],  /* output position vector */
		 LIGOTimeGPS *time0, /* time at which position is calculated */
		 LALDetector  detector, /* detector*/
		 EphemerisData *edat);

/* ****************************************************** */

/*@}*/

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif  /* end of double inclusion protection */

