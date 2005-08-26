/*  
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
 * 
 */



/**
 * DriveHoughFStat.h 
 * \author Badri Krishnan
 * Date : August 2005
 * \brief
 * Header file for DriveHoughFStat.c 
 * 
 ****/




#ifndef _DRIVEHOUGHCOLOR_H
#define _DRIVEHOUGHCOLOR_H

/* standard includes */
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glob.h>
#include <time.h>
#include <errno.h> 

/* lal includes */
#include <lal/UserInput.h>
#include <lal/LALStdlib.h>
#include <lal/PulsarDataTypes.h>
#include <lal/SFTfileIO.h>
#include <lal/AVFactories.h>
#include <lal/RngMedBias.h>
#include <lal/LALComputeAM.h>
#include <lal/ComputeSky.h>
#include <lal/LALInitBarycenter.h>
#include <lal/Velocity.h>
#include <lal/LALDemod.h>
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/Date.h>
#include <lal/LUT.h> 

#include <lalapps.h>

/******************************************************
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif


/******************************************************
 *  Assignment of Id string using NRCSID()
 */

NRCSID( DRIVEHOUGHFSTATH, "$Id$" );

/******************************************************
 *  Error codes and messages.
 */
 
#define DRIVEHOUGHFSTAT_ENORM 0
#define DRIVEHOUGHFSTAT_ESUB  1
#define DRIVEHOUGHFSTAT_EARG  2
#define DRIVEHOUGHFSTAT_EBAD  3
#define DRIVEHOUGHFSTAT_EFILE 4
#define DRIVEHOUGHFSTAT_ENULL 5

#define DRIVEHOUGHFSTAT_MSGENORM "Normal exit"
#define DRIVEHOUGHFSTAT_MSGESUB  "Subroutine failed"
#define DRIVEHOUGHFSTAT_MSGEARG  "Error parsing arguments"
#define DRIVEHOUGHFSTAT_MSGEBAD  "Bad argument values"
#define DRIVEHOUGHFSTAT_MSGEFILE "Could not create output file"
#define DRIVEHOUGHFSTAT_MSGENULL "Null pointer"



/* ******************************************************************
 *  Structure, enum, union, etc., typdefs.
 */

/** structure containing Fstat values in a frequency range for a single alpha, delta, fdot */
  typedef struct tagFstatVector {
    REAL8 alpha;   /**< right ascension */
    REAL8 delta;   /**< declination */
    INT4  msp;     /**< largest spindown order */
    REAL8 *fdot;   /**< vector of length msp -- first, second and higher spindowns */
    INT8  length;  /**< number of frequency values */
    REAL8 fStart;  /**< starting frequency in Hz */
    REAL8 deltaF;  /**< frequency resolution in Hz*/
    REAL8 fBand;   /**< frequency band in Hz */
    REAL8 *data;   /**< Fstat values for the corresponding frequency values */ 
  } FstatVector;


  /** structure containing dtector velocity and position for set of timestamps */
  typedef struct tagTimeVelPosVector {
    INT4 length;     /**< number of time stamps */
    LIGOTimeGPS *ts; /**< time stamps */
    REAL8 *velx;     /**< x-component of velocity in equatorial coordinates */
    REAL8 *vely;     /**< y-component of velocity in equatorial coordinates */
    REAL8 *velz;     /**< z-component of velocity in equatorial coordinates */
    REAL8 *posx;     /**< x-component of position in equatorial coordinates */
    REAL8 *posy;     /**< y-component of position in equatorial coordinates */
    REAL8 *posz;     /**< z-component of position in equatorial coordinates */
  } TimeVelPosVector; 
    
#ifdef  __cplusplus
}                /* Close C++ protection */
#endif


#endif     /* Close double-include protection _DRIVEHOUGHFSTAT_H */
