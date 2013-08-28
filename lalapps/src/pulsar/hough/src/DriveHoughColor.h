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

/**
 * \file DriveHoughColor.h
 * \author Alicia Sintes, Badri Krishnan
 * \brief Header file for non-demodulated Hough search
 *
 * File Name: DRIVEHOUGHCOLOR.h
 *
 * Authors: Sintes, A.M., Krishnan, B.
 *
 * History:   Created by Sintes June 16, 2003
 * to test part of the Hough-Driver code.
 *
 * -----------------------------------------------------------------------
 */
 
/*
 *   Protection against double inclusion (include-loop protection)
 *     Note the naming convention!
 */

#ifndef _DRIVEHOUGHCOLOR_H
#define _DRIVEHOUGHCOLOR_H

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glob.h>
#include <time.h>
#include <errno.h> 

#include <lal/Date.h>
#include <lal/DetectorSite.h>
#include <lal/LALDatatypes.h>
#include <lal/LALHough.h>
#include <lal/RngMedBias.h>
#include <lal/LALRunningMedian.h>
#include <lal/Velocity.h>
#include <lal/Statistics.h>
#include <lal/ComputeFstat.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/LALInitBarycenter.h>
#include <lal/SFTClean.h>
#include <lalapps.h>
#include <gsl/gsl_cdf.h>

#include "./PeakSelect.h"

/******************************************************
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif


/******************************************************
 *  Error codes and messages.
 */
 
#define DRIVEHOUGHCOLOR_ENORM 0
#define DRIVEHOUGHCOLOR_ESUB  1
#define DRIVEHOUGHCOLOR_EARG  2
#define DRIVEHOUGHCOLOR_EBAD  3
#define DRIVEHOUGHCOLOR_EFILE 4
#define DRIVEHOUGHCOLOR_EDIR 4
#define DRIVEHOUGHCOLOR_ENULL 5
#define DRIVEHOUGHCOLOR_ENONULL 5

#define DRIVEHOUGHCOLOR_MSGENORM "Normal exit"
#define DRIVEHOUGHCOLOR_MSGESUB  "Subroutine failed"
#define DRIVEHOUGHCOLOR_MSGEARG  "Error parsing arguments"
#define DRIVEHOUGHCOLOR_MSGEBAD  "Bad argument values"
#define DRIVEHOUGHCOLOR_MSGEFILE "Could not create output file"
#define DRIVEHOUGHCOLOR_MSGEDIR  "Could not create directory"
#define DRIVEHOUGHCOLOR_MSGENULL "Null pointer"
#define DRIVEHOUGHCOLOR_MSGENONULL "Not a Null pointer"

/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    XLALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, "$Id$", statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    XLALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              "$Id$", (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( DRIVEHOUGHCOLOR_ESUB, DRIVEHOUGHCOLOR_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return DRIVEHOUGHCOLOR_ESUB;                                  \
  }                                                                  \
} while (0)
/******************************************************************/

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

#define PIXELFACTOR  2 


/* ******************************************************************
 *  Structure, enum, union, etc., typdefs.
 */

  typedef struct tagREAL8Cart3CoorVector{
    UINT4   	  length; /**< number of elements */
    REAL8Cart3Coor  *data; /**< x.y.z */
  } REAL8Cart3CoorVector;
  
  typedef struct tagHoughSignificantEvent{
    REAL8  nStar;  /**< most significant number count in a skypatch*/
    REAL8  nStarSignificance; /**< significance of number count nStar */
    REAL8  freqStar;  /**< frequency of nStar */
    REAL8  alphaStar; /**< right-ascension of nStar */
    REAL8  deltaStar; /**< declination of nStar */
    REAL8  fdotStar; /**< value of first spindown parameter */
  } HoughSignificantEvent;

  typedef struct tagHoughSignificantEventVector{
    INT4  length;
    HoughSignificantEvent *event;
  } HoughSignificantEventVector;


  typedef struct tagHoughSkyPatchesInfo{
    UINT4 numSkyPatches;
    REAL8 *alpha;
    REAL8 *delta;
    REAL8 *alphaSize;
    REAL8 *deltaSize;
  } HoughSkyPatchesInfo;
    
/**
 * struct fo storing all the variables affected by the
 * selection of a subset of SFTs
 */
  typedef struct tagBestVariables{
    UINT4 length;   /**< the number of SFTs to be selected */
    REAL8Vector *weightsV; /**< noise and AM weights */
    REAL8Vector *timeDiffV; /**< the vector of time diffs */
    REAL8Cart3CoorVector *velV; /**< vector of detector velocities */
    HOUGHPeakGramVector *pgV; /**< the vector of peakgrams */
  } BestVariables;


/*
 *  Functions Declarations (i.e., prototypes).
 */

void Periodo2PSDrng (LALStatus  *status,
		     REAL8Periodogram1    *psd,
		     REAL8Periodogram1    *peri,
		     UINT2                *blocksRNG);


/* ****************************************************** */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif


#endif     /* Close double-include protection _DRIVEHOUGHCOLOR_H */
