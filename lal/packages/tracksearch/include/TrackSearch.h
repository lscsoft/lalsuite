/*
*  Copyright (C) 2007 Cristina Valeria Torres, Jolien Creighton
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

/*-----------------------------------------------------------------------
 *
 * File Name: LALTrackSearch.h
 *
 * Current Developer: Torres, Cristina V.  (LLO)
 * Original Developer: Balasubramanian, R. (Cardiff University, UK)
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------

 *
 * NAME
 * TrackSearch.h
 *
 * SYNOPSIS
 * #include <lal/TrackSearch.h>
 *
 * DESCRIPTION
 * Defines Structures function prototypes and macros for use by the
 * curve tracking algorithm to find curves in time frequency maps.
 * This header file and the other source files for curve tracking are adapted
 * from Carsten Steger's detect-lines package.
 *
 * DIAGNOSTICS
 *
 *-----------------------------------------------------------------------
 */

#ifndef _TRACKSEARCH_H
#define _TRACKSEARCH_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/TimeFreq.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID (TRACKSEARCHH, "$Id$");

/* Mathematical constants */
#ifndef SQRT_2_PI_INV
#define SQRT_2_PI_INV 0.39894228040143272L  /* (1.0 / sqrt(2.0 * pi))  */
#endif
#define MAX_ANGLE_DIFFERENCE (LAL_PI/6.0) /* the max angle differene between 2 adjacent line points */
  /*
   * Need to study if we should loosen this angle differenc
   * requirement to try and building longer lines that are less
   * straight so that chirps are easier to follow!
   * Currently 30Deg consider considering relax to 45Deg
   */
#define MASK_SIZE 4      /* size of the Gaussian mask in units of sigma (std. deviation) */
#define PIXEL_BOUNDARY (0.6)
<<<<<<< HEAD:lal/packages/tracksearch/include/TrackSearch.h
#define MAX_ANGLE_DIFF (LAL_PI/6.0) /* maximum one sided angle within which to search 
=======
  /* Permanent Swap out of below Wed-Dec-17-2008:200812171543 Used to
     use 4.0*/
#define MAX_ANGLE_DIFF (LAL_PI/6.0) /* maximum one sided angle within which to search
>>>>>>> master:lal/packages/tracksearch/include/TrackSearch.h
				   for a continuation of the line */
#define LENGTH_THRESHOLD 3  /* A LOWER threshold on the length of the curves in pixels to process*/
#define MAX_CURVE_LENGTH 16384 /* the maximum length allowed */
#define MAX_NUMBER_OF_CURVES 65536 /* the maximum number of curves allowed */


/* Error codes and strings */

#define TS_NULL_POINTER 1
#define TS_NON_NULL_POINTER 2
#define TS_ALLOC 3
#define TS_UNITIALIZED_MASKS 4
#define TS_MAP_DIMENSION 5
#define TS_ILLEGAL_PARAMS 6
#define TS_LINE_START 7
#define TS_ARRAY_OVERFLOW 8
#define TS_TOO_MANY_CURVES 9
#define TS_OCTANT 10
#define TS_MASKSIZE 11
#define TS_SIGMASIZE 12
#define TS_ALLOC_ERROR 13

#define TS_MSGNON_NULL_POINTER  "Allocation of memory to Non null pointer"
#define TS_MSG_NULL_POINTER     "Null pointer passed when non NUll pointer expected"
#define TS_MSG_ALLOC            "flag for allocating storage space does not contain a valid value (0,1,2)"
#define TS_MSG_UNINITIALIZED_MASKS "Gaussian masks appear not to have been initialized"
#define TS_MSG_MAP_DIMENSION     "Dimensions of Map have to be each greater than 8"
#define TS_MSG_ILLEGAL_PARAMS   "Illegal parameter values ( high and low should be positive and high>low)"
#define TS_MSG_LINE_START "The number of line start points is not consistent with previous module"
#define MSG_TS_ARRAY_OVERFLOW  "Array bounds can be crossed "
#define MSG_TS_TOO_MANY_CURVES " Too many curves found in map "
#define MSG_TS_OCTANT "octant cannot be greater than 4 "
#define MSG_TS_MASKSIZE "maskSize too small"
#define MSG_TS_SIGMASIZE "sigma value to small to be useful"
#define MSG_TS_ALLOC_ERROR "error trying to allocate or deallocate memory"
/* structures */

typedef struct tagCurve
{
  INT4 n; /* number of points in the curve */
  CHAR junction; /* =1 if the curve has a junction */
  CHAR trash; /* Linked list postprocessing K-keep D-drop */
  INT4 *row; /* the row coordinates of the n points */
  INT4 *col; /* the column coordinates of the n points */
  LIGOTimeGPS *gpsStamp;  /* real gps timestate variable*/
  REAL4 *fBinHz; /* frequency value not bin label */
  REAL4 *depth; /* the "height" of the pixel in the TF map corresponding
		  to (col[i],row[i]) */
  REAL4 totalPower; /* resulting numerical intergration along ridge */
  REAL4 snrEstimate; /*estimate of trigger SNR */
} Curve;

typedef struct tagTrackSearchStore  /* Structure for storage space for the algorithm */
{
  INT4 height;            /* height of the TF map (increasing frequency direction) */
  INT4 width;             /* width of the TF map (increasing time direction ) */
  CHAR **isLine;          /* 2D map used to flag potential line points */
  REAL4 **k[5];           /* arrays used to store the image convolved with the first
			     and second derivatives of the Gaussian kernel */
  REAL4 **imageTemp;      /* A temporary 2D array to contain intermediate transformed images;
			   Also used to store the maximum eigenValue*/
  REAL8 *gaussMask[3];    /* The Gaussian mask and its first and second derivatives. */
  REAL4 *eigenVec;        /* Array to store vector in maximal direction each point at each pixel point*/
  INT4  numLStartPoints; /* a variable which contains the number of possible starting points for curves */
  INT4  numLPoints;      /* A variable which contains the number of possible line points */
} TrackSearchStore;

/*
 * Enum type for the possible operations for curve cutting
 */
typedef enum tagTrackSearchCut
  {
    none, PandL, nPandnL, PandnL, nPandL
  }TrackSearchCut;

/*
 * Output structure for the tracksearch algorithm
 */
typedef struct tagTrackSearchOut
{
  INT4 numberOfCurves;          /* the number of curves found*/
  Curve *curves;                /* a pointer to a array of numberOfCurves curves */
  TrackSearchStore store;       /* a pointer to the temporary storage space */
  REAL4 minPowerCut;            /* copy of applied threshold value */
  REAL4 minSNRCut;              /* copy of applied min SNR value */
  REAL4 minLengthCut;           /* copy of applied threshold value */
  REAL4 startThreshCut;         /* copy of applied Lh threshold */
  REAL4 linePThreshCut;         /* copy of applied Ll threshold */
  TrackSearchCut thresholdLogic;/* Logic operation to apply with mPC and mLC*/
} TrackSearchOut;

/*
 * Structure to help determine the correlation of integer map indices
 * and real gps times and frequencies
 */
typedef struct tagTrackSearchMapMarkingParams
{
  REAL4 deltaT;/*The time bin resolution dur/#bins */
  REAL4 dataDeltaT;/*Sampling rate of data used*/
  LIGOTimeGPS mapStartGPS;
  LIGOTimeGPS mapStopGPS;
  INT4 mapTimeBins;
  INT4 mapFreqBins;
} TrackSearchMapMarkingParams;

/*
 * Structure with params specific to the LAL library function
 */
typedef struct tagTrackSearchParams
{
  INT4  height;   /* height of Map */
  INT4  width;    /* width of map */
  INT4  allocFlag;/* flag for allocating/deallocating the temporary storage structure
		  flag must be = 1 the first time the function is called
		  flag=0 -- no allocation of space
		  flag=1 -- allocate memory and initialize store
		  flag=2 -- free store and exit the function
		  */
  REAL4 sigma; /* the width of the lines to search for (in pixels) */
  REAL4 high;  /* the upper threshold on the 2nd derivative (along maximal direction)*/
  REAL4 low;   /* the lower threshold on the 2nd derivative (along maximal direction)*/
  BOOLEAN autoL; /* Bool flag to call auto Lambda code if needed */
  REAL4 rLow;  /* the lower threshold relative to upper threshold
		*  (only used if using auto threshold adjustments)
		*/
  REAL4 autoLambda; /* given Z score value determine best auto lambda for run */
} TrackSearchParams;

/*
 * Function Prototypes
 */
void
LALSignalTrackSearch(LALStatus *,
		     TrackSearchOut *,
		     const TimeFreqRep *,
		     TrackSearchParams *);

void
LALTrackSearchInsertMarkers(LALStatus *,
			    TrackSearchOut *,
			    TrackSearchMapMarkingParams *);


#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif






