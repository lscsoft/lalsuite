
/*----------------------------------------------------------------------- 
 * 
 * File Name: LALTrackSearch.h
 * 
 * Origin: Balasubramanian, R. (Cardiff University, UK)
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
#define MASK_SIZE 4      /* size of the Gaussian mask in units of sigma (std. deviation) */
#define PIXEL_BOUNDARY (0.6)
#define MAX_ANGLE_DIFF (LAL_PI/6.0) /* maximum one sided angle within which to search 
				   for a continuation of the line */
#define LENGTH_THRESHOLD 3  /* A LOWER threshold on the length of the curves in pixels to process*/
#define MAX_CURVE_LENGTH 16384 /* the maximum length allowed */
#define MAX_NUMBER_OF_CURVES 2048 /* the maximum number of curves allowed */


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
/* structures */
typedef struct tagCurve
{
  INT4 n; /* number of points in the curve */
  CHAR junction; /* =1 if the curve has a junction and =0 if no junction */ 
  INT4 *row; /* the row coordinates of the n points */
  INT4 *col; /* the column coordinates of the n points */
  REAL4 *depth; /* the "height" of the pixel in the TF map corresponding 
		  to (col[i],row[i]) */
  REAL4 totalPower; /* resulting numerical intergration along ridge */
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


typedef struct tagTrackSearchOut /* Output Structure for the Track Search algorithm */
{
  INT4 numberOfCurves;           /* the number of curves found*/
  Curve *curves;                 /* a pointer to a array of numberOfCurves curves */
  TrackSearchStore store;       /* a pointer to the temporary storage space */
} TrackSearchOut;

   
typedef struct tagTrackSearchParams /* Parameter structure for the Track search algorithm*/
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
} TrackSearchParams;


/* function Prototypes */

void 
LALSignalTrackSearch(LALStatus *, TrackSearchOut *, const TimeFreqRep *, TrackSearchParams *);

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif






