/*-----------------------------------------------------------------------
 *
 * File Name: Overlap.h
 *
 * Author: J.D. Romano
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * Overlap.h
 *
 * SYNOPSIS
 * #include "Overlap.h"
 *
 * DESCRIPTION
 * Error codes, typedefs, and protypes for the function Overlap().
 *
 * DIAGNOSTICS
 *
 *-----------------------------------------------------------------------
 */

#ifndef _OVERLAP_H
#define _OVERLAP_H

#include "LALStdlib.h"

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (OVERLAPH, "$Id$");

#define OVERLAP_ENULLIP    1
#define OVERLAP_ESITE      2
#define OVERLAP_ESIZE      4
#define OVERLAP_EDELTAF    8
#define OVERLAP_ENULLOP    16
#define OVERLAP_ESIZEMM    32
#define OVERLAP_ENULLD     64

#define OVERLAP_MSGENULLIP "Pointer to input parameters must be non-null"
#define OVERLAP_MSGESITE   "Site ID must have legitimate enum type"
#define OVERLAP_MSGESIZE   "Specified length of output vector must be > 0"
#define OVERLAP_MSGEDELTAF "Frequency spacing must be > 0"
#define OVERLAP_MSGENULLOP "Pointer to output vector must be non-null"
#define OVERLAP_MSGESIZEMM "Length of output vector must agree with length specified in input parameters"
#define OVERLAP_MSGENULLD  "Pointer to data member of output vector must be non-null"

#define C_LIGHT (2.99792458e10) /* speed of light [cm/sec] */

typedef enum {
  LHO,      /* LIGO Hanford Observatory    */
  LLO,      /* LIGO Livingston Observatory */
  VIRGO,
  GEO600,
  TAMA300
}
IFOsite;

#define SITENAMELIST  { "LHO", "LLO", "VIRGO", "GEO600", "TAMA300" }
#define NUMBEROFSITES ((int)(TAMA300 + 1))

typedef struct tagSiteParameters {
  IFOsite   siteID;
  CHAR*     siteName;
  REAL4     vertexNorth;
  REAL4     vertexWest;
  REAL4     arm1Orientation;
  REAL4     arm2Orientation;
  REAL4     armLengthCM;
}
SiteParameters;

typedef struct tagSiteCoordinates {
  REAL4     vertex[3];
  REAL4     arm1[3];
  REAL4     arm2[3];
}
SiteCoordinates;

typedef struct tagOverlapParameters {
  IFOsite   site1ID;  /* ID number for interferometer site #1 */
  IFOsite   site2ID;  /* ID number for interferometer site #2 */
  INT4      length;   /* length of vector containing overlap red. function */
  REAL4     deltaF;   /* frequency spacing for overlap reduction function */
}
OverlapParameters;

void
Overlap ( Status *, REAL4Vector *, OverlapParameters * );


#ifdef  __cplusplus
}
#endif

#endif /* _OVERLAP_H */
