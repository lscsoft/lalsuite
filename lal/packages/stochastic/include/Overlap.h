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
 * Error codes, typedefs, and protypes for the function overlap().
 *
 * DIAGNOSTICS
 *
 *-----------------------------------------------------------------------
 */

#ifndef _OVERLAP_H
#define _OVERLAP_H

#ifndef _LALSTDLIB_H
#include "LALStdlib.h"
#ifndef _LALSTDLIB_H
#define _LALstdlib_H
#endif
#endif

NRCSID (OVERLAPH, "$Id$");

#define OVERLAP_ENULLP 1
#define OVERLAP_ENULLV 2
#define OVERLAP_ESIZE  4
#define OVERLAP_EDFREQ 8
#define OVERLAP_ESITE  16
#define OVERLAP_ESZMM  32
#define OVERLAP_ENULLD 64

#define OVERLAP_MSGENULLP "Address of parameter structure must be non-null"
#define OVERLAP_MSGENULLV "Address of vector pointer must be non-null"
#define OVERLAP_MSGESIZE  "Length of output vector must be positive"
#define OVERLAP_MSGEDFREQ "Frequency spacing must be positive"
#define OVERLAP_MSGESITE  "Site ID must have legitimate enum type"
#define OVERLAP_MSGESZMM  "Vector lengths must agree"
#define OVERLAP_MSGENULLD "Data area of vector must be non-null"

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

#endif
