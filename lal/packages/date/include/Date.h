/*-----------------------------------------------------------------------
 *
 * File Name: Date.h
 *
 * Author: David Chin <dwchin@umich.edu>
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * Date.h
 *
 * SYNOPSIS
 * #include "Date.h"
 *
 * DESCRIPTION
 * Data type definitions for date and time manipulation.
 * Function prototypes for date and time manipulation functions.
 *
 * DIAGNOSTICS
 *
 *-----------------------------------------------------------------------
 */

#ifndef _DATE_H
#define _DATE_H

#ifndef _REENTRANT
#   define _REENTRANT
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "LALRCSID.h"

#include "LALConstants.h"
#include "LALAtomicDatatypes.h"
#include "LALDatatypes.h"
#include "LALStatusMacros.h"
#include "LALStdlib.h"

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (DATEH, "$Id$");

/*
 * Julian.c
 */
#define JULIAN_ENULLINPUT    1
#define JULIAN_ENULLOUTPUT   2
#define JULIAN_EDATETOOEARLY 3

#define JULIAN_MSGENULLINPUT "Input is NULL"
#define JULIAN_MSGENULLOUTPUT "Output is NULL"
#define JULIAN_MSGEDATETOOEARLY "Date too early: Julian Day can only be computed for dates >= 1900-Mar"

/*
 * UtoGPS.c
 */
#define UTOGPS_ENULLINPUT   1
#define UTOGPS_ENULLOUTPUT  2
    
#define UTOGPS_MSGENULLINPUT "Input is NULL"
#define UTOGPS_MSGENULLOUTPUT "Output is NULL"
    
/*
 * Utime.c
 */
#define UTIME_ENULLINPUT  1
#define UTIME_ENULLOUTPUT 2
#define UTIME_ERANGE      3

#define UTIME_MSGENULLINPUT "Input is NULL"
#define UTIME_MSGENULLOUTPUT "Output is NULL"
#define UTIME_MSGERANGE "Input time out of range: 0 <= utc_seconds <= 946684823"

/*
 * DateString.c
 */
#define DATESTRING_ENULLINPUT    1
#define DATESTRING_ENULLOUTPUT   2
#define DATESTRING_EBUFFTOOSMALL 3

#define DATESTRING_MSGENULLINPUT "Input is NULL"
#define DATESTRING_MSGENULLOUTPUT "Output is NULL"
#define DATESTRING_MSGEBUFFTOOSMALL "Output timestamp string too small: min. size = 26"
    
/*
 * LMST1.c
 */
#define LMST1_ENULLINPUT  1
#define LMST1_ENULLOUTPUT 2
    
#define GMST1_ENULLINPUT  1
#define GMST1_ENULLOUTPUT 2
    
#define LMST1_MSGENULLINPUT "Input is NULL"
#define LMST1_MSGENULLOUTPUT "Output is NULL"
    
#define GMST1_MSGENULLINPUT "Input is NULL"
#define GMST1_MSGENULLOUTPUT "Output is NULL"

#define MST_SEC 0
#define MST_HRS 1
#define MST_DEG 2
#define MST_RAD 3

    
/*
 * SecsToLALDate.c
 */
#define SECSTOLALDATE_ENULLOUTPUT 1
    
#define SECSTOLALDATE_MSGENULLOUTPUT "Output is NULL"

    
/* The standard Unix tm structure */
typedef struct
tm
LALUnixDate;

/* This time object is exactly like LIGOTimeGPS, except for the name */
typedef struct
tagLIGOTimeUnix
{
    INT4 unixSeconds;
    INT4 unixNanoSeconds;
}
LIGOTimeUnix;

/* Encode timezone information */
typedef struct
tagLALTimezone
{
    INT4 secondsWest; /* seconds West of UTC */
    INT4 dst;         /* Daylight Savings Time correction to apply */
}
LALTimezone;    

/* Date and time structure */
typedef struct
tagLALDate
{
    LALUnixDate unixDate;
    INT4        residualNanoSeconds; /* residual nanoseconds */
    LALTimezone timezone;            /* timezone information */
}
LALDate;

/*
 * 9. Functions Declarations (i.e., prototypes).
 */
void LALJulianDay (LALStatus*, INT4*, const LALDate*);
void LALModJulianDay(LALStatus*, REAL8*, const LALDate*);
void LALJulianDate(LALStatus*, REAL8*, const LALDate*);
void LALModJulianDate(LALStatus*, REAL8*, const LALDate*);
void LALUtoGPS(LALStatus*, LIGOTimeGPS*, const LIGOTimeUnix*);
void LALGPStoU(LALStatus*, LIGOTimeUnix*, const LIGOTimeGPS*);
void LALUtime(LALStatus*, LALDate*, const LIGOTimeUnix*);
void LALDateString(LALStatus*, CHARVector*, const LALDate*);
void LALGMST1(LALStatus*, REAL8*, const LALDate*, INT4);
void LALLMST1(LALStatus*, REAL8*, const LALDate*, REAL8, INT4);
void LALSecsToLALDate(LALStatus*, LALDate*, REAL8);


#ifdef  __cplusplus
}
#endif

#endif /* _DATE_H */
