/*-----------------------------------------------------------------------
 *
 * File Name: Validation1.h
 *
 * Authors: Sintes, A.M., Krishnan, B.
 *
 * Revision: $Id$
 *
 * History:   Created by Sintes June 16, 2003
 *    to test part of the Hough-Driver code.
 *    No input from SFT data yet implemented here.
 *-----------------------------------------------------------------------
 */
 
/*
 *   Protection against double inclusion (include-loop protection)
 *     Note the naming convention!
 */

#ifndef _VALIDATION1_H
#define _VALIDATION1_H

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

#include "./SFTbin.h"
#include "./Velocity.h"
#include "./PeakSelect.h"

/******************************************************
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif


/******************************************************
 *  Assignment of Id string using NRCSID()
 */

NRCSID (VALIDATION1H, "$Id$");

/******************************************************
 *  Error codes and messages.
 */
 
#define VALIDATION1_ENORM 0
#define VALIDATION1_ESUB  1
#define VALIDATION1_EARG  2
#define VALIDATION1_EBAD  3
#define VALIDATION1_EFILE 4

#define VALIDATION1_MSGENORM "Normal exit"
#define VALIDATION1_MSGESUB  "Subroutine failed"
#define VALIDATION1_MSGEARG  "Error parsing arguments"
#define VALIDATION1_MSGEBAD  "Bad argument values"
#define VALIDATION1_MSGEFILE "Could not create output file"


/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, VALIDATION1H, statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              VALIDATION1H, (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( VALIDATION1_ESUB, VALIDATION1_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return VALIDATION1_ESUB;                                  \
  }                                                                  \
} while (0)
/******************************************************************/

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif


/* ******************************************************************
 *  Structure, enum, union, etc., typdefs.
 */

typedef struct tagREAL8Cart3CoorVector{
  UINT4   	  length; /* number of elements */
  REAL8Cart3Coor  *data; /* x.y.z */
} REAL8Cart3CoorVector;

typedef struct tagLIGOTimeGPSVector{
  UINT4        length; /* number of elements */
  LIGOTimeGPS  *time; /* the collection of times */
} LIGOTimeGPSVector;

/* ****************************************************** */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif


#endif     /* Close double-include protection _VALIDATION1_H */
