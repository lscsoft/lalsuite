/*-----------------------------------------------------------------------
 *
 * File Name: LALDatatypes.h
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#ifndef _LALATOMICDATATYPES_H
#define _LALATOMICDATATYPES_H

#include "LALConfig.h"
#include "LALRCSID.h"

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( LALATOMICDATATYPESH, "$Id$" );

typedef char CHAR;
typedef unsigned char UCHAR;
typedef unsigned char BOOLEAN;

/* Integer types */

#if SIZEOF_SHORT == 2
  typedef short INT2;
  typedef unsigned short UINT2;
#elif SIZEOF_INT == 2
  typedef int INT2;
  typedef unsigned int UINT2;
#else
# error "ERROR: NO 2 BYTE INTEGER FOUND"
#endif

#if SIZEOF_INT == 4
  typedef int INT4;
  typedef unsigned int UINT4;
#elif SIZEOF_LONG == 4
  typedef long INT4;
  typedef unsigned long UINT4;
#else
# error "ERROR: NO 4 BYTE INTEGER FOUND"
#endif

#if SIZEOF_LONG == 8
  typedef long INT8;
  typedef unsigned long UINT8;
#elif SIZEOF_LONG_LONG == 8
  typedef long long INT8;
  typedef unsigned long long UINT8;
#else
# error "ERROR: NO 8 BYTE INTEGER FOUND"
#endif

/* Real types */

#if SIZEOF_FLOAT == 4
  typedef float REAL4;
#else
# error "ERROR: NO 4 BYTE REAL FOUND"
#endif

#if SIZEOF_DOUBLE == 8
  typedef double REAL8;
#else
# error "ERROR: NO 8 BYTE REAL FOUND"
#endif

/* Complex types */

typedef struct
tagCOMPLEX8
{
  REAL4 re;
  REAL4 im;
}
COMPLEX8;

typedef struct
tagCOMPLEX16
{
  REAL8 re;
  REAL8 im;
}
COMPLEX16;

#ifdef  __cplusplus
}
#endif

#endif /* _LALATOMICDATATYPES_H */
