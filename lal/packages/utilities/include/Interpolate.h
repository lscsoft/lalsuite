/*----------------------------------------------------------------------- 
 * 
 * File Name: Interpolate.h
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _INTERPOLATE_H
#define _INTERPOLATE_H

#ifndef _LALDATATYPES_H
#include "LALDatatypes.h"
#ifndef _LALDATATYPES_H
#define _LALDATATYPES_H
#endif
#endif

NRCSID (INTERPOLATEH, "$Id: Interpolate.h");

#define INTERPOLATE_ENULL 1
#define INTERPOLATE_ESIZE 2
#define INTERPOLATE_EZERO 4

#define INTERPOLATE_MSGENULL "Null pointer"
#define INTERPOLATE_MSGESIZE "Invalid size" 
#define INTERPOLATE_MSGEZERO "Zero divide"

typedef struct
tagInterpolateOut
{
  REAL4  y;
  REAL4 dy;
}
InterpolateOut;

typedef struct
tagInterpolatePar
{
  INT4   n;
  REAL4 *x;
  REAL4 *y;
}
InterpolatePar;

void
PolynomialInterpolation (
    Status         *status,
    InterpolateOut *output,
    REAL4           target,
    InterpolatePar *params
    );

#endif /* _INTERPOLATE_H */
