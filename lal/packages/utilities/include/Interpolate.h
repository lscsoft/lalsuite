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

#include "LALDatatypes.h"

#ifdef __cplusplus
extern "C" {
#endif


NRCSID (INTERPOLATEH, "$Id: Interpolate.h");

#define INTERPOLATE_ENULL 1
#define INTERPOLATE_ESIZE 2
#define INTERPOLATE_EZERO 4

#define INTERPOLATE_MSGENULL "Null pointer"
#define INTERPOLATE_MSGESIZE "Invalid size" 
#define INTERPOLATE_MSGEZERO "Zero divide"

typedef struct
tagSInterpolateOut
{
  REAL4  y;
  REAL4 dy;
}
SInterpolateOut;

typedef struct
tagDInterpolateOut
{
  REAL8  y;
  REAL8 dy;
}
DInterpolateOut;

typedef struct
tagSInterpolatePar
{
  UINT4  n;
  REAL4 *x;
  REAL4 *y;
}
SInterpolatePar;

typedef struct
tagDInterpolatePar
{
  UINT4  n;
  REAL8 *x;
  REAL8 *y;
}
DInterpolatePar;

void
LALSPolynomialInterpolation (
    LALStatus          *status,
    SInterpolateOut *output,
    REAL4            target,
    SInterpolatePar *params
    );

void
LALDPolynomialInterpolation (
    LALStatus          *status,
    DInterpolateOut *output,
    REAL8            target,
    DInterpolatePar *params
    );


#ifdef __cplusplus
}
#endif

#endif /* _INTERPOLATE_H */
