/*----------------------------------------------------------------------- 
 * 
 * File Name: Interpolate.c
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _MATH_H
#include <math.h>
#ifndef _MATH_H
#define _MATH_H
#endif
#endif

#ifndef _STRING_H
#include <string.h>
#ifndef _STRING_H
#define _STRING_H
#endif
#endif

#ifndef _LALSTDLIB_H
#include "LALStdlib.h"
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _INTERPOLATE_H
#include "Interpolate.h"
#ifndef _INTERPOLATE_H
#define _INTERPOLATE_H
#endif
#endif


NRCSID (INTERPOLATEC, "$Id$");

void
PolynomialInterpolation (
    Status         *status,
    InterpolateOut *output,
    REAL4           target,
    InterpolatePar *params
    )
{
  REAL4 *dn;   /* difference in a step down */
  REAL4 *up;   /* difference in a step up   */
  REAL4  diff;
  INT4   near;
  INT4   order;
  INT4   n;
  INT4   i;

  INITSTATUS (status, INTERPOLATEC);

  ASSERT (output, status, INTERPOLATE_ENULL, INTERPOLATE_MSGENULL);
  ASSERT (params, status, INTERPOLATE_ENULL, INTERPOLATE_MSGENULL);
  ASSERT (params->x, status, INTERPOLATE_ENULL, INTERPOLATE_MSGENULL);
  ASSERT (params->y, status, INTERPOLATE_ENULL, INTERPOLATE_MSGENULL);

  n = params->n;
  ASSERT (n > 1, status, INTERPOLATE_ESIZE, INTERPOLATE_MSGESIZE);

  dn = (REAL4 *) LALMalloc (n*sizeof(REAL4));
  ASSERT (dn, status, INTERPOLATE_ENULL, INTERPOLATE_MSGENULL);

  up = (REAL4 *) LALMalloc (n*sizeof(REAL4));
  ASSERT (up, status, INTERPOLATE_ENULL, INTERPOLATE_MSGENULL);


  /*
   *
   * Initialize dn[] and up[] and find element of domain nearest the target.
   *
   */


  memcpy (dn, params->y, n*sizeof(REAL4));
  memcpy (up, params->y, n*sizeof(REAL4));
  diff = target - params->x[near = 0];
  for (i = 1; diff > 0 && i < n; ++i)
  {
    REAL4 tmp = target - params->x[i];
    diff = (fabs(tmp) < fabs(diff) ? near = i, tmp : diff);
  }
  output->y = params->y[near];


  /*
   *
   * Recompute dn[] and up[] for each polynomial order.
   *
   */


  for (order = 1; order < n; ++order)
  {
    INT4 imax = n - order;
    for (i = 0; i < imax; ++i)
    {
      REAL4 xdn = params->x[i];
      REAL4 xup = params->x[i + order];
      REAL4 den = xdn - xup;
      REAL4 fac;
      ASSERT (den != 0, status, INTERPOLATE_EZERO, INTERPOLATE_MSGEZERO);
      fac   = (dn[i + 1] - up[i])/den;
      dn[i] = fac*(xdn - target);
      up[i] = fac*(xup - target);
    }

    /* go down unless impossible */
    output->y += (near < imax ? dn[near] : up[--near]);
  }

  output->dy = fabs(dn[0]) < fabs(up[0]) ? fabs(dn[0]) : fabs(up[0]);

  LALFree (dn);
  LALFree (up);
  RETURN  (status);
}
