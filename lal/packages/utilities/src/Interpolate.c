/*----------------------------------------------------------------------- 
 * 
 * File Name: Interpolate.c
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include <string.h>
#include "LALStdlib.h"
#include "Interpolate.h"

NRCSID (INTERPOLATEC, "$Id$");

void
SPolynomialInterpolation (
    Status          *status,
    SInterpolateOut *output,
    REAL4            target,
    SInterpolatePar *params
    )
{
  REAL4 *dn;   /* difference in a step down */
  REAL4 *up;   /* difference in a step up   */
  REAL4  diff;
  UINT4  near;
  UINT4  order;
  UINT4  n;
  UINT4  i;

  INITSTATUS (status, "SPolynomialInterpolation", INTERPOLATEC);

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


void
DPolynomialInterpolation (
    Status          *status,
    DInterpolateOut *output,
    REAL8            target,
    DInterpolatePar *params
    )
{
  REAL8 *dn;   /* difference in a step down */
  REAL8 *up;   /* difference in a step up   */
  REAL8  diff;
  UINT4  near;
  UINT4  order;
  UINT4  n;
  UINT4  i;

  INITSTATUS (status, "DPolynomialInterpolation", INTERPOLATEC);

  ASSERT (output, status, INTERPOLATE_ENULL, INTERPOLATE_MSGENULL);
  ASSERT (params, status, INTERPOLATE_ENULL, INTERPOLATE_MSGENULL);
  ASSERT (params->x, status, INTERPOLATE_ENULL, INTERPOLATE_MSGENULL);
  ASSERT (params->y, status, INTERPOLATE_ENULL, INTERPOLATE_MSGENULL);

  n = params->n;
  ASSERT (n > 1, status, INTERPOLATE_ESIZE, INTERPOLATE_MSGESIZE);

  dn = (REAL8 *) LALMalloc (n*sizeof(REAL8));
  ASSERT (dn, status, INTERPOLATE_ENULL, INTERPOLATE_MSGENULL);

  up = (REAL8 *) LALMalloc (n*sizeof(REAL8));
  ASSERT (up, status, INTERPOLATE_ENULL, INTERPOLATE_MSGENULL);


  /*
   *
   * Initialize dn[] and up[] and find element of domain nearest the target.
   *
   */


  memcpy (dn, params->y, n*sizeof(REAL8));
  memcpy (up, params->y, n*sizeof(REAL8));
  diff = target - params->x[near = 0];
  for (i = 1; diff > 0 && i < n; ++i)
  {
    REAL8 tmp = target - params->x[i];
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
      REAL8 xdn = params->x[i];
      REAL8 xup = params->x[i + order];
      REAL8 den = xdn - xup;
      REAL8 fac;
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
