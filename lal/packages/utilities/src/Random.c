/*----------------------------------------------------------------------- 
 * 
 * File Name: Random.c
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _TIME_H
#include <time.h>
#ifndef _TIME_H
#define _TIME_H
#endif
#endif

#ifndef _MATH_H
#include <math.h>
#ifndef _MATH_H
#define _MATH_H
#endif
#endif

#ifndef _LALSTDLIB_H
#include "LALStdlib.h"
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _RANDOM_H
#include "Random.h"
#ifndef _RANDOM_H
#define _RANDOM_H
#endif
#endif

NRCSID (RANDOMC, "$Id$");

static const INT4 a = 16807;
static const INT4 m = 2147483647;
static const INT4 q = 127773;
static const INT4 r = 2836;

static const REAL4 eps = 1.2e-7;

static
INT4
BasicRandom (INT4 i)
{
  INT4 k;
  k = i/q;
  i = a*(i - k*q) - r*k;
  if (i < 0)
    i += m;
  return i;
}

void
CreateRandomParams (
    Status        *status,
    RandomParams **params,
    INT4           seed
    )
{
  INT4 n;

  INITSTATUS (status, RANDOMC);

  ASSERT (params, status, RANDOM_ENULL, RANDOM_MSGENULL);
  ASSERT (!*params, status, RANDOM_ENNUL, RANDOM_MSGENNUL);

  *params = (RandomParams *) LALMalloc (sizeof(RandomParams));
  ASSERT (*params, status, RANDOM_ENULL, RANDOM_MSGENULL);

  while (seed == 0)
  {
    seed = time (NULL);
  }

  if (seed < 0)
  {
    seed = -seed;
  }

  (*params)->i = seed;

  for (n = 0; n < 8; ++n)
  {
    (*params)->i = BasicRandom ((*params)->i);
  }

  for (n = 0; n < sizeof((*params)->v)/sizeof(INT4); ++n)
  {
    (*params)->v[n] = (*params)->i = BasicRandom ((*params)->i);
  }

  (*params)->y = (*params)->v[0];

  RETURN (status);
}

void
DestroyRandomParams (
    Status        *status,
    RandomParams **params
    )
{
  INITSTATUS (status, RANDOMC);

  ASSERT (params, status, RANDOM_ENULL, RANDOM_MSGENULL);
  ASSERT (*params, status, RANDOM_ENULL, RANDOM_MSGENULL);

  LALFree (*params);
  *params = NULL;

  RETURN (status);
}


void
UniformDeviate (
    Status       *status,
    REAL4        *deviate,
    RandomParams *params
    )
{
  INT4 ndiv;
  INT4 n;

  INITSTATUS (status, RANDOMC);

  ASSERT (deviate, status, RANDOM_ENULL, RANDOM_MSGENULL);
  ASSERT (params, status, RANDOM_ENULL, RANDOM_MSGENULL);

  ndiv = 1 + (m - 1)/(sizeof(params->v)/sizeof(INT4));
  n    = params->y/ndiv;

  params->y = params->v[n];
  params->v[n] = params->i = BasicRandom (params->i);

  *deviate = params->y/(REAL4)m;
  if (*deviate > 1 - eps)
  {
    *deviate = 1 - eps;
  }

  RETURN (status);
}


void
NormalDeviates (
    Status       *status,
    REAL4Vector  *deviates,
    RandomParams *params
    )
{
  REAL4 *data;
  INT4   half;

  INITSTATUS (status, RANDOMC);
  ATTATCHSTATUSPTR (status);

  ASSERT (params, status, RANDOM_ENULL, RANDOM_MSGENULL);
  ASSERT (deviates, status, RANDOM_ENULL, RANDOM_MSGENULL);
  ASSERT (deviates->data, status, RANDOM_ENULL, RANDOM_MSGENULL);
  ASSERT (deviates->length > 0, status, RANDOM_ESIZE, RANDOM_MSGESIZE);

  data = deviates->data;
  half = deviates->length/2;

  while (half-- > 0)
  {
    REAL4 u;
    REAL4 v;
    REAL4 x;
    REAL4 y;
    REAL4 rsq;
    REAL4 fac;

    do {
      UniformDeviate (status->statusPtr, &u, params);
      CHECKSTATUSPTR (status);
      UniformDeviate (status->statusPtr, &v, params);
      CHECKSTATUSPTR (status);
      x   = 2*u - 1;
      y   = 2*v - 1;
      rsq = x*x + y*y;
    }
    while (rsq >= 1 || rsq == 0);

    fac     = sqrt(-2*log(rsq)/rsq);
    *data++ = fac*x;
    *data++ = fac*y;
  }

  /* do it again if there is an odd amount of data */
  if (deviates->length % 2)
  {
    REAL4 u;
    REAL4 v;
    REAL4 x;
    REAL4 y;
    REAL4 rsq;
    REAL4 fac;

    do {
      UniformDeviate (status->statusPtr, &u, params);
      CHECKSTATUSPTR (status);
      UniformDeviate (status->statusPtr, &v, params);
      CHECKSTATUSPTR (status);
      x   = 2*u - 1;
      y   = 2*v - 1;
      rsq = x*x + y*y;
    }
    while (rsq >= 1 || rsq == 0);

    fac   = sqrt(-2*log(rsq)/rsq);
    *data = fac*x;
    /* throw away y */
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}
