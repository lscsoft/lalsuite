/*
*  Copyright (C) 2007 Alexander Dietz, Jolien Creighton, Kipp Cannon
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <time.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/Random.h>
#include <lal/Sequence.h>
#include <lal/XLALError.h>

/**
 * \defgroup Random_c Module Random.c
 * \ingroup Random_h
 *
 * \brief Functions for generating random numbers.
 *
 * \heading{Description}
 *
 * The routines <tt>LALCreateRandomParams()</tt> and <tt>LALDestroyRandomParams()</tt>
 * create and destroy a parameter structure for the generation of random
 * variables.  The creation routine requires a random number seed \c seed.
 * If the seed is zero then a seed is generated using the current time.
 *
 * The routine <tt>LALUniformDeviate()</tt> returns a single random deviate
 * distributed uniformly between zero and unity.
 *
 * The routine <tt>LALNormalDeviates()</tt> fills a vector with normal (Gaussian)
 * deviates with zero mean and unit variance, whereas the function\c XLALNormalDeviate just returns one normal distributed random number.
 *
 * \heading{Operating Instructions}
 *
 * \code
 * static LALStatus     status;
 * static RandomParams *params;
 * static REAL4Vector  *vector;
 * UINT4 i;
 * INT4 seed = 0;
 *
 * LALCreateVector( &status, &vector, 9999 );
 * LALCreateRandomParams( &status, &params, seed );
 *
 * /\* fill vector with uniform deviates *\/
 * for ( i = 0; i < vector->length; ++i )
 * {
 * LALUniformDeviate( &status, vector->data + i, params );
 * }
 *
 * /\* fill vector with normal deviates *\/
 * LALNormalDeviates( &status, vector, params );
 *
 * LALDestroyRandomParams( &status, &params );
 * LALDestroyVector( &status, &vector );
 * \endcode
 *
 * \heading{Algorithm}
 *
 * This is an implementation of the random number generators \c ran1 and
 * \c gasdev described in Numerical Recipes [\ref ptvf1992].
 *
 */
/*@{*/

static const INT4 a = 16807;
static const INT4 m = 2147483647;
static const INT4 q = 127773;
static const INT4 r = 2836;

static const REAL4 eps = 1.2e-7;

/*
 *
 * XLAL Routines.
 *
 */

INT4 XLALBasicRandom( INT4 i )
{
  INT4 k;
  k = i/q;
  i = a*(i - k*q) - r*k;
  if (i < 0)
    i += m;
  return i;
}

RandomParams * XLALCreateRandomParams( INT4 seed )
{
  RandomParams *params;

  params = XLALMalloc( sizeof( *params) );
  if ( ! params )
    XLAL_ERROR_NULL( XLAL_ENOMEM );

  while ( seed == 0 ) /* use system clock to get seed */
    seed = time( NULL );

  if ( seed < 0 )
    seed = -seed;

  XLALResetRandomParams( params, seed );

  return params;
}


void XLALResetRandomParams( RandomParams *params, INT4 seed )
{
  UINT4 n;

  params->i = seed;
  for ( n = 0; n < 8; ++n ) /* shuffle 8 times */
    params->i = XLALBasicRandom( params->i );

  /* populate vector of random numbers */
  for ( n = 0; n < sizeof( params->v )/sizeof( *params->v ); ++n )
    params->v[n] = params->i = XLALBasicRandom( params->i );

  /* first random number is the 0th element of v */
  params->y = params->v[0];

}


void XLALDestroyRandomParams( RandomParams *params )
{
  XLALFree( params );
}


REAL4 XLALUniformDeviate( RandomParams *params )
{
  REAL4 ans;
  INT4 ndiv;
  INT4 n;

  if ( ! params )
    XLAL_ERROR_REAL4( XLAL_EFAULT );

  /* randomly choose which element of the vector of random numbers to use */
  ndiv = 1 + (m - 1)/(sizeof(params->v)/sizeof(*params->v));
  n = params->y/ndiv;
  params->y = params->v[n];

  /* repopulate this element */
  params->v[n] = params->i = XLALBasicRandom( params->i );

  ans = params->y/(REAL4)m;
  if ( ans > 1 - eps ) /* make sure it is not exactly 1 */
    ans = 1 - eps;

  return ans;
}


int XLALNormalDeviates( REAL4Vector *deviates, RandomParams *params )
{
  REAL4 *data;
  INT4   half;

  if ( ! deviates || ! deviates->data || ! params )
    XLAL_ERROR( XLAL_EFAULT );
  if ( ! deviates->length )
    XLAL_ERROR( XLAL_EBADLEN );

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
      u = XLALUniformDeviate( params );
      v = XLALUniformDeviate( params );
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
      u = XLALUniformDeviate( params );
      v = XLALUniformDeviate( params );
      x   = 2*u - 1;
      y   = 2*v - 1;
      rsq = x*x + y*y;
    }
    while (rsq >= 1 || rsq == 0);

    fac   = sqrt(-2*log(rsq)/rsq);
    *data = fac*x;
    /* throw away y */
  }

  return XLAL_SUCCESS;
}

REAL4 XLALNormalDeviate( RandomParams *params )
{
  REAL4Sequence *deviates;
  REAL4 deviate;

  if ( ! params )
    XLAL_ERROR_REAL4( XLAL_EFAULT );

  /* create a vector */
  deviates = XLALCreateREAL4Sequence(1);
  if(!deviates)
    XLAL_ERROR_REAL4( XLAL_EFUNC );

  /* call the actual function */
  XLALNormalDeviates( deviates, params );
  deviate = deviates->data[0];

  /* destroy the vector */
  XLALDestroyREAL4Sequence(deviates);

  return deviate;
}

/*
 *
 * LAL Routines.
 *
 */



void
LALCreateRandomParams (
    LALStatus     *status,
    RandomParams **params,
    INT4           seed
    )
{
  INITSTATUS(status);

  ASSERT (params, status, RANDOMH_ENULL, RANDOMH_MSGENULL);
  ASSERT (!*params, status, RANDOMH_ENNUL, RANDOMH_MSGENNUL);

  *params = XLALCreateRandomParams( seed );
  if ( ! params )
  {
    XLALClearErrno();
    ABORT( status, RANDOMH_ENULL, RANDOMH_MSGENULL );
  }

  RETURN (status);
}



void
LALDestroyRandomParams (
    LALStatus     *status,
    RandomParams **params
    )
{
  INITSTATUS(status);

  ASSERT (params, status, RANDOMH_ENULL, RANDOMH_MSGENULL);
  ASSERT (*params, status, RANDOMH_ENULL, RANDOMH_MSGENULL);

  XLALDestroyRandomParams( *params );
  *params = NULL;

  RETURN (status);
}



void
LALUniformDeviate (
    LALStatus    *status,
    REAL4        *deviate,
    RandomParams *params
    )
{
  INITSTATUS(status);

  ASSERT (deviate, status, RANDOMH_ENULL, RANDOMH_MSGENULL);
  ASSERT (params, status, RANDOMH_ENULL, RANDOMH_MSGENULL);

  *deviate = XLALUniformDeviate( params );
  if ( XLAL_IS_REAL4_FAIL_NAN( *deviate ) )
  {
    XLALClearErrno();
    ABORT( status, RANDOMH_ENULL, RANDOMH_MSGENULL );
  }

  RETURN (status);
}



void
LALNormalDeviates (
    LALStatus    *status,
    REAL4Vector  *deviates,
    RandomParams *params
    )
{
  INITSTATUS(status);

  ASSERT (params, status, RANDOMH_ENULL, RANDOMH_MSGENULL);
  ASSERT (deviates, status, RANDOMH_ENULL, RANDOMH_MSGENULL);
  ASSERT (deviates->data, status, RANDOMH_ENULL, RANDOMH_MSGENULL);
  ASSERT (deviates->length > 0, status, RANDOMH_ESIZE, RANDOMH_MSGESIZE);

  if ( XLALNormalDeviates( deviates, params ) != XLAL_SUCCESS )
  {
    int errnum = xlalErrno;
    XLALClearErrno();
    switch ( errnum )
    {
      case XLAL_EFAULT:
        ABORT( status, RANDOMH_ENULL, RANDOMH_MSGENULL );
      case XLAL_EBADLEN:
        ABORT( status, RANDOMH_ESIZE, RANDOMH_MSGESIZE );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN (status);
}
/*@}*/
