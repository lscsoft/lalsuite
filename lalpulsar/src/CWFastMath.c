// Copyright (C) 2014 Karl Wette
// Copyright (C) 2011, 2014 David Keitel
// Copyright (C) 2011 Reinhard Prix
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
//

/*---------- INCLUDES ----------*/
#include "CWFastMath.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_vector.h>

/*---------- local DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)

/*----- Macros ----- */

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/

/* empty initializers  */

/*---------- internal prototypes ----------*/

/* ----- module-local fast lookup-table handling of logarithms ----- */

/**
 * Lookup-table for logarithms log(x)
 * Holds an array 'data' of 'length' for values log(x) for x in the range (0, xmax]
 */
#define LOGLUT_XMAX 	3.0	// LUT for range (0,numDetectors+1), currently numDetectors = 2 FIXME: get this dynamically
#define LOGLUT_LENGTH 	2000	// number of LUT values to pre-compute
static gsl_vector *logLUT = NULL; 	/**< module-global lookup-table for logarithms log(x) */
#define LOGLUT_DXINV  ((LOGLUT_LENGTH)/(LOGLUT_XMAX))	// 1/dx with dx = xmax/length

static int XLALCreateLogLUT ( void );	/* only ever used internally, destructor is in exported API */

/* ----- module-local fast lookup-table handling of negative exponentials ----- */

/**
 * Lookup-table for negative exponentials e^(-x)
 * Holds an array 'data' of 'length' for values e^(-x) for x in the range [0, xmax]
 */
#define EXPLUT_XMAX 	20.0	// LUT down to e^(-20) = 2.0612e-09
#define EXPLUT_LENGTH 	2000	// number of LUT values to pre-compute
static gsl_vector *expLUT = NULL; 	/**< module-global lookup-table for negative exponentials e^(-x) */
#define EXPLUT_DXINV  ((EXPLUT_LENGTH)/(EXPLUT_XMAX))	// 1/dx with dx = xmax/length

static int XLALCreateExpLUT ( void );	/* only ever used internally, destructor is in exported API */

/*==================== FUNCTION DEFINITIONS ====================*/

/** Generate a lookup-table logLUT for log(x) over the interval x in (0, xmax], using 'length' points. */
int
XLALCreateLogLUT ( void )
{
  /* create empty output LUT */
  gsl_vector *ret;
  XLAL_CHECK ( ( ret = gsl_vector_alloc ( LOGLUT_LENGTH + 1)) != NULL, XLAL_ENOMEM, "Failed call to gsl_vector_alloc (%s).", LOGLUT_LENGTH +1 );

  /* fill output LUT */
  REAL8 dx = LOGLUT_XMAX / LOGLUT_LENGTH;
  UINT4 i;
  for ( i=0; i <= LOGLUT_LENGTH; i ++ ) {
    REAL8 xi = i * dx;
    gsl_vector_set ( ret, i, log( xi ) );
  } /* for i < length() */

  /* 'return' this by setting the global vector */
  logLUT = ret;

  return XLAL_SUCCESS;

} /* XLALCreateLogLUT() */


/**
 * Destructor function for logLUT_t lookup table
 */
void
XLALDestroyLogLUT ( void )
{

  if ( !logLUT ) {
    return;
  }

  gsl_vector_free ( logLUT );

  logLUT = NULL;

  return;

} /* XLALDestroyLogLUT() */


/**
 * Fast logarithmic function log(x) using lookup-table (LUT).
 * We need to compute log(x) for x in (0,xmax], typically in a B-stat
 * integral of the form int e^-x dx: this means that small values e^(-x)
 * will not contribute much to the integral and are less important than
 * values close to 1. Therefore we pre-compute a LUT of e^(-x) for x in [0, xmax],
 * in Npoints points, and set e^(-x) = 0 for x < xmax.
 *
 * NOTE: if module-global logLUT=NULL, we create it here
 * NOTE: if argument is outside of (0,xmax], we use math-lib log(x) instead of LUT
 */
REAL8
XLALFastLog ( REAL8 x )
{

  if ( x > LOGLUT_XMAX ) { /* for values bigger than xmax, use normal log function */
    return log(x);
  }

  XLAL_CHECK_REAL8 ( x >= 0, XLAL_EDOM, "Negative argument: x=%f.", x );

  /* if lookup table doesn't exist yet: generate it now */
  XLAL_CHECK_REAL8 ( logLUT || ( XLALCreateLogLUT() == XLAL_SUCCESS), XLAL_EFUNC, "Failed call to XLALCreateLogLUT()." );

  /* find index of closest point xp in LUT to xm */
  UINT4 i0 = (UINT4) ( x * LOGLUT_DXINV + 0.5 );

  return gsl_vector_get ( logLUT, i0 );

} /* XLALFastLog() */


/**
 * Generate an exponential lookup-table expLUT for e^(-x)
 * over the interval x in [0, xmax], using 'length' points.
 */
int
XLALCreateExpLUT ( void )
{
  /* create empty output LUT */
  gsl_vector *ret;
  if ( ( ret = gsl_vector_alloc ( EXPLUT_LENGTH + 1)) == NULL ) {
    XLALPrintError ("%s: failed to gsl_vector_alloc (%s)\n", __func__, EXPLUT_LENGTH +1 );
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  /* fill output LUT */
  REAL8 dx = EXPLUT_XMAX / EXPLUT_LENGTH;
  UINT4 i;
  for ( i=0; i <= EXPLUT_LENGTH; i ++ )
    {
      REAL8 xi = i * dx;

      gsl_vector_set ( ret, i, exp( - xi ) );

    } /* for i < length() */

  /* 'return' this by setting the global vector */
  expLUT = ret;

  return XLAL_SUCCESS;

} /* XLALCreateExpLUT() */


/**
 * Destructor function for expLUT_t lookup table
 */
void
XLALDestroyExpLUT ( void )
{
  if ( !expLUT )
    return;

  gsl_vector_free ( expLUT );

  expLUT = NULL;

  return;

} /* XLALDestroyExpLUT() */


/**
 * Fast exponential function e^-x using lookup-table (LUT).
 * We need to compute exp(-x) for x >= 0, typically in a B-stat
 * integral of the form int e^-x dx: this means that small values e^(-x)
 * will not contribute much to the integral and are less important than
 * values close to 1. Therefore we pre-compute a LUT of e^(-x) for x in [0, xmax],
 * in Npoints points, and set e^(-x) = 0 for x < xmax.
 *
 * NOTE: if module-global expLUT=NULL, we create it here
 * NOTE: if argument is negative, we use math-lib exp(-x) instead of LUT
 */
REAL8
XLALFastNegExp ( REAL8 mx )
{
  if ( mx > EXPLUT_XMAX )	/* for values smaller than e^(-xmax) we truncate to 0 */
    return 0.0;

  if ( mx < 0 )
    return exp ( - mx  );

  /* if lookup table doesn't exist yet: generate it now */
  if ( !expLUT && ( XLALCreateExpLUT() != XLAL_SUCCESS) ) {
    XLAL_ERROR_REAL8 ( XLAL_EFUNC );
  }

  /* find index of closest point xp in LUT to xm */
  UINT4 i0 = (UINT4) ( mx * EXPLUT_DXINV + 0.5 );

  return gsl_vector_get ( expLUT, i0 );

} /* XLALFastNegExp() */
