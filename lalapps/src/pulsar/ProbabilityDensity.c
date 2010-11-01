/*
 * Copyright (C) 2010 Reinhard Prix
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

/*********************************************************************************/
/** \author R. Prix
 * \file
 * \brief
 * Module implementing a probability density function (pdf) object,
 * and useful methods to operate on such objects.
 *
 *
 *********************************************************************************/
#include "config.h"

/* System includes */
#include <math.h>

/* gsl includes */
#include <gsl/gsl_math.h>

/* LAL-includes */
#include <lal/XLALError.h>
#include <lal/AVFactories.h>
#include <lal/LALMalloc.h>
#include <lal/LALStdlib.h>

#include "ProbabilityDensity.h"

/* ----- MACRO definitions ---------- */

/* ---------- internal type definitions ---------- */

/* ---------- internal prototypes ---------- */
/* empty struct initializers */


/* ==================== function definitions ==================== */
/** Function to generate random samples drawn from the given pdf(x)
 *
 * NOTE: if the 'sampling' field is NULL, it will be set the first call to this function.
 */
REAL8
XLALDrawFromPDF1D ( pdf1D_t *pdf,	/**< [in] probability density to sample from */
                    const gsl_rng *rng	/**< random-number generator */
                    )
{
  const char *fn = __func__;

  /* check input consistency */
  if ( !pdf || !rng ) {
    XLALPrintError ("%s: NULL input 'pdf = %p' or 'rng = %p'\n", fn, pdf, rng );
    XLAL_ERROR_REAL8 ( fn, XLAL_EINVAL );
  }
  if ( XLALCheckValidPDF1D ( pdf ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: invalid pdf.\n", fn );
    XLAL_ERROR_REAL8 ( fn, XLAL_EFUNC );
  }

  /* ----- special case 1: single value with certainty */
  if ( pdf->xTics->length == 1 )
    return pdf->xTics->data[0];

  /* ----- special case 2: uniform pdf in [xMin, xMax] */
  if ( pdf->xTics->length == 2 )
    return gsl_ran_flat ( rng, pdf->xTics->data[0], pdf->xTics->data[1] );

  /* ----- general case: draw from discretized pdf ----- */


  // ----- buffer gsl_ran_discrete_t if not computed previously
  if ( pdf->sampling == NULL )
    {
      // first need to prepare discrete *probabilities* prob[i] from probability *density* function pdf[i]
      // namely simply using the respective bin-sizes xBin[i]: prob[i] = pdf[i] * xBin[i]
      UINT4 numBins = pdf->xTics->length - 1;
      REAL8 *prob;
      if ( (prob = XLALMalloc ( numBins * sizeof(*prob) )) == NULL ) {
        XLALPrintError ("%s: failed to XLALMalloc %d-array of REAL8s\n", fn, numBins );
        XLAL_ERROR_REAL8 ( fn, XLAL_ENOMEM );
      }
      UINT4 i;
      for ( i = 0; i < numBins; i ++ )
        {
          REAL8 dx_i = pdf->xTics->data[i+1] - pdf->xTics->data[i];
          prob[i] = pdf->probDens->data[i] * dx_i;
        }

      if ( (pdf->sampling = gsl_ran_discrete_preproc ( numBins, prob ) ) == NULL ) {
        XLALPrintError ("%s: gsl_ran_discrete_preproc() failed\n", fn );
        XLAL_ERROR ( fn, XLAL_EFAILED );
      }

      XLALFree ( prob );
    } /* if pdf->sampling == NULL */

  // ----- draw an index from pdf->probDens
  UINT4 ind = gsl_ran_discrete (rng, pdf->sampling );
  // get the corresponding bin-boundaries of bin[i] = [x[i], x[i+1]]
  REAL8 x0 = pdf->xTics->data[ind];
  REAL8 x1 = pdf->xTics->data[ind+1];

  // and do another uniform draw from [x0, x1]	(thanks to Karl for that suggestion ;)
  // this could possibly be further smoothed by drawing from a linear pdf interpolating the neighbouring bins..
  return gsl_ran_flat ( rng, x0, x1 );

} /* XLALDrawFromPDF1D() */

/** Checks internal consistency of pdf1D object.
 *
 * If lalDebugLevel > 0, also checks normalization of pdf if it claims to be normalized.
 *
 *
 * Return: XLAL_SUCCESS if pdf seems OK, XLAL-error otherwise
 */
int
XLALCheckValidPDF1D ( const pdf1D_t *pdf )
{
  const char *fn = __func__;

  /* check input consistency */
  if ( !pdf ) {
    XLALPrintError ("%s: NULL input 'pdf'\n", fn );
    XLAL_ERROR_REAL8 ( fn, XLAL_EINVAL );
  }

  /* check presence of required xTics array */
  if ( (pdf->xTics == NULL) || (pdf->xTics->length==0) || (pdf->xTics->data==NULL) ) {
    XLALPrintError ("%s: invalid pdf->xTics = NULL, length 0 or NULL data field \n", fn );
    XLAL_ERROR ( fn, XLAL_EDOM );
  }

  /* ----- allowed special case 1: single value with certainty */
  if ( (pdf->xTics->length == 1) && (pdf->probDens == NULL) )
    return XLAL_SUCCESS;

  /* ----- allowed special case 2: uniform pdf in [xMin, xMax] */
  if ( (pdf->xTics->length == 2) && (pdf->probDens == NULL) )
    return XLAL_SUCCESS;

  /* ----- general case: discretized pdf ----- */
  UINT4 numBins = pdf->xTics->length - 1;

  /* check valid pdf array */
  if ( pdf->probDens == NULL ) {
    XLALPrintError ("%s: invalid NULL pdf->probDens for required numBins=%d\n", numBins );
    XLAL_ERROR ( fn, XLAL_EDOM );
  }
  if ( pdf->probDens->length != numBins ) {
    XLALPrintError ("%s: invalid length pdf->probDens->length=%d but should be %d\n", pdf->probDens->length, numBins );
    XLAL_ERROR ( fn, XLAL_EDOM );
  }
  if ( pdf->probDens->data == NULL ) {
    XLALPrintError ("%s: invalid pdf->probDens->data = NULL\n", fn );
    XLAL_ERROR ( fn, XLAL_EDOM );
  }

  /* ----- if pdf claims to be normalized, check that ----- */
  if ( (lalDebugLevel > 0) && pdf->isNormalized )
    {
      REAL8 norm = 0;
      UINT4 i;
      for ( i=0; i < numBins; i ++ )
        {
          REAL8 dx_i = pdf->xTics->data[i+1] - pdf->xTics->data[i];
          norm += pdf->probDens->data[i] * dx_i;	// sum up probability in bin[i]: prob[i] = pdf[i] * dx[i]
        } /* for i < numBins */
      REAL8 relErr = 1e-12;	// generous, given double precision
      if ( gsl_fcmp (norm, 1.0, relErr) != 0 ) {
        XLALPrintError ("%s: pdf claims to be normalized, but norm = %.16f differs from 1.0 by more than %g relative error\n", fn, norm, relErr );
        XLAL_ERROR ( fn, XLAL_EDOM );
      }
    } // if pdf->isNormalized

  /* we found no problem, so we assume it's OK */
  return XLAL_SUCCESS;

} /* XLALCheckValidPDF1D() */

/** Destructor function for 1-D pdf
 */
void
XLALDestroyPDF1D ( pdf1D_t *pdf )
{
  if ( !pdf )
    return;

  if ( pdf->xTics )
    XLALDestroyREAL8Vector ( pdf->xTics );
  if ( pdf->probDens )
    XLALDestroyREAL8Vector ( pdf->probDens );

  if ( pdf->sampling )
    gsl_ran_discrete_free ( pdf->sampling );


  XLALFree ( pdf );

  return;

} /* XLALDestroyPDF1D() */


/** Creator function for a 'singular' 1D pdf, containing a single value with certainty, ie P(x0)=1, and P(x!=x0)=0
 *
 * This is encoded as an xTics array containing just one value: x0, and prob=NULL, sampling=NULL
 */
pdf1D_t *
XLALCreateSingularPDF1D ( REAL8 x0	/**< domain of pdf is a single point: x0 */
                          )
{
  const char *fn = __func__;

  /* allocate memory for output pdf */
  pdf1D_t *ret;

  if ( ( ret = XLALCalloc ( 1, sizeof(*ret) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc ( 1, %d)\n", fn, sizeof(*ret) );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  if ( ( ret->xTics = XLALCreateREAL8Vector ( 1 )) == NULL ) {
    XLALPrintError ("%s: surprisingly, XLALCreateREAL8Vector(1) failed!\n", fn );
    XLALFree ( ret );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  ret->xTics->data[0] = x0;	/* only value required: P[x0]=1 */

  return ret;

} /* XLALCreateSingularPDF1D() */

/** Creator function for a uniform 1D pdf over [xMin, xMax]
 *
 * This is encoded as an xTics array containing just two values: x[0]=xMin, x[1]=xMax,
 * and prob=NULL, sampling=NULL {not required to draw from this pdf}
 */
pdf1D_t *
XLALCreateUniformPDF1D ( REAL8 xMin,	/**< lower boundary of domain interval */
                         REAL8 xMax	/**< upper boundary of domain interval */
                         )
{
  const char *fn = __func__;

  /* check input */
  if ( xMax < xMin ) {
    XLALPrintError ("%s: invalid input, xMax=%f must be > xMin = %f\n", fn, xMax, xMin );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  /* allocate memory for output pdf */
  pdf1D_t *ret;

  if ( ( ret = XLALCalloc ( 1, sizeof(*ret) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc ( 1, %d)\n", fn, sizeof(*ret) );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  if ( ( ret->xTics = XLALCreateREAL8Vector ( 2 )) == NULL ) {
    XLALPrintError ("%s: surprisingly, XLALCreateREAL8Vector(2) failed!\n", fn );
    XLALFree ( ret );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  ret->xTics->data[0] = xMin;
  ret->xTics->data[1] = xMax;

  return ret;

} /* XLALCreateUniformPDF1D() */


/** Creator function for a generic discrete 1D pdf over [xMin, xMax], discretized into numBins bins
 *
 * NOTE: generates a uniform sampling of the domain [xMin, xMax] in numBins
 * NOTE2: returns the P[i] array 'prob' initialized to 0, so after calling this function
 * the user still needs to feed in the correct values for the probabilities P[i] of x in [x[i],x[i+1]]
 */
pdf1D_t *
XLALCreateDiscretePDF1D ( REAL8 xMin,	/**< lower boundary of domain interval */
                          REAL8 xMax,	/**< upper boundary of domain interval */
                          UINT4 numBins /**< number of bins to discretize PDF into */
                         )
{
  const char *fn = __func__;

  /* check input */
  if ( xMax < xMin ) {
    XLALPrintError ("%s: invalid input, xMax=%f must be > xMin = %f\n", fn, xMax, xMin );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }
  if ( numBins == 0 ) {
    XLALPrintError ("%s: invalid input, numBins must be positive!\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  /* allocate memory for output pdf */
  pdf1D_t *ret;

  if ( ( ret = XLALCalloc ( 1, sizeof(*ret) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc ( 1, %d)\n", fn, sizeof(*ret) );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  if ( ( ret->xTics = XLALCreateREAL8Vector ( numBins + 1 )) == NULL ) {
    XLALPrintError ("%s: surprisingly, XLALCreateREAL8Vector(%d) failed!\n", fn, numBins + 1 );
    XLALFree ( ret );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }
  if ( ( ret->probDens = XLALCreateREAL8Vector ( numBins )) == NULL ) {
    XLALPrintError ("%s: surprisingly, XLALCreateREAL8Vector(%d) failed!\n", fn, numBins );
    XLALDestroyREAL8Vector ( ret->xTics );
    XLALFree ( ret );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  /* initialize N bins uniformly spaced over [xMin, xMax], ie. N+1 'tics' */
  UINT4 i;
  REAL8 dx = ( xMax - xMin ) / numBins;
  for (i = 0; i < numBins + 1; i ++ )
    {
      REAL8 xi = xMin + i * dx;
      ret->xTics->data[i] = xi;

    } /* for i < numBins+1 */

  /* initialized pdf bins to zero */
  memset ( ret->probDens->data, 0, numBins * sizeof( *ret->probDens->data ) );

  return ret;

} /* XLALCreateDiscretePDF1D() */

/** Method to normalize the given pdf1D.
 * Only does something if necessary, ie if pdf isn't normalized already
 */
int
XLALNormalizePDF1D ( pdf1D_t *pdf )
{
  const char *fn = __func__;

  /* sanity checks */
  if ( !pdf  ) {
    XLALPrintError ("%s: invalid NULL input 'pdf'\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( XLALCheckValidPDF1D ( pdf ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: invalid pdf\n", fn );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }

  /* ----- special case 1: single value with certainty: nothing to be done */
  if ( (pdf->xTics->length == 1) && (pdf->probDens == NULL) )
    return XLAL_SUCCESS;

  /* ----- allowed special case 2: uniform pdf in [xMin, xMax] */
  if ( (pdf->xTics->length == 2) && (pdf->probDens == NULL) )
    return XLAL_SUCCESS;

  /* is normalized already? nothing to be done */
  if ( pdf->isNormalized )
    return XLAL_SUCCESS;

  /* ----- general case: compute norm = int pdf(x) dx and divide pdf(x) to properly normalize */
  UINT4 numBins = pdf->xTics->length - 1;
  UINT4 i;
  REAL8 norm = 0;
  for ( i=0; i < numBins; i ++ )
    {
      REAL8 dx_i = pdf->xTics->data[i+1] - pdf->xTics->data[i];
      norm += pdf->probDens->data[i] * dx_i;	// sum up probability in bin[i]: prob[i] = pdf[i] * dx[i]
    } /* for i < numBins */

  REAL8 invNorm = 1.0 / norm;
  for ( i=0; i < numBins; i ++ )
    {
      pdf->probDens->data[i] *= invNorm;
    }

  pdf->isNormalized = 1;

  return XLAL_SUCCESS;

} /* XLALNormalizePDF1D() */

