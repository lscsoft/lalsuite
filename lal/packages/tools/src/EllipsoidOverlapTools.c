/*
*  Copyright (C) 2007 Anand Sengupta, Craig Robinson
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

/*-----------------------------------------------------------------------
 *
 * File Name: EllipsoidOverlapTools.c
 *
 * Author: Anand S. Sengupta and Craig Robinson
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="EllipsoidOverlapToolsCV">
Author: Anand S. Sengupta and Craig Robinson
$Id$
</lalVerbatim>
#endif

#include <lal/EllipsoidOverlapTools.h>

NRCSID( ELLIPSOIDOVERLAPTOOLSC,
        "$Id$" );


static REAL8 fContact (REAL8 x, void *params);

/* ---------------------------------------------------------------------------
 * This function return the contact function of two ellipsoids as defined in
 * the paper by Perram & Wertheim, JCP, v58, pp 409-416, eq 3.7 with a
 * negative sign. This is minimised to figure out if the ellipsoids
 * intersect.
 * --------------------------------------------------------------------------*/



static REAL8 fContact (REAL8 x, void *params)
{
    fContactWorkSpace *p
            = (fContactWorkSpace *)params;

    gsl_permutation   *p1    = p->p1;
    gsl_vector        *tmpV  = p->tmpV;
    gsl_matrix        *C     = p->C;
    gsl_matrix        *A     = p->tmpA;
    gsl_matrix        *B     = p->tmpB;

    INT4     s1;
    REAL8    result;

    gsl_matrix_memcpy ( A, p->invQ1);
    gsl_matrix_memcpy ( B, p->invQ2);

    gsl_matrix_scale (B, x);
    gsl_matrix_scale (A, (1.0L-x));

    gsl_matrix_add (A, B);

    gsl_linalg_LU_decomp( A, p1, &s1 );
    gsl_linalg_LU_invert( A, p1, C );

    /* Evaluate the product C x r_AB */
    gsl_blas_dsymv (CblasUpper, 1.0, C, p->r_AB, 0.0, tmpV);

    /* Now evaluate transpose(r_AB) x (C x r_AB) */
    gsl_blas_ddot (p->r_AB, tmpV, &result);

    result *= (x*(1.0L-x));

    return  (-result);
}

/*
-------------------------------------------------------------------------
 * This function allocates and initialises the memory and parameters for
 * the workspace used for determining whether the two ellipsoids intersect.
 * If the shape matrices a and b are not null, these will be used in the
 * workspace; otherwise they will point to null, and will have to be
 * set by hand.
 * ------------------------------------------------------------------------*/
fContactWorkSpace * XLALInitFContactWorkSpace(
                       UINT4                          n,
                       const gsl_matrix              *a,
                       const gsl_matrix              *b,
                       const gsl_min_fminimizer_type *T,
                       REAL8                          conv
                                             )
{

  static const char *func = "XLALInitFContactWorkSpace";
  fContactWorkSpace *workSpace = NULL;

  /* Check the parameters passed in are sensible */
  if ( n <= 0 )
    XLAL_ERROR_NULL( func, XLAL_EINVAL );

  if ( conv <= 0.0 )
    XLAL_ERROR_NULL( func, XLAL_EINVAL );

  if ( !T )
    XLAL_ERROR_NULL( func, XLAL_EFAULT );

  /* The matrices a and b should be both supplied or both null */
  if ((!( a && b )) && ( a || b ))
        XLAL_ERROR_NULL( func, XLAL_EINVAL );

  /* Check the matrices conform to the expected sizes */
  if ( a )
  {
    if ( a->size1 != n || a->size2 != n || b->size1 != n || b->size2 != n )
    {
      XLAL_ERROR_NULL( func, XLAL_EBADLEN );
    }
  }


  /* Allocate the workspace */
  workSpace = (fContactWorkSpace *) LALCalloc( 1, sizeof(fContactWorkSpace) );
  if ( !workSpace )
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );

  if ( a )
  {
    workSpace->invQ1         = a;
    workSpace->invQ2         = b;
  }

  XLAL_CALLGSL( workSpace->tmpA = gsl_matrix_calloc( n, n ) );
  XLAL_CALLGSL( workSpace->tmpB = gsl_matrix_calloc( n, n ) );
  XLAL_CALLGSL( workSpace->C    = gsl_matrix_calloc( n, n ) );
  XLAL_CALLGSL( workSpace->p1   = gsl_permutation_alloc( n ) );
  XLAL_CALLGSL( workSpace->tmpV = gsl_vector_calloc( n ) );
  XLAL_CALLGSL( workSpace->r_AB = gsl_vector_calloc( n ) );
  XLAL_CALLGSL( workSpace->s    = gsl_min_fminimizer_alloc( T ) );

  /* Check all of the above were allocated properly */
  if (!workSpace->tmpA || !workSpace->tmpB || !workSpace->C ||
      !workSpace->p1 || !workSpace->tmpV || !workSpace->r_AB
      || !workSpace->s )
  {
    XLALFreeFContactWorkSpace( workSpace );
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );
  }

  /* Now set the rest of the parameters to the correct values */
  workSpace->T         = T;
  workSpace->convParam = conv;
  workSpace->n         = n;

  return workSpace;
}

/*
-------------------------------------------------------------------------
 * This function minimises fContact defined above using the
 * Brent method. It returns the minima with a negative sign (which then
 * becomes the maxima of the actual contact function. This can be compared
 * to 1 to check if two ellipsoids indeed overlap.
 * ------------------------------------------------------------------------*/
REAL8 XLALCheckOverlapOfEllipsoids (
        const gsl_vector   *ra,
        const gsl_vector   *rb,
        fContactWorkSpace  *workSpace )
{
    static const char *func = "XLALCheckOverlapOfEllipsoids";

    gsl_function        F;
    INT4                min_status;
    INT4                iter = 0, max_iter = 100;
    REAL8               m = 0.6180339887;
    REAL8               a = 0.0L, b = 1.0L;
    gsl_min_fminimizer  *s = workSpace->s;

    /* Sanity check on input arguments */
    if ( !ra || !rb || !workSpace )
      XLAL_ERROR_REAL8( func, XLAL_EFAULT );

    if ( ra->size != rb->size || ra->size != workSpace->n )
      XLAL_ERROR_REAL8( func, XLAL_EBADLEN);


    /* Set r_AB to be rb - ra */
    XLAL_CALLGSL( gsl_vector_memcpy( workSpace->r_AB, rb) );
    XLAL_CALLGSL( gsl_vector_sub (workSpace->r_AB, ra) );

    if ( gsl_vector_isnull( workSpace->r_AB ))
    {
      XLALPrintWarning("Position vectors ra and rb are identical.\n");
      return 0;
    }

    F.function = &fContact;
    F.params   = workSpace;

    XLAL_CALLGSL( min_status = gsl_min_fminimizer_set (s, &F, m, a, b) );
    if ( min_status != GSL_SUCCESS )
      XLAL_ERROR_REAL8( func, XLAL_EFUNC );

    do
    {
        iter++;
        XLAL_CALLGSL( min_status = gsl_min_fminimizer_iterate (s) );
        if (min_status != GSL_SUCCESS )
        {
          if (min_status == GSL_EBADFUNC)
            XLAL_ERROR_REAL8( func, XLAL_EFUNC | XLAL_EFPINVAL );
          else if (min_status == GSL_EZERODIV)
            XLAL_ERROR_REAL8( func, XLAL_EFUNC | XLAL_EFPDIV0 );
          else
            XLAL_ERROR_REAL8( func, XLAL_EFUNC );
        }

        m = gsl_min_fminimizer_x_minimum (s);
        a = gsl_min_fminimizer_x_lower (s);
        b = gsl_min_fminimizer_x_upper (s);

        XLAL_CALLGSL( min_status = gsl_min_test_interval (a, b, workSpace->convParam, 0.0) );
        if (min_status != GSL_CONTINUE && min_status != GSL_SUCCESS )
          XLAL_ERROR_REAL8( func, XLAL_EFUNC );
    }
    while (min_status == GSL_CONTINUE && iter < max_iter );
    /* End of minimization routine */

    /* Throw an error if max iterations would have been exceeded */
    if ( iter == max_iter && min_status == GSL_CONTINUE )
    {
      XLAL_ERROR_REAL8( func, XLAL_EMAXITER );
    }

    return ( -(s->f_minimum) );
}


/*
-------------------------------------------------------------------------
 * This function frees the memory allocated using XLALInitFContactWorkSpace.
 * ------------------------------------------------------------------------*/
void XLALFreeFContactWorkSpace( fContactWorkSpace *workSpace )
{

  static const char *func = "XLALFreeFContactWorkSpace";

  if ( !workSpace )
    XLAL_ERROR_VOID( func, XLAL_EFAULT );

  /* Free all the allocated memory */
  if (workSpace->tmpA) XLAL_CALLGSL( gsl_matrix_free( workSpace->tmpA ));
  if (workSpace->tmpB) XLAL_CALLGSL( gsl_matrix_free( workSpace->tmpB ));
  if (workSpace->C)    XLAL_CALLGSL( gsl_matrix_free( workSpace->C ));
  if (workSpace->p1)   XLAL_CALLGSL( gsl_permutation_free( workSpace->p1 ));
  if (workSpace->tmpV) XLAL_CALLGSL( gsl_vector_free( workSpace->tmpV ));
  if (workSpace->r_AB) XLAL_CALLGSL( gsl_vector_free( workSpace->r_AB ));
  if (workSpace->s)    XLAL_CALLGSL( gsl_min_fminimizer_free( workSpace->s ));

  LALFree( workSpace );
  return;
}
