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

/* -------------------------------------------------------------------------
 * This function minimises fContact defined above using the
 * Brent method. It returns the minima with a negative sign (which then
 * becomes the maxima of the actual contact function. This can be compared
 * to 1 to check if two ellipsoids indeed overlap.
 * ------------------------------------------------------------------------*/
REAL8 XLALCheckOverlapOfEllipsoids (
        gsl_vector         *ra, 
        gsl_vector         *rb,
        fContactWorkSpace  *workSpace )
{
    gsl_function        F;
    INT4                min_status;
    INT4                iter = 0, max_iter = 100;
    REAL8               m = 0.6180339887;
    REAL8               a = 0.0L, b = 1.0L;
    gsl_min_fminimizer  *s = workSpace->s;

    gsl_vector_memcpy( workSpace->r_AB, rb);
    gsl_vector_sub (workSpace->r_AB, ra);

    F.function = &fContact;
    F.params   = workSpace;

    gsl_min_fminimizer_set (s, &F, m, a, b);

    do
    { 
        iter++;
        min_status = gsl_min_fminimizer_iterate (s);

        m = gsl_min_fminimizer_x_minimum (s);
        a = gsl_min_fminimizer_x_lower (s);
        b = gsl_min_fminimizer_x_upper (s);

        min_status = gsl_min_test_interval (a, b, workSpace->convParam, 0.0);
    }
    while (min_status == GSL_CONTINUE && iter < max_iter && s->f_minimum > -1.0L  );
    /* End of minimization routine */


    return ( -(s->f_minimum) );
}





