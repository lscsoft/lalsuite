#include    <math.h>
#include    <lal/LALStdlib.h>
#include    <lal/LALGSL.h>
#include    <lal/LALError.h>
#include    <lal/LIGOMetadataUtils.h>
#include    <gsl/gsl_errno.h>
#include    <gsl/gsl_math.h>
#include    <gsl/gsl_min.h>
#include    <gsl/gsl_vector.h>
#include    <gsl/gsl_matrix.h>
#include    <gsl/gsl_blas.h>
#include    <gsl/gsl_linalg.h>

#ifndef _ELLIPSOIDOVERLAPTOOLSH
#define _ELLIPSOIDOVERLAPTOOLSH

typedef struct tagfContactWorkSpace
{
    /* Dimension of the matrices & vectors */
    INT4              n;

    /* Vector r_AB = r_B - r_A and           */
    /* shape matrices for ellipsoid centered */
    /* at points A and B                     */
    gsl_vector        *r_AB;
    const gsl_matrix  *invQ1, *invQ2;

    /* Parameters for minimizing the contact function */
    REAL8             convParam;
    const             gsl_min_fminimizer_type *T;
    gsl_min_fminimizer  *s;

    /* Temporary workspace variables */
    gsl_matrix        *tmpA, *tmpB, *C;
    gsl_vector        *tmpV;
    gsl_permutation   *p1;
}
fContactWorkSpace;


/* Function Prototypes */
REAL8 XLALCheckOverlapOfEllipsoids (
        const gsl_vector   *ra, 
        const gsl_vector   *rb,
        fContactWorkSpace  *workSpace );


fContactWorkSpace * XLALInitFContactWorkSpace(
                       INT4                           n,
                       const gsl_matrix              *a,
                       const gsl_matrix              *b,
                       const gsl_min_fminimizer_type *T,
                       REAL8                          conv
                                             );

void XLALFreeFContactWorkSpace( fContactWorkSpace *workSpace );

/* Functions for generating the error matrix and position vectors for triggers */
gsl_matrix * XLALGetErrorMatrixFromSnglInspiral(SnglInspiralTable *event);

gsl_vector * XLALGetPositionFromSnglInspiral( SnglInspiralTable *table );

#endif
