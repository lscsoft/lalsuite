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
 * File Name: EllipsoidOverlapTools.h
 *
 * Author: Robinson, C. A., and Sengupta, A.
 *
 *
 *-----------------------------------------------------------------------
 */

#if 0
Author: Robinson, C. A. and Sengupta, A. S.
#endif

#ifndef _ELLIPSOIDOVERLAPTOOLS_H
#define _ELLIPSOIDOVERLAPTOOLS_H


#include    <math.h>
#include    <lal/LALStdlib.h>
#include    <lal/LALGSL.h>
#include    <lal/LALError.h>

#include    <gsl/gsl_errno.h>
#include    <gsl/gsl_math.h>
#include    <gsl/gsl_min.h>
#include    <gsl/gsl_vector.h>
#include    <gsl/gsl_matrix.h>
#include    <gsl/gsl_blas.h>
#include    <gsl/gsl_linalg.h>

#ifdef  __cplusplus
extern "C" {
#endif

#ifdef SWIG /* SWIG interface directives */
SWIGLAL(IMMUTABLE_MEMBERS(tagfContactWorkSpace, invQ1, invQ2));
#endif /* SWIG */
typedef struct tagfContactWorkSpace
{
    /* Dimension of the matrices & vectors */
    UINT4             n;

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
                       UINT4                          n,
                       const gsl_matrix              *a,
                       const gsl_matrix              *b,
                       const gsl_min_fminimizer_type *T,
                       REAL8                          conv
                                             );

void XLALFreeFContactWorkSpace( fContactWorkSpace *workSpace );

#ifdef  __cplusplus
}
#endif

#endif   /* _ELLIPSOIDOVERLAPTOOLS_H */
