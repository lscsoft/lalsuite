/*
*  Copyright (C) 2007 Alexander Dietz, Drew Keppel, Craig Robinson
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
#ifndef _COINCINSPIRALELLIPSOID_H
#define _COINCINSPIRALELLIPSOID_H

#ifdef  __cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------
 *
 * File Name: CoincInspiralEllipsoid.h
 *
 * Author: Robinson, C. A.
 *
 *-----------------------------------------------------------------------
 */

/**
 * \defgroup CoincInspiralEllipsoid_h Header CoincInspiralEllipsoid.h
 * \ingroup pkg_CBC_NEW
 * \author Robinson, C. A.
 *
 * \brief Provides function definitions for performing inspiral coincidence analysis
 * using error ellipsoids.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/CoincInspiralEllipsoid.h>
 * \endcode
 *
 */
/*@{*/

#include    <math.h>
#include    <lal/LALStdlib.h>
#include    <lal/LALGSL.h>
#include    <lal/LALError.h>
#include    <lal/LIGOMetadataTables.h>
#include    <lal/LIGOMetadataInspiralUtils.h>
#include    <lal/EllipsoidOverlapTools.h>

#include    <gsl/gsl_errno.h>
#include    <gsl/gsl_math.h>
#include    <gsl/gsl_min.h>
#include    <gsl/gsl_vector.h>
#include    <gsl/gsl_matrix.h>
#include    <gsl/gsl_blas.h>
#include    <gsl/gsl_linalg.h>

/**
 * The \c TriggerErrorList is a linked list used within e-thinca. It
 * contains pointers to the \c SnglInspiralTable for a given trigger,
 * and its associated error matrix and position vector.
 */
typedef struct tagTriggerErrorList
{
  SnglInspiralTable          *trigger;
  gsl_matrix                 *err_matrix;
  gsl_vector                 *position;
  struct tagTriggerErrorList *next;
}
TriggerErrorList;


/* Functions for performing coincident analysis */
void
LALCreateTwoIFOCoincListEllipsoid(
    LALStatus                  *status,
    CoincInspiralTable        **coincOutput,
    SnglInspiralTable          *snglInput,
    InspiralAccuracyList       *accuracyParams
    );

/* Functions for checking for coincidence between inspiral events */
INT2 XLALCompareInspiralsEllipsoid(
      TriggerErrorList              *aPtr,
      TriggerErrorList              *bPtr,
      fContactWorkSpace             *workSpace,
      InspiralAccuracyList          *params
      );

void
XLALSnglInspiralCoincTestEllipsoid(
    CoincInspiralTable         *coincInspiral,
    SnglInspiralTable          *snglInspiral,
    InspiralAccuracyList       *accuracyParams
    );

/* Functions for generating the error matrix and position vectors for triggers */
gsl_matrix * XLALGetErrorMatrixFromSnglInspiral(
     SnglInspiralTable *event,
     REAL8              eMatch
     );

int XLALSetErrorMatrixFromSnglInspiral(gsl_matrix        *shape,
                                       SnglInspiralTable *event,
                                       REAL8              eMatch
                                       );

gsl_vector * XLALGetPositionFromSnglInspiral( SnglInspiralTable *table );

int XLALSetTimeInPositionVector( gsl_vector *position,
                                 REAL8       timeShift
                               );

/* Function for calculating the e-thinca parameter between coincident events */
REAL8 XLALCalculateEThincaParameter(
          SnglInspiralTable *table1,
          SnglInspiralTable *table2,
          InspiralAccuracyList* accuracyParams
);

REAL8 XLALEThincaParameterForInjection(
                    SimInspiralTable  *injection,
                    SnglInspiralTable *trigger
                    );

/*@}*/ /* end:CoincInspiralEllipsoid_h */

#ifdef  __cplusplus
}
#endif

#endif   /* _COINCINSPIRALELLIPSOID_H */
