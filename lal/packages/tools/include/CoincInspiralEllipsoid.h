/*-----------------------------------------------------------------------
 *
 * File Name: CoincInspiralEllipsoid.h
 *
 * Author: Robinson, C. A.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="CoincInspiralEllipsoidHV">
Author: Robinson, C. A.
$Id$
</lalVerbatim>
#endif

#ifndef _COINCINSPIRALELLIPSOID_H
#define _COINCINSPIRALELLIPSOID_H


#include    <math.h>
#include    <lal/LALStdlib.h>
#include    <lal/LALGSL.h>
#include    <lal/LALError.h>
#include    <lal/LIGOMetadataTables.h>
#include    <lal/LIGOMetadataUtils.h>
#include    <lal/EllipsoidOverlapTools.h>

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

NRCSID( COINCINSPIRALELLIPSOIDH, "$Id$" );


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

gsl_vector * XLALGetPositionFromSnglInspiral( SnglInspiralTable *table );

int XLALSetTimeInPositionVector( gsl_vector *position,
                                 REAL8       time
                               );

/* Function for calculating the e-thinca parameter between coincident events */
REAL8 XLALCalculateEThincaParameter(
          SnglInspiralTable *table1,
          SnglInspiralTable *table2
             );

#ifdef  __cplusplus
}
#endif

#endif   /* _COINCINSPIRALELLIPSOID_H */
