/*----------------------------------------------------------------------- 
 * 
 * File Name: Integrate.h
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _INTEGRATE_H
#define _INTEGRATE_H

#include "LALDatatypes.h"

#ifdef __cplusplus
extern "C" {
#endif

NRCSID (INTEGRATEH, "$Id$");


#define INTEGRATE_ENULL 1
#define INTEGRATE_ETYPE 2
#define INTEGRATE_EIDOM 4
#define INTEGRATE_EMXIT 8

#define INTEGRATE_MSGENULL "Null pointer"
#define INTEGRATE_MSGETYPE "Unknown integral type"
#define INTEGRATE_MSGEIDOM "Invalid domain"
#define INTEGRATE_MSGEMXIT "Maximum iterations exceeded"

typedef void (REAL4LALFunction) (LALStatus *s, REAL4 *y, REAL4 x, void *p);
typedef void (REAL8LALFunction) (LALStatus *s, REAL8 *y, REAL8 x, void *p);


typedef enum
{
  ClosedInterval,     /* evaluate integral on a closed interval             */
  OpenInterval,       /* evaluate integral on an open interval              */
  SingularLowerLimit, /* integrate an inv sqrt singularity at lower limit   */
  SingularUpperLimit, /* integrate an inv sqrt singularity at upper limit   */
  InfiniteDomainPow,  /* integrate infinite domain with power-law falloff   */
  InfiniteDomainExp   /* integrate infinite domain with exponential falloff */
}
IntegralType;


typedef struct
tagSIntegrateIn
{
  REAL4LALFunction *function;
  REAL4             xmax;
  REAL4             xmin;
  IntegralType      type;
}
SIntegrateIn;


typedef struct
tagDIntegrateIn
{
  REAL8LALFunction *function;
  REAL8             xmax;
  REAL8             xmin;
  IntegralType      type;
}
DIntegrateIn;


void
LALSRombergIntegrate (
    LALStatus       *status,
    REAL4        *result,
    SIntegrateIn *input,
    void         *params
    );


void
LALDRombergIntegrate (
    LALStatus       *status,
    REAL8        *result,
    DIntegrateIn *input,
    void         *params
    );


#ifdef __cplusplus
}
#endif

#endif /* _INTEGRATE_H */
