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

#ifndef _LALDATATYPES_H
#include "LALDatatypes.h"
#ifndef _LALDATATYPES_H
#define _LALDATATYPES_H
#endif
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

typedef void (REAL4LALFunction) (Status *s, REAL4 *y, REAL4 x, void *p);


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
tagIntegrateIn
{
  REAL4LALFunction *function;
  REAL4             xmax;
  REAL4             xmin;
  IntegralType      type;
}
IntegrateIn;


void
RombergIntegrate (
    Status      *status,
    REAL4       *result,
    IntegrateIn *input,
    void        *params
    );

#endif /* _INTEGRATE_H */
