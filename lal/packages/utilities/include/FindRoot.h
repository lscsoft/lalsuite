/*----------------------------------------------------------------------- 
 * 
 * File Name: FindRoot.h
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _FINDROOT_H
#define _FINDROOT_H

#include "LALDatatypes.h"

#ifdef __cplusplus
extern "C" {
#endif


NRCSID (FINDROOTH, "$Id$");

#define FINDROOT_ENULL 1
#define FINDROOT_EIDOM 2
#define FINDROOT_EMXIT 4
#define FINDROOT_EBRKT 8

#define FINDROOT_MSGENULL "Null pointer"
#define FINDROOT_MSGEIDOM "Invalid initial domain"
#define FINDROOT_MSGEMXIT "Maximum iterations exceeded"
#define FINDROOT_MSGEBRKT "Root not bracketed"

typedef void (REAL4LALFunction) (LALStatus *s, REAL4 *y, REAL4 x, void *p);
typedef void (REAL8LALFunction) (LALStatus *s, REAL8 *y, REAL8 x, void *p);

typedef struct
tagSFindRootIn
{
  REAL4LALFunction *function;
  REAL4             xmax;
  REAL4             xmin;
  REAL4             xacc;
}
SFindRootIn;

typedef struct
tagDFindRootIn
{
  REAL8LALFunction *function;
  REAL8             xmax;
  REAL8             xmin;
  REAL8             xacc;
}
DFindRootIn;

void
LALSBracketRoot (
    LALStatus      *status,
    SFindRootIn *inout,
    void        *params
    );

void
LALDBracketRoot (
    LALStatus      *status,
    DFindRootIn *inout,
    void        *params
    );

void
LALSBisectionFindRoot (
    LALStatus      *status,
    REAL4       *root,
    SFindRootIn *input,
    void        *params
    );

void
LALDBisectionFindRoot (
    LALStatus      *status,
    REAL8       *root,
    DFindRootIn *input,
    void        *params
    );


#ifdef __cplusplus
}
#endif

#endif /* _FINDROOT_H */
