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

#ifndef _LALDATATYPES_H
#include "LALDatatypes.h"
#ifndef _LALDATATYPES_H
#define _LALDATATYPES_H
#endif
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

typedef void (REAL4LALFunction) (Status *s, REAL4 *y, REAL4 x, void *p);

typedef struct
tagFindRootIn
{
  REAL4LALFunction *function;
  REAL4             xmax;
  REAL4             xmin;
  REAL4             xacc;
}
FindRootIn;

void
BracketRoot (
    Status     *status,
    FindRootIn *inout,
    void       *params
    );

void
BisectionFindRoot (
    Status     *status,
    REAL4      *root,
    FindRootIn *input,
    void       *params
    );

#endif
