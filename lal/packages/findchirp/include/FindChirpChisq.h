/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpChisq.h
 *
 * Author: Anderson, W. G., and Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _FINDCHIRPCHISQ_H
#define _FINDCHIRPCHISQ_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (SPFINDCHIRPCHISQH, "$Id$");

#define FINDCHIRPCHISQ_ENULL 1
#define FINDCHIRPCHISQ_ENNUL 2
#define FINDCHIRPCHISQ_ECHIZ 3

#define FINDCHIRPCHISQ_MSGENULL "Null pointer"
#define FINDCHIRPCHISQ_MSGENNUL "Non-null pointer"
#define FINDCHIRPCHISQ_MSGECHIZ "Number of chisq bins is zero or negative"


typedef struct
tagFindChirpChisqInput
{
  COMPLEX8Vector               *qtildeVec;
  COMPLEX8Vector               *qVec;
}
FindChirpChisqInput;


typedef struct
tagFindChirpChisqParams
{
  UINT4Vector                  *chisqBinVec;
  ComplexFFTPlan               *plan;
  REAL4                         chisqNorm;
}
FindChirpChisqParams;


void
LALFindChirpChisqVeto (
    LALStatus                  *status,
    REAL4Vector                *chisqVec,
    FindChirpChisqInput        *input,
    FindChirpChisqParams       *params
    );

#ifdef  __cplusplus
}
#endif

#endif _FINDCHIRPCHISQ_H
