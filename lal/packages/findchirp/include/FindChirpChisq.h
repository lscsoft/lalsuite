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
#define FINDCHIRPCHISQ_ENUMZ 3
#define FINDCHIRPCHISQ_ECHIZ 4
#define FINDCHIRPCHISQ_EALOC 5

#define FINDCHIRPCHISQ_MSGENULL "Null pointer"
#define FINDCHIRPCHISQ_MSGENNUL "Non-null pointer"
#define FINDCHIRPCHISQ_MSGENUMZ "Number of points is zero or negative"
#define FINDCHIRPCHISQ_MSGECHIZ "Number of chisq bins is zero or negative"
#define FINDCHIRPCHISQ_MSGEALOC "Memory allocation error"


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
  REAL4                         chisqNorm;
  UINT4Vector                  *chisqBinVec;
  ComplexFFTPlan               *plan;
  COMPLEX8Vector               *qtildeBinVec;
  COMPLEX8Vector              **qBinVecPtr;
}
FindChirpChisqParams;

void
LALFindChirpChisqVetoInit (
    LALStatus                  *status,
    FindChirpChisqParams       *params,
    UINT4                       numChisqBins,
    UINT4                       numPoints
    );

void
LALFindChirpChisqVetoFinalize (
    LALStatus                  *status,
    FindChirpChisqParams       *params,
    UINT4                       numChisqBins
    );

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
