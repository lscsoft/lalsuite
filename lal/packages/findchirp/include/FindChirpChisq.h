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

#if 0
<lalVerbatim file="FindChirpChisqHV">
Author: Allen, B., Brown, D. A. and Creighton, J. D. E.
$Id$
</lalVerbatim> 

<lalLaTeX>

\section{Header \texttt{FindChirpChisq.h}}
\label{s:FindChirp.h}

Provides routines to perform chisq veton on binary inspiral chirps.

</lalLaTeX>
#endif

#ifndef _FINDCHIRPCHISQH_H
#define _FINDCHIRPCHISQH_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (SPFINDCHIRPCHISQH, "$Id$");

#if 0
<lalLaTeX> 
\subsection*{Error codes} 
</lalLaTeX>
#endif
/* <lalErrTable> */
#define FINDCHIRPCHISQH_ENULL 1
#define FINDCHIRPCHISQH_ENNUL 2
#define FINDCHIRPCHISQH_ENUMZ 3
#define FINDCHIRPCHISQH_ECHIZ 4
#define FINDCHIRPCHISQH_EALOC 5
#define FINDCHIRPCHISQH_MSGENULL "Null pointer"
#define FINDCHIRPCHISQH_MSGENNUL "Non-null pointer"
#define FINDCHIRPCHISQH_MSGENUMZ "Number of points is zero or negative"
#define FINDCHIRPCHISQH_MSGECHIZ "Number of chisq bins is zero or negative"
#define FINDCHIRPCHISQH_MSGEALOC "Memory allocation error"
/* </lalErrTable> */


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

#endif /* _FINDCHIRPCHISQH_H */
