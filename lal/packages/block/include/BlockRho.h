/****************** <lalVerbatim file="BlockRhoHV">
Author: Tibbits, M. M.
$Id$
********************************* </lalVerbatim> */

/******************************* <lalLaTeX file="BlockRhoH">
\section{Header \texttt{BlockRho.h}}
\label{s:BlockRho.h}

Provides structures and definitions global to the Rho algorithms,
notably \texttt{BlockSearchParams}.

\subsection*{Error conditions}

\input{BlockRhoHErrTab}

\subsection*{Structures and Unions}

*************************************************** </lalLaTeX> */

#ifndef _BLOCK_RHO_H_
#define _BLOCK_RHO_H_

#include <math.h>
#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALConstants.h>
#include <lal/LALMoment.h>
#include <lal/Comm.h>
#include <lal/Matrix.h>

/* C++ protection */
#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( BLOCKRHOH, "$Id$");


/******** <lalErrTable file="BlockRhoHErrTab"> ********/
#define BLOCKRHOH_ENULL 1
#define BLOCKRHOH_ENNUL 2
#define BLOCKRHOH_EALOC 3
#define BLOCKRHOH_ENUMZ 4
#define	BLOCKRHOH_EARG  5
#define BLOCKRHOH_EDATA 6

#define BLOCKRHOH_MSGENULL "Null pointer"
#define BLOCKRHOH_MSGENNUL "Non-null pointer"
#define BLOCKRHOH_MSGEALOC "Memory allocation error"
#define BLOCKRHOH_MSGENUMZ "Data segment length is zero"
#define	BLOCKRHOH_MSGEARG  "Error parsing command-line arguments"
#define BLOCKRHOH_MSGEDATA "Too few input data points to define a Rho statistic"

/******** </lalErrTable> ********/

void LALBlockRho2 (
 LALStatus		*status,
 REAL8			*result,
 REAL8			*rpeak,
 INT4			*myindex,
 REAL8Sequence		*data,
 UINT4			*marginOfExclusion			
);

void LALBlockRho3 (
 LALStatus               *status,
 REAL8                   *result,
 REAL8                   *rpeak,
 INT4                    *myindex,
 REAL8Sequence           *data
);


/* C++ protection */
#ifdef  __cplusplus
}
#endif

/******************************* <lalLaTeX file="BlockRhoH">
  \vfill{\footnotesize\input{BlockRhoHV}}
  \newpage\input{BlockRho2C}
  %\newpage\input{BlockRho3C}
*************************************************** </lalLaTeX> */

#endif /* _BLOCK_RHO_H_ */
