/************************************ <lalVerbatim file="LALRunningMedianHV">
Author: B. Machenschalk
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALRunningMedian.h}}
\label{s:LALRunningMedian.h}

Provides routines to efficiently calculate the running median

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALRunningMedian.h>
\end{verbatim}

\noindent This header covers routines to efficiently calculate the
running median of REAL4 and REAL8 sequences

</lalLaTeX> */


#ifndef _LALRUNNINGMEDIAN_H
#define _LALRUNNINGMEDIAN_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (LALRUNNINGMEDIANH, "$Id$");

/* <lalLaTeX>

\subsection*{Error conditions}
\input{LALRunningMedianHErrTab}

</lalLaTeX> */

/*
<lalErrTable file="LALRunningMedianHErrTab">
*/

#define LALRUNNINGMEDIANH_EMALOC1 1
#define LALRUNNINGMEDIANH_EMALOC2 2
#define LALRUNNINGMEDIANH_EMALOC3 3
#define LALRUNNINGMEDIANH_EMALOC4 4
#define LALRUNNINGMEDIANH_EMALOC5 5
#define LALRUNNINGMEDIANH_EMALOC6 6
#define LALRUNNINGMEDIANH_ECV 7
#define LALRUNNINGMEDIANH_ENULL 8
#define LALRUNNINGMEDIANH_EZERO 9
#define LALRUNNINGMEDIANH_ELARGE 10
#define LALRUNNINGMEDIANH_EIMED 11

#define LALRUNNINGMEDIANH_MSGEMALOC1 "Could not allocate indexblock"
#define LALRUNNINGMEDIANH_MSGEMALOC2 "Could not allocate checks"
#define LALRUNNINGMEDIANH_MSGEMALOC3 "Could not allocate checks4shift"
#define LALRUNNINGMEDIANH_MSGEMALOC4 "Could not allocate nodeaddresses"
#define LALRUNNINGMEDIANH_MSGEMALOC5 "Could not aloocate first node"
#define LALRUNNINGMEDIANH_MSGEMALOC6 "Could not allocate node"
#define LALRUNNINGMEDIANH_MSGECV     "Could not create output vector (LALCreateVector failed)"
#define LALRUNNINGMEDIANH_MSGENULL   "Invalid input: NULL pointer."
#define LALRUNNINGMEDIANH_MSGEZERO   "Invalid input: block length must be >2"
#define LALRUNNINGMEDIANH_MSGELARGE  "Invalid input: block length larger than imput length"
#define LALRUNNINGMEDIANH_MSGEIMED   "Invalid input: wrong size of median array"

/*
</lalErrTable>
*/


/* Structures. */

/* <lalLaTeX>
\subsection*{Structures}
This is the parameter structure for the LALRunningMedian functions.
Currently the only parameter supported is the blocksize, the number
of elements a single median is calculated from. The current
implementation requires the blocksize to be $< 2$.
\begin{verbatim}
typedef struct tagLALRunningMedianPar
{
  UINT4 blocksize;
}
LALRunningMedianPar;
\end{verbatim}
</lalLaTeX> */

typedef struct tagLALRunningMedianPar
{
  UINT4 blocksize;
}
LALRunningMedianPar;


/* Function prototypes. */

void
LALDRunningMedian( LALStatus *status,
		   REAL8Sequence *medians,
		   const REAL8Sequence *input,
		   LALRunningMedianPar param);

void
LALSRunningMedian( LALStatus *status,
		   REAL4Sequence *medians,
		   const REAL4Sequence *input,
		   LALRunningMedianPar param);

void
LALDRunningMedian2( LALStatus *status,
		    REAL8Sequence *medians,
		    const REAL8Sequence *input,
		    LALRunningMedianPar param);

/* <lalLaTeX>
\newpage\input{LALRunningMedianC}
</lalLaTeX> */


/* Test program. */
/* <lalLaTeX>
\newpage\input{LALRunningMedianTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _LALRUNNINGMEDIAN_H */
