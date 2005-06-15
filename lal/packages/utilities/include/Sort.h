/*----------------------------------------------------------------------- 
 * 
 * File Name: Sort.h
 * 
 * Author: Creighton, T. D.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------*/

/*

<lalVerbatim file="SortHV">
$Id$
</lalVerbatim>

<lalLaTeX>

\section{Header \texttt{Sort.h}}
\label{s:Sort.h}

Provides routines for sorting, indexing, and ranking real vector
elements.

\subsection{Synopsis}
\begin{verbatim}
#include <lal/Sort.h>
\end{verbatim}

</lalLaTeX> */

#ifndef _SORT_H
#define _SORT_H

#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID(SORTH,"$Id$");

/* <lalLaTeX>

\subsection*{Error conditions}
\input{SortHErrTab}
</lalLaTeX>

<lalErrTable file="SortHErrTab"> */

#define SORTH_ENUL 1
#define SORTH_ELEN 2
#define SORTH_EMEM 3

#define SORTH_MSGENUL "Null pointer"
#define SORTH_MSGELEN "Length mismatch"
#define SORTH_MSGEMEM "Memory allocation error"

/* </lalErrTable>

<lalLaTeX>
\subsection*{Structures}
</lalLaTeX> */

/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{HeapSortC}
</lalLaTeX> */
void LALSHeapSort(LALStatus      *status,
	       REAL4Vector *vector);

void LALSHeapIndex(LALStatus      *status,
		INT4Vector  *indx,
		REAL4Vector *vector);

void LALSHeapRank(LALStatus      *status,
	       INT4Vector  *rank,
	       REAL4Vector *vector);

void LALDHeapSort(LALStatus      *status,
	       REAL8Vector *vector);

void LALDHeapIndex(LALStatus      *status,
		INT4Vector  *indx,
		REAL8Vector *vector);

void LALDHeapRank(LALStatus      *status,
	       INT4Vector  *rank,
	       REAL8Vector *vector);

/* <lalLaTeX>
\newpage\input{SortTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _SORT_H */
