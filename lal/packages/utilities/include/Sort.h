/*----------------------------------------------------------------------- 
 * 
 * File Name: Sort.h
 * 
 * Author: Creighton, T. D.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------*/

/* <lalLaTeX>

\section{Header \texttt{Sort.h}}

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

\subsection{Error conditions}
\begin{tabular}{|c|l|l|}
\hline
status & status                      & Explanation                      \\
 code  & description                 &                                  \\
\hline
\tt 1  & \tt Null pointer            & Missing a required pointer.      \\
\tt 2  & \tt Length mismatch         & Vectors are of different length. \\
\tt 3  & \tt Memory allocation error & Could not allocate memory.       \\
\hline
\end{tabular}

</lalLaTeX> */

#define SORT_ENUL 1
#define SORT_ELEN 2
#define SORT_EMEM 3

#define SORT_MSGENUL "Null pointer"
#define SORT_MSGELEN "Length mismatch"
#define SORT_MSGEMEM "Memory allocation error"

/* <lalLaTeX>
\subsection{Structures}
</lalLaTeX> */

/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{HeapSortC}
</lalLaTeX> */
void LALSHeapSort(LALStatus      *stat,
	       REAL4Vector *vector);

void LALSHeapIndex(LALStatus      *stat,
		INT4Vector  *indx,
		REAL4Vector *vector);

void LALSHeapRank(LALStatus      *stat,
	       INT4Vector  *rank,
	       REAL4Vector *vector);

void LALDHeapSort(LALStatus      *stat,
	       REAL8Vector *vector);

void LALDHeapIndex(LALStatus      *stat,
		INT4Vector  *indx,
		REAL8Vector *vector);

void LALDHeapRank(LALStatus      *stat,
	       INT4Vector  *rank,
	       REAL8Vector *vector);

/* <lalLaTeX>
\newpage\input{SortTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _SORT_H */
