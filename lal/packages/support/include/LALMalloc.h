/************************************ <lalVerbatim file="LALMallocHV">
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALMalloc.h}}
\label{s:LALMalloc.h}

Provides standard LAL memory allocation/deallocation routines.

\subsection*{Synopsis}
\begin{verbatim}
#include "LALMalloc.h"
\end{verbatim}

\noindent This header covers routines that replace the standard \verb+malloc+,
\verb+calloc+, \verb+realloc+, and \verb+free+.  All memory allocation and
deallocation in LAL should use these replacement functions.

\vfill{\footnotesize\input{LALMallocHV}}
\newpage\input{LALMallocC}
\newpage\input{LALMallocTestC}

</lalLaTeX> */



#ifndef _LALMALLOC_H
#define _LALMALLOC_H

#include "LALRCSID.h"

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID( LALMALLOCH, "$Id$" );

void *
LALMalloc( size_t n );

void
LALFree( void *p );

void *
LALCalloc( size_t m, size_t n );

void *
LALRealloc( void *p, size_t n );

void
LALCheckMemoryLeaks( void );


#ifdef  __cplusplus
}
#endif

#endif /* _LALMALLOC_H */
