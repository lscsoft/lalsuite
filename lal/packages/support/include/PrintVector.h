/************************************ <lalVerbatim file="PrintVectorHV">
Author: Allen, B.
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{PrintVector.h}}
\label{s:PrintVector.h}

This is a simple utility to print vectors into a file.

\subsection*{Synopsis}
\begin{verbatim}
#include "PrintVector.h"
\end{verbatim}

\noindent Defines PrintVector interface.

\vfill{\footnotesize\input{PrintVectorHV}}
\newpage\input{PrintVectorC}

</lalLaTeX> */

#ifndef _PRINTVECTOR_H
#define _PRINTVECTOR_H

#include "LALRCSID.h"

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID( PRINTVECTORH, "$Id$" );

void PrintVector( REAL4Vector *vector );

#ifdef  __cplusplus
}
#endif

#endif /* _PRINTVECTOR_H */
