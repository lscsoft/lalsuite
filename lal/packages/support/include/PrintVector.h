/************************************ <lalVerbatim file="PrintVectorHV">
Author: Allen, B.; generalized by J.T. Whelan
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{PrintVector.h}}
\label{s:PrintVector.h}

This is a simple utility to print vectors into a file.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/PrintVector.h>
\end{verbatim}

\noindent Contains the prototypes for the
LAL$\langle\mbox{DT}\rangle$PrintVector functions

\vfill{\footnotesize\input{PrintVectorHV}}
\newpage\input{PrintVectorC}
\newpage\input{PrintVectorTestC}

</lalLaTeX> */

#ifndef _PRINTVECTOR_H
#define _PRINTVECTOR_H

#include <lal/LALRCSID.h>

#ifndef _LALSTDLIB_H
#include <lal/LALStdlib.h>
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID( PRINTVECTORH, "$Id$" );

void LALCHARPrintVector( CHARVector *vector );
void LALI2PrintVector( INT2Vector *vector );
void LALI4PrintVector( INT4Vector *vector );
void LALI8PrintVector( INT8Vector *vector );
void LALU2PrintVector( UINT2Vector *vector );
void LALU4PrintVector( UINT4Vector *vector );
void LALU8PrintVector( UINT8Vector *vector );
void LALPrintVector( REAL4Vector *vector );
void LALSPrintVector( REAL4Vector *vector );
void LALDPrintVector( REAL8Vector *vector );
void LALCPrintVector( COMPLEX8Vector *vector );
void LALZPrintVector( COMPLEX16Vector *vector );

#ifdef  __cplusplus
}
#endif

#endif /* _PRINTVECTOR_H */
