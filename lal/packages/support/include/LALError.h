/************************************ <lalVerbatim file="LALErrorHV">
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALError.h}}
\label{s:LALError.h}

Provides routines to report errors.  Should not generally be used in production
code.

\subsection*{Synopsis}
\begin{verbatim}
#include "LALError.h"
\end{verbatim}

\noindent This header covers routines that print error messages and abort.

\vfill{\footnotesize\input{LALErrorHV}}
\newpage\input{LALErrorC}

</lalLaTeX> */

#ifndef _LALERROR_H
#define _LALERROR_H

#include "LALRCSID.h"

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID( LALERRORH, "$Id$" );

int
LALPrintError( const char *fmt, ... );

void
LALAbort( const char *msg );


#ifdef  __cplusplus
}
#endif

#endif /* _LALERROR_H */
