#if 0 /* autodoc block */

<lalVerbatim file="LALStdioHV">
$Id$
</lalVerbatim>

<lalLaTeX>

\section{Header \texttt{LALStdio.h}}
\label{s:LALStdio.h}

Provides LAL functions similar to the non-file functions in \verb+<stdio.h>+.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
\end{verbatim}

\noindent This header provides the LALsnprintf function.

\vfill{\footnotesize\input{LALStdioHV}}
\newpage\input{LALStdC}
</lalLaTeX>
#endif /* autodoc block */

#ifndef _LALSTDIO_H
#define _LALSTDIO_H

#include <stdio.h>
#include <stdarg.h>
#include <lal/LALRCSID.h>

#ifdef __cplusplus
extern "C" {
#endif

NRCSID( LALSTDIOH, "$Id$" );

#define LALFopen fopen
#define LALFclose fclose

int
LALSnprintf( char *, size_t, const char *, ... );

int
LALVsnprintf( char *, size_t, const char *, va_list );

#ifdef __cplusplus
}
#endif

#endif /* _LALSTDIO_H */
