#if 0 /* autodoc block */

<lalVerbatim file="LALStdioHV">
$Id$
</lalVerbatim>

<lalLaTeX>

\section{Header \texttt{LALStdio.h}}
\label{s:LALStdio.h}

Provides standard LAL IO functions.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
\end{verbatim}

\noindent These headers cover the quasi-LAL IO functions and the LALsnprintf
functions.  Only use \texttt{FileIO.h} in test code that links to
the \texttt{lalsupport} library.

\vfill{\footnotesize\input{LALStdioHV}}
\newpage\input{LALStdC}
\newpage\input{FileIOC}

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
