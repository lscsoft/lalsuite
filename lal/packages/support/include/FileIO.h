/* <lalVerbatim file="FileIOHV">
$Id$
</lalVerbatim> */

/* <lalLaTeX>
\section{Header \texttt{FileIO.h}}
\label{s:FileIO.h}

Provides standard LAL support IO functions.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
\end{verbatim}
\noindent Only use \texttt{FileIO.h} in test code that links to
the \texttt{lalsupport} library.

\vfill{\footnotesize\input{FileIOHV}}
\newpage\input{FileIOC}
</lalLaTeX> */ 

#ifndef _FILEIO_H
#define _FILEIO_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <lal/LALRCSID.h>

NRCSID( FILEIOH, "$Id$" );

#ifndef LALFopen
#define LALFopen fopen
#endif

#ifndef LALFclose
#define LALFclose fclose
#endif

FILE *
LALOpenDataFile( const char * );

#ifdef __cplusplus
}
#endif

#endif /* _FILEIO_H */
