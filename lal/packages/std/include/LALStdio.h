/*
*  Copyright (C) 2007 Jolien Creighton
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

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
