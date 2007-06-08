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
