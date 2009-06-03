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

/************************************ <lalVerbatim file="LALStdlibHV">
Author: J. D. E. Creighton, T. D. Creighton
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALConfig.h}}
\label{s:LALConfig.h}

Defines configuration macro constants.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALConfig.h>
\end{verbatim}

\noindent This header (which is not technically in the \texttt{std} package;
rather it is generated directly in the \texttt{include/lal} directory during
configuration) is included in essentially every other header file.  It
contains macro constants that are defined at configuration time.  They are:

\idx[Constant]{LAL\_VERSION}
\idx[Constant]{LAL\_VERSION\_MAJOR}
\idx[Constant]{LAL\_VERSION\_MINOR}
\idx[Constant]{LAL\_CONFIGURE\_ARGS}
\idx[Constant]{LAL\_CONFIGURE\_DATE}
\idx[Constant]{LAL\_CVS\_TAG}
\idx[Constant]{LAL\_SIZEOF\_DOUBLE}
\idx[Constant]{LAL\_SIZEOF\_FLOAT}
\idx[Constant]{LAL\_SIZEOF\_INT}
\idx[Constant]{LAL\_SIZEOF\_LONG}
\idx[Constant]{LAL\_SIZEOF\_LONG\_LONG}
\idx[Constant]{LAL\_SIZEOF\_SHORT}
\idx[Constant]{LAL\_NDEBUG}
\idx[Constant]{NOLALMACROS}
\idx[Constant]{LAL\_PTHREAD\_LOCK}
\idx[Constant]{LAL\_FRAME\_ENABLED}
\idx[Constant]{LAL\_MPI\_ENABLED}
\begin{description}
\item[\texttt{LAL\_VERSION}] Constant string containing the version of LAL.
\item[\texttt{LAL\_VERSION\_MAJOR}] Integer representing the major version
  number of LAL.
\item[\texttt{LAL\_VERSION\_MINOR}] Integer representing the minor version
  number of LAL.
\item[\texttt{LAL\_CONFIGURE\_ARGS}] Constant string containing the arguments
  given to the configure script.
\item[\texttt{LAL\_CONFIGURE\_DATE}] Constant string containing the date
  when LAL was configured.
\item[\texttt{LAL\_CVS\_TAG}] Constant string containing the CVS tag used to
  checkout LAL (blank if none).
\item[\texttt{LAL\_SIZEOF\_DOUBLE}] Integer representing the size of a
  double precision floating point number.
\item[\texttt{LAL\_SIZEOF\_FLOAT}] Integer representing the size of a
  single precision floating point number.
\item[\texttt{LAL\_SIZEOF\_INT}] Integer representing the size of a
  integer.
\item[\texttt{LAL\_SIZEOF\_LONG}] Integer representing the size of a
  long integer.
\item[\texttt{LAL\_SIZEOF\_LONG\_LONG}] Integer representing the size of a
  long long integer.
\item[\texttt{LAL\_SIZEOF\_SHORT}] Integer representing the size of a
  short integer.
\item[\texttt{LAL\_NDEBUG}] Defined if debugging is turned off (use the
  configure argument \texttt{--disable-debug} to do this).
\item[\texttt{NOLALMACROS}] Defined if status macros are replaced by functions
  (where possible) (use the configure argument \texttt{--disable-macros} to do
  this).
\item[\texttt{LAL\_PTHREAD\_LOCK}] Defined if POSIX thread mutex locking is
  to be used for threadsafety (use the configure argument
  \texttt{--enable-pthread-lock} to do this).
\item[\texttt{LAL\_FRAME\_ENABLED}] Defined if LAL frame-format data reading
  routines will be compiled (use the configure argument
  \texttt{--enable-frame} to do this).
\item[\texttt{LAL\_MPI\_ENABLED}] Defined if LAL MPI routines will be compiled
  (use the configure argument \texttt{--enable-mpi} to do this).
\end{description}


\newpage
\section{Header \texttt{LALStdlib.h}}
\label{s:LALStdlib.h}

Includes the standard LAL header files.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALStdlib.h>
\end{verbatim}

\noindent This header is the overall header for the \verb@std@
package.  It provides the datatypes, constants, and macros required by
most LAL functions, by including the following header files in the
\verb@std@ package:

\vspace{1ex}

</lalLaTeX> */

#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H

/* <lalVerbatim> */
#include <lal/LALRCSID.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStatusMacros.h>
/* </lalVerbatim>
<lalLaTeX>

\noindent\verb@LALStdlib.h@ also includes function prototype headers
for certain standard modules used by many LAL routines:

\vspace{1ex}

</lalLaTeX>
<lalVerbatim> */
#include <stdio.h>
#include <stdarg.h>
#include <lal/LALMalloc.h>
/* </lalVerbatim>

<lalLaTeX>
\vfill{\footnotesize\input{LALStdlibHV}}
</lalLaTeX> */

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (LALSTDLIBH, "$Id$");

/* These are non-ANSI standard routines that will be allowed in LAL */
int getopt( int, char * const *, const char * );
FILE *popen( const char *, const char * );
int pclose( FILE * );

#ifdef  __cplusplus
}
#endif

#endif /* _LALSTDLIB_H */
