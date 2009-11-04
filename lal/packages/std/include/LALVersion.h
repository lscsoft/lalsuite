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

<lalVerbatim file="LALVersionHV">
$Id$
</lalVerbatim>

<lalLaTeX>

\section{Header \texttt{LALVersion.h}}
\label{s:LALVersion.h}

Provides routines for reporting the LAL version.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALVersion.h>
\end{verbatim}

\noindent This header covers the routines for reporting the LAL version.


\subsection*{Global variables}
\idx[Constant]{lalVersion}
\idx[Constant]{lalVersionMajor}
\idx[Constant]{lalVersionMinor}
\idx[Constant]{lalConfigureArgs}
\idx[Constant]{lalConfigureDate}
\idx[Constant]{lalCVSTag}
\begin{verbatim}
extern const char *lalVersion;
extern const int   lalVersionMajor;
extern const int   lalVersionMinor;
extern const char *lalConfigureArgs;
extern const char *lalConfigureDate;
extern const char *lalCVSTag;
\end{verbatim}

These constant variables are set at compile time and included into the LAL
library.  They contain the information about the version of LAL and the
configuration information.

\subsection*{Error conditions}
\input{LALVersionHErrTab}

\newpage\input{LALVersionC}
\newpage\input{LALVersionTestC}

</lalLaTeX>
#endif /* autodoc block */

#ifndef _LALVERSION_H
#define _LALVERSION_H

#include <lal/LALDatatypes.h>

#ifdef __cplusplus
extern "C" {
#endif

NRCSID( LALVERSIONH, "$Id$" );

/* <lalErrTable file="LALVersionHErrTab"> */

#define LALVERSIONH_ENULL 1
#define LALVERSIONH_ESIZE 2
#define LALVERSIONH_ESPRN 4
#define LALVERSIONH_ESHRT 8

#define LALVERSIONH_MSGENULL "Null string pointer."
#define LALVERSIONH_MSGESIZE "Zero string size."
#define LALVERSIONH_MSGESPRN "Error in snprintf."
#define LALVERSIONH_MSGESHRT "String too short."

/* </lalErrTable> */

extern const char *lalVersion;
extern const int   lalVersionMajor;
extern const int   lalVersionMinor;
extern const int   lalVersionMicro;
extern const char *lalBuildDate;
extern const char *lalConfigureArgs;
extern const char *lalConfigureDate;
extern const char *lalCVSTag;


void
LALVersion( LALStatus *status, CHAR *message, UINT4 size, INT4 verbose );

#ifdef __cplusplus
}
#endif

#endif /* _LALVERSION_H */
