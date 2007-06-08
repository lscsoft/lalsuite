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

/************************************ <lalVerbatim file="LALSampleHV">
Author: Creighton, T. D.
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\section{Header \texttt{LALSample.h}}
\label{s:LALSample.h}

Example header for LAL.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALSample.h>
\end{verbatim}

\noindent This header provides two trivial functions to divide real
numbers.  It exists primarily to demonstrate documentation and coding
standards for LAL headers.

******************************************************* </lalLaTeX> */

#ifndef _LALSAMPLE_H  /* Double-include protection. */
#define _LALSAMPLE_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID( LALSAMPLEH, "$Id$" );

/********************************************************** <lalLaTeX>
\subsection*{Error codes}
</lalLaTeX>
***************************************************** <lalErrTable> */
#define LALSAMPLEH_ENULL 1
#define LALSAMPLEH_EDIV0 2
#define LALSAMPLEH_MSGENULL "Arguments contained an unexpected null pointer"
#define LALSAMPLEH_MSGEDIV0 "Division by zero"
/*************************************************** </lalErrTable> */

/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{LALSampleHV}}
\newpage\input{LALSampleC}
******************************************************* </lalLaTeX> */

/* Function prototypes */

void
LALREAL8Invert( LALStatus *status, REAL8 *output, REAL8 input );

void
LALREAL8Divide( LALStatus *status, REAL8 *output, REAL8 numer, REAL8 denom);

/********************************************************** <lalLaTeX>
\newpage\input{LALSampleTestC}
******************************************************* </lalLaTeX> */

#ifdef  __cplusplus
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. */
