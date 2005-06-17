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
