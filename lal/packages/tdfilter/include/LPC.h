#ifndef _LPC_H
#define _LPC_H

#include <lal/LALDatatypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <lal/LALMalloc.h>

/*<lalLaTeX file="LPCH">
\section{Header \texttt{LPC.h}}
\label{s:LPC.h}

Functions for linear predictor filters.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LPC.h>
\end{verbatim}

\subsection*{Error Conditions}
\input{LPCHErrTab}
</lalLaTeX>*/

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif


NRCSID (LPCH, "$Id$");

  /******************************** <lalErrTable file="LPCHErrTab"> */
#define LPCH_EIN 1
#define LPCH_EMEM 2
#define LPCH_ENUM 3

#define LPCH_MSGEIN "invalid input"
#define LPCH_MSGEMEM "memory error"
#define LPCH_MSGENUM "numerical error"
/*************************************************** </lalErrTable> */

void LALPolystab(LALStatus *status,
		 REAL4Vector *a);

void LALLPC(LALStatus *status,
	    REAL4Vector *aout,    /* returned filter coefficients */
	    REAL4Vector *x,    /* training data */
	    UINT4 p            /* filter order */
	    );

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif
