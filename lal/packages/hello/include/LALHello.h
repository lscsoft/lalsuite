/************************************ <lalVerbatim file="LALHelloHV">
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALHello.h}}
\label{s:LALHello.h}

Provides routines to print ``hello, LSC!''

\subsection*{Synopsis}
\begin{verbatim}
#include "LALHello.h"
\end{verbatim}

\noindent This header covers the routine to print the greeting message.

</lalLaTeX> */


#ifndef _LALHELLO_H
#define _LALHELLO_H

#include "LALDatatypes.h"

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (LALHELLOH, "$Id$");

/* <lalLaTeX>

\subsection*{Error conditions}
\input{LALHelloHErrTab}

</lalLaTeX> */

/*
<lalErrTable file="LALHelloHErrTab">
*/

#define LALHELLOH_ENMSG  1
#define LALHELLOH_EOPEN  2
#define LALHELLOH_EWRITE 4
#define LALHELLOH_ECLOSE 8
#define LALHELLOH_EFLUSH 16

#define LALHELLOH_MSGENMSG  "Null output message string."
#define LALHELLOH_MSGEOPEN  "Could not open file."
#define LALHELLOH_MSGEWRITE "Error in writing to file."
#define LALHELLOH_MSGECLOSE "Error in closing file."
#define LALHELLOH_MSGEFLUSH "Error in flushing stdout."

/*
</lalErrTable>
*/


/* Structures. */

/* <lalLaTeX>
\subsection*{Structures}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{LALHelloC}
</lalLaTeX> */

void
LALHello( Status *status, const CHAR *fileName );


/* Test program. */

/* <lalLaTeX>
\newpage\input{LALHelloTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _LALHELLO_H */
