/* <lalVerbatim file="LIGOLwXMLReadHV">
$Id$
</lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LIGOLwXMLRead.h}}
\label{s:LIGOLwXMLRead.h}

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LIGOLwXMLRead.h>
\end{verbatim}

</lalLaTeX> */


#ifndef _LIGOLWXMLREAD_H
#define _LIGOLWXMLREAD_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID( LIGOLWXMLREADH, "$Id$" );

/* <lalLaTeX>

\subsection*{Error conditions}
\input{LIGOLwXMLReadHE}

</lalLaTeX> */

/* <lalErrTable file="LIGOLwXMLReadHE"> */
#define LALHELLOH_ENULL  1
#define LALHELLOH_MSGENULL  "Null pointer"
/* </lalErrTable> */

/* <lalLaTeX>
\newpage\input{LIGOLwXMLReadC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _LIGOLWXMLREAD_H */
