/********************************** <lalVerbatim file="LALInitBarycenterHV">
Author: Cutler, C.
$Id$
*********************************** </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALInitBarycenter.h}}
\label{s:LALInitBarycenter.h}

Provides a routine for reading Earth and Sun position information from
data files.
 

\subsection*{Synopsis}

\begin{verbatim}
#include "LALInitBarycenter.h"
\end{verbatim}

\noindent This header covers the routine \verb@LALInitBarycenter.c@.
Since it involves file I/O, it is placed in the \verb@support@
package, and included in the \verb@lalsupport@ library.

</lalLaTeX> */



#ifndef _LALINITBARYCENTER_H    /* Protect against double-inclusion */
#define _LALINITBARYCENTER_H

#include <lal/LALBarycenter.h>
#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (LALINITBARYCENTERH,"$Id$");

/* <lalErrTable file="LALInitBarycenterHErrorTable"> */
#define LALINITBARYCENTERH_EOPEN 1

#define LALINITBARYCENTERH_EEPHFILE 32

#define LALINITBARYCENTERH_MSGEOPEN "Could not open ephemeris file."
#define LALINITBARYCENTERH_MSGEEPHFILE "Error in reading an ephemeris file."
/* </lalErrTable> */

/* <lalLaTeX>
\subsection*{Error conditions}
\vspace{0.1in}
\input{LALInitBarycenterHErrorTable}
</lalLaTeX> */

/* Function prototypes. */

void LALInitBarycenter(LALStatus *, EphemerisData *);


/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{LALInitBarycenterHV}}
\newpage\input{LALInitBarycenterC}
******************************************************* </lalLaTeX> */

#ifdef  __cplusplus
}
#endif      /* Close C++ protection */

#endif      /* Close double-include protection */























