/************************************ <lalVerbatim file="LALStdlibHV">
Author: Finn, L. S.
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

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
#include <lal/LALMalloc.h>
/* </lalVerbatim> 

<lalLaTeX>
\vfill{\footnotesize\input{LALStdlibHV}}
</lalLaTeX> */

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (LALSTDLIBH, "$Id$");

#define LALFopen  fopen    /* will be replaced by custom unit */
#define LALFclose fclose   /* will be replaced by custom unit */

#ifdef  __cplusplus
}
#endif

#endif /* _LALSTDLIB_H */
