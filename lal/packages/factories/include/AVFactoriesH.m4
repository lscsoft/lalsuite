/*----------------------------------------------------------------------- 

File Name: AVFactories.h

<lalVerbatim file="AVFactoriesHV">
Revision: $Id$
</lalVerbatim>

-------------------------------------------------------------------------*/

/* <lalLaTeX>

\section{Header \texttt{AVFactories.h}}
\label{s:AVFactories.h}

Provides prototype and status code information for use of CreateVector,
CreateArray, ResizeVector, ResizeArray, DestroyVector and DestroyArray

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/AVFactories.h>
\end{verbatim}

</lalLaTeX> */

#ifndef _AVFACTORIES_H
#define _AVFACTORIES_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (AVFACTORIESH, "$Id$");

/* <lalLaTeX>

\subsection*{Error conditions}
\input{AVFactoriesHErrTab}

</lalLaTeX> */

/*
<lalErrTable file="AVFactoriesHErrTab">
*/
#define AVFACTORIESH_ELENGTH 1
#define AVFACTORIESH_EVPTR   2
#define AVFACTORIESH_EUPTR   4
#define AVFACTORIESH_EDPTR   8
#define AVFACTORIESH_EMALLOC 16
#define AVFACTORIESH_MSGELENGTH  "Illegal length."
#define AVFACTORIESH_MSGEVPTR    "Null vector/array handle."
#define AVFACTORIESH_MSGEUPTR    "Non-null vector/array pointer."
#define AVFACTORIESH_MSGEDPTR    "Null vector/array data."
#define AVFACTORIESH_MSGEMALLOC  "Malloc failure."
/*
</lalErrTable>
*/

/* Function prototypes. */
/* <lalLaTeX>
\newpage\input{VectorFactoriesC}
\newpage\input{ArrayFactoriesC}
</lalLaTeX> */

dnl There is no CHARArray type
define(`TYPECODE',`CHAR')
include(`VFactoriesBaseH.m4')

define(`TYPECODE',`I2')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`I4')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`I8')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`U2')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`U4')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`U8')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`S')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`D')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`C')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`Z')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

define(`TYPECODE',`')
include(`VFactoriesBaseH.m4')
include(`AFactoriesBaseH.m4')

/* Test program. */

/* <lalLaTeX>
\newpage\input{VectorFactoriesTestC}
\newpage\input{ArrayFactoriesTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _AVFACTORIES_H */
