/*-----------------------------------------------------------------------

File Name: ArrayFactories.c

<lalVerbatim file="ArrayFactoriesCV">
Revision: $Id$
</lalVerbatim>

-------------------------------------------------------------------------*/

/* <lalLaTeX>

\subsection{Module \texttt{ArrayFactories.c}}
\label{ss:ArrayFactories.c}

Create/destroy $\langle\mbox{datatype}\rangle$Array objects.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{ArrayFactoriesD}
\idx{LALZCreateArray()}
\idx{LALCCreateArray()}
\idx{LALDCreateArray()}
\idx{LALSCreateArray()}
\idx{LALI2CreateArray()}
\idx{LALI4CreateArray()}
\idx{LALI8CreateArray()}
\idx{LALU2CreateArray()}
\idx{LALU4CreateArray()}
\idx{LALU8CreateArray()}
\idx{LALCreateArray()}
\idx{LALZDestroyArray()}
\idx{LALCDestroyArray()}
\idx{LALDDestroyArray()}
\idx{LALSDestroyArray()}
\idx{LALI2DestroyArray()}
\idx{LALI4DestroyArray()}
\idx{LALI8DestroyArray()}
\idx{LALU2DestroyArray()}
\idx{LALU4DestroyArray()}
\idx{LALU8DestroyArray()}
\idx{LALDestroyArray()}

\subsubsection*{Description}

The \texttt{CreateArray} family of functions create a
$\langle\mbox{datatype}\rangle$\texttt{Array} of the appropriate dimensions.

The \texttt{DestroyArray} family of functions return the storage allocated by
the \texttt{CreateArray} functions to the system.

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()
LALFree()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{ArrayFactoriesCV}}

</lalLaTeX> */

#include <string.h>
#include "LALStdlib.h"
#include "AVFactories.h"

/* <lalVerbatim file="ArrayFactoriesNRCSID"> */
NRCSID( ARRAYFACTORIESC, "$Id$" );
/* </lalVerbatim> */

define(`TYPECODE',`Z')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`C')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`D')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`S')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`I2')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`I4')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`I8')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`U2')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`U4')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`U8')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')
