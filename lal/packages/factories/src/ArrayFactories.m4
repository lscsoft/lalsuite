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
\index{\texttt{LALZCreateArray()}}
\index{\texttt{LALCCreateArray()}}
\index{\texttt{LALDCreateArray()}}
\index{\texttt{LALSCreateArray()}}
\index{\texttt{LALI2CreateArray()}}
\index{\texttt{LALI4CreateArray()}}
\index{\texttt{LALI8CreateArray()}}
\index{\texttt{LALU2CreateArray()}}
\index{\texttt{LALU4CreateArray()}}
\index{\texttt{LALU8CreateArray()}}
\index{\texttt{LALCreateArray()}}
\index{\texttt{LALZDestroyArray()}}
\index{\texttt{LALCDestroyArray()}}
\index{\texttt{LALDDestroyArray()}}
\index{\texttt{LALSDestroyArray()}}
\index{\texttt{LALI2DestroyArray()}}
\index{\texttt{LALI4DestroyArray()}}
\index{\texttt{LALI8DestroyArray()}}
\index{\texttt{LALU2DestroyArray()}}
\index{\texttt{LALU4DestroyArray()}}
\index{\texttt{LALU8DestroyArray()}}
\index{\texttt{LALDestroyArray()}}

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
include(`DestroyArray.m4')

define(`TYPECODE',`C')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`D')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`S')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`I2')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`I4')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`I8')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`U2')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`U4')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`U8')
include(`CreateArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`')
include(`CreateArray.m4')
include(`DestroyArray.m4')
