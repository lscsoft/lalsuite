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
\index{\verb&ZCreateArray()&}
\index{\verb&CCreateArray()&}
\index{\verb&DCreateArray()&}
\index{\verb&SCreateArray()&}
\index{\verb&I2CreateArray()&}
\index{\verb&I4CreateArray()&}
\index{\verb&I8CreateArray()&}
\index{\verb&U2CreateArray()&}
\index{\verb&U4CreateArray()&}
\index{\verb&U8CreateArray()&}
\index{\verb&CreateArray()&}
\index{\verb&ZDestroyArray()&}
\index{\verb&CDestroyArray()&}
\index{\verb&DDestroyArray()&}
\index{\verb&SDestroyArray()&}
\index{\verb&I2DestroyArray()&}
\index{\verb&I4DestroyArray()&}
\index{\verb&I8DestroyArray()&}
\index{\verb&U2DestroyArray()&}
\index{\verb&U4DestroyArray()&}
\index{\verb&U8DestroyArray()&}
\index{\verb&DestroyArray()&}

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
