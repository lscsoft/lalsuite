/*-----------------------------------------------------------------------

File Name: VectorFactories.c

<lalVerbatim file="VectorFactoriesCV">
Revision: $Id$
</lalVerbatim>

-------------------------------------------------------------------------*/

/* <lalLaTeX>

\subsection{Module \texttt{VectorFactories.c}}
\label{ss:VectorFactories.c}

Create/destroy $\langle\mbox{datatype}\rangle$Vector objects.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{VectorFactoriesD}
\idx{LALZCreateVector()}
\idx{LALCCreateVector()}
\idx{LALDCreateVector()}
\idx{LALSCreateVector()}
\idx{LALI2CreateVector()}
\idx{LALI4CreateVector()}
\idx{LALI8CreateVector()}
\idx{LALU2CreateVector()}
\idx{LALU4CreateVector()}
\idx{LALU8CreateVector()}
\idx{LALCHARCreateVector()}
\idx{LALCreateVector()}
\idx{LALZDestroyVector()}
\idx{LALCDestroyVector()}
\idx{LALDDestroyVector()}
\idx{LALSDestroyVector()}
\idx{LALI2DestroyVector()}
\idx{LALI4DestroyVector()}
\idx{LALI8DestroyVector()}
\idx{LALU2DestroyVector()}
\idx{LALU4DestroyVector()}
\idx{LALU8DestroyVector()}
\idx{LALCHARDestroyVector()}
\idx{LALDestroyVector()}

\subsubsection*{Description}

The \texttt{CreateVector} family of functions create a
$\langle\mbox{datatype}\rangle$\texttt{Vector} of the appropriate
dimensions.

The \texttt{ResizeVector} family of functions changes the amount of
storage allocated by the \texttt{CreateVector} functions.

The \texttt{DestroyVector} family of functions return the storage allocated by
the \texttt{CreateVector} functions to the system.

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()
LALFree()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{VectorFactoriesCV}}

</lalLaTeX> */

#include "LALStdlib.h"
#include "AVFactories.h"

/* <lalVerbatim file="VectorFactoriesNRCSID"> */
NRCSID( VECTORFACTORIESC, "$Id$" );
/* </lalVerbatim> */

define(`TYPECODE',`Z')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`C')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`D')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`S')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`I2')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`I4')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`I8')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`U2')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`U4')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`U8')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`CHAR')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')
