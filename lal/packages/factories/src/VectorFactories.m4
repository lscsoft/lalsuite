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
\index{\texttt{LALZCreateVector()}}
\index{\texttt{LALCCreateVector()}}
\index{\texttt{LALDCreateVector()}}
\index{\texttt{LALSCreateVector()}}
\index{\texttt{LALI2CreateVector()}}
\index{\texttt{LALI4CreateVector()}}
\index{\texttt{LALI8CreateVector()}}
\index{\texttt{LALU2CreateVector()}}
\index{\texttt{LALU4CreateVector()}}
\index{\texttt{LALU8CreateVector()}}
\index{\texttt{LALCHARCreateVector()}}
\index{\texttt{LALCreateVector()}}
\index{\texttt{LALZDestroyVector()}}
\index{\texttt{LALCDestroyVector()}}
\index{\texttt{LALDDestroyVector()}}
\index{\texttt{LALSDestroyVector()}}
\index{\texttt{LALI2DestroyVector()}}
\index{\texttt{LALI4DestroyVector()}}
\index{\texttt{LALI8DestroyVector()}}
\index{\texttt{LALU2DestroyVector()}}
\index{\texttt{LALU4DestroyVector()}}
\index{\texttt{LALU8DestroyVector()}}
\index{\texttt{LALCHARDestroyVector()}}
\index{\texttt{LALDestroyVector()}}

\subsubsection*{Description}

The \texttt{CreateVector} family of functions create a
$\langle\mbox{datatype}\rangle$\texttt{Vector} of the appropriate dimensions.

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
include(`DestroyVector.m4')

define(`TYPECODE',`C')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`D')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`S')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`I2')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`I4')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`I8')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`U2')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`U4')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`U8')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`CHAR')
include(`CreateVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`')
include(`CreateVector.m4')
include(`DestroyVector.m4')
