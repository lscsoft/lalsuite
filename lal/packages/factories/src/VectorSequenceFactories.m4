/*----------------------------------------------------------------------- 

File Name: VectorSequenceFactories.c

<lalVerbatim file="VectorSequenceFactoriesCV">
Revision: $Id$
</lalVerbatim>

-------------------------------------------------------------------------*/

/* <lalLaTeX>

\subsection{Module \texttt{VectorSequenceFactories.c}}
\label{ss:VectorSequenceFactories.c}

Create/destroy $\langle\mbox{datatype}\rangle$VectorSequence objects. 

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{VectorSequenceFactoriesD}
\index{\texttt{LALZCreateVectorSequence()}}
\index{\texttt{LALCCreateVectorSequence()}}
\index{\texttt{LALDCreateVectorSequence()}}
\index{\texttt{LALSCreateVectorSequence()}}
\index{\texttt{LALI2CreateVectorSequence()}}
\index{\texttt{LALI4CreateVectorSequence()}}
\index{\texttt{LALI8CreateVectorSequence()}}
\index{\texttt{LALU2CreateVectorSequence()}}
\index{\texttt{LALU4CreateVectorSequence()}}
\index{\texttt{LALU8CreateVectorSequence()}}
\index{\texttt{LALCHARCreateVectorSequence()}}
\index{\texttt{LALCreateVectorSequence()}}
\index{\texttt{LALZDestroyVectorSequence()}}
\index{\texttt{LALCDestroyVectorSequence()}}
\index{\texttt{LALDDestroyVectorSequence()}}
\index{\texttt{LALSDestroyVectorSequence()}}
\index{\texttt{LALI2DestroyVectorSequence()}}
\index{\texttt{LALI4DestroyVectorSequence()}}
\index{\texttt{LALI8DestroyVectorSequence()}}
\index{\texttt{LALU2DestroyVectorSequence()}}
\index{\texttt{LALU4DestroyVectorSequence()}}
\index{\texttt{LALU8DestroyVectorSequence()}}
\index{\texttt{LALCHARDestroyVectorSequence()}}
\index{\texttt{LALDestroyVectorSequence()}}

\subsubsection*{Description}

The \texttt{CreateVectorSequence} family of functions create a
$\langle\mbox{datatype}\rangle$\texttt{VectorSequence} of the
appropriate dimensions.

The \texttt{DestroyVectorSequence} family of functions return the storage
allocated by the \texttt{CreateVectorSequence} functions to the system.

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()
LALFree()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{VectorSequenceFactoriesCV}}

</lalLaTeX> */

#include "LALStdlib.h"
#include "SeqFactories.h"

/* <lalVerbatim file="VectorSequenceFactoriesNRCSID"> */
NRCSID( VECTORSEQUENCEFACTORIESC, "$Id$" );
/* </lalVerbatim> */

define(`TYPECODE',`Z')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`C')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`D')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`S')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`I2')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`I4')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`I8')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`U2')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`U4')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`U8')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`CHAR')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')
