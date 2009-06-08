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
\idx{LALZCreateVectorSequence()}
\idx{LALCCreateVectorSequence()}
\idx{LALDCreateVectorSequence()}
\idx{LALSCreateVectorSequence()}
\idx{LALI2CreateVectorSequence()}
\idx{LALI4CreateVectorSequence()}
\idx{LALI8CreateVectorSequence()}
\idx{LALU2CreateVectorSequence()}
\idx{LALU4CreateVectorSequence()}
\idx{LALU8CreateVectorSequence()}
\idx{LALCHARCreateVectorSequence()}
\idx{LALCreateVectorSequence()}
\idx{LALZDestroyVectorSequence()}
\idx{LALCDestroyVectorSequence()}
\idx{LALDDestroyVectorSequence()}
\idx{LALSDestroyVectorSequence()}
\idx{LALI2DestroyVectorSequence()}
\idx{LALI4DestroyVectorSequence()}
\idx{LALI8DestroyVectorSequence()}
\idx{LALU2DestroyVectorSequence()}
\idx{LALU4DestroyVectorSequence()}
\idx{LALU8DestroyVectorSequence()}
\idx{LALCHARDestroyVectorSequence()}
\idx{LALDestroyVectorSequence()}

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
