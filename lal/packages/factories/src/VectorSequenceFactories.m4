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
\index{\verb&ZCreateVectorSequence()&}
\index{\verb&CCreateVectorSequence()&}
\index{\verb&DCreateVectorSequence()&}
\index{\verb&SCreateVectorSequence()&}
\index{\verb&I2CreateVectorSequence()&}
\index{\verb&I4CreateVectorSequence()&}
\index{\verb&I8CreateVectorSequence()&}
\index{\verb&U2CreateVectorSequence()&}
\index{\verb&U4CreateVectorSequence()&}
\index{\verb&U8CreateVectorSequence()&}
\index{\verb&CHARCreateVectorSequence()&}
\index{\verb&CreateVectorSequence()&}
\index{\verb&ZDestroyVectorSequence()&}
\index{\verb&CDestroyVectorSequence()&}
\index{\verb&DDestroyVectorSequence()&}
\index{\verb&SDestroyVectorSequence()&}
\index{\verb&I2DestroyVectorSequence()&}
\index{\verb&I4DestroyVectorSequence()&}
\index{\verb&I8DestroyVectorSequence()&}
\index{\verb&U2DestroyVectorSequence()&}
\index{\verb&U4DestroyVectorSequence()&}
\index{\verb&U8DestroyVectorSequence()&}
\index{\verb&CHARDestroyVectorSequence()&}
\index{\verb&DestroyVectorSequence()&}

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
