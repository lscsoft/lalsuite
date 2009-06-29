/*-----------------------------------------------------------------------

File Name: ArraySequenceFactories.c

<lalVerbatim file="ArraySequenceFactoriesCV">
Revision: $Id$
</lalVerbatim>

-------------------------------------------------------------------------*/

/* <lalLaTeX>

\subsection{Module \texttt{ArraySequenceFactories.c}}
\label{ss:ArraySequenceFactories.c}

Create/destroy $\langle\mbox{datatype}\rangle$ArraySequence objects.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{ArraySequenceFactoriesD}
\idx{LALZCreateArraySequence()}
\idx{LALCCreateArraySequence()}
\idx{LALDCreateArraySequence()}
\idx{LALSCreateArraySequence()}
\idx{LALI2CreateArraySequence()}
\idx{LALI4CreateArraySequence()}
\idx{LALI8CreateArraySequence()}
\idx{LALU2CreateArraySequence()}
\idx{LALU4CreateArraySequence()}
\idx{LALU8CreateArraySequence()}
\idx{LALCreateArraySequence()}
\idx{LALZDestroyArraySequence()}
\idx{LALCDestroyArraySequence()}
\idx{LALDDestroyArraySequence()}
\idx{LALSDestroyArraySequence()}
\idx{LALI2DestroyArraySequence()}
\idx{LALI4DestroyArraySequence()}
\idx{LALI8DestroyArraySequence()}
\idx{LALU2DestroyArraySequence()}
\idx{LALU4DestroyArraySequence()}
\idx{LALU8DestroyArraySequence()}
\idx{LALDestroyArraySequence()}

\subsubsection*{Description}

The \texttt{CreateArraySequence} family of functions create a
$\langle\mbox{datatype}\rangle$\texttt{ArraySequence} of the
appropriate dimensions.

The \texttt{DestroyArraySequence} family of functions return the storage
allocated by the \texttt{CreateArraySequence} functions to the system.

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()
LALFree()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{ArraySequenceFactoriesCV}}

</lalLaTeX> */

#include "LALStdlib.h"
#include "AVFactories.h"
#include "SeqFactories.h"

/* <lalVerbatim file="ArraySequenceFactoriesNRCSID"> */
NRCSID( ARRAYSEQUENCEFACTORIESC, "$Id$" );
/* </lalVerbatim> */

define(`TYPECODE',`Z')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`C')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`D')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`S')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`I2')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`I4')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`I8')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`U2')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`U4')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`U8')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')
