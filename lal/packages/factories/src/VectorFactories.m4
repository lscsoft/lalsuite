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
\index{\verb&ZCreateVector()&}
\index{\verb&CCreateVector()&}
\index{\verb&DCreateVector()&}
\index{\verb&SCreateVector()&}
\index{\verb&I2CreateVector()&}
\index{\verb&I4CreateVector()&}
\index{\verb&I8CreateVector()&}
\index{\verb&U2CreateVector()&}
\index{\verb&U4CreateVector()&}
\index{\verb&U8CreateVector()&}
\index{\verb&CHARCreateVector()&}
\index{\verb&CreateVector()&}
\index{\verb&ZDestroyVector()&}
\index{\verb&CDestroyVector()&}
\index{\verb&DDestroyVector()&}
\index{\verb&SDestroyVector()&}
\index{\verb&I2DestroyVector()&}
\index{\verb&I4DestroyVector()&}
\index{\verb&I8DestroyVector()&}
\index{\verb&U2DestroyVector()&}
\index{\verb&U4DestroyVector()&}
\index{\verb&U8DestroyVector()&}
\index{\verb&CHARDestroyVector()&}
\index{\verb&DestroyVector()&}

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
