/*----------------------------------------------------------------------- 

File Name: SeqFactories.h

<lalVerbatim file="SeqFactoriesHV">
Revision: $Id$
</lalVerbatim>

-------------------------------------------------------------------------*/

/* <lalLaTeX>

\section{Header \texttt{SeqFactories.h}}
\label{s:SeqFactories.h}

Provides prototype and status code information for use of CreateVectorSequence
and DestroyVectorSequence.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/SeqFactories.h>
\end{verbatim}

</lalLaTeX> */

#ifndef _SEQFACTORIES_H
#define _SEQFACTORIES_H

#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (SEQFACTORIESH, "$Id$");

/* <lalLaTeX>

\subsection*{Error conditions}
\input{SeqFactoriesHErrTab}

</lalLaTeX> */

/*
<lalErrTable file="SeqFactoriesHErrTab">
*/

#define SEQFACTORIESH_ESLENGTH  1
#define SEQFACTORIESH_EVLENGTH  2
#define SEQFACTORIESH_EALENGTH  4
#define SEQFACTORIESH_EVPTR     8
#define SEQFACTORIESH_EUPTR    16
#define SEQFACTORIESH_EDPTR    32
#define SEQFACTORIESH_EINPTR   64
#define SEQFACTORIESH_EMALLOC 128

#define SEQFACTORIESH_MSGESLENGTH "Illegal sequence length."
#define SEQFACTORIESH_MSGEVLENGTH "Illegal vector length."
#define SEQFACTORIESH_MSGEALENGTH "Illegal array dimension."
#define SEQFACTORIESH_MSGEVPTR    "Null sequence handle."
#define SEQFACTORIESH_MSGEUPTR    "Non-null sequence pointer."
#define SEQFACTORIESH_MSGEDPTR    "Null sequence data."
#define SEQFACTORIESH_MSGEINPTR   "Null input pointer."
#define SEQFACTORIESH_MSGEMALLOC  "Malloc failure."

/*
</lalErrTable>
*/


/* Structures. */
/* <lalLaTeX>

\subsection*{Structures}
\begin{verbatim}
CreateVectorSequenceIn
\end{verbatim}
\idx[Type]{CreateVectorSequenceIn}

\noindent This structure stores the input required for creating a vector
sequence.  This input includes the length of the sequence (i.e., the number of
vectors) and the length of each vector.  The fields are:

\begin{description}
\item[\texttt{UINT4 length}] The sequence length.
\item[\texttt{UINT4 vectorLength}] The length of each vector in the sequence.
\end{description}

</lalLaTeX> */

typedef struct tagCreateVectorSequenceIn {
  UINT4 length;
  UINT4 vectorLength;
} CreateVectorSequenceIn;


/* <lalLaTeX>

\begin{verbatim}
CreateArraySequenceIn
\end{verbatim}
\idx[Type]{CreateArraySequenceIn}

\noindent This structure stores the input required for creating an array
sequence.  This input includes the length of the sequence (i.e., the number of
array) and the dimensions of each array index.  The fields are:

\begin{description}
\item[\texttt{UINT4 length}] The sequence length.
\item[\texttt{UINT4Vector *dimLength}] The dimensions of each array
index (the same for every array in the sequence).
\end{description}

</lalLaTeX> */

typedef struct tagCreateArraySequenceIn {
  UINT4 length;
  UINT4Vector *dimLength;
} CreateArraySequenceIn;


/* Function prototypes. */
/* <lalLaTeX>
\newpage\input{VectorSequenceFactoriesC}
</lalLaTeX> */


void LALCreateSequence(LALStatus *, REAL4Sequence **, UINT4);
void LALDestroySequence(LALStatus *, REAL4Sequence **);

void LALCreateVectorSequence(LALStatus *, REAL4VectorSequence **,
                             CreateVectorSequenceIn *);
void LALDestroyVectorSequence(LALStatus *, REAL4VectorSequence**);

void LALCreateArraySequence(LALStatus *, REAL4ArraySequence **,
                            CreateArraySequenceIn *);
void LALDestroyArraySequence(LALStatus *, REAL4ArraySequence **);
                       
define(`TYPECODE',`CHAR')
include(`SeqFactoriesBaseH.m4')

define(`TYPECODE',`I2')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

define(`TYPECODE',`I4')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

define(`TYPECODE',`I8')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

define(`TYPECODE',`U2')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

define(`TYPECODE',`U4')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

define(`TYPECODE',`U8')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

define(`TYPECODE',`S')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

define(`TYPECODE',`D')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

define(`TYPECODE',`C')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')

define(`TYPECODE',`Z')
include(`SeqFactoriesBaseH.m4')
include(`ArraySeqFactoriesBaseH.m4')


/* <lalLaTeX>
\newpage\input{VectorSequenceFactoriesTestC}
</lalLaTeX> */

/* <lalLaTeX>
\newpage\input{ArraySequenceFactoriesTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _SEQFACTORIES_H */
