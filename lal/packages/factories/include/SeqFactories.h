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
#define SEQFACTORIESH_EVPTR     4
#define SEQFACTORIESH_EUPTR     8
#define SEQFACTORIESH_EDPTR    16
#define SEQFACTORIESH_EINPTR   32
#define SEQFACTORIESH_EMALLOC  64

#define SEQFACTORIESH_MSGESLENGTH "Illegal sequence length."
#define SEQFACTORIESH_MSGEVLENGTH "Illegal vector length."
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
\index{\texttt{CreateVectorSequenceIn}}

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



/* Function prototypes. */
/* <lalLaTeX>
\newpage\input{VectorSequenceFactoriesC}
</lalLaTeX> */

void LALCreateSequence(LALStatus *, REAL4Sequence **, UINT4);
void LALCHARCreateSequence(LALStatus *, CHARSequence **, UINT4);
void LALI2CreateSequence(LALStatus *, INT2Sequence **, UINT4);
void LALI4CreateSequence(LALStatus *, INT4Sequence **, UINT4);
void LALI8CreateSequence(LALStatus *, INT8Sequence **, UINT4);
void LALU2CreateSequence(LALStatus *, UINT2Sequence **, UINT4);
void LALU4CreateSequence(LALStatus *, UINT4Sequence **, UINT4);
void LALU8CreateSequence(LALStatus *, UINT8Sequence **, UINT4);
void LALSCreateSequence(LALStatus *, REAL4Sequence **, UINT4);
void LALDCreateSequence(LALStatus *, REAL8Sequence **, UINT4);
void LALCCreateSequence(LALStatus *, COMPLEX8Sequence **, UINT4);
void LALZCreateSequence(LALStatus *, COMPLEX16Sequence **, UINT4);

void LALDestroySequence(LALStatus *, REAL4Sequence **);
void LALCHARDestroySequence(LALStatus *, CHARSequence **);
void LALI2DestroySequence(LALStatus *, INT2Sequence **);
void LALI4DestroySequence(LALStatus *, INT4Sequence **);
void LALI8DestroySequence(LALStatus *, INT8Sequence **);
void LALU2DestroySequence(LALStatus *, UINT2Sequence **);
void LALU4DestroySequence(LALStatus *, UINT4Sequence **);
void LALU8DestroySequence(LALStatus *, UINT8Sequence **);
void LALSDestroySequence(LALStatus *, REAL4Sequence **);
void LALDDestroySequence(LALStatus *, REAL8Sequence **);
void LALCDestroySequence(LALStatus *, COMPLEX8Sequence **);
void LALZDestroySequence(LALStatus *, COMPLEX16Sequence **);

void LALCreateVectorSequence(LALStatus *, 
                             REAL4VectorSequence **,
			     CreateVectorSequenceIn *);
void LALCHARCreateVectorSequence(LALStatus *, 
                             CHARVectorSequence **,
			     CreateVectorSequenceIn *);
void LALI2CreateVectorSequence(LALStatus *, 
			     INT2VectorSequence **,
			     CreateVectorSequenceIn *);
void LALI4CreateVectorSequence(LALStatus *, 
			     INT4VectorSequence **,
			     CreateVectorSequenceIn *);
void LALI8CreateVectorSequence(LALStatus *, 
			     INT8VectorSequence **,
			     CreateVectorSequenceIn *);
void LALU2CreateVectorSequence(LALStatus *, 
			     UINT2VectorSequence **,
			     CreateVectorSequenceIn *);
void LALU4CreateVectorSequence(LALStatus *, 
			     UINT4VectorSequence **,
			     CreateVectorSequenceIn *);
void LALU8CreateVectorSequence(LALStatus *, 
			     UINT8VectorSequence **,
			     CreateVectorSequenceIn *);
void LALSCreateVectorSequence(LALStatus *, 
			     REAL4VectorSequence **,
			     CreateVectorSequenceIn *);
void LALDCreateVectorSequence(LALStatus *, 
			     REAL8VectorSequence **,
			     CreateVectorSequenceIn *);
void LALCCreateVectorSequence(LALStatus *, 
			     COMPLEX8VectorSequence **, 
			     CreateVectorSequenceIn *);
void LALZCreateVectorSequence(LALStatus *, 
			     COMPLEX16VectorSequence **, 
			     CreateVectorSequenceIn *);

void LALDestroyVectorSequence (LALStatus *, 
                             REAL4VectorSequence **);
void LALCHARDestroyVectorSequence (LALStatus *, 
                             CHARVectorSequence **);
void LALI2DestroyVectorSequence(LALStatus *, 
			     INT2VectorSequence **);
void LALI4DestroyVectorSequence(LALStatus *, 
			     INT4VectorSequence **);
void LALI8DestroyVectorSequence(LALStatus *, 
			     INT8VectorSequence **);
void LALU2DestroyVectorSequence(LALStatus *, 
			     UINT2VectorSequence **);
void LALU4DestroyVectorSequence(LALStatus *, 
			     UINT4VectorSequence **);
void LALU8DestroyVectorSequence(LALStatus *, 
			     UINT8VectorSequence **);
void LALSDestroyVectorSequence(LALStatus *, 
			     REAL4VectorSequence **);
void LALDDestroyVectorSequence(LALStatus *, 
			     REAL8VectorSequence **);
void LALCDestroyVectorSequence(LALStatus *, 
			     COMPLEX8VectorSequence **);
void LALZDestroyVectorSequence(LALStatus *, 
			     COMPLEX16VectorSequence **);

/* Test program. */

/* <lalLaTeX>
\newpage\input{VectorSequenceFactoriesTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _SEQFACTORIES_H */
