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
#include "SeqFactories.h"
\end{verbatim}

</lalLaTeX> */

#ifndef _SEQFACTORIES_H
#define _SEQFACTORIES_H

#include "LALDatatypes.h"
#include "AVFactories.h"

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
\index{\verb&CreateVectorSequenceIn&}

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

void CreateSequence(Status *, REAL4Sequence **, UINT4);
void CHARCreateSequence(Status *, CHARSequence **, UINT4);
void I2CreateSequence(Status *, INT2Sequence **, UINT4);
void I4CreateSequence(Status *, INT4Sequence **, UINT4);
void I8CreateSequence(Status *, INT8Sequence **, UINT4);
void U2CreateSequence(Status *, UINT2Sequence **, UINT4);
void U4CreateSequence(Status *, UINT4Sequence **, UINT4);
void U8CreateSequence(Status *, UINT8Sequence **, UINT4);
void SCreateSequence(Status *, REAL4Sequence **, UINT4);
void DCreateSequence(Status *, REAL8Sequence **, UINT4);
void CCreateSequence(Status *, COMPLEX8Sequence **, UINT4);
void ZCreateSequence(Status *, COMPLEX16Sequence **, UINT4);

void DestroySequence(Status *, REAL4Sequence **);
void CHARDestroySequence(Status *, CHARSequence **);
void I2DestroySequence(Status *, INT2Sequence **);
void I4DestroySequence(Status *, INT4Sequence **);
void I8DestroySequence(Status *, INT8Sequence **);
void U2DestroySequence(Status *, UINT2Sequence **);
void U4DestroySequence(Status *, UINT4Sequence **);
void U8DestroySequence(Status *, UINT8Sequence **);
void SDestroySequence(Status *, REAL4Sequence **);
void DDestroySequence(Status *, REAL8Sequence **);
void CDestroySequence(Status *, COMPLEX8Sequence **);
void ZDestroySequence(Status *, COMPLEX16Sequence **);

void CreateVectorSequence(Status *, 
                             REAL4VectorSequence **,
			     CreateVectorSequenceIn *);
void CHARCreateVectorSequence(Status *, 
                             CHARVectorSequence **,
			     CreateVectorSequenceIn *);
void I2CreateVectorSequence(Status *, 
			     INT2VectorSequence **,
			     CreateVectorSequenceIn *);
void I4CreateVectorSequence(Status *, 
			     INT4VectorSequence **,
			     CreateVectorSequenceIn *);
void I8CreateVectorSequence(Status *, 
			     INT8VectorSequence **,
			     CreateVectorSequenceIn *);
void U2CreateVectorSequence(Status *, 
			     UINT2VectorSequence **,
			     CreateVectorSequenceIn *);
void U4CreateVectorSequence(Status *, 
			     UINT4VectorSequence **,
			     CreateVectorSequenceIn *);
void U8CreateVectorSequence(Status *, 
			     UINT8VectorSequence **,
			     CreateVectorSequenceIn *);
void SCreateVectorSequence(Status *, 
			     REAL4VectorSequence **,
			     CreateVectorSequenceIn *);
void DCreateVectorSequence(Status *, 
			     REAL8VectorSequence **,
			     CreateVectorSequenceIn *);
void CCreateVectorSequence(Status *, 
			     COMPLEX8VectorSequence **, 
			     CreateVectorSequenceIn *);
void ZCreateVectorSequence(Status *, 
			     COMPLEX16VectorSequence **, 
			     CreateVectorSequenceIn *);

void DestroyVectorSequence (Status *, 
                             REAL4VectorSequence **);
void CHARDestroyVectorSequence (Status *, 
                             CHARVectorSequence **);
void I2DestroyVectorSequence(Status *, 
			     INT2VectorSequence **);
void I4DestroyVectorSequence(Status *, 
			     INT4VectorSequence **);
void I8DestroyVectorSequence(Status *, 
			     INT8VectorSequence **);
void U2DestroyVectorSequence(Status *, 
			     UINT2VectorSequence **);
void U4DestroyVectorSequence(Status *, 
			     UINT4VectorSequence **);
void U8DestroyVectorSequence(Status *, 
			     UINT8VectorSequence **);
void SDestroyVectorSequence(Status *, 
			     REAL4VectorSequence **);
void DDestroyVectorSequence(Status *, 
			     REAL8VectorSequence **);
void CDestroyVectorSequence(Status *, 
			     COMPLEX8VectorSequence **);
void ZDestroyVectorSequence(Status *, 
			     COMPLEX16VectorSequence **);

/* Test program. */

/* <lalLaTeX>
\newpage\input{VectorSequenceFactoriesTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _SEQFACTORIES_H */
