/*----------------------------------------------------------------------- 

File Name: AVFactories.h

<lalVerbatim file="AVFactoriesHV">
Revision: $Id$
</lalVerbatim>

-------------------------------------------------------------------------*/

/* <lalLaTeX>

\section{Header \texttt{AVFactories.h}}
\label{s:AVFactories.h}

Provides prototype and status code information for use of CreateVector,
CreateArray, LALDestroyVector and DestroyArray

\subsection*{Synopsis}
\begin{verbatim}
#include "AVFactories.h"
\end{verbatim}

</lalLaTeX> */

#ifndef _AVFACTORIES_H
#define _AVFACTORIES_H

#include "LALDatatypes.h"

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (AVFACTORIESH, "$Id$");

/* <lalLaTeX>

\subsection*{Error conditions}
\input{AVFactoriesHErrTab}

</lalLaTeX> */

/*
<lalErrTable file="AVFactoriesHErrTab">
*/
#define AVFACTORIESH_ELENGTH 1
#define AVFACTORIESH_EVPTR   2
#define AVFACTORIESH_EUPTR   4
#define AVFACTORIESH_EDPTR   8
#define AVFACTORIESH_EMALLOC 16
#define AVFACTORIESH_MSGELENGTH  "Illegal length."
#define AVFACTORIESH_MSGEVPTR    "Null vector/array handle."
#define AVFACTORIESH_MSGEUPTR    "Non-null vector/array pointer."
#define AVFACTORIESH_MSGEDPTR    "Null vector/array data."
#define AVFACTORIESH_MSGEMALLOC  "Malloc failure."
/*
</lalErrTable>
*/

/* Function prototypes. */
/* <lalLaTeX>
\newpage\input{VectorFactoriesC}
\newpage\input{ArrayFactoriesC}
</lalLaTeX> */

void LALCHARCreateVector(LALStatus *, CHARVector **, UINT4);
void LALI2CreateVector(LALStatus *, INT2Vector **, UINT4);
void LALI4CreateVector(LALStatus *, INT4Vector **, UINT4);
void LALI8CreateVector(LALStatus *, INT8Vector **, UINT4);
void LALU2CreateVector(LALStatus *, UINT2Vector **, UINT4);
void LALU4CreateVector(LALStatus *, UINT4Vector **, UINT4);
void LALU8CreateVector(LALStatus *, UINT8Vector **, UINT4);
void LALCreateVector(LALStatus *, REAL4Vector **, UINT4);
void LALSCreateVector(LALStatus *, REAL4Vector **, UINT4);
void LALDCreateVector(LALStatus *, REAL8Vector **, UINT4);
void LALCCreateVector(LALStatus *, COMPLEX8Vector **, UINT4);
void LALZCreateVector(LALStatus *, COMPLEX16Vector **, UINT4);

void LALI2CreateArray(LALStatus *, INT2Array **, UINT4Vector *);
void LALI4CreateArray(LALStatus *, INT4Array **, UINT4Vector *);
void LALI8CreateArray(LALStatus *, INT8Array **, UINT4Vector *);
void LALU2CreateArray(LALStatus *, UINT2Array **, UINT4Vector *);
void LALU4CreateArray(LALStatus *, UINT4Array **, UINT4Vector *);
void LALU8CreateArray(LALStatus *, UINT8Array **, UINT4Vector *);
void LALCreateArray(LALStatus *, REAL4Array **, UINT4Vector *);
void LALSCreateArray(LALStatus *, REAL4Array **, UINT4Vector *);
void LALDCreateArray(LALStatus *, REAL8Array **, UINT4Vector *);
void LALCCreateArray(LALStatus *, COMPLEX8Array **, UINT4Vector *);
void LALZCreateArray(LALStatus *, COMPLEX16Array **, UINT4Vector *);
 
void LALCHARDestroyVector(LALStatus *, CHARVector **);
void LALI2DestroyVector(LALStatus *, INT2Vector **);
void LALI4DestroyVector(LALStatus *, INT4Vector **);
void LALI8DestroyVector(LALStatus *, INT8Vector **);
void LALU2DestroyVector(LALStatus *, UINT2Vector **);
void LALU4DestroyVector(LALStatus *, UINT4Vector **);
void LALU8DestroyVector(LALStatus *, UINT8Vector **);
void LALDestroyVector(LALStatus *, REAL4Vector **);
void LALSDestroyVector(LALStatus *, REAL4Vector **);
void LALDDestroyVector(LALStatus *, REAL8Vector **);
void LALCDestroyVector(LALStatus *, COMPLEX8Vector **);
void LALZDestroyVector(LALStatus *, COMPLEX16Vector **);

void LALI2DestroyArray(LALStatus *, INT2Array **);
void LALI4DestroyArray(LALStatus *, INT4Array **);
void LALI8DestroyArray(LALStatus *, INT8Array **);
void LALU2DestroyArray(LALStatus *, UINT2Array **);
void LALU4DestroyArray(LALStatus *, UINT4Array **);
void LALU8DestroyArray(LALStatus *, UINT8Array **);
void LALDestroyArray(LALStatus *, REAL4Array **);
void LALSDestroyArray(LALStatus *, REAL4Array **);
void LALDDestroyArray(LALStatus *, REAL8Array **);
void LALCDestroyArray(LALStatus *, COMPLEX8Array **);
void LALZDestroyArray(LALStatus *, COMPLEX16Array **);

/* Test program. */

/* <lalLaTeX>
\newpage\input{VectorFactoriesTestC}
\newpage\input{ArrayFactoriesTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _AVFACTORIES_H */
