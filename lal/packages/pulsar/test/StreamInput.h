/********************************** <lalVerbatim file="StreamInputHV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\section{Header \texttt{StreamInput.h}}
\label{s:StreamInput.h}

Provides routines to read data from an open stream and store it in LAL
data structures.

\subsection*{Synopsis}
\begin{verbatim}
#include "StreamInput.h"
\end{verbatim}

\noindent This header provides prototypes for routines that construct
LAL data structures using the data from a file (or other I/O) stream.
The routines do not provide a system-level interface to create files
and open or close file streams; they simply assume that they have been
passed an open, readable stream.

According to the LAL specification in place at the time of writing
these modules, routines that perform file I/O (that is, which invoke
the C \verb@FILE@ type) should not be part of the LAL, but may be part
of standalone or test programs.  These routines, while quite generic,
do involve reading from file (or other input) streams, and so are not
included in the library; they are distributed as separate modules
within the \verb@test@ directory of the \verb@pulsar@ package.

However, these routines are compartmentalized in such a way that they
can easily be converted if the LAL specification is later expanded to
include I/O streams.  In partucular, the only file I/O commands used
are \verb@fgets()@ and \verb@feof()@.  Thus the upgrade would involve
only the following global changes:
\begin{enumerate}
\item Replace all occurrences of \verb@FILE *@ with the name of the
LAL I/O stream pointer type.
\item Replace all occurrences of \verb@fgets()@ and \verb@feof()@ with
equivalent LAL functions.
\end{enumerate}
In particular, there is no need to translate routines such as
\verb@fscanf()@; one should simply read data into a LAL
\verb@CHARVector@ and then use \verb@sscanf()@ to format the input.
This is the philosophy adopted in the following modules.

******************************************************* </lalLaTeX> */

#ifndef _STREAMINPUT_H
#define _STREAMINPUT_H

#include <lal/LALStdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

NRCSID(STREAMINPUTH,"$Id$");

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define STREAMINPUTH_ENUL 1
#define STREAMINPUTH_EOUT 2
#define STREAMINPUTH_EMEM 3
#define STREAMINPUTH_ELEN 4

#define STREAMINPUTH_MSGENUL "Unexpected null pointer in arguments"
#define STREAMINPUTH_MSGEOUT "Output handle points to a non-null pointer"
#define STREAMINPUTH_MSGEMEM "Memory allocation error"
#define STREAMINPUTH_MSGELEN "No numbers were read"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Structures}

******************************************************* </lalLaTeX> */

/* <lalLaTeX>
\vfill{\footnotesize\input{StreamInputHV}}
</lalLaTeX> */

/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{StreamVectorInputC}
</lalLaTeX> */
void
LALCHARReadVector( LALStatus  *stat, CHARVector **vector, FILE *stream );

void
LALI2ReadVector( LALStatus  *stat, INT2Vector **vector, FILE *stream );

void
LALI4ReadVector( LALStatus  *stat, INT4Vector **vector, FILE *stream );

void
LALI8ReadVector( LALStatus  *stat, INT8Vector **vector, FILE *stream );

void
LALU2ReadVector( LALStatus  *stat, UINT2Vector **vector, FILE *stream );

void
LALU4ReadVector( LALStatus  *stat, UINT4Vector **vector, FILE *stream );

void
LALU8ReadVector( LALStatus  *stat, UINT8Vector **vector, FILE *stream );

void
LALI2ReadVector( LALStatus  *stat, INT2Vector **vector, FILE *stream );

void
LALI4ReadVector( LALStatus  *stat, INT4Vector **vector, FILE *stream );

void
LALSReadVector( LALStatus  *stat, REAL4Vector **vector, FILE *stream );

void
LALDReadVector( LALStatus  *stat, REAL8Vector **vector, FILE *stream );

/* <lalLaTeX>
\newpage\input{StreamVectorSequenceInputC}
</lalLaTeX> */
void
LALCHARReadVectorSequence( LALStatus  *stat,
			   CHARVectorSequence **sequence,
			   FILE *stream );

void
LALI2ReadVectorSequence( LALStatus  *stat,
			 INT2VectorSequence **sequence,
			 FILE *stream );

void
LALI4ReadVectorSequence( LALStatus  *stat,
			 INT4VectorSequence **sequence,
			 FILE *stream );

void
LALI8ReadVectorSequence( LALStatus  *stat,
			 INT8VectorSequence **sequence,
			 FILE *stream );

void
LALU2ReadVectorSequence( LALStatus  *stat,
			 UINT2VectorSequence **sequence,
			 FILE *stream );

void
LALU4ReadVectorSequence( LALStatus  *stat,
			 UINT4VectorSequence **sequence,
			 FILE *stream );

void
LALU8ReadVectorSequence( LALStatus  *stat,
			 UINT8VectorSequence **sequence,
			 FILE *stream );

void
LALSReadVectorSequence( LALStatus  *stat,
			REAL4VectorSequence **sequence,
			FILE *stream );

void
LALDReadVectorSequence( LALStatus  *stat,
			REAL8VectorSequence **sequence,
			FILE *stream );

/* <lalLaTeX>
\newpage\input{StreamInputTestC}
</lalLaTeX> */

#ifdef __cplusplus
}
#endif

#endif /* _STREAMINPUT_H */
