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
passed an open, readable stream.  Nonetheless, because they involve
I/O stream manipulation, these routines are placed in the
\verb@lalsupport@ library rather than in \verb@lal@ proper.

The routines in \verb@StreamVectorInput.c@ and
\verb@StreamVectorSequenceInput.c@ are compartmentalized in such a way
that they can easily be converted if the LAL specification later
changes the way in which I/O streams are handled.  In partucular, the
only file I/O commands used are \verb@fgets()@ and \verb@feof()@.
Thus the upgrade would involve only the following global changes:
\begin{enumerate}
\item Replace all occurrences of \verb@FILE *@ with the name of the
LAL I/O stream pointer type.
\item Replace all occurrences of \verb@fgets()@ and \verb@feof()@ with
equivalent LAL functions.
\end{enumerate}
In particular, there is no need to translate routines such as
\verb@fscanf()@; one should simply read data into a LAL
\verb@CHARVector@ and then use \verb@sscanf()@ to format the input.
This is the approach used in the numerical input routines in
\verb@StreamVectorInput.c@ and \verb@StreamVectorSequenceInput.c@.

The routines in \verb@StreamSequenceInput.c@ are less robust but much
more efficient: they use \verb@fscanf()@ to parse the input stream
directly.  They are intended primarily for test programs that may need
to read large datafiles of undetermined length.


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
#define STREAMINPUTH_ENUL  1
#define STREAMINPUTH_EOUT  2
#define STREAMINPUTH_EMEM  3
#define STREAMINPUTH_ELEN  4
#define STREAMINPUTH_ESLEN 5
#define STREAMINPUTH_EVLEN 6
#define STREAMINPUTH_EDLEN 7
#define STREAMINPUTH_EDIM  8
#define STREAMINPUTH_EFMT  9

#define STREAMINPUTH_MSGENUL  "Unexpected null pointer in arguments"
#define STREAMINPUTH_MSGEOUT  "Output handle points to a non-null pointer"
#define STREAMINPUTH_MSGEMEM  "Memory allocation error"
#define STREAMINPUTH_MSGELEN  "No numbers were read"
#define STREAMINPUTH_MSGESLEN "Not enough numbers read to fill sequence"
#define STREAMINPUTH_MSGEVLEN "Could not determine complex vectorLength"
#define STREAMINPUTH_MSGEDLEN "No dimLength given"
#define STREAMINPUTH_MSGEDIM  "Inconsistent or non-positive arrayDim value"
#define STREAMINPUTH_MSGEFMT  "Badly formatted number"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Types}

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
LALI2ReadVector( LALStatus  *stat, INT2Vector **vector, FILE *stream, BOOLEAN strict );

void
LALI4ReadVector( LALStatus  *stat, INT4Vector **vector, FILE *stream, BOOLEAN strict );

void
LALI8ReadVector( LALStatus  *stat, INT8Vector **vector, FILE *stream, BOOLEAN strict );

void
LALU2ReadVector( LALStatus  *stat, UINT2Vector **vector, FILE *stream, BOOLEAN strict );

void
LALU4ReadVector( LALStatus  *stat, UINT4Vector **vector, FILE *stream, BOOLEAN strict );

void
LALU8ReadVector( LALStatus  *stat, UINT8Vector **vector, FILE *stream, BOOLEAN strict );

void
LALSReadVector( LALStatus  *stat, REAL4Vector **vector, FILE *stream, BOOLEAN strict );

void
LALDReadVector( LALStatus  *stat, REAL8Vector **vector, FILE *stream, BOOLEAN strict );

/* <lalLaTeX>
\newpage\input{StreamVectorSequenceInputC}
</lalLaTeX> */
void
LALCHARReadVectorSequence( LALStatus  *stat, CHARVectorSequence **sequence, FILE *stream );

void
LALI2ReadVectorSequence( LALStatus  *stat, INT2VectorSequence **sequence, FILE *stream );

void
LALI4ReadVectorSequence( LALStatus  *stat, INT4VectorSequence **sequence, FILE *stream );

void
LALI8ReadVectorSequence( LALStatus  *stat, INT8VectorSequence **sequence, FILE *stream );

void
LALU2ReadVectorSequence( LALStatus  *stat, UINT2VectorSequence **sequence, FILE *stream );

void
LALU4ReadVectorSequence( LALStatus  *stat, UINT4VectorSequence **sequence, FILE *stream );

void
LALU8ReadVectorSequence( LALStatus  *stat, UINT8VectorSequence **sequence, FILE *stream );

void
LALSReadVectorSequence( LALStatus  *stat, REAL4VectorSequence **sequence, FILE *stream );

void
LALDReadVectorSequence( LALStatus  *stat, REAL8VectorSequence **sequence, FILE *stream );

/* <lalLaTeX>
\newpage\input{StreamSequenceInputC}
</lalLaTeX> */
void
LALCHARReadSequence( LALStatus *stat, CHARSequence **sequence, FILE *stream );

void
LALI2ReadSequence( LALStatus *stat, INT2Sequence **sequence, FILE *stream );

void
LALI4ReadSequence( LALStatus *stat, INT4Sequence **sequence, FILE *stream );

void
LALI8ReadSequence( LALStatus *stat, INT8Sequence **sequence, FILE *stream );

void
LALU2ReadSequence( LALStatus *stat, UINT2Sequence **sequence, FILE *stream );

void
LALU4ReadSequence( LALStatus *stat, UINT4Sequence **sequence, FILE *stream );

void
LALU8ReadSequence( LALStatus *stat, UINT8Sequence **sequence, FILE *stream );

void
LALSReadSequence( LALStatus *stat, REAL4Sequence **sequence, FILE *stream );

void
LALDReadSequence( LALStatus *stat, REAL8Sequence **sequence, FILE *stream );

void
LALCReadSequence( LALStatus *stat, COMPLEX8Sequence **sequence, FILE *stream );

void
LALZReadSequence( LALStatus *stat, COMPLEX16Sequence **sequence, FILE *stream );

/* <lalLaTeX>
\newpage\input{StreamSeriesInputC}
</lalLaTeX> */
void
LALI2ReadTSeries( LALStatus *stat, INT2TimeSeries *series, FILE *stream );
void
LALI4ReadTSeries( LALStatus *stat, INT4TimeSeries *series, FILE *stream );
void
LALI8ReadTSeries( LALStatus *stat, INT8TimeSeries *series, FILE *stream );
void
LALU2ReadTSeries( LALStatus *stat, UINT2TimeSeries *series, FILE *stream );
void
LALU4ReadTSeries( LALStatus *stat, UINT4TimeSeries *series, FILE *stream );
void
LALU8ReadTSeries( LALStatus *stat, UINT8TimeSeries *series, FILE *stream );
void
LALSReadTSeries( LALStatus *stat, REAL4TimeSeries *series, FILE *stream );
void
LALDReadTSeries( LALStatus *stat, REAL8TimeSeries *series, FILE *stream );
void
LALCReadTSeries( LALStatus *stat, COMPLEX8TimeSeries *series, FILE *stream );
void
LALZReadTSeries( LALStatus *stat, COMPLEX16TimeSeries *series, FILE *stream );

void
LALI2ReadTVectorSeries( LALStatus *stat, INT2TimeVectorSeries *series, FILE *stream );
void
LALI4ReadTVectorSeries( LALStatus *stat, INT4TimeVectorSeries *series, FILE *stream );
void
LALI8ReadTVectorSeries( LALStatus *stat, INT8TimeVectorSeries *series, FILE *stream );
void
LALU2ReadTVectorSeries( LALStatus *stat, UINT2TimeVectorSeries *series, FILE *stream );
void
LALU4ReadTVectorSeries( LALStatus *stat, UINT4TimeVectorSeries *series, FILE *stream );
void
LALU8ReadTVectorSeries( LALStatus *stat, UINT8TimeVectorSeries *series, FILE *stream );
void
LALSReadTVectorSeries( LALStatus *stat, REAL4TimeVectorSeries *series, FILE *stream );
void
LALDReadTVectorSeries( LALStatus *stat, REAL8TimeVectorSeries *series, FILE *stream );
void
LALCReadTVectorSeries( LALStatus *stat, COMPLEX8TimeVectorSeries *series, FILE *stream );
void
LALZReadTVectorSeries( LALStatus *stat, COMPLEX16TimeVectorSeries *series, FILE *stream );

void
LALI2ReadTArraySeries( LALStatus *stat, INT2TimeArraySeries *series, FILE *stream );
void
LALI4ReadTArraySeries( LALStatus *stat, INT4TimeArraySeries *series, FILE *stream );
void
LALI8ReadTArraySeries( LALStatus *stat, INT8TimeArraySeries *series, FILE *stream );
void
LALU2ReadTArraySeries( LALStatus *stat, UINT2TimeArraySeries *series, FILE *stream );
void
LALU4ReadTArraySeries( LALStatus *stat, UINT4TimeArraySeries *series, FILE *stream );
void
LALU8ReadTArraySeries( LALStatus *stat, UINT8TimeArraySeries *series, FILE *stream );
void
LALSReadTArraySeries( LALStatus *stat, REAL4TimeArraySeries *series, FILE *stream );
void
LALDReadTArraySeries( LALStatus *stat, REAL8TimeArraySeries *series, FILE *stream );
void
LALCReadTArraySeries( LALStatus *stat, COMPLEX8TimeArraySeries *series, FILE *stream );
void
LALZReadTArraySeries( LALStatus *stat, COMPLEX16TimeArraySeries *series, FILE *stream );

void
LALI2ReadFSeries( LALStatus *stat, INT2FrequencySeries *series, FILE *stream );
void
LALI4ReadFSeries( LALStatus *stat, INT4FrequencySeries *series, FILE *stream );
void
LALI8ReadFSeries( LALStatus *stat, INT8FrequencySeries *series, FILE *stream );
void
LALU2ReadFSeries( LALStatus *stat, UINT2FrequencySeries *series, FILE *stream );
void
LALU4ReadFSeries( LALStatus *stat, UINT4FrequencySeries *series, FILE *stream );
void
LALU8ReadFSeries( LALStatus *stat, UINT8FrequencySeries *series, FILE *stream );
void
LALSReadFSeries( LALStatus *stat, REAL4FrequencySeries *series, FILE *stream );
void
LALDReadFSeries( LALStatus *stat, REAL8FrequencySeries *series, FILE *stream );
void
LALCReadFSeries( LALStatus *stat, COMPLEX8FrequencySeries *series, FILE *stream );
void
LALZReadFSeries( LALStatus *stat, COMPLEX16FrequencySeries *series, FILE *stream );

/* <lalLaTeX>
%\newpage\input{StreamInputTestC}
</lalLaTeX> */

#ifdef __cplusplus
}
#endif

#endif /* _STREAMINPUT_H */
