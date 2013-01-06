/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#ifndef _STREAMINPUT_H
#define _STREAMINPUT_H

#include <lal/LALStdlib.h>
#include <lal/Grid.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
   \addtogroup StreamInput_h
\author Creighton, T. D.

   \brief Provides routines to read data from an open stream and store it in LAL data structures.

\heading{Synopsis}
\code
#include "StreamInput.h"
\endcode

This header provides prototypes for routines that construct
LAL data structures using the data from a file (or other I/O) stream.
The routines do not provide a system-level interface to create files
and open or close file streams; they simply assume that they have been
passed an open, readable stream.  Nonetheless, because they involve
I/O stream manipulation, these routines are placed in the
\c lalsupport library rather than in \c lal proper.

The routines in \ref StreamVectorInput_c and
\ref StreamVectorSequenceInput_c are compartmentalized in such a way
that they can easily be converted if the LAL specification later
changes the way in which I/O streams are handled.  In partucular, the
only file I/O commands used are <tt>fgets()</tt> and <tt>feof()</tt>.
Thus the upgrade would involve only the following global changes:
<ol>
<li> Replace all occurrences of <tt>FILE *</tt> with the name of the
LAL I/O stream pointer type.</li>
<li> Replace all occurrences of <tt>fgets()</tt> and <tt>feof()</tt> with
equivalent LAL functions.</li>
</ol>
In particular, there is no need to translate routines such as
<tt>fscanf()</tt>; one should simply read data into a LAL
\c CHARVector and then use <tt>sscanf()</tt> to format the input.
This is the approach used in the numerical input routines in
\ref StreamVectorInput_c and \ref StreamVectorSequenceInput_c.

The routines in \ref StreamSequenceInput_c are less robust but much
more efficient: they use <tt>fscanf()</tt> to parse the input stream
directly.  They are intended primarily for test programs that may need
to read large datafiles of undetermined length.  The routines in
\ref StreamSeriesInput_c and \ref StreamGridInput_c also parse the
input stream directly using <tt>fscanf()</tt>, to avoid potentially
crippling computational overhead.

*/ /*@{*/

/** \name Error Codes */ /*@{*/
#define STREAMINPUTH_ENUL  1    /**< Unexpected null pointer in arguments */
#define STREAMINPUTH_EOUT  2    /**< Output handle points to a non-null pointer */
#define STREAMINPUTH_EMEM  3    /**< Memory allocation error */
#define STREAMINPUTH_ELEN  4    /**< No numbers were read */
#define STREAMINPUTH_ESLEN 5    /**< Not enough numbers read to fill sequence */
#define STREAMINPUTH_EVLEN 6    /**< Could not determine complex vectorLength */
#define STREAMINPUTH_EDLEN 7    /**< Dimension lengths inconsistent or not given */
#define STREAMINPUTH_EDIM  8    /**< Inconsistent or non-positive arrayDim value */
#define STREAMINPUTH_EFMT  9    /**< Badly formatted number */
#define STREAMINPUTH_EBUF 10    /**< BUFFSIZE not a multiple of largest complex type size */
/*@}*/
/*@}*/

#define STREAMINPUTH_MSGENUL  "Unexpected null pointer in arguments"
#define STREAMINPUTH_MSGEOUT  "Output handle points to a non-null pointer"
#define STREAMINPUTH_MSGEMEM  "Memory allocation error"
#define STREAMINPUTH_MSGELEN  "No numbers were read"
#define STREAMINPUTH_MSGESLEN "Not enough numbers read to fill sequence"
#define STREAMINPUTH_MSGEVLEN "Could not determine complex vectorLength"
#define STREAMINPUTH_MSGEDLEN "Dimension lengths inconsistent or not given"
#define STREAMINPUTH_MSGEDIM  "Inconsistent or non-positive arrayDim value"
#define STREAMINPUTH_MSGEFMT  "Badly formatted number"
#define STREAMINPUTH_MSGEBUF  "BUFFSIZE not a multiple of largest complex type size"


/* Function prototypes. */


/**
   \defgroup StreamVectorInput_c Module StreamVectorInput.c
   \ingroup StreamInput_h
   \author Creighton, T. D.

   \brief Reads data from a single line in an input stream.

\heading{Description}

These routines read ASCII data from the I/O stream <tt>*stream</tt>
until a newline or the end-of-input is reached.  (The line can be of
arbitrary length; the data is temporarily stored in a linked list of
buffers.)  Once read, a LAL vector structure <tt>**vector</tt> is
created and the data stored in it.  The routine passes back a pointer
to the new structure.  For the numerical routines, the \c strict
parameter determines whether the routine will do strict error checking
based on the contents of the input stream (see below).

The basic routine in this module is <tt>LALCHARReadVector()</tt>, which
simply stores bytes read from <tt>*stream</tt> until the next newline
character <tt>'\n'</tt>, null character <tt>'\0'</tt>, or the end of the
input as determined by the <tt>feof()</tt> function.  The vector
includes the newline (if present), and also an explicit <tt>'\0'</tt> at
the end, if one was not already present.  This routine should
\e not be used to read a binary data stream, which are not
logically divided into ``lines''.  Unless it aborts due to invalid
arguments or failed memory allocation, <tt>LALCHARReadVector()</tt> will
always return successfully regardless of the contents of the input
stream; <tt>*vector</tt> will created containing at least a single
<tt>'\0'</tt> terminator, if nothing else.

The other routines in this module use <tt>LALCHARReadVector()</tt> to
read a line, and then parse it into numerical datatypes using the
corresponding routine in the \ref StringConvert.c.
Conversion stops when the routine encounters a character that cannot
be parsed as part of a number.  If \c strict is 0, the routine
will fail only due to invalid arguments or memory allocation failure,
not from a poorly-formatted input stream; if no numbers are read,
<tt>*vector</tt> will remain \c NULL, but no error will be reported.
(In this mode, the calling routine should always test the output
before trying to dereference it, in order to avoid segmentation
violations.)  If \c strict is nonzero, the routine will report an
error if the input stream was poorly formatted, either an \c ELEN
error if no numbers were read, or \c EFMT if a character was
encountered that was neither part of a parseable number nor
whitespace.

Note that \c strict=0 allows an input stream to contain blank
lines or comments.  A comment begins with any character that cannot
occur in a valid number, which will cause the numerical parser to skip
the rest of the line.  The usual comment delimiters are <tt>'#'</tt> and
<tt>'%'</tt>, but any character except <tt>'+'</tt> <tt>'-'</tt>,
<tt>'e'</tt>, <tt>'E'</tt>, <tt>'.'</tt>, digits, and whitespace will work.

*/
/*@{*/
void LALCHARReadVector( LALStatus  *status, CHARVector **vector, FILE *stream );
void LALI2ReadVector( LALStatus  *status, INT2Vector **vector, FILE *stream, BOOLEAN strict );
void LALI4ReadVector( LALStatus  *status, INT4Vector **vector, FILE *stream, BOOLEAN strict );
void LALI8ReadVector( LALStatus  *status, INT8Vector **vector, FILE *stream, BOOLEAN strict );
void LALU2ReadVector( LALStatus  *status, UINT2Vector **vector, FILE *stream, BOOLEAN strict );
void LALU4ReadVector( LALStatus  *status, UINT4Vector **vector, FILE *stream, BOOLEAN strict );
void LALU8ReadVector( LALStatus  *status, UINT8Vector **vector, FILE *stream, BOOLEAN strict );
void LALSReadVector( LALStatus  *status, REAL4Vector **vector, FILE *stream, BOOLEAN strict );
void LALDReadVector( LALStatus  *status, REAL8Vector **vector, FILE *stream, BOOLEAN strict );
/*@}*/


/**
   \defgroup StreamVectorSequenceInput_c Module StreamVectorSequenceInput.c
   \ingroup StreamInput_h
   \author Creighton, T. D.

   \brief Reads the entire contents of an input stream into a vector sequence.

\heading{Description}

These routines read data from the I/O stream <tt>*stream</tt> until the
end-of-input is reached.  Each line is stored as a data vector, and
the vectors are combined into a LAL vector sequence structure
<tt>**sequence</tt>.  Each line vector is padded with zeros to match the
length of the longest line.  The routine passes back a pointer to the
new structure.

The routine <tt>LALCHARReadVectorSequence()</tt> essentially stores an
image of the I/O stream as a sequence of lines padded with <tt>'\0'</tt>
characters.  However, it will skip over any empty lines, which occur,
for instance, when the end-of-input or a null character <tt>'\0'</tt>
occurs immediately following a newline character <tt>'\n'</tt>.  The
numeric routines will additionally skip blank lines, comment lines, or
other input lines that have no parseable numbers in them.  (As with
the routines in \ref StreamVectorInput.c, comment in sindicated by a
<tt>#</tt> sign at the beginning of a line or a <tt>%</tt> sign anywhere
in the line, signifying that the remainder of the line is to be
ignored.)  However, if an input line contains \e any parseable
data, then the corresponding vector in the vector sequence will be
allocated (and padded with zeros, if it is shorter than the longest
line).

\heading{Algorithm}

These functions first create a linked list of vectors, using the
routines in \ref StreamVectorInput.c to read them in.  Once the list
is complete, the longest vector length is determined, and the vector
sequence is created and filled.

The numeric routines skip over blank, comment, or otherwise
unparseable lines by catching and handling the \c LEN error code
generated by the vector input routine.  However, it is worth pointing
out that the vector input routine will have generated an error message
if the error reporting bit in \c lalDebugLevel was set.  The
vector sequence input routines will therefore generate a followup
messages indicating that the preceding error was successfully dealt
with.  So you may see pairs of \c ABORT: and \c CONTINUE:
error messages when reading files containing blank or comment lines.

*/
/*@{*/
void LALCHARReadVectorSequence( LALStatus  *status, CHARVectorSequence **sequence, FILE *stream );
void LALI2ReadVectorSequence( LALStatus  *status, INT2VectorSequence **sequence, FILE *stream );
void LALI4ReadVectorSequence( LALStatus  *status, INT4VectorSequence **sequence, FILE *stream );
void LALI8ReadVectorSequence( LALStatus  *status, INT8VectorSequence **sequence, FILE *stream );
void LALU2ReadVectorSequence( LALStatus  *status, UINT2VectorSequence **sequence, FILE *stream );
void LALU4ReadVectorSequence( LALStatus  *status, UINT4VectorSequence **sequence, FILE *stream );
void LALU8ReadVectorSequence( LALStatus  *status, UINT8VectorSequence **sequence, FILE *stream );
void LALSReadVectorSequence( LALStatus  *status, REAL4VectorSequence **sequence, FILE *stream );
void LALDReadVectorSequence( LALStatus  *status, REAL8VectorSequence **sequence, FILE *stream );
/*@}*/


/**
   \defgroup StreamSequenceInput_c Module StreamSequenceInput.c
   \ingroup StreamInput_h
   \author Creighton, T. D.

   \brief Converts an input stream into a data sequence.

\heading{Description}

These routines read data from the I/O stream <tt>*stream</tt> until the
end-of-input is reached.  (The input can be of arbitrary length; the
data is temporarily stored in a linked list of buffers.)  Once read, a
LAL sequence structure <tt>**sequence</tt> is created and the data
stored in it.  The routine passes back a pointer to the new structure.

The routine <tt>LALCHARReadSequence()</tt> simply stores the entire
remaining contents of the I/O stream in a \c CHARSequence,
including whitespace, newline <tt>'\n'</tt>, null <tt>'\0'</tt>, or other
special characters.  (It can in principle be used to read and store
binary data as a sequence of bytes.  Note that the end-of-transmission
byte <tt>'\004'</tt> does \e not necessarily mark the end-of-input,
which is instead determined using the <tt>feof()</tt> function.)

The other routines in this module interpret the input as a sequence of
whitespace-separated numbers, which are parsed directly from the I/O
stream using <tt>fscanf()</tt>.  The sequence is terminated at the
end-of-input or at any point where <tt>fscanf()</tt> is unable to parse
the input.

For the complex input routines <tt>LALCReadSequence()</tt> and
<tt>LALZReadSequence()</tt>, each pair of numbers read are interpreted
as the real and imaginary parts of a complex number.  The usual input
format is for each line to contain a pair of numbers, but
<tt>fscanf()</tt> does not distinguish between newline and other
whitespace characters, so neither do these routines.

Unlike the numerical routines in other \ref StreamInput.h modules,
these routines have no mechanism to deal with comments; every
whitespace-delimited substring will be treated as a number.

\heading{Algorithm}

These routines read data into a linked list of buffers, to allow
memory allocation to occur in batches for improved efficiency.  The
numerical routines also use <tt>fscanf()</tt> directly on the I/O stream
to avoid the inefficiency of storing and parsing intermediate
character strings, as is done by the corresponding vector sequence
input routines.  This reduces robustness and versatility (as
indicated, for instance, by the inability of dealing with comments),
and increases the number of potential points-of-failure (by requiring
a consistent implementation across platforms of <tt>getc()</tt> and
<tt>fscanf()</tt>, rather than the single function <tt>fgets()</tt> used
by other stream input routines).  However, these sacrifices are
necessary to allow LAL applications to ingest large quantities of
numerical data efficiently.

*/
/*@{*/
int XLALCHARReadSequence( CHARSequence **sequence, FILE *stream );
void LALCHARReadSequence( LALStatus *status, CHARSequence **sequence, FILE *stream );
void LALI2ReadSequence( LALStatus *status, INT2Sequence **sequence, FILE *stream );
void LALI4ReadSequence( LALStatus *status, INT4Sequence **sequence, FILE *stream );
void LALI8ReadSequence( LALStatus *status, INT8Sequence **sequence, FILE *stream );
void LALU2ReadSequence( LALStatus *status, UINT2Sequence **sequence, FILE *stream );
void LALU4ReadSequence( LALStatus *status, UINT4Sequence **sequence, FILE *stream );
void LALU8ReadSequence( LALStatus *status, UINT8Sequence **sequence, FILE *stream );
void LALSReadSequence( LALStatus *status, REAL4Sequence **sequence, FILE *stream );
void LALDReadSequence( LALStatus *status, REAL8Sequence **sequence, FILE *stream );
void LALCReadSequence( LALStatus *status, COMPLEX8Sequence **sequence, FILE *stream );
void LALZReadSequence( LALStatus *status, COMPLEX16Sequence **sequence, FILE *stream );
/*@}*/



void
LALI2ReadTSeries( LALStatus *status, INT2TimeSeries *series, FILE *stream );
void
LALI4ReadTSeries( LALStatus *status, INT4TimeSeries *series, FILE *stream );
void
LALI8ReadTSeries( LALStatus *status, INT8TimeSeries *series, FILE *stream );
void
LALU2ReadTSeries( LALStatus *status, UINT2TimeSeries *series, FILE *stream );
void
LALU4ReadTSeries( LALStatus *status, UINT4TimeSeries *series, FILE *stream );
void
LALU8ReadTSeries( LALStatus *status, UINT8TimeSeries *series, FILE *stream );
void
LALSReadTSeries( LALStatus *status, REAL4TimeSeries *series, FILE *stream );
void
LALDReadTSeries( LALStatus *status, REAL8TimeSeries *series, FILE *stream );
void
LALCReadTSeries( LALStatus *status, COMPLEX8TimeSeries *series, FILE *stream );
void
LALZReadTSeries( LALStatus *status, COMPLEX16TimeSeries *series, FILE *stream );

void
LALI2ReadTVectorSeries( LALStatus *status, INT2TimeVectorSeries *series, FILE *stream );
void
LALI4ReadTVectorSeries( LALStatus *status, INT4TimeVectorSeries *series, FILE *stream );
void
LALI8ReadTVectorSeries( LALStatus *status, INT8TimeVectorSeries *series, FILE *stream );
void
LALU2ReadTVectorSeries( LALStatus *status, UINT2TimeVectorSeries *series, FILE *stream );
void
LALU4ReadTVectorSeries( LALStatus *status, UINT4TimeVectorSeries *series, FILE *stream );
void
LALU8ReadTVectorSeries( LALStatus *status, UINT8TimeVectorSeries *series, FILE *stream );
void
LALSReadTVectorSeries( LALStatus *status, REAL4TimeVectorSeries *series, FILE *stream );
void
LALDReadTVectorSeries( LALStatus *status, REAL8TimeVectorSeries *series, FILE *stream );
void
LALCReadTVectorSeries( LALStatus *status, COMPLEX8TimeVectorSeries *series, FILE *stream );
void
LALZReadTVectorSeries( LALStatus *status, COMPLEX16TimeVectorSeries *series, FILE *stream );

void
LALI2ReadTArraySeries( LALStatus *status, INT2TimeArraySeries *series, FILE *stream );
void
LALI4ReadTArraySeries( LALStatus *status, INT4TimeArraySeries *series, FILE *stream );
void
LALI8ReadTArraySeries( LALStatus *status, INT8TimeArraySeries *series, FILE *stream );
void
LALU2ReadTArraySeries( LALStatus *status, UINT2TimeArraySeries *series, FILE *stream );
void
LALU4ReadTArraySeries( LALStatus *status, UINT4TimeArraySeries *series, FILE *stream );
void
LALU8ReadTArraySeries( LALStatus *status, UINT8TimeArraySeries *series, FILE *stream );
void
LALSReadTArraySeries( LALStatus *status, REAL4TimeArraySeries *series, FILE *stream );
void
LALDReadTArraySeries( LALStatus *status, REAL8TimeArraySeries *series, FILE *stream );
void
LALCReadTArraySeries( LALStatus *status, COMPLEX8TimeArraySeries *series, FILE *stream );
void
LALZReadTArraySeries( LALStatus *status, COMPLEX16TimeArraySeries *series, FILE *stream );

void
LALI2ReadFSeries( LALStatus *status, INT2FrequencySeries *series, FILE *stream );
void
LALI4ReadFSeries( LALStatus *status, INT4FrequencySeries *series, FILE *stream );
void
LALI8ReadFSeries( LALStatus *status, INT8FrequencySeries *series, FILE *stream );
void
LALU2ReadFSeries( LALStatus *status, UINT2FrequencySeries *series, FILE *stream );
void
LALU4ReadFSeries( LALStatus *status, UINT4FrequencySeries *series, FILE *stream );
void
LALU8ReadFSeries( LALStatus *status, UINT8FrequencySeries *series, FILE *stream );
void
LALSReadFSeries( LALStatus *status, REAL4FrequencySeries *series, FILE *stream );
void
LALDReadFSeries( LALStatus *status, REAL8FrequencySeries *series, FILE *stream );
void
LALCReadFSeries( LALStatus *status, COMPLEX8FrequencySeries *series, FILE *stream );
void
LALZReadFSeries( LALStatus *status, COMPLEX16FrequencySeries *series, FILE *stream );




void
LALI2ReadGrid( LALStatus *status, INT2Grid **grid, FILE *stream );
void
LALI4ReadGrid( LALStatus *status, INT4Grid **grid, FILE *stream );
void
LALI8ReadGrid( LALStatus *status, INT8Grid **grid, FILE *stream );
void
LALU2ReadGrid( LALStatus *status, UINT2Grid **grid, FILE *stream );
void
LALU4ReadGrid( LALStatus *status, UINT4Grid **grid, FILE *stream );
void
LALU8ReadGrid( LALStatus *status, UINT8Grid **grid, FILE *stream );
void
LALSReadGrid( LALStatus *status, REAL4Grid **grid, FILE *stream );
void
LALDReadGrid( LALStatus *status, REAL8Grid **grid, FILE *stream );
void
LALCReadGrid( LALStatus *status, COMPLEX8Grid **grid, FILE *stream );
void
LALZReadGrid( LALStatus *status, COMPLEX16Grid **grid, FILE *stream );






#ifdef __cplusplus
}
#endif

#endif /* _STREAMINPUT_H */
