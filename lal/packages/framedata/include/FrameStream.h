/**** <lalVerbatim file="FrameStreamHV">
 * Author: Jolien D. E. Creighton
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \section{Header \texttt{FrameStream.h}}
 *
 * Low-level routines for manupulating frame data streams.
 *
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <stdio.h>
 * #include <lal/FrameStream.h>
 * \end{verbatim}
 *
 * A frame stream is like a file stream except that it streams along the set
 * of frames in a set of frame files.  These routines are low-level routines
 * that allow you to extract frame data.
 *
 **** </lalLaTeX> */

#include <lal/LALDatatypes.h>

#ifndef _FRAMESTREAM_H
#define _FRAMESTREAM_H

#ifdef __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( FRAMESTREAMH, "$Id$" );

/**** <lalLaTeX>
 *
 * \subsection*{Error conditions}
 *
 **** </lalLaTeX> */
/**** <lalErrTable> */
#define FRAMESTREAMH_ENULL 1
#define FRAMESTREAMH_ENNUL 2
#define FRAMESTREAMH_EALOC 4
#define FRAMESTREAMH_EPIPE 8
#define FRAMESTREAMH_EFILE 16
#define FRAMESTREAMH_EOPEN 32
#define FRAMESTREAMH_EREAD 64
#define FRAMESTREAMH_ECHAN 128
#define FRAMESTREAMH_ETYPE 256
#define FRAMESTREAMH_ERROR 512

#define FRAMESTREAMH_MSGENULL "Null pointer"
#define FRAMESTREAMH_MSGENNUL "Non-null pointer"
#define FRAMESTREAMH_MSGEALOC "Memory allocation error"
#define FRAMESTREAMH_MSGEPIPE "Pipe open error"
#define FRAMESTREAMH_MSGEFILE "Frame data files not found"
#define FRAMESTREAMH_MSGEOPEN "Frame file open error"
#define FRAMESTREAMH_MSGEREAD "Frame file read error"
#define FRAMESTREAMH_MSGECHAN "Could not find ADC channel"
#define FRAMESTREAMH_MSGETYPE "Invalid ADC type"
#define FRAMESTREAMH_MSGERROR "Frame stream error"
/**** </lalErrTable> */

/**** <lalLaTeX>
 *
 * \subsection*{Structures}
 * \idx[Type]{FrameStream}
 * \idx[Type]{FrameStreamPos}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct tagFrameStream FrameStream;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 *
 * This structure details the state of the frame stream.  The contents are
 * private; you should not tamper with them!
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagFrameStreamPos
{
  LIGOTimeGPS gpstime;
  UINT4       filenum;
  UINT4       frnum;
}
FrameStreamPos;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 * 
 * This structure contains a record of the state of a frame stream; this
 * record can be used to restore the stream to the state when the record
 * was made (provided the stream has not been closed).  The fields are:
 * \begin{description}
 * \item[\texttt{gpstime}] the GPS time of the open frame when the record
 *     was made.
 * \item[\texttt{filenum}] the file number of a list of frame files that was
 *     open when the record was made.
 * \item[\texttt{filenum}] the frame number of the frames within the open
 *     frame file that was open when the record was made.
 * \end{description}
 *
 * \vfill{\footnotesize\input{FrameStreamHV}}
 * \newpage\input{FrameStreamC}
 * \newpage\input{FrameStreamTestC}
 *
 **** </lalLaTeX> */


void
LALOpenFrameStream(
    LALStatus    *status,
    FrameStream **stream,
    const CHAR   *dirname,
    const CHAR   *headname
    );

void
LALCloseFrameStream(
    LALStatus    *status,
    FrameStream **stream
    );

void
LALSFrameReadADCTimeSeries(
    LALStatus        *status,
    REAL4TimeSeries **series,
    const CHAR       *channel,
    FrameStream      *stream
    );

void
LALI2FrameReadADCTimeSeries(
    LALStatus       *status,
    INT2TimeSeries **series,
    const CHAR      *channel,
    FrameStream     *stream
    );

void
LALNextFrame(
    LALStatus   *status,
    INT4        *error,
    FrameStream *stream
    );

void
LALFrameStreamError(
    LALStatus   *status,
    INT4        *error,
    FrameStream *stream
    );

void
LALFrameStreamGetPos(
    LALStatus      *status,
    FrameStreamPos *position,
    FrameStream    *stream
    );

void
LALFrameStreamSetPos(
    LALStatus      *status,
    FrameStreamPos *position,
    FrameStream    *stream
    );

#ifdef __cplusplus
#pragma {
}
#endif

#endif /* _FRAMESTREAM_H */
