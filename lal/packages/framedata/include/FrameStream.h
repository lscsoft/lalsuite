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
 * that allow you to extract frame data.  Many of these routines have names
 * similar to the standard C file stream manipulation routines and perform
 * similar functions.
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
#define FRAMESTREAMH_ENULL 00001
#define FRAMESTREAMH_ENNUL 00002
#define FRAMESTREAMH_EALOC 00004
#define FRAMESTREAMH_EFILE 00010
#define FRAMESTREAMH_EOPEN 00020
#define FRAMESTREAMH_EREAD 00040
#define FRAMESTREAMH_ETIME 00100
#define FRAMESTREAMH_ESIZE 00200
#define FRAMESTREAMH_ECHAN 00400
#define FRAMESTREAMH_ETYPE 01000
#define FRAMESTREAMH_ERROR 02000
#define FRAMESTREAMH_EDONE 04000

#define FRAMESTREAMH_MSGENULL "Null pointer"
#define FRAMESTREAMH_MSGENNUL "Non-null pointer"
#define FRAMESTREAMH_MSGEALOC "Memory allocation error"
#define FRAMESTREAMH_MSGEFILE "Frame data files not found"
#define FRAMESTREAMH_MSGEOPEN "Frame file open error"
#define FRAMESTREAMH_MSGEREAD "Frame file read error"
#define FRAMESTREAMH_MSGETIME "Invalid ADC offset time"
#define FRAMESTREAMH_MSGESIZE "Invalid vector length"
#define FRAMESTREAMH_MSGECHAN "Could not find ADC channel"
#define FRAMESTREAMH_MSGETYPE "Invalid ADC type"
#define FRAMESTREAMH_MSGERROR "Frame stream error"
#define FRAMESTREAMH_MSGEDONE "End of frame data"
/**** </lalErrTable> */

/**** <lalLaTeX>
 *
 * \subsection*{Structures}
 * \idx[Type]{FrStream}
 * \idx[Type]{FrPos}
 * \idx[Type]{ChannelType}
 * \idx[Type]{FrChanIn}
 * \idx[Type]{FrOutPar}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct tagFrStream FrStream;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 *
 * This structure details the state of the frame stream.  The contents are
 * private; you should not tamper with them!
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagFrPos
{
  LIGOTimeGPS epoch;
  UINT4       filenum;
  UINT4       frnum;
}
FrPos;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 * 
 * This structure contains a record of the state of a frame stream; this
 * record can be used to restore the stream to the state when the record
 * was made (provided the stream has not been closed).  The fields are:
 * \begin{description}
 * \item[\texttt{epoch}] the GPS time of the open frame when the record
 *     was made.
 * \item[\texttt{filenum}] the file number of a list of frame files that was
 *     open when the record was made.
 * \item[\texttt{frnum}] the frame number of the frames within the open
 *     frame file that was open when the record was made.
 * \end{description}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef enum
tagChannelType
{ ProcDataChannel, ADCDataChannel, SimDataChannel }
ChannelType;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 * 
 * These are the various types of channel that can be specified for read/write.
 * They are ``post-processed data'' (\texttt{ProcDataChannel}), ``ADC data''
 * (\texttt{ADCDataChannel}), and ``simulated data'' (\texttt{SimDataChannel}).
 *
 **** </lalLaTeX> */

/**** <lalVerbatim> */
typedef struct
tagFrChanIn
{
  const CHAR *name;
  ChannelType type;
}
FrChanIn;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 * 
 * This structure specifies the channel to read as input.  The fields are:
 * \begin{description}
 * \item[\texttt{name}] the name of the channel.
 * \item[\texttt{type}] the channel type.
 * \end{description}
 *
 **** </lalLaTeX> */

/**** <lalVerbatim> */
typedef struct
tagFrOutPar
{
  const CHAR *prefix;
  ChannelType type;
  UINT4 nframes;
  UINT4 frame;
  UINT4 run;
}
FrOutPar;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 * 
 * This structure specifies the parameters for output of data to a frame.
 * The fields are:
 * \begin{description}
 * \item[\texttt{prefix}] the prefix to attach to the output frame file name.
 * \item[\texttt{type}] the type of channel to create in the output frames.
 * \item[\texttt{nframes}] the number of frames to output in the frame file.
 * \item[\texttt{frame}] the number the first frame of output.
 * \item[\texttt{run}] the number this data run.
 * \end{description}
 *
 * \vfill{\footnotesize\input{FrameStreamHV}}
 * \newpage\input{FrameStreamC}
 * \newpage\input{FrameStreamTestC}
 *
 **** </lalLaTeX> */

void
LALFrOpen(
    LALStatus    *status,
    FrStream    **stream,
    const CHAR   *dirname,
    const CHAR   *pattern
    );

void
LALFrClose(
    LALStatus  *status,
    FrStream  **stream
    );

void
LALFrEnd(
    LALStatus *status,
    INT4      *end,
    FrStream  *stream
    );

void
LALFrNext(
    LALStatus *status,
    FrStream  *stream
    );

void
LALFrRewind( 
    LALStatus *status,
    FrStream  *stream
    );

void
LALFrSeek(
    LALStatus   *status,
    LIGOTimeGPS *epoch,
    FrStream    *stream
    );

void
LALFrTell(
    LALStatus   *status,
    LIGOTimeGPS *epoch,
    FrStream    *stream
    );

void
LALFrGetPos(
    LALStatus *status,
    FrPos     *position,
    FrStream  *stream
    );

void
LALFrSetPos(
    LALStatus *status,
    FrPos     *position,
    FrStream  *stream
    );

void
LALFrGetINT2TimeSeries(
    LALStatus      *status,
    INT2TimeSeries *series,
    FrChanIn       *chanin,
    FrStream       *stream
    );

void
LALFrGetREAL4TimeSeries(
    LALStatus       *status,
    REAL4TimeSeries *series,
    FrChanIn        *chanin,
    FrStream        *stream
    );

void
LALFrWriteREAL4TimeSeries(
    LALStatus       *status,
    REAL4TimeSeries *series,
    FrOutPar        *params
    );

#ifdef __cplusplus
#pragma {
}
#endif

#endif /* _FRAMESTREAM_H */
