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
#include <lal/FrameCache.h>

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
#define FRAMESTREAMH_ETREQ 010000
#define FRAMESTREAMH_EDGAP 020000

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
#define FRAMESTREAMH_MSGETREQ "No data at time requested"
#define FRAMESTREAMH_MSGEDGAP "Gap in the data"
/**** </lalErrTable> */

/**** <lalLaTeX>
 *
 * \subsection*{Structures}
 * \idx[Type]{FrState}
 * \idx[Type]{FrFileInfo}
 * \idx[Type]{FrStream}
 * \idx[Type]{FrPos}
 * \idx[Type]{ChannelType}
 * \idx[Type]{FrChanIn}
 * \idx[Type]{FrOutPar}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef enum
{
  LAL_FR_OK  = 0,  /* nominal */
  LAL_FR_ERR = 1,  /* error in frame stream */
  LAL_FR_END = 2,  /* end of frame stream */
  LAL_FR_GAP = 4,  /* gap in frame stream */
  LAL_FR_URL = 8,  /* error opening frame URL */
  LAL_FR_TOC = 16  /* error reading frame TOC */
}
FrState;
typedef enum
{
  LAL_FR_SILENT_MODE     = 0,
  LAL_FR_TIMEWARN_MODE   = 1,  /* display warning for invalid time requests */
  LAL_FR_GAPINFO_MODE    = 2,  /* display info for gaps in data */
  LAL_FR_VERBOSE_MODE    = 3,  /* display warnings and info */
  LAL_FR_IGNOREGAP_MODE  = 4,  /* ignore gaps in data */
  LAL_FR_IGNORETIME_MODE = 8,  /* ignore invalid times requested */
  LAL_FR_DEFAULT_MODE    = 15  /* ignore time/gaps but report warnings & info */
}
FrMode;
struct FrFile;
typedef struct tagFrFileInfo
{
  INT4  ind;
  CHAR *url;
  INT4  t0;
  INT4  dt;
}
FrFileInfo;
typedef struct tagFrStream
{
  FrState        state;
  INT4           mode;
  LIGOTimeGPS    epoch;
  UINT4          nfile;
  FrFileInfo    *flist;
  UINT4          fnum;
  struct FrFile *file;
  INT4           pos;
}
FrStream;
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
  UINT4       fnum;
  INT4        pos;
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
 * \item[\texttt{fnum}] the file number of a list of frame files that was
 *     open when the record was made.
 * \item[\texttt{pos}] the position within the
 *     frame file that was open when the record was made.
 * \end{description}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef enum
{ LAL_ADC_CHAN, LAL_SIM_CHAN, LAL_PROC_CHAN }
FrChanType;
/* for backwards compatability... */
#define ChannelType FrChanType
#define ProcDataChannel LAL_PROC_CHAN
#define ADCDataChannel  LAL_ADC_CHAN
#define SimDataChannel  LAL_SIM_CHAN
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
  const CHAR *source;
  const CHAR *description;
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
 * \item[\texttt{source}] the source identifier to attach to the output
 *   frame file name.
 * \item[\texttt{description}] the description identifier to attach to the
 *   output frame file name.
 * \item[\texttt{type}] the type of channel to create in the output frames.
 * \item[\texttt{nframes}] the number of frames to output in the frame file.
 * \item[\texttt{frame}] the number the first frame of output.
 * \item[\texttt{run}] the number this data run.
 * \end{description}
 * The output frame file name will be
 * $\langle\mbox{source}\rangle$\verb+-+$\langle\mbox{description}\rangle$%
 * \verb+-+$\langle\mbox{GPS start time}\rangle$\verb+-+%
 * $\langle\mbox{duration}\rangle$\verb+.gwf+.
 *
 * \vfill{\footnotesize\input{FrameStreamHV}}
 * \newpage\input{FrameStreamC}
 * \newpage\input{FrameSeriesC}
 * \newpage\input{FrameStreamTestC}
 *
 **** </lalLaTeX> */

FrStream * XLALFrCacheOpen( FrCache *cache );
FrStream * XLALFrOpen( const char *dirname, const char *pattern );
int XLALFrClose( FrStream *stream );
int XLALFrSetMode( FrStream *stream, int mode );
int XLALFrState( FrStream *stream );
int XLALFrClearErr( FrStream *stream );
int XLALFrRewind( FrStream *stream );
int XLALFrNext( FrStream *stream );
int XLALFrSeek( FrStream *stream, const LIGOTimeGPS *epoch );
int XLALFrTell( LIGOTimeGPS *epoch, FrStream *stream );
int XLALFrGetpos( FrPos *position, FrStream *stream );
int XLALFrSetpos( FrStream *stream, FrPos *position );
int XLALFrGetTimeSeriesType( const char *channel, FrStream *stream );
int XLALFrGetINT2TimeSeries( INT2TimeSeries *series, FrStream *stream );
int XLALFrGetINT4TimeSeries( INT4TimeSeries *series, FrStream *stream );
int XLALFrGetINT8TimeSeries( INT8TimeSeries *series, FrStream *stream );
int XLALFrGetREAL4TimeSeries( REAL4TimeSeries *series, FrStream *stream );
int XLALFrGetREAL8TimeSeries( REAL8TimeSeries *series, FrStream *stream );
int XLALFrGetCOMPLEX8TimeSeries( COMPLEX8TimeSeries *series, FrStream *stream );
int XLALFrGetCOMPLEX16TimeSeries( COMPLEX16TimeSeries *series, FrStream *stream );
int XLALFrGetINT2TimeSeriesMetadata( INT2TimeSeries *series, FrStream *stream );
int XLALFrGetINT4TimeSeriesMetadata( INT4TimeSeries *series, FrStream *stream );
int XLALFrGetINT8TimeSeriesMetadata( INT8TimeSeries *series, FrStream *stream );
int XLALFrGetREAL4TimeSeriesMetadata( REAL4TimeSeries *series, FrStream *stream );
int XLALFrGetREAL8TimeSeriesMetadata( REAL8TimeSeries *series, FrStream *stream );
int XLALFrGetCOMPLEX8TimeSeriesMetadata( COMPLEX8TimeSeries *series, FrStream *stream );
int XLALFrGetCOMPLEX16TimeSeriesMetadata( COMPLEX16TimeSeries *series, FrStream *stream );
int XLALFrGetINT2FrequencySeries( INT2FrequencySeries *series, FrStream *stream );
int XLALFrGetINT4FrequencySeries( INT4FrequencySeries *series, FrStream *stream );
int XLALFrGetINT8FrequencySeries( INT8FrequencySeries *series, FrStream *stream );
int XLALFrGetREAL4FrequencySeries( REAL4FrequencySeries *series, FrStream *stream );
int XLALFrGetREAL8FrequencySeries( REAL8FrequencySeries *series, FrStream *stream );
int XLALFrGetCOMPLEX8FrequencySeries( COMPLEX8FrequencySeries *series, FrStream *stream );
int XLALFrGetCOMPLEX16FrequencySeries( COMPLEX16FrequencySeries *series, FrStream *stream );


/*
 *
 * LAL Routines.
 *
 */

void
LALFrCacheOpen(
    LALStatus  *status,
    FrStream  **output,
    FrCache    *cache
    );

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
LALFrSetMode(
    LALStatus *status,
    INT4       mode,
    FrStream  *stream
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
    LALStatus         *status,
    const LIGOTimeGPS *epoch,
    FrStream          *stream
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
LALFrGetTimeSeriesType(
    LALStatus   *status,
    LALTYPECODE *output,
    FrChanIn    *chanin,
    FrStream    *stream
    );

void
LALFrGetINT2TimeSeries(
    LALStatus      *status,
    INT2TimeSeries *series,
    FrChanIn       *chanin,
    FrStream       *stream
    );

void
LALFrGetINT2TimeSeriesMetadata(
    LALStatus      *status,
    INT2TimeSeries *series,
    FrChanIn       *chanin,
    FrStream       *stream
    );

void
LALFrGetINT4TimeSeries(
    LALStatus      *status,
    INT4TimeSeries *series,
    FrChanIn       *chanin,
    FrStream       *stream
    );

void
LALFrGetINT4TimeSeriesMetadata(
    LALStatus      *status,
    INT4TimeSeries *series,
    FrChanIn       *chanin,
    FrStream       *stream
    );

void
LALFrGetINT8TimeSeries(
    LALStatus      *status,
    INT8TimeSeries *series,
    FrChanIn       *chanin,
    FrStream       *stream
    );

void
LALFrGetINT8TimeSeriesMetadata(
    LALStatus      *status,
    INT8TimeSeries *series,
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
LALFrGetREAL4TimeSeriesMetadata(
    LALStatus       *status,
    REAL4TimeSeries *series,
    FrChanIn        *chanin,
    FrStream        *stream
    );

void
LALFrGetREAL8TimeSeries(
    LALStatus       *status,
    REAL8TimeSeries *series,
    FrChanIn        *chanin,
    FrStream        *stream
    );

void
LALFrGetREAL8TimeSeriesMetadata(
    LALStatus       *status,
    REAL8TimeSeries *series,
    FrChanIn        *chanin,
    FrStream        *stream
    );

void
LALFrGetCOMPLEX8TimeSeries(
    LALStatus          *status,
    COMPLEX8TimeSeries *series,
    FrChanIn           *chanin,
    FrStream           *stream
    );

void
LALFrGetCOMPLEX8TimeSeriesMetadata(
    LALStatus          *status,
    COMPLEX8TimeSeries *series,
    FrChanIn           *chanin,
    FrStream           *stream
    );

void
LALFrGetCOMPLEX16TimeSeries(
    LALStatus           *status,
    COMPLEX16TimeSeries *series,
    FrChanIn            *chanin,
    FrStream            *stream
    );

void
LALFrGetCOMPLEX16TimeSeriesMetadata(
    LALStatus           *status,
    COMPLEX16TimeSeries *series,
    FrChanIn            *chanin,
    FrStream            *stream
    );


void
LALFrGetINT2FrequencySeries(
    LALStatus      *status,
    INT2FrequencySeries *series,
    FrChanIn       *chanin,
    FrStream       *stream
    );

void
LALFrGetINT4FrequencySeries(
    LALStatus      *status,
    INT4FrequencySeries *series,
    FrChanIn       *chanin,
    FrStream       *stream
    );

void
LALFrGetINT8FrequencySeries(
    LALStatus      *status,
    INT8FrequencySeries *series,
    FrChanIn       *chanin,
    FrStream       *stream
    );

void
LALFrGetREAL4FrequencySeries(
    LALStatus       *status,
    REAL4FrequencySeries *series,
    FrChanIn        *chanin,
    FrStream        *stream
    );

void
LALFrGetREAL8FrequencySeries(
    LALStatus       *status,
    REAL8FrequencySeries *series,
    FrChanIn        *chanin,
    FrStream        *stream
    );

void
LALFrGetCOMPLEX8FrequencySeries(
    LALStatus          *status,
    COMPLEX8FrequencySeries *series,
    FrChanIn           *chanin,
    FrStream           *stream
    );

void
LALFrGetCOMPLEX16FrequencySeries(
    LALStatus           *status,
    COMPLEX16FrequencySeries *series,
    FrChanIn            *chanin,
    FrStream            *stream
    );


void
LALFrWriteINT2TimeSeries(
    LALStatus       *status,
    INT2TimeSeries  *series,
    FrOutPar        *params
    );

void
LALFrWriteINT2TimeSeries(
    LALStatus       *status,
    INT2TimeSeries  *series,
    FrOutPar        *params
    );

void
LALFrWriteINT4TimeSeries(
    LALStatus       *status,
    INT4TimeSeries  *series,
    FrOutPar        *params
    );

void
LALFrWriteINT8TimeSeries(
    LALStatus       *status,
    INT8TimeSeries  *series,
    FrOutPar        *params
    );

void
LALFrWriteREAL4TimeSeries(
    LALStatus       *status,
    REAL4TimeSeries *series,
    FrOutPar        *params
    );

void
LALFrWriteREAL8TimeSeries(
    LALStatus       *status,
    REAL8TimeSeries *series,
    FrOutPar        *params
    );

void
LALFrWriteCOMPLEX8TimeSeries(
    LALStatus          *status,
    COMPLEX8TimeSeries *series,
    FrOutPar           *params
    );

void
LALFrWriteCOMPLEX16TimeSeries(
    LALStatus           *status,
    COMPLEX16TimeSeries *series,
    FrOutPar            *params
    );


void
LALFrWriteINT2FrequencySeries(
    LALStatus       *status,
    INT2FrequencySeries  *series,
    FrOutPar        *params,
    INT4             subtype
    );

void
LALFrWriteINT2FrequencySeries(
    LALStatus       *status,
    INT2FrequencySeries  *series,
    FrOutPar        *params,
    INT4             subtype
    );

void
LALFrWriteINT4FrequencySeries(
    LALStatus       *status,
    INT4FrequencySeries  *series,
    FrOutPar        *params,
    INT4             subtype
    );

void
LALFrWriteINT8FrequencySeries(
    LALStatus       *status,
    INT8FrequencySeries  *series,
    FrOutPar        *params,
    INT4             subtype
    );

void
LALFrWriteREAL4FrequencySeries(
    LALStatus       *status,
    REAL4FrequencySeries *series,
    FrOutPar        *params,
    INT4             subtype
    );

void
LALFrWriteREAL8FrequencySeries(
    LALStatus       *status,
    REAL8FrequencySeries *series,
    FrOutPar        *params,
    INT4             subtype
    );

void
LALFrWriteCOMPLEX8FrequencySeries(
    LALStatus          *status,
    COMPLEX8FrequencySeries *series,
    FrOutPar           *params,
    INT4             subtype
    );

void
LALFrWriteCOMPLEX16FrequencySeries(
    LALStatus           *status,
    COMPLEX16FrequencySeries *series,
    FrOutPar            *params,
    INT4             subtype
    );

#ifdef __cplusplus
#pragma {
}
#endif

#endif /* _FRAMESTREAM_H */
