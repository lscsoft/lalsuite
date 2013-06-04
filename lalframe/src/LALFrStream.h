/*
*  Copyright (C) 2007 Jolien Creighton, Kipp Cannon
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
#include <lal/LALDatatypes.h>
#include <lal/LALCache.h>

#ifndef _LALFRSTREAM_H
#define _LALFRSTREAM_H

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \defgroup FrameStream_h Header FrameStream.h
 * \ingroup pkg_framedata
 *
 * \author Jolien D. E. Creighton
 *
 * \brief Low-level routines for manupulating frame data streams.
 *
 * \heading{Synopsis}
 * \code
 * #include <stdio.h>
 * #include <lal/LALFrStream.h>
 * \endcode
 *
 * A frame stream is like a file stream except that it streams along the set
 * of frames in a set of frame files.  These routines are low-level routines
 * that allow you to extract frame data.  Many of these routines have names
 * similar to the standard C file stream manipulation routines and perform
 * similar functions.
 *
*/
/*@{*/

/**\name Error Codes */
/*@{*/
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
/*@}*/

typedef enum
{
  LAL_FR_STREAM_OK  = 0,  /* nominal */
  LAL_FR_STREAM_ERR = 1,  /* error in frame stream */
  LAL_FR_STREAM_END = 2,  /* end of frame stream */
  LAL_FR_STREAM_GAP = 4,  /* gap in frame stream */
  LAL_FR_STREAM_URL = 8,  /* error opening frame URL */
  LAL_FR_STREAM_TOC = 16  /* error reading frame TOC */
}
FrState;
typedef enum
{
  LAL_FR_STREAM_SILENT_MODE     = 0,
  LAL_FR_STREAM_TIMEWARN_MODE   = 1,  /* display warning for invalid time requests */
  LAL_FR_STREAM_GAPINFO_MODE    = 2,  /* display info for gaps in data */
  LAL_FR_STREAM_VERBOSE_MODE    = 3,  /* display warnings and info */
  LAL_FR_STREAM_IGNOREGAP_MODE  = 4,  /* ignore gaps in data */
  LAL_FR_STREAM_IGNORETIME_MODE = 8,  /* ignore invalid times requested */
  LAL_FR_STREAM_DEFAULT_MODE    = 15, /* ignore time/gaps but report warnings & info */
  LAL_FR_STREAM_CHECKSUM_MODE   = 16  /* ensure that file checksums are OK */
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

/**
 *
 * This structure details the state of the frame stream.  The contents are
 * private; you should not tamper with them!
 *
*/
typedef struct tagLALFrStream
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
LALFrStream;


/**
 *
 * This structure contains a record of the state of a frame stream; this
 * record can be used to restore the stream to the state when the record
 * was made (provided the stream has not been closed).  The fields are:
 * <dl>
 * <dt>epoch</dt><dd> the GPS time of the open frame when the record  was made.</dd>
 * <dt>fnum</dt><dd> the file number of a list of frame files that was open when the record was made.</dd>
 * <dt>pos</dt><dd> the position within the frame file that was open when the record was made.</dd>
 * </dl>
 *
*/
typedef struct
tagLALFrStreamPos
{
  LIGOTimeGPS epoch;
  UINT4       fnum;
  INT4        pos;
}
LALFrStreamPos;

typedef enum
{ LAL_ADC_CHAN, LAL_SIM_CHAN, LAL_PROC_CHAN }
FrChanType;
/* for backwards compatability... */
#define ChannelType FrChanType
#define ProcDataChannel LAL_PROC_CHAN
#define ADCDataChannel  LAL_ADC_CHAN
#define SimDataChannel  LAL_SIM_CHAN

/**
 *
 * These are the various types of channel that can be specified for read/write.
 * They are "post-processed data" (\c ProcDataChannel), "ADC data"
 * (\c ADCDataChannel), and "simulated data" (\c SimDataChannel).
 *
*/


/**
 *
 * This structure specifies the channel to read as input.  The fields are:
 * <dl>
 * <dt>name</dt><dd> the name of the channel.
 * </dd><dt>type</dt><dd> the channel type.
 * </dd></dl>
 *
*/
#ifdef SWIG /* SWIG interface directives */
SWIGLAL(STRUCT_IMMUTABLE(tagFrChanIn, name));
#endif /* SWIG */
typedef struct
tagFrChanIn
{
  const CHAR *name;
  ChannelType type;
}
FrChanIn;


/**
 *
 * This structure specifies the parameters for output of data to a frame.
 * The fields are:
 * <dl>
 * <dt>source</dt><dd> the source identifier to attach to the output frame file name.</dd>
 * <dt>description</dt><dd> the description identifier to attach to the output frame file name.</dd>
 * <dt>type</dt><dd> the type of channel to create in the output frames.</dd>
 * <dt>nframes</dt><dd> the number of frames to output in the frame file.</dd>
 * <dt>frame</dt><dd> the number the first frame of output.</dd>
 * <dt>run</dt><dd> the number this data run.</dd>
 * </dl>
 * The output frame file name will be
 * \f$\langle\mbox{source}\rangle\f$<tt>-</tt>\f$\langle\mbox{description}\rangle\f$%
 * <tt>-</tt>\f$\langle\mbox{GPS start time}\rangle\f$<tt>-</tt>%
 * \f$\langle\mbox{duration}\rangle\f$<tt>.gwf</tt>.
 *
*/
#ifdef SWIG /* SWIG interface directives */
SWIGLAL(STRUCT_IMMUTABLE(tagFrOutPar, source, description));
#endif /* SWIG */
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


LALFrStream * XLALFrStreamCacheOpen( LALCache *cache );
LALFrStream * XLALFrStreamOpen( const char *dirname, const char *pattern );
void XLALFrStreamClose( LALFrStream *stream );
int XLALFrStreamSetMode( LALFrStream *stream, int mode );
int XLALFrStreamGetState( LALFrStream *stream );
int XLALFrStreamClearErr( LALFrStream *stream );
int XLALFrStreamRewind( LALFrStream *stream );
int XLALFrStreamNext( LALFrStream *stream );
int XLALFrStreamSeek( LALFrStream *stream, const LIGOTimeGPS *epoch );
int XLALFrStreamTell( LIGOTimeGPS *epoch, LALFrStream *stream );
int XLALFrStreamGetpos( LALFrStreamPos *position, LALFrStream *stream );
int XLALFrStreamSetpos( LALFrStream *stream, LALFrStreamPos *position );
int XLALFrStreamGetTimeSeriesType( const char *channel, LALFrStream *stream );

int XLALFrStreamGetINT2TimeSeries( INT2TimeSeries *series, LALFrStream *stream );
int XLALFrStreamGetINT4TimeSeries( INT4TimeSeries *series, LALFrStream *stream );
int XLALFrStreamGetINT8TimeSeries( INT8TimeSeries *series, LALFrStream *stream );
int XLALFrStreamGetREAL4TimeSeries( REAL4TimeSeries *series, LALFrStream *stream );
int XLALFrStreamGetREAL8TimeSeries( REAL8TimeSeries *series, LALFrStream *stream );
int XLALFrStreamGetCOMPLEX8TimeSeries( COMPLEX8TimeSeries *series, LALFrStream *stream );
int XLALFrStreamGetCOMPLEX16TimeSeries( COMPLEX16TimeSeries *series, LALFrStream *stream );

int XLALFrStreamGetINT2TimeSeriesMetadata( INT2TimeSeries *series, LALFrStream *stream );
int XLALFrStreamGetINT4TimeSeriesMetadata( INT4TimeSeries *series, LALFrStream *stream );
int XLALFrStreamGetINT8TimeSeriesMetadata( INT8TimeSeries *series, LALFrStream *stream );
int XLALFrStreamGetREAL4TimeSeriesMetadata( REAL4TimeSeries *series, LALFrStream *stream );
int XLALFrStreamGetREAL8TimeSeriesMetadata( REAL8TimeSeries *series, LALFrStream *stream );
int XLALFrStreamGetCOMPLEX8TimeSeriesMetadata( COMPLEX8TimeSeries *series, LALFrStream *stream );
int XLALFrStreamGetCOMPLEX16TimeSeriesMetadata( COMPLEX16TimeSeries *series, LALFrStream *stream );

int XLALFrStreamGetINT2FrequencySeries( INT2FrequencySeries *series, LALFrStream *stream );
int XLALFrStreamGetINT4FrequencySeries( INT4FrequencySeries *series, LALFrStream *stream );
int XLALFrStreamGetINT8FrequencySeries( INT8FrequencySeries *series, LALFrStream *stream );
int XLALFrStreamGetREAL4FrequencySeries( REAL4FrequencySeries *series, LALFrStream *stream );
int XLALFrStreamGetREAL8FrequencySeries( REAL8FrequencySeries *series, LALFrStream *stream );
int XLALFrStreamGetCOMPLEX8FrequencySeries( COMPLEX8FrequencySeries *series, LALFrStream *stream );
int XLALFrStreamGetCOMPLEX16FrequencySeries( COMPLEX16FrequencySeries *series, LALFrStream *stream );
int XLALFrStreamGetVectorLength ( CHAR *name, LALFrStream *stream );

INT2TimeSeries *XLALFrStreamReadINT2TimeSeries( LALFrStream *stream, const char *chname, const LIGOTimeGPS *start, REAL8 duration, size_t lengthlimit );
INT4TimeSeries *XLALFrStreamReadINT4TimeSeries( LALFrStream *stream, const char *chname, const LIGOTimeGPS *start, REAL8 duration, size_t lengthlimit );
INT8TimeSeries *XLALFrStreamReadINT8TimeSeries( LALFrStream *stream, const char *chname, const LIGOTimeGPS *start, REAL8 duration, size_t lengthlimit );
REAL4TimeSeries *XLALFrStreamReadREAL4TimeSeries( LALFrStream *stream, const char *chname, const LIGOTimeGPS *start, REAL8 duration, size_t lengthlimit );
REAL8TimeSeries *XLALFrStreamReadREAL8TimeSeries( LALFrStream *stream, const char *chname, const LIGOTimeGPS *start, REAL8 duration, size_t lengthlimit );
COMPLEX8TimeSeries *XLALFrStreamReadCOMPLEX8TimeSeries( LALFrStream *stream, const char *chname, const LIGOTimeGPS *start, REAL8 duration, size_t lengthlimit );
COMPLEX16TimeSeries *XLALFrStreamReadCOMPLEX16TimeSeries( LALFrStream *stream, const char *chname, const LIGOTimeGPS *start, REAL8 duration, size_t lengthlimit );

REAL8TimeSeries * XLALFrStreamInputREAL8TimeSeries( LALFrStream *stream, const char *channel, const LIGOTimeGPS *start, REAL8 duration, size_t lengthlimit );


/*
 *
 * LAL Routines.
 *
 */

void
LALFrCacheOpen(
    LALStatus  *status,
    LALFrStream  **output,
    LALCache    *cache
    );

void
LALFrOpen(
    LALStatus    *status,
    LALFrStream    **stream,
    const CHAR   *dirname,
    const CHAR   *pattern
    );

void
LALFrClose(
    LALStatus  *status,
    LALFrStream  **stream
    );

void
LALFrSetMode(
    LALStatus *status,
    INT4       mode,
    LALFrStream  *stream
    );

void
LALFrEnd(
    LALStatus *status,
    INT4      *end,
    LALFrStream  *stream
    );

void
LALFrNext(
    LALStatus *status,
    LALFrStream  *stream
    );

void
LALFrRewind(
    LALStatus *status,
    LALFrStream  *stream
    );

void
LALFrSeek(
    LALStatus         *status,
    const LIGOTimeGPS *epoch,
    LALFrStream          *stream
    );

void
LALFrTell(
    LALStatus   *status,
    LIGOTimeGPS *epoch,
    LALFrStream    *stream
    );

void
LALFrGetPos(
    LALStatus *status,
    LALFrStreamPos     *position,
    LALFrStream  *stream
    );

void
LALFrSetPos(
    LALStatus *status,
    LALFrStreamPos     *position,
    LALFrStream  *stream
    );

void
LALFrGetTimeSeriesType(
    LALStatus   *status,
    LALTYPECODE *output,
    FrChanIn    *chanin,
    LALFrStream    *stream
    );

void
LALFrGetINT2TimeSeries(
    LALStatus      *status,
    INT2TimeSeries *series,
    FrChanIn       *chanin,
    LALFrStream       *stream
    );

void
LALFrGetINT2TimeSeriesMetadata(
    LALStatus      *status,
    INT2TimeSeries *series,
    FrChanIn       *chanin,
    LALFrStream       *stream
    );

void
LALFrGetINT4TimeSeries(
    LALStatus      *status,
    INT4TimeSeries *series,
    FrChanIn       *chanin,
    LALFrStream       *stream
    );

void
LALFrGetINT4TimeSeriesMetadata(
    LALStatus      *status,
    INT4TimeSeries *series,
    FrChanIn       *chanin,
    LALFrStream       *stream
    );

void
LALFrGetINT8TimeSeries(
    LALStatus      *status,
    INT8TimeSeries *series,
    FrChanIn       *chanin,
    LALFrStream       *stream
    );

void
LALFrGetINT8TimeSeriesMetadata(
    LALStatus      *status,
    INT8TimeSeries *series,
    FrChanIn       *chanin,
    LALFrStream       *stream
    );

void
LALFrGetREAL4TimeSeries(
    LALStatus       *status,
    REAL4TimeSeries *series,
    FrChanIn        *chanin,
    LALFrStream        *stream
    );

void
LALFrGetREAL4TimeSeriesMetadata(
    LALStatus       *status,
    REAL4TimeSeries *series,
    FrChanIn        *chanin,
    LALFrStream        *stream
    );

void
LALFrGetREAL8TimeSeries(
    LALStatus       *status,
    REAL8TimeSeries *series,
    FrChanIn        *chanin,
    LALFrStream        *stream
    );

void
LALFrGetREAL8TimeSeriesMetadata(
    LALStatus       *status,
    REAL8TimeSeries *series,
    FrChanIn        *chanin,
    LALFrStream        *stream
    );

void
LALFrGetCOMPLEX8TimeSeries(
    LALStatus          *status,
    COMPLEX8TimeSeries *series,
    FrChanIn           *chanin,
    LALFrStream           *stream
    );

void
LALFrGetCOMPLEX8TimeSeriesMetadata(
    LALStatus          *status,
    COMPLEX8TimeSeries *series,
    FrChanIn           *chanin,
    LALFrStream           *stream
    );

void
LALFrGetCOMPLEX16TimeSeries(
    LALStatus           *status,
    COMPLEX16TimeSeries *series,
    FrChanIn            *chanin,
    LALFrStream            *stream
    );

void
LALFrGetCOMPLEX16TimeSeriesMetadata(
    LALStatus           *status,
    COMPLEX16TimeSeries *series,
    FrChanIn            *chanin,
    LALFrStream            *stream
    );


void
LALFrGetINT2FrequencySeries(
    LALStatus      *status,
    INT2FrequencySeries *series,
    FrChanIn       *chanin,
    LALFrStream       *stream
    );

void
LALFrGetINT4FrequencySeries(
    LALStatus      *status,
    INT4FrequencySeries *series,
    FrChanIn       *chanin,
    LALFrStream       *stream
    );

void
LALFrGetINT8FrequencySeries(
    LALStatus      *status,
    INT8FrequencySeries *series,
    FrChanIn       *chanin,
    LALFrStream       *stream
    );

void
LALFrGetREAL4FrequencySeries(
    LALStatus       *status,
    REAL4FrequencySeries *series,
    FrChanIn        *chanin,
    LALFrStream        *stream
    );

void
LALFrGetREAL8FrequencySeries(
    LALStatus       *status,
    REAL8FrequencySeries *series,
    FrChanIn        *chanin,
    LALFrStream        *stream
    );

void
LALFrGetCOMPLEX8FrequencySeries(
    LALStatus          *status,
    COMPLEX8FrequencySeries *series,
    FrChanIn           *chanin,
    LALFrStream           *stream
    );

void
LALFrGetCOMPLEX16FrequencySeries(
    LALStatus           *status,
    COMPLEX16FrequencySeries *series,
    FrChanIn            *chanin,
    LALFrStream            *stream
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

/*@}*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALFRSTREAM_H */
