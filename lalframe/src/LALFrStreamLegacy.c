/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Kipp Cannon
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

/* LEGACY CODE */

#include <math.h>
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/LALCache.h>
#include <lal/LALFrStream.h>

void
LALFrCacheOpen(LALStatus * status, LALFrStream ** output, LALCache * cache)
{
    LALFrStream *stream;

    XLALPrintDeprecationWarning(__func__, "XLALFrStreamCacheOpen");
    INITSTATUS(status);
    ASSERT(cache, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    ASSERT(output, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    ASSERT(!*output, status, FRAMESTREAMH_ENNUL, FRAMESTREAMH_MSGENNUL);

    stream = *output = XLALFrStreamCacheOpen(cache);
    if (!stream) {
        int errnum = xlalErrno;
        XLALClearErrno();
        switch (errnum) {
        case XLAL_ENOMEM:
            ABORT(status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC);
        case XLAL_EIO:
            ABORT(status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN);
        default:
            ABORTXLAL(status);
        }
    }

    RETURN(status);
}

void
LALFrOpen(LALStatus * status,
    LALFrStream ** stream, const CHAR * dirname, const CHAR * pattern)
{
    XLALPrintDeprecationWarning(__func__, "XLALFrStreamOpen");
    INITSTATUS(status);
    ATTATCHSTATUSPTR(status);
    ASSERT(stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    ASSERT(!*stream, status, FRAMESTREAMH_ENNUL, FRAMESTREAMH_MSGENNUL);

    *stream = XLALFrStreamOpen(dirname, pattern);

    DETATCHSTATUSPTR(status);
    RETURN(status);
}

void LALFrClose(LALStatus * status, LALFrStream ** stream)
{
    XLALPrintDeprecationWarning(__func__, "XLALFrStreamClose");
    INITSTATUS(status);
    ASSERT(stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    ASSERT(*stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    XLALFrStreamClose(*stream);
    *stream = NULL;
    RETURN(status);
}

void LALFrSetMode(LALStatus * status, INT4 mode, LALFrStream * stream)
{
    XLALPrintDeprecationWarning(__func__, "XLALFrStreamSetMode");
    INITSTATUS(status);
    ASSERT(stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    stream->mode = mode;
    RETURN(status);
}

void LALFrEnd(LALStatus * status, INT4 * end, LALFrStream * stream)
{
    XLALPrintDeprecationWarning(__func__, "XLALFrStreamEnd");
    INITSTATUS(status);
    ASSERT(stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    ASSERT(end, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    *end = XLALFrStreamState(stream) & LAL_FR_STREAM_END;
    RETURN(status);
}

void LALFrRewind(LALStatus * status, LALFrStream * stream)
{
    XLALPrintDeprecationWarning(__func__, "XLALFrStreamRewind");
    INITSTATUS(status);
    ASSERT(stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    if (XLALFrStreamRewind(stream)) {
        XLALClearErrno();
        if (stream->state & LAL_FR_STREAM_URL) {        /* problem was in opening a file */
            ABORT(status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN);
        }
        if (stream->state & LAL_FR_STREAM_TOC) {        /* problem was in reading a file */
            ABORT(status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD);
        }
    }
    RETURN(status);
}

void LALFrNext(LALStatus * status, LALFrStream * stream)
{
    CHAR frErrMsg[1024];
    int code;

    XLALPrintDeprecationWarning(__func__, "XLALFrStreamNext");
    INITSTATUS(status);
    ASSERT(stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);

    if (stream->state & LAL_FR_STREAM_ERR) {
        ABORT(status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR);
    }
    if (stream->state & LAL_FR_STREAM_END) {
        ABORT(status, FRAMESTREAMH_EDONE, FRAMESTREAMH_MSGEDONE);
    }

    code = XLALFrStreamNext(stream);
    if (code < 0) {
        XLALClearErrno();
        if (stream->state & LAL_FR_STREAM_ERR) {
            if (stream->state & LAL_FR_STREAM_URL) {    /* must have failed to open a file */
                snprintf(frErrMsg, sizeof(frErrMsg) / sizeof(*frErrMsg),
                    "Could not open URL %s\n",
                    stream->cache->list[stream->fnum].url);
                LALError(status, frErrMsg);
                ABORT(status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN);
            }
            if (stream->state & LAL_FR_STREAM_TOC) {    /* must have failed to read a file */
                snprintf(frErrMsg, sizeof(frErrMsg) / sizeof(*frErrMsg),
                    "Could not read TOC from %s\n",
                    stream->cache->list[stream->fnum].url);
                LALError(status, frErrMsg);
                ABORT(status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD);
            }
        } else {        /* must be a gap error */

            ABORT(status, FRAMESTREAMH_EDGAP, FRAMESTREAMH_MSGEDGAP);
        }
    }

    RETURN(status);
}

void LALFrSeek(LALStatus * status, const LIGOTimeGPS * epoch,
    LALFrStream * stream)
{
    CHAR frErrMsg[1024];
    int code;

    XLALPrintDeprecationWarning(__func__, "XLALFrStreamSeek");
    INITSTATUS(status);
    ASSERT(stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    ASSERT(epoch, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    if (stream->state & LAL_FR_STREAM_ERR) {
        ABORT(status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR);
    }

    code = XLALFrStreamSeek(stream, epoch);
    if (code < 0) {
        XLALClearErrno();
        if (stream->state & LAL_FR_STREAM_ERR) {        /* a file error */
            if (stream->state & LAL_FR_STREAM_URL) {    /* must have failed to open a file */
                snprintf(frErrMsg, sizeof(frErrMsg) / sizeof(*frErrMsg),
                    "Could not open URL %s\n",
                    stream->cache->list[stream->fnum].url);
                LALError(status, frErrMsg);
                ABORT(status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN);
            }
            if (stream->state & LAL_FR_STREAM_TOC) {    /* must have failed to read a file */
                snprintf(frErrMsg, sizeof(frErrMsg) / sizeof(*frErrMsg),
                    "Could not read TOC from %s\n",
                    stream->cache->list[stream->fnum].url);
                LALError(status, frErrMsg);
                ABORT(status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD);
            }
        } else {        /* must be too early, too late, or in a gap */

            ABORT(status, FRAMESTREAMH_ETREQ, FRAMESTREAMH_MSGETREQ);
        }
    }

    RETURN(status);
}

void LALFrTell(LALStatus * status, LIGOTimeGPS * epoch, LALFrStream * stream)
{
    XLALPrintDeprecationWarning(__func__, "XLALFrStreamTell");
    INITSTATUS(status);
    ASSERT(stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    ASSERT(epoch, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    if (stream->state & LAL_FR_STREAM_ERR) {
        ABORT(status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR);
    }
    XLALFrStreamTell(epoch, stream);
    RETURN(status);
}

void
LALFrGetPos(LALStatus * status, LALFrStreamPos * position,
    LALFrStream * stream)
{
    XLALPrintDeprecationWarning(__func__, "XLALFrStreamGetpos");
    INITSTATUS(status);
    ASSERT(position, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    ASSERT(stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    if (stream->state & LAL_FR_STREAM_ERR) {
        ABORT(status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR);
    }
    XLALFrStreamGetpos(position, stream);
    RETURN(status);
}

void
LALFrSetPos(LALStatus * status, LALFrStreamPos * position,
    LALFrStream * stream)
{
    XLALPrintDeprecationWarning(__func__, "XLALFrStreamSetpos");
    INITSTATUS(status);
    ASSERT(position, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    ASSERT(stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
    if (stream->state & LAL_FR_STREAM_ERR) {
        ABORT(status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR);
    }
    if (XLALFrStreamSetpos(stream, position)) {
        XLALClearErrno();
        if (stream->state & LAL_FR_STREAM_ERR) {
            if (stream->state & LAL_FR_STREAM_URL) {    /* must have failed to open a file */
                ABORT(status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN);
            }
            if (stream->state & LAL_FR_STREAM_TOC) {    /* must have failed to read a file */
                ABORT(status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD);
            }
        }
    }
    RETURN(status);
}

void LALFrGetTimeSeriesType(LALStatus * status, LALTYPECODE * output,
    FrChanIn * chanin, LALFrStream * stream)
{
    XLALPrintDeprecationWarning(__func__, "XLALFrStreamGetTimeSeriesType");
    INITSTATUS(status);
    *output = XLALFrStreamGetTimeSeriesType(chanin->name, stream);
    RETURN(status);
}

/* GET SERIES FUNCTIONS */

#define DEFINE_LAL_GET_TS_FUNCTION(laltype) \
    void LALFrGet ## laltype ## TimeSeries(LALStatus *status, laltype ## TimeSeries *series, FrChanIn *chanin, LALFrStream *stream) \
    { \
        int errnum; \
        int code; \
        XLALPrintDeprecationWarning(__func__, "XLALFrStreamGet" #laltype "TimeSeries"); \
        INITSTATUS(status); \
        strcpy(series->name, chanin->name); \
        XLAL_TRY(code = XLALFrStreamGet ## laltype ## TimeSeries(series, stream), errnum); \
        if ((code < 0) || errnum) { \
            if (stream->state & LAL_FR_STREAM_END) { \
                ABORT(status, FRAMESTREAMH_EDONE, FRAMESTREAMH_MSGEDONE); \
            } \
            ABORT(status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD); \
        } \
        RETURN(status); \
    }

#define DEFINE_LAL_GET_TSM_FUNCTION(laltype) \
    void LALFrGet ## laltype ## TimeSeriesMetadata(LALStatus *status, laltype ## TimeSeries *series, FrChanIn *chanin, LALFrStream *stream) \
    { \
        int code; \
        XLALPrintDeprecationWarning(__func__, "XLALFrStreamGet" #laltype "TimeSeriesMetadata"); \
        INITSTATUS(status); \
        strcpy(series->name, chanin->name); \
        code = XLALFrStreamGet ## laltype ## TimeSeriesMetadata(series, stream); \
        if (code < 0) { \
            ABORT(status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD); \
        } \
        RETURN(status); \
    }

#define DEFINE_LAL_GET_FS_FUNCTION(laltype) \
    void LALFrGet ## laltype ## FrequencySeries(LALStatus *status, laltype ## FrequencySeries *series, FrChanIn *chanin, LALFrStream *stream) \
    { \
        int code; \
        XLALPrintDeprecationWarning(__func__, "XLALFrStreamGet" #laltype "FrequencySeries"); \
        INITSTATUS(status); \
        strcpy(series->name, chanin->name); \
        code = XLALFrStreamGet ## laltype ## FrequencySeries(series, stream); \
        if (code < 0) { \
            ABORT(status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD); \
        } \
        RETURN(status); \
    }

/* *INDENT-OFF* */
DEFINE_LAL_GET_TS_FUNCTION(INT2)
DEFINE_LAL_GET_TS_FUNCTION(INT4)
DEFINE_LAL_GET_TS_FUNCTION(INT8)
DEFINE_LAL_GET_TS_FUNCTION(REAL4)
DEFINE_LAL_GET_TS_FUNCTION(REAL8)
DEFINE_LAL_GET_TS_FUNCTION(COMPLEX8)
DEFINE_LAL_GET_TS_FUNCTION(COMPLEX16)

DEFINE_LAL_GET_TSM_FUNCTION(INT2)
DEFINE_LAL_GET_TSM_FUNCTION(INT4)
DEFINE_LAL_GET_TSM_FUNCTION(INT8)
DEFINE_LAL_GET_TSM_FUNCTION(REAL4)
DEFINE_LAL_GET_TSM_FUNCTION(REAL8)
DEFINE_LAL_GET_TSM_FUNCTION(COMPLEX8)
DEFINE_LAL_GET_TSM_FUNCTION(COMPLEX16)

DEFINE_LAL_GET_FS_FUNCTION(REAL4)
DEFINE_LAL_GET_FS_FUNCTION(REAL8)
DEFINE_LAL_GET_FS_FUNCTION(COMPLEX8)
DEFINE_LAL_GET_FS_FUNCTION(COMPLEX16)
/* *INDENT-ON* */

/* WRITE SERIES FUNCTIONS */

/* FIXME: now only supports one frame per file and ProcData channels only */
#define DEFINE_LAL_WRITE_TS_FUNCTION(laltype) \
    void LALFrWrite ## laltype ## TimeSeries(LALStatus *status, laltype ## TimeSeries *series, FrOutPar *params) \
    { \
        LALFrameH *frame; \
        char fname[FILENAME_MAX]; \
        double duration; \
        int t0, dt; \
        XLALPrintDeprecationWarning(__func__, "XLALFrWrite" #laltype "TimeSeries"); \
        INITSTATUS(status); \
        duration = series->deltaT * series->data->length; \
        t0 = series->epoch.gpsSeconds; \
        dt = (int)ceil(XLALGPSGetREAL8(&series->epoch)+duration) - t0; \
        snprintf(fname, sizeof(fname), "%s-%s-%d-%d.gwf", \
             params->source ? params->source : "F", \
             params->description ? params->description : "UNKNOWN", \
             t0, dt); \
        frame = XLALFrameNew(&series->epoch, duration, "LAL", params->run, params->frame, 0); \
        XLALFrameAdd ## laltype ## TimeSeriesProcData(frame, series); \
        XLALFrameWrite(frame, fname); \
        XLALFrameFree(frame); \
        RETURN(status); \
    }

#define DEFINE_LAL_WRITE_FS_FUNCTION(laltype) \
    void LALFrWrite ## laltype ## FrequencySeries(LALStatus *status, laltype ## FrequencySeries *series, FrOutPar *params, int subtype) \
    { \
        LALFrameH *frame; \
        char fname[FILENAME_MAX]; \
        double duration; \
        int t0, dt; \
        XLALPrintDeprecationWarning(__func__, "XLALFrWrite" #laltype "FrequencySeries"); \
        INITSTATUS(status); \
        duration = series->deltaF ? 1.0 / series->deltaF : 1.0; \
        t0 = series->epoch.gpsSeconds; \
        dt = (int)ceil(XLALGPSGetREAL8(&series->epoch)+duration) - t0; \
        dt = dt < 1 ? 1 : dt; \
        snprintf(fname, sizeof(fname), "%s-%s-%d-%d.gwf", \
             params->source ? params->source : "F", \
             params->description ? params->description : "UNKNOWN", \
             t0, dt); \
        frame = XLALFrameNew(&series->epoch, duration, "LAL", params->run, params->frame, 0); \
        XLALFrameAdd ## laltype ## FrequencySeriesProcData(frame, series, subtype); \
        XLALFrameWrite(frame, fname); \
        XLALFrameFree(frame); \
        RETURN(status); \
    }

/* *INDENT-OFF* */
DEFINE_LAL_WRITE_TS_FUNCTION(INT2)
DEFINE_LAL_WRITE_TS_FUNCTION(INT4)
DEFINE_LAL_WRITE_TS_FUNCTION(INT8)
DEFINE_LAL_WRITE_TS_FUNCTION(REAL4)
DEFINE_LAL_WRITE_TS_FUNCTION(REAL8)
DEFINE_LAL_WRITE_TS_FUNCTION(COMPLEX8)
DEFINE_LAL_WRITE_TS_FUNCTION(COMPLEX16)

DEFINE_LAL_WRITE_FS_FUNCTION(REAL4)
DEFINE_LAL_WRITE_FS_FUNCTION(REAL8)
DEFINE_LAL_WRITE_FS_FUNCTION(COMPLEX8)
DEFINE_LAL_WRITE_FS_FUNCTION(COMPLEX16)
/* *INDENT-ON* */
