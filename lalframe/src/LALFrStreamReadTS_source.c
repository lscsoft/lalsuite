/*
*  Copyright (C) 2013 Jolien Creighton
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

#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define STYPE CONCAT2(TYPE,TimeSeries)

#define CREATESERIES CONCAT2(XLALCreate,STYPE)
#define DESTROYSERIES CONCAT2(XLALDestroy,STYPE)
#define RESIZESERIES CONCAT2(XLALResize,STYPE)

#define READSERIES CONCAT2(XLALFrFileRead,STYPE)
#define READSERIESMETA CONCAT3(XLALFrFileRead,STYPE,Metadata)
#define STREAMGETSERIES CONCAT2(XLALFrStreamGet,STYPE)
#define STREAMGETSERIESMETA CONCAT3(XLALFrStreamGet,STYPE,Metadata)
#define STREAMREADSERIES CONCAT2(XLALFrStreamRead,STYPE)

int STREAMGETSERIES(STYPE * series, LALFrStream * stream)
{
    const REAL8 fuzz = 0.1 / 16384.0;   /* smallest discernable time */
    const size_t size = sizeof(TYPE);
    size_t noff;
    size_t need;
    size_t ncpy;
    TYPE *dest;
    STYPE *buffer;
    LIGOTimeGPS tend;
    INT8 tnow;
    INT8 tbeg;
    int gap = 0;

    XLAL_CHECK(!(stream->state & LAL_FR_STREAM_END), XLAL_EIO);
    XLAL_CHECK(!(stream->state & LAL_FR_STREAM_ERR), XLAL_EIO);

    /* if series does not have allocation for data,
     * we are to return metadata only, so we don't
     * need to load data in the next call */
    if (series->data && series->data->length)
        buffer = READSERIES(stream->file, series->name, stream->pos);
    else
        buffer = READSERIESMETA(stream->file, series->name, stream->pos);
    if (!buffer)
        XLAL_ERROR(XLAL_EFUNC);

    tnow = XLALGPSToINT8NS(&stream->epoch);
    tbeg = XLALGPSToINT8NS(&buffer->epoch);

    /* Make sure that we aren't requesting data
     * that comes before the current frame.
     * Allow 1 millisecond padding to account
     * for double precision */
    if (tnow + 1000 < tbeg) {
        DESTROYSERIES(buffer);
        XLAL_ERROR(XLAL_ETIME);
    }

    /* compute number of points offset very carefully:
     * if current time is within fuzz of a sample, get
     * that sample; otherwise get the sample just after
     * the requested time */
    noff = ceil((1e-9 * (tnow - tbeg) - fuzz) / buffer->deltaT);

    /* adjust current time to be exactly the first sample
     * (rounded to nearest nanosecond) */
    tnow = tbeg + floor(1e9 * noff * buffer->deltaT + 0.5);
    XLALINT8NSToGPS(&series->epoch, tnow);
    series->deltaT = buffer->deltaT;
    series->sampleUnits = buffer->sampleUnits;

    /* end here if all you want is metadata */
    if (!series->data || !series->data->length) {
        DESTROYSERIES(buffer);
        return 0;
    }

    /* the rest of this function is to get the required
     * amount of data and copy it into the series */

    dest = series->data->data;  /* pointer to where to put the data */
    need = series->data->length;        /* number of points that are needed */

    if (noff > buffer->data->length) {
        /* invalid time offset */
        DESTROYSERIES(buffer);
        XLAL_ERROR(XLAL_ETIME);
    }

    /* copy as much of the buffer is needed */
    ncpy =
        (buffer->data->length - noff) <
        need ? buffer->data->length - noff : need;
    memcpy(dest, buffer->data->data + noff, ncpy * size);
    dest += ncpy;
    need -= ncpy;

    DESTROYSERIES(buffer);

    /* continue while data is required */
    while (need) {

        /* goto next frame */
        if (XLALFrStreamNext(stream) < 0)
            XLAL_ERROR(XLAL_EFUNC);
        if (stream->state & LAL_FR_STREAM_END)
            XLAL_ERROR(XLAL_EIO,
                "End of frame stream while %zd points remain to be read",
                need);

        /* load more data */
        buffer = READSERIES(stream->file, series->name, stream->pos);
        if (!buffer)
            XLAL_ERROR(XLAL_EFUNC);

        if (stream->state & LAL_FR_STREAM_GAP) {
            /* gap in data: reset dest and need and set epoch */
            dest = series->data->data;
            need = series->data->length;
            series->epoch = buffer->epoch;
            gap = 1;
        }

        /* copy data */
        ncpy = buffer->data->length < need ? buffer->data->length : need;
        memcpy(dest, buffer->data->data, ncpy * size);
        dest += ncpy;
        need -= ncpy;
        DESTROYSERIES(buffer);
    }

    /* update stream start time so that it corresponds to the
     * exact time of the next sample to be read */
    stream->epoch = series->epoch;
    XLALGPSAdd(&stream->epoch, series->data->length * series->deltaT);

    /* are we still within the current frame? */
    XLALFrFileQueryGTime(&tend, stream->file, stream->pos);
    XLALGPSAdd(&tend, XLALFrFileQueryDt(stream->file, stream->pos));
    if (XLALGPSCmp(&tend, &stream->epoch) <= 0) {
        /* advance a frame... note that failure here is
         * benign so we suppress gap warnings: these will
         * be triggered on the next read (if one is done) */
        int savemode = stream->mode;
        LIGOTimeGPS saveepoch = stream->epoch;
        stream->mode |= LAL_FR_STREAM_IGNOREGAP_MODE;   /* ignore gaps for now */
        if (XLALFrStreamNext(stream) < 0) {
            stream->mode = savemode;
            XLAL_ERROR(XLAL_EFUNC);
        }
        if (!(stream->state & LAL_FR_STREAM_GAP))       /* no gap: reset epoch */
            stream->epoch = saveepoch;
        stream->mode = savemode;
    }

    /* make sure to set the gap flag in the stream state
     * if a gap had been encountered during the reading */
    if (gap)
        stream->state |= LAL_FR_STREAM_GAP;
    /* FIXME:
     * does this need to cause a failure if mode is set to fail on gaps? */

    /* if the stream state is an error then fail */
    if (stream->state & LAL_FR_STREAM_ERR)
        XLAL_ERROR(XLAL_EIO);

    return 0;
}

int STREAMGETSERIESMETA(STYPE * series, LALFrStream * stream)
{
    STYPE *buffer;
    XLAL_CHECK(!(stream->state & LAL_FR_STREAM_ERR), XLAL_EIO);
    XLAL_CHECK(!(stream->state & LAL_FR_STREAM_END), XLAL_EIO);
    buffer = READSERIESMETA(stream->file, series->name, stream->pos);
    if (!buffer)
        XLAL_ERROR(XLAL_EFUNC);
    series->epoch = buffer->epoch;
    series->deltaT = buffer->deltaT;
    series->sampleUnits = buffer->sampleUnits;
    DESTROYSERIES(buffer);
    return 0;
}

STYPE *STREAMREADSERIES(LALFrStream * stream, const char *chname,
    const LIGOTimeGPS * start, double duration, size_t lengthlimit)
{
    STYPE *series;
    size_t length;

    /* create and initialize a zero-length time series vector */
    series = CREATESERIES(chname, start, 0.0, 0.0, &lalADCCountUnit, 0);
    if (!series)
        XLAL_ERROR_NULL(XLAL_EFUNC);

    /* get the time series meta-data */
    if (STREAMGETSERIES(series, stream)) {
        DESTROYSERIES(series);
        XLAL_ERROR_NULL(XLAL_EFUNC);
    }

    /* resize the time series to the correct number of samples */
    length = duration / series->deltaT;
    if (lengthlimit && (lengthlimit < length))
        length = lengthlimit;
    if (!RESIZESERIES(series, 0, length)) {
        DESTROYSERIES(series);
        XLAL_ERROR_NULL(XLAL_EFUNC);
    }

    /* read the data */
    if (XLALFrStreamSeek(stream, start) || STREAMGETSERIES(series, stream)) {
        DESTROYSERIES(series);
        XLAL_ERROR_NULL(XLAL_EFUNC);
    }

    return series;
}

#undef STYPE
#undef CREATESERIES
#undef DESTROYSERIES
#undef RESIZESERIES
#undef READSERIES
#undef READSERIESMETA
#undef STREAMGETSERIES
#undef STREAMREADSERIES

#undef CONCAT2x
#undef CONCAT2
#undef CONCAT3x
#undef CONCAT3
#undef STRING
