#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/LALFrameIO.h>
#include <lal/LALFrStream.h>

int XLALFrStreamGetVectorLength(const char *chname, LALFrStream * stream)
{
    return XLALFrFileQueryChanVectorLength(stream->file, chname, stream->pos);
}

LALTYPECODE XLALFrStreamGetTimeSeriesType(const char *chname,
    LALFrStream * stream)
{
    return XLALFrFileQueryChanType(stream->file, chname, stream->pos);
}

#define TYPE INT2
#include "LALFrStreamReadTS_source.c"
#undef TYPE

#define TYPE INT4
#include "LALFrStreamReadTS_source.c"
#undef TYPE

#define TYPE INT8
#include "LALFrStreamReadTS_source.c"
#undef TYPE

#define TYPE UINT2
#include "LALFrStreamReadTS_source.c"
#undef TYPE

#define TYPE UINT4
#include "LALFrStreamReadTS_source.c"
#undef TYPE

#define TYPE UINT8
#include "LALFrStreamReadTS_source.c"
#undef TYPE

#define TYPE REAL4
#include "LALFrStreamReadTS_source.c"
#include "LALFrStreamReadFS_source.c"
#undef TYPE

#define TYPE REAL8
#include "LALFrStreamReadTS_source.c"
#include "LALFrStreamReadFS_source.c"
#undef TYPE

#define TYPE COMPLEX8
#include "LALFrStreamReadTS_source.c"
#include "LALFrStreamReadFS_source.c"
#undef TYPE

#define TYPE COMPLEX16
#include "LALFrStreamReadTS_source.c"
#include "LALFrStreamReadFS_source.c"
#undef TYPE

/* scalar-to-scalar copy */
#define COPY_S2S(dest, orig, n) \
    do { \
        size_t i; \
        for (i = 0; i < n; ++i) (dest)[i] = (orig)[i]; \
    } while (0)

/* scalar-to-complex copy */
#define COPY_S2C(dest, orig, n) \
    do { \
        size_t i; \
        for (i = 0; i < n; ++i) (dest)[i].re = (orig)[i]; \
    } while (0)

/* complex-to-complex copy */
#define COPY_C2C(dest, orig, n) \
    do { \
        size_t i; \
        for (i = 0; i < n; ++i) { \
            (dest)[i].re = (orig)[i].re; \
            (dest)[i].im = (orig)[i].im; \
        } \
    } while (0)

#define INPUTTS(series, desttype, origtype, promotion, stream, chname, start, duration, lengthlimit) \
    do { \
        origtype ## TimeSeries *origin; \
        origin = XLALFrStreamRead##origtype##TimeSeries((stream),(chname),(start),(duration),(lengthlimit)); \
        series = XLALCreate##desttype##TimeSeries((chname),(start),origin->f0,origin->deltaT,&origin->sampleUnits,origin->data->length); \
        if (!origin || !series) { \
            XLALDestroy##origtype##TimeSeries(origin); \
            XLALDestroy##desttype##TimeSeries(series); \
            XLAL_ERROR_NULL(XLAL_EFUNC); \
        } \
        COPY_##promotion(series->data->data, origin->data->data, origin->data->length); \
        XLALDestroy##origtype##TimeSeries(origin); \
    } while(0)

#define INPUTFS(series, desttype, origtype, promotion, stream, chname, epoch) \
    do { \
        origtype ## FrequencySeries *origin; \
        origin = XLALFrStreamRead##origtype##FrequencySeries((stream),(chname),(epoch)); \
        series = XLALCreate##desttype##FrequencySeries((chname),(epoch),origin->f0,origin->deltaF,&origin->sampleUnits,origin->data->length); \
        if (!origin || !series) { \
            XLALDestroy##origtype##FrequencySeries(origin); \
            XLALDestroy##desttype##FrequencySeries(series); \
            XLAL_ERROR_NULL(XLAL_EFUNC); \
        } \
        COPY_##promotion(series->data->data, origin->data->data, origin->data->length); \
        XLALDestroy##origtype##FrequencySeries(origin); \
    } while(0)

REAL8TimeSeries *XLALFrStreamInputREAL8TimeSeries(LALFrStream * stream,
    const char *chname, const LIGOTimeGPS * start, double duration,
    size_t lengthlimit)
{
    REAL8TimeSeries *series;
    LALTYPECODE typecode;

    if (XLALFrStreamSeek(stream, start))
        XLAL_ERROR_NULL(XLAL_EFUNC);

    typecode = XLALFrFileQueryChanType(stream->file, chname, stream->pos);
    switch (typecode) {
    case LAL_I2_TYPE_CODE:
        INPUTTS(series, REAL8, INT2, S2S, stream, chname, start, duration,
            lengthlimit);
        break;
    case LAL_I4_TYPE_CODE:
        INPUTTS(series, REAL8, INT4, S2S, stream, chname, start, duration,
            lengthlimit);
        break;
    case LAL_I8_TYPE_CODE:
        INPUTTS(series, REAL8, INT8, S2S, stream, chname, start, duration,
            lengthlimit);
        break;
    case LAL_U2_TYPE_CODE:
        INPUTTS(series, REAL8, UINT2, S2S, stream, chname, start, duration,
            lengthlimit);
        break;
    case LAL_U4_TYPE_CODE:
        INPUTTS(series, REAL8, UINT4, S2S, stream, chname, start, duration,
            lengthlimit);
        break;
    case LAL_U8_TYPE_CODE:
        INPUTTS(series, REAL8, UINT8, S2S, stream, chname, start, duration,
            lengthlimit);
        break;
    case LAL_S_TYPE_CODE:
        INPUTTS(series, REAL8, REAL4, S2S, stream, chname, start, duration,
            lengthlimit);
        break;
    case LAL_D_TYPE_CODE:
        series =
            XLALFrStreamReadREAL8TimeSeries(stream, chname, start,
            duration, lengthlimit);
        if (!series)
            XLAL_ERROR_NULL(XLAL_EFUNC);
        break;
    case LAL_C_TYPE_CODE:
    case LAL_Z_TYPE_CODE:
        XLAL_PRINT_ERROR("Cannot convert complex type to float type");
    default:
        XLAL_ERROR_NULL(XLAL_ETYPE);
    }
    return series;
}

COMPLEX16TimeSeries *XLALFrStreamInputCOMPLEX16TimeSeries(LALFrStream *
    stream, const char *chname, const LIGOTimeGPS * start, double duration,
    size_t lengthlimit)
{
    COMPLEX16TimeSeries *series;
    LALTYPECODE typecode;

    if (XLALFrStreamSeek(stream, start))
        XLAL_ERROR_NULL(XLAL_EFUNC);

    typecode = XLALFrFileQueryChanType(stream->file, chname, stream->pos);
    switch (typecode) {
    case LAL_I2_TYPE_CODE:
        INPUTTS(series, COMPLEX16, INT2, S2S, stream, chname, start,
            duration, lengthlimit);
        break;
    case LAL_I4_TYPE_CODE:
        INPUTTS(series, COMPLEX16, INT4, S2S, stream, chname, start,
            duration, lengthlimit);
        break;
    case LAL_I8_TYPE_CODE:
        INPUTTS(series, COMPLEX16, INT8, S2S, stream, chname, start,
            duration, lengthlimit);
        break;
    case LAL_U2_TYPE_CODE:
        INPUTTS(series, COMPLEX16, UINT2, S2S, stream, chname, start,
            duration, lengthlimit);
        break;
    case LAL_U4_TYPE_CODE:
        INPUTTS(series, COMPLEX16, UINT4, S2S, stream, chname, start,
            duration, lengthlimit);
        break;
    case LAL_U8_TYPE_CODE:
        INPUTTS(series, COMPLEX16, UINT8, S2S, stream, chname, start,
            duration, lengthlimit);
        break;
    case LAL_S_TYPE_CODE:
        INPUTTS(series, COMPLEX16, REAL4, S2S, stream, chname, start,
            duration, lengthlimit);
        break;
    case LAL_D_TYPE_CODE:
        INPUTTS(series, COMPLEX16, REAL8, S2S, stream, chname, start,
            duration, lengthlimit);
        break;
    case LAL_C_TYPE_CODE:
        INPUTTS(series, COMPLEX16, COMPLEX8, S2S, stream, chname, start,
            duration, lengthlimit);
    case LAL_Z_TYPE_CODE:
        series = XLALFrStreamReadCOMPLEX16TimeSeries(stream, chname, start,
            duration, lengthlimit);
        if (!series)
            XLAL_ERROR_NULL(XLAL_EFUNC);
    default:
        XLAL_ERROR_NULL(XLAL_ETYPE);
    }
    return series;
}

REAL8FrequencySeries *XLALFrStreamInputREAL8FrequencySeries(LALFrStream *
    stream, const char *chname, const LIGOTimeGPS * epoch)
{
    REAL8FrequencySeries *series;
    LALTYPECODE typecode;

    if (XLALFrStreamSeek(stream, epoch))
        XLAL_ERROR_NULL(XLAL_EFUNC);

    typecode = XLALFrFileQueryChanType(stream->file, chname, stream->pos);
    switch (typecode) {
    case LAL_S_TYPE_CODE:
        INPUTFS(series, REAL8, REAL4, S2S, stream, chname, epoch);
        break;
    case LAL_D_TYPE_CODE:
        series = XLALFrStreamReadREAL8FrequencySeries(stream, chname, epoch);
        if (!series)
            XLAL_ERROR_NULL(XLAL_EFUNC);
        break;
    case LAL_C_TYPE_CODE:
    case LAL_Z_TYPE_CODE:
        XLAL_PRINT_ERROR("Cannot convert complex type to float type");
    default:
        XLAL_ERROR_NULL(XLAL_ETYPE);
    }
    return series;
}

COMPLEX16FrequencySeries
    *XLALFrStreamInputCOMPLEX16FrequencySeries(LALFrStream * stream,
    const char *chname, const LIGOTimeGPS * epoch)
{
    COMPLEX16FrequencySeries *series;
    LALTYPECODE typecode;

    if (XLALFrStreamSeek(stream, epoch))
        XLAL_ERROR_NULL(XLAL_EFUNC);

    typecode = XLALFrFileQueryChanType(stream->file, chname, stream->pos);
    switch (typecode) {
    case LAL_S_TYPE_CODE:
        INPUTFS(series, COMPLEX16, REAL4, S2S, stream, chname, epoch);
        break;
    case LAL_D_TYPE_CODE:
        INPUTFS(series, COMPLEX16, REAL8, S2S, stream, chname, epoch);
        break;
    case LAL_C_TYPE_CODE:
        INPUTFS(series, COMPLEX16, COMPLEX8, S2S, stream, chname, epoch);
    case LAL_Z_TYPE_CODE:
        series =
            XLALFrStreamReadCOMPLEX16FrequencySeries(stream, chname, epoch);
        if (!series)
            XLAL_ERROR_NULL(XLAL_EFUNC);
    default:
        XLAL_ERROR_NULL(XLAL_ETYPE);
    }
    return series;
}
