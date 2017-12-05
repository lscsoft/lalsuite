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
/**
 * @addtogroup LALFrStreamRead_c
 * @brief Provides routines for reading data from a #LALFrStream.
 * @details
 * The routines described below perform either random-access or sequential
 * reading of time-series or frequency-series data from #LALFrStream streams.
 * @{
 */

#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/LALFrameIO.h>
#include <lal/LALFrStream.h>

/**
 * @brief Returns the number of data points in channel @p chname in the
 * current frame in frame stream @p stream.
 * @param chname String containing the name of the channel.
 * @param stream Pointer to a #LALFrStream structure.
 * @returns The length of the data vector of the channel in the current frame.
 * @retval -1 Failure.
 */
int XLALFrStreamGetVectorLength(const char *chname, LALFrStream * stream)
{
    return XLALFrFileQueryChanVectorLength(stream->file, chname, stream->pos);
}

/** 
 * @brief Returns the type code for the data type of channel @p chname in the
 * current frame in frame stream @p stream.
 * @param chname String containing the name of the channel.
 * @param stream Pointer to a #LALFrStream structure.
 * @returns The #LALTYPECODE value of the data type of the channel.
 * @retval LAL_CHAR_TYPE_CODE Channel is an array of type char.
 * @retval LAL_I2_TYPE_CODE Channel is an array of type int16_t.
 * @retval LAL_I4_TYPE_CODE Channel is an array of type int32_t.
 * @retval LAL_I8_TYPE_CODE Channel is an array of type int64_t.
 * @retval LAL_UCHAR_TYPE_CODE Channel is an array of type unsigned char.
 * @retval LAL_U2_TYPE_CODE Channel is an array of type uint16_t.
 * @retval LAL_U4_TYPE_CODE Channel is an array of type uint32_t.
 * @retval LAL_U8_TYPE_CODE Channel is an array of type uint64_t.
 * @retval LAL_S_TYPE_CODE Channel is an array of type float.
 * @retval LAL_D_TYPE_CODE Channel is an array of type double.
 * @retval LAL_C_TYPE_CODE Channel is an array of type float complex.
 * @retval LAL_Z_TYPE_CODE Channel is an array of type double complex.
 * @retval -1 Failure.
 */
LALTYPECODE XLALFrStreamGetTimeSeriesType(const char *chname, LALFrStream * stream)
{
    return XLALFrFileQueryChanType(stream->file, chname, stream->pos);
}

/** @cond */

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

/** @endcond */


/**
 * @name Random-Access Time Series Reading Routines
 *
 * These routines allow the user to read channel data for a specified
 * period of time.  The following is an example of how to read 16 seconds of
 * REAL8 data from channel `STRAIN` beginning at GPS time 123456789 from
 * frame files in the current directory.
 *
 * @code
 * #include <lal/LALFrStream.h>
 * ...
 * const char *dirname = ".";
 * const char *pattern = "*.gwf";
 * const char *chname = "STRAIN";
 * LALFrStream *stream;
 * REAL8TimeSeries *series;
 * LIGOTimeGPS start = { 123456789 };
 * REAL8 duration = 16.0;
 * ...
 * stream = XLALFrStreamOpen(dirname, pattern);
 * series = XLALFrStreamInputREAL8TimeSeries(stream, chname, &start, duration, 0);
 * @endcode
 *
 * Note that the XLALFrStreamInputREAL8TimeSeries() will convert channel data
 * of a different datatype to REAL8 data (so, for example, if the data stored
 * in the frame file is REAL4 data, then it is cast as REAL8 data).  If this
 * behaviour is not wanted, use XLALFrStreamReadREAL8TimeSeries() instead.
 *
 * @{
 */

/**
 * @fn INT2TimeSeries *XLALFrStreamReadINT2TimeSeries(LALFrStream * stream, const char *chname, const LIGOTimeGPS * start, double duration, size_t lengthlimit)
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified start time and duration.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * beginning at a specified point in time and lasting a specified duration.
 * If there is a gap in the data, this routine skips to the next contiguous
 * set of data of the required duration.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param start Pointer to a LIGOTimeGPS structure specifying the start time.
 * @param duration The duration of the data to read, in seconds.
 * @param lengthlimit The maximum number of points to read or 0 for unlimited.
 * @returns Pointer to a new INT2TimeSeries containing the specified data.
 * @retval NULL Failure.
 */

/**
 * @fn INT4TimeSeries *XLALFrStreamReadINT4TimeSeries(LALFrStream * stream, const char *chname, const LIGOTimeGPS * start, double duration, size_t lengthlimit)
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified start time and duration.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * beginning at a specified point in time and lasting a specified duration.
 * If there is a gap in the data, this routine skips to the next contiguous
 * set of data of the required duration.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param start Pointer to a LIGOTimeGPS structure specifying the start time.
 * @param duration The duration of the data to read, in seconds.
 * @param lengthlimit The maximum number of points to read or 0 for unlimited.
 * @returns Pointer to a new INT4TimeSeries containing the specified data.
 * @retval NULL Failure.
 */

/**
 * @fn INT8TimeSeries *XLALFrStreamReadINT8TimeSeries(LALFrStream * stream, const char *chname, const LIGOTimeGPS * start, double duration, size_t lengthlimit)
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified start time and duration.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * beginning at a specified point in time and lasting a specified duration.
 * If there is a gap in the data, this routine skips to the next contiguous
 * set of data of the required duration.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param start Pointer to a LIGOTimeGPS structure specifying the start time.
 * @param duration The duration of the data to read, in seconds.
 * @param lengthlimit The maximum number of points to read or 0 for unlimited.
 * @returns Pointer to a new INT8TimeSeries containing the specified data.
 * @retval NULL Failure.
 */

/**
 * @fn REAL8TimeSeries *XLALFrStreamReadREAL8TimeSeries(LALFrStream * stream, const char *chname, const LIGOTimeGPS * start, double duration, size_t lengthlimit)
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified start time and duration.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * beginning at a specified point in time and lasting a specified duration.
 * If there is a gap in the data, this routine skips to the next contiguous
 * set of data of the required duration.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param start Pointer to a LIGOTimeGPS structure specifying the start time.
 * @param duration The duration of the data to read, in seconds.
 * @param lengthlimit The maximum number of points to read or 0 for unlimited.
 * @returns Pointer to a new REAL8TimeSeries containing the specified data.
 * @retval NULL Failure.
 */

/**
 * @fn UINT2TimeSeries *XLALFrStreamReadUINT2TimeSeries(LALFrStream * stream, const char *chname, const LIGOTimeGPS * start, double duration, size_t lengthlimit)
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified start time and duration.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * beginning at a specified point in time and lasting a specified duration.
 * If there is a gap in the data, this routine skips to the next contiguous
 * set of data of the required duration.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param start Pointer to a LIGOTimeGPS structure specifying the start time.
 * @param duration The duration of the data to read, in seconds.
 * @param lengthlimit The maximum number of points to read or 0 for unlimited.
 * @returns Pointer to a new UINT2TimeSeries containing the specified data.
 * @retval NULL Failure.
 */

/**
 * @fn UINT4TimeSeries *XLALFrStreamReadUINT4TimeSeries(LALFrStream * stream, const char *chname, const LIGOTimeGPS * start, double duration, size_t lengthlimit)
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified start time and duration.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * beginning at a specified point in time and lasting a specified duration.
 * If there is a gap in the data, this routine skips to the next contiguous
 * set of data of the required duration.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param start Pointer to a LIGOTimeGPS structure specifying the start time.
 * @param duration The duration of the data to read, in seconds.
 * @param lengthlimit The maximum number of points to read or 0 for unlimited.
 * @returns Pointer to a new UINT4TimeSeries containing the specified data.
 * @retval NULL Failure.
 */

/**
 * @fn UINT8TimeSeries *XLALFrStreamReadUINT8TimeSeries(LALFrStream * stream, const char *chname, const LIGOTimeGPS * start, double duration, size_t lengthlimit)
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified start time and duration.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * beginning at a specified point in time and lasting a specified duration.
 * If there is a gap in the data, this routine skips to the next contiguous
 * set of data of the required duration.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param start Pointer to a LIGOTimeGPS structure specifying the start time.
 * @param duration The duration of the data to read, in seconds.
 * @param lengthlimit The maximum number of points to read or 0 for unlimited.
 * @returns Pointer to a new UINT8TimeSeries containing the specified data.
 * @retval NULL Failure.
 */

/**
 * @fn REAL4TimeSeries *XLALFrStreamReadREAL4TimeSeries(LALFrStream * stream, const char *chname, const LIGOTimeGPS * start, double duration, size_t lengthlimit)
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified start time and duration.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * beginning at a specified point in time and lasting a specified duration.
 * If there is a gap in the data, this routine skips to the next contiguous
 * set of data of the required duration.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param start Pointer to a LIGOTimeGPS structure specifying the start time.
 * @param duration The duration of the data to read, in seconds.
 * @param lengthlimit The maximum number of points to read or 0 for unlimited.
 * @returns Pointer to a new REAL4TimeSeries containing the specified data.
 * @retval NULL Failure.
 */

/**
 * @fn REAL8TimeSeries *XLALFrStreamReadREAL8TimeSeries(LALFrStream * stream, const char *chname, const LIGOTimeGPS * start, double duration, size_t lengthlimit)
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified start time and duration.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * beginning at a specified point in time and lasting a specified duration.
 * If there is a gap in the data, this routine skips to the next contiguous
 * set of data of the required duration.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param start Pointer to a LIGOTimeGPS structure specifying the start time.
 * @param duration The duration of the data to read, in seconds.
 * @param lengthlimit The maximum number of points to read or 0 for unlimited.
 * @returns Pointer to a new REAL8TimeSeries containing the specified data.
 * @retval NULL Failure.
 */

/**
 * @fn COMPLEX8TimeSeries *XLALFrStreamReadCOMPLEX8TimeSeries(LALFrStream * stream, const char *chname, const LIGOTimeGPS * start, double duration, size_t lengthlimit)
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified start time and duration.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * beginning at a specified point in time and lasting a specified duration.
 * If there is a gap in the data, this routine skips to the next contiguous
 * set of data of the required duration.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param start Pointer to a LIGOTimeGPS structure specifying the start time.
 * @param duration The duration of the data to read, in seconds.
 * @param lengthlimit The maximum number of points to read or 0 for unlimited.
 * @returns Pointer to a new COMPLEX8TimeSeries containing the specified data.
 * @retval NULL Failure.
 */

/**
 * @fn COMPLEX16TimeSeries *XLALFrStreamReadCOMPLEX16TimeSeries(LALFrStream * stream, const char *chname, const LIGOTimeGPS * start, double duration, size_t lengthlimit)
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified start time and duration.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * beginning at a specified point in time and lasting a specified duration.
 * If there is a gap in the data, this routine skips to the next contiguous
 * set of data of the required duration.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param start Pointer to a LIGOTimeGPS structure specifying the start time.
 * @param duration The duration of the data to read, in seconds.
 * @param lengthlimit The maximum number of points to read or 0 for unlimited.
 * @returns Pointer to a new COMPLEX16TimeSeries containing the specified data.
 * @retval NULL Failure.
 */

/**
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified start time and duration, and performs any needed type conversion.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * beginning at a specified point in time and lasting a specified duration.
 * If there is a gap in the data, this routine skips to the next contiguous
 * set of data of the required duration.  If the channel being read is not
 * REAL8, the data is converted to type REAL8.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param start Pointer to a LIGOTimeGPS structure specifying the start time.
 * @param duration The duration of the data to read, in seconds.
 * @param lengthlimit The maximum number of points to read or 0 for unlimited.
 * @returns Pointer to a new REAL8TimeSeries containing the specified data.
 * @retval NULL Failure.
 */
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
#if __GNUC__ >= 7
	__attribute__ ((fallthrough));
#endif
    default:
        XLAL_ERROR_NULL(XLAL_ETYPE);
    }
    return series;
}

/**
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified start time and duration, and performs any needed type conversion.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * beginning at a specified point in time and lasting a specified duration.
 * If there is a gap in the data, this routine skips to the next contiguous
 * set of data of the required duration.  If the channel being read is not
 * COMPLEX16, the data is converted to type COMPLEX16.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param start Pointer to a LIGOTimeGPS structure specifying the start time.
 * @param duration The duration of the data to read, in seconds.
 * @param lengthlimit The maximum number of points to read or 0 for unlimited.
 * @returns Pointer to a new COMPLEX16TimeSeries containing the specified data.
 * @retval NULL Failure.
 */
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
        break;
    case LAL_Z_TYPE_CODE:
        series = XLALFrStreamReadCOMPLEX16TimeSeries(stream, chname, start,
            duration, lengthlimit);
        if (!series)
            XLAL_ERROR_NULL(XLAL_EFUNC);
        break;
    default:
        XLAL_ERROR_NULL(XLAL_ETYPE);
    }
    return series;
}

/** @} */

/**
 * @name Random-Access Frequency Series Reading Routines
 *
 * These routines allow the user to read channel data for a specified
 * period of time.  The following is an example of how to read a frequency
 * series REAL8 data from channel `STRAIN_PSD` beginning at GPS time 123456789
 * from frame files in the current directory.
 *
 * @code
 * #include <lal/LALFrStream.h>
 * ...
 * const char *dirname = ".";
 * const char *pattern = "*.gwf";
 * const char *chname = "STRAIN_PSD";
 * LALFrStream *stream;
 * REAL8FrequencySeries *series;
 * LIGOTimeGPS start = { 123456789 };
 * ...
 * stream = XLALFrStreamOpen(dirname, pattern);
 * series = XLALFrStreamInputREAL8FrequencySeries(stream, chname, &start);
 * ...
 * @endcode
 *
 * Note that the XLALFrStreamInputREAL8FrequencySeries() will convert channel
 * data of a different datatype to REAL8 data (so, for example, if the data
 * stored in the frame file is REAL4 data, then it is cast as REAL8 data).  If
 * this behaviour is not wanted, use XLALFrStreamReadREAL8Frequencyeries()
 * instead.
 * @{
 */

/**
 * @fn REAL4FrequencySeries *XLALFrStreamReadREAL4FrequencySeries(LALFrStream * stream, const char *chname, const LIGOTimeGPS * epoch)
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified epoch.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * contained in a frame spanning the specified epoch.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param epoch Pointer to a LIGOTimeGPS structure specifying epoch.
 * @returns Pointer to a new REAL4FrequencySeries containing the specified
 * data.
 * @retval NULL Failure.
 */

/**
 * @fn REAL8FrequencySeries *XLALFrStreamReadREAL8FrequencySeries(LALFrStream * stream, const char *chname, const LIGOTimeGPS * epoch)
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified epoch.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * contained in a frame spanning the specified epoch.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param epoch Pointer to a LIGOTimeGPS structure specifying epoch.
 * @returns Pointer to a new REAL8FrequencySeries containing the specified
 * data.
 * @retval NULL Failure.
 */

/**
 * @fn COMPLEX8FrequencySeries *XLALFrStreamReadCOMPLEX8FrequencySeries(LALFrStream * stream, const char *chname, const LIGOTimeGPS * epoch)
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified epoch.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * contained in a frame spanning the specified epoch.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param epoch Pointer to a LIGOTimeGPS structure specifying epoch.
 * @returns Pointer to a new COMPLEX8FrequencySeries containing the specified
 * data.
 * @retval NULL Failure.
 */

/**
 * @fn COMPLEX16FrequencySeries *XLALFrStreamReadCOMPLEX16FrequencySeries(LALFrStream * stream, const char *chname, const LIGOTimeGPS * epoch)
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified epoch.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * contained in a frame spanning the specified epoch.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param epoch Pointer to a LIGOTimeGPS structure specifying epoch.
 * @returns Pointer to a new COMPLEX16FrequencySeries containing the specified
 * data.
 * @retval NULL Failure.
 */

/**
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified epoch, and performs any needed type conversion.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * contained in a frame spanning the specified epoch.  If the channel being
 * read is not REAL8, the data is converted to type REAL8.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param epoch Pointer to a LIGOTimeGPS structure specifying epoch.
 * @returns Pointer to a new REAL8FrequencySeries containing the specified
 * data.
 * @retval NULL Failure.
 */
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
#if __GNUC__ >= 7
	__attribute__ ((fallthrough));
#endif
    default:
        XLAL_ERROR_NULL(XLAL_ETYPE);
    }
    return series;
}

/**
 * @brief Reads a time series channel from a #LALFrStream stream with a
 * specified epoch, and performs any needed type conversion.
 * @details
 * This routine reads data from a specified channel from a #LALFrStream 
 * contained in a frame spanning the specified epoch.  If the channel being
 * read is not COMPLEX16, the data is converted to type COMPLEX16.
 * @param stream Pointer to the #LALFrStream stream.
 * @param chname String with the channel name to read.
 * @param epoch Pointer to a LIGOTimeGPS structure specifying epoch.
 * @returns Pointer to a new COMPLEX16FrequencySeries containing the specified
 * data.
 * @retval NULL Failure.
 */
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
        break;
    case LAL_Z_TYPE_CODE:
        series =
            XLALFrStreamReadCOMPLEX16FrequencySeries(stream, chname, epoch);
        if (!series)
            XLAL_ERROR_NULL(XLAL_EFUNC);
        break;
    default:
        XLAL_ERROR_NULL(XLAL_ETYPE);
    }
    return series;
}

/** @} */

/**
 * @name Sequential Time Series Reading Routines
 *
 * These routines are useful for reading sequential blocks of data from
 * a #LALFrStream stream.  The following is an example of how to read
 * REAL8 data from channel `STRAIN` in 1 second blocks starting at the
 * beginning of a stream from frame files in the current directory.
 *
 * @code
 * #include <lal/TimeSeries.h>
 * #include <lal/LALFrStream.h>
 * ...
 * const char *dirname = ".";
 * const char *pattern = "*.gwf";
 * const char *chname = "STRAIN";
 * LIGOTimeGPS epoch = LIGOTIMEGPSZERO; // value is unimportant
 * LALFrStream *stream;
 * REAL8TimeSeries *series;
 * REAL8 duration = 1.0;
 * ...
 * stream = XLALFrStreamOpen(dirname, pattern);
 *
 * // initial values of metadata other than the name are not important
 * series = XLALCreateREAL8TimeSeries(chname, &epoch, 0.0, 0.0, &lalDimensionlessUnit, 0);
 *
 * // update series metadata, particularly the sample rate deltaT
 * XLALFrStreamGetREAL8TimeSeriesMetadata(series, stream);
 *
 * // allocate the required amount of data to the series
 * series = XLALResizeREAL8TimeSeries(series, 0, duration / series->deltaT);
 *
 * // read data from stream one second at a time
 * while (! XLALFrStreamEnd(stream)) {
 *     // get the next second of data
 *     XLALFrStreamGetREAL8TimeSeries(series, stream);
 *     ...
 * }
 * ...
 * @endcode
 * @{
 */

/**
 * @fn int XLALFrStreamGetINT2TimeSeries(INT2TimeSeries * series, LALFrStream * stream)
 * @brief Reads time series data from the current position in a #LALFrStream stream.
 * @details
 * This routine is provided a time series @p series that will be populated with
 * data read from @p stream.  The channel name of the data to be read is given
 * by `stream->name` and the amount of data to be read is given by
 * `stream->data->length`.  The remaining metadata quantities in @p stream
 * are populated: `stream->epoch`, `stream->deltaT`, `stream->sampleUnits`. 
 * If `stream->data` is NULL or `stream->data->length` is 0, only the
 * metadata is read; otherwise `stream->data->data` is populated with the
 * channel data.  The data begins at the current position of the stream,
 * but if a gap is encountered, the stream skips forward until a gap-less
 * stretch of data that is long enough is found.  If this happens, the
 * #LAL_FR_STREAM_GAP flag is set on the LALFrStreamState state of the stream.
 * @param[in, out] series Pointer to the time series containing the name of the channel to be read and an allocated sequence of the length of data to be read.
 * @param stream Pointer to a #LALFrStream stream.
 * @retval 0 Success.
 * @retval <0 Read error or end of stream encountered.
 */

/**
 * @fn int XLALFrStreamGetINT4TimeSeries(INT4TimeSeries * series, LALFrStream * stream)
 * @copydoc XLALFrStreamGetINT2TimeSeries()
 */

/**
 * @fn int XLALFrStreamGetINT8TimeSeries(INT8TimeSeries * series, LALFrStream * stream)
 * @copydoc XLALFrStreamGetINT2TimeSeries()
 */

/**
 * @fn int XLALFrStreamGetUINT2TimeSeries(UINT2TimeSeries * series, LALFrStream * stream)
 * @copydoc XLALFrStreamGetINT2TimeSeries()
 */

/**
 * @fn int XLALFrStreamGetUINT4TimeSeries(UINT4TimeSeries * series, LALFrStream * stream)
 * @copydoc XLALFrStreamGetINT2TimeSeries()
 */

/**
 * @fn int XLALFrStreamGetUINT8TimeSeries(UINT8TimeSeries * series, LALFrStream * stream)
 * @copydoc XLALFrStreamGetINT2TimeSeries()
 */

/**
 * @fn int XLALFrStreamGetREAL4TimeSeries(REAL4TimeSeries * series, LALFrStream * stream)
 * @copydoc XLALFrStreamGetINT2TimeSeries()
 */

/**
 * @fn int XLALFrStreamGetREAL8TimeSeries(REAL8TimeSeries * series, LALFrStream * stream)
 * @copydoc XLALFrStreamGetINT2TimeSeries()
 */

/**
 * @fn int XLALFrStreamGetCOMPLEX8TimeSeries(COMPLEX8TimeSeries * series, LALFrStream * stream)
 * @copydoc XLALFrStreamGetINT2TimeSeries()
 */

/**
 * @fn int XLALFrStreamGetCOMPLEX16TimeSeries(COMPLEX16TimeSeries * series, LALFrStream * stream)
 * @copydoc XLALFrStreamGetINT2TimeSeries()
 */

/**
 * @fn int XLALFrStreamGetINT2TimeSeriesMetadata(INT2TimeSeries * series, LALFrStream * stream)
 * @brief Reads time series metadata from the current position in a #LALFrStream stream.
 * @details
 * This routine is provided a time series @p series that will be populated with
 * metadata read from @p stream.  The metadata is associated with the channel
 * whose name is specified by by `stream->name`.  The remaining metadata
 * quantities in @p stream are populated: `stream->epoch`, `stream->deltaT`,
 * `stream->sampleUnits`.  No data is actually read by this routine.
 * @param[in, out] series Pointer to the time series containing the name of the channel to be read and an allocated sequence of the length of data to be read.
 * @param stream Pointer to a #LALFrStream stream.
 * @retval 0 Success.
 * @retval <0 Failure.
 */

/**
 * @fn int XLALFrStreamGetINT4TimeSeriesMetadata(INT4TimeSeries * series, LALFrStream * stream)
 * @copydoc XLALFrStreamGetINT2TimeSeriesMetadata()
 */

/**
 * @fn int XLALFrStreamGetINT8TimeSeriesMetadata(INT8TimeSeries * series, LALFrStream * stream)
 * @copydoc XLALFrStreamGetINT2TimeSeriesMetadata()
 */

/**
 * @fn int XLALFrStreamGetUINT2TimeSeriesMetadata(UINT2TimeSeries * series, LALFrStream * stream)
 * @copydoc XLALFrStreamGetINT2TimeSeriesMetadata()
 */

/**
 * @fn int XLALFrStreamGetUINT4TimeSeriesMetadata(UINT4TimeSeries * series, LALFrStream * stream)
 * @copydoc XLALFrStreamGetINT2TimeSeriesMetadata()
 */

/**
 * @fn int XLALFrStreamGetUINT8TimeSeriesMetadata(UINT8TimeSeries * series, LALFrStream * stream)
 * @copydoc XLALFrStreamGetINT2TimeSeriesMetadata()
 */

/**
 * @fn int XLALFrStreamGetREAL4TimeSeriesMetadata(REAL4TimeSeries * series, LALFrStream * stream)
 * @copydoc XLALFrStreamGetINT2TimeSeriesMetadata()
 */

/**
 * @fn int XLALFrStreamGetREAL8TimeSeriesMetadata(REAL8TimeSeries * series, LALFrStream * stream)
 * @copydoc XLALFrStreamGetINT2TimeSeriesMetadata()
 */

/**
 * @fn int XLALFrStreamGetCOMPLEX8TimeSeriesMetadata(COMPLEX8TimeSeries * series, LALFrStream * stream)
 * @copydoc XLALFrStreamGetINT2TimeSeriesMetadata()
 */

/**
 * @fn int XLALFrStreamGetCOMPLEX16TimeSeriesMetadata(COMPLEX16TimeSeries * series, LALFrStream * stream)
 * @copydoc XLALFrStreamGetINT2TimeSeriesMetadata()
 */

/** @} */

/**
 * @name Sequential Frequency Series Reading Routines
 *
 * These routines are useful for reading sequential blocks of data from
 * a #LALFrStream stream.  The following is an example of how to read
 * REAL8 data from channel `STRAIN_PSD` in from each frame starting at the
 * beginning of a stream from frame files in the current directory.
 *
 * @code
 * #include <stdio.h>
 * #include <lal/Sequence.h>
 * #include <lal/LALFrStream.h>
 * ...
 * const char *dirname = ".";
 * const char *pattern = "*.gwf";
 * const char *chname = "STRAIN_PSD";
 * LALFrStream *stream;
 * UINT4 length;
 * ...
 * stream = XLALFrStreamOpen(dirname, pattern);
 *
 * // read frequency series from each frame in stream sequentially
 * while (! XLALFrStreamEnd(stream)) {
 *     REAL8FrequencySeries series;
 *     snprintf(series.name, sizeof(series.name), "%s", chname);
 *     XLALFrStreamGetREAL8FrequencySeries(&series, stream);
 *     XLALFrStreamNext(stream);
 *     ...
 *     XLALDestroyREAL8Sequence(&series.data);
 * }
 * ...
 * @endcode
 * @{
 */

/**
 * @fn int XLALFrStreamGetREAL4FrequencySeries(REAL4FrequencySeries * series, LALFrStream * stream);
 * @brief Reads frequency series data from the current frame in a #LALFrStream
 * stream.
 * @details
 * This routine is provided a frequency series @p series that will be populated 
 * with data read from the current frame in @p stream.  The channel name of
 * the data to read is given by `series->name`. The metadata contained in
 * @p series is over-written with the metadata of the named frequency series in
 * @p stream, and `series->data` will be allocated and set to the frequency
 * series data contained in the frame.
 *
 * @attention
 * `series->data` is over-written by this routine; it should be unallocated
 * prior to calling this routine.
 *
 * @param[in,out] series Pointer to the frequency series containing the name
 * of the channel to be read.
 * @param stream Pointer to a #LALFrStream stream.
 *
 * @retval 0 Success.
 * @retval <0 Failure.
 */

/**
 * @fn int XLALFrStreamGetREAL8FrequencySeries(REAL8FrequencySeries * series, LALFrStream * stream);
 * @copydoc XLALFrStreamGetREAL4FrequencySeries()
 */

/**
 * @fn int XLALFrStreamGetCOMPLEX8FrequencySeries(COMPLEX8FrequencySeries * series, LALFrStream * stream);
 * @copydoc XLALFrStreamGetREAL4FrequencySeries()
 */

/**
 * @fn int XLALFrStreamGetCOMPLEX16FrequencySeries(COMPLEX16FrequencySeries * series, LALFrStream * stream);
 * @copydoc XLALFrStreamGetREAL4FrequencySeries()
 */

/** @} */

/** @} */
