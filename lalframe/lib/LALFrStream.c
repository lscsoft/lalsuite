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
/**
 * @addtogroup LALFrStream_c
 * @brief Provides routines for opening, closing, positioning, and performing
 * other manipulations on a #LALFrStream.
 *
 * @details
 * A frame stream is like a file stream except that it streams along the set
 * of frames in a set of frame files.  These routines are high-level routines
 * that allow you to extract frame data.  Many of these routines have names
 * similar to the standard C file stream manipulation routines and perform
 * similar functions.
 *
 * The routines XLALFrStreamOpen() and XLALFrStreamClose() are used to open
 * and close a frame stream.  The stream is created by XLALFrStreamOpen(),
 * and must be a pointer to NULL before it is opened.  It must have
 * been created prior to calling XLALFrStreamClose(), and after this call,
 * the stream will be a pointer to  NULL.  The routine
 * XLALFrStreamOpen() requires the user to specify the directory name of the
 * frame files and the head names.  If the directory is  NULL, the
 * routine uses the current director (.).  The head names specifies
 * which files are the wanted files in the specified directory.  Wildcards are
 * allowed.  For example, to get LLO frames only, the head names could be set
 * to L-*.gwf.  If the head name is  NULL, the default value
 * *.gwf is used.  The routine XLALFrStreamCacheOpen() is like
 * XLALFrStreamOpen() except that the list of frame files is taken from a
 * frame file cache.  [In fact, XLALFrStreamOpen() simply uses
 * XLALFrCacheGenerate() and XLALFrStreamCacheOpen() to create the
 * stream.]
 *
 * The routine XLALFrStreamSetMode() is used to change the operating mode
 * of a frame stream, which determines how the routines try to accomodate
 * gaps in data and requests for times when there is no data (e.g., before
 * the beginning of the data, after the end of the data, or in some missing
 * data).
 * The default mode, which is given the value
 * #LAL_FR_STREAM_DEFAULT_MODE, prints warnings if a time requested
 * corresponds to a time when there is no data (but then skips to the first
 * avaliable data) and prints an info message when a gap in the data occurs
 * (but then skips beyond the gap).  This default mode is equal to the
 * combination
 * #LAL_FR_STREAM_VERBOSE_MODE | #LAL_FR_STREAM_IGNOREGAP_MODE | #LAL_FR_STREAM_IGNORETIME_MODE
 * where  #LAL_FR_STREAM_VERBOSE_MODE is equal to the combination
 * #LAL_FR_STREAM_TIMEWARN_MODE | #LAL_FR_STREAM_GAPINFO_MODE.  Use
 * #LAL_FR_STREAM_VERBOSE_MODE to print out warnings when requesting times
 * with no data and print out an info message when a gap in the data is
 * encountered.  Unless the mode is supplemented with
 * #LAL_FR_STREAM_IGNOREGAP_MODE, gaps encountered in the data will cause
 * a routine to exit with a non-zero status code; similarly,
 * #LAL_FR_STREAM_IGNORETIME_MODE prevents routines from failing if a time
 * when there is not data is requested.  Set the mode to
 * #LAL_FR_STREAM_SILENT_MODE to suppress the warning and info messages but
 * still cause routines to fail when data is not available.
 * Note: the default value  #LAL_FR_STREAM_DEFAULT_MODE is assumed initially,
 * but this is not necessarily the recommended mode --- it is adopted for
 * compatibility reasons.
 *
 * The routine XLALFrStreamEnd() determines if the end-of-frame-data flag for
 * the data stream has been set.
 *
 * The routine XLALFrStreamNext() advances the frame stream to the
 * beginning of the next frame.
 *
 * The routine XLALFrStreamRewind() rewinds the frame stream to the first
 * frame.
 *
 * The routine XLALFrStreamSeek() sets the frame stream to a specified time,
 * or the earliest time after the specified time if that time is not available
 * (e.g., if it is before the beginning of the frame stream or if it is in a
 * gap in the frame data).  The routine XLALFrStreamTell() returns the
 * current time within the frame stream.
 *
 * The routine XLALFrStreamGetpos() returns a structure containing the
 * current frame stream position.  The frame stream can later be restored to
 * this position using XLALFrStreamSetpos().
 *
 * @{
 */

#include <config.h>
#ifndef HAVE_GETHOSTNAME_PROTOTYPE
int gethostname(char *name, int len);
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/Date.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/LALCache.h>
#include <lal/LALFrameIO.h>
#include <lal/LALFrStream.h>

/* INTERNAL ROUTINES */
/** @cond */

static int XLALFrStreamFileClose(LALFrStream * stream)
{
    XLALFrFileClose(stream->file);
    stream->file = NULL;
    stream->pos = 0;
    return 0;
}

static int XLALFrStreamFileOpen(LALFrStream * stream, UINT4 fnum)
{
    if (!stream->cache || !stream->cache->list)
        XLAL_ERROR(XLAL_EINVAL, "No files in stream file cache");
    if (fnum >= stream->cache->length)
        XLAL_ERROR(XLAL_EINVAL, "File index too large");
    if (stream->file)
        XLALFrStreamFileClose(stream);
    stream->pos = 0;
    stream->fnum = fnum;
    stream->file = XLALFrFileOpenURL(stream->cache->list[fnum].url);
    if (!stream->file) {
        stream->state |= LAL_FR_STREAM_ERR | LAL_FR_STREAM_URL;
        XLAL_ERROR(XLAL_EFUNC);
    }
    if (stream->mode & LAL_FR_STREAM_CHECKSUM_MODE) {
        if (!XLALFrFileCksumValid(stream->file)) {
            stream->state |= LAL_FR_STREAM_ERR;
            XLALFrStreamFileClose(stream);
            XLAL_ERROR(XLAL_EIO, "Invalid checksum in file %s",
                stream->cache->list[fnum].url);
        }
    }
    XLALFrFileQueryGTime(&stream->epoch, stream->file, 0);
    return 0;
}

/** @endcond */

/* EXPORTED ROUTINES */

/**
 * @name Routines to Open, Close and Get/Set Modes of a LALFrStream
 * @{
 */

/**
 * @brief Closes a LALFrStream
 * @details
 * This routine closes all file pointers and deallocates memory associated
 * with a given #LALFrStream.  It performs no action if @p stream is NULL.
 * @param stream Pointer to the #LALFrStream structure to be closed.
 * @retval 0 Success.
 * @retval <0 Failure.
 */
int XLALFrStreamClose(LALFrStream * stream)
{
    if (stream) {
        XLALDestroyCache(stream->cache);
        XLALFrStreamFileClose(stream);
        LALFree(stream);
    }
    return 0;
}

/**
 * @brief Opens a LALFrStream associated with a LALCache
 * @details
 * This routine creates a #LALFrStream that is a stream associated with
 * the frame files contained in a LALCache.
 * @param cache Pointer to a LALCache structure describing the frame files to stream.
 * @returns Pointer to a newly created #LALFrStream structure.
 * @retval NULL Failure.
 */
LALFrStream *XLALFrStreamCacheOpen(LALCache * cache)
{
    LALFrStream *stream;
    size_t i;

    if (!cache)
        XLAL_ERROR_NULL(XLAL_EFAULT);

    stream = LALCalloc(1, sizeof(*stream));
    if (!stream)
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    stream->cache = XLALCacheDuplicate(cache);

    /* check cache entries for t0 and dt; if these are not set then read
     * the framefile to try to get them */
    for (i = 0; i < stream->cache->length; ++i) {
        if (stream->cache->list[i].t0 == 0 || stream->cache->list[i].dt == 0) {
            LIGOTimeGPS end;
            size_t nFrame;
            if (XLALFrStreamFileOpen(stream, i) < 0) {
                XLALFrStreamClose(stream);
                XLAL_ERROR_NULL(XLAL_EIO);
            }
            nFrame = XLALFrFileQueryNFrame(stream->file);
            stream->cache->list[i].t0 = stream->epoch.gpsSeconds;
            XLALFrFileQueryGTime(&end, stream->file, nFrame - 1);
            XLALGPSAdd(&end, XLALFrFileQueryDt(stream->file, nFrame - 1));
            stream->cache->list[i].dt =
                ceil(XLALGPSGetREAL8(&end)) - stream->cache->list[i].t0;
            XLALFrStreamFileClose(stream);
        }
    }

    /* sort and uniqify the cache */
    if (XLALCacheSort(stream->cache) || XLALCacheUniq(stream->cache)) {
        XLALFrStreamClose(stream);
        XLAL_ERROR_NULL(XLAL_EFUNC);
    }

    stream->mode = LAL_FR_STREAM_DEFAULT_MODE;

    /* open up the first file */
    if (XLALFrStreamFileOpen(stream, 0) < 0) {
        XLALFrStreamClose(stream);
        XLAL_ERROR_NULL(XLAL_EFUNC);
    }
    return stream;
}

/**
 * @brief Opens a LALFrStream for specified frame files.
 * @details
 * This routine creates a #LALFrStream that is a stream associated with
 * the frame files in the specified directory matching the specified pattern.
 * The directory containing the frame files is specified by the parameter
 * @p dirname, or the current directory (@p .) if NULL.  Files matching
 * @p pattern will included in the stream.  Wildcards are allowed.  For
 * example, to get LLO frames only, @p pattern could be set to `L-*.gwf`.
 * If @p pattern is NULL, the default value `*.gwf` is used.
 * @param dirname String containing the director name containing the frame files
 * @param pattern Pattern matching the desired frame files.
 * @returns Pointer to a newly created #LALFrStream structure.
 * @retval NULL Failure.
 */
LALFrStream *XLALFrStreamOpen(const char *dirname, const char *pattern)
{
    LALFrStream *stream;
    LALCache *cache;

    cache = XLALCacheGlob(dirname, pattern);
    if (!cache)
        XLAL_ERROR_NULL(XLAL_EFUNC);

    stream = XLALFrStreamCacheOpen(cache);
    if (!stream)
        XLAL_ERROR_NULL(XLAL_EFUNC);

    XLALDestroyCache(cache);
    return stream;
}

/**
 * @brief Returns the current operating mode of a LALFrStream
 * @details
 * The operating mode of a #LALFrStream determines how routines try
 * to accommodate gaps in data and requests for times when there is
 * no data (e.g., before the beginning of the data, after the end of
 * the data, or in some period of missing data).  See XLALFrStreamSetMode()
 * for a description of the #LALFrStreamMode modes.
 * @param stream Pointer to a #LALFrStream structure whose mode will be
 * determined.
 * @returns
 * The current #LALFrStreamMode mode, which is a bit field flag of indicating
 * the current operating modes of the #LALFrStream.
 */
int XLALFrStreamGetMode(LALFrStream * stream)
{
    return stream->mode;
}

/**
 * @brief Change the operating mode of a LALFrStream
 * @details
 * The operating mode of a #LALFrStream determines how routines try
 * to accommodate gaps in data and requests for times when there is
 * no data (e.g., before the beginning of the data, after the end of
 * the data, or in some period of missing data).
 *
 * The default #LALFrStreamMode mode, which is given the value
 * #LAL_FR_STREAM_DEFAULT_MODE, prints warnings if a time requested
 * corresponds to a time when there is no data (but then skips to the first
 * avaliable data) and prints an info message when a gap in the data occurs
 * (but then skips beyond the gap).  This default mode is equal to the
 * combination
 * #LAL_FR_STREAM_VERBOSE_MODE | #LAL_FR_STREAM_IGNOREGAP_MODE | #LAL_FR_STREAM_IGNORETIME_MODE
 * where  #LAL_FR_STREAM_VERBOSE_MODE is equal to the combination
 * #LAL_FR_STREAM_TIMEWARN_MODE | #LAL_FR_STREAM_GAPINFO_MODE.  Use
 * #LAL_FR_STREAM_VERBOSE_MODE to print out warnings when requesting times
 * with no data and print out an info message when a gap in the data is
 * encountered.  Unless the mode is supplemented with
 * #LAL_FR_STREAM_IGNOREGAP_MODE, gaps encountered in the data will cause
 * a routine to exit with a non-zero status code; similarly,
 * #LAL_FR_STREAM_IGNORETIME_MODE prevents routines from failing if a time
 * when there is not data is requested.  Set the mode to
 * #LAL_FR_STREAM_SILENT_MODE to suppress the warning and info messages but
 * still cause routines to fail when data is not available.
 * To enable frame file checksum checking, set the #LAL_FR_STREAM_CHECKSUM_MODE
 * bit.
 *
 * @note The default value  #LAL_FR_STREAM_DEFAULT_MODE is assumed initially,
 * but this is not necessarily the recommended mode --- it is adopted for
 * compatibility reasons.
 *
 * @param stream Pointer to a #LALFrStream structure whose mode will be changed.
 * @param mode Bit lag field specifying the operating modes.
 * @retval 0 Success.
 * @retval <0 Current file does not pass frame file checksum.
 */
int XLALFrStreamSetMode(LALFrStream * stream, int mode)
{
    stream->mode = mode;
    /* if checksum mode is turned on, do checksum on current file */
    if ((mode & LAL_FR_STREAM_CHECKSUM_MODE) && (stream->file))
        return XLALFrFileCksumValid(stream->file) ? 0 : -1;
    return 0;
}

/** @} */

/**
 * @name Routines for Positioning and Manipulating the State of a LALFrStream
 * @{
 */

/**
 * @brief Gets the current state of a LALFrStream
 * @details
 * Gets the #LALFrStreamState state value in a #LALFrStream structure.  The
 * value returned is a bit field specifying the combination of states described
 * by the #LALFrStreamState bits.
 * @param stream Pointer to a #LALFrStream structure whose state will be returned.
 * @returns
 * The #LALFrStreamState state of the #LALFrStream.
 */
int XLALFrStreamState(LALFrStream * stream)
{
    return stream->state;
}

/**
 * @brief Checks to see if a LALFrStream is at the end of the stream
 * @details
 * Determines if the #LAL_FR_STREAM_END bit is set in the #LALFrStreamState
 * state value in a #LALFrStream structure, indicating that the stream is
 * at the end.
 * @param stream Pointer to a #LALFrStream structure.
 * @retval 0 The #LAL_FR_STREAM_END bit is not set in @p stream.
 * @retval 1 The #LAL_FR_STREAM_END bit is set in @p stream.
 */ 
int XLALFrStreamEnd(LALFrStream * stream)
{
    return stream->state & LAL_FR_STREAM_END;
}

/**
 * @brief Checks to see if a LALFrStream has encountered an error
 * @details
 * Determines if the #LAL_FR_STREAM_ERR bit is set in the #LALFrStreamState
 * state value in a #LALFrStream structure, indicating that there has been
 * an error in the stream.
 * @param stream Pointer to a #LALFrStream structure.
 * @retval 0 The #LAL_FR_STREAM_ERR bit is not set in @p stream.
 * @retval 1 The #LAL_FR_STREAM_ERR bit is set in @p stream.
 */ 
int XLALFrStreamError(LALFrStream * stream)
{
    return stream->state & LAL_FR_STREAM_ERR;
}

/**
 * @brief Resets the state of a LALFrStream
 * @details
 * Sets the #LALFrStreamState state of a #LALFrStream structure to the
 * value #LAL_FR_STREAM_OK.
 * @param stream Pointer to a #LALFrStream structure.
 * @retval 0 Success.
 */ 
int XLALFrStreamClearErr(LALFrStream * stream)
{
    stream->state = LAL_FR_STREAM_OK;
    return 0;
}

/**
 * @brief Rewinds a LALFrStream stream
 * @details
 * Rewinds a #LALFrStream stream to the beginning and resets the state
 * to #LAL_FR_STREAM_OK.
 * @param stream Pointer to a #LALFrStream structure.
 * @retval 0 Success.
 * @retval <0 Failure.
 */
int XLALFrStreamRewind(LALFrStream * stream)
{
    XLALFrStreamFileClose(stream);
    stream->state = LAL_FR_STREAM_OK;
    if (XLALFrStreamFileOpen(stream, 0) < 0)
        XLAL_ERROR(XLAL_EFUNC);
    return 0;
}

/**
 * @brief Advance a LALFrStream stream to the beginning of the next frame
 * @details
 * The position of a LALFrStream is advanced so that the next read will
 * be at the next frame.  If the stream is at the end, the #LAL_FR_STREAM_END
 * bit of the LALFrStreamState state is set, and the routine returns the
 * return code 1.  If there is a gap in the data before the next frame,
 * the #LAL_FR_STREAM_GAP bit of the LALFrStreamState state is set, and the
 * routine returns the return code 2.  If, however, the
 * #LAL_FR_STREAM_IGNOREGAP_MODE bit is not set in the LALFrStreamMode mode
 * then the routine produces an error if a gap is encountered.
 * @param stream Pointer to a #LALFrStream structure.
 * @retval 2 Gap in the data is encountered.
 * @retval 1 End of stream encountered.
 * @retval 0 Normal success.
 * @retval <0 Failure.
 */
int XLALFrStreamNext(LALFrStream * stream)
{
    /* timing accuracy: tenth of a sample interval for a 16kHz fast channel */
    const INT8 tacc = (INT8) floor(0.1 * 1e9 / 16384.0);
    const char *url1;
    const char *url2;
    int pos1;
    int pos2;
    INT8 texp = 0;
    INT8 tact;

    if (stream->state & LAL_FR_STREAM_END)
        return 1;       /* end code */

    /* turn off gap bit */
    stream->state &= ~LAL_FR_STREAM_GAP;

    url2 = url1 = stream->cache->list[stream->fnum].url;
    pos2 = pos1 = stream->pos;

    /* open a new file if necessary */
    if (!stream->file) {
        if (stream->fnum >= stream->cache->length) {
            stream->state |= LAL_FR_STREAM_END;
            return 1;
        }
        if (XLALFrStreamFileOpen(stream, stream->fnum) < 0)
            XLAL_ERROR(XLAL_EFUNC);
    }
    if (stream->file) {
        INT4 nFrame = XLALFrFileQueryNFrame(stream->file);
        if (stream->pos < nFrame) {
            LIGOTimeGPS gpstime;
            XLALGPSToINT8NS(XLALFrFileQueryGTime(&gpstime, stream->file,
                  stream->pos));
            texp =
                XLALGPSToINT8NS(XLALGPSAdd(&gpstime,
                    XLALFrFileQueryDt(stream->file, stream->pos)));
            ++stream->pos;
        }
        if (stream->pos >= nFrame) {
            XLALFrStreamFileClose(stream);
            ++stream->fnum;
        }
        pos2 = stream->pos;
    }
    /* open a new file if necessary */
    if (!stream->file) {
        if (stream->fnum >= stream->cache->length) {
            stream->state |= LAL_FR_STREAM_END;
            return 1;
        }
        if (XLALFrStreamFileOpen(stream, stream->fnum) < 0)
            XLAL_ERROR(XLAL_EFUNC);
        url2 = stream->cache->list[stream->fnum].url;
        pos2 = stream->pos;
    }
    /* compute actual start time of this new frame */
    tact =
        XLALGPSToINT8NS(XLALFrFileQueryGTime(&stream->epoch, stream->file,
            stream->pos));

    /* INT8 is platform dependent, cast to long long for llabs() call */
    if (llabs((long long)(texp - tact)) > tacc) { /* there is a gap */
        stream->state |= LAL_FR_STREAM_GAP;
        if (stream->mode & LAL_FR_STREAM_GAPINFO_MODE) {
            XLAL_PRINT_INFO("Gap in frame data between times %.6f and %.6f",
                1e-9 * texp, 1e-9 * tact);
        }
        if (!(stream->mode & LAL_FR_STREAM_IGNOREGAP_MODE)) {
            XLAL_PRINT_ERROR("Gap in frame data");
            XLAL_PRINT_ERROR("Time %.6f is end of frame %d of file %s",
                1e-9 * texp, pos1, url1);
            XLAL_PRINT_ERROR("Time %.6f is start of frame %d of file %s",
                1e-9 * tact, pos2, url2);
            XLAL_ERROR(XLAL_ETIME);
        }
        return 2;       /* gap code */
    }
    return 0;
}

/**
 * @brief Seeks a LALFrStream stream to data at a given time
 * @details
 * The position of a LALFrStream is set so that the next read will
 * be at the specified time.  #LAL_FR_STREAM_END and #LAL_FR_STREAM_GAP
 * bits are turned off in the #LALFrStreamState state.  If the time is before
 * the beginning of the stream, the stream position is set to the beginning of
 * the stream and the routine returns with code 1.  If the time is after the
 * end of the stream, the #LAL_FR_STREAM_END bit is set in the
 * #LALFrStreamState state, and the routine returns with code 2.  If the time
 * is in a gap in the data, the #LAL_FR_STREAM_GAP bit is set in the
 * #LALFrStreamState state, the position is advanced to the next data, and the
 * routine returns with code 3.  If, however, the
 * #LAL_FR_STREAM_IGNORETIME_MODE bit is not set in the LALFrStreamMode mode
 * then these conditions result in an error.
 * @param stream Pointer to a #LALFrStream structure.
 * @param epoch The LIGOTimeGPS time of the next data to read.
 * @retval 3 Time requested is in a gap in the data.
 * @retval 2 Time requested is after the end of the stream.
 * @retval 1 Time requested is before the beginning of the stream.
 * @retval 0 Normal success.
 * @retval <0 Failure.
 */
int XLALFrStreamSeek(LALFrStream * stream, const LIGOTimeGPS * epoch)
{
    double twant = XLALGPSGetREAL8(epoch);
    LALCacheEntry *entry;

    /* close file if one is open */
    XLALFrStreamFileClose(stream);

    /* clear EOF or GAP states; preserve ERR state */
    if (stream->state & LAL_FR_STREAM_ERR)
        stream->state = LAL_FR_STREAM_ERR;
    else
        stream->state = LAL_FR_STREAM_OK;

    /* is epoch before first file? */
    if (epoch->gpsSeconds < stream->cache->list->t0) {
        XLALFrStreamRewind(stream);
        stream->state |= LAL_FR_STREAM_GAP;
        /* is this reported as an error? */
        if (!(stream->mode & LAL_FR_STREAM_IGNORETIME_MODE)) {
            /* FIXME:  if this is an error, should the stream state say so? */
            /* stream->state |= LAL_FR_STREAM_ERR; */
            XLAL_ERROR(XLAL_ETIME);
        }
        if (stream->mode & LAL_FR_STREAM_TIMEWARN_MODE)
            XLAL_PRINT_WARNING("Requested time %d before first frame",
                epoch->gpsSeconds);
        return 1;       /* before first file code */
    }

    /* seek for the time in the cache */
    entry = XLALCacheEntrySeek(stream->cache, twant);
    if (!entry) {       /* seek failed: only happens if time is past end of cache */
        stream->fnum = stream->cache->length;
        stream->epoch = *epoch;
        stream->state |= LAL_FR_STREAM_END;
        /* is this reported as an error? */
        if (!(stream->mode & LAL_FR_STREAM_IGNORETIME_MODE)) {
            /* FIXME:  if this is an error, should the stream state say so? */
            /* stream->state |= LAL_FR_STREAM_ERR; */
            XLAL_ERROR(XLAL_ETIME);
        }
        if (stream->mode & LAL_FR_STREAM_TIMEWARN_MODE)
            XLAL_PRINT_WARNING("Requested time %d after last frame",
                epoch->gpsSeconds);
        return 2;       /* after last file code */
    }

    /* now we must find the position within the frame file */
    for (stream->fnum = entry - stream->cache->list;
        stream->fnum < stream->cache->length; ++stream->fnum) {
        /* check the file contents to determine the position that matches */
        size_t nFrame;
        if (XLALFrStreamFileOpen(stream, stream->fnum) < 0)
            XLAL_ERROR(XLAL_EFUNC);
        if (epoch->gpsSeconds < stream->cache->list[stream->fnum].t0) {
            /* detect a gap between files */
            stream->state |= LAL_FR_STREAM_GAP;
            break;
        }
        nFrame = XLALFrFileQueryNFrame(stream->file);
        for (stream->pos = 0; stream->pos < (int)nFrame; ++stream->pos) {
            LIGOTimeGPS start;
            int cmp;
            XLALFrFileQueryGTime(&start, stream->file, stream->pos);
            cmp = XLALGPSCmp(epoch, &start);
            if (cmp >= 0
                && XLALGPSDiff(epoch,
                    &start) < XLALFrFileQueryDt(stream->file, stream->pos))
                break;  /* this is the frame! */
            if (cmp < 0) {
                /* detect a gap between frames within a file */
                stream->state |= LAL_FR_STREAM_GAP;
                break;
            }
        }
        if (stream->pos < (int)nFrame)  /* we've found the frame */
            break;
        /* oops... not in this frame file, go on to the next one */
        /* probably the frame file was mis-named.... */
        XLALFrStreamFileClose(stream);
    }

    if (stream->fnum >= stream->cache->length) {
        /* we've gone right to the end without finding it! */
        stream->fnum = stream->cache->length;
        stream->epoch = *epoch;
        stream->state |= LAL_FR_STREAM_END;
        /* is this reported as an error? */
        if (!(stream->mode & LAL_FR_STREAM_IGNORETIME_MODE)) {
            /* FIXME:  if this is an error, should the stream state say so? */
            /* stream->state |= LAL_FR_STREAM_ERR; */
            XLAL_ERROR(XLAL_ETIME);
        }
        if (stream->mode & LAL_FR_STREAM_TIMEWARN_MODE)
            XLAL_PRINT_WARNING("Requested time %d after last frame",
                epoch->gpsSeconds);
        return 2;       /* after last file code */
    }

    /* set the time of the stream */
    if (stream->state & LAL_FR_STREAM_GAP) {
        XLALFrFileQueryGTime(&stream->epoch, stream->file, stream->pos);
        if (stream->mode & LAL_FR_STREAM_TIMEWARN_MODE)
            XLAL_PRINT_WARNING("Requested time %.6f in gap in frame data",
                twant);
        if (!(stream->mode & LAL_FR_STREAM_IGNORETIME_MODE))
            XLAL_ERROR(XLAL_ETIME);
        return 3;       /* in a gap code */
    }
    stream->epoch = *epoch;
    return 0;
}

/**
 * @brief Seeks a LALFrStream stream by a time offset
 * @details
 * The position of a LALFrStream is set so that the next read will
 * be at the specified time offset.  The offset @p dt is a number of
 * seconds relative to the @p whence postion, which can be
 * @p SEEK_SET to seek relative to the beginning of the stream,
 * @p SEEK_CUR to seek relative to the current position of the stream,
 * or @p SEEK_END to seek relative to the end of the stream.
 * The return codes and conditions are the same as XLALFrStreamSeek().
 * @param stream Pointer to a #LALFrStream structure.
 * @param dt The offset time in seconds.
 * @param whence The position whence to seek: one of @p SEEK_SET, @p SEEK_CUR,
 * or @p SEEK_END.
 * @retval 3 Time requested is in a gap in the data.
 * @retval 2 Time requested is after the end of the stream.
 * @retval 1 Time requested is before the beginning of the stream.
 * @retval 0 Normal success.
 * @retval <0 Failure.
 */
int XLALFrStreamSeekO(LALFrStream * stream, double dt, int whence)
{
    LIGOTimeGPS epoch;
    switch (whence) {
    case SEEK_SET:
        if (XLALFrStreamRewind(stream) < 0)
            XLAL_ERROR(XLAL_EFUNC);
        /* FALL THROUGH */
    case SEEK_CUR:
        epoch = stream->epoch;
        break;
    case SEEK_END:
        /* go to the last frame */
        XLALFrStreamFileClose(stream);
        if (XLALFrStreamFileOpen(stream, stream->cache->length - 1) < 0)
            XLAL_ERROR(XLAL_EFUNC);
        if ((stream->pos = XLALFrFileQueryNFrame(stream->file) - 1) < 0)
            XLAL_ERROR(XLAL_EFUNC);
        if (XLALFrFileQueryGTime(&epoch, stream->file, stream->pos) == NULL)
            XLAL_ERROR(XLAL_EFUNC);
        /* add duration of last frame to dt */
        dt += XLALFrFileQueryDt(stream->file, stream->pos);
        break;
    default:
        XLAL_ERROR(XLAL_EINVAL,
            "Invalid whence value: use SEEK_SET, SEEK_CUR, or SEEK_END");
    }
    XLALGPSAdd(&epoch, dt);
    if (XLALFrStreamSeek(stream, &epoch) < 0)
        XLAL_ERROR(XLAL_EFUNC);
    return 0;
}

/**
 * @brief Tells the current time of the current position of a LALFrStream
 * stream
 * @param[out] epoch Pointer to a LIGOTimeGPS structure that will contain the
 * time of the current position of the stream.
 * @param[in] stream Pointer to a #LALFrStream structure.
 * @retval 0 Success.
 */
int XLALFrStreamTell(LIGOTimeGPS * epoch, LALFrStream * stream)
{
    *epoch = stream->epoch;
    return 0;
}

/**
 * @brief Gets the current position of a LALFrStream stream
 * @details
 * The XLALFrStreamGetpos() and XLALFrStreamSetpos() provide the ability
 * to save the position of a #LALFrStream stream and to return the stream
 * to that previously saved position.  This can be useful, e.g., when reading
 * several different channels from the data files.
 * @param[out] position Pointer to a #LALFrStreamPos structure that will save
 * the current position.
 * @param[in] stream Pointer to a #LALFrStream structure.
 * time of the current position of the stream.
 * @retval 0 Success.
 */
int XLALFrStreamGetpos(LALFrStreamPos * position, LALFrStream * stream)
{
    position->epoch = stream->epoch;
    position->fnum = stream->fnum;
    position->pos = stream->pos;
    return 0;
}

/**
 * @brief Sets the current position of a LALFrStream stream
 * @details
 * The XLALFrStreamGetpos() and XLALFrStreamSetpos() provide the ability
 * to save the position of a #LALFrStream stream and to return the stream
 * to that previously saved position.  This can be useful, e.g., when reading
 * several different channels from the data files.
 * @param[in, out] stream Pointer to a #LALFrStream structure.
 * @param[in] position Pointer to a #LALFrStreamPos structure that has a
 * previously saved position.
 * @retval 0 Success.
 */
int XLALFrStreamSetpos(LALFrStream * stream, const LALFrStreamPos * position)
{
    /* clear EOF or GAP states; preserve ERR state */
    if (stream->state & LAL_FR_STREAM_ERR)
        stream->state = LAL_FR_STREAM_ERR;
    else
        stream->state = LAL_FR_STREAM_OK;

    if (stream->fnum != position->fnum) {
        XLALFrStreamFileClose(stream);
        if (position->fnum >= stream->fnum) {
            stream->fnum = stream->cache->length;
            stream->state |= LAL_FR_STREAM_END;
            XLAL_ERROR(XLAL_EINVAL);
        }
        if (XLALFrStreamFileOpen(stream, position->fnum) < 0)
            XLAL_ERROR(XLAL_EFUNC);
    }
    stream->epoch = position->epoch;
    stream->pos = position->pos;
    if (stream->pos > (INT4) XLALFrFileQueryNFrame(stream->file)) {
        stream->state |= LAL_FR_STREAM_ERR;
        XLAL_ERROR(XLAL_EINVAL);
    }
    return 0;
}

/** @} */

/** @} */
