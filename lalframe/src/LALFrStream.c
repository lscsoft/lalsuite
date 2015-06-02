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

/* EXPORTED ROUTINES */

int XLALFrStreamClose(LALFrStream * stream)
{
    if (stream) {
        XLALDestroyCache(stream->cache);
        XLALFrStreamFileClose(stream);
        LALFree(stream);
    }
    return 0;
}

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

int XLALFrStreamGetMode(LALFrStream * stream)
{
    return stream->mode;
}

int XLALFrStreamSetMode(LALFrStream * stream, int mode)
{
    stream->mode = mode;
    /* if checksum mode is turned on, do checksum on current file */
    if ((mode & LAL_FR_STREAM_CHECKSUM_MODE) && (stream->file))
        return XLALFrFileCksumValid(stream->file) ? 0 : -1;
    return 0;
}

int XLALFrStreamState(LALFrStream * stream)
{
    return stream->state;
}

int XLALFrStreamEnd(LALFrStream * stream)
{
    return stream->state & LAL_FR_STREAM_END;
}

int XLALFrStreamError(LALFrStream * stream)
{
    return stream->state & LAL_FR_STREAM_ERR;
}

int XLALFrStreamClearErr(LALFrStream * stream)
{
    stream->state = LAL_FR_STREAM_OK;
    return 0;
}

int XLALFrStreamRewind(LALFrStream * stream)
{
    XLALFrStreamFileClose(stream);
    stream->state = LAL_FR_STREAM_OK;
    if (XLALFrStreamFileOpen(stream, 0) < 0)
        XLAL_ERROR(XLAL_EFUNC);
    return 0;
}

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

int XLALFrStreamTell(LIGOTimeGPS * epoch, LALFrStream * stream)
{
    *epoch = stream->epoch;
    return 0;
}

int XLALFrStreamGetpos(LALFrStreamPos * position, LALFrStream * stream)
{
    position->epoch = stream->epoch;
    position->fnum = stream->fnum;
    position->pos = stream->pos;
    return 0;
}

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
