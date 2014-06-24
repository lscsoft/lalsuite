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

#include "config.h"
#include <stdlib.h>
#include <string.h>

#include <lal/XLALError.h>
#include <lal/LALFrameU.h>

enum {
    LAL_FRAMEU_FRAME_LIBRARY_UNAVAILABLE,
    LAL_FRAMEU_FRAME_LIBRARY_FRAMEL,
    LAL_FRAMEU_FRAME_LIBRARY_FRAMEC,
};

/* enable FrameL support if available */
#if defined HAVE_FRAMEL_H && defined HAVE_LIBFRAME
#   include "LALFrameUFrameL.h"
#   define CASE_FRAMEL(errval, function, ...) case LAL_FRAMEU_FRAME_LIBRARY_FRAMEL: return function ## _FrameL_ (__VA_ARGS__)
#   ifndef LAL_FRAMEU_FRAME_LIBRARY_DEFAULT
#       define LAL_FRAMEU_FRAME_LIBRARY_DEFAULT LAL_FRAMEU_FRAME_LIBRARY_FRAMEL
#   endif
#else
#   define CASE_FRAMEL(errval, function, ...) case LAL_FRAMEU_FRAME_LIBRARY_FRAMEL: XLAL_ERROR_VAL(errval, XLAL_EERR, "FrameL library unavailable")
#endif

/* enable FrameC support if available */
#if defined HAVE_FRAMECPPC_FRAMEC_H && defined HAVE_LIBFRAMECPPC
#   include "LALFrameUFrameC.h"
#   define CASE_FRAMEC(errval, function, ...) case LAL_FRAMEU_FRAME_LIBRARY_FRAMEC: return function ## _FrameC_ (__VA_ARGS__)
#   ifndef LAL_FRAMEU_FRAME_LIBRARY_DEFAULT
#       define LAL_FRAMEU_FRAME_LIBRARY_DEFAULT LAL_FRAMEU_FRAME_LIBRARY_FRAMEC
#   endif
#else
#   define CASE_FRAMEC(errval, function, ...) case LAL_FRAMEU_FRAME_LIBRARY_FRAMEC: XLAL_ERROR_VAL(errval, XLAL_EERR, "FrameC library unavailable")
#endif

/* fall-back: no frame library available */
/* actually, just kill compilation... */
#ifndef LAL_FRAMEU_FRAME_LIBRARY_DEFAULT
#define LAL_FRAMEU_FRAME_LIBRARY_DEFAULT LAL_FRAMEU_FRAME_LIBRARY_UNAVAILABLE
#error No frame library available
#endif

#define FRAME_LIBRARY_SELECT_VAL(errval, function, ...) \
    do { \
        switch (XLALFrameLibrary()) { \
        CASE_FRAMEL(errval, function, __VA_ARGS__); \
        CASE_FRAMEC(errval, function, __VA_ARGS__); \
        default: \
            XLAL_ERROR_VAL(errval, XLAL_EERR, "No frame library available"); \
        } \
    } while (0)

#define FRAME_LIBRARY_SELECT_VOID(function, ...) FRAME_LIBRARY_SELECT_VAL(/*void*/, function, __VA_ARGS__)
#define FRAME_LIBRARY_SELECT_NULL(function, ...) FRAME_LIBRARY_SELECT_VAL(NULL, function, __VA_ARGS__)
#define FRAME_LIBRARY_SELECT_REAL8(function, ...) FRAME_LIBRARY_SELECT_VAL(XLAL_REAL8_FAIL_NAN, function, __VA_ARGS__)
#define FRAME_LIBRARY_SELECT(function, ...) FRAME_LIBRARY_SELECT_VAL(XLAL_FAILURE, function, __VA_ARGS__)

/* 
 * Routine that returns selected frame library:
 * if LAL_FRAME_LIBRARY is set, use the value from that environment;
 * otherwise use the default value.
 * Note: this is NOT threadsafe, but I doubt the frame libraries are
 * threadsafe in any case....
 */
static int XLALFrameLibrary(void)
{
    static int lalFrameLibrary = -1;
    if (lalFrameLibrary < 0) {
        const char *env = getenv("LAL_FRAME_LIBRARY");
        if (env) {
            if (strcmp(env, "FrameL") == 0) {
#if defined HAVE_FRAMEL_H && defined HAVE_LIBFRAME
                lalFrameLibrary = LAL_FRAMEU_FRAME_LIBRARY_FRAMEL;
#else
                XLAL_ERROR_VAL(LAL_FRAMEU_FRAME_LIBRARY_UNAVAILABLE, XLAL_ESYS,
                    "LAL_FRAME_LIBRARY=%s: FrameL frame library not available", env);
#endif
            } else if (strcmp(env, "FrameC") == 0) {
#if defined HAVE_FRAMECPPC_FRAMEC_H && defined HAVE_LIBFRAMECPPC
                lalFrameLibrary = LAL_FRAMEU_FRAME_LIBRARY_FRAMEC;
#else
                XLAL_ERROR_VAL(LAL_FRAMEU_FRAME_LIBRARY_UNAVAILABLE, XLAL_ESYS,
                    "LAL_FRAME_LIBRARY=%s: FrameC frame library not available", env);
#endif
            } else {
                XLAL_PRINT_WARNING
                    ("LAL_FRAME_LIBRARY=%s: No valid frame library specified [expect: \"FrameL\" or \"FrameC\"]");
                lalFrameLibrary = LAL_FRAMEU_FRAME_LIBRARY_DEFAULT;
            }
        } else {
            lalFrameLibrary = LAL_FRAMEU_FRAME_LIBRARY_DEFAULT;
        }
        switch (lalFrameLibrary) {
        case LAL_FRAMEU_FRAME_LIBRARY_FRAMEL:
            XLAL_PRINT_INFO("Using the FrameL frame library");
            break;
        case LAL_FRAMEU_FRAME_LIBRARY_FRAMEC:
            XLAL_PRINT_INFO("Using the FrameC frame library");
            break;
        default:
            XLAL_ERROR_VAL(LAL_FRAMEU_FRAME_LIBRARY_UNAVAILABLE, XLAL_EERR, "No frame library available");
        }
    }
    return lalFrameLibrary;
}

void XLALFrameUFrFileClose(LALFrameUFrFile * stream)
{
    FRAME_LIBRARY_SELECT_VOID(XLALFrameUFrFileClose, stream);
}

LALFrameUFrFile *XLALFrameUFrFileOpen(const char *filename, const char *mode)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrFileOpen, filename, mode);
}

int XLALFrameUFileCksumValid(LALFrameUFrFile * stream)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFileCksumValid, stream);
}

void XLALFrameUFrTOCFree(LALFrameUFrTOC * toc)
{
    FRAME_LIBRARY_SELECT_VOID(XLALFrameUFrTOCFree, toc);
}

LALFrameUFrTOC *XLALFrameUFrTOCRead(LALFrameUFrFile * stream)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrTOCRead, stream);
}

size_t XLALFrameUFrTOCQueryNFrame(const LALFrameUFrTOC * toc)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrTOCQueryNFrame, toc);
}

double XLALFrameUFrTOCQueryGTimeModf(double *iptr, const LALFrameUFrTOC * toc, size_t pos)
{
    FRAME_LIBRARY_SELECT_REAL8(XLALFrameUFrTOCQueryGTimeModf, iptr, toc, pos);
}

double XLALFrameUFrTOCQueryDt(const LALFrameUFrTOC * toc, size_t pos)
{
    FRAME_LIBRARY_SELECT_REAL8(XLALFrameUFrTOCQueryDt, toc, pos);
}

size_t XLALFrameUFrTOCQueryAdcN(const LALFrameUFrTOC * toc)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrTOCQueryAdcN, toc);
}

const char *XLALFrameUFrTOCQueryAdcName(const LALFrameUFrTOC * toc, size_t adc)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrTOCQueryAdcName, toc, adc);
}

size_t XLALFrameUFrTOCQuerySimN(const LALFrameUFrTOC * toc)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrTOCQuerySimN, toc);
}

const char *XLALFrameUFrTOCQuerySimName(const LALFrameUFrTOC * toc, size_t sim)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrTOCQuerySimName, toc, sim);
}

size_t XLALFrameUFrTOCQueryProcN(const LALFrameUFrTOC * toc)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrTOCQueryProcN, toc);
}

const char *XLALFrameUFrTOCQueryProcName(const LALFrameUFrTOC * toc, size_t proc)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrTOCQueryProcName, toc, proc);
}

size_t XLALFrameUFrTOCQueryDetectorN(const LALFrameUFrTOC * toc)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrTOCQueryDetectorN, toc);
}

const char *XLALFrameUFrTOCQueryDetectorName(const LALFrameUFrTOC * toc, size_t det)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrTOCQueryDetectorName, toc, det);
}

void XLALFrameUFrameHFree(LALFrameUFrameH * frame)
{
    FRAME_LIBRARY_SELECT_VOID(XLALFrameUFrameHFree, frame);
}

LALFrameUFrameH *XLALFrameUFrameHAlloc(const char *name, double start, double dt, int frnum)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrameHAlloc, name, start, dt, frnum);
}

LALFrameUFrameH *XLALFrameUFrameHRead(LALFrameUFrFile * stream, int pos)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrameHRead, stream, pos);
}

int XLALFrameUFrameHWrite(LALFrameUFrFile * stream, LALFrameUFrameH * frame)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrameHWrite, stream, frame);
}

int XLALFrameUFrameHFrChanAdd(LALFrameUFrameH * frame, LALFrameUFrChan * channel)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrameHFrChanAdd, frame, channel);
}

int XLALFrameUFrameHFrDetectorAdd(LALFrameUFrameH * frame, LALFrameUFrDetector * detector)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrameHFrDetectorAdd, frame, detector);
}

int XLALFrameUFrameHFrHistoryAdd(LALFrameUFrameH * frame, LALFrameUFrHistory * history)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrameHFrHistoryAdd, frame, history);
}

const char *XLALFrameUFrameHQueryName(const LALFrameUFrameH * frame)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrameHQueryName, frame);
}

int XLALFrameUFrameHQueryRun(const LALFrameUFrameH * frame)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrameHQueryRun, frame);
}

int XLALFrameUFrameHQueryFrame(const LALFrameUFrameH * frame)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrameHQueryFrame, frame);
}

int XLALFrameUFrameHQueryDataQuality(const LALFrameUFrameH * frame)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrameHQueryDataQuality, frame);
}

double XLALFrameUFrameHQueryGTimeModf(double *iptr, const LALFrameUFrameH * frame)
{
    FRAME_LIBRARY_SELECT_REAL8(XLALFrameUFrameHQueryGTimeModf, iptr, frame);
}

int XLALFrameUFrameHQueryULeapS(const LALFrameUFrameH * frame)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrameHQueryULeapS, frame);
}

double XLALFrameUFrameHQueryDt(const LALFrameUFrameH * frame)
{
    FRAME_LIBRARY_SELECT_REAL8(XLALFrameUFrameHQueryDt, frame);
}

int XLALFrameUFrameHSetRun(LALFrameUFrameH * frame, int run)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrameHSetRun, frame, run);
}

void XLALFrameUFrChanFree(LALFrameUFrChan * channel)
{
    FRAME_LIBRARY_SELECT_VOID(XLALFrameUFrChanFree, channel);
}

LALFrameUFrChan *XLALFrameUFrChanRead(LALFrameUFrFile * stream, const char *name, size_t pos)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrChanRead, stream, name, pos);
}

LALFrameUFrChan *XLALFrameUFrAdcChanAlloc(const char *name, int dtype, size_t ndata)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrAdcChanAlloc, name, dtype, ndata);
}

LALFrameUFrChan *XLALFrameUFrSimChanAlloc(const char *name, int dtype, size_t ndata)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrSimChanAlloc, name, dtype, ndata);
}

LALFrameUFrChan *XLALFrameUFrProcChanAlloc(const char *name, int type, int subtype, int dtype, size_t ndata)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrProcChanAlloc, name, type, subtype, dtype, ndata);
}

const char *XLALFrameUFrChanQueryName(const LALFrameUFrChan * channel)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrChanQueryName, channel);
}

double XLALFrameUFrChanQueryTimeOffset(const LALFrameUFrChan * channel)
{
    FRAME_LIBRARY_SELECT_REAL8(XLALFrameUFrChanQueryTimeOffset, channel);
}

int XLALFrameUFrChanSetSampleRate(LALFrameUFrChan * channel, double sampleRate)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrChanSetSampleRate, channel, sampleRate);
}

int XLALFrameUFrChanSetTimeOffset(LALFrameUFrChan * channel, double timeOffset)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrChanSetTimeOffset, channel, timeOffset);
}

int XLALFrameUFrChanVectorAlloc(LALFrameUFrChan * channel, int dtype, size_t ndata)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrChanVectorAlloc, channel, dtype, ndata);
}

int XLALFrameUFrChanVectorCompress(LALFrameUFrChan * channel, int compressLevel)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrChanVectorCompress, channel, compressLevel);
}

int XLALFrameUFrChanVectorExpand(LALFrameUFrChan * channel)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrChanVectorExpand, channel);
}

const char *XLALFrameUFrChanVectorQueryName(const LALFrameUFrChan * channel)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrChanVectorQueryName, channel);
}

int XLALFrameUFrChanVectorQueryCompress(const LALFrameUFrChan * channel)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrChanVectorQueryCompress, channel);
}

int XLALFrameUFrChanVectorQueryType(const LALFrameUFrChan * channel)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrChanVectorQueryType, channel);
}

void *XLALFrameUFrChanVectorQueryData(const LALFrameUFrChan * channel)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrChanVectorQueryData, channel);
}

size_t XLALFrameUFrChanVectorQueryNBytes(const LALFrameUFrChan * channel)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrChanVectorQueryNBytes, channel);
}

size_t XLALFrameUFrChanVectorQueryNData(const LALFrameUFrChan * channel)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrChanVectorQueryNData, channel);
}

size_t XLALFrameUFrChanVectorQueryNDim(const LALFrameUFrChan * channel)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrChanVectorQueryNDim, channel);
}

size_t XLALFrameUFrChanVectorQueryNx(const LALFrameUFrChan * channel, size_t dim)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrChanVectorQueryNx, channel, dim);
}

double XLALFrameUFrChanVectorQueryDx(const LALFrameUFrChan * channel, size_t dim)
{
    FRAME_LIBRARY_SELECT_REAL8(XLALFrameUFrChanVectorQueryDx, channel, dim);
}

double XLALFrameUFrChanVectorQueryStartX(const LALFrameUFrChan * channel, size_t dim)
{
    FRAME_LIBRARY_SELECT_REAL8(XLALFrameUFrChanVectorQueryStartX, channel, dim);
}

const char *XLALFrameUFrChanVectorQueryUnitX(const LALFrameUFrChan * channel, size_t dim)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrChanVectorQueryUnitX, channel, dim);
}

const char *XLALFrameUFrChanVectorQueryUnitY(const LALFrameUFrChan * channel)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrChanVectorQueryUnitY, channel);
}

int XLALFrameUFrChanVectorSetName(LALFrameUFrChan * channel, const char *name)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrChanVectorSetName, channel, name);
}

int XLALFrameUFrChanVectorSetDx(LALFrameUFrChan * channel, double dx)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrChanVectorSetDx, channel, dx);
}

int XLALFrameUFrChanVectorSetStartX(LALFrameUFrChan * channel, double x0)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrChanVectorSetStartX, channel, x0);
}

int XLALFrameUFrChanVectorSetUnitX(LALFrameUFrChan * channel, const char *unit)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrChanVectorSetUnitX, channel, unit);
}

int XLALFrameUFrChanVectorSetUnitY(LALFrameUFrChan * channel, const char *unit)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrChanVectorSetUnitY, channel, unit);
}

void XLALFrameUFrDetectorFree(LALFrameUFrDetector * detector)
{
    FRAME_LIBRARY_SELECT_VOID(XLALFrameUFrDetectorFree, detector);
}

LALFrameUFrDetector *XLALFrameUFrDetectorRead(LALFrameUFrFile * stream, const char *name)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrDetectorRead, stream, name);
}

LALFrameUFrDetector *XLALFrameUFrDetectorAlloc(const char *name, const char *prefix, double latitude, double longitude,
    double elevation, double azimuthX, double azimuthY, double altitudeX, double altitudeY, double midpointX, double midpointY,
    int localTime)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrDetectorAlloc, name, prefix, latitude, longitude, elevation, azimuthX, azimuthY,
        altitudeX, altitudeY, midpointX, midpointY, localTime);
}

const char *XLALFrameUFrDetectorQueryName(const LALFrameUFrDetector * detector)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrDetectorQueryName, detector);
}

const char *XLALFrameUFrDetectorQueryPrefix(const LALFrameUFrDetector * detector)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrDetectorQueryPrefix, detector);
}

double XLALFrameUFrDetectorQueryLongitude(const LALFrameUFrDetector * detector)
{
    FRAME_LIBRARY_SELECT_REAL8(XLALFrameUFrDetectorQueryLongitude, detector);
}

double XLALFrameUFrDetectorQueryLatitude(const LALFrameUFrDetector * detector)
{
    FRAME_LIBRARY_SELECT_REAL8(XLALFrameUFrDetectorQueryLatitude, detector);
}

double XLALFrameUFrDetectorQueryElevation(const LALFrameUFrDetector * detector)
{
    FRAME_LIBRARY_SELECT_REAL8(XLALFrameUFrDetectorQueryElevation, detector);
}

double XLALFrameUFrDetectorQueryArmXAzimuth(const LALFrameUFrDetector * detector)
{
    FRAME_LIBRARY_SELECT_REAL8(XLALFrameUFrDetectorQueryArmXAzimuth, detector);
}

double XLALFrameUFrDetectorQueryArmYAzimuth(const LALFrameUFrDetector * detector)
{
    FRAME_LIBRARY_SELECT_REAL8(XLALFrameUFrDetectorQueryArmYAzimuth, detector);
}

double XLALFrameUFrDetectorQueryArmXAltitude(const LALFrameUFrDetector * detector)
{
    FRAME_LIBRARY_SELECT_REAL8(XLALFrameUFrDetectorQueryArmXAltitude, detector);
}

double XLALFrameUFrDetectorQueryArmYAltitude(const LALFrameUFrDetector * detector)
{
    FRAME_LIBRARY_SELECT_REAL8(XLALFrameUFrDetectorQueryArmYAltitude, detector);
}

double XLALFrameUFrDetectorQueryArmXMidpoint(const LALFrameUFrDetector * detector)
{
    FRAME_LIBRARY_SELECT_REAL8(XLALFrameUFrDetectorQueryArmXMidpoint, detector);
}

double XLALFrameUFrDetectorQueryArmYMidpoint(const LALFrameUFrDetector * detector)
{
    FRAME_LIBRARY_SELECT_REAL8(XLALFrameUFrDetectorQueryArmYMidpoint, detector);
}

int XLALFrameUFrDetectorQueryLocalTime(const LALFrameUFrDetector * detector)
{
    FRAME_LIBRARY_SELECT(XLALFrameUFrDetectorQueryLocalTime, detector);
}

void XLALFrameUFrHistoryFree(LALFrameUFrHistory * history)
{
    FRAME_LIBRARY_SELECT_VOID(XLALFrameUFrHistoryFree, history);
}

LALFrameUFrHistory *XLALFrameUFrHistoryAlloc(const char *name, double gpssec, const char *comment)
{
    FRAME_LIBRARY_SELECT_NULL(XLALFrameUFrHistoryAlloc, name, gpssec, comment);
}
