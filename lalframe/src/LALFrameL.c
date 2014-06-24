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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/Date.h>
#include <lal/XLALError.h>
#include <lal/LALFrameL.h>

#include <FrIO.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* resolve incomplete types with FrameL structs and custom structs */
#define tagLALFrameUFrameH		FrameH
#define tagLALFrameUFrTOC		FrTOC
#define tagLALFrameUFrHistory		FrHistory
#include <lal/LALFrameU.h>
#include "LALFrameUFrameL.h"

enum { XLAL_FRAMEU_FR_FILE_MODE_R, XLAL_FRAMEU_FR_FILE_MODE_W };
struct tagLALFrameUFrFile {
    int mode;
    FrFile *handle;
};

enum { XLAL_FRAMEU_FR_CHAN_TYPE_ADC, XLAL_FRAMEU_FR_CHAN_TYPE_PROC,
    XLAL_FRAMEU_FR_CHAN_TYPE_SIM
};
struct tagLALFrameUFrChan {
    int type;
    union {
        FrAdcData *adc;
        FrProcData *proc;
        FrSimData *sim;
    } handle;
};

/* for safety, have a string prefix that is nul-terminated as well */
struct tagLALFrameUFrDetector {
    FrDetector *handle;
    char prefix[3];
};

/* local helper functions */

/*
 * Mac OS X 10.9 (Mavericks) displays the following warning regarding
 * the usage of tempnam():
 *
 * warning: 'tempnam' is deprecated: This function is provided for
 * compatibility reasons only. Due to security concerns inherent in the
 * design of tempnam(3), it is highly recommended that you use mkstemp(3)
 * instead. [-Wdeprecated-declarations]
 *
 * wrap mkstemp() as a drop in replacement.
 */
static int mytmpfd(char *tmpfname)
{
    snprintf(tmpfname, L_tmpnam, "%s/tmp.lal.XXXXXX", P_tmpdir);
    return mkstemp(tmpfname);
}

/* TODO: get XLAL error macros going! */

/*
 * FrFile functions
 */

void XLALFrameUFrFileClose_FrameL_(LALFrameUFrFile * stream)
{
    if (stream) {
        if (stream->mode == XLAL_FRAMEU_FR_FILE_MODE_R)
            FrFileIEnd(stream->handle);
        else
            FrFileOEnd(stream->handle);
        LALFree(stream);
    }
}

LALFrameUFrFile *XLALFrameUFrFileOpen_FrameL_(const char *filename, const char *mode)
{
    union {
        const char *cs;
        char *s;
    } fname = {
    filename};
    LALFrameUFrFile *stream;

    if (*mode != 'r' && *mode != 'w')
        XLAL_ERROR_NULL(XLAL_EINVAL);

    stream = LALCalloc(1, sizeof(*stream));
    if (!stream)
        XLAL_ERROR_NULL(XLAL_ENOMEM);

    if (*mode == 'r') {
        if (filename == NULL || strcmp(filename, "-") == 0) {
            /* hack for reading from stdin: dump to a tmpfile */
            char tmpfname[L_tmpnam + 1];
            FILE *tmpfp;
#define BUFSZ 16384
            char buf[BUFSZ];
#undef BUFSZ
            tmpfp = fdopen(mytmpfd(tmpfname), "w");
            while (fwrite(buf, 1, fread(buf, 1, sizeof(buf), stdin),
                    tmpfp) == sizeof(buf)) ;
            fclose(tmpfp);
            stream->handle = FrFileINew(tmpfname);
            unlink(tmpfname);   /* remove tempfile when closed */
        } else
            stream->handle = FrFileINew(fname.s);
        if (!stream->handle) {
            LALFree(stream);
            XLAL_ERROR_NULL(XLAL_EIO, "FrFileINew failed.");
        }
        stream->mode = XLAL_FRAMEU_FR_FILE_MODE_R;
    } else {
        if (filename == NULL || strcmp(filename, "-") == 0) {
            /* output to stdout */
            struct FrIO *frfd;
            frfd = malloc(sizeof(*frfd));
#ifdef FRIOCFILE
            frfd->fp = stdout;
#else
            frfd->fd = 1;
#endif
            stream->handle = FrFileONewFd(frfd, -1);    /* -1: don't alter compression */
        } else
            stream->handle = FrFileONew(fname.s, -1);   /* -1: don't alter compression */
        if (!stream->handle) {
            LALFree(stream);
            XLAL_ERROR_NULL(XLAL_EIO, "FrFileONew failed.");
        }
        stream->handle->noTOCts = FR_NO;
        stream->mode = XLAL_FRAMEU_FR_FILE_MODE_W;
    }

    return stream;
}

/* FIXME: WARNING: this value might need to change in the future */
#define FR_FILE_HEADER_SIZE 40  /* size of frame file header in bytes */

int XLALFrameUFileCksumValid_FrameL_(LALFrameUFrFile * stream)
{
    FrameH *frame = NULL;
    FRBOOL chkSumFiFlag = stream->handle->chkSumFiFlag;
    FRBOOL chkSumFrFlag = stream->handle->chkSumFrFlag;
    stream->handle->chkSumFiFlag = FR_YES;
    stream->handle->chkSumFrFlag = FR_NO;

    /* sequential read */

    /* FIXME: WARNING: HACK: 
     * the following line may need to be updated in future */
    FrIOSet(stream->handle->frfd, FR_FILE_HEADER_SIZE);

    while ((frame = FrameReadRecycle(stream->handle, frame))) ;
    stream->handle->chkSumFiFlag = chkSumFiFlag;
    stream->handle->chkSumFrFlag = chkSumFrFlag;
    if (stream->handle->error != FR_OK)
        XLAL_ERROR_VAL(0, XLAL_EFAILED, "%s", FrErrorGetHistory());
    if (stream->handle->chkTypeFiRead == 0)
        XLAL_ERROR_VAL(0, XLAL_EFAILED, "Missing checksum");
    if (stream->handle->chkSumFiRead != stream->handle->chkSumFi)
        XLAL_ERROR_VAL(0, XLAL_EFAILED, "Bad checksum");
    FrFileIRewind(stream->handle);
    return 1;
}

/*
 * FrTOC functions
 */

void XLALFrameUFrTOCFree_FrameL_(UNUSED LALFrameUFrTOC * toc)
{
    /* Note: toc is owned by frame file */
    return;
}

LALFrameUFrTOC *XLALFrameUFrTOCRead_FrameL_(LALFrameUFrFile * stream)
{
    return FrTOCReadFull(stream->handle);
}

size_t XLALFrameUFrTOCQueryNFrame_FrameL_(const LALFrameUFrTOC * toc)
{
    return toc->nFrame;
}

double XLALFrameUFrTOCQueryGTimeModf_FrameL_(double *iptr, const LALFrameUFrTOC * toc,
    size_t pos)
{
    if (pos < (size_t) toc->nFrame) {
        LIGOTimeGPS epoch;
        XLALGPSSet(&epoch, toc->GTimeS[pos], toc->GTimeN[pos]);
        return XLALGPSModf(iptr, &epoch);
    }
    XLAL_ERROR_REAL8(XLAL_EINVAL);
}

double XLALFrameUFrTOCQueryDt_FrameL_(const LALFrameUFrTOC * toc, size_t pos)
{
    if (pos >= (size_t)toc->nFrame)
        XLAL_ERROR_REAL8(XLAL_EINVAL, "pos = %zu out of range", pos);
    return toc->dt[pos];
}

size_t XLALFrameUFrTOCQueryAdcN_FrameL_(const LALFrameUFrTOC * toc)
{
    FrTOCts *ts;
    size_t nAdc = 0;
    for (ts = toc->adc; ts != NULL; ts = ts->next)
        ++nAdc;
    return nAdc;
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrTOCQueryAdcName_FrameL_(const LALFrameUFrTOC * toc,
    size_t adc)
{
    FrTOCts *ts = toc->adc;
    size_t i = 0;
    while (ts)
        if (adc == i++)
            return ts->name;
        else
            ts = ts->next;
    XLAL_ERROR_NULL(XLAL_EINVAL, "adc = %zu out of range", adc);
}

size_t XLALFrameUFrTOCQuerySimN_FrameL_(const LALFrameUFrTOC * toc)
{
    FrTOCts *ts;
    size_t nSim = 0;
    for (ts = toc->sim; ts != NULL; ts = ts->next)
        ++nSim;
    return nSim;
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrTOCQuerySimName_FrameL_(const LALFrameUFrTOC * toc,
    size_t sim)
{
    FrTOCts *ts = toc->sim;
    size_t i = 0;
    while (ts)
        if (sim == i++)
            return ts->name;
        else
            ts = ts->next;
    XLAL_ERROR_NULL(XLAL_EINVAL, "sim = %zu out of range", sim);
}

size_t XLALFrameUFrTOCQueryProcN_FrameL_(const LALFrameUFrTOC * toc)
{
    FrTOCts *ts;
    int nProc = 0;
    for (ts = toc->proc; ts != NULL; ts = ts->next)
        ++nProc;
    return nProc;
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrTOCQueryProcName_FrameL_(const LALFrameUFrTOC * toc,
    size_t proc)
{
    FrTOCts *ts = toc->proc;
    size_t i = 0;
    while (ts)
        if (proc == i++)
            return ts->name;
        else
            ts = ts->next;
    XLAL_ERROR_NULL(XLAL_EINVAL, "proc = %zu out of range", proc);
}

size_t XLALFrameUFrTOCQueryDetectorN_FrameL_(const LALFrameUFrTOC * toc)
{
    FrTOCdet *d;
    int nDetector = 0;
    for (d = toc->detector; d != NULL; d = d->next)
        ++nDetector;
    return nDetector;
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrTOCQueryDetectorName_FrameL_(const LALFrameUFrTOC * toc,
    size_t det)
{
    FrTOCdet *d = toc->detector;
    size_t i = 0;
    while (d)
        if (det == i++)
            return d->name;
        else
            d = d->next;
    XLAL_ERROR_NULL(XLAL_EINVAL, "det = %zu out of range", det);
}

/*
 * FrameH functions
 */

void XLALFrameUFrameHFree_FrameL_(LALFrameUFrameH * frame)
{
    FrameFree(frame);
    return;
}

LALFrameUFrameH *XLALFrameUFrameHAlloc_FrameL_(const char *name, double start,
    double dt, int frnum)
{
    LALFrameUFrameH *frame;
    frame = calloc(1, sizeof(*frame));
    if (!frame)
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    frame->classe = FrameHDef();
    frame->name = strdup(name);
    frame->GTimeS = floor(start);
    frame->GTimeN = floor(0.5 + 1e9 * (start - frame->GTimeS));
    frame->frame = frnum;
    frame->dt = dt;
    frame->ULeapS = XLALLeapSeconds(frame->GTimeS);
    return frame;
}

LALFrameUFrameH *XLALFrameUFrameHRead_FrameL_(LALFrameUFrFile * stream, int pos)
{
    LALFrameUFrameH *frame;
    LALFrameUFrameH *copy;

    /* make sure the TOC is read */
    if (stream->handle->toc == NULL)
        if (FrTOCReadFull(stream->handle) == NULL)
            XLAL_ERROR_NULL(XLAL_EIO, "FrTOCReadFull failed");

    /* go to the right position */
    if (FrFileIOSet(stream->handle,
            stream->handle->toc->positionH[pos]) == -1)
        XLAL_ERROR_NULL(XLAL_EIO, "FrFileIOSet failed");

    /* get the frame */
    frame = FrameRead(stream->handle);
    if (!frame)
        XLAL_ERROR_NULL(XLAL_EIO, "FrameRead failed");

    copy = FrameHCopy(frame);
    FrSetIni(stream->handle);
    stream->handle->curFrame = NULL;
    FrameFree(frame);

    return copy;
}

int XLALFrameUFrameHWrite_FrameL_(LALFrameUFrFile * stream, LALFrameUFrameH * frame)
{
    if (FrameWrite(frame, stream->handle) != FR_OK)
        XLAL_ERROR(XLAL_EIO, "FrameWrite failed");
    return 0;
}

/* function to add a channel to a frame */

/* need a bunch of helper routines to copy channels */
/* note: these use system memory management routines */

static FrHistory *XLALFrHistoryCopy(const FrHistory * original)
{
    FrHistory *copy;
    copy = calloc(1, sizeof(*copy));
    if (!copy)
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    memcpy(copy, original, sizeof(*copy));
    /* in case of failure, set to zero/null the things we need
     * to allocate separately */
    copy->name = NULL;
    copy->comment = NULL;
    copy->next = NULL;
    if (original->name) {
        copy->name = strdup(original->name);
        if (!copy->name)
            goto failure;
    }
    if (original->comment) {
        copy->comment = strdup(original->comment);
        if (!copy->comment)
            goto failure;
    }
    if (original->next) {
        copy->next = XLALFrHistoryCopy(original->next);
        if (!copy->next)
            goto failure;
    }
    /* successful return */
    return copy;
  failure:     /* unsuccessful return */
    FrHistoryFree(copy);
    XLAL_ERROR_NULL(XLAL_ENOMEM);
}

static FrAdcData *XLALFrAdcDataCopy(const FrAdcData * original)
{
    FrAdcData *copy;
    copy = calloc(1, sizeof(*copy));
    if (!copy)
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    memcpy(copy, original, sizeof(*copy));
    /* in case of failure, set to zero/null the things we need
     * to allocate separately */
    copy->name = NULL;
    copy->comment = NULL;
    copy->units = NULL;
    copy->data = NULL;
    copy->aux = NULL;
    copy->next = NULL;
    if (original->name) {
        copy->name = strdup(original->name);
        if (!copy->name)
            goto failure;
    }
    if (original->comment) {
        copy->comment = strdup(original->comment);
        if (!copy->comment)
            goto failure;
    }
    if (original->units) {
        copy->units = strdup(original->units);
        if (!copy->units)
            goto failure;
    }
    if (original->data) {
        copy->data = FrVectCopy(original->data);
        if (!copy->data)
            goto failure;
    }
    if (original->aux) {
        copy->aux = FrVectCopy(original->aux);
        if (!copy->aux)
            goto failure;
    }
    if (original->next) {
        copy->next = XLALFrAdcDataCopy(original->next);
        if (!copy->next)
            goto failure;
    }
    /* successful return */
    return copy;
  failure:     /* unsuccessful return */
    FrAdcDataFree(copy);
    XLAL_ERROR_NULL(XLAL_ENOMEM);
}

static FrProcData *XLALFrProcDataCopy(const FrProcData * original)
{
    FrProcData *copy;
    copy = calloc(1, sizeof(*copy));
    if (!copy)
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    memcpy(copy, original, sizeof(*copy));
    /* in case of failure, set to zero/null the things we need
     * to allocate separately */
    copy->name = NULL;
    copy->comment = NULL;
    copy->nAuxParam = 0;
    copy->auxParam = NULL;
    copy->auxParamNames = NULL;
    copy->data = NULL;
    copy->aux = NULL;
    copy->table = NULL;
    copy->history = NULL;
    copy->next = NULL;
    if (original->name) {
        copy->name = strdup(original->name);
        if (!copy->name)
            goto failure;
    }
    if (original->comment) {
        copy->comment = strdup(original->comment);
        if (!copy->comment)
            goto failure;
    }
    if (original->nAuxParam > 0) {
        int i;
        copy->auxParam = calloc(original->nAuxParam, sizeof(*copy->auxParam));
        copy->auxParamNames =
            calloc(original->nAuxParam, sizeof(*copy->auxParamNames));
        if (!copy->auxParam || !copy->auxParamNames)
            goto failure;
        copy->nAuxParam = original->nAuxParam;  /* now set value of nAuxParam */
        memcpy(copy->auxParam, original->auxParam,
            copy->nAuxParam * sizeof(*copy->auxParam));
        for (i = 0; i < copy->nAuxParam; ++i)
            copy->auxParamNames[i] = strdup(original->auxParamNames[i]);
    }
    if (original->data) {
        copy->data = FrVectCopy(original->data);
        if (!copy->data)
            goto failure;
    }
    if (original->aux) {
        copy->aux = FrVectCopy(original->aux);
        if (!copy->aux)
            goto failure;
    }
    if (original->table) {
        copy->table = FrTableCopy(original->table);
        if (!copy->table)
            goto failure;
    }
    if (original->history) {
        copy->history = XLALFrHistoryCopy(original->history);
        if (!copy->history)
            goto failure;
    }
    if (original->next) {
        copy->next = XLALFrProcDataCopy(original->next);
        if (!copy->next)
            goto failure;
    }
    /* successful return */
    return copy;
  failure:     /* unsuccessful return */
    FrProcDataFree(copy);
    XLAL_ERROR_NULL(XLAL_ENOMEM);
}

static FrSimData *XLALFrSimDataCopy(const FrSimData * original)
{
    FrSimData *copy;
    copy = calloc(1, sizeof(*copy));
    if (!copy)
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    memcpy(copy, original, sizeof(*copy));
    /* in case of failure, set to zero/null the things we need
     * to allocate separately */
    copy->name = NULL;
    copy->comment = NULL;
    copy->data = NULL;
    copy->input = NULL;
    copy->table = NULL;
    copy->next = NULL;
    if (original->name) {
        copy->name = strdup(original->name);
        if (!copy->name)
            goto failure;
    }
    if (original->comment) {
        copy->comment = strdup(original->comment);
        if (!copy->comment)
            goto failure;
    }
    if (original->data) {
        copy->data = FrVectCopy(original->data);
        if (!copy->data)
            goto failure;
    }
    if (original->input) {
        copy->input = FrVectCopy(original->input);
        if (!copy->input)
            goto failure;
    }
    if (original->table) {
        copy->table = FrTableCopy(original->table);
        if (!copy->table)
            goto failure;
    }
    if (original->next) {
        copy->next = XLALFrSimDataCopy(original->next);
        if (!copy->next)
            goto failure;
    }
    /* successful return */
    return copy;
  failure:     /* unsuccessful return */
    FrSimDataFree(copy);
    XLAL_ERROR_NULL(XLAL_ENOMEM);
}

int XLALFrameUFrameHFrChanAdd_FrameL_(LALFrameUFrameH * frame,
    LALFrameUFrChan * channel)
{
    LALFrameUFrChan copy;
    switch (channel->type) {
    case XLAL_FRAMEU_FR_CHAN_TYPE_ADC:
        copy.handle.adc = XLALFrAdcDataCopy(channel->handle.adc);
        if (!copy.handle.adc)
            XLAL_ERROR(XLAL_EFUNC);
        if (frame->rawData == NULL)
            FrRawDataNew(frame);
        copy.handle.adc->next = frame->rawData->firstAdc;
        frame->rawData->firstAdc = copy.handle.adc;
        break;
    case XLAL_FRAMEU_FR_CHAN_TYPE_PROC:
        copy.handle.proc = XLALFrProcDataCopy(channel->handle.proc);
        if (!copy.handle.proc)
            XLAL_ERROR(XLAL_EFUNC);
        copy.handle.proc->next = frame->procData;
        frame->procData = copy.handle.proc;
        break;
    case XLAL_FRAMEU_FR_CHAN_TYPE_SIM:
        copy.handle.sim = XLALFrSimDataCopy(channel->handle.sim);
        if (!copy.handle.sim)
            XLAL_ERROR(XLAL_EFUNC);
        copy.handle.sim->next = frame->simData;
        frame->simData = copy.handle.sim;
        break;
    default:
        XLAL_ERROR(XLAL_EINVAL, "Invalid channel type");
    }
    return 0;
}

int XLALFrameUFrameHFrDetectorAdd_FrameL_(LALFrameUFrameH * frame,
    LALFrameUFrDetector * detector)
{
    FrDetector *copy;
    copy = calloc(1, sizeof(*copy));
    if (!copy)
        XLAL_ERROR(XLAL_ENOMEM);

    copy->classe = FrDetectorDef();
    if (detector->handle->name) {
        copy->name = strdup(detector->handle->name);
        if (!copy->name) {
            free(copy);
            XLAL_ERROR(XLAL_ENOMEM);
        }
    }
    if (detector->handle->prefix)
        memcpy(copy->prefix, detector->handle->prefix, 2);

    copy->longitude = detector->handle->longitude;
    copy->latitude = detector->handle->latitude;
    copy->elevation = detector->handle->elevation;
    copy->armXazimuth = detector->handle->armXazimuth;
    copy->armYazimuth = detector->handle->armYazimuth;
    copy->armXaltitude = detector->handle->armXaltitude;
    copy->armYaltitude = detector->handle->armYaltitude;
    copy->armXmidpoint = detector->handle->armXmidpoint;
    copy->armYmidpoint = detector->handle->armYmidpoint;
    copy->localTime = detector->handle->localTime;

    /* now affix it to frame */
    copy->next = frame->detectProc;
    frame->detectProc = copy;
    return 0;
}

int XLALFrameUFrameHFrHistoryAdd_FrameL_(LALFrameUFrameH * frame,
    LALFrameUFrHistory * history)
{
    FrHistory *copy;
    copy =
        FrHistoryNew(history->name ? history->name : frame->name,
        history->time, history->comment);
    if (!copy)
        XLAL_ERROR(XLAL_EIO, "FrHistoryNew failed");
    copy->next = frame->history;
    frame->history = copy;
    return 0;
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrameHQueryName_FrameL_(const LALFrameUFrameH * frame)
{
    return frame->name;
}

int XLALFrameUFrameHQueryRun_FrameL_(const LALFrameUFrameH * frame)
{
    return frame->run;
}

int XLALFrameUFrameHQueryFrame_FrameL_(const LALFrameUFrameH * frame)
{
    return frame->frame;
}

int XLALFrameUFrameHQueryDataQuality_FrameL_(const LALFrameUFrameH * frame)
{
    return frame->dataQuality;
}

double XLALFrameUFrameHQueryGTimeModf_FrameL_(double *iptr,
    const LALFrameUFrameH * frame)
{
    LIGOTimeGPS epoch;
    XLALGPSSet(&epoch, frame->GTimeS, frame->GTimeN);
    return XLALGPSModf(iptr, &epoch);
}

int XLALFrameUFrameHQueryULeapS_FrameL_(const LALFrameUFrameH * frame)
{
    return frame->ULeapS;
}

double XLALFrameUFrameHQueryDt_FrameL_(const LALFrameUFrameH * frame)
{
    return frame->dt;
}

/* functions to set frame metadata */

int XLALFrameUFrameHSetRun_FrameL_(LALFrameUFrameH * frame, int run)
{
    frame->run = run;
    return 0;
}

/* 
 * FrChan functions
 */

void XLALFrameUFrChanFree_FrameL_(LALFrameUFrChan * channel)
{
    if (channel) {
        switch (channel->type) {
        case XLAL_FRAMEU_FR_CHAN_TYPE_ADC:
            FrAdcDataFree(channel->handle.adc);
            break;
        case XLAL_FRAMEU_FR_CHAN_TYPE_PROC:
            FrProcDataFree(channel->handle.proc);
            break;
        case XLAL_FRAMEU_FR_CHAN_TYPE_SIM:
            FrSimDataFree(channel->handle.sim);
            break;
        default:
            break;
        }
        LALFree(channel);
    }
    return;
}

LALFrameUFrChan *XLALFrameUFrChanRead_FrameL_(LALFrameUFrFile * stream,
    const char *name, size_t pos)
{
    LALFrameUFrChan *channel;
    FrTOCts *ts;
    FrTOC *toc;
    double gtime;

    toc = stream->handle->toc;
    if (!toc)
        toc = FrTOCReadFull(stream->handle);
    if (!toc)
        XLAL_ERROR_NULL(XLAL_EIO, "FrTOCReadFull failed");
    if (pos >= (size_t) toc->nFrame)
        XLAL_ERROR_NULL(XLAL_EINVAL, "pos = %zu out of range");

    /* the gps time of the frame at position pos */
    gtime = toc->GTimeS[pos] + 1e-9 * toc->GTimeN[pos];

    /* allocate memory for channel */
    channel = LALCalloc(1, sizeof(*channel));
    if (!channel)
        XLAL_ERROR_NULL(XLAL_ENOMEM);

    /* scan adc channels */
    stream->handle->relocation = FR_NO;
    for (ts = toc->adc; ts != NULL; ts = ts->next)
        if (strcmp(name, ts->name) == 0) {
            FrAdcData *adc;
            if (FrTOCSetPos(stream->handle, ts->position[pos]))
                goto failure;
            adc = FrAdcDataRead(stream->handle);
            if (!adc)
                goto failure;
            gtime += adc->timeOffset;
            //adc->next = NULL;
            adc->data = FrVectReadNext(stream->handle, gtime, adc->name);
            channel->handle.adc = adc;
            channel->type = XLAL_FRAMEU_FR_CHAN_TYPE_ADC;
            return channel;
        }

    /* scan proc channels */
    for (ts = toc->proc; ts != NULL; ts = ts->next)
        if (strcmp(name, ts->name) == 0) {
            FrProcData *proc;
            if (FrTOCSetPos(stream->handle, ts->position[pos]))
                goto failure;
            proc = FrProcDataRead(stream->handle);
            if (!proc)
                goto failure;
            gtime += proc->timeOffset;
            proc->next = NULL;
            proc->data = FrVectReadNext(stream->handle, gtime, proc->name);
            channel->handle.proc = proc;
            channel->type = XLAL_FRAMEU_FR_CHAN_TYPE_PROC;
            return channel;
        }

    /* scan sim channels */
    for (ts = toc->sim; ts != NULL; ts = ts->next)
        if (strcmp(name, ts->name) == 0) {
            FrSimData *sim;
            if (FrTOCSetPos(stream->handle, ts->position[pos]))
                goto failure;
            sim = FrSimDataRead(stream->handle);
            if (!sim)
                goto failure;
            gtime += sim->timeOffset;
            sim->next = NULL;
            sim->data = FrVectReadNext(stream->handle, gtime, sim->name);
            channel->handle.sim = sim;
            channel->type = XLAL_FRAMEU_FR_CHAN_TYPE_SIM;
            return channel;
        }

    /* couldn't find channel */
  failure:
    LALFree(channel);
    XLAL_ERROR_NULL(XLAL_ENAME, "Channel %s not found", name);
}

/* stripped-down copy of FrAdcDataNewF */
static FrAdcData *XLALFrameUFrAdcDataNew(const char *name, int type)
{
    FrAdcData *adcData;
    adcData = calloc(1, sizeof(*adcData));
    if (!adcData)
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    adcData->classe = FrAdcDataDef();
    adcData->name = strdup(name);
    if (!adcData->name) {
        FrAdcDataFree(adcData);
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    }
    switch (type) {
    case FR_VECT_4S:
        adcData->nBits = 32;
        break;
    case FR_VECT_2S:
        adcData->nBits = 16;
        break;
    case FR_VECT_C:
        adcData->nBits = 8;
        break;
    case FR_VECT_4R:
        adcData->nBits = 32;
        break;
    case FR_VECT_8R:
        adcData->nBits = 64;
        break;
    default:   /* invalid type */
        FrAdcDataFree(adcData);
        XLAL_ERROR_NULL(XLAL_ETYPE);
    }
    return adcData;
}

/* stripped-down copy of FrSimDataNew */
static FrSimData *XLALFrameUFrSimDataNew(const char *name)
{
    FrSimData *simData;
    simData = calloc(1, sizeof(*simData));
    if (!simData)
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    simData->classe = FrSimDataDef();
    simData->name = strdup(name);
    if (!simData->name) {
        FrSimDataFree(simData);
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    }
    return simData;
}

/* stripped-down copy of FrProcDataNew */
static FrProcData *XLALFrameUFrProcDataNew(const char *name)
{
    FrProcData *procData;
    procData = calloc(1, sizeof(*procData));
    if (!procData)
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    procData->classe = FrProcDataDef();
    procData->name = strdup(name);
    if (!procData->name) {
        FrProcDataFree(procData);
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    }
    return procData;
}

/* helper function to create new channel of appropriate type */
static LALFrameUFrChan *XLALFrameUFrChanAlloc(const char *name, int chanType,
    int dataType, size_t ndata)
{
    LALFrameUFrChan *channel;
    channel = LALCalloc(1, sizeof(*channel));
    if (!channel)
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    channel->type = chanType;
    switch (chanType) {
    case XLAL_FRAMEU_FR_CHAN_TYPE_ADC:
        channel->handle.adc = XLALFrameUFrAdcDataNew(name, dataType);
        if (!channel->handle.adc) {
            LALFree(channel);
            XLAL_ERROR_NULL(XLAL_EFUNC);
        }
        break;
    case XLAL_FRAMEU_FR_CHAN_TYPE_SIM:
        channel->handle.sim = XLALFrameUFrSimDataNew(name);
        if (!channel->handle.sim) {
            LALFree(channel);
            XLAL_ERROR_NULL(XLAL_EFUNC);
        }
        break;
    case XLAL_FRAMEU_FR_CHAN_TYPE_PROC:
        channel->handle.proc = XLALFrameUFrProcDataNew(name);
        if (!channel->handle.proc) {
            LALFree(channel);
            XLAL_ERROR_NULL(XLAL_EFUNC);
        }
        break;
    default:   /* unrecognized channel type */
        LALFree(channel);
        XLAL_ERROR_NULL(XLAL_ETYPE);
    }
    if (XLALFrameUFrChanVectorAlloc(channel, dataType, ndata) < 0) {
        XLALFrameUFrChanFree(channel);
        XLAL_ERROR_NULL(XLAL_EFUNC);
    }
    return channel;
}

LALFrameUFrChan *XLALFrameUFrAdcChanAlloc_FrameL_(const char *name, int dtype,
    size_t ndata)
{
    return XLALFrameUFrChanAlloc(name, XLAL_FRAMEU_FR_CHAN_TYPE_ADC, dtype,
        ndata);
}

LALFrameUFrChan *XLALFrameUFrSimChanAlloc_FrameL_(const char *name, int dtype,
    size_t ndata)
{
    return XLALFrameUFrChanAlloc(name, XLAL_FRAMEU_FR_CHAN_TYPE_SIM, dtype,
        ndata);
}

LALFrameUFrChan *XLALFrameUFrProcChanAlloc_FrameL_(const char *name, int type,
    int subtype, int dtype, size_t ndata)
{
    LALFrameUFrChan *channel;
    channel =
        XLALFrameUFrChanAlloc(name, XLAL_FRAMEU_FR_CHAN_TYPE_PROC, dtype,
        ndata);
    if (!channel)
        XLAL_ERROR_NULL(XLAL_EFUNC);
    /* set type and subtype metadata */
    channel->handle.proc->type = type;
    channel->handle.proc->subType = subtype;
    return channel;
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrChanQueryName_FrameL_(const LALFrameUFrChan * channel)
{
    switch (channel->type) {
    case XLAL_FRAMEU_FR_CHAN_TYPE_ADC:
        return channel->handle.adc->name;
        break;
    case XLAL_FRAMEU_FR_CHAN_TYPE_SIM:
        return channel->handle.sim->name;
        break;
    case XLAL_FRAMEU_FR_CHAN_TYPE_PROC:
        return channel->handle.proc->name;
        break;
    default:   /* unrecognized channel type */
        XLAL_ERROR_NULL(XLAL_ETYPE);
    }
}

double XLALFrameUFrChanQueryTimeOffset_FrameL_(const LALFrameUFrChan * channel)
{
    switch (channel->type) {
    case XLAL_FRAMEU_FR_CHAN_TYPE_ADC:
        return channel->handle.adc->timeOffset;
        break;
    case XLAL_FRAMEU_FR_CHAN_TYPE_SIM:
        return channel->handle.sim->timeOffset;
        break;
    case XLAL_FRAMEU_FR_CHAN_TYPE_PROC:
        return channel->handle.proc->timeOffset;
        break;
    default:   /* unrecognized channel type */
        XLAL_ERROR_REAL8(XLAL_ETYPE);
    }
}

int XLALFrameUFrChanSetSampleRate_FrameL_(LALFrameUFrChan * channel,
    double sampleRate)
{
    switch (channel->type) {
    case XLAL_FRAMEU_FR_CHAN_TYPE_ADC:
        channel->handle.adc->sampleRate = sampleRate;
        return 0;
    case XLAL_FRAMEU_FR_CHAN_TYPE_SIM:
        channel->handle.sim->sampleRate = sampleRate;
        return 0;
    case XLAL_FRAMEU_FR_CHAN_TYPE_PROC:        /* does not support setting sample rate */
    default:   /* unrecognized channel type */
        XLAL_ERROR(XLAL_ETYPE);
    }
}

int XLALFrameUFrChanSetTimeOffset_FrameL_(LALFrameUFrChan * channel,
    double timeOffset)
{
    if (timeOffset < 0)         /* timeOffset must be non-negative */
        XLAL_ERROR(XLAL_EINVAL, "Time offset must be non-negative");
    switch (channel->type) {
    case XLAL_FRAMEU_FR_CHAN_TYPE_ADC:
        channel->handle.adc->timeOffset = timeOffset;
        return 0;
    case XLAL_FRAMEU_FR_CHAN_TYPE_SIM:
        channel->handle.sim->timeOffset = timeOffset;
        return 0;
    case XLAL_FRAMEU_FR_CHAN_TYPE_PROC:
        channel->handle.sim->timeOffset = timeOffset;
	return 0;
    default:   /* unrecognized channel type */
        XLAL_ERROR(XLAL_ETYPE);
    }
}

/*
 * FrVect functions
 */

/* helper function to allocate vector in channel */
int XLALFrameUFrChanVectorAlloc_FrameL_(LALFrameUFrChan * channel, int dtype,
    size_t ndata)
{
    switch (channel->type) {
    case XLAL_FRAMEU_FR_CHAN_TYPE_ADC:
        /* determine bits */
        switch (dtype) {
        case FR_VECT_4S:
            channel->handle.adc->nBits = 32;
            break;
        case FR_VECT_2S:
            channel->handle.adc->nBits = 16;
            break;
        case FR_VECT_C:
            channel->handle.adc->nBits = 8;
            break;
        case FR_VECT_4R:
            channel->handle.adc->nBits = 32;
            break;
        case FR_VECT_8R:
            channel->handle.adc->nBits = 64;
            break;
        default:       /* invalid vector type for adc data */
            XLAL_ERROR(XLAL_ETYPE);
        }
        /* discard old vector (if present) */
        FrVectFree(channel->handle.adc->data);
        /* create new vector; note negative type means use malloc, not calloc */
        channel->handle.adc->data =
            FrVectNew1D(channel->handle.adc->name, -dtype, ndata, 0.0, NULL,
            NULL);
        if (!channel->handle.adc->data)
            XLAL_ERROR(XLAL_EIO, "FrVectNew1D failed");
        break;
    case XLAL_FRAMEU_FR_CHAN_TYPE_SIM:
        /* make sure type is allowed */
        switch (dtype) {
        case FR_VECT_4S:
        case FR_VECT_2S:
        case FR_VECT_C:
        case FR_VECT_4R:
        case FR_VECT_8R:
            break;
        default:       /* invalid vector type for sim data */
            XLAL_ERROR(XLAL_ETYPE);
        }
        /* discard old vector (if present) */
        FrVectFree(channel->handle.sim->data);
        /* create new vector; note negative type means use malloc, not calloc */
        channel->handle.sim->data =
            FrVectNew1D(channel->handle.sim->name, -dtype, ndata, 0.0, NULL,
            NULL);
        if (!channel->handle.sim->data)
            XLAL_ERROR(XLAL_EIO, "FrVectNew1D failed");
        break;
    case XLAL_FRAMEU_FR_CHAN_TYPE_PROC:
        /* discard old vector (if present) */
        FrVectFree(channel->handle.proc->data);
        /* create new vector; note negative type means use malloc, not calloc */
        channel->handle.proc->data =
            FrVectNew1D(channel->handle.proc->name, -dtype, ndata, 0.0, NULL,
            NULL);
        if (!channel->handle.proc->data)
            XLAL_ERROR(XLAL_EIO, "FrVectNew1D failed");
        break;
    default:   /* invalid channel type */
        XLAL_ERROR(XLAL_ETYPE);
    }
    return 0;
}

/* helper function that gets the FrVect structure from a channel */
static FrVect *XLALFrameUFrChanVectorPtr(const LALFrameUFrChan * channel)
{
    switch (channel->type) {
    case XLAL_FRAMEU_FR_CHAN_TYPE_ADC:
        return channel->handle.adc->data;
    case XLAL_FRAMEU_FR_CHAN_TYPE_SIM:
        return channel->handle.sim->data;
    case XLAL_FRAMEU_FR_CHAN_TYPE_PROC:
        return channel->handle.proc->data;
    default:   /* unrecognized type */
        XLAL_ERROR_NULL(XLAL_ETYPE);
    }
}

int XLALFrameUFrChanVectorCompress_FrameL_(LALFrameUFrChan * channel,
    int compressLevel)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR(XLAL_EFUNC);
    FrVectCompress(vect, compressLevel, 0);
    return 0;
}

int XLALFrameUFrChanVectorExpand_FrameL_(LALFrameUFrChan * channel)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR(XLAL_EFUNC);
    FrVectExpand(vect);
    return 0;
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrChanVectorQueryName_FrameL_(const LALFrameUFrChan * channel)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR_NULL(XLAL_EFUNC);
    return vect->name;
}

int XLALFrameUFrChanVectorQueryCompress_FrameL_(const LALFrameUFrChan * channel)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR(XLAL_EFUNC);
    return vect->compress;
}

int XLALFrameUFrChanVectorQueryType_FrameL_(const LALFrameUFrChan * channel)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR(XLAL_EFUNC);
    return vect->type;
}

/* retrieves a handle to the data vector in the FrVect structure */
void *XLALFrameUFrChanVectorQueryData_FrameL_(const LALFrameUFrChan * channel)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR_NULL(XLAL_EFUNC);
    return vect->data;
}

size_t XLALFrameUFrChanVectorQueryNBytes_FrameL_(const LALFrameUFrChan * channel)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR(XLAL_EFUNC);
    return vect->nBytes;
}

size_t XLALFrameUFrChanVectorQueryNData_FrameL_(const LALFrameUFrChan * channel)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR(XLAL_EFUNC);
    return vect->nData;
}

size_t XLALFrameUFrChanVectorQueryNDim_FrameL_(const LALFrameUFrChan * channel)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR(XLAL_EFUNC);
    return vect->nDim;
}

size_t XLALFrameUFrChanVectorQueryNx_FrameL_(const LALFrameUFrChan * channel,
    size_t dim)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR(XLAL_EFUNC);
    if (dim >= vect->nDim)
        XLAL_ERROR(XLAL_EINVAL, "dim = %zu out of range", dim);
    return vect->nx[dim];
}

double XLALFrameUFrChanVectorQueryDx_FrameL_(const LALFrameUFrChan * channel,
    size_t dim)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR_REAL8(XLAL_EFUNC);
    if (dim >= vect->nDim)
        XLAL_ERROR_REAL8(XLAL_EINVAL, "dim = %zu out of range", dim);
    return vect->dx[dim];
}

double XLALFrameUFrChanVectorQueryStartX_FrameL_(const LALFrameUFrChan * channel,
    size_t dim)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR_REAL8(XLAL_EFUNC);
    if (dim >= vect->nDim)
        XLAL_ERROR_REAL8(XLAL_EINVAL, "dim = %zu out of range", dim);
    return vect->startX[dim];
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrChanVectorQueryUnitX_FrameL_(const LALFrameUFrChan * channel,
    size_t dim)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR_NULL(XLAL_EFUNC);
    if (dim >= vect->nDim)
        XLAL_ERROR_NULL(XLAL_EINVAL, "dim = %zu out of range", dim);
    return vect->unitX[dim];
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrChanVectorQueryUnitY_FrameL_(const LALFrameUFrChan * channel)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR_NULL(XLAL_EFUNC);
    return vect->unitY;
}

int XLALFrameUFrChanVectorSetName_FrameL_(LALFrameUFrChan * channel, const char *name)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR(XLAL_EFUNC);
    vect->name = strdup(name);
    return 0;
}

int XLALFrameUFrChanVectorSetUnitY_FrameL_(LALFrameUFrChan * channel,
    const char *unit)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR(XLAL_EFUNC);
    vect->unitY = strdup(unit);
    return 0;
}

/* NOTE: only support 1-dimensional vectors */

int XLALFrameUFrChanVectorSetDx_FrameL_(LALFrameUFrChan * channel, double dx)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR(XLAL_EFUNC);
    vect->dx[0] = dx;
    return 0;
}

int XLALFrameUFrChanVectorSetStartX_FrameL_(LALFrameUFrChan * channel, double x0)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR(XLAL_EFUNC);
    vect->startX[0] = x0;
    return 0;
}

int XLALFrameUFrChanVectorSetUnitX_FrameL_(LALFrameUFrChan * channel,
    const char *unit)
{
    FrVect *vect;
    vect = XLALFrameUFrChanVectorPtr(channel);
    if (!vect)
        XLAL_ERROR(XLAL_EFUNC);
    vect->unitX[0] = strdup(unit);
    return 0;
}

/*
 * FrDetector functions
 */

void XLALFrameUFrDetectorFree_FrameL_(LALFrameUFrDetector * detector)
{
    if (detector) {
        FrDetectorFree(detector->handle);
        LALFree(detector);
    }
    return;
}

LALFrameUFrDetector *XLALFrameUFrDetectorRead_FrameL_(LALFrameUFrFile * stream,
    const char *name)
{
    LALFrameUFrDetector *detector;
    FrTOCdet *d;
    for (d = stream->handle->toc->detector; d != NULL; d = d->next)
        if (strcmp(d->name, name) == 0) {
            char prefix[3];
            FrDetector *det;
            if (FrTOCSetPos(stream->handle, d->position) != 0)
                XLAL_ERROR_NULL(XLAL_EIO, "FrTOCSetPos failed");
            det = FrDetectorRead(stream->handle);
            prefix[0] = det->prefix[0];
            prefix[1] = det->prefix[1];
            prefix[2] = 0;
            detector =
                XLALFrameUFrDetectorAlloc(det->name, prefix, det->latitude,
                det->longitude, det->elevation, det->armXazimuth,
                det->armYazimuth, det->armXaltitude, det->armYaltitude,
                det->armXmidpoint, det->armYmidpoint, det->localTime);
            return detector;
        }
    /* didn't find a detector of that name */
    XLAL_ERROR_NULL(XLAL_ENAME, "Detector %s not found", name);
}

LALFrameUFrDetector *XLALFrameUFrDetectorAlloc_FrameL_(const char *name,
    const char *prefix, double latitude, double longitude, double elevation,
    double azimuthX, double azimuthY, double altitudeX, double altitudeY,
    double midpointX, double midpointY, int localTime)
{
    LALFrameUFrDetector *detector;
    detector = LALCalloc(1, sizeof(*detector));
    if (!detector)
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    detector->handle = calloc(1, sizeof(*detector->handle));
    if (!detector->handle)
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    detector->handle->classe = FrDetectorDef();
    detector->handle->name = strdup(name);
    if (!detector->handle->name) {
        XLALFrameUFrDetectorFree(detector);
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    }
    if (prefix) {
        memcpy(detector->prefix, prefix, 2);
        memcpy(detector->handle->prefix, prefix, 2);
    }

    detector->handle->longitude = longitude;    /* longitude (east of greenwich) in radians      */
    detector->handle->latitude = latitude;      /* latitude (north of equator) in radians        */
    detector->handle->elevation = elevation;    /* detector altitude (meter)                     */
    detector->handle->armXazimuth = azimuthX;   /* orientation of X arm in radians CW from North */
    detector->handle->armYazimuth = azimuthY;   /* orientation of Y arm in radians CW from North */
    /* Azimuth values should be in the range 0 to 2pi */
    detector->handle->armXaltitude = altitudeX; /* altitude (pitch) of the X arm                 */
    detector->handle->armYaltitude = altitudeY; /* altitude (pitch) of the Y arm                 */
    detector->handle->armXmidpoint = midpointX; /* vertex to  middle of the X arm distance       */
    detector->handle->armYmidpoint = midpointY; /* vertex to  middle of the Y arm distance       */
    detector->handle->localTime = localTime;    /* local time - UTC time (second)                */
    return detector;
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrDetectorQueryName_FrameL_(const LALFrameUFrDetector *
    detector)
{
    return detector->handle->name;
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrDetectorQueryPrefix_FrameL_(const LALFrameUFrDetector *
    detector)
{
    return detector->prefix;
}

double XLALFrameUFrDetectorQueryLongitude_FrameL_(const LALFrameUFrDetector *
    detector)
{
    return detector->handle->longitude;
}

double XLALFrameUFrDetectorQueryLatitude_FrameL_(const LALFrameUFrDetector * detector)
{
    return detector->handle->latitude;
}

double XLALFrameUFrDetectorQueryElevation_FrameL_(const LALFrameUFrDetector *
    detector)
{
    return detector->handle->elevation;
}

double XLALFrameUFrDetectorQueryArmXAzimuth_FrameL_(const LALFrameUFrDetector *
    detector)
{
    return detector->handle->armXazimuth;
}

double XLALFrameUFrDetectorQueryArmYAzimuth_FrameL_(const LALFrameUFrDetector *
    detector)
{
    return detector->handle->armYazimuth;
}

double XLALFrameUFrDetectorQueryArmXAltitude_FrameL_(const LALFrameUFrDetector *
    detector)
{
    return detector->handle->armXaltitude;
}

double XLALFrameUFrDetectorQueryArmYAltitude_FrameL_(const LALFrameUFrDetector *
    detector)
{
    return detector->handle->armYaltitude;
}

double XLALFrameUFrDetectorQueryArmXMidpoint_FrameL_(const LALFrameUFrDetector *
    detector)
{
    return detector->handle->armXmidpoint;
}

double XLALFrameUFrDetectorQueryArmYMidpoint_FrameL_(const LALFrameUFrDetector *
    detector)
{
    return detector->handle->armYmidpoint;
}

int XLALFrameUFrDetectorQueryLocalTime_FrameL_(const LALFrameUFrDetector * detector)
{
    return detector->handle->localTime;
}

/*
 * FrHistory routines
 */

void XLALFrameUFrHistoryFree_FrameL_(LALFrameUFrHistory * history)
{
    FrHistoryFree(history);
    return;
}

LALFrameUFrHistory *XLALFrameUFrHistoryAlloc_FrameL_(const char *name, double gpssec,
    const char *comment)
{
    FrHistory *history;
    history = calloc(1, sizeof(*history));
    if (!history)
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    history->classe = FrHistoryDef();
    history->time = floor(gpssec);
    if (name && !(history->name = strdup(name))) {
        XLALFrameUFrHistoryFree(history);
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    }
    if (!(history->comment = strdup(comment))) {
        XLALFrameUFrHistoryFree(history);
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    }
    return history;
}
