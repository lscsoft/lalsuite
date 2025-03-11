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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#define _GNU_SOURCE   /* for mkstemp() */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/Date.h>

#ifndef P_tmpdir
#define P_tmpdir "/tmp"
#endif

/* suppress warnings from FrameC headers; remove this when headers are fixed */
#if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 402
#define GCC_DIAG_STR(s) #s
#define GCC_DIAG_JOINSTR(x,y) GCC_DIAG_STR(x ## y)
# define GCC_DIAG_DO_PRAGMA(x) _Pragma (#x)
# define GCC_DIAG_PRAGMA(x) GCC_DIAG_DO_PRAGMA(GCC diagnostic x)
# if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406
#  define GCC_DIAG_OFF(x) GCC_DIAG_PRAGMA(push) \
	GCC_DIAG_PRAGMA(ignored GCC_DIAG_JOINSTR(-W,x))
#  define GCC_DIAG_ON(x) GCC_DIAG_PRAGMA(pop)
# else
#  define GCC_DIAG_OFF(x) GCC_DIAG_PRAGMA(ignored GCC_DIAG_JOINSTR(-W,x))
#  define GCC_DIAG_ON(x)  GCC_DIAG_PRAGMA(warning GCC_DIAG_JOINSTR(-W,x))
# endif
#else
# define GCC_DIAG_OFF(x)
# define GCC_DIAG_ON(x)
#endif

GCC_DIAG_OFF(strict-prototypes)

#ifdef CHAR
#undef CHAR
#endif
#define CHAR FRAMEC_CHAR
#include <framecppc/FrameC.h>
#include <framecppc/FrameH.h>
#include <framecppc/Stream.h>
#include <framecppc/FrTOC.h>
#undef CHAR
#define CHAR CHAR

GCC_DIAG_ON(strict-prototypes)

/* resolve incomplete types with FrameC structs and custom structs */
#define tagLALFrameUFrameH		frame_h
#define tagLALFrameUFrChan		fr_chan
#define tagLALFrameUFrTOC		fr_toc_t_
#define tagLALFrameUFrHistory		fr_history

#ifndef FRAMEC_HANDLES_STDIO_STDOUT
struct tagLALFrameUFrFile {
    fr_file_t *handle;
    fr_file_mode_t mode;
    FILE *tmpfp;
};
#define FRFILE(stream) ((stream)->handle)
#else
#define tagLALFrameUFrFile		fr_file
#define FRFILE(stream) (stream)
#endif

#include <lal/LALFrameU.h>
#include "LALFrameUFrameC.h"

/* for safety, have a string prefix that is nul-terminated as well */
struct tagLALFrameUFrDetector {
    fr_detector_t *handle;
    char prefix[3];
};

/* calls a FrameC function and assigns the return value to retval */
#define CALL_FRAMEC_FUNCTION_RETVAL(retval, err, f, ...) \
    do { \
        FrameCError *error = NULL; \
        retval = f(&error, __VA_ARGS__); \
        if (error) { \
            XLAL_PRINT_ERROR("FrameC function %s() failed with code %d: %s", #f, (error)->s_errno, (error)->s_message); \
            free(error); \
            err = 1; \
        } \
        else \
            err = 0; \
    } while (0)

/* calls a FrameC function having no return value */
#define CALL_FRAMEC_FUNCTION(err, f, ...) \
    do { \
        FrameCError *error = NULL; \
        f(&error, __VA_ARGS__); \
        if (error) { \
            XLAL_PRINT_ERROR("FrameC function %s() failed with code %d: %s", #f, (error)->s_errno, (error)->s_message); \
            free(error); \
            err = 1; \
        } \
        else \
            err = 0; \
    } while (0)

/* calls a FrameC and returns with 0 on success or -1 on failure */
#define TRY_FRAMEC_FUNCTION(f, ...) \
    do { \
        int err; \
        CALL_FRAMEC_FUNCTION(err, f, __VA_ARGS__); \
        if (err) { \
            XLAL_ERROR(XLAL_EFUNC); \
        } \
        return 0; \
    } while (0)

/* calls a FrameC and returns a specified retval or errval on failure */
#define TRY_FRAMEC_FUNCTION_VAL(retval, errval, f, ...) \
    do { \
        int err; \
        CALL_FRAMEC_FUNCTION(err, f, __VA_ARGS__); \
        if (err) { \
            XLAL_ERROR_VAL(errval, XLAL_EFUNC); \
        } \
        return retval; \
    } while (0)

/* calls a FrameC and returns its integer return value or -1 on failure */
#define TRY_FRAMEC_FUNCTION_INT(f, ...) \
    do { \
        int retval; \
        int err; \
        CALL_FRAMEC_FUNCTION_RETVAL(retval, err, f, __VA_ARGS__); \
        if (err) { \
            XLAL_ERROR(XLAL_EFUNC); \
        } \
        return retval; \
    } while (0)

/* calls a FrameC and returns its pointer return value or NULL on failure */
#define TRY_FRAMEC_FUNCTION_PTR(f, ...) \
    do { \
        void *retval; \
        int err; \
        CALL_FRAMEC_FUNCTION_RETVAL(retval, err, f, __VA_ARGS__); \
        if (err) { \
            XLAL_ERROR_NULL(XLAL_EFUNC); \
        } \
        return retval; \
    } while (0)

/* calls a FrameC and returns */
#define TRY_FRAMEC_FUNCTION_VOID(f, ...) \
    do { \
        int err; \
        CALL_FRAMEC_FUNCTION(err, f, __VA_ARGS__); \
        if (err) { \
            XLAL_ERROR_VOID(XLAL_EFUNC); \
        } \
        return; \
    } while (0)

/*
 * FrFile functions
 */

#ifndef FRAMEC_HANDLES_STDIO_STDOUT

#define BUFSZ 16384

void XLALFrameUFrFileClose_FrameC_(LALFrameUFrFile * stream)
{
    int err;
    if (stream) {
        CALL_FRAMEC_FUNCTION(err, FrameCFileClose, FRFILE(stream));
        if (stream->tmpfp) {
            /* must dump temporary file contents to stdout */
            char buf[BUFSZ];
            while (fwrite(buf, 1, fread(buf, 1, sizeof(buf), stream->tmpfp), stdout) == sizeof(buf)) ;
            fclose(stream->tmpfp);
        }
        LALFree(stream);
        if (err)
            XLAL_ERROR_VOID(XLAL_EFUNC);
    }
    return;
}

LALFrameUFrFile *XLALFrameUFrFileOpen_FrameC_(const char *filename, const char *mode)
{
    char tmpfname[L_tmpnam + 1] = "";
    LALFrameUFrFile *stream;
    int err;
    stream = LALCalloc(1, sizeof(*stream));
    if (*mode == 'r')
        stream->mode = FRAMEC_FILE_MODE_INPUT;
    else if (*mode == 'w')
        stream->mode = FRAMEC_FILE_MODE_OUTPUT;
    else {
        LALFree(stream);
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }
    if (filename == NULL || strcmp(filename, "-") == 0) {
        /* hack to allow reading from stdin */
        filename = tmpfname;
        snprintf(tmpfname, L_tmpnam, "%s/tmp.lal.XXXXXX", P_tmpdir);
        if (stream->mode == FRAMEC_FILE_MODE_INPUT) {
            /* dump stdin to a temporary file */
            char buf[BUFSZ];
            stream->tmpfp = fdopen(mkstemp(tmpfname), "w");
            while (fwrite(buf, 1, fread(buf, 1, sizeof(buf), stdin), stream->tmpfp) == sizeof(buf)) ;
            fclose(stream->tmpfp);
            stream->tmpfp = NULL;
        } else {
            stream->tmpfp = fdopen(mkstemp(tmpfname), "r");
        }
    }
    CALL_FRAMEC_FUNCTION_RETVAL(stream->handle, err, FrameCFileOpen, filename, stream->mode);
    if (*tmpfname) {
        unlink(tmpfname);       /* remove temporary file when closed */
    }
    if (err) {
        LALFree(stream);
        XLAL_ERROR_NULL(XLAL_EFUNC);
    }
    return stream;
}

int XLALFrameUFileCksumValid_FrameC_(LALFrameUFrFile * stream)
{
    TRY_FRAMEC_FUNCTION_INT(FrameCFileCksumValid, FRFILE(stream));
}

#undef BUFSZ

#else /* defined FRAMEC_HANDLES_STDIN_STDOUT */

void XLALFrameUFrFileClose_FrameC_(LALFrameUFrFile * stream)
{
    TRY_FRAMEC_FUNCTION_VOID(FrameCFileClose, stream);
}

LALFrameUFrFile *XLALFrameUFrFileOpen_FrameC_(const char *filename, const char *mode)
{
    fr_file_mode_t m;
    if (*mode == 'r')
        m = FRAMEC_FILE_MODE_INPUT;
    else if (*mode == 'w')
        m = FRAMEC_FILE_MODE_OUTPUT;
    else
        return NULL;
    /* FIXME: does this work? */
    if (filename == NULL)
        filename = "-"; /* indicates we should use stdin or stdout */
    TRY_FRAMEC_FUNCTION_PTR(FrameCFileOpen, filename, m);
}

int XLALFrameUFileCksumValid_FrameC_(LALFrameUFrFile * stream)
{
    TRY_FRAMEC_FUNCTION_INT(FrameCFileCksumValid, stream);
}

#endif /* defined FRAMEC_HANDLES_STDIN_STDOUT */

/*
 * FrTOC functions
 */

void XLALFrameUFrTOCFree_FrameC_(LALFrameUFrTOC * toc)
{
    TRY_FRAMEC_FUNCTION_VOID(FrameCFrTOCFree, toc);
}

LALFrameUFrTOC *XLALFrameUFrTOCRead_FrameC_(LALFrameUFrFile * stream)
{
    TRY_FRAMEC_FUNCTION_PTR(FrameCFrTOCRead, FRFILE(stream));
}

size_t XLALFrameUFrTOCQueryNFrame_FrameC_(const LALFrameUFrTOC * toc)
{
    fr_toc_nframe_t nframe;
    TRY_FRAMEC_FUNCTION_VAL(nframe, -1, FrameCFrTOCQuery, toc, FR_TOC_FIELD_NFRAME, &nframe, FR_TOC_FIELD_LAST);
}

double XLALFrameUFrTOCQueryGTimeModf_FrameC_(double *iptr, const LALFrameUFrTOC * toc, size_t pos)
{
    fr_toc_nframe_t nframe;
    fr_toc_t0_t *gtime;
    LIGOTimeGPS epoch;
    int err;
    CALL_FRAMEC_FUNCTION(err, FrameCFrTOCQuery, toc, FR_TOC_FIELD_NFRAME, &nframe, FR_TOC_FIELD_LAST);
    if (err || nframe <= pos)
        XLAL_ERROR_REAL8(XLAL_EINVAL);
    gtime = LALCalloc(nframe, sizeof(*gtime));
    if (!gtime)
        XLAL_ERROR_REAL8(XLAL_EINVAL);
    CALL_FRAMEC_FUNCTION(err, FrameCFrTOCQuery, toc, FR_TOC_FIELD_GTIME, gtime, nframe, FR_TOC_FIELD_LAST);
    if (err) {
        LALFree(gtime);
        XLAL_ERROR_REAL8(XLAL_EINVAL);
    }
    XLALGPSSet(&epoch, gtime[pos].sec, gtime[pos].nan);
    LALFree(gtime);
    return XLALGPSModf(iptr, &epoch);
}

double XLALFrameUFrTOCQueryDt_FrameC_(const LALFrameUFrTOC * toc, size_t pos)
{
    fr_toc_nframe_t nframe;
    fr_toc_dt_t dt;
    double retval;
    int err;
    CALL_FRAMEC_FUNCTION(err, FrameCFrTOCQuery, toc, FR_TOC_FIELD_NFRAME, &nframe, FR_TOC_FIELD_LAST);
    if (err || nframe <= pos)
        return -1.0;
    dt = LALCalloc(nframe, sizeof(*dt));
    if (!dt)
        return -1.0;
    CALL_FRAMEC_FUNCTION(err, FrameCFrTOCQuery, toc, FR_TOC_FIELD_DT, dt, nframe, FR_TOC_FIELD_LAST);
    retval = err ? -1.0 : dt[pos];
    LALFree(dt);
    return retval;
}

size_t XLALFrameUFrTOCQueryAdcN_FrameC_(const LALFrameUFrTOC * toc)
{
    fr_toc_adc_n_t n;
    TRY_FRAMEC_FUNCTION_VAL(n, -1, FrameCFrTOCQuery, toc, FR_TOC_FIELD_ADC_N, &n, FR_TOC_FIELD_LAST);
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrTOCQueryAdcName_FrameC_(const LALFrameUFrTOC * toc, size_t adc)
{
    fr_toc_adc_name_t *names;
    fr_toc_adc_n_t n;
    const char *name;
    int err;
    CALL_FRAMEC_FUNCTION(err, FrameCFrTOCQuery, toc, FR_TOC_FIELD_ADC_N, &n, FR_TOC_FIELD_LAST);
    if (err || n <= adc)
        return NULL;
    names = LALCalloc(n, sizeof(*names));
    if (!names)
        return NULL;
    CALL_FRAMEC_FUNCTION(err, FrameCFrTOCQuery, toc, FR_TOC_FIELD_ADC_NAMES, names, n, FR_TOC_FIELD_LAST);
    name = err ? NULL : names[adc];
    LALFree(names);
    return name;
}

size_t XLALFrameUFrTOCQuerySimN_FrameC_(const LALFrameUFrTOC * toc)
{
    fr_toc_sim_n_t n;
    TRY_FRAMEC_FUNCTION_VAL(n, -1, FrameCFrTOCQuery, toc, FR_TOC_FIELD_SIM_N, &n, FR_TOC_FIELD_LAST);
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrTOCQuerySimName_FrameC_(const LALFrameUFrTOC * toc, size_t sim)
{
    fr_toc_sim_name_t *names;
    fr_toc_sim_n_t n;
    const char *name;
    int err;
    CALL_FRAMEC_FUNCTION(err, FrameCFrTOCQuery, toc, FR_TOC_FIELD_SIM_N, &n, FR_TOC_FIELD_LAST);
    if (err || n <= sim)
        return NULL;
    names = LALCalloc(n, sizeof(*names));
    if (!names)
        return NULL;
    CALL_FRAMEC_FUNCTION(err, FrameCFrTOCQuery, toc, FR_TOC_FIELD_SIM_NAMES, names, n, FR_TOC_FIELD_LAST);
    name = err ? NULL : names[sim];
    LALFree(names);
    return name;
}

size_t XLALFrameUFrTOCQueryProcN_FrameC_(const LALFrameUFrTOC * toc)
{
    fr_toc_proc_n_t n;
    TRY_FRAMEC_FUNCTION_VAL(n, -1, FrameCFrTOCQuery, toc, FR_TOC_FIELD_PROC_N, &n, FR_TOC_FIELD_LAST);
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrTOCQueryProcName_FrameC_(const LALFrameUFrTOC * toc, size_t proc)
{
    fr_toc_proc_name_t *names;
    fr_toc_proc_n_t n;
    const char *name;
    int err;
    CALL_FRAMEC_FUNCTION(err, FrameCFrTOCQuery, toc, FR_TOC_FIELD_PROC_N, &n, FR_TOC_FIELD_LAST);
    if (err || n <= proc)
        return NULL;
    names = LALCalloc(n, sizeof(*names));
    if (!names)
        return NULL;
    CALL_FRAMEC_FUNCTION(err, FrameCFrTOCQuery, toc, FR_TOC_FIELD_PROC_NAMES, names, n, FR_TOC_FIELD_LAST);
    name = err ? NULL : names[proc];
    LALFree(names);
    return name;
}

size_t XLALFrameUFrTOCQueryDetectorN_FrameC_(const LALFrameUFrTOC * toc)
{
    fr_toc_detector_n_t n;
    TRY_FRAMEC_FUNCTION_VAL(n, -1, FrameCFrTOCQuery, toc, FR_TOC_FIELD_DETECTOR_N, &n, FR_TOC_FIELD_LAST);
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrTOCQueryDetectorName_FrameC_(const LALFrameUFrTOC * toc, size_t det)
{
    fr_toc_detector_name_t *names;
    fr_toc_detector_n_t n;
    const char *name;
    int err;
    CALL_FRAMEC_FUNCTION(err, FrameCFrTOCQuery, toc, FR_TOC_FIELD_DETECTOR_N, &n, FR_TOC_FIELD_LAST);
    if (err || n <= det)
        return NULL;
    names = LALCalloc(n, sizeof(*names));
    if (!names)
        return NULL;
    CALL_FRAMEC_FUNCTION(err, FrameCFrTOCQuery, toc, FR_TOC_FIELD_DETECTOR_NAMES, names, n, FR_TOC_FIELD_LAST);
    name = err ? NULL : names[det];
    LALFree(names);
    return name;
}

/*
 * FrameH functions
 */

void XLALFrameUFrameHFree_FrameC_(LALFrameUFrameH * frame)
{
    TRY_FRAMEC_FUNCTION_VOID(FrameCFrameHFree, frame);
}

LALFrameUFrameH *XLALFrameUFrameHAlloc_FrameC_(const char *name, double start1, double start2, double dt, int frnum)
{
    double fp, fp1, fp2;
    double ip, ip1, ip2;
    frame_h_gtime_t gtime;

    /* break start time into integer and fractional parts */
    fp1 = modf(start1, &ip1);
    fp2 = modf(start2, &ip2);
    fp = modf(fp1 + fp2, &ip);
    ip += ip1 + ip2;
    if (fp < 0.0) { /* make sure fractional part is positive */
        fp += 1.0;
        ip -= 1.0;
    }

    gtime.sec = ip;
    gtime.nan = floor(0.5 + 1e9 * fp);
    if (gtime.nan >= 1000000000) { /* handle round-up corner case */
        gtime.nan -= 1000000000;
        gtime.sec += 1;
    }

    TRY_FRAMEC_FUNCTION_PTR(FrameCFrameHAlloc, name, gtime, dt, frnum);
}

LALFrameUFrameH *XLALFrameUFrameHRead_FrameC_(LALFrameUFrFile * stream, int pos)
{
    TRY_FRAMEC_FUNCTION_PTR(FrameCFrameHRead, FRFILE(stream), pos);
}

int XLALFrameUFrameHWrite_FrameC_(LALFrameUFrFile * stream, LALFrameUFrameH * frame)
{
    TRY_FRAMEC_FUNCTION(FrameCFrameHWrite, FRFILE(stream), frame);
}

/* function to add a channel to a frame */

int XLALFrameUFrameHFrChanAdd_FrameC_(LALFrameUFrameH * frame, LALFrameUFrChan * channel)
{
    TRY_FRAMEC_FUNCTION(FrameCFrameHFrChanAdd, frame, channel);
}

/* function to add a detector to a frame */

int XLALFrameUFrameHFrDetectorAdd_FrameC_(LALFrameUFrameH * frame, LALFrameUFrDetector * detector)
{
    TRY_FRAMEC_FUNCTION(FrameCFrameHFrDetectorAdd, frame, detector->handle);
}

/* function to add a detector to a frame */

int XLALFrameUFrameHFrHistoryAdd_FrameC_(LALFrameUFrameH * frame, LALFrameUFrHistory * history)
{
    TRY_FRAMEC_FUNCTION(FrameCFrameHFrHistoryAdd, frame, history);
}

/* functions to query frame metadata */

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrameHQueryName_FrameC_(const LALFrameUFrameH * frame)
{
    frame_h_name_t name;
    TRY_FRAMEC_FUNCTION_VAL(name, NULL, FrameCFrameHQuery, frame, FRAME_H_FIELD_NAME, &name, FRAME_H_FIELD_LAST);
}

int XLALFrameUFrameHQueryRun_FrameC_(const LALFrameUFrameH * frame)
{
    frame_h_run_t run;
    TRY_FRAMEC_FUNCTION_VAL(run, -1, FrameCFrameHQuery, frame, FRAME_H_FIELD_RUN, &run, FRAME_H_FIELD_LAST);
}

int XLALFrameUFrameHQueryFrame_FrameC_(const LALFrameUFrameH * frame)
{
    frame_h_frame_t frnum;
    TRY_FRAMEC_FUNCTION_VAL(frnum, -1, FrameCFrameHQuery, frame, FRAME_H_FIELD_FRAME, &frnum, FRAME_H_FIELD_LAST);
}

int XLALFrameUFrameHQueryDataQuality_FrameC_(const LALFrameUFrameH * frame)
{
    frame_h_data_quality_t dq;
    TRY_FRAMEC_FUNCTION_VAL(dq, -1, FrameCFrameHQuery, frame, FRAME_H_FIELD_DATA_QUALITY, &dq, FRAME_H_FIELD_LAST);
}

double XLALFrameUFrameHQueryGTimeModf_FrameC_(double *iptr, const LALFrameUFrameH * frame)
{
    frame_h_gtime_t start;
    LIGOTimeGPS epoch;
    TRY_FRAMEC_FUNCTION_VAL(XLALGPSModf(iptr, XLALGPSSet(&epoch, start.sec,
                start.nan)), XLAL_REAL8_FAIL_NAN, FrameCFrameHQuery, frame, FRAME_H_FIELD_GTIME, &start, FRAME_H_FIELD_LAST);
}

int XLALFrameUFrameHQueryULeapS_FrameC_(const LALFrameUFrameH * frame)
{
    frame_h_uleaps_t uleaps;
    TRY_FRAMEC_FUNCTION_VAL(uleaps, -1, FrameCFrameHQuery, frame, FRAME_H_FIELD_ULEAPS, &uleaps, FRAME_H_FIELD_LAST);
}

double XLALFrameUFrameHQueryDt_FrameC_(const LALFrameUFrameH * frame)
{
    frame_h_dt_t dt;
    TRY_FRAMEC_FUNCTION_VAL(dt, -1, FrameCFrameHQuery, frame, FRAME_H_FIELD_DT, &dt, FRAME_H_FIELD_LAST);
}

/* functions to set frame metadata */

int XLALFrameUFrameHSetRun_FrameC_(LALFrameUFrameH * frame, int run)
{
    frame_h_run_t run_ = run;
    TRY_FRAMEC_FUNCTION(FrameCFrameHSet, frame, FRAME_H_FIELD_RUN, run_, FRAME_H_FIELD_LAST);
}

/*
 * FrChan functions
 */

void XLALFrameUFrChanFree_FrameC_(LALFrameUFrChan * channel)
{
    TRY_FRAMEC_FUNCTION_VOID(FrameCFrChanFree, channel);
}

LALFrameUFrChan *XLALFrameUFrChanRead_FrameC_(LALFrameUFrFile * stream, const char *name, size_t pos)
{
    TRY_FRAMEC_FUNCTION_PTR(FrameCFrChanRead, FRFILE(stream), name, pos);
}

LALFrameUFrChan *XLALFrameUFrAdcChanAlloc_FrameC_(const char *name, int dtype, size_t ndata)
{
    TRY_FRAMEC_FUNCTION_PTR(FrameCFrChanAlloc, name, FR_ADC_CHAN_TYPE, dtype, ndata);
}

LALFrameUFrChan *XLALFrameUFrSimChanAlloc_FrameC_(const char *name, int dtype, size_t ndata)
{
    TRY_FRAMEC_FUNCTION_PTR(FrameCFrChanAlloc, name, FR_SIM_CHAN_TYPE, dtype, ndata);
}

static fr_proc_type_t XLALFrProcType(int type)
{
    switch (type) {
    case 0:
        return FR_PROC_TYPE_UNKNOWN;
    case 1:
        return FR_PROC_TYPE_TIME_SERIES;
    case 2:
        return FR_PROC_TYPE_FREQUENCY_SERIES;
    case 3:
        return FR_PROC_TYPE_OTHER_1D_SERIES;
    case 4:
        return FR_PROC_TYPE_TIME_FREQUENCY;
    case 5:
        return FR_PROC_TYPE_WAVELET;
    case 6:
        /* NOTE: lower case r in Fr */
        return Fr_PROC_TYPE_MULTI_DIMENSIONAL;
    default:
        return -1;
    }
}

static fr_proc_sub_type_t XLALFrProcSubType(int subtype)
{
    switch (subtype) {
    case 0:
        return FR_PROC_SUB_TYPE_UNKNOWN;
    case 1:
        return FR_PROC_SUB_TYPE_DFT;
    case 2:
        return FR_PROC_SUB_TYPE_AMPLITUDE_SPECTRAL_DENSITY;
    case 3:
        return FR_PROC_SUB_TYPE_POWER_SPECTRAL_DENSITY;
    case 4:
        return FR_PROC_SUB_TYPE_CROSS_SPECTRAL_DENSITY;
    case 5:
        return FR_PROC_SUB_TYPE_COHERENCE;
    case 6:
        return FR_PROC_SUB_TYPE_TRANSFER_FUNCTION;
    default:
        return -1;
    }
}

LALFrameUFrChan *XLALFrameUFrProcChanAlloc_FrameC_(const char *name, int type, int subtype, int dtype, size_t ndata)
{
    TRY_FRAMEC_FUNCTION_PTR(FrameCFrProcChanAlloc, name, XLALFrProcType(type), XLALFrProcSubType(subtype), dtype, ndata);
}

/* channel query functions */

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrChanQueryName_FrameC_(const LALFrameUFrChan * channel)
{
    fr_chan_name_t name;
    TRY_FRAMEC_FUNCTION_VAL(name, NULL, FrameCFrChanQuery, channel, FR_CHAN_FIELD_NAME, &name, FR_CHAN_FIELD_LAST);
}

double XLALFrameUFrChanQueryTimeOffset_FrameC_(const LALFrameUFrChan * channel)
{
    fr_chan_time_offset_t offset;
    TRY_FRAMEC_FUNCTION_VAL(offset, 0, FrameCFrChanQuery, channel, FR_CHAN_FIELD_TIME_OFFSET, &offset, FR_CHAN_FIELD_LAST);
}

/* channel set functions */

int XLALFrameUFrChanSetSampleRate_FrameC_(LALFrameUFrChan * channel, double sampleRate)
{
    fr_chan_sample_rate_t srate = sampleRate;
    TRY_FRAMEC_FUNCTION(FrameCFrChanSet, channel, FR_CHAN_FIELD_SAMPLE_RATE, srate, FR_CHAN_FIELD_LAST);
}

int XLALFrameUFrChanSetTimeOffset_FrameC_(LALFrameUFrChan * channel, double timeOffset)
{
    fr_chan_time_offset_t offset = timeOffset;
    TRY_FRAMEC_FUNCTION(FrameCFrChanSet, channel, FR_CHAN_FIELD_TIME_OFFSET, offset, FR_CHAN_FIELD_LAST);
}

int XLALFrameUFrChanSetTRange_FrameC_(LALFrameUFrChan * channel, double tRange)
{
    fr_chan_t_range_t dt = tRange;
    TRY_FRAMEC_FUNCTION(FrameCFrChanSet, channel, FR_CHAN_FIELD_T_RANGE, dt, FR_CHAN_FIELD_LAST);
}

/*
 * FrVect functions
 */

int XLALFrameUFrChanVectorAlloc_FrameC_(LALFrameUFrChan * channel, int dtype, size_t ndata)
{
    TRY_FRAMEC_FUNCTION(FrameCFrChanVectorAlloc, channel, dtype, ndata);
}

/* get the FrameC compression scheme for a particular compressLevel */
/* UNUSED */
/*
static fr_vect_compression_schemes_t XLALFrVectCompressionScheme(int compressLevel)
{
    switch (compressLevel) {
    case 0:
    case 256:
        return FR_VECT_COMPRESS_RAW;
    case 1:
    case 257:
        return FR_VECT_COMPRESS_GZIP;
    case 3:
    case 259:
        return FR_VECT_COMPRESS_DIFF_GZIP;
    case 5:
    case 261:
        return FR_VECT_COMPRESS_ZERO_SUPPRESS_WORD_2;
    case 8:
    case 264:
        return FR_VECT_COMPRESS_ZERO_SUPPRESS_WORD_4;
    default:
        return FR_VECT_COMPRESS_UNKNOWN;
    }
}
*/

/* functions to compress and expand FrVect structures within a channel */

int XLALFrameUFrChanVectorCompress_FrameC_(LALFrameUFrChan * channel, int compressLevel)
{
    /* Work around bug in FrameC FrameCFrChanVectorCompress() by disabling
     * compression.  Revert this once FrameCFrChanVectorCompress() works. */
    (void)channel;
    (void)compressLevel;
    XLAL_PRINT_WARNING("Compression not currently implemented with FrameC");
    return 0;
    /*
    TRY_FRAMEC_FUNCTION(FrameCFrChanVectorCompress, channel, XLALFrVectCompressionScheme(compressLevel));
    */
}

int XLALFrameUFrChanVectorExpand_FrameC_(LALFrameUFrChan * channel)
{
    TRY_FRAMEC_FUNCTION(FrameCFrChanVectorExpand, channel);
}

/* functions to query FrVect structures within a channel */

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrChanVectorQueryName_FrameC_(const LALFrameUFrChan * channel)
{
    fr_vect_name_t name;
    TRY_FRAMEC_FUNCTION_VAL(name, NULL, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_NAME, &name, FR_VECT_FIELD_LAST);
}

int XLALFrameUFrChanVectorQueryCompress_FrameC_(const LALFrameUFrChan * channel)
{
    fr_vect_compress_t compress;
    TRY_FRAMEC_FUNCTION_VAL(compress, -1, FrameCFrChanVectorQuery, channel,
        FR_VECT_FIELD_COMPRESS, &compress, FR_VECT_FIELD_LAST);
}

int XLALFrameUFrChanVectorQueryType_FrameC_(const LALFrameUFrChan * channel)
{
    fr_vect_type_t type;
    TRY_FRAMEC_FUNCTION_VAL(type, -1, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_TYPE, &type, FR_VECT_FIELD_LAST);
}

/* retrieves a handle to the data vector in the FrVect structure */
void *XLALFrameUFrChanVectorQueryData_FrameC_(const LALFrameUFrChan * channel)
{
    fr_vect_data_t data;
    TRY_FRAMEC_FUNCTION_VAL(data, NULL, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_DATA, &data, FR_VECT_FIELD_LAST);
}

size_t XLALFrameUFrChanVectorQueryNBytes_FrameC_(const LALFrameUFrChan * channel)
{
    fr_vect_nbytes_t nbytes;
    TRY_FRAMEC_FUNCTION_VAL(nbytes, -1, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_NBYTES, &nbytes, FR_VECT_FIELD_LAST);
}

size_t XLALFrameUFrChanVectorQueryNData_FrameC_(const LALFrameUFrChan * channel)
{
    fr_vect_ndata_t ndata;
    TRY_FRAMEC_FUNCTION_VAL(ndata, -1, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_NDATA, &ndata, FR_VECT_FIELD_LAST);
}

size_t XLALFrameUFrChanVectorQueryNDim_FrameC_(const LALFrameUFrChan * channel)
{
    fr_vect_ndim_t ndim;
    TRY_FRAMEC_FUNCTION_VAL(ndim, -1, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_NDIM, &ndim, FR_VECT_FIELD_LAST);
}

size_t XLALFrameUFrChanVectorQueryNx_FrameC_(const LALFrameUFrChan * channel, size_t dim)
{
    fr_vect_ndim_t ndim;
    fr_vect_nx_t *nx;
    int err;
    size_t retval;
    CALL_FRAMEC_FUNCTION(err, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_NDIM, &ndim, FR_VECT_FIELD_LAST);
    if (err || ndim <= dim)
        return -1;
    nx = LALCalloc(ndim, sizeof(*nx));
    if (!nx)
        return -1;
    CALL_FRAMEC_FUNCTION(err, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_NX, nx, ndim, FR_VECT_FIELD_LAST);
    retval = err ? (size_t) (-1) : nx[dim];
    LALFree(nx);
    return retval;
}

double XLALFrameUFrChanVectorQueryDx_FrameC_(const LALFrameUFrChan * channel, size_t dim)
{
    fr_vect_ndim_t ndim;
    fr_vect_dx_t *dx;
    int err;
    double retval;
    CALL_FRAMEC_FUNCTION(err, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_NDIM, &ndim, FR_VECT_FIELD_LAST);
    if (err || ndim <= dim)
        return -1.0;
    dx = LALCalloc(ndim, sizeof(*dx));
    if (!dx)
        return -1.0;
    CALL_FRAMEC_FUNCTION(err, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_DX, dx, ndim, FR_VECT_FIELD_LAST);
    retval = err ? -1.0 : dx[dim];
    LALFree(dx);
    return retval;
}

double XLALFrameUFrChanVectorQueryStartX_FrameC_(const LALFrameUFrChan * channel, size_t dim)
{
    fr_vect_ndim_t ndim;
    fr_vect_startx_t *startX;
    int err;
    double retval;
    CALL_FRAMEC_FUNCTION(err, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_NDIM, &ndim, FR_VECT_FIELD_LAST);
    if (err || ndim <= dim)
        return -1.0;
    startX = LALCalloc(ndim, sizeof(*startX));
    if (!startX)
        return -1.0;
    CALL_FRAMEC_FUNCTION(err, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_START_X, startX, ndim, FR_VECT_FIELD_LAST);
    retval = err ? -1.0 : startX[dim];
    LALFree(startX);
    return retval;
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrChanVectorQueryUnitX_FrameC_(const LALFrameUFrChan * channel, size_t dim)
{
    fr_vect_ndim_t ndim;
    fr_vect_unit_x_t *unitX;
    int err;
    const char *retval;
    CALL_FRAMEC_FUNCTION(err, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_NDIM, &ndim, FR_VECT_FIELD_LAST);
    if (err || ndim <= dim)
        return NULL;
    unitX = LALCalloc(ndim, sizeof(*unitX));
    if (!unitX)
        return NULL;
    CALL_FRAMEC_FUNCTION(err, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_UNIT_X, unitX, ndim, FR_VECT_FIELD_LAST);
    retval = err ? NULL : unitX[dim];
    LALFree(unitX);
    return retval;
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrChanVectorQueryUnitY_FrameC_(const LALFrameUFrChan * channel)
{
    fr_vect_unit_y_t unitY;
    TRY_FRAMEC_FUNCTION_VAL(unitY, NULL, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_UNIT_Y, &unitY, FR_VECT_FIELD_LAST);
}

/* functions to set FrVect structures within a channel */

int XLALFrameUFrChanVectorSetName_FrameC_(LALFrameUFrChan * channel, const char *name)
{
    TRY_FRAMEC_FUNCTION(FrameCFrChanVectorSet, channel, FR_VECT_FIELD_NAME, (fr_vect_name_t) (name), FR_VECT_FIELD_LAST);
}

int XLALFrameUFrChanVectorSetUnitY_FrameC_(LALFrameUFrChan * channel, const char *unit)
{
    TRY_FRAMEC_FUNCTION(FrameCFrChanVectorSet, channel, FR_VECT_FIELD_UNIT_Y,
        (fr_vect_unit_y_t) (unit), 1, FR_VECT_FIELD_LAST);
}

/* NOTE: only support 1-dimensional vectors */

int XLALFrameUFrChanVectorSetDx_FrameC_(LALFrameUFrChan * channel, double dx)
{
    fr_vect_dx_t dx_ = dx;
    TRY_FRAMEC_FUNCTION(FrameCFrChanVectorSet, channel, FR_VECT_FIELD_DX, &dx_, 1, FR_VECT_FIELD_LAST);
}

int XLALFrameUFrChanVectorSetStartX_FrameC_(LALFrameUFrChan * channel, double x0)
{
    fr_vect_startx_t x0_ = x0;
    TRY_FRAMEC_FUNCTION(FrameCFrChanVectorSet, channel, FR_VECT_FIELD_START_X, &x0_, 1, FR_VECT_FIELD_LAST);
}

int XLALFrameUFrChanVectorSetUnitX_FrameC_(LALFrameUFrChan * channel, const char *unit)
{
    fr_vect_unit_x_t unit_ = unit;
    TRY_FRAMEC_FUNCTION(FrameCFrChanVectorSet, channel, FR_VECT_FIELD_UNIT_X, &unit_, 1, FR_VECT_FIELD_LAST);
}

/*
 * FrDetector functions
 */

void XLALFrameUFrDetectorFree_FrameC_(LALFrameUFrDetector * detector)
{
    if (detector) {
        int err;
        CALL_FRAMEC_FUNCTION(err, FrameCFrDetectorFree, detector->handle);
        LALFree(detector);
        if (err)
            XLAL_ERROR_VOID(XLAL_EFUNC);
    }
    return;
}

LALFrameUFrDetector *XLALFrameUFrDetectorRead_FrameC_(LALFrameUFrFile * stream, const char *name)
{
    fr_detector_prefix_t prefix;
    LALFrameUFrDetector *detector;
    int err;
    detector = LALCalloc(1, sizeof(*detector));
    if (!detector)
        return NULL;
    CALL_FRAMEC_FUNCTION(err, detector->handle = FrameCFrDetectorRead, FRFILE(stream), name);
    if (err || !detector->handle) {
        XLALFrameUFrDetectorFree_FrameC_(detector);
        return NULL;
    }
    /* read the prefix and use it to set the detector's prefix */
    CALL_FRAMEC_FUNCTION(err, FrameCFrDetectorQuery, detector->handle,
        FR_DETECTOR_FIELD_PREFIX, &prefix, FR_DETECTOR_FIELD_LAST);
    if (err) {
        XLALFrameUFrDetectorFree_FrameC_(detector);
        return NULL;
    }
    memcpy(detector->prefix, prefix.prefix, 2);
    return detector;
}

LALFrameUFrDetector *XLALFrameUFrDetectorAlloc_FrameC_(const char *name,
    const char *prefix, double latitude, double longitude, double elevation,
    double azimuthX, double azimuthY, double altitudeX, double altitudeY, double midpointX, double midpointY, int localTime)
{
    LALFrameUFrDetector *detector;
    int err;
    detector = LALCalloc(1, sizeof(*detector));
    if (!detector)
        return NULL;
    if (prefix)
        memcpy(detector->prefix, prefix, 2);
    CALL_FRAMEC_FUNCTION_RETVAL(detector->handle, err, FrameCFrDetectorAlloc,
        name, prefix, latitude, longitude, elevation, azimuthX, azimuthY,
        altitudeX, altitudeY, midpointX, midpointY, localTime);
    if (err || !detector->handle) {
        XLALFrameUFrDetectorFree_FrameC_(detector);
        return NULL;
    }
    return detector;
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrDetectorQueryName_FrameC_(const LALFrameUFrDetector * detector)
{
    fr_detector_name_t name;
    TRY_FRAMEC_FUNCTION_VAL(name, NULL, FrameCFrDetectorQuery,
        detector->handle, FR_DETECTOR_FIELD_NAME, &name, FR_DETECTOR_FIELD_LAST);
}

/* WARNING: returns pointer to memory that is lost when frame is freed */
const char *XLALFrameUFrDetectorQueryPrefix_FrameC_(const LALFrameUFrDetector * detector)
{
    return detector->prefix;
}

double XLALFrameUFrDetectorQueryLongitude_FrameC_(const LALFrameUFrDetector * detector)
{
    fr_detector_longitude_t longitude;
    TRY_FRAMEC_FUNCTION_VAL(longitude, -1.0, FrameCFrDetectorQuery,
        detector->handle, FR_DETECTOR_FIELD_LONGITUDE, &longitude, FR_DETECTOR_FIELD_LAST);
}

double XLALFrameUFrDetectorQueryLatitude_FrameC_(const LALFrameUFrDetector * detector)
{
    fr_detector_latitude_t latitude;
    TRY_FRAMEC_FUNCTION_VAL(latitude, -1.0, FrameCFrDetectorQuery,
        detector->handle, FR_DETECTOR_FIELD_LATITUDE, &latitude, FR_DETECTOR_FIELD_LAST);
}

double XLALFrameUFrDetectorQueryElevation_FrameC_(const LALFrameUFrDetector * detector)
{
    fr_detector_elevation_t elevation;
    TRY_FRAMEC_FUNCTION_VAL(elevation, -1.0, FrameCFrDetectorQuery,
        detector->handle, FR_DETECTOR_FIELD_ELEVATION, &elevation, FR_DETECTOR_FIELD_LAST);
}

double XLALFrameUFrDetectorQueryArmXAzimuth_FrameC_(const LALFrameUFrDetector * detector)
{
    fr_detector_arm_x_azimuth_t armXazimuth;
    TRY_FRAMEC_FUNCTION_VAL(armXazimuth, -1.0, FrameCFrDetectorQuery,
        detector->handle, FR_DETECTOR_FIELD_ARM_X_AZIMUTH, &armXazimuth, FR_DETECTOR_FIELD_LAST);
}

double XLALFrameUFrDetectorQueryArmYAzimuth_FrameC_(const LALFrameUFrDetector * detector)
{
    fr_detector_arm_y_azimuth_t armYazimuth;
    TRY_FRAMEC_FUNCTION_VAL(armYazimuth, -1.0, FrameCFrDetectorQuery,
        detector->handle, FR_DETECTOR_FIELD_ARM_Y_AZIMUTH, &armYazimuth, FR_DETECTOR_FIELD_LAST);
}

double XLALFrameUFrDetectorQueryArmXAltitude_FrameC_(const LALFrameUFrDetector * detector)
{
    fr_detector_arm_x_altitude_t armXaltitude;
    TRY_FRAMEC_FUNCTION_VAL(armXaltitude, -1.0, FrameCFrDetectorQuery,
        detector->handle, FR_DETECTOR_FIELD_ARM_X_ALTITUDE, &armXaltitude, FR_DETECTOR_FIELD_LAST);
}

double XLALFrameUFrDetectorQueryArmYAltitude_FrameC_(const LALFrameUFrDetector * detector)
{
    fr_detector_arm_y_altitude_t armYaltitude;
    TRY_FRAMEC_FUNCTION_VAL(armYaltitude, -1.0, FrameCFrDetectorQuery,
        detector->handle, FR_DETECTOR_FIELD_ARM_Y_ALTITUDE, &armYaltitude, FR_DETECTOR_FIELD_LAST);
}

double XLALFrameUFrDetectorQueryArmXMidpoint_FrameC_(const LALFrameUFrDetector * detector)
{
    fr_detector_arm_x_midpoint_t armXmidpoint;
    TRY_FRAMEC_FUNCTION_VAL(armXmidpoint, -1.0, FrameCFrDetectorQuery,
        detector->handle, FR_DETECTOR_FIELD_ARM_X_MIDPOINT, &armXmidpoint, FR_DETECTOR_FIELD_LAST);
}

double XLALFrameUFrDetectorQueryArmYMidpoint_FrameC_(const LALFrameUFrDetector * detector)
{
    fr_detector_arm_y_midpoint_t armYmidpoint;
    TRY_FRAMEC_FUNCTION_VAL(armYmidpoint, -1.0, FrameCFrDetectorQuery,
        detector->handle, FR_DETECTOR_FIELD_ARM_Y_MIDPOINT, &armYmidpoint, FR_DETECTOR_FIELD_LAST);
}

int XLALFrameUFrDetectorQueryLocalTime_FrameC_(const LALFrameUFrDetector * detector)
{
    fr_detector_localtime_t localTime;
    TRY_FRAMEC_FUNCTION_VAL(localTime, -1.0, FrameCFrDetectorQuery,
        detector->handle, FR_DETECTOR_FIELD_LOCAL_TIME, &localTime, FR_DETECTOR_FIELD_LAST);
}

/*
 * FrHistory routines
 */

void XLALFrameUFrHistoryFree_FrameC_(LALFrameUFrHistory * history)
{
    TRY_FRAMEC_FUNCTION_VOID(FrameCFrHistoryFree, history);
}

LALFrameUFrHistory *XLALFrameUFrHistoryAlloc_FrameC_(const char *name, double gpssec, const char *comment)
{
    TRY_FRAMEC_FUNCTION_PTR(FrameCFrHistoryAlloc, name, (fr_history_time_t) floor(gpssec), comment);
}
