#include <math.h>
#include <stdio.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/Date.h>

#ifdef CHAR
#undef CHAR
#endif
#define CHAR FRAMEC_CHAR
#include <framec/FrameC.h>
#include <framec/FrameH.h>
#include <framec/Stream.h>
#include <framec/FrTOC.h>
#undef CHAR
#define CHAR CHAR

/* resolve incomplete types with FrameL structs and custom structs */
#define tagLALFrameUFrameH		frame_h
#define tagLALFrameUFrFile		fr_file
#define tagLALFrameUFrChan		fr_chan
#define tagLALFrameUFrTOC		fr_toc_t_
#define tagLALFrameUFrHistory		fr_history
#include <lal/LALFrameU.h>

/* for safety, have a string prefix that is nul-terminated as well */
struct tagLALFrameUFrDetector {
    fr_detector_t *handle;
    char prefix[3];
};

/* TODO: do something with XLAL error macros;
 * TODO: also, free the error handle */

#define CALL_FRAMEC_FUNCTION(err, f, ...) do { \
		FrameCError *error = NULL; \
		f(&error, __VA_ARGS__); \
		if (error) { \
			fprintf(stderr, "XLAL Error - %s (%s:%d): function %s failed with code %d: %s\n", __func__, __FILE__, __LINE__, #f, (error)->s_errno, (error)->s_message); \
			err = 1; \
		}\
		else err = 0; \
	} while (0)

#define TRY_FRAMEC_FUNCTION(f, ...) do { \
		FrameCError *error = NULL; \
		f(&error, __VA_ARGS__); \
		if (error) { \
			fprintf(stderr, "XLAL Error - %s (%s:%d): function %s failed with code %d: %s\n", __func__, __FILE__, __LINE__, #f, (error)->s_errno, (error)->s_message); \
			return -1; \
		}\
		return 0; \
	} while (0)

#define TRY_FRAMEC_FUNCTION_VAL(retval, errval, f, ...) do { \
		FrameCError *error = NULL; \
		f(&error, __VA_ARGS__); \
		if (error) { \
			fprintf(stderr, "XLAL Error - %s (%s:%d): function %s failed with code %d: %s\n", __func__, __FILE__, __LINE__, #f, (error)->s_errno, (error)->s_message); \
			return errval; \
		}\
		return retval; \
	} while (0)

#define TRY_FRAMEC_FUNCTION_PTR(f, ...) do { \
		FrameCError *error = NULL; \
		void *ptr = f(&error, __VA_ARGS__); \
		if (error) { \
			fprintf(stderr, "XLAL Error - %s (%s:%d): function %s failed with code %d: %s\n", __func__, __FILE__, __LINE__, #f, (error)->s_errno, (error)->s_message); \
			return NULL; \
		}\
		return ptr; \
	} while (0)

#define TRY_FRAMEC_FUNCTION_VOID(f, ...) do { \
		FrameCError *error = NULL; \
		f(&error, __VA_ARGS__); \
		if (error) { \
			fprintf(stderr, "XLAL Error - %s (%s:%d): function %s failed with code %d: %s\n", __func__, __FILE__, __LINE__, #f, (error)->s_errno, (error)->s_message); \
		}\
		return; \
	} while (0)



/* Stream functions */

int XLALFrameUFrFileClose(LALFrameUFrFile * stream)
{
    TRY_FRAMEC_FUNCTION(FrameCFileClose, stream);
}


LALFrameUFrFile *XLALFrameUFrFileOpen(const char *filename, const char *mode)
{
    fr_file_mode_t m;
    if (*mode == 'r')
        m = FRAMEC_FILE_MODE_INPUT;
    else if (*mode == 'w')
        m = FRAMEC_FILE_MODE_OUTPUT;
    else
        return NULL;
    TRY_FRAMEC_FUNCTION_PTR(FrameCFileOpen, filename, m);
}

int XLALFrameUFileCksumValid(LALFrameUFrFile * stream)
{
    /* TODO Enable this when routine becomes available! */
    // TRY_FRAMEC_FUNCTION( FrameCFileCksumValid, stream );
    stream = NULL;
    return 0;
}


/* TOC functions */

void XLALFrameUFrTOCFree(LALFrameUFrTOC * toc)
{
    TRY_FRAMEC_FUNCTION_VOID(FrameCFrTOCFree, toc);
}

LALFrameUFrTOC *XLALFrameUFrTOCRead(LALFrameUFrFile * stream)
{
    TRY_FRAMEC_FUNCTION_PTR(FrameCFrTOCRead, stream);
}

size_t XLALFrameUFrTOCQueryNFrame(const LALFrameUFrTOC * toc)
{
    fr_toc_nframe_t nframe;
    TRY_FRAMEC_FUNCTION_VAL(nframe, -1, FrameCFrTOCQuery, toc, FR_TOC_FIELD_NFRAME, &nframe, FR_TOC_FIELD_LAST);
}

double XLALFrameUFrTOCQueryGTimeModf(double *iptr, const LALFrameUFrTOC * toc, size_t pos)
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

double XLALFrameUFrTOCQueryDt(const LALFrameUFrTOC * toc, size_t pos)
{
    fr_toc_nframe_t nframe;
    fr_toc_dt_element_t *dt;
    double retval;
    int err;
    CALL_FRAMEC_FUNCTION(err, FrameCFrTOCQuery, toc, FR_TOC_FIELD_NFRAME, &nframe, FR_TOC_FIELD_LAST);
    if (err || nframe <= pos)
        return -1.0;
    dt = LALCalloc(nframe, sizeof(*dt));
    if (!dt)
        return -1.0;
    CALL_FRAMEC_FUNCTION(err, FrameCFrTOCQuery, toc, FR_TOC_FIELD_GTIME, dt, nframe, FR_TOC_FIELD_LAST);
    retval = err ? -1.0 : dt[pos];
    LALFree(dt);
    return retval;
}

size_t XLALFrameUFrTOCQueryAdcN(const LALFrameUFrTOC * toc)
{
    fr_toc_adc_n_t n;
    TRY_FRAMEC_FUNCTION_VAL(n, -1, FrameCFrTOCQuery, toc, FR_TOC_FIELD_ADC_N, &n, FR_TOC_FIELD_LAST);
}

const char *XLALFrameUFrTOCQueryAdcName(const LALFrameUFrTOC * toc, size_t adc)
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

size_t XLALFrameUFrTOCQuerySimN(const LALFrameUFrTOC * toc)
{
    fr_toc_sim_n_t n;
    TRY_FRAMEC_FUNCTION_VAL(n, -1, FrameCFrTOCQuery, toc, FR_TOC_FIELD_SIM_N, &n, FR_TOC_FIELD_LAST);
}

const char *XLALFrameUFrTOCQuerySimName(const LALFrameUFrTOC * toc, size_t sim)
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

size_t XLALFrameUFrTOCQueryProcN(const LALFrameUFrTOC * toc)
{
    fr_toc_proc_n_t n;
    TRY_FRAMEC_FUNCTION_VAL(n, -1, FrameCFrTOCQuery, toc, FR_TOC_FIELD_PROC_N, &n, FR_TOC_FIELD_LAST);
}

const char *XLALFrameUFrTOCQueryProcName(const LALFrameUFrTOC * toc, size_t proc)
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

size_t XLALFrameUFrTOCQueryDetectorN(const LALFrameUFrTOC * toc)
{
    fr_toc_detector_n_t n;
    TRY_FRAMEC_FUNCTION_VAL(n, -1, FrameCFrTOCQuery, toc, FR_TOC_FIELD_DETECTOR_N, &n, FR_TOC_FIELD_LAST);
}

const char *XLALFrameUFrTOCQueryDetectorName(const LALFrameUFrTOC * toc, size_t det)
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

/* FrameH functions */

/* frame allocation/free functions */

void XLALFrameUFrameHFree(LALFrameUFrameH * frame)
{
    TRY_FRAMEC_FUNCTION_VOID(FrameCFrameHFree, frame);
}


LALFrameUFrameH *XLALFrameUFrameHAlloc(const char *name, double start, double dt, int frnum)
{
    frame_h_gtime_t gtime;
    gtime.sec = (int) start;
    gtime.nan = (int) floor(0.5 + 1e9 * (start - gtime.sec));
    TRY_FRAMEC_FUNCTION_PTR(FrameCFrameHAlloc, name, gtime, dt, frnum);
}

/* I/O functions for frame */

LALFrameUFrameH *XLALFrameUFrameHRead(LALFrameUFrFile * stream, int pos)
{
    TRY_FRAMEC_FUNCTION_PTR(FrameCFrameHRead, stream, pos);
}

int XLALFrameUFrameHWrite(LALFrameUFrFile * stream, LALFrameUFrameH * frame)
{
    TRY_FRAMEC_FUNCTION(FrameCFrameHWrite, stream, frame);
}

/* function to add a channel to a frame */

int XLALFrameUFrameHFrChanAdd(LALFrameUFrameH * frame, LALFrameUFrChan * channel)
{
    TRY_FRAMEC_FUNCTION(FrameCFrameHFrChanAdd, frame, channel);
}

/* function to add a detector to a frame */

int XLALFrameUFrameHFrDetectorAdd(LALFrameUFrameH * frame, LALFrameUFrDetector * detector)
{
    TRY_FRAMEC_FUNCTION(FrameCFrameHFrDetectorAdd, frame, detector->handle);
}

/* function to add a detector to a frame */

int XLALFrameUFrameHFrHistoryAdd(LALFrameUFrameH * frame, LALFrameUFrHistory * history)
{
    TRY_FRAMEC_FUNCTION(FrameCFrameHFrHistoryAdd, frame, history);
}

/* functions to query frame metadata */

/* note: returns a pointer to internal memory that will be lost when frame is freed */
const char *XLALFrameUFrameHQueryName(const LALFrameUFrameH * frame)
{
    frame_h_name_t name;
    TRY_FRAMEC_FUNCTION_VAL(name, NULL, FrameCFrameHQuery, frame, FRAME_H_FIELD_NAME, &name, FRAME_H_FIELD_LAST);
}

int XLALFrameUFrameHQueryRun(const LALFrameUFrameH * frame)
{
    frame_h_run_t run;
    TRY_FRAMEC_FUNCTION_VAL(run, -1, FrameCFrameHQuery, frame, FRAME_H_FIELD_RUN, &run, FRAME_H_FIELD_LAST);
}

int XLALFrameUFrameHQueryFrame(const LALFrameUFrameH * frame)
{
    frame_h_frame_t frnum;
    TRY_FRAMEC_FUNCTION_VAL(frnum, -1, FrameCFrameHQuery, frame, FRAME_H_FIELD_FRAME, &frnum, FRAME_H_FIELD_LAST);
}

int XLALFrameUFrameHQueryDataQuality(const LALFrameUFrameH * frame)
{
    frame_h_data_quality_t dq;
    TRY_FRAMEC_FUNCTION_VAL(dq, -1, FrameCFrameHQuery, frame, FRAME_H_FIELD_DATA_QUALITY, &dq, FRAME_H_FIELD_LAST);
}

double XLALFrameUFrameHQueryGTimeModf(double *iptr, const LALFrameUFrameH * frame)
{
    frame_h_gtime_t start;
    LIGOTimeGPS epoch;
    TRY_FRAMEC_FUNCTION_VAL(XLALGPSModf(iptr, XLALGPSSet(&epoch, start.sec, start.nan)), XLAL_REAL8_FAIL_NAN, FrameCFrameHQuery, frame, FRAME_H_FIELD_GTIME, &start, FRAME_H_FIELD_LAST);
}

int XLALFrameUFrameHQueryULeapS(const LALFrameUFrameH * frame)
{
    frame_h_uleaps_t uleaps;
    TRY_FRAMEC_FUNCTION_VAL(uleaps, -1, FrameCFrameHQuery, frame, FRAME_H_FIELD_ULEAPS, &uleaps, FRAME_H_FIELD_LAST);
}

double XLALFrameUFrameHQueryDt(const LALFrameUFrameH * frame)
{
    frame_h_dt_t dt;
    TRY_FRAMEC_FUNCTION_VAL(dt, -1, FrameCFrameHQuery, frame, FRAME_H_FIELD_DT, &dt, FRAME_H_FIELD_LAST);
}

/* functions to set frame metadata */

int XLALFrameUFrameHSetRun(LALFrameUFrameH * frame, int run)
{
    frame_h_run_t run_ = run;
    TRY_FRAMEC_FUNCTION(FrameCFrameHSet, frame, FRAME_H_FIELD_RUN, run_, FRAME_H_FIELD_LAST);
}


/* Channel functions */

void XLALFrameUFrChanFree(LALFrameUFrChan * channel)
{
    TRY_FRAMEC_FUNCTION_VOID(FrameCFrChanFree, channel);
}

LALFrameUFrChan *XLALFrameUFrChanRead(LALFrameUFrFile * stream, const char *name, size_t pos)
{
    TRY_FRAMEC_FUNCTION_PTR(FrameCFrChanRead, stream, name, pos);
}

LALFrameUFrChan *XLALFrameUFrAdcChanAlloc(const char *name, int dtype, size_t ndata)
{
    TRY_FRAMEC_FUNCTION_PTR(FrameCFrChanAlloc, name, FR_ADC_CHAN_TYPE, dtype, ndata);
}

LALFrameUFrChan *XLALFrameUFrSimChanAlloc(const char *name, int dtype, size_t ndata)
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

LALFrameUFrChan *XLALFrameUFrProcChanAlloc(const char *name, int type, int subtype, int dtype, size_t ndata)
{
    TRY_FRAMEC_FUNCTION_PTR(FrameCFrProcChanAlloc, name, XLALFrProcType(type), XLALFrProcSubType(subtype), dtype, ndata);
}


/* Channel query functions */

const char *XLALFrameUFrChanQueryName(const LALFrameUFrChan * channel)
{
    fr_chan_name_t name;
    TRY_FRAMEC_FUNCTION_VAL(name, NULL, FrameCFrChanQuery, channel, FR_CHAN_FIELD_NAME, &name, FR_CHAN_FIELD_LAST);
}

double XLALFrameUFrChanQueryTimeOffset(const LALFrameUFrChan * channel)
{
    /*  TODO: implement when functionality is made available */
    /*
       fr_chan_time_offset_t offset;
       TRY_FRAMEC_FUNCTION_VAL( offset, 0, FrameCFrChanQuery, channel, FR_CHAN_FIELD_TIME_OFFSET, &offset, FR_CHAN_FIELD_LAST );
     */
    channel = NULL;
    return 0;
}

/* Channel set functions */

int XLALFrameUFrChanSetSampleRate(LALFrameUFrChan * channel, double sampleRate)
{
    fr_chan_sample_rate_t srate = sampleRate;
    TRY_FRAMEC_FUNCTION(FrameCFrChanSet, channel, FR_CHAN_FIELD_SAMPLE_RATE, srate, FR_CHAN_FIELD_LAST);
}


/* ChanVector functions */

int XLALFrameUFrChanVectorAlloc(LALFrameUFrChan * channel, int dtype, size_t ndata)
{
    TRY_FRAMEC_FUNCTION(FrameCFrChanVectorAlloc, channel, dtype, ndata);
}


/* get the FrameC compression scheme for a particular Frame Spec compressLevel */
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

/* functions to compress and expand FrVect structures within a channel */

int XLALFrameUFrChanVectorCompress(LALFrameUFrChan * channel, int compressLevel)
{
    TRY_FRAMEC_FUNCTION(FrameCFrChanVectorCompress, channel, XLALFrVectCompressionScheme(compressLevel));
}

int XLALFrameUFrChanVectorExpand(LALFrameUFrChan * channel)
{
    TRY_FRAMEC_FUNCTION(FrameCFrChanVectorExpand, channel);
}

/* functions to query FrVect structures within a channel */

/* note: returns a pointer to internal memory which is lost when channel is freed */
const char *XLALFrameUFrChanVectorQueryName(const LALFrameUFrChan * channel)
{
    fr_vect_name_t name;
    TRY_FRAMEC_FUNCTION_VAL(name, NULL, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_DATA, &name, FR_VECT_FIELD_LAST);
}

int XLALFrameUFrChanVectorQueryCompress(const LALFrameUFrChan * channel)
{
    fr_vect_compress_t compress;
    TRY_FRAMEC_FUNCTION_VAL(compress, -1, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_COMPRESS, &compress, FR_VECT_FIELD_LAST);
}

int XLALFrameUFrChanVectorQueryType(const LALFrameUFrChan * channel)
{
    fr_vect_type_t type;
    TRY_FRAMEC_FUNCTION_VAL(type, -1, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_TYPE, &type, FR_VECT_FIELD_LAST);
}

/* retrieves a handle to the data vector in the FrVect structure */
void *XLALFrameUFrChanVectorQueryData(const LALFrameUFrChan * channel)
{
    fr_vect_data_t data;
    TRY_FRAMEC_FUNCTION_VAL(data, NULL, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_DATA, &data, FR_VECT_FIELD_LAST);
}

size_t XLALFrameUFrChanVectorQueryNBytes(const LALFrameUFrChan * channel)
{
    fr_vect_nbytes_t nbytes;
    TRY_FRAMEC_FUNCTION_VAL(nbytes, -1, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_NBYTES, &nbytes, FR_VECT_FIELD_LAST);
}

size_t XLALFrameUFrChanVectorQueryNData(const LALFrameUFrChan * channel)
{
    fr_vect_ndata_t ndata;
    TRY_FRAMEC_FUNCTION_VAL(ndata, -1, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_NDATA, &ndata, FR_VECT_FIELD_LAST);
}

size_t XLALFrameUFrChanVectorQueryNDim(const LALFrameUFrChan * channel)
{
    fr_vect_ndim_t ndim;
    TRY_FRAMEC_FUNCTION_VAL(ndim, -1, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_NDIM, &ndim, FR_VECT_FIELD_LAST);
}

size_t XLALFrameUFrChanVectorQueryNx(const LALFrameUFrChan * channel, size_t dim)
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

double XLALFrameUFrChanVectorQueryDx(const LALFrameUFrChan * channel, size_t dim)
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

double XLALFrameUFrChanVectorQueryStartX(const LALFrameUFrChan * channel, size_t dim)
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

const char *XLALFrameUFrChanVectorQueryUnitX(const LALFrameUFrChan * channel, size_t dim)
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

const char *XLALFrameUFrChanVectorQueryUnitY(const LALFrameUFrChan * channel)
{
    fr_vect_unit_y_t unitY;
    TRY_FRAMEC_FUNCTION_VAL(unitY, NULL, FrameCFrChanVectorQuery, channel, FR_VECT_FIELD_UNIT_Y, &unitY, FR_VECT_FIELD_LAST);
}


/* functions to set FrVect structures within a channel */

int XLALFrameUFrChanVectorSetName(LALFrameUFrChan * channel, const char *name)
{
    TRY_FRAMEC_FUNCTION(FrameCFrChanVectorSet, channel, FR_VECT_FIELD_NAME, (fr_vect_name_t) (name), FR_VECT_FIELD_LAST);
}


/*
int XLALFrameUFrChanVectorSetData(LALFrameUFrChan *channel, void *data, size_t ndata)
{
	fr_vect_ndim_t ndim_ = ndim;
	fr_vect_dx_t *dx_;
	size_t i;
	int err;
	dx_ = LALCalloc(ndim, sizeof(*dx_));
	for (i = 0; i < ndim; ++i)
		dx_[i] = dx[i];
	CALL_FRAMEC_FUNCTION( err, FrameCFrChanVectorSet, channel, FR_VECT_FIELD_DX, dx_, ndim_, FR_VECT_FIELD_LAST );
	LALFree(dx_);
	return err ? -1 : 0;
}
*/


/* TODO:
XLALFrameUFrChanVectorSetDx
XLALFrameUFrChanVectorSetStartX
XLALFrameUFrChanVectorSetUnitX
XLALFrameUFrChanVectorSetUnitY
*/

/*
int XLALFrameUFrChanVectorSetDx(LALFrameUFrChan *channel, const double *dx, size_t ndim)
{
	fr_vect_ndim_t ndim_ = ndim;
	fr_vect_dx_t *dx_;
	size_t i;
	int err;
	dx_ = LALCalloc(ndim, sizeof(*dx_));
	for (i = 0; i < ndim; ++i)
		dx_[i] = dx[i];
	CALL_FRAMEC_FUNCTION( err, FrameCFrChanVectorSet, channel, FR_VECT_FIELD_DX, dx_, ndim_, FR_VECT_FIELD_LAST );
	LALFree(dx_);
	return err ? -1 : 0;
}
*/


/* Detector functions */

void XLALFrameUFrDetectorFree(LALFrameUFrDetector * detector)
{
    if (detector) {
        int err;
        CALL_FRAMEC_FUNCTION(err, FrameCFrDetectorFree, detector->handle);
        LALFree(detector);
    }
    return;
}

LALFrameUFrDetector *XLALFrameUFrDetectorRead(LALFrameUFrFile * stream, const char *name)
{
    fr_detector_prefix_t prefix;
    LALFrameUFrDetector *detector;
    int err;
    detector = LALCalloc(1, sizeof(*detector));
    if (!detector)
        return NULL;
    CALL_FRAMEC_FUNCTION(err, detector->handle = FrameCFrDetectorRead, stream, name);
    if (err || !detector->handle) {
        XLALFrameUFrDetectorFree(detector);
        return NULL;
    }
    /* read the prefix and use it to set the detector's prefix */
    CALL_FRAMEC_FUNCTION(err, FrameCFrDetectorQuery, detector->handle, FR_DETECTOR_FIELD_PREFIX, &prefix, FR_DETECTOR_FIELD_LAST);
    if (err) {
        XLALFrameUFrDetectorFree(detector);
        return NULL;
    }
    memcpy(detector->prefix, prefix.prefix, 2);
    return detector;
}

LALFrameUFrDetector *XLALFrameUFrDetectorAlloc(const char *name, const char *prefix, double latitude, double longitude,
                                               double elevation, double azimuthX, double azimuthY, double altitudeX, double altitudeY, double midpointX, double midpointY, int localTime)
{
    LALFrameUFrDetector *detector;
    int err;
    detector = LALCalloc(1, sizeof(*detector));
    if (!detector)
        return NULL;
    if (prefix)
        memcpy(detector->prefix, prefix, 2);
    CALL_FRAMEC_FUNCTION(err, detector->handle = FrameCFrDetectorAlloc, name, prefix, latitude, longitude, elevation, azimuthX, azimuthY, altitudeX, altitudeY, midpointX, midpointY, localTime);
    if (err || !detector->handle) {
        XLALFrameUFrDetectorFree(detector);
        return NULL;
    }
    return detector;
}

const char *XLALFrameUFrDetectorQueryName(const LALFrameUFrDetector * detector)
{
    fr_detector_name_t name;
    TRY_FRAMEC_FUNCTION_VAL(name, NULL, FrameCFrDetectorQuery, detector->handle, FR_DETECTOR_FIELD_NAME, &name, FR_DETECTOR_FIELD_LAST);
}

const char *XLALFrameUFrDetectorQueryPrefix(const LALFrameUFrDetector * detector)
{
    return detector->prefix;
}

double XLALFrameUFrDetectorQueryLongitude(const LALFrameUFrDetector * detector)
{
    fr_detector_longitude_t longitude;
    TRY_FRAMEC_FUNCTION_VAL(longitude, -1.0, FrameCFrDetectorQuery, detector->handle, FR_DETECTOR_FIELD_LONGITUDE, &longitude, FR_DETECTOR_FIELD_LAST);
}

double XLALFrameUFrDetectorQueryLatitude(const LALFrameUFrDetector * detector)
{
    fr_detector_latitude_t latitude;
    TRY_FRAMEC_FUNCTION_VAL(latitude, -1.0, FrameCFrDetectorQuery, detector->handle, FR_DETECTOR_FIELD_LATITUDE, &latitude, FR_DETECTOR_FIELD_LAST);
}

double XLALFrameUFrDetectorQueryElevation(const LALFrameUFrDetector * detector)
{
    fr_detector_elevation_t elevation;
    TRY_FRAMEC_FUNCTION_VAL(elevation, -1.0, FrameCFrDetectorQuery, detector->handle, FR_DETECTOR_FIELD_ELEVATION, &elevation, FR_DETECTOR_FIELD_LAST);
}

double XLALFrameUFrDetectorQueryArmXAzimuth(const LALFrameUFrDetector * detector)
{
    fr_detector_arm_x_azimuth_t armXazimuth;
    TRY_FRAMEC_FUNCTION_VAL(armXazimuth, -1.0, FrameCFrDetectorQuery, detector->handle, FR_DETECTOR_FIELD_ARM_X_AZIMUTH, &armXazimuth, FR_DETECTOR_FIELD_LAST);
}

double XLALFrameUFrDetectorQueryArmYAzimuth(const LALFrameUFrDetector * detector)
{
    fr_detector_arm_y_azimuth_t armYazimuth;
    TRY_FRAMEC_FUNCTION_VAL(armYazimuth, -1.0, FrameCFrDetectorQuery, detector->handle, FR_DETECTOR_FIELD_ARM_Y_AZIMUTH, &armYazimuth, FR_DETECTOR_FIELD_LAST);
}

double XLALFrameUFrDetectorQueryArmXAltitude(const LALFrameUFrDetector * detector)
{
    fr_detector_arm_x_altitude_t armXaltitude;
    TRY_FRAMEC_FUNCTION_VAL(armXaltitude, -1.0, FrameCFrDetectorQuery, detector->handle, FR_DETECTOR_FIELD_ARM_X_ALTITUDE, &armXaltitude, FR_DETECTOR_FIELD_LAST);
}

double XLALFrameUFrDetectorQueryArmYAltitude(const LALFrameUFrDetector * detector)
{
    fr_detector_arm_y_altitude_t armYaltitude;
    TRY_FRAMEC_FUNCTION_VAL(armYaltitude, -1.0, FrameCFrDetectorQuery, detector->handle, FR_DETECTOR_FIELD_ARM_Y_ALTITUDE, &armYaltitude, FR_DETECTOR_FIELD_LAST);
}

double XLALFrameUFrDetectorQueryArmXMidpoint(const LALFrameUFrDetector * detector)
{
    fr_detector_arm_x_midpoint_t armXmidpoint;
    TRY_FRAMEC_FUNCTION_VAL(armXmidpoint, -1.0, FrameCFrDetectorQuery, detector->handle, FR_DETECTOR_FIELD_ARM_X_MIDPOINT, &armXmidpoint, FR_DETECTOR_FIELD_LAST);
}

double XLALFrameUFrDetectorQueryArmYMidpoint(const LALFrameUFrDetector * detector)
{
    fr_detector_arm_y_midpoint_t armYmidpoint;
    TRY_FRAMEC_FUNCTION_VAL(armYmidpoint, -1.0, FrameCFrDetectorQuery, detector->handle, FR_DETECTOR_FIELD_ARM_Y_MIDPOINT, &armYmidpoint, FR_DETECTOR_FIELD_LAST);
}

int XLALFrameUFrDetectorQueryLocalTime(const LALFrameUFrDetector * detector)
{
    fr_detector_localtime_t localTime;
    TRY_FRAMEC_FUNCTION_VAL(localTime, -1.0, FrameCFrDetectorQuery, detector->handle, FR_DETECTOR_FIELD_LOCAL_TIME, &localTime, FR_DETECTOR_FIELD_LAST);
}


/* History routines */

void XLALFrameUFrHistoryFree(LALFrameUFrHistory * history)
{
    TRY_FRAMEC_FUNCTION_VOID(FrameCFrHistoryFree, history);
}

LALFrameUFrHistory *XLALFrameUFrHistoryAlloc(const char *name, double gpssec, const char *comment)
{
    TRY_FRAMEC_FUNCTION_PTR(FrameCFrHistoryAlloc, name, (fr_history_time_t) floor(gpssec), comment);
}
