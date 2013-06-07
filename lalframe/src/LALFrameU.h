#include <stddef.h>
#include <time.h>

#ifndef _LALFRAMEU_H
#define _LALFRAMEU_H

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif
struct tagLALFrameUFrameH;
struct tagLALFrameUFrFile;
struct tagLALFrameUFrTOC;
struct tagLALFrameUFrChan;
struct tagLALFrameUFrDetector;
struct tagLALFrameUFrHistory;

typedef struct tagLALFrameUFrameH LALFrameUFrameH;
typedef struct tagLALFrameUFrFile LALFrameUFrFile;
typedef struct tagLALFrameUFrTOC LALFrameUFrTOC;
typedef struct tagLALFrameUFrChan LALFrameUFrChan;
typedef struct tagLALFrameUFrDetector LALFrameUFrDetector;
typedef struct tagLALFrameUFrHistory LALFrameUFrHistory;

/* from Appendix B */
enum {
    LAL_FRAMEU_FR_VECT_COMPRESS_RAW = 0,
    LAL_FRAMEU_FR_VECT_COMPRESS_GZIP = 1,
    LAL_FRAMEU_FR_VECT_COMPRESS_DIFF_GZIP = 3,
    LAL_FRAMEU_FR_VECT_COMPRESS_ZERO_SUPPRESS_WORD_2 = 5,
    LAL_FRAMEU_FR_VECT_COMPRESS_ZERO_SUPPRESS_WORD_4 = 8
};

/* from Appendix C */
enum {
    LAL_FRAMEU_FR_VECT_C = 0,
    LAL_FRAMEU_FR_VECT_2S = 1,
    LAL_FRAMEU_FR_VECT_8R = 2,
    LAL_FRAMEU_FR_VECT_4R = 3,
    LAL_FRAMEU_FR_VECT_4S = 4,
    LAL_FRAMEU_FR_VECT_8S = 5,
    LAL_FRAMEU_FR_VECT_8C = 6,
    LAL_FRAMEU_FR_VECT_16C = 7,
    LAL_FRAMEU_FR_VECT_STRING = 8,
    LAL_FRAMEU_FR_VECT_2U = 9,
    LAL_FRAMEU_FR_VECT_4U = 10,
    LAL_FRAMEU_FR_VECT_8U = 11,
    LAL_FRAMEU_FR_VECT_1U = 12,
};

/* From 4.3.2.11 FrProcData */
/* FrProcData types */
enum {
    LAL_FRAMEU_FR_PROC_TYPE_UNKNOWN = 0,
    LAL_FRAMEU_FR_PROC_TYPE_TIME_SERIES = 1,
    LAL_FRAMEU_FR_PROC_TYPE_FREQUENCY_SERIES = 2,
    LAL_FRAMEU_FR_PROC_TYPE_OTHER_1D_SERIES = 3,
    LAL_FRAMEU_FR_PROC_TYPE_TIME_FREQUENCY = 4,
    LAL_FRAMEU_FR_PROC_TYPE_WAVELET = 5,
    LAL_FRAMEU_Fr_PROC_TYPE_MULTI_DIMENSIONAL = 6
};

/* FrProcData sub-types for frequency series */
enum {
    LAL_FRAMEU_FR_PROC_SUB_TYPE_UNKNOWN = 0,
    LAL_FRAMEU_FR_PROC_SUB_TYPE_DFT = 1,
    LAL_FRAMEU_FR_PROC_SUB_TYPE_AMPLITUDE_SPECTRAL_DENSITY = 2,
    LAL_FRAMEU_FR_PROC_SUB_TYPE_POWER_SPECTRAL_DENSITY = 3,
    LAL_FRAMEU_FR_PROC_SUB_TYPE_CROSS_SPECTRAL_DENSITY = 4,
    LAL_FRAMEU_FR_PROC_SUB_TYPE_COHERENCE = 5,
    LAL_FRAMEU_FR_PROC_SUB_TYPE_TRANSFER_FUNCTION = 6
};

/* TODO: add routines:
int XLALFrameUFrFileIGWDVersion(LALFrameUFrFile *stream);
*/

void XLALFrameUFrFileClose(LALFrameUFrFile * stream);
LALFrameUFrFile *XLALFrameUFrFileOpen(const char *filename, const char *mode);
int XLALFrameUFileCksumValid(LALFrameUFrFile * stream);

void XLALFrameUFrTOCFree(LALFrameUFrTOC * toc);
LALFrameUFrTOC *XLALFrameUFrTOCRead(LALFrameUFrFile * stream);

size_t XLALFrameUFrTOCQueryNFrame(const LALFrameUFrTOC * toc);
double XLALFrameUFrTOCQueryGTimeModf(double *iptr, const LALFrameUFrTOC * toc,
    size_t pos);
double XLALFrameUFrTOCQueryDt(const LALFrameUFrTOC * toc, size_t pos);
size_t XLALFrameUFrTOCQueryAdcN(const LALFrameUFrTOC * toc);
const char *XLALFrameUFrTOCQueryAdcName(const LALFrameUFrTOC * toc,
    size_t adc);
size_t XLALFrameUFrTOCQuerySimN(const LALFrameUFrTOC * toc);
const char *XLALFrameUFrTOCQuerySimName(const LALFrameUFrTOC * toc,
    size_t sim);
size_t XLALFrameUFrTOCQueryProcN(const LALFrameUFrTOC * toc);
const char *XLALFrameUFrTOCQueryProcName(const LALFrameUFrTOC * toc,
    size_t proc);
size_t XLALFrameUFrTOCQueryDetectorN(const LALFrameUFrTOC * toc);
const char *XLALFrameUFrTOCQueryDetectorName(const LALFrameUFrTOC * toc,
    size_t det);

void XLALFrameUFrameHFree(LALFrameUFrameH * frame);
LALFrameUFrameH *XLALFrameUFrameHAlloc(const char *name, double start,
    double dt, int frnum);
LALFrameUFrameH *XLALFrameUFrameHRead(LALFrameUFrFile * stream, int pos);
int XLALFrameUFrameHWrite(LALFrameUFrFile * stream, LALFrameUFrameH * frame);
int XLALFrameUFrameHFrChanAdd(LALFrameUFrameH * frame,
    LALFrameUFrChan * channel);
int XLALFrameUFrameHFrDetectorAdd(LALFrameUFrameH * frame,
    LALFrameUFrDetector * detector);
int XLALFrameUFrameHFrHistoryAdd(LALFrameUFrameH * frame,
    LALFrameUFrHistory * history);

const char *XLALFrameUFrameHQueryName(const LALFrameUFrameH * frame);
int XLALFrameUFrameHQueryRun(const LALFrameUFrameH * frame);
int XLALFrameUFrameHQueryFrame(const LALFrameUFrameH * frame);
int XLALFrameUFrameHQueryDataQuality(const LALFrameUFrameH * frame);
double XLALFrameUFrameHQueryGTimeModf(double *iptr,
    const LALFrameUFrameH * frame);
int XLALFrameUFrameHQueryULeapS(const LALFrameUFrameH * frame);
double XLALFrameUFrameHQueryDt(const LALFrameUFrameH * frame);

int XLALFrameUFrameHSetRun(LALFrameUFrameH * frame, int run);

void XLALFrameUFrChanFree(LALFrameUFrChan * channel);
LALFrameUFrChan *XLALFrameUFrChanRead(LALFrameUFrFile * stream,
    const char *name, size_t pos);
LALFrameUFrChan *XLALFrameUFrAdcChanAlloc(const char *name, int dtype,
    size_t ndata);
LALFrameUFrChan *XLALFrameUFrSimChanAlloc(const char *name, int dtype,
    size_t ndata);
LALFrameUFrChan *XLALFrameUFrProcChanAlloc(const char *name, int type,
    int subtype, int dtype, size_t ndata);

const char *XLALFrameUFrChanQueryName(const LALFrameUFrChan * channel);
double XLALFrameUFrChanQueryTimeOffset(const LALFrameUFrChan * channel);

int XLALFrameUFrChanSetSampleRate(LALFrameUFrChan * channel,
    double sampleRate);

int XLALFrameUFrChanVectorAlloc(LALFrameUFrChan * channel, int dtype,
    size_t ndata);
int XLALFrameUFrChanVectorCompress(LALFrameUFrChan * channel,
    int compressLevel);
int XLALFrameUFrChanVectorExpand(LALFrameUFrChan * channel);

const char *XLALFrameUFrChanVectorQueryName(const LALFrameUFrChan * channel);
int XLALFrameUFrChanVectorQueryCompress(const LALFrameUFrChan * channel);
int XLALFrameUFrChanVectorQueryType(const LALFrameUFrChan * channel);
void *XLALFrameUFrChanVectorQueryData(const LALFrameUFrChan * channel);
size_t XLALFrameUFrChanVectorQueryNBytes(const LALFrameUFrChan * channel);
size_t XLALFrameUFrChanVectorQueryNData(const LALFrameUFrChan * channel);
size_t XLALFrameUFrChanVectorQueryNDim(const LALFrameUFrChan * channel);
size_t XLALFrameUFrChanVectorQueryNx(const LALFrameUFrChan * channel,
    size_t dim);
double XLALFrameUFrChanVectorQueryDx(const LALFrameUFrChan * channel,
    size_t dim);
double XLALFrameUFrChanVectorQueryStartX(const LALFrameUFrChan * channel,
    size_t dim);
const char *XLALFrameUFrChanVectorQueryUnitX(const LALFrameUFrChan * channel,
    size_t dim);
const char *XLALFrameUFrChanVectorQueryUnitY(const LALFrameUFrChan * channel);

int XLALFrameUFrChanVectorSetName(LALFrameUFrChan * channel,
    const char *name);
int XLALFrameUFrChanVectorSetDx(LALFrameUFrChan * channel, double dx);
int XLALFrameUFrChanVectorSetStartX(LALFrameUFrChan * channel, double x0);
int XLALFrameUFrChanVectorSetUnitX(LALFrameUFrChan * channel,
    const char *unit);
int XLALFrameUFrChanVectorSetUnitY(LALFrameUFrChan * channel,
    const char *unit);

/* TODO: a bunch more things to set coming up!!! */

void XLALFrameUFrDetectorFree(LALFrameUFrDetector * detector);
LALFrameUFrDetector *XLALFrameUFrDetectorRead(LALFrameUFrFile * stream,
    const char *name);
LALFrameUFrDetector *XLALFrameUFrDetectorAlloc(const char *name,
    const char *prefix, double latitude, double longitude, double elevation,
    double azimuthX, double azimuthY, double altitudeX, double altitudeY,
    double midpointX, double midpointY, int localTime);

const char *XLALFrameUFrDetectorQueryName(const LALFrameUFrDetector *
    detector);
const char *XLALFrameUFrDetectorQueryPrefix(const LALFrameUFrDetector *
    detector);
double XLALFrameUFrDetectorQueryLongitude(const LALFrameUFrDetector *
    detector);
double XLALFrameUFrDetectorQueryLatitude(const LALFrameUFrDetector *
    detector);
double XLALFrameUFrDetectorQueryElevation(const LALFrameUFrDetector *
    detector);
double XLALFrameUFrDetectorQueryArmXAzimuth(const LALFrameUFrDetector *
    detector);
double XLALFrameUFrDetectorQueryArmYAzimuth(const LALFrameUFrDetector *
    detector);
double XLALFrameUFrDetectorQueryArmXAltitude(const LALFrameUFrDetector *
    detector);
double XLALFrameUFrDetectorQueryArmYAltitude(const LALFrameUFrDetector *
    detector);
double XLALFrameUFrDetectorQueryArmXMidpoint(const LALFrameUFrDetector *
    detector);
double XLALFrameUFrDetectorQueryArmYMidpoint(const LALFrameUFrDetector *
    detector);
int XLALFrameUFrDetectorQueryLocalTime(const LALFrameUFrDetector * detector);

void XLALFrameUFrHistoryFree(LALFrameUFrHistory * history);
LALFrameUFrHistory *XLALFrameUFrHistoryAlloc(const char *name, double gpssec,
    const char *comment);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif /* _LALFRAMEU_H */
