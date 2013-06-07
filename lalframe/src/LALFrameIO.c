#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <lal/LALDatatypes.h>
#include <lal/LALDetectors.h>
#include <lal/LALString.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/Date.h>

#include <lal/LALFrameU.h>
#include <lal/LALFrameIO.h>

struct tagLALFrFile {
    LALFrameUFrFile *file;
    LALFrameUFrTOC *toc;
};

int XLALFrFileClose(LALFrFile * frfile)
{
    if (frfile) {
        if (frfile->file) {
            XLALFrameUFrFileClose(frfile->file);
            frfile->file = NULL;
        }
        if (frfile->toc) {
            XLALFrameUFrTOCFree(frfile->toc);
            frfile->toc = NULL;
        }
        LALFree(frfile);
    }
    return 0;
}

LALFrFile *XLALFrFileOpenURL(const char *url)
{
    LALFrFile *frfile;
    char prot[FILENAME_MAX] = "";
    char host[FILENAME_MAX] = "";
    char path[FILENAME_MAX] = "";
    int n;

    XLAL_CHECK_NULL(url, XLAL_EFAULT);
    XLAL_CHECK_NULL(strlen(url) < FILENAME_MAX, XLAL_EBADLEN,
        "url %s is too long", url);

    n = sscanf(url, "%[^:]://%[^/]%s", prot, host, path);
    if (n != 3) {       /* perhaps the hostname has been omitted */
        XLALStringCopy(host, "localhost", sizeof(host));
        if (n != 2) {   /* assume the whole thing is a file path */
            XLALStringCopy(prot, "file", sizeof(prot));
            XLALStringCopy(path, url, sizeof(path));
        }
    }

    /* only file urls are supported at this time */
    if (strcmp(prot, "file"))   /* not a file url */
        XLAL_ERROR_NULL(XLAL_EINVAL, "Unsupported protocol %s", prot);

    /* only files on localhost are supported at this time */
    if (strcmp(host, "localhost")) {    /* not explicitly localhost */
        /* make sure the host *is* the localhost */
        char localhost[FILENAME_MAX];
        gethostname(localhost, FILENAME_MAX - 1);
        if (strcmp(host, localhost))    /* not localhost */
            XLAL_ERROR_NULL(XLAL_EINVAL, "Cannot read files from host %s",
                host);
    }

    /* open frame file in read mode */
    /*
     * frfile = XLALFrameUFrFileOpen(path, "r");
     * if (!frfile)
     * XLAL_ERROR_NULL(XLAL_EIO, "Could not open frame file %s", path);
     */
    frfile = LALMalloc(sizeof(*frfile));
    if (!frfile)
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    frfile->file = XLALFrameUFrFileOpen(path, "r");
    if (!frfile->file) {
        LALFree(frfile);
        XLAL_ERROR_NULL(XLAL_EIO, "Could not open frame file %s", path);
    }
    frfile->toc = XLALFrameUFrTOCRead(frfile->file);
    if (!frfile->toc) {
        XLALFrameUFrFileClose(frfile->file);
        LALFree(frfile);
        XLAL_ERROR_NULL(XLAL_EIO, "Could not open TOC for frame file %s",
            path);
    }

    return frfile;
}

size_t XLALFrFileQueryNFrame(const LALFrFile * frfile)
{
    return XLALFrameUFrTOCQueryNFrame(frfile->toc);
}

LIGOTimeGPS *XLALFrFileQueryGTime(LIGOTimeGPS * start,
    const LALFrFile * frfile, size_t pos)
{
    double ip, fp;      /* integer part and fraction part */
    fp = XLALFrameUFrTOCQueryGTimeModf(&ip, frfile->toc, pos);
    return XLALGPSSet(start, ip, XLAL_BILLION_REAL8 * fp);
}

double XLALFrFileQueryDt(const LALFrFile * frfile, size_t pos)
{
    return XLALFrameUFrTOCQueryDt(frfile->toc, pos);
}

LALTYPECODE XLALFrFileQueryChanType(const LALFrFile * frfile,
    const char *chname, size_t pos)
{
    LALFrameUFrChan *channel;
    int type;
    channel = XLALFrameUFrChanRead(frfile->file, chname, pos);
    if (!channel)
        XLAL_ERROR(XLAL_ENAME);
    type = XLALFrameUFrChanVectorQueryType(channel);
    XLALFrameUFrChanFree(channel);
    switch (type) {
    case LAL_FRAMEU_FR_VECT_C:
        return LAL_CHAR_TYPE_CODE;
    case LAL_FRAMEU_FR_VECT_2S:
        return LAL_I2_TYPE_CODE;
    case LAL_FRAMEU_FR_VECT_8R:
        return LAL_D_TYPE_CODE;
    case LAL_FRAMEU_FR_VECT_4R:
        return LAL_S_TYPE_CODE;
    case LAL_FRAMEU_FR_VECT_4S:
        return LAL_I4_TYPE_CODE;
    case LAL_FRAMEU_FR_VECT_8S:
        return LAL_I8_TYPE_CODE;
    case LAL_FRAMEU_FR_VECT_8C:
        return LAL_C_TYPE_CODE;
    case LAL_FRAMEU_FR_VECT_16C:
        return LAL_Z_TYPE_CODE;
    case LAL_FRAMEU_FR_VECT_2U:
        return LAL_U2_TYPE_CODE;
    case LAL_FRAMEU_FR_VECT_4U:
        return LAL_U4_TYPE_CODE;
    case LAL_FRAMEU_FR_VECT_8U:
        return LAL_U8_TYPE_CODE;
    case LAL_FRAMEU_FR_VECT_1U:
        return LAL_UCHAR_TYPE_CODE;
    case LAL_FRAMEU_FR_VECT_STRING:
        XLAL_ERROR(XLAL_ETYPE,
            "No LAL typecode equivalent to FrVect type STRING");
    default:
        XLAL_ERROR(XLAL_ETYPE, "Unrecognized FrVect type %d", type);
    }
    return -1;  /* never get here anyway... */
}

size_t XLALFrFileQueryChanVectorLength(const LALFrFile * frfile,
    const char *chname, size_t pos)
{
    LALFrameUFrChan *channel;
    size_t length;
    channel = XLALFrameUFrChanRead(frfile->file, chname, pos);
    if (!channel)
        XLAL_ERROR(XLAL_ENAME);
    length = XLALFrameUFrChanVectorQueryNData(channel);
    XLALFrameUFrChanFree(channel);
    return length;
}

int XLALFrFileCksumValid(LALFrFile * frfile)
{
    int result;
    /* this process might mess up the TOC so need to reread it afterwards */
    XLALFrameUFrTOCFree(frfile->toc);
    result = XLALFrameUFileCksumValid(frfile->file);
    frfile->toc = XLALFrameUFrTOCRead(frfile->file);
    return result;
}

#define TDOM 1
#define FDOM 2

#define DOM TDOM

#define TYPE INT2
#define VEXT 2S
#include "LALFrameIO_source.c"
#undef VEXT
#undef TYPE

#define TYPE INT4
#define VEXT 4S
#include "LALFrameIO_source.c"
#undef VEXT
#undef TYPE

#define TYPE INT8
#define VEXT 8S
#include "LALFrameIO_source.c"
#undef VEXT
#undef TYPE

#define TYPE UINT2
#define VEXT 2U
#include "LALFrameIO_source.c"
#undef VEXT
#undef TYPE

#define TYPE UINT4
#define VEXT 4U
#include "LALFrameIO_source.c"
#undef VEXT
#undef TYPE

#define TYPE UINT8
#define VEXT 8U
#include "LALFrameIO_source.c"
#undef VEXT
#undef TYPE

#define TYPE REAL4
#define VEXT 4R
#include "LALFrameIO_source.c"
#undef VEXT
#undef TYPE

#define TYPE REAL8
#define VEXT 8R
#include "LALFrameIO_source.c"
#undef VEXT
#undef TYPE

#define TYPE COMPLEX8
#define VEXT 8C
#include "LALFrameIO_source.c"
#undef VEXT
#undef TYPE

#define TYPE COMPLEX16
#define VEXT 16C
#include "LALFrameIO_source.c"
#undef VEXT
#undef TYPE

#undef DOM
#define DOM FDOM

#define TYPE REAL4
#define VEXT 4R
#include "LALFrameIO_source.c"
#undef VEXT
#undef TYPE

#define TYPE REAL8
#define VEXT 8R
#include "LALFrameIO_source.c"
#undef VEXT
#undef TYPE

#define TYPE COMPLEX8
#define VEXT 8C
#include "LALFrameIO_source.c"
#undef VEXT
#undef TYPE

#define TYPE COMPLEX16
#define VEXT 16C
#include "LALFrameIO_source.c"
#undef VEXT
#undef TYPE

#undef DOM
#undef TDOM
#undef FDOM

int XLALFrameAddFrHistory(LALFrameH * frame, const char *name,
    const char *comment)
{
    LALFrameUFrHistory *history = NULL;
    LIGOTimeGPS now;

    /* get current time */
    if (!XLALGPSTimeNow(&now))
        XLAL_ERROR(XLAL_EFUNC);

    history = XLALFrameUFrHistoryAlloc(name, now.gpsSeconds, comment);
    if (!history)
        XLAL_ERROR(XLAL_EFUNC);

    if (XLALFrameUFrameHFrHistoryAdd(frame, history)) {
        XLALFrameUFrHistoryFree(history);
        XLAL_ERROR(XLAL_EFUNC);
    }

    XLALFrameUFrHistoryFree(history);
    return 0;
}

static time_t XLALSecondsSinceUnixEpoch(struct tm *tm)
{
    return tm->tm_sec + tm->tm_min * 60 + tm->tm_hour * 3600 +
        tm->tm_yday * 86400 + (tm->tm_year - 70) * 31536000 +
        ((tm->tm_year - 69) / 4) * 86400 -
        ((tm->tm_year - 1) / 100) * 86400 +
        ((tm->tm_year + 299) / 400) * 86400;
}

/* determine local seasonal time - UTC (seconds) */
static int XLALLocalTime(char site, int gpssec)
{
    struct tm utc;
    struct tm loc;
    time_t tutc;
    time_t tloc;
    const char *zone;
    char *orig = NULL;

    XLALGPSToUTC(&utc, gpssec);
    tutc = XLALSecondsSinceUnixEpoch(&utc);

    switch (site) {
    case 'G':
        zone = "Europe/Berlin";
        break;
    case 'H':
        zone = "PST8PDT";
        break;
    case 'L':
        zone = "CST6CDT";
        break;
    case 'T':
        zone = "Japan";
        break;
    case 'V':
        zone = "Europe/Rome";
        break;
    default:   /* use UTC */
        return 0;
    }
    orig = getenv("TZ");
    if (orig)
        orig = strdup(orig);
    setenv("TZ", zone, 1);
    tzset();

    loc = *localtime(&tutc);
    tloc = XLALSecondsSinceUnixEpoch(&loc);

    if (orig) {
        setenv("TZ", orig, 1);
        tzset();
        free(orig);
    }

    return (int)floor(0.5 + difftime(tloc, tutc));
}

int XLALFrameAddFrDetector(LALFrameH * frame, const LALFrDetector * detector)
{
    LALFrameUFrDetector *d;
    double gpssec;
    int localTime;
    XLALFrameUFrameHQueryGTimeModf(&gpssec, frame);
    localTime = XLALLocalTime(detector->prefix[0], gpssec);
    d = XLALFrameUFrDetectorAlloc(detector->name, detector->prefix,
        detector->vertexLatitudeRadians,
        detector->vertexLongitudeRadians,
        detector->vertexElevation,
        detector->xArmAzimuthRadians,
        detector->yArmAzimuthRadians,
        detector->xArmAltitudeRadians,
        detector->yArmAltitudeRadians,
        detector->xArmMidpoint, detector->yArmMidpoint, localTime);
    XLALFrameUFrameHFrDetectorAdd(frame, d);
    XLALFrameUFrDetectorFree(d);
    return 0;
}

void XLALFrameFree(LALFrameH * frame)
{
    XLALFrameUFrameHFree(frame);
    return;
}

LALFrameH *XLALFrameNew(const LIGOTimeGPS * epoch, double duration,
    const char *project, int run, int frnum, int detectorFlags)
{
    LALFrameH *frame = NULL;
    int detind;

    /* allocate frame and set run */
    frame =
        XLALFrameUFrameHAlloc(project, XLALGPSGetREAL8(epoch), duration,
        frnum);
    if (!frame)
        XLAL_ERROR_NULL(XLAL_EFUNC);
    if (XLALFrameUFrameHSetRun(frame, run) < 0) {
        XLALFrameFree(frame);
        XLAL_ERROR_NULL(XLAL_EFUNC);
    }

    /* add detectors */
    for (detind = 0; detind < LAL_NUM_DETECTORS; ++detind) {
        int detflg = 1 << 2 * detind;
        if ((detflg & detectorFlags))   /* yes, one ampersand! */
            XLALFrameAddFrDetector(frame,
                &lalCachedDetectors[detind].frDetector);
    }

    return frame;
}

/*
 * FIXME: How Do I Set timeOffset ????
 * TODO: add compressLevel 
 */

#define DEFINE_FR_CHAN_ADD_TS_FUNCTION(chantype, laltype, vectype, compress) \
	int XLALFrameAdd ## laltype ## TimeSeries ## chantype ## Data(LALFrameH *frame, const laltype ## TimeSeries *series) \
	{ \
		/* const char unitX[] = "s"; */ \
		LALFrameUFrChan *channel = NULL; \
		void *data = NULL; \
		channel = XLALFrameUFrAdcChanAlloc(series->name, LAL_FRAMEU_FR_VECT_ ## vectype, series->data->length); \
		if (!channel) \
			goto failure; \
		data = XLALFrameUFrChanVectorQueryData(channel); \
		if (!data) \
			goto failure; \
		memcpy(data, series->data->data, series->data->length * sizeof(*series->data->data)); \
		/* TODO: set a bunch more metadata */ \
		XLALFrameUFrChanSetSampleRate(channel, 1.0/series->deltaT); \
		XLALFrameUFrChanVectorCompress(channel, LAL_FRAMEU_FR_VECT_COMPRESS_ ## compress); \
		XLALFrameUFrameHFrChanAdd(frame, channel); \
		XLALFrameUFrChanFree(channel); \
		return 0; \
	failure: /* unsuccessful exit */ \
		XLALFrameUFrChanFree(channel); \
		XLAL_ERROR(XLAL_EFUNC); \
	}

#define DEFINE_FR_PROC_CHAN_ADD_TS_FUNCTION(laltype, vectype, compress) \
int XLALFrameAdd ## laltype ## TimeSeriesProcData(LALFrameH *frame, const laltype ## TimeSeries *series) \
	{ \
		/* const char unitX[] = "s"; */ \
		LALFrameUFrChan *channel = NULL; \
		void *data = NULL; \
		channel = XLALFrameUFrProcChanAlloc(series->name, LAL_FRAMEU_FR_PROC_TYPE_TIME_SERIES, LAL_FRAMEU_FR_PROC_SUB_TYPE_UNKNOWN, LAL_FRAMEU_FR_VECT_ ## vectype, series->data->length); \
		if (!channel) \
			goto failure; \
		data = XLALFrameUFrChanVectorQueryData(channel); \
		if (!data) \
			goto failure; \
		memcpy(data, series->data->data, series->data->length * sizeof(*series->data->data)); \
		/* TODO: set a bunch more metadata */ \
		XLALFrameUFrChanVectorCompress(channel, LAL_FRAMEU_FR_VECT_COMPRESS_ ## compress); \
		XLALFrameUFrameHFrChanAdd(frame, channel); \
		XLALFrameUFrChanFree(channel); \
		return 0; \
	failure: /* unsuccessful exit */ \
		XLALFrameUFrChanFree(channel); \
		XLAL_ERROR(XLAL_EFUNC); \
	}

#define DEFINE_FR_PROC_CHAN_ADD_FS_FUNCTION(laltype, vectype, compress) \
int XLALFrameAdd ## laltype ## FrequencySeriesProcData(LALFrameH *frame, const laltype ## FrequencySeries *series, int subtype) \
	{ \
		/* const char unitX[] = "s^-1"; */ \
		LALFrameUFrChan *channel = NULL; \
		void *data = NULL; \
		channel = XLALFrameUFrProcChanAlloc(series->name, LAL_FRAMEU_FR_PROC_TYPE_FREQUENCY_SERIES, subtype, LAL_FRAMEU_FR_VECT_ ## vectype, series->data->length); \
		if (!channel) \
			goto failure; \
		data = XLALFrameUFrChanVectorQueryData(channel); \
		if (!data) \
			goto failure; \
		memcpy(data, series->data->data, series->data->length * sizeof(*series->data->data)); \
		/* TODO: set a bunch more metadata */ \
		XLALFrameUFrChanVectorCompress(channel, LAL_FRAMEU_FR_VECT_COMPRESS_ ## compress); \
		XLALFrameUFrameHFrChanAdd(frame, channel); \
		XLALFrameUFrChanFree(channel); \
		return 0; \
	failure: /* unsuccessful exit */ \
		XLALFrameUFrChanFree(channel); \
		XLAL_ERROR(XLAL_EFUNC); \
	}

/* *INDENT-OFF* */
DEFINE_FR_CHAN_ADD_TS_FUNCTION(Adc, INT2,  2S, ZERO_SUPPRESS_WORD_2)
DEFINE_FR_CHAN_ADD_TS_FUNCTION(Adc, INT4,  4S, ZERO_SUPPRESS_WORD_4)
DEFINE_FR_CHAN_ADD_TS_FUNCTION(Adc, REAL4, 4R, DIFF_GZIP)
DEFINE_FR_CHAN_ADD_TS_FUNCTION(Adc, REAL8, 8R, DIFF_GZIP)

DEFINE_FR_CHAN_ADD_TS_FUNCTION(Sim, INT2,  2S, ZERO_SUPPRESS_WORD_2)
DEFINE_FR_CHAN_ADD_TS_FUNCTION(Sim, INT4,  4S, ZERO_SUPPRESS_WORD_4)
DEFINE_FR_CHAN_ADD_TS_FUNCTION(Sim, REAL4, 4R, DIFF_GZIP)
DEFINE_FR_CHAN_ADD_TS_FUNCTION(Sim, REAL8, 8R, DIFF_GZIP)

DEFINE_FR_PROC_CHAN_ADD_TS_FUNCTION(INT2,       2S, ZERO_SUPPRESS_WORD_2)
DEFINE_FR_PROC_CHAN_ADD_TS_FUNCTION(INT4,       4S, ZERO_SUPPRESS_WORD_4)
DEFINE_FR_PROC_CHAN_ADD_TS_FUNCTION(INT8,       8S, DIFF_GZIP)
DEFINE_FR_PROC_CHAN_ADD_TS_FUNCTION(UINT2,      2U, ZERO_SUPPRESS_WORD_2)
DEFINE_FR_PROC_CHAN_ADD_TS_FUNCTION(UINT4,      4U, ZERO_SUPPRESS_WORD_4)
DEFINE_FR_PROC_CHAN_ADD_TS_FUNCTION(UINT8,      8U, DIFF_GZIP)
DEFINE_FR_PROC_CHAN_ADD_TS_FUNCTION(REAL4,      4R, DIFF_GZIP)
DEFINE_FR_PROC_CHAN_ADD_TS_FUNCTION(REAL8,      8R, DIFF_GZIP)
DEFINE_FR_PROC_CHAN_ADD_TS_FUNCTION(COMPLEX8,   8C, DIFF_GZIP)
DEFINE_FR_PROC_CHAN_ADD_TS_FUNCTION(COMPLEX16, 16C, DIFF_GZIP)

DEFINE_FR_PROC_CHAN_ADD_FS_FUNCTION(REAL4,      4R, DIFF_GZIP)
DEFINE_FR_PROC_CHAN_ADD_FS_FUNCTION(REAL8,      8R, DIFF_GZIP)
DEFINE_FR_PROC_CHAN_ADD_FS_FUNCTION(COMPLEX8,   8C, DIFF_GZIP)
DEFINE_FR_PROC_CHAN_ADD_FS_FUNCTION(COMPLEX16, 16C, DIFF_GZIP)
/* *INDENT-ON* */

int XLALFrameWrite(LALFrameH * frame, const char *fname)
{
    // LALFrFile *frfile = NULL;
    LALFrameUFrFile *frfile = NULL;
    char tmpfname[FILENAME_MAX];

    /* open temporary file */
    snprintf(tmpfname, sizeof(tmpfname), "%s.tmp", fname);
    frfile = XLALFrameUFrFileOpen(tmpfname, "w");
    if (!frfile)
        goto failure;

    /* write frame */
    XLALFrameUFrameHWrite(frfile, frame);

    /* close temporary file */
    XLALFrameUFrFileClose(frfile);

    /* rename file */
    rename(tmpfname, fname);
    return 0;

  failure:     /* unsuccessful exit */
    XLALFrameUFrFileClose(frfile);
    /* TODO: remove tempfile */
    return -1;
}

static int charcmp(const void *c1, const void *c2)
{
    char a = *(const char *)c1;
    char b = *(const char *)c2;
    return (a > b) - (a < b);
}

/*
 * Based on a channel name, format a standard frame filename;
 * in the process, determine recognized detectors and sites;
 * returns detector flags, if recognized detectors are found.
 */
static int XLALFrameFileName(char *fname, size_t size, const char *chname,
    const LIGOTimeGPS * epoch, double duration)
{
    char site[LAL_NUM_DETECTORS + 1] = "";
    char *desc;
    const char *cs;
    char *s;
    int detflgs = 0;
    int t0;
    int dt;

    /* parse chname to get identified sites and detectors */
    /* strip out detectors from "XmYn...:"-style prefix */
    for (cs = chname; *cs; cs += 2) {
        int d;
        /* when you get to a colon, you're done! */
        if (*cs == ':')
            break;
        /* see if this is an unexpected format */
        if (strlen(cs) <= 2 || !isupper(cs[0]) || !isdigit(cs[1])) {
            /* parse error so reset detflgs and site */
            detflgs = 0;
            site[0] = 0;
            break;
        }
        /* try to find this detector */
        for (d = 0; d < LAL_NUM_DETECTORS; ++d)
            if (0 == strncmp(cs, lalCachedDetectors[d].frDetector.prefix, 2)) {
                /* found it: put it in sites and detflgs */
                detflgs |= 1 << 2 * d;
                strncat(site, cs, 1);
            }
    }

    /* sort and uniqify sites */
    qsort(site, strlen(site), 1, charcmp);
    cs = s = site;
    for (cs = s = site; *s; ++cs)
        if (*s != *cs)
            *++s = *cs;

    /* description is a modified version of chname */
    /* replace invalid description char with '_' */
    desc = XLALStringDuplicate(chname);
    for (s = desc; *s; ++s)
        if (!isalnum(*s))
            *s = '_';

    /* determine start time field and duration field */
    t0 = epoch->gpsSeconds;
    dt = (int)ceil(XLALGPSGetREAL8(epoch) + duration) - t0;

    /* now format the file name */
    snprintf(fname, size, "%s-%s-%d-%d.gwf", *site ? site : "X", desc, t0,
        dt);

    LALFree(desc);
    return detflgs;
}

#define DEFINE_FR_WRITE_TS_FUNCTION(laltype) \
    int XLALFrWrite ## laltype ## TimeSeries(const laltype ## TimeSeries *series, int frnum) \
    { \
        LALFrameH *frame; \
        double duration; \
        char fname[FILENAME_MAX]; \
        int detflgs; \
        duration = series->deltaT * series->data->length; \
        detflgs = XLALFrameFileName(fname, sizeof(fname), series->name, &series->epoch, duration); \
        frame = XLALFrameNew(&series->epoch, duration, "LAL", 0, frnum, detflgs); \
        XLALFrameAdd ## laltype ## TimeSeriesProcData(frame, series); \
        XLALFrameWrite(frame, fname); \
        XLALFrameFree(frame); \
        return 0; \
    }

/* *INDENT-OFF* */
DEFINE_FR_WRITE_TS_FUNCTION(INT2)
DEFINE_FR_WRITE_TS_FUNCTION(INT4)
DEFINE_FR_WRITE_TS_FUNCTION(INT8)
DEFINE_FR_WRITE_TS_FUNCTION(REAL4)
DEFINE_FR_WRITE_TS_FUNCTION(REAL8)
DEFINE_FR_WRITE_TS_FUNCTION(COMPLEX8)
DEFINE_FR_WRITE_TS_FUNCTION(COMPLEX16)
/* *INDENT-ON* */
#define DEFINE_FR_WRITE_FS_FUNCTION(laltype) \
    int XLALFrWrite ## laltype ## FrequencySeries(const laltype ## FrequencySeries *series, int frnum, int subtype) \
    { \
    	LALFrameH *frame; \
    	double duration; \
    	char fname[FILENAME_MAX]; \
    	int detflgs; \
    	duration = series->deltaF > 0.0 ? 1.0 / series->deltaF : 1.0; \
    	detflgs = XLALFrameFileName(fname, sizeof(fname), series->name, &series->epoch, duration); \
    	frame = XLALFrameNew(&series->epoch, duration, "LAL", 0, frnum, detflgs); \
    	XLALFrameAdd ## laltype ## FrequencySeriesProcData(frame, series, subtype); \
    	XLALFrameWrite(frame, fname); \
    	return 0; \
    }
/* *INDENT-OFF* */
DEFINE_FR_WRITE_FS_FUNCTION(REAL4)
DEFINE_FR_WRITE_FS_FUNCTION(REAL8)
DEFINE_FR_WRITE_FS_FUNCTION(COMPLEX8)
DEFINE_FR_WRITE_FS_FUNCTION(COMPLEX16)
/* *INDENT-ON* */

#if 0
int main(void)
{
    const LIGOTimeGPS epoch = { 1000000000, 0 };
    const double duration = 16.0;
    const char *project = "LIGO";
    const int run = 0;
    const int frnum = 0;
    const int detectorFlags =
        LAL_LHO_4K_DETECTOR_BIT | LAL_LLO_4K_DETECTOR_BIT;
    const double srate = 16384.0;
    double deltaT = 1.0 / srate;
    size_t length = duration * srate;
    size_t j;

    LALFrameH *frame;
    REAL4TimeSeries *series;

    series =
        XLALCreateREAL4TimeSeries("dummy", &epoch, 0.0, deltaT,
        &lalStrainUnit, length);
    for (j = 0; j < length; ++j)
        series->data->data[j] = j % 256;

    frame =
        XLALFrameNew(&epoch, duration, project, run, frnum, detectorFlags);

    XLALFrameAddREAL4TimeSeriesProcData(frame, series);

    XLALDestroyREAL4TimeSeries(series);

    XLALFrameWrite(frame, "dummy.gwf");

    XLALFrameFree(frame);

    return 0;
}
#endif
