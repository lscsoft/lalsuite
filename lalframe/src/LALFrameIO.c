/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton, Robert Adam Mercer, Xavier Siemens
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
#include <unistd.h>
#ifndef HAVE_GETHOSTNAME_PROTOTYPE
int gethostname(char *name, int len);
#endif

#include <math.h>

#include <lal/LALDetectors.h>
#include <lal/LALStdio.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>

#include <lal/LALString.h>
#include <lal/LALCache.h>
#include <lal/LALFrameIO.h>

/* FIXME: WARNING: this value might need to change in the future */
#define FR_FILE_HEADER_SIZE 40  /* size of frame file header in bytes */

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

LALFrFile *XLALFrFileOpenURL(const char *url)
{
    struct FrFile *frfile = NULL;
    char prot[FILENAME_MAX] = "";
    char host[FILENAME_MAX] = "";
    char path[FILENAME_MAX] = "";
    int n;

    if (!url)
        XLAL_ERROR_NULL(XLAL_EFAULT);
    if (strlen(url) >= FILENAME_MAX) {
        XLALPrintError("XLAL Error - %s: URL too long: %s\n", __func__,
                       url);
        XLAL_ERROR_NULL(XLAL_EBADLEN);
    }

    n = sscanf(url, "%[^:]://%[^/]%s", prot, host, path);
    if (n != 3) {       /* perhaps the hostname has been omitted */
        strncpy(host, "localhost", sizeof(host) - 1);
        n = sscanf(url, "%[^:]://%s", prot, path);
        if (n != 2) {   /* assume the whole thing is a file path */
            strncpy(prot, "file", sizeof(prot) - 1);
            strncpy(path, url, sizeof(path) - 1);
        }
    }

    if (strcmp(prot, "file")) { /* not a file URL */
        XLALPrintError("XLAL Error - %s: unsupported protocol %s\n",
                       __func__, prot);
        XLAL_ERROR_NULL(XLAL_EINVAL);
    } else {    /* OK: this is a file URL */

        if (strcmp(host, "localhost")) {        /* not explicitly localhost *//* make sure that the host is the localhost */
            char localhost[FILENAME_MAX];
            gethostname(localhost, FILENAME_MAX - 1);
            if (strcmp(host, localhost)) {      /* not localhost */
                XLALPrintError
                    ("XLAL Error - %s: cannot read files from remote host %s\n",
                     __func__, host);
                XLAL_ERROR_NULL(XLAL_EINVAL);
            }
        }
        frfile = FrFileINew(path);
        if (!frfile) {
            XLALPrintError
                ("XLAL Error - %s: could not open frame file %s\n",
                 __func__, path);
            XLAL_ERROR_NULL(XLAL_EIO);
        }
    }

    return frfile;
}


/* code taken from the FrCheck program in the FrameL library */
int XLALFrFileCksumValid(const LALFrFile * frfile)
{
    union {
        const LALFrFile *cf;
        LALFrFile *f;
    } u = {
    frfile};
    LALFrFile *iFile = u.f;
    LALFrameH *frame = NULL;
    int retval = 0;
    FRBOOL chkSumFiFlag = iFile->chkSumFiFlag;
    FRBOOL chkSumFrFlag = iFile->chkSumFrFlag;
    iFile->chkSumFiFlag = FR_YES;
    iFile->chkSumFrFlag = FR_NO;
    /* sequential read */
    /* FIXME: WARNING: HACK: following line may need to be updated in future */
    FrIOSet(iFile->frfd, FR_FILE_HEADER_SIZE);
    while ((frame = FrameReadRecycle(iFile, frame)));
    iFile->chkSumFiFlag = chkSumFiFlag;
    iFile->chkSumFrFlag = chkSumFrFlag;
    if (iFile->error != FR_OK) {
        XLALPrintError("XLAL Error - %s: %s\n", __func__,
                       FrErrorGetHistory());
        XLAL_ERROR(XLAL_EFAILED);
        return -1;
    }
    if (iFile->chkTypeFiRead == 0) {
        XLALPrintWarning("XLAL Warning - %s: missing checksum\n",
                         __func__);
        retval = 1;     /* missing checksum */
    } else if (iFile->chkSumFiRead != iFile->chkSumFi) {
        XLALPrintError("XLAL Error - %s: bad checksum\n", __func__);
        retval = -1;    /* bad checksum */
    } else
        retval = 0;
    FrFileIRewind(iFile);
    return retval;
}


int XLALFrameAddFrHistory(LALFrameH * frame, const char *name,
                          const char *comment)
{
    union {
        const char *cs;
        char *s;
    } namecnvrt;        /* get rid of const qual */
    union {
        const char *cs;
        char *s;
    } commentcnvrt;     /* get rid of const qual */
    LIGOTimeGPS now;
    FrHistory *history;

    /* get current time */
    if (!XLALGPSTimeNow(&now))
        XLAL_ERROR(XLAL_EFUNC);

    /* this nonsense is to convert const char * to char * ... don't worry,
     * the frame library just copies the string anyway */
    namecnvrt.cs = name;
    commentcnvrt.cs = comment;

    /* now create the history */
    history = FrHistoryNew(namecnvrt.s, now.gpsSeconds, commentcnvrt.s);
    if (!history)
        XLAL_ERROR(XLAL_EERR);  /* "internal" error */

    /* attach history to the frame structure */
    if (frame) {
        /* behaviour is identical to FrHistoryAdd if name is NULL */
        if (!name)
            FrStrCpy(&history->name, frame->name);
        history->next = frame->history;
        frame->history = history;
    }

    return 0;
}

static FrDetector *XLALFrDetectorNew(int detector)
{
    const LALDetector *lalDetector;
    FrDetector *frDetector;
    char *detectorName;

    if (detector < 0 || detector >= LAL_NUM_DETECTORS)
        XLAL_ERROR_NULL(XLAL_EINVAL);
    lalDetector = lalCachedDetectors + detector;

    detectorName = XLALStringDuplicate(lalDetector->frDetector.name);
    frDetector = FrDetectorNew(detectorName);
    LALFree(detectorName);
    if (!frDetector)
        XLAL_ERROR_NULL(XLAL_EERR);     /* "internal" error */

    memcpy(frDetector->prefix, lalDetector->frDetector.prefix, 2);
    frDetector->longitude = lalDetector->frDetector.vertexLongitudeRadians;
    frDetector->latitude = lalDetector->frDetector.vertexLatitudeRadians;
    frDetector->elevation = lalDetector->frDetector.vertexElevation;
    frDetector->armXazimuth = lalDetector->frDetector.xArmAzimuthRadians;
    frDetector->armYazimuth = lalDetector->frDetector.yArmAzimuthRadians;
    frDetector->armXaltitude = lalDetector->frDetector.xArmAltitudeRadians;
    frDetector->armYaltitude = lalDetector->frDetector.yArmAltitudeRadians;
    frDetector->armXmidpoint = lalDetector->frDetector.xArmMidpoint;
    frDetector->armYmidpoint = lalDetector->frDetector.yArmMidpoint;
    frDetector->localTime = 0;

    return frDetector;
}

int XLALFrameAddFrDetector(LALFrameH * frame, const LALDetector * detector)
{
    FrDetector *d;
    char *detectorName;
    detectorName = XLALStringDuplicate(detector->frDetector.name);
    d = FrDetectorNew(detectorName);
    LALFree(detectorName);
    if (!d)
        XLAL_ERROR(XLAL_EERR);  /* "internal" error */

    memcpy(d->prefix, detector->frDetector.prefix, 2);
    d->longitude = detector->frDetector.vertexLongitudeRadians;
    d->latitude = detector->frDetector.vertexLatitudeRadians;
    d->elevation = detector->frDetector.vertexElevation;
    d->armXazimuth = detector->frDetector.xArmAzimuthRadians;
    d->armYazimuth = detector->frDetector.yArmAzimuthRadians;
    d->armXaltitude = detector->frDetector.xArmAltitudeRadians;
    d->armYaltitude = detector->frDetector.yArmAltitudeRadians;
    d->armXmidpoint = detector->frDetector.xArmMidpoint;
    d->armYmidpoint = detector->frDetector.yArmMidpoint;
    d->localTime = 0;

    /* add this detector */
    d->next = frame->detectProc;
    frame->detectProc = d;
    return 0;
}


void XLALFrameFree(LALFrameH * frame)
{
    FrameFree(frame);
}


LALFrameH *XLALFrameNew(LIGOTimeGPS * epoch, double duration,
                        const char *project, int run, int frnum,
                        int detectorFlags)
{
    static char histidname[] = __FILE__ " Id";
    static char histtagname[] = __FILE__ " Tag";
  /** \deprecated FIXME: the following code uses obsolete CVS ID tags.
   *  It should be modified to use git version information. */
    static char rcsname[] = "$Name$";
    static char rcsid[] = "$Id$";
    int detector;
    char *proj;
    LALFrameH *frame;

    proj = XLALStringDuplicate(project);
    if (!proj)
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    frame = FrameHNew(proj);
    LALFree(proj);
    if (!frame)
        XLAL_ERROR_NULL(XLAL_EERR);     /* "internal" error */

    frame->run = run;
    frame->frame = frnum;
    frame->GTimeS = epoch->gpsSeconds;
    frame->GTimeN = epoch->gpsNanoSeconds;
    frame->ULeapS = XLALLeapSeconds(epoch->gpsSeconds);
    frame->dt = duration;
    /* TODO: what about dataQuality ? */

    for (detector = 0; detector < LAL_NUM_DETECTORS; ++detector) {
        int detector_flag = 1 << 2 * detector;
        if ((detector_flag & detectorFlags)) {  /* yes, one ampersand! */
            /* add this detector */
            FrDetector *d;
            d = XLALFrDetectorNew(detector);
            if (!d)
                XLAL_ERROR_NULL(XLAL_EFUNC);
            d->next = frame->detectProc;
            frame->detectProc = d;
        }
    }

    /* add history: name of history field is this function's name */
    if (XLALFrameAddFrHistory(frame, histidname, rcsid) < 0)
        XLAL_ERROR_NULL(XLAL_EFUNC);
    if (XLALFrameAddFrHistory(frame, histtagname, rcsname) < 0)
        XLAL_ERROR_NULL(XLAL_EFUNC);
    return frame;
}


static FrVect *XLALFrVectINT4TimeSeries(INT4TimeSeries * series)
{
    char seconds[LALUnitTextSize] = "s";
    char units[LALUnitTextSize];
    FrVect *vect;

    if (NULL ==
        XLALUnitAsString(units, sizeof(units), &series->sampleUnits))
        XLAL_ERROR_NULL(XLAL_EFUNC);

    vect =
        FrVectNew1D(series->name, FR_VECT_4S, series->data->length,
                    series->deltaT, seconds, units);
    if (!vect)
        XLAL_ERROR_NULL(XLAL_EERR);     /* "internal" error */
    vect->startX[0] = 0.0;

    memcpy(vect->data, series->data->data,
           series->data->length * sizeof(*series->data->data));

    FrVectCompress(vect, 8, 0);
    if (vect->compress == 0)
        FrVectCompress(vect, 6, 1);

    return vect;
}


static FrVect *XLALFrVectREAL4TimeSeries(REAL4TimeSeries * series)
{
    char seconds[LALUnitTextSize] = "s";
    char units[LALUnitTextSize];
    FrVect *vect;

    if (NULL ==
        XLALUnitAsString(units, sizeof(units), &series->sampleUnits))
        XLAL_ERROR_NULL(XLAL_EFUNC);

    vect =
        FrVectNew1D(series->name, FR_VECT_4R, series->data->length,
                    series->deltaT, seconds, units);
    if (!vect)
        XLAL_ERROR_NULL(XLAL_EERR);     /* "internal" error */
    vect->startX[0] = 0.0;

    memcpy(vect->data, series->data->data,
           series->data->length * sizeof(*series->data->data));

    FrVectCompress(vect, 8, 0);
    if (vect->compress == 0)
        FrVectCompress(vect, 6, 1);

    return vect;
}


static FrVect *XLALFrVectREAL8TimeSeries(REAL8TimeSeries * series)
{
    char seconds[LALUnitTextSize] = "s";
    char units[LALUnitTextSize];
    FrVect *vect;

    if (NULL ==
        XLALUnitAsString(units, sizeof(units), &series->sampleUnits))
        XLAL_ERROR_NULL(XLAL_EFUNC);

    vect =
        FrVectNew1D(series->name, FR_VECT_8R, series->data->length,
                    series->deltaT, seconds, units);
    if (!vect)
        XLAL_ERROR_NULL(XLAL_EERR);     /* "internal" error */
    vect->startX[0] = 0.0;

    memcpy(vect->data, series->data->data,
           series->data->length * sizeof(*series->data->data));

    FrVectCompress(vect, 6, 1);

    return vect;
}


int XLALFrameAddREAL8TimeSeriesProcData(LALFrameH * frame,
                                        REAL8TimeSeries * series)
{
    LIGOTimeGPS frameEpoch;
    FrProcData *proc;
    FrVect *vect;
    REAL8 duration;

    duration = series->deltaT * series->data->length;

    vect = XLALFrVectREAL8TimeSeries(series);
    if (!vect)
        XLAL_ERROR(XLAL_EFUNC);

    proc = FrProcDataNewV(frame, vect);
    if (!proc) {
        FrVectFree(vect);
        XLAL_ERROR(XLAL_EERR);
    }

    /* time offset: compute this from frame time */
    frameEpoch.gpsSeconds = frame->GTimeS;
    frameEpoch.gpsNanoSeconds = frame->GTimeN;
    proc->timeOffset = XLALGPSDiff(&series->epoch, &frameEpoch);

    /* remaining metadata */
    proc->type = 1;
    proc->subType = 0;
    proc->tRange = duration;
    proc->fShift = 0.0;
    proc->phase = 0.0;
    proc->BW = 0.0;

    return 0;
}

int XLALFrameAddREAL4TimeSeriesProcData(LALFrameH * frame,
                                        REAL4TimeSeries * series)
{
    LIGOTimeGPS frameEpoch;
    FrProcData *proc;
    FrVect *vect;
    REAL8 duration;

    duration = series->deltaT * series->data->length;

    vect = XLALFrVectREAL4TimeSeries(series);
    if (!vect)
        XLAL_ERROR(XLAL_EFUNC);

    proc = FrProcDataNewV(frame, vect);
    if (!proc) {
        FrVectFree(vect);
        XLAL_ERROR(XLAL_EERR);
    }

    /* time offset: compute this from frame time */
    frameEpoch.gpsSeconds = frame->GTimeS;
    frameEpoch.gpsNanoSeconds = frame->GTimeN;
    proc->timeOffset = XLALGPSDiff(&series->epoch, &frameEpoch);

    /* remaining metadata */
    proc->type = 1;
    proc->subType = 0;
    proc->tRange = duration;
    proc->fShift = 0.0;
    proc->phase = 0.0;
    proc->BW = 0.0;

    return 0;
}

int XLALFrameAddINT4TimeSeriesProcData(LALFrameH * frame,
                                       INT4TimeSeries * series)
{
    LIGOTimeGPS frameEpoch;
    FrProcData *proc;
    FrVect *vect;
    REAL8 duration;

    duration = series->deltaT * series->data->length;

    vect = XLALFrVectINT4TimeSeries(series);
    if (!vect)
        XLAL_ERROR(XLAL_EFUNC);

    proc = FrProcDataNewV(frame, vect);
    if (!proc) {
        FrVectFree(vect);
        XLAL_ERROR(XLAL_EERR);
    }

    /* time offset: compute this from frame time */
    frameEpoch.gpsSeconds = frame->GTimeS;
    frameEpoch.gpsNanoSeconds = frame->GTimeN;
    proc->timeOffset = XLALGPSDiff(&series->epoch, &frameEpoch);

    /* remaining metadata */
    proc->type = 1;
    proc->subType = 0;
    proc->tRange = duration;
    proc->fShift = 0.0;
    proc->phase = 0.0;
    proc->BW = 0.0;

    return 0;
}


int XLALFrameAddREAL4TimeSeriesSimData(LALFrameH * frame,
                                       REAL4TimeSeries * series)
{
    LIGOTimeGPS frameEpoch;
    FrSimData *sim;
    FrVect *vect;

    vect = XLALFrVectREAL4TimeSeries(series);
    if (!vect)
        XLAL_ERROR(XLAL_EFUNC);

    sim =
        FrSimDataNew(frame, series->name, 1. / series->deltaT,
                     series->data->length, -32);
    if (!sim) {
        FrVectFree(vect);
        XLAL_ERROR(XLAL_EERR);
    }
    FrVectFree(sim->data);
    sim->data = vect;

    /* time offset: compute this from frame time */
    frameEpoch.gpsSeconds = frame->GTimeS;
    frameEpoch.gpsNanoSeconds = frame->GTimeN;
    sim->timeOffset = XLALGPSDiff(&series->epoch, &frameEpoch);

    /* remaining metadata */
    sim->fShift = 0;
    sim->phase = 0;

    return 0;
}

int XLALFrameAddREAL8TimeSeriesSimData(LALFrameH * frame,
                                       REAL8TimeSeries * series)
{
    LIGOTimeGPS frameEpoch;
    FrSimData *sim;
    FrVect *vect;

    vect = XLALFrVectREAL8TimeSeries(series);
    if (!vect)
        XLAL_ERROR(XLAL_EFUNC);

    sim =
        FrSimDataNew(frame, series->name, 1. / series->deltaT,
                     series->data->length, -64);
    if (!sim) {
        FrVectFree(vect);
        XLAL_ERROR(XLAL_EERR);
    }
    FrVectFree(sim->data);
    sim->data = vect;

    /* time offset: compute this from frame time */
    frameEpoch.gpsSeconds = frame->GTimeS;
    frameEpoch.gpsNanoSeconds = frame->GTimeN;
    sim->timeOffset = XLALGPSDiff(&series->epoch, &frameEpoch);

    /* remaining metadata */
    sim->fShift = 0;
    sim->phase = 0;

    return 0;
}

int XLALFrameAddREAL4TimeSeriesAdcData(LALFrameH * frame,
                                       REAL4TimeSeries * series)
{
    FrAdcData *adc;
    int i;

    adc =
        FrAdcDataNew(frame, series->name, 1. / series->deltaT,
                     series->data->length, -32);
    if (!adc) {
        XLAL_ERROR(XLAL_EERR);
    }

    for (i = 0; i < (int) series->data->length; i++) {
        adc->data->dataF[i] = series->data->data[i];
        /*    fprintf(stdout,"%d %f %f\n",i, adc->data->dataF[i], series->data->data[i]); */
    }

    FrVectCompress(adc->data, 8, 0);
    if (adc->data->compress == 0)
        FrVectCompress(adc->data, 6, 1);

    return 0;
}

int XLALFrameWrite(LALFrameH * frame, const char *fname, int compressLevel)
{
    LALFrFile *frfile;
    char tmpfname[FILENAME_MAX];
    int c;

    /* set temporary filename */
    c = snprintf(tmpfname, sizeof(tmpfname), "%s.tmp", fname);
    if (c < 0 || c > (int) sizeof(tmpfname) - 2)
        XLAL_ERROR(XLAL_ENAME);

    /* open temporary file */
    frfile = FrFileONew(tmpfname, compressLevel);
    if (!frfile)
        XLAL_ERROR(XLAL_EIO);

    /* write frame to temporary filename */
    FrameWrite(frame, frfile);
    FrFileOEnd(frfile);

    /* rename tmpfile */
    if (rename(tmpfname, fname) < 0)
        XLAL_ERROR(XLAL_ESYS);

    return 0;
}
