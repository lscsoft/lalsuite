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

#include <stdlib.h>
#include <string.h>
#include <lal/LALFrameU.h>
#include "utils.h"

LALFrameUFrameH *framecpy(LALFrameUFrFile * frfile, size_t pos)
{
    LALFrameUFrameH *iframe;
    LALFrameUFrameH *oframe;
    iframe = XLALFrameUFrameHRead(frfile, pos);
    oframe = framedup(iframe);
    XLALFrameUFrameHFree(iframe);
    return oframe;
}

LALFrameUFrameH *framedup(LALFrameUFrameH * frame)
{
    LALFrameUFrameH *output;
    const char *name;
    double tfp, tip, dt;
    int frnum, run;

    if (!frame)
        return NULL;

    name = XLALFrameUFrameHQueryName(frame);
    tfp = XLALFrameUFrameHQueryGTimeModf(&tip, frame);
    dt = XLALFrameUFrameHQueryDt(frame);
    frnum = XLALFrameUFrameHQueryFrame(frame);
    run = XLALFrameUFrameHQueryRun(frame);

    output = XLALFrameUFrameHAlloc(name, tip + tfp, dt, frnum);
    XLALFrameUFrameHSetRun(output, run);

    return output;
}

int copydetectors(LALFrameUFrameH * frame, LALFrameUFrFile * frfile)
{
    LALFrameUFrTOC *toc;
    size_t ndet, det;
    toc = XLALFrameUFrTOCRead(frfile);

    /* loop over detectors in input file */
    ndet = XLALFrameUFrTOCQueryDetectorN(toc);
    for (det = 0; det < ndet; ++det) {
        const char *name;
        LALFrameUFrDetector *detector;
        name = XLALFrameUFrTOCQueryDetectorName(toc, det);
        detector = XLALFrameUFrDetectorRead(frfile, name);
        copydetector(frame, detector);
        XLALFrameUFrDetectorFree(detector);
    }

    XLALFrameUFrTOCFree(toc);
    return 0;
}

int copydetector(LALFrameUFrameH * frame, LALFrameUFrDetector * detector)
{
    const char *name;
    const char *prefix;
    double longitude;
    double latitude;
    double elevation;
    double azimuthx;
    double azimuthy;
    double altitudex;
    double altitudey;
    double midpointx;
    double midpointy;
    int loctime;

    LALFrameUFrDetector *detectorcopy;

    if (!frame || !detector)
        return -1;

    name = XLALFrameUFrDetectorQueryName(detector);
    prefix = XLALFrameUFrDetectorQueryPrefix(detector);
    loctime = XLALFrameUFrDetectorQueryLocalTime(detector);
    longitude = XLALFrameUFrDetectorQueryLongitude(detector);
    latitude = XLALFrameUFrDetectorQueryLatitude(detector);
    elevation = XLALFrameUFrDetectorQueryElevation(detector);
    azimuthx = XLALFrameUFrDetectorQueryArmXAzimuth(detector);
    altitudex = XLALFrameUFrDetectorQueryArmXAltitude(detector);
    midpointx = XLALFrameUFrDetectorQueryArmXMidpoint(detector);
    azimuthy = XLALFrameUFrDetectorQueryArmYAzimuth(detector);
    altitudey = XLALFrameUFrDetectorQueryArmYAltitude(detector);
    midpointy = XLALFrameUFrDetectorQueryArmYMidpoint(detector);

    detectorcopy =
        XLALFrameUFrDetectorAlloc(name, prefix, latitude, longitude,
        elevation, azimuthx, azimuthy, altitudex, altitudey, midpointx,
        midpointy, loctime);
    XLALFrameUFrameHFrDetectorAdd(frame, detectorcopy);
    XLALFrameUFrDetectorFree(detectorcopy);
    return 0;
}

int copychannels(LALFrameUFrameH * frame, LALFrameUFrFile * frfile,
    size_t pos, const char *match)
{
    LALFrameUFrTOC *toc;
    size_t nadc, adc;
    size_t nproc, proc;
    size_t nsim, sim;

    toc = XLALFrameUFrTOCRead(frfile);

    /* loop over channels in input file */

    nadc = XLALFrameUFrTOCQueryAdcN(toc);
    for (adc = 0; adc < nadc; ++adc) {
        const char *name;
        LALFrameUFrChan *channel;
        name = XLALFrameUFrTOCQueryAdcName(toc, adc);
        if (match && strcmp(name, match))
            continue;   /*does not match */
        channel = XLALFrameUFrChanRead(frfile, name, pos);
        copychannel(frame, channel, ADC_CHAN_TYPE);
        XLALFrameUFrChanFree(channel);
    }

    nproc = XLALFrameUFrTOCQueryProcN(toc);
    for (proc = 0; proc < nproc; ++proc) {
        const char *name;
        LALFrameUFrChan *channel;
        name = XLALFrameUFrTOCQueryProcName(toc, proc);
        if (match && strcmp(name, match))
            continue;   /*does not match */
        channel = XLALFrameUFrChanRead(frfile, name, pos);
        copychannel(frame, channel, PROC_CHAN_TYPE);
        XLALFrameUFrChanFree(channel);
    }

    nsim = XLALFrameUFrTOCQuerySimN(toc);
    for (sim = 0; sim < nsim; ++sim) {
        const char *name;
        LALFrameUFrChan *channel;
        name = XLALFrameUFrTOCQuerySimName(toc, sim);
        if (match && strcmp(name, match))
            continue;   /*does not match */
        channel = XLALFrameUFrChanRead(frfile, name, pos);
        copychannel(frame, channel, SIM_CHAN_TYPE);
        XLALFrameUFrChanFree(channel);
    }

    XLALFrameUFrTOCFree(toc);
    return 0;
}

int copychannel(LALFrameUFrameH * frame, LALFrameUFrChan * channel,
    int chantype)
{
    const char *channame;
    const char *vectname;
    int dtype;
    size_t nbytes;
    size_t ndata;
    void *dataorig;
    void *datacopy;
    const char *unity;
    const char *unitx;
    double x0;
    double dx;
    LALFrameUFrChan *chancopy;

    if (!frame || !channel)
        return -1;

    channame = XLALFrameUFrChanQueryName(channel);
    vectname = XLALFrameUFrChanVectorQueryName(channel);

    XLALFrameUFrChanVectorExpand(channel);
    dtype = XLALFrameUFrChanVectorQueryType(channel);
    nbytes = XLALFrameUFrChanVectorQueryNBytes(channel);
    ndata = XLALFrameUFrChanVectorQueryNData(channel);
    x0 = XLALFrameUFrChanVectorQueryStartX(channel, 0);
    dx = XLALFrameUFrChanVectorQueryDx(channel, 0);
    unitx = XLALFrameUFrChanVectorQueryUnitX(channel, 0);
    unity = XLALFrameUFrChanVectorQueryUnitY(channel);
    dataorig = XLALFrameUFrChanVectorQueryData(channel);

    switch (chantype) {
    case ADC_CHAN_TYPE:
        chancopy = XLALFrameUFrAdcChanAlloc(channame, dtype, ndata);
        break;
    case PROC_CHAN_TYPE:
        chancopy =
            XLALFrameUFrProcChanAlloc(channame,
            LAL_FRAMEU_FR_PROC_TYPE_TIME_SERIES,
            LAL_FRAMEU_FR_PROC_SUB_TYPE_UNKNOWN, dtype, ndata);
        break;
    case SIM_CHAN_TYPE:
        chancopy = XLALFrameUFrSimChanAlloc(channame, dtype, ndata);
        break;
    default:
        abort();
    }

    XLALFrameUFrChanVectorAlloc(chancopy, dtype, ndata);
    XLALFrameUFrChanVectorSetName(chancopy, vectname);
    XLALFrameUFrChanVectorSetUnitY(chancopy, unity);
    XLALFrameUFrChanVectorSetUnitX(chancopy, unitx);
    XLALFrameUFrChanVectorSetStartX(chancopy, x0);
    XLALFrameUFrChanVectorSetDx(chancopy, dx);
    datacopy = XLALFrameUFrChanVectorQueryData(chancopy);
    memcpy(datacopy, dataorig, nbytes);
    XLALFrameUFrameHFrChanAdd(frame, chancopy);
    XLALFrameUFrChanFree(chancopy);
    return 0;
}
