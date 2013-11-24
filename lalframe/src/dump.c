/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Robert Adam Mercer
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
#include <lal/LALFrameU.h>

#define FAILURE(...) do { fprintf(stderr, __VA_ARGS__); exit(1); } while (0)

int quotelvl = 0;
int indent(int);
#define indent indent(quotelvl)

int dumpdetector(LALFrameUFrFile * frfile, size_t det);
int dumpframe(LALFrameUFrFile * frfile, size_t pos);
int dumpchannel(LALFrameUFrChan * channel, size_t chan, const char *chantype);
const char *typestr(int type);
const char *compressstr(int type);

int main(int argc, char *argv[])
{
    char stdio[] = "-";
    char *defaultfilev[1] = { stdio };
    int filec = (argc == 1 ? 1 : argc - 1);
    char **filev = (argc == 1 ? defaultfilev : &argv[1]);
    int f;

    for (f = 0; f < filec; ++f) {
        char *fname = filev[f];
        LALFrameUFrFile *frfile = NULL;
        LALFrameUFrTOC *toc = NULL;
        size_t nframe = 0;
        size_t ndet = 0;
        size_t det = 0;
        size_t pos = 0;

        frfile = XLALFrameUFrFileOpen(fname, "r");
        if (!frfile)
            FAILURE("file %s not found\n", fname);

        if (f > 0)
            printf("\n");
        if (filec > 1)
            printf("==> %s <==\n", fname);

        toc = XLALFrameUFrTOCRead(frfile);
        if (!toc)
            FAILURE("no TOC found\n");

        /* loop over detectors in file */
        ndet = XLALFrameUFrTOCQueryDetectorN(toc);
        for (det = 0; det < ndet; ++det)
            dumpdetector(frfile, det);

        /* loop over frames in file */
        nframe = XLALFrameUFrTOCQueryNFrame(toc);
        if ((int)(nframe) <= 0)
            FAILURE("no frames found\n");
        for (pos = 0; pos < nframe; ++pos)
            dumpframe(frfile, pos);

        XLALFrameUFrTOCFree(toc);
        XLALFrameUFrFileClose(frfile);
    }

    return 0;
}

int dumpdetector(LALFrameUFrFile * frfile, size_t det)
{
    LALFrameUFrTOC *toc;
    LALFrameUFrDetector *detector;
    const char *name;
    const char *prefix;
    /*
     * double longitude;
     * double latitude;
     * double elevation;
     * double azimuthx;
     * double azimuthy;
     * double altitudex;
     * double altitudey;
     * double midpointx;
     * double midpointy;
     */
    int loctime;

    if (!frfile)
        return -1;

    toc = XLALFrameUFrTOCRead(frfile);
    name = XLALFrameUFrTOCQueryDetectorName(toc, det);
    detector = XLALFrameUFrDetectorRead(frfile, name);
    if (!detector) {
        --quotelvl;
        return -1;
    }

    name = XLALFrameUFrDetectorQueryName(detector);
    prefix = XLALFrameUFrDetectorQueryPrefix(detector);
    loctime = XLALFrameUFrDetectorQueryLocalTime(detector);
    /*
     * longitude = XLALFrameUFrDetectorQueryLongitude(detector);
     * latitude = XLALFrameUFrDetectorQueryLatitude(detector);
     * elevation = XLALFrameUFrDetectorQueryElevation(detector);
     * azimuthx = XLALFrameUFrDetectorQueryArmXAzimuth(detector);
     * altitudex = XLALFrameUFrDetectorQueryArmXAltitude(detector);
     * midpointx = XLALFrameUFrDetectorQueryArmXMidpoint(detector);
     * azimuthy = XLALFrameUFrDetectorQueryArmYAzimuth(detector);
     * altitudey = XLALFrameUFrDetectorQueryArmYAltitude(detector);
     * midpointy = XLALFrameUFrDetectorQueryArmYMidpoint(detector);
     */

    indent;
    printf("FrDetector %zu %s (%s): local time = %d\n", det, name, prefix,
        loctime);

    XLALFrameUFrTOCFree(toc);
    XLALFrameUFrDetectorFree(detector);
    return 0;
}

int dumpframe(LALFrameUFrFile * frfile, size_t pos)
{
    LALFrameUFrameH *frame;
    LALFrameUFrTOC *toc;
    const char *name;
    int run;
    int frnum;
    int dq;
    double tip;
    double tfp;
    int leaps;
    double dt;
    size_t nadc;
    size_t nsim;
    size_t nproc;
    size_t adc;
    size_t sim;
    size_t proc;

    if (!frfile)
        return -1;

    /* read frame */
    frame = XLALFrameUFrameHRead(frfile, pos);
    if (!frame)
        FAILURE("unable to open frame %zu\n", pos);

    /* get frame metadata */
    name = XLALFrameUFrameHQueryName(frame);
    run = XLALFrameUFrameHQueryRun(frame);
    frnum = XLALFrameUFrameHQueryFrame(frame);
    dq = XLALFrameUFrameHQueryDataQuality(frame);
    tfp = XLALFrameUFrameHQueryGTimeModf(&tip, frame);
    leaps = XLALFrameUFrameHQueryULeapS(frame);
    dt = XLALFrameUFrameHQueryDt(frame);

    /* print frame metadata */
    indent;
    printf("FrameH %zu", pos);
    if (name && *name)
        printf(" %s", name);
    printf(" run %d, frame %d", run, frnum);
    printf(": dq = %d, ", dq);
    if (tfp == 0.0)
        printf("t0 = %d s, ", (int)tip);
    else
        printf("t0 = %f s, ", tip + tfp);
    printf("dt = %g s, ", dt);
    printf("TAI-UTC = %d\n", leaps);

    /* get channel counts */
    toc = XLALFrameUFrTOCRead(frfile);
    nadc = XLALFrameUFrTOCQueryAdcN(toc);
    nsim = XLALFrameUFrTOCQuerySimN(toc);
    nproc = XLALFrameUFrTOCQueryProcN(toc);

    /* print adc channel info */
    for (adc = 0; adc < nadc; ++adc) {
        LALFrameUFrChan *channel;
        name = XLALFrameUFrTOCQueryAdcName(toc, adc);
        channel = XLALFrameUFrChanRead(frfile, name, pos);
        dumpchannel(channel, adc, "Adc");
        XLALFrameUFrChanFree(channel);
    }

    /* print sim channel info */
    for (sim = 0; sim < nsim; ++sim) {
        LALFrameUFrChan *channel;
        name = XLALFrameUFrTOCQuerySimName(toc, sim);
        channel = XLALFrameUFrChanRead(frfile, name, pos);
        dumpchannel(channel, sim, "Sim");
        XLALFrameUFrChanFree(channel);
    }

    /* print proc channel info */
    for (proc = 0; proc < nproc; ++proc) {
        LALFrameUFrChan *channel;
        name = XLALFrameUFrTOCQueryProcName(toc, proc);
        channel = XLALFrameUFrChanRead(frfile, name, pos);
        dumpchannel(channel, proc, "Proc");
        XLALFrameUFrChanFree(channel);
    }

    XLALFrameUFrTOCFree(toc);
    XLALFrameUFrameHFree(frame);
    return 0;
}

int dumpchannel(LALFrameUFrChan * channel, size_t chan, const char *chantype)
{
    const char *channame;
    double toffset;
    const char *vectname;
    int compress;
    size_t nbytes;
    size_t ndata;
    int type;
    size_t ndim;
    size_t dim;
    const char *unity;

    if (!channel)
        return -1;

    ++quotelvl;

    channame = XLALFrameUFrChanQueryName(channel);
    toffset = XLALFrameUFrChanQueryTimeOffset(channel);

    vectname = XLALFrameUFrChanVectorQueryName(channel);
    compress = XLALFrameUFrChanVectorQueryCompress(channel);
    nbytes = XLALFrameUFrChanVectorQueryNBytes(channel);
    ndata = XLALFrameUFrChanVectorQueryNData(channel);
    type = XLALFrameUFrChanVectorQueryType(channel);
    ndim = XLALFrameUFrChanVectorQueryNDim(channel);
    unity = XLALFrameUFrChanVectorQueryUnitY(channel);

    /* print channel information */
    indent;
    printf("Fr%sData %zu %s", chantype, chan, channame);
    if (toffset > 0.0)
        printf(", offset = %g s", toffset);
    if (strcmp(channame, vectname))
        printf(", FrVect %s", vectname);
    printf(": ");
    if (nbytes / 1073741824)
        printf("%2gGi", floor(0.5 + nbytes / 1073741824.0));
    else if (nbytes / 1048576)
        printf("%2gMi", floor(0.5 + nbytes / 1048576.0));
    else if (nbytes / 1024)
        printf("%2gKi", floor(0.5 + nbytes / 1024.0));
    else
        printf("%zuB", nbytes);
    printf(" (%s)", compressstr(compress));
    printf(", %zu %s pts", ndata, typestr(type));
    if (unity && *unity)
        printf(", yunits = %s", unity);
    if (ndim == 1) {
        double dx;
        double x0;
        const char *unitx;
        dx = XLALFrameUFrChanVectorQueryDx(channel, 0);
        x0 = XLALFrameUFrChanVectorQueryStartX(channel, 0);
        unitx = XLALFrameUFrChanVectorQueryUnitX(channel, 0);
        printf(", x0 = %g, dx = %g", x0, dx);
        if (unitx && *unitx)
            printf(", xunits = %s", unitx);
        printf("\n");
    } else {
        printf(", %zu dim\n", ndim);
        ++quotelvl;
        for (dim = 0; dim < ndim; ++dim) {
            size_t nx;
            double dx;
            double x0;
            const char *unitx;

            nx = XLALFrameUFrChanVectorQueryNx(channel, dim);
            dx = XLALFrameUFrChanVectorQueryDx(channel, dim);
            x0 = XLALFrameUFrChanVectorQueryStartX(channel, dim);
            unitx = XLALFrameUFrChanVectorQueryUnitX(channel, dim);
            indent;
            printf("nx[%zu] = %zu", dim, nx);
            printf(", x0[%zu] = %g", dim, x0);
            printf(", dx[%zu] = %g", dim, dx);
            if (unitx && *unitx)
                printf(", xunits[%zu] = %s", dim, unitx);
            printf("\n");
        }
        --quotelvl;
    }
    --quotelvl;
    return 0;
}

const char *typestr(int type)
{
    switch (type) {
    case LAL_FRAMEU_FR_VECT_C:
        return "char";
    case LAL_FRAMEU_FR_VECT_2S:
        return "int16_t";
    case LAL_FRAMEU_FR_VECT_4S:
        return "int32_t";
    case LAL_FRAMEU_FR_VECT_8S:
        return "int64_t";
    case LAL_FRAMEU_FR_VECT_1U:
        return "uint8_t";
    case LAL_FRAMEU_FR_VECT_2U:
        return "uint16_t";
    case LAL_FRAMEU_FR_VECT_4U:
        return "uint32_t";
    case LAL_FRAMEU_FR_VECT_8U:
        return "uint64_t";
    case LAL_FRAMEU_FR_VECT_4R:
        return "float";
    case LAL_FRAMEU_FR_VECT_8R:
        return "double";
    case LAL_FRAMEU_FR_VECT_8C:
        return "float complex";
    case LAL_FRAMEU_FR_VECT_16C:
        return "double complex";
    case LAL_FRAMEU_FR_VECT_STRING:
        return "string";
    default:
        break;
    }
    return "unknown";
}

const char *compressstr(int compress)
{
    switch (compress) {
    case LAL_FRAMEU_FR_VECT_COMPRESS_RAW:
        return "RAW";
    case LAL_FRAMEU_FR_VECT_COMPRESS_GZIP:
        return "GZIP compression";
    case LAL_FRAMEU_FR_VECT_COMPRESS_DIFF_GZIP:
        return "DIFF_GZIP compression";
    case LAL_FRAMEU_FR_VECT_COMPRESS_ZERO_SUPPRESS_WORD_2:
        return "ZERO_SUPPRESS_WORD_2 compression";
    case LAL_FRAMEU_FR_VECT_COMPRESS_ZERO_SUPPRESS_WORD_4:
        return "ZERO_SUPPRESS_WORD_4 compression";
    default:
        break;
    }
    return "unknown compression";
}

/* output quote formatting */
#undef indent
int indent(int n)
{
    const char *tab = "  ";
    switch (n) {
    default:
    case 3:
        fputs(tab, stdout);
    case 2:
        fputs(tab, stdout);
        fputs(tab, stdout);
        fputs("... ", stdout);
        break;
    case 1:
        fputs(tab, stdout);
        fputs("- ", stdout);
        break;
    case 0:
        fputs("- ", stdout);
        break;
    }
    return 0;
}
