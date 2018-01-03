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

/**
 * @defgroup lalfr_print lalfr-print
 * @ingroup lalframe_programs
 *
 * @brief Prints channel data from frame files
 *
 * ### Synopsis
 *
 *     lalfr-print [file ...]
 *
 * ### Description
 *
 * The `lalfr-print` utility reads the contents of the channels from each
 * `file` and prints the data to the standard output.  If `file` is a single
 * dash (`-`) or absent, `lalfr-print` reads from the standard input.
 *
 * For each channel contained in `file`, the output is written in two-column
 * format where the first column is the GPS time of each sample and the second
 * column contains the corresponding sample values.  The columns are separated
 * by a tab character (`\t`) and each line is separated by a newline character
 * (`\n`).
 *
 * Each channel in `file` is written sequentially and are separated by a line
 * containing a separator consisting of the character `#` followed by the name
 * of the next channel to be written.  If more than one file argument is
 * present then the separate files are processed sequentially and separated in
 * the output by a line containing the separator `==> file <==` where `file` is
 * the name of the current file being processed.
 *
 * ### Examples
 *
 * The command:
 *
 *     lalfr-print file.gwf
 *
 * prints to standard output all the channels contained in `file.gwf`.  If
 * there are more than one channel present in that file, they can be split into
 * separate files, each containing a single channel's data, with the command:
 *
 *     lalfr-print file.gwf | split -p '#'
 *
 * and the resulting files `xaa`, `xab`, etc., contain the data from each of
 * the channels in `file.gwf`.
 *
 * @sa @ref lalfr_dump, @ref lalfr_fmt
 */

#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALFrameU.h>

#define FS "\t" /* field separator */
#define RS "\n" /* record separator */

#define FAILURE(...) do { fprintf(stderr, __VA_ARGS__); exit(1); } while (0)

int printchannel(LALFrameUFrChan * channel, double x0, int fdom);
int printval(void *data, size_t i, int dtype);

int main(int argc, char *argv[])
{
    char stdio[] = "-";
    char *defaultfilev[1] = { stdio };
    int filec = (argc == 1 ? 1 : argc - 1);
    char **filev = (argc == 1 ? defaultfilev : &argv[1]);
    int f;

    for (f = 0; f < filec; ++f) {
        char *fname = filev[f];
        LALFrameUFrFile *frfile;
        LALFrameUFrTOC *toc;
        size_t nadc;
        size_t nsim;
        size_t nproc;
        size_t nframe;
        size_t chan;
        size_t pos;

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

        nframe = XLALFrameUFrTOCQueryNFrame(toc);
        if (!toc)
            FAILURE("no frames found\n");

        nadc = XLALFrameUFrTOCQueryAdcN(toc);
        nsim = XLALFrameUFrTOCQuerySimN(toc);
        nproc = XLALFrameUFrTOCQueryProcN(toc);

        for (chan = 0; chan < nadc; ++chan) {
            LALFrameUFrChan *channel;
            const char *name;
            name = XLALFrameUFrTOCQueryAdcName(toc, chan);
            printf("# time (s)\t%s", name);
            for (pos = 0; pos < nframe; ++pos) {
                double tip;
                double tfp;
                tfp = XLALFrameUFrTOCQueryGTimeModf(&tip, toc, pos);
                channel = XLALFrameUFrChanRead(frfile, name, pos);
                if (pos == 0) {
                    const char *unit = XLALFrameUFrChanVectorQueryUnitY(channel);
                    if (unit && strlen(unit))
                        printf(" (%s)", unit);
                    printf("\n");
                }
                printchannel(channel, tip + tfp, 0);
                XLALFrameUFrChanFree(channel);
            }
        }

        for (chan = 0; chan < nsim; ++chan) {
            LALFrameUFrChan *channel;
            const char *name;
            name = XLALFrameUFrTOCQuerySimName(toc, chan);
            printf("# time (s)\t%s", name);
            for (pos = 0; pos < nframe; ++pos) {
                double tip;
                double tfp;
                tfp = XLALFrameUFrTOCQueryGTimeModf(&tip, toc, pos);
                channel = XLALFrameUFrChanRead(frfile, name, pos);
                if (pos == 0) {
                    const char *unit = XLALFrameUFrChanVectorQueryUnitY(channel);
                    if (unit && strlen(unit))
                        printf(" (%s)", unit);
                    printf("\n");
                }
                printchannel(channel, tip + tfp, 0);
                XLALFrameUFrChanFree(channel);
            }
        }

        for (chan = 0; chan < nproc; ++chan) {
            LALFrameUFrChan *channel;
            const char *name;
            int fdom = 0;
            name = XLALFrameUFrTOCQueryProcName(toc, chan);
            for (pos = 0; pos < nframe; ++pos) {
                double tip;
                double tfp;
                tfp = XLALFrameUFrTOCQueryGTimeModf(&tip, toc, pos);
                channel = XLALFrameUFrChanRead(frfile, name, pos);
                if (pos == 0) {
                    const char *unit;
                    unit = XLALFrameUFrChanVectorQueryUnitX(channel, 0);
                    if (unit && strcmp(unit, "s") == 0)
                        printf("# time (s)\t%s", name);
                    else if (unit && strcmp(unit, "s^-1") == 0) {
                        printf("# freq (s^-1)\t%s", name);
                        fdom = 1;
                    }
                    else if (unit && strlen(unit)) {
                        printf("# sample (%s)\t%s", unit, name);
                        fdom = 1;
                    }
                    else {
                        printf("# sample\t%s", name);
                        fdom = 1;
                    }
                    unit = XLALFrameUFrChanVectorQueryUnitY(channel);
                    if (unit && strlen(unit))
                        printf(" (%s)", unit);
                    printf("\n");
                }
                if (fdom)
                    printchannel(channel, 0.0, fdom);
                else
                    printchannel(channel, tip + tfp, 0);
                XLALFrameUFrChanFree(channel);
            }
        }

        XLALFrameUFrTOCFree(toc);
        XLALFrameUFrFileClose(frfile);
    }

    return 0;
}

int printchannel(LALFrameUFrChan * channel, double x0, int fdom)
{
    /* const char *name; */
    double x0ip, x0fp;
    double dx;
    void *data;
    int dtype;
    size_t ndata;
    size_t i;

    if (!channel)
        return -1;

    XLALFrameUFrChanVectorExpand(channel);

    /* name = XLALFrameUFrChanQueryName(channel); */
    if (!fdom)
        x0 += XLALFrameUFrChanQueryTimeOffset(channel);

    ndata = XLALFrameUFrChanVectorQueryNData(channel);
    x0 += XLALFrameUFrChanVectorQueryStartX(channel, 0);
    dx = XLALFrameUFrChanVectorQueryDx(channel, 0);
    data = XLALFrameUFrChanVectorQueryData(channel);
    dtype = XLALFrameUFrChanVectorQueryType(channel);

    x0fp = modf(x0, &x0ip);
    for (i = 0; i < ndata; ++i) {
        double xip, xfp;
        xfp = modf(x0fp + i * dx, &xip);
        xip += x0ip;
        printf("%ld.%09ld", (long)xip, (long)fabs(1e9 * xfp));
        fputs(FS, stdout);
        printval(data, i, dtype);
        fputs(RS, stdout);
    }

    return 0;
}

int printval(void *data, size_t i, int dtype)
{
    switch (dtype) {
    case LAL_FRAMEU_FR_VECT_C:
        return printf("%c", ((char *)data)[i]);
    case LAL_FRAMEU_FR_VECT_2S:
        return printf("%" PRIi16, ((int16_t *) data)[i]);
    case LAL_FRAMEU_FR_VECT_4S:
        return printf("%" PRIi32, ((int32_t *) data)[i]);
    case LAL_FRAMEU_FR_VECT_8S:
        return printf("%" PRIi64, ((int64_t *) data)[i]);
    case LAL_FRAMEU_FR_VECT_1U:
        return printf("%" PRIu8, ((uint8_t *) data)[i]);
    case LAL_FRAMEU_FR_VECT_2U:
        return printf("%" PRIu16, ((uint16_t *) data)[i]);
    case LAL_FRAMEU_FR_VECT_4U:
        return printf("%" PRIu32, ((uint32_t *) data)[i]);
    case LAL_FRAMEU_FR_VECT_8U:
        return printf("%" PRIu64, ((uint64_t *) data)[i]);
    case LAL_FRAMEU_FR_VECT_4R:
        return printf("%e", (double)((float *)data)[i]);
    case LAL_FRAMEU_FR_VECT_8R:
        return printf("%e", ((double *)data)[i]);
    case LAL_FRAMEU_FR_VECT_8C:
        return printf("(%e,%e)", (double)((float *)data)[2 * i],
            (double)((float *)data)[2 * i + 1]);
    case LAL_FRAMEU_FR_VECT_16C:
        return printf("(%e,%e)", ((float *)data)[2 * i],
            ((float *)data)[2 * i + 1]);
    case LAL_FRAMEU_FR_VECT_STRING:
        return printf("%s", ((char **)data)[i]);
    default:
        break;
    }
    return -1;
}
