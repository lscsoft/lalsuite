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

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALFrameU.h>

#define FAILURE(...) do { fprintf(stderr, __VA_ARGS__); exit(99); } while (0)

int readdata(size_t * ndim, size_t * size, double **data, FILE * fp);

int main(int argc, char *argv[])
{
    char stdio[] = "-";
    char *defaultfilev[1] = { stdio };
    int filec = (argc == 1 ? 1 : argc - 1);
    char **filev = (argc == 1 ? defaultfilev : &argv[1]);
    LALFrameUFrFile *output;
    int f;

    output = XLALFrameUFrFileOpen(stdio, "w");
    if (!output)
        FAILURE("could not create output file\n");

    for (f = 0; f < filec; ++f) {
        char *fname = filev[f];
        LALFrameUFrameH *frame;
        LALFrameUFrChan *channel;
        double srate = 16384.0; /* default sample rate */
        double t0 = 0.0;        /* default gps start time */
        double dt;
        int type = LAL_FRAMEU_FR_PROC_TYPE_TIME_SERIES;
        int subtype = LAL_FRAMEU_FR_PROC_SUB_TYPE_UNKNOWN;
        int dtype = LAL_FRAMEU_FR_VECT_8R;
        size_t ndim;
        size_t dim;
        size_t dimlen;
        double *data;
        double *vect;
        FILE *fp;

        if (strcmp(fname, "-") == 0)
            fp = stdin;
        else
            fp = fopen(fname, "r");
        if (!output)
            FAILURE("could not open file %s\n", fname);

        readdata(&ndim, &dimlen, &data, fp);

        if (ndim == 1)  /* use default sample rate and start time */
            dt = dimlen / srate;
        else {  /* first dimension is gps time */
            /* assume sample rate in Hz is an integer */
            srate =
                floor(0.5 + (dimlen - 1) / (data[(dimlen - 1) * ndim] -
                    data[0]));
            dt = dimlen / srate;
            t0 = data[0];
        }

        frame = XLALFrameUFrameHAlloc(fname, t0, dt, f);
        XLALFrameUFrameHSetRun(frame, 0);

        for (dim = (ndim == 1 ? 0 : 1); dim < ndim; ++dim) {
            char channame[256];
            size_t i;

            snprintf(channame, sizeof(channame), "C%02zu", dim);
            channel =
                XLALFrameUFrProcChanAlloc(channame, type, subtype, dtype,
                dimlen);
            XLALFrameUFrChanVectorAlloc(channel, dtype, dimlen);
            XLALFrameUFrChanVectorSetName(channel, channame);
            XLALFrameUFrChanVectorSetDx(channel, 1.0 / srate);
            XLALFrameUFrChanVectorSetUnitX(channel, "s");

            vect = XLALFrameUFrChanVectorQueryData(channel);
            for (i = 0; i < dimlen; ++i)
                vect[i] = data[i * ndim + dim];

            XLALFrameUFrameHFrChanAdd(frame, channel);
            XLALFrameUFrChanFree(channel);
        }
        free(data);

        XLALFrameUFrameHWrite(output, frame);
        XLALFrameUFrameHFree(frame);
    }

    XLALFrameUFrFileClose(output);
    return 0;
}

int readdata(size_t * ndim, size_t * dimlen, double **data, FILE * fp)
{
    char line[LINE_MAX];
    size_t block = 1024;
    size_t bufsz = 0;
    size_t rows = 0;
    size_t cols = 0;
    size_t col;

    *data = NULL;
    while (fgets(line, sizeof(line), fp)) {
        char *s;
        char *endp;
        if (*line == '#')       /* ignore lines beginning with '#' */
            continue;
        if (cols == 0) {        /* count columns on first line */
            endp = line;
            while (1) {
                s = endp;
                /* work around bug in glibc < 2.16
                 * http://sourceware.org/bugzilla/show_bug.cgi?id=13970 */
                double v = strtod(s, &endp);
                (void)v;
                if (s == endp || *endp == '\0')
                    break;
                ++cols;
            }
            if (cols == 0)
                FAILURE("format error on input line %zu\n", rows + 1);
        }
        endp = line;

        if (rows == bufsz) {
            bufsz += block;
            *data = realloc(*data, bufsz * cols * sizeof(**data));
        }
        for (col = 0; col < cols; ++col) {
            s = endp;
            (*data)[rows * cols + col] = strtod(s, &endp);
            if (s == endp || *endp == '\0')
                FAILURE("format error on input line %zu\n", rows + 1);
        }
        ++rows;
    }

    *data = realloc(*data, rows * cols * sizeof(**data));
    *ndim = cols;
    *dimlen = rows;

    return 0;
}
