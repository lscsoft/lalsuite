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
 * @defgroup lalfr_fmt lalfr-fmt
 * @ingroup lalframe_programs
 *
 * @brief Formats multicolumn text files into channels in a frame file
 *
 * ### Synopsis
 *
 *     lalfr-fmt [file ...]
 *
 * ### Description
 *
 * The `lalfr-fmt` utility reads each multicolumn input `file` file and formats
 * the columns into different channels that are output as a frame file to the
 * standard output.  The `file` operands are processed in command-line order.
 * If `file` is a single dash (`-`) or absent, `lalfr-fmt` reads from the
 * standard input.
 * 
 * If `file` contains only one column of data, each row is interpreted as a
 * sample of a single channel having a sample rate of 16384 Hz.  If there are
 * more than one column in `file` then the first column is interpreted as the
 * GPS time of each sample.  The remaining columns then describe the samples
 * for separate channels.  The output channels are named `C01`, `C02`, etc.
 *
 * If the first line of `file` begins with the `#` character followed
 * tab-separated list of column names with units in parentheses then these
 * are used for the channel names and the sample units.
 *
 * ### Examples
 *
 * The command:
 *
 *     seq 0 0.02 99.98 | paste - - - - - | lalfr-fmt > file.gwf
 *
 * produces a frame file `file.gwf` containing four channels named `C01`,
 * `C02`, `C03`, and `C04`, each with a sample rate of 10 Hz.
 *
 * The command:
 *
 *     lalsim-detector-noise -s 1000000000 -t 64 --aligo-zerodet-highpower | lalfr-fmt > noise.gwf
 *
 * produces a frame file `noise.gwf` containing 64 seconds of simulated aLIGO
 * noise.
 *
 * @sa @ref lalfr_print
 */

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALString.h>
#include <lal/LALFrameU.h>

#define FAILURE(...) do { fprintf(stderr, __VA_ARGS__); exit(99); } while (0)

char *mystrdup(const char *s);
int readdata(size_t * ndim, size_t * size, double **data, char ***channames, char ***chanunits, FILE * fp);

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
        char **channames;
        char **chanunits;
        FILE *fp;

        if (strcmp(fname, "-") == 0)
            fp = stdin;
        else
            fp = fopen(fname, "r");
        if (!fp)
            FAILURE("could not open file %s\n", fname);

        readdata(&ndim, &dimlen, &data, &channames, &chanunits, fp);

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
            size_t i;

            channel = XLALFrameUFrProcChanAlloc(channames[dim], type, subtype, dtype, dimlen);
            XLALFrameUFrChanVectorAlloc(channel, dtype, dimlen);
            XLALFrameUFrChanVectorSetName(channel, channames[dim]);
            XLALFrameUFrChanVectorSetDx(channel, 1.0 / srate);
            XLALFrameUFrChanVectorSetUnitX(channel, "s");
            XLALFrameUFrChanVectorSetUnitY(channel, chanunits[dim]);

            vect = XLALFrameUFrChanVectorQueryData(channel);
            for (i = 0; i < dimlen; ++i)
                vect[i] = data[i * ndim + dim];

            XLALFrameUFrameHFrChanAdd(frame, channel);
            XLALFrameUFrChanFree(channel);
            free(channames[dim]);
            free(chanunits[dim]);
        }
        free(data);
        free(channames);
        free(chanunits);

        XLALFrameUFrameHWrite(output, frame);
        XLALFrameUFrameHFree(frame);
    }

    XLALFrameUFrFileClose(output);
    return 0;
}

int readdata(size_t * ndim, size_t * dimlen, double **data, char ***channames, char ***chanunits, FILE * fp)
{
    char header[LINE_MAX] = "";
    char line[LINE_MAX];
    char *hdr = header;
    char *tok;
    size_t block = 1024;
    size_t bufsz = 0;
    size_t rows = 0;
    size_t cols = 0;
    size_t col;
	int c;

    /* determine if the first line contains channel names and units */
    if ((c = fgetc(fp)) == '#')
        fgets(header, sizeof(header), fp); /* save header line */
    else
        ungetc(c, fp);

    *data = NULL;
    *channames = NULL;
    *chanunits = NULL;

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
    *channames = calloc(cols, sizeof(**channames));
    *chanunits = calloc(cols, sizeof(**chanunits));

    /* now scan header for channel names and units */
    for (col = 0; col < cols; ++col) {
        char name[LINE_MAX] = "";
        char unit[LINE_MAX] = "";
        int n = 0;
        if (hdr && (tok = XLALStringToken(&hdr, "\t", 0)))
            n = sscanf(tok, "%s (%[^)])", name, unit);
        if (n < 1)
            snprintf(name, sizeof(name), "C%02zu", col);
        (*channames)[col] = mystrdup(name);
        (*chanunits)[col] = mystrdup(unit);
    }

    return 0;
}

/* strdup is not C99... */
char *mystrdup(const char *str)
{
    size_t len = strlen(str);
    char *dup = malloc(len + 1);
    strcpy(dup, str);
    dup[len] = '\0';
    return dup;
}
