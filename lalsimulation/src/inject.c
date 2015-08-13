/*
*  Copyright (C) 2015 Jolien Creighton
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
 * @defgroup lalsim_inject lalsim-inject
 * @ingroup lalsimulation_programs
 *
 * @brief Injects an induced gravitational wave strain into detector data
 *
 * ### Synopsis
 *
 *     lalsim-inject [options] targetfile injectfile1 [injectfile2 ...]
 *
 * ### Description
 *
 * The `lalsim-inject` utility takes gravitational wave detector data
 * contained in `targetfile` and adds to it the gravitational wave
 * induced strain data in `injectfile1` ....  The result is written
 * to standard output.  All input and output is two-column ascii format
 * data where the first column contains the GPS timestamp of each sample
 * and the second column contains the detector strain value.
 *
 * ### Options
 *
 * <DL>
 * <DT>`-h`, `--help`    <DD>print a help message and exit</DD>
 * <DT>`-v`, `--verbose` <DD>verbose output</DD>
 * </DL>
 *
 * ### Environment
 *
 * The `LAL_DEBUG_LEVEL` can used to control the error and warning reporting of
 * `lalsim-inject`.  Common values are: `LAL_DEBUG_LEVEL=0` which suppresses
 * error messages, `LAL_DEBUG_LEVEL=1`  which prints error messages alone,
 * `LAL_DEBUG_LEVEL=3` which prints both error messages and warning messages,
 * and `LAL_DEBUG_LEVEL=7` which additionally prints informational messages.
 *
 * ### Exit Status
 *
 * The `lalsim-inject` utility exits 0 on success, and >0 if an error occurs.
 *
 * ### Example
 *
 * The following set of commands produces 16 seconds of simulated detector
 * noise for the LHO detector starting at GPS time 1000000000; produces a
 * synthetic binary neutron star signal in the LHO detector that has a
 * geocentric end time of 1000000008; and adds the signal to the noise:
 *
 *
 *     lalsim-detector-noise --aligo-zerodet-highpower -s 1000000000 -t 16 > noise
 *     lalsim-inspiral | lalsim-detector-strain -D H1 -a 1:23:45 -d 45.0 -p 30.0 -t 1000000008 > signal
 *     lalsim-inject noise signal > output
 *
 * The resulting file `output` contains the simulated signal contained in file
 * `signal` injected into the simulated noise contained in file `noise`.
 */


#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/LALSimulation.h>

char *program;
int usage(void);
int isoption(char *arg, const char *opt);
REAL8TimeSeries *readdata(FILE * fp);

int main(int argc, char *argv[])
{
    char tstr[32];      // string to hold GPS time -- 31 characters is enough
    REAL8TimeSeries **h;
    int verbose = 0;
    int first_stdin = -1;
    int n;
    int c;
    size_t j;

    program = argv[0];

    /* parse options */
    for (c = 1; c < argc; ++c) {
        char *arg = argv[c];
        if (*arg != '-' || strcmp(arg, "-") == 0)       /* end of options */
            break;
        if (strcmp(arg, "--") == 0) {   /* stop option parsing */
            ++c;        /* go on to next argument */
            break;
        }
        if (isoption(arg, "-help") || isoption(arg, "--help")) {
            usage();
            return 0;
        }
        if (isoption(arg, "-verbose") || isoption(arg, "--verbose"))
            verbose = 1;
    }

    /* parse remaining arguments, reading the various time series */

    if (argc - c < 2) {
        fprintf(stderr, "%s: must provide one target file and at least one injection file\n", program);
        usage();
        return 1;
    }

    h = calloc(argc - c, sizeof(*h));
    for (n = 0; c < argc; ++c, ++n) {
        char *fname = argv[c];
        FILE *fp = stdin;
        if (verbose)
            fprintf(stderr, "%s: reading %s data ", program, n ? "injection" : "target");
        if (strcmp(fname, "-") != 0) {
            if (verbose)
                fprintf(stderr, "from file %s\n", fname);
            if (!(fp = fopen(fname, "r"))) {
                fprintf(stderr, "%s: could not open file %s\n", program, fname);
                return 1;
            }
        } else if (verbose)
            fprintf(stderr, "from <stdin>\n");

        if (fp == stdin) {
            if (first_stdin < 0) {
                h[n] = readdata(fp);
                first_stdin = n;
            } else      /* already read stdin: copy data */
                h[n] = XLALCutREAL8TimeSeries(h[first_stdin], 0, h[first_stdin]->data->length);
        } else
            h[n] = readdata(fp);

        if (!h[n]) {
            fprintf(stderr, "%s: failed to read data from file %s\n", program, fname);
        }
        if (verbose)
            fprintf(stderr, "%s: %d points of strain data read\n", program, (int)h[n]->data->length);
    }

    /* add injections to target */
    for (c = 1; c < n; ++c) {
        /* sample rates deduced from file might have roundoff
         * errors: if an injection series' deltaT is close enough 
         * the target series' deltaT, make it equal */
        if (fabs(1.0 / h[c]->deltaT - 1.0 / h[0]->deltaT) < 0.1)
            h[c]->deltaT = h[0]->deltaT;
        else {
            fprintf(stderr, "%s: incorrect sample rate for injection %d -- must match target\n", program, c);
            return 1;
        }

        /* now add */
        if (verbose)
            fprintf(stderr, "%s: adding injection %d to target\n", program, c);
        if (XLALSimAddInjectionREAL8TimeSeries(h[0], h[c], NULL) < 0) {
            fprintf(stderr, "%s: failed to add injection %d to target\n", program, c);
            return 1;
        }
    }

    /* output results */
    fprintf(stdout, "# time (s)\tSTRAIN (strain)\n");
    for (j = 0; j < h[0]->data->length; ++j) {
        LIGOTimeGPS t = h[0]->epoch;
        fprintf(stdout, "%s\t%e\n", XLALGPSToStr(tstr, XLALGPSAdd(&t, j * h[0]->deltaT)), h[0]->data->data[j]);
    }

    /* cleanup and exit */
    while (n--)
        XLALDestroyREAL8TimeSeries(h[n]);
    free(h);
    LALCheckMemoryLeaks();
    return 0;
}

/* determine if argument matches option name */
int isoption(char *arg, const char *opt)
{
    return strstr(opt, arg) == opt;
}

/* read a file containing timeseries data */
REAL8TimeSeries *readdata(FILE * fp)
{
    REAL8TimeSeries *h;
    const size_t block = 1024;
    LIGOTimeGPS start;
    LIGOTimeGPS end;
    double dt;
    double *data = NULL;
    size_t bufsz = 0;
    size_t n, l;
    char line[LINE_MAX];
    char t0[LINE_MAX];
    char t1[LINE_MAX];

    for (l = 0, n = 0; fgets(line, sizeof(line), fp); ++l) {
        int c;
        if (*line == '#')
            continue;
        if (n == bufsz) {       /* allocate more memory */
            bufsz += block;
            data = realloc(data, bufsz * sizeof(*data));
        }
        c = sscanf(line, "%s %le", n ? t1 : t0, data + n);
        if (c != 2) {
            fprintf(stderr, "%s: format error on line %zd: %s\n", program, l, line);
            exit(1);
        }
        ++n;
    }
    data = realloc(data, n * sizeof(*data));

    XLALStrToGPS(&start, t0, NULL);
    XLALStrToGPS(&end, t1, NULL);
    dt = XLALGPSDiff(&end, &start) / (n - 1);
    h = XLALCreateREAL8TimeSeries("strain", &start, 0.0, dt, &lalStrainUnit, n);
    memcpy(h->data->data, data, n * sizeof(*data));

    free(data);
    return h;
}

int usage(void)
{
    /* *INDENT-OFF* */
    fprintf(stderr, "\
usage: %s [options] targetfile injectfile1 [injectfile2 ...]\n\
options:\n\
        -h, --help     print this message and exit\n\
        -v, --verbose  verbose output\n\
", program);
    /* *INDENT-ON* */
    return 0;
}
