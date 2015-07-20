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

/**
 * @defgroup lalfr_vis lalfr-vis
 * @ingroup lalframe_programs
 *
 * @brief Visualize frame data
 *
 * ### Synopsis
 *
 *     lalfr-vis --channel=channel --frame-cache=cachefile --start-time=tstart --duration=deltat [--output=outfile] [--highpass=minfreq] [--lowpass=maxfreq] [--pad=padding] [--resample=srate] [--spectrum=resolution]
 *
 *     lalfr-vis --channel=channel --frame-glob=globstring --start-time=tstart --duration=deltat [--output=outfile] [--highpass=minfreq] [--lowpass=maxfreq] [--pad=padding] [--resample=srate] [--spectrum=resolution]
 *
 *
 * ### Description
 *
 * The `lalfr-vis` utility reads a requested interval
 * [`tstart`,`tstart+deltat`) of `channel` data from frame files that are
 * either indexed in the `cachefile` or matching the pattern `globstring` as
 * described by `glob(3)`.  The output is written to `outfile` and the format
 * of the output is deter- mined by extension of `outfile` as described below.
 * If `outfile` is not specified, the output is written to the standard output
 * in two-column ascii format data.
 *
 * The `lalfr-vis` can optionally perform certain manipulations of the data
 * that is read, including:
 *
 * * High-pass filtering of the data.
 * * Low-pass filtering of the data.
 * * Resampling of the data.
 * * Computing the power spectrum of the data.
 *
 * If any of the filtering or resampling operations are performed, it is
 * recommended additional `padding` is used.  This additional data, before and
 * after the requested interval, will be discarded before output (or before
 * computing the power spectrum) which will remove filter transients.  Note
 * that data will be read for the entire interval
 * [`tstart-padding`,`tstart+deltat+padding`) and must be available.
 * 
 *
 * ### Options
 *
 * <DL>
 * <DT>`-h`, `--help`</DT>
 * <DD>Prints the help message.</DD>
 * <DT>`-c channel`, `--channel=channel`</DT>
 * <DD>The channel name that is to be read.</DD>
 * <DT>`-f cachefile`, `--frame-cache=cachefile`</DT>
 * <DD>The cachefile indexing the frame files to be used.</DD>
 * <DT>`-g globstring `, `--frame-glob=globstring`</DT>
 * <DD>The globstring identifying the frame files to be used.</DD>
 * <DT>`-o outfile`, `--output=outfile`</DT>
 * <DD>The  `outfile` to use.  The extension of `outfile` is used to determine
 * the format of the output.  The output formats are described below. </DD>
 * <DT>`-s tstart`, `--start-time=tstart`</DT>
 * <DD>The time `tstart` GPS seconds of the data to read.  If padding is
 * specified with the `-P` option or the `--pad` option, an additional amount
 * of data preceding `tstart` is also read, and discarded after any requested
 * filtering of the data is complete.</DD>
 * <DT>`-t deltat`, `--duration=deltat`</DT>
 * <DD>The duration `deltat` in seconds of data to read.  If padding is
 * specified with the `-P` option or the `--pad` option, an additional amount
 * of data is also read, and discarded after  any  requested filtering of the
 * data is complete.</DD>
 * <DT>`-H minfreq`, `--highpass=minfreq`</DT>
 * <DD>High-pass filter the data at `minfreq` Hertz.  An additional amount of
 * data, which will be discarded after the filtering, should be specified with
 * the `-P` option or the `--pad` option to remove filter transients.</DD>
 * <DT>`-L maxfreq`, `--lowpass=maxfreq`</DT>
 * <DD>Low-pass filter the data at `maxfreq` Hertz.  An additional amount of
 * data, which will be discarded after the filtering, should be specified with
 * the `-P` option or the `--pad` option to remove filter transients.</DD>
 * <DT>`-P padding`, `--pad=padding`</DT>
 * <DD>Read padding additional seconds of data before and after the requested
 * interval.  This data is then dropped after all requested filtering has been
 * performed in order to remove filter transients.</DD> 
 * <DT>`-R srate`, `--resample=srate`</DT>
 * <DD>Resample the data to sampling rate `srate` Hertz.  An additional amount
 * of data, which will be discarded after the filtering, should be specified
 * with the `-P` option or the `--pad` option to remove filter transients.</DD> 
 * <DT>`-S resolution`, `--spectrum=resolution`</DT>
 * <DD>Compute the power spectrum of the data at the requested `resolution` in
 * Hertz.  Depending on the output format, either the amplitude spectral
 * density or the power spectral density is output.</DD>
 * </DL>
 *
 * ### Output Formats
 *
 * If the -o or --output option is used to specify an output file, the
 * extension  of that file is used to determine the output format.  Supported
 * output formats include:
 * <DL>
 * <DT>`.au`</DT>
 * <DD>Au audio file format.  This format is for time-series data only and may
 * not be used if the `-S` or the `--spectrum` option is used.</DD>
 * <DT>`.eps`</DT>
 * <DD>Encapsulated PostScript (EPS) graphical file format.  A plot of the data
 * values as a function of time is produced using the `gnuplot(1)` program, if
 * available.  If the `-S` or the `--spectrum` option is used, a log-log plot
 * of the amplitude spectrum is produced instead.</DD>
 * <DT>`.gif`</DT>
 * <DD>Graphics Interchange Format (GIF) graphical file format.  A plot of the
 * data values as a function of time is produced using the `gnuplot(1)`
 * program, if available.  If the `-S` or the `--spectrum` option is used, a
 * log-log plot of the amplitude spectrum is produced instead.</DD>
 * <DT>`.jpg`</DT>
 * <DD>JPEG graphical file format.  A plot of the data values as a function of
 * time is produced using the `gnuplot(1)` program, if available.  If the `-S`
 * or the `--spectrum` option is used, a log-log plot of the amplitude spectrum
 * is produced instead.</DD>
 * <DT>`.pdf`</DT>
 * <DD>Portable Document Format (PDF) graphical file format.  A plot of the
 * data values as a function of time is produced using the `gnuplot(1)`
 * program, if available.  If the `-S` or the `--spectrum` option is used, a
 * log-log plot of the amplitude spectrum is produced instead.</DD>
 * <DT>`.png`</DT>
 * <DD>Portable Network Graphics (PNG) graphical file format.  A plot of the
 * data values as a function of time is produced using the `gnuplot(1)`
 * program, if available.  If the `-S` or the `--spectrum` option is used, a
 * log-log plot of the amplitude spectrum is produced instead.</DD>
 * <DT>`.ps`</DT>
 * <DD>PostScript (PS) graphical file format.  A plot of the data values as a
 * function of time is produced using the `gnuplot(1)` program, if available.
 * If the `-S` or the `--spectrum` option is used, a log-log plot of the
 * amplitude spectrum is produced instead.</DD>
 * <DT>`.svg`</DT>
 * <DD>Scalable Vector Graphics (SVG) graphical file format.  A plot of the
 * data values as a function of time is produced using the `gnuplot(1)`
 * program, if available.  If the `-S` or the `--spectrum` optione used, a
 * log-log plot of the amplitude spectrum is produced instead.</DD>
 * <DT>`.wav`</DT>
 * <DD>Waveform Audio File Format (WAVE) audio file format.  This format is for
 * time-series data only and may not be used if the `-S` or the `--spectrum`
 * option is used.</DD>
 * <DT>`.xml`</DT>
 * <DD>XML-based LIGO-lightweight (LIGOLw) file format.  If the `-S` or the
 * `--spectrum` option is used, the power spectral density data is written to
 * the file.</DD>
 * </DL>
 *
 * If none of these extensions are used then the output will be in two-column
 * ascii format.  The first column will be the GPS time of each sample of data
 * and the second column will be the sample values.  However, if the `-S` or
 * the `--spectrum` option is used, the first column will be the frequency of
 * each sample of the spectrum and the second column will be the value of the
 * power spectral density at that frequency.
 * 
 *
 * ### Environment
 * 
 * The `LAL_DEBUG_LEVEL` can used to control the error and warning reporting of
 * `lalfr-vis`.  Common values are: `LAL_DEBUG_LEVEL=0` which suppresses error
 * messages, `LAL_DEBUG_LEVEL=1`  which prints error messages alone,
 * `LAL_DEBUG_LEVEL=3` which prints both error messages and warning messages,
 * and `LAL_DEBUG_LEVEL=7` which additionally prints informational messages.
 *
 *
 * ### Exit Status
 *
 * The `lalfr-vis` utility exits 0 on success, and >0 if an error occurs.
 *
 * ### Examples
 * 
 * The command:
 * 
 *     lalfr-vis -c H1:LSC-STRAIN -g "H-*.gwf" -s 1000000000 -t 16 -o out.wav
 *
 * will read 16 seconds beginning at GPS time 1000000000 of `H1:LSC-STRAIN`
 * data from frame files matching `H-*.gwf` in the current directory and output
 * the data as a WAVE audio file `out.wav`.
 *
 * The command:
 *
 *     lalfr-vis -c L1:LSC-STRAIN -f LLO.cache -s 1000000001 -t 64 -R 2048 -H 10 -L 1000 -P 1 -S 0.25 -o out.png
 *
 * will read 66 seconds beginning at GPS time 1000000000 of `L1:LSC-STRAIN`
 * data from frame files indexed in `LLO.cache`, and the following
 * manipulations will be performed: the data will be resampled to a sampling
 * rate of 2048 Hz, the data will be high-pass filtered at 10 Hz, the data will
 * be low-pass filtered at 1000 Hz, the first and last 1 second of data will be
 * dropped (to remove filter transients), a power spectrum will be computed
 * with 0.25 Hz resolution, and a PNG file displaying a log-log plot of the
 * amplitude spectral density will output in file `out.pnd`.
 *
 * @sa lalfr_stream
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <lal/LALStdlib.h>
#include <lal/LALgetopt.h>
#include <lal/LALString.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/Window.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Audio.h>
#include <lal/LALFrStream.h>
#include <lal/Units.h>

#define FS "\t" /* field separator */
#define RS "\n" /* record separator */
#define MAXDUR 16.0     /* max segment duration (s) */
#define MAXLEN 16384    /* max number of points in buffer */

/* globals */
LALCache *cache;
char *channel;
char *outfile;
double minfreq;
double maxfreq;
double srate;
double pad;
double df;
double t0;
double dt;

int usage(const char *program);
int parseargs(int argc, char **argv);
void output_ts(const char *fname, REAL8TimeSeries *series);
void gnuplot_output_ts(const char *fname, const char *fmt, REAL8TimeSeries *series);
void output_fs(const char *fname, REAL8FrequencySeries *series);
void gnuplot_output_fs(const char *fname, const char *fmt, REAL8FrequencySeries *series);

void xml_output_ts(const char *fname, REAL8TimeSeries *series);
void xml_output_fs(const char *fname, REAL8FrequencySeries *series);

int main(int argc, char *argv[])
{
    LALFrStream *stream;
    REAL8TimeSeries *series;
    LIGOTimeGPS start;

    XLALSetErrorHandler(XLALAbortErrorHandler);

    parseargs(argc, argv);

    /* get the data */
    stream = XLALFrStreamCacheOpen(cache);
    XLALGPSSetREAL8(&start, t0 - pad);
    series = XLALFrStreamInputREAL8TimeSeries(stream, channel, &start, dt + 2.0 * pad, 0);
    XLALFrStreamClose(stream);

    /* manipulate the data */
    if (srate > 0)
        XLALResampleREAL8TimeSeries(series, 1.0 / srate);
    if (minfreq > 0)
        XLALHighPassREAL8TimeSeries(series, minfreq, 0.9, 8);
    if (maxfreq > 0)
        XLALLowPassREAL8TimeSeries(series, maxfreq, 0.9, 8);
    if (pad > 0)
        series = XLALResizeREAL8TimeSeries(series, pad / series->deltaT, dt / series->deltaT);

    if (df > 0) { /* we are computing a spectrum */
        REAL8FrequencySeries *spectrum;
        REAL8FFTPlan *plan;
        REAL8Window *window;
        size_t seglen = 1.0 / (df * series->deltaT);

        /* make sure that the time series length is commensurate with seglen */
        if (((2 * series->data->length) % seglen) != 0) {
            size_t newlen = ((2 * series->data->length) / seglen) * seglen;
            series = XLALResizeREAL8TimeSeries(series, 0, newlen);
        }

        spectrum = XLALCreateREAL8FrequencySeries(series->name, &series->epoch, 0.0, df, &lalDimensionlessUnit, seglen/2 + 1);
        plan = XLALCreateForwardREAL8FFTPlan(seglen, 0);
        window = XLALCreateHannREAL8Window(seglen);
        XLALREAL8AverageSpectrumWelch(spectrum, series, seglen, seglen/2, window, plan);
        if (minfreq > 0 || maxfreq > 0) {
            size_t first = minfreq / spectrum->deltaF;
            size_t last = maxfreq > 0 ? maxfreq / spectrum->deltaF : spectrum->data->length;
            spectrum = XLALResizeREAL8FrequencySeries(spectrum, first, last - first);
        }
        output_fs(outfile, spectrum);
        XLALDestroyREAL8Window(window);
        XLALDestroyREAL8FFTPlan(plan);
        XLALDestroyREAL8FrequencySeries(spectrum);
    } else { /* we are outputting a time series */
        output_ts(outfile, series);
    }

    XLALDestroyREAL8TimeSeries(series);
    return 0;
}

int parseargs(int argc, char **argv)
{
    struct LALoption long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"channel", required_argument, 0, 'c'},
        {"frame-cache", required_argument, 0, 'f'},
        {"frame-glob", required_argument, 0, 'g'},
        {"output", required_argument, 0, 'o'},
        {"start-time", required_argument, 0, 's'},
        {"duration", required_argument, 0, 't'},
        {"highpass", required_argument, 0, 'H'},
        {"lowpass", required_argument, 0, 'L'},
        {"pad", required_argument, 0, 'P'},
        {"resample", required_argument, 0, 'R'},
        {"spectrum", required_argument, 0, 'S'},
        {0, 0, 0, 0}
    };
    char args[] = "hc:f:g:o:s:t:H:L:P:R:S:";
    while (1) {
        int option_index = 0;
        int c;

        c = LALgetopt_long_only(argc, argv, args, long_options, &option_index);
        if (c == -1)    /* end of options */
            break;

        switch (c) {
        case 0:        /* if option set a flag, nothing else to do */
            if (long_options[option_index].flag)
                break;
            else {
                fprintf(stderr, "error parsing option %s with argument %s\n",
                    long_options[option_index].name, LALoptarg);
                exit(1);
            }
        case 'h':      /* help */
            usage(argv[0]);
            exit(0);
        case 'c':      /* channel */
            channel = XLALStringDuplicate(LALoptarg);
            break;
        case 'f':      /* frame-cache */
            cache = XLALCacheImport(LALoptarg);
            break;
        case 'g':      /* frame-glob */
            cache = XLALCacheGlob(NULL, LALoptarg);
            break;
        case 'o':      /* output */
            outfile = XLALStringDuplicate(LALoptarg);
            break;
        case 's':      /* start-time */
            t0 = atof(LALoptarg);
            break;
        case 't':      /* duration */
            dt = atof(LALoptarg);
            break;
        case 'H':      /* highpass */
            minfreq = atof(LALoptarg);
            break;
        case 'L':      /* lowpass */
            maxfreq = atof(LALoptarg);
            break;
        case 'P':      /* pad */
            pad = atof(LALoptarg);
            break;
        case 'R':      /* start-time */
            srate = atof(LALoptarg);
            break;
        case 'S':      /* spectrum */
            df = atof(LALoptarg);
            break;
        case '?':
        default:
            fprintf(stderr, "unknown error while parsing options\n");
            exit(1);
        }
    }

    if (LALoptind < argc) {
        fprintf(stderr, "extraneous command line arguments:\n");
        while (LALoptind < argc)
            fprintf(stderr, "%s\n", argv[LALoptind++]);
        exit(1);
    }

    /* sanity check parameters */

    if (!channel) {
        fprintf(stderr, "must specify a channel\n");
        usage(argv[0]);
        exit(1);
    }
    if (!cache) {
        fprintf(stderr, "must specify a frame cache or frame files\n");
        usage(argv[0]);
        exit(1);
    }

    return 0;
}

int usage(const char *program)
{
    /* basename */
    program = strrchr(program, '/') ? strrchr(program, '/') + 1 : program;
    fprintf(stderr, "usage: %s [options]\n", program);
    fprintf(stderr, "options:\n");
    fprintf(stderr, "\t-h, --help                   \tprint this message and exit\n");
    fprintf(stderr, "\t-c CHAN, --channel=CHAN      \tchannel name CHAN\n");
    fprintf(stderr, "\t-f CACHE, --frame-cache=CACHE\tframe cache file CACHE\n");
    fprintf(stderr, "\t-g GLOB, --frame-glob=GLOB   \tframe file glob string GLOB\n");
    fprintf(stderr, "\t-o OUTFILE, --output=OUTFILE \toutput to file OUTFILE\n");
    fprintf(stderr, "\t-s T0, --start-time=T0       \tGPS start time T0 (s)\n");
    fprintf(stderr, "\t-t DT, --duration=DT         \tduration DT (s)\n");
    fprintf(stderr, "\t-H FMIN, --highpass=FMIN     \thighpass filter at frequency FMIN (Hz)\n");
    fprintf(stderr, "\t-L FMAX, --lowpass=FMAX      \tlowpass filter at frequency FMAX (Hz)\n");
    fprintf(stderr, "\t-P PAD, --pad=PAD            \tadd PAD data at start and end (s)\n");
    fprintf(stderr, "\t-R SRATE, --resample=SRATE   \tresample to rate SRATE (Hz)\n");
    fprintf(stderr, "\t-S DF, --spectrum=DF         \tcompute spectrum with resolution DF (Hz)\n");
    fprintf(stderr, "\nOUTPUT FORMATS:\n");
    fprintf(stderr, "\t");
    fprintf(stderr, ".au, ");
    fprintf(stderr, ".eps, ");
    fprintf(stderr, ".gif, ");
    fprintf(stderr, ".jpg, ");
    fprintf(stderr, ".pdf, ");
    fprintf(stderr, ".png, ");
    fprintf(stderr, ".ps, ");
    fprintf(stderr, ".svg, ");
    fprintf(stderr, ".wav, ");
    fprintf(stderr, ".xml\n");
    fprintf(stderr, "\nEXAMPLES:\n");
    fprintf(stderr, "\n\tOutput 64s of strain data to a WAV file:\n");
    fprintf(stderr, "\t\t%s -c \"H1:LSC-STRAIN\" -g H-H1_RDS_C03_L2-864903135-128.gwf -s 864903136 -t 64 -P 1 -H 10 -L 1000 -R 2048 -o data.wav\n", program);
    fprintf(stderr, "\n\tOutput amplitude spectrum to a PDF file:\n");
    fprintf(stderr, "\t\t%s -c \"H1:LSC-STRAIN\" -g H-H1_RDS_C03_L2-864903135-128.gwf -s 864903136 -t 64 -P 1 -H 10 -L 1000 -R 2048 -S 0.25 -o spec.pdf\n", program);
    return 0;
}

void output_ts(const char *fname, REAL8TimeSeries *series)
{
    const char *ext;
    FILE *fp;
    /* determine output type from extension */
    ext = strrchr(fname ? fname : "", '.');
    ext = ext ? ext + 1 : "";
    if (XLALStringCaseCompare(ext, "wav") == 0) {
        fp = fopen(fname, "w");
        XLALAudioWAVRecordREAL8TimeSeries(fp, series);
        fclose(fp);
    }
    else if (XLALStringCaseCompare(ext, "au") == 0) {
        fp = fopen(fname, "w");
        XLALAudioAURecordREAL8TimeSeries(fp, series);
        fclose(fp);
    }
    else if (
        XLALStringCaseCompare(ext, "eps") == 0
        || XLALStringCaseCompare(ext, "jpeg") == 0
        || XLALStringCaseCompare(ext, "jpg") == 0
        || XLALStringCaseCompare(ext, "gif") == 0
        || XLALStringCaseCompare(ext, "png") == 0
        || XLALStringCaseCompare(ext, "ps") == 0
        || XLALStringCaseCompare(ext, "pdf") == 0
        || XLALStringCaseCompare(ext, "svg") == 0
        )
        gnuplot_output_ts(fname, ext, series);
    else if (XLALStringCaseCompare(ext, "xml") == 0)
        xml_output_ts(fname, series);
    else { /* default is an ascii file */
        size_t i;
        fp = fname ? fopen(fname, "w") : stdout;
        for (i = 0; i < series->data->length; ++i) {
            char tstr[32];
            LIGOTimeGPS t = series->epoch;
            XLALGPSToStr(tstr, XLALGPSAdd(&t, i * series->deltaT));
            fprintf(fp, "%s\t%.15g\n", tstr, series->data->data[i]);
        }
        if (fname)
            fclose(fp);
    }
    return;
}

void output_fs(const char *fname, REAL8FrequencySeries *series)
{
    const char *ext;
    FILE *fp;
    /* determine output type from extension */
    ext = strrchr(fname ? fname : "", '.');
    ext = ext ? ext + 1 : "";
    if (XLALStringCaseCompare(ext, "wav") == 0) {
        fprintf(stderr, "cannot output a spectrum to an audio file\n");
        exit(1);
    }
    else if (XLALStringCaseCompare(ext, "au") == 0) {
        fprintf(stderr, "cannot output a spectrum to an audio file\n");
        exit(1);
    }
    else if (
        XLALStringCaseCompare(ext, "eps") == 0
        || XLALStringCaseCompare(ext, "jpeg") == 0
        || XLALStringCaseCompare(ext, "jpg") == 0
        || XLALStringCaseCompare(ext, "gif") == 0
        || XLALStringCaseCompare(ext, "png") == 0
        || XLALStringCaseCompare(ext, "ps") == 0
        || XLALStringCaseCompare(ext, "pdf") == 0
        || XLALStringCaseCompare(ext, "svg") == 0
        )
        gnuplot_output_fs(fname, ext, series);
    else if (XLALStringCaseCompare(ext, "xml") == 0)
        xml_output_fs(fname, series);
    else { /* default is an ascii file */
        size_t i;
        fp = fname ? fopen(fname, "w") : stdout;
        for (i = 0; i < series->data->length; ++i)
            fprintf(fp, "%f\t%.15g\n", i * series->deltaF, series->data->data[i]);
        if (fname)
            fclose(fp);
    }
    return;
}

void gnuplot_output_ts(const char *fname, const char *fmt, REAL8TimeSeries *series)
{
    size_t i;
    char *units;
    FILE *gp = popen("gnuplot -persistent", "w");

    if (!gp) {
        fprintf(stderr, "require program gnuplot to output .%s files", fmt);
        exit(1);
    }

    /* get the sample units as a string */
    units = XLALUnitToString(&series->sampleUnits);
    if (units == NULL || *units == '\0') {
        LALFree(units);
        units = XLALStringDuplicate("unknown units");
    }

    /* modify fmt as required by gnuplot */
    if (strcmp(fmt, "jpg") == 0)
        fmt = "jpeg";
    else if (strcmp(fmt, "ps") == 0)
        fmt = "postscript landscape";
    else if (strcmp(fmt, "eps") == 0)
        fmt = "postscript eps";

    /* issue gnuplot commands */
    fprintf(gp, "set terminal %s\n", fmt);
    fprintf(gp, "set output '%s'\n", fname);
    fprintf(gp, "set key off\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set xlabel 'time (s)'\n");
    fprintf(gp, "set ylabel 'value (%s)'\n", units);
    fprintf(gp, "set title '%s @ %d.%09d\n", series->name, series->epoch.gpsSeconds, series->epoch.gpsNanoSeconds);
    fprintf(gp, "plot '-' with lines\n");
    for (i = 0; i < series->data->length; ++i)
        fprintf(gp, "%.9f\t%e\n", i * series->deltaT, series->data->data[i]);
    fprintf(gp, "e");

    pclose(gp);
    LALFree(units);
    return;
}

void gnuplot_output_fs(const char *fname, const char *fmt, REAL8FrequencySeries *series)
{
    LALUnit asdunit;
    char *units;
    size_t i;
    FILE *gp = popen("gnuplot -persistent", "w");

    if (!gp) {
        fprintf(stderr, "require program gnuplot to output .%s files", fmt);
        exit(1);
    }

    /* get the sample units as a string */
    XLALUnitSqrt(&asdunit, &series->sampleUnits);
    units = XLALUnitToString(&asdunit);
    if (units == NULL || *units == '\0') {
        LALFree(units);
        units = XLALStringDuplicate("unknown units");
    }

    /* modify fmt as required by gnuplot */
    if (strcmp(fmt, "jpg") == 0)
        fmt = "jpeg";
    else if (strcmp(fmt, "ps") == 0)
        fmt = "postscript landscape";
    else if (strcmp(fmt, "eps") == 0)
        fmt = "postscript eps";

    /* issue gnuplot commands */
    fprintf(gp, "set terminal %s\n", fmt);
    fprintf(gp, "set output '%s'\n", fname);
    fprintf(gp, "set key off\n");
    fprintf(gp, "set grid xtics mxtics ytics\n");
    fprintf(gp, "set xlabel 'frequency (Hz)'\n");
    fprintf(gp, "set ylabel 'amplitude spectral density (%s)'\n", units);
    fprintf(gp, "set title '%s @ %d.%09d\n", series->name, series->epoch.gpsSeconds, series->epoch.gpsNanoSeconds);
    fprintf(gp, "set logscale\n");
    fprintf(gp, "plot '-' with lines\n");
    for (i = 0; i < series->data->length; ++i)
        fprintf(gp, "%.9f\t%e\n", series->f0 + i * series->deltaF, sqrt(series->data->data[i]));
    fprintf(gp, "e");

    pclose(gp);
    LALFree(units);
    return;
}



void xml_begin_xml(FILE *fp);
void xml_end_xml(FILE *fp);

void xml_begin_freq_series(FILE *fp);
void xml_begin_time_series(FILE *fp);
void xml_end_series(FILE *fp);

void xml_put_gps(LIGOTimeGPS *epoch, FILE *fp);
void xml_put_f0(double f0, FILE *fp);
void xml_put_det(const char *name, FILE *fp);

void xml_begin_freq_array(REAL8FrequencySeries *series, FILE *fp);
void xml_begin_time_array(REAL8TimeSeries *series, FILE *fp);
void xml_end_array(FILE *fp);

void xml_put_stream(REAL8Vector *vector, double dx, FILE *fp);

void xml_begin_xml(FILE *fp)
{
    fputs("<?xml version='1.0' encoding='utf-8'?>\n", fp);
    fputs("<!DOCTYPE LIGO_LW SYSTEM \"http://ldas-sw.ligo.caltech.edu/doc/ligolwAPI/html/ligolw_dtd.txt\">\n", fp);
    fputs("<LIGO_LW>\n", fp);
    return;
}
void xml_end_xml(FILE *fp)
{
    fputs("</LIGO_LW>\n", fp);
    return;
}

void xml_begin_freq_series(FILE *fp)
{
    fputs("\t<LIGO_LW Name=\"REAL8FrequencySeries\">\n", fp);
    return;
}
void xml_begin_time_series(FILE *fp)
{
    fputs("\t<LIGO_LW Name=\"REAL8TimeSeries\">\n", fp);
    return;
}
void xml_end_series(FILE *fp)
{
    fputs("\t</LIGO_LW>\n", fp);
    return;
}

void xml_put_gps(LIGOTimeGPS *epoch, FILE *fp)
{
    char tstr[32];
    XLALGPSToStr(tstr, epoch);
    fprintf(fp, "\t\t<Time Type=\"GPS\" Name=\"epoch\">%s</Time>\n", tstr);
    return;
}
void xml_put_f0(double f0, FILE *fp)
{
    fprintf(fp, "\t\t<Param Type=\"real_8\" Name=\"f0:param\" Unit=\"s^-1\">"
        "%g</Param>\n", f0);
    return;
}
void xml_put_det(const char *name, FILE *fp)
{
    const char *dets[] = {"G1", "H1", "H2", "K1", "L1", "T1", "V1"};
    size_t ndet = XLAL_NUM_ELEM(dets);
    size_t d;
    fputs("\t\t<Param Type=\"lstring\" Name=\"instrument:param\">", fp);
    for (d = 0; d < ndet; ++d)
        if (strncmp(name, dets[d], 2) == 0) {
            fputs(dets[d], fp);
            break;
        }
    if (d == ndet)
        fputs("??", fp);
    fputs("</Param>\n", fp);
    return;
}

void xml_begin_freq_array(REAL8FrequencySeries *series, FILE *fp)
{
    fprintf(fp, "\t\t<Array Type=\"real_8\" Name=\"%s:array\" Unit=\"%s\">\n",
        series->name, XLALUnitToString(&series->sampleUnits));
    fprintf(fp, "\t\t\t<Dim Start=\"%g\" Scale=\"%.15g\" "
        "Name=\"Frequency\" Unit=\"s^-1\">%u</Dim>\n", series->f0,
        series->deltaF, series->data->length);
    fputs("\t\t\t<Dim Name=\"Frequency,Real\">2</Dim>\n", fp);
    return;
}
void xml_begin_time_array(REAL8TimeSeries *series, FILE *fp)
{
    fprintf(fp, "\t\t<Array Type=\"real_8\" Name=\"%s:array\" Unit=\"%s\">\n",
        series->name, XLALUnitToString(&series->sampleUnits));
    fprintf(fp, "\t\t\t<Dim Start=\"0\" Scale=\"%.15g\" "
        "Name=\"Time\" Unit=\"s\">%u</Dim>\n", series->deltaT,
        series->data->length);
    fputs("\t\t\t<Dim Name=\"Time,Real\">2</Dim>\n", fp);
    return;
}
void xml_end_array(FILE *fp)
{
    fputs("\t\t</Array>\n", fp);
    return;
}

void xml_put_stream(REAL8Vector *vector, double dx, FILE *fp)
{
#define DELIM ","
    size_t i;
    fputs("\t\t\t<Stream Delimiter=\"" DELIM "\" Type=\"Local\">\n", fp);
    for (i = 0; i < vector->length; ++i)
        fprintf(fp, "\t\t\t\t%.15g" DELIM "%.15g" DELIM "\n", i * dx,
            vector->data[i]);
    fputs("\t\t\t</Stream>\n", fp);
    return;
#undef DELIM
}

void xml_output_ts(const char *fname, REAL8TimeSeries *series)
{
    FILE *fp = fopen(fname, "w");
    xml_begin_xml(fp);

    xml_begin_time_series(fp);

    xml_put_gps(&series->epoch, fp);
    xml_put_f0(series->f0, fp);

    xml_begin_time_array(series, fp);

    xml_put_stream(series->data, series->deltaT, fp);

    xml_end_array(fp);

    xml_put_det(series->name, fp);

    xml_end_series(fp);

    xml_end_xml(fp);

    fclose(fp);
    return;
}

void xml_output_fs(const char *fname, REAL8FrequencySeries *series)
{
    FILE *fp = fopen(fname, "w");
    xml_begin_xml(fp);

    xml_begin_freq_series(fp);

    xml_put_gps(&series->epoch, fp);
    xml_put_f0(series->f0, fp);

    xml_begin_freq_array(series, fp);

    xml_put_stream(series->data, series->deltaF, fp);

    xml_end_array(fp);

    xml_put_det(series->name, fp);

    xml_end_series(fp);

    xml_end_xml(fp);

    fclose(fp);
    return;
}
