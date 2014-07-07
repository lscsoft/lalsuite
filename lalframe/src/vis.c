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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <getopt.h>

#include <lal/LALStdlib.h>
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
    struct option long_options[] = {
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

        c = getopt_long_only(argc, argv, args, long_options, &option_index);
        if (c == -1)    /* end of options */
            break;

        switch (c) {
        case 0:        /* if option set a flag, nothing else to do */
            if (long_options[option_index].flag)
                break;
            else {
                fprintf(stderr, "error parsing option %s with argument %s\n",
                    long_options[option_index].name, optarg);
                exit(1);
            }
        case 'h':      /* help */
            usage(argv[0]);
            exit(0);
        case 'c':      /* channel */
            channel = strdup(optarg);
            break;
        case 'f':      /* frame-cache */
            cache = XLALCacheImport(optarg);
            break;
        case 'g':      /* frame-glob */
            cache = XLALCacheGlob(NULL, optarg);
            break;
        case 'o':      /* output */
            outfile = strdup(optarg);
            break;
        case 's':      /* start-time */
            t0 = atof(optarg);
            break;
        case 't':      /* duration */
            dt = atof(optarg);
            break;
        case 'H':      /* highpass */
            minfreq = atof(optarg);
            break;
        case 'L':      /* lowpass */
            maxfreq = atof(optarg);
            break;
        case 'P':      /* pad */
            pad = atof(optarg);
            break;
        case 'R':      /* start-time */
            srate = atof(optarg);
            break;
        case 'S':      /* spectrum */
            df = atof(optarg);
            break;
        case '?':
        default:
            fprintf(stderr, "unknown error while parsing options\n");
            exit(1);
        }
    }

    if (optind < argc) {
        fprintf(stderr, "extraneous command line arguments:\n");
        while (optind < argc)
            fprintf(stderr, "%s\n", argv[optind++]);
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
    if (strcasecmp(ext, "wav") == 0) {
        fp = fopen(fname, "w");
        XLALAudioWAVRecordREAL8TimeSeries(fp, series);
        fclose(fp);
    }
    else if (strcasecmp(ext, "au") == 0) {
        fp = fopen(fname, "w");
        XLALAudioAURecordREAL8TimeSeries(fp, series);
        fclose(fp);
    }
    else if (
        strcasecmp(ext, "eps") == 0
        || strcasecmp(ext, "jpeg") == 0
        || strcasecmp(ext, "jpg") == 0
        || strcasecmp(ext, "gif") == 0
        || strcasecmp(ext, "png") == 0
        || strcasecmp(ext, "ps") == 0
        || strcasecmp(ext, "pdf") == 0
        || strcasecmp(ext, "svg") == 0
        )
        gnuplot_output_ts(fname, ext, series);
    else if (strcasecmp(ext, "xml") == 0)
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
    if (strcasecmp(ext, "wav") == 0) {
        fprintf(stderr, "cannot output a spectrum to an audio file\n");
        exit(1);
    }
    else if (strcasecmp(ext, "au") == 0) {
        fprintf(stderr, "cannot output a spectrum to an audio file\n");
        exit(1);
    }
    else if (
        strcasecmp(ext, "eps") == 0
        || strcasecmp(ext, "jpeg") == 0
        || strcasecmp(ext, "jpg") == 0
        || strcasecmp(ext, "gif") == 0
        || strcasecmp(ext, "png") == 0
        || strcasecmp(ext, "ps") == 0
        || strcasecmp(ext, "pdf") == 0
        || strcasecmp(ext, "svg") == 0
        )
        gnuplot_output_fs(fname, ext, series);
    else if (strcasecmp(ext, "xml") == 0)
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
    size_t ndet = sizeof(dets)/sizeof(*dets);
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
