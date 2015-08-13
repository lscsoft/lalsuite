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

/*
 * The following command injects waveform.dat directly overhead in Hanford:
 *
 * lalsim-detector-strain --detector-prefix H1 --gps-time 1079467252.35 --ra 0 --dec 46.45515 --psi 324 --verbose waveform.dat
 * and this injects it directly below:
 *
 * lalsim-detector-strain --detector-prefix H1 --gps-time 1079467252.35 --ra 180 --dec -46.45515 --psi -324 --verbose waveform.dat
 */

/**
 * @defgroup lalsim_detector_strain lalsim-detector-strain
 * @ingroup lalsimulation_programs
 *
 * @brief Computes the strain on a detector given a gravitational waveform
 *
 * ### Synopsis
 *
 *     lalsim-detector-strain [options] [file]
 *
 * ### Description
 *
 * The `lalsim-detector-strain` utility converts a gravitational waveform
 * in `file` or standard input if `file` is absent into an induced detector
 * strain on a specified detector.  The input data should be in a
 * three-column ascii format with the first column being the time of each
 * sample, the second column being the plus-polarization of the gravitational
 * waveform, and the third column being the cross-polarization of the
 * gravitational waveform.  The timestamps on the input file should be
 * centered time 0 (conventions vary: some waveforms will end near time 0;
 * others will be centered on time 0).  The output is written to standard
 * output in two column ascii format where the first column is the GPS
 * timestamp of each sample and the second column is the strain induced on
 * the detector.
 *
 * ### Options
 *
 * <DL>
 * <DT>`-h`, `--help`</DT>
 * <DD>print a help message and exit</DD>
 * <DT>`-v`, `--verbose`</DT>
 * <DD>verbose output</DD>
 * <DT>`-r`, `--radians`</DT>
 * <DD>use radians rather than decimal degrees</DD>
 * <DT>`-O`, `--overhead`</DT>
 * <DD>signal from directly overhead</DD>
 * <DT>`-D` PREFIX, `--detector-prefix=`PREFIX</DT>
 * <DD>(required unless overhead) detector prefix (e.g., 'H1', 'L1', 'V1')</DD>
 * <DT>`-t` EPOCH, `--gps-time=`EPOCH</DT>
 * <DD>(required) time of arrival at earth geocenter (or at detector if
 * overhead): this is added to the timestamp of the input data, which should be
 * an waveform about time = 0</DD>
 * <DT>`-a` RA, `--right-ascension=`RA</DT>
 * <DD>(required unless overhead) right ascension in H:M:S format or decimal
 * degrees</DD>
 * <DT>`-d` DEC, `--declination=`DEC</DT>
 * <DD>(required unless overhead) declination in D:M:S format or decimal
 * degrees</DD>
 * <DT>`-p` PSI,` --polarization-angle=`PSI</DT>
 * <DD>(required) polarization angle in degrees</DD>
 * </DL>
 *
 * ### Environment
 *
 * The `LAL_DEBUG_LEVEL` can used to control the error and warning reporting of
 * `lalsim-detector-strain`.  Common values are: `LAL_DEBUG_LEVEL=0` which
 * suppresses error messages, `LAL_DEBUG_LEVEL=1`  which prints error messages
 * alone, `LAL_DEBUG_LEVEL=3` which prints both error messages and warning
 * messages, and `LAL_DEBUG_LEVEL=7` which additionally prints informational
 * messages.
 *
 * ### Exit Status
 *
 * The `lalsim-detector-strain` utility exits 0 on success, and >0 if an error
 * occurs.
 *
 * ### Example
 *
 * The command:
 *
 *     lalsim-inspiral | lalsim-detector-strain -D H1 -a 1:23:45 -d 45.0 -p 30.0 -t 1000000000
 *
 * outputs to standard output in two-column ascii format the strain induced
 * on the LHO observatory detector from a 1.4 solar mass + 1.4 solar mass
 * binary inspiral at 1 Mpc distance originating from source at
 * right-ascension 1h 23m 45s, declination 45 degrees, and polarization angle
 * 30 degrees that arrives at the geocenter at GPS time 1000000000.  The
 * first column contains the GPS time of each sample and the second column
 * contains the induced detector strain at that time.
 *
 * The command:
 *
 *     lalsim-inspiral | lalsim-detector-strain -O -p 0.0 -t 1000000000
 *
 * produces a similar output, but now for a signal coming from directly
 * overhead of a arbitrary detector with polarization angle 0 (optimally
 * oriented); the GPS arrival time is now at the detector location.
 */


#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <lal/LALConstants.h>
#include <lal/LALgetopt.h>
#include <lal/LALDatatypes.h>
#include <lal/LALDetectors.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/LALSimulation.h>

#define INVALID_EPOCH { .gpsSeconds = LAL_INT4_MAX, .gpsNanoSeconds = LAL_INT4_MAX }

/*
 * Timezones for certain detectors. 
 *
 * PSTPDT: Pacific Time
 * PST (standard): UTC-8
 * PDT (daylight): UTC-7
 * Daylight begin: Second Sunday of March at 02:00 PST
 * Daylight end: Last Sunday of October at 02:00 PDT
 *
 * CSTCDT: Pacific Time
 * PST (standard): UTC-6
 * PDT (daylight): UTC-5
 * Daylight begin: Second Sunday of March at 02:00 PST
 * Daylight end: Last Sunday of October at 02:00 PDT
 *
 * CETCEST: Central Europe Time
 * CET (standard): UTC+1
 * CEST (daylight): UTC+2
 * Daylight begin: Last Sunday of March at 02:00 CET = 01:00 UTC
 * Daylight end: Last Sunday of October at 03:00 CEST = 01:00 UTC
 *
 * JST: Japanese Standard Time
 * JST: UTC+9
 * No daylight savings time
 */
#define TZ_PSTPDT  "PST08:00:00PDT07:00:00,M3.2.0/02:00:00,M11.1.0/02:00:00"
#define TZ_CSTCDT  "CST06:00:00CDT05:00:00,M3.2.0/02:00:00,M11.1.0/02:00:00"
#define TZ_CETCEST "CET-01:00:00CEST-02:00:00,M3.5.0/02:00:00,M10.5.0/03:00:00"
#define TZ_JST     "JST-09:00:00"

#define TZ_LHO   TZ_PSTPDT
#define TZ_LLO   TZ_CSTCDT
#define TZ_GEO   TZ_CETCEST
#define TZ_VIRGO TZ_CETCEST
#define TZ_KAGRA TZ_JST
#define TZ_TAMA  TZ_JST

struct params {
    int verbose;
    FILE *fp;
    const LALDetector *detector;
    int overhead;
    LIGOTimeGPS epoch;
    double ra;
    double dec;
    double psi;
};
int printparams(struct params p);
char *sitetime(char *s, size_t size, time_t * timer, int site);
double strtora(const char *s, int degrees);
double strtodec(const char *s, int degrees);
int fputra(double ra, FILE * fp);
int fputdec(double ra, FILE * fp);
struct params parseargs(int argc, char **argv);
int usage(const char *program);
int readdata(REAL8TimeSeries ** hplus, REAL8TimeSeries ** hcross, FILE * fp);

int main(int argc, char *argv[])
{
    char tstr[32];      // string to hold GPS time -- 31 characters is enough
    REAL8TimeSeries *hplus = NULL;
    REAL8TimeSeries *hcross = NULL;
    REAL8TimeSeries *h;
    struct params p;
    size_t j;

    XLALSetErrorHandler(XLALBacktraceErrorHandler);

    p = parseargs(argc, argv);
    printparams(p);

    readdata(&hplus, &hcross, p.fp);

    /* shift injection to required epoch */
    XLALGPSAddGPS(&hplus->epoch, &p.epoch);
    XLALGPSAddGPS(&hcross->epoch, &p.epoch);

    /* compute detector strain */
    if (p.overhead) {
        double cos2psi = cos(2.0 * p.psi);
        double sin2psi = sin(2.0 * p.psi);
        /* ra and dec are ignored parameters; psi is used to mix plus and
         * cross polarizations; note sign convention has F+ = -1 and Fx = 0
	 * for psi = 0, and F+ = 0, Fx = 1 for psi = 45 degrees:
	 * see Sec. "The Gravitational Wave Tensor in the Earth Fixed Frame",
	 * LIGO-T010110 https://dcc.ligo.org/LIGO-T010110/public */
        h = XLALCreateREAL8TimeSeries("STRAIN", &hplus->epoch, hplus->f0, hplus->deltaT, &hplus->sampleUnits, hplus->data->length);
        for (j = 0; j < h->data->length; ++j)
            h->data->data[j] = -hplus->data->data[j] * cos2psi + hcross->data->data[j] * sin2psi;
        fprintf(stdout, "# time (s)\tSTRAIN (strain)\n");
    } else {
        h = XLALSimDetectorStrainREAL8TimeSeries(hplus, hcross, p.ra, p.dec, p.psi, p.detector);
        fprintf(stdout, "# time (s)\t%s:STRAIN (strain)\n", p.detector->frDetector.prefix);
    }

    for (j = 0; j < h->data->length; ++j) {
        LIGOTimeGPS t = h->epoch;
        fprintf(stdout, "%s\t%e\n", XLALGPSToStr(tstr, XLALGPSAdd(&t, j * h->deltaT)), h->data->data[j]);
    }

    XLALDestroyREAL8TimeSeries(h);
    XLALDestroyREAL8TimeSeries(hcross);
    XLALDestroyREAL8TimeSeries(hplus);
    LALCheckMemoryLeaks();
    return 0;
}

int printparams(struct params p)
{
    if (p.verbose) {
        char tstr[32];  // string to hold GPS time -- 31 characters is enough
        struct tm utc = {.tm_sec = 0 };
        time_t timer;
        char date[64];
        double gmst_rad = XLALGreenwichMeanSiderealTime(&p.epoch);
        if (p.overhead) {
            XLALGPSToUTC(&utc, round(XLALGPSGetREAL8(&p.epoch)));
            timer = round(XLALGPSGetREAL8(&p.epoch)) + XLAL_EPOCH_UNIX_GPS;
            strftime(date, sizeof(date), "%c %Z", &utc);
            fputs("detector arrival time:        GPS = ", stderr);
            fputs(XLALGPSToStr(tstr, &p.epoch), stderr);
            fputs("\n", stderr);
            fputs("coordinated universal time:   ", stderr);
            fprintf(stderr, "%s", date);
            fputs("\n", stderr);
            fputs("Greenwich mean sidereal time: HH:MM:SS = ", stderr);
            fputra(gmst_rad, stderr);
            fprintf(stderr, " (%g deg, %g rad)", fmod(gmst_rad, 2.0 * LAL_PI) / LAL_PI_180, fmod(gmst_rad, 2.0 * LAL_PI));
            fputs("\n", stderr);
            fputs("polarization angle:           ", stderr);
            fprintf(stderr, "%g deg, %g rad", p.psi / LAL_PI_180, p.psi);
            fputs("\n", stderr);
        } else {
            LIGOTimeGPS epoch = p.epoch;
            double longitude = p.detector->frDetector.vertexLongitudeRadians;
            double latitude = p.detector->frDetector.vertexLatitudeRadians;
            double lmst_rad = gmst_rad + longitude;
            double ha = lmst_rad - p.ra;
            double altitude = asin(sin(p.dec) * sin(latitude) + cos(p.dec) * cos(latitude) * cos(ha));
            double azimuth = -atan2(cos(p.dec) * sin(ha), sin(p.dec) * cos(latitude) - cos(p.dec) * sin(latitude) * cos(ha));
            double fplus, fcross;
            double dt = XLALTimeDelayFromEarthCenter(p.detector->location, p.ra, p.dec, &p.epoch);
            XLALGPSAdd(&epoch, dt);
            XLALGPSToUTC(&utc, round(XLALGPSGetREAL8(&p.epoch)));
            timer = round(XLALGPSGetREAL8(&p.epoch)) + XLAL_EPOCH_UNIX_GPS;
            XLALComputeDetAMResponse(&fplus, &fcross, p.detector->response, p.ra, p.dec, p.psi, gmst_rad);
            strftime(date, sizeof(date), "%c %Z", &utc);
            fputs("geocentric arrival time:      GPS = ", stderr);
            fputs(XLALGPSToStr(tstr, &p.epoch), stderr);
            fputs("\n", stderr);
            fputs("coordinated universal time:   ", stderr);
            fprintf(stderr, "%s", date);
            fputs("\n", stderr);
            fputs("Greenwich mean sidereal time: HH:MM:SS = ", stderr);
            fputra(gmst_rad, stderr);
            fprintf(stderr, " (%g deg, %g rad)", fmod(gmst_rad, 2.0 * LAL_PI) / LAL_PI_180, fmod(gmst_rad, 2.0 * LAL_PI));
            fputs("\n", stderr);
            fputs("right ascension:              HH:MM:SS = ", stderr);
            fputra(p.ra, stderr);
            fprintf(stderr, " (%g deg, %g rad)", p.ra / LAL_PI_180, p.ra);
            fputs("\n", stderr);
            fputs("declination:                  DD:MM:SS = ", stderr);
            fputdec(p.dec, stderr);
            fprintf(stderr, " (%g deg, %g rad)", p.dec / LAL_PI_180, p.dec);
            fputs("\n", stderr);
            fputs("polarization angle:           ", stderr);
            fprintf(stderr, "%g deg, %g rad", p.psi / LAL_PI_180, p.psi);
            fputs("\n", stderr);
            fputs("detector:                     ", stderr);
            fprintf(stderr, "%s", p.detector->frDetector.name);
            fputs("\n", stderr);
            fputs("time delay from Earth center: ", stderr);
            fprintf(stderr, "%+f s", dt);
            fputs("\n", stderr);
            fputs("arrival time at detector:     GPS = ", stderr);
            fputs(XLALGPSToStr(tstr, &epoch), stderr);
            fputs("\n", stderr);
            fputs("local time:                   ", stderr);
            fprintf(stderr, "%s", sitetime(tstr, sizeof(tstr), &timer, *p.detector->frDetector.prefix));
            fputs("\n", stderr);
            fputs("latitude North:               DD:MM:SS = ", stderr);
            fputdec(latitude, stderr);
            fprintf(stderr, " (%g deg, %g rad)", latitude / LAL_PI_180, latitude);
            fputs("\n", stderr);
            fputs("longitude East:               HH:MM:SS = ", stderr);
            fputra(longitude, stderr);
            fprintf(stderr, " (%g deg, %g rad)", longitude / LAL_PI_180, longitude);
            fputs("\n", stderr);
            fputs("local mean sidereal time:     HH:MM:SS = ", stderr);
            fputra(lmst_rad, stderr);
            fprintf(stderr, " (%g deg, %g rad)", fmod(lmst_rad, 2.0 * LAL_PI) / LAL_PI_180, fmod(lmst_rad, 2.0 * LAL_PI));
            fputs("\n", stderr);
            fputs("hour angle:                   HH:MM:SS = ", stderr);
            fputra(ha, stderr);
            fprintf(stderr, " (%g deg, %g rad)", fmod(ha, 2.0 * LAL_PI) / LAL_PI_180, fmod(ha, 2.0 * LAL_PI));
            fputs("\n", stderr);
            fputs("altitude:                     ", stderr);
            fprintf(stderr, "%g deg, %g rad", altitude / LAL_PI_180, altitude);
            fputs("\n", stderr);
            fputs("azimuth West of North:        ", stderr);
            fprintf(stderr, "%g deg, %g rad", azimuth / LAL_PI_180, azimuth);
            fputs("\n", stderr);
            fputs("detector antenna response:    ", stderr);
            fprintf(stderr, "F+ = %g, Fx = %g", fplus, fcross);
            fputs("\n", stderr);
        }
    }
    return 0;
}

char *sitetime(char *s, size_t size, time_t * timer, int site)
{
    static char tz_lho[] = TZ_LHO;
    static char tz_llo[] = TZ_LLO;
    static char tz_geo[] = TZ_GEO;
    static char tz_virgo[] = TZ_VIRGO;
    static char tz_kagra[] = TZ_KAGRA;
    static char tz_tama[] = TZ_TAMA;
    char *tz_orig = getenv("TZ");
    char *tz_site = tz_orig;
    switch (site) {
    case 'H':
        tz_site = tz_lho;
        break;
    case 'L':
        tz_site = tz_llo;
        break;
    case 'G':
        tz_site = tz_geo;
        break;
    case 'V':
        tz_site = tz_virgo;
        break;
    case 'K':
        tz_site = tz_kagra;
        break;
    case 'T':
        tz_site = tz_tama;
        break;
    default:
        fprintf(stderr, "warning: unknown time zone - using localtime\n");
        tz_site = tz_orig;
        break;
    }
    setenv("TZ", tz_site, 1);
    strftime(s, size, "%c %Z", localtime(timer));
    if (tz_orig)
        setenv("TZ", tz_orig, 1);
    else
        unsetenv("TZ");
    return s;
}

int readdata(REAL8TimeSeries ** hplus, REAL8TimeSeries ** hcross, FILE * fp)
{
    const size_t block = 1024;
    LIGOTimeGPS start;
    LIGOTimeGPS end;
    double dt;
    double *hp = NULL;
    double *hc = NULL;
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
            hp = realloc(hp, bufsz * sizeof(*hp));
            hc = realloc(hc, bufsz * sizeof(*hc));
        }
        c = sscanf(line, "%s %le %le", n ? t1 : t0, hp + n, hc + n);
        if (c != 3) {
            fprintf(stderr, "error: format error on line %zd: %s\n", l, line);
            exit(1);
        }
        ++n;
    }
    hp = realloc(hp, n * sizeof(*hp));
    hc = realloc(hc, n * sizeof(*hp));

    XLALStrToGPS(&start, t0, NULL);
    XLALStrToGPS(&end, t1, NULL);
    dt = XLALGPSDiff(&end, &start) / (n - 1);
    *hplus = XLALCreateREAL8TimeSeries("h_plus", &start, 0.0, dt, &lalStrainUnit, n);
    *hcross = XLALCreateREAL8TimeSeries("h_cross", &start, 0.0, dt, &lalStrainUnit, n);
    memcpy((*hplus)->data->data, hp, n * sizeof(*hp));
    memcpy((*hcross)->data->data, hc, n * sizeof(*hc));

    free(hp);
    free(hc);
    return 0;
}

double strtora(const char *string, int degrees)
{
    double h = 0, m = 0, s = 0;
    double ra;
    int c;

    /* try H:M:S format first */
    if (strchr(string, ':'))
        c = sscanf(string, "%lf:%lf:%lf", &h, &m, &s);
    else {
        ra = atof(string);
        if (degrees)
            ra *= LAL_PI_180;
        goto done;
    }

    if (c == 0) {
        fprintf(stderr, "error parsing option --right-ascension with argument %s\n", string);
        exit(1);
    }
    m = copysign(m, h) / 60.0;
    s = copysign(s, h) / 3600.0;
    ra = (h + m + s) * 15.0 * LAL_PI_180;

  done:
    ra = fmod(ra, 2.0 * LAL_PI);
    if (ra < 0.0)
        ra += 2.0 * LAL_PI;
    return ra;
}

double strtodec(const char *string, int degrees)
{
    double d = 0, m = 0, s = 0;
    double dec;
    int c;

    /* try D:M:S format first */
    if (strchr(string, ':'))
        c = sscanf(string, "%lf:%lf:%lf", &d, &m, &s);
    else {
        dec = atof(string);
        if (degrees)
            dec *= LAL_PI_180;
        goto done;
    }

    if (c == 0) {
        fprintf(stderr, "error parsing option --declination with argument %s\n", string);
        exit(1);
    }
    m = copysign(m, d) / 60.0;
    s = copysign(s, d) / 3600.0;
    dec = (d + m + s) * LAL_PI_180;
  done:
    dec = fmod(dec + LAL_PI, 2.0 * LAL_PI) - LAL_PI;
    return dec;
}

/* output right ascension in hms format */
int fputra(double ra, FILE * fp)
{
    double h, m, s;
    ra = fmod(ra, 2.0 * LAL_PI);
    if (ra < 0.0)
        ra += 2.0 * LAL_PI;
    ra /= 15.0 * LAL_PI_180;    /* convert to hours */
    m = 60.0 * modf(ra, &h);
    s = 60.0 * modf(m, &s);
    return fprintf(fp, "%02d:%02d:%04.1f", (int)h, (int)fabs(m), fabs(s));
}

/* output right ascension in hms format */
int fputdec(double dec, FILE * fp)
{
    double d, m, s;
    dec = fmod(dec + LAL_PI, 2.0 * LAL_PI) - LAL_PI;
    dec /= LAL_PI_180;  /* convert to degrees */
    m = 60.0 * modf(dec, &d);
    s = 60.0 * modf(m, &s);
    return fprintf(fp, "%c%02d:%02d:%04.1fs", dec < 0 ? '-' : '+', (int)fabs(d), (int)fabs(m), fabs(s));
}

struct params parseargs(int argc, char **argv)
{
    LIGOTimeGPS invalid_epoch = INVALID_EPOCH;
    char *ra_string = NULL;
    char *dec_string = NULL;
    char *psi_string = NULL;
    int degrees = 1;
    struct params p = {
        .verbose = 0,
        .fp = stdin,
        .overhead = 0,
        .detector = NULL,
        .epoch = invalid_epoch,
        .ra = HUGE_VAL,
        .dec = HUGE_VAL,
        .psi = HUGE_VAL,
    };
    struct LALoption long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"verbose", no_argument, 0, 'v'},
        {"radians", no_argument, 0, 'r'},
        {"overhead", no_argument, 0, 'O'},
        {"detector-prefix", required_argument, 0, 'D'},
        {"gps-time", required_argument, 0, 't'},
        {"right-ascension", required_argument, 0, 'a'},
        {"ra", required_argument, 0, 'a'},
        {"alpha", required_argument, 0, 'a'},
        {"declination", required_argument, 0, 'd'},
        {"dec", required_argument, 0, 'd'},
        {"delta", required_argument, 0, 'd'},
        {"polarization-angle", required_argument, 0, 'p'},
        {"psi", required_argument, 0, 'p'},
        {0, 0, 0, 0}
    };
    char args[] = "hvrOD:t:a:d:p:";
    int fail = 0;
    int d;
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
                fprintf(stderr, "error parsing option %s with argument %s\n", long_options[option_index].name, LALoptarg);
                exit(1);
            }
        case 'h':      /* help */
            usage(argv[0]);
            exit(0);
        case 'v':      /* verbose */
            p.verbose = 1;
            break;
        case 'r':      /* radians */
            degrees = 0;
            break;
        case 'O':      /* overhead */
            p.overhead = 1;
            break;
        case 'D':      /* detector-prefix */
            for (d = 0; d < LAL_NUM_DETECTORS; ++d) {
                if (strcmp(LALoptarg, lalCachedDetectors[d].frDetector.prefix) == 0) {
                    p.detector = lalCachedDetectors + d;
                    break;
                }
            }
            break;
        case 't':      /* gps-time */
            XLALStrToGPS(&p.epoch, LALoptarg, NULL);
            break;
        case 'a':      /* right-ascension */
            ra_string = LALoptarg;
            break;
        case 'd':      /* right-ascension */
            dec_string = LALoptarg;
            break;
        case 'p':      /* polarization-angle */
            psi_string = LALoptarg;
            break;
        case '?':
        default:
            fprintf(stderr, "unknown error while parsing options\n");
            exit(1);
        }
    }
    switch (argc - LALoptind) {
    case 0:
        break;
    case 1:    /* the input file name */
        p.fp = fopen(argv[LALoptind], "r");
        if (!p.fp) {
            fprintf(stderr, "error: could not open file %s\n", argv[LALoptind]);
            exit(1);
        }

        break;
    default:
        fprintf(stderr, "extraneous command line arguments:\n");
        while (++LALoptind < argc)
            fprintf(stderr, "%s\n", argv[LALoptind]);
        exit(1);
    }

    /* make sure all required arguments have been specified */
    if (XLALGPSCmp(&p.epoch, &invalid_epoch) == 0) {
        fprintf(stderr, "error: unspecified arrival time\n");
        fail = 1;
    }

    if (psi_string) {
            p.psi = atof(psi_string);
            if (degrees)
                p.psi *= LAL_PI_180;
    } else {
        fprintf(stderr, "error: unspecified polarization angle\n");
        fail = 1;
    }

    if (p.overhead) {
        /* remaining parameters are ignored */
        if (p.detector)
            fprintf(stderr, "warning: ignoring detector for overhead injection\n");
        if (ra_string)
            fprintf(stderr, "warning: ignoring right ascension for overhead injection\n");
        if (dec_string)
            fprintf(stderr, "warning: ignoring declination for overhead injection\n");
        goto done;
    }

    if (p.detector == NULL) {
        fprintf(stderr, "error: unspecified detector prefix\n");
        fprintf(stderr, "recognized detector prefixes are:\n");
        for (d = 0; d < LAL_NUM_DETECTORS; ++d)
            fprintf(stderr, "\t%s: %s\n", lalCachedDetectors[d].frDetector.prefix, lalCachedDetectors[d].frDetector.name);
        fail = 1;
    }

    if (ra_string)
        p.ra = strtora(ra_string, degrees);
    else {
        fprintf(stderr, "error: unspecified right ascension\n");
        fail = 1;
    }
    
    if (dec_string)
        p.dec = strtodec(dec_string, degrees);
    else {
        fprintf(stderr, "error: unspecified declination\n");
        fail = 1;
    }

done:
    if (fail)
        exit(1);

    return p;
}

int usage(const char *program)
{
    /* *INDENT-OFF* */
    fprintf(stderr, "\
usage: %s [options] [file]\n\
options:\n\
        -h, --help      print this message and exit\n\
        -v, --verbose   verbose output\n\
        -r, --radians   use radians rather than decimal degrees\n\
        -O, --overhead  signal from directly overhead\n\
        -D PREFIX, --detector-prefix=PREFIX     (required unless overhead)\n\
                detector prefix (e.g., 'H1', 'L1', 'V1')\n\
        -t EPOCH, --gps-time=EPOCH              (required)\n\
                time of arrival at earth geocenter (or at detector if overhead):\n\
                this is added to the timestamp of the input data,\n\
                which should be an waveform about time = 0\n\
        -a RA, --right-ascension=RA             (required unless overhead)\n\
                right ascension in H:M:S format or decimal degrees\n\
        -d DEC, --declination=DEC               (required unless overhead)\n\
                declination in D:M:S format or decimal degrees\n\
        -p PSI, --polarization-angle=PSI        (required)\n\
                polarization angle in degrees\n\
", program);
    /* *INDENT-ON* */
    return 0;
}
