#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_rng.h>

#include <lal/LALStdlib.h>
#include <lal/LALgetopt.h>
#include <lal/LALConstants.h>
#include <lal/LALDetectors.h>
#include <lal/LALSimulation.h>
#include <lal/TimeSeries.h>
#include <lal/Date.h>
#include <lal/Audio.h>
#include <lal/LALSimBurst.h>

/* global variables for parameters (with default values) */
LALDetector detector;
LIGOTimeGPS epoch;
double hrss = 1.0;
double f_min;
double f_max;
double deltaT = 1.0/16384.0;
double V = 1.0;
char output[FILENAME_MAX];

int usage(const char *program);
int parseargs(int argc, char **argv);
int fprintgps(FILE *fp, LIGOTimeGPS *t);

#define CAT(a,b) a ## b
#define FLT(i) CAT(i,.)
#define HMS2RAD(h,m,s) (LAL_PI * ((h) + ((m) + (s) / 60.) / 60.) / 12.0)
#define DMS2RAD(d,m,s) ((signbit(FLT(d)) ? -LAL_PI : LAL_PI) * (abs(d) + ((m) + (s) / 60.) / 60.) / 180.0)

int main(int argc, char *argv[])
{
	/* 
	 * RA and DEC are centered on Mon R2 IRS 3
	 * See: J. Giannakopoulou et al. 1997 ApJ 487 346 doi:10.1086/304574
	 */
	double ra  = HMS2RAD( 6,  7, 47.8);
	double dec = DMS2RAD(-6, 22, 55);
	double psi = 0.0;
	REAL8TimeSeries *h, *hplus, *hcross;
	gsl_rng *rng;
	FILE *fp = stdout;

	/* XLALSetErrorHandler(XLALBacktraceErrorHandler); */
	XLALSetErrorHandler(XLALAbortErrorHandler);

	parseargs(argc, argv);

	gsl_rng_env_setup();
	rng = gsl_rng_alloc(gsl_rng_default);

	XLALSimUnicorn(&hplus, &hcross, f_min, f_max, V, hrss, deltaT, rng);
	XLALGPSAddGPS(&hplus->epoch, &epoch);
	XLALGPSAddGPS(&hcross->epoch, &epoch);

	h = XLALSimDetectorStrainREAL8TimeSeries(hplus, hcross, ra, dec, psi, &detector);

	if (*output)
		fp = fopen(output, "w");

	/* determine type of output from filename */
	if (*output && 0 == strcmp(strrchr(output, '.'), ".wav"))
		XLALAudioWAVRecordREAL8TimeSeries(fp, h);
	else if (*output && 0 == strcmp(strrchr(output, '.'), ".au"))
		XLALAudioAURecordREAL8TimeSeries(fp, h);
	else { /* ascii file */
		size_t j;
		fprintf(fp, "# time (s)\t%s:STRAIN (strain)\n", detector.frDetector.prefix);
		for (j = 0; j < h->data->length; ++j) {
			LIGOTimeGPS t = h->epoch;
			XLALGPSAdd(&t, j * h->deltaT);
			fprintgps(fp, &t);
			fprintf(fp, "\t%e\n", h->data->data[j]);
		}
	}

	if (*output)
		fclose(fp);

	XLALDestroyREAL8TimeSeries(h);
	XLALDestroyREAL8TimeSeries(hcross);
	XLALDestroyREAL8TimeSeries(hplus);
	gsl_rng_free(rng);
	LALCheckMemoryLeaks();

	return 0;
}

int parseargs(int argc, char **argv)
{
	struct LALoption long_options[] = {
			{ "help", no_argument, 0, 'h' },
			{ "detector-site", required_argument, 0, 'd' },
			{ "min-frequency", required_argument, 0, 'f' },
			{ "max-frequency", required_argument, 0, 'F' },
			{ "hrss", required_argument, 0, 'H' },
			{ "output-file", required_argument, 0, 'o' },
			{ "sample-rate", required_argument, 0, 's' },
			{ "gps-start-time", required_argument, 0, 't' },
			{ "time-freq-volume", required_argument, 0, 'V' },
			{ 0, 0, 0, 0 }
		};
	char args[] = "hd:f:F:H:o:s:t:V:";
	while (1) {
		int option_index = 0;
		int c;
	
		c = LALgetopt_long_only(argc, argv, args, long_options, &option_index);
		if (c == -1) /* end of options */
			break;
	
		switch (c) {
			case 0: /* if option set a flag, nothing else to do */
				if (long_options[option_index].flag)
					break;
				else {
					fprintf(stderr, "error parsing option %s with argument %s\n", long_options[option_index].name, LALoptarg);
					exit(1);
				}
			case 'h': /* help */
				usage(argv[0]);
				exit(0);
			case 'd': /* detector-site */
				switch (*LALoptarg) {
					case 'G':
						detector = lalCachedDetectors[LAL_GEO_600_DETECTOR];
						break;
					case 'H':
						detector = lalCachedDetectors[LAL_LHO_4K_DETECTOR];
						break;
					case 'L':
						detector = lalCachedDetectors[LAL_LLO_4K_DETECTOR];
						break;
					case 'V':
						detector = lalCachedDetectors[LAL_VIRGO_DETECTOR];
						break;
					default:
						fprintf(stderr, "unrecognized detector site %s - must be 'G', 'H', 'L' or 'V'\n", LALoptarg);
						exit(1);
				}
			case 'f': /* min-frequency */
				f_min = atof(LALoptarg);
				break;
			case 'F': /* max-frequency */
				f_max = atof(LALoptarg);
				break;
			case 'H': /* hrss */
				hrss = atof(LALoptarg);
				break;
			case 'o': /* output-file */
				strncpy(output, LALoptarg, sizeof(output) - 1);
				break;
			case 's': /* sample-rate */
				deltaT = 1.0 / atof(LALoptarg);
				break;
			case 't': /* gps-start-time */
				XLALGPSSetREAL8(&epoch, atof(LALoptarg));
				break;
			case 'V': /* time-freq-volume */
				V = atof(LALoptarg);
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

	if (V < LAL_2_PI) {
		fprintf(stderr, "error: time-frequency volume must be at least 2/pi\n");
		usage(argv[0]);
		exit(1);
	}

	if (strlen(detector.frDetector.prefix) == 0) {
		fprintf(stderr, "error: must specify a detector site\n");
		usage(argv[0]);
		exit(1);
	}

	if (f_max == 0.0) /* default value: Nyquist */
		f_max = 0.5 / deltaT;
	else if (f_max < f_min) {
		fprintf(stderr, "error: the maximum frequency must be greater than the minimum frequency\n");
		usage(argv[0]);
		exit(1);
	}
	else if (f_max > 0.5 / deltaT) {
		fprintf(stderr, "error: the maximum frequency must be less than the Nyquist frequency\n");
		usage(argv[0]);
		exit(1);
	}

	return 0;
}
	
int usage( const char *program )
{
	fprintf(stderr, "usage: %s [options]\n", program);
	fprintf(stderr, "options:\n" );
	fprintf(stderr, "\t-h, --help     \tprint this message and exit\n");
	fprintf(stderr, "\t-d detectorSite\t(required) detector site (G|H|L|V)\n");
	fprintf(stderr, "\t-f minFrequency\t(default=0) minimum frequency (Hz)\n");
	fprintf(stderr, "\t-F maxFrequency\t(default=Nyquist) maximum frequency (Hz)\n");
	fprintf(stderr, "\t-H hrss        \t(default=1) root-sum-squared strain hrss (s)\n");
	fprintf(stderr, "\t-o outfile     \t(default=stdout) output filename\n");
	fprintf(stderr, "\t-s sampleRate  \t(default=16384) sample rate (Hz)\n");
	fprintf(stderr, "\t-t GPSStartTime\t(default=0) start time relative to GPS epoch (s)\n");
	fprintf(stderr, "\t-V TimeFreqVol \t(default=1) pixel time-frequency volume\n");
	fprintf(stderr, "environment:\n" );
	fprintf(stderr, "\tGSL_RNG_SEED   \trandom number generator seed\n");
	fprintf(stderr, "\tGSL_RNG_TYPE   \trandom number generator type\n");
	return 0;
}

int fprintgps(FILE *fp, LIGOTimeGPS *t)
{
	int s = t->gpsSeconds;
	int ns = t->gpsNanoSeconds;
	int sgn = signbit(s + ns);
	return fprintf(fp, "%s%d.%09d", sgn ? "-" : "", abs(s), abs(ns));
}
