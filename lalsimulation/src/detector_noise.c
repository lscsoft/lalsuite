/*
*  Copyright (C) 2007 Jolien Creighton
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
#include <getopt.h>

#include <gsl/gsl_rng.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/LALSimNoise.h>

double (*psdfunc)(double);
int (*opsdfunc)(REAL8FrequencySeries *, double);
double srate = 16384; // sampling rate in Hertz
double segdur = 4; // duration of a segment in seconds
double tstart;
double duration;
double overrideflow;
double flow;
int official;
int psdonly;
const char *detector;

int usage(const char *program);
int parseargs(int argc, char **argv);

int main(int argc, char *argv[])
{
	size_t length;
	size_t stride;
	size_t n;
	REAL8FrequencySeries *psd = NULL;
	REAL8TimeSeries *seg = NULL;
	LIGOTimeGPS epoch;
	gsl_rng *rng;

	XLALSetErrorHandler(XLALAbortErrorHandler);

	parseargs(argc, argv);
	if (overrideflow > 0.0)
		flow = overrideflow;
	length = segdur * srate;
	stride = length / 2;

	XLALGPSSetREAL8(&epoch, tstart);
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(gsl_rng_default);
	psd = XLALCreateREAL8FrequencySeries(detector, &epoch, 0.0, srate/length, &lalSecondUnit, length/2 + 1);
	if (official && opsdfunc)
		opsdfunc(psd, flow);
	else
		XLALSimNoisePSD(psd, flow, psdfunc);
	if (psdonly) { // output PSD and exit
		size_t klow = flow / psd->deltaF;
		size_t k;
		for (k = klow; k < length/2 - 1; ++k)
			fprintf(stdout, "%e\t%e\n", k * psd->deltaF, sqrt(psd->data->data[k]));
		goto end;
	}

	n = duration * srate;
	seg = XLALCreateREAL8TimeSeries("STRAIN", &epoch, 0.0, 1.0/srate, &lalStrainUnit, length);
	XLALSimNoise(seg, 0, psd, rng); // first time to initialize
	while (1) { // infinite loop
		double t0 = XLALGPSGetREAL8(&seg->epoch);
		size_t j;
		for (j = 0; j < stride; ++j, --n) { // output first stride points
			if (n == 0) // check if we're done
				goto end;
			printf("%.9f\t%e\n", t0 + j * seg->deltaT, seg->data->data[j]);
		}
		XLALSimNoise(seg, stride, psd, rng); // make more data
	}

end:
	XLALDestroyREAL8TimeSeries(seg);
	XLALDestroyREAL8FrequencySeries(psd);
	LALCheckMemoryLeaks();

	return 0;
}

int parseargs( int argc, char **argv )
{
	struct option long_options[] = {
			{ "help", no_argument, 0, 'h' },
			{ "aligo-nosrm", no_argument, 0, 'A' },
			{ "aligo-zerodet-lowpower", no_argument, 0, 'B' },
			{ "aligo-zerodet-highpower", no_argument, 0, 'C' },
			{ "aligo-nsnsopt", no_argument, 0, 'D' },
			{ "aligo-bhbh20deg", no_argument, 0, 'E' },
			{ "aligo-highfreq", no_argument, 0, 'F' },
			{ "iligo-srd", no_argument, 0, 'I' },
			{ "virgo", no_argument, 0, 'v' },
			{ "advvirgo", no_argument, 0, 'V' },
			{ "geo", no_argument, 0, 'g' },
			{ "tama", no_argument, 0, 'T' },
			{ "kagra", no_argument, 0, 'K' },
			{ "official", no_argument, 0, 'O' },
			{ "psd-only", no_argument, 0, 'P' },
			{ "start-time", required_argument, 0, 's' },
			{ "duration", required_argument, 0, 't' },
			{ "sample-rate", required_argument, 0, 'r' },
			{ "segment-duration", required_argument, 0, 'd' },
			{ "low-frequency", required_argument, 0, 'f' },
			{ 0, 0, 0, 0 }
		};
	char args[] = "hIABCDEFOPvVgTKs:t:r:d:f:";
	while (1) {
		int option_index = 0;
		int c;

		c = getopt_long_only(argc, argv, args, long_options, &option_index);
		if (c == -1) /* end of options */
			break;

		switch (c) {
			case 0: /* if option set a flag, nothing else to do */
				if (long_options[option_index].flag)
					break;
				else {
					fprintf(stderr, "error parsing option %s with argument %s\n", long_options[option_index].name, optarg);
					exit(1);
				}
			case 'h': /* help */
				usage(argv[0]);
				exit(0);
			case 'A': /* aligo-nosrm */
				psdfunc = XLALSimNoisePSDaLIGONoSRMLowPower;
				opsdfunc = XLALSimNoisePSDaLIGONoSRMLowPowerGWINC;
				flow = 9.0;
				detector = "aLIGO";
				break;
			case 'B': /* aligo-zerodet-lowpower */
				psdfunc = XLALSimNoisePSDaLIGOZeroDetLowPower;
				opsdfunc = XLALSimNoisePSDaLIGOZeroDetLowPowerGWINC;
				flow = 9.0;
				detector = "aLIGO";
				break;
			case 'C': /* aligo-zerodet-highpower */
				psdfunc = XLALSimNoisePSDaLIGOZeroDetHighPower;
				opsdfunc = XLALSimNoisePSDaLIGOZeroDetHighPowerGWINC;
				flow = 9.0;
				detector = "aLIGO";
				break;
			case 'D': /* aligo-nsnsopt */
				psdfunc = XLALSimNoisePSDaLIGONSNSOpt;
				opsdfunc = XLALSimNoisePSDaLIGONSNSOptGWINC;
				flow = 9.0;
				detector = "aLIGO";
				break;
			case 'E': /* aligo-bhbh20deg */
				psdfunc = XLALSimNoisePSDaLIGOBHBH20Deg;
				opsdfunc = XLALSimNoisePSDaLIGOBHBH20DegGWINC;
				flow = 9.0;
				detector = "aLIGO";
				break;
			case 'F': /* aligo-highfreq */
				psdfunc = XLALSimNoisePSDaLIGOHighFrequency;
				opsdfunc = XLALSimNoisePSDaLIGOHighFrequencyGWINC;
				flow = 9.0;
				detector = "aLIGO";
				break;
			case 'I': /* iligo-srd */
				psdfunc = XLALSimNoisePSDiLIGOSRD;
				flow = 30.0;
				detector = "LIGO SRD";
				break;
			case 'v': /* initial Virgo */
				psdfunc = XLALSimNoisePSDVirgo;
				flow = 5.0;
				detector = "Virgo";
				break;
			case 'V': /* Advanced Virgo */
				psdfunc = XLALSimNoisePSDAdvVirgo;
				flow = 1.0;
				detector = "AdvVirgo";
				break;
			case 'g': /* GEO600 */
				psdfunc = XLALSimNoisePSDGEO;
				flow = 30.0;
				detector = "GEO600";
				break;
			case 'T': /* TAMA300 */
				psdfunc = XLALSimNoisePSDTAMA;
				flow = 30.0;
				detector = "TAMA300";
				break;
			case 'K': /* KAGRA (formerly LCGT) */
				psdfunc = XLALSimNoisePSDKAGRA;
				flow = 5.0;
				detector = "KAGRA";
				break;
			case 'O': /* official */
				official = 1;
				break;
			case 'P': /* start-time */
				psdonly = 1;
				break;
			case 's': /* start-time */
				tstart = atof(optarg);
				break;
			case 't': /* duration */
				duration = atof(optarg);
				break;
			case 'r': /* duration */
				srate = atof(optarg);
				break;
			case 'd': /* segment duration */
				segdur = atof(optarg);
				break;
			case 'f': /* low frequency */
				overrideflow = atof(optarg);
				break;
			case '?':
			default:
				fprintf(stderr, "unknown error while parsing options\n");
				exit(1);
		}
	}

	if ( optind < argc ) {
		fprintf(stderr, "extraneous command line arguments:\n");
		while (optind < argc)
			fprintf(stderr, "%s\n", argv[optind++]);
		exit(1);
	}

	if (! psdfunc || (!psdonly && duration == 0.0)) {
		fprintf(stderr, "must select a noise model and a duration\n");
		usage(argv[0]);
		exit(1);
	}

	return 0;
}

int usage( const char *program )
{
	fprintf(stderr, "usage: %s [options]\n", program);
	fprintf(stderr, "options:\n" );
	fprintf(stderr, "\t-h, --help                   \tprint this message and exit\n");
	fprintf(stderr, "\t-A, --aligo-nosrm            \taLIGO no SRM noise\n");
	fprintf(stderr, "\t-B, --aligo-zerodet-lowpower \taLIGO zero detuning low power noise\n");
	fprintf(stderr, "\t-C, --aligo-zerodet-highpower\taLIGO zero detuning high power noise\n");
	fprintf(stderr, "\t-D, --aligo-nsnsopt          \taLIGO NSNS optimized noise\n");
	fprintf(stderr, "\t-E, --aligo-bhbh20deg        \taLIGO BHBH optimized 20 deg detuning noise\n");
	fprintf(stderr, "\t-F, --aligo-highfreq         \taLIGO kHz narrowband noise\n");
	fprintf(stderr, "\t-I, --iligo-srd              \tiLIGO SRD noise power\n");
	fprintf(stderr, "\t-v, --virgo                  \tinitial Virgo noise power\n");
	fprintf(stderr, "\t-V, --advvirgo               \tAdvanced Virgo noise power\n");
	fprintf(stderr, "\t-g, --geo                    \tGEO600 noise power\n");
	fprintf(stderr, "\t-T, --tama                   \tTAMA300 noise power\n");
	fprintf(stderr, "\t-K, --kagra                  \tKAGRA noise power\n");
	fprintf(stderr, "\t-O, --official               \tuse official data files\n");
	fprintf(stderr, "\t-P, --psd-only               \toutput PSD only\n");
	fprintf(stderr, "\t-s, --start-time             \tGPS start time (s)\n");
	fprintf(stderr, "\t-t, --duration               \t(required) duration of data to produce (s)\n");
	fprintf(stderr, "\t-r, --sample-rate            \tsample rate (Hz) [16384]\n");
	fprintf(stderr, "\t-d, --segment-duration       \tsegment duration (s) [4]\n");
	fprintf(stderr, "\t-f, --low-frequency          \toverride default low frequency (Hz)\n");
	return 0;
}
