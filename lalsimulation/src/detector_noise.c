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

/**
 * @defgroup lalsim_detector_noise lalsim-detector-noise
 * @ingroup lalsimulation_programs
 *
 * @brief Simulates detector noise
 *
 * ### Synopsis
 *
 *     lalsim-detector-noise [options]
 *
 * ### Description
 *
 * The `lalsim-detector-noise` utility produces a continuous stream of
 * simulated detector noise for a specified interval of time and for
 * a specified noise PSD.  Alternatively, `lalsim-detector-noise` outputs
 * a requested noise PSD.  The output is written to the standard output
 * in two-column ascii format data in which the first column contains either
 * the GPS times of each sample or the frequency of each PSD component,
 * and the second column contains the value of that sample.
 *
 * ### Options
 *
 * <DL>
 * <DT>`-h`, `--help`</DT>
 * <DD>print this message and exit</DD>
 * <DT>`-0`, `--0noise`</DT>
 * <DD>no noise (generates zeros)</DD>
 * <DT>`-A`, `--aligo-nosrm`</DT>
 * <DD>aLIGO no SRM noise</DD>
 * <DT>`-B`, `--aligo-zerodet-lowpower`</DT>
 * <DD>aLIGO zero detuning low power noise</DD>
 * <DT>`-C`, `--aligo-zerodet-highpower`</DT>
 * <DD>aLIGO zero detuning high power noise</DD>
 * <DT>`-D`, `--aligo-nsnsopt`</DT>
 * <DD>aLIGO NSNS optimized noise</DD>
 * <DT>`-E`, `--aligo-bhbh20deg`</DT>
 * <DD>aLIGO BHBH optimized 20 deg detuning noise</DD>
 * <DT>`-F`, `--aligo-highfreq`</DT>
 * <DD>aLIGO kHz narrowband noise</DD>
 * <DT>`-I`, `--iligo-srd`</DT>
 * <DD>iLIGO SRD noise power</DD>
 * <DT>`-v`, `--virgo`</DT>
 * <DD>initial Virgo noise power</DD>
 * <DT>`-V`, `--advvirgo`</DT>
 * <DD>Advanced Virgo noise power</DD>
 * <DT>`-g`, `--geo`</DT>
 * <DD>GEO600 noise power</DD>
 * <DT>`-G`, `--geohf`</DT>
 * <DD>GEO-HF noise power</DD>
 * <DT>`-T`, `--tama`</DT>
 * <DD>TAMA300 noise power</DD>
 * <DT>`-K`, `--kagra`</DT>
 * <DD>KAGRA noise power</DD>
 * <DT>`-O`, `--official`</DT>
 * <DD>use official data files</DD>
 * <DT>`-P`, `--psd-only`</DT>
 * <DD>output PSD only</DD>
 * <DT>`-s`, `--start-time` GPSSTART</DT>
 * <DD>GPS start time (s)</DD>
 * <DT>`-t`, `--duration` DURATION</DT>
 * <DD>(required) duration of data to produce (s)</DD>
 * <DT>`-r`, `--sample-rate` SRATE</DT>
 * <DD>sample rate (Hz) [16384]</DD>
 * <DT>`-d`, `--segment-duration` SEGDUR</DT>
 * <DD>segment duration (s) [4]</DD>
 * <DT>`-f`, `--low-frequency` FLOW</DT>
 * <DD>override default low frequency (Hz)</DD>
 * </DL>
 *
 * ### Environment
 *
 * The `LAL_DEBUG_LEVEL` can used to control the error and warning reporting of
 * `lalsim-detector-noise`.  Common values are: `LAL_DEBUG_LEVEL=0` which
 * suppresses error messages, `LAL_DEBUG_LEVEL=1`  which prints error messages
 * alone, `LAL_DEBUG_LEVEL=3` which prints both error messages and warning
 * messages, and `LAL_DEBUG_LEVEL=7` which additionally prints informational
 * messages.
 *
 * The `GSL_RNG_SEED` and `GSL_RNG_TYPE` environment variables can be used
 * to set the random number generator seed and type respectively.
 *
 * ### Exit Status
 *
 * The `lalsim-detector-noise` utility exits 0 on success, and >0 if an error
 * occurs.
 *
 * ### Example
 *
 * The command:
 *
 *     lalsim-detector-noise --aligo-zerodet-highpower -s 1000000000 -t 1000
 *
 * will stream 1000 seconds of aLIGO zero detuning high power noise
 * beginning at GPS time 1000000000.
 *
 * The command:
 *
 *     lalsim-detector-noise --iligo-srd -P
 *
 * outputs the Initial LIGO PSD.
 *
 * The command:
 *
 *     lalsim-detector-noise -0 -s 1000000000 -t 1000
 *
 * will stream 1000 seconds of zero-noise beginning at GPS time 1000000000.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>

#include <lal/LALStdlib.h>
#include <lal/LALgetopt.h>
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
LIGOTimeGPS tstart = LIGOTIMEGPSZERO;
double duration;
double overrideflow;
double flow;
int official;
int psdonly;
const char *detector;
const char *prefix;

int usage(const char *program);
int parseargs(int argc, char **argv);
double zeronoise(double f);

int main(int argc, char *argv[])
{
	char tstr[32]; // string to hold GPS time -- 31 characters is enough
	size_t length;
	size_t stride;
	size_t n;
	REAL8FrequencySeries *psd = NULL;
	REAL8TimeSeries *seg = NULL;
	gsl_rng *rng;

	XLALSetErrorHandler(XLALAbortErrorHandler);

	parseargs(argc, argv);
	if (overrideflow > 0.0)
		flow = overrideflow;
	length = segdur * srate;
	stride = length / 2;

	/* handle 0noise case first */
	if (strcmp(detector, "0noise") == 0) {
		/* just print out a bunch of zeros */
		if (psdonly) {
			double deltaF = srate / length;
			size_t klow = flow / deltaF;
			size_t k;
			fprintf(stdout, "# freq (s^-1)\tPSD (strain^2 s)\n");
			for (k = klow; k < length/2 - 1; ++k)
				fprintf(stdout, "%e\t%e\n", k * deltaF, 0.0);
		} else {
			size_t j;
			fprintf(stdout, "# time (s)\tNOISE (strain)\n");
			n = duration * srate;
			for (j = 0; j < n; ++j) { 
				LIGOTimeGPS t = tstart;
				fprintf(stdout, "%s\t%e\n", XLALGPSToStr(tstr, XLALGPSAdd(&t, j/srate)), 0.0);
			}
		}
		return 0;
	}

	gsl_rng_env_setup();
	rng = gsl_rng_alloc(gsl_rng_default);
	psd = XLALCreateREAL8FrequencySeries(detector, &tstart, 0.0, srate/length, &lalSecondUnit, length/2 + 1);
	if (official && opsdfunc)
		opsdfunc(psd, flow);
	else
		XLALSimNoisePSD(psd, flow, psdfunc);
	if (psdonly) { // output PSD and exit
		size_t klow = flow / psd->deltaF;
		size_t k;
		fprintf(stdout, "# freq (s^-1)\tPSD (strain^2 s)\n");
		for (k = klow; k < length/2 - 1; ++k)
			fprintf(stdout, "%e\t%e\n", k * psd->deltaF, sqrt(psd->data->data[k]));
		goto end;
	}

	n = duration * srate;
	seg = XLALCreateREAL8TimeSeries("STRAIN", &tstart, 0.0, 1.0/srate, &lalStrainUnit, length);
	XLALSimNoise(seg, 0, psd, rng); // first time to initialize
	fprintf(stdout, "# time (s)\tNOISE (strain)\n");
	while (1) { // infinite loop
		size_t j;
		for (j = 0; j < stride; ++j, --n) { // output first stride points
			LIGOTimeGPS t = seg->epoch;
			if (n == 0) // check if we're done
				goto end;
			fprintf(stdout, "%s\t%e\n", XLALGPSToStr(tstr, XLALGPSAdd(&t, j * seg->deltaT)), seg->data->data[j]);
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
	struct LALoption long_options[] = {
			{ "help", no_argument, 0, 'h' },
			{ "0noise", no_argument, 0, '0' },
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
			{ "geohf", no_argument, 0, 'G' },
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
	char args[] = "hI0ABCDEFOPvVgGTKs:t:r:d:f:";
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
			case '0': /* 0noise */
				/* psdfunc and opsdfunc are ignored so just choose anything */
				psdfunc = XLALSimNoisePSDaLIGONoSRMLowPower;
				opsdfunc = XLALSimNoisePSDaLIGONoSRMLowPowerGWINC;
				flow = 9.0;
				detector = "0noise";
				break;
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
			case 'G': /* GEO-HF */
				psdfunc = XLALSimNoisePSDGEOHF;
				flow = 50.0;
				detector = "GEOHF";
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
				{
					char *endp = NULL;
					if (XLALStrToGPS(&tstart, LALoptarg, &endp) < 0 || strlen(endp)) {
						fprintf(stderr, "could not parse GPS string `%s'\n", LALoptarg);
						exit(1);
					}
					break;
				}
			case 't': /* duration */
				duration = atof(LALoptarg);
				break;
			case 'r': /* sample-rate */
				srate = atof(LALoptarg);
				break;
			case 'd': /* segment duration */
				segdur = atof(LALoptarg);
				break;
			case 'f': /* low frequency */
				overrideflow = atof(LALoptarg);
				break;
			case '?':
			default:
				fprintf(stderr, "unknown error while parsing options\n");
				exit(1);
		}
	}

	if ( LALoptind < argc ) {
		fprintf(stderr, "extraneous command line arguments:\n");
		while (LALoptind < argc)
			fprintf(stderr, "%s\n", argv[LALoptind++]);
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
	fprintf(stderr, "\t-0, --0noise                 \tno noise (generates zeros)\n");
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
        fprintf(stderr, "\t-G, --geohf                  \tGEO-HF noise power\n");
	fprintf(stderr, "\t-T, --tama                   \tTAMA300 noise power\n");
	fprintf(stderr, "\t-K, --kagra                  \tKAGRA noise power\n");
	fprintf(stderr, "\t-O, --official               \tuse official data files\n");
	fprintf(stderr, "\t-P, --psd-only               \toutput PSD only\n");
	fprintf(stderr, "\t-s, --start-time GPSSTART    \tGPS start time (s)\n");
	fprintf(stderr, "\t-t, --duration DURATION      \t(required) duration of data to produce (s)\n");
	fprintf(stderr, "\t-r, --sample-rate SRATE      \tsample rate (Hz) [16384]\n");
	fprintf(stderr, "\t-d, --segment-duration SEGDUR\tsegment duration (s) [4]\n");
	fprintf(stderr, "\t-f, --low-frequency FLOW     \toverride default low frequency (Hz)\n");
	return 0;
}
