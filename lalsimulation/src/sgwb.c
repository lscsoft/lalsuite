/*
*  Copyright (C) 2011 Jolien Creighton
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
#include <lal/LALSimSGWB.h>

double tstart;
double duration;
double flow = 10.0;
double Omega0;
LALDetector detectors[LAL_NUM_DETECTORS];
size_t numDetectors;

int usage(const char *program);
int parseargs(int argc, char **argv);

int main(int argc, char *argv[])
{
	const double H0 = 0.72 * LAL_H0FAC_SI; // Hubble's constant in seconds
	const double srate = 16384.0; // sampling rate in Hertz
	const size_t length = 65536; // number of points in a segment
	const size_t stride = length / 2; // number of points in a stride
	size_t i, n;
	REAL8FrequencySeries *OmegaGW = NULL;
	REAL8TimeSeries **seg = NULL;
	LIGOTimeGPS epoch;
	gsl_rng *rng;

	XLALSetErrorHandler(XLALAbortErrorHandler);

	parseargs(argc, argv);

	XLALGPSSetREAL8(&epoch, tstart);
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(gsl_rng_default);
	OmegaGW = XLALSimSGWBOmegaGWFlatSpectrum(Omega0, flow, srate/length, length/2 + 1);

	n = duration * srate;
	seg = LALCalloc(numDetectors, sizeof(*seg));
	for (i = 0; i < numDetectors; ++i)
		seg[i] = XLALCreateREAL8TimeSeries("STRAIN", &epoch, 0.0, 1.0/srate, &lalStrainUnit, length);

	XLALSimSGWB(seg, detectors, numDetectors, 0, OmegaGW, H0, rng); // first time to initilize
	while (1) { // infinite loop
		double t0 = XLALGPSGetREAL8(&seg[0]->epoch);
		size_t j;
		for (j = 0; j < stride; ++j, --n) { // output first stride points
			if (n == 0) // check if we're done
				goto end;
			printf("%.9f", t0 + j * seg[0]->deltaT);
			for (i = 0; i < numDetectors; ++i)
				printf("\t%e", seg[i]->data->data[j]);
			printf("\n");
		}
		XLALSimSGWB(seg, detectors, numDetectors, stride, OmegaGW, H0, rng); // make more data
	}

end:
	for (i = 0; i < numDetectors; ++i)
		XLALDestroyREAL8TimeSeries(seg[i]);
	XLALFree(seg);
	XLALDestroyREAL8FrequencySeries(OmegaGW);
	LALCheckMemoryLeaks();

	return 0;
}

int parseargs( int argc, char **argv )
{
	struct option long_options[] = {
			{ "help", no_argument, 0, 'h' },
			{ "geo", no_argument, 0, 'G' },
			{ "hanford", no_argument, 0, 'H' },
			{ "livingston", no_argument, 0, 'L' },
			{ "virgo", no_argument, 0, 'V' },
			{ "start-time", required_argument, 0, 's' },
			{ "duration", required_argument, 0, 't' },
			{ "Omega0", required_argument, 0, 'W' },
			{ "low-frequency", required_argument, 0, 'f' },
			{ 0, 0, 0, 0 }
		};
	char args[] = "hGHLVs:t:W:f:";
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
			case 'G': /* geo */
				detectors[numDetectors++] = lalCachedDetectors[LAL_GEO_600_DETECTOR];
				break;
			case 'H': /* hanford */
				detectors[numDetectors++] = lalCachedDetectors[LAL_LHO_4K_DETECTOR];
				break;
			case 'L': /* livingston */
				detectors[numDetectors++] = lalCachedDetectors[LAL_LLO_4K_DETECTOR];
				break;
			case 'V': /* livingston */
				detectors[numDetectors++] = lalCachedDetectors[LAL_VIRGO_DETECTOR];
				break;
			case 's': /* start-time */
				tstart = atof(optarg);
				break;
			case 't': /* duration */
				duration = atof(optarg);
				break;
			case 'W': /* Omega0 */
				Omega0 = atof(optarg);
				break;
			case 'f': /* low-frequency */
				flow = atof(optarg);
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

	if (duration == 0.0 || Omega0 == 0.0) {
		fprintf(stderr, "must select a duration and a value of Omega0\n");
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
	fprintf(stderr, "\t-G, --geo                    \tinclude GEO\n");
	fprintf(stderr, "\t-H, --hanford                \tinclude LHO\n");
	fprintf(stderr, "\t-L, --livingston             \tinclude LLO\n");
	fprintf(stderr, "\t-V, --virgo                  \tinclude Virgo\n");
	fprintf(stderr, "\t-s, --start-time GPSSTART    \tGPS start time (s)\n");
	fprintf(stderr, "\t-t, --duration   DURATION    \t(required) duration of data to produce (s)\n");
	fprintf(stderr, "\t-W, --Omega0     OMEGA0      \t(required) flat spectral energy density\n");
	fprintf(stderr, "\t-f, --low-frequency FLOW     \tlow frequency cutoff (Hz) (default = 10 Hz)\n");
	return 0;
}
