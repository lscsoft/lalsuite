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

/**
 * @defgroup lalsim_sgwb lalsim-sgwb
 * @ingroup lalsimulation_programs
 *
 * @brief Simulates a stochastic gravitational wave background
 *
 * ### Synopsis
 *
 *     lalsim-sgwb [options]
 *
 * ### Description
 *
 * The `lalsim-sgwb` utility produces a continuous stream of simulated detector
 * stochastic gravitational wave background noise for a specified interval of
 * time, for the specified detectors, and for a specified flat spectral energy
 * density.  The output is written to the standard output in multi-column ascii
 * format data in which the first column contains the GPS times of each sample
 * and the remaining columns contain the (correlated) noise strain values
 * in the various detectors.
 *
 * ### Options
 *
 * <DL>
 * <DT>`-h`, `--help`</DT>                   	<DD>print a help message and exit</DD>
 * <DT>`-G`, `--geo`</DT>                    	<DD>include GEO</DD>
 * <DT>`-H`, `--hanford`</DT>                	<DD>include LHO</DD>
 * <DT>`-L`, `--livingston`</DT>             	<DD>include LLO</DD>
 * <DT>`-V`, `--virgo`</DT>                  	<DD>include Virgo</DD>
 * <DT>`-s`, `--start-time` GPSSTART</DT>    	<DD>GPS start time (s)</DD>
 * <DT>`-t`, `--duration`   DURATION</DT>    	<DD>(required) duration of data to produce (s)</DD>
 * <DT>`-r`, `--sample-rate` SRATE</DT>            	<DD>sample rate (Hz) [16384]</DD>
 * <DT>`-W`, `--Omega0`     OMEGA0</DT>      	<DD>(required) flat spectral energy density</DD>
 * <DT>`-f`, `--low-frequency` FLOW</DT>     	<DD>low frequency cutoff (Hz) (default = 10 Hz)</DD>
 * </DL>
 *
 * ### Environment
 *
 * The `LAL_DEBUG_LEVEL` can used to control the error and warning reporting of
 * `lalsim-sgwb`.  Common values are: `LAL_DEBUG_LEVEL=0` which suppresses
 * error messages, `LAL_DEBUG_LEVEL=1`  which prints error messages alone,
 * `LAL_DEBUG_LEVEL=3` which prints both error messages and warning messages,
 * and `LAL_DEBUG_LEVEL=7` which additionally prints informational messages.
 *
 * The `GSL_RNG_SEED` and `GSL_RNG_TYPE` environment variables can be used
 * to set the random number generator seed and type respectively.
 *
 * ### Exit Status
 *
 * The `lalsim-sgwb` utility exits 0 on success, and >0 if an error occurs.
 *
 * ### Example
 *
 * The command:
 *
 *     lalsim-sgwb -H -L -V -W 1e-6 -s 1000000000 -t 1000
 *
 * will stream 1000 seconds of stochastic gravitational-wave background
 * strain noise in the LHO, LLO, and Virgo detectors having
 * \f$ {\Omega_0=10^{-6}} \f$ .
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
#include <lal/LALSimSGWB.h>

double srate = 16384.0; // sampling rate in Hertz
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
	char tstr[32]; // string to hold GPS time -- 31 characters is enough
	const double H0 = 0.72 * LAL_H0FAC_SI; // Hubble's constant in seconds
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
	printf("# time (s)");
	for (i = 0; i < numDetectors; ++i) {
		char name[LALNameLength];
		snprintf(name, sizeof(name), "%s:STRAIN", detectors[i].frDetector.prefix);
		seg[i] = XLALCreateREAL8TimeSeries(name, &epoch, 0.0, 1.0/srate, &lalStrainUnit, length);
		printf("\t%s (strain)", name);
	}
	printf("\n");

	XLALSimSGWB(seg, detectors, numDetectors, 0, OmegaGW, H0, rng); // first time to initilize

	while (1) { // infinite loop
		size_t j;
		for (j = 0; j < stride; ++j, --n) { // output first stride points
			LIGOTimeGPS t = seg[0]->epoch;
			if (n == 0) // check if we're done
				goto end;
			printf("%s", XLALGPSToStr(tstr, XLALGPSAdd(&t, j * seg[0]->deltaT)));
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
	struct LALoption long_options[] = {
			{ "help", no_argument, 0, 'h' },
			{ "geo", no_argument, 0, 'G' },
			{ "hanford", no_argument, 0, 'H' },
			{ "livingston", no_argument, 0, 'L' },
			{ "virgo", no_argument, 0, 'V' },
			{ "start-time", required_argument, 0, 's' },
			{ "duration", required_argument, 0, 't' },
			{ "sample-rate", required_argument, 0, 'r' },
			{ "Omega0", required_argument, 0, 'W' },
			{ "low-frequency", required_argument, 0, 'f' },
			{ 0, 0, 0, 0 }
		};
	char args[] = "hGHLVs:t:r:W:f:";
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
				tstart = atof(LALoptarg);
				break;
			case 't': /* duration */
				duration = atof(LALoptarg);
				break;
			case 'r': /* sample-rate */
				srate = atof(LALoptarg);
				break;
			case 'W': /* Omega0 */
				Omega0 = atof(LALoptarg);
				break;
			case 'f': /* low-frequency */
				flow = atof(LALoptarg);
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
	fprintf(stderr, "\t-s, --start-time    GPSSTART \tGPS start time (s)\n");
	fprintf(stderr, "\t-t, --duration      DURATION \t(required) duration of data to produce (s)\n");
	fprintf(stderr, "\t-r, --sample-rate   SRATE    \tsample rate (Hz) [16384]\n");
	fprintf(stderr, "\t-W, --Omega0        OMEGA0   \t(required) flat spectral energy density\n");
	fprintf(stderr, "\t-f, --low-frequency FLOW     \tlow frequency cutoff (Hz) (default = 10 Hz)\n");
	return 0;
}
