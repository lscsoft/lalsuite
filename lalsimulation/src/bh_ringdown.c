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
 * @defgroup lalsim_bh_ringdown lalsim-bh-ringdown
 * @ingroup lalsimulation_programs
 *
 * @brief Simulates a gravitational waveform from black hole ringdown.
 *
 * ### Synopsis
 *
 *     lalsim-bh-ringdown [-h] -M Msolar -a a -r distanceMpc -e fracEnergy -i inclination [-q azimuth] -l l -m m
 *
 * ### Description
 *
 * The `lalsim-bh-ringdown` utility produces a stream of a simulated
 * gravitational waveform for ringdown radiation from a quasinormal
 * mode of a Kerr black hole with mode numbers @p l and @p m.  The
 * dimensionless Kerr spin parameter @p a, black hole mass in solar
 * masses @p Msolar, fraction of mass lost in ringdown radiation
 * @p e, distance to the observer in Mpc @p distanceMpc, and inclination of the
 * observer relative to the black hole's spin axis @p inclination must be
 * specified.  The output is written to standard output in three-column ascii
 * format.  The first column gives the time corresponding to each sample, the
 * second column gives the value of the plus-polarization of the waveform, and
 * the third column gives the value of the cross-polarization of the waveform.
 * 
 * ### Options
 *
 * <DL>
 * <DT>`-h`, `--help`</DT>
 * <DD>print a help message and exit</DD>
 * <DT>`-M` Msolar</DT>
 * <DD>(required) set black hole mass (solar masses)</DD>
 * <DT>`-a` a</DT>
 * <DD>(required) set value of dimensionless spin parameter a/M, |a/M|<1 (Leaver: |a/M|<0.5)</DD>
 * <DT>`-r` distanceMpc</DT>
 * <DD>(required) set distance (Mpc)</DD>
 * <DT>`-e` fracEnergy</DT>
 * <DD>(required) set energy radiated (fraction of mass)</DD>
 * <DT>`-i` inclination</DT>
 * <DD>(required) set inclination angle (degrees)</DD>
 * <DT>`-q` azimuth</DT>
 * <DD>(optional: default=0) set azimuth angle (degrees)</DD>
 * <DT>`-l` l</DT>
 * <DD>(required) set value of mode number l, l>=0</DD>
 * <DT>`-m` m</DT>
 * <DD>(required) set value of mode number m, abs(m)<=l</DD>
 * </DL>
 *
 * ### Environment
 *
 * The `LAL_DEBUG_LEVEL` can used to control the error and warning reporting of
 * `lalsim-bh-ringdown`.  Common values are: `LAL_DEBUG_LEVEL=0` which
 * suppresses error messages, `LAL_DEBUG_LEVEL=1`  which prints error messages
 * alone, `LAL_DEBUG_LEVEL=3` which prints both error messages and warning
 * messages, and `LAL_DEBUG_LEVEL=7` which additionally prints informational
 * messages.
 *
 * ### Exit Status
 *
 * The `lalsim-bh-ringdown` utility exits 0 on success, and >0 if an error
 * occurs.
 *
 * ### Example
 *
 * The command:
 *
 *     lalsim-bh-ringdown -M 10 -a 0.97 -r 1.0 -e 0.01 -i 45.0 -l 2 -m 2
 *
 * produces a three-column ascii output to standard output; the rows are
 * samples (at the rate of 16384 Hz), and the three columns are 1. the
 * time of each sample, 2. the plus-polarization strain, and 3. the
 * cross-polarization strain.  The waveform produced is for a 10 solar
 * mass black hole spinning with Kerr parameter a/M = 0.97 at a distance
 * of 1 Mpc and inclination of 45 degrees that radiates 1% of its mass in the l
 * = 2, m = 2 quasinormal mode.
 */

#include <complex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <lal/LALStdlib.h>
#include <lal/LALgetopt.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/LALSimBlackHoleRingdown.h>

#define a_invalid -100.0 /* invalid */
#define i_invalid -361.0 /* invalid */
#define l_invalid (INT_MIN + 1) /* invalid */
#define m_invalid (INT_MAX) /* invalid */
double dt = 1.0/16384.0;
double a = a_invalid;
double M = 0.0;
double r = 0.0;
double e = 0.0;
double i = i_invalid;
double q = 0.0;
int l = l_invalid;
int m = m_invalid;

int usage(const char *program);
int parseargs(int argc, char **argv);

int main(int argc, char *argv[])
{
  	LIGOTimeGPS epoch = {0, 0};
	REAL8TimeSeries *hplus;
	REAL8TimeSeries *hcross;
	size_t j;

	XLALSetErrorHandler(XLALBacktraceErrorHandler);

	parseargs(argc, argv);

	XLALSimBlackHoleRingdown(&hplus, &hcross, &epoch, q, dt, M, a, e, r, i, l, m);

	fprintf(stdout, "# time (s)\th_plus (strain)\th_cross (strain)\n");
	for (j = 0; j < hplus->data->length; ++j)
		fprintf(stdout, "%.9f\t%e\t%e\n", j*dt, hplus->data->data[j], hcross->data->data[j]);

	XLALDestroyREAL8TimeSeries(hcross);
	XLALDestroyREAL8TimeSeries(hplus);
	LALCheckMemoryLeaks();

	return 0;
}

int parseargs( int argc, char **argv )
{
	struct LALoption long_options[] = {
			{ "help", no_argument, 0, 'h' },
			{ "leaver", no_argument, 0, 'L' },
			{ "mass", required_argument, 0, 'M' },
			{ "spin", required_argument, 0, 'a' },
			{ "inclination", required_argument, 0, 'i' },
			{ "azimuth", required_argument, 0, 'q' },
			{ "energy", required_argument, 0, 'e' },
			{ "distance", required_argument, 0, 'r' },
			{ "l", required_argument, 0, 'l' },
			{ "m", required_argument, 0, 'm' },
			{ 0, 0, 0, 0 }
		};
	char args[] = "hM:a:i:q:e:r:l:m:";
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
			case 'M': /* mass */
				M = LAL_MSUN_SI * atof(LALoptarg);
				break;
			case 'a': /* spin */
				a = atof(LALoptarg);
				break;
			case 'i': /* inclination */
				i = LAL_PI_180 * atof(LALoptarg);
				break;
			case 'q': /* azimuth */
				q = LAL_PI_180 * atof( LALoptarg );
				break;
			case 'e': /* energy */
				e = atof(LALoptarg);
				break;
			case 'r': /* distance */
				r = 1e6 * LAL_PC_SI * atof(LALoptarg);
				break;
			case 'l':
				l = atoi(LALoptarg);
				break;
			case 'm':
				m = atoi(LALoptarg);
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
	
	if (a == a_invalid || l == l_invalid || m == m_invalid || M <= 0.0 || e <= 0.0 || e >= 1.0 || r <= 0.0 || i == i_invalid) {
		fprintf(stderr, "must specify mass, spin, distance, frac. energy loss, l, m\n");
		usage(argv[0]);
		exit(1);
	}

	if (l < 2) {
		fprintf(stderr, "must specify l >= 2\n");
		exit(1);
	}

	if (abs(m) > l) {
		fprintf(stderr, "must specify m <= l\n");
		exit(1);
	}
	
	return 0;
}
	
int usage( const char *program )
{
	fprintf(stderr, "usage: %s [options]\n", program);
	fprintf(stderr, "options:\n" );
	fprintf(stderr, "\t-h, --help     \tprint this message and exit\n");
	fprintf(stderr, "\t-M Msolar      \t(required) set black hole mass (solar masses)\n");
	fprintf(stderr, "\t-a a           \t(required) set value of a, -1<a<1\n");
	fprintf(stderr, "\t-r distanceMpc \t(required) set distance (Mpc)\n");
	fprintf(stderr, "\t-e fracEnergy  \t(required) set energy radiated (fraction of M)\n");
	fprintf(stderr, "\t-i inclination \t(required) set inclination angle (degrees)\n");
	fprintf(stderr, "\t-q azimuth     \t(default=0) set azimuth angle (degrees)\n");
	fprintf(stderr, "\t-l l           \t(required) set value of l, l>=2\n");
	fprintf(stderr, "\t-m m           \t(required) set value of m, abs(m)<=l\n");
	return 0;
}
