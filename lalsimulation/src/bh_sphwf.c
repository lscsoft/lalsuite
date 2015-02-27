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

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <lal/LALStdlib.h>
#include <lal/LALgetopt.h>
#include <lal/LALConstants.h>
#include <lal/LALSimBlackHoleRingdown.h>

#define theta_invalid -100.0 /* invalid */
#define a_invalid -100.0 /* invalid */
#define l_invalid (INT_MIN + 1) /* invalid */
#define m_invalid (INT_MAX) /* invalid */
#define s_invalid (INT_MAX) /* invalid */
double theta = theta_invalid;
double a = a_invalid;
int l = l_invalid;
int m = m_invalid;
int s = s_invalid;
int spherical = 0;
int plotfmt = 0;

int usage( const char *program );
int parseargs( int argc, char **argv );

int main( int argc, char *argv[] )
{
	XLALSetErrorHandler(XLALBacktraceErrorHandler);

	parseargs(argc, argv);

	if (theta == theta_invalid) { /* make a table */
		fprintf(stdout, "# theta(deg)\t    Re(S)   \t    Im(S)\n");
		for (theta = 0.0; theta <= 180.0; theta += 10.0) {
			COMPLEX16 sphwf = XLALSimBlackHoleRingdownSpheroidalWaveFunction(LAL_PI_180 * theta, a, l, m, s);
			fprintf(stdout, "%8g\t%e\t%e\n", theta, creal(sphwf), cimag(sphwf));
		}
	} else { /* evaluate at specified value */
		COMPLEX16 sphwf = XLALSimBlackHoleRingdownSpheroidalWaveFunction(LAL_PI_180 * theta, a, l, m, s);
		fprintf(stdout, "Spheroidal wave function (s)S(l,m)(cos(theta),a):\n");
		fprintf(stdout, "(%d)S(%d,%d)(cos(%g deg),%g) = %g + %g i\n", s, l, m, theta, a, creal(sphwf), cimag(sphwf));
	}
	return 0;
}

int parseargs(int argc, char **argv)
{
	struct LALoption long_options[] = {
			{ "help", no_argument, 0, 'h' },
			{ "spin", required_argument, 0, 'a' },
			{ "inclination", required_argument, 0, 'i' },
			{ "l", required_argument, 0, 'l' },
			{ "m", required_argument, 0, 'm' },
			{ "s", required_argument, 0, 's' },
			{ 0, 0, 0, 0 }
		};
	char args[] = "ha:i:l:m:s:";
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
			case 'a': /* spin */
				a = atof(LALoptarg);
				break;
			case 'i': /* inclination */
				theta = atof(LALoptarg);
				break;
			case 'l':
				l = atoi(LALoptarg);
				break;
			case 'm':
				m = atoi(LALoptarg);
				break;
			case 's':
				s = atoi(LALoptarg);
				break;
			case '?':
			default:
				fprintf(stderr, "unknown error while parsing options\n");
				exit(1);
		}
	}

	if (LALoptind < argc) {
		fprintf( stderr, "extraneous command line arguments:\n" );
		while (LALoptind < argc)
			fprintf(stderr, "%s\n", argv[LALoptind++]);
		exit(1);
	}

	if (a == a_invalid || l == l_invalid || m == m_invalid || s == s_invalid) {
		fprintf(stderr, "must specify a, l, m, and s\n");
		usage(argv[0]);
		exit(1);
	}

	if (fabs(a) >= 1.0) {
		fprintf(stderr, "must specify |a| < 1\n");
		exit(1);
	}

	if (l < abs(s)) {
		fprintf(stderr, "must specify l >= |s|\n");
		exit(1);
	}

	if (abs(m) > l) {
		fprintf(stderr, "must specify |m| <= l\n");
		exit(1);
	}

	return 0;
}

int usage(const char *program)
{
	fprintf(stderr, "usage: %s [options]\n", program);
	fprintf(stderr, "options:\n" );
	fprintf(stderr, "\t-h, --help     \tprint this message and exit\n");
	fprintf(stderr, "\t-i theta       \t(optional) inclination (polar) angle theta (deg)\n");
	fprintf(stderr, "\t-a a           \t(required) set value of a, -1<a<1\n");
	fprintf(stderr, "\t-l l           \t(required) set value of l, l>=0\n");
	fprintf(stderr, "\t-m m           \t(required) set value of m, abs(m)<=l\n");
	fprintf(stderr, "\t-s s           \t(required) set value of s, s<=0\n");
	return 0;
}
