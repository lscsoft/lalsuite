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
 * @defgroup lalsim_bh_sphwf lalsim-bh-sphwf
 * @ingroup lalsimulation_programs
 *
 * @brief Evaluates a spin-weighted spheroidal wave function
 *
 * ### Synopsis
 *
 *     lalsim-bh-sphwf [-h] [-i theta] -a a -l l -m m -s s
 *
 * ### Description
 *
 * The `lalsim-bh-sphwf` utility prints the value of the spin-weighted
 * spheroidal wave function \f$ {}_{s}S_{\ell,m}(\cos\theta, a) \f$
 * for the specified dimensionless spin parameter @p a, spin weight @p s
 * mode numbers @p l and @p m, and polar angle @p theta.  If the parameter
 * @p theta is not given, the utility prints a table of the values of
 * the spin-weighted spheroidal wave function.
 * 
 * ### Options
 *
 * <DL>
 * <DT>`-h`, `--help`</DT>
 * <DD>print a help message and exit</DD>
 * <DT>`-i` theta</DT>
 * <DD>(optional) set the polar angle (degrees)</DD>
 * <DT>`-a` a</DT>
 * <DD>(required) set value of dimensionless spin parameter a/M, |a/M|<1</DD>
 * <DT>`-l` l</DT>
 * <DD>(required) set value of mode number l, l>=0</DD>
 * <DT>`-m` m</DT>
 * <DD>(required) set value of mode number m, abs(m)<=l</DD>
 * <DT>`-s` s</DT>
 * <DD>(required) set value of spin weight s, s<=0</DD>
 * </DL>
 *
 * ### Environment
 *
 * The `LAL_DEBUG_LEVEL` can used to control the error and warning reporting of
 * `lalsim-bh-sphwf`.  Common values are: `LAL_DEBUG_LEVEL=0` which suppresses
 * error messages, `LAL_DEBUG_LEVEL=1`  which prints error messages alone,
 * `LAL_DEBUG_LEVEL=3` which prints both error messages and warning messages,
 * and `LAL_DEBUG_LEVEL=7` which additionally prints informational messages.
 *
 * ### Exit Status
 *
 * The `lalsim-bh-sphwf` utility exits 0 on success, and >0 if an error
 * occurs.
 *
 * ### Example
 *
 * The command:
 *
 *     lalsim-bh-sphwf -a 0.97 -l 2 -m 2 -s -2
 *
 * prints a table of the spin-weighted spheroidal harmonic for
 * spin weight -2 (gravitational perturbations) in the l = m = 2 quasinormal
 * mode for black hole with Kerr spin parameter a/M = 0.97:
 *
@verbatim
# theta(deg)	    Re(S)   	    Im(S)
       0	1.872344e+00	-6.463929e-02
      10	1.827823e+00	-6.206959e-02
      20	1.694394e+00	-5.225983e-02
      30	1.504102e+00	-4.281635e-02
      40	1.273225e+00	-3.208053e-02
      50	1.027770e+00	-2.170445e-02
      60	7.909086e-01	-1.294164e-02
      70	5.808472e-01	-6.804991e-03
      80	4.043798e-01	-2.360517e-03
      90	2.668276e-01	-8.673617e-19
     100	1.656352e-01	9.846400e-04
     110	9.604459e-02	1.063797e-03
     120	5.104255e-02	8.091283e-04
     130	2.420241e-02	4.838700e-04
     140	9.739511e-03	2.305645e-04
     150	3.034816e-03	8.061289e-05
     160	5.924425e-04	1.697014e-05
     170	3.675324e-05	1.099041e-06
     180	0.000000e+00	0.000000e+00
@endverbatim
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
