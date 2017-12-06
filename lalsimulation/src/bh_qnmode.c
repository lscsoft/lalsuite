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
#include <getopt.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALSimBlackHoleRingdown.h>

#define a_invalid -100.0 /* invalid */
#define l_invalid (INT_MIN + 1) /* invalid */
#define m_invalid (INT_MAX) /* invalid */
#define s_invalid (INT_MAX) /* invalid */
REAL8 M = 0.0;
REAL8 a = a_invalid;
int l = l_invalid;
int m = m_invalid;
int s = s_invalid;
int leaver = 0;

int usage(const char *program);
int parseargs(int argc, char **argv);
int output_mode_table(void);
int output_mode(void);


int output_mode_table(void)
{
	double avec[] = {0.0,0.1,0.2,0.3,0.4,0.45,0.49,0.4999};
	int asign;
	size_t numa = sizeof(avec) / sizeof(*avec);
	size_t i;
	
	fprintf(stdout, "# quasinormal mode table for l=%d m=%d s=%d ", l, m, s);
	if (leaver)
		fprintf(stdout, "(Leaver's conventions)\n");
	else
		fprintf(stdout, "(standard conventions)\n");
	if (M == 0.0)
		fprintf(stdout, "#  a   \t        %s        \t        omega\n", leaver ? "A" : "E");
	else
		fprintf(stdout, "#  a   \t  frequency (Hz)  \t   quality   \n");

	/* do positive spins first, then do negative spins */
	for (asign = 1; abs(asign) == 1; asign -= 2) {
		for (i = 0; i < numa; ++i) {
			COMPLEX16 A, omega;
    			a = avec[i];
    			XLALSimBlackHoleRingdownModeEigenvaluesLeaver(&A, &omega, asign*a, l, m, s);
    			if (!leaver) { /* convert to standard conventions */
				a *= 2.0;
				omega *= 0.5;
				A += s*(s+1);
    			}
			if (M == 0.0)
				fprintf(stdout, "%.4f \t(%.5f,%+.5f)\t(%+.6f,%.6f)\n", a, creal(A), cimag(A), creal(omega), cimag(omega));
			else {
				double f, Q;
				f = fabs(creal(omega) / (LAL_TWOPI * M * LAL_MTSUN_SI));
				Q = fabs(creal(omega)) / (-2.0 * cimag(omega));
				fprintf(stdout, "%.4f \t %12.3f \t\t %8.3f\n", a, f, Q);
			}
		}
		if (asign == 1)
			fprintf(stdout, "\n");
	}
	return 0;
}

int output_mode(void)
{
	COMPLEX16 A, omega;

	if (!leaver) /* spin was not specified in Leaver conventions */
		a *= 0.5; /* change to Leaver conventions */

    	XLALSimBlackHoleRingdownModeEigenvaluesLeaver(&A, &omega, a, l, m, s);

	if (!leaver) { /* change from Leaver conventions */
		omega *= 0.5;
		a *= 2.0;
	}

	/*
	 * if mass has been specified, return f and Q
	 * otherwise return omega in either leaver or standard conventions
	 */

	fprintf(stdout, "using %s conventions\n", leaver ? "Leaver's" : "standard");
	fprintf(stdout, "mode l = %d, m = %d, s = %d\n", l, m, s);
	fprintf(stdout, "spin a = %g (dimensionless)\n", a);
	fprintf(stdout, "M * omega = (%+.6f,%+.6f)\n", creal(omega), cimag(omega));
	if (M != 0.0) {
		double f, Q;
		f = fabs(creal(omega) / (LAL_TWOPI * M * LAL_MTSUN_SI));
		Q = fabs(creal(omega)) / (-2.0 * cimag(omega));

		fprintf(stdout, "mass M = %g solar masses\n", M);
		fprintf(stdout, "frequency f = %g Hz\n", f);
		fprintf(stdout, "quality Q = %g\n", Q);
	}
	
	return 0;
}


int main(int argc, char *argv[])
{
	XLALSetErrorHandler(XLALBacktraceErrorHandler);

	parseargs(argc, argv);

	if (a == a_invalid)
		output_mode_table();
	else
		output_mode();

	LALCheckMemoryLeaks();
	return 0;
}

int parseargs( int argc, char **argv )
{
	struct option long_options[] = {
			{ "help", no_argument, 0, 'h' },
			{ "leaver", no_argument, 0, 'L' },
			{ "mass", required_argument, 0, 'M' },
			{ "spin", required_argument, 0, 'a' },
			{ "l", required_argument, 0, 'l' },
			{ "m", required_argument, 0, 'm' },
			{ "s", required_argument, 0, 's' },
			{ 0, 0, 0, 0 }
		};
	char args[] = "hLM:a:l:m:s:";
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
		case 'L': /* leaver */
			leaver = 1;
			break;
		case 'M': /* mass */
			M = atof(optarg);
			break;
		case 'a': /* spin */
			a = atof(optarg);
			break;
		case 'l':
			l = atoi(optarg);
			break;
		case 'm':
			m = atoi(optarg);
			break;
		case 's':
			s = atoi(optarg);
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

	if (l == l_invalid || m == m_invalid || s == s_invalid) {
		fprintf(stderr, "must specify l, m, and s\n");
		usage(argv[0]);
		exit(1);
	}

	if (!(a == a_invalid || (leaver && fabs(a) < 0.5) || (!leaver && fabs(a) < 1.0))) {
		fprintf(stderr, "must specify |a| < 1 (or |a| < 0.5 with Leaver's conventions)\n");
		exit(1);
	}

	if (leaver && M != 0.0) {
		fprintf(stderr, "do not use both --leaver and --mass options\n");
		usage(argv[0]);
		exit(1);
	}

	return 0;
}

int usage(const char *program)
{
	fprintf(stderr, "usage: %s [options]\n", program);
	fprintf(stderr, "options:\n");
	fprintf(stderr, "\t-h, --help     \tprint this message and exit\n");
	fprintf(stderr, "\t-L, --leaver   \tuse Leaver's conventions\n");
	fprintf(stderr, "\t-M Msolar      \t(optional) set black hole mass (solar masses)\n");
	fprintf(stderr, "\t-a a           \t(optional) set value of a, |a|<1 (Leaver: |a|<0.5)\n");
	fprintf(stderr, "\t-l l           \t(required) set value of l, l>=0\n");
	fprintf(stderr, "\t-m m           \t(required) set value of m, abs(m)<=l\n");
	fprintf(stderr, "\t-s s           \t(required) set value of s, s<=0\n");
	fprintf(stderr, "description:\n");
	fprintf(stderr, "\tprints the eigenvalues for given a if a is specified\n");
	fprintf(stderr, "\tprints a table of mode eigenvalues if a is not specified\n");
	fprintf(stderr, "\tprints frequency and quality factor if mass is specified\n");
	fprintf(stderr, "\tuses Leaver's conventions (G = c = 2M = 1) if --leaver is used\n");
	return 0;
}
