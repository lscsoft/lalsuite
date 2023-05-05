/*
 * Copyright (C) 2013 J. Creighton, B. Lackey
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <lal/LALConstants.h>
#include <lal/LALgetopt.h>
#include <lal/LALStdlib.h>
#include <lal/LALSimNeutronStar.h>

int usage(const char *program);
int parseargs(int argc, char *argv[]);
int output(const char *fmt, double c, double m, double r, double w, double z, double k2); 

/* global variables  */

LALSimNeutronStarEOS *global_eos = NULL;

const char *const default_fmt = "%r\\t%m\\n";
const char *const default_fmt_U = "%r\\t%m\\t%w\\t%z\\n";
const char *global_fmt;

double global_epsrel = 0;

int global_virial = 0;
double global_rho_c, global_e_c = 0;

const long default_npts = 100;
long global_npts;

int main(int argc, char *argv[])
{
    double logpmin = 75.5;
    double logpmax;
    double dlogp;
    long i;

    XLALSetErrorHandler(XLALAbortErrorHandler);

    global_npts = default_npts;
    global_fmt = default_fmt;
    parseargs(argc, argv);

    if (global_rho_c == 0 && global_e_c != 0)
    {
        logpmin = log(XLALSimNeutronStarEOSPressureOfEnergyDensity(global_e_c, global_eos));
    }

    if (global_rho_c != 0 && global_e_c == 0)
    {
        logpmin = log(XLALSimNeutronStarEOSPressureOfRestMassDensity(global_rho_c, global_eos));
    }

    logpmax = log(XLALSimNeutronStarEOSMaxPressure(global_eos));
    dlogp = (logpmax - logpmin) / global_npts;

    for (i = 0; i < global_npts; ++i) {
        double pc = exp(logpmin + (0.5 + i) * dlogp);
        double c, m, r, k2, I1, I2, I3, J1, J2, J3, w, z;
        if (global_virial == 0)
        {
            if (global_epsrel == 0)
            {
                XLALSimNeutronStarTOVODEIntegrate(&r, &m, &k2, pc, global_eos);
            }
            else
            {
                XLALSimNeutronStarTOVODEIntegrateWithTolerance(&r, &m, &k2, pc, global_eos, global_epsrel);
            }
            /* convert units */
            m /= LAL_MSUN_SI;       /* mass in solar masses */
            c = m * LAL_MRSUN_SI / r;       /* compactness (dimensionless) */
            r /= 1000.0;    /* radius in km */
            w = z = 0.0 / 0.0; /* GRV2 and GRV3 (dimesionless) NaN since -U is not used */

            output(global_fmt, c, m, r, w, z, k2);
        }
        else if (global_virial == 1)
        {
            if (global_epsrel == 0)
            {
                XLALSimNeutronStarVirialODEIntegrate(&r, &m, &I1, &I2, &I3, &J1, &J2, &J3, &k2, pc, global_eos);
            }
            else
            {
                XLALSimNeutronStarVirialODEIntegrateWithTolerance(&r, &m, &I1, &I2, &I3, &J1, &J2, &J3, &k2, pc, global_eos, global_epsrel);
            }
            /* convert units */
            m /= LAL_MSUN_SI;       /* mass in solar masses */
            c = m * LAL_MRSUN_SI / r;       /* compactness (dimensionless) */
            r /= 1000.0;    /* radius in km */
            w = fabs(I1/(I2 + I3) - 1.0); /* GRV2 (dimesionless) */
            z = fabs(J1/(J2 + J3) - 1.0); /* GRV2 (dimesionless) */

            output(global_fmt, c, m, r, w, z, k2);
        }   
    }

    XLALDestroySimNeutronStarEOS(global_eos);
    LALCheckMemoryLeaks();
    return 0;
}

int output(const char *fmt, double c, double m, double r, double w, double z, double k2)
{
    int i;
    for (i = 0; fmt[i]; ++i) {
        switch (fmt[i]) {
        case '%':      /* conversion character */
            switch (fmt[i + 1]) {
            case '%':
                fputc('%', stdout);
                break;
            case 'c':
                fprintf(stdout, "%.18e", c);
                break;
            case 'k':
                fprintf(stdout, "%.18e", k2);
                break;
            case 'l':
                fprintf(stdout, "%.18e", (2.0 / 3.0) * k2 / pow(c, 5));
                break;
            case 'm':
                fprintf(stdout, "%.18e", m);
                break;
            case 'r':
                fprintf(stdout, "%.18e", r);
                break;
            case 'w':
                fprintf(stdout, "%.5e", w);
                break;
            case 'z':
                fprintf(stdout, "%.5e", z);
                break;
            default:
                fprintf(stderr,
                    "unrecognized conversion %%%c in output format\n",
                    fmt[i + 1]);
                exit(1);
                break;
            }
            ++i;
            break;
        case '\\':     /* escaped character */
            switch (fmt[i + 1]) {
            case '\\':
                fputc('\\', stdout);
                break;
            case 'a':
                fputc('\a', stdout);
                break;
            case 'b':
                fputc('\b', stdout);
                break;
            case 'f':
                fputc('\f', stdout);
                break;
            case 'n':
                fputc('\n', stdout);
                break;
            case 'r':
                fputc('\r', stdout);
                break;
            case 't':
                fputc('\t', stdout);
                break;
            case 'v':
                fputc('\v', stdout);
                break;
            case '\0': /* impossible! */
                fprintf(stderr,
                    "impossible escaped NUL character in format\n");
                exit(1);
            default:
                fputc(fmt[i + 1], stdout);
                break;
            }
            ++i;
            break;
        default:
            fputc(fmt[i], stdout);
            break;
        }
    }
    return 0;
}

int parseargs(int argc, char **argv)
{
    struct LALoption long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"mass-radius-virial", no_argument, 0, 'U'},
        {"epsilon_relative", required_argument, 0, 'E'},
        {"eos-file", required_argument, 0, 'f'},
        {"eos-name", required_argument, 0, 'n'},
        {"central-rest-mass-density", required_argument, 0, 'd'},
        {"central-energy-density", required_argument, 0, 'e'},
        {"format", required_argument, 0, 'F'},
        {"npts", required_argument, 0, 'N'},
        {"polytrope", no_argument, 0, 'P'},
        {"gamma", required_argument, 0, 'G'},
        {"pressure", required_argument, 0, 'p'},
        {"density", required_argument, 0, 'r'},
        {"piecewisepolytrope", no_argument, 0, 'Q'},
        {"logp1", required_argument, 0, 'q'},
        {"gamma1", required_argument, 0, '1'},
        {"gamma2", required_argument, 0, '2'},
        {"gamma3", required_argument, 0, '3'},
        {"4parameterspectraldecomposition", no_argument, 0, 'S'},
        {"SDgamma0", required_argument, 0, 'w'},
        {"SDgamma1", required_argument, 0, 'x'},
        {"SDgamma2", required_argument, 0, 'y'},
        {"SDgamma3", required_argument, 0, 'z'},
        {0, 0, 0, 0}
    };
    char args[] = "hUE:f:n:d:e:F:N:PG:p:r:Qq:1:2:3:Sw:x:y:z:";

    /* quantities for 1-piece polytrope: */
    int polytropeFlag = 0;
    double Gamma = 0, reference_pressure_si = 0, reference_density_si = 0;

    /* quantities for 4-parameter piecewise polytrope: */
    int piecewisePolytropeFlag = 0;
    double logp1_si = 0, gamma1 = 0, gamma2 = 0, gamma3 = 0;

    /* quatities for 4-coeff. spectral decompositions */
    int spectralFlag = 0;
    double SDgamma0 = 0, SDgamma1 = 0, SDgamma2 = 0, SDgamma3 = 0;

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
                fprintf(stderr, "error parsing option %s with argument %s\n",
                    long_options[option_index].name, LALoptarg);
                exit(1);
            }
        case 'h':      /* help */
            usage(argv[0]);
            exit(0);
        case 'U':
            if (global_fmt == default_fmt)
                global_fmt = default_fmt_U;
            global_virial = 1;
            break;
        case 'E':
            global_epsrel = atof(LALoptarg);
            if (global_epsrel < 0 || global_epsrel == 0 || global_epsrel > 1) {
                fprintf(stderr, "invalid value of TOV solver routine relative error\n");
                exit(1);
            }
            break;
        case 'f':      /* eos-file */
            global_eos = XLALSimNeutronStarEOSFromFile(LALoptarg);
            break;
        case 'n':      /* eos-name */
            global_eos = XLALSimNeutronStarEOSByName(LALoptarg);
            break;
         case 'd':
            global_rho_c = atof(LALoptarg);
            if (global_rho_c < 0) {
                fprintf(stderr, "invalid value of central rest-mass density\n");
                exit(1);
            }
            break;
        case 'e':
            global_e_c = atof(LALoptarg);
            if (global_e_c < 0) {
                fprintf(stderr, "invalid value of central energy density\n");
                exit(1);
            }
            break;
        case 'F':      /* format */
            global_fmt = LALoptarg;
            break;
        case 'N':      /* npts */
            global_npts = strtol(LALoptarg, NULL, 0);
            if (global_npts < 1) {
                fprintf(stderr, "invalid number of points\n");
                exit(1);
            }
            break;

            /* using a 1-piece polytrope */
        case 'P':
            polytropeFlag = 1;
            break;
        case 'G':
            Gamma = atof(LALoptarg);
            break;
        case 'p':
            reference_pressure_si = atof(LALoptarg);
            break;
        case 'r':
            reference_density_si = atof(LALoptarg);
            break;

            /* using a 4-piece polytrope */
        case 'Q':
            piecewisePolytropeFlag = 1;
            break;
        case 'q':
            logp1_si = atof(LALoptarg);
            break;
        case '1':
            gamma1 = atof(LALoptarg);
            break;
        case '2':
            gamma2 = atof(LALoptarg);
            break;
        case '3':
            gamma3 = atof(LALoptarg);
            break;

           /* using 4-coeff. spectral decomposition */
        case 'S':
            spectralFlag = 1;
            break;
        case 'w':
            SDgamma0 = atof(LALoptarg);
            break;
        case 'x':
            SDgamma1 = atof(LALoptarg);
            break;
        case 'y':
            SDgamma2 = atof(LALoptarg);
            break;
        case 'z':
            SDgamma3 = atof(LALoptarg);
            break;

        default:
            fprintf(stderr, "unknown error while parsing options\n");
            exit(1);
        }
    }

    /* set eos to 1-piece polytrope */
    if (polytropeFlag == 1)
        global_eos =
            XLALSimNeutronStarEOSPolytrope(Gamma, reference_pressure_si,
            reference_density_si);

    /* set eos to 4-parameter piecewise polytrope */
    if (piecewisePolytropeFlag == 1)
        global_eos =
            XLALSimNeutronStarEOS4ParameterPiecewisePolytrope(logp1_si,
            gamma1, gamma2, gamma3);

    /* set eos to 4-coeff. spectral decomposition */
    if (spectralFlag == 1)
        global_eos =
            XLALSimNeutronStarEOS4ParameterSpectralDecomposition(SDgamma0,
            SDgamma1, SDgamma2, SDgamma3);

    if (LALoptind < argc) {
        fprintf(stderr, "extraneous command line arguments:\n");
        while (LALoptind < argc)
            fprintf(stderr, "%s\n", argv[LALoptind++]);
        exit(1);
    }

    if (!global_eos) {
        fprintf(stderr, "error: no equation of state selected\n");
        usage(argv[0]);
        exit(1);
    }

    return 0;
}

int usage(const char *program)
{
    int e, c;
    fprintf(stderr, "usage: %s [options]\n", program);
    fprintf(stderr,
        "\t-h, --help                   \tprint this message and exit\n");
    fprintf(stderr,
        "\t-f FILE, --eos-file=FILE     \tuse EOS with from data filename FILE\n");
    fprintf(stderr,
        "\t-n NAME, --eos-name=NAME     \tuse EOS with name NAME\n");
    fprintf(stderr, "\n");
    fprintf(stderr,
        "\t-U, --calculate-virial-theorem      \tcalculate and print Virial theorem besides M, R\n");
    fprintf(stderr, "\n");
    fprintf(stderr,
        "\t-E, --choose-relative-error      \tchoose a value for the TOV solver routine relative error (default epsrel=1e-6)\n");
    fprintf(stderr, "\n");
    fprintf(stderr,
        "\tYou may choose a central value for the rest-mass density (in kg/m^3) or the energy density (in J/m^3) as follows\n");
    fprintf(stderr,
        "\t-d, --central-rest-mass-density          \tcentral rest-mass density\n");
    fprintf(stderr,
        "\t-e, --central-energy-density             \tcentral energy density\n");
    fprintf(stderr,
        "\tIf none are chosen the the program will use a standard fixed value\n");
    fprintf(stderr, "\n");
    fprintf(stderr,
        "\t-P, --polytrope                  \tuse single polytrope\n");
    fprintf(stderr, "\t-G GAMMA, --gamma=GAMMA          \tadiabatic index\n");
    fprintf(stderr,
        "\t-p PRESSURE, --pressure=PRESSURE \tpressure at reference density\n");
    fprintf(stderr,
        "\t-r DENSITY, --density=DENSITY    \treference density\n");
    fprintf(stderr, "\n");
    fprintf(stderr,
        "\t-Q, --piecewisepolytrope         \tuse 4-parameter piecewise polytrope (PRD 79, 124032 (2009))\n");
    fprintf(stderr,
        "\t-q log(p_1), --logp1=log(p_1)    \tlog of pressure at rho_1=10^17.7 kg/m^3\n");
    fprintf(stderr,
        "\t-1 Gamma_1, --gamma1=Gamma_1     \tadiabatic index <10^17.7 kg/m^3\n");
    fprintf(stderr,
        "\t-2 Gamma_2, --gamma2=Gamma_2     \tadiabatic index 10^17.7--10^18 kg/m^3\n");
    fprintf(stderr,
        "\t-3 Gamma_3, --gamma3=Gamma_3     \tadiabatic index >10^18.0 kg/m^3\n");
    fprintf(stderr, "\n");
    fprintf(stderr,
        "\t-S --4paramspectraldecomp        \tuse 4-parameter spectral decomposition (PRD 82, 103011 (2010))\n");
    fprintf(stderr,
        "\t-w SDgamma0, --SDgamma0=SDgamma0 \tadiabatic index spectral decomposition coefficient 1 0.2--2.0\n");
    fprintf(stderr,
        "\t-x SDgamma1, --SDgamma1=SDgamma1 \tadiabatic index spectral decomposition coefficient 2 -1.6--1.7\n");
    fprintf(stderr,
        "\t-y SDgamma2, --SDgamma2=SDgamma2 \tadiabatic index spectral decomposition coefficient 3 -0.8--0.6\n");
    fprintf(stderr,
        "\t-z SDgamma3, --SDgamma3=SDgamma3 \tadiabatic index spectral decomposition coefficient 4 -0.2--0.2\n");
    fprintf(stderr, "\n");
    fprintf(stderr,
        "\t-N NPTS, --npts=NPTS         \toutput NPTS points [%ld]\n",
        default_npts);
    fprintf(stderr,
        "\t-F FORMAT, --format=FORMAT   \toutput format FORMAT [\"%s\"]\n",
        default_fmt);
    fprintf(stderr, "format string conversions:\n");
    fprintf(stderr,
        "\t%%c\t is replaced by the dimensionless compactness M/R\n");
    fprintf(stderr,
        "\t%%l\t is replaced by the dimensionless tidal parameter\n");
    fprintf(stderr, "\t%%k\t is replaced by the tidal love number k2\n");
    fprintf(stderr, "\t%%m\t is replaced by the mass in solar masses\n");
    fprintf(stderr, "\t%%r\t is replaced by the radius in kilometers\n");
    fprintf(stderr, "recognized tabulated equations of state:");
	for (e = 0, c = 0; e < (int)(sizeof(lalSimNeutronStarEOSNames)/sizeof(*lalSimNeutronStarEOSNames)); ++e) {
		c += fprintf(stderr, "%s%s", c ? ", " : "\n\t", lalSimNeutronStarEOSNames[e]);
		if (c > 50)
			c = 0;
	}
    fprintf(stderr, "\n");
    return 0;
}
