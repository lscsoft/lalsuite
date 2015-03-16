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
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
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
int output(const char *fmt, double c, double m, double r, double k2);

/* global variables  */

LALSimNeutronStarEOS *global_eos = NULL;
double global_mass;

int main(int argc, char *argv[])
{
    LALSimNeutronStarFamily *fam;
    XLALSetErrorHandler(XLALAbortErrorHandler);

    parseargs(argc, argv);

    fam = XLALCreateSimNeutronStarFamily(global_eos);
    printf("Equation of State: %s\n", XLALSimNeutronStarEOSName(global_eos));
    printf("Maximum Mass (solar) = %g\n",
        XLALSimNeutronStarMaximumMass(fam) / LAL_MSUN_SI);
    if (global_mass != 0.0) {
        double m = global_mass * LAL_MSUN_SI;
        double p = XLALSimNeutronStarCentralPressure(m, fam);
        double r = XLALSimNeutronStarRadius(m, fam);
        double k = XLALSimNeutronStarLoveNumberK2(m, fam);
        double c = global_mass * LAL_MRSUN_SI / r;
        double l = (2.0 / 3.0) * k / pow(c, 5);
        printf("Parameters for Neutron Star Mass (solar) = %g\n", global_mass);
        printf("- Central Pressure (Pa) = %g\n", p);
        printf("- Radius (km) = %g\n", r / 1000.0);
        printf("- Compactness (dimensionless) = %g\n", c);
        printf("- Love Number (dimensionless) = %g\n", k);
        printf("- Tidal Parameter (dimensionless) = %g\n", l);
    }

    XLALDestroySimNeutronStarFamily(fam);
    XLALDestroySimNeutronStarEOS(global_eos);
    LALCheckMemoryLeaks();
    return 0;
}

int parseargs(int argc, char **argv)
{
    struct LALoption long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"mass", required_argument, 0, 'm'},
        {"eos-file", required_argument, 0, 'f'},
        {"eos-name", required_argument, 0, 'n'},
        {"polytrope", no_argument, 0, 'P'},
        {"gamma", required_argument, 0, 'G'},
        {"pressure", required_argument, 0, 'p'},
        {"density", required_argument, 0, 'r'},
        {"piecewisepolytrope", no_argument, 0, 'Q'},
        {"logp1", required_argument, 0, 'q'},
        {"gamma1", required_argument, 0, '1'},
        {"gamma2", required_argument, 0, '2'},
        {"gamma3", required_argument, 0, '3'},
        {0, 0, 0, 0}
    };
    char args[] = "hm:f:n:PG:p:r:Qq:1:2:3:";

    /* quantities for 1-piece polytrope: */
    int polytropeFlag = 0;
    double Gamma = 0, reference_pressure_si = 0, reference_density_si = 0;

    /* quantities for 4-parameter piecewise polytrope: */
    int piecewisePolytropeFlag = 0;
    double logp1_si = 0, gamma1 = 0, gamma2 = 0, gamma3 = 0;

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
        case 'm':      /* mass */
            global_mass = atof(LALoptarg);
            break;
        case 'f':      /* eos-file */
            global_eos = XLALSimNeutronStarEOSFromFile(LALoptarg);
            break;
        case 'n':      /* eos-name */
            global_eos = XLALSimNeutronStarEOSByName(LALoptarg);
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
    fprintf(stderr, "usage: %s [options]\n", program);
    fprintf(stderr,
        "\t-h, --help                   \tprint this message and exit\n");
    fprintf(stderr,
        "\t-m MASS, --mass=MASS         \tparameters of a neutron star of mass MASS (solar)\n");
    fprintf(stderr,
        "\t-f FILE, --eos-file=FILE     \tuse EOS with from data filename FILE\n");
    fprintf(stderr,
        "\t-n NAME, --eos-name=NAME     \tuse EOS with name NAME\n");
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
    return 0;
}
