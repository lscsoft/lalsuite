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

#include <stdio.h>
#include <stdlib.h>

#include <lal/LALConstants.h>
#include <lal/LALgetopt.h>
#include <lal/LALStdlib.h>
#include <lal/LALSimNeutronStar.h>

int usage(const char *program);
int parseargs(int argc, char *argv[]);
int output(const char *fmt, double h, double p, double epsilon, double rho,
    double dedp, double vsound);

/* global variables  */
int global_cgs = 0;
int global_geom = 0;
LALSimNeutronStarEOS *global_eos = NULL;

const char *const default_fmt = "%p\\t%r\\n";
const char *global_fmt;

const long default_npts = 100;
long global_npts;

/* unit conversion constants */
#define CENTIMETER_SI 0.01      /* 1 cm in m */
#define DYNE_PER_CENTIMETER_SQUARED_SI 0.1      /* 1 dyn/cm^2 in Pa */
#define GRAM_PER_CENTIMETER_CUBED_SI 1000.0     /* 1 g/cm^3 in kg/m^3 */

int main(int argc, char *argv[])
{
    double hmax;
    long i;

    XLALSetErrorHandler(XLALAbortErrorHandler);

    global_npts = default_npts;
    global_fmt = default_fmt;
    parseargs(argc, argv);

    hmax = XLALSimNeutronStarEOSMaxPseudoEnthalpy(global_eos);
    for (i = 0; i < global_npts; i++) {
        double h = hmax * (0.5 + i) / global_npts;
        double p;
        double epsilon;
        double rho;
        double dedp;
        double vsound;
        if (global_geom) {      /* output in geometric units */
            p = XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(h,
                global_eos);
            epsilon =
                XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized
                (h, global_eos);
            rho =
                XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometerized
                (h, global_eos);
            dedp =
                XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometerized
                (p, global_eos);
            vsound =
                XLALSimNeutronStarEOSSpeedOfSoundGeometerized(h, global_eos);
        } else {        /* output in SI units by default */
            p = XLALSimNeutronStarEOSPressureOfPseudoEnthalpy(h, global_eos);
            epsilon =
                XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpy(h,
                global_eos);
            rho =
                XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpy(h,
                global_eos);
            dedp =
                XLALSimNeutronStarEOSEnergyDensityDerivOfPressure(p,
                global_eos);
            vsound = XLALSimNeutronStarEOSSpeedOfSound(h, global_eos);
            if (global_cgs) {   /* output in cgs units by converting from SI units */
                p /= DYNE_PER_CENTIMETER_SQUARED_SI;
                epsilon /= DYNE_PER_CENTIMETER_SQUARED_SI;      /* same as ERG_PER_CENTIMETER_CUBED_SI */
                rho /= GRAM_PER_CENTIMETER_CUBED_SI;
                /* dedp is dimensionless */
                vsound /= CENTIMETER_SI;
            }
        }
        output(global_fmt, h, p, epsilon, rho, dedp, vsound);
    }

    XLALDestroySimNeutronStarEOS(global_eos);
    LALCheckMemoryLeaks();
    return 0;
}

int output(const char *fmt, double h, double p, double epsilon, double rho,
    double dedp, double vsound)
{
    int i;
    for (i = 0; fmt[i]; ++i) {
        switch (fmt[i]) {
        case '%':      /* conversion character */
            switch (fmt[i + 1]) {
            case '%':
                fputc('%', stdout);
                break;
            case 'h':
                fprintf(stdout, "%e", h);
                break;
            case 'p':
                fprintf(stdout, "%e", p);
                break;
            case 'e':
                fprintf(stdout, "%e", epsilon);
                break;
            case 'r':
                fprintf(stdout, "%e", rho);
                break;
            case 'd':
                fprintf(stdout, "%e", dedp);
                break;
            case 'v':
                fprintf(stdout, "%e", vsound);
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
        {"cgs", no_argument, 0, 'c'},
        {"geom", no_argument, 0, 'g'},
        {"eos-file", required_argument, 0, 'f'},
        {"eos-name", required_argument, 0, 'n'},
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
        {0, 0, 0, 0}
    };
    char args[] = "hcgf:n:F:N:PG:p:r:Qq:1:2:3:";

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
        case 'c':      /* cgs */
            global_cgs = 1;
            break;
        case 'g':      /* geom */
            global_geom = 1;
            break;
        case 'f':      /* eos-file */
            global_eos = XLALSimNeutronStarEOSFromFile(LALoptarg);
            break;
        case 'n':      /* eos-name */
            global_eos = XLALSimNeutronStarEOSByName(LALoptarg);
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
        "\t-h, --help                       \tprint this message and exit\n");
    fprintf(stderr,
        "\t-f FILE, --eos-file=FILE         \tuse EOS with from data filename FILE\n");
    fprintf(stderr,
        "\t-n NAME, --eos-name=NAME         \tuse EOS with name NAME\n");
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
    fprintf(stderr, "\t-c, --cgs                        \tuse CGS units\n");
    fprintf(stderr,
        "\t-g, --geom                       \tuse geometerized units of length (G=c=1)\n");
    fprintf(stderr,
        "\t-N NPTS, --npts=NPTS             \toutput NPTS points [%ld]\n",
        default_npts);
    fprintf(stderr,
        "\t-F FORMAT, --format=FORMAT       \toutput format FORMAT [\"%s\"]\n",
        default_fmt);
    fprintf(stderr, "Format string conversions:\n");
    fprintf(stderr, "\t%%h\t is replaced by the pseudo-enthalpy\n");
    fprintf(stderr, "\t%%p\t is replaced by the pressure\n");
    fprintf(stderr, "\t%%e\t is replaced by the energy density\n");
    fprintf(stderr, "\t%%r\t is replaced by the rest-mass density\n");
    fprintf(stderr, "\t%%v\t is replaced by the speed of sound\n");

    return 0;
}
