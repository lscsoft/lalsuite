#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/LALSimNeutronStar.h>

int usage(const char *program);
int parseargs(int argc, char *argv[]);
int output(const char *fmt, double c, double m, double r, double k2);

/* global variables  */

LALSimNeutronStarEOS *global_eos = NULL;

const char *const default_fmt = "%r\\t%m\\n";
const char *global_fmt;

const long default_npts = 100;
long global_npts;

int main(int argc, char *argv[])
{
    const double logpmin = 75.5;
    double logpmax;
    double dlogp;
    long i;

    XLALSetErrorHandler(XLALAbortErrorHandler);

    global_npts = default_npts;
    global_fmt = default_fmt;
    parseargs(argc, argv);

    logpmax = log(XLALSimNeutronStarEOSMaxPressure(global_eos));
    /* THIS ONLY GIVES global_npts-1 POINTS: */
    dlogp = (logpmax - logpmin) / global_npts;

    for (i = 1; i < global_npts; ++i) { /* is this npts-1 or npts points? */
        double pc = exp(logpmin + i * dlogp);
        double c, m, r, k2;
        XLALSimNeutronStarTOVODEIntegrate(&r, &m, &k2, pc, global_eos);
        /* convert units */
        m /= LAL_MSUN_SI;       /* mass in solar masses */
        c = m * LAL_MRSUN_SI / r;       /* compactness (dimensionless) */
        r /= 1000.0;    /* radius in km */
        output(global_fmt, c, m, r, k2);
    }

    XLALDestroySimNeutronStarEOS(global_eos);
    LALCheckMemoryLeaks();
    return 0;
}

int output(const char *fmt, double c, double m, double r, double k2)
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
                fprintf(stdout, "%e", c);
                break;
            case 'k':
                fprintf(stdout, "%e", k2);
                break;
            case 'l':
                fprintf(stdout, "%e", (2.0 / 3.0) * k2 / pow(c, 5));
                break;
            case 'm':
                fprintf(stdout, "%e", m);
                break;
            case 'r':
                fprintf(stdout, "%e", r);
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
    struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
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
    char args[] = "hf:n:F:N:PG:p:r:Qq:1:2:3:";

    /* quantities for 1-piece polytrope: */
    int polytropeFlag = 0;
    double Gamma = 0, reference_pressure_si = 0, reference_density_si = 0;

    /* quantities for 4-parameter piecewise polytrope: */
    int piecewisePolytropeFlag = 0;
    double logp1_si = 0, gamma1 = 0, gamma2 = 0, gamma3 = 0;

    while (1) {
        int option_index = 0;
        int c;

        c = getopt_long_only(argc, argv, args, long_options, &option_index);
        if (c == -1)    /* end of options */
            break;

        switch (c) {
        case 0:        /* if option set a flag, nothing else to do */
            if (long_options[option_index].flag)
                break;
            else {
                fprintf(stderr, "error parsing option %s with argument %s\n",
                    long_options[option_index].name, optarg);
                exit(1);
            }
        case 'h':      /* help */
            usage(argv[0]);
            exit(0);
        case 'f':      /* eos-file */
            global_eos = XLALSimNeutronStarEOSFromFile(optarg);
            break;
        case 'n':      /* eos-name */
            global_eos = XLALSimNeutronStarEOSByName(optarg);
            break;
        case 'F':      /* format */
            global_fmt = optarg;
            break;
        case 'N':      /* npts */
            global_npts = strtol(optarg, NULL, 0);
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
            Gamma = atof(optarg);
            break;
        case 'p':
            reference_pressure_si = atof(optarg);
            break;
        case 'r':
            reference_density_si = atof(optarg);
            break;

            /* using a 4-piece polytrope */
        case 'Q':
            piecewisePolytropeFlag = 1;
            break;
        case 'q':
            logp1_si = atof(optarg);
            break;
        case '1':
            gamma1 = atof(optarg);
            break;
        case '2':
            gamma2 = atof(optarg);
            break;
        case '3':
            gamma3 = atof(optarg);
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

    if (optind < argc) {
        fprintf(stderr, "extraneous command line arguments:\n");
        while (optind < argc)
            fprintf(stderr, "%s\n", argv[optind++]);
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
    fprintf(stderr, "\n");
    fprintf(stderr,
        "\t-N NPTS, --npts=NPTS         \toutput NPTS points [%ld]\n",
        default_npts);
    fprintf(stderr,
        "\t-F FORMAT, --format=FORMAT   \toutput format FORMAT [\"%s\"]\n",
        default_fmt);
    fprintf(stderr, "Format string conversions:\n");
    fprintf(stderr,
        "\t%%c\t is replaced by the dimensionless compactness M/R\n");
    fprintf(stderr,
        "\t%%l\t is replaced by the dimensionless tidal parameter\n");
    fprintf(stderr, "\t%%k\t is replaced by the tidal love number k2\n");
    fprintf(stderr, "\t%%m\t is replaced by the mass in solar masses\n");
    fprintf(stderr, "\t%%r\t is replaced by the radius in kilometers\n");
    return 0;
}
