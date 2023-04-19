#include <lal/LALStdlib.h>
#include <lal/LALDict.h>

/* UNITS ASSOCIATED WITH PARAMETERS */

static struct {const char *name; const char *unit;} lalSimInspiralREAL8WaveformParams[] = {
        /* length quantities */
        {"distance",                 "m"},
        /* mass quantities */
        {"mass1",             "kg"},
        {"mass2",             "kg"},
        {"total_mass",        "kg"},
        {"chirp_mass",        "kg"},
        {"mass_difference",   "kg"},
        {"reduced_mass",      "kg"},
        /* time quantities */
        {"deltaT",                "s"},
        /* frequency quantities */
        {"f22_start",             "Hz"},
        {"f22_ref",          "Hz"},
        {"f_max",                 "Hz"},
        {"deltaF",                "Hz"},
        {"f_ecc",                 "Hz"},
        /* angle quantities */
        {"inclination",           "rad"},
        {"phi_ref",             "rad"},
        {"longAscNodes",          "rad"},
        {"meanPerAno",            "rad"},
        {"spin1_tilt",            "rad"},
        {"spin1_phi",             "rad"},
        {"spin2_tilt",            "rad"},
        {"spin2_phi",             "rad"},
        /* adimensional quantities */
	{"sym_mass_ratio",           ""},
        {"mass_ratio",               ""},
        {"spin1_norm",               ""},
        {"spin2_norm",               ""},
        {"spin1x",                   ""},
        {"spin1y",                   ""},
        {"spin1z",                   ""},
        {"spin2x",                   ""},
        {"spin2y",                   ""},
        {"spin2z",                   ""},
        {"eccentricity",             ""},
	{"lambda1", 		     ""},
	{"lambda2", 		     ""},
	{"dQuadMon1",		     ""},
	{"dQuadMon2",		     ""},
    };
