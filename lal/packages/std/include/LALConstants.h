/*----------------------------------------------------------------------- 
 * 
 * File Name: LALConstants.h
 * 
 * Author: Creighton, T. D.
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * LALConstants.h
 * 
 * SYNOPSIS 
 * #include "LALConstants.h" 
 *
 * DESCRIPTION
 * A number of useful constants for LAL.  Computational constants are
 * taken from the IEEE standard 754 for binary arithmetic.
 * Mathematical constants are from the GNU C math.h header file.
 * Physical constants, astrophysical parameters, and cosmological
 * parameters are taken from the paper by Particle Data Group,
 * R. M. Barnett et al., Phys. Rev. D v.54 p.1, 1996
 * 
 * DIAGNOSTICS 
 *
 * CALLS
 * 
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _LALCONSTANTS_H
#define _LALCONSTANTS_H

#include "LALRCSID.h"

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (LALCONSTANTSH, "$Id$");

/* Computational constants (dimensionless) */

#define LAL_REAL4_MANT 24 /* Bits of precision in the mantissa of a REAL4 */
#define LAL_REAL4_MAX 3.40282347e+38 /* Largest REAL4 */
#define LAL_REAL4_MIN 1.17549435e-38 /* Smallest nonzero REAL4 */
#define LAL_REAL4_EPS 1.19209290e-07 /* 0.5^(LAL_REAL4_MANT-1) */
/* I.e. the difference between 1 and the next resolveable REAL4 */
#define LAL_REAL8_MANT 53 /* Bits of precision in the mantissa of a REAL8 */
#define LAL_REAL8_MAX 1.7976931348623157e+308 /* Largest REAL8 */
#define LAL_REAL8_MIN 2.2250738585072014e-308 /* Smallest nonzero REAL8 */
#define LAL_REAL8_EPS 2.2204460492503131e-16  /* 0.5^(LAL_REAL8_MANT-1) */
/* I.e. the difference between 1 and the next resolveable REAL8 */


/* Mathematical constants (dimensionless) */

#define LAL_E         2.7182818284590452353602874713526625L  /* e */
#define LAL_LOG2E     1.4426950408889634073599246810018922L  /* log_2 e */
#define LAL_LOG10E    0.4342944819032518276511289189166051L  /* log_10 e */
#define LAL_LN2       0.6931471805599453094172321214581766L  /* log_e 2 */
#define LAL_LN10      2.3025850929940456840179914546843642L  /* log_e 10 */
#define LAL_SQRT2     1.4142135623730950488016887242096981L  /* sqrt(2) */
#define LAL_SQRT1_2   0.7071067811865475244008443621048490L  /* 1/sqrt(2) */
#define LAL_GAMMA     0.5772156649015328606065120900824024L  /* gamma */
/* Assuming we're not near a black hole or in Tennessee... */
#define LAL_PI        3.1415926535897932384626433832795029L  /* pi */
#define LAL_TWOPI     6.2831853071795864769252867665590058L  /* 2*pi */
#define LAL_PI_2      1.5707963267948966192313216916397514L  /* pi/2 */
#define LAL_PI_4      0.7853981633974483096156608458198757L  /* pi/4 */
#define LAL_1_PI      0.3183098861837906715377675267450287L  /* 1/pi */
#define LAL_2_PI      0.6366197723675813430755350534900574L  /* 2/pi */
#define LAL_2_SQRTPI  1.1283791670955125738961589031215452L  /* 2/sqrt(pi) */

/* Physical constants, defined (SI) */

#define LAL_C_SI      299792458 /* Speed of light in vacuo, m s^-1 */
#define LAL_EPSILON0_SI  8.8541878176203898505365630317107503e-12
/* Permittivity of free space, C^2 N^-1 m^-2 */
#define LAL_MU0_SI    1.2566370614359172953850573533118012e-6
/* Permeability of free space, N A^-2 */
#define LAL_GEARTH_SI 9.80665 /* Standard gravity, m s^-2 */


/* Physical constants, measured (SI or dimensionless) */

#define LAL_G_SI      6.67259e-11    /* Gravitational constant, N m^2 kg^-2 */
#define LAL_H_SI      6.6260755e-34  /* Planck constant, J s */
#define LAL_HBAR_SI   1.05457266e-34 /* Reduced Planck constant, J s */
#define LAL_MPL_SI    2.17671e-8     /* Planck mass, kg */
#define LAL_LPL_SI    1.61605e-35    /* Planck length, m */
#define LAL_TPL_SI    5.39056e-44    /* Planck time, s */
#define LAL_K_SI      1.380658e-23   /* Boltzmann constant, J K^-1 */
#define LAL_R_SI      8.314511       /* Ideal gas constant, J K^-1 */
#define LAL_MOL       6.0221367e23   /* Avogadro constant, dimensionless */
#define LAL_BWIEN_SI  2.897756e-3    /* Wien displacement law constant, m K */
#define LAL_SIGMA_SI  5.67051e-8  /* Stefan-Boltzmann constant, W m^-2 K^-4 */
#define LAL_AMU_SI    1.6605402e-27  /* Atomic mass unit, kg */
#define LAL_MP_SI     1.6726231e-27  /* Proton mass, kg */
#define LAL_ME_SI     9.1093897e-31  /* Electron mass, kg */
#define LAL_QE_SI     1.60217733e-19 /* Electron charge, C */
#define LAL_ALPHA  7.297354677e-3 /* Fine structure constant, dimensionless */
#define LAL_RE_SI     2.81794092e-15 /* Classical electron radius, m */
#define LAL_LAMBDAE_SI 3.86159323e-13 /* Electron Compton wavelength, m */
#define LAL_AB_SI     5.29177249e-11 /* Bohr radius, m */
#define LAL_MUB_SI    9.27401543e-24 /* Bohr magneton, J T^-1 */
#define LAL_MUN_SI    5.05078658e-27 /* Nuclear magneton, J T^-1 */


/* Astrophysical parameters (SI) */

#define LAL_REARTH_SI 6.378140e6      /* Earth equatorial radius, m */
#define LAL_MEARTH_SI 5.97370e24      /* Earth mass, kg */
#define LAL_RSUN_SI   6.960e8         /* Solar equatorial radius, m */
#define LAL_MSUN_SI   1.98892e30      /* Solar mass, kg */
#define LAL_MRSUN_SI  1.47662504e3    /* Geometrized solar mass, m */
#define LAL_MTSUN_SI  4.92549095e-6   /* Geometrized solar mass, s */
#define LAL_LSUN_SI   3.846e26        /* Solar luminosity, W */
#define LAL_AU_SI     1.4959787066e11 /* Astronomical unit, m */
#define LAL_PC_SI     3.0856775807e16 /* Parsec, m */
#define LAL_YRTROP_SI 31556925.2      /* Tropical year (1994), s */
#define LAL_YRSID_SI  31558149.8      /* Sidereal year (1994), s */
#define LAL_DAYSID_SI 86164.09053     /* Mean sidereal day, s */
#define LAL_LYR_SI    9.46052817e15   /* ``Tropical'' lightyear (1994), m */


/* Cosmological parameters (SI) */

#define LAL_H0FAC_SI  3.2407792903e-18 /* Hubble constant prefactor, s^-1 */
#define LAL_H0_SI     2e-18            /* Approximate Hubble constant, s^-1 */
/* Hubble constant H0 = h0*HOFAC, where h0 is around 0.65 */
#define LAL_RHOCFAC_SI 1.68860e-9   /* Critical density prefactor, J m^-3 */
#define LAL_RHOC_SI   7e-10         /* Approximate critical density, J m^-3 */
/* Critical density RHOC = h0*h0*RHOCFAC, where h0 is around 0.65 */
#define LAL_TCBR_SI   2.726   /* Cosmic background radiation temperature, K */
#define LAL_VCBR_SI   3.695e5 /* Solar velocity with respect to CBR, m s^-1 */
#define LAL_RHOCBR_SI 4.177e-14 /* Energy density of CBR, J m^-3 */
#define LAL_NCBR_SI   4.109e8   /* Number density of CBR photons, m^-3 */
#define LAL_SCBR_SI   3.993e-14 /* Entropy density of CBR, J K^-1 m^-3 */


#ifdef  __cplusplus
}
#endif

#endif /* _LALCONSTANTS_H */
