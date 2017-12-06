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
 * \addtogroup LALConstants_h
 * \author Creighton, T. D.
 * \brief Provides standard numerical constants for LAL.
 *
 * This header defines a number of useful numerical constants
 * for use in LAL routines.  These constants come in three basic
 * flavours: arithmetic and mathematical constants, fundamental (or
 * defined) physical constants, and measured astrophysical and
 * cosmological parameters.
 *
 * Note that this header is not included automatically by the header
 * <tt>LALStdlib.h</tt>.  Include it explicitly if you need any of these
 * constants.
 */
/*@{*/

#ifndef _LALCONSTANTS_H
#define _LALCONSTANTS_H

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \name Floating-point constants
 * The following constants define the precision and range of
 * floating-point arithmetic in LAL.  They are taken from the IEEE
 * standard 754 for binary arithmetic.  All numbers are dimensionless.
 */
/*@{*/
#define LAL_REAL4_MANT 24 /**< Bits of precision in the mantissa of a REAL4 */
#define LAL_REAL4_MAX 3.40282347e+38 /**< Largest REAL4 */
#define LAL_REAL4_MIN 1.17549435e-38 /**< Smallest nonzero REAL4 */
#define LAL_REAL4_EPS 1.19209290e-07 /**< 0.5^(LAL_REAL4_MANT-1), ie the difference between 1 and the next resolveable REAL4 */
#define LAL_REAL8_MANT 53 /**< Bits of precision in the mantissa of a REAL8 */
#define LAL_REAL8_MAX 1.7976931348623157e+308 /**< Largest REAL8 */
#define LAL_REAL8_MIN 2.2250738585072014e-308 /**< Smallest nonzero REAL8 */
#define LAL_REAL8_EPS 2.2204460492503131e-16  /**< 0.5^(LAL_REAL8_MANT-1), ie the difference between 1 and the next resolveable REAL8 */
/*@}*/


/**
 * \name Integer constants
 * Extremal integer values, all expressed as unsigned long long.
 */
/*@{*/
#define LAL_UINT8_MAX   LAL_UINT8_C(18446744073709551615)
#define LAL_UINT4_MAX   LAL_UINT8_C(4294967295)
#define LAL_UINT2_MAX   LAL_UINT8_C(65535)
#define LAL_INT8_MAX    LAL_UINT8_C(9223372036854775807)
#define LAL_INT4_MAX    LAL_UINT8_C(2147483647)
#define LAL_INT2_MAX    LAL_UINT8_C(32767)
/*@}*/


/**
 * \name Mathematical constants
 * The following are fundamental mathematical constants.  They are mostly
 * taken from the GNU C <tt>math.h</tt> header (with the exception of
 * <tt>LAL_TWOPI</tt>, which was computed using Maple).  All numbers are
 * dimensionless. The value of exp(gamma) is taken from
 * http://www.research.att.com/~njas/sequences/A073004
 */
/*@{*/
#define LAL_E         2.7182818284590452353602874713526625  /**< e */
#define LAL_LOG2E     1.4426950408889634073599246810018922  /**< log_2 e */
#define LAL_LOG10E    0.4342944819032518276511289189166051  /**< log_10 e */
#define LAL_LN2       0.6931471805599453094172321214581766  /**< log_e 2 */
#define LAL_LN10      2.3025850929940456840179914546843642  /**< log_e 10 */
#define LAL_SQRT2     1.4142135623730950488016887242096981  /**< sqrt(2) */
#define LAL_SQRT1_2   0.7071067811865475244008443621048490  /**< 1/sqrt(2) */
#define LAL_GAMMA     0.5772156649015328606065120900824024  /**< gamma */
#define LAL_EXPGAMMA  1.7810724179901979852365041031071795  /**< exp(gamma) */
/* Assuming we're not near a black hole or in Tennessee... */
#define LAL_PI        3.1415926535897932384626433832795029  /**< pi */
#define LAL_TWOPI     6.2831853071795864769252867665590058  /**< 2*pi */
#define LAL_PI_2      1.5707963267948966192313216916397514  /**< pi/2 */
#define LAL_PI_4      0.7853981633974483096156608458198757  /**< pi/4 */
#define LAL_1_PI      0.3183098861837906715377675267450287  /**< 1/pi */
#define LAL_2_PI      0.6366197723675813430755350534900574  /**< 2/pi */
#define LAL_2_SQRTPI  1.1283791670955125738961589031215452  /**< 2/sqrt(pi) */
#define LAL_PI_180    1.7453292519943295769236907684886127e-2 /**< pi/180 */
#define LAL_180_PI    57.295779513082320876798154814105170 /**< 180/pi */
/*@}*/

/**
 * \name Exact physical constants
 * The following physical constants are defined to have exact values.
 * The values of \f$c\f$ and \f$g\f$ are taken from \cite Barnet_1996,
 * \f$p_\mathrm{atm}\f$ is from \cite Lang_1992, while \f$\epsilon_0\f$ and
 * \f$\mu_0\f$ are computed from \f$c\f$ using exact formulae.  The use
 * of a Julian year (365.25 days) as standard is specified by the IAU.
 * They are given in the SI units shown.
 */
/*@{*/
#define LAL_C_SI      299792458 /**< Speed of light in vacuo, m s^-1 */
#define LAL_EPSILON0_SI  8.8541878176203898505365630317107503e-12 /**< Permittivity of free space, C^2 N^-1 m^-2 */
#define LAL_MU0_SI    1.2566370614359172953850573533118012e-6 /**< Permeability of free space, N A^-2 */
#define LAL_GEARTH_SI 9.80665 /**< Standard gravity, m s^-2 */
#define LAL_PATM_SI 101325 /**< Standard atmosphere, Pa */
#define LAL_YRJUL_SI 31557600 /**< Julian year, s */
#define LAL_LYR_SI 9.4607304725808e15 /**< (Julian) Lightyear, m */
/*@}*/

/**
 * \name Physical constants
 * The following are measured fundamental physical constants, with values
 * given in \cite Barnet_1996.  When not dimensionless, they are given
 * in the SI units shown.
 */
/*@{*/
#define LAL_G_SI      6.67259e-11    /**< Gravitational constant, N m^2 kg^-2 */
#define LAL_H_SI      6.6260755e-34  /**< Planck constant, J s */
#define LAL_HBAR_SI   1.0545726691251019773669079307477023e-34 /**< Reduced Planck constant, J s.  = LAL_H_SI / LAL_TWOPI */
#define LAL_MPL_SI    2.1767140835297016797409334934257335e-8  /**< Planck mass, kg.  = sqrt{LAL_HBAR_SI * LAL_C_SI / LAL_G_SI} */
#define LAL_LPL_SI    1.6160486159348859434398861412879278e-35 /**< Planck length, m.  = sqrt{LAL_HBAR_SI * LAL_G_SI / LAL_C_SI^3} */
#define LAL_TPL_SI    5.3905579437054615411301846068720240e-44 /**< Planck time, s.  = sqrt{LAL_HBAR_SI * LAL_G_SI / LAL_C_SI^5} */
#define LAL_K_SI      1.380658e-23   /**< Boltzmann constant, J K^-1 */
#define LAL_R_SI      8.314511       /**< Ideal gas constant, J K^-1 */
#define LAL_MOL       6.0221367e23   /**< Avogadro constant, dimensionless */
#define LAL_BWIEN_SI  2.897756e-3    /**< Wien displacement law constant, m K */
#define LAL_SIGMA_SI  5.67051e-8  /**< Stefan-Boltzmann constant, W m^-2 K^-4 */
#define LAL_AMU_SI    1.6605402e-27  /**< Atomic mass unit, kg */
#define LAL_MP_SI     1.6726231e-27  /**< Proton mass, kg */
#define LAL_ME_SI     9.1093897e-31  /**< Electron mass, kg */
#define LAL_QE_SI     1.60217733e-19 /**< Electron charge, C */
#define LAL_ALPHA  7.297354677e-3 /**< Fine structure constant, dimensionless */
#define LAL_RE_SI     2.81794092e-15 /**< Classical electron radius, m */
#define LAL_LAMBDAE_SI 3.86159323e-13 /**< Electron Compton wavelength, m */
#define LAL_AB_SI     5.29177249e-11 /**< Bohr radius, m */
#define LAL_MUB_SI    9.27401543e-24 /**< Bohr magneton, J T^-1 */
#define LAL_MUN_SI    5.05078658e-27 /**< Nuclear magneton, J T^-1 */
/*@}*/

/**
 * \name Astrophysical parameters
 * The following parameters are derived from measured properties of the
 * Earth and Sun.  The values are taken from \cite Barnet_1996, except
 * for the obliquity of the ecliptic plane and the eccentricity of
 * Earth's orbit, which are taken from \cite Lang_1992.  All values are
 * given in the SI units shown.  Note that the ``year'' and
 * ``light-year'' have exactly defined values, and appear under
 * ``Exact physical constants''.
 */
/*@{*/
#define LAL_REARTH_SI 6.378140e6      /**< Earth equatorial radius, m */
#define LAL_AWGS84_SI 6.378137e6      /**< Semimajor axis of WGS-84 Reference Ellipsoid, m */
#define LAL_BWGS84_SI 6.356752314e6   /**< Semiminor axis of WGS-84 Reference Ellipsoid, m */
#define LAL_MEARTH_SI 5.97370e24      /**< Earth mass, kg */
#define LAL_IEARTH    0.409092804     /**< Earth inclination (2000), radians */
#define LAL_COSIEARTH 0.91748206215761919815    /**< Cosine of Earth inclination (2000) */
#define LAL_SINIEARTH 0.39777715572793088957    /**< Sine of Earth inclination (2000) */
#define LAL_EEARTH    0.0167          /**< Earth orbital eccentricity */
#define LAL_RSUN_SI   6.960e8         /**< Solar equatorial radius, m */
#define LAL_MSUN_SI   1.98892e30      /**< Solar mass, kg */
#define LAL_MRSUN_SI  1.4766254500421874513093320107664308e3  /**< Geometrized solar mass, m.  = LAL_MSUN_SI / LAL_MPL_SI * LAL_LPL_SI */
#define LAL_MTSUN_SI  4.9254923218988636432342917247829673e-6 /**< Geometrized solar mass, s.  = LAL_MSUN_SI / LAL_MPL_SI * LAL_TPL_SI */
#define LAL_LSUN_SI   3.846e26        /**< Solar luminosity, W */
#define LAL_AU_SI     1.4959787066e11 /**< Astronomical unit, m */
#define LAL_PC_SI     3.0856775807e16 /**< Parsec, m */
#define LAL_YRTROP_SI 31556925.2      /**< Tropical year (1994), s */
#define LAL_YRSID_SI  31558149.8      /**< Sidereal year (1994), s */
#define LAL_DAYSID_SI 86164.09053     /**< Mean sidereal day, s */
/*@}*/

/**
 * \name Cosmological parameters
 * The following cosmological parameters are derived from measurements of
 * the Hubble expansion rate and of the cosmic background radiation
 * (CBR).  Data are taken from \cite Barnet_1996.  In what follows, the
 * normalized Hubble constant \f$h_0\f$ is equal to the actual Hubble
 * constant \f$H_0\f$ divided by \f$\langle H
 * \rangle=100\,\mathrm{km}\,\mathrm{s}^{-1}\mathrm{Mpc}^{-1}\f$.  Thus the
 * Hubble constant can be written as:
 * \f$H_0 = \langle H \rangle h_0\f$.
 * Similarly, the critical energy density \f$\rho_c\f$ required for spatial
 * flatness is given by: \f$\rho_c = \langle\rho\rangle h_0^2\f$.
 * Current estimates give \f$h_0\f$ a value of around 0.65, which is what is
 * assumed below.  All values are in the SI units shown.
 */
/*@{*/
#define LAL_H0FAC_SI  3.2407792903e-18 /**< Hubble constant prefactor, s^-1 */
#define LAL_H0_SI     2e-18            /**< Approximate Hubble constant, s^-1 */
/* Hubble constant H0 = h0*HOFAC, where h0 is around 0.65 */
#define LAL_RHOCFAC_SI 1.68860e-9   /**< Critical density prefactor, J m^-3 */
#define LAL_RHOC_SI   7e-10         /**< Approximate critical density, J m^-3 */
/* Critical density RHOC = h0*h0*RHOCFAC, where h0 is around 0.65 */
#define LAL_TCBR_SI   2.726   /**< Cosmic background radiation temperature, K */
#define LAL_VCBR_SI   3.695e5 /**< Solar velocity with respect to CBR, m s^-1 */
#define LAL_RHOCBR_SI 4.177e-14 /**< Energy density of CBR, J m^-3 */
#define LAL_NCBR_SI   4.109e8   /**< Number density of CBR photons, m^-3 */
#define LAL_SCBR_SI   3.993e-14 /**< Entropy density of CBR, J K^-1 m^-3 */
/*@}*/

/*@}*/
#ifdef  __cplusplus
}
#endif
#endif /* _LALCONSTANTS_H */
