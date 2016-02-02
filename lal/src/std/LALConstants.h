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
 * @defgroup LALConstants_h Header LALConstants.h
 * @ingroup lal_std
 * @author Creighton, T. D.
 * @brief Provides standard numerical constants for LAL.
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
 *//*@{*/

/*@{*/
#ifndef _LALCONSTANTS_H
#define _LALCONSTANTS_H

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * @name Floating-point constants
 * The following constants define the precision and range of
 * floating-point arithmetic in LAL.  They are taken from the IEEE
 * standard 754 for binary arithmetic.  All numbers are dimensionless.
 * @see http://dx.doi.org/10.1109/IEEESTD.2008.4610935
 */
/*@{*/
#if __STDC_VERSION__ >= 199901L
#define LAL_REAL4_MANT 24 /**< Bits of precision in the mantissa of a REAL4 */
#define LAL_REAL4_MAX 0x1.fffffe0000000p+127 /**< Largest normalized REAL4 number (2-2^-23)*2^127 */
#define LAL_REAL4_MIN 0x1.0000000000000p-126 /**< Smallest normalized REAL4 number 2^-126 */
#define LAL_REAL4_EPS 0x1.0000000000000p-23 /**< Difference between 1 and the next resolvable REAL4 2^-23 */
#define LAL_REAL8_MANT 53 /**< Bits of precision in the mantissa of a REAL8 */
#define LAL_REAL8_MAX 0x1.fffffffffffffp+1023 /**< Largest normalized REAL8 number (2-2^-52)*2^1023 */
#define LAL_REAL8_MIN 0x1.0000000000000p-1022 /**< Smallest normalized REAL8 number 2^-1022 */
#define LAL_REAL8_EPS 0x1.0000000000000p-52 /**< Difference between 1 and the next resolvable REAL8 2^-52 */
#else
#define LAL_REAL4_MANT 24
#define LAL_REAL4_MAX 3.4028234663852886e+38
#define LAL_REAL4_MIN 1.1754943508222875e-38
#define LAL_REAL4_EPS 1.1920928955078125e-07
#define LAL_REAL8_MANT 53
#define LAL_REAL8_MAX 1.797693134862315708145274237317043568e308
#define LAL_REAL8_MIN 2.225073858507201383090232717332404064e-308
#define LAL_REAL8_EPS 2.220446049250313080847263336181640625e-16
#endif
/*@}*/

/**
 * @name Integer constants
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
 * @name Mathematical constants
 * All are dimensionless.
 *
 * @def LAL_E
 * @brief Euler's constant, e
 * @see https://oeis.org/A001113
 *
 * @def LAL_LOG2E
 * @brief base-2 logarithm of e, log_2(e)
 * @see http://oeis.org/A007525
 *
 * @def LAL_LOG10E
 * @brief common logarithm of e, log_10(e)
 * @see http://oeis.org/A002285
 *
 * @def LAL_LN2
 * @brief natural log of 2, ln(2)
 * @see http://oeis.org/A002162
 *
 * @def LAL_LN10
 * @brief natural log of 10, ln(10)
 * @see http://oeis.org/A002392
 *
 * @def LAL_SQRT2
 * @brief Pythagoras's constant, sqrt(2)
 * @see https://oeis.org/A002193
 *
 * @def LAL_SQRT1_2
 * @brief 1/sqrt(2)
 * @see https://oeis.org/A010503
 *
 * @def LAL_GAMMA
 * @brief Euler-Mascheroni constant, gamma
 * @see https://oeis.org/A001620
 *
 * @def LAL_EXPGAMMA
 * @brief exp(gamma)
 * @see https://oeis.org/A073004
 *
 * @def LAL_PI
 * @brief Archimedes's constant, pi
 * @see https://oeis.org/A000796
 *
 * @def LAL_TWOPI
 * @brief 2*pi is circumference of a circle divided by its radius
 * @see https://oeis.org/A019692
 *
 * @def LAL_PI_2
 * @brief pi/2
 * @see http://oeis.org/A019669
 *
 * @def LAL_PI_4
 * @brief pi/4 is the least positive solution to sin(x) = cos(x)
 * @see http://oeis.org/A003881
 *
 * @def LAL_1_PI
 * @brief 1/pi is the ratio of the volume of a regular octahedron
 * to the volume of the circumscribed sphere
 * @see http://oeis.org/A049541
 *
 * @def LAL_2_PI
 * @brief 2/pi is Buffon's constant
 * @see http://oeis.org/A060294
 *
 * @def LAL_2_SQRTPI
 * @brief 2/sqrt(pi)
 * @see http://oeis.org/A190732
 *
 * @def LAL_PI_180
 * @brief 180/pi is the number of degrees in one radian
 * @see http://oeis.org/A072097
 *
 * @def LAL_180_PI
 * @brief pi/180 is the number of radians in one degree
 * @see http://oeis.org/A019685
 *
 * @def LAL_LNPI
 * @brief natural log of pi, ln(pi)
 * @see http://oeis.org/A053510
 */
/*@{*/
#define LAL_E         2.718281828459045235360287471352662498
#define LAL_LOG2E     1.442695040888963407359924681001892137
#define LAL_LOG10E    0.434294481903251827651128918916605082
#define LAL_LN2       0.693147180559945309417232121458176568
#define LAL_LN10      2.302585092994045684017991454684364208
#define LAL_SQRT2     1.414213562373095048801688724209698079
#define LAL_SQRT1_2   0.707106781186547524400844362104849039
#define LAL_GAMMA     0.577215664901532860606512090082402431
#define LAL_EXPGAMMA  1.781072417990197985236504103107179549
/* Assuming we're not near a black hole or in Tennessee... */
#define LAL_PI        3.141592653589793238462643383279502884
#define LAL_TWOPI     6.283185307179586476925286766559005768
#define LAL_PI_2      1.570796326794896619231321691639751442
#define LAL_PI_4      0.785398163397448309615660845819875721
#define LAL_1_PI      0.318309886183790671537767526745028724
#define LAL_2_PI      0.636619772367581343075535053490057448
#define LAL_2_SQRTPI  1.128379167095512573896158903121545172
#define LAL_PI_180    1.745329251994329576923690768488612713e-2
#define LAL_180_PI    57.295779513082320876798154814105170332
#define LAL_LNPI      1.144729885849400174143427351353058712
/*@}*/

/**
 * @name Exact physical constants
 * The following physical constants are defined to have exact values.
 * The dimensions in SI units are as shown.
 * @see http://dx.doi.org/10.1103/RevModPhys.84.1527
 */
/*@{*/
#define LAL_C_SI 299792458e0 /**< Speed of light in vacuo, m s^-1 */
#define LAL_EPSILON0_SI 8.854187817620389850536563031710750261e-12 /**< Permittivity of free space, C^2 N^-1 m^-2 */
#define LAL_MU0_SI 1.256637061435917295385057353311801154e-6 /**< Permeability of free space, N A^-2 */
#define LAL_GEARTH_SI 9.80665 /**< Standard gravity, m s^-2 */
#define LAL_PATM_SI 101325e0 /**< Standard atmosphere, Pa */
/*@}*/

/**
 * @name Primary physical constants
 * These physical constants are given to the precision
 * to which they are known.  Other physical constants
 * derived from these are given in the next section.
 * @see http://dx.doi.org/10.1103/RevModPhys.84.1527
 */
/*@{*/
#define LAL_QE_SI 1.60217656535e-19 /**< Electron charge, C */
#define LAL_ALPHA 0.007297352569824 /**< Fine structure constant, dimensionless */
#define LAL_RYD_SI 10973731.56853955 /**< Rydberg constant, m^-1 */
#define LAL_MP_ME 1836.1526724575 /**< Proton-electron mass ratio, dimensionless */
#define LAL_MP_AMU 1.00727646681290 /**< Proton molar mass, kg mol^-1 */
#define LAL_R_SI 8.314462175 /**< Molar gas constant, J mol^-1 K^-1 */
#define LAL_G_SI 6.67384e-11 /**< Gravitational constant, N m^2 kg^-2 */
/*@}*/

/**
 * @name Derived physical constants
 * The following constants are derived from the primary
 * physical constants.  When not dimensionless, they are
 * given in the SI units shown.  Precision beyond the
 * accuracy is retained for these constants in order
 * that equivalent combinations yield the same value.
 */
/*@{*/

/**
 * @brief Reduced Planck constant, J s
 * @details
 * LAL_HBAR_SI = ((LAL_QE_SI * LAL_QE_SI) / (4 * LAL_PI * LAL_EPSILON0_SI * LAL_C_SI * LAL_ALPHA))
 */
#define LAL_HBAR_SI 1.054571726620627613731476059278130752e-34

/**
 * @brief Planck constant, J s
 * @details
 * LAL_H_SI = LAL_TWOPI * LAL_HBAR_SI
 */
#define LAL_H_SI 6.626069578069735006764896984166285474e-34

/**
 * @brief Planck mass, kg
 * @details
 * LAL_MPL_SI =	sqrt(LAL_HBAR_SI * LAL_C_SI / LAL_G_SI)
 */
#define LAL_MPL_SI 2.176509253770671995475098872904844844e-8

/**
 * @brief Planck length, m
 * @details
 * LAL_LPL_SI =	(LAL_HBAR_SI / (LAL_MPL_SI * LAL_C_SI))
 */
#define LAL_LPL_SI 1.616199256687500316071519963813388523e-35

/**
 * @brief Planck time, s
 * @details
 * LAL_TPL_SI =	(LAL_LPL_SI / LAL_C_SI)
 */
#define LAL_TPL_SI 5.391060427168919359777623104225619054e-44

/**
 * @brief Planck luminosity, J s^-1
 * @details
 * LAL_LUMPL_SI = (LAL_C_SI * LAL_C_SI * LAL_C_SI * LAL_C_SI * LAL_C_SI) / (LAL_G_SI)
 */
#define LAL_LUMPL_SI 3.628504984913064522721519179529402840e52

/**
 * @brief Electron mass, kg
 * @details
 * LAL_ME_SI = ((2 * LAL_RYD_SI * LAL_H_SI) / (LAL_C_SI * LAL_ALPHA * LAL_ALPHA))
 */
#define LAL_ME_SI 9.109382914246707433083507609082565768e-31

/**
 * @brief Proton mass, kg
 * @details
 * LAL_MP_SI = (LAL_ME_SI * LAL_MP_ME)
 */
#define LAL_MP_SI 1.672621778243278140372641234618498933e-27

/**
 * @brief Atomic mass uint, kg
 * @details
 * LAL_AMU_SI = (LAL_MP_SI / LAL_MP_AMU)
 */
#define LAL_AMU_SI 1.660538921886641249862170970849963401e-27

/**
 * @brief Bohr radius, m
 * @details
 * LAL_AB_SI = (LAL_ALPHA / (4 * LAL_PI * LAL_RYD_SI))
 */
#define LAL_AB_SI 5.291772109231583273539324387807930183e-11

/**
 * @brief Electron Compton wavelength, m
 * @details
 * LAL_LAMBDAE_SI = (2 * LAL_PI * LAL_ALPHA * LAL_AB_SI)
 */
#define LAL_LAMBDAE_SI 2.426310238942902554172796303895891665e-12

/**
 * @brief Classical electron radius, m
 * @details
 * LAL_RE_SI = (LAL_ALPHA * LAL_ALPHA * LAL_AB_SI)
 */
#define LAL_RE_SI 2.817940326717505398482365681324267402e-15

/**
 * @brief Bohr magneton, J T^-1
 * LAL_MUB_SI = (LAL_LAMBDAE_SI * LAL_C_SI * LAL_QE_SI / (4 * LAL_PI))
 * @details
 */
#define LAL_MUB_SI 9.274009681982817185200341325664841126e-24

/**
 * @brief Nuclear magneton, J T^-1
 * @details
 * LAL_MUN_SI =	(LAL_MUB_SI / LAL_MP_ME)
 */
#define LAL_MUN_SI 5.050783532924044193470852682533321512e-27

/**
 * @brief Avogadro constant, dimensionless
 * @details
 * LAL_MOL = (0.001 * LAL_MP_AMU / LAL_MP_SI)
 */
#define LAL_MOL 6.022141286901230737569270163033480013e23

/**
 * @brief Boltzmann constant, J K^-1
 * @details
 * LAL_K_SI = (LAL_R_SI / LAL_MOL)
 */
#define LAL_K_SI 1.380648805614175830977374450051504830e-23

/**
 * @brief Stefan-Boltzmann constant, W m^-2 K^-4
 * @details
 * LAL_SIGMA_SI = ((LAL_PI * LAL_PI * LAL_K_SI * LAL_K_SI * LAL_K_SI * LAL_K_SI) / (60 * LAL_HBAR_SI * LAL_HBAR_SI * LAL_HBAR_SI * LAL_C_SI * LAL_C_SI))
 */
#define LAL_SIGMA_SI 5.670372694104351560189523450147062255e-8

/**
 * @brief Second radiation constant, m K
 * @details
 * LAL_C2RAD_SI = (LAL_H_SI * LAL_C_SI / LAL_K_SI)
 */
#define LAL_C2RAD_SI 1.438776955885524212375359308606994213e-2

/**
 * @brief Wien displacement law constant, m K
 * @details
 * LAL_BWIEN_SI = (LAL_C2RAD_SI / X)
 *
 * where the factor X satisfies
 *
 * X * exp(X) = 5 * (exp(X) - 1)
 */
#define LAL_BWIEN_SI 2.897772113049799237638435257604006684e-3
/*@}*/

/**
 * @name Exact astrophysical parameters
 * The following astrophysical constants are defined to have exact values.
 * The dimensions in SI units are as shown.
 * @see http://asa.usno.navy.mil/static/files/2015/Astronomical_Constants_2015.pdf
 */
/*@{*/
#define LAL_SOL_SID 1.002737909350795 /**< Ratio of mean solar day to sidereal day, dimensionless */
#define LAL_DAYJUL_SI 86400e0 /**< Julian day, s */
#define LAL_YRJUL_SI 31557600e0 /**< Julian year, s */
#define LAL_LYR_SI 9460730472580800 /**< (Julian) Lightyear, m */
#define LAL_AU_SI 149597870700e0 /**< Astronomical unit, m */
#define LAL_PC_SI 3.085677581491367278913937957796471611e16 /**< Parsec, m */
/*@}*/

/**
 * @name Primary astrophysical parameters
 * These astrophysical constants are given to the precision
 * to which they are known.  Other physical constants
 * derived from these are given in the next section.
 */
/*@{*/

/**
 * @brief Earth equatorial radius, m
 * @see http://asa.usno.navy.mil/static/files/2015/Astronomical_Constants_2015.pdf
 */
#define LAL_REARTH_SI 6378136.6

/**
 * @brief Semimajor axis of WGS-84 Reference Ellipsoid, m
 * @see Department of Defense World Geodedic System 1984 http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf
 */
#define LAL_AWGS84_SI 6378137e0

/**
 * @brief Semiminor axis of WGS-84 Reference Ellipsoid, m
 * @details
 * LAL_BWGS84_SI = LAL_AWGS84_SI * (1 - f)
 *
 * where f is the flattening
 *
 * 1/f = 298.257223563
 *
 * @note This constant is not given to full precision for compatibility.
 * @see Department of Defense World Geodedic System 1984 http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf
 */
#define LAL_BWGS84_SI 6356752.314

/**
 * @brief Earth inclination (2000), radians
 * @details
 * This is the measured value of the mean obliquity of the
 * ecliptic, 84381.406 arcseconds, converted to radians.
 * @see http://asa.usno.navy.mil/static/files/2015/Astronomical_Constants_2015.pdf
 */
#define LAL_IEARTH 0.409092600600582871467239393761915655

/**
 * @brief Earth orbital eccentricity, dimensionless
 * @see E. Myles Standish and James G. Williams, Orbital Ephemerides of the Sun, Moon, and Planets ftp://ssd.jpl.nasa.gov/pub/eph/planets/ioms/ExplSupplChap8.pdf
 */
#define LAL_EEARTH 0.0167

/**
 * @brief Geocentric gravitational constant, m^3 s^-2 (TCB)
 * @see http://asa.usno.navy.mil/static/files/2015/Astronomical_Constants_2015.pdf
 */
#define LAL_GMEARTH_SI 3.986004418e+14

/**
 * @brief Solar equatorial radius, m
 * @see http://dx.doi.org/10.1088/0004-637X/750/2/135
 */
#define LAL_RSUN_SI 6.96342e+8

/**
 * @brief Solar luminosity, W
 * @see http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
 */
#define LAL_LSUN_SI 3.846e+26

/**
 * @brief Solar mass parameter, m^3 s^-2 (TCB)
 * @see http://asa.usno.navy.mil/static/files/2015/Astronomical_Constants_2015.pdf
 */
#define LAL_GMSUN_SI 1.32712442099e+20

/**
 * @brief Tropical year (2000), s
 * @see Borkowski, K. M., The Tropical Year and Solar Calendar, Journal of the Royal Astronomical Society of Canada, Vol. 85, NO. 3/JUN, P.121, 1991 http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1991JRASC..85..121B&data_type=PDF_HIGH&whole_paper=YES&type=PRINTER&filetype=.pdf
 */
#define LAL_YRTROP_SI 31556925.1874707200

/**
 * @brief Sidereal year (2000), s
 * @see http://hpiers.obspm.fr/eop-pc/models/constants.html
 */
#define LAL_YRSID_SI 31558149.763545600
/*@}*/

/**
 * @name Derived astrophysical parameters
 * The following constants are derived from the primary
 * astrophysical constants.  When not dimensionless, they are
 * given in the SI units shown.  Precision beyond the
 * accuracy is retained for these constants in order
 * that equivalent combinations yield the same value.
 */
/*@{*/

/**
 * @brief Cosine of Earth inclination (2000)
 * @details
 * LAL_COSIEARTH = cos(LAL_IEARTH)
 */
#define LAL_COSIEARTH 0.917482143065241841533315838574859003

/**
 * @brief Sine of Earth inclination (2000)
 * LAL_SINIEARTH = sin(LAL_IEARTH)
 */
#define LAL_SINIEARTH 0.397776969112605992551264763661918798

/**
 * @brief Earth mass, kg
 * @details
 * LAL_MEARTH_SI = LAL_GMEARTH_SI / LAL_G_SI
 */
#define LAL_MEARTH_SI 5.972580130779281493113409970871342436e24

/**
 * @brief Solar mass, kg
 * @details
 * LAL_MSUN_SI = LAL_GMSUN_SI / LAL_G_SI
 */
#define LAL_MSUN_SI 1.988546954961461467461011951140572744e30

/**
 * @brief Geometrized solar mass, m
 * @details
 * LAL_MRSUN_SI = LAL_GMSUN_SI / (LAL_C_SI * LAL_C_SI)
 */
#define LAL_MRSUN_SI 1.476625061404649406193430731479084713e3

/**
 * @brief Geometrized solar mass, s
 * @details
 * LAL_MTSUN_SI = LAL_GMSUN_SI / (LAL_C_SI * LAL_C_SI * LAL_C_SI)
 */
#define LAL_MTSUN_SI 4.925491025543575903411922162094833998e-6

/**
 * @brief Mean sidereal day, s
 * @details
 * LAL_DAYSID_SI = LAL_DAYJUL_SI / LAL_SOL_SID
 */
#define LAL_DAYSID_SI 86164.090530832885726908802698342011678582
/*@}*/

/**
 * @name Cosmological parameters
 * The following cosmological parameters are derived from measurements of
 * the Hubble expansion rate and of the cosmic microwave background radiation.
 * In what follows, the normalized Hubble constant \f$h_0\f$ is equal to the
 * actual Hubble constant \f$H_0\f$ divided by
 * \f$\langle H \rangle=100\,\mathrm{km}\,\mathrm{s}^{-1}\mathrm{Mpc}^{-1}\f$.
 * Thus the Hubble constant can be written as:
 * \f$H_0 = \langle H \rangle * h_0\f$.  Similarly, the critical energy density
 * \f$\rho_c\f$ required for spatial flatness is given by:
 * \f$\rho_c = \langle\rho\rangle h_0^2\f$.
 * Current estimates give \f$h_0\f$ a value of around 0.69 
 * which is what is assumed below.
 * All values are in the SI units shown.
 * @see http://arxiv.org/abs/1303.5062
 * @see http://dx.doi.org/10.1088/0067-0049/208/2/20
 */
/*@{*/

/**
 * @brief Hubble constant prefactor, s^-1
 * @details
 * LAL_H0FAC_SI = 100 km s^-1 Mpc^-1
 */
#define LAL_H0FAC_SI 3.240779289444365023237687716352957261e-18

/**
 * @brief Approximate Hubble constant, s^-1
 * @details
 * LAL_H0_SI = h0 * LAL_H0FAC_SI
 *
 * where h0 is approximately 0.69 (the value adopted here).
 * @see http://arxiv.org/abs/1303.5062
 * @see http://dx.doi.org/10.1088/0067-0049/208/2/20
 */
#define LAL_H0_SI (0.69 * LAL_H0FAC_SI)

/**
 * @brief Critical energy density prefactor, J m^-3
 * @details
 * LAL_RHOCFAC_SI = 3 * (LAL_H0FAC_SI * LAL_C_SI)^2 / (8 * LAL_PI * LAL_G_SI)
 */
#define LAL_RHOCFAC_SI 1.688285614132491133606558221644153036e-9

/**
 * @brief Approximate critical energy density, J m^-3
 * @details
 * LAL_RHOC_SI = h0^2 * LAL_RHOCFAC_SI
 *
 * where h0 is approximately 0.69 (the value adopted here).
 * @see http://arxiv.org/abs/1303.5062
 * @see http://dx.doi.org/10.1088/0067-0049/208/2/20
 */
#define LAL_RHOC_SI (0.69 * 0.69 * LAL_RHOCFAC_SI)

/**
 * @brief Cosmic microwave background radiation temperature, K
 * @see http://dx.doi.org/10.1088/0004-637X/707/2/916
 */
#define LAL_TCMB_SI 2.72548

/**
 * @brief Solar velocity with respect to the cosmic microwave background radiation, m s^-1
 * @details
 * Adopted value is v/c = 0.0012338
 * @see http://dx.doi.org/10.1088/0004-637X/707/2/916
 */
#define LAL_VCMB_SI 369883.9346804

/**
 * @brief Number density of cosmic microwave background radiation photons, m^-3
 * @details
 * LAL_NCMB_SI = 16 * zeta(3) * LAL_PI * (LAL_K_SI * LAL_TCMB_SI / (LAL_C_SI * LAL_H_SI))^3
 *
 * where zeta is the Riemann zeta function and zeta(3) is Apery's constant.
 * @see http://oeis.org/A002117
 */
#define LAL_NCMB_SI 4.107177389981626568308052789272762276e8

/**
 * @brief Entropy density of cosmic microwave background radiation, J K^-1 m^-3
 * @details
 * LAL_SCMB_SI = 4 * LAL_PI^2 * LAL_K_SI * (LAL_K_SI * LAL_TCMB_SI / (LAL_C_SI * LAL_HBAR_SI))^3 / 45
 */
#define LAL_SCMB_SI 2.042295723200401840880167645233175756e-14
/*@}*/

/*@}*/
/*@}*/
#ifdef  __cplusplus
}
#endif
#endif /* _LALCONSTANTS_H */ 
