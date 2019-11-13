import decimal

D = decimal.Decimal

precision = 36 # should be sufficient for quadruple-precision
extra_precision = 5 # some extra precision for correct rounding
quantize = D('1.' + '0'*precision)
decimal.getcontext().prec = precision + extra_precision
decimal.getcontext().capitals = 0


def as_str(x):
	expn = x.log10().quantize(D('1.'), rounding=decimal.ROUND_DOWN)
	if expn == 0:
		return str(x.quantize(quantize))
	elif expn < 0:
		expn -= D(1)
	return str(x.scaleb(-expn).quantize(quantize)) + 'e' + str(expn)

def wien_factor_function(x):
	return x * x.exp() - D(5) * x.exp() + D(5)

def wien_factor_function_deriv(x):
	return (x - D(4)) * x.exp()

def compute_wien_factor():
	""" Use Newton's method. """
	x0 = D(5)
	while True:
		x = x0 - wien_factor_function(x0) / wien_factor_function_deriv(x0)
		if  x.quantize(quantize) == x0.quantize(quantize):
			return x
		x0 = x

def compute_pi():
	pi0 = D('0')
	k = 0
	while True:
		term = D(4) / D(8*k + 1)
		term -= D(2) / D(8*k + 4)
		term -= D(1) / D(8*k + 5)
		term -= D(1) / D(8*k + 6)
		term *= D(16) ** D(-k)
		pi = pi0 + term
		if  pi.quantize(quantize) == pi0.quantize(quantize):
			return pi
		pi0 = pi
		k = k + 1

# computes Apery's constant = Riemann zeta function zeta(3)
def compute_zeta3():
	n, lasts, s, fact2n, factn, sign = 0, D(1), D(0), D(1), D(1), -1
	while s != lasts:
		lasts = s
		n += 1
		factn *= D(n)
		fact2n *= D(2*n) * D(2*n - 1)
		bincoef = fact2n / (factn * factn)
		sign *= -1
		s += sign / (bincoef * D(n)**D(3))
	return D('2.5')*s


def cos(x):
	i, lasts, s, fact, num, sign = 0, 0, 1, 1, 1, 1
	while s != lasts:
		lasts = s
		i += 2
		fact *= i * (i-1)
		num *= x * x
		sign *= -1
		s += num / fact * sign
	return +s

def sin(x):
	i, lasts, s, fact, num, sign = 1, 0, x, 1, x, 1
	while s != lasts:
		lasts = s
		i += 2
		fact *= i * (i-1)
		num *= x * x
		sign *= -1
		s += num / fact * sign
	return +s


print r'''/*
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
 *//*@{*/'''

print r'''
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
 */'''
print '/*@{*/'

LAL_REAL4_BITS = 32 # bits in REAL4
LAL_REAL4_MANT = 24 # bits in REAL4 mantissa
LAL_REAL4_EXPN = LAL_REAL4_BITS - LAL_REAL4_MANT - 1 # bits in REAL4 exponent
LAL_REAL4_MAX  = (2 - 2**(-LAL_REAL4_MANT+1)) * 2**(2**LAL_REAL4_EXPN - 1)
LAL_REAL4_MIN  = 2**(2 - 2**LAL_REAL4_EXPN)
LAL_REAL4_EPS  = 2**(1 - LAL_REAL4_MANT)

LAL_REAL8_BITS = 64 # bits in REAL8
LAL_REAL8_MANT = 53 # bits in REAL8 mantissa
LAL_REAL8_EXPN = LAL_REAL8_BITS - LAL_REAL8_MANT - 1 # bits in REAL8 exponent
LAL_REAL8_MAX  = (2 - 2**(-LAL_REAL8_MANT+1)) * 2**(2**LAL_REAL8_EXPN - 1)
LAL_REAL8_MIN  = 2**(2 - 2**LAL_REAL8_EXPN)
LAL_REAL8_EPS  = 2**(1 - LAL_REAL8_MANT)

print '#if __STDC_VERSION__ >= 199901L'

print '#define LAL_REAL4_MANT',
print LAL_REAL4_MANT,
print '/**< Bits of precision in the mantissa of a REAL4 */'

print '#define LAL_REAL4_MAX',
print float(LAL_REAL4_MAX).hex(),
print '/**< Largest normalized REAL4 number (2-2^-23)*2^127 */'

print '#define LAL_REAL4_MIN',
print float(LAL_REAL4_MIN).hex(),
print '/**< Smallest normalized REAL4 number 2^-126 */'

print '#define LAL_REAL4_EPS',
print float(LAL_REAL4_EPS).hex(),
print '/**< Difference between 1 and the next resolvable REAL4 2^-23 */'

print '#define LAL_REAL8_MANT',
print LAL_REAL8_MANT,
print '/**< Bits of precision in the mantissa of a REAL8 */'

print '#define LAL_REAL8_MAX',
print float(LAL_REAL8_MAX).hex(),
print '/**< Largest normalized REAL8 number (2-2^-52)*2^1023 */'

print '#define LAL_REAL8_MIN',
print float(LAL_REAL8_MIN).hex(),
print '/**< Smallest normalized REAL8 number 2^-1022 */'

print '#define LAL_REAL8_EPS',
print float(LAL_REAL8_EPS).hex(),
print '/**< Difference between 1 and the next resolvable REAL8 2^-52 */'

print '#else'
print '#define LAL_REAL4_MANT', LAL_REAL4_MANT
print '#define LAL_REAL4_MAX', repr(LAL_REAL4_MAX)
print '#define LAL_REAL4_MIN', repr(LAL_REAL4_MIN)
print '#define LAL_REAL4_EPS', repr(LAL_REAL4_EPS)
print '#define LAL_REAL8_MANT', LAL_REAL8_MANT
print '#define LAL_REAL8_MAX', as_str(D(LAL_REAL8_MAX))
print '#define LAL_REAL8_MIN', as_str(D(LAL_REAL8_MIN))
print '#define LAL_REAL8_EPS', as_str(D(LAL_REAL8_EPS))
print '#endif'

print '/*@}*/'

print r'''
/**
 * @name Integer constants
 * Extremal integer values, all expressed as unsigned long long.
 */'''
print '/*@{*/'
print '#define LAL_UINT8_MAX   LAL_UINT8_C(%u)' % (2**(8*8)-1)
print '#define LAL_UINT4_MAX   LAL_UINT8_C(%u)' % (2**(8*4)-1)
print '#define LAL_UINT2_MAX   LAL_UINT8_C(%u)' % (2**(8*2)-1)
print '#define LAL_INT8_MAX    LAL_UINT8_C(%u)' % (2**(8*8-1)-1)
print '#define LAL_INT4_MAX    LAL_UINT8_C(%u)' % (2**(8*4-1)-1)
print '#define LAL_INT2_MAX    LAL_UINT8_C(%u)' % (2**(8*2-1)-1)
print '/*@}*/'



print r'''
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
 */'''

LAL_EXPGAMMA = D('1.7810724179901979852365041031071795491696452143034302053')

LAL_E = D('1').exp()
LAL_LN2 = D('2').ln()
LAL_LN10 = D('10').ln()
LAL_LOG2E = D('1') / LAL_LN2
LAL_LOG10E = LAL_E.log10()
LAL_SQRT2 = D('2').sqrt()
LAL_SQRT1_2 = D('1') / LAL_SQRT2
LAL_GAMMA = LAL_EXPGAMMA.ln()

print '/*@{*/'

print '#define LAL_E        ',
print str(LAL_E.quantize(quantize))

print '#define LAL_LOG2E    ',
print str(LAL_LOG2E.quantize(quantize))

print '#define LAL_LOG10E   ',
print str(LAL_LOG10E.quantize(quantize))

print '#define LAL_LN2      ',
print str(LAL_LN2.quantize(quantize))

print '#define LAL_LN10     ',
print str(LAL_LN10.quantize(quantize))

print '#define LAL_SQRT2    ',
print str(LAL_SQRT2.quantize(quantize))

print '#define LAL_SQRT1_2  ',
print str(LAL_SQRT1_2.quantize(quantize))

print '#define LAL_GAMMA    ',
print str(LAL_GAMMA.quantize(quantize))

print '#define LAL_EXPGAMMA ',
print str(LAL_EXPGAMMA.quantize(quantize))

LAL_PI = compute_pi()
LAL_TWOPI = D('2') * LAL_PI
LAL_PI_2 = LAL_PI / D('2')
LAL_PI_4 = LAL_PI / D('4')
LAL_1_PI = D('1') / LAL_PI
LAL_2_PI = D('2') / LAL_PI
LAL_2_SQRTPI = D('2') / LAL_PI.sqrt()
LAL_PI_180 = LAL_PI / D('180')
LAL_180_PI = D('180') / LAL_PI
LAL_LNPI = LAL_PI.ln()

print "/* Assuming we're not near a black hole or in Tennessee... */"

print '#define LAL_PI       ',
print str(LAL_PI.quantize(quantize))

print '#define LAL_TWOPI    ',
print str(LAL_TWOPI.quantize(quantize))

print '#define LAL_PI_2     ',
print str(LAL_PI_2.quantize(quantize))

print '#define LAL_PI_4     ',
print str(LAL_PI_4.quantize(quantize))

print '#define LAL_1_PI     ',
print str(LAL_1_PI.quantize(quantize))

print '#define LAL_2_PI     ',
print str(LAL_2_PI.quantize(quantize))

print '#define LAL_2_SQRTPI ',
print str(LAL_2_SQRTPI.quantize(quantize))

print '#define LAL_PI_180   ',
print as_str(LAL_PI_180)

print '#define LAL_180_PI   ',
print str(LAL_180_PI.quantize(quantize))

print '#define LAL_LNPI     ',
print str(LAL_LNPI.quantize(quantize))

print '/*@}*/'

print r'''
/**
 * @name Exact physical constants
 * The following physical constants are defined to have exact values.
 * The dimensions in SI units are as shown.
 * @see http://dx.doi.org/10.1103/RevModPhys.84.1527
 */'''

print '/*@{*/'
LAL_C_SI = D('299792458')
LAL_MU0_SI = D('4') * LAL_PI * D('1E-7')
LAL_EPSILON0_SI = D('1') / (LAL_MU0_SI * LAL_C_SI * LAL_C_SI)
LAL_GEARTH_SI = D('9.80665')
LAL_PATM_SI = D('101325')

print '#define LAL_C_SI',
print str(LAL_C_SI) + 'e0',
print '/**< Speed of light in vacuo, m s^-1 */'

print '#define LAL_EPSILON0_SI',
print as_str(LAL_EPSILON0_SI),
print '/**< Permittivity of free space, C^2 N^-1 m^-2 */'

print '#define LAL_MU0_SI',
print as_str(LAL_MU0_SI),
print '/**< Permeability of free space, N A^-2 */'

print '#define LAL_GEARTH_SI',
print LAL_GEARTH_SI,
print '/**< Standard gravity, m s^-2 */'

print '#define LAL_PATM_SI',
print str(LAL_PATM_SI) + 'e0',
print '/**< Standard atmosphere, Pa */'

print '/*@}*/'

print r'''
/**
 * @name Primary physical constants
 * These physical constants are given to the precision
 * to which they are known.  Other physical constants
 * derived from these are given in the next section.
 * @see http://dx.doi.org/10.1103/RevModPhys.84.1527
 */'''

print '/*@{*/'

LAL_QE_SI = 			D('1.60217656535E-19')
print '#define LAL_QE_SI',
print str(LAL_QE_SI),
print '/**< Electron charge, C */'

LAL_ALPHA = 			D('7.297352569824E-3')
print '#define LAL_ALPHA',
print str(LAL_ALPHA),
print '/**< Fine structure constant, dimensionless */'

LAL_RYD_SI = 			D('10973731.56853955')
print '#define LAL_RYD_SI',
print str(LAL_RYD_SI),
print '/**< Rydberg constant, m^-1 */'

LAL_MP_ME = 			D('1836.1526724575')
print '#define LAL_MP_ME',
print str(LAL_MP_ME),
print '/**< Proton-electron mass ratio, dimensionless */'

LAL_MP_AMU = 			D('1.00727646681290')
print '#define LAL_MP_AMU',
print str(LAL_MP_AMU),
print '/**< Proton molar mass, kg mol^-1 */'

LAL_R_SI = 			D('8.314462175')
print '#define LAL_R_SI',
print str(LAL_R_SI),
print '/**< Molar gas constant, J mol^-1 K^-1 */'

LAL_G_SI = 			D('6.67384e-11')
print '#define LAL_G_SI',
print str(LAL_G_SI),
print '/**< Gravitational constant, N m^2 kg^-2 */'

print '/*@}*/'

print r'''
/**
 * @name Derived physical constants
 * The following constants are derived from the primary
 * physical constants.  When not dimensionless, they are
 * given in the SI units shown.  Precision beyond the
 * accuracy is retained for these constants in order
 * that equivalent combinations yield the same value.
 */'''

print '/*@{*/'

LAL_HBAR_SI = ((LAL_QE_SI * LAL_QE_SI) / (D('4') * LAL_PI * LAL_EPSILON0_SI * LAL_C_SI * LAL_ALPHA))
print '''
/**
 * @brief Reduced Planck constant, J s
 * @details
 * LAL_HBAR_SI = ((LAL_QE_SI * LAL_QE_SI) / (4 * LAL_PI * LAL_EPSILON0_SI * LAL_C_SI * LAL_ALPHA))
 */'''
print '#define LAL_HBAR_SI',
print as_str(LAL_HBAR_SI)

LAL_H_SI = LAL_TWOPI * LAL_HBAR_SI
print '''
/**
 * @brief Planck constant, J s
 * @details
 * LAL_H_SI = LAL_TWOPI * LAL_HBAR_SI
 */'''
print '#define LAL_H_SI',
print as_str(LAL_H_SI)

LAL_MPL_SI =			(LAL_HBAR_SI * LAL_C_SI / LAL_G_SI).sqrt()
print '''
/**
 * @brief Planck mass, kg
 * @details
 * LAL_MPL_SI =	sqrt(LAL_HBAR_SI * LAL_C_SI / LAL_G_SI)
 */'''
print '#define LAL_MPL_SI',
print as_str(LAL_MPL_SI)

LAL_LPL_SI =			(LAL_HBAR_SI / (LAL_MPL_SI * LAL_C_SI))
print '''
/**
 * @brief Planck length, m
 * @details
 * LAL_LPL_SI =	(LAL_HBAR_SI / (LAL_MPL_SI * LAL_C_SI))
 */'''
print '#define LAL_LPL_SI',
print as_str(LAL_LPL_SI)

LAL_TPL_SI =			(LAL_LPL_SI / LAL_C_SI)
print '''
/**
 * @brief Planck time, s
 * @details
 * LAL_TPL_SI =	(LAL_LPL_SI / LAL_C_SI)
 */'''
print '#define LAL_TPL_SI',
print as_str(LAL_TPL_SI)

LAL_ME_SI = 			((D('2') * LAL_RYD_SI * LAL_H_SI) / (LAL_C_SI * LAL_ALPHA * LAL_ALPHA))
print '''
/**
 * @brief Electron mass, kg
 * @details
 * LAL_ME_SI = ((2 * LAL_RYD_SI * LAL_H_SI) / (LAL_C_SI * LAL_ALPHA * LAL_ALPHA))
 */'''
print '#define LAL_ME_SI',
print as_str(LAL_ME_SI)

LAL_MP_SI =			(LAL_ME_SI * LAL_MP_ME)
print '''
/**
 * @brief Proton mass, kg
 * @details
 * LAL_MP_SI = (LAL_ME_SI * LAL_MP_ME)
 */'''
print '#define LAL_MP_SI',
print as_str(LAL_MP_SI)

LAL_AMU_SI =			(LAL_MP_SI / LAL_MP_AMU)
print '''
/**
 * @brief Atomic mass uint, kg
 * @details
 * LAL_AMU_SI = (LAL_MP_SI / LAL_MP_AMU)
 */'''
print '#define LAL_AMU_SI',
print as_str(LAL_AMU_SI)

LAL_AB_SI =			(LAL_ALPHA / (D('4') * LAL_PI * LAL_RYD_SI))
print '''
/**
 * @brief Bohr radius, m
 * @details
 * LAL_AB_SI = (LAL_ALPHA / (4 * LAL_PI * LAL_RYD_SI))
 */'''
print '#define LAL_AB_SI',
print as_str(LAL_AB_SI)

LAL_LAMBDAE_SI =		(D('2') * LAL_PI * LAL_ALPHA * LAL_AB_SI)
print '''
/**
 * @brief Electron Compton wavelength, m
 * @details
 * LAL_LAMBDAE_SI = (2 * LAL_PI * LAL_ALPHA * LAL_AB_SI)
 */'''
print '#define LAL_LAMBDAE_SI',
print as_str(LAL_LAMBDAE_SI)

LAL_RE_SI = 			(LAL_ALPHA * LAL_ALPHA * LAL_AB_SI)
print '''
/**
 * @brief Classical electron radius, m
 * @details
 * LAL_RE_SI = (LAL_ALPHA * LAL_ALPHA * LAL_AB_SI)
 */'''
print '#define LAL_RE_SI',
print as_str(LAL_RE_SI)

LAL_MUB_SI = 			(LAL_LAMBDAE_SI * LAL_C_SI * LAL_QE_SI / (D('4') * LAL_PI))
print '''
/**
 * @brief Bohr magneton, J T^-1
 * LAL_MUB_SI = (LAL_LAMBDAE_SI * LAL_C_SI * LAL_QE_SI / (4 * LAL_PI))
 * @details
 */'''
print '#define LAL_MUB_SI',
print as_str(LAL_MUB_SI)

LAL_MUN_SI =			(LAL_MUB_SI / LAL_MP_ME)
print '''
/**
 * @brief Nuclear magneton, J T^-1
 * @details
 * LAL_MUN_SI =	(LAL_MUB_SI / LAL_MP_ME)
 */'''
print '#define LAL_MUN_SI',
print as_str(LAL_MUN_SI)

LAL_MOL =			(D('0.001') * LAL_MP_AMU / LAL_MP_SI)
print '''
/**
 * @brief Avogadro constant, dimensionless
 * @details
 * LAL_MOL = (0.001 * LAL_MP_AMU / LAL_MP_SI)
 */'''
print '#define LAL_MOL',
print as_str(LAL_MOL)

LAL_K_SI =			(LAL_R_SI / LAL_MOL)
print '''
/**
 * @brief Boltzmann constant, J K^-1
 * @details
 * LAL_K_SI = (LAL_R_SI / LAL_MOL)
 */'''
print '#define LAL_K_SI',
print as_str(LAL_K_SI)

LAL_SIGMA_SI =			((LAL_PI * LAL_PI * LAL_K_SI * LAL_K_SI * LAL_K_SI * LAL_K_SI) / (D('60') * LAL_HBAR_SI * LAL_HBAR_SI * LAL_HBAR_SI * LAL_C_SI * LAL_C_SI))
print '''
/**
 * @brief Stefan-Boltzmann constant, W m^-2 K^-4
 * @details
 * LAL_SIGMA_SI = ((LAL_PI * LAL_PI * LAL_K_SI * LAL_K_SI * LAL_K_SI * LAL_K_SI) / (60 * LAL_HBAR_SI * LAL_HBAR_SI * LAL_HBAR_SI * LAL_C_SI * LAL_C_SI))
 */'''
print '#define LAL_SIGMA_SI',
print as_str(LAL_SIGMA_SI)

# Second radiation constant (m K)
LAL_C2RAD_SI =			(LAL_H_SI * LAL_C_SI / LAL_K_SI)
print '''
/**
 * @brief Second radiation constant, m K
 * @details
 * LAL_C2RAD_SI = (LAL_H_SI * LAL_C_SI / LAL_K_SI)
 */'''
print '#define LAL_C2RAD_SI',
print as_str(LAL_C2RAD_SI)

# Factor in Wein displacement law: x e^x = 5 (e^x - 1)
LAL_C2RAD_BWIEN = compute_wien_factor()
#print '#define LAL_C2RAD_BWIEN',
#print as_str(LAL_C2RAD_BWIEN),

# Wein displacement law constant (m K)
LAL_BWIEN_SI = 			(LAL_C2RAD_SI / LAL_C2RAD_BWIEN)
print '''
/**
 * @brief Wien displacement law constant, m K
 * @details
 * LAL_BWIEN_SI = (LAL_C2RAD_SI / X)
 *
 * where the factor X satisfies
 *
 * X * exp(X) = 5 * (exp(X) - 1)
 */'''
print '#define LAL_BWIEN_SI',
print as_str(LAL_BWIEN_SI)

print '/*@}*/'

print r'''
/**
 * @name Exact astrophysical parameters
 * The following astrophysical constants are defined to have exact values.
 * The dimensions in SI units are as shown.
 * @see http://asa.usno.navy.mil/static/files/2015/Astronomical_Constants_2015.pdf
 */'''

print '/*@{*/'

LAL_SOL_SID = D('1.002737909350795')
LAL_DAYJUL_SI = D('86400')
LAL_YRJUL_SI = D('31557600')
LAL_LYR_SI = LAL_YRJUL_SI * LAL_C_SI
LAL_AU_SI = D('149597870700')
LAL_PC_SI = LAL_AU_SI * D('3600') * LAL_180_PI

print '#define LAL_SOL_SID',
print LAL_SOL_SID,
print '/**< Ratio of mean solar day to sidereal day, dimensionless */'

print '#define LAL_DAYJUL_SI',
print str(LAL_DAYJUL_SI) + 'e0',
print '/**< Julian day, s */'

print '#define LAL_YRJUL_SI',
print str(LAL_YRJUL_SI) + 'e0',
print '/**< Julian year, s */'

print '#define LAL_LYR_SI',
print LAL_LYR_SI,
print '/**< (Julian) Lightyear, m */'

print '#define LAL_AU_SI',
print str(LAL_AU_SI) + 'e0',
print '/**< Astronomical unit, m */'

print '#define LAL_PC_SI',
print as_str(LAL_PC_SI),
print '/**< Parsec, m */'

print '/*@}*/'

print r'''
/**
 * @name Primary astrophysical parameters
 * These astrophysical constants are given to the precision
 * to which they are known.  Other physical constants
 * derived from these are given in the next section.
 */'''

print '/*@{*/'

LAL_GMSUN_SI = D('1.32712442099E20')
LAL_GMEARTH_SI = D('3.986004418E14')
LAL_REARTH_SI = D('6378136.6')
LAL_AWGS84_SI = D('6378137')
LAL_FWGS84 = D('1') / D('298.257223563')
LAL_BWGS84_SI = LAL_AWGS84_SI * (D('1') - LAL_FWGS84)
LAL_IEARTH = D('84381.406') * LAL_PI_180 / D('3600')
LAL_EEARTH = D('0.0167')
LAL_RSUN_SI = D('696342E3')
LAL_LSUN_SI = D('3.846E26')
LAL_YRTROP_SI = D('365.2421896698') * LAL_DAYJUL_SI
LAL_YRSID_SI = D('365.256363004') * LAL_DAYJUL_SI

print '''
/**
 * @brief Earth equatorial radius, m
 * @see http://asa.usno.navy.mil/static/files/2015/Astronomical_Constants_2015.pdf
 */'''
print '#define LAL_REARTH_SI',
print LAL_REARTH_SI

print '''
/**
 * @brief Semimajor axis of WGS-84 Reference Ellipsoid, m
 * @see Department of Defense World Geodedic System 1984 http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf
 */'''
print '#define LAL_AWGS84_SI',
print str(LAL_AWGS84_SI) + 'e0'

print '''
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
 */'''
print '#define LAL_BWGS84_SI',
print LAL_BWGS84_SI.quantize(D('1000000.000'))

print '''
/**
 * @brief Earth inclination (2000), radians
 * @details
 * This is the measured value of the mean obliquity of the
 * ecliptic, 84381.406 arcseconds, converted to radians.
 * @see http://asa.usno.navy.mil/static/files/2015/Astronomical_Constants_2015.pdf
 */'''
print '#define LAL_IEARTH',
print str(LAL_IEARTH.quantize(quantize))

print '''
/**
 * @brief Earth orbital eccentricity, dimensionless
 * @see E. Myles Standish and James G. Williams, Orbital Ephemerides of the Sun, Moon, and Planets ftp://ssd.jpl.nasa.gov/pub/eph/planets/ioms/ExplSupplChap8.pdf
 */'''
print '#define LAL_EEARTH',
print LAL_EEARTH

print '''
/**
 * @brief Geocentric gravitational constant, m^3 s^-2 (TCB)
 * @see http://asa.usno.navy.mil/static/files/2015/Astronomical_Constants_2015.pdf
 */'''
print '#define LAL_GMEARTH_SI',
print LAL_GMEARTH_SI

print '''
/**
 * @brief Solar equatorial radius, m
 * @see http://dx.doi.org/10.1088/0004-637X/750/2/135
 */'''
print '#define LAL_RSUN_SI',
print LAL_RSUN_SI

print '''
/**
 * @brief Solar luminosity, W
 * @see http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
 */'''
print '#define LAL_LSUN_SI',
print LAL_LSUN_SI

print '''
/**
 * @brief Solar mass parameter, m^3 s^-2 (TCB)
 * @see http://asa.usno.navy.mil/static/files/2015/Astronomical_Constants_2015.pdf
 */'''
print '#define LAL_GMSUN_SI',
print LAL_GMSUN_SI

print '''
/**
 * @brief Tropical year (2000), s
 * @see Borkowski, K. M., The Tropical Year and Solar Calendar, Journal of the Royal Astronomical Society of Canada, Vol. 85, NO. 3/JUN, P.121, 1991 http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1991JRASC..85..121B&data_type=PDF_HIGH&whole_paper=YES&type=PRINTER&filetype=.pdf
 */'''
print '#define LAL_YRTROP_SI',
print LAL_YRTROP_SI

print '''
/**
 * @brief Sidereal year (2000), s
 * @see http://hpiers.obspm.fr/eop-pc/models/constants.html
 */'''
print '#define LAL_YRSID_SI',
print LAL_YRSID_SI

print '/*@}*/'

print r'''
/**
 * @name Derived astrophysical parameters
 * The following constants are derived from the primary
 * astrophysical constants.  When not dimensionless, they are
 * given in the SI units shown.  Precision beyond the
 * accuracy is retained for these constants in order
 * that equivalent combinations yield the same value.
 */'''

print '/*@{*/'

LAL_COSIEARTH = cos(LAL_IEARTH)
LAL_SINIEARTH = sin(LAL_IEARTH)
LAL_MEARTH_SI = LAL_GMEARTH_SI / LAL_G_SI
LAL_MSUN_SI = LAL_GMSUN_SI / LAL_G_SI
LAL_MRSUN_SI = LAL_GMSUN_SI / (LAL_C_SI * LAL_C_SI)
LAL_MTSUN_SI = LAL_GMSUN_SI / (LAL_C_SI * LAL_C_SI * LAL_C_SI)
LAL_DAYSID_SI = LAL_DAYJUL_SI / LAL_SOL_SID

print '''
/**
 * @brief Cosine of Earth inclination (2000)
 * @details
 * LAL_COSIEARTH = cos(LAL_IEARTH)
 */'''
print '#define LAL_COSIEARTH',
print str(LAL_COSIEARTH.quantize(quantize))

print '''
/**
 * @brief Sine of Earth inclination (2000)
 * LAL_SINIEARTH = sin(LAL_IEARTH)
 */'''
print '#define LAL_SINIEARTH',
print str(LAL_SINIEARTH.quantize(quantize))

print '''
/**
 * @brief Earth mass, kg
 * @details
 * LAL_MEARTH_SI = LAL_GMEARTH_SI / LAL_G_SI
 */'''
print '#define LAL_MEARTH_SI',
print as_str(LAL_MEARTH_SI)

print '''
/**
 * @brief Solar mass, kg
 * @details
 * LAL_MSUN_SI = LAL_GMSUN_SI / LAL_G_SI
 */'''
print '#define LAL_MSUN_SI',
print as_str(LAL_MSUN_SI)

print '''
/**
 * @brief Geometrized solar mass, m
 * @details
 * LAL_MRSUN_SI = LAL_GMSUN_SI / (LAL_C_SI * LAL_C_SI)
 */'''
print '#define LAL_MRSUN_SI',
print as_str(LAL_MRSUN_SI)

print '''
/**
 * @brief Geometrized solar mass, s
 * @details
 * LAL_MTSUN_SI = LAL_GMSUN_SI / (LAL_C_SI * LAL_C_SI * LAL_C_SI)
 */'''
print '#define LAL_MTSUN_SI',
print as_str(LAL_MTSUN_SI)

print '''
/**
 * @brief Mean sidereal day, s
 * @details
 * LAL_DAYSID_SI = LAL_DAYJUL_SI / LAL_SOL_SID
 */'''
print '#define LAL_DAYSID_SI',
print LAL_DAYSID_SI.quantize(quantize)

print '/*@}*/'

h0 = 0.69 # current reasonable guess for the normalized Hubble constant
print r'''
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
 * Current estimates give \f$h_0\f$ a value of around''',
print h0,
print r'''
 * which is what is assumed below.
 * All values are in the SI units shown.
 * @see http://arxiv.org/abs/1303.5062
 * @see http://dx.doi.org/10.1088/0067-0049/208/2/20
 */'''
print '/*@{*/'


print '''
/**
 * @brief Hubble constant prefactor, s^-1
 * @details
 * LAL_H0FAC_SI = 100 km s^-1 Mpc^-1
 */'''
LAL_H0FAC_SI = D('1') / ((D('10') * LAL_PC_SI))
print '#define LAL_H0FAC_SI',
print as_str(LAL_H0FAC_SI)

print '''
/**
 * @brief Approximate Hubble constant, s^-1
 * @details
 * LAL_H0_SI = h0 * LAL_H0FAC_SI
 *
 * where h0 is approximately %g (the value adopted here).
 * @see http://arxiv.org/abs/1303.5062
 * @see http://dx.doi.org/10.1088/0067-0049/208/2/20
 */''' % h0
print '#define LAL_H0_SI (%g * LAL_H0FAC_SI)' % h0

print '''
/**
 * @brief Critical energy density prefactor, J m^-3
 * @details
 * LAL_RHOCFAC_SI = 3 * (LAL_H0FAC_SI * LAL_C_SI)^2 / (8 * LAL_PI * LAL_G_SI)
 */'''
LAL_RHOCFAC_SI = D('3') * (LAL_H0FAC_SI * LAL_C_SI)**(D('2')) / (D('8') * LAL_PI * LAL_G_SI)
print '#define LAL_RHOCFAC_SI',
print as_str(LAL_RHOCFAC_SI)

print '''
/**
 * @brief Approximate critical energy density, J m^-3
 * @details
 * LAL_RHOC_SI = h0^2 * LAL_RHOCFAC_SI
 *
 * where h0 is approximately %g (the value adopted here).
 * @see http://arxiv.org/abs/1303.5062
 * @see http://dx.doi.org/10.1088/0067-0049/208/2/20
 */''' % h0
print '#define LAL_RHOC_SI (%g * %g * LAL_RHOCFAC_SI)' % (h0,h0)

print '''
/**
 * @brief Cosmic microwave background radiation temperature, K
 * @see http://dx.doi.org/10.1088/0004-637X/707/2/916
 */'''
LAL_TCMB_SI = D('2.72548')
print '#define LAL_TCMB_SI',
print LAL_TCMB_SI

print '''
/**
 * @brief Solar velocity with respect to the cosmic microwave background radiation, m s^-1
 * @details
 * Adopted value is v/c = 0.0012338
 * @see http://dx.doi.org/10.1088/0004-637X/707/2/916
 */'''
LAL_VCMB_SI = D('0.0012338') * LAL_C_SI
print '#define LAL_VCMB_SI',
print LAL_VCMB_SI

print '''
/**
 * @brief Number density of cosmic microwave background radiation photons, m^-3
 * @details
 * LAL_NCMB_SI = 16 * zeta(3) * LAL_PI * (LAL_K_SI * LAL_TCMB_SI / (LAL_C_SI * LAL_H_SI))^3
 *
 * where zeta is the Riemann zeta function and zeta(3) is Apery's constant.
 * @see http://oeis.org/A002117
 */'''
zeta3 = compute_zeta3()
LAL_NCMB_SI = D(16) * zeta3 * LAL_PI * (LAL_K_SI * LAL_TCMB_SI / (LAL_C_SI * LAL_H_SI))**3
print '#define LAL_NCMB_SI',
print as_str(LAL_NCMB_SI)

print '''
/**
 * @brief Entropy density of cosmic microwave background radiation, J K^-1 m^-3
 * @details
 * LAL_SCMB_SI = 4 * LAL_PI^2 * LAL_K_SI * (LAL_K_SI * LAL_TCMB_SI / (LAL_C_SI * LAL_HBAR_SI))^3 / 45
 */'''
zeta3 = compute_zeta3()
LAL_SCMB_SI = D('4') * LAL_PI**D('2') * LAL_K_SI * (LAL_K_SI * LAL_TCMB_SI / (LAL_C_SI * LAL_HBAR_SI))**D('3') / D('45')
print '#define LAL_SCMB_SI',
print as_str(LAL_SCMB_SI)

print '/*@}*/'

print r'''
/*@}*/
/*@}*/
#ifdef  __cplusplus
}
#endif
#endif /* _LALCONSTANTS_H */ '''
