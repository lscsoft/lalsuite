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

/** \file
 * \ingroup std
 * \author Creighton, T. D.
 * \date $Id$
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

/********************************* <lalVerbatim file="LALConstantsHV">
Author: Creighton, T. D.
$Id$
********************************** </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALConstants.h}}
\label{s:LALConstants.h}

Provides standard numerical constants for LAL.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALConstants.h>
\end{verbatim}

\noindent This header defines a number of useful numerical constants
for use in LAL routines.  These constants come in three basic
flavours: arithmetic and mathematical constants, fundamental (or
defined) physical constants, and measured astrophysical and
cosmological parameters.

Note that, unlike the other headers in the \verb@std@ package, this
header is \emph{not} included automatically by the header
\verb@LALStdlib.h@.  Include it explicitly if you need any of these
constants.

</lalLaTeX> */

#ifndef _LALCONSTANTS_H
#define _LALCONSTANTS_H

#include <lal/LALRCSID.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (LALCONSTANTSH, "$Id$");

/* <lalLaTeX>

\subsection*{Mathematical constants}
\idx[Constant]{LAL\_REAL4\_MANT}
\idx[Constant]{LAL\_REAL4\_MAX}
\idx[Constant]{LAL\_REAL4\_MIN}
\idx[Constant]{LAL\_REAL4\_EPS}
\idx[Constant]{LAL\_REAL8\_MANT}
\idx[Constant]{LAL\_REAL8\_MAX}
\idx[Constant]{LAL\_REAL8\_MIN}
\idx[Constant]{LAL\_REAL8\_EPS}
\idx[Constant]{LAL\_E}
\idx[Constant]{LAL\_LOG2E}
\idx[Constant]{LAL\_LOG10E}
\idx[Constant]{LAL\_LN2}
\idx[Constant]{LAL\_LN10}
\idx[Constant]{LAL\_SQRT2}
\idx[Constant]{LAL\_SQRT1\_2}
\idx[Constant]{LAL\_GAMMA}
\idx[Constant]{LAL\_PI}
\idx[Constant]{LAL\_TWOPI}
\idx[Constant]{LAL\_PI\_2}
\idx[Constant]{LAL\_PI\_4}
\idx[Constant]{LAL\_1\_PI}
\idx[Constant]{LAL\_2\_PI}
\idx[Constant]{LAL\_2\_SQRTPI}
\idx[Constant]{LAL\_PI\_180}
\idx[Constant]{LAL\_180\_PI}

The following constants define the precision and range of
floating-point arithmetic in LAL.  They are taken from the IEEE
standard~754 for binary arithmetic.  All numbers are dimensionless.

\begin{center}
\begin{tabular}{|lll|}
\hline
Name & Value & Description \\
\hline
\tt LAL\_REAL4\_MANT & 24 &
	Bits in {\tt REAL4} mantissa \\
\tt LAL\_REAL4\_MAX  & $3.40282347\times10^{38}$ &
	Largest {\tt REAL4} \\
\tt LAL\_REAL4\_MIN  & $1.17549435\times10^{-38}$ &
	Smallest positive {\tt REAL4} \\
\tt LAL\_REAL4\_EPS  & $1.19209290\times10^{-7}$ &
	$2^{-(\mathtt{LAL\_REAL4\_MANT}-1)}$ \\
\hline
\tt LAL\_REAL8\_MANT & 53 &
	Bits in {\tt REAL8} mantissa \\
\tt LAL\_REAL8\_MAX  & $1.7976931348623157\times10^{308}$ &
	Largest {\tt REAL8} \\
\tt LAL\_REAL8\_MIN  & $2.2250738585072014\times10^{-308}$ &
	Smallest positive {\tt REAL8} \\
\tt LAL\_REAL8\_EPS  & $2.2204460492503131\times10^{-16}$ &
	$2^{-(\mathtt{LAL\_REAL8\_MANT}-1)}$ \\
\hline
\end{tabular}
\end{center}

\noindent\verb@LAL_REAL4_EPS@ and \verb@LAL_REAL8_EPS@ can be thought
of as the difference between 1 and the next representable \verb@REAL4@
or \verb@REAL8@ number.

\vspace{3ex}

</lalLaTeX> */

/** \name Floating-point constants
 * The following constants define the precision and range of
 * floating-point arithmetic in LAL.  They are taken from the IEEE
 * standard 754 for binary arithmetic.  All numbers are dimensionless. */
/*@{*/
#define LAL_REAL4_MANT 24 /**< Bits of precision in the mantissa of a REAL4 */
#define LAL_REAL4_MAX 3.40282347e+38 /**< Largest REAL4 */
#define LAL_REAL4_MIN 1.17549435e-38 /**< Smallest nonzero REAL4 */
#define LAL_REAL4_EPS 1.19209290e-07 /**< 0.5^(LAL_REAL4_MANT-1), i.e. the difference between 1 and the next resolveable REAL4 */
#define LAL_REAL8_MANT 53 /**< Bits of precision in the mantissa of a REAL8 */
#define LAL_REAL8_MAX 1.7976931348623157e+308 /**< Largest REAL8 */
#define LAL_REAL8_MIN 2.2250738585072014e-308 /**< Smallest nonzero REAL8 */
#define LAL_REAL8_EPS 2.2204460492503131e-16  /**< 0.5^(LAL_REAL8_MANT-1), i.e. the difference between 1 and the next resolveable REAL8 */
/*@}*/

/* <lalLaTeX>

The following are fundamental mathematical constants.  They are mostly
taken from the GNU C \verb@math.h@ header (with the exception of
\verb@LAL_TWOPI@, which was computed using Maple).  All numbers are
dimensionless.

\begin{center}
\begin{tabular}{|llc|}
\hline
Name & Value & Expression \\
\hline
\tt LAL\_E         & 2.7182818284590452353602874713526625 & $e$ \\
\tt LAL\_LOG2E     & 1.4426950408889634073599246810018922 & $\log_2 e$ \\
\tt LAL\_LOG10E    & 0.4342944819032518276511289189166051 & $\log_{10} e$ \\
\tt LAL\_LN2       & 0.6931471805599453094172321214581766 & $\log_e 2$ \\
\tt LAL\_LN10      & 2.3025850929940456840179914546843642 & $\log_e 10$ \\
\tt LAL\_SQRT2     & 1.4142135623730950488016887242096981 & $\sqrt{2}$ \\
\tt LAL\_SQRT1\_2  & 0.7071067811865475244008443621048490 & $1/\sqrt{2}$ \\
\tt LAL\_GAMMA     & 0.5772156649015328606065120900824024 & $\gamma$ \\
\tt LAL\_PI        & 3.1415926535897932384626433832795029 & $\pi$ \\
\tt LAL\_TWOPI     & 6.2831853071795864769252867665590058 & $2\pi$ \\
\tt LAL\_PI\_2     & 1.5707963267948966192313216916397514 & $\pi/2$ \\
\tt LAL\_PI\_4     & 0.7853981633974483096156608458198757 & $\pi/4$ \\
\tt LAL\_1\_PI     & 0.3183098861837906715377675267450287 & $1/\pi$ \\
\tt LAL\_2\_PI     & 0.6366197723675813430755350534900574 & $2/\pi$ \\
\tt LAL\_2\_SQRTPI & 1.1283791670955125738961589031215452 & $2/\sqrt{\pi}$ \\
\tt LAL\_PI\_180   & 1.7453292519943295769236907684886127$\times10^{-2}$ &
  $\pi/180$ \\
\tt LAL\_180\_PI   & 57.295779513082320876798154814105170 & $180/\pi$ \\
\hline
\end{tabular}
\end{center}

</lalLaTeX> */

/** \name Mathematical constants
 * The following are fundamental mathematical constants.  They are mostly
 * taken from the GNU C <tt>math.h</tt> header (with the exception of
 * <tt>LAL_TWOPI</tt>, which was computed using Maple).  All numbers are
 * dimensionless. The value of exp(gamma) is taken from
 * http://www.research.att.com/~njas/sequences/A073004 */
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

/* <lalLaTeX>

\subsection*{Physical constants}
\idx[Constant]{LAL\_C\_SI}
\idx[Constant]{LAL\_EPSILON0\_SI}
\idx[Constant]{LAL\_MU0\_SI}
\idx[Constant]{LAL\_GEARTH\_SI}
\idx[Constant]{LAL\_PATM\_SI}
\idx[Constant]{LAL\_YRJUL\_SI}
\idx[Constant]{LAL\_LYR\_SI}
\idx[Constant]{LAL\_G\_SI}
\idx[Constant]{LAL\_H\_SI}
\idx[Constant]{LAL\_HBAR\_SI}
\idx[Constant]{LAL\_MPL\_SI}
\idx[Constant]{LAL\_LPL\_SI}
\idx[Constant]{LAL\_TPL\_SI}
\idx[Constant]{LAL\_K\_SI}
\idx[Constant]{LAL\_R\_SI}
\idx[Constant]{LAL\_MOL}
\idx[Constant]{LAL\_BWIEN\_SI}
\idx[Constant]{LAL\_SIGMA\_SI}
\idx[Constant]{LAL\_AMU\_SI}
\idx[Constant]{LAL\_MP\_SI}
\idx[Constant]{LAL\_ME\_SI}
\idx[Constant]{LAL\_QP\_SI}
\idx[Constant]{LAL\_ALPHA}
\idx[Constant]{LAL\_RE\_SI}
\idx[Constant]{LAL\_LAMBDAE\_SI}
\idx[Constant]{LAL\_AB\_SI}
\idx[Constant]{LAL\_MUB\_SI}
\idx[Constant]{LAL\_MUN\_SI}

The following physical constants are defined to have exact values.
The values of $c$ and $g$ are taken from~\cite{Barnet:1996},
$p_\mathrm{atm}$ is from~\cite{Lang:1992}, while $\epsilon_0$ and
$\mu_0$ are computed from $c$ using exact formulae.  The use of a
Julian year (365.25 days) as standard is specified by the IAU.
They are given in the SI units shown.

\begin{center}
\begin{tabular}{|lll|}
\hline
Name & Value & Description \\
\hline
\tt LAL\_C\_SI        & $299\,792\,458\,\mathrm{m}\,\mathrm{s}^{-1}$ &
	Speed of light $c$ in free space \\
\tt LAL\_EPSILON0\_SI & \multicolumn{2}{l|}{
	$8.8541878176203898505365630317107503\times10^{-12}\,
	\mathrm{C}^2\mathrm{N}^{-1}\mathrm{m}^{-2}$} \\
& & Permittivity $\epsilon_0$ of free space \\
\tt LAL\_MU0\_SI      & \multicolumn{2}{l|}{
	$1.2566370614359172953850573533118012\times10^{-6}\,
	\mathrm{N}\,\mathrm{A}^{-2}$} \\
& & Permeability $\mu_0$ of free space \\
\tt LAL\_GEARTH\_SI   & $9.80665\,\mathrm{m}\,\mathrm{s}^{-2}$ &
	Standard gravity $g$ \\
\tt LAL\_PATM\_SI     & $101\,325\,\mathrm{Pa}$ &
	Standard atmospheric pressure $p_\mathrm{atm}$ \\
\tt LAL\_YRJUL\_SI    & $31\,557\,600\,\mathrm{s}$ &
	(Julian) year \\
\tt LAL\_LYR\_SI    & $9.4607304725808\times10^{15}\,\mathrm{m}$ &
	$c\times$(Julian) year\\
\hline
\end{tabular}
\end{center}

</lalLaTeX> */

/** \name Exact physical constants
 * The following physical constants are defined to have exact values.
 * The values of \f$c\f$ and \f$g\f$ are taken from Barnet (1996),
 * \f$p_\mathrm{atm}\f$ is from Lang (1992), while \f$\epsilon_0\f$ and
 * \f$\mu_0\f$ are computed from \f$c\f$ using exact formulae.  The use
 * of a Julian year (365.25 days) as standard is specified by the IAU.
 * They are given in the SI units shown. */
/*@{*/
#define LAL_C_SI      299792458 /**< Speed of light in vacuo, m s^-1 */
#define LAL_EPSILON0_SI  8.8541878176203898505365630317107503e-12 /**< Permittivity of free space, C^2 N^-1 m^-2 */
#define LAL_MU0_SI    1.2566370614359172953850573533118012e-6 /**< Permeability of free space, N A^-2 */
#define LAL_GEARTH_SI 9.80665 /**< Standard gravity, m s^-2 */
#define LAL_PATM_SI 101325 /**< Standard atmosphere, Pa */
#define LAL_YRJUL_SI 31557600 /**< Julian year, s */
#define LAL_LYR_SI 9.4607304725808e15 /**< (Julian) Lightyear, m */
/*@}*/

/* <lalLaTeX>

The following are measured fundamental physical constants, with values
given in \cite{Barnet:1996}.  When not dimensionless, they are given
in the SI units shown.

\begin{center}
\begin{tabular}{|lll|}
\hline
Name & Value & Description \\
\hline
\tt LAL\_G\_SI     & $6.67259\times10^{-11}\,\mathrm{N}\,\mathrm{m}^{2}
	\mathrm{kg}^{-2}$ & Gravitational constant $G$ \\
\tt LAL\_H\_SI     & $6.6260755\times10^{-34}\,\mathrm{J}\,\mathrm{s}$ &
	Planck constant $h$ \\
\tt LAL\_HBAR\_SI  & $1.05457266\times10^{-34}\,\mathrm{J}\,\mathrm{s}$ &
	Reduced Planck constant $\hbar$ \\
\tt LAL\_MPL\_SI   & $2.17671\times10^{-8}\,\mathrm{kg}$ & Planck mass \\
\tt LAL\_LPL\_SI   & $1.61605\times10^{-35}\,\mathrm{m}$ & Planck length \\
\tt LAL\_TPL\_SI   & $5.39056\times10^{-44}\,\mathrm{s}$ & Planck time \\
\tt LAL\_K\_SI     & $1.380658\times10^{-23}\,\mathrm{J}\,\mathrm{K}^{-1}$ &
	Boltzmann constant $k$ \\
\tt LAL\_R\_SI     & $8.314511\,\mathrm{J}\,\mathrm{K}^{-1}$ &
	Ideal gas constant $R$ \\
\tt LAL\_MOL       & $6.0221367\times10^{23}$ & Avogadro constant \\
\tt LAL\_BWIEN\_SI & $2.897756\times10^{-3}\,\mathrm{m}\,\mathrm{K}$ &
	Wien displacement law constant $b$ \\
\tt LAL\_SIGMA\_SI & $5.67051\times10^{-8}\,\mathrm{W}\,\mathrm{m}^{-2}
	\mathrm{K}^{-4}$ & Stefan-Boltzmann constant $\sigma$ \\
\tt LAL\_AMU\_SI   & $1.6605402\times10^{-27}\,\mathrm{kg}$ &
	Atomic mass unit \\
\tt LAL\_MP\_SI    & $1.6726231\times10^{-27}\,\mathrm{kg}$ & Proton mass \\
\tt LAL\_ME\_SI    & $9.1093897\times10^{-31}\,\mathrm{kg}$ & Electron mass \\
\tt LAL\_QP\_SI    & $1.60217733\times10^{-19}\,\mathrm{C}$ & Proton charge \\
\tt LAL\_ALPHA     & $7.297354677\times10^{-3}$ & Fine structure constant \\
\tt LAL\_RE\_SI    & $2.81794092\times10^{-15}\,\mathrm{m}$ &
	Classical electron radius $r_e$ \\
\tt LAL\_LAMBDAE\_SI & $3.86159323\times10^{-13}\,\mathrm{m}$ &
	Electron Compton wavelength $\lambda_e$ \\
\tt LAL\_AB\_SI    & $5.29177249\times10^{-11}\,\mathrm{m}$ & Bohr radius $a$\\
\tt LAL\_MUB\_SI   & $9.27401543\times10^{-24}\,\mathrm{J}\,\mathrm{T}^{-1}$ &
	Bohr magneton $\mu_B$ \\
\tt LAL\_MUN\_SI   & $5.05078658\times10^{-27}\,\mathrm{J}\,\mathrm{T}^{-1}$ &
	Nuclear magneton $\mu_N$ \\
\hline
\end{tabular}
\end{center}
</lalLaTeX> */

/** \name Physical constants
 * The following are measured fundamental physical constants, with values
 * given in Barnet (1996).  When not dimensionless, they are given
 * in the SI units shown. */
/*@{*/
#define LAL_G_SI      6.67259e-11    /**< Gravitational constant, N m^2 kg^-2 */
#define LAL_H_SI      6.6260755e-34  /**< Planck constant, J s */
#define LAL_HBAR_SI   1.05457266e-34 /**< Reduced Planck constant, J s */
#define LAL_MPL_SI    2.17671e-8     /**< Planck mass, kg */
#define LAL_LPL_SI    1.61605e-35    /**< Planck length, m */
#define LAL_TPL_SI    5.39056e-44    /**< Planck time, s */
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

/* <lalLaTeX>

\subsection*{Astrophysical parameters}
\idx[Constant]{LAL\_REARTH\_SI}
\idx[Constant]{LAL\_AWGS84\_SI}
\idx[Constant]{LAL\_BWGS84\_SI}
\idx[Constant]{LAL\_MEARTH\_SI}
\idx[Constant]{LAL\_IEARTH}
\idx[Constant]{LAL\_EEARTH}
\idx[Constant]{LAL\_RSUN\_SI}
\idx[Constant]{LAL\_MSUN\_SI}
\idx[Constant]{LAL\_MRSUN\_SI}
\idx[Constant]{LAL\_MTSUN\_SI}
\idx[Constant]{LAL\_LSUN\_SI}
\idx[Constant]{LAL\_AU\_SI}
\idx[Constant]{LAL\_PC\_SI}
\idx[Constant]{LAL\_YRTROP\_SI}
\idx[Constant]{LAL\_YRSID\_SI}
\idx[Constant]{LAL\_DAYSID\_SI}
\idx[Constant]{LAL\_H0\_SI}
\idx[Constant]{LAL\_H0FAC\_SI}
\idx[Constant]{LAL\_RHOC\_SI}
\idx[Constant]{LAL\_RHOCFAC\_SI}
\idx[Constant]{LAL\_TCBR\_SI}
\idx[Constant]{LAL\_VCBR\_SI}
\idx[Constant]{LAL\_RHOCBR\_SI}
\idx[Constant]{LAL\_NCBR\_SI}
\idx[Constant]{LAL\_SCBR\_SI}

The following parameters are derived from measured properties of the
Earth and Sun.  The values are taken from~\cite{Barnet:1996}, except
for the obliquity of the ecliptic plane and the eccentricity of
Earth's orbit, which are taken from~\cite{Lang:1992}.  All values are
given in the SI units shown.  Note that the ``year'' and ``light-year''
have exactly defined values, and appear under ``Exact physical constants''.

\begin{center}
\begin{tabular}{|lll|}
\hline
Name & Value & Description \\
\hline
\tt LAL\_REARTH\_SI & $6.378140\times10^6\,\mathrm{m}$ &
	Earth equatorial radius \\
\tt LAL\_AWGS84\_SI & $6.378137\times10^6\,\mathrm{m}$ &
        Semimajor axis of WGS-84 Reference Ellipsoid \\
\tt LAL\_BWGS84\_SI & $6.356752314\times10^6\,\mathrm{m}$ &
        Semiminor axis of WGS-84 Reference Ellipsoid \\
\tt LAL\_MEARTH\_SI & $5.97370\times10^{24}\,\mathrm{kg}$ & Earth mass \\
\tt LAL\_IEARTH     & $0.409092804\,\mathrm{rad}$ &
	Obliquity of the ecliptic (2000) \\
\tt LAL\_EEARTH     & 0.0167 & Earth orbital eccentricity \\
\tt LAL\_RSUN\_SI   & $6.960\times10^8\,\mathrm{m}$ & Solar equatorial radius\\
\tt LAL\_MSUN\_SI   & $1.98892\times10^{30}\,\mathrm{kg}$ & Solar mass \\
\tt LAL\_MRSUN\_SI  & $1.47662504\times10^3\,\mathrm{m}$ &
	Geometrized solar mass (length) \\
\tt LAL\_MTSUN\_SI  & $4.92549095\times10^{-6}\,\mathrm{s}$ &
	Geometrized solar mass (time) \\
\tt LAL\_LSUN\_SI   & $3.846\times10^{26}\,\mathrm{W}$ & Solar luminosity \\
\tt LAL\_AU\_SI     & $1.4959787066\times10^{11}\,\mathrm{m}$ &
	Astronomical unit \\
\tt LAL\_PC\_SI     & $3.0856775807\times10^{16}\,\mathrm{m}$ & Parsec \\
\tt LAL\_YRTROP\_SI & $31\,556\,925.2\,\mathrm{s}$ & Tropical year (1994) \\
\tt LAL\_YRSID\_SI  & $31\,558\,149.8\,\mathrm{s}$ & Sidereal year (1994) \\
\tt LAL\_DAYSID\_SI & $86\,164.09053\,\mathrm{s}$ & Mean sidereal day \\
\hline
\end{tabular}
\end{center}

</lalLaTeX> */

/** \name Astrophysical parameters
 * The following parameters are derived from measured properties of the
 * Earth and Sun.  The values are taken from Barnet (1996), except
 * for the obliquity of the ecliptic plane and the eccentricity of
 * Earth's orbit, which are taken from Lang (1992).  All values are
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
#define LAL_EEARTH    0.0167          /**< Earth orbital eccentricity */
#define LAL_RSUN_SI   6.960e8         /**< Solar equatorial radius, m */
#define LAL_MSUN_SI   1.98892e30      /**< Solar mass, kg */
#define LAL_MRSUN_SI  1.47662504e3    /**< Geometrized solar mass, m */
#define LAL_MTSUN_SI  4.92549095e-6   /**< Geometrized solar mass, s */
#define LAL_LSUN_SI   3.846e26        /**< Solar luminosity, W */
#define LAL_AU_SI     1.4959787066e11 /**< Astronomical unit, m */
#define LAL_PC_SI     3.0856775807e16 /**< Parsec, m */
#define LAL_YRTROP_SI 31556925.2      /**< Tropical year (1994), s */
#define LAL_YRSID_SI  31558149.8      /**< Sidereal year (1994), s */
#define LAL_DAYSID_SI 86164.09053     /**< Mean sidereal day, s */
/*@}*/

/* <lalLaTeX>

The following cosmological parameters are derived from measurements of
the Hubble expansion rate and of the cosmic background radiation
(CBR).  Data are taken from~\cite{Barnet:1996}.  In what follows, the
normalized Hubble constant $h_0$ is equal to the actual Hubble
constant $H_0$ divided by $\langle H
\rangle=100\,\mathrm{km}\,\mathrm{s}^{-1}\mathrm{Mpc}^{-1}$.  Thus the
Hubble constant can be written as:
$$
H_0 = \langle H \rangle h_0 \; .
$$
Similarly, the critical energy density $\rho_c$ required for spatial
flatness is given by:
$$
\rho_c = \langle\rho\rangle h_0^2 \; .
$$
Current estimates give $h_0$ a value of around 0.65, which is what is
assumed below.  All values are in the SI units shown.

\begin{center}
\begin{tabular}{|lll|}
\hline
Name & Value & Description \\
\hline
\tt LAL\_H0\_SI      & $2\times10^{-18}\,\mathrm{s}^{-1}$ &
	Approx.\ Hubble constant $H_0$ \\
\tt LAL\_H0FAC\_SI   & $3.2407792903\times10^{-18}\,\mathrm{s}^{-1}$ &
	$H_0/h_0$ \\
\tt LAL\_RHOC\_SI    & $7\times10^{-10}\,\mathrm{J}\,\mathrm{m}^{-3}$ &
	Approx.\ critical energy density $\rho_c$ \\
\tt LAL\_RHOCFAC\_SI & $1.68860\times10^{-9}\,\mathrm{J}\,\mathrm{m}^{-3}$ &
	$\rho_c/h_0^2$ \\
\tt LAL\_TCBR\_SI    & $2.726 \mathrm{K}$ &
	CBR temperature \\
\tt LAL\_VCBR\_SI    & $3.695\times10^5\,\mathrm{m}\,\mathrm{s}^{-1}$ &
	Solar velocity with respect to CBR \\
\tt LAL\_RHOCBR\_SI  & $4.177\times10^{-14}\,\mathrm{J}\,\mathrm{m}^{-3}$ &
	Energy density of CBR \\
\tt LAL\_NCBR\_SI    & $4.109\times10^8\,\mathrm{m}^{-3}$ &
	Number density of CBR photons \\
\tt LAL\_SCBR\_SI    & $3.993\times10^{-14}\,\mathrm{J}\,\mathrm{K}^{-1}
	\mathrm{m}^{-3}$ & Entropy density of CBR \\
\hline
\end{tabular}
\end{center}

</lalLaTeX> */

/** \name Cosmological parameters
 * The following cosmological parameters are derived from measurements of
 * the Hubble expansion rate and of the cosmic background radiation
 * (CBR).  Data are taken from Barnet (1996).  In what follows, the
 * normalized Hubble constant \f$h_0\f$ is equal to the actual Hubble
 * constant \f$H_0\f$ divided by \f$\langle H
 * \rangle=100\,\mathrm{km}\,\mathrm{s}^{-1}\mathrm{Mpc}^{-1}\f$.  Thus the
 * Hubble constant can be written as:
 * \f$H_0 = \langle H \rangle h_0\f$.
 * Similarly, the critical energy density \f$\rho_c\f$ required for spatial
 * flatness is given by: \f$\rho_c = \langle\rho\rangle h_0^2\f$.
 * Current estimates give \f$h_0\f$ a value of around 0.65, which is what is
 * assumed below.  All values are in the SI units shown. */
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


/* <lalLaTeX>

\vfill{\footnotesize\input{LALConstantsHV}}

</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _LALCONSTANTS_H */
