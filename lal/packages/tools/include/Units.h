/*
*  Copyright (C) 2007 Jolien Creighton, Kipp Cannon, John Whelan
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

/********************************* <lalVerbatim file="UnitsHV">
Author: J. T. Whelan <john.whelan@ligo.org>
$Id$
********************************** </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{Units.h}}
\label{tools:s:Units.h}

Provides prototypes for manipulation of units and declares
\texttt{extern} constants for the basic and derived SI units.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/Units.h>
\end{verbatim}

\noindent This header provides prototypes for functions to manipulate
the \texttt{LALUnit} structure.  It also defines \texttt{extern}
constants for a set of predefined units, which are designed to make
the structure easier to use.  For instance, to determine whether a
quantity has units of strain per root hertz, one constructs the unit
"strain per root hertz" from the predefined \texttt{lalStrainUnit} and
\texttt{lalHertzUnit} constant structures using the
\texttt{LALUnitRaise()} and \texttt{LALUnitMultiply()} functions, then
compares that to the unit structure in question using the
\texttt{LALUnitCompare()} function.

The LALUnit datatype itself is included in the header
\texttt{LALDatatypes.h}, and defines a unit in terms of an integer
power of ten multiplier along with rational powers of the basic SI
units (meters, kilograms, seconds, Amperes, and Kelvins) and two
custom units (strain and ADC counts).

\subsection*{Error conditions}
\input{UnitsHErrTable}

\subsection*{Structures}

\subsubsection*{LALUnitPair}
\idx[Type]{LALUnitPair}
Consists of a pair of unit structures; used as an input structure for
the \texttt{LALUnitCompare()} and \texttt{LALUnitMultiply()} functions.
The fields are:
\begin{description}
\item[\texttt{LALUnit *unitOne}] The first unit.
\item[\texttt{LALUnit *unitTwo}] The second unit.
\end{description}

\subsubsection*{RAT4}
\idx[Type]{RAT4}
A four-byte rational number, used as a parameter structure for
\texttt{LALUnitRaise()}.  The fields are:
\begin{description}
\item[\texttt{INT2 numerator}] The numerator.
\item[\texttt{UINT2 denominatorMinusOne}] One less than the denominator.
\end{description}

\vfill{\footnotesize\input{UnitsHV}}

\newpage\input{UnitDefsC}
\newpage\input{UnitNormalizeC}
\newpage\input{UnitRaiseC}
\newpage\input{UnitMultiplyC}
\newpage\input{UnitCompareC}

\newpage\subsection{XLAL Functions}

\subsubsection*{Synopsis}
\begin{verbatim}
#include <lal/Units.h>

char * XLALUnitAsString( char *string, UINT4 length, const LALUnit *input );
LALUnit * XLALParseUnitString( LALUnit *output, const char *string );
int XLALUnitNormalize( LALUnit *unit );
int XLALUnitCompare( const LALUnit *unit1, const LALUnit *unit2 );
LALUnit * XLALUnitMultiply( LALUnit *output, const LALUnit *unit1,
    const LALUnit *unit2 );
LALUnit * XLALUnitRaiseRAT4( LALUnit *output, const LALUnit *input,
    const RAT4 *power );
LALUnit * XLALUnitRaiseINT2( LALUnit *output, const LALUnit *input,
    INT2 power );
LALUnit * XLALUnitSquare( LALUnit *output, const LALUnit *input );
LALUnit * XLALUnitSqrt( LALUnit *output, const LALUnit *input );
\end{verbatim}
\idx{XLALUnitAsString}
\idx{XLALParseUnitString}
\idx{XLALUnitNormalize}
\idx{XLALUnitCompare}
\idx{XLALUnitMultiply}
\idx{XLALUnitRaiseRAT4}
\idx{XLALUnitRaiseINT2}
\idx{XLALUnitSquare}
\idx{XLALUnitSqrt}

\subsubsection*{Description}

\verb+XLALUnitAsString+ converts a \verb+LALUnit+ structure into a character
string of maximum length \verb+length+ (including NUL termination)
representation of the units.  The inverse function, \verb+XLALParseUnitString+
parses a character string to produce a \verb+LALUnit+ structure; if
\verb+output+ is \verb+NULL+, memory for the output is allocated.  If the input
\verb+string+ is \verb+NULL+ or is empty then the output units are
dimensionless: \verb+lalDimensionlessUnit+.

\verb+XLALUnitNormalize+ puts a \verb+LALUnit+ structure into normal form
by simplifying all unit exponent fractions to their simplest form.

\verb+XLALUnitCompare+ compares two \verb+LALUnit+ structures: they are the
same if their normal forms are identical.

\verb+XLALUnitMultiply+ multiplies two \verb+LALUnit+ structures.  The result
is put into normal form.

\verb+XLALUnitRaiseRAT4+ raises a \verb+LALUnit+ structure to a rational
power given by the \verb+RAT4+ structure \verb+power+.
\verb+XLALUnitRaiseINT2+ raises a \verb+LALUnit+ structure to an integer
power \verb+power+.
\verb+XLALUnitSquare+ produces the square of a \verb+LALUnit+ structure.
\verb+XLALUnitSqrt+ produces the square-root of a \verb+LALUnit+ structure.

\subsubsection*{Return Values}

\verb+XLALUnitAsString+ returns the pointer to the input \verb+string+, which
is populated with the unit string if successful.  If there is a failure,
\verb+XLALUnitAsString+ returns a \verb+NULL+ pointer and \verb+xlalErrno+
is set to one of the following values:  \verb+XLAL_EFAULT+ if one of the
input pointers is \verb+NULL+ or \verb+XLAL_EBADLEN+ if the length of the
string is insufficent for the unit string.

\verb+XLALParseUnitString+ returns the pointer \verb+output+ upon return
or a pointer to newly allocated memory if \verb+output+ was \verb+NULL+;
on failure, \verb+XLALParseUnitString+ returns \verb+NULL+ and sets
\verb+xlalErrno+ to one of the following values:  \verb+XLAL_ENOMEM+
if the routine was unable to allocate memory for the output or
\verb+XLAL_EFAILED+ if the routine was unable to parse the unit string.

\verb+XLALUnitNormalize+ returns 0 upon success or \verb+XLAL_FAILURE+
if the input pointer is \verb+NULL+, in which case \verb+xlalErrno+
is set to \verb+XLAL_EFAULT+

\verb+XLALUnitCompare+ returns 0 if the the normal form of the two unit
structures are the same or > 0 if they are different.  It returns
\verb+XLAL_FAILURE+ and \verb+xlalErrno+ is set to \verb+XLAL_EFAULT+ if
one of the input pointers is \verb+NULL+.

\verb+XLALUnitMultiply+
\verb+XLALUnitRaiseRAT4+
\verb+XLALUnitRaiseINT2+
\verb+XLALUnitSquare+ and
\verb+XLALUnitSqrt+ all return a pointer to the output unit structure
\verb+output+ upon success or \verb+NULL+ upon failure.  If there is
a failure, \verb+xlalErrno+ is set to one of the following values:
\verb+XLAL_EFAULT+ if one of the input pointers is \verb+NULL+,
\verb+XLAL_ERANGE+ if one of the unit powers exceeds the allowed range,
or \verb+XLAL_EINVAL+ (for the raise functions only) if the unit power
would not be an integer.

\newpage\input{UnitsTestC}

</lalLaTeX> */

#ifndef _UNITS_H
#define _UNITS_H

#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (UNITSH, "$Id$");

/******************************** <lalErrTable file="UnitsHErrTable"> */

#define UNITSH_ENULLPIN         1
#define UNITSH_ENULLPOUT        2
#define UNITSH_ENULLPD          3
#define UNITSH_ENULLPPARAM      4
#define UNITSH_ESTRINGSIZE      5
#define UNITSH_EOVERFLOW        6
#define UNITSH_ENONINT          7
#define UNITSH_EPARSE           8

#define UNITSH_MSGENULLPIN      "Null pointer to input"
#define UNITSH_MSGENULLPOUT     "Null pointer to output"
#define UNITSH_MSGENULLPD       "Null pointer to data member of vector"
#define UNITSH_MSGENULLPPARAM   "Null pointer to parameters"
#define UNITSH_MSGESTRINGSIZE   "Output string too short"
#define UNITSH_MSGEOVERFLOW     "Exponent outside of (U)INT2 bounds"
#define UNITSH_MSGENONINT       "Non-integer power of ten"
#define UNITSH_MSGEPARSE        "Error parsing unit string"

/************************************ </lalErrTable> */


/* The parameter structure for LALUnitRaise contains the numerator and
 * denominator-minus-one of the rational power.
 */

typedef struct
tagRAT4
{
  INT2 numerator;
  UINT2 denominatorMinusOne;
} RAT4;


/*********************************************************
 *                                                       *
 *       Functions to manipulate unit structures         *
 *                                                       *
 *********************************************************/


/* XLAL routines */

char * XLALUnitAsString( char *string, UINT4 length, const LALUnit *input );
LALUnit * XLALParseUnitString( LALUnit *output, const char *string );
int XLALUnitNormalize( LALUnit *unit );
int XLALUnitCompare( const LALUnit *unit1, const LALUnit *unit2 );
LALUnit * XLALUnitMultiply( LALUnit *output, const LALUnit *unit1,
    const LALUnit *unit2 );
LALUnit * XLALUnitDivide( LALUnit *output, const LALUnit *unit1,
    const LALUnit *unit2 );
LALUnit * XLALUnitRaiseRAT4( LALUnit *output, const LALUnit *input,
    const RAT4 *power );
LALUnit * XLALUnitRaiseINT2( LALUnit *output, const LALUnit *input,
    INT2 power );
LALUnit * XLALUnitSquare( LALUnit *output, const LALUnit *input );
LALUnit * XLALUnitSqrt( LALUnit *output, const LALUnit *input );
LALUnit * XLALUnitInvert( LALUnit *output, const LALUnit *input );
REAL8 XLALUnitPrefactor(const LALUnit *unit);
int XLALUnitIsDimensionless(const LALUnit *unit);
REAL8 XLALUnitRatio(const LALUnit *unit1, const LALUnit *unit2);



/* LALUnitNormalize will reduce the rational powers in the basic unit
 * exponents, e.g. converting 2/2 to 1/1 and 3/6 to 1/2.
 */

void LALUnitNormalize (LALStatus *status, LALUnit *output,
		       const LALUnit *input);

/* Several functions take two Unit variables as input, so need the
 * following structure.
 */

typedef struct
tagLALUnitPair
{
  const LALUnit   *unitOne;
  const LALUnit   *unitTwo;
}
LALUnitPair;

/* LALUnitMultiply will multiply together two Unit variables and
 * output their product; it will call LALUnitNormalize.
 */

void LALUnitMultiply (LALStatus *status, LALUnit *output,
		      const LALUnitPair *input);

/* LALUnitCompare will compare two Unit variables and return true if
 * they are the equivalent (the same power of ten offset as well as
 * equivalent ratioanl powers of the base units).
 */

void LALUnitCompare (LALStatus *status, BOOLEAN *output,
		      const LALUnitPair *input);

/* LALUnitRaise will raise a unit structure to a rational power; the
 * most common choices will presumably be -1, 2, and 1/2.  An error
 * occurs if input->powerOfTen is not evenly divisible by the
 * denominator of the specified power
 */

void LALUnitRaise (LALStatus *status, LALUnit *output,
		   const LALUnit *input, const RAT4 *power);

/* LALUnitAsString will convert an LALUnit structure into a
 * human-readable text form.
 */

void LALUnitAsString (LALStatus *status, CHARVector *output,
		      const LALUnit *input);

void
LALParseUnitString ( LALStatus *status,
		     LALUnit *output,
		     const CHARVector *input );

enum { LALUnitNameSize = sizeof("strain") };
enum { LALUnitTextSize = sizeof("10^-32768 m^-32768/32767 kg^-32768/32767 "
				"s^-32768/32767 A^-32768/32767 "
				"K^-32768/32767 strain^-32768/32767 "
				"count^-32768/32767") };

extern const CHAR lalUnitName[LALNumUnits][LALUnitNameSize];

/*********************************************************
 *                                                       *
 *                 Predefined units                      *
 *                                                       *
 *********************************************************/

/* Predefined constant units make it easier for programmers to specify
 * and compare (using LALUnitCompare) units more easily.  Those given
 * here are an example; more can be added.
 */

/* LALUnitsTest.c will verify the definitions of the derived units,
 * for example using LALUnitRaise, LALUnitMultiply and LALUnitCompare
 * to show that 1 Farad = 1 Coulomb Volt^-1
 */

extern const LALUnit lalDimensionlessUnit;

/* Basic Units */
extern const LALUnit lalMeterUnit    ;
extern const LALUnit lalKiloGramUnit ;
extern const LALUnit lalSecondUnit   ;
extern const LALUnit lalAmpereUnit   ;
extern const LALUnit lalKelvinUnit   ;
extern const LALUnit lalStrainUnit   ;
extern const LALUnit lalADCCountUnit ;

/* Derived Mechanical Units */
extern const LALUnit lalHertzUnit     ;
extern const LALUnit lalNewtonUnit    ;
extern const LALUnit lalJouleUnit     ;
extern const LALUnit lalWattUnit      ;
extern const LALUnit lalPascalUnit    ;

/* Derived Electromagnetic Units */
extern const LALUnit lalCoulombUnit   ;
extern const LALUnit lalVoltUnit      ;
extern const LALUnit lalOhmUnit       ;
extern const LALUnit lalFaradUnit     ;
extern const LALUnit lalWeberUnit     ;
extern const LALUnit lalTeslaUnit     ;
extern const LALUnit lalHenryUnit     ;

/* Powers of Ten */
extern const LALUnit lalYottaUnit     ;
extern const LALUnit lalZettaUnit     ;
extern const LALUnit lalExaUnit       ;
extern const LALUnit lalPetaUnit      ;
extern const LALUnit lalTeraUnit      ;
extern const LALUnit lalGigaUnit      ;
extern const LALUnit lalMegaUnit      ;
extern const LALUnit lalKiloUnit      ;
extern const LALUnit lalHectoUnit     ;
extern const LALUnit lalDekaUnit      ;
extern const LALUnit lalDeciUnit      ;
extern const LALUnit lalCentiUnit     ;
extern const LALUnit lalMilliUnit     ;
extern const LALUnit lalMicroUnit     ;
extern const LALUnit lalNanoUnit      ;
extern const LALUnit lalPicoUnit      ;
extern const LALUnit lalFemtoUnit     ;
extern const LALUnit lalAttoUnit      ;
extern const LALUnit lalZeptoUnit     ;
extern const LALUnit lalYoctoUnit     ;

/* Convenient Scaled Units */
extern const LALUnit lalGramUnit      ;
extern const LALUnit lalAttoStrainUnit;
extern const LALUnit lalPicoFaradUnit ;


/*********************************************************
 *                                                       *
 *    Functions to manipulate quantities with units      *
 *                                                       *
 *********************************************************/

/* The various RescaleUnits routines will change the power of ten
 * offset in the unit structure of a structured datatype to the
 * specified value and multiply each element of the structure by the
 * appropriate factor, so that for instance {15, 0, .24} x 10^-3 kg
 * becomes {.015, 0, .00024} kg if newPowerOfTen is 0.
 */
/* There will be routines for each type of series; examples are given
 * below.
 */
/*
void LALSFRescaleUnits (LALStatus *status, REAL4FrequencySeries *output,
			const REAL4FrequencySeries *input,
			const INT2 *newPowerOfTen);

void LALCTRescaleUnits (LALStatus *status, COMPLEX8TimeSeries *output,
			const COMPLEX8TimeSeries *input,
			const INT2 *newPowerOfTen);

void LALI2TVRescaleUnits (LALStatus *status, INT2TimeVectorSeries *output,
			  const INT2TimeVectorSeries *input,
			  const INT2 *newPowerOfTen);

void LALU2TARescaleUnits (LALStatus *status, UINT2TimeArraySeries *output,
			  const UINT2TimeArraySeries *input,
			  const INT2 *newPowerOfTen);
*/

#ifdef  __cplusplus
}
#endif

#endif /* _UNITS_H */

