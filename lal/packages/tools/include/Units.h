/********************************* <lalVerbatim file="UnitsHV">
Author: Whelan, J. T.
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
\index{\texttt{LALUnitPair}}
Consists of a pair of unit structures; used as an input structure for
the \texttt{LALUnitCompare()} and \texttt{LALUnitMultiply()} functions.
The fields are:
\begin{description}
\item[\texttt{LALUnit unitOne}] The first unit.
\item[\texttt{LALUnit unitTwo}] The second unit.
\end{description}

\subsubsection*{RAT4}
\index{\texttt{RAT4}}
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


/*********************************************************
 *                                                       *
 *       Functions to manipulate unit structures         *
 *                                                       *
 *********************************************************/

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
  LALUnit   unitOne;
  LALUnit   unitTwo;
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


/* The parameter structure for LALUnitRaise contains the numerator and
 * denominator-minus-one of the rational power.  
 */

typedef struct
tagRAT4 
{
  INT2 numerator;
  UINT2 denominatorMinusOne;
} RAT4;

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

