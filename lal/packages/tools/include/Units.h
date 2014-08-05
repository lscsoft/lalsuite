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

#ifndef _UNITS_H
#define _UNITS_H

#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \addtogroup Units_h
 * \author J. T. Whelan <john.whelan@ligo.org>
 *
 * \brief Provides prototypes for manipulation of units and declares
 * \c extern constants for the basic and derived SI units.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/Units.h>
 * \endcode
 *
 * This header provides prototypes for functions to manipulate
 * the \c LALUnit structure.  It also defines \c extern
 * constants for a set of predefined units, which are designed to make
 * the structure easier to use.  For instance, to determine whether a
 * quantity has units of strain per root hertz, one constructs the unit
 * &quot;strain per root hertz&quot; from the predefined \c lalStrainUnit and
 * \c lalHertzUnit constant structures using the
 * LALUnitRaise() and LALUnitMultiply() functions, then
 * compares that to the unit structure in question using the
 * LALUnitCompare() function.
 *
 * The LALUnit datatype itself is included in the header
 * \ref LALDatatypes.h, and defines a unit in terms of an integer
 * power of ten multiplier along with rational powers of the basic SI
 * units (meters, kilograms, seconds, Amperes, and Kelvins) and two
 * custom units (strain and ADC counts).
 *
 * ### XLAL interface to Units.h functions ###
 *
 * XLALUnitAsString() converts a ::LALUnit structure into a character
 * string of maximum length \c length (including NUL termination)
 * representation of the units.  The inverse function, XLALParseUnitString()
 * parses a character string to produce a \c LALUnit structure; if
 * \c output is \c NULL, memory for the output is allocated.  If the input
 * \c string is \c NULL or is empty then the output units are
 * dimensionless: lalDimensionlessUnit.
 *
 * XLALUnitNormalize() puts a ::LALUnit structure into normal form
 * by simplifying all unit exponent fractions to their simplest form.
 *
 * XLALUnitCompare() compares two ::LALUnit structures: they are the
 * same if their normal forms are identical.
 *
 * XLALUnitMultiply() multiplies two ::LALUnit structures.  The result
 * is put into normal form.
 *
 * XLALUnitRaiseRAT4() raises a ::LALUnit structure to a rational
 * power given by the ::RAT4 structure \c power.
 * XLALUnitRaiseINT2() raises a ::LALUnit structure to an integer
 * power \c power.
 * XLALUnitSquare() produces the square of a ::LALUnit structure.
 * XLALUnitSqrt() produces the square-root of a ::LALUnit structure.
 *
 * ### Return Values ###
 *
 * XLALUnitAsString() returns the pointer to the input \c string, which
 * is populated with the unit string if successful.  If there is a failure,
 * XLALUnitAsString() returns a \c NULL pointer and \c ::xlalErrno
 * is set to one of the following values:  \c #XLAL_EFAULT if one of the
 * input pointers is \c NULL or \c #XLAL_EBADLEN if the length of the
 * string is insufficent for the unit string.
 *
 * XLALParseUnitString() returns the pointer \c output upon return
 * or a pointer to newly allocated memory if \c output was \c NULL;
 * on failure, \c XLALParseUnitString returns \c NULL and sets
 * \c ::xlalErrno to one of the following values:  \c #XLAL_ENOMEM
 * if the routine was unable to allocate memory for the output or
 * \c #XLAL_EFAILED if the routine was unable to parse the unit string.
 *
 * XLALUnitNormalize() returns 0 upon success or \c #XLAL_FAILURE
 * if the input pointer is \c NULL, in which case \c xlalErrno
 * is set to \c #XLAL_EFAULT
 *
 * XLALUnitCompare() returns 0 if the the normal form of the two unit
 * structures are the same or \> 0 if they are different.  It returns
 * \c #XLAL_FAILURE and \c ::xlalErrno is set to \c #XLAL_EFAULT if
 * one of the input pointers is \c NULL.
 *
 * XLALUnitMultiply(), XLALUnitRaiseRAT4(), XLALUnitRaiseINT2(), XLALUnitSquare() and
 * XLALUnitSqrt() all return a pointer to the output unit structure
 * \c output upon success or \c NULL upon failure.  If there is
 * a failure, \c ::xlalErrno is set to one of the following values:
 * \c #XLAL_EFAULT if one of the input pointers is \c NULL,
 * \c #XLAL_ERANGE if one of the unit powers exceeds the allowed range,
 * or \c #XLAL_EINVAL (for the raise functions only) if the unit power
 * would not be an integer.
 *
 * @{
 * \defgroup UnitDefs_c 		Module UnitDefs.c
 * \defgroup UnitNormalize_c 	Module UnitNormalize.c
 * \defgroup UnitRaise_c 		Module UnitRaise.c
 * \defgroup UnitMultiply_c 	Module UnitMultiply.c
 * \defgroup UnitCompare_c 	Module UnitCompare.c
 * @}
 */
/*@{*/

/** \name Error Codes */
/*@{*/
#define UNITSH_ENULLPIN         1	/**< Null pointer to input */
#define UNITSH_ENULLPOUT        2	/**< Null pointer to output */
#define UNITSH_ENULLPD          3	/**< Null pointer to data member of vector */
#define UNITSH_ENULLPPARAM      4	/**< Null pointer to parameters */
#define UNITSH_ESTRINGSIZE      5	/**< Output string too short */
#define UNITSH_EOVERFLOW        6	/**< Exponent outside of (U)INT2 bounds */
#define UNITSH_ENONINT          7	/**< Non-integer power of ten */
#define UNITSH_EPARSE           8	/**< Error parsing unit string */
/*@}*/

/** \cond DONT_DOXYGEN */
#define UNITSH_MSGENULLPIN      "Null pointer to input"
#define UNITSH_MSGENULLPOUT     "Null pointer to output"
#define UNITSH_MSGENULLPD       "Null pointer to data member of vector"
#define UNITSH_MSGENULLPPARAM   "Null pointer to parameters"
#define UNITSH_MSGESTRINGSIZE   "Output string too short"
#define UNITSH_MSGEOVERFLOW     "Exponent outside of (U)INT2 bounds"
#define UNITSH_MSGENONINT       "Non-integer power of ten"
#define UNITSH_MSGEPARSE        "Error parsing unit string"
/** \endcond */

/**
 * A four-byte rational number, used as a parameter structure for
 * LALUnitRaise().
 */
typedef struct
tagRAT4
{
  INT2 numerator;		/**< The numerator */
  UINT2 denominatorMinusOne;	/**< One less than the denominator */
} RAT4;

/**
 * Consists of a pair of unit structures; used as an input structure for
 * the LALUnitCompare() and LALUnitMultiply() functions.
 */
#ifdef SWIG /* SWIG interface directives */
SWIGLAL(IMMUTABLE_MEMBERS(tagLALUnitPair, unitOne, unitTwo));
#endif /* SWIG */
typedef struct
tagLALUnitPair
{
  const LALUnit   *unitOne;	/**< The first unit */
  const LALUnit   *unitTwo;	/**< The second unit */
}
LALUnitPair;

/*@}*/
/* end: Units_h */

/*********************************************************
 *                                                       *
 *       Functions to manipulate unit structures         *
 *                                                       *
 *********************************************************/

#ifndef SWIG /* exclude from SWIG interface */

/* XLAL routines */
char * XLALUnitAsString( char *string, UINT4 length, const LALUnit *input );
char * XLALUnitToString( const LALUnit *input );
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

/* Obsolete LAL-interface functions */
void LALUnitNormalize (LALStatus *status, LALUnit *output, const LALUnit *input);
void LALUnitMultiply (LALStatus *status, LALUnit *output, const LALUnitPair *input);
void LALUnitCompare (LALStatus *status, BOOLEAN *output, const LALUnitPair *input);
void LALUnitRaise (LALStatus *status, LALUnit *output, const LALUnit *input, const RAT4 *power);
void LALUnitAsString (LALStatus *status, CHARVector *output, const LALUnit *input);
void LALParseUnitString ( LALStatus *status, LALUnit *output, const CHARVector *input );

enum enumLALUnitNameSize {
  LALUnitNameSize = sizeof("strain")
};
enum enumLALUnitTextSize {
  LALUnitTextSize = sizeof("10^-32768 m^-32768/32767 kg^-32768/32767 "
                           "s^-32768/32767 A^-32768/32767 "
                           "K^-32768/32767 strain^-32768/32767 "
                           "count^-32768/32767")
};

extern const CHAR lalUnitName[LALNumUnits][LALUnitNameSize];

#endif /* SWIG */

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
