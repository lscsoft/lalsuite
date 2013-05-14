/*
*  Copyright (C) 2007 Jolien Creighton, John Whelan
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

#include <lal/LALStdlib.h>
#include <string.h>
#include <ctype.h>
#include <lal/Units.h>

#define UNITDEFSC_TEMPSIZE 20


/**
\author J. T. Whelan <john.whelan@ligo.org>
\addtogroup UnitDefs_c

\brief Defines basic and derived SI units and a function to produce a text
string corresponding to a unit structure.

LALUnitAsString() converts the unit structure
<tt>*input</tt> into a text string which is stored in the character
vector <tt>*output</tt>.  Note that the resulting text string is
expressed solely in terms of the basic units (m, kg, s, A,
K, strain and counts), and is thus not necessarily the most
convenient way to check the units of a quantity.  A better method is
to construct a unit structure containing the expected units, then
compare that to the actual units using LALUnitCompare().

LALParseUnitString() reconstructs the original
\c LALUnit structure from the string output by
LALUnitAsString().  It is very sensitive to the exact format
of the string and is not intended for use in parsing user-entered
strings.

\heading{Algorithm}

LALUnitAsString() moves through the unit structure, appending
the appropriate text to the string as it goes along.

LALParseUnitString() moves through the input string, one
character at a time, building an ::LALUnit structure as a it
goes along, so long as it encounters precisely the syntax expected.

\heading{Notes}

This file also defines a number of \c constant unit structures
(declared \c extern in \ref Units_h).  Zeroth is
\c lalDimensionlessUnit, which is simply a ::LALUnit
structure to be associated with a unitless quantity.
First, the relevant fundamental SI units and two custom units of use in
gravitational wave detection:

<table><tr><th>Constant</th><th>Name</th><th>Abbr.</th><th>Physical Quantity</th></tr>
<tr><td>\c #lalMeterUnit</td><td>meter</td><td>m</td><td>length</td></tr>
<tr><td>\c #lalKiloGramUnit</td><td>kilogram</td><td>kg</td><td>mass</td></tr>
<tr><td>\c #lalSecondUnit</td><td>second</td><td>s</td><td>time</td></tr>
<tr><td>\c #lalAmpereUnit</td><td>ampere</td><td>A</td><td>electric current</td></tr>
<tr><td>\c #lalKelvinUnit</td><td>kelvin</td><td>K</td><td>thermodynamic temperature</td></tr>
<tr><td>\c #lalStrainUnit</td><td>strain</td><td>\f$\epsilon\f$</td><td>gravitational strain</td></tr>
<tr><td>\c #lalADCCountUnit</td><td>ADC count</td><td>count</td><td>A-to-D converter counts</td></tr>
</table>

Next, the named derived units in the SI [\ref Halliday_2001]:

<table><tr><th>Constant</th><th>Name</th><th>Abbr.</th><th>Physical Quantity</th><th>Def.</th><th>Fundamental</th></tr>
<tr><td>\c #lalHertzUnit</td><td>hertz</td><td>Hz</td><td>frequency</td><td> s\f$^{-1}\f$</td><td>s\f$^{-1}\f$</td></tr>
<tr><td>\c #lalNewtonUnit</td><td>newton</td><td>N</td><td>force</td><td> kg\f$\cdot\f$ m/s\f$^2\f$</td><td>m kg s\f$^{-2}\f$</td></tr>
<tr><td>\c #lalPascalUnit</td><td>pascal</td><td>Pa</td><td>pressure</td><td> N/m\f$^2\f$</td><td>m\f$^{-1}\f$ kg s\f$^{-2}\f$</td></tr>
<tr><td>\c #lalJouleUnit</td><td>joule</td><td>J</td><td>energy</td><td> N\f$\cdot\f$m</td><td>m\f$^2\f$ kg s\f$^{-2}\f$</td></tr>
<tr><td>\c #lalWattUnit</td><td>watt</td><td>W</td><td>power</td><td> J/s</td><td>m\f$^2\f$ kg s\f$^{-3}\f$</td></tr>
<tr><td>\c #lalCoulombUnit</td><td>coulomb</td><td>C</td><td>electric charge</td><td>A\f$\cdot\f$s</td><td>s A</td></tr>
<tr><td>\c #lalVoltUnit</td><td>volt</td><td>V</td><td>potential</td><td> W/A</td><td>m\f$^2\f$ kg s\f$^{-3}\f$ A\f$^{-1}\f$</td></tr>
<tr><td>\c #lalOhmUnit</td><td>ohm</td><td>\f$\Omega\f$</td><td>resistance</td><td> V/A</td><td>m\f$^2\f$ kg s\f$^{-3}\f$ A\f$^{-2}\f$</td></tr>
<tr><td>\c #lalFaradUnit</td><td>farad</td><td>F</td><td>capacitance</td><td> C/V</td><td>m\f$^{-2}\f$ kg\f$^{-1}\f$ s\f$^4\f$ A\f$^2\f$</td></tr>
<tr><td>\c #lalWeberUnit</td><td>weber</td><td>Wb</td><td>magnetic flux</td><td> V\f$\cdot\f$s</td><td>m\f$^2\f$ kg s\f$^{-2}\f$ A\f$^{-1}\f$</td></tr>
<tr><td>\c #lalHenryUnit</td><td>henry</td><td>H</td><td>inductance</td><td> V\f$\cdot\f$s/A</td><td>m\f$^2\f$ kg s\f$^{-2}\f$ A\f$^{-2}\f$</td></tr>
<tr><td>\c #lalTeslaUnit</td><td>tesla</td><td>T</td><td>magnetic flux density</td><td> Wb/m\f$^2\f$</td><td>kg s\f$^{-2}\f$ A\f$^{-1}\f$</td></tr>
</table>

The powers of ten (SI prefixes)

<table><tr><th>Constant</th><th>Prefix</th><th>Abbr.</th><th>Value</th></tr>
<tr><td>\c #lalYottaUnit</td><td>yotta</td><td>Y</td><td>\f$10^{ 24}\f$</td></tr>
<tr><td>\c #lalZettaUnit</td><td>zetta</td><td>Z</td><td>\f$10^{ 21}\f$</td></tr>
<tr><td>\c #lalExaUnit</td><td>exa</td><td>E</td><td>\f$10^{ 18}\f$</td></tr>
<tr><td>\c #lalPetaUnit</td><td>peta</td><td>P</td><td>\f$10^{ 15}\f$</td></tr>
<tr><td>\c #lalTeraUnit</td><td>tera</td><td>T</td><td>\f$10^{ 12}\f$</td></tr>
<tr><td>\c #lalGigaUnit</td><td>giga</td><td>G</td><td>\f$10^{  9}\f$</td></tr>
<tr><td>\c #lalMegaUnit</td><td>mega</td><td>M</td><td>\f$10^{  6}\f$</td></tr>
<tr><td>\c #lalKiloUnit</td><td>kilo</td><td>k</td><td>\f$10^{  3}\f$</td></tr>
<tr><td>\c #lalHectoUnit</td><td>hecto</td><td>h</td><td>\f$10^{  2}\f$</td></tr>
<tr><td>\c #lalDekaUnit</td><td>deka</td><td>da</td><td>\f$10^{  1}\f$</td></tr>
<tr><td>\c #lalDeciUnit</td><td>deci</td><td>d</td><td>\f$10^{ -1}\f$</td></tr>
<tr><td>\c #lalCentiUnit</td><td>centi</td><td>c</td><td>\f$10^{ -2}\f$</td></tr>
<tr><td>\c #lalMilliUnit</td><td>milli</td><td>m</td><td>\f$10^{ -3}\f$</td></tr>
<tr><td>\c #lalMicroUnit</td><td>micro</td><td>\f$\mu\f$</td><td>\f$10^{ -6}\f$</td></tr>
<tr><td>\c #lalNanoUnit</td><td>nano</td><td>n</td><td>\f$10^{ -9}\f$</td></tr>
<tr><td>\c #lalPicoUnit</td><td>pico</td><td>p</td><td>\f$10^{-12}\f$</td></tr>
<tr><td>\c #lalFemtoUnit</td><td>femto</td><td>f</td><td>\f$10^{-15}\f$</td></tr>
<tr><td>\c #lalAttoUnit</td><td>atto</td><td>a</td><td>\f$10^{-18}\f$</td></tr>
<tr><td>\c #lalZeptoUnit</td><td>zepto</td><td>z</td><td>\f$10^{-21}\f$</td></tr>
<tr><td>\c #lalYoctoUnit</td><td>yocto</td><td>y</td><td>\f$10^{-24}\f$</td></tr>
</table>

And finally a couple of convenient scaled units:

<table><tr><th>Constant</th><th>Name</th><th>Abbr.\</th><th>Def.</th><th>Fundamental</th></tr>
<tr><td>\c #lalGramUnit</td><td>gram</td><td>g</td><td> \f$10^{-3}\f$ kg</td><td>\f$10^{-3}\f$ kg</td></tr>
<tr><td>\c #lalAttoStrainUnit</td><td>attostrain</td><td>a\f$\epsilon\f$</td><td>\f$10^{-18} \epsilon\f$</td><td>\f$10^{-18} \epsilon\f$</td></tr>
<tr><td>\c #lalPicoFaradUnit</td><td>picofarad</td><td>pF</td><td>\f$10^{-12}\f$ F</td><td>\f$10^{-12}\f$ m\f$^{-2}\f$ kg\f$^{-1}\f$ s\f$^4\f$ A\f$^2\f$</td></tr>
</table>

*/
/*@{*/

/** To convert a units structure to a string repesentation, we need to
 * define the names of the basic units.
 */
const CHAR lalUnitName[LALNumUnits][LALUnitNameSize] =
{
  "m", "kg", "s", "A", "K", "strain", "count"
};

/* ********************************************************
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

const LALUnit lalDimensionlessUnit = {  0, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< dimensionless units */

/** \name Basic Units */
/*@{*/
const LALUnit lalMeterUnit         = {  0, { 1, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< meter [m] */
const LALUnit lalKiloGramUnit      = {  0, { 0, 1, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< kilogram [kg]*/
const LALUnit lalSecondUnit        = {  0, { 0, 0, 1, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< second [s] */
const LALUnit lalAmpereUnit        = {  0, { 0, 0, 0, 1, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Ampere [A] */
const LALUnit lalKelvinUnit        = {  0, { 0, 0, 0, 0, 1, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Kelvin [K] */
const LALUnit lalStrainUnit        = {  0, { 0, 0, 0, 0, 0, 1, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Strain [1] */
const LALUnit lalADCCountUnit      = {  0, { 0, 0, 0, 0, 0, 0, 1}, { 0, 0, 0, 0, 0, 0, 0} };	/**< ADC count [count] */
/*@}*/

/** \name Derived Mechanical Units */
/*@{*/
const LALUnit lalHertzUnit         = {  0, { 0, 0,-1, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Hertz [Hz] */
const LALUnit lalNewtonUnit        = {  0, { 1, 1,-2, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Newton [N] */
const LALUnit lalPascalUnit        = {  0, {-1, 1,-2, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Pascal [Pa] */
const LALUnit lalJouleUnit         = {  0, { 2, 1,-2, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Joule [J] */
const LALUnit lalWattUnit          = {  0, { 2, 1,-3, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Watt [W ] */
/*@}*/

/** \name Derived Electromagnetic Units */
/*@{*/
const LALUnit lalCoulombUnit       = {  0, { 0, 0, 1, 1, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Coulomb [C] */
const LALUnit lalVoltUnit          = {  0, { 2, 1,-3,-1, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Volt [V] */
const LALUnit lalOhmUnit           = {  0, { 2, 1,-3,-2, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Ohm [\f$\Omega\f$] */
const LALUnit lalFaradUnit         = {  0, {-2,-1, 4, 2, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Farad [F] */
const LALUnit lalWeberUnit         = {  0, { 2, 1,-2,-1, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Weber [Wb] */
const LALUnit lalHenryUnit         = {  0, { 2, 1,-2,-2, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Henry [H] */
const LALUnit lalTeslaUnit         = {  0, { 0, 1,-2,-1, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Tesla [T] */
/*@}*/

/** \name Powers of Ten */
/*@{*/
const LALUnit lalYottaUnit         = { 24, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Yotta [1e24] */
const LALUnit lalZettaUnit         = { 21, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Zetta [1e21] */
const LALUnit lalExaUnit           = { 18, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Exa [1e18] */
const LALUnit lalPetaUnit          = { 15, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Peta [1e15] */
const LALUnit lalTeraUnit          = { 12, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Tera [1e12] */
const LALUnit lalGigaUnit          = {  9, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Giga [1e9] */
const LALUnit lalMegaUnit          = {  6, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Mega [1e6] */
const LALUnit lalKiloUnit          = {  3, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Kilo [1e3] */
const LALUnit lalHectoUnit         = {  2, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Hecto [1e2] */
const LALUnit lalDekaUnit          = {  1, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Deka [1e1] */
const LALUnit lalDeciUnit          = { -1, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Deci [1e-1] */
const LALUnit lalCentiUnit         = { -2, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Centi [1e-2] */
const LALUnit lalMilliUnit         = { -3, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Milli [1e-3] */
const LALUnit lalMicroUnit         = { -6, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Micro [1e-6] */
const LALUnit lalNanoUnit          = { -9, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Nano [1e-9] */
const LALUnit lalPicoUnit          = {-12, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Pico [1e-12] */
const LALUnit lalFemtoUnit         = {-15, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Femto [1e-15] */
const LALUnit lalAttoUnit          = {-18, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Atto [1e-18] */
const LALUnit lalZeptoUnit         = {-21, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Zepto [1e-21] */
const LALUnit lalYoctoUnit         = {-24, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Yocto [1e-24] */
/*@}*/

/** \name Convenient Scaled Units */
/*@{*/
const LALUnit lalGramUnit          = { -3, { 0, 1, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< Gram [1e-3] */
const LALUnit lalAttoStrainUnit    = {-18, { 0, 0, 0, 0, 0, 1, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< AttoStrain [1e-18] */
const LALUnit lalPicoFaradUnit     = {-12, {-2,-1, 2, 2, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };	/**< PicoFarad [1e-12 F] */
/*@}*/

/* Static function to read a number into a character array */
/* returns 0 on success, 1 on failure */
/* leaves *charPtrPtr pointing to first non-digit */
static int readNumber( char temp[], const char **charPtrPtr )
{
  CHAR *tempPtr, *tempStopPtr;

  tempPtr = temp;
  /* make sure we don't fall off end of temporary array */
  tempStopPtr = temp + UNITDEFSC_TEMPSIZE;

  if ( ! isdigit(**charPtrPtr) ) return 1;

  do
  {
    *tempPtr = **charPtrPtr;
    ++tempPtr, ++*charPtrPtr;
    if (tempPtr >= tempStopPtr) return 1;
  }
  while ( isdigit(**charPtrPtr) );
  *tempPtr = '\0';
  return 0;
}

/* Static function to read a string into a character array */
/* returns 0 on success, 1 on failure */
/* leaves *charPtrPtr pointing to first non-letter */
static int readString( char temp[UNITDEFSC_TEMPSIZE], const char **charPtrPtr )
{
  CHAR *tempPtr, *tempStopPtr;

  tempPtr = temp;
  /* make sure we don't fall off end of temporary array */
  tempStopPtr = temp + UNITDEFSC_TEMPSIZE;

  if ( ! isalpha(**charPtrPtr) ) return 1;

  do
  {
    *tempPtr = **charPtrPtr;
    ++tempPtr, ++*charPtrPtr;
    if (tempPtr >= tempStopPtr) return 1;
  }
  while ( isalpha(**charPtrPtr) );
  *tempPtr = '\0';
  return 0;
}

/** Returns the pointer to the input \c string, which
 * is populated with the unit string if successful.  If there is a failure,
 * XLALUnitAsString() returns a \c NULL pointer and \c xlalErrno
 * is set to one of the following values:  \c XLAL_EFAULT if one of the
 * input pointers is \c NULL or \c XLAL_EBADLEN if the length of the
 * string is insufficent for the unit string.
 */
char * XLALUnitAsString( char *string, UINT4 length, const LALUnit *input )
{
  UINT2        i;
  CHAR         temp[UNITDEFSC_TEMPSIZE];
  INT2         numer;
  CHAR         *charPtr, *charStopPtr;

  if ( ! string || ! input )
    XLAL_ERROR_NULL( XLAL_EFAULT );
  if ( ! length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );

  charPtr = string;
  charStopPtr = string + length;
  /* Points one past the end of the array! */

  *charPtr = '\0';

  if (input->powerOfTen != 0)
  {
    sprintf(temp, "10^%d", input->powerOfTen);
    if ( charPtr + strlen(temp) >= charStopPtr)
      XLAL_ERROR_NULL( XLAL_EBADLEN );
    strncpy(charPtr, temp, charStopPtr - charPtr);
    charPtr += strlen(temp);
  } /* if (input->powerOfTen != 0) */

  for (i=0; i<LALNumUnits; ++i)
  {
    numer = input->unitNumerator[i];
    if (numer != 0)
    {
      if (charPtr != string)
      {
	*charPtr = ' ';
	++charPtr;
      } /* if (charPtr != output->data) */
      if (input->unitDenominatorMinusOne[i] == 0)
      {
	if (numer == 1)
	{
	  sprintf(temp, "%s", lalUnitName[i]);
	} /* if (numer == 1) */
	else
	{
	  sprintf(temp, "%s^%d", lalUnitName[i], numer);
	}
      } /* if (input->unitDenominatorMinusOne[i] == 0) */
      else {
	sprintf(temp, "%s^%d/%d", lalUnitName[i], numer,
		 input->unitDenominatorMinusOne[i] + 1);
      }
      if ( charPtr + strlen(temp) >= charStopPtr)
        XLAL_ERROR_NULL( XLAL_EBADLEN );
      strncpy(charPtr, temp, charStopPtr - charPtr);
      charPtr += strlen(temp);
    } /* if (numer != 0) */
  }  /* for (i=0; i<LALNumUnits; ++i) */

  return string;
}

/** Allocates and returns a new string, which is populated with the unit
 * string.  If there is a failure, returns a \c NULL pointer and \c xlalErrno
 * is set to one of the error values of \c XLALUnitAsString or \c XLALMalloc.
 * Caller is responsible for freeing return value with \c XLALFree.
 */
char * XLALUnitToString( const LALUnit *input )
{
  char *output = NULL, *buf = XLALMalloc(LALUnitTextSize);
  if (buf)
  {
    output = XLALUnitAsString(buf, LALUnitTextSize, input);
    if (output)
      output = XLALRealloc(buf, strlen(buf) + 1);
    if (!output)
      XLALFree(buf);
  }
  return output;
}

/** \deprecated Use XLALUnitAsString() instead.
 */
void
LALUnitAsString( LALStatus *status,
		 CHARVector *output,
		 const LALUnit *input )

{
  INITSTATUS(status);
  /* ATTATCHSTATUSPTR (status); */

  ASSERT( input != NULL, status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );

  ASSERT( output != NULL, status, UNITSH_ENULLPOUT, UNITSH_MSGENULLPOUT );

  ASSERT( output->data != NULL, status, UNITSH_ENULLPD,
	  UNITSH_MSGENULLPD );

  ASSERT( output->length > 0, status,
	  UNITSH_ESTRINGSIZE, UNITSH_MSGESTRINGSIZE );

  if ( ! XLALUnitAsString( output->data, output->length, input ) )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EFAULT: /* a NULL pointer was passed to the XLAL function */
        if ( ! input )
        {
          ABORT( status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );
        }
        else /* must have been a NULL output data pointer */
        {
          ABORT( status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );
        }
      default: /* otherwise must have almost overwritten the string */
        ABORT( status, UNITSH_ESTRINGSIZE, UNITSH_MSGESTRINGSIZE );
    }
  }

  RETURN(status);
}

/** Returns the pointer \c output upon return
 * or a pointer to newly allocated memory if \c output was \c NULL;
 * on failure, \c XLALParseUnitString() returns \c NULL and sets
 * ::xlalErrno to one of the following values: \c #XLAL_ENOMEM
 * if the routine was unable to allocate memory for the output or
 * \c #XLAL_EFAILED if the routine was unable to parse the unit string.
 */
LALUnit * XLALParseUnitString( LALUnit *output, const char *string )
{
  UINT2        i;
  INT2         sign;
  CHAR         temp[20];
  int outputAllocated = 0;

  if ( ! output )
  {
    output = LALMalloc( sizeof( *output ) );
    outputAllocated = 1;
    if ( ! output )
      XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  /* Start with dimensionless (all zeros) and fill in from there */
  *output = lalDimensionlessUnit;

  /* If the string is NULL, it represents dimensionless */
  if ( ! string )
    return output;

  /* Strip leading whitespace */
  string += strspn(string, "\t\n\v\f\r ");

  /* If the string is empty, it represents dimensionless */
  if ( ! *string )
    return output;

  /* Look for power of ten; note LALUnitsAsString is set up to say
   * "10^1" rather than "10", so need not allow for missing '^'
   */
  if (*string == '1' && *(string+1) == '0' && *(string+2) == '^')
  {
    string += 3;
    /* now pointing at first digit of power of ten (or minus sign) */

    if ( *string == '-'  )
    {
      sign = -1;
      ++string;
    }
    else
    {
      sign = 1;
    }

    /* read power of ten into temp[]; return value of 1 means failure */
    if ( readNumber( temp, &string ) )
    {
      if ( outputAllocated )
        LALFree( output );
      XLAL_ERROR_NULL( XLAL_EFAILED );
    }
    /* string now points to one after end of power of ten */

    output->powerOfTen = sign*atoi(temp);

    /* If the power of ten was all there was, return */
    if (*string == '\0')
      return output;

    if ( *string != ' ')
    {
      if ( outputAllocated )
        LALFree( output );
      XLAL_ERROR_NULL( XLAL_EFAILED );
    }

    ++string;
  } /* if (*string == '1' && *(string+1) == '0' && *(string+2) == '^') */

  /* string now points to start of first unit */

  /* Read units and exponents, one unit per pass of the following do loop */
  do
  {
    /* read unit name into temp[]; return value of 1 means failure */
    if ( readString( temp, &string ) )
    {
      if ( outputAllocated )
        LALFree( output );
      XLAL_ERROR_NULL( XLAL_EFAILED );
    }

    /* string now points to one after end of unit name */

    /* find which unit name this matches */
    for (i=0; strcmp(temp,lalUnitName[i]); ++i)
    {
      if (i>=LALNumUnits) /* didn't find it */
      {
        if ( outputAllocated )
          LALFree( output );
        XLAL_ERROR_NULL( XLAL_EFAILED );
      }
    }

    /* Make sure we haven't already read in this unit */
    if ( output->unitNumerator[i] || output->unitDenominatorMinusOne[i] )
    {
      if ( outputAllocated )
        LALFree( output );
      XLAL_ERROR_NULL( XLAL_EFAILED );
    }

    if ( *string == ' ' || *string == '\0' )
    { /* We have only one power of the unit */
      output->unitNumerator[i] = 1;
    }
    else if ( *string == '^' )
    {
      ++string;
      /* now points to the first digit of the exponent, or minus sign */

      if ( *string == '-'  )
      {
	sign = -1;
	++string;
      }
      else
      {
	sign = 1;
      }

      /* read exponent numerator into temp[];
	 return value of 1 means failure */
      if ( readNumber( temp, &string ) )
      {
        if ( outputAllocated )
          LALFree( output );
        XLAL_ERROR_NULL( XLAL_EFAILED );
      }
      output->unitNumerator[i] = sign * atoi(temp);

      if ( *string == '/' )
      {
	++string;
	/* now points to first digit of denominator */

	/* read exponent denominator into temp[];
	   return value of 1 means failure */
	if ( readNumber( temp, &string ) || temp[0] == '0')
        {
          if ( outputAllocated )
            LALFree( output );
          XLAL_ERROR_NULL( XLAL_EFAILED );
        }
	output->unitDenominatorMinusOne[i] = atoi(temp) - 1;
      } /* if ( *string == '/' ) */
    } /* else if ( *string == '^' ) */
    else
    {
      if ( outputAllocated )
        LALFree( output );
      XLAL_ERROR_NULL( XLAL_EFAILED );
    }

    if ( *string == ' ') ++string;

  }
  while ( *string != '\0' );

  return output;
}


/** \deprecated Use XLALParseUnitString() instead.
 */
void
LALParseUnitString ( LALStatus *status,
		     LALUnit *output,
		     const CHARVector *input )

{
  CHAR         *charPtr, *charStopPtr;

  INITSTATUS(status);

  ASSERT( input != NULL, status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );

  ASSERT( output != NULL, status, UNITSH_ENULLPOUT, UNITSH_MSGENULLPOUT );

  ASSERT( input->data != NULL, status, UNITSH_ENULLPD, UNITSH_MSGENULLPD );

  /* ensure that there's a '\0' within the input CHARVector */
  charPtr = input->data;
  charStopPtr = charPtr + strlen(input->data);
  /* Should point to first '\0' in string */
  if (charStopPtr >= charPtr + input->length)
  {
    ABORT( status, UNITSH_EPARSE, UNITSH_MSGEPARSE );
  }

  /* call the XLAL function */
  output = XLALParseUnitString( output, charPtr );
  if ( ! output ) /* there was a parse error */
  {
    XLALClearErrno();
    ABORT( status, UNITSH_EPARSE, UNITSH_MSGEPARSE );
  }

  RETURN(status);
}
/*@}*//* end: UnitDefs_c */
