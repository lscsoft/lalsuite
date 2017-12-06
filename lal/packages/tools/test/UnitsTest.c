/*
*  Copyright (C) 2014 Karl Wette (XLALification)
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

/**
 * \author J. T. Whelan <john.whelan@ligo.org>
 * \file
 * \ingroup Units_h
 *
 * \brief Test Suite for unit manipulation programs
 *
 * ### Description ###
 *
 * This program tests the various units-manipulation routines, as well as
 * the pre-defined units defined in \ref Units.h.  For each
 * successful test, it prints "PASS" to standard output.
 *
 * ### Uses ###
 *
 * \code
 * XLALMalloc()
 * XLALUnitAsString()
 * XLALParseUnitString()
 * XLALFree()
 * XLALUnitMultiply()
 * XLALUnitRaiseRAT4()
 * XLALUnitNormalize()
 * XLALUnitCompare()
 * \endcode
 *
 */

/** \cond  DONT_DOXYGEN */

#include <config.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>

int main(void) {

  int errnum;
  LALUnit unit1, unit2, unit3, *punit;
  RAT4 power;
  CHAR *string, buffer[LALUnitTextSize];
  UINT4 length;

  printf("Checking input validation of XLALUnitAsString:\n");

  length = sizeof("m^2 kg s^-3 A^-1") - 1;
  string = XLALMalloc(length);
  XLAL_CHECK( string != NULL, XLAL_EFAILED );

  XLAL_TRY( XLALUnitAsString( NULL, length, &unit1 ), errnum );
  XLAL_CHECK( errnum == XLAL_EFAULT, XLAL_EFAILED );
  printf("  PASS: Null pointer to output string\n");

  XLAL_TRY( XLALUnitAsString( string, 0, &unit1 ), errnum );
  XLAL_CHECK( errnum == XLAL_EBADLEN, XLAL_EFAILED );
  printf("  PASS: Bad output string length\n");

  XLAL_TRY( XLALUnitAsString( string, length, NULL ), errnum );
  XLAL_CHECK( errnum == XLAL_EFAULT, XLAL_EFAILED );
  printf("  PASS: Null pointer to input LALUnit\n");

  XLAL_TRY( XLALUnitAsString( string, length, &lalVoltUnit ), errnum );
  XLAL_CHECK( errnum == XLAL_EBADLEN, XLAL_EFAILED );
  printf("  PASS: Output string too short\n");

  XLALFree(string);

  printf("Testing response of XLALUnitAsString to valid data:\n");

  length = sizeof("m^2 kg s^-3 A^-1");
  string = XLALMalloc(length);
  XLAL_CHECK( string != NULL, XLAL_EFAILED );

  XLAL_CHECK( XLALUnitAsString( string, length, &lalVoltUnit ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( strcmp(string, "m^2 kg s^-3 A^-1") == 0, XLAL_EFAILED );
  printf("  PASS: String representation of Volt is '%s'\n", string);

  XLALFree(string);

  length = sizeof("10^-12 m^-2 kg^-1 s^2 A^2");
  string = XLALMalloc(length);
  XLAL_CHECK( string != NULL, XLAL_EFAILED );

  XLAL_CHECK( XLALUnitAsString( string, length, &lalPicoFaradUnit ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( strcmp(string, "10^-12 m^-2 kg^-1 s^2 A^2") == 0, XLAL_EFAILED );
  printf("  PASS: String representation of PicoFarad is '%s'\n", string);

  XLAL_CHECK( XLALUnitAsString( string, length, &lalMeterUnit ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( strcmp(string, "m") == 0, XLAL_EFAILED );
  printf("  PASS: String representation of Meter is '%s'\n", string);

  XLAL_CHECK( XLALUnitAsString( string, length, &lalCoulombUnit ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( strcmp(string, "s A") == 0, XLAL_EFAILED );
  printf("  PASS: String representation of Coulomb is '%s'\n", string);

  XLAL_CHECK( XLALUnitAsString( string, length, &lalMegaUnit ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( strcmp(string, "10^6") == 0, XLAL_EFAILED );
  printf("  PASS: String representation of Mega is '%s'\n", string);

  XLAL_CHECK( XLALUnitAsString( string, length, &lalAttoStrainUnit ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( strcmp(string, "10^-18 strain") == 0, XLAL_EFAILED );
  printf("  PASS: String representation of AttoStrain is '%s'\n", string);

  XLAL_CHECK( XLALUnitAsString( string, length, &lalDimensionlessUnit ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( strcmp(string, "") == 0, XLAL_EFAILED );
  printf("  PASS: String representation of dimensionless is '%s'\n", string);

  XLALFree(string);

  printf("Checking input validation of XLALParseUnitString:\n");

  punit = XLALParseUnitString( NULL, "" );
  XLAL_CHECK( punit != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( punit, &lalDimensionlessUnit) == 0, XLAL_EFAILED );
  XLALFree(punit);
  printf("  PASS: Null pointer to output LALUnit allocates new one\n");

  punit = XLALParseUnitString( &unit1, NULL );
  XLAL_CHECK( punit == &unit1, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalDimensionlessUnit) == 0, XLAL_EFAILED );
  printf("  PASS: Null pointer to string returns dimensionless\n");

  memset(buffer, 0, LALUnitTextSize);

  strncpy( buffer, "10^", LALUnitTextSize );
  XLAL_TRY( XLALParseUnitString( &unit1, buffer ), errnum );
  XLAL_CHECK( errnum == XLAL_EFAILED, XLAL_EFAILED );
  printf("  PASS: Expected failure to parse '%s'\n", buffer);

  strncpy( buffer, "s^3 s^2", LALUnitTextSize );
  XLAL_TRY( XLALParseUnitString( &unit1, buffer ), errnum );
  XLAL_CHECK( errnum == XLAL_EFAILED, XLAL_EFAILED );
  printf("  PASS: Expected failure to parse '%s'\n", buffer);

  strncpy( buffer, "V", LALUnitTextSize );
  XLAL_TRY( XLALParseUnitString( &unit1, buffer ), errnum );
  XLAL_CHECK( errnum == XLAL_EFAILED, XLAL_EFAILED );
  printf("  PASS: Expected failure to parse '%s'\n", buffer);

  printf("Testing response of XLALParseUnitString to valid data:\n");

  XLAL_CHECK( XLALUnitAsString( buffer, LALUnitTextSize, &lalVoltUnit ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALParseUnitString( &unit1, buffer ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalVoltUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: Volt <-> string\n");

  XLAL_CHECK( XLALUnitAsString( buffer, LALUnitTextSize, &lalPicoFaradUnit ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALParseUnitString( &unit1, buffer ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalPicoFaradUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: PicoFarad <-> string\n");

  XLAL_CHECK( XLALUnitAsString( buffer, LALUnitTextSize, &lalMeterUnit ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALParseUnitString( &unit1, buffer ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalMeterUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: Meter <-> string\n");

  XLAL_CHECK( XLALUnitAsString( buffer, LALUnitTextSize, &lalAttoStrainUnit ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALParseUnitString( &unit1, buffer ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalAttoStrainUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: AttoStrain <-> string\n");

  XLAL_CHECK( XLALUnitAsString( buffer, LALUnitTextSize, &lalCoulombUnit ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALParseUnitString( &unit1, buffer ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalCoulombUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: Coulomb <-> string\n");

  XLAL_CHECK( XLALUnitAsString( buffer, LALUnitTextSize, &lalMegaUnit ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALParseUnitString( &unit1, buffer ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalMegaUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: Mega <-> string\n");

  XLAL_CHECK( XLALUnitAsString( buffer, LALUnitTextSize, &lalDimensionlessUnit ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALParseUnitString( &unit1, buffer ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalDimensionlessUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: Dimensionless <-> string\n");

  unit1 = unit2 = lalDimensionlessUnit;

  printf("Checking input validation of LALUnitMultiply:\n");

  XLAL_TRY( XLALUnitMultiply( NULL, &unit1, &unit2 ), errnum );
  XLAL_CHECK( errnum == XLAL_EFAULT, XLAL_EFAILED );
  printf("  PASS: Null pointer to output LALUnit\n");

  XLAL_TRY( XLALUnitMultiply( &unit3, NULL, NULL ), errnum );
  XLAL_CHECK( errnum == XLAL_EFAULT, XLAL_EFAILED );
  printf("  PASS: Null pointer to input LALUnits\n");

  unit1.powerOfTen = 20000;
  unit2.powerOfTen = 30000;
  XLAL_TRY( XLALUnitMultiply( &unit3, &unit1, &unit2 ), errnum );
  XLAL_CHECK( errnum == XLAL_ERANGE, XLAL_EFAILED );
  printf("  PASS: Exponent outside of (U)INT2 bounds\n");

  unit1.powerOfTen = 0;
  unit2.powerOfTen = 0;
  unit1.unitNumerator[2] = 12345;
  unit1.unitDenominatorMinusOne[2] = 23456;
  unit2.unitNumerator[2] = 23456;
  unit2.unitDenominatorMinusOne[2] = 12345;
  XLAL_TRY( XLALUnitMultiply( &unit3, &unit1, &unit2 ), errnum );
  XLAL_CHECK( errnum == XLAL_ERANGE, XLAL_EFAILED );
  printf("  PASS: Product outside of (U)INT2 bounds\n");

  printf("Checking input validation of XLALUnitRaiseRAT4:\n");

  XLAL_TRY( XLALUnitRaiseRAT4( NULL, &unit1, &power ), errnum );
  XLAL_CHECK( errnum == XLAL_EFAULT, XLAL_EFAILED );
  printf("  PASS: Null pointer to output LALUnit\n");

  XLAL_TRY( XLALUnitRaiseRAT4( &unit3, NULL, &power ), errnum );
  XLAL_CHECK( errnum == XLAL_EFAULT, XLAL_EFAILED );
  printf("  PASS: Null pointer to input LALUnit\n");

  XLAL_TRY( XLALUnitRaiseRAT4( &unit3, &unit1, NULL ), errnum );
  XLAL_CHECK( errnum == XLAL_EFAULT, XLAL_EFAILED );
  printf("  PASS: Null pointer to input RAT4\n");

  unit1 = lalKiloUnit;
  power.numerator = 1;
  power.denominatorMinusOne = 1;
  XLAL_TRY( XLALUnitRaiseRAT4( &unit2, &unit1, &power ), errnum );
  XLAL_CHECK( errnum == XLAL_EINVAL, XLAL_EFAILED );
  printf("  PASS: Non-integer power of ten\n");

  unit1.unitNumerator[2] = 20000;
  unit1.unitNumerator[2] += 20000;
  power.numerator = 2;
  power.denominatorMinusOne = 0;
  XLAL_TRY( XLALUnitRaiseRAT4( &unit2, &unit1, &power ), errnum );
  XLAL_CHECK( errnum == XLAL_ERANGE, XLAL_EFAILED );
  printf("  PASS: Exponent outside of (U)INT2 bounds\n");

  printf("Testing response of XLALUnitNormalize to valid data:\n");

  unit1 = lalDimensionlessUnit;
  unit1.unitNumerator[0] = 2;
  unit1.unitDenominatorMinusOne[0] = 5;
  XLAL_CHECK( XLALUnitNormalize( &unit1 ) == XLAL_SUCCESS, XLAL_EFAILED );
  unit2 = lalDimensionlessUnit;
  unit2.unitNumerator[0] = 1;
  unit2.unitDenominatorMinusOne[0] = 2;
  XLAL_CHECK( XLALUnitCompare( &unit1, &unit2 ) == 0, XLAL_EFAILED );
  printf("  PASS: 2/6 reduces to 1/3\n");

  printf("Checking input validation of XLALUnitCompare:\n");

  XLAL_TRY( XLALUnitCompare( &unit1, NULL ), errnum );
  XLAL_CHECK( errnum == XLAL_EFAULT, XLAL_EFAILED );
  XLAL_TRY( XLALUnitCompare( NULL, &unit2 ), errnum );
  XLAL_CHECK( errnum == XLAL_EFAULT, XLAL_EFAILED );
  printf("  PASS: Null pointer to input LALUnits\n");

  printf("Testing response of XLALUnitCompare to valid data:\n");

  XLAL_CHECK( XLALUnitCompare( &lalMeterUnit, &lalMeterUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: m = m\n");

  XLAL_CHECK( XLALUnitCompare( &lalMeterUnit, &lalKiloGramUnit ) != 0, XLAL_EFAILED );
  printf("  PASS: m != kg\n");

  XLAL_CHECK( XLALUnitCompare( &lalGramUnit, &lalKiloGramUnit ) != 0, XLAL_EFAILED );
  printf("  PASS: g != kg\n");

  printf("Testing definitions of basic units:\n");

  unit1 = lalDimensionlessUnit;
  unit1.unitNumerator[LALUnitIndexMeter] = 1;
  unit2 = lalMeterUnit;
  XLAL_CHECK( XLALUnitCompare( &unit1, &unit2 ) == 0, XLAL_EFAILED );
  XLAL_CHECK( unit2.unitNumerator[LALUnitIndexMeter] == 1, XLAL_EFAILED );
  printf("  PASS: meter\n");

  unit1 = lalDimensionlessUnit;
  unit1.unitNumerator[LALUnitIndexKiloGram] = 1;
  unit2 = lalKiloGramUnit;
  XLAL_CHECK( XLALUnitCompare( &unit1, &unit2 ) == 0, XLAL_EFAILED );
  XLAL_CHECK( unit2.unitNumerator[LALUnitIndexKiloGram] == 1, XLAL_EFAILED );
  printf("  PASS: kilogram\n");

  unit1 = lalDimensionlessUnit;
  unit1.unitNumerator[LALUnitIndexSecond] = 1;
  unit2 = lalSecondUnit;
  XLAL_CHECK( XLALUnitCompare( &unit1, &unit2 ) == 0, XLAL_EFAILED );
  XLAL_CHECK( unit2.unitNumerator[LALUnitIndexSecond] == 1, XLAL_EFAILED );
  printf("  PASS: second\n");

  unit1 = lalDimensionlessUnit;
  unit1.unitNumerator[LALUnitIndexAmpere] = 1;
  unit2 = lalAmpereUnit;
  XLAL_CHECK( XLALUnitCompare( &unit1, &unit2 ) == 0, XLAL_EFAILED );
  XLAL_CHECK( unit2.unitNumerator[LALUnitIndexAmpere] == 1, XLAL_EFAILED );
  printf("  PASS: Ampere\n");

  unit1 = lalDimensionlessUnit;
  unit1.unitNumerator[LALUnitIndexKelvin] = 1;
  unit2 = lalKelvinUnit;
  XLAL_CHECK( XLALUnitCompare( &unit1, &unit2 ) == 0, XLAL_EFAILED );
  XLAL_CHECK( unit2.unitNumerator[LALUnitIndexKelvin] == 1, XLAL_EFAILED );
  printf("  PASS: Kelvin\n");

  unit1 = lalDimensionlessUnit;
  unit1.unitNumerator[LALUnitIndexStrain] = 1;
  unit2 = lalStrainUnit;
  XLAL_CHECK( XLALUnitCompare( &unit1, &unit2 ) == 0, XLAL_EFAILED );
  XLAL_CHECK( unit2.unitNumerator[LALUnitIndexStrain] == 1, XLAL_EFAILED );
  printf("  PASS: strain\n");

  unit1 = lalDimensionlessUnit;
  unit1.unitNumerator[LALUnitIndexADCCount] = 1;
  unit2 = lalADCCountUnit;
  XLAL_CHECK( XLALUnitCompare( &unit1, &unit2 ) == 0, XLAL_EFAILED );
  XLAL_CHECK( unit2.unitNumerator[LALUnitIndexADCCount] == 1, XLAL_EFAILED );
  printf("  PASS: ADC Count\n");

  printf("Testing definitions of derived units:\n");

  power.numerator = -1;
  power.denominatorMinusOne = 0;
  XLAL_CHECK( XLALUnitRaiseRAT4( &unit1, &lalSecondUnit, &power ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalHertzUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: Hz is s^-1\n");

  power.numerator = -2;
  power.denominatorMinusOne = 0;
  XLAL_CHECK( XLALUnitRaiseRAT4( &unit2, &lalSecondUnit, &power ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitMultiply( &unit1, &lalMeterUnit, &unit2 ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitMultiply( &unit2, &lalKiloGramUnit, &unit1 ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit2, &lalNewtonUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: N is kg m s^-2\n");

  power.numerator = -2;
  power.denominatorMinusOne = 0;
  XLAL_CHECK( XLALUnitRaiseRAT4( &unit2, &lalMeterUnit, &power ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitMultiply( &unit1, &lalNewtonUnit, &unit2 ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalPascalUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: Pa is N m^-2\n");

  XLAL_CHECK( XLALUnitMultiply( &unit1, &lalNewtonUnit, &lalMeterUnit ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalJouleUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: J is N m\n");

  power.numerator = -1;
  power.denominatorMinusOne = 0;
  XLAL_CHECK( XLALUnitRaiseRAT4( &unit2, &lalSecondUnit, &power ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitMultiply( &unit1, &lalJouleUnit, &unit2 ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalWattUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: W is J/s\n");

  XLAL_CHECK( XLALUnitMultiply( &unit1, &lalAmpereUnit, &lalSecondUnit ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalCoulombUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: C is A s\n");

  power.numerator = -1;
  power.denominatorMinusOne = 0;
  XLAL_CHECK( XLALUnitRaiseRAT4( &unit2, &lalAmpereUnit, &power ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitMultiply( &unit1, &lalWattUnit, &unit2 ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalVoltUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: V is W/A\n");

  power.numerator = -1;
  power.denominatorMinusOne = 0;
  XLAL_CHECK( XLALUnitRaiseRAT4( &unit2, &lalAmpereUnit, &power ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitMultiply( &unit1, &lalVoltUnit, &unit2 ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalOhmUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: Ohm is V/A\n");

  power.numerator = -1;
  power.denominatorMinusOne = 0;
  XLAL_CHECK( XLALUnitRaiseRAT4( &unit2, &lalVoltUnit, &power ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitMultiply( &unit1, &lalCoulombUnit, &unit2 ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalFaradUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: F is C/V\n");

  XLAL_CHECK( XLALUnitMultiply( &unit1, &lalVoltUnit, &lalSecondUnit ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalWeberUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: Wb is V s\n");

  power.numerator = -1;
  power.denominatorMinusOne = 0;
  XLAL_CHECK( XLALUnitRaiseRAT4( &unit2, &lalAmpereUnit, &power ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitMultiply( &unit1, &lalSecondUnit, &unit2 ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitMultiply( &unit1, &lalVoltUnit, &unit1 ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalHenryUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: H is V s/A\n");

  power.numerator = -2;
  power.denominatorMinusOne = 0;
  XLAL_CHECK( XLALUnitRaiseRAT4( &unit2, &lalMeterUnit, &power ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitMultiply( &unit1, &lalWeberUnit, &unit2 ) != NULL, XLAL_EFAILED );
  XLAL_CHECK( XLALUnitCompare( &unit1, &lalTeslaUnit ) == 0, XLAL_EFAILED );
  printf("  PASS: T is Wb m^-2\n");

  LALCheckMemoryLeaks();

  printf("PASS: All tests\n");

  return 0;

}

/** \endcond */
