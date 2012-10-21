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

/**
\author J. T. Whelan <john.whelan@ligo.org>
\file
\ingroup Units_h

\brief Test Suite for unit manipulation programs

\heading{Usage}
\code
UnitsTest [options]
Options:
  -h         print help
  -q         quiet: run silently
  -v         verbose: print extra information
  -d level   set lalDebugLevel to level
\endcode

\heading{Description}

This program tests the various units-manipulation routines, as well as
the pre-defined units defined in \ref Units.h.  For each
successful test, it prints "PASS" to standard output.

\heading{Uses}
\code
LALCHARCreateVector()
LALUnitAsString()
LALParseUnitString()
LALCHARDestroyVector()
LALUnitMultiply()
LALUnitRaise()
LALUnitNormalize()
LALUnitCompare()
\endcode

\heading{Notes}

*/
/** @{ */
/**\name Error Codes */ /*@{*/
#define UNITSTESTC_ENOM 0
#define UNITSTESTC_ECHK 1
#define UNITSTESTC_EFLS 2
#define UNITSTESTC_MSGENOM "Nominal exit"
#define UNITSTESTC_MSGECHK "Error checking failed to catch bad data"
#define UNITSTESTC_MSGEFLS "Incorrect answer for valid data"
/*@}*/
/** @} */

/** \cond  DONT_DOXYGEN */
#include <config.h>

#include <stdlib.h>

#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#define CODES_(x) #x
#define CODES(x) CODES_(x)

extern char *optarg;
extern int   optind;

extern int lalDebugLevel;
int verbose    = 0;

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

static void
TestStatus (LALStatus *status, const char *expectedCodes, int exitCode);

/* The main function */
int main( int argc, char *argv[] )
{
  static LALStatus   status;
  LALUnit            unit1, unit2;
  LALUnitPair        unitPair;
  RAT4               power;
  CHARVector         dummy;
  CHARVector         *string;
  BOOLEAN            answer;

  lalDebugLevel = LALMSGLVL1;

  ParseOptions( argc, argv );

  printf("Checking Input Validation of LALUnitAsString:\n");

  string = NULL;

  LALCHARCreateVector(&status, &string, sizeof("m^2 kg s^-3 A^-1")-1);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
#ifndef LAL_NDEBUG

  if ( ! lalNoDebug )
  {
    LALUnitAsString( &status, NULL, &unit1 );
    TestStatus(&status, CODES(UNITSH_ENULLPOUT), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGENULLPOUT);

    dummy.data = NULL;

    LALUnitAsString( &status, &dummy, &unit1 );
    TestStatus(&status, CODES(UNITSH_ENULLPD), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGENULLPD);

    LALUnitAsString( &status, string, NULL );
    TestStatus(&status, CODES(UNITSH_ENULLPIN), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGENULLPIN);
  }
#else
  (void)dummy;
#endif /* LAL_NDEBUG */

  LALUnitAsString( &status, string, &lalVoltUnit );
  TestStatus(&status, CODES(UNITSH_ESTRINGSIZE), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGESTRINGSIZE);

  LALCHARDestroyVector(&status, &string);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  printf("Testing response of LALUnitAsString to valid data:\n");

  LALCHARCreateVector(&status, &string, sizeof("m^2 kg s^-3 A^-1"));
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  LALUnitAsString(&status, string, &lalVoltUnit);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (strcmp(string->data, "m^2 kg s^-3 A^-1"))
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: string representation of Volt\n");

  LALCHARDestroyVector(&status, &string);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  LALCHARCreateVector(&status, &string, sizeof("10^-12 m^-2 kg^-1 s^2 A^2"));
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  LALUnitAsString(&status, string, &lalPicoFaradUnit);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (strcmp(string->data, "10^-12 m^-2 kg^-1 s^2 A^2"))
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: string representation of PicoFarad is '%s'\n",
	 string->data);

  LALUnitAsString(&status, string, &lalMeterUnit);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (strcmp(string->data, "m"))
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: string representation of Meter is '%s'\n", string->data);

  LALUnitAsString(&status, string, &lalCoulombUnit);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (strcmp(string->data, "s A"))
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: string representation of Coulomb is '%s'\n", string->data);

  LALUnitAsString(&status, string, &lalMegaUnit);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (strcmp(string->data, "10^6"))
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: string representation of Mega is '%s'\n", string->data);

  LALUnitAsString(&status, string, &lalAttoStrainUnit);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (strcmp(string->data, "10^-18 strain"))
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: string representation of AttoStrain is '%s'\n",
	 string->data);

  LALUnitAsString(&status, string, &lalDimensionlessUnit);

  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (strcmp(string->data, ""))
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: string representation of dimensionless is '%s'\n",
	 string->data);
  LALCHARDestroyVector(&status, &string);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  printf("Checking Input Validation of LALParseUnitString:\n");

  string = NULL;

  LALCHARCreateVector(&status, &string, LALUnitTextSize);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
#ifndef LAL_NDEBUG

  if ( ! lalNoDebug )
  {
    LALParseUnitString( &status,  &unit1, NULL );
    TestStatus(&status, CODES(UNITSH_ENULLPIN), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGENULLPIN);

    dummy.data = NULL;

    LALParseUnitString( &status, &unit1, &dummy );
    TestStatus(&status, CODES(UNITSH_ENULLPD), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGENULLPD);

    LALParseUnitString( &status, NULL, string );
    TestStatus(&status, CODES(UNITSH_ENULLPOUT), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGENULLPOUT);
  }
#endif /* LAL_NDEBUG */

  strncpy( string->data, "10^", string->length );

  LALParseUnitString( &status, &unit1, string );
  TestStatus(&status, CODES(UNITSH_EPARSE), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGEPARSE);

  strncpy( string->data, "s^3 s^2", string->length );

  LALParseUnitString( &status, &unit1, string );
  TestStatus(&status, CODES(UNITSH_EPARSE), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGEPARSE);

  strncpy( string->data, "V", string->length );

  LALParseUnitString( &status, &unit1, string );
  TestStatus(&status, CODES(UNITSH_EPARSE), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGEPARSE);

  printf("Testing response of LALParseUnitString to valid data:\n");

  unitPair.unitOne = &unit1;
  unitPair.unitTwo = &lalVoltUnit;

  LALUnitAsString(&status, string, &lalVoltUnit);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  LALParseUnitString(&status, &unit1, string);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (!answer)
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: Volt <-> string\n");

  unitPair.unitTwo = &lalPicoFaradUnit;

  LALUnitAsString(&status, string, unitPair.unitTwo);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  LALParseUnitString(&status, &unit1, string);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (!answer)
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: PicoFarad <-> string\n");

  unitPair.unitTwo = &lalMeterUnit;

  LALUnitAsString(&status, string, unitPair.unitTwo);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  LALParseUnitString(&status, &unit1, string);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (!answer)
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: Meter <-> string\n");

  unitPair.unitTwo = &lalAttoStrainUnit;

  LALUnitAsString(&status, string, unitPair.unitTwo);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  LALParseUnitString(&status, &unit1, string);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (!answer)
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: AttoStrain <-> string\n");

  unitPair.unitTwo = &lalCoulombUnit;

  LALUnitAsString(&status, string, unitPair.unitTwo);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  LALParseUnitString(&status, &unit1, string);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (!answer)
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: Coulomb <-> string\n");

  unitPair.unitTwo = &lalMegaUnit;

  LALUnitAsString(&status, string, unitPair.unitTwo);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  LALParseUnitString(&status, &unit1, string);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (!answer)
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: Mega <-> string\n");

  unitPair.unitTwo = &lalDimensionlessUnit;

  LALUnitAsString(&status, string, unitPair.unitTwo);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  LALParseUnitString(&status, &unit1, string);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (!answer)
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: Dimensionless <-> string\n");


  LALCHARDestroyVector(&status, &string);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);


#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    unit1 = unit2 = lalDimensionlessUnit;

    unitPair.unitTwo = &unit2;

    printf("Checking Input Validation of LALUnitMultiply:\n");

    LALUnitMultiply( &status, NULL, &unitPair );
    TestStatus(&status, CODES(UNITSH_ENULLPOUT), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGENULLPOUT);

    LALUnitMultiply( &status, &unit1, NULL );
    TestStatus(&status, CODES(UNITSH_ENULLPIN), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGENULLPIN);

    unit1.powerOfTen = 20000;
    unit2.powerOfTen = 30000;

    LALUnitMultiply( &status, &unit1, &unitPair );
    TestStatus(&status, CODES(UNITSH_EOVERFLOW), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGEOVERFLOW);

    unit1.powerOfTen = 0;
    unit2.powerOfTen = 0;
    unit1.unitNumerator[2] = 12345;
    unit1.unitDenominatorMinusOne[2] = 23456;
    unit2.unitNumerator[2] = 23456;
    unit2.unitDenominatorMinusOne[2] = 12345;

    LALUnitMultiply( &status, &unit1, &unitPair );
    TestStatus(&status, CODES(UNITSH_EOVERFLOW), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGEOVERFLOW);

    printf("Checking Input Validation of LALUnitRaise:\n");

    LALUnitRaise( &status, NULL, &unit1, &power );
    TestStatus(&status, CODES(UNITSH_ENULLPOUT), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGENULLPOUT);

    LALUnitRaise( &status, &unit1, NULL, &power );
    TestStatus(&status, CODES(UNITSH_ENULLPIN), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGENULLPIN);

    LALUnitRaise( &status, &unit1, &unit1, NULL );
    TestStatus(&status, CODES(UNITSH_ENULLPPARAM), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGENULLPPARAM);

    unit1 = lalKiloUnit;
    power.numerator = 1;
    power.denominatorMinusOne = 1;
    LALUnitRaise( &status, &unit1, &unit1, &power );
    TestStatus(&status, CODES(UNITSH_ENONINT), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGENONINT);

    unit1.unitNumerator[2] = 20000;
    unit1.unitNumerator[2] += 20000;
    power.numerator = 2;
    power.denominatorMinusOne = 0;
    LALUnitRaise( &status, &unit1, &unit1, &power );
    TestStatus(&status, CODES(UNITSH_EOVERFLOW), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGEOVERFLOW);

    printf("Checking Input Validation of LALNormalize:\n");

    LALUnitRaise( &status, NULL, &unit1, &power );
    TestStatus(&status, CODES(UNITSH_ENULLPOUT), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGENULLPOUT);

    LALUnitRaise( &status, &unit1, NULL, &power );
    TestStatus(&status, CODES(UNITSH_ENULLPIN), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGENULLPIN);
  }
#endif /* LAL_NDEBUG */

  printf("Testing response of LALUnitNormalize to valid data:\n");

  unit2 = lalDimensionlessUnit;
  unit2.unitNumerator[0] = 2;
  unit2.unitDenominatorMinusOne[0] = 5;

  LALUnitNormalize( &status, &unit1, &unit2 );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  unit2 = lalDimensionlessUnit;
  unit2.unitNumerator[0] = 1;
  unit2.unitDenominatorMinusOne[0] = 2;

  unitPair.unitOne = &unit1;
  unitPair.unitTwo = &unit2;

  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (!answer)
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: 2/6 reduces to 1/3\n");


#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    printf("Checking Input Validation of LALUnitCompare:\n");

    LALUnitCompare( &status, NULL, &unitPair );
    TestStatus(&status, CODES(UNITSH_ENULLPOUT), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGENULLPOUT);

    LALUnitCompare( &status, &answer, NULL );
    TestStatus(&status, CODES(UNITSH_ENULLPIN), UNITSTESTC_ECHK);
    printf("  PASS: %s\n", UNITSH_MSGENULLPIN);
  }
#endif /* LAL_NDEBUG */

  printf("Testing response of LALUnitCompare to valid data:\n");

  unitPair.unitOne = &lalMeterUnit;
  unitPair.unitTwo = &lalMeterUnit;

  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: m = m\n");

  unitPair.unitTwo = &lalKiloGramUnit;

  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( answer )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: m != kg\n");

  unitPair.unitOne = &lalGramUnit;

  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( answer )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: g != kg\n");

  printf("Testing definitions of basic units:\n");

  unitPair.unitOne = &unit1;

  unit1 = lalDimensionlessUnit;
  unit1.unitNumerator[LALUnitIndexMeter] = 1;
  unitPair.unitTwo = &lalMeterUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer || unitPair.unitTwo->unitNumerator[LALUnitIndexMeter] != 1 )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: meter\n");

  unit1 = lalDimensionlessUnit;
  unit1.unitNumerator[LALUnitIndexKiloGram] = 1;
  unitPair.unitTwo = &lalKiloGramUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer || unitPair.unitTwo->unitNumerator[LALUnitIndexKiloGram] != 1 )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: kilogram\n");

  unit1 = lalDimensionlessUnit;
  unit1.unitNumerator[LALUnitIndexSecond] = 1;
  unitPair.unitTwo = &lalSecondUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer || unitPair.unitTwo->unitNumerator[LALUnitIndexSecond] != 1 )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: second\n");

  unit1 = lalDimensionlessUnit;
  unit1.unitNumerator[LALUnitIndexAmpere] = 1;
  unitPair.unitTwo = &lalAmpereUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer || unitPair.unitTwo->unitNumerator[LALUnitIndexAmpere] != 1 )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: Ampere\n");

  unit1 = lalDimensionlessUnit;
  unit1.unitNumerator[LALUnitIndexKelvin] = 1;
  unitPair.unitTwo = &lalKelvinUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer || unitPair.unitTwo->unitNumerator[LALUnitIndexKelvin] != 1 )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: Kelvin\n");

  unit1 = lalDimensionlessUnit;
  unit1.unitNumerator[LALUnitIndexStrain] = 1;
  unitPair.unitTwo = &lalStrainUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer || unitPair.unitTwo->unitNumerator[LALUnitIndexStrain] != 1 )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: strain\n");

  unit1 = lalDimensionlessUnit;
  unit1.unitNumerator[LALUnitIndexADCCount] = 1;
  unitPair.unitTwo = &lalADCCountUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer || unitPair.unitTwo->unitNumerator[LALUnitIndexADCCount] != 1 )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: ADC Count\n");

  printf("Testing definitions of derived units:\n");

  power.numerator = -1;
  power.denominatorMinusOne = 0;
  LALUnitRaise( &status, &unit1, &lalSecondUnit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitTwo = &lalHertzUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: Hz is s^-1\n");

  power.numerator = -2;
  power.denominatorMinusOne = 0;
  LALUnitRaise( &status, &unit2, &lalSecondUnit, &power);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &lalMeterUnit;
  unitPair.unitTwo = &unit2;
  LALUnitMultiply( &status, &unit1, &unitPair); /* unit is now m/s^2 */
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &lalKiloGramUnit;
  unitPair.unitTwo = &unit1;
  LALUnitMultiply( &status, &unit2, &unitPair);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &unit2;
  unitPair.unitTwo = &lalNewtonUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: N is kg m s^-2\n");

  power.numerator = -2;
  power.denominatorMinusOne = 0;
  LALUnitRaise( &status, &unit2, &lalMeterUnit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &lalNewtonUnit;
  unitPair.unitTwo = &unit2;
  LALUnitMultiply( &status, &unit1, &unitPair);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &unit1;
  unitPair.unitTwo = &lalPascalUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: Pa is N m^-2\n");

  unitPair.unitOne = &lalNewtonUnit;
  unitPair.unitTwo = &lalMeterUnit;
  LALUnitMultiply( &status, &unit1, &unitPair);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &unit1;
  unitPair.unitTwo = &lalJouleUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: J is N m\n");

  power.numerator = -1;
  power.denominatorMinusOne = 0;
  LALUnitRaise( &status, &unit2, &lalSecondUnit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &lalJouleUnit;
  unitPair.unitTwo = &unit2;
  LALUnitMultiply( &status, &unit1, &unitPair);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &unit1;
  unitPair.unitTwo = &lalWattUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: W is J/s\n");

  unitPair.unitOne = &lalAmpereUnit;
  unitPair.unitTwo = &lalSecondUnit;
  LALUnitMultiply( &status, &unit1, &unitPair);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &unit1;
  unitPair.unitTwo = &lalCoulombUnit;

  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: C is A s\n");

  power.numerator = -1;
  power.denominatorMinusOne = 0;
  LALUnitRaise( &status, &unit2, &lalAmpereUnit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &lalWattUnit;
  unitPair.unitTwo = &unit2;
  LALUnitMultiply( &status, &unit1, &unitPair);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &unit1;
  unitPair.unitTwo = &lalVoltUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: V is W/A\n");

  power.numerator = -1;
  power.denominatorMinusOne = 0;
  LALUnitRaise( &status, &unit2, &lalAmpereUnit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &lalVoltUnit;
  unitPair.unitTwo = &unit2;
  LALUnitMultiply( &status, &unit1, &unitPair);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &unit1;
  unitPair.unitTwo = &lalOhmUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: Ohm is V/A\n");

  power.numerator = -1;
  power.denominatorMinusOne = 0;
  LALUnitRaise( &status, &unit2, &lalVoltUnit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &lalCoulombUnit;
  unitPair.unitTwo = &unit2;
  LALUnitMultiply( &status, &unit1, &unitPair);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &unit1;
  unitPair.unitTwo = &lalFaradUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: F is C/V\n");

  unitPair.unitOne = &lalVoltUnit;
  unitPair.unitTwo = &lalSecondUnit;
  LALUnitMultiply( &status, &unit1, &unitPair);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &unit1;
  unitPair.unitTwo = &lalWeberUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: Wb is V s\n");

  power.numerator = -1;
  power.denominatorMinusOne = 0;
  LALUnitRaise( &status, &unit2, &lalAmpereUnit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &lalSecondUnit;
  unitPair.unitTwo = &unit2;
  LALUnitMultiply( &status, &unit1, &unitPair);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &lalVoltUnit;
  unitPair.unitTwo = &unit1;
  LALUnitMultiply( &status, &unit1, &unitPair);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &unit1;
  unitPair.unitTwo = &lalHenryUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: H is V s/A\n");

  power.numerator = -2;
  power.denominatorMinusOne = 0;
  LALUnitRaise( &status, &unit2, &lalMeterUnit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  unitPair.unitOne = &lalWeberUnit;
  unitPair.unitTwo = &unit2;
  LALUnitMultiply( &status, &unit1, &unitPair);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = &unit1;
  unitPair.unitTwo = &lalTeslaUnit;

  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer )
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: T is Wb m^-2\n");

  LALCheckMemoryLeaks();

  printf("PASS: All tests\n");

  return UNITSTESTC_ENOM;
}


/*----------------------------------------------------------------------*/


/*
 * TestStatus ()
 *
 * Routine to check that the status code status->statusCode agrees with one of
 * the codes specified in the space-delimited string ignored; if not,
 * exit to the system with code exitcode.
 *
 */
static void
TestStatus (LALStatus *status, const char *ignored, int exitcode)
{
  char  str[64];
  char *tok;

  if (verbose)
  {
    REPORTSTATUS (status);
  }

  if (strncpy (str, ignored, sizeof (str)))
  {
    if ((tok = strtok (str, " ")))
    {
      do
      {
        if (status->statusCode == atoi (tok))
        {
          return;
        }
      }
      while ((tok = strtok (NULL, " ")));
    }
    else
    {
      if (status->statusCode == atoi (tok))
      {
        return;
      }
    }
  }

  fprintf (stderr, "\nExiting to system with code %d\n", exitcode);
  exit (exitcode);
}

/*
 * Usage ()
 *
 * Prints a usage message for program program and exits with code exitcode.
 *
 */
static void
Usage (const char *program, int exitcode)
{
  fprintf (stderr, "Usage: %s [options]\n", program);
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "  -h         print this message\n");
  fprintf (stderr, "  -q         quiet: run silently\n");
  fprintf (stderr, "  -v         verbose: print extra information\n");
  fprintf (stderr, "  -d level   set lalDebugLevel to level\n");
  exit (exitcode);
}


/*
 * ParseOptions ()
 *
 * Parses the argc - 1 option strings in argv[].
 *
 */
static void
ParseOptions (int argc, char *argv[])
{
  FILE *fp;

  while (1)
  {
    int c = -1;

    c = getopt (argc, argv, "hqvd:");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'd': /* set debug level */
        lalDebugLevel = atoi (optarg);
        break;

      case 'v': /* verbose */
        ++verbose;
        break;

      case 'q': /* quiet: run silently (ignore error messages) */
        fp = freopen ("/dev/null", "w", stderr);
        if (fp == NULL)
        {
          fprintf(stderr, "Error: unable to open /dev/null\n");
          exit(1);
        }
        fp = freopen ("/dev/null", "w", stdout);
        if (fp == NULL)
        {
          fprintf(stderr, "Error: unable to open /dev/null\n");
          exit(1);
        }
        break;

      case 'h':
        Usage (argv[0], 0);
        break;

      default:
        Usage (argv[0], 1);
    }

  }

  if (optind < argc)
  {
    Usage (argv[0], 1);
  }

  return;
}
/** \endcond */
