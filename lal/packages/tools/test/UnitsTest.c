/******************************** <lalVerbatim file="UnitsTestCV">
Author: Whelan, J. T.
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{UnitsTest.c}}
\label{ss:UnitsTest.c}

Test Suite for unit manipulation programs

\subsubsection*{Usage}
\begin{verbatim}

\subsubsection*{Usage}
\begin{verbatim}
UnitsTest [options]
Options:
  -h         print help
  -q         quiet: run silently
  -v         verbose: print extra information
  -d level   set lalDebugLevel to level
\end{verbatim}

\subsubsection*{Description}

This program tests the various units-manipulation routines, as well as
the pre-defined units defined in \texttt{Units.h}.  For each
successful test, it prints ``PASS'' to standard output.

\subsubsection*{Exit codes}
\input{UnitsTestCE}

\subsubsection*{Uses}
\begin{verbatim}
LALCHARCreateVector()
LALUnitAsString()
LALCHARDestroyVector()
LALUnitMultiply()
LALUnitRaise()
LALUnitNormalize()
LALUnitCompare()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{UnitsTestCV}}
******************************************************* </lalLaTeX> */

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

int lalDebugLevel = 1;
int verbose    = 0;

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

static void
TestStatus (LALStatus *status, const char *expectedCodes, int exitCode);

static void
ClearStatus (LALStatus *status);

NRCSID( UNITSTESTC, "$Id$" );

/***************************** <lalErrTable file="UnitsTestCE"> */
#define UNITSTESTC_ENOM 0
#define UNITSTESTC_ECHK 1
#define UNITSTESTC_EFLS 2
#define UNITSTESTC_MSGENOM "Nominal exit"
#define UNITSTESTC_MSGECHK "Error checking failed to catch bad data"
#define UNITSTESTC_MSGEFLS "Incorrect answer for valid data"
/***************************** </lalErrTable> */

/* The main function */
int main( int argc, char *argv[] )
{
  static LALStatus   status;
  static LALUnit     unit;
  static LALUnit     unity;
  LALUnitPair        unitPair;
  RAT4               power;
  CHARVector         dummy;
  CHARVector         *string;
  BOOLEAN            answer;

  ParseOptions( argc, argv );

  printf("Checking Input Validation of LALUnitAsString:\n");

  string = NULL;

  LALCHARCreateVector(&status, &string, sizeof("m^2 kg s^-3 A^-1")-1);
#ifndef NDEBUG

  LALUnitAsString( &status, NULL, &unit );
  TestStatus(&status, CODES(UNITSH_ENULLPOUT), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGENULLPOUT);

  dummy.data = NULL;

  LALUnitAsString( &status, &dummy, &unit );
  TestStatus(&status, CODES(UNITSH_ENULLPDOUT), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGENULLPDOUT);

  LALUnitAsString( &status, string, NULL );
  TestStatus(&status, CODES(UNITSH_ENULLPIN), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGENULLPIN);
#endif /* NDEBUG */

  unit = lalVoltUnit;

  LALUnitAsString( &status, string, &unit );
  TestStatus(&status, CODES(UNITSH_ESTRINGSIZE), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGESTRINGSIZE);

  LALCHARDestroyVector(&status, &string);

  printf("Testing response of LALUnitAsString to valid data:\n");

  LALCHARCreateVector(&status, &string, sizeof("m^2 kg s^-3 A^-1"));

  LALUnitAsString(&status, string, &unit);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (strcmp(string->data, "m^2 kg s^-3 A^-1"))
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: string representation of Volt\n");

  LALCHARDestroyVector(&status, &string);

  LALCHARCreateVector(&status, &string, sizeof("10^-12 m^-2 kg^-1 s^2 A^2"));
 
  unit = lalPicoFaradUnit;
  LALUnitAsString(&status, string, &unit);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (strcmp(string->data, "10^-12 m^-2 kg^-1 s^2 A^2"))
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: string representation of PicoFarad is '%s'\n",
	 string->data);
 
  unit = lalMeterUnit;
  LALUnitAsString(&status, string, &unit);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (strcmp(string->data, "m"))
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: string representation of Meter is '%s'\n", string->data);
 
  unit = lalCoulombUnit;
  LALUnitAsString(&status, string, &unit);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (strcmp(string->data, "s A"))
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: string representation of Coulomb is '%s'\n", string->data);
 
  unit = lalMegaUnit;
  LALUnitAsString(&status, string, &unit);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (strcmp(string->data, "10^6"))
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: string representation of Mega is '%s'\n", string->data);
 
  unit = lalAttoStrainUnit;
  LALUnitAsString(&status, string, &unit);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (strcmp(string->data, "10^-18 strain"))
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: string representation of AttoStrain is '%s'\n",
	 string->data);
 
  unit = unity;
  LALUnitAsString(&status, string, &unit);

  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (strcmp(string->data, ""))
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: string representation of unity is '%s'\n", string->data);
  LALCHARDestroyVector(&status, &string);

#ifndef NDEBUG

  printf("Checking Input Validation of LALUnitMultiply:\n");

  LALUnitMultiply( &status, NULL, &unitPair );
  TestStatus(&status, CODES(UNITSH_ENULLPOUT), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGENULLPOUT);

  LALUnitMultiply( &status, &unit, NULL );
  TestStatus(&status, CODES(UNITSH_ENULLPIN), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGENULLPIN);

  unitPair.unitOne.powerOfTen = 20000;
  unitPair.unitTwo.powerOfTen = 30000;

  LALUnitMultiply( &status, &unit, &unitPair );
  TestStatus(&status, CODES(UNITSH_EOVERFLOW), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGEOVERFLOW);

  unitPair.unitOne.powerOfTen = 0;
  unitPair.unitTwo.powerOfTen = 0;
  unitPair.unitOne.unitNumerator[2] = 12345;
  unitPair.unitOne.unitDenominatorMinusOne[2] = 23456;
  unitPair.unitTwo.unitNumerator[2] = 23456;
  unitPair.unitTwo.unitDenominatorMinusOne[2] = 12345;

  LALUnitMultiply( &status, &unit, &unitPair );
  TestStatus(&status, CODES(UNITSH_EOVERFLOW), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGEOVERFLOW);

  printf("Checking Input Validation of LALUnitRaise:\n");

  LALUnitRaise( &status, NULL, &unit, &power );
  TestStatus(&status, CODES(UNITSH_ENULLPOUT), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGENULLPOUT);

  LALUnitRaise( &status, &unit, NULL, &power );
  TestStatus(&status, CODES(UNITSH_ENULLPIN), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGENULLPIN);

  LALUnitRaise( &status, &(unitPair.unitOne), &unit, NULL );
  TestStatus(&status, CODES(UNITSH_ENULLPPARAM), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGENULLPPARAM);

  unit = lalKiloUnit;
  power.numerator = 1;
  power.denominatorMinusOne = 1;
  LALUnitRaise( &status, &(unitPair.unitOne), &unit, &power );
  TestStatus(&status, CODES(UNITSH_ENONINT), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGENONINT);

  unit.unitNumerator[2] = 40000;
  power.numerator = 2;
  power.denominatorMinusOne = 0;
  LALUnitRaise( &status, &(unitPair.unitOne), &unit, &power );
  TestStatus(&status, CODES(UNITSH_EOVERFLOW), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGEOVERFLOW);

  printf("Checking Input Validation of LALNormalize:\n");

  LALUnitRaise( &status, NULL, &unit, &power );
  TestStatus(&status, CODES(UNITSH_ENULLPOUT), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGENULLPOUT);

  LALUnitRaise( &status, &unit, NULL, &power );
  TestStatus(&status, CODES(UNITSH_ENULLPIN), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGENULLPIN);
  
#endif /* NDEBUG */

  printf("Testing response of LALUnitNormalize to valid data:\n");

  unit = unity;
  unit.unitNumerator[0] = 2;
  unit.unitDenominatorMinusOne[0] = 5;
  unitPair.unitTwo = unity;
  unitPair.unitTwo.unitNumerator[0] = 1;
  unitPair.unitTwo.unitDenominatorMinusOne[0] = 2;

  LALUnitNormalize( &status, &(unitPair.unitOne), &unit );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if (!answer)
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: 2/6 reduces to 1/3\n");  


#ifndef NDEBUG

  printf("Checking Input Validation of LALUnitCompare:\n");

  LALUnitCompare( &status, NULL, &unitPair );
  TestStatus(&status, CODES(UNITSH_ENULLPOUT), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGENULLPOUT);

  LALUnitCompare( &status, &answer, NULL );
  TestStatus(&status, CODES(UNITSH_ENULLPIN), UNITSTESTC_ECHK);
  printf("  PASS: %s\n", UNITSH_MSGENULLPIN);

#endif /* NDEBUG */

  printf("Testing response of LALUnitCompare to valid data:\n");

  unitPair.unitOne = lalMeterUnit;
  unitPair.unitTwo = lalMeterUnit;
  
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: m = m\n");

  unitPair.unitTwo = lalKiloGramUnit;
  
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( answer ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: m != kg\n");

  unitPair.unitOne = lalGramUnit;

  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( answer ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: g != kg\n");

  printf("Testing definitions of basic units:\n");

  unitPair.unitOne = unity;
  unitPair.unitOne.unitNumerator[LALUnitIndexMeter] = 1;
  unitPair.unitTwo = lalMeterUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer || unitPair.unitTwo.unitNumerator[LALUnitIndexMeter] != 1 ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: meter\n");

  unitPair.unitOne = unity;
  unitPair.unitOne.unitNumerator[LALUnitIndexKiloGram] = 1;
  unitPair.unitTwo = lalKiloGramUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer || unitPair.unitTwo.unitNumerator[LALUnitIndexKiloGram] != 1 ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: kilogram\n");

  unitPair.unitOne = unity;
  unitPair.unitOne.unitNumerator[LALUnitIndexSecond] = 1;
  unitPair.unitTwo = lalSecondUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer || unitPair.unitTwo.unitNumerator[LALUnitIndexSecond] != 1 ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: second\n");

  unitPair.unitOne = unity;
  unitPair.unitOne.unitNumerator[LALUnitIndexAmpere] = 1;
  unitPair.unitTwo = lalAmpereUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer || unitPair.unitTwo.unitNumerator[LALUnitIndexAmpere] != 1 ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: Ampere\n");

  unitPair.unitOne = unity;
  unitPair.unitOne.unitNumerator[LALUnitIndexKelvin] = 1;
  unitPair.unitTwo = lalKelvinUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer || unitPair.unitTwo.unitNumerator[LALUnitIndexKelvin] != 1 ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: Kelvin\n");

  unitPair.unitOne = unity;
  unitPair.unitOne.unitNumerator[LALUnitIndexStrain] = 1;
  unitPair.unitTwo = lalStrainUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer || unitPair.unitTwo.unitNumerator[LALUnitIndexStrain] != 1 ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: strain\n");

  unitPair.unitOne = unity;
  unitPair.unitOne.unitNumerator[LALUnitIndexADCCount] = 1;
  unitPair.unitTwo = lalADCCountUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer || unitPair.unitTwo.unitNumerator[LALUnitIndexADCCount] != 1 ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: ADC Count\n");

  printf("Testing definitions of derived units:\n");
  
  power.numerator = -1;
  power.denominatorMinusOne = 0;
  unit = lalSecondUnit;
  LALUnitRaise( &status, &(unitPair.unitOne), &unit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitTwo = lalHertzUnit;
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
  unit = lalSecondUnit;
  LALUnitRaise( &status, &unitPair.unitTwo, &unit, &power);
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = lalMeterUnit;
  LALUnitMultiply( &status, &unit, &unitPair); /* unit is now m/s^2 */
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = lalKiloGramUnit;
  unitPair.unitTwo = unit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = unit;
  unitPair.unitTwo = lalNewtonUnit;
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
  LALUnitRaise( &status, &(unitPair.unitTwo), &lalMeterUnit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = lalNewtonUnit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = unit;
  unitPair.unitTwo = lalPascalUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: Pa is N m^-2\n");

  unitPair.unitOne = lalNewtonUnit;
  unitPair.unitTwo = lalMeterUnit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = unit;
  unitPair.unitTwo = lalJouleUnit;
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
  LALUnitRaise( &status, &(unitPair.unitTwo), &lalSecondUnit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = lalJouleUnit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = unit;
  unitPair.unitTwo = lalWattUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: W is J/s\n");

  unitPair.unitOne = lalAmpereUnit;
  unitPair.unitTwo = lalSecondUnit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = unit;
  unitPair.unitTwo = lalCoulombUnit;

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
  LALUnitRaise( &status, &(unitPair.unitTwo), &lalAmpereUnit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = lalWattUnit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = unit;
  unitPair.unitTwo = lalVoltUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: V is W/A\n");

  /*
  unitPair.unitOne = lalCoulombUnit;
  unitPair.unitTwo = lalVoltUnit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = unit;
  unitPair.unitTwo = lalJouleUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: J is C V\n");
  */

  power.numerator = -1;
  power.denominatorMinusOne = 0;
  LALUnitRaise( &status, &(unitPair.unitTwo), &lalAmpereUnit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = lalVoltUnit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = unit;
  unitPair.unitTwo = lalOhmUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: Ohm is V/A\n");

  /*
  unitPair.unitOne = lalAmpereUnit;
  unitPair.unitTwo = lalOhmUnit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = unit;
  unitPair.unitTwo = lalVoltUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: V is A Ohm\n");
  */

  /*
  unitPair.unitOne = lalFaradUnit;
  unitPair.unitTwo = lalVoltUnit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = unit;
  unitPair.unitTwo = lalCoulombUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: C is F V\n");
  */

  power.numerator = -1;
  power.denominatorMinusOne = 0;
  LALUnitRaise( &status, &(unitPair.unitTwo), &lalVoltUnit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = lalCoulombUnit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = unit;
  unitPair.unitTwo = lalFaradUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: F is C/V\n");

  /*
  power.numerator = 2;
  power.denominatorMinusOne = 0;
  LALUnitRaise( &status, &(unitPair.unitOne), &lalOhmUnit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = lalAmpereUnit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = unit;
  unitPair.unitTwo = lalWattUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: W is A^2 Ohm\n");
  */

  unitPair.unitOne = lalVoltUnit;
  unitPair.unitTwo = lalSecondUnit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = unit;
  unitPair.unitTwo = lalWeberUnit;
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
  LALUnitRaise( &status, &(unitPair.unitTwo), &lalAmpereUnit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = lalSecondUnit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = lalVoltUnit;
  unitPair.unitTwo = unit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = unit;
  unitPair.unitTwo = lalHenryUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: H is V s/A\n");
  
  /*
  power.numerator = -1;
  power.denominatorMinusOne = 0;
  LALUnitRaise( &status, &(unitPair.unitTwo), &lalSecondUnit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = lalAmpereUnit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = lalHenryUnit;
  unitPair.unitTwo = unit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = unit;
  unitPair.unitTwo = lalVoltUnit;
  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: V is H A/s\n");
  */

  power.numerator = -2;
  power.denominatorMinusOne = 0;
  LALUnitRaise( &status, &(unitPair.unitTwo), &lalMeterUnit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);

  unitPair.unitOne = lalWeberUnit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = unit;
  unitPair.unitTwo = lalTeslaUnit;

  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: T is Wb m^-2\n");

  /*
  unitPair.unitOne = lalAmpereUnit;
  unitPair.unitTwo = lalMeterUnit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  power.numerator = -1;
  power.denominatorMinusOne = 0;
  LALUnitRaise( &status, &(unitPair.unitTwo), &unit, &power );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = lalNewtonUnit;
  LALUnitMultiply( &status, &unit, &unitPair);  
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  unitPair.unitOne = unit;
  unitPair.unitTwo = lalTeslaUnit;

  LALUnitCompare( &status, &answer, &unitPair );
  TestStatus(&status, CODES(0), UNITSTESTC_EFLS);
  if ( !answer ) 
  {
    fprintf (stderr, "\nExiting to system with code %d\n", UNITSTESTC_EFLS);
    return UNITSTESTC_EFLS;
  }
  printf("  PASS: T is N/(A m)\n");
  */

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
 *
 * ClearStatus ()
 *
 * Recursively applies DETATCHSTATUSPTR() to status structure to destroy
 * linked list of statuses.
 *
 */
void
ClearStatus (LALStatus *status)
{
  if (status->statusPtr)
  {
    ClearStatus      (status->statusPtr);
    DETATCHSTATUSPTR (status);
  }
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
        freopen ("/dev/null", "w", stderr);
        freopen ("/dev/null", "w", stdout);
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
