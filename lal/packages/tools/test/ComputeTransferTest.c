/**** <lalVerbatim file="ComputeTransferTestCV">
 * Author: Patrick R. Brady
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * \subsection{Program \texttt{ComputeTransferTest.c}}
 * 
 * Tests the computation of transfer functions.
 *
 * \subsubsection*{Usage}
 *
 * \begin{verbatim}
 * ComputeTransferTest
 * \end{verbatim}
 *
 * \subsubsection*{Description}
 *
 **** </lalLaTeX> */

#include <stdio.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/PrintFTSeries.h>
#include <lal/Calibration.h>

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS( pstat ); return 1; } \
  else ((void)0)

#define TESTSTATUSERR( pstat, code ) \
  if ( (pstat)->statusCode != code ) { REPORTSTATUS( pstat ); return 1; } \
  else ((void)0)

const char *usage = "Usage: %s [options]\nOptions:\n"
    "\t-h\t\tprint this message and exit\n"
    "\t-d lvl\t\tset debug level to lvl\n"
    "\t-v\t\tverbose output\n";


int lalDebugLevel = LALMEMDBG;
int verbose = 0;
const char *program;

int main( int argc, char *argv[] )
{
  enum { TransferLength = 16384 };
  const REAL8 df = 1.0;
  static LALStatus status;
  static CalibrationRecord calrec;
  static LIGOTimeGPS now;
  static LALUnit units;
  REAL8 poles[] = { 0.3, 0.5 };
  REAL8 zeros[] = { 0.1, 0.2, 0.4 };
  REAL8Vector pvec;
  REAL8Vector zvec;
  COMPLEX8 fdata[TransferLength];
  COMPLEX8Vector fvec;
  COMPLEX8FrequencySeries fser;

  /* initialization */
  pvec.length = sizeof( poles ) / sizeof( *poles );
  pvec.data   = poles;
  zvec.length = sizeof( zeros ) / sizeof( *zeros );
  zvec.data   = zeros;
  fvec.length = sizeof( fdata ) / sizeof( *fdata );
  fvec.data   = fdata;
  fser.epoch  = now;
  fser.f0     = 0;
  fser.deltaF = df;
  fser.data   = &fvec;
  fser.sampleUnits = units;
  strncpy( fser.name, "transfer", sizeof( fser.name ) );

  /* parse options */
  program = *argv;
  while ( ++argv && --argc > 0 )
  {
    if ( ! strcmp( *argv, "-d" ) )
    {
      lalDebugLevel = atoi( *++argv );
      --argc;
    }
    else if ( ! strcmp( *argv, "-v" ) )
    {
      verbose = 1;
    }
    else if ( ! strcmp( *argv, "-h" ) )
    {
      fprintf( stderr, usage, program );
      return 0;
    }
    else
    {
      fprintf( stderr, "invalid option %s\n", *argv );
      fprintf( stderr, usage, program );
    }
  }

  calrec.type     = CalibrationZPG;
  calrec.zeros    = &zvec;
  calrec.poles    = &pvec;
  calrec.transfer = &fser;

  /* check response to invalid input */
# ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    CalibrationRecord badcalrec;

    LALComputeTransfer( &status, NULL );
    TESTSTATUSERR( &status, CALIBRATIONH_ENULL );

    badcalrec = calrec;
    badcalrec.zeros = NULL;
    LALComputeTransfer( &status, &badcalrec );
    TESTSTATUSERR( &status, CALIBRATIONH_ENULL );

    badcalrec = calrec;
    badcalrec.poles = NULL;
    LALComputeTransfer( &status, &badcalrec );
    TESTSTATUSERR( &status, CALIBRATIONH_ENULL );

    badcalrec = calrec;
    badcalrec.transfer = NULL;
    LALComputeTransfer( &status, &badcalrec );
    TESTSTATUSERR( &status, CALIBRATIONH_ENULL );

    fputs( "PASS: Test response to invalid data\n", stderr );
  }
# endif /* LAL_NDEBUG */

  LALComputeTransfer( &status, &calrec );
  TESTSTATUS( &status );
  LALCPrintFrequencySeries( calrec.transfer, "transfer.out" );

  fputs( "PASS: Test response to valid data\n", stderr );

  LALCheckMemoryLeaks();
  fputs( "PASS: All tests\n", stderr );
  return 0;
}
