/**** <lalVerbatim file="TimeFreqFFTTestCV">
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * 
 * \subsection{Program \texttt{TimeFreqFFTTest.c}}
 * \label{ss:TimeFreqFFTTest.c}
 * 
 * Tests the routines in \verb+TimeFreqFFT.h+.
 * 
 * \subsection*{Usage}
 * \begin{verbatim}
 * TimeFreqFFTTest [options]
 * Options:
 *   -h         print this message
 *   -q         quiet: run silently
 *   -v         verbose: print extra information
 *   -d level   set lalDebugLevel to level
 * \end{verbatim}
 * 
 * \subsubsection*{Description}
 * \subsubsection*{Exit codes}
 * \begin{tabular}{|c|l|}
 * \hline
 * Code & Explanation                   \\
 * \hline
 * \tt 0 & Success, normal exit.         \\
 * \tt 1 & Subroutine failed.            \\
 * \hline
 * \end{tabular}
 * 
 * \subsubsection*{Uses}
 * \subsubsection*{Notes}
 * 
 * \vfill{\footnotesize\input{TimeFreqFFTTestCV}}
 * 
 **** </lalLaTeX> */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALStdio.h>
#include <lal/Units.h>
#include <lal/Random.h>
#include <lal/PrintFTSeries.h>
#include <lal/TimeFreqFFT.h>

#define CODES_(x) #x
#define CODES(x) CODES_(x)

extern char *optarg;
extern int   optind;

int lalDebugLevel = 0;
int verbose       = 0;

static void
Usage( const char *program, int exitflag );

static void
ParseOptions( int argc, char *argv[] );

static void
TestStatus( LALStatus *status, const char *expectedCodes, int exitCode );

static void
ClearStatus( LALStatus *status );


int main( int argc, char *argv[] )
{
  const UINT4 n  = 65536;
  const REAL4 dt = 1.0 / 16384.0;
  static LALStatus status;

  static REAL4TimeSeries         x;
  static COMPLEX8FrequencySeries X;

  static COMPLEX8TimeSeries      z;
  static COMPLEX8FrequencySeries Z;

  RealFFTPlan    *fwdRealPlan    = NULL;
  RealFFTPlan    *revRealPlan    = NULL;
  ComplexFFTPlan *fwdComplexPlan = NULL;
  ComplexFFTPlan *revComplexPlan = NULL;
  RandomParams   *rand           = NULL;

  UINT4 j;

  ParseOptions( argc, argv );

  LALSCreateVector( &status, &x.data, n );
  TestStatus( &status, CODES( 0 ), 1 );
  LALCCreateVector( &status, &X.data, n / 2 + 1 );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCCreateVector( &status, &z.data, n );
  TestStatus( &status, CODES( 0 ), 1 );
  LALCCreateVector( &status, &Z.data, n );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCreateForwardRealFFTPlan( &status, &fwdRealPlan, n, 0 );
  TestStatus( &status, CODES( 0 ), 1 );
  LALCreateReverseRealFFTPlan( &status, &revRealPlan, n, 0 );
  TestStatus( &status, CODES( 0 ), 1 );
  LALCreateForwardComplexFFTPlan( &status, &fwdComplexPlan, n, 0 );
  TestStatus( &status, CODES( 0 ), 1 );
  LALCreateReverseComplexFFTPlan( &status, &revComplexPlan, n, 0 );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCreateRandomParams( &status, &rand, 0 );
  TestStatus( &status, CODES( 0 ), 1 );


  /*
   *
   * Try the real transform.
   *
   */
  

  x.f0 = 0;
  x.deltaT = dt;
  x.sampleUnits = lalMeterUnit;
  LALSnprintf( x.name, sizeof( x.name ), "x" );
  LALNormalDeviates( &status, x.data, rand );
  TestStatus( &status, CODES( 0 ), 1 );
  for ( j = 0; j < n; ++j ) /* add a 60 Hz line */
  {
    REAL4 t = j * dt;
    x.data->data[j] += 0.1 * cos( LAL_TWOPI * 60.0 * t );
  }
  LALSPrintTimeSeries( &x, "x.out" );

  LALSnprintf( X.name, sizeof( X.name ), "X" );
  LALTimeFreqRealFFT( &status, &X, &x, fwdRealPlan );
  TestStatus( &status, CODES( 0 ), 1 );
  LALCPrintFrequencySeries( &X, "X.out" );

  LALFreqTimeRealFFT( &status, &x, &X, revRealPlan );
  TestStatus( &status, CODES( 0 ), 1 );
  LALSPrintTimeSeries( &x, "xx.out" );


  /*
   *
   * Try the complex transform.
   *
   */


  z.f0 = 0;
  z.deltaT = dt;
  z.sampleUnits = lalVoltUnit;
  LALSnprintf( z.name, sizeof( z.name ), "z" );
  { /* dirty hack */
    REAL4Vector tmp;
    tmp.length = 2 * z.data->length;
    tmp.data   = (REAL4 *)z.data->data;
    LALNormalDeviates( &status, &tmp, rand );
  }
  for ( j = 0; j < n; ++j ) /* add a 50 Hz line and a 500 Hz ringdown */
  {
    REAL4 t = j * dt;
    z.data->data[j].re += 0.2 * cos( LAL_TWOPI * 50.0 * t );
    z.data->data[j].im += exp( -t ) * sin( LAL_TWOPI * 500.0 * t );
  }
  LALCPrintTimeSeries( &z, "z.out" );
  TestStatus( &status, CODES( 0 ), 1 );

  LALSnprintf( Z.name, sizeof( Z.name ), "Z" );
  LALTimeFreqComplexFFT( &status, &Z, &z, fwdComplexPlan );
  TestStatus( &status, CODES( 0 ), 1 );
  LALCPrintFrequencySeries( &Z, "Z.out" );

  LALFreqTimeComplexFFT( &status, &z, &Z, revComplexPlan );
  TestStatus( &status, CODES( 0 ), 1 );
  LALCPrintTimeSeries( &z, "zz.out" );

  LALDestroyRandomParams( &status, &rand );
  TestStatus( &status, CODES( 0 ), 1 );

  LALDestroyRealFFTPlan( &status, &fwdRealPlan );
  TestStatus( &status, CODES( 0 ), 1 );
  LALDestroyRealFFTPlan( &status, &revRealPlan );
  TestStatus( &status, CODES( 0 ), 1 );
  LALDestroyComplexFFTPlan( &status, &fwdComplexPlan );
  TestStatus( &status, CODES( 0 ), 1 );
  LALDestroyComplexFFTPlan( &status, &revComplexPlan );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCDestroyVector( &status, &Z.data );
  TestStatus( &status, CODES( 0 ), 1 );
  LALCDestroyVector( &status, &z.data );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCDestroyVector( &status, &X.data );
  TestStatus( &status, CODES( 0 ), 1 );
  LALSDestroyVector( &status, &x.data );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCheckMemoryLeaks();
  return 0;
}


/*
 *
 * Override ordinary lalsupport version of these routines so that they put
 * a hash at the beginning of the comment lines.
 *
 */

void LALCPrintFrequencySeries(
    COMPLEX8FrequencySeries *series,
    const CHAR *filename
    ) 
{
  REAL8 f;
  COMPLEX8 *data;
  FILE *fp;
  UINT4 i;
  static LALStatus status;
  CHARVector *unitString;

  if (series==NULL) return;

  /* *(series->data) is a COMPLEX8Sequence */
  /* series->data->data is a pointer to COMPLEX8 */

  /* open output file */
  fp=LALFopen(filename,"w");
  fprintf(fp,"# %s\n",series->name);
  if (series->epoch.gpsSeconds && series->epoch.gpsNanoSeconds) {
    fprintf(fp,"# Epoch is %d seconds, %d nanoseconds\n",
            series->epoch.gpsSeconds,series->epoch.gpsNanoSeconds);
  }
  else {
    fprintf(fp,"#\n");
  }
  unitString = NULL;
  LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
  LALUnitAsString(&status, unitString, &(series->sampleUnits));
  fprintf(fp,"# Units are (%s)\n",unitString->data);
  fprintf(fp,"# Freq (Hz)\tRe(Value)\tIm(Value)\n");
  LALCHARDestroyVector(&status, &unitString);
  for ( i = 0; i < series->data->length; ++i )
  {
    f = series->f0 + i * series->deltaF;
    data = &(series->data->data[i]);
    fprintf(fp,"%e\t%e\t%e\n",f,data->re,data->im);
  }

  LALFclose(fp);

  return;
}


void LALCPrintTimeSeries( COMPLEX8TimeSeries *series, const CHAR *filename ) 
{
  REAL8 t;
  COMPLEX8 *data;
  FILE *fp;
  UINT4 i;
  static LALStatus status;
  CHARVector *unitString;

  if (series==NULL) return;

  /* *(series->data) is a COMPLEX8Sequence */
  /* series->data->data is a pointer to COMPLEX8 */

  /* Make a COMPLEX8 pointer which points to the first memory address not
   * belonging to the sequence
   */

  /* open output file */
  fp=LALFopen(filename,"w");
  fprintf(fp,"# %s\n",series->name);
  if (series->f0) {
     fprintf(fp,"# Heterodyned at %g Hz\n",series->f0);
  }
  else {
    fprintf(fp,"#\n");
  }
  fprintf(fp,"# Epoch is %d seconds, %d nanoseconds\n",
          series->epoch.gpsSeconds,series->epoch.gpsNanoSeconds);
  unitString = NULL;
  LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
  LALUnitAsString(&status, unitString, &(series->sampleUnits));
  fprintf(fp,"# Units are (%s)\n",unitString->data);
  fprintf(fp,"# Seconds since epoch\tRe(Value)\tIm(Value)\n");
  LALCHARDestroyVector(&status, &unitString);
  for ( i = 0; i < series->data->length; ++i )
  {
    t = i * series->deltaT;
    data = &(series->data->data[i]);
    fprintf(fp,"%e\t%e\t%e\n",t,data->re,data->im);
  }     

  LALFclose(fp);

  return;
}

void LALSPrintTimeSeries( REAL4TimeSeries *series, const CHAR *filename ) 
{
  REAL8 t;
  REAL4 *data;
  FILE *fp;
  UINT4 i;
  static LALStatus status;
  CHARVector *unitString;

  if (series==NULL) return;

  /* *(series->data) is a REAL4Sequence */
  /* series->data->data is a pointer to REAL4 */

  /* Make a REAL4 pointer which points to the first memory address not
   * belonging to the sequence
   */

  /* open output file */
  fp=LALFopen(filename,"w");
  fprintf(fp,"# %s\n",series->name);
  if (series->f0) {
     fprintf(fp,"# Heterodyned at %g Hz\n",series->f0);
  }
  else {
    fprintf(fp,"#\n");
  }
  fprintf(fp,"# Epoch is %d seconds, %d nanoseconds\n",
          series->epoch.gpsSeconds,series->epoch.gpsNanoSeconds);
  unitString = NULL;
  LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
  LALUnitAsString(&status, unitString, &(series->sampleUnits));
  fprintf(fp,"# Units are (%s)\n",unitString->data);
  fprintf(fp,"# Seconds since epoch\tValue\n");
  LALCHARDestroyVector(&status, &unitString);
  for ( i = 0; i < series->data->length; ++i )
  {
    t = i * series->deltaT;
    data = &(series->data->data[i]);
    fprintf(fp,"%e\t%e\n",t,*data);
  }     

  LALFclose(fp);

  return;
}


/*
 * TestStatus()
 *
 * Routine to check that the status code status->statusCode agrees with one of
 * the codes specified in the space-delimited string ignored; if not,
 * exit to the system with code exitcode.
 *
 */
static void
TestStatus( LALStatus *status, const char *ignored, int exitcode )
{
  char  str[64];
  char *tok;

  if ( verbose )
  {
    REPORTSTATUS( status );
  }

  if ( strncpy( str, ignored, sizeof( str ) ) )
  {
    if ( (tok = strtok( str, " " ) ) )
    {
      do
      {
        if ( status->statusCode == atoi( tok ) )
        {
          return;
        }
      }
      while ( ( tok = strtok( NULL, " " ) ) );
    }
    else
    {
      if ( status->statusCode == atoi( tok ) )
      {
        return;
      }
    }
  }

  fprintf( stderr, "\nExiting to system with code %d\n", exitcode );
  exit( exitcode );
}


/*
 *
 * ClearStatus()
 *
 * Recursively applies DETATCHSTATUSPTR() to status structure to destroy
 * linked list of statuses.
 *
 */
void
ClearStatus( LALStatus *status )
{
  if ( status->statusPtr )
  {
    ClearStatus( status->statusPtr );
    DETATCHSTATUSPTR( status );
  }
}


/*
 *
 * Usage()
 *
 * Prints a usage message for program program and exits with code exitcode.
 *
 */
static void
Usage( const char *program, int exitcode )
{
  fprintf( stderr, "Usage: %s [options]\n", program );
  fprintf( stderr, "Options:\n" );
  fprintf( stderr, "  -h         print this message\n" );
  fprintf( stderr, "  -q         quiet: run silently\n" );
  fprintf( stderr, "  -v         verbose: print extra information\n" );
  fprintf( stderr, "  -d level   set lalDebugLevel to level\n" );
  exit( exitcode );
}


/*
 * ParseOptions()
 *
 * Parses the argc - 1 option strings in argv[].
 *
 */
static void
ParseOptions( int argc, char *argv[] )
{
  while ( 1 )
  {
    int c = -1;

    c = getopt( argc, argv, "hqvd:" );
    if ( c == -1 )
    {
      break;
    }

    switch ( c )
    {
      case 'd': /* set debug level */
        lalDebugLevel = atoi( optarg );
        break;

      case 'v': /* verbose */
        ++verbose;
        break;

      case 'q': /* quiet: run silently (ignore error messages) */
        freopen( "/dev/null", "w", stderr );
        freopen( "/dev/null", "w", stdout );
        break;

      case 'h':
        Usage( argv[0], 0 );
        break;

      default:
        Usage( argv[0], 1 );
    }

  }

  if ( optind < argc )
  {
    Usage( argv[0], 1 );
  }

  return;
}
