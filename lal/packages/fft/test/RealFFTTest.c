#if 0 /* autodoc block */

<lalVerbatim file="RealFFTTestCV">
$Id$
</lalVerbatim>

<lalLaTeX>

\subsection{Program \texttt{RealFFTTest.c}}
\label{ss:RealFFTTest.c}

Tests the routines in \verb+RealFFT.h+.

\subsection*{Usage}
\begin{verbatim}
RealFFTTest [options]
Options:
  -h         print this message
  -q         quiet: run silently
  -v         verbose: print extra information
  -d level   set lalDebugLevel to level
\end{verbatim}

\subsubsection*{Description}
\subsubsection*{Exit codes}
\begin{tabular}{|c|l|}
\hline
 Code & Explanation                   \\
\hline
\tt 0 & Success, normal exit.         \\
\tt 1 & Subroutine failed.            \\
\hline
\end{tabular}

\subsubsection*{Uses}
\subsubsection*{Notes}

\vfill{\footnotesize\input{RealFFTTestCV}}

</lalLaTeX>

#endif /* autodoc block */

#include <stdio.h>
#include <math.h>
#include <lal/LALConfig.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/VectorOps.h>

#define CODES_(x) #x
#define CODES(x) CODES_(x)

NRCSID (MAIN, "$Id$");

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
  const UINT4 m = 4;   /* example length of sequence of vectors */
  const UINT4 n = 32;  /* example vector length */

  static LALStatus status; 

  RealFFTPlan            *pfwd = NULL;
  RealFFTPlan            *pinv = NULL;
  REAL4Vector            *hvec = NULL;
  COMPLEX8Vector         *Hvec = NULL;
  REAL4Vector            *Pvec = NULL;
  REAL4VectorSequence    *hseq = NULL;
  COMPLEX8VectorSequence *Hseq = NULL;
  REAL4VectorSequence    *Pseq = NULL;
  CreateVectorSequenceIn  seqinp;
  CreateVectorSequenceIn  seqout;

  UINT4 i;
  UINT4 j;
  FILE *fp;

  ParseOptions( argc, argv );

  fp = verbose ? stdout : NULL ;  

  /* create vectors and sequences */

  seqinp.length       = m;
  seqinp.vectorLength = n;
  seqout.length       = m;
  seqout.vectorLength = n/2 + 1;

  LALSCreateVector( &status, &hvec, n );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCCreateVector( &status, &Hvec, n/2 + 1 );
  TestStatus( &status, CODES( 0 ), 1 );

  LALSCreateVector( &status, &Pvec, n/2 + 1 );
  TestStatus( &status, CODES( 0 ), 1 );

  LALSCreateVectorSequence( &status, &hseq, &seqinp );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCCreateVectorSequence( &status, &Hseq, &seqout );
  TestStatus( &status, CODES( 0 ), 1 );

  LALSCreateVectorSequence( &status, &Pseq, &seqout );
  TestStatus( &status, CODES( 0 ), 1 );

  /* create plans */

  LALEstimateFwdRealFFTPlan( &status, &pfwd, n );
  TestStatus( &status, CODES( 0 ), 1 );

  LALEstimateInvRealFFTPlan( &status, &pinv, n );
  TestStatus( &status, CODES( 0 ), 1 );

  LALDestroyRealFFTPlan( &status, &pfwd );
  TestStatus( &status, CODES( 0 ), 1 );

  LALDestroyRealFFTPlan( &status, &pinv );
  TestStatus( &status, CODES( 0 ), 1 );

  LALMeasureFwdRealFFTPlan( &status, &pfwd, n );
  TestStatus( &status, CODES( 0 ), 1 );

  LALMeasureInvRealFFTPlan( &status, &pinv, n );
  TestStatus( &status, CODES( 0 ), 1 );

  /* assign data ... */
  fp ? fprintf( fp, "\nInitial data:\n" ) : 0;
  for ( i = 0; i < n; ++i )
  {
    REAL4 data = cos( 7*i ) - 1;
    hvec->data[i] = data;
    for ( j = 0; j < m; ++j )
    {
      hseq->data[i + j*n] = data;
      fp ? fprintf( fp, "% 9.3f\t", data ) : 0;
      data += 1;
    }
    fp ? fprintf( fp, "\n" ) : 0;
  }

  /* perform FFTs */

  fp ? fprintf( fp, "\nSingle Forward FFT:\n" ) : 0;
  LALFwdRealFFT( &status, Hvec, hvec, pfwd );
  TestStatus( &status, CODES( 0 ), 1 );
  for ( i = 0; i < Hvec->length; ++i )
    fp ? fprintf( fp, "(% 9.3f, % 9.3f)\n", Hvec->data[i].re, Hvec->data[i].im ) : 0;

  fp ? fprintf( fp, "\nSingle Forward FFT Power:\n" ) : 0;
  LALRealPowerSpectrum( &status, Pvec, hvec, pfwd );
  TestStatus( &status, CODES( 0 ), 1 );
  for ( i = 0; i < Pvec->length; ++i )
    fp ? fprintf( fp, "%12.3f\n", Pvec->data[i] ) : 0;

  fp ? fprintf( fp, "\nSingle Inverse FFT:\n") : 0;
  LALInvRealFFT (&status, hvec, Hvec, pinv);
  TestStatus( &status, CODES( 0 ), 1 );
  for (i = 0; i < hvec->length; ++i)
    fp ? fprintf( fp, "% 9.3f\n", hvec->data[i]/n ) : 0;

  fp ? fprintf( fp, "\nMultiple Forward FFT:\n") : 0;
  LALFwdRealSequenceFFT (&status, Hseq, hseq, pfwd);
  TestStatus( &status, CODES( 0 ), 1 );
  for ( i = 0; i < Hseq->vectorLength; ++i )
  {
    for ( j = 0; j < Hseq->length; ++j )
      fp ? fprintf( fp, "(% 9.3f, % 9.3f)\t",
                    Hseq->data[i + j*Hseq->vectorLength].re,
                    Hseq->data[i + j*Hseq->vectorLength].im) : 0;
    fp ? fprintf( fp, "\n" ) : 0;
  }

  fp ? fprintf( fp, "\nMultiple Forward FFT Power:\n" ) : 0;
  LALRealSequencePowerSpectrum( &status, Pseq, hseq, pfwd );
  TestStatus( &status, CODES( 0 ), 1 );
  for ( i = 0; i < Pseq->vectorLength; ++i )
  {
    for ( j = 0; j < Pseq->length; ++j )
      fp ? fprintf( fp, "%12.3f\t", Pseq->data[i + j*Pseq->vectorLength] ) : 0;
    fp ? fprintf( fp, "\n" ) : 0;
  }

  fp ? fprintf( fp, "\nMultiple Inverse FFT:\n" ) : 0;
  LALInvRealSequenceFFT( &status, hseq, Hseq, pinv );
  TestStatus( &status, CODES( 0 ), 1 );
  for ( i = 0; i < hseq->vectorLength; ++i )
  {
    for ( j = 0; j < hseq->length; ++j )
      fp ? fprintf( fp, "% 9.3f\t", hseq->data[i + j*hseq->vectorLength]/n ) : 0;
    fp ? fprintf( fp, "\n" ) : 0;
  }

  /* destroy plans, vectors, and sequences */
  LALDestroyRealFFTPlan( &status, &pfwd );
  TestStatus( &status, CODES( 0 ), 1 );

  LALDestroyRealFFTPlan( &status, &pinv );
  TestStatus( &status, CODES( 0 ), 1 );

  LALSDestroyVector( &status, &hvec );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCDestroyVector( &status, &Hvec );
  TestStatus( &status, CODES( 0 ), 1 );

  LALSDestroyVector( &status, &Pvec );
  TestStatus( &status, CODES( 0 ), 1 );

  LALSDestroyVectorSequence( &status, &hseq );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCDestroyVectorSequence( &status, &Hseq );
  TestStatus( &status, CODES( 0 ), 1 );

  LALSDestroyVectorSequence( &status, &Pseq );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCheckMemoryLeaks();
  return 0;
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
