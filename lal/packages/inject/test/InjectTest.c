/*********************************** <lalVerbatim file="InjectTestCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{InjectTest.c}}
\label{ss:InjectTest.c}

Injects an inspiral signal into detector noise.

\subsubsection*{Usage}
\begin{verbatim}
InjectTest [-s sourcefile] [-t transferfile] [-i infile | -o outfile] [-d debuglevel] [-r randomseed]
\end{verbatim}

\subsubsection*{Description}

This program generates Galactic inspiral waveform signals and injects
them into ADC data.  The following option flags are accepted:
\begin{itemize}
\item[\texttt{-s}] Reads Galactic source positions and masses from the
file \verb@sourcefile@, whose format is given below.  If not
specified, masses are chosen randomly between 1 and 2$M_\odot$, and
are placed at the Galactic core.
\item[\texttt{-t}] Specifies an index file \verb@transferfile@ giving
the names and epochs of data files containing the instrument transfer
function (ADC counts per unit strain).  If not specified, a fake
internal whitening filter is used.
\item[\texttt{-i}] Specifies an index file \verb@infile@ giving the
names and timespans of data files containing the instrument output.
If not specified, a single stretch of Gaussian noise will be generated
internally and injected.
\item[\texttt{-o}] Used only if a \verb@-i@ option was \emph{not}
given, this option specifies an output file \verb@outfile@ where the
simulated data will be written.  If the \verb@-i@ option \emph{was}
given, the injected data will be written to files with names derived
from the input data files.  If neither \verb@-i@ nor \verb@-o@ is
specified, the injection routines are exercised but the simulated data
is discarded.
\item[\texttt{-d}] Sets the debug level to \verb@debuglevel@.  If not
specified, level 0 is assumed.
\item[\texttt{-r}] Sets the random number seed to \verb@randomseed@.
If not specified, the seed is gerenated from the current time.
\end{itemize}

\paragraph{Format for \texttt{sourcefile}:} The source file consists
of four columns of floating-point numbers, representing the distance
$\rho$ of the system from the Galactic axis in kpc, the height $z$
above (or below) the Galactic plane in kpc, and the two binary masses
in Solar masses.  All are read in as \verb@REAL4@.

\paragraph{Format for \texttt{infile}:} The input file index consists
of four columns.  The first is a character string specifying the name
of a file containing ADC data, the second is an \verb@INT8@ specifying
the GPS start time (epoch) of the data set in nanoseconds, the third
is a \verb@UINT4@ specifying the length of the data file, and the
fourth is a \verb@REAL8@ specifying the sampling interval in seconds.

The ADC data files pointed to by \verb@infile@ are a single column of
\verb@INT2@ data giving the ADC output at each time sample.

\paragraph{Format for \texttt{transferfile}:} The transfer function
file index consists of four columns.  The first is a character string
specifying the name of a file containing a transfer function, the
second is an \verb@INT8@ specifying the GPS epoch when the transfer
function was taken (or when it is most valid) in nanoseconds, the
third is a \verb@REAL8@ specifying the starting frequency of the
transfer function data in hertz, and the fourth is a \verb@REAL8@
specifying the frequency sampling interval in hertz.

The transfer function data files pointed to by \verb@transferfile@
consist of two columns of \verb@REAL4@ data, giving the real and
imaginary parts of the transfer function (in ADC counts per unit
strain) at each frequency sample.

\subsubsection*{Exit codes}
****************************************** </lalLaTeX><lalErrTable> */
#define INJECTTESTC_ENORM 0
#define INJECTTESTC_ESUB  1
#define INJECTTESTC_EARG  2
#define INJECTTESTC_EVAL  3
#define INJECTTESTC_EFILE 4
#define INJECTTESTC_EMEM  5

#define INJECTTESTC_MSGENORM "Normal exit"
#define INJECTTESTC_MSGESUB  "Subroutine failed"
#define INJECTTESTC_MSGEARG  "Error parsing arguments"
#define INJECTTESTC_MSGEVAL  "Input argument out of valid range"
#define INJECTTESTC_MSGEFILE "Could not open file"
#define INJECTTESTC_MSGEMEM  "Out of memory"
/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Algorithm}

The general structure of the routine is as follows, assuming
\emph{all} input file options \verb@-s@, \verb@-t@, and \verb@-i@ have
been specified.
\begin{enumerate}
\item The file \verb@infile@ is read and the earliest start time and
latest stop time are determined.  The start time for the first
injected signal is offset \emph{backwards} from the earliest data
start time by a random amount.
\item The file \verb@sourcefile@ is read, and the first set of source
parameters is sent to \verb@LALGetInspiralParams()@ and
\verb@GeneratePPNInspiral@ to generate a waveform.  The total timespan
of the waveform is determined.
\item The file \verb@transferfile@ is read, and the transfer function
file nearest to the midpoint of the waveform is determined.  That file
is then read in, and used with \verb@LALSimulateCoherentGW()@ to
generate a simulated detector output.
\item From the list in \verb@infile@, the set of data files that
overlap with the waveform is determined.  These files are read in and
injected using the function \verb@LALSI2InjectTimeSeries()@.
\item The injected data are written to output files whose names are
derived from the input data files by adding a \verb@.sim@ extension.
The last of the datasets is kept in memory and not written, since
subsequent inspiral injections may also overlap with it.
\item The start time for the next signal is determined by adding a
random amount to the end of the previous signal, and steps 2--5 are
repeated.  This continues until either the latest stop time is passed,
or the sourcelist in \verb@sourcefile@ is exhausted.
\item The last dataset remaining in memory is written to its output
file.
\end{enumerate}

The negative time offset for the first signal is a random number
between 0 and 5~minutes.  The time offset between signals is a random
number between 1 and 2~minutes.  These are set by \verb@#define@d
constants.

If no \verb@infile@ is specified, a 1024~second 1024~Hz ADC datastream
is generated internally and populated with white Gaussian noise having
a variance $\sigma^2=10$.  These numbers can also be changed by
changing \verb@#define@d constants.  The injected data is written to a
file only if an \verb@outfile@ was specified with the \verb@-o@
option.

If no \verb@sourcefile@ is specified, an unlimited number of source
parameters are generated randomly.  Each source is positioned at the
Galactic core (but with random inclination, phase, and polarization to
ensure a varying amplitude), with masses uniformly and randomly
distributed between 1 and 2$M_\odot$.

If no \verb@transferfile@ is specified, a simple white bandpass from
40~Hz to 600~Hz is used.  Obviously this is not partiucularly
realistic, and I may improve upon it, but you should really just feed
in your favourite transfer function with the \verb@-t@ option.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALPrintError()                 LALCheckMemoryLeaks()
LALMalloc()                     LALFree()
\end{verbatim}

\subsubsection*{Notes}

Under development.  At present it doesn't do anything; it just exits.

\vfill{\footnotesize\input{InjectTestCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/Random.h>
#include <lal/Inject.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GeneratePPNInspiral.h>
#include "StreamInput.h"

NRCSID( INJECTTESTC, "$Id$" );

/* Default parameter settings. */
int lalDebugLevel = 0;

/* Range of random initial offset (nanoseconds). */
#define INIT_OFFSET_MIN   (-300000000000L)
#define INIT_OFFSET_RANGE (300000000000L)

/* Range of random offset between signals (nanoseconds). */
#define OFFSET_MIN   (60000000000)
#define OFFSET_RANGE (60000000000)

/* Parameters of internally-generated data segments. */
#define EPOCH   (662342400000000000L) /* start of third millenium, UTC */
#define DELTAT  (0.0009765625)
#define NPTS_T  (1048576)
#define SIGMASQ (10.0)

/* Parameters of internally-generated transfer function. */
#define F0     (40.0)
#define DELTAF (560)
#define NPTS_F (2)
#define TRANS  (1.0)

/* Usage format string. */
#define USAGE "Usage: %s [-s sourcefile] [-d debuglevel]\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, INJECTTESTC, statement ? statement :      \
                 "", (msg) );                                        \
}                                                                    \
while (0)

#define INFO( statement )                                            \
do                                                                   \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 INJECTTESTC, (statement) );                         \
}                                                                    \
while (0)

#define SUB( func, statusptr )                                       \
do                                                                   \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( INJECTTESTC_ESUB, INJECTTESTC_MSGESUB,                      \
         "Function call \"" #func "\" failed:" );                    \
  return INJECTTESTC_ESUB;                                           \
}                                                                    \
while (0)

#define CHECKVAL( val, lower, upper )                                \
do                                                                   \
if ( ( (val) <= (lower) ) || ( (val) > (upper) ) )                   \
{                                                                    \
  ERROR( INJECTTESTC_EVAL, INJECTTESTC_MSGEVAL,                      \
         "Value of " #val " out of range:" );                        \
  LALPrintError( #val " = %f, range = (%f,%f]\n", (REAL8)(val),      \
                 (REAL8)(lower), (REAL8)(upper) );                   \
  return INJECTTESTC_EVAL;                                           \
}                                                                    \
while (0)

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

int
main(int argc, char **argv)
{
  /* Command-line parsing variables. */
  int arg;                   /* command-line argument counter */
  static LALStatus stat;     /* status structure */
  CHAR *sourcefile = NULL;   /* name of sourcefile */
  CHAR *transferfile = NULL; /* name of transferfile */
  CHAR *infile = NULL;       /* name of infile */
  CHAR *outfile = NULL;      /* name of outfile */
  INT4 seed = 0;             /* random number seed */

  /* Input file reading variables. */
  FILE *fp;             /* generic file pointer */
  UINT4 i;              /* generic index over file lines */
  CHARVectorSequence *inImage = NULL; /* character image of infile */
  UINT4 inLength;       /* number of lines of infile */
  REAL8 *deltaT;        /* array of sampling rates for each input file */
  INT8 *tStart, *tStop; /* arrays of start and stop times */
  CHARVectorSequence *transferImage = NULL; /* image of transferfile */
  UINT4 transferLength; /* number of lines of transferfile */
  INT8 *epoch;          /* array of epochs for each transfer function file */
  REAL8 *f0, *deltaF;   /* arrays giving frequency samples for each file */
  REAL4VectorSequence *sourceList = NULL; /* contents of sourcefile */

  /* Other global variables. */
  RandomParams *params = NULL; /* parameters of pseudorandom sequence */
  INT8 tInitial, tFinal; /* earliest and latest possible inject times */
  INT8 t;                /* current injection time */


  /*******************************************************************
   *                                                                 *
   * Parse command-line arguments.  arg stores the current position. *
   *                                                                 *
   *******************************************************************/

  arg = 1;
  while ( arg < argc ) {
    /* Parse source file option. */
    if ( !strcmp( argv[arg], "-s" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	sourcefile = argv[arg++];
      }else{
	ERROR( INJECTTESTC_EARG, INJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INJECTTESTC_EARG;
      }
    }
    /* Parse transfer file option. */
    else if ( !strcmp( argv[arg], "-t" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	transferfile = argv[arg++];
      }else{
	ERROR( INJECTTESTC_EARG, INJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INJECTTESTC_EARG;
      }
    }
    /* Parse input file option. */
    else if ( !strcmp( argv[arg], "-i" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	infile = argv[arg++];
      }else{
	ERROR( INJECTTESTC_EARG, INJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INJECTTESTC_EARG;
      }
    }
    /* Parse output file option. */
    else if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	outfile = argv[arg++];
      }else{
	ERROR( INJECTTESTC_EARG, INJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INJECTTESTC_EARG;
      }
    }
    /* Parse debug level option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	lalDebugLevel = atoi( argv[arg++] );
      }else{
	ERROR( INJECTTESTC_EARG, INJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INJECTTESTC_EARG;
      }
    }
    /* Parse random seed option. */
    else if ( !strcmp( argv[arg], "-r" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	seed = atoi( argv[arg++] );
      }else{
	ERROR( INJECTTESTC_EARG, INJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INJECTTESTC_EARG;
      }
    }
    /* Check for unrecognized options. */
    else if ( argv[arg][0] == '-' ) {
      ERROR( INJECTTESTC_EARG, INJECTTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return INJECTTESTC_EARG;
    }
  } /* End of argument parsing loop. */


  /*******************************************************************
   *                                                                 *
   * Read input files, and store data to arrays.                     *
   *                                                                 *
   *******************************************************************/

  /* Read infile and store the data in structures. */
  if ( infile ) {
    UINT4 vectorLength; /* length of each line */

    if ( ( fp = fopen( infile, "r" ) ) == NULL ) {
      ERROR( INJECTTESTC_EFILE, INJECTTESTC_MSGEFILE, infile );
    }
    SUB( LALCHARReadVectorSequence( &stat, &inImage, fp ), &stat );
    fclose( fp );

    /* Create lists of start, stop, and sample times. */
    inLength = inImage->length;
    vectorLength = inImage->vectorLength;
    tStart = (INT8 *)LALMalloc( inLength*sizeof(INT8) );
    tStop = (INT8 *)LALMalloc( inLength*sizeof(INT8) );
    deltaT = (REAL8 *)LALMalloc( inLength*sizeof(REAL8) );
    if ( !( tStart && tStop && deltaT ) ) {
      ERROR( INJECTTESTC_EMEM, INJECTTESTC_MSGEMEM, 0 );
    }
    for ( i = 0; i < inLength; i++ ) {
      UINT4 length; /* number of points in the data file */
      if ( 3 > sscanf( inImage->data + i*vectorLength,
		       "%*s %Li %u %lf\n", tStart + i, &length,
		       deltaT + i ) )
	inLength = i;
      else
	tStop[i] = tStart[i] + (INT8)( length*deltaT[i] );
    }

    /* Determine the earliest start time and latest stop time. */
    tInitial = tStart[0];
    tFinal = tStop[0];
    for ( i = 0; i < inLength; i++ ) {
      if ( tStart[i] < tInitial )
	tInitial = tStart[i];
      if ( tStop[i] > tFinal )
	tFinal = tStop[i];
    }
  }

  /* If there's no infile, set parameters from #defined constants. */
  else {
    inLength = 1;
    tStart = (INT8 *)LALMalloc( sizeof(INT8) );
    tStop = (INT8 *)LALMalloc( sizeof(INT8) );
    deltaT = (REAL8 *)LALMalloc( sizeof(REAL8) );
    if ( !( tStart && tStop && deltaT ) ) {
      ERROR( INJECTTESTC_EMEM, INJECTTESTC_MSGEMEM, 0 );
    }
    *tStart = EPOCH;
    *tStop = EPOCH + (INT8)( NPTS_T*DELTAT );
    *deltaT = DELTAT;
  }


  /* Read transferfile and store the data in structures. */
  if ( transferfile ) {
    UINT4 vectorLength; /* length of each line */

    if ( ( fp = fopen( transferfile, "r" ) ) == NULL ) {
      ERROR( INJECTTESTC_EFILE, INJECTTESTC_MSGEFILE, transferfile );
    }
    SUB( LALCHARReadVectorSequence( &stat, &transferImage, fp ),
	 &stat );
    fclose( fp );

    /* Create lists of epochs, and start and sample frequencies. */
    transferLength = transferImage->length;
    vectorLength = inImage->vectorLength;
    epoch = (INT8 *)LALMalloc( transferLength*sizeof(INT8) );
    f0 = (REAL8 *)LALMalloc( transferLength*sizeof(REAL8) );
    deltaF = (REAL8 *)LALMalloc( transferLength*sizeof(REAL8) );
    if ( !( epoch && f0 && deltaF ) ) {
      ERROR( INJECTTESTC_EMEM, INJECTTESTC_MSGEMEM, 0 );
    }
    for ( i = 0; i < transferLength; i++ )
      if ( 3 > sscanf( transferImage->data + i*vectorLength,
		       "%*s %Li %lf %lf\n", epoch + i, f0 + i,
		       deltaF + i ) )
	transferLength = i;
  }

  /* If there's no transferfile, set parameters from constants. */
  else {
    transferLength = 1;
    epoch = (INT8 *)LALMalloc( sizeof(INT8) );
    f0 = (REAL8 *)LALMalloc( sizeof(REAL8) );
    deltaF = (REAL8 *)LALMalloc( sizeof(REAL8) );
    if ( !( epoch && f0 && deltaF ) ) {
      ERROR( INJECTTESTC_EMEM, INJECTTESTC_MSGEMEM, 0 );
    }
    *epoch = EPOCH;
    *f0 = F0;
    *deltaF = DELTAF;
  }


  /* Read sourcelist to an array. */
  if ( sourcefile ) {
    if ( ( fp = fopen( sourcefile, "r" ) ) == NULL ) {
      ERROR( INJECTTESTC_EFILE, INJECTTESTC_MSGEFILE, sourcefile );
    }
    SUB( LALSReadVectorSequence( &stat, &sourceList, fp ), &stat );
    fclose( fp );
  }

  /* Generate random parameters and initial start time. */
  {
    REAL4 frac; /* random number between 0 and 1 */
    SUB( LALCreateRandomParams( &stat, &params, seed ), &stat );
    SUB( LALUniformDeviate( &stat, &frac, params ), &stat );
    t = tInitial + INIT_OFFSET_MIN + (INT8)( INIT_OFFSET_RANGE*frac );
  }





  /*LALCheckMemoryLeaks();*/
  INFO( INJECTTESTC_MSGENORM );
  return INJECTTESTC_ENORM;
}
