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
#define INJECTTESTC_EOVER 6

#define INJECTTESTC_MSGENORM "Normal exit"
#define INJECTTESTC_MSGESUB  "Subroutine failed"
#define INJECTTESTC_MSGEARG  "Error parsing arguments"
#define INJECTTESTC_MSGEVAL  "Input argument out of valid range"
#define INJECTTESTC_MSGEFILE "Could not open file"
#define INJECTTESTC_MSGEMEM  "Out of memory"
#define INJECTTESTC_MSGEOVER "Detector ADC files overlap in time"
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
#include <lal/LALConstants.h>
#include <lal/SeqFactories.h>
#include <lal/Sort.h>
#include <lal/Random.h>
#include <lal/Inject.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GeneratePPNInspiral.h>
#include "StreamInput.h"

NRCSID( INJECTTESTC, "$Id$" );

/* Default parameter settings. */
int lalDebugLevel = 0;

/* Some global constants. */
#define NAMELENGTH (1024) /* maximum length of a file name */
#define NPARAM (4)        /* number of source parameters in sourcefile */

/* Range of random initial offset (nanoseconds). */
#define INIT_OFFSET_MIN   (-30000000000L)
#define INIT_OFFSET_RANGE (30000000000L)

/* Range of random offset between signals (nanoseconds). */
#define OFFSET_MIN   (30000000000)
#define OFFSET_RANGE (30000000000)

/* Range of randomly-generated masses. */
#define M_MIN   (1.0)
#define M_RANGE (1.0)

/* Parameters of internally-generated data segments. */
/*#define EPOCH   (662342400000000000L)  start of third millenium, UTC */
#define EPOCH   (315187200000000000L) /* about Jan. 1, 1990 */
#define DELTAT  (0.0009765625)
#define NPTS_T  (65536)
#define SIGMASQ (10.0)

/* Parameters of internally-generated transfer function. */
#define F0     (50.0)
#define DELTAF (400.0)
#define NPTS_F (2)
#define TRANS  (1.0e17)

/* Other waveform parameters. */
#define FSTART   (40.0)
#define FSTOP    (600.0)
#define DTSAMPLE (0.01)

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

/* A function to convert INT8 nanoseconds to LIGOTimeGPS. */
void
I8ToLIGOTimeGPS( LIGOTimeGPS *output, INT8 input );


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
  UINT4 sourceLength;   /* number of lines of sourcefile */

  /* Other global variables. */
  INT4Vector *idx = NULL;      /* index ordering input files */
  INT4 *id;                    /* pointer to idx->data */
  RandomParams *params = NULL; /* parameters of pseudorandom sequence */
  INT8 tInitial, tFinal; /* earliest and latest possible inject times */
  INT8 t, tEnd;          /* start and stop time of current signal */
  UINT4 nOut, nSource;   /* current source and output file numbers */
  DetectorResponse detector; /* the detector in question */
  INT2TimeSeries output;     /* detector ACD output */
  CHAR outname[NAMELENGTH];  /* name of current output file */


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

  /* First, set up fields of output and transfer functions. */
  output.data = NULL;
  detector.transfer = (COMPLEX8FrequencySeries *)
    LALMalloc( sizeof(COMPLEX8FrequencySeries) );
  if ( !(detector.transfer) ) {
    ERROR( INJECTTESTC_EMEM, INJECTTESTC_MSGEMEM, 0 );
    return INJECTTESTC_EMEM;
  }
  detector.transfer->data = NULL;


  /* Read infile and store the data in structures. */
  if ( infile ) {
    UINT4 vectorLength; /* length of each line */

    if ( ( fp = fopen( infile, "r" ) ) == NULL ) {
      ERROR( INJECTTESTC_EFILE, INJECTTESTC_MSGEFILE, infile );
      return INJECTTESTC_EFILE;
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
      return INJECTTESTC_EMEM;
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

    /* Sort start and stop times, and make sure they don't overlap. */
    {
      REAL8Vector *start = NULL;
      SUB( LALDCreateVector( &stat, &start, inLength ), &stat );
      for ( i = 0; i < inLength; i++ )
	start->data[i] = (REAL8)( tStart[i] );
      SUB( LALI4CreateVector( &stat, &idx, inLength ), &stat );
      SUB( LALDHeapIndex( &stat, idx, start ), &stat );
      SUB( LALDDestroyVector( &stat, &start ), &stat );
      id = idx->data;
      for ( i = 0; i < inLength - 1; i++ )
	if ( tStart[id[i+1]] < tStop[id[i]] ) {
	  ERROR( INJECTTESTC_EOVER, INJECTTESTC_MSGEOVER, 0 );
	  return INJECTTESTC_EOVER;
	}
      tInitial = tStart[id[0]];
      tFinal = tStop[id[inLength-1]];
    }

    /* Read first input file. */
    {
      CHAR fname[NAMELENGTH];           /* input file name */
      INT2VectorSequence *input = NULL; /* input data */
      sscanf( inImage->data + id[0]*inImage->vectorLength, "%s",
	      fname );
      if ( ( fp = fopen( fname, "r" ) ) == NULL ) {
	ERROR( INJECTTESTC_EFILE, INJECTTESTC_MSGEFILE, fname );
	return INJECTTESTC_EFILE;
      }
      SUB( LALI2ReadVectorSequence( &stat, &input, fp ), &stat );
      fclose( fp );
      SUB( LALI2CreateVector( &stat, &(output.data), input->length ),
	   &stat );
      memcpy( output.data->data, input->data,
	      input->length*sizeof(INT2) );
      SUB( LALI2DestroyVectorSequence( &stat, &input ), &stat );
      I8ToLIGOTimeGPS( &(output.epoch), tStart[id[0]] );
      output.deltaT = deltaT[id[0]];
      sprintf( outname, "%s.sim", fname );
    }
  }

  /* If there's no infile, set parameters from #defined constants. */
  else {
    inLength = 1;
    id = (UINT4 *)LALMalloc( sizeof(UINT4) );
    tStart = (INT8 *)LALMalloc( sizeof(INT8) );
    tStop = (INT8 *)LALMalloc( sizeof(INT8) );
    if ( !( tStart && tStop ) ) {
      ERROR( INJECTTESTC_EMEM, INJECTTESTC_MSGEMEM, 0 );
      return INJECTTESTC_EMEM;
    }
    *id = 0;
    *tStart = tInitial = EPOCH;
    *tStop = tFinal = EPOCH + (1000000000L)*(INT8)( NPTS_T*DELTAT );
  }

  /* Read transferfile and store the data in structures. */
  if ( transferfile ) {
    UINT4 vectorLength; /* length of each line */

    if ( ( fp = fopen( transferfile, "r" ) ) == NULL ) {
      ERROR( INJECTTESTC_EFILE, INJECTTESTC_MSGEFILE, transferfile );
      return INJECTTESTC_EFILE;
    }
    SUB( LALCHARReadVectorSequence( &stat, &transferImage, fp ),
	 &stat );
    fclose( fp );

    /* Create lists of epochs, and start and sample frequencies. */
    transferLength = transferImage->length;
    vectorLength = transferImage->vectorLength;
    epoch = (INT8 *)LALMalloc( transferLength*sizeof(INT8) );
    f0 = (REAL8 *)LALMalloc( transferLength*sizeof(REAL8) );
    deltaF = (REAL8 *)LALMalloc( transferLength*sizeof(REAL8) );
    if ( !( epoch && f0 && deltaF ) ) {
      ERROR( INJECTTESTC_EMEM, INJECTTESTC_MSGEMEM, 0 );
      return INJECTTESTC_EMEM;
    }
    for ( i = 0; i < transferLength; i++ )
      if ( 3 > sscanf( transferImage->data + i*vectorLength,
		       "%*s %Li %lf %lf\n", epoch + i, f0 + i,
		       deltaF + i ) )
	transferLength = i;
  }


  /* Read sourcelist to an array. */
  if ( sourcefile ) {
    if ( ( fp = fopen( sourcefile, "r" ) ) == NULL ) {
      ERROR( INJECTTESTC_EFILE, INJECTTESTC_MSGEFILE, sourcefile );
      return INJECTTESTC_EFILE;
    }
    SUB( LALSReadVectorSequence( &stat, &sourceList, fp ), &stat );
    fclose( fp );
    sourceLength = sourceList->length;
  } else {
    sourceLength = (UINT4)( -1 );
  }


  /* Generate random parameters and initial start time. */
  {
    REAL4 frac; /* random number between 0 and 1 */
    SUB( LALCreateRandomParams( &stat, &params, seed ), &stat );
    SUB( LALUniformDeviate( &stat, &frac, params ), &stat );
    t = tInitial + INIT_OFFSET_MIN + (INT8)( INIT_OFFSET_RANGE*frac );
  }
  /* Generate random initial data, if necessary. */
  if ( !infile ) {
    REAL4Vector *deviates = NULL; /* Gaussian random input */

    SUB( LALI2CreateVector( &stat, &(output.data), NPTS_T ), &stat );
    SUB( LALSCreateVector( &stat, &deviates, NPTS_T ), &stat );
    SUB( LALNormalDeviates( &stat, deviates, params ), &stat );
    for ( i = 0; i < NPTS_T; i++ )
      output.data->data[i] = (INT2)( SIGMASQ*deviates->data[i] );
    SUB( LALSDestroyVector( &stat, &deviates ), &stat );
    I8ToLIGOTimeGPS( &(output.epoch), *tStart );
    output.deltaT = DELTAT;
    if ( outfile )
      sprintf( outname, "%s", outfile );
    else
      outname[0] = '\0';
  }


  /* Stub to specify the detector. */
  detector.latitude = 0.0;
  detector.longitude = 0.0;


  /*******************************************************************
   *                                                                 *
   * Start injecting sources.                                        *
   *                                                                 *
   *******************************************************************/

  nSource = 0;
  nOut = 0;
  while ( ( nOut < inLength ) && ( nSource < sourceLength ) &&
	  ( t < tFinal ) ) {
    CoherentGW waveform;              /* gravitational waveform */
    REAL4TimeVectorSeries a, phi;     /* fields of waveform */
    REAL4TimeSeries signal;           /* detector response */

    /* Compute the waveform. */
    {
      GalacticInspiralParamStruc position; /* location of source */
      PPNParamStruc ppnParams;             /* parameters of source */

      /* Get the source position. */
      if ( sourcefile ) {
	position.rho = sourceList->data[NPARAM*nSource];
	position.z = sourceList->data[NPARAM*nSource + 1];
	SUB( LALUniformDeviate( &stat, &(position.lGal), params ),
	     &stat );
	position.lGal *= LAL_TWOPI;
	position.m1 = sourceList->data[NPARAM*nSource + 2];
	position.m2 = sourceList->data[NPARAM*nSource + 3];
      } else {
	position.rho = 0.0;
	position.z = 0.0;
	position.lGal = 0.0;
	SUB( LALUniformDeviate( &stat, &(position.m1), params ),
	     &stat );
	SUB( LALUniformDeviate( &stat, &(position.m2), params ),
	     &stat );
	position.m1 = M_MIN + M_RANGE*position.m1;
	position.m2 = M_MIN + M_RANGE*position.m2;
      }

      /* Set up input parameters. */
      SUB( LALGetInspiralParams( &stat, &ppnParams, &position,
				 params ), &stat );
      ppnParams.fStartIn = FSTART;
      ppnParams.fStopIn = FSTOP;
      ppnParams.ppn = NULL;
      ppnParams.lengthIn = 0;
      a.deltaT = phi.deltaT = DTSAMPLE;
      a.data = phi.data = NULL;
      I8ToLIGOTimeGPS( &(a.epoch), t );
      phi.epoch = a.epoch;
      waveform.h = NULL;
      waveform.a = &a;
      waveform.phi = &phi;


      /* DEBUG:
      printf( "t=%4i m1=%f m2=%f mTot=%f eta=%f\n",
	      (INT4)( ( t - EPOCH )/1000000000 ), position.m1,
	      position.m2, ppnParams.mTot, ppnParams.eta );
      */


      /* Generate waveform. */
      SUB( LALGeneratePPNInspiral( &stat, &waveform, &ppnParams ),
	   &stat );
      if ( ppnParams.dfdt > 4.0 ) {
	INFO( "Waveform sampling interval may be too large." );
      }

      /* Determine end of generated waveform. */
      tEnd = t + (INT8)( 1.0e9 * waveform.a->deltaT
			 * waveform.a->data->length ) + 1;
    }


    /* Get the nearest transfer function. */
    if ( transferfile ) {
      CHAR fname[NAMELENGTH]; /* name of nearest transfer fn file */
      INT8 tDiff1 = abs( t - epoch[0] ); /* epoch of closest file */
      INT4 k = 0;                        /* index of nearest file */
      REAL4VectorSequence *function = NULL; /* transfer function */
      for ( i = 1; i < transferLength; i++ ) {
	INT8 tDiff2 = abs( t - epoch[i] );
	if ( tDiff2 < tDiff1 ) {
	  tDiff1 = tDiff2;
	  k = i;
	}
      }
      sscanf( transferImage->data + k*transferImage->vectorLength,
	      "%s", fname );
      if ( ( fp = fopen( fname, "r" ) ) == NULL ) {
	ERROR( INJECTTESTC_EFILE, INJECTTESTC_MSGEFILE, fname );
	return INJECTTESTC_EFILE;
      }
      SUB( LALSReadVectorSequence( &stat, &function, fp ), &stat );
      fclose( fp );
      SUB( LALCCreateVector( &stat, &(detector.transfer->data),
			     function->length ), &stat );
      memcpy( detector.transfer->data->data, function->data,
	      2*function->length*sizeof(REAL4) );
      SUB( LALSDestroyVectorSequence( &stat, &function ), &stat );
      I8ToLIGOTimeGPS( &(detector.transfer->epoch), epoch[k] );
      detector.transfer->f0 = f0[k];
      detector.transfer->deltaF = deltaF[k];
    }

    /* If there is no transferfile, make up a transfer function. */
    else {
      I8ToLIGOTimeGPS( &(detector.transfer->epoch), EPOCH );
      detector.transfer->f0 = F0;
      detector.transfer->deltaF = DELTAF;
      SUB( LALCCreateVector( &stat, &(detector.transfer->data),
			     NPTS_F ), &stat );
      for ( i = 0; i < NPTS_F; i++ ) {
	detector.transfer->data->data[i].re = TRANS;
	detector.transfer->data->data[i].im = 0.0;
      }
    }


    /* DEBUG:
    if ( nSource == 1 ) {
      FILE *fptest = fopen( "test.wav", "w" );
      for ( i = 0; i < waveform.a->data->length; i++ )
	fprintf( fptest, "%10.3e %10.3e %10.3e\n",
		 waveform.a->data->data[2*i],
		 waveform.a->data->data[2*i+1],
		 waveform.phi->data->data[i] );
      fclose( fptest );
    }
    */


    /* Compute detector response. */
    signal.epoch = waveform.a->epoch;
    signal.deltaT = DELTAT;
    signal.data = NULL;
    SUB( LALSCreateVector( &stat, &(signal.data), (UINT4)
			   (waveform.a->data->length*DTSAMPLE/DELTAT) + 1 ),
	 &stat );
    SUB( LALSimulateCoherentGW( &stat, &signal, &waveform, &detector ),
	 &stat );
    SUB( LALCDestroyVector( &stat, &(detector.transfer->data) ),
	 &stat );
    SUB( LALSDestroyVectorSequence( &stat, &(waveform.a->data) ),
	 &stat );
    SUB( LALSDestroyVectorSequence( &stat, &(waveform.phi->data) ),
	 &stat );



    /* DEBUG:
    if ( nSource == 1 ) {
      FILE *fptest = fopen( "test.sig", "w" );
      for ( i = 0; i < signal.data->length; i++ )
	fprintf( fptest, "%10.3e\n", signal.data->data[i] );
      fclose( fptest );
    }
    */


    /* Inject signal into detector output. */
    SUB( LALSI2InjectTimeSeries( &stat, &output, &signal, params ),
	 &stat );
    while ( ( nOut < inLength ) && ( tEnd > tStop[id[nOut]] ) ) {

      /* Write data to output file. */
      if ( outname[0] != '\0' ) {
	if ( ( fp = fopen( outname, "w" ) ) == NULL ) {
	  ERROR( INJECTTESTC_EFILE, INJECTTESTC_MSGEFILE, outname );
	  return INJECTTESTC_EFILE;
	}
	for ( i = 0; i < output.data->length; i++ )
	  fprintf( fp, "%6i\n", output.data->data[i] );
	fclose( fp );
	outname[0] = '\0';
      }
      nOut++;

      /* Read data from input file. */
      if ( nOut < inLength ) {
	CHAR fname[NAMELENGTH];           /* input file name */
	INT2VectorSequence *input = NULL; /* input data */
	SUB( LALI2DestroyVector( &stat, &(output.data) ), &stat );
	sscanf( inImage->data + id[nOut]*inImage->vectorLength, "%s",
		fname );
	if ( ( fp = fopen( fname, "r" ) ) == NULL ) {
	  ERROR( INJECTTESTC_EFILE, INJECTTESTC_MSGEFILE, fname );
	  return INJECTTESTC_EFILE;
	}
	SUB( LALI2ReadVectorSequence( &stat, &input, fp ), &stat );
	fclose( fp );
	SUB( LALI2CreateVector( &stat, &(output.data),
				input->length ), &stat );
	memcpy( output.data->data, input->data,
		input->length*sizeof(INT2) );
	SUB( LALI2DestroyVectorSequence( &stat, &input ), &stat );
	I8ToLIGOTimeGPS( &(output.epoch), tStart[id[nOut]] );
	output.deltaT = deltaT[id[nOut]];
	sprintf( outname, "%s.sim", fname );
      }
    }

    /* Signal injection is done, so move on to the next source. */
    SUB( LALSDestroyVector( &stat, &(signal.data) ), &stat );
    nSource++;
    {
      REAL4 frac; /* random number between 0 and 1 */
      SUB( LALUniformDeviate( &stat, &frac, params ), &stat );
      t = tEnd + OFFSET_MIN + (INT8)( OFFSET_RANGE*frac );
    }
  }


  /*******************************************************************
   *                                                                 *
   * Clean up.                                                       *
   *                                                                 *
   *******************************************************************/

  /* Print remaining output data. */
  if ( outname[0] != '\0' ) {
    if ( ( fp = fopen( outname, "w" ) ) == NULL ) {
      ERROR( INJECTTESTC_EFILE, INJECTTESTC_MSGEFILE, outname );
      return INJECTTESTC_EFILE;
    }
    for ( i = 0; i < output.data->length; i++ )
      fprintf( fp, "%6i\n", output.data->data[i] );
    fclose( fp );
  }

  /* Destroy remaining memory. */
  SUB( LALDestroyRandomParams( &stat, &params ), &stat );
  SUB( LALI2DestroyVector( &stat, &(output.data) ), &stat );
  LALFree( detector.transfer );

  if ( infile ) {
    SUB( LALCHARDestroyVectorSequence( &stat, &inImage ), &stat );
    SUB( LALI4DestroyVector( &stat, &idx ), &stat );
    LALFree( deltaT );
  } else {
    LALFree( id );
  }
  LALFree( tStart );
  LALFree( tStop );

  if ( transferfile ) {
    SUB( LALCHARDestroyVectorSequence( &stat, &transferImage ),
	 &stat );
    LALFree( epoch );
    LALFree( f0 );
    LALFree( deltaF );
  }

  if ( sourcefile ) {
    SUB( LALSDestroyVectorSequence( &stat, &sourceList ), &stat );
  }

  LALCheckMemoryLeaks();
  INFO( INJECTTESTC_MSGENORM );
  return INJECTTESTC_ENORM;
}


/* A function to convert INT8 nanoseconds to LIGOTimeGPS. */
void
I8ToLIGOTimeGPS( LIGOTimeGPS *output, INT8 input )
{
  INT4 s = input / 1000000000;
  output->gpsSeconds = s;
  output->gpsNanoSeconds = input - 1000000000*s;
  return;
}
