/******************************** <lalVerbatim file="IIRFilterTestCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{IIRFilterTest.c}}
\label{s:IIRFilterTest.c}

Tests the routines in \verb@IIRFilter.h@.

\subsubsection*{Usage}
\begin{verbatim}
IIRFilterTest [-h] [-o] [-t] [-d debuglevel] [-n npts] [-r reps] [-w freq]
\end{verbatim}

\subsubsection*{Description}

This program generates a time series vector, and passes it through a
third-order Butterworth low-pass filter.  By default, running this
program with no arguments simply passes an impulse function to the
filter routines, producing no output.  All filter parameters are set
from \verb@#define@d constants.

The \verb@-h@ option causes a usage message to print to \verb@stderr@.

The \verb@-o@ option tells the program to print the input and output
vectors to data files: \verb@out.0@ stores the initial impulse,
\verb@out.1@ the response computed using \verb@SIIRFilter()@,
\verb@out.2@ the response computed using \verb@IIRFilterREAL4()@,
\verb@out.3@ the response from \verb@IIRFilterREAL4Vector()@, and
\verb@out.4@ the response from \verb@IIRFilterREAL4VectorR()@.

The \verb@-t@ option causes \verb@IIRFilterTest@ to fill the time
series vector with Gaussian random deviates, and prints execution
times for the various filter subroutines to \verb@stdout@.  To
generate useful timing data, the default size of the time vector is
increased (unless explicitly set, below).

The \verb@-d@ option changes the debug level from 0 to the specified
value \verb@debuglevel@.

The \verb@-n@ option sets the size of the time vectors to
\verb@npts@.

The \verb@-r@ option applies each filter to the data \verb@reps@
times.

The \verb@-w@ option sets the characteristic frequency of the filter
to \verb@freq@ in the $w$-plane (described in \verb@ZPGFilter.h@).

\subsubsection*{Exit codes}
****************************************** </lalLaTeX><lalErrTable> */
#define IIRFILTERTESTC_ENORM 0
#define IIRFILTERTESTC_ESUB  1
#define IIRFILTERTESTC_EARG  2
#define IIRFILTERTESTC_EBAD  3
#define IIRFILTERTESTC_EFILE 4

#define IIRFILTERTESTC_MSGENORM "Normal exit"
#define IIRFILTERTESTC_MSGESUB  "Subroutine failed"
#define IIRFILTERTESTC_MSGEARG  "Error parsing arguments"
#define IIRFILTERTESTC_MSGEBAD  "Bad argument value"
#define IIRFILTERTESTC_MSGEFILE "Could not create output file"
/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Algorithm}

A third-order Butterworth low-pass filter is defined by the following
power response function:
$$
|T(w)|^2 = \frac{1}{1+(w/w_c)^6}\;,
$$
where the frequency parameter $w=\tan(\pi f\Delta t)$ maps the Nyquist
interval onto the entire real axis, and $w_c$ is the (transformed)
characteristic frequency set by the \verb@-w@ option above.  A stable
time-domain filter has poles with $|z|<1$ in the $z$-plane
representation, which means that the $w$-plane poles should have
$\Im(w)>0$.  We construct a transfer function using only the
positive-imaginary poles of the power response function:
$$
T(w) = \frac{iw_c^3}{\prod_{k=0}^2
	\left[w-w_c e^{(2k+1)i\pi/6}\right]}\;,
$$
where we have chosen a phase coefficient $i$ in the numerator in order
to get a purely real DC response.  This ensures that the $z$-plane
representation of the filter will have a real gain, resulting in a
physically-realizable IIR filter.

The poles and gain of the transfer function $T(w)$ are simply read off
of the equation above, and are stored in a \verb@COMPLEX8ZPGFilter@.
This is transformed from the $w$-plane to the $z$-plane representation
using \verb@LALWToZCOMPLEX8ZPGFilter()@, and then used to create an
IIR filter with \verb@LALCreateREAL4IIRFilter()@.  This in turn is
used by the routines \verb@LALSIIRFilter()@,
\verb@LALIIRFilterREAL4()@, \verb@LALIIRFilterREAL4Vector()@, and
\verb@LALIIRFilterREAL4VectorR()@ to filter a data vector containing
either a unit impulse or white Gaussian noise (for more useful timing
information).

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALPrintError()
LALSCreateVector()              LALSDestroyVector()
LALCreateRandomParams()         LALDestroyRandomParams()
LALCreateCOMPLEX8ZPGFilter()    LALDestroyCOMPLEX8ZPGFilter()
LALCreateREAL4IIRFilter()       LALDestroyREAL4IIRFilter()
LALNormalDeviates()
LALWToZCOMPLEX8ZPGFilter()
LALSIIRFilter()
LALIIRFilterREAL4()
LALIIRFilterREAL4Vector()
LALIIRFilterREAL4VectorR()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{BandPassTestCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/IIRFilter.h>
#include <lal/ZPGFilter.h>

NRCSID(IIRFILTERTESTC,"$Id$");

/* Default parameters. */
INT4 lalDebugLevel=0;
#define NPTS   4096    /* Default length of time series */
#define NPTS_T 4194304 /* Length of time series for timing runs */
#define WC     0.01    /* Characteristic frequency in w-plane */
#define REPS   1       /* Number of repeated filterings */

/* Output filenames. */
#define OUTFILE0 "out.0"
#define OUTFILE1 "out.1"
#define OUTFILE2 "out.2"
#define OUTFILE3 "out.3"
#define OUTFILE4 "out.4"

/* Mathematical constant. */
#define SQRT3_2 0.8660254037844386467637231707529361L

/* Usage format string. */
#define USAGE "Usage: %s [-d debuglevel] [-n npts] [-r reps] [-w freq] [-o] [-t]\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
		   "        %s %s\n", (code), *argv, __FILE__,       \
		   __LINE__, IIRFILTERTESTC, statement ? statement : \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
		   "        %s\n", *argv, __FILE__, __LINE__,        \
		   IIRFILTERTESTC, (statement) );                    \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( IIRFILTERTESTC_ESUB, IIRFILTERTESTC_MSGESUB,              \
           "Function call \"" #func "\" failed:" );                  \
    return IIRFILTERTESTC_ESUB;                                      \
  }                                                                  \
} while (0)

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

/* Local static functions. */
static void
PrintVector(FILE *fp, REAL4Vector *vector);

int
main(int argc, char **argv)
{
  static LALStatus stat; /* Status pointer for subroutines */
  BOOLEAN print = 0;     /* Whether output will be printed */
  BOOLEAN doTime = 0;    /* Whether filters will be timed */
  BOOLEAN points = 0;    /* Whether number of points is given */
  INT4 npts = NPTS;      /* Number of points in time series */
  INT4 reps = REPS;      /* Number of repeated filterings */
  INT4 arg;              /* argument counter */
  UINT4 i, j;            /* Index counters */
  REAL4 wc = WC;         /* Characteristic w-plane frequency */
  REAL4 *data1;          /* Time series data */
  REAL4 *data2;          /* More time series data */
  REAL4Vector *input1 = NULL; /* A time series input vector */
  REAL4Vector *input2 = NULL; /* Another time series input vector */
  REAL4Vector *output = NULL; /* A time series output vector */
  REAL4IIRFilter *iirFilter = NULL; /* The IIR filter */
  FILE *fp = NULL; /* The output file */
  clock_t start;   /* Clock time before starting to filter */
  clock_t stop;    /* Clock time after filtering */

  /* Parse argument list.  i stores the current position. */
  arg = 1;
  while ( arg < argc ) {
    /* Parse debuglevel option. */
    if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        lalDebugLevel = atoi( argv[arg++] );
      } else {
	ERROR( IIRFILTERTESTC_EARG, IIRFILTERTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return IIRFILTERTESTC_EARG;
      }
    }
    /* Parse number of points. */
    else if ( !strcmp( argv[arg], "-n" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
	npts = atoi( argv[arg++] );
	points = 1;
      } else {
	ERROR( IIRFILTERTESTC_EARG, IIRFILTERTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return IIRFILTERTESTC_EARG;
      }
    }
    /* Parse number of repetitions. */
    else if ( !strcmp( argv[arg], "-r" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
	reps = strtoul( argv[arg++], NULL, 10 );
      } else {
	ERROR( IIRFILTERTESTC_EARG, IIRFILTERTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return IIRFILTERTESTC_EARG;
      }
    }
    /* Parse characteristic frequency. */
    else if ( !strcmp( argv[arg], "-w" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
	wc = atof( argv[arg++] );
      } else {
	ERROR( IIRFILTERTESTC_EARG, IIRFILTERTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return IIRFILTERTESTC_EARG;
      }
    }
    /* Parse timing and output flags. */
    else if ( !strcmp( argv[arg], "-t" ) ) {
      arg++;
      doTime = 1;
    }
    else if ( !strcmp( argv[arg], "-o" ) ) {
      arg++;
      print = 1;
    }
    /* Unrecognized option. */
    else {
      ERROR( IIRFILTERTESTC_EARG, IIRFILTERTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return IIRFILTERTESTC_EARG;
    }
  } /* End of argument parsing loop. */

  /* Check argument values. */
  if ( doTime && !points )
    npts = NPTS_T;
  if ( npts <= 0 ) {
    ERROR( IIRFILTERTESTC_EBAD, IIRFILTERTESTC_MSGEBAD, "npts <= 0:" );
    LALPrintError( USAGE, *argv );
    return IIRFILTERTESTC_EBAD;
  }
  if ( reps <= 0 ) {
    ERROR( IIRFILTERTESTC_EBAD, IIRFILTERTESTC_MSGEBAD, "reps <= 0:" );
    LALPrintError( USAGE, *argv );
    return IIRFILTERTESTC_EBAD;
  }
  if ( wc <= 0.0 ) {
    ERROR( IIRFILTERTESTC_EBAD, IIRFILTERTESTC_MSGEBAD, "freq <= 0:" );
    LALPrintError( USAGE, *argv );
    return IIRFILTERTESTC_EBAD;
  }

  /* Create the time-domain filter. */
  {
    /* First create ZPG filter used to define IIR filter. */
    COMPLEX8ZPGFilter *zpgFilter = NULL;
    SUB( LALCreateCOMPLEX8ZPGFilter( &stat, &zpgFilter, 0, 3 ),
	 &stat );
    zpgFilter->poles->data[0].re = wc*SQRT3_2;
    zpgFilter->poles->data[0].im = wc*0.5;
    zpgFilter->poles->data[1].re = 0.0;
    zpgFilter->poles->data[1].im = wc;
    zpgFilter->poles->data[2].re = -wc*SQRT3_2;
    zpgFilter->poles->data[2].im = wc*0.5;
    zpgFilter->gain.re = 0.0;
    zpgFilter->gain.im = wc*wc*wc;

    /* Now create IIR filter and destroy the ZPG filter. */
    SUB( LALWToZCOMPLEX8ZPGFilter( &stat, zpgFilter ), &stat );
    SUB( LALCreateREAL4IIRFilter( &stat, &iirFilter, zpgFilter ),
	 &stat );
    SUB( LALDestroyCOMPLEX8ZPGFilter( &stat, &zpgFilter ), &stat );
  }

  /* Allocate memory for the time series. */
  SUB( LALSCreateVector( &stat, &input1, npts ), &stat );
  SUB( LALSCreateVector( &stat, &input2, npts ), &stat );
  SUB( LALSCreateVector( &stat, &output, npts ), &stat );

  /* Create the input time series. */
  if ( doTime ) {
    RandomParams *params = NULL; /* Params for random generator */
    SUB( LALCreateRandomParams( &stat, &params, 0 ), &stat );
    INFO( "Begining generation of Gaussian random deviates.\n" );
    SUB( LALNormalDeviates( &stat, input1, params ), &stat );
    INFO( "Finished generating random deviates.\n" );
    SUB( LALDestroyRandomParams( &stat, &params ), &stat );
  } else {
    memset( input1->data, 0, npts*sizeof(REAL4) );
    input1->data[npts/2] = 1.0;
  }
  memcpy( input2->data, input1->data, npts*sizeof(REAL4) );
  if ( print ) {
    if ( !( fp = fopen( OUTFILE0, "w" ) ) ) {
      ERROR( IIRFILTERTESTC_EFILE, IIRFILTERTESTC_MSGEFILE, 0 );
      return IIRFILTERTESTC_EFILE;
    }
    PrintVector(fp,input1);
  }

  /* Filter the time series using SIIRFilter(). */
  if ( doTime ) {
    if ( reps == 1 )
      fprintf( stdout, "Filtering %i points:\n", npts );
    else
      fprintf( stdout, "Filtering %i points %i times:\n", npts,
	       reps );
  }
  j = reps;
  start = clock();
  while ( j-- ) {
    data1 = input1->data;
    data2 = output->data;
    i = npts;
    while ( i-- )
      *(data2++) = LALSIIRFilter( *(data1++), iirFilter );
  }
  stop = clock();
  if ( doTime )
    fprintf( stdout, "Elapsed time for SIIRFilter():            %.2f"
	     " s\n", (double)( stop - start )/CLOCKS_PER_SEC );
  if ( print ) {
    if ( !( fp = fopen( OUTFILE1, "w" ) ) ) {
      ERROR( IIRFILTERTESTC_EFILE, IIRFILTERTESTC_MSGEFILE, 0 );
      return IIRFILTERTESTC_EFILE;
    }
    PrintVector( fp, output );
  }

  /* Filter the time series using IIRFilterREAL4(). */
  j = reps;
  start = clock();
  while ( j-- ) {
    data1 = input1->data;
    data2 = output->data;
    i = npts;
    while ( i-- )
      SUB( LALIIRFilterREAL4( &stat, data2++, *(data1++), iirFilter ),
	   &stat );
  }
  stop = clock();
  if ( doTime )
    fprintf( stdout, "Elapsed time for IIRFilterREAL4():        %.2f"
	     " s\n", (double)( stop - start )/CLOCKS_PER_SEC );
  if( print ) {
    if ( !( fp = fopen( OUTFILE2, "w" ) ) ) {
      ERROR( IIRFILTERTESTC_EFILE, IIRFILTERTESTC_MSGEFILE, 0 );
      return IIRFILTERTESTC_EFILE;
    }
    PrintVector( fp, output );
  }

  /* Filter the time series using IIRFilterREAL4Vector(). */
  j = reps;
  start = clock();
  while ( j-- )
    SUB( LALIIRFilterREAL4Vector( &stat, input1, iirFilter ), &stat );
  stop=clock();
  if( doTime )
    fprintf( stdout, "Elapsed time for IIRFilterREAL4Vector():  %.2f"
	     " s\n", (double)( stop - start )/CLOCKS_PER_SEC );
  if ( print ) {
    if( !( fp = fopen( OUTFILE3, "w" ) ) ) {
      ERROR( IIRFILTERTESTC_EFILE, IIRFILTERTESTC_MSGEFILE, 0 );
      return IIRFILTERTESTC_EFILE;
    }
    PrintVector( fp, input1 );
  }

  /* Filter the time series using IIRFilterREAL4VectorR(). */
  j = reps;
  start = clock();
  while ( j-- )
    SUB( LALIIRFilterREAL4VectorR( &stat, input2, iirFilter ), &stat );
  stop=clock();
  if ( doTime )
    fprintf( stdout, "Elapsed time for IIRFilterREAL4VectorR(): %.2f"
	     " s\n", (double)( stop - start )/CLOCKS_PER_SEC );
  if( print ) {
    if( !( fp = fopen( OUTFILE4, "w" ) ) ) {
      ERROR( IIRFILTERTESTC_EFILE, IIRFILTERTESTC_MSGEFILE, 0 );
      return IIRFILTERTESTC_EFILE;
    }
    PrintVector( fp, input2 );
  }

  /* Free memory and exit. */
  SUB( LALSDestroyVector( &stat, &input1 ), &stat );
  SUB( LALSDestroyVector( &stat, &input2 ), &stat );
  SUB( LALSDestroyVector( &stat, &output ), &stat );
  SUB( LALDestroyREAL4IIRFilter( &stat, &iirFilter ), &stat );
  LALCheckMemoryLeaks();
  INFO( IIRFILTERTESTC_MSGENORM );
  return IIRFILTERTESTC_ENORM;
}


static void
PrintVector(FILE *fp, REAL4Vector *vector)
{
  INT4 i=vector->length;
  REAL4 *data=vector->data;

  while ( i-- )
    fprintf( fp, "%10.3e\n", *(data++) );

  return;
}
