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
IIRFilterTest [-f filtertag] [-o] [-t] [-d debuglevel] [-n npts] [-r reps] [-w freq]
\end{verbatim}

\subsubsection*{Description}

This program generates a time series vector, and passes it through a
third-order Butterworth low-pass filter.  By default, running this
program with no arguments simply passes an impulse function to the
filter routines, producing no output.  All filter parameters are set
from \verb@#define@d constants.  The following option flags are
accepted:
\begin{itemize}
\item[\texttt{-f}] Specifies which filter(s) to be used:
\verb@filtertag@ is a token containing one or more character codes
from \verb@a@ to \verb@f@ and/or from \verb@A@ to \verb@D@, each
corresponding to a different filter:
\begin{center}
\begin{tabular}{rcl@{\qquad\qquad}rcl}
\verb@a@ &=& \verb@LALSIIRFilter()@             &
\verb@A@ &=& \verb@LALDIIRFilter()@             \\
\verb@b@ &=& \verb@LALIIRFilterREAL4()@         &
\verb@B@ &=& \verb@LALIIRFilterREAL8()@         \\
\verb@c@ &=& \verb@LALIIRFilterREAL4Vector()@   &
\verb@C@ &=& \verb@LALIIRFilterREAL8Vector()@   \\
\verb@d@ &=& \verb@LALIIRFilterREAL4VectorR()@  &
\verb@D@ &=& \verb@LALIIRFilterREAL8VectorR()@  \\
\verb@e@ &=& \verb@LALDIIRFilterREAL4Vector()@  &&& \\
\verb@f@ &=& \verb@LALDIIRFilterREAL4VectorR()@ &&&
\end{tabular}
\end{center}
If not specified, \verb@-f abcd@ is assumed.

\item[\texttt{-o}] Prints the input and output vectors to data files:
\verb@out.0@ stores the initial impulse, \verb@out.@$c$ the response
computed using the filter with character code $c$ (above).  If not
specified, the routines are exercised, but no output is written.

\item[\texttt{-t}] Causes \verb@IIRFilterTest@ to fill the time series
vector with Gaussian random deviates, and prints execution times for
the various filter subroutines to \verb@stdout@.  To generate useful
timing data, the default size of the time vector is increased (unless
explicitly set, below).

\item[\texttt{-d}] Changes the debug level from 0 to the specified
value \verb@debuglevel@.

\item[\texttt{-n}] Sets the size of the time vectors to \verb@npts@.
If not specified, 4096 points are used (4194304 if the \verb@-t@
option was also given).

\item[\texttt{-r}] Applies each filter to the data \verb@reps@ times
instead of just once.

\item[\texttt{-w}] Sets the characteristic frequency of the filter to
\verb@freq@ in the $w$-plane (described in \verb@ZPGFilter.h@).  If
not specified, \verb@-w 0.01@ is assumed (i.e.\ a characteristic
frequency of 2\% of Nyquist).
\end{itemize}

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

\newpage

\subsubsection*{Algorithm}

\begin{wrapfigure}{r}{0.6\textwidth}
\vspace{-7ex}
\begin{center}
\resizebox{0.55\textwidth}{!}{\includegraphics{tdfilter_iirfiltertest}}
\\ \parbox{0.55\textwidth}{\caption{\label{fig:iirfiltertest} Impulse
response functions in the frequency domain (computed as a windowed
FFT) for various filtering routines, run using the \texttt{-r~5}
option.}}
\end{center}
\vspace{-4ex}
\end{wrapfigure}
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

\subsubsection*{Sample output}

Running this program on a 1.3~GHz Intel machine with no optimization
produced the following typical timing information:
\begin{verbatim}
> IIRFilterTest -r 5 -t -f abcdefABCD
Filtering 4194304 points 5 times:
Elapsed time for LALSIIRFilter():             1.39 s
Elapsed time for LALDIIRFilter():             1.79 s
Elapsed time for LALIIRFilterREAL4():         2.86 s
Elapsed time for LALIIRFilterREAL8():         3.25 s
Elapsed time for LALIIRFilterREAL4Vector():   1.52 s
Elapsed time for LALIIRFilterREAL8Vector():   2.13 s
Elapsed time for LALIIRFilterREAL4VectorR():  1.33 s
Elapsed time for LALIIRFilterREAL8VectorR():  1.96 s
Elapsed time for LALDIIRFilterREAL4Vector():  1.12 s
Elapsed time for LALDIIRFilterREAL4VectorR(): 1.06 s
\end{verbatim}
From these results it is clear that the mixed-precision vector
filtering routines are the most efficient, outperforming even the
purely single-precision vector filtering routines by 20\%--30\%.  This
was unanticipated; by my count the main inner loop of the
single-precision routines contain $2M+2N+1$ dereferences and $2M+2N-3$
floating-point operations per vector element, where $M$ and $N$ are
the direct and recursive filter orders, whereas the mixed-precision
routines contain $2\max\{M,N\}+2M+2$ dereferences and $2M+2N-1$
floating-point operations per element.  However, most of the
dereferences in the mixed-precision routines are to short internal
arrays rather than to the larger data vector, which might cause some
speedup.

Running the same command with the \verb@-t@ flag replaced with
\verb@-o@ generates files containing the impulse response of the
filters.  The frequency-domain impulse response is shown in
Fig.~\ref{fig:iirfiltertest}.  This shows the steady improvement in
truncation error from single- to mixed- to double-precision filtering.


\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALSCreateVector()              LALDCreateVector()
LALSDestroyVector()             LALDDestroyVector()
LALCreateRandomParams()         LALDestroyRandomParams()
LALNormalDeviates()             LALPrintError()
LALCreateCOMPLEX8ZPGFilter()    LALCreateCOMPLEX16ZPGFilter()
LALDestroyCOMPLEX8ZPGFilter()   LALDestroyCOMPLEX16ZPGFilter()
LALCreateREAL4IIRFilter()       LALCreateREAL8IIRFilter()
LALDestroyREAL4IIRFilter()      LALDestroyREAL8IIRFilter()
LALWToZCOMPLEX8ZPGFilter()      LALWToZCOMPLEX16ZPGFilter()
LALSIIRFilter()                 LALDIIRFilter()
LALIIRFilterREAL4()             LALIIRFilterREAL8()
LALIIRFilterREAL4Vector()       LALIIRFilterREAL8Vector()
LALIIRFilterREAL4VectorR()      LALIIRFilterREAL8VectorR()
LALDIIRFilterREAL4Vector()      LALDIIRFilterREAL4VectorR()
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
#define TAG    "abcd"  /* Filter codes */

/* Mathematical constant. */
#define SQRT3_2 0.8660254037844386467637231707529361L

/* Usage format string. */
#define USAGE "Usage: %s [-d debuglevel] [-f filtertag] [-n npts]\n" \
"\t[-r reps] [-w freq] [-o] [-t]\n"

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

/* A script to assign the zpgFilter, which may be single or double
   precision. */
#define ASSIGNFILTER                                                 \
do {                                                                 \
  zpgFilter->poles->data[0].re = wc*SQRT3_2;                         \
  zpgFilter->poles->data[0].im = wc*0.5;                             \
  zpgFilter->poles->data[1].re = 0.0;                                \
  zpgFilter->poles->data[1].im = wc;                                 \
  zpgFilter->poles->data[2].re = -wc*SQRT3_2;                        \
  zpgFilter->poles->data[2].im = wc*0.5;                             \
  zpgFilter->gain.re = 0.0;                                          \
  zpgFilter->gain.im = wc*wc*wc;                                     \
} while (0)

/* A script to print a data vector, which may be single or double
   precision. */
#define PRINTDATA( outfile )                                         \
do {                                                                 \
  FILE *fp = fopen( (outfile), "w" );                                \
  if ( !fp ) {                                                       \
    ERROR( IIRFILTERTESTC_EFILE, IIRFILTERTESTC_MSGEFILE, 0 );       \
    return IIRFILTERTESTC_EFILE;                                     \
  }                                                                  \
  i = npts;                                                          \
  while ( i-- )                                                      \
    fprintf( fp, "%23.16e\n", *(data++) );                           \
} while (0)


/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

int
main(int argc, char **argv)
{
  static LALStatus stat; /* Status pointer for subroutines */
  BOOLEAN print = 0;     /* Whether output will be printed */
  BOOLEAN doTime = 0;    /* Whether filters will be timed */
  BOOLEAN points = 0;    /* Whether number of points is given */
  const CHAR *tag = TAG; /* List of filter character codes */
  INT4 npts = NPTS;      /* Number of points in time series */
  INT4 reps = REPS;      /* Number of repeated filterings */
  INT4 arg;              /* argument counter */
  INT4 i, j;             /* Index counters */
  REAL8 wc = WC;         /* Characteristic w-plane frequency */
  REAL4Vector *sInput = NULL;  /* REAL4 time series input vector */
  REAL8Vector *dInput = NULL;  /* REAL8 time series input vector */
  REAL4Vector *sOutput = NULL; /* REAL4 time series output vector */
  REAL8Vector *dOutput = NULL; /* REAL8 time series output vector */
  REAL4IIRFilter *sIIRFilter = NULL; /* The REAL4 IIR filter */
  REAL8IIRFilter *dIIRFilter = NULL; /* The REAL8 IIR filter */
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
    /* Parse filtertag option. */
    else if ( !strcmp( argv[arg], "-f" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        tag = argv[arg++];
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

  /* Create the time-domain filter(s). */
  if ( strchr( tag, 'a' ) || strchr( tag, 'b' ) ||
       strchr( tag, 'c' ) || strchr( tag, 'd' ) ) {
    COMPLEX8ZPGFilter *zpgFilter = NULL;
    SUB( LALCreateCOMPLEX8ZPGFilter( &stat, &zpgFilter, 0, 3 ),
	 &stat );
    ASSIGNFILTER;
    SUB( LALWToZCOMPLEX8ZPGFilter( &stat, zpgFilter ), &stat );
    SUB( LALCreateREAL4IIRFilter( &stat, &sIIRFilter, zpgFilter ),
	 &stat );
    SUB( LALDestroyCOMPLEX8ZPGFilter( &stat, &zpgFilter ), &stat );
  }
  if ( strchr( tag, 'A' ) || strchr( tag, 'B' ) ||
       strchr( tag, 'C' ) || strchr( tag, 'D' ) ||
       strchr( tag, 'e' ) || strchr( tag, 'f' ) ) {
    COMPLEX16ZPGFilter *zpgFilter = NULL;
    SUB( LALCreateCOMPLEX16ZPGFilter( &stat, &zpgFilter, 0, 3 ),
	 &stat );
    ASSIGNFILTER;
    SUB( LALWToZCOMPLEX16ZPGFilter( &stat, zpgFilter ), &stat );
    SUB( LALCreateREAL8IIRFilter( &stat, &dIIRFilter, zpgFilter ),
	 &stat );
    SUB( LALDestroyCOMPLEX16ZPGFilter( &stat, &zpgFilter ), &stat );
  }
  if ( !sIIRFilter && !dIIRFilter ) {
    ERROR( IIRFILTERTESTC_EBAD, IIRFILTERTESTC_MSGEBAD, "filtertag:" );
    LALPrintError( USAGE, *argv );
    return IIRFILTERTESTC_EBAD;
  }

  /* Create the input time series. */
  SUB( LALSCreateVector( &stat, &sInput, npts ), &stat );
  if ( doTime ) {
    RandomParams *params = NULL; /* Params for random generator */
    SUB( LALCreateRandomParams( &stat, &params, 0 ), &stat );
    INFO( "Begining generation of Gaussian random deviates.\n" );
    SUB( LALNormalDeviates( &stat, sInput, params ), &stat );
    INFO( "Finished generating random deviates.\n" );
    SUB( LALDestroyRandomParams( &stat, &params ), &stat );
  } else {
    memset( sInput->data, 0, npts*sizeof(REAL4) );
    sInput->data[npts/2] = 1.0;
  }
  if ( print ) {
    REAL4 *data = sInput->data;
    PRINTDATA( "out.0" );
  }

  /* Create other time series. */
  if ( strchr( tag, 'a' ) || strchr( tag, 'b' ) ||
       strchr( tag, 'c' ) || strchr( tag, 'd' ) ||
       strchr( tag, 'e' ) || strchr( tag, 'f' ) ) {
    SUB( LALSCreateVector( &stat, &sOutput, npts ), &stat );
  }
  if ( strchr( tag, 'A' ) || strchr( tag, 'B' ) ||
       strchr( tag, 'C' ) || strchr( tag, 'D' ) ) {
    SUB( LALDCreateVector( &stat, &dInput, npts ), &stat );
    SUB( LALDCreateVector( &stat, &dOutput, npts ), &stat );
    for ( i = 0; i < npts; i++ )
      dInput->data[i] = (REAL8)( sInput->data[i] );
  }

  /* Filter the time series. */
  if ( doTime ) {
    if ( reps == 1 )
      fprintf( stdout, "Filtering %i points:\n", npts );
    else
      fprintf( stdout, "Filtering %i points %i times:\n", npts,
	       reps );
  }

  /* Using LALSIIRFilter(): */
  if ( strchr( tag, 'a' ) ) {
    memcpy( sOutput->data, sInput->data, npts*sizeof(REAL4) );
    j = reps;
    start = clock();
    while ( j-- ) {
      REAL4 *data = sOutput->data;
      i = npts;
      while ( i-- ) {
	*data = LALSIIRFilter( *data, sIIRFilter );
	data++;
      }
    }
    stop = clock();
    if ( doTime )
      fprintf( stdout, "Elapsed time for LALSIIRFilter():             %.2f"
	       " s\n", (double)( stop - start )/CLOCKS_PER_SEC );
    if ( print ) {
      REAL4 *data = sOutput->data;
      PRINTDATA( "out.a" );
    }
  }

  /* Using LALDIIRFilter(): */
  if ( strchr( tag, 'A' ) ) {
    memcpy( dOutput->data, dInput->data, npts*sizeof(REAL8) );
    j = reps;
    start = clock();
    while ( j-- ) {
      REAL8 *data = dOutput->data;
      i = npts;
      while ( i-- ) {
	*data = LALDIIRFilter( *data, dIIRFilter );
	data++;
      }
    }
    stop = clock();
    if ( doTime )
      fprintf( stdout, "Elapsed time for LALDIIRFilter():             %.2f"
	       " s\n", (double)( stop - start )/CLOCKS_PER_SEC );
    if ( print ) {
      REAL8 *data = dOutput->data;
      PRINTDATA( "out.A" );
    }
  }

  /* Using LALIIRFilterREAL4(). */
  if ( strchr( tag, 'b' ) ) {
    memset( sIIRFilter->history->data, 0,
	    sIIRFilter->history->length*sizeof(REAL4) );
    memcpy( sOutput->data, sInput->data, npts*sizeof(REAL4) );
    j = reps;
    start = clock();
    while ( j-- ) {
      REAL4 *data = sOutput->data;
      i = npts;
      while ( i-- ) {
	SUB( LALIIRFilterREAL4( &stat, data, *data, sIIRFilter ),
	     &stat );
	data++;
      }
    }
    stop = clock();
    if ( doTime )
      fprintf( stdout, "Elapsed time for LALIIRFilterREAL4():         %.2f"
	       " s\n", (double)( stop - start )/CLOCKS_PER_SEC );
    if( print ) {
      REAL4 *data = sOutput->data;
      PRINTDATA( "out.b" );
    }
  }

  /* Using LALIIRFilterREAL8(). */
  if ( strchr( tag, 'B' ) ) {
    memset( dIIRFilter->history->data, 0,
	    dIIRFilter->history->length*sizeof(REAL8) );
    memcpy( dOutput->data, dInput->data, npts*sizeof(REAL8) );
    j = reps;
    start = clock();
    while ( j-- ) {
      REAL8 *data = dOutput->data;
      i = npts;
      while ( i-- ) {
	SUB( LALIIRFilterREAL8( &stat, data, *data, dIIRFilter ),
	     &stat );
	data++;
      }
    }
    stop = clock();
    if ( doTime )
      fprintf( stdout, "Elapsed time for LALIIRFilterREAL8():         %.2f"
	       " s\n", (double)( stop - start )/CLOCKS_PER_SEC );
    if( print ) {
      REAL8 *data = dOutput->data;
      PRINTDATA( "out.B" );
    }
  }

  /* Using LALIIRFilterREAL4Vector(). */
  if ( strchr( tag, 'c' ) ) {
    memset( sIIRFilter->history->data, 0,
	    sIIRFilter->history->length*sizeof(REAL4) );
    memcpy( sOutput->data, sInput->data, npts*sizeof(REAL4) );
    j = reps;
    start = clock();
    while ( j-- )
      SUB( LALIIRFilterREAL4Vector( &stat, sOutput, sIIRFilter ),
	   &stat );
    stop=clock();
    if( doTime )
      fprintf( stdout, "Elapsed time for LALIIRFilterREAL4Vector():   %.2f"
	       " s\n", (double)( stop - start )/CLOCKS_PER_SEC );
    if ( print ) {
      REAL4 *data = sOutput->data;
      PRINTDATA( "out.c" );
    }
  }

  /* Using LALIIRFilterREAL8Vector(). */
  if ( strchr( tag, 'C' ) ) {
    memset( dIIRFilter->history->data, 0,
	    dIIRFilter->history->length*sizeof(REAL8) );
    memcpy( dOutput->data, dInput->data, npts*sizeof(REAL8) );
    j = reps;
    start = clock();
    while ( j-- )
      SUB( LALIIRFilterREAL8Vector( &stat, dOutput, dIIRFilter ),
	   &stat );
    stop=clock();
    if( doTime )
      fprintf( stdout, "Elapsed time for LALIIRFilterREAL8Vector():   %.2f"
	       " s\n", (double)( stop - start )/CLOCKS_PER_SEC );
    if ( print ) {
      REAL8 *data = dOutput->data;
      PRINTDATA( "out.C" );
    }
  }

  /* Using LALIIRFilterREAL4VectorR(). */
  if ( strchr( tag, 'd' ) ) {
    memcpy( sOutput->data, sInput->data, npts*sizeof(REAL4) );
    j = reps;
    start = clock();
    while ( j-- )
      SUB( LALIIRFilterREAL4VectorR( &stat, sOutput, sIIRFilter ),
	   &stat );
    stop=clock();
    if ( doTime )
      fprintf( stdout, "Elapsed time for LALIIRFilterREAL4VectorR():  %.2f"
	       " s\n", (double)( stop - start )/CLOCKS_PER_SEC );
    if( print ) {
      REAL4 *data = sOutput->data;
      PRINTDATA( "out.d" );
    }
  }

  /* Using LALIIRFilterREAL8VectorR(). */
  if ( strchr( tag, 'D' ) ) {
    memcpy( dOutput->data, dInput->data, npts*sizeof(REAL8) );
    j = reps;
    start = clock();
    while ( j-- )
      SUB( LALIIRFilterREAL8VectorR( &stat, dOutput, dIIRFilter ),
	   &stat );
    stop=clock();
    if ( doTime )
      fprintf( stdout, "Elapsed time for LALIIRFilterREAL8VectorR():  %.2f"
	       " s\n", (double)( stop - start )/CLOCKS_PER_SEC );
    if( print ) {
      REAL8 *data = dOutput->data;
      PRINTDATA( "out.D" );
    }
  }

  /* Using LALDIIRFilterREAL4Vector(). */
  if ( strchr( tag, 'e' ) ) {
    memset( dIIRFilter->history->data, 0,
	    dIIRFilter->history->length*sizeof(REAL8) );
    memcpy( sOutput->data, sInput->data, npts*sizeof(REAL4) );
    j = reps;
    start = clock();
    while ( j-- )
      SUB( LALDIIRFilterREAL4Vector( &stat, sOutput, dIIRFilter ),
	   &stat );
    stop=clock();
    if( doTime )
      fprintf( stdout, "Elapsed time for LALDIIRFilterREAL4Vector():  %.2f"
	       " s\n", (double)( stop - start )/CLOCKS_PER_SEC );
    if ( print ) {
      REAL4 *data = sOutput->data;
      PRINTDATA( "out.e" );
    }
  }

  /* Using LALDIIRFilterREAL4VectorR(). */
  if ( strchr( tag, 'f' ) ) {
    memcpy( sOutput->data, sInput->data, npts*sizeof(REAL4) );
    j = reps;
    start = clock();
    while ( j-- )
      SUB( LALDIIRFilterREAL4VectorR( &stat, sOutput, dIIRFilter ),
	   &stat );
    stop=clock();
    if( doTime )
      fprintf( stdout, "Elapsed time for LALDIIRFilterREAL4VectorR(): %.2f"
	       " s\n", (double)( stop - start )/CLOCKS_PER_SEC );
    if ( print ) {
      REAL4 *data = sOutput->data;
      PRINTDATA( "out.f" );
    }
  }

  /* Free memory and exit. */
  SUB( LALSDestroyVector( &stat, &sInput ), &stat );
  if ( sOutput )
    SUB( LALSDestroyVector( &stat, &sOutput ), &stat );
  if ( dInput )
    SUB( LALDDestroyVector( &stat, &dInput ), &stat );
  if ( dOutput )
    SUB( LALDDestroyVector( &stat, &dOutput ), &stat );
  if ( sIIRFilter )
    SUB( LALDestroyREAL4IIRFilter( &stat, &sIIRFilter ), &stat );
  if ( dIIRFilter )
    SUB( LALDestroyREAL8IIRFilter( &stat, &dIIRFilter ), &stat );
  LALCheckMemoryLeaks();
  INFO( IIRFILTERTESTC_MSGENORM );
  return IIRFILTERTESTC_ENORM;
}
