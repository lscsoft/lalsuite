/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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

/** \cond DONT_DOXYGEN */
/** \endcond */
#include <complex.h>
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


/**
 * \author Creighton, T. D.
 * \file
 * \ingroup IIRFilter_h
 *
 * \brief Tests the routines in \ref IIRFilter_h.
 *
 * ### Usage ###
 *
 * \code
 * IIRFilterTest [-f filtertag] [-o] [-t] [-d debuglevel] [-n npts] [-r reps] [-w freq]
 * \endcode
 *
 * ### Description ###
 *
 * This program generates a time series vector, and passes it through a
 * third-order Butterworth low-pass filter.  By default, running this
 * program with no arguments simply passes an impulse function to the
 * filter routines, producing no output.  All filter parameters are set
 * from <tt>\#define</tt>d constants.  The following option flags are
 * accepted:
 * <ul>
 * <li>[<tt>-f</tt>] Specifies which filter(s) to be used:
 * \c filtertag is a token containing one or more character codes
 * from \c a to \c f and/or from \c A to \c D, each
 * corresponding to a different filter:
 *
 * <table>
 * <tr><td>\c a</td><td>=</td><td><tt>LALSIIRFilter()</tt></td><td>\c A</td><td>=</td><td><tt>LALDIIRFilter()</tt></td></tr>
 * <tr><td>\c b</td><td>=</td><td><tt>LALIIRFilterREAL4()</tt></td><td>\c B</td><td>=</td><td><tt>LALIIRFilterREAL8()</tt></td></tr>
 * <tr><td>\c c</td><td>=</td><td><tt>LALIIRFilterREAL4Vector()</tt></td><td>\c C</td><td>=</td><td><tt>LALIIRFilterREAL8Vector()</tt></td></tr>
 * <tr><td>\c d</td><td>=</td><td><tt>LALIIRFilterREAL4VectorR()</tt></td><td>\c D</td><td>=</td><td><tt>LALIIRFilterREAL8VectorR()</tt></td></tr>
 * <tr><td>\c e</td><td>=</td><td><tt>LALDIIRFilterREAL4Vector()</tt></td><td></td><td></td><td></td></tr>
 * <tr><td>\c f</td><td>=</td><td><tt>LALDIIRFilterREAL4VectorR()</tt></td><td></td><td></td><td></td></tr>
 * </table>
 *
 * If not specified, <tt>-f abcd</tt> is assumed.</li>
 *
 * <li>[<tt>-o</tt>] Prints the input and output vectors to data files:
 * <tt>out.0</tt> stores the initial impulse, <tt>out.</tt>\f$c\f$ the response
 * computed using the filter with character code \f$c\f$ (above).  If not
 * specified, the routines are exercised, but no output is written.</li>
 *
 * <li>[<tt>-t</tt>] Causes \c IIRFilterTest to fill the time series
 * vector with Gaussian random deviates, and prints execution times for
 * the various filter subroutines to \c stdout.  To generate useful
 * timing data, the default size of the time vector is increased (unless
 * explicitly set, below).</li>
 *
 * <li>[<tt>-d</tt>] Changes the debug level from 0 to the specified
 * value \c debuglevel.</li>
 *
 * <li>[<tt>-n</tt>] Sets the size of the time vectors to \c npts.
 * If not specified, 4096 points are used (4194304 if the <tt>-t</tt>
 * option was also given).</li>
 *
 * <li>[<tt>-r</tt>] Applies each filter to the data \c reps times
 * instead of just once.</li>
 *
 * <li>[<tt>-w</tt>] Sets the characteristic frequency of the filter to
 * \c freq in the \f$w\f$-plane (described in \ref ZPGFilter.h).  If
 * not specified, <tt>-w 0.01</tt> is assumed (i.e.\ a characteristic
 * frequency of 2\% of Nyquist).</li>
 * </ul>
 *
 * ### Algorithm ###
 *
 * \image html tdfilter_iirfiltertest.png "Fig.[fig_iirfiltertest]: Impulse response functions in the frequency domain (computed as a windowed FFT) for various filtering routines, run using the <tt>-r 5</tt> option."
 * \image latex tdfilter_iirfiltertest.pdf "Impulse response functions in the frequency domain (computed as a windowed FFT) for various filtering routines, run using the <tt>-r 5</tt> option." width=0.55\textwidth
 *
 * A third-order Butterworth low-pass filter is defined by the following
 * power response function:
 * \f[
 * |T(w)|^2 = \frac{1}{1+(w/w_c)^6}\;,
 * \f]
 * where the frequency parameter \f$w=\tan(\pi f\Delta t)\f$ maps the Nyquist
 * interval onto the entire real axis, and \f$w_c\f$ is the (transformed)
 * characteristic frequency set by the <tt>-w</tt> option above.  A stable
 * time-domain filter has poles with \f$|z|<1\f$ in the \f$z\f$-plane
 * representation, which means that the \f$w\f$-plane poles should have
 * \f$\Im(w)>0\f$.  We construct a transfer function using only the
 * positive-imaginary poles of the power response function:
 * \f[
 * T(w) = \frac{iw_c^3}{\prod_{k=0}^2
 * \left[w-w_c e^{(2k+1)i\pi/6}\right]}\;,
 * \f]
 * where we have chosen a phase coefficient \f$i\f$ in the numerator in order
 * to get a purely real DC response.  This ensures that the \f$z\f$-plane
 * representation of the filter will have a real gain, resulting in a
 * physically-realizable IIR filter.
 *
 * The poles and gain of the transfer function \f$T(w)\f$ are simply read off
 * of the equation above, and are stored in a \c COMPLEX8ZPGFilter.
 * This is transformed from the \f$w\f$-plane to the \f$z\f$-plane representation
 * using <tt>LALWToZCOMPLEX8ZPGFilter()</tt>, and then used to create an
 * IIR filter with <tt>LALCreateREAL4IIRFilter()</tt>.  This in turn is
 * used by the routines <tt>LALSIIRFilter()</tt>,
 * <tt>LALIIRFilterREAL4()</tt>, <tt>LALIIRFilterREAL4Vector()</tt>, and
 * <tt>LALIIRFilterREAL4VectorR()</tt> to filter a data vector containing
 * either a unit impulse or white Gaussian noise (for more useful timing
 * information).
 *
 * ### Sample output ###
 *
 * Running this program on a 1.3 GHz Intel machine with no optimization
 * produced the following typical timing information:
 * \code
 * > IIRFilterTest -r 5 -t -f abcdefABCD
 * Filtering 4194304 points 5 times:
 * Elapsed time for LALSIIRFilter():             1.39 s
 * Elapsed time for LALDIIRFilter():             1.79 s
 * Elapsed time for LALIIRFilterREAL4():         2.86 s
 * Elapsed time for LALIIRFilterREAL8():         3.25 s
 * Elapsed time for LALIIRFilterREAL4Vector():   1.52 s
 * Elapsed time for LALIIRFilterREAL8Vector():   2.13 s
 * Elapsed time for LALIIRFilterREAL4VectorR():  1.33 s
 * Elapsed time for LALIIRFilterREAL8VectorR():  1.96 s
 * Elapsed time for LALDIIRFilterREAL4Vector():  1.12 s
 * Elapsed time for LALDIIRFilterREAL4VectorR(): 1.06 s
 * \endcode
 * From these results it is clear that the mixed-precision vector
 * filtering routines are the most efficient, outperforming even the
 * purely single-precision vector filtering routines by 20\%--30\%.  This
 * was unanticipated; by my count the main inner loop of the
 * single-precision routines contain \f$2M+2N+1\f$ dereferences and \f$2M+2N-3\f$
 * floating-point operations per vector element, where \f$M\f$ and \f$N\f$ are
 * the direct and recursive filter orders, whereas the mixed-precision
 * routines contain \f$2\max\{M,N\}+2M+2\f$ dereferences and \f$2M+2N-1\f$
 * floating-point operations per element.  However, most of the
 * dereferences in the mixed-precision routines are to short internal
 * arrays rather than to the larger data vector, which might cause some
 * speedup.
 *
 * Running the same command with the <tt>-t</tt> flag replaced with
 * <tt>-o</tt> generates files containing the impulse response of the
 * filters.  The frequency-domain impulse response is shown in
 * Fig.\figref{fig_iirfiltertest}.  This shows the steady improvement in
 * truncation error from single- to mixed- to double-precision filtering.
 *
 */

/** \name Error Codes */
/*@{*/
#define IIRFILTERTESTC_ENORM 0	/**< Normal exit */
#define IIRFILTERTESTC_ESUB  1	/**< Subroutine failed */
#define IIRFILTERTESTC_EARG  2	/**< Error parsing arguments */
#define IIRFILTERTESTC_EBAD  3	/**< Bad argument value */
#define IIRFILTERTESTC_EFILE 4	/**< Could not create output file */
/*@}*/

/** \cond DONT_DOXYGEN */
#define IIRFILTERTESTC_MSGENORM "Normal exit"
#define IIRFILTERTESTC_MSGESUB  "Subroutine failed"
#define IIRFILTERTESTC_MSGEARG  "Error parsing arguments"
#define IIRFILTERTESTC_MSGEBAD  "Bad argument value"
#define IIRFILTERTESTC_MSGEFILE "Could not create output file"

/* Default parameters. */
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
		   __LINE__, "$Id$", statement ? statement : \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
		   "        %s\n", *argv, __FILE__, __LINE__,        \
		   "$Id$", (statement) );                    \
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
  zpgFilter->poles->data[0] = wc*SQRT3_2;                         \
  zpgFilter->poles->data[0] += I*wc*0.5;                             \
  zpgFilter->poles->data[1] = 0.0;                                \
  zpgFilter->poles->data[1] += I*wc;                                 \
  zpgFilter->poles->data[2] = -wc*SQRT3_2;                        \
  zpgFilter->poles->data[2] += I*wc*0.5;                             \
  zpgFilter->gain = 0.0;                                          \
  zpgFilter->gain = I*wc*wc*wc;                                     \
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
/** \endcond */
