/********************************* <lalVerbatim file="BandPassTestCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{BandPassTest.c}}
\label{s:BandPassTest.c}

Tests time-domain high- and low-pass filters.

\subsubsection*{Usage}
\begin{verbatim}
BandPassTest [-d debuglevel] [-o outfile] [-f f1 f2 a1 a2 order] [-n npts offset]
\end{verbatim}

\subsubsection*{Description}

This program generates a time series with an impulse in it, and passes
it through a time-domain low-pass or high-pass filter.  By default,
running this program with no arguments simply tests the subroutines,
producing no output.  All filter parameters are set from
\verb@#define@d constants.

The \verb@-d@ option sets the debug level to the specified value
\verb@debuglevel@.  The \verb@-o@ flag tells the program to print the
impulse response to the specified data file \verb@outfile@.  The
\verb@-f@ option sets the filter to have power attenuations \verb@a1@
and \verb@a2@ at the frequencies \verb@f1@ and \verb@f2@ (in units of
the sampling frequency).  The \verb@-n@ option sets the length of the
time series to \verb@npts@ and places the impulse a number of samples
\verb@offset@ into it.

\subsubsection*{Exit codes}
****************************************** </lalLaTeX><lalErrTable> */
#define BANDPASSTESTC_ENORM 0
#define BANDPASSTESTC_ESUB  1
#define BANDPASSTESTC_EARG  2
#define BANDPASSTESTC_EBAD  3
#define BANDPASSTESTC_EFILE 4

#define BANDPASSTESTC_MSGENORM "Normal exit"
#define BANDPASSTESTC_MSGESUB  "Subroutine failed"
#define BANDPASSTESTC_MSGEARG  "Error parsing arguments"
#define BANDPASSTESTC_MSGEBAD  "Bad argument values"
#define BANDPASSTESTC_MSGEFILE "Could not create output file"
/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALPrintError()
LALSCreateVector()              LALSDestroyVector()
LALButterworthREAL4TimeSeries()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{BandPassTestCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <lal/AVFactories.h>
#include <lal/BandPassTimeSeries.h>

NRCSID(BANDPASSTESTC,"$Id$");

/* Default parameters. */
INT4 lalDebugLevel=0;
#define F1 0.01     /* Lower frequency of transition band. */
#define F2 0.015    /* Upper frequency of transition band. */
#define A1 0.9      /* Desired attenuation at F1. */
#define A2 0.1      /* Desired attenuation at F2. */
#define ORDER 20    /* Maximum filter order. */
#define NPTS 4096   /* Length of time series. */
#define OFFSET 1024 /* Offset of the impulse from the start. */

/* Usage format string. */
#define USAGE "Usage: %s [-d debuglevel] [-o outfile] [-f f1 f2 a1 a2 order] [-n npts offset]\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
		   "        %s %s\n", (code), *argv, __FILE__,       \
		   __LINE__, BANDPASSTESTC, statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
		   "        %s\n", *argv, __FILE__, __LINE__,        \
		   BANDPASSTESTC, (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( BANDPASSTESTC_ESUB, BANDPASSTESTC_MSGESUB,                \
           "Function call \"" #func "\" failed:" );                  \
    return BANDPASSTESTC_ESUB;                                       \
  }                                                                  \
} while (0)

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

int
main(int argc, char **argv)
{
  static LALStatus stat;            /* LALStatus pointer */
  CHAR *fname = NULL;               /* The output filename */
  INT4 arg;                         /* Argument counter */
  UINT4 i;                          /* Index counter */
  UINT4 npts = NPTS;                /* Num. of points in time series */
  UINT4 offset = OFFSET;            /* Position of delta function */
  static REAL4TimeSeries series;    /* Time series */
  REAL4 *data;                      /* Time series data */
  static PassBandParamStruc params; /* Filter parameters */
  FILE *fp=NULL;                    /* Output file */

  /* Set up the default filter parameters. */
  params.f1 = F1;
  params.f2 = F2;
  params.a1 = A1;
  params.a2 = A2;
  params.nMax = ORDER;

  /* Parse argument list.  i stores the current position. */
  arg = 1;
  while ( arg < argc ) {
    /* Parse debuglevel option. */
    if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        lalDebugLevel = atoi( argv[arg++] );
      } else {
	ERROR( BANDPASSTESTC_EARG, BANDPASSTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return BANDPASSTESTC_EARG;
      }
    }
    /* Parse output file option. */
    else if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        fname = argv[arg++];
      } else {
	ERROR( BANDPASSTESTC_EARG, BANDPASSTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return BANDPASSTESTC_EARG;
      }
    }
    /* Parse filter options. */
    else if ( !strcmp( argv[arg], "-f" ) ) {
      if ( argc > arg + 5 ) {
        arg++;
	params.f1=atof(argv[arg++]);
	params.f2=atof(argv[arg++]);
	params.a1=atof(argv[arg++]);
	params.a2=atof(argv[arg++]);
	params.nMax=atoi(argv[arg++]);
      } else {
	ERROR( BANDPASSTESTC_EARG, BANDPASSTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return BANDPASSTESTC_EARG;
      }
    }
    /* Parse time series options. */
    else if ( !strcmp( argv[arg], "-n" ) ) {
      if ( argc > arg + 2 ) {
        arg++;
	npts=atoi(argv[arg++]);
	offset=atoi(argv[arg++]);
      } else {
	ERROR( BANDPASSTESTC_EARG, BANDPASSTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return BANDPASSTESTC_EARG;
      }
    }
    /* Unrecognized option. */
    else {
      ERROR( BANDPASSTESTC_EARG, BANDPASSTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return BANDPASSTESTC_EARG;
    }
  } /* End of argument parsing loop. */

  if ( offset >= npts ) {
    ERROR( BANDPASSTESTC_EBAD, BANDPASSTESTC_MSGEBAD, 0 );
    LALPrintError( "\toffset=%i must be less than npts=%i\n", offset,
		   npts );
    return BANDPASSTESTC_EBAD;
  }

  /* Create the time series. */
  LALSnprintf( series.name, LALNameLength, "%s", "Impulse" );
  series.deltaT = 1.0;
  SUB( LALSCreateVector( &stat, &(series.data), npts ), &stat );
  memset( series.data->data, 0, npts*sizeof(REAL4) );
  series.data->data[offset] = 1.0;

  /* Filter the time series. */
  SUB( LALButterworthREAL4TimeSeries( &stat, &series, &params ),
       &stat );

  /* Print the output, if the -o option was given. */
  if ( fname ) {
    fp = fopen( fname, "w" );
    if ( !fp ){
      ERROR( BANDPASSTESTC_EFILE, BANDPASSTESTC_MSGEFILE, 0 );
      return BANDPASSTESTC_EFILE;
    }
    for ( data = series.data->data, i=0; i<npts; data++, i++ )
      fprintf( fp, "%8.3e\n", *data );
    fclose( fp );
  }

  /* Free memory and exit. */
  SUB( LALSDestroyVector( &stat, &(series.data) ), &stat );
  LALCheckMemoryLeaks();
  INFO( BANDPASSTESTC_MSGENORM );
  return BANDPASSTESTC_ENORM;
}
