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
BandPassTest [-d debuglevel] [-i infile | -n npts dt offset] [-o outfile]
             [-f f1 f2 a1 a2 order]
\end{verbatim}

\subsubsection*{Description}

This program applies a Butterworth time-domain low-pass or high-pass
filter to a time series, using the routine
\verb@LALDButterworthREAL4TimeSeries()@.  The following option flags
are accepted:
\begin{itemize}
\item[\texttt{-d}] Changes the default debug level from 0 to
\verb@debuglevel@.
\item[\texttt{-i}] Reads the input time series from \verb@infile@
using the routine \verb@LALSReadTSeries()@; see \verb@StreamInput.h@
for a description of the file format.
\item[\texttt{-n}] Generates an input time series of length
\verb@npts@ and sampling interval \verb@dt@, containing just an
impulse at sample index \verb@offset@.  If the \verb@-i@ option is
also given, it overrides this option.  If neither are given,
\verb@-n 4096 1.0 1024@ is assumed.
\item[\texttt{-o}] Writes the output time series to \verb@outfile@,
using the routine \verb@LALSWriteTSeries()@; see \verb@StreamOutput.h@
for a description of the file format.  If not specified, the routines
are exercised, but no output is written.
\item[\texttt{-f}] Sets the filter to have attenuation \verb@a1@ and
\verb@a2@ at frequencies \verb@f1@ and \verb@f2@, with a maximum
filter order of \verb@order@; see \verb@ButterworthTimeSeries.c@ for a
description of how these values are interpreted.  If not specified,
\verb@-f 0.01 0.015 0.9 0.1 20@ is assumed.
\end{itemize}

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
#define BANDPASSTESTC_MSGEFILE "Could not open file"
/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALPrintError()                         LALSnprintf()
LALSCreateVector()                      LALSDestroyVector()
LALSReadTSeries()                       LALSWriteTSeries()
LALDButterworthREAL4TimeSeries()        LALCheckMemoryLeaks()
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
#include <lal/StreamInput.h>
#include <lal/StreamOutput.h>

NRCSID(BANDPASSTESTC,"$Id$");

/* Default parameters. */
INT4 lalDebugLevel=0;
#define NPTS 4096   /* Length of time series. */
#define DT 1.0      /* Sampling interval of time series. */
#define OFFSET 1024 /* Offset of the impulse from the start. */
#define F1 0.01     /* Lower frequency of transition band. */
#define F2 0.015    /* Upper frequency of transition band. */
#define A1 0.9      /* Desired attenuation at F1. */
#define A2 0.1      /* Desired attenuation at F2. */
#define ORDER 20    /* Maximum filter order. */

/* Usage format string. */
#define USAGE "Usage: %s [-d debuglevel] [-i infile | -n npts dt offset]\n" \
"\t[-o outfile] [-f f1 f2 a1 a2 order]\n"

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
  static LALStatus stat; /* LALStatus pointer */
  CHAR *infile = NULL;   /* The input filename */
  CHAR *outfile = NULL;  /* The output filename */
  INT4 arg;              /* Argument counter */
  UINT4 npts = NPTS;     /* Number of points in time series */
  UINT4 offset = OFFSET; /* Position of delta function */
  REAL8 dt = DT;         /* Sampling interval. */
  static REAL4TimeSeries series;    /* Time series */
  static PassBandParamStruc params; /* Filter parameters */

  XLALSetErrorHandler( XLALAbortErrorHandler );

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
    /* Parse input file option. */
    else if ( !strcmp( argv[arg], "-i" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        infile = argv[arg++];
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
        outfile = argv[arg++];
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
      if ( argc > arg + 3 ) {
        arg++;
	npts=atoi(argv[arg++]);
	dt=atof(argv[arg++]);
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

  /* Check input values. */
  if ( !infile ) {
    if ( offset >= npts ) {
      ERROR( BANDPASSTESTC_EBAD, BANDPASSTESTC_MSGEBAD, 0 );
      LALPrintError( "\toffset=%i must be less than npts=%i\n", offset,
		     npts );
      return BANDPASSTESTC_EBAD;
    }
  }

  /* Create the time series. */
  if ( infile ) {
    FILE *fp = fopen( infile, "r" );
    if ( !fp ) {
      ERROR( BANDPASSTESTC_EFILE, BANDPASSTESTC_MSGEFILE, infile );
      return BANDPASSTESTC_EFILE;
    }
    SUB( LALSReadTSeries( &stat, &series, fp ), &stat );
    fclose( fp );
  } else {
    LALSnprintf( series.name, LALNameLength, "%s", "Impulse" );
    series.deltaT = dt;
    SUB( LALSCreateVector( &stat, &(series.data), npts ), &stat );
    memset( series.data->data, 0, npts*sizeof(REAL4) );
    series.data->data[offset] = 1.0;
  }

  /* Filter the time series. */
  SUB( LALDButterworthREAL4TimeSeries( &stat, &series, &params ),
       &stat );

  /* Print the output, if the -o option was given. */
  if ( outfile ) {
    FILE *fp = fopen( outfile, "w" );
    if ( !fp ){
      ERROR( BANDPASSTESTC_EFILE, BANDPASSTESTC_MSGEFILE, outfile );
      return BANDPASSTESTC_EFILE;
    }
    SUB( LALSWriteTSeries( &stat, fp, &series ), &stat );
    fclose( fp );
  }

  /* Free memory and exit. */
  SUB( LALSDestroyVector( &stat, &(series.data) ), &stat );
  LALCheckMemoryLeaks();
  INFO( BANDPASSTESTC_MSGENORM );
  return BANDPASSTESTC_ENORM;
}
