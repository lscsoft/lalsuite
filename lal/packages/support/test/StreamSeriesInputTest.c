/************************ <lalVerbatim file="StreamSeriesInputTestCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{StreamSeriesInputTest.c}}
\label{ss:StreamSeriesInputTest.c}

Reads a time or frequency series from a file, and writes it to another
file.

\subsubsection*{Usage}
\begin{verbatim}
StreamSeriesInputTest [-o outfile] [-i infile stype dtype] [-d debuglevel]
\end{verbatim}

\subsubsection*{Description}

This test program parses data from an input file or from \verb@stdin@,
using the routines in \verb@StreamSeriesInput.c@, and possibly
generating output using the routines in \verb@StreamSeriesOutput.c@.
The following option flags are accepted:
\begin{itemize}
\item[\texttt{-o}] Writes the output to \verb@outfile@.  If
\verb@outfile@ is given as \verb@stdout@, the data is written to
standard output (\emph{not} to a file named \verb@stdout@).  If the
\verb@-o@ flag is not given, the input routines are exercised, but no
output is written.
\item[\texttt{-i}] Specifies the input file name \verb@infile@, series
type \verb@stype@, and base datatype \verb@dtype@.  Series type is a
single character: either \verb@t@ (time series), \verb@v@ (time vector
series), \verb@a@ (time array series), or \verb@f@ (frequency series).
Base datatype may be \verb@i2@ (\verb@INT2@), \verb@i4@ (\verb@INT4@),
\verb@i8@ (\verb@INT8@), \verb@u2@ (\verb@UINT2@), \verb@u4@
(\verb@UINT4@), \verb@u8@ (\verb@UINT8@), \verb@s@ (\verb@REAL4@),
\verb@d@ (\verb@REAL8@), \verb@c@ (\verb@COMPLEX8@), or \verb@z@
(\verb@COMPLEX16@).  If the \verb@-i@ flag is not given,
\verb@-i StreamSeriesInput.dat f s@ is assumed (this file is provided
with the distribution so that running the code with no arguments, \'a
la \verb@make check@, will perform a nontrivial test of the
algorithm).
\item[\texttt{-d}] Sets the debug level to \verb@debuglevel@; if
absent, \verb@-d 0@ is assumed.
\end{itemize}

See the documentation in \verb@StreamSeriesInput.c@ and
\verb@StreamSeriesOutput.c@ for discussion of the input and output
data file formats.

\subsubsection*{Exit codes}
****************************************** </lalLaTeX><lalErrTable> */
#define STREAMSERIESINPUTTESTC_ENORM 0
#define STREAMSERIESINPUTTESTC_ESUB  1
#define STREAMSERIESINPUTTESTC_EARG  2
#define STREAMSERIESINPUTTESTC_EFILE 3

#define STREAMSERIESINPUTTESTC_MSGENORM "Normal exit"
#define STREAMSERIESINPUTTESTC_MSGESUB  "Subroutine failed"
#define STREAMSERIESINPUTTESTC_MSGEARG  "Error parsing arguments"
#define STREAMSERIESINPUTTESTC_MSGEFILE "Could not open file"
/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel                           LALPrintError()
LALOpenDataFile()                       LALCheckMemoryLeaks()
LAL<typecode>ReadTSeries()              LAL<typecode>WriteTSeries()
LAL<typecode>ReadTVectorSeries()        LAL<typecode>WriteTVectorSeries()
LAL<typecode>ReadTArraySeries()         LAL<typecode>WriteTArraySeries()
LAL<typecode>ReadFSeries()              LAL<typecode>WriteFSeries()
LAL<typecode>DestroyVector()
LAL<typecode>DestroyVectorSequence()
LAL<typecode>DestroyArraySequence()
\end{verbatim}
where \verb@<typecode>@ is any of \verb@I2@, \verb@I4@, \verb@I8@,
\verb@U2@, \verb@U4@, \verb@U8@, \verb@S@, \verb@D@, \verb@C@,
\verb@Z@.

\subsubsection*{Notes}

\vfill{\footnotesize\input{StreamSeriesInputTestCV}}

******************************************************* </lalLaTeX> */

#include <stdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/StreamInput.h>
#include <lal/StreamOutput.h>

NRCSID( STREAMSERIESINPUTTESTC, "$Id$" );

/* Default parameter settings. */
int lalDebugLevel = 0;
#define INFILE "StreamSeriesInput.data"

/* Usage format string. */
#define USAGE "Usage: %s [-o outfile] [-d debuglevel]\n"        \
"\t[-i infile {t | v | a | f} {i2 | i4 | i8 | u2 | u4 | u8 | s | d | c | z}]\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, STREAMSERIESINPUTTESTC,                   \
		 statement ? statement : "", (msg) );                \
}                                                                    \
else (void)(0)

#define INFO( statement )                                            \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 STREAMSERIESINPUTTESTC, (statement) );              \
}                                                                    \
else (void)(0)

#define SUB( func, statusptr )                                       \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( STREAMSERIESINPUTTESTC_ESUB, STREAMSERIESINPUTTESTC_MSGESUB,\
         "Function call \"" #func "\" failed:" );                    \
  return STREAMSERIESINPUTTESTC_ESUB;                                \
}                                                                    \
else (void)(0)

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

int
main(int argc, char **argv)
{
  static LALStatus stat;       /* top-level status structure */
  INT4 arg;                    /* index over command-line options */
  CHAR *outfile = NULL;        /* name of output file */
  const CHAR *infile = INFILE; /* name of input file */
  const CHAR *dtype = "s";     /* data type tag */
  const CHAR *stype = "f";     /* series type tag */
  FILE *fpIn = NULL;           /* input file pointer */
  FILE *fpOut = NULL;          /* output file pointer */

  /* Parse argument list.  arg stores the current position. */
  arg = 1;
  while ( arg < argc ) {
    /* Parse output file option. */
    if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	outfile = argv[arg++];
      } else {
	ERROR( STREAMSERIESINPUTTESTC_EARG,
	       STREAMSERIESINPUTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return STREAMSERIESINPUTTESTC_EARG;
      }
    }
    /* Parse input file option. */
    else if ( !strcmp( argv[arg], "-i" ) ) {
      if ( argc > arg + 3 ) {
	arg++;
	infile = argv[arg++];
	stype = argv[arg++];
	dtype = argv[arg++];
      } else {
	ERROR( STREAMSERIESINPUTTESTC_EARG,
	       STREAMSERIESINPUTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return STREAMSERIESINPUTTESTC_EARG;
      }
    }
    /* Parse debug level option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	lalDebugLevel = atoi( argv[arg++] );
      } else {
	ERROR( STREAMSERIESINPUTTESTC_EARG,
	       STREAMSERIESINPUTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return STREAMSERIESINPUTTESTC_EARG;
      }
    }
    /* Check for unrecognized options. */
    else {
      ERROR( STREAMSERIESINPUTTESTC_EARG,
	     STREAMSERIESINPUTTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return STREAMSERIESINPUTTESTC_EARG;
    }
  } /* End of argument parsing loop. */

  /* Open input and output files. */
  if ( strcmp( infile, "stdin" ) ) {
    if ( !( fpIn = LALOpenDataFile( infile ) ) ) {
      ERROR( STREAMSERIESINPUTTESTC_EFILE, "- "
	     STREAMSERIESINPUTTESTC_MSGEFILE, infile );
      return STREAMSERIESINPUTTESTC_EFILE;
    }
  } else
    fpIn = stdin;
  if ( outfile ) {
    if ( strcmp( outfile, "stdout" ) ) {
      if ( !( fpOut = fopen( outfile, "w" ) ) ) {
	ERROR( STREAMSERIESINPUTTESTC_EFILE, "- "
	       STREAMSERIESINPUTTESTC_MSGEFILE, outfile );
	return STREAMSERIESINPUTTESTC_EFILE;
      }
    } else
      fpOut = stdout;
  }

  /* Read (and possibly write) the data as a time series. */
  if ( !strcmp( stype, "t" ) ) {
    if ( !strcmp( dtype, "i2" ) ) {
      static INT2TimeSeries series;
      SUB( LALI2ReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI2WriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI2DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "i4" ) ) {
      static INT4TimeSeries series;
      SUB( LALI4ReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI4WriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI4DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "i8" ) ) {
      static INT8TimeSeries series;
      SUB( LALI8ReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI8WriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI8DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u2" ) ) {
      static UINT2TimeSeries series;
      SUB( LALU2ReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU2WriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU2DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u4" ) ) {
      static UINT4TimeSeries series;
      SUB( LALU4ReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU4WriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU4DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u8" ) ) {
      static UINT8TimeSeries series;
      SUB( LALU8ReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU8WriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU8DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "s" ) ) {
      static REAL4TimeSeries series;
      SUB( LALSReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALSWriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALSDestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "d" ) ) {
      static REAL8TimeSeries series;
      SUB( LALDReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALDWriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALDDestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "c" ) ) {
      static COMPLEX8TimeSeries series;
      SUB( LALCReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALCWriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALCDestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "z" ) ) {
      static COMPLEX16TimeSeries series;
      SUB( LALZReadTSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALZWriteTSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALZDestroyVector( &stat, &(series.data) ), &stat );
    } else {
      ERROR( STREAMSERIESINPUTTESTC_EARG,
	     STREAMSERIESINPUTTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return STREAMSERIESINPUTTESTC_EARG;
    }
  }

  /* Read (and possibly write) the data as a time vector series. */
  else if ( !strcmp( stype, "v" ) ) {
    if ( !strcmp( dtype, "i2" ) ) {
      static INT2TimeVectorSeries series;
      SUB( LALI2ReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI2WriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI2DestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "i4" ) ) {
      static INT4TimeVectorSeries series;
      SUB( LALI4ReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI4WriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI4DestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "i8" ) ) {
      static INT8TimeVectorSeries series;
      SUB( LALI8ReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI8WriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI8DestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u2" ) ) {
      static UINT2TimeVectorSeries series;
      SUB( LALU2ReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU2WriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU2DestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u4" ) ) {
      static UINT4TimeVectorSeries series;
      SUB( LALU4ReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU4WriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU4DestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u8" ) ) {
      static UINT8TimeVectorSeries series;
      SUB( LALU8ReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU8WriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU8DestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "s" ) ) {
      static REAL4TimeVectorSeries series;
      SUB( LALSReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALSWriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALSDestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "d" ) ) {
      static REAL8TimeVectorSeries series;
      SUB( LALDReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALDWriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALDDestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "c" ) ) {
      static COMPLEX8TimeVectorSeries series;
      SUB( LALCReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALCWriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALCDestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "z" ) ) {
      static COMPLEX16TimeVectorSeries series;
      SUB( LALZReadTVectorSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALZWriteTVectorSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALZDestroyVectorSequence( &stat, &(series.data) ), &stat );
    } else {
      ERROR( STREAMSERIESINPUTTESTC_EARG,
	     STREAMSERIESINPUTTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return STREAMSERIESINPUTTESTC_EARG;
    }
  }

  /* Read (and possibly write) the data as a time array series. */
  else if ( !strcmp( stype, "a" ) ) {
    if ( !strcmp( dtype, "i2" ) ) {
      static INT2TimeArraySeries series;
      SUB( LALI2ReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI2WriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI2DestroyArraySequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "i4" ) ) {
      static INT4TimeArraySeries series;
      SUB( LALI4ReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI4WriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI4DestroyArraySequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "i8" ) ) {
      static INT8TimeArraySeries series;
      SUB( LALI8ReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI8WriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI8DestroyArraySequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u2" ) ) {
      static UINT2TimeArraySeries series;
      SUB( LALU2ReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU2WriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU2DestroyArraySequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u4" ) ) {
      static UINT4TimeArraySeries series;
      SUB( LALU4ReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU4WriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU4DestroyArraySequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u8" ) ) {
      static UINT8TimeArraySeries series;
      SUB( LALU8ReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU8WriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU8DestroyArraySequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "s" ) ) {
      static REAL4TimeArraySeries series;
      SUB( LALSReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALSWriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALSDestroyArraySequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "d" ) ) {
      static REAL8TimeArraySeries series;
      SUB( LALDReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALDWriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALDDestroyArraySequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "c" ) ) {
      static COMPLEX8TimeArraySeries series;
      SUB( LALCReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALCWriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALCDestroyArraySequence( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "z" ) ) {
      static COMPLEX16TimeArraySeries series;
      SUB( LALZReadTArraySeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALZWriteTArraySeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALZDestroyArraySequence( &stat, &(series.data) ), &stat );
    } else {
      ERROR( STREAMSERIESINPUTTESTC_EARG,
	     STREAMSERIESINPUTTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return STREAMSERIESINPUTTESTC_EARG;
    }
  }

  /* Read (and possibly write) the data as a frequency series. */
  else if ( !strcmp( stype, "f" ) ) {
    if ( !strcmp( dtype, "i2" ) ) {
      static INT2FrequencySeries series;
      SUB( LALI2ReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI2WriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI2DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "i4" ) ) {
      static INT4FrequencySeries series;
      SUB( LALI4ReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI4WriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI4DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "i8" ) ) {
      static INT8FrequencySeries series;
      SUB( LALI8ReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALI8WriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALI8DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u2" ) ) {
      static UINT2FrequencySeries series;
      SUB( LALU2ReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU2WriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU2DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u4" ) ) {
      static UINT4FrequencySeries series;
      SUB( LALU4ReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU4WriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU4DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "u8" ) ) {
      static UINT8FrequencySeries series;
      SUB( LALU8ReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALU8WriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALU8DestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "s" ) ) {
      static REAL4FrequencySeries series;
      SUB( LALSReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALSWriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALSDestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "d" ) ) {
      static REAL8FrequencySeries series;
      SUB( LALDReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALDWriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALDDestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "c" ) ) {
      static COMPLEX8FrequencySeries series;
      SUB( LALCReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALCWriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALCDestroyVector( &stat, &(series.data) ), &stat );
    } else if ( !strcmp( dtype, "z" ) ) {
      static COMPLEX16FrequencySeries series;
      SUB( LALZReadFSeries( &stat, &series, fpIn ), &stat );
      if ( fpOut ) {
	SUB( LALZWriteFSeries( &stat, fpOut, &series ), &stat );
      }
      SUB( LALZDestroyVector( &stat, &(series.data) ), &stat );
    } else {
      ERROR( STREAMSERIESINPUTTESTC_EARG,
	     STREAMSERIESINPUTTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return STREAMSERIESINPUTTESTC_EARG;
    }
  }

  /* Unrecognized series type. */
  else {
    ERROR( STREAMSERIESINPUTTESTC_EARG,
	   STREAMSERIESINPUTTESTC_MSGEARG, 0 );
    LALPrintError( USAGE, *argv );
    return STREAMSERIESINPUTTESTC_EARG;
  }

  /* Close files, cleanup, and exit. */
  if ( infile && strcmp( infile, "stdin" ) )
    fclose( fpIn );
  if ( outfile && strcmp( outfile, "stdout" ) )
    fclose( fpOut );
  LALCheckMemoryLeaks();
  INFO( STREAMSERIESINPUTTESTC_MSGENORM );
  return STREAMSERIESINPUTTESTC_ENORM;
}
