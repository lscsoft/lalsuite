/************************ <lalVerbatim file="StreamSeriesInputTestCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{StreamSeriesInputTest.c}}
\label{ss:StreamSeriesInputTest.c}

Reads a \verb@REAL4VectorSequence@ from a file.

\subsubsection*{Usage}
\begin{verbatim}
StreamSeriesInputTest [-o outfile] [-d debuglevel] [-t]
        [{-s | -v | -a | -f} {i2 | i4 | i8 | u2 | u4 | u8 | s | d | c | z} infile]
\end{verbatim}

\subsubsection*{Description}

This test program parses data from an input file or from \verb@stdin@.
The following option flags are accepted:
\begin{itemize}
\item[\texttt{-o}] Writes the output to \verb@outfile@.  If
\verb@outfile@ is given as \verb@stdout@, the data is written to
standard output (\emph{not} to a file named \verb@stdout@).  If the
\verb@-o@ flag is not given, the routines are exercised, but no output
is written.
\item[\texttt{-d}] Sets the debug level to \verb@debuglevel@; if
absent, \verb@-d 0@ is assumed.
\item[\texttt{-t}] Writes to \verb@stderr@ the system time required to
read the file.
\item[\texttt{-v}] Reads the contents of \verb@infile@ as a sequence
of vectors to be parsed by the routines
\verb@LAL<datatype>ReadVectorSequence()@, where \verb@<datatype>@ is
determined by the argument immediately following the \verb@-v@ option
flag.  If \verb@infile@ is given as \verb@stdin@, the data is read
from standard input (\emph{not} from a file named \verb@stdin@).
\item[\texttt{-s}] As \verb@-v@, above, except that the file contents
are parsed by the routines \verb@LAL<datatype>ReadSequence()@.  If
neither \verb@-v@ nor \verb@-s@ is specified,
\verb@-v s StreamSeriesInput.dat@ is assumed (this file is provided with the
distribution so that running the code with no arguments, \'a la
\verb@make check@, will perform a nontrivial test of the algorithm).
\end{itemize}

For data read in as a character vector sequences, the output will
consist of a number of lines equal to the length of the sequence, with
each line being the length of the vector; all non-graphical characters
in the vector (including the various types of whitespace) will be
replaced with single spaces.  For character sequences, the output will
essentially be a copy of the input.  For numerical vector sequences,
the output will consist of separate lines for each vector of the
sequence, with each line printing the components of the vector in some
type-dependent format.  For numerical sequences, each line of output
contains a single number, or, in the case of complex datatypes, two
numbers representing the real and imaginary components, again in some
type-dependent format.

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
lalDebugLevel                   LALPrintError()
LALOpenDataFile()               LALCheckMemoryLeaks()
LALCHARReadVectorSequence()     LALCHARDestroyVectorSequence()
LALI2ReadVectorSequence()       LALI2DestroyVectorSequence()
LALI4ReadVectorSequence()       LALI4DestroyVectorSequence()
LALI8ReadVectorSequence()       LALI8DestroyVectorSequence()
LALU2ReadVectorSequence()       LALU2DestroyVectorSequence()
LALU4ReadVectorSequence()       LALU4DestroyVectorSequence()
LALU8ReadVectorSequence()       LALU8DestroyVectorSequence()
LALSReadVectorSequence()        LALSDestroyVectorSequence()
LALDReadVectorSequence()        LALDDestroyVectorSequence()
LALCHARReadSequence()           LALCHARDestroyVector()
LALI2ReadSequence()             LALI2DestroyVector()
LALI4ReadSequence()             LALI4DestroyVector()
LALI8ReadSequence()             LALI8DestroyVector()
LALU2ReadSequence()             LALU2DestroyVector()
LALU4ReadSequence()             LALU4DestroyVector()
LALU8ReadSequence()             LALU8DestroyVector()
LALSReadSequence()              LALSDestroyVector()
LALDReadSequence()              LALDDestroyVector()
LALCReadSequence()              LALCDestroyVector()
LALZReadSequence()              LALZDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{StreamSeriesInputTestCV}}

******************************************************* </lalLaTeX> */

#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/StreamInput.h>
#include <lal/StreamOutput.h>

NRCSID(STREAMSERIESINPUTTESTC,"$Id$");

/* Default parameter settings. */
int lalDebugLevel = 0;
#define INFILE "StreamSeriesInput.data"

/* Usage format string. */
#define USAGE "Usage: %s [-o outfile] [-d debuglevel] [-t]\n"        \
"\t[-v {ch | i2 | i4 | i8 | u2 | u4 | u8 | s | d} infile]\n"         \
"\t[-s {ch | i2 | i4 | i8 | u2 | u4 | u8 | s | d | c | z} infile]\n"

/* Valid datatypes for vectors or sequences. */
#define VTYPES "ch i2 i4 i8 u2 u4 u8 s d"
#define STYPES "ch i2 i4 i8 u2 u4 u8 s d c z"

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

/* Macro to determine the print format for integers. */
#define GETINTFORMAT                                                 \
do {                                                                 \
  if ( fpOut ) {                                                     \
    UINT4 nTot = values->length*dim;                                 \
    REAL8 max = 0.0;                                                 \
    BOOLEAN neg = 0;                                                 \
    for ( i = 0; i < nTot; i++ ) {                                   \
      REAL8 x = (REAL8)( values->data[i] );                          \
      REAL8 y = fabs( x );                                           \
      if ( x > max )                                                 \
        max = x;                                                     \
      if ( x != y )                                                  \
        neg = 1;                                                     \
    }                                                                \
    max += 0.5;                                                      \
    if ( neg )                                                       \
      max *= 10.0;                                                   \
    sprintf( format, "%%%ui", (UINT4)( log( max )/log( 10 ) ) + 2 ); \
  }                                                                  \
} while (0)

#define GETUINTFORMAT                                                \
do {                                                                 \
  if ( fpOut ) {                                                     \
    UINT4 nTot = values->length*dim;                                 \
    REAL8 max = 0.0;                                                 \
    BOOLEAN neg = 0;                                                 \
    for ( i = 0; i < nTot; i++ ) {                                   \
      REAL8 x = (REAL8)( values->data[i] );                          \
      REAL8 y = fabs( x );                                           \
      if ( x > max )                                                 \
        max = x;                                                     \
      if ( x != y )                                                  \
        neg = 1;                                                     \
    }                                                                \
    max += 0.5;                                                      \
    if ( neg )                                                       \
      max *= 10.0;                                                   \
    sprintf( format, "%%%uu", (UINT4)( log( max )/log( 10 ) ) + 1 ); \
  }                                                                  \
} while (0)

#define GETREALFORMAT                                                \
do {                                                                 \
  if ( fpOut ) {                                                     \
    UINT4 nTot = values->length*dim;                                 \
    BOOLEAN neg = 0;                                                 \
    for ( i = 0; i < nTot; i++ )                                     \
      if ( values->data[i] < 0.0 )                                   \
	neg = 1;                                                     \
    if ( neg )                                                       \
      sprintf( format, "%%10.3e" );                                  \
    else                                                             \
      sprintf( format, "%%9.3e" );                                   \
  }                                                                  \
} while (0)

#define GETCOMPLEXFORMAT                                             \
do {                                                                 \
  if ( fpOut ) {                                                     \
    UINT4 nTot = values->length*dim;                                 \
    BOOLEAN neg = 0;                                                 \
    for ( i = 0; i < nTot; i++ )                                     \
      if ( values->data[i].re < 0.0 || values->data[i].im < 0.0 )    \
	neg = 1;                                                     \
    if ( neg )                                                       \
      sprintf( format, "%%10.3e" );                                  \
    else                                                             \
      sprintf( format, "%%9.3e" );                                   \
  }                                                                  \
} while (0)

/* Macros for printing sequences and vector sequences. */
#define PRINTVECTORSEQUENCE                                          \
do {                                                                 \
  if ( fpOut )                                                       \
    for ( i = 0; i < values->length; i++ ) {                         \
      fprintf( fpOut, format, values->data[i*dim] );                 \
      for ( j = 1; j < dim; j++ ) {                                  \
        fprintf( fpOut, " " );                                       \
        fprintf( fpOut, format, values->data[i*dim+j] );             \
      }                                                              \
      fprintf( fpOut, "\n" );                                        \
    }                                                                \
} while (0)

#define PRINTSEQUENCE                                                \
do {                                                                 \
  if ( fpOut )                                                       \
    for ( i = 0; i < values->length; i++ ) {                         \
      fprintf( fpOut, format, values->data[i] );                     \
      fprintf( fpOut, "\n" );                                        \
    }                                                                \
} while (0)

#define PRINTCOMPLEXSEQUENCE                                         \
do {                                                                 \
  if ( fpOut )                                                       \
    for ( i = 0; i < values->length; i++ ) {                         \
      fprintf( fpOut, format, values->data[i].re );                  \
      fprintf( fpOut, " " );                                         \
      fprintf( fpOut, format, values->data[i].im );                  \
      fprintf( fpOut, "\n" );                                        \
    }                                                                \
} while (0)


/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

int
main(int argc, char **argv)
{
#if 1
  static LALStatus stat;       /* top-level status structure */
  FILE *fp;
  REAL8FrequencySeries series;

  lalDebugLevel = 7;

  if ( argc != 3 ) {
    fprintf( stderr, "%s: Does nothing at present\n", *argv );
    return 0;
  }
  memset( &series, 0, sizeof(REAL8FrequencySeries) );
  fp = fopen( argv[1], "r" );
  SUB( LALDReadFSeries( &stat, &series, fp ), &stat );
  fclose( fp );
  fp = fopen( argv[2], "w" );
  SUB( LALDWriteFSeries( &stat, fp, &series ), &stat );
  fclose( fp );
  SUB( LALDDestroyVector( &stat, &(series.data) ), &stat );
  LALCheckMemoryLeaks();
  return 0;

#else
  static LALStatus stat;       /* top-level status structure */
  INT4 arg;                    /* index over command-line options */
  UINT4 i, j;                  /* indecies */
  UINT4 dim = 1;               /* dimension of vector data */
  BOOLEAN timing = 0;          /* whether to perform timing tests */
  BOOLEAN vector = 1;          /* input is a sequence of vectors */
  CHAR format[256];            /* integer output format */
  CHAR *outfile = NULL;        /* name of output file */
  const CHAR *infile = INFILE; /* name of input file */
  const CHAR *datatype = "s";  /* data type tag */
  FILE *fpIn = NULL;           /* input file pointer */
  FILE *fpOut = NULL;          /* output file pointer */
  clock_t start = 0, stop = 0; /* data input timestamps */




  /* Parse argument list.  arg stores the current position. */
  arg = 1;
  while ( arg < argc ) {
    /* Parse output file option. */
    if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	outfile = argv[arg++];
      } else {
	ERROR( STREAMSERIESINPUTTESTC_EARG, STREAMSERIESINPUTTESTC_MSGEARG, 0 );
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
	ERROR( STREAMSERIESINPUTTESTC_EARG, STREAMSERIESINPUTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return STREAMSERIESINPUTTESTC_EARG;
      }
    }
    /* Parse timing option. */
    else if ( !strcmp( argv[arg], "-t" ) ) {
      arg++;
      timing = 1;
    }
    /* Parse vector sequence input option. */
    else if ( !strcmp( argv[arg], "-v" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	datatype = argv[arg++];
	infile = argv[arg++];
	vector = 1;
      } else {
	ERROR( STREAMSERIESINPUTTESTC_EARG, STREAMSERIESINPUTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return STREAMSERIESINPUTTESTC_EARG;
      }
    }
    /* Parse plain sequence input option. */
    else if ( !strcmp( argv[arg], "-s" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	datatype = argv[arg++];
	infile = argv[arg++];
	vector = 0;
      } else {
	ERROR( STREAMSERIESINPUTTESTC_EARG, STREAMSERIESINPUTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return STREAMSERIESINPUTTESTC_EARG;
      }
    }
    /* Check for unrecognized options. */
    else {
      ERROR( STREAMSERIESINPUTTESTC_EARG, STREAMSERIESINPUTTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return STREAMSERIESINPUTTESTC_EARG;
    }
  } /* End of argument parsing loop. */

  /* Open input and output files. */
  if ( strcmp( infile, "stdin" ) ) {
    if ( !( fpIn = LALOpenDataFile( infile ) ) ) {
      ERROR( STREAMSERIESINPUTTESTC_EFILE, "- " STREAMSERIESINPUTTESTC_MSGEFILE,
	     infile );
      return STREAMSERIESINPUTTESTC_EFILE;
    }
  } else
    fpIn = stdin;
  if ( outfile ) {
    if ( strcmp( outfile, "stdout" ) ) {
      if ( !( fpOut = fopen( outfile, "w" ) ) ) {
	ERROR( STREAMSERIESINPUTTESTC_EFILE, "- " STREAMSERIESINPUTTESTC_MSGEFILE,
	       outfile );
	return STREAMSERIESINPUTTESTC_EFILE;
      }
    } else
      fpOut = stdout;
  }

  /* Read the data as a vector sequence, and print it according to the
     selected datatype. */
  if ( vector ) {
    if ( !strcmp( datatype, "ch" ) ) {
      CHARVectorSequence *values = NULL;
      start = clock();
      SUB( LALCHARReadVectorSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      dim = values->vectorLength;
      PRINTCHARVECTORSEQUENCE;
      SUB( LALCHARDestroyVectorSequence( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "i2" ) ) {
      INT2VectorSequence *values = NULL;
      start = clock();
      SUB( LALI2ReadVectorSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      dim = values->vectorLength;
      GETINTFORMAT;
      PRINTVECTORSEQUENCE;
      SUB( LALI2DestroyVectorSequence( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "i4" ) ) {
      INT4VectorSequence *values = NULL;
      start = clock();
      SUB( LALI4ReadVectorSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      dim = values->vectorLength;
      GETINTFORMAT;
      PRINTVECTORSEQUENCE;
      SUB( LALI4DestroyVectorSequence( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "i8" ) ) {
      INT8VectorSequence *values = NULL;
      start = clock();
      SUB( LALI8ReadVectorSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      dim = values->vectorLength;
      GETINTFORMAT;
      PRINTVECTORSEQUENCE;
      SUB( LALI8DestroyVectorSequence( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "u2" ) ) {
      UINT2VectorSequence *values = NULL;
      start = clock();
      SUB( LALU2ReadVectorSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      dim = values->vectorLength;
      GETUINTFORMAT;
      PRINTVECTORSEQUENCE;
      SUB( LALU2DestroyVectorSequence( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "u4" ) ) {
      UINT4VectorSequence *values = NULL;
      start = clock();
      SUB( LALU4ReadVectorSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      dim = values->vectorLength;
      GETUINTFORMAT;
      PRINTVECTORSEQUENCE;
      SUB( LALU4DestroyVectorSequence( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "u8" ) ) {
      UINT8VectorSequence *values = NULL;
      start = clock();
      SUB( LALU8ReadVectorSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      dim = values->vectorLength;
      GETUINTFORMAT;
      PRINTVECTORSEQUENCE;
      SUB( LALU8DestroyVectorSequence( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "s" ) ) {
      REAL4VectorSequence *values = NULL;
      start = clock();
      SUB( LALSReadVectorSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      dim = values->vectorLength;
      GETREALFORMAT;
      PRINTVECTORSEQUENCE;
      SUB( LALSDestroyVectorSequence( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "d" ) ) {
      REAL8VectorSequence *values = NULL;
      start = clock();
      SUB( LALDReadVectorSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      dim = values->vectorLength;
      GETREALFORMAT;
      PRINTVECTORSEQUENCE;
      SUB( LALDDestroyVectorSequence( &stat, &values ), &stat );
    } else {
      ERROR( -1, "Internal consistency error!", 0 );
    }
  }

  /* Read the data as a plain sequence, and print it according to the
     selected datatype. */
  else {
    if ( !strcmp( datatype, "ch" ) ) {
      CHARVector *values = NULL;
      start = clock();
      SUB( LALCHARReadSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      PRINTCHARSEQUENCE;
      SUB( LALCHARDestroyVector( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "i2" ) ) {
      INT2Vector *values = NULL;
      start = clock();
      SUB( LALI2ReadSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      GETINTFORMAT;
      PRINTSEQUENCE;
      SUB( LALI2DestroyVector( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "i4" ) ) {
      INT4Vector *values = NULL;
      start = clock();
      SUB( LALI4ReadSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      GETINTFORMAT;
      PRINTSEQUENCE;
      SUB( LALI4DestroyVector( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "i8" ) ) {
      INT8Vector *values = NULL;
      start = clock();
      SUB( LALI8ReadSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      GETINTFORMAT;
      PRINTSEQUENCE;
      SUB( LALI8DestroyVector( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "u2" ) ) {
      UINT2Vector *values = NULL;
      start = clock();
      SUB( LALU2ReadSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      GETUINTFORMAT;
      PRINTSEQUENCE;
      SUB( LALU2DestroyVector( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "u4" ) ) {
      UINT4Vector *values = NULL;
      start = clock();
      SUB( LALU4ReadSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      GETUINTFORMAT;
      PRINTSEQUENCE;
      SUB( LALU4DestroyVector( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "u8" ) ) {
      UINT8Vector *values = NULL;
      start = clock();
      SUB( LALU8ReadSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      GETUINTFORMAT;
      PRINTSEQUENCE;
      SUB( LALU8DestroyVector( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "s" ) ) {
      REAL4Vector *values = NULL;
      start = clock();
      SUB( LALSReadSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      GETREALFORMAT;
      PRINTSEQUENCE;
      SUB( LALSDestroyVector( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "d" ) ) {
      REAL8Vector *values = NULL;
      start = clock();
      SUB( LALDReadSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      GETREALFORMAT;
      PRINTSEQUENCE;
      SUB( LALDDestroyVector( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "c" ) ) {
      COMPLEX8Vector *values = NULL;
      start = clock();
      SUB( LALCReadSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      GETCOMPLEXFORMAT;
      PRINTCOMPLEXSEQUENCE;
      SUB( LALCDestroyVector( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "z" ) ) {
      COMPLEX16Vector *values = NULL;
      start = clock();
      SUB( LALZReadSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      GETCOMPLEXFORMAT;
      PRINTCOMPLEXSEQUENCE;
      SUB( LALZDestroyVector( &stat, &values ), &stat );
    } else {
      ERROR( -1, "Internal consistency error!", 0 );
    }
  }

  /* Print timing information, if requested. */
  if ( timing )
    fprintf( stderr, "CPU time required to read data: %.2f s\n",
	     (double)( stop - start )/CLOCKS_PER_SEC );

  /* Close files, cleanup, and exit. */
  if ( infile && strcmp( infile, "stdin" ) )
    fclose( fpIn );
  if ( outfile && strcmp( outfile, "stdout" ) )
    fclose( fpOut );
  LALCheckMemoryLeaks();
  INFO( STREAMSERIESINPUTTESTC_MSGENORM );
  return STREAMSERIESINPUTTESTC_ENORM;

#endif

}
