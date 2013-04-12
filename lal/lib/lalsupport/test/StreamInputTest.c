/*
*  Copyright (C) 2007 Teviet Creighton
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

/**
   \file
   \ingroup StreamInput_h
   \author Creighton, T. D.

\brief Reads a sequence or vector sequence from a file.

\heading{Usage}
\code
StreamInputTest [-o outfile] [-d debuglevel] [-t]
                [-v {ch | i2 | i4 | i8 | u2 | u4 | u8 | s | d} infile]
                [-s {ch | i2 | i4 | i8 | u2 | u4 | u8 | s | d | c | z} infile]
\endcode

\heading{Description}

This test program parses data from an input file or from \c stdin.
The following option flags are accepted:
<ul>
<li>[<tt>-o</tt>] Writes the output to \c outfile.  If
\c outfile is given as \c stdout, the data is written to
standard output (\e not to a file named \c stdout).  If the
<tt>-o</tt> flag is not given, the routines are exercised, but no output
is written.</li>
<li>[<tt>-d</tt>] Sets the debug level to \c debuglevel; if
absent, <tt>-d 0</tt> is assumed.</li>
<li>[<tt>-t</tt>] Writes to \c stderr the system time required to
read the file.</li>
<li>[<tt>-v</tt>] Reads the contents of \c infile as a sequence
of vectors to be parsed by the routines
<tt>LAL\<datatype\>ReadVectorSequence()</tt>, where <tt>\<datatype\></tt> is
determined by the argument immediately following the <tt>-v</tt> option
flag.  If \c infile is given as \c stdin, the data is read
from standard input (\e not from a file named \c stdin).</li>
<li>[<tt>-s</tt>] As <tt>-v</tt>, above, except that the file contents
are parsed by the routines <tt>LAL\<datatype\>ReadSequence()</tt>.  If
neither <tt>-v</tt> nor <tt>-s</tt> is specified,
<tt>-v s StreamInput.dat</tt> is assumed (this file is provided with the
distribution so that running the code with no arguments, \'a la
<tt>make check</tt>, will perform a nontrivial test of the algorithm).</li>
</ul>

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
*/

/** \name Error Codes */ /*@{*/
#define STREAMINPUTTESTC_ENORM 0        /**< Normal exit */
#define STREAMINPUTTESTC_ESUB  1        /**< Subroutine failed */
#define STREAMINPUTTESTC_EARG  2        /**< Error parsing arguments */
#define STREAMINPUTTESTC_EFILE 3        /**< Could not open file */
/*@}*/

/** \cond DONT_DOXYGEN */

#define STREAMINPUTTESTC_MSGENORM "Normal exit"
#define STREAMINPUTTESTC_MSGESUB  "Subroutine failed"
#define STREAMINPUTTESTC_MSGEARG  "Error parsing arguments"
#define STREAMINPUTTESTC_MSGEFILE "Could not open file"

#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/StreamInput.h>

/* Default parameter settings. */
extern int lalDebugLevel;
#define INFILE DATADIR "StreamInput.data"

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
		 __LINE__, "$Id$", statement ? statement : \
                 "", (msg) );                                        \
}                                                                    \
else (void)(0)

#define INFO( statement )                                            \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 "$Id$", (statement) );                    \
}                                                                    \
else (void)(0)

#define SUB( func, statusptr )                                       \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( STREAMINPUTTESTC_ESUB, STREAMINPUTTESTC_MSGESUB,            \
         "Function call \"" #func "\" failed:" );                    \
  return STREAMINPUTTESTC_ESUB;                                      \
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

#define GETCOMPLEX8FORMAT                                            \
do {                                                                 \
  if ( fpOut ) {                                                     \
    UINT4 nTot = values->length*dim;                                 \
    BOOLEAN neg = 0;                                                 \
    for ( i = 0; i < nTot; i++ )                                     \
      if ( crealf(values->data[i]) < 0.0 || cimagf(values->data[i]) < 0.0 ) \
	neg = 1;                                                     \
    if ( neg )                                                       \
      sprintf( format, "%%10.3e" );                                  \
    else                                                             \
      sprintf( format, "%%9.3e" );                                   \
  }                                                                  \
} while (0)

#define GETCOMPLEX16FORMAT                                           \
do {                                                                 \
  if ( fpOut ) {                                                     \
    UINT4 nTot = values->length*dim;                                 \
    BOOLEAN neg = 0;                                                 \
    for ( i = 0; i < nTot; i++ )                                     \
      if ( creal(values->data[i]) < 0.0 || cimag(values->data[i]) < 0.0 ) \
	neg = 1;                                                     \
    if ( neg )                                                       \
      sprintf( format, "%%10.3e" );                                  \
    else                                                             \
      sprintf( format, "%%9.3e" );                                   \
  }                                                                  \
} while (0)

/* Macros for printing sequences and vector sequences. */
#define PRINTCHARVECTORSEQUENCE                                      \
do {                                                                 \
  if ( fpOut )                                                       \
    for ( i = 0; i < values->length; i++ ) {                         \
      for ( j = 0; j < dim; j++ ) {                                  \
	int c = (int)( values->data[i*dim+j] );                      \
	if ( isgraph( c ) )                                          \
	  fputc( c, fpOut );                                         \
	else                                                         \
	  fputc( (int)(' '), fpOut );                                \
      }                                                              \
      fputc( (int)('\n'), fpOut );                                   \
    }                                                                \
} while (0)

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

#define PRINTCHARSEQUENCE                                            \
do {                                                                 \
  if ( fpOut )                                                       \
    for ( i = 0; i < values->length; i++ )                           \
      fputc( (int)( values->data[i] ), fpOut );                      \
} while (0)

#define PRINTSEQUENCE                                                \
do {                                                                 \
  if ( fpOut )                                                       \
    for ( i = 0; i < values->length; i++ ) {                         \
      fprintf( fpOut, format, values->data[i] );                     \
      fprintf( fpOut, "\n" );                                        \
    }                                                                \
} while (0)

#define PRINTCOMPLEX8SEQUENCE                                        \
do {                                                                 \
  if ( fpOut )                                                       \
    for ( i = 0; i < values->length; i++ ) {                         \
      fprintf( fpOut, format, crealf(values->data[i]) );             \
      fprintf( fpOut, " " );                                         \
      fprintf( fpOut, format, cimagf(values->data[i]) );             \
      fprintf( fpOut, "\n" );                                        \
    }                                                                \
} while (0)

#define PRINTCOMPLEX16SEQUENCE                                       \
do {                                                                 \
  if ( fpOut )                                                       \
    for ( i = 0; i < values->length; i++ ) {                         \
      fprintf( fpOut, format, creal(values->data[i]) );              \
      fprintf( fpOut, " " );                                         \
      fprintf( fpOut, format, cimag(values->data[i]) );              \
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

  lalDebugLevel = 0;

  /* Parse argument list.  arg stores the current position. */
  arg = 1;
  while ( arg < argc ) {
    /* Parse output file option. */
    if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	outfile = argv[arg++];
      } else {
	ERROR( STREAMINPUTTESTC_EARG, STREAMINPUTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return STREAMINPUTTESTC_EARG;
      }
    }
    /* Parse debug level option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	lalDebugLevel = atoi( argv[arg++] );
      } else {
	ERROR( STREAMINPUTTESTC_EARG, STREAMINPUTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return STREAMINPUTTESTC_EARG;
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
	ERROR( STREAMINPUTTESTC_EARG, STREAMINPUTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return STREAMINPUTTESTC_EARG;
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
	ERROR( STREAMINPUTTESTC_EARG, STREAMINPUTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return STREAMINPUTTESTC_EARG;
      }
    }
    /* Check for unrecognized options. */
    else {
      ERROR( STREAMINPUTTESTC_EARG, STREAMINPUTTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return STREAMINPUTTESTC_EARG;
    }
  } /* End of argument parsing loop. */

  /* Open input and output files. */
  if ( strcmp( infile, "stdin" ) ) {
    if ( !( fpIn = LALOpenDataFile( infile ) ) ) {
      ERROR( STREAMINPUTTESTC_EFILE, "- " STREAMINPUTTESTC_MSGEFILE,
	     infile );
      return STREAMINPUTTESTC_EFILE;
    }
  } else
    fpIn = stdin;
  if ( outfile ) {
    if ( strcmp( outfile, "stdout" ) ) {
      if ( !( fpOut = fopen( outfile, "w" ) ) ) {
	ERROR( STREAMINPUTTESTC_EFILE, "- " STREAMINPUTTESTC_MSGEFILE,
	       outfile );
	return STREAMINPUTTESTC_EFILE;
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
      GETCOMPLEX8FORMAT;
      PRINTCOMPLEX8SEQUENCE;
      SUB( LALCDestroyVector( &stat, &values ), &stat );
    } else if ( !strcmp( datatype, "z" ) ) {
      COMPLEX16Vector *values = NULL;
      start = clock();
      SUB( LALZReadSequence( &stat, &values, fpIn ), &stat );
      stop = clock();
      GETCOMPLEX16FORMAT;
      PRINTCOMPLEX16SEQUENCE;
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
  INFO( STREAMINPUTTESTC_MSGENORM );
  return STREAMINPUTTESTC_ENORM;
}
/** \endcond */
