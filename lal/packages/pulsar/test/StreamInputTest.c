/********************************* <lalVerbatim file="StreamInputTestCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{StreamInputTest.c}}
\label{ss:StreamInputTest.c}

Reads a \verb@REAL4VectorSequence@ from a file.

\subsubsection*{Usage}
\begin{verbatim}
StreamInputTest [-o outfile] [-d debuglevel] [infile]
\end{verbatim}

\subsubsection*{Description}

This test program reads a sequence of vectors of floating-point
numbers from an input file and writes them to \verb@stdout@.  The
following option flags are accepted:
\begin{itemize}
\item[\texttt{-o}] Writes the output to \verb@outfile@ instead of
\verb@stdout@.
\item[\texttt{-d}] Sets the debug level to \verb@debuglevel@.
\end{itemize}
Once the above options are processed, the remaining command-line
argument is the name of the input file.  If this is \verb@stdin@, the
data is read from standard input (\emph{not} from a file named
\verb@stdin@).  If the input file argument is missing, the data is
read from a file named \verb@StreamInput.dat@ provided with the
distribution; this is so that running the code with no arguments (\'a
la \verb@make check@) will perform a nontrivial test of the algorithm.

\subsubsection*{Exit codes}
****************************************** </lalLaTeX><lalErrTable> */
#define STREAMINPUTTESTC_ENORM 0
#define STREAMINPUTTESTC_ESUB  1
#define STREAMINPUTTESTC_EARG  2
#define STREAMINPUTTESTC_EFILE 3

#define STREAMINPUTTESTC_MSGENORM "Normal exit"
#define STREAMINPUTTESTC_MSGESUB  "Subroutine failed"
#define STREAMINPUTTESTC_MSGEARG  "Error parsing arguments"
#define STREAMINPUTTESTC_MSGEFILE "Could not open file"
/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALPrintError()                 LALCheckMemoryLeaks()
LALSReadVectorSequence()        LALDestroyVectorSequence()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{StreamInputTestCV}}

******************************************************* </lalLaTeX> */

#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include "StreamInput.h"

NRCSID(STREAMINPUTTESTC,"$Id$");

/* Default parameter settings. */
int lalDebugLevel = 0;
#define INFILE "StreamInput.dat"

/* Usage format string. */
#define USAGE "Usage: %s [-o outfile] [-d debuglevel] [infile]\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, STREAMINPUTTESTC, statement ? statement : \
                 "", (msg) );                                        \
}                                                                    \
else (void)(0)

#define INFO( statement )                                            \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 STREAMINPUTTESTC, (statement) );                    \
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

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

int
main(int argc, char **argv)
{
  static LALStatus stat; /* top-level status structure */
  INT4 arg;
  UINT4 i;               /* an index */
  CHAR *outfile = NULL;  /* name of output file */
  const CHAR *infile = INFILE; /* name of input file */
  REAL4VectorSequence *numbers = NULL; /* data read */
  FILE *fp;              /* input/output file pointer */

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
    /* Check for unrecognized options. */
    else if ( argv[arg][0] == '-' ) {
      ERROR( STREAMINPUTTESTC_EARG, STREAMINPUTTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return STREAMINPUTTESTC_EARG;
    }
    /* Parse remaining parameters. */
    else {
      if ( argc == arg + 1 ) {
	if ( !strcmp( argv[arg++], "stdin" ) ) {
	  infile = NULL;
	} else {
	  infile = argv[arg++];
	}
      } else {
	ERROR( STREAMINPUTTESTC_EARG, STREAMINPUTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return STREAMINPUTTESTC_EARG;
      }
    }
  } /* End of argument parsing loop. */

  /* Open input file. */
  if ( infile ) {
    if ( !( fp = fopen( infile, "r" ) ) ) {
      ERROR( STREAMINPUTTESTC_EFILE, "- " STREAMINPUTTESTC_MSGEFILE,
	     infile );
      return STREAMINPUTTESTC_EFILE;
    }
  } else
    fp = stdin;

  /* Read data and close file. */
  SUB( LALSReadVectorSequence( &stat, &numbers, fp ), &stat );
  if ( infile )
    fclose( fp );

  /* Open output file. */
  if ( outfile ) {
    if ( !( fp = fopen( outfile, "r" ) ) ) {
      ERROR( STREAMINPUTTESTC_EFILE, "- " STREAMINPUTTESTC_MSGEFILE,
	     outfile );
      return STREAMINPUTTESTC_EFILE;
    }
  } else
    fp = stdout;

  /* Write data and close file. */
  i = numbers->length;
  {
    REAL4 *data = numbers->data;
    while ( i-- ) {
      UINT4 j = numbers->vectorLength;
      while ( j-- )
	fprintf( fp, "%10.3e ", *(data++) );
      fprintf( fp, "\n" );
    }
  }
  if ( outfile )
    fclose( fp );

  /* Clean up memory and exit. */
  SUB( LALSDestroyVectorSequence( &stat, &numbers ), &stat );
  LALCheckMemoryLeaks();
  INFO( STREAMINPUTTESTC_MSGENORM );
  return STREAMINPUTTESTC_ENORM;
}
