/******************************** <lalVerbatim file="LALSampleTestCV">
Author: Creighton, T. D.
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{LALSampleTest.c}}
\label{ss:LALSampleTest.c}

Example program for LAL.

\subsubsection*{Usage}
\begin{verbatim}
LALSampleTest [numer denom [lalDebugLevel]]
\end{verbatim}

\subsubsection*{Description}

This program demonstrates LAL coding and documentation standards for
test programs.  It reads two numbers \verb@numer@ and \verb@denom@
from the command line, computes their quotient using the function
\verb@LALREAL8Divide()@, and prints the result to \verb@stdout@.  It is
run with \verb@lalDebugLevel@ = 0, unless set by the optional third
argument.

\subsubsection*{Exit codes}
\input{LALSampleTestCE}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALPrintError()
LALREAL8Divide()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALSampleTestCV}}
******************************************************* </lalLaTeX> */

#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALSample.h>

NRCSID( LALSAMPLETESTC, "$Id$" );

/***************************** <lalErrTable file="LALSampleTestCE"> */
#define LALSAMPLETESTC_ENOM 0
#define LALSAMPLETESTC_EARG 1
#define LALSAMPLETESTC_ESUB 2
#define LALSAMPLETESTC_MSGENOM "Nominal exit"
#define LALSAMPLETESTC_MSGEARG "Error parsing command-line arguments"
#define LALSAMPLETESTC_MSGESUB "Subroutine returned error"
/***************************** </lalErrTable> */

/* Declare and set the default lalDebugLevel */
int lalDebugLevel = 0;

/* A local macro for printing error messages */
#define EXIT( code, program, message )                                \
  do {                                                                \
    if (( lalDebugLevel & LALERROR ) && (code))                          \
      LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"\
                     "        %s\n", (code), (program), __FILE__,     \
                     __LINE__, LALSAMPLETESTC, (message) );           \
    else if ( lalDebugLevel & LALINFO )                                  \
      LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"    \
                     "        %s\n", (program), __FILE__, __LINE__,   \
                     LALSAMPLETESTC, (message) );                     \
    return (code);                                                    \
  } while (0)

/* The main function */
int
main( int argc, char **argv )
{
  static LALStatus stat;
  REAL8 ratio;

  /* Parse input line. */
  if ( argc == 1 )
    EXIT( LALSAMPLETESTC_ENOM, argv[0], LALSAMPLETESTC_MSGENOM );
  else if ( argc == 4 )
    lalDebugLevel = atoi( argv[3] );
  else if ( argc != 3 )
    {
      LALPrintError( "Usage: %s [numer denom [ lalDebugLevel ]]\n",
		     argv[0] );
      EXIT( LALSAMPLETESTC_EARG, argv[0], LALSAMPLETESTC_MSGEARG );
    }

  /* Compute ratio. */
  LALREAL8Divide( &stat, &ratio, atof( argv[1] ), atof( argv[2] ) );
  if ( stat.statusCode )
    EXIT( LALSAMPLETESTC_ESUB, argv[0], LALSAMPLETESTC_MSGESUB );

  /* Print result. */
  printf( "Ratio: %f / %f = %f\n", atof( argv[1] ), atof( argv[2] ),
	  ratio );
  EXIT( LALSAMPLETESTC_ENOM, argv[0], LALSAMPLETESTC_MSGENOM );
}
