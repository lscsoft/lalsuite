/*********************************** <lalVerbatim file="InjectTestCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{InjectTest.c}}
\label{ss:InjectTest.c}

Injects an inspiral signal into Gaussian noise.

\subsubsection*{Usage}
\begin{verbatim}
InjectTest [-s sourcefile] [-d debuglevel]
\end{verbatim}

\subsubsection*{Description}

This test program does nothing at present.  The following option flags
are accepted (and then ignored):
\begin{itemize}
\item[\texttt{-s}] Reads source positions from the file
\verb@sourcefile@.
\item[\texttt{-d}] Sets the debug level to \verb@debuglevel@.
\end{itemize}

\subsubsection*{Exit codes}
****************************************** </lalLaTeX><lalErrTable> */
#define INJECTTESTC_ENORM 0
#define INJECTTESTC_ESUB  1
#define INJECTTESTC_EARG  2
#define INJECTTESTC_EVAL  3

#define INJECTTESTC_MSGENORM "Normal exit"
#define INJECTTESTC_MSGESUB  "Subroutine failed"
#define INJECTTESTC_MSGEARG  "Error parsing arguments"
#define INJECTTESTC_MSGEVAL  "Input argument out of valid range"
/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALPrintError()                 LALCheckMemoryLeaks()
LALMalloc()                     LALFree()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{InjectTestCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/Inject.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GeneratePPNInspiral.h>

NRCSID( INJECTTESTC, "$Id$" );

/* Default parameter settings. */
int lalDebugLevel = 0;

/* Usage format string. */
#define USAGE "Usage: %s [-s sourcefile] [-d debuglevel]\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, INJECTTESTC, statement ? statement :      \
                 "", (msg) );                                        \
}                                                                    \
while (0)

#define INFO( statement )                                            \
do                                                                   \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 INJECTTESTC, (statement) );                         \
}                                                                    \
while (0)

#define SUB( func, statusptr )                                       \
do                                                                   \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( INJECTTESTC_ESUB, INJECTTESTC_MSGESUB,                      \
         "Function call \"" #func "\" failed:" );                    \
  return INJECTTESTC_ESUB;                                           \
}                                                                    \
while (0)

#define CHECKVAL( val, lower, upper )                                \
do                                                                   \
if ( ( (val) <= (lower) ) || ( (val) > (upper) ) )                   \
{                                                                    \
  ERROR( INJECTTESTC_EVAL, INJECTTESTC_MSGEVAL,                      \
         "Value of " #val " out of range:" );                        \
  LALPrintError( #val " = %f, range = (%f,%f]\n", (REAL8)(val),      \
                 (REAL8)(lower), (REAL8)(upper) );                   \
  return INJECTTESTC_EVAL;                                           \
}                                                                    \
while (0)

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

int
main(int argc, char **argv)
{
  int arg;
  static LALStatus stat;
  CHAR *sourcefile = NULL;

  /* Parse argument list.  arg stores the current position. */
  arg = 1;
  while ( arg < argc ) {
    /* Parse source file option. */
    if ( !strcmp( argv[arg], "-s" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	sourcefile = argv[arg++];
      }else{
	ERROR( INJECTTESTC_EARG, INJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INJECTTESTC_EARG;
      }
    }
    /* Parse debug level option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	lalDebugLevel = atoi( argv[arg++] );
      }else{
	ERROR( INJECTTESTC_EARG, INJECTTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return INJECTTESTC_EARG;
      }
    }
    /* Check for unrecognized options. */
    else if ( argv[arg][0] == '-' ) {
      ERROR( INJECTTESTC_EARG, INJECTTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return INJECTTESTC_EARG;
    }
  } /* End of argument parsing loop. */


  /* I do absolutely nothing! */


  LALCheckMemoryLeaks();
  INFO( INJECTTESTC_MSGENORM );
  return INJECTTESTC_ENORM;
}
