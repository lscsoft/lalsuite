/**** <lalVerbatim file="BlockTestCV">
 * Author: Matthew M. Tibbits
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * \subsection{Program \texttt{BlockTest.c}}
 * 
 * \subsubsection*{Usage}
 *
 * \begin{verbatim}
 * BlockTest
 * \end{verbatim}
 *
 * \subsubsection*{Description}
 *
 * \vfill{\footnotesize\input{BlockTestCV}}
 **** </lalLaTeX> */

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/BlockRho.h>

#define TestStatus( ps ) \
  if ( (ps)->statusCode ) { \
    fprintf( stderr, "Failed LAL routine near line %d\n", __LINE__ ); \
    exit( 1 ); \
  } else ((void)0)

int lalDebugLevel = LALMSGLVL3;

int main( void )
{
  return 77;
}
