/* <lalVerbatim> */
#include "LALStdlib.h"
#include "LALPrimer.h"

NRCSID( LALPRIMERC, "$Id$" );

void
REAL4Invert( Status *stat, REAL4 *output, REAL4 input )
     /* Computes the inverse of a REAL4 number. */
{
  INITSTATUS( stat, "REAL4Invert", LALPRIMERC );

  /* This traps coding errors in the calling routine. */
  ASSERT( output != NULL, stat, LALPRIMERH_ENULL, LALPRIMERH_MSGENULL );

  /* This traps runtime errors. */
  if ( input == 0.0 )
    ABORT( stat, LALPRIMERH_EDIV0, LALPRIMERH_MSGEDIV0 );

  *output = 1.0/input;
  RETURN( stat );
}


void
REAL4Divide( Status *stat, REAL4 *output, REAL4 numer, REAL4 denom )
     /* Computes the ratio of two REAL4 numbers. */
{
  INITSTATUS( stat, "REAL4Divide", LALPRIMERC );
  ATTATCHSTATUSPTR( stat );

  TRY( REAL4Invert( stat->statusPtr, output, denom ), stat );
  *output *= numer;

  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
/* </lalVerbatim> */
