/* <lalVerbatim> */
#include <stdlib.h>
#include "LALStdlib.h"
#include "LALPrimer.h"

NRCSID( LALPRIMERTESTC, "$Id$" );

int lalDebugLevel = 0;

int
main( int argc, char **argv )
     /* Divides two numbers given on the command input line. */
{
  static LALStatus stat;
  REAL4 ratio;

  /* Parse input line. */
  if ( argc == 4 )
    lalDebugLevel = atoi( argv[3] );
  else if ( argc != 3 )
    {
      fprintf( stderr, "Usage: %s numer denom [ lalDebugLevel ]\n",
	       argv[0] );
      return 0; /* so that test script won't fail */
    }

  /* Compute ratio. */
  REAL4Divide( &stat, &ratio, atof( argv[1] ), atof( argv[2] ) );
  if ( stat.statusCode && ( lalDebugLevel > 0 ) )
    fprintf( stderr,
	     "Error[0] 1: program %s, file %s, line %i, %s\n"
	     "         Function REAL4Divide() failed\n",
	     argv[0], __FILE__, __LINE__, LALPRIMERTESTC );

  /* Print result. */
  if ( !stat.statusCode )
    fprintf( stdout, "\nRatio = %f\n", ratio );
  REPORTSTATUS( &stat );
  return stat.statusCode;
}
/* </lalVerbatim> */
