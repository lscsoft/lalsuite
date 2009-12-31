/*
*  Copyright (C) 2007 Jolien Creighton
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

/* <lalVerbatim> */
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALPrimer.h>

NRCSID( LALPRIMERTESTC, "$Id$" );

extern int lalDebugLevel;

int
main( int argc, char **argv )
     /* Divides two numbers given on the command input line. */
{
  static LALStatus stat;
  REAL4 ratio;

  lalDebugLevel = 0;

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
