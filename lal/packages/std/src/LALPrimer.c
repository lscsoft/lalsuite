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
#include <lal/LALStdlib.h>
#include <lal/LALPrimer.h>

void
REAL4Invert( LALStatus *stat, REAL4 *output, REAL4 input )
     /* Computes the inverse of a REAL4 number. */
{
  INITSTATUS(stat);

  /* This traps coding errors in the calling routine. */
  ASSERT( output != NULL, stat, LALPRIMERH_ENULL, LALPRIMERH_MSGENULL );

  /* This traps runtime errors. */
  if ( input == 0.0 )
    ABORT( stat, LALPRIMERH_EDIV0, LALPRIMERH_MSGEDIV0 );

  *output = 1.0/input;
  RETURN( stat );
}


void
REAL4Divide( LALStatus *stat, REAL4 *output, REAL4 numer, REAL4 denom )
     /* Computes the ratio of two REAL4 numbers. */
{
  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  TRY( REAL4Invert( stat->statusPtr, output, denom ), stat );
  *output *= numer;

  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
/* </lalVerbatim> */
