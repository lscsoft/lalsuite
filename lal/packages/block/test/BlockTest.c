/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton
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

#include <lal/LALRCSID.h>
NRCSID (BLOCKTESTC,"$Id$");

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
