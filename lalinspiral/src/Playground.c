/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Stephen Fairhurst
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

/*-----------------------------------------------------------------------
 *
 * File Name: Playground.c
 *
 * Author: Brown, D. A.
 *
 *-----------------------------------------------------------------------
 */

#include <lal/LALInspiral.h>
#include <lal/Date.h>


/** \ingroup Date_h
    \author Brown D. A.

    \brief Determines if a given time is playground data.

\heading{Algorithm}

The playground algorithm is given in LIGO techincal document T030020-01.
Briefly, \f$t\f$ is playground if
\f{equation}{
t - 729273613 \% 6370 < 600.
\f}

*/
int
XLALINT8NanoSecIsPlayground (
    INT8         ns
    )

{
  const INT8 start = 729273613 * XLAL_BILLION_INT8;
  const INT8 interval = 6370 * XLAL_BILLION_INT8;
  const INT8 length = 600 * XLAL_BILLION_INT8;

  return (ns - start) % interval < length;
}
