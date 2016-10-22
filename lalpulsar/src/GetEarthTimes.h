/*
*  Copyright (C) 2007 David M. Whitbeck, Jolien Creighton, Ian Jones, Reinhard Prix, Teviet Creighton
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


#ifndef _GETEARTHTIMES_H
#define _GETEARTHTIMES_H

#include <lal/LALStdlib.h>
#include <lal/LALBarycenter.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup GetEarthTimes_h  Header GetEarthTimes.h
 * \author Creighton, T. D.
 * \ingroup lalpulsar_general
 * \brief Computes the next sidereal midnight and autumnal equinox.
 */
/*@{*/

int XLALGetEarthTimes( const LIGOTimeGPS *tepoch, REAL8 *tMidnight, REAL8 *tAutumn );

/*@}*/

#ifdef __cplusplus
}
#endif

#endif
