/*
*  Copyright (C) 2007 Karsten Wiesner, Duncan Brown, Eirini Messaritaki, Jolien Creighton
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
 * File Name: Chisq_CPU.h
 *
 * Author: Wiesner, K. Anderson, W. G., and Brown, D. A.
 *
 *
 *-----------------------------------------------------------------------
 */
#ifndef _CHISQ_CPU_H
#define _CHISQ_CPU_H

#include <lal/FindChirp.h>

void Chisq_CPU (REAL4* chisq, COMPLEX8* q, COMPLEX8* qtilde, FindChirpChisqParams* params, 
		UINT4 numPoints, UINT4 numChisqBins, REAL4 chisqNorm, LALStatus *status);

#endif /* _CHISQ_CPU_H */
