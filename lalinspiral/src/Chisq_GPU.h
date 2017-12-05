/*
*  Copyright (C) 2010 Karsten Wiesner, Duncan Brown, Eirini Messaritaki, Jolien Creighton
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
 * File Name: Chisq_GPU.h
 *
 * Author: Wiesner, K., Anderson, W. G., and Brown, D. A. 
 *
 *
 *-----------------------------------------------------------------------
 */

/* We only want to generate SWIG wrappers if LALINSPIRAL_CUDA_ENABLED. */
#if !defined(SWIG) || defined(LALINSPIRAL_CUDA_ENABLED)

#ifndef _CHISQ_GPU_H
#define _CHISQ_GPU_H

void Chisq_GPU(REAL4* chisq, COMPLEX8* q, COMPLEX8* qtilde,  UINT4* chisqBin,
	       UINT4 numPoints, UINT4 numChisqBins, REAL4 chisqNorm);

#endif /* _CHISQ_GPU_H */
#endif /* !defined(SWIG) || defined(LALINSPIRAL_CUDA_ENABLED) */
