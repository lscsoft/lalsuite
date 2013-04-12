/*
*  Copyright (C) 2010 Karsten Wiesner, Stas Babak, Duncan Brown, Eirini Messaritaki, Jolien Creighton, Reinhard Prix, Craig Robinson
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
 * File Name: Chisq_CPU.c
 *
 * Author: Wiesner, K., Anderson, W. G., and Brown, D. A.
 *
 *
 *-----------------------------------------------------------------------
 */
#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpChisq.h>
#include "Chisq_CPU.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

void Chisq_CPU (REAL4* chisq, COMPLEX8* q, COMPLEX8* qtilde, FindChirpChisqParams* params, 
    UINT4 numPoints, UINT4 numChisqBins, REAL4 chisqNorm, LALStatus UNUSED *status)
{
  
  UINT4* chisqBin     = params->chisqBinVec->data;
  COMPLEX8* qtildeBin = params->qtildeBinVec->data;
  int rc;

  for ( UINT4 l = 0; l < numChisqBins; ++l )
    {
      memset( qtildeBin, 0, numPoints * sizeof(COMPLEX8) );
      
      memcpy( qtildeBin + chisqBin[l], qtilde + chisqBin[l],
	      (chisqBin[l+1] - chisqBin[l]) * sizeof(COMPLEX8) );
      
      rc = XLALCOMPLEX8VectorFFT(params->qBinVecPtr[l], \
          params->qtildeBinVec, params->plan);
      if (rc != 0)
        XLAL_ERROR_VOID(rc);
    }
  
  memset( chisq, 0, numPoints * sizeof(REAL4) );
  
  for ( UINT4 j = 0; j < numPoints; ++j )
    {
      for ( UINT4 l = 0; l < numChisqBins; ++l )
	{
	  REAL4 Xl = crealf(params->qBinVecPtr[l]->data[j]);
	  REAL4 Yl = cimagf(params->qBinVecPtr[l]->data[j]);
	  REAL4 deltaXl = chisqNorm * Xl -
	    (chisqNorm * crealf(q[j]) / (REAL4) (numChisqBins));
	  REAL4 deltaYl = chisqNorm * Yl -
	    (chisqNorm * cimagf(q[j]) / (REAL4) (numChisqBins));
	  
	  chisq[j] += deltaXl * deltaXl + deltaYl * deltaYl;
	}
    }
  
}
