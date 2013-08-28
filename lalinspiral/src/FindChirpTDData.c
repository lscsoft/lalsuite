/*
*  Copyright (C) 2007 Duncan Brown, Gareth Jones, Jolien Creighton, Craig Robinson
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
 * File Name: FindChirpTDData.c
 *
 * Author: Brown D. A., and Creighton, J. D. E.
 *
 *-----------------------------------------------------------------------
 */

/**
 * \author Brown, D. A., and Creighton, J. D. E.
 * \file
 * \ingroup FindChirpTD_h
 *
 * \brief Time domain filtering code.
 *
 * \heading{Description}
 *
 * \heading{Algorithm}
 *
 * \heading{Uses}
 *
 * \heading{Notes}
 *
 */

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpTD.h>
#include <lal/FindChirpSP.h>

void
LALFindChirpTDData (
    LALStatus                  *status,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec,
    FindChirpDataParams        *params
    )

{
  Approximant approx;
  UINT4 i, k;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( params, status,
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );
  ASSERT( params->ampVec, status,
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );
  ASSERT( params->ampVec->data, status,
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );

  /* check the approximant */
  switch ( params->approximant )
  {
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case GeneratePPN:
    case PadeT1:
    case EOB:
    case EOBNR:
    case EOBNRv2:
    case BCVSpin:
    case FindChirpPTF:
    case IMRPhenomB:
    case AmpCorPPN:
      /* store the input approximant */
      approx = params->approximant;
      break;

    default:
      ABORT( status, FINDCHIRPTDH_EMAPX, FINDCHIRPTDH_MSGEMAPX );
      break;
  }

  /* set the amplitude vector to unity */
  for ( k = 0; k < params->ampVec->length; ++k )
  {
    params->ampVec->data[k] = 1.0;
  }

  /* call the stationary phase conditioning to do the real work */
  params->approximant = FindChirpSP;
  LALFindChirpSPData ( status->statusPtr, fcSegVec, dataSegVec, params );
  CHECKSTATUSPTR( status );
  params->approximant = approx;

  for ( i = 0; i < fcSegVec->length; ++i )
  {
    FindChirpSegment *fcSeg = fcSegVec->data + i;

    /* store the waveform approximant in the data segment */
    fcSeg->approximant = params->approximant;

    /* zero the tmpltPower and segNorm vectors as they are incorrect */
    memset( fcSeg->tmpltPowerVec->data, 0,
        fcSeg->tmpltPowerVec->length * sizeof(REAL4) );
    memset( fcSeg->segNorm->data, 0,
        fcSeg->segNorm->length * sizeof(REAL4) );
  }

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
