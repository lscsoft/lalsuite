/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpTDData.c
 *
 * Author: Brown D. A., and Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0 
<lalVerbatim file="FindChirpTDDataCV">
Author: Brown, D. A., and Creighton, J. D. E.
$Id$
</lalVerbatim> 

<lalLaTeX>
\subsection{Module \texttt{FindChirpTDData.c}}
\label{ss:FindChirpTDData.c}

Time domain filtering code.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpTDDataCP}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpTDDataCV}}
</lalLaTeX> 
#endif

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpTD.h>


NRCSID (FINDCHIRPTDDATAC, "$Id$");


/* <lalVerbatim file="FindChirpTDDataCP"> */
void
LALFindChirpTDData (
    LALStatus                  *status,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec,
    FindChirpDataParams        *params
    )
/* </lalVerbatim> */
{
  Approximant approx;
  UINT4 i, k;

  INITSTATUS( status, "LALFindChirpTDData", FINDCHIRPTDDATAC );
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
  params->approximant = TaylorF2;
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
