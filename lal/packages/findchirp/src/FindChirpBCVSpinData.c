/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpBCVSpinData.c
 *
 * Author: Brown D. A., Spinning BCV-Modifications: Jones, G
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0 
<lalVerbatim file="FindChirpBCVSpinDataCV">
Author: Brown, D. A., Spinning BCV-Modifications: Jones, G.
$Id$
</lalVerbatim> 

<lalLaTeX>
\subsection{Module \texttt{FindChirpBCVSpinData.c}}
\label{ss:FindChirpBCVSpinData.c}

\input{FindChirpBCVSpinDataCDoc}

\vfill{\footnotesize\input{FindChirpBCVSpinDataCV}}
</lalLaTeX> 
#endif

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>


NRCSID (FINDCHIRPBCVSPINDATAC, "$Id$");

/*documenation later*/
void
LALFindChirpBCVSpinData (
    LALStatus                  *status,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec,
    FindChirpSPDataParams      *params
    )

{
  UINT4                 i, k;
  UINT4                 cut;

  REAL4                *w;
  REAL4                *amp;
  REAL4                *ampBCV; 
  COMPLEX8             *wtilde;
  REAL4		       *tmpltPower;
  REAL4		       *tmpltPowerBCV;

  REAL4Vector          *dataVec;
  REAL4                *spec;
  COMPLEX8             *resp;

  COMPLEX8             *outputData;
  COMPLEX8             *outputData1;
  COMPLEX8             *outputData2;
  COMPLEX8             *outputData3;
  COMPLEX8             *outputDataBCV;

  UINT4                *chisqBin    = NULL;
  UINT4                *chisqBinBCV = NULL;
  UINT4                 numChisqBins;
  UINT4                 chisqPt;
  REAL4                 increment;
  REAL4                 nextBin;
  REAL4                 partSum;
  REAL4                 Power  = 0.0 ;
  REAL4                 PowerBCV = 0.0 ;
  REAL4                 I = 0.0;
  REAL4                 J = 0.0;
  REAL4                 K = 0.0;
  REAL4                 L = 0.0;
  REAL4                 M = 0.0;

  REAL4                 Beta; /* Spin parameter, initialise? */  
  REAL4                 denominator;
  REAL4                 denominator1;

  FindChirpSegment     *fcSeg;
  DataSegment          *dataSeg;



  /*declaration*/
  INITSTATUS( status, "LALFindChirpBCVSpinData", FINDCHIRPBCVSPINDATAC );
  ATTATCHSTATUSPTR( status );


  /*code*/

  /* compute strain */
  for ( k = 0; k < fcSeg->data->data->length; ++k )
  {
    REAL4 p = outputData[k].re;
    REAL4 q = outputData[k].im;
    REAL4 x = resp[k].re * params->dynRange;
    REAL4 y = resp[k].im * params->dynRange;

    outputData[k].re =  p*x - q*y;
    outputData[k].im =  p*y + q*x;
  }

  /* 
   * moments necessary for the calculation of
   * the BCVSpin normalization parameters
   */

  for ( k = 1; k < fcSeg->data->data->length; ++k )
  {
    fcSeg->segNorm += amp[k] * amp[k] * wtilde[k].re ; /* for std-candle */
    I += 4.0 * amp[k] * amp[k] * wtilde[k].re ;
    J += 4.0 * amp[k] * amp[k] * wtilde[k].re * 
      cos(Beta * amp[k] / ampBCV[k]); /* amp[k]/ampBCV[k] = f^-2/3  */
    K += 4.0 * amp[k] * amp[k] * wtilde[k].re * 
      sin(Beta * amp[k] / ampBCV[k]);
    L += 4.0 * amp[k] * amp[k] * wtilde[k].re * 
      sin(2 * Beta * amp[k] / ampBCV[k]);
    M += 4.0 * amp[k] * amp[k] * wtilde[k].re * 
      cos(2 * Beta * amp[k] / ampBCV[k]);
  }

  /* 
   * the calculation of the orthonormalised 
   * amplitude vectors A1, A2, A3 (using Eirini's
   * a1, b1, b2 for time being)
   *
   */   

  denominator = I*M  +  0.5*pow(I,2) - pow(J,2);
  denominator1 = sqrt ( 0.25 * pow(I,3) + M*(pow(J,2) - 
        pow(K,2)) - 0.5*(pow(J,2) + pow(K,2)) - I*(pow(L,2) + 
        pow(M,2)) + 2*J*K*L );

  fcSeg->a1 = 1.0 / sqrt(I) ;

  fcSeg->b1 = 1.0 /sqrt(denominator) * (I * cos(Beta * amp[k]) -  J);

  fcSeg->b2 = 1.0 /denominator1 * ( sin(Beta * amp[k]) - 
      (I*L - J*K)*cos(Beta * amp[k])/denominator + 
      (J*L - K*M + 0.5*I*K)/denominator );
  /* note that these quntities need to be multiplied by amp[k] to match */
  /* definitions in documentation */

  /*
   *
   * conditioning data for FindChirpFilter.c
   * note lack of exponential terms
   * Also, amp[k] factor now included
   *
   */

  outputData1[k].re *= fcSeg->a1 * amp[k] * wtilde[k].re;
  outputData2[k].re *= fcSeg->b1 * amp[k] * wtilde[k].re;
  outputData3[k].re *= fcSeg->b2 * amp[k] * wtilde[k].re;


  DETATCHSTATUSPTR( status );
  RETURN( status );

}
