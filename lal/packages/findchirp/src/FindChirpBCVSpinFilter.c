 /*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpBCVSpinFilter.c
 *
 * Author: Brown D. A., Spinning BCV-Modifications: Jones, G
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0 
<lalVerbatim file="FindChirpBCVSpinFilterCV">
Author: Brown, D. A., Spinning BCV-Modifications: Jones, G.
$Id$
</lalVerbatim> 

<lalLaTeX>
\subsection{Module \texttt{FindChirpBCVSpinFilter.c}}
\label{ss:FindChirpBCVSpinFilter.c}

\input{FindChirpBCVSpinFilterCDoc}

\vfill{\footnotesize\input{FindChirpBCVSpinFilterCV}}
</lalLaTeX> 
#endif

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>


NRCSID (FINDCHIRPBCVSPINFILTERC, "$Id$");

/*documenation later*/
void
LALFindChirpBCVSpinFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params,             
    FindChirpSPDataParams      *inputParams,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec
  )

{
  UINT4                 j, k;
  UINT4                 numPoints;
  UINT4                 deltaEventIndex;
  UINT4                 ignoreIndex;
  REAL4                 deltaT;
  REAL4                 norm;
  REAL4                 modqsqThresh;
  REAL4                 rhosqThresh;
  REAL4                 mismatch;
  REAL4                 chisqThreshFac;
  REAL4                 modChisqThresh;
  UINT4                 numChisqBins;
  UINT4                 eventStartIdx = 0;
  REAL4                 chirpTime     = 0;
  BOOLEAN               haveChisq     = 0;
  COMPLEX8             *qtilde        = NULL; 
  COMPLEX8             *qtildeBCV     = NULL; 
  COMPLEX8             *q             = NULL; 
  COMPLEX8             *qBCV          = NULL;
  COMPLEX8             *inputData     = NULL;
  COMPLEX8             *inputDataBCV  = NULL;
  COMPLEX8             *tmpltSignal   = NULL;
  SnglInspiralTable    *thisEvent     = NULL;
  LALMSTUnitsAndAcc     gmstUnits;
  FindChirpSegment      *fcSeg;
  DataSegment           *dataSeg; 
  REAL4                 templateNorm;
  REAL4                 modqsq;
  COMPLEX8              *wtilde;  /* need new pointer name? */
  REAL4                 *amp;
  REAL4                 *ampBCV;
  REAL4                 I = 0.0;
  REAL4                 J = 0.0;
  REAL4                 K = 0.0;
  REAL4                 L = 0.0;
  REAL4                 M = 0.0;
  REAL4                 Beta; /* Spin parameter, value from bank or external loop */  
  REAL4                 denominator;
  REAL4                 denominator1;
  REAL4                 a1;
  REAL4                 a2;                  
  REAL4                 a3;     
  COMPLEX8             *inputData1;
  COMPLEX8             *inputData2;
  COMPLEX8             *inputData3;
  FindChirpChisqInput  *chisqInput;
  FindChirpChisqInput  *chisqInputBCV;

  INITSTATUS( status, "LALFindChirpBCVSpinFilter", FINDCHIRPBCVSPINFILTERC );
  ATTATCHSTATUSPTR( status );


  /*
   *    
   * check that the arguments are reasonable
   * may need to remove asserts regarding chisq
   *          
   */

  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( eventList, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( !*eventList, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the filter parameters are reasonable */
  ASSERT( params->deltaT > 0, status,
      FINDCHIRPH_EDTZO, FINDCHIRPH_MSGEDTZO );
  ASSERT( params->rhosqThresh > 0, status,
      FINDCHIRPH_ERHOT, FINDCHIRPH_MSGERHOT );
  ASSERT( params->chisqThresh > 0, status,
      FINDCHIRPH_ECHIT, FINDCHIRPH_MSGECHIT );

  /* check that the fft plan exists */
  ASSERT( params->invPlan, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the workspace vectors exist */
  ASSERT(params->qVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qVec->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVec->data,status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL);
  ASSERT(params->qVecBCV, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qVecBCV->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVecBCV, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVecBCV->data,status, FINDCHIRPH_ENULL, 
	  FINDCHIRPH_MSGENULL);
  
  /* check that the chisq parameter and input structures exist */
  ASSERT( params->chisqParams, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->chisqInput,   status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->chisqInputBCV,status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* if a rhosqVec vector has been created, check we can store data in it */
  if ( params->rhosqVec )
  {
    ASSERT( params->rhosqVec->data->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
    ASSERT( params->rhosqVec->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* if a chisqVec vector has been created, check we can store data in it */
  if ( params->chisqVec )
  {
    ASSERT( params->chisqVec->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* make sure that the input structure exists */
  ASSERT( input, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the input structure contains some input */
  ASSERT( input->tmplt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->fcTmplt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->segment, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the template and the segment are both BCV */
  ASSERT( input->fcTmplt->approximant == BCVSpin, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  ASSERT( input->segment->approximant == BCVSpin, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );

/*
   *
   * point local pointers to input and output pointers
   * Check that I actually need all this
   *
   */


  /* workspace vectors */
  q    = params->qVec->data;
  qBCV = params->qVecBCV->data;
  qtilde    = params->qtildeVec->data;
  qtildeBCV = params->qtildeVecBCV->data;
  
  
  /* template and data */
  inputData     = input->segment->data->data->data;
  inputDataBCV  = input->segment->dataBCV->data->data;
  tmpltSignal   = input->fcTmplt->data->data;
  templateNorm  = input->fcTmplt->tmpltNorm;   /* this is expPsi */
  deltaT        = params->deltaT;


/*
 *
 * code which prev. would be in data function
 *
 */




  amp        = inputParams->ampVec->data;
  ampBCV     = inputParams->ampVecBCV->data;
  wtilde     = inputParams->wtildeVec->data;

  for ( k = 1; k < fcSeg->data->data->length; ++k )
  {
    I += 4.0 * amp[k] * amp[k] * wtilde[k].re ;
    J += 4.0 * amp[k] * amp[k] * wtilde[k].re * 
      cos(Beta * amp[k] / ampBCV[k]);                
    K += 4.0 * amp[k] * amp[k] * wtilde[k].re * 
      sin(Beta * amp[k] / ampBCV[k]);
    L += 4.0 * amp[k] * amp[k] * wtilde[k].re * 
      sin(2 * Beta * amp[k] / ampBCV[k]);
    M += 4.0 * amp[k] * amp[k] * wtilde[k].re * 
      cos(2 * Beta * amp[k] / ampBCV[k]);
  }

 denominator = I*M  +  0.5*pow(I,2) - pow(J,2);
 denominator1 = sqrt ( 0.25 * pow(I,3) + M*(pow(J,2) - 
        pow(K,2)) - 0.5*(pow(J,2) + pow(K,2)) - I*(pow(L,2) + 
        pow(M,2)) + 2*J*K*L );

  /* 
   * the calculation of the orthonormalised 
   * amplitude vectors a1, a2, a3, put in a loop
   *
   */   


 for ( k = 1; k < fcSeg->data->data->length; ++k )
  {
  a1 = 1.0 * amp[k]/ sqrt(I) ;

  a2 = 1.0 * amp[k]/sqrt(denominator) * (I * cos(Beta * amp[k]) -  J);

  a3 = 1.0 * amp[k]/denominator1 * ( sin(Beta * amp[k]) - 
      (I*L - J*K)*cos(Beta * amp[k])/denominator + 
      (J*L - K*M + 0.5*I*K)/denominator );
  }
  
 

  /*
   * initialising outputData vectors to
   * calibrated detector output as calc in LAL..Data()
   * note lack of exponential terms, these are
   * calc in LALFindChirpBCVSpinTemplate()
   */



  inputData1 = fcSeg->data->data->data;
  /*  inputData2 = fcSeg->data->data->data;
      inputData3 = fcSeg->data->data->data; i dont think ill need these   */
  
  /*
   *
   * compute qtilde, qtildeBCV, and q, qBCV 
   * need to create qtildeBCVSpin, qBCVSpin
   *
   * needs close checking
   *
   */


  memset( qtilde,    0, numPoints * sizeof(COMPLEX8) );
  memset( qtildeBCV, 0, numPoints * sizeof(COMPLEX8) );

  /* qtilde positive frequency, not DC or nyquist */
  for ( k = 1; k < numPoints/2; ++k )
  {
    REAL4 r        = inputData1[k].re;
    REAL4 s        = 0.0 - inputData1[k].im;    /* note complex conjugate */
    /* REAL4 rBCV     = inputData2[k].re;
       REAL4 sBCV     = 0.0 - inputData2[k].im; */   /* note complex conjugate */
    /* REAL4 rBCVSpin = inputData3[k].re;
       REAL4 sBCVSpin = 0.0 - inputData3[k].im;  */  /* note complex conjugate */

    REAL4 x = tmpltSignal[k].re;
    REAL4 y = tmpltSignal[k].im;     
 
 
    qtilde[k].re        = r * x - s * y ;
    qtilde[k].im        = r * y + s * x ;
    /* qtildeBCV[k].re     = rBCV * x - sBCV * y ;
       qtildeBCV[k].im     = rBCV * y + sBCV * x ; */
    /*    qtildeBCVSpin[k].re = rBCVSpin * x - sBCVSpin * y ;
	  qtildeBCVSpin[k].im = rBCVSpin * y + sBCVSpin * x ; */
    
    
    qtildeBCV[k]     = qtilde[k];
    /*qtildeBCVSpin[k] = qtilde[k]; */

    qtilde[k].re        *= a1;
    qtildeBCV[k].re     *= a2;
    /*qtildeBCVSpin[k].re *= a3; */

    qtilde[k].im        *= a1;
    qtildeBCV[k].im     *= a2;
    /*qtildeBCVSpin[k].im *= a3; */

  }

  /* qtilde negative frequency only: not DC or nyquist */
  if ( params->computeNegFreq )
  {
    for ( k = numPoints/2 + 2; k < numPoints - 1; ++k )
    {
      REAL4 r        = inputData1[k].re;
      REAL4 s        = 0.0 - inputData1[k].im;    /* note complex conjugate */
      /* REAL4 rBCV     = inputData2[k].re;
	 REAL4 sBCV     = 0.0 - inputData2[k].im; */   /* note complex conjugate */
      /* REAL4 rBCVSpin = inputData3[k].re;
	 REAL4 sBCVSpin = 0.0 - inputData3[k].im;  */  /* note complex conjugate */

      REAL4 x = tmpltSignal[k].re;
      REAL4 y = tmpltSignal[k].im;

      qtilde[k].re = r * x - s * y ;
      qtilde[k].im = r * y + s * x ;
      /*qtildeBCV[k].re = rBCV * x - sBCV * y ;
	qtildeBCV[k].im = rBCV * y + sBCV * x ;*/
      /*qtildeBCVSpin[k].re = rBCVSpin * x - sBCVSpin * y ;
	qtildeBCVSpin[k].im = rBCVSpin * y + sBCVSpin * x ; */

      qtildeBCV[k]     = qtilde[k];
      /*qtildeBCVSpin[k] = qtilde[k]; */

     
      qtilde[k].re        *= a1;
      qtildeBCV[k].re     *= a2;
      /*qtildeBCVSpin[k].re *= a3; */

      qtilde[k].im        *= a1;
      qtildeBCV[k].im     *= a2;
      /*qtildeBCVSpin[k].im *= a3; */

    }
   }
 
  

 
 



/*code*/




  DETATCHSTATUSPTR( status );
  RETURN( status );
}
