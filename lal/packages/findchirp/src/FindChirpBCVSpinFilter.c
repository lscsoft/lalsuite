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
  UINT4                 eventStartIdx  = 0;
  REAL4                 chirpTime      = 0;
  BOOLEAN               haveChisq      = 0;
  COMPLEX8             *qtilde         = NULL; 
  COMPLEX8             *qtildeBCVSpin1 = NULL; 
  COMPLEX8             *qtildeBCVSpin2 = NULL; 
  COMPLEX8             *q              = NULL; 
  COMPLEX8             *qBCVSpin1      = NULL;
  COMPLEX8             *qBCVSpin2      = NULL;
  COMPLEX8             *inputData      = NULL;
  COMPLEX8             *inputDataBCV   = NULL;
  COMPLEX8             *tmpltSignal    = NULL;
  SnglInspiralTable    *thisEvent      = NULL;
  LALMSTUnitsAndAcc     gmstUnits;
  FindChirpSegment      *fcSeg;
  DataSegment           *dataSeg; 
  REAL4                 templateNorm;
  REAL4                 modqsq;
  COMPLEX8              *wtilde;  /* need new pointer name? */
  REAL4                 *amp;
  REAL4                 *ampBCV;
  REAL4                 *ampBCVSpin1;
  REAL4                 *ampBCVSpin2;
  REAL4                 I = 0.0;
  REAL4                 J = 0.0;
  REAL4                 K = 0.0;
  REAL4                 L = 0.0;
  REAL4                 M = 0.0;
  REAL4                 Beta; /* Spin parameter, value from bank or external loop */  
  REAL4                 rootI;
  REAL4                 denominator;
  REAL4                 rootDenominator;
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
  q         = params->qVec->data;
  qBCVSpin1 = params->qVecBCVSpin1->data;
  qBCVSpin2 = params->qVecBCVSpin2->data; 

  qtilde         = params->qtildeVec->data;
  qtildeBCVSpin1 = params->qtildeVecBCVSpin1->data;
  qtildeBCVSpin2 = params->qtildeVecBCVSpin2->data; 
  
  
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




  amp         = inputParams->ampVec->data;
  ampBCV      = inputParams->ampVecBCV->data;
  ampBCVSpin1 = inputParams->ampVecBCVSpin1->data;
  ampBCVSpin2 = inputParams->ampVecBCVSpin2->data;

  
  wtilde     = inputParams->wtildeVec->data;

  for ( k = 1; k < fcSeg->data->data->length; ++k )
  {
    I += amp[k] * amp[k] * wtilde[k].re ;
    J += amp[k] * amp[k] * wtilde[k].re * 
      cos(Beta * ampBCVSpin2[k]);                
    K += amp[k] * amp[k] * wtilde[k].re * 
      sin(Beta * ampBCVSpin2[k]);
    L += amp[k] * amp[k] * wtilde[k].re * 
      sin(2 * Beta * ampBCVSpin2[k]);
    M += amp[k] * amp[k] * wtilde[k].re * 
      cos(2 * Beta * ampBCVSpin2[k]);
  }

  /* Taking multiplucation outside loop lessens cost */

  I *= 4;
  J *= 4;
  K *= 4;
  L *= 4;
  M *= 4;

  /* Expensive or well used quantities calc before loop */

 rootI           = sqrt(I);
 denominator     = I*M  +  0.5*pow(I,2) - pow(J,2);
 rootDenominator = sqrt(denominator);
 denominator1    = sqrt ( 0.25 * pow(I,3) 
	+ M*(pow(J,2) - pow(K,2)) 
	- 0.5 *(pow(J,2) + pow(K,2)) 
	- I * (pow(L,2) + pow(M,2)) 
	+ 2*J*K*L );


 
 

  /*
   * initialising outputData vectors to
   * calibrated detector output as calc in LAL..Data()
   * note lack of exponential terms, these are
   * calc in LALFindChirpBCVSpinTemplate()
   */



  inputData1 = fcSeg->data->data->data;
 
  
  /*
   *
   * compute qtilde, qtildeBCVSpin1 and qtildeBCVSpin2
   *
   * needs close checking
   *
   */


  memset( qtilde,         0, numPoints * sizeof(COMPLEX8) );
  memset( qtildeBCVSpin1, 0, numPoints * sizeof(COMPLEX8) );
  memset( qtildeBCVSpin2, 0, numPoints * sizeof(COMPLEX8) );

  /* qtilde positive frequency, not DC or nyquist */
  for ( k = 1; k < numPoints/2; ++k )
  {
    REAL4 r        = inputData1[k].re;
    REAL4 s        = 0.0 - inputData1[k].im;       /* note complex conjugate */
   
    REAL4 x = tmpltSignal[k].re;
    REAL4 y = tmpltSignal[k].im;     
 
 
    qtilde[k].re        = r * x - s * y ;
    qtilde[k].im        = r * y + s * x ;
     
    
    qtildeBCVSpin1[k] = qtilde[k];
    qtildeBCVSpin2[k] = qtilde[k]; 

    /* real parts */
    qtilde[k].re         *= amp[k]/ rootI ;                        /* A1 */

    qtildeBCVSpin1[k].re *= amp[k]/ (rootDenominator * rootI);              
    qtildeBCVSpin1[k].re *=(I * cos(Beta * ampBCVSpin2[k]) -  J);  /* A2 */

    qtildeBCVSpin2[k].re *= (amp[k]/denominator1); 
    qtildeBCVSpin2[k].re *= (sin(Beta * ampBCVSpin2[k])
	    - ( (I*L - J*K) * cos(Beta * ampBCVSpin2[k]) / denominator)
	    + ( (J*L - K*M + 0.5*I*K) / denominator ) );           /* A3 */


    /* imaginary parts */
    qtilde[k].im         *= amp[k]/ rootI ;                        /* A1 */

    qtildeBCVSpin1[k].im *= amp[k]/(rootDenominator * rootI);  
    qtildeBCVSpin1[k].im *= (I * cos(Beta * ampBCVSpin2[k]) -  J); /* A2 */
 
    qtildeBCVSpin2[k].im *= (amp[k]/denominator1); 
    qtildeBCVSpin2[k].im *= (sin(Beta * ampBCVSpin2[k])
	    - ( (I*L - J*K) * cos(Beta * ampBCVSpin2[k]) / denominator)
	    + ( (J*L - K*M + 0.5*I*K) / denominator ) );           /* A3 */
  }

  /* qtilde negative frequency only: not DC or nyquist */
  if ( params->computeNegFreq )
  {
    for ( k = numPoints/2 + 2; k < numPoints - 1; ++k )
    {
      REAL4 r        = inputData1[k].re;
      REAL4 s        = 0.0 - inputData1[k].im;    /* note complex conjugate */
     

      REAL4 x = tmpltSignal[k].re;
      REAL4 y = tmpltSignal[k].im;

      qtilde[k].re = r * x - s * y ;
      qtilde[k].im = r * y + s * x ;
     

      qtildeBCVSpin1[k] = qtilde[k];
      qtildeBCVSpin2[k] = qtilde[k]; 

      /* real parts */
      qtilde[k].re         *= amp[k]/ rootI ;                        /* A1 */

      qtildeBCVSpin1[k].re *= amp[k]/(rootDenominator * rootI);              
      qtildeBCVSpin1[k].re *=(I * cos(Beta * ampBCVSpin2[k]) -  J);  /* A2 */

      qtildeBCVSpin2[k].re *= (amp[k]/denominator1); 
      qtildeBCVSpin2[k].re *= (sin(Beta * ampBCVSpin2[k])
	    - ( (I*L - J*K) * cos(Beta * ampBCVSpin2[k]) / denominator)
	    + ( (J*L - K*M + 0.5*I*K) / denominator ) );             /* A3 */

      

      /* imaginary parts */
      qtilde[k].im         *= amp[k]/ rootI ;                        /* A1 */

      qtildeBCVSpin1[k].im *= amp[k]/(rootDenominator * rootI);  
      qtildeBCVSpin1[k].im *= (I * cos(Beta * ampBCVSpin2[k]) -  J); /* A2 */
 
      qtildeBCVSpin2[k].im *= (amp[k]/denominator1); 
      qtildeBCVSpin2[k].im *= (sin(Beta * ampBCVSpin2[k])
	    - ( (I*L - J*K) * cos(Beta * ampBCVSpin2[k]) / denominator)
	    + ( (J*L - K*M + 0.5*I*K) / denominator ) );             /* A3 */
     
    }
   }
 
   /* 
    *
    * inverse fft to get q, qBCVSpin1 and qBCVSpin2
    *    
    */

   LALCOMPLEX8VectorFFT( status->statusPtr, params->qVec, 
		   params->qtildeVec, params->invPlan );
   CHECKSTATUSPTR( status );

   LALCOMPLEX8VectorFFT( status->statusPtr, params->qVecBCVSpin1, 
		   params->qtildeVecBCVSpin1, params->invPlan );
   CHECKSTATUSPTR( status );

   LALCOMPLEX8VectorFFT( status->statusPtr, params->qVecBCVSpin2, 
		   params->qtildeVecBCVSpin2, params->invPlan );
		   CHECKSTATUSPTR( status );

   /* 
    *
    * calculate signal to noise squared 
    *
    */

   /* square and add terms calc above */

   /* limits ?, normalisation */

   /*  for ( j = 0; j < numPoints; ++j)
       { 
   rhosq[j]   = q[j].re * q[j].re + q[j].im * q[j].im; 
   rhosq[j]  += qBCVSpin1[j].re * qBCVSpin1[j].re + qBCVSpin1[j].im * qBCVSpin1[j].im; 
   rhosq[j]  += qBCVSpin2[j].re * qBCVSpin2[j].re + qBCVSpin2[j].im * qBCVSpin2[j].im; 
   } */

/*code*/




  DETATCHSTATUSPTR( status );
  RETURN( status );
}
