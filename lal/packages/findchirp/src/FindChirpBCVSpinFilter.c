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

Provides functions to filter data for spinning BCV templates.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpBCVSpinFilterCP}
\idx{LALFindChirpBCVSpinFilter()}

The function \texttt{LALFindChirpBCVSpinFilter()} filters data for 
spinning BCV templates as described by the algorithm below.

\subsubsection*{Algorithm}

Blah.

\subsubsection*{Uses}
\begin{verbatim}
LALCalloc()
LALFree()
LALCreateVector()
LALDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

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
#include <lal/FindChirpBCVSpin.h>


NRCSID (FINDCHIRPBCVSPINFILTERC, "$Id$");

/* <lalVerbatim file="FindChirpBCVSpinFilterCP"> */
void
LALFindChirpBCVSpinFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params,             
    FindChirpDataParams        *fcDataParams,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec
  )
/* </lalVerbatim> */
{
  UINT4                 i, j, k;
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

  REAL4			rhosq;

  REAL4 		rho;
  REAL4                 invMaxRho;

  REAL4			maxRho;
  UINT4			maxRhoCount;

  REAL4   		alpha1;
  REAL4                 alpha2;
  REAL4                 alpha3;
  REAL4                 alpha4;
  REAL4                 alpha5;
  REAL4                 alpha6;
  REAL4			wvff;
  REAL4                 wvft;              /*params used to produce BCVSpin waveform */	
  
  FILE     		*fpRho =  NULL;
  FILE                  *fpWvtIm =  NULL;
  FILE                  *fpWvtRe =  NULL;
  FILE			*fpWvfIm =  NULL;
  FILE			*fpWvfRe =  NULL;
  FILE                  *fpwtilde = NULL;
  FILE                  *fpf23      = NULL;
 

  RealFFTPlan           *prev = NULL;
  REAL4Vector           *hVec = NULL;
  REAL4Vector           *HVec = NULL;

  INITSTATUS( status, "LALFindChirpBCVSpinFilter", FINDCHIRPBCVSPINFILTERC );
  ATTATCHSTATUSPTR( status );


  /*
   *    
   * check that the arguments are reasonable
   * may need to remove asserts regarding chisq
   *          
   */

 fprintf (stdout, "before asserts in FindChirpBCVSpinFilter \n");

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
 /* ASSERT(params->qVecBCV, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qVecBCV->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVecBCV, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVecBCV->data,status, FINDCHIRPH_ENULL, 
	  FINDCHIRPH_MSGENULL);*/
  
  /* check that the chisq parameter and input structures exist */
  ASSERT( params->chisqParams, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->chisqInput,   status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  /*ASSERT( params->chisqInputBCV,status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );*/

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

  /* make sure that the template and the segment are both BCVSpin */
  ASSERT( input->fcTmplt->approximant == BCVSpin, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  ASSERT( input->segment->approximant == BCVSpin, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );

/*fprintf (stdout, "after asserts in FindChirpBCVSpinFilter \n");*/

/*
*
*check sense of this later
*
*/


fpRho   = fopen ("rho.dat", "w");
fpWvtIm = fopen ("wvtIm.dat", "w");
fpWvtRe = fopen ("wvtRe.dat", "w");
fpWvfIm = fopen ("wvfIm.dat", "w");
fpWvfRe = fopen ("wvfRe.dat", "w");
fpwtilde = fopen ("wtilde.dat", "w");
fpf23     = fopen("f23.dat", "w");


for (i = 0; i < dataSegVec->length; ++i)
{
fcSeg = & (fcSegVec->data[i]);
}
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
  
  numPoints = params->qVec->length; 
 
  /* template and data */
  inputData     = input->segment->data->data->data;
 /* inputDataBCV  = input->segment->dataBCV->data->data;*/
  tmpltSignal   = input->fcTmplt->data->data;
  templateNorm  = input->fcTmplt->tmpltNorm;   /* this is expPsi */
  deltaT        = params->deltaT;


/*
 *
 * code which prev. would be in data function
 *
 */




  amp         = fcDataParams->ampVec->data;
 /* ampBCV      = fcDataParams->ampVecBCV->data;*/
  ampBCVSpin1 = fcDataParams->ampVecBCVSpin1->data;
  ampBCVSpin2 = fcDataParams->ampVecBCVSpin2->data;

/*fprintf (stdout, "before moments calc. in FindChirpBCVSpinFilter \n");*/

  
  wtilde     = fcDataParams->wtildeVec->data;

/*fprintf (stdout, "before moments calc. in FindChirpBCVSpinFilter1 \n");*/

fprintf (stdout, "fcSeg->data->data->length: %d \n", fcSeg->data->data->length);

Beta = 0.1;

fprintf (stdout, "Beta: %e \n", Beta);

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

fprintf (fpwtilde, "%d\t%e\n", k, wtilde[k].re);
fprintf (fpf23, "%d\t%e\n", k, ampBCVSpin2[k]);

  }

  /* Taking multiplucation outside loop lessens cost */

  I *= 4;
  J *= 4;
  K *= 4;
  L *= 2;
  M *= 2;

fprintf (stdout, "I = %e \n", I); 
fprintf (stdout, "J = %e \n", J);
fprintf (stdout, "K = %e \n", K);
fprintf (stdout, "L = %e \n", L);
fprintf (stdout, "M = %e \n", M);

/*fprintf (stdout, "after moments calc. in FindChirpBCVSpinFilter \n");*/


  /* Expensive or well used quantities calc before loop */

 rootI           = sqrt(I);
 denominator     = I*M  +  0.5*pow(I,2) - pow(J,2);
 rootDenominator = sqrt(denominator);
 denominator1    = sqrt ( 0.25 * pow(I,3) 
	+ M*(pow(J,2) - pow(K,2)) 
	- 0.5 * I * (pow(J,2) + pow(K,2)) 
	- I * (pow(L,2) + pow(M,2)) 
	+ 2*J*K*L );

/*fprintf (stdout, "rootI           = %e \n", rootI);
fprintf (stdout, "denominator     = %e \n", denominator);
fprintf (stdout, "rootDenominator = %e \n", rootDenominator);
fprintf (stdout, "denominator1    = %e \n", denominator1);*/


 
 

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


/*fprintf (stdout, "before qtilde calc. in FindChirpBCVSpinFilter \n");*/


  memset( qtilde,         0, numPoints * sizeof(COMPLEX8) );
  memset( qtildeBCVSpin1, 0, numPoints * sizeof(COMPLEX8) );
  memset( qtildeBCVSpin2, 0, numPoints * sizeof(COMPLEX8) );

/*fprintf (stdout, "numPoints = %d \n", numPoints);*/


  /* qtilde positive frequency, not DC or nyquist */
  for ( k = 1; k < numPoints/2; ++k )
  {
    REAL4 r        = inputData1[k].re;
    REAL4 s        = 0.0 - inputData1[k].im;       /* note complex conjugate */
   
    REAL4 x = tmpltSignal[k].re;
    REAL4 y = tmpltSignal[k].im;     
 
 
    qtilde[k].re        = r * x - s * y ;
    qtilde[k].im        = r * y + s * x ;

/* is  this correct? */
                                                                                                                             
      qtilde[k].re *= wtilde[k].re;
      qtilde[k].im *= wtilde[k].re;

/*fprintf (stdout, "loop         = %d \n", k);
fprintf (stdout, "x            = %e \n", x);
fprintf (stdout, "y            = %e \n", y);

fprintf (stdout, "qtilde[k].re = %e \n", qtilde[k].re);
fprintf (stdout, "qtilde[k].im = %e \n", qtilde[k].im);*/
    

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

/*fprintf (stdout, "qtilde[k].re         = %e \n", qtilde[k].re);
fprintf (stdout, "qtildeBCVSpin1[k].re = %e \n", qtildeBCVSpin1[k].re);
fprintf (stdout, "qtildeBCVSpin2[k].re = %e \n", qtildeBCVSpin2[k].re);*/




    /* imaginary parts */
    qtilde[k].im         *= amp[k]/ rootI ;                        /* A1 */

    qtildeBCVSpin1[k].im *= amp[k]/(rootDenominator * rootI);  
    qtildeBCVSpin1[k].im *= (I * cos(Beta * ampBCVSpin2[k]) -  J); /* A2 */
 
    qtildeBCVSpin2[k].im *= (amp[k]/denominator1); 
    qtildeBCVSpin2[k].im *= (sin(Beta * ampBCVSpin2[k])
	    - ( (I*L - J*K) * cos(Beta * ampBCVSpin2[k]) / denominator)
	    + ( (J*L - K*M + 0.5*I*K) / denominator ) );           /* A3 */

/*fprintf (stdout, "qtilde[k].im         = %e \n", qtilde[k].im);
fprintf (stdout, "qtildeBCVSpin1[k].im = %e \n", qtildeBCVSpin1[k].im);
fprintf (stdout, "qtildeBCVSpin2[k].im = %e \n", qtildeBCVSpin2[k].im);*/


  }

/*fprintf (stdout, "just after +ve  freq loop \n");*/

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
    
/* is  this correct? */

      qtilde[k].re *= wtilde[k].re;
      qtilde[k].im *= wtilde[k].re;
 

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
 
/*fprintf (stdout, "after qtilde calc. in FindChirpBCVSpinFilter \n");*/


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


/*fprintf (stdout, "after FFT calc. in FindChirpBCVSpinFilter \n");

fprintf (stdout, "length of q            = %d \n", params->qVec->length);
fprintf (stdout, "length of qVecBCVSpin1 = %d \n", params->qVecBCVSpin1->length);
fprintf (stdout, "length of qVecBCVSpin2 = %d \n", params->qVecBCVSpin2->length);*/

/* for ( j = 0; j < numPoints; ++j)
       {
fprintf (stdout, "q[j]            = %e \n", q[j]);
fprintf (stdout, "qBCVSpin1[j]    = %e \n", qBCVSpin1[j]);
fprintf (stdout, "qBCVSpin2[j]    = %e \n", qBCVSpin2[j]);

fprintf (stdout, "q[j].re            = %e \n", q[j].re);
fprintf (stdout, "qBCVSpin1[j].re    = %e \n", qBCVSpin1[j].re);
fprintf (stdout, "qBCVSpin2[j].re    = %e \n", qBCVSpin2[j].re);

fprintf (stdout, "q[j].im            = %e \n", q[j].im);
fprintf (stdout, "qBCVSpin1[j].im    = %e \n", qBCVSpin1[j].im);
fprintf (stdout, "qBCVSpin2[j].im    = %e \n", qBCVSpin2[j].im);
       }*/

   /* 
    *
    * calculate signal to noise squared 
    *
    */

   /* square and add terms calc above */

   /* limits ?, normalisation */

memset (params->rhosqVec->data->data, 0, numPoints * sizeof (REAL4) );

maxRho = 0;
maxRhoCount = 0;

   for ( j = 0; j < numPoints; ++j)
       {


   rhosq   = q[j].re * q[j].re + q[j].im * q[j].im; 
   rhosq  += qBCVSpin1[j].re * qBCVSpin1[j].re + qBCVSpin1[j].im * qBCVSpin1[j].im; 
   rhosq  += qBCVSpin2[j].re * qBCVSpin2[j].re + qBCVSpin2[j].im * qBCVSpin2[j].im;
 
   /* Inverse FFT normalisation, division by numPoints squared */

 /*  rhosq  = rhosq / pow(numPoints,2);*/

/*fprintf (stdout, "rhosq            = %e \n", rhosq);*/

/* to account for 2's  in   eqn 22,23 */
    
  rhosq *= 4;

   params->rhosqVec->data->data[j] = rhosq;  /*copied from Eirini's code */

/*fprintf (stdout, "rhosqVec            = %e \n",  params->rhosqVec->data->data[j] );*/

   rho    = pow(rhosq, 0.5);
   
fprintf (fpRho, "%d\t%e\n", j, rho);

   /* finding max value of rho in time */
if (rho > maxRho)
	{
	maxRho = rho;
        maxRhoCount = j;
	}

/*fprintf (stdout, "maxRho       = %e \n", maxRho);
fprintf (stdout, "maxRhoCount  = %d \n", maxRhoCount);*/
 

       } 

 invMaxRho = 1 / maxRho;
                                                                                                                             
/* calc alpha values in time domain */
                                                                                                                             
   alpha1 = q[maxRhoCount].re * invMaxRho;
   alpha2 = q[maxRhoCount].im * invMaxRho;
   alpha3 = qBCVSpin1[maxRhoCount].re * invMaxRho;
   alpha4 = qBCVSpin1[maxRhoCount].im * invMaxRho;
   alpha5 = qBCVSpin2[maxRhoCount].re * invMaxRho;
   alpha6 = qBCVSpin2[maxRhoCount].im * invMaxRho;

fprintf (stdout, "alpha1       = %e \n", alpha1);
fprintf (stdout, "alpha2       = %e \n", alpha2);
fprintf (stdout, "alpha3       = %e \n", alpha3);
fprintf (stdout, "alpha4       = %e \n", alpha4);
fprintf (stdout, "alpha5       = %e \n", alpha5);
fprintf (stdout, "alpha6       = %e \n", alpha6);



                                                                                                                  
 /*calc freq domain waveform, store in qtilde[k]  */

  for ( k = 0; k < (numPoints/2)+1; ++k )
  {                                                                                                                          
    REAL4 x = tmpltSignal[k].re;
    REAL4 y = tmpltSignal[k].im;
                                                                                                                             
/* check here */                                                                                                                             
    qtilde[k].re         = (alpha1 * x - alpha2 * y) * amp[k]/ rootI ;
    qtilde[k].re        += (alpha3 * x - alpha4 * y) 
				*  amp[k]/(rootDenominator * rootI)
				* (I * cos(Beta * ampBCVSpin2[k]) -  J) ;
    qtilde[k].re        += (alpha5 * x - alpha6 * y) 
				* (amp[k]/denominator1)
				*  (sin(Beta * ampBCVSpin2[k])
            - ( (I*L - J*K) * cos(Beta * ampBCVSpin2[k]) / denominator)
            + ( (J*L - K*M + 0.5*I*K) / denominator ) ) ;     


    qtilde[k].im         = (alpha1 * y + alpha2 * x) * amp[k]/ rootI ;
    qtilde[k].im        += (alpha3 * y + alpha4 * x)
                                *  amp[k]/(rootDenominator * rootI)
                                * (I * cos(Beta * ampBCVSpin2[k]) -  J) ;
    qtilde[k].im        += (alpha5 * y + alpha6 * x)
                                * (amp[k]/denominator1)
                                *  (sin(Beta * ampBCVSpin2[k])
            - ( (I*L - J*K) * cos(Beta * ampBCVSpin2[k]) / denominator)
            + ( (J*L - K*M + 0.5*I*K) / denominator ) ) ;
    }
  
fprintf (stdout, "numPoints         = %d \n", numPoints);



/* inverse FFT to find time domain waveform, store in hVec */


        LALCreateReverseRealFFTPlan(status->statusPtr, &prev, numPoints,0);
        LALSCreateVector(status->statusPtr, &hVec, numPoints);
        LALSCreateVector(status->statusPtr, &HVec, numPoints);
/*
COMPLEX8Vector *HVec = NULL;
*/

        for (i=1; i<=(numPoints/2); i++)
        {
                HVec->data[i] = qtilde[i].re;
                HVec->data[numPoints-i] = -qtilde[i].im;  /* - to correct time reversal */
	}       
                HVec->data[0] = 0.;
                HVec->data[numPoints/2] = 0.;

 	
 for (i=1; i<(numPoints/2)+1; i++)
        {
	 fprintf (fpWvfIm, "%d %e %e\n", i, HVec->data[i], HVec->data[numPoints-i]);
	}



/*
LALReverseRealFFT( status->statusPtr, hVec,
                   HVec, prev );
   CHECKSTATUSPTR( status );
*/

fprintf (stdout, "after FFT \n");

LALREAL4VectorFFT( status->statusPtr, hVec, HVec, prev );
CHECKSTATUSPTR( status );

/*
LALCOMPLEX8VectorFFT( status->statusPtr, params->qVec,
                   params->qtildeVec, params->invPlan );
   CHECKSTATUSPTR( status );*/

for ( j = 0; j < numPoints; ++j)
       
	{
/*	fprintf (stdout, "q[j]            = %e \n", q[j]);
	fprintf (stdout, "q[j].re         = %e \n", q[j].re);
	fprintf (stdout, "q[j].im         = %e \n", q[j].im);*/

/*	fprintf (fpWvtIm, "%d\t%e\n", j, q[j].im);*/
/*        fprintf (fpWvtRe, "%d\t%e\n", j, q[j].re);*/
 fprintf (fpWvtRe, "%d\t%e\n", j, hVec->data[j]);


       }

fclose (fpRho);
fclose (fpWvfIm);
fclose (fpWvfRe);
fclose (fpWvtIm);
fclose (fpWvtRe);
fclose (fpwtilde);
fclose (fpf23);

LALSDestroyVector ( status->statusPtr, &hVec);
LALSDestroyVector ( status->statusPtr, &HVec);

LALDestroyRealFFTPlan ( status->statusPtr, &prev);

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
