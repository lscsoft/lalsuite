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
  REAL4                 deltaF;
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
  REAL8                 *ampBCVSpin1;
  REAL8                 *ampBCVSpin2;
  REAL8                 I = 0.0;
  REAL8                 J = 0.0;
  REAL8                 K = 0.0;
  REAL8                 L = 0.0;
  REAL8                 M = 0.0;
  REAL4                 Beta; /* Spin parameter, value from bank or external loop */  
  REAL8                 rootI;
  REAL8                 denominator;
  REAL8                 rootDenominator;
  REAL8                 denominator1;
  REAL8                 numerator1;
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
  REAL4                 invRho;

  REAL4			maxRho;
  UINT4			maxRhoCount;

  REAL4                 normFac;

  REAL4   		alpha1;
  REAL4                 alpha2;
  REAL4                 alpha3;
  REAL4                 alpha4;
  REAL4                 alpha5;
  REAL4                 alpha6;
  REAL4                 alphaSumSq;
  REAL4			wvff;
  REAL4                 wvft;              /*params used to produce BCVSpin waveform */	
  
  REAL4                 Sevenby6 = 7.0/6.0;
  REAL4			Twoby3   = 2.0/3.0;

  REAL8                 deltaTto7by6;
  REAL8		        deltaTto2by3;
     

  FILE     		*fpRho =  NULL;
  FILE                  *fpWvtIm =  NULL;
  FILE                  *fpRecon =  NULL;
  FILE			*fpWvfIm =  NULL;
  FILE			*fpWvfRe =  NULL;
  FILE                  *fpwtilde = NULL;
  FILE                  *fpf23      = NULL;
  FILE                  *fpamp =  NULL;
  FILE                  *fpampBCV =  NULL;
  FILE                  *fpampBCVSpin1 =  NULL;
  FILE                  *fpampBCVSpin2 =  NULL;
  FILE                  *fpA1;
  FILE                  *fpA2;
  FILE                  *fpA3;
  FILE                  *fpalphaSumSq;
  
  FILE                  *fpqtilde = NULL;
  FILE                  *fpqtildeBCVSpin1 = NULL;
  FILE                  *fpqtildeBCVSpin2 = NULL;

  FILE                  *fpq = NULL;
  FILE                  *fpqBCVSpin1 = NULL;
  FILE                  *fpqBCVSpin2 = NULL;
 
  FILE                  *fpStrain1Re = NULL;
  FILE                  *fpStrain1Im = NULL;

  FILE                  *fpDataNorm;

  RealFFTPlan           *prev = NULL;
  REAL4Vector           *hVec = NULL;
  REAL4Vector           *HVec = NULL;
  /*
  RealFFTPlan           *prev1 = NULL;
  REAL4Vector           *qNewVec = NULL;
  REAL4Vector           *qtildeNewVec = NULL;
  REAL4Vector           *qNewVecBCVSpin1 = NULL;
  REAL4Vector           *qtildeNewVecBCVSpin1 = NULL;
  REAL4Vector           *qNewVecBCVSpin2 = NULL;
  REAL4Vector           *qtildeNewVecBCVSpin2 = NULL;
  */

  RealFFTPlan           *prev1 = NULL;
  REAL4Vector           *dVec = NULL;
  REAL4Vector           *DVec = NULL;


  REAL8Vector           *A1Vec = NULL;
  REAL8Vector           *A2Vec = NULL;
  REAL8Vector           *A3Vec = NULL;
 
  REAL4                 factor;
  REAL4                 T0; 

  REAL4                 normData;
  REAL4                 invRootNormData;

  REAL4                 invNumPoints;

  REAL4                 A1A1;
  REAL4                 A2A2;
  REAL4                 A3A3;
  REAL4                 A1A2;
  REAL4                 A1A3;
  REAL4                 A2A3;


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
fpRecon = fopen ("Recon.dat", "w");
fpWvfIm = fopen ("wvfIm.dat", "w");
fpWvfRe = fopen ("wvfRe.dat", "w");
fpwtilde = fopen ("wtilde.dat", "w");
fpf23     = fopen("f23.dat", "w");
/*fpamp = fopen("amp.dat","w");
fpampBCV = fopen("ampBCV.dat","w");*/
fpampBCVSpin1 = fopen("ampBCVSpin1.dat","w");
fpampBCVSpin2 = fopen("ampBCVSpin2.dat","w");
fpA1 = fopen ("A1.dat","w");
fpA2 = fopen ("A2.dat","w");
fpA3 = fopen ("A3.dat","w");
fpalphaSumSq = fopen("alphaSumSq.dat","w");

fpqtilde = fopen ("qtilde.dat","w");
fpqtildeBCVSpin1 = fopen ("qtildeBCVSpin1.dat","w");
fpqtildeBCVSpin2 = fopen ("qtildeBCVSpin2.dat","w");

fpq = fopen ("q.dat","w");
fpqBCVSpin1 = fopen ("qBCVSpin1.dat","w");
fpqBCVSpin2 = fopen ("qBCVSpin2.dat","w");

fpStrain1Re = fopen("Strain1Re.dat","w");
fpStrain1Im = fopen("Strain1Im.dat","w");

fpDataNorm = fopen("DataNorm.dat","w");


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
fprintf (stdout, "Filter  code, numPoints = %d\n", numPoints); 
  /* template and data */
  inputData     = input->segment->data->data->data;
 /* inputDataBCV  = input->segment->dataBCV->data->data;*/
  tmpltSignal   = input->fcTmplt->data->data;  /* this is expPsi */
  templateNorm  = input->fcTmplt->tmpltNorm;   
  deltaT        = params->deltaT;

  normFac       = 4./numPoints;


/*
 *
 * code which prev. would be in data function
 *
 */




  amp         = fcDataParams->ampVec->data;
 /* ampBCV      = fcDataParams->ampVecBCV->data;*/
  ampBCVSpin1 = fcDataParams->ampVecBCVSpin1->data;
  ampBCVSpin2 = fcDataParams->ampVecBCVSpin2->data;

fprintf (stdout, "before moments calc. in FindChirpBCVSpinFilter \n");

  
  wtilde     = fcDataParams->wtildeVec->data;

fprintf (stdout, "before moments calc. in FindChirpBCVSpinFilter1 \n");

fprintf (stdout, "fcSeg->data->data->length: %d \n", fcSeg->data->data->length);

Beta = 100.0;

fprintf (stdout, "Beta: %e \n", Beta);
fprintf (stdout, "deltaT: %e \n", deltaT);

/* for ( k = 1; k < fcSeg->data->data->length; ++k )
  {
    fprintf  (stdout, "%d\t%e\n", k,amp[k]);
  }*/

fprintf( stdout, "dt 7/6 %e\n",  pow(deltaT, Sevenby6  )); 

/* deltaTto7by6 = pow(deltaT, Sevenby6);*/
deltaTto2by3 = pow(deltaT, Twoby3);

/* we now calc. true value for f^(-2/3) from the ampBCVSpin2 vector, we need to
calc. the true value since f^(-2/3) is used inside sin and cos fns. */

 for ( k = 1; k < fcSeg->data->data->length; ++k )
  {
 /*   ampBCVSpin1[k] *= deltaTto7by6; */  /* converting amp to f */  
    ampBCVSpin2[k] *= deltaTto2by3;       
  }

 for ( k = 1; k < fcSeg->data->data->length; ++k )
  {
	fprintf(fpampBCVSpin1, "%d\t%e\n", k, ampBCVSpin1[k]);
        fprintf(fpampBCVSpin2, "%d\t%e\n", k, ampBCVSpin2[k]);
  }


fprintf (stdout, "after amp calcs");


  for ( k = 1; k < fcSeg->data->data->length; ++k )
  {
    I += ampBCVSpin1[k] * ampBCVSpin1[k] * wtilde[k].re ;
    J += ampBCVSpin1[k] * ampBCVSpin1[k] * wtilde[k].re * 
      cos(Beta * ampBCVSpin2[k]);                
    K += ampBCVSpin1[k] * ampBCVSpin1[k] * wtilde[k].re * 
      sin(Beta * ampBCVSpin2[k]);
    L += ampBCVSpin1[k] * ampBCVSpin1[k] * wtilde[k].re * 
      sin(2 * Beta * ampBCVSpin2[k]);
    M += ampBCVSpin1[k] * ampBCVSpin1[k] * wtilde[k].re * 
      cos(2 * Beta * ampBCVSpin2[k]);

/*fprintf (stdout, "k=  %d\t I= %e\n ", k, I);*/
fprintf (fpwtilde, "%d\t%e\t%e\n", k, wtilde[k].re, wtilde[k].im);
fprintf (fpf23, "%d\t%e\n", k, ampBCVSpin2[k]);

  }

 /* Taking multiplucation outside loop lessens cost */

deltaF  =  1.0/((REAL4)numPoints * deltaT);
fprintf(stdout, "deltaF %e\n", deltaF);

  I *= 4*deltaF;
  J *= 4*deltaF;
  K *= 4*deltaF;
  L *= 2*deltaF;
  M *= 2*deltaF;

/* !!! To find absolute values of these moments multiply by (deltaT)^(-7/3) !!! */ 


fprintf (stdout, "I = %e \n", I); 
fprintf (stdout, "J = %e \n", J);
fprintf (stdout, "K = %e \n", K);
fprintf (stdout, "L = %e \n", L);
fprintf (stdout, "M = %e \n", M);

fprintf (stdout, "after moments calc. in FindChirpBCVSpinFilter \n");


  /* Expensive or well used quantities calc before loop */

 rootI           = sqrt(I);
 denominator     = I*M  +  0.5*pow(I,2) - pow(J,2);
 rootDenominator = sqrt(denominator);
 /*denominator1    = sqrt ( (0.25 * pow(I,3)) 
	+ (M*(pow(J,2) - pow(K,2))) 
	- (0.5 * I * (pow(J,2) + pow(K,2))) 
	- (I * (pow(L,2) + pow(M,2))) 
	+ (2*J*K*L) );*/
numerator1     = (I*L)-(J*K);
denominator1   =  sqrt( (0.5*pow(I,2)) -(I*M) - pow(K,2) -  (pow(numerator1,2)/denominator) );


fprintf (stdout, "rootI           = %e \n", rootI);
fprintf (stdout, "denominator     = %e \n", denominator);
fprintf (stdout, "rootDenominator = %e \n", rootDenominator);
fprintf (stdout, "denominator1    = %e \n", denominator1);
fprintf (stdout, "numerator1    = %e \n", numerator1);

  
LALDCreateVector(status->statusPtr, &A1Vec, (numPoints/2)+1);
LALDCreateVector(status->statusPtr, &A2Vec, (numPoints/2)+1);
LALDCreateVector(status->statusPtr, &A3Vec, (numPoints/2)+1);

A1Vec->data[0] = 0;
A2Vec->data[0] = 0;
A3Vec->data[0] = 0;

if (Beta == 0.0)
{
 fprintf(stdout,"inside beta = 0 loop \n");

 for ( k = 1; k < fcSeg->data->data->length; ++k )
	{
 	 A1Vec->data[k] = ampBCVSpin1[k] / rootI;
         A2Vec->data[k] = 0.0;
	 A3Vec->data[k] = 0.0;

         fprintf (fpA1, "%d\t%e\n", k, A1Vec->data[k]);
         fprintf (fpA2, "%d\t%e\n", k, A2Vec->data[k]);
         fprintf (fpA3, "%d\t%e\n", k, A3Vec->data[k]);
       }
}
else
{
 fprintf(stdout,"inside beta not = 0 loop \n");

 for ( k = 1; k < fcSeg->data->data->length; ++k )
  	{
    	A1Vec->data[k] = ampBCVSpin1[k] / rootI;
    	A2Vec->data[k] = ampBCVSpin1[k] * (   ( cos(Beta * ampBCVSpin2[k]) )    -  (J/I) ) * rootI
     	                    / rootDenominator ;
    	A3Vec->data[k] = (ampBCVSpin1[k]/denominator1) * 
        	                    ( sin(Beta * ampBCVSpin2[k]) 
                	              - (K/I)
                        	      - (numerator1 * ( cos(Beta * ampBCVSpin2[k]) - (J/I) )/denominator )  
                        	    ) * rootI;          

 	 fprintf (fpA1, "%d\t%e\n", k, A1Vec->data[k]);
 	 fprintf (fpA2, "%d\t%e\n", k, A2Vec->data[k]);
         fprintf (fpA3, "%d\t%e\n", k, A3Vec->data[k]);
  	 }
}  

/* checking orthonormalisation of A vectors */

A1A1 = 0.;
A2A2 = 0.;
A3A3 = 0.;
A1A2 = 0.;
A1A3 = 0.;
A2A3 = 0.;

for (k=0; k < (numPoints/2) + 1; ++k)
	{ 
         A1A1 += A1Vec->data[k] * A1Vec->data[k] * wtilde[k].re;
	 A2A2 += A2Vec->data[k] * A2Vec->data[k] * wtilde[k].re;
         A3A3 += A3Vec->data[k] * A3Vec->data[k] * wtilde[k].re;
         A1A2 += A1Vec->data[k] * A2Vec->data[k] * wtilde[k].re;  
         A1A3 += A1Vec->data[k] * A3Vec->data[k] * wtilde[k].re;       
         A2A3 += A2Vec->data[k] * A3Vec->data[k] * wtilde[k].re;  
       }

A1A1 *= 4 * deltaF;
A2A2 *= 4 * deltaF;
A3A3 *= 4 * deltaF;
A1A2 *= 4 * deltaF;
A1A3 *= 4 * deltaF;
A2A3 *= 4 * deltaF;

fprintf (stdout, "A1 cross A1 %e\n", A1A1);  
fprintf (stdout, "A2 cross A2 %e\n", A2A2);
fprintf (stdout, "A3 cross A3 %e\n", A3A3);
fprintf (stdout, "A1 cross A2 %e\n", A1A2);
fprintf (stdout, "A1 cross A3 %e\n", A1A3);
fprintf (stdout, "A2 cross A3 %e\n", A2A3);


  /*
   * initialising outputData vectors to
   * calibrated detector output as calc in LAL..Data()
   * note lack of exponential terms, these are
   * calc in LALFindChirpBCVSpinTemplate()
   */



  inputData1 = fcSeg->data->data->data;
 
  for (k = 0; k < (numPoints/2)+1; ++k )
{
  fprintf (fpStrain1Re, "%d\t%e\n", k, inputData1[k].re);
  fprintf (fpStrain1Im, "%d\t%e\n", k, inputData1[k].im);
} 

 
  /*
   *
   * compute qtilde, qtildeBCVSpin1 and qtildeBCVSpin2
   *
   * needs close checking
   *
   */


/* finding cross produc of data with itself,  to be used for normalisation later  */

normData = 0.;

for (k = 0; k < (numPoints/2)+1; ++k )
{
normData += ((inputData1[k].re * inputData1[k].re) + (inputData1[k].im * inputData1[k].im))
             * wtilde[k].re;
}

 normData *= deltaT * normFac;   
/* normData *= 2.0/numPoints; MAYBE CHANGE BACK TO THIS!!!*/

fprintf (stdout, "input data normData = %e\n", normData);

invRootNormData = pow(normData,-0.5);
/*
for (k = 0; k < (numPoints/2)+1; ++k )
{
inputData1[k].re *= invRootNormData;
inputData1[k].im *= invRootNormData;  
}*/
  /*  UNCOMMENT LOOP TO NORMALISE INPUT  */

        LALCreateReverseRealFFTPlan(status->statusPtr, &prev1, numPoints,0);
        LALSCreateVector(status->statusPtr, &dVec, numPoints);
        LALSCreateVector(status->statusPtr, &DVec, numPoints);
                                                                                                                             
        for (i=1; i<(numPoints/2); i++)
                                                                                                                             
        {
		DVec->data[i] = inputData1[i].re;
                DVec->data[numPoints-i] = inputData1[i].im;
        }
                DVec->data[0] = 0.;
                DVec->data[numPoints/2] = 0.;
                                                                                                                                                                                                                                                         
fprintf (stdout, "before FFT \n");
                                                                                                                             
LALREAL4VectorFFT( status->statusPtr, dVec, DVec, prev1 );
CHECKSTATUSPTR( status );

for (i=0; i < numPoints; i++)
	{       
                dVec->data[i] *= 1./numPoints; 
		fprintf (fpDataNorm,  "%e\n", dVec->data[i]);
	}

normData = 0.;
                                                                                                                             
for (k = 0; k < (numPoints/2)+1; ++k )
{
normData += ((inputData1[k].re * inputData1[k].re) + (inputData1[k].im * inputData1[k].im))
             * wtilde[k].re;
}
                                                                                                                             
normData *= deltaT * normFac;
                                                                                                                             
fprintf (stdout, "input data normData = %e\n", normData);

/*

A1norm = 0.;
                                                                                                                             
for (k = 0; k < (numPoints/2)+1; ++k )
{
A1norm += ((tmpltSignal[k].re * tmpltSignal[k].re) + (tmpltSignal[k].im * tmpltSignal[k].im))
             * wtilde[k].re * A1Vec->data[k] * A1Vec->data[k] ;
}
                                                                                                                             
A1norm *= 2./numPoints;
                                                                                                                             
fprintf (stdout, "A1norm = %e\n", A1norm);


invRootNormData = pow(normData,-0.5);
                                                                                                                             
for (k = 0; k < (numPoints/2)+1; ++k )
{
tmpltSignal[k].re *= invRootNormData;
tmpltSignal[k].im *= invRootNormData;
}
                                                                                                                             
normData = 0.;
                                                                                                                             
for (k = 0; k < (numPoints/2)+1; ++k )
{
normData += ((tmpltSignal[k].re * tmpltSignal[k].re) + (tmpltSignal[k].im * tmpltSignal[k].im))
             * wtilde[k].re;
}
                                                                                                                             
normData *= 4 * deltaF;
                                                                                                                             
fprintf (stdout, "template normData = %e\n", normData);


*/


normData = 0.;
                                                                                                                             
for (k = 0; k < (numPoints/2)+1; ++k )
{
normData += wtilde[k].re * ( pow(tmpltSignal[k].re,2) + pow(tmpltSignal[k].im,2) ) ;
}
                                                                                                                             
normData *= 4 * deltaF;
                                                                                                                             
fprintf (stdout, " 4 int df over snf   = %e\n", normData);


/*
invRootNormData = pow(normData,-0.5);
                                                                                                                             
for (k = 0; k < (numPoints/2)+1; ++k )
{
tmpltSignal[k].re *= invRootNormData;
tmpltSignal[k].im *= invRootNormData;
}
*/                                                                                                                             
                                                                                                                             







/* IFFT inputdata  back to time domain, check abs values */

/* LALCOMPLEX8VectorFFT( status->statusPtr, input,
                   , params->invPlan );
   CHECKSTATUSPTR( status );*/

fprintf (stdout, "before qtilde calc. in FindChirpBCVSpinFilter \n");


  memset( qtilde,         0, numPoints * sizeof(COMPLEX8) );
  memset( qtildeBCVSpin1, 0, numPoints * sizeof(COMPLEX8) );
  memset( qtildeBCVSpin2, 0, numPoints * sizeof(COMPLEX8) );

fprintf (stdout, "numPoints = %d \n", numPoints);


  /* qtilde positive frequency, not DC or nyquist */
  for ( k = 1; k < numPoints/2; ++k )
  {
    REAL4 r        =   inputData1[k].re; 
    REAL4 s        =   inputData1[k].im;      
   
    REAL4 x =  tmpltSignal[k].re;
    REAL4 y =  0. - tmpltSignal[k].im;     
 
    qtilde[k].re        = r * x + s * y ;
    qtilde[k].im        = s * x - r * y ;
   


   /* qtilde[k].re        = r * x - s * y ;
    qtilde[k].im        = r * y + s * x ;*/

/* fprintf(stdout, "%e \t %e  \n ", qtilde[k].re, qtilde[k].im);*/


/* is  this correct? */
                                                                                                                             
      qtilde[k].re *= wtilde[k].re;
      qtilde[k].im *= wtilde[k].re; 

/*fprintf (stdout, "loop         = %d \n", k);
fprintf (stdout, "x            = %e \n", x);
fprintf (stdout, "y            = %e \n", y);

fprintf (stdout, "qtilde[k].re = %e \n", qtilde[k].re);
fprintf (stdout, "qtilde[k].im = %e \n", qtilde[k].im);*/
    

    qtildeBCVSpin1[k] =  qtilde[k]; 
    qtildeBCVSpin2[k] =  qtilde[k];   

    /* real parts */

     qtilde[k].re         *= A1Vec->data[k];
     qtildeBCVSpin1[k].re *= A2Vec->data[k];
     qtildeBCVSpin2[k].re *= A3Vec->data[k];   /*   UNCOMMENT THIS!!!! */


 /*  qtilde[k].re         *= amp[k]/ rootI ;                        /* A1 */

 /*   qtildeBCVSpin1[k].re *= amp[k]/ (rootDenominator * rootI);              
    qtildeBCVSpin1[k].re *=(I * cos(Beta * ampBCVSpin2[k]) -  J);  /* A2 */

 /*   qtildeBCVSpin2[k].re *= (amp[k]/denominator1); 
    qtildeBCVSpin2[k].re *= (sin(Beta * ampBCVSpin2[k])
	    - ( (I*L - J*K) * cos(Beta * ampBCVSpin2[k]) / denominator)
	    + ( (J*L - K*M + 0.5*I*K) / denominator ) );           /* A3 */   

/*fprintf (stdout, "qtilde[k].re         = %e \n", qtilde[k].re);
fprintf (stdout, "qtildeBCVSpin1[k].re = %e \n", qtildeBCVSpin1[k].re);
fprintf (stdout, "qtildeBCVSpin2[k].re = %e \n", qtildeBCVSpin2[k].re);*/




    /* imaginary parts */

     qtilde[k].im         *= A1Vec->data[k];
     qtildeBCVSpin1[k].im *= A2Vec->data[k];
     qtildeBCVSpin2[k].im *= A3Vec->data[k];



/*    qtilde[k].im         *= amp[k]/ rootI ;                        /* A1 */

/*    qtildeBCVSpin1[k].im *= amp[k]/(rootDenominator * rootI);  
    qtildeBCVSpin1[k].im *= (I * cos(Beta * ampBCVSpin2[k]) -  J); /* A2 */
 
/*    qtildeBCVSpin2[k].im *= (amp[k]/denominator1); 
    qtildeBCVSpin2[k].im *= (sin(Beta * ampBCVSpin2[k])
	    - ( (I*L - J*K) * cos(Beta * ampBCVSpin2[k]) / denominator)
	    + ( (J*L - K*M + 0.5*I*K) / denominator ) );           /* A3 */

fprintf (fpqtilde,         "%d\t%e\t%e\n", k, qtilde[k].re,         qtilde[k].im);
fprintf (fpqtildeBCVSpin1, "%d\t%e\t%e\n", k, qtildeBCVSpin1[k].re, qtildeBCVSpin1[k].im);
fprintf (fpqtildeBCVSpin2, "%d\t%e\t%e\n", k, qtildeBCVSpin2[k].re, qtildeBCVSpin2[k].im);


  }

fprintf (stdout, "just after +ve  freq loop \n");

  /* qtilde negative frequency only: not DC or nyquist */
  if ( params->computeNegFreq )
  {
    for ( k = numPoints/2 + 2; k < numPoints - 1; ++k )
    {
      REAL4 r        = inputData1[k].re;
      REAL4 s        = inputData1[k].im;    /* note complex conjugate */
     

      REAL4 x = tmpltSignal[k].re;
      REAL4 y =  - tmpltSignal[k].im;

      qtilde[k].re = r * x - s * y ;
      qtilde[k].im = r * y + s * x ;
    
/* is  this correct? */

      qtilde[k].re *= wtilde[k].re;
      qtilde[k].im *= wtilde[k].re;
 

      qtildeBCVSpin1[k] = qtilde[k];
      qtildeBCVSpin2[k] = qtilde[k]; 

/*  LALSCreateVector(status->statusPtr, &A1, numPoints);
  LALSCreateVector(status->statusPtr, &A2, numPoints);
  LALSCreateVector(status->statusPtr, &A3, numPoints);*/




      /* real parts */

     qtilde[k].re         *= A1Vec->data[k];
     qtildeBCVSpin1[k].re *= A2Vec->data[k];
     qtildeBCVSpin2[k].re *= A3Vec->data[k];



/*      qtilde[k].re         *= amp[k]/ rootI ;                        /* A1 */

 /*     qtildeBCVSpin1[k].re *= amp[k]/(rootDenominator * rootI);              
      qtildeBCVSpin1[k].re *=(I * cos(Beta * ampBCVSpin2[k]) -  J);  /* A2 */

 /*     qtildeBCVSpin2[k].re *= (amp[k]/denominator1); 
      qtildeBCVSpin2[k].re *= (sin(Beta * ampBCVSpin2[k])
	    - ( (I*L - J*K) * cos(Beta * ampBCVSpin2[k]) / denominator)
	    + ( (J*L - K*M + 0.5*I*K) / denominator ) );             /* A3 */

      

      /* imaginary parts */
   
     qtilde[k].im         *= A1Vec->data[k];
     qtildeBCVSpin1[k].im *= A2Vec->data[k];
     qtildeBCVSpin2[k].im *= A3Vec->data[k];



/*   qtilde[k].im         *= amp[k]/ rootI ;                        /* A1 */

/*      qtildeBCVSpin1[k].im *= amp[k]/(rootDenominator * rootI);  
      qtildeBCVSpin1[k].im *= (I * cos(Beta * ampBCVSpin2[k]) -  J); /* A2 */
 
/*      qtildeBCVSpin2[k].im *= (amp[k]/denominator1); 
      qtildeBCVSpin2[k].im *= (sin(Beta * ampBCVSpin2[k])
	    - ( (I*L - J*K) * cos(Beta * ampBCVSpin2[k]) / denominator)
	    + ( (J*L - K*M + 0.5*I*K) / denominator ) );             /* A3 */
     
    }
   }
 


/*fprintf (stdout, "after qtilde calc. in FindChirpBCVSpinFilter \n");*/
 for ( k = 0; k < numPoints; ++k )
{
fprintf (fpWvtIm, "%d\t%e\n", k, qtilde[k].re);
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

fprintf (stdout, "overlap  %e\t%e\n", q[0].re,q[0].im);

/*
        LALCreateReverseRealFFTPlan(status->statusPtr, &prev1, numPoints,0);
        LALSCreateVector(status->statusPtr, &qNewVec, numPoints);
        LALSCreateVector(status->statusPtr, &qtildeNewVec, numPoints);
        LALSCreateVector(status->statusPtr, &qNewVecBCVSpin1, numPoints);
        LALSCreateVector(status->statusPtr, &qtildeNewVecBCVSpin1, numPoints);
        LALSCreateVector(status->statusPtr, &qNewVecBCVSpin2, numPoints);
        LALSCreateVector(status->statusPtr, &qtildeNewVecBCVSpin2, numPoints);

    for (i=0; i<(numPoints); i++)
        {
                qtildeNewVec->data[i] = 0.;
                qtildeNewVecBCVSpin1->data[i] = 0.;
                qtildeNewVecBCVSpin1->data[i] = 0.;
        }


    for (i=1; i<(numPoints/2); i++)
        {
                qtildeNewVec->data[i] = qtilde[i].re;
                qtildeNewVec->data[numPoints-i] = qtilde[i].im;

                qtildeNewVecBCVSpin1->data[i] = qtildeBCVSpin1[i].re;
                qtildeNewVecBCVSpin1->data[numPoints-i] = qtildeBCVSpin1[i].im;

                qtildeNewVecBCVSpin2->data[i] = qtildeBCVSpin2[i].re;
                qtildeNewVecBCVSpin2->data[numPoints-i] = qtildeBCVSpin2[i].im;


        }
                qtildeNewVec->data[0]           = 0.;
                qtildeNewVec->data[numPoints/2] = 0.;
               
                qtildeNewVecBCVSpin1->data[0]           = 0.;
                qtildeNewVecBCVSpin1->data[numPoints/2] = 0.;

                qtildeNewVecBCVSpin2->data[0]           = 0.;
                qtildeNewVecBCVSpin2->data[numPoints/2] = 0.;

fprintf (stdout, "before FFT \n");
                                                                                                                             
LALREAL4VectorFFT( status->statusPtr, qNewVec, qtildeNewVec, prev1 );
CHECKSTATUSPTR( status );

LALREAL4VectorFFT( status->statusPtr, qNewVecBCVSpin1, qtildeNewVecBCVSpin1, prev1 );
CHECKSTATUSPTR( status );

LALREAL4VectorFFT( status->statusPtr, qNewVecBCVSpin2, qtildeNewVecBCVSpin2, prev1 );
CHECKSTATUSPTR( status );

fprintf (stdout, "after FFT \n"); */




/*fprintf (stdout, "after FFT calc. in FindChirpBCVSpinFilter \n");

fprintf (stdout, "length of q            = %d \n", params->qVec->length);
fprintf (stdout, "length of qVecBCVSpin1 = %d \n", params->qVecBCVSpin1->length);
fprintf (stdout, "length of qVecBCVSpin2 = %d \n", params->qVecBCVSpin2->length);*/

for ( j = 0; j < numPoints; ++j)
       {
fprintf (fpq,         "%d\t%e\t%e\n", j, q[j].re,         q[j].im);
fprintf (fpqBCVSpin1, "%d\t%e\t%e\n", j, qBCVSpin1[j].re, qBCVSpin1[j].im);
fprintf (fpqBCVSpin2, "%d\t%e\t%e\n", j, qBCVSpin2[j].re, qBCVSpin2[j].im);
       }

   /* 
    *
    * calculate signal to noise squared 
    *
    */
    
memset (params->rhosqVec->data->data, 0, numPoints * sizeof (REAL4) );

maxRho = 0;
maxRhoCount = 0;
                                                                                                                           
for ( j = 0; j < numPoints; ++j)
	{
	q[j].re *= normFac;  /*4 * invNumPoints; */    /*  2 * invNumPoints;*/  
        q[j].im *= normFac; /*- 4 * invNumPoints; */    /* -2 * invNumPoints;*/
        qBCVSpin1[j].re *= normFac;/*  4 * invNumPoints; */  /*   2 * invNumPoints;*/
        qBCVSpin1[j].im *= normFac;/* -4 * invNumPoints; */ /*  -2 * invNumPoints;*/
        qBCVSpin2[j].re *= normFac;  /*  4 * invNumPoints;*/  /*   2 * invNumPoints;*/
       	qBCVSpin2[j].im *= normFac;/* -4 * invNumPoints;*/  /*  -2 * invNumPoints;*/
   /*    }
   for ( j = 0; j < numPoints; ++j)
       {     */
   	rhosq   = 0.;
   	rhosq   = q[j].re * q[j].re + q[j].im * q[j].im; 
   	rhosq  += qBCVSpin1[j].re * qBCVSpin1[j].re + qBCVSpin1[j].im * qBCVSpin1[j].im; 
   	rhosq  += qBCVSpin2[j].re * qBCVSpin2[j].re + qBCVSpin2[j].im * qBCVSpin2[j].im;

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

   	invRho = 1/rho;

	/* choice of 2 different naming conventions for alpha's */
	/* alpha1 = q[j].re * invRho;
   	alpha2 = q[j].im * invRho;
   	alpha3 = qBCVSpin1[j].re * invRho;
   	alpha4 = qBCVSpin1[j].im * invRho;
   	alpha5 = qBCVSpin2[j].re * invRho;
   	alpha6 = qBCVSpin2[j].im * invRho;*/

   	alpha1 = q[j].re * invRho;
   	alpha4 = q[j].im * invRho;
  	alpha2 = qBCVSpin1[j].re * invRho;
   	alpha5 = qBCVSpin1[j].im * invRho;
   	alpha3 = qBCVSpin2[j].re * invRho;
   	alpha6 = qBCVSpin2[j].im * invRho;

        /* alphaSumSq calc for checking purposes, should = 1 */
        alphaSumSq  =  pow(alpha1,2) + pow(alpha2,2) + pow(alpha3,2)
                    +  pow(alpha4,2) + pow(alpha5,2) + pow(alpha6,2);

	fprintf (fpalphaSumSq, "%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", j, alphaSumSq,
                                 alpha1, alpha2, alpha3, alpha4, alpha5, alpha6);

       } 

       /*
	*
	* Reconstructing waveform corresponding to maximum rho event
	*
	*/

invMaxRho = 1 / maxRho;                                                    
fprintf (stdout, "maxRho       = %e \n", maxRho);
fprintf (stdout, "maxRhoCount  = %d \n", maxRhoCount);
                                                                     
/* calc alpha values in time domain */
                                      
/* choice of 2 diff naming conventions for alpha's */                                                                        /*   alpha1 = q[maxRhoCount].re * invMaxRho;
alpha2 = q[maxRhoCount].im * invMaxRho;
alpha3 = qBCVSpin1[maxRhoCount].re * invMaxRho;
alpha4 = qBCVSpin1[maxRhoCount].im * invMaxRho;
alpha5 = qBCVSpin2[maxRhoCount].re * invMaxRho;
alpha6 = qBCVSpin2[maxRhoCount].im * invMaxRho;*/

alpha1 = q[maxRhoCount].re * invMaxRho;
alpha4 = q[maxRhoCount].im * invMaxRho;
alpha2 = qBCVSpin1[maxRhoCount].re * invMaxRho;
alpha5 = qBCVSpin1[maxRhoCount].im * invMaxRho;
alpha3 = qBCVSpin2[maxRhoCount].re * invMaxRho;
alpha6 = qBCVSpin2[maxRhoCount].im * invMaxRho;

fprintf (stdout, "alpha1       = %e \n", alpha1);
fprintf (stdout, "alpha2       = %e \n", alpha2);
fprintf (stdout, "alpha3       = %e \n", alpha3);
fprintf (stdout, "alpha4       = %e \n", alpha4);
fprintf (stdout, "alpha5       = %e \n", alpha5);
fprintf (stdout, "alpha6       = %e \n", alpha6);

alphaSumSq  =  pow(alpha1,2) + pow(alpha2,2) + pow(alpha3,2) 
            +  pow(alpha4,2) + pow(alpha5,2) + pow(alpha6,2);

fprintf (stdout, "alphaSumSq       = %e \n", alphaSumSq);

/*calc freq domain waveform, store in qtilde[k]  */

/* setting DC frequencies to 0 */
qtilde[0].re = 0;
qtilde[0].im = 0;

/* setting NYQUIST frequencies to 0 */
qtilde[numPoints/2].re = 0;
qtilde[numPoints/2].im = 0;

/* creates time domain offset for reconstructed waveform */ 
factor = LAL_TWOPI * 0.5;
fprintf(stdout, "factor=%e\n", factor);

for ( k = 1; k < (numPoints/2); ++k )
	{                                                                                                                   
	REAL4 x = tmpltSignal[k].re; 
        REAL4 y = - tmpltSignal[k].im;
  
	/* corresponds to other alpha naming convention */          
        /*  qtilde[k].re         = (alpha1 * x - alpha2 * y) * A1Vec->data[k] ;
        qtilde[k].re        += (alpha3 * x - alpha4 * y) * A2Vec->data[k] ;
  	qtilde[k].re        += (alpha5 * x - alpha6 * y) * A3Vec->data[k] ;     

    	qtilde[k].im         = (alpha1 * y + alpha2 * x) * A1Vec->data[k] ;
        qtilde[k].im        += (alpha3 * y + alpha4 * x) * A2Vec->data[k] ;
        qtilde[k].im        += (alpha5 * y + alpha6 * x) * A3Vec->data[k] ;*/
   
	qtilde[k].re         = (alpha1 * x - alpha4 * y) * A1Vec->data[k] ;
    	qtilde[k].re        += (alpha2 * x - alpha5 * y) * A2Vec->data[k] ;
        qtilde[k].re        += (alpha3 * x - alpha6 * y) * A3Vec->data[k] ;
                                                                                                                             
        qtilde[k].im         = (alpha1 * y + alpha4 * x) * A1Vec->data[k] ;
	qtilde[k].im        += (alpha2 * y + alpha5 * x) * A2Vec->data[k] ;
	qtilde[k].im        += (alpha3 * y + alpha6 * x) * A3Vec->data[k] ;
	}

/* calculate overlap of reconstructed waveform with itself, should  = 1 */  
normData =  0.;

for (k = 0; k < (numPoints/2)+1; ++k )
	{
	normData += ((qtilde[k].re * qtilde[k].re) + (qtilde[k].im * qtilde[k].im))
                   * wtilde[k].re;
	}
                                                                                                                             
normData *= 4 * deltaF;
                                                                                                                             
fprintf (stdout, "normData = %e\n", normData);

/* inverse FFT to find time domain reconstructed waveform, store in hVec */

LALCreateReverseRealFFTPlan(status->statusPtr, &prev, numPoints,0);
LALSCreateVector(status->statusPtr, &hVec, numPoints);
LALSCreateVector(status->statusPtr, &HVec, numPoints);

/* setting DC and NYQUIST parts to 0 */
HVec->data[0] = 0.;
HVec->data[numPoints/2] = 0.;

for (i=1; i<(numPoints/2); i++)
	{
	HVec->data[i] = qtilde[i].re *  cos((REAL4)i*factor) - qtilde[i].im *  sin((REAL4)i*factor);
	HVec->data[numPoints-i] = qtilde[i].im *  cos((REAL4)i*factor) + qtilde[i].re *  sin((REAL4)i*factor);  
	/* no - sign needed to correct time reversal */
	}       
                
/*for (i=1; i<(numPoints/2)+1; i++)
        {
	fprintf (fpWvfIm, "%d %e %e\n", i, HVec->data[i], HVec->data[numPoints-i]);
	}*/

fprintf (stdout, "before FFT \n");

LALREAL4VectorFFT( status->statusPtr, hVec, HVec, prev );
CHECKSTATUSPTR( status );

fprintf (stdout, "after FFT \n");

for ( j = 0; j < numPoints; ++j)
	{
/*	fprintf (stdout, "q[j]            = %e \n", q[j]);
	fprintf (stdout, "q[j].re         = %e \n", q[j].re);
	fprintf (stdout, "q[j].im         = %e \n", q[j].im);
	fprintf (fpWvtIm, "%d\t%e\n", j, q[j].im);
        fprintf (fpWvtRe, "%d\t%e\n", j, q[j].re);*/
 
        hVec->data[j] *= deltaF;         /* trying to normalise  reconstructed wvf */

        /* hVec->data[j] *= maxRho; */   /* sets Recon wvf equal in ampl. to input wave */

	fprintf (fpRecon, "%d\t%e\n", j, hVec->data[j]);


       }

fprintf (stdout, "before closing output files \n");


fclose (fpRho);
fclose (fpWvfIm);
fclose (fpWvfRe);
fclose (fpWvtIm);
fclose (fpRecon);
fclose (fpwtilde);
fclose (fpf23);

/*  fclose (fpamp);
fclose (fpampBCV); */
fclose (fpampBCVSpin1);
fclose (fpampBCVSpin2);
fclose (fpA1);
fclose (fpA2);
fclose (fpA3);
fclose (fpalphaSumSq);

fclose (fpqtilde);
fclose (fpqtildeBCVSpin1);
fclose (fpqtildeBCVSpin2);

fclose (fpq);
fclose (fpqBCVSpin1);
fclose (fpqBCVSpin2);

fclose (fpStrain1Re);
fclose (fpStrain1Im);

fclose (fpDataNorm);

LALDDestroyVector ( status->statusPtr, &A1Vec);
LALDDestroyVector ( status->statusPtr, &A2Vec);
LALDDestroyVector ( status->statusPtr, &A3Vec);
LALSDestroyVector ( status->statusPtr, &hVec);
LALSDestroyVector ( status->statusPtr, &HVec);
LALDestroyRealFFTPlan ( status->statusPtr, &prev);
LALSDestroyVector ( status->statusPtr, &dVec);
LALSDestroyVector ( status->statusPtr, &DVec);
LALDestroyRealFFTPlan ( status->statusPtr, &prev1);

fprintf (stdout, "just before end  of FindChirpBCVSpinFilter \n");

DETATCHSTATUSPTR( status );
RETURN( status );
}
