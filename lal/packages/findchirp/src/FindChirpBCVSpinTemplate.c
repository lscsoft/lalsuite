/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpSBCVTemplate.c
 *
 * Author: Brown D. A., Spinning BCV-Modifications: Jones, G
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpBCVSpin.h>

NRCSID (FINDCHIRPSBCVTEMPLATEC, "$Id$");

#if 0 
<lalVerbatim file="FindChirpBCVSpinTemplateCV">
Author: Brown, D. A., Spinning BCV-Modifications: Jones, G
$Id$
</lalVerbatim> 

<lalLaTeX>
\subsection{Module \texttt{FindChirpBCVSpinTemplate.c}}
\label{ss:FindChirpBCVSpinTemplate.c}

Provides functions to create spinning BCV detection templates in a form that
can be used by the \texttt{FindChirpBCVSpinFilter()} function.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpBCVSpinTemplateCP}
\idx{LALFindChirpBCVSpinTemplate()}

The function \texttt{LALFindChirpBCVSpinTemplate()} creates the 
spinning BCV template as described by the algorithm below.

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

\vfill{\footnotesize\input{FindChirpBCVSpinTemplateCV}}
</lalLaTeX> 
#endif

/* <lalVerbatim file="FindChirpBCVSpinTemplateCP"> */
void
LALFindChirpBCVSpinTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpTmpltParams       *params,
    FindChirpDataParams        *fcDataParams	
    )
/* </lalVerbatim> */
{
  UINT4        		numPoints        = 0;
  REAL4        		deltaF           = 0.0;
  COMPLEX8    	       *expPsi           = NULL;
  REAL4       	       *xfac             = NULL;
  REAL4        		x1               = 0.0;  
  REAL4        		psi0             = 0.0;
  REAL4        		psi00            = 0.0;
  REAL4        		psi05            = 0.0;
  REAL4        		psi10            = 0.0;
  REAL4        		psi15            = 0.0;
  REAL4        		psi20            = 0.0;        
  REAL4       		fFinal           = 0.0;  
  INT4      	        k                = 0;
  INT4         		kmin             = 0;     
  INT4         		kmax             = 0;    
  REAL8                *ampBCVSpin1;
  REAL8                *ampBCVSpin2;
  COMPLEX8             *wtilde;
  REAL8                 I                = 0.0;
  REAL8                 J                = 0.0;
  REAL8                 K                = 0.0;
  REAL8                 L                = 0.0;
  REAL8                 M                = 0.0;
  REAL4                 beta; 
  REAL8                 rootI;
  REAL8                 denominator;
  REAL8                 rootDenominator;
  REAL8                 denominator1;
  REAL8                 numerator1;
  REAL4                 Twoby3           = 2.0/3.0;
  REAL8                 deltaTto2by3;
  REAL8                *A1Vec            = NULL;
  REAL8                *A2Vec            = NULL;
  REAL8                *A3Vec            = NULL;
  REAL4                 deltaT;
  REAL4                 fLow;
  REAL4                 A1A1             = 0.0;
  REAL4                 A2A2             = 0.0;
  REAL4                 A3A3 		 = 0.0;
  REAL4                 A1A2 		 = 0.0;
  REAL4                 A1A3 		 = 0.0;
  REAL4                 A2A3 		 = 0.0;
  REAL4                 mTotal           = 0.0;
  REAL4                 rLSOto3by2       = 0.0;
  int		        doTest;
  
  FILE                 *fpTmpltIm  	 = NULL;
  FILE                 *fpTmpltRe  	 = NULL;
  FILE                 *fpwtilde   	 = NULL;
  FILE                 *fpA1       	 = NULL;
  FILE                 *fpA2       	 = NULL;
  FILE                 *fpA3             = NULL;
  
  FILE                 *fpTemplate       = NULL;

  REAL4                 factor;
  INT4                  i;
  RealFFTPlan          *prev             = NULL;
  REAL4Vector          *hVec             = NULL;
  REAL4Vector          *HVec             = NULL;
  REAL4                 invNumPoints;
  
  INITSTATUS( status, "LALFindChirpBCVSpinTemplate", FINDCHIRPSBCVTEMPLATEC );
  ATTATCHSTATUSPTR( status );

  /* check that the output structures exist */
  ASSERT( fcTmplt, status, 
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
  ASSERT( fcTmplt->data, status, 
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
  ASSERT( fcTmplt->data->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );

  /* check that the parameter structure exists */         
  ASSERT( params, status, FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );

  /* check that the parameter structure is set */
  /* to the correct waveform approximant       */
  ASSERT( params->approximant == BCVSpin, status, 
      FINDCHIRPBCVSPINH_EMAPX, FINDCHIRPBCVSPINH_MSGEMAPX );

  /* check that the timestep is positive */
  ASSERT( params->deltaT > 0, status,                        
      FINDCHIRPBCVSPINH_EDELT, FINDCHIRPBCVSPINH_MSGEDELT );

  /* check that the input exists */
  ASSERT( tmplt, status, FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );

  /*
   * Choose level of output
   */
  
  doTest = 0; /* Set to 1 for lots of output, useful for testing */
  
  if (doTest ==1)
  {	  
  fpTemplate       = fopen ("TemplateBCVSpin.dat","w");	  
  fpTmpltIm        = fopen ("tmpltIm.dat","w");
  fpTmpltRe        = fopen ("tmpltRe.dat","w");
  fpwtilde         = fopen ("wtilde.dat", "w");
  fpA1             = fopen ("A1.dat","w");
  fpA2             = fopen ("A2.dat","w");
  fpA3             = fopen ("A3.dat","w");
  }

  /*
   *
   * compute the BCVSpin template
   *
   */

  /* set up pointers */
  expPsi      = fcTmplt->data->data;
  numPoints   = fcTmplt->data->length;
  xfac        = params->xfacVec->data;
  xfac        = params->xfacVec->data;
  
  /* Reading in data needed to calculate moments */
  wtilde      = fcDataParams->wtildeVec->data;
  ampBCVSpin1 = fcDataParams->ampVecBCVSpin1->data;
  ampBCVSpin2 = fcDataParams->ampVecBCVSpin2->data;
  
  /* store the waveform approximant */
  fcTmplt->approximant = BCVSpin;

  /* zero output */
  memset( expPsi, 0, numPoints * sizeof(COMPLEX8) );

  /* psi coefficients */
  psi00 = tmplt->psi0;            /* BCV only uses psi0, psi15:            */
  psi05 = 0.0; /*tmplt->psi1;*/   /* -> psi1,2,4 don't exist in tmplt      */
  psi10 = 0.0; /*tmplt->psi2;*/   /* -> use if statements to define these? */
  psi15 = tmplt->psi3;            /* & which name convention to use?       */
  psi20 = 0.0; /*tmplt->psi4;*/
  /* XXX work needed here... */
  
  /* parameters */
  deltaT        = params->deltaT;
  deltaTto2by3  = pow(deltaT, Twoby3);
  deltaF        = 1.0 / ( (REAL4) params->deltaT * (REAL4) numPoints ); 
  x1            = pow( deltaF, -1.0/3.0 );
  fLow          = params->fLow;
  fFinal        = tmplt->fFinal;
  kmin = params->fLow / deltaF > 1 ? params->fLow / deltaF : 1;
  kmax = fFinal / deltaF < numPoints/2 ? fFinal / deltaF : numPoints/2;
  beta = tmplt->beta; 

  /* Preliminary BCVSpin bank does not populate fFinal */
  /* will estimate fFinal form psi0, psi3 as quick fix */
  
  if (fFinal == 0.0)
  {	  
  	fprintf (stdout, "BCVSpin bank has not populated fFinal so we shall"); 
        fprintf (stdout, "estimate its values using psi0 and psi3 \n" );
        
        rLSOto3by2 = 14.69693846; /* 6 to 3by2) */	 

	fFinal = (-psi00 * 16 * LAL_PI) / (psi15 * rLSOto3by2);
  }
 
  /* since we have redefined fFinal we must redefine kmax */
  
  kmax = fFinal / deltaF < numPoints/2 ? fFinal / deltaF : numPoints/2;
   

  
  if (doTest ==1)
  {	  
  	fprintf (stdout, "psi0      = %e \n", psi00);
  	fprintf (stdout, "psi3      = %e \n", psi15);
  	fprintf (stdout, "fFinal    = %e \n", fFinal);
  	fprintf (stdout, "fLow      = %e \n", fLow);
	fprintf (stdout, "numPoints = %d \n", numPoints);
	fprintf (stdout, "deltaT    = %e \n", deltaT);
	fprintf (stdout, "deltaF    = %e \n", deltaF);
  	fprintf (stdout, "kmin      = %d \n", kmin);
  	fprintf (stdout, "kmax      = %d \n", kmax);
  	fprintf (stdout, "beta      = %e \n", beta);
  }
 
  /* compute psi0: used in range reduction */
  {
    REAL4 x    = x1 * xfac[kmin];
    REAL4 psi  = 
    psi20 + (x * x) * ( psi15 + x * ( psi10 + x * ( psi05 + x * ( psi00 ))));
    psi0 = -2 * LAL_PI * ( floor ( 0.5 * psi / LAL_PI ) );

    if (doTest ==1)
    {	    
    	fprintf (stdout, "xfac[kmin] = %e \n", xfac[kmin]);
    	fprintf (stdout, "x = %e \n", x);
    	fprintf (stdout, "psi = %e \n", psi);
    	fprintf (stdout, "psi0 = %e \n", psi0);
    }		
  }
  /* XXX work needed here... check psi */

  /*
   *
   * calculate the stationary phase chirp
   *
   */

   /*fprintf (stdout, "just before loop in template code \n");*/

    for ( k = kmin; k < kmax ; ++k )
    {
      REAL4 x    = x1 * xfac[k];
      REAL4 psi  = 
      psi20 + (x * x) * ( psi15 + x * ( psi10 + x * ( psi05 + x * ( psi00 ))));
      REAL4 psi1 = psi + psi0;

      /* XXX work needed here... check psi */
      /* leaving psi calc here, needs to be inside template loop 
       * but not inside loop over data segments 
       */

      /* range reduction of psi1 */
      while ( psi1 < -LAL_PI )
	{
	  psi1 += 2 * LAL_PI;
	  psi0 += 2 * LAL_PI;
	}
      while ( psi1 > LAL_PI )
	{
	  psi1 -= 2 * LAL_PI;
	  psi0 -= 2 * LAL_PI;
	}

      /* compute approximate sine and cosine of psi1 */
      /* XXX The sign of this is different than the SP filtering
       * because the data is conjugated instead of the template in the
       * BCV code */
      expPsi[k].im =   -sin(psi1);
      expPsi[k].re =   cos(psi1);
        
      if (doTest ==1)
      {	      
     		fprintf (fpTmpltIm, "%d\t%e\n",k, expPsi[k].im);
      	 	fprintf (fpTmpltRe, "%d\t%e\n",k, expPsi[k].re);
      }

    }

        /* Uncomment following code to output template in time domain  */
       
   if (doTest == 1)
   {	   
        LALCreateReverseRealFFTPlan(status->statusPtr, &prev, numPoints,0);
        LALSCreateVector(status->statusPtr, &hVec, numPoints);
        LALSCreateVector(status->statusPtr, &HVec, numPoints);

       for (i=0; i<(numPoints); i++)
        {
           HVec->data[i] = 0.0;
        }

      factor = LAL_TWOPI * 0.5;


       for (i=kmin; i<kmax; i++)
        {
                HVec->data[i]           =  expPsi[i].re *  cos((REAL4)i*factor) - expPsi[i].im *  sin((REAL4)i*factor);
                HVec->data[numPoints-i] =  expPsi[i].im *  cos((REAL4)i*factor) + expPsi[i].re *  sin((REAL4)i*factor);
        }
                HVec->data[0] = 0.;
                HVec->data[numPoints/2] = 0.;

        /*for ( i = 0; i < numPoints; ++i)
        {
           fprintf (fpTemplate, "%d\t%e\n", i, HVec->data[i]);
        }*/


	fprintf (stdout, "before FFT \n");
                                                                                                                           
	LALREAL4VectorFFT( status->statusPtr, hVec, HVec, prev );
	CHECKSTATUSPTR( status );
                                                                                                                             
	fprintf (stdout, "after FFT \n");

	invNumPoints = 1.0 / ((REAL4) numPoints); 
	fprintf (stdout, "invNumPoints  %e\n", invNumPoints );

	for ( i = 0; i < numPoints; ++i)
                                                                                                                             
        {
           hVec->data[i] *= invNumPoints;
           fprintf (fpTemplate, "%d\t%e\n", i, hVec->data[i]);
        }

	fclose (fpTemplate);

	LALDestroyRealFFTPlan ( status->statusPtr, &prev);
	LALSDestroyVector ( status->statusPtr, &hVec);
	LALSDestroyVector ( status->statusPtr, &HVec);
   }

  /*
   *
   * Calculating the amplitude vectors A1Vec, A2Vec, A3Vec
   *
   */
    
  for ( k = kmin; k < kmax; ++k )
  {
          I += ampBCVSpin1[k] * ampBCVSpin1[k] * wtilde[k].re ;
          J += ampBCVSpin1[k] * ampBCVSpin1[k] * wtilde[k].re *
                  cos(beta * ampBCVSpin2[k] * deltaTto2by3);
          K += ampBCVSpin1[k] * ampBCVSpin1[k] * wtilde[k].re *
                  sin(beta * ampBCVSpin2[k] * deltaTto2by3);
          L += ampBCVSpin1[k] * ampBCVSpin1[k] * wtilde[k].re *
                  sin(2 * beta * ampBCVSpin2[k] * deltaTto2by3);
          M += ampBCVSpin1[k] * ampBCVSpin1[k] * wtilde[k].re *
                  cos(2 * beta * ampBCVSpin2[k] * deltaTto2by3);
                                                                                                                             
          if (doTest == 1)
          {
                  fprintf (fpwtilde,"%d\t%e\t%e\n",k,wtilde[k].re, wtilde[k].im);
          }
  }
                                                                                                                             
  /* Taking multiplucation outside loop lessens cost */
                                                                                                                             
  I *= 4*deltaF;
  J *= 4*deltaF;
  K *= 4*deltaF;
  L *= 2*deltaF;
  M *= 2*deltaF;
                                                                                                                         
  /* To find absolute values of these moments multiply by (deltaT)^(-7/3)  */
	  
  if (doTest == 1)
  {
          fprintf (stdout, "Moments of the noise: \n");
          fprintf (stdout, "I = %e \n", I);
          fprintf (stdout, "J = %e \n", J);
          fprintf (stdout, "K = %e \n", K);
          fprintf (stdout, "L = %e \n", L);
          fprintf (stdout, "M = %e \n\n", M);
  }
                                                                                                                             
  /* Expensive or well used quantities calc before loop */
 
  
  rootI           = sqrt(I);
  denominator     = I*M  +  0.5*pow(I,2) - pow(J,2);
  rootDenominator = sqrt(denominator);
  numerator1      = (I*L)-(J*K);
  denominator1    =  sqrt( (0.5*pow(I,2)) -(I*M) - pow(K,2)
          -  (pow(numerator1,2)/denominator) );
                           
  fcTmplt->momentI               = I;
  fcTmplt->momentJ               = J;
  fcTmplt->momentK               = K;
  
  fcTmplt->rootMomentI           = rootI;
  fcTmplt->numFactor             = denominator; 
  fcTmplt->numFactor1            = rootDenominator;
  fcTmplt->numFactor2            = numerator1; 
  fcTmplt->numFactor3            = denominator1; 
  
  if (doTest == 1)
  {
	fprintf (stdout, "rootI           = %e \n", rootI);
      	fprintf (stdout, "denominator     = %e \n", denominator);
	fprintf (stdout, "rootDenominator = %e \n", rootDenominator);
	fprintf (stdout, "denominator1    = %e \n", denominator1);
	fprintf (stdout, "numerator1      = %e \n\n", numerator1);    
  }
  
  A1Vec = fcTmplt->A1BCVSpin->data;
  A2Vec = fcTmplt->A2BCVSpin->data;
  A3Vec = fcTmplt->A3BCVSpin->data;
  
  /*LALDCreateVector(status->statusPtr, &A1Vec, (numPoints/2)+1);
  LALDCreateVector(status->statusPtr, &A2Vec, (numPoints/2)+1);
  LALDCreateVector(status->statusPtr, &A3Vec, (numPoints/2)+1);*/
  
  /* IMPROVE THIS */
  memset( A1Vec, 0, ((numPoints/2)+1) * sizeof(REAL4) );
  memset( A2Vec, 0, ((numPoints/2)+1) * sizeof(REAL4) );
  memset( A3Vec, 0, ((numPoints/2)+1) * sizeof(REAL4) );
  
  A1Vec[0] = 0;  
  A2Vec[0] = 0;
  A3Vec[0] = 0;
  
  if (beta == 0.0)
  {
  	/* fprintf(stdout,"Inside beta = 0 loop \n"); */
	                                                                                                                 
        /*for ( k = 1; k < ((numPoints/2)+1); ++k )*/
        for ( k = kmin; k < kmax; ++k )
	{
		A1Vec[k] = ampBCVSpin1[k] / rootI;
		A2Vec[k] = 0.0;
		A3Vec[k] = 0.0;
		                                                        
		if (doTest == 1)
		{
			fprintf (fpA1, "%d\t%e\n", k, A1Vec[k]);
			fprintf (fpA2, "%d\t%e\n", k, A2Vec[k]);
			fprintf (fpA3, "%d\t%e\n", k, A3Vec[k]);
	        }
	 }
  }
  else
    {
  	/* fprintf(stdout,"Inside beta not = 0 loop \n"); */

        /*for ( k = 1; k < ((numPoints/2)+1); ++k )*/
 	for ( k = kmin; k < kmax; ++k )
  	{
    		A1Vec[k] = ampBCVSpin1[k] / rootI;
    		A2Vec[k] = ampBCVSpin1[k] 
                        * (   ( cos(beta * ampBCVSpin2[k] * deltaTto2by3) )    
			-  (J/I) ) * rootI
     	                / rootDenominator ;
    		A3Vec[k] = (ampBCVSpin1[k]/denominator1) * 
        	        ( sin(beta * ampBCVSpin2[k]  * deltaTto2by3) 
                	- (K/I)
                        - (numerator1 * ( cos(beta * ampBCVSpin2[k] 
                        * deltaTto2by3) - (J/I) )/denominator )  ) 
                        * rootI;          		
	
		if (doTest ==1)
		{
 	 		fprintf (fpA1, "%d\t%e\n", k, A1Vec[k]);
 	 		fprintf (fpA2, "%d\t%e\n", k, A2Vec[k]);
         		fprintf (fpA3, "%d\t%e\n", k, A3Vec[k]);
		}
  	 }
  }  

  /* checking orthonormalisation of A vectors */

  if (doTest == 1)
  {	
	fprintf (stdout, "Checking orthonormalisation of amplitude vectors \n");
	  
  	/*for (k=0; k < (numPoints/2) + 1; ++k)*/
  	for (k=kmin; k < kmax; ++k)
  	{ 
  		A1A1 += A1Vec[k] * A1Vec[k] * wtilde[k].re;
		A2A2 += A2Vec[k] * A2Vec[k] * wtilde[k].re;
        	A3A3 += A3Vec[k] * A3Vec[k] * wtilde[k].re;
        	A1A2 += A1Vec[k] * A2Vec[k] * wtilde[k].re;  
        	A1A3 += A1Vec[k] * A3Vec[k] * wtilde[k].re;       
        	A2A3 += A2Vec[k] * A3Vec[k] * wtilde[k].re;  
  	}

  	A1A1 *= 4 * deltaF;
  	A2A2 *= 4 * deltaF;
  	A3A3 *= 4 * deltaF;
  	A1A2 *= 4 * deltaF;
  	A1A3 *= 4 * deltaF;
  	A2A3 *= 4 * deltaF;

  	fprintf (stdout, "A1hat cross A1hat %e\n", A1A1);  
  	fprintf (stdout, "A2hat cross A2hat %e\n", A2A2);
  	fprintf (stdout, "A3hat cross A3hat %e\n", A3A3);
  	fprintf (stdout, "A1hat cross A2hat %e\n", A1A2);
  	fprintf (stdout, "A1hat cross A3hat %e\n", A1A3);
  	fprintf (stdout, "A2hat cross A3hat %e\n\n", A2A3);
  }


  if (doTest ==1)
  {
        fclose (fpTmpltIm);
        fclose (fpTmpltRe);
        fclose (fpwtilde);
        fclose (fpA1);
        fclose (fpA2);
        fclose (fpA3);
  }


  DETATCHSTATUSPTR( status );
  RETURN( status );
}


