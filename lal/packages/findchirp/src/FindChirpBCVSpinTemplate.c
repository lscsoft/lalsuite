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
    FindChirpTmpltParams       *params
    )
/* </lalVerbatim> */
{
  UINT4        numPoints  = 0;
  REAL4        deltaF     = 0.0;
  COMPLEX8    *expPsi     = NULL;
  REAL4       *xfac       = NULL;
  REAL4        x1         = 0.0;  
  REAL4        psi0       = 0.0;
  REAL4        psi00      = 0.0;
  REAL4        psi05      = 0.0;
  REAL4        psi10      = 0.0;
  REAL4        psi15      = 0.0;
  REAL4        psi20      = 0.0;        
  REAL4        fHi        = 0.0;  
  INT4         k          = 0;
  /*REAL4        factor;*/
  INT4         kmin       = 0;     
  INT4         kmax       = 0;    
 
  FILE        *fpTmpltIm      = NULL;
  FILE        *fpTmpltRe      = NULL;
  /* FILE        *fpTemplate     = NULL; */

  /*RealFFTPlan *prev           = NULL;
  REAL4Vector *hVec           = NULL;
  REAL4Vector *HVec           = NULL;
  REAL4        invNumPoints;*/

  int		doTest;
  
  INITSTATUS( status, "LALFindChirpBCVSpinTemplate", FINDCHIRPSBCVTEMPLATEC );
  ATTATCHSTATUSPTR( status );

 /*
   *
   * check that the arguments are reasonable
   * same as Eirini's except for aproximant
   *
   */
  

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
  fpTmpltIm  =  fopen("tmpltIm.dat","w");
  fpTmpltRe  =  fopen("tmpltRe.dat","w");
  /* fpTemplate =  fopen("Template.dat","w"); */
  }

  /*
   *
   * compute the BCVSpin template
   *
   */


  /* set up pointers */
  expPsi    = fcTmplt->data->data;
  xfac      = params->xfacVec->data;
  numPoints = fcTmplt->data->length;
  
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
  deltaF = 1.0 / ( (REAL4) params->deltaT * (REAL4) numPoints ); 

  x1 = pow( deltaF, -1.0/3.0 );
  fHi  = tmplt->fFinal;
  
  kmin = params->fLow / deltaF > 1 ? params->fLow / deltaF : 1;
  kmax = fHi / deltaF < numPoints/2 ? fHi / deltaF : numPoints/2;

  if (doTest ==1)
  {	  
  	fprintf (stdout, "psi0    = %e \n", psi00);
  	fprintf (stdout, "psi3    = %e \n", psi15);
  	fprintf (stdout, "fHi     = %e \n", fHi);
  	fprintf (stdout, "fLow    = %e \n", params->fLow);
  	fprintf (stdout, "deltaF  = %e \n", deltaF);
  	fprintf (stdout, "kmin    = %d \n", kmin);
  	fprintf (stdout, "kmax    = %d \n", kmax);
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
      expPsi[k].im =   sin(psi1);
      expPsi[k].re =   cos(psi1);
        
      if (doTest ==1)
      {	      
     		fprintf (fpTmpltIm, "%d\t%e\n",k, expPsi[k].im);
      	 	fprintf (fpTmpltRe, "%d\t%e\n",k, expPsi[k].re);
      }

    }

if (doTest ==1)
{	
	fclose (fpTmpltIm);
	fclose (fpTmpltRe);
}
        /* Uncomment following code to output template in time domain  */
        /*
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
                HVec->data[numPoints-i] = expPsi[i].im *  cos((REAL4)i*factor) - expPsi[i].re *  sin((REAL4)i*factor);*/
        /* - to correct time reversal */
        /*}
                HVec->data[0] = 0.;
                HVec->data[numPoints/2] = 0.;*/

        /*for ( i = 0; i < numPoints; ++i)
        {
           fprintf (fpTemplate, "%d\t%e\n", i, HVec->data[i]);
        }*/


/*fprintf (stdout, "before FFT \n");*/
/*                                                                                                                           
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
*/  

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


