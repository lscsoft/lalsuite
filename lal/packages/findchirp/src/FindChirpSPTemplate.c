/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpSPTemplate.c
 *
 * Author: Brown D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0 
<lalVerbatim file="FindChirpSPTemplateCV">
Author: Brown, D. A.
$Id$
</lalVerbatim> 

<lalLaTeX>
\subsection{Module \texttt{FindChirpSPTemplate.c}}
\label{ss:FindChirpSPTemplat.c}

Provides functions to create stationary phase inspiral templates in a
form that can be used by the \texttt{FindChirpFilter()} function.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpSPTemplateCP}
\idx{LALFindChirpSPTemplateInit()}
\idx{LALFindChirpSPTemplateFinalize()}
\idx{LALFindChirpSPTemplate()}

The function \texttt{LALFindChirpSPTemplateInit()} takes as input the address
of a structure of type \texttt{FindChirpInitParams} containing the correct
values to intialize a search. It creates a structure of type
\texttt{FindChirpSPTmpltParams} as described above and returns its address.

The function \texttt{LALFindChirpSPTemplateFinalize()} takes as the address
of a structure of type \texttt{FindChirpSPTmpltParams} destroys this 
structure and sets the address to NULL.

The function \texttt{LALFindChirpSPTemplate()} creates the stationary phase
template as described by the algorithm below.

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

\vfill{\footnotesize\input{FindChirpSPDataCV}}
</lalLaTeX> 
#endif

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>


NRCSID (FINDCHIRPSPTEMPLATEC, "$Id$");

/* <lalVerbatim file="FindChirpSPTemplateCP"> */
void
LALFindChirpSPTemplateInit (
    LALStatus                  *status,
    FindChirpSPTmpltParams    **output,
    FindChirpInitParams        *params
    )
/* </lalVerbatim> */
{
  UINT4                         k;
  FindChirpSPTmpltParams       *outputPtr;
  REAL4                        *xfac = NULL;
  const REAL4                   exponent = -1.0/3.0;

  INITSTATUS( status, "LALFindChirpSPTemplateInit", FINDCHIRPSPTEMPLATEC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */

  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( output, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );  
  ASSERT( !*output, status, FINDCHIRPSPH_ENNUL, FINDCHIRPSPH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  
  /* make sure that the number of points is positive */
  ASSERT( params->numPoints > 0, status, 
      FINDCHIRPSPH_ENUMZ, FINDCHIRPSPH_MSGENUMZ );


  /*
   *
   * create tmplt generation parameters structure
   *
   */


  /* create the output structure */
  outputPtr = *output = (FindChirpSPTmpltParams *)
    LALCalloc( 1, sizeof(FindChirpSPTmpltParams) );
  if ( ! outputPtr )
  {
    ABORT( status, FINDCHIRPSPH_EALOC, FINDCHIRPSPH_MSGEALOC );
  }

  /* create the vector to store x^(-7/6) */
  LALCreateVector( status->statusPtr, &(outputPtr->xfacVec), 
      params->numPoints/2 + 1 );
  BEGINFAIL( status )
  {
    LALFree( outputPtr );
    *output = NULL;
  }
  ENDFAIL( status );

  xfac = outputPtr->xfacVec->data;
  memset( xfac, 0, outputPtr->xfacVec->length * sizeof(REAL4) );

  for (k = 1; k < outputPtr->xfacVec->length; ++k) 
    xfac[k] = pow( (REAL4) k, exponent );
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN (status);
}



/* <lalVerbatim file="FindChirpSPTemplateCP"> */
void
LALFindChirpSPTemplateFinalize (
    LALStatus                  *status,
    FindChirpSPTmpltParams    **output
    )
/* </lalVerbatim> */
{
  FindChirpSPTmpltParams       *outputPtr;

  INITSTATUS( status, "LALFindChirpSPTemplateFinalize", FINDCHIRPSPTEMPLATEC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* make sure handle is non-null and points to a non-null pointer */
  ASSERT( output, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( *output, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  

  /*
   *
   * destroy tmplt generation parameters structure
   *
   */


  /* local pointer to output */
  outputPtr = *output;

  /* destroy the vector of x^(-7/6) */
  LALDestroyVector( status->statusPtr, &(outputPtr->xfacVec) );
  CHECKSTATUSPTR( status );

  /* free the structure */
  LALFree( outputPtr );
  *output = NULL;


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}



/* <lalVerbatim file="FindChirpSPTemplateCP"> */
void
LALFindChirpSPTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpSPTmpltParams     *params
    )
/* </lalVerbatim> */
{
  UINT4         numPoints  = 0;
  REAL4         deltaF     = 0.0;
  REAL4         m          = 0.0;
  REAL4         eta        = 0.0;
  REAL4         mu         = 0.0;
  COMPLEX8     *expPsi     = NULL;
  REAL4        *xfac       = NULL;
  REAL4         x1         = 0.0;
  REAL4         psi0       = 0.0;
  REAL4         fHi        = 0.0;
  INT4          k          = 0;
  INT4          kmin       = 0;
  INT4          kmax       = 0;

  REAL4         distNorm;
  const REAL4   cannonDist = 1.0; /* Mpc */

  /* pn constants */
  REAL4 c0;
  REAL4 c10;
  REAL4 c15;
  REAL4 c20;

  /* taylor coefficents */
  const REAL4 s2 = -0.16605;
  const REAL4 s4 =  0.00761;
  const REAL4 c2 = -0.49670;
  const REAL4 c4 =  0.03705;
  
  
  INITSTATUS( status, "LALFindChirpSPTemplate", FINDCHIRPSPTEMPLATEC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */
  

  /* check that the output structures exist */
  ASSERT( fcTmplt, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( fcTmplt->data, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( fcTmplt->data->data, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  /* check that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  /* check that the timestep is positive */
  ASSERT( params->deltaT > 0, status, 
      FINDCHIRPSPH_EDELT, FINDCHIRPSPH_MSGEDELT );

  /* check that the input exists */
  ASSERT( tmplt, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );


  /*
   *
   * compute the stationary phase template
   *
   */


  /* set up pointers */
  expPsi = fcTmplt->data->data;
  xfac = params->xfacVec->data;
  numPoints = fcTmplt->data->length;

  /* zero output */
  memset( expPsi, 0, numPoints * sizeof(COMPLEX8) );

  /* parameters */
  deltaF = 1.0 / ( (REAL4) params->deltaT * (REAL4) numPoints );
  m      = (REAL4) tmplt->totalMass;
  eta    = (REAL4) tmplt->eta;
  mu     = (REAL4) tmplt->mu;

  /* template dependent normalisation */
  distNorm = 2.0 * LAL_MRSUN_SI / (cannonDist * 1.0e6 * LAL_PC_SI);
  distNorm *= params->dynRange;

  fcTmplt->tmpltNorm = sqrt( (5.0*mu) / 96.0 ) *
    pow( m / (LAL_PI*LAL_PI) , 1.0/3.0 ) *
    pow( LAL_MTSUN_SI / (REAL4) params->deltaT, -1.0/6.0 );

  fcTmplt->tmpltNorm *= fcTmplt->tmpltNorm;

  fcTmplt->tmpltNorm *= distNorm * distNorm;

  /* pN constants */
  c0  = 3.0/(eta*128.0);
  c10 = 3715.0/756.0 + eta*55.0/9.0;
  c15 = -16*LAL_PI;
  c20 = 15293365.0/508032.0 + eta*(27145.0/504.0 + eta*3085.0/72.0);

  /* x1 */
  x1 = pow( LAL_PI * m * LAL_MTSUN_SI * deltaF, -1.0/3.0 );

  /* frequency cutoffs */
  fHi = 1.0 / (6.0 * sqrt(6.0) * LAL_PI * m * LAL_MTSUN_SI);
  kmin = params->fLow / deltaF > 1 ? params->fLow / deltaF : 1;
  kmax = fHi / deltaF < numPoints/2 ? fHi / deltaF : numPoints/2;

  /* compute psi0: used in range reduction */
  {
    REAL4 x = x1 * xfac[kmin];
    REAL4 psi = c0 * x * ( c20 + x * ( c15 + x * (c10 + x * x ) ) );
    psi0 = -2 * LAL_PI * ( floor ( 0.5 * psi / LAL_PI ) );
  }


  /*
   *
   * calculate the stationary phase chirp
   *
   */


  for ( k = kmin; k < kmax ; ++k )
  {
    REAL4 x = x1 * xfac[k];
    REAL4 psi = c0 * x * ( c20 + x * ( c15 + x * (c10 + x * x ) ) );
    REAL4 psi1 = psi + psi0;
    REAL4 psi2;

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
    if ( psi1 < -LAL_PI/2 )
    {
      psi1 = -LAL_PI - psi1;
      psi2 = psi1 * psi1;
      /* XXX minus sign added because of new sign convention for fft */
      expPsi[k].im = - psi1 * ( 1 + psi2 * ( s2 + psi2 * s4 ) );
      expPsi[k].re = -1 - psi2 * ( c2 + psi2 * c4 );
    }
    else if ( psi1 > LAL_PI/2 )
    {
      psi1 = LAL_PI - psi1;
      psi2 = psi1 * psi1;
      /* XXX minus sign added because of new sign convention for fft */
      expPsi[k].im = - psi1 * ( 1 + psi2 * ( s2 + psi2 * s4 ) );
      expPsi[k].re = -1 - psi2 * ( c2 + psi2 * c4 );
    }
    else
    {
      psi2 = psi1 * psi1;
      /* XXX minus sign added because of new sign convention for fft */
      expPsi[k].im = - psi1 * ( 1 + psi2 * ( s2 + psi2 * s4 ) );
      expPsi[k].re = 1 + psi2 * ( c2 + psi2 * c4 );
    }

  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}



/* <lalVerbatim file="FindChirpBCVTemplateCP"> */
void
LALFindChirpBCVTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpSPTmpltParams     *params
    )
/* </lalVerbatim> */
{
  UINT4        numPoints  = 0;
  REAL4        deltaF     = 0.0;
  REAL4        m          = 0.0;  
  REAL4        chirpMass  = 0.0;
  REAL4        eta        = 0.0; 
  REAL4        mu         = 0.0;  /* now only used in normalisation */
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
  INT4         kmin       = 0;     
  INT4         kmax       = 0;    
  REAL4        distNorm;
  const REAL4  cannonDist = 1.0; /* Mpc */
  
  INITSTATUS( status, "LALFindChirpSPTemplate", FINDCHIRPSPTEMPLATEC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */
  

  /* check that the output structures exist */
  ASSERT( fcTmplt, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( fcTmplt->data, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( fcTmplt->data->data, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  /* check that the parameter structure exists */         
  ASSERT( params, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  /* check that the timestep is positive */
  ASSERT( params->deltaT > 0, status,                        
      FINDCHIRPSPH_EDELT, FINDCHIRPSPH_MSGEDELT );

  /* check that the input exists */
  ASSERT( tmplt, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );


  /*
   *
   * compute the stationary phase template
   *
   */


  /* set up pointers */
  expPsi    = fcTmplt->data->data;
  xfac      = params->xfacVec->data;
  numPoints = fcTmplt->data->length;

  /* zero output */
  memset( expPsi, 0, numPoints * sizeof(COMPLEX8) );

  /* psi coefficients */
  psi00 = tmplt->psi0;        /* BCV only uses psi0, psi15:            */
  psi05 = 0.0; /*tmplt->psi1;*/ /* -> psi1,2,4 don't exist in tmplt      */
  psi10 = 0.0; /*tmplt->psi2;*/ /* -> use if statements to define these? */
  psi15 = tmplt->psi3;        /* & which name convention to use?       */
  psi20 = 0.0; /*tmplt->psi4;*/
/* work needed here... */

  /* parameters */
  deltaF = 1.0 / ( (REAL4) params->deltaT * (REAL4) numPoints ); /* not defined in tmplt */
  m      = - psi15 / ( 16 * LAL_PI * LAL_PI * psi00 );
  eta    = 3 / ( 128 * psi00 * pow( LAL_PI * m, 5.0/3.0 ) );       
  mu     = eta * m;
/* work needed here... check definitions are correct */

  /* defining chirp mass (to the power of 5/3 => rename) */
  chirpMass = pow( 1.0 / LAL_PI, 5.0/3.0) * ( 3 / (128 * psi00));
/* work needed here... required? */

  /* template dependent normalisation */
  distNorm = 2.0 * LAL_MRSUN_SI / (cannonDist * 1.0e6 * LAL_PC_SI);
  distNorm *= params->dynRange;

  fcTmplt->tmpltNorm = sqrt( (5.0*mu) / 96.0 ) *
    pow( m / (LAL_PI*LAL_PI) , 1.0/3.0 ) *
    pow( LAL_MTSUN_SI / (REAL4) params->deltaT, -1.0/6.0 );
  
  fcTmplt->tmpltNorm *= fcTmplt->tmpltNorm;
  
  fcTmplt->tmpltNorm *= distNorm * distNorm;

  
  
  /* x1 */   /* does this explanation suffice? */
  x1 = pow( deltaF, -1.0/3.0 );
/* work needed here ... check x1 */

  /* frequency cutoffs */
  fHi  = tmplt->fendBCV;
  kmin = params->fLow / deltaF > 1 ? params->fLow / deltaF : 1;
  kmax = fHi / deltaF < numPoints/2 ? fHi / deltaF : numPoints/2;
  
  /* compute psi0: used in range reduction */
  {
    REAL4 x    = x1 * xfac[kmin];
    REAL4 psi  = psi20 + x * ( psi15 + x * ( psi10 + x * ( psi05 + x * ( psi00 ))));
          psi0 = -2 * LAL_PI * ( floor ( 0.5 * psi / LAL_PI ) );
  }
/* work needed here... check psi */


  /*
   *
   * calculate the stationary phase chirp
   *
   */


  for ( k = kmin; k < kmax ; ++k )
    {
      REAL4 x    = x1 * xfac[k];
      REAL4 psi  = psi20 + x * ( psi15 + x * ( psi10 + x * ( psi05 + x * ( psi00 ))));
      REAL4 psi1 = psi + psi0;
      REAL4 psi2;  /* defining psi2 every time through the loop necessary? */
/* work needed here... check psi */  

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
      expPsi[k].im = - sin(psi1);
      expPsi[k].re =   cos(psi1);
/* work needed here... expensive computation method */
    }

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


