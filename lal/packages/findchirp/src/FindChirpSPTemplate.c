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
\label{ss:FindChirpSPTemplate.c}

Provides functions to create stationary phase inspiral templates in a
form that can be used by the \texttt{FindChirpFilter()} function.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpSPTemplateCP}
\idx{LALFindChirpSPTemplate()}

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

\vfill{\footnotesize\input{FindChirpSPTemplateCV}}
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
LALFindChirpSPTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpTmpltParams       *params
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
  INT4          k          = 0;
  INT4          kmin       = 0;
  INT4          kmax       = 0;

  REAL4         distNorm;
  const REAL4   cannonDist = 1.0; /* Mpc */

  /* pn constants */
  REAL4 c0, c10, c15, c20; 

  /* variables used to compute chirp time */
  REAL4 c0T, c2T, c3T, c4T, xT, x2T, x3T, x4T, x8T;

  /* chebychev coefficents for expansion of sin and cos */
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
  ASSERT( fcTmplt, status, 
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( fcTmplt->data, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( fcTmplt->data->data, status, 
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  /* check that the parameter structure exists */
  ASSERT( params, status, 
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( params->xfacVec, status, 
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( params->xfacVec->data, status, 
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  /* check that the timestep is positive */
  ASSERT( params->deltaT > 0, status, 
      FINDCHIRPSPH_EDELT, FINDCHIRPSPH_MSGEDELT );

  /* check that the input exists */
  ASSERT( tmplt, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  /* check that the parameter structure is set */
  /* to the correct waveform approximant       */
  if ( params->approximant != TaylorF2 )
    ABORT( status, FINDCHIRPSPH_EMAPX, FINDCHIRPSPH_MSGEMAPX );


  /*
   *
   * compute the stationary phase template
   *
   */


  /* set up pointers */
  expPsi = fcTmplt->data->data;
  xfac = params->xfacVec->data;
  numPoints = fcTmplt->data->length;

  /* store the waveform approximant */
  fcTmplt->tmplt.approximant = params->approximant;

  /* zero output */
  memset( expPsi, 0, numPoints * sizeof(COMPLEX8) );

  /* parameters */
  deltaF = 1.0 / ( (REAL4) params->deltaT * (REAL4) numPoints );
  m      = (REAL4) tmplt->totalMass;
  eta    = (REAL4) tmplt->eta;
  mu     = (REAL4) tmplt->mu;

  if ( m <= 0 || eta <= 0 || mu <= 0 )
  {
    ABORT( status, FINDCHIRPH_EMASS, FINDCHIRPH_MSGEMASS );
  }

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
  tmplt->fCutoff = 1.0 / (6.0 * sqrt(6.0) * LAL_PI * m * LAL_MTSUN_SI);
  kmin = params->fLow / deltaF > 1 ? params->fLow / deltaF : 1;
  kmax = tmplt->fCutoff / deltaF < numPoints/2 ? 
    tmplt->fCutoff / deltaF : numPoints/2;

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


  /*
   *
   * compute the length of the stationary phase chirp
   *
   */


  c0T = 5 * m * LAL_MTSUN_SI / (256 * eta);
  c2T = 743.0/252.0 + eta * 11.0/3.0;
  c3T = -32 * LAL_PI/3;
  c4T = 3058673.0/508032.0 + eta * (5429.0/504.0 + eta * 617.0/72.0);
  xT  = pow( LAL_PI * m * LAL_MTSUN_SI * params->fLow, 1.0/3.0);
  x2T = xT * xT;
  x3T = xT * x2T;
  x4T = x2T * x2T;
  x8T = x4T * x4T;
  tmplt->tC = c0T * (1 + c2T * x2T + c3T * x3T + c4T * x4T) / x8T;

  /* copy the template parameters to the finchirp template structure */
  memcpy( &(fcTmplt->tmplt), tmplt, sizeof(InspiralTemplate) );

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
