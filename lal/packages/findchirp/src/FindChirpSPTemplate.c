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

#include <lal/LALStdlib.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>


NRCSID (FINDCHIRPSPTEMPLATEC, "$Id$");

void
LALFindChirpSPTemplateInit (
    LALStatus                  *status,
    FindChirpSPTmpltParams    **output,
    FindChirpInitParams        *params
    )
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



void
LALFindChirpSPTemplateFinalize (
    LALStatus                  *status,
    FindChirpSPTmpltParams    **output
    )
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



void
LALFindChirpSPTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpSPTmpltParams     *params
    )
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
   * caclulate the stationary phase chirp
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
      expPsi[k].im = psi1 * ( 1 + psi2 * ( s2 + psi2 * s4 ) );
      expPsi[k].re = -1 - psi2 * ( c2 + psi2 * c4 );
    }
    else if ( psi1 > LAL_PI/2 )
    {
      psi1 = LAL_PI - psi1;
      psi2 = psi1 * psi1;
      expPsi[k].im = psi1 * ( 1 + psi2 * ( s2 + psi2 * s4 ) );
      expPsi[k].re = -1 - psi2 * ( c2 + psi2 * c4 );
    }
    else
    {
      psi2 = psi1 * psi1;
      expPsi[k].im = psi1 * ( 1 + psi2 * ( s2 + psi2 * s4 ) );
      expPsi[k].re = 1 + psi2 * ( c2 + psi2 * c4 );
    }

  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}




