/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpSBCVTemplate.c
 *
 * Author: Brown D. A., Spinning BCV modifications by Jones, G
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
#include <lal/FindChirpSP.h>

NRCSID (FINDCHIRPSBCVTEMPLATEC, "$Id$");

/*documenation later*/
void
LALFindChirpBCVSpinTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpSPTmpltParams     *params
		)
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
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( fcTmplt->data, status, 
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( fcTmplt->data->data, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  /* check that the parameter structure exists */         
  ASSERT( params, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  /* check that the parameter structure is set */
  /* to the correct waveform approximant       */
  ASSERT( params->approximant == BCVSpin, status, 
      FINDCHIRPSPH_EMAPX, FINDCHIRPSPH_MSGEMAPX );

  /* check that the timestep is positive */
  ASSERT( params->deltaT > 0, status,                        
      FINDCHIRPSPH_EDELT, FINDCHIRPSPH_MSGEDELT );

  /* check that the input exists */
  ASSERT( tmplt, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

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


  /*code*/

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

