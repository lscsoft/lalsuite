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
 
  FILE        *fpTmplt      =  NULL;
 
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

fprintf (stdout, "expPsi    = %e \n", expPsi);
fprintf (stdout, "xfac      = %e \n", xfac);
fprintf (stdout, "numPoints = %d \n", numPoints);


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

fprintf (stdout, "psi00 = %e \n", psi00);
fprintf (stdout, "psi15 = %e \n", psi15);



 /* parameters */
  deltaF = 1.0 / ( (REAL4) params->deltaT * (REAL4) numPoints ); 
  /* XXX not defined in tmplt */
  m      = - psi15 / ( 16 * LAL_PI * LAL_PI * psi00 );       /* total mass */
  eta    = 3 / ( 128 * psi00 * pow( LAL_PI * m, 5.0/3.0 ) ); /* symmetric mass ratio */      
  mu     = eta * m;                                          /* ? */
  /* defns checked against Bank documentation */

fprintf (stdout, "m     = %e \n", m);
fprintf (stdout, "eta   = %e \n", eta);
fprintf (stdout, "mu    = %e \n", mu);


  /* defining chirp mass (to the power of 5/3 => rename) */
  chirpMass = pow( 1.0 / LAL_PI, 5.0/3.0) * ( 3 / (128 * psi00));
  /*  chirpMass = m * pow(eta,3/5); */

  /*
   *
   * template dependent normalisation
   * i need to check this & make necessary changes
   *
   */


  /* template dependent normalisation */
  distNorm = 2.0 * LAL_MRSUN_SI / (cannonDist * 1.0e6 * LAL_PC_SI);
  distNorm *= params->dynRange;

  fcTmplt->tmpltNorm = sqrt( (5.0*mu) / 96.0 ) *
    pow( m / (LAL_PI*LAL_PI) , 1.0/3.0 ) *
    pow( LAL_MTSUN_SI / (REAL4) params->deltaT, -1.0/6.0 );
  
  fcTmplt->tmpltNorm *= fcTmplt->tmpltNorm;
  
  fcTmplt->tmpltNorm *= distNorm * distNorm;
  /* is this the same for BCVSpin? */
  
  
  /* x1 */   /* does this explanation suffice? */
  x1 = pow( deltaF, -1.0/3.0 );
  /* XXX work needed here ... check x1 */

fprintf (stdout, "x1 = %e \n", x1);


  /* frequency cutoffs */
  fHi  = tmplt->fFinal;
  kmin = params->fLow / deltaF > 1 ? params->fLow / deltaF : 1;
  kmax = fHi / deltaF < numPoints/2 ? fHi / deltaF : numPoints/2;
 
fprintf (stdout, "fHi     = %e \n", fHi);
fprintf (stdout, "fLow    = %e \n", params->fLow);
fprintf (stdout, "deltaF  = %e \n", deltaF);

fprintf (stdout, "kmin    = %d \n", kmin);
fprintf (stdout, "kmax    = %d \n", kmax);

 
  /* compute psi0: used in range reduction */
  {
    REAL4 x    = x1 * xfac[kmin];
    REAL4 psi  = 
      psi20 + (x * x) * ( psi15 + x * ( psi10 + x * ( psi05 + x * ( psi00 ))));
    psi0 = -2 * LAL_PI * ( floor ( 0.5 * psi / LAL_PI ) );

fprintf (stdout, "xfac[kmin] = %e \n", xfac[kmin]);

fprintf (stdout, "x = %e \n", x);

fprintf (stdout, "psi = %e \n", psi);
fprintf (stdout, "psi0 = %e \n", psi0);
 
 }
  /* XXX work needed here... check psi */



  /*
   *
   * calculate the stationary phase chirp
   *
   */

fprintf (stdout, "just before loop in template code \n");

fpTmplt  =  fopen("tmplt.dat","w");

  for ( k = kmin; k < kmax ; ++k )
    {
      REAL4 x    = x1 * xfac[k];
      REAL4 psi  = 
        psi20 + (x * x) * ( psi15 + x * ( psi10 + x * ( psi05 + x * ( psi00 ))));
      REAL4 psi1 = psi + psi0;
      REAL4 psi2;  /* defining psi2 every time through the loop necessary? */
                   /* where is psi2 used? */

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

	fprintf (fpTmplt, "%d\t%e\t%e\n",k, expPsi[k].im, expPsi[k].re);

      /* XXX work needed here... expensive computation method */

fprintf (stdout, "expPsi[k].im = %e \n", expPsi[k].im);
fprintf (stdout, "expPsi[k].re = %e \n", expPsi[k].re);



    }

fclose (fpTmplt);
  /*code*/

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

