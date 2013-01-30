/*
*  Copyright (C) 2007 Darren Woods, Drew Keppel, Duncan Brown, Gareth Jones,
*             Jolien Creighton, Patrick Brady, Thomas Cokelaer, Evan Ochsner
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/*-----------------------------------------------------------------------
 *
 * File Name: FindChirpSPTemplate.c
 *
 * Author: Brown D. A.
 *
 *-----------------------------------------------------------------------
 */

/**

\author Brown, D. A.
\file
\ingroup FindChirpSP_h

\brief Provides functions to create stationary phase inspiral templates in a
form that can be used by the <tt>FindChirpFilter()</tt> function.

\heading{Prototypes}

The function <tt>LALFindChirpSPTemplate()</tt> creates the stationary phase
template as described by the algorithm below.

\heading{Algorithm}

Blah.

\heading{Uses}
\code
LALCalloc()
LALFree()
LALCreateVector()
LALDestroyVector()
\endcode

\heading{Notes}

*/

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>

double
XLALFindChirpChirpTime (double m1,
			double m2,
			double fLower,
			int order)
{

    /* variables used to compute chirp time */
    double c0T, c2T, c3T, c4T, c5T, c6T, c6LogT, c7T;
    double xT, x2T, x3T, x4T, x5T, x6T, x7T, x8T;
    double m = m1 + m2;
    double eta = m1 * m2 / m / m;

    c0T = c2T = c3T = c4T = c5T = c6T = c6LogT = c7T = 0.;

    /* Switch on PN order, set the chirp time coeffs for that order */
    switch (order)
    {
    case 8:
    case 7:
      c7T = LAL_PI * (14809.0 * eta * eta / 378.0 - 75703.0 * eta / 756.0 - 15419335.0 / 127008.0);
    case 6:
      c6T = LAL_GAMMA * 6848.0 / 105.0 - 10052469856691.0 / 23471078400.0 + LAL_PI * LAL_PI * 128.0 / 3.0 + eta * (3147553127.0 / 3048192.0 - LAL_PI * LAL_PI * 451.0 / 12.0) - eta * eta * 15211.0 / 1728.0 + eta * eta * eta * 25565.0 / 1296.0 + log (4.0) * 6848.0 / 105.0;
      c6LogT = 6848.0 / 105.0;
    case 5:
      c5T = 13.0 * LAL_PI * eta / 3.0 - 7729.0 * LAL_PI / 252.0;
    case 4:
      c4T = 3058673.0 / 508032.0 + eta * (5429.0 / 504.0 + eta * 617.0 / 72.0);
      c3T = -32.0 * LAL_PI / 5.0;
      c2T = 743.0 / 252.0 + eta * 11.0 / 3.0;
      c0T = 5.0 * m * LAL_MTSUN_SI / (256.0 * eta);
      break;
    default:
      fprintf (stderr, "ERROR!!!\n");
      break;
    }

    /* This is the PN parameter v evaluated at the lower freq. cutoff */
    xT = pow (LAL_PI * m * LAL_MTSUN_SI * fLower, 1.0 / 3.0);
    x2T = xT * xT;
    x3T = xT * x2T;
    x4T = x2T * x2T;
    x5T = x2T * x3T;
    x6T = x3T * x3T;
    x7T = x3T * x4T;
    x8T = x4T * x4T;

    /* Computes the chirp time as tC = t(v_low)    */
    /* tC = t(v_low) - t(v_upper) would be more    */
    /* correct, but the difference is negligble.   */

    /* This formula works for any PN order, because */
    /* higher order coeffs will be set to zero.     */

    return c0T * (1 + c2T * x2T + c3T * x3T + c4T * x4T + c5T * x5T + (c6T + c6LogT * log (xT)) * x6T + c7T * x7T) / x8T;
}




void
LALFindChirpSPTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpTmpltParams       *params
    )

{
  UINT4         numPoints  = 0;
  REAL4         deltaF     = 0.0;
  REAL4         m          = 0.0;
  REAL4         eta        = 0.0;
  REAL4         mu         = 0.0;
  REAL4         S1z        = 0.0;
  REAL4         S2z        = 0.0;
  REAL4         mass_delta = 0.0;
  REAL4         chis       = 0.0;
  REAL4         chia       = 0.0;
  REAL4         chi1       = 0.0;
  REAL4         chi2       = 0.0;
  REAL4         qm_def1    = 0.0;
  REAL4         qm_def2    = 0.0;
  REAL4         pn_beta    = 0.0;
  REAL4         pn_sigma   = 0.0;
  REAL4         pn_gamma   = 0.0;
  COMPLEX8     *expPsi     = NULL;
  REAL4        *xfac       = NULL;
  REAL4         x1         = 0.0;
  REAL4         psi        = 0.0;
  REAL4         psi0       = 0.0;
  INT4          k          = 0;
  INT4          f          = 0;
  INT4          kmin       = 0;
  INT4          kmax       = 0;
  REAL4         fLow       = -1;
  CHAR          infomsg[512];

  REAL4         distNorm;
  const REAL4   cannonDist = 1.0; /* Mpc */

  /* pn constants */
  REAL4 c0, c10, c15, c20, c25, c25Log, c30, c30Log, c35, c40P;
  REAL4 x;

  /* chebychev coefficents for expansion of sin and cos */
  const REAL4 s2 = -0.16605;
  const REAL4 s4 =  0.00761;
  const REAL4 c2 = -0.49670;
  const REAL4 c4 =  0.03705;


  INITSTATUS(status);
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
  if ( params->approximant != FindChirpSP )
  {
    ABORT( status, FINDCHIRPSPH_EMAPX, FINDCHIRPSPH_MSGEMAPX );
  }
  LALInfo( status, "Generating template using FindChirpSP" );


  /*
   *
   * compute the stationary phase template
   *
   */


  /* set up pointers */
  expPsi = fcTmplt->data->data;
  xfac = params->xfacVec->data;
  numPoints = 2 * (fcTmplt->data->length - 1);

  /* set the waveform approximant */
  tmplt->approximant = params->approximant;

  /* set the pN order of the template */
  tmplt->order = params->order;

  /* zero output */
  memset( expPsi, 0, fcTmplt->data->length * sizeof(COMPLEX8) );

  /* parameters */
  deltaF = 1.0 / ( (REAL4) params->deltaT * (REAL4) numPoints );
  m      = (REAL4) tmplt->totalMass;
  eta    = (REAL4) tmplt->eta;
  mu     = (REAL4) tmplt->mu;
  S1z    = tmplt->spin1[2];
  S2z    = tmplt->spin2[2];
  mass_delta = (tmplt->mass1 - tmplt->mass2) / (m);
  chis   = 0.5 * (tmplt->spin1[2] + tmplt->spin2[2]);
  chia   = 0.5 * (tmplt->spin1[2] - tmplt->spin2[2]);
  chi1 = tmplt->mass1 / m;
  chi2 = tmplt->mass2 / m;
  qm_def1 = 1; /* The QM deformability parameters */
  qm_def2 = 1; /* This is 1 for black holes and larger for neutron stars */

  /* Eq. (6.23) in arXiv:0810.5336 */
  pn_beta = (113./12.- 19./3. * eta) * chis + 113./12. * mass_delta * chia;
  
  /* See Eq. (6.24) in arXiv:0810.5336 */
  /* 9b,c,d in arXiv:astro-ph/0504538 */
  pn_sigma = eta * (721./48. *S1z*S2z-247./48.*S1z*S2z);
  pn_sigma += (720*qm_def1 - 1)/96.0 * (chi1*chi1*S1z*S1z);
  pn_sigma += (720*qm_def2 - 1)/96.0 * (chi2*chi2*S2z*S2z);
  pn_sigma -= (240*qm_def1 - 7)/96.0 * (chi1*chi1*S1z*S1z);
  pn_sigma -= (240*qm_def2 - 7)/96.0 * (chi2*chi2*S2z*S2z);

  /* See Eq. (6.25) in arXiv:0810.5336 */
  pn_gamma = (732985./2268. - 24260./81. * eta - 340./9. * eta * eta ) * chis;
  pn_gamma += (732985./2268. +140./9.0 * eta) * chia * mass_delta;

  fprintf(stderr,"%e %e %e \n",pn_beta,pn_sigma,pn_gamma);

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

  if ( lalDebugLevel & LALINFO )
  {
    snprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg),
        "tmpltNorm = %e\n", fcTmplt->tmpltNorm );
    LALInfo( status, infomsg );
  }

  /* Initialize all PN phase coeffs to zero. */
  c0 = c10 = c15 = c20 = c25 = c25Log = 0.;
  c30 = c30Log = c35 = c40P = 0.;

  /* Switch on PN order, set the appropriate phase coeffs for that order */
  switch( params->order )
  {
    case LAL_PNORDER_PSEUDO_FOUR:
      c40P = 3923.0;
    case LAL_PNORDER_THREE_POINT_FIVE:
      c35 = LAL_PI*(77096675.0/254016.0 + eta*378515.0/1512.0
            - eta*eta*74045.0/756.0);
    case LAL_PNORDER_THREE:
      c30 = 11583231236531.0/4694215680.0 - LAL_GAMMA*6848.0/21.0
            - LAL_PI*LAL_PI*640.0/3.0 + eta*(LAL_PI*LAL_PI*2255.0/12.0
            - 15737765635.0/3048192.0) + eta*eta*76055.0/1728.0
            - eta*eta*eta*127825.0/1296.0 - 6848.0*log(4.0)/21.0;
      c30Log = -6848.0/21.0;
    case LAL_PNORDER_TWO_POINT_FIVE:
      c25 = LAL_PI*38645.0/756.0 - LAL_PI*eta*65.0/9.0 - pn_gamma;
      c25Log = 3*c25;
    case LAL_PNORDER_TWO:
      c20 = 15293365.0/508032.0 + eta*(27145.0/504.0 + eta*3085.0/72.0);
      c20 -= 10. * pn_sigma;
      c15 = -16*LAL_PI + 4.*pn_beta;
      c10 = 3715.0/756.0 + eta*55.0/9.0;
      c0  = 3.0/(eta*128.0);
      break;
    default:
      ABORT( status, FINDCHIRPSPH_EORDR, FINDCHIRPSPH_MSGEORDR );
      break;
  }

  /* x1 */
  x1 = pow( LAL_PI * m * LAL_MTSUN_SI * deltaF, -1.0/3.0 );

  /* frequency cutoffs */
  if (params->dynamicTmpltFlow)
  {
    /* Dynamic lower cutoff
     * Work out longest length for template
     * Keep a few extra sample points for safety */
    REAL4 currTime;
    REAL4 maxT = ((REAL4) params->deltaT * (REAL4) (numPoints-12))/4.;
    maxT -= 0.5 * (REAL4) params->invSpecTrunc * (REAL4) params->deltaT;
    fLow = -1;
    for (f=1; f < 100; f++)
    {
      currTime = XLALFindChirpChirpTime( tmplt->mass1,
                                        tmplt->mass2,
                                        (double) f,
                                        params->order);
      if (currTime < maxT)
      {
        fLow = (REAL4) f;
        break;
      }
    }   
    /* If nothing passed then fail */
    if ( fLow < 0)
    {
      ABORT( status, FINDCHIRPH_EFLOX, FINDCHIRPH_MSGEFLOX );
    }
  }
  else
  {
    fLow = params->fLow;
  }

  kmin = fLow / deltaF > 1 ? fLow / deltaF : 1;
  kmax = tmplt->fFinal / deltaF < numPoints/2 ?
    tmplt->fFinal / deltaF : numPoints/2;

  /* compute psi0: used in range reduction */

  /* This formula works for any PN order, because */
  /* higher order coeffs will be set to zero.     */

    x = x1 * xfac[kmin];
    psi = c0 * ( x * ( c20 + x * ( c15 + x * (c10 + x * x ) ) )
                + c25 - c25Log * log(x) + (1.0/x)
                * ( c30 - c30Log * log(x) + (1.0/x) * ( c35 - (1.0/x)
                * c40P * log(x) ) ) );
    psi0 = -2 * LAL_PI * ( floor ( 0.5 * psi / LAL_PI ) );


  /*
   *
   * calculate the stationary phase chirp
   *
   */

  /* This formula works for any PN order, because */
  /* higher order coeffs will be set to zero.     */

    for ( k = kmin; k < kmax ; ++k )
    {
      REAL4 x_0 = x1 * xfac[k];
      REAL4 psi_0 = c0 * ( x_0 * ( c20 + x_0 * ( c15 + x_0 * (c10 + x_0 * x_0 ) ) )
                  + c25 - c25Log * log(x_0) + (1.0/x_0) * ( c30 - c30Log * log(x_0)
                  + (1.0/x_0) * ( c35 - (1.0/x_0) * c40P * log(x_0) ) ) );
      REAL4 psi1 = psi_0 + psi0;
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

      /* if reverse chirp bank option selected, switch sign of imag. part */
      if ( params->reverseChirpBank )
      {
        expPsi[k].im = - expPsi[k].im;
      }

    }


  /*
   *
   * compute the length of the stationary phase chirp
   *
   */

    tmplt->tC = XLALFindChirpChirpTime( tmplt->mass1,
					tmplt->mass2,
					fLow,
					params->order);

  /* copy the template parameters to the findchirp template structure */
  memcpy( &(fcTmplt->tmplt), tmplt, sizeof(InspiralTemplate) );

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
