/*
*  Copyright (C) 2007 Duncan Brown, Gareth Jones
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
 * File Name: FindChirpSBCVTemplate.c
 *
 * Author: Brown D. A., Spinning BCV-Modifications: Jones, G
 *
 *-----------------------------------------------------------------------
 */

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpBCVSpin.h>

/**

\author Brown, D. A., Spinning BCV-Modifications: Jones, G
\file
\ingroup FindChirpBCVSpin_h

\brief Provides functions to create spinning BCV detection templates in a form that
can be used by the <tt>FindChirpBCVSpinFilter()</tt> function.

\heading{Prototypes}

The function <tt>LALFindChirpBCVSpinTemplate()</tt> creates the
spinning BCV template as described by the algorithm below.

\heading{Algorithm}

This code calculates a number of quantities required by the
<tt>LALFindChirpBCVSpinFilterSegment()</tt> function.
To improve efficiency we calculate every template dependent
quantity that does not require detector data
in this function.
Since the template bank does not provide values for \f$f_{final}\f$ we calculate
it here. For every combination of \f$\psi_0\f$ and \f$\psi_3\f$ values provided by
the template bank we find the frequency at the last stable orbit using
\f$f_{final} = (-16 \pi \psi_0)/(\psi_3 r_{LSO}^{3/2})\f$
where \f$r_{LSO} = 6M_{\odot}\f$ is the separation of the binaries components.
We then calculate the complex phase of the template
\f$\psi_{NM}(f) = \psi_{initial} + f^{-5/3}(\psi_0 + f \psi_3)\f$
between \f$f_{low}\f$ and \f$f_{final}\f$.
Next we calculate 5 moments which are required to construct the template:
\f{eqnarray}{
I &=& 4 \int_{0}^{\infty} f^{-7/3} \frac{df} {S_{n} (f) } \nonumber\\
J &=& 4 \int_{0}^{\infty} f^{-7/3} \cos ( \beta f^{-2/3})
\frac{df} {S_{n} (f) } \nonumber\\
K &=& 4 \int_{0}^{\infty} f^{-7/3} \sin ( \beta f^{-2/3})
\frac{df} {S_{n} (f) } \nonumber\\
L &=& 2 \int_{0}^{\infty} f^{-7/3} \sin (2\beta f^{-2/3})
\frac{df} {S_{n} (f) } \nonumber\\
M &=& 2 \int_{0}^{\infty} f^{-7/3} \cos (2\beta f^{-2/3})
\frac{df} {S_{n} (f) }
\f}
In practise we integrate between our lowest non-zero frequency sample
point <tt>k=1</tt> (a division-by-zero error would
occur at \f$0 Hz\f$) and the Nyquist frequency <tt>k=numPoints/2</tt>. From
these moments we then find the
orthonormalised amplitude vectors:
\f{eqnarray}{
\mathcal{\widehat{A}}_1(f)   =  \frac { f^{-7/6} }  { I^{1/2}}
\nonumber
\f}

\f{eqnarray}{
\mathcal{\widehat{A}}_2(f)   =  \frac{ f^{-7/6} \bigg
[ \cos(\beta f^{-2/3}) - \frac{J}{I} \bigg ] I^{1/2} }
{ \bigg[ IM + \frac{I^{2}}{2}
- J^{2} \bigg] ^{1/2} }
\nonumber
\f}

\f{eqnarray}{
\mathcal{\widehat{A}}_3(f) & = & \frac{ f^{-7/6} \bigg [
\sin(\beta f^{-2/3})
- \frac{K}{I}
- \frac{IL - JK}{IM + \frac{I^{2}}{2} - J^{2}}
\big[\cos(\beta f^{-2/3}) -\frac{J}{I} \big ]
\bigg ] I^{1/2} }
{ \bigg [
\frac{I^{2}}{2} - IM - K^{2} - \frac{ (IL - JK)^{2}}
{IM + \frac{I^{2}}{2} - J^{2} }
\bigg ] ^{1/2} }
\f}
where \f$\beta\f$ is provided by the template bank code and the
\f$f^{-7/6}\f$ and \f$f^{-2/3}\f$ vectors were calculated previously
in <tt>LALFindChirpDataInit()</tt>. To avoid division-by-zero
errors we explicitly set
\f$\mathcal{\widehat{A}}_2(f) = \mathcal{\widehat{A}}_3(f) = 0\f$
when \f$\beta = 0\f$.


\heading{Uses}
\code
LALCalloc()
LALFree()
LALCreateVector()
LALDestroyVector()
\endcode

\heading{Notes}



*/


void
LALFindChirpBCVSpinTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpTmpltParams       *params,
    FindChirpDataParams        *fcDataParams
    )

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
  REAL4                 rLSOto3by2       = 0.0;

  INITSTATUS(status);
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

  /*
   *
   * compute the BCVSpin template
   *
   */

  /* set up pointers */
  expPsi      = fcTmplt->data->data;
  numPoints   = 2 * (fcTmplt->data->length - 1);
  xfac        = params->xfacVec->data;
  xfac        = params->xfacVec->data;

  /* Reading in data needed to calculate moments */
  wtilde      = fcDataParams->wtildeVec->data;
  ampBCVSpin1 = fcDataParams->ampVecBCVSpin1->data;
  ampBCVSpin2 = fcDataParams->ampVecBCVSpin2->data;

  /* store the waveform approximant */
  tmplt->approximant = BCVSpin;

  /* zero output */
  memset( expPsi, 0, fcTmplt->data->length * sizeof(COMPLEX8) );

  /* psi coefficients */
  psi00 = tmplt->psi0;            /* BCV only uses psi0, psi15:            */
  psi05 = 0.0; /*tmplt->psi1;*/   /* -> psi1,2,4 don't exist in tmplt      */
  psi10 = 0.0; /*tmplt->psi2;*/   /* -> use if statements to define these? */
  psi15 = tmplt->psi3;            /* & which name convention to use?       */
  psi20 = 0.0; /*tmplt->psi4;*/

  /* parameters */
  deltaT        = params->deltaT;
  deltaTto2by3  = pow(deltaT, Twoby3);
  deltaF        = 1.0 / ( (REAL4) params->deltaT * (REAL4) numPoints );
  x1            = pow( deltaF, -1.0/3.0 );
  fFinal        = tmplt->fFinal;
  kmin = params->fLow / deltaF > 1 ? params->fLow / deltaF : 1;
  kmax = fFinal / deltaF < numPoints/2 ? fFinal / deltaF : numPoints/2;
  beta = tmplt->beta;

  /* Preliminary BCVSpin bank does not populate fFinal */
  /* will estimate fFinal form psi0, psi3 as quick fix */

  if (fFinal == 0.0)
  {
        rLSOto3by2 = 14.69693846; /* 6 to 3by2) */

	fFinal = (-psi00 * 16 * LAL_PI) / (psi15 * rLSOto3by2);

        tmplt->fFinal = fFinal;
  }



  /* since we have redefined fFinal we must redefine kmax */

  kmax = fFinal / deltaF < numPoints/2 ? fFinal / deltaF : numPoints/2;

  /* compute psi0: used in range reduction */
  {
    REAL4 x    = x1 * xfac[kmin];
    REAL4 psi  =
    psi20 + (x * x) * ( psi15 + x * ( psi10 + x * ( psi05 + x * ( psi00 ))));
    psi0 = -2 * LAL_PI * ( floor ( 0.5 * psi / LAL_PI ) );
  }

  /*
   *
   * calculate the stationary phase chirp
   *
   */

    for ( k = kmin; k < kmax ; ++k )
    {
      REAL4 x    = x1 * xfac[k];
      REAL4 psi  =
      psi20 + (x * x) * ( psi15 + x * ( psi10 + x * ( psi05 + x * ( psi00 ))));
      REAL4 psi1 = psi + psi0;

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
      expPsi[k].imagf_FIXME =   -sin(psi1);
      expPsi[k].realf_FIXME =   cos(psi1);

    }

  /*
   *
   * Calculating the amplitude vectors A1Vec, A2Vec, A3Vec
   *
   */

  for ( k = kmin; k < kmax; ++k )
  {
          I += ampBCVSpin1[k] * ampBCVSpin1[k] * crealf(wtilde[k]) ;
          J += ampBCVSpin1[k] * ampBCVSpin1[k] * crealf(wtilde[k]) *
                  cos(beta * ampBCVSpin2[k] * deltaTto2by3);
          K += ampBCVSpin1[k] * ampBCVSpin1[k] * crealf(wtilde[k]) *
                  sin(beta * ampBCVSpin2[k] * deltaTto2by3);
          L += ampBCVSpin1[k] * ampBCVSpin1[k] * crealf(wtilde[k]) *
                  sin(2 * beta * ampBCVSpin2[k] * deltaTto2by3);
          M += ampBCVSpin1[k] * ampBCVSpin1[k] * crealf(wtilde[k]) *
                  cos(2 * beta * ampBCVSpin2[k] * deltaTto2by3);
  }

  /* Taking multiplucation outside loop lessens cost */

  I *= 4*deltaF;
  J *= 4*deltaF;
  K *= 4*deltaF;
  L *= 2*deltaF;
  M *= 2*deltaF;

  /* To find absolute values of these moments multiply by (deltaT)^(7/3)  */

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

  A1Vec = fcTmplt->A1BCVSpin->data;
  A2Vec = fcTmplt->A2BCVSpin->data;
  A3Vec = fcTmplt->A3BCVSpin->data;

  memset( A1Vec, 0, ((numPoints/2)+1) * sizeof(REAL8) );
  memset( A2Vec, 0, ((numPoints/2)+1) * sizeof(REAL8) );
  memset( A3Vec, 0, ((numPoints/2)+1) * sizeof(REAL8) );

  A1Vec[0] = 0;
  A2Vec[0] = 0;
  A3Vec[0] = 0;

  if (beta == 0.0)
  {
        for ( k = kmin; k < kmax; ++k )
	{
		A1Vec[k] = ampBCVSpin1[k] / rootI;
		A2Vec[k] = 0.0;
		A3Vec[k] = 0.0;
	 }
  }
  else
    {
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
  	 }
  }

  /* copy the template parameters to the finchirp template structure */
  memcpy( &(fcTmplt->tmplt), tmplt, sizeof(InspiralTemplate) );


  DETATCHSTATUSPTR( status );
  RETURN( status );
}
