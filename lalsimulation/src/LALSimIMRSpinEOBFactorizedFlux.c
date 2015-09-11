/*
*  Copyright (C) 2010 Craig Robinson, Yi Pan, Prayush Kumar (minor changes)
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


/**
 * \author Craig Robinson, Yi Pan
 *
 * \brief Function to compute the factorized flux as uses in the SEOBNRv1
 * model. Flux function given in
 * Taracchini et al. ( PRD 86, 024011 (2012), arXiv 1202.0790 ).
 * All equation numbers in this file refer to equations of this paper,
 * unless otherwise specified.
 */

#ifndef _LALSIMIMRSPINEOBFACTORIZEDFLUX_C
#define _LALSIMIMRSPINEOBFACTORIZEDFLUX_C

#include <complex.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMREOBNRv2.h"
#include "LALSimIMRSpinEOB.h"

#include "LALSimIMRSpinEOBAuxFuncs.c"
#include "LALSimIMREOBNQCCorrection.c"
#include "LALSimIMRSpinEOBFactorizedWaveform.c"


int UsePrec = 0;

/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */

static REAL8 XLALInspiralSpinFactorizedFlux(
                      REAL8Vector           *values,
                      EOBNonQCCoeffs        *nqcCoeffs,
                      const REAL8           omega,
                      SpinEOBParams         *ak,
                      const REAL8            H,
                      const INT4             lMax,
                      const UINT4            SpinAlignedEOBversion
                     );

/*------------------------------------------------------------------------------------------
 *
 *          Defintions of functions.
 *
 *------------------------------------------------------------------------------------------
 */

/**
 * This function calculates the spin factorized-resummed GW energy flux
 * for given dynamical variables.
 */

static REAL8 XLALInspiralSpinFactorizedFlux(
                      REAL8Vector           *values,    /**< dynamical variables */
                      EOBNonQCCoeffs        *nqcCoeffs, /**< pre-computed NQC coefficients */
                      const REAL8           omega,      /**< orbital frequency */
                      SpinEOBParams         *ak,        /**< physical parameters */
                      const REAL8            H,         /**< real Hamiltonian */
                      const INT4             lMax,      /**< upper limit of the summation over l */
                      const UINT4            SpinAlignedEOBversion /**< 1 for SEOBNRv1, 2 for SEOBNRv2 */
                     )

{

  REAL8 flux = 0.0;
  REAL8 v;
  REAL8 omegaSq;
  COMPLEX16 hLM;
  INT4 l, m;

  //EOBNonQCCoeffs nqcCoeffs;

#ifndef LAL_NDEBUG
  if ( !values || !ak )
  {
    XLAL_ERROR_REAL8( XLAL_EFAULT );
  }
#endif

  if ( lMax < 2 )
  {
    XLAL_ERROR_REAL8( XLAL_EINVAL );
  }

  /* Omegs is the derivative of phi */
  omegaSq = omega*omega;

  v = cbrt( omega );

  /* Update the factorized multipole coefficients, w.r.t. new spins */
  if ( UsePrec )
  {
	/* Assume that initial conditions are available at this point, to
	 * compute the chiS and chiA parameters.
     * Calculate the values of chiS and chiA, as given in Eq.16 of
	 * Precessing EOB paper. Assuming \vec{L} to be pointing in the
     * direction of \vec{r}\times\vec{p} */
    /* TODO: Check the mass scaling of spins */
	REAL8 rcrossp[3], rcrosspMag, s1dotL, s2dotL;
	REAL8 chiS, chiA, tplspin;

	rcrossp[0] = values->data[1]*values->data[5] - values->data[2]*values->data[4];
	rcrossp[1] = values->data[2]*values->data[3] - values->data[0]*values->data[5];
	rcrossp[2] = values->data[0]*values->data[4] - values->data[1]*values->data[3];
	rcrosspMag = sqrt(rcrossp[0]*rcrossp[0] + rcrossp[1]*rcrossp[1] +
        rcrossp[2]*rcrossp[2]);

	rcrossp[0] /= rcrosspMag;
	rcrossp[1] /= rcrosspMag;
	rcrossp[2] /= rcrosspMag;

	s1dotL = values->data[6]*rcrossp[0] + values->data[7]*rcrossp[1]
			+ values->data[8]*rcrossp[2];
	s2dotL = values->data[9]*rcrossp[0] + values->data[10]*rcrossp[1]
			+ values->data[11]*rcrossp[2];

	chiS = 0.5 * (s1dotL + s2dotL);
	chiA = 0.5 * (s1dotL - s2dotL);

	/* Compute the test-particle limit spin of the deformed-Kerr background */
	/* TODO: Check this is actually the way it works in LAL */
	switch ( SpinAlignedEOBversion )
	{
     case 1:
       tplspin = 0.0;
       break;
     case 2:
       tplspin = (1.-2.* ak->eobParams->eta) * chiS + (ak->eobParams->m1
			- ak->eobParams->m2)/(ak->eobParams->m1 + ak->eobParams->m2) * chiA;
       break;
     default:
       XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
       XLAL_ERROR( XLAL_EINVAL );
       break;
	}

	/* ************************************************* */
	/* Re-Populate the Waveform structures               */
	/* ************************************************* */

	/* Re-compute the spinning coefficients for hLM */
	if ( XLALSimIMREOBCalcSpinFacWaveformCoefficients( ak->eobParams->hCoeffs,
			ak->eobParams->m1, ak->eobParams->m2, ak->eobParams->eta,
			tplspin, chiS, chiA, SpinAlignedEOBversion ) == XLAL_FAILURE )
	{
		XLALDestroyREAL8Vector( values );
		XLAL_ERROR( XLAL_EFUNC );
	}
  }

//  printf( "v = %.16e\n", v );
  for ( l = 2; l <= lMax; l++ )
  {
    for ( m = 1; m <= l; m++ )
    {

      if ( XLALSimIMRSpinEOBFluxGetSpinFactorizedWaveform( &hLM, values, v, H,
            l, m, ak ) == XLAL_FAILURE )
      {
        XLAL_ERROR_REAL8( XLAL_EFUNC );
      }
      /* For the 2,2 mode, we apply NQC correction to the flux */
      if ( l == 2 && m == 2 )
      {
        COMPLEX16 hNQC;
        /*switch ( SpinAlignedEOBversion )
        {
          case 1:
            XLALSimIMRGetEOBCalibratedSpinNQC( &nqcCoeffs, l, m, ak->eobParams->eta, ak->a );
            break;
          case 2:
            // XLALSimIMRGetEOBCalibratedSpinNQCv2( &nqcCoeffs, l, m, ak->eobParams->eta, ak->a );
            XLALSimIMRGetEOBCalibratedSpinNQC3D( &nqcCoeffs, l, m, ak->eobParams->eta, ak->a, (ak->chi1 - ak->chi2)/2. );
            break;
          default:
            XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
            XLAL_ERROR( XLAL_EINVAL );
            break;
        }*/
        XLALSimIMREOBNonQCCorrection( &hNQC, values, omega, nqcCoeffs );
        /* Eq. 16 */
        hLM *= hNQC;
      }
      //printf( "l = %d, m = %d, mag(hLM) = %.17e, omega = %.16e\n", l, m, sqrt(creal(hLM)*creal(hLM)+cimag(hLM)*cimag(hLM)), omega );
      /* Eq. 13 */
      flux += (REAL8)(m * m) * omegaSq * ( creal(hLM)*creal(hLM) + cimag(hLM)*cimag(hLM) );
    }
  }
  return flux * LAL_1_PI / 8.0;
}

#endif /* _LALSIMIMRSPINEOBFACTORIZEDFLUX_C */
