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


static int UsePrec = 0;

/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */

static REAL8 XLALInspiralSpinFactorizedFlux (REAL8Vector * values,
					     EOBNonQCCoeffs * nqcCoeffs,
					     const REAL8 omega,
					     SpinEOBParams * ak,
					     const REAL8 H,
					     const UINT4 lMax,
					     const UINT4
					     SpinAlignedEOBversion);

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

static REAL8
XLALInspiralSpinFactorizedFlux (REAL8Vector * values,	/**< dynamical variables */
				EOBNonQCCoeffs * nqcCoeffs,
							/**< pre-computed NQC coefficients */
				const REAL8 omega,	/**< orbital frequency */
				SpinEOBParams * ak,	/**< physical parameters */
				const REAL8 H,		/**< real Hamiltonian */
				const UINT4 lMax,	/**< upper limit of the summation over l */
				UNUSED const UINT4 SpinAlignedEOBversion  /**< 1 for SEOBNRv1, 2 for SEOBNRv2, 4 for SEOBNRv4 */
  )
{

  if ( nqcCoeffs==NULL ) {
      XLAL_ERROR_REAL8 (XLAL_EINVAL);
    }
  REAL8 flux = 0.0;
  REAL8 v;
  REAL8 omegaSq;
  COMPLEX16 hLM;
  INT4 l, m;

  //EOBNonQCCoeffs nqcCoeffs;

#ifndef LAL_NDEBUG
  if (!values || !ak)
    {
      XLAL_ERROR_REAL8 (XLAL_EFAULT);
    }
#endif

  if (lMax < 2)
    {
      XLAL_ERROR_REAL8 (XLAL_EINVAL);
    }

  /* Omegs is the derivative of phi */
  omegaSq = omega * omega;

  v = cbrt (omega);

    COMPLEX16 hT= 0.;
//  printf( "v = %.16e\n", v );
  for (l = 2; l <= (INT4) lMax; l++)
    {
      for (m = 1; m <= l; m++)
	{
	  INT4 use_optimized_v2 = 0;
        if ( (ak->seobCoeffs->tidal1->lambda2Tidal != 0. && ak->seobCoeffs->tidal1->omega02Tidal != 0.) || (ak->seobCoeffs->tidal2->lambda2Tidal != 0. && ak->seobCoeffs->tidal2->omega02Tidal != 0.) ) {
            if (XLALSimIMRSpinEOBGetSpinFactorizedWaveform
                (&hLM, values, v, H, l, m, ak, use_optimized_v2
                 ) == XLAL_FAILURE)
            {
                XLAL_ERROR_REAL8 (XLAL_EFUNC);
            }
            if (XLALSimIMRSpinEOBWaveformTidal
                (&hT, values, v, l, m, ak
                 ) == XLAL_FAILURE)
            {
                XLAL_ERROR_REAL8 (XLAL_EFUNC);
            }
        }
        else {
            if (XLALSimIMRSpinEOBFluxGetSpinFactorizedWaveform
                (&hLM, values, v, H, l, m, ak, use_optimized_v2,
                 NULL) == XLAL_FAILURE)
            {
                XLAL_ERROR_REAL8 (XLAL_EFUNC);
            }
        }

	  /* For the 2,2 mode, we apply NQC correction to the flux */
	  if (l == 2 && m == 2)
	    {
	      COMPLEX16 hNQC;
	      XLALSimIMREOBNonQCCorrection (&hNQC, values, omega, nqcCoeffs);
	      /* Eq. 16 */
	      hLM *= hNQC;
	    }
        hLM += hT;

	  //printf( "l = %d, m = %d, mag(hLM) = %.17e, omega = %.16e\n", l, m, sqrt(creal(hLM)*creal(hLM)+cimag(hLM)*cimag(hLM)), omega );
	  /* Eq. 13 */
	  flux +=
	    (REAL8) (m * m) * omegaSq * (creal (hLM) * creal (hLM) +
					 cimag (hLM) * cimag (hLM));
	}
    }
  return flux * LAL_1_PI / 8.0;
}

#endif /* _LALSIMIMRSPINEOBFACTORIZEDFLUX_C */
