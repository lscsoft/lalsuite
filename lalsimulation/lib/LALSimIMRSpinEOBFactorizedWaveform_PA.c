/*
*  Copyright (C) 2010 Craig Robinson, Yi Pan, Prayush Kumar, Andrea Taracchini
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
 * \author Craig Robinson, Yi Pan, Prayush Kumar, Andrea Taracchini
 *
 * \brief Function to compute the factorized waveform as uses in the SEOBNRv1 model.
 * Waveform expressions are given by
 * Taracchini et al. ( PRD 86, 024011 (2012), arXiv 1202.0790 ).
 * All equation numbers in this file refer to equations of this paper,
 * unless otherwise specified.
 * Coefficients of the so-called "deltalm" terms are given by
 * Damour et al. PRD 79, 064004 (2009) [arXiv:0811.2069] and Pan et al. PRD 83, 064003 (2011) [arXiv:1006.0431],
 * henceforth DIN and PBFRT.
 *
 * Functions to compute the factorized waveform for the purpose of computing the
 * flux, i.e. returning only the absolute value of the multipoles. The tail term
 * Tlm is used in its resummed form, given by Eq. (42) of Damour, Nagar and
 * Bernuzzi, PRD 87 084035 (2013) [arxiv:1212.4357], called DNB here.
 *
 */

#ifndef _LALSIMIMRSPINEOBFACTORIZEDWAVEFORMPA_C
#define _LALSIMIMRSPINEOBFACTORIZEDWAVEFORMPA_C

#include <complex.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMREOBNRv2.h"
#include "LALSimIMRSpinEOB.h"

#include "LALSimIMREOBNewtonianMultipole.c"
#include "LALSimIMRSpinEOBAuxFuncs.c"
#include "LALSimIMREOBNQCCorrection.c"
#include "LALSimIMRSpinEOBHamiltonian.c"

/* OPTIMIZED */
#include "LALSimIMRSpinEOBHamiltonianOptimized.c"
/* END OPTIMIZED */

/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */


UNUSED static INT4 XLALSimIMRSpinEOBGetSpinFactorizedWaveform_PA (COMPLEX16 * restrict hlm, REAL8Vector * restrict values, const REAL8 v, const REAL8 Hreal, const INT4 l, const INT4 m, SpinEOBParams * restrict params, REAL8 vPhi
								      /**< Use optimized v2? */
  );




/**
 * This function calculates hlm mode factorized-resummed waveform
 * for given dynamical variables.
 * Eq. 17 and the entire Appendix of the paper.
 */
static INT4
XLALSimIMRSpinEOBGetSpinFactorizedWaveform_PA (COMPLEX16 * restrict hlm,
						      /**< OUTPUT, hlm waveforms */
					    REAL8Vector * restrict values,
						      /**< dyanmical variables */
					    const REAL8 v,
						      /**< velocity */
					    const REAL8 Hreal,
						      /**< real Hamiltonian */
					    const INT4 l,
						      /**< l mode index */
					    const INT4 m,
						      /**< m mode index */
					    SpinEOBParams * restrict params,
						       /**< Spin EOB parameters */
					    REAL8 vPhi
						      /**< NonKepler coeff */
  )
{
  /* Status of function calls */
  INT4 status;
  INT4 i;
  INT4 use_hm = 0;


  REAL8 eta;
  REAL8 r, pp, Omega, v2, vh, vh3, k, hathatk, eulerlogxabs;	//pr
  REAL8 Slm, deltalm, rholm;
  COMPLEX16 auxflm = 0.0;
  COMPLEX16 Tlm, rholmPwrl;
  COMPLEX16 hNewton;
  gsl_sf_result lnr1, arg1, z2;

  /* Non-Keplerian velocity */
  REAL8 vPhi2;

  /* Pre-computed coefficients */
  FacWaveformCoeffs *hCoeffs = params->eobParams->hCoeffs;

  if (abs (m) > l)
    {
      XLAL_ERROR (XLAL_EINVAL);
    }
  if (m == 0)
    {
      XLAL_ERROR (XLAL_EINVAL);
    }
  use_hm = params->use_hm;
  eta = params->eobParams->eta;

  /* Check our eta was sensible */
  if (eta > 0.25 && eta < 0.25 +1e-4) {
      eta = 0.25;
  }
  if (eta > 0.25)
    {
      XLALPrintError
	("XLAL Error - %s: Eta = %.16f seems to be > 0.25 - this isn't allowed!\n",
	 __func__,eta);
      XLAL_ERROR (XLAL_EINVAL);
    }
  /*else if ( eta == 0.25 && m % 2 )
     {
     // If m is odd and dM = 0, hLM will be zero
     memset( hlm, 0, sizeof( COMPLEX16 ) );
     return XLAL_SUCCESS;
     } */

  r = values->data[0];
  //pr    = values->data[2];
  pp = values->data[3];

  v2 = v * v;
  Omega = v2 * v;
  vh3 = Hreal * Omega;
  vh = cbrt (vh3);
  eulerlogxabs = LAL_GAMMA + log (2.0 * (REAL8) m * v);

  /* Calculate the non-Keplerian velocity */
  if (params->alignedSpins)
    {
     
      vPhi = r * cbrt (vPhi);
      vPhi *= Omega;
      vPhi2 = vPhi * vPhi;
    }
  else
    {
      vPhi = v;
      vPhi2 = v2;
    }
  /* Calculate the newtonian multipole, 1st term in Eq. 17, given by Eq. A1 */
  // YP: !!!!! SEOBNRv3devel temporary change !!!!!
  status = XLALSimIMRSpinEOBCalculateNewtonianMultipole (&hNewton, vPhi2, r,
							 values->data[1],
							 (UINT4) l, m,
							 params->eobParams);
#if YPPrecWave
  status = XLALSimIMRSpinEOBCalculateNewtonianMultipole (&hNewton, vPhi2, r,
							 values->data[12] +
							 values->data[13],
							 (UINT4) l, m,
							 params->eobParams);
#endif
  // YP: !!!!! SEOBNRv3devel temporary change !!!!!
  if (status == XLAL_FAILURE)
    {
      XLAL_ERROR (XLAL_EFUNC);
    }
  /* Calculate the source term, 2nd term in Eq. 17, given by Eq. A5 */
  if (((l + m) % 2) == 0)
    {
      Slm = (Hreal * Hreal - 1.) / (2. * eta) + 1.;
    }
  else
    {
      Slm = v * pp;
    }
  //printf( "Hreal = %e, Slm = %e, eta = %e\n", Hreal, Slm, eta );

  /* Calculate the Tail term, 3rd term in Eq. 17, given by Eq. A6 */
  k = m * Omega;
  hathatk = Hreal * k;
  XLAL_CALLGSL (status =
		gsl_sf_lngamma_complex_e (l + 1.0, -2.0 * hathatk, &lnr1,
					  &arg1));
  if (status != GSL_SUCCESS)
    {
      XLALPrintError ("XLAL Error - %s: Error in GSL function: %s\n",
                      __func__, gsl_strerror (status));
      XLAL_ERROR (XLAL_EFUNC);
    }
  XLAL_CALLGSL (status = gsl_sf_fact_e (l, &z2));
  if (status != GSL_SUCCESS)
    {
      XLALPrintError ("XLAL Error - %s: Error in GSL function: %s\n",
                      __func__, gsl_strerror (status));
      XLAL_ERROR (XLAL_EFUNC);
    }
  Tlm =
    cexp ((lnr1.val + LAL_PI * hathatk) +
	  I * (arg1.val + 2.0 * hathatk * log (4.0 * k / sqrt (LAL_E))));
  Tlm /= z2.val;


  /* Calculate the residue phase and amplitude terms */
  /* deltalm is the 4th term in Eq. 17, delta 22 given by Eq. A15, others  */
  /* rholm is the 5th term in Eq. 17, given by Eqs. A8 - A14 */
  /* auxflm is a special part of the 5th term in Eq. 17, given by Eq. A15 */
  /* Actual values of the coefficients are defined in the next function of this file */
  switch (l)
    {
    case 2:
      switch (abs (m))
	{
	case 2:
	  deltalm = vh3 * (hCoeffs->delta22vh3 + vh3 * (hCoeffs->delta22vh6
							+
							vh * vh *
							(hCoeffs->delta22vh9 *
							 vh))) +
	    hCoeffs->delta22v5 * v * v2 * v2 +
	    hCoeffs->delta22v6 * v2 * v2 * v2 +
	    hCoeffs->delta22v8 * v2 * v2 * v2 * v2;
	  rholm =
	    1. + v2 * (hCoeffs->rho22v2 +
		       v * (hCoeffs->rho22v3 +
			    v * (hCoeffs->rho22v4 +
				 v * (hCoeffs->rho22v5 +
				      v * (hCoeffs->rho22v6 +
					   hCoeffs->rho22v6l * eulerlogxabs +
					   v * (hCoeffs->rho22v7 +
						v * (hCoeffs->rho22v8 +
						     hCoeffs->rho22v8l *
						     eulerlogxabs +
						     (hCoeffs->rho22v10 +
						      hCoeffs->rho22v10l *
						      eulerlogxabs) *
						     v2)))))));
	  break;
	case 1:
	  {
	    deltalm = vh3 * (hCoeffs->delta21vh3 + vh3 * (hCoeffs->delta21vh6
							  +
							  vh *
							  (hCoeffs->
							   delta21vh7 +
							   (hCoeffs->
							    delta21vh9) * vh *
							   vh))) +
	      hCoeffs->delta21v5 * v * v2 * v2 +
	      hCoeffs->delta21v7 * v2 * v2 * v2 * v;
	    rholm =
	      1. + v * (hCoeffs->rho21v1 +
			v * (hCoeffs->rho21v2 +
			     v * (hCoeffs->rho21v3 +
				  v * (hCoeffs->rho21v4 +
				       v * (hCoeffs->rho21v5 +
					    v * (hCoeffs->rho21v6 +
						 hCoeffs->rho21v6l *
						 eulerlogxabs +
						 v * (hCoeffs->rho21v7 +
						      hCoeffs->rho21v7l *
						      eulerlogxabs +
						      v * (hCoeffs->rho21v8 +
							   hCoeffs->rho21v8l *
							   eulerlogxabs +
							   (hCoeffs->
							    rho21v10 +
							    hCoeffs->
							    rho21v10l *
							    eulerlogxabs) *
							   v2))))))));
                 if (use_hm){
                   //RC: This terms are in Eq.A11 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
                   auxflm = v * (hCoeffs->f21v1 + v2 * (hCoeffs->f21v3 + v * hCoeffs->f21v4 + v2 * (hCoeffs->f21v5 + v * hCoeffs->f21v6 + v2*hCoeffs->f21v7c)));
                 }
                 else{
                   auxflm = v * hCoeffs->f21v1;
                 }
               }
	  break;
	default:
	  XLAL_ERROR (XLAL_EINVAL);
	  break;
	}
      break;
    case 3:
      switch (m)
	{
	case 3:
    deltalm =
        vh3 * (hCoeffs->delta33vh3 +
             vh3 * (hCoeffs->delta33vh6 + hCoeffs->delta33vh9 * vh3)) +
             hCoeffs->delta33v5 * v * v2 * v2;
    if(use_hm){
      //RC: This terms are in Eq.A6 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
      rholm =
          1. + v2 * (hCoeffs->rho33v2 +
                     v * (hCoeffs->rho33v3 +
                          v * (hCoeffs->rho33v4 +
                               v * (hCoeffs->rho33v5 +
                                    v * (hCoeffs->rho33v6 +
                                         hCoeffs->rho33v6l * eulerlogxabs +
                                         v * (hCoeffs->rho33v7 +
                                              v * (hCoeffs->rho33v8 +
                                               hCoeffs->rho33v8l *
                                               eulerlogxabs +
                                             v2*(hCoeffs->rho33v10 + hCoeffs->rho33v10l*eulerlogxabs))))))));
       //RC: This terms are in Eq.A10 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
        auxflm = v * (v2 * (hCoeffs->f33v3 + v * (hCoeffs->f33v4 + v * (hCoeffs->f33v5  + v * hCoeffs->f33v6)))) + _Complex_I * vh3 * vh3 * hCoeffs->f33vh6;
            }
            else
            {
            rholm =
                1. + v2 * (hCoeffs->rho33v2 +
                           v * (hCoeffs->rho33v3 +
                                v * (hCoeffs->rho33v4 +
                                     v * (hCoeffs->rho33v5 +
                                          v * (hCoeffs->rho33v6 +
                                               hCoeffs->rho33v6l * eulerlogxabs +
                                               v * (hCoeffs->rho33v7 +
                                                    v * (hCoeffs->rho33v8 +
                                                     hCoeffs->rho33v8l *
                                                     eulerlogxabs)))))));
            auxflm = v * v2 * hCoeffs->f33v3;
            }
	  break;
	case 2:
	  deltalm =
	    vh3 * (hCoeffs->delta32vh3 +
		   vh * (hCoeffs->delta32vh4 +
			 vh * vh * (hCoeffs->delta32vh6 +
				    hCoeffs->delta32vh9 * vh3)));
	  rholm =
	    1. + v * (hCoeffs->rho32v +
		      v * (hCoeffs->rho32v2 +
			   v * (hCoeffs->rho32v3 +
				v * (hCoeffs->rho32v4 +
				     v * (hCoeffs->rho32v5 +
					  v * (hCoeffs->rho32v6 +
					       hCoeffs->rho32v6l *
					       eulerlogxabs +
					       (hCoeffs->rho32v8 +
						hCoeffs->rho32v8l *
						eulerlogxabs) * v2))))));
	  break;
	case 1:
	  deltalm = vh3 * (hCoeffs->delta31vh3 + vh3 * (hCoeffs->delta31vh6
							+
							vh *
							(hCoeffs->delta31vh7 +
							 hCoeffs->delta31vh9 *
							 vh * vh))) +
	    hCoeffs->delta31v5 * v * v2 * v2;
	  rholm =
	    1. + v2 * (hCoeffs->rho31v2 +
		       v * (hCoeffs->rho31v3 +
			    v * (hCoeffs->rho31v4 +
				 v * (hCoeffs->rho31v5 +
				      v * (hCoeffs->rho31v6 +
					   hCoeffs->rho31v6l * eulerlogxabs +
					   v * (hCoeffs->rho31v7 +
						(hCoeffs->rho31v8 +
						 hCoeffs->rho31v8l *
						 eulerlogxabs) * v))))));
	  auxflm = v * v2 * hCoeffs->f31v3;
	  break;
	default:
	  XLAL_ERROR (XLAL_EINVAL);
	  break;
	}
      break;
    case 4:
      switch (m)
	{
	case 4:
    if(use_hm){
      //RC: This terms are in Eq.A15 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
      deltalm = vh3 * (hCoeffs->delta44vh3 + vh3 * (hCoeffs->delta44vh6 + vh3 * hCoeffs->delta44vh9))
        + hCoeffs->delta44v5 * v2 * v2 * v;
      //RC: This terms are in Eq.A8 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
      rholm = 1. + v2 * (hCoeffs->rho44v2
               + v * (hCoeffs->rho44v3 + v * (hCoeffs->rho44v4
                      +
                      v *
                      (hCoeffs->
                       rho44v5 +
                       v * (hCoeffs->
                        rho44v6 +
                        hCoeffs->
                        rho44v6l *
                        eulerlogxabs +
                      v2 *( hCoeffs->rho44v8 + hCoeffs->rho44v8l * eulerlogxabs
                        +v2 * (hCoeffs->rho44v10 + hCoeffs->rho44v10l * eulerlogxabs) ) )))));
    }
    else{
	  deltalm = vh3 * (hCoeffs->delta44vh3 + hCoeffs->delta44vh6 * vh3)
	    + hCoeffs->delta44v5 * v2 * v2 * v;

    rholm = 1. + v2 * (hCoeffs->rho44v2
             + v * (hCoeffs->rho44v3 + v * (hCoeffs->rho44v4
                    +
                    v *
                    (hCoeffs->
                     rho44v5 +
                     (hCoeffs->
                      rho44v6 +
                      hCoeffs->
                      rho44v6l *
                      eulerlogxabs) *
                     v))));
    }

	  break;
	case 3:
	  deltalm = vh3 * (hCoeffs->delta43vh3 + vh * (hCoeffs->delta43vh4
						       +
						       hCoeffs->delta43vh6 *
						       vh * vh));
	  rholm =
	    1. + v * (hCoeffs->rho43v +
		      v * (hCoeffs->rho43v2 +
			   v2 * (hCoeffs->rho43v4 +
				 v * (hCoeffs->rho43v5 +
				      (hCoeffs->rho43v6 +
				       hCoeffs->rho43v6l * eulerlogxabs) *
				      v))));
	  auxflm = v * hCoeffs->f43v;
	  break;
	case 2:
	  deltalm = vh3 * (hCoeffs->delta42vh3 + hCoeffs->delta42vh6 * vh3);
	  rholm = 1. + v2 * (hCoeffs->rho42v2
			     + v * (hCoeffs->rho42v3 +
				    v * (hCoeffs->rho42v4 +
					 v * (hCoeffs->rho42v5 +
					      (hCoeffs->rho42v6 +
					       hCoeffs->rho42v6l *
					       eulerlogxabs) * v))));
	  break;
	case 1:
	  deltalm = vh3 * (hCoeffs->delta41vh3 + vh * (hCoeffs->delta41vh4
						       +
						       hCoeffs->delta41vh6 *
						       vh * vh));
	  rholm =
	    1. + v * (hCoeffs->rho41v +
		      v * (hCoeffs->rho41v2 +
			   v2 * (hCoeffs->rho41v4 +
				 v * (hCoeffs->rho41v5 +
				      (hCoeffs->rho41v6 +
				       hCoeffs->rho41v6l * eulerlogxabs) *
				      v))));
	  auxflm = v * hCoeffs->f41v;
	  break;
	default:
	  XLAL_ERROR (XLAL_EINVAL);
	  break;
	}
      break;
    case 5:
      switch (m)
	{
	case 5:
  if(use_hm){
    //RC: This terms are in Eq.A16 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
    deltalm =
       vh3 *(hCoeffs->delta55vh3 +vh3*(hCoeffs->delta55vh6 +vh3 *(hCoeffs->delta55vh9)))
       + hCoeffs->delta55v5 * v2 * v2 * v;
    //RC: This terms are in Eq.A9 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
    rholm =
      1. + v2 * (hCoeffs->rho55v2 +
           v * (hCoeffs->rho55v3 +
          v * (hCoeffs->rho55v4 +
         v * (hCoeffs->rho55v5 +
            v * (hCoeffs->rho55v6 + hCoeffs->rho55v6l * eulerlogxabs +
          v2 * (hCoeffs->rho55v8 + hCoeffs->rho55v8l * eulerlogxabs +
        v2 * (hCoeffs->rho55v10 + hCoeffs->rho55v10l * eulerlogxabs )))))));
    //RC: This terms are in Eq.A12 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028
    auxflm = v2 * v * (hCoeffs->f55v3 + v * (hCoeffs->f55v4 + v * (hCoeffs->f55v5c) ));
  }
  else
  {
	  deltalm =
	    hCoeffs->delta55vh3 * vh3 + hCoeffs->delta55v5 * v2 * v2 * v;
	  rholm =
	    1. + v2 * (hCoeffs->rho55v2 +
		       v * (hCoeffs->rho55v3 +
			    v * (hCoeffs->rho55v4 +
				 v * (hCoeffs->rho55v5 +
				      hCoeffs->rho55v6 * v))));
            }
	  break;
	case 4:
	  deltalm = vh3 * (hCoeffs->delta54vh3 + hCoeffs->delta54vh4 * vh);
	  rholm = 1. + v2 * (hCoeffs->rho54v2 + v * (hCoeffs->rho54v3
						     + hCoeffs->rho54v4 * v));
	  break;
	case 3:
	  deltalm = hCoeffs->delta53vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho53v2
			     + v * (hCoeffs->rho53v3 +
				    v * (hCoeffs->rho53v4 +
					 hCoeffs->rho53v5 * v)));
	  break;
	case 2:
	  deltalm = vh3 * (hCoeffs->delta52vh3 + hCoeffs->delta52vh4 * vh);
	  rholm = 1. + v2 * (hCoeffs->rho52v2 + v * (hCoeffs->rho52v3
						     + hCoeffs->rho52v4 * v));
	  break;
	case 1:
	  deltalm = hCoeffs->delta51vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho51v2
			     + v * (hCoeffs->rho51v3 +
				    v * (hCoeffs->rho51v4 +
					 hCoeffs->rho51v5 * v)));
	  break;
	default:
	  XLAL_ERROR (XLAL_EINVAL);
	  break;
	}
      break;
    case 6:
      switch (m)
	{
	case 6:
	  deltalm = hCoeffs->delta66vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho66v2 + v * (hCoeffs->rho66v3
						     + hCoeffs->rho66v4 * v));
	  break;
	case 5:
	  deltalm = hCoeffs->delta65vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho65v2 + hCoeffs->rho65v3 * v);
	  break;
	case 4:
	  deltalm = hCoeffs->delta64vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho64v2 + v * (hCoeffs->rho64v3
						     + hCoeffs->rho64v4 * v));
	  break;
	case 3:
	  deltalm = hCoeffs->delta63vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho63v2 + hCoeffs->rho63v3 * v);
	  break;
	case 2:
	  deltalm = hCoeffs->delta62vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho62v2 + v * (hCoeffs->rho62v3
						     + hCoeffs->rho62v4 * v));
	  break;
	case 1:
	  deltalm = hCoeffs->delta61vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho61v2 + hCoeffs->rho61v3 * v);
	  break;
	default:
	  XLAL_ERROR (XLAL_EINVAL);
	  break;
	}
      break;
    case 7:
      switch (m)
	{
	case 7:
	  deltalm = hCoeffs->delta77vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
	  break;
	case 6:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho76v2 * v2;
	  break;
	case 5:
	  deltalm = hCoeffs->delta75vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho75v2 + hCoeffs->rho75v3 * v);
	  break;
	case 4:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho74v2 * v2;
	  break;
	case 3:
	  deltalm = hCoeffs->delta73vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
	  break;
	case 2:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho72v2 * v2;
	  break;
	case 1:
	  deltalm = hCoeffs->delta71vh3 * vh3;
	  rholm = 1. + v2 * (hCoeffs->rho71v2 + hCoeffs->rho71v3 * v);
	  break;
	default:
	  XLAL_ERROR (XLAL_EINVAL);
	  break;
	}
      break;
    case 8:
      switch (m)
	{
	case 8:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho88v2 * v2;
	  break;
	case 7:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho87v2 * v2;
	  break;
	case 6:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho86v2 * v2;
	  break;
	case 5:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho85v2 * v2;
	  break;
	case 4:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho84v2 * v2;
	  break;
	case 3:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho83v2 * v2;
	  break;
	case 2:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho82v2 * v2;
	  break;
	case 1:
	  deltalm = 0.0;
	  rholm = 1. + hCoeffs->rho81v2 * v2;
	  break;
	default:
	  XLAL_ERROR (XLAL_EINVAL);
	  break;
	}
      break;
    default:
      XLAL_ERROR (XLAL_EINVAL);
      break;
    }

  /* Raise rholm to the lth power */
  rholmPwrl = 1.0;
  for(i = 0 ; i < l ; i++) rholmPwrl *= rholm;
  /* In the equal-mass odd m case, there is no contribution from nonspin terms,
   * and the only contribution comes from the auxflm term that is proportional to chiA (asymmetric spins).
   * In this case, we must ignore the nonspin terms directly, since the leading term defined by
   * CalculateThisMultipolePrefix in LALSimIMREOBNewtonianMultipole.c is not zero (see comments there).
   */
  if (eta == 0.25 && m % 2)
    {
      rholmPwrl = auxflm;
    }
  else
    {
      rholmPwrl += auxflm;
    }
  /* Put all factors in Eq. 17 together */
  *hlm = Tlm * cexp (I * deltalm) * Slm * rholmPwrl;
  *hlm *= hNewton;
  /*if (r > 8.5)
     {
     printf("YP::FullWave: Reh = %.16e, Imh = %.16e, hAmp = %.16e, hPhi = %.16e\n",creal(*hlm),cimag(*hlm),cabs(*hlm),carg(*hlm));
     } */
  return XLAL_SUCCESS;
}


#endif /* _LALSIMIMRSPINEOBFACTORIZEDWAVEFORMPA */
