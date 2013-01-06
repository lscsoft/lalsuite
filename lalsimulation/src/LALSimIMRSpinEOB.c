/*
*  Copyright (C) 2011 Craig Robinson, Enrico Barausse
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
 * \author Craig Robinson
 *
 * \ file 
 *
 * \brief Functions for producing EOB waveforms for 
 * spinning binaries, as described in Barausse and Buonanno ( arXiv 0912.3517 ).
 */

#include <complex.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/LALAdaptiveRungeKutta4.h>
#include <gsl/gsl_sf_gamma.h>

#include <math.h>

#include "LALSimIMREOBNRv2.h"
#include "LALSimIMRSpinEOB.h"

/* Include static functions */
#include "LALSimIMREOBFactorizedWaveform.c" 
#include "LALSimIMREOBNewtonianMultipole.c"
#include "LALSimIMRSpinEOBInitialConditions.c"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

static int
XLALEOBSpinStopCondition(double UNUSED t,
                           const double values[],
                           double dvalues[],
                           void *funcParams
                          )
{

  SpinEOBParams *params = (SpinEOBParams *)funcParams;
  double omega_x, omega_y, omega_z, omega;
  double r2;

  omega_x = values[1]*dvalues[2] - values[2]*dvalues[1];
  omega_y = values[2]*dvalues[0] - values[0]*dvalues[2];
  omega_z = values[0]*dvalues[1] - values[1]*dvalues[0];

  r2 = values[0]*values[0] + values[1]*values[1] + values[2]*values[2];
  omega = sqrt( omega_x*omega_x + omega_y*omega_y + omega_z*omega_z )/r2;

  //if ( omega < params->eobParams->omega )
  if ( r2 < 36. && omega < params->eobParams->omega )
  {
    return 1;
  }

  params->eobParams->omega = omega;
  return GSL_SUCCESS;
}


int XLALSimIMRSpinEOBCalculateSigmaKerr( REAL8Vector *sigmaKerr,
                                   REAL8        mass1,
                                   REAL8        mass2,
                                   REAL8Vector *s1,
                                   REAL8Vector *s2 )
{

  UINT4 i;
  REAL8 totalMass = mass1 + mass2;

  for ( i = 0; i < 3; i++ )
  {
    sigmaKerr->data[i] = (s1->data[i] + s2->data[i])/(totalMass*totalMass);
  }

  return XLAL_SUCCESS;
}


int XLALSimIMRSpinEOBCalculateSigmaStar( REAL8Vector *sigmaStar,
                                   REAL8        mass1,
                                   REAL8        mass2,
                                   REAL8Vector *s1,
                                   REAL8Vector *s2 )
{

  UINT4 i;
  REAL8 totalMass = mass1 + mass2;

  for ( i = 0; i < 3; i++ )
  {
    sigmaStar->data[i] = (mass2/mass1 * s1->data[i] + mass1/mass2 * s2->data[i])/(totalMass*totalMass);
  }

  return XLAL_SUCCESS;
}


INT4 XLALSimIMRSpinEOBGetSpinFactorizedWaveform( COMPLEX16         * restrict hlm,
				REAL8Vector           * restrict values,
                                const REAL8           v,
                                const REAL8           Hreal,
                                const INT4            l,
                                const INT4            m,
                                SpinEOBParams        * restrict params
                                )
{
        /* Status of function calls */
        INT4 status;
        INT4 i;

        REAL8 eta;	
	REAL8 r, pr, pp, Omega, v2, vh, vh3, k, hathatk, eulerlogxabs;
	REAL8 Slm, deltalm, rholm, rholmPwrl;
	COMPLEX16 Tlm;
        COMPLEX16 hNewton;
	gsl_sf_result lnr1, arg1, z2;

        /* Non-Keplerian velocity */
        REAL8 vPhi, vPhi2;

        /* Pre-computed coefficients */
        FacWaveformCoeffs *hCoeffs = params->eobParams->hCoeffs;

	if ( abs(m) > (INT4) l )
	{
	  XLAL_ERROR( XLAL_EINVAL );
	}
	

        eta = params->eobParams->eta;

        /* Check our eta was sensible */
        if ( eta > 0.25 )
        {
          XLALPrintError("XLAL Error - %s: Eta seems to be > 0.25 - this isn't allowed!\n", __func__ );
          XLAL_ERROR( XLAL_EINVAL );
        }
        else if ( eta == 0.25 && m % 2 )
        {
          /* If m is odd and dM = 0, hLM will be zero */
          memset( hlm, 0, sizeof( COMPLEX16 ) );
          return XLAL_SUCCESS;
        }
        
	r	= values->data[0];
        pr	= values->data[2];
	pp	= values->data[3];

	v2	= v * v;
        Omega   = v2 * v;
        vh3     = Hreal * Omega;
	vh	= cbrt(vh3);
	eulerlogxabs = LAL_GAMMA + log( 2.0 * (REAL8)m * v );

        /* Calculate the non-Keplerian velocity */
        if ( params->alignedSpins )
        {
          vPhi = XLALSimIMRSpinAlignedEOBNonKeplerCoeff( values->data, params );

          if ( XLAL_IS_REAL8_FAIL_NAN( vPhi ) )
          {
            XLAL_ERROR( XLAL_EFUNC );
          }

          vPhi  = r * cbrt(vPhi);
          vPhi *= Omega;
          vPhi2 = vPhi*vPhi;
        }
        else
        {
          vPhi2 = v2;
        }

        /* Calculate the newtonian multipole */
        status = XLALSimIMREOBCalculateNewtonianMultipole( &hNewton, vPhi2, r,
                         values->data[1], (UINT4)l, m, params->eobParams );
        if ( status == XLAL_FAILURE )
        {
          XLAL_ERROR( XLAL_EFUNC );
        }

        /* Calculate the source term */
	if ( ( (l+m)%2 ) == 0)
	{ 
	  Slm = (Hreal*Hreal - 1.)/(2.*eta) + 1.;
	}
	else
	{
	  Slm = v * pp;
	}
        //printf( "Hreal = %e, Slm = %e, eta = %e\n", Hreal, Slm, eta );

        /* Calculate the Tail term */	
	k	= m * Omega;
	hathatk = Hreal * k;
	XLAL_CALLGSL( status = gsl_sf_lngamma_complex_e( l+1.0, -2.0*hathatk, &lnr1, &arg1 ) );
	if (status != GSL_SUCCESS)
	{
	  XLALPrintError("XLAL Error - %s: Error in GSL function\n", __func__ );
	  XLAL_ERROR( XLAL_EFUNC );
	}
	XLAL_CALLGSL( status = gsl_sf_fact_e( l, &z2 ) );
	if ( status != GSL_SUCCESS)
	{
	  XLALPrintError("XLAL Error - %s: Error in GSL function\n", __func__ );
	  XLAL_ERROR( XLAL_EFUNC );
	}
	Tlm = XLALCOMPLEX16Exp( XLALCOMPLEX16Rect( lnr1.val + LAL_PI * hathatk, 
				arg1.val + 2.0 * hathatk * log(4.0*k/sqrt(LAL_E)) ) );
	Tlm = XLALCOMPLEX16DivReal( Tlm, z2.val );


        /* Calculate the residue phase and amplitude terms */
	switch( l )
	{
	  case 2:
	    switch( abs(m) )
	    {
	      case 2:
	        deltalm = vh3*(hCoeffs->delta22vh3 + vh3*(hCoeffs->delta22vh6 
			+ vh*vh*(hCoeffs->delta22vh8 + hCoeffs->delta22vh9*vh))) 
			+ hCoeffs->delta22v5 *v*v2*v2;
		rholm	= 1. + v2*(hCoeffs->rho22v2 + v*(hCoeffs->rho22v3
			+ v*(hCoeffs->rho22v4
			+ v*(hCoeffs->rho22v5 + v*(hCoeffs->rho22v6 
			+ hCoeffs->rho22v6l*eulerlogxabs + v*(hCoeffs->rho22v7 
			+ v*(hCoeffs->rho22v8 + hCoeffs->rho22v8l*eulerlogxabs 
			+ (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs)*v2)))))));
	        break;
	      case 1:
                {
	        deltalm = vh3*(hCoeffs->delta21vh3 + vh3*(hCoeffs->delta21vh6
			+ vh*(hCoeffs->delta21vh7 + (hCoeffs->delta21vh9)*vh*vh))) 
			+ hCoeffs->delta21v5*v*v2*v2 + hCoeffs->delta21v7*v2*v2*v2*v;
		rholm	= 1. + v*(hCoeffs->rho21v1
			+ v*( hCoeffs->rho21v2 + v*(hCoeffs->rho21v3 + v*(hCoeffs->rho21v4 
			+ v*(hCoeffs->rho21v5 + v*(hCoeffs->rho21v6 + hCoeffs->rho21v6l*eulerlogxabs 
			+ v*(hCoeffs->rho21v7 + hCoeffs->rho21v7l * eulerlogxabs 
			+ v*(hCoeffs->rho21v8 + hCoeffs->rho21v8l * eulerlogxabs 
			+ (hCoeffs->rho21v10 + hCoeffs->rho21v10l * eulerlogxabs)*v2))))))));
                }
	        break;
	      default:
                XLAL_ERROR( XLAL_EINVAL );
                break;
	    }
	    break;
	  case 3:
	    switch (m)
	    {
	      case 3:
	        deltalm = vh3*(hCoeffs->delta33vh3 + vh3*(hCoeffs->delta33vh6 + hCoeffs->delta33vh9*vh3)) 
                        + hCoeffs->delta33v5*v*v2*v2 + hCoeffs->delta33v7*v2*v2*v2*v;
		rholm	= 1. + v2*(hCoeffs->rho33v2 + v*(hCoeffs->rho33v3 + v*(hCoeffs->rho33v4 
			+ v*(hCoeffs->rho33v5 + v*(hCoeffs->rho33v6 + hCoeffs->rho33v6l*eulerlogxabs
			+ v*(hCoeffs->rho33v7 + (hCoeffs->rho33v8 + hCoeffs->rho33v8l*eulerlogxabs)*v))))));
	        break;
	      case 2:
		deltalm = vh3*(hCoeffs->delta32vh3 + vh*(hCoeffs->delta32vh4 + vh*vh*(hCoeffs->delta32vh6
			+ hCoeffs->delta32vh9*vh3)));
		rholm	= 1. + v*(hCoeffs->rho32v 
			+ v*(hCoeffs->rho32v2 + v*(hCoeffs->rho32v3 + v*(hCoeffs->rho32v4 + v*(hCoeffs->rho32v5
			+ v*(hCoeffs->rho32v6 + hCoeffs->rho32v6l*eulerlogxabs
			+ (hCoeffs->rho32v8 + hCoeffs->rho32v8l*eulerlogxabs)*v2))))));
		break;
	      case 1:
		deltalm = vh3*(hCoeffs->delta31vh3 + vh3*(hCoeffs->delta31vh6
			+ vh*(hCoeffs->delta31vh7 + hCoeffs->delta31vh9*vh*vh))) 
			+ hCoeffs->delta31v5*v*v2*v2;
		rholm	= 1. + v2*(hCoeffs->rho31v2 + v*(hCoeffs->rho31v3 + v*(hCoeffs->rho31v4 
			+ v*(hCoeffs->rho31v5 + v*(hCoeffs->rho31v6 + hCoeffs->rho31v6l*eulerlogxabs 
			+ v*(hCoeffs->rho31v7 + (hCoeffs->rho31v8 + hCoeffs->rho31v8l*eulerlogxabs)*v))))));
		break;
              default:
                XLAL_ERROR( XLAL_EINVAL );
                break;
	    }
	    break;
	  case 4:
	    switch (m)
	    {
	      case 4:
                
	        deltalm = vh3*(hCoeffs->delta44vh3 + hCoeffs->delta44vh6 *vh3)
                          + hCoeffs->delta44v5*v2*v2*v;
		rholm	= 1. + v2*(hCoeffs->rho44v2
			+ v*( hCoeffs->rho44v3 + v*(hCoeffs->rho44v4
			+ v*(hCoeffs->rho44v5 + (hCoeffs->rho44v6
			+ hCoeffs->rho44v6l*eulerlogxabs)*v))));
	        break;
	      case 3:
	        deltalm = vh3*(hCoeffs->delta43vh3 + vh*(hCoeffs->delta43vh4 
			+ hCoeffs->delta43vh6*vh*vh));
		rholm	= 1. + v*(hCoeffs->rho43v
			+ v*(hCoeffs->rho43v2
			+ v2*(hCoeffs->rho43v4 + v*(hCoeffs->rho43v5
			+ (hCoeffs->rho43v6 + hCoeffs->rho43v6l*eulerlogxabs)*v))));
	        break;
	      case 2:
		deltalm = vh3*(hCoeffs->delta42vh3 + hCoeffs->delta42vh6*vh3);
		rholm	= 1. + v2*(hCoeffs->rho42v2
			+ v*(hCoeffs->rho42v3 + v*(hCoeffs->rho42v4 + v*(hCoeffs->rho42v5
			+ (hCoeffs->rho42v6 + hCoeffs->rho42v6l*eulerlogxabs)*v))));
		break;
	      case 1:
		deltalm = vh3*(hCoeffs->delta41vh3 + vh*(hCoeffs->delta41vh4
			+ hCoeffs->delta41vh6*vh*vh));
		rholm	= 1. + v*(hCoeffs->rho41v 
			+ v*(hCoeffs->rho41v2
			+ v2*(hCoeffs->rho41v4 + v*(hCoeffs->rho41v5 
			+ (hCoeffs->rho41v6 +  hCoeffs->rho41v6l*eulerlogxabs)*v))));
		break;
	      default:
                XLAL_ERROR( XLAL_EINVAL );
                break;
	    }
	    break;
	  case 5:
	    switch (m)
	    {
	      case 5:
	        deltalm = hCoeffs->delta55vh3*vh3 + hCoeffs->delta55v5*v2*v2*v;
		rholm	= 1. + v2*( hCoeffs->rho55v2 
			+ v*(hCoeffs->rho55v3 + v*(hCoeffs->rho55v4 
                        + v*(hCoeffs->rho55v5 + hCoeffs->rho55v6*v))));
	        break;
	      case 4:
		deltalm = vh3*(hCoeffs->delta54vh3 + hCoeffs->delta54vh4*vh);
		rholm	= 1. + v2*(hCoeffs->rho54v2 + v*(hCoeffs->rho54v3
			+ hCoeffs->rho54v4*v));
		break;
	      case 3:
	        deltalm = hCoeffs->delta53vh3 * vh3;
		rholm	= 1. + v2*(hCoeffs->rho53v2 
			+ v*(hCoeffs->rho53v3 + v*(hCoeffs->rho53v4 + hCoeffs->rho53v5*v)));
	        break;
	      case 2:
		deltalm = vh3*(hCoeffs->delta52vh3 + hCoeffs->delta52vh4*vh);
		rholm	= 1. + v2*(hCoeffs->rho52v2 + v*(hCoeffs->rho52v3
			+ hCoeffs->rho52v4*v));
		break;
	      case 1:
		deltalm = hCoeffs->delta51vh3 * vh3;
		rholm	= 1. + v2*(hCoeffs->rho51v2 
			+ v*(hCoeffs->rho51v3 + v*(hCoeffs->rho51v4 + hCoeffs->rho51v5*v)));
		break;
	      default:
                XLAL_ERROR( XLAL_EINVAL );
                break;
	    }
	    break;
	  case 6:
	    switch (m)
	    {
	      case 6:
	        deltalm = hCoeffs->delta66vh3*vh3;
		rholm	= 1. + v2*(hCoeffs->rho66v2 + v*(hCoeffs->rho66v3
			+ hCoeffs->rho66v4*v));
	        break;
	      case 5:
		deltalm = hCoeffs->delta65vh3*vh3;
		rholm	= 1. + v2*(hCoeffs->rho65v2 + hCoeffs->rho65v3*v);
		break;
	      case 4:
		deltalm = hCoeffs->delta64vh3 * vh3;
		rholm	= 1. + v2*(hCoeffs->rho64v2 + v*(hCoeffs->rho64v3
			+ hCoeffs->rho64v4*v));
		break;
	      case 3:
	        deltalm = hCoeffs->delta63vh3 * vh3;
		rholm	= 1. + v2*(hCoeffs->rho63v2 + hCoeffs->rho63v3*v);
	        break;
	      case 2:
		deltalm = hCoeffs->delta62vh3 * vh3;
		rholm	= 1. + v2*(hCoeffs->rho62v2 + v*(hCoeffs->rho62v3
			+ hCoeffs->rho62v4 * v));
		break;
	      case 1:
		deltalm = hCoeffs->delta61vh3 * vh3;
		rholm	= 1. + v2*(hCoeffs->rho61v2 + hCoeffs->rho61v3*v);
		break;
	      default:
                XLAL_ERROR( XLAL_EINVAL );
                break;
	    }
	    break;
	  case 7:
	    switch (m)
	    {
	      case 7:
	        deltalm = hCoeffs->delta77vh3 * vh3;
		rholm	= 1. + v2*(hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
	        break;
	      case 6:
	        deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho76v2 * v2;
	        break;
	      case 5:
		deltalm = hCoeffs->delta75vh3 * vh3;
		rholm	= 1. + v2*(hCoeffs->rho75v2 + hCoeffs->rho75v3*v);
		break;
	      case 4:
		deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho74v2 * v2;
		break;
	      case 3:
	        deltalm = hCoeffs->delta73vh3 *vh3;
		rholm	= 1. + v2*(hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
	        break;
	      case 2:
		deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho72v2 * v2;
		break;
	      case 1:
		deltalm = hCoeffs->delta71vh3 * vh3;
		rholm	= 1. + v2*(hCoeffs->rho71v2 +hCoeffs->rho71v3 * v);
		break;
	      default:
                XLAL_ERROR( XLAL_EINVAL );
                break;
	    }
	    break;
	  case 8:
	    switch (m)
	    {
	      case 8:
	        deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho88v2 * v2;
	        break;
	      case 7:
		deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho87v2 * v2;
		break;
	      case 6:
	        deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho86v2 * v2;
	        break;
	      case 5:
		deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho85v2 * v2;
		break;
	      case 4:
		deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho84v2 * v2;
		break;
	      case 3:
	        deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho83v2 * v2;
	        break;
	      case 2:
		deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho82v2 * v2;
		break;
	      case 1:
		deltalm = 0.0;
		rholm	= 1. + hCoeffs->rho81v2 * v2;
		break;
	      default:
                XLAL_ERROR( XLAL_EINVAL );
                break;
	    }
	    break;
	  default:
            XLAL_ERROR( XLAL_EINVAL );
            break; 
	}

        /* Raise rholm to the lth power */
        rholmPwrl = 1.0;
        i = l;
        while ( i-- )
        {
          rholmPwrl *= rholm;
        }

        //printf( "rholm^l = %.16e, Tlm = %.16e + i %.16e, \nSlm = %.16e, hNewton = %.16e + i %.16e, delta = %.16e\n", rholmPwrl, Tlm.re, Tlm.im, Slm, hNewton.re, hNewton.im, deltalm );

	*hlm = XLALCOMPLEX16MulReal( XLALCOMPLEX16Mul( Tlm, XLALCOMPLEX16Polar( 1.0, deltalm) ), 
				     Slm*rholmPwrl );
        *hlm = XLALCOMPLEX16Mul( *hlm, hNewton );

	return XLAL_SUCCESS;
}


int XLALSimIMRSpinEOBWaveform(
        REAL8TimeSeries **hplus,
        REAL8TimeSeries **hcross,
        LIGOTimeGPS     *tc,
        const REAL8     UNUSED phiC,
        const REAL8     deltaT,
        const REAL8     m1,
        const REAL8     m2,
        const REAL8     fMin,
        const REAL8     r,
        const REAL8     inc,
        const REAL8     spin1[],
        const REAL8     spin2[]
     )
{

  int i;

  REAL8Vector *values = NULL;

  /* Parameters of the system */
  REAL8 mTotal, eta, mTScaled;
  REAL8 amp0, amp;

  /* Variables for the integrator */
  ark4GSLIntegrator       *integrator = NULL;
  REAL8Array              *dynamics   = NULL;
  INT4                    retLen;
  /*REAL8  UNUSED           tMax;*/

  /* Accuracies of adaptive Runge-Kutta integrator */
  const REAL8 EPS_ABS = 1.0e-9;
  const REAL8 EPS_REL = 1.0e-8;

  /* Spins not scaled by the mass */
  REAL8 mSpin1[3], mSpin2[3];

  /* Parameter structures containing important parameters for the model */
  SpinEOBParams           seobParams;
  SpinEOBHCoeffs          seobCoeffs;
  EOBParams               eobParams;
  FacWaveformCoeffs       hCoeffs;
  NewtonMultipolePrefixes prefixes;

  /* Initialize parameters */
  mTotal = m1 + m2;
  mTScaled = mTotal * LAL_MTSUN_SI;
  eta    = m1 * m2 / (mTotal*mTotal);

  amp0 = 4. * mTotal * LAL_MRSUN_SI * eta / r;

  /* TODO: Insert potentially necessary checks on the arguments */

  /* Allocate the values vector to contain the ICs */
  /* For this model, it contains 12 dynamical variables: */
  /* values[0-2]  - x (Cartesian separation vector) */
  /* values[3-5]  - p (Cartesian momentum) */
  /* values[6-8]  - spin of body 1 */
  /* values[9-11] - spin of body 2 */
  if ( !(values = XLALCreateREAL8Vector( 14 )) )
  {
    XLAL_ERROR(  XLAL_ENOMEM );
  }
  memset( values->data, 0, values->length * sizeof( REAL8 ));

  /* Set up structures and calculate necessary PN parameters */
  /* Due to precession, these need to get calculated in every step */
  /* TODO: Only calculate non-spinning parts once */
  memset( &seobParams, 0, sizeof(seobParams) );
  memset( &seobCoeffs, 0, sizeof(seobCoeffs) );
  memset( &eobParams, 0, sizeof(eobParams) );
  memset( &hCoeffs, 0, sizeof( hCoeffs ) );
  memset( &prefixes, 0, sizeof( prefixes ) );

  seobParams.seobCoeffs   = &seobCoeffs;
  seobParams.eobParams    = &eobParams;
  eobParams.hCoeffs       = &hCoeffs;
  eobParams.prefixes      = &prefixes;

  eobParams.m1  = m1;
  eobParams.m2  = m2;
  eobParams.eta = eta;

  memcpy( mSpin1, spin1, sizeof( mSpin1 ) );
  memcpy( mSpin2, spin2, sizeof( mSpin2 ) );

  for ( i = 0; i < 3; i++ )
  {
    mSpin1[i] *= m1*m1;
    mSpin2[i] *= m2*m2;
  }

  /* Populate the initial structures */
  if ( XLALSimIMREOBCalcFacWaveformCoefficients( &hCoeffs, eta) == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  if ( XLALSimIMREOBComputeNewtonMultipolePrefixes( &prefixes, eobParams.m1, eobParams.m2 )
         == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* TODO: Set the initial conditions */

  if ( XLALSimIMRSpinEOBInitialConditions( values, m1, m2, fMin, inc, mSpin1, mSpin2, &seobParams ) == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }
  //exit(0);
#if 0
  values->data[0] = 0.;
  values->data[1] = 12.845228155660482;
  values->data[2] = -4.553894189373296;
  values->data[3] = -0.006165987975074341 /*/ eta*/;
  values->data[4] = 0.10049046440176972 /*/ eta*/;
  values->data[5] = 0.28877341851636174 /*/ eta*/;
  values->data[6] = 4.*m1*m1 * 0.1125;
  values->data[7] = 4.*m1*m1 * -0.09742785792574934;
  values->data[8] = 4.*m1*m1 * -0.16875;
  values->data[9] = 4.*m2*m2 *-0.1060660171779821;
  values->data[10] =4.*m2*m2 * 6.938893903907228e-18;
  values->data[11] =4.*m2*m2 * -0.10606601717798211;
#endif
#if 0
  values->data[0] = 0.;
  values->data[1] = 0.;
  values->data[2] = 50.;
  values->data[3] = 0.;
  values->data[4] = -0.1453043197514489;
  values->data[5] = -6.644349740254723e-06;
  values->data[6] = 0.5 * m1*m1;
  values->data[7] = 0.;
  values->data[8] = 0.;
  values->data[9] = 0.5 * m2*m2;
  values->data[10] = 0.;
  values->data[11] = 0.;
#endif

  /* Initialize the GSL integrator */
  if (!(integrator = XLALAdaptiveRungeKutta4Init(14, XLALSpinHcapNumericalDerivative, XLALEOBSpinStopCondition, EPS_ABS, EPS_REL)))
  {
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  integrator->stopontestonly = 1;

  retLen = XLALAdaptiveRungeKutta4( integrator, &seobParams, values->data, 0., 20./mTScaled, deltaT/mTScaled, &dynamics );
  if ( retLen == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  printf("To be the man, you've got to beat the man! Woooooooo!!!!\n" );

  REAL8 *posVecx = dynamics->data+retLen;
  REAL8 *posVecy = dynamics->data+2*retLen;
  REAL8 *posVecz = dynamics->data+3*retLen;
  REAL8 *momVecx = dynamics->data+4*retLen;
  REAL8 *momVecy = dynamics->data+5*retLen;
  REAL8 *momVecz = dynamics->data+6*retLen;
  REAL8 *s1Vecx = dynamics->data+7*retLen;
  REAL8 *s1Vecy = dynamics->data+8*retLen;
  REAL8 *s1Vecz = dynamics->data+9*retLen;
  REAL8 *s2Vecx = dynamics->data+10*retLen;
  REAL8 *s2Vecy = dynamics->data+11*retLen;
  REAL8 *s2Vecz = dynamics->data+12*retLen;
  REAL8 *vphi   = dynamics->data+13*retLen;

  FILE *out = fopen( "seobDynamics.dat", "w" );

  for ( i = 0; i < retLen; i++ )
  {
    fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", i*deltaT/mTScaled, posVecx[i], posVecy[i], posVecz[i], momVecx[i], momVecy[i], momVecz[i],
              s1Vecx[i]/(4.*m1*m1), s1Vecy[i]/(4.*m1*m1), s1Vecz[i]/(4.*m1*m1), s2Vecx[i]/(4.*m2*m2), s2Vecy[i]/(4.*m2*m2), s2Vecz[i]/(4.*m2*m2) );
  }
  fclose( out );

  /* We can now calculate the waveform */
  REAL8 vX, vY, vZ, rCrossV_x, rCrossV_y, rCrossV_z, omega, vOmega;
  REAL8 magPosVec, LNhx, LNhy, LNhz, magL, alpha;

  REAL8TimeSeries *hPlusTS  = XLALCreateREAL8TimeSeries( "H_PLUS", tc, 0.0, deltaT, &lalStrainUnit, retLen );
  REAL8TimeSeries *hCrossTS = XLALCreateREAL8TimeSeries( "H_CROSS", tc, 0.0, deltaT, &lalStrainUnit, retLen );

  for ( i = 0; i < retLen; i++ )
  {
    for ( unsigned int j = 0; j < values->length; j++ )
    {
      values->data[j] = dynamics->data[(j+1)*retLen + i];
    }

    vX = XLALSpinHcapNumDerivWRTParam( 3, values->data, &seobParams );
    vY = XLALSpinHcapNumDerivWRTParam( 4, values->data, &seobParams );
    vZ = XLALSpinHcapNumDerivWRTParam( 5, values->data, &seobParams );

    rCrossV_x = posVecy[i] * vZ - posVecz[i] * vY;
    rCrossV_y = posVecz[i] * vX - posVecx[i] * vZ;
    rCrossV_z = posVecx[i] * vY - posVecy[i] * vX;

    magPosVec = sqrt(posVecx[i]*posVecx[i] + posVecy[i]*posVecy[i] + posVecz[i]*posVecz[i] );

    omega = sqrt(rCrossV_x*rCrossV_x + rCrossV_y*rCrossV_y + rCrossV_z*rCrossV_z ) / (magPosVec*magPosVec);
    vOmega = cbrt( omega );

    amp = amp0 * vOmega * vOmega;

    LNhx = posVecy[i] * momVecz[i] - posVecz[i] * momVecy[i];
    LNhy = posVecz[i] * momVecx[i] - posVecx[i] * momVecz[i];
    LNhz = posVecx[i] * momVecy[i] - posVecy[i] * momVecx[i];

    magL = sqrt(LNhx*LNhx + LNhy*LNhy + LNhz*LNhz);

    LNhx = LNhx / magL;
    LNhy = LNhy / magL;
    LNhz = LNhz / magL;

    alpha = atan2( LNhy, LNhx );

    printf( "alpha = %.16e, omega = %.16e, LNhz = %.16e, vphi = %.16e\n", alpha, omega, LNhz, vphi[i] );
 
    hPlusTS->data->data[i]  = - 0.5 * amp * cos( 2.*vphi[i]) * cos(2.*alpha) * (1. + LNhz*LNhz) 
                            + amp * sin(2.*vphi[i]) * sin(2.*alpha)*LNhz;

    hCrossTS->data->data[i] = - 0.5 * amp * cos( 2.*vphi[i]) * sin(2.*alpha) * (1. + LNhz*LNhz)
                            - amp * sin(2.*vphi[i]) * cos(2.*alpha) * LNhz;

  }

  /* Point the output pointers to the relevant time series and return */
  (*hplus)  = hPlusTS;
  (*hcross) = hCrossTS;


  return XLAL_SUCCESS;
}
