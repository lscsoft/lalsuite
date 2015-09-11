    /*
*  Copyright (C) 2011 Craig Robinson, Enrico Barausse, Yi Pan, 
*                2014 Prayush Kumar (Precessing EOB)
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
 * @addtogroup LALSimIMRSpinAlignedEOB_c
 *
 * @author Craig Robinson, Yi Pan, Prayush Kumar, Stas Babak, Andrea Taracchini
 *
 * @brief Functions for producing SEOBNRv1 waveforms for
 * spinning binaries, as described in
 * Taracchini et al. ( PRD 86, 024011 (2012), arXiv 1202.0790 ).
 * All equation numbers in this file refer to equations of this paper,
 * unless otherwise specified.
 * 
 * @review SEOBNRv1 has been reviewd by Riccardo Sturani, B. Sathyaprakash and
 * Prayush Kumar.
 * The review concluded in fall 2012.  
 *
 * \brief Functions for producing SEOBNRv2 waveforms for
 * spinning binaries, as described in
 * Taracchini et al. ( arXiv 1311.2544 ).
 * 
 * @review SEOBNRv2 has been reviewed by Riccardo Sturani, Prayush Kumar and
 * Stas Babak.
 * The review concluded with git hash 
 * 5bc6bb861de2eb72ca403b9e0f529d83080490fe (August 2014).
 *
 * \brief Functions for producing SEOBNRv3 waveforms for 
 * precessing binaries of spinning compact objects, as described
 * in Pan et al. ( arXiv 1307.6232 ) == YPP.
 */

#include <math.h>
#include <complex.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/LALAdaptiveRungeKutta4.h>
#include <lal/SphericalHarmonics.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>

#include "LALSimIMREOBNRv2.h"
#include "LALSimIMRSpinEOB.h"
#include "LALSimInspiralPrecess.h"
#include "LALSimBlackHoleRingdown.h"
#include "LALSimFindAttachTime.h"

/* Include all the static function files we need */
#include "LALSimIMREOBNQCCorrection.c"
#include "LALSimInspiraldEnergyFlux.c"
#include "LALSimIMREOBNewtonianMultipole.c"
#include "LALSimIMREOBFactorizedWaveform.c"
#include "LALSimIMREOBHybridRingdown.c"
#include "LALSimIMRSpinEOBAuxFuncs.c"
#include "LALSimIMRSpinEOBFactorizedWaveformCoefficients.c"
#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBFactorizedWaveform.c"
#include "LALSimIMRSpinEOBFactorizedFlux.c"
#include "LALSimIMRSpinEOBHcapNumericalDerivative.c"
#include "LALSimIMRSpinAlignedEOBHcapDerivative.c"
#include "LALSimIMRSpinEOBInitialConditions.c"


#define debugOutput 0

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

static double f_alphadotcosi( double x, void * inparams )
{
	PrecEulerAnglesIntegration* params = (PrecEulerAnglesIntegration*) inparams;
	
	REAL8 alphadot = gsl_spline_eval_deriv( params->alpha_spline, x, params->alpha_acc );
	REAL8 beta = gsl_spline_eval( params->beta_spline, x, params->beta_acc );
	
	return -1. * alphadot * cos(beta);	
	
}

static int UNUSED
XLALEOBSpinStopCondition(double UNUSED t,
                           const double values[],
                           double dvalues[],
                           void *funcParams
                          )
{

  SpinEOBParams *params = (SpinEOBParams *)funcParams;
  REAL8 omega, omega_xyz[3];
  double r2;

  r2 = values[0]*values[0] + values[1]*values[1] + values[2]*values[2];
  cross_product( values, dvalues, omega_xyz );
  omega = sqrt(inner_product( omega_xyz, omega_xyz )) / r2;

  /* Terminate when omega reaches peak, and separation is < 6M */
  //if ( omega < params->eobParams->omega )
  if ( r2 < 36. && omega < params->eobParams->omega )
  {
    return 1;
  }

  params->eobParams->omega = omega;
  return GSL_SUCCESS;
}

UNUSED static int
XLALEOBSpinStopConditionBasedOnPR(double UNUSED t,
                           const double values[],
                           double dvalues[],
                           void UNUSED *funcParams
                          )
{
  int debugPK = 0; int debugPKverbose = 0;
  INT4 i;
  SpinEOBParams UNUSED *params = (SpinEOBParams *)funcParams;
  
  REAL8 r2, pDotr = 0;
  REAL8 p[3], r[3], pdotVec[3], rdotVec[3];
  REAL8 omega, omega_xyz[3], L[3], dLdt1[3], dLdt2[3];

  memcpy( r, values, 3*sizeof(REAL8));
  memcpy( p, values+3, 3*sizeof(REAL8));
  memcpy( rdotVec, dvalues, 3*sizeof(REAL8));
  memcpy( pdotVec, dvalues+3, 3*sizeof(REAL8));

  r2 = inner_product(r,r);
  cross_product( values, dvalues, omega_xyz );
  omega = sqrt(inner_product( omega_xyz, omega_xyz )) / r2;
  pDotr = inner_product( p, r ) / sqrt(r2);
  
  REAL8 rdot;
  rdot = inner_product(rdotVec, r) / sqrt(r2);
  double prDot = - inner_product( p, r )*rdot/r2
                + inner_product( pdotVec, r )/sqrt(r2)
                + inner_product( rdotVec, p )/sqrt(r2);

  cross_product( r, pdotVec, dLdt1 );
  cross_product( rdotVec, p, dLdt2 );
  cross_product( r, p, L );

  /* ********************************************************** */
  /* *******  Different termination conditions Follow  ******** */  
  /* ********************************************************** */
  
  /* Terminate if any derivative is Nan */
  for( i = 0; i < 12; i++ )
  {
	  if( isnan(dvalues[i]) || isnan(values[i]) )
	  {          
		  if(debugPK){printf("\n  isnan reached. r2 = %f\n", r2); fflush(NULL); }
          XLALPrintError( "XLAL Error - %s: nan reached at r2 = %f \n", __func__, r2);
          XLAL_ERROR( XLAL_EINVAL );

		  return 1;
	  }
  }
  
  /* ********************************************************** */
  /* *******  Unphysical orbital conditions  ******** */  
  /* ********************************************************** */
  
  /* Terminate if p_r points outwards */  
  if ( r2 < 16 && pDotr >= 0  )
  {
    if(debugPK){
      printf("\n Integration stopping, p_r pointing outwards -- out-spiraling!\n");
      fflush(NULL); 
    }
    return 1;
  }

  /* Terminate if rdot is >0 (OUTspiraling) */  
  if ( r2 < 16 && rdot >= 0  )
  {
    if(debugPK){ 
      printf("\n Integration stopping, dr/dt>0 -- out-spiraling!\n"); 
      fflush(NULL);
    }
    return 1;
  }
  
  /* Terminate if dp_R/dt > 0, i.e. radial momentum is increasing */
  if(r2 < 4. && prDot > 0. )
  {
    if(debugPK){
      printf("\n Integration stopping as prDot = %lf at r = %lf\n",
            prDot, sqrt(r2));
      fflush(NULL);
    }
    return 1;
  }
    
//  /* Terminate if dL/dt > 0, i.e. angular momentum is increasing */
//  if(r2 < 16. && LMagdot > 0. )
//  {
//    if(debugPK){
//        printf("\n Integration stopping as d|L|/dt = %lf at r = %lf\n",
//                LMagdot, sqrt(r2));
//        fflush(NULL);
//    }
//    return 1;
//  }
    
//    if(r2 < 16. && ( dvalues[12] < 0. || dvalues[13] < 0.)  ) {
//        if(debugPK)printf("\n Integration stopping, psiDot<0.01\n");
//        return 1;
//    }

  /* **************************************************************** */
  /*                         Omega related                            */
  /* **************************************************************** */
  
  /* Terminate when omega reaches peak, and separation is < 4M */
  if ( r2 < 16. && omega < params->eobParams->omega )
    params->eobParams->omegaPeaked = 1;

  /* If omega has gone through a second extremum, break */
  if ( r2 < 4. && params->eobParams->omegaPeaked == 1
                && omega > params->eobParams->omega ) 
  {
    if(debugPK) {
      printf("\n Integration stopping, omega reached second extremum\n");
      fflush(NULL);
    }
    return 1;
  }
  
  /* If Momega did not evolve above 0.01 even though r < 4 or omega<0.14 for r<2, break */
  if((r2 < 16. && omega < 0.04) || (r2 < 4. && omega < 0.14 && params->eobParams->omegaPeaked == 1  ) )
  {
    if(debugPK){ 
      printf("\n Integration stopping for omega below threshold, omega=%f at r = %f\n", omega, sqrt(r2));
      fflush(NULL);
    }
    return 1;
  }

    if(r2 < 16. && omega > 1. )
    {
        if(debugPK){
            printf("\n Integration stopping, omega>1 at r = %f\n", sqrt(r2));
            fflush(NULL);
        }
        return 1;
    }
  params->eobParams->omega = omega;
  
  /* **************************************************************** */
  /*              related to Numerical values of x/p/derivatives      */
  /* **************************************************************** */

  /* If momentum derivatives are too large numerically, break */
  if ( r2 < 25 && (fabs(dvalues[3]) > 10 || fabs(dvalues[4]) > 10 || fabs(dvalues[5]) > 10) )
  {
    if(debugPK){
      printf("\n Integration stopping, dpdt > 10 -- too large!\n");
      fflush(NULL);}
    return 1;
  }
  
  /* If p_\Phi is too large numerically, break */
  if( r2 < 16. && values[5] > 10 )
  {
    if(debugPK){ 
      printf("Integration stopping, Pphi > 10 now\n\n"); 
      fflush(NULL); 
    }
    return 1;
  }

  /* **************************************************************** */
  /*              Last resort conditions                              */
  /* **************************************************************** */
  
//  if(r2 < 1.747*1.747) {
//    if(debugPK) {
//      printf("\n Integration stopping, r<2M\n");
//      fflush(NULL);
//    }
//    return 1;
//  }
//  
  /* Very verbose output */
  if(debugPKverbose && r2 < 16.) {
//    printf("\n Integration continuing, r = %f, omega = %f, rdot = %f, pr = %f, prDot = %f\n", 
//      sqrt(r2), omega, rdot, pDotr, prDot);
    printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
      t, values[0], values[1], values[2], values[3], values[4], values[5],
      values[6], values[7], values[8], values[9], values[10], values[11], 
      values[12], values[13], omega);
//    printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", 
//    dvalues[0], dvalues[1], dvalues[2], dvalues[3], dvalues[4], dvalues[5],
//    dvalues[6],  dvalues[7], dvalues[8], dvalues[9], dvalues[10], dvalues[11],
//    dvalues[12], dvalues[13]);
  }
  
  return GSL_SUCCESS;
}


/**
 * Stopping condition for the regular resolution EOB orbital evolution
 * -- stop when reaching max orbital frequency in strong field.
 * At each test,
 * if omega starts to decrease, return 1 to stop evolution;
 * if not, update omega with current value and return GSL_SUCCESS to continue evolution.
 */
static int
XLALEOBSpinAlignedStopCondition(double UNUSED t,  /**< UNUSED */
                           const double values[], /**< dynamical variable values */
                           double dvalues[],      /**< dynamical variable time derivative values */
                           void *funcParams       /**< physical parameters */
                          )
{

  REAL8 omega, r;
  SpinEOBParams *params = (SpinEOBParams *)funcParams;

  r     = values[0];
  omega = dvalues[1];

  //if ( omega < params->eobParams->omega )
  if ( r < 6. && omega < params->eobParams->omega )
  {
    return 1;
  }

  params->eobParams->omega = omega;
  return GSL_SUCCESS;
}

/**
 * Stopping condition for the high resolution EOB orbital evolution
 * -- stop when reaching a minimum radius 0.3M out of the EOB horizon (Eqs. 9b, 37)
 * or when getting nan in any of the four ODE equations
 * At each test,
 * if conditions met, return 1 to stop evolution;
 * if not, return GSL_SUCCESS to continue evolution.
 */
static int
XLALSpinAlignedHiSRStopCondition(double UNUSED t,  /**< UNUSED */
                           const double UNUSED values[], /**< dynamical variable values */
                           double dvalues[],      /**< dynamical variable time derivative values */
                           void UNUSED *funcParams       /**< physical parameters */
                          )
{
  if ( dvalues[2] >= 0. || isnan( dvalues[3] ) || isnan (dvalues[2]) || isnan (dvalues[1]) || isnan (dvalues[0]) )
  {
      //XLALPrintError( "XLAL Error - %s: nan in dvalues \n", __func__);
      //XLAL_ERROR( XLAL_EINVAL );
      return 1;
  }
  return GSL_SUCCESS;
}


/**
 * This function returns the frequency at which the peak amplitude occurs
 * in SEOBNRv(x)
 *
 */
double XLALSimIMRSpinAlignedEOBPeakFrequency(
    REAL8 m1SI,                 /**< mass of companion 1 (kg) */
    REAL8 m2SI,                 /**< mass of companion 2 (kg) */
    const REAL8 spin1z,         /**< z-component of the dimensionless spin of object 1 */
    const REAL8 spin2z,         /**< z-component of the dimensionless spin of object 2 */
    UINT4 SpinAlignedEOBversion /**< 1 for SEOBNRv1, 2 for SEOBNRv2 */
    )
{

  /* The return variable, frequency in Hz */
  double retFreq;

  REAL8 nrOmega;

  /* The mode to use; currently, only l=2, m=2 is supported */
  INT4  ll = 2;
  INT4  mm = 2;

  /* Mass parameters (we'll work in units of solar mass) */
  REAL8  m1 = m1SI / LAL_MSUN_SI;
  REAL8  m2 = m2SI / LAL_MSUN_SI;
  REAL8  Mtotal = m1+m2;
  REAL8  eta = m1*m2 / (Mtotal * Mtotal);

  /* We need spin vectors for the SigmaKerr function */
  REAL8Vector *sigmaKerr = XLALCreateREAL8Vector( 3 );
  int ii;
  REAL8  aa;

  REAL8Vector s1Vec, s2Vec;
  s1Vec.length = s2Vec.length = 3;
  REAL8 spin1[3] = {0., 0., spin1z};
  REAL8 spin2[3] = {0., 0., spin2z};
  s1Vec.data = spin1;
  s2Vec.data = spin2;
  /* the SigmaKerr function uses spins, not chis */
  for( ii = 0; ii < 3; ii++ )
  {
    s1Vec.data[ii] *= m1*m1;
    s2Vec.data[ii] *= m2*m2;
  }

  /* Calculate the normalized spin of the deformed-Kerr background */ 
  if ( XLALSimIMRSpinEOBCalculateSigmaKerr( sigmaKerr, m1, m2, &s1Vec, &s2Vec ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLAL_ERROR( XLAL_EFUNC );
  }

  aa = sigmaKerr->data[2];

 /* Now get the frequency at the peak amplitude */
  switch ( SpinAlignedEOBversion )
  {
   case 1:
     nrOmega = GetNRSpinPeakOmega( ll, mm, eta, aa );
     break;
   case 2:
     nrOmega = GetNRSpinPeakOmegav2( ll, mm, eta, aa );
     break;
   default:
     XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
     XLAL_ERROR( XLAL_EINVAL );
     break;
  }

  retFreq = nrOmega / (2 * LAL_PI * Mtotal * LAL_MTSUN_SI);

  /* Free memory */
  XLALDestroyREAL8Vector( sigmaKerr );

  return retFreq;
}

/**
 * This function generates spin-aligned SEOBNRv1 waveforms h+ and hx.
 * Currently, only the h22 harmonic is available.
 * STEP 0) Prepare parameters, including pre-computed coefficients
 * for EOB Hamiltonian, flux and waveform
 * STEP 1) Solve for initial conditions
 * STEP 2) Evolve EOB trajectory until reaching the peak of orbital frequency
 * STEP 3) Step back in time by tStepBack and volve EOB trajectory again
 * using high sampling rate, stop at 0.3M out of the "EOB horizon".
 * STEP 4) Locate the peak of orbital frequency for NQC and QNM calculations
 * STEP 5) Calculate NQC correction using hi-sampling data
 * STEP 6) Calculate QNM excitation coefficients using hi-sampling data
 * STEP 7) Generate full inspiral waveform using desired sampling frequency
 * STEP 8) Generate full IMR modes -- attaching ringdown to inspiral
 * STEP 9) Generate full IMR hp and hx waveforms
 */
int XLALSimIMRSpinAlignedEOBWaveform(
        REAL8TimeSeries **hplus,     /**<< OUTPUT, +-polarization waveform */
        REAL8TimeSeries **hcross,    /**<< OUTPUT, x-polarization waveform */
        const REAL8     phiC,        /**<< coalescence orbital phase (rad) */ 
        REAL8           deltaT,      /**<< sampling time step */
        const REAL8     m1SI,        /**<< mass-1 in SI unit */ 
        const REAL8     m2SI,        /**<< mass-2 in SI unit */
        const REAL8     fMin,        /**<< starting frequency (Hz) */
        const REAL8     r,           /**<< distance in SI unit */
        const REAL8     inc,         /**<< inclination angle */
        const REAL8     spin1z,      /**<< z-component of spin-1, dimensionless */
        const REAL8     spin2z,       /**<< z-component of spin-2, dimensionless */
        UINT4           SpinAlignedEOBversion /**<< 1 for SEOBNRv1, 2 for SEOBNRv2 */
     )
{
  /* If the EOB version flag is neither 1 nor 2, exit */
  if (SpinAlignedEOBversion != 1 && SpinAlignedEOBversion != 2)
  {
    XLALPrintError("XLAL Error - %s: SEOBNR version flag incorrectly set to %u\n",
        __func__, SpinAlignedEOBversion);
    XLAL_ERROR( XLAL_EERR );
  }

  Approximant SpinAlignedEOBapproximant = (SpinAlignedEOBversion == 1) ? SEOBNRv1 : SEOBNRv2;

  /* 
   * Check spins
   */
  if ( spin1z < -1.0 || spin2z < -1.0 )
  {
    XLALPrintError( "XLAL Error - %s: Component spin less than -1!\n", __func__);
    XLAL_ERROR( XLAL_EINVAL );
  }
  /* If either spin > 0.6, model not available, exit */
  if ( SpinAlignedEOBversion == 1 && ( spin1z > 0.6 || spin2z > 0.6 ) )
  {
    XLALPrintError( "XLAL Error - %s: Component spin larger than 0.6!\nSEOBNRv1 is only available for spins in the range -1 < a/M < 0.6.\n", __func__);
    XLAL_ERROR( XLAL_EINVAL );
  }
 /* For v2 the upper bound is 0.99 */
  if ( SpinAlignedEOBversion == 2 && ( spin1z > 0.99 || spin2z > 0.99 )) 
  {
    XLALPrintError( "XLAL Error - %s: Component spin larger than 0.99!\nSEOBNRv2 is only available for spins in the range -1 < a/M < 0.99.\n", __func__);
    XLAL_ERROR( XLAL_EINVAL );
  }

  INT4 i;

  REAL8Vector *values = NULL;

  /* EOB spin vectors used in the Hamiltonian */
  REAL8Vector *sigmaStar = NULL;
  REAL8Vector *sigmaKerr = NULL;
  REAL8       a, tplspin;
  REAL8       chiS, chiA;

  /* Wrapper spin vectors used to calculate sigmas */
  REAL8Vector s1Vec, s1VecOverMtMt;
  REAL8Vector s2Vec, s2VecOverMtMt;
  REAL8       spin1[3] = {0, 0, spin1z};
  REAL8       spin2[3] = {0, 0, spin2z};
  REAL8       s1Data[3], s2Data[3], s1DataNorm[3], s2DataNorm[3];

  /* Parameters of the system */
  REAL8 m1, m2, mTotal, eta, mTScaled;
  REAL8 amp0;
  REAL8 sSub = 0.0;
  LIGOTimeGPS tc = LIGOTIMEGPSZERO;

  /* Dynamics of the system */
  REAL8Vector rVec, phiVec, prVec, pPhiVec;
  REAL8       omega, v, ham;

  /* Cartesian vectors needed to calculate Hamiltonian */
  REAL8Vector cartPosVec, cartMomVec;
  REAL8       cartPosData[3], cartMomData[3];

  /* Signal mode */
  COMPLEX16   hLM;
  REAL8Vector *sigReVec = NULL, *sigImVec = NULL;

  /* Non-quasicircular correction */
  EOBNonQCCoeffs nqcCoeffs;
  COMPLEX16      hNQC;
  REAL8Vector    *ampNQC = NULL, *phaseNQC = NULL;

  /* Ringdown freq used to check the sample rate */
  COMPLEX16Vector modefreqVec;
  COMPLEX16      modeFreq;

  /* Spin-weighted spherical harmonics */
  COMPLEX16  MultSphHarmP;
  COMPLEX16  MultSphHarmM;

  /* We will have to switch to a high sample rate for ringdown attachment */
  REAL8 deltaTHigh;
  UINT4 resampFac;
  UINT4 resampPwr;
  REAL8 resampEstimate;

  /* How far will we have to step back to attach the ringdown? */
  REAL8 tStepBack;
  INT4  nStepBack;

  /* Dynamics and details of the high sample rate part used to attach the ringdown */
  UINT4 hiSRndx;
  REAL8Vector timeHi, rHi, phiHi, prHi, pPhiHi;
  REAL8Vector *sigReHi = NULL, *sigImHi = NULL;
  REAL8Vector *omegaHi = NULL;

  /* Indices of peak frequency and final point */
  /* Needed to attach ringdown at the appropriate point */
  UINT4 peakIdx = 0, finalIdx = 0;

  /* (2,2) and (2,-2) spherical harmonics needed in (h+,hx) */
  REAL8 y_1, y_2, z1, z2;

  /* Variables for the integrator */
  LALAdaptiveRungeKutta4Integrator       *integrator = NULL;
  REAL8Array              *dynamics   = NULL;
  REAL8Array              *dynamicsHi = NULL;
  INT4                    retLen;
  REAL8  UNUSED           tMax;

  /* Accuracies of adaptive Runge-Kutta integrator */
  const REAL8 EPS_ABS = 1.0e-10;
  const REAL8 EPS_REL = 1.0e-9;

  /*
   * STEP 0) Prepare parameters, including pre-computed coefficients 
   *         for EOB Hamiltonian, flux and waveform  
   */

  /* Parameter structures containing important parameters for the model */
  SpinEOBParams           seobParams;
  SpinEOBHCoeffs          seobCoeffs;
  EOBParams               eobParams;
  FacWaveformCoeffs       hCoeffs;
  NewtonMultipolePrefixes prefixes;

  /* Initialize parameters */
  m1 = m1SI / LAL_MSUN_SI;
  m2 = m2SI / LAL_MSUN_SI;
  mTotal = m1 + m2;
  mTScaled = mTotal * LAL_MTSUN_SI;
  eta    = m1 * m2 / (mTotal*mTotal);
    
    /* For v2 the upper bound is mass ratio 100 */
    if ( SpinAlignedEOBversion == 2 && eta < 100./101./101.)
    {
        XLALPrintError( "XLAL Error - %s: Mass ratio larger than 100!\nSEOBNRv2 is only available for mass ratios up to 100.\n", __func__);
        XLAL_ERROR( XLAL_EINVAL );
    }

  amp0 = mTotal * LAL_MRSUN_SI / r;

  if (pow(LAL_PI*fMin*mTScaled,-2./3.) < 10.0)
  {
    XLAL_PRINT_WARNING("Waveform generation may fail due to high starting frequency. The starting frequency corresponds to a small initial radius of %.2fM. We recommend a lower starting frequency that corresponds to an estimated starting radius > 10M.", pow(LAL_PI*fMin*mTScaled,-2.0/3.0));
  }
 
  /* TODO: Insert potentially necessary checks on the arguments */

  /* Calculate the time we will need to step back for ringdown */
  tStepBack = 100. * mTScaled;
  nStepBack = ceil( tStepBack / deltaT );

  /* Calculate the resample factor for attaching the ringdown */
  /* We want it to be a power of 2 */
  /* If deltaT > Mtot/50, reduce deltaT by the smallest power of two for which deltaT < Mtot/50 */
  resampEstimate = 50. * deltaT / mTScaled;
  resampFac = 1;
  //resampFac = 1 << (UINT4)ceil(log2(resampEstimate));
  
  if ( resampEstimate > 1. )
  {
    resampPwr = (UINT4)ceil( log2( resampEstimate ) );
    while ( resampPwr-- )
    {
      resampFac *= 2u;
    }
  }


  /* Allocate the values vector to contain the initial conditions */
  /* Since we have aligned spins, we can use the 4-d vector as in the non-spin case */
  if ( !(values = XLALCreateREAL8Vector( 4 )) )
  {
    XLAL_ERROR( XLAL_ENOMEM );
  }
  memset ( values->data, 0, values->length * sizeof( REAL8 ));

  /* Set up structures and calculate necessary PN parameters */
  /* Unlike the general case, we only need to calculate these once */
  memset( &seobParams, 0, sizeof(seobParams) );
  memset( &seobCoeffs, 0, sizeof(seobCoeffs) );
  memset( &eobParams, 0, sizeof(eobParams) );
  memset( &hCoeffs, 0, sizeof( hCoeffs ) );
  memset( &prefixes, 0, sizeof( prefixes ) );

  /* Before calculating everything else, check sample freq is high enough */
  modefreqVec.length = 1;
  modefreqVec.data   = &modeFreq;

  if ( XLALSimIMREOBGenerateQNMFreqV2( &modefreqVec, m1, m2, spin1, spin2, 2, 2, 1, SpinAlignedEOBapproximant ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* If Nyquist freq < 220 QNM freq, exit */
  if ( deltaT > LAL_PI / creal(modeFreq) )
  {
    XLALPrintError( "XLAL Error - %s: Ringdown frequency > Nyquist frequency!\nAt present this situation is not supported.\n", __func__);
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EINVAL );
  }

  if ( !(sigmaStar = XLALCreateREAL8Vector( 3 )) )
  {
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_ENOMEM );
  }

  if ( !(sigmaKerr = XLALCreateREAL8Vector( 3 )) )
  {
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_ENOMEM );
  }

  seobParams.alignedSpins = 1;
  seobParams.tortoise     = 1;
  seobParams.sigmaStar    = sigmaStar;
  seobParams.sigmaKerr    = sigmaKerr;
  seobParams.seobCoeffs   = &seobCoeffs;
  seobParams.eobParams    = &eobParams;
  seobParams.nqcCoeffs    = &nqcCoeffs;
  eobParams.hCoeffs       = &hCoeffs;
  eobParams.prefixes      = &prefixes;

  eobParams.m1  = m1;
  eobParams.m2  = m2;
  eobParams.eta = eta;

  s1Vec.length = s2Vec.length = 3;
  s1VecOverMtMt.length = s2VecOverMtMt.length = 3;
  s1Vec.data   = s1Data;
  s2Vec.data   = s2Data;
  s1VecOverMtMt.data   = s1DataNorm;
  s2VecOverMtMt.data   = s2DataNorm;

  /* copy the spins into the appropriate vectors, and scale them by the mass */
  memcpy( s1Data, spin1, sizeof( s1Data ) );
  memcpy( s2Data, spin2, sizeof( s2Data ) );
  memcpy( s1DataNorm, spin1, sizeof( s1DataNorm ) );
  memcpy( s2DataNorm, spin2, sizeof( s2DataNorm ) );

  /* Calculate chiS and chiA */

  chiS = 0.5 * (spin1[2] + spin2[2]);
  chiA = 0.5 * (spin1[2] - spin2[2]);

  for( i = 0; i < 3; i++ )
  {
    s1Data[i] *= m1*m1;
    s2Data[i] *= m2*m2;
  }
 for ( i = 0; i < 3; i++ )
  {
    s1DataNorm[i] = s1Data[i]/mTotal/mTotal;
    s2DataNorm[i] = s2Data[i]/mTotal/mTotal;
  }
  seobParams.s1Vec    = &s1VecOverMtMt;
  seobParams.s2Vec    = &s2VecOverMtMt;

  cartPosVec.length = cartMomVec.length = 3;
  cartPosVec.data = cartPosData;
  cartMomVec.data = cartMomData;
  memset( cartPosData, 0, sizeof( cartPosData ) );
  memset( cartMomData, 0, sizeof( cartMomData ) );

  /* Populate the initial structures */
  if ( XLALSimIMRSpinEOBCalculateSigmaStar( sigmaStar, m1, m2, &s1Vec, &s2Vec ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  if ( XLALSimIMRSpinEOBCalculateSigmaKerr( sigmaKerr, m1, m2, &s1Vec, &s2Vec ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Calculate the value of a */
  /* XXX I am assuming that, since spins are aligned, it is okay to just use the z component XXX */
  /* TODO: Check this is actually the way it works in LAL */
  switch ( SpinAlignedEOBversion )
  {
     case 1:
       tplspin = 0.0;
       break;
     case 2:
       tplspin = (1.-2.*eta) * chiS + (m1 - m2)/(m1 + m2) * chiA;
       break;
     default:
       XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
       XLAL_ERROR( XLAL_EINVAL );
       break;
  }
  /*for ( i = 0; i < 3; i++ )
  {
    a += sigmaKerr->data[i]*sigmaKerr->data[i];
  }
  a = sqrt( a );*/
  seobParams.a = a = sigmaKerr->data[2];
  seobParams.chi1 = spin1[2];
  seobParams.chi2 = spin2[2];

  /* Now compute the spinning H coefficients and store them in seobCoeffs */
  if ( XLALSimIMRCalculateSpinEOBHCoeffs( &seobCoeffs, eta, a, SpinAlignedEOBversion ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  if ( XLALSimIMREOBCalcSpinFacWaveformCoefficients( &hCoeffs, m1, m2, eta, tplspin, chiS, chiA, SpinAlignedEOBversion ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  if ( XLALSimIMREOBComputeNewtonMultipolePrefixes( &prefixes, eobParams.m1, eobParams.m2 )
         == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  switch ( SpinAlignedEOBversion )
  {
     case 1:
       if ( XLALSimIMRGetEOBCalibratedSpinNQC( &nqcCoeffs, 2, 2, eta, a ) == XLAL_FAILURE )
       {
         XLAL_ERROR( XLAL_EFUNC );
       }
       break;
     case 2:
       if ( XLALSimIMRGetEOBCalibratedSpinNQC3D( &nqcCoeffs, 2, 2, m1, m2, a, chiA ) == XLAL_FAILURE )
       {
         XLAL_ERROR( XLAL_EFUNC );
       }
       break;
     default:
       XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
       XLAL_ERROR( XLAL_EINVAL );
       break;
  }

  /*
   * STEP 1) Solve for initial conditions
   */

  /* Set the initial conditions. For now we use the generic case */
  /* Can be simplified if spin-aligned initial conditions solver available. The cost of generic code is negligible though. */
  REAL8Vector *tmpValues = XLALCreateREAL8Vector( 14 );
  if ( !tmpValues )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_ENOMEM );
  }

  memset( tmpValues->data, 0, tmpValues->length * sizeof( REAL8 ) );

  /* We set inc zero here to make it easier to go from Cartesian to spherical coords */
  /* No problem setting inc to zero in solving spin-aligned initial conditions. */
  /* inc is not zero in generating the final h+ and hx */
  if ( XLALSimIMRSpinEOBInitialConditionsV2( tmpValues, m1, m2, fMin, 0, s1Data, s2Data, &seobParams ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( tmpValues );
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /*fprintf( stderr, "ICs = %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", tmpValues->data[0], tmpValues->data[1], tmpValues->data[2],
      tmpValues->data[3], tmpValues->data[4], tmpValues->data[5], tmpValues->data[6], tmpValues->data[7], tmpValues->data[8],
      tmpValues->data[9], tmpValues->data[10], tmpValues->data[11] );*/

  /* Taken from Andrea's code */
/*  memset( tmpValues->data, 0, tmpValues->length*sizeof(tmpValues->data[0]));*/
#if 0
  tmpValues->data[0] = 4.84796591370541;
  tmpValues->data[3] = -0.0186210227887131;
  tmpValues->data[4] = 0.5499544739926464;//tmpValues->data[0]; // q=1
  fprintf( stderr, "ICs = %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", tmpValues->data[0], tmpValues->data[1], tmpValues->data[2],
      tmpValues->data[3], tmpValues->data[4], tmpValues->data[5], tmpValues->data[6], tmpValues->data[7], tmpValues->data[8],
      tmpValues->data[9], tmpValues->data[10], tmpValues->data[11] );

#endif
#if 0
  tmpValues->data[0] = 19.9982539582;
  tmpValues->data[3] = -0.000390702473305;
  tmpValues->data[4] = 4.71107185264/tmpValues->data[0]; // q=1, chi1=chi2=0.98
#endif
#if 0
  tmpValues->data[0] = 19.996332305;
  tmpValues->data[3] = -0.000176807206312;
  tmpValues->data[4] = 4.84719922687/tmpValues->data[0]; // q=8
#endif
#if 0
  tmpValues->data[0] = 6.22645094958;
  tmpValues->data[3] = -0.00851784427559;
  tmpValues->data[4] = 3.09156589713/tmpValues->data[0]; // q=8 chi1=0.5 TEST DYNAMICS
#endif
#if 0
  tmpValues->data[0] = 19.9996712714;
  tmpValues->data[3] = -0.00016532905477;
  tmpValues->data[4] = 4.77661989696/tmpValues->data[0]; // q=8 chi1=0.5
#endif

  /* Now convert to Spherical */
  /* The initial conditions code returns Cartesian components of four vectors x, p, S1 and S2,
   * in the special case that the binary starts on the x-axis and the two spins are aligned
   * with the orbital angular momentum along the z-axis.
   * Therefore, in spherical coordinates the initial conditions are
   * r = x; phi = 0.; pr = px; pphi = r * py.
   */
  values->data[0] = tmpValues->data[0];
  values->data[1] = 0.;
  values->data[2] = tmpValues->data[3];
  values->data[3] = tmpValues->data[0] * tmpValues->data[4];

  //fprintf( stderr, "Spherical initial conditions: %e %e %e %e\n", values->data[0], values->data[1], values->data[2], values->data[3] );

  /*
   * STEP 2) Evolve EOB trajectory until reaching the peak of orbital frequency
   */

  /* Now we have the initial conditions, we can initialize the adaptive integrator */
  if (!(integrator = XLALAdaptiveRungeKutta4Init(4, XLALSpinAlignedHcapDerivative, XLALEOBSpinAlignedStopCondition, EPS_ABS, EPS_REL)))
  {
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  integrator->stopontestonly = 1;
  integrator->retries = 1;

  retLen = XLALAdaptiveRungeKutta4( integrator, &seobParams, values->data, 0., 20./mTScaled, deltaT/mTScaled, &dynamics );
  if ( retLen == XLAL_FAILURE || dynamics == NULL )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Set up pointers to the dynamics */
  rVec.length = phiVec.length = prVec.length = pPhiVec.length = retLen;
  rVec.data    = dynamics->data+retLen;
  phiVec.data  = dynamics->data+2*retLen;
  prVec.data   = dynamics->data+3*retLen;
  pPhiVec.data = dynamics->data+4*retLen;

  //printf( "We think we hit the peak at time %e\n", dynamics->data[retLen-1] );

  /* TODO : Insert high sampling rate / ringdown here */
  #if debugOutput
  FILE *out = fopen( "saDynamics.dat", "w" );
  for ( i = 0; i < retLen; i++ )
  {
    fprintf( out, "%.16e %.16e %.16e %.16e %.16e\n", dynamics->data[i], rVec.data[i], phiVec.data[i], prVec.data[i], pPhiVec.data[i] );
  }
  fclose( out );
  #endif

  if (tStepBack > retLen*deltaT)
  {
    tStepBack = 0.5*retLen*deltaT; //YPnote: if 100M of step back > actual time of evolution, step back 50% of the later
    nStepBack = ceil( tStepBack / deltaT );
  }
 
  /*
   * STEP 3) Step back in time by tStepBack and volve EOB trajectory again 
   *         using high sampling rate, stop at 0.3M out of the "EOB horizon".
   */

  /* Set up the high sample rate integration */
  hiSRndx = retLen - nStepBack;
  deltaTHigh = deltaT / (REAL8)resampFac;

  #if debugOutput
  fprintf( stderr, "Stepping back %d points - we expect %d points at high SR\n", nStepBack, nStepBack*resampFac );
  fprintf( stderr, "Commencing high SR integration... from %.16e %.16e %.16e %.16e %.16e\n",
     (dynamics->data)[hiSRndx],rVec.data[hiSRndx], phiVec.data[hiSRndx], prVec.data[hiSRndx], pPhiVec.data[hiSRndx] );
  #endif

  values->data[0] = rVec.data[hiSRndx];
  values->data[1] = phiVec.data[hiSRndx];
  values->data[2] = prVec.data[hiSRndx];
  values->data[3] = pPhiVec.data[hiSRndx];
  /* For HiSR evolution, we stop at a radius 0.3M from the deformed Kerr singularity, 
   * or when any derivative of Hamiltonian becomes nan */
  integrator->stop = XLALSpinAlignedHiSRStopCondition;

  retLen = XLALAdaptiveRungeKutta4( integrator, &seobParams, values->data, 0., 20./mTScaled, deltaTHigh/mTScaled, &dynamicsHi );
  if ( retLen == XLAL_FAILURE || dynamicsHi == NULL )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  //fprintf( stderr, "We got %d points at high SR\n", retLen );

  /* Set up pointers to the dynamics */
  rHi.length = phiHi.length = prHi.length = pPhiHi.length = timeHi.length = retLen;
  timeHi.data = dynamicsHi->data;
  rHi.data    = dynamicsHi->data+retLen;
  phiHi.data  = dynamicsHi->data+2*retLen;
  prHi.data   = dynamicsHi->data+3*retLen;
  pPhiHi.data = dynamicsHi->data+4*retLen;

  #if debugOutput
  out = fopen( "saDynamicsHi.dat", "w" );
  for ( i = 0; i < retLen; i++ )
  {
    fprintf( out, "%.16e %.16e %.16e %.16e %.16e\n", timeHi.data[i], rHi.data[i], phiHi.data[i], prHi.data[i], pPhiHi.data[i] );
  }
  fclose( out );
  #endif

  /* Allocate the high sample rate vectors */
  sigReHi  = XLALCreateREAL8Vector( retLen + (UINT4)ceil( 20 / ( cimag(modeFreq) * deltaTHigh )) );
  sigImHi  = XLALCreateREAL8Vector( retLen + (UINT4)ceil( 20 / ( cimag(modeFreq) * deltaTHigh )) );
  omegaHi  = XLALCreateREAL8Vector( retLen + (UINT4)ceil( 20 / ( cimag(modeFreq) * deltaTHigh )) );
  ampNQC   = XLALCreateREAL8Vector( retLen );
  phaseNQC = XLALCreateREAL8Vector( retLen );

  if ( !sigReHi || !sigImHi || !omegaHi || !ampNQC || !phaseNQC )
  {
    XLAL_ERROR( XLAL_ENOMEM );
  }

  memset( sigReHi->data, 0, sigReHi->length * sizeof( sigReHi->data[0] ));
  memset( sigImHi->data, 0, sigImHi->length * sizeof( sigImHi->data[0] ));

  /* Populate the high SR waveform */
  REAL8 omegaOld = 0.0;
  INT4  phaseCounter = 0;
  for ( i = 0; i < retLen; i++ )
  {
    values->data[0] = rHi.data[i];
    values->data[1] = phiHi.data[i];
    values->data[2] = prHi.data[i];
    values->data[3] = pPhiHi.data[i];

    omega = XLALSimIMRSpinAlignedEOBCalcOmega( values->data, &seobParams );
    if (omega < 1.0e-15) omega = 1.0e-9; //YPnote: make sure omega>0 during very-late evolution when numerical errors are huge.
    omegaHi->data[i] = omega;            //YPnote: omega<0 is extremely rare and had only happenned after relevant time interval.  
    v = cbrt( omega );

    /* Calculate the value of the Hamiltonian */
    cartPosVec.data[0] = values->data[0];
    cartMomVec.data[0] = values->data[2];
    cartMomVec.data[1] = values->data[3] / values->data[0];

    ham = XLALSimIMRSpinEOBHamiltonian( eta, &cartPosVec, &cartMomVec, &s1VecOverMtMt, &s2VecOverMtMt, sigmaKerr, sigmaStar, seobParams.tortoise, &seobCoeffs );

    if ( XLALSimIMRSpinEOBGetSpinFactorizedWaveform( &hLM, values, v, ham, 2, 2, &seobParams )
           == XLAL_FAILURE )
    {
      /* TODO: Clean-up */
      XLAL_ERROR( XLAL_EFUNC );
    }

    ampNQC->data[i]  = cabs( hLM );
    sigReHi->data[i] = (REAL4)(amp0 * creal(hLM));
    sigImHi->data[i] = (REAL4)(amp0 * cimag(hLM));
    phaseNQC->data[i]= carg( hLM ) + phaseCounter * LAL_TWOPI;

    if ( i && phaseNQC->data[i] > phaseNQC->data[i-1] )
    {
      phaseCounter--;
      phaseNQC->data[i] -= LAL_TWOPI;
    }

    if ( omega <= omegaOld && !peakIdx )
    {
      //printf( "Have we got the peak? omegaOld = %.16e, omega = %.16e\n", omegaOld, omega );
      peakIdx = i;
    }
    omegaOld = omega;
  }
  //printf( "We now think the peak is at %d\n", peakIdx );
  finalIdx = retLen - 1;

  /*
   * STEP 4) Locate the peak of orbital frequency for NQC and QNM calculations
   */

  /* Stuff to find the actual peak time */
  gsl_spline    *spline = NULL;
  gsl_interp_accel *acc = NULL;
  REAL8 omegaDeriv1; //, omegaDeriv2;
  REAL8 time1, time2;
  REAL8 timePeak, timewavePeak = 0., omegaDerivMid;
  REAL8 sigAmpSqHi = 0., oldsigAmpSqHi = 0.;
  INT4  peakCount = 0;

  spline = gsl_spline_alloc( gsl_interp_cspline, retLen );
  acc    = gsl_interp_accel_alloc();

  time1 = dynamicsHi->data[peakIdx];

  gsl_spline_init( spline, dynamicsHi->data, omegaHi->data, retLen );
  omegaDeriv1 = gsl_spline_eval_deriv( spline, time1, acc );
  if ( omegaDeriv1 > 0. )
  {
    time2 = dynamicsHi->data[peakIdx+1];
    //omegaDeriv2 = gsl_spline_eval_deriv( spline, time2, acc );
  }
  else
  {
    //omegaDeriv2 = omegaDeriv1;
    time2 = time1;
    time1 = dynamicsHi->data[peakIdx-1];
    peakIdx--;
    omegaDeriv1 = gsl_spline_eval_deriv( spline, time1, acc );
  }

  do
  {
    timePeak = ( time1 + time2 ) / 2.;
    omegaDerivMid = gsl_spline_eval_deriv( spline, timePeak, acc );

    if ( omegaDerivMid * omegaDeriv1 < 0.0 )
    {
      //omegaDeriv2 = omegaDerivMid;
      time2 = timePeak;
    }
    else
    {
      omegaDeriv1 = omegaDerivMid;
      time1 = timePeak;
    }
  }
  while ( time2 - time1 > 1.0e-5 );

  /*gsl_spline_free( spline );
  gsl_interp_accel_free( acc );
  */

  //XLALPrintInfo( "Estimation of the peak is now at time %.16e\n", timePeak );

  /* Having located the peak of orbital frequency, we set time and phase of coalescence */
  XLALGPSAdd( &tc, -mTScaled * (dynamics->data[hiSRndx] + timePeak));
  gsl_spline_init( spline, dynamicsHi->data, phiHi.data, retLen );
  sSub = gsl_spline_eval( spline, timePeak, acc ) - phiC;
  gsl_spline_free( spline );
  gsl_interp_accel_free( acc );
  /* Apply phiC to hi-sampling waveforms */
  REAL8 thisReHi, thisImHi;
  REAL8 csSub2 = cos(2.0 * sSub);
  REAL8 ssSub2 = sin(2.0 * sSub);
  for ( i = 0; i < retLen; i++)
  {
    thisReHi = sigReHi->data[i];
    thisImHi = sigImHi->data[i];
    sigReHi->data[i] =   thisReHi * csSub2 - thisImHi * ssSub2;
    sigImHi->data[i] =   thisReHi * ssSub2 + thisImHi * csSub2; 
  }

  /*
   * STEP 5) Calculate NQC correction using hi-sampling data
   */

  /* Calculate nonspin and amplitude NQC coefficients from fits and interpolation table */
  /*switch ( SpinAlignedEOBversion )
  {
     case 1:
       if ( XLALSimIMRGetEOBCalibratedSpinNQC( &nqcCoeffs, 2, 2, eta, a ) == XLAL_FAILURE )
       {
         XLAL_ERROR( XLAL_EFUNC );
       }
       break;
     case 2:
       if ( XLALSimIMRGetEOBCalibratedSpinNQC3D( &nqcCoeffs, 2, 2, eta, a, chiA ) == XLAL_FAILURE )
       {
         XLAL_ERROR( XLAL_EFUNC );
       }
       break;
     default:
       XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
       XLAL_ERROR( XLAL_EINVAL );
       break;
  }*/

  /* Calculate phase NQC coefficients */
  if ( XLALSimIMRSpinEOBCalculateNQCCoefficients( ampNQC, phaseNQC, &rHi, &prHi, omegaHi,
          2, 2, timePeak, deltaTHigh/mTScaled, m1, m2, a, chiA, chiS, &nqcCoeffs, SpinAlignedEOBversion ) == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Calculate the time of amplitude peak. Despite the name, this is in fact the shift in peak time from peak orb freq time */
  switch ( SpinAlignedEOBversion )
  {
     case 1:
     timewavePeak = XLALSimIMREOBGetNRSpinPeakDeltaT(2, 2, eta,  a);
       break;
     case 2:
     timewavePeak = XLALSimIMREOBGetNRSpinPeakDeltaTv2(2, 2, m1, m2, spin1z, spin2z );
       break;
     default:
       XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
       XLAL_ERROR( XLAL_EINVAL );
       break;
  }

  /* Apply to the high sampled part */
  #if debugOutput
  out = fopen( "saWavesHi.dat", "w" );
  #endif 
  for ( i = 0; i < retLen; i++ )
  {
    values->data[0] = rHi.data[i];
    values->data[1] = phiHi.data[i] - sSub;
    values->data[2] = prHi.data[i];
    values->data[3] = pPhiHi.data[i];

    if ( XLALSimIMREOBNonQCCorrection( &hNQC, values, omegaHi->data[i], &nqcCoeffs ) == XLAL_FAILURE )
    {
      XLAL_ERROR( XLAL_EFUNC );
    }

    hLM = sigReHi->data[i];
    hLM += I * sigImHi->data[i];
    #if debugOutput
    fprintf( out, "%.16e %.16e %.16e %.16e %.16e\n", timeHi.data[i], creal(hLM), cimag(hLM), creal(hNQC), cimag(hNQC) );
    #endif

    hLM *= hNQC;
    sigReHi->data[i] = (REAL4) creal(hLM);
    sigImHi->data[i] = (REAL4) cimag(hLM);
    sigAmpSqHi = creal(hLM)*creal(hLM)+cimag(hLM)*cimag(hLM);
    if (sigAmpSqHi < oldsigAmpSqHi && peakCount == 0 && (i-1)*deltaTHigh/mTScaled < timePeak - timewavePeak) 
    {
      timewavePeak = (i-1)*deltaTHigh/mTScaled;
      peakCount += 1;
    }
    oldsigAmpSqHi = sigAmpSqHi;
  }
  #if debugOutput
  fclose(out);
  printf("NQCs entering hNQC: %f, %f, %f, %f, %f, %f\n", nqcCoeffs.a1, nqcCoeffs.a2,nqcCoeffs.a3, nqcCoeffs.a3S, nqcCoeffs.a4, nqcCoeffs.a5 );
  printf("NQCs entering hNQC: %f, %f, %f, %f\n", nqcCoeffs.b1, nqcCoeffs.b2,nqcCoeffs.b3, nqcCoeffs.b4 );
  #endif
  if (timewavePeak < 1.0e-16 || peakCount == 0)
  {
    //printf("YP::warning: could not locate mode peak, use calibrated time shift of amplitude peak instead.\n");
    /* NOTE: instead of looking for the actual peak, use the calibrated value,    */
    /*       ignoring the error in using interpolated NQC instead of iterated NQC */
    timewavePeak = timePeak - timewavePeak;
  }
  
  /*
   * STEP 6) Calculate QNM excitation coefficients using hi-sampling data
   */

  /*out = fopen( "saInspWaveHi.dat", "w" );
  for ( i = 0; i < retLen; i++ )
  {
    fprintf( out, "%.16e %.16e %.16e\n", timeHi.data[i], sigReHi->data[i], sigImHi->data[i] );
  }
  fclose( out );*/
  
  /* Attach the ringdown at the time of amplitude peak */
  REAL8 combSize = 7.5; /* Eq. 34 */
  REAL8 chi = (spin1[2] + spin2[2]) / 2. + ((spin1[2] - spin2[2]) / 2.) * ((m1 - m2)/(m1+m2)) / (1. - 2. * eta);

  /* Modify the combsize for SEOBNRv2 */
  /* If chi1=chi2=0, comb = 11. if chi < 0.8, comb = 12. if chi >= 0.8, comb =
   * 13.5 */
  if( SpinAlignedEOBversion == 2 )
  {
    combSize = (spin1[2] == 0. && spin2[2] == 0.) ? 11. : (( eta > 10./121. && chi >= 0.8 ) ? 8.5 : 12.);
    if ( (eta > 30./31./31. && eta <= 10./121. && chi >= 0.8) || (eta <= 30./31./31. && chi >= 0.8 && chi < 0.9) )
      combSize = 13.5;
  }

  REAL8 timeshiftPeak;
  timeshiftPeak = timePeak - timewavePeak;
  if ( SpinAlignedEOBversion == 2)
  {
    timeshiftPeak = (timePeak - timewavePeak) > 0. ? (timePeak - timewavePeak) : 0.;
  }

  /*printf("YP::timePeak and timewavePeak: %.16e and %.16e\n",timePeak,timewavePeak);
  printf("YP::timeshiftPeak and combSize: %.16e and %.16e\n",timeshiftPeak,combSize);
  printf("PK::chi and SpinAlignedEOBversion: %.16e and %u\n\n", chi,SpinAlignedEOBversion);*/

  REAL8Vector *rdMatchPoint = XLALCreateREAL8Vector( 3 );
  if ( !rdMatchPoint )
  {
    XLAL_ERROR( XLAL_ENOMEM );
  }

  if ( combSize > timePeak - timeshiftPeak )
  {
    XLALPrintError( "The comb size looks to be too big!!!\n" );
  }

  rdMatchPoint->data[0] = combSize < timePeak - timeshiftPeak ? timePeak - timeshiftPeak - combSize : 0;
  rdMatchPoint->data[1] = timePeak - timeshiftPeak;
  rdMatchPoint->data[2] = dynamicsHi->data[finalIdx];
  #if debugOutput
  printf("YP::comb range: %f, %f\n",rdMatchPoint->data[0],rdMatchPoint->data[1]);
  #endif
  rdMatchPoint->data[0] -= fmod( rdMatchPoint->data[0], deltaTHigh/mTScaled );
  rdMatchPoint->data[1] -= fmod( rdMatchPoint->data[1], deltaTHigh/mTScaled );
  if ( XLALSimIMREOBHybridAttachRingdown( sigReHi, sigImHi, 2, 2,
              deltaTHigh, m1, m2, spin1[0], spin1[1], spin1[2], spin2[0], spin2[1], spin2[2],
              &timeHi, rdMatchPoint, SpinAlignedEOBapproximant, 1.0 )
          == XLAL_FAILURE ) 
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  /*
   * STEP 7) Generate full inspiral waveform using desired sampling frequency
   */

  /* Now create vectors at the correct sample rate, and compile the complete waveform */
  sigReVec = XLALCreateREAL8Vector( rVec.length + ceil( sigReHi->length / resampFac ) );
  sigImVec = XLALCreateREAL8Vector( sigReVec->length );

  memset( sigReVec->data, 0, sigReVec->length * sizeof( REAL8 ) );
  memset( sigImVec->data, 0, sigImVec->length * sizeof( REAL8 ) );
 
  /* Generate full inspiral waveform using desired sampling frequency */
  /* TODO - Check vectors were allocated */
  for ( i = 0; i < (INT4)rVec.length; i++ )
  {
    values->data[0] = rVec.data[i];
    values->data[1] = phiVec.data[i] - sSub;
    values->data[2] = prVec.data[i];
    values->data[3] = pPhiVec.data[i];

    omega = XLALSimIMRSpinAlignedEOBCalcOmega( values->data, &seobParams );
    v = cbrt( omega );

    /* Calculate the value of the Hamiltonian */
    cartPosVec.data[0] = values->data[0];
    cartMomVec.data[0] = values->data[2];
    cartMomVec.data[1] = values->data[3] / values->data[0];

    ham = XLALSimIMRSpinEOBHamiltonian( eta, &cartPosVec, &cartMomVec, &s1VecOverMtMt, &s2VecOverMtMt, sigmaKerr, sigmaStar, seobParams.tortoise, &seobCoeffs );

    if ( XLALSimIMRSpinEOBGetSpinFactorizedWaveform( &hLM, values, v, ham, 2, 2, &seobParams )
           == XLAL_FAILURE )
    {
      /* TODO: Clean-up */
      XLAL_ERROR( XLAL_EFUNC );
    }

    if ( XLALSimIMREOBNonQCCorrection( &hNQC, values, omega, &nqcCoeffs ) == XLAL_FAILURE )
    {
      XLAL_ERROR( XLAL_EFUNC );
    }

    hLM *= hNQC;

    sigReVec->data[i] = amp0 * creal(hLM);
    sigImVec->data[i] = amp0 * cimag(hLM);
  }

  /*
   * STEP 8) Generate full IMR modes -- attaching ringdown to inspiral
   */

  /* Attach the ringdown part to the inspiral */
  for ( i = 0; i < (INT4)(sigReHi->length / resampFac); i++ )
  {
    sigReVec->data[i+hiSRndx] = sigReHi->data[i*resampFac];
    sigImVec->data[i+hiSRndx] = sigImHi->data[i*resampFac];
  }

  /*
   * STEP 9) Generate full IMR hp and hx waveforms
   */
  
  /* For now, let us just try to create a waveform */
  REAL8TimeSeries *hPlusTS  = XLALCreateREAL8TimeSeries( "H_PLUS", &tc, 0.0, deltaT, &lalStrainUnit, sigReVec->length );
  REAL8TimeSeries *hCrossTS = XLALCreateREAL8TimeSeries( "H_CROSS", &tc, 0.0, deltaT, &lalStrainUnit, sigImVec->length );

  /* TODO change to using XLALSimAddMode function to combine modes */
  /* For now, calculate -2Y22 * h22 + -2Y2-2 * h2-2 directly (all terms complex) */
  /* Compute spin-weighted spherical harmonics and generate waveform */
  REAL8 coa_phase = 0.0;

  MultSphHarmP = XLALSpinWeightedSphericalHarmonic( inc, coa_phase, -2, 2, 2 );
  MultSphHarmM = XLALSpinWeightedSphericalHarmonic( inc, coa_phase, -2, 2, -2 );

  y_1 =   creal(MultSphHarmP) + creal(MultSphHarmM);
  y_2 =   cimag(MultSphHarmM) - cimag(MultSphHarmP);
  z1 = - cimag(MultSphHarmM) - cimag(MultSphHarmP);
  z2 =   creal(MultSphHarmM) - creal(MultSphHarmP);

  for ( i = 0; i < (INT4)sigReVec->length; i++ )
  {
    REAL8 x1 = sigReVec->data[i];
    REAL8 x2 = sigImVec->data[i];

    hPlusTS->data->data[i]  = (x1 * y_1) + (x2 * y_2);
    hCrossTS->data->data[i] = (x1 * z1) + (x2 * z2);
  }

  /* Point the output pointers to the relevant time series and return */
  (*hplus)  = hPlusTS;
  (*hcross) = hCrossTS;

  /* Free memory */
  XLALDestroyREAL8Vector( tmpValues );
  XLALDestroyREAL8Vector( sigmaKerr );
  XLALDestroyREAL8Vector( sigmaStar );
  XLALDestroyREAL8Vector( values );
  XLALDestroyREAL8Vector( rdMatchPoint );
  XLALDestroyREAL8Vector( ampNQC );
  XLALDestroyREAL8Vector( phaseNQC );
  XLALDestroyREAL8Vector( sigReVec );
  XLALDestroyREAL8Vector( sigImVec );
  XLALAdaptiveRungeKutta4Free( integrator );
  XLALDestroyREAL8Array( dynamics );
  XLALDestroyREAL8Array( dynamicsHi );
  XLALDestroyREAL8Vector( sigReHi );
  XLALDestroyREAL8Vector( sigImHi );
  XLALDestroyREAL8Vector( omegaHi );

  return XLAL_SUCCESS;
}








/** ********************************************************************
 *  THE FOLLOWING HAS FUNCTIONS FOR THE PRECESSING EOB MODEL
 *  ********************************************************************
 * */
/**
 * This function generates precessing spinning SEOBNRv3 waveforms h+ and hx.
 * Currently, only h2m harmonics will be generated. 
 * 
 * Input conventions:
 * Cartesian coordinate system: initial \vec{L} is along the z-axis
 * phiC       : in radians
 * deltaT     : in SI units (Hz)
 * m1SI, m2SI : in SI units (kg)
 * fMin       : in SI units (Hz)
 * r          : in SI units (m)
 * inc        : in radians
 * INspin{1,2}: in dimensionless units of m{1,2}^2
 * 
 * Evolution conventions:
 * values[0-2]: r in units of total mass
 * values[3-5]: dimensionless units
 * values[6-8]: 
 * values[9-11]:
 * 
 * STEP 0) Prepare parameters, including pre-computed coefficients
 * for EOB Hamiltonian, flux and waveform
 * STEP 1) Solve for initial conditions
 * STEP 2) Evolve EOB trajectory until reaching the peak of orbital frequency
 * STEP 3) Step back in time by tStepBack and volve EOB trajectory again
 * using high sampling rate, stop at UNDECIDED.
 * STEP 4) Locate the peak of orbital frequency for NQC and QNM calculations
 * STEP 4.5) Rotation of waveform as per the alpha(t) and iota(t)
 * 
 * (Following under construction)
 * STEP 5) Calculate NQC correction using hi-sampling data
 * STEP 6) Calculate QNM excitation coefficients using hi-sampling data
 * STEP 7) Generate full inspiral waveform using desired sampling frequency
 * STEP 8) Generate full IMR modes -- attaching ringdown to inspiral
 * STEP 9) Generate full IMR hp and hx waveforms
 */

int XLALSimIMRSpinEOBWaveform(
        REAL8TimeSeries **hplus,
         REAL8TimeSeries **hcross,
        //LIGOTimeGPS     *tc,
        const REAL8      phiC,
        const REAL8     deltaT,
        const REAL8     m1SI,
        const REAL8     m2SI,
        const REAL8     fMin,
        const REAL8     r,
        const REAL8     inc,
        const REAL8     INspin1[],
        const REAL8     INspin2[]
     )

{
    REAL8Vector   *dynamicsHi = NULL;
    REAL8Vector  *AttachPars = NULL; 
    SphHarmTimeSeries *hlmPTSHi = NULL;
    SphHarmTimeSeries *hlmPTSout = NULL;
    SphHarmTimeSeries *hIMRlmJTSHi = NULL;

    XLALSimIMRSpinEOBWaveformAll(hplus, hcross, &dynamicsHi, &hlmPTSout, &hlmPTSHi, &hIMRlmJTSHi, &AttachPars, 
                        phiC, deltaT, m1SI, m2SI, fMin, r, inc, INspin1[0],
                        INspin1[1], INspin1[2], INspin2[0], INspin2[1],
                        INspin2[2]);


    //int i;
    
    //printf("Stas: checking \n");
    //for (i=0; i<5; i++){
    //    printf("AttachPars: %d, %f \n", i, AttachPars->data[i]);
    //    printf("dyn: %d, %f \n", i, dynamicsHi->data[i]);
    //}
    

    //printf("cleaning memory... \n");
    //if (dynamicsHi == NULL){
    //    printf("dy is already null\n");
    //}else{
    //    printf("Stas: dynamics is cleaned \n");
    //}
    //if (AttachPars == NULL){
    //    printf("att pars is already null\n");
    //}else{
    //    printf("Stas: attach pars is cleaned \n");
    //}
    XLALDestroyREAL8Vector( dynamicsHi );
    
    XLALDestroyREAL8Vector( AttachPars );
    
    XLALDestroySphHarmTimeSeries(hlmPTSout);
    //printf("Stas: Pwave is cleaned \n");

   
    XLALDestroySphHarmTimeSeries(hlmPTSHi);
    //printf("Stas:  J wave is cleaned \n");

    XLALDestroySphHarmTimeSeries(hIMRlmJTSHi);
    //printf("Stas: J wave IMR is cleaned \n");

    return XLAL_SUCCESS;

}


int XLALSimIMRSpinEOBWaveformAll(
        REAL8TimeSeries **hplus,
        REAL8TimeSeries **hcross,
        REAL8Vector      **dynHi, /**<< Here we store and return the seob dynamics for high sampling (end of inspiral) */
        SphHarmTimeSeries **hlmPTSoutput, /**<< Here we store and return the PWave (high sampling) */
        SphHarmTimeSeries **hlmPTSHiOutput, /**<< Here we store and return the JWave (high sampling) */
        SphHarmTimeSeries **hIMRlmJTSHiOutput, /**<< Here we store and return the JWaveIMR (high sampling) */
        REAL8Vector     **AttachPars,   /**<< Parameters of RD attachment: */ 
        //LIGOTimeGPS     *tc,
        const REAL8      phiC,
        const REAL8     deltaT,
        const REAL8     m1SI,
        const REAL8     m2SI,
        const REAL8     fMin,
        const REAL8     r,
        const REAL8     inc,
        const REAL8     INspin1x,
        const REAL8     INspin1y,
        const REAL8     INspin1z,
        const REAL8     INspin2x,
        const REAL8     INspin2y,
        const REAL8     INspin2z
     )

{
  /* TODO: Insert potentially necessary checks on the arguments */

  /* FIXME: Moved this definition out of prototype to allow SWIG interaction*/
  REAL8 INspin1[3], INspin2[3];
  INspin1[0] = INspin1x;
  INspin1[1] = INspin1y;
  INspin1[2] = INspin1z;
  INspin2[0] = INspin2x;
  INspin2[1] = INspin2y;
  INspin2[2] = INspin2z;

  INT4 UNUSED ret;
  INT4 debugPK = 0, debugCustomIC = 0, debugNoNQC = 0;
  INT4 debugRD = 0;
  FILE *out = NULL;
  INT4 i=0;
  INT4 k=0;
  UINT4 j=0;
  INT4 UNUSED status;
  LIGOTimeGPS tc = LIGOTIMEGPSZERO;

  REAL8Vector *AttachParams;
  REAL8Array  *dynamicsHi;
  
  /* SEOBNRv3 model */
  Approximant spinEOBApproximant = SEOBNRv3;
  /* FIXME: The underlying aligned spin EOB model is hard-coded here. 
   * This should be reflected in the name, e.g. SPEOBNRv2 ? */
  INT4 SpinAlignedEOBversion =2;

  /* Vector to store the initial spins */
  REAL8 spin1[3] = {0,0,0}, spin2[3] = {0,0,0};
  memcpy( spin1, INspin1, 3*sizeof(REAL8));
  memcpy( spin2, INspin2, 3*sizeof(REAL8));
  
  /* Check the initial misalignment of component spins, w.r.t. the z-axis
   * theta{1,2}Ini measure the angle w.r.t. the orbital ang. momentum
   * If these angles are < 1.e-3 radians, call the aligned-spin model directly
   * */ 
  REAL8 theta1Ini = 0, theta2Ini = 0;
  REAL8 spin1Norm = -1, spin2Norm = -1;
  if( INspin1[0] != 0 || INspin1[1] != 0 ) {
    spin1Norm = sqrt( INspin1[0]*INspin1[0] + INspin1[1]*INspin1[1] +                          INspin1[2]*INspin1[2] );
    theta1Ini = acos( INspin1[2] / spin1Norm );
  
    if ( theta1Ini <= 1.0e-5) {
      spin1[0] = 0.;
      spin1[1] = 0.;
      spin1[2] = spin1Norm * INspin1[2] / fabs(INspin1[2]);
    }
  }
  if( INspin2[0] != 0 || INspin2[1] != 0 ) {
    spin2Norm = sqrt( INspin2[0]*INspin2[0] + INspin2[1]*INspin2[1] +                          INspin2[2]*INspin2[2] );
    theta2Ini = acos( INspin2[2] / spin2Norm );
  
    if ( theta2Ini <= 1.0e-5) {
      spin2[0] = 0.;
      spin2[1] = 0.;
      spin2[2] = spin2Norm * INspin2[2] / fabs(INspin2[2]);
    }
  }
  if (theta1Ini < 1.0e-3 && theta2Ini < 1.0e-3) {
    ret = XLALSimIMRSpinAlignedEOBWaveform(
          hplus, hcross, phiC, deltaT, m1SI, m2SI, fMin, r, inc,
          spin1Norm, spin2Norm, SpinAlignedEOBversion);
    return ret;
  }
  if ( debugPK ) {
    printf( "theta1Ini, theta2Ini =  %3.10f, %3.10f\n", theta1Ini, theta2Ini );
    printf( "INspin1 = {%3.10f, %3.10f, %3.10f}\n", 
            INspin1[0], INspin1[1], INspin1[2] );
    printf( "INspin2 = {%3.10f, %3.10f, %3.10f}\n", 
            INspin2[0], INspin2[1], INspin2[2] );
    printf( "spin1 = {%3.10f, %3.10f, %3.10f}\n", spin1[0], spin1[1], spin1[2] );
    printf( "spin2 = {%3.10f, %3.10f, %3.10f}\n", spin2[0], spin2[1], spin2[2] );
  }
  
  /* *******************************************************************/
  /* ********************** Memory Allocation **************************/
  /* *******************************************************************/
  
  /* Allocate the values vector to contain the ICs            */
  /* Throughout the code, there are 12 dynamical variables:   */
  /* values[0-2]  - x (Cartesian tortoise separation vector)  */
  /* values[3-5]  - p (Cartesian momentum)                    */
  /* values[6-8]  - spin of body 1                            */
  /* values[9-11] - spin of body 2                            */
  /* Two additional orbital quantities are evolved            */
  /* values[12]   - orbital phase                             */
  /* values[13]   - alpha dot cos(inclination)                */
  REAL8Vector *values = NULL;
  if ( !(values = XLALCreateREAL8Vector( 14 )) )
  {
    XLAL_ERROR(  XLAL_ENOMEM );
  }
  memset( values->data, 0, values->length * sizeof( REAL8 ));

  /* Vector to contain derivatives of ICs */
  REAL8Vector *dvalues = NULL; 
  if ( !(dvalues = XLALCreateREAL8Vector( 14 )) )
  {
    XLAL_ERROR(  XLAL_ENOMEM );
  }
  REAL8       rdotvec[3] = {0,0,0};
 
  /* EOB spin vectors used in the Hamiltonian */
  REAL8       a = 0, tplspin = 0;
  REAL8       chiS = 0, chiA = 0;
  REAL8       spinNQC = 0;
  REAL8Vector *sigmaStar = NULL;
  REAL8Vector *sigmaKerr = NULL;
  if ( !(sigmaStar = XLALCreateREAL8Vector( 3 )) )
  {
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_ENOMEM );
  }
  if ( !(sigmaKerr = XLALCreateREAL8Vector( 3 )) )
  {
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_ENOMEM );
  }
  
  /* Spins not scaled by the mass */
  REAL8 mSpin1[3] = {0,0,0}, mSpin2[3] = {0,0,0};

  /* Wrapper spin vectors used to calculate sigmas */
  REAL8Vector s1Vec, s1VecOverMtMt;
  REAL8Vector s2Vec, s2VecOverMtMt;
  REAL8       s1Data[3]     = {0,0,0}, s2Data[3]     = {0,0,0},
              s1DataNorm[3] = {0,0,0}, s2DataNorm[3] = {0,0,0};

  /* Parameters of the system */
  REAL8 m1 = 0, m2 = 0, mTotal = 0, eta = 0, mTScaled = 0;
  REAL8 amp0 = 0;
  REAL8 UNUSED amp = 0;

  /* Dynamics of the system */
  REAL8Vector tVec, phiDMod, phiMod;
  REAL8Vector posVecx, posVecy, posVecz, momVecx, momVecy, momVecz;
  REAL8Vector s1Vecx, s1Vecy, s1Vecz, s2Vecx, s2Vecy, s2Vecz;
  REAL8       omega = 0, v = 0, ham = 0;

  /* Cartesian vectors needed to calculate Hamiltonian */
  REAL8Vector cartPosVec, cartMomVec;
  REAL8 cartPosData[3] = {0,0,0}, cartMomData[3] = {0,0,0};
  REAL8 rcrossrdotNorm = 0, rvec[3]    = {0,0,0}, rcrossrdot[3] = {0,0,0};
  REAL8 pvec[3] = {0,0,0},  rcrossp[3] = {0,0,0}, rcrosspMag    = 0;
  REAL8 s1dotLN = 0, s2dotLN = 0;
  REAL8 UNUSED s1dotL = 0, s2dotL = 0;

  /* Polar vectors needed for waveform modes calculation */
  REAL8Vector UNUSED polarDynamics, cartDynamics;
  REAL8 polData[4] = {0,0,0,0};  
  
  /* Signal mode */
  COMPLEX16   hLM = 0 + 0j;
  /* COMPLEX16 *sigReVec = NULL, *sigImVec = NULL;*/

  /* Non-quasicircular correction */
  COMPLEX16      hNQC = 0 + 0j;
  EOBNonQCCoeffs nqcCoeffs;
  memset( &nqcCoeffs, 0, sizeof( nqcCoeffs ) );

  /* Ringdown freq used to check the sample rate */
  COMPLEX16Vector  modefreqVec;
  COMPLEX16        modeFreq;

  /* Spin-weighted spherical harmonics */
  COMPLEX16  Y22 = 0 + 0j, Y2m2 = 0 + 0j, Y21 = 0 + 0j, Y2m1 = 0 + 0j;
  COMPLEX16  Y20 = 0 + 0j;
  /* (2,2) and (2,-2) spherical harmonics needed in (h+,hx) */
  /*REAL8  y_1, y_2, z1, z2;*/

  /*Parameters for the portion of the orbit we integrate at high sampling rate*/
  REAL8 deltaTHigh = 0, resampEstimate = 0;
  UINT4 resampFac  = 0, resampPwr      = 0;

  /* How far will we have to step back to attach the ringdown? */
  REAL8 tStepBack = 0, HiSRstart = 0;
  INT4  nStepBack = 0;
  UINT4 hiSRndx   = 0;

  /* Dynamics and details of the high sample rate part used to attach the
   * ringdown */
  REAL8Vector timeHi,    phiDModHi, phiModHi;
  REAL8Vector posVecxHi, posVecyHi, posVeczHi, momVecxHi, momVecyHi, momVeczHi;
  REAL8Vector s1VecxHi,  s1VecyHi,  s1VeczHi,  s2VecxHi,  s2VecyHi,  s2VeczHi;
  REAL8Vector *sigReHi = NULL, *sigImHi = NULL;
  //REAL8Vector * = NULL;

  /* Indices of peak frequency and final point           */
  /* Needed to attach ringdown at the appropriate point  */
  UINT4 UNUSED peakIdx = 0, finalIdx = 0;

  /* Set up structures and calculate necessary PN parameters */
  /* Due to precession, these need to get calculated in every step */
  /* TODO: Only calculate non-spinning parts once */
  SpinEOBParams           seobParams;
  SpinEOBHCoeffs          seobCoeffs;
  EOBParams               eobParams;
  FacWaveformCoeffs       hCoeffs;
  NewtonMultipolePrefixes prefixes;
  
  memset( &seobParams, 0, sizeof(seobParams) );
  memset( &seobCoeffs, 0, sizeof(seobCoeffs) );
  memset( &eobParams,  0, sizeof(eobParams) );
  memset( &hCoeffs,    0, sizeof( hCoeffs ) );
  memset( &prefixes,   0, sizeof( prefixes ) );

  /* Miscellaneous memory for Ringdown attachment        */
  REAL8 tPeakOmega = 0, tAttach = 0, combSize = 0,/*longCombSize,*/ deltaNQC =0;
  REAL8 UNUSED sh  = 0;
  REAL8 vX = 0, vY = 0, vZ = 0;//, rCrossV_x = 0, rCrossV_y = 0, rCrossV_z = 0;
  REAL8 UNUSED vOmega = 0, omegasav = 0, omegasav2 = 0;
  REAL8 magR = 0, Lx   = 0, Ly   = 0, Lz = 0, magL = 0, magJ = 0, 
        LNhx = 0, LNhy = 0, LNhz = 0, magLN =0, Jx = 0, Jy   = 0, Jz = 0;
  REAL8 aI2P = 0, bI2P = 0, gI2P = 0, aP2J = 0, bP2J = 0, gP2J = 0;
  REAL8 chi1J= 0, chi2J= 0, chiJ = 0; 
  REAL8 UNUSED kappaJL = 0;
  REAL8 JLN = 0.0;
  REAL8 JframeEx[3] = {0,0,0}, JframeEy[3] = {0,0,0}, JframeEz[3] = {0,0,0};
  REAL8 LframeEx[3] = {0,0,0}, LframeEy[3] = {0,0,0}, LframeEz[3] = {0,0,0};
  
  /* Variables for the integrator */
  LALAdaptiveRungeKutta4Integrator       *integrator = NULL;
  REAL8Array              *dynamics   = NULL;//, *dynamicsHi = NULL;
  INT4                    retLen = 0, retLenLow = 0, retLenHi = 0;
  INT4                    retLenRDPatch = 0, retLenRDPatchLow = 0;
  //REAL8                   tMax;
  
  /* Interpolation spline */
  gsl_spline    *spline = NULL;
  gsl_interp_accel *acc = NULL;

  /* Accuracies of adaptive Runge-Kutta integrator */
  const REAL8 EPS_ABS = 1.0e-9;
  const REAL8 EPS_REL = 1.0e-8;

  /* Memory for the calculation of the alpha(t) and beta(t) angles */
  REAL8 tmpR[3], tmpRdot[3];
  INT4 phaseCounterA = 0, phaseCounterB = 0;
  
  REAL8Vector *LN_x = NULL;  REAL8Vector *LN_y = NULL; REAL8Vector *LN_z = NULL;
  REAL8Vector *Alpha = NULL; REAL8Vector *Beta = NULL; REAL8Vector *Gamma = NULL;
  REAL8Vector *LN_xHi = NULL; REAL8Vector *LN_yHi = NULL;
  REAL8Vector *LN_zHi = NULL;
  REAL8Vector *AlphaHi = NULL; REAL8Vector *BetaHi = NULL;
  REAL8Vector *GammaHi = NULL;
  
  REAL8 precEulerresult = 0, precEulererror = 0;
  gsl_integration_workspace * precEulerw = gsl_integration_workspace_alloc (1000);
  gsl_function precEulerF;
  PrecEulerAnglesIntegration precEulerparams;

    
  /* Stuff to find the actual peak time */
  //REAL8 omegaDeriv = 0, time1 = 0, time2 = 0;
  UNUSED REAL8  timewavePeak = 0.;
  //REAL8  omegaDerivMid=0.0;
  REAL8 UNUSED sigAmpSqHi = 0., oldsigAmpSqHi = 0.;
  /*INT4   peakCount = 0;*/

  /* Memory to help with rotation of dynamics */
  REAL8 JExnorm = 0;
  SphHarmTimeSeries  *hlmPTS = NULL;
  SphHarmTimeSeries  *hIMRlmJTS = NULL;
  REAL8Sequence *tlist        = NULL;
  REAL8Sequence *tlistRDPatch = NULL;

  SphHarmTimeSeries *hlmPTSHi = NULL;
  SphHarmTimeSeries *hlmPTSout = NULL;
  SphHarmTimeSeries *hIMRlmJTSHi = NULL;
  
  REAL8Sequence *tlistHi        = NULL;
  REAL8Sequence *tlistRDPatchHi = NULL;
  
  /* Memory for ringdown attachment */
  REAL8Vector *rdMatchPoint = XLALCreateREAL8Vector( 3 );  
  REAL8 alJtoI = 0, betJtoI = 0, gamJtoI = 0;
  SphHarmTimeSeries *hIMRlmITS = NULL;
  COMPLEX16 x11 = 0 + 0j;

  /* *******************************************************************/
  /* ********************** Memory Initialization **********************/
  /* *******************************************************************/
  
  /* Initialize mass parameters */
  m1       = m1SI / LAL_MSUN_SI;
  m2       = m2SI / LAL_MSUN_SI;
  mTotal   = m1 + m2;
  mTScaled = mTotal * LAL_MTSUN_SI;
  eta      = m1 * m2 / (mTotal*mTotal);
  
  /* Initialize amplitude scaling parameter */
  amp0 = mTotal * LAL_MRSUN_SI / r;
  //amp0 = 4. * mTotal * LAL_MRSUN_SI * eta / r;
  
  
  if (debugPK) {
    printf("Stas, here is the passes functions\n");
   printf("Inputs: m1 = %.16e, m2 = %.16e, fMin = %.16e, inclination = %.16e\n", 
          m1, m2, (double) fMin, (double) inc );
    printf("Mtotal = %.16e, eta = %.16e \n", mTotal, eta);
    printf("Inputs: spin1 = {%.16e, %.16e, %.16e}\n",  
            spin1[0], spin1[1], spin1[2]);
    printf("Inputs: spin2 = {%.16e, %.16e, %.16e}\n", 
            spin2[0], spin2[1], spin2[2]);
  }  

  /* Calculate the time we will need to step back for ringdown */
  tStepBack = 150. * mTScaled;
  nStepBack = ceil( tStepBack / deltaT );

  /* Calculate the resample factor for attaching the ringdown                 */
  /* We want it to be a power of 2                                            */
  /* If deltaT > Mtot/50, reduce deltaT by the smallest power of two for which*/
  /* deltaT < Mtot/50                                                         */
  resampEstimate = 50. * deltaT / mTScaled;
  resampFac = 1;
  if ( resampEstimate > 1. ) {
    resampPwr = (UINT4)ceil( log2( resampEstimate ) );
    while( resampPwr-- ) { resampFac *= 2u; }
  }

  /* Wrapper spin vectors used to calculate sigmas                 */
  /* Note: spin{1,2} have INspin{1,2} (\chi vectors) at this point */
  s1Vec.length = s2Vec.length = 3;
  s1Vec.data   = s1Data;
  s2Vec.data   = s2Data;
  
  s1VecOverMtMt.length = s2VecOverMtMt.length = 3;
  s1VecOverMtMt.data   = s1DataNorm;
  s2VecOverMtMt.data   = s2DataNorm;

  for( i = 0; i < 3; i++ )
  {
    s1Data[i]     = mSpin1[i] = spin1[i] * m1 * m1;
    s2Data[i]     = mSpin2[i] = spin2[i] * m2 * m2;
    s1DataNorm[i] = s1Data[i] / mTotal / mTotal;
    s2DataNorm[i] = s2Data[i] / mTotal / mTotal;
  }
  
  if (debugPK){
    printf("Stas, here is the passes functions\n");
   printf("Inputs: m1 = %.16e, m2 = %.16e, fMin = %.16e, inclination = %.16e\n", 
            m1, m2, (double) fMin, (double) inc );
    printf("Inputs: spin1 = {%.16e, %.16e, %.16e}\n", 
            spin1[0], spin1[1], spin1[2]);
    printf("Inputs: spin2 = {%.16e, %.16e, %.16e}\n",
            spin2[0], spin2[1], spin2[2]);
    printf("Inputs: s1V = {%.16e, %.16e, %.16e}\n",
            s1Vec.data[0], s1Vec.data[1], s1Vec.data[2]);
    printf("Inputs: s2V = {%.16e, %.16e, %.16e}\n",
            s2Vec.data[0], s2Vec.data[1], s2Vec.data[2]);
  }

  /* ************************************************* */
  /* Populate the initial structures                   */
  /* ************************************************* */
  
  /* ************************************************* */
  /* Spin parameters                                   */
  if ( XLALSimIMRSpinEOBCalculateSigmaStar( 
          sigmaStar, m1, m2, &s1Vec, &s2Vec ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }
  if ( XLALSimIMRSpinEOBCalculateSigmaKerr( 
          sigmaKerr, m1, m2, &s1Vec, &s2Vec ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Calculate the value of a */
  seobParams.a = a = sqrt( inner_product(sigmaKerr->data, sigmaKerr->data) );
  seobParams.s1Vec = &s1VecOverMtMt;
  seobParams.s2Vec = &s2VecOverMtMt;

  if (debugPK){
      printf("Inputs: sigma = {%.16e, %.16e, %.16e}\n", 
        sigmaKerr->data[0], sigmaKerr->data[1], sigmaKerr->data[2]);
      printf("Inputs: star = {%.16e, %.16e, %.16e}\n",  
        sigmaStar->data[0], sigmaStar->data[1], sigmaStar->data[2]);
      printf("Inputs: a = %.16e\n", a);
      fflush(NULL);
  }

  if ( !(AttachParams = XLALCreateREAL8Vector( 5 )) )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_ENOMEM );
  }
  *AttachPars = AttachParams;
  /* ************************************************* */
  /* Cartesian position and momentum vectors           */
  
  /* Cartesian vectors needed to calculate Hamiltonian */
  cartPosVec.length = cartMomVec.length = 3;
  cartPosVec.data = cartPosData;
  cartMomVec.data = cartMomData;
  memset( cartPosData, 0, sizeof( cartPosData ) );
  memset( cartMomData, 0, sizeof( cartMomData ) );
 

  /* ************************************************* */
  /* Waveform parameter structures                     */
  
  /* Spin-EOB parameters */
  seobParams.alignedSpins = 0;
  seobParams.tortoise     = 1;
  seobParams.sigmaStar    = sigmaStar;
  seobParams.sigmaKerr    = sigmaKerr;
  seobParams.seobCoeffs   = &seobCoeffs;
  seobParams.eobParams    = &eobParams;
  seobParams.nqcCoeffs    = &nqcCoeffs;
  /* Non-Spin-EOB parameters */
  eobParams.hCoeffs       = &hCoeffs;
  eobParams.prefixes      = &prefixes;
  seobCoeffs.SpinAlignedEOBversion = SpinAlignedEOBversion;
  eobParams.m1  = m1;
  eobParams.m2  = m2;
  eobParams.eta = eta;

  /* ************************************************* */
  /* Pre-compute the Hamiltonian coefficients          */
  if ( XLALSimIMRCalculateSpinEOBHCoeffs( &seobCoeffs, eta, a,
                          SpinAlignedEOBversion ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Pre-compute the coefficients for the Newtonian factor of hLM */
  if ( XLALSimIMREOBComputeNewtonMultipolePrefixes( &prefixes, eobParams.m1,
			eobParams.m2 ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* ************************************************* */
  /* ***** Get    the INITIAL CONDITIONS               */
  /* ************************************************* */
  
  if( debugPK )
  {
    printf("Calling the XLALSimIMRSpinEOBInitialConditions function!\n");
    printf(
      "Inputs: m1 = %.16e, m2 = %.16e, fMin = %.16e, inclination = %.16e\n",
                      m1, m2, (double) fMin, (double) inc );
    printf("Inputs: mSpin1 = {%.16e, %.16e, %.16e}\n",  
                      mSpin1[0], mSpin1[1], mSpin1[2]);
    printf("Inputs: mSpin2 = {%.16e, %.16e, %.16e}\n",  
                      mSpin2[0], mSpin2[1], mSpin2[2]);
    fflush(NULL);
  }

  REAL8Vector* tmpValues2 = NULL;
  tmpValues2 = XLALCreateREAL8Vector( 14 ); 

  if ( XLALSimIMRSpinEOBInitialConditions( tmpValues2, m1, m2, fMin, inc,
       	mSpin1, mSpin2, &seobParams ) == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }
  tmpValues2->data[12] = 0.;
  tmpValues2->data[13] = 0.;

  if(debugCustomIC) {
//    //v2 0.8 .06
//    tmpValues2->data[0]=  2.4947503693442151e+01;
//    tmpValues2->data[1]= 0.0000000000000000e+00;
//    tmpValues2->data[2]= 0.0000000000000000e+00;
//    tmpValues2->data[3]= -1.1443986771509659e-04;
//    tmpValues2->data[4]= 2.0974358101454393e-01;
//    tmpValues2->data[5]= 0.0000000000000000e+00;
//    tmpValues2->data[6]= 0.0000000000000000e+00;
//    tmpValues2->data[7]=  0.0000000000000000e+00;
//    tmpValues2->data[8]= 5.5555555555555558e-01;
//    tmpValues2->data[9]= 0.0000000000000000e+00;
//    tmpValues2->data[10]= 0.0000000000000000e+00;
//    tmpValues2->data[11]= 1.6666666666666666e-02;
//
//    //v1 0.5 0.5
//    tmpValues2->data[0]=  2.4977926627506903e+01;
//    tmpValues2->data[1]= 0.0000000000000000e+00;
//    tmpValues2->data[2]= 0.0000000000000000e+00;
//    tmpValues2->data[3]= -1.1624761163076212e-04;
//    tmpValues2->data[4]= 2.1082424860474935e-01;
//    tmpValues2->data[5]= 0.0000000000000000e+00;
//    tmpValues2->data[6]= 0.0000000000000000e+00;
//    tmpValues2->data[7]=  0.0000000000000000e+00;
//    tmpValues2->data[8]= 3.4722222222222221e-01;
//    tmpValues2->data[9]= 0.0000000000000000e+00;
//    tmpValues2->data[10]= 0.0000000000000000e+00;
//    tmpValues2->data[11]= 1.3888888888888888e-02;
  
//    v2devel 0.8 0.6
//    tmpValues2->data[0]=  2.4947499317612962e+01;
//    tmpValues2->data[1]= 0.0000000000000000e+00;
//    tmpValues2->data[2]= 0.0000000000000000e+00;
//    tmpValues2->data[3]= -9.7392171179022500e-05;
//    tmpValues2->data[4]= 2.0974360118734747e-01;
//    tmpValues2->data[5]= 0.0000000000000000e+00;
//    tmpValues2->data[6]= 0.0000000000000000e+00;
//    tmpValues2->data[7]=  0.0000000000000000e+00;
//    tmpValues2->data[8]= 5.5555555555555558e-01;
//    tmpValues2->data[9]= 0.0000000000000000e+00;
//    tmpValues2->data[10]= 0.0000000000000000e+00;
//    tmpValues2->data[11]= 1.6666666666666666e-02;
//  v1devel 0.5 0.5
//    tmpValues2->data[0]=  2.4977922258283822e+01;
//    tmpValues2->data[1]= 0.0000000000000000e+00;
//    tmpValues2->data[2]= 0.0000000000000000e+00;
//    tmpValues2->data[3]= -9.8650292751988489e-05;
//    tmpValues2->data[4]= 2.1082427042465196e-01;
//    tmpValues2->data[5]= 0.0000000000000000e+00;
//    tmpValues2->data[6]= 0.0000000000000000e+00;
//    tmpValues2->data[7]=  0.0000000000000000e+00;
//    tmpValues2->data[8]= 3.4722222222222221e-01;
//    tmpValues2->data[9]= 0.0000000000000000e+00;
//    tmpValues2->data[10]= 0.0000000000000000e+00;
//    tmpValues2->data[11]= 1.3888888888888888e-02;
//
//    //  v1devel 0.5 0.5 short
//    tmpValues2->data[0]=  1.4210453570976361e+01;
//    tmpValues2->data[1]= 0.0000000000000000e+00;
//    tmpValues2->data[2]= 0.0000000000000000e+00;
//    tmpValues2->data[3]= -4.9684240542886924e-04;
//    tmpValues2->data[4]= 2.9021826002997486e-01;
//    tmpValues2->data[5]= 0.0000000000000000e+00;
//    tmpValues2->data[6]= 0.0000000000000000e+00;
//    tmpValues2->data[7]=  0.0000000000000000e+00;
//    tmpValues2->data[8]= 3.4722222222222221e-01;
//    tmpValues2->data[9]= 0.0000000000000000e+00;
//    tmpValues2->data[10]= 0.0000000000000000e+00;
//    tmpValues2->data[11]= 1.3888888888888888e-02;
 
    //  v1devel 0.5 0.5 long
//    tmpValues2->data[0]=  3.5925771624136843e+01;
//    tmpValues2->data[1]= 0.0000000000000000e+00;
//    tmpValues2->data[2]= 0.0000000000000000e+00;
//    tmpValues2->data[3]= -3.4725194158512486e-05;
//    tmpValues2->data[4]= 1.7311431867161625e-01;
//    tmpValues2->data[5]= 0.0000000000000000e+00;
//    tmpValues2->data[6]= 0.0000000000000000e+00;
//    tmpValues2->data[7]=  0.0000000000000000e+00;
//    tmpValues2->data[8]= 3.4722222222222221e-01;
//    tmpValues2->data[9]= 0.0000000000000000e+00;
//    tmpValues2->data[10]= 0.0000000000000000e+00;
//    tmpValues2->data[11]= 1.3888888888888888e-02;
//
//    //  v1devel 0.5 0.5 medium
//    tmpValues2->data[0]=  2.2611498964597168e+01;
//    tmpValues2->data[1]= 0.0000000000000000e+00;
//    tmpValues2->data[2]= 0.0000000000000000e+00;
//    tmpValues2->data[3]= -1.3187256837094701e-04;
//    tmpValues2->data[4]= 2.2274218670264384e-01;
//    tmpValues2->data[5]= 0.0000000000000000e+00;
//    tmpValues2->data[6]= 0.0000000000000000e+00;
//    tmpValues2->data[7]=  0.0000000000000000e+00;
//    tmpValues2->data[8]= 3.4722222222222221e-01;
//    tmpValues2->data[9]= 0.0000000000000000e+00;
//    tmpValues2->data[10]= 0.0000000000000000e+00;
//    tmpValues2->data[11]= 1.3888888888888888e-02;


    tmpValues2->data[0]=  1.8656971023349701e+01;
    tmpValues2->data[1]= 0;
    tmpValues2->data[2]= 0;
    tmpValues2->data[3]= -4.8895110550351657e-04;
    tmpValues2->data[4]= 2.4619856470732523e-01;
    tmpValues2->data[5]= -2.6870560011919208e-05;
    tmpValues2->data[6]= 2.9451156124040437e-02;
    tmpValues2->data[7]= -1.1434246556844330e-02;
    tmpValues2->data[8]= 2.3068072871192638e-01;
    tmpValues2->data[9]= 6.2954182743516019e-03;
    tmpValues2->data[10]= 8.8246456372334629e-03;
    tmpValues2->data[11]= 1.5200695012151960e-01;

/*
    tmpValues2->data[0]=  18.692546;
    tmpValues2->data[1]= -3.081570006533801e-23;
    tmpValues2->data[2]= -2.655985321845305e-19;
    tmpValues2->data[3]= -0.0005055333254955916;
    tmpValues2->data[4]= 0.247049159617696;
    tmpValues2->data[5]= -5.732707053941175e-05;
    tmpValues2->data[6]= -2.3036527238752252e-02;
    tmpValues2->data[7]= -3.3429958143641740e-02;
    tmpValues2->data[8]= 2.3933425361333752e-01;
    tmpValues2->data[9]= 7.5408416145236912e-02;
    tmpValues2->data[10]= 1.7433761720072755e-01;
    tmpValues2->data[11]= -2.3113050136189490e-02;
*/

  //tmpValues2->data[3] *= 1.02;
  //tmpValues2->data[4] *= 1.02;
  //tmpValues2->data[5] *= -0.01;
  }
  
  if(debugPK)
  {
    printf("Setting up initial conditions, returned values are:\n");
	  for( j=0; j < tmpValues2->length; j++)
      printf("%.16le\n", tmpValues2->data[j]);
  }

  /* Copy over the IC to values */
  for (j = 0; j < tmpValues2->length; j++)
    values->data[j] = tmpValues2->data[j];
  
  /* ************************************************************************ */
  /* Calculate the values of chiS and chiA, as given in Eq.16 in YPP
   * Assuming \vec{L} to be pointing in the direction of \vec{r}\times\vec{p} */
  
  if(debugPK)printf("\nReached the point where LN is to be calculated\n");

  memset( dvalues->data, 0, 14 * sizeof(REAL8) );
  if( XLALSpinHcapRvecDerivative( 0, values->data, dvalues->data,
          (void*) &seobParams) == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }
  if(debugPK)printf("\nCalculated Rdot\n");

  /* Calculate r cross rDot */ 
  memcpy( rdotvec, dvalues->data, 3*sizeof(REAL8) );
  memcpy( rvec,    values->data,  3*sizeof(REAL8) );
  memcpy( pvec,    values->data+3,3*sizeof(REAL8) );
  
  cross_product( rvec, rdotvec, rcrossrdot );
  rcrossrdotNorm = sqrt(inner_product( rcrossrdot, rcrossrdot ));
  for( i = 0; i < 3; i++ ) { rcrossrdot[i] /= rcrossrdotNorm; }
   
  s1dotLN = inner_product( spin1, rcrossrdot );
  s2dotLN = inner_product( spin2, rcrossrdot );

  /* Why are we calculating s{1,2} dot L ? */
  cross_product( rvec, pvec, rcrossp );
  rcrosspMag = sqrt(inner_product(rcrossp, rcrossp));
  for( i = 0; i < 3; i++ ) { rcrossp[i] /= rcrosspMag; }
    
  s1dotL = inner_product( spin1, rcrossrdot );
  s2dotL = inner_product( spin2, rcrossrdot );

  if(debugPK) {
    printf("rXp = %3.10f %3.10f %3.10f\n", rcrossp[0], rcrossp[1], rcrossp[2]);
    printf("rXrdot = %3.10f %3.10f %3.10f\n", rcrossrdot[0], rcrossrdot[1], rcrossrdot[2]);
    fflush(NULL);
  }

  chiS = 0.5 * (s1dotLN + s2dotLN);
  chiA = 0.5 * (s1dotLN - s2dotLN);

  /* Compute the test-particle limit spin of the deformed-Kerr background */
  switch ( SpinAlignedEOBversion )
  {
    case 1:
      tplspin = 0.0;
      break;
    case 2:
      tplspin = (1.-2.*eta) * chiS + (m1 - m2)/(m1 + m2) * chiA;
      break;
    default:
      XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
      XLAL_ERROR( XLAL_EINVAL );
      break;
  }

  /* ************************************************* */
  /* Populate the Waveform initial structures          */
  /* ************************************************* */
  /* Pre-compute the non-spinning and spinning coefficients for hLM factors */
  if( debugPK ) {
    printf("tplspin = %.12e, chiS = %.12e, chiA = %.12e\n", tplspin, chiS, chiA);
    fflush(NULL);
  }
  
  if ( XLALSimIMREOBCalcSpinFacWaveformCoefficients( &hCoeffs, m1, m2, eta,
        tplspin, chiS, chiA, SpinAlignedEOBversion ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* ************************************************* */
  /* Compute the coefficients for the NQC Corrections  */
  /* ************************************************* */
  /* FIXME check NQC version used */
  spinNQC = (1.-2.*eta) * chiS + (m1 - m2)/(m1 + m2) * chiA;

  //Stas: Check if initial frequency is too high
  //
  REAL8 NRPeakOmega22 = GetNRSpinPeakOmegav2(2, 2, eta, spinNQC) / mTScaled;
  REAL8 freqMinRad = pow(10.0, -1.5)/(LAL_PI*mTScaled);
  REAL8 signOfa = spinNQC/fabs(spinNQC);
  REAL8 spn2 = spinNQC*spinNQC;
  REAL8 Z1 = 1.0 + pow(1.0 - spn2, 1./3.)*(pow(1.0+spinNQC, 1./3.) + pow(1.0-spinNQC, 1./3.) );
  REAL8 Z2 = sqrt(3.0*spn2 + Z1*Z1);
  REAL8 rISCO = 3.0 + Z2 - signOfa*sqrt( (3.-Z1)*(3.+Z1 + 2.*Z2));
  REAL8 fISCO = pow(rISCO, -1.5)/(LAL_PI*mTScaled);

  REAL8 f_star =  freqMinRad;


  if (debugPK){
      printf("Stas - spin = %4.10f \n", spinNQC);
      printf("Stas - NRPeakOmega22 =  %4.10f,   %4.10f \n",  GetNRSpinPeakOmegav2(2, 2, eta, spinNQC) / mTotal,  GetNRSpinPeakOmegav2(2, 2, eta, spinNQC));
      printf("Stas ---- check for fmin NRPeakOmega22 = %4.10f, freqMinRad = %4.10f \n", NRPeakOmega22, freqMinRad);
      printf("Stas -- minf freq is min( %4.10f, %4.10f )\n", NRPeakOmega22*0.1, freqMinRad);
      printf("Stas -- initial radius (apr) %4.10f \n", pow(LAL_PI*fMin*mTScaled,-2./3.) );
      printf("Stas -- omega_orb_0 = %4.10f \n",  rcrossrdotNorm/inner_product( rvec, rvec )/mTScaled);
      printf("Stas -- rISCO  = %4.10f , freq_ISCO = %4.10f \n", rISCO, fISCO);

  }

  if (pow(LAL_PI*fMin*mTScaled,-2./3.) < 10.0)
  {
      XLAL_PRINT_WARNING("Waveform generation may fail due to high starting frequency. The starting frequency corresponds to a small initial radius of %.2fM. We recommend a lower starting frequency that corresponds to an estimated starting radius > 10M.", pow(LAL_PI*fMin*mTScaled,-2.0/3.0));
  }
  if (fMin > f_star){
      XLALPrintError("XLAL Error - %s: Intitial frequency is too high, the limit is %4.10f \n", __func__, f_star);
      XLALDestroyREAL8Vector( sigmaKerr );
      XLALDestroyREAL8Vector( sigmaStar );
      XLALDestroyREAL8Vector( values );
      XLALDestroyREAL8Vector( rdMatchPoint );

      XLAL_ERROR (XLAL_EINVAL );
  } 

  switch ( SpinAlignedEOBversion )
  {
	  case 1:
	    if(debugPK)printf("\t NQC: spin used = %.12e\n", spinNQC);
	    XLALSimIMRGetEOBCalibratedSpinNQC( &nqcCoeffs, 2, 2, eta, spinNQC );
	    break;
	  case 2:
	    if(debugPK)printf("\t NQC: spins used = %.12e, %.12e\n", spinNQC, chiA);
	    XLALSimIMRGetEOBCalibratedSpinNQC3D( &nqcCoeffs, 2, 2, m1, m2, spinNQC, chiA );
	    break;
	  default:
	    XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
	    XLAL_ERROR( XLAL_EINVAL );
	    break;
  }
  
  /* FIXME NOTE the lines below, NQC coeffs are put to zero for debugging */
  if (debugNoNQC) {
    nqcCoeffs.a1 = nqcCoeffs.a2 = nqcCoeffs.a3 = nqcCoeffs.a3S = nqcCoeffs.a4 = 
      nqcCoeffs.a5 = nqcCoeffs.b1 = nqcCoeffs.b2 = nqcCoeffs.b3 = nqcCoeffs.b4 = 
      0;
  }
 
  if ( debugPK )
  {
	  /* Print out all NQC coefficients */
    printf(" NQC: a1 = %.16e, a2 = %.16e,\n a3 = %.16e, a3S = %.16e,\n a4 = %.16e, a5 = %.16e\n b1 = %.16e, b2 = %.16e,\n b3 = %.16e, b4 = %.16e\n",
        nqcCoeffs.a1, nqcCoeffs.a2, nqcCoeffs.a3, 
        nqcCoeffs.a3S, nqcCoeffs.a4, nqcCoeffs.a5, 
        nqcCoeffs.b1, nqcCoeffs.b2, nqcCoeffs.b3, nqcCoeffs.b4 );

	  /* Print out all mass parameters */
	  printf("m1SI = %lf, m2SI = %lf, m1 = %lf, m2 = %lf\n",
			(double) m1SI, (double) m2SI, (double) m1, (double) m2 );
	  printf("mTotal = %lf, mTScaled = %lf, eta = %lf\n",
			(double) mTotal, (double) mTScaled, (double) eta );
	  /* Print out all spin parameters */
	  printf("spin1 = {%lf,%lf,%lf}, spin2 = {%lf,%lf,%lf}\n",
			(double) spin1[0], (double) spin1[1], (double) spin1[2],
			(double) spin2[0], (double) spin2[1], (double) spin2[2]);
	  printf("mSpin1 = {%lf,%lf,%lf}, mSpin2 = {%lf,%lf,%lf}\n",
			(double) mSpin1[0], (double) mSpin1[1], (double) mSpin1[2],
			(double) mSpin2[0], (double) mSpin2[1], (double) mSpin2[2]);
	
      printf("sigmaStar = {%lf,%lf,%lf}, sigmaKerr = {%lf,%lf,%lf}\n",
			(double) sigmaStar->data[0], (double) sigmaStar->data[1],
			(double) sigmaStar->data[2], (double) sigmaKerr->data[0],
			(double) sigmaKerr->data[1], (double) sigmaKerr->data[2]);
	  printf("a = %lf, tplspin = %lf, chiS = %lf, chiA = %lf\n",
			(double) a, (double) tplspin, (double) chiS, (double) chiA);
	  printf("s1Vec = {%lf,%lf,%lf}, s2Vec = {%lf,%lf,%lf}\n",
			(double) seobParams.s1Vec->data[0], (double) seobParams.s1Vec->data[1],
			(double) seobParams.s1Vec->data[2], (double) seobParams.s2Vec->data[0],
			(double) seobParams.s2Vec->data[1], (double) seobParams.s2Vec->data[2]);
	  printf("a is used to compute Hamiltonian coefficients,\n tplspin and chiS and chiA for the multipole coefficients\n");
	  
    printf("s1VecOverM = {%.12e,%.12e,%.12e}\n", (double) s1VecOverMtMt.data[0], 
	  (double) s1VecOverMtMt.data[1], (double) s1VecOverMtMt.data[2]);
	  printf("s2VecOverM = {%.12e,%.12e,%.12e}\n", (double) s2VecOverMtMt.data[0],
	  (double) s2VecOverMtMt.data[1], (double) s2VecOverMtMt.data[2]);
      
    double StasS1 = sqrt(spin1[0]*spin1[0] + spin1[1]*spin1[1] +spin1[2]*spin1[2]);
    double StasS2 = sqrt(spin2[0]*spin2[0] + spin2[1]*spin2[1] +spin2[2]*spin2[2]);
    printf("Stas: amplitude of spin1 = %.16e, amplitude of spin2 = %.16e, theta1 = %.16e , theta2 = %.16e, phi1 = %.16e, phi2 = %.16e  \n", 
            StasS1, StasS2, acos(spin1[2]/StasS1)/LAL_PI, acos(spin2[2]/StasS2)/LAL_PI, 
            atan2(spin1[1], spin1[0])/LAL_PI, atan2(spin2[1], spin2[0])/LAL_PI );
    fflush(NULL);
  }

  /* ***************************************************************** */
  /* ******** Integrate the Hamiltonian Equations ******************** */
  /* ***************************************************************** */
  
  /* ************************************** */
  /* Low Sampling-rate integration          */
  
  /* Initialize the GSL integrator */
  if (!(integrator = XLALAdaptiveRungeKutta4Init(14, 
          XLALSpinHcapNumericalDerivative, XLALEOBSpinStopConditionBasedOnPR,
          EPS_ABS, EPS_REL)))
  {
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Ensure that integration stops ONLY when the stopping condition is True */
  integrator->stopontestonly = 1;
 
  if( debugPK ){
    printf("\n r = {%f,%f,%f}\n",
      values->data[0], values->data[1], values->data[2]);
    fflush(NULL);
  }

  if(debugPK) { printf("\n\n BEGINNING THE EVOLUTION\n\n"); fflush(NULL); }
  
  /* Call the integrator */
  retLen = XLALAdaptiveRungeKutta4( integrator, &seobParams, values->data,
				0., 20./mTScaled, deltaT/mTScaled, &dynamics );
  retLenLow = retLen;
  
  if ( retLen == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }
  
  if(debugPK) { printf("\n\n FINISHED THE EVOLUTION\n\n"); fflush(NULL); }  
  
  /* ************************************** */
  /* Store the dynamics in separate vectors */
  
  tVec.length = posVecx.length = posVecy.length = posVecz.length = 
  momVecx.length = momVecy.length = momVecz.length = 
  s1Vecx.length = s1Vecy.length = s1Vecz.length = 
  s2Vecx.length = s2Vecy.length = s2Vecz.length = 
  phiDMod.length = phiMod.length = retLen;
  
  tVec.data   = dynamics->data;
  posVecx.data = dynamics->data+retLen;
  posVecy.data = dynamics->data+2*retLen;
  posVecz.data = dynamics->data+3*retLen;
  momVecx.data = dynamics->data+4*retLen;
  momVecy.data = dynamics->data+5*retLen;
  momVecz.data = dynamics->data+6*retLen;
  s1Vecx.data = dynamics->data+7*retLen;
  s1Vecy.data = dynamics->data+8*retLen;
  s1Vecz.data = dynamics->data+9*retLen;
  s2Vecx.data = dynamics->data+10*retLen;
  s2Vecy.data = dynamics->data+11*retLen;
  s2Vecz.data = dynamics->data+12*retLen;
  phiDMod.data= dynamics->data+13*retLen;
  phiMod.data = dynamics->data+14*retLen;

  if (debugPK) {
    /* Write the dynamics to file */
    out = fopen( "seobDynamics.dat", "w" );
    for ( i = 0; i < retLen; i++ )
    {
       //YP: output orbital phase and phase modulation separately, instead of their sum
       fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", 
       tVec.data[i], 
       posVecx.data[i], posVecy.data[i], posVecz.data[i], 
       momVecx.data[i], momVecy.data[i], momVecz.data[i],
       s1Vecx.data[i], s1Vecy.data[i], s1Vecz.data[i], 
       s2Vecx.data[i], s2Vecy.data[i], s2Vecz.data[i],
       phiDMod.data[i], phiMod.data[i] );
    }
    fclose( out );
    
    printf("AT We are here %d %d\n",retLen, nStepBack);
    fflush(NULL);
  }
  
  /* ************************************** */
  /* High Sampling-rate integration         */

  /* If 150M of step back > actual time of evolution, step back 50% of that */
  if (tStepBack > retLen*deltaT)
  {
    tStepBack = 0.5*retLen*deltaT; 
    nStepBack = ceil( tStepBack / deltaT );
  }

  /* Step back in time by tStepBack and volve EOB trajectory again 
   * using high sampling rate. */
  hiSRndx    = retLen - nStepBack;
  deltaTHigh = deltaT / (REAL8)resampFac;
  HiSRstart  = tVec.data[hiSRndx];
  
  if (debugPK) {
    printf("AT We are here %d\n",nStepBack);
    printf("Stas: start HighSR at %.16e \n", HiSRstart);
    printf( "Stepping back %d points - we expect %d points at high SR\n", nStepBack, nStepBack*resampFac );
    fflush(NULL);
  }
  
  /* Copy over the dynamics at "hiSRndx" over as "initial values" */
  values->data[0] = posVecx.data[hiSRndx];
  values->data[1] = posVecy.data[hiSRndx];
  values->data[2] = posVecz.data[hiSRndx];
  values->data[3] = momVecx.data[hiSRndx];
  values->data[4] = momVecy.data[hiSRndx];
  values->data[5] = momVecz.data[hiSRndx];
  values->data[6] = s1Vecx.data[hiSRndx];
  values->data[7] = s1Vecy.data[hiSRndx];
  values->data[8] = s1Vecz.data[hiSRndx];
  values->data[9] = s2Vecx.data[hiSRndx];
  values->data[10]= s2Vecy.data[hiSRndx];
  values->data[11]= s2Vecz.data[hiSRndx];
  values->data[12]= phiDMod.data[hiSRndx];
  values->data[13]= phiMod.data[hiSRndx];

  if (debugPK){
    printf( "Commencing high SR integration.. From: \n" );
    for( i=0; i<12; i++)printf("%.16e\n", values->data[i]);
    fflush(NULL);
  }
  
  /* For HiSR evolution, we stop at FIXME */
  integrator->stop = XLALEOBSpinStopConditionBasedOnPR;

  /* Call the integrator */
  retLen = XLALAdaptiveRungeKutta4( integrator, &seobParams, values->data, 
									0., 20./mTScaled, deltaTHigh/mTScaled, &dynamicsHi );
  retLenHi = retLen;

  if ( retLen == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  if(debugPK) { printf( "Finished high SR integration... \n" ); fflush(NULL); }

  /* Set up pointers to the dynamics */
  timeHi.length = posVecxHi.length = posVecyHi.length = posVeczHi.length = 
  momVecxHi.length = momVecyHi.length = momVeczHi.length = 
  s1VecxHi.length = s1VecyHi.length = s1VeczHi.length = 
  s2VecxHi.length = s2VecyHi.length = s2VeczHi.length = 
  phiDModHi.length = phiModHi.length = retLen;

  timeHi.data   = dynamicsHi->data;
  posVecxHi.data = dynamicsHi->data+retLen;
  posVecyHi.data = dynamicsHi->data+2*retLen;
  posVeczHi.data = dynamicsHi->data+3*retLen;
  momVecxHi.data = dynamicsHi->data+4*retLen;
  momVecyHi.data = dynamicsHi->data+5*retLen;
  momVeczHi.data = dynamicsHi->data+6*retLen;
  s1VecxHi.data = dynamicsHi->data+7*retLen;
  s1VecyHi.data = dynamicsHi->data+8*retLen;
  s1VeczHi.data = dynamicsHi->data+9*retLen;
  s2VecxHi.data = dynamicsHi->data+10*retLen;
  s2VecyHi.data = dynamicsHi->data+11*retLen;
  s2VeczHi.data = dynamicsHi->data+12*retLen;
  phiDModHi.data= dynamicsHi->data+13*retLen;
  phiModHi.data = dynamicsHi->data+14*retLen;

  if (debugPK){
    out = fopen( "seobDynamicsHi.dat", "w" );
    for ( i = 0; i < retLen; i++ )
    {
      fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", 
      timeHi.data[i], 
      posVecxHi.data[i], posVecyHi.data[i], posVeczHi.data[i], 
      momVecxHi.data[i], momVecyHi.data[i], momVeczHi.data[i],
      s1VecxHi.data[i], s1VecyHi.data[i], s1VeczHi.data[i], 
      s2VecxHi.data[i], s2VecyHi.data[i], s2VeczHi.data[i],
      phiDModHi.data[i], phiModHi.data[i] );
    }
    fclose( out );
  }

  /* ************************************************************************ */
  /* **************     Compute alpha (t) and beta (t).  ******************** */
  /* ************************************************************************ */
  
  /* Interpolate trajectories to compute L_N (t) in order to get alpha (t) and
   * beta (t). */
  if(debugPK){
    fprintf( stderr, "Generating Alpha and Beta angle timeseries at low SR\n" );
    fflush(NULL);
  }

  LN_x = XLALCreateREAL8Vector( retLenLow );
  LN_y = XLALCreateREAL8Vector( retLenLow );
  LN_z = XLALCreateREAL8Vector( retLenLow );
  
  Alpha = XLALCreateREAL8Vector( retLenLow );
  Beta   = XLALCreateREAL8Vector( retLenLow );
  
  gsl_spline *x_spline = gsl_spline_alloc( gsl_interp_cspline, retLenLow );
  gsl_spline *y_spline = gsl_spline_alloc( gsl_interp_cspline, retLenLow );
  gsl_spline *z_spline = gsl_spline_alloc( gsl_interp_cspline, retLenLow );
  
  gsl_interp_accel *x_acc    = gsl_interp_accel_alloc();
  gsl_interp_accel *y_acc    = gsl_interp_accel_alloc();
  gsl_interp_accel *z_acc    = gsl_interp_accel_alloc();
  
  gsl_spline_init( x_spline, tVec.data, posVecx.data, retLenLow );
  gsl_spline_init( y_spline, tVec.data, posVecy.data, retLenLow );
  gsl_spline_init( z_spline, tVec.data, posVecz.data, retLenLow );
  
  if (debugPK){
    fprintf( stderr, "WRiting Alpha and Beta angle timeseries at low SR to alphaANDbeta.dat\n" );
    fflush(NULL);
    out = fopen( "alphaANDbeta.dat","w");
  }
  for( i=0; i < retLenLow; i++ )
  {
    tmpR[0] = posVecx.data[i]; tmpR[1] = posVecy.data[i]; tmpR[2] = posVecz.data[i];
		tmpRdot[0] = gsl_spline_eval_deriv( x_spline, tVec.data[i], x_acc );
		tmpRdot[1] = gsl_spline_eval_deriv( y_spline, tVec.data[i], y_acc );
		tmpRdot[2] = gsl_spline_eval_deriv( z_spline, tVec.data[i], z_acc );
		
		LN_x->data[i] = tmpR[1] * tmpRdot[2] - tmpR[2] * tmpRdot[1];
		LN_y->data[i] = tmpR[2] * tmpRdot[0] - tmpR[0] * tmpRdot[2];
		LN_z->data[i] = tmpR[0] * tmpRdot[1] - tmpR[1] * tmpRdot[0];
		
		magLN = sqrt(LN_x->data[i] * LN_x->data[i] + LN_y->data[i] * LN_y->data[i] 
              + LN_z->data[i] * LN_z->data[i]);
		LN_x->data[i] /= magLN; LN_y->data[i] /= magLN; LN_z->data[i] /= magLN;
		
		/* Unwrap the two angles */
    if (fabs(LN_x->data[i]) <= 1.e-7 && fabs(LN_y->data[i]) <=1.e-7){
      Alpha->data[i] = 0.0;
    } else {
      Alpha->data[i] = atan2( LN_y->data[i], LN_x->data[i] ) 
                      + phaseCounterA * LAL_TWOPI;
    }
    
    if( i && Alpha->data[i] - Alpha->data[i-1] > 5. )
		{
			phaseCounterA--; 
			Alpha->data[i] -= LAL_TWOPI;
		} else if ( i && Alpha->data[i] - Alpha->data[i-1] < -5. )
		{
			phaseCounterA++;
			Alpha->data[i] += LAL_TWOPI;
		}
    
    /* FIXME: Why is there a "0" multiplying phaseCounterB ?? */
		Beta->data[i] = acos( LN_z->data[i] ) + 0*phaseCounterB * LAL_TWOPI;
    
		if( i && Beta->data[i] > Beta->data[i-1] )
		{
			phaseCounterB--;
			//Beta->data[i] -= LAL_TWOPI;
		}
    
    if(debugPK) {
      fprintf( out, "%.16e %.16e %.16e %d %d %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", 
      tVec.data[i], Alpha->data[i], Beta->data[i], 
      phaseCounterA, phaseCounterB, tmpR[0], tmpR[1], tmpR[2], 
      tmpRdot[0], tmpRdot[1], tmpRdot[2],
      LN_x->data[i], LN_y->data[i], LN_z->data[i] );
    }
	}
  if (debugPK) fclose(out);

  /* Integrate \dot{\alpha} \cos{\beta} to get the final Euler angle*/
  //gsl_spline_free(x_spline);
  //gsl_interp_accel_free(x_acc);
  //x_spline = gsl_spline_alloc( gsl_interp_cspline, retLenLow );
  //x_acc = gsl_interp_accel_alloc();  
  gsl_spline_init( x_spline, tVec.data, Alpha->data, retLenLow );
  
  //y_spline = gsl_spline_alloc( gsl_interp_cspline, retLenLow );
  //y_acc = gsl_interp_accel_alloc();
  gsl_spline_init( y_spline, tVec.data, Beta->data, retLenLow );
  
  Gamma = XLALCreateREAL8Vector( retLenLow );

  precEulerparams.alpha_spline = x_spline;
  precEulerparams.alpha_acc    = x_acc;
  precEulerparams.beta_spline  = y_spline;
  precEulerparams.beta_acc     = y_acc;
   
  precEulerF.function = &f_alphadotcosi;
  precEulerF.params   = &precEulerparams;
   
  if (debugPK) {
    fprintf( stderr,"WRiting Gamma angle timeseries at low SR to gamma.dat\n");
    fflush(NULL);
    out = fopen( "gamma.dat","w");  
  }
  
  for( i = 0; i < retLenLow; i++ )
  {
    if( i==0 ) { Gamma->data[i] = 0.; }
    else
		{
      gsl_integration_qags (&precEulerF, tVec.data[i-1], tVec.data[i], 1e-9, 1e-9, 1000, precEulerw, &precEulerresult, &precEulererror);
			Gamma->data[i] = Gamma->data[i-1] + precEulerresult;
    }
		 			
		if (debugPK) 
      fprintf( out, "%.16e %.16e %.16e %.16e\n", tVec.data[i], Gamma->data[i],
                                              precEulerresult, precEulererror);
  }  
  if (debugPK) fclose(out);

/*
 * STEP 4.2) Interpolate trajectories to compute L_N (t) in order to get alpha (t) and beta (t)
 *
 *fprintf( stderr, "Generating Alpha and Beta angle timeseries at high SR\n" );*/


  LN_xHi = XLALCreateREAL8Vector( retLenHi );
  LN_yHi = XLALCreateREAL8Vector( retLenHi );
  LN_zHi = XLALCreateREAL8Vector( retLenHi );
  
  AlphaHi = XLALCreateREAL8Vector( retLenHi );
  BetaHi   = XLALCreateREAL8Vector( retLenHi );
  
  gsl_spline_free(x_spline);
  gsl_interp_accel_free(x_acc);
  gsl_spline_free(y_spline);
  gsl_interp_accel_free(y_acc);
  gsl_spline_free(z_spline);
  gsl_interp_accel_free(z_acc);
  
  x_spline = gsl_spline_alloc( gsl_interp_cspline, retLenHi );
  y_spline = gsl_spline_alloc( gsl_interp_cspline, retLenHi );
  z_spline = gsl_spline_alloc( gsl_interp_cspline, retLenHi );
  
  x_acc    = gsl_interp_accel_alloc();
  y_acc    = gsl_interp_accel_alloc();
  z_acc    = gsl_interp_accel_alloc();
  
  gsl_spline_init( x_spline, timeHi.data, posVecxHi.data, retLenHi );
  gsl_spline_init( y_spline, timeHi.data, posVecyHi.data, retLenHi );
  gsl_spline_init( z_spline, timeHi.data, posVeczHi.data, retLenHi );
 
  if (debugPK){
    fprintf( stderr, "WRiting Alpha and Beta angle timeseries at High SR to alphaANDbetaHi.dat\n" );
    fflush(NULL);
    out = fopen( "alphaANDbetaHi.dat","w");
  }
  for( i=0; i < retLenHi; i++ )
  {
		tmpR[0] = posVecxHi.data[i]; tmpR[1] = posVecyHi.data[i]; tmpR[2] = posVeczHi.data[i];
		tmpRdot[0] = gsl_spline_eval_deriv( x_spline, timeHi.data[i], x_acc );
		tmpRdot[1] = gsl_spline_eval_deriv( y_spline, timeHi.data[i], y_acc );
		tmpRdot[2] = gsl_spline_eval_deriv( z_spline, timeHi.data[i], z_acc );
		
		LN_xHi->data[i] = tmpR[1] * tmpRdot[2] - tmpR[2] * tmpRdot[1];
		LN_yHi->data[i] = tmpR[2] * tmpRdot[0] - tmpR[0] * tmpRdot[2];
		LN_zHi->data[i] = tmpR[0] * tmpRdot[1] - tmpR[1] * tmpRdot[0];
		
		magLN = sqrt(LN_xHi->data[i] * LN_xHi->data[i] + LN_yHi->data[i] * LN_yHi->data[i] + LN_zHi->data[i] * LN_zHi->data[i]);
		LN_xHi->data[i] /= magLN; LN_yHi->data[i] /= magLN; LN_zHi->data[i] /= magLN;
		
		/* Unwrap the two angles */
    if (fabs(LN_xHi->data[i]) <= 1.e-7 && fabs(LN_yHi->data[i]) <=1.e-7){
      AlphaHi->data[i] = 0.0;
    } else {
      AlphaHi->data[i] = atan2( LN_yHi->data[i], LN_xHi->data[i] ) + phaseCounterA * LAL_TWOPI;
    }
    if( i && AlphaHi->data[i] - AlphaHi->data[i-1] > 5. )
		{
			phaseCounterA--; 
			AlphaHi->data[i] -= LAL_TWOPI;
		}
		else if( i && AlphaHi->data[i] - AlphaHi->data[i-1] < -5. )
		{
			phaseCounterA++;
			AlphaHi->data[i] += LAL_TWOPI;
		}
		
		/* Make sure that Alpha agrees initially with the low SR value */
		AlphaHi->data[i] -= (AlphaHi->data[0] - Alpha->data[hiSRndx]);
		
    /* FIXME: Why is phaseCounterB multiplied by 0?? */
		BetaHi->data[i] = acos( LN_zHi->data[i] ) + 0*phaseCounterB * LAL_TWOPI;
		if( i && BetaHi->data[i] > BetaHi->data[i-1] )
		{
			phaseCounterB--;
			//Beta->data[i] -= LAL_TWOPI;
		}
    
    if (debugPK){
      fprintf( out, "%.16e %.16e %.16e %d %d %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", 
      tVec.data[hiSRndx]+timeHi.data[i], 
		  AlphaHi->data[i], BetaHi->data[i], 
      phaseCounterA, phaseCounterB, 
      tmpR[0], tmpR[1], tmpR[2], tmpRdot[0], tmpRdot[1], tmpRdot[2],
			LN_xHi->data[i], LN_yHi->data[i], LN_zHi->data[i] );
    }
	}
  if (debugPK) fclose(out);

  /* Integrate \dot{\alpha} \cos{\beta} to get the final Euler angle*/
  //x_spline = gsl_spline_alloc( gsl_interp_cspline, retLenHi );
  //x_acc = gsl_interp_accel_alloc();  
  gsl_spline_init( x_spline, timeHi.data, AlphaHi->data, retLenHi );
  
  //y_spline = gsl_spline_alloc( gsl_interp_cspline, retLenHi );
  //y_acc = gsl_interp_accel_alloc();
  gsl_spline_init( y_spline, timeHi.data, BetaHi->data, retLenHi );
  
  GammaHi = XLALCreateREAL8Vector( retLenHi );

  //PrecEulerAnglesIntegration precEulerparams;
  precEulerparams.alpha_spline = x_spline;
  precEulerparams.alpha_acc    = x_acc;
  precEulerparams.beta_spline  = y_spline;
  precEulerparams.beta_acc     = y_acc;
   
  /*gsl_integration_workspace * */
  gsl_integration_workspace_free( precEulerw );
  precEulerw = gsl_integration_workspace_alloc (1000);
  precEulerresult = 0, precEulererror = 0;
   
  //gsl_function precEulerF;
  precEulerF.function = &f_alphadotcosi;
  precEulerF.params = &precEulerparams;
   
  if(debugPK){
    fprintf( stderr, 
      "WRiting Gamma angle timeseries at High SR to gammaHi.dat\n" );
    out = fopen( "gammaHi.dat","w");  
  }
  
  for( i = 0; i < retLenHi; i++ )
  {
    if( i==0 ) { GammaHi->data[i] = Gamma->data[hiSRndx]; }
    else
		{
		 gsl_integration_qags(&precEulerF, timeHi.data[i-1], timeHi.data[i],
        1e-9, 1e-9, 1000, precEulerw, &precEulerresult, &precEulererror);
		 GammaHi->data[i] = GammaHi->data[i-1] + precEulerresult;
		}
		 			
		if (debugPK) 
      fprintf( out, "%.16e %.16e %.16e %.16e\n", 
        tVec.data[hiSRndx]+timeHi.data[i], GammaHi->data[i], 
        precEulerresult, precEulererror);
	}  
  if (debugPK) fclose(out);

  gsl_integration_workspace_free( precEulerw );
  gsl_spline_free( x_spline );
  gsl_spline_free( y_spline );
  gsl_spline_free( z_spline );
  gsl_interp_accel_free( x_acc );
  gsl_interp_accel_free( y_acc );
  gsl_interp_accel_free( z_acc );
   
  /* WaveStep 1.5: moved to here  */
  modefreqVec.length = 1;
  modefreqVec.data   = &modeFreq;
  
  if ( XLALSimIMREOBGenerateQNMFreqV2( &modefreqVec, m1, m2, spin1, spin2, 2, 2, 1, spinEOBApproximant ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  retLenRDPatch = (UINT4)ceil( 40 / ( cimag(modeFreq) * deltaTHigh ));
  retLenRDPatchLow = (UINT4)ceil( 40 / ( cimag(modeFreq) * deltaT ));
  //printf("Stas modeFreq = %f, retLenRDPatch = %d EOB_RD_EFOLDS = %f \n", cimag(modeFreq)*deltaTHigh, retLenRDPatch, EOB_RD_EFOLDS);
  
  /* Allocate the high sample rate vectors */
  //sigReHi  = XLALCreateREAL8Vector( retLen + retLenRDPatch );
  //sigImHi  = XLALCreateREAL8Vector( retLen + retLenRDPatch );
  //omegaHi  = XLALCreateREAL8Vector( retLenHi + retLenRDPatch);        
  //if (!omegaHi )
  //{
    //XLAL_ERROR( XLAL_ENOMEM );
  //}

  /*memset( sigReHi->data, 0, sigReHi->length * sizeof( sigReHi->data[0] ));
  memset( sigImHi->data, 0, sigImHi->length * sizeof( sigImHi->data[0] ));*/



  /* ************************************************************************ */
  /* **************     Waveform Generation  ******************************** */
  /* ************************************************************************ */
  
  /* WaveStep 1
   * Locate merger point (max omega), calculate J, chi and kappa at merger, and construct final J frame
   */
  /* WaveStep 1.1: locate merger point */
    // Find tAmpMax to determin the tAttachment
    REAL8 *radiusVec;
    REAL8 radiusData[retLenHi];
    radiusVec = &radiusData[0];
//    REAL8Vector timeHiV;
//    timeHiV.length = retLenHi;
//    timeHiV.data = dynamicsHi->data;
//    printf("(retLen, retLenHi, dt, tTot, tTotHi)=(%d,%d,%f,%f,%f) \n", retLen, retLenHi, timeHiV.data[1]-timeHiV.data[0], ( timeHiV.data[1]-timeHiV.data[0])*retLen,( timeHiV.data[1]-timeHiV.data[0])*retLenHi);

    for ( i = 0; i < retLenHi; i++ )
    {
        for ( j = 0; j < 3; j++ )
        {
            values->data[j] = dynamicsHi->data[(j+1)*retLenHi + i];
        }
        radiusData[i] = sqrt(values->data[0]*values->data[0] + values->data[1]*values->data[1] + values->data[2]*values->data[2]);
    }
    if (debugPK){
        printf("Andrea the attachment time based on Omega is %f \n", tAttach);
    }

    
    int found = 0;
    if (debugPK) {
        printf("Stas searching for maxima in omega .... \n");
    }
    REAL8 tMaxOmega;
    tPeakOmega = XLALSimLocateOmegaTime(dynamicsHi, values->length, retLenHi, seobParams, seobCoeffs, m1, m2, radiusVec, &found, &tMaxOmega);

    if(tPeakOmega == 0.0 || found==0){
        if (debugPK){
            printf("maximum of omega and A(r) is not found, looking for amplitude maximum\n");
        }
    }
    else{
        if (debugPK){
            printf("The omega-related time is found and it's %f \n", tPeakOmega);
        }
    }
    if(debugPK) {
      printf( "Estimation of the peak is now at time %.16e, %.16e \n", 
            tPeakOmega, tPeakOmega+HiSRstart);
      fflush(NULL);
    }
  
    //exit(0);
  
  /* WaveStep 1.2: calculate J at merger */
  spline = gsl_spline_alloc( gsl_interp_cspline, retLen );
  acc    = gsl_interp_accel_alloc();

  /* Obtain the dynamics at the time where Omega = d\phi/dt peaks */
  for ( j = 0; j < values->length; j++ )
  {
    gsl_spline_init( spline, 
                    dynamicsHi->data, dynamicsHi->data+(j+1)*retLen, retLen );
    values->data[j] = gsl_spline_eval( spline, tPeakOmega, acc );
  }

  /* Calculate dr/dt */
  memset( dvalues->data, 0, 14*sizeof(dvalues->data[0]));
  if( XLALSpinHcapRvecDerivative( 0, values->data, dvalues->data, 
        &seobParams) != XLAL_SUCCESS )
  {
    printf( " Calculation of dr/dt at t = tPeakOmega \n");
    XLAL_ERROR( XLAL_EFUNC );
  }
    
  /* Calculare Omega = r x dr/dt */
  memcpy( rdotvec, dvalues->data, 3*sizeof(REAL8));
  memcpy( rvec,    values->data,  3*sizeof(REAL8));
  memcpy( pvec,    values->data+3,3*sizeof(REAL8));
  cross_product( rvec, rdotvec, rcrossrdot );
  cross_product( rvec, pvec,    rcrossp );
  /* Calculate L at OmegaPeak */
  Lx = rcrossp[0]; Ly = rcrossp[1]; Lz = rcrossp[2];
  magL = sqrt(inner_product(rcrossp, rcrossp));
  /* Calculate J at OmegaPeak */
  Jx = eta*Lx + values->data[6] + values->data[9];
  Jy = eta*Ly + values->data[7] + values->data[10];
  Jz = eta*Lz + values->data[8] + values->data[11];
  magJ = sqrt( Jx*Jx + Jy*Jy + Jz*Jz );
  
  if(debugPK){ 
    printf("J at merger: %e, %e, %e (mag = %e)\n", Jx, Jy, Jz, magJ); 
    fflush(NULL); 
  }

  /* WaveStep 1.3: calculate chi and kappa at merger */
  chi1J = values->data[6]*Jx + values->data[7] *Jy + values->data[8] *Jz;
  chi2J = values->data[9]*Jx + values->data[10]*Jy + values->data[11]*Jz;
  chi1J/= magJ*m1*m1/mTotal/mTotal;
  chi2J/= magJ*m2*m2/mTotal/mTotal;
  /*chiJ = (chi1J+chi2J)/2. + (chi1J-chi2J)/2.*sqrt(1. - 4.*eta)/(1. - 2.*eta);*/
  chiJ = (chi1J+chi2J)/2. + (chi1J-chi2J)/2.*((m1-m2)/(m1+m2))/(1. - 2.*eta);
  kappaJL = (Lx*Jx + Ly*Jy + Lz*Jz) / magL / magJ;

  magLN = sqrt(inner_product(rcrossrdot, rcrossrdot));
  JLN =  (Jx*rcrossrdot[0] + Jy*rcrossrdot[1] + Jz*rcrossrdot[2])/magLN/magJ;

  
  if(debugPK) { 
      printf("chiJ = %3.10f\n", chiJ); fflush(NULL); fflush(NULL); 
      printf("J.L = %4.11f \n", kappaJL); fflush(NULL); 
      printf("J.LN = %4.11f \n", JLN);
      printf("L.LN = %4.11f \n", (Lx*rcrossrdot[0] + Ly*rcrossrdot[1] + Lz*rcrossrdot[2])/magL/magLN);
  }


  sh = 0.0;
  /* WaveStep 1.4: calculate combsize and deltaNQC */
  switch ( SpinAlignedEOBversion )
  {                     
    case 1:
      combSize = 7.5;
      deltaNQC = XLALSimIMREOBGetNRSpinPeakDeltaT(2, 2, eta,  chiJ);
      if ( debugPK ) {
        printf("v1 RD prescriptions are used! %3.10f %3.10f\n", 
                combSize, deltaNQC); 
        fflush(NULL);
      }
      break;
    case 2:
      combSize = 12.;       
      if ( chi1J == 0. && chi2J == 0. ) combSize = 11.;
      if (chiJ >= 0.8 && eta >30.0/(31.*31) && eta <10.0/121.) combSize = 13.5;
      if (chiJ >= 0.9 && eta < 30./(31.*31.)) combSize = 12.0;
      if (chiJ >= 0.8 && eta > 10./121.) combSize = 8.5;       
      deltaNQC = XLALSimIMREOBGetNRSpinPeakDeltaTv2(2, 2, m1, m2, chi1J, chi2J );
      // FIXME !!!!! 
      if ( debugPK ) {       
        printf("v2 RD prescriptions are used! %3.10f %3.10f\n", 
                combSize, deltaNQC);
        //if (deltaNQC > 2.5) {
        //    deltaNQC = 1.0;
        //    printf("WARNING!!!!!!! ----->> DELTANQC WAS SET TO 2.5 !!!!!\n");
        //}
        fflush(NULL);
      }
      break;
    default:
      XLALPrintError( "XLAL Error - %s: wrong SpinAlignedEOBversion value, must be 1 or 2!\n", __func__ );
      XLAL_ERROR( XLAL_EINVAL );
      break;
  }
  
  /*** Here is the old attempt of Yi to modify parameters according to the 
       opening angle  
  combSize    *= 1.0 + 9.0 * (1.0 - fabs(kappaJL));
  //longCombSize = combSize;
  deltaNQC    += 10.0 * (1.0 - fabs(kappaJL));*/
  
  // (Stas) !!! NOTE: tAttach is further modified by small shift "sh" computed and applied in XLALSimIMREOBHybridAttachRingdown !!!
  // FIXME
  tAttach = tPeakOmega - deltaNQC; //- 1.0;
  
  if (! found){
     tAttach = tPeakOmega;
  }

  if (debugPK){
    printf("For RD: DeltaNQC = %3.10f, comb = %3.10f \n", deltaNQC, combSize);
    printf("NOTE! that additional shift (sh) is computed and added in XLALSimIMREOBHybridAttachRingdown\n");
	  fflush(NULL);
  } 

  // FIXME
  //combSize = 40.0; 

  /* WaveStep 1.5: get  */
  /* WaveStep 1.6: construct J-frame */
  JframeEz[0] = Jx / magJ;
  JframeEz[1] = Jy / magJ;
  JframeEz[2] = Jz / magJ;
  if ( 1.-JframeEz[2] < 1.0e-13 )
  { JframeEx[0] = 1.; JframeEx[1] = 0.; JframeEx[2] = 0.; }
  else {
    JframeEx[0] = JframeEz[1];
    JframeEx[1] = -JframeEz[0];
    JframeEx[2] = 0.;
  }
  // FIXME try to change the sign of x:
  //if (JframeEx[0] < 0.0){
  //  JframeEx[0] = -JframeEz[1];
  //  JframeEx[1] = JframeEz[0];
  //}
  
  JExnorm = sqrt(JframeEx[0]*JframeEx[0] + JframeEx[1]*JframeEx[1]);
  //JframeEx[0] /= sqrt( JframeEz[0]*JframeEz[0] + JframeEz[1]*JframeEz[1] );
  //JframeEx[1] /= sqrt( JframeEz[0]*JframeEz[0] + JframeEz[1]*JframeEz[1] );
    
  JframeEx[0] /= JExnorm;
  JframeEx[1] /= JExnorm;
  JframeEy[0] = JframeEz[1]*JframeEx[2] - JframeEz[2]*JframeEx[1];
  JframeEy[1] = JframeEz[2]*JframeEx[0] - JframeEz[0]*JframeEx[2];
  JframeEy[2] = JframeEz[0]*JframeEx[1] - JframeEz[1]*JframeEx[0];
  

  // FIXME try to change the sign of x:
  //JframeEy[0] = - JframeEy[0]; 
  //JframeEy[1] = - JframeEy[1]; 
  //JframeEy[2] = - JframeEy[2]; 

  if(debugPK) { 
    printf("J-frameEx = [%e\t%e\t%e]\n", JframeEx[0], JframeEx[1], JframeEx[2]);printf("J-frameEy = [%e\t%e\t%e]\n", JframeEy[0], JframeEy[1], JframeEy[2]);
    printf("J-frameEz = [%e\t%e\t%e]\n", JframeEz[0], JframeEz[1], JframeEz[2]);fflush(NULL);
  }

  /* WaveStep 2
   * Calculate quasi-nonprecessing waveforms
   */
  /* WaveStep 2.1: create time-series containers for euler angles and hlm harmonics */
  retLen = retLenLow;
  
  REAL8TimeSeries *alphaI2PTS = XLALCreateREAL8TimeSeries( "alphaI2P", 
                                    &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  REAL8TimeSeries  *betaI2PTS = XLALCreateREAL8TimeSeries(  "betaI2P", 
                                    &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  REAL8TimeSeries *gammaI2PTS = XLALCreateREAL8TimeSeries( "gammaI2P", 
                                    &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  REAL8TimeSeries *alphaP2JTS = XLALCreateREAL8TimeSeries( "alphaP2J", 
                                    &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  REAL8TimeSeries  *betaP2JTS = XLALCreateREAL8TimeSeries(  "betaP2J", 
                                    &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  REAL8TimeSeries *gammaP2JTS = XLALCreateREAL8TimeSeries( "gammaP2J", 
                                    &tc, 0.0, deltaT, &lalStrainUnit, retLen );

  COMPLEX16TimeSeries *h22TS   = XLALCreateCOMPLEX16TimeSeries( "H_22",  
                                  &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h21TS   = XLALCreateCOMPLEX16TimeSeries( "H_21",   
                                  &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h20TS   = XLALCreateCOMPLEX16TimeSeries( "H_20",   
                                  &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h2m1TS  = XLALCreateCOMPLEX16TimeSeries( "H_2m1", 
                                  &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h2m2TS  = XLALCreateCOMPLEX16TimeSeries( "H_2m2",  
                                  &tc, 0.0, deltaT, &lalStrainUnit, retLen );


  COMPLEX16TimeSeries *h22PTS  = NULL;  
  COMPLEX16TimeSeries *h21PTS  = NULL;
  COMPLEX16TimeSeries *h20PTS  = NULL;
  COMPLEX16TimeSeries *h2m1PTS = NULL;
  COMPLEX16TimeSeries *h2m2PTS = NULL;
  
  COMPLEX16TimeSeries *h22JTS  = NULL; 
  COMPLEX16TimeSeries *h21JTS  = NULL;
  COMPLEX16TimeSeries *h20JTS  = NULL;
  COMPLEX16TimeSeries *h2m1JTS = NULL;
  COMPLEX16TimeSeries *h2m2JTS = NULL;
  
  COMPLEX16TimeSeries *hJTS    = NULL;
  
  COMPLEX16TimeSeries *hIMR22JTS  = NULL;
  COMPLEX16TimeSeries *hIMR21JTS  = NULL;
  COMPLEX16TimeSeries *hIMR20JTS  = NULL;
  COMPLEX16TimeSeries *hIMR2m1JTS = NULL;
  COMPLEX16TimeSeries *hIMR2m2JTS = NULL;
  
  //COMPLEX16TimeSeries *h22PTS  = XLALCreateCOMPLEX16TimeSeries( "HP_22",  
  //                                &tc, 0.0, deltaT, &lalStrainUnit, retLen );

  //COMPLEX16TimeSeries *h22JTS  = XLALCreateCOMPLEX16TimeSeries( "HJ_22",  
  //                                &tc, 0.0, deltaT, &lalStrainUnit, retLen );
  
  //COMPLEX16TimeSeries *hJTS    = XLALCreateCOMPLEX16TimeSeries( "HJ",     
  //                                &tc, 0.0, deltaT, &lalStrainUnit, retLen );

  //COMPLEX16TimeSeries *hIMR22JTS  = XLALCreateCOMPLEX16TimeSeries( "HIMRJ_22",
  //            &tc, 0.0, deltaT, &lalStrainUnit, retLen + retLenRDPatchLow );
  
  COMPLEX16TimeSeries *hIMRJTS    = XLALCreateCOMPLEX16TimeSeries( "HIMRJ", 
              &tc, 0.0, deltaT, &lalStrainUnit, retLen + retLenRDPatchLow );

  
  REAL8TimeSeries  *hPlusTS  = XLALCreateREAL8TimeSeries( "H_PLUS",
              &tc, 0.0, deltaT, &lalStrainUnit, retLen + retLenRDPatchLow );
  REAL8TimeSeries  *hCrossTS = XLALCreateREAL8TimeSeries( "H_CROSS", 
              &tc, 0.0, deltaT, &lalStrainUnit, retLen + retLenRDPatchLow );
  
  if ( !(tlist = XLALCreateREAL8Vector( retLen )) 
    || !(tlistRDPatch = XLALCreateREAL8Vector( retLen + retLenRDPatchLow )) )
  {
    XLAL_ERROR(  XLAL_ENOMEM );
  }
  memset( tlist->data, 0, tlist->length * sizeof( REAL8 ));
  memset( tlistRDPatch->data, 0, tlistRDPatch->length * sizeof( REAL8 ));
  for ( i = 0; i < retLen; i++ )
  {
    tlist->data[i] = i * deltaT/mTScaled;
    tlistRDPatch->data[i] = i * deltaT/mTScaled;
  }
  for ( i = retLen; i < retLen + retLenRDPatchLow; i++ )
  {
    tlistRDPatch->data[i] = i * deltaT/mTScaled;
  }

  /* WaveStep 2.3.1: main loop for quasi-nonprecessing waveform generation */
  // Generating modes for coarsely sampled portion 
  if(debugPK){
    printf("Generating precessing-frame modes for coarse dynamics\n");
    fflush(NULL); 
    out = fopen( "rotDynamics.dat", "w" );
  }
  
  for ( i = 0; i < retLen; i++ )
  {
    for ( j = 0; j < values->length; j++ )
    {
      values->data[j] = dynamics->data[(j+1)*retLen + i];
    }
    
    /* Calculate dr/dt */
    memset( dvalues->data, 0, 14*sizeof(dvalues->data[0]));
    if( XLALSpinHcapRvecDerivative( 0, values->data, dvalues->data, 
          &seobParams) != XLAL_SUCCESS )
    {
      printf( " Calculation of dr/dt at t = tPeakOmega \n");
      XLAL_ERROR( XLAL_EFUNC );
    }
    
    /* Calculate omega */
    memcpy( rdotvec, dvalues->data, 3*sizeof(REAL8));
    rvec[0] = posVecx.data[i]; rvec[1] = posVecy.data[i];
    rvec[2] = posVecz.data[i];
    cross_product( rvec, rdotvec, rcrossrdot );
    magR = sqrt(inner_product(rvec, rvec));
    omega = sqrt(inner_product(rcrossrdot, rcrossrdot)) / (magR*magR);
    vOmega = v = cbrt( omega );
    amp = amp0 * vOmega * vOmega;
    
    /* Cartesian vectors needed to calculate Hamiltonian */
    cartPosVec.length = cartMomVec.length = 3;
    cartPosVec.data = cartPosData;
    cartMomVec.data = cartMomData;
    memset( cartPosData, 0, sizeof( cartPosData ) );
    memset( cartMomData, 0, sizeof( cartMomData ) );

    LNhx  = rcrossrdot[0];
    LNhy  = rcrossrdot[1];
    LNhz  = rcrossrdot[2];
    magLN = sqrt(LNhx*LNhx + LNhy*LNhy + LNhz*LNhz);
    LNhx  = LNhx / magLN;
    LNhy  = LNhy / magLN;
    LNhz  = LNhz / magLN;

    /*aI2P = atan2( LNhy, LNhx );
    bI2P = acos( LNhz );
    gI2P = -phiMod.data[i]; */ /* this one is defined w.r.t. L not LN*/
    aI2P = Alpha->data[i];
    bI2P = Beta->data[i];
    gI2P = Gamma->data[i];

    LframeEx[0] =  cos(aI2P)*cos(bI2P)*cos(gI2P) - sin(aI2P)*sin(gI2P);
    LframeEx[1] =  sin(aI2P)*cos(bI2P)*cos(gI2P) + cos(aI2P)*sin(gI2P);
    LframeEx[2] = -sin(bI2P)*cos(gI2P);
    LframeEy[0] = -cos(aI2P)*cos(bI2P)*sin(gI2P) - sin(aI2P)*cos(gI2P);
    LframeEy[1] = -sin(aI2P)*cos(bI2P)*sin(gI2P) + cos(aI2P)*cos(gI2P);
    LframeEy[2] =  sin(bI2P)*sin(gI2P);
    LframeEz[0] =  LNhx;
    LframeEz[1] =  LNhy;
    LframeEz[2] =  LNhz;
    aP2J = atan2(JframeEz[0]*LframeEy[0]+JframeEz[1]*LframeEy[1]+JframeEz[2]*LframeEy[2],
                 JframeEz[0]*LframeEx[0]+JframeEz[1]*LframeEx[1]+JframeEz[2]*LframeEx[2]);
    bP2J = acos( JframeEz[0]*LframeEz[0]+JframeEz[1]*LframeEz[1]+JframeEz[2]*LframeEz[2]);
    gP2J = atan2(  JframeEy[0]*LframeEz[0]+JframeEy[1]*LframeEz[1]+JframeEy[2]*LframeEz[2],
                 -(JframeEx[0]*LframeEz[0]+JframeEx[1]*LframeEz[1]+JframeEx[2]*LframeEz[2]));

/*if (i==0||i==1900) printf("{{%f,%f,%f},{%f,%f,%f},{%f,%f,%f}}\n",JframeEx[0],JframeEx[1],JframeEx[2],JframeEy[0],JframeEy[1],JframeEy[2],JframeEz[0],JframeEz[1],JframeEz[2]);
if (i==0||i==1900) printf("{{%f,%f,%f},{%f,%f,%f},{%f,%f,%f}}\n",LframeEx[0],LframeEx[1],LframeEx[2],LframeEy[0],LframeEy[1],LframeEy[2],LframeEz[0],LframeEz[1],LframeEz[2]);
if (i==0||i==1900) printf("YP: study time = %f\n",i*deltaT/mTScaled);
if (i==1900) printf("YP: gamma: %f, %f, %f, %f\n", JframeEy[0]*LframeEz[0]+JframeEy[1]*LframeEz[1]+JframeEy[2]*LframeEz[2], JframeEx[0]*LframeEz[0]+JframeEx[1]*LframeEz[1]+JframeEx[2]*LframeEz[2], gP2J, atan2(-0.365446,-0.378524));
    */
    /* I2P Euler angles are stored only for debugging purposes */
    alphaI2PTS->data->data[i] = -aI2P;
     betaI2PTS->data->data[i] = bI2P;
    gammaI2PTS->data->data[i] = -gI2P;
    alphaP2JTS->data->data[i] = -aP2J;
    /*betaP2JTS->data->data[i] = 0.*LAL_PI/2.-bP2J;*/
    betaP2JTS->data->data[i] = bP2J;
    gammaP2JTS->data->data[i] = -gP2J;

    /* Calculate the value of the Hamiltonian */
    memcpy( cartPosVec.data, values->data,   3*sizeof(REAL8) );
    memcpy( cartMomVec.data, values->data+3, 3*sizeof(REAL8) );

    /*if (i == 287)
    {
      printf("%f, %f %f %f, %f %f %f\n",eta, cartPosVec.data[0], cartPosVec.data[1], cartPosVec.data[2], cartMomVec.data[0], cartMomVec.data[1], cartMomVec.data[2]);
    }*/
    
    /* Update Hamiltonian coefficients as per |Skerr| */
    for( k = 0; k < 3; k++ )
    {
			s1DataNorm[k] = values->data[k+6];
			s2DataNorm[k] = values->data[k+9];
			s1Data[k] = s1DataNorm[k] * mTotal * mTotal;
			s2Data[k] = s2DataNorm[k] * mTotal * mTotal;
    }
    s1Vec.data = s1Data;
    s2Vec.data = s2Data;
    
    s1VecOverMtMt.data = s1DataNorm;
    s2VecOverMtMt.data = s2DataNorm;
    
    seobParams.s1Vec = &s1VecOverMtMt;
    seobParams.s2Vec = &s2VecOverMtMt;
    
    XLALSimIMRSpinEOBCalculateSigmaStar( sigmaStar, m1, m2, &s1Vec, &s2Vec );
    XLALSimIMRSpinEOBCalculateSigmaKerr( sigmaKerr, m1, m2, &s1Vec, &s2Vec );
    
    seobParams.a = a = sqrt(inner_product(sigmaKerr->data, sigmaKerr->data));
    
    rcrossrdot[0] = LNhx; 
    rcrossrdot[1] = LNhy; 
    rcrossrdot[2] = LNhz;
    s1dotLN = inner_product(s1Data, rcrossrdot) / (m1*m1);
    s2dotLN = inner_product(s2Data, rcrossrdot) / (m2*m2);
    
    chiS = 0.5 * (s1dotLN + s2dotLN);
    chiA = 0.5 * (s1dotLN - s2dotLN);

    switch ( SpinAlignedEOBversion )
    {
     case 1:
       tplspin = 0.0;
       break;
     case 2:
       tplspin = (1.-2.*eta) * chiS + (m1 - m2)/(m1 + m2) * chiA;
       break;
     default:
       XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
       XLAL_ERROR( XLAL_EINVAL );
       break;
    }

    if ( XLALSimIMRCalculateSpinEOBHCoeffs( &seobCoeffs, eta, a, 
                          SpinAlignedEOBversion ) == XLAL_FAILURE )
      printf("\nSomething went wrong evaluating XLALSimIMRCalculateSpinEOBHCoeffs in step %d of coarse dynamics\n", 
			i );

    /* Update hlm coefficients */
    if ( XLALSimIMREOBCalcSpinFacWaveformCoefficients( &hCoeffs, m1, m2, eta, 
        tplspin, chiS, chiA, SpinAlignedEOBversion ) == XLAL_FAILURE )
      printf("\nSomething went wrong evaluating XLALSimIMRCalculateSpinEOBHCoeffs in step %d of coarse dynamics\n", 
			i );
    
    ham = XLALSimIMRSpinEOBHamiltonian( eta, &cartPosVec, &cartMomVec,
                  &s1VecOverMtMt, &s2VecOverMtMt,
                  sigmaKerr, sigmaStar, seobParams.tortoise, &seobCoeffs );

    /* Calculate the polar data */
    polarDynamics.length = 4;
    polarDynamics.data   = polData;
    
    /* Calculate the orbital angular momentum */
    cross_product( values->data, values->data+3, rcrossp );
    magL = sqrt(inner_product(rcrossp, rcrossp));

    polData[0] = sqrt(inner_product(values->data, values->data));
    polData[1] = phiMod.data[i] + phiDMod.data[i];
    polData[2] = inner_product(values->data, values->data+3) / polData[0];
    polData[3] = magL;

    if ( XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform( &hLM, 
          &polarDynamics, values, v, ham, 2, 2, &seobParams ) == XLAL_FAILURE )
    {
      XLAL_ERROR( XLAL_EFUNC );
    }
    if ( XLALSimIMRSpinEOBNonQCCorrection( &hNQC, 
          values, omega, &nqcCoeffs ) == XLAL_FAILURE )
    {
      XLAL_ERROR( XLAL_EFUNC );
    }
    hLM *= hNQC;
    h22TS->data->data[i]  = hLM;
    h2m2TS->data->data[i] = conjl(hLM);

    if (debugPK) {
      fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e ",
             i*deltaT/mTScaled, aI2P, bI2P, gI2P, aP2J, bP2J, gP2J );
      //fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
      //     vX, vY, vZ, LNhx, LNhy, LNhz, creal(hLM), cimag(hLM) );
      fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
                      polData[0], polData[1], polData[2], polData[3], ham, v, 
                      creal(hLM), cimag(hLM), creal(hLM/hNQC), cimag(hLM/hNQC) );
      fflush(NULL);
    }

    if ( XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform( &hLM, 
        &polarDynamics, values, v, ham, 2, 1, &seobParams ) == XLAL_FAILURE )
    {
      XLAL_ERROR( XLAL_EFUNC );
    }
    h21TS->data->data[i]  = hLM;
    h2m1TS->data->data[i] = conjl(hLM);
    h20TS->data->data[i]  = 0.0;

    /*if (i == 95 )
    {
      printf("%.16e %.16e %.16e\n",ham, omega, v);
      printf("%.16e %.16e %.16e\n",values->data[0],values->data[1],values->data[2]);
      printf("%.16e %.16e %.16e\n",values->data[3],values->data[4],values->data[5]);
      printf("%.16e %.16e %.16e\n",values->data[6],values->data[7],values->data[8]);
      printf("%.16e %.16e %.16e\n",values->data[9],values->data[10],values->data[11]);
      printf("%.16e %.16e %.16e %.16e %.16e %.16e\n",nqcCoeffs.a1,nqcCoeffs.a2,nqcCoeffs.a3,nqcCoeffs.a3S,nqcCoeffs.a4,nqcCoeffs.a5);
      printf("%.16e %.16e %.16e %.16e\n",nqcCoeffs.b1,nqcCoeffs.b2,nqcCoeffs.b3,nqcCoeffs.b4);
      printf("%.16e %.16e, %.16e %.16e\n",creal(hLM),cimag(hLM),creal(hNQC),cimag(hNQC));
    }*/

    /*hPlusTS->data->data[i]  = - 0.5 * amp * cos( 2.*vphi[i]) * cos(2.*alpha) * (1. + LNhz*LNhz)
                            + amp * sin(2.*vphi[i]) * sin(2.*alpha)*LNhz;

    hCrossTS->data->data[i] = - 0.5 * amp * cos( 2.*vphi[i]) * sin(2.*alpha) * (1. + LNhz*LNhz)
                            - amp * sin(2.*vphi[i]) * cos(2.*alpha) * LNhz;*/

  }
  if(debugPK) {
      fclose( out );
      printf("YP: quasi-nonprecessing modes generated.\n"); fflush(NULL);
  }

  /* WaveStep 2.4.1: add quasi-nonprecessing spherical harmonic modes to the SphHarmTimeSeries structure */
  hlmPTS = XLALSphHarmTimeSeriesAddMode( hlmPTS, h22TS, 2, 2 );
  hlmPTS = XLALSphHarmTimeSeriesAddMode( hlmPTS, h21TS, 2, 1 );
  hlmPTS = XLALSphHarmTimeSeriesAddMode( hlmPTS, h20TS, 2, 0 );
  hlmPTS = XLALSphHarmTimeSeriesAddMode( hlmPTS, h2m1TS, 2, -1 );
  hlmPTS = XLALSphHarmTimeSeriesAddMode( hlmPTS, h2m2TS, 2, -2 );
  XLALSphHarmTimeSeriesSetTData( hlmPTS, tlist );

  h22PTS  = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, 2 );
  h21PTS  = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, 1 );
  h20PTS  = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, 0 );
  h2m1PTS = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, -1);
  h2m2PTS = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, -2);
  
  /* Write waveforms in precessing frame */
  if (debugPK) {
    printf("YP: SphHarmTS structures populated.\n"); fflush(NULL);
    out = fopen( "PWaves.dat", "w" );
    for ( i = 0; i < retLen; i++ )
    {
      fprintf( out, 
      "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
          i*deltaT/mTScaled, 
          creal(h22PTS->data->data[i]), cimag(h22PTS->data->data[i]),
          creal(h21PTS->data->data[i]), cimag(h21PTS->data->data[i]),
          creal(h20PTS->data->data[i]), cimag(h20PTS->data->data[i]),
          creal(h2m1PTS->data->data[i]), cimag(h2m1PTS->data->data[i]),
          creal(h2m2PTS->data->data[i]), cimag(h2m2PTS->data->data[i]) );
    }
    fclose( out );
    printf("YP: P-frame waveforms written to file.\n");
    fflush(NULL);
  }


  /* WaveStep 2.3.2: main loop for quasi-nonprecessing waveform generation -- HIGH SR */
  retLen = retLenHi;
  
  REAL8TimeSeries *alphaI2PTSHi = XLALCreateREAL8TimeSeries( "alphaI2P", 
                                &tc, 0.0, deltaTHigh, &lalStrainUnit, retLen );
  REAL8TimeSeries  *betaI2PTSHi = XLALCreateREAL8TimeSeries(  "betaI2P", 
                                &tc, 0.0, deltaTHigh, &lalStrainUnit, retLen );
  REAL8TimeSeries *gammaI2PTSHi = XLALCreateREAL8TimeSeries( "gammaI2P", 
                                &tc, 0.0, deltaTHigh, &lalStrainUnit, retLen );
  REAL8TimeSeries *alphaP2JTSHi = XLALCreateREAL8TimeSeries( "alphaP2J", 
                                &tc, 0.0, deltaTHigh, &lalStrainUnit, retLen );
  REAL8TimeSeries  *betaP2JTSHi = XLALCreateREAL8TimeSeries(  "betaP2J", 
                                &tc, 0.0, deltaTHigh, &lalStrainUnit, retLen );
  REAL8TimeSeries *gammaP2JTSHi = XLALCreateREAL8TimeSeries( "gammaP2J", 
                                &tc, 0.0, deltaTHigh, &lalStrainUnit, retLen );

  COMPLEX16TimeSeries *h22TSHi   = XLALCreateCOMPLEX16TimeSeries( "H_22",   
                                &tc, 0.0, deltaTHigh, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h21TSHi   = XLALCreateCOMPLEX16TimeSeries( "H_21",   
                                &tc, 0.0, deltaTHigh, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h20TSHi   = XLALCreateCOMPLEX16TimeSeries( "H_20",   
                                &tc, 0.0, deltaTHigh, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h2m1TSHi  = XLALCreateCOMPLEX16TimeSeries( "H_2m1",  
                                &tc, 0.0, deltaTHigh, &lalStrainUnit, retLen );
  COMPLEX16TimeSeries *h2m2TSHi  = XLALCreateCOMPLEX16TimeSeries( "H_2m2",  
                                &tc, 0.0, deltaTHigh, &lalStrainUnit, retLen );
  
  
  COMPLEX16TimeSeries *h22PTSHi  = NULL; 
  COMPLEX16TimeSeries *h21PTSHi  = NULL;
  COMPLEX16TimeSeries *h20PTSHi  = NULL;
  COMPLEX16TimeSeries *h2m1PTSHi = NULL;
  COMPLEX16TimeSeries *h2m2PTSHi = NULL;
  
  COMPLEX16TimeSeries *h22JTSHi  = NULL;
  COMPLEX16TimeSeries *h21JTSHi  = NULL;
  COMPLEX16TimeSeries *h20JTSHi  = NULL;
  COMPLEX16TimeSeries *h2m1JTSHi = NULL;
  COMPLEX16TimeSeries *h2m2JTSHi = NULL;
  COMPLEX16TimeSeries *hJTSHi    = NULL;
  
  COMPLEX16TimeSeries *hIMR22JTSHi  = NULL;
  COMPLEX16TimeSeries *hIMR21JTSHi  = NULL;
  COMPLEX16TimeSeries *hIMR20JTSHi  = NULL;
  COMPLEX16TimeSeries *hIMR2m1JTSHi = NULL;
  COMPLEX16TimeSeries *hIMR2m2JTSHi = NULL;
  
  //COMPLEX16TimeSeries *h22PTSHi  = XLALCreateCOMPLEX16TimeSeries( "HP_22",  
  //                              &tc, 0.0, deltaTHigh, &lalStrainUnit, retLen );
  
  //COMPLEX16TimeSeries *h22JTSHi  = XLALCreateCOMPLEX16TimeSeries( "HJ_22",  
  //                              &tc, 0.0, deltaTHigh, &lalStrainUnit, retLen );

  //COMPLEX16TimeSeries *hIMR22JTSHi  = XLALCreateCOMPLEX16TimeSeries(
  // "HIMRJ_22Hi",  &tc, 0.0, deltaTHigh, &lalStrainUnit, retLen + retLenRDPatch);

  COMPLEX16TimeSeries *hIMRJTSHi    = XLALCreateCOMPLEX16TimeSeries(
   "HIMRJHi",     &tc, 0.0, deltaTHigh, &lalStrainUnit, retLen + retLenRDPatch);


   
  if ( !(tlistHi = XLALCreateREAL8Vector( retLen )) || !(tlistRDPatchHi = XLALCreateREAL8Vector( retLen + retLenRDPatch )) )
  {
    XLAL_ERROR(  XLAL_ENOMEM );
  }
  memset( tlistHi->data, 0, tlistHi->length * sizeof( REAL8 ));
  memset( tlistRDPatchHi->data, 0, tlistRDPatchHi->length * sizeof( REAL8 ));
  for ( i = 0; i < retLen; i++ )
  {
    tlistHi->data[i] = i * deltaTHigh/mTScaled + HiSRstart;
    tlistRDPatchHi->data[i] = i * deltaTHigh/mTScaled + HiSRstart;
  }
  for ( i = retLen; i < retLen + retLenRDPatch; i++ )
  {
    tlistRDPatchHi->data[i] = i * deltaTHigh/mTScaled + HiSRstart;
  }

  
  // Generating modes for finely sampled portion 
  if(debugPK) {
    printf("Generating precessing-frame modes for fine dynamics\n");
    fflush(NULL);
  }
  if (debugPK) out = fopen( "rotDynamicsHi.dat", "w" );
  for ( i = 0; i < retLen; i++ )
  {
    for ( j = 0; j < values->length; j++ )
    {
      values->data[j] = dynamicsHi->data[(j+1)*retLen + i];
      /*if ( i == 0 ) printf("value (length=%d) %d: %.16e\n", values->length, j,values->data[j]);*/
    }

    /* Calculate dr/dt */
    memset( dvalues->data, 0, 14*sizeof(dvalues->data[0]));
    if( XLALSpinHcapRvecDerivative( 0, values->data, dvalues->data, 
          &seobParams) != XLAL_SUCCESS )
    {
      printf( " Calculation of dr/dt at t = tPeakOmega \n");
      XLAL_ERROR( XLAL_EFUNC );
    }
    
    /* Calculate omega */
    rvec[0] = posVecxHi.data[i]; rvec[1] = posVecyHi.data[i];
    rvec[2] = posVeczHi.data[i];
    memcpy( rdotvec, dvalues->data, 3*sizeof(REAL8));
    vX = rdotvec[0]; vY = rdotvec[1]; vZ = rdotvec[2]; 
    cross_product( rvec, rdotvec, rcrossrdot );
    
    magR   = sqrt(inner_product(rvec, rvec));
    omega  = sqrt(inner_product(rcrossrdot, rcrossrdot)) / (magR*magR);
    vOmega = v = cbrt( omega );
    amp = amp0 * vOmega * vOmega;

    /* Cartesian vectors needed to calculate Hamiltonian */
    cartPosVec.length = cartMomVec.length = 3;
    cartPosVec.data = cartPosData;
    cartMomVec.data = cartMomData;
    memset( cartPosData, 0, sizeof( cartPosData ) );
    memset( cartMomData, 0, sizeof( cartMomData ) );

    magLN = sqrt(inner_product(rcrossrdot, rcrossrdot));
    LNhx  = rcrossrdot[0] / magLN;
    LNhy  = rcrossrdot[1] / magLN;
    LNhz  = rcrossrdot[2] / magLN;

    /*aI2P = atan2( LNhy, LNhx );
    bI2P = acos( LNhz );
    gI2P = -phiModHi.data[i]; */  
    aI2P = AlphaHi->data[i];
    bI2P = BetaHi->data[i];
    gI2P = GammaHi->data[i];
    LframeEx[0] =  cos(aI2P)*cos(bI2P)*cos(gI2P) - sin(aI2P)*sin(gI2P);
    LframeEx[1] =  sin(aI2P)*cos(bI2P)*cos(gI2P) + cos(aI2P)*sin(gI2P);
    LframeEx[2] = -sin(bI2P)*cos(gI2P);
    LframeEy[0] = -cos(aI2P)*cos(bI2P)*sin(gI2P) - sin(aI2P)*cos(gI2P);
    LframeEy[1] = -sin(aI2P)*cos(bI2P)*sin(gI2P) + cos(aI2P)*cos(gI2P);
    LframeEy[2] =  sin(bI2P)*sin(gI2P);
    LframeEz[0] =  LNhx;
    LframeEz[1] =  LNhy;
    LframeEz[2] =  LNhz;
    //if(debugPK) { 
    //   printf("L-frameEx = [%e\t%e\t%e]\n", LframeEx[0], LframeEx[1], LframeEx[2]);printf("L-frameEy = [%e\t%e\t%e]\n", LframeEy[0], LframeEy[1], LframeEy[2]);
    //   printf("L-frameEz = [%e\t%e\t%e]\n", LframeEz[0], LframeEz[1], LframeEz[2]);fflush(NULL);
    //}

    aP2J = atan2(JframeEz[0]*LframeEy[0]+JframeEz[1]*LframeEy[1]+JframeEz[2]*LframeEy[2],
                 JframeEz[0]*LframeEx[0]+JframeEz[1]*LframeEx[1]+JframeEz[2]*LframeEx[2]);
    bP2J = acos( JframeEz[0]*LframeEz[0]+JframeEz[1]*LframeEz[1]+JframeEz[2]*LframeEz[2]);
    gP2J = atan2(  JframeEy[0]*LframeEz[0]+JframeEy[1]*LframeEz[1]+JframeEy[2]*LframeEz[2],
                 -(JframeEx[0]*LframeEz[0]+JframeEx[1]*LframeEz[1]+JframeEx[2]*LframeEz[2]));

/*if (i==0||i==1900) printf("{{%f,%f,%f},{%f,%f,%f},{%f,%f,%f}}\n",JframeEx[0],JframeEx[1],JframeEx[2],JframeEy[0],JframeEy[1],JframeEy[2],JframeEz[0],JframeEz[1],JframeEz[2]);
if (i==0||i==1900) printf("{{%f,%f,%f},{%f,%f,%f},{%f,%f,%f}}\n",LframeEx[0],LframeEx[1],LframeEx[2],LframeEy[0],LframeEy[1],LframeEy[2],LframeEz[0],LframeEz[1],LframeEz[2]);
if (i==0||i==1900) printf("YP: study time = %f\n",i*deltaTHigh/mTScaled);
if (i==1900) printf("YP: gamma: %f, %f, %f, %f\n", JframeEy[0]*LframeEz[0]+JframeEy[1]*LframeEz[1]+JframeEy[2]*LframeEz[2], JframeEx[0]*LframeEz[0]+JframeEx[1]*LframeEz[1]+JframeEx[2]*LframeEz[2], gP2J, atan2(-0.365446,-0.378524));*/
    /* I2P Euler angles are stored only for debugging purposes */
    alphaI2PTSHi->data->data[i] = -aI2P;
     betaI2PTSHi->data->data[i] = bI2P;
    gammaI2PTSHi->data->data[i] = -gI2P;
    alphaP2JTSHi->data->data[i] = -aP2J;
     betaP2JTSHi->data->data[i] = bP2J;
    gammaP2JTSHi->data->data[i] = -gP2J;

    /* Calculate the value of the Hamiltonian */
    memcpy( cartPosVec.data, values->data,   3*sizeof(REAL8) );
    memcpy( cartMomVec.data, values->data+3, 3*sizeof(REAL8) );

    /* Update Hamiltonian coefficients as per |Skerr| */
    for( k = 0; k < 3; k++ )
    {
			s1DataNorm[k] = values->data[k+6];
			s2DataNorm[k] = values->data[k+9];
			s1Data[k] = s1DataNorm[k] * mTotal * mTotal;
			s2Data[k] = s2DataNorm[k] * mTotal * mTotal;
    }
    s1Vec.data = s1Data;
    s2Vec.data = s2Data;
    
    s1VecOverMtMt.data = s1DataNorm;
    s2VecOverMtMt.data = s2DataNorm;
    
    seobParams.s1Vec = &s1VecOverMtMt;
    seobParams.s2Vec = &s2VecOverMtMt;
    
    XLALSimIMRSpinEOBCalculateSigmaStar( sigmaStar, m1, m2, &s1Vec, &s2Vec );
    XLALSimIMRSpinEOBCalculateSigmaKerr( sigmaKerr, m1, m2, &s1Vec, &s2Vec );
    
    seobParams.a = a = sqrt(inner_product(sigmaKerr->data, sigmaKerr->data));
    
    rcrossrdot[0] = LNhx; 
    rcrossrdot[1] = LNhy; 
    rcrossrdot[2] = LNhz;
    s1dotLN = inner_product(s1Data, rcrossrdot) / (m1*m1);
    s2dotLN = inner_product(s2Data, rcrossrdot) / (m2*m2);
    
    chiS = 0.5 * (s1dotLN + s2dotLN);
    chiA = 0.5 * (s1dotLN - s2dotLN);

    switch ( SpinAlignedEOBversion )
    {
     case 1:
       tplspin = 0.0;
       break;
     case 2:
       tplspin = (1.-2.*eta) * chiS + (m1 - m2)/(m1 + m2) * chiA;
       break;
     default:
       XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
       XLAL_ERROR( XLAL_EINVAL );
       break;
    }

    if ( XLALSimIMRCalculateSpinEOBHCoeffs( &seobCoeffs, eta, a, 
                          SpinAlignedEOBversion ) == XLAL_FAILURE )
      printf("\nSomething went wrong evaluating XLALSimIMRCalculateSpinEOBHCoeffs in step %d of coarse dynamics\n", 
			i );

    /* Update hlm coefficients */
    if ( XLALSimIMREOBCalcSpinFacWaveformCoefficients( &hCoeffs, m1, m2, eta, 
        tplspin, chiS, chiA, SpinAlignedEOBversion ) == XLAL_FAILURE )
      printf("\nSomething went wrong evaluating XLALSimIMRCalculateSpinEOBHCoeffs in step %d of coarse dynamics\n", 
			i );
    
    ham = XLALSimIMRSpinEOBHamiltonian( eta, &cartPosVec, &cartMomVec,
                  &s1VecOverMtMt, &s2VecOverMtMt,
                  sigmaKerr, sigmaStar, seobParams.tortoise, &seobCoeffs );

    /* Calculate the polar data */
    polarDynamics.length = 4;
    polarDynamics.data   = polData;
    
    /* Calculate the orbital angular momentum */
    cross_product(values->data, values->data+3, rcrossp);

    polData[0] = sqrt(inner_product(values->data, values->data));
    polData[1] = phiModHi.data[i] + phiDModHi.data[i];
    polData[2] = inner_product(values->data, values->data+3) / polData[0];
    polData[3] = sqrt(inner_product(rcrossp, rcrossp));


    if ( XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform( &hLM, &polarDynamics, values, v, ham, 2, 2, &seobParams )
           == XLAL_FAILURE )
    {
      XLAL_ERROR( XLAL_EFUNC );
    }
    if ( XLALSimIMRSpinEOBNonQCCorrection( &hNQC, values, omega, &nqcCoeffs ) == XLAL_FAILURE )
    {
      XLAL_ERROR( XLAL_EFUNC );
    }
    hLM *= hNQC;
    h22TSHi->data->data[i]  = hLM;
    h2m2TSHi->data->data[i] = conjl(hLM);

    if (debugPK){
       fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e ",
             i*deltaTHigh/mTScaled + HiSRstart, aI2P, bI2P, gI2P, aP2J, bP2J, gP2J );
       fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e, %.16e\n",
             vX, vY, vZ, LNhx, LNhy, LNhz, creal(hLM), cimag(hLM), polData[0]);
    }

    if ( XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform( &hLM, &polarDynamics, values, v, ham, 2, 1, &seobParams )
           == XLAL_FAILURE )
    {
      XLAL_ERROR( XLAL_EFUNC );
    }
    h21TSHi->data->data[i]  = hLM;
    h2m1TSHi->data->data[i] = conjl(hLM);
    h20TSHi->data->data[i]  = 0.0;

    /*if (i == 95 )
    {
      printf("%.16e %.16e %.16e\n",ham, omega, v);
      printf("%.16e %.16e %.16e\n",values->data[0],values->data[1],values->data[2]);
      printf("%.16e %.16e %.16e\n",values->data[3],values->data[4],values->data[5]);
      printf("%.16e %.16e %.16e\n",values->data[6],values->data[7],values->data[8]);
      printf("%.16e %.16e %.16e\n",values->data[9],values->data[10],values->data[11]);
      printf("%.16e %.16e %.16e %.16e %.16e %.16e\n",nqcCoeffs.a1,nqcCoeffs.a2,nqcCoeffs.a3,nqcCoeffs.a3S,nqcCoeffs.a4,nqcCoeffs.a5);
      printf("%.16e %.16e %.16e %.16e\n",nqcCoeffs.b1,nqcCoeffs.b2,nqcCoeffs.b3,nqcCoeffs.b4);
      printf("%.16e %.16e, %.16e %.16e\n\n",creal(hLM),cimag(hLM),creal(hNQC),cimag(hNQC));
    }*/

    /*hPlusTS->data->data[i]  = - 0.5 * amp * cos( 2.*vphi[i]) * cos(2.*alpha) * (1. + LNhz*LNhz)
                            + amp * sin(2.*vphi[i]) * sin(2.*alpha)*LNhz;

    hCrossTS->data->data[i] = - 0.5 * amp * cos( 2.*vphi[i]) * sin(2.*alpha) * (1. + LNhz*LNhz)
                            - amp * sin(2.*vphi[i]) * cos(2.*alpha) * LNhz;*/

  }
  if (debugPK){
      fclose( out );
      printf("YP: quasi-nonprecessing modes generated.for High SR\n");
  }

  /* WaveStep 2.4.2: add quasi-nonprecessing spherical harmonic modes to the SphHarmTimeSeries structure */
  hlmPTSHi = XLALSphHarmTimeSeriesAddMode( hlmPTSHi, h22TSHi, 2, 2 );
  hlmPTSHi = XLALSphHarmTimeSeriesAddMode( hlmPTSHi, h21TSHi, 2, 1 );
  hlmPTSHi = XLALSphHarmTimeSeriesAddMode( hlmPTSHi, h20TSHi, 2, 0 );
  hlmPTSHi = XLALSphHarmTimeSeriesAddMode( hlmPTSHi, h2m1TSHi, 2, -1 );
  hlmPTSHi = XLALSphHarmTimeSeriesAddMode( hlmPTSHi, h2m2TSHi, 2, -2 );
  XLALSphHarmTimeSeriesSetTData( hlmPTSHi, tlistHi );
  
  hlmPTSout = XLALSphHarmTimeSeriesAddMode( hlmPTSout, h22TSHi, 2, 2 );
  hlmPTSout = XLALSphHarmTimeSeriesAddMode( hlmPTSout, h21TSHi, 2, 1 );
  hlmPTSout = XLALSphHarmTimeSeriesAddMode( hlmPTSout, h20TSHi, 2, 0 );
  hlmPTSout = XLALSphHarmTimeSeriesAddMode( hlmPTSout, h2m1TSHi, 2, -1 );
  hlmPTSout = XLALSphHarmTimeSeriesAddMode( hlmPTSout, h2m2TSHi, 2, -2 );
  //XLALSphHarmTimeSeriesSetTData( hlmPTSout, tlistHi );
  *hlmPTSoutput = hlmPTSout;

  h22PTSHi  = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, 2 );
  h21PTSHi  = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, 1 );
  h20PTSHi  = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, 0 );
  h2m1PTSHi = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, -1);
  h2m2PTSHi = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, -2);
  
  if (debugPK){ 
    printf("YP: SphHarmTS structures populated for HIgh SR.\n"); 

    out = fopen( "PWavesHi.dat", "w" );
    for ( i = 0; i < retLen; i++ )
    {
      fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
        i*deltaTHigh/mTScaled + HiSRstart, creal(h22PTSHi->data->data[i]),
        cimag(h22PTSHi->data->data[i]),
        creal(h21PTSHi->data->data[i]), cimag(h21PTSHi->data->data[i]),
        creal(h20PTSHi->data->data[i]), cimag(h20PTSHi->data->data[i]),
        creal(h2m1PTSHi->data->data[i]), cimag(h2m1PTSHi->data->data[i]),
        creal(h2m2PTSHi->data->data[i]), cimag(h2m2PTSHi->data->data[i]));
     }
     fclose( out );
     printf("YP: P-frame waveforms written to file for High SR.\n");
     fflush(NULL);
  }
  
//  // Find tAmpMax to determin the tAttachment
//    REAL8 *radiusVec;
//    REAL8 radiusData[retLen];
//    radiusVec = &radiusData[0];
//    for ( i = 0; i < retLen; i++ )
//    {
//        for ( j = 0; j < 3; j++ )
//        {
//            values->data[j] = dynamicsHi->data[(j+1)*retLen + i];
//        }
//        radiusData[i] = sqrt(values->data[0]*values->data[0] + values->data[1]*values->data[1] + values->data[2]*values->data[2]);
//    }
//    if (debugPK){
//        printf("Andrea the attachment time based on Omega is %f \n", tAttach);
//    }
  int foundAmp = 0;
    REAL8 tMaxAmp;
  double tAmpMax =  XLALSimLocateAmplTime(&timeHi, h22PTSHi->data, radiusVec, &foundAmp, &tMaxAmp);
  
  if(foundAmp==0){
      if (debugPK){
          printf("maximum of Amplitude or dot{Ampl} are not found \n");
      }
  }
  else{
      if (debugPK){
          printf("The Amplitude-related time is found and it's %f \n", tAmpMax);
      }
  }
  if( foundAmp==0 && found == 0){
      XLALPrintError("Houston, we've got a problem SOS, SOS, SOS, cannot find the RD attachment point...\n");
      XLAL_ERROR( XLAL_EINVAL );
      
  }
   
  
  if (tAmpMax < tAttach){
      if (debugPK){
          printf("Stas the attachment time will be modified from %f to %f \n", tAttach, tAmpMax);
      }      
      tAttach = tAmpMax;
  }
  // FIXME 
  //tAttach = 136.3;
//  tAttach  = tAttach - 6.0;


    //for (i=0; i<5; i++){
    //    printf("AttachPars stas: %d, %f \n", i, AttachParams->data[i]);
    //}
   

  if (debugPK){
      printf("Stas: the final decision on the attachment time is %f \n", tAttach);
  }
  //exit(0);
  
  /*if(debugPK) {
    printf( "Estimation of the peak is now at time %.16e, %.16e \n", 
          tPeakOmega, tPeakOmega+HiSRstart);
    fflush(NULL);
  }*/

  /* Changing retLen back to the Low SR value */
  retLen = retLenLow;

  /* WaveStep 3
   * Generate IMR waveforms in the merger J (~ final spin) -frame
   */
  if ( XLALSimInspiralPrecessionRotateModes( hlmPTS, alphaP2JTS, betaP2JTS, gammaP2JTS ) == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }
  h22JTS  = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, 2 );
  h21JTS  = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, 1 );
  h20JTS  = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, 0 );
  h2m1JTS = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, -1);
  h2m2JTS = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, -2);

  
  if (debugPK){
    printf("YP: PtoJ rotation done.\n");

    out = fopen( "JWaves.dat", "w" );
    for ( i = 0; i < retLen; i++ )
    {
      fprintf( out, 
      "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
        i*deltaT/mTScaled, 
        creal(h22JTS->data->data[i]), cimag(h22JTS->data->data[i]),
        creal(h21JTS->data->data[i]), cimag(h21JTS->data->data[i]),
        creal(h20JTS->data->data[i]), cimag(h20JTS->data->data[i]),
        creal(h2m1JTS->data->data[i]), cimag(h2m1JTS->data->data[i]),
        creal(h2m2JTS->data->data[i]), cimag(h2m2JTS->data->data[i]) );
    }
    fclose( out );
    printf("YP: P-frame waveforms written to file.\n");
    fflush(NULL);
  }

  /* Stas: Rotating the high sampling part */  
  if ( XLALSimInspiralPrecessionRotateModes( hlmPTSHi, 
        alphaP2JTSHi, betaP2JTSHi, gammaP2JTSHi ) == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }
  h22JTSHi  = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, 2 );
  h21JTSHi  = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, 1 );
  h20JTSHi  = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, 0 );
  h2m1JTSHi = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, -1);
  h2m2JTSHi = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, -2);

  *hlmPTSHiOutput = hlmPTSHi;

  if (debugPK) {
    out = fopen( "JWavesHi.dat", "w" );
    for ( i = 0; i < retLenHi; i++ )
    {
      fprintf( out, 
        "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
        /*timeHi.data[i]+HiSRstart, creal(h22JTSHi->data->data[i]), cimag(h22JTSHi->data->data[i]),*/
        timeHi.data[i]+HiSRstart, 
        creal(h22JTSHi->data->data[i]), cimag(h22JTSHi->data->data[i]),
        creal(h21JTSHi->data->data[i]), cimag(h21JTSHi->data->data[i]),
        creal(h20JTSHi->data->data[i]), cimag(h20JTSHi->data->data[i]),
        creal(h2m1JTSHi->data->data[i]), cimag(h2m1JTSHi->data->data[i]),
        creal(h2m2JTSHi->data->data[i]), cimag(h2m2JTSHi->data->data[i]) );
     }
     fclose( out );
  }

  sigReHi  = XLALCreateREAL8Vector( retLenHi + retLenRDPatch );
  sigImHi  = XLALCreateREAL8Vector( retLenHi + retLenRDPatch );
  if ( !sigReHi || !sigImHi || !rdMatchPoint )
  {
    XLAL_ERROR( XLAL_ENOMEM );
  }
  memset( sigReHi->data, 0, sigReHi->length * sizeof( sigReHi->data[0] ));
  memset( sigImHi->data, 0, sigImHi->length * sizeof( sigImHi->data[0] ));

  if ( combSize > tAttach )
  {
      XLALPrintError( "Function: %s, The comb size looks to be too big!!!\n", __func__ );
      XLAL_ERROR( XLAL_EINVAL );
  }
  rdMatchPoint->data[0] = combSize < tAttach ? tAttach - combSize : 0;
  rdMatchPoint->data[1] = tAttach;
  rdMatchPoint->data[2] = (retLenHi-1)*deltaTHigh/mTScaled;
  if (debugPK){
    printf("YP::comb range: %f, %f\n", 
        rdMatchPoint->data[0],rdMatchPoint->data[1]);
    printf("Stas, tAttach = %f, comb range = %f to %f \n", 
         tAttach, rdMatchPoint->data[0]+HiSRstart,
                  rdMatchPoint->data[1]+HiSRstart);
    fflush(NULL);
  }

  rdMatchPoint->data[0] -= fmod( rdMatchPoint->data[0], deltaTHigh/mTScaled );
  rdMatchPoint->data[1] -= fmod( rdMatchPoint->data[1], deltaTHigh/mTScaled );

  if (debugPK){
      printf("Stas: matching points (again) %f, %f or %f, %f \n", 
      rdMatchPoint->data[0],rdMatchPoint->data[1],
      rdMatchPoint->data[0]+HiSRstart,rdMatchPoint->data[1]+HiSRstart);
    fflush(NULL);
  }


  // FIXME Here I'll try to make a loop to check if tAttach is good or not
  //
  int CheckRDAttachment = 0;
  double S1amp = sqrt(spin1[0]*spin1[0] + spin1[1]*spin1[1] +spin1[2]*spin1[2]);


  if (CheckRDAttachment || S1amp >= 0.0){

      if (debugPK) printf("Stas: checking the RD attachment point....\n");

      int pass = 0;
      rdMatchPoint->data[0] = combSize < tAttach ? tAttach - combSize : 0;
      rdMatchPoint->data[1] = tAttach;
      rdMatchPoint->data[2] = (retLenHi-1)*deltaTHigh/mTScaled;
      rdMatchPoint->data[0] -= fmod( rdMatchPoint->data[0], deltaTHigh/mTScaled );
      rdMatchPoint->data[1] -= fmod( rdMatchPoint->data[1], deltaTHigh/mTScaled );

      REAL8 thr = 1.;
      REAL8 ratio22 = 1.0;
      REAL8 ratio2m2 = 1.0;
 
      //hJTSHi  = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, 2 );

      //printf("Stas check tmieHi-> %d, sigReHi -> %d \n", timeHi.length, sigReHi->length);
      for ( i = 0; i < retLenHi; i++ )
        {
          sigReHi->data[i] = creal(h22JTSHi->data->data[i]);
          sigImHi->data[i] = cimag(h22JTSHi->data->data[i]);
        }
      if( XLALSimCheckRDattachment(sigReHi, sigImHi, &ratio22, tAttach, 2, 2,
                    deltaTHigh, m1, m2, 0.0, 0.0, chi1J, 0.0, 0.0, chi2J,
                    &timeHi, rdMatchPoint, spinEOBApproximant, kappaJL ) == XLAL_FAILURE )
      {
          XLAL_ERROR( XLAL_EFUNC );
      }
      
      memset( sigReHi->data, 0, sigReHi->length * sizeof( sigReHi->data[0] ));
      memset( sigImHi->data, 0, sigImHi->length * sizeof( sigImHi->data[0] ));
      
      //hJTSHi  = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, -2 );

      for ( i = 0; i < retLenHi; i++ )
        {
          sigReHi->data[i] = creal(h2m2JTSHi->data->data[i]);
          sigImHi->data[i] = cimag(h2m2JTSHi->data->data[i]);
        }
      if( XLALSimCheckRDattachment(sigReHi, sigImHi, &ratio2m2, tAttach, 2, -2,
                    deltaTHigh, m1, m2, 0.0, 0.0, chi1J, 0.0, 0.0, chi2J,
                    &timeHi, rdMatchPoint, spinEOBApproximant, kappaJL ) == XLAL_FAILURE )
      {
          XLAL_ERROR( XLAL_EFUNC );
      }

      if(ratio22 <= thr && ratio2m2 <=thr){
          pass = 1;
      }

      if (pass == 0){

           if (debugPK) printf("Adjusting RD attachment point... \n");
           memset( sigReHi->data, 0, sigReHi->length * sizeof( sigReHi->data[0] ));
           memset( sigImHi->data, 0, sigImHi->length * sizeof( sigImHi->data[0] ));
            
           int found_att = XLALSimAdjustRDattachmentTime( sigReHi, sigImHi, h22JTSHi, h2m2JTSHi,  
                    &ratio22, &ratio2m2, &tAttach, thr,
                    deltaTHigh, m1, m2, 0.0, 0.0, chi1J, 0.0, 0.0, chi2J,
                    &timeHi, rdMatchPoint, spinEOBApproximant, kappaJL, combSize, tMaxOmega, tMaxAmp);
             
           if (debugPK){
             if (found_att == 1){
                 printf("we have found new attachment point tAtt = %f, with ratios %f, %f \n",
                       tAttach, ratio22, ratio2m2);
             }
               else if (found_att == 2) {
                   printf("we haven't found proper attachment point, we picked the best point at tAtt = %f\n",
                          tAttach);

               }
               else{
                 printf("we haven't found proper attachment point, best ratios are %f, %f at tAtt = %f \n",
                       ratio22, ratio2m2, tAttach);
             }
           }

      }    
      memset( sigReHi->data, 0, sigReHi->length * sizeof( sigReHi->data[0] ));
      memset( sigImHi->data, 0, sigImHi->length * sizeof( sigImHi->data[0] ));

     
  } //end of RD attachment if 

  if (debugRD){
     out = fopen( "tAttach.dat", "w" );
     fprintf( out, "%.16e    %.16e    %.16e   %.16e \n", tPeakOmega, deltaNQC, tAmpMax, tAttach); 
     fclose(out);
 
  }
  AttachParams->data[0] = tPeakOmega; 
  AttachParams->data[1] = deltaNQC; 
  AttachParams->data[2] = tAmpMax;
  AttachParams->data[3] = tAttach;
  AttachParams->data[4] = HiSRstart;

  rdMatchPoint->data[0] = combSize < tAttach ? tAttach - combSize : 0;
  rdMatchPoint->data[1] = tAttach;
  rdMatchPoint->data[2] = (retLenHi-1)*deltaTHigh/mTScaled;
  rdMatchPoint->data[0] -= fmod( rdMatchPoint->data[0], deltaTHigh/mTScaled );
  rdMatchPoint->data[1] -= fmod( rdMatchPoint->data[1], deltaTHigh/mTScaled );

  
  /*** Stas Let's try to attach RD to 2,2 mode: ***/
  for ( k = 2; k > -3; k-- )
  {
    hJTSHi  = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, k );
    for ( i = 0; i < retLenHi; i++ )
    {
      sigReHi->data[i] = creal(hJTSHi->data->data[i]);
      sigImHi->data[i] = cimag(hJTSHi->data->data[i]);
    }
    if ( XLALSimIMREOBHybridAttachRingdown( sigReHi, sigImHi, 2, k,
                deltaTHigh, m1, m2, 0.0, 0.0, chi1J, 0.0, 0.0, chi2J,
                &timeHi, rdMatchPoint, spinEOBApproximant, kappaJL )
                //&timeHi, rdMatchPoint, spinEOBApproximant, JLN )
            == XLAL_FAILURE )
    {
      XLAL_ERROR( XLAL_EFUNC );
    }

    for ( i = 0; i < (int)sigReHi->length; i++ )
    {
      hIMRJTSHi->data->data[i] = (sigReHi->data[i]) + I * (sigImHi->data[i]);
    }
    hIMRlmJTSHi = XLALSphHarmTimeSeriesAddMode( hIMRlmJTSHi, hIMRJTSHi, 2, k );
  }
  XLALSphHarmTimeSeriesSetTData( hIMRlmJTSHi, tlistRDPatchHi );
  if (debugPK){ printf("YP: J wave RD attachment done.\n"); fflush(NULL); }

  hIMR22JTSHi  = XLALSphHarmTimeSeriesGetMode( hIMRlmJTSHi, 2, 2 );
  hIMR21JTSHi  = XLALSphHarmTimeSeriesGetMode( hIMRlmJTSHi, 2, 1 );
  hIMR20JTSHi  = XLALSphHarmTimeSeriesGetMode( hIMRlmJTSHi, 2, 0 );
  hIMR2m1JTSHi = XLALSphHarmTimeSeriesGetMode( hIMRlmJTSHi, 2, -1);
  hIMR2m2JTSHi = XLALSphHarmTimeSeriesGetMode( hIMRlmJTSHi, 2, -2);
 
  *hIMRlmJTSHiOutput = hIMRlmJTSHi;
  if (debugPK){
    out = fopen( "JIMRWavesHi.dat", "w" );
    for ( i = 0; i < retLenHi + retLenRDPatch; i++ )
    {
      fprintf( out, 
        "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
        tlistRDPatchHi->data[i], 
        creal(hIMR22JTSHi->data->data[i]), cimag(hIMR22JTSHi->data->data[i]),
        creal(hIMR21JTSHi->data->data[i]), cimag(hIMR21JTSHi->data->data[i]),
        creal(hIMR20JTSHi->data->data[i]), cimag(hIMR20JTSHi->data->data[i]),
        creal(hIMR2m1JTSHi->data->data[i]), cimag(hIMR2m1JTSHi->data->data[i]),
        creal(hIMR2m2JTSHi->data->data[i]), cimag(hIMR2m2JTSHi->data->data[i]) );
     }
     fclose( out );
  
     printf("Stas: down sampling and make the complete waveform in J-frame\n");
     fflush(NULL);
  }

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  /*** attach hi sampling part and resample  */
  COMPLEX16TimeSeries *hIMRJTS2mHi    = NULL;
  for ( k = 2; k > -3; k-- )
  {
     hIMRJTS2mHi = XLALSphHarmTimeSeriesGetMode( hIMRlmJTSHi, 2, k );
     for ( i = 0; i < (int)sigReHi->length; i++ )
     {
      sigReHi->data[i] = creal(hIMRJTS2mHi->data->data[i]);
      sigImHi->data[i] = cimag(hIMRJTS2mHi->data->data[i]);
     }
     /* recycling h20PTS */
     hJTS = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, k );
     for (i = 0; i< retLenLow; i++){
         if (i*deltaT/mTScaled > rdMatchPoint->data[1]+HiSRstart) break;//Andrea
         hIMRJTS->data->data[i] = hJTS->data->data[i];
     }
      int idxRD = i; //Andrea
     spline = gsl_spline_alloc( gsl_interp_cspline, retLenHi + retLenRDPatch );
     acc    = gsl_interp_accel_alloc();   
     gsl_spline_init( spline, tlistRDPatchHi->data, sigReHi->data, retLenHi + retLenRDPatch );
//     for (i = retLenLow-4; i< retLenLow+retLenRDPatchLow; i++){
//          hIMRJTS->data->data[i] = gsl_spline_eval( spline, tlistRDPatch->data[i], acc );
      //Andrea
      for (i = idxRD; i< retLenLow+retLenRDPatchLow; i++){
          if( i*deltaT/mTScaled <= tlistRDPatchHi->data[ retLenHi + retLenRDPatch - 1]) {
              hIMRJTS->data->data[i] = gsl_spline_eval( spline, tlistRDPatch->data[i], acc );
          }
          else {
              hIMRJTS->data->data[i] = 0.;
          }
      }
     gsl_spline_free(spline);
     gsl_interp_accel_free(acc);

     spline = gsl_spline_alloc( gsl_interp_cspline, retLenHi + retLenRDPatch );
     acc    = gsl_interp_accel_alloc();     
     gsl_spline_init( spline, 
              tlistRDPatchHi->data, sigImHi->data, retLenHi + retLenRDPatch );
     for (i = idxRD; i< retLenLow+retLenRDPatchLow; i++){
         if( i*deltaT/mTScaled <= tlistRDPatchHi->data[ retLenHi + retLenRDPatch - 1]) {
             hIMRJTS->data->data[i] += I * gsl_spline_eval( spline, tlistRDPatch->data[i], acc );
         }
         else {
             hIMRJTS->data->data[i] += I*0.;
         }
     }

     gsl_spline_free(spline);
     gsl_interp_accel_free(acc);

     hIMRlmJTS = XLALSphHarmTimeSeriesAddMode( hIMRlmJTS, hIMRJTS, 2, k );
     /*if (k==1) exit(0);*/
  }
  XLALSphHarmTimeSeriesSetTData( hIMRlmJTS, tlistRDPatch );
  if (debugPK){ printf("Stas: J-wave with RD  generated.\n"); fflush(NULL); }

  hIMR22JTS  = XLALSphHarmTimeSeriesGetMode( hIMRlmJTS, 2, 2 );
  hIMR21JTS  = XLALSphHarmTimeSeriesGetMode( hIMRlmJTS, 2, 1 );
  hIMR20JTS  = XLALSphHarmTimeSeriesGetMode( hIMRlmJTS, 2, 0 );
  hIMR2m1JTS = XLALSphHarmTimeSeriesGetMode( hIMRlmJTS, 2, -1);
  hIMR2m2JTS = XLALSphHarmTimeSeriesGetMode( hIMRlmJTS, 2, -2);

  if (debugPK){
     out = fopen( "JIMRWaves.dat", "w" );
     for ( i = 0; i < retLenLow + retLenRDPatchLow; i++ )
     {
        fprintf( out, 
          "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
          i*deltaT/mTScaled, 
          creal(hIMR22JTS->data->data[i]), cimag(hIMR22JTS->data->data[i]),
          creal(hIMR21JTS->data->data[i]), cimag(hIMR21JTS->data->data[i]),
          creal(hIMR20JTS->data->data[i]), cimag(hIMR20JTS->data->data[i]),
          creal(hIMR2m1JTS->data->data[i]), cimag(hIMR2m1JTS->data->data[i]),
          creal(hIMR2m2JTS->data->data[i]), cimag(hIMR2m2JTS->data->data[i]) );
     }
     fclose( out );
     printf("YP: IMR J wave written to file.\n");
     fflush(NULL);
  }

  /**** First atempt to compute Euler Angles and rotate to I frame */  
  gamJtoI = atan2(JframeEz[1], -JframeEz[0]);
  betJtoI = acos(JframeEz[2]);
  alJtoI = atan2(JframeEy[2], -JframeEx[2]);

  if (debugPK){
    printf("Stas: J->I EA = %.16e, %.16e, %.16e \n", alJtoI, betJtoI, gamJtoI);
    fflush(NULL); 
  }
 
  /*COMPLEX16TimeSeries *hITS    = XLALCreateCOMPLEX16TimeSeries( "HJ",     &tc, 0.0, deltaT, &lalStrainUnit, retLen );*/
  COMPLEX16TimeSeries *hIMR22ITS  = NULL;
  COMPLEX16TimeSeries *hIMR21ITS  = NULL;
  COMPLEX16TimeSeries *hIMR20ITS  = NULL;
  COMPLEX16TimeSeries *hIMR2m1ITS = NULL;
  COMPLEX16TimeSeries *hIMR2m2ITS = NULL;


  //COMPLEX16TimeSeries *hIMR22ITS  = XLALCreateCOMPLEX16TimeSeries( "HIMRI_22",
  //  &tc, 0.0, deltaT, &lalStrainUnit, retLenLow + retLenRDPatchLow );
  REAL8TimeSeries *alpI        = XLALCreateREAL8TimeSeries( "alphaJ2I",
        &tc, 0.0, deltaT, &lalStrainUnit, retLenLow + retLenRDPatchLow );
  REAL8TimeSeries *betI        = XLALCreateREAL8TimeSeries( "betaaJ2I",
        &tc, 0.0, deltaT, &lalStrainUnit, retLenLow + retLenRDPatchLow );
  REAL8TimeSeries *gamI        = XLALCreateREAL8TimeSeries( "gammaJ2I",
        &tc, 0.0, deltaT, &lalStrainUnit, retLenLow + retLenRDPatchLow );

  for (i=0; i< retLenLow + retLenRDPatchLow; i++){
      alpI->data->data[i] = -alJtoI;
      betI->data->data[i] = betJtoI;
      gamI->data->data[i] = -gamJtoI;
  }

  hIMRlmITS = XLALSphHarmTimeSeriesAddMode( hIMRlmITS, hIMR22JTS, 2, 2 );
  hIMRlmITS = XLALSphHarmTimeSeriesAddMode( hIMRlmITS, hIMR21JTS, 2, 1 );
  hIMRlmITS = XLALSphHarmTimeSeriesAddMode( hIMRlmITS, hIMR20JTS, 2, 0 );
  hIMRlmITS = XLALSphHarmTimeSeriesAddMode( hIMRlmITS, hIMR2m1JTS, 2, -1 );
  hIMRlmITS = XLALSphHarmTimeSeriesAddMode( hIMRlmITS, hIMR2m2JTS, 2, -2 );
  XLALSphHarmTimeSeriesSetTData( hIMRlmITS, tlistRDPatch );

  if (debugPK){ printf("Rotation to inertial I-frame\n"); fflush(NULL); }
  if ( XLALSimInspiralPrecessionRotateModes( hIMRlmITS, 
                alpI, betI, gamI ) == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }
  hIMR22ITS  = XLALSphHarmTimeSeriesGetMode( hIMRlmITS, 2, 2 );
  hIMR21ITS  = XLALSphHarmTimeSeriesGetMode( hIMRlmITS, 2, 1 );
  hIMR20ITS  = XLALSphHarmTimeSeriesGetMode( hIMRlmITS, 2, 0 );
  hIMR2m1ITS = XLALSphHarmTimeSeriesGetMode( hIMRlmITS, 2, -1);
  hIMR2m2ITS = XLALSphHarmTimeSeriesGetMode( hIMRlmITS, 2, -2);

  if (debugPK){
    out = fopen( "IWaves.dat", "w" );
    for ( i = 0; i < retLenLow + retLenRDPatchLow; i++ )
    {
      fprintf( out, 
        "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
        /*timeHi.data[i]+HiSRstart, creal(h22JTSHi->data->data[i]), cimag(h22JTSHi->data->data[i]),*/
        tlistRDPatch->data[i], 
        creal(hIMR22ITS->data->data[i]), cimag(hIMR22ITS->data->data[i]),
        creal(hIMR21ITS->data->data[i]), cimag(hIMR21ITS->data->data[i]),
        creal(hIMR20ITS->data->data[i]), cimag(hIMR20ITS->data->data[i]),
        creal(hIMR2m1ITS->data->data[i]), cimag(hIMR2m1ITS->data->data[i]),
        creal(hIMR2m2ITS->data->data[i]), cimag(hIMR2m2ITS->data->data[i]) );
    }
    fclose( out );
  }

  Y22 = XLALSpinWeightedSphericalHarmonic( inc, phiC, -2, 2, 2 );
  Y2m2 = XLALSpinWeightedSphericalHarmonic( inc, phiC, -2, 2, -2 );
  Y21 = XLALSpinWeightedSphericalHarmonic( inc, phiC, -2, 2, 1 );
  Y2m1 = XLALSpinWeightedSphericalHarmonic( inc, phiC, -2, 2, -1 );
  Y20 = XLALSpinWeightedSphericalHarmonic( inc, phiC, -2, 2, 0 );
 
  for ( i = 0; i < (INT4)hIMR22ITS->data->length; i++ )
  {
    x11 = Y22*hIMR22ITS->data->data[i] + Y21*hIMR21ITS->data->data[i] 
          + Y20*hIMR20ITS->data->data[i] + Y2m1*hIMR2m1ITS->data->data[i] 
          + Y2m2*hIMR2m2ITS->data->data[i];    
    hPlusTS->data->data[i]  = creal(x11);
    hCrossTS->data->data[i] = -cimag(x11);     
  }
    
  /* Point the output pointers to the relevant time series and return */
  (*hplus)  = hPlusTS;
  (*hcross) = hCrossTS;
     
  if (debugPK){
    printf("plus and cross are computed, freeing the memory\n"); 
    fflush(NULL); 
  }
  /*hplus =  XLALSphHarmTimeSeriesGetMode( hIMRlmJTS, 2, 2 );
  hcross =  XLALSphHarmTimeSeriesGetMode( hIMRlmJTS, 2, 2 );*/
  
  XLALDestroyREAL8TimeSeries(alphaI2PTS);
  XLALDestroyREAL8TimeSeries(betaI2PTS);
  XLALDestroyREAL8TimeSeries(gammaI2PTS);
  XLALDestroyREAL8TimeSeries(alphaP2JTS);
  XLALDestroyREAL8TimeSeries(betaP2JTS);
  XLALDestroyREAL8TimeSeries(gammaP2JTS);

  if(debugPK){ printf("Memory cleanup 1 done.\n"); fflush(NULL); }
  XLALDestroyREAL8TimeSeries(alphaI2PTSHi);
  XLALDestroyREAL8TimeSeries(betaI2PTSHi);
  XLALDestroyREAL8TimeSeries(gammaI2PTSHi);
  XLALDestroyREAL8TimeSeries(alphaP2JTSHi);
  XLALDestroyREAL8TimeSeries(betaP2JTSHi);
  XLALDestroyREAL8TimeSeries(gammaP2JTSHi);
 

  ////////////////////////////////////////


  XLALDestroyCOMPLEX16TimeSeries(h22TS); 
  XLALDestroyCOMPLEX16TimeSeries(h21TS); 
  XLALDestroyCOMPLEX16TimeSeries(h20TS);  
  XLALDestroyCOMPLEX16TimeSeries(h2m1TS); 
  XLALDestroyCOMPLEX16TimeSeries(h2m2TS); 
  
  //XLALDestroyCOMPLEX16TimeSeries(h22PTS); 
  //XLALDestroyCOMPLEX16TimeSeries(h21PTS);
  //XLALDestroyCOMPLEX16TimeSeries(h20PTS);
  //XLALDestroyCOMPLEX16TimeSeries(h2m1PTS);
  //XLALDestroyCOMPLEX16TimeSeries(h2m2PTS);
  
  //printf("Memory cleanup 2 done.\n"); fflush(NULL); 
  XLALDestroySphHarmTimeSeries(hlmPTS);
  // printf("Memory cleanup 2.5 done.\n"); fflush(NULL); 
  XLALDestroySphHarmTimeSeries(hIMRlmJTS);
  //printf("Memory cleanup 2.6 done.\n"); fflush(NULL); 
  //XLALDestroySphHarmTimeSeries(hIMRlmJTS);
  //XLALDestroySphHarmTimeSeries(hIMRlmITS);
  //XLALDestroyCOMPLEX16TimeSeries(hIMR22ITS); 
  //XLALDestroyCOMPLEX16TimeSeries(hIMR21ITS); 
  //XLALDestroyCOMPLEX16TimeSeries(hIMR20ITS); 
  //XLALDestroyCOMPLEX16TimeSeries(hIMR2m1ITS); 
  //XLALDestroyCOMPLEX16TimeSeries(hIMR2m2ITS); 
  
  XLALDestroyCOMPLEX16TimeSeries(XLALSphHarmTimeSeriesGetMode( hIMRlmITS, 2, 2 ));
  XLALDestroyCOMPLEX16TimeSeries(XLALSphHarmTimeSeriesGetMode( hIMRlmITS, 2, 1 ));
  XLALDestroyCOMPLEX16TimeSeries(XLALSphHarmTimeSeriesGetMode( hIMRlmITS, 2, 0 ));
  XLALDestroyCOMPLEX16TimeSeries(XLALSphHarmTimeSeriesGetMode( hIMRlmITS, 2, -1 ));
  XLALDestroyCOMPLEX16TimeSeries(XLALSphHarmTimeSeriesGetMode( hIMRlmITS, 2, -2 ));
  //printf("Memory cleanup 2.7 done.\n"); fflush(NULL); 
 
  //XLALDestroyCOMPLEX16TimeSeries(hIMR22JTS); 
  //XLALDestroyCOMPLEX16TimeSeries(hIMR21JTS); 
  //XLALDestroyCOMPLEX16TimeSeries(hIMR20JTS); 
  //XLALDestroyCOMPLEX16TimeSeries(hIMR2m1JTS); 
  //XLALDestroyCOMPLEX16TimeSeries(hIMR2m2JTS); 
  
  XLALDestroyCOMPLEX16TimeSeries(hIMRJTS);
  if(debugPK){ printf("Memory cleanup 2 done.\n"); fflush(NULL); }
  //printf("Memory cleanup 2.8 done.\n"); fflush(NULL); 
  


  XLALDestroyCOMPLEX16TimeSeries(h22TSHi);  
  XLALDestroyCOMPLEX16TimeSeries(h21TSHi); 
  XLALDestroyCOMPLEX16TimeSeries(h20TSHi); 
  XLALDestroyCOMPLEX16TimeSeries(h2m1TSHi); 
  XLALDestroyCOMPLEX16TimeSeries(h2m2TSHi);
  
  //XLALDestroyCOMPLEX16TimeSeries(h22PTSHi); 
  //XLALDestroyCOMPLEX16TimeSeries(h21PTSHi);
  //XLALDestroyCOMPLEX16TimeSeries(h20PTSHi);
  //XLALDestroyCOMPLEX16TimeSeries(h2m1PTSHi);
  //XLALDestroyCOMPLEX16TimeSeries(h2m2PTSHi);
  //XLALDestroyCOMPLEX16TimeSeries(hIMR22JTSHi); 
  //XLALDestroyCOMPLEX16TimeSeries(hIMR21JTSHi); 
  //XLALDestroyCOMPLEX16TimeSeries(hIMR20JTSHi); 
  //XLALDestroyCOMPLEX16TimeSeries(hIMR2m1JTSHi); 
  //XLALDestroyCOMPLEX16TimeSeries(hIMR2m2JTSHi); 
  
  XLALDestroyCOMPLEX16TimeSeries(hIMRJTSHi);
    

  XLALDestroyREAL8Vector( values );
  XLALDestroyREAL8Vector( dvalues );
  XLALDestroyREAL8Vector( sigmaStar );
  XLALDestroyREAL8Vector( sigmaKerr );
  XLALDestroyREAL8Vector( sigReHi );
  XLALDestroyREAL8Vector( sigImHi );
  //XLALDestroyREAL8Vector( omegaHi );
    
  if(debugPK){ printf("Memory cleanup 3 done.\n"); fflush(NULL); }
  //printf("Memory cleanup 3 done.\n"); fflush(NULL); 
  XLALAdaptiveRungeKutta4Free(integrator);
  XLALDestroyREAL8Array( dynamics );
  //XLALDestroyREAL8Array( dynamicsHi );
    
  XLALDestroyREAL8Vector( LN_x );
  XLALDestroyREAL8Vector( LN_y );
  XLALDestroyREAL8Vector( LN_z );
  XLALDestroyREAL8Vector( LN_xHi );
  XLALDestroyREAL8Vector( LN_yHi );
  XLALDestroyREAL8Vector( LN_zHi );
  XLALDestroyREAL8Vector( Alpha );
  XLALDestroyREAL8Vector( Beta );
  XLALDestroyREAL8Vector( AlphaHi );
  XLALDestroyREAL8Vector( BetaHi );
  XLALDestroyREAL8Vector( Gamma );
  XLALDestroyREAL8Vector( GammaHi );
  XLALDestroyREAL8TimeSeries( alpI);
  XLALDestroyREAL8TimeSeries( betI);
  XLALDestroyREAL8TimeSeries( gamI);


  XLALDestroyREAL8Vector( rdMatchPoint );
  //XLALDestroyREAL8Vector( AttachParams );
  XLALDestroyREAL8Vector( tmpValues2 );
  //printf("Memory cleanup 4 done.\n"); fflush(NULL); 
    
    
  //XLALDestroyREAL8Vector( tlist );
  //XLALDestroyREAL8Vector( tlistHi );
  //XLALDestroyREAL8Vector( tlistRDPatch );
  //XLALDestroyREAL8Vector( tlistRDPatchHi );
   
  //printf("Memory cleanup ALL done.\n"); fflush(NULL); 
  if(debugPK){ printf("Memory cleanup ALL done.\n"); fflush(NULL); }

  /* FIXME: Temporary code to convert REAL8Array to REAL8Vector because SWIG
   *        doesn't seem to like REAL8Array */
  REAL8Vector *tmp_vec;
  tmp_vec = XLALCreateREAL8Vector(dynamicsHi->dimLength->data[0] * dynamicsHi->dimLength->data[1]);
  UINT4 tmp_idx_ii;
  for (tmp_idx_ii=0; tmp_idx_ii < tmp_vec->length; tmp_idx_ii++)
  {
      tmp_vec->data[tmp_idx_ii] = dynamicsHi->data[tmp_idx_ii];
  }
  *dynHi = tmp_vec;
  XLALDestroyREAL8Array( dynamicsHi );

  return XLAL_SUCCESS;
}


