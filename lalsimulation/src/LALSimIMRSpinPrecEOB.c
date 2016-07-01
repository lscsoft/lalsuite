/*
*  Copyright (C) 2011 Craig Robinson, Enrico Barausse, Yi Pan,
*                2014 Prayush Kumar, Stas Babak, Andrea Taracchini (Precessing EOB)
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
#ifndef _LALSIMIMRSPINPRECEOB_C
#define _LALSIMIMRSPINPRECEOB_C
/**
 * @addtogroup LALSimIMRSpinPrecEOB_c
 *
 * @author Craig Robinson, Yi Pan, Prayush Kumar, Stas Babak, Andrea Taracchini
 *
 * \brief Functions for producing SEOBNRv3 waveforms for
 * precessing binaries of spinning compact objects, as described
 * in Pan et al. ( PRD 89, 084006 (2014) ) == YPP.
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
#include "LALSimBlackHoleRingdownPrec.c"
#include "LALSimFindAttachTime.h"

/* Include all the static function files we need */
#include "LALSimIMREOBNQCCorrection.c"
#include "LALSimInspiraldEnergyFlux.c"
#include "LALSimIMREOBNewtonianMultipole.c"
#include "LALSimIMRSpinEOBFactorizedWaveformPrec.c"
#include "LALSimIMREOBHybridRingdownPrec.c"
#include "LALSimIMRSpinEOBAuxFuncs.c"
#include "LALSimIMRSpinEOBFactorizedWaveformCoefficientsPrec.c"
#include "LALSimIMRSpinEOBHamiltonianPrec.c"
#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBFactorizedFluxPrec.c"
#include "LALSimIMRSpinEOBHcapNumericalDerivativePrec.c"
#include "LALSimIMRSpinAlignedEOBHcapDerivative.c"
#include "LALSimIMRSpinEOBInitialConditions.c"
#include "LALSimIMRSpinEOBInitialConditionsPrec.c"


#define debugOutput 0

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define FREE_EVERYTHING \
if ( hlmPTS != NULL ) \
XLALDestroySphHarmTimeSeries(hlmPTS); \
if ( hIMRlmJTS != NULL ) \
XLALDestroySphHarmTimeSeries(hIMRlmJTS); \
if ( hIMRlmITS != NULL ) \
XLALDestroySphHarmTimeSeries(hIMRlmITS); \
\
if ( dynamics != NULL ) \
XLALDestroyREAL8Array( dynamics ); \
if ( dynamicsV2 != NULL ) \
XLALDestroyREAL8Array( dynamicsV2 ); \
if ( dynamicsV2Hi != NULL ) \
XLALDestroyREAL8Array( dynamicsV2Hi ); \
if ( dynamicsHi != NULL ) \
XLALDestroyREAL8Array( dynamicsHi ); \
\
XLALDestroyREAL8Vector( valuesV2 ); \
XLALDestroyREAL8Vector( values ); \
XLALDestroyREAL8Vector( dvalues ); \
\
XLALDestroyREAL8Vector( sigmaStar ); \
XLALDestroyREAL8Vector( sigmaKerr ); \
\
XLALDestroyREAL8Vector( rdMatchPoint ); \
\
XLALDestroyREAL8Vector( Alpha ); \
XLALDestroyREAL8Vector( Beta ); \
XLALDestroyREAL8Vector( Gamma ); \
XLALDestroyREAL8Vector( AlphaHi ); \
XLALDestroyREAL8Vector( BetaHi ); \
XLALDestroyREAL8Vector( GammaHi ); \
\
XLALDestroyREAL8Vector( sigReHi ); \
XLALDestroyREAL8Vector( sigImHi ); \
\
XLALDestroyREAL8TimeSeries(alphaI2PTS); \
XLALDestroyREAL8TimeSeries(betaI2PTS); \
XLALDestroyREAL8TimeSeries(gammaI2PTS); \
\
XLALDestroyREAL8TimeSeries(alphaI2PTSHi); \
XLALDestroyREAL8TimeSeries(betaI2PTSHi); \
XLALDestroyREAL8TimeSeries(gammaI2PTSHi); \
\
XLALDestroyREAL8TimeSeries(alphaP2JTS); \
XLALDestroyREAL8TimeSeries(betaP2JTS); \
XLALDestroyREAL8TimeSeries(gammaP2JTS); \
\
XLALDestroyREAL8TimeSeries(alphaP2JTSHi); \
XLALDestroyREAL8TimeSeries(betaP2JTSHi); \
XLALDestroyREAL8TimeSeries(gammaP2JTSHi); \
\
XLALDestroyREAL8TimeSeries(alpI); \
XLALDestroyREAL8TimeSeries(betI); \
XLALDestroyREAL8TimeSeries(gamI); \
\
XLALDestroyCOMPLEX16TimeSeries(h22TS); \
XLALDestroyCOMPLEX16TimeSeries(h21TS); \
XLALDestroyCOMPLEX16TimeSeries(h20TS); \
XLALDestroyCOMPLEX16TimeSeries(h2m1TS); \
XLALDestroyCOMPLEX16TimeSeries(h2m2TS); \
\
XLALDestroyCOMPLEX16TimeSeries(h22TSHi); \
XLALDestroyCOMPLEX16TimeSeries(h21TSHi); \
XLALDestroyCOMPLEX16TimeSeries(h20TSHi); \
XLALDestroyCOMPLEX16TimeSeries(h2m1TSHi); \
XLALDestroyCOMPLEX16TimeSeries(h2m2TSHi); \
\
XLALDestroyCOMPLEX16TimeSeries(hIMRJTS); \
XLALDestroyCOMPLEX16TimeSeries(hIMRJTSHi); \
\
XLALAdaptiveRungeKutta4Free(integrator); \
\
XLALDestroyREAL8TimeSeries(hPlusTS); \
XLALDestroyREAL8TimeSeries(hCrossTS);


#define FREE_SPHHARM \
XLALDestroyREAL8Vector( tlist ); \
XLALDestroyREAL8Vector( tlistHi ); \
XLALDestroyREAL8Vector( timeJFull ); \
XLALDestroyREAL8Vector( timeIFull ); \
XLALDestroyREAL8Vector( tlistRDPatch ); \
XLALDestroyREAL8Vector( tlistRDPatchHi );

#define PRINT_PARAMS do { \
XLALPrintError("--approximant SEOBNRv3 --f-min %.16e --m1 %.16e --m2 %.16e --spin1x %.16e --spin1y %.16e --spin1z %.16e  --spin2x %.16e --spin2y %.16e --spin2z %.16e --inclination %.16e --distance %.16e --phiRef %.16e --sample-rate %.16e\n",\
               fMin, m1SI/LAL_MSUN_SI, m2SI/LAL_MSUN_SI, INspin1x, INspin1y, INspin1z, INspin2x, INspin2y, INspin2z, inc, r/(1e6 * LAL_PC_SI), phiC, 1./deltaT);\
} while(0);

/**
 * Computes RHS of ODE for gamma. Eq. 10 of PRD 89, 084006 (2014)
 */
static double f_alphadotcosi( double x, void * inparams )
{
	PrecEulerAnglesIntegration* params = (PrecEulerAnglesIntegration*) inparams;

	REAL8 alphadot = gsl_spline_eval_deriv( params->alpha_spline, x, params->alpha_acc );
	REAL8 beta = gsl_spline_eval( params->beta_spline, x, params->beta_acc );

	return -1. * alphadot * cos(beta);

}

/**
 * Given the trajectory in an inertial frame, this computes Euler angles
 * of the roation from the inertial frame to the minimal-rotation frame
 * that co-precesses with LN(t) = rvec(t) x rdotvec(t)
 */
static int EulerAnglesI2P(REAL8Vector *Alpha, /**<< output: alpha Euler angle */
                 REAL8Vector *Beta, /**<< output: beta Euler angle */
                 REAL8Vector *Gamma, /**<< output: gamma Euler angle */
                 INT4 *phaseCounterA, /**<< output: counter for unwrapping of alpha */
                 INT4 *phaseCounterB, /**<< output: counter for unwrapping of beta */
                 const REAL8Vector tVec, /**<< time series */
                 const REAL8Vector posVecx, /**<< x time series */
                 const REAL8Vector posVecy, /**<< y time series */
                 const REAL8Vector posVecz, /**<< z time series */
                 const UINT4 retLenLow, /**<< Array length of the trajectory */
                 const REAL8 InitialAlpha, /**<< Initial alpha (used only if flag_highSR=1) */
                 const REAL8 InitialGamma, /**<< Initial gamma */
                 UINT4 flag_highSR /**<< Flag to indicate whether one is analyzing the high SR trajectory */) {
    UINT4 i = 0;
    REAL8Vector *LN_x = NULL, *LN_y = NULL, *LN_z = NULL;
    REAL8 tmpR[3], tmpRdot[3], magLN;
    REAL8 precEulerresult = 0, precEulererror = 0;
    gsl_integration_workspace * precEulerw = gsl_integration_workspace_alloc (1000);
    gsl_function precEulerF;
    PrecEulerAnglesIntegration precEulerparams;
    REAL8 inGamma = InitialGamma; 
    
    LN_x = XLALCreateREAL8Vector( retLenLow );
    LN_y = XLALCreateREAL8Vector( retLenLow );
    LN_z = XLALCreateREAL8Vector( retLenLow );
    
    gsl_spline *x_spline = gsl_spline_alloc( gsl_interp_cspline, retLenLow );
    gsl_spline *y_spline = gsl_spline_alloc( gsl_interp_cspline, retLenLow );
    gsl_spline *z_spline = gsl_spline_alloc( gsl_interp_cspline, retLenLow );
    
    gsl_interp_accel *x_acc    = gsl_interp_accel_alloc();
    gsl_interp_accel *y_acc    = gsl_interp_accel_alloc();
    gsl_interp_accel *z_acc    = gsl_interp_accel_alloc();
    
    gsl_spline_init( x_spline, tVec.data, posVecx.data, retLenLow );
    gsl_spline_init( y_spline, tVec.data, posVecy.data, retLenLow );
    gsl_spline_init( z_spline, tVec.data, posVecz.data, retLenLow );

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

        /*  Eq. 19 of PRD 89, 084006 (2014) */
        /*  Also unwrap the two angles */
        if (fabs(LN_x->data[i]) <= 1.e-10 && fabs(LN_y->data[i]) <=1.e-10){
            Alpha->data[i] = 0.0;
            inGamma = 0.0;
        } else {
            Alpha->data[i] = atan2( LN_y->data[i], LN_x->data[i] )
                             +  *phaseCounterA * LAL_TWOPI;
            if (i==0 && flag_highSR != 1){
                inGamma = -Alpha->data[i];
            }
        }
        
        if( i>0 && Alpha->data[i] - Alpha->data[i-1] > 5. ) {
            *phaseCounterA = *phaseCounterA - 1;
            Alpha->data[i] -= LAL_TWOPI;
        }
        else if ( i && Alpha->data[i] - Alpha->data[i-1] < -5. ) {
            *phaseCounterA = *phaseCounterA + 1;
            Alpha->data[i] += LAL_TWOPI;
        }
        if (LN_z->data[i] >1.) {
            LN_z->data[i] = 1.;
        }
        if (LN_z->data[i] <-1.) {
            LN_z->data[i] = -1.;
        }
        if ( flag_highSR == 1) {
            Alpha->data[i] -= (Alpha->data[0] - InitialAlpha);
        }

        if (fabs(1.0 - LN_z->data[i]) < 1.e-12){
            REAL8 LN_xy;
            LN_xy = sqrt(LN_x->data[i]*LN_x->data[i] + LN_y->data[i]*LN_y->data[i]);
            //LN_z->data[i] = sqrt(1.0 - LN_xy*LN_xy);
            Beta->data[i] = atan2(LN_xy, LN_z->data[i]);
            //printf("here   ");
        }else{
            Beta->data[i] = acos( LN_z->data[i] );
        }
        if( i>0 && Beta->data[i] > Beta->data[i-1] ) {
            *phaseCounterB = *phaseCounterB - 1;
        }
    }
    /* Integrate \dot{\alpha} \cos{\beta} to get the final Euler angle
     Eq. 20 of PRD 89, 084006 (2014) */
    gsl_spline_init( x_spline, tVec.data, Alpha->data, retLenLow );
    gsl_spline_init( y_spline, tVec.data, Beta->data, retLenLow );
        
    precEulerparams.alpha_spline = x_spline;
    precEulerparams.alpha_acc    = x_acc;
    precEulerparams.beta_spline  = y_spline;
    precEulerparams.beta_acc     = y_acc;
    
    precEulerF.function = &f_alphadotcosi;
    precEulerF.params   = &precEulerparams;
    
    for( i = 0; i < retLenLow; i++ )
    {
        //if( i==0 ) { Gamma->data[i] = InitialGamma; }
        if( i==0 ) { Gamma->data[i] = inGamma; }
        else
        {
            gsl_integration_qags (&precEulerF, tVec.data[i-1], tVec.data[i], 1e-9, 1e-9, 1000, precEulerw, &precEulerresult, &precEulererror);
            Gamma->data[i] = Gamma->data[i-1] + precEulerresult;
        }
    }
    gsl_integration_workspace_free( precEulerw );
    gsl_spline_free( x_spline );
    gsl_spline_free( y_spline );
    gsl_spline_free( z_spline );
    gsl_interp_accel_free( x_acc );
    gsl_interp_accel_free( y_acc );
    gsl_interp_accel_free( z_acc );
    XLALDestroyREAL8Vector( LN_x );
    XLALDestroyREAL8Vector( LN_y );
    XLALDestroyREAL8Vector( LN_z );
    return XLAL_SUCCESS;
}

/**
 * Given Euler angles to go from initial inertial frame to precessing frama
 * and the LNhat vector, this functions computes the Euler angles to 
 * go from the precessing frame to the frame of the total angular
 * momentum
 */
static void EulerAnglesP2J(
                REAL8 *aP2J, /**<< alpha Euler angle from precessing to final-J frame */
                REAL8 *bP2J, /**<< beta Euler angle from precessing to final-J frame */
                REAL8 *gP2J, /**<< gamma Euler angle from precessing to final-J frame */
                const REAL8 aI2P, /**<< alpha Euler angle from inertial to precessing frame */
                const REAL8 bI2P, /**<< beta Euler angle from inertial to precessing frame */
                const REAL8 gI2P, /**<< gamma Euler angle from inertial to precessing frame */
                const REAL8 LNhx, /**<< x component of LNhat */
                const REAL8 LNhy, /**<< y component of LNhat */
                const REAL8 LNhz, /**<< z component of LNhat */
                const REAL8 JframeEx[], /**<< x-axis of the total-angular-momentum frame */
                const REAL8 JframeEy[], /**<< y-axis of the total-angular-momentum frame */
                const REAL8 JframeEz[] /**<< z-axis of the total-angular-momentum frame */
) {
    REAL8 LframeEx[3] = {0,0,0}, LframeEy[3] = {0,0,0}, LframeEz[3] = {0,0,0};
    LframeEx[0] =  cos(aI2P)*cos(bI2P)*cos(gI2P) - sin(aI2P)*sin(gI2P);
    LframeEx[1] =  sin(aI2P)*cos(bI2P)*cos(gI2P) + cos(aI2P)*sin(gI2P);
    LframeEx[2] = -sin(bI2P)*cos(gI2P);
    LframeEy[0] = -cos(aI2P)*cos(bI2P)*sin(gI2P) - sin(aI2P)*cos(gI2P);
    LframeEy[1] = -sin(aI2P)*cos(bI2P)*sin(gI2P) + cos(aI2P)*cos(gI2P);
    LframeEy[2] =  sin(bI2P)*sin(gI2P);
    LframeEz[0] =  LNhx;
    LframeEz[1] =  LNhy;
    LframeEz[2] =  LNhz;
    REAL8 normJ, normLz;
    normJ = JframeEz[0]*JframeEz[0]+JframeEz[1]*JframeEz[1]+JframeEz[2]*JframeEz[2];
    normLz = LframeEz[0]*LframeEz[0]+LframeEz[1]*LframeEz[1]+LframeEz[2]*LframeEz[2];
    *aP2J = atan2(JframeEz[0]*LframeEy[0]+JframeEz[1]*LframeEy[1]+JframeEz[2]*LframeEy[2],
                 JframeEz[0]*LframeEx[0]+JframeEz[1]*LframeEx[1]+JframeEz[2]*LframeEx[2]);
    *bP2J = acos( JframeEz[0]*LframeEz[0]+JframeEz[1]*LframeEz[1]+JframeEz[2]*LframeEz[2]);
    if (*bP2J < 1.e-4){
        *bP2J = acos((JframeEz[0]*LframeEz[0]+JframeEz[1]*LframeEz[1]+JframeEz[2]*LframeEz[2])/sqrt(normJ*normLz));

    }

    *gP2J = atan2(  JframeEy[0]*LframeEz[0]+JframeEy[1]*LframeEz[1]+JframeEy[2]*LframeEz[2],
                 -(JframeEx[0]*LframeEz[0]+JframeEx[1]*LframeEz[1]+JframeEx[2]*LframeEz[2]));
    
    /* I2P Euler angles are stored only for debugging purposes */
    if ( fabs(*bP2J-LAL_PI) < 1.e-10){
        *gP2J = 0.0;
        *aP2J = atan2( JframeEx[1], JframeEx[0]);
    }
    
    if ( fabs(*bP2J) < 1.e-10){
        *gP2J = 0.0;
        *aP2J = atan2( JframeEx[1], JframeEx[0]);
    }
}


/**
 * Stopping conditions for dynamics integration for SEOBNRv3
 */
UNUSED static int
XLALEOBSpinPrecStopConditionBasedOnPR(double UNUSED t,
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
    if (debugPK){  XLAL_PRINT_INFO("XLALEOBSpinPrecStopConditionBasedOnPR:: r = %e %e\n", sqrt(r2), omega);}
    if (debugPK){  XLAL_PRINT_INFO("XLALEOBSpinPrecStopConditionBasedOnPR:: values = %e %e %e %e %e %e\n", values[6], values[7], values[8], values[9], values[10], values[11]);}
    if (debugPK){  XLAL_PRINT_INFO("XLALEOBSpinPrecStopConditionBasedOnPR:: dvalues = %e %e %e %e %e %e\n",dvalues[6], dvalues[7], dvalues[8], dvalues[9], dvalues[10], dvalues[11]);}
  REAL8 rdot;
  // this is d(r)/dt obtained by differentiating r2 (see above)
  rdot = inner_product(rdotVec, r) / sqrt(r2);
  // This is d/dt(pDotr) see pDotr above.
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
		  if(debugPK){XLAL_PRINT_INFO("\n  isnan reached. r2 = %f\n", r2); fflush(NULL); }
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
      XLAL_PRINT_INFO("\n Integration stopping, p_r pointing outwards -- out-spiraling!\n");
      fflush(NULL);
    }
    return 1;
  }

  /* Terminate if rdot is >0 (OUTspiraling) for separation <4M */
  if ( r2 < 16 && rdot >= 0  )
  {
    if(debugPK){
      XLAL_PRINT_INFO("\n Integration stopping, dr/dt>0 -- out-spiraling!\n");
      fflush(NULL);
    }
    return 1;
  }

  /* Terminate if dp_R/dt > 0, i.e. radial momentum is increasing for separation <2M */
  if(r2 < 4. && prDot > 0. )
  {
    if(debugPK){
      XLAL_PRINT_INFO("\n Integration stopping as prDot = %lf at r = %lf\n",
            prDot, sqrt(r2));
      fflush(NULL);
    }
    return 1;
  }

    if(r2 < 16. && ( sqrt(values[3]*values[3] + values[4]*values[4] + values[5]*values[5]) > 10. )) {
        if(debugPK)XLAL_PRINT_INFO("\n Integration stopping |pvec|> 10\n");
        fflush(NULL);
        return 1;
    }

    if(r2 < 16. && ( sqrt(values[3]*values[3] + values[4]*values[4] + values[5]*values[5]) < 1.e-10 )) {
        if(debugPK)XLAL_PRINT_INFO("\n Integration stopping |pvec|<1e-10\n");
        fflush(NULL);
        return 1;
    }

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
      XLAL_PRINT_INFO("\n Integration stopping, omega reached second extremum\n");
      fflush(NULL);
    }
    return 1;
  }

  /* If Momega did not evolve above 0.01 even though r < 4 or omega<0.14 for r<2, break */
  if((r2 < 16. && omega < 0.04) || (r2 < 4. && omega < 0.14 && params->eobParams->omegaPeaked == 1  ) )
  {
    if(debugPK){
      XLAL_PRINT_INFO("\n Integration stopping for omega below threshold, omega=%f at r = %f\n", omega, sqrt(r2));
      fflush(NULL);
    }
    return 1;
  }

    if(r2 < 16. && omega > 1. )
    {
        if(debugPK){
            XLAL_PRINT_INFO("\n Integration stopping, omega>1 at r = %f\n", sqrt(r2));
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
      XLAL_PRINT_INFO("\n Integration stopping, dpdt > 10 -- too large!\n");
      fflush(NULL);}
    return 1;
  }

  /* If p_\Phi is too large numerically, break */
  if( r2 < 16. && values[5] > 10 )
  {
    if(debugPK){
      XLAL_PRINT_INFO("Integration stopping, Pphi > 10 now\n\n");
      fflush(NULL);
    }
    return 1;
  }
  /* If rdot inclreases, break */
   if ( r2 < 9. && rdot>params->prev_dr ) {
      if(debugPK){
          XLAL_PRINT_INFO("\n Integration stopping, dr/dt increasing!\n");
          fflush(NULL);
        }
      params->prev_dr=rdot;
      return 1;
    }
    params->prev_dr=rdot;

  /* **************************************************************** */
  /*              Last resort conditions                              */
  /* **************************************************************** */

  /* Very verbose output */
  if(debugPKverbose && r2 < 16.) {
    XLAL_PRINT_INFO("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
      t, values[0], values[1], values[2], values[3], values[4], values[5],
      values[6], values[7], values[8], values[9], values[10], values[11],
      values[12], values[13], omega);
  }

  return GSL_SUCCESS;
}


/**
 * Stopping condition for the regular resolution SEOBNRv1/2 orbital evolution
 * -- stop when reaching max orbital frequency in strong field.
 * At each test,
 * if omega starts to decrease, return 1 to stop evolution;
 * if not, update omega with current value and return GSL_SUCCESS to continue evolution.
 */
static int
XLALEOBSpinPrecAlignedStopCondition(double UNUSED t,  /**< UNUSED */
                           const double values[], /**< dynamical variable values */
                           double dvalues[],      /**< dynamical variable time derivative values */
                           void *funcParams       /**< physical parameters */
                          )
{
  int debugPK = 0;
  REAL8 omega, r;
  SpinEOBParams *params = (SpinEOBParams *)funcParams;

  r     = values[0];
  omega = dvalues[1];
    if (debugPK){  XLAL_PRINT_INFO("XLALEOBSpinPrecAlignedStopCondition:: r = %e\n", r);}

  if ( r < 6. && omega < params->eobParams->omega )
  {
    return 1;
  }

  params->eobParams->omega = omega;
  return GSL_SUCCESS;
}

/**
 * Stopping condition for the high resolution SEOBNRv1/2 orbital evolution
 * -- stop when reaching a minimum radius 0.3M out of the EOB horizon (Eqs. 9b, 37)
 * or when getting nan in any of the four ODE equations
 * At each test,
 * if conditions met, return 1 to stop evolution;
 * if not, return GSL_SUCCESS to continue evolution.
 */
static int
XLALSpinPrecAlignedHiSRStopCondition(double UNUSED t,  /**< UNUSED */
                           const double UNUSED values[], /**< dynamical variable values */
                           double dvalues[],      /**< dynamical variable time derivative values */
                           void UNUSED *funcParams       /**< physical parameters */
                          )
{
    int debugPK = 0;
    if (debugPK){
        XLAL_PRINT_INFO("XLALSpinPrecAlignedHiSRStopCondition:: r = %e\n", values[0]);
        XLAL_PRINT_INFO("values[0], values[1], values[2], values[3], dvalues[0], dvalues[1], dvalues[2], dvalues[3] = %e %e %e %e %e %e %e %e\n", values[0], values[1], values[2], values[3], dvalues[0], dvalues[1], dvalues[2], dvalues[3]);
    }

  if ( dvalues[2] >= 0. || isnan( dvalues[3] ) || isnan (dvalues[2]) || isnan (dvalues[1]) || isnan (dvalues[0]) )
  {
      return 1;
  }
  return GSL_SUCCESS;
}


/**
 * Standard interface for SEOBNRv3 waveform generator: calls XLALSimIMRSpinEOBWaveformAll
 */
int XLALSimIMRSpinEOBWaveform(
        REAL8TimeSeries **hplus, /**<< OUTPUT, +-polarization waveform */
         REAL8TimeSeries **hcross, /**<< OUTPUT, x-polarization waveform */
        const REAL8      phiC, /**<< coalescence orbital phase (rad) */
        const REAL8     deltaT, /**<< sampling time step */
        const REAL8     m1SI, /**<< mass-1 in SI unit (kg) */
        const REAL8     m2SI, /**<< mass-2 in SI unit (kg) 8*/
        const REAL8     fMin, /**<< starting frequency (Hz) */
        const REAL8     r, /**<< luminosity distance in SI unit (m) */
        const REAL8     inc, /**<< inclination angle */
        const REAL8     INspin1[],  /**<< spin1 */
        const REAL8     INspin2[] /**<< spin2 */
     )

{
    int ret;
    REAL8Vector   *dynamicsHi = NULL;
    REAL8Vector  *AttachPars = NULL;
    /** This time series contains harmonics in precessing (P) frame, no RD, for the end of the signal (high samling part)*/
    SphHarmTimeSeries *hlmPTSHi = NULL;
    /** This stores harmonics in J-frame, no RD, for the end of the signal (high samling part) */ 
    SphHarmTimeSeries *hlmPTSout = NULL;
    /** Stores harmonics in J-frame with RD,  for the end of the signal (high samling part)  */
    SphHarmTimeSeries *hIMRlmJTSHi = NULL;
    /** stores harmonics of the full waveform in I-frame */
    SphHarmTimeSeries *hIMR = NULL;

    ret = XLALSimIMRSpinEOBWaveformAll(hplus, hcross, &dynamicsHi, &hlmPTSout, &hlmPTSHi, &hIMRlmJTSHi, &hIMR, &AttachPars,
                        phiC, deltaT, m1SI, m2SI, fMin, r, inc, INspin1[0], INspin1[1], INspin1[2], INspin2[0], INspin2[1], INspin2[2]);
    if (ret == XLAL_SUCCESS){
        if (*hplus == NULL || *hcross == NULL){
             XLALPrintError("Houston-2, we've got a problem SOS, SOS, SOS, the waveform generator returns NULL!!!... m1 = %.18e, m2 = %.18e, fMin = %.18e, inclination = %.18e,   spin1 = {%.18e, %.18e, %.18e},   spin2 = {%.18e, %.18e, %.18e} \n", 
                       m1SI/LAL_MSUN_SI, m2SI/LAL_MSUN_SI, (double)fMin, (double)inc,  INspin1[0], INspin1[1], INspin1[2], INspin2[0], INspin2[1], INspin2[2]);
             XLAL_ERROR( XLAL_ENOMEM );
        }
        if ((*hplus)->data == NULL || (*hcross)->data == NULL){
             XLALPrintError("Houston-3, we've got a problem SOS, SOS, SOS, the waveform generator returns NULL!!!... m1 = %.18e, m2 = %.18e, fMin = %.18e, inclination = %.18e,   spin1 = {%.18e, %.18e, %.18e},   spin2 = {%.18e, %.18e, %.18e} \n", 
                       m1SI/LAL_MSUN_SI, m2SI/LAL_MSUN_SI, (double)fMin, (double)inc,  INspin1[0], INspin1[1], INspin1[2], INspin2[0], INspin2[1], INspin2[2]);
             XLAL_ERROR( XLAL_ENOMEM );
        }
    }

    if(dynamicsHi)
          XLALDestroyREAL8Vector( dynamicsHi );
    if(AttachPars)
          XLALDestroyREAL8Vector( AttachPars );
    if(hlmPTSout)
          XLALDestroySphHarmTimeSeries(hlmPTSout);
    if(hlmPTSHi)
          XLALDestroySphHarmTimeSeries(hlmPTSHi);
    if(hIMRlmJTSHi)
          XLALDestroySphHarmTimeSeries(hIMRlmJTSHi);
    if(hIMR)
          XLALDestroySphHarmTimeSeries(hIMR);
    return ret;
}

/**
 * This function generates precessing spinning SEOBNRv3 waveforms h+ and hx.
 * Currently, only h2m harmonics will be generated.
 *
 * Input conventions:
 * Cartesian coordinate system: initial \f$\vec{L}_N\f$ is in the xz plane, rotated away from
 * the z-axis by an angle inc
 * phiC       : in radians
 * deltaT     : in SI units (s)
 * m1SI, m2SI : in SI units (kg)
 * fMin       : in SI units (Hz)
 * r          : in SI units (m)
 * inc        : in radians
 * INspin{1,2}: in dimensionless units of m{1,2}^2
 *
 * Evolution conventions:
 * values[0-2]: r vector in units of total mass
 * values[3-5]: pstar vector in units of reduced mass
 * values[6-8]: S1 vector in units of (total mass)^2
 * values[9-11]: S2 vector in units of (total mass)^2
 * values[12-13]: phases (carrier and precession (Thomas)) in rads
 *
 * Note that when the initial opening angles of the spins w.r.t. the initial Newtonian
 * angular momentum are small, the aligned-spin SEOBNRv2 dynamics is used.
 * However, the waveforms are then generated according the the SEOBNRv3 model:
 * for example, they include the (2,1) mode.
 *
 * STEP 0) Prepare parameters, including pre-computed coefficients
 * for EOB Hamiltonian, flux and waveform
 * STEP 1) Solve for initial conditions
 * STEP 2) Evolve EOB trajectory both at low and high sampling rate
 * STEP 3) Compute Euler angles to go from initial inertial frame to
 * precessing frame
 * STEP 4) Locate merger point and at that time calculate J, chi and kappa,
 * and construct final J frame
 * STEP 5) Generate quasi-nonprecessing waveforms in precessing frame
 * STEP 6) Rotate quasi-nonprecessing waveforms from precessing to final-J-frame
 * STEP 7) Attach ringdown to final-J-frame modes
 * STEP 8) Rotate modes from final final-J-frame to initial inertial frame
 * STEP 9) Compute h+, hx
 */
int XLALSimIMRSpinEOBWaveformAll(
        REAL8TimeSeries **hplus,  /**<< output: hplus GW polarization */
        REAL8TimeSeries **hcross, /**<< output: hcross GW polarization */
        REAL8Vector      **dynHi, /**<< Here we store and return the seob dynamics for high sampling (end of inspiral) */
        SphHarmTimeSeries **hlmPTSoutput, /**<< Here we store and return the PWave (high sampling) */
        SphHarmTimeSeries **hlmPTSHiOutput, /**<< Here we store and return the JWave (high sampling) */
        SphHarmTimeSeries **hIMRlmJTSHiOutput, /**<< Here we store and return the JWaveIMR (high sampling) */
        SphHarmTimeSeries **hIMRoutput,       /**<< Here we store and return the IWave (full) */
        REAL8Vector     **AttachPars,   /**<< Parameters of RD attachment: */
        const REAL8      phiC,      /**<< intitial orbital phase */
        const REAL8     deltaT,     /**<< sampling time step */
        const REAL8     m1SI,       /**<< mass of first object in SI */
        const REAL8     m2SI,       /**<< mass of second object in SI */
        const REAL8     fMin,       /**<< fMin */
        const REAL8     r,          /**<< luminosity distance in SI */
        const REAL8     inc,        /**<< inclination */
        const REAL8     INspin1x,   /**<< spin1 x-component */
        const REAL8     INspin1y,   /**<< spin1 y-component */
        const REAL8     INspin1z,   /**<< spin1 z-component */
        const REAL8     INspin2x,   /**<< spin2 x-component */
        const REAL8     INspin2y,   /**<< spin2 y-component */
        const REAL8     INspin2z    /**<< spin2 z-component */
     )

{
  REAL8 INspin1[3], INspin2[3];
  INspin1[0] = INspin1x;
  INspin1[1] = INspin1y;
  INspin1[2] = INspin1z;
  INspin2[0] = INspin2x;
  INspin2[1] = INspin2y;
  INspin2[2] = INspin2z;

  INT4 debugPK = 0, debugCustomIC = 0, debugNoNQC = 0;
  INT4 debugRD = 0;
  FILE *out = NULL;
  INT4 i=0;
  INT4 k=0;
  UINT4 j=0;
  LIGOTimeGPS tc = LIGOTIMEGPSZERO;
  //REAL8 coa_phase_offset = 0;

  REAL8Vector *AttachParams = NULL;

  Approximant spinEOBApproximant = SEOBNRv3;
  /* The underlying aligned-spin EOB model is hard-coded here */
  INT4 SpinAlignedEOBversion = 2;

  /* Vector to store the initial spins */
  //REAL8 spin1[3] = {0,0,0}, spin2[3] = {0,0,0}, InitLhat[3] = {sin(inc),0.,cos(inc)};
  REAL8 spin1[3] = {0,0,0}, spin2[3] = {0,0,0}, InitLhat[3] = {0.0,0.0,1.0};
  // FIXME
  memcpy( spin1, INspin1, 3*sizeof(REAL8));
  memcpy( spin2, INspin2, 3*sizeof(REAL8));

  /* Check the initial misalignment of component spins, w.r.t. the z-axis
   * theta{1,2}Ini measure the angle w.r.t. the orbital ang. momentum
   * */
  REAL8 theta1Ini = 0, theta2Ini = 0;
  REAL8 spin1Norm = -1, spin2Norm = -1;
  spin1Norm = sqrt( INspin1[0]*INspin1[0] + INspin1[1]*INspin1[1] + INspin1[2]*INspin1[2] );
  spin2Norm = sqrt( INspin2[0]*INspin2[0] + INspin2[1]*INspin2[1] + INspin2[2]*INspin2[2] );
  if ( INspin1[0] == 0. && INspin1[1] == 0. && INspin1[2] == 0. ) {
        spin1Norm = 0.;
  }
  if ( INspin2[0] == 0. && INspin2[1] == 0. && INspin2[2] == 0. ) {
        spin2Norm = 0.;
  }

  /** NOTE: if the spin amgnitude is less than 1.e-5 we put them explicitely to zero! */ 
  if (spin1Norm <= 1.0e-5) {INspin1[0]=0.;INspin1[1]=0.;INspin1[2]=0.;}
  if (spin2Norm <= 1.0e-5) {INspin2[0]=0.;INspin2[1]=0.;INspin2[2]=0.;}
  REAL8 acosarg1 = 0, acosarg2 = 0;
  if( spin1Norm > 1.0e-5){
      acosarg1 = (InitLhat[0]*INspin1[0] + InitLhat[1]*INspin1[1] + InitLhat[2]*INspin1[2])/ spin1Norm ;
      if (acosarg1 > 1.) {
          theta1Ini = 0.;
      }
      else if (acosarg1 < -1.) {
          theta1Ini = LAL_PI;
      }
      else {
          theta1Ini = acos( acosarg1 );
      }
    }

  if( spin2Norm > 1.0e-5){
      acosarg2 = (InitLhat[0]*INspin2[0] + InitLhat[1]*INspin2[1] + InitLhat[2]*INspin2[2])/ spin2Norm ;
      if (acosarg2 > 1.) {
          theta2Ini = 0.;
      }
      else if (acosarg2 < -1.) {
          theta2Ini = LAL_PI;
      }
      else {
          theta2Ini = acos( acosarg2 );
      }
  }

    /** If the misalignment angle of spins with orbital angular momentum is less than <1.e-4 we treat them as aligned for evolution of the 
     * orbital dynamics  */
    REAL8 EPS_ALIGN = 1.0e-4, incA = inc;
    bool SpinsAlmostAligned = ( fabs(theta1Ini) <= EPS_ALIGN  || fabs(theta1Ini) >= LAL_PI - EPS_ALIGN) && ( fabs(theta2Ini) <= EPS_ALIGN || fabs(theta2Ini) >= LAL_PI - EPS_ALIGN);
    if ( SpinsAlmostAligned ) {
        if (debugPK) {XLAL_PRINT_INFO("Both spins are almost aligned/antialigned to initial LNhat: we use V2 dynamics\n");
            XLAL_PRINT_INFO("spin1Norm spin2Norm %e %e\n", spin1Norm, spin2Norm);
        }
        spin1[0] = 0.;
        spin1[1] = 0.;
        spin1[2] = spin1Norm*copysign(1., cos(theta1Ini));
        spin2[0] = 0.;
        spin2[1] = 0.;
        spin2[2] = spin2Norm*copysign(1., cos(theta2Ini));
        incA = 0.;
    }

    if ( debugPK ) {
        XLAL_PRINT_INFO( "InitLhat = {%3.10f, %3.10f, %3.10f}\n",
               InitLhat[0], InitLhat[1], InitLhat[2] );

    XLAL_PRINT_INFO( "theta1Ini, theta2Ini =  %3.10f, %3.10f\n", theta1Ini, theta2Ini );
    XLAL_PRINT_INFO( "INspin1 = {%3.10f, %3.10f, %3.10f}\n",
            INspin1[0], INspin1[1], INspin1[2] );
    XLAL_PRINT_INFO( "INspin2 = {%3.10f, %3.10f, %3.10f}\n",
            INspin2[0], INspin2[1], INspin2[2] );
    XLAL_PRINT_INFO( "spin1 = {%3.10f, %3.10f, %3.10f}\n", spin1[0], spin1[1], spin1[2] );
    XLAL_PRINT_INFO( "spin2 = {%3.10f, %3.10f, %3.10f}\n", spin2[0], spin2[1], spin2[2] );
  }

  /* *******************************************************************/
  /* ********************** Memory Allocation **************************/
  /* *******************************************************************/

  /* Allocate the values vector to contain the ICs            */
  /* Throughout the code, there are 12 dynamical variables:   */
  /* values[0-2]  is x (Cartesian tortoise separation vector)  */
  /* values[3-5]  is p (Cartesian momentum)                    */
  /* values[6-8]  is spin of body 1                            */
  /* values[9-11] is spin of body 2                            */
  /* Two additional orbital quantities are evolved            */
  /* values[12]   is orbital phase                             */
  /* values[13]   is alpha dot cos(inclination)                */
    REAL8Vector *values = NULL;
  if ( !(values = XLALCreateREAL8Vector( 14 )) )
  {
    XLALPrintError("Failed to allocate REAL8Vector values!\n");
    PRINT_PARAMS
    XLAL_ERROR(  XLAL_ENOMEM );
  }
  memset( values->data, 0, values->length * sizeof( REAL8 ));

  /* Vector to contain derivatives of ICs */
  REAL8Vector *dvalues = NULL;
  if ( !(dvalues = XLALCreateREAL8Vector( 14 )) )
  {
    XLALDestroyREAL8Vector( values );
    XLALPrintError("Failed to allocate REAL8Vector dvalues!\n");
    PRINT_PARAMS
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
    XLALDestroyREAL8Vector( dvalues );
    XLALPrintError("Failed to allocate REAL8Vector sigmaStar!\n");
    PRINT_PARAMS
    XLAL_ERROR( XLAL_ENOMEM );
  }
  if ( !(sigmaKerr = XLALCreateREAL8Vector( 3 )) )
  {
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLALDestroyREAL8Vector( dvalues );
    XLALPrintError("Failed to allocate REAL8Vector sigmaKerr!\n");
    PRINT_PARAMS
    XLAL_ERROR( XLAL_ENOMEM );
  }

  /* Spins not scaled by the mass */
  REAL8 mSpin1[3] = {0,0,0}, mSpin2[3] = {0,0,0};

  /* Wrapper spin vectors used to calculate sigmas */
  REAL8Vector s1Vec, s1VecOverMtMt;
  REAL8Vector s2Vec, s2VecOverMtMt;
   /* In s1Data and s2Data will store chi_i * m_i^2 in geomtric units  */
   /* In s1DataNorm and s2DataNorm will store chi_i * m_i^2/M^2 */
  REAL8       s1Data[3]     = {0,0,0}, s2Data[3]     = {0,0,0},
              s1DataNorm[3] = {0,0,0}, s2DataNorm[3] = {0,0,0};

  /* Parameters of the system */
  REAL8 m1 = 0, m2 = 0, mTotal = 0, eta = 0, mTScaled = 0;
  REAL8 amp0 = 0;

  /* Dynamics of the system */
  REAL8Vector tVec, phiDMod, phiMod;
  REAL8Vector posVecx, posVecy, posVecz, momVecx, momVecy, momVecz;
  REAL8Vector s1Vecx, s1Vecy, s1Vecz, s2Vecx, s2Vecy, s2Vecz;
  REAL8       omega = 0, v = 0, ham = 0;

  /* Cartesian vectors needed to calculate Hamiltonian */
  REAL8Vector cartPosVec, cartMomVec;
  REAL8 cartPosData[3] = {0,0,0}, cartMomData[3] = {0,0,0};
  REAL8 rcrossrdotNorm = 0, rvec[3]    = {0,0,0}, rcrossrdot[3] = {0,0,0};
  REAL8 pvec[3] = {0,0,0},  rcrossp[3] = {0,0,0};//, rcrosspMag    = 0;
  REAL8 s1dotLN = 0, s2dotLN = 0;

  /* Polar vectors needed for waveform modes calculation */
  REAL8Vector  polarDynamics;
  REAL8 polData[4] = {0,0,0,0};

  /* Signal mode */
  COMPLEX16   hLM = 0.0 + 0.0*j;

  /* Non-quasicircular correction */
  COMPLEX16      hNQC = 0.0 + 0.0*j;
  EOBNonQCCoeffs nqcCoeffs;
  memset( &nqcCoeffs, 0, sizeof( nqcCoeffs ) );

  /* Ringdown freq used to check the sample rate */
  COMPLEX16Vector  modefreqVec;
  COMPLEX16        modeFreq;

  /* Spin-weighted spherical harmonics */
  COMPLEX16  Y22 = 0.0 + 0.0j, Y2m2 = 0.0 + 0.0j, Y21 = 0.0 + 0.0j, Y2m1 = 0.0 + 0.0j;
  COMPLEX16  Y20 = 0.0 + 0.0j;
  /* (2,2) and (2,-2) spherical harmonics needed in (h+,hx) */

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


  /* Set up structures and calculate necessary PN parameters */
    /* SpinEOBParams contains:
        EOBParams               *eobParams; // see below
        SpinEOBHCoeffs       *seobCoeffs; // see below
        EOBNonQCCoeffs     *nqcCoeffs; // see below
        REAL8Vector             *s1Vec; // spin1
        REAL8Vector             *s2Vec; // spin2
        REAL8Vector             *sigmaStar; // Eq. 5.3 of Barausse and Buonanno PRD 81, 084024 (2010)
        REAL8Vector             *sigmaKerr; // Eq. 5.2 of Barausse and Buonanno PRD 81, 084024 (2010)
        REAL8                        a; // |sigmaKerr|
        REAL8                        chi1; // spin1.LNhat/m1^2
        REAL8                        chi2; // spin2.LNhat/m2^2
        REAL8                        prev_dr; // stores dr/dt for stopping condition purposes
        int                              alignedSpins; // flag to indicate whther the binary is precessing or not
        int                              tortoise; // flag to switch on/off tortoise coordinates when calling Hamiltonian
        int i                            gnoreflux;  // flag to swith off radiation reaction when calling the ODEs via XLALSpinPrecHcapNumericalDerivative
        */
    SpinEOBParams           seobParams;
    /* SpinEOBHCoeffs contains:
        double KK; // nonspinning calibration in Hamiltonian (for SEOBNRv2: 1.712 − 1.804eta − 39:77eta^2 + 103.2eta^3)
        double k0; // Delta_i coefficients in the Delta_u potential Eq. 8 of PRD 86, 024011 (2012) and https://dcc.ligo.org/T1400476
        double k1;
        double k2;
        double k3;
        double k4;
        double k5;
        double k5l;
        double b3; // omega_{fd}^1 frame dragging parameter in Hamiltonian (unused)
        double bb3; // omega_{fd}^2 frame dragging parameter in Hamiltonian (unused)
        double d1; // spin-orbit calibration in Hamiltonian for SEOBNRv1
        double d1v2; // spin-orbit calibration in Hamiltonian for SEOBNRv1
        double dheffSS; // spin-spin calibration in Hamiltonian for SEOBNRv1
        double dheffSSv2; // spin-spin calibration in Hamiltonian for SEOBNRv1
        UINT4    SpinAlignedEOBversion;
        int      updateHCoeffs;
    */
    SpinEOBHCoeffs          seobCoeffs;
    /* EOBParams contains parameters common to nonspin and spin EOBNR models,
        including mass ratio, masses, pre-computed coefficients for potential, flux and waveforms,
        NQC coefficients and Newtonian multiple prefixes */
    EOBParams               eobParams;
    /*FacWaveformCoeffscontaining the coefficients for calculating the factorized waveform.
      The coefficients are precomputed in the function  XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients */
    FacWaveformCoeffs       hCoeffs;
    /* NewtonMultipolePrefixes contains all the terms of the Newtonian multipole which
       are constant over the course of the evolution, and can therefore be
       pre-computed. They are stored in a two-dimensional array, which is
       indexed as values[l][m]. This is filled by XLALSimIMREOBComputeNewtonMultipolePrefixes */
    NewtonMultipolePrefixes prefixes;

  memset( &seobParams, 0, sizeof(seobParams) );
  memset( &seobCoeffs, 0, sizeof(seobCoeffs) );
  memset( &eobParams,  0, sizeof(eobParams) );
  memset( &hCoeffs,    0, sizeof( hCoeffs ) );
  memset( &prefixes,   0, sizeof( prefixes ) );
    
    REAL8TimeSeries *alphaI2PTS = NULL, *betaI2PTS = NULL, *gammaI2PTS = NULL, *alphaP2JTS = NULL, *betaP2JTS = NULL, *gammaP2JTS = NULL;
    REAL8TimeSeries *alphaI2PTSHi = NULL, *betaI2PTSHi = NULL, *gammaI2PTSHi = NULL, *alphaP2JTSHi = NULL, *betaP2JTSHi = NULL, *gammaP2JTSHi = NULL;
    COMPLEX16TimeSeries *h22TS = NULL, *h21TS = NULL, *h20TS = NULL, *h2m1TS = NULL, *h2m2TS = NULL;
    COMPLEX16TimeSeries *h22TSHi = NULL, *h21TSHi = NULL, *h20TSHi = NULL, *h2m1TSHi = NULL, *h2m2TSHi = NULL;
    
    COMPLEX16TimeSeries *hIMRJTSHi  = NULL;
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
    
    COMPLEX16TimeSeries *hIMRJTS2mHi    = NULL;
    COMPLEX16TimeSeries *hIMR22ITS  = NULL;
    COMPLEX16TimeSeries *hIMR21ITS  = NULL;
    COMPLEX16TimeSeries *hIMR20ITS  = NULL;
    COMPLEX16TimeSeries *hIMR2m1ITS = NULL;
    COMPLEX16TimeSeries *hIMR2m2ITS = NULL;
    REAL8TimeSeries *alpI = NULL, *betI = NULL, *gamI = NULL;
    
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
    REAL8TimeSeries  *hPlusTS  = NULL;
    REAL8TimeSeries  *hCrossTS = NULL;
    COMPLEX16TimeSeries *hIMRJTS = NULL;
    
    
  /* Miscellaneous memory for Ringdown attachment        */
  REAL8 tPeakOmega = 0, tAttach = 0, combSize = 0,/*longCombSize,*/ deltaNQC =0;
  REAL8  sh  = 0;
  REAL8 magR = 0, Lx   = 0, Ly   = 0, Lz = 0, magL = 0, magJ = 0,
        LNhx = 0, LNhy = 0, LNhz = 0, magLN =0, Jx = 0, Jy   = 0, Jz = 0;
  REAL8 aI2P = 0, bI2P = 0, gI2P = 0;
  REAL8 aP2J = 0, bP2J = 0, gP2J = 0;
  REAL8 chi1J= 0, chi2J= 0, chiJ = 0;
  REAL8 chi1L= 0.0, chi2L = 0.0, chiL = 0.0;
  REAL8 kappaJL = 0;
  REAL8 JLN = 0.0;
  REAL8 JframeEx[3] = {0,0,0}, JframeEy[3] = {0,0,0}, JframeEz[3] = {0,0,0};

  /* Variables for the integrator */
  LALAdaptiveRungeKutta4Integrator       *integrator = NULL;
  REAL8Array              *dynamics   = NULL;
  REAL8Array              *dynamicsHi   = NULL;
  REAL8Array              *dynamicsV2   = NULL;
  REAL8Array              *dynamicsV2Hi   = NULL;
  INT4                    retLen = 0, retLenLow = 0, retLenHi = 0;
  INT4                    retLenRDPatchHi= 0, retLenRDPatchLow = 0;

  /* Interpolation spline */
  gsl_spline    *spline = NULL;
  gsl_interp_accel *acc = NULL;

  /* Accuracies of adaptive Runge-Kutta integrator */
    /* Note that this accuracies are lower than those used in SEOBNRv2: they allow reasonable runtimes for precessing systems */
    /* These accuracies can be adjusted according to desired accuracy and runtime */
   REAL8 EPS_ABS = 1.0e-8;
  const REAL8 EPS_REL = 1.0e-8;
  /* Relax abs accuracy in case of highly symmetric case that would otherwise slow down significantly */
  if (sqrt((INspin1[0] + INspin2[0])*(INspin1[0] + INspin2[0]) + (INspin1[1] + INspin2[1])*(INspin1[1] + INspin2[1])) < 1.0e-10 && !SpinsAlmostAligned)
  {
      if (debugPK) XLAL_PRINT_INFO("EPS_ABS is decreased!\n");
      EPS_ABS = 1.0e-4;
  }

    /* Memory for the calculation of the alpha(t) and beta(t) angles */
  REAL8Vector *Alpha = NULL, *Beta = NULL, *Gamma = NULL;
  REAL8Vector *AlphaHi = NULL,  *BetaHi = NULL, *GammaHi = NULL;
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
  REAL8Sequence *timeJFull = NULL;
  REAL8Sequence *timeIFull = NULL;

  /* Memory for ringdown attachment */
  REAL8Vector *rdMatchPoint = XLALCreateREAL8Vector( 3 );
  if ( !rdMatchPoint )
  {
    XLALPrintError("Failed to allocate REAL8Vector rdMatchPoint!\n");
    PRINT_PARAMS
    XLAL_ERROR( XLAL_ENOMEM );
  }
  REAL8 alJtoI = 0, betJtoI = 0, gamJtoI = 0;
  SphHarmTimeSeries *hIMRlmITS = NULL;
  COMPLEX16 x11 = 0.0 + 0.0j;

  /* *******************************************************************/
  /* ********************** Memory Initialization **********************/
  /* *******************************************************************/

  /* Initialize mass parameters */
  m1       = m1SI / LAL_MSUN_SI;
  m2       = m2SI / LAL_MSUN_SI;
  mTotal   = m1 + m2;
  mTScaled = mTotal * LAL_MTSUN_SI;
  eta      = m1 * m2 / (mTotal*mTotal);

  if (eta > 0.25){
      m2=m1;
      eta = 0.25;
  }

  /* Initialize amplitude scaling parameter */
  amp0 = mTotal * LAL_MRSUN_SI / r;

  if (debugPK) {
    XLAL_PRINT_INFO("Stas, here is the passes functions\n");
   XLAL_PRINT_INFO("Inputs: m1 = %.16e, m2 = %.16e, fMin = %.16e, inclination = %.16e\n",
          m1, m2, (double) fMin, (double) inc );
    XLAL_PRINT_INFO("Mtotal = %.16e, eta = %.16e \n", mTotal, eta);
    XLAL_PRINT_INFO("Inputs: spin1 = {%.16e, %.16e, %.16e}\n",
            spin1[0], spin1[1], spin1[2]);
    XLAL_PRINT_INFO("Inputs: spin2 = {%.16e, %.16e, %.16e}\n",
            spin2[0], spin2[1], spin2[2]);
  }

  /* Calculate the time we will need to step back for ringdown */
    /* Stepping back 150M in time is very conservative: the attachment point usually lies in the last 20M of the trajectory */
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
    s1Vec.data[i]     = mSpin1[i] = spin1[i] * m1 * m1;
    s2Vec.data[i]     = mSpin2[i] = spin2[i] * m2 * m2;
    s1VecOverMtMt.data[i] = s1Vec.data[i] / mTotal / mTotal;
    s2VecOverMtMt.data[i] = s2Vec.data[i] / mTotal / mTotal;
  }

  if (debugPK){
    XLAL_PRINT_INFO("Stas, here is the passes functions\n");
   XLAL_PRINT_INFO("Inputs: m1 = %.16e, m2 = %.16e, fMin = %.16e, inclination = %.16e\n",
            m1, m2, (double) fMin, (double) inc );
    XLAL_PRINT_INFO("Inputs: spin1 = {%.16e, %.16e, %.16e}\n",
            spin1[0], spin1[1], spin1[2]);
    XLAL_PRINT_INFO("Inputs: spin2 = {%.16e, %.16e, %.16e}\n",
            spin2[0], spin2[1], spin2[2]);
    XLAL_PRINT_INFO("Inputs: s1V = {%.16e, %.16e, %.16e}\n",
            s1Vec.data[0], s1Vec.data[1], s1Vec.data[2]);
    XLAL_PRINT_INFO("Inputs: s2V = {%.16e, %.16e, %.16e}\n",
            s2Vec.data[0], s2Vec.data[1], s2Vec.data[2]);
  }

/* *********************************************************************************
 * *********************************************************************************
 * STEP 0) Prepare parameters, including pre-computed coefficients
 * for EOB Hamiltonian, flux and waveform
 * *********************************************************************************
 * ********************************************************************************* */
  /* ************************************************* */
  /* Populate the initial structures                   */
  /* ************************************************* */
  /* Spin parameters */
  if ( XLALSimIMRSpinEOBCalculateSigmaStar(
          sigmaStar, m1, m2, &s1Vec, &s2Vec ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLALDestroyREAL8Vector( dvalues );
    XLALDestroyREAL8Vector( rdMatchPoint );
    XLALPrintError("XLALSimIMRSpinEOBCalculateSigmaStar failed!\n");
    PRINT_PARAMS
    XLAL_ERROR( XLAL_EDOM );
  }
  if ( XLALSimIMRSpinEOBCalculateSigmaKerr(
          sigmaKerr, m1, m2, &s1Vec, &s2Vec ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLALDestroyREAL8Vector( dvalues );
    XLALDestroyREAL8Vector( rdMatchPoint );
    XLALPrintError("XLALSimIMRSpinEOBCalculateSigmaKerr failed!\n");
    PRINT_PARAMS
    XLAL_ERROR( XLAL_EDOM );
  }

  /* Calculate the value of a, that is magnitude of Eq. 31 in PRD 86, 024011 (2012) */
  seobParams.a = a = sqrt( inner_product(sigmaKerr->data, sigmaKerr->data) );
  seobParams.s1Vec = &s1VecOverMtMt;
  seobParams.s2Vec = &s2VecOverMtMt;

  if (debugPK){
      XLAL_PRINT_INFO("Inputs: sigma = {%.16e, %.16e, %.16e}\n",
        sigmaKerr->data[0], sigmaKerr->data[1], sigmaKerr->data[2]);
      XLAL_PRINT_INFO("Inputs: star = {%.16e, %.16e, %.16e}\n",
        sigmaStar->data[0], sigmaStar->data[1], sigmaStar->data[2]);
      XLAL_PRINT_INFO("Inputs: a = %.16e\n", a);
      fflush(NULL);
  }

  if ( !(AttachParams = XLALCreateREAL8Vector( 5 )) )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLALDestroyREAL8Vector( dvalues );
    XLALDestroyREAL8Vector( rdMatchPoint );
    XLALPrintError("Failed to allocate REAL8Vector AttachParams!\n");
    PRINT_PARAMS
    XLAL_ERROR( XLAL_ENOMEM );
  }
  *AttachPars = AttachParams;

   /* Cartesian vectors needed to calculate Hamiltonian */
  cartPosVec.length = cartMomVec.length = 3;
  cartPosVec.data = cartPosData;
  cartMomVec.data = cartMomData;
  memset( cartPosData, 0, sizeof( cartPosData ) );
  memset( cartMomData, 0, sizeof( cartMomData ) );


  /* ************************************************* */
  /* Waveform parameter structures              */
  /* ************************************************* */
  /* Spin-EOB parameters */
  seobParams.alignedSpins = 0;
  seobParams.tortoise     = 1;
  seobParams.ignoreflux = 0;
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

  /* Pre-compute the Hamiltonian coefficients
    double KK; // nonspinning calibration in Hamiltonian (for SEOBNRv2: 1.712 − 1.804eta − 39:77eta^2 + 103.2eta^3)
    double k0; // Delta_i coefficients in the Delta_u potential Eq. 8 of PRD 86, 024011 (2012) and https://dcc.ligo.org/T1400476
    double k1;
    double k2;
    double k3;
    double k4;
    double k5;
    double k5l;
    double b3; // omega_{fd}^1 frame dragging parameter in Hamiltonian (unused)
    double bb3; // omega_{fd}^2 frame dragging parameter in Hamiltonian (unused)
    double d1; // spin-orbit calibration in Hamiltonian for SEOBNRv1
    double d1v2; // spin-orbit calibration in Hamiltonian for SEOBNRv1
    double dheffSS; // spin-spin calibration in Hamiltonian for SEOBNRv1
    double dheffSSv2; // spin-spin calibration in Hamiltonian for SEOBNRv1
   */
  if ( XLALSimIMRCalculateSpinPrecEOBHCoeffs( &seobCoeffs, eta, a,
                          SpinAlignedEOBversion ) == XLAL_FAILURE ) /* This function returns XLAL_SUCCESS or calls XLAL_ERROR( XLAL_EINVAL ) */
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLALDestroyREAL8Vector( dvalues );
    XLALDestroyREAL8Vector( rdMatchPoint );
    XLALPrintError("XLALSimIMRCalculateSpinPrecEOBHCoeffs failed!\n");
    PRINT_PARAMS
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Pre-compute the coefficients for the Newtonian factor of hLM
     Eq. A1 of PRD 86, 024011 (2012)  */
  if ( XLALSimIMREOBComputeNewtonMultipolePrefixes( &prefixes, eobParams.m1,
			eobParams.m2 ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLALDestroyREAL8Vector( dvalues );
    XLALDestroyREAL8Vector( rdMatchPoint );
    XLALPrintError("XLALSimIMREOBComputeNewtonMultipolePrefixes failed!\n");
    PRINT_PARAMS
    XLAL_ERROR( XLAL_EDOM );
  }

/* *********************************************************************************
 * *********************************************************************************
 * STEP 1) Solve for initial conditions, according to Sec. IV A of
 * PRD 74, 104005 (2006)
 * *********************************************************************************
 * ********************************************************************************* */
  if( debugPK )
  {
    XLAL_PRINT_INFO("Calling the XLALSimIMRSpinEOBInitialConditionsPrec function!\n");
    XLAL_PRINT_INFO(
      "Inputs: m1 = %.16e, m2 = %.16e, fMin = %.16e, inclination = %.16e\n",
                      m1, m2, (double) fMin, (double) inc );
    XLAL_PRINT_INFO("Inputs: mSpin1 = {%.16e, %.16e, %.16e}\n",
                      mSpin1[0], mSpin1[1], mSpin1[2]);
    XLAL_PRINT_INFO("Inputs: mSpin2 = {%.16e, %.16e, %.16e}\n",
                      mSpin2[0], mSpin2[1], mSpin2[2]);
    fflush(NULL);
  }

  //values = XLALCreateREAL8Vector( 14 );
  REAL8 incl_temp = 0.0;  // !!!! For comparison with C++ and NR we need inc = 0 for initial conditions
  
  //incl_temp = inc; //FIXME
  /* The initial condition construction is based on PRD 74, 104005 (2006) */
  /* If the initial spin opening angles are small, then use SEOBNRv2 (aligned-spin) dyamics */
  if ( SpinsAlmostAligned ) {
        seobParams.alignedSpins = 1;
        seobParams.chi1 = copysign( spin1Norm, cos(theta1Ini) );
        seobParams.chi2 = copysign( spin2Norm, cos(theta2Ini) );
        chiS = 0.5*(seobParams.chi1 + seobParams.chi2);
        chiA = 0.5*(seobParams.chi1 - seobParams.chi2);
        tplspin = (1.-2.*eta) * chiS + (m1 - m2)/(m1 + m2) * chiA;
        if ( XLALSimIMREOBCalcSpinFacWaveformCoefficients( seobParams.eobParams->hCoeffs, m1, m2, eta, tplspin, chiS, chiA, SpinAlignedEOBversion ) == XLAL_FAILURE ) /* This function returns XLAL_SUCCESS or calls XLAL_ERROR( XLAL_EINVAL ) */
      {
            XLALDestroyREAL8Vector( sigmaKerr );
            XLALDestroyREAL8Vector( sigmaStar );
            XLALDestroyREAL8Vector( values );
            XLALDestroyREAL8Vector( dvalues );
            XLALDestroyREAL8Vector( rdMatchPoint );
            XLALPrintError("XLALSimIMREOBCalcSpinFacWaveformCoefficients failed!\n");
            PRINT_PARAMS
            XLAL_ERROR( XLAL_EFUNC );
        }

        if ( XLALSimIMRSpinEOBInitialConditions( values, m1, m2, fMin, incA,
                                                mSpin1, mSpin2, &seobParams, /* use_optimized_v2 = */ 0 ) == XLAL_FAILURE ) /* This function returns XLAL_SUCCESS or calls XLAL_ERROR with XLAL_EINVAL, XLAL_ENOMEM, or XLAL_EMAXITER */
        {
            XLALDestroyREAL8Vector( sigmaKerr );
            XLALDestroyREAL8Vector( sigmaStar );
            XLALDestroyREAL8Vector( values );
            XLALDestroyREAL8Vector( dvalues );
            XLALDestroyREAL8Vector( rdMatchPoint );
            XLALPrintError("XLALSimIMRSpinEOBInitialConditions failed!\n");
            PRINT_PARAMS
            XLAL_ERROR( XLAL_EFUNC );
        }
        seobParams.alignedSpins = 0;
    }
    else {
        if ( XLALSimIMRSpinEOBInitialConditionsPrec( values, m1, m2, fMin, incl_temp,
                                                    mSpin1, mSpin2, &seobParams ) == XLAL_FAILURE ) /* This function returns XLAL_SUCCESS or calls XLAL_ERROR with XLAL_EINVAL, XLAL_ENOMEM, or XLAL_EMAXITER */
        {
            XLALDestroyREAL8Vector( sigmaKerr );
            XLALDestroyREAL8Vector( sigmaStar );
            XLALDestroyREAL8Vector( values );
            XLALDestroyREAL8Vector( dvalues );
            XLALDestroyREAL8Vector( rdMatchPoint );
            XLALPrintError("XLALSimIMRSpinEOBInitialConditionsPrec failed!\n");
            PRINT_PARAMS
            XLAL_ERROR( XLAL_EFUNC );
        }
  }
  /* Initial phases */
  values->data[12] = 0.;
  values->data[13] = 0.;

  if(debugCustomIC) {
      /* Hardcode initial conditions for debugging purposes here */
      values->data[0]= 15.4898001256;
      values->data[1]= 0.;
      values->data[2]= -4.272074867808132e-19;
      values->data[3]= -0.0009339635043526475;
      values->data[4]= 0.2802596444562164;
      values->data[5]= -0.0001262371378125648;
      values->data[6]= 0.03950299224160406;
      values->data[7]= 0.1033495061350166;
      values->data[8]= 0.02382287037711037;
      values->data[9]= 0.07463668902857602;
      values->data[10]= 0.001769731591445356;
      values->data[11]= 0.04303525354405329;
  }

  if(debugPK)
  {
    XLAL_PRINT_INFO("Setting up initial conditions, returned values are:\n");
	  for( j=0; j < values->length; j++)
      XLAL_PRINT_INFO("%.16le\n", values->data[j]);
  }

  if(debugPK)XLAL_PRINT_INFO("\nReached the point where LN is to be calculated\n");

  memset( dvalues->data, 0, 14 * sizeof(REAL8) );
  if( XLALSpinPrecHcapRvecDerivative( 0, values->data, dvalues->data,
          (void*) &seobParams) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLALDestroyREAL8Vector( dvalues );
    XLALDestroyREAL8Vector( rdMatchPoint );
    XLALPrintError("XLALSpinPrecHcapRvecDerivative failed!\n");
    PRINT_PARAMS
    XLAL_ERROR( XLAL_EDOM );
  }
  if(debugPK)XLAL_PRINT_INFO("\nCalculated Rdot\n");

  /* Calculate r cross rDot */
  memcpy( rdotvec, dvalues->data, 3*sizeof(REAL8) );
  memcpy( rvec,    values->data,  3*sizeof(REAL8) );
  memcpy( pvec,    values->data+3,3*sizeof(REAL8) );

  cross_product( rvec, rdotvec, rcrossrdot );
  rcrossrdotNorm = sqrt(inner_product( rcrossrdot, rcrossrdot ));
  for( i = 0; i < 3; i++ ) { rcrossrdot[i] /= rcrossrdotNorm; }

/* Calculate the values of chiS and chiA, as given in Eq. 17 of
 * PRD 89, 084006 (2014) */
  s1dotLN = inner_product( spin1, rcrossrdot );
  s2dotLN = inner_product( spin2, rcrossrdot );

///* An alternative is to project the spins onto L = rXp */
// cross_product( rvec, pvec, rcrossp );
//  rcrosspMag = sqrt(inner_product(rcrossp, rcrossp));
//  for( i = 0; i < 3; i++ ) { rcrossp[i] /= rcrosspMag; }
//
//  s1dotL = inner_product( spin1, rcrossrdot );
//  s2dotL = inner_product( spin2, rcrossrdot );

  if(debugPK) {
    XLAL_PRINT_INFO("rXp = %3.10f %3.10f %3.10f\n", rcrossp[0], rcrossp[1], rcrossp[2]);
    XLAL_PRINT_INFO("rXrdot = %3.10f %3.10f %3.10f\n", rcrossrdot[0], rcrossrdot[1], rcrossrdot[2]);
    fflush(NULL);
  }

  chiS = 0.5 * (s1dotLN + s2dotLN);
  chiA = 0.5 * (s1dotLN - s2dotLN);

  /* Compute the test-particle limit spin of the deformed-Kerr background: (S1+S2).LNhat/Mtot^2 */
  switch ( SpinAlignedEOBversion )
  {
    case 1:
       /* See below Eq. 17 of PRD 86, 041011 (2012) */
      tplspin = 0.0;
      break;
    case 2:
      /* See below Eq. 4 of PRD 89, 061502(R) (2014)*/
      tplspin = (1.-2.*eta) * chiS + (m1 - m2)/(m1 + m2) * chiA;
      break;
    default:
        XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
        XLALDestroyREAL8Vector( sigmaKerr );
        XLALDestroyREAL8Vector( sigmaStar );
        XLALDestroyREAL8Vector( values );
        XLALDestroyREAL8Vector( dvalues );
        XLALDestroyREAL8Vector( rdMatchPoint );
        PRINT_PARAMS
        XLAL_ERROR( XLAL_EINVAL );
      break;
  }

  /* ************************************************* */
  /* Populate the Waveform initial structures  */
  /* ************************************************* */
  /* Pre-compute the non-spinning and spinning coefficients for hLM factors */
  if( debugPK ) {
    XLAL_PRINT_INFO("tplspin = %.12e, chiS = %.12e, chiA = %.12e\n", tplspin, chiS, chiA);
    fflush(NULL);
  }

  if ( XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients( eobParams.hCoeffs, m1, m2, eta,
        tplspin, chiS, chiA, 3 ) == XLAL_FAILURE ) /* This function returns XLAL_SUCCESS or calls XLAL_ERROR( XLAL_EINVAL ) */
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLALDestroyREAL8Vector( dvalues );
    XLALDestroyREAL8Vector( rdMatchPoint );
    XLALPrintError("XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients failed!\n");
    PRINT_PARAMS
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* ************************************************* */
  /* Compute the coefficients for the NQC Corrections  */
  /* ************************************************* */
    /* See below Eq. 4 of PRD 89, 061502(R) (2014)*/
  spinNQC = (1.-2.*eta) * chiS + (m1 - m2)/(m1 + m2) * chiA;

  /* Check if initial frequency is too high: we choose an initial minimum separation
   *  of 10M as a compromise between reliability of initial conditions and length
   *  of the waveform */
  REAL8 NRPeakOmega22 = GetNRSpinPeakOmegav2(2, 2, eta, spinNQC) / mTScaled;
  REAL8 freqMinRad = pow(10.0, -1.5)/(LAL_PI*mTScaled);
  REAL8 signOfa = copysign(1., spinNQC);
  REAL8 spn2 = spinNQC*spinNQC;
   /* Kerr ISCO radius from Eq. 2.21 of Astrophysical Journal, Vol. 178, pp. 347-370 (1972) */
  REAL8 Z1 = 1.0 + pow(1.0 - spn2, 1./3.)*(pow(1.0+spinNQC, 1./3.) + pow(1.0-spinNQC, 1./3.) );
  REAL8 Z2 = sqrt(3.0*spn2 + Z1*Z1);
  REAL8 rISCO = 3.0 + Z2 - signOfa*sqrt( (3.-Z1)*(3.+Z1 + 2.*Z2));
  REAL8 fISCO = pow(rISCO, -1.5)/(LAL_PI*mTScaled);

  if (debugPK){
      XLAL_PRINT_INFO("Stas - spin = %4.10f \n", spinNQC);
      XLAL_PRINT_INFO("Stas - NRPeakOmega22 =  %4.10f,   %4.10f \n",  GetNRSpinPeakOmegav2(2, 2, eta, spinNQC) / mTotal,  GetNRSpinPeakOmegav2(2, 2, eta, spinNQC));
      XLAL_PRINT_INFO("Stas ---- check for fmin NRPeakOmega22 = %4.10f, freqMinRad = %4.10f \n", NRPeakOmega22, freqMinRad);
      XLAL_PRINT_INFO("Stas -- minf freq is min( %4.10f, %4.10f )\n", NRPeakOmega22*0.1, freqMinRad);
      XLAL_PRINT_INFO("Stas -- initial radius (apr) %4.10f \n", pow(LAL_PI*fMin*mTScaled,-2./3.) );
      XLAL_PRINT_INFO("Stas -- omega_orb_0 = %4.10f \n",  rcrossrdotNorm/inner_product( rvec, rvec )/mTScaled);
      XLAL_PRINT_INFO("Stas -- rISCO  = %4.10f , freq_ISCO = %4.10f \n", rISCO, fISCO);

  }

    /* Check if initial frequency is too high: we choose an initial minimum separation
     *  of 10M as a compromise between reliability of initial conditions ans length
     *  of the waveform. We use Newtonian Kepler's law */
    if (pow(LAL_PI*fMin*mTScaled,-2./3.) < 10.0)
  {
      XLAL_PRINT_WARNING("Waveform generation may fail due to high starting frequency. The starting frequency corresponds to a small initial radius of %.2fM. We recommend a lower starting frequency that corresponds to an estimated starting radius > 10M.", pow(LAL_PI*fMin*mTScaled,-2.0/3.0));
  }
  if (fMin > freqMinRad){
        XLALPrintError("XLAL Error - %s: Intitial frequency is too high, the limit is %4.10f \n", __func__, freqMinRad);
        XLALDestroyREAL8Vector( sigmaKerr );
        XLALDestroyREAL8Vector( sigmaStar );
        XLALDestroyREAL8Vector( values );
        XLALDestroyREAL8Vector( dvalues );
        XLALDestroyREAL8Vector( rdMatchPoint );
        PRINT_PARAMS
        XLAL_ERROR (XLAL_EDOM );
  }

  switch ( SpinAlignedEOBversion )
  {
	  case 1:
	    if(debugPK)XLAL_PRINT_INFO("\t NQC: spin used = %.12e\n", spinNQC);
	    XLALSimIMRGetEOBCalibratedSpinNQC( &nqcCoeffs, 2, 2, eta, spinNQC );
	    break;
	  case 2:
	    if(debugPK)XLAL_PRINT_INFO("\t NQC: spins used = %.12e, %.12e\n", spinNQC, chiA);
	    XLALSimIMRGetEOBCalibratedSpinNQC3D( &nqcCoeffs, 2, 2, m1, m2, spinNQC, chiA );
	    break;
	  default:
        XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
        XLALDestroyREAL8Vector( sigmaKerr );
        XLALDestroyREAL8Vector( sigmaStar );
        XLALDestroyREAL8Vector( values );
        XLALDestroyREAL8Vector( dvalues );
        XLALDestroyREAL8Vector( rdMatchPoint );
        PRINT_PARAMS
        XLAL_ERROR( XLAL_EINVAL );
	    break;
  }

  /* NQC coeffs can be put to zero for debugging here */
  if (debugNoNQC) {
    nqcCoeffs.a1 = nqcCoeffs.a2 = nqcCoeffs.a3 = nqcCoeffs.a3S = nqcCoeffs.a4 =
      nqcCoeffs.a5 = nqcCoeffs.b1 = nqcCoeffs.b2 = nqcCoeffs.b3 = nqcCoeffs.b4 =
      0;
  }

  if ( debugPK )
  {
	  /* Print out all NQC coefficients */
    XLAL_PRINT_INFO(" NQC: a1 = %.16e, a2 = %.16e,\n a3 = %.16e, a3S = %.16e,\n a4 = %.16e, a5 = %.16e\n b1 = %.16e, b2 = %.16e,\n b3 = %.16e, b4 = %.16e\n",
        nqcCoeffs.a1, nqcCoeffs.a2, nqcCoeffs.a3,
        nqcCoeffs.a3S, nqcCoeffs.a4, nqcCoeffs.a5,
        nqcCoeffs.b1, nqcCoeffs.b2, nqcCoeffs.b3, nqcCoeffs.b4 );

	  /* Print out all mass parameters */
	  XLAL_PRINT_INFO("m1SI = %lf, m2SI = %lf, m1 = %lf, m2 = %lf\n",
			(double) m1SI, (double) m2SI, (double) m1, (double) m2 );
	  XLAL_PRINT_INFO("mTotal = %lf, mTScaled = %lf, eta = %lf\n",
			(double) mTotal, (double) mTScaled, (double) eta );
	  /* Print out all spin parameters */
	  XLAL_PRINT_INFO("spin1 = {%lf,%lf,%lf}, spin2 = {%lf,%lf,%lf}\n",
			(double) spin1[0], (double) spin1[1], (double) spin1[2],
			(double) spin2[0], (double) spin2[1], (double) spin2[2]);
	  XLAL_PRINT_INFO("mSpin1 = {%lf,%lf,%lf}, mSpin2 = {%lf,%lf,%lf}\n",
			(double) mSpin1[0], (double) mSpin1[1], (double) mSpin1[2],
			(double) mSpin2[0], (double) mSpin2[1], (double) mSpin2[2]);

      XLAL_PRINT_INFO("sigmaStar = {%lf,%lf,%lf}, sigmaKerr = {%lf,%lf,%lf}\n",
			(double) sigmaStar->data[0], (double) sigmaStar->data[1],
			(double) sigmaStar->data[2], (double) sigmaKerr->data[0],
			(double) sigmaKerr->data[1], (double) sigmaKerr->data[2]);
	  XLAL_PRINT_INFO("a = %lf, tplspin = %lf, chiS = %lf, chiA = %lf\n",
			(double) a, (double) tplspin, (double) chiS, (double) chiA);
	  XLAL_PRINT_INFO("s1Vec = {%lf,%lf,%lf}, s2Vec = {%lf,%lf,%lf}\n",
			(double) seobParams.s1Vec->data[0], (double) seobParams.s1Vec->data[1],
			(double) seobParams.s1Vec->data[2], (double) seobParams.s2Vec->data[0],
			(double) seobParams.s2Vec->data[1], (double) seobParams.s2Vec->data[2]);
	  XLAL_PRINT_INFO("a is used to compute Hamiltonian coefficients,\n tplspin and chiS and chiA for the multipole coefficients\n");

    XLAL_PRINT_INFO("s1VecOverM = {%.12e,%.12e,%.12e}\n", (double) s1VecOverMtMt.data[0],
	  (double) s1VecOverMtMt.data[1], (double) s1VecOverMtMt.data[2]);
	  XLAL_PRINT_INFO("s2VecOverM = {%.12e,%.12e,%.12e}\n", (double) s2VecOverMtMt.data[0],
	  (double) s2VecOverMtMt.data[1], (double) s2VecOverMtMt.data[2]);

    double StasS1 = sqrt(spin1[0]*spin1[0] + spin1[1]*spin1[1] +spin1[2]*spin1[2]);
    double StasS2 = sqrt(spin2[0]*spin2[0] + spin2[1]*spin2[1] +spin2[2]*spin2[2]);
    XLAL_PRINT_INFO("Stas: amplitude of spin1 = %.16e, amplitude of spin2 = %.16e, theta1 = %.16e , theta2 = %.16e, phi1 = %.16e, phi2 = %.16e  \n",
            StasS1, StasS2, acos(spin1[2]/StasS1)/LAL_PI, acos(spin2[2]/StasS2)/LAL_PI,
            atan2(spin1[1], spin1[0])/LAL_PI, atan2(spin2[1], spin2[0])/LAL_PI );
    fflush(NULL);
  }

/* *********************************************************************************
 * *********************************************************************************
 * STEP 2) Evolve EOB trajectory both at low and high sampling rate
 * *********************************************************************************
 * ********************************************************************************* */
  /* Low Sampling-rate integration */
  /* Initialize the GSL integrator */
    REAL8Vector *valuesV2 = NULL;
    if ( !(valuesV2 = XLALCreateREAL8Vector( 4 )) )
    {
        XLALDestroyREAL8Vector( sigmaKerr );
        XLALDestroyREAL8Vector( sigmaStar );
        XLALDestroyREAL8Vector( values );
        XLALDestroyREAL8Vector( dvalues );
        XLALDestroyREAL8Vector( rdMatchPoint );
        XLALPrintError("Failed to allocate REAL8Vector valuesV2!\n");
        PRINT_PARAMS
        XLAL_ERROR(  XLAL_ENOMEM );
    }
    memset( valuesV2->data, 0, valuesV2->length * sizeof( REAL8 ));
    /* If spins are almost aligned with LNhat, use SEOBNRv2 dynamics */
    if ( SpinsAlmostAligned ) {
        /* In SEOBNRv2 the dynamical variables are r, phi, p_r^*, p_phi */
        valuesV2->data[0] = values->data[0];
        valuesV2->data[1] = 0.;
        valuesV2->data[2] = values->data[3];
        valuesV2->data[3] = values->data[0] * values->data[4];

        seobParams.alignedSpins = 1;
        seobParams.chi1 = spin1Norm*copysign(1., cos(theta1Ini));
        seobParams.chi2 = spin2Norm*copysign(1., cos(theta2Ini));
        if (!(integrator = XLALAdaptiveRungeKutta4Init(4, XLALSpinAlignedHcapDerivative, XLALEOBSpinPrecAlignedStopCondition, EPS_ABS, EPS_REL)))
        {
            XLALDestroyREAL8Vector( sigmaKerr );
            XLALDestroyREAL8Vector( sigmaStar );
            XLALDestroyREAL8Vector( values );
            XLALDestroyREAL8Vector( dvalues );
            XLALDestroyREAL8Vector( rdMatchPoint );
            XLALDestroyREAL8Vector( valuesV2 );
            XLALPrintError("XLALAdaptiveRungeKutta4Init failed!\n");
            PRINT_PARAMS
            XLAL_ERROR( XLAL_EDOM );
        }
        seobParams.alignedSpins = 0;
    }
    else {
        if (!(integrator = XLALAdaptiveRungeKutta4Init(14,
                                                       XLALSpinPrecHcapNumericalDerivative, XLALEOBSpinPrecStopConditionBasedOnPR,
                                                       EPS_ABS, EPS_REL)))
        {
            XLALDestroyREAL8Vector( sigmaKerr );
            XLALDestroyREAL8Vector( sigmaStar );
            XLALDestroyREAL8Vector( values );
            XLALDestroyREAL8Vector( dvalues );
            XLALDestroyREAL8Vector( rdMatchPoint );
            XLALDestroyREAL8Vector( valuesV2 );
            XLALPrintError("XLALAdaptiveRungeKutta4Init failed!\n");
            PRINT_PARAMS
            XLAL_ERROR( XLAL_EDOM );
        }
    }

  /* Ensure that integration stops ONLY when the stopping condition is True */
  integrator->stopontestonly = 1;
  /* When this option is set to 0, the integration can be exceeddingly slow for spin-aligned systems */
  integrator->retries = 1;

  if( debugPK ){
    XLAL_PRINT_INFO("\n r = {%f,%f,%f}\n",
      values->data[0], values->data[1], values->data[2]);
    fflush(NULL);
  }

  if(debugPK) { XLAL_PRINT_INFO("\n\n BEGINNING THE EVOLUTION\n\n"); fflush(NULL); }

    REAL8Vector rVec, phiVec, prVec, pPhiVec;
    rVec.length = phiVec.length = prVec.length = pPhiVec.length = 0;
    rVec.data = 0; phiVec.data = 0; prVec.data = 0; pPhiVec.data = 0;
  /* Call the integrator */
    /* If spins are almost aligned with LNhat, use SEOBNRv2 dynamics */
    if ( SpinsAlmostAligned ) {
        seobParams.alignedSpins = 1;
        seobParams.chi1 = spin1Norm*copysign(1.,cos(theta1Ini));
        seobParams.chi2 = spin2Norm*copysign(1.,cos(theta2Ini));

        retLen = XLALAdaptiveRungeKutta4( integrator, &seobParams, valuesV2->data,
                                         0., 20./mTScaled, deltaT/mTScaled, &dynamicsV2 );
        seobParams.alignedSpins = 0;
        retLenLow = retLen;
        if ( retLenLow == XLAL_FAILURE )
        {
            XLALDestroyREAL8Vector( sigmaKerr );
            XLALDestroyREAL8Vector( sigmaStar );
            XLALDestroyREAL8Vector( values );
            XLALDestroyREAL8Vector( dvalues );
            XLALDestroyREAL8Vector( rdMatchPoint );
            XLALDestroyREAL8Vector( valuesV2 );
            if (dynamicsV2 != NULL){
                XLALDestroyREAL8Array( dynamicsV2 );
            }
            XLALPrintError("XLALAdaptiveRungeKutta4 failed!\n");
            PRINT_PARAMS
            XLAL_ERROR( XLAL_EFAILED );
        }
        tVec.length=rVec.length = phiVec.length = prVec.length = pPhiVec.length = retLenLow;
        tVec.data   = dynamicsV2->data;
        rVec.data    = dynamicsV2->data+retLenLow;
        phiVec.data  = dynamicsV2->data+2*retLenLow;
        prVec.data   = dynamicsV2->data+3*retLenLow;
        pPhiVec.data = dynamicsV2->data+4*retLenLow;
        dynamics = XLALCreateREAL8ArrayL( 2, 15, (UINT4)retLenLow );
        for (i = 0; i < retLenLow; i++) {
            dynamics->data[i] = tVec.data[i];
            dynamics->data[retLenLow + i] = rVec.data[i]*cos(phiVec.data[i]);
            dynamics->data[2*retLenLow + i] = rVec.data[i]*sin(phiVec.data[i]);
            dynamics->data[3*retLenLow + i] = 0.;
            dynamics->data[4*retLenLow + i] = prVec.data[i]*cos(phiVec.data[i]) - pPhiVec.data[i]/rVec.data[i]*sin(phiVec.data[i]);
            dynamics->data[5*retLenLow + i] = prVec.data[i]*sin(phiVec.data[i]) + pPhiVec.data[i]/rVec.data[i]*cos(phiVec.data[i]);
            dynamics->data[6*retLenLow + i] = 0.;
            dynamics->data[7*retLenLow + i] = 0.;
            dynamics->data[8*retLenLow + i] = 0.;
            dynamics->data[9*retLenLow + i] = spin1[2]*(m1*m1/mTotal/mTotal);
            dynamics->data[10*retLenLow + i] = 0.;
            dynamics->data[11*retLenLow + i] = 0.;
            dynamics->data[12*retLenLow + i] = spin2[2]*(m2*m2/mTotal/mTotal);
            dynamics->data[13*retLenLow + i]= phiVec.data[i];
            dynamics->data[14*retLenLow + i]= 0.;
        }
    }
    else {
        retLen = XLALAdaptiveRungeKutta4( integrator, &seobParams, values->data,
                                         0., 20./mTScaled, deltaT/mTScaled, &dynamics );
        retLenLow = retLen;
        if ( retLenLow == XLAL_FAILURE )
        {
            XLALPrintError("XLALAdaptiveRungeKutta4 failed!\n");
            PRINT_PARAMS
            XLAL_ERROR( XLAL_EFAILED );
        }
        tVec.length = retLenLow;
        tVec.data   = dynamics->data;
    }

    /* This is common to all cases, even when we use SEOBNRv2 dynamics */
        posVecx.data = dynamics->data+retLenLow;
        posVecy.data = dynamics->data+2*retLenLow;
        posVecz.data = dynamics->data+3*retLenLow;
        momVecx.data = dynamics->data+4*retLenLow;
        momVecy.data = dynamics->data+5*retLenLow;
        momVecz.data = dynamics->data+6*retLenLow;
        s1Vecx.data = dynamics->data+7*retLenLow;
        s1Vecy.data = dynamics->data+8*retLenLow;
        s1Vecz.data = dynamics->data+9*retLenLow;
        s2Vecx.data = dynamics->data+10*retLenLow;
        s2Vecy.data = dynamics->data+11*retLenLow;
        s2Vecz.data = dynamics->data+12*retLenLow;
        phiDMod.data= dynamics->data+13*retLenLow;
        phiMod.data = dynamics->data+14*retLenLow;
     if(debugPK) { XLAL_PRINT_INFO("\n\n FINISHED THE EVOLUTION\n\n"); fflush(NULL);  }

  if (debugPK) {
    /* Write the dynamics to file */
    out = fopen( "seobDynamics.dat", "w" );
    for ( i = 0; i < retLenLow; i++ )
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

    XLAL_PRINT_INFO("AT We are here %d %d\n",retLenLow, nStepBack);
    fflush(NULL);
  }

  /* High Sampling-rate integration  */
  /* Stepback by 150M wrt the ending time of the low-sampling-rate trajecotry */
  /* If 150M of step back > actual time of evolution, step back 50% of that */
  if (tStepBack > retLenLow*deltaT)
  {
    tStepBack = 0.5*retLenLow*deltaT;
    nStepBack = ceil( tStepBack / deltaT );
  }

  /* Step back in time by tStepBack and volve EOB trajectory again
   * using high sampling rate. */
  hiSRndx    = retLenLow - nStepBack;
  deltaTHigh = deltaT / (REAL8)resampFac;
  HiSRstart  = tVec.data[hiSRndx];

  if (debugPK) {
    XLAL_PRINT_INFO("AT We are here %d\n",nStepBack);
    XLAL_PRINT_INFO("Stas: start HighSR at %.16e \n", HiSRstart);
    XLAL_PRINT_INFO( "Stepping back %d points - we expect %d points at high SR\n", nStepBack, nStepBack*resampFac );
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
    XLAL_PRINT_INFO( "Commencing high SR integration.. From: \n" );
    for( i=0; i<12; i++)XLAL_PRINT_INFO("%.16e\n", values->data[i]);
    fflush(NULL);
  }
  seobParams.prev_dr=0.;
  /* For HiSR evolution, we stop whenver XLALEOBSpinPrecStopConditionBasedOnPR is triggered */
  integrator->stop = XLALEOBSpinPrecStopConditionBasedOnPR;

  REAL8Vector rVecHi, phiVecHi, prVecHi, pPhiVecHi;
  rVecHi.length = phiVecHi.length = prVecHi.length = pPhiVecHi.length = 0;
  rVecHi.data = 0; phiVecHi.data = 0; prVecHi.data = 0; pPhiVecHi.data = 0;

    /** If spins are aligned/antialigned with LNhat to within 1e-4 rads, then use SEOBNRv2 dynamics */
  if ( SpinsAlmostAligned ) {
        seobParams.alignedSpins = 1;
        seobParams.chi1 = spin1Norm*copysign(1., cos(theta1Ini));
        seobParams.chi2 = spin2Norm*copysign(1.,cos(theta2Ini));

        valuesV2->data[0] = rVec.data[hiSRndx];
        valuesV2->data[1] = phiVec.data[hiSRndx];
        valuesV2->data[2] = prVec.data[hiSRndx];
        valuesV2->data[3] = pPhiVec.data[hiSRndx];
        if (debugPK) {XLAL_PRINT_INFO("Start high SR integration at: valuesV2->data[0], valuesV2->data[1], valuesV2->data[2], valuesV2->data[3] %e %e %e %e\n", valuesV2->data[0], valuesV2->data[1], valuesV2->data[2], valuesV2->data[3]);}

        integrator->stop = XLALSpinPrecAlignedHiSRStopCondition;

        retLen = XLALAdaptiveRungeKutta4( integrator, &seobParams, valuesV2->data, 0., tStepBack/mTScaled, deltaTHigh/mTScaled, &dynamicsV2Hi );

        seobParams.alignedSpins = 0;
        retLenHi = retLen;
        if ( retLenHi == XLAL_FAILURE )
        {
            XLALDestroyREAL8Vector( sigmaKerr );
            XLALDestroyREAL8Vector( sigmaStar );
            XLALDestroyREAL8Vector( values );
            XLALDestroyREAL8Vector( dvalues );
            XLALDestroyREAL8Vector( rdMatchPoint );
            XLALDestroyREAL8Vector( valuesV2 );
            if (dynamicsV2 != NULL){
                XLALDestroyREAL8Array( dynamicsV2 );
            }
            if (dynamicsV2Hi != NULL){
                XLALDestroyREAL8Array( dynamicsV2Hi );
            }
            XLALPrintError("XLALAdaptiveRungeKutta4 failed!\n");
            PRINT_PARAMS
            XLAL_ERROR( XLAL_EFAILED );
        }
        timeHi.length = phiVecHi.length = prVecHi.length = pPhiVecHi.length = retLenHi;
        timeHi.data   = dynamicsV2Hi->data;
        rVecHi.data    = dynamicsV2Hi->data+retLenHi;
        phiVecHi.data  = dynamicsV2Hi->data+2*retLenHi;
        prVecHi.data   = dynamicsV2Hi->data+3*retLenHi;
        pPhiVecHi.data = dynamicsV2Hi->data+4*retLenHi;
      dynamicsHi = XLALCreateREAL8ArrayL( 2, 15, (UINT4)retLenHi );
      for (i = 0; i < retLenHi; i++) {
          dynamicsHi->data[i] = timeHi.data[i];
          dynamicsHi->data[retLenHi + i] = rVecHi.data[i]*cos(phiVecHi.data[i]);
          dynamicsHi->data[2*retLenHi + i]  = rVecHi.data[i]*sin(phiVecHi.data[i]);
          dynamicsHi->data[3*retLenHi + i] = 0.;
          dynamicsHi->data[4*retLenHi + i] = prVecHi.data[i]*cos(phiVecHi.data[i]) - pPhiVecHi.data[i]/rVecHi.data[i]*sin(phiVecHi.data[i]);
          dynamicsHi->data[5*retLenHi + i] = prVecHi.data[i]*sin(phiVecHi.data[i]) + pPhiVecHi.data[i]/rVecHi.data[i]*cos(phiVecHi.data[i]);
          dynamicsHi->data[6*retLenHi + i] = 0.;
          dynamicsHi->data[7*retLenHi + i] = 0.;
          dynamicsHi->data[8*retLenHi + i] = 0.;
          dynamicsHi->data[9*retLenHi + i] = spin1[2]*(m1*m1/mTotal/mTotal);
          dynamicsHi->data[10*retLenHi + i] = 0.;
          dynamicsHi->data[11*retLenHi + i] = 0.;
          dynamicsHi->data[12*retLenHi + i] = spin2[2]*(m2*m2/mTotal/mTotal);
          dynamicsHi->data[13*retLenHi + i]= phiVecHi.data[i];
          dynamicsHi->data[14*retLenHi + i]  = 0.;
      }
     }
    else {
        retLen = XLALAdaptiveRungeKutta4( integrator, &seobParams, values->data,
                                         0., tStepBack/mTScaled, deltaTHigh/mTScaled, &dynamicsHi );
        retLenHi = retLen;
        if ( retLenHi == XLAL_FAILURE )
        {
            XLALDestroyREAL8Vector( sigmaKerr );
            XLALDestroyREAL8Vector( sigmaStar );
            XLALDestroyREAL8Vector( values );
            XLALDestroyREAL8Vector( dvalues );
            XLALDestroyREAL8Vector( rdMatchPoint );
            XLALDestroyREAL8Vector( valuesV2 );
            if (dynamicsV2 != NULL){
                XLALDestroyREAL8Array( dynamicsV2 );
            }
            if (dynamicsV2Hi != NULL){
                XLALDestroyREAL8Array( dynamicsV2Hi );
            }
            if (dynamics != NULL){
                XLALDestroyREAL8Array( dynamics );
            }
            if (dynamicsHi != NULL){
                XLALDestroyREAL8Array( dynamicsHi );
            }
            XLALPrintError("XLALAdaptiveRungeKutta4 failed!\n");
            PRINT_PARAMS
            XLAL_ERROR( XLAL_EFAILED );
        }
        timeHi.length = retLenHi;
        timeHi.data   = dynamicsHi->data;
    }
    retLenHi = retLen;
    if(debugPK){XLAL_PRINT_INFO("retLenHi = %d\n", retLenHi);}
    posVecxHi.length = posVecyHi.length = posVeczHi.length =
    momVecxHi.length = momVecyHi.length = momVeczHi.length =
    s1VecxHi.length = s1VecyHi.length = s1VeczHi.length =
    s2VecxHi.length = s2VecyHi.length = s2VeczHi.length =
    phiDModHi.length = phiModHi.length = retLenHi;

    /* This is common to all cases, even when we use SEOBNRv2 dynamics */
        posVecxHi.data = dynamicsHi->data+retLenHi;
        posVecyHi.data = dynamicsHi->data+2*retLenHi;
        posVeczHi.data = dynamicsHi->data+3*retLenHi;
        momVecxHi.data = dynamicsHi->data+4*retLenHi;
        momVecyHi.data = dynamicsHi->data+5*retLenHi;
        momVeczHi.data = dynamicsHi->data+6*retLenHi;
        s1VecxHi.data = dynamicsHi->data+7*retLenHi;
        s1VecyHi.data = dynamicsHi->data+8*retLenHi;
        s1VeczHi.data = dynamicsHi->data+9*retLenHi;
        s2VecxHi.data = dynamicsHi->data+10*retLenHi;
        s2VecyHi.data = dynamicsHi->data+11*retLenHi;
        s2VeczHi.data = dynamicsHi->data+12*retLenHi;
        phiDModHi.data= dynamicsHi->data+13*retLenHi;
        phiModHi.data = dynamicsHi->data+14*retLenHi;


  if(debugPK) { XLAL_PRINT_INFO( "Finished high SR integration... \n" ); fflush(NULL); }

   if (debugPK){
    out = fopen( "seobDynamicsHi.dat", "w" );
    for ( i = 0; i < retLenHi; i++ )
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


/* *********************************************************************************
 * *********************************************************************************
 * STEP 3) Compute Euler angles to go from initial inertial frame to
 * precessing frame: Eqs. 19-20 of PRD 89, 084006 (2014)
* *********************************************************************************
* **********************************************************************************/
/* Interpolate trajectories to compute L_N (t) in order to get alpha(t) and
   * beta(t). */
    INT4 phaseCounterA = 0, phaseCounterB = 0;
    Alpha = XLALCreateREAL8Vector( retLenLow );
    Beta   = XLALCreateREAL8Vector( retLenLow );
    Gamma = XLALCreateREAL8Vector( retLenLow );
    EulerAnglesI2P(Alpha, Beta, Gamma, &phaseCounterA, &phaseCounterB, tVec, posVecx, posVecy, posVecz, retLenLow, 0., 0.5*LAL_PI, 0);
    //EulerAnglesI2P(Alpha, Beta, Gamma, &phaseCounterA, &phaseCounterB, tVec, posVecx, posVecy, posVecz, retLenLow, 0., 0.0, 0);

    if (debugPK){
        XLAL_PRINT_INFO("Writing Alpha and Beta angle timeseries at low SR to alphaANDbeta.dat\n" );
        fflush(NULL);
        out = fopen( "alphaANDbeta.dat","w");
        for( i = 0; i < retLenLow; i++ ) {
            fprintf( out, "%.16e %.16e %.16e\n", tVec.data[i], Alpha->data[i], Beta->data[i]);
        }
        fclose(out);
        
        XLAL_PRINT_INFO("Writing Gamma angle timeseries at low SR to gamma.dat\n");
        fflush(NULL);
        out = fopen( "gamma.dat","w");
        for( i = 0; i < retLenLow; i++ ) {
            fprintf( out, "%.16e %.16e\n", tVec.data[i], Gamma->data[i]);
        }
        fclose(out);
    }

    AlphaHi = XLALCreateREAL8Vector( retLenHi );
    BetaHi   = XLALCreateREAL8Vector( retLenHi );
    GammaHi = XLALCreateREAL8Vector( retLenHi );
    EulerAnglesI2P(AlphaHi, BetaHi, GammaHi, &phaseCounterA, &phaseCounterB, timeHi, posVecxHi, posVecyHi, posVeczHi, retLenHi, Alpha->data[hiSRndx], Gamma->data[hiSRndx], 1);
    if (debugPK){
        XLAL_PRINT_INFO("Writing Alpha and Beta angle timeseries at high SR to alphaANDbetaHi.dat\n" );
        fflush(NULL);
        out = fopen( "alphaANDbetaHi.dat","w");
        for( i = 0; i < retLenHi; i++ ) {
            fprintf( out, "%.16e %.16e %.16e\n", tVec.data[hiSRndx] + timeHi.data[i], AlphaHi->data[i], BetaHi->data[i]);
        }
        fclose(out);
        
        XLAL_PRINT_INFO("Writing Gamma angle timeseries at high SR to gammaHi.dat\n");
        fflush(NULL);
        out = fopen( "gammaHi.dat","w");
        for( i = 0; i < retLenHi; i++ ) {
            fprintf( out, "%.16e %.16e\n", tVec.data[hiSRndx] + timeHi.data[i], GammaHi->data[i]);
        }
        fclose(out);
    }

  /* Compute leading QNM to determine how long the ringdown will be */
  modefreqVec.length = 1;
  modefreqVec.data   = &modeFreq;

  if ( XLALSimIMREOBGenerateQNMFreqV2Prec( &modefreqVec, m1, m2, spin1, spin2, 2, 2, 1, spinEOBApproximant ) == XLAL_FAILURE ) /* This function returns XLAL_SUCCESS or calls XLAL_ERROR( XLAL_EINVAL ). Internally, it calls XLALSimIMREOBFinalMassSpinPrec which has the same behavior */ {
     FREE_EVERYTHING
     FREE_SPHHARM
      XLALPrintError("XLALSimIMREOBGenerateQNMFreqV2Prec failed!\n");
      PRINT_PARAMS
      XLAL_ERROR( XLAL_EFUNC );
  }

  retLenRDPatchHi= (UINT4)ceil( 40 / ( cimag(modeFreq) * deltaTHigh ));
  retLenRDPatchLow = (UINT4)ceil( 40 / ( cimag(modeFreq) * deltaT ));

/* *********************************************************************************
 * *********************************************************************************
 * STEP 4) Locate merger point and at that time calculate J, chi and kappa,
 * and construct final J frame
 * *********************************************************************************
 * **********************************************************************************/
  /* Locate merger point */
    // Find tAmpMax to determin the tAttachment
    REAL8Vector *radiusVec = NULL;
    if ( !(radiusVec = XLALCreateREAL8Vector( retLenHi )) ) {
        FREE_EVERYTHING
        FREE_SPHHARM
        XLALPrintError("Failed to allocate REAL8Vector radiusVec!\n");
        PRINT_PARAMS
        XLAL_ERROR( XLAL_ENOMEM );
    }

    for ( i = 0; i < retLenHi; i++ )
    {
        for ( j = 0; j < 3; j++ )
        {
            values->data[j] = dynamicsHi->data[(j+1)*retLenHi + i];
        }
        radiusVec->data[i] = sqrt(values->data[0]*values->data[0] + values->data[1]*values->data[1] + values->data[2]*values->data[2]);
    }


    int foundPeakOmega = 0;
    if (debugPK) {
        XLAL_PRINT_INFO("Stas searching for maxima in omega .... \n");
    }
    REAL8 tMaxOmega;
    tPeakOmega = XLALSimLocateOmegaTime(dynamicsHi, values->length, retLenHi, seobParams, seobCoeffs, m1, m2, radiusVec, &foundPeakOmega, &tMaxOmega);

    if(tPeakOmega == 0.0 || foundPeakOmega==0){
        if (debugPK){
            XLAL_PRINT_INFO("maximum of omega and A(r) is not found, looking for amplitude maximum\n");
        }
    }
    else{
        if (debugPK){
            XLAL_PRINT_INFO("The omega-related time is found and it's %f \n", tPeakOmega);
        }
    }
    if(debugPK) {
      XLAL_PRINT_INFO( "Estimation of the peak is now at time %.16e, %.16e \n",
            tPeakOmega, tPeakOmega+HiSRstart);
      fflush(NULL);
    }


  /* Calculate J at merger */
  spline = gsl_spline_alloc( gsl_interp_cspline, retLenHi );
  acc    = gsl_interp_accel_alloc();

  /* Obtain the dynamics at the time where Omega = d\phi/dt peaks */
  for ( j = 0; j < values->length; j++ )
  {
    gsl_spline_init( spline,
                    dynamicsHi->data, dynamicsHi->data+(j+1)*retLenHi, retLenHi );
    values->data[j] = gsl_spline_eval( spline, tPeakOmega, acc );
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  /* Get the phase offset required to ensure that the orbital phase = phiC at
   * tPeakOmega */
  //coa_phase_offset = values->data[12] - phiC;

  /* Calculate dr/dt */
  memset( dvalues->data, 0, 14*sizeof(dvalues->data[0]));
  if( XLALSpinPrecHcapRvecDerivative( 0, values->data, dvalues->data,
        &seobParams) != XLAL_SUCCESS )  {
      XLAL_PRINT_INFO( " Calculation of dr/dt at t = tPeakOmega \n");
      FREE_EVERYTHING
      FREE_SPHHARM
      XLALPrintError("XLALSpinPrecHcapRvecDerivative failed!\n");
      PRINT_PARAMS
      XLAL_ERROR( XLAL_EDOM );
  }

  /* Calculate Omega = r x dr/dt */
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
    XLAL_PRINT_INFO("J at merger: %e, %e, %e (mag = %e)\n", Jx, Jy, Jz, magJ);
    fflush(NULL);
  }

  /* Calculate chi and kappa at merger: Eq. 26 and 27 of PRD 89, 084006 (2014).
   Here we project the spins onto J */
  chi1J = values->data[6]*Jx + values->data[7] *Jy + values->data[8] *Jz;
  chi2J = values->data[9]*Jx + values->data[10]*Jy + values->data[11]*Jz;
  chi1J/= magJ*m1*m1/mTotal/mTotal;
  chi2J/= magJ*m2*m2/mTotal/mTotal;
  chiJ = (chi1J+chi2J)/2. + (chi1J-chi2J)/2.*((m1-m2)/(m1+m2))/(1. - 2.*eta);
  kappaJL = (Lx*Jx + Ly*Jy + Lz*Jz) / magL / magJ;
  magLN = sqrt(inner_product(rcrossrdot, rcrossrdot));
  JLN =  (Jx*rcrossrdot[0] + Jy*rcrossrdot[1] + Jz*rcrossrdot[2])/magLN/magJ;
  // Stas: here is the second option which we will use: projection of the spin on
  // L at l.r. projection on J doesn't work for anti-aligned system
  chi1L = values->data[6]*Lx + values->data[7] *Ly + values->data[8] *Lz;
  chi2L = values->data[9]*Lx + values->data[10]*Ly + values->data[11]*Lz;
  chi1L /=  magL*m1*m1/mTotal/mTotal;
  chi2L /=  magL*m2*m2/mTotal/mTotal;
  chiL = (chi1L+chi2L)/2. + (chi1L-chi2L)/2.*((m1-m2)/(m1+m2))/(1. - 2.*eta);


  if(debugPK) {
      XLAL_PRINT_INFO("chi1J,chi2J,chiJ = %3.10f %3.10f %3.10f\n", chi1J,chi2J,chiJ); fflush(NULL); fflush(NULL);
      XLAL_PRINT_INFO("chi1L,chi2L,chiL = %3.10f %3.10f %3.10f\n", chi1L,chi2L,chiL); fflush(NULL); fflush(NULL);
      XLAL_PRINT_INFO("J.L = %4.11f \n", kappaJL); fflush(NULL);
      XLAL_PRINT_INFO("J.LN = %4.11f \n", JLN);
      XLAL_PRINT_INFO("L.LN = %4.11f \n", (Lx*rcrossrdot[0] + Ly*rcrossrdot[1] + Lz*rcrossrdot[2])/magL/magLN);
  }


  sh = 0.0;
  /* Calculate combsize and deltaNQC */
  switch ( SpinAlignedEOBversion )
  {
    /* Eqs. 33-34 of PRD 86, 024011 (2012) */
    case 1:
      combSize = 7.5;
      deltaNQC = XLALSimIMREOBGetNRSpinPeakDeltaT(2, 2, eta,  chiJ);
      if ( debugPK ) {
        XLAL_PRINT_INFO("v1 RD prescriptions are used! %3.10f %3.10f\n",
                combSize, deltaNQC);
        fflush(NULL);
      }
      break;
     /* The following if's are documented in sections "DeltaNQC" and "￼PseudoQNM prescriptions" of  https://dcc.ligo.org/T1400476 */
    case 2:
      combSize = 12.;
      if ( chi1L == 0. && chi2L == 0. ) combSize = 11.;
      if (chiL >= 0.8 && eta >30.0/(31.*31) && eta <10.0/121.) combSize = 13.5;
      if (chiL >= 0.9 && eta < 30./(31.*31.)) combSize = 12.0;
      if (chiL >= 0.8 && eta > 10./121.) combSize = 8.5;
      deltaNQC = XLALSimIMREOBGetNRSpinPeakDeltaTv2(2, 2, m1, m2, chi1L, chi2L );
      if ( debugPK ) {
        XLAL_PRINT_INFO("v2 RD prescriptions are used! %3.10f %3.10f\n",
                combSize, deltaNQC);
        fflush(NULL);
      }
      break;
    default:
      XLALPrintError( "XLAL Error - %s: wrong SpinAlignedEOBversion value, must be 1 or 2!\n", __func__ );
      FREE_EVERYTHING
      FREE_SPHHARM
      PRINT_PARAMS
      XLAL_ERROR( XLAL_EINVAL );
      break;
  }


  /* NOTE: tAttach is further modified by small shift "sh" computed and applied in XLALSimIMREOBHybridAttachRingdownPrec !!!
    Manually choose tAttach here */
  tAttach = tPeakOmega - deltaNQC;

  /* tPeakOmega is initialized to 0 by default, so even when found is false, tPeakOmega has numerical value */
  if ( ! foundPeakOmega ){
     tAttach = tPeakOmega;
  }

  if (debugPK){
     XLAL_PRINT_INFO("tPeakOmega - deltaNQC, tPeakOmega, deltaNQC %e %e %e\n", tPeakOmega - deltaNQC + HiSRstart, tPeakOmega + HiSRstart, deltaNQC);
    XLAL_PRINT_INFO("For RD: DeltaNQC = %3.10f, comb = %3.10f \n", deltaNQC, combSize);
    XLAL_PRINT_INFO("NOTE! that additional shift (sh) is computed and added in XLALSimIMREOBHybridAttachRingdownPrec\n");
	  fflush(NULL);
  }

  /* Manually choose combSize */
  //combSize = 40.0;

  /* Construct J-frame */
  JframeEz[0] = Jx / magJ;
  JframeEz[1] = Jy / magJ;
  JframeEz[2] = Jz / magJ;

  if ( fabs(1.+ JframeEz[2]) <= 1.0e-13 )
    { JframeEx[0] = -1.; JframeEx[1] = 0.; JframeEx[2] = 0.; } // anti-aligned
   //{ JframeEx[0] = 0.0; JframeEx[1] = -1.0; JframeEx[2] = 0.; } // anti-aligned
  else {
      if ( fabs(1. - JframeEz[2]) <= 1.0e-13 )
          //{ JframeEx[0] = 1.; JframeEx[1] = 0.; JframeEx[2] = 0.; }  // aligned
          { JframeEx[0] = 1.0; JframeEx[1] = 0.0; JframeEx[2] = 0.; }  // aligned
      else {
            JframeEx[0] = JframeEz[1];
            JframeEx[1] = -JframeEz[0];
            JframeEx[2] = 0.;
      }
  }

  JExnorm = sqrt(JframeEx[0]*JframeEx[0] + JframeEx[1]*JframeEx[1]);

  JframeEx[0] /= JExnorm;
  JframeEx[1] /= JExnorm;
  JframeEy[0] = JframeEz[1]*JframeEx[2] - JframeEz[2]*JframeEx[1];
  JframeEy[1] = JframeEz[2]*JframeEx[0] - JframeEz[0]*JframeEx[2];
  JframeEy[2] = JframeEz[0]*JframeEx[1] - JframeEz[1]*JframeEx[0];


  if(debugPK) {
    XLAL_PRINT_INFO("J-frameEx = [%.16e\t%.16e\t%.16e]\n", JframeEx[0], JframeEx[1], JframeEx[2]);XLAL_PRINT_INFO("J-frameEy = [%.16e\t%.16e\t%.16e]\n", JframeEy[0], JframeEy[1], JframeEy[2]);
    XLAL_PRINT_INFO("J-frameEz = [%.16e\t%.16e\t%.16e]\n", JframeEz[0], JframeEz[1], JframeEz[2]);fflush(NULL);
  }

/* *********************************************************************************
 * *********************************************************************************
 * STEP 5) Generate quasi-nonprecessing waveforms in precessing frame
 * *********************************************************************************
 * **********************************************************************************/
  /* Create time-series containers for euler angles and hlm harmonics */
    if ( !(alphaI2PTS = XLALCreateREAL8TimeSeries( "alphaI2P", &tc, 0.0, deltaT, &lalStrainUnit, retLenLow )) + !(betaI2PTS = XLALCreateREAL8TimeSeries(  "betaI2P", &tc, 0.0, deltaT, &lalStrainUnit, retLenLow )) + !(gammaI2PTS = XLALCreateREAL8TimeSeries( "gammaI2P", &tc, 0.0, deltaT, &lalStrainUnit, retLenLow )) + !(alphaP2JTS = XLALCreateREAL8TimeSeries( "alphaP2J", &tc, 0.0, deltaT, &lalStrainUnit, retLenLow )) + !(betaP2JTS = XLALCreateREAL8TimeSeries(  "betaP2J", &tc, 0.0, deltaT, &lalStrainUnit, retLenLow )) + !(gammaP2JTS = XLALCreateREAL8TimeSeries( "gammaP2J", &tc, 0.0, deltaT, &lalStrainUnit, retLenLow )) ) {
        FREE_EVERYTHING
        FREE_SPHHARM
        XLALPrintError("Failed to allocate alphaI2PTS or betaI2PTS or gammaI2PTS or alphaP2JTS or betaP2JTS or gammaP2JTS!\n");
        PRINT_PARAMS
        XLAL_ERROR(  XLAL_ENOMEM );
    }
    
    if ( !(h22TS   = XLALCreateCOMPLEX16TimeSeries( "H_22",  &tc, 0.0, deltaT, &lalStrainUnit, retLenLow )) + !(h21TS   = XLALCreateCOMPLEX16TimeSeries( "H_21", &tc, 0.0, deltaT, &lalStrainUnit, retLenLow )) + !(h20TS   = XLALCreateCOMPLEX16TimeSeries( "H_20",  &tc, 0.0, deltaT, &lalStrainUnit, retLenLow )) + !(h2m1TS  = XLALCreateCOMPLEX16TimeSeries( "H_2m1", &tc, 0.0, deltaT, &lalStrainUnit, retLenLow )) + !(h2m2TS  = XLALCreateCOMPLEX16TimeSeries( "H_2m2", &tc, 0.0, deltaT, &lalStrainUnit, retLenLow ))) {
        FREE_EVERYTHING
        FREE_SPHHARM
        XLALPrintError("Failed to allocate h22TS or h21TS or h20TS or h2m1TS or h2m2TS!\n");
        PRINT_PARAMS
        XLAL_ERROR(  XLAL_ENOMEM );
     }

  hIMRJTS    = XLALCreateCOMPLEX16TimeSeries( "HIMRJ",
              &tc, 0.0, deltaT, &lalStrainUnit, retLenLow + retLenRDPatchLow );


  if ( !(tlist = XLALCreateREAL8Vector( retLenLow ))
    || !(tlistRDPatch = XLALCreateREAL8Vector( retLenLow + retLenRDPatchLow )) )
  {
    XLALPrintError("Failed to allocate REAL8Vector tlistRDPatch!\n");
    PRINT_PARAMS
    XLAL_ERROR(  XLAL_ENOMEM );
  }
  // timeJFullHi will be a copy  of tlistRDPatchHi which will be  returned together with J-frame modes
  if ( !(timeJFull = XLALCreateREAL8Vector( retLenLow + retLenRDPatchLow )) )
  {
    XLALPrintError("Failed to allocate REAL8Vector timeJFull!\n");
    PRINT_PARAMS
    XLAL_ERROR(  XLAL_ENOMEM );
  }
  // timeIFullHi will be a copy  of tlistRDPatchHi which will be  returned together with I-frame modes
  if ( !(timeIFull = XLALCreateREAL8Vector( retLenLow + retLenRDPatchLow )) )
  {
    XLALPrintError("Failed to allocate REAL8Vector timeIFull!\n");
    PRINT_PARAMS
    XLAL_ERROR(  XLAL_ENOMEM );
  }
  memset( tlist->data, 0, tlist->length * sizeof( REAL8 ));
  memset( tlistRDPatch->data, 0, tlistRDPatch->length * sizeof( REAL8 ));
  memset( timeJFull->data, 0, timeJFull->length * sizeof( REAL8 ));
  memset( timeIFull->data, 0, timeIFull->length * sizeof( REAL8 ));
  for ( i = 0; i < retLenLow; i++ )
  {
    tlist->data[i] = i * deltaT/mTScaled;
    tlistRDPatch->data[i] = i * deltaT/mTScaled;
  }
  for ( i = retLenLow; i < retLenLow + retLenRDPatchLow; i++ )
  {
    tlistRDPatch->data[i] = i * deltaT/mTScaled;
  }

  /* Main loop for quasi-nonprecessing waveform generation */
  // Generating modes for coarsely sampled portion
  //FILE *out2 = fopen("Lframe.dat", "w");  //FIXME
  if(debugPK){
    XLAL_PRINT_INFO("Generating precessing-frame modes for coarse dynamics\n");
    fflush(NULL);
    out = fopen( "rotDynamics.dat", "w" );
  }

  for ( i = 0; i < retLenLow; i++ )
  {
    for ( j = 0; j < values->length; j++ )
    {
      values->data[j] = dynamics->data[(j+1)*retLenLow + i];
    }

    /* Calculate dr/dt */
    memset( dvalues->data, 0, 14*sizeof(dvalues->data[0]));
    if( XLALSpinPrecHcapRvecDerivative( 0, values->data, dvalues->data,
          &seobParams) != XLAL_SUCCESS )
    {
      XLAL_PRINT_INFO( " Calculation of dr/dt at t = tPeakOmega \n");
        FREE_EVERYTHING
        FREE_SPHHARM
        XLALPrintError("XLALSpinPrecHcapRvecDerivative failed!\n");
        PRINT_PARAMS
        XLAL_ERROR( XLAL_EDOM );
    }

    /* Calculate omega */
    memcpy( rdotvec, dvalues->data, 3*sizeof(REAL8));
    rvec[0] = posVecx.data[i]; rvec[1] = posVecy.data[i];
    rvec[2] = posVecz.data[i];
    cross_product( rvec, rdotvec, rcrossrdot );
    magR = sqrt(inner_product(rvec, rvec));
    omega = sqrt(inner_product(rcrossrdot, rcrossrdot)) / (magR*magR);
    v = cbrt( omega );

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

    /* this one is defined w.r.t. L not LN*/
    aI2P = Alpha->data[i];
    bI2P = Beta->data[i];
    gI2P = Gamma->data[i];

    //if (debugPK){
        //REAL8 LframeEx[3] = {0,0,0}, LframeEy[3] = {0,0,0}, LframeEz[3] = {0,0,0};
        //LframeEx[0] =  cos(aI2P)*cos(bI2P)*cos(gI2P) - sin(aI2P)*sin(gI2P);
       // LframeEx[1] =  sin(aI2P)*cos(bI2P)*cos(gI2P) + cos(aI2P)*sin(gI2P);
        //LframeEx[2] = -sin(bI2P)*cos(gI2P);
        //LframeEy[0] = -cos(aI2P)*cos(bI2P)*sin(gI2P) - sin(aI2P)*cos(gI2P);
        //LframeEy[1] = -sin(aI2P)*cos(bI2P)*sin(gI2P) + cos(aI2P)*cos(gI2P);
        //LframeEy[2] =  sin(bI2P)*sin(gI2P);
        //LframeEz[0] =  LNhx;
        //LframeEz[1] =  LNhy;
        //LframeEz[2] =  LNhz;
        
        //fprintf(out2, "%.16e   %.16e  %.16e  %.16e    %.16e  %.16e  %.16e     %.16e  %.16e  %.16e \n",
        //        i*deltaT/mTScaled, LframeEx[0],  LframeEx[1],  LframeEx[2],  LframeEy[0],  LframeEy[1], 
        //         LframeEy[2],  LframeEz[0],  LframeEz[1],  LframeEz[2]);

    //}

    EulerAnglesP2J(&aP2J, &bP2J, &gP2J, aI2P, bI2P, gI2P, LNhx, LNhy, LNhz, JframeEx, JframeEy, JframeEz);
                     
      /** Euler angles to go from precessing to J-frame. Note that we follow in this code passive rotation Z-Y-Z and notations
       * of Arun et al. arXiv 0810.5336, however the Wiegner D-matrix and mode rotation coded up in LAL for active rotation 
       * We make the angle transformation here in order to match two conventions  */
    alphaI2PTS->data->data[i] = -aI2P;
    betaI2PTS->data->data[i] = bI2P;
    gammaI2PTS->data->data[i] = -gI2P;
    alphaP2JTS->data->data[i] = -aP2J;
    betaP2JTS->data->data[i] = bP2J;
    gammaP2JTS->data->data[i] = -gP2J;


    /* Calculate the value of the Hamiltonian */
    memcpy( cartPosVec.data, values->data,   3*sizeof(REAL8) );
    memcpy( cartMomVec.data, values->data+3, 3*sizeof(REAL8) );

    /* Update Hamiltonian coefficients as per |Skerr| */
      s1Vec.data = s1Data;
      s2Vec.data = s2Data;
      
      s1VecOverMtMt.data = s1DataNorm;
      s2VecOverMtMt.data = s2DataNorm;
    for( k = 0; k < 3; k++ )
    {
			s1VecOverMtMt.data[k] = values->data[k+6];
			s2VecOverMtMt.data[k] = values->data[k+9];
			s1Vec.data[k] = s1VecOverMtMt.data[k] * mTotal * mTotal;
			s2Vec.data[k] = s2VecOverMtMt.data[k] * mTotal * mTotal;
    }

    seobParams.s1Vec = &s1VecOverMtMt;
    seobParams.s2Vec = &s2VecOverMtMt;

    XLALSimIMRSpinEOBCalculateSigmaStar( sigmaStar, m1, m2, &s1Vec, &s2Vec );
    XLALSimIMRSpinEOBCalculateSigmaKerr( sigmaKerr, m1, m2, &s1Vec, &s2Vec );

    seobParams.a = a = sqrt(inner_product(sigmaKerr->data, sigmaKerr->data));

    rcrossrdot[0] = LNhx;
    rcrossrdot[1] = LNhy;
    rcrossrdot[2] = LNhz;
    /* Eq. 16 of PRD 89, 084006 (2014): it's S_{1,2}/m_{1,2}^2.LNhat */
    s1dotLN = inner_product(s1Vec.data, rcrossrdot) / (m1*m1);
    s2dotLN = inner_product(s2Vec.data, rcrossrdot) / (m2*m2);

    chiS = 0.5 * (s1dotLN + s2dotLN);
    chiA = 0.5 * (s1dotLN - s2dotLN);

    switch ( SpinAlignedEOBversion )
    {
     case 1:
        /* See below Eq. 17 of PRD 86, 041011 (2012) */
       tplspin = 0.0;
       break;
     case 2:
       /* See below Eq. 4 of PRD 89, 061502(R) (2014)*/
       tplspin = (1.-2.*eta) * chiS + (m1 - m2)/(m1 + m2) * chiA;
       break;
     default:
       XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
       FREE_EVERYTHING
       FREE_SPHHARM
       XLAL_ERROR( XLAL_EINVAL );
       break;
    }

    if ( XLALSimIMRCalculateSpinPrecEOBHCoeffs( &seobCoeffs, eta, a,
                          SpinAlignedEOBversion ) == XLAL_FAILURE )
      XLAL_PRINT_INFO("\nSomething went wrong evaluating XLALSimIMRCalculateSpinPrecEOBHCoeffs in step %d of coarse dynamics\n",
			i );

    /* Update hlm coefficients */
    if ( XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients( eobParams.hCoeffs, m1, m2, eta,
                                                          tplspin, chiS, chiA, 3 ) == XLAL_FAILURE ) /* This function returns XLAL_SUCCESS or calls XLAL_ERROR( XLAL_EINVAL ) */
      {
      XLAL_PRINT_INFO("\nSomething went wrong evaluating XLALSimIMRCalculateSpinPrecEOBHCoeffs in step %d of coarse dynamics\n",
			i );
        XLALPrintError("XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients failed!\n");
        PRINT_PARAMS
        XLAL_ERROR( XLAL_EFUNC );
    }

    ham = XLALSimIMRSpinPrecEOBHamiltonian( eta, &cartPosVec, &cartMomVec,
                  &s1VecOverMtMt, &s2VecOverMtMt,
                  sigmaKerr, sigmaStar, seobParams.tortoise, &seobCoeffs );

    /* Calculate the polar data */
    polarDynamics.length = 4;
    polarDynamics.data   = polData;

    /* Calculate the orbital angular momentum */
    cross_product( values->data, values->data+3, rcrossp );
    magL = sqrt(inner_product(rcrossp, rcrossp));

    polarDynamics.data[0] = sqrt(inner_product(values->data, values->data));
    polarDynamics.data[1] = phiMod.data[i] + phiDMod.data[i];
    polarDynamics.data[2] = inner_product(values->data, values->data+3) / polarDynamics.data[0];
    polarDynamics.data[3] = magL;

    if ( XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform( &hLM,
          &polarDynamics, values, v, ham, 2, 2, &seobParams ) == XLAL_FAILURE )
    {
        FREE_EVERYTHING
        FREE_SPHHARM
        XLALPrintError("XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform failed!\n");
        PRINT_PARAMS
        XLAL_ERROR( XLAL_EDOM );
    }
    if ( XLALSimIMRSpinEOBNonQCCorrection( &hNQC,
          values, omega, &nqcCoeffs ) == XLAL_FAILURE )
    {
        FREE_EVERYTHING
        FREE_SPHHARM
        XLALPrintError("XLALSimIMRSpinEOBNonQCCorrection failed!\n");
        PRINT_PARAMS
        XLAL_ERROR( XLAL_EDOM );
    }
    hLM *= hNQC;
    h22TS->data->data[i]  = hLM;
    h2m2TS->data->data[i] = conj(hLM);

    if (debugPK) {
      fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e ",
             i*deltaT/mTScaled, aI2P, bI2P, gI2P, aP2J, bP2J, gP2J );
      fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
                      polData[0], polData[1], polData[2], polData[3], ham, v,
                      creal(hLM), cimag(hLM), creal(hLM/hNQC), cimag(hLM/hNQC) );
      fflush(NULL);
    }

    if ( XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform( &hLM,
        &polarDynamics, values, v, ham, 2, 1, &seobParams ) == XLAL_FAILURE )
    {
        FREE_EVERYTHING
        FREE_SPHHARM
        XLALPrintError("XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform failed!\n");
        PRINT_PARAMS
        XLAL_ERROR( XLAL_EDOM );
    }
    h21TS->data->data[i]  = hLM;
    h2m1TS->data->data[i] = conj(hLM);
    h20TS->data->data[i]  = 0.0;
  }
  if(debugPK) {
      fclose( out );
      //fclose( out2 );
      XLAL_PRINT_INFO("YP: quasi-nonprecessing modes generated.\n"); fflush(NULL);
  }

  /* Add quasi-nonprecessing spherical harmonic modes to the SphHarmTimeSeries structure */
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
    XLAL_PRINT_INFO("YP: SphHarmTS structures populated.\n"); fflush(NULL);
    out = fopen( "PWaves.dat", "w" );
    for ( i = 0; i < retLenLow; i++ )
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
    XLAL_PRINT_INFO("YP: P-frame waveforms written to file.\n");
    fflush(NULL);
  }


  /* Main loop for quasi-nonprecessing waveform generation -- HIGH SR */
    if ( !(alphaI2PTSHi = XLALCreateREAL8TimeSeries( "alphaI2P", &tc, 0.0, deltaT, &lalStrainUnit, retLenHi )) + !(betaI2PTSHi = XLALCreateREAL8TimeSeries(  "betaI2P", &tc, 0.0, deltaT, &lalStrainUnit, retLenHi )) + !(gammaI2PTSHi = XLALCreateREAL8TimeSeries( "gammaI2P", &tc, 0.0, deltaT, &lalStrainUnit, retLenHi )) + !(alphaP2JTSHi = XLALCreateREAL8TimeSeries( "alphaP2J", &tc, 0.0, deltaT, &lalStrainUnit, retLenHi )) + !(betaP2JTSHi = XLALCreateREAL8TimeSeries(  "betaP2J", &tc, 0.0, deltaT, &lalStrainUnit, retLenHi )) + !(gammaP2JTSHi = XLALCreateREAL8TimeSeries( "gammaP2J", &tc, 0.0, deltaT, &lalStrainUnit, retLenHi )) )
    {
        FREE_EVERYTHING
        XLALDestroyREAL8Vector( tlistHi );
        XLALDestroyREAL8Vector( timeJFull );
        XLALDestroyREAL8Vector( timeIFull );
        XLALDestroyREAL8Vector( tlistRDPatch );
        XLALDestroyREAL8Vector( tlistRDPatchHi );
        XLALPrintError("Failed to allocate alphaI2PTSHi or betaI2PTSHi or gammaI2PTSHi or alphaP2JTSHi or betaP2JTSHi or gammaP2JTSHi!\n");
        PRINT_PARAMS
        XLAL_ERROR(  XLAL_ENOMEM );
    }
    
    if ( !(h22TSHi   = XLALCreateCOMPLEX16TimeSeries( "H_22",  &tc, 0.0, deltaT, &lalStrainUnit, retLenHi )) + !(h21TSHi   = XLALCreateCOMPLEX16TimeSeries( "H_21", &tc, 0.0, deltaT, &lalStrainUnit, retLenHi )) + !(h20TSHi   = XLALCreateCOMPLEX16TimeSeries( "H_20",  &tc, 0.0, deltaT, &lalStrainUnit, retLenHi )) + !(h2m1TSHi  = XLALCreateCOMPLEX16TimeSeries( "H_2m1", &tc, 0.0, deltaT, &lalStrainUnit, retLenHi )) + !(h2m2TSHi  = XLALCreateCOMPLEX16TimeSeries( "H_2m2", &tc, 0.0, deltaT, &lalStrainUnit, retLenHi ))) {
        FREE_EVERYTHING
        XLALDestroyREAL8Vector( tlistHi );
        XLALDestroyREAL8Vector( timeJFull );
        XLALDestroyREAL8Vector( timeIFull );
        XLALDestroyREAL8Vector( tlistRDPatch );
        XLALDestroyREAL8Vector( tlistRDPatchHi );
        XLALPrintError("Failed to allocate h22TSHi or h21TSHi or h20TSHi or h2m1TSHi or h2m2TSHi!\n");
        PRINT_PARAMS
        XLAL_ERROR(  XLAL_ENOMEM );
    }

   hIMRJTSHi    = XLALCreateCOMPLEX16TimeSeries(
   "HIMRJHi",     &tc, 0.0, deltaTHigh, &lalStrainUnit, retLenHi + retLenRDPatchHi);

  if ( !(tlistHi = XLALCreateREAL8Vector( retLenHi )) || !(tlistRDPatchHi = XLALCreateREAL8Vector( retLenHi + retLenRDPatchHi)) )
  {
    XLALPrintError("Failed to allocate REAL8Vector tlistHi!\n");
    PRINT_PARAMS
    XLAL_ERROR(  XLAL_ENOMEM );
  }
  memset( tlistHi->data, 0, tlistHi->length * sizeof( REAL8 ));
  memset( tlistRDPatchHi->data, 0, tlistRDPatchHi->length * sizeof( REAL8 ));
  for ( i = 0; i < retLenHi; i++ )
  {
    tlistHi->data[i] = i * deltaTHigh/mTScaled + HiSRstart;
    tlistRDPatchHi->data[i] = i * deltaTHigh/mTScaled + HiSRstart;
  }
  for ( i = retLenHi; i < retLenHi + retLenRDPatchHi; i++ )
  {
    tlistRDPatchHi->data[i] = i * deltaTHigh/mTScaled + HiSRstart;
  }


  // Generating modes for finely sampled portion
  if(debugPK) {
    XLAL_PRINT_INFO("Generating precessing-frame modes for fine dynamics\n");
    fflush(NULL);
  }
  if (debugPK) out = fopen( "rotDynamicsHi.dat", "w" );
  for ( i = 0; i < retLenHi; i++ )
  {
    for ( j = 0; j < values->length; j++ )
    {
      values->data[j] = dynamicsHi->data[(j+1)*retLenHi + i];
    }

    /* Calculate dr/dt */
    memset( dvalues->data, 0, 14*sizeof(dvalues->data[0]));
    if( XLALSpinPrecHcapRvecDerivative( 0, values->data, dvalues->data,
          &seobParams) != XLAL_SUCCESS )
    {
        XLAL_PRINT_INFO( " Calculation of dr/dt at t = tPeakOmega \n");
        FREE_EVERYTHING
        XLALDestroyREAL8Vector( tlistHi );
        XLALDestroyREAL8Vector( timeJFull );
        XLALDestroyREAL8Vector( timeIFull );
        XLALDestroyREAL8Vector( tlistRDPatch );
        XLALDestroyREAL8Vector( tlistRDPatchHi );
        XLALPrintError("XLALSpinPrecHcapRvecDerivative failed!\n");
        PRINT_PARAMS
        XLAL_ERROR( XLAL_EDOM );
    }

    /* Calculate omega */
    rvec[0] = posVecxHi.data[i];
    rvec[1] = posVecyHi.data[i];
    rvec[2] = posVeczHi.data[i];
    memcpy( rdotvec, dvalues->data, 3*sizeof(REAL8));
    cross_product( rvec, rdotvec, rcrossrdot );

    magR   = sqrt(inner_product(rvec, rvec));
    omega  = sqrt(inner_product(rcrossrdot, rcrossrdot)) / (magR*magR);
    v = cbrt( omega );

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

    aI2P = AlphaHi->data[i];
    bI2P = BetaHi->data[i];
    gI2P = GammaHi->data[i];
    
    EulerAnglesP2J(&aP2J, &bP2J, &gP2J, aI2P, bI2P, gI2P, LNhx, LNhy, LNhz, JframeEx, JframeEy, JframeEz);

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
    s1Vec.data = s1Data;
    s2Vec.data = s2Data;

    s1VecOverMtMt.data = s1DataNorm;
    s2VecOverMtMt.data = s2DataNorm;
    for( k = 0; k < 3; k++ )
    {
        s1VecOverMtMt.data[k] = values->data[k+6];
        s2VecOverMtMt.data[k] = values->data[k+9];
        s1Vec.data[k] = s1VecOverMtMt.data[k] * mTotal * mTotal;
        s2Vec.data[k] = s2VecOverMtMt.data[k] * mTotal * mTotal;
    }

    seobParams.s1Vec = &s1VecOverMtMt;
    seobParams.s2Vec = &s2VecOverMtMt;

    XLALSimIMRSpinEOBCalculateSigmaStar( sigmaStar, m1, m2, &s1Vec, &s2Vec );
    XLALSimIMRSpinEOBCalculateSigmaKerr( sigmaKerr, m1, m2, &s1Vec, &s2Vec );

    seobParams.a = a = sqrt(inner_product(sigmaKerr->data, sigmaKerr->data));

    rcrossrdot[0] = LNhx;
    rcrossrdot[1] = LNhy;
    rcrossrdot[2] = LNhz;
    s1dotLN = inner_product(s1Vec.data, rcrossrdot) / (m1*m1);
    s2dotLN = inner_product(s2Vec.data, rcrossrdot) / (m2*m2);

    chiS = 0.5 * (s1dotLN + s2dotLN);
    chiA = 0.5 * (s1dotLN - s2dotLN);

    switch ( SpinAlignedEOBversion )
    {
     case 1:
       /* See below Eq. 17 of PRD 86, 041011 (2012) */
       tplspin = 0.0;
       break;
     case 2:
       /* See below Eq. 4 of PRD 89, 061502(R) (2014)*/
       tplspin = (1.-2.*eta) * chiS + (m1 - m2)/(m1 + m2) * chiA;
       break;
     default:
       XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
       FREE_EVERYTHING
        XLALDestroyREAL8Vector( tlistHi );
        XLALDestroyREAL8Vector( timeJFull );
        XLALDestroyREAL8Vector( timeIFull );
        XLALDestroyREAL8Vector( tlistRDPatch );
        XLALDestroyREAL8Vector( tlistRDPatchHi );
        PRINT_PARAMS
        XLAL_ERROR( XLAL_EINVAL );
       break;
    }

    if ( XLALSimIMRCalculateSpinPrecEOBHCoeffs( &seobCoeffs, eta, a,
                          SpinAlignedEOBversion ) == XLAL_FAILURE )
      XLAL_PRINT_INFO("\nSomething went wrong evaluating XLALSimIMRCalculateSpinPrecEOBHCoeffs in step %d of coarse dynamics\n",
			i );

    /* Update hlm coefficients */
    if ( XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients( eobParams.hCoeffs, m1, m2, eta,
                                                          tplspin, chiS, chiA, 3 ) == XLAL_FAILURE ) {
        XLAL_PRINT_INFO("\nSomething went wrong evaluating XLALSimIMRCalculateSpinPrecEOBHCoeffs in step %d of coarse dynamics\n",
			i );
        FREE_EVERYTHING
        XLALDestroyREAL8Vector( tlistHi );
        XLALDestroyREAL8Vector( timeJFull );
        XLALDestroyREAL8Vector( timeIFull );
        XLALDestroyREAL8Vector( tlistRDPatch );
        XLALDestroyREAL8Vector( tlistRDPatchHi );
        XLALPrintError("XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients failed!\n");
        PRINT_PARAMS
        XLAL_ERROR( XLAL_EDOM );
    }

    ham = XLALSimIMRSpinPrecEOBHamiltonian( eta, &cartPosVec, &cartMomVec,
                  &s1VecOverMtMt, &s2VecOverMtMt,
                  sigmaKerr, sigmaStar, seobParams.tortoise, &seobCoeffs );

    /* Calculate the polar data */
    polarDynamics.length = 4;
    polarDynamics.data   = polData;

    /* Calculate the orbital angular momentum */
    cross_product(values->data, values->data+3, rcrossp);

    polarDynamics.data[0] = sqrt(inner_product(values->data, values->data));
    polarDynamics.data[1] = phiModHi.data[i] + phiDModHi.data[i];
    polarDynamics.data[2] = inner_product(values->data, values->data+3) / polarDynamics.data[0];
    polarDynamics.data[3] = sqrt(inner_product(rcrossp, rcrossp));


    if ( XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform( &hLM, &polarDynamics, values, v, ham, 2, 2, &seobParams )
           == XLAL_FAILURE )
    {
        FREE_EVERYTHING
        XLALDestroyREAL8Vector( tlistHi );
        XLALDestroyREAL8Vector( timeJFull );
        XLALDestroyREAL8Vector( timeIFull );
        XLALDestroyREAL8Vector( tlistRDPatch );
        XLALDestroyREAL8Vector( tlistRDPatchHi );
        XLALPrintError("XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform failed!\n");
        PRINT_PARAMS
        XLAL_ERROR( XLAL_EDOM );
    }
    if ( XLALSimIMRSpinEOBNonQCCorrection( &hNQC, values, omega, &nqcCoeffs ) == XLAL_FAILURE )
    {
        FREE_EVERYTHING
        XLALDestroyREAL8Vector( tlistHi );
        XLALDestroyREAL8Vector( timeJFull );
        XLALDestroyREAL8Vector( timeIFull );
        XLALDestroyREAL8Vector( tlistRDPatch );
        XLALDestroyREAL8Vector( tlistRDPatchHi );
        XLALPrintError("XLALSimIMRSpinEOBNonQCCorrection failed!\n");
        PRINT_PARAMS
        XLAL_ERROR( XLAL_EDOM );
    }
    hLM *= hNQC;
    h22TSHi->data->data[i]  = hLM;
    h2m2TSHi->data->data[i] = conj(hLM);

    if (debugPK){
       fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e ",
             i*deltaTHigh/mTScaled + HiSRstart, aI2P, bI2P, gI2P, aP2J, bP2J, gP2J );
       fprintf( out, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e, %.16e\n",
             rdotvec[0], rdotvec[1], rdotvec[2], LNhx, LNhy, LNhz, creal(hLM), cimag(hLM), polData[0]);
    }

    if ( XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform( &hLM, &polarDynamics, values, v, ham, 2, 1, &seobParams )
           == XLAL_FAILURE )
    {
        FREE_EVERYTHING
        XLALDestroyREAL8Vector( tlistHi );
        XLALDestroyREAL8Vector( timeJFull );
        XLALDestroyREAL8Vector( timeIFull );
        XLALDestroyREAL8Vector( tlistRDPatch );
        XLALDestroyREAL8Vector( tlistRDPatchHi );
        XLALPrintError("XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform failed!\n");
        PRINT_PARAMS
        XLAL_ERROR( XLAL_EDOM );
    }
    h21TSHi->data->data[i]  = hLM;
    h2m1TSHi->data->data[i] = conj(hLM);
    h20TSHi->data->data[i]  = 0.0;
  }
  if (debugPK){
      fclose( out );
      XLAL_PRINT_INFO("YP: quasi-nonprecessing modes generated.for High SR\n");
  }

  /* Add quasi-nonprecessing spherical harmonic modes to the SphHarmTimeSeries structure */
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
  *hlmPTSoutput = hlmPTSout;

  h22PTSHi  = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, 2 );
  h21PTSHi  = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, 1 );
  h20PTSHi  = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, 0 );
  h2m1PTSHi = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, -1);
  h2m2PTSHi = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, -2);

  if (debugPK){
    XLAL_PRINT_INFO("YP: SphHarmTS structures populated for HIgh SR.\n");

    out = fopen( "PWavesHi.dat", "w" );
    for ( i = 0; i < retLenHi; i++ )
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
     XLAL_PRINT_INFO("YP: P-frame waveforms written to file for High SR.\n");
     fflush(NULL);
  }

  int foundAmp = 0;
    REAL8 tMaxAmp;
  double tAmpMax =  XLALSimLocateAmplTime(&timeHi, h22PTSHi->data, radiusVec, &foundAmp, &tMaxAmp);

  if(foundAmp==0){
      if (debugPK){
          XLAL_PRINT_INFO("maximum of Amplitude or dot{Ampl} are not found \n");
      }
  }
  else{
      if (debugPK){
          XLAL_PRINT_INFO("The Amplitude-related time is found and it's %f \n", tAmpMax);
      }
  }

  if( foundAmp==0 && foundPeakOmega == 0){
      if (debugPK){
          XLALPrintError("Houston, we've got a problem SOS, SOS, SOS, cannot find the RD attachment point... m1 = %.16e, m2 = %.16e, fMin = %.16e, inclination = %.16e,  Mtotal = %.16e, eta = %.16e, spin1 = {%.16e, %.16e, %.16e},   spin2 = {%.16e, %.16e, %.16e} \n", 
                   m1, m2, (double)fMin, (double)inc, mTotal, eta, spin1[0], spin1[1], spin1[2], spin2[0], spin2[1], spin2[2]);
      }
      
      tAttach = tAttach -2.0;    
      //FREE_EVERYTHING
      //XLALDestroyREAL8Vector( timeJFull );
      //XLALDestroyREAL8Vector( timeIFull );
      //XLALDestroyREAL8Vector( tlistRDPatch );
      //XLALDestroyREAL8Vector( tlistRDPatchHi );
      //XLAL_ERROR( XLAL_EINVAL );
  }


  if (tAmpMax < tAttach){
      if (debugPK){
          XLAL_PRINT_INFO("Stas the attachment time will be modified from %f to %f \n", tAttach, tAmpMax);
      }
      tAttach = tAmpMax;
  }

  if (debugPK){
      XLAL_PRINT_INFO("Stas: the final decision on the attachment time is %f \n", tAttach);
  }

  /* Changing retLen back to the Low SR value */
//  retLen = retLenLow;

/* *********************************************************************************
 * *********************************************************************************
 * STEP 6) Rotate quasi-nonprecessing waveforms from precessing to final J-frame
 * *********************************************************************************
 * **********************************************************************************/
  if ( XLALSimInspiralPrecessionRotateModes( hlmPTS, alphaP2JTS, betaP2JTS, gammaP2JTS ) == XLAL_FAILURE )
  {
      FREE_EVERYTHING
      XLALDestroyREAL8Vector( timeJFull );
      XLALDestroyREAL8Vector( timeIFull );
      XLALDestroyREAL8Vector( tlistRDPatch );
      XLALDestroyREAL8Vector( tlistRDPatchHi );
      XLALPrintError("XLALSimInspiralPrecessionRotateModes failed!\n");
      PRINT_PARAMS
      XLAL_ERROR( XLAL_EFAILED );
  }
  h22JTS  = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, 2 );
  h21JTS  = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, 1 );
  h20JTS  = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, 0 );
  h2m1JTS = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, -1);
  h2m2JTS = XLALSphHarmTimeSeriesGetMode( hlmPTS, 2, -2);


  if (debugPK){
    XLAL_PRINT_INFO("YP: PtoJ rotation done.\n");

    out = fopen( "JWaves.dat", "w" );
    for ( i = 0; i < retLenLow; i++ )
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
    XLAL_PRINT_INFO("YP: P-frame waveforms written to file.\n");
    fflush(NULL);
  }

  /* Stas: Rotating the high sampling part */
  if ( XLALSimInspiralPrecessionRotateModes( hlmPTSHi,
        alphaP2JTSHi, betaP2JTSHi, gammaP2JTSHi ) == XLAL_FAILURE )
  {
      FREE_EVERYTHING
      XLALDestroyREAL8Vector( timeJFull );
      XLALDestroyREAL8Vector( timeIFull );
      XLALDestroyREAL8Vector( tlistRDPatch );
      XLALDestroyREAL8Vector( tlistRDPatchHi );
      XLALPrintError("XLALSimInspiralPrecessionRotateModes failed!\n");
      PRINT_PARAMS
      XLAL_ERROR( XLAL_EFAILED );
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
        timeHi.data[i]+HiSRstart,
        creal(h22JTSHi->data->data[i]), cimag(h22JTSHi->data->data[i]),
        creal(h21JTSHi->data->data[i]), cimag(h21JTSHi->data->data[i]),
        creal(h20JTSHi->data->data[i]), cimag(h20JTSHi->data->data[i]),
        creal(h2m1JTSHi->data->data[i]), cimag(h2m1JTSHi->data->data[i]),
        creal(h2m2JTSHi->data->data[i]), cimag(h2m2JTSHi->data->data[i]) );
     }
     fclose( out );
  }

  sigReHi  = XLALCreateREAL8Vector( retLenHi + retLenRDPatchHi);
  sigImHi  = XLALCreateREAL8Vector( retLenHi + retLenRDPatchHi);
  if ( !sigReHi || !sigImHi )
  {
    XLALPrintError("Failed to allocate REAL8Vector sigReHi or sigImHi!\n");
    PRINT_PARAMS
    XLAL_ERROR( XLAL_ENOMEM );
  }
  memset( sigReHi->data, 0, sigReHi->length * sizeof( sigReHi->data[0] ));
  memset( sigImHi->data, 0, sigImHi->length * sizeof( sigImHi->data[0] ));


/* *********************************************************************************
 * *********************************************************************************
 * STEP 7) Attach ringdown to J-frame modes
 * *********************************************************************************
 * **********************************************************************************/
  if ( combSize > tAttach )
  {
      XLALPrintError( "Function: %s, The comb size looks to be too big!!!\n", __func__ );
      FREE_EVERYTHING
      XLALDestroyREAL8Vector( timeJFull );
      XLALDestroyREAL8Vector( timeIFull );
      XLALDestroyREAL8Vector( tlistRDPatch );
      XLALDestroyREAL8Vector( tlistRDPatchHi );
      PRINT_PARAMS
      XLAL_ERROR( XLAL_EINVAL );
  }
  rdMatchPoint->data[0] = combSize < tAttach ? tAttach - combSize : 0;
  rdMatchPoint->data[1] = tAttach;
  rdMatchPoint->data[2] = (retLenHi-1)*deltaTHigh/mTScaled;
  if (debugPK){
    XLAL_PRINT_INFO("YP::comb range: %f, %f\n",
        rdMatchPoint->data[0],rdMatchPoint->data[1]);
    XLAL_PRINT_INFO("Stas, tAttach = %f, comb range = %f to %f \n",
         tAttach, rdMatchPoint->data[0]+HiSRstart,
                  rdMatchPoint->data[1]+HiSRstart);
    fflush(NULL);
  }

  rdMatchPoint->data[0] = trunc(rdMatchPoint->data[0] / (deltaTHigh/mTScaled) ) * (deltaTHigh/mTScaled);
  rdMatchPoint->data[1] = trunc(rdMatchPoint->data[1] / (deltaTHigh/mTScaled) ) * (deltaTHigh/mTScaled);

  if (debugPK){
      XLAL_PRINT_INFO("Stas: matching points (again) %f, %f or %f, %f \n",
      rdMatchPoint->data[0],rdMatchPoint->data[1],
      rdMatchPoint->data[0]+HiSRstart,rdMatchPoint->data[1]+HiSRstart);
    fflush(NULL);
  }

    sh = 0.;
    REAL8 chi = tplspin/(1. - 2.*eta);
    if (chi >= 0.7 && chi < 0.8) {
        sh = -9. * (eta - 0.25);
    }
    if ((eta > 30. / 31. / 31. && eta <= 10. / 121. && chi >= 0.8) || (eta <= 30. / 31. / 31. && chi >= 0.8 && chi < 0.9)) {
        //This is case 4 in T1400476 - v3
        sh = -9. * (eta - 0.25) * (1. + 2. * exp(-(chi - 0.85) * (chi - 0.85) / 0.05 / 0.05)) * (1. + 1. / (1. + exp((eta - 0.01) / 0.001)));
    }
    if (eta < 30. / 31. / 31. && chi >= 0.9) {
        //This is case 5 in T1400476 - v3
        sh = 0.55 - 9. * (eta - 0.25) * (1. + 2. * exp(-(chi - 0.85) * (chi - 0.85) / 0.05 / 0.05)) * (1. + 1. / (1. + exp((eta - 0.01) / 0.001)));
    }
    if (eta > 10. / 121. && chi >= 0.8) {
        //This is case 6 in T1400476 - v3
        sh = 1. - 9. * (eta - 0.25) * (1. + 2. * exp(-(chi - 0.85) * (chi - 0.85) / 0.05 / 0.05)) * (1. + 1. / (1. + exp((eta - 0.01) / 0.001)));
    }


  /* Self-adjusting ringdown attachment. See https://dcc.ligo.org/T1500548 */
  int CheckRDAttachment = 1;
  if (CheckRDAttachment ){

      if (debugPK) XLAL_PRINT_INFO("Stas: checking the RD attachment point....\n");fflush(NULL);

      int pass = 0;
      rdMatchPoint->data[0] = combSize < tAttach ? tAttach - combSize : 0;
      rdMatchPoint->data[1] = tAttach;
      rdMatchPoint->data[0] = rdMatchPoint->data[0] - sh;
      rdMatchPoint->data[1] = rdMatchPoint->data[1] - sh;
      rdMatchPoint->data[2] = (retLenHi-1)*deltaTHigh/mTScaled;
      rdMatchPoint->data[0] = trunc(rdMatchPoint->data[0] / (deltaTHigh/mTScaled) ) * (deltaTHigh/mTScaled);
      rdMatchPoint->data[1] = trunc(rdMatchPoint->data[1] / (deltaTHigh/mTScaled) ) * (deltaTHigh/mTScaled);

      REAL8 thr = 1.01;
      REAL8 ratio22 = 1.0;
      REAL8 ratio2m2 = 1.0;
      REAL8 timediff = 0.;

      for ( i = 0; i < retLenHi; i++ )
        {
          sigReHi->data[i] = creal(h22JTSHi->data->data[i]);
          sigImHi->data[i] = cimag(h22JTSHi->data->data[i]);
        }

      if( XLALSimCheckRDattachment(sigReHi, sigImHi, &ratio22, tAttach, 2, 2,
                    deltaTHigh, m1, m2, 0.0, 0.0, chi1L, 0.0, 0.0, chi2L,
                    &timeHi, rdMatchPoint, spinEOBApproximant, kappaJL, &timediff) == XLAL_FAILURE )
      {
          FREE_EVERYTHING
          XLALDestroyREAL8Vector( timeJFull );
          XLALDestroyREAL8Vector( timeIFull );
          XLALDestroyREAL8Vector( tlistRDPatch );
          XLALDestroyREAL8Vector( tlistRDPatchHi );
          XLALPrintError("XLALSimCheckRDattachment failed!\n");
          PRINT_PARAMS
          XLAL_ERROR( XLAL_EDOM );
      }

      memset( sigReHi->data, 0, sigReHi->length * sizeof( sigReHi->data[0] ));
      memset( sigImHi->data, 0, sigImHi->length * sizeof( sigImHi->data[0] ));

      for ( i = 0; i < retLenHi; i++ )
        {
          sigReHi->data[i] = creal(h2m2JTSHi->data->data[i]);
          sigImHi->data[i] = cimag(h2m2JTSHi->data->data[i]);
        }

      if( XLALSimCheckRDattachment(sigReHi, sigImHi, &ratio2m2, tAttach, 2, -2,
                    deltaTHigh, m1, m2, 0.0, 0.0, chi1L, 0.0, 0.0, chi2L,
                    &timeHi, rdMatchPoint, spinEOBApproximant, kappaJL, &timediff ) == XLAL_FAILURE ) {
          FREE_EVERYTHING
          XLALDestroyREAL8Vector( timeJFull );
          XLALDestroyREAL8Vector( timeIFull );
          XLALDestroyREAL8Vector( tlistRDPatch );
          XLALDestroyREAL8Vector( tlistRDPatchHi );
          XLALPrintError("XLALSimCheckRDattachment failed!\n");
          PRINT_PARAMS
          XLAL_ERROR( XLAL_EDOM );
      }

      if (debugPK){
          XLAL_PRINT_INFO("At first check: ratio22 ratio2m2 %e %e\n", ratio22, ratio2m2);fflush(NULL);
      }
      if(ratio22 <= thr && ratio2m2 <=thr){
          pass = 1;
      }
      if (pass == 0){
          thr = 1.;
           if (debugPK){
               XLAL_PRINT_INFO("Adjusting RD attachment point... \n");fflush(NULL);
               XLAL_PRINT_INFO("initial ratios: %f,  %f \n", ratio22,  ratio2m2);
               fflush(NULL);
           }
               
           memset( sigReHi->data, 0, sigReHi->length * sizeof( sigReHi->data[0] ));
           memset( sigImHi->data, 0, sigImHi->length * sizeof( sigImHi->data[0] ));

           int found_att = XLALSimAdjustRDattachmentTime( sigReHi, sigImHi, h22JTSHi, h2m2JTSHi,
                    &ratio22, &ratio2m2, &tAttach, thr,
                    deltaTHigh, m1, m2, 0.0, 0.0, chi1L, 0.0, 0.0, chi2L,
                    &timeHi, rdMatchPoint, spinEOBApproximant, kappaJL, combSize, tMaxOmega, tMaxAmp);

           if (debugPK){
             if (found_att == 1){
                 XLAL_PRINT_INFO("we have found new attachment point tAtt = %f, with ratios %f, %f \n",
                       tAttach, ratio22, ratio2m2);fflush(NULL);
             }
               else if (found_att == 2) {
                   XLAL_PRINT_INFO("we haven't found proper attachment point, we picked the best point at tAtt = %f with ratios %f, %f\n", tAttach, ratio22, ratio2m2);fflush(NULL);

               }
               else{
                 XLAL_PRINT_INFO("we haven't found proper attachment point, best ratios are %f, %f at tAtt = %f \n",
                       ratio22, ratio2m2, tAttach);fflush(NULL);
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
  rdMatchPoint->data[0] = rdMatchPoint->data[0] - sh;
  rdMatchPoint->data[1] = rdMatchPoint->data[1] - sh;
  rdMatchPoint->data[2] = (retLenHi-1)*deltaTHigh/mTScaled;
  rdMatchPoint->data[0] = trunc(rdMatchPoint->data[0] / (deltaTHigh/mTScaled) ) * (deltaTHigh/mTScaled);
  rdMatchPoint->data[1] = trunc(rdMatchPoint->data[1] / (deltaTHigh/mTScaled) ) * (deltaTHigh/mTScaled);


  /*** Stas Let's try to attach RD to 2,2 mode: ***/
  INT4 m = 0;
  for ( m = -2; m < 3; m++ )
  {
    hJTSHi  = XLALSphHarmTimeSeriesGetMode( hlmPTSHi, 2, m );
    for ( i = 0; i < retLenHi; i++ )
    {
      sigReHi->data[i] = creal(hJTSHi->data->data[i]);
      sigImHi->data[i] = cimag(hJTSHi->data->data[i]);
    }
    if ( XLALSimIMREOBHybridAttachRingdownPrec( sigReHi, sigImHi, 2, m,
                deltaTHigh, m1, m2, 0.0, 0.0, chi1L, 0.0, 0.0, chi2L,
                &timeHi, rdMatchPoint, spinEOBApproximant, kappaJL )
            == XLAL_FAILURE ) {
        FREE_EVERYTHING
        XLALDestroyREAL8Vector( timeJFull );
        XLALDestroyREAL8Vector( timeIFull );
        XLALDestroyREAL8Vector( tlistRDPatch );
        XLALDestroyREAL8Vector( tlistRDPatchHi );
        XLALPrintError("XLALSimIMREOBHybridAttachRingdownPrec failed!\n");
        PRINT_PARAMS
        XLAL_ERROR( XLAL_EDOM );
    }

    for ( i = 0; i < (int)sigReHi->length; i++ )
    {
      hIMRJTSHi->data->data[i] = (sigReHi->data[i]) + I * (sigImHi->data[i]);
    }
    hIMRlmJTSHi = XLALSphHarmTimeSeriesAddMode( hIMRlmJTSHi, hIMRJTSHi, 2, m );
  }
  XLALSphHarmTimeSeriesSetTData( hIMRlmJTSHi, tlistRDPatchHi );
  if (debugPK){ XLAL_PRINT_INFO("YP: J wave RD attachment done.\n"); fflush(NULL); }

  hIMR22JTSHi  = XLALSphHarmTimeSeriesGetMode( hIMRlmJTSHi, 2, 2 );
  hIMR21JTSHi  = XLALSphHarmTimeSeriesGetMode( hIMRlmJTSHi, 2, 1 );
  hIMR20JTSHi  = XLALSphHarmTimeSeriesGetMode( hIMRlmJTSHi, 2, 0 );
  hIMR2m1JTSHi = XLALSphHarmTimeSeriesGetMode( hIMRlmJTSHi, 2, -1);
  hIMR2m2JTSHi = XLALSphHarmTimeSeriesGetMode( hIMRlmJTSHi, 2, -2);
 
  /** We search for maximum of the frame-invariant amplitude and set its time as time of
   * coalescence (or reference time) */ 
  // find tc which is accosiated with max of invariant amplitude:
  REAL8Vector *invAmp = XLALCreateREAL8Vector( hIMR22JTSHi->data->length ); 
  for (i=0; i < retLenHi + retLenRDPatchHi; i++ ){
      invAmp->data[i] =  creal(hIMR22JTSHi->data->data[i])*creal(hIMR22JTSHi->data->data[i]) + 
          cimag(hIMR22JTSHi->data->data[i]) * cimag(hIMR22JTSHi->data->data[i])  + 
           creal(hIMR21JTSHi->data->data[i])*creal(hIMR21JTSHi->data->data[i]) + 
          cimag(hIMR21JTSHi->data->data[i]) * cimag(hIMR21JTSHi->data->data[i])  + 
           creal(hIMR20JTSHi->data->data[i])*creal(hIMR20JTSHi->data->data[i]) + 
          cimag(hIMR20JTSHi->data->data[i]) * cimag(hIMR20JTSHi->data->data[i])  + 
           creal(hIMR2m1JTSHi->data->data[i])*creal(hIMR2m1JTSHi->data->data[i]) + 
          cimag(hIMR2m1JTSHi->data->data[i]) * cimag(hIMR2m1JTSHi->data->data[i])  + 
           creal(hIMR2m2JTSHi->data->data[i])*creal(hIMR2m2JTSHi->data->data[i]) + 
          cimag(hIMR2m2JTSHi->data->data[i]) * cimag(hIMR2m2JTSHi->data->data[i]);
  }
  REAL8 invAmpmax = invAmp->data[0];
  int i_maxiA = 0;
  for (i=1; i < retLenHi + retLenRDPatchHi - 1; i++){
      if ( invAmp->data[i] >= invAmp->data[i-1] && invAmp->data[i] > invAmp->data[i+1] && invAmp->data[i] > invAmpmax){
          i_maxiA = i;
          invAmpmax = invAmp->data[i];
      }
  } 
  if(debugPK){  
      XLAL_PRINT_INFO("We set tc = %.16e, %.16e \n", (tlistRDPatchHi->data[i_maxiA] ), -mTScaled * (tlistRDPatchHi->data[i_maxiA] ));
  }


  /* Having found the time of peak, we set the time of coalescence */
  //XLALGPSAdd(&tc, -mTScaled * (tPeakOmega + HiSRstart) );
  XLALGPSAdd( &tc, -mTScaled * (tlistRDPatchHi->data[i_maxiA] ));
  XLALDestroyREAL8Vector( invAmp ); // we don't need it anymore

  hPlusTS  = XLALCreateREAL8TimeSeries( "H_PLUS",
              &tc, 0.0, deltaT, &lalStrainUnit, retLenLow + retLenRDPatchLow );
  hCrossTS = XLALCreateREAL8TimeSeries( "H_CROSS",
              &tc, 0.0, deltaT, &lalStrainUnit, retLenLow + retLenRDPatchLow );



  *hIMRlmJTSHiOutput = hIMRlmJTSHi;
  if (debugPK){
    out = fopen( "JIMRWavesHi.dat", "w" );
    for ( i = 0; i < retLenHi + retLenRDPatchHi; i++ )
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

     XLAL_PRINT_INFO("Stas: down sampling and make the complete waveform in J-frame\n");
     fflush(NULL);
  }


  /*** attach hi sampling part and resample  */
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
         if (i*deltaT/mTScaled > HiSRstart) break;//Andrea
         hIMRJTS->data->data[i] = hJTS->data->data[i];
     }
      int idxRD = i; //Andrea
     spline = gsl_spline_alloc( gsl_interp_cspline, retLenHi + retLenRDPatchHi);
     acc    = gsl_interp_accel_alloc();
     gsl_spline_init( spline, tlistRDPatchHi->data, sigReHi->data, retLenHi + retLenRDPatchHi);
      for (i = idxRD; i< retLenLow+retLenRDPatchLow; i++){
          if( i*deltaT/mTScaled <= tlistRDPatchHi->data[ retLenHi + retLenRDPatchHi- 1]) {
              hIMRJTS->data->data[i] = gsl_spline_eval( spline, tlistRDPatch->data[i], acc );
          }
          else {
              hIMRJTS->data->data[i] = 0.;
          }
      }
     gsl_spline_free(spline);
     gsl_interp_accel_free(acc);

     spline = gsl_spline_alloc( gsl_interp_cspline, retLenHi + retLenRDPatchHi);
     acc    = gsl_interp_accel_alloc();
     gsl_spline_init( spline,
              tlistRDPatchHi->data, sigImHi->data, retLenHi + retLenRDPatchHi);
     for (i = idxRD; i< retLenLow+retLenRDPatchLow; i++){
         if( i*deltaT/mTScaled <= tlistRDPatchHi->data[ retLenHi + retLenRDPatchHi- 1]) {
             hIMRJTS->data->data[i] += I * gsl_spline_eval( spline, tlistRDPatch->data[i], acc );
         }
         else {
             hIMRJTS->data->data[i] += I*0.;
         }
     }
     gsl_spline_free(spline);
     gsl_interp_accel_free(acc);

     hIMRlmJTS = XLALSphHarmTimeSeriesAddMode( hIMRlmJTS, hIMRJTS, 2, k );
  }
  for (i=0; i<(int)tlistRDPatch->length; i++){
      timeJFull->data[i] = tlistRDPatch->data[i];    
  }
  XLALSphHarmTimeSeriesSetTData( hIMRlmJTS, timeJFull );
  if (debugPK){ XLAL_PRINT_INFO("Stas: J-wave with RD  generated.\n"); fflush(NULL); }

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
     XLAL_PRINT_INFO("YP: IMR J wave written to file.\n");
     fflush(NULL);
  }


/* *********************************************************************************
 * *********************************************************************************
 * STEP 8) Rotate modes from final final-J-frame to initial inertial frame
 * *********************************************************************************
 * **********************************************************************************/
 /* Euler Angles to go from final-J-frame to initial inertial frame  */
  gamJtoI = atan2(JframeEz[1], -JframeEz[0]);
  betJtoI = acos(JframeEz[2]);
  alJtoI = atan2(JframeEy[2], JframeEx[2]);
  if (fabs(betJtoI - LAL_PI) < 1.e-10){
      gamJtoI = 0.0;
      alJtoI = atan2(JframeEx[1], JframeEx[0]);
  }
  if (fabs(betJtoI) < 1.e-10){
        gamJtoI = 0.0;
        alJtoI = atan2(JframeEx[1], JframeEx[0]);
    }

  if (debugPK){
    XLAL_PRINT_INFO("Stas: J->I EA = %.16e, %.16e, %.16e \n", alJtoI, betJtoI, gamJtoI);
    fflush(NULL);
  }

   alpI        = XLALCreateREAL8TimeSeries( "alphaJ2I", &tc, 0.0, deltaT, &lalStrainUnit, retLenLow + retLenRDPatchLow );
   betI        = XLALCreateREAL8TimeSeries( "betaJ2I", &tc, 0.0, deltaT, &lalStrainUnit, retLenLow + retLenRDPatchLow );
   gamI        = XLALCreateREAL8TimeSeries( "gammaJ2I", &tc, 0.0, deltaT, &lalStrainUnit, retLenLow + retLenRDPatchLow );

  /** NOTE we are making again transformation between notations of Arun et. al adopted here and 
   * Wiegner matrix (active rotation) coded up in LAL */
  for (i=0; i< retLenLow + retLenRDPatchLow; i++){
      alpI->data->data[i] = -alJtoI;
      betI->data->data[i] = betJtoI;
      gamI->data->data[i] = -gamJtoI;
  }
  for (i=0; i<(int)tlistRDPatch->length; i++){
      timeIFull->data[i] = tlistRDPatch->data[i];    
  }

  hIMRlmITS = XLALSphHarmTimeSeriesAddMode( hIMRlmITS, hIMR22JTS, 2, 2 );
  hIMRlmITS = XLALSphHarmTimeSeriesAddMode( hIMRlmITS, hIMR21JTS, 2, 1 );
  hIMRlmITS = XLALSphHarmTimeSeriesAddMode( hIMRlmITS, hIMR20JTS, 2, 0 );
  hIMRlmITS = XLALSphHarmTimeSeriesAddMode( hIMRlmITS, hIMR2m1JTS, 2, -1 );
  hIMRlmITS = XLALSphHarmTimeSeriesAddMode( hIMRlmITS, hIMR2m2JTS, 2, -2 );
  XLALSphHarmTimeSeriesSetTData( hIMRlmITS, timeIFull );

  /* Rotate J-frame modes */
  if (debugPK){ XLAL_PRINT_INFO("Rotation to inertial I-frame\n"); fflush(NULL); }
  if ( XLALSimInspiralPrecessionRotateModes( hIMRlmITS,
                alpI, betI, gamI ) == XLAL_FAILURE ) {
      FREE_EVERYTHING
      XLALPrintError("XLALSimInspiralPrecessionRotateModes failed!\n");
      PRINT_PARAMS
      XLAL_ERROR( XLAL_EFAILED );
  }
  hIMR22ITS  = XLALSphHarmTimeSeriesGetMode( hIMRlmITS, 2, 2 );
  hIMR21ITS  = XLALSphHarmTimeSeriesGetMode( hIMRlmITS, 2, 1 );
  hIMR20ITS  = XLALSphHarmTimeSeriesGetMode( hIMRlmITS, 2, 0 );
  hIMR2m1ITS = XLALSphHarmTimeSeriesGetMode( hIMRlmITS, 2, -1);
  hIMR2m2ITS = XLALSphHarmTimeSeriesGetMode( hIMRlmITS, 2, -2);

  *hIMRoutput = XLALSphHarmTimeSeriesAddMode( *hIMRoutput, hIMR22ITS, 2, 2 );
  *hIMRoutput = XLALSphHarmTimeSeriesAddMode( *hIMRoutput, hIMR21ITS, 2, 1 );
  *hIMRoutput = XLALSphHarmTimeSeriesAddMode( *hIMRoutput, hIMR20ITS, 2, 0 );
  *hIMRoutput = XLALSphHarmTimeSeriesAddMode( *hIMRoutput, hIMR2m1ITS, 2, -1 );
  *hIMRoutput = XLALSphHarmTimeSeriesAddMode( *hIMRoutput, hIMR2m2ITS, 2, -2 );
  XLALSphHarmTimeSeriesSetTData( *hIMRoutput, tlistRDPatch );

  if (debugPK){
    out = fopen( "IWaves.dat", "w" );
    for ( i = 0; i < retLenLow + retLenRDPatchLow; i++ )
    {
      fprintf( out,
        "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
        tlistRDPatch->data[i],
        creal(hIMR22ITS->data->data[i]), cimag(hIMR22ITS->data->data[i]),
        creal(hIMR21ITS->data->data[i]), cimag(hIMR21ITS->data->data[i]),
        creal(hIMR20ITS->data->data[i]), cimag(hIMR20ITS->data->data[i]),
        creal(hIMR2m1ITS->data->data[i]), cimag(hIMR2m1ITS->data->data[i]),
        creal(hIMR2m2ITS->data->data[i]), cimag(hIMR2m2ITS->data->data[i]) );
    }
    fclose( out );
  }

/* *********************************************************************************
 * *********************************************************************************
 * STEP 9) Compute h+, hx
 * ********************************************************************************
 * **********************************************************************************/
    //incl_temp = inc;
    //incl_temp = 0.0;
    /** NOTE that we have use  now different convention: the phi_ref (or one might call it phi_c) is now
     * the azimuthal phase of the observer in the source (I-)frame. Together with inclination angle it defines 
     * the position of the observer in the (I-)frame associated with the source at t=0 */

    Y22 = XLALSpinWeightedSphericalHarmonic( inc, -phiC, -2, 2, 2 );
    Y2m2 = XLALSpinWeightedSphericalHarmonic( inc, -phiC, -2, 2, -2 );
    Y21 = XLALSpinWeightedSphericalHarmonic( inc, -phiC, -2, 2, 1 );
    Y2m1 = XLALSpinWeightedSphericalHarmonic( inc, -phiC, -2, 2, -1 );
    Y20 = XLALSpinWeightedSphericalHarmonic( inc, -phiC, -2, 2, 0 );
    /*if ( SpinsAlmostAligned ) {
        Y22 = XLALSpinWeightedSphericalHarmonic( inc, coa_phase_offset, -2, 2, 2 );
        Y2m2 = XLALSpinWeightedSphericalHarmonic( inc, coa_phase_offset, -2, 2, -2 );
        Y21 = XLALSpinWeightedSphericalHarmonic( inc, coa_phase_offset, -2, 2, 1 );
        Y2m1 = XLALSpinWeightedSphericalHarmonic( inc, coa_phase_offset, -2, 2, -1 );
        Y20 = XLALSpinWeightedSphericalHarmonic( inc, coa_phase_offset, -2, 2, 0 );
    }
    else {
        Y22 = XLALSpinWeightedSphericalHarmonic( incl_temp, coa_phase_offset, -2, 2, 2 );
        Y2m2 = XLALSpinWeightedSphericalHarmonic( incl_temp, coa_phase_offset, -2, 2, -2 );
        Y21 = XLALSpinWeightedSphericalHarmonic( incl_temp, coa_phase_offset, -2, 2, 1 );
        Y2m1 = XLALSpinWeightedSphericalHarmonic( incl_temp, coa_phase_offset, -2, 2, -1 );
        Y20 = XLALSpinWeightedSphericalHarmonic( incl_temp, coa_phase_offset, -2, 2, 0 );
    }*/

    if(debugPK){ XLAL_PRINT_INFO("Ylm %e %e %e %e %e %e %e %e %e %e \n", creal(Y22), cimag(Y22), creal(Y2m2), cimag(Y2m2), creal(Y21), cimag(Y21), creal(Y2m1), cimag(Y2m1), creal(Y20), cimag (Y20)); fflush(NULL); }

  for ( i = 0; i < (INT4)hIMR22ITS->data->length; i++ )
  {
    x11 = Y22*hIMR22ITS->data->data[i] + Y21*hIMR21ITS->data->data[i]
          + Y20*hIMR20ITS->data->data[i] + Y2m1*hIMR2m1ITS->data->data[i]
          + Y2m2*hIMR2m2ITS->data->data[i];
    hPlusTS->data->data[i]  = amp0*creal(x11);
    hCrossTS->data->data[i] = -amp0*cimag(x11);
  }

  /* Point the output pointers to the relevant time series and return */
  (*hplus)  = hPlusTS;
  (*hcross) = hCrossTS;

  if (debugPK){
    XLAL_PRINT_INFO("plus and cross are computed, freeing the memory\n");
    fflush(NULL);
  }
    XLALDestroyREAL8TimeSeries(alphaI2PTS);
    XLALDestroyREAL8TimeSeries(betaI2PTS);
    XLALDestroyREAL8TimeSeries(gammaI2PTS);
    XLALDestroyREAL8TimeSeries(alphaP2JTS);
    XLALDestroyREAL8TimeSeries(betaP2JTS);
    XLALDestroyREAL8TimeSeries(gammaP2JTS);
    XLALDestroyREAL8TimeSeries(alphaI2PTSHi);
    XLALDestroyREAL8TimeSeries(betaI2PTSHi);
    XLALDestroyREAL8TimeSeries(gammaI2PTSHi);
    XLALDestroyREAL8TimeSeries(alphaP2JTSHi);
    XLALDestroyREAL8TimeSeries(betaP2JTSHi);
    XLALDestroyREAL8TimeSeries(gammaP2JTSHi);
    XLALDestroyCOMPLEX16TimeSeries(h22TS);
    XLALDestroyCOMPLEX16TimeSeries(h21TS);
    XLALDestroyCOMPLEX16TimeSeries(h20TS);
    XLALDestroyCOMPLEX16TimeSeries(h2m1TS);
    XLALDestroyCOMPLEX16TimeSeries(h2m2TS);
    XLALDestroySphHarmTimeSeries(hlmPTS);
    XLALDestroySphHarmTimeSeries(hIMRlmJTS);
    XLALDestroySphHarmTimeSeries(hIMRlmITS);
    XLALDestroyCOMPLEX16TimeSeries(hIMRJTS);
    XLALDestroyCOMPLEX16TimeSeries(h22TSHi);
    XLALDestroyCOMPLEX16TimeSeries(h21TSHi);
    XLALDestroyCOMPLEX16TimeSeries(h20TSHi);
    XLALDestroyCOMPLEX16TimeSeries(h2m1TSHi);
    XLALDestroyCOMPLEX16TimeSeries(h2m2TSHi);
    XLALDestroyCOMPLEX16TimeSeries(hIMRJTSHi);
    XLALDestroyREAL8Vector( values );
    XLALDestroyREAL8Vector( dvalues );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigReHi );
    XLALDestroyREAL8Vector( sigImHi );
    XLALAdaptiveRungeKutta4Free(integrator);
    XLALDestroyREAL8Array( dynamics );
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
    XLALDestroyREAL8Vector( radiusVec  );
    
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
    XLALDestroyREAL8Vector( valuesV2 );
    if (dynamicsV2 != NULL){
        XLALDestroyREAL8Array( dynamicsV2 );
    }
    if (dynamicsV2Hi != NULL){
        XLALDestroyREAL8Array( dynamicsV2Hi );
    }
  if(debugPK){ XLAL_PRINT_INFO("Memory cleanup ALL done.\n"); fflush(NULL); }

  return XLAL_SUCCESS;
}
#undef FREE_EVERYTHING
#endif
