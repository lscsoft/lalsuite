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


#include <math.h>
#include <complex.h>
#include <lal/LALSimIMR.h>
#include <lal/LALSimInspiral.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/LALAdaptiveRungeKutta4.h>
#include <lal/SphericalHarmonics.h>
#include <gsl/gsl_sf_gamma.h>

#include "LALSimIMREOBNRv2.h"
#include "LALSimIMRSpinEOB.h"

/* Include all the static function files we need */
#include "LALSimIMREOBHybridRingdown.c"
#include "LALSimIMREOBFactorizedWaveform.c"
#include "LALSimIMREOBNewtonianMultipole.c"
#include "LALSimIMREOBNQCCorrection.c"
#include "LALSimIMRSpinEOBInitialConditions.c"
#include "LALSimIMRSpinEOBAuxFuncs.c"
#include "LALSimIMRSpinAlignedEOBHcapDerivative.c"
#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBFactorizedWaveform.c"
#include "LALSimIMRSpinEOBFactorizedFlux.c"

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
  int status;
  INT4 SpinAlignedEOBversion = 1;

  REAL8Vector *values = NULL;

  /* Parameters of the system */
  REAL8 mTotal, eta, mTScaled;
  REAL8 amp0, amp;
  REAL8 a, tplspin;
  REAL8 chiS, chiA;
  REAL8Vector *sigmaStar = NULL;
  REAL8Vector *sigmaKerr = NULL;

  /* Spins not scaled by the mass */
  REAL8 mSpin1[3], mSpin2[3];

  /* Parameter structures containing important parameters for the model */
  SpinEOBParams           seobParams;
  SpinEOBHCoeffs          seobCoeffs;
  EOBParams               eobParams;
  FacWaveformCoeffs       hCoeffs;
  NewtonMultipolePrefixes prefixes;

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

  /* Calculate chiS and chiA */
  /* Assuming we are in the minimally-rotating frame, such that the orbital
   * angular momentum is along the z-axis at the initial time. */
  chiS = 0.5 * (spin1[2] + spin2[2]);
  chiA = 0.5 * (spin1[2] - spin2[2]);

  /* Wrapper spin vectors used to calculate sigmas */
  REAL8Vector s1Vec, s1VecOverMtMt;
  REAL8Vector s2Vec, s2VecOverMtMt;
  REAL8       s1Data[3], s2Data[3], s1DataNorm[3], s2DataNorm[3];
  REAL8       omega, v, ham;

  s1VecOverMtMt.data   = s1DataNorm;
  s2VecOverMtMt.data   = s2DataNorm;

  memcpy( s1Data, spin1, sizeof(s1Data) );
  memcpy( s2Data, spin2, sizeof(s2Data) );
  memcpy( s1DataNorm, spin1, sizeof( s1DataNorm ) );
  memcpy( s2DataNorm, spin2, sizeof( s2DataNorm ) );

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

  s1Vec.length = s2Vec.length = 3;
  s1Vec.data   = s1Data;
  s2Vec.data   = s2Data;

  /* Populate the initial structures */
  if ( XLALSimIMRSpinEOBCalculateSigmaStar( sigmaStar, m1, m2,
                              &s1Vec, &s2Vec ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  if ( XLALSimIMRSpinEOBCalculateSigmaKerr( sigmaKerr, m1, m2,
                              &s1Vec, &s2Vec ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  seobParams.a = a = sigmaKerr->data[2];
  s1VecOverMtMt.length = s2VecOverMtMt.length = 3;
  seobParams.s1Vec    = &s1VecOverMtMt;
  seobParams.s2Vec    = &s2VecOverMtMt;

  /* Variables for the integrator */
  ark4GSLIntegrator       *integrator = NULL;
  REAL8Array              *dynamics   = NULL;
  INT4                    retLen;
  /*REAL8  UNUSED           tMax;*/

  /* Accuracies of adaptive Runge-Kutta integrator */
  const REAL8 EPS_ABS = 1.0e-9;
  const REAL8 EPS_REL = 1.0e-8;

  /* Initialize parameters */
  mTotal = m1 + m2;
  mTScaled = mTotal * LAL_MTSUN_SI;
  eta    = m1 * m2 / (mTotal*mTotal);

  amp0 = 4. * mTotal * LAL_MRSUN_SI * eta / r;

  /* Cartesian vectors needed to calculate Hamiltonian */
  REAL8Vector cartPosVec, cartMomVec;
  REAL8       cartPosData[3], cartMomData[3];

  cartPosVec.length = cartMomVec.length = 3;
  cartPosVec.data = cartPosData;
  cartMomVec.data = cartMomData;
  memset( cartPosData, 0, sizeof( cartPosData ) );
  memset( cartMomData, 0, sizeof( cartMomData ) );

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

  seobParams.tortoise     = 1;
  seobParams.sigmaStar    = sigmaStar;
  seobParams.sigmaKerr    = sigmaKerr;
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

  /* Populate the initial structures */
  if ( XLALSimIMRCalculateSpinEOBHCoeffs( &seobCoeffs, eta, a,
                          SpinAlignedEOBversion ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  if ( XLALSimIMREOBCalcPrecNoSpinFacWaveformCoefficients( &hCoeffs, m1, m2, eta,
        tplspin, chiS, chiA, SpinAlignedEOBversion ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( sigmaKerr );
    XLALDestroyREAL8Vector( sigmaStar );
    XLALDestroyREAL8Vector( values );
    XLAL_ERROR( XLAL_EFUNC );
  }

  if ( XLALSimIMREOBCalcPrecSpinFacWaveformCoefficients( &hCoeffs, m1, m2, eta,
        tplspin, chiS, chiA, SpinAlignedEOBversion ) == XLAL_FAILURE )
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

  /* TODO: Set the initial conditions */

  if ( XLALSimIMRSpinEOBInitialConditions( values, m1, m2, fMin, inc, mSpin1, mSpin2, &seobParams ) == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }
  //exit(0);
  //YP::{x,y,z,px,py,pz,s1x,s1y,s1z,s2x,s2y,s2z} =
  //{15.87, 0, 0, -0.000521675194648, 0.278174373488, -0.00012666165246, -0.270452950188, -0.216802131414, 0.00133043857763, 0, 0, 0};
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
#if 1
  values->data[0] = 15.87;
  values->data[1] = 0.;
  values->data[2] = 0.;
  values->data[3] = -0.000521675194648;
  values->data[4] = 0.278174373488;
  values->data[5] = -0.00012666165246;
  values->data[6] = -0.270452950188;
  values->data[7] = -0.216802131414;
  values->data[8] = 0.00133043857763;
  values->data[9] = 0.;
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
  REAL8 vX, vY, vZ, rCrossV_x, rCrossV_y, rCrossV_z, vOmega;
  REAL8 magPosVec, LNhx, LNhy, LNhz, magL, alpha;

  COMPLEX16   hLM;
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

  /* Cartesian vectors needed to calculate Hamiltonian */
  cartPosVec.length = cartMomVec.length = 3;
  cartPosVec.data = cartPosData;
  cartMomVec.data = cartMomData;
  memset( cartPosData, 0, sizeof( cartPosData ) );
  memset( cartMomData, 0, sizeof( cartMomData ) );

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

    printf( "alpha = %.16e, omega = %.16e, LNhz = %.16e, vphi = %.16e\n",
             alpha, omega, LNhz, vphi[i] );

    /* Calculate the value of the Hamiltonian */
    cartPosVec.data[0] = values->data[0];
    cartMomVec.data[0] = values->data[2];
    cartMomVec.data[1] = values->data[3] / values->data[0];

    omega = XLALSimIMRSpinAlignedEOBCalcOmega( values->data, &seobParams );
    v = cbrt( omega );

    ham = XLALSimIMRSpinEOBHamiltonian( eta, &cartPosVec, &cartMomVec,
                  &s1VecOverMtMt, &s2VecOverMtMt,
                  sigmaKerr, sigmaStar, seobParams.tortoise, &seobCoeffs );

    status = XLALSimIMRSpinEOBGetSpinFactorizedWaveform( &hLM, values, v,
                  ham, 2, 2, &seobParams );

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
