/*
*  Copyright (C) 2007 Stas Babak, David Churches, Duncan Brown, David Chin, Jolien Creighton,
*                     B.S. Sathyaprakash, Craig Robinson , Thomas Cokelaer, Evan Ochsner
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

\author Craig Robinson

\file

\brief Functions to generate the EOBNRv2 waveforms, as defined in
Pan et al, PRD84, 124052(2011).

*/

#include <complex.h>
#include <lal/Units.h>
#include <lal/LALAdaptiveRungeKutta4.h>
#include <lal/FindRoot.h>
#include <lal/SeqFactories.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <gsl/gsl_sf_gamma.h>
#include "LALSimIMREOBNRv2.h"
#include "LALSimBlackHoleRingdown.h"

/* Include all the static function files we need */
#include "LALSimIMREOBFactorizedWaveform.c"
#include "LALSimIMREOBFactorizedFlux.c"
#include "LALSimIMREOBNQCCorrection.c"
#include "LALSimIMREOBNewtonianMultipole.c"
#include "LALSimIMREOBHybridRingdown.c"
#include "LALSimInspiraldEnergyFlux.c"

/** 
 * The maximum number of modes available to us in this model
 */
static const int EOBNRV2_NUM_MODES_MAX = 5;

/**
 *  The following declarations are so that the compiler quits. It appears that they are here because
 *  some of these functions are called before their definitions in the following code. 
 *  Perhaps this can be made cleaner by just moving the definitions, but for now I'm keeping with what
 *  the original code did.
**/

static int
XLALSimIMREOBNRv2SetupFlux(
             expnCoeffsdEnergyFlux *ak, /**<<PN expansion coefficients (only relevant fields will be populated)*/
             REAL8                 eta  /**<< Symmetric mass ratio */
             );

static REAL8
omegaofrP4PN (
             const REAL8 r,
             const REAL8 eta,
             void *params);

static
int LALHCapDerivativesP4PN(   double t,
                               const double values[],
                               double dvalues[],
                               void   *funcParams);

static
REAL8 XLALCalculateOmega (   REAL8 eta,
                             REAL8 r,
                             REAL8 pr,
                             REAL8 pPhi,
                             EOBACoefficients *aCoeffs );

static
int XLALFirstStoppingCondition(double t,
                               const double values[],
                               double dvalues[],
                               void   *funcParams);

static
int XLALHighSRStoppingCondition(double t,
                                const double values[],
                                double dvalues[],
                                void   *funcParams);

static
REAL8 XLALprInitP4PN(REAL8 p, void  *params);

static
REAL8 XLALpphiInitP4PN( const REAL8 r,
                       EOBACoefficients * restrict coeffs );

static
REAL8 XLALrOfOmegaP4PN( REAL8 r, void *params);

static
REAL8 XLALvrP4PN(const REAL8 r, const REAL8 omega, pr3In *params);

static
size_t find_instant_freq(const REAL8TimeSeries *hp, const REAL8TimeSeries *hc, const REAL8 target, const size_t start, const int fsign);

static
size_t find_instant_freq_hlm(const COMPLEX16TimeSeries *hlm,
    const REAL8 target, const size_t start);

/*-------------------------------------------------------------------*/
/*                      pseudo-4PN functions                         */
/*-------------------------------------------------------------------*/

/**
 * Calculates the initial orbital momentum.
 */
REAL8
XLALpphiInitP4PN(
            const REAL8 r,                     /**<< Initial orbital separation */
            EOBACoefficients * restrict coeffs /**<< Pre-computed EOB A coefficients */
            )
{
  REAL8 u, u2, A, dA;

  u  = 1. / r;
  u2 = u*u;

  /* Comparing this expression with the old one, I *think* I will need to multiply */
  /* dA by - r^2. TODO: Check that this should be the case */

  A  = XLALCalculateEOBA( r, coeffs );
  dA = - r*r * XLALCalculateEOBdAdr( r, coeffs );
  return sqrt(-dA/(2.*u*A + u2 * dA));
/* why is it called phase? This is initial j!? */
}


/**
 * The name of this function is slightly misleading, as it does
 * not computs the initial pr directly. It instead calculated
 * dH/dpr - vr for the given values of vr and p. This function is
 * then passed to a root finding function, which determines the 
 * value of pr which makes this function zero. That value is then
 * the initial pr.
 */
REAL8
XLALprInitP4PN(
             REAL8 p,     /**<< The pr value we are currently testing */
             void *params /**<< The pr3In structure containing necessary parameters */
             )
{
  REAL8  pr;
  REAL8  u, u2, p2, p3, p4, A;
  REAL8  onebyD, AbyD, Heff, HReal, etahH;
  REAL8 eta, z3, r, vr, q;
  pr3In *ak;

  /* TODO: Change this to use tortoise coord */

  ak = (pr3In *) params;

  eta = ak->eta;
  vr = ak->vr;
  r = ak->r;
  q = ak->q;

   p2 = p*p;
   p3 = p2*p;
   p4 = p2*p2;
   u = 1./ r;
   u2 = u*u;
   z3 = 2. * (4. - 3. * eta) * eta;

   A = XLALCalculateEOBA( r, ak->aCoeffs );
   onebyD = 1. / XLALCalculateEOBD( r, eta );
   AbyD = A * onebyD;

   Heff = sqrt(A*(1. + AbyD * p2 + q*q * u2 + z3 * p4 * u2));
   HReal = sqrt(1. + 2.*eta*(Heff - 1.)) / eta;
   etahH = eta*Heff*HReal;

/* This sets pr = dH/dpr - vr, calls rootfinder,
   gets value of pr s.t. dH/pr = vr */
   pr = -vr +  A*(AbyD*p + 2. * z3 * u2 * p3)/etahH;

   return pr;
}


/*-------------------------------------------------------------------*/

/**
 * This function calculates the initial omega for a given value of r.
 */
static REAL8
omegaofrP4PN (
             const REAL8 r,   /**<< Initial orbital separation (in units of total mass M */
             const REAL8 eta, /**<< Symmetric mass ratio */
             void *params     /**<< The pre-computed A coefficients */
             )
{
   REAL8 u, u3, A, dA, omega;

   EOBACoefficients *coeffs;
   coeffs = (EOBACoefficients *) params;

   u = 1./r;
   u3 = u*u*u;

   A  = XLALCalculateEOBA( r, coeffs );
   dA = XLALCalculateEOBdAdr( r, coeffs );

   /* Comparing this expression with the old one, I *think* I will need to multiply */
   /* dA by - r^2. TODO: Check that this should be the case */
   dA = - r * r * dA;

   omega = sqrt(u3) * sqrt ( -0.5 * dA /(1. + 2.*eta * (A/sqrt(A+0.5 * u*dA)-1.)));
   return omega;
}

/*-------------------------------------------------------------------*/

/**
 * Function called within a root-finding algorithm used to determine
 * the initial radius for a given value of omega.
 */
REAL8
XLALrOfOmegaP4PN(
            REAL8 r,     /**<< Test value of the initial radius */
            void *params /**<< pr3In structure, containing useful parameters, including omega */
            )
{
  REAL8  omega1,omega2;
  pr3In *pr3in;

#ifndef LAL_NDEBUG
  if ( !params )
    XLAL_ERROR_REAL8( XLAL_EFAULT );
#endif

  pr3in = (pr3In *) params;

  omega1 = pr3in->omega;
  omega2 = omegaofrP4PN(r, pr3in->eta, pr3in->aCoeffs);
  return ( -omega1 + omega2 );
}

/*-------------------------------------------------------------------*/

/**
 * This function computes the derivatives of the EOB Hamiltonian w.r.t. the dynamical
 * variables, and therefore the derivatives of the dynamical variables w.r.t. time.
 * As such this gets called in the Runge-Kutta integration of the orbit.
 */
int
LALHCapDerivativesP4PN( double UNUSED t,        /**<< Current time (GSL requires it to be a parameter, but it's irrelevant) */
                        const REAL8 values[],   /**<< The dynamics r, phi, pr, pphi */
                        REAL8       dvalues[],  /**<< The derivatives dr/dt, dphi/dt. dpr/dt and dpphi/dt */
                        void        *funcParams /**<< Structure containing all the necessary parameters */
                      )
{

  EOBParams *params = NULL;

  /* Max l to sum up to in the factorized flux */
  const INT4 lMax = 8;

  REAL8 eta, omega;

  double r, p, q;
  REAL8 r2, p2, p4, q2;
  REAL8 u, u2, u3;
  REAL8 A, AoverSqrtD, dAdr, Heff, Hreal;
  REAL8 HeffHreal;
  REAL8 flux;
  REAL8 z3;

  /* Factorized flux function takes a REAL8Vector */
  /* We will wrap our current array for now */
  REAL8Vector valuesVec;

  params = (EOBParams *) funcParams;

  valuesVec.length = 4;
  memcpy( &(valuesVec.data), &(values), sizeof(REAL8 *) );

  eta = params->eta;

  z3   = 2. * ( 4. - 3. * eta ) * eta;

  r = values[0];
  p = values[2];
  q = values[3];

  u  = 1.0 / r;
  u2 = u * u;
  u3 = u2 * u;
  r2 = r * r;
  p2 = p*p;
  p4 = p2 * p2;
  q2 = q * q;

  A          = XLALCalculateEOBA(r, params->aCoeffs );
  dAdr       = XLALCalculateEOBdAdr(r, params->aCoeffs );
  AoverSqrtD = A / sqrt( XLALCalculateEOBD(r, eta) );

  /* Note that Hreal as given here is missing a factor of 1/eta */
  /* This is because it only enters into the derivatives in     */
  /* the combination eta*Hreal*Heff, so the eta would get       */
  /* cancelled out anyway. */

  Heff  = XLALEffectiveHamiltonian( eta, r, p, q, params->aCoeffs );
  Hreal = sqrt( 1. + 2.*eta*(Heff - 1.) );

  HeffHreal = Heff * Hreal;

  /* rDot */
  dvalues[0] = AoverSqrtD * u2 * p * (r2 + 2. * p2 * z3 * A ) / HeffHreal;

  /* sDot */
  dvalues[1] = omega = q * A * u2 / HeffHreal;

  /* Note that the only field of dvalues used in the flux is dvalues->data[1] */
  /* which we have already calculated. */
  flux = XLALSimIMREOBFactorizedFlux( &valuesVec, omega, params, lMax );

  /* pDot */
  dvalues[2] = 0.5 * AoverSqrtD * u3 * (  2.0 * ( q2 + p4 * z3) * A
                     - r * ( q2 + r2 + p4 * z3 ) * dAdr ) / HeffHreal
                     - (p / q) * (flux / (eta * omega));

  /* qDot */
  dvalues[3] = - flux / (eta * omega);

  return GSL_SUCCESS;
}

/**
 * Function which calculates omega = dphi/dt
 */
static
REAL8 XLALCalculateOmega(   REAL8 eta,                /**<< Symmetric mass ratio */
                            REAL8 r,                  /**<< Orbital separation (units of total mass M) */
                            REAL8 pr,                 /**<< Tortoise co-ordinate */
                            REAL8 pPhi,               /**<< Momentum pphi */
                            EOBACoefficients *aCoeffs /**<< Pre-computed EOB A coefficients */
                        )
{

  REAL8 A, Heff, Hreal, HeffHreal;

  A    = XLALCalculateEOBA( r, aCoeffs );
  Heff = XLALEffectiveHamiltonian( eta, r, pr, pPhi, aCoeffs );
  Hreal = sqrt( 1. + 2.*eta*(Heff - 1.) );

  HeffHreal = Heff * Hreal;
  return pPhi * A / (HeffHreal*r*r);
}

/**
 * Function which will determine whether to stop the evolution for the initial,
 * user-requested sample rate. We stop in this case when we have reached the peak
 * orbital frequency. 
 */
static int
XLALFirstStoppingCondition(double UNUSED t,              /**<< Current time (required by GSL) */
                           const double UNUSED values[], /**<< Current dynamics (required by GSL) */
                           double dvalues[],             /**<< Derivatives of dynamics w.r.t. time */
                           void *funcParams              /**<< Structure containing necessary parameters */
                          )
{

  EOBParams *params = (EOBParams *)funcParams;
  double omega = dvalues[1];

  if ( omega < params->omega )
  {
    return 1;
  }

  params->omega = omega;
  return GSL_SUCCESS;
}

/**
 * Function which will determine whether to stop the evolution for the high sample rate.
 * In this case, the data obtained will be used to attach the ringdown, so to make sure 
 * we won't be interpolating data too near the final points, we push this integration
 * as far as is feasible.
 */
static int
XLALHighSRStoppingCondition(double UNUSED t,       /**<< Current time (required by GSL) */
                           const double values[],  /**<< Current dynamics  */
                           double dvalues[],       /**<< Derivatives of dynamics w.r.t. time */
                           void UNUSED *funcParams /**<< Structure containing necessary parameters */
                          )
{
  EOBParams *params = (EOBParams *)funcParams;
  REAL8 rstop;
  if ( params->eta > 0.1 )
  {
    rstop = 1.25 - params->eta;
  }
  else
  {
    rstop = 2.1 - 10.0 * params->eta;
  }

  if ( values[0] <= rstop || isnan(dvalues[3]) || isnan (dvalues[2]) || isnan (dvalues[1]) || isnan (dvalues[0]) )
  {
    return 1;
  }
  return 0;
}

/*-------------------------------------------------------------------*/

/**
 * Calculates the initial radial velocity
 */
REAL8 XLALvrP4PN( const REAL8 r,    /**<< Orbital separation (in units of total mass M) */
                 const REAL8 omega, /**<< Orbital frequency (dimensionless: M*omega)*/
                 pr3In *params      /**<< pr3In structure containing some necessary parameters */
                )
{
  REAL8 A, dAdr, d2Adr2, dA, d2A, NA, DA, dDA, dNA, d2DA, d2NA;
  REAL8 r2, r3, r4, r5, u, u2, v, x1;
  REAL8 eta, FDIS;
  REAL8 twoUAPlusu2dA;

  expnCoeffsdEnergyFlux ak;

  EOBACoefficients *aCoeffs = params->aCoeffs;

  eta = params->eta;
  r2 = r*r;
  r3 = r2*r;
  r4 = r2*r2;
  r5 = r2*r3;

  u = 1./ r;

  u2 = u*u;

  NA = r4 * aCoeffs->n4
     + r5 * aCoeffs->n5;

  DA = aCoeffs->d0
     + r  * aCoeffs->d1
     + r2 * aCoeffs->d2
     + r3 * aCoeffs->d3
     + r4 * aCoeffs->d4
     + r5 * aCoeffs->d5;

  dNA = 4. * aCoeffs->n4 * r3
      + 5. * aCoeffs->n5 * r4;

  dDA = aCoeffs->d1
      + 2. * aCoeffs->d2 * r
      + 3. * aCoeffs->d3 * r2
      + 4. * aCoeffs->d4 * r3
      + 5. * aCoeffs->d5 * r4;

  d2NA = 12. * aCoeffs->n4 * r2
       + 20. * aCoeffs->n5 * r3;

  d2DA = 2. * aCoeffs->d2
       + 6. * aCoeffs->d3 * r
       + 12. * aCoeffs->d4 * r2
       + 20. * aCoeffs->d5 * r3;

  A      = NA/DA;
  dAdr   = ( dNA * DA - dDA * NA ) / (DA*DA);
  d2Adr2 = (d2NA*DA - d2DA*NA - 2.*dDA*DA*dAdr)/(DA*DA);

  /* The next derivatives are with respect to u, not r */

  dA  = - r2 * dAdr;
  d2A =  r4 * d2Adr2 + 2.*r3 * dAdr;

  v = cbrt(omega);

  XLALSimIMREOBNRv2SetupFlux( &ak, eta);

  FDIS = - XLALSimInspiralFp8PP(v, &ak)/(eta*omega);

  twoUAPlusu2dA = 2.* u * A + u2 * dA;
  x1 = -r2 * sqrt (-dA * twoUAPlusu2dA * twoUAPlusu2dA * twoUAPlusu2dA )
                / (2.* u * dA * dA + A*dA - u * A * d2A);
  return (FDIS * x1);
}

/**
 * Calculates the time window over which the ringdown attachment takes
 * place. These values were calibrated to numerical relativity simulations,
 * and come from Pan et al, PRD84, 124052(2011). 
 * The time returned is in units of M.
 */
static REAL8
GetRingdownAttachCombSize( 
                         INT4 l, /**<< Mode l */
                         INT4 m  /**<< Mode m */
                         )
{

   switch ( l )
   {
     case 2:
       switch ( m )
       {
         case 2:
           return 5.;
           break;
         case 1:
           return 8.;
           break;
         default:
           XLAL_ERROR_REAL8( XLAL_EINVAL );
           break;
        }
        break;
     case 3:
       if ( m == 3 )
       {
         return 12.;
       }
       else
       {
         XLAL_ERROR_REAL8( XLAL_EINVAL );
       }
       break;
     case 4:
       if ( m == 4 )
       {
         return 9.;
       }
       else
       {
         XLAL_ERROR_REAL8( XLAL_EINVAL );
       }
       break;
     case 5:
       if ( m == 5 )
       {
         return 8.;
       }
       else
       {
         XLAL_ERROR_REAL8( XLAL_EINVAL );
       }
       break;
     default:
       XLAL_ERROR_REAL8( XLAL_EINVAL );
       break;
  }

  /* It should not be possible to get to this point */
  /* Put an return path here to avoid compiler warning */
  XLALPrintError( "We shouldn't ever reach this point!\n" );
  XLAL_ERROR_REAL8( XLAL_EINVAL );

}

/**
 * Sets up the various PN coefficients which are needed to calculate
 * the flux. This is only used in setting the initial conditions, as in
 * the waveform evolution, the flux comes from the waveform itself.
 * Some of these values are only really needed as intermediate steps on
 * the way to getting the Pade coefficients, but since they are contained
 * within the structure, we may as well set them.
 */
static int
XLALSimIMREOBNRv2SetupFlux(
             expnCoeffsdEnergyFlux *ak, /**<<PN expansion coefficients (only relevant fields will be populated)*/
             REAL8                 eta  /**<< Symmetric mass ratio */
             )
{

  REAL8 a1, a2, a3, a4, a5, a6, a7, a8;
  
  /* Powers of a */
  REAL8 a12, a22, a32, a42, a52, a62, a72, a23, a33, a43, a53, a34, a44;

  REAL8 c1, c2, c3, c4, c5, c6, c7, c8;
  
  REAL8 vlso, vpole;

  memset( ak, 0, sizeof( *ak ) );

  /* Taylor coefficients of flux. */
  ak->fPaN = ak->FTaN = 32.*eta*eta/5.;
  ak->FTa1 = 0.;
  ak->FTa2 = -1247./336.-35.*eta/12.;
  ak->FTa3 = 4.*LAL_PI;
  ak->FTa4 = -44711./9072.+9271.*eta/504.+65.*eta*eta/18.;
  ak->FTa5 = -(8191./672.+583./24.*eta)*LAL_PI;
  ak->FTl6 = -1712./105.;
  ak->FTa6 = 6643739519./69854400. + 16.*LAL_PI*LAL_PI/3. + ak->FTl6 * log (4.L)
           - 1712./105.*LAL_GAMMA + (-134543./7776. + 41.*LAL_PI*LAL_PI/48.)*eta
           - 94403./3024. *eta*eta - 775./324. * eta*eta*eta;
  ak->FTa7 = LAL_PI * (-16285./504. + 214745./1728.*eta
              + 193385./3024.*eta*eta);
  ak->FTa8 = - 117.5043907226773;
  ak->FTl8 =   52.74308390022676;

  /* For the EOBPP model, vpole and vlso are tuned to NR */
  vpole = ak->vpolePP = 0.85;
  vlso  = ak->vlsoPP  = 1.0;

  /* Pade coefficients of f(v);  assumes that a0=1 => c0=1 */
  a1 = ak->FTa1 - 1./vpole;
  a2 = ak->FTa2 - ak->FTa1/vpole;
  a3 = ak->FTa3 - ak->FTa2/vpole;
  a4 = ak->FTa4 - ak->FTa3/vpole;
  a5 = ak->FTa5 - ak->FTa4/vpole;
  a6 = ak->FTa6 - ak->FTa5/vpole + ak->FTl6*log(vlso);
  a7 = ak->FTa7 - ( ak->FTa6 + ak->FTl6*log(vlso))/vpole;
  a8 = ak->FTa8 - ak->FTa7/vpole + ak->FTl8*log(vlso);

  c1 = -a1;
  c2 = -(c1*c1 - a2)/c1;
  c3 = -(c1*pow(c2+c1,2.) + a3)/(c1*c2);
  c4 = -(c1*pow(c2+c1,3.) + c1*c2*c3*(c3+2*c2+2*c1) - a4)/(c1*c2*c3);
  c5 = -(c1*pow(pow(c1+c2,2.)+c2*c3,2.) + c1*c2*c3*pow(c1+c2+c3+c4,2.) + a5)
     /(c1*c2*c3*c4);
  c6 = -(c1*pow(c1+c2, 5.)
     + c1*c2*c3*(pow(c1+c2+c3,3.) + 3.*pow(c1+c2,3.)
     + 3.*c2*c3*(c1+c2) + c3*c3*(c2-c1))
     + c1*c2*c3*c4*(pow(c1+c2+c3+c4,2.)
     + 2.*pow(c1+c2+c3,2.) + c3*(c4-2.*c1))
     + c1*c2*c3*c4*c5*(c5 + 2.*(c1+c2+c3+c4))
     - a6)/(c1*c2*c3*c4*c5);

  a12 = a1*a1;
  a22 = a2*a2;
  a32 = a3*a3;
  a42 = a4*a4;
  a52 = a5*a5;
  a62 = a6*a6;
  a72 = a7*a7;
  a23 = a22*a2;
  a33 = a32*a3;
  a43 = a42*a4;
  a53 = a52*a5;
  a34 = a33*a3;
  a44 = a43*a4;

  c7 = ((a23 + a32 + a12*a4 - a2*(2*a1*a3 + a4)) *
     (-a44 - a32*a52 + a1*a53 - a22*a62 + a33*a7
     + a22*a5*a7 + a42*(3*a3*a5 + 2*a2*a6 + a1*a7)
     + a3*(2*a2*a5*a6 + a1*a62 - a1*a5*a7)
     - 2*a4*((a32 + a1*a5)*a6 + a2*(a52 + a3*a7))))/
     ((a33 + a1*a42 + a22*a5 - a3*(2*a2*a4 + a1*a5))*
     (a34 + a22*a42 - a43 - 2*a1*a2*a4*a5 + a12*a52
     - a2*a52 - a23*a6 - a12*a4*a6 + a2*a4*a6
     - a32*(3*a2*a4 + 2*a1*a5 + a6)
     + 2*a3*((a22 + a4)*a5 + a1*(a42 + a2*a6))));

  c8 = -(((a33 + a1*a42 + a22*a5 - a3*(2*a2*a4 + a1*a5)) *
     (pow(a4,5) + pow(a5,4) - 2*a33*a5*a6 + 2*a1*a3*a52*a6
     + 2*a3*a5*a62 + a12*a62*a6 - 2*a3*a52*a7
     + 2*a1*a32*a6*a7 - 2*a12*a5*a6*a7 + a32*a72
     + pow(a3,4)*a8 - 2*a1*a32*a5*a8 + a12*a52*a8
     - a32*a6*a8 - a43*(4*a3*a5 + 3*a2*a6 + 2*a1*a7 + a8)
     + a42*(3*a2*a52 + 3*a32*a6 + 4*a1*a5*a6 + a62
     + 4*a2*a3*a7 + 2*a5*a7 + a22*a8 + 2*a1*a3*a8)
     + a22*(a52*a6 - 2*a3*a6*a7 + 2*a3*a5*a8)
     + a22*a2*(a72 - a6*a8) - a4*(3*a52*a6 - 2*a22*a62 + 2*a33*a7
     + 4*a22*a5*a7 + a2*a72 - a2*a6*a8 - 3*a32*(a52 - a2*a8)
     + 2*a3*(a2*a5*a6 + 2*a1*a62 + a6*a7 - a5*a8)
     + 2*a1*(a52*a5 - a2*a6*a7 + a2*a5*a8) + a12*(-a72 + a6*a8))
     + a2*(-a62*a6 + 2*a5*a6*a7 + 2*a1*a5*(-a62 + a5*a7)
     + a32*(a62 + 2*a5*a7) - a52*a8
     - 2*a3*(a52*a5 + a1*(a72 - a6*a8)))))
     / ((pow(a3,4) + a22*a42 - a43 - 2*a1*a2*a4*a5
     + a12*a52 - a2*a52 - a22*a2*a6 - a12*a4*a6 + a2*a4*a6
     - a32*(3*a2*a4 + 2*a1*a5 + a6)
     + 2*a3*((a22 + a4)*a5 + a1*(a42 + a2*a6)))
     * (-pow(a4,4) - a32*a52 + a1*a52*a5 - a22*a62 + a33*a7
     + a22*a5*a7 + a42*(3*a3*a5 + 2*a2*a6 + a1*a7)
     + a3*(2*a2*a5*a6 + a1*a62 - a1*a5*a7)
     - 2*a4*((a32 + a1*a5)*a6 + a2*(a52 + a3*a7)))));

  ak->fPa1 = c1;
  ak->fPa2 = c2;
  ak->fPa3 = c3;
  ak->fPa4 = c4;
  ak->fPa5 = c5;
  ak->fPa6 = c6;
  ak->fPa7 = c7;
  ak->fPa8 = c8;

  return XLAL_SUCCESS;
}

/* return the index before the instantaneous frequency rises past target */
static size_t find_instant_freq(const REAL8TimeSeries *hp, const REAL8TimeSeries *hc, const REAL8 target, const size_t start, const int fsign) {
  size_t k = start + 1;
  const size_t n = hp->data->length - 1;

  /* Use second order differencing to find the instantaneous frequency as
   * h = A e^(2 pi i f t) ==> f = d/dt(h) / (2*pi*h) */
  for (; k < n; k++) {
    const REAL8 hpDot = (hp->data->data[k+1] - hp->data->data[k-1]) / (2 * hp->deltaT);
    const REAL8 hcDot = (hc->data->data[k+1] - hc->data->data[k-1]) / (2 * hc->deltaT);
    REAL8 f = hcDot * hp->data->data[k] - hpDot * hc->data->data[k];
    f /= LAL_TWOPI;
    f /= hp->data->data[k] * hp->data->data[k] + hc->data->data[k] * hc->data->data[k];
    if (fsign != 0) f = -f;
//printf("this f: %f\n",f);
    if (f >= target) return k - 1;
  }
//printf("target f: %f\n",target);
  printf("Error: initial frequency too high, no waveform generated");
  XLAL_ERROR(XLAL_EDOM);
}

static size_t find_instant_freq_hlm(const COMPLEX16TimeSeries *hlm,
    const REAL8 target, const size_t start) {
  size_t k = start + 1;
  const size_t n = hlm->data->length - 1;
  COMPLEX16 *dataPtr = hlm->data->data;
  /* Use second order differencing to find the instantaneous frequency as
   * h_lm = A e^(-2 pi i f t) ==> f = - d/dt(h_lm) / (2*pi*h_lm) */
  for (; k < n; k++) {
    const REAL8 hreDot = (creal(dataPtr[k+1]) - creal(dataPtr[k-1]))
        / (2 * hlm->deltaT);
    const REAL8 himDot = (cimag(dataPtr[k+1]) - cimag(dataPtr[k-1]))
        / (2 * hlm->deltaT);
    REAL8 f = hreDot * cimag(dataPtr[k]) - himDot * creal(dataPtr[k]);
    f /= LAL_TWOPI;
    f /= creal(dataPtr[k]) * creal(dataPtr[k])
        + cimag(dataPtr[k]) * cimag(dataPtr[k]);
    if (f >= target) return k - 1;
  }
  printf("Error: initial frequency too high, no waveform generated");
  XLAL_ERROR(XLAL_EDOM);
}

/* Engine function of EOBNRv2 */
static int
XLALSimIMREOBNRv2Generator(
              REAL8TimeSeries   **hplus,
              REAL8TimeSeries   **hcross,
              SphHarmTimeSeries **h_lms,
              const REAL8       phiC,
              const REAL8       deltaT,
              const REAL8       m1SI,
              const REAL8       m2SI,
              const REAL8       fLower,
              const REAL8       distance,
              const REAL8       inclination,
              const int         higherModeFlag
              )
{
   UINT4                   count, nn=4, hiSRndx=0;
   UINT4                   flag_fLower_extend = 0;
   size_t                  cut_ind = 0.;

   /* Vector containing the current signal mode */
   COMPLEX16TimeSeries     *sigMode = NULL;

   REAL8Vector             rVec, phiVec, prVec, pPhiVec;
   REAL8Vector             rVecHi, phiVecHi, prVecHi, pPhiVecHi, tVecHi;

   /* Masses in solar masses */
   REAL8                   mass1, mass2, totalMass;
   REAL8                   eta, m, r, s, p, q, dt, t, v, omega, f, ampl0;
   REAL8                   omegaOld;

   REAL8Vector             *values, *dvalues;

   /* Time to use as the epoch of the returned time series */
   LIGOTimeGPS             epoch = LIGOTIMEGPSZERO;

   /* Variables for the integrator */
   ark4GSLIntegrator       *integrator = NULL;
   REAL8Array              *dynamics   = NULL;
   REAL8Array              *dynamicsHi = NULL;
   INT4                    retLen;
   REAL8                   tMax;

   pr3In                   pr3in;

   /* Stuff for pre-computed EOB values */
   EOBParams eobParams;
   EOBACoefficients        aCoeffs;
   FacWaveformCoeffs       hCoeffs;
   NewtonMultipolePrefixes prefixes;

   /* Variables to allow the waveform to be generated */
   /* from a specific fLower */
   REAL8      sSub = 0.0;  /* Initial phase, and phase to subtract */

   REAL8      rmin = 20;        /* Smallest value of r at which to generate the waveform */
   COMPLEX16  hLM;             /* Factorized waveform */
   COMPLEX16  hNQC;            /* Non-quasicircular correction */
   UINT4      i, j, modeL;
   INT4       modeM;
   INT4       nModes;         /* number of modes required */

   /* Used for EOBNR */
   COMPLEX16Vector *modefreqs;
   UINT4 resampFac;
   UINT4 resampPwr; /* Power of 2 for resampling */
   REAL8 resampEstimate;

   /* For checking XLAL return codes */
   INT4 xlalStatus;

   /* Accuracy of root finding algorithms */
   const REAL8 xacc = 1.0e-12;

   /* Min and max rInit and prInit to consider in root finding */
   const REAL8 rInitMin = 3.;
   const REAL8 rInitMax = 1000.;
   
   const REAL8 prInitMin = -10.;
   const REAL8 prInitMax = 5.;

   /* Accuracies of adaptive Runge-Kutta integrator */
   const REAL8 EPS_ABS = 1.0e-12;
   const REAL8 EPS_REL = 1.0e-10;

   REAL8 tStepBack; /* We need to step back 20M to attach ringdown */
   UINT4 nStepBack; /* Num points to step back */

   /* Stuff at higher sample rate */
   /* Note that real and imaginary parts of signal mode are in separate vectors */
   /* This is necessary for the ringdown attachment code */
   REAL8Vector             *sigReHi, *sigImHi;
   REAL8Vector             *phseHi, *omegaHi;
   UINT4                   lengthHiSR;

   /* Used in the calculation of the non-quasicircular correctlon */
   REAL8Vector             *ampNQC, *q1, *q2, *q3, *p1, *p2;
   EOBNonQCCoeffs           nqcCoeffs;

   /* Inidices of the ringdown matching points */
   /* peakIdx is the index where omega is a maximum */
   /* finalIdx is the index of the last point before the */
   /* integration breaks */
   /* startIdx is the index at which the waveform crosses fLower */
   REAL8Vector             *rdMatchPoint;
   UINT4                   peakIdx  = 0;
   UINT4                   finalIdx = 0;
   UINT4                   startIdx = 0;
   /* Counter used in unwrapping the waveform phase for NQC correction */
   INT4 phaseCounter;

   /* The list of available modes */
   const INT4 lmModes[5][2] = {{2, 2},
                               {2, 1},
                               {3, 3},
                               {4, 4},
                               {5, 5}};

   INT4 currentMode;

   /* Checks on input */
   if ( (!hplus || !hcross) && !h_lms )
   {
     XLAL_ERROR( XLAL_EFAULT );
   }

   if ( hplus && hcross ) // If given polarization double pointers **hp, **hc...
   {
     if ( *hplus || *hcross ) // ...Make sure *hp, *hc are NULL
     {
       XLALPrintError( "(*hplus) and (*hcross) are expected to be NULL; got %p and %p\n",
           *hplus, *hcross);
       XLAL_ERROR( XLAL_EFAULT );
     }
   }
   if ( h_lms ) // If given Sph. Harm. double pointer **h_lms...
   {
     if ( *h_lms ) // ...Make sure *h_lms is NULL
     {
       XLALPrintError( "(*h_lms) is expected to be NULL; got %p\n",
           *h_lms);
       XLAL_ERROR( XLAL_EFAULT );
     }
   }

   if ( distance <= 0.0 )
   {
     XLALPrintError( "Distance must be > 0.\n" );
     XLAL_ERROR( XLAL_EINVAL );
   }

   if ( m1SI <= 0. || m2SI <= 0. )
   {
     XLALPrintError( "Component masses must be > zero!\n" );
     XLAL_ERROR( XLAL_EFUNC );
   }

   /* Allocate some memory */
   values    = XLALCreateREAL8Vector( nn );
   dvalues   = XLALCreateREAL8Vector( nn );

   if ( !values || !dvalues )
   {
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
     XLAL_ERROR( XLAL_EFUNC );
   }

   /* From this point on, we will want to use masses in solar masses */
   mass1 = m1SI / LAL_MSUN_SI;
   mass2 = m2SI / LAL_MSUN_SI;

   totalMass = mass1 + mass2;

   eta = mass1 * mass2 / (totalMass*totalMass);

   /* We will also need the total mass in units of seconds */
   m = totalMass * LAL_MTSUN_SI;

   dt  = deltaT;

   /* The Runge-Kutta integrator needs an estimate of the length of the waveform */
   /* It doesn't really matter, as the length will be determined by the stopping condition */
   tMax = 5.;

   /* Set the amplitude depending on whether the distance is given */
   ampl0 = totalMass * LAL_MRSUN_SI/ distance;

   /* Check we get a sensible answer */
   if ( ampl0 == 0.0 )
   {
     XLALPrintWarning( "Generating waveform of zero amplitude!!\n" );
   }

   /* Set the number of modes depending on whether the user wants higher order modes */
   if ( higherModeFlag == 0 )
   {
     nModes = 1;
   }
   else if ( higherModeFlag == 1 )
   {
     nModes = EOBNRV2_NUM_MODES_MAX;
   }
   else
   {
     XLALPrintError( "Higher mode flag appears to be uninitialised " 
         "(expected 0 or 1, but got %d\n)", higherModeFlag );
     XLAL_ERROR( XLAL_EINVAL );
   }

   /* Check that the (l,m,0) QNM freq. is less than the Nyquist freq. */
   modefreqs = XLALCreateCOMPLEX16Vector( 3 );
   for ( currentMode = 0; currentMode < nModes; currentMode++ )
   {
     count = 0;

     modeL = lmModes[currentMode][0];
     modeM = lmModes[currentMode][1];
     /* Get QNM frequencies */
     xlalStatus = XLALSimIMREOBGenerateQNMFreqV2( modefreqs, mass1, mass2, NULL, NULL, modeL, modeM, 3, EOBNRv2);
     if ( xlalStatus != XLAL_SUCCESS )
     {
       XLALDestroyCOMPLEX16Vector( modefreqs );
       XLALDestroyREAL8Vector( values );
       XLALDestroyREAL8Vector( dvalues );
       XLAL_ERROR( XLAL_EFUNC );
     }
     /* If Nyquist freq. <  (l,m,0) QNM freq., exit */
     /* Note that we cancelled a factor of 2 occuring on both sides */
     if ( LAL_PI / creal(modefreqs->data[0]) < dt )
     {
       XLALDestroyCOMPLEX16Vector( modefreqs );
       XLALDestroyREAL8Vector( values );
       XLALDestroyREAL8Vector( dvalues );
       XLALPrintError( "(%d,%d) mode ringdown freq greater than Nyquist freq. "
             "Increase sample rate or consider using EOB approximant.\n",modeL,modeM );
       XLAL_ERROR( XLAL_EINVAL );
     }
   }

   /* Calculate the time we will need to step back for ringdown */
   tStepBack = 20.0 * m;
   nStepBack = ceil( tStepBack / dt );

   /* Set up structures for pre-computed EOB coefficients */
   memset( &eobParams, 0, sizeof(eobParams) );
   eobParams.eta = eta;
   eobParams.m1  = mass1;
   eobParams.m2  = mass2;
   eobParams.aCoeffs   = &aCoeffs;
   eobParams.hCoeffs   = &hCoeffs;
   eobParams.nqcCoeffs = &nqcCoeffs;
   eobParams.prefixes  = &prefixes;

   if ( XLALCalculateEOBACoefficients( &aCoeffs, eta ) == XLAL_FAILURE 
   ||   XLALSimIMREOBCalcFacWaveformCoefficients( &hCoeffs, eta) == XLAL_FAILURE 
   ||   XLALSimIMREOBComputeNewtonMultipolePrefixes( &prefixes, eobParams.m1, eobParams.m2 )
         == XLAL_FAILURE )
   {
     XLALDestroyCOMPLEX16Vector( modefreqs );
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
     XLAL_ERROR( XLAL_EFUNC );
   }

   /* For the dynamics, we need to use preliminary calculated versions   */
   /*of the NQC coefficients for the (2,2) mode. We calculate them here. */
   if ( XLALSimIMREOBGetCalibratedNQCCoeffs( &nqcCoeffs, 2, 2, eta ) == XLAL_FAILURE )
   {
     XLALDestroyCOMPLEX16Vector( modefreqs );
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
     XLAL_ERROR( XLAL_EFUNC );
   }

   /* Calculate the resample factor for attaching the ringdown */
   /* We want it to be a power of 2 */
   /* Of course, we only want to do this if the required SR > current SR... */
   /* The form chosen for the resampleEstimate will essentially set */
   /* deltaT = M / 20. ( or less taking into account the power of 2 stuff */
   resampEstimate = 20.* dt / m;
   resampFac      = 1;

   if ( resampEstimate > 1 )
   {
     resampPwr = (UINT4)ceil(log2(resampEstimate));
     while ( resampPwr-- )
     {
       resampFac *= 2u;
     }
   }

   /* The length of the vectors at the higher sample rate will be */
   /* the step back time plus the ringdown */
   lengthHiSR = ( nStepBack + (UINT4)(2. * EOB_RD_EFOLDS / cimag(modefreqs->data[0]) / dt) ) * resampFac;

   /* Double it for good measure */
   lengthHiSR *= 2;

   /* We are now done with the ringdown modes - destroy them here */
   XLALDestroyCOMPLEX16Vector( modefreqs );
   modefreqs = NULL;

   /* Find the initial velocity given the lower frequency */
   f     = fLower;
   omega = f * LAL_PI * m;
   v     = cbrt( omega );

   /* initial r as a function of omega - where to start evolution */
   pr3in.eta = eta;
   pr3in.aCoeffs = &aCoeffs;

   /* We will be changing the starting r if it is less than rmin */
   /* Therefore, we should reset pr3in.omega later if necessary. */
   /* For now we need it so that we can see what the initial r is. */

   pr3in.omega = omega;

   /* if ( XLALrOfOmegaP4PN(rInitMin, &pr3in) < 0.)
   {
     XLALPrintError( "Initial orbital frequency too high. The corresponding initial radius < %fM\n", rInitMin);
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
     XLAL_ERROR( XLAL_EFUNC );
   } */
   /* If initial frequency too high and initial r < 10M, start at r = 10M, set flag and remove low freq waveform at the end */
   if ( XLALrOfOmegaP4PN(10., &pr3in) < 0.)
   {
     flag_fLower_extend = 1;
     omega = pow(10., -1.5);
     pr3in.omega = omega;
     v = cbrt( omega );
     f = omega / LAL_PI / m;
   }
   if ( XLALrOfOmegaP4PN(rInitMax, &pr3in) > 0.)
   {
     XLALPrintError( "Initial orbital frequency too low. The corresponding initial radius > %fM\n", rInitMax);
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
     XLAL_ERROR( XLAL_EFUNC );
   }
   r = XLALDBisectionFindRoot( XLALrOfOmegaP4PN, rInitMin,
              rInitMax, xacc, &pr3in);
   if ( XLAL_IS_REAL8_FAIL_NAN( r ) )
   {
     XLALPrintError( "Failed solving the initial radius. The desired initial orbital frequency is %e\n", omega );
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
     XLAL_ERROR( XLAL_EFUNC );
   }

   /* We want the waveform to generate from a point which won't cause */
   /* problems with the initial conditions. Therefore we force the code */
   /* to start at least at r = rmin (in units of M). */

   r = (r<rmin) ? rmin : r;
   pr3in.r = r;

   /* Now that r is changed recompute omega corresponding */
   /* to that r and only then compute initial pr and pphi */

   omega = omegaofrP4PN( r, eta, &aCoeffs );
   pr3in.omega = omega;
   q = XLALpphiInitP4PN(r, &aCoeffs );
   /* first we compute vr (we need coeef->Fp6) */
   pr3in.q = q;
   pr3in.vr = XLALvrP4PN(r, omega, &pr3in);
   /* then we compute the initial value of p */
   p = XLALDBisectionFindRoot( XLALprInitP4PN, prInitMin, prInitMax, xacc, &pr3in);
   if ( XLAL_IS_REAL8_FAIL_NAN( p ) )
   {
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
     XLAL_ERROR( XLAL_EFUNC );
   }
   /* We need to change P to be the tortoise co-ordinate */
   /* TODO: Change prInit to calculate this directly */
   p = p * XLALCalculateEOBA(r, &aCoeffs);
   p = p / sqrt( XLALCalculateEOBD( r, eta ) );

   s = 0.0;

   values->data[0] = r;
   values->data[1] = s;
   values->data[2] = p;
   values->data[3] = q;

   /* And their higher sample rate counterparts */
   /* Allocate memory for temporary arrays */
   sigReHi  = XLALCreateREAL8Vector ( lengthHiSR );
   sigImHi  = XLALCreateREAL8Vector ( lengthHiSR );
   phseHi  = XLALCreateREAL8Vector ( lengthHiSR );
   omegaHi = XLALCreateREAL8Vector ( lengthHiSR );

   /* Allocate NQC vectors */
   ampNQC = XLALCreateREAL8Vector ( lengthHiSR );
   q1     = XLALCreateREAL8Vector ( lengthHiSR );
   q2     = XLALCreateREAL8Vector ( lengthHiSR );
   q3     = XLALCreateREAL8Vector ( lengthHiSR );
   p1     = XLALCreateREAL8Vector ( lengthHiSR );
   p2     = XLALCreateREAL8Vector ( lengthHiSR );

   if ( !sigReHi || !sigImHi || !phseHi )
   {
     XLALDestroyREAL8Vector( sigReHi );
     XLALDestroyREAL8Vector( sigImHi );
     XLALDestroyREAL8Vector( phseHi );
     XLALDestroyREAL8Vector( omegaHi );
     XLALDestroyREAL8Vector( ampNQC );
     XLALDestroyREAL8Vector( q1 );
     XLALDestroyREAL8Vector( q2 );
     XLALDestroyREAL8Vector( q3 );
     XLALDestroyREAL8Vector( p1 );
     XLALDestroyREAL8Vector( p2 );
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
     XLAL_ERROR( XLAL_ENOMEM );
   }

   memset(sigReHi->data, 0, sigReHi->length * sizeof( REAL8 ));
   memset(sigImHi->data, 0, sigImHi->length * sizeof( REAL8 ));
   memset(phseHi->data, 0, phseHi->length * sizeof( REAL8 ));
   memset(omegaHi->data, 0, omegaHi->length * sizeof( REAL8 ));
   memset(ampNQC->data, 0, ampNQC->length * sizeof( REAL8 ));
   memset(q1->data, 0, q1->length * sizeof( REAL8 ));
   memset(q2->data, 0, q2->length * sizeof( REAL8 ));
   memset(q3->data, 0, q3->length * sizeof( REAL8 ));
   memset(p1->data, 0, p1->length * sizeof( REAL8 ));
   memset(p2->data, 0, p2->length * sizeof( REAL8 ));

   /* Initialize the GSL integrator */
   if (!(integrator = XLALAdaptiveRungeKutta4Init(nn, LALHCapDerivativesP4PN, XLALFirstStoppingCondition, EPS_ABS, EPS_REL)))
   {
     XLALDestroyREAL8Vector( sigReHi );
     XLALDestroyREAL8Vector( sigImHi );
     XLALDestroyREAL8Vector( phseHi );
     XLALDestroyREAL8Vector( omegaHi );
     XLALDestroyREAL8Vector( ampNQC );
     XLALDestroyREAL8Vector( q1 );
     XLALDestroyREAL8Vector( q2 );
     XLALDestroyREAL8Vector( q3 );
     XLALDestroyREAL8Vector( p1 );
     XLALDestroyREAL8Vector( p2 );
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
     XLAL_ERROR( XLAL_EFUNC );
   }

   integrator->stopontestonly = 1;

   count = 0;

   /* Use the new adaptive integrator */
   /* TODO: Implement error checking */
   retLen = XLALAdaptiveRungeKutta4( integrator, &eobParams, values->data, 0., tMax/m, dt/m, &dynamics );

   /* We should have integrated to the peak of the frequency by now */
   hiSRndx = retLen - nStepBack;

   /* Set up the vectors, and re-initialize everything for the high sample rate */
   rVec.length  = phiVec.length = prVec.length = pPhiVec.length = retLen;
   rVec.data    = dynamics->data+retLen;
   phiVec.data  = dynamics->data+2*retLen;
   prVec.data   = dynamics->data+3*retLen;
   pPhiVec.data = dynamics->data+4*retLen;

   dt = dt/(REAL8)resampFac;
   values->data[0] = rVec.data[hiSRndx];
   values->data[1] = phiVec.data[hiSRndx];
   values->data[2] = prVec.data[hiSRndx];
   values->data[3] = pPhiVec.data[hiSRndx];

   /* We want to use a different stopping criterion for the higher sample rate */
   integrator->stop = XLALHighSRStoppingCondition;

   retLen = XLALAdaptiveRungeKutta4( integrator, &eobParams, values->data,
     0, (lengthHiSR-1)*dt/m, dt/m, &dynamicsHi );

   rVecHi.length  = phiVecHi.length = prVecHi.length = pPhiVecHi.length = tVecHi.length = retLen;
   rVecHi.data    = dynamicsHi->data+retLen;
   phiVecHi.data  = dynamicsHi->data+2*retLen;
   prVecHi.data   = dynamicsHi->data+3*retLen;
   pPhiVecHi.data = dynamicsHi->data+4*retLen;
   tVecHi.data    = dynamicsHi->data;

   /* We are now finished with the adaptive RK, so we can free its resources */
   XLALAdaptiveRungeKutta4Free( integrator );
   integrator = NULL;

   /* Now we have the dynamics, we tweak the factorized coefficients for the waveform */
  if ( XLALSimIMREOBModifyFacWaveformCoefficients( &hCoeffs, eta) == XLAL_FAILURE )
  {
     XLALDestroyREAL8Vector( sigReHi );
     XLALDestroyREAL8Vector( sigImHi );
     XLALDestroyREAL8Vector( phseHi );
     XLALDestroyREAL8Vector( omegaHi );
     XLALDestroyREAL8Vector( ampNQC );
     XLALDestroyREAL8Vector( q1 );
     XLALDestroyREAL8Vector( q2 );
     XLALDestroyREAL8Vector( q3 );
     XLALDestroyREAL8Vector( p1 );
     XLALDestroyREAL8Vector( p2 );
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
     XLAL_ERROR( XLAL_EFUNC );
  }

  /* We can now prepare to output the waveform */
  /* We want to start outputting when the 2,2 mode crosses the user-requested fLower */
  /* Find the point where we reach the low frequency cutoff */
  REAL8 lfCut = fLower * LAL_PI*m;
  if (flag_fLower_extend == 1)
  {
    lfCut = pow(10., -1.5);
  }

  i = 0;
  /* TODO: discrete search for lfCut generates discontinuity w.r.t change in physical parameter. 
           Implement a continuous search. */
  while ( i < hiSRndx )
  {
    omega = XLALCalculateOmega( eta, rVec.data[i], prVec.data[i], pPhiVec.data[i], &aCoeffs );
    if ( omega > lfCut || fabs( omega - lfCut ) < 1.0e-5 )
    {
      break;
    }
    i++;
  }

  if ( i == hiSRndx )
  {
     XLALDestroyREAL8Vector( sigReHi );
     XLALDestroyREAL8Vector( sigImHi );
     XLALDestroyREAL8Vector( phseHi );
     XLALDestroyREAL8Vector( omegaHi );
     XLALDestroyREAL8Vector( ampNQC );
     XLALDestroyREAL8Vector( q1 );
     XLALDestroyREAL8Vector( q2 );
     XLALDestroyREAL8Vector( q3 );
     XLALDestroyREAL8Vector( p1 );
     XLALDestroyREAL8Vector( p2 );
     XLALDestroyREAL8Vector( values );
     XLALDestroyREAL8Vector( dvalues );
    XLALPrintError( "Low frequency cut-off is too close to coalescence frequency.\n" );
    XLAL_ERROR( XLAL_EFAILED );
  }

  startIdx = i;

  omegaOld = 0.0;
  phaseCounter = 0;
  for ( i=0; i < (UINT4)retLen; i++ )
  {
    omega = XLALCalculateOmega( eta, rVecHi.data[i], prVecHi.data[i], pPhiVecHi.data[i], &aCoeffs );
    omegaHi->data[i] = omega;
    /* For now we re-populate values - there may be a better way to do this */
    values->data[0] = r = rVecHi.data[i];
    values->data[1] = s = phiVecHi.data[i] - sSub;
    values->data[2] = p = prVecHi.data[i];
    values->data[3] = q = pPhiVecHi.data[i];

    if ( omega <= omegaOld && !peakIdx )
    {
      peakIdx = i-1;
    }
    omegaOld = omega;
  }
  finalIdx = retLen - 1;

  /* Stuff to find the actual peak time */
  gsl_spline    *spline = NULL;
  gsl_interp_accel *acc = NULL;
  REAL8 omegaDeriv1;
  REAL8 time1, time2;   
  REAL8 timePeak, omegaDerivMid;

  spline = gsl_spline_alloc( gsl_interp_cspline, retLen );
  acc    = gsl_interp_accel_alloc();

  time1 = dynamicsHi->data[peakIdx];

  gsl_spline_init( spline, dynamicsHi->data, omegaHi->data, retLen );
  omegaDeriv1 = gsl_spline_eval_deriv( spline, time1, acc );
  if ( omegaDeriv1 > 0. )
  {
    time2 = dynamicsHi->data[peakIdx+1];
  }
  else
  {
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
      time2 = timePeak;
    }
    else
    {
      omegaDeriv1 = omegaDerivMid;
      time1 = timePeak;
    }
  }
  while ( time2 - time1 > 1.0e-5 );

  /* XLALPrintInfo( "Estimation of the peak is now at time %e\n", timePeak ); */

  /* Set the coalescence time and phase */
  /* It is not easy to define an exact coalescence time and phase for an IMR time-domain model */
  /* Therefore we set them at the time when the orbital frequency reaches maximum */
  /* Note that the coalescence phase is defined for the ORBITAL phase, not the GW phasae */
  /* With PN corrections in the GW modes, GW phase is not exactly m times orbital phase */
  /* In brief, at the highest orbital frequency, the orbital phase is phiC/2 */ 
  //t = m * (dynamics->data[hiSRndx] + dynamicsHi->data[peakIdx] - dynamics->data[startIdx]);
  t = m * (dynamics->data[hiSRndx] + timePeak - dynamics->data[startIdx]);
  gsl_spline_init( spline, dynamicsHi->data, phiVecHi.data, retLen );
  /* sSub = phiVecHi.data[peakIdx] - phiC/2.; */
  sSub = gsl_spline_eval( spline, timePeak, acc ) - phiC;

  gsl_spline_free( spline );
  gsl_interp_accel_free( acc );

  XLALGPSAdd( &epoch, -t);

  /* Allocate vectors for polarizations if they are being output */
  /* Their length should be the length of the inspiral + the merger/ringdown */
  if( hplus && hcross )
  {
    *hplus = XLALCreateREAL8TimeSeries( "H_PLUS", &epoch, 0.0, deltaT,
        &lalStrainUnit, hiSRndx - startIdx + lengthHiSR / resampFac );
    *hcross = XLALCreateREAL8TimeSeries( "H_CROSS", &epoch, 0.0, deltaT,
        &lalStrainUnit, (*hplus)->data->length );

    if ( !(*hplus) || !(*hcross) ) // Check hp, hc allocated properly
    {
      if ( *hplus )  
      { 
        XLALDestroyREAL8TimeSeries( *hplus ); 
        *hplus = NULL;
      }
      if ( *hcross )  
      { 
        XLALDestroyREAL8TimeSeries( *hcross ); 
        *hcross = NULL;
      }
      XLALDestroyREAL8Vector( sigReHi );
      XLALDestroyREAL8Vector( sigImHi );
      XLALDestroyREAL8Vector( phseHi );
      XLALDestroyREAL8Vector( omegaHi );
      XLALDestroyREAL8Vector( ampNQC );
      XLALDestroyREAL8Vector( q1 );
      XLALDestroyREAL8Vector( q2 );
      XLALDestroyREAL8Vector( q3 );
      XLALDestroyREAL8Vector( p1 );
      XLALDestroyREAL8Vector( p2 );
      XLALDestroyREAL8Vector( values );
      XLALDestroyREAL8Vector( dvalues );
    }
    memset( (*hplus)->data->data, 0, (*hplus)->data->length * sizeof(REAL8) );
    memset( (*hcross)->data->data, 0, (*hcross)->data->length * sizeof(REAL8) );
  }

  /* We can now start calculating things for NQCs, and hiSR waveform */

  for ( currentMode = 0; currentMode < nModes; currentMode++ )
  {
     sigMode = XLALCreateCOMPLEX16TimeSeries( "H_MODE", &epoch, 0.0, deltaT,
         &lalStrainUnit, hiSRndx - startIdx + lengthHiSR / resampFac );

     memset(sigMode->data->data, 0, sigMode->data->length * sizeof(COMPLEX16) );
     count = 0;

     modeL = lmModes[currentMode][0];
     modeM = lmModes[currentMode][1];

     /* If we have an equal mass system, some modes will be zero */
     if ( eta == 0.25 && modeM % 2 )
     {
       continue;
     }

    phaseCounter = 0;
    for ( i=0; i < (UINT4)retLen; i++ )
    {
      omega = XLALCalculateOmega( eta, rVecHi.data[i], prVecHi.data[i], pPhiVecHi.data[i], &aCoeffs );
      omegaHi->data[i] = omega;
      /* For now we re-populate values - there may be a better way to do this */
      values->data[0] = r = rVecHi.data[i];
      values->data[1] = s = phiVecHi.data[i] - sSub;
      values->data[2] = p = prVecHi.data[i];
      values->data[3] = q = pPhiVecHi.data[i];

      v = cbrt( omega );

      xlalStatus = XLALSimIMREOBGetFactorizedWaveform( &hLM, values, v, modeL, modeM, &eobParams );

      ampNQC->data[i] = cabs( hLM );
      sigReHi->data[i] = (REAL8) ampl0 * creal(hLM);
      sigImHi->data[i] = (REAL8) ampl0 * cimag(hLM);
      phseHi->data[i] = carg( hLM ) + phaseCounter * LAL_TWOPI;
      if ( i && phseHi->data[i] > phseHi->data[i-1] )
      {
        phaseCounter--;
        phseHi->data[i] -= LAL_TWOPI;
      }
      q1->data[i] = p*p / (r*r*omega*omega);
      q2->data[i] = q1->data[i] / r;
      q3->data[i] = q2->data[i] / sqrt(r);
      p1->data[i] = p / ( r*omega );
      p2->data[i] = p1->data[i] * p*p;
    }

    /* Calculate the NQC correction */
    XLALSimIMREOBCalculateNQCCoefficients( &nqcCoeffs, ampNQC, phseHi, q1,q2,q3,p1,p2, modeL, modeM, timePeak, dt/m, eta );

    /* Now create the (low-sampled) part of the waveform */
    i = startIdx;
    while ( i < hiSRndx )
    {
       omega = XLALCalculateOmega( eta, rVec.data[i], prVec.data[i], pPhiVec.data[i], &aCoeffs );
       /* For now we re-populate values - there may be a better way to do this */
       values->data[0] = r = rVec.data[i];
       values->data[1] = s = phiVec.data[i] - sSub;
       values->data[2] = p = prVec.data[i];
       values->data[3] = q = pPhiVec.data[i];

       v = cbrt( omega );

       xlalStatus = XLALSimIMREOBGetFactorizedWaveform( &hLM, values, v, modeL, modeM, &eobParams );

       xlalStatus = XLALSimIMREOBNonQCCorrection( &hNQC, values, omega, &nqcCoeffs );
       hLM *= hNQC;

       sigMode->data->data[count] = hLM * ampl0;

       count++;
       i++;
    }

    /* Now apply the NQC correction to the high sample part */
    for ( i = 0; i <= finalIdx; i++ )
    {
      omega = XLALCalculateOmega( eta, rVecHi.data[i], prVecHi.data[i], pPhiVecHi.data[i], &aCoeffs );

      /* For now we re-populate values - there may be a better way to do this */
      values->data[0] = r = rVecHi.data[i];
      values->data[1] = s = phiVecHi.data[i] - sSub;
      values->data[2] = p = prVecHi.data[i];
      values->data[3] = q = pPhiVecHi.data[i];

      xlalStatus = XLALSimIMREOBNonQCCorrection( &hNQC, values, omega, &nqcCoeffs );

      hLM = sigReHi->data[i];
      hLM += I * sigImHi->data[i];

      hLM *= hNQC;
      sigReHi->data[i] = creal(hLM);
      sigImHi->data[i] = cimag(hLM);
    }

     /*--------------------------------------------------------------
      * Attach the ringdown waveform to the end of inspiral
       -------------------------------------------------------------*/
     rdMatchPoint = XLALCreateREAL8Vector( 3 );
     if ( !rdMatchPoint )
     {
       XLALDestroyCOMPLEX16TimeSeries( sigMode );
       XLALDestroyREAL8TimeSeries( *hplus );  *hplus  = NULL;
       XLALDestroyREAL8TimeSeries( *hcross ); *hcross = NULL;
       XLALDestroySphHarmTimeSeries( *h_lms ); *h_lms = NULL;
       XLALDestroyREAL8Vector( sigReHi );
       XLALDestroyREAL8Vector( sigImHi );
       XLALDestroyREAL8Vector( phseHi );
       XLALDestroyREAL8Vector( omegaHi );
       XLALDestroyREAL8Vector( ampNQC );
       XLALDestroyREAL8Vector( q1 );
       XLALDestroyREAL8Vector( q2 );
       XLALDestroyREAL8Vector( q3 );
       XLALDestroyREAL8Vector( p1 );
       XLALDestroyREAL8Vector( p2 );
       XLALDestroyREAL8Vector( values );
       XLALDestroyREAL8Vector( dvalues );
       XLAL_ERROR( XLAL_ENOMEM );
     }

     /* Check the first matching point is sensible */
     if ( ceil( tStepBack / ( 2.0 * dt ) ) > peakIdx )
     {
       XLALPrintError( "Invalid index for first ringdown matching point.\n" );
       XLALDestroyCOMPLEX16TimeSeries( sigMode );
       XLALDestroyREAL8TimeSeries( *hplus );  *hplus  = NULL;
       XLALDestroyREAL8TimeSeries( *hcross ); *hcross = NULL;
       XLALDestroySphHarmTimeSeries( *h_lms ); *h_lms = NULL;
       XLALDestroyREAL8Vector( sigReHi );
       XLALDestroyREAL8Vector( sigImHi );
       XLALDestroyREAL8Vector( phseHi );
       XLALDestroyREAL8Vector( omegaHi );
       XLALDestroyREAL8Vector( ampNQC );
       XLALDestroyREAL8Vector( q1 );
       XLALDestroyREAL8Vector( q2 );
       XLALDestroyREAL8Vector( q3 );
       XLALDestroyREAL8Vector( p1 );
       XLALDestroyREAL8Vector( p2 );
       XLALDestroyREAL8Vector( values );
       XLALDestroyREAL8Vector( dvalues );
       XLAL_ERROR( XLAL_EFAILED );
     }

     REAL8 combSize = GetRingdownAttachCombSize( modeL, modeM );
     REAL8 nrPeakDeltaT = XLALSimIMREOBGetNRPeakDeltaT( modeL, modeM, eta );

     if ( combSize > timePeak )
     {
       XLALPrintWarning( "Comb size not as big as it should be\n" );
     }

     rdMatchPoint->data[0] = combSize < timePeak + nrPeakDeltaT ? timePeak + nrPeakDeltaT - combSize : 0;
     rdMatchPoint->data[1] = timePeak + nrPeakDeltaT;
     rdMatchPoint->data[2] = dynamicsHi->data[finalIdx];

     xlalStatus = XLALSimIMREOBHybridAttachRingdown(sigReHi, sigImHi,
                   modeL, modeM, dt, mass1, mass2, 0, 0, 0, 0, 0, 0, &tVecHi, rdMatchPoint, EOBNRv2 );
     if (xlalStatus != XLAL_SUCCESS )
     {
       XLALDestroyREAL8Vector( rdMatchPoint );
       XLALDestroyCOMPLEX16TimeSeries( sigMode );
       XLALDestroyREAL8TimeSeries( *hplus );  *hplus  = NULL;
       XLALDestroyREAL8TimeSeries( *hcross ); *hcross = NULL;
       XLALDestroySphHarmTimeSeries( *h_lms ); *h_lms = NULL;
       XLALDestroyREAL8Vector( sigReHi );
       XLALDestroyREAL8Vector( sigImHi );
       XLALDestroyREAL8Vector( phseHi );
       XLALDestroyREAL8Vector( omegaHi );
       XLALDestroyREAL8Vector( ampNQC );
       XLALDestroyREAL8Vector( q1 );
       XLALDestroyREAL8Vector( q2 );
       XLALDestroyREAL8Vector( q3 );
       XLALDestroyREAL8Vector( p1 );
       XLALDestroyREAL8Vector( p2 );
       XLALDestroyREAL8Vector( values );
       XLALDestroyREAL8Vector( dvalues );
       XLAL_ERROR( XLAL_EFUNC );
     }

     XLALDestroyREAL8Vector( rdMatchPoint );

     for(j=0; j<sigReHi->length; j+=resampFac)
     {
       sigMode->data->data[count] = sigReHi->data[j];
       sigMode->data->data[count] += I * sigImHi->data[j];
       count++;
     }

     /* Add mode to final output h+ and hx */
     if( hplus && hcross)
     {
       XLALSimAddMode( *hplus, *hcross, sigMode, inclination, 0., modeL, modeM, 1 );
     }
     else if( !(*h_lms) ) // Create a new SphHarmTimeSeries to hold 1st mode
     {
       if (flag_fLower_extend == 1) // Find where 1st mode (2,2) hits fLower
       {
         cut_ind = find_instant_freq_hlm(sigMode, fLower, 1);
         sigMode = XLALResizeCOMPLEX16TimeSeries(sigMode, cut_ind,
             sigMode->data->length - cut_ind);
         if (!sigMode )
           XLAL_ERROR(XLAL_EFUNC);
       }
       *h_lms = XLALSphHarmTimeSeriesAddMode(sigMode, modeL, modeM, NULL);
     }
     else // Append additional modes into existing SphHarmTimeSeries
     {
       if (flag_fLower_extend == 1) // Resize all modes to length of (2,2)
       {
         sigMode = XLALResizeCOMPLEX16TimeSeries(sigMode, cut_ind,
             sigMode->data->length - cut_ind);
         if (!sigMode )
           XLAL_ERROR(XLAL_EFUNC);
       }
       *h_lms = XLALSphHarmTimeSeriesAddMode( sigMode, modeL, modeM, *h_lms);
     }

  } /* End loop over modes */

  /* clip the parts below f_min */
  if (flag_fLower_extend == 1)
  {
    if ( cos(inclination) < 0.0 )
    {
      cut_ind = find_instant_freq(*hplus, *hcross, fLower, 1, 1); 
    }
    else
    {
      cut_ind = find_instant_freq(*hplus, *hcross, fLower, 1, 0); 
    }
    *hplus = XLALResizeREAL8TimeSeries(*hplus, cut_ind, (*hplus)->data->length - cut_ind);
    *hcross = XLALResizeREAL8TimeSeries(*hcross, cut_ind, (*hcross)->data->length - cut_ind);
    if (!(*hplus) || !(*hcross))
      XLAL_ERROR(XLAL_EFUNC);
  }

  /* Clean up */
  XLALDestroyREAL8Vector( values );
  XLALDestroyREAL8Vector( dvalues );
  XLALDestroyREAL8Array( dynamics );
  XLALDestroyREAL8Array( dynamicsHi );
  XLALDestroyREAL8Vector ( sigReHi );
  XLALDestroyREAL8Vector ( sigImHi );
  XLALDestroyCOMPLEX16TimeSeries( sigMode );
  XLALDestroyREAL8Vector ( phseHi );
  XLALDestroyREAL8Vector ( omegaHi );
  XLALDestroyREAL8Vector( ampNQC );
  XLALDestroyREAL8Vector( q1 );
  XLALDestroyREAL8Vector( q2 );
  XLALDestroyREAL8Vector( q3 );
  XLALDestroyREAL8Vector( p1 );
  XLALDestroyREAL8Vector( p2 );

  return XLAL_SUCCESS;
}

/**
 * This function generates the plus and cross polarizations for the dominant
 * (2,2) mode of the EOBNRv2 approximant. This model is defined in Pan et al,
 * PRD84, 124052(2011).
 */
int
XLALSimIMREOBNRv2DominantMode(
              REAL8TimeSeries **hplus,      /**<< The +-polarization waveform (returned) */
              REAL8TimeSeries **hcross,     /**<< The x-polarization waveform (returned) */
              const REAL8       phiC,       /**<< The phase at the coalescence time (twice the orbital phase at the max orbital frequency moment) */
              const REAL8       deltaT,     /**<< Sampling interval (in seconds) */
              const REAL8       m1SI,       /**<< First component mass (in kg) */
              const REAL8       m2SI,       /**<< Second component mass (in kg) */
              const REAL8       fLower,     /**<< Starting frequency (in Hz) */
              const REAL8       distance,   /**<< Distance to source (in metres) */
              const REAL8       inclination /**<< Inclination of the source (in radians) */
              )
{

  if ( XLALSimIMREOBNRv2Generator(hplus, hcross, NULL, phiC, deltaT, m1SI, m2SI,
              fLower, distance, inclination, 0 ) == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  return XLAL_SUCCESS;
}

/**
 * This function generates the plus and cross polarizations for the EOBNRv2 approximant
 * with all available modes included. This model is defined in Pan et al,
 * PRD84, 124052(2011).
 */
int
XLALSimIMREOBNRv2AllModes(
              REAL8TimeSeries **hplus,      /**<< The +-polarization waveform (returned) */
              REAL8TimeSeries **hcross,     /**<< The x-polarization waveform (returned) */
              const REAL8       phiC,       /**<< The orbital phase at the coalescence time */
              const REAL8       deltaT,     /**<< Sampling interval (in seconds) */
              const REAL8       m1SI,       /**<< First component mass (in kg) */
              const REAL8       m2SI,       /**<< Second component mass (in kg) */
              const REAL8       fLower,     /**<< Starting frequency (in Hz) */
              const REAL8       distance,   /**<< Distance to source (in metres) */
              const REAL8       inclination /**<< Inclination of the source (in radians) */
              )
{

  if ( XLALSimIMREOBNRv2Generator(hplus, hcross, NULL, phiC, deltaT, m1SI, m2SI,
              fLower, distance, inclination, 1 ) == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  return XLAL_SUCCESS;
}

/**
 * Wrapper function to generate the -2 spin-weighted spherical harmonic modes
 * (as opposed to generating the polarizations). This model is defined in
 * Pan et al, PRD84, 124052(2011). Returns the (2,2), (2,1), (3,3), (4,4), (5,5)
 * SWSH modes in a SphHarmTimeSeries struct.
 */
SphHarmTimeSeries *XLALSimIMREOBNRv2Modes(
        const REAL8 phiRef,  /**< Orbital phase at coalescence (radians) */
        const REAL8 deltaT,  /**< Sampling interval (s) */
        const REAL8 m1,      /**< First component mass (kg) */
        const REAL8 m2,      /**< Second component mass (kg) */
        const REAL8 fLower,  /**< Starting GW frequency (Hz) */
        const REAL8 distance /**< Distance to sources (m) */
        )
{
  SphHarmTimeSeries *hlms = NULL;
  if ( XLALSimIMREOBNRv2Generator(NULL, NULL, &hlms, phiRef, deltaT, m1, m2,
              fLower, distance, 0., 1) == XLAL_FAILURE )
  {
    XLAL_ERROR_NULL( XLAL_EFUNC );
  }

  return hlms;
}
