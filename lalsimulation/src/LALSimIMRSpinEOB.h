#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMREOBNRv2.h"

#ifndef _LALSIMIMRSPINEOB_H
#define _LALSIMIMRSPINEOB_H

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * Parameters for the spinning EOB model, used in calculating the Hamiltonian.
 * The Hamiltonian is given in Barausse and Buonanno (http://arxiv.org/pdf/0912.3517)
 * The parameters correspond to the following as found in the paper:
 * KK - K found in the equations for \Delta_is (Eqn 5.77-5.81)
 * k0 - \Delta_0 given in Eqn 5.77
 * k1 - \Delta_1 given in Eqn 5.78
 * k2 - \Delta_2 given in Eqn 5.79
 * k3 - \Delta_3 given in Eqn 5.80
 * k4 - \Delta_4 given in Eqn 5.81
 * b3 - \omega^{fd}_2 given in Eqn 5.40
 * bb3 - \omega^{fd}_1 given in Eqn 5.40
 */

typedef struct
tagSpinEOBHCoeffs
{
  double KK;
  double k0;
  double k1;
  double k2;
  double k3;
  double k4;
  double b3;
  double bb3;
}
SpinEOBHCoeffs;

typedef struct
tagSpinEOBParams
{
  EOBParams               *eobParams;
  SpinEOBHCoeffs          *seobCoeffs;
  REAL8Vector             *sigmaStar;
  REAL8Vector             *sigmaKerr;
  REAL8                   a;
  int                     alignedSpins;
  int                     tortoise;
}
SpinEOBParams;

/* We need to encapsulate the data for the GSL derivative function */
typedef
struct tagHcapDerivParams
{
   const REAL8   *values;
   SpinEOBParams *params;
   UINT4         varyParam;
}
HcapDerivParams;

/* We need to encapsulate the data for calculating spherical 2nd derivatives */
typedef
struct tagHcapSphDeriv2Params
{
  const REAL8     *sphValues;
  SpinEOBParams   *params;
  UINT4           varyParam1;
  UINT4           varyParam2;
}
HcapSphDeriv2Params;

/**
 * This function calculates the function \f$\Delta_r(r)\f$ which appears in the spinning EOB
 * potential function.
 */
REAL8 XLALSimIMRSpinEOBHamiltonianDeltaR(
        SpinEOBHCoeffs *coeffs, /**<< Pre-computed coefficients which appear in the function */
        const REAL8    r,       /**<< Current orbital radius (in units of total mass) */
        const REAL8    eta,     /**<< Symmetric mass ratio */
        const REAL8    a        /**<< Normalized deformed Kerr spin */
        );

REAL8 XLALSimIMRSpinEOBHamiltonian(
               const REAL8    eta,
               REAL8Vector    * restrict x,
               REAL8Vector    * restrict p,
               REAL8Vector    * restrict sigmaKerr,
               REAL8Vector    * restrict sigmaStar,
               int                       tortoise,
               SpinEOBHCoeffs *coeffs);

int XLALCalculateSpinEOBCoeffs(
        SpinEOBHCoeffs *params,
        const REAL8    eta,
        const REAL8    a
        );

int XLALSpinAlignedHcapDerivative(
                          double                t,
                          const REAL8           values[],
                          REAL8                 dvalues[],
                          void                  *funcParams
                               );

REAL8 XLALSimIMRSpinAlignedEOBCalcOmega(
                      const REAL8          values[],
                      SpinEOBParams        *funcParams
                      );

REAL8 XLALSimIMRSpinAlignedEOBNonKeplerCoeff(
                      const REAL8           values[],
                      SpinEOBParams         *funcParams
                      );

int XLALSpinHcapNumericalDerivative(
                          double                t,
                          const REAL8           values[],
                          REAL8                 dvalues[],
                          void                  *funcParams
                               );


REAL8 XLALSpinHcapNumDerivWRTParam(
                       const INT4 paramIdx,
                       const REAL8 values[],
                       SpinEOBParams *params
                       );

int XLALSimIMRSpinEOBCalculateSigmaKerr( REAL8Vector *sigmaKerr,
                                   REAL8        mass1,
                                   REAL8        mass2,
                                   REAL8Vector *s1,
                                   REAL8Vector *s2 );

int XLALSimIMRSpinEOBCalculateSigmaStar( REAL8Vector *sigmaStar,
                                   REAL8        mass1,
                                   REAL8        mass2,
                                   REAL8Vector *s1,
                                   REAL8Vector *s2 );

int XLALSimIMRCalculateSpinEOBHCoeffs(
        SpinEOBHCoeffs *coeffs,
        const REAL8    eta,
        const REAL8    a
        );

int XLALSimIMRSpinEOBWaveform(
        REAL8TimeSeries **hplus,
        REAL8TimeSeries **hcross,
        LIGOTimeGPS     *tc,
        const REAL8     phiC,
        const REAL8     deltaT,
        const REAL8     m1,
        const REAL8     m2,
        const REAL8     fMin,
        const REAL8     r,
        const REAL8     inc,
        const REAL8     spin1[],
        const REAL8     spin2[]
     );

int XLALSimIMRSpinAlignedEOBWaveform(
        REAL8TimeSeries **hplus,
        REAL8TimeSeries **hcross,
        LIGOTimeGPS     *tc,
        const REAL8     phiC,
        const REAL8     deltaT,
        const REAL8     m1,
        const REAL8     m2,
        const REAL8     fMin,
        const REAL8     r,
        const REAL8     inc,
        const REAL8     spin1[],
        const REAL8     spin2[]
     );

REAL8 XLALInspiralSpinFactorizedFlux(
                      REAL8Vector           *values,
                      const REAL8           omega,
                      SpinEOBParams         *ak,
                      const REAL8            H,
                      const INT4             lMax
                     );

INT4 XLALSimIMRSpinEOBGetSpinFactorizedWaveform( COMPLEX16             * restrict hlm,
                                REAL8Vector           * restrict values,
                                const REAL8           v,
                                const REAL8           Hreal,
                                const INT4            l,
                                const INT4            m,
                                SpinEOBParams         * restrict params
                                );

int XLALSimIMRSpinEOBInitialConditions(
                      REAL8Vector *initConds,
                      const REAL8  mass1,
                      const REAL8  mass2,
                      const REAL8  fMin,
                      const REAL8  inc,
                      const REAL8  spin1[],
                      const REAL8  spin2[],
                      SpinEOBParams *params
                      );

#endif /* _LALSIMIMRSPINEOB_H */
