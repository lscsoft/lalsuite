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

#endif /* _LALSIMIMRSPINEOB_H */
