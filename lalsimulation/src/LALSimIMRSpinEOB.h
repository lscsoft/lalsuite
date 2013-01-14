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
 * KK - K found in the equations for \f$\Delta_i\f$ (Eqn 5.77-5.81)
 * k0 - \f$\Delta_0\f$ given in Eqn 5.77
 * k1 - \f$\Delta_1\f$ given in Eqn 5.78
 * k2 - \f$\Delta_2\f$ given in Eqn 5.79
 * k3 - \f$\Delta_3\f$ given in Eqn 5.80
 * k4 - \f$\Delta_4\f$ given in Eqn 5.81
 * b3 - \f$\omega^{fd}_2\f$ given in Eqn 5.40
 * bb3 - \f$\omega^{fd}_1\f$ given in Eqn 5.40
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

/**
 * Parameters for the spinning EOB model.
 * 1) eobParams contains parameters common to nonspin and spin EOBNR models,
 *    including mass ratio, masses, pre-computed coefficients for potential, flux and waveforms,
 *    NQC coefficients and Newtonian multiple prefixes.
 * 2) seobCoeffs contans parameters for calculating the spin EOB Hamiltonian.
 * 3) sigmaStar and sigmaKerr are effective spins of the test-particle and background.
 * 4) a is the spin value being used for test-particle limit spin terms.
 * 5) alignedSpins and tortoise are controling flags.
 */

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

#endif /* _LALSIMIMRSPINEOB_H */
