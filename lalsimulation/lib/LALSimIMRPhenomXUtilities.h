#ifndef _LALSIM_IMR_PHENOMX_UTILITIES_H
#define _LALSIM_IMR_PHENOMX_UTILITIES_H

/*
 * Copyright (C) 2018 Geraint Pratten
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */


/*
 * \author Geraint Pratten
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/LALSimInspiral.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
  //#include <complex.h>

/* ********************** NUMERICAL UTILITY FUNCTIONS ********************* */
size_t NextPow2(const size_t n);
bool IMRPhenomX_StepFuncBool(const double t, const double t1);
void IMRPhenomX_InternalNudge(REAL8 x, REAL8 X, REAL8 epsilon);
bool IMRPhenomX_ApproxEqual(REAL8 x, REAL8 y, REAL8 epsilon);

double pow_2_of(double number);
double pow_3_of(double number);
double pow_4_of(double number);
double pow_5_of(double number);
double pow_6_of(double number);
double pow_7_of(double number);
double pow_8_of(double number);
double pow_9_of(double number);

REAL8 XLALSimIMRPhenomXsign(REAL8 x);

REAL8 XLALSimIMRPhenomXatan2tol(REAL8 a, REAL8 b, REAL8 tol);


/* We want to expose the functions below - they are useful */
/* ********************** FREQUENCY CONVERSIONS ********************* */
REAL8 XLALSimIMRPhenomXUtilsMftoHz(REAL8 Mf, REAL8 Mtot_Msun);
REAL8 XLALSimIMRPhenomXUtilsHztoMf(REAL8 fHz, REAL8 Mtot_Msun);

REAL8 XLALSimIMRPhenomXLina(REAL8 eta, REAL8 S, REAL8 dchi, REAL8 delta);
REAL8 XLALSimIMRPhenomXLinb(REAL8 eta, REAL8 S, REAL8 dchi, REAL8 delta);
REAL8 XLALSimIMRPhenomXPsi4ToStrain(double eta, double S, double dchi);

/* ********************** MECO, ISCO, ETC ********************* */
REAL8 XLALSimIMRPhenomXfMECO(REAL8 eta, REAL8 chi1L, REAL8 chi2L);
REAL8 XLALSimIMRPhenomXfISCO(REAL8 chif);

/* ********************** SPIN PARAMETERISATIONS ********************* */
REAL8 XLALSimIMRPhenomXchiPN(REAL8 eta, REAL8 chi1l, REAL8 chi2l);
REAL8 XLALSimIMRPhenomXchiPNHat(REAL8 eta, REAL8 chi1l, REAL8 chi2l);
REAL8 XLALSimIMRPhenomXchiEff(REAL8 eta, REAL8 chi1l, REAL8 chi2l);
REAL8 XLALSimIMRPhenomXdchi(REAL8 chi1l, REAL8 chi2l);
REAL8 XLALSimIMRPhenomXSTotR(REAL8 eta, REAL8 chi1l, REAL8 chi2l);

/* ********************** MASS PARAMETERISATIONS ********************* */

/* ********************** FINAL STATE ********************* */
REAL8 XLALSimIMRPhenomXFinalSpin2017(REAL8 eta, REAL8 chi1L, REAL8 chi2L);
REAL8 XLALSimIMRPhenomXFinalMass2017(REAL8 eta, REAL8 chi1L, REAL8 chi2L);
REAL8 XLALSimIMRPhenomXErad2017(REAL8 eta, REAL8 chi1L, REAL8 chi2L);
REAL8 XLALSimIMRPhenomXPrecessingFinalSpin2017(REAL8 eta, REAL8 chi1L, REAL8 chi2L, REAL8 chi_inplane);

/* Check masses and spins */
INT4 XLALIMRPhenomXPCheckMassesAndSpins(REAL8 *m1, REAL8 *m2, REAL8 *chi1x, REAL8 *chi1y, REAL8 *chi1z, REAL8 *chi2x, REAL8 *chi2y, REAL8 *chi2z);

/* ********************** ANALYTICAL MODEL WRAPPERS *********************  */
REAL8 XLALSimIMRPhenomXIntermediatePhase22AnsatzAnalytical(REAL8 ff, REAL8 fRD, REAL8 fDA, REAL8 a0, REAL8 a1, REAL8 a2, REAL8 a3, REAL8 a4, REAL8 aL);
REAL8 XLALSimIMRPhenomXIntermediateAmplitude22AnsatzAnalytical(REAL8 ff, REAL8 ff7o6, REAL8 a0, REAL8 a1, REAL8 a2, REAL8 a3, REAL8 a4, REAL8 a5);
REAL8 XLALSimIMRPhenomXRingdownPhase22AnsatzAnalytical(REAL8 ff, REAL8 fRD, REAL8 fDA, REAL8 a0, REAL8 a1, REAL8 a2, REAL8 a4, REAL8 aL);
REAL8 XLALSimIMRPhenomXRingdownPhaseDeriv22AnsatzAnalytical(REAL8 ff, REAL8 fRD, REAL8 fDA, REAL8 a0, REAL8 a1, REAL8 a2, REAL8 a4, REAL8 aL);
REAL8 XLALSimIMRPhenomXRingdownAmplitude22AnsatzAnalytical(REAL8 ff, REAL8 fRD, REAL8 fDA, REAL8 gamma1, REAL8 gamma2, REAL8 gamma3);

REAL8 XLALSimIMRPhenomXAmp22Prefactor(REAL8 eta);




/*
  The functions below are XLAL exposed of the QNM ringdown and damping frequency used for the
  IMRPhenomX model: https://arxiv.org/abs/2001.11412.

  See:
  https://arxiv.org/src/2001.10914v1/anc/QNMs/CoefficientStatsfring22.m
  https://arxiv.org/src/2001.10914v1/anc/QNMs/CoefficientStatsfdamp22.m
*/
REAL8 XLALSimIMRPhenomXfring22(const REAL8 af);
REAL8 XLALSimIMRPhenomXfdamp22(const REAL8 af);




#ifdef __cplusplus
}
#endif

#endif /* _LALSIM_IMR_PHENOMX_UTILITIES_H */
