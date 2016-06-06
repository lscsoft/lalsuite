#ifndef _LALSIM_IMR_PHENOMD_INTERNALS_H
#define _LALSIM_IMR_PHENOMD_INTERNALS_H

/*
 * Copyright (C) 2015 Michael Puerrer, Sebastian Khan, Frank Ohme, Ofek Birnholtz, Lionel London
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
 * \author Michael Puerrer, Sebastian Khan, Frank Ohme, Ofek Birnholtz, Lionel London
 *
 * \file
 *
 * \brief Internal function for IMRPhenomD phenomenological waveform model.
 * See \ref LALSimIMRPhenom_c for more details.
 *
 */

/*
This waveform uses the TaylorF2 coefficients for it's inspiral phase augmented
by higher order phenomenological terms tuned to SEOBv2-Hybrid waveforms.
Below are lines copied from LALSimInspiralPNCoefficients.c which are the TaylorF2
phase coefficients we have used.
We document them here in case changes to that file changes the behaviour
of this waveform.

    const REAL8 mtot = m1 + m2;
    const REAL8 d = (m1 - m2) / (m1 + m2);
    const REAL8 eta = m1*m2/mtot/mtot;
    const REAL8 m1M = m1/mtot;
    const REAL8 m2M = m2/mtot;
    // Use the spin-orbit variables from arXiv:1303.7412, Eq. 3.9
    // We write dSigmaL for their (\delta m/m) * \Sigma_\ell
    // There's a division by mtotal^2 in both the energy and flux terms
    // We just absorb the division by mtotal^2 into SL and dSigmaL

    const REAL8 SL = m1M*m1M*chi1L + m2M*m2M*chi2L;
    const REAL8 dSigmaL = d*(m2M*chi2L - m1M*chi1L);

    const REAL8 pfaN = 3.L/(128.L * eta);
    //Non-spin phasing terms - see arXiv:0907.0700, Eq. 3.18
    pfa->v[0] = 1.L;
    pfa->v[2] = 5.L*(743.L/84.L + 11.L * eta)/9.L;
    pfa->v[3] = -16.L*LAL_PI;
    pfa->v[4] = 5.L*(3058.673L/7.056L + 5429.L/7.L * eta
                     + 617.L * eta*eta)/72.L;
    pfa->v[5] = 5.L/9.L * (7729.L/84.L - 13.L * eta) * LAL_PI;
    pfa->vlogv[5] = 5.L/3.L * (7729.L/84.L - 13.L * eta) * LAL_PI;
    pfa->v[6] = (11583.231236531L/4.694215680L
                     - 640.L/3.L * LAL_PI * LAL_PI - 6848.L/21.L*LAL_GAMMA)
                 + eta * (-15737.765635L/3.048192L
                     + 2255./12. * LAL_PI * LAL_PI)
                 + eta*eta * 76055.L/1728.L
                 - eta*eta*eta * 127825.L/1296.L;
    pfa->v[6] += (-6848.L/21.L)*log(4.);
    pfa->vlogv[6] = -6848.L/21.L;
    pfa->v[7] = LAL_PI * ( 77096675.L/254016.L
                     + 378515.L/1512.L * eta - 74045.L/756.L * eta*eta);

    // Spin-orbit terms - can be derived from arXiv:1303.7412, Eq. 3.15-16
    const REAL8 pn_gamma = (554345.L/1134.L + 110.L*eta/9.L)*SL + (13915.L/84.L - 10.L*eta/3.)*dSigmaL;
    switch( spinO )
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
            pfa->v[7] += (-8980424995.L/762048.L + 6586595.L*eta/756.L - 305.L*eta*eta/36.L)*SL - (170978035.L/48384.L - 2876425.L*eta/672.L - 4735.L*eta*eta/144.L) * dSigmaL;
*/


#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <lal/LALStdlib.h>
#include <lal/LALSimIMR.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/LALSimInspiral.h>

#include "LALSimIMRPhenomD.h"

// NOTE: At the moment we have separate functions for each Phenom coefficient;
// these could be collected together

/**
  * Structure holding all coefficients for the amplitude
  */
typedef struct tagIMRPhenomDAmplitudeCoefficients {
  double eta;         // symmetric mass-ratio
  double chi1, chi2;  // dimensionless aligned spins, convention m1 >= m2.
  double q;           // asymmetric mass-ratio (q>=1)
  double chi;         // PN reduced spin parameter
  double fRD;         // ringdown frequency
  double fDM;         // imaginary part of the ringdown frequency (damping time)

  double fmaxCalc;    // frequency at which the mrerger-ringdown amplitude is maximum

  // Phenomenological inspiral amplitude coefficients
  double rho1;
  double rho2;
  double rho3;

  // Phenomenological intermediate amplitude coefficients
  double delta0;
  double delta1;
  double delta2;
  double delta3;
  double delta4;

  // Phenomenological merger-ringdown amplitude coefficients
  double gamma1;
  double gamma2;
  double gamma3;

  // Coefficients for collocation method. Used in intermediate amplitude model
  double f1, f2, f3;
  double v1, v2, v3;
  double d1, d2;

  // Transition frequencies for amplitude
  // We don't *have* to store them, but it may be clearer.
  double fInsJoin;    // Ins = Inspiral
  double fMRDJoin;    // MRD = Merger-Ringdown
}
IMRPhenomDAmplitudeCoefficients;

/**
  * Structure holding all coefficients for the phase
  */
typedef struct tagIMRPhenomDPhaseCoefficients {
  double eta;         // symmetric mass-ratio
  double chi1, chi2;  // dimensionless aligned spins, convention m1 >= m2.
  double q;           // asymmetric mass-ratio (q>=1)
  double chi;         // PN reduced spin parameter
  double fRD;         // ringdown frequency
  double fDM;         // imaginary part of the ringdown frequency (damping time)

  // Phenomenological inspiral phase coefficients
  double sigma1;
  double sigma2;
  double sigma3;
  double sigma4;
  double sigma5;

  // Phenomenological intermediate phase coefficients
  double beta1;
  double beta2;
  double beta3;

  // Phenomenological merger-ringdown phase coefficients
  double alpha1;
  double alpha2;
  double alpha3;
  double alpha4;
  double alpha5;

  // C1 phase connection coefficients
  double C1Int;
  double C2Int;
  double C1MRD;
  double C2MRD;

  // Transition frequencies for phase
  double fInsJoin;    // Ins = Inspiral
  double fMRDJoin;    // MRD = Merger-Ringdown
}
IMRPhenomDPhaseCoefficients;


 /**
   * Structure holding all additional coefficients needed for the delta amplitude functions.
   */
typedef struct tagdeltaUtility {
  double f12;
  double f13;
  double f14;
  double f15;
  double f22;
  double f23;
  double f24;
  double f32;
  double f33;
  double f34;
  double f35;
} DeltaUtility;

/*
 *
 * Internal function prototypes; f stands for geometric frequency "Mf"
 *
 */

////////////////////////////// Miscellaneous functions //////////////////////////////

static double chiPN(double eta, double chi1, double chi2);
static size_t NextPow2(const size_t n);
// static double StepFunc(const double t, const double t1);
static bool StepFunc_boolean(const double t, const double t1);

static inline double pow_2_of(double number);
static inline double pow_3_of(double number);
static inline double pow_4_of(double number);

static double Subtract3PNSS(double m1, double m2, double M, double chi1, double chi2);

/******************************* Constants to save floating-point pow calculations *******************************/

/**
 * useful powers in GW waveforms: 1/6, 1/3, 2/3, 4/3, 5/3, 7/3, 8/3
 * calculated using only one invocation of 'pow', the rest are just multiplications and divisions
 */
typedef struct tagUsefulPowers
{
    REAL8 sixth;
    REAL8 third;
    REAL8 two_thirds;
    REAL8 four_thirds;
    REAL8 five_thirds;
	REAL8 two;
    REAL8 seven_thirds;
    REAL8 eight_thirds;
} UsefulPowers;

/**
 * must be called before the first usage of *p
 */
static int init_useful_powers(UsefulPowers * p, REAL8 number);

/**
 * useful powers of LAL_PI, calculated once and kept constant - to be initied with a call to
 * init_useful_powers(&powers_of_pi, LAL_PI);
 *
 * only declared here, defined in LALSIMIMRPhenomD.c (because this c file is "included" like an h file)
 */
extern UsefulPowers powers_of_pi;

/**
 * used to cache the recurring (frequency-independant) prefactors of AmpInsAnsatz. Must be inited with a call to
 * init_amp_ins_prefactors(&prefactors, p);
 */
typedef struct tagAmpInsPrefactors
{
	double two_thirds;
	double one;
	double four_thirds;
	double five_thirds;
	double two;
	double seven_thirds;
	double eight_thirds;
	double three;

	double amp0;
} AmpInsPrefactors;

/**
 * must be called before the first usage of *prefactors
 */
static int init_amp_ins_prefactors(AmpInsPrefactors * prefactors, IMRPhenomDAmplitudeCoefficients* p);

/**
 * used to cache the recurring (frequency-independant) prefactors of PhiInsAnsatzInt. Must be inited with a call to
 * init_phi_ins_prefactors(&prefactors, p, pn);
 */
typedef struct tagPhiInsPrefactors
{
	double initial_phasing;
	double third;
	double third_with_logv;
	double two_thirds;
	double one;
	double four_thirds;
	double five_thirds;
	double two;
	double logv;
	double minus_third;
	double minus_two_thirds;
	double minus_one;
	double minus_five_thirds;
} PhiInsPrefactors;

/**
 * must be called before the first usage of *prefactors
 */
static int init_phi_ins_prefactors(PhiInsPrefactors * prefactors, IMRPhenomDPhaseCoefficients* p, PNPhasingSeries *pn);


/******************************* integer powers floating-point pow calculations *******************************/

/**
 * calc square of number without floating point 'pow'
 */
static inline double pow_2_of(double number)
{
	return (number*number);
}

/**
 * calc cube of number without floating point 'pow'
 */
static inline double pow_3_of(double number)
{
	return (number*number*number);
}

/**
 * calc fourth power of number without floating point 'pow'
 */
static inline double pow_4_of(double number)
{
	double pow2 = pow_2_of(number);
	return pow2 * pow2;
}

//////////////////////// Final spin, final mass, fring, fdamp ///////////////////////

static double FinalSpin0815_s(double eta, double s);
UNUSED static double FinalSpin0815(double eta, double chi1, double chi2);
static double EradRational0815_s(double eta, double s);
static double EradRational0815(double eta, double chi1, double chi2);
static double fring(double eta, double chi1, double chi2, double finalspin);
static double fdamp(double eta, double chi1, double chi2, double finalspin);

/******************************* Amplitude functions *******************************/

static double amp0Func(double eta);

///////////////////////////// Amplitude: Inspiral functions /////////////////////////

static double rho1_fun(double eta, double chiPN);
static double rho2_fun(double eta, double chiPN);
static double rho3_fun(double eta, double chiPN);
static double AmpInsAnsatz(double Mf, UsefulPowers * powers_of_Mf, AmpInsPrefactors * prefactors);
static double DAmpInsAnsatz(double Mf, IMRPhenomDAmplitudeCoefficients* p);

////////////////////////// Amplitude: Merger-Ringdown functions //////////////////////

static double gamma1_fun(double eta, double chiPN);
static double gamma2_fun(double eta, double chiPN);
static double gamma3_fun(double eta, double chiPN);
static double AmpMRDAnsatz(double f, IMRPhenomDAmplitudeCoefficients* p);
static double DAmpMRDAnsatz(double f, IMRPhenomDAmplitudeCoefficients* p);
static double fmaxCalc(IMRPhenomDAmplitudeCoefficients* p);

//////////////////////////// Amplitude: Intermediate functions ///////////////////////

static double AmpIntAnsatz(double f, IMRPhenomDAmplitudeCoefficients* p);
static double AmpIntColFitCoeff(double eta, double chiPN); //this is the v2 value
static double delta0_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d);
static double delta1_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d);
static double delta2_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d);
static double delta3_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d);
static double delta4_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d);
static void ComputeDeltasFromCollocation(IMRPhenomDAmplitudeCoefficients* p);

///////////////////////////// Amplitude: glueing function ////////////////////////////

static IMRPhenomDAmplitudeCoefficients* ComputeIMRPhenomDAmplitudeCoefficients(double eta, double chi1, double chi2, double finspin);
static double IMRPhenDAmplitude(double f, IMRPhenomDAmplitudeCoefficients *p, UsefulPowers *powers_of_f, AmpInsPrefactors * prefactors);

/********************************* Phase functions *********************************/

/////////////////////////////// Phase: Ringdown functions ////////////////////////////

static double alpha1Fit(double eta, double chiPN);
static double alpha2Fit(double eta, double chiPN);
static double alpha3Fit(double eta, double chiPN);
static double alpha4Fit(double eta, double chiPN);
static double alpha5Fit(double eta, double chiPN);
static double PhiMRDAnsatzInt(double f, IMRPhenomDPhaseCoefficients *p);
static double DPhiMRD(double f, IMRPhenomDPhaseCoefficients *p);

/////////////////////////// Phase: Intermediate functions ///////////////////////////

static double beta1Fit(double eta, double chiPN);
static double beta2Fit(double eta, double chiPN);
static double beta3Fit(double eta, double chiPN);
static double PhiIntAnsatz(double f, IMRPhenomDPhaseCoefficients *p);
static double DPhiIntAnsatz(double f, IMRPhenomDPhaseCoefficients *p);
static double DPhiIntTemp(double ff, IMRPhenomDPhaseCoefficients *p);

///////////////////////////// Phase: Inspiral functions /////////////////////////////

static double sigma1Fit(double eta, double chiPN);
static double sigma2Fit(double eta, double chiPN);
static double sigma3Fit(double eta, double chiPN);
static double sigma4Fit(double eta, double chiPN);
static double PhiInsAnsatzInt(double f, UsefulPowers * powers_of_Mf, PhiInsPrefactors * prefactors, IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn);
static double DPhiInsAnsatzInt(double ff, IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn);

////////////////////////////// Phase: glueing function //////////////////////////////

static IMRPhenomDPhaseCoefficients* ComputeIMRPhenomDPhaseCoefficients(double eta, double chi1, double chi2, double finspin, const LALSimInspiralTestGRParam *extraParams);
static void ComputeIMRPhenDPhaseConnectionCoefficients(IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn, PhiInsPrefactors * prefactors);
static double IMRPhenDPhase(double f, IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn, UsefulPowers *powers_of_f, PhiInsPrefactors * prefactors);

#endif	// of #ifndef _LALSIM_IMR_PHENOMD_INTERNALS_H
