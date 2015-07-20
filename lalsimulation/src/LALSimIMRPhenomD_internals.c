/*
 * Copyright (C) 2015 Michael Puerrer, Sebastian Khan, Frank Ohme
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
 * \author Michael Puerrer, Sebastian Khan, Frank Ohme
 *
 * \file
 *
 * \brief C code for IMRPhenomD phenomenological waveform model.
 * See ... for details.
 *
 * This is an aligned-spin frequency domain model.
 *
 * @note The model was calibrated to mass-ratios [1:1,1:4,1:8,1:18].
 * * Along the mass-ratio 1:1 line it was calibrated to spins  [-0.95, +0.98].
 * * Along the mass-ratio 1:4 line it was calibrated to spins  [-0.75, +0.75].
 * * Along the mass-ratio 1:8 line it was calibrated to spins  [-0.85, +0.85].
 * * Along the mass-ratio 1:18 line it was calibrated to spins [-0.8, +0.4].
 * The calibration points will be given in forthcoming papers.
 *
 * @note The model is usable outside this parameter range,
 * and in tests to date gives sensible physical results,
 * but conclusive statements on the physical fidelity of
 * the model for these parameters await comparisons against further
 * numerical-relativity simulations.
 *
 */

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

// Constants in Mathematica CForm expressions
const double Pi = LAL_PI;
const double EulerGamma = LAL_GAMMA;

// NOTE: At the moment we have separate functions for each Phenom coefficient;
// these could be collected together

typedef struct tagIMRPhenomDAmplitudeCoefficients {
  double eta;         // symmetric mass-ratio
  double chi1, chi2;  // dimensionless aligned spins, convention m1 >= m2.
  double q;           // asymmetric mass-ratio (q>=1)
  double chi;         // PN reduced spin parameter
  double fRD;         // ringdown frequency
  double fDM;         // imaginary part of the ringdown frequency (damping time)

  double fmaxCalc;

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

  // Coefficients for collocation method
  double f1, f2, f3;
  double v1, v2, v3;
  double d1, d2;

  // Transition frequencies for amplitude
  // We don't *have* to store them, but it may be clearer.
  double fInsJoin;
  double fMRDJoin;
}
IMRPhenomDAmplitudeCoefficients;

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
  double fInsJoin;
  double fMRDJoin;
}
IMRPhenomDPhaseCoefficients;

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

/**
 *
 * Internal function prototypes; f stands for geometric frequency "Mf"
 *
 */

////////////////////////////// Miscellaneous functions //////////////////////////////

static double chiPN(double eta, double chi1, double chi2);
static double PlanckTaper(const double t, const double t1, const double t2);
static size_t NextPow2(const size_t n);

//////////////////////// Final spin, final mass, fring, fdamp ///////////////////////

static double FinalSpin0714_s(double eta, double s);
static double FinalSpin0714(double eta, double chi1, double chi2);
static double EradRational_s(double eta, double s);
static double EradRational(double eta, double chi1, double chi2);
static double fring(double eta, double chi1, double chi2);
static double fdamp(double eta, double chi1, double chi2);

/******************************* Amplitude functions *******************************/

static double amp0Func(double eta);

///////////////////////////// Amplitude: Inspiral functions /////////////////////////

static double rho1_fun(double eta, double chiPN);
static double rho2_fun(double eta, double chiPN);
static double rho3_fun(double eta, double chiPN);
static double AmpInsAnsatz(double f, IMRPhenomDAmplitudeCoefficients* p);
static double DAmpInsAnsatz(double f, IMRPhenomDAmplitudeCoefficients* p);

////////////////////////// Amplitude: Merger-Ringdown functions //////////////////////

static double gamma1_fun(double eta, double chiPN);
static double gamma2_fun(double eta, double chiPN);
static double gamma3_fun(double eta, double chiPN);
static double AmpMRDAnsatz(double f, IMRPhenomDAmplitudeCoefficients* p);
static double DAmpMRDAnsatz(double f, IMRPhenomDAmplitudeCoefficients* p);
static double fmaxCalc(IMRPhenomDAmplitudeCoefficients* p);

//////////////////////////// Amplitude: Intermediate functions ///////////////////////

static double AmpIntAnsatz(double f, IMRPhenomDAmplitudeCoefficients* p);
static double AmpIntColFitCoeff(double eta, double chiPN);
static double delta0_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d);
static double delta1_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d);
static double delta2_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d);
static double delta3_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d);
static double delta4_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d);
static void ComputeDeltasFromCollocation(IMRPhenomDAmplitudeCoefficients* p);

///////////////////////////// Amplitude: glueing function ////////////////////////////

static IMRPhenomDAmplitudeCoefficients* ComputeIMRPhenomDAmplitudeCoefficients(double eta, double chi1, double chi2);
static double IMRPhenDAmplitude(double f, IMRPhenomDAmplitudeCoefficients *p);

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
static double PhiInsAnsatzInt(double f, IMRPhenomDPhaseCoefficients *p);
static double DPhiInsAnsatzInt(double ff, IMRPhenomDPhaseCoefficients *p);

////////////////////////////// Phase: glueing function //////////////////////////////

static IMRPhenomDPhaseCoefficients* ComputeIMRPhenomDPhaseCoefficients(double eta, double chi1, double chi2);
static void ComputeIMRPhenDPhaseConnectionCoefficients(IMRPhenomDPhaseCoefficients *p);
static double IMRPhenDPhase(double f, IMRPhenomDPhaseCoefficients *p);

/**
 *
 * Internal function definitions
 *
 * */

////////////////////////////// Miscellaneous functions //////////////////////////////

// PN reduced spin parameter
// See Eq 5.9 in http://arxiv.org/pdf/1107.1267v2.pdf
static double chiPN(double eta, double chi1, double chi2) {
  // Convention m1 >= m2
  double q = (1.0 + sqrt(1.0 - 4.0*eta) - 2.0*eta) / (2.0*eta);
  double M = 1; // only used for delta; value is irrelevant
  double m1 = M*q/(1.0+q);
  double m2 = M*1.0/(1.0+q);
  double chi_s = (chi1 + chi2) / 2.0;
  double chi_a = (chi1 - chi2) / 2.0;
  double delta = (m1 - m2) / M;
  return chi_s * (1.0 - 76.0/113.0*eta) + delta*chi_a;
}

// Planck taper function. See http://arxiv.org/abs/1003.2939
static double PlanckTaper(const double t, const double t1, const double t2) {
  if (t <= t1)
    return 0.0;
  else if (t >= t2)
    return 1.0;
  else
    return 1.0 / (exp((t2 - t1)/(t - t1) + (t2 - t1)/(t - t2)) + 1.0);
}

// Return the closest higher power of 2
static size_t NextPow2(const size_t n) {
  return 1 << (size_t) ceil(log2(n));
}

//////////////////////// Final spin, final mass, fring, fdamp ////////////////////////

// TODO: Add documentation of these functions

static double FinalSpin0714_s(double eta, double s) {
  double eta2 = eta*eta;
  double eta3 = eta2*eta;
  double s2 = s*s;
  double s3 = s2*s;
  double s4 = s3*s;

  return 3.4641016151377544*eta - 3.896044039898422*eta2 + 4.133657521035006*eta3 +
   (1 - 2.6507474704915883*eta + 0.5478819823269401*eta2 - 4.032593676598959*eta3)*s +
   (0.12443985124277726*eta - 0.9295544432946068*eta2 - 0.10106384443604632*eta3)*s2 +
   (-0.5142453933108723*eta + 1.8929952741242566*eta2)*s3 + (-0.6162184851271666*eta + 2.384721226741833*eta2)*s4;
}

static double FinalSpin0714(double eta, double chi1, double chi2) {
  // Convention m1 >= m2
  double m1 = 0.5 * (1.0 + sqrt(1.0 - 4.0*eta));
  double m2 = 0.5 * (1.0 - sqrt(1.0 - 4.0*eta));
  double m1s = m1*m1;
  double m2s = m2*m2;
  double s = (m1s * chi1 + m2s * chi2) / (m1s + m2s); // CHECK: this is different in EradRational: Is this correct?

  return FinalSpin0714_s(eta, s);
}

static double EradRational_s(double eta, double s) {
  double eta2 = eta*eta;
  double eta3 = eta2*eta;

  return ((0.0731529149096712*eta + 0.26674934921953924*eta2 + 0.8576638056091429*eta3)*
     (1. + (0.7306635209598181 - 5.097805627396959*eta + 7.881058808305713*eta2)*s))/(1. + (-0.15735564839807778 - 4.051664713223574*eta)*s);
}

static double EradRational(double eta, double chi1, double chi2) {
  // Convention m1 >= m2
  double m1 = 0.5 * (1.0 + sqrt(1.0 - 4.0*eta));
  double m2 = 0.5 * (1.0 - sqrt(1.0 - 4.0*eta));
  double m1s = m1*m1;
  double m2s = m2*m2;
  double s = m1s * chi1 + m2s * chi2; // CHECK: this is different in FinalSpin0714: Is this correct?

  return EradRational_s(eta, s);
}

static double fring(double eta, double chi1, double chi2) {
  double return_val;

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *iFring = gsl_spline_alloc(gsl_interp_cspline, QNMData_length);
  gsl_spline_init(iFring, QNMData_a, QNMData_fring, QNMData_length);

  return_val = gsl_spline_eval(iFring, FinalSpin0714(eta, chi1, chi2), acc) / (1.0 - EradRational(eta, chi1, chi2));

  gsl_spline_free(iFring);
  gsl_interp_accel_free(acc);
  return return_val;
}

static double fdamp(double eta, double chi1, double chi2) {
  double return_val;
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *iFdamp = gsl_spline_alloc(gsl_interp_cspline, QNMData_length);
  gsl_spline_init(iFdamp, QNMData_a, QNMData_fdamp, QNMData_length);

  return_val = gsl_spline_eval(iFdamp, FinalSpin0714(eta, chi1, chi2), acc) / (1.0 - EradRational(eta, chi1, chi2));

  gsl_spline_free(iFdamp);
  gsl_interp_accel_free(acc);
  return return_val;
}

/******************************* Amplitude functions *******************************/

static double amp0Func(double eta) {
  return (sqrt(0.6666666666666666)*sqrt(eta))/pow(Pi,0.16666666666666666);
}

///////////////////////////// Amplitude: Inspiral functions /////////////////////////

// Phenom coefficients rho1, ..., rho3 from direct fit
// AmpInsDFFitCoeffChiPNFunc[eta, chiPN]

static double rho1_fun(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 9971.492472945103 - 17037.888891610917*eta
    + (12668.61511289507 + 346089.83163884433*eta - 1.2149167739883522e6*eta2)*xi
    + (-67553.74726851289 + 1.3861298880085826e6*eta - 3.963246463568903e6*eta2)*xi2
    + (-60370.53480844165 + 804019.1829911621*eta - 2.0902056443125196e6*eta2)*xi3;
}

static double rho2_fun(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return -66310.52250429206 + 109014.21151257299*eta
    + (-17171.567397449413 - 3.490892198148065e6*eta + 1.1374537779220121e7*eta2)*xi
    + (741885.6177097048 - 1.3100372065215468e7*eta + 3.643880366528771e7*eta2)*xi2
    + (598326.6388874187 - 7.43208895148914e6*eta + 1.8925268527639613e7*eta2)*xi3;
}

static double rho3_fun(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 121291.86682475785 - 185262.388853303*eta
    + (-152202.97634608237 + 8.743489674196277e6*eta - 2.6916230268901654e7*eta2)*xi
    + (-1.971833989372485e6 + 3.091520754702774e7*eta - 8.390584129432291e7*eta2)*xi2
    + (-1.4569827161973754e6 + 1.7071365122228015e7*eta - 4.2745390843862176e7*eta2)*xi3;
}

// The Newtonian term in LAL is fine and we should use exactly the same (either hardcoded or call).
// The higher order amplitude corrections in LAL are wrong.
// So, we just use the Mathematica expression for convencience.
static double AmpInsAnsatz(double Mf, IMRPhenomDAmplitudeCoefficients* p) {
  double eta = p->eta;
  double chi1 = p->chi1;
  double chi2 = p->chi2;
  double rho1 = p->rho1;
  double rho2 = p->rho2;
  double rho3 = p->rho3;

  double chi12 = chi1*chi1;
  double chi13 = chi12*chi1;
  double chi22 = chi2*chi2;
  double eta2 = eta*eta;
  double eta3 = eta*eta2;
  double Mf2 = Mf*Mf;
  double Mf3 = Mf*Mf2;
  double Pi2 = Pi*Pi;
  double Seta = sqrt(1.0 - 4.0*eta);

  // optimized expression
  return 1 + ((-969 + 1804*eta)*pow(Pi*Mf,0.6666666666666666))/672.
  + ((chi1*(81*(1 + Seta) - 44*eta) + chi2*(81 - 81*Seta - 44*eta))*Mf*Pi)/48.
  + ((-27312085 - 10287648*chi22 - 10287648*chi12*(1 + Seta) + 10287648*chi22*Seta
  + 24*(-1975055 + 857304*chi12 - 994896*chi1*chi2 + 857304*chi22)*eta + 35371056*eta2)
  *pow(Pi*Mf,1.3333333333333333))/8.128512e6
  + (pow(Pi*Mf,1.6666666666666667)*(-6048*chi13*(-1 - Seta + (3 + Seta)*eta)
  + chi1*(287213*(1 + Seta) - 4*(93414 + 2083*Seta)*eta - 35632*eta2)
  + chi2*(-((287213 + 6048*chi22)*(-1 + Seta)) + 4*(-93414 + 1512*chi22*(-3 + Seta) + 2083*Seta)*eta
  - 35632*eta2) + 42840*(-1 + 4*eta)*Pi))/32256.
  - (Mf2*Pi2*(-336*(-3248849057 + 1809550512*chi12 - 2954929824*chi1*chi2 + 1809550512*chi22)*eta2
  - 324322727232*eta3 + 7*(177520268561 + 29362199328*chi22 + 29362199328*chi12*(1 + Seta)
  - 29362199328*chi22*Seta + 12160253952*(chi1 + chi2 + chi1*Seta - chi2*Seta)*Pi)
  + 12*eta*(-545384828789 + 49568837472*chi1*chi2 - 12312458928*chi22
  + 77616*chi12*(-158633 + 282718*Seta) - 21943440288*chi22*Seta - 8345272320*(chi1 + chi2)*Pi
  + 21384760320*Pi2)))/6.0085960704e10 + pow(Mf,2.333333333333333)*rho1
  + pow(Mf,2.6666666666666665)*rho2 + Mf3*rho3;
}

static double DAmpInsAnsatz(double f, IMRPhenomDAmplitudeCoefficients* p) {
  double eta = p->eta;
  double chi1 = p->chi1;
  double chi2 = p->chi2;
  double rho1 = p->rho1;
  double rho2 = p->rho2;
  double rho3 = p->rho3;

  return ((-969 + 1804*eta)*pow(Pi,0.6666666666666666))/(1008.*pow(f,0.3333333333333333)) +
   ((chi1*(81*(1 + sqrt(1 - 4*eta)) - 44*eta) + chi2*(81 - 81*sqrt(1 - 4*eta) - 44*eta))*Pi)/48. +
   ((-27312085 - 10287648*pow(chi2,2) - 10287648*pow(chi1,2)*(1 + sqrt(1 - 4*eta)) + 10287648*pow(chi2,2)*sqrt(1 - 4*eta) +
        24*(-1975055 + 857304*pow(chi1,2) - 994896*chi1*chi2 + 857304*pow(chi2,2))*eta + 35371056*pow(eta,2))*pow(f,0.3333333333333333)*
      pow(Pi,1.3333333333333333))/6.096384e6 + (5*pow(f,0.6666666666666666)*pow(Pi,1.6666666666666667)*
      (-6048*pow(chi1,3)*(-1 - sqrt(1 - 4*eta) + (3 + sqrt(1 - 4*eta))*eta) +
        chi1*(287213*(1 + sqrt(1 - 4*eta)) - 4*(93414 + 2083*sqrt(1 - 4*eta))*eta - 35632*pow(eta,2)) +
        chi2*(-((287213 + 6048*pow(chi2,2))*(-1 + sqrt(1 - 4*eta))) +
           4*(-93414 + 1512*pow(chi2,2)*(-3 + sqrt(1 - 4*eta)) + 2083*sqrt(1 - 4*eta))*eta - 35632*pow(eta,2)) + 42840*(-1 + 4*eta)*Pi))/96768. -
   (f*pow(Pi,2)*(-336*(-3248849057 + 1809550512*pow(chi1,2) - 2954929824*chi1*chi2 + 1809550512*pow(chi2,2))*pow(eta,2) -
        324322727232*pow(eta,3) + 7*(177520268561 + 29362199328*pow(chi2,2) + 29362199328*pow(chi1,2)*(1 + sqrt(1 - 4*eta)) -
           29362199328*pow(chi2,2)*sqrt(1 - 4*eta) + 12160253952*(chi1 + chi2 + chi1*sqrt(1 - 4*eta) - chi2*sqrt(1 - 4*eta))*Pi) +
        12*eta*(-545384828789 + 49568837472*chi1*chi2 - 12312458928*pow(chi2,2) + 77616*pow(chi1,2)*(-158633 + 282718*sqrt(1 - 4*eta)) -
           21943440288*pow(chi2,2)*sqrt(1 - 4*eta) - 8345272320*(chi1 + chi2)*Pi + 21384760320*pow(Pi,2))))/3.0042980352e10 +
   2.333333333333333*pow(f,1.333333333333333)*rho1 + 2.6666666666666665*pow(f,1.6666666666666665)*rho2 + 3*pow(f,2)*rho3;
}

/////////////////////////// Amplitude: Merger-Ringdown functions ///////////////////////

// Phenom coefficients gamma1, ..., gamma3
// AmpMRDAnsatzFunc[]

static double gamma1_fun(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 0.006929325756444066 + 0.028178602456554704*eta
    + (0.004712463375341761 - 0.10042540687248233*eta + 0.1913418949251075*eta2)*xi
    + (0.000992288691457378 - 0.0735538521653164*eta + 0.1623598588976758*eta2)*xi2
    + (-0.00012994126836641466 - 0.014868851443782867*eta + 0.03420784550184175*eta2)*xi3;
}

static double gamma2_fun(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 1.0338973605758845 - 0.0314363028397228*eta
    + (0.6461456789695995 - 8.452571007081147*eta + 25.039441014539896*eta2)*xi
    + (0.6170831015808893 - 13.59833956349621*eta + 42.618165831477675*eta2)*xi2
    + (0.20878227999538593 - 4.991235292948315*eta + 15.879490996567085*eta2)*xi3;
}

static double gamma3_fun(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 1.3638671714681179 - 0.2648078572888508*eta
    + (0.17346232045631538 - 2.7148000635403484*eta + 7.242203410070453*eta2)*xi
    + (0.2115831718297191 - 3.616857536014833*eta + 10.649274664019698*eta2)*xi2
    + (0.07706402964475789 - 1.275335316106863*eta + 3.874082717359707*eta2)*xi3;
}

static double AmpMRDAnsatz(double f, IMRPhenomDAmplitudeCoefficients* p) {
  double fRD = p->fRD;
  double fDM = p->fDM;
  double gamma1 = p->gamma1;
  double gamma2 = p->gamma2;
  double gamma3 = p->gamma3;

  return exp( -(f - fRD)*gamma2 / (fDM*gamma3) )
    * (fDM*gamma3*fabs(gamma1)) / (pow(f - fRD,2) + pow(fDM,2)*pow(gamma3,2));
}

static double DAmpMRDAnsatz(double f, IMRPhenomDAmplitudeCoefficients* p) {
  double fRD = p->fRD;
  double fDM = p->fDM;
  double gamma1 = p->gamma1;
  double gamma2 = p->gamma2;
  double gamma3 = p->gamma3;

  return (-2*fDM*(f - fRD)*gamma3*fabs(gamma1)) / ( exp(((f - fRD)*gamma2)/(fDM*gamma3)) * pow(pow(f - fRD,2) + pow(fDM,2)*pow(gamma3,2),2)) -
   (gamma2*fabs(gamma1)) / ( exp(((f - fRD)*gamma2)/(fDM*gamma3)) * (pow(f - fRD,2) + pow(fDM,2)*pow(gamma3,2)));
}

static double fmaxCalc(IMRPhenomDAmplitudeCoefficients* p) {
  double fRD = p->fRD;
  double fDM = p->fDM;
  double gamma2 = p->gamma2;
  double gamma3 = p->gamma3;

  // NOTE: There's a problem with this expression becoming imaginary if the spin is large.
  //return fabs(fRD + (fDM*(-1 + sqrt(1 - pow(gamma2,2)))*gamma3)/gamma2);
  //return fabs(fRD + (fDM*(-1 + sqrt( fabs(1 - pow(gamma2,2))))*gamma3)/gamma2);

  // Possible fix: if gamma2 >= 1 then set the square root term to zero.
  if (gamma2 <= 1)
    return fabs(fRD + (fDM*(-1 + sqrt(1 - pow(gamma2,2)))*gamma3)/gamma2);
  else
    return fabs(fRD + (fDM*(-1)*gamma3)/gamma2);
}

///////////////////////////// Amplitude: Intermediate functions ////////////////////////

// Phenom coefficients delta0, ..., delta4 determined from collocation method
// (constraining 3 values and 2 derivatives)
// AmpIntAnsatzFunc[]


static double AmpIntAnsatz(double Mf, IMRPhenomDAmplitudeCoefficients* p) {
  double Mf2 = Mf*Mf;
  double Mf3 = Mf*Mf2;
  double Mf4 = Mf*Mf3;
  return p->delta0 + p->delta1*Mf + p->delta2*Mf2 + p->delta3*Mf3 + p->delta4*Mf4;
}

// AmpIntColFitCoeff\[Chi]PNFunc[eta, chiPN][[2]] // CForm
static double AmpIntColFitCoeff(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 0.8193644340977788 + 2.519365484276569*eta
    + (1.1630617232064524 - 2.302683026216893*eta + 6.086868143445912*eta2)*xi
    + (0.7516886002177731 - 2.5295728215368083*eta + 6.06319783479259*eta2)*xi2
    + (0.1736858389317789 - 0.7143291487179297*eta + 1.7339015693455153*eta2)*xi3;
}

static double delta0_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d) {
  double f1 = p->f1;
  double f2 = p->f2;
  double f3 = p->f3;
  double v1 = p->v1;
  double v2 = p->v2;
  double v3 = p->v3;
  double d1 = p->d1;
  double d2 = p->d2;

  double f12 = d->f12;
  double f13 = d->f13;
  double f14 = d->f14;
  double f15 = d->f15;
  double f22 = d->f22;
  double f23 = d->f23;
  double f24 = d->f24;
  double f32 = d->f32;
  double f33 = d->f33;
  double f34 = d->f34;
  double f35 = d->f35;

  return -((d2*f15*f22*f3 - 2*d2*f14*f23*f3 + d2*f13*f24*f3 - d2*f15*f2*f32 + d2*f14*f22*f32
  - d1*f13*f23*f32 + d2*f13*f23*f32 + d1*f12*f24*f32 - d2*f12*f24*f32 + d2*f14*f2*f33
  + 2*d1*f13*f22*f33 - 2*d2*f13*f22*f33 - d1*f12*f23*f33 + d2*f12*f23*f33 - d1*f1*f24*f33
  - d1*f13*f2*f34 - d1*f12*f22*f34 + 2*d1*f1*f23*f34 + d1*f12*f2*f35 - d1*f1*f22*f35
  + 4*f12*f23*f32*v1 - 3*f1*f24*f32*v1 - 8*f12*f22*f33*v1 + 4*f1*f23*f33*v1 + f24*f33*v1
  + 4*f12*f2*f34*v1 + f1*f22*f34*v1 - 2*f23*f34*v1 - 2*f1*f2*f35*v1 + f22*f35*v1 - f15*f32*v2
  + 3*f14*f33*v2 - 3*f13*f34*v2 + f12*f35*v2 - f15*f22*v3 + 2*f14*f23*v3 - f13*f24*v3
  + 2*f15*f2*f3*v3 - f14*f22*f3*v3 - 4*f13*f23*f3*v3 + 3*f12*f24*f3*v3 - 4*f14*f2*f32*v3
  + 8*f13*f22*f32*v3 - 4*f12*f23*f32*v3) / (pow(f1 - f2,2)*pow(f1 - f3,3)*pow(-f2 + f3,2)));
}

static double delta1_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d) {
  double f1 = p->f1;
  double f2 = p->f2;
  double f3 = p->f3;
  double v1 = p->v1;
  double v2 = p->v2;
  double v3 = p->v3;
  double d1 = p->d1;
  double d2 = p->d2;

  double f12 = d->f12;
  double f13 = d->f13;
  double f14 = d->f14;
  double f15 = d->f15;
  double f22 = d->f22;
  double f23 = d->f23;
  double f24 = d->f24;
  double f32 = d->f32;
  double f33 = d->f33;
  double f34 = d->f34;
  double f35 = d->f35;

  return -((-(d2*f15*f22) + 2*d2*f14*f23 - d2*f13*f24 - d2*f14*f22*f3 + 2*d1*f13*f23*f3
  + 2*d2*f13*f23*f3 - 2*d1*f12*f24*f3 - d2*f12*f24*f3 + d2*f15*f32 - 3*d1*f13*f22*f32
  - d2*f13*f22*f32 + 2*d1*f12*f23*f32 - 2*d2*f12*f23*f32 + d1*f1*f24*f32 + 2*d2*f1*f24*f32
  - d2*f14*f33 + d1*f12*f22*f33 + 3*d2*f12*f22*f33 - 2*d1*f1*f23*f33 - 2*d2*f1*f23*f33
  + d1*f24*f33 + d1*f13*f34 + d1*f1*f22*f34 - 2*d1*f23*f34 - d1*f12*f35 + d1*f22*f35
  - 8*f12*f23*f3*v1 + 6*f1*f24*f3*v1 + 12*f12*f22*f32*v1 - 8*f1*f23*f32*v1 - 4*f12*f34*v1
  + 2*f1*f35*v1 + 2*f15*f3*v2 - 4*f14*f32*v2 + 4*f12*f34*v2 - 2*f1*f35*v2 - 2*f15*f3*v3
  + 8*f12*f23*f3*v3 - 6*f1*f24*f3*v3 + 4*f14*f32*v3 - 12*f12*f22*f32*v3 + 8*f1*f23*f32*v3)
  / (pow(f1 - f2,2)*pow(f1 - f3,3)*pow(-f2 + f3,2)));
}

static double delta2_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d) {
  double f1 = p->f1;
  double f2 = p->f2;
  double f3 = p->f3;
  double v1 = p->v1;
  double v2 = p->v2;
  double v3 = p->v3;
  double d1 = p->d1;
  double d2 = p->d2;

  double f12 = d->f12;
  double f13 = d->f13;
  double f14 = d->f14;
  double f15 = d->f15;
  double f23 = d->f23;
  double f24 = d->f24;
  double f32 = d->f32;
  double f33 = d->f33;
  double f34 = d->f34;
  double f35 = d->f35;

  return -((d2*f15*f2 - d1*f13*f23 - 3*d2*f13*f23 + d1*f12*f24 + 2*d2*f12*f24 - d2*f15*f3
  + d2*f14*f2*f3 - d1*f12*f23*f3 + d2*f12*f23*f3 + d1*f1*f24*f3 - d2*f1*f24*f3 - d2*f14*f32
  + 3*d1*f13*f2*f32 + d2*f13*f2*f32 - d1*f1*f23*f32 + d2*f1*f23*f32 - 2*d1*f24*f32 - d2*f24*f32
  - 2*d1*f13*f33 + 2*d2*f13*f33 - d1*f12*f2*f33 - 3*d2*f12*f2*f33 + 3*d1*f23*f33 + d2*f23*f33
  + d1*f12*f34 - d1*f1*f2*f34 + d1*f1*f35 - d1*f2*f35 + 4*f12*f23*v1 - 3*f1*f24*v1 + 4*f1*f23*f3*v1
  - 3*f24*f3*v1 - 12*f12*f2*f32*v1 + 4*f23*f32*v1 + 8*f12*f33*v1 - f1*f34*v1 - f35*v1 - f15*v2
  - f14*f3*v2 + 8*f13*f32*v2 - 8*f12*f33*v2 + f1*f34*v2 + f35*v2 + f15*v3 - 4*f12*f23*v3 + 3*f1*f24*v3
  + f14*f3*v3 - 4*f1*f23*f3*v3 + 3*f24*f3*v3 - 8*f13*f32*v3 + 12*f12*f2*f32*v3 - 4*f23*f32*v3)
  / (pow(f1 - f2,2)*pow(f1 - f3,3)*pow(-f2 + f3,2)));
}

static double delta3_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d) {
  double f1 = p->f1;
  double f2 = p->f2;
  double f3 = p->f3;
  double v1 = p->v1;
  double v2 = p->v2;
  double v3 = p->v3;
  double d1 = p->d1;
  double d2 = p->d2;

  double f12 = d->f12;
  double f13 = d->f13;
  double f14 = d->f14;
  double f22 = d->f22;
  double f24 = d->f24;
  double f32 = d->f32;
  double f33 = d->f33;
  double f34 = d->f34;

  return -((-2*d2*f14*f2 + d1*f13*f22 + 3*d2*f13*f22 - d1*f1*f24 - d2*f1*f24 + 2*d2*f14*f3
  - 2*d1*f13*f2*f3 - 2*d2*f13*f2*f3 + d1*f12*f22*f3 - d2*f12*f22*f3 + d1*f24*f3 + d2*f24*f3
  + d1*f13*f32 - d2*f13*f32 - 2*d1*f12*f2*f32 + 2*d2*f12*f2*f32 + d1*f1*f22*f32 - d2*f1*f22*f32
  + d1*f12*f33 - d2*f12*f33 + 2*d1*f1*f2*f33 + 2*d2*f1*f2*f33 - 3*d1*f22*f33 - d2*f22*f33
  - 2*d1*f1*f34 + 2*d1*f2*f34 - 4*f12*f22*v1 + 2*f24*v1 + 8*f12*f2*f3*v1 - 4*f1*f22*f3*v1
  - 4*f12*f32*v1 + 8*f1*f2*f32*v1 - 4*f22*f32*v1 - 4*f1*f33*v1 + 2*f34*v1 + 2*f14*v2
  - 4*f13*f3*v2 + 4*f1*f33*v2 - 2*f34*v2 - 2*f14*v3 + 4*f12*f22*v3 - 2*f24*v3 + 4*f13*f3*v3
  - 8*f12*f2*f3*v3 + 4*f1*f22*f3*v3 + 4*f12*f32*v3 - 8*f1*f2*f32*v3 + 4*f22*f32*v3)
  / (pow(f1 - f2,2)*pow(f1 - f3,3)*pow(-f2 + f3,2)));
}

static double delta4_fun(IMRPhenomDAmplitudeCoefficients* p, DeltaUtility* d) {
  double f1 = p->f1;
  double f2 = p->f2;
  double f3 = p->f3;
  double v1 = p->v1;
  double v2 = p->v2;
  double v3 = p->v3;
  double d1 = p->d1;
  double d2 = p->d2;

  double f12 = d->f12;
  double f13 = d->f13;
  double f22 = d->f22;
  double f23 = d->f23;
  double f32 = d->f32;
  double f33 = d->f33;

  return -((d2*f13*f2 - d1*f12*f22 - 2*d2*f12*f22 + d1*f1*f23 + d2*f1*f23 - d2*f13*f3 + 2*d1*f12*f2*f3
  + d2*f12*f2*f3 - d1*f1*f22*f3 + d2*f1*f22*f3 - d1*f23*f3 - d2*f23*f3 - d1*f12*f32 + d2*f12*f32
  - d1*f1*f2*f32 - 2*d2*f1*f2*f32 + 2*d1*f22*f32 + d2*f22*f32 + d1*f1*f33 - d1*f2*f33 + 3*f1*f22*v1
  - 2*f23*v1 - 6*f1*f2*f3*v1 + 3*f22*f3*v1 + 3*f1*f32*v1 - f33*v1 - f13*v2 + 3*f12*f3*v2 - 3*f1*f32*v2
  + f33*v2 + f13*v3 - 3*f1*f22*v3 + 2*f23*v3 - 3*f12*f3*v3 + 6*f1*f2*f3*v3 - 3*f22*f3*v3)
  / (pow(f1 - f2,2)*pow(f1 - f3,3)*pow(-f2 + f3,2)));
}

// Calculates delta_i's
static void ComputeDeltasFromCollocation(IMRPhenomDAmplitudeCoefficients* p) {

  // Four evenly spaced collocation points in the interval [f1,f4].
  double f1 = 0.014;
  double f3 = p->fmaxCalc;
  double dfx = (f3 - f1)/2.0;
  double f2 = f1 + dfx;

  // v1 is inspiral model evaluated at f1
  // d1 is derivative of inspiral model evaluated at f1
  double v1 = AmpInsAnsatz(f1, p); // FIXME: depends on rho_i's
  double d1 = DAmpInsAnsatz(f1, p); // FIXME: depends on rho_i's

  // v4??? is merger-ringdown model evaluated at f3
  // d2 is derivative of merger-ringdown model evaluated at f3
  //
  double v3 = AmpMRDAnsatz(f3, p); // FIXME: depends on gamma_i's
  double d2 = DAmpMRDAnsatz(f3, p); // FIXME: depends on gamma_i's

  // v2 is the value of the amplitude evaluated at f2
  // they come from the fit of the collocation points in the intermediate region
  double v2 = AmpIntColFitCoeff(p->eta, p->chi);

  p->f1 = f1;
  p->f2 = f2;
  p->f3 = f3;
  p->v1 = v1;
  p->v2 = v2;
  p->v3 = v3;
  p->d1 = d1;
  p->d2 = d2;

  // Now compute the delta_i's from the collocation coefficients
  // Precompute common quantities here and pass along to delta functions.
  DeltaUtility d;
  d.f12 = f1*f1;
  d.f13 = f1*d.f12;
  d.f14 = f1*d.f13;
  d.f15 = f1*d.f14;
  d.f22 = f2*f2;
  d.f23 = f2*d.f22;
  d.f24 = f2*d.f23;
  d.f32 = f3*f3;
  d.f33 = f3*d.f32;
  d.f34 = f3*d.f33;
  d.f35 = f3*d.f34;
  p->delta0 = delta0_fun(p, &d);
  p->delta1 = delta1_fun(p, &d);
  p->delta2 = delta2_fun(p, &d);
  p->delta3 = delta3_fun(p, &d);
  p->delta4 = delta4_fun(p, &d);
}


///////////////////////////// Amplitude: glueing function ////////////////////////////

static IMRPhenomDAmplitudeCoefficients* ComputeIMRPhenomDAmplitudeCoefficients(double eta, double chi1, double chi2) {
  IMRPhenomDAmplitudeCoefficients *p = (IMRPhenomDAmplitudeCoefficients *) XLALMalloc(sizeof(IMRPhenomDAmplitudeCoefficients));

  p->eta = eta;
  p->chi1 = chi1;
  p->chi2 = chi2;

  p->q = (1.0 + sqrt(1.0 - 4.0*eta) - 2.0*eta) / (2.0*eta);
  p->chi = chiPN(eta, chi1, chi2);

  p->fRD = fring(eta, chi1, chi2);
  p->fDM = fdamp(eta, chi1, chi2);

  // Compute gamma_i's, rho_i's first then delta_i's
  p->gamma1 = gamma1_fun(eta, p->chi);
  p->gamma2 = gamma2_fun(eta, p->chi);
  p->gamma3 = gamma3_fun(eta, p->chi);

  p->fmaxCalc = fmaxCalc(p); // NOTE: needs the MR gamma_i coeffs already calculated in p!

  p->rho1 = rho1_fun(eta, p->chi);
  p->rho2 = rho2_fun(eta, p->chi);
  p->rho3 = rho3_fun(eta, p->chi);

  // compute delta_i's
  ComputeDeltasFromCollocation(p);

  return p;
}

// Call ComputeIMRPhenomDAmplitudeCoefficients() first!
static double IMRPhenDAmplitude(double f, IMRPhenomDAmplitudeCoefficients *p) {
  // The inspiral, intermediate and merger-ringdown amplitude parts
  // FIXME: could avoid evaluating two of those when we are far enough away from one of the transition frequencies

  // Transition frequencies
  p->fInsJoin = 0.014;
  p->fMRDJoin = p->fmaxCalc;

  double AmpPreFac = amp0Func(p->eta) * pow(f, -7.0/6.0);
  double AmpIns = AmpPreFac * AmpInsAnsatz(f, p);
  double AmpInt = AmpPreFac * AmpIntAnsatz(f, p);
  double AmpMRD = AmpPreFac * AmpMRDAnsatz(f, p);

  double df = 0.001; // frequency window length for Planck taper
  double PT1 = PlanckTaper(f, p->fInsJoin, p->fInsJoin + df);
  double PT2 = PlanckTaper(f, p->fMRDJoin, p->fMRDJoin + df);

  return (1.0 - PT1) * AmpIns + PT1 * AmpInt * (1.0 - PT2) + PT2 * AmpMRD;
}

/********************************* Phase functions *********************************/

////////////////////////////// Phase: Ringdown functions ///////////////////////////

// alpha_i i=1,2,3 are the phenomenological intermediate coefficients depending on eta and chiPN
// PhiRingdownAnsatz is the ringdown phasing in terms of the alpha_i coefficients

static double alpha1Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 48.268391497848555 + 611.7245082528094*eta
    + (-22.52115933057628 + 2432.0201299154805*eta - 6046.652884557147*eta2)*xi
    + (-54.362482633751085 + 3011.8243021641515*eta - 9368.252485211071*eta2)*xi2
    + (-19.68404253845273 + 1002.5518994042775*eta - 3362.9647886372677*eta2)*xi3;
}

static double alpha2Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return -0.0673058164684716 - 0.1830627533897372*eta
    + (-0.18049467080235737 + 1.1177063913620364*eta - 2.9518177020311915*eta2)*xi
    + (-0.16615004046059806 + 1.7202871344838704*eta - 4.7026420500401205*eta2)*xi2
    + (-0.04865897674674206 + 0.6093863685081735*eta - 1.7372061063559312*eta2)*xi3;
}

static double alpha3Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 6.6200489338307085 - 380.80778662617956*eta
    + (7.990338219795906 - 1554.1611800117628*eta + 3689.0921582110873*eta2)*xi
    + (19.552335169497248 - 1775.9434000999295*eta + 5262.501198452011*eta2)*xi2
    + (8.970735925928109 - 574.9441185334387*eta + 1841.0002392874653*eta2)*xi3;
}

static double alpha4Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return -0.03717271633074854 + 1.4261077580289205*eta
    + (-0.10487146504736497 + 1.0535635483450132*eta - 0.20431576108236085*eta2)*xi
    + (-0.09032850971624289 + 0.8708538997789368*eta + 0.02651866135290333*eta2)*xi2
    + (-0.026912018722678448 + 0.26839255388193417*eta - 0.05860768899388999*eta2)*xi3;
}

static double alpha5Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 0.9929836969310366 + 0.003010359365293129*eta
  + (-0.11319782076414776 + 1.9896055092224239*eta - 6.024331489827541*eta2)*xi
  + (-0.13174717365370286 + 2.6346685770150726*eta - 8.278336740553794*eta2)*xi2
  + (-0.044075227973024274 + 0.9030871993168443*eta - 2.860732570370078*eta2)*xi3;
}

static double PhiMRDAnsatzInt(double f, IMRPhenomDPhaseCoefficients *p) {
  return -(p->alpha2/f) + (4*p->alpha3*pow(f,0.75))/3. + p->alpha1*f + p->alpha4*atan((f - p->alpha5*p->fRD)/p->fDM);
}

static double DPhiMRD(double f, IMRPhenomDPhaseCoefficients *p) {
  return (p->alpha1 + p->alpha2/pow(f,2) + p->alpha3/pow(f,0.25) + p->alpha4/(p->fDM*(1 + pow(f - p->alpha5*p->fRD,2)/pow(p->fDM,2)))) / p->eta;
}

///////////////////////////// Phase: Intermediate functions /////////////////////////////

// beta_i i=1,2,3 are the phenomenological intermediate coefficients depending on eta and chiPN
// PhiIntAnsatz is the intermediate phasing in terms of the beta_i coefficients


// \[Beta]1Fit = PhiIntFitCoeff\[Chi]PNFunc[\[Eta], \[Chi]PN][[1]]

static double beta1Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 97.45117083987036 - 41.31509513566375*eta
    + (148.39440087200862 - 1358.654044408998*eta + 2610.1336681598723*eta2)*xi
    + (131.6877588301213 - 1346.3392425741838*eta + 2636.2940288867408*eta2)*xi2
    + (38.57251427350905 - 392.6296853455257*eta + 770.252899105536*eta2)*xi3;
}

static double beta2Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return -3.2804577873466485 - 9.06185751848799*eta
    + (-12.414871795376841 + 55.49086816210587*eta - 106.12827149099643*eta2)*xi
    + (-11.957149752470606 + 76.85798857225546*eta - 155.4568382340443*eta2)*xi2
    + (-3.414799615653319 + 25.58891148247438*eta - 54.44356409784692*eta2)*xi3;
}

static double beta3Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return -0.000025165874934707377 + 0.00001979585325937374*eta
    + (-0.000018380049120626495 + 0.000021869365408771506*eta + 0.00008271002269128104*eta2)*xi
    + (7.162864746777629e-6 - 0.000055897084288544904*eta + 0.00019176840794129924*eta2)*xi2
    + (5.4509322526231396e-6 - 0.00003224698772106724*eta + 0.00007984258587881925*eta2)*xi3;
}

static double PhiIntAnsatz(double Mf, IMRPhenomDPhaseCoefficients *p) {
  return  p->beta1*Mf - p->beta3/(3.*pow(Mf,3)) + p->beta2*log(Mf);
}

static double DPhiIntAnsatz(double Mf, IMRPhenomDPhaseCoefficients *p) {
  return (p->beta1 + p->beta3/pow(Mf,4) + p->beta2/Mf) / p->eta;
}

static double DPhiIntTemp(double ff, IMRPhenomDPhaseCoefficients *p) {
  double eta = p->eta;
  double beta1 = p->beta1;
  double beta2 = p->beta2;
  double beta3 = p->beta3;
  double C2Int = p->C2Int;

  return C2Int + (beta1 + beta3/pow(ff,4) + beta2/ff)/eta;
}

///////////////////////////// Phase: Inspiral functions /////////////////////////////

// sigma_i i=1,2,3,4 are the phenomenological inspiral coefficients depending on eta and chiPN
// PhiInsAnsatzInt is a souped up TF2 phasing which depends on the sigma_i coefficients

static double sigma1Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 2096.155693064352 + 1464.8669504930606*eta
    + (1307.4284854073637 + 18366.582650093205*eta - 43679.41361741789*eta2)*xi
    + (-840.5072830340943 + 32136.410168522973*eta - 108834.93282891942*eta2)*xi2
    + (449.7311701253021 + 8385.049840384685*eta - 44612.56586161423*eta2)*xi3;
}

static double sigma2Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return -10114.056476682446 - 44631.01107848535*eta
    + (-6541.308770527784 - 266959.2342301296*eta + 686328.3232181441*eta2)*xi
    + (3405.6372141416923 - 437507.720946488*eta + 1.6318171313072182e6*eta2)*xi2
    + (-7462.6485633155735 - 114585.25183484034*eta + 674402.4691673126*eta2)*xi3;
}

static double sigma3Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 22933.658292557113 + 230960.008073928*eta
    + (14961.084015996932 + 1.1940181344003242e6*eta - 3.1042239706327254e6*eta2)*xi
    + (-3038.1665950966894 + 1.8720322854868704e6*eta - 7.309145014733551e6*eta2)*xi2
    + (42738.22871637542 + 467502.01890796213*eta - 3.0648534997008806e6*eta2)*xi3;
}

static double sigma4Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return -14621.715252090613 - 377812.8578198768*eta
    + (-9608.682697135695 - 1.7108925259726832e6*eta + 4.332924603450101e6*eta2)*xi
    + (-22366.683297514493 - 2.5019716395126972e6*eta + 1.0274495906301545e7*eta2)*xi2
    + (-85360.30079323666 - 570025.3446157051*eta + 4.396844348659892e6*eta2)*xi3;
}

static double PhiInsAnsatzInt(double Mf, IMRPhenomDPhaseCoefficients *p) {
  double eta = p->eta;
  double chi1 = p->chi1;
  double chi2 = p->chi2;
  double sigma1 = p->sigma1;
  double sigma2 = p->sigma2;
  double sigma3 = p->sigma3;
  double sigma4 = p->sigma4;

  // Obtain LAL TaylorF2 phasing coefficients which are tested in ../test/PNCoefficients.c.
  // FIXME: MP: we probably only want to call XLALSimInspiralTaylorF2AlignedPhasing() once per waveform
  // and save the coefficients somewhere.
  // E.g., we could just store the pointer to the PNPhasingSeries struct in IMRPhenomDPhaseCoefficients!
  // PNPhasingSeries *pn = NULL;
  // // Convention m1 >= m2
  // double m1 = 0.5 * (1.0 + sqrt(1.0 - 4.0*eta));
  // double m2 = 0.5 * (1.0 - sqrt(1.0 - 4.0*eta));
  // XLALSimInspiralTaylorF2AlignedPhasing(&pn, m1, m2, chi1, chi2, 1, 1, LAL_SIM_INSPIRAL_SPIN_ORDER_35PN);

  // Assemble PN phasing series
  const double v = cbrt(Pi*Mf);
  const double logv = log(v);
  const double v2 = v * v;
  const double v3 = v * v2;
  const double v4 = v * v3;
  const double v5 = v * v4;
  const double v6 = v * v5;
  // const double v7 = v * v6;

  // double phasing = 0.0;
  // phasing += pn->v[7] * v7;
  // phasing += (pn->v[6] + pn->vlogv[6] * logv) * v6;
  // phasing += (pn->v[5] + pn->vlogv[5] * logv) * v5;
  // phasing += pn->v[4] * v4;
  // phasing += pn->v[3] * v3;
  // phasing += pn->v[2] * v2;
  // phasing += pn->v[0]; // * v^0
  // phasing /= v5;
  // phasing -= LAL_PI_4;
  // LALFree(pn);

  // Now add higher order terms that were calibrated for PhenomD
  // phasing += (
  //           (sigma1/Pi) * v3
  //         + (3.0/(4.0*pow(Pi,4.0/3.0))*sigma2) * v4
  //         + (3.0/(5.0*pow(Pi,5.0/3.0))*sigma3) * v5
  //         + (1.0/(2.0*Pi*Pi)*sigma4) * v6
  //         ) / eta;

  // Note: MP: this is the CForm of the TF2 phase from Mathematica
  // Ultimately, we want to use the expression above using LAL TF2 coefficients.
  // FIXME: optimize
  double eta2 = eta*eta;
  double Pi2 = Pi*Pi;
  double phasing = (6065*(chi1 + chi2))/1728. - (35*(chi1 - chi2)*sqrt(1 - 4*eta))/192. \
- (732985*(chi1 + chi2))/(193536.*eta) - (732985*(chi1 - chi2)*sqrt(1 \
- 4*eta))/(193536.*eta) + (85*(chi1 + chi2)*eta)/192. - (161*Pi)/384. \
+ (38645*Pi)/(32256.*eta) + 3/(128.*eta*v5) + 55/(384.*v3) + \
3715/(32256.*eta*v3) - (19*(chi1 + chi2))/(64.*v2) + (113*(chi1 + \
chi2))/(256.*eta*v2) + (113*(chi1 - chi2)*sqrt(1 - \
4*eta))/(256.*eta*v2) - (3*Pi)/(8.*eta*v2) + 27145/(21504.*v) + \
(75*pow(chi1 - chi2,2))/(64.*v) + (15*pow(chi1 + \
chi2,2))/(1024.*v) + 15293365/(2.1676032e7*eta*v) - (1215*pow(chi1 \
- chi2,2))/(4096.*eta*v) - (1215*pow(chi1 + chi2,2))/(4096.*eta*v) \
- (1215*(chi1 - chi2)*(chi1 + chi2)*sqrt(1 - 4*eta))/(2048.*eta*v) + \
(3085*eta)/(3072.*v) - (15737765635*v)/1.30056192e8 + \
(11583231236531*v)/(2.0028653568e11*eta) + (76055*eta*v)/73728. - \
(127825*eta2*v)/55296. - (107*EulerGamma*v)/(14.*eta) - (195*(chi1 + \
chi2)*Pi*v)/32. + (1135*(chi1 + chi2)*Pi*v)/(128.*eta) + (1135*(chi1 \
- chi2)*sqrt(1 - 4*eta)*Pi*v)/(128.*eta) + (2255*Pi2*v)/512. - \
(5*Pi2*v)/eta + (10566655595*(chi1 + chi2)*v2)/6.5028096e7 + \
(26804935*(chi1 - chi2)*sqrt(1 - 4*eta)*v2)/516096. - \
(25150083775*(chi1 + chi2)*v2)/(2.60112384e8*eta) - \
(25150083775*(chi1 - chi2)*sqrt(1 - 4*eta)*v2)/(2.60112384e8*eta) - \
(1042165*(chi1 + chi2)*eta*v2)/258048. - (1985*(chi1 - chi2)*sqrt(1 - \
4*eta)*eta*v2)/4096. + (5345*(chi1 + chi2)*eta2*v2)/3072. + \
(378515*Pi*v2)/64512. + (77096675*Pi*v2)/(1.0838016e7*eta) - \
(74045*eta*Pi*v2)/32256. + ((sigma1*v3)/Pi + \
(3*sigma2*v4)/(4.*pow(Pi,1.3333333333333333)) + \
(3*sigma3*v5)/(5.*pow(Pi,1.6666666666666667)) + \
(sigma4*v6)/(2.*Pi2))/eta - (107*v*log(2))/(7.*eta) + (6065*(chi1 + \
chi2)*logv)/576. - (35*(chi1 - chi2)*sqrt(1 - 4*eta)*logv)/64. - \
(732985*(chi1 + chi2)*logv)/(64512.*eta) - (732985*(chi1 - \
chi2)*sqrt(1 - 4*eta)*logv)/(64512.*eta) + (85*(chi1 + \
chi2)*eta*logv)/64. - (65*Pi*logv)/128. + \
(38645*Pi*logv)/(10752.*eta) - (107*v*logv)/(14.*eta);

  return phasing;
}

static double DPhiInsAnsatzInt(double Mf, IMRPhenomDPhaseCoefficients *p) {
  double eta = p->eta;
  double chi1 = p->chi1;
  double chi2 = p->chi2;
  double sigma1 = p->sigma1;
  double sigma2 = p->sigma2;
  double sigma3 = p->sigma3;
  double sigma4 = p->sigma4;

  // Obtain LAL TaylorF2 phasing coefficients which are tested in ../test/PNCoefficients.c.
  // FIXME: MP: we probably only want to call XLALSimInspiralTaylorF2AlignedPhasing() once per waveform
  // and save the coefficients somewhere.
  // E.g., we could just store the pointer to the PNPhasingSeries struct in IMRPhenomDPhaseCoefficients!
  // PNPhasingSeries *pn = NULL;
  // // Convention m1 >= m2
  // double m1 = 0.5 * (1.0 + sqrt(1.0 - 4.0*eta));
  // double m2 = 0.5 * (1.0 - sqrt(1.0 - 4.0*eta));
  // XLALSimInspiralTaylorF2AlignedPhasing(&pn, m1, m2, chi1, chi2, 1, 1, LAL_SIM_INSPIRAL_SPIN_ORDER_35PN);

  // Assemble PN phasing series
  const double v = cbrt(Pi*Mf);
  const double logv = log(v);
  const double v2 = v * v;
  const double v3 = v * v2;
  const double v4 = v * v3;
  const double v5 = v * v4;
  const double v6 = v * v5;
  const double v7 = v * v6;
  const double v8 = v * v7;

  // Apply the correct prefactors to LAL phase coefficients to get the
  // phase derivative dphi / dMf = dphi/dv * dv/dMf
  // double Dphasing = 0.0;
  // Dphasing += +2.0 * pn->v[7] * v7;
  // Dphasing += (pn->v[6] + pn->vlogv[6] * (1.0 + logv)) * v6;
  // Dphasing += pn->vlogv[5] * v5;
  // Dphasing += -1.0 * pn->v[4] * v4;
  // Dphasing += -2.0 * pn->v[3] * v3;
  // Dphasing += -3.0 * pn->v[2] * v2;
  // Dphasing += -5.0 * pn->v[0];
  // Dphasing /= v8 * 3.0/Pi;
  // LALFree(pn);

  // Now add higher order terms that were calibrated for PhenomD
  // Dphasing += (
  //         sigma1
  //       + (pow(Pi, -1.0/3.0)*sigma2) * v
  //       + (pow(Pi, -2.0/3.0)*sigma3) * v2
  //       + (sigma4/Pi) * v3
  //       ) / eta;

  // Note: MP: this is the CForm of the TF2 phase derivative from Mathematica
  // Ultimately, we want to use the expression above using LAL TF2 coefficients.
  // optimized CForm
  double eta2 = eta*eta;
  double Pi2 = Pi*Pi;
  double Pi3 = Pi*Pi2;
  double Dphasing = (-5*Pi)/(128.*eta*v8) - (55*Pi)/(384.*v6) - \
(3715*Pi)/(32256.*eta*v6) + (19*(chi1 + chi2)*Pi)/(96.*v5) - \
(113*(chi1 + chi2)*Pi)/(384.*eta*v5) - (113*(chi1 - chi2)*sqrt(1 - \
4*eta)*Pi)/(384.*eta*v5) + Pi2/(4.*eta*v5) - (27145*Pi)/(64512.*v4) - \
(25*pow(chi1 - chi2,2)*Pi)/(64.*v4) - (5*pow(chi1 + \
chi2,2)*Pi)/(1024.*v4) - (15293365*Pi)/(6.5028096e7*eta*v4) + \
(405*pow(chi1 - chi2,2)*Pi)/(4096.*eta*v4) + (405*pow(chi1 + \
chi2,2)*Pi)/(4096.*eta*v4) + (405*(chi1 - chi2)*(chi1 + chi2)*sqrt(1 \
- 4*eta)*Pi)/(2048.*eta*v4) - (3085*eta*Pi)/(9216.*v4) + (6065*(chi1 \
+ chi2)*Pi)/(1728.*v3) - (35*(chi1 - chi2)*sqrt(1 - \
4*eta)*Pi)/(192.*v3) - (732985*(chi1 + chi2)*Pi)/(193536.*eta*v3) - \
(732985*(chi1 - chi2)*sqrt(1 - 4*eta)*Pi)/(193536.*eta*v3) + \
(85*(chi1 + chi2)*eta*Pi)/(192.*v3) - (65*Pi2)/(384.*v3) + \
(38645*Pi2)/(32256.*eta*v3) - (15737765635*Pi)/(3.90168576e8*v2) + \
(10052469856691*Pi)/(6.0085960704e11*eta*v2) + \
(76055*eta*Pi)/(221184.*v2) - (127825*eta2*Pi)/(165888.*v2) - \
(107*EulerGamma*Pi)/(42.*eta*v2) - (65*(chi1 + chi2)*Pi2)/(32.*v2) + \
(1135*(chi1 + chi2)*Pi2)/(384.*eta*v2) + (1135*(chi1 - chi2)*sqrt(1 - \
4*eta)*Pi2)/(384.*eta*v2) + (2255*Pi3)/(1536.*v2) - \
(5*Pi3)/(3.*eta*v2) + (10566655595*(chi1 + chi2)*Pi)/(9.7542144e7*v) \
+ (26804935*(chi1 - chi2)*sqrt(1 - 4*eta)*Pi)/(774144.*v) - \
(25150083775*(chi1 + chi2)*Pi)/(3.90168576e8*eta*v) - \
(25150083775*(chi1 - chi2)*sqrt(1 - 4*eta)*Pi)/(3.90168576e8*eta*v) - \
(1042165*(chi1 + chi2)*eta*Pi)/(387072.*v) - (1985*(chi1 - \
chi2)*sqrt(1 - 4*eta)*eta*Pi)/(6144.*v) + (5345*(chi1 + \
chi2)*eta2*Pi)/(4608.*v) + (378515*Pi2)/(96768.*v) + \
(77096675*Pi2)/(1.6257024e7*eta*v) - (74045*eta*Pi2)/(48384.*v) + \
(sigma1 + (sigma2*v)/pow(Pi,0.3333333333333333) + \
(sigma3*v2)/pow(Pi,0.6666666666666666) + (sigma4*v3)/Pi)/eta - \
(107*Pi*log(2))/(21.*eta*v2) - (107*Pi*logv)/(42.*eta*v2);

  return Dphasing;
}

///////////////////////////////// Phase: glueing function ////////////////////////////////

static IMRPhenomDPhaseCoefficients* ComputeIMRPhenomDPhaseCoefficients(double eta, double chi1, double chi2) {
  IMRPhenomDPhaseCoefficients *p = (IMRPhenomDPhaseCoefficients *) XLALMalloc(sizeof(IMRPhenomDPhaseCoefficients));

  // Convention m1 >= m2
  p->eta = eta;
  p->chi1 = chi1;
  p->chi2 = chi2;

  p->q = (1.0 + sqrt(1.0 - 4.0*eta) - 2.0*eta) / (2.0*eta);
  p->chi = chiPN(eta, chi1, chi2);

  p->sigma1 = sigma1Fit(eta, p->chi);
  p->sigma2 = sigma2Fit(eta, p->chi);
  p->sigma3 = sigma3Fit(eta, p->chi);
  p->sigma4 = sigma4Fit(eta, p->chi);

  p->beta1 = beta1Fit(eta, p->chi);
  p->beta2 = beta2Fit(eta, p->chi);
  p->beta3 = beta3Fit(eta, p->chi);

  p->alpha1 = alpha1Fit(eta, p->chi);
  p->alpha2 = alpha2Fit(eta, p->chi);
  p->alpha3 = alpha3Fit(eta, p->chi);
  p->alpha4 = alpha4Fit(eta, p->chi);
  p->alpha5 = alpha5Fit(eta, p->chi);

  p->fRD = fring(eta, chi1, chi2);
  p->fDM = fdamp(eta, chi1, chi2);

  return p;
}

static void ComputeIMRPhenDPhaseConnectionCoefficients(IMRPhenomDPhaseCoefficients *p) {
  double eta = p->eta;

  // Transition frequencies
  p->fInsJoin=0.018;
  p->fMRDJoin=0.5*p->fRD;

  // Compute C1Int and C2Int coeffs
  // Equations to solve for to get C(1) continuous join
  // PhiIns (f)  =   PhiInt (f) + C1Int + C2Int f
  // Joining at fInsJoin
  // PhiIns (fInsJoin)  =   PhiInt (fInsJoin) + C1Int + C2Int fInsJoin
  // PhiIns'(fInsJoin)  =   PhiInt'(fInsJoin) + C2Int
  double DPhiIns = DPhiInsAnsatzInt(p->fInsJoin, p);
  double DPhiInt = DPhiIntAnsatz(p->fInsJoin, p);
  p->C2Int = DPhiIns - DPhiInt;
  p->C1Int = PhiInsAnsatzInt(p->fInsJoin, p)
    - 1.0/eta * PhiIntAnsatz(p->fInsJoin, p) - p->C2Int * p->fInsJoin;

  // Compute C1MRD and C2MRD coeffs
  // Equations to solve for to get C(1) continuous join
  // PhiInsInt (f)  =   PhiMRD (f) + C1MRD + C2MRD f
  // Joining at fMRDJoin
  // Where \[Phi]InsInt(f) is the \[Phi]Ins+\[Phi]Int joined function
  // PhiInsInt (fMRDJoin)  =   PhiMRD (fMRDJoin) + C1MRD + C2MRD fMRDJoin
  // PhiInsInt'(fMRDJoin)  =   PhiMRD'(fMRDJoin) + C2MRD
  // temporary Intermediate Phase function to Join up the Merger-Ringdown
  double PhiIntTempVal = 1.0/eta * PhiIntAnsatz(p->fMRDJoin, p) + p->C1Int + p->C2Int*p->fMRDJoin;
  double DPhiIntTempVal = DPhiIntTemp(p->fMRDJoin, p);
  double DPhiMRDVal = DPhiMRD(p->fMRDJoin, p);
  p->C2MRD = DPhiIntTempVal - DPhiMRDVal;
  p->C1MRD = PhiIntTempVal - 1.0/eta * PhiMRDAnsatzInt(p->fMRDJoin, p) - p->C2MRD*p->fMRDJoin;
}

static double IMRPhenDPhase(double f, IMRPhenomDPhaseCoefficients *p) {
  // The inspiral, intermendiate and merger-ringdown phase parts
  // FIXME: could avoid evaluating two of those when we are far enough away from one of the transition frequencies
  double PhiIns = PhiInsAnsatzInt(f, p);
  double PhiInt = 1.0/p->eta * PhiIntAnsatz(f, p) + p->C1Int + p->C2Int * f;
  double PhiMRD = 1.0/p->eta * PhiMRDAnsatzInt(f, p) + p->C1MRD + p->C2MRD * f;

  double df1 = 0.0001; // frequency window length for Planck taper
  double df2 = 0.005; // frequency window length for Planck taper
  double PT1 = PlanckTaper(f, p->fInsJoin, p->fInsJoin + df1);
  double PT2 = PlanckTaper(f, p->fMRDJoin, p->fMRDJoin + df2);

  return (1.0 - PT1) * PhiIns + PT1 * PhiInt * (1.0 - PT2) + PT2 * PhiMRD;
}
