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
static double StepFunc(const double t, const double t1);

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
static double AmpInsAnsatz(double Mf, IMRPhenomDAmplitudeCoefficients* p);
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
static double PhiInsAnsatzInt(double f, IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn);
static double DPhiInsAnsatzInt(double ff, IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn);

////////////////////////////// Phase: glueing function //////////////////////////////

static IMRPhenomDPhaseCoefficients* ComputeIMRPhenomDPhaseCoefficients(double eta, double chi1, double chi2, double finspin, const LALSimInspiralTestGRParam *extraParams);
static void ComputeIMRPhenDPhaseConnectionCoefficients(IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn);
static double IMRPhenDPhase(double f, IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn);

/*
 *
 * Internal function definitions
 *
 */

////////////////////////////// Miscellaneous functions //////////////////////////////

/**
 * PN reduced spin parameter
 * See Eq 5.9 in http://arxiv.org/pdf/1107.1267v2.pdf
 */
static double chiPN(double eta, double chi1, double chi2) {
  // Convention m1 >= m2 and chi1 is the spin on m1
  double delta = sqrt(1.0 - 4.0*eta);
  double chi_s = (chi1 + chi2) / 2.0;
  double chi_a = (chi1 - chi2) / 2.0;
  return chi_s * (1.0 - eta*76.0/113.0) + delta*chi_a;
}

/**
 * Return the closest higher power of 2
 */
static size_t NextPow2(const size_t n) {
  return 1 << (size_t) ceil(log2(n));
}

/**
 * Step function
 */
static double StepFunc(const double t, const double t1) {
  if (t < t1)
    return 0.0;
  else
    return 1.0;
}

//////////////////////// Final spin, final mass, fring, fdamp ////////////////////////

// Final Spin and Radiated Energy formulas described in 1508.07250

/**
 * Formula to predict the final spin. Equation 3.6 arXiv:1508.07250
 * s defined around Equation 3.6.
 */
static double FinalSpin0815_s(double eta, double s) {
  double eta2 = eta*eta;
  double eta3 = eta2*eta;
  double eta4 = eta3*eta;
  double s2 = s*s;
  double s3 = s2*s;
  double s4 = s3*s;

return 3.4641016151377544*eta - 4.399247300629289*eta2 + 
   9.397292189321194*eta3 - 13.180949901606242*eta4 + 
   (1 - 0.0850917821418767*eta - 5.837029316602263*eta2)*s + 
   (0.1014665242971878*eta - 2.0967746996832157*eta2)*s2 + 
   (-1.3546806617824356*eta + 4.108962025369336*eta2)*s3 + 
   (-0.8676969352555539*eta + 2.064046835273906*eta2)*s4;
}

/**
 * Wrapper function for FinalSpin0815_s.
 */
static double FinalSpin0815(double eta, double chi1, double chi2) {
  // Convention m1 >= m2
  double m1 = 0.5 * (1.0 + sqrt(1.0 - 4.0*eta));
  double m2 = 0.5 * (1.0 - sqrt(1.0 - 4.0*eta));
  double m1s = m1*m1;
  double m2s = m2*m2;
  // s defined around Equation 3.6 arXiv:1508.07250
  double s = (m1s * chi1 + m2s * chi2);
  return FinalSpin0815_s(eta, s);
}

/**
 * Formula to predict the total radiated energy. Equation 3.7 and 3.8 arXiv:1508.07250.
 * Input parameter s defined around Equation 3.7 and 3.8.
 */
static double EradRational0815_s(double eta, double s) {
  double eta2 = eta*eta;
  double eta3 = eta2*eta;
  double eta4 = eta3*eta;

  return ((0.055974469826360077*eta + 0.5809510763115132*eta2 - 0.9606726679372312*eta3 + 3.352411249771192*eta4)*
    (1. + (-0.0030302335878845507 - 2.0066110851351073*eta + 7.7050567802399215*eta2)*s))/(1. + (-0.6714403054720589 - 1.4756929437702908*eta + 7.304676214885011*eta2)*s);
}

/**
 * Wrapper function for EradRational0815_s.
 */
static double EradRational0815(double eta, double chi1, double chi2) {
  // Convention m1 >= m2
  double m1 = 0.5 * (1.0 + sqrt(1.0 - 4.0*eta));
  double m2 = 0.5 * (1.0 - sqrt(1.0 - 4.0*eta));
  double m1s = m1*m1;
  double m2s = m2*m2;
  // arXiv:1508.07250
  double s = (m1s * chi1 + m2s * chi2) / (m1s + m2s);

  return EradRational0815_s(eta, s);
}

/**
 * fring is the real part of the ringdown frequency
 * 1508.07250 figure 9
 */
static double fring(double eta, double chi1, double chi2, double finspin) {
  double return_val;

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *iFring = gsl_spline_alloc(gsl_interp_cspline, QNMData_length);
  gsl_spline_init(iFring, QNMData_a, QNMData_fring, QNMData_length);

  return_val = gsl_spline_eval(iFring, finspin, acc) / (1.0 - EradRational0815(eta, chi1, chi2));

  gsl_spline_free(iFring);
  gsl_interp_accel_free(acc);
  return return_val;
}

/**
 * fdamp is the complex part of the ringdown frequency
 * 1508.07250 figure 9
 */
static double fdamp(double eta, double chi1, double chi2, double finspin) {
  double return_val;
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *iFdamp = gsl_spline_alloc(gsl_interp_cspline, QNMData_length);
  gsl_spline_init(iFdamp, QNMData_a, QNMData_fdamp, QNMData_length);

  return_val = gsl_spline_eval(iFdamp, finspin, acc) / (1.0 - EradRational0815(eta, chi1, chi2));

  gsl_spline_free(iFdamp);
  gsl_interp_accel_free(acc);
  return return_val;
}

/******************************* Amplitude functions *******************************/

/**
 * amplitude scaling factor defined by Eq.17 in 1508.07253.
 */
static double amp0Func(double eta) {
  return (sqrt(2.0/3.0)*sqrt(eta))/pow(LAL_PI,1.0/6.0);
}

///////////////////////////// Amplitude: Inspiral functions /////////////////////////

// Phenom coefficients rho1, ..., rho3 from direct fit
// AmpInsDFFitCoeffChiPNFunc[eta, chiPN]

/**
 * rho_1 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double rho1_fun(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 3931.8979897196696 - 17395.758706812805*eta
  + (3132.375545898835 + 343965.86092361377*eta - 1.2162565819981997e6*eta2)*xi
  + (-70698.00600428853 + 1.383907177859705e6*eta - 3.9662761890979446e6*eta2)*xi2
  + (-60017.52423652596 + 803515.1181825735*eta - 2.091710365941658e6*eta2)*xi3;
}

/**
 * rho_2 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double rho2_fun(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return -40105.47653771657 + 112253.0169706701*eta
  + (23561.696065836168 - 3.476180699403351e6*eta + 1.137593670849482e7*eta2)*xi
  + (754313.1127166454 - 1.308476044625268e7*eta + 3.6444584853928134e7*eta2)*xi2
  + (596226.612472288 - 7.4277901143564405e6*eta + 1.8928977514040343e7*eta2)*xi3;
}

/**
 * rho_3 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double rho3_fun(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 83208.35471266537 - 191237.7264145924*eta +
  (-210916.2454782992 + 8.71797508352568e6*eta - 2.6914942420669552e7*eta2)*xi
  + (-1.9889806527362722e6 + 3.0888029960154563e7*eta - 8.390870279256162e7*eta2)*xi2
  + (-1.4535031953446497e6 + 1.7063528990822166e7*eta - 4.2748659731120914e7*eta2)*xi3;
}

// The Newtonian term in LAL is fine and we should use exactly the same (either hardcoded or call).
// We just use the Mathematica expression for convenience.
/**
 * Inspiral amplitude plus rho phenom coefficents. rho coefficients computed
 * in rho1_fun, rho2_fun, rho3_fun functions.
 * Amplitude is a re-expansion. See 1508.07253 and Equation 29, 30 and Appendix B arXiv:1508.07253 for details
 */
static double AmpInsAnsatz(double Mf, IMRPhenomDAmplitudeCoefficients* p) {
  double eta = p->eta;
  double chi1 = p->chi1;
  double chi2 = p->chi2;
  double rho1 = p->rho1;
  double rho2 = p->rho2;
  double rho3 = p->rho3;

  double chi12 = chi1*chi1;
  double chi22 = chi2*chi2;
  double eta2 = eta*eta;
  double eta3 = eta*eta2;
  double Mf2 = Mf*Mf;
  double Mf3 = Mf*Mf2;
  double Pi = LAL_PI;
  double Pi2 = Pi*Pi;
  double Seta = sqrt(1.0 - 4.0*eta);

  return 1 + ((-969 + 1804*eta)*pow(Pi*Mf,2.0/3.0))/672.
  + ((chi1*(81*(1 + Seta) - 44*eta) + chi2*(81 - 81*Seta - 44*eta))*Mf*Pi)/48.
  + ((-27312085 - 10287648*chi22 - 10287648*chi12*(1 + Seta) + 10287648*chi22*Seta
  + 24*(-1975055 + 857304*chi12 - 994896*chi1*chi2 + 857304*chi22)*eta + 35371056*eta2)
  *pow(Pi*Mf,4.0/3.0))/8.128512e6
  + (pow(Pi*Mf,5.0/3.0)*(chi2*(-285197*(-1 + Seta) + 4*(-91902 + 1579*Seta)*eta - 35632*eta2)
  + chi1*(285197*(1 + Seta) - 4*(91902 + 1579*Seta)*eta - 35632*eta2)
  + 42840*(-1 + 4*eta)*Pi))/32256. - (Mf2*Pi2*(-336*(-3248849057 + 2943675504*chi12
  - 3339284256*chi1*chi2 + 2943675504*chi22)*eta2 - 324322727232*eta3
  - 7*(-177520268561 + 107414046432*chi22 + 107414046432*chi12*(1 + Seta)
  - 107414046432*chi22*Seta + 11087290368*(chi1 + chi2 + chi1*Seta
  - chi2*Seta)*Pi) + 12*eta*(-545384828789 - 176491177632*chi1*chi2
  + 202603761360*chi22 + 77616*chi12*(2610335 + 995766*Seta)
  - 77287373856*chi22*Seta + 5841690624*(chi1 + chi2)*Pi + 21384760320*Pi2)))/6.0085960704e10
  + pow(Mf,7.0/3.0)*rho1 + pow(Mf,8.0/3.0)*rho2 + Mf3*rho3;
}

/**
 * Take the AmpInsAnsatz expression and compute the first derivative
 * with respect to frequency to get the expression below.
 */
static double DAmpInsAnsatz(double Mf, IMRPhenomDAmplitudeCoefficients* p) {
  double eta = p->eta;
  double chi1 = p->chi1;
  double chi2 = p->chi2;
  double rho1 = p->rho1;
  double rho2 = p->rho2;
  double rho3 = p->rho3;

  double chi12 = chi1*chi1;
  double chi22 = chi2*chi2;
  double eta2 = eta*eta;
  double eta3 = eta*eta2;
  double Mf2 = Mf*Mf;
  double Pi = LAL_PI;
  double Pi2 = Pi*Pi;
  double Seta = sqrt(1.0 - 4.0*eta);

   return ((-969 + 1804*eta)*pow(Pi,2.0/3.0))/(1008.*pow(Mf,1.0/3.0))
   + ((chi1*(81*(1 + Seta) - 44*eta) + chi2*(81 - 81*Seta - 44*eta))*Pi)/48.
   + ((-27312085 - 10287648*chi22 - 10287648*chi12*(1 + Seta)
   + 10287648*chi22*Seta + 24*(-1975055 + 857304*chi12 - 994896*chi1*chi2 + 857304*chi22)*eta
   + 35371056*eta2)*pow(Mf,1.0/3.0)*pow(Pi,4.0/3.0))/6.096384e6
   + (5*pow(Mf,2.0/3.0)*pow(Pi,5.0/3.0)*(chi2*(-285197*(-1 + Seta)
   + 4*(-91902 + 1579*Seta)*eta - 35632*eta2) + chi1*(285197*(1 + Seta)
   - 4*(91902 + 1579*Seta)*eta - 35632*eta2) + 42840*(-1 + 4*eta)*Pi))/96768.
   - (Mf*Pi2*(-336*(-3248849057 + 2943675504*chi12 - 3339284256*chi1*chi2 + 2943675504*chi22)*eta2 - 324322727232*eta3
   - 7*(-177520268561 + 107414046432*chi22 + 107414046432*chi12*(1 + Seta) - 107414046432*chi22*Seta
   + 11087290368*(chi1 + chi2 + chi1*Seta - chi2*Seta)*Pi)
   + 12*eta*(-545384828789 - 176491177632*chi1*chi2 + 202603761360*chi22 + 77616*chi12*(2610335 + 995766*Seta)
   - 77287373856*chi22*Seta + 5841690624*(chi1 + chi2)*Pi + 21384760320*Pi2)))/3.0042980352e10
   + (7.0/3.0)*pow(Mf,4.0/3.0)*rho1 + (8.0/3.0)*pow(Mf,5.0/3.0)*rho2 + 3*Mf2*rho3;
}

/////////////////////////// Amplitude: Merger-Ringdown functions ///////////////////////

// Phenom coefficients gamma1, ..., gamma3
// AmpMRDAnsatzFunc[]

/**
 * gamma 1 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double gamma1_fun(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 0.006927402739328343 + 0.03020474290328911*eta
  + (0.006308024337706171 - 0.12074130661131138*eta + 0.26271598905781324*eta2)*xi
  + (0.0034151773647198794 - 0.10779338611188374*eta + 0.27098966966891747*eta2)*xi2
  + (0.0007374185938559283 - 0.02749621038376281*eta + 0.0733150789135702*eta2)*xi3;
}

/**
 * gamma 2 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double gamma2_fun(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 1.010344404799477 + 0.0008993122007234548*eta
  + (0.283949116804459 - 4.049752962958005*eta + 13.207828172665366*eta2)*xi
  + (0.10396278486805426 - 7.025059158961947*eta + 24.784892370130475*eta2)*xi2
  + (0.03093202475605892 - 2.6924023896851663*eta + 9.609374464684983*eta2)*xi3;
}

/**
 * gamma 3 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double gamma3_fun(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 1.3081615607036106 - 0.005537729694807678*eta
  + (-0.06782917938621007 - 0.6689834970767117*eta + 3.403147966134083*eta2)*xi
  + (-0.05296577374411866 - 0.9923793203111362*eta + 4.820681208409587*eta2)*xi2 
  + (-0.006134139870393713 - 0.38429253308696365*eta + 1.7561754421985984*eta2)*xi3;
}

/**
 * Ansatz for the merger-ringdown amplitude. Equation 19 arXiv:1508.07253
 */
static double AmpMRDAnsatz(double f, IMRPhenomDAmplitudeCoefficients* p) {
  double fRD = p->fRD;
  double fDM = p->fDM;
  double gamma1 = p->gamma1;
  double gamma2 = p->gamma2;
  double gamma3 = p->gamma3;
  return exp( -(f - fRD)*gamma2 / (fDM*gamma3) )
    * (fDM*gamma3*gamma1) / (pow(f - fRD,2) + pow(fDM,2)*pow(gamma3,2));
}

/**
 * first frequency derivative of AmpMRDAnsatz
 */
static double DAmpMRDAnsatz(double f, IMRPhenomDAmplitudeCoefficients* p) {
  double fRD = p->fRD;
  double fDM = p->fDM;
  double gamma1 = p->gamma1;
  double gamma2 = p->gamma2;
  double gamma3 = p->gamma3;
  return (-2*fDM*(f - fRD)*gamma3*gamma1) / ( exp(((f - fRD)*gamma2)/(fDM*gamma3)) * pow(pow(f - fRD,2) + pow(fDM,2)*pow(gamma3,2),2)) -
   (gamma2*gamma1) / ( exp(((f - fRD)*gamma2)/(fDM*gamma3)) * (pow(f - fRD,2) + pow(fDM,2)*pow(gamma3,2)));
}

/**
 * Equation 20 arXiv:1508.07253 (called f_peak in paper)
 * analytic location of maximum of AmpMRDAnsatz
 */
static double fmaxCalc(IMRPhenomDAmplitudeCoefficients* p) {
  double fRD = p->fRD;
  double fDM = p->fDM;
  double gamma2 = p->gamma2;
  double gamma3 = p->gamma3;

  // NOTE: There's a problem with this expression from the paper becoming imaginary if gamma2>=1
  // Fix: if gamma2 >= 1 then set the square root term to zero.
  if (gamma2 <= 1)
    return fabs(fRD + (fDM*(-1 + sqrt(1 - pow(gamma2,2)))*gamma3)/gamma2);
  else
    return fabs(fRD + (fDM*(-1)*gamma3)/gamma2);
}

///////////////////////////// Amplitude: Intermediate functions ////////////////////////

// Phenom coefficients delta0, ..., delta4 determined from collocation method
// (constraining 3 values and 2 derivatives)
// AmpIntAnsatzFunc[]

/**
 * Ansatz for the intermediate amplitude. Equation 21 arXiv:1508.07253
 */
static double AmpIntAnsatz(double Mf, IMRPhenomDAmplitudeCoefficients* p) {
  double Mf2 = Mf*Mf;
  double Mf3 = Mf*Mf2;
  double Mf4 = Mf*Mf3;
  return p->delta0 + p->delta1*Mf + p->delta2*Mf2 + p->delta3*Mf3 + p->delta4*Mf4;
}

/**
 * The function name stands for 'Amplitude Intermediate Collocation Fit Coefficient'
 * This is the 'v2' value in Table 5 of arXiv:1508.07253
 */
static double AmpIntColFitCoeff(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 0.8149838730507785 + 2.5747553517454658*eta
  + (1.1610198035496786 - 2.3627771785551537*eta + 6.771038707057573*eta2)*xi
  + (0.7570782938606834 - 2.7256896890432474*eta + 7.1140380397149965*eta2)*xi2
  + (0.1766934149293479 - 0.7978690983168183*eta + 2.1162391502005153*eta2)*xi3;
}

  /**
  * The following functions (delta{0,1,2,3,4}_fun) were derived
  * in mathematica according to
  * the constraints detailed in arXiv:1508.07253,
  * section 'Region IIa - intermediate'.
  * These are not given in the paper.
  * Can be rederived by solving Equation 21 for the constraints
  * given in Equations 22-26 in arXiv:1508.07253
  */
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

/**
 * Calculates delta_i's
 * Method described in arXiv:1508.07253 section 'Region IIa - intermediate'
 */
static void ComputeDeltasFromCollocation(IMRPhenomDAmplitudeCoefficients* p) {
  // Three evenly spaced collocation points in the interval [f1,f3].
  double f1 = AMP_fJoin_INS;
  double f3 = p->fmaxCalc;
  double dfx = (f3 - f1)/2.0;
  double f2 = f1 + dfx;

  // v1 is inspiral model evaluated at f1
  // d1 is derivative of inspiral model evaluated at f1
  double v1 = AmpInsAnsatz(f1, p);
  double d1 = DAmpInsAnsatz(f1, p);

  // v3 is merger-ringdown model evaluated at f3
  // d2 is derivative of merger-ringdown model evaluated at f3
  double v3 = AmpMRDAnsatz(f3, p);
  double d2 = DAmpMRDAnsatz(f3, p);

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

/**
 * A struct containing all the parameters that need to be calculated
 * to compute the phenomenological amplitude
 */
static IMRPhenomDAmplitudeCoefficients* ComputeIMRPhenomDAmplitudeCoefficients(double eta, double chi1, double chi2, double finspin) {
  IMRPhenomDAmplitudeCoefficients *p = (IMRPhenomDAmplitudeCoefficients *) XLALMalloc(sizeof(IMRPhenomDAmplitudeCoefficients));

  p->eta = eta;
  p->chi1 = chi1;
  p->chi2 = chi2;

  p->q = (1.0 + sqrt(1.0 - 4.0*eta) - 2.0*eta) / (2.0*eta);
  p->chi = chiPN(eta, chi1, chi2);

  p->fRD = fring(eta, chi1, chi2, finspin);
  p->fDM = fdamp(eta, chi1, chi2, finspin);

  // Compute gamma_i's, rho_i's first then delta_i's
  p->gamma1 = gamma1_fun(eta, p->chi);
  p->gamma2 = gamma2_fun(eta, p->chi);
  p->gamma3 = gamma3_fun(eta, p->chi);

  p->fmaxCalc = fmaxCalc(p);

  p->rho1 = rho1_fun(eta, p->chi);
  p->rho2 = rho2_fun(eta, p->chi);
  p->rho3 = rho3_fun(eta, p->chi);

  // compute delta_i's
  ComputeDeltasFromCollocation(p);

  return p;
}

// Call ComputeIMRPhenomDAmplitudeCoefficients() first!
/**
 * This function computes the IMR amplitude given phenom coefficients.
 * Defined in VIII. Full IMR Waveforms arXiv:1508.07253
 */
static double IMRPhenDAmplitude(double f, IMRPhenomDAmplitudeCoefficients *p) {
  // Defined in VIII. Full IMR Waveforms arXiv:1508.07253
  // The inspiral, intermediate and merger-ringdown amplitude parts

  // Transition frequencies
  p->fInsJoin = AMP_fJoin_INS;
  p->fMRDJoin = p->fmaxCalc;

  double AmpPreFac = amp0Func(p->eta) * pow(f, -7.0/6.0);
  double AmpIns = AmpPreFac * AmpInsAnsatz(f, p);
  double AmpInt = AmpPreFac * AmpIntAnsatz(f, p);
  double AmpMRD = AmpPreFac * AmpMRDAnsatz(f, p);

  double SF1 = StepFunc(f, p->fInsJoin);
  double SF2 = StepFunc(f, p->fMRDJoin);

  return (1.0 - SF1) * AmpIns + SF1 * AmpInt * (1.0 - SF2) + SF2 * AmpMRD;  
}

/********************************* Phase functions *********************************/

////////////////////////////// Phase: Ringdown functions ///////////////////////////

// alpha_i i=1,2,3,4,5 are the phenomenological intermediate coefficients depending on eta and chiPN
// PhiRingdownAnsatz is the ringdown phasing in terms of the alpha_i coefficients

/**
 * alpha 1 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double alpha1Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 43.31514709695348 + 638.6332679188081*eta
    + (-32.85768747216059 + 2415.8938269370315*eta - 5766.875169379177*eta2)*xi
    + (-61.85459307173841 + 2953.967762459948*eta - 8986.29057591497*eta2)*xi2
    + (-21.571435779762044 + 981.2158224673428*eta - 3239.5664895930286*eta2)*xi3;
}

/**
 * alpha 2 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double alpha2Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return -0.07020209449091723 - 0.16269798450687084*eta
  + (-0.1872514685185499 + 1.138313650449945*eta - 2.8334196304430046*eta2)*xi
  + (-0.17137955686840617 + 1.7197549338119527*eta - 4.539717148261272*eta2)*xi2
  + (-0.049983437357548705 + 0.6062072055948309*eta - 1.682769616644546*eta2)*xi3;
}

/**
 * alpha 3 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double alpha3Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 9.5988072383479 - 397.05438595557433*eta
  + (16.202126189517813 - 1574.8286986717037*eta + 3600.3410843831093*eta2)*xi
  + (27.092429659075467 - 1786.482357315139*eta + 5152.919378666511*eta2)*xi2
  + (11.175710130033895 - 577.7999423177481*eta + 1808.730762932043*eta2)*xi3;
}

/**
 * alpha 4 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double alpha4Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return -0.02989487384493607 + 1.4022106448583738*eta
  + (-0.07356049468633846 + 0.8337006542278661*eta + 0.2240008282397391*eta2)*xi
  + (-0.055202870001177226 + 0.5667186343606578*eta + 0.7186931973380503*eta2)*xi2
  + (-0.015507437354325743 + 0.15750322779277187*eta + 0.21076815715176228*eta2)*xi3;
}

/**
 * alpha 5 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double alpha5Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 0.9974408278363099 - 0.007884449714907203*eta
  + (-0.059046901195591035 + 1.3958712396764088*eta - 4.516631601676276*eta2)*xi
  + (-0.05585343136869692 + 1.7516580039343603*eta - 5.990208965347804*eta2)*xi2
  + (-0.017945336522161195 + 0.5965097794825992*eta - 2.0608879367971804*eta2)*xi3;
}

/**
 * Ansatz for the merger-ringdown phase Equation 14 arXiv:1508.07253
 */
static double PhiMRDAnsatzInt(double f, IMRPhenomDPhaseCoefficients *p) {
  return -(p->alpha2/f) + (4*p->alpha3*pow(f,0.75))/3. + p->alpha1*f + p->alpha4*atan((f - p->alpha5*p->fRD)/p->fDM);
}

/**
 * First frequency derivative of PhiMRDAnsatzInt
 */
static double DPhiMRD(double f, IMRPhenomDPhaseCoefficients *p) {
  return (p->alpha1 + p->alpha2/pow(f,2) + p->alpha3/pow(f,0.25) + p->alpha4/(p->fDM*(1 + pow(f - p->alpha5*p->fRD,2)/pow(p->fDM,2)))) / p->eta;
}

///////////////////////////// Phase: Intermediate functions /////////////////////////////

// beta_i i=1,2,3 are the phenomenological intermediate coefficients depending on eta and chiPN
// PhiIntAnsatz is the intermediate phasing in terms of the beta_i coefficients


// \[Beta]1Fit = PhiIntFitCoeff\[Chi]PNFunc[\[Eta], \[Chi]PN][[1]]

/**
 * beta 1 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double beta1Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 97.89747327985583 - 42.659730877489224*eta
  + (153.48421037904913 - 1417.0620760768954*eta + 2752.8614143665027*eta2)*xi
  + (138.7406469558649 - 1433.6585075135881*eta + 2857.7418952430758*eta2)*xi2
  + (41.025109467376126 - 423.680737974639*eta + 850.3594335657173*eta2)*xi3;
}

/**
 * beta 2 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double beta2Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return -3.282701958759534 - 9.051384468245866*eta
  + (-12.415449742258042 + 55.4716447709787*eta - 106.05109938966335*eta2)*xi
  + (-11.953044553690658 + 76.80704618365418*eta - 155.33172948098394*eta2)*xi2
  + (-3.4129261592393263 + 25.572377569952536*eta - 54.408036707740465*eta2)*xi3;
}

/**
 * beta 3 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double beta3Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return -0.000025156429818799565 + 0.000019750256942201327*eta
  + (-0.000018370671469295915 + 0.000021886317041311973*eta + 0.00008250240316860033*eta2)*xi
  + (7.157371250566708e-6 - 0.000055780000112270685*eta + 0.00019142082884072178*eta2)*xi2
  + (5.447166261464217e-6 - 0.00003220610095021982*eta + 0.00007974016714984341*eta2)*xi3;
}

/**
 * ansatz for the intermediate phase defined by Equation 16 arXiv:1508.07253
 */
static double PhiIntAnsatz(double Mf, IMRPhenomDPhaseCoefficients *p) {
  // 1./eta in paper omitted and put in when need in the functions:
  // ComputeIMRPhenDPhaseConnectionCoefficients
  // IMRPhenDPhase
  return  p->beta1*Mf - p->beta3/(3.*pow(Mf,3)) + p->beta2*log(Mf);
}

/**
 * First frequency derivative of PhiIntAnsatz
 * (this time with 1/eta explicitly factored in)
 */
static double DPhiIntAnsatz(double Mf, IMRPhenomDPhaseCoefficients *p) {
  return (p->beta1 + p->beta3/pow(Mf,4) + p->beta2/Mf) / p->eta;
}

/**
 * temporary instance of DPhiIntAnsatz used when computing
 * coefficients to make the phase C(1) continuous between regions.
 */
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

/**
 * sigma 1 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double sigma1Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 2096.551999295543 + 1463.7493168261553*eta
  + (1312.5493286098522 + 18307.330017082117*eta - 43534.1440746107*eta2)*xi
  + (-833.2889543511114 + 32047.31997183187*eta - 108609.45037520859*eta2)*xi2
  + (452.25136398112204 + 8353.439546391714*eta - 44531.3250037322*eta2)*xi3;
}

/**
 * sigma 2 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double sigma2Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return -10114.056472621156 - 44631.01109458185*eta
  + (-6541.308761668722 - 266959.23419307504*eta + 686328.3229317984*eta2)*xi
  + (3405.6372187679685 - 437507.7208209015*eta + 1.6318171307344697e6*eta2)*xi2
  + (-7462.648563007646 - 114585.25177153319*eta + 674402.4689098676*eta2)*xi3;
}

/**
 * sigma 3 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double sigma3Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return 22933.658273436497 + 230960.00814979506*eta
  + (14961.083974183695 + 1.1940181342318142e6*eta - 3.1042239693052764e6*eta2)*xi
  + (-3038.166617199259 + 1.8720322849093592e6*eta - 7.309145012085539e6*eta2)*xi2
  + (42738.22871475411 + 467502.018616601*eta - 3.064853498512499e6*eta2)*xi3;
}

/**
 * sigma 4 phenom coefficient. See corresponding row in Table 5 arXiv:1508.07253
 */
static double sigma4Fit(double eta, double chi) {
  double xi = -1 + chi;
  double xi2 = xi*xi;
  double xi3 = xi2*xi;
  double eta2 = eta*eta;

  return -14621.71522218357 - 377812.8579387104*eta
  + (-9608.682631509726 - 1.7108925257214056e6*eta + 4.332924601416521e6*eta2)*xi
  + (-22366.683262266528 - 2.5019716386377467e6*eta + 1.0274495902259542e7*eta2)*xi2
  + (-85360.30079034246 - 570025.3441737515*eta + 4.396844346849777e6*eta2)*xi3;
}

/**
 * Ansatz for the inspiral phase.
 * We call the LAL TF2 coefficients here.
 * The exact values of the coefficients used are given
 * as comments in the top of this file
 * Defined by Equation 27 and 28 arXiv:1508.07253
 */
static double PhiInsAnsatzInt(double Mf, IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn) {
  double sigma1 = p->sigma1;
  double sigma2 = p->sigma2;
  double sigma3 = p->sigma3;
  double sigma4 = p->sigma4;
  double Pi = LAL_PI;

  // Assemble PN phasing series
  const double v = cbrt(Pi*Mf);
  const double logv = log(v);
  const double v2 = v * v;
  const double v3 = v * v2;
  const double v4 = v * v3;
  const double v5 = v * v4;
  const double v6 = v * v5;
  const double v7 = v * v6;

  double phasing = 0.0;
  phasing += pn->v[7] * v7;
  phasing += (pn->v[6] + pn->vlogv[6] * logv) * v6;
  phasing += (pn->v[5] + pn->vlogv[5] * logv) * v5;
  phasing += pn->v[4] * v4;
  phasing += pn->v[3] * v3;
  phasing += pn->v[2] * v2;
  phasing += pn->v[0]; // * v^0
  phasing /= v5;
  phasing -= LAL_PI_4;

  // Now add higher order terms that were calibrated for PhenomD
  phasing += (
            (sigma1/Pi) * v3
          + (3.0/(4.0*pow(Pi,4.0/3.0))*sigma2) * v4
          + (3.0/(5.0*pow(Pi,5.0/3.0))*sigma3) * v5
          + (1.0/(2.0*Pi*Pi)*sigma4) * v6
          ) / p->eta;

  return phasing;
}

/**
 * First frequency derivative of PhiInsAnsatzInt
 */
static double DPhiInsAnsatzInt(double Mf, IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn) {
  double sigma1 = p->sigma1;
  double sigma2 = p->sigma2;
  double sigma3 = p->sigma3;
  double sigma4 = p->sigma4;
  double Pi = LAL_PI;

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
  double Dphasing = 0.0;
  Dphasing += +2.0 * pn->v[7] * v7;
  Dphasing += (pn->v[6] + pn->vlogv[6] * (1.0 + logv)) * v6;
  Dphasing += pn->vlogv[5] * v5;
  Dphasing += -1.0 * pn->v[4] * v4;
  Dphasing += -2.0 * pn->v[3] * v3;
  Dphasing += -3.0 * pn->v[2] * v2;
  Dphasing += -5.0 * pn->v[0];
  Dphasing /= v8 * 3.0/Pi;

  // Now add higher order terms that were calibrated for PhenomD
  Dphasing += (
          sigma1
        + (pow(Pi, -1.0/3.0)*sigma2) * v
        + (pow(Pi, -2.0/3.0)*sigma3) * v2
        + (sigma4/Pi) * v3
        ) / p->eta;

  return Dphasing;
}

///////////////////////////////// Phase: glueing function ////////////////////////////////

/**
 * A struct containing all the parameters that need to be calculated
 * to compute the phenomenological phase
 */
static IMRPhenomDPhaseCoefficients* ComputeIMRPhenomDPhaseCoefficients(double eta, double chi1, double chi2, double finspin, const LALSimInspiralTestGRParam *extraParams) {

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

  p->fRD = fring(eta, chi1, chi2, finspin);
  p->fDM = fdamp(eta, chi1, chi2, finspin);

  if (extraParams!=NULL)
    {
      if (XLALSimInspiralTestGRParamExists(extraParams,"dsigma1")) p->sigma1*=(1.0+XLALSimInspiralGetTestGRParam(extraParams,"dsigma1"));
      if (XLALSimInspiralTestGRParamExists(extraParams,"dsigma2")) p->sigma2*=(1.0+XLALSimInspiralGetTestGRParam(extraParams,"dsigma2"));
      if (XLALSimInspiralTestGRParamExists(extraParams,"dsigma3")) p->sigma3*=(1.0+XLALSimInspiralGetTestGRParam(extraParams,"dsigma3"));
      if (XLALSimInspiralTestGRParamExists(extraParams,"dsigma4")) p->sigma4*=(1.0+XLALSimInspiralGetTestGRParam(extraParams,"dsigma4"));
      if (XLALSimInspiralTestGRParamExists(extraParams,"dbeta1")) p->beta1*=(1.0+XLALSimInspiralGetTestGRParam(extraParams,"dbeta1"));
      if (XLALSimInspiralTestGRParamExists(extraParams,"dbeta2")) p->beta2*=(1.0+XLALSimInspiralGetTestGRParam(extraParams,"dbeta2"));
      if (XLALSimInspiralTestGRParamExists(extraParams,"dbeta3")) p->beta3*=(1.0+XLALSimInspiralGetTestGRParam(extraParams,"dbeta3"));
      if (XLALSimInspiralTestGRParamExists(extraParams,"dalpha1")) p->alpha1*=(1.0+XLALSimInspiralGetTestGRParam(extraParams,"dalpha1"));
      if (XLALSimInspiralTestGRParamExists(extraParams,"dalpha2")) p->alpha2*=(1.0+XLALSimInspiralGetTestGRParam(extraParams,"dalpha2"));
      if (XLALSimInspiralTestGRParamExists(extraParams,"dalpha3")) p->alpha3*=(1.0+XLALSimInspiralGetTestGRParam(extraParams,"dalpha3"));
      if (XLALSimInspiralTestGRParamExists(extraParams,"dalpha4")) p->alpha4*=(1.0+XLALSimInspiralGetTestGRParam(extraParams,"dalpha4"));
      if (XLALSimInspiralTestGRParamExists(extraParams,"dalpha5")) p->alpha5*=(1.0+XLALSimInspiralGetTestGRParam(extraParams,"dalpha5"));
    }

  return p;
}

/**
 * This function aligns the three phase parts (inspiral, intermediate and merger-rindown)
 * such that they are c^1 continuous at the transition frequencies
 * Defined in VIII. Full IMR Waveforms arXiv:1508.07253
 */
static void ComputeIMRPhenDPhaseConnectionCoefficients(IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn) {
  double eta = p->eta;

  // Transition frequencies
  // Defined in VIII. Full IMR Waveforms arXiv:1508.07253
  p->fInsJoin=PHI_fJoin_INS;
  p->fMRDJoin=0.5*p->fRD;

  // Compute C1Int and C2Int coeffs
  // Equations to solve for to get C(1) continuous join
  // PhiIns (f)  =   PhiInt (f) + C1Int + C2Int f
  // Joining at fInsJoin
  // PhiIns (fInsJoin)  =   PhiInt (fInsJoin) + C1Int + C2Int fInsJoin
  // PhiIns'(fInsJoin)  =   PhiInt'(fInsJoin) + C2Int
  double DPhiIns = DPhiInsAnsatzInt(p->fInsJoin, p, pn);
  double DPhiInt = DPhiIntAnsatz(p->fInsJoin, p);
  p->C2Int = DPhiIns - DPhiInt;
  p->C1Int = PhiInsAnsatzInt(p->fInsJoin, p, pn)
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

/**
 * This function computes the IMR phase given phenom coefficients.
 * Defined in VIII. Full IMR Waveforms arXiv:1508.07253
 */
static double IMRPhenDPhase(double f, IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn) {
  // Defined in VIII. Full IMR Waveforms arXiv:1508.07253
  // The inspiral, intermendiate and merger-ringdown phase parts
  double PhiIns = PhiInsAnsatzInt(f, p, pn);
  double PhiInt = 1.0/p->eta * PhiIntAnsatz(f, p) + p->C1Int + p->C2Int * f;
  double PhiMRD = 1.0/p->eta * PhiMRDAnsatzInt(f, p) + p->C1MRD + p->C2MRD * f;

  double SF1 = StepFunc(f, p->fInsJoin);
  double SF2 = StepFunc(f, p->fMRDJoin);

  return (1.0 - SF1) * PhiIns + SF1 * PhiInt * (1.0 - SF2) + SF2 * PhiMRD;
}
