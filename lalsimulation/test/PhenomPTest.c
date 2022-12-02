/*
 *  Copyright (C) 2013 Michael Puerrer
 *  Test code for LALSimIMRPhenomP(v1)
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <lal/Units.h>
#include <lal/LALAdaptiveRungeKuttaIntegrator.h>
#include <lal/LALConstants.h>
#include <lal/FindRoot.h>
#include <lal/SeqFactories.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include <lal/LALSimNoise.h>
#include <lal/ComplexFFT.h>

#include <lal/ComplexFFT.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>

#include <LALSimIMRPhenomP.c> /* Include source directly so we can carry out unit tests for internal functions */

#define MYUNUSED(expr) do { (void)(expr); } while (0)

void prC(const char *name, COMPLEX16 x);
void prC(const char *name, COMPLEX16 x) {
  printf("%s: %g + I %g\n", name, creal(x), cimag(x));
}

REAL8 sqr(REAL8 x);
REAL8 sqr(REAL8 x) { return x*x; }

REAL8 MatchSI(COMPLEX16FrequencySeries **htilde1, COMPLEX16FrequencySeries **htilde2, REAL8 fMin, REAL8 fMax, REAL8 df);
REAL8 MatchSI(COMPLEX16FrequencySeries **htilde1, COMPLEX16FrequencySeries **htilde2, REAL8 fMin, REAL8 fMax, REAL8 df) {
  // Frequencies are in Hz!

  // Assume that waveforms use same frequency points
  int len1 = (*htilde1)->data->length;
  int len2 = (*htilde1)->data->length;
  if (len1 != len2) {
    XLALPrintError("Length of waveforms differs!\n");
    XLAL_ERROR(XLAL_EDOM); // FIXME
  }
  int n = len1;

  COMPLEX16Vector *integrand = XLALCreateCOMPLEX16Vector(n);
  REAL8 PSDfact = XLALSimNoisePSDaLIGOZeroDetHighPower((fMax-fMin)/4.);
  REAL8 tableS;
  COMPLEX16 h1,h2;

  REAL8 norm1 = 0.0;
  REAL8 norm2 = 0.0;
  int iStart = (int)(fMin / df);
  int iStop = (int)(fMax / df);
  iStart = (iStart < 1) ? 1 : iStart;
  iStop = (iStop > n) ? n : iStop;
  for (int i=iStart; i<iStop; i++) {
    REAL8 f = i*df;
    tableS = PSDfact / XLALSimNoisePSDaLIGOZeroDetHighPower(f);
    h1 = ((*htilde1)->data->data)[i];
    h2 = ((*htilde2)->data->data)[i];
    integrand->data[i] = h1 * conj(h2) * tableS;
    norm1 += sqr(cabs(h1)) * tableS;
    norm2 += sqr(cabs(h2)) * tableS;
    // printf("f = %g\tnoise(f) = %g\ttableS[i] = %g\tintegrand[i] = %g\n", f, XLALSimNoisePSDaLIGOZeroDetHighPower(f), tableS, creal(integrand[i]));
    // printf("{norm1, norm2} = {%g,%g}\t{%g,%g}\n", norm1,norm2, cabs(h1)*tableS, cabs(h2)*tableS);
  }

  n = iStop-iStart;
  int zpf = 10;
  int m = n + 2*zpf*n; // zero-pad on both sides
  //printf("Total length %d\n", m);

  COMPLEX16FFTPlan *myplan = XLALCreateForwardCOMPLEX16FFTPlan(m, 0);
  COMPLEX16Vector *array_in = XLALCreateCOMPLEX16Vector(m);
  COMPLEX16Vector *array_out = XLALCreateCOMPLEX16Vector(m);

  // fill input array
  for (int i=0; i<zpf*n; i++)
    array_in->data[i] = 0.0;
  for (int i=zpf*n; i<(zpf*n + n); i++)
    array_in->data[i] = integrand->data[iStart + i-zpf*n];
  for (int i=zpf*n + n; i<m; i++)
    array_in->data[i] = 0.0;

  XLALCOMPLEX16VectorFFT(array_out, array_in, myplan);

  REAL8 match=0; REAL8 val;
  for (int i=0; i<m; i++) {
    val = cabs(array_out->data[i]);
    if (val > match)
      match = val;
  }

  XLALFree(array_in);
  XLALFree(array_out);
  XLALDestroyCOMPLEX16FFTPlan(myplan);

  return match / sqrt(norm1*norm2);
}

void dump_file(const char *filename, COMPLEX16FrequencySeries *hptilde, COMPLEX16FrequencySeries *hctilde, REAL8 M);
void dump_file(const char *filename, COMPLEX16FrequencySeries *hptilde, COMPLEX16FrequencySeries *hctilde, REAL8 M) {
  COMPLEX16 hp;
  COMPLEX16 hc;
  REAL8 f;
  FILE *out;
  REAL8 f_max_prime = 0;
  int wflen = hptilde->data->length;
  REAL8 deltaF = hptilde->deltaF;
  out = fopen(filename, "w");
  for (int i=1;i<wflen;i++) {
    f = i*deltaF;
    hp = (hptilde->data->data)[i];
    hc = (hctilde->data->data)[i];
    if (!(creal(hp) == 0 && cimag(hp) == 0 && creal(hc) == 0 && cimag(hc) == 0)) {
      f_max_prime = (f > f_max_prime) ? f : f_max_prime;
      fprintf(out, "%.15g\t%.15g\t%.15g\t%.15g\t%.15g\n", f*LAL_MTSUN_SI*M, creal(hp), cimag(hp), creal(hc), cimag(hc));
    }
  }
  fclose(out);
}

bool approximatelyEqual(REAL8 a, REAL8 b, REAL8 epsilon);
bool approximatelyEqualC(COMPLEX16 a, COMPLEX16 b, REAL8 epsilon);

bool approximatelyEqual(REAL8 a, REAL8 b, REAL8 epsilon) {
  return !gsl_fcmp(a, b, epsilon); // gsl_fcmp() returns 0 if the numbers a, b are approximately equal to a relative accuracy epsilon.
}

bool approximatelyEqualC(COMPLEX16 a, COMPLEX16 b, REAL8 epsilon) {
  return approximatelyEqual(creal(a), creal(b), epsilon) && approximatelyEqual(cimag(a), cimag(b), epsilon);
}

void print_difference(const char *name, REAL8 u, REAL8 u_expected);
void print_difference(const char *name, REAL8 u, REAL8 u_expected) {
  printf("%-8s: %-20.17g\t%-20.17g\t%-20.17g\n", name, u, u_expected, u - u_expected);
}

static void Test_alpha_epsilon(void);
static void Test_alpha_epsilon(void) {
  printf("\n** Test_alpha_epsilon: **\n");
  const REAL8 f = 0.01;
  const REAL8 q = 4;
  const REAL8 chil = 0.5625;
  const REAL8 chip = 0.18;

  NNLOanglecoeffs angcoeffs;
  ComputeNNLOanglecoeffs(&angcoeffs,q,chil,chip);

  const REAL8 omega = LAL_PI * f;
  const REAL8 logomega = log(omega);
  const REAL8 omega_cbrt = cbrt(omega);
  const REAL8 omega_cbrt2 = omega_cbrt*omega_cbrt;
  const REAL8 alpha = (angcoeffs.alphacoeff1/omega
                    + angcoeffs.alphacoeff2/omega_cbrt2
                    + angcoeffs.alphacoeff3/omega_cbrt
                    + angcoeffs.alphacoeff4*logomega
                    + angcoeffs.alphacoeff5*omega_cbrt);

  const REAL8 epsilon = (angcoeffs.epsiloncoeff1/omega
                      + angcoeffs.epsiloncoeff2/omega_cbrt2
                      + angcoeffs.epsiloncoeff3/omega_cbrt
                      + angcoeffs.epsiloncoeff4*logomega
                      + angcoeffs.epsiloncoeff5*omega_cbrt);

  const REAL8 alpha_expected = -11.8195574;
  const REAL8 epsilon_expected = -11.9359726;

  print_difference("alpha", alpha, alpha_expected);
  print_difference("epsilon", epsilon, epsilon_expected);

  const REAL8 eps = 1e-5;

  XLAL_CHECK_EXIT(
       approximatelyEqual(alpha,    alpha_expected,   eps)
    && approximatelyEqual(epsilon,  epsilon_expected, eps)
    && "Test_alpha_epsilon()"
  );
}

static void Test_XLALSimIMRPhenomPCalculateModelParameters(void);
static void Test_XLALSimIMRPhenomPCalculateModelParameters(void) {
  printf("\n** Test_XLALSimIMRPhenomPCalculateModelParameters: **\n");

  REAL8 chi1_l, chi2_l, chi_eff, chip, thetaJ, alpha0;

  REAL8 m1_SI = 10 * LAL_MSUN_SI;
  REAL8 m2_SI = 40 * LAL_MSUN_SI;
  REAL8 s1x = 0.3;
  REAL8 s1y = 0;
  REAL8 s1z = 0.45;
  REAL8 s2x = 0;
  REAL8 s2y = 0;
  REAL8 s2z = 0.45;
  REAL8 lnhatx = sin(0.4);
  REAL8 lnhaty = 0;
  REAL8 lnhatz = cos(0.4);
  REAL8 f_min = 20;
  REAL8 f_ref = f_min;
  IMRPhenomP_version_type version = IMRPhenomPv2_V;

  XLALSimIMRPhenomPCalculateModelParametersOld(
      &chi1_l,            /**< Output: aligned spin on companion 1 */
      &chi2_l,            /**< Output: aligned spin on companion 2 */
      &chip,              /**< Output: Effective spin in the orbital plane */
      &thetaJ,            /**< Output: Angle between J0 and line of sight (z-direction) */
      &alpha0,            /**< Output: Initial value of alpha angle (azimuthal precession angle) */
      m1_SI,              /**< Mass of companion 1 (kg) */
      m2_SI,              /**< Mass of companion 2 (kg) */
      f_ref,              /**< Reference GW frequency (Hz) */
      lnhatx,             /**< Initial value of LNhatx: orbital angular momentum unit vector */
      lnhaty,             /**< Initial value of LNhaty */
      lnhatz,             /**< Initial value of LNhatz */
      s1x,                /**< Initial value of s1x: dimensionless spin of larger BH */
      s1y,                /**< Initial value of s1y: dimensionless spin of larger BH */
      s1z,                /**< Initial value of s1z: dimensionless spin of larger BH */
      s2x,                /**< Initial value of s2x: dimensionless spin of larger BH */
      s2y,                /**< Initial value of s2y: dimensionless spin of larger BH */
      s2z,                /**< Initial value of s2z: dimensionless spin of larger BH */
      version);

  chi_eff = (m1_SI*chi1_l + m2_SI*chi2_l) / (m1_SI + m2_SI); /* Effective aligned spin */

  REAL8 chi_eff_expected = 0.4378425478398173;
  REAL8 chip_expected = 0.1752382540388927;
  REAL8 thetaJ_expected = 0.29197409372473093;
  REAL8 alpha0_expected = LAL_PI;

  print_difference("chi_eff", chi_eff, chi_eff_expected);
  print_difference("chip", chip, chip_expected);
  print_difference("thetaJ", thetaJ, thetaJ_expected);
  print_difference("alpha0", alpha0, alpha0_expected);

  //const REAL8 eps = DBL_EPSILON;
  const REAL8 eps = 1e-5;

  XLAL_CHECK_EXIT(approximatelyEqual(chi_eff,  chi_eff_expected, eps)
    && approximatelyEqual(chip,     chip_expected, eps)
    && approximatelyEqual(thetaJ,   thetaJ_expected, eps)
    && approximatelyEqual(alpha0,   alpha0_expected, eps)
    && "Test_XLALSimIMRPhenomPCalculateModelParameters()"
  );
}


#if 0
static void Test_PhenomC(void);
static void Test_PhenomC(void) {
  printf("\n** Test_PhenomC: **\n");
  const REAL8 fHz = 40.6051;
  const REAL8 eta = 0.16;
  const REAL8 m1 = 10;
  const REAL8 m2 = 40;
  const REAL8 chi = 0.45;
  const REAL8 chip = 0.18;
  const REAL8 M = (m1+m2);
  const REAL8 distance = 100 * 1e6 * LAL_PC_SI;
  const REAL8 phic = 0;
  REAL8 phPhenomC = 0.0;
  REAL8 aPhenomC  = 0.0;

  BBHPhenomCParams *PCparams = ComputeIMRPhenomCParamsRDmod(
    m1,    /**< mass of companion 1 (solar masses) */
    m2,    /**< mass of companion 2 (solar masses) */
    chi,   /**< Reduced aligned spin of the binary chi = (m1*chi1 + m2*chi2)/M */
    chip, /**< Dimensionless spin in the orbial plane */
    NULL); /**< No testing GR parameters */

  int errcode = IMRPhenomCGenerateAmpPhase( &aPhenomC, &phPhenomC, fHz, eta, PCparams );
  if( errcode != XLAL_SUCCESS )
    exit(-1);
  phPhenomC -= 2.*phic; /* Note: phic is orbital phase */
  REAL8 amp0 = 2. * sqrt(5. / (64.*LAL_PI)) * M * LAL_MRSUN_SI * M * LAL_MTSUN_SI / distance;
  printf("test %g\n", M* LAL_MRSUN_SI * M* LAL_MTSUN_SI / distance);
  COMPLEX16 hPC = amp0 * aPhenomC * cexp(-I*phPhenomC); /* Assemble PhenomC waveform. */
  printf("phPhenomC: %g\n", phPhenomC);
  printf("aPhenomC: %g\tamp0: %g\n", aPhenomC, amp0);
  printf("LAL_MRSUN_SI, LAL_MTSUN_SI, LAL_PC_SI: %g\t%g\t%g\n",LAL_MRSUN_SI, LAL_MTSUN_SI, LAL_PC_SI);
  prC("hPC", hPC);

  COMPLEX16 hPC_expected = -4.08291e-23 - I * 8.89596e-23;

  const REAL8 eps = 1e-5;

  XLAL_CHECK_EXIT(
       approximatelyEqualC(hPC,      hPC_expected, eps)
    && "Test_PhenomC()"
  );
}
#endif

#if 0
static void Test_PhenomPCore(void);
static void Test_PhenomPCore(void) {
  printf("\n** Test_PhenomPCore: **\n");
  BBHPhenomCParams *PCparams = ComputeIMRPhenomCParamsRDmod(10, 40, 0.45, 0.18, NULL);
  REAL8 q = 4;
  REAL8 chi_eff = 0.45;
  const REAL8 chil = (1.0+q)/q * chi_eff; /* dimensionless aligned spin of the largest BH */
  printf("chil %g\n", chil);
  REAL8 chip = 0.18;
  NNLOanglecoeffs angcoeffs;
  ComputeNNLOanglecoeffs(&angcoeffs,q,chil,chip);
  const REAL8 m_sec = 50 * LAL_MTSUN_SI;   /* Total mass in seconds */
  const REAL8 piM = LAL_PI * m_sec;
  const REAL8 omega_ref = piM * 20;
  const REAL8 logomega_ref = log(omega_ref);
  const REAL8 omega_ref_cbrt = cbrt(omega_ref);
  const REAL8 omega_ref_cbrt2 = omega_ref_cbrt*omega_ref_cbrt;
  const REAL8 alphaNNLOoffset = (angcoeffs.alphacoeff1/omega_ref
                              + angcoeffs.alphacoeff2/omega_ref_cbrt2
                              + angcoeffs.alphacoeff3/omega_ref_cbrt
                              + angcoeffs.alphacoeff4*logomega_ref
                              + angcoeffs.alphacoeff5*omega_ref_cbrt);
  printf("alphaNNLOoffset %g\n", alphaNNLOoffset);

  SpinWeightedSphericalHarmonic_l2 Y2m;
  const REAL8 ytheta  = 0.279523;
  const REAL8 yphi    = alphaNNLOoffset - 0.0234206;
  printf("yphi %g\n", alphaNNLOoffset - 0.0234206);
  Y2m.Y2m2 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2, -2);
  Y2m.Y2m1 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2, -1);
  Y2m.Y20  = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2,  0);
  Y2m.Y21  = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2,  1);
  Y2m.Y22  = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2,  2);

  COMPLEX16 hp, hc;
  REAL8 phasing;
  REAL8 fHz = 40.6051; // Mf = 0.01 for M=50Msun
  const UINT4 version = 1;
  IMRPhenomDAmplitudeCoefficients *pAmp = NULL;
  IMRPhenomDPhaseCoefficients *pPhi = NULL;
  PNPhasingSeries *PNparams = NULL;
  AmpInsPrefactors * amp_prefactors = NULL;
  PhiInsPrefactors * phi_prefactors = NULL;

  int ret = PhenomPCoreOneFrequency(
    fHz,                     /**< frequency (Hz) */
    0.16,                    /**< symmetric mass ratio */
    0.45,                    /**< dimensionless effective total aligned spin */
    0.18,                    /**< dimensionless spin in the orbial plane */
    100 * 1e6 * LAL_PC_SI,   /**< distance of source (m) */
    50,                      /**< total mass (Solar masses) */
    0,                       /**< orbital coalescence phase (rad) */
    pAmp,                    /**< Internal IMRPhenomD amplitude coefficients */
    pPhi,                    /**< Internal IMRPhenomD phase coefficients */
    PCparams,                /**< internal PhenomC parameters */
    PNparams,                /**< PN inspiral phase coefficients */
    &angcoeffs,              /**< struct with PN coeffs for the NNLO angles */
    &Y2m,                    /**< struct of l=2 spherical harmonics of spin weight -2 */
    0,0,
    &hp,                     /**< output: \f$\tilde h_+\f$ */
    &hc,                     /**< output: \f$\tilde h_+\f$ */
    &phasing,                /**< Output: overall phasing */
    version,
    amp_prefactors,
    phi_prefactors
);

  MYUNUSED(ret);
  prC("hp", hp);
  prC("hc", hc);

  COMPLEX16 hp_expected = 2.06975e-23 - I * 9.29353e-23;
  COMPLEX16 hc_expected = -9.29441e-23 - I * 2.06616e-23;
  const REAL8 eps = 1e-5;

  XLAL_CHECK_EXIT(
       approximatelyEqualC(hp, hp_expected, eps)
    && approximatelyEqualC(hc, hc_expected, eps)
    && "Test_PhenomPCore()"
  );
}
#endif

static void Test_XLALSimIMRPhenomP(void);
static void Test_XLALSimIMRPhenomP(void) {
  printf("\n** Test_XLALSimIMRPhenomP: **\n");

  REAL8 chi1_l, chi2_l, chip, thetaJ, alpha0;

  REAL8 m1 = 10;
  REAL8 m2 = 40;
  REAL8 m1_SI = m1 * LAL_MSUN_SI;
  REAL8 m2_SI = m2 * LAL_MSUN_SI;
  REAL8 s1x = 0.3;
  REAL8 s1y = 0;
  REAL8 s1z = 0.45;
  REAL8 s2x = 0;
  REAL8 s2y = 0;
  REAL8 s2z = 0.45;
  REAL8 lnhatx = sin(0.4);
  REAL8 lnhaty = 0;
  REAL8 lnhatz = cos(0.4);
  REAL8 f_min = 20;
  REAL8 f_ref = f_min;
  IMRPhenomP_version_type version = IMRPhenomPv2_V;
  NRTidal_version_type NRTv = NoNRT_V;

  XLALSimIMRPhenomPCalculateModelParametersOld(
      &chi1_l,            /**< Output: aligned spin on companion 1 */
      &chi2_l,            /**< Output: aligned spin on companion 2 */
      &chip,              /**< Output: Effective spin in the orbital plane */
      &thetaJ,            /**< Output: Angle between J0 and line of sight (z-direction) */
      &alpha0,            /**< Output: Initial value of alpha angle (azimuthal precession angle) */
      m1_SI,              /**< Mass of companion 1 (kg) */
      m2_SI,              /**< Mass of companion 2 (kg) */
      f_ref,              /**< Reference GW frequency (Hz) */
      lnhatx,             /**< Initial value of LNhatx: orbital angular momentum unit vector */
      lnhaty,             /**< Initial value of LNhaty */
      lnhatz,             /**< Initial value of LNhatz */
      s1x,                /**< Initial value of s1x: dimensionless spin of larger BH */
      s1y,                /**< Initial value of s1y: dimensionless spin of larger BH */
      s1z,                /**< Initial value of s1z: dimensionless spin of larger BH */
      s2x,                /**< Initial value of s2x: dimensionless spin of larger BH */
      s2y,                /**< Initial value of s2y: dimensionless spin of larger BH */
      s2z,                /**< Initial value of s2z: dimensionless spin of larger BH */
      version);

  COMPLEX16FrequencySeries *hptilde = NULL;
  COMPLEX16FrequencySeries *hctilde = NULL;
  REAL8 phic = 0;
  REAL8 deltaF = 0.06;
  REAL8 f_max = 0; // 8000;
  REAL8 distance = 100 * 1e6 * LAL_PC_SI;

  int ret = XLALSimIMRPhenomP(
    &hptilde,                 /**< Frequency-domain waveform h+ */
    &hctilde,                 /**< Frequency-domain waveform hx */
    chi1_l,                   /**< aligned spin on companion 1 */
    chi2_l,                   /**< aligned spin on companion 2 */
    chip,                     /**< Effective spin in the orbital plane */
    thetaJ,                   /**< Angle between J0 and line of sight (z-direction) */
    m1_SI,                    /**< mass of companion1 (kg) */
    m2_SI,                    /**< mass of companion1 (kg) */
    distance,                 /**< Distance of source (m) */
    alpha0,                   /**< Initial value of alpha angle */
    phic,                     /**< Orbital coalescence phase (rad) */
    deltaF,                   /**< Sampling frequency (Hz) */
    f_min,                    /**< Starting GW frequency (Hz) */
    f_max,                    /**< End frequency; 0 defaults to ringdown cutoff freq */
    f_ref,                    /**< Reference frequency */
    version,
    NRTv,
    NULL);                   /**<linked list containing the extra testing GR parameters */

  dump_file("PhenomP_Test1.dat", hptilde, hctilde, m1+m2);
  MYUNUSED(ret);
  COMPLEX16 hp = (hptilde->data->data)[1000];
  COMPLEX16 hc = (hctilde->data->data)[1000];
  prC("hp", hp);
  prC("hc", hc);

  COMPLEX16 hp_expected = -5.17642e-23 + I * 2.60463e-23;
  COMPLEX16 hc_expected =  2.6046e-23 + I * 5.17592e-23;
  const REAL8 eps = 1e-5;

  XLAL_CHECK_EXIT(
       approximatelyEqualC(hp, hp_expected, eps)
    && approximatelyEqualC(hc, hc_expected, eps)
    && "XLALSimIMRPhenomP()"
  );

}

#if 0
static void Test_PhenomC_PhenomP(void);
static void Test_PhenomC_PhenomP(void) {
  printf("\n** Test_PhenomC_PhenomP: **\n");
  REAL8 eta = 0.16;
  REAL8 chi = 0.45;
  REAL8 q = (1 + sqrt(1 - 4*eta) - 2*eta)/(2.*eta);
  REAL8 M = 50;

  // Parameters for XLALSimIMRPhenomPCalculateModelParameters
  COMPLEX16FrequencySeries *hptilde = NULL;
  COMPLEX16FrequencySeries *hctilde = NULL;
  REAL8 phic = 0;
  REAL8 deltaF = 0.06;
  REAL8 m1_SI = M *   q/(1+q) * LAL_MSUN_SI;
  REAL8 m2_SI = M * 1.0/(1+q) * LAL_MSUN_SI;
  REAL8 f_min = 10;
  REAL8 f_ref = f_min;
  REAL8 f_max = 0; // 8000;
  REAL8 distance = 100 * 1e6 * LAL_PC_SI;
  REAL8 s1x = 0;
  REAL8 s1y = 0;
  REAL8 s1z = chi;
  REAL8 s2x = 0;
  REAL8 s2y = 0;
  REAL8 s2z = chi;
  REAL8 lnhatx = 0;
  REAL8 lnhaty = 0;
  REAL8 lnhatz = 1;
  const UINT4 version = 1;
  const LALSimInspiralTestGRParam *nonGR = NULL;

  REAL8 chi1_l, chi2_l, chip, thetaJ, alpha0;

  XLALSimIMRPhenomPCalculateModelParameters(
      &chi1_l,            /**< Output: aligned spin on companion 1 */
      &chi2_l,            /**< Output: aligned spin on companion 2 */
      &chip,              /**< Output: Effective spin in the orbital plane */
      &thetaJ,            /**< Output: Angle between J0 and line of sight (z-direction) */
      &alpha0,            /**< Output: Initial value of alpha angle (azimuthal precession angle) */
      m1_SI,              /**< Mass of companion 1 (kg) */
      m2_SI,              /**< Mass of companion 2 (kg) */
      f_ref,              /**< Starting GW frequency (Hz) */
      lnhatx,             /**< Initial value of LNhatx: orbital angular momentum unit vector */
      lnhaty,             /**< Initial value of LNhaty */
      lnhatz,             /**< Initial value of LNhatz */
      s1x,                /**< Initial value of s1x: dimensionless spin of larger BH */
      s1y,                /**< Initial value of s1y: dimensionless spin of larger BH */
      s1z,                /**< Initial value of s1z: dimensionless spin of larger BH */
      s2x,                /**< Initial value of s2x: dimensionless spin of larger BH */
      s2y,                /**< Initial value of s2y: dimensionless spin of larger BH */
      s2z);               /**< Initial value of s2z: dimensionless spin of larger BH */

//  printf("chi_eff = %g\n", chi_eff);
  printf("chip = %g\n", chip);
  printf("eta = %g\n", eta);
  printf("thetaJ = %g\n", thetaJ);

  int ret = XLALSimIMRPhenomP(
    &hptilde,                 /**< Frequency-domain waveform h+ */
    &hctilde,                 /**< Frequency-domain waveform hx */
    chi1_l,                   /**< aligned spin on companion 1 */
    chi2_l,                   /**< aligned spin on companion 2 */
    chip,                     /**< Effective spin in the orbital plane */
    thetaJ,                   /**< Angle between J0 and line of sight (z-direction) */
    m1_SI,                    /**< Mass of companion 1 (kg) */
    m2_SI,                    /**< Mass of companion 2 (kg) */
    distance,                 /**< Distance of source (m) */
    alpha0,                   /**< Initial value of alpha angle (azimuthal precession angle) */
    phic,                     /**< Orbital coalescence phase (rad) */
    deltaF,                   /**< Sampling frequency (Hz) */
    f_min,                    /**< Starting GW frequency (Hz) */
    f_max,                    /**< End frequency; 0 defaults to ringdown cutoff freq */
    f_ref,                    /**< Reference frequency */
    version,
    nonGR);

  int wflen = hptilde->data->length;
  REAL8 f_max_prime = 0;
  FILE *out;
  out = fopen("XLALSimIMRPhenomP.dat", "w");
  for (int i=1;i<wflen;i++) {
    REAL8 f = i*deltaF;
    COMPLEX16 hp = (hptilde->data->data)[i];
    COMPLEX16 hc = (hctilde->data->data)[i];
    if (!(creal(hp) == 0 && cimag(hp) == 0 && creal(hc) == 0 && cimag(hc) == 0)) {
      f_max_prime = (f > f_max_prime) ? f : f_max_prime;
      fprintf(out, "%.15g\t%.15g\t%.15g\t%.15g\t%.15g\n", f*LAL_MTSUN_SI*M, creal(hp), cimag(hp), creal(hc), cimag(hc));
    }
  }
  fclose(out);
  printf("f_max_prime = %g\n", f_max_prime);


  COMPLEX16FrequencySeries *htildePC = NULL;
  ret = XLALSimIMRPhenomCGenerateFD(
      &htildePC,             /**< FD waveform */
      phic,                  /**< orbital phase at peak (rad) */
      deltaF,                /**< sampling interval (Hz) */
      m1_SI,                 /**< mass of companion 1 (kg) */
      m2_SI,                 /**< mass of companion 2 (kg) */
      chi,                   /**< mass-weighted aligned-spin parameter */
      f_min,                 /**< starting GW frequency (Hz) */
      f_max,                 /**< end frequency; 0 defaults to ringdown cutoff freq */
      distance,               /**< distance of source (m) */
      nonGR
  );

  out = fopen("XLALSimIMRPhenomC.dat", "w");
  wflen = hptilde->data->length;
  for (int i=1;i<wflen;i++) {
    REAL8 f = i*deltaF;
    COMPLEX16 hp = (htildePC->data->data)[i];
    if (!(creal(hp) == 0 && cimag(hp) == 0))
      fprintf(out, "%.15g\t%.15g\t%.15g\n", f*LAL_MTSUN_SI*M, creal(hp), cimag(hp));
  }
  fclose(out);

  // Now compute match between PhenomC and PhenomP for this aligned configuration
  REAL8 match = MatchSI(&hptilde, &htildePC, f_min, f_max_prime, deltaF);
  REAL8 match_expected = 0.999465;

  const REAL8 eps = 1e-5;

  XLAL_CHECK_EXIT(
       approximatelyEqualC(match, match_expected, eps)
    && "Test_PhenomC_PhenomP()"
  );

  printf("match(PhenomP_aligned, PhenomC) = %g\n", match);
  MYUNUSED(ret);
}
#endif

static void Test_XLALSimIMRPhenomP_f_ref(void);
static void Test_XLALSimIMRPhenomP_f_ref(void) {
  // For aligned spins f_ref should not change the waveform
  printf("\n** Test_XLALSimIMRPhenomP_f_ref: **\n");

  REAL8 chi1_l, chi2_l, chip, thetaJ, alpha0;

  REAL8 m1 = 10;
  REAL8 m2 = 40;
  REAL8 m1_SI = m1 * LAL_MSUN_SI;
  REAL8 m2_SI = m2 * LAL_MSUN_SI;
  REAL8 s1x = 0.45*sin(0.4);
  REAL8 s1y = 0;
  REAL8 s1z = 0.45*cos(0.4);
  REAL8 s2x = 0.45*sin(0.4);
  REAL8 s2y = 0;
  REAL8 s2z = 0.45*cos(0.4);
  REAL8 lnhatx = sin(0.4);
  REAL8 lnhaty = 0;
  REAL8 lnhatz = cos(0.4);
  REAL8 f_min = 20;
  REAL8 f_ref = f_min;
  IMRPhenomP_version_type version = IMRPhenomPv2_V;
  NRTidal_version_type NRTv = NoNRT_V;

  XLALSimIMRPhenomPCalculateModelParametersOld(
      &chi1_l,            /**< Output: aligned spin on companion 1 */
      &chi2_l,            /**< Output: aligned spin on companion 2 */
      &chip,              /**< Output: Effective spin in the orbital plane */
      &thetaJ,            /**< Output: Angle between J0 and line of sight (z-direction) */
      &alpha0,            /**< Output: Initial value of alpha angle (azimuthal precession angle) */
      m1_SI,              /**< Mass of companion 1 (kg) */
      m2_SI,              /**< Mass of companion 2 (kg) */
      f_ref,              /**< Reference GW frequency (Hz) */
      lnhatx,             /**< Initial value of LNhatx: orbital angular momentum unit vector */
      lnhaty,             /**< Initial value of LNhaty */
      lnhatz,             /**< Initial value of LNhatz */
      s1x,                /**< Initial value of s1x: dimensionless spin of larger BH */
      s1y,                /**< Initial value of s1y: dimensionless spin of larger BH */
      s1z,                /**< Initial value of s1z: dimensionless spin of larger BH */
      s2x,                /**< Initial value of s2x: dimensionless spin of larger BH */
      s2y,                /**< Initial value of s2y: dimensionless spin of larger BH */
      s2z,                /**< Initial value of s2z: dimensionless spin of larger BH */
      version);

  COMPLEX16FrequencySeries *hptilde = NULL;
  COMPLEX16FrequencySeries *hctilde = NULL;
  REAL8 phic = 0;
  REAL8 deltaF = 0.06;
  REAL8 f_max = 0; // 8000;
  REAL8 distance = 100 * 1e6 * LAL_PC_SI;

  int ret = XLALSimIMRPhenomP(
    &hptilde,                 /**< Frequency-domain waveform h+ */
    &hctilde,                 /**< Frequency-domain waveform hx */
    chi1_l,                   /**< aligned spin on companion 1 */
    chi2_l,                   /**< aligned spin on companion 2 */
    chip,                     /**< Effective spin in the orbital plane */
    thetaJ,                   /**< Angle between J0 and line of sight (z-direction) */
    m1_SI,                    /**< mass of companion1 (kg) */
    m2_SI,                    /**< mass of companion1 (kg) */
    distance,                 /**< Distance of source (m) */
    alpha0,                   /**< Initial value of alpha angle (azimuthal precession angle) */
    phic,                     /**< Orbital coalescence phase (rad) */
    deltaF,                   /**< Sampling frequency (Hz) */
    f_min,                    /**< Starting GW frequency (Hz) */
    f_max,                    /**< End frequency; 0 defaults to ringdown cutoff freq */
    f_ref,                    /**< Reference frequency */
    version,
    NRTv,
    NULL);

  dump_file("PhenomP_Test_f_ref1.dat", hptilde, hctilde, m1+m2);
  MYUNUSED(ret);
  COMPLEX16 hp = (hptilde->data->data)[1000];
  COMPLEX16 hc = (hctilde->data->data)[1000];
  printf("f_ref = %g\n", f_ref);
  prC("hp", hp);
  prC("hc", hc);

  // Now repeat for a different f_ref
  f_ref = 5;

  XLALSimIMRPhenomPCalculateModelParametersOld(
      &chi1_l,            /**< Output: aligned spin on companion 1 */
      &chi2_l,            /**< Output: aligned spin on companion 2 */
      &chip,              /**< Output: Effective spin in the orbital plane */
      &thetaJ,            /**< Output: Angle between J0 and line of sight (z-direction) */
      &alpha0,            /**< Output: Initial value of alpha angle (azimuthal precession angle) */
      m1_SI,              /**< Mass of companion 1 (kg) */
      m2_SI,              /**< Mass of companion 2 (kg) */
      f_ref,              /**< Reference GW frequency (Hz) */
      lnhatx,             /**< Initial value of LNhatx: orbital angular momentum unit vector */
      lnhaty,             /**< Initial value of LNhaty */
      lnhatz,             /**< Initial value of LNhatz */
      s1x,                /**< Initial value of s1x: dimensionless spin of larger BH */
      s1y,                /**< Initial value of s1y: dimensionless spin of larger BH */
      s1z,                /**< Initial value of s1z: dimensionless spin of larger BH */
      s2x,                /**< Initial value of s2x: dimensionless spin of larger BH */
      s2y,                /**< Initial value of s2y: dimensionless spin of larger BH */
      s2z,                /**< Initial value of s2z: dimensionless spin of larger BH */
      version);

  COMPLEX16FrequencySeries *hptilde2 = NULL;
  COMPLEX16FrequencySeries *hctilde2 = NULL;

  ret = XLALSimIMRPhenomP(
    &hptilde2,                /**< Frequency-domain waveform h+ */
    &hctilde2,                /**< Frequency-domain waveform hx */
    chi1_l,                   /**< aligned spin on companion 1 */
    chi2_l,                   /**< aligned spin on companion 2 */
    chip,                     /**< Effective spin in the orbital plane */
    thetaJ,                   /**< Angle between J0 and line of sight (z-direction) */
    m1_SI,                    /**< mass of companion1 (kg) */
    m2_SI,                    /**< mass of companion1 (kg) */
    distance,                 /**< Distance of source (m) */
    alpha0,                   /**< Initial value of alpha angle (azimuthal precession angle) */
    phic,                     /**< Orbital coalescence phase (rad) */
    deltaF,                   /**< Sampling frequency (Hz) */
    f_min,                    /**< Starting GW frequency (Hz) */
    f_max,                    /**< End frequency; 0 defaults to ringdown cutoff freq */
    f_ref,                    /**< Reference frequency */
    version,
    NRTv,
    NULL);


  dump_file("PhenomP_Test_f_ref2.dat", hptilde2, hctilde2, m1+m2);
  MYUNUSED(ret);
  COMPLEX16 hp2 = (hptilde2->data->data)[1000];
  COMPLEX16 hc2 = (hctilde2->data->data)[1000];
  printf("f_ref = %g\n", f_ref);
  prC("hp", hp2);
  prC("hc", hc2);

  const REAL8 eps = 1e-5;

  XLAL_CHECK_EXIT(
       approximatelyEqualC(hp, hp2, eps)
    && approximatelyEqualC(hc, hc2, eps)
    && "XLALSimIMRPhenomP_f_ref()"
  );
}

int main(int argc, char *argv[]) {
  MYUNUSED(argc);
  MYUNUSED(argv);

#ifndef _OPENMP
  Test_alpha_epsilon();
  Test_XLALSimIMRPhenomPCalculateModelParameters();
  //Test_PhenomC();
  //Test_PhenomPCore();
  Test_XLALSimIMRPhenomP();
  //Test_PhenomC_PhenomP();
  Test_XLALSimIMRPhenomP_f_ref();
#else
  MYUNUSED(Test_alpha_epsilon);
  MYUNUSED(Test_XLALSimIMRPhenomPCalculateModelParameters);
  //MYUNUSED(Test_PhenomC);
  //MYUNUSED(Test_PhenomPCore);
  MYUNUSED(Test_XLALSimIMRPhenomP);
  //MYUNUSED(Test_PhenomC_PhenomP);
  MYUNUSED(Test_XLALSimIMRPhenomP_f_ref);
#endif

  printf("\nAll done!\n");
  return 0;
}
