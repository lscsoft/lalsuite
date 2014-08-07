/*
 *  Copyright (C) 2013,2014 Michael Puerrer, Alejandro Bohe
 *  Reuses code found in:
 *    - LALSimIMRPhenomC
 *    - LALSimInspiralSpinTaylorF2
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
 * \author Michael Puerrer, Alejandro Bohe
 *
 * \file
 *
 * \brief Functions for producing PhenomP waveforms for precessing binaries,
 * as described in Hannam et al., arXiv:1308.3271 [gr-qc].
 */

#include <stdlib.h>
#include <math.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <lal/Units.h>
#include <lal/XLALError.h>
#include <lal/SphericalHarmonics.h>
#include "LALSimIMR.h"
#include "LALSimIMRPhenomC_internals.c" /* This is ugly, but allows us to reuse internal PhenomC functions without making those functions XLAL */

/**************************** PhenomP internal function prototypes *****************************/

/* PhenomC parameters for modified ringdown: Uses final spin formula of Barausse & Rezzolla, Astrophys.J.Lett.704:L40-L44, 2009 */
BBHPhenomCParams *ComputeIMRPhenomCParamsRDmod(
  const REAL8 m1,   /**< Mass of companion 1 (solar masses) */
  const REAL8 m2,   /**< Mass of companion 2 (solar masses) */
  const REAL8 chi,  /**< Reduced aligned spin of the binary chi = (m1*chi1 + m2*chi2)/M */
  const REAL8 chip  /**< Dimensionless spin in the orbital plane */
);

typedef struct tagNNLOanglecoeffs {
    REAL8 alphacoeff1; /* Coefficient of omega^(-1)   in alphaNNLO */
    REAL8 alphacoeff2; /* Coefficient of omega^(-2/3) in alphaNNLO */
    REAL8 alphacoeff3; /* Coefficient of omega^(-1/3) in alphaNNLO */
    REAL8 alphacoeff4; /* Coefficient of log(omega)   in alphaNNLO */
    REAL8 alphacoeff5; /* Coefficient of omega^(1/3)  in alphaNNLO */

    REAL8 epsiloncoeff1; /* Coefficient of omega^(-1)   in epsilonNNLO */
    REAL8 epsiloncoeff2; /* Coefficient of omega^(-2/3) in epsilonNNLO */
    REAL8 epsiloncoeff3; /* Coefficient of omega^(-1/3) in epsilonNNLO */
    REAL8 epsiloncoeff4; /* Coefficient of log(omega)   in epsilonNNLO */
    REAL8 epsiloncoeff5; /* Coefficient of omega^(1/3)  in epsilonNNLO */
} NNLOanglecoeffs;

static void ComputeNNLOanglecoeffs(
  NNLOanglecoeffs *angcoeffs,  /**< Output: Structure to store results */
  const REAL8 q,               /**< Mass-ratio (convention q>1) */
  const REAL8 chil,            /**< Dimensionless aligned spin of the largest BH */
  const REAL8 chip             /**< Dimensionless spin component in the orbital plane */
);

typedef struct tagSpinWeightedSphericalHarmonic_l2 {
  COMPLEX16 Y2m2, Y2m1, Y20, Y21, Y22;
} SpinWeightedSphericalHarmonic_l2;

/* Internal core function to calculate PhenomP polarizations for a single frequency. */
int PhenomPCore(
  const REAL8 fHz,                        /**< Frequency (Hz) */
  const REAL8 eta,                        /**< Symmetric mass ratio */
  const REAL8 chi_eff,                    /**< Dimensionless effective total aligned spin */
  const REAL8 chip,                       /**< Dimensionless spin in the orbital plane */
  const REAL8 distance,                   /**< Distance of source (m) */
  const REAL8 M,                          /**< Total mass (Solar masses) */
  const REAL8 phic,                       /**< Orbital phase at the peak of the underlying non precessing model (rad) */
  BBHPhenomCParams *PCparams,             /**< Internal PhenomC parameters */
  NNLOanglecoeffs *angcoeffs,             /**< Struct with PN coeffs for the NNLO angles */
  SpinWeightedSphericalHarmonic_l2 *Y2m,  /**< Struct of l=2 spherical harmonics of spin weight -2 */
  const REAL8 alphaoffset,                /**< f_ref dependent offset for alpha angle */
  const REAL8 epsilonoffset,              /**< f_ref dependent offset for epsilon angle */
  COMPLEX16 *hp,                          /**< Output: \tilde h_+ */
  COMPLEX16 *hc                           /**< Output: \tilde h_+ */
);

/* Simple 2PN version of L, without any spin terms expressed as a function of v */
REAL8 L2PNR(
  const REAL8 v,   /**< Cubic root of (Pi * Frequency (geometric)) */
  const REAL8 eta  /**< Symmetric mass-ratio */
);

void WignerdCoefficients(
  const REAL8 v,          /**< Cubic root of (Pi * Frequency (geometric)) */
  const REAL8 SL,         /**< Dimensionfull aligned spin */
  const REAL8 eta,        /**< Symmetric mass-ratio */
  const REAL8 Sp,         /**< Dimensionfull spin component in the orbital plane */
  REAL8 *cos_beta_half,   /**< Output: cos(beta/2) */
  REAL8 *sin_beta_half    /**< Output: sin(beta/2) */
);

REAL8 FinalSpinBarausse2009_all_spin_on_larger_BH(
  const REAL8 nu,     /**< Symmetric mass-ratio */
  const REAL8 chi,    /**< Reduced aligned spin of the binary chi = (m1*chi1 + m2*chi2)/M */
  const REAL8 chip    /**< Dimensionless spin in the orbital plane */
);

REAL8 FinalSpinBarausse2009_aligned_spin_equally_distributed(
  const REAL8 nu,     /**< Symmetric mass-ratio */
  const REAL8 chi,    /**< Reduced aligned spin of the binary chi = (m1*chi1 + m2*chi2)/M */
  const REAL8 chip    /**< Dimensionless spin in the orbital plane */
);

REAL8 FinalSpinBarausse2009(  /* Barausse & Rezzolla, Astrophys.J.Lett.704:L40-L44, 2009, arXiv:0904.2577 */
  const REAL8 nu,               /**< Symmetric mass-ratio */
  const REAL8 a1,               /**< |a_1| norm of dimensionless spin vector for BH 1 */
  const REAL8 a2,               /**< |a_2| norm of dimensionless spin vector for BH 2 */
  const REAL8 cos_alpha,        /**< cos(alpha) = \hat a_1 . \hat a_2 (Eq. 7) */
  const REAL8 cos_beta_tilde,   /**< cos(\tilde beta)  = \hat a_1 . \hat L (Eq. 9) */
  const REAL8 cos_gamma_tilde   /**< cos(\tilde gamma) = \hat a_2 . \hat L (Eq. 9)*/
);

static size_t NextPow2(const size_t n); /* Return the closest higher power of 2. */


/*************************************** Implementation ****************************************/

int XLALSimIMRPhenomPCalculateModelParameters(
    REAL8 *chi_eff,                 /**< Output: Effective aligned spin */
    REAL8 *chip,                    /**< Output: Effective spin in the orbital plane */
    REAL8 *eta,                     /**< Output: Symmetric mass-ratio */
    REAL8 *thetaJ,                  /**< Output: Angle between J0 and line of sight (z-direction) */
    REAL8 *phiJ,                    /**< Output: Angle of J0 in the plane of the sky */
    REAL8 *alpha0,                  /**< Output: Initial value of alpha angle */
    const REAL8 m1_SI,              /**< Mass of companion 1 (kg) */
    const REAL8 m2_SI,              /**< Mass of companion 2 (kg) */
    const REAL8 f_ref,              /**< Reference GW frequency (Hz) */
    const REAL8 lnhatx,             /**< Initial value of LNhatx: orbital angular momentum unit vector */
    const REAL8 lnhaty,             /**< Initial value of LNhaty */
    const REAL8 lnhatz,             /**< Initial value of LNhatz */
    const REAL8 s1x,                /**< Initial value of s1x: dimensionless spin of BH 1 */
    const REAL8 s1y,                /**< Initial value of s1y: dimensionless spin of BH 1 */
    const REAL8 s1z,                /**< Initial value of s1z: dimensionless spin of BH 1 */
    const REAL8 s2x,                /**< Initial value of s2x: dimensionless spin of BH 2 */
    const REAL8 s2y,                /**< Initial value of s2y: dimensionless spin of BH 2 */
    const REAL8 s2z)                /**< Initial value of s2z: dimensionless spin of BH 2 */
{
  /* Check inputs for sanity */
  if (!chi_eff)  XLAL_ERROR(XLAL_EFAULT);
  if (!chip)     XLAL_ERROR(XLAL_EFAULT);
  if (!eta)      XLAL_ERROR(XLAL_EFAULT);
  if (!thetaJ)   XLAL_ERROR(XLAL_EFAULT);
  if (!phiJ)     XLAL_ERROR(XLAL_EFAULT);
  if (!alpha0)   XLAL_ERROR(XLAL_EFAULT);

  if (f_ref <= 0) {
    XLALPrintError("Reference frequency must be positive.\n");
    XLAL_ERROR(XLAL_EDOM);
  }

  const REAL8 m1 = m1_SI / LAL_MSUN_SI;   /* Masses in solar masses */
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;
  const REAL8 M = m1+m2;
  const REAL8 m1_2 = m1*m1;
  const REAL8 m2_2 = m2*m2;
  *eta = m1 * m2 / (M*M);    /* Symmetric mass-ratio */

  /* Aligned spins */
  const REAL8 chi1_l = lnhatx*s1x + lnhaty*s1y + lnhatz*s1z; /* Dimensionless aligned spin on BH 1 */
  const REAL8 chi2_l = lnhatx*s2x + lnhaty*s2y + lnhatz*s2z; /* Dimensionless aligned spin on BH 2 */

  /* Spin components orthogonal to lnhat */
  const REAL8 S1_perp_x = (s1x - chi1_l*lnhatx) * m1_2;
  const REAL8 S1_perp_y = (s1y - chi1_l*lnhaty) * m1_2;
  const REAL8 S1_perp_z = (s1z - chi1_l*lnhatz) * m1_2;
  const REAL8 S2_perp_x = (s2x - chi2_l*lnhatx) * m2_2;
  const REAL8 S2_perp_y = (s2y - chi2_l*lnhaty) * m2_2;
  const REAL8 S2_perp_z = (s2z - chi2_l*lnhatz) * m2_2;
  const REAL8 S1_perp = sqrt(S1_perp_x*S1_perp_x + S1_perp_y*S1_perp_y + S1_perp_z*S1_perp_z);
  const REAL8 S2_perp = sqrt(S2_perp_x*S2_perp_x + S2_perp_y*S2_perp_y + S2_perp_z*S2_perp_z);

  *chi_eff = (m1*chi1_l + m2*chi2_l) / M; /* Effective aligned spin */

  const REAL8 A1 = 2 + (3*m2) / (2*m1);
  const REAL8 A2 = 2 + (3*m1) / (2*m2);
  const REAL8 ASp1 = A1*S1_perp;
  const REAL8 ASp2 = A2*S2_perp;
  const REAL8 num = (ASp2 > ASp1) ? ASp2 : ASp1;
  const REAL8 den = (m2 > m1) ? A2*m2_2 : A1*m1_2;
  *chip = num / den; /*  chip = max(A1 Sp1, A2 Sp2) / (A_i m_i^2) for i index of larger BH */

  /* Compute L, J0 and orientation angles */
  const REAL8 m_sec = M * LAL_MTSUN_SI;   /* Total mass in seconds */
  const REAL8 piM = LAL_PI * m_sec;
  const REAL8 v_ref = cbrt(piM * f_ref);

  const REAL8 L0 = M*M * L2PNR(v_ref, *eta); /* Use 2PN approximation for L. */
  const REAL8 Jx0 = L0 * lnhatx + m1_2*s1x + m2_2*s2x;
  const REAL8 Jy0 = L0 * lnhaty + m1_2*s1y + m2_2*s2y;
  const REAL8 Jz0 = L0 * lnhatz + m1_2*s1z + m2_2*s2z;
  const REAL8 J0 = sqrt(Jx0*Jx0 + Jy0*Jy0 + Jz0*Jz0);

  *thetaJ = acos(Jz0 / J0); /* Angle between J0 and line of sight (z-direction) */
  if (Jx0 == 0.0 && Jy0 == 0.0)
    *phiJ = 0;
  else
    *phiJ = atan2(Jy0, Jx0); /* Angle of J0 in the plane of the sky */
    /* Note: Compared to the similar code in SpinTaylorF2 we have defined phiJ as the angle between the positive
    (rather than the negative) x-axis and the projection of J0, since this is a more natural definition of the angle.
    We have also renamed the angle from psiJ to phiJ. */

  /* Rotate Lnhat back to frame where J is along z, to figure out initial alpha */
  /* Note: We have flipped the sign of the sin(phiJ) terms to rotate by -phiJ in agreement with our definition of phiJ above. */
  const REAL8 rotLx = lnhatx*cos(*thetaJ)*cos(*phiJ) + lnhaty*cos(*thetaJ)*sin(*phiJ) + lnhatz*sin(*thetaJ);
  const REAL8 rotLy = -lnhatx*sin(*phiJ) + lnhaty*cos(*phiJ);
  /* AB: this rotation sends the line of sight to the xz plane so that, if we want to interpret the
     value of alpha0 that we calculate below as an initial orientation in this new frame we have to
     compute our Ylms at yphi=0. In addition, in XLALSimIMRPhenomP I use the fact that the waveform
     only depends on yphi-alpha0 and shortcut the alpha0 dependence by putting it in the Ylm. */
  if (rotLx == 0.0 && rotLy == 0.0)
    *alpha0 = 0.0;
  else
    *alpha0 = atan2(rotLy, rotLx);

  return XLAL_SUCCESS;
}

int XLALSimIMRPhenomP(
  COMPLEX16FrequencySeries **hptilde,   /**< Output: Frequency-domain waveform h+ */
  COMPLEX16FrequencySeries **hctilde,   /**< Output: Frequency-domain waveform hx */
  const REAL8 chi_eff,                  /**< Effective aligned spin */
  const REAL8 chip,                     /**< Effective spin in the orbital plane */
  const REAL8 eta,                      /**< Symmetric mass-ratio */
  const REAL8 thetaJ,                   /**< Angle between J0 and line of sight (z-direction) */
  const REAL8 phiJ,                     /**< Angle of J0 in the plane of the sky */
  const REAL8 Mtot_SI,                  /**< Total mass of binary (kg) */
  const REAL8 distance,                 /**< Distance of source (m) */
  const REAL8 alpha0,                   /**< Initial value of alpha angle */
  const REAL8 phic,                     /**< Orbital phase at the peak of the underlying non precessing model (rad) */
  const REAL8 deltaF,                   /**< Sampling frequency (Hz) */
  const REAL8 f_min,                    /**< Starting GW frequency (Hz) */
  const REAL8 f_max,                    /**< End frequency; 0 defaults to ringdown cutoff freq */
  const REAL8 f_ref)                    /**< Reference frequency */
{
  const REAL8 M = Mtot_SI / LAL_MSUN_SI;  /* External units: SI; internal units: solar masses */
  const REAL8 m_sec = M * LAL_MTSUN_SI;   /* Total mass in seconds */
  const REAL8 q = (1.0 + sqrt(1.0 - 4.0*eta) - 2.0*eta)/(2.0*eta); /* Mass-ratio */
  const REAL8 m1 = M * 1.0 / (1+q);
  const REAL8 m2 = M * q / (1+q);
  const REAL8 piM = LAL_PI * m_sec;
  const REAL8 v0 = cbrt(piM * f_ref);

  static LIGOTimeGPS ligotimegps_zero = {0, 0};

  /* Check inputs for sanity */
  if (*hptilde)       XLAL_ERROR(XLAL_EFAULT);
  if (*hctilde)       XLAL_ERROR(XLAL_EFAULT);
  if (deltaF <= 0)    XLAL_ERROR(XLAL_EDOM);
  if (M <= 0)         XLAL_ERROR(XLAL_EDOM);
  if (f_min <= 0)     XLAL_ERROR(XLAL_EDOM);
  if (f_max < 0)      XLAL_ERROR(XLAL_EDOM);
  if (distance <= 0)  XLAL_ERROR(XLAL_EDOM);

  if (f_ref <= 0) {
      XLALPrintError("Reference frequency must be positive.\n");
      XLAL_ERROR(XLAL_EDOM);
  }

  if (eta < 0.0453515){ /* q = 20 */
      XLALPrintError("Mass ratio is way outside the calibration range. m1/m2 should be <= 20.\n");
      XLAL_ERROR(XLAL_EDOM);
  }
  else if (eta < 0.16) { /* q = 4 */
      XLALPrintWarning("Warning: The model is only calibrated for m1/m2 <= 4.\n");
  }

  if (fabs(chi_eff) > 1) XLAL_ERROR(XLAL_EDOM);
  if (fabs(chip) > 1) XLAL_ERROR(XLAL_EDOM);

  /* If spins are above 0.9 or below -0.9, throw an error. */
  /* The rationale behind this is given at this page: https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/WaveformsReview/IMRPhenomCdevel-SanityCheck01 */
  if (chi_eff > 0.9 || chi_eff < -0.9){
      XLALPrintError("Spins outside the range [-0.9,0.9] are not supported\n");
      XLAL_ERROR(XLAL_EDOM);
  }

  const REAL8 chil = (1.0+q)/q * chi_eff; /* dimensionless aligned spin of the largest BH */
  NNLOanglecoeffs angcoeffs;
  ComputeNNLOanglecoeffs(&angcoeffs,q,chil,chip);

  /* Compute the offsets due to the choice of integration constant in alpha and epsilon PN formula */
  const REAL8 omega_ref = piM * f_ref;
  const REAL8 logomega_ref = log(omega_ref);
  const REAL8 omega_ref_cbrt = v0;
  const REAL8 omega_ref_cbrt2 = omega_ref_cbrt*omega_ref_cbrt;
  const REAL8 alphaNNLOoffset = (angcoeffs.alphacoeff1/omega_ref
                              + angcoeffs.alphacoeff2/omega_ref_cbrt2
                              + angcoeffs.alphacoeff3/omega_ref_cbrt
                              + angcoeffs.alphacoeff4*logomega_ref
                              + angcoeffs.alphacoeff5*omega_ref_cbrt);

  const REAL8 epsilonNNLOoffset = (angcoeffs.epsiloncoeff1/omega_ref
                                + angcoeffs.epsiloncoeff2/omega_ref_cbrt2
                                + angcoeffs.epsiloncoeff3/omega_ref_cbrt
                                + angcoeffs.epsiloncoeff4*logomega_ref
                                + angcoeffs.epsiloncoeff5*omega_ref_cbrt);

  /* Compute Ylm's only once and pass them to PhenomPCore() below. */
  SpinWeightedSphericalHarmonic_l2 Y2m;
  const REAL8 ytheta  = thetaJ;
  const REAL8 yphi    = 0;
  Y2m.Y2m2 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2, -2);
  Y2m.Y2m1 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2, -1);
  Y2m.Y20  = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2,  0);
  Y2m.Y21  = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2,  1);
  Y2m.Y22  = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2,  2);


  /* Phenomenological parameters */
  BBHPhenomCParams *PCparams;
  //PCparams = ComputeIMRPhenomCParams(m1, m2, chi_eff); /* original PhenomC */
  PCparams = ComputeIMRPhenomCParamsRDmod(m1, m2, chi_eff, chip); /* PhenomC with ringdown using Barausse 2009 formula for final spin */
  if (!PCparams) XLAL_ERROR(XLAL_EFUNC);
  if (PCparams->fCut <= f_min) {
      XLALPrintError("(fCut = 0.15M) <= f_min\n");
      XLAL_ERROR(XLAL_EDOM);
  }

  /* Default f_max to params->fCut */
  REAL8 f_max_prime = f_max ? f_max : PCparams->fCut;
  f_max_prime = (f_max_prime > PCparams->fCut) ? PCparams->fCut : f_max_prime;
  if (f_max_prime <= f_min) {
      XLALPrintError("f_max <= f_min\n");
      XLAL_ERROR(XLAL_EDOM);
  }
  XLALPrintInfo("f_max_prime = %g\tfcut = %g\tv = %g\n", f_max_prime, PCparams->fCut, cbrt(piM * f_max_prime));


  /* Allocate hp, hc */
  XLALPrintInfo("f_max / deltaF = %g\n", f_max_prime / deltaF);
  size_t n = NextPow2(f_max_prime / deltaF) + 1; /* Note: Should explain why the length is one plus a power of 2. */
  XLALPrintInfo("n = %d\n", (int)(n));
  *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &ligotimegps_zero, 0.0, deltaF, &lalStrainUnit, n);
  *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &ligotimegps_zero, 0.0, deltaF, &lalStrainUnit, n);
  memset((*hptilde)->data->data, 0, n * sizeof(COMPLEX16));
  memset((*hctilde)->data->data, 0, n * sizeof(COMPLEX16));
  XLALUnitMultiply(&((*hptilde)->sampleUnits), &((*hptilde)->sampleUnits), &lalSecondUnit);
  XLALUnitMultiply(&((*hctilde)->sampleUnits), &((*hctilde)->sampleUnits), &lalSecondUnit);
  if (!(*hptilde) || !(*hctilde))
    XLAL_ERROR(XLAL_EFUNC);

  /* Test output */
  XLALPrintInfo("eta: %g\n", eta);
  XLALPrintInfo("m1: %g\n", m1);
  XLALPrintInfo("m2: %g\n", m2);
  XLALPrintInfo("chi: %g\n", chi_eff);
  XLALPrintInfo("chip: %g\n", chip);
  XLALPrintInfo("thetaJ: %g\n", thetaJ);
  XLALPrintInfo("phiJ: %g\n", phiJ);
  XLALPrintInfo("ind_mix = (int)(f_mix / deltaF) = %d\n", (int)(f_min / deltaF));
  XLALPrintInfo("ind_max = (int)(f_max / deltaF) = %d\n", (int)(f_max_prime / deltaF));

  /* Note: there will usually be zero data at the beginning and end of the frequency series  */
  size_t i_min = (size_t) (f_min / deltaF);
  size_t i_max = (size_t) (f_max_prime / deltaF);
  int errcode = XLAL_SUCCESS;
  /*
    We can't call XLAL_ERROR() directly with OpenMP on.
    Keep track of return codes for each thread and in addition use flush to get out of
    the parallel for loop as soon as possible if something went wrong in any thread.
  */
  #pragma omp parallel for
  for (size_t i = i_min; i < i_max; i++) {
    COMPLEX16 hp_val, hc_val;
    REAL8 f = i * deltaF;
    int per_thread_errcode;

    #pragma omp flush(errcode)
    if (errcode != XLAL_SUCCESS)
      goto skip;

    /* Generate the waveform */
    per_thread_errcode = PhenomPCore(f, eta, chi_eff, chip, distance, M, phic,
                              PCparams, &angcoeffs, &Y2m,
                              alphaNNLOoffset - alpha0, epsilonNNLOoffset,
                              &hp_val, &hc_val);

    if (per_thread_errcode != XLAL_SUCCESS) {
      errcode = per_thread_errcode;
      #pragma omp flush(errcode)
    }


    ((*hptilde)->data->data)[i] = hp_val;
    ((*hctilde)->data->data)[i] = hc_val;

    skip: /* this statement intentionally left blank */;
  }

  if( errcode != XLAL_SUCCESS )
    XLAL_ERROR(errcode);
  else
    return XLAL_SUCCESS;
}

/******************************* PhenomP internal functions *********************************/

/* Internal core function to calculate PhenomP polarizations for a single frequency. */
int PhenomPCore(
  const REAL8 fHz,                        /**< Frequency (Hz) */
  const REAL8 eta,                        /**< Symmetric mass ratio */
  const REAL8 chi_eff,                    /**< Dimensionless effective total aligned spin */
  const REAL8 chip,                       /**< Dimensionless spin in the orbital plane */
  const REAL8 distance,                   /**< Distance of source (m) */
  const REAL8 M,                          /**< Total mass (Solar masses) */
  const REAL8 phic,                       /**< Orbital phase at the peak of the underlying non precessing model (rad) */
  BBHPhenomCParams *PCparams,             /**< Internal PhenomC parameters */
  NNLOanglecoeffs *angcoeffs,  		  /**< Struct with PN coeffs for the NNLO angles */
  SpinWeightedSphericalHarmonic_l2 *Y2m,  /**< Struct of l=2 spherical harmonics of spin weight -2 */
  const REAL8 alphaoffset,                /**< f_ref dependent offset for alpha angle */
  const REAL8 epsilonoffset,              /**< f_ref dependent offset for epsilon angle */
  COMPLEX16 *hp,                          /**< Output: \tilde h_+ */
  COMPLEX16 *hc)                          /**< Output: \tilde h_+ */
{

  /* Calculate PhenomC amplitude and phase for a given frequency. */
  REAL8 phPhenomC = 0.0;
  REAL8 aPhenomC  = 0.0;
  int errcode = IMRPhenomCGenerateAmpPhase( &aPhenomC, &phPhenomC, fHz, eta, PCparams );
  if( errcode != XLAL_SUCCESS )
    XLAL_ERROR(XLAL_EFUNC);
  phPhenomC -= 2.*phic; /* Note: phic is orbital phase */
  REAL8 amp0 = M * LAL_MRSUN_SI * M * LAL_MTSUN_SI / distance;
  COMPLEX16 hPC = amp0 * aPhenomC * cexp(-I*phPhenomC); /* Assemble PhenomC waveform. */

  REAL8 f = fHz*LAL_MTSUN_SI*M; /* Frequency in geometric units */

  /*
   * We put all spin on the larger BH. Here m2 >= m1.
   * chi_eff = (m1*0 + m2*chi2)/M
   * chil = chi2 = chi_eff / m2 (for M=1)
   * SL = m2^2 chi2 = m2*M*chi = q/(1+q) * chi (for M=1)
   * Sp = chip*m2^2
   */

  const REAL8 q = (1.0 + sqrt(1.0 - 4.0*eta) - 2.0*eta)/(2.0*eta);
  const REAL8 m2 = q/(1.0+q);         /* Compute the mass of the larger BH for unit total mass M=1. */
  const REAL8 SL = chi_eff*m2;        /* Dimensionfull aligned spin of the largest BH. */
  const REAL8 Sperp = chip*(m2*m2);   /* Dimensionfull spin component in the orbital plane. S_perp = S_2_perp */


  /*
   * Calculate plus and cross polarizations of the PhenomP model for individual m's.
   * The general expression for the modes h^P_{2m}(t) is given by Eq. 1 of arXiv:1308.3271.
   * We calculate the frequency domain l=2 plus and cross polarizations separately for each m = -2, ... , 2.
   * The expression of the polarizations times the Ylm's in code notation are:
   *    \tilde (h_2m)_+ = e^{-2i epsilon} (e^{-i m alpha} d^2_{-2,m} (-2Y_2m) + e^{+i m alpha} d^2_{2,m} (-2Y_2m)^*) * hPC / 2
   *    \tilde (h_2m)_x = e^{-2i epsilon} (e^{-i m alpha} d^2_{-2,m} (-2Y_2m) - e^{+i m alpha} d^2_{2,m} (-2Y_2m)^*) * hPC / 2
   * where the d^l_{m',m} are Wigner d-matrices evaluated at - beta, and hPC is the PhenomC frequency domain model hPC(f) = PCAmp(f) e^{-i PCPhase(f)}.
   *
   */

  /* Compute PN NNLO angles */
  const REAL8 omega = LAL_PI * f;
  const REAL8 logomega = log(omega);
  const REAL8 omega_cbrt = cbrt(omega);
  const REAL8 omega_cbrt2 = omega_cbrt*omega_cbrt;

  REAL8 alpha = (angcoeffs->alphacoeff1/omega
              + angcoeffs->alphacoeff2/omega_cbrt2
              + angcoeffs->alphacoeff3/omega_cbrt
              + angcoeffs->alphacoeff4*logomega
              + angcoeffs->alphacoeff5*omega_cbrt) - alphaoffset;

  REAL8 epsilon = (angcoeffs->epsiloncoeff1/omega
                + angcoeffs->epsiloncoeff2/omega_cbrt2
                + angcoeffs->epsiloncoeff3/omega_cbrt
                + angcoeffs->epsiloncoeff4*logomega
                + angcoeffs->epsiloncoeff5*omega_cbrt) - epsilonoffset;

  /* Calculate intermediate expressions cos(beta/2), sin(beta/2) and powers thereof for Wigner d's. */
  REAL8 cBetah, sBetah; /* cos(beta/2), sin(beta/2) */
  WignerdCoefficients(omega_cbrt, SL, eta, Sperp, &cBetah, &sBetah);
  const REAL8 cBetah2 = cBetah*cBetah;
  const REAL8 cBetah3 = cBetah2*cBetah;
  const REAL8 cBetah4 = cBetah3*cBetah;
  const REAL8 sBetah2 = sBetah*sBetah;
  const REAL8 sBetah3 = sBetah2*sBetah;
  const REAL8 sBetah4 = sBetah3*sBetah;

  /* Compute Wigner d coefficients
    The expressions below agree with refX and Mathematica
    d2  = Table[WignerD[{2, mp, 2}, 0, -\[Beta], 0], {mp, -2, 2}]
    dm2 = Table[WignerD[{2, mp, -2}, 0, -\[Beta], 0], {mp, -2, 2}]
  */
  COMPLEX16 d2[5]   = {sBetah4, 2*cBetah*sBetah3, sqrt(6.0)*sBetah2*cBetah2, 2*cBetah3*sBetah, cBetah4};
  COMPLEX16 dm2[5]  = {d2[4], -d2[3], d2[2], -d2[1], d2[0]}; /* Exploit symmetry d^2_{-2,-m} = (-1)^m d^2_{2,m} */

  COMPLEX16 Y2mA[5] = {Y2m->Y2m2, Y2m->Y2m1, Y2m->Y20, Y2m->Y21, Y2m->Y22};
  COMPLEX16 hp_sum = 0;
  COMPLEX16 hc_sum = 0;

  /* Sum up contributions to \tilde h+ and \tilde hx */
  /* Precompute powers of e^{i m alpha} */
  COMPLEX16 cexp_i_alpha = cexp(+I*alpha);
  COMPLEX16 cexp_2i_alpha = cexp_i_alpha*cexp_i_alpha;
  COMPLEX16 cexp_mi_alpha = 1.0/cexp_i_alpha;
  COMPLEX16 cexp_m2i_alpha = cexp_mi_alpha*cexp_mi_alpha;
  COMPLEX16 cexp_im_alpha[5] = {cexp_m2i_alpha, cexp_mi_alpha, 1.0, cexp_i_alpha, cexp_2i_alpha};
  for(int m=-2; m<=2; m++) {
    COMPLEX16 T2m   = cexp_im_alpha[-m+2] * dm2[m+2] *      Y2mA[m+2];  /*  = cexp(-I*m*alpha) * dm2[m+2] *      Y2mA[m+2] */
    COMPLEX16 Tm2m  = cexp_im_alpha[m+2]  * d2[m+2]  * conj(Y2mA[m+2]); /*  = cexp(+I*m*alpha) * d2[m+2]  * conj(Y2mA[m+2]) */
    hp_sum +=     T2m + Tm2m;
    hc_sum += +I*(T2m - Tm2m);
  }

  COMPLEX16 eps_phase_hPC = cexp(-2*I*epsilon) * hPC / 2.0;
  *hp = eps_phase_hPC * hp_sum;
  *hc = eps_phase_hPC * hc_sum;

  return XLAL_SUCCESS;
}


void ComputeNNLOanglecoeffs(
  NNLOanglecoeffs *angcoeffs, /**< Output: Structure to store results */
  const REAL8 q,              /**< Mass-ratio (convention q>1) */
  const REAL8 chil,           /**< Dimensionless aligned spin of the largest BH */
  const REAL8 chip)           /**< Dimensionless spin component in the orbital plane */
{
  const REAL8 m2 = q/(1. + q);
  const REAL8 m1 = 1./(1. + q);
  const REAL8 dm = m1 - m2;
  const REAL8 mtot = 1.;
  const REAL8 eta = m1*m2; /* mtot = 1 */
  const REAL8 eta2 = eta*eta;
  const REAL8 eta3 = eta2*eta;
  const REAL8 eta4 = eta3*eta;
  const REAL8 mtot2 = mtot*mtot;
  const REAL8 mtot4 = mtot2*mtot2;
  const REAL8 mtot6 = mtot4*mtot2;
  const REAL8 mtot8 = mtot6*mtot2;
  const REAL8 chil2 = chil*chil;
  const REAL8 chip2 = chip*chip;
  const REAL8 chip4 = chip2*chip2;
  const REAL8 dm2 = dm*dm;
  const REAL8 dm3 = dm2*dm;
  const REAL8 m2_2 = m2*m2;
  const REAL8 m2_3 = m2_2*m2;
  const REAL8 m2_4 = m2_3*m2;
  const REAL8 m2_5 = m2_4*m2;
  const REAL8 m2_6 = m2_5*m2;
  const REAL8 m2_7 = m2_6*m2;
  const REAL8 m2_8 = m2_7*m2;


  angcoeffs->alphacoeff1 = (-0.18229166666666666 - (5*dm)/(64.*m2));

  angcoeffs->alphacoeff2 = ((-15*dm*m2*chil)/(128.*mtot2*eta) - (35*m2_2*chil)/(128.*mtot2*eta));

  angcoeffs->alphacoeff3 = (-1.7952473958333333 - (4555*dm)/(7168.*m2) -
        (15*chip2*dm*m2_3)/(128.*mtot4*eta2) -
        (35*chip2*m2_4)/(128.*mtot4*eta2) - (515*eta)/384. - (15*dm2*eta)/(256.*m2_2) -
        (175*dm*eta)/(256.*m2));

  angcoeffs->alphacoeff4 = - (35*LAL_PI)/48. - (5*dm*LAL_PI)/(16.*m2) +
     (5*dm2*chil)/(16.*mtot2) + (5*dm*m2*chil)/(3.*mtot2) +
     (2545*m2_2*chil)/(1152.*mtot2) -
     (5*chip2*dm*m2_5*chil)/(128.*mtot6*eta3) -
     (35*chip2*m2_6*chil)/(384.*mtot6*eta3) + (2035*dm*m2*chil)/(21504.*mtot2*eta) +
     (2995*m2_2*chil)/(9216.*mtot2*eta);

  angcoeffs->alphacoeff5 = (4.318908476114694 + (27895885*dm)/(2.1676032e7*m2) -
        (15*chip4*dm*m2_7)/(512.*mtot8*eta4) -
        (35*chip4*m2_8)/(512.*mtot8*eta4) -
        (485*chip2*dm*m2_3)/(14336.*mtot4*eta2) +
        (475*chip2*m2_4)/(6144.*mtot4*eta2) +
        (15*chip2*dm2*m2_2)/(256.*mtot4*eta) + (145*chip2*dm*m2_3)/(512.*mtot4*eta) +
        (575*chip2*m2_4)/(1536.*mtot4*eta) + (39695*eta)/86016. + (1615*dm2*eta)/(28672.*m2_2) -
        (265*dm*eta)/(14336.*m2) + (955*eta2)/576. + (15*dm3*eta2)/(1024.*m2_3) +
        (35*dm2*eta2)/(256.*m2_2) + (2725*dm*eta2)/(3072.*m2) - (15*dm*m2*LAL_PI*chil)/(16.*mtot2*eta) -
        (35*m2_2*LAL_PI*chil)/(16.*mtot2*eta) + (15*chip2*dm*m2_7*chil2)/(128.*mtot8*eta4) +
        (35*chip2*m2_8*chil2)/(128.*mtot8*eta4) +
        (375*dm2*m2_2*chil2)/(256.*mtot4*eta) + (1815*dm*m2_3*chil2)/(256.*mtot4*eta) +
        (1645*m2_4*chil2)/(192.*mtot4*eta));

  angcoeffs->epsiloncoeff1 = (-0.18229166666666666 - (5*dm)/(64.*m2));

  angcoeffs->epsiloncoeff2 = ((-15*dm*m2*chil)/(128.*mtot2*eta) - (35*m2_2*chil)/(128.*mtot2*eta));

  angcoeffs->epsiloncoeff3 = (-1.7952473958333333 - (4555*dm)/(7168.*m2) - (515*eta)/384. -
        (15*dm2*eta)/(256.*m2_2) - (175*dm*eta)/(256.*m2));

  angcoeffs->epsiloncoeff4 = - (35*LAL_PI)/48. - (5*dm*LAL_PI)/(16.*m2) +
     (5*dm2*chil)/(16.*mtot2) + (5*dm*m2*chil)/(3.*mtot2) +
     (2545*m2_2*chil)/(1152.*mtot2) + (2035*dm*m2*chil)/(21504.*mtot2*eta) +
     (2995*m2_2*chil)/(9216.*mtot2*eta);

  angcoeffs->epsiloncoeff5 = (4.318908476114694 + (27895885*dm)/(2.1676032e7*m2) + (39695*eta)/86016. +
        (1615*dm2*eta)/(28672.*m2_2) - (265*dm*eta)/(14336.*m2) + (955*eta2)/576. +
        (15*dm3*eta2)/(1024.*m2_3) + (35*dm2*eta2)/(256.*m2_2) +
        (2725*dm*eta2)/(3072.*m2) - (15*dm*m2*LAL_PI*chil)/(16.*mtot2*eta) - (35*m2_2*LAL_PI*chil)/(16.*mtot2*eta) +
        (375*dm2*m2_2*chil2)/(256.*mtot4*eta) + (1815*dm*m2_3*chil2)/(256.*mtot4*eta) +
        (1645*m2_4*chil2)/(192.*mtot4*eta));
}


REAL8 L2PNR(
  const REAL8 v,   /**< Cubic root of (Pi * Frequency (geometric)) */
  const REAL8 eta) /**< Symmetric mass-ratio */
{
  const REAL8 mu = eta; /* M=1 */
  const REAL8 v2 = v*v;
  const REAL8 v3 = v2*v;
  const REAL8 v4 = v3*v;
  const REAL8 eta2 = eta*eta;

  /* Simple 2PN version of L, without any spin terms expressed as a function of v:
  [Kidder, Phys. Rev. D 52, 821–847 (1995), Eq. 2.9] */

  return mu*sqrt((1 - ((3 - eta)*v2)/3. + (4.75 + eta/9.)*eta*v4)/v2)*
    (1 + ((1 - 3*eta)*v2)/2. + (3*(1 - 7*eta + 13*eta2)*v4)/8. +
      ((14 - 41*eta + 4*eta2)*v4)/(4.*pow(1 - ((3 - eta)*v2)/3. + (4.75 + eta/9.)*eta*v4,2)) +
      ((3 + eta)*v2)/(1 - ((3 - eta)*v2)/3. + (4.75 + eta/9.)*eta*v4) +
      ((7 - 10*eta - 9*eta2)*v4)/(2.*(1 - ((3 - eta)*v2)/3. + (4.75 + eta/9.)*eta*v4)));
}

void WignerdCoefficients(
  const REAL8 v,        /**< Cubic root of (Pi * Frequency (geometric)) */
  const REAL8 SL,       /**< Dimensionfull aligned spin */
  const REAL8 eta,      /**< Symmetric mass-ratio */
  const REAL8 Sp,       /**< Dimensionfull spin component in the orbital plane */
  REAL8 *cos_beta_half, /**< Output: cos(beta/2) */
  REAL8 *sin_beta_half) /**< Output: sin(beta/2) */
{
  /* cos(beta) = \hat J . \hat L = (1 + (Sp / (L + SL)) )^(-1/2) */
  /* We compute
      cos(beta/2) = ((1 + cos(beta))/2)^(1/2)
      sin(beta/2) = ((1 - cos(beta))/2)^(1/2)
    by Taylor expanding in the parameter s (see below) up to third order.
  */
  REAL8 s = Sp / (L2PNR(v, eta) + SL);  /* s := Sp / (L + SL) */
  REAL8 s2 = s*s;
  REAL8 s3 = s2*s;
  *cos_beta_half = 1.0 - s2/8.0;            /* cos(beta/2) */
  *sin_beta_half = s/2.0 - (3.0*s3)/16.0;   /* sin(beta/2) */
}

REAL8 FinalSpinBarausse2009_all_spin_on_larger_BH(
  const REAL8 nu,     /**< Symmetric mass-ratio */
  const REAL8 chi,    /**< Reduced aligned spin of the binary chi = (m1*chi1 + m2*chi2)/M */
  const REAL8 chip)   /**< Dimensionless spin in the orbital plane */
{
  /* Use convention m1>m2 as in arXiv:0904.2577 */
  /* Put all spin on larger BH: a1 = (chip, 0, chi), a2 = (0,0,0), L = (0,0,1) */
  const REAL8 a1_x = chip;
  const REAL8 a1_y = 0;
  const REAL8 a1_z = chi;
  const REAL8 a2_x = 0;
  const REAL8 a2_y = 0;
  const REAL8 a2_z = 0;

  const REAL8 a1 = sqrt(a1_x*a1_x + a1_y*a1_y + a1_z*a1_z);
  const REAL8 a2 = sqrt(a2_x*a2_x + a2_y*a2_y + a2_z*a2_z);

  const REAL8 cos_alpha = (a1*a2 == 0) ? 0.0 : a1_z*a2_z/(a1*a2); /* cos(alpha) = \hat a1 . \hat a2 (Eq. 7) */
  const REAL8 cos_beta_tilde  = (a1 == 0) ? 0.0 : a1_z/a1;  /* \cos(\tilde \beta)  = \hat a1 . \hat L  (Eq. 9) */
  const REAL8 cos_gamma_tilde = (a2 == 0) ? 0.0 : a2_z/a2;  /* \cos(\tilde \gamma) = \hat a2 . \hat L (Eq. 9) */

  return FinalSpinBarausse2009(nu, a1, a2, cos_alpha, cos_beta_tilde, cos_gamma_tilde);
}

REAL8 FinalSpinBarausse2009_aligned_spin_equally_distributed(
  const REAL8 nu,     /**< Symmetric mass-ratio */
  const REAL8 chi,    /**< Reduced aligned spin of the binary chi = (m1*chi1 + m2*chi2)/M */
  const REAL8 chip)   /**< Dimensionless spin in the orbital plane */
{
  /* Use convention m1>m2 as in arXiv:0904.2577 */
  /* Put all spin on larger BH: a1 = (chip, 0, chi), a2 = (0,0,0), L = (0,0,1) */
  const REAL8 a1_x = chip;
  const REAL8 a1_y = 0;
  const REAL8 a1_z = chi;
  const REAL8 a2_x = 0;
  const REAL8 a2_y = 0;
  const REAL8 a2_z = chi;

  const REAL8 a1 = sqrt(a1_x*a1_x + a1_y*a1_y + a1_z*a1_z);
  const REAL8 a2 = sqrt(a2_x*a2_x + a2_y*a2_y + a2_z*a2_z);

  const REAL8 cos_alpha = (a1*a2 == 0) ? 0.0 : a1_z*a2_z/(a1*a2); /* cos(alpha) = \hat a1 . \hat a2 (Eq. 7) */
  const REAL8 cos_beta_tilde  = (a1 == 0) ? 0.0 : a1_z/a1;  /* \cos(\tilde \beta)  = \hat a1 . \hat L  (Eq. 9) */
  const REAL8 cos_gamma_tilde = (a2 == 0) ? 0.0 : a2_z/a2;  /* \cos(\tilde \gamma) = \hat a2 . \hat L (Eq. 9) */

  return FinalSpinBarausse2009(nu, a1, a2, cos_alpha, cos_beta_tilde, cos_gamma_tilde);
}


REAL8 FinalSpinBarausse2009(  /* Barausse & Rezzolla, Astrophys.J.Lett.704:L40-L44, 2009, arXiv:0904.2577 */
  const REAL8 nu,               /**< Symmetric mass-ratio */
  const REAL8 a1,               /**< |a_1| norm of dimensionless spin vector for BH 1 */
  const REAL8 a2,               /**< |a_2| norm of dimensionless spin vector for BH 2 */
  const REAL8 cos_alpha,        /**< cos(alpha) = \hat a_1 . \hat a_2 (Eq. 7) */
  const REAL8 cos_beta_tilde,   /**< cos(\tilde beta)  = \hat a_1 . \hat L (Eq. 9) */
  const REAL8 cos_gamma_tilde)  /**< cos(\tilde gamma) = \hat a_2 . \hat L (Eq. 9)*/
{
  REAL8 q = (2*nu)/(1 + sqrt(1 - 4*nu) - 2*nu); /* Use convention q = m2/m1 <= 1 as in arXiv:0904.2577 */

  /* These parameters are defined in eq. 3. */
  const REAL8 s4 = -0.1229;
  const REAL8 s5 = 0.4537;
  const REAL8 t0 = -2.8904;
  const REAL8 t2 = -3.5171;
  const REAL8 t3 = 2.5763;

  /* shorthands */
  const REAL8 nu2 = nu*nu;
  const REAL8 q2 = q*q;
  const REAL8 q4 = q2*q2;
  const REAL8 q2p = 1 + q2;
  const REAL8 q2p2 = q2p*q2p;
  const REAL8 qp = 1 + q;
  const REAL8 qp2 = qp*qp;
  const REAL8 a1_2 = a1*a1;
  const REAL8 a2_2 = a2*a2;

  /* l = \tilde l/(m1*m2), where \tilde l = S_fin - (S1 + S2) = L - J_rad (Eq. 4) */
  const REAL8 l = 2*sqrt(3.0) + t2*nu + t3*nu2
                + (s4 / q2p2) * (a1_2 + a2_2*q4 + 2*a1*a2*q2*cos_alpha)
                + ((s5*nu + t0 + 2)/q2p) * (a1*cos_beta_tilde + a2*cos_gamma_tilde*q2); /* |l| (Eq. 10) */
  const REAL8 l2 = l*l;

  /* a_fin = S_fin/M^2  (Eq. 6) */
  const REAL8 a_fin = (1.0 / qp2) * sqrt(a1_2 + a2_2*q4 + 2*a1*a2*q2*cos_alpha + 2*(a1*cos_beta_tilde + a2*q2*cos_gamma_tilde)*l*q + l2*q2);
  return a_fin;
}


/* PhenomC parameters for modified ringdown: Uses final spin formula of Barausse & Rezzolla, Astrophys.J.Lett.704:L40-L44, 2009 */
BBHPhenomCParams *ComputeIMRPhenomCParamsRDmod(
  const REAL8 m1,   /**< Mass of companion 1 (solar masses) */
  const REAL8 m2,   /**< Mass of companion 2 (solar masses) */
  const REAL8 chi,  /**< Reduced aligned spin of the binary chi = (m1*chi1 + m2*chi2)/M */
  const REAL8 chip) /**< Dimensionless spin in the orbital plane */
{

  BBHPhenomCParams *p = NULL;
  p = ComputeIMRPhenomCParams(m1, m2, chi); /* populate parameters with the original PhenomC setup */
  if( !p )
    XLAL_ERROR_NULL(XLAL_EFUNC);

  const REAL8 M = m1 + m2;
  const REAL8 eta = m1 * m2 / (M * M);

  REAL8 finspin = FinalSpinBarausse2009_all_spin_on_larger_BH(eta, chi, chip);
  if( fabs(finspin) > 1.0 )
    XLAL_ERROR_NULL( XLAL_EDOM );

  p->afin = finspin;

  /* Get the Ringdown frequency */
  REAL8 prefac = (1./(2.*LAL_PI)) * LAL_C_SI * LAL_C_SI * LAL_C_SI / (LAL_G_SI * M * LAL_MSUN_SI);
  REAL8 k1 = 1.5251;
  REAL8 k2 = -1.1568;
  REAL8 k3 = 0.1292;

  p->fRingDown = (prefac * (k1 + k2 * pow(1. - fabs(finspin), k3)));
  p->MfRingDown = p->m_sec * p->fRingDown;

  /* Get the quality factor of ring-down, using Eq (5.6) of Main paper (arxiv.org/pdf/1005.3306v3.pdf) */
  p->Qual = (0.7000 + (1.4187 * pow(1.0 - fabs(finspin), -0.4990)) );

  /* Get the transition frequencies, at which the model switches phase and
   * amplitude prescriptions, as used in Eq.(5.9), (5.13) of the Main paper */
  p->f1 = 0.1 * p->fRingDown;
  p->f2 = p->fRingDown;
  p->f0 = 0.98 * p->fRingDown;
  p->d1 = 0.005;
  p->d2 = 0.005;
  p->d0 = 0.015;

  /* Get the coefficients beta1, beta2, defined in Eq 5.7 of the main paper */
  REAL8 Mfrd = p->MfRingDown;

  p->b2 = ((-5./3.)* p->a1 * pow(Mfrd,(-8./3.)) - p->a2/(Mfrd*Mfrd) -
      (p->a3/3.)*pow(Mfrd,(-4./3.)) + (2./3.)* p->a5 * pow(Mfrd,(-1./3.)) + p->a6)/eta;

  REAL8 psiPMrd = IMRPhenomCGeneratePhasePM( p->fRingDown, eta, p );

  p->b1 = psiPMrd - (p->b2 * Mfrd);

  /* Taking the upper cut-off frequency as 0.15M */
  /*p->fCut = (1.7086 * eta * eta - 0.26592 * eta + 0.28236) / p->piM;*/
  p->fCut = 0.15 / p->m_sec;

  return p;
}
