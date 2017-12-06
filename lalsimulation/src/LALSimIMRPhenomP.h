#ifndef _LALSIM_IMR_PHENOMP_H
#define _LALSIM_IMR_PHENOMP_H

/*
 *  Copyright (C) 2013,2014,2015 Michael Puerrer, Alejandro Bohe
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

#include <lal/LALStdlib.h>
#include <lal/LALSimIMR.h>
#include <lal/LALConstants.h>

#include "LALSimIMRPhenomC_internals.h"
#include "LALSimIMRPhenomD_internals.h"

#include <lal/FrequencySeries.h>
#include <lal/LALSimInspiral.h>


/* ************************** PhenomP internal function prototypes *****************************/


/* PhenomC parameters for modified ringdown: Uses final spin formula of Barausse & Rezzolla, Astrophys.J.Lett.704:L40-L44, 2009 */
static BBHPhenomCParams *ComputeIMRPhenomCParamsRDmod(
  const REAL8 m1,   /**< Mass of companion 1 (solar masses) */
  const REAL8 m2,   /**< Mass of companion 2 (solar masses) */
  const REAL8 chi,  /**< Reduced aligned spin of the binary chi = (m1*chi1 + m2*chi2)/M */
  const REAL8 chip,  /**< Dimensionless spin in the orbital plane */
  const LALSimInspiralTestGRParam *extraParams /**< linked list containing the extra testing GR parameters */
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

/* Internal core function to calculate PhenomP polarizations for a sequence of frequences. */
static int PhenomPCore(
  COMPLEX16FrequencySeries **hptilde,   /**< Output: Frequency-domain waveform h+ */
  COMPLEX16FrequencySeries **hctilde,   /**< Output: Frequency-domain waveform hx */
  const REAL8 chi1_l_in,                /**< Dimensionless aligned spin on companion 1 */
  const REAL8 chi2_l_in,                /**< Dimensionless aligned spin on companion 2 */
  const REAL8 chip,                     /**< Effective spin in the orbital plane */
  const REAL8 thetaJ,                   /**< Angle between J0 and line of sight (z-direction) */
  const REAL8 m1_SI_in,                 /**< Mass of companion 1 (kg) */
  const REAL8 m2_SI_in,                 /**< Mass of companion 2 (kg) */
  const REAL8 distance,                 /**< Distance of source (m) */
  const REAL8 alpha0,                   /**< Initial value of alpha angle (azimuthal precession angle) */
  const REAL8 phic,                     /**< Orbital phase at the peak of the underlying non precessing model (rad) */
  const REAL8 f_ref,                    /**< Reference frequency */
  const REAL8Sequence *freqs,           /**< Frequency points at which to evaluate the waveform (Hz) */
  double deltaF,                        /**< Sampling frequency (Hz).
   * If deltaF > 0, the frequency points given in freqs are uniformly spaced with
   * spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
   * Then we will use deltaF = 0 to create the frequency series we return. */
  IMRPhenomP_version_type IMRPhenomP_version, /**< IMRPhenomPv1 uses IMRPhenomC, IMRPhenomPv2 uses IMRPhenomD */
  const LALSimInspiralTestGRParam *extraParams /**< linked list containing the extra testing GR parameters */
);

/* Internal core function to calculate PhenomP polarizations for a single frequency. */
static int PhenomPCoreOneFrequency(
  const REAL8 fHz,                        /**< Frequency (Hz) */
  const REAL8 eta,                        /**< Symmetric mass ratio */
  const REAL8 chi1_l,                     /**< Dimensionless aligned spin on companion 1 */
  const REAL8 chi2_l,                     /**< Dimensionless aligned spin on companion 2 */
  const REAL8 chip,                       /**< Dimensionless spin in the orbital plane */
  const REAL8 distance,                   /**< Distance of source (m) */
  const REAL8 M,                          /**< Total mass (Solar masses) */
  const REAL8 phic,                       /**< Orbital phase at the peak of the underlying non precessing model (rad) */
  IMRPhenomDAmplitudeCoefficients *pAmp,  /**< Internal IMRPhenomD amplitude coefficients */
  IMRPhenomDPhaseCoefficients *pPhi,      /**< Internal IMRPhenomD phase coefficients */
  BBHPhenomCParams *PCparams,             /**< Internal PhenomC parameters */
  PNPhasingSeries *PNparams,              /**< PN inspiral phase coefficients */
  NNLOanglecoeffs *angcoeffs,             /**< Struct with PN coeffs for the NNLO angles */
  SpinWeightedSphericalHarmonic_l2 *Y2m,  /**< Struct of l=2 spherical harmonics of spin weight -2 */
  const REAL8 alphaoffset,                /**< f_ref dependent offset for alpha angle (azimuthal precession angle) */
  const REAL8 epsilonoffset,              /**< f_ref dependent offset for epsilon angle */
  COMPLEX16 *hp,                          /**< Output: tilde h_+ */
  COMPLEX16 *hc,                          /**< Output: tilde h_+ */
  REAL8 *phasing,                         /**< Output: overall phasing */
  const UINT4 IMRPhenomP_version,         /**< Version number: 1 uses IMRPhenomC, 2 uses IMRPhenomD */
  AmpInsPrefactors *amp_prefactors,       /**< pre-calculated (cached for saving runtime) coefficients for amplitude. See LALSimIMRPhenomD_internals.c*/
  PhiInsPrefactors *phi_prefactors        /**< pre-calculated (cached for saving runtime) coefficients for phase. See LALSimIMRPhenomD_internals.*/
);

/* Simple 2PN version of L, without any spin terms expressed as a function of v */
static REAL8 L2PNR(
  const REAL8 v,   /**< Cubic root of (Pi * Frequency (geometric)) */
  const REAL8 eta  /**< Symmetric mass-ratio */
);

static REAL8 L2PNR_v1(
  const REAL8 v,   /**< Cubic root of (Pi * Frequency (geometric)) */
  const REAL8 eta  /**< Symmetric mass-ratio */
);

static void WignerdCoefficients(
  REAL8 *cos_beta_half,   /**< Output: cos(beta/2) */
  REAL8 *sin_beta_half,   /**< Output: sin(beta/2) */
  const REAL8 v,          /**< Cubic root of (Pi * Frequency (geometric)) */
  const REAL8 SL,         /**< Dimensionfull aligned spin */
  const REAL8 eta,        /**< Symmetric mass-ratio */
  const REAL8 Sp          /**< Dimensionfull spin component in the orbital plane */
);

static void WignerdCoefficients_SmallAngleApproximation(
  REAL8 *cos_beta_half, /**< Output: cos(beta/2) */
  REAL8 *sin_beta_half, /**< Output: sin(beta/2) */
  const REAL8 v,        /**< Cubic root of (Pi * Frequency (geometric)) */
  const REAL8 SL,       /**< Dimensionfull aligned spin */
  const REAL8 eta,      /**< Symmetric mass-ratio */
  const REAL8 Sp        /**< Dimensionfull spin component in the orbital plane */
);

static void CheckMaxOpeningAngle(
  const REAL8 m1,     /**< Mass of companion 1 (solar masses) */
  const REAL8 m2,     /**< Mass of companion 2 (solar masses) */
  const REAL8 chi1_l, /**< Aligned spin of BH 1 */
  const REAL8 chi2_l, /**< Aligned spin of BH 2 */
  const REAL8 chip    /**< Dimensionless spin in the orbital plane */
);

static REAL8 FinalSpinIMRPhenomD_all_in_plane_spin_on_larger_BH(
  const REAL8 m1,     /**< Mass of companion 1 (solar masses) */
  const REAL8 m2,     /**< Mass of companion 2 (solar masses) */
  const REAL8 chi1_l, /**< Aligned spin of BH 1 */
  const REAL8 chi2_l, /**< Aligned spin of BH 2 */
  const REAL8 chip    /**< Dimensionless spin in the orbital plane */
);

static REAL8 FinalSpinBarausse2009_all_spin_on_larger_BH(
  const REAL8 nu,     /**< Symmetric mass-ratio */
  const REAL8 chi,    /**< Effective aligned spin of the binary:  chi = (m1*chi1 + m2*chi2)/M  */
  const REAL8 chip    /**< Dimensionless spin in the orbital plane */
);

static REAL8 FinalSpinBarausse2009(  /* Barausse & Rezzolla, Astrophys.J.Lett.704:L40-L44, 2009, arXiv:0904.2577 */
  const REAL8 nu,               /**< Symmetric mass-ratio */
  const REAL8 a1,               /**< |a_1| norm of dimensionless spin vector for BH 1 */
  const REAL8 a2,               /**< |a_2| norm of dimensionless spin vector for BH 2 */
  const REAL8 cos_alpha,        /**< cos(alpha) = \\hat a_1 . \\hat a_2 (Eq. 7) */
  const REAL8 cos_beta_tilde,   /**< cos(\\tilde beta)  = \\hat a_1 . \\hat L (Eq. 9) */
  const REAL8 cos_gamma_tilde   /**< cos(\\tilde gamma) = \\hat a_2 . \\hat L (Eq. 9)*/
);

static bool approximately_equal(REAL8 x, REAL8 y, REAL8 epsilon);
static void nudge(REAL8 *x, REAL8 X, REAL8 epsilon);

#endif	// of #ifndef _LALSIM_IMR_PHENOMP_H
