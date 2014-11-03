/*
 * Copyright (C) 2011 J. Creighton
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

#ifndef _LALSIMBLACKHOLERINGDOWN_H
#define _LALSIMBLACKHOLERINGDOWN_H

#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/* LOW-LEVEL ROUTINES (USE LEAVER'S CONVENSIONS G = c = 2M = 1) */


/**
 * Low-level routine that computes the black hole quasinormal mode
 * eigenefrequency, omega, and angular separation constant A for a given
 * (l,m) mode and spin-weight s (s=-2 for gravitational perturbations).
 *
 * Implements Leaver's method by simultaneously
 * solving the continued fraction equations Eq. (21) and Eq. (27)
 * of Leaver (1985):
 * E. W. Leaver "An analyitic representation for the quasi-normal
 * modes of Kerr black holes", Proc. R. Soc. Lond. A 402 285-298 (1985).
 *
 * \warning The variables are represented in Leaver's conventions
 * in which G = c = 2M = 1.  In particular this means, |a| < 0.5.
 *
 * \todo Extend so that overtones can be computed too.
 */
int XLALSimBlackHoleRingdownModeEigenvaluesLeaver(
	COMPLEX16 *A,		/**< angular separation constant [returned] */
	COMPLEX16 *omega,		/**< eigenfrequency [returned] */
	double a,		/**< spin parameter (note: |a| < 0.5) */
	int l,			/**< mode value l */
	int m,			/**< mode value m */
	int s			/**< spin weight (s = -2 for gravitational perturbations) */
);


/**
 * Low-level routine that evaluates the spheroidal wave function at a
 * specified value of mu = cos(theta) for a given (l,m) mode and
 * spin-weight s (s=-2 for gravitational perturbations).
 * Also requires the angular separation constant A and eigenvalue
 * omega for that mode, which are calculated by the routine
 * XLALSimBlackHoleRingdownModeEigenvaluesLeaver().
 *
 * Implements Leaver's method by simultaneously
 * solving the continued fraction equations Eq. (21) and Eq. (27)
 * of Leaver (1985):
 * E. W. Leaver "An analyitic representation for the quasi-normal
 * modes of Kerr black holes", Proc. R. Soc. Lond. A 402 285-298 (1985).
 *
 * \warning The variables are represented in Leaver's conventions
 * in which G = c = 2M = 1.  In particular this means, |a| < 0.5.
 *
 * \todo Extend so that overtones can be computed too.
 */
COMPLEX16 XLALSimBlackHoleRingdownSpheroidalWaveFunctionLeaver(
	double mu,		/**< cosine of polar angle */
	double a,		/**< spin parameter (note: |a| < 0.5) */
	int l,			/**< mode value l */
	int m,			/**< mode value m */
	int s,			/**< spin weight (s = -2 for gravitational perturbations) */
	COMPLEX16 A,		/**< angular separation constant */
	COMPLEX16 omega		/**< eigenfrequency */
);


/* HIGH-LEVEL ROUTINES */


/**
 * Computes the frequency and quality factor of a specified quasinormal
 * mode (l,m) of spin weight s perturbations (s=-2 for gravitational
 * perturbations) of a black hole of a specified mass and spin.
 *
 * Uses the method of Leaver (1985):
 * E. W. Leaver "An analyitic representation for the quasi-normal
 * modes of Kerr black holes", Proc. R. Soc. Lond. A 402 285-298 (1985).
 *
 * \note The dimensionless spin assumes values between -1 and 1.
 *
 * \todo Extend so that overtones can be computed too.
 */
int XLALSimBlackHoleRingdownMode(
	double *frequency,		/**< mode frequency (Hz) [returned] */
	double *quality,		/**< mode quality factor [returned] */
	double mass,			/**< black hole mass (kg) */
	double dimensionless_spin,	/**< black hole dimensionless spin parameter (-1,+1) */
	int l,				/**< polar mode number */
	int m,				/**< azimuthal mode number */
	int s				/**< spin weight (s=-2 for gravitational radiation) */
);


/**
 * Evaluates the value of spheroidal wave function at a given
 * polar angle theta for a specified mode (l,m) and spin weight s
 * (s=-2 for gravitational perturbations) and
 * dimensionless spin parameter.
 *
 * Uses the method of Leaver (1985):
 * E. W. Leaver "An analyitic representation for the quasi-normal
 * modes of Kerr black holes", Proc. R. Soc. Lond. A 402 285-298 (1985).
 *
 * \note The dimensionless spin assumes values between -1 and 1.
 *
 * \todo Extend so that overtones can be computed too.
 */
COMPLEX16 XLALSimBlackHoleRingdownSpheroidalWaveFunction(
	double theta,			/**< polar angle (radians) */
	double dimensionless_spin,	/**< black hole dimensionless spin parameter */
	int l,				/**< polar mode number */
	int m,				/**< azimuthal mode number */
	int s				/**< spin weight (s=-2 for gravitational radiation) */
);


/**
 * Computes the waveform for the ringdown of a black hole
 * quasinormal mode (l,m).
 *
 * \note The dimensionless spin assumes values between -1 and 1.
 *
 * \todo Extend so that overtones can be computed too.
 */
int XLALSimBlackHoleRingdown(
	REAL8TimeSeries **hplus,	/**< plus-polarization waveform [returned] */
	REAL8TimeSeries **hcross,	/**< cross-polarization waveform [returned] */
	const LIGOTimeGPS *t0,		/**< start time of ringdown */
	double phi0,			/**< initial phase of ringdown (rad) */
	double deltaT,			/**< sampling interval (s) */
	double mass,			/**< black hole mass (kg) */
	double dimensionless_spin,	/**< black hole dimensionless spin parameter */
	double fractional_mass_loss,	/**< fraction of mass radiated in this mode */
	double distance,		/**< distance to source (m) */
	double inclination,		/**< inclination of source's spin axis (rad) */
	int l,				/**< polar mode number */
	int m				/**< azimuthal mode number */
);

/**
 * Computes the final mass and spin of the black hole resulting from merger.
 * They are given by fittings of NR simulations results. Specifically,
 * for EOBNR, Table I of Buonanno et al. PRD76, 104049;
 * for EOBNRv2 and EOBNRv2HM, Eqs. 29a and 29b of Pan et al. PRD84, 124052;
 * for SEOBNRv1, Eq. 8 of Tichy and Marronetti PRD78, 081501 and
 * Eqs. 1 and 3 of Barausse and Rezzolla ApJ704, L40.
 */
INT4 XLALSimIMREOBFinalMassSpin(
	REAL8    *finalMass,  /**<< OUTPUT, the final mass (scaled by original total mass) */
	REAL8    *finalSpin,  /**<< OUTPUT, the final spin (scaled by final mass) */
  const REAL8     mass1,      /**<< The mass of the 1st component of the system */
  const REAL8     mass2,      /**<< The mass of the 2nd component of the system */
  const REAL8     spin1[3],   /**<< The spin of the 1st object; only needed for spin waveforms */
  const REAL8     spin2[3],   /**<< The spin of the 2nd object; only needed for spin waveforms */
  Approximant     approximant /**<< The waveform approximant being used */
);

/**
 * This function generates the quasinormal mode frequencies for a black
 * hole ringdown. At present, this function works for the 22, 21, 33, 44
 * and 55 modes, and includes 8 overtones. The final frequencies are
 * computed by interpolating the data found on the webpage of
 * Vitor Cardoso, http://centra.ist.utl.pt/~vitor/?page=ringdown
 * In this page, frequecy data are given for positive final spins only.
 * For a negative final spin chi<0 case, the (l,m) mode frequency is given by
 * the (l,-m) mode frequency of the positive final spin -chi case.
 */
INT4 XLALSimIMREOBGenerateQNMFreqV2(
  COMPLEX16Vector *modefreqs, /**<< OUTPUT, complex freqs of overtones (scaled by total mass) */
  const REAL8      mass1,     /**<< The mass of the 1st component (in Solar masses) */
  const REAL8      mass2,     /**<< The mass of the 2nd component (in Solar masses) */
  const REAL8      spin1[3],  /**<< The spin of the 1st object; only needed for spin waveforms */
  const REAL8      spin2[3],  /**<< The spin of the 2nd object; only needed for spin waveforms */
  UINT4            l,         /**<< The l value of the mode in question */
  INT4             m,         /**<< The m value of the mode in question */
  UINT4            nmodes,    /**<< The number of overtones that should be included (max 8) */
  Approximant      approximant/**<< The waveform approximant being used */
  );

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMBLACKHOLERINGDOWN_H */
