/*
 * Copyright (C) 2010 Reinhard Prix
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

/*********************************************************************************/
/** \author R. Prix
 * \file
 * \brief
 * Generate samples of various statistics (F-stat, F-atoms, B-stat,...) drawn from their
 * respective distributions, assuming Gaussian noise, and drawing signal params from
 * their (given) priors
 *
 * This is based on synthesizeBstat, and is mostly meant to be used for efficient
 * Monte-Carlos studies, ROC curves etc
 *
 *********************************************************************************/

#ifndef _SYNTHESIZE_CW_DRAWS_H  /* Double-include protection. */
#define _SYNTHESIZE_CW_DRAWS_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/*---------- INCLUDES ----------*/

/* GSL includes */
#include <gsl/gsl_rng.h>

/* LAL includes */
#include <lal/SkyCoordinates.h>
#include <lal/LALComputeAM.h>

#include <lal/ProbabilityDensity.h>
#include <lal/TransientCW_utils.h>

/*---------- exported types ----------*/

/** Enumeration of allowed amplitude-prior types
 */
typedef enum {
  AMP_PRIOR_TYPE_PHYSICAL = 0,	/**< 'physical' priors: isotropic pdf{cosi,psi,phi0} AND flat pdf(h0) */
  AMP_PRIOR_TYPE_CANONICAL,	/**< 'canonical' priors: uniform in A^mu up to h_max */
  AMP_PRIOR_TYPE_LAST
} AmpPriorType_t;

/** Signal (amplitude) parameter ranges
 */
typedef struct tagAmplitudePrior_t {
  pdf1D_t *pdf_h0Nat;	/**< pdf for h0/sqrt{Sn} */
  REAL8 fixedSNR;	/**< alternative 1: adjust h0 to fix the optimal SNR of the signal */
  BOOLEAN fixRhohMax;	/**< alternative 2: draw h0 with fixed rhohMax = h0Max * (detM)^(1/8) <==> canonical Fstat prior */

  pdf1D_t *pdf_cosi;	/**< pdf(cosi) */
  pdf1D_t *pdf_psi;	/**< pdf(psi) */
  pdf1D_t *pdf_phi0;	/**< pdf(phi0) */
} AmplitudePrior_t;

/** struct for buffering of AM-coeffs, if signal for same sky-position is injected
 */
typedef struct tagmultiAMBuffer_t {
  SkyPosition skypos;		/**< sky-position for which we have AM-coeffs computed already */
  MultiAMCoeffs *multiAM;;	/**< pre-computed AM-coeffs for skypos */
} multiAMBuffer_t;

/** Hold all (generally) randomly drawn injection parameters: skypos, amplitude-params, M_mu_nu, transient-window, SNR
 */
typedef struct tagInjParams_t
{
  SkyPosition skypos;
  PulsarAmplitudeParams ampParams;
  PulsarAmplitudeVect ampVect;
  AntennaPatternMatrix M_mu_nu;
  transientWindow_t transientWindow;
  REAL8 SNR;
  REAL8 detM1o8;	// (detMp)^(1/8): rescale param between h0, and rhoh = h0 * (detMp)^(1/8)
} InjParams_t;


/*---------- Global variables ----------*/

/* empty struct initializers */
extern multiAMBuffer_t empty_multiAMBuffer;
extern InjParams_t empty_InjParams_t;

/*---------- exported prototypes [API] ----------*/
int XLALDrawCorrelatedNoise ( PulsarAmplitudeVect n_mu, const gsl_matrix *L, gsl_rng * rng );

// ----- API to synthesize F-stat atoms for transient CW searches
FstatAtomVector* XLALGenerateFstatAtomVector ( const DetectorStateSeries *detStates, const AMCoeffs *amcoeffs );
MultiFstatAtomVector*XLALGenerateMultiFstatAtomVector ( const MultiDetectorStateSeries *multiDetStates, const MultiAMCoeffs *multiAM );

int XLALAddNoiseToFstatAtomVector ( FstatAtomVector *atoms, gsl_rng * rng );
int XLALAddNoiseToMultiFstatAtomVector ( MultiFstatAtomVector *multiAtoms, gsl_rng * rng );

REAL8 XLALAddSignalToFstatAtomVector ( FstatAtomVector* atoms, AntennaPatternMatrix *M_mu_nu, const PulsarAmplitudeVect A_Mu, transientWindow_t transientWindow );
REAL8 XLALAddSignalToMultiFstatAtomVector ( MultiFstatAtomVector* multiAtoms, AntennaPatternMatrix *M_mu_nu, const PulsarAmplitudeVect A_Mu, transientWindow_t transientWindow, INT4 lineX );

int XLALRescaleMultiFstatAtomVector ( MultiFstatAtomVector* multiAtoms,	REAL8 rescale );
int write_InjParams_to_fp ( FILE * fp, const InjParams_t *par, UINT4 dataStartGPS );

MultiFstatAtomVector *
XLALSynthesizeTransientAtoms ( InjParams_t *injParamsOut,
                               SkyPosition skypos,
                               AmplitudePrior_t AmpPrior,
                               transientWindowRange_t transientInjectRange,
                               const MultiDetectorStateSeries *multiDetStates,
                               BOOLEAN SignalOnly,
                               multiAMBuffer_t *multiAMBuffer,
                               gsl_rng *rng,
                               INT4 lineX
                               );



#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
