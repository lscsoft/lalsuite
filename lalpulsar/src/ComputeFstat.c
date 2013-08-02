//
// Copyright (C) 2012, 2013 Karl Wette
// Copyright (C) 2007 Chris Messenger
// Copyright (C) 2006 John T. Whelan, Badri Krishnan
// Copyright (C) 2005, 2006, 2007, 2010 Reinhard Prix
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
//

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

#include <lal/ComputeFstat.h>
#include <lal/ExtrapolatePulsarSpins.h>

#define MYSIGN(x) ( ((x) < 0) ? (-1.0):(+1.0) )
#define SQ(x) ( (x) * (x) )

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

static const LALStatus empty_status;

///// Internal function prototypes /////

// Performs some common input sanity checks, and allocate a FstatInputData* struct
static FstatInputData*
SetupFstat_Common(
  MultiSFTVector **multiSFTs,
  MultiNoiseWeights **multiWeights,
  const EphemerisData *edat,
  const SSBprecision SSBprec
  );

///// Internal struct definitions /////

// Common input data for F-statistic algorithms
typedef struct {
  UINT4 numDetectors;					// Number of detectors
  CHAR detectorNames[PULSAR_MAX_DETECTORS][3];		// Names of detectors
  REAL8 Fnorm;						// F-statistic normalisation factor
  MultiNoiseWeights *multiWeights;			// Multi-detector noise weights
} FstatInputData_Common;

// Input data specific to F-statistic algorithms
typedef struct tagFstatInputData_Demod FstatInputData_Demod;
typedef struct tagFstatInputData_Resamp FstatInputData_Resamp;

// Internal definition of input data structure
struct tagFstatInputData {
  FstatInputData_Common common;				// Common input data
  FstatInputData_Demod* demod;				// Demodulation input data
  FstatInputData_Resamp* resamp;			// Resampling input data
};

///// Include F-statistic algorithm implementations /////

#include "ComputeFstat_Demod.c"
#include "ComputeFstat_Resamp.c"

///// Function definitions /////

static FstatInputData*
SetupFstat_Common(
  MultiSFTVector **multiSFTs,
  MultiNoiseWeights **multiWeights,
  const EphemerisData *edat,
  const SSBprecision SSBprec
  )
{

  // Check input
  XLAL_CHECK_NULL(multiSFTs != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(*multiSFTs != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(multiWeights != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(edat != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(SSBprec < SSBPREC_LAST, XLAL_EINVAL);

  // Check number of SFTs
  XLAL_CHECK_NULL((*multiSFTs)->length > 0, XLAL_EINVAL, "Found no SFTs!");
  XLAL_CHECK_NULL((*multiSFTs)->length <= PULSAR_MAX_DETECTORS, XLAL_EINVAL, "Supports only up to PULSAR_MAX_DETECTORS=%u detectors", PULSAR_MAX_DETECTORS);
  for (UINT4 X = 0; X < (*multiSFTs)->length; ++X) {

    // Check number of SFTs for each detector
    XLAL_CHECK_NULL((*multiSFTs)->data[X]->length > 0, XLAL_EINVAL, "Found no SFTs from detector %u", X);
    for (UINT4 alpha = 0; alpha < (*multiSFTs)->data[X]->length; ++alpha) {

      // Check length of SFTs
      XLAL_CHECK_NULL((*multiSFTs)->data[X]->data[alpha].data->length > 0, XLAL_EINVAL,
                      "Found zero-length SFT from detector %u, position %u", X, alpha);

    }

  }

  // If noise weights were supplied ...
  if (*multiWeights != NULL) {

    // Check numbers of noise weights match SFTs
    XLAL_CHECK_NULL((*multiWeights)->length == (*multiSFTs)->length, XLAL_EINVAL,
                    "Number of noise weight detectors does not match SFTS: %u != %u",
                    (*multiWeights)->length, (*multiSFTs)->length);
    for (UINT4 X = 0; X < (*multiSFTs)->length; ++X) {

      // Check number of noise weights for each detector
      XLAL_CHECK_NULL((*multiWeights)->data[X]->length == (*multiSFTs)->data[X]->length, XLAL_EINVAL,
                      "Number of noise weights from detector %u does not match SFTS: %u != %u",
                      X, (*multiWeights)->data[X]->length, (*multiSFTs)->data[X]->length);

    }

  }

  // Allocate input data struct
  FstatInputData* input = XLALCalloc(1, sizeof(*input));
  XLAL_CHECK_NULL(input != NULL, XLAL_ENOMEM);

  // Save number of detectors, and copy name of each detector
  input->common.numDetectors = (*multiSFTs)->length;
  for (UINT4 X = 0; X < (*multiSFTs)->length; ++X) {
    strncpy(input->common.detectorNames[X], (*multiSFTs)->data[X]->data[0].name, 2);
  }

  // If no noise weights were supplied ...
  if (*multiWeights == NULL) {

    // Correction to F-statistic quantities computed by XLALComputeFstat() without noise weights
    const REAL8 Tsft = 1.0 / (*multiSFTs)->data[0]->data[0].deltaF;
    input->common.Fnorm = 1.0 / sqrt( 0.5 * Tsft );

  } else {
    input->common.Fnorm = 0.0;
  }

  // Save pointer to input noise weights, set supplied pointer to NULL
  input->common.multiWeights = *multiWeights;
  *multiWeights = NULL;

  return input;

}

FstatInputDataVector*
XLALCreateFstatInputDataVector(
  const UINT4 length
  )
{

  // Allocate and initialise vector container
  FstatInputDataVector* inputs = XLALCalloc(1, sizeof(*inputs));
  XLAL_CHECK_NULL(inputs != NULL, XLAL_ENOMEM);
  inputs->length = length;

  // Allocate and initialise vector data
  if (inputs->length > 0) {
    inputs->data = XLALCalloc(inputs->length, sizeof(inputs->data[0]));
    XLAL_CHECK_NULL(inputs->data != NULL, XLAL_ENOMEM);
  }

  return inputs;

}

void
XLALDestroyFstatInputDataVector(
  FstatInputDataVector* inputs
  )
{
  if (inputs != NULL) {
    for (UINT4 i = 0; i < inputs->length; ++i) {
      XLALDestroyFstatInputData(inputs->data[i]);
    }
    XLALFree(inputs->data);
    XLALFree(inputs);
  }
}

int
XLALComputeFstat(
  FstatResults **Fstats,
  FstatInputData *input,
  const PulsarDopplerParams *doppler,
  const REAL8 dFreq,
  const UINT4 numFreqBins,
  const FstatQuantities whatToCompute
  )
{

  // Check input
  XLAL_CHECK(Fstats != NULL, XLAL_EFAULT);
  XLAL_CHECK(input != NULL, XLAL_EFAULT);
  XLAL_CHECK(doppler != NULL, XLAL_EFAULT);
  XLAL_CHECK(doppler->orbit == NULL, XLAL_EINVAL, "Binary parameters are currently not supported!");
  XLAL_CHECK(dFreq > 0 || (numFreqBins == 1 && dFreq >= 0), XLAL_EINVAL);
  XLAL_CHECK(numFreqBins > 0, XLAL_EINVAL);
  XLAL_CHECK(0 < whatToCompute && whatToCompute < FSTATQ_LAST, XLAL_EINVAL);

  // Allocate results struct, if needed
  if (*Fstats == NULL) {
    *Fstats = XLALCalloc(1, sizeof(**Fstats));
    XLAL_CHECK(*Fstats != NULL, XLAL_ENOMEM);
  }

  // Get constant pointer to common input data
  const FstatInputData_Common *common = &(input->common);
  const UINT4 numDetectors = common->numDetectors;

  // Enlarge result arrays if they are too small
  const bool moreFreqBins = (numFreqBins > (*Fstats)->internalalloclen);
  const bool moreDetectors = (numDetectors > (*Fstats)->numDetectors);
  if (moreFreqBins || moreDetectors) {

    // Enlarge multi-detector 2F array
    if ((whatToCompute & FSTATQ_2F) && moreFreqBins) {
      (*Fstats)->twoF = XLALRealloc((*Fstats)->twoF, numFreqBins*sizeof((*Fstats)->twoF[0]));
      XLAL_CHECK((*Fstats)->twoF != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->twoF to length %u", numFreqBins);
    }

    // Enlarge multi-detector Fa & Fb array
    if ((whatToCompute & FSTATQ_FAFB) && moreFreqBins) {
      (*Fstats)->FaFb = XLALRealloc((*Fstats)->FaFb, numFreqBins*sizeof((*Fstats)->FaFb[0]));
      XLAL_CHECK((*Fstats)->FaFb != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->FaFb to length %u", numFreqBins);
    }

    // Enlarge 2F per detector arrays
    if ((whatToCompute & FSTATQ_2F_PER_DET) && (moreFreqBins || moreDetectors)) {
      for (UINT4 X = 0; X < numDetectors; ++X) {
        (*Fstats)->twoFPerDet[X] = XLALRealloc((*Fstats)->twoFPerDet[X], numFreqBins*sizeof((*Fstats)->twoFPerDet[X][0]));
        XLAL_CHECK((*Fstats)->twoFPerDet[X] != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->twoFPerDet[%u] to length %u", X, numFreqBins);
      }
    }

    // Enlarge Fa & Fb per detector arrays
    if ((whatToCompute & FSTATQ_FAFB_PER_DET) && (moreFreqBins || moreDetectors)) {
      for (UINT4 X = 0; X < numDetectors; ++X) {
        (*Fstats)->FaFbPerDet[X] = XLALRealloc((*Fstats)->FaFbPerDet[X], numFreqBins*sizeof((*Fstats)->FaFbPerDet[X][0]));
        XLAL_CHECK((*Fstats)->FaFbPerDet[X] != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->FaFbPerDet[%u] to length %u", X, numFreqBins);
      }
    }

    // Enlarge F-atoms per detector arrays, and initialise to NULL
    if ((whatToCompute & FSTATQ_ATOMS_PER_DET) && (moreFreqBins || moreDetectors)) {
      for (UINT4 X = 0; X < numDetectors; ++X) {
        (*Fstats)->multiFatoms = XLALRealloc((*Fstats)->multiFatoms, numFreqBins*sizeof((*Fstats)->multiFatoms[0]));
        XLAL_CHECK((*Fstats)->multiFatoms != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->multiFatoms to length %u", numFreqBins);

        // If more detectors are needed, destroy multi-F-atom vectors so they can be re-allocated later
        if (moreDetectors) {
          for (UINT4 k = 0; k < numFreqBins; ++k) {
            XLALDestroyMultiFstatAtomVector((*Fstats)->multiFatoms[k]);
            (*Fstats)->multiFatoms[k] = NULL;
          }
        } else {
          for (UINT4 k = (*Fstats)->internalalloclen; k < numFreqBins; ++k) {
            (*Fstats)->multiFatoms[k] = NULL;
          }
        }

      }
    }

    // Update allocated length of arrays
    (*Fstats)->internalalloclen = numFreqBins;

  } // if (moreFreqBins || moreDetectors)

  // Initialise result struct parameters
  (*Fstats)->doppler = *doppler;
  (*Fstats)->dFreq = dFreq;
  (*Fstats)->numFreqBins = numFreqBins;
  (*Fstats)->numDetectors = numDetectors;
  memcpy((*Fstats)->detectorNames, common->detectorNames, sizeof(common->detectorNames));
  (*Fstats)->whatWasComputed = whatToCompute;

  // Call the appropriate algorithm function to compute the F-statistic
  if (input->demod != NULL) {
    XLAL_CHECK(ComputeFstat_Demod(*Fstats, common, input->demod) == XLAL_SUCCESS, XLAL_EFUNC);
  } else if (input->resamp != NULL) {
    XLAL_CHECK(ComputeFstat_Resamp(*Fstats, common, input->resamp) == XLAL_SUCCESS, XLAL_EFUNC);
  } else {
    XLAL_ERROR(XLAL_EFAILED, "Invalid FstatInputData struct passed to %s()", __func__);
  }

  // Correct F-statistic quantities when no noise weights are given
  if (common->Fnorm != 0.0) {
    const REAL8 Fnorm = common->Fnorm;
    const REAL8 Fnorm_sqr = Fnorm * Fnorm;

    // Correct antenna pattern matrix
    (*Fstats)->Mmunu.Sinv_Tsft = 2.0 / Fnorm_sqr;   // equivalent to Tsft

    // Correct multi-detector 2F array
    if (whatToCompute & FSTATQ_2F) {
      for (UINT4 k = 0; k < (*Fstats)->numFreqBins; ++k) {
        (*Fstats)->twoF[k] *= Fnorm_sqr;
        (*Fstats)->twoF[k] += 4;
      }
    }

    // Correct multi-detector F-parts array
    if (whatToCompute & FSTATQ_FAFB) {
      for (UINT4 k = 0; k < (*Fstats)->numFreqBins; ++k) {
        (*Fstats)->FaFb[k].Fa *= Fnorm;
        (*Fstats)->FaFb[k].Fb *= Fnorm;
      }
    }

    // Correct 2F per detector arrays
    if (whatToCompute & FSTATQ_2F_PER_DET) {
      for (UINT4 X = 0; X < numDetectors; ++X) {
        for (UINT4 k = 0; k < (*Fstats)->numFreqBins; ++k) {
          (*Fstats)->twoFPerDet[X][k] *= Fnorm_sqr;
          (*Fstats)->twoFPerDet[X][k] += 4;
        }
      }
    }

    // Correct F-parts per detector arrays
    if (whatToCompute & FSTATQ_FAFB_PER_DET) {
      for (UINT4 X = 0; X < numDetectors; ++X) {
        for (UINT4 k = 0; k < (*Fstats)->numFreqBins; ++k) {
          (*Fstats)->FaFbPerDet[X][k].Fa *= Fnorm;
          (*Fstats)->FaFbPerDet[X][k].Fb *= Fnorm;
        }
      }
    }

    // Correct F-atoms per detector arrays, and initialise to NULL
    if (whatToCompute & FSTATQ_ATOMS_PER_DET) {
      for (UINT4 k = 0; k < (*Fstats)->numFreqBins; ++k) {
        for (UINT4 X = 0; X < numDetectors; ++X) {
          FstatAtomVector *atomX = (*Fstats)->multiFatoms[k]->data[X];
          for (UINT4 alpha = 0; alpha < atomX->length; ++alpha) {
            atomX->data[alpha].Fa_alpha *= Fnorm;
            atomX->data[alpha].Fb_alpha *= Fnorm;
          }
        }
      }
    }

  } // if (common->Fnorm != 0.0)

  return XLAL_SUCCESS;

}

void
XLALDestroyFstatInputData(
  FstatInputData* input
  )
{
  if (input != NULL) {
    XLALDestroyMultiNoiseWeights(input->common.multiWeights);
    if (input->demod != NULL) {
      DestroyFstatInputData_Demod(input->demod);
    } else if (input->resamp != NULL) {
      DestroyFstatInputData_Resamp(input->resamp);
    }
    XLALFree(input);
  }
}

void
XLALDestroyFstatResults(
  FstatResults* Fstats
  )
{
  if (Fstats != NULL) {
    XLALFree(Fstats->twoF);
    XLALFree(Fstats->FaFb);
    for (UINT4 X = 0; X < PULSAR_MAX_DETECTORS; ++X) {
      XLALFree(Fstats->twoFPerDet[X]);
      XLALFree(Fstats->FaFbPerDet[X]);
      if (Fstats->multiFatoms != NULL) {
        for (UINT4 n = 0; n < Fstats->internalalloclen; ++n) {
          XLALDestroyMultiFstatAtomVector(Fstats->multiFatoms[n]);
        }
        XLALFree(Fstats->multiFatoms);
      }
    }
    XLALFree(Fstats);
  }
}

int
XLALEstimatePulsarAmplitudeParams(
  PulsarCandidate *pulsarParams,
  const LIGOTimeGPS* FaFb_refTime,
  const COMPLEX16 Fa,
  const COMPLEX16 Fb,
  const CmplxAntennaPatternMatrix *Mmunu
  )
{

  REAL8 A1h, A2h, A3h, A4h;
  REAL8 Ad, Bd, Cd, Dd, Ed;
  REAL8 normAmu;
  REAL8 A1check, A2check, A3check, A4check;

  REAL8 Asq, Da, disc;
  REAL8 aPlus, aCross;
  REAL8 Ap2, Ac2;
  REAL8 beta;
  REAL8 phi0, psi;
  REAL8 b1, b2, b3;
  REAL8 h0, cosi;

  REAL8 cosphi0, sinphi0, cos2psi, sin2psi;

  REAL8 tolerance = LAL_REAL4_EPS;

  gsl_vector *x_mu, *A_Mu;
  gsl_matrix *M_Mu_Nu;
  gsl_matrix *Jh_Mu_nu;
  gsl_permutation *permh;
  gsl_matrix *tmp, *tmp2;
  int signum;

  XLAL_CHECK ( pulsarParams != NULL, XLAL_EFAULT );
  XLAL_CHECK ( Mmunu != NULL, XLAL_EFAULT );

  Ad = Mmunu->Ad;
  Bd = Mmunu->Bd;
  Cd = Mmunu->Cd;
  Ed = Mmunu->Ed;
  Dd = Ad * Bd - Cd * Cd - Ed * Ed;

  normAmu = 2.0 / sqrt(2.0 * Mmunu->Sinv_Tsft);	/* generally *very* small!! */

  /* ----- GSL memory allocation ----- */
  XLAL_CHECK ( ( x_mu = gsl_vector_calloc (4) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( A_Mu = gsl_vector_calloc (4) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( M_Mu_Nu = gsl_matrix_calloc (4, 4) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( Jh_Mu_nu = gsl_matrix_calloc (4, 4) ) != NULL, XLAL_ENOMEM );

  XLAL_CHECK ( ( permh = gsl_permutation_calloc ( 4 ) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( tmp = gsl_matrix_calloc (4, 4) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( tmp2 = gsl_matrix_calloc (4, 4) ) != NULL, XLAL_ENOMEM );

  /* ----- fill vector x_mu */
  gsl_vector_set (x_mu, 0,   creal(Fa) );	/* x_1 */
  gsl_vector_set (x_mu, 1,   creal(Fb) );        /* x_2 */
  gsl_vector_set (x_mu, 2, - cimag(Fa) );	/* x_3 */
  gsl_vector_set (x_mu, 3, - cimag(Fb) );	/* x_4 */

  /* ----- fill matrix M^{mu,nu} [symmetric: use UPPER HALF ONLY!!]*/
  gsl_matrix_set (M_Mu_Nu, 0, 0,   Bd / Dd );
  gsl_matrix_set (M_Mu_Nu, 1, 1,   Ad / Dd );
  gsl_matrix_set (M_Mu_Nu, 0, 1, - Cd / Dd );

  gsl_matrix_set (M_Mu_Nu, 0, 3, - Ed / Dd );
  gsl_matrix_set (M_Mu_Nu, 1, 2,   Ed / Dd );

  gsl_matrix_set (M_Mu_Nu, 2, 2,   Bd / Dd );
  gsl_matrix_set (M_Mu_Nu, 3, 3,   Ad / Dd );
  gsl_matrix_set (M_Mu_Nu, 2, 3, - Cd / Dd );

  /* get (un-normalized) MLE's for amplitudes A^mu  = M^{mu,nu} x_nu */

  /* GSL-doc: int gsl_blas_dsymv (CBLAS_UPLO_t Uplo, double alpha, const gsl_matrix * A,
   *                              const gsl_vector * x, double beta, gsl_vector * y )
   *
   * compute the matrix-vector product and sum: y = alpha A x + beta y
   * for the symmetric matrix A. Since the matrix A is symmetric only its
   * upper half or lower half need to be stored. When Uplo is CblasUpper
   * then the upper triangle and diagonal of A are used, and when Uplo
   * is CblasLower then the lower triangle and diagonal of A are used.
   */
  XLAL_CHECK( gsl_blas_dsymv (CblasUpper, 1.0, M_Mu_Nu, x_mu, 0.0, A_Mu) == 0, XLAL_EERR );

  A1h = gsl_vector_get ( A_Mu, 0 );
  A2h = gsl_vector_get ( A_Mu, 1 );
  A3h = gsl_vector_get ( A_Mu, 2 );
  A4h = gsl_vector_get ( A_Mu, 3 );

  Asq = SQ(A1h) + SQ(A2h) + SQ(A3h) + SQ(A4h);
  Da = A1h * A4h - A2h * A3h;
  disc = sqrt ( SQ(Asq) - 4.0 * SQ(Da) );

  Ap2  = 0.5 * ( Asq + disc );
  aPlus = sqrt(Ap2);		/* not yet normalized */

  Ac2 = 0.5 * ( Asq - disc );
  aCross = sqrt( Ac2 );
  aCross *= MYSIGN ( Da );      /* not yet normalized */

  beta = aCross / aPlus;

  b1 =   A4h - beta * A1h;
  b2 =   A3h + beta * A2h;
  b3 = - A1h + beta * A4h ;

  psi  = 0.5 * atan ( b1 /  b2 );	/* in [-pi/4,pi/4] (gauge used also by TDS) */
  phi0 =       atan ( b2 / b3 );	/* in [-pi/2,pi/2] */

  /* Fix remaining sign-ambiguity by checking sign of reconstructed A1 */
  A1check = aPlus * cos(phi0) * cos(2.0*psi) - aCross * sin(phi0) * sin(2*psi);
  if ( A1check * A1h <  0 )
    phi0 += LAL_PI;

  cosphi0 = cos(phi0);
  sinphi0 = sin(phi0);
  cos2psi = cos(2*psi);
  sin2psi = sin(2*psi);

  /* check numerical consistency of estimated Amu and reconstructed */
  A1check =   aPlus * cosphi0 * cos2psi - aCross * sinphi0 * sin2psi;
  A2check =   aPlus * cosphi0 * sin2psi + aCross * sinphi0 * cos2psi;
  A3check = - aPlus * sinphi0 * cos2psi - aCross * cosphi0 * sin2psi;
  A4check = - aPlus * sinphi0 * sin2psi + aCross * cosphi0 * cos2psi;

  if ( ( fabs( (A1check - A1h)/A1h ) > tolerance ) ||
       ( fabs( (A2check - A2h)/A2h ) > tolerance ) ||
       ( fabs( (A3check - A3h)/A3h ) > tolerance ) ||
       ( fabs( (A4check - A4h)/A4h ) > tolerance ) )
  {
    if ( lalDebugLevel )
      XLALPrintError ( "WARNING %s(): Difference between estimated and reconstructed Amu exceeds tolerance of %g\n",
                       __func__, tolerance );
  }

  /* translate A_{+,x} into {h_0, cosi} */
  h0 = aPlus + sqrt ( disc );  /* not yet normalized ! */
  cosi = aCross / h0;

  /* ========== Estimate the errors ========== */

  /* ----- compute derivatives \partial A^\mu / \partial B^\nu, where
   * we consider the output-variables B^\nu = (h0, cosi, phi0, psi)
   * where aPlus = 0.5 * h0 * (1 + cosi^2)  and aCross = h0 * cosi
   */
  { /* Ahat^mu is defined as A^mu with the replacements: A_+ --> A_x, and A_x --> h0 */
    REAL8 A1hat =   aCross * cosphi0 * cos2psi - h0 * sinphi0 * sin2psi;
    REAL8 A2hat =   aCross * cosphi0 * sin2psi + h0 * sinphi0 * cos2psi;
    REAL8 A3hat = - aCross * sinphi0 * cos2psi - h0 * cosphi0 * sin2psi;
    REAL8 A4hat = - aCross * sinphi0 * sin2psi + h0 * cosphi0 * cos2psi;

    /* ----- A1 =   aPlus * cosphi0 * cos2psi - aCross * sinphi0 * sin2psi; ----- */
    gsl_matrix_set (Jh_Mu_nu, 0, 0,   A1h / h0 );	/* dA1/h0 */
    gsl_matrix_set (Jh_Mu_nu, 0, 1,   A1hat );          /* dA1/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 0, 2,   A3h );		/* dA1/dphi0 */
    gsl_matrix_set (Jh_Mu_nu, 0, 3, - 2.0 * A2h );	/* dA1/dpsi */

    /* ----- A2 =   aPlus * cosphi0 * sin2psi + aCross * sinphi0 * cos2psi; ----- */
    gsl_matrix_set (Jh_Mu_nu, 1, 0,   A2h / h0 );	/* dA2/h0 */
    gsl_matrix_set (Jh_Mu_nu, 1, 1,   A2hat );          /* dA2/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 1, 2,   A4h );		/* dA2/dphi0 */
    gsl_matrix_set (Jh_Mu_nu, 1, 3,   2.0 * A1h );	/* dA2/dpsi */

    /* ----- A3 = - aPlus * sinphi0 * cos2psi - aCross * cosphi0 * sin2psi; ----- */
    gsl_matrix_set (Jh_Mu_nu, 2, 0,   A3h / h0 );	/* dA3/h0 */
    gsl_matrix_set (Jh_Mu_nu, 2, 1,   A3hat );          /* dA3/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 2, 2, - A1h );		/* dA3/dphi0 */
    gsl_matrix_set (Jh_Mu_nu, 2, 3, - 2.0 * A4h );	/* dA3/dpsi */

    /* ----- A4 = - aPlus * sinphi0 * sin2psi + aCross * cosphi0 * cos2psi; ----- */
    gsl_matrix_set (Jh_Mu_nu, 3, 0,   A4h / h0 );	/* dA4/h0 */
    gsl_matrix_set (Jh_Mu_nu, 3, 1,   A4hat );          /* dA4/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 3, 2, - A2h );		/* dA4/dphi0 */
    gsl_matrix_set (Jh_Mu_nu, 3, 3,   2.0 * A3h );	/* dA4/dpsi */
  }

  /* ----- compute inverse matrices Jh^{-1} by LU-decomposition ----- */
  XLAL_CHECK( gsl_linalg_LU_decomp (Jh_Mu_nu, permh, &signum ) == 0, XLAL_EERR );

  /* inverse matrix */
  XLAL_CHECK( gsl_linalg_LU_invert (Jh_Mu_nu, permh, tmp ) == 0, XLAL_EERR );
  gsl_matrix_memcpy ( Jh_Mu_nu, tmp );

  /* ----- compute Jh^-1 . Minv . (Jh^-1)^T ----- */

  /* GSL-doc: gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha,
   *                          const gsl_matrix *A, const gsl_matrix *B, double beta, gsl_matrix *C)
   * These functions compute the matrix-matrix product and sum
   * C = \alpha op(A) op(B) + \beta C
   * where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans
   * and similarly for the parameter TransB.
   */

  /* first tmp = Minv . (Jh^-1)^T */
  XLAL_CHECK( gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, M_Mu_Nu, Jh_Mu_nu, 0.0, tmp ) == 0, XLAL_EERR );
  /* then J^-1 . tmp , store result in tmp2 */
  XLAL_CHECK( gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Jh_Mu_nu, tmp, 0.0, tmp2 ) == 0, XLAL_EERR );
  gsl_matrix_memcpy ( Jh_Mu_nu, tmp2 );

  /* ===== debug-output resulting matrices ===== */
  /* propagate initial-phase from Fstat-reference-time to refTime of Doppler-params */
  XLAL_CHECK( XLALExtrapolatePulsarPhase ( &phi0, pulsarParams->Doppler.fkdot, pulsarParams->Doppler.refTime, phi0, *FaFb_refTime )
              == XLAL_SUCCESS, XLAL_EFUNC );

  if ( phi0 < 0 )             /* make sure phi0 in [0, 2*pi] */
    phi0 += LAL_TWOPI;
  phi0 = fmod ( phi0, LAL_TWOPI );

  /* fill candidate-struct with the obtained signal-parameters and error-estimations */
  pulsarParams->Amp.h0     = normAmu * h0;
  pulsarParams->Amp.cosi   = cosi;
  pulsarParams->Amp.phi0   = phi0;
  pulsarParams->Amp.psi    = psi;

  /* read out principal estimation-errors from diagonal elements of inverse Fisher-matrix*/
  pulsarParams->dAmp.h0     = normAmu * sqrt( gsl_matrix_get (Jh_Mu_nu, 0, 0 ) );
  pulsarParams->dAmp.cosi   = sqrt( gsl_matrix_get (Jh_Mu_nu, 1, 1 ) );
  pulsarParams->dAmp.phi0   = sqrt( gsl_matrix_get (Jh_Mu_nu, 2, 2 ) );
  pulsarParams->dAmp.psi    = sqrt( gsl_matrix_get (Jh_Mu_nu, 3, 3 ) );
  /* return also the full Amplitude-Fisher matrix: invert Jh_Mu_nu */
  XLAL_CHECK( gsl_linalg_LU_decomp (Jh_Mu_nu, permh, &signum ) == 0, XLAL_EERR );
  XLAL_CHECK( gsl_linalg_LU_invert (Jh_Mu_nu, permh, tmp ) == 0, XLAL_EERR );
  pulsarParams->AmpFisherMatrix = tmp;

  /* ----- free GSL memory ----- */
  gsl_vector_free ( x_mu );
  gsl_vector_free ( A_Mu );
  gsl_matrix_free ( M_Mu_Nu );
  gsl_matrix_free ( Jh_Mu_nu );
  gsl_permutation_free ( permh );
  gsl_matrix_free ( tmp2 );

  return XLAL_SUCCESS;

}
