/*
 * Copyright (C) 2010 Reinhard Prix
 * Copyright (C) 2011, 2014 David Keitel
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

/*---------- INCLUDES ----------*/

/* GSL includes */
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <lal/SynthesizeCWDraws.h>

/*---------- DEFINES ----------*/
#define SQ(x) ((x)*(x))

/*----- SWITCHES -----*/
/*---------- internal types ----------*/

/*---------- Global variables ----------*/

/*---------- internal prototypes ----------*/

/*==================== FUNCTION DEFINITIONS ====================*/

/**
 * Generate 4 random-noise draws n_mu = {n_1, n_2, n_3, n_4} with correlations according to
 * the matrix M = L L^T, which is passed in as input.
 *
 * Note: you need to pass a pre-allocated 4-vector n_mu.
 * Note2: this function is meant as a lower-level noise-generation utility, called
 * from a higher-level function to translate the antenna-pattern functions into pre-factorized Lcor
 */
int
XLALDrawCorrelatedNoise ( PulsarAmplitudeVect n_mu,	/**< [out] 4d vector of noise-components {n_mu}, with correlation L * L^T */
                          const gsl_matrix *L,		/**< [in] correlator matrix to get n_mu = L_mu_nu * norm_nu from 4 uncorr. unit variates norm_nu */
                          gsl_rng * rng			/**< gsl random-number generator */
                          )
{
  int gslstat;

  /* ----- check input arguments ----- */
  if ( !L || (L->size1 != 4) || (L->size2 != 4) ) {
    XLALPrintError ( "%s: Invalid correlator matrix, n_mu must be pre-allocate 4x4 matrix", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( !rng ) {
    XLALPrintError ("%s: invalid NULL input as gsl random-number generator!\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  /* ----- generate 4 normal-distributed, uncorrelated random numbers ----- */
  PulsarAmplitudeVect n;

  n[0] = gsl_ran_gaussian ( rng, 1.0 );
  n[1] = gsl_ran_gaussian ( rng, 1.0 );
  n[2] = gsl_ran_gaussian ( rng, 1.0 );
  n[3] = gsl_ran_gaussian ( rng, 1.0 );

  /* use four normal-variates {norm_nu} with correlator matrix L to get: n_mu = L_{mu nu} norm_nu
   * which gives {n_\mu} satisfying cov(n_mu,n_nu) = (L L^T)_{mu nu} = M_{mu nu}
   */
  gsl_vector_view n_view = gsl_vector_view_array ( n, 4 );
  gsl_vector_view n_mu_view = gsl_vector_view_array ( n_mu, 4 );

  /* int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)
   * compute the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H
   * for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
   */
  if ( (gslstat = gsl_blas_dgemv (CblasNoTrans, 1.0, L, &n_view.vector, 0.0, &n_mu_view.vector)) != 0 ) {
    XLALPrintError ( "%s: gsl_blas_dgemv(L * norm) failed: %s\n", __func__, gsl_strerror (gslstat) );
    XLAL_ERROR ( XLAL_EFAILED );
  }

  return XLAL_SUCCESS;

} /* XLALDrawCorrelatedNoise() */

/**
 * Generate an FstatAtomVector for given antenna-pattern functions.
 * Simply creates FstatAtomVector and initializes with antenna-pattern function.
 */
FstatAtomVector*
XLALGenerateFstatAtomVector ( const DetectorStateSeries *detStates,	/**< input detector-state series, only used for timestamps */
                              const AMCoeffs *amcoeffs			/**< input antenna-pattern functions {a_i, b_i} */
                              )
{
  /* check input consistency */
  if ( !detStates || !detStates->data ) {
    XLALPrintError ("%s: invalid NULL input in detStates=%p or detStates->data=%p\n", __func__, detStates, detStates->data );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }
  if ( !amcoeffs || !amcoeffs->a || !amcoeffs->b ) {
    XLALPrintError ("%s: invalid NULL input in amcoeffs=%p or amcoeffs->a=%p, amcoeffs->b=%p\n", __func__, amcoeffs, amcoeffs->a, amcoeffs->b );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }
  UINT4 numAtoms = detStates->length;
  if ( numAtoms != amcoeffs->a->length || numAtoms != amcoeffs->b->length ) {
    XLALPrintError ("%s: inconsistent lengths numAtoms=%d amcoeffs->a = %d, amecoeffs->b = %d\n", __func__, numAtoms, amcoeffs->a->length, amcoeffs->b->length );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /* prepare output vector */
  FstatAtomVector *atoms;
  if ( ( atoms = XLALCreateFstatAtomVector ( numAtoms ) ) == NULL ) {
    XLALPrintError ("%s: XLALCreateFstatAtomVector(%d) failed.\n", __func__, numAtoms );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  atoms->TAtom = detStates->deltaT;
  UINT4 alpha;
  for ( alpha=0; alpha < numAtoms; alpha ++ )
    {
      REAL8 a = amcoeffs->a->data[alpha];
      REAL8 b = amcoeffs->b->data[alpha];

      atoms->data[alpha].timestamp = detStates->data[alpha].tGPS.gpsSeconds - 0.5*detStates->deltaT;	/* shift back to SFT start-time */
      atoms->data[alpha].a2_alpha = a * a;
      atoms->data[alpha].b2_alpha = b * b;
      atoms->data[alpha].ab_alpha = a * b;

      /* Fa,Fb are zero-initialized from XLALCreateFstatAtomVector() */

    } /* for alpha < numAtoms */

  /* return result */
  return atoms;

} /* XLALGenerateFstatAtomVector() */


/**
 * Generate a multi-FstatAtomVector for given antenna-pattern functions.
 * Simply creates MultiFstatAtomVector and initializes with antenna-pattern function.
 */
MultiFstatAtomVector*
XLALGenerateMultiFstatAtomVector ( const MultiDetectorStateSeries *multiDetStates, 	/**< [in] multi-detector state series, only used for timestamps */
                                   const MultiAMCoeffs *multiAM				/**< input antenna-pattern functions {a_i, b_i} */
                                   )
{
  /* check input consistency */
  if ( !multiDetStates || !multiDetStates->data ) {
    XLALPrintError ("%s: invalid NULL input in 'multiDetStates'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }
  if ( !multiAM || !multiAM->data || !multiAM->data[0] ) {
    XLALPrintError ("%s: invalid NULL input in 'mutiAM'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  UINT4 numDet = multiDetStates->length;
  if ( numDet != multiAM->length ) {
    XLALPrintError ("%s: inconsistent number of detectors in multiDetStates (%d) and multiAM (%d)\n", __func__, multiDetStates->length, multiAM->length );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /* create multi-atoms vector */
  MultiFstatAtomVector *multiAtoms;
  if ( ( multiAtoms = XLALCalloc ( 1, sizeof(*multiAtoms) )) == NULL ) {
    XLALPrintError ("%s: XLALCalloc ( 1, %d) failed.\n", __func__, sizeof(*multiAtoms) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  multiAtoms->length = numDet;
  if ( ( multiAtoms->data = XLALCalloc ( numDet, sizeof(*multiAtoms->data) ) ) == NULL ) {
    XLALPrintError ("%s: XLALCalloc ( %d, %d) failed.\n", __func__, numDet, sizeof(*multiAtoms->data) );
    XLALFree ( multiAtoms );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  /* loop over detectors and generate each atoms-vector individually */
  UINT4 X;
  for ( X=0; X < numDet; X ++ )
    {
      if ( ( multiAtoms->data[X] = XLALGenerateFstatAtomVector ( multiDetStates->data[X], multiAM->data[X] )) == NULL ) {
        XLALPrintError ("%s: XLALGenerateFstatAtomVector() failed.\n", __func__ );
        XLALDestroyMultiFstatAtomVector ( multiAtoms );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }

    } /* for X < numDet */

  /* return result */
  return multiAtoms;

} /* XLALGenerateMultiFstatAtomVector() */

/**
 * Add Gaussian-noise components to given FstatAtomVector
 */
int
XLALAddNoiseToFstatAtomVector ( FstatAtomVector *atoms,	/**< input atoms-vector, noise will be added to this */
                                gsl_rng * rng		/**< random-number generator */
                                )
{
  /* check input consistency */
  if ( !atoms || !rng ) {
    XLALPrintError ("%s: invalid NULL input for 'atoms'=%p or random-number generator 'rng'=%p\n", __func__, atoms, rng );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  UINT4 numAtoms = atoms->length;

  /* prepare gsl-matrix for correlator L = 1/2 * [ a, a ; b , b ] */
  gsl_matrix *Lcor;
  if ( (Lcor = gsl_matrix_calloc ( 4, 4 )) == NULL ) {
    XLALPrintError ("%s: gsl_matrix_calloc ( 4, 4 ) failed.\n", __func__ );
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  /* prepare placeholder for 4 n_mu noise draws */
  PulsarAmplitudeVect n_mu;

  /* ----- step through atoms and synthesize noise ----- */
  UINT4 alpha;
  for ( alpha=0; alpha < numAtoms; alpha ++ )
    {
      /* unfortunately we need {a,b} here, but
       * the atoms only store {a^2, b^2, ab }
       * so we need to invert this [module arbitrary relative sign)
       */
      REAL8 a2 = atoms->data[alpha].a2_alpha;
      REAL8 b2 = atoms->data[alpha].b2_alpha;
      REAL8 ab = atoms->data[alpha].ab_alpha;

      REAL8 a_by_2 = 0.5 * sqrt(a2);
      REAL8 b_by_2 = 0.5 * sqrt(b2);
      /* convention: always set sign on b */
      if ( ab < 0 )
        b_by_2 = - b_by_2;

      /* upper-left block */
      gsl_matrix_set ( Lcor, 0, 0, a_by_2 );
      gsl_matrix_set ( Lcor, 0, 1, a_by_2 );
      gsl_matrix_set ( Lcor, 1, 0, b_by_2);
      gsl_matrix_set ( Lcor, 1, 1, b_by_2);
      /* lower-right block: +2 on all components */
      gsl_matrix_set ( Lcor, 2, 2, a_by_2 );
      gsl_matrix_set ( Lcor, 2, 3, a_by_2 );
      gsl_matrix_set ( Lcor, 3, 2, b_by_2);
      gsl_matrix_set ( Lcor, 3, 3, b_by_2);

      if ( XLALDrawCorrelatedNoise ( n_mu, Lcor, rng ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: failed to XLALDrawCorrelatedNoise().\n", __func__ );
        XLAL_ERROR ( XLAL_EFUNC );
      }
      REAL8 x1,x2,x3,x4;
      x1 = n_mu[0];
      x2 = n_mu[1];
      x3 = n_mu[2];
      x4 = n_mu[3];

      /* add this to Fstat-atom */
      /* relation Fa,Fb <--> x_mu: see Eq.(72) in CFSv2-LIGO-T0900149-v2.pdf */
      atoms->data[alpha].Fa_alpha += crectf( x1, - x3 );
      atoms->data[alpha].Fb_alpha += crectf( x2, - x4 );

    } /* for i < numAtoms */

  /* free internal memory */
  gsl_matrix_free ( Lcor );

  return XLAL_SUCCESS;

} /* XLALAddNoiseToFstatAtomVector() */


/**
 * Add Gaussian-noise components to given multi-FstatAtomVector
 */
int
XLALAddNoiseToMultiFstatAtomVector ( MultiFstatAtomVector *multiAtoms,	/**< input multi atoms-vector, noise will be added to this */
                                     gsl_rng * rng			/**< random-number generator */
                                     )
{
  /* check input consistency */
  if ( !multiAtoms || !multiAtoms->data ) {
    XLALPrintError ("%s: invalid NULL input in 'multiAtoms'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( !rng ) {
    XLALPrintError ("%s: invalid NULL input for random-number generator 'rng'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  UINT4 numDetectors = multiAtoms->length;

  UINT4 X;
  for ( X=0; X < numDetectors; X ++ )
    {
      if ( XLALAddNoiseToFstatAtomVector ( multiAtoms->data[X], rng ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: XLALAddNoiseToFstatAtomVector() failed.\n", __func__ );
        XLAL_ERROR ( XLAL_EFUNC );
      }

    } /* for X < numDetectors */

  return XLAL_SUCCESS;

} /* XLALSynthesizeMultiFstatAtomVector4Noise() */


/**
 * Add signal s_mu = M_mu_nu A^nu within the given transient-window
 * to given atoms.
 *
 * RETURN: SNR^2 of the injected signal
 * and the effective AntennaPatternMatrix M_mu_nu for this signal.
 */
REAL8
XLALAddSignalToFstatAtomVector ( FstatAtomVector* atoms,	 /**< [in/out] atoms vectors containing antenna-functions and possibly noise {Fa,Fb} */
                                 AntennaPatternMatrix *M_mu_nu,	 /**< [out] effective antenna-pattern matrix for the injected signal */
                                 const PulsarAmplitudeVect A_Mu, /**< [in] input canonical amplitude vector A^mu = {A1,A2,A3,A4} */
                                 transientWindow_t transientWindow /**< transient signal window */
                                 )
{
  int gslstat;

  /* check input consistency */
  if ( !atoms || !atoms->data ) {
    XLALPrintError ( "%s: Invalid NULL input 'atoms'\n", __func__ );
    XLAL_ERROR_REAL8 ( XLAL_EINVAL );
  }
  if ( !M_mu_nu ) {
    XLALPrintError ( "%s: Invalid NULL input 'M_mu_nu'\n", __func__ );
    XLAL_ERROR_REAL8 ( XLAL_EINVAL );
  }

  /* prepare transient-window support */
  UINT4 t0, t1;
  if ( XLALGetTransientWindowTimespan ( &t0, &t1, transientWindow ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALGetTransientWindowTimespan() failed.\n", __func__ );
    XLAL_ERROR_REAL8 ( XLAL_EFUNC );
  }

  /* prepare gsl-matrix for Mh_mu_nu = [ a^2, a*b ; a*b , b^2 ] */
  gsl_matrix *Mh_mu_nu;
  if ( (Mh_mu_nu = gsl_matrix_calloc ( 4, 4 )) == NULL ) {
    XLALPrintError ("%s: gsl_matrix_calloc ( 4, 4 ) failed.\n", __func__ );
    XLAL_ERROR_REAL8 ( XLAL_ENOMEM );
  }

  gsl_vector_const_view A_Mu_view = gsl_vector_const_view_array ( A_Mu, 4 );

  REAL8 TAtom = atoms->TAtom;
  UINT4 numAtoms = atoms->length;
  UINT4 alpha;
  REAL8 Ad = 0, Bd = 0, Cd = 0;		// usual non-transient antenna-pattern functions
  REAL8 Ap = 0, Bp = 0, Cp = 0;		// "effective" transient-CW antenna-pattern matrix M'_mu_nu

  for ( alpha=0; alpha < numAtoms; alpha ++ )
    {
      UINT4 ti = atoms->data[alpha].timestamp;
      REAL8 win = XLALGetTransientWindowValue ( ti, t0, t1, transientWindow.tau, transientWindow.type );

      if ( win == 0 )
        continue;

      /* compute sh_mu = sqrt(gamma/2) * Mh_mu_nu A^nu * win, where Mh_mu_nu is now just
       * the per-atom block matrix [a^2,  ab; ab, b^2 ]
       * where Sn=1, so gamma = Sinv*TAtom = TAtom
       */
      // NOTE: for sh_mu: only LINEAR in window-function, NOT quadratic! -> see notes
      REAL8 a2 = win * atoms->data[alpha].a2_alpha;
      REAL8 b2 = win * atoms->data[alpha].b2_alpha;
      REAL8 ab = win * atoms->data[alpha].ab_alpha;

      Ad += a2;
      Bd += b2;
      Cd += ab;

      // we also compute M'_mu_nu, which will be used to estimate optimal SNR
      // NOTE: M'_mu_nu is QUADRATIC in window-function!, so we multiply by win again
      Ap += win * a2;
      Bp += win * b2;
      Cp += win * ab;

      /* upper-left block */
      gsl_matrix_set ( Mh_mu_nu, 0, 0, a2 );
      gsl_matrix_set ( Mh_mu_nu, 1, 1, b2 );
      gsl_matrix_set ( Mh_mu_nu, 0, 1, ab );
      gsl_matrix_set ( Mh_mu_nu, 1, 0, ab );
      /* lower-right block: +2 on all components */
      gsl_matrix_set ( Mh_mu_nu, 2, 2, a2 );
      gsl_matrix_set ( Mh_mu_nu, 3, 3, b2 );
      gsl_matrix_set ( Mh_mu_nu, 2, 3, ab );
      gsl_matrix_set ( Mh_mu_nu, 3, 2, ab );

      /* placeholder for 4-vector sh_mu */
      PulsarAmplitudeVect sh_mu = {0,0,0,0};
      gsl_vector_view sh_mu_view = gsl_vector_view_array ( sh_mu, 4 );

      /* int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)
       * compute the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H
       * for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
       *
       * sh_mu = sqrt(gamma/2) * Mh_mu_nu A^nu, where here gamma = TAtom, as Sinv=1 for multi-IFO value, and weights for SinvX!=1 have already been absorbed in atoms through XLALComputeMultiAMCoeffs()
       */
      REAL8 norm_s = sqrt(TAtom / 2.0);
      if ( (gslstat = gsl_blas_dgemv (CblasNoTrans, norm_s, Mh_mu_nu, &A_Mu_view.vector, 0.0, &sh_mu_view.vector)) != 0 ) {
        XLALPrintError ( "%s: gsl_blas_dgemv(L * norm) failed: %s\n", __func__, gsl_strerror (gslstat) );
        XLAL_ERROR_REAL8 ( XLAL_EFAILED );
      }

      /* add this signal to the atoms, using the relation Fa,Fb <--> x_mu: see Eq.(72) in CFSv2-LIGO-T0900149-v2.pdf */
      REAL8 s1,s2,s3,s4;
      s1 = sh_mu[0];
      s2 = sh_mu[1];
      s3 = sh_mu[2];
      s4 = sh_mu[3];

      atoms->data[alpha].Fa_alpha += crectf( s1, - s3 );
      atoms->data[alpha].Fb_alpha += crectf( s2, - s4 );

    } /* for alpha < numAtoms */

  /* compute optimal SNR^2 expected for this signal,
   * using rho2 = A^mu M'_mu_nu A^nu = T/Sn( A' [A1^2+A3^2] + 2C' [A1A2 +A3A4] + B' [A2^2+A4^2])
   * NOTE: here we use the "effective" transient-CW antenna-pattern matrix M'_mu_nu
   */
  REAL8 A1 = A_Mu[0];
  REAL8 A2 = A_Mu[1];
  REAL8 A3 = A_Mu[2];
  REAL8 A4 = A_Mu[3];

  REAL8 rho2 = TAtom  * ( Ap * ( SQ(A1) + SQ(A3) ) + 2.0*Cp * ( A1*A2 + A3*A4 ) + Bp * ( SQ(A2) + SQ(A4) ) );

  /* return "effective" transient antenna-pattern matrix */
  M_mu_nu->Sinv_Tsft = TAtom;	/* everything here in units of Sn, so effectively Sn=1 */
  M_mu_nu->Ad = Ap;
  M_mu_nu->Bd = Bp;
  M_mu_nu->Cd = Cp;
  M_mu_nu->Dd = Ap * Bp - Cp * Cp;

  /* free memory */
  gsl_matrix_free ( Mh_mu_nu );

  /* return SNR^2 */
  return rho2;

} /* XLALAddSignalToFstatAtomVector() */


/**
 * Add given signal s_mu = M_mu_nu A^nu within the given transient-window
 * to multi-IFO noise-atoms.
 *
 * RETURN: SNR^2 of the injected signal
 * and the effective AntennaPatternMatrix M_mu_nu for this signal.
 */
REAL8
XLALAddSignalToMultiFstatAtomVector ( MultiFstatAtomVector* multiAtoms,	 /**< [in/out] multi atoms vectors containing antenna-functions and possibly noise {Fa,Fb} */
                                      AntennaPatternMatrix *M_mu_nu,	 /**< [out] effective multi-IFO antenna-pattern matrix for the injected signal */
                                      const PulsarAmplitudeVect A_Mu, 	/**< [in] input canonical amplitude vector A^mu = {A1,A2,A3,A4} */
                                      transientWindow_t transientWindow, /**< [in] transient signal window */
                                      INT4 lineX	 		/**< [in] if >= 0: generate signal only for detector 'lineX': must be within 0,...(Ndet-1) */
                                      )
{
  /* check input consistency */
  if ( !multiAtoms || !multiAtoms->data ) {
    XLALPrintError ( "%s: Invalid NULL input 'multiAtoms'\n", __func__ );
    XLAL_ERROR_REAL8 ( XLAL_EINVAL );
  }
  if ( !M_mu_nu ) {
    XLALPrintError ( "%s: Invalid NULL input 'M_mu_nu'\n", __func__ );
    XLAL_ERROR_REAL8 ( XLAL_EINVAL );
  }

  UINT4 numDet = multiAtoms->length;
  XLAL_CHECK ( (lineX < 0) || ((UINT4)lineX < numDet), XLAL_EINVAL, "Inconsistent input of lineX = %d, not within 0 ... Ndet-1 (= %d)\n", lineX, numDet );

  UINT4 X;
  REAL8 rho2 = 0;
  XLAL_INIT_MEM( (*M_mu_nu));

  for ( X=0; X < numDet; X ++ )
    {
      REAL8 rho2X;
      AntennaPatternMatrix M_mu_nu_X;
      PulsarAmplitudeVect A0_Mu = {0,0,0,0};	// zero amplitude signal for simulating single-IFO line

      if ( (lineX >= 0) && ((UINT4)lineX != X) ) {
        rho2X = XLALAddSignalToFstatAtomVector ( multiAtoms->data[X], &M_mu_nu_X, A0_Mu, transientWindow );	// zero-signal injection
      }
      else {
        rho2X = XLALAddSignalToFstatAtomVector ( multiAtoms->data[X], &M_mu_nu_X, A_Mu, transientWindow );	// actual signal injection
      }
      XLAL_CHECK_REAL8 ( xlalErrno == 0, XLAL_EFUNC );

      rho2 += rho2X;			/* multi-IFO SNR^2 = sum_X SNR_X^2 */
      M_mu_nu->Ad += M_mu_nu_X.Ad;	/* multi-IFO M_mu_nu = sum_X M_mu_nu_X */
      M_mu_nu->Bd += M_mu_nu_X.Bd;
      M_mu_nu->Cd += M_mu_nu_X.Cd;

      M_mu_nu->Sinv_Tsft += M_mu_nu_X.Sinv_Tsft;	/* noise adds harmonically 1/S = sum_X (1/S_X) */

    } /* for X < numDet */

  /* update sub-determinant */
  M_mu_nu->Dd = M_mu_nu->Ad * M_mu_nu->Bd - SQ(M_mu_nu->Cd);

  /* return SNR^2 */
  return rho2;

} /* XLALAddSignalToMultiFstatAtomVector() */

/**
 * Generates a multi-Fstat atoms vector for given parameters, drawing random parameters wherever required.
 *
 * Input: detector states, signal sky-pos (or allsky), amplitudes (or range), transient window range
 *
 */
MultiFstatAtomVector *
XLALSynthesizeTransientAtoms ( InjParams_t *injParamsOut,			/**< [out] return summary of injected signal parameters (can be NULL) */
                               SkyPosition skypos,				/**< (Alpha,Delta,system). Use Alpha < 0 to signal 'allsky' */
                               AmplitudePrior_t AmpPrior,			/**< [in] amplitude-parameter priors to draw signals from */
                               transientWindowRange_t transientInjectRange,	/**< transient-window range for injections (flat priors) */
                               const MultiDetectorStateSeries *multiDetStates, 	/**< [in] multi-detector state series covering observation time */
                               BOOLEAN SignalOnly,				/**< [in] switch to generate signal draws without noise */
                               multiAMBuffer_t *multiAMBuffer,			/**< [in/out] buffer for AM-coefficients if re-using same skyposition (must be !=NULL) */
                               gsl_rng *rng,					/**< [in/out] gsl random-number generator */
                               INT4 lineX,					/**< [in] if >= 0: generate signal only for detector 'lineX': must be within 0,...(Ndet-1) */
                               const MultiNoiseWeights *multiNoiseWeights	/** [in] per-detector noise weights SX^-1/S^-1, no per-SFT variation (can be NULL for unit weights) */
                               )
{
  /* check input */
  XLAL_CHECK_NULL ( rng && multiAMBuffer && multiDetStates, XLAL_EINVAL, "Invalid NULL input!\n");
  XLAL_CHECK_NULL ( !multiNoiseWeights || multiNoiseWeights->data, XLAL_EINVAL, "Invalid NULL input for multiNoiseWeights->data!\n" );

  /* -----  if Alpha < 0 ==> draw skyposition isotropically from all-sky */
  if ( skypos.longitude < 0 )
    {
      skypos.longitude = gsl_ran_flat ( rng, 0, LAL_TWOPI );	/* alpha uniform in [ 0, 2pi ] */
      skypos.latitude  = acos ( gsl_ran_flat ( rng, -1, 1 ) ) - LAL_PI_2;	/* cos(delta) uniform in [ -1, 1 ] */
      skypos.system    = COORDINATESYSTEM_EQUATORIAL;
      /* never re-using buffered AM-coeffs here, as we randomly draw new skypositions */
      if ( multiAMBuffer->multiAM ) XLALDestroyMultiAMCoeffs ( multiAMBuffer->multiAM );
      multiAMBuffer->multiAM = NULL;
    } /* if random skypos to draw */
  else /* otherwise: re-use AM-coeffs if already computed, or initialize them if for the first time */
    {
      if ( multiAMBuffer->multiAM )
        if ( multiAMBuffer->skypos.longitude != skypos.longitude ||
             multiAMBuffer->skypos.latitude  != skypos.latitude  ||
             multiAMBuffer->skypos.system    != skypos.system )
          {
            XLALDestroyMultiAMCoeffs ( multiAMBuffer->multiAM );
            multiAMBuffer->multiAM = NULL;
          }
    } /* if single skypos given */

  /* ----- generate antenna-pattern functions for this sky-position */
  if ( !multiAMBuffer->multiAM && (multiAMBuffer->multiAM = XLALComputeMultiAMCoeffs ( multiDetStates, multiNoiseWeights, skypos )) == NULL ) {
    XLALPrintError ( "%s: XLALComputeMultiAMCoeffs() failed with xlalErrno = %d\n", __func__, xlalErrno );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  multiAMBuffer->skypos = skypos; /* store buffered skyposition */

  /* ----- generate a pre-initialized F-stat atom vector containing only the antenna-pattern coefficients */
  MultiFstatAtomVector *multiAtoms;
  if ( (multiAtoms = XLALGenerateMultiFstatAtomVector ( multiDetStates, multiAMBuffer->multiAM )) == NULL ) {
    XLALPrintError ( "%s: XLALGenerateMultiFstatAtomVector() failed with xlalErrno = %d\n", __func__, xlalErrno );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  /* ----- draw amplitude vector A^mu from given ranges in {h0, cosi, psi, phi0} */
  PulsarAmplitudeParams Amp;
  if ( AmpPrior.fixedSNR > 0 )	/* special treatment of fixed-SNR: use h0=1, later rescale signal */
    Amp.h0 = 1.0;
  else if ( AmpPrior.fixedSNR == 0 )/* same as setting h0 = 0 */
    Amp.h0 = 0;
  else					/* otherwise, draw from h0-prior */
    Amp.h0 = XLALDrawFromPDF1D ( AmpPrior.pdf_h0Nat, rng );

  Amp.cosi = XLALDrawFromPDF1D ( AmpPrior.pdf_cosi, rng );
  Amp.psi  = XLALDrawFromPDF1D ( AmpPrior.pdf_psi,  rng );
  Amp.phi0 = XLALDrawFromPDF1D ( AmpPrior.pdf_phi0, rng );

  /* convert amplitude params to 'canonical' vector coordinates */
  PulsarAmplitudeVect A_Mu;
  if ( XLALAmplitudeParams2Vect ( A_Mu, Amp ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALAmplitudeParams2Vect() failed with xlalErrno = %d\n", __func__, xlalErrno );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  /* ----- draw transient-window parameters from given ranges using flat priors */
  transientWindow_t XLAL_INIT_DECL(injectWindow);
  injectWindow.type = transientInjectRange.type;
  if ( injectWindow.type != TRANSIENT_NONE )	/* nothing to be done if no window */
    {
      injectWindow.t0  = (UINT4) gsl_ran_flat ( rng, transientInjectRange.t0, transientInjectRange.t0 + transientInjectRange.t0Band );
      injectWindow.tau = (UINT4) gsl_ran_flat ( rng, transientInjectRange.tau, transientInjectRange.tau + transientInjectRange.tauBand );
    }

  /* ----- add transient signal to the Fstat atoms */
  AntennaPatternMatrix M_mu_nu;
  REAL8 rho2 = XLALAddSignalToMultiFstatAtomVector ( multiAtoms, &M_mu_nu, A_Mu, injectWindow, lineX );
  if ( xlalErrno ) {
    XLALPrintError ( "%s: XLALAddSignalToMultiFstatAtomVector() failed with xlalErrno = %d\n", __func__, xlalErrno );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  /* ----- special signal rescaling if 1) fixedSNR OR 2) fixed rhohMax */
  REAL8 detM1o8 = sqrt ( M_mu_nu.Sinv_Tsft ) * pow ( M_mu_nu.Dd, 0.25 );	// (detM)^(1/8) = sqrt(Tsft/Sn) * (Dp)^(1/4)

  /* 1) if fixedSNR signal is requested: rescale everything to the desired SNR now */
  if ( (AmpPrior.fixedSNR > 0) || AmpPrior.fixRhohMax )
    {
      if ( (AmpPrior.fixedSNR > 0) && AmpPrior.fixRhohMax ) { /* double-check consistency: only one is allowed */
        XLALPrintError ("%s: Something went wrong: both [fixedSNR = %f > 0] and [fixedRhohMax==true] are not allowed!\n", __func__, AmpPrior.fixedSNR );
        XLAL_ERROR_NULL ( XLAL_EDOM );
      }

      REAL8 rescale = 1.0;

      if ( AmpPrior.fixedSNR > 0 )
        rescale = AmpPrior.fixedSNR / sqrt(rho2);	// rescale atoms by this factor, such that SNR = fixedSNR
      if ( AmpPrior.fixRhohMax )
        rescale = 1.0 / detM1o8;	// we drew h0 in [0, rhohMax], so we now need to rescale as h0Max = rhohMax/(detM)^(1/8)

      if ( XLALRescaleMultiFstatAtomVector ( multiAtoms, rescale ) != XLAL_SUCCESS ) {	      /* rescale atoms */
        XLALPrintError ( "%s: XLALRescaleMultiFstatAtomVector() failed with xlalErrno = %d\n", __func__, xlalErrno );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }

      Amp.h0 *= rescale;	      /* rescale amplitude-params for consistency */
      UINT4 i; for (i=0; i < 4; i ++) A_Mu[i] *= rescale;

      rho2 *= SQ(rescale);	      /* rescale reported optimal SNR */

    } /* if fixedSNR > 0 OR fixedRhohMax */

  /* ----- add noise to the Fstat atoms, unless --SignalOnly was specified */
  if ( !SignalOnly )
    if ( XLALAddNoiseToMultiFstatAtomVector ( multiAtoms, rng ) != XLAL_SUCCESS ) {
      XLALPrintError ("%s: XLALAddNoiseToMultiFstatAtomVector() failed with xlalErrno = %d\n", __func__, xlalErrno );
      XLAL_ERROR_NULL ( XLAL_EFUNC );
    }

  /* ----- if requested: return all inject signal parameters */
  if ( injParamsOut )
    {
      injParamsOut->skypos = skypos;
      injParamsOut->ampParams = Amp;
      UINT4 i; for (i=0; i < 4; i ++) injParamsOut->ampVect[i] = A_Mu[i];
      injParamsOut->multiAM = *multiAMBuffer->multiAM;
      injParamsOut->transientWindow = injectWindow;
      injParamsOut->SNR = sqrt(rho2);
      injParamsOut->detM1o8 = detM1o8;
    } /* if injParamsOut */


  return multiAtoms;

} /* XLALSynthesizeTransientAtoms() */

/**
 * Rescale a given multi-Fstat atoms vector {Fa,Fb} by given scalar factor.
 * This is used to rescale signal atoms to desired fixed SNR.
 */
int
XLALRescaleMultiFstatAtomVector ( MultiFstatAtomVector* multiAtoms,	/**< [in/out] multi atoms vectors containing a signal in {Fa,Fb} to be rescaled */
                                  REAL8 rescale				/**< rescale factor: Fa' = rescale * Fa, and Fb'= rescale * Fb */
                                  )
{
  /* check input */
  if ( !multiAtoms ) {
    XLALPrintError ("%s: invalid NULL input 'multiAtoms'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  UINT4 numDet = multiAtoms->length;
  UINT4 X;

  for ( X=0; X < numDet; X ++ )
    {
      FstatAtomVector *atoms = multiAtoms->data[X];
      UINT4 numAtoms = atoms->length;
      UINT4 i;
      for ( i=0; i < numAtoms; i ++ )
        {
          FstatAtom *thisAtom = &atoms->data[i];

          thisAtom->Fa_alpha *= ((REAL4) rescale);

          thisAtom->Fb_alpha *= ((REAL4) rescale);

        } /* for i < numAtoms */

    } /* for X < numDet */

  return XLAL_SUCCESS;

} /* XLALRescaleMultiFstatAtomVector() */


/**
 * Write an injection-parameters structure to the given file-pointer,
 * adding one line with the injection parameters
 */
int
write_InjParams_to_fp ( FILE * fp,			/**< [in] file-pointer to output file */
                        const InjParams_t *par,		/**< [in] injection params to write. NULL means write header-comment line */
                        const UINT4 dataStartGPS,	/**< [in] data start-time in GPS seconds (used to turn window 't0' into offset from dataStartGPS) */
                        const BOOLEAN outputMmunuX,	/**< [in] write per-IFO antenna pattern matrices? */
                        const UINT4 numDetectors	/**< [in] number of detectors, only needed to construct M_mu_nu_X_header_string */
                        )
{
  /* input consistency */
  if ( ! fp ) {
    XLALPrintError ("%s: invalid NULL input 'fp'\n", __func__);
    XLAL_ERROR ( XLAL_EINVAL );
  }

  int ret;
  /* if requested, write header-comment line */
  if ( par == NULL ) {
    char M_mu_nu_X_header_string[256] = "";
    if ( outputMmunuX ) {
      char buf0[256];
      for ( UINT4 X = 0; X < numDetectors ; X ++ ) {
        snprintf ( buf0, sizeof(buf0), "      AdX[%d]     BdX[%d]     CdX[%d]     DdX[%d]", X, X, X, X );
        strcat ( M_mu_nu_X_header_string, buf0 );
      }
    }
    ret = fprintf(fp, "%%%%Alpha Delta      SNR       h0   cosi    psi   phi0        A1        A2        A3        A4         Ad         Bd         Cd         Dd        t0[d]      tau[d]   type (detMp)^(1/8)%s\n", M_mu_nu_X_header_string);
    if ( ret < 0 ) {
      XLALPrintError ("%s: failed to fprintf() to given file-pointer 'fp'.\n", __func__ );
      XLAL_ERROR ( XLAL_EIO );
    }

    return XLAL_SUCCESS;	/* we're done here */

  } /* if par == NULL */

  REAL8 t0_d = 1.0 * ( par->transientWindow.t0 - dataStartGPS ) / DAY24;
  REAL8 tau_d = 1.0 * par->transientWindow.tau / DAY24;

  char M_mu_nu_X_string[256] = "";
  if ( outputMmunuX ) {
    if ( par->multiAM.length != numDetectors ) {
      XLALPrintError ("%s: length of multiAM different than numDetectors (%d!=%d).\n", __func__, par->multiAM.length, numDetectors );
      XLAL_ERROR ( XLAL_EINVAL );
    }
    strcat ( M_mu_nu_X_string, "        ");
    char buf0[256];
    for ( UINT4 X = 0; X < numDetectors ; X ++ ) {
      snprintf ( buf0, sizeof(buf0), " %10.5g %10.5g %10.5g %10.5g", par->multiAM.data[X]->A, par->multiAM.data[X]->B, par->multiAM.data[X]->C, par->multiAM.data[X]->D );
      strcat ( M_mu_nu_X_string, buf0 );
    }
  }

  /* if injParams given, output them to the file */
  ret = fprintf ( fp, " %5.3f %6.3f   %6.3f  %7.3g %6.3f %6.3f %6.3f % 9.3g % 9.3g % 9.3g % 9.3g %10.5g %10.5g %10.5g %10.5g    %9.3f  %9.3f    %1d %9.3g%s\n",
                  par->skypos.longitude, par->skypos.latitude,						/* skypos */
                  par->SNR,										/* SNR */
                  par->ampParams.h0, par->ampParams.cosi, par->ampParams.psi, par->ampParams.phi0,	/* amplitude params {h0,cosi,psi,phi0}*/
                  par->ampVect[0], par->ampVect[1], par->ampVect[2], par->ampVect[3],			/* ampltiude vector A^mu */
                  par->multiAM.Mmunu.Ad, par->multiAM.Mmunu.Bd, par->multiAM.Mmunu.Cd, par->multiAM.Mmunu.Dd,	/* antenna-pattern matrix components */
                  t0_d, tau_d, par->transientWindow.type,		/* transient-window params */
                  par->detM1o8,										/* rescale parameter (detMp)^(1/8) */
                  M_mu_nu_X_string										/* optional output string for per-IFO antenna pattern matrices */
                  );
  if ( ret < 0 ) {
    XLALPrintError ("%s: failed to fprintf() to given file-pointer 'fp'.\n", __func__ );
    XLAL_ERROR ( XLAL_EIO );
  }

 return XLAL_SUCCESS;

} /* write_InjParams_to_fp() */

