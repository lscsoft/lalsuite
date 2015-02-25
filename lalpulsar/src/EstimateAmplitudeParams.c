//
// Copyright (C) 2014 Reinhard Prix
// Copyright (C) 2012, 2013, 2014 Karl Wette
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

// ---------- includes ----------
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

#include <lal/ExtrapolatePulsarSpins.h>

#include <lal/EstimateAmplitudeParams.h>

// ---------- local macro definitions ----------
#define SQ(x) ( (x) * (x) )
#define MYSIGN(x) ( ((x) < 0) ? (-1.0):(+1.0) )


// ==================== function definitions ====================

///
/// Estimate the amplitude parameters of a pulsar CW signal, given its phase parameters,
/// constituent parts of the \f$\mathcal{F}\f$-statistic, and antenna pattern matrix.
///
/// \note Parameter-estimation based on large parts on Yousuke's notes and implemention (in CFSv1),
/// extended for error-estimation.
///
int
XLALEstimatePulsarAmplitudeParams ( PulsarCandidate *pulsarParams,	///< [in,out] Pulsar candidate parameters.
                                    const LIGOTimeGPS* FaFb_refTime,	///< [in] Reference time of \f$F_a\f$ and \f$F_b\f$, may differ from pulsar candidate reference time.
                                    const COMPLEX8 Fa,			///< [in] Complex \f$\mathcal{F}\f$-statistic amplitude \f$F_a\f$.
                                    const COMPLEX8 Fb,			///< [in] Complex \f$\mathcal{F}\f$-statistic amplitude \f$F_b\f$.
                                    const AntennaPatternMatrix *Mmunu	///< [in] Antenna pattern matrix \f$M_{\mu\nu}\f$.
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

  XLAL_CHECK ( pulsarParams != NULL, XLAL_EINVAL );
  XLAL_CHECK ( FaFb_refTime != NULL, XLAL_EINVAL );
  XLAL_CHECK ( isfinite(creal(Fa)) && isfinite(cimag(Fb)), XLAL_EDOM,
               "Fa = (%g, %g) is not finite", creal(Fa), cimag(Fa) );
  XLAL_CHECK ( isfinite(creal(Fb)) && isfinite(cimag(Fb)), XLAL_EDOM,
               "Fb = (%g, %g) is not finite", creal(Fb), cimag(Fb) );
  XLAL_CHECK ( Mmunu != NULL, XLAL_EINVAL );

  Ad = Mmunu->Ad;
  Bd = Mmunu->Bd;
  Cd = Mmunu->Cd;
  Ed = Mmunu->Ed;
  Dd = Ad * Bd - Cd * Cd - Ed * Ed;

  normAmu = 2.0 / sqrt(2.0 * Mmunu->Sinv_Tsft); // generally *very* small!!

  // ----- GSL memory allocation -----
  XLAL_CHECK ( ( x_mu = gsl_vector_calloc (4) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( A_Mu = gsl_vector_calloc (4) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( M_Mu_Nu = gsl_matrix_calloc (4, 4) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( Jh_Mu_nu = gsl_matrix_calloc (4, 4) ) != NULL, XLAL_ENOMEM );

  XLAL_CHECK ( ( permh = gsl_permutation_calloc ( 4 ) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( tmp = gsl_matrix_calloc (4, 4) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( tmp2 = gsl_matrix_calloc (4, 4) ) != NULL, XLAL_ENOMEM );

  // ----- fill vector x_mu
  gsl_vector_set (x_mu, 0,   creal(Fa) );       // x_1
  gsl_vector_set (x_mu, 1,   creal(Fb) );       // x_2
  gsl_vector_set (x_mu, 2, - cimag(Fa) );       // x_3
  gsl_vector_set (x_mu, 3, - cimag(Fb) );       // x_4

  // ----- fill matrix M^{mu,nu} [symmetric: use UPPER HALF ONLY!!]
  gsl_matrix_set (M_Mu_Nu, 0, 0,   Bd / Dd );
  gsl_matrix_set (M_Mu_Nu, 1, 1,   Ad / Dd );
  gsl_matrix_set (M_Mu_Nu, 0, 1, - Cd / Dd );

  gsl_matrix_set (M_Mu_Nu, 0, 3, - Ed / Dd );
  gsl_matrix_set (M_Mu_Nu, 1, 2,   Ed / Dd );

  gsl_matrix_set (M_Mu_Nu, 2, 2,   Bd / Dd );
  gsl_matrix_set (M_Mu_Nu, 3, 3,   Ad / Dd );
  gsl_matrix_set (M_Mu_Nu, 2, 3, - Cd / Dd );

  // get (un-normalized) MLE's for amplitudes A^mu  = M^{mu,nu} x_nu

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
  aPlus = sqrt(Ap2);            // not yet normalized

  Ac2 = 0.5 * ( Asq - disc );
  aCross = sqrt( Ac2 );
  aCross *= MYSIGN ( Da );      // not yet normalized

  beta = aCross / aPlus;

  b1 =   A4h - beta * A1h;
  b2 =   A3h + beta * A2h;
  b3 = - A1h + beta * A4h ;

  psi  = 0.5 * atan ( b1 /  b2 );       // in [-pi/4,pi/4] (gauge used also by TDS)
  phi0 =       atan ( b2 / b3 );        // in [-pi/2,pi/2]

  // Fix remaining sign-ambiguity by checking sign of reconstructed A1
  A1check = aPlus * cos(phi0) * cos(2.0*psi) - aCross * sin(phi0) * sin(2*psi);
  if ( A1check * A1h <  0 )
    phi0 += LAL_PI;

  cosphi0 = cos(phi0);
  sinphi0 = sin(phi0);
  cos2psi = cos(2*psi);
  sin2psi = sin(2*psi);

  // check numerical consistency of estimated Amu and reconstructed
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

  // translate A_{+,x} into {h_0, cosi}
  h0 = aPlus + sqrt ( disc );  // not yet normalized !
  cosi = aCross / h0;

  // ========== Estimate the errors ==========

  // ----- compute derivatives \partial A^\mu / \partial B^\nu, where
  // we consider the output-variables B^\nu = (h0, cosi, phi0, psi)
  // where aPlus = 0.5 * h0 * (1 + cosi^2)  and aCross = h0 * cosi
  { // Ahat^mu is defined as A^mu with the replacements: A_+ --> A_x, and A_x --> h0
    REAL8 A1hat =   aCross * cosphi0 * cos2psi - h0 * sinphi0 * sin2psi;
    REAL8 A2hat =   aCross * cosphi0 * sin2psi + h0 * sinphi0 * cos2psi;
    REAL8 A3hat = - aCross * sinphi0 * cos2psi - h0 * cosphi0 * sin2psi;
    REAL8 A4hat = - aCross * sinphi0 * sin2psi + h0 * cosphi0 * cos2psi;

    // ----- A1 =   aPlus * cosphi0 * cos2psi - aCross * sinphi0 * sin2psi; -----
    gsl_matrix_set (Jh_Mu_nu, 0, 0,   A1h / h0 );       /* dA1/h0 */
    gsl_matrix_set (Jh_Mu_nu, 0, 1,   A1hat );          /* dA1/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 0, 2,   A3h );            /* dA1/dphi0 */
    gsl_matrix_set (Jh_Mu_nu, 0, 3, - 2.0 * A2h );      /* dA1/dpsi */

    // ----- A2 =   aPlus * cosphi0 * sin2psi + aCross * sinphi0 * cos2psi; -----
    gsl_matrix_set (Jh_Mu_nu, 1, 0,   A2h / h0 );       /* dA2/h0 */
    gsl_matrix_set (Jh_Mu_nu, 1, 1,   A2hat );          /* dA2/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 1, 2,   A4h );            /* dA2/dphi0 */
    gsl_matrix_set (Jh_Mu_nu, 1, 3,   2.0 * A1h );      /* dA2/dpsi */

    // ----- A3 = - aPlus * sinphi0 * cos2psi - aCross * cosphi0 * sin2psi; -----
    gsl_matrix_set (Jh_Mu_nu, 2, 0,   A3h / h0 );       /* dA3/h0 */
    gsl_matrix_set (Jh_Mu_nu, 2, 1,   A3hat );          /* dA3/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 2, 2, - A1h );            /* dA3/dphi0 */
    gsl_matrix_set (Jh_Mu_nu, 2, 3, - 2.0 * A4h );      /* dA3/dpsi */

    // ----- A4 = - aPlus * sinphi0 * sin2psi + aCross * cosphi0 * cos2psi; -----
    gsl_matrix_set (Jh_Mu_nu, 3, 0,   A4h / h0 );       /* dA4/h0 */
    gsl_matrix_set (Jh_Mu_nu, 3, 1,   A4hat );          /* dA4/dcosi */
    gsl_matrix_set (Jh_Mu_nu, 3, 2, - A2h );            /* dA4/dphi0 */
    gsl_matrix_set (Jh_Mu_nu, 3, 3,   2.0 * A3h );      /* dA4/dpsi */
  }

  // ----- compute inverse matrices Jh^{-1} by LU-decomposition -----
  XLAL_CHECK( gsl_linalg_LU_decomp (Jh_Mu_nu, permh, &signum ) == 0, XLAL_EERR );

  // inverse matrix
  XLAL_CHECK( gsl_linalg_LU_invert (Jh_Mu_nu, permh, tmp ) == 0, XLAL_EERR );
  gsl_matrix_memcpy ( Jh_Mu_nu, tmp );

  // ----- compute Jh^-1 . Minv . (Jh^-1)^T -----

  /* GSL-doc: gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha,
   *                          const gsl_matrix *A, const gsl_matrix *B, double beta, gsl_matrix *C)
   * These functions compute the matrix-matrix product and sum
   * C = \alpha op(A) op(B) + \beta C
   * where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans
   * and similarly for the parameter TransB.
   */

  // first tmp = Minv . (Jh^-1)^T
  XLAL_CHECK( gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, M_Mu_Nu, Jh_Mu_nu, 0.0, tmp ) == 0, XLAL_EERR );
  // then J^-1 . tmp , store result in tmp2
  XLAL_CHECK( gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Jh_Mu_nu, tmp, 0.0, tmp2 ) == 0, XLAL_EERR );
  gsl_matrix_memcpy ( Jh_Mu_nu, tmp2 );

  // ===== debug-output resulting matrices =====
  // propagate initial-phase from Fstat-reference-time to refTime of Doppler-params
  // XLALExtrapolatePulsarPhase() guarantees propagated phi0 is in [0, 2*pi]
  const REAL8 dtau = XLALGPSDiff( &pulsarParams->Doppler.refTime, FaFb_refTime );
  XLAL_CHECK( XLALExtrapolatePulsarPhase( &phi0, pulsarParams->Doppler.fkdot, phi0, dtau ) == XLAL_SUCCESS, XLAL_EFUNC );

  // fill candidate-struct with the obtained signal-parameters and error-estimations
  pulsarParams->Amp.h0     = normAmu * h0;
  pulsarParams->Amp.cosi   = cosi;
  pulsarParams->Amp.phi0   = phi0;
  pulsarParams->Amp.psi    = psi;

  // read out principal estimation-errors from diagonal elements of inverse Fisher-matrix
  pulsarParams->dAmp.h0     = normAmu * sqrt( gsl_matrix_get (Jh_Mu_nu, 0, 0 ) );
  pulsarParams->dAmp.cosi   = sqrt( gsl_matrix_get (Jh_Mu_nu, 1, 1 ) );
  pulsarParams->dAmp.phi0   = sqrt( gsl_matrix_get (Jh_Mu_nu, 2, 2 ) );
  pulsarParams->dAmp.psi    = sqrt( gsl_matrix_get (Jh_Mu_nu, 3, 3 ) );
  // return also the full Amplitude-Fisher matrix: invert Jh_Mu_nu
  XLAL_CHECK( gsl_linalg_LU_decomp (Jh_Mu_nu, permh, &signum ) == 0, XLAL_EERR );
  XLAL_CHECK( gsl_linalg_LU_invert (Jh_Mu_nu, permh, tmp ) == 0, XLAL_EERR );
  pulsarParams->AmpFisherMatrix = tmp;

  // ----- free GSL memory -----
  gsl_vector_free ( x_mu );
  gsl_vector_free ( A_Mu );
  gsl_matrix_free ( M_Mu_Nu );
  gsl_matrix_free ( Jh_Mu_nu );
  gsl_permutation_free ( permh );
  gsl_matrix_free ( tmp2 );

  return XLAL_SUCCESS;

} // XLALEstimatePulsarAmplitudeParams()

///
/// Convert amplitude params from 'physical' coordinates \f$(h_0, \cos\iota, \psi, \phi_0)\f$ into
/// 'canonical' coordinates \f$A^\mu = (A_1, A_2, A_3, A_4)\f$. The equations can be found in
/// \cite JKS98 or \cite Prix07 Eq.(2).
///
int
XLALAmplitudeParams2Vect ( PulsarAmplitudeVect A_Mu,		///< [out] Canonical amplitude coordinates \f$A^\mu = (A_1, A_2, A_3, A_4)\f$.
                           const PulsarAmplitudeParams Amp	///< [in] Physical amplitude params \f$(h_0, \cos\iota, \psi, \phi_0)\f$.
                           )
{

  REAL8 aPlus = 0.5 * Amp.h0 * ( 1.0 + SQ(Amp.cosi) );
  REAL8 aCross = Amp.h0 * Amp.cosi;
  REAL8 cos2psi = cos ( 2.0 * Amp.psi );
  REAL8 sin2psi = sin ( 2.0 * Amp.psi );
  REAL8 cosphi0 = cos ( Amp.phi0 );
  REAL8 sinphi0 = sin ( Amp.phi0 );

  XLAL_CHECK( A_Mu != NULL, XLAL_EINVAL );

  A_Mu[0] =  aPlus * cos2psi * cosphi0 - aCross * sin2psi * sinphi0;
  A_Mu[1] =  aPlus * sin2psi * cosphi0 + aCross * cos2psi * sinphi0;
  A_Mu[2] = -aPlus * cos2psi * sinphi0 - aCross * sin2psi * cosphi0;
  A_Mu[3] = -aPlus * sin2psi * sinphi0 + aCross * cos2psi * cosphi0;

  return XLAL_SUCCESS;

} // XLALAmplitudeParams2Vect()

///
/// Compute amplitude params \f$(h_0, \cos\iota, \psi, \phi_0)\f$ from amplitude-vector \f$A^\mu = (A_1, A_2, A_3, A_4)\f$.
/// Adapted from algorithm in XLALEstimatePulsarAmplitudeParams().
///
int
XLALAmplitudeVect2Params ( PulsarAmplitudeParams *Amp,		///< [out] Physical amplitude params \f$(h_0, \cos\iota, \psi, \phi_0)\f$.
                           const PulsarAmplitudeVect A_Mu	///< [in] Canonical amplitude coordinates \f$A^\mu = (A_1, A_2, A_3, A_4)\f$.
                           )
{

  REAL8 h0Ret, cosiRet, psiRet, phi0Ret;

  REAL8 A1, A2, A3, A4, Asq, Da, disc;
  REAL8 Ap2, Ac2, aPlus, aCross;
  REAL8 beta, b1, b2, b3;

  XLAL_CHECK( A_Mu != NULL, XLAL_EINVAL );
  XLAL_CHECK( Amp != NULL, XLAL_EINVAL );

  A1 = A_Mu[0];
  A2 = A_Mu[1];
  A3 = A_Mu[2];
  A4 = A_Mu[3];

  Asq = SQ(A1) + SQ(A2) + SQ(A3) + SQ(A4);
  Da = A1 * A4 - A2 * A3;

  disc = sqrt ( SQ(Asq) - 4.0 * SQ(Da) );

  Ap2  = 0.5 * ( Asq + disc );
  aPlus = sqrt(Ap2);

  Ac2 = 0.5 * ( Asq - disc );
  aCross = MYSIGN(Da) * sqrt( Ac2 );

  beta = aCross / aPlus;

  b1 =   A4 - beta * A1;
  b2 =   A3 + beta * A2;
  b3 = - A1 + beta * A4 ;

  // amplitude params in LIGO conventions
  psiRet  = 0.5 * atan2 ( b1,  b2 );  /* [-pi/2,pi/2] */
  phi0Ret =       atan2 ( b2,  b3 );  /* [-pi, pi] */

  // Fix remaining sign-ambiguity by checking sign of reconstructed A1
  REAL8 A1check = aPlus * cos(phi0Ret) * cos(2.0*psiRet) - aCross * sin(phi0Ret) * sin(2*psiRet);
  if ( A1check * A1 < 0 ) {
    phi0Ret += LAL_PI;
  }

  h0Ret = aPlus + sqrt ( disc );
  cosiRet = aCross / h0Ret;

  // make unique by fixing the gauge to be psi in [-pi/4, pi/4], phi0 in [0, 2*pi]
  while ( psiRet > LAL_PI_4 ) {
    psiRet  -= LAL_PI_2;
    phi0Ret -= LAL_PI;
  }
  while ( psiRet < - LAL_PI_4 ) {
    psiRet  += LAL_PI_2;
    phi0Ret += LAL_PI;
  }
  while ( phi0Ret < 0 ) {
    phi0Ret += LAL_TWOPI;
  }

  while ( phi0Ret > LAL_TWOPI ) {
    phi0Ret -= LAL_TWOPI;
  }

  // Return final answer
  Amp->h0   = h0Ret;
  Amp->cosi = cosiRet;
  Amp->psi  = psiRet;
  Amp->phi0 = phi0Ret;

  return XLAL_SUCCESS;

} // XLALAmplitudeVect2Params()
