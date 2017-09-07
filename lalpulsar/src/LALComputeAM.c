/*
*  Copyright (C) 2007 John Whelan, Reinhard Prix
*  Copyright (C) 2007 Jolien Creighton, Maria Alessandra Papa, Steve Berukoff, Xavier Siemens
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
 * \author S.J. Berukoff, Reinhard Prix, John Whelan
 * \date 2007
 * \addtogroup LALComputeAM_h
 * \brief Computes quantities for amplitude demodulation.
 *
 * This routine computes the quantities \f$a(t)\f$ and \f$b(t)\f$ as defined in
 * Jaranowski, Krolak, and Schutz \cite JKS98 , hereafter JKS.  These
 * functions quantify the dependence of the detector output on the
 * beam-pattern functions \f$F_{+}\f$ and \f$F_{\times}\f$; in fact, \f$a(t)\f$ and
 * \f$b(t)\f$ <i>are</i> the beam-pattern functions, without the dependence
 * on polarization angle and detector arm angle.  Since the
 * LALDemod() suite is an attempt to compute an optimal statistic,
 * it is necessary to include these quantities in the computation.
 * Otherwise, the motion of the Earth as it revolves about its axis will
 * smear the signal into several neighboring bins centered about the
 * search frequency, consequently losing valuable SNR.
 *
 * ### Algorithm ###
 *
 * The routine is really simple.  From JKS,
 * \f{eqnarray*}
 * F_{+} &=& \sin\zeta [ a(t) \cos 2 \psi + b(t) \sin 2 \psi ] \\
 * F_{\times} &=& \sin\zeta [ b(t) \cos 2 \psi - a(t) \sin 2 \psi ]
 * \f}
 * We use the routine LALComputeDetAMResponse() to calculate
 * \f$F_{+}\f$ and \f$F_{\times}\f$ for a given polarization angle, and then
 * extract \f$a(t)\f$ and \f$b(t)\f$, once for each timestamp \f$t\f$.  Additionally,
 * computation of the optimal statistic requires that we compute inner
 * products of these two quantities for later use.
 *
 */

/*---------- INCLUDES ----------*/
#include <lal/LALComputeAM.h>
#include <lal/SinCosLUT.h>

/*---------- local DEFINES and macros ----------*/

#define SQ(x) (x) * (x)

/*---------- internal types ----------*/

/*---------- internal prototypes ----------*/

/*==================== FUNCTION DEFINITIONS ====================*/

/**
 * Compute single time-stamp antenna-pattern coefficients a(t), b(t)
 * Note: this function uses REAL8 precision, so this can be used in
 * high-precision integration of the F-metric
 *
 */
int
XLALComputeAntennaPatternCoeffs ( REAL8 *ai,   			/**< [out] antenna-pattern function a(t) */
				  REAL8 *bi,			/**< [out] antenna-pattern function b(t) */
				  const SkyPosition *skypos,	/**< [in] skyposition {alpha, delta} */
				  const LIGOTimeGPS *tGPS,	/**< [in] GPS time t */
				  const LALDetector *site,	/**< [in] detector */
				  const EphemerisData *edat	/**< [in] ephemeris-data */
				  )
{
  LALStatus XLAL_INIT_DECL(status);
  EarthState XLAL_INIT_DECL(earth);

  if ( !ai || !bi || !skypos || !tGPS || !site || !edat) {
    XLAL_ERROR( XLAL_EINVAL );
  }

  XLAL_CHECK( XLALBarycenterEarth (&earth, tGPS, edat ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* ---------- compute the detector tensor ---------- */
  REAL8 sinG, cosG, sinGcosG, sinGsinG, cosGcosG;
  SymmTensor3d detT;

  sinG = sin ( earth.gmstRad );
  cosG = cos ( earth.gmstRad );

  sinGsinG = sinG * sinG;
  sinGcosG = sinG * cosG;
  cosGcosG = cosG * cosG;

  detT.d11 = site->response[0][0] * cosGcosG
        - 2 * site->response[0][1] * sinGcosG
            + site->response[1][1] * sinGsinG;
  detT.d22 = site->response[0][0] * sinGsinG
        + 2 * site->response[0][1] * sinGcosG
            + site->response[1][1] * cosGcosG;
  detT.d12 = (site->response[0][0] - site->response[1][1]) * sinGcosG
        + site->response[0][1] * (cosGcosG - sinGsinG);
  detT.d13 = site->response[0][2] * cosG
            - site->response[1][2] * sinG;
  detT.d23 = site->response[0][2] * sinG
            + site->response[1][2] * cosG;
  detT.d33 = site->response[2][2];


  /*---------- We write components of xi and eta vectors in SSB-fixed coords */
  REAL8 delta, alpha;
  REAL8 sin1delta, cos1delta;
  REAL8 sin1alpha, cos1alpha;

  REAL8 xi1, xi2;
  REAL8 eta1, eta2, eta3;


  alpha = skypos->longitude;
  delta = skypos->latitude;

  sin1delta = sin(delta);
  cos1delta = cos(delta);

  sin1alpha = sin(alpha);
  cos1alpha = cos(alpha);

  xi1 = - sin1alpha;
  xi2 =  cos1alpha;
  eta1 = sin1delta * cos1alpha;
  eta2 = sin1delta * sin1alpha;
  eta3 = - cos1delta;

  /*---------- Compute the a(t_i) and b(t_i) ---------- */
  (*ai) = detT.d11 * ( xi1 * xi1 - eta1 * eta1 )
    + 2 * detT.d12 * ( xi1*xi2 - eta1*eta2 )
    - 2 * detT.d13 *             eta1 * eta3
    +     detT.d22 * ( xi2*xi2 - eta2*eta2 )
    - 2 * detT.d23 *             eta2 * eta3
    -     detT.d33 *             eta3*eta3;

  (*bi) = detT.d11 * 2 * xi1 * eta1
    + 2 * detT.d12 *   ( xi1 * eta2 + xi2 * eta1 )
    + 2 * detT.d13 *     xi1 * eta3
    +     detT.d22 * 2 * xi2 * eta2
    + 2 * detT.d23 *     xi2 * eta3;


  return XLAL_SUCCESS;

} /* XLALComputeAntennaPatternCoeffs() */


/**
 * <b>Replace</b> AM-coeffs by weighted AM-coeffs, i.e.
 * \f$a_{X\alpha} \rightarrow \widehat{a}_{X\alpha} \equiv a_{X\alpha} \sqrt{w_{X\alpha}}\f$, and
 * \f$b_{X\alpha} \rightarrow \widehat{b}_{X\alpha} \equiv a_{X\alpha} \sqrt{w_{X\alpha}}\f$,
 * where \f$w_{X\alpha}\f$ are the \a multiWeights for SFT \f$\alpha\f$ and detector \f$X\f$.
 *
 * Also compute the resulting per-detector \f$X\f$ antenna-pattern matrix coefficients
 * \f$\widehat{A}_X \equiv \sum_{\alpha} \widehat{a}^2_{X\alpha}\f$,
 * \f$\widehat{B}_X \equiv \sum_{\alpha} \widehat{b}^2_{X\alpha}\f$,
 * \f$\widehat{C}_X \equiv \sum_{\alpha} \widehat{a}_{X\alpha}\widehat{b}_{X\alpha}\f$,
 *
 * and corresponding multi-detector antenna-pattern matrix coefficients
 * \f$\widehat{A} \equiv \sum_{X} \widehat{A}_X\f$,
 * \f$\widehat{B} \equiv \sum_{X} \widehat{B}_X\f$,
 * \f$\widehat{C} \equiv \sum_{X} \widehat{C}_X\f$.
 *
 * See Sec.4.1 in CFSv2.pdf notes https://dcc.ligo.org/cgi-bin/DocDB/ShowDocument?docid=T0900149&version=4
 *
 * \note *) this function modifies the input AMCoeffs->{a,b,c} *in place* !
 * \note *) if multiWeights = NULL, we assume unit-weights.
 */
int
XLALWeightMultiAMCoeffs (  MultiAMCoeffs *multiAMcoef, const MultiNoiseWeights *multiWeights )
{

  /* ----- input sanity checks ----- */
  if ( !multiAMcoef ) {
    XLALPrintError ("%s: illegal NULL input received in 'multiAMcoefs'.\n", __func__ );
    XLAL_ERROR( XLAL_EINVAL );
  }
  UINT4 numDetectors = multiAMcoef->length;
  /* make sure identical number of detectors in amCoefs and weights */
  if ( multiWeights && (multiWeights->length != numDetectors) ) {
    XLALPrintError("%s: multiWeights must be NULL or have the same number of detectors (numDet=%d) as mulitAMcoef (numDet=%d)!\n", __func__, multiWeights->length, numDetectors );
    XLAL_ERROR( XLAL_EINVAL );
  }
  /* make sure identical number of timesteps in a_X, b_X and (if given) weights w_X, respectively */
  UINT4 X;
  for ( X = 0; X < numDetectors; X ++ )
    {
      UINT4 numStepsX = multiAMcoef->data[X]->a->length;
      if ( numStepsX != multiAMcoef->data[X]->b->length ) {
        XLALPrintError ("%s: per-SFT antenna-pattern series have different length: a_alpha (len=%d), b_alpha (len=%d)\n", __func__, numStepsX, multiAMcoef->data[X]->b->length );
        XLAL_ERROR ( XLAL_EINVAL );
      }
      if ( multiWeights && (multiWeights->data[X]->length != numStepsX )) {
        XLALPrintError("%s: multiWeights[X=%d] must be NULL or have the same length (len=%d) as mulitAMcoef[X] (len=%d)!\n", __func__, X, multiWeights->data[X]->length, numStepsX );
        XLAL_ERROR( XLAL_EINVAL );
      }
    } // for X < numDetectors

  REAL4 Ad = 0, Bd = 0, Cd = 0, Ed = 0;	// multi-IFO values
  /* ---------- main loop over detectors X ---------- */
  for ( X=0; X < numDetectors; X ++)
    {
      AMCoeffs *amcoeX = multiAMcoef->data[X];
      UINT4 numStepsX = amcoeX->a->length;

      /* ----- if given, apply noise-weights to all Antenna-pattern coefficients from detector X ----- */
      if ( multiWeights )
        {
          REAL8Vector *weightsX = multiWeights->data[X];
          UINT4 alpha;	// SFT-index
          REAL8 ooSinv_Tsft = 1.0 / multiWeights->Sinv_Tsft;
          for(alpha = 0; alpha < numStepsX; alpha++)
            {
              REAL8 weight =  multiWeights->isNotNormalized ? weightsX->data[alpha]*ooSinv_Tsft : weightsX->data[alpha];
              REAL8 Sqwi = sqrt ( weight );
              /* apply noise-weights, *replace* original a, b by noise-weighed version! */
	      amcoeX->a->data[alpha] *= Sqwi;
	      amcoeX->b->data[alpha] *= Sqwi;
            } // for alpha < numSteps
        } // if weights

      UINT4 alpha;	// SFT-index
      REAL4 AdX = 0, BdX = 0, CdX = 0, EdX = 0;	// single-IFO values
      /* compute single-IFO antenna-pattern coefficients AX,BX,CX, by summing over time-steps 'alpha' */
      for(alpha = 0; alpha < numStepsX; alpha++)
        {
          REAL4 ahat = amcoeX->a->data[alpha];
          REAL4 bhat = amcoeX->b->data[alpha];

          AdX += ahat * ahat;
          BdX += bhat * bhat;
          CdX += ahat * bhat;
          // EdX = 0; // trivial for real-values a,b
        } /* for alpha < numStepsX */

      /* store those */
      amcoeX->A = AdX;
      amcoeX->B = BdX;
      amcoeX->C = CdX;
      amcoeX->D = AdX * BdX - CdX * CdX - EdX * EdX;

      // in the unlikely event of a degenerate M-matrix with D = sqrt ( det M ) <= 0,
      // we set D->inf, in order to set the corresponding F-value to zero rather than >>1
      // By setting 'D=inf', we also allow upstream catching/filtering on such singular cases
      if ( amcoeX->D <= 0 ) {
	amcoeX->D = INFINITY;
      }
      /* compute multi-IFO antenna-pattern coefficients A,B,C,E by summing over IFOs X */
      Ad += AdX;
      Bd += BdX;
      Cd += CdX;
      // Ed = 0; // trivial for real-valued a,b
    } /* for X < numDetectors */

  multiAMcoef->Mmunu.Ad = Ad;
  multiAMcoef->Mmunu.Bd = Bd;
  multiAMcoef->Mmunu.Cd = Cd;
  multiAMcoef->Mmunu.Dd = Ad * Bd - Cd * Cd - Ed * Ed;

  // in the unlikely event of a degenerate M-matrix with D = sqrt(det M_munu) <= 0,
  // we set D->inf, in order to set the corresponding F-value to zero rather than >>1
  // By setting 'D=inf', we also allow upstream catching/filtering on such singular cases
  if ( multiAMcoef->Mmunu.Dd <= 0 ) {
    multiAMcoef->Mmunu.Dd = INFINITY;
  }


  if ( multiWeights ) {
    multiAMcoef->Mmunu.Sinv_Tsft = multiWeights->Sinv_Tsft;
  }

  return XLAL_SUCCESS;

} /* XLALWeightMultiAMCoeffs() */


/**
 * Compute the 'amplitude coefficients' \f$a(t)\sin\zeta\f$,
 * \f$b(t)\sin\zeta\f$ as defined in \cite JKS98 for a series of
 * timestamps.
 *
 * The input consists of the DetectorState-timeseries, which contains
 * the detector-info and the LMST's corresponding to the different times.
 *
 * \note This implementation is based on the geometrical definition of
 * \f$a\sin\zeta\f$ and \f$b\sin\zeta\f$ as detector response
 * coefficients in a preferred polarization basis.  (It is thereby
 * more general than the JKS expressions and could be used e.g., with
 * the response tensor of a bar detector with no further modification
 * needed.)
 *
 * \note The fields AMCoeffs->{A, B, C, D} are not computed by this function,
 * as they require correct SFT noise-weights. These fields would be computed, for
 * example by XLALWeightMultiAMCoeffs().
 *
 */
AMCoeffs *
XLALComputeAMCoeffs ( const DetectorStateSeries *DetectorStates,	/**< timeseries of detector states */
                      SkyPosition skypos				/**< {alpha,delta} of the source */
                      )
{
  /* ---------- check input consistency ---------- */
  if ( !DetectorStates ) {
    XLALPrintError ("%s: invalid NULL input 'DetectorStates'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /* currently requires sky-pos to be in equatorial coordinates (FIXME) */
  if ( skypos.system != COORDINATESYSTEM_EQUATORIAL ) {
    XLALPrintError ("%s: only equatorial coordinates currently supported in 'skypos'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /*---------- We write components of xi and eta vectors in SSB-fixed coords */
  REAL4 alpha = skypos.longitude;
  REAL4 delta = skypos.latitude;

  REAL4 sin1delta, cos1delta;
  REAL4 sin1alpha, cos1alpha;
  XLAL_CHECK_NULL( XLALSinCosLUT (&sin1delta, &cos1delta, delta ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_NULL( XLALSinCosLUT (&sin1alpha, &cos1alpha, alpha ) == XLAL_SUCCESS, XLAL_EFUNC );

  REAL4 xi1 = - sin1alpha;
  REAL4 xi2 =  cos1alpha;
  REAL4 eta1 = sin1delta * cos1alpha;
  REAL4 eta2 = sin1delta * sin1alpha;
  REAL4 eta3 = - cos1delta;

  /* prepare output vector */
  UINT4 numSteps = DetectorStates->length;
  AMCoeffs *coeffs;
  if ( ( coeffs = XLALCreateAMCoeffs ( numSteps ) ) == NULL ) {
    XLALPrintError ("%s: XLALCreateAMCoeffs(%d) failed\n", __func__, numSteps );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  /*---------- Compute the a(t_i) and b(t_i) ---------- */
  UINT4 i;
  for ( i=0; i < numSteps; i++ )
    {
      REAL4 ai, bi;

      SymmTensor3 *d = &(DetectorStates->data[i].detT);

      ai =    d->d11 * ( xi1 * xi1 - eta1 * eta1 )
	+ 2 * d->d12 * ( xi1*xi2 - eta1*eta2 )
	- 2 * d->d13 *             eta1 * eta3
	+     d->d22 * ( xi2*xi2 - eta2*eta2 )
	- 2 * d->d23 *             eta2 * eta3
	-     d->d33 *             eta3*eta3;

      bi =    d->d11 * 2 * xi1 * eta1
	+ 2 * d->d12 *   ( xi1 * eta2 + xi2 * eta1 )
	+ 2 * d->d13 *     xi1 * eta3
	+     d->d22 * 2 * xi2 * eta2
	+ 2 * d->d23 *     xi2 * eta3;

      coeffs->a->data[i] = ai;
      coeffs->b->data[i] = bi;

    } /* for i < numSteps */

  /* return the result */
  return coeffs;

} /* XLALComputeAMCoeffs() */

/**
 * Multi-IFO version of XLALComputeAMCoeffs().
 * Computes noise-weighted combined multi-IFO antenna pattern functions.
 *
 * This function applies
 * the noise-weights and computes  the multi-IFO antenna-pattern matrix components
 * {A, B, C}, and single-IFO matrix components {A_X,B_X,C_X} for detector X.
 *
 * Therefore: DONT use XLALWeightMultiAMCoeffs() on the result!
 *
 * \note *) an input of multiWeights = NULL corresponds to unit-weights
 */
MultiAMCoeffs *
XLALComputeMultiAMCoeffs ( const MultiDetectorStateSeries *multiDetStates, 	/**< [in] detector-states at timestamps t_i */
                           const MultiNoiseWeights *multiWeights,		/**< [in] noise-weigths at timestamps t_i (can be NULL) */
                           SkyPosition skypos					/**< source sky-position [in equatorial coords!] */
                           )
{
  /* check input consistency */
  if ( !multiDetStates ) {
    XLALPrintError ("%s: invalid NULL input argument 'multiDetStates'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  UINT4 numDetectors = multiDetStates->length;

  /* prepare output vector */
  MultiAMCoeffs *ret;
  if ( ( ret = XLALCalloc( 1, sizeof( *ret ) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc( 1, %zu)\n", __func__, sizeof( *ret ) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  ret->length = numDetectors;
  if ( ( ret->data = XLALCalloc ( numDetectors, sizeof ( *ret->data ) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc(%d, %zu)\n", __func__, numDetectors, sizeof ( *ret->data ) );
    XLALFree ( ret );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  /* loop over detectors and generate AMCoeffs for each one */
  UINT4 X;
  for ( X=0; X < numDetectors; X ++ )
    {
      if ( (ret->data[X] = XLALComputeAMCoeffs ( multiDetStates->data[X], skypos )) == NULL ) {
        XLALPrintError ("%s: call to XLALComputeAMCoeffs() failed with xlalErrno = %d\n", __func__, xlalErrno );
        XLALDestroyMultiAMCoeffs ( ret );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }

    } /* for X < numDetectors */

  /* apply noise-weights and compute antenna-pattern matrix {A,B,C} */
  if ( XLALWeightMultiAMCoeffs (  ret, multiWeights ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: call to XLALWeightMultiAMCoeffs() failed with xlalErrno = %d\n", __func__, xlalErrno );
    XLALDestroyMultiAMCoeffs ( ret );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  /* return result */
  return ret;

} /* XLALComputeMultiAMCoeffs() */


/* ---------- creators/destructors for AM-coeffs -------------------- */
/**
 * Create an AMCeoffs vector for given number of timesteps
 */
AMCoeffs *
XLALCreateAMCoeffs ( UINT4 numSteps )
{
  AMCoeffs *ret;

  if ( ( ret = XLALCalloc ( 1, sizeof (*ret) ) ) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc ( 1, %zu )\n", __func__, sizeof (*ret) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  if ( ( ret->a = XLALCreateREAL4Vector ( numSteps )) == NULL ) {
    XLALPrintError ("%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numSteps );
    XLALDestroyAMCoeffs ( ret );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  if ( ( ret->b = XLALCreateREAL4Vector ( numSteps )) == NULL ) {
    XLALPrintError ("%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numSteps );
    XLALDestroyAMCoeffs ( ret );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  return ret;

} /* XLALCreateAMCoeffs() */


/**
 * Destroy a MultiAMCoeffs structure.
 *
 * Note, this is "NULL-robust" in the sense that it will not crash
 * on NULL-entries anywhere in this struct, so it can be used
 * for failure-cleanup even on incomplete structs
 */
void
XLALDestroyMultiAMCoeffs ( MultiAMCoeffs *multiAMcoef )
{
  UINT4 X;

  if ( ! multiAMcoef )
    return;

  if ( multiAMcoef->data )
    {
      for ( X=0; X < multiAMcoef->length; X ++ )
	{
	  XLALDestroyAMCoeffs ( multiAMcoef->data[X] );
	} /* for X < numDetectors */
      LALFree ( multiAMcoef->data );
    }
  LALFree ( multiAMcoef );

  return;

} /* XLALDestroyMultiAMCoeffs() */

/**
 * Destroy a AMCoeffs structure.
 *
 * \note This function is "NULL-robust" in the sense that it will not crash
 * on NULL-entries anywhere in this struct, so it can be used
 * for failure-cleanup even on incomplete structs
 */
void
XLALDestroyAMCoeffs ( AMCoeffs *amcoef )
{
  if ( ! amcoef )
    return;

  if ( amcoef->a )
    XLALDestroyREAL4Vector ( amcoef->a );
  if ( amcoef->b )
    XLALDestroyREAL4Vector ( amcoef->b );

  LALFree ( amcoef );

  return;

} /* XLALDestroyAMCoeffs() */

