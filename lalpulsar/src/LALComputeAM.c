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
 * Jaranowski, Krolak, and Schutz \ref JKS98, hereafter JKS.  These
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
#include <lal/ComputeFstat.h>

/*---------- local DEFINES and macros ----------*/

#define SQ(x) (x) * (x)

/*---------- internal types ----------*/

/*---------- Global variables ----------*/
/* empty initializers  */
const AMCoeffs empty_AMCoeffs;
const MultiAMCoeffs empty_MultiAMCoeffs;
const AntennaPatternMatrix empty_AntennaPatternMatrix;

static const LALStatus empty_LALStatus;
static const EarthState empty_EarthState;

/*---------- internal prototypes ----------*/


/*==================== FUNCTION DEFINITIONS ====================*/

/**
 * Compute the 'amplitude coefficients' \f$a(t), b(t)\f$ as defined in
 * \ref JKS98 for a series of timestamps.
 *
 * The input consists of the DetectorState-timeseries, which contains
 * the detector-info and the LMST's corresponding to the different times.
 *
 * In order to allow re-using the output-structure AMCoeffs for subsequent
 * calls, we require the REAL4Vectors a and b to be allocated already and
 * to have the same length as the DetectoStates-timeseries.
 *
 * \note This is an alternative implementation to LALComputeAM() with
 * the aim to be both simpler and faster.
 * The difference being that we don't implicitly re-derive the final expression
 * here but simply try to implement the final expressions (12), (13) in \ref JKS98
 * in the most economical way possible.
 */
void
LALGetAMCoeffs(LALStatus *status,				/**< [in/out] LAL status structure pointer */
	       AMCoeffs *coeffs,				/**< [out] amplitude-coeffs {a(t_i), b(t_i)} */
	       const DetectorStateSeries *DetectorStates,	/**< timeseries of detector states */
	       SkyPosition skypos				/**< {alpha,delta} of the source */
	       )
{
  REAL4 ah1, ah2, ah3, ah4, ah5;
  REAL4 a1, a2, a3, a4, a5;

  REAL4 bh1, bh2, bh3, bh4;
  REAL4 b1, b2, b3, b4;

  REAL4 delta, alpha;
  REAL4 sin1delta, cos1delta, sin2delta, cos2delta;

  REAL4 gam, lambda;
  REAL4 norm;
  UINT4 i, numSteps;

  INITSTATUS(status);

  /*---------- check input ---------- */
  ASSERT ( DetectorStates, status, LALCOMPUTEAMH_ENULL, LALCOMPUTEAMH_MSGENULL);

  numSteps = DetectorStates->length;

  /* require the coeffients-vectors to be allocated and consistent with timestamps */
  ASSERT ( coeffs, status, LALCOMPUTEAMH_ENULL, LALCOMPUTEAMH_MSGENULL);
  ASSERT ( coeffs->a && coeffs->b, status, LALCOMPUTEAMH_ENULL, LALCOMPUTEAMH_MSGENULL);
  ASSERT ( (coeffs->a->length == numSteps) && (coeffs->b->length == numSteps), status,
	   LALCOMPUTEAMH_EINPUT,  LALCOMPUTEAMH_MSGEINPUT);

  /* require sky-pos to be in equatorial coordinates */
  ASSERT ( skypos.system == COORDINATESYSTEM_EQUATORIAL, status,
	   SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );

  /*---------- detector paramters: lambda, L, gamma */
  {
    /* FIXME: put into DetectorStateSeries */
    /* orientation of detector arms */
    REAL8 xAzi = DetectorStates->detector.frDetector.xArmAzimuthRadians;
    REAL8 yAzi = DetectorStates->detector.frDetector.yArmAzimuthRadians;

    /* get detector orientation gamma */
    gam = LAL_PI_2 - 0.5 * (xAzi + yAzi);
    /* get detector position latitude (lambda) */
    lambda = DetectorStates->detector.frDetector.vertexLatitudeRadians;
    /*
    printf ("IFO = %s: sin(zeta) = %f\n", DetectorStates->detector.frDetector.name, sin( xAzi - yAzi ) );
    */
  }

  /*---------- coefficient ahN, bhN dependent ONLY on detector-position  ---------- */
  /* FIXME: put these coefficients into DetectorStateSeries */
  {
    REAL4 sin2gamma, cos2gamma;
    REAL4 sin1lambda, cos1lambda;
    REAL4 sin2lambda, cos2lambda;

    sin_cos_LUT (&sin2gamma, &cos2gamma, 2.0f * gam );
    sin_cos_LUT (&sin1lambda, &cos1lambda, lambda );

    sin2lambda = 2.0f * sin1lambda * cos1lambda;
    cos2lambda = cos1lambda * cos1lambda - sin1lambda * sin1lambda;

    /* coefficients for a(t) */
    ah1 = 0.0625f * sin2gamma * (3.0f - cos2lambda);	/* 1/16 = 0.0625 */
    ah2 = - 0.25f * cos2gamma * sin1lambda;
    ah3 =   0.25f * sin2gamma * sin2lambda;
    ah4 =  -0.5f  * cos2gamma * cos1lambda;
    ah5 =  0.75f  * sin2gamma * cos1lambda * cos1lambda;

    /* coefficients for b(t) */
    bh1 =           cos2gamma * sin1lambda;
    bh2 =   0.25f * sin2gamma * (3.0f - cos2lambda);
    bh3 =           cos2gamma * cos1lambda;
    bh4 =   0.5f  * sin2gamma * sin2lambda;
  }

  /*---------- coefficients aN, bN dependent ONLY on {ahN, bhN} and source-latitude delta */
  alpha = skypos.longitude;
  delta = skypos.latitude;

  sin_cos_LUT (&sin1delta, &cos1delta, delta );
  sin2delta = 2.0f * sin1delta * cos1delta;
  cos2delta = cos1delta * cos1delta - sin1delta * sin1delta;

  /* coefficients for a(t) */
  a1 = ah1 * ( 3.0f - cos2delta );
  a2 = ah2 * ( 3.0f - cos2delta );
  a3 = ah3 * sin2delta;
  a4 = ah4 * sin2delta;
  a5 = ah5 * cos1delta * cos1delta;

  /* coefficients for b(t) */
  b1 = bh1 * sin1delta;
  b2 = bh2 * sin1delta;
  b3 = bh3 * cos1delta;
  b4 = bh4 * cos1delta;


  /*---------- Compute the a(t_i) and b(t_i) ---------- */
  coeffs->A = 0;
  coeffs->B = 0;
  coeffs->C = 0;
  coeffs->D = 0;
  for ( i=0; i < numSteps; i++ )
    {
      REAL4 ah;
      REAL4 cos1ah, sin1ah, cos2ah, sin2ah;
      REAL4 ai, bi;

      ah = alpha - DetectorStates->data[i].LMST;

      sin_cos_LUT ( &sin1ah, &cos1ah, ah );
      sin2ah = 2.0f * sin1ah * cos1ah;
      cos2ah = cos1ah * cos1ah - sin1ah * sin1ah;

      ai = a1 * cos2ah + a2 * sin2ah + a3 * cos1ah + a4 * sin1ah + a5;
      bi = b1 * cos2ah + b2 * sin2ah + b3 * cos1ah + b4 * sin1ah;
      coeffs->a->data[i] = ai;
      coeffs->b->data[i] = bi;

      /* sum A, B, C on the fly */
      coeffs->A += ai * ai;
      coeffs->B += bi * bi;
      coeffs->C += ai * bi;

    } /* for i < numSteps */

  /* finish calculation of A,B,C, D */
  norm = 2.0f / numSteps;
  coeffs->A *= norm;
  coeffs->B *= norm;
  coeffs->C *= norm;

  coeffs->D = coeffs->A * coeffs->B - coeffs->C * coeffs->C;

  RETURN(status);

} /* LALGetAMCoeffs() */

/**
 * Compute the 'amplitude coefficients' \f$a(t)\sin\zeta\f$,
 * \f$b(t)\sin\zeta\f$ as defined in \ref JKS98 for a series of
 * timestamps.
 *
 * The input consists of the DetectorState-timeseries, which contains
 * the detector-info and the LMST's corresponding to the different times.
 *
 * In order to allow re-using the output-structure AMCoeffs for subsequent
 * calls, we require the REAL4Vectors a and b to be allocated already and
 * to have the same length as the DetectoStates-timeseries.
 *
 * \note This is an alternative implementation to both LALComputeAM()
 * and LALGetAMCoeffs(), which uses the geometrical definition of
 * \f$a\sin\zeta\f$ and \f$b\sin\zeta\f$ as detector response
 * coefficients in a preferred polarization basis.  (It is thereby
 * more general than the JKS expressions and could be used e.g., with
 * the response tensor of a bar detector with no further modification
 * needed.)
 */
void
LALNewGetAMCoeffs(LALStatus *status,			/**< [in/out] LAL status structure pointer */
	       AMCoeffs *coeffs,			/**< [out] amplitude-coeffs {a(t_i), b(t_i)} */
	       const DetectorStateSeries *DetectorStates,/**< timeseries of detector states */
	       SkyPosition skypos			/**< {alpha,delta} of the source */
	       )
{
  REAL4 delta, alpha;
  REAL4 sin1delta, cos1delta;
  REAL4 sin1alpha, cos1alpha;

  REAL4 xi1, xi2;
  REAL4 eta1, eta2, eta3;
  REAL4 norm;
  UINT4 i, numSteps;

  INITSTATUS(status);

  /*---------- check input ---------- */
  ASSERT ( DetectorStates, status, LALCOMPUTEAMH_ENULL, LALCOMPUTEAMH_MSGENULL);

  numSteps = DetectorStates->length;

  /* require the coeffients-vectors to be allocated and consistent with timestamps */
  ASSERT ( coeffs, status, LALCOMPUTEAMH_ENULL, LALCOMPUTEAMH_MSGENULL);
  ASSERT ( coeffs->a && coeffs->b, status, LALCOMPUTEAMH_ENULL, LALCOMPUTEAMH_MSGENULL);
  ASSERT ( (coeffs->a->length == numSteps) && (coeffs->b->length == numSteps), status,
	   LALCOMPUTEAMH_EINPUT,  LALCOMPUTEAMH_MSGEINPUT);

  /* require sky-pos to be in equatorial coordinates */
  ASSERT ( skypos.system == COORDINATESYSTEM_EQUATORIAL, status,
	   SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );

  /*---------- We write components of xi and eta vectors in SSB-fixed coords */
  alpha = skypos.longitude;
  delta = skypos.latitude;

  sin_cos_LUT (&sin1delta, &cos1delta, delta );
  sin_cos_LUT (&sin1alpha, &cos1alpha, alpha );
  // see Eq.(17) in CFSv2 notes (version v3):
  // https://dcc.ligo.org/cgi-bin/private/DocDB/ShowDocument?docid=1665&version=3
  xi1 =   sin1alpha;
  xi2 =  -cos1alpha;
  eta1 = -sin1delta * cos1alpha;
  eta2 = -sin1delta * sin1alpha;
  eta3 = cos1delta;

  /*---------- Compute the a(t_i) and b(t_i) ---------- */
  coeffs->A = 0;
  coeffs->B = 0;
  coeffs->C = 0;
  coeffs->D = 0;
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

      /*
      printf("xi = (%f,%f)\n",xi1,xi2);
      printf("eta = (%f,%f,%f)\n",eta1,eta2,eta3);
      printf("d = (%f %f %f\n",d->d11,d->d12,d->d13);
      printf("     %f %f %f\n",d->d12,d->d22,d->d23);
      printf("     %f %f %f)\n",d->d13,d->d23,d->d33);
      */

      coeffs->a->data[i] = ai;
      coeffs->b->data[i] = bi;

      /* sum A, B, C on the fly */
      coeffs->A += ai * ai;
      coeffs->B += bi * bi;
      coeffs->C += ai * bi;

    } /* for i < numSteps */

  /* finish calculation of A,B,C, D */
  norm = 2.0f / numSteps;
  coeffs->A *= norm;
  coeffs->B *= norm;
  coeffs->C *= norm;

  coeffs->D = coeffs->A * coeffs->B - coeffs->C * coeffs->C;

  RETURN(status);

} /* LALNewGetAMCoeffs() */



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
  LALStatus status = empty_LALStatus;
  EarthState earth = empty_EarthState;

  if ( !ai || !bi || !skypos || !tGPS || !site || !edat) {
    XLAL_ERROR( XLAL_EINVAL );
  }

  LALBarycenterEarth (&status, &earth, tGPS, edat );
  if ( status.statusCode ) {
    XLALPrintError ("%s: Call to LALBarycenterEarth() failed. statusCode=%d\n", __func__, status.statusCode );
    XLAL_ERROR( XLAL_EFUNC );
  }

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
 * Multi-IFO version of LALGetAMCoeffs().
 * Get all antenna-pattern coefficients for all input detector-series.
 *
 * NOTE: contrary to LALGetAMCoeffs(), this functions *allocates* the output-vector,
 * use XLALDestroyMultiAMCoeffs() to free this.
 */
void
LALGetMultiAMCoeffs (LALStatus *status,			/**< [in/out] LAL status structure pointer */
		     MultiAMCoeffs **multiAMcoef,	/**< [out] AM-coefficients for all input detector-state series */
		     const MultiDetectorStateSeries *multiDetStates, /**< [in] detector-states at timestamps t_i */
		     SkyPosition skypos			/**< source sky-position [in equatorial coords!] */
		     )
{
  UINT4 X, numDetectors;
  MultiAMCoeffs *ret = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* check input */
  ASSERT (multiDetStates, status,LALCOMPUTEAMH_ENULL, LALCOMPUTEAMH_MSGENULL);
  ASSERT (multiDetStates->length, status,LALCOMPUTEAMH_ENULL, LALCOMPUTEAMH_MSGENULL);
  ASSERT (multiAMcoef, status,LALCOMPUTEAMH_ENULL, LALCOMPUTEAMH_MSGENULL);
  ASSERT ( *multiAMcoef == NULL, status,LALCOMPUTEAMH_ENONULL, LALCOMPUTEAMH_MSGENONULL);
  ASSERT ( skypos.system == COORDINATESYSTEM_EQUATORIAL, status, LALCOMPUTEAMH_EINPUT, LALCOMPUTEAMH_MSGEINPUT );

  numDetectors = multiDetStates->length;

  if ( ( ret = LALCalloc( 1, sizeof( *ret ) )) == NULL ) {
    ABORT (status, LALCOMPUTEAMH_EMEM, LALCOMPUTEAMH_MSGEMEM);
  }
  ret->length = numDetectors;
  if ( ( ret->data = LALCalloc ( numDetectors, sizeof ( *ret->data ) )) == NULL ) {
    LALFree ( ret );
    ABORT (status, LALCOMPUTEAMH_EMEM, LALCOMPUTEAMH_MSGEMEM);
  }

  for ( X=0; X < numDetectors; X ++ )
    {
      AMCoeffs *amcoeX = NULL;
      UINT4 numStepsX = multiDetStates->data[X]->length;

      ret->data[X] = LALCalloc ( 1, sizeof ( *(ret->data[X]) ) );
      amcoeX = ret->data[X];
      amcoeX->a = XLALCreateREAL4Vector ( numStepsX );
      if ( (amcoeX->b = XLALCreateREAL4Vector ( numStepsX )) == NULL ) {
	LALPrintError ("\nOut of memory!\n\n");
	goto failed;
      }

      /* LALGetAMCoeffs (status->statusPtr, amcoeX, multiDetStates->data[X], skypos ); */
      LALNewGetAMCoeffs (status->statusPtr, amcoeX, multiDetStates->data[X], skypos );
      if ( status->statusPtr->statusCode )
	{
	  LALPrintError ( "\nCall to LALNewGetAMCoeffs() has failed ... \n\n");
	  goto failed;
	}

    } /* for X < numDetectors */

  goto success;

 failed:
  /* free all memory allocated so far */
  XLALDestroyMultiAMCoeffs ( ret );
  ABORT ( status, -1, "LALGetMultiAMCoeffs() failed" );

 success:
  (*multiAMcoef) = ret;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALGetMultiAMCoeffs() */


/**
 * Original antenna-pattern function by S Berukoff
 */
void LALComputeAM (LALStatus          *status,
		   AMCoeffs           *coe,
		   LIGOTimeGPS        *ts,
		   AMCoeffsParams     *params)
{

  REAL4 zeta;                  /* sine of angle between detector arms        */
  INT4 i;                      /* temporary loop index                       */
  LALDetAMResponse response;   /* output of LALComputeDetAMResponse          */

  REAL4 sumA2=0.0;
  REAL4 sumB2=0.0;
  REAL4 sumAB=0.0;             /* variables to store scalar products         */
  INT4 length=coe->a->length;  /* length of input time series                */

  REAL4 cos2psi;
  REAL4 sin2psi;               /* temp variables                             */

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  /* Must put an ASSERT checking that ts vec and coe vec are same length */

  /* Compute the angle between detector arms, then the reciprocal */
  {
    LALFrDetector det = params->das->pDetector->frDetector;
    zeta = 1.0/(sin(det.xArmAzimuthRadians - det.yArmAzimuthRadians));
    if(params->das->pDetector->type == LALDETECTORTYPE_CYLBAR) zeta=1.0;
  }

  cos2psi = cos(2.0 * params->polAngle);
  sin2psi = sin(2.0 * params->polAngle);

  /* Note the length is the same for the b vector */
  for(i=0; i<length; ++i)
    {
      REAL4 *a = coe->a->data;
      REAL4 *b = coe->b->data;

      /* Compute F_plus, F_cross */
      LALComputeDetAMResponse(status->statusPtr, &response, params->das, &ts[i]);

      /*  Compute a, b from JKS eq 10,11
       *  a = zeta * (F_plus*cos(2\psi)-F_cross*sin(2\psi))
       *  b = zeta * (F_cross*cos(2\psi)+Fplus*sin(2\psi))
       */
      a[i] = zeta * (response.plus*cos2psi-response.cross*sin2psi);
      b[i] = zeta * (response.cross*cos2psi+response.plus*sin2psi);

      /* Compute scalar products */
      sumA2 += SQ(a[i]);                       /*  A  */
      sumB2 += SQ(b[i]);                       /*  B  */
      sumAB += (a[i]) * (b[i]);                /*  C  */
    }

  {
    /* Normalization factor */
    REAL8 L = 2.0/(REAL8)length;

    /* Assign output values and normalise */
    coe->A = L*sumA2;
    coe->B = L*sumB2;
    coe->C = L*sumAB;
    coe->D = ( coe->A * coe->B - SQ(coe->C) );
    /* protection against case when AB=C^2 */
    if(coe->D == 0) coe->D=1.0e-9;
  }

  /* Normal exit */

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* LALComputeAM() */

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

  REAL8 Ad = 0, Bd = 0, Cd = 0;	// multi-IFO values
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
          for(alpha = 0; alpha < numStepsX; alpha++)
            {
              REAL8 Sqwi = sqrt ( weightsX->data[alpha] );
              /* apply noise-weights, *replace* original a, b by noise-weighed version! */
	      amcoeX->a->data[alpha] *= Sqwi;
	      amcoeX->b->data[alpha] *= Sqwi;
            } // for alpha < numSteps
        } // if weights

      UINT4 alpha;	// SFT-index
      REAL8 AdX = 0, BdX = 0, CdX = 0;	// single-IFO values
      /* compute single-IFO antenna-pattern coefficients AX,BX,CX, by summing over time-steps 'alpha' */
      for(alpha = 0; alpha < numStepsX; alpha++)
        {
          REAL8 ahat = amcoeX->a->data[alpha];
          REAL8 bhat = amcoeX->b->data[alpha];

          AdX += ahat * ahat;
          BdX += bhat * bhat;
          CdX += ahat * bhat;
        } /* for alpha < numStepsX */

      /* store those */
      amcoeX->A = AdX;
      amcoeX->B = BdX;
      amcoeX->C = CdX;
      amcoeX->D = AdX * BdX - CdX * CdX;

      /* compute multi-IFO antenna-pattern coefficients A,B,C by summing over IFOs X */
      Ad += AdX;
      Bd += BdX;
      Cd += CdX;

    } /* for X < numDetectors */

  multiAMcoef->Mmunu.Ad = Ad;
  multiAMcoef->Mmunu.Bd = Bd;
  multiAMcoef->Mmunu.Cd = Cd;
  multiAMcoef->Mmunu.Dd = Ad * Bd - Cd * Cd;

  if ( multiWeights )
    multiAMcoef->Mmunu.Sinv_Tsft = multiWeights->Sinv_Tsft;

  return XLAL_SUCCESS;

} /* XLALWeightMultiAMCoeffs() */


/**
 * Compute the 'amplitude coefficients' \f$a(t)\sin\zeta\f$,
 * \f$b(t)\sin\zeta\f$ as defined in \ref JKS98 for a series of
 * timestamps.
 *
 * The input consists of the DetectorState-timeseries, which contains
 * the detector-info and the LMST's corresponding to the different times.
 *
 * \note This is an equivalent implementation to LALNewGetAMCoeffs(),
 * only using the XLAL interface instead.
 * This implementation is based on the geometrical definition of
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
  sin_cos_LUT (&sin1delta, &cos1delta, delta );
  sin_cos_LUT (&sin1alpha, &cos1alpha, alpha );

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
 * \note *) contrary to LALGetMultiAMCoeffs(), and XLALComputeAMCoeffs(), this function applies
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
    XLALPrintError ("%s: failed to XLALCalloc( 1, %d)\n", __func__, sizeof( *ret ) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  ret->length = numDetectors;
  if ( ( ret->data = XLALCalloc ( numDetectors, sizeof ( *ret->data ) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc(%d, %d)\n", __func__, numDetectors, sizeof ( *ret->data ) );
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
    XLALPrintError ("%s: failed to XLALCalloc ( 1, %d )\n", __func__, sizeof (*ret) );
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

