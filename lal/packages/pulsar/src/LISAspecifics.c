/*
 * Copyright (C) 2006 Reinhard Prix
 * Copyright (C) 2006 John T. Whelan
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

/** \author R. Prix, J. T. Whelan
 * \file
 * \brief
 * LISA-specific implementations for Fstat/continuous-wave searches on LISA TDI observables.
 *
 */

/*---------- INCLUDES ----------*/
#include <math.h>
#include <string.h>

#include <lal/LALError.h>
#include <ComputeFstat.h>

#include "LISAspecifics.h"


NRCSID( LISASPECIFICSC, "$Id$");

/*---------- local DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)

/*----- Macros ----- */
/** Simple Euklidean scalar product for two 3-dim vectors in cartesian coords */
#define SCALAR(u,v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])
#define SQ(x) ( (x) * (x) )
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )
/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- empty initializers ---------- */
const LALDetector empty_LALDetector;

/*---------- Global variables ----------*/

/*---------- internal prototypes ----------*/
int XLALgetLISAtwoArmRAAIFO ( CmplxDetectorTensor *detT, const DetectorArm *armA, const DetectorArm *armB, PulsarDopplerParams doppler);

/*==================== FUNCTION DEFINITIONS ====================*/


/** Set up the \em LALDetector struct representing LISA X, Y, Z TDI observables.
 * INPUT: channelNum = '1', '2', '3', '4', '5', '6': detector-tensor corresponding to TDIs X, Y, Z, Y-Z, Z-X, X-Y respectively.
 * return -1 on ERROR, 0 if OK
 */
int
XLALcreateLISA (LALDetector *Detector,	/**< [out] LALDetector */
		CHAR channelNum		/**< [in] which TDI observable: '1' = X, '2'= Y, '3' = Z, '4' = Y-Z, '5' = Z-X, '6' = X-Y */
		)
{
  LALDetector Detector1 = empty_LALDetector;

  if ( !Detector )
    return -1;

  switch ( channelNum )
    {
    case '1':
      strcpy ( Detector1.frDetector.name, "Z1: LISA TDI X" );
      strcpy ( Detector1.frDetector.prefix, "Z1");
      break;
    case '2':
      strcpy ( Detector1.frDetector.name, "Z2: LISA TDI Y" );
      strcpy ( Detector1.frDetector.prefix, "Z2");
      break;
    case '3':
      strcpy ( Detector1.frDetector.name, "Z3: LISA TDI Z" );
      strcpy ( Detector1.frDetector.prefix, "Z3");
      break;
    case '4':
      strcpy ( Detector1.frDetector.name, "Z4: LISA TDI Y-Z" );
      strcpy ( Detector1.frDetector.prefix, "Z4");
      break;
    case '5':
      strcpy ( Detector1.frDetector.name, "Z5: LISA TDI Z-X" );
      strcpy ( Detector1.frDetector.prefix, "Z5");
      break;
    case '6':
      strcpy ( Detector1.frDetector.name, "Z6: LISA TDI X-Y" );
      strcpy ( Detector1.frDetector.prefix, "Z6");
      break;
    default:
      LALPrintError ("\nIllegal LISA TDI index '%c': must be one of {'1', '2', '3', '4', '5', '6'}.\n\n", channelNum );
      return -1;
      break;
    } /* switch (detIndex) */

  /* fill frDetector with dummy numbers: meaningless for LISA */
  Detector1.frDetector.vertexLongitudeRadians = 0;
  Detector1.frDetector.vertexLatitudeRadians = 0;
  Detector1.frDetector.vertexElevation = 0;
  Detector1.frDetector.xArmAltitudeRadians = 0;
  Detector1.frDetector.xArmAzimuthRadians = 0;
  Detector1.frDetector.yArmAltitudeRadians = 0;
  Detector1.frDetector.yArmAzimuthRadians = 0;

  Detector1.type = LALDETECTORTYPE_ABSENT;

  /* however: need to be careful to put some non-zero numbers for location
   * otherwise LALBarycenter() will spit out NaNs ...
   */
  Detector1.location[0] = 1;
  Detector1.location[1] = 1;
  Detector1.location[2] = 1;

  (*Detector) = Detector1;

  return 0;

} /* XLALcreateLISA() */


/** Precompute the arm-geometry for LISA, which is later used
 * for assembling the RAA detector-tensor (which depends on
 * frequency and skyposition
 */
int
XLALprecomputeLISAarms ( DetectorState *detState )
{
  REAL4 x[3], y[3], z[3];		/* position x,y,z of spacecraft n */

  enum { SC1 = 0, SC2, SC3 };		/* translate TDI-conventions to C-indexing */
  enum { I1 = 0, I2, I3 };		/* convenience indices */

  UINT4 send_l[3] = { SC3, SC1, SC2 };	/* canonical sender-spacecraft 's' for arm 'l' */
  UINT4 rec_l [3] = { SC2, SC3, SC1 };	/* canonical receiver-spacecraft 'r' for arm 'l' */

  UINT4 arm;

  if ( !detState )
    return -1;

  { /* determine space-craft positions [MLDC Challenge 1] */
    REAL4 e = 0.00965;				/* LISA eccentricity */
    REAL4 ae = LAL_AU_SI * e / LAL_C_SI;	/* in units of time */
    REAL4 kappa = 0, lambda = 0;		/* MLDC default */
    REAL4 sin_alpha, cos_alpha;

    REAL8 Om = LAL_TWOPI /  LAL_YRSID_SI;
    REAL8 alpha_t;
    UINT4 i;
    REAL8 tGPS = GPS2REAL8( detState->tGPS ) - LISA_TIME_ORIGIN;

    alpha_t = Om * tGPS + kappa;
    sin_cos_LUT (&sin_alpha, &cos_alpha, alpha_t);

    for ( i = 0; i < 3; i ++ )
      {
	REAL4 beta = 2.0 * i * LAL_PI / 3.0  + lambda;	/* relative orbital phase in constellation */
	REAL4 sin_beta, cos_beta;
	sin_cos_LUT ( &sin_beta, &cos_beta, beta );

	x[i] = detState->earthState.posNow[0] + ae * ( sin_alpha * cos_alpha * sin_beta  - ( 1.0f + SQ(sin_alpha)) * cos_beta );
	y[i] = detState->earthState.posNow[1] + ae * ( sin_alpha * cos_alpha * cos_beta  - ( 1.0f + SQ(cos_alpha)) * sin_beta );
	z[i] = detState->earthState.posNow[2] - sqrt(3.0f) * ae * ( cos_alpha * cos_beta + sin_alpha * sin_beta );
      } /* for i=[0:2] */

  } /* determine spacecraft positions */

  for ( arm = 0; arm < 3; arm ++ )
    {
      UINT4 send, rec;
      REAL4 n[3], L, invL;

      /* get canonical (ie following TDI conventions) senders and receivers this arm */
      send = send_l[arm];
      rec  = rec_l [arm];

      /* get un-normalized arm-vectors first */
      n[0] = x[rec] - x[send];
      n[1] = y[rec] - y[send];
      n[2] = z[rec] - z[send];

      /* get armlength in seconds, and normalize */
      detState->detArms[arm].L_c = L = sqrt ( SCALAR(n,n) );
      invL = 1.0f / L;

      n[0] *= invL;
      n[1] *= invL;
      n[2] *= invL;

      /* pre-compute the "basis tensor" n x n  for this arm */
      XLALTensorSquareVector ( &(detState->detArms[arm].basisT), n[0], n[1], n[2] );

      /* store arm unit-vector */
      detState->detArms[arm].n[0] = n[0];
      detState->detArms[arm].n[1] = n[1];
      detState->detArms[arm].n[2] = n[2];

    } /* for arm = 0, 1, 2 */

  return 0;

} /* XLALprecomputeLISAarms() */


/* Construct the long-wavelength-limit (LWL) detector tensor for LISA, given the prefix
 * and the precomputed detector arms.
 *
 * RETURN 0 = OK, -1 = ERROR
 */
int
XLALgetLISADetectorTensorLWL ( DetectorTensor *detT, 		/**< [out]: LISA LWL detector-tensor */
			       const Detector3Arms detArms,	/**< [in] precomputed detector-arms */
			       CHAR channelNum )		/**< channel-number (as a char)  '1', '2', '3' .. */
{
  LISAarmT armA, armB;
  CHAR chan1 = 0;
  CHAR chan2 = 0;
  DetectorTensor detT1, detT2;

  if ( !detT )
    return -1;

  /* we distinuish (currently) 3 different TDI observables: X (channel=1), Y (channel=2), Z (channel=3) */
  switch ( channelNum )
    {
    case '1': 	/* TDI observable 'X' */
      armA = LISA_ARM2; armB = LISA_ARM3;
      break;
    case '2':		/* TDI observable 'Y' */
      armA = LISA_ARM3; armB = LISA_ARM1;
      break;
    case '3':		/* TDI observable 'Z' */
      armA = LISA_ARM1; armB = LISA_ARM2;
      break;
    case '4':
      chan1 = '2'; chan2 = '3';
      break;
    case '5':
      chan1 = '3'; chan2 = '1';
      break;
    case '6':
      chan1 = '1'; chan2 = '2';
      break;
    default:	/* unknown */
      LALPrintError ("\nInvalid channel-number '%c' for LISA \n\n", channelNum );
      xlalErrno = XLAL_EINVAL;
      return -1;
      break;
    } /* switch channel[1] */

  if (chan1 == 0)
    {
      XLALSubtractDetectorTensors ( detT, &(detArms[armA].basisT), &(detArms[armB].basisT) );
      XLALMultiplyDetectorTensor ( detT, detT, 0.5f );	/* (1/2)*(nA x nA - nB x nB ) */
    }
  else
    {
      if ( XLALgetLISADetectorTensorLWL ( &detT1, detArms, chan1 ) != 0 ) {
	LALPrintError ("\nXLALgetLISADetectorTensorLWL() failed !\n\n");
	xlalErrno = XLAL_EINVAL;
	return -1;
      }
      if ( XLALgetLISADetectorTensorLWL ( &detT2, detArms, chan2 ) != 0 ) {
	LALPrintError ("\nXLALgetLISADetectorTensorLWL() failed !\n\n");
	xlalErrno = XLAL_EINVAL;
	return -1;
      }

      XLALSubtractDetectorTensors ( detT, &detT1, &detT2 );	/* d_X - d_Y etc */
    } /* multi-channel "detector" such as 'X-Y' etc */

  return 0;

} /* XLALgetLISADetectorTensorLWL() */


/* Construct the rigid-adiabatic-approximation (RAA) detector tensor for LISA, given the prefix,
 * the precomputed arms, and the Frequency and skyposition.
 *
 * RETURN 0 = OK, -1 = ERROR
 */
int
XLALgetLISADetectorTensorRAA ( CmplxDetectorTensor *detT, 	/**< [out]: LISA LWL detector-tensor */
			       const Detector3Arms detArms,	/**< [in] precomputed detector-arms */
			       PulsarDopplerParams doppler,   	/**< [in] doppler parameters including frequency and sky position */
			       CHAR channelNum)			/**< [in] channel-number (as a char)  '1', '2', '3' .. */
{
  LISAarmT armA, armB;
  CHAR chan1 = 0;
  CHAR chan2 = 0;
  CmplxDetectorTensor detT1, detT2;

  if ( !detT )
    return -1;

  /* we distinuish (currently) 3 different TDI observables: X (channel=1), Y (channel=2), Z (channel=3) */
  switch ( channelNum )
    {
    case '1': 	/* TDI observable 'X' */
      armA = LISA_ARM2; armB = LISA_ARM3;
      break;
    case '2':		/* TDI observable 'Y' */
      armA = LISA_ARM3; armB = LISA_ARM1;
      break;
    case '3':		/* TDI observable 'Z' */
      armA = LISA_ARM1; armB = LISA_ARM2;
      break;
    case '4':
      chan1 = '2'; chan2 = '3';
      break;
    case '5':
      chan1 = '3'; chan2 = '1';
      break;
    case '6':
      chan1 = '1'; chan2 = '2';
      break;
    default:	/* unknown */
      LALPrintError ("\nInvalid channel-number '%c' for LISA \n\n", channelNum );
      xlalErrno = XLAL_EINVAL;
      return -1;
      break;
    } /* switch channel[1] */

  if (chan1 == 0) {
    if ( XLALgetLISAtwoArmRAAIFO ( detT, &(detArms[armA]), &(detArms[armB]), doppler ) != 0 ) {
      LALPrintError ("\nXLALgetLISAtwoArmRAAIFO() failed !\n\n");
      xlalErrno = XLAL_EINVAL;
      return -1;
    }
  } else {
    if ( XLALgetLISADetectorTensorRAA ( &detT1, detArms, doppler, chan1 ) != 0 ) {
      LALPrintError ("\nXLALgetCmplxLISADetectorTensorRAA() failed !\n\n");
      xlalErrno = XLAL_EINVAL;
      return -1;
    }
    if ( XLALgetLISADetectorTensorRAA ( &detT2, detArms, doppler, chan2 ) != 0 ) {
      LALPrintError ("\nXLALgetCmplxLISADetectorTensorRAA() failed !\n\n");
      xlalErrno = XLAL_EINVAL;
      return -1;
    }
    XLALSubtractDetectorTensors ( &detT->re, &detT1.re, &detT2.re );
    XLALSubtractDetectorTensors ( &detT->im, &detT1.im, &detT2.im );
  }

  return 0;

} /* XLALgetCmplxLISADetectorTensor() */


/** return a rigid-adiabatic-approximation (RAA) two-arm IFO detector tensor for LISA, given 'armA' and 'armB'
 * This implements LISA RAA-tensor using spacecraft-orbits described by Eq.(2.1) in LISA-MLCD
 * 'challenge1.pdf' document, see http://astrogravs.nasa.gov/docs/mldc/round1.html
 *
 * RETURN 0 = OK, -1 = ERROR
 *
 * \note: armA=0 means "arm 1" in TDI notation, therefore the enums LISA_ARM1, LISA_ARM2, LISA_ARM3
 * should be used to avoid confusion.
 */
int
XLALgetLISAtwoArmRAAIFO ( CmplxDetectorTensor *detT, 	/**< [out]: two-arm IFO detector-tensor */
			  const DetectorArm *detArmA,	/**< [in] precomputed detector-arm 'A' */
			  const DetectorArm *detArmB,	/**< [in] precomputed detector-arm 'B' */
			  PulsarDopplerParams doppler   /**< [in] doppler parameters including frequency and sky position */
			  )
{

  REAL4 pifL_c = LAL_PI_2 * doppler.fkdot[0] * (detArmA->L_c + detArmB->L_c);
  REAL4 pifL_3c = pifL_c / 3.0f;
  REAL4 sin1Delta, cos1Delta;
  REAL4 sin1Alpha, cos1Alpha;
  REAL4 kdotnA, kdotnB;
  REAL4 eta, sinpha, cospha;
  REAL4 sinc_eta, sineta, coseta;
  REAL4 k[3];
  COMPLEX8 coeffAA, coeffBB;

  if ( !detT || !detArmA || !detArmB )
    return -1;

  sin_cos_LUT (&sin1Delta, &cos1Delta, doppler.Delta );
  sin_cos_LUT (&sin1Alpha, &cos1Alpha, doppler.Alpha );

  k[0] = - cos1Delta * cos1Alpha;
  k[1] = - cos1Delta * sin1Alpha;
  k[2] = - sin1Delta;

  kdotnA = SCALAR ( k, (detArmA->n) );
  kdotnB = SCALAR ( k, (detArmB->n) );

  /* ----- calculate complex coefficient in front of basis-tensor (1/2)*(nA x nA) ----- */

  /* first term */
  sin_cos_LUT (&sinpha, &cospha, pifL_3c * ( 3.0f - (kdotnA + 2.0f * kdotnB) ) );

  eta = pifL_c * (1.0f + kdotnA);
  sin_cos_LUT ( &sineta, &coseta, eta );
  sinc_eta = sineta / eta;

  coeffAA.re = 0.5f * cospha * sinc_eta;
  coeffAA.im = 0.5f * sinpha * sinc_eta;

  /* second term */
  sin_cos_LUT (&sinpha, &cospha, pifL_3c * ( - 3.0f - (kdotnA + 2.0f * kdotnB) ) );

  eta = pifL_c * (1.0f - kdotnA);
  sin_cos_LUT ( &sineta, &coseta, eta );
  sinc_eta = sineta / eta;

  coeffAA.re += 0.5f * cospha * sinc_eta;
  coeffAA.im += 0.5f * sinpha * sinc_eta;

  /* ----- calculate coefficient in front of (1/2)*(nB x nB) */

  /* first term */
  sin_cos_LUT (&sinpha, &cospha, pifL_3c * ( 3.0f + (2.0f * kdotnA + kdotnB) ) );

  eta = pifL_c * (1.0f - kdotnB);
  sin_cos_LUT ( &sineta, &coseta, eta );
  sinc_eta = sineta / eta;

  coeffBB.re = 0.5f * cospha * sinc_eta;
  coeffBB.im = 0.5f * sinpha * sinc_eta;

  /* second term */
  sin_cos_LUT (&sinpha, &cospha, pifL_3c * ( - 3.0f + (2.0f * kdotnA + kdotnB) ) );

  eta = pifL_c * (1.0f + kdotnB);
  sin_cos_LUT ( &sineta, &coseta, eta );
  sinc_eta = sineta / eta;

  coeffBB.re += 0.5f * cospha * sinc_eta;
  coeffBB.im += 0.5f * sinpha * sinc_eta;

  /* now we can express the detector tensor in the rigid adiabatic approximation (RAA): 
   * detT = coeffAA * basisA - coeffBB * basisB
   */
  {
    DetectorTensor tmpA, tmpB;
    XLALMultiplyDetectorTensor ( &tmpA, &detArmA->basisT, coeffAA.re );
    XLALMultiplyDetectorTensor ( &tmpB, &detArmB->basisT, coeffBB.re );
    XLALSubtractDetectorTensors( &detT->re, &tmpA, &tmpB );

    XLALMultiplyDetectorTensor ( &tmpA, &detArmA->basisT, coeffAA.im );
    XLALMultiplyDetectorTensor ( &tmpB, &detArmB->basisT, coeffBB.im );
    XLALSubtractDetectorTensors( &detT->im, &tmpA, &tmpB );
  }

  return 0;

} /* XLALgetLISAtwoArmRAAIFO() */
