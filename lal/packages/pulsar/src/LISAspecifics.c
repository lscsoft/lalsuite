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

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- empty initializers ---------- */
const LALDetector empty_LALDetector;

/*---------- Global variables ----------*/

/*---------- internal prototypes ----------*/
int XLALgetLISAtwoArmIFO ( DetectorTensor *detT, LIGOTimeGPS tGPS, LISAarmT armA, LISAarmT armB );

int XLALgetLISAtwoArmRAAIFO ( CmplxDetectorTensor *detT, const DetectorState *detState, PulsarDopplerParams doppler, LISAarmT armA, LISAarmT armB );

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

/* Construct the long-wavelength-limit (LWL) detector tensor for LISA, given the prefix
 *
 * RETURN 0 = OK, -1 = ERROR
 */
int
XLALgetLISADetectorTensor ( DetectorTensor *detT, 	/**< [out]: LISA LWL detector-tensor */
			    LIGOTimeGPS tGPS,		/**< [in] GPS time to compute IFO at */
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

  if (chan1 == 0) {
    if ( XLALgetLISAtwoArmIFO ( detT, tGPS, armA, armB ) != 0 ) {
      LALPrintError ("\nXLALgetLISAtwoArmIFO() failed !\n\n");
      xlalErrno = XLAL_EINVAL;
      return -1;
    }
  } else {
    if ( XLALgetLISADetectorTensor ( &detT1, tGPS, chan1 ) != 0 ) {
      LALPrintError ("\nXLALgetLISADetectorTensor() failed !\n\n");
      xlalErrno = XLAL_EINVAL;
      return -1;
    }
    if ( XLALgetLISADetectorTensor ( &detT2, tGPS, chan2 ) != 0 ) {
      LALPrintError ("\nXLALgetLISADetectorTensor() failed !\n\n");
      xlalErrno = XLAL_EINVAL;
      return -1;
    }
    detT->d11 = detT1.d11 - detT2.d11;
    detT->d12 = detT1.d12 - detT2.d12;
    detT->d13 = detT1.d13 - detT2.d13;
    detT->d22 = detT1.d22 - detT2.d22;
    detT->d23 = detT1.d23 - detT2.d23;
    detT->d33 = detT1.d33 - detT2.d33;
  }

  return 0;

} /* XLALgetLISADetectorTensor() */


/** return a long-wavelength-limit (LWL) two-arm IFO detector tensor for LISA, given 'armA' and 'armB'
 * This implements LISA LWL-tensor using spacecraft-orbits described by Eq.(2.1) in LISA-MLCD
 * 'challenge1.pdf' document, see http://astrogravs.nasa.gov/docs/mldc/round1.html
 *
 * RETURN 0 = OK, -1 = ERROR
 *
 * \note: armA=0 means "arm 1" in TDI notation, therefore the enums LISA_ARM1, LISA_ARM2, LISA_ARM3
 * should be used to avoid confusion.
 */
int
XLALgetLISAtwoArmIFO ( DetectorTensor *detT, 	/**< [out]: two-arm IFO detector-tensor */
		       LIGOTimeGPS tGPS,	/**< [in] GPS time to compute IFO at */
		       LISAarmT armA, 		/**< [in] first arm */
		       LISAarmT armB )		/**< [in] second arm */
{
  UINT4 n;				/* space-craft index 1,2,3 */
  REAL4 x[3], y[3], z[3];		/* position x,y,z of spacecraft n */

  enum { SC1 = 0, SC2, SC3 };		/* translate TDI-conventions to C-indexing */
  enum { I1 = 0, I2, I3 };		/* translate spacial indices 1,2,3 to C-indexing  */

  UINT4 send_l[3] = { SC3, SC1, SC2 };	/* sender-spacecraft 's' for arm 'l' */
  UINT4 rec_l [3] = { SC2, SC3, SC1 };	/* receiver-spacecraft 'r' for arm 'l' */

  UINT4 sendA, recA, sendB, recB; 	/* respective senders/receivers for arms 'A' and 'B' */
  REAL4 nA[3], LA, nB[3], LB;		/* unit-vectors and armlength of the two arms involved */

  { /* determine space-craft positions [MLDC Challenge 1] */
    REAL4 a = LAL_AU_SI;
    REAL4 e = 0.00965;		/* eccentricity */
    REAL4 kappa = 0, lambda = 0;	/* MLDC default */
    REAL8 ti;
    REAL8 Om = LAL_TWOPI /  LAL_YRSID_SI;
    REAL8 alpha_ti;
    REAL4 sin_alpha, cos_alpha;
    REAL4 sqrt3;

    sqrt3 = sqrt(3.0f);
    ti = XLALGPSGetREAL8( &tGPS ) - LISA_TIME_ORIGIN;
    alpha_ti = Om * ti + kappa;

    sin_cos_LUT (&sin_alpha, &cos_alpha, alpha_ti);
    for ( n = 1; n <= 3; n ++ )
      {
	REAL4 beta = 2.0 * ( n - 1.0 ) * LAL_PI / 3.0  + lambda;	/* relative orbital phase in constellation */
	REAL4 sin_beta, cos_beta;
	sin_cos_LUT ( &sin_beta, &cos_beta, beta );

	x[n-1] = a * cos_alpha  + a * e * ( sin_alpha * cos_alpha * sin_beta  - ( 1.0 + SQ(sin_alpha)) * cos_beta );
	y[n-1] = a * sin_alpha  + a * e * ( sin_alpha * cos_alpha * cos_beta  - ( 1.0 + SQ(cos_alpha)) * sin_beta );
	z[n-1] = - sqrt3 * a * e * ( cos_alpha * cos_beta + sin_alpha * sin_beta );
      } /* for n = 1,2,3 */
  } /* determine spacecraft positions x[n-1], y[n-1], z[n-1] */

  /* get corresponding senders and receivers for the arms involved */
  sendA = send_l[armA];
  recA  = rec_l [armA];

  sendB = send_l[armB];
  recB  = rec_l [armB];

  /* get un-normalized arm-vectors first */
  nA[I1] = x[recA] - x[sendA];
  nA[I2] = y[recA] - y[sendA];
  nA[I3] = z[recA] - z[sendA];

  nB[I1] = x[recB] - x[sendB];
  nB[I2] = y[recB] - y[sendB];
  nB[I3] = z[recB] - z[sendB];

  /* get lenths and normalize */
  LA = sqrt ( SCALAR(nA,nA) );
  LB = sqrt ( SCALAR(nB,nB) );
  nA[I1] /= LA;	nA[I2] /= LA; 	nA[I3] /= LA;
  nB[I1] /= LB; 	nB[I2] /= LB; 	nB[I3] /= LB;

  /* now we can express the detector tensor in the long-wavelength limit (LWL): */
  detT->d11 =  0.5 * ( nA[I1] * nA[I1] - nB[I1] * nB[I1] );
  detT->d12 =  0.5 * ( nA[I1] * nA[I2] - nB[I1] * nB[I2] );
  detT->d13 =  0.5 * ( nA[I1] * nA[I3] - nB[I1] * nB[I3] );

  detT->d22 =  0.5 * ( nA[I2] * nA[I2] - nB[I2] * nB[I2] );
  detT->d23 =  0.5 * ( nA[I2] * nA[I3] - nB[I2] * nB[I3] );

  detT->d33 =  0.5 * ( nA[I3] * nA[I3] - nB[I3] * nB[I3] );

  return 0;

} /* XLALgetLISAtwoArmIFO() */


/* Construct the rigid-adiabatic-approximation (RAA) detector tensor for LISA, given the prefix
 *
 * RETURN 0 = OK, -1 = ERROR
 */
int
XLALgetCmplxLISADetectorTensor ( CmplxDetectorTensor *detT, 	/**< [out]: LISA LWL detector-tensor */
				 const DetectorState *detState, /**< [in] detector-state info */
				 PulsarDopplerParams doppler,   /**< [in] doppler parameters including frequency and sky position */
				 CHAR channelNum)		/**< [in] channel-number (as a char)  '1', '2', '3' .. */
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
    if ( XLALgetLISAtwoArmRAAIFO ( detT, detState, doppler, armA, armB ) != 0 ) {
      LALPrintError ("\nXLALgetLISAtwoArmRAAIFO() failed !\n\n");
      xlalErrno = XLAL_EINVAL;
      return -1;
    }
  } else {
    if ( XLALgetCmplxLISADetectorTensor ( &detT1, detState, doppler, chan1 ) != 0 ) {
      LALPrintError ("\nXLALgetCmplxLISADetectorTensor() failed !\n\n");
      xlalErrno = XLAL_EINVAL;
      return -1;
    }
    if ( XLALgetCmplxLISADetectorTensor ( &detT2, detState, doppler, chan2 ) != 0 ) {
      LALPrintError ("\nXLALgetCmplxLISADetectorTensor() failed !\n\n");
      xlalErrno = XLAL_EINVAL;
      return -1;
    }
    detT->d11.re = detT1.d11.re - detT2.d11.re;
    detT->d11.im = detT1.d11.im - detT2.d11.im;
    detT->d12.re = detT1.d12.re - detT2.d12.re;
    detT->d12.im = detT1.d12.im - detT2.d12.im;
    detT->d13.re = detT1.d13.re - detT2.d13.re;
    detT->d13.im = detT1.d13.im - detT2.d13.im;
    detT->d22.re = detT1.d22.re - detT2.d22.re;
    detT->d22.im = detT1.d22.im - detT2.d22.im;
    detT->d23.re = detT1.d23.re - detT2.d23.re;
    detT->d23.im = detT1.d23.im - detT2.d23.im;
    detT->d33.re = detT1.d33.re - detT2.d33.re;
    detT->d33.im = detT1.d33.im - detT2.d33.im;
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
			  const DetectorState *detState, /**< [in] detector-state info */
			  PulsarDopplerParams doppler,   /**< [in] doppler parameters including frequency and sky position */
			  LISAarmT armA, 		/**< [in] first arm */
			  LISAarmT armB )		/**< [in] second arm */
{
  UINT4 n;				/* space-craft index 1,2,3 */
  REAL4 x[3], y[3], z[3];		/* position x,y,z of spacecraft n */

  enum { SC1 = 0, SC2, SC3 };		/* translate TDI-conventions to C-indexing */
  enum { I1 = 0, I2, I3 };		/* translate spacial indices 1,2,3 to C-indexing  */

  UINT4 send_l[3] = { SC3, SC1, SC2 };	/* sender-spacecraft 's' for arm 'l' */
  UINT4 rec_l [3] = { SC2, SC3, SC1 };	/* receiver-spacecraft 'r' for arm 'l' */

  UINT4 sendA, recA, sendB, recB; 	/* respective senders/receivers for arms 'A' and 'B' */
  REAL4 nA[3], LA, nB[3], LB;		/* unit-vectors and armlength of the two arms involved */
  REAL4 pifL_c;
  REAL4 delta, alpha;
  REAL4 sin1delta, cos1delta;
  REAL4 sin1alpha, cos1alpha;
  REAL4 k[3];
  REAL4 kdotnA, kdotnB;
  REAL4 sinpha, cospha, eta, sinceta;
  COMPLEX8 coeffAA, coeffBB;

  { /* determine space-craft positions [MLDC Challenge 1] */
    REAL4 a = LAL_AU_SI;
    REAL4 e = 0.00965;		/* eccentricity */
    REAL4 kappa = 0, lambda = 0;	/* MLDC default */
    REAL8 ti;
    REAL8 Om = LAL_TWOPI /  LAL_YRSID_SI;
    REAL8 alpha_ti;
    REAL4 sin_alpha, cos_alpha;
    REAL4 sqrt3;

    sqrt3 = sqrt(3.0f);
    ti = XLALGPSGetREAL8( &detState->tGPS ) - LISA_TIME_ORIGIN;
    alpha_ti = Om * ti + kappa;

    sin_cos_LUT (&sin_alpha, &cos_alpha, alpha_ti);
    for ( n = 1; n <= 3; n ++ )
      {
	REAL4 beta = 2.0 * ( n - 1.0 ) * LAL_PI / 3.0  + lambda;	/* relative orbital phase in constellation */
	REAL4 sin_beta, cos_beta;
	sin_cos_LUT ( &sin_beta, &cos_beta, beta );

	x[n-1] = a * cos_alpha  + a * e * ( sin_alpha * cos_alpha * sin_beta  - ( 1.0 + SQ(sin_alpha)) * cos_beta );
	y[n-1] = a * sin_alpha  + a * e * ( sin_alpha * cos_alpha * cos_beta  - ( 1.0 + SQ(cos_alpha)) * sin_beta );
	z[n-1] = - sqrt3 * a * e * ( cos_alpha * cos_beta + sin_alpha * sin_beta );
      } /* for n = 1,2,3 */
  } /* determine spacecraft positions x[n-1], y[n-1], z[n-1] */

  /* get corresponding senders and receivers for the arms involved */
  sendA = send_l[armA];
  recA  = rec_l [armA];

  sendB = send_l[armB];
  recB  = rec_l [armB];

  /* get un-normalized arm-vectors first */
  nA[I1] = x[recA] - x[sendA];
  nA[I2] = y[recA] - y[sendA];
  nA[I3] = z[recA] - z[sendA];

  nB[I1] = x[recB] - x[sendB];
  nB[I2] = y[recB] - y[sendB];
  nB[I3] = z[recB] - z[sendB];

  /* get lenths and normalize */
  LA = sqrt ( SCALAR(nA,nA) );
  LB = sqrt ( SCALAR(nB,nB) );
  nA[I1] /= LA;	nA[I2] /= LA; 	nA[I3] /= LA;
  nB[I1] /= LB; 	nB[I2] /= LB; 	nB[I3] /= LB;

  pifL_c = doppler.fkdot[0] * (LA + LB) * (LAL_PI_2 / LAL_C_SI);

  alpha = doppler.Alpha;
  delta = doppler.Delta;

  sin_cos_LUT (&sin1delta, &cos1delta, delta );
  sin_cos_LUT (&sin1alpha, &cos1alpha, alpha );

  k[I1] = - cos1delta * cos1alpha;
  k[I2] = - cos1delta * sin1alpha;
  k[I3] = - sin1delta;

  kdotnA = SCALAR(k,nA);
  kdotnB = SCALAR(k,nB);

  /* calculate coefficient of nA x nA */
  sin_cos_LUT (&sinpha, &cospha,
	       ( pifL_c / 3.0 ) * ( 3.0 - (kdotnA + 2.0*kdotnB) )
	       );
  eta = pifL_c * (1.0 + kdotnA);
  sinceta = sin(eta) / eta;
  coeffAA.re = cospha * sinceta / 4.0;
  coeffAA.im = sinpha * sinceta / 4.0;

  sin_cos_LUT (&sinpha, &cospha,
	       ( pifL_c / 3.0 ) * ( - 3.0 - (kdotnA + 2.0*kdotnB) )
	       );
  eta = pifL_c * (1.0 - kdotnA);
  sinceta = sin(eta) / eta;
  coeffAA.re += cospha * sinceta / 4.0;
  coeffAA.im += sinpha * sinceta / 4.0;

  /* calculate coefficient of nB x nB */
  sin_cos_LUT (&sinpha, &cospha,
	       ( pifL_c / 3.0 ) * ( 3.0 + (2.0*kdotnA + kdotnB) )
	       );
  eta = pifL_c * (1.0 - kdotnB);
  sinceta = sin(eta) / eta;
  coeffBB.re = cospha * sinceta / 4.0;
  coeffBB.im = sinpha * sinceta / 4.0;

  sin_cos_LUT (&sinpha, &cospha,
	       ( pifL_c / 3.0 ) * ( - 3.0 + (2.0*kdotnA + kdotnB) )
	       );
  eta = pifL_c * (1.0 + kdotnB);
  sinceta = sin(eta) / eta;
  coeffBB.re += cospha * sinceta / 4.0;
  coeffBB.im += sinpha * sinceta / 4.0;

  /* now we can express the detector tensor in the rigid adiabatic approximation (RAA): */
  detT->d11.re =  coeffAA.re * nA[I1] * nA[I1] - coeffBB.re * nB[I1] * nB[I1];
  detT->d12.re =  coeffAA.re * nA[I1] * nA[I2] - coeffBB.re * nB[I1] * nB[I2];
  detT->d13.re =  coeffAA.re * nA[I1] * nA[I3] - coeffBB.re * nB[I1] * nB[I3];

  detT->d22.re =  coeffAA.re * nA[I2] * nA[I2] - coeffBB.re * nB[I2] * nB[I2];
  detT->d23.re =  coeffAA.re * nA[I2] * nA[I3] - coeffBB.re * nB[I2] * nB[I3];

  detT->d33.re =  coeffAA.re * nA[I3] * nA[I3] - coeffBB.re * nB[I3] * nB[I3];

  detT->d11.im =  coeffAA.im * nA[I1] * nA[I1] - coeffBB.im * nB[I1] * nB[I1];
  detT->d12.im =  coeffAA.im * nA[I1] * nA[I2] - coeffBB.im * nB[I1] * nB[I2];
  detT->d13.im =  coeffAA.im * nA[I1] * nA[I3] - coeffBB.im * nB[I1] * nB[I3];

  detT->d22.im =  coeffAA.im * nA[I2] * nA[I2] - coeffBB.im * nB[I2] * nB[I2];
  detT->d23.im =  coeffAA.im * nA[I2] * nA[I3] - coeffBB.im * nB[I2] * nB[I3];

  detT->d33.im =  coeffAA.im * nA[I3] * nA[I3] - coeffBB.im * nB[I3] * nB[I3];

  return 0;

} /* XLALgetLISAtwoArmRAAIFO() */
