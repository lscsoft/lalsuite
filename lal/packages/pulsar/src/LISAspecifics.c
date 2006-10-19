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
#define __USE_ISOC99 1
#include <math.h>

#include "ComputeFstat.h"

NRCSID( LISASPECIFICS, "$Id$");

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

/*---------- Global variables ----------*/

/*---------- internal prototypes ----------*/

/*==================== FUNCTION DEFINITIONS ====================*/

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
    REAL8 ti = XLALGPSGetREAL8( &tGPS );
    REAL8 Om = LAL_TWOPI /  LAL_YRSID_SI;
    REAL8 alpha_ti  = Om * ti + kappa;
    REAL4 sin_alpha, cos_alpha;
    REAL4 sqrt3 = sqrt(3.0f);

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
  nA[0] = x[recA] - x[sendA]; 
  nA[1] = y[recA] - y[sendA]; 
  nA[2] = z[recA] - z[sendA];
      
  nB[0] = x[recB] - x[sendB]; 
  nB[1] = y[recB] - y[sendB]; 
  nB[2] = z[recB] - z[sendB];
  
  /* get lenths and normalize */
  LA = sqrt ( SCALAR(nA,nA) );
  LB = sqrt ( SCALAR(nB,nB) );
  nA[0] /= LA;	nA[1] /= LA; 	nA[2] /= LA;
  nB[0] /= LB; 	nB[1] /= LB; 	nB[2] /= LB;

  /* now we can express the detector tensor in the long-wavelength limit (LWL): */
  detT->d11 =  0.5 * ( nA[I1] * nA[I1] - nB[I1] * nB[I1] );
  detT->d12 =  0.5 * ( nA[I1] * nA[I2] - nB[I1] * nB[I2] );
  detT->d13 =  0.5 * ( nA[I1] * nA[I3] - nB[I1] * nB[I3] );
  
  detT->d22 =  0.5 * ( nA[I2] * nA[I2] - nB[I2] * nB[I2] );
  detT->d23 =  0.5 * ( nA[I2] * nA[I3] - nB[I2] * nB[I3] );
  
  detT->d33 =  0.5 * ( nA[I3] * nA[I3] - nB[I3] * nB[I3] );

  return 0;

} /* XLALgetLISAtwoArmIFO() */

