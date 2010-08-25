/*
*  Copyright (C) 2007 Jolien Creighton, Maria Alessandra Papa, Steve Berukoff
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
 * \author Berukoff, S.J.
 * \date 2005
 * \ingroup pulsarAntenna
 * \file
 * \brief Header-file for computing antenna-pattern components for amplitude demodulation.
 *
 * <tt>#include <lal/LALComputeAM.h></tt>
 *
 * In order to compute the optimal statistic for pulsar searches, one must take account of the
 * various modulations that change the emitted, (fairly) simple sinusoid into a non-trivial function
 * of parameters.  The frequency evolution of the signal (spindown effects, Doppler modulation, etc.)
 * have already been accounted for; this routine filters the amplitude modulation effects.
 */

#ifndef _LALCOMPUTEAM_H
#define _LALCOMPUTEAM_H

#include <math.h>
#include <lal/DetResponse.h>
#include <lal/DetectorSite.h>
#include <lal/LALBarycenter.h>

#ifdef __cplusplus
extern "C" {
#endif

NRCSID (LALCOMPUTEAMH, "$Id: LALComputeAM.h");

/** \name Error codes */
/*@{*/
#define LALCOMPUTEAMH_ENOTS        1
#define LALCOMPUTEAMH_EBCERR       2
#define LALCOMPUTEAMH_EESERR       3
#define LALCOMPUTEAMH_EEPH         4
#define LALCOMPUTEAMH_EDAS         5
#define LALCOMPUTEAMH_EFRD         6

#define LALCOMPUTEAMH_MSGENOTS    "Input LIGOTimeGPS Vector is wrong size or NULL"
#define LALCOMPUTEAMH_MSGEBCERR   "Baryinput pointer is invalid"
#define LALCOMPUTEAMH_MSGEESERR   "EarthState structure invalid, or pointer NULL"
#define LALCOMPUTEAMH_MSGEEPH     "Ephemeris Table invalid, or pointer NULL"
#define LALCOMPUTEAMH_MSGEDAS     "Detector and source information invalid, or pointer NULL"
#define LALCOMPUTEAMH_MSGEFRD     "Detector geometry information invalid, or pointer NULL"
/*@}*/

/** This structure contains the output of the routine: a(t), b(t),
 * and the scalar products therein.  That is:
 */
typedef struct AMCoeffsTag
{
  REAL4Vector     *a;          /**< the function a(t)         */
  REAL4Vector     *b;          /**< the function b(t)         */
  REAL4           A;           /**< the scalar product (a||a) */
  REAL4           B;           /**< the scalar product (b||b) */
  REAL4           C;           /**< the scalar product (a||b) */
  REAL4           D;           /**< the quantity AB-C^2       */
} AMCoeffs;

typedef struct CmplxAMCoeffsTag
{
  COMPLEX8Vector     *a;          /**< the a coefficient evaluated at the relevant times */
  COMPLEX8Vector     *b;          /**< the b coefficient evaluated at the relevant times  */
} CmplxAMCoeffs;

/** This structure contains the parameters for the routine.  They include:
 */
typedef struct AMCoeffsParamsTag
{
  BarycenterInput      *baryinput;  /**< data from Barycentring routine */
  EarthState           *earth;      /**< from LALBarycenter()           */
  EphemerisData        *edat;       /**< the ephemerides                */
  LALDetAndSource      *das;        /**< det and source information     */
  LALFrDetector        *det;        /**< detector geometry              */
  REAL4                polAngle;    /**< polarization angle             */
} AMCoeffsParams;


/* exported prototypes */
void LALComputeAM (LALStatus          *status,
		   AMCoeffs           *coe,
		   LIGOTimeGPS        *ts,
		   AMCoeffsParams     *params);

#ifdef __cplusplus
}
#endif

#endif /* _LALCOMPUTEAM_H */
