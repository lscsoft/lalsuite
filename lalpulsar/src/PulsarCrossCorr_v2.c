/*
 *  Copyright (C) 2012 John Whelan, Shane Larson and Badri Krishnan
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

#include <lal/PulsarCrossCorr_v2.h>

#define SQUARE(x) ((x)*(x))

/** Calculate the Doppler-shifted frequency associate with each SFT in a list */
/* This is according to Eqns 2.11 and 2.12 of Dhurandhar et al 2008 */
int XLALGetDopplerShiftedFrequencyInfo
  (
   REAL8Vector         *shiftedFreqs, /**< Output list of shifted frequencies */
   UINT4Vector         *lowestBins,   /**< Output list of bin indices */
   REAL8Vector         *kappaValues,  /**< Output list of bin offsets */
   UINT4               numBins,       /**< Number of frequency bins to use */
   PulsarDopplerParams *dopp,         /**< Doppler parameters for signal */
   LIGOTimeGPSVector   *timestamps,   /**< List of SFT timestamps */
   REAL8VectorSequence *vByC          /**< List of dimensionless detector velocities*/
  )
{
  UINT8 numSFTs;
  UINT8 indI;
  REAL8 khat[3];
  UINT4 k;
  REAL8 timeDiff, factor, fhat;

  numSFTs = timestamps->length;
  if ( shiftedFreqs->length !=numSFTs
       || lowestBins->length !=numSFTs
       || kappaValues->length !=numSFTs
       || vByC->length !=numSFTs ) {
    XLALPrintError("Lengths of SFT-indexed lists don't match!");
    XLAL_ERROR(XLAL_EBADLEN );
  }

  if ( vByC->vectorLength != 3 ) {
    XLALPrintError("Velocity vector must be three-dimensional!");
    XLAL_ERROR(XLAL_EBADLEN );
  }

  if ( numBins < 1 ) {
    XLALPrintError("Must specify a positive number of bins to use!");
    XLAL_ERROR(XLAL_EBADLEN );
  }

  if ( dopp->orbit ) {
    XLALPrintError("Doppler shifting for binary orbit not yet implemented!");
    XLAL_ERROR(XLAL_EINVAL);
   
  }

  /* unit vector along wave propagation direction */
  khat[0] = - cos(dopp->Delta) * cos(dopp->Alpha);
  khat[1] = - cos(dopp->Delta) * sin(dopp->Alpha);
  khat[2] = - sin(dopp->Delta);

  /* now calculate the intrinsic signal frequency in the SFT */
  /* fhat = f_0 + f_1(t-t0) + f_2(t-t0)^2/2 + ... */

  /* this is the sft reference time  - the pulsar reference time */
  for (indI=0; indI < numSFTs; indI++) {
    timeDiff = XLALGPSDiff( &(timestamps->data[indI]), &(dopp->refTime));
    fhat = dopp->fkdot[0]; /* initialization */
    factor = 1.0;
    for (k = 1;  k < PULSAR_MAX_SPINS; k++) {
      factor *= timeDiff / k;
      fhat += dopp->fkdot[k] * factor;
    }
    factor = 1.;
    for (k = 0; k < 3; k++) {
      factor += khat[k] * vByC->data[indI*3+k];
    }
    shiftedFreqs->data[indI] = factor * fhat;
    lowestBins->data[indI]
      = ceil(shiftedFreqs->data[indI] * timestamps->deltaT - 0.5*numBins);
    kappaValues->data[indI] = lowestBins->data[indI]
      - shiftedFreqs->data[indI] * timestamps->deltaT;
  }

  return XLAL_SUCCESS;

}
