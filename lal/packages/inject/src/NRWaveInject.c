/*
 * Copyright (C) 2006 S.Fairhurst, B. Krishnan, L.Santamaria
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

/** \file NRWaveIO.c
 *  \ingroup NRWaveIO
 *  \author S.Fairhurst, B.Krishnan, L.Santamaria
 * 
 *  \brief Functions for reading/writing numerical relativity waveforms
 *
 * $Id$ 
 *
 */

#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/NRWaveIO.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/Date.h>
#include <lal/LALConstants.h>
#include <lal/NRWaveInject.h>
#include <lal/SphericalHarmonics.h>


/** Takes a (sky averaged) numerical relativity waveform and returns the
 * waveform appropriate for given coalescence phase and inclination angles */
/** for the moment only mode (2,2) implemented */
REAL4TimeVectorSeries *
XLALOrientNRWave( 
    REAL4TimeVectorSeries *strain,         /**< sky average h+, hx data */ 
    UINT4                  modeL,          /**< L                       */
    INT4                   modeM,          /**< M                       */
    REAL4                  inclination,    /**< binary inclination      */
    REAL4                  coa_phase       /**< binary coalescence phase*/)
{
    COMPLEX16  MultSphHarm;
    REAL4      tmp1, tmp2;
    UINT4      vecLength, k;

    vecLength = strain->data->vectorLength;

/* Calculating the (2,2) Spherical Harmonic */
    MultSphHarm = SphHarm( modeL, modeM, inclination, coa_phase );

/* Filling the data vector with the data multiplied by the Harmonic */
    for ( k = 0; k < vecLength; k++)
    {
	tmp1 = strain->data->data[k];
	tmp2 = strain->data->data[vecLength + k];

	strain->data->data[k] = 
	    (tmp1 * MultSphHarm.re) - 
	    (tmp2 * MultSphHarm.im);

	strain->data->data[vecLength + k] = 
	    (tmp2 * MultSphHarm.re) +
	    (tmp1 * MultSphHarm.im);
    }
  return( strain );
}


REAL4TimeSeries *
XLALCalculateNRStrain( REAL4TimeVectorSeries *strain, /**< h+, hx time series data*/
		       SimInspiralTable      *inj,    /**< injection details      */
		       CHAR                  *ifo,    /**< interferometer */
		       INT4            sampleRate     /**< sample rate of time series */)
{
  LALDetector            det;
  double                 fplus;
  double                 fcross;
  double                 tDelay;
  REAL4TimeSeries       *htData = NULL;
  UINT4                  vecLength, k;
  InterferometerNumber   ifoNumber = LAL_UNKNOWN_IFO;

  /* get the detector information */
  memset( &det, 0, sizeof(LALDetector) );
  ifoNumber = XLALIFONumber( ifo );
  XLALReturnDetector( &det, ifoNumber );

  /* compute detector response */
  XLALComputeDetAMResponse(&fplus, &fcross, det.response, inj->longitude, 
      inj->latitude, inj->polarization, inj->end_time_gmst);

  /* calculate the time delay */
  tDelay = XLALTimeDelayFromEarthCenter( det.location, inj->longitude,
      inj->latitude, &(inj->geocent_end_time) );

  /* create htData */
  htData = LALCalloc(1, sizeof(*htData));
  if (!htData) 
  {
    XLAL_ERROR_NULL( "XLALCalculateNRStrain", XLAL_ENOMEM );
  }
  vecLength = strain->data->vectorLength;
  htData->data = XLALCreateREAL4Vector( vecLength );
  if ( ! htData->data )
  {
    XLAL_ERROR_NULL( "XLALCalculateNRStrain", XLAL_ENOMEM );
  }

  /* store the htData */
  htData->epoch = *XLALGPSAdd( &(inj->geocent_end_time), tDelay );
  htData->deltaT = strain->deltaT;

  for ( k = 0; k < vecLength; ++k )
  {
    htData->data->data[k] = fplus * strain->data->data[k]  + 
      fcross * strain->data->data[vecLength + k];
  }

  /*interpolate to given sample rate */
  htData = XLALInterpolateNRWave( htData, sampleRate);

  return( htData );
}



/** Function for interpolating time series to a given sampling rate. 
    Input vector is destroyed and a new vector is allocated.
*/
REAL4TimeSeries *
XLALInterpolateNRWave( REAL4TimeSeries *in,           /**< input strain time series */
		       INT4            sampleRate     /**< sample rate of time series */)
{

  REAL4TimeSeries *ret=NULL;
  REAL8 deltaTin, deltaTout, r, y1, y2;
  REAL8 tObs; /* duration of signal */
  UINT4 k, lo, numPoints;
  
  deltaTin = in->deltaT;
  tObs = deltaTin * in->data->length;  

  /* length of output vector */ 
  numPoints = (UINT4) (sampleRate * tObs);


  /* allocate memory */
  ret = LALCalloc(1, sizeof(*ret));
  if (!ret) 
  {
    XLAL_ERROR_NULL( "XLALCalculateNRStrain", XLAL_ENOMEM );
  }

  ret->data = XLALCreateREAL4Vector( numPoints );
  if (! ret->data) 
  {
    XLAL_ERROR_NULL( "XLALCalculateNRStrain", XLAL_ENOMEM );
  }

  ret->deltaT = 1./sampleRate;
  deltaTout = ret->deltaT;

  /* copy stuff from in which should be the same */
  ret->epoch = in->epoch;
  ret->f0 = in->f0;
  ret->sampleUnits = in->sampleUnits;
  strcpy(ret->name, in->name);

  /* go over points of output vector and interpolate linearly
     using closest points of input */
  for (k = 0; k < numPoints; k++) {

    lo = (UINT4)( k*deltaTout / deltaTin);

    /* y1 and y2 are the input values at x1 and x2 */
    /* here we need to make sure that we don't exceed
       bounds of input vector */
    y1 = in->data->data[lo];
    y2 = in->data->data[lo+1];

    /* we want to calculate y2*r + y1*(1-r) where
       r = (x-x1)/(x2-x1) */
    r = k*deltaTout / deltaTin - lo;
    
    ret->data->data[k] = y2 * r + y1 * (1 - r);
  }

  /* destroy input vector */
  XLALDestroyREAL4Vector ( in->data);
  LALFree(in);

  return ret;
}



/** For given inspiral parameters, find nearest waveform in 
    catalog of numerical relativity waveforms.  At the moment, only
    the mass ratio is considered.  
*/
CHAR *
XLALFindNRFile( NRWaveCatalog *nrCatalog,   /**< input  NR wave catalog  */
		SimInspiralTable      *inj, /**< injection details  */
		INT4  modeL,                /**< mode index l*/
		INT4  modeM                 /**< mode index m*/)
{

  REAL8 massRatioIn, massRatio, diff, newDiff;
  UINT4 k, best;
  CHAR *ret=NULL;

  massRatioIn = inj->mass1/inj->mass2;
  massRatio = nrCatalog->data[0].massRatio;

  diff = fabs(massRatio - massRatioIn);

  /* look over catalog and fimd waveform closest in mass ratio */
  for (best = 0, k = 0; k < nrCatalog->length; k++) {

    massRatio = nrCatalog->data[k].massRatio;
    newDiff = fabs(massRatio - massRatioIn);

    if ( diff > newDiff) {
      diff = newDiff;
      best = k;
    }    
  } 

  ret = LALMalloc( LALNameLength*sizeof(CHAR) );

  strcpy(ret, nrCatalog->data[best].filename);

  return ret;

}
