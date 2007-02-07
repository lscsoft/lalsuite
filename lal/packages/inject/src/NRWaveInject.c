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


/** Spherical Harmonic for the l=2 mode */
/* This function is for the moment not used in the code */
/* It'll be useful only when XLALOrientNRWave is allowed to take L, M as input */
COMPLEX16 SphHarm ( 
    UINT4   L,      /* value of L */
    INT4    M,      /* value of M */
    REAL4   theta,  /* angle with respect to the z axis */
    REAL4   phi    /* angle with respect to the x axis */)

{
    COMPLEX16  out; /* complex number */
    REAL4      deptheta; /** dependency on theta */

    if (L == 2)
    {
	switch ( M )
	{
	    case -2:
		deptheta = sqrt( 5.0 / ( 64.0 * LAL_PI ) ) * ( 1.0 - cos( theta ))*( 1.0 - cos( theta ));
		out.re = deptheta * cos( -2.0*phi );
		out.im = deptheta * sin( -2.0*phi );
		break;

	    case -1:
		deptheta = sqrt( 5.0 / ( 16.0 * LAL_PI ) ) * sin( theta )*( 1.0 - cos( theta ));
		out.re = deptheta * cos( -phi );
		out.im = deptheta * sin( -phi );
		break;

	    case 0:
		deptheta = sqrt( 15.0 / ( 32.0 * LAL_PI ) ) * sin( theta )*sin( theta );
		out.re = deptheta;
		out.im = deptheta;
		break;

	    case 1:
		deptheta = sqrt( 5.0 / ( 16.0 * LAL_PI ) ) * sin( theta )*( 1.0 + cos( theta ));
		out.re = deptheta * cos( phi );
		out.im = deptheta * sin( phi );
		break;
		
	    case 2:
		deptheta = sqrt( 5.0 / ( 64.0 * LAL_PI ) ) * ( 1.0 + cos( theta ))*( 1.0 + cos( theta ));
		out.re = deptheta * cos( 2.0*phi );
		out.im = deptheta * sin( 2.0*phi );
		break;	   
	    
	    default:
		/* Error message informing that the chosen M is incompatible with L=2*/
		break;
	}
    }
    else 
    {
	/* Error message informing that L!=2 is not yet implemented*/
    }
    
    return( out );
}


/** Spherical Harmonic for the l=2, m=2 mode */
/* Eventually all these functions should be placed in their own library */
COMPLEX16 SphHarm22 ( 
    REAL4   theta,  /* angle with respect to the z axis */
    REAL4   phi    /* angle with respect to the x axis */)

{
    COMPLEX16  out; /* complex number */
    REAL4      deptheta; /** dependency on theta */

    deptheta = sqrt( 5.0 / ( 64.0 * LAL_PI ) ) * ( 1.0 + cos( theta ))*( 1.0 + cos( theta ));

    out.re = deptheta * cos( 2.0*phi );
    out.im = deptheta * sin( 2.0*phi );

    return( out );
}


/** Takes a (sky averaged) numerical relativity waveform and returns the
 * waveform appropriate for given coalescence phase and inclination angles */
/** for the moment only mode (2,2) implemented */
REAL4TimeVectorSeries *
XLALOrientNRWave( 
    REAL4TimeVectorSeries *strain,         /**< sky average h+, hx data */ 
    REAL4                  inclination,    /**< binary inclination      */
    REAL4                  coa_phase       /**< binary coalescence phase*/)
{
    COMPLEX16  MultSphHarm;
    REAL4      tmp1, tmp2;
    UINT4      vecLength, k;

    vecLength = strain->data->vectorLength;

/* Calculating the (2,2) Spherical Harmonic */
    MultSphHarm = SphHarm22( inclination, coa_phase );

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
XLALCalculateNRStrain( 
    REAL4TimeVectorSeries *strain, /**< h+, hx time series data*/
    SimInspiralTable      *inj,    /**< injection details      */
    CHAR                  *ifo,    /**< interferometer */
    INT4            sampleRate     /**< sample rate of time series */)
{
  LALDetector            det;
  double                 fplus;
  double                 fcross;
  double                 tDelay;
  REAL4TimeSeries       *htData = NULL;
  int                    k;

  memset( &det, 0, sizeof(LALDetector) );
  ifoNumber = XLALIFONumber( ifo );
  XLALReturnDetector( &det, ifoNumber );

  XLALComputeDetAMResponse(&fplus, &fcross, det.response, inj->longitude, 
      inj->latitude, inj->polarization, inj->end_time_gmst);

  tDelay = XLALTimeDelayFromEarthCenter( det.location, inj->longitude,
      inj->latitude, &(inj->geocent_end_time) );

  /* store the ht data */
  htData = LALCalloc(1, sizeof(*htData));
  if (!htData) 
  {
    XLAL_ERROR_NULL( "XLALCalculateNRStrain", XLAL_ENOMEM );
  }
  htData->epoch = *XLALGPSAdd( &(inj->geocent_end_time), tDelay );
  htData->deltaT = strain->deltaT;
  htData->data = XLALCreateREAL4Vector( strain->data->vectorLength );
  if ( ! htData->data )
  {
    XLAL_ERROR_NULL( "XLALCalculateNRStrain", XLAL_ENOMEM );
  }

  for ( k = 0; k < strain->data->vectorLength; ++k )
  {
    htData->data->data[k] = fplus * strain->data->data[k]  + 
      fcross * strain->data->data[strain->data->vectorLength + k];
  }

  /* XXX TO DO: interpolate data to the required sample rate XXX */

  return( htData );
}

