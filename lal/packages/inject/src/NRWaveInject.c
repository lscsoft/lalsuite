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

/** Spherical Harmonic for the l=2 mode */
COMPLEX16 SphHarm ( 
    UINT4   L,      /**< value of L */
    INT4    M,      /**< value of M */
    REAL4   theta,  /**< angle with respect to the z axis */
    REAL4   phi     /**< angle with respect to the x axis */)

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
		/* Error message informing that the chosen M is incompatible with L*/
		printf("Sorry, the value chosen for m is not compatible with l");
		break;
	}
    }

    else if (L == 3)
    {
	switch ( M )
	{
	    case -3:
		deptheta = Sqrt(21./(2.*LAL_PI))*Cos(theta/2.)*Power(Sin(theta/2.),5);
		out.re = deptheta * cos( -3.0*phi );
		out.im = deptheta * sin( -3.0*phi );
		break;

	    case -2:
		deptheta = Sqrt(7./4.*LAL_PI)*(2 + 3*Cos(theta))*Power(Sin(theta/2.),4);
		out.re = deptheta * cos( -2.0*phi );
		out.im = deptheta * sin( -2.0*phi );
		break;

	    case -1:
		deptheta = Sqrt(35./(2.*LAL_PI))*(Sin(theta) + 4*Sin(2*theta) - 3*Sin(3*theta))/32.;
		out.re = deptheta * cos( -phi );
		out.im = deptheta * sin( -phi );
		break;

	    case 0:
		deptheta = (Sqrt(105./(2.*LAL_PI))*Cos(theta)*Power(Sin(theta),2))/4.;
		out.re = deptheta;
		out.im = deptheta;
		break;

	    case 1:
		deptheta = -Sqrt(35./(2.*LAL_PI))*(Sin(theta) - 4*Sin(2*theta) - 3*Sin(3*theta))/32.;
		out.re = deptheta * cos( phi );
		out.im = deptheta * sin( phi );
		break;
		
	    case 2:
		deptheta = Sqrt(7./LAL_PI)*Power(Cos(theta/2.),4)*(-2 + 3*Cos(theta))/2.;
		out.re = deptheta * cos( 2.0*phi );
		out.im = deptheta * sin( 2.0*phi );
		break;	   
	    
	    case 3:
		deptheta = -Sqrt(21./(2.*LAL_PI))*Power(Cos(theta/2.),5)*Sin(theta/2.);
		out.re = deptheta * cos( 3.0*phi );
		out.im = deptheta * sin( 3.0*phi );
		break;	   

	    default:
		/* Error message informing that the chosen M is incompatible with L*/
		printf("Sorry, the value chosen for m is not compatible with l\n");
		break;
	}
    }    

    else if (L == 4)
    {
	switch ( M )
	{

	    case -4:
		deptheta = 3.*Sqrt(7./LAL_PI)*Power(Cos(theta/2.),2)*Power(Sin(theta/2.),6);
		out.re = deptheta * cos( -4.0*phi );
		out.im = deptheta * sin( -4.0*phi );
		break;

	    case -3:
		deptheta = 3.*Sqrt(7./(2.*LAL_PI))*Cos(theta/2.)*(1 + 2*Cos(theta))*Power(Sin(theta/2.),5);
		out.re = deptheta * cos( -3.0*phi );
		out.im = deptheta * sin( -3.0*phi );
		break;

	    case -2:
		deptheta = (3*(9. + 14.*Cos(theta) + 7.*Cos(2*theta))*Power(Sin(theta/2.),4))/(4.*Sqrt(LAL_PI));
		out.re = deptheta * cos( -2.0*phi );
		out.im = deptheta * sin( -2.0*phi );
		break;

	    case -1:
		deptheta = (3.*(3.*Sin(theta) + 2.*Sin(2*theta) + 7.*Sin(3*theta) - 7.*Sin(4*theta)))/(32.*Sqrt(2*LAL_PI));
		out.re = deptheta * cos( -phi );
		out.im = deptheta * sin( -phi );
		break;

	    case 0:
		deptheta = (3.*Sqrt(5./(2.*LAL_PI))*(5. + 7.*Cos(2*theta))*Power(Sin(theta),2))/16.;
		out.re = deptheta;
		out.im = deptheta;
		break;

	    case 1:
		deptheta = (3.*(3.*Sin(theta) - 2.*Sin(2*theta) + 7.*Sin(3*theta) + 7.*Sin(4*theta)))/(32.*Sqrt(2*LAL_PI));
		out.re = deptheta * cos( phi );
		out.im = deptheta * sin( phi );
		break;
		
	    case 2:
		deptheta = (3.*Power(Cos(theta/2.),4)*(9. - 14.*Cos(theta) + 7.*Cos(2*theta)))/(4.*Sqrt(LAL_PI));
		out.re = deptheta * cos( 2.0*phi );
		out.im = deptheta * sin( 2.0*phi );
		break;	   
	    
	    case 3:
		deptheta = -3.*Sqrt(7./(2.*LAL_PI))*Power(Cos(theta/2.),5)*(-1. + 2.*Cos(theta))*Sin(theta/2.);
		out.re = deptheta * cos( 3.0*phi );
		out.im = deptheta * sin( 3.0*phi );
		break;	   

	    case 4:
		deptheta = 3.*Sqrt(7./LAL_PI)*Power(Cos(theta/2.),6)*Power(Sin(theta/2.),2);
		out.re = deptheta * cos( 4.0*phi );
		out.im = deptheta * sin( 4.0*phi );
		break;	   

	    default:
		/* Error message informing that the chosen M is incompatible with L*/
		printf("Sorry, the value chosen for m is not compatible with l\n");
		break;
	}
    }    

    else if (L == 5)
    {
	switch ( M )
	{

	    case -5:
		deptheta = Sqrt(330./LAL_PI)*Power(Cos(theta/2.),3)*Power(Sin(theta/2.),7);
		out.re = deptheta * cos( -5.0*phi );
		out.im = deptheta * sin( -5.0*phi );
		break;

	    case -4:
		deptheta = Sqrt(33./LAL_PI)*Power(Cos(theta/2.),2)*(2. + 5.*Cos(theta))*Power(Sin(theta/2.),6);
		out.re = deptheta * cos( -4.0*phi );
		out.im = deptheta * sin( -4.0*phi );
		break;

	    case -3:
		deptheta = (Sqrt(33./(2.*LAL_PI))*Cos(theta/2.)*(17. + 24.*Cos(theta) + 15.*Cos(2.*theta))*Power(Sin(theta/2.),5))/4.;
		out.re = deptheta * cos( -3.0*phi );
		out.im = deptheta * sin( -3.0*phi );
		break;

	    case -2:
		deptheta = (Sqrt(11./LAL_PI)*(32. + 57.*Cos(theta) + 36.*Cos(2.*theta) + 15.*Cos(3.*theta))*Power(Sin(theta/2.),4))/8.;
		out.re = deptheta * cos( -2.0*phi );
		out.im = deptheta * sin( -2.0*phi );
		break;

	    case -1:
		deptheta = (Sqrt(77./LAL_PI)*(2.*Sin(theta) + 8.*Sin(2.*theta) + 3.*Sin(3.*theta) + 12.*Sin(4.*theta) - 15.*Sin(5.*theta)))/256.;
		out.re = deptheta * cos( -phi );
		out.im = deptheta * sin( -phi );
		break;

	    case 0:
		deptheta = (Sqrt(1155./(2.*LAL_PI))*(5.*Cos(theta) + 3.*Cos(3.*theta))*Power(Sin(theta),2))/32.;
		out.re = deptheta;
		out.im = deptheta;
		break;

	    case 1:
		deptheta = Sqrt(77./LAL_PI)*(-2.*Sin(theta) + 8.*Sin(2.*theta) - 3.*Sin(3.*theta) + 12.*Sin(4.*theta) + 15.*Sin(5.*theta)))/256.;
		out.re = deptheta * cos( phi );
		out.im = deptheta * sin( phi );
		break;
		
	    case 2:
		deptheta = Sqrt(11./LAL_PI)*Power(Cos(theta/2.),4)*(-32. + 57.*Cos(theta) - 36.*Cos(2.*theta) + 15.*Cos(3.*theta)))/8.;
		out.re = deptheta * cos( 2.0*phi );
		out.im = deptheta * sin( 2.0*phi );
		break;	   
	    
	    case 3:
		deptheta = -Sqrt(33./(2.*LAL_PI))*Power(Cos(theta/2.),5)*(17. - 24.*Cos(theta) + 15.*Cos(2.*theta))*Sin(theta/2.))/4.;
		out.re = deptheta * cos( 3.0*phi );
		out.im = deptheta * sin( 3.0*phi );
		break;	   

	    case 4:
		deptheta = Sqrt(33./LAL_PI)*Power(Cos(theta/2.),6)*(-2. + 5.*Cos(theta))*Power(Sin(theta/2.),2);
		out.re = deptheta * cos( 4.0*phi );
		out.im = deptheta * sin( 4.0*phi );
		break;	   

	    case 5:
		deptheta = -Sqrt(330./LAL_PI)*Power(Cos(theta/2.),7)*Power(Sin(theta/2.),3);
		out.re = deptheta * cos( 5.0*phi );
		out.im = deptheta * sin( 5.0*phi );
		break;	   

	    default:
		/* Error message informing that the chosen M is incompatible with L*/
		printf("Sorry, the value chosen for m is not compatible with l\n");
		break;
	}
    }    


    else 
    {
	/* Error message informing that L!=2 is not yet implemented*/
	fprintf(stderr, "Sorry, for the moment we haven't implemented anything other than l=2");
    }
    
    return( out );
}
