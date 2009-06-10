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

/** \file NRWaveInject.c
 *  \ingroup NRWaveInject
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
#include <lal/LIGOMetadataUtils.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/LALConstants.h>
#include <lal/NRWaveInject.h>
#include <lal/Random.h>
#include <lal/Inject.h>
#include <lal/LALSimulation.h>
#include <lal/LALDetectors.h>
#include <lal/LALComplex.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/SkyCoordinates.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Window.h>

#include <gsl/gsl_heapsort.h>

/* macro to "use" unused function parameters */
#define UNUSED(expr) do { (void)(expr); } while (0)

NRCSID (NRWAVEINJECTC, "$Id$");

int compare_abs_float(const void *a, const void *b);
int compare_abs_double(const void *a, const void *b);

/** Takes a strain of h+ and hx data and stores it in a temporal
 *  strain in order to perform the sum over l and m modes **/
REAL4TimeVectorSeries *
XLALSumStrain(
    REAL4TimeVectorSeries *tempstrain,     /**< storing variable */
    REAL4TimeVectorSeries *strain          /**< variable to add  */)
{
    UINT4      vecLength, length, k;

    vecLength = strain->data->vectorLength;
    length = strain->data->length;

    for ( k = 0; k < vecLength*length; k++)
      {
	tempstrain->data->data[k] += strain->data->data[k];
      }
    return( tempstrain );
}

/* REAL8 version */
REAL8TimeVectorSeries *
XLALSumStrainREAL8(
    REAL8TimeVectorSeries *tempstrain,     /**< storing variable */
    REAL8TimeVectorSeries *strain          /**< variable to add  */)
{
    UINT4      vecLength, length, k;

    vecLength = strain->data->vectorLength;
    length = strain->data->length;

    for ( k = 0; k < vecLength*length; k++)
      {
	tempstrain->data->data[k] += strain->data->data[k];
      }
    return( tempstrain );
}


/** Takes a (sky averaged) numerical relativity waveform and returns the
 * waveform appropriate for given coalescence phase and inclination angles */
/* REAL4TimeVectorSeries */
INT4
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
    /* need some error checking */
    XLALSphHarm( &MultSphHarm, modeL, modeM, inclination, coa_phase );

    /* Filling the data vector with the data multiplied by the Harmonic */
    for ( k = 0; k < vecLength; k++)
    {
	tmp1 = strain->data->data[k];
	tmp2 = strain->data->data[vecLength + k];

	strain->data->data[k] =
	    (tmp1 * MultSphHarm.re) +
	    (tmp2 * MultSphHarm.im);

	strain->data->data[vecLength + k] =
	    (tmp2 * MultSphHarm.re) -
	    (tmp1 * MultSphHarm.im);
    }

    return 0;
/*   return( strain ); */
}


/** Takes a (sky averaged) numerical relativity waveform and returns the
 * waveform appropriate for given coalescence phase and inclination angles */
REAL8TimeVectorSeries *
XLALOrientNRWaveREAL8(
    REAL8TimeVectorSeries *strain,         /**< sky average h+, hx data */
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
    /* need some error checking */
    XLALSphHarm( &MultSphHarm, modeL, modeM, inclination, coa_phase );

    /* Filling the data vector with the data multiplied by the Harmonic */
    for ( k = 0; k < vecLength; k++)
    {
	tmp1 = strain->data->data[k];
	tmp2 = strain->data->data[vecLength + k];

	strain->data->data[k] =
	    (tmp1 * MultSphHarm.re) +
	    (tmp2 * MultSphHarm.im);

	strain->data->data[vecLength + k] =
	    (tmp2 * MultSphHarm.re) -
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

  if( ifoNumber == LAL_UNKNOWN_IFO )
  {
    XLALPrintWarning( "Unknown ifo! Injecting signal overhead!\n" );
    fplus  = 1.0;
    fcross = 0.0;
    tDelay = 0.0;
  }
  else
  {
    /* compute detector response */
    XLALComputeDetAMResponse(&fplus, &fcross, det.response, inj->longitude,
          inj->latitude, inj->polarization, inj->end_time_gmst);

    /* calculate the time delay */
    tDelay = XLALTimeDelayFromEarthCenter( det.location, inj->longitude,
        inj->latitude, &(inj->geocent_end_time) );
  }

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
    LALFree( htData );
    XLAL_ERROR_NULL( "XLALCalculateNRStrain", XLAL_ENOMEM );
  }

  /* store the htData */

  /* add timedelay to inj->geocent_end_time */
  {
    REAL8 tEnd;
    tEnd =  XLALGPSGetREAL8(&(inj->geocent_end_time) );
    tEnd += tDelay;
    XLALGPSSetREAL8(&(htData->epoch), tEnd );
  }
  /* Using XLALGPSAdd is bad because that changes the GPS time! */
  /*  htData->epoch = *XLALGPSAdd( &(inj->geocent_end_time), tDelay ); */

  htData->deltaT = strain->deltaT;
  htData->sampleUnits = lalADCCountUnit;

  for ( k = 0; k < vecLength; ++k )
  {
    htData->data->data[k] = (fplus * strain->data->data[k]  +
			     fcross * strain->data->data[vecLength + k]) / inj->distance;
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
  REAL8 deltaTin, deltaTout, r, y_1, y_2;
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

    /* y_1 and y_2 are the input values at x1 and x2 */
    /* here we need to make sure that we don't exceed
       bounds of input vector */
    if ( lo < in->data->length - 1) {
      y_1 = in->data->data[lo];
      y_2 = in->data->data[lo+1];

      /* we want to calculate y_2*r + y_1*(1-r) where
	 r = (x-x1)/(x2-x1) */
      r = k*deltaTout / deltaTin - lo;

      ret->data->data[k] = y_2 * r + y_1 * (1 - r);
    }
    else {
      ret->data->data[k] = 0.0;
    }
  }

  /* destroy input vector */
  XLALDestroyREAL4Vector ( in->data);
  LALFree(in);

  return ret;
}



int compare_abs_float(const void *a, const void *b){

  const REAL4 *af, *bf;

  af = (const REAL4 *)a;
  bf = (const REAL4 *)b;

  if ( fabs(*af) > fabs(*bf))
    return 1;
  else if  ( fabs(*af) < fabs(*bf))
    return -1;
  else
    return 0;
}


int compare_abs_double(const void *a, const void *b){

  const REAL8 *af, *bf;

  af = (const REAL8 *)a;
  bf = (const REAL8 *)b;

  if ( fabs(*af) > fabs(*bf))
    return 1;
  else if  ( fabs(*af) < fabs(*bf))
    return -1;
  else
    return 0;
}


/** Function for calculating the coalescence time (defined to be the peak) of a NR wave */
INT4
XLALFindNRCoalescenceTime(REAL8 *tc,
			  const REAL4TimeVectorSeries *in   /**< input strain time series */)
{

  size_t *ind=NULL;
  size_t len;
  REAL4 *sumSquare=NULL;
  UINT4 k;

  len = in->data->vectorLength;
  ind = LALCalloc(1,len*sizeof(*ind));

  sumSquare = LALCalloc(1, len*sizeof(*sumSquare));

  for (k=0; k < len; k++) {
    sumSquare[k] = in->data->data[k]*in->data->data[k] +
      in->data->data[k + len]*in->data->data[k + len];
  }

  gsl_heapsort_index( ind, sumSquare, len, sizeof(REAL4), compare_abs_float);

  *tc = ind[len-1] * in->deltaT;

  LALFree(ind);
  LALFree(sumSquare);

  return 0;
}


/** Function for calculating the coalescence time (defined to be the
    peak) of a NR wave
    This uses the peak of h(t)
*/
INT4
XLALFindNRCoalescenceTimeFromhoft(REAL8 *tc,
			  const REAL4TimeSeries *in   /**< input strain time series */)
{

  size_t *ind=NULL;
  size_t len;

  len = in->data->length;
  ind = LALCalloc(1,len*sizeof(*ind));

/*   gsl_heapsort_index( ind, in->data->data, len, sizeof(REAL4), compare_abs_float); */

  *tc = ind[len-1] * in->deltaT;

  LALFree(ind);

  return 0;
}


/** Function for calculating the coalescence time (defined to be the peak) of a NR wave */
INT4
XLALFindNRCoalescenceTimeREAL8(REAL8 *tc,
			       const REAL8TimeSeries *in   /**< input strain time series */)
{

  size_t *ind=NULL;
  size_t len;

  len = in->data->length;
  ind = LALCalloc(len, sizeof(*ind));

  gsl_heapsort_index( ind, in->data->data, len, sizeof(REAL8), compare_abs_double);

  *tc = ind[len-1] * in->deltaT;

  LALFree(ind);

  return 0;
}



/** For given inspiral parameters, find nearest waveform in
    catalog of numerical relativity waveforms.  At the moment, only
    the mass ratio is considered.
*/
INT4
XLALFindNRFile( NRWaveMetaData   *out,       /**< output wave data */
	        NRWaveCatalog    *nrCatalog, /**< input  NR wave catalog  */
		const SimInspiralTable *inj,       /**< injection details  */
		INT4  modeL,                 /**< mode index l*/
		INT4  modeM                  /**< mode index m*/)
{

  REAL8 massRatioIn, massRatio, diff, newDiff;
  UINT4 k, best=0;

  /* check arguments are sensible */
  if ( !out ) {
    LALPrintError ("\nOutput pointer is NULL !\n\n");
    XLAL_ERROR ( "XLALFindNRFile", XLAL_EINVAL);
  }

  if ( !nrCatalog ) {
    LALPrintError ("\n NR Catalog pointer is NULL !\n\n");
    XLAL_ERROR ( "XLALFindNRFile", XLAL_EINVAL);
  }

  if ( !inj ) {
    LALPrintError ("\n SimInspiralTable pointer is NULL !\n\n");
    XLAL_ERROR ( "XLALFindNRFile", XLAL_EINVAL);
  }

  /* we want to check for the highest mass before calculating the mass ratio */
  if (inj->mass2 > inj->mass1) {
     massRatioIn = inj->mass2/inj->mass1;
  }
  else {
     massRatioIn = inj->mass1/inj->mass2;
  }

  /*   massRatio = nrCatalog->data[0].massRatio; */

  /*   diff = fabs(massRatio - massRatioIn); */

  /* look over catalog and fimd waveform closest in mass ratio */
  /* initialize diff */
  diff = -1;
  for (k = 0; k < nrCatalog->length; k++) {

    /* check catalog data is not null */
    if ( nrCatalog->data + k == NULL ) {
      LALPrintError ("\n NR Catalog data is NULL !\n\n");
      XLAL_ERROR ( "XLALFindNRFile", XLAL_EINVAL);
    }

    /* look for waveforms with correct mode values */
    if ((modeL == nrCatalog->data[k].mode[0]) && (modeM == nrCatalog->data[k].mode[1])) {

      massRatio = nrCatalog->data[k].massRatio;
      newDiff = fabs(massRatio - massRatioIn);

      if ( (diff < 0) || (diff > newDiff)) {
	diff = newDiff;
	best = k;
      }
    } /* if (modeL == ...) */
  } /* loop over waveform catalog */

  /* error checking if waveforms with input mode values were not found */
  if ( diff < 0) {
    LALPrintError ("\n Input mode numbers not found !\n\n");
    XLAL_ERROR ( "XLALFindNRFile", XLAL_EINVAL);
  }

  /* copy best match to output */
  memcpy(out, nrCatalog->data + best, sizeof(NRWaveMetaData));

  return 0;

}

/** Spin 2 weighted spherical Harmonic */
INT4 XLALSphHarm ( COMPLEX16 *out, /**< output */
		   UINT4   L,      /**< value of L */
		   INT4    M,      /**< value of M */
		   REAL4   theta,  /**< angle with respect to the z axis */
		   REAL4   phi     /**< angle with respect to the x axis */)

{
    REAL4      deptheta; /** dependency on theta */

  /* check arguments are sensible */
    if ( !out ) {
      LALPrintError ("\nOutput pointer is NULL !\n\n");
      XLAL_ERROR ( "XLALSphHarm", XLAL_EINVAL);
    }

    if (L == 2)
    {
      switch ( M ){
      case -2:
	deptheta = sqrt( 5.0 / ( 64.0 * LAL_PI ) ) * ( 1.0 - cos( theta ))*( 1.0 - cos( theta ));
	out->re = deptheta * cos( -2.0*phi );
	out->im = deptheta * sin( -2.0*phi );
	break;

      case -1:
	deptheta = sqrt( 5.0 / ( 16.0 * LAL_PI ) ) * sin( theta )*( 1.0 - cos( theta ));
	out->re = deptheta * cos( -phi );
	out->im = deptheta * sin( -phi );
	break;

      case 0:
	deptheta = sqrt( 15.0 / ( 32.0 * LAL_PI ) ) * sin( theta )*sin( theta );
	out->re = deptheta;
	out->im = 0.0;
	break;

      case 1:
	deptheta = sqrt( 5.0 / ( 16.0 * LAL_PI ) ) * sin( theta )*( 1.0 + cos( theta ));
	out->re = deptheta * cos( phi );
	out->im = deptheta * sin( phi );
	break;

      case 2:
	deptheta = sqrt( 5.0 / ( 64.0 * LAL_PI ) ) * ( 1.0 + cos( theta ))*( 1.0 + cos( theta ));
	out->re = deptheta * cos( 2.0*phi );
	out->im = deptheta * sin( 2.0*phi );
	break;

      default:
	/* Error message informing that the chosen M is incompatible with L*/
	LALPrintError ("\n Inconsistent (L, M) values \n\n");
	XLAL_ERROR ( "XLALSphHarm", XLAL_EINVAL);
	break;
      }  /* switch (M) */
    }  /* L==2*/

    else if (L == 3)
      {
	switch ( M ) {
	case -3:
	  deptheta = sqrt(21./(2.*LAL_PI))*cos(theta/2.)*pow(sin(theta/2.),5);
	  out->re = deptheta * cos( -3.0*phi );
	  out->im = deptheta * sin( -3.0*phi );
	  break;

	case -2:
	  deptheta = sqrt(7./(4.*LAL_PI))*(2 + 3*cos(theta))*pow(sin(theta/2.),4);
	  out->re = deptheta * cos( -2.0*phi );
	  out->im = deptheta * sin( -2.0*phi );
	  break;

	case -1:
	  deptheta = sqrt(35./(2.*LAL_PI))*(sin(theta) + 4*sin(2*theta) - 3*sin(3*theta))/32.;
	  out->re = deptheta * cos( -phi );
	  out->im = deptheta * sin( -phi );
	  break;

	case 0:
	  deptheta = (sqrt(105./(2.*LAL_PI))*cos(theta)*pow(sin(theta),2))/4.;
	  out->re = deptheta;
	  out->im = 0.0;
	  break;

	case 1:
	  deptheta = -sqrt(35./(2.*LAL_PI))*(sin(theta) - 4*sin(2*theta) - 3*sin(3*theta))/32.;
	  out->re = deptheta * cos( phi );
	  out->im = deptheta * sin( phi );
	  break;

	case 2:
	  deptheta = sqrt(7./LAL_PI)*pow(cos(theta/2.),4)*(-2 + 3*cos(theta))/2.;
	  out->re = deptheta * cos( 2.0*phi );
	  out->im = deptheta * sin( 2.0*phi );
	  break;

	case 3:
	  deptheta = -sqrt(21./(2.*LAL_PI))*pow(cos(theta/2.),5)*sin(theta/2.);
	  out->re = deptheta * cos( 3.0*phi );
	  out->im = deptheta * sin( 3.0*phi );
	  break;

	default:
	  /* Error message informing that the chosen M is incompatible with L*/
	  LALPrintError ("\n Inconsistent (L, M) values \n\n");
	  XLAL_ERROR ( "XLALSphHarm", XLAL_EINVAL);
	  break;
	}
      }   /* L==3 */

    else if (L == 4)
    {
      switch ( M )	{

      case -4:
	deptheta = 3.*sqrt(7./LAL_PI)*pow(cos(theta/2.),2)*pow(sin(theta/2.),6);
	out->re = deptheta * cos( -4.0*phi );
	out->im = deptheta * sin( -4.0*phi );
	break;

      case -3:
	deptheta = 3.*sqrt(7./(2.*LAL_PI))*cos(theta/2.)*(1 + 2*cos(theta))*pow(sin(theta/2.),5);
	out->re = deptheta * cos( -3.0*phi );
	out->im = deptheta * sin( -3.0*phi );
	break;

      case -2:
	deptheta = (3*(9. + 14.*cos(theta) + 7.*cos(2*theta))*pow(sin(theta/2.),4))/(4.*sqrt(LAL_PI));
	out->re = deptheta * cos( -2.0*phi );
	out->im = deptheta * sin( -2.0*phi );
	break;

      case -1:
	deptheta = (3.*(3.*sin(theta) + 2.*sin(2*theta) + 7.*sin(3*theta) - 7.*sin(4*theta)))/(32.*sqrt(2*LAL_PI));
	out->re = deptheta * cos( -phi );
	out->im = deptheta * sin( -phi );
	break;

      case 0:
	deptheta = (3.*sqrt(5./(2.*LAL_PI))*(5. + 7.*cos(2*theta))*pow(sin(theta),2))/16.;
	out->re = deptheta;
	out->im = 0.0;
	break;

      case 1:
	deptheta = (3.*(3.*sin(theta) - 2.*sin(2*theta) + 7.*sin(3*theta) + 7.*sin(4*theta)))/(32.*sqrt(2*LAL_PI));
	out->re = deptheta * cos( phi );
	out->im = deptheta * sin( phi );
	break;

      case 2:
	deptheta = (3.*pow(cos(theta/2.),4)*(9. - 14.*cos(theta) + 7.*cos(2*theta)))/(4.*sqrt(LAL_PI));
	out->re = deptheta * cos( 2.0*phi );
	out->im = deptheta * sin( 2.0*phi );
	break;

      case 3:
	deptheta = -3.*sqrt(7./(2.*LAL_PI))*pow(cos(theta/2.),5)*(-1. + 2.*cos(theta))*sin(theta/2.);
	out->re = deptheta * cos( 3.0*phi );
	out->im = deptheta * sin( 3.0*phi );
	break;

      case 4:
	deptheta = 3.*sqrt(7./LAL_PI)*pow(cos(theta/2.),6)*pow(sin(theta/2.),2);
	out->re = deptheta * cos( 4.0*phi );
	out->im = deptheta * sin( 4.0*phi );
	break;

      default:
	/* Error message informing that the chosen M is incompatible with L*/
	LALPrintError ("\n Inconsistent (L, M) values \n\n");
	XLAL_ERROR ( "XLALSphHarm", XLAL_EINVAL);
	break;
      }
    }    /* L==4 */

    else if (L == 5)
    {
	switch ( M )	{

	case -5:
	  deptheta = sqrt(330./LAL_PI)*pow(cos(theta/2.),3)*pow(sin(theta/2.),7);
	  out->re = deptheta * cos( -5.0*phi );
	  out->im = deptheta * sin( -5.0*phi );
	  break;

	case -4:
	  deptheta = sqrt(33./LAL_PI)*pow(cos(theta/2.),2)*(2. + 5.*cos(theta))*pow(sin(theta/2.),6);
	  out->re = deptheta * cos( -4.0*phi );
	  out->im = deptheta * sin( -4.0*phi );
	  break;

	case -3:
	  deptheta = (sqrt(33./(2.*LAL_PI))*cos(theta/2.)*(17. + 24.*cos(theta) + 15.*cos(2.*theta))*pow(sin(theta/2.),5))/4.;
	  out->re = deptheta * cos( -3.0*phi );
	  out->im = deptheta * sin( -3.0*phi );
	  break;

	case -2:
	  deptheta = (sqrt(11./LAL_PI)*(32. + 57.*cos(theta) + 36.*cos(2.*theta) + 15.*cos(3.*theta))*pow(sin(theta/2.),4))/8.;
	  out->re = deptheta * cos( -2.0*phi );
	  out->im = deptheta * sin( -2.0*phi );
	  break;

	case -1:
	  deptheta = (sqrt(77./LAL_PI)*(2.*sin(theta) + 8.*sin(2.*theta) + 3.*sin(3.*theta) + 12.*sin(4.*theta) - 15.*sin(5.*theta)))/256.;
	  out->re = deptheta * cos( -phi );
	  out->im = deptheta * sin( -phi );
	  break;

	case 0:
	  deptheta = (sqrt(1155./(2.*LAL_PI))*(5.*cos(theta) + 3.*cos(3.*theta))*pow(sin(theta),2))/32.;
	  out->re = deptheta;
	  out->im = 0.0;
	  break;

	case 1:
	  deptheta = sqrt(77./LAL_PI)*(-2.*sin(theta) + 8.*sin(2.*theta) - 3.*sin(3.*theta) + 12.*sin(4.*theta) + 15.*sin(5.*theta))/256.;
	  out->re = deptheta * cos( phi );
	  out->im = deptheta * sin( phi );
	  break;

	case 2:
	  deptheta = sqrt(11./LAL_PI)*pow(cos(theta/2.),4)*(-32. + 57.*cos(theta) - 36.*cos(2.*theta) + 15.*cos(3.*theta))/8.;
	  out->re = deptheta * cos( 2.0*phi );
	  out->im = deptheta * sin( 2.0*phi );
	  break;

	case 3:
	  deptheta = -sqrt(33./(2.*LAL_PI))*pow(cos(theta/2.),5)*(17. - 24.*cos(theta) + 15.*cos(2.*theta))*sin(theta/2.)/4.;
	  out->re = deptheta * cos( 3.0*phi );
	  out->im = deptheta * sin( 3.0*phi );
	  break;

	case 4:
	  deptheta = sqrt(33./LAL_PI)*pow(cos(theta/2.),6)*(-2. + 5.*cos(theta))*pow(sin(theta/2.),2);
	  out->re = deptheta * cos( 4.0*phi );
	  out->im = deptheta * sin( 4.0*phi );
	  break;

	case 5:
	  deptheta = -sqrt(330./LAL_PI)*pow(cos(theta/2.),7)*pow(sin(theta/2.),3);
	  out->re = deptheta * cos( 5.0*phi );
	  out->im = deptheta * sin( 5.0*phi );
	  break;

	default:
	  /* Error message informing that the chosen M is incompatible with L*/
	  LALPrintError ("\n Inconsistent (L, M) values \n\n");
	  XLAL_ERROR ( "XLALSphHarm", XLAL_EINVAL);
	  break;
	}
    }  /* L==5 */


    else
    {
      /* Error message informing that L!=2 is not yet implemented*/
      LALPrintError ("\n These (L, M) values not implemented yet \n\n");
      XLAL_ERROR ( "XLALSphHarm", XLAL_EINVAL);
    }

    return 0;
}


void LALInjectStrainGW( LALStatus                 *status,
			REAL4TimeSeries           *injData,
			REAL4TimeVectorSeries     *strain,
			SimInspiralTable          *thisInj,
			CHAR                      *ifo,
			REAL8                     dynRange)
{

  REAL8 sampleRate;
  REAL4TimeSeries *htData = NULL;
  UINT4  k;
  REAL8 offset;
  InspiralApplyTaper taper = INSPIRAL_TAPER_NONE;

  INITSTATUS (status, "LALNRInject",  NRWAVEINJECTC);
  ATTATCHSTATUSPTR (status);

  /* sampleRate = 1.0/strain->deltaT;   */
  /* use the sample rate required for the output time series */
  sampleRate = 1.0/injData->deltaT;

  /*compute strain for given sky location*/
  htData = XLALCalculateNRStrain( strain, thisInj, ifo, sampleRate );
  if ( !htData )
  {
    ABORTXLAL( status );
  }

  /* multiply the input data by dynRange */
  for ( k = 0 ; k < htData->data->length ; ++k )
  {
    htData->data->data[k] *= dynRange;
  }

  XLALFindNRCoalescenceTime( &offset, strain);

  XLALGPSAdd( &(htData->epoch), -offset);


  /* Taper the signal if required */
  if ( strcmp( thisInj->taper, "TAPER_NONE" ) )
  {

    if ( ! strcmp( "TAPER_START", thisInj->taper ) )
    {
      taper = INSPIRAL_TAPER_START;
    }
    else if (  ! strcmp( "TAPER_END", thisInj->taper ) )
    {
      taper = INSPIRAL_TAPER_END;
    }
    else if (  ! strcmp( "TAPER_STARTEND", thisInj->taper ) )
    {
      taper = INSPIRAL_TAPER_STARTEND;
    }
    else
    {
      XLALPrintError( "Unsupported tapering type specified: %s\n", thisInj->taper );
      XLALDestroyREAL4Vector ( htData->data);
      LALFree(htData);
      ABORT( status, NRWAVEINJECT_EVAL, NRWAVEINJECT_MSGEVAL );
    }
    if ( XLALInspiralWaveTaper( htData->data, taper ) == XLAL_FAILURE )
    {
      XLALClearErrno();
      XLALDestroyREAL4Vector ( htData->data);
      LALFree(htData);
      ABORTXLAL( status );
    }
  }

  /* Band-passing probably not as important for NR-type waveforms */
  /* TODO: Implement if required at a later date */
  if ( thisInj->bandpass )
  {
    XLALPrintError( "Band-passing not yet implemented for InjectStrainGW.\n" );
    XLALDestroyREAL4Vector ( htData->data);
    LALFree(htData);
    ABORTXLAL( status );
  }

  /* inject the htData into injection time stream */
  LALSSInjectTimeSeries( status->statusPtr, injData, htData );
  BEGINFAIL( status )
  {
    XLALDestroyREAL4Vector ( htData->data);
    LALFree( htData );
  }
  ENDFAIL( status );

  /* set channel name */
  LALSnprintf( injData->name, LIGOMETA_CHANNEL_MAX * sizeof( CHAR ),
    "%s:STRAIN", ifo );

  XLALDestroyREAL4Vector ( htData->data);
  LALFree(htData);

  DETATCHSTATUSPTR(status);
  RETURN(status);

}


/** REAL8 version of above but using Jolien's new functions */
void LALInjectStrainGWREAL8( LALStatus                 *status,
			     REAL8TimeSeries           *injData,
			     REAL8TimeVectorSeries     *strain,
			     SimInspiralTable          *thisInj,
			     CHAR                      *ifo,
			     REAL8                     dynRange)
{

  REAL8 sampleRate;
  REAL8TimeSeries *htData = NULL;
  REAL8TimeSeries *hplus = NULL;
  REAL8TimeSeries *hcross = NULL;
  UINT4  k, len;
  REAL8 offset;
  InterferometerNumber  ifoNumber = LAL_UNKNOWN_IFO;
  LALDetector det;

  INITSTATUS (status, "LALNRInject",  NRWAVEINJECTC);
  ATTATCHSTATUSPTR (status);

  /* dynRange is unused in this function */
  UNUSED(dynRange);

  /* get the detector information */
  memset( &det, 0, sizeof(LALDetector) );
  ifoNumber = XLALIFONumber( ifo );
  XLALReturnDetector( &det, ifoNumber );

  /* Band-passing and tapering not yet implemented for REAL8 data */
  if ( strcmp( thisInj->taper, "TAPER_NONE" ) || thisInj->bandpass )
  {
    XLALPrintError( "Tapering/band-passing of REAL8 injection currently unsupported\n" );
    ABORT( status, NRWAVEINJECT_EVAL, NRWAVEINJECT_MSGEVAL );
  }

  /* sampleRate = 1.0/strain->deltaT;   */
  /* use the sample rate required for the output time series */
  sampleRate = 1.0/injData->deltaT;
  len = strain->data->vectorLength;

  hplus = XLALCreateREAL8TimeSeries ( strain->name, &strain->epoch, strain->f0,
				      strain->deltaT, &strain->sampleUnits, len);
  hcross = XLALCreateREAL8TimeSeries ( strain->name, &strain->epoch, strain->f0,
				       strain->deltaT, &strain->sampleUnits, len);

  for ( k = 0; k < len; k++) {
    hplus->data->data[k] = strain->data->data[k];
    hplus->data->data[k] = strain->data->data[k + len];
  }

  htData = XLALSimDetectorStrainREAL8TimeSeries( hplus, hcross, thisInj->longitude,
						 thisInj->latitude, thisInj->polarization,
						 &det);

  XLALFindNRCoalescenceTimeREAL8( &offset, htData);
  XLALGPSAdd( &(htData->epoch), -offset);

  XLALSimAddInjectionREAL8TimeSeries( injData, htData, NULL);

  /* set channel name */
  LALSnprintf( injData->name, LIGOMETA_CHANNEL_MAX * sizeof( CHAR ),
    "%s:STRAIN", ifo );

  XLALDestroyREAL8TimeSeries ( htData);
  XLALDestroyREAL8TimeSeries ( hplus);
  XLALDestroyREAL8TimeSeries ( hcross);

  DETATCHSTATUSPTR(status);
  RETURN(status);

}


/** construct the channel name corresponding to a particular mode
    and polarization in frame file containing nr data */
CHAR* XLALGetNinjaChannelName(CHAR *polarisation, UINT4 l, INT4 m)
{
  /* variables */
  CHAR sign;
  CHAR *channel=NULL;

  if ( !((strncmp(polarisation, "plus", 4) == 0) || (strncmp(polarisation, "cross", 5) == 0))) {
    XLAL_ERROR_NULL( "XLALGetNinjaChannelName",   XLAL_EINVAL );
  }

  /* allocate memory for channel */
  channel = (CHAR *)LALCalloc(1, LIGOMETA_CHANNEL_MAX * sizeof(CHAR));

  /* get sign of m */
  if (m < 0)
  {
    /* negative */
    strncpy(&sign, "n", 1);
  }
  else
  {
    /* positive */
    strncpy(&sign, "p", 1);
  }

  /* set channel name */
  LALSnprintf(channel, LIGOMETA_CHANNEL_MAX, "h%s_l%d_m%c%d", polarisation, l, sign, abs(m));

  /* return channel name */
  return channel;
}


/** Function for parsing numrel group name and converting it into a enum element.
    This needs to be robust enough to be able to handle the information as submitted
    by the groups. Is there a cleaner way to do this?
    add or modify the group names as required
*/
NumRelGroup XLALParseNumRelGroupName( CHAR *name)
{

  NumRelGroup ret=NINJA_GROUP_LAST;

  if ( !name ) {
    return ret;
  }

  if ( strstr( name, "Caltech") || strstr( name, "CIT") )
    ret = NINJA_GROUP_CIT;
  else if ( strstr( name, "AEI") || strstr(name, "Potsdam")  )
    ret = NINJA_GROUP_AEI;
  else if ( strstr( name, "Jena") )
    ret = NINJA_GROUP_JENA;
  else if ( strstr( name, "Pretorius") || strstr( name, "Princeton"))
    ret = NINJA_GROUP_PRINCETON;
  else if ( strstr( name, "Cornell") )
    ret = NINJA_GROUP_CORNELL;
  else if ( strstr( name, "PSU") || strstr(name, "Penn State") )
    ret = NINJA_GROUP_PSU;
  else if ( strstr( name, "LSU") )
    ret = NINJA_GROUP_LSU;
  else if ( strstr( name, "RIT") || strstr( name, "Rochester") )
    ret = NINJA_GROUP_RIT;
  else if ( strstr( name, "FAU") || strstr( name, "Florida Atlantic") )
    ret = NINJA_GROUP_FAU;
  else if ( strstr( name, "UTB") )
    ret = NINJA_GROUP_UTB;
  else if ( strstr( name, "UIUC") || strstr(name, "Urbana Champaign") )
    ret = NINJA_GROUP_UIUC;

  return ret;
}
