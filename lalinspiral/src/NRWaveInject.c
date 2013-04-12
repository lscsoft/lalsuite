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

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/NRWaveIO.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataInspiralUtils.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/LALConstants.h>
#include <lal/NRWaveInject.h>
#include <lal/Random.h>
#include <lal/Inject.h>
#include <lal/LALSimulation.h>
#include <lal/LALDetectors.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/SkyCoordinates.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Window.h>
#include <lal/SphericalHarmonics.h>

#include <gsl/gsl_heapsort.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

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
      (tmp1 * creal(MultSphHarm)) +
      (tmp2 * MultSphHarm.im);

    strain->data->data[vecLength + k] =
      (tmp2 * creal(MultSphHarm)) -
      (tmp1 * MultSphHarm.im);
  }

  return 0;
  /*   return( strain ); */
}

/** Takes a (sky averaged) numerical relativity waveform and returns the
 * waveform appropriate for given coalescence phase and inclination angles */
void
XLALOrientNRWaveTimeSeriesREAL8(
    REAL8TimeSeries        *plus,	   /**< NEEDS DOCUMENTATION */
    REAL8TimeSeries        *cross,	   /**< NEEDS DOCUMENTATION */
    UINT4                  modeL,          /**< L                       */
    INT4                   modeM,          /**< M                       */
    REAL4                  inclination,    /**< binary inclination      */
    REAL4                  coa_phase       /**< binary coalescence phase*/)
{
  COMPLEX16  MultSphHarm;
  REAL4      tmp1, tmp2;
  UINT4      vecLength, k;

  vecLength = plus->data->length;

  /* Calculating the (2,2) Spherical Harmonic */
  /* need some error checking */
  XLALSphHarm( &MultSphHarm, modeL, modeM, inclination, coa_phase );

  /* Filling the data vector with the data multiplied by the Harmonic */
  for ( k = 0; k < vecLength; k++)
  {
    tmp1 = plus->data->data[k];
    tmp2 = cross->data->data[k];

    plus->data->data[k] =
      (tmp1 * creal(MultSphHarm)) +
      (tmp2 * MultSphHarm.im);

    cross->data->data[k] =
      (tmp2 * creal(MultSphHarm)) -
      (tmp1 * MultSphHarm.im);
  }

  return;
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
      (tmp1 * creal(MultSphHarm)) +
      (tmp2 * MultSphHarm.im);

    strain->data->data[vecLength + k] =
      (tmp2 * creal(MultSphHarm)) -
      (tmp1 * MultSphHarm.im);
  }
  return( strain );
}


  REAL4TimeSeries *
XLALCalculateNRStrain( REAL4TimeVectorSeries *strain, /**< h+, hx time series data*/
    SimInspiralTable      *inj,    /**< injection details      */
    const CHAR            *ifo,    /**< interferometer */
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
  ifoNumber = (InterferometerNumber) XLALIFONumber( ifo );
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
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }
  vecLength = strain->data->vectorLength;
  htData->data = XLALCreateREAL4Vector( vecLength );
  if ( ! htData->data )
  {
    LALFree( htData );
    XLAL_ERROR_NULL( XLAL_ENOMEM );
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
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  ret->data = XLALCreateREAL4Vector( numPoints );
  if (! ret->data)
  {
    XLAL_ERROR_NULL( XLAL_ENOMEM );
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
    /* Copy the end point if we are exactly at it */
    else if(lo==in->data->length-1){
      ret->data->data[k]=in->data->data[lo];
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


/** Function for interpolating time series to a given sampling rate.
  Input vector is destroyed and a new vector is allocated.
  */
  REAL8TimeSeries *
XLALInterpolateNRWaveREAL8( REAL8TimeSeries *in,           /**< input strain time series */
    INT4            sampleRate     /**< sample rate of time series */)
{

  REAL8TimeSeries *ret=NULL;
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
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  ret->data = XLALCreateREAL8Vector( numPoints );
  if (! ret->data)
  {
    XLAL_ERROR_NULL( XLAL_ENOMEM );
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
    }    /* Copy the end point if we are exactly at it */
    else if(lo==in->data->length-1){
      ret->data->data[k]=in->data->data[lo];
    }
    else {
      ret->data->data[k] = 0.0;
    }
  }

  /* destroy input vector */
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
XLALFindNRCoalescenceTime(REAL8 *tc,  /**< FIXME: !TO BE DOCUMENTED! */
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


  INT4
XLALFindNRCoalescencePlusCrossREAL8(REAL8 *tc,  /**< FIXME: !TO BE DOCUMENTED! */
    const REAL8TimeSeries *plus,   /**< input strain plus time series */
    const REAL8TimeSeries *cross   /**< input strain cross time series */)
{

  size_t *ind=NULL;
  size_t len;
  REAL8 *sumSquare=NULL;
  UINT4 k;

  len = plus->data->length;
  ind = LALCalloc(1,len*sizeof(*ind));

  sumSquare = LALCalloc(1, len*sizeof(*sumSquare));

  for (k=0; k < len; k++) {
    sumSquare[k] = plus->data->data[k]*plus->data->data[k] +
      cross->data->data[k]*cross->data->data[k];
  }

  gsl_heapsort_index( ind, sumSquare, len, sizeof(REAL8), compare_abs_double);

  *tc = ind[len-1] * plus->deltaT;

  LALFree(ind);
  LALFree(sumSquare);

  return 0;
}


/** Function for calculating the coalescence time (defined to be the
  peak) of a NR wave
  This uses the peak of h(t)
  */
  INT4
XLALFindNRCoalescenceTimeFromhoft(REAL8 *tc,   /**< FIXME: !TO BE DOCUMENTED! */
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
XLALFindNRCoalescenceTimeREAL8(REAL8 *tc,  /**< FIXME: !TO BE DOCUMENTED! */
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
    XLAL_ERROR ( XLAL_EINVAL);
  }

  if ( !nrCatalog ) {
    LALPrintError ("\n NR Catalog pointer is NULL !\n\n");
    XLAL_ERROR ( XLAL_EINVAL);
  }

  if ( !inj ) {
    LALPrintError ("\n SimInspiralTable pointer is NULL !\n\n");
    XLAL_ERROR ( XLAL_EINVAL);
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
      XLAL_ERROR ( XLAL_EINVAL);
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
    XLAL_ERROR ( XLAL_EINVAL);
  }

  /* copy best match to output */
  memcpy(out, nrCatalog->data + best, sizeof(NRWaveMetaData));

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
  LALSimInspiralApplyTaper taper = LAL_SIM_INSPIRAL_TAPER_NONE;

  INITSTATUS(status);
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
      taper = LAL_SIM_INSPIRAL_TAPER_START;
    }
    else if (  ! strcmp( "TAPER_END", thisInj->taper ) )
    {
      taper = LAL_SIM_INSPIRAL_TAPER_END;
    }
    else if (  ! strcmp( "TAPER_STARTEND", thisInj->taper ) )
    {
      taper = LAL_SIM_INSPIRAL_TAPER_STARTEND;
    }
    else
    {
      XLALPrintError( "Unsupported tapering type specified: %s\n", thisInj->taper );
      XLALDestroyREAL4Vector ( htData->data);
      LALFree(htData);
      ABORT( status, NRWAVEINJECT_EVAL, NRWAVEINJECT_MSGEVAL );
    }
    if ( XLALSimInspiralREAL4WaveTaper( htData->data, taper ) == XLAL_FAILURE )
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

  /* Cast to REAL8 for injection */
  REAL8TimeSeries *injData8=XLALCreateREAL8TimeSeries("temporary", &(injData->epoch), injData->f0, injData->deltaT, &(injData->sampleUnits), injData->data->length);
  if(!injData8) XLAL_ERROR_VOID(XLAL_ENOMEM,"Unable to allocate injection buffer\n");
  REAL8TimeSeries *htData8=XLALCreateREAL8TimeSeries("signal, temp", &(htData->epoch), htData->f0, htData->deltaT, &(htData->sampleUnits), htData->data->length);
  if(!htData8) XLAL_ERROR_VOID(XLAL_ENOMEM,"Unable to allocate signal buffer\n");

  for(UINT4 i=0;i<htData->data->length;i++) htData8->data->data[i]=(REAL8)htData->data->data[i];
  for(UINT4 i=0;i<injData->data->length;i++) injData8->data->data[i]=(REAL8)injData->data->data[i];

  /* inject the htData into injection time stream */
  int retcode = XLALSimAddInjectionREAL8TimeSeries(injData8, htData8, NULL);
  if(retcode!=XLAL_SUCCESS)
  {
    XLALDestroyREAL8TimeSeries(htData8);
    XLALDestroyREAL8TimeSeries(injData8);
    XLALDestroyREAL4Vector(htData->data);
    LALFree(htData);
    ABORTXLAL(status);
  }

  for(UINT4 i=0;i<injData8->data->length;i++) injData->data->data[i]=(REAL4)injData8->data->data[i];

  XLALDestroyREAL8TimeSeries(injData8);
  XLALDestroyREAL8TimeSeries(htData8);


  /* set channel name */
  snprintf( injData->name, sizeof(injData->name),
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
    REAL8                     UNUSED dynRange)
{

  REAL8TimeSeries *htData = NULL;
  REAL8TimeSeries *hplus = NULL;
  REAL8TimeSeries *hcross = NULL;
  UINT4  k, len;
  REAL8 offset;
  InterferometerNumber  ifoNumber = LAL_UNKNOWN_IFO;
  LALDetector det;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* get the detector information */
  memset( &det, 0, sizeof(LALDetector) );
  ifoNumber = (InterferometerNumber) XLALIFONumber( ifo );
  XLALReturnDetector( &det, ifoNumber );

  /* Band-passing and tapering not yet implemented for REAL8 data */
  if ( strcmp( thisInj->taper, "TAPER_NONE" ) || thisInj->bandpass )
  {
    XLALPrintError( "Tapering/band-passing of REAL8 injection currently unsupported\n" );
    ABORT( status, NRWAVEINJECT_EVAL, NRWAVEINJECT_MSGEVAL );
  }

  /* use the sample rate required for the output time series */
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

  int retcode=XLALSimAddInjectionREAL8TimeSeries( injData, htData, NULL);
  XLAL_CHECK_VOID(retcode==XLAL_SUCCESS, XLAL_EFUNC, "Unable to add injection to data\n");

  /* set channel name */
  snprintf( injData->name, sizeof(injData->name),
      "%s:STRAIN", ifo );

  XLALDestroyREAL8TimeSeries ( htData);
  XLALDestroyREAL8TimeSeries ( hplus);
  XLALDestroyREAL8TimeSeries ( hcross);

  DETATCHSTATUSPTR(status);
  RETURN(status);

}


/** construct the channel name corresponding to a particular mode
  and polarization in frame file containing nr data */
CHAR* XLALGetNinjaChannelName(const CHAR *polarisation, UINT4 l, INT4 m)
{
  /* variables */
  CHAR sign;
  CHAR *channel=NULL;

  if ( !((strncmp(polarisation, "plus", 4) == 0) || (strncmp(polarisation, "cross", 5) == 0))) {
    XLAL_ERROR_NULL( XLAL_EINVAL );
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
  snprintf(channel, LIGOMETA_CHANNEL_MAX, "h%s_l%d_m%c%d", polarisation, l, sign, abs(m));

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
