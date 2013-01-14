/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton
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

#include <complex.h>
#include <math.h>
#include <string.h>

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/VectorOps.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/LALCalibration.h>

void XLALDestroyCalData( LALCalData *caldata )
{
  XLALDestroyCOMPLEX8FrequencySeries( caldata->digitalFilterReference );
  XLALDestroyCOMPLEX8FrequencySeries( caldata->actuationReference );
  XLALDestroyCOMPLEX8FrequencySeries( caldata->openLoopGainReference );
  XLALDestroyCOMPLEX8FrequencySeries( caldata->cavityGainReference );
  XLALDestroyCOMPLEX8FrequencySeries( caldata->responseReference );
  XLALDestroyREAL4TimeSeries( caldata->openLoopFactors );
  XLALDestroyREAL4TimeSeries( caldata->cavityFactors );
  XLALFree( caldata );
  return;
}

int XLALAverageCalibrationFactors(
    REAL8 *cal_alpha,    /**< returned value of alpha factor */
    REAL8 *cal_gamma,    /**< returned value of gamma factor */
    LIGOTimeGPS *epoch,  /**< epoch to begin averaging */
    REAL8 duration,      /**< duration of averaging (0 = use just one point) */
    LALCalData *caldata  /**< calibration reference data */
    )
{
  const REAL4 tiny = 1e-6;
  INT4 offset;
  INT4 npts;
  INT4 i;

  /* check for invalid input */
  if ( ! cal_alpha || ! cal_gamma || ! epoch || ! caldata )
    XLAL_ERROR( XLAL_EFAULT );
  if (! caldata->cavityFactors || ! caldata->openLoopFactors || duration < 0.0)
    XLAL_ERROR( XLAL_EINVAL );

  /* make sure cal_alpha and cal_gamma factors have consistent interval */

  /* cal_alpha and cal_gamma factors have to start at same time */
  if ( XLALGPSCmp( &caldata->cavityFactors->epoch, &caldata->openLoopFactors->epoch ) )
    XLAL_ERROR( XLAL_ETIME );
  /* cal_alpha and cal_gamma factors have to have same time step */
  if ( fabs( caldata->cavityFactors->deltaT - caldata->openLoopFactors->deltaT ) > tiny )
    XLAL_ERROR( XLAL_ETIME );
  /* cal_alpha and cal_gamma factors have to have same length */
  if ( caldata->cavityFactors->data->length != caldata->openLoopFactors->data->length )
    XLAL_ERROR( XLAL_EBADLEN );

  /* compute offset for averaging cal_alpha and cal_gamma */
  /* first point at or before the requested time */
  offset = floor( XLALGPSDiff( epoch, &caldata->cavityFactors->epoch ) / caldata->cavityFactors->deltaT );
  if ( offset < 0 )
    XLAL_ERROR( XLAL_ETIME );

  /* figure out how many points of cal_alpha and cal_gamma to average */
  npts = floor( duration / caldata->cavityFactors->deltaT + 0.5 );
  if ( npts < 1 )
    npts = 1;

  /* make sure there is enough factors data */
  if ( offset + npts > (INT4)caldata->cavityFactors->data->length )
    XLAL_ERROR( XLAL_ESIZE );

  /* compute average cal_alpha and cal_gamma */
  *cal_alpha = *cal_gamma = 0;
  for ( i = 0; i < npts; ++i )
  {
    REAL4 alp = caldata->cavityFactors->data->data[i+offset];
    REAL4 gam = caldata->openLoopFactors->data->data[i+offset];
    if ( alp < tiny || gam < tiny )
      XLALPrintWarning( "XLAL Warning - %s: Zero or negative factor found\n\talpha=%g gamma=%g\n", __func__, alp, gam );
    *cal_alpha += alp;
    *cal_gamma += gam;
  }
  *cal_alpha /= npts;
  *cal_gamma /= npts;

  return 0;
}

COMPLEX8FrequencySeries * XLALUpdateReferenceResponse(
    LALCalData *caldata, /**< calibration reference data */
    REAL8 cal_alpha,     /**< value of the cal_alpha factor */
    REAL8 cal_gamma      /**< value of the cal_gamma factor */
    )
{
  enum { DARM_ERR, AS_Q } type;
  COMPLEX8FrequencySeries *response;
  const REAL4 tiny = 1e-6;
  UINT4 k;

  if ( ! caldata )
    XLAL_ERROR_NULL( XLAL_EFAULT );
  if (!caldata->responseReference || !caldata->responseReference->data ||
      !caldata->cavityGainReference || !caldata->cavityGainReference->data ||
      !caldata->openLoopGainReference || !caldata->openLoopGainReference->data )
    XLAL_ERROR_NULL( XLAL_EINVAL );

  /* reference functions have to have same epoch */
  if ( XLALGPSCmp( &caldata->cavityGainReference->epoch, &caldata->responseReference->epoch ) || XLALGPSCmp( &caldata->openLoopGainReference->epoch, &caldata->responseReference->epoch ) )
    XLAL_ERROR_NULL( XLAL_ETIME );
  /* reference functions have to have same frequencies */
  if ( fabs( caldata->cavityGainReference->deltaF - caldata->responseReference->deltaF ) > tiny || fabs( caldata->openLoopGainReference->deltaF - caldata->responseReference->deltaF ) > tiny )
    XLAL_ERROR_NULL( XLAL_EFREQ );
  /* reference functions have to have same length */
  if ( caldata->cavityGainReference->data->length != caldata->responseReference->data->length || caldata->openLoopGainReference->data->length != caldata->responseReference->data->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );

  /* determine if we are obtaining a response function for AS_Q or DARM_ERR */
  if ( strstr( caldata->responseReference->name, "DARM_ERR" ) )
    type = DARM_ERR;
  else if ( strstr( caldata->responseReference->name, "AS_Q" ) )
    type = AS_Q;
  else
    XLAL_ERROR_NULL( XLAL_ENAME );

  /* copy reference response */
  response = XLALCutCOMPLEX8FrequencySeries( caldata->responseReference, 0, caldata->responseReference->data->length );
  /* NOTE: epoch of returned response is the same as the reference response
   * since there is no input here as to which epoch the response function
   * should have had. */

  /* update reference response with specified factors */
  for ( k = 0; k < response->data->length; ++k )
  {
    COMPLEX8 C = caldata->cavityGainReference->data->data[k];
    COMPLEX8 G = caldata->openLoopGainReference->data->data[k];
    G *= cal_gamma;
    C *= ((type == DARM_ERR) ? cal_gamma : cal_alpha);
    response->data->data[k] = ((C == 0.0) ? 0.0 : (G + 1.0) / C);
  }

  return response;
}

/* functional macros for procedures that get reused */

#define GET_POLAR_RESPONSE( rad, phi, epoch, duration, caldata ) \
  do { \
    COMPLEX8FrequencySeries *rawresponse; \
    REAL8 cal_alpha, cal_gamma; \
    UINT4 k; \
    rad = phi = NULL; \
    if ( 0 > XLALAverageCalibrationFactors( &cal_alpha, &cal_gamma, epoch, duration, caldata ) ) break; \
    XLALPrintInfo( "XLAL Info - %s: Using calibration factors alpha=%g gamma=%g\n", __func__, cal_alpha, cal_gamma ); \
    if ( ! ( rawresponse = XLALUpdateReferenceResponse( caldata, cal_alpha, cal_gamma ) ) ) break; \
    rad = XLALCreateREAL8FrequencySeries( "response magnitude", epoch, rawresponse->f0, rawresponse->deltaF, &rawresponse->sampleUnits, rawresponse->data->length ); \
    phi = XLALCreateREAL8FrequencySeries( "response phase", epoch, rawresponse->f0, rawresponse->deltaF, &lalDimensionlessUnit, rawresponse->data->length ); \
    if ( ! rad || ! phi ) break; \
    for ( k = 0; k < rawresponse->data->length; ++k ) { \
      rad->data->data[k] = cabsf( rawresponse->data->data[k] ); \
      phi->data->data[k] = cargf( rawresponse->data->data[k] ); \
    } \
    XLALDestroyCOMPLEX8FrequencySeries( rawresponse ); \
  } while (0)


#define POLAR_INTERPOLATE_RESPONSE( response, rad, phi ) \
  do { \
    UINT4 j, k; \
    XLALREAL8VectorUnwrapAngle( phi->data, phi->data ); \
    for ( k = 0; k < response->data->length; ++k ) { \
      REAL8 r, p, x; \
      j = floor( x = k * response->deltaF / rad->deltaF ); \
      if ( j > rad->data->length - 2 ) j = rad->data->length - 2; \
      x -= j; \
      r = (1.0 - x)*rad->data->data[j] + x*rad->data->data[j+1]; \
      p = (1.0 - x)*phi->data->data[j] + x*phi->data->data[j+1]; \
      response->data->data[k] = r*cos(p) + I*r*sin(p); \
    } \
  } while (0)


COMPLEX8FrequencySeries * XLALCreateCOMPLEX8Response(
    LIGOTimeGPS *epoch,    /**< epoch of response function */
    REAL8        duration, /**< duration of averaging (0 = use one point) */
    REAL8        deltaF,   /**< frequency resolution of requested response */
    UINT4        length,   /**< number of frequency bins in response */
    LALCalData  *caldata   /**< calibration reference data */
    )
{
  COMPLEX8FrequencySeries *response;
  REAL8FrequencySeries *rad;
  REAL8FrequencySeries *phi;

  GET_POLAR_RESPONSE( rad, phi, epoch, duration, caldata );
  if ( ! rad || ! phi )
    XLAL_ERROR_NULL( XLAL_EFUNC );

  response = XLALCreateCOMPLEX8FrequencySeries( "response", epoch, 0.0, deltaF, &rad->sampleUnits, length );

  POLAR_INTERPOLATE_RESPONSE( response, rad, phi );

  XLALDestroyREAL8FrequencySeries( phi );
  XLALDestroyREAL8FrequencySeries( rad );

  return response;
}


COMPLEX16FrequencySeries * XLALCreateCOMPLEX16Response(
    LIGOTimeGPS *epoch,    /**< epoch of response function */
    REAL8        duration, /**< duration of averaging (0 = use one point) */
    REAL8        deltaF,   /**< frequency resolution of requested response */
    UINT4        length,   /**< number of frequency bins in response */
    LALCalData  *caldata   /**< calibration reference data */
    )
{
  COMPLEX16FrequencySeries *response;
  REAL8FrequencySeries *rad;
  REAL8FrequencySeries *phi;

  GET_POLAR_RESPONSE( rad, phi, epoch, duration, caldata );
  if ( ! rad || ! phi )
    XLAL_ERROR_NULL( XLAL_EFUNC );

  response = XLALCreateCOMPLEX16FrequencySeries( "response", epoch, 0.0, deltaF, &rad->sampleUnits, length );

  POLAR_INTERPOLATE_RESPONSE( response, rad, phi );

  XLALDestroyREAL8FrequencySeries( phi );
  XLALDestroyREAL8FrequencySeries( rad );

  return response;
}


int XLALUpdateResponse(
    COMPLEX8FrequencySeries *response,  /**< response function to return */
    REAL8 duration,                     /**< duration for averaging factors */
    LALCalData *caldata                 /**< calibration reference data */
    )
{
  REAL8FrequencySeries *rad;
  REAL8FrequencySeries *phi;

  GET_POLAR_RESPONSE( rad, phi, &response->epoch, duration, caldata );
  if ( ! rad || ! phi )
    XLAL_ERROR( XLAL_EFUNC );

  response->sampleUnits = rad->sampleUnits;

  POLAR_INTERPOLATE_RESPONSE( response, rad, phi );

  XLALDestroyREAL8FrequencySeries( phi );
  XLALDestroyREAL8FrequencySeries( rad );

  return 0;
}
