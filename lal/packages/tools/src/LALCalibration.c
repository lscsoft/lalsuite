#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/LALCalibration.h>


void XLALDestroyCalData( LALCalData *caldata )
{
  XLALDestroyCOMPLEX8FrequencySeries( caldata->openLoopGainReference );
  XLALDestroyCOMPLEX8FrequencySeries( caldata->cavityGainReference );
  XLALDestroyCOMPLEX8FrequencySeries( caldata->responseReference );
  XLALDestroyREAL4TimeSeries( caldata->openLoopFactors );
  XLALDestroyREAL4TimeSeries( caldata->cavityFactors );
  XLALFree( caldata );
  return;
}


int XLALUpdateResponse( 
    COMPLEX8FrequencySeries *response,  /**< response function to return */
    REAL8 duration,                     /**< duration for averaging factors */
    LALCalData *caldata                 /**< calibration reference data */
    )
{
  static const char *func = "XLALUpdateResponse";
  const REAL4 tiny = 1e-6;
  INT4 offset;
  INT4 npts;
  INT4 i;
  UINT4 k;

  REAL8 alpha;
  REAL8 gamma;

  if ( ! response || ! caldata )
    XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! caldata->cavityFactors || ! caldata->openLoopFactors 
      || ! caldata->responseReference || ! caldata->cavityGainReference
      || ! caldata->openLoopGainReference )
    XLAL_ERROR( func, XLAL_EINVAL );

  /* make sure alpha and gamma factors have consistent interval */

  /* alpha and gamma factors have to start at same time */
  if ( XLALGPSCmp( &caldata->cavityFactors->epoch, &caldata->openLoopFactors->epoch ) )
    XLAL_ERROR( func, XLAL_ETIME );
  /* alpha and gamma factors have to have same time step */
  if ( fabs( caldata->cavityFactors->deltaT - caldata->openLoopFactors->deltaT ) > tiny )
    XLAL_ERROR( func, XLAL_ETIME );
  /* alpha and gamma factors have to have same length */
  if ( caldata->cavityFactors->data->length != caldata->openLoopFactors->data->length )
    XLAL_ERROR( func, XLAL_EBADLEN );

  /* reference functions have to have same epoch */
  if ( XLALGPSCmp( &caldata->cavityGainReference->epoch, &caldata->responseReference->epoch ) || XLALGPSCmp( &caldata->openLoopGainReference->epoch, &caldata->responseReference->epoch ) )
    XLAL_ERROR( func, XLAL_ETIME );
  /* reference functions have to have same frequencies */
  if ( fabs( caldata->cavityGainReference->deltaF - caldata->responseReference->deltaF ) > tiny || fabs( caldata->openLoopGainReference->deltaF - caldata->responseReference->deltaF ) > tiny )
    XLAL_ERROR( func, XLAL_EFREQ );
  /* reference functions have to have same length */
  if ( caldata->cavityGainReference->data->length != caldata->responseReference->data->length || caldata->openLoopGainReference->data->length != caldata->responseReference->data->length )
    XLAL_ERROR( func, XLAL_EBADLEN );


  /* requested epoch must be after reference function validity */
  if ( XLALGPSCmp(&response->epoch,&caldata->responseReference->epoch) < 0 )
    XLAL_ERROR( func, XLAL_ETIME );


  /* figure out how many points of alpha and gamma to average */
  if ( duration == 0 )
    npts = 1;
  else
    npts = floor( duration / caldata->cavityFactors->deltaT + 0.5 );
  if ( npts < 1 )
    XLAL_ERROR( func, XLAL_ETIME );


  /* compute offset for averaging alpha and gamma */
  /* first point at or before the requested time */
  offset = floor( XLALGPSDiff( &response->epoch, &caldata->cavityFactors->epoch ) / caldata->cavityFactors->deltaT );
  if ( offset < 0 )
    XLAL_ERROR( func, XLAL_ETIME );

  /* make sure there is enough factors data */
  if ( offset + npts > (INT4)caldata->cavityFactors->data->length )
    XLAL_ERROR( func, XLAL_ESIZE );

  /* compute average alpha and gamma */
  alpha = gamma = 0;
  for ( i = 0; i < npts; ++i )
  {
    REAL4 alp = caldata->cavityFactors->data->data[i+offset];
    REAL4 gam = caldata->openLoopFactors->data->data[i+offset];
    if ( alp < tiny || gam < tiny )
      XLALPrintWarning( "XLAL Warning - %s: Zero or negative factor found\n\talpha=%g gamma=%g\n", func, alp, gam );
    alpha += alp;
    gamma += gam;
  }
  alpha /= npts;
  gamma /= npts;
  XLALPrintInfo( "XLAL Info - %s: Using calibration factors alpha=%g gamma=%g\n", func, alpha, gamma );

  /* interpolate to requested frequencies */
  for ( k = 0; k < response->data->length; ++k )
  {
    COMPLEX8 G;
    COMPLEX8 C;
    COMPLEX8 R1;
    COMPLEX8 R2;
    REAL8 x;
    UINT4 j;

    /* compute response at two points for interpolation */
    x = k * response->deltaF / caldata->openLoopGainReference->deltaF;
    j = floor( x );
    if ( j > caldata->openLoopGainReference->data->length - 2 )
      j = caldata->openLoopGainReference->data->length - 2;
    x -= j;


    /*
     *
     * First Point
     *
     */

    C.re = caldata->cavityGainReference->data->data[j].re;
    C.im = caldata->cavityGainReference->data->data[j].im;
    G.re = caldata->openLoopGainReference->data->data[j].re;
    G.im = caldata->openLoopGainReference->data->data[j].im;

    /* update using averaged factors */
    C.re *= alpha;
    C.im *= alpha;
    G.re *= gamma;
    G.im *= gamma;

    /* now compute the response function */
    /* first add one to G and then divide by C */
    if ( C.re == 0 && C.im == 0 )
      R1.re = R1.im = 0;
    else
    {
      G.re += 1.0;
      if ( fabs( C.re ) >= fabs( C.im ) )
      {
        REAL4 rat = C.im / C.re;
        REAL4 den = C.re + rat * C.im;
        R1.re = ( G.re + rat * G.im ) / den;
        R1.im = ( G.im - rat * G.re ) / den;
      }
      else
      {
        REAL4 rat = C.re / C.im;
        REAL4 den = C.im + rat * C.re;
        R1.re = ( G.re * rat + G.im ) / den;
        R1.im = ( G.im * rat - G.re ) / den;
      }
    }


    /*
     *
     * Second Point
     *
     */

    C.re = caldata->cavityGainReference->data->data[j+1].re;
    C.im = caldata->cavityGainReference->data->data[j+1].im;
    G.re = caldata->openLoopGainReference->data->data[j+1].re;
    G.im = caldata->openLoopGainReference->data->data[j+1].im;

    /* update using averaged factors */
    C.re *= alpha;
    C.im *= alpha;
    G.re *= gamma;
    G.im *= gamma;

    /* now compute the response function */
    /* first add one to G and then divide by C */
    if ( C.re == 0 && C.im == 0 )
      R2.re = R2.im = 0;
    else
    {
      G.re += 1.0;
      if ( fabs( C.re ) >= fabs( C.im ) )
      {
        REAL4 rat = C.im / C.re;
        REAL4 den = C.re + rat * C.im;
        R2.re = ( G.re + rat * G.im ) / den;
        R2.im = ( G.im - rat * G.re ) / den;
      }
      else
      {
        REAL4 rat = C.re / C.im;
        REAL4 den = C.im + rat * C.re;
        R2.re = ( G.re * rat + G.im ) / den;
        R2.im = ( G.im * rat - G.re ) / den;
      }
    }


    /*
     *
     * Now Interpolate
     *
     */

    {
      REAL8 rad1 = sqrt( R1.re * R1.re + R1.im * R1.im );
      REAL8 rad2 = sqrt( R2.re * R2.re + R2.im * R2.im );
      REAL8 arg1 = atan2( R1.im, R1.re );
      REAL8 arg2 = atan2( R2.im, R2.re );
      REAL8 rad;
      REAL8 arg;
      rad = rad1 + x * (rad2 - rad1);
      arg = arg1 + x * (arg2 - arg1);
      response->data->data[k].re = rad * cos( arg );
      response->data->data[k].im = rad * sin( arg );
    }

  }

  return 0;
}

