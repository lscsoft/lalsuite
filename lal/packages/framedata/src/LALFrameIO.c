#include <math.h>

#include <lal/LALDetectors.h>
#include <lal/LALStdio.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>

#include <lal/LALString.h>
#include <lal/LALFrameIO.h>
#include <lal/LALCalibration.h>


FrDetector * XLALFrDetectorNew( int detector )
{
  static const char *func = "XLALFrDetectorNew";
  const LALDetector *lalDetector;
  FrDetector *frDetector;
  char *detectorName;

  if ( detector < 0 || detector >= LAL_NUM_DETECTORS )
    XLAL_ERROR_NULL( func, XLAL_EINVAL );
  lalDetector = lalCachedDetectors + detector;

  detectorName = XLALStringDuplicate( lalDetector->frDetector.name );
  frDetector = FrDetectorNew( detectorName );
  LALFree( detectorName );
  if ( ! frDetector )
    XLAL_ERROR_NULL( func, XLAL_EERR ); /* "internal" error */

  memcpy( frDetector->prefix, lalDetector->frDetector.prefix, 2 );
  frDetector->longitude    = lalDetector->frDetector.vertexLongitudeRadians;
  frDetector->latitude     = lalDetector->frDetector.vertexLatitudeRadians;
  frDetector->elevation    = lalDetector->frDetector.vertexElevation;
  frDetector->armXazimuth  = lalDetector->frDetector.xArmAzimuthRadians;
  frDetector->armYazimuth  = lalDetector->frDetector.yArmAzimuthRadians;
  frDetector->armXaltitude = lalDetector->frDetector.xArmAltitudeRadians;
  frDetector->armYaltitude = lalDetector->frDetector.yArmAltitudeRadians;
  frDetector->armXmidpoint = lalDetector->frDetector.xArmMidpoint;
  frDetector->armYmidpoint = lalDetector->frDetector.yArmMidpoint;
  frDetector->localTime    = 0;

  return frDetector;
}


FrameH * XLALFrameNew( LIGOTimeGPS *epoch, double duration, 
    const char *project, int run, int frnum, int detectorFlags )
{
  static const char *func = "XLALFrameNew";
  static char rcsid[] = "$Id$";
  int detector;
  char *proj;
  FrameH *frame;

  proj = XLALStringDuplicate( project );
  if ( ! proj )
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );
  frame = FrameHNew( proj );
  LALFree( proj );
  if ( ! frame )
    XLAL_ERROR_NULL( func, XLAL_EERR ); /* "internal" error */

  frame->run    = run;
  frame->frame  = frnum;
  frame->GTimeS = epoch->gpsSeconds;
  frame->GTimeN = epoch->gpsNanoSeconds;
  frame->ULeapS = XLALLeapSeconds( epoch->gpsSeconds );
  frame->dt     = duration;
  /* TODO: what about dataQuality ? */

  for ( detector = 0; detector < LAL_NUM_DETECTORS; ++detector )
  {
    int detector_flag = 1 << 2 * detector;
    if ( ( detector_flag & detectorFlags ) ) /* yes, one ampersand! */
    {
      /* add this detector */
      FrDetector *d;
      d = XLALFrDetectorNew( detector );
      if ( ! d )
        XLAL_ERROR_NULL( func, XLAL_EFUNC );
      d->next = frame->detectProc;
      frame->detectProc = d;
    }
  }

  FrHistoryAdd( frame, rcsid );
  return frame;
}


FrVect * XLALFrVectREAL4TimeSeries( REAL4TimeSeries *series )
{
  static const char *func = "XLALFrVectREAL4TimeSeries";
  char seconds[LALUnitTextSize] = "s";
  char units[LALUnitTextSize];
  FrVect *vect;

  if ( NULL == XLALUnitAsString( units, sizeof( units ), &series->sampleUnits ) )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  vect = FrVectNew1D( series->name, FR_VECT_4R, series->data->length, series->deltaT, seconds, units );
  if ( ! vect )
    XLAL_ERROR_NULL( func, XLAL_EERR ); /* "internal" error */
  vect->startX[0] = 0.0;

  memcpy( vect->data, series->data->data, series->data->length * sizeof( *series->data->data ) );

  return vect;
}


FrVect * XLALFrVectREAL8TimeSeries( REAL8TimeSeries *series )
{
  static const char *func = "XLALFrVectREAL8TimeSeries";
  char seconds[LALUnitTextSize] = "s";
  char units[LALUnitTextSize];
  FrVect *vect;

  if ( NULL == XLALUnitAsString( units, sizeof( units ), &series->sampleUnits ) )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  vect = FrVectNew1D( series->name, FR_VECT_8R, series->data->length, series->deltaT, seconds, units );
  if ( ! vect )
    XLAL_ERROR_NULL( func, XLAL_EERR ); /* "internal" error */
  vect->startX[0] = 0.0;

  memcpy( vect->data, series->data->data, series->data->length * sizeof( *series->data->data ) );

  return vect;
}


FrVect * XLALFrVectCOMPLEX8TimeSeries( COMPLEX8TimeSeries *series )
{
  static const char *func = "XLALFrVectCOMPLEX8TimeSeries";
  char seconds[LALUnitTextSize] = "s";
  char units[LALUnitTextSize];
  FrVect *vect;

  if ( NULL == XLALUnitAsString( units, sizeof( units ), &series->sampleUnits ) )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  vect = FrVectNew1D( series->name, FR_VECT_8C, series->data->length, series->deltaT, seconds, units );
  if ( ! vect )
    XLAL_ERROR_NULL( func, XLAL_EERR ); /* "internal" error */
  vect->startX[0] = 0.0;

  memcpy( vect->data, series->data->data, series->data->length * sizeof( *series->data->data ) );

  return vect;
}


/*  FIXME: Should f0 be startX in the vect or fShift in the proc data ??? */

FrVect * XLALFrVectCOMPLEX16TimeSeries( COMPLEX16TimeSeries *series )
{
  static const char *func = "XLALFrVectCOMPLEX8TimeSeries";
  char seconds[LALUnitTextSize] = "s";
  char units[LALUnitTextSize];
  FrVect *vect;

  if ( NULL == XLALUnitAsString( units, sizeof( units ), &series->sampleUnits ) )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  vect = FrVectNew1D( series->name, FR_VECT_16C, series->data->length, series->deltaT, seconds, units );
  if ( ! vect )
    XLAL_ERROR_NULL( func, XLAL_EERR ); /* "internal" error */
  vect->startX[0] = 0.0;

  memcpy( vect->data, series->data->data, series->data->length * sizeof( *series->data->data ) );

  return vect;
}


FrVect * XLALFrVectREAL4FrequencySeries( REAL4FrequencySeries *series )
{
  static const char *func = "XLALFrVectREAL4FrequencySeries";
  char hertz[LALUnitTextSize] = "s^-1";
  char units[LALUnitTextSize];
  FrVect *vect;

  if ( NULL == XLALUnitAsString( units, sizeof( units ), &series->sampleUnits ) )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  vect = FrVectNew1D( series->name, FR_VECT_4R, series->data->length, series->deltaF, hertz, units );
  if ( ! vect )
    XLAL_ERROR_NULL( func, XLAL_EERR ); /* "internal" error */
  vect->startX[0] = series->f0;

  memcpy( vect->data, series->data->data, series->data->length * sizeof( *series->data->data ) );

  return vect;
}


FrVect * XLALFrVectREAL8FrequencySeries( REAL8FrequencySeries *series )
{
  static const char *func = "XLALFrVectREAL8FrequencySeries";
  char hertz[LALUnitTextSize] = "s^-1";
  char units[LALUnitTextSize];
  FrVect *vect;

  if ( NULL == XLALUnitAsString( units, sizeof( units ), &series->sampleUnits ) )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  vect = FrVectNew1D( series->name, FR_VECT_8R, series->data->length, series->deltaF, hertz, units );
  if ( ! vect )
    XLAL_ERROR_NULL( func, XLAL_EERR ); /* "internal" error */
  vect->startX[0] = series->f0;

  memcpy( vect->data, series->data->data, series->data->length * sizeof( *series->data->data ) );

  return vect;
}


FrVect * XLALFrVectCOMPLEX8FrequencySeries( COMPLEX8FrequencySeries *series )
{
  static const char *func = "XLALFrVectCOMPLEX8FrequencySeries";
  char hertz[LALUnitTextSize] = "s^-1";
  char units[LALUnitTextSize];
  FrVect *vect;

  if ( NULL == XLALUnitAsString( units, sizeof( units ), &series->sampleUnits ) )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  vect = FrVectNew1D( series->name, FR_VECT_8C, series->data->length, series->deltaF, hertz, units );
  if ( ! vect )
    XLAL_ERROR_NULL( func, XLAL_EERR ); /* "internal" error */
  vect->startX[0] = series->f0;

  memcpy( vect->data, series->data->data, series->data->length * sizeof( *series->data->data ) );

  return vect;
}


FrVect * XLALFrVectCOMPLEX16FrequencySeries( COMPLEX16FrequencySeries *series )
{
  static const char *func = "XLALFrVectCOMPLEX16FrequencySeries";
  char hertz[LALUnitTextSize] = "s^-1";
  char units[LALUnitTextSize];
  FrVect *vect;

  if ( NULL == XLALUnitAsString( units, sizeof( units ), &series->sampleUnits ) )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  vect = FrVectNew1D( series->name, FR_VECT_8C, series->data->length, series->deltaF, hertz, units );
  if ( ! vect )
    XLAL_ERROR_NULL( func, XLAL_EERR ); /* "internal" error */
  vect->startX[0] = series->f0;

  memcpy( vect->data, series->data->data, series->data->length * sizeof( *series->data->data ) );

  return vect;
}


int XLALFrameAddCalRef( FrameH *frame, COMPLEX8FrequencySeries *series, int version, double duration )
{
  static const char *func = "XLALFrameAddCalRef";
  char representation[] = "calibration";
  char comment[] = "$Id$";
  char prefix[3];
  FrDetector *detector;
  FrStatData *sdat;
  FrVect *vect;
  int tstart;
  int tend;

  /* To find the detector, look for the detector prefix stored as the first
   * two characters of the channel name. */
  prefix[0] = series->name[0];
  prefix[1] = series->name[1];
  prefix[2] = 0;

  detector = FrameFindDetector( frame, prefix );
  if ( ! detector )
  {
    XLALPrintError( "XLAL Error - %s: no detector associated with prefix %s in frame\n", func, prefix );
    XLAL_ERROR( func, XLAL_EINVAL );
  }

  vect = XLALFrVectCOMPLEX8FrequencySeries( series );
  if ( ! vect )
    XLAL_ERROR( func, XLAL_EFUNC );

  tstart = series->epoch.gpsSeconds;
  tend = (int)ceil( series->epoch.gpsSeconds + 1e-9*series->epoch.gpsNanoSeconds + duration );
  sdat = FrStatDataNew( series->name, comment, representation, tstart, tend, version, vect, NULL );
  if ( ! sdat )
  {
    FrVectFree( vect );
    XLAL_ERROR( func, XLAL_EERR );
  }

  FrStatDataAdd( detector, sdat );

  return 0;
}


COMPLEX8FrequencySeries * XLALFrameGetCalRef( LIGOTimeGPS *validUntil, LIGOTimeGPS *epoch, const char *channel, FrameH *frame )
{
  static const char *func = "XLALFrameGetCalRef";
  char prefix[3];
  LALUnit unit;
  COMPLEX8FrequencySeries *series;
  FrStatData *sdat;
  char *chan;

  /* To find the detector, look for the detector prefix stored as the first
   * two characters of the channel name. */
  prefix[0] = channel[0];
  prefix[1] = channel[1];
  prefix[2] = 0;

  chan = XLALStringDuplicate( channel );
  sdat = FrameFindStatData( frame, prefix, chan, epoch->gpsSeconds );
  LALFree( chan );
  if ( ! sdat )
  {
    XLALPrintError( "XLAL Error - %s: no stat data channel %s for GPS time %d in frame\n", func, channel, epoch->gpsSeconds );
    XLAL_ERROR_NULL( func, XLAL_EINVAL );
  }

  /* parse units */
  if ( NULL == XLALParseUnitString( &unit, sdat->data->unitY ) )
  {
    XLALPrintError( "XLAL Error - %s: could not parse unit string %s\n", func, sdat->data->unitY );
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  }

  /* use validUntil as dummy structure */
  validUntil->gpsSeconds     = sdat->timeStart;
  validUntil->gpsNanoSeconds = 0;

  series = XLALCreateCOMPLEX8FrequencySeries( sdat->data->name, validUntil, sdat->data->startX[0], sdat->data->dx[0], &unit, sdat->data->nData );
  if ( ! series )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  memcpy( series->data->data, sdat->data->data, series->data->length * sizeof( *series->data->data ) );

  /* report end time of validity in validUntil */
  validUntil->gpsSeconds = sdat->timeEnd;

  return series;
}


int XLALFrameAddCalFac( FrameH *frame, REAL4TimeSeries *series )
{
  static const char *func = "XLALFrameAddCalFac";
  char comment[] = "$Id$";
  LIGOTimeGPS frameEpoch;
  FrProcData *proc;
  FrVect *vect;

  vect = XLALFrVectREAL4TimeSeries( series );
  if ( ! vect )
    XLAL_ERROR( func, XLAL_EFUNC );

  proc = FrProcDataNewV( frame, vect );
  if ( ! proc )
  {
    FrVectFree( vect );
    XLAL_ERROR( func, XLAL_EERR );
  }

  /* comment is rcs id */
  FrStrCpy( &proc->comment, comment );

  /* time offset: compute this from frame time */
  frameEpoch.gpsSeconds     = frame->GTimeS;
  frameEpoch.gpsNanoSeconds = frame->GTimeN;
  proc->timeOffset = XLALGPSDiff( &series->epoch, &frameEpoch );

  /* remaining metadata */
  proc->type    = 1;
  proc->subType = 0;
  proc->tRange  = series->data->length * series->deltaT;
  proc->fShift  = 0.0;
  proc->phase   = 0.0;
  proc->BW      = 0.0;

  return 0;
}


REAL4TimeSeries * XLALFrameGetCalFac( const char *channel, FrameH *frame )
{
  static const char *func = "XLALFrameGetCalFac";
  LIGOTimeGPS epoch;
  REAL4TimeSeries *series;
  FrProcData *proc;
  char *chan;

  chan = XLALStringDuplicate( channel );
  proc = FrProcDataFind( frame, chan );
  LALFree( chan );
  if ( ! proc )
  {
    XLALPrintError( "XLAL Error - %s: no proc data channel %s in frame\n", func, channel );
    XLAL_ERROR_NULL( func, XLAL_EINVAL );
  }

  epoch.gpsSeconds     = frame->GTimeS;
  epoch.gpsNanoSeconds = frame->GTimeN;
  XLALGPSAdd( &epoch, proc->timeOffset );

  series = XLALCreateREAL4TimeSeries( proc->data->name, &epoch, proc->fShift, proc->data->dx[0], &lalDimensionlessUnit, proc->data->nData );
  if ( ! series )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  memcpy( series->data->data, proc->data->data, series->data->length * sizeof( *series->data->data ) );

  return series;
}


LALCalData * XLALFrGetCalData( LIGOTimeGPS *epoch, const char *readoutChannel, const char *fname )
{
  static const char *func = "XLALFrGetCalData";
  char responseName[LALNameLength]  = "Xn:CAL-RESPONSE_";
  char cavgainName[LALNameLength]   = "Xn:CAL-CAV_GAIN_";
  char oloopgainName[LALNameLength] = "Xn:CAL-OLOOP_GAIN";
  char cavfacName[LALNameLength]    = "Xn:CAL-CAV_FAC";
  char oloopfacName[LALNameLength]  = "Xn:CAL-OLOOP_FAC";
  char readoutPoint[LALNameLength]; /* DARM_ERR or AS_Q */
  char ifo[3];
  char *fileName;
  LIGOTimeGPS tend;
  LALCalData *caldata;
  FrFile *frfile;
  FrameH *frame;

  /* setup channel names based on readout channel */
  XLALStringCopy( ifo, readoutChannel, sizeof( ifo ) );
  XLALStringCopy( readoutPoint, readoutChannel + strlen("Xn:LSC-"), sizeof( readoutPoint ) );
  XLALStringConcatenate( responseName, readoutPoint, sizeof( responseName ) );
  XLALStringConcatenate( cavgainName, readoutPoint, sizeof( cavgainName ) );
  memcpy( responseName, ifo, 2 );
  memcpy( cavgainName, ifo, 2 );
  memcpy( oloopgainName, ifo, 2 );
  memcpy( cavfacName, ifo, 2 );
  memcpy( oloopfacName, ifo, 2 );

  if ( strcmp( readoutPoint, "AS_Q" ) && strcmp( readoutPoint, "DARM_ERR" ) )
    XLAL_ERROR_NULL( func, XLAL_ENAME );
  

  caldata = LALCalloc( 1, sizeof( *caldata ) );
  if ( ! caldata )
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );

  fileName = XLALStringDuplicate( fname );
  if ( ! fileName )
  {
    XLALDestroyCalData( caldata );
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  }
  frfile = FrFileINew( fileName );
  XLALFree( fileName );
  if ( ! frfile )
  {
    XLALDestroyCalData( caldata );
    XLAL_ERROR_NULL( func, XLAL_EERR );
  }
  frame = FrameRead( frfile );
  if ( ! frame )
  {
    FrFileIEnd( frfile );
    XLALDestroyCalData( caldata );
    XLAL_ERROR_NULL( func, XLAL_EERR );
  }

  caldata->cavityFactors = XLALFrameGetCalFac( cavfacName, frame );
  if ( ! caldata->cavityFactors )
  {
    FrFileIEnd( frfile );
    XLALDestroyCalData( caldata );
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  }
  caldata->openLoopFactors = XLALFrameGetCalFac( oloopfacName, frame );
  if ( ! caldata->openLoopFactors )
  {
    FrFileIEnd( frfile );
    XLALDestroyCalData( caldata );
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  }
  caldata->responseReference = XLALFrameGetCalRef( &tend, epoch, responseName, frame );
  if ( ! caldata->responseReference )
  {
    FrFileIEnd( frfile );
    XLALDestroyCalData( caldata );
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  }
  caldata->cavityGainReference = XLALFrameGetCalRef( &tend, epoch, cavgainName, frame );
  if ( ! caldata->cavityGainReference )
  {
    FrFileIEnd( frfile );
    XLALDestroyCalData( caldata );
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  }
  caldata->openLoopGainReference = XLALFrameGetCalRef( &tend, epoch, oloopgainName, frame );
  if ( ! caldata->openLoopGainReference )
  {
    FrFileIEnd( frfile );
    XLALDestroyCalData( caldata );
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  }

  FrFileIEnd( frfile );
  return caldata;
}
