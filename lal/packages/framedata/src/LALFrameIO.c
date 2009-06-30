/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton, Robert Adam Mercer, Xavier Siemens
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

#include <config.h>
#include <unistd.h>
#ifndef HAVE_GETHOSTNAME_PROTOTYPE
int gethostname(char *name, int len);
#endif

#include <math.h>

#include <lal/LALDetectors.h>
#include <lal/LALStdio.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>

#include <lal/LALString.h>
#include <lal/FrameCache.h>
#include <lal/LALFrameIO.h>
#include <lal/LALCalibration.h>

#include <lal/LALRCSID.h>
NRCSID (LALFRAMEIOC,"$Id$");

/* FIXME: WARNING: this value might need to change in the future */
#define FR_FILE_HEADER_SIZE 40 /* size of frame file header in bytes */

struct FrFile * XLALFrOpenURL( const char *url )
{
  static const char *func = "XLALFrOpenURL";
  struct FrFile *frfile = NULL;
  char prot[FILENAME_MAX] = "";
  char host[FILENAME_MAX] = "";
  char path[FILENAME_MAX] = "";
  int n;

  if ( ! url )
    XLAL_ERROR_NULL( func, XLAL_EFAULT );
  if ( strlen( url ) >= FILENAME_MAX )
  {
    XLALPrintError( "XLAL Error - %s: URL too long: %s\n", func, url );
    XLAL_ERROR_NULL( func, XLAL_EBADLEN );
  }

  n = sscanf( url, "%[^:]://%[^/]%s", prot, host, path );
  if ( n != 3 ) /* perhaps the hostname has been omitted */
  {
    strncpy( host, "localhost", sizeof( host ) - 1 );
    n = sscanf( url, "%[^:]://%s", prot, path );
    if ( n != 2 ) /* assume the whole thing is a file path */
    {
      strncpy( prot, "file", sizeof( prot ) - 1 );
      strncpy( path, url, sizeof( path ) - 1 );
    }
  }

  if ( strcmp( prot, "file" ) ) /* not a file URL */
  {
    XLALPrintError( "XLAL Error - %s: unsupported protocol %s\n", func, prot );
    XLAL_ERROR_NULL( func, XLAL_EINVAL );
  }
  else /* OK: this is a file URL */
  {
    if ( strcmp( host, "localhost" ) ) /* not explicitly localhost */
    { /* make sure that the host is the localhost */
      char localhost[FILENAME_MAX];
      gethostname( localhost, FILENAME_MAX - 1 );
      if ( strcmp( host, localhost ) ) /* not localhost */
      {
        XLALPrintError( "XLAL Error - %s: cannot read files from remote host %s\n", func, host );
        XLAL_ERROR_NULL( func, XLAL_EINVAL );
      }
    }
    frfile = FrFileINew( path );
    if ( ! frfile )
    {
      XLALPrintError( "XLAL Error - %s: could not open frame file %s\n", func, path );
      XLAL_ERROR_NULL( func, XLAL_EIO );
    }
  }

  return frfile;
}


/* code taken from the FrCheck program in the FrameL library */
int XLALFrFileCheckSum( FrFile *iFile )
{
  static const char *func = "XLALFrFileCheckSum";
  FrameH *frame = NULL;
  FRBOOL chkSumFiFlag = iFile->chkSumFiFlag;
  FRBOOL chkSumFrFlag = iFile->chkSumFrFlag;
  iFile->chkSumFiFlag = FR_YES;
  iFile->chkSumFrFlag = FR_NO;
  /* sequential read */
  /* FIXME: WARNING: HACK: following line may need to be updated in future */
  FrIOSet(iFile->frfd, FR_FILE_HEADER_SIZE);
  while ((frame = FrameReadRecycle(iFile, frame)))
    ;
  iFile->chkSumFiFlag = chkSumFiFlag;
  iFile->chkSumFrFlag = chkSumFrFlag;
  if ( iFile->error != FR_OK ) {
    XLALPrintError( "XLAL Error - %s: %s\n", func, FrErrorGetHistory() );
    XLAL_ERROR( func, XLAL_EFAILED );
    return -1;
  }
  if ( iFile->chkTypeFiRead == 0 ) {
    XLALPrintWarning( "XLAL Warning - %s: missing checksum\n", func );
    return 1;
  }
  if ( iFile->chkSumFiRead != iFile->chkSumFi ) {
    XLALPrintError( "XLAL Error - %s: bad checksum\n", func );
    return -1;
  }
  return 0;
}


/* This routine is essentially the same as FrHistoryNew but the name of the
 * history can be set too. */
FrHistory * XLALFrHistoryAdd( FrameH *frame, const char *name, const char *comment )
{
  static const char *func = "XLALFrHistoryAdd";
  union { const char *cs; char *s; } namecnvrt; /* get rid of const qual */
  union { const char *cs; char *s; } commentcnvrt; /* get rid of const qual */
  LIGOTimeGPS now;
  FrHistory *history;

  /* get current time */
  if ( ! XLALGPSTimeNow( &now ) )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  /* this nonsense is to convert const char * to char * ... don't worry,
   * the frame library just copies the string anyway */
  namecnvrt.cs    = name;
  commentcnvrt.cs = comment;

  /* now create the history */
  history = FrHistoryNew( namecnvrt.s, now.gpsSeconds, commentcnvrt.s );
  if ( ! history )
    XLAL_ERROR_NULL( func, XLAL_EERR ); /* "internal" error */

  /* attach history to the frame structure */
  if ( frame )
  {
    /* behaviour is identical to FrHistoryAdd if name is NULL */
    if ( ! name )
      FrStrCpy( &history->name, frame->name );
    history->next = frame->history;
    frame->history = history;
  }

  return history;
}

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
  static char histidname[] = __FILE__ " Id";
  static char histtagname[] = __FILE__ " Tag";
  static char rcsname[] = "$Name$";
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

  /* add history: name of history field is this function's name */
  if ( ! XLALFrHistoryAdd( frame, histidname, rcsid ) )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  if ( ! XLALFrHistoryAdd( frame, histtagname, rcsname ) )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return frame;
}


FrVect * XLALFrVectINT4TimeSeries( INT4TimeSeries *series )
{
  static const char *func = "XLALFrVectINT4TimeSeries";
  char seconds[LALUnitTextSize] = "s";
  char units[LALUnitTextSize];
  FrVect *vect;

  if ( NULL == XLALUnitAsString( units, sizeof( units ), &series->sampleUnits ) )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  vect = FrVectNew1D( series->name, FR_VECT_4S, series->data->length, series->deltaT, seconds, units );
  if ( ! vect )
    XLAL_ERROR_NULL( func, XLAL_EERR ); /* "internal" error */
  vect->startX[0] = 0.0;

  memcpy( vect->data, series->data->data, series->data->length * sizeof( *series->data->data ) );

  FrVectCompress (vect, 8, 0);
  if (vect->compress == 0) FrVectCompress (vect, 6, 1);

  return vect;
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

  FrVectCompress (vect, 8, 0);
  if (vect->compress == 0) FrVectCompress (vect, 6, 1);

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

  FrVectCompress (vect, 6, 1);

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

  FrVectCompress (vect, 6, 1);

  return vect;
}



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

  FrVectCompress (vect, 6, 1);

  return vect;
}

/*  FIXME: Should f0 be startX in the vect or fShift in the proc data ??? */

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

  FrVectCompress (vect, 8, 0);
  if (vect->compress == 0) FrVectCompress (vect, 6, 1);

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

  FrVectCompress (vect, 6, 1);

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

  FrVectCompress (vect, 6, 1);

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

  FrVectCompress (vect, 6, 1);

  return vect;
}


int XLALFrameAddCalRef( FrameH *frame, COMPLEX8FrequencySeries *series, int version, double duration )
{
  static const char *func = "XLALFrameAddCalRef";
  char representation[] = "freq_series";
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
  if ( ! sdat || sdat->timeStart > (UINT4)epoch->gpsSeconds || sdat->timeEnd < (UINT4)epoch->gpsSeconds )
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


int XLALFrameAddCalFac( FrameH *frame, REAL4TimeSeries *series, int version )
{
  static const char *func = "XLALFrameAddCalFac";
  char representation[] = "time_series";
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

  vect = XLALFrVectREAL4TimeSeries( series );
  if ( ! vect )
    XLAL_ERROR( func, XLAL_EFUNC );

  tstart = series->epoch.gpsSeconds;
  tend = (int)ceil( series->epoch.gpsSeconds + 1e-9*series->epoch.gpsNanoSeconds + series->data->length * series->deltaT );
  sdat = FrStatDataNew( series->name, comment, representation, tstart, tend, version, vect, NULL );
  if ( ! sdat )
  {
    FrVectFree( vect );
    XLAL_ERROR( func, XLAL_EERR );
  }

  FrStatDataAdd( detector, sdat );

  return 0;
}


REAL4TimeSeries * XLALFrameGetCalFac( LIGOTimeGPS *epoch, const char *channel, FrameH *frame )
{
  static const char *func = "XLALFrameGetCalFac";
  char prefix[3];
  LIGOTimeGPS tmpEpoch;
  REAL4TimeSeries *series;
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
  if ( ! sdat || sdat->timeStart > (UINT4)epoch->gpsSeconds || sdat->timeEnd < (UINT4)epoch->gpsSeconds )
  {
    XLALPrintError( "XLAL Error - %s: no stat data channel %s for GPS time %d in frame\n", func, channel, epoch->gpsSeconds );
    XLAL_ERROR_NULL( func, XLAL_EINVAL );
  }

  tmpEpoch.gpsSeconds     = sdat->timeStart;
  tmpEpoch.gpsNanoSeconds = 0;
  series = XLALCreateREAL4TimeSeries( sdat->data->name, &tmpEpoch, 0.0, sdat->data->dx[0], &lalDimensionlessUnit, sdat->data->nData );
  if ( ! series )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  memcpy( series->data->data, sdat->data->data, series->data->length * sizeof( *series->data->data ) );

  return series;
}

int XLALFrameAddREAL8TimeSeriesProcData( FrameH *frame, REAL8TimeSeries *series )
{
	static const char * func = "XLALFrameAddREAL8TimeSeriesProcData";
	char rcsinfo[]  = "$Id$" "$Name$";
	LIGOTimeGPS frameEpoch;
	FrProcData *proc;
	FrVect *vect;
	REAL8 duration;

	duration = series->deltaT * series->data->length;

	vect = XLALFrVectREAL8TimeSeries( series );
	if ( ! vect )
		XLAL_ERROR( func, XLAL_EFUNC );

	proc = FrProcDataNewV( frame, vect );
	if ( ! proc ) {
		FrVectFree( vect );
		XLAL_ERROR( func, XLAL_EERR );
	}

	/* comment is rcs info of this routine */
/* 	FrStrCpy( &proc->comment, rcsinfo ); */

	/* time offset: compute this from frame time */
	frameEpoch.gpsSeconds     = frame->GTimeS;
	frameEpoch.gpsNanoSeconds = frame->GTimeN;
	proc->timeOffset = XLALGPSDiff( &series->epoch, &frameEpoch );

	/* remaining metadata */
	proc->type    = 1;
	proc->subType = 0;
	proc->tRange  = duration;
	proc->fShift  = 0.0;
	proc->phase   = 0.0;
	proc->BW      = 0.0;

	return 0;
}

int XLALFrameAddREAL4TimeSeriesProcData( FrameH *frame, REAL4TimeSeries *series )
{
	static const char * func = "XLALFrameAddREAL4TimeSeriesProcData";
	char rcsinfo[]  = "$Id$" "$Name$";
	LIGOTimeGPS frameEpoch;
	FrProcData *proc;
	FrVect *vect;
	REAL8 duration;

	duration = series->deltaT * series->data->length;

	vect = XLALFrVectREAL4TimeSeries( series );
	if ( ! vect )
		XLAL_ERROR( func, XLAL_EFUNC );

	proc = FrProcDataNewV( frame, vect );
	if ( ! proc ) {
		FrVectFree( vect );
		XLAL_ERROR( func, XLAL_EERR );
	}

	/* comment is rcs info of this routine */
/* 	FrStrCpy( &proc->comment, rcsinfo ); */

	/* time offset: compute this from frame time */
	frameEpoch.gpsSeconds     = frame->GTimeS;
	frameEpoch.gpsNanoSeconds = frame->GTimeN;
	proc->timeOffset = XLALGPSDiff( &series->epoch, &frameEpoch );

	/* remaining metadata */
	proc->type    = 1;
	proc->subType = 0;
	proc->tRange  = duration;
	proc->fShift  = 0.0;
	proc->phase   = 0.0;
	proc->BW      = 0.0;

	return 0;
}



int XLALFrameAddINT4TimeSeriesProcData( FrameH *frame, INT4TimeSeries *series )
{
	static const char * func = "XLALFrameAddINT4TimeSeriesProcData";
	char rcsinfo[]  = "$Id$" "$Name$";
	LIGOTimeGPS frameEpoch;
	FrProcData *proc;
	FrVect *vect;
	REAL8 duration;

	duration = series->deltaT * series->data->length;

	vect = XLALFrVectINT4TimeSeries( series );
	if ( ! vect )
		XLAL_ERROR( func, XLAL_EFUNC );

	proc = FrProcDataNewV( frame, vect );
	if ( ! proc ) {
		FrVectFree( vect );
		XLAL_ERROR( func, XLAL_EERR );
	}

	/* comment is rcs info of this routine */
/* 	FrStrCpy( &proc->comment, rcsinfo ); */

	/* time offset: compute this from frame time */
	frameEpoch.gpsSeconds     = frame->GTimeS;
	frameEpoch.gpsNanoSeconds = frame->GTimeN;
	proc->timeOffset = XLALGPSDiff( &series->epoch, &frameEpoch );

	/* remaining metadata */
	proc->type    = 1;
	proc->subType = 0;
	proc->tRange  = duration;
	proc->fShift  = 0.0;
	proc->phase   = 0.0;
	proc->BW      = 0.0;

	return 0;
}


int XLALFrameAddREAL4TimeSeriesSimData( FrameH *frame, REAL4TimeSeries *series )
{
	static const char * func = "XLALFrameAddREAL4TimeSeriesSimData";
	char rcsinfo[]  = "$Id$" "$Name$";
  LIGOTimeGPS frameEpoch;
	FrSimData *sim;
	FrVect *vect;

	vect = XLALFrVectREAL4TimeSeries( series );
	if ( ! vect )
		XLAL_ERROR( func, XLAL_EFUNC );

	sim = FrSimDataNew( frame, series->name, series->deltaT, series->data->length, -32 );
	if ( ! sim ) {
		FrVectFree( vect );
		XLAL_ERROR( func, XLAL_EERR );
	}
  FrVectFree( sim->data );
  sim->data = vect;

	/* comment is rcs info of this routine */
/* 	FrStrCpy( &proc->comment, rcsinfo ); */

  /* time offset: compute this from frame time */
  frameEpoch.gpsSeconds     = frame->GTimeS;
  frameEpoch.gpsNanoSeconds = frame->GTimeN;
  sim->timeOffset = XLALGPSDiff( &series->epoch, &frameEpoch );

  /* remaining metadata */
  sim->fShift = 0;
  sim->phase  = 0;

	return 0;
}

int XLALFrameAddREAL4TimeSeriesAdcData( FrameH *frame, REAL4TimeSeries *series )
{
	static const char * func = "XLALFrameAddREAL4TimeSeriesAdcData";
	char rcsinfo[]  = "$Id$" "$Name$";
	LIGOTimeGPS frameEpoch;
	FrAdcData *adc;
	int i;

 	adc = FrAdcDataNew( frame, series->name, 1./series->deltaT, series->data->length, -32 );
/* 	adc = FrAdcDataNewF (frame, series->name, rcsinfo, 1, 1, -32, 0, 1, "Counts", 1.0/series->deltaT, series->data->length); */
	if ( ! adc ) {
	  XLAL_ERROR( func, XLAL_EERR );
	}

	for(i=0; i < (int)series->data->length; i++) {
	  adc->data->dataF[i] = series->data->data[i];
/* 	  fprintf(stdout,"%d %f %f\n",i, adc->data->dataF[i], series->data->data[i]); */
	}

	FrVectCompress (adc->data, 8, 0);
	if (adc->data->compress == 0) FrVectCompress (adc->data, 6, 1);

	return 0;
}

#if 0
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
#endif


LALCalData * XLALFrameGetCalData( LIGOTimeGPS *epoch, const char *readoutChannel, FrameH *frame )
{
  static const char *func = "XLALFrameGetCalData";
  LALCalData *caldata;
  char oloopgainName[LALNameLength]     = "Xn:CAL-OLOOP_GAIN";
  char cavgainName[LALNameLength]       = "Xn:CAL-CAV_GAIN_";
  char responseName[LALNameLength]      = "Xn:CAL-RESPONSE_";
  char actuationName[LALNameLength]     = "Xn:CAL-ACTUATION";
  char digfltName[LALNameLength]        = "Xn:CAL-DIGFLT_";
  char cavfacName[LALNameLength]        = "Xn:CAL-CAV_FAC";
  char oloopfacName[LALNameLength]      = "Xn:CAL-OLOOP_FAC";
  char readoutPoint[LALNameLength]; /* DARM_ERR or AS_Q */
  char ifo[3];
  LIGOTimeGPS tend;
  int errnum;

  /* setup channel names based on readout channel */
  XLALStringCopy( ifo, readoutChannel, sizeof( ifo ) );
  XLALStringCopy( readoutPoint, readoutChannel + strlen("Xn:LSC-"), sizeof( readoutPoint ) );
  XLALStringConcatenate( responseName, readoutPoint, sizeof( responseName ) );
  XLALStringConcatenate( cavgainName, readoutPoint, sizeof( cavgainName ) );
  XLALStringConcatenate( digfltName, readoutPoint, sizeof( digfltName ) );
  memcpy( oloopgainName, ifo, 2 );
  memcpy( actuationName, ifo, 2 );
  memcpy( responseName, ifo, 2 );
  memcpy( cavgainName, ifo, 2 );
  memcpy( digfltName, ifo, 2 );
  memcpy( cavfacName, ifo, 2 );
  memcpy( oloopfacName, ifo, 2 );

  if ( strcmp( readoutPoint, "AS_Q" ) && strcmp( readoutPoint, "DARM_ERR" ) )
    XLAL_ERROR_NULL( func, XLAL_ENAME );

  caldata = LALCalloc( 1, sizeof( *caldata ) );
  if ( ! caldata )
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );

  /* caldata->cavityFactors = XLALFrameGetCalFac( cavfacName, frame ); */
  caldata->cavityFactors = XLALFrameGetCalFac( epoch, cavfacName, frame );
  if ( ! caldata->cavityFactors )
  {
    XLALDestroyCalData( caldata );
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  }
  /* caldata->openLoopFactors = XLALFrameGetCalFac( oloopfacName, frame ); */
  caldata->openLoopFactors = XLALFrameGetCalFac( epoch, oloopfacName, frame );
  if ( ! caldata->openLoopFactors )
  {
    XLALDestroyCalData( caldata );
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  }
  caldata->responseReference = XLALFrameGetCalRef( &tend, epoch, responseName, frame );
  if ( ! caldata->responseReference )
  {
    XLALDestroyCalData( caldata );
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  }
  caldata->cavityGainReference = XLALFrameGetCalRef( &tend, epoch, cavgainName, frame );
  if ( ! caldata->cavityGainReference )
  {
    XLALDestroyCalData( caldata );
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  }
  caldata->openLoopGainReference = XLALFrameGetCalRef( &tend, epoch, oloopgainName, frame );
  if ( ! caldata->openLoopGainReference )
  {
    XLALDestroyCalData( caldata );
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  }
  XLAL_TRY( caldata->actuationReference = XLALFrameGetCalRef( &tend, epoch, actuationName, frame ), errnum );
  XLAL_TRY( caldata->digitalFilterReference = XLALFrameGetCalRef( &tend, epoch, digfltName, frame ), errnum );

  return caldata;
}


LALCalData * XLALFrGetCalData( LIGOTimeGPS *epoch, const char *readoutChannel, const char *fname )
{
  static const char *func = "XLALFrGetCalData";
  LALCalData *caldata;
  char *fileName;
  FrFile *frfile;
  FrameH *frame;

  fileName = XLALStringDuplicate( fname );
  if ( ! fileName )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  frfile = FrFileINew( fileName );
  XLALFree( fileName );
  if ( ! frfile )
    XLAL_ERROR_NULL( func, XLAL_EERR );
  frame = FrameRead( frfile );
  if ( ! frame )
  {
    FrFileIEnd( frfile );
    XLAL_ERROR_NULL( func, XLAL_EERR );
  }
  caldata = XLALFrameGetCalData( epoch, readoutChannel, frame );
  if ( ! caldata )
  {
    FrFileIEnd( frfile );
    XLALDestroyCalData( caldata );
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  }

  FrFileIEnd( frfile );
  return caldata;
}


LALCalData * XLALFrCacheGetCalData( LIGOTimeGPS *epoch, const char *readoutChannel, FrCache *cache )
{
  static const char *func = "XLALFrCacheGetCalData";
  LALCalData *caldata;
  char srcRegEx[] = "X";
  char dscRegEx[] = "Xn_CAL_";
  FrCacheSieve sieve;
  FrFile *frfile;
  FrameH *frame;

  memset( &sieve, 0, sizeof(sieve) );
  memcpy( srcRegEx, readoutChannel, 1 ); /* copy site */
  memcpy( dscRegEx, readoutChannel, 2 ); /* copy ifo */

  /* first sieve: get all cal frames containing the specified time */
  sieve.srcRegEx = srcRegEx; /* those files for the correct site */
  sieve.dscRegEx = dscRegEx; /* those files for the correct ifo */
  /* only want frame files that contain the specified time */
  sieve.earliestTime = epoch->gpsSeconds;
  sieve.latestTime = epoch->gpsSeconds;
  cache = XLALFrSieveCache( cache, &sieve );
  if ( ! cache )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  if ( cache->numFrameFiles < 1 )
  {
    XLALPrintError( "XLAL Error - %s: no matching calibration frame files found in cache %s\n", func );
    XLALFrDestroyCache( cache );
    XLAL_ERROR_NULL( func, XLAL_EIO );
  }

  /*
   * We want to open the LAST of the remaining frames!
   * This frame will have the highest version number and
   * it will be the latest one to contain the specified
   * time (in case the files overlap).
   */
  frfile = XLALFrOpenURL( cache->frameFiles[cache->numFrameFiles - 1].url );
  XLALFrDestroyCache( cache ); /* don't need this anymore */
  if ( ! frfile )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  frame = FrameRead( frfile );
  if ( ! frame )
  {
    FrFileIEnd( frfile );
    XLAL_ERROR_NULL( func, XLAL_EERR );
  }
  caldata = XLALFrameGetCalData( epoch, readoutChannel, frame );
  if ( ! caldata )
  {
    FrFileIEnd( frfile );
    XLALDestroyCalData( caldata );
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  }

  FrFileIEnd( frfile );
  return caldata;
}


int XLALFrameWrite(FrameH *frame, const char *fname, int compressLevel)
{
  static const char *func = "XLALFrameWrite";
  FrFile *frfile;
  char tmpfname[FILENAME_MAX];
  int c;

  /* set temporary filename */
  c = LALSnprintf( tmpfname, sizeof( tmpfname ), "%s.tmp", fname );
  if ( c < 0 || c > (int)sizeof(tmpfname) - 2 )
    XLAL_ERROR( func, XLAL_ENAME );

  /* open temporary file */
  frfile = FrFileONew( tmpfname, compressLevel );
  if ( !frfile )
    XLAL_ERROR( func, XLAL_EIO );

  /* write frame to temporary filename */
  FrameWrite( frame, frfile );
  FrFileOEnd( frfile );

  /* rename tmpfile */
  if (rename(tmpfname, fname) < 0)
    XLAL_ERROR( func, XLAL_ESYS );

  return 0;
}
