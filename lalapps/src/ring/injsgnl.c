#include <string.h>

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/GenerateBurst.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/Units.h>

#include "lalapps.h"
#include "injsgnl.h"
#include "getresp.h"
#include "errutil.h"

RCSID( "$Id$" );

#define FILENAME_LENGTH 255

int inject_signal( REAL4TimeSeries *series, int injectSignalType, 
    const char *injectFile, const char *calCacheFile, REAL4 responseScale )
{
  /* note: duration is only used for response, and can be relatively coarse */
  const REAL8 duration = 16; /* determines deltaF=1/dataDuration Hz*/
  LALStatus                status     = blank_status;
  COMPLEX8FrequencySeries *response   = NULL;
  SimBurstTable           *injectList = NULL;
  SimBurstTable           *thisInject;
  char                     injFile[FILENAME_LENGTH + 1];
  LIGOTimeGPS              epoch;
  UINT4                    numInject;
  INT4                     startSec;
  INT4                     stopSec;
  int                      strainData;
  char                     ifoName[3];
  char                     name[LALNameLength];
  REAL8                    sampleRate;

  /* copy injectFile to injFile (to get rid of const qual) */
  strncpy( injFile, injectFile, sizeof( injFile ) - 1 );
  LALSnprintf( name, sizeof( name ), "%s_INJ", series->name );
  strncpy( ifoName, series->name, 2 );
  ifoName[2] = 0;

  /* get list of injections for this data epoch */
  verbose( "reading simulated-burst tables from file %s\n", injFile );
  startSec = series->epoch.gpsSeconds;
  stopSec  = startSec + ceil( 1e-9 * series->epoch.gpsNanoSeconds
      + series->deltaT * series->data->length );

  switch ( injectSignalType )
  {
    case burst_inject:
      LAL_CALL( LALSimBurstTableFromLIGOLw( &status, &injectList, injFile,
            startSec, stopSec ), &status );
      break;
    default:
      error( "unrecognized injection signal type\n" );
  }

  /* count the number of injections */
  numInject  = 0;
  for ( thisInject = injectList; thisInject; thisInject = thisInject->next )
    ++numInject;

  if ( numInject == 0 )
    verbose( "no injections to perform in this data epoch\n" );
  else /* perform the injections */
  {
    /* get a representative response function */
    epoch.gpsSeconds     = startSec;
    epoch.gpsNanoSeconds = 0;
    verbose( "getting response function for GPS time %d.%09d\n",
        epoch.gpsSeconds, epoch.gpsNanoSeconds );

    /* determine if this is strain data */
    strainData = XLALUnitCompare( &series->sampleUnits, &lalStrainUnit );

    /* determine sample rate of data (needed for response) */
    sampleRate = 1.0/series->deltaT;

    /* this gets an impulse response if we have strain data */
    response = get_response( calCacheFile, ifoName, &epoch, duration,
        sampleRate, responseScale, strainData );

    /* units must be counts for inject; reset below if they were strain */
    series->sampleUnits = lalADCCountUnit;

    /* inject the signals */
    verbose( "injecting %u signal%s into time series\n", numInject,
        numInject == 1 ? "" : "s" );
    switch ( injectSignalType )
    {
      case burst_inject:
        LAL_CALL( LALBurstInjectSignals(&status, series, injectList, response),
            &status );
        break;
      default:
        error( "unrecognized injection signal type\n" );
    }

    /* correct the name */
    strncpy( series->name, name, sizeof( series->name ) - 1 );

    /* reset units if necessary */
    if ( strainData )
      series->sampleUnits = lalStrainUnit;

    /* free memory */
    while ( injectList )
    {
      thisInject = injectList;
      injectList = injectList->next;
      LALFree( thisInject );
    }

    XLALDestroyCOMPLEX8Vector( response->data );
    LALFree( response );
  }

  return 0;
}
