/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Lisa M. Goggin, Patrick Brady
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

#include <string.h>

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/GenerateRing.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Units.h>
#include <lal/FindChirp.h>

#include "lalapps.h"
#include "injsgnl.h"
#include "getresp.h"
#include "errutil.h"

RCSID( "$Id$" );

/* maximum length of filename */
#define FILENAME_LENGTH 255

/* routine to inject a signal with parameters read from a LIGOLw-format file */
int inject_signal( 
    REAL4TimeSeries   *series, 
    int                injectSignalType, 
    const char        *injectFile, 
    const char        *calCacheFile, 
    REAL4              responseScale, 
    const char        *channel_name 
    )
{
  /* note: duration is only used for response, and can be relatively coarse */
  const      REAL8 duration = 16; /* determines deltaF=1/dataDuration Hz*/
  LALStatus                status     = blank_status;
  COMPLEX8FrequencySeries *response   = NULL;
  SimInspiralTable        *injectList = NULL;
  SimInspiralTable        *thisInject;
  SimRingdownTable        *ringList = NULL;
  char                     injFile[FILENAME_LENGTH + 1];
  LIGOTimeGPS              epoch;
  UINT4                    numInject = 0;
  INT4                     startSec;
  INT4                     stopSec;
  int                      strainData;
  char                     ifoName[3];
  char                     name[LALNameLength];
  REAL8                    sampleRate;
  INT4                     calType=0;

  /* xml output data */
  CHAR                  fname[FILENAME_MAX];
  MetadataTable         ringinjections;
  LIGOLwXMLStream       xmlfp;
  
  /* copy injectFile to injFile (to get rid of const qual) */
  strncpy( injFile, injectFile, sizeof( injFile ) - 1 );
  snprintf( name, sizeof( name ), "%s_INJ", series->name );
  strncpy( ifoName, series->name, 2 );
  ifoName[2] = 0;

  /* get list of injections for this data epoch */
  verbose( "reading simulated-ring tables from file %s\n", injFile );
  startSec = series->epoch.gpsSeconds;
  stopSec  = startSec + ceil( 1e-9 * series->epoch.gpsNanoSeconds
      + series->deltaT * series->data->length );

/* call the approprate LAL injection routine */
  switch ( injectSignalType )
  {
    case ring_inject:
      ringList = 
      XLALSimRingdownTableFromLIGOLw( injFile, startSec, stopSec );
      numInject  = 0;
      for ( thisInject = ringList; thisInject; thisInject = thisInject->next )
            ++numInject;
    break;
    case imr_inject: case imr_ring_inject: case EOBNR_inject: case Phenom_inject:
      numInject = 
        SimInspiralTableFromLIGOLw( &injectList, injFile, startSec, stopSec );
      break;
    default:
      error( "unrecognized injection signal type\n" );
  }

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
    strainData = !XLALUnitCompare( &series->sampleUnits, &lalStrainUnit );

    /* determine sample rate of data (needed for response) */
    sampleRate = 1.0/series->deltaT;

    /* this gets an impulse response if we have strain data */
    response = get_response( calCacheFile, ifoName, &epoch, duration,
        sampleRate, responseScale, strainData, channel_name );

    /* units must be counts for inject; reset below if they were strain */
    series->sampleUnits = lalADCCountUnit;
   
   /* inject the signals */
    verbose( "injecting %u signal%s into time series\n", numInject,
        numInject == 1 ? "" : "s" );

    switch ( injectSignalType )
    {
      case ring_inject:
        LAL_CALL( LALRingInjectSignals(&status, series, ringList, response, calType),
            &status );
        break;
      case imr_inject: case imr_ring_inject:
        ringList = (SimRingdownTable *) XLALCalloc( 1, sizeof(SimRingdownTable) );
        LAL_CALL( LALFindChirpInjectIMR( &status, series, injectList, ringList, 
              response, injectSignalType ), &status );
        break;        
      case EOBNR_inject: case Phenom_inject:
        LAL_CALL( LALFindChirpInjectSignals( &status, series, injectList, response ), &status );
        break;
      default:
        error( "unrecognized injection signal type\n" );
    }

    /* correct the name */
    strncpy( series->name, name, sizeof( series->name ) - 1 );

    /* reset units if necessary */
    if ( strainData )
      series->sampleUnits = lalStrainUnit;

    switch ( injectSignalType )
    {
      case imr_inject: case imr_ring_inject:
        /* write output to LIGO_LW XML file    */
        ringinjections.simRingdownTable = NULL;
        ringinjections.simRingdownTable = ringList;
                 
        /* create the output file name */
        snprintf( fname, sizeof(fname), "HL-INJECTIONS_0-%d-%d.xml", 
            startSec, stopSec - startSec );
        fprintf( stdout, "Writing the injection details to %s\n", fname);
    
        /* open the xml file */
        memset( &xmlfp, 0, sizeof(LIGOLwXMLStream) );
        LALOpenLIGOLwXMLFile( &status, &xmlfp, fname);
    
        /* write the sim_ringdown table */
        if ( ringinjections.simRingdownTable )
        {
          LALBeginLIGOLwXMLTable( &status, &xmlfp, sim_ringdown_table );
          LALWriteLIGOLwXMLTable( &status, &xmlfp, ringinjections,
            sim_ringdown_table );
          LALEndLIGOLwXMLTable ( &status, &xmlfp );
        }
    
        while ( ringinjections.simRingdownTable->next )
        {
          ringList=ringinjections.simRingdownTable;
          ringinjections.simRingdownTable = ringinjections.simRingdownTable->next;
          LALFree( ringList );
        }
    
        /* close the injection file */
        LALCloseLIGOLwXMLFile ( &status, &xmlfp );
    
        /* free memory */

        while ( injectList )
        {
          thisInject = injectList;
          injectList = injectList->next;
          LALFree( thisInject );
        }
        break;
    }
    
    /* free memory */
    while ( ringList )
    {
      ringinjections.simRingdownTable = ringList;
      ringList = ringList->next;
      LALFree( ringinjections.simRingdownTable );
    }
    
    XLALDestroyCOMPLEX8Vector( response->data );
    LALFree( response );
  }
  
  return 0;
}
