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

/**
 * \author Jolien D. E. Creighton
 * \file
 *
 * \heading{Program \ref Read40mData.c}
 *
 * Tests the low-level frame stream routines by reading Caltech 40m-prototype
 * data from Nov1994.
 *
 * \heading{Usage}
 *
 * \code
 * Read40mData
 * \endcode
 *
 * \heading{Description}
 *
 * This program reads the channels \c IFO_DMRO and \c IFO_Lock from
 * all the frames in the directory set in the environment \c LAL_FRAME_PATH
 * (or the current directory if this environment is not set) and prints them
 * to new frame files containing only the in-lock \c IFO_DMRO channel.
 *
*/

#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/PrintFTSeries.h>
#include <lal/FrameStream.h>

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 1; } else ((void)0)

#define CHANNEL "IFO_DMRO"


int main( void )
{
  static LALStatus status;
  const float duration = 60.0; /* seconds of data to read at a time */
  const short lockmax  = 10;
  const short lockmin  = 1;
  REAL4TimeSeries tser;
  INT2TimeSeries  dmro;
  INT2TimeSeries  lock;
  FrChanIn  dmroin = { CHANNEL, ADCDataChannel };
  FrChanIn  lockin = { "IFO_Lock", ADCDataChannel };
  FrStream *stream = NULL;
  FrOutPar  outpar = { "C1:" CHANNEL, "REDUCED", ADCDataChannel, 6, 0, 0 };
  char *dirname = getenv( "LAL_FRAME_PATH" );
  int  locklost = 1;


  /*
   *
   * Open the frame stream.
   *
   */

  LALFrOpen( &status, &stream, dirname, "C1-*.F" );
  TESTSTATUS( &status );


  /*
   *
   * Set up time series to read DMRO and Lock channels.
   *
   */

  dmro.data = NULL;
  lock.data = NULL;
  LALFrGetINT2TimeSeries( &status, &dmro, &dmroin, stream );
  TESTSTATUS( &status );
  LALFrGetINT2TimeSeries( &status, &lock, &lockin, stream );
  TESTSTATUS( &status );
  LALI2CreateVector( &status, &dmro.data, duration / dmro.deltaT );
  TESTSTATUS( &status );
  LALI2CreateVector( &status, &lock.data, duration / lock.deltaT );
  TESTSTATUS( &status );
  memcpy( &tser, &dmro, sizeof( tser ) );
  tser.data = NULL;
  LALSCreateVector( &status, &tser.data, duration / tser.deltaT );
  TESTSTATUS( &status );


  /*
   *
   * Scan frame data writing locked data (after 3 minutes into lock).
   *
   */

  while ( 1 )
  {
    FrPos frpos;
    UINT4 i;
    INT8  tacc = 0.1 * 1e9 / 16384;
    INT8  texp;
    INT8  tact;
    int   locked;


    /* find a segment that is locked */
    do
    {
      locked = 1;

      /* save the position of the start of this segment */
      LALFrGetPos( &status, &frpos, stream );
      TESTSTATUS( &status );

      /* read the lock channel */
      LALFrGetINT2TimeSeries( &status, &lock, &lockin, stream );
      if ( status.statusCode == FRAMESTREAMH_EDONE )
      {
        goto end;
      }
      TESTSTATUS( &status );

      /* see if there were any periods out of lock */
      for ( i = 0; i < lock.data->length; ++i )
      {
        if ( lock.data->data[i] < lockmin || lock.data->data[i] > lockmax )
        {
          locked = 0;
          locklost = 1;
        }
      }
      if ( locklost )
      {
        printf( "\nsearching for lock..." );
      }
    }
    while ( ! locked );


    /* if lock has been lost, skip ahead 3 minutes and re-check lock */
    if ( locklost )
    { /* skip 3 minutes into locked segment */
      LIGOTimeGPS epoch;
      puts( " acquired: skipping three minutes" );
      LALFrTell( &status, &epoch, stream );
      TESTSTATUS( &status );
      epoch.gpsSeconds += 3 * 60;
      LALFrSeek( &status, &epoch, stream );
      TESTSTATUS( &status );
      locklost = 0;
      continue;
    }


    /* reset stream position to start of locked segment */
    LALFrSetPos( &status, &frpos, stream );
    TESTSTATUS( &status );


    /* compute expected stop time for this segment (assuming no gaps) */
    texp  = (INT8)1000000000 * (INT8)dmro.epoch.gpsSeconds;
    texp += (INT8)dmro.epoch.gpsNanoSeconds;
    texp += (INT8)( 1e9 * dmro.data->length * dmro.deltaT );


    /* read the DMRO channel */
    LALFrGetINT2TimeSeries( &status, &dmro, &dmroin, stream );
    if ( status.statusCode == FRAMESTREAMH_EDONE )
    {
      goto end;
    }
    TESTSTATUS( &status );


    /* compute the actual stop time for this segment */
    tact  = (INT8)1000000000 * (INT8)dmro.epoch.gpsSeconds;
    tact += (INT8)dmro.epoch.gpsNanoSeconds;


    /* warn if a gap in the data is found */
    if ( abs( texp - tact ) > tacc && ! locklost )
      puts( "Gap in frame data!" );


    /* copy data to tser and output it */
    printf( "%s-%s-%d-%d.gwf\n", outpar.source, outpar.description,
        tser.epoch.gpsSeconds,
        (int)ceil( 1e-9 * tser.epoch.gpsNanoSeconds
          + tser.data->length * tser.deltaT ) );
    tser.epoch = dmro.epoch;
    for ( i = 0; i < tser.data->length; ++i )
    {
      tser.data->data[i] = dmro.data->data[i];
    }
    LALFrWriteREAL4TimeSeries( &status, &tser, &outpar );
    TESTSTATUS( &status );
  }


  /*
   *
   * Cleanup and exit.
   *
   */

end:

  LALFrClose( &status, &stream );
  TESTSTATUS( &status );

  LALI2DestroyVector( &status, &dmro.data );
  TESTSTATUS( &status );

  LALI2DestroyVector( &status, &lock.data );
  TESTSTATUS( &status );

  LALSDestroyVector( &status, &tser.data );
  TESTSTATUS( &status );

  LALCheckMemoryLeaks();
  return 0;
}
