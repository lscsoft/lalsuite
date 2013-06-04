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
 * \heading{Program \ref ReadLHOData.c}
 *
 * Tests the low-level frame stream routines by reading LHO E5 frame data.
 *
 * \heading{Usage}
 *
 * \code
 * ReadLHOData
 * \endcode
 *
 * \heading{Description}
 *
 * This program extracts the channel <tt>H2:LSC-AS_Q</tt> from all the frames in
 * the directory set in the environment \c LAL_FRAME_PATH (or the current
 * directory if this environment is not set) and writes it to new frame files.
 *
*/

#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/PrintFTSeries.h>
#include <lal/LALFrStream.h>

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 1; } else ((void)0)

#ifndef CHANNEL
#define CHANNEL "H2:LSC-AS_Q"
#endif


int main( void )
{
  static LALStatus status;
  const float duration = 60.0; /* seconds of data to read at a time */
  REAL4TimeSeries chan;
  FrChanIn  chanin = { CHANNEL, ADCDataChannel };
  LALFrStream *stream = NULL;
  FrOutPar  outpar = { CHANNEL, "TEST", ADCDataChannel, 6, 0, 0 };
  char *dirname = getenv( "LAL_FRAME_PATH" );

  /* open the frame stream */
  LALFrOpen( &status, &stream, dirname, "H-*.gwf" );
  TESTSTATUS( &status );

  /* get channel info */
  chan.data = NULL;
  LALFrGetREAL4TimeSeries( &status, &chan, &chanin, stream );
  TESTSTATUS( &status );

  /* allocate data storage for channel */
  LALSCreateVector( &status, &chan.data, duration / chan.deltaT );
  TESTSTATUS( &status );

  /* loop until all data has been read */
  while ( 1 )
  {
    INT8 tacc = 0.1 * 1e9 / 16384;
    INT8 texp;
    INT8 tact;

    texp  = (INT8)1000000000 * (INT8)chan.epoch.gpsSeconds;
    texp += (INT8)chan.epoch.gpsNanoSeconds;
    texp += (INT8)( 1e9 * chan.data->length * chan.deltaT );

    LALFrGetREAL4TimeSeries( &status, &chan, &chanin, stream );
    if ( status.statusCode == FRAMESTREAMH_EDONE )
    {
      break;
    }
    TESTSTATUS( &status );

    tact  = (INT8)1000000000 * (INT8)chan.epoch.gpsSeconds;
    tact += (INT8)chan.epoch.gpsNanoSeconds;

    if ( abs( texp - tact ) > tacc )
      puts( "Gap in frame data!" );

    printf( "%s-%s-%d-%d.gwf\n", outpar.source, outpar.description,
        chan.epoch.gpsSeconds,
        (int)ceil( 1e-9 * chan.epoch.gpsNanoSeconds
          + chan.data->length * chan.deltaT ) );
    LALFrWriteREAL4TimeSeries( &status, &chan, &outpar );
    TESTSTATUS( &status );
  }

  LALFrClose( &status, &stream );
  TESTSTATUS( &status );

  LALSDestroyVector( &status, &chan.data );
  TESTSTATUS( &status );

  LALCheckMemoryLeaks();
  return 0;
}
