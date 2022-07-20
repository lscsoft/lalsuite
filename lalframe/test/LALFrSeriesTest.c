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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/**
 * \author Jolien D. E. Creighton
 * \file
 *
 * ### Program FrameStreamTest.c ###
 *
 * Tests the low-level frame stream routines.
 *
 * ### Usage ###
 *
 * \code
 * FrameStreamTest
 * \endcode
 *
 * ### Description ###
 *
 * This program reads the channels <tt>H1:LSC-AS_Q</tt> from all the fake frames
 * <tt>F-TEST-*.gwf</tt> in the directory TEST_DATA_DIR, and prints them to files.
 *
 */

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/PrintFTSeries.h>
#include <lal/LALFrStream.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#ifndef CHANNEL
#define CHANNEL "H1:LSC-AS_Q"
#endif


int main( void )
{
  const UINT4 npts = 200001; /* read data in weirdly sized blocks */
  LALFrStream *stream;
  LALFrStreamPos frpos;
  INT4TimeSeries *chan;
  LIGOTimeGPS epoch;
  INT4 file;

  XLALSetErrorHandler(XLALAbortErrorHandler);

  stream = XLALFrStreamOpen( TEST_DATA_DIR, "F-TEST-*.gwf" );
  if ( !stream )
    return 1;

  if ( XLALFrStreamSetMode( stream, LAL_FR_STREAM_VERBOSE_MODE | LAL_FR_STREAM_CHECKSUM_MODE ) )
    return 1;

  /* seek to some initial time */
  epoch.gpsSeconds     = 600000071;
  epoch.gpsNanoSeconds = 123456789;
  if ( XLALFrStreamSeek( stream, &epoch ) )
    return 1;

  /* save this position */
  if ( XLALFrStreamGetpos( &frpos, stream ) )
    return 1;

  /* only channel name and npts are used.  other metadata will be filled in
   * later */
  chan = XLALCreateINT4TimeSeries( CHANNEL, &epoch, 0.0, 0.0, &lalDimensionlessUnit, npts );
  if ( !chan )
    return 1;

  for ( file = 0; file < 8; file++ )
  {
    CHAR fname[256];
    LALTYPECODE typecode = XLALFrStreamGetTimeSeriesType( CHANNEL, stream );
    if ( typecode != LAL_I4_TYPE_CODE )
    {
      if ( typecode < 0 )
        fprintf( stderr, "failure\n" );
      else
        fprintf( stderr, "Wrong data type!\n" );
      return 1;
    }
    if ( XLALFrStreamGetINT4TimeSeries( chan, stream ) )
      return 1;
    sprintf( fname, CHANNEL ".%03d", file );
    LALI4PrintTimeSeries( chan, fname );
  }

  /* go back to saved time */
  if ( XLALFrStreamSetpos( stream, &frpos ) )
    return 1;

  if ( XLALFrStreamGetINT4TimeSeries( chan, stream ) )
    return 1;

  LALI4PrintTimeSeries( chan, CHANNEL ".999" );

  XLALFrStreamClose( stream );

  XLALDestroyINT4TimeSeries( chan );

  return 0;
}
