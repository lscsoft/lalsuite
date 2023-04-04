/*
 * Copyright (C) 2023 Karl Wette
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA
 */

/*
 * Author: Karl Wette
 *
 * Read a REAL8 channel from frame files, and output new frame files with a
 * variety of channels types, to test whether MakeSFTs can handle them.
 */

#include <stdio.h>

#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/LALFrStream.h>

int main( int argc, char *argv[] )
{

  if ( argc != 4 ) {
    printf( "usage: %s <framecache> <channel> <outdir>\n", argv[0] );
    return 1;
  }

  /* Open frame stream from cache */
  LALCache *framecache = XLALCacheImport( argv[1] );
  XLAL_CHECK_MAIN( framecache != NULL, XLAL_EFUNC );
  LALFrStream *framestream = XLALFrStreamCacheOpen( framecache );
  XLAL_CHECK_MAIN( framestream != NULL, XLAL_EFUNC );

  /* Set up reading time series from REAL8 channel */
  REAL8TimeSeries *time_series_REAL8 = XLALCreateREAL8TimeSeries( argv[2], NULL, 0.0, 0.0, &lalDimensionlessUnit, 0 );
  XLAL_CHECK_MAIN( time_series_REAL8 != NULL, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALFrStreamGetREAL8TimeSeriesMetadata( time_series_REAL8, framestream ) == XLAL_SUCCESS, XLAL_EFUNC );
  time_series_REAL8 = XLALResizeREAL8TimeSeries( time_series_REAL8, 0, 1800.0 / time_series_REAL8->deltaT );
  XLAL_CHECK_MAIN( time_series_REAL8 != NULL, XLAL_EFUNC );

  /* Until end of frame stream */
  while ( !XLALFrStreamEnd( framestream ) ) {

    /* Read next chunk of time series */
    strcpy( time_series_REAL8->name, argv[2] );
    XLAL_CHECK_MAIN( XLALFrStreamGetREAL8TimeSeries( time_series_REAL8, framestream ) == XLAL_SUCCESS, XLAL_EFUNC );

    /* Create output frame */
    LALFrameH *outframe = XLALFrameNew( &time_series_REAL8->epoch, 1800.0, "MultiChFrame", 1, 0, 0 );
    XLAL_CHECK_MAIN( outframe != NULL, XLAL_EFUNC );

    /* Add REAL8 channels */
    snprintf( time_series_REAL8->name, sizeof( time_series_REAL8->name ), "%c%c:AdcREAL8", argv[2][0], argv[2][1] );
    XLAL_CHECK_MAIN( XLALFrameAddREAL8TimeSeriesAdcData( outframe, time_series_REAL8 ) == XLAL_SUCCESS, XLAL_EFUNC );
    snprintf( time_series_REAL8->name, sizeof( time_series_REAL8->name ), "%c%c:ProcREAL8", argv[2][0], argv[2][1] );
    XLAL_CHECK_MAIN( XLALFrameAddREAL8TimeSeriesProcData( outframe, time_series_REAL8 ) == XLAL_SUCCESS, XLAL_EFUNC );
    snprintf( time_series_REAL8->name, sizeof( time_series_REAL8->name ), "%c%c:SimREAL8", argv[2][0], argv[2][1] );
    XLAL_CHECK_MAIN( XLALFrameAddREAL8TimeSeriesSimData( outframe, time_series_REAL8 ) == XLAL_SUCCESS, XLAL_EFUNC );

    /* Make REAL4 copy of REAL8 time series */
    REAL4TimeSeries *time_series_REAL4 = XLALConvertREAL8TimeSeriesToREAL4 ( time_series_REAL8 );
    XLAL_CHECK_MAIN( time_series_REAL4 != NULL, XLAL_EFUNC );

    /* Add REAL4 channels */
    snprintf( time_series_REAL4->name, sizeof( time_series_REAL4->name ), "%c%c:AdcREAL4", argv[2][0], argv[2][1] );
    XLAL_CHECK_MAIN( XLALFrameAddREAL4TimeSeriesAdcData( outframe, time_series_REAL4 ) == XLAL_SUCCESS, XLAL_EFUNC );
    snprintf( time_series_REAL4->name, sizeof( time_series_REAL4->name ), "%c%c:ProcREAL4", argv[2][0], argv[2][1] );
    XLAL_CHECK_MAIN( XLALFrameAddREAL4TimeSeriesProcData( outframe, time_series_REAL4 ) == XLAL_SUCCESS, XLAL_EFUNC );
    snprintf( time_series_REAL4->name, sizeof( time_series_REAL4->name ), "%c%c:SimREAL4", argv[2][0], argv[2][1] );
    XLAL_CHECK_MAIN( XLALFrameAddREAL4TimeSeriesSimData( outframe, time_series_REAL4 ) == XLAL_SUCCESS, XLAL_EFUNC );

    /* Write output frame */
    char *fname = XLALStringAppendFmt( NULL, "%s/%c-%c%c_%s-%d-%d.gwf",
                                       argv[3],
                                       argv[2][0],
                                       argv[2][0], argv[2][1],
                                       "MultiChFrame",
                                       time_series_REAL8->epoch.gpsSeconds,
                                       1800 );
    XLAL_CHECK_MAIN( fname != NULL, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALFrameWrite( outframe, fname ) == XLAL_SUCCESS, XLAL_EFUNC );

    /* Free memory */
    XLALDestroyREAL4TimeSeries( time_series_REAL4 );
    XLALFree( fname );

  }

  /* Free memory */
  XLALDestroyREAL8TimeSeries( time_series_REAL8 );
  XLALFrStreamClose( framestream );
  XLALDestroyCache( framecache );

  LALCheckMemoryLeaks();

  return 0;

}
