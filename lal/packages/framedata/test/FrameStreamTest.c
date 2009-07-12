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

/**** <lalVerbatim file="FrameStreamTestCV">
 * Author: Jolien D. E. Creighton
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * \subsection{Program \texttt{FrameStreamTest.c}}
 *
 * Tests the low-level frame stream routines.
 *
 * \subsubsection*{Usage}
 *
 * \begin{verbatim}
 * FrameStreamTest
 * \end{verbatim}
 *
 * \subsubsection*{Description}
 *
 * This program reads the channels \verb+H1:LSC-AS_Q+ from all the fake frames
 * \verb+F-TEST-*.gwf+ in the directory set in the environment
 * \verb+LAL_FRAME_PATH+ * (or the current directory if this environment is not
 * set) and prints them to files.
 *
 **** </lalLaTeX> */

#include <FrameL.h>

#include <stdio.h>
#include <unistd.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/PrintFTSeries.h>
#include <lal/FrameStream.h>

#include <lal/LALRCSID.h>
NRCSID (FRAMESTREAMTESTC,"$Id$");

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { \
    fprintf( stderr, "Failure on line %d\n", __LINE__ ); \
    REPORTSTATUS(pstat); return 1; \
  } else ((void)0)

#ifndef CHANNEL
#define CHANNEL "H1:LSC-AS_Q"
#endif

INT4 lalDebugLevel = LALWARNING | LALINFO;

int main( void )
{
  static LALStatus status;
  const UINT4 npts = 200001;
  FrChanIn  chanin = { CHANNEL, ADCDataChannel };
  FrStream *stream = NULL;
  FrPos     frpos;
  static INT4TimeSeries chan; /* must zero the f0 field */
  LIGOTimeGPS epoch;
  CHAR *dirname = getenv( "LAL_FRAME_PATH" );
  INT4 file = 0;

  /* test files are version 4 frames */
  /* ignore this test if using earlier version of frame library */
  if ( FRAMELIB_VERSION < 4 )
    return 77;

  chan.data = NULL;
  LALI4CreateVector( &status, &chan.data, npts );
  TESTSTATUS( &status );

  LALFrOpen( &status, &stream, dirname, "F-TEST-*.gwf" );
  TESTSTATUS( &status );

  if ( XLALFrSetMode( stream, LAL_FR_VERBOSE_MODE | LAL_FR_CHECKSUM_MODE ) )
    return 1;

  /* seek to some initial time */
  epoch.gpsSeconds     = 600000071;
  epoch.gpsNanoSeconds = 123456789;
  LALFrSeek( &status, &epoch, stream );
  TESTSTATUS( &status );

  /* save this position */
  LALFrGetPos( &status, &frpos, stream );
  TESTSTATUS( &status );

  while ( 1 )
  {
    LALTYPECODE typecode;
    CHAR fname[256];
    LALFrGetTimeSeriesType( &status, &typecode, &chanin, stream );
    if ( status.statusCode == FRAMESTREAMH_EDONE )
    {
      break;
    }
    TESTSTATUS( &status );
    if ( typecode != LAL_I4_TYPE_CODE )
    {
      fprintf( stderr, "Wrong data type!\n" );
      return 1;
    }
    LALFrGetINT4TimeSeries( &status, &chan, &chanin, stream );
    if ( status.statusCode == FRAMESTREAMH_EDONE )
    {
      break;
    }
    TESTSTATUS( &status );
    sprintf( fname, CHANNEL ".%03d", file++ );
    LALI4PrintTimeSeries( &chan, fname );
  }

  /* go back to saved time */
  LALFrSetPos( &status, &frpos, stream );
  TESTSTATUS( &status );

  LALFrGetINT4TimeSeries( &status, &chan, &chanin, stream );
  TESTSTATUS( &status );

  LALI4PrintTimeSeries( &chan, CHANNEL ".999" );

  LALFrClose( &status, &stream );
  TESTSTATUS( &status );

  LALI4DestroyVector( &status, &chan.data );
  TESTSTATUS( &status );

  LALCheckMemoryLeaks();
  return 0;
}
