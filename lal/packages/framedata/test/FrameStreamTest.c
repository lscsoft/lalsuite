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
 * \verb+F-*.F+ in the directory set in the environment \verb+LAL_FRAME_PATH+
 * (or the current directory if this environment is not set) and prints them
 * to files.
 *
 **** </lalLaTeX> */

#include <stdio.h>
#include <unistd.h>
#include <FrameL.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/PrintFTSeries.h>
#include <lal/FrameStream.h>

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 1; } else ((void)0)

#ifndef CHANNEL
#define CHANNEL "H1:LSC-AS_Q"
#endif

INT4 lalDebugLevel = LALMEMDBG;

int main( void )
{
  static LALStatus status;
  const unsigned int npts = 1048576;
  FrChanIn  chanin = { CHANNEL, ADCDataChannel };
  FrStream *stream = NULL;
  REAL4TimeSeries chan;
  char *dirname = getenv( "LAL_FRAME_PATH" );
  int file = 0;

  /* test files are version 4 frames */
  /* ignore this test if using earlier version of frame library */
  if ( FRAMELIB_VERSION < 4 )
    return 77;

  chan.data = NULL;
  LALCreateVector( &status, &chan.data, npts );
  TESTSTATUS( &status );

  LALFrOpen( &status, &stream, dirname, "F-*.F" );
  TESTSTATUS( &status );

  while ( 1 )
  {
    char fname[256];
    LALFrGetREAL4TimeSeries( &status, &chan, &chanin, stream );
    if ( status.statusCode == FRAMESTREAMH_EDONE )
    {
      break;
    }
    TESTSTATUS( &status );
    sprintf( fname, CHANNEL ".%03d", file++ );
    puts( fname );
    LALSPrintTimeSeries( &chan, fname );
  }

  LALFrClose( &status, &stream );
  TESTSTATUS( &status );

  LALDestroyVector( &status, &chan.data );
  TESTSTATUS( &status );

  LALCheckMemoryLeaks();
  return 0;
}
