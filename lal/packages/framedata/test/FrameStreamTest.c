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
 * This program reads the channels \verb+LSC-AS_Q+ and \verb+LSC-AS_I+ from
 * all the frames in the directory set in the environment \verb+LAL_FRAME_PATH+
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

static void LALSDestroyTimeSeries( LALStatus *status, REAL4TimeSeries **series )
{
  INITSTATUS( status, "LALSDestroyTimeSeries", "$Id$" );
  ATTATCHSTATUSPTR( status );
  ASSERT( series, status, 1, "Null pointer" );
  ASSERT( *series, status, 1, "Null pointer" );
  TRY( LALSDestroyVector( status->statusPtr, &(*series)->data ), status );
  LALFree( *series );
  *series = NULL;
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

INT4 lalDebugLevel = LALMSGLVL3;

int main( void )
{
  static LALStatus  status;
  FrameStream      *stream = NULL;
  FrameStreamPos    frpos;
  INT4              error;
  UINT4             count;
  CHAR             *dirname;

  /* test files are version 4 frames */
  /* ignore this test if using earlier version of frame library */
  if ( FRAMELIB_VERSION < 4 )
    return 77;

  dirname = getenv( "LAL_FRAME_PATH" );
  LALOpenFrameStream( &status, &stream, dirname, NULL );
  TESTSTATUS( &status );

  LALFrameStreamGetPos( &status, &frpos, stream );
  TESTSTATUS( &status );

  error = count = 0;
  while ( error >= 0 && count < 10 )
  {
    CHAR fname[64];
    REAL4TimeSeries *series = NULL;

    LALSFrameReadADCTimeSeries( &status, &series, "LSC-AS_Q", stream );
    TESTSTATUS( &status );

    sprintf( fname, "LSC-AS_Q-%d.out", count++ );
    puts( fname );
    LALSPrintTimeSeries( series, fname );

    LALSDestroyTimeSeries( &status, &series );
    TESTSTATUS( &status );

    LALNextFrame( &status, &error, stream );
    TESTSTATUS( &status );
    if ( error > 0 )
      puts( "gap in data" );
  }

  LALFrameStreamSetPos( &status, &frpos, stream );
  TESTSTATUS( &status );

  error = count = 0;
  while ( error >= 0 && count < 10 )
  {
    CHAR fname[64];
    REAL4TimeSeries *series = NULL;

    LALSFrameReadADCTimeSeries( &status, &series, "LSC-AS_I", stream );
    TESTSTATUS( &status );

    sprintf( fname, "LSC-AS_I-%d.out", count++ );
    puts( fname );
    LALSPrintTimeSeries( series, fname );

    LALSDestroyTimeSeries( &status, &series );
    TESTSTATUS( &status );

    LALNextFrame( &status, &error, stream );
    TESTSTATUS( &status );
    if ( error > 0 )
      puts( "gap in data" );
  }

  LALCloseFrameStream( &status, &stream );
  TESTSTATUS( &status );

  LALCheckMemoryLeaks();
  return 0;
}
