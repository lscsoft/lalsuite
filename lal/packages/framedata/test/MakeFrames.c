/**** <lalVerbatim file="MakeFramesCV">
 * Author: Jolien D. E. Creighton
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * \subsection{Program \texttt{MakeFrames.c}}
 *
 * Make some frames with random Gaussian noise.
 * 
 * \subsubsection*{Usage}
 *
 * \begin{verbatim}
 * MakeFrames
 * \end{verbatim}
 *
 * \subsubsection*{Description}
 *
 * This program makes some frames with one ADC channel containing random
 * Gaussian noise.
 *
 **** </lalLaTeX> */

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/Units.h>
#include <lal/FrameStream.h>

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 1; } else ((void)0)

#ifndef CHANNEL
#define CHANNEL "H1:LSC-AS_Q"
#endif

INT4 lalDebugLevel = LALMSGLVL3;

int main( void )
{
  static LALStatus       status;
  static REAL4TimeSeries series;
  const UINT4   time = 60;          /* duration of frame file */
  const INT4    seed = 10;          /* random number seed     */
  const REAL4   rate = 16384;       /* sample rate (Hz)       */
  const UINT4   npts = time * rate; /* number of points       */
  FrOutPar      opar = { "F", ADCDataChannel, 6, 0, 0 };
  RandomParams *rpar = NULL;

  /* initialize */

  LALCreateRandomParams( &status, &rpar, seed );
  TESTSTATUS( &status );

  strncpy( series.name, CHANNEL, sizeof( series.name ) );
  series.epoch.gpsSeconds = 600000000;
  series.sampleUnits = lalADCCountUnit;
  series.deltaT = 1 / rate;
  LALCreateVector( &status, &series.data, npts );
  TESTSTATUS( &status );

  /* generate first frame file worth of data and write it */

  LALNormalDeviates( &status, series.data, rpar );
  TESTSTATUS( &status );

  LALFrWriteREAL4TimeSeries( &status, &series, &opar );
  TESTSTATUS( &status );

  /* generate second frame file worth of data and write it */

  series.epoch.gpsSeconds += time;
  LALNormalDeviates( &status, series.data, rpar );
  TESTSTATUS( &status );

  LALFrWriteREAL4TimeSeries( &status, &series, &opar );
  TESTSTATUS( &status );

  /* generate third frame file worth of data and write it */

  series.epoch.gpsSeconds += time;
  LALNormalDeviates( &status, series.data, rpar );
  TESTSTATUS( &status );

  LALFrWriteREAL4TimeSeries( &status, &series, &opar );
  TESTSTATUS( &status );

  /* cleanup */

  LALDestroyVector( &status, &series.data );
  TESTSTATUS( &status );

  LALDestroyRandomParams( &status, &rpar );
  TESTSTATUS( &status );

  LALCheckMemoryLeaks();
  return 0;
}
