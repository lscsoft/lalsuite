/**** <lalVerbatim file="RingTestCV">
 * Author: Jolien Creighton
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * \subsection{Program \texttt{RingTest.c}}
 * 
 * Calls the various routines in \verb+Ring.h+.
 *
 * \subsubsection*{Usage}
 *
 * \begin{verbatim}
 * RingTest
 * \end{verbatim}
 *
 * \subsubsection*{Description}
 *
 * The program creates a ring template, which is output to the file
 * \verb+ring.out+, and a black hole ring waveform, which is output to the
 * file \verb+bhring.out+.  The program then generates a template bank,
 * which is written to the file \verb+bank.out+.
 * 
 * \vfill{\footnotesize\input{RingTestCV}}
 **** </lalLaTeX> */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/Ring.h>

#define TestStatus( ps ) \
  if ( (ps)->statusCode ) { \
    fprintf( stderr, "Failed LAL routine near line %d\n", __LINE__ ); \
    exit( 1 ); \
  } else ((void)0)

int lalDebugLevel = LALMSGLVL3;

int main( void )
{
  const UINT4 npts = 4096;
  const REAL4 srate = 2048;
  static LALStatus status;
  static REAL4TimeSeries ring;
  static RingTemplateInput tmplt;
  static BlackHoleRingInput bhring;
  static RingTemplateBank *bank;
  static RingTemplateBankInput bankin;

  UINT4 i;
  FILE *fp;


  LALCreateVector( &status, &ring.data, npts );
  TestStatus( &status );
  ring.deltaT = 1 / srate;

  tmplt.frequency = 10;
  tmplt.quality = 10;

  LALComputeRingTemplate( &status, &ring, &tmplt );
  TestStatus( &status );
  
  fp = fopen( "ring.out", "w" );
  for ( i = 0; i < ring.data->length; ++i )
    fprintf( fp, "%e\t%e\t%e\n", i * ring.deltaT, ring.data->data[i],
        sqrt( 2 * LAL_PI )
        * exp( - LAL_PI * tmplt.frequency * ring.deltaT * i / tmplt.quality )
        * cos( 2 * LAL_PI * tmplt.frequency * ring.deltaT * i ) );
  fclose( fp );

  bhring.solarMasses = 100;
  bhring.dimensionlessSpin = 0.98;
  bhring.percentMassLoss = 1;
  bhring.distanceMpc = 1;

  LALComputeBlackHoleRing( &status, &ring, &bhring );
  TestStatus( &status );

  fp = fopen( "bhring.out", "w" );
  for ( i = 0; i < ring.data->length; ++i )
    fprintf( fp, "%e\t%e\n", i * ring.deltaT, ring.data->data[i] );
  fclose( fp );

  bankin.minFrequency = 100;
  bankin.maxFrequency = 500;
  bankin.minQuality = 2;
  bankin.maxQuality = 20;
  bankin.maxMismatch = 0.03;

  LALCreateRingTemplateBank( &status, &bank, &bankin );
  TestStatus( &status );
  fp = fopen( "bank.out", "w" );
  for ( i = 0; i < bank->numTmplt; ++i )
    fprintf( fp, "%e\t%e\n", bank->tmplt[i].frequency, bank->tmplt[i].quality );
  fclose( fp );

  LALDestroyRingTemplateBank( &status, &bank );
  TestStatus( &status );
  LALDestroyVector( &status, &ring.data );
  TestStatus( &status );

  LALCheckMemoryLeaks();
  return 0;
}
