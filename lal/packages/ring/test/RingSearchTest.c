/**** <lalVerbatim file="RingSearchTestCV">
 * Author: Jolien Creighton
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * \subsection{Program \texttt{RingSearchTest.c}}
 * 
 * Calls the various routines in \verb+RingSearch.h+.
 *
 * \subsubsection*{Usage}
 *
 * \begin{verbatim}
 * RingSearchTest
 * \end{verbatim}
 *
 * \subsubsection*{Description}
 *
 * This program filters some random Gaussian noise with an injected ringdown
 * over a very small template bank and prints the events found.
 *
 * \vfill{\footnotesize\input{RingSearchTestCV}}
 **** </lalLaTeX> */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/PrintFTSeries.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/Ring.h>
#include <lal/RingSearch.h>

#define TestStatus( ps ) \
  if ( (ps)->statusCode ) { \
    fprintf( stderr, "Failed LAL routine near line %d\n", __LINE__ ); \
    exit( 1 ); \
  } else ((void)0)

int lalDebugLevel = LALMSGLVL3;

int main( void )
{
  const float ampinj = 10;
  const float srate  = 1024;
  const int   nseg   = 3;

  const char *argv[] = { "filterparams", "-segsz", "65536", "-speclen", "4096",
      "-flow", "40", "-fmin", "150", "-fmax", "200", "-qmin", "2",
      "-qmax", "10", "-maxmm", "0.1", "-thresh", "6", "-scale", "2000", NULL };
  int argc = sizeof( argv ) / sizeof( *argv ) - 1;

  static LALStatus                status;
  static COMPLEX8FrequencySeries  response;
  static REAL4FrequencySeries     spectrum;
  static REAL4TimeSeries          channel;
  static REAL4TimeSeries          inject;
  static RingTemplateInput        injTmplt;
  static RandomParams            *randpar;

  static RingSearchInput          input;
  static RingSearchData           data;
  static RingEventList           *events;
  static RingSearchParams        *params;

  RAT4 minusOne = { -1, 0 };
  LALUnitPair unitPair;
  LALUnit unit;
  UINT4 i;
  UINT4 j;
  UINT4 k;

  
  /*
   *
   * InitSearch
   *
   */

  LALRingSearchInit( &status, &params, argv, argc );
  TestStatus( &status );

  /*
   *
   * ConditionData
   *
   */


  /* first make up the data */
  LALSCreateVector( &status, &channel.data, nseg * params->segmentSize );
  TestStatus( &status );
  LALSCreateVector( &status, &spectrum.data, params->segmentSize / 2 + 1 );
  TestStatus( &status );
  LALCCreateVector( &status, &response.data, params->segmentSize / 2 + 1 );
  TestStatus( &status );
  LALCreateRandomParams( &status, &randpar, 1 );
  TestStatus( &status );

  data.channel  = &channel;
  data.spectrum = &spectrum;
  /* data.spectrum = NULL; */
  data.response = &response;

  channel.deltaT      = 1 / srate;
  channel.sampleUnits = lalADCCountUnit;
  /*
  memset( channel.data->data, 0,
      channel.data->length * sizeof( *channel.data->data ) );
   */
  LALNormalDeviates( &status, channel.data, randpar );
  TestStatus( &status );

  injTmplt.frequency = 150;
  injTmplt.quality   = 2;
  inject.deltaT      = 1 / srate;
  inject.sampleUnits = lalADCCountUnit;
  LALSCreateVector( &status, &inject.data, params->segmentSize / 4 );
  TestStatus( &status );
  LALComputeRingTemplate( &status, &inject, &injTmplt );
  TestStatus( &status );
  for ( j = 0; j < inject.data->length; ++j )
  {
    channel.data->data[j + params->segmentSize / 2] +=
      ampinj * inject.data->data[j];
  }
  LALSDestroyVector( &status, &inject.data );
  TestStatus( &status );



  LALSPrintTimeSeries( &channel, "channel.dat" );

  spectrum.deltaF = srate / params->segmentSize;
  unitPair.unitOne = &lalADCCountUnit;
  unitPair.unitTwo = &lalADCCountUnit;
  LALUnitMultiply( &status, &spectrum.sampleUnits, &unitPair );
  TestStatus( &status );
  unitPair.unitOne = &spectrum.sampleUnits;
  unitPair.unitTwo = &lalSecondUnit;
  LALUnitMultiply( &status, &spectrum.sampleUnits, &unitPair );
  TestStatus( &status );
  for ( k = 0; k < spectrum.data->length; ++k )
  {
    spectrum.data->data[k] = 2 / srate;
  }

  response.deltaF = srate / params->segmentSize;
  LALUnitRaise( &status, &unit, &lalADCCountUnit, &minusOne );
  TestStatus( &status );
  unitPair.unitOne = &lalStrainUnit;
  unitPair.unitTwo = &unit;
  LALUnitMultiply( &status, &response.sampleUnits, &unitPair );
  TestStatus( &status );
  for ( k = 0; k < response.data->length; ++k )
  {
    response.data->data[k].re = 1;
    response.data->data[k].im = 0;
  }

  LALCPrintFrequencySeries( &response, "response.dat" );
  
  /* OK, actually do the data conditioning here */
  LALRingSearchConditionData( &status, params, &data );
  TestStatus( &status );

  LALSPrintFrequencySeries( params->invSpectrum, "invspec.dat" );
  for ( i = 0; i < params->numSegments; ++i )
  {
    char fname[64];
    sprintf( fname, "segment-%03d.dat", i );
    LALCPrintFrequencySeries( params->dataSegment + i, fname );
  }

  LALDestroyRandomParams( &status, &randpar );
  TestStatus( &status );
  LALCDestroyVector( &status, &response.data );
  TestStatus( &status );
  LALSDestroyVector( &status, &spectrum.data );
  TestStatus( &status );
  LALSDestroyVector( &status, &channel.data );
  TestStatus( &status );

  /*
   *
   * FilterData
   *
   */

  /* do just a couple of templates */
  input.startTemplate = 0;
  input.templatesToDo = params->templateBank->numTmplt;
  input.templatesToDo = input.templatesToDo > 2 ? 2 : input.templatesToDo;

  /* params->maximizeEvents = 0; */
  params->keepResults = 1;
  LALRingSearch( &status, &events, &input, params );
  TestStatus( &status );

  while ( events )
  {
    RingEventList *next = events->next;
    fprintf( stdout, "%9d.%09d %e %e %e %e\n",
        (int)( events->startTimeNS / 1000000000 ),
        (int)( events->startTimeNS % 1000000000 ),
        events->snr, events->amplitude, events->frequency, events->quality );
    LALFree( events );
    events = next;
  }

  for ( i = 0; i < params->numResults; ++i )
  {
    LALSPrintTimeSeries( params->result + i, params->result[i].name );
  }


  /*
   *
   * FinalizeSearch
   *
   */

  LALRingSearchFini( &status, &params );
  TestStatus( &status );

  LALCheckMemoryLeaks();
  return 0;
}
