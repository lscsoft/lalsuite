#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Units.h>
#include <lal/Random.h>
#include <lal/PrintFTSeries.h>
#include <lal/RingSearch.h>

#define TestStatus( ps ) \
  if ( (ps)->statusCode ) { \
    fprintf( stderr, "Failed LAL routine near line %d\n", __LINE__ ); \
    exit( 1 ); \
  } else ((void)0)

int lalDebugLevel = LALMSGLVL3;

int main( void )
{
  static LALStatus status;
  REAL4FrequencySeries spec;
  REAL4TimeSeries data;
  AvgSpecParams spar;
  RandomParams *rpar = NULL;
  UINT4 segsz = 1024;
  UINT4 nseg = 20;

  data.epoch.gpsSeconds     = 700000000;
  data.epoch.gpsNanoSeconds = 0;
  data.sampleUnits = lalADCCountUnit;
  data.deltaT = 1.0 / 1024.0;
  data.data = NULL;
  LALCreateVector( &status, &data.data, nseg * segsz );
  TestStatus( &status );

  spec.data = NULL;
  LALCreateVector( &status, &spec.data, segsz / 2 + 1 );
  TestStatus( &status );

  LALCreateRandomParams( &status, &rpar, 1 );
  TestStatus( &status );
  LALNormalDeviates( &status, data.data, rpar );
  TestStatus( &status );
  LALDestroyRandomParams( &status, &rpar );
  TestStatus( &status );

  spar.segsize = segsz;
  spar.fwdplan = NULL;
  LALCreateForwardRealFFTPlan( &status, &spar.fwdplan, segsz, 0 );

  LALMedianSpectrum( &status, &spec, &data, &spar );
  TestStatus( &status );

  LALSPrintFrequencySeries( &spec, "medspec.dat" );

  LALDestroyRealFFTPlan( &status, &spar.fwdplan );
  TestStatus( &status );

  LALDestroyVector( &status, &spec.data );
  TestStatus( &status );

  LALDestroyVector( &status, &data.data );
  TestStatus( &status );

  LALCheckMemoryLeaks();

  return 0;
}
