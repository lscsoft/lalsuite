#include <lal/LALStdlib.h>
#include <lal/FrameCache.h>

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 1; } else ((void)0)

int lalDebugLevel = 1;

int main( void )
{
  static LALStatus     status;
  static FrCache      *input;
  static FrCache      *cache;
  static FrCacheSieve  sieve;

  LALFrCacheImport( &status, &input, "catalog.test" );
  TESTSTATUS( &status );

  sieve.srcRegEx = "[F]";
  LALFrCacheSieve( &status, &cache, input, &sieve );
  TESTSTATUS( &status );

  LALFrCacheExport( &status, cache, "catalog.out" );
  TESTSTATUS( &status );
  LALDestroyFrCache( &status, &cache );
  TESTSTATUS( &status );
  LALDestroyFrCache( &status, &input );
  TESTSTATUS( &status );

  LALCheckMemoryLeaks();
  return 0;
}
