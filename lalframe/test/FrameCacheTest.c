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

#include <lal/LALStdlib.h>
#include <lal/FrameCache.h>

#include <lal/LALRCSID.h>
NRCSID (FRAMECACHETESTC,"$Id$");

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
