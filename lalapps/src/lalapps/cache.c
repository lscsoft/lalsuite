/*
*  Copyright (C) 2007 Jolien Creighton
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALCache.h>
#include <lal/FileIO.h>

int main( int argc, char *argv[] )
{
	LALFILE *outfile = NULL;
	LALCache *cache;
	int arg = 1;
	XLALSetErrorHandler( XLALExitErrorHandler );
	if ( argc > 1 && ! strcmp(argv[1],"-o") ) {
		outfile = XLALFileOpen( argv[2], "w" );
		arg += 2;
	}
	cache = XLALCacheGlob( NULL, argc == 1 ? NULL : argv[arg] );
	for ( ; arg < argc; ++arg ) {
		LALCache *tmp = cache;
		LALCache *add;
		add = XLALCacheGlob( NULL, argv[arg] );
		cache = XLALCacheMerge( tmp, add );
		XLALDestroyCache( add );
		XLALDestroyCache( tmp );
	}
	XLALCacheFileWrite( outfile ? outfile : LALSTDOUT, cache );
	XLALFileClose( outfile );
	XLALDestroyCache( cache );
	LALCheckMemoryLeaks();
	return 0;
}
