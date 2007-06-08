/*
*  Copyright (C) 2007 Duncan Brown, Stephen Fairhurst
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

#include <config.h>
#include <stdio.h>
#include <FrameL.h>

int main( int argc, char *argv[] )
{
  char history[] = "Created by " PACKAGE "-" VERSION ".";
  struct FrFile *frfileout;

  if ( argc < 3 )
  {
    fprintf( stderr, "usage: %s outfile infiles\n", argv[0] );
    return 1;
  }

  if ( !strcmp(argv[1],"--output") )
  { 
    if ( argc < 4  )
    {
      fprintf( stderr, "alternative usage: %s --output outfile infiles\n", 
          argv[0] );
      return 1;
    }
    --argc;
    ++argv;
  }

  FrLibIni( NULL, stderr, 0 );

  if ( ! ( frfileout = FrFileONewH( *++argv, 0, history ) ) )
    return fputs( "could not open output file!\n", stderr ), 1;

  --argc;
  while ( --argc > 0 )
  {
    struct FrFile *frfilein;
    struct FrameH *frame;

    if ( ! ( frfilein = FrFileINew( *++argv ) ) )
      return fprintf( stderr, "input file %s not found!\n", *argv ), 1;
    while ( ( frame = FrameRead( frfilein ) ) )
    {
      FrameWrite( frame, frfileout );
      FrameFree( frame );
    }
    FrFileIEnd( frfilein );
  }

  FrFileOEnd( frfileout );

  return 0;
}
