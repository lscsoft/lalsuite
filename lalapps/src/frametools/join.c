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
