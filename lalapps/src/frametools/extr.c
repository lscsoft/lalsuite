#include <config.h>
#include <stdio.h>
#include <FrameL.h>

int main( int argc, char *argv[] )
{
  char history[] = "Created by " PACKAGE "-" VERSION ".";
  struct FrFile *frfileout;
  struct FrFile *frfilein;
  struct FrameH *frameout;
  struct FrameH *framein;
  char  *fnameout;
  char  *fnamein;
  char **channels;
  int    nchannels;

  if ( argc < 4 )
  {
    fprintf( stderr, "usage: %s outfile infile channels\n", argv[0] );
    return 1;
  }

  fnameout  = argv[1];
  fnamein   = argv[2];
  channels  = &argv[3];
  nchannels = argc - 3;

  FrLibIni( NULL, stderr, 0 );

  if ( ! ( frfilein = FrFileINew( fnamein ) ) )
    return fputs( "input file not found!\n", stderr ), 1;
  if ( ! ( frfileout = FrFileONewH( fnameout, 0, history ) ) )
    return fputs( "could not open output file!\n", stderr ), 1;

  while ( ( framein = FrameRead( frfilein ) ) )
  {
    int chan;
    if ( ! ( frameout = FrameHNew( framein->name ) ) )
      return fputs( "allocation error!\n", stderr ), 1;
    frameout->run = framein->run;
    frameout->frame = framein->frame;
    frameout->dataQuality = framein->dataQuality;
    frameout->GTimeS = framein->GTimeS;
    frameout->GTimeN = framein->GTimeN;
    frameout->ULeapS = framein->ULeapS;
    /* frameout->localTime = framein->localTime; */
    frameout->dt = framein->dt;
    for ( chan = 0; chan < nchannels; ++chan )
    {
      struct FrAdcData *adc;
      if ( ! ( adc = FrAdcDataFind( framein, channels[chan] ) ) )
        return fprintf( stderr, "channel %s not found!\n", channels[chan] ), 1;
      if ( ! FrAdcDataCopy( adc, frameout ) )
        return fputs( "allocation error!\n", stderr ), 1;
    }
    FrameWrite( frameout, frfileout );    
    FrameFree( frameout );
    FrameFree( framein );
  }

  FrFileOEnd( frfileout );
  FrFileIEnd( frfilein );

  return 0;
}
