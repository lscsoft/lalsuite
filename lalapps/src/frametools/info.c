#include <config.h>
#include <stdio.h>
#include <FrameL.h>

const char *typestr( int type );
const char *typestr( int type )
{
  switch ( type )
  {
    case FR_VECT_C:
      return "char";
    case FR_VECT_2S:
      return "int_2";
    case FR_VECT_4S:
      return "int_4";
    case FR_VECT_8S:
      return "int_8";
    case FR_VECT_1U:
      return "unsigned char";
    case FR_VECT_2U:
      return "unsigned int_2";
    case FR_VECT_4U:
      return "unsigned int_4";
    case FR_VECT_8U:
      return "unsigned int_8";
    case FR_VECT_4R:
      return "float";
    case FR_VECT_8R:
      return "double";
    case FR_VECT_C8:
      return "complex";
    case FR_VECT_C16:
      return "double complex";
    default:
      return "unknown";
  }
  return "unknown";
}


int main( int argc, char *argv[] )
{
  if ( argc == 1 )
  {
    fprintf( stderr, "usage: %s framefiles\n", argv[0] );
    return 1;
  }

  FrLibIni( NULL, stderr, 0 );  

  while ( --argc > 0 )
  {
    struct FrStatData *sdata;
    struct FrFile *frfile;
    struct FrameH *frame;
    char *fname = *++argv;
    if ( ! ( frfile = FrFileINew( fname ) ) )
      return fprintf( stderr, "file %s not found!\n", fname ), 1;
    fprintf( stdout, "\n>>> Info for frame file %s:\n", fname );
    sdata = frfile->sDataCur;
    while ( ( frame = FrameRead( frfile ) ) )
    {
      fprintf( stdout, "\n>> %s, run %u, frame %u:\n", frame->name, frame->run,
          frame->frame );
      fprintf( stdout, "GPS time (s)   = %u.%09u\n", frame->GTimeS,
          frame->GTimeN );
      /* fprintf( stdout, "local time (s) = %d\n", frame->localTime ); */
      fprintf( stdout, "leap seconds   = %u\n", frame->ULeapS );
      fprintf( stdout, "duration (s)   = %f\n", frame->dt );
      if ( frame->rawData )
      {
        struct FrAdcData *adc = frame->rawData->firstAdc;
        fprintf( stdout, "adc channels:\n" );
        while ( adc )
        {
          fprintf( stdout, "\t%s: (crate %u, channel %u)",
              adc->name, adc->channelGroup, adc->channelNumber );
          if ( adc->data )
            fprintf( stdout, ", %ld %s points @ %f Hz", (long)adc->data->nData,
                typestr( adc->data->type ), adc->sampleRate );
          fprintf( stdout, "\n" );
          adc = adc->next;
        }
      }
      if ( frame->procData )
      {
        struct FrProcData *proc = frame->procData;
        while ( proc )
        {
#if defined FR_VERS && FR_VERS < 5000
          fprintf( stdout, "\t%s: srate = %f Hz,", proc->name,
              proc->sampleRate );
#else
          fprintf( stdout, "\t%s:", proc->name );
#endif
          if ( proc->data )
          {
            int dim;
            fprintf( stdout, " %ld %s points [%s]", (long)proc->data->nData,
                typestr( proc->data->type ), proc->data->unitY );
            for ( dim = 0; dim < (int)proc->data->nDim; ++dim )
              fprintf( stdout, ", nx(%d) = %ld dx(%d) = %f %s", dim,
                  (long)proc->data->nx[dim], dim, proc->data->dx[dim],
                  proc->data->unitX[dim] );
          }
          fprintf( stdout, "\n" );
          proc = proc->next;
        }
      }
      FrameFree( frame );
    }
    FrFileIEnd( frfile );
  }

  return 0;
}
