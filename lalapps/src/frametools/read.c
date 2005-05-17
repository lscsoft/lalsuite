#include <config.h>
#include <stdio.h>
#include <FrameL.h>

static const char *typestr( int type )
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

static size_t typesize( int type )
{
  switch ( type )
  {
    case FR_VECT_C:
      return sizeof( char );
    case FR_VECT_2S:
      return sizeof( short );
    case FR_VECT_4S:
      return sizeof( int );
    case FR_VECT_8S:
      return sizeof( FRLONG );
    case FR_VECT_1U:
      return sizeof( unsigned char );
    case FR_VECT_2U:
      return sizeof( unsigned short );
    case FR_VECT_4U:
      return sizeof( unsigned int );
    case FR_VECT_8U:
      return sizeof( FRULONG );
    case FR_VECT_4R:
      return sizeof( float );
    case FR_VECT_8R:
      return sizeof( double );
    case FR_VECT_C8:
      return 2 * sizeof( float );
    case FR_VECT_C16:
      return 2 * sizeof( double );
    default:
      return 0;
  }
  return 0;
}

static int writeFrVectData( FILE *fp, double x0, double dx, void *data,
    size_t nobj, int type )
{
  size_t i;
  switch ( type )
  {
    case FR_VECT_C:
      for ( i = 0; i < nobj; ++i )
      {
        char *val = data;
        fprintf( fp, "%e\t%hd\n", x0 + i * dx, val[i] );
      }
      break;
    case FR_VECT_2S:
      for ( i = 0; i < nobj; ++i )
      {
        short *val = data;
        fprintf( fp, "%e\t%d\n", x0 + i * dx, val[i] );
      }
      break;
    case FR_VECT_4S:
      for ( i = 0; i < nobj; ++i )
      {
        int *val = data;
        fprintf( fp, "%e\t%d\n", x0 + i * dx, val[i] );
      }
      break;
    case FR_VECT_8S:
      for ( i = 0; i < nobj; ++i )
      {
        FRLONG *val = data;
        fprintf( fp, "%e\t%ld\n", x0 + i * dx, (long)val[i] );
      }
      break;
    case FR_VECT_1U:
      for ( i = 0; i < nobj; ++i )
      {
        unsigned char *val = data;
        fprintf( fp, "%e\t%c\n", x0 + i * dx, val[i] );
      }
      break;
    case FR_VECT_2U:
      for ( i = 0; i < nobj; ++i )
      {
        unsigned short *val = data;
        fprintf( fp, "%e\t%hu\n", x0 + i * dx, val[i] );
      }
      break;
    case FR_VECT_4U:
      for ( i = 0; i < nobj; ++i )
      {
        unsigned int *val = data;
        fprintf( fp, "%e\t%u\n", x0 + i * dx, val[i] );
      }
      break;
      return sizeof( unsigned int );
    case FR_VECT_8U:
      for ( i = 0; i < nobj; ++i )
      {
        FRULONG *val = data;
        fprintf( fp, "%e\t%lu\n", x0 + i * dx, (unsigned long)val[i] );
      }
      break;
    case FR_VECT_4R:
      for ( i = 0; i < nobj; ++i )
      {
        float *val = data;
        fprintf( fp, "%e\t%e\n", x0 + i * dx, val[i] );
      }
      break;
    case FR_VECT_8R:
      for ( i = 0; i < nobj; ++i )
      {
        double *val = data;
        fprintf( fp, "%e\t%e\n", x0 + i * dx, val[i] );
      }
      break;
    case FR_VECT_C8:
      for ( i = 0; i < nobj; ++i )
      {
        float *val = data;
        fprintf( fp, "%e\t%e\t%e\n", x0 + i * dx, val[2*i], val[2*i+1] );
      }
      break;
    case FR_VECT_C16:
      for ( i = 0; i < nobj; ++i )
      {
        double *val = data;
        fprintf( fp, "%e\t%e\t%e\n", x0 + i * dx, val[2*i], val[2*i+1] );
      }
      break;
    default:
      return -1;
  }
  return nobj;
}

static int writeFrVect( FILE *fp, struct FrVect *v )
{
  size_t size;
  unsigned int dim;
  char *ptr;
  fprintf( fp, "# name  = %s\n", v->name ? v->name : "" );
  fprintf( fp, "# type  = %s\n", typestr( v->type ) );
  fprintf( fp, "# nData = %ld\n", (long)v->nData );
  fprintf( fp, "# nDim  = %d\n", v->nDim );
  fprintf( fp, "# dims  = %ld", (long)v->nx[0] );
  for ( dim = 1; dim < v->nDim; ++dim )
    fprintf( fp, ", %ld", (long)v->nx[dim] );
  fprintf( fp, "\n" );
  fprintf( fp, "# unitX = %s", v->unitX[0] ? v->unitX[0] : "" );
  for ( dim = 1; dim < v->nDim; ++dim )
    fprintf( fp, ", %s", v->unitX[dim] ? v->unitX[dim] : "" );
  fprintf( fp, "\n" );
  fprintf( fp, "# unitY = %s\n", v->unitY ? v->unitY : "" );
  size = typesize( v->type );
  ptr = v->data;
  writeFrVectData( fp, v->startX[0], v->dx[0], ptr, v->nx[0], v->type );
  ptr += size * v->nx[0];
  for ( dim = 1; dim < v->nDim; ++dim )
  {
    fprintf( fp, "&\n" );
    writeFrVectData( fp, v->startX[dim], v->dx[dim], ptr, v->nx[dim], v->type );
    ptr += size * v->nx[dim];
  }
  return 1;
}


int main( int argc, char *argv[] )
{
  struct FrFile *frfile;
  struct FrameH *frame;
  char *fname;
  char **channels;
  int nchannels;

  if ( argc < 3 )
  {
    fprintf( stderr, "usage: %s infilename channels\n", argv[0] );
    return 1;
  }

  fname     = argv[1];
  channels  = &argv[2];
  nchannels = argc - 2;

  FrLibIni( NULL, stderr, 0 );

  if ( ! ( frfile = FrFileINew( fname ) ) )
    return fputs( "file not found!\n", stderr ), 1;

  while ( ( frame = FrameRead( frfile ) ) )
  {
    int chan;
    for ( chan = 0; chan < nchannels; ++chan )
    {
      struct FrAdcData *adc;
      struct FrProcData *proc;
      char ofname[256];
      FILE *fp;
      if ( ( adc = FrAdcDataFind( frame, channels[chan] ) ) )
      {
        if ( frame->dt < 1 )
          sprintf( ofname, "%s-%03u.dat", channels[chan], frame->frame );
        else
          sprintf( ofname, "%s-%u.dat", channels[chan], frame->GTimeS );
        if ( ! ( fp = fopen( ofname, "w" ) ) )
          return fprintf( stderr, "could not open output file %s!", fname ), 1;
        fprintf( stderr, "writing to file %s\n", ofname );
        fprintf( fp, "## name             = %s\n", adc->name ? adc->name : "" );
        fprintf( fp, "## comment          = %s\n", adc->comment?adc->comment:"" );
        fprintf( fp, "## crate            = %u\n", adc->channelGroup );
        fprintf( fp, "## channel          = %u\n", adc->channelNumber );
        fprintf( fp, "## bias             = %f\n", adc->bias );
        fprintf( fp, "## slope            = %f\n", adc->slope );
        fprintf( fp, "## units            = %s\n", adc->units ? adc->units : "" );
        fprintf( fp, "## sample rate (Hz) = %f\n", adc->sampleRate );
        fprintf( fp, "## GPS time (s)     = %u.%09u\n", frame->GTimeS,
            frame->GTimeN );
        /* fprintf( fp, "## time offset (s)  = %u.%09u\n", adc->timeOffsetS,
           adc->timeOffsetN ); */
        fprintf( fp, "## freq shift (Hz)  = %f\n", adc->fShift );
        fprintf( fp, "## data valid       = %s\n",
            adc->dataValid ? "no" : "yes" );
        if ( adc->data )
        {
          writeFrVect( fp, adc->data );
        }
      }
      else if ( ( proc = FrProcDataFind( frame, channels[chan] ) ) )
      {
        if ( frame->dt < 1 )
          sprintf( ofname, "%s-%03u.dat", channels[chan], frame->frame );
        else
          sprintf( ofname, "%s-%u.dat", channels[chan], frame->GTimeS );
        if ( ! ( fp = fopen( ofname, "w" ) ) )
          return fprintf( stderr, "could not open output file %s!", fname ), 1;
        fprintf( stderr, "writing to file %s\n", ofname );
        fprintf( fp, "## name             = %s\n",
            proc->name ? proc->name : "" );
        fprintf( fp, "## comment          = %s\n",
            proc->comment ? proc->comment : "" );
        fprintf( fp, "## GPS time (s)     = %u.%09u\n", frame->GTimeS,
            frame->GTimeN );
        fprintf( fp, "## freq shift (Hz)  = %f\n", proc->fShift );
        if ( proc->data )
        {
          writeFrVect( fp, proc->data );
        }
      }
      else
      {
        return fprintf( stderr, "channel %s not found!\n", channels[chan] ), 1;
      }
    }
    FrameFree( frame );
  }

  FrFileIEnd( frfile );

  return 0;
}
