int isnan(double value);
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <FrameL.h>
#include "series.h"

int lalDebugLevel = 0;

char *get_next_line( char *line, size_t size, FILE *fp )
{
  char *s;
  do
    s = fgets( line, size, fp );
  while ( ( line[0] == '#' || line[0] == '%' ) && s );
  return s;
}

int read_freq_series( struct series *ser, const char *fname )
{
  float a, b;
  char line[256];
  int n;
  FILE *fp;
  fp = fopen( fname, "r" );
  n = 0;

  ser->dom = Trans;
  ser->type = FR_VECT_8C;
  
  /* count lines */
  while ( get_next_line( line, sizeof( line ), fp ) )
      ++n;
  rewind( fp );

  ser->size = n;
  ser->data = calloc( 2 * n, sizeof( *ser->data ) );

  get_next_line( line, sizeof( line ), fp );
  sscanf( line, "%*le %e %e\n", &a, &b );
  if ( isnan( a ) || isnan( b ) )
  {
    ser->data[0] = 1;
    ser->data[1] = 0;
  }
  else
  {
    ser->data[0] = a * cos( b );
    ser->data[1] = a * sin( b );
  }

  get_next_line( line, sizeof( line ), fp );
  sscanf( line, "%le %e %e\n", &ser->step, &a, &b );
  if ( isnan( a ) || isnan( b ) )
  {
    ser->data[2] = 1;
    ser->data[3] = 0;
  }
  else
  {
    ser->data[2] = a * cos( b );
    ser->data[3] = a * sin( b );
  }

  n = 2;
  while ( get_next_line( line, sizeof( line ), fp ) )
  {
    sscanf( line, "%*e %e %e\n", &a, &b );
    if ( isnan( a ) || isnan( b ) )
    {
      ser->data[2*n]   = 1;
      ser->data[2*n+1] = 1;
    }
    else
    {
      ser->data[2*n]   = a * cos( b );
      ser->data[2*n+1] = a * sin( b );
    }
    ++n;
  }

  fclose( fp );
  return n;
}

#define CALURL "http://blue.ligo-wa.caltech.edu/engrun/Calib_Home/html/cal_home.html"
#define C_CHANNEL "CAL-CAV_GAIN"
#define R_CHANNEL "CAL-RESPONSE"
#define USAGE( s ) do { \
  fprintf( stderr, "Usage: %s -run run -time gpssec -ifo 'H1'|'H2'|'L1'", s ); \
  fprintf( stderr, " -ver VXX R-file C-file \n\t [ [ -ifo" );\
  fprintf( stderr, " 'H1'|'H2'|'L1' -ver VXX R-file2 C-file2 ] ...]\n" ); \
  fprintf( stderr, "Calibration files found at URL:\n" CALURL "\n" ); \
  exit( 1 ); } while ( 0 )

int main( int argc, char *argv[] )
{
  static char   site;
  static size_t size;
  static double step;
  struct FrFile *frfile = NULL;
  struct FrameH *frame  = NULL;
  const char *run = NULL;
  const char *ifo = NULL;
  const char *ver = NULL;
  int sec = 0;
  int arg;

  /* parse arguments */
  if ( argc == 1 )
  {
    USAGE( argv[0] );
  }
  for ( arg = 1; arg < argc; ++arg )
  {
    char Rname[64];
    char Rilwd[64];
    char Cname[64];
    char Cilwd[64];
    if ( strstr( argv[arg], "-run" ) )
    {
      if ( run )
      {
        fprintf( stderr, "Error: run \"%s\" already specified\n", run );
        USAGE( argv[0] );
      }
      run = argv[++arg];
    }
    else if ( strstr( argv[arg], "-time" ) )
    {
      if ( sec )
      {
        fprintf( stderr, "Error: GPS time %d already specified\n", sec );
        USAGE( argv[0] );
      }
      sec = atoi( argv[++arg] );
    }
    else if ( strstr( argv[arg], "-ifo" ) )
    {
      if ( ifo )
      {
        fprintf( stderr, "Error: ifo \"%s\" already specified\n", ifo );
        USAGE( argv[0] );
      }
      if ( ! run )
      {
        fprintf( stderr, "Error: run not specified\n" );
        USAGE( argv[0] );
      }
      ifo = argv[++arg];
      if ( step && step != *ifo )
      {
        fprintf( stderr, "Error: ifos must all be at same site\n" );
        return 1;
      }
      site = *ifo;
      sprintf( Rname, "%s\\:" R_CHANNEL, ifo );
      sprintf( Rilwd, "%s-%s-" R_CHANNEL ".ilwd", run, ifo );
      sprintf( Cname, "%s\\:" C_CHANNEL, ifo );
      sprintf( Cilwd, "%s-%s-" C_CHANNEL ".ilwd", run, ifo );
    }
    else if ( strstr( argv[arg], "-ver" ) )
    {
      if ( ! run || ! ifo )
      {
        fprintf( stderr, "Error: run and/or ifo not specified\n" );
        USAGE( argv[0] );
      }
      ver = argv[++arg];
    }
    else if ( strstr( argv[arg], "-h" ) )
    {
      USAGE( argv[0] );
    }
    else
    {
      struct series R;
      struct series C;
      int code;
      if ( ! ifo || ! run || ! ver )
      {
        fprintf( stderr, "Error: ifo, run or version not specified\n" );
        USAGE( argv[0] );
      }
      if ( ! sec )
      {
        fprintf( stderr, "Error: GPS time not specified\n" );
        USAGE( argv[0] );
      }
      /* get R and C data and metadata */
      R.name = Rname;
      R.unit = "strain/ct";
      C.name = Cname;
      C.unit = "ct/strain";
      R.tbeg.gpsSeconds = C.tbeg.gpsSeconds = sec;
      R.tbeg.gpsNanoSeconds = C.tbeg.gpsNanoSeconds = 0;
      code = read_freq_series( &R, argv[arg] );
      if ( ! code )
      {
        fprintf( stderr, "Error: could not read file %s\n", argv[arg] );
        return 1;
      }
      if ( ! size )
      {
        size = R.size;
        step = R.step;
      }
      if ( size != R.size || step != R.step )
      {
        fprintf( stderr, "Error: data domain mismatch\nAll data files must" );
        fprintf( stderr, " have same length and frequency interval\n" );
        return 1;
      }
      code = read_freq_series( &C, argv[++arg] );
      if ( ! code )
      {
        fprintf( stderr, "Error: could not read file %s\n", argv[arg] );
        return 1;
      }
      if ( size != C.size || step != C.step )
      {
        fprintf( stderr, "Error: data domain mismatch\nAll data files must" );
        fprintf( stderr, " have same length and frequency interval\n" );
        return 1;
      }
      epoch_add( &R.tend, &R.tbeg, (R.size - 1)/(R.step * R.size) );
      epoch_add( &C.tend, &C.tbeg, (C.size - 1)/(C.step * C.size) );
#if 0
      code = write_ilwd( Rilwd, &R );
      if ( code )
      {
        fprintf( stderr, "Error: could not write file %s\n", Rilwd );
        return 1;
      }
      code = write_ilwd( Cilwd, &C );
      if ( code )
      {
        fprintf( stderr, "Error: could not write file %s\n", Cilwd );
        return 1;
      }
#endif
      if ( ! frfile )
      {
        char fname[256];
        int dt = (int)ceil( 1e-9 * R.tbeg.gpsNanoSeconds + 1.0 / R.step );
        sprintf( fname, "%c-CAL_REF_%s_%s-%d-%d.gwf", *ifo, ver, ifo, 
            R.tbeg.gpsSeconds, dt );
        frfile = FrFileONew( fname, 0 );
      }
      /* don't mangle the channel names for frames */
      R.name = Rname;
      C.name = Cname;
      sprintf( Rname, "%s:" R_CHANNEL, ifo );
      sprintf( Cname, "%s:" C_CHANNEL, ifo );
      frame = fr_add_proc_data( frame, &R );
      frame = fr_add_proc_data( frame, &C );
      FrameWrite( frame, frfile );
    }
  }

  FrFileOEnd( frfile );
  return 0;
}
