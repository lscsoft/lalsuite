int isnan(double value);
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <FrameL.h>
#include <series.h>

int lalDebugLevel = 0;

char *get_next_line( char *line, size_t size, FILE *fp )
{
  char *s;
  do
    s = fgets( line, size, fp );
  while ( ( line[0] == '#' || line[0] == '%' ) && s );
  return s;
}

int read_time_series( struct series *aser, struct series *abser,
    const char *fname )
{
  char line[256];
  int n;
  int t0;
  int t1;
  int dt;
  FILE *fp;

  aser->dom = Time;
  aser->type = FR_VECT_8C;
  abser->dom = Time;
  abser->type = FR_VECT_8C;

  fp = fopen( fname, "r" );

  get_next_line( line, sizeof( line ), fp );
  sscanf( line, "%d", &t0 );
  get_next_line( line, sizeof( line ), fp );
  sscanf( line, "%d", &dt );
  dt -= t0;

  aser->tbeg.gpsSeconds  = t0;
  aser->tbeg.gpsNanoSeconds  = 0;
  abser->tbeg.gpsSeconds = t0;
  abser->tbeg.gpsNanoSeconds = 0;
  aser->step  = dt;
  abser->step = dt;

  /* scan to end to get final time */
  while ( get_next_line( line, sizeof( line ), fp ) )
    ;
  rewind( fp );
  sscanf( line, "%d", &t1 );
  aser->tend.gpsSeconds  = t1 - dt;
  aser->tend.gpsNanoSeconds  = 0;
  abser->tend.gpsSeconds = t1 - dt;
  abser->tend.gpsNanoSeconds = 0;
  n = 1 + ( t1 - t0 ) / dt;
  aser->size  = n;
  abser->size = n;

  aser->data  = calloc( 2 * n, sizeof( *aser->data  ) );
  abser->data = calloc( 2 * n, sizeof( *abser->data ) );

  while ( get_next_line( line, sizeof( line ), fp ) )
  {
    float a;
    float ab;
    int t;
    sscanf( line, "%d %f %f", &t, &ab, &a );
    if ( ! isnan( a ) && ! isnan( ab ) )
    {
      int i;
      i = ( t - t0 ) / dt;
      aser->data[2*i]  = a;
      abser->data[2*i] = ab;
    }
  }

  fclose( fp );
  return n;
}

#define CALURL "http://blue.ligo-wa.caltech.edu/engrun/Calib_Home/html/cal_home.html"

#define USAGE( s ) do { \
  fprintf( stderr, "Usage: %s -run run -ifo 'H1'|'H2'|'L1' file", s );\
  fprintf( stderr, " [ [ -ifo 'H1'|'H2'|'L1' file2 ] ...]\n" ); \
  fprintf( stderr, "Calibration files found at URL:\n" CALURL "\n" ); \
  exit( 1 ); } while ( 0 )

#define A_CHANNEL "CAL-CAV_FAC"
#define AB_CHANNEL "CAL-OLOOP_FAC"

int main( int argc, char *argv[] )
{
  static LIGOTimeGPS  tbeg;
  static size_t size;
  static double step;
  static char   site;
  struct FrFile *frfile = NULL;
  struct FrameH *frame  = NULL;
  const char *run = NULL;
  const char *ifo = NULL;
  int arg;

  /* parse arguments */
  if ( argc == 1 )
  {
    USAGE( argv[0] );
  }
  for ( arg = 1; arg < argc; ++arg )
  {
    char aname[64];
    char ailwd[64];
    char abname[64];
    char abilwd[64];
    if ( strstr( argv[arg], "-run" ) )
    {
      if ( run )
      {
        fprintf( stderr, "Error: run \"%s\" already specified\n", run );
        USAGE( argv[0] );
      }
      run = argv[++arg];
    }
    else if ( strstr( argv[arg], "-ifo" ) )
    {
      if ( ! run )
      {
        fprintf( stderr, "Error: run not specified\n" );
        USAGE( argv[0] );
      }
      ifo = argv[++arg];
      if ( site && site != *ifo )
      {
        fprintf( stderr, "Error: ifos must all be at same site\n" );
        return 1;
      }
      site = *ifo;
      sprintf( aname, "%s\\:" A_CHANNEL, ifo );
      sprintf( ailwd, "%s-%s-" A_CHANNEL ".ilwd", run, ifo );
      sprintf( abname, "%s\\:" AB_CHANNEL, ifo );
      sprintf( abilwd, "%s-%s-" AB_CHANNEL ".ilwd", run, ifo );
    }
    else if ( strstr( argv[arg], "-h" ) )
    {
      USAGE( argv[0] );
    }
    else
    {
      struct series a;
      struct series ab;
      int code;
      if ( ! ifo )
      {
        fprintf( stderr, "Error: ifo not specified\n" );
        USAGE( argv[0] );
      }
      /* get a and ab data and metadata */
      a.name  = aname;
      a.unit  = "none";
      ab.name = abname;
      ab.unit = "none";
      code = read_time_series( &a, &ab, argv[arg] );
      if ( ! code )
      {
        fprintf( stderr, "Error: could not read file %s\n", argv[arg] );
        return 1;
      }
      if ( size )
      {
        if ( size != a.size || step != a.step 
            || tbeg.gpsSeconds != a.tbeg.gpsSeconds || 
            tbeg.gpsNanoSeconds != a.tbeg.gpsNanoSeconds )
        {
          fprintf( stderr, "Error: data domain mismatch\n" );
          fprintf( stderr, "Error: all series must have same time range\n" );
          return 1;
        }
      }
      else
      {
        size = a.size;
        step = a.step;
        tbeg = a.tbeg;
      }
      code = write_ilwd( ailwd, &a );
      if ( code )
      {
        fprintf( stderr, "Error: could not write file %s\n", ailwd );
        return 1;
      }
      code = write_ilwd( abilwd, &ab );
      if ( code )
      {
        fprintf( stderr, "Error: could not write file %s\n", abilwd );
        return 1;
      }
      /* correct the end times */
      epoch_add( &a.tend, &a.tbeg, a.step * a.size );
      epoch_add( &ab.tend, &ab.tbeg, ab.step * ab.size );
      if ( ! frfile )
      {
        char fname[256];
        int dt = (int)ceil( epoch_diff( &a.tend, &a.tbeg ) );
        sprintf( fname, "%c-CAL_FAC-%d-%d.gwf", *ifo, a.tbeg.gpsSeconds, dt );
        frfile = FrFileONew( fname, 0 );
      }
      /* don't mangle the channel names for frames */
      a.name  = aname;
      ab.name = abname;
      sprintf( aname, "%s:" A_CHANNEL, ifo );
      sprintf( abname, "%s:" AB_CHANNEL, ifo );
      frame = fr_add_proc_data( frame, &a );
      frame = fr_add_proc_data( frame, &ab );
    }
  }

  FrameWrite( frame, frfile );
  FrFileOEnd( frfile );
  return 0;
}
