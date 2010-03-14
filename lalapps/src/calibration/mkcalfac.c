/*
*  Copyright (C) 2007 Alexander Dietz, Duncan Brown, Jolien Creighton
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

int isnan(double value);
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALFrameL.h>
#include <series.h>

RCSID( "$Id$" );

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "mkcalfac"

extern int lalDebugLevel;
static int sensemon_format;
static int skip_first_line;
int deltaT =0;
const char *run = NULL;
const char *channel = NULL;

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
  int t1 = 0;
  int dt;
  FILE *fp;
  float a,b,ab;
  int t, tmp_t1;

  aser->dom = Time;
  aser->type = FR_VECT_8C;
  abser->dom = Time;
  abser->type = FR_VECT_8C;

  fp = fopen( fname, "r" );

  if ( skip_first_line )
  {
    /* skip the first line as sensmon seems to bugger this up */
    get_next_line( line, sizeof( line ), fp );
  }

  get_next_line( line, sizeof( line ), fp );
  sscanf( line, "%d", &t0 );
  dt=deltaT;

  aser->tbeg.gpsSeconds  = t0;
  aser->tbeg.gpsNanoSeconds  = 0;
  abser->tbeg.gpsSeconds = t0;
  abser->tbeg.gpsNanoSeconds = 0;
  aser->step  = dt;
  abser->step = dt;

  /* scan to end to get final time */
  while ( get_next_line( line, sizeof( line ), fp ) )
  {

    sscanf( line, "%d", &tmp_t1 );
    if ( ( (tmp_t1 - t0) % dt ) == 0 )
    {
      t1 = tmp_t1;
    }

  }
  rewind( fp );
  aser->tend.gpsSeconds  = t1 - dt;
  aser->tend.gpsNanoSeconds  = 0;
  abser->tend.gpsSeconds = t1 - dt;
  abser->tend.gpsNanoSeconds = 0;
  n = 1 + ( t1 - t0 ) / dt;
  aser->size  = n;
  abser->size = n;

  aser->data  = calloc( 2 * n, sizeof( *aser->data  ) );
  abser->data = calloc( 2 * n, sizeof( *abser->data ) );

  if ( skip_first_line )
  {
    /* skip the first line as sensmon seems to bugger this up */
    get_next_line( line, sizeof( line ), fp );
  }

  while ( get_next_line( line, sizeof( line ), fp ) )
  {


    if ( sensemon_format )
    {
      /* sensemon format is: Time Range Line Alpha Beta */
      sscanf( line, "%d %*f %*f %f %f", &t, &a, &b );
      ab = a * b;
    }
    else
    {
      if ( !strcmp(run, "S3")) {
	/* gaby format is: Time alpha*beta alpha beta */
	sscanf( line, "%d %f %f", &t, &ab, &a );
      } else {
	sscanf( line, "%d %f %f", &t, &a, &ab );
      }

    }

    if ( (t - t0) % deltaT )
    {
      fprintf( stderr, "warning: skipping line\n\t%s\n", line );
    }
    else
    {
      if ( ! isnan( a ) && ! isnan( ab ) )
      {
        int i;
        i = ( t - t0 ) / dt;
        aser->data[2*i]  = a;
        abser->data[2*i] = ab;
      }
    }
  }

  fclose( fp );
  return n;
}

#define CALURL "http://blue.ligo-wa.caltech.edu/engrun/Calib_Home/html/cal_home.html"

void usage ( char *s )
{
  fprintf( stdout, "Usage: %s [options] [factorfile]\n", s );
  fprintf( stdout, "\nOptions:\n" );
  fprintf( stdout, "  --help                print this message\n");
  fprintf( stdout, "  --ifo IFO             set the frame ifo name to IFO (L1, H1, etc.)\n" );
  fprintf( stdout, "  --run RUN             set the frame run name to RUN (E11, S2, etc.)\n" );
  fprintf( stdout, "  --channel CHANNEL     set the channel name to CHANNEL\n");
  fprintf( stdout, "  --time TIME           set the used time spacing (001, 060 etc.)\n" );
  fprintf( stdout, "  --version VER         set the frame version name to to RUN (V01, V02, etc.)\n" );
  fprintf( stdout, "  --sensemon-format     read the text file in sensemon format\n" );
  fprintf( stdout, "  --skip-first-line     skip the first line of the file\n" );
  fprintf( stdout, "\nFactor File:\n" );
  fprintf( stdout,
      "The last argument must be a calibration factor file. This must be an ASCII\n"
      "file in the format:\n"
      "\n"
      "GPStime         alpha*beta      alpha           beta\n"
      "\n"
      "If the option --sensemon-format is given, the ASCII file must be in the\n"
      "format:\n"
      "\n"
      "GPStime         range           LineAmp         alpha           beta\n"
      "\n"
      "Any comment lines must begin with %% or #.\n" );
  return;
}

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
  const char *ver = NULL;
  const char *ifo = NULL;
  const char *timeSpace = NULL;
  int arg;
  int done = 0;
  char aname[64];
  char abname[64];
  struct series a;
  struct series ab;
  int code;
  char fname[256];
  char achannel[256];
  char bchannel[256];
  int dt;

  lalDebugLevel = 0;

  sensemon_format = 0;
  skip_first_line = 0;

  /* parse arguments */
  if ( argc == 1 )
  {
    usage( argv[0] );
    exit( 1 );
  }
  for ( arg = 1; arg < argc; ++arg )
  {
    if ( strstr( argv[arg], "--run" ) )
    {
      run = argv[++arg];
    }
    else if ( strstr( argv[arg], "--channel" ) )
    {
      channel = argv[++arg];
    }
    else if ( strstr( argv[arg], "--ifo" ) )
    {
      ifo = argv[++arg];
      site = *ifo;
    }
    else if ( strstr( argv[arg], "--time" ) )
    {
      timeSpace= argv[++arg];
      deltaT = atoi ( timeSpace );
    }
    else if ( strstr( argv[arg], "--version" ) )
    {
      ver = argv[++arg];
    }
    else if ( strstr( argv[arg], "--help" ) )
    {
      usage( argv[0] );
      exit( 0 );
    }
    else if ( strstr( argv[arg], "-h" ) )
    {
      usage( argv[0] );
      exit( 0 );
    }
    else if ( strstr( argv[arg], "--sensemon-format" ) )
    {
      sensemon_format = 1;
    }
    else if ( strstr( argv[arg], "--skip-first-line" ) )
    {
      skip_first_line = 1;
    }
    else if (strstr( argv[arg], "-V" ) )
      {
	/* print version information and exit */
        fprintf( stdout, "CVS Version: " CVS_ID_STRING "\n"
		 "CVS Tag: " CVS_NAME_STRING "\n" );
	exit(0);
      }
    else
    {
      if ( !ifo || ! run || ! ver || !channel )
      {
        fprintf( stderr, "Error: ifo, run, channel or version not specified\n" );
        exit( 1 );
      }
      if ( done )
      {
        fprintf( stderr,
            "Error: only one calibration can be generated at a time\n" );
        exit( 1 );
      }
      done = 1;
      /* output file names */
      sprintf( aname, "%s\\:" A_CHANNEL, ifo );
      sprintf( abname, "%s\\:" AB_CHANNEL, ifo );

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

      /* correct the end times */
      epoch_add( &a.tend, &a.tbeg, a.step * a.size );
      epoch_add( &ab.tend, &ab.tbeg, ab.step * ab.size );
      if ( ! frfile )
      {
        dt = (int)ceil( epoch_diff( &a.tend, &a.tbeg ) );
        sprintf( fname, "%c-%s_CAL_FAC_%s_%s_%s_%s-%d-%d.gwf",
		 *ifo, ifo, channel,
		 run, ver, timeSpace,
		 a.tbeg.gpsSeconds, dt );
        frfile = FrFileONew( fname, 0 );
      }

      /* set the correct channel names */
      if (strstr(run, "S5")) {
	sprintf( achannel,  A_CHANNEL "_%s", channel );
	sprintf( bchannel, AB_CHANNEL "_%s", channel );
      } else {
	sprintf( achannel,  A_CHANNEL  );
	sprintf( bchannel, AB_CHANNEL  );
      }


      /* don't mangle the channel names for frames */
      a.name  = aname;
      ab.name = abname;
      sprintf( aname, "%s:%s", ifo, achannel );
      sprintf( abname, "%s:%s",ifo, bchannel );
      frame = fr_add_proc_data( frame, &a );
      frame = fr_add_proc_data( frame, &ab );
    }
  }

  FrameWrite( frame, frfile );
  FrFileOEnd( frfile );
  return 0;
}
