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
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/*
 * simplereadcalframe.c
 *
 * Short program that will seek to a given GPS time and dump the relevant
 * calibration data that is valid at that time.
 *
 *   usage: simplereadcalframe -f framefile -e gpstime -c channel -d detector
 *       framefile (string): path to a calibration frame file
 *       gpstime (integer):  epoch (gps seconds) of calibration data to recover
 *       channel (string):   either "AS_Q" or "DARM_ERR"
 *       detector (string):  one of "L1", "H1", or "H2"
 *
 * This is sample code.  The intent is to be short and illustrative.
 *
 * Compile with:
 *
 *      cc -I/opt/lscsoft/libframe/include -o simplereadcalframe simplereadcalframe.c -L/opt/lscsoft/libframe/lib -lFrame -lm
 *
 */


#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lal/LALFrameL.h>

char *program;
void print_reference( FrStatData *sdata );
void print_factors( FrStatData *sdata );
void usage( int errcode );
void fail( const char *fmt, ... );

int main( int argc, char *argv[] )
{
  char *filename   = NULL;
  char *channel    = NULL;
  char *detector   = NULL;
  char oloopgain[] = "Xn:CAL-OLOOP_GAIN";
  char cavgain[]   = "Xn:CAL-CAV_GAIN_????????";
  char response[]  = "Xn:CAL-RESPONSE_????????";
  char oloopfac[]  = "Xn:CAL-OLOOP_FAC";
  char cavfac[]    = "Xn:CAL-CAV_FAC";
  FrFile     *frfile;
  FrameH     *frame;
  FrStatData *sdata;
  int gpstime = -1;
  int arg;

  program = argv[0];

  /* parse args */
  for ( arg = 1; arg < argc; ++arg )
  {
    if ( 0 == strcmp( argv[arg], "-h" ) )
      usage( 0 );
    else if ( 0 == strcmp( argv[arg], "-f" ) )
      filename = argv[++arg];
    else if ( 0 == strcmp( argv[arg], "-e" ) )
      gpstime = atoi( argv[++arg] );
    else if ( 0 == strcmp( argv[arg], "-c" ) )
      channel = argv[++arg];
    else if ( 0 == strcmp( argv[arg], "-d" ) )
      detector = argv[++arg];
    else
      usage( 1 );
  }
  /* make sure that required data has been found */
  if ( ! filename || ! channel || gpstime < 0 )
    usage( 1 );
  /* put appropriate readout points into channel names */
  if ( 0 == strcmp( channel, "AS_Q" ) || 0 == strcmp( channel, "DARM_ERR" ) )
  {
    strcpy( strchr( cavgain, '?' ), channel );
    strcpy( strchr( response, '?' ), channel );
  }
  else
    usage( 1 );
  /* put appropriate detector prefixes into channel names */
  if ( 0 == strcmp( detector, "L1" )
      || 0 == strcmp( detector, "H1" )
      || 0 == strcmp( detector, "H2" ) )
  {
    memcpy( oloopgain, detector, 2 );
    memcpy( cavgain, detector, 2 );
    memcpy( response, detector, 2 );
    memcpy( oloopfac, detector, 2 );
    memcpy( cavfac, detector, 2 );
  }
  else
    usage( 1 );

  /* open frame file and read frame */
  if ( ! ( frfile = FrFileINew( filename ) ) )
    fail( "could not open frame file %s\n", filename );
  if ( ! ( frame = FrameRead( frfile ) ) )
    fail( "could read frame from file %s\n", filename );

  /* write out open loop gain */
  if ( ! ( sdata = FrameFindStatData( frame, detector, oloopgain, gpstime ) ) )
    fail( "could not find static data %s for epoch %d in file %s\n",
        oloopgain, gpstime, filename );
  print_reference( sdata );

  /* write out cavity gain */
  if ( ! ( sdata = FrameFindStatData( frame, detector, cavgain, gpstime ) ) )
    fail( "could not find static data %s for epoch %d in file %s\n",
        cavgain, gpstime, filename );
  print_reference( sdata );

  /* write out response */
  if ( ! ( sdata = FrameFindStatData( frame, detector, response, gpstime ) ) )
    fail( "could not find static data %s for epoch %d in file %s\n",
        response, gpstime, filename );
  print_reference( sdata );

  /* write out open loop factors (gamma) */
  if ( ! ( sdata = FrameFindStatData( frame, detector, oloopfac, gpstime ) ) )
    fail( "could not find static data %s for epoch %d in file %s\n",
        oloopfac, gpstime, filename );
  print_factors( sdata );

  /* write out cavity factors (alpha) */
  if ( ! ( sdata = FrameFindStatData( frame, detector, cavfac, gpstime ) ) )
    fail( "could not find static data %s for epoch %d in file %s\n",
        cavfac, gpstime, filename );
  print_factors( sdata );

  FrFileIEnd( frfile );
  return 0;
}

void usage( int exitcode )
{
  fprintf( stderr, "usage: %s -f framefile -e gpstime -c channel -d detector\n",
      program );
  fprintf( stderr, "    framefile (string): path to a calibration frame file\n" );
  fprintf( stderr, "    gpstime (integer):  epoch (gps seconds) of calibration data to recover \n" );
  fprintf( stderr, "    channel (string):   either \"AS_Q\" or \"DARM_ERR\"\n" );
  fprintf( stderr, "    detector (string):  one of \"L1\", \"H1\", or \"H2\"\n" );
  exit( exitcode );
}

void fail( const char *fmt, ... )
{
  va_list ap;
  va_start( ap, fmt );
  vfprintf( stderr, fmt, ap );
  va_end( ap );
  exit( 1 );
}

char * generate_filename( FrStatData *sData )
{
  static char filename[FILENAME_MAX];
  char channel[FILENAME_MAX];
  char site;
  char *c;
  int ver;
  int t0;
  int dt;

  site = sData->detector->prefix[0];

  /* copy and transform channel name */
  strncpy( channel, sData->name, sizeof(channel) - 1 );
  while ( ( c = strpbrk( channel, "-: " ) ) )
    *c = '_';

  /* version, start, and end time */
  ver = sData->version;
  t0 = sData->timeStart;
  dt = sData->timeEnd - sData->timeStart;

  snprintf( filename, sizeof(filename), "%c-%s_V%d-%d-%d.dat", site, channel, ver, t0, dt );

  return filename;
}

void print_reference( FrStatData *sData )
{
  char *filename;
  FILE *fp;
  size_t i;
  filename = generate_filename( sData );
  fp = fopen( filename, "w" );
  for ( i = 1; i < sData->data->nData; ++i )
    fprintf( fp, "%.7e\t%.7e\t% .7e\n", i * sData->data->dx[0],
        hypot( sData->data->dataF[2*i+1], sData->data->dataF[2*i] ),
        atan2( sData->data->dataF[2*i+1], sData->data->dataF[2*i] ) );
  fclose( fp );
  return;
}

void print_factors( FrStatData *sData )
{
  char *filename;
  FILE *fp;
  size_t i;
  filename = generate_filename( sData );
  fp = fopen( filename, "w" );
  for ( i = 0; i < sData->data->nData; ++i )
    fprintf( fp, "%10d\t%.9f\n",
        (int)(sData->timeStart + sData->data->startX[0] + i*sData->data->dx[0]),
        sData->data->dataF[i] );
  fclose( fp );
  return;
}
