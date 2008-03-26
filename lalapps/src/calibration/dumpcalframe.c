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
 * dumpcalframe.c
 *
 * Dumps the contents of a calibration frame file to ascii files.
 *
 * Compile using:
 *
 *      cc -I/opt/lscsoft/libframe/include -o dumpcalframe dumpcalframe.c -L/opt/lscsoft/libframe/lib -lFrame -lm
 *
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <FrameL.h>


char * generate_filename( FrStatData *sData );
void print_reference( FrStatData *sData );
void print_factors( FrStatData *sData );

int main( int argc, char *argv[] )
{
  int arg;

  if ( argc == 1 )
  {
    fprintf( stderr, "usage: %s file1 [file2 ...]\n", argv[0] );
    return 1;
  }
  
  for ( arg = 1; arg < argc; ++arg )
  {
    FrFile *frfile = FrFileINew( argv[arg] );
    FrameH *frame;
    fprintf( stderr, "File %s\n", argv[arg] );
    while ( ( frame = FrameRead( frfile ) ) )
    {
      FrDetector *det;
      fprintf( stderr, "\tFrame %d\n", frame->frame );
      for ( det = frame->detectProc; det != NULL; det = det->next )
      {
        FrStatData *sData;
        fprintf( stderr, "\t\tDetector %s\n", det->name );
        for ( sData = det->sData; sData != NULL; sData = sData->next )
        {
          if ( 0 == strcmp( sData->representation, "CAL_REF" ) )
          {
            fprintf( stderr, "\t\t\tCalibration Data %s\n", sData->name );
            fprintf( stderr, "\t\t\t\tVersion %d\n", sData->version );
            fprintf( stderr, "\t\t\t\tTime %d-%d\n", sData->timeStart, sData->timeEnd );
            print_reference( sData );
          }
          else if ( 0 == strcmp( sData->representation, "CAL_FAC" ) )
          {
            fprintf( stderr, "\t\t\tCalibration Data %s\n", sData->name );
            fprintf( stderr, "\t\t\t\tVersion %d\n", sData->version );
            fprintf( stderr, "\t\t\t\tTime %d-%d\n", sData->timeStart, sData->timeEnd );
            print_factors( sData );
          }
          else if ( 0 == strcmp( sData->representation, "time_series" ) && 0 !=  strcmp( sData->comment, "time domain FIR filter" ))
          {
            fprintf( stderr, "\t\t\tCalibration Data %s\n", sData->name );
            fprintf( stderr, "\t\t\t\tVersion %d\n", sData->version );
            fprintf( stderr, "\t\t\t\tTime %d-%d\n", sData->timeStart, sData->timeEnd );
            print_factors( sData );
          }
        }
      }
      FrameFree( frame );
    }
    FrFileIEnd( frfile );
  }

  return 0;
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
