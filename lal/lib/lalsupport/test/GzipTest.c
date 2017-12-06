/*
*  Copyright (C) 2013 Matthew Pitkin
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

#include <getopt.h>

#include <lal/LALStdlib.h>
#include <lal/XLALError.h>
#include <lal/LALString.h>
#include <lal/FileIO.h>

#define USAGE \
"Usage: %s [options]\n\n"\
" --help (-h)              display this message\n"\
" --file (-f)              name of ascii text file to gzip/unzip\n"\
" --gzip (-g)              gzip (compress) the text file\n"\
" --gunzip (-u)            gunzip (decompress) the gzipped text file\n"\
"\n"

int main(int argc, char **argv){
  CHAR *filename = NULL;
  INT4 gzip = 0, gunzip = 0;

  struct option long_options[] =
  {
    { "help",      no_argument,        0, 'h' },
    { "file",      required_argument,  0, 'f' },
    { "gzip",      no_argument,     NULL, 'g' },
    { "guzip",     no_argument,     NULL, 'u' },
    { 0, 0, 0, 0 }
  };

  CHAR args[] = "hf:gu";
  CHAR *program = argv[0];

  /* get input arguments */
  while(1){
    int option_index = 0;
    int c;

    c = getopt_long( argc, argv, args, long_options, &option_index );
    if ( c == -1 ) /* end of options */
      break;

    switch(c){
      case 0: /* if option set a flag, nothing else to do */
        if ( long_options[option_index].flag )
          break;
        else
          fprintf(stderr, "Error parsing option %s with argument %s\n", long_options[option_index].name, optarg );
      case 'h': /* help message */
        fprintf(stderr, USAGE, program);
        exit(0);
      case 'f': /* input file */
        filename = XLALStringDuplicate( optarg );
        break;
      case 'g': /* gzip the file */
        gzip = 1;
        break;
      case 'u': /* gunzip the file */
        gunzip = 1;
        break;
      case '?':
        fprintf(stderr, "Unknown error while parsing options\n" );
        exit(0);
      default:
        fprintf(stderr, "Unknown error while parsing options\n" );
        exit(0);
    }
  }

  if ( filename == NULL ){
    fprintf(stderr, "Must specify an input file\n");
    fprintf(stderr, USAGE, program);
    exit(0);
  }

  if ( ( !gzip && !gunzip ) || ( gzip && gunzip ) ){
    fprintf(stderr, "Must specify whether you want to either gzip (-g) or gunzip (-u) the input file.\n");
    fprintf(stderr, USAGE, program);
    exit(0);
  }

  if ( gzip ){
    /* zip the file */
    if ( XLALGzipTextFile( filename ) != XLAL_SUCCESS ){
      fprintf(stderr, "Gzip of %s has failed\n", filename);
      exit(1);
    }
  }
  if ( gunzip ){
    /* unzip the file */
    if ( XLALGunzipTextFile( filename ) != XLAL_SUCCESS ){
      fprintf(stderr, "Guzip of %s has failed\n", filename);
      exit(1);
    }
  }

  return 0;
}
