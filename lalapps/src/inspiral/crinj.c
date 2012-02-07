/*
*  Copyright (C) 2007 Chad Hanna
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

/*----------------------------------------------------------------------- 
 * 
 * File Name: crinj.c
 *
 * Author: Ravi Kumar and Chad Hanna
 * 
 * 
 *-----------------------------------------------------------------------
 */


#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/TimeDelay.h>


#define USAGE \
  "lalapps_inspinj [options]\n"\
"\nDefaults are shown in brackets\n\n" \
"  --help                      display this message\n"\
"  --source-file FILE          read source parameters from FILE\n"\
"  --RA-bin-size DEGREES       input right ascension bin size in degrees\n"\
"  --cosDEC-bin-size (0,1]     input Cos(Declination) bin size from +0,1\n"\
"  --logDIST-bin-size (0,3]    input log10(Distance) bin size from +0,3\n"\
"\n"

#define KPC ( 1e3 * LAL_PC_SI )
#define MPC ( 1e6 * LAL_PC_SI )
#define GPC ( 1e9 * LAL_PC_SI )

struct {
  char   name[16];
  double ra;
  double dec;
  double dist;
  double lum;
  double fudge;
} *source_data;

typedef struct {
  char NAME[16];
  double RA;
  double cosDEC;
  double logDIST;
  double DEC;
  double DIST;
  double LUM;
  double FUDGE;
  int RAbin;
  int cosDECbin;
  int logDISTbin;
  unsigned long long ID;
} *binsourcedata;


int num_source = 0;

void read_source_data( char *sourceFileName )
{
  char line[256];
  FILE *fp;
  int i;
  fp = fopen( sourceFileName, "r" );


  if ( ! fp )
  {
    perror( "read_source_data" );
    printf( "Could not find file %s\n", sourceFileName );
    exit( 1 );
  }

  num_source = 0;
  while ( fgets( line, sizeof( line ), fp ) )
    if ( line[0] == '#' )
      continue;
    else 
      ++num_source;
  rewind( fp );
  
  source_data = calloc( num_source, sizeof( *source_data ) );

  if ( ! source_data )
  {
    printf("alloc error\n" );
    exit( 1 );
  }

  i = 0;
  while ( fgets( line, sizeof( line ), fp ) )
    if ( line[0] == '#' )
      continue;
    else
    {
      char ra_sgn, dec_sgn;
      double ra_h, ra_m, dec_d, dec_m;
      int c;

      c = sscanf( line, "%s %c%le:%le %c%le:%le %le %le %le",
          source_data[i].name, &ra_sgn, &ra_h, &ra_m, &dec_sgn, &dec_d, &dec_m,
          &source_data[i].dist, &source_data[i].lum, &source_data[i].fudge );

      if ( c != 10 )
      {
        fprintf( stderr, "error parsing source datafile %s line %d\n", sourceFileName,i );
        exit( 1 );
      }

      /* by convention, overall sign is carried only on hours/degrees entry */
      source_data[i].ra  = ( ra_h + ra_m / 60.0 ) * LAL_PI / 12.0;
      source_data[i].dec = ( dec_d + dec_m / 60.0 ) * LAL_PI / 180.0;

      if ( ra_sgn == '-' )
        source_data[i].ra *= -1;
      if ( dec_sgn == '-' )
        source_data[i].dec *= -1;
      ++i;
    }
  fclose( fp );
}

int main( int argc, char *argv[] )
{
  /* Declare some local variables */

  /* set up inital debugging values */
  float RAbinsize = 0;
  float cosDECbinsize = 0;
  float logDISTbinsize = 0;
  FILE *CFP = NULL;                             /* pointer to the cluster file */
  char *sourceFileName = NULL;
  int i = 0;
  int n = 0;
  int m = 0;
  binsourcedata bin_source_data = NULL;
  double totalLUM = 0;
  double sumDIST =0;
  double sumRA = 0;
  double sumcosDEC =0;
  
 
  char   NAME[16];
  double RA = 0;
  double DEC = 0;
  double DIST =0;
  double cosDEC = 0;
  double logDIST = 0;
  double LUM = 0;
  double RAbin = 0;
  double cosDECbin = 0;
  double logDISTbin = 0;
  unsigned long long ID = 0;

  char   cNAME[16];
  double cRA = 0;
  double cDEC = 0;
  double cDIST = 0;
  double cLUM = 0;
  char bigNAME[16];
 
  int num_c = 0;
  

/*
  struct bin_source_data *bin_source_data = NULL;
  struct clusteredsourcedata *clustered_source_data = NULL;
*/
  /* getopt arguments */
  struct option long_options[] =
  {
    {"help",                          no_argument, 0,                'h'},
    {"source-file",             required_argument, 0,                'f'},
    {"RA-bin-size",             required_argument, 0,                'r'},
    {"cosDEC-bin-size",         required_argument, 0,                'c'},
    {"logDIST-bin-size",        required_argument, 0,                'l'},
    {0, 0, 0, 0}
  };
  int c;


  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "hf:r:c:l:", long_options, &option_index );

    /* detect the end of the options */
    if ( c == - 1 )
    {
      break;
    }

    switch ( c )
    {
      case 0:
        /* if this option set a flag, do nothing else now */
        if ( long_options[option_index].flag != 0 )
        {
          break;
        }
        else
        {
          fprintf( stderr, "error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;

      case 'f':
        optarg_len = strlen( optarg ) + 1;
        sourceFileName = calloc( 1, optarg_len * sizeof(char) );
        memcpy( sourceFileName, optarg, optarg_len * sizeof(char) );
        break;

      case 'r':
        RAbinsize = atof( optarg );
        RAbinsize *= LAL_PI/180.0;
        break;

      case 'c':
        cosDECbinsize = atof( optarg );
        break;

      case 'l':
        logDISTbinsize = atof( optarg ); 
        break;

      case 'h':
        fprintf( stderr, USAGE );
        exit( 0 );
        break;

      case '?':
        fprintf( stderr, USAGE );
        exit( 1 );
        break;

      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        fprintf( stderr, USAGE );
        exit( 1 );
        break;
    }
  }

  read_source_data( sourceFileName );

  printf("finished reading sourcelist...numsource = %d\n",num_source); 
  bin_source_data = calloc(num_source, sizeof(*bin_source_data));
  

  /* binning code */
   i=0; 
  /* bin the source file */
  for (i=0;i < num_source; i++){
    
  memcpy(bin_source_data[i].NAME,source_data[i].name,sizeof(source_data[i].name));

    bin_source_data[i].LUM = source_data[i].lum;
    bin_source_data[i].RA = source_data[i].ra;
    bin_source_data[i].DEC = source_data[i].dec;
    bin_source_data[i].DIST = source_data[i].dist;
    if (source_data[i].dec < 0)
      bin_source_data[i].cosDEC = -1.0*cos(source_data[i].dec);
    else
      bin_source_data[i].cosDEC = cos(source_data[i].dec);
    bin_source_data[i].logDIST = log10(source_data[i].dist);
    bin_source_data[i].RAbin = floor(bin_source_data[i].RA/RAbinsize);
    bin_source_data[i].cosDECbin = floor((1.0+bin_source_data[i].cosDEC)/cosDECbinsize);
    bin_source_data[i].logDISTbin = floor(bin_source_data[i].logDIST/logDISTbinsize);
    bin_source_data[i].ID = 10000*bin_source_data[i].RAbin + 
                         100*bin_source_data[i].cosDECbin +
                         1000000*bin_source_data[i].logDISTbin;

    }

  /* Sort the list */
    
    for(n = 0; n < num_source; n++){
      for(m = (n+1); m < num_source; m++){
        if (bin_source_data[m].ID < bin_source_data[n].ID){
          memcpy(NAME, bin_source_data[m].NAME,sizeof(bin_source_data[m].NAME));
          LUM = bin_source_data[m].LUM;
          RA = bin_source_data[m].RA;
          DEC = bin_source_data[m].DEC;
          DIST = bin_source_data[m].DIST;
          cosDEC = bin_source_data[m].cosDEC;
          logDIST = bin_source_data[m].logDIST;
          RAbin = bin_source_data[m].RAbin;
          cosDECbin = bin_source_data[m].cosDECbin;
          logDISTbin = bin_source_data[m].logDISTbin;
          ID = bin_source_data[m].ID;
          memcpy(bin_source_data[m].NAME, bin_source_data[n].NAME,sizeof(bin_source_data[n].NAME));
          bin_source_data[m].LUM = bin_source_data[n].LUM;
          bin_source_data[m].RA = bin_source_data[n].RA;
          bin_source_data[m].DEC = bin_source_data[n].DEC;
          bin_source_data[m].DIST = bin_source_data[n].DIST;
          bin_source_data[m].cosDEC = bin_source_data[n].cosDEC;
          bin_source_data[m].logDIST = bin_source_data[n].logDIST;
          bin_source_data[m].RAbin = bin_source_data[n].RAbin;
          bin_source_data[m].cosDECbin = bin_source_data[n].cosDECbin;
          bin_source_data[m].logDISTbin = bin_source_data[n].logDISTbin;
          bin_source_data[m].ID = bin_source_data[n].ID;
          memcpy(bin_source_data[n].NAME, NAME,sizeof(NAME));
          bin_source_data[n].LUM = LUM;
          bin_source_data[n].RA = RA;
          bin_source_data[n].DEC = DEC;
          bin_source_data[n].DIST = DIST;
          bin_source_data[n].cosDEC = cosDEC;
          bin_source_data[n].logDIST = logDIST;
          bin_source_data[n].RAbin = RAbin;
          bin_source_data[n].cosDECbin = cosDECbin;
          bin_source_data[n].logDISTbin = logDISTbin;
          bin_source_data[n].ID = ID;

        }
      }
    }

  /* Clustering and write it to a file called ClusterList.dat*/

 
  n = 0;
  CFP = fopen("ClusterList.dat", "w");
  if (CFP == NULL){
    printf("error opening ClusterList.dat for writing\n");
    return 1;
    }

  
  fprintf(CFP, "# Clustered galaxy list\n");
  fprintf(CFP, "# NAME\t\t RA\t\t DEC\t\t DIST\t\t LUM\t\t FUDGE\n");

  totalLUM = bin_source_data[0].LUM;
  sumDIST = bin_source_data[0].DIST * bin_source_data[0].LUM;
  sumRA = bin_source_data[0].RA * bin_source_data[0].LUM;
  sumcosDEC = bin_source_data[0].cosDEC * bin_source_data[0].LUM;
  memcpy(bigNAME, bin_source_data[0].NAME,sizeof(bin_source_data[0].NAME));
  
  for(n=1; n<num_source; n++){
      if (bin_source_data[n].ID != bin_source_data[n-1].ID){
        num_c++;
        cLUM = totalLUM;
        cRA = sumRA/totalLUM;
        /*Check Sign and bin over cos(DEC)*/
        if( sumcosDEC < 0)
          cDEC = -acos(sumcosDEC/totalLUM);
        else
          cDEC = acos(sumcosDEC/totalLUM);
        cDIST = sumDIST/totalLUM;
        memcpy(cNAME, bigNAME,sizeof(bigNAME));
        fprintf(CFP, "c%-10s\t %+02d:%0.2f\t %+02d:%0.2f\t %8.1f\t %8.3f\t %.2f\n", 
                cNAME,
                (int) floor(cRA*12.0/LAL_PI), 
                (fmod(cRA*12.0/LAL_PI,1.0)*60.0),
                (int) floor(cDEC*180.0/LAL_PI), 
                (fabs(fmod(cDEC*180.0/LAL_PI,1.0)*60.0)), 
                cDIST, cLUM, 1.00 );
        /* initialize */
        memcpy(bigNAME, bin_source_data[n].NAME,sizeof(bin_source_data[n].NAME));
        totalLUM = bin_source_data[n].LUM;
        sumDIST = bin_source_data[n].DIST*bin_source_data[n].LUM;
        sumRA = bin_source_data[n].RA*bin_source_data[n].LUM;
        sumcosDEC = bin_source_data[n].cosDEC*bin_source_data[n].LUM;
      }
      else
        { 
        if (bin_source_data[n].LUM > bin_source_data[n-1].LUM)
          memcpy(bigNAME, bin_source_data[n].NAME,sizeof(bin_source_data[n].NAME));
        totalLUM += bin_source_data[n].LUM;
        sumDIST += bin_source_data[n].DIST*bin_source_data[n].LUM;
        sumRA += bin_source_data[n].RA*bin_source_data[n].LUM;
        sumcosDEC += bin_source_data[n].cosDEC*bin_source_data[n].LUM;
      }
  }

printf("Clustered list written to ClusterList.dat\n");
printf("Clustered list size = %d\n", num_c);
printf("Fraction of sources = %f\n", ((float) num_c) / ((float)num_source));


fclose(CFP);
return 0;
}
