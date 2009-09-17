/*
*  Copyright (C) 2007 Duncan Brown
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
 * File Name: inspxmlinj.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <config.h>
#include <lalapps.h>
#include <processtable.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/Date.h>
#include <lal/TimeDelay.h>
#include <lal/lalGitID.h>
#include <lalappsGitID.h>

#define USAGE \
  "lalapps_inspxmlinj [options]\n"\
"\nDefaults are shown in brackets\n\n" \
"  --help                  display this message\n"\
"  --end-time SEC          geocentric end time seconds (700000000)\n"\
"  --end-time-ns SEC       geocentric end time nanoseconds (0)\n"\
"  --mass-1 MASS1          first binary component mass in solar masses (1.4)\n" \
"  --mass-2 MASS2          second binary component mass in solar masses (1.4)\n" \
"  --distance MPC          distance from center of earth to binary in Mpc (1)\n"\
"  --longitude LONG        longitude in specified coordinate system (0) \n"\
"  --latitude LAT          latitude in specified coordinate system (0)\n"\
"  --coordinate-system C   sky coordinates [horizon|galactic|geocentric](horizon)\n"\
"  --inclination IOTA      binary inclination angle (0)\n"\
"  --coalescence-phase PHI binary coalescence phase (0)\n"\
"  --polarization PSI      binary polarization angle (0)\n"\
"  --waveform NAME         set waveform type to NAME (GeneratePPNtwoPN)\n"\
"  --user-tag STRING       set the usertag to STRING\n"\
"\n"


RCSID( "$Id$" );

#define MPC ( 1e6 * LAL_PC_SI )

#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "inspxmlinj"

ProcessParamsTable *next_process_param( const char *name, const char *type,
    const char *fmt, ... )
{
  ProcessParamsTable *pp;
  va_list ap;
  pp = calloc( 1, sizeof( *pp ) );
  if ( ! pp )
  {
    perror( "next_process_param" );
    exit( 1 );
  }
  strncpy( pp->program, PROGRAM_NAME, LIGOMETA_PROGRAM_MAX );
  snprintf( pp->param, LIGOMETA_PARAM_MAX, "--%s", name );
  strncpy( pp->type, type, LIGOMETA_TYPE_MAX );
  va_start( ap, fmt );
  vsnprintf( pp->value, LIGOMETA_VALUE_MAX, fmt, ap );
  va_end( ap );
  return pp;
}

int main ( int argc, char *argv[] )
{
  /* xml output data */
  CHAR                  fname[256];
  CHAR                 *userTag = NULL;
  LALStatus             status = blank_status;
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         injections;
  SimInspiralTable      injParams;
  ProcessParamsTable   *this_proc_param;
  LIGOLwXMLStream       xmlfp;

  /* default sky position */
  SkyPosition           binaryPosition;

  struct option long_options[] =
  {
    {"help",                    no_argument,       0,                'H'},
    {"end-time",                required_argument, 0,                'a'},
    {"end-time-ns",             required_argument, 0,                'b'},
    {"mass-1",                  required_argument, 0,                'c'},
    {"mass-2",                  required_argument, 0,                'd'},
    {"distance",                required_argument, 0,                'e'},
    {"longitude",               required_argument, 0,                'f'},
    {"latitude",                required_argument, 0,                'g'},
    {"coordinate-system",       required_argument, 0,                'h'},
    {"inclinaton",              required_argument, 0,                'i'},
    {"coalescence-phase",       required_argument, 0,                'j'},
    {"polarization",            required_argument, 0,                'k'},
    {"waveform,"                required_argument, 0,                'l'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {0, 0, 0, 0}
  };
  int c;

  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "LALMSGLVL2" );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) 
    calloc( 1, sizeof(ProcessTable) );
  XLALGPSTimeNow(&(proctable.processTable->start_time));
  if (strcmp(CVS_REVISION,"$Revi" "sion$"))
    {
      LAL_CALL( populate_process_table( &status, proctable.processTable, 
                                        PROGRAM_NAME, CVS_REVISION,
                                        CVS_SOURCE, CVS_DATE ), &status );
    }
  else
    {
      LAL_CALL( populate_process_table( &status, proctable.processTable, 
                                        PROGRAM_NAME, lalappsGitCommitID,
                                        lalappsGitGitStatus,
                                        lalappsGitCommitDate ), &status );
    }
  snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );

  /* clear the waveform field */
  memset( &(injParams.waveform), 0, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR) );

  /* set the default values */
  injParams.mass1 = 1.4;
  injParams.mass2 = 1.4;
  injParams.eta = 0.25;
  injParams.distance = 1.0;
  injParams.inclination = 0;
  injParams.coa_phase = 0;
  injParams.polarization = 0;

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    long int gpsinput;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "Ha:b:c:d:e:f:g:h:i:j:k:l:Z:", long_options, &option_index );

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

      case 'a':
        {
          long int gendsec = atol( optarg );
          if ( gendsec < 441417609 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS start time is prior to " 
                "Jan 01, 1994  00:00:00 UTC:\n"
                "(%ld specified)\n",
                long_options[option_index].name, gendsec );
            exit( 1 );
          }
          if ( gendsec > 999999999 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS start time is after " 
                "Sep 14, 2011  01:46:26 UTC:\n"
                "(%ld specified)\n", 
                long_options[option_index].name, gendsec );
            exit( 1 );
          }
          injParams.geocent_end_time.gpsSeconds = (INT4) gendsec;
          this_proc_param = this_proc_param->next = 
            next_process_param( long_options[option_index].name, "int", 
                "%ld", gendsec );
        }
        break;

      case 'b':
        {
          long int gendnansec = atol( optarg );
          if ( gendnansec < 0 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS start time nanoseconds is negative\n",
                long_options[option_index].name );
            exit( 1 );
          }
          if ( gendnansec > 999999999 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS start time nanoseconds is greater than unity:\n" 
                "Must be <= 999999999 (%ld specified)\n", 
                long_options[option_index].name, gendnansec );
            exit( 1 );
          }
          injParams.geocent_end_time.gpsNanoSeconds = (INT4) gendnansec;
          this_proc_param = this_proc_param->next = 
            next_process_param( long_options[option_index].name, "int", 
                "%ld", gendnansec );
        }
        break;

      case 'c':
        injParams.mass1 = (REAL4) atof( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "real_4", 
              "%d", rand_seed );
        break;

      case 'd':
        injParams.mass2 = (REAL4) atof( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "real_4", 
              "%d", rand_seed );
        break;

      case 'e':
        injParams.distance = (REAL4) atof( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "real_4", 
              "%d", rand_seed );
        break;

      case 'f':
        abort();
        break;

      case 'g':
        abort();
        break;

      case 'h':
        abort();
        break;

      case 'i':
        injParams.inclination = (REAL4) atof( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "real_4", 
              "%d", rand_seed );
        break;

      case 'j':
        injParams.coa_phase= (REAL4) atof( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "real_4", 
              "%d", rand_seed );
        break;

      case 'k':
        injParams.polarization = (REAL4) atof( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "real_4", 
              "%d", rand_seed );
        break;

      case 'l':
        snprintf( &(injParams.waveform), 
                  LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "%s", optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "string", 
              "%s", optarg );
        break;

      case 'Z':
        /* create storage for the usertag */
        optarg_len = strlen( optarg ) + 1;
        userTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( userTag, optarg, optarg_len );

        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
                  PROGRAM_NAME );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "-userTag" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
                  optarg );
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
    }
  }

  if ( ! *waveform )
  {
    /* default to Tev's GeneratePPNInspiral as used in */
    snprintf( waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), 
              "GeneratePPNtwoPN" );
  }





