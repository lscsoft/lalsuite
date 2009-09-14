/*
*  Copyright (C) 2007 Drew Keppel
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
 * File Name: injcut.c
 *
 * Author: Keppel, D.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <config.h>

#if !defined HAVE_GSL_GSL_FFT_REAL_H || !defined HAVE_LIBGSL
int main( void ) { fprintf( stderr, "no gsl: disabled\n" ); return 77; }
#else

#include <math.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <lalapps.h>
#include <processtable.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOLwXML.h>
#include <lal/Date.h>
#include <lal/lalGitID.h>
#include <lalappsGitID.h>

RCSID( "$Id$" );
#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "injcut"

#define USAGE \
"lalapps_injcut [options]\n"\
"\nDefaults are shown in brackets\n\n" \
" [--help]                    display this message\n"\
" [--version]                 print version information and exit\n"\
" [--verbose]                 print progress information\n"\
"  --injection-file INJ_FILE  read injection parameters from INJ_FILE\n"\
"  --output OUTPUT            write output data to file: OUTPUT\n"\
" [--mass-cut] MASS_TYPE      keep only triggers in mass range of type\n"\
"                             MASS_TYPE (mtotal|mchirp|mcomp)\n"\
"                             (if MASS_TYPE = mcomp, mass2 is the smaller mass)\n"\
" [--mass-range-low] MIN      set the minimum mass to MIN\n"\
" [--mass-range-high] MAX     set the maximum mass to MAX\n"\
" [--mass2-range-low] MIN     set the minimum mass2 to MIN\n"\
" [--mass2-range-high] MAX    set the maximum mass2 to MAX\n"\
"\n"


/* all units are in kpc since this is what GalacticInspiralParamStruc expects */
static ProcessParamsTable *next_process_param( 
        const char *name, 
        const char *type,
            const char *fmt, ... );

extern int vrbflg;


ProcessParamsTable *next_process_param( 
        const char *name, 
        const char *type,
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
  if ( ! strcmp( name, "userTag" ) || ! strcmp( name, "user-tag" ) )
    snprintf( pp->param, LIGOMETA_PARAM_MAX, "-userTag" );
  else
    snprintf( pp->param, LIGOMETA_PARAM_MAX, "--%s", name );
  strncpy( pp->type, type, LIGOMETA_TYPE_MAX );
  va_start( ap, fmt );
  vsnprintf( pp->value, LIGOMETA_VALUE_MAX, fmt, ap );
  va_end( ap );

  return pp;
}


int main( int argc, char *argv[] )
{
  LALStatus             status = blank_status;
  /* command line options */
  CHAR         *massCut = NULL;     /* mass cut type */
  REAL4         minMass = -1;       /* minimum mass */
  REAL4         maxMass = -1;       /* maximum mass */
  REAL4         minMass2 = -1;       /* minimum mass2 */
  REAL4         maxMass2 = -1;       /* maximum mass2 */
  CHAR         *injectFileName = NULL;
  CHAR         *outputFileName = NULL;

  /* program variables */
  INT4   numSimEvents = 0;

  /* xml input data */
  SimInspiralTable     *simEventHead = NULL;

  /* xml output data */
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         injections;
  ProcessParamsTable   *this_proc_param;
  SimInspiralTable     *this_inj = NULL;
  LIGOLwXMLStream       xmlfp;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"help",                    no_argument,       0,                'h'},
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"version",                 no_argument,       0,                'V'},
    {"debug-level",             required_argument, 0,                'z'},
    {"injection-file",          required_argument, 0,                'I'},
    {"mass-cut",                required_argument, 0,                'M'},
    {"mass-range-low",          required_argument, 0,                't'},
    {"mass-range-high",         required_argument, 0,                'T'},
    {"mass2-range-low",         required_argument, 0,                'r'},
    {"mass2-range-high",        required_argument, 0,                'R'},
    {"output",                  required_argument, 0,                'o'},
    {0, 0, 0, 0}
  };
  int c;

  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );

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


  /*
   *
   * parse command line arguments
   *
   */

     
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "c:C:hI:m:M:o:t:T:Vz:Z:", long_options, &option_index );

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

      case 'I':
        /* create storage for the injection file name */
        optarg_len = strlen( optarg ) + 1;
        injectFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( injectFileName, optarg, optarg_len );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "%s", injectFileName );
        break;

      case 'M':
        /* create storage for the missed injection file name */
        optarg_len = strlen( optarg ) + 1;
        massCut = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( massCut, optarg, optarg_len );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "%s", massCut );
        break;

      case 'o':
        /* create storage for the output file name */
        optarg_len = strlen( optarg ) + 1;
        outputFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( outputFileName, optarg, optarg_len );
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "%s", outputFileName );
        break;

      case 't':
        /* minimum total mass */
        minMass = (REAL4) atof( optarg );
        if ( minMass <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "miniumum mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, minMass );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%e", minMass );
        break;

      case 'T':
        /* maximum total mass */
        maxMass = (REAL4) atof( optarg );
        if ( maxMass <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "maxiumum mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, maxMass );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%e", maxMass );
        break;

      case 'r':
        /* minimum total mass2 */
        minMass2 = (REAL4) atof( optarg );
        if ( minMass2 <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "miniumum mass2 must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, minMass2 );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%e", minMass2 );
        break;

      case 'R':
        /* maximum total mass2 */
        maxMass2 = (REAL4) atof( optarg );
        if ( maxMass2 <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "maxiumum mass2 must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, maxMass2 );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%e", maxMass2 );
        break;

      case 'z':
        set_debug_level( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
            "string", "%s", optarg );
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "INJection CUTter routine\n" 
            "Drew Keppel\n"
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n" );
        fprintf( stdout, lalappsGitID );
        exit( 0 );
        break;
        
      case 'h':
      case '?':
        fprintf( stderr, USAGE );
        exit( 0 );
        break;

      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        fprintf( stderr, USAGE );
        exit( 1 );
    }
  }

  if ( optind < argc )
  {
    fprintf( stderr, "extraneous command line arguments:\n" );
    while ( optind < argc )
    {
      fprintf ( stderr, "%s\n", argv[optind++] );
    }
    exit( 1 );
  }

  /* check that the input and output file names have been specified */
  if ( ! injectFileName )
  {
    fprintf( stderr, "--injection-file must be specified\n" );
    exit( 1 );
  }
  if ( ! outputFileName )
  {
    fprintf( stderr, "--output must be specified\n" );
    exit( 1 );
  }

  /* check that the mass inputs make sense */
  if ( ( massCut || minMass >= 0 || maxMass >= 0 ) &&
       ! ( massCut && minMass >= 0 && maxMass >= 0 ) )
  {
    fprintf( stderr, "--mass-cut, --mass-range-low, and --mass-rang-high "
        "must all be used together\n" );
    exit( 1 );
  }

  if ( massCut && ( minMass >= maxMass ) )
  {
    fprintf( stderr, "--mass-range-low must be less than "
        "--mass-range-high\n" );
    exit( 1 );
  }

  if ( massCut && ( ! strcmp( "mcomp", massCut ) ) &&
       ( minMass2 < 0 || maxMass2 <0 ) )
  {
    fprintf( stderr, "--mass2-range-low and --mass2-rang-high \n"
        "must be specified if using --mass-cut mcomp\n" );
    exit( 1 );
  }

  if ( massCut && ( ! strcmp( "mcomp", massCut ) ) &&
       minMass2 >= maxMass2 )
  {
    fprintf( stderr, "--mass2-range-low must be less than "
        "--mass2-range-high\n" );
    exit( 1 );
  }

  /* null out the head of the linked list */
  injections.simInspiralTable = NULL;

  /*
   *
   * read in the injection XML file
   *
   */

  if ( vrbflg )
  {
    fprintf( stdout, "reading injections from %s... ", injectFileName );
  }

  numSimEvents = SimInspiralTableFromLIGOLw( &simEventHead, 
      injectFileName, 0, 0 );

  if ( vrbflg ) fprintf( stdout, "got %d injections\n", numSimEvents );

  if ( numSimEvents < 0 )
  {
    fprintf( stderr, "error: unable to read sim_inspiral table from %s\n",
        injectFileName );
    exit( 1 );
  }

  /* keep only injections in chirp mass range */
  if ( ( minMass >= 0 ) && ( maxMass >= 0 ) && ! strcmp( "mchirp", massCut ) )
  {
    numSimEvents = XLALSimInspiralChirpMassCut( &simEventHead,
        minMass, maxMass );
    if ( vrbflg )
    {
      fprintf( stdout, "kept %d injections chirp mass range %f to %f\n",
          numSimEvents, minMass, maxMass );
    }
  }

  /* keep only injections in component mass range */
  if ( ( minMass >= 0 ) && ( maxMass >= 0 ) &&
       ( minMass2 >= 0 ) && ( maxMass2 >= 0 ) && ! strcmp( "mcomp", massCut ) )
  {
    numSimEvents = XLALSimInspiralCompMassCut( &simEventHead,
        minMass, maxMass, minMass2, maxMass2 );
    if ( vrbflg )
    {
      fprintf( stdout, "kept %d injections component mass range %f to %f\n",
          numSimEvents, minMass, maxMass );
    }
  }

  /* keep only injections in total mass range */
  if ( ( minMass >= 0 ) && ( maxMass >= 0 ) && ! strcmp( "mtotal", massCut ) )
  {
    numSimEvents = XLALSimInspiralTotalMassCut( &simEventHead,
        minMass, maxMass );
    if ( vrbflg )
    {
      fprintf( stdout, "kept %d injections total mass range %f to %f\n",
          numSimEvents, minMass, maxMass );
    }
  }

  injections.simInspiralTable = simEventHead;

  /* open the xml file */
  if ( vrbflg ) fprintf( stdout, "writing output xml files... " );
  memset( &xmlfp, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlfp, outputFileName ),
      &status );

  /* write the process table */
  snprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "H1H2L1" );
  XLALGPSTimeNow(&(proctable.processTable->end_time));
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, process_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, proctable, 
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
  free( proctable.processTable );

  /* free the unused process param entry */
  this_proc_param = procparams.processParamsTable;
  procparams.processParamsTable = procparams.processParamsTable->next;
  free( this_proc_param );

  /* write the process params table */
  if ( procparams.processParamsTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, process_params_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, procparams, 
          process_params_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
    while( procparams.processParamsTable )
    {
      this_proc_param = procparams.processParamsTable;
      procparams.processParamsTable = this_proc_param->next;
      free( this_proc_param );
    }
  }

  /* write the sim_inspiral table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, sim_inspiral_table ), 
        &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, injections, 
          sim_inspiral_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
  while ( injections.simInspiralTable )
  {
    this_inj = injections.simInspiralTable;
    injections.simInspiralTable = injections.simInspiralTable->next;
    LALFree( this_inj );
  }

  /* close the injection file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlfp ), &status );

  /* check for memory leaks and exit */
  LALCheckMemoryLeaks();
  return 0;
}
#endif
