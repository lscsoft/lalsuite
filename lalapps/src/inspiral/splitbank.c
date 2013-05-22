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
 * File Name: splitbank.c
 *
 * Author: Brown, D. A.
 * 
 * 
 *-----------------------------------------------------------------------
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <time.h>

#include <lalapps.h>
#include <processtable.h>

#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/Date.h>

#include <LALAppsVCSInfo.h>

#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "splitbank"

#define USAGE \
"Usage: %s [options] [LIGOLW XML input files]\n\n"\
"  --help                    display this message\n"\
"  --verbose                 print progress information\n"\
"  --version                 print version information\n"\
"  --user-tag STRING         set the process_params usertag to STRING\n"\
"  --comment STRING          set the process table comment to STRING\n"\
"\n"\
"  --bank-file FILE          read template bank parameters from FILE\n"\
"  --number-of-banks N       split template bank into N files\n"\
"  --minimal-match M         set minimal match of triggered bank to M\n"\

extern int vrbflg;                      /* verbocity of lal function    */


int main ( int argc, char *argv[] )
{
  /* lal function variables */
  LALStatus             status = blank_status;

  /* template bank generation parameters */
  CHAR   *bankFileName = NULL;
  INT4    numOutBanks = 0;
  REAL4   minMatch = -1;

  /* output data */
  MetadataTable         inputBank;
  MetadataTable         outputBank;
  MetadataTable         proctable;
  MetadataTable         procparams;
  ProcessParamsTable   *this_proc_param = NULL;
  LIGOLwXMLStream       xmlStream;

  /* counters and other variables */
  INT4 i, j;
  INT4 numTmplts = 0;
  INT4 numTmpltsWritten = 0;
  INT4 numPerFile = 0;
  CHAR *gpsHyphen;
  char outBankFileName[FILENAME_MAX];
  CHAR bankFileNameHead[FILENAME_MAX];
  CHAR bankFileNameTail[FILENAME_MAX];
  CHAR comment[LIGOMETA_COMMENT_MAX];  
  CHAR *userTag = NULL;
  SnglInspiralTable *thisTmplt = NULL;
  SnglInspiralTable *tmpTmplt = NULL;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"version",                 no_argument,       0,                'V'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {"comment",                 required_argument, 0,                's'},    
    {"help",                    no_argument,       0,                'h'}, 
    {"bank-file",               required_argument, 0,                'v'},
    {"number-of-banks",         required_argument, 0,                'n'},
    {"minimal-match",           required_argument, 0,                'M'},
    {0, 0, 0, 0}
  };
  int c;


  /*
   * 
   * initialize things
   *
   */


  lal_errhandler = LAL_ERR_EXIT;
  setvbuf( stdout, NULL, _IONBF, 0 );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  XLALGPSTimeNow(&(proctable.processTable->start_time));
  XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME, LALAPPS_VCS_IDENT_ID,
      LALAPPS_VCS_IDENT_STATUS, LALAPPS_VCS_IDENT_DATE, 0);
  this_proc_param = procparams.processParamsTable = 
    (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );


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
        "i:n:VZ:hs:M:", 
        long_options, &option_index );

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

      case 'v':
        optarg_len = strlen( optarg ) + 1;
        bankFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( bankFileName, optarg, optarg_len );
        snprintf( procparams.processParamsTable->program, 
            LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
        snprintf( procparams.processParamsTable->type, 
            LIGOMETA_TYPE_MAX, "string" );
        snprintf( procparams.processParamsTable->param, 
            LIGOMETA_PARAM_MAX, "--%s", long_options[option_index].name );
        snprintf( procparams.processParamsTable->value, 
            LIGOMETA_VALUE_MAX, "%s", optarg );
        break;

      case 'n':
        numOutBanks = (INT4) atoi( optarg );
        if ( numOutBanks < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "Number of output banks must be greater than zero:" 
              "(%d specified)\n",
              long_options[option_index].name, numOutBanks );
          exit( 1 );
        }
        else if ( numOutBanks > 99 )
        {
          fprintf( stderr, 
              "Warning: generating more than 99 banks is not reccomended!\n" );
        }
        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
            "%s", PROGRAM_NAME );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "int" );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
            "--%s", long_options[option_index].name );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%d", 
            numOutBanks );
        break;

      case 's':
        if ( strlen( optarg ) > LIGOMETA_COMMENT_MAX - 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "comment must be less than %d characters\n",
              long_options[option_index].name, LIGOMETA_COMMENT_MAX );
          exit( 1 );
        }
        else
        {
          snprintf( comment, LIGOMETA_COMMENT_MAX, "%s", optarg );
        }
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
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "-userTag" );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            optarg );
        break;

      case 'M':
        minMatch = (REAL4) atof( optarg );
        if ( minMatch <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "minimal match of bank must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, minMatch );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
            PROGRAM_NAME );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "float" );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s",
            long_options[option_index].name );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%e",
            minMatch );
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "Inspiral Template Bank Splitter\n" 
            "Duncan Brown <duncan@gravity.phys.uwm.edu>\n");
        XLALOutputVersionString(stderr, 0);
        exit( 0 );
        break;

      case '?':
        fprintf( stderr, USAGE, argv[0] );
        exit( 1 );
        break;

      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        fprintf( stderr, USAGE, argv[0] );
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

  /* check the values of the arguments */
  if ( ! bankFileName )
  {
    fprintf( stderr, "Error: --bank-file must be specified\n" );
    exit( 1 );
  }

  if ( ! numOutBanks )
  {
    fprintf( stderr, "Error: --number-of-banks must be specified\n" );
    exit( 1 );
  }

  if ( minMatch < 0 )
  {
    fprintf( stderr, "Error: --minimal-match must be specified\n" );
    exit( 1 );
  }


  /*
   *
   * read in the template bank from the input file
   *
   */


  /* read in the template bank from a ligo lw xml file */
  inputBank.snglInspiralTable = NULL;
  numTmplts = LALSnglInspiralTableFromLIGOLw( &(inputBank.snglInspiralTable), 
      bankFileName, 0, -1 );
  if ( numTmplts < 0 )
  {
    fprintf( stderr, "error: unable to read templates from %s\n", 
        bankFileName );
    exit( 1 );
  }
  
  if ( vrbflg ) fprintf( stdout, "read %d templates from %s\n", 
      numTmplts, bankFileName );

  /* find the hypen just before the GPS start time of the bank */
  gpsHyphen = NULL;
  gpsHyphen = strstr( bankFileName, "-" );
  if ( ! gpsHyphen )
  {
    fprintf( stderr, "Error: could not find first hypen in file name %s\n",
        bankFileName );
    exit( 1 );
  }
  gpsHyphen = strstr( gpsHyphen + 1, "-" );
  if ( ! gpsHyphen )
  {
    fprintf( stderr, "Error: could not find second hypen in file name %s\n",
        bankFileName );
    exit( 1 );
  }

  /* store the name of the template bank file */
  memcpy( bankFileNameHead, bankFileName, 
      (size_t) gpsHyphen - (size_t) bankFileName < FILENAME_MAX ? 
      (gpsHyphen - bankFileName) * sizeof(CHAR) : FILENAME_MAX * sizeof(CHAR) );
  strncpy( bankFileNameTail, gpsHyphen + 1, FILENAME_MAX );

  if ( vrbflg )
  {
    fprintf( stdout, "head of bank file name is %s\n", bankFileNameHead );
    fprintf( stdout, "tail of bank file name is %s\n", bankFileNameTail );
  }

  
  /*
   *
   * write out the individual tempate bank files
   *
   */
  

  /* compute the number of templates per output file */
  numPerFile = floor( ( numTmplts - 0.5 )/ numOutBanks + 1 );
  thisTmplt = inputBank.snglInspiralTable;
  if ( vrbflg ) fprintf( stdout, "writing around %d templates per file\n", 
      numPerFile );

  for ( i = 0; i < numOutBanks; ++i )
  {
    /* open the output xml file */
    memset( outBankFileName, 0, FILENAME_MAX * sizeof(CHAR) );
    snprintf( outBankFileName, FILENAME_MAX, "%s_%2.2d-%s", 
        bankFileNameHead, i, bankFileNameTail );
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
    LAL_CALL( LALOpenLIGOLwXMLFile( &status , &xmlStream, outBankFileName), 
        &status );

    if ( vrbflg ) 
      fprintf( stdout, "writing templates to %s... ", outBankFileName );

    /* write process table */
    XLALGPSTimeNow(&(proctable.processTable->end_time));
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, proctable, 
          process_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
    
    /* write process_params table */
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
          process_params_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, procparams, 
          process_params_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

    /* write the templates to the file */
    outputBank.snglInspiralTable = thisTmplt;
    numTmpltsWritten = 0;

    if ( thisTmplt )
    {
      LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream, 
            sngl_inspiral_table), &status );

      for ( j = 0; j < numPerFile - 1 && thisTmplt->next; ++j )
      {
        thisTmplt = thisTmplt->next;
      }
      tmpTmplt = thisTmplt->next;
      thisTmplt->next = NULL;
      thisTmplt = tmpTmplt;

      LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputBank,
            sngl_inspiral_table), &status );
      LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream), &status );
    }

    while ( outputBank.snglInspiralTable )
    {
      ++numTmpltsWritten;
      tmpTmplt = outputBank.snglInspiralTable;
      outputBank.snglInspiralTable = outputBank.snglInspiralTable->next;
      LALFree( tmpTmplt );
    }

    LAL_CALL( LALCloseLIGOLwXMLFile( &status, &xmlStream), &status );

    if ( vrbflg ) fprintf( stdout, "%d templates\n", numTmpltsWritten );
  }

  LALCheckMemoryLeaks();
  exit( 0 );
}
