/*
*  Copyright (C) 2007 Lisa M. Goggin
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
   * File Name: collatesim.c
   *
   * Author: Goggin, L. M. 
   *
   *-----------------------------------------------------------------------
   */

#include <stdio.h>
#include <stdlib.h>
#include <config.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <glob.h>
#include <lalapps.h>
#include <processtable.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXMLRingdownRead.h>
#include <lal/LIGOLwXML.h>
#include <lal/Date.h>

#include <LALAppsVCSInfo.h>

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "collasimm"

#define USAGE \
  "lalapps_rinj [options]\n"\
"\nDefaults are shown in brackets\n\n" \
"  --help                   display this message\n"\
"  --verbose                turn verbose flag on\n"\
"\n"\
"Input data source:\n"\
"  --glob GLOB                  use pattern GLOB to determine the input files\n"\
"  --input FILE                 read list of input XML files from FILE\n"\
"\n"\
"Output data destination:\n"\
"  --output FILE                write output data to FILE\n"

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
calloc( 1, sizeof(ProcessParamsTable) ); \
snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        PROGRAM_NAME ); \
snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
        long_options[option_index].name ); \
snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

#define MAX_PATH 4096

/* function to read the next line of data from the input file list */
char *get_next_line( char *line, size_t size, FILE *fp );
char *get_next_line( char *line, size_t size, FILE *fp )
{
  char *s;
  do
    s = fgets( line, size, fp );
  while ( ( line[0] == '#' || line[0] == '%' ) && s );
  return s;
}

extern int vrbflg;

int main( int argc, char *argv[] )
{
  /* lal initialization variables */
  LALStatus status = blank_status;

  /*  program option variables */
  CHAR comment[LIGOMETA_COMMENT_MAX];
  char *inputGlob = NULL;
  char *inputFileName = NULL;
  char *outputFileName = NULL;
  char line[MAX_PATH];
  FILE *fp = NULL;
  glob_t globbedFiles;
  int numInFiles = 0;
  char **inFileNameList;
  int j;
  int  errnum;

  MetadataTable         proctable;
  MetadataTable         procparams;
  ProcessParamsTable   *this_proc_param;
 
  SimRingdownTable   **eventHandle = NULL;
  SimRingdownTable    *eventHead = NULL;
  SimRingdownTable    *thisEvent = NULL;
  SimRingdownTable    *prevEvent = NULL;

  LIGOLwXMLStream       xmlStream;
  MetadataTable         outputTable;

  UINT4                 numEvents = 0;
  UINT4                 numEventsKept = 0;
  INT8 ta=0.0;
  INT8 tb;
  
 /*
  *
  * initialization
  *
  */
  
  
  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  
  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *)
    calloc( 1, sizeof(ProcessTable) );
  XLALGPSTimeNow(&(proctable.processTable->start_time));

  XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME,
      lalAppsVCSIdentId, lalAppsVCSIdentStatus, lalAppsVCSIdentDate, 0);

  this_proc_param = procparams.processParamsTable =
    (ProcessParamsTable *)
    calloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );
  
    
  /*
   *
   * parse command line arguments
   *
   */
  
  
  while (1)
  {
    /* getopt arguments */
    static struct option long_options[] =
    {
      {"verbose",             no_argument,        &vrbflg,              1 },
      {"help",                    no_argument,            0,           'h'},
      {"comment",                 required_argument,      0,           'c'},
      {"version",                 no_argument,            0,           'V'},
      {"glob",                    required_argument,      0,           'g'},
      {"input",                   required_argument,      0,           'i'},
      {"output",                  required_argument,      0,           'o'}
    };
    int c;

    /* getopt_long stores the option index here. */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only ( argc, argv, "h:g:i:o:V", long_options, &option_index );

    /* detect the end of the options */
    if ( c == - 1 )
      break;

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

      case 'h':
         fprintf( stdout, USAGE );
         exit( 0 );
         break;

      case 'g':
         /* create storage for the input file glob */
         optarg_len = strlen( optarg ) + 1;
         inputGlob = (CHAR *) calloc( optarg_len, sizeof(CHAR));
         memcpy( inputGlob, optarg, optarg_len );
         ADD_PROCESS_PARAM( "string", "'%s'", optarg );
         break;

      case 'i':
         /* create storage for the input file name */
         optarg_len = strlen( optarg ) + 1;
         inputFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
         memcpy( inputFileName, optarg, optarg_len );
         ADD_PROCESS_PARAM( "string", "%s", optarg );
         break;

      case 'o':
         /* create storage for the output file name */
         optarg_len = strlen( optarg ) + 1;
         outputFileName = (CHAR *) calloc( optarg_len,sizeof(CHAR));
         memcpy( outputFileName, optarg, optarg_len );
         ADD_PROCESS_PARAM( "string", "%s", optarg );
         break;

      default:
         fprintf( stderr, "unknown error while parsing options\n" );
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

  if ( inputGlob )
  {
    /* use glob() to get a list of the input file names */
    if ( glob( inputGlob, GLOB_ERR, NULL, &globbedFiles ) )
    {
      fprintf( stderr, "error globbing files from %s\n", inputGlob );
      perror( "error:" );
      exit( 1 );
    }
    
    numInFiles = globbedFiles.gl_pathc;
    inFileNameList = (char **) LALCalloc( numInFiles, sizeof(char *) );
    
    for ( j = 0; j < numInFiles; ++j )
    {
      inFileNameList[j] = globbedFiles.gl_pathv[j];
    }
  }
  else if ( inputFileName )
  {
    /* read the list of input filenames from a file */
    fp = fopen( inputFileName, "r" );
    if ( ! fp )
    {
      fprintf( stderr, "could not open file containing list of xml files\n" );
      perror( "error:" );
      exit( 1 );
    }
    
    /* count the number of lines in the file */
    while ( get_next_line( line, sizeof(line), fp ) )
    {
      ++numInFiles;
    }
    rewind( fp );
    
    /* allocate memory to store the input file names */
    inFileNameList = (char **) LALCalloc( numInFiles, sizeof(char *) );
    
    /* read in the input file names */
    for ( j = 0; j < numInFiles; ++j )
    {
      inFileNameList[j] = (char *) LALCalloc( MAX_PATH, sizeof(char) );
      get_next_line( line, sizeof(line), fp );
      strncpy( inFileNameList[j], line, strlen(line) - 1);
    }
    fclose( fp );
  }
  else
  {
    fprintf( stderr, "no input file mechanism specified\n" );
    exit( 1 );
  }
  
  if ( vrbflg )
  {
    fprintf( stdout, "reading input triggers from:\n" );
    for ( j = 0; j < numInFiles; ++j )
    {
      fprintf( stdout, "%s\n", inFileNameList[j] );
    }
  }
 
  for ( j = 0; j < numInFiles; ++j )  
  {
    if ( ! prevEvent )
    {
        eventHandle = &thisEvent;
    }
    else
    {
      eventHandle = &(prevEvent->next);
    }
  
    /* read the events from the file into a temporary list */

    XLAL_TRY( *eventHandle = XLALSimRingdownTableFromLIGOLw( inFileNameList[j], 0, 0 ), errnum);
    if ( ! *eventHandle )
      switch ( errnum )
      {
        case XLAL_EDATA:
          XLALPrintError("Unable to read sngl_ringdown table from %s\n", 
              inFileNameList[j] );
          /*LALFree(thisInputFile);*/
          XLALClearErrno();
          break;
        default:
          XLALSetErrno( errnum);
          XLAL_ERROR( XLAL_EFUNC );
      }
      thisEvent = *eventHandle;
      numEvents++;
      if ( ! eventHead ) eventHead = thisEvent;

      if ( ! prevEvent )
      {
        ta = XLALGPSToINT8NS( &(thisEvent->h_start_time) );
        prevEvent = thisEvent;
        thisEvent = thisEvent->next;
        ++numEventsKept;
      }
      else
      {
        tb = XLALGPSToINT8NS( &(thisEvent->h_start_time) );
        if( ta!=tb)
        {
          prevEvent = thisEvent;
          thisEvent = thisEvent->next;
          ++numEventsKept;
          ta=tb;
        }
      }
  }
  
      
  /*
   *
   * write output data
   *
   */


  /* write the main output file containing found injections */
  if ( vrbflg ) fprintf( stdout, "writing output xml files... " );
  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlStream, outputFileName ), &status );

  /* write out the process and process params tables */
  if ( vrbflg ) fprintf( stdout, "process... " );
  XLALGPSTimeNow(&(proctable.processTable->end_time));
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ),
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, proctable,
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
  free( proctable.processTable );
  
  /* write the process params table */
  if ( vrbflg ) fprintf( stdout, "process_params... " );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream,
        process_params_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, procparams,
        process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

  /* Write the results to the ringdown table */
  if ( eventHead )
  {
    if ( vrbflg ) fprintf( stdout, "sim_ringdown... " );
    outputTable.simRingdownTable = eventHead;
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream,
          sim_ringdown_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable,
          sim_ringdown_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status);
  }

  /* close the output file */
  LAL_CALL( LALCloseLIGOLwXMLFile(&status, &xmlStream), &status);
  if ( vrbflg ) fprintf( stdout, "done\n" );

  /*
   * 
   *   free memory and exit
   *
   */
  
  
  /* free the ringdown events we saved */
  while ( eventHead )
  {
    thisEvent = eventHead;
    eventHead = eventHead->next;
    LALFree( thisEvent );
  }
}
