/*----------------------------------------------------------------------- 
 * 
 * File Name: trigger2tmplt.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <lal/LIGOLwXML.h>
#include <lalapps.h>
#include "ligolwbank.h"

RCSID( "$Id$" );

#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"


/*
 *
 * compare function used by qsort
 *
 */


int compareTmplts ( const void *a, const void *b )
{
  SnglInspiralTable *aPtr, *bPtr;

  aPtr = (SnglInspiralTable *) a;
  bPtr = (SnglInspiralTable *) b;

  if ( aPtr->mass1 > bPtr->mass1 )
  {
    return 1;
  }
  else if ( aPtr->mass1 < bPtr->mass1 )
  {
    return -1;
  }
  if ( aPtr->mass2 > bPtr->mass2 )
  {
    return 1;
  }
  else if ( aPtr->mass2 < bPtr->mass2 )
  {
    return -1;
  }
  else
  {
    return 0;
  }
}


/*
 *
 * variables that control program behaviour
 *
 */



int main ( int argc, char *argv[] )
{
  extern int vrbflg;            /* verbocity of lal function    */
  char *inputFileName = NULL;   /* input file name              */
  char *outputFileName = NULL;  /* output file name             */
  struct option long_options[] =
  {
    /* these options set a flag */
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"input-file",              required_argument, 0,                'i'},
    {"output-file",             required_argument, 0,                'o'},
    {"debug-level",             required_argument, 0,                'z'},
    {0, 0, 0, 0}
  };
  int i;
  int numEvents;
  REAL4 epsilon = 1e-6;
  SnglInspiralTable    *eventHead = NULL;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;
  SnglInspiralTable   **eventHandle = NULL;
  MetadataTable         outputTable;
  LIGOLwXMLStream       results;
  LALStatus             status = blank_status;
  

  /*
   *
   * parse the command line arguments
   *
   */


  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t namelen;

    i = getopt_long_only( argc, argv, "i:o:z:", long_options, &option_index );

    /* detect the end of the options */
    if ( i == - 1 )
    {
      break;
    }

    switch ( i )
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

      case 'i':
        {
          /* create storage for the input file name name */
          namelen = strlen( optarg ) + 1;
          inputFileName = (CHAR *) LALCalloc( namelen, sizeof(CHAR));
          memcpy( inputFileName, optarg, namelen );
        }
        break;

      case 'o':
        {
          /* create storage for the output file name name */
          namelen = strlen( optarg ) + 1;
          outputFileName = (CHAR *) LALCalloc( namelen, sizeof(CHAR));
          memcpy( outputFileName, optarg, namelen );
        }
        break;
        
      case 'z':
        set_debug_level( optarg );
        break;
      case '?':
        exit( 1 );
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

  
  /* 
   *
   * check that the input and output file names have been specified 
   *
   */


  if ( ! inputFileName )
  {
    fprintf( stderr, "--input-file must be specified\n" );
    exit( 1 );
  }
  if ( ! outputFileName )
  {
    fprintf( stderr, "--output-file must be specified\n" );
    exit( 1 );
  }


  /*
   *
   * read the entire sngl_inspiral table from the input file
   *
   */


  numEvents = SnglInspiralTableFromLIGOLw( &eventHead, inputFileName, 0, 0 );
  if ( numEvents < 0 )
  {
    fprintf( stderr, "error: unable to read sngl_inspiral table from %s\n", 
        inputFileName );
    exit( 1 );
  }
  else if ( numEvents == 0 )
  {
    fprintf( stdout, "no sngl_iinspiral events found in file: %s\n"
        "exiting cleanly\n" , inputFileName );
    LALFree( inputFileName );
    LALFree( outputFileName );
    exit( 0 );
  }

  if ( vrbflg ) fprintf( stdout, "parsed %d sngl_inspiral rows from %s\n", 
      numEvents, inputFileName );

  
  /*
   *
   * order the templates and discard duplicates
   *
   */


  /* create an array of ptrs and qsort it */
  eventHandle = (SnglInspiralTable **) 
    LALCalloc( numEvents, sizeof(SnglInspiralTable *) );

  for ( i = 0, thisEvent = eventHead; i < numEvents; 
      ++i, thisEvent = thisEvent->next )
  {
    eventHandle[i] = thisEvent;
  }

  if ( vrbflg ) fprintf( stdout, "sorting events... " );
  qsort( eventHandle, numEvents, sizeof(SnglInspiralTable *), compareTmplts );
  if ( vrbflg ) fprintf( stdout, "done\n" );

  /* remove duplicate templates from the sorted list */
  prevEvent = eventHead;
  for ( i = 1; i < numEvents; ++i )
  {
    thisEvent = eventHandle[i];
    if ( (fabs( prevEvent->mass1 - thisEvent->mass1 ) < epsilon) &&
        (fabs( prevEvent->mass2 - thisEvent->mass2 ) < epsilon) )
    {
      /* remove this event from the linked list */
      if ( vrbflg ) 
        fprintf( stdout, "removing event %d from linked list\n", i );
      prevEvent->next = thisEvent->next;
      LALFree( thisEvent );
      thisEvent = prevEvent->next;
    }
    prevEvent = thisEvent;
  }


  /*
   *
   * write the sorted, unique list out to an XML file
   *
   */


  outputTable.snglInspiralTable = eventHead;
  
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &results, outputFileName),
      &status );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, sngl_inspiral_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, outputTable, 
        sngl_inspiral_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &results ), &status );


  /*
   *
   * free memory and exit
   *
   */


  LALFree( eventHandle );
  while ( eventHead )
  {
    thisEvent = eventHead;
    eventHead = eventHead->next;
    LALFree( eventHead );
  }
  LALFree( inputFileName );
  LALFree( outputFileName );

  LALCheckMemoryLeaks();
  
  return 0;
}
