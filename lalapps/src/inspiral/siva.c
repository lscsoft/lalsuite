/*--------------------------------------------------------------------------
 * File name: siva.c  -  Simple Inspiral Veto Application
 * Author: Peter Shawhan (borrowing heavily from sire.c)
 * Revision: $Id $
 *--------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/SegmentsIO.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXMLRead.h>
#include <lalapps.h>

RCSID("$Id$");

#define USAGE \
  "Usage: lalapps_siva INPUTFILE VETOFILE OUTPUTFILE\n"\
"\n"\
"  INPUTFILE    XML file of inspiral triggers to read\n"\
"  VETOFILE     Segment list file of time intervals to be vetoed\n"\
"  OUTPUTFILE   Name of XML file to create with surviving inspiral triggers\n"\
"\n"\
"This is the Simple Inspiral Veto Application.  It is Simple in that it\n"\
"only handles one input trigger file and one input veto file, and it\n"\
"completely ignores all LIGO_LW tables in the input file except for the\n"\
"sngl_inspiral table.  The output XML file contains only a sngl_inspiral table.\n"

int main( int argc, char *argv[] )
{
  /* lal initialization variables */
  LALStatus stat = blank_status;

  /* Argument pointers */
  CHAR *inFile = NULL;
  CHAR *vetoFile = NULL;
  CHAR *outFile = NULL;

  /* Program operation variables */
  LALSegList vetoSegs;
  SnglInspiralTable *events = NULL;
  INT4 numEvents = 0;
  INT4 numEventsKept = 0;
  SnglInspiralTable    *eventHead = NULL;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *tmpEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;

  LIGOLwXMLStream       xmlStream;
  MetadataTable         outputTable;

  /*------ Beginning of code ------*/

  /*-- Check command-line arguments --*/
  if ( argc != 4 ) {
    printf( USAGE );
    exit( 0 );
  }
  inFile = argv[1];
  vetoFile = argv[2];
  outFile = argv[3];

  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "33" );


  /*-- Initialize the veto segment list, and read the veto file --*/

  XLALSegListInit( &vetoSegs );
  LAL_CALL( LALSegListRead( &stat, &vetoSegs, vetoFile, "" ), &stat );
  /*-- Make sure the list of veto segments is coalesced for fast searching --*/
  XLALSegListCoalesce( &vetoSegs );


  /*-- Read the inspiral events from the file --*/

  numEvents = LALSnglInspiralTableFromLIGOLw( &events, inFile, 0, -1 );
  if ( numEvents < 0 ) {
    fprintf( stderr, "error: unable to read sngl_inspiral table from %s\n", 
	     inFile );
    exit( 1 );
  }

  /*-- Report the number of events read in --*/
  printf( "Read %d events from input file\n", numEvents );


  /*-- Loop over inspiral triggers --*/

  thisEvent = events;
  while ( thisEvent ) {

    /*-- Check the time of this event against the veto segment list --*/

    if ( XLALSegListSearch( &vetoSegs, &(thisEvent->end_time) ) == NULL ) {
      /* This inspiral trigger does not fall within any veto segment */

      /* keep the trigger and increment the count of triggers */
      if ( ! eventHead ) eventHead = thisEvent;
      prevEvent = thisEvent;
      thisEvent = thisEvent->next;
      ++numEventsKept;

    } else {
      /*-- This event's end_time falls within one of the veto segments --*/

      /* discard the trigger and move to the next one */
      if ( prevEvent ) prevEvent->next = thisEvent->next;
      tmpEvent = thisEvent;
      thisEvent = thisEvent->next;
      LAL_CALL ( LALFreeSnglInspiral ( &stat, &tmpEvent ), &stat);

    }

  }

  /* make sure that the linked list is properly terminated */
  if ( prevEvent && prevEvent->next ) prevEvent->next->next = NULL;

  /*-- Report the number of events kept --*/
  printf( "Kept %d events\n", numEventsKept );


  /*-- Write out the surviving triggers --*/

  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &stat, &xmlStream, outFile ), &stat );

  if ( eventHead )
  {
    outputTable.snglInspiralTable = eventHead;
    LAL_CALL( LALBeginLIGOLwXMLTable( &stat, &xmlStream, 
          sngl_inspiral_table ), &stat );
    LAL_CALL( LALWriteLIGOLwXMLTable( &stat, &xmlStream, outputTable, 
          sngl_inspiral_table ), &stat );
    LAL_CALL( LALEndLIGOLwXMLTable( &stat, &xmlStream ), &stat);
  }

  /* close the output file */
  LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);


  /*-- Free memory and exit --*/

  while ( eventHead )
  {
    thisEvent = eventHead;
    eventHead = eventHead->next;
    LAL_CALL ( LALFreeSnglInspiral ( &stat, &thisEvent ), &stat);
  }

  LALCheckMemoryLeaks();
  exit( 0 );
}
