/*----------------------------------------------------------------------- 
 * 
 * File Name: inspiral.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>

#include <lalapps.h>
#include <processtable.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>
#include <lal/Date.h>

#include "inspiral.h"

#define RESULT_FILE "results.xml"

RCSID( "$Id$" );

#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"

/* global debugging options */
extern int vrbflg;                      /* verbocity of lal function    */

int main( int argc, char *argv[] )
{
  LALStatus             status = blank_status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;
  MetadataTable         proctable;
  LIGOLwXMLStream       results;


  /*
   *
   * initialization
   *
   */


  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );

  /* initialize the process table with entries and start time */
  proctable.processTable = (ProcessTable *) LALMalloc( sizeof(ProcessTable) );

  LAL_CALL( 
      populate_process_table( &status, proctable.processTable, 
       PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ),
      &status );
  LAL_CALL( 
      LALGPSTimeNow ( &status, &(proctable.processTable->start_time),
        &accuracy ), &status );


  /*
   *
   * parse command line options and write process_params
   *
   */

  /* inspiral_initialize( argc, argv, procparams ); */


  /*
   *
   * do something
   *
   */


  sleep( 2 );


  /*
   *
   * close the output file, write the result file and exit sucessfully
   *
   */


  LAL_CALL( 
      LALGPSTimeNow ( &status, &(proctable.processTable->end_time),
        &accuracy ), &status );

  sprintf( proctable.processTable->ifos, "L1" );

  memset( &results, 0, sizeof(LIGOLwXMLStream) );

  LAL_CALL( 
      LALOpenLIGOLwXMLFile( &status, &results, RESULT_FILE ), &status );

  LAL_CALL( 
      LALBeginLIGOLwXMLTable( &status, &results, process_table ), &status );

  LAL_CALL( 
      LALWriteLIGOLwXMLTable( &status, &results, 
        proctable, process_table ), &status );

  LAL_CALL( 
      LALEndLIGOLwXMLTable ( &status, &results ), &status );

  LAL_CALL( 
      LALCloseLIGOLwXMLFile ( &status, &results ), &status );

  LALFree( proctable.processTable );

  LALCheckMemoryLeaks();

  /* exit sucessfully */
  exit( 0 );
}
