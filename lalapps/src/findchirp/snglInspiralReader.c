#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOLwXML.h>
#include <lalapps.h>
#include "event_utils.h"
RCSID("$Id$");

#define MAXSTR 2048

/* Usage format string. */
#define USAGE "Usage: %s --input infile --table tablename --outfile filename \
    [--snrstar snrstar] [--noplayground] [--help]\n"

#define SNGLINSPIRALREADER_EARG   1
#define SNGLINSPIRALREADER_EROW   2
#define SNGLINSPIRALREADER_EFILE  3

#define SNGLINSPIRALREADER_MSGEARG   "Error parsing arguments"
#define SNGLINSPIRALREADER_MSGROW    "Error reading row from XML table"
#define SNGLINSPIRALREADER_MSGEFILE  "Could not open file"

#define TRUE  1
#define FALSE 0

static int getline(char *line, int max, FILE *fpin)
{
    int i;
    CHAR tmpline[MAXSTR];

    for (i=0;i<MAXSTR;i++) { *(line+i) = '\0' ; }
    if (fgets(tmpline, max, fpin) == NULL){
        return 0;
    }
    else{
        strncpy(line,tmpline,strlen(tmpline)-1);
        return strlen(line);
    }
}


/****************************************************************************
 * 
 * FUNCTION TESTS IF THE FILE CONTAINS ANY PLAYGROUND DATA
 * 
 ***************************************************************************/
static int isPlayground(INT4 gpsStart, INT4 gpsEnd){
    INT4 runStart=729273613;
    INT4 playInterval=6370;
    INT4 playLength=600;
    INT4 segStart,segEnd,segMiddle;

    segStart = (gpsStart - runStart)%playInterval;
    segEnd   = (gpsEnd - runStart)%playInterval;
    segMiddle = gpsStart + (INT4) (0.5 * (gpsEnd - gpsStart));
    segMiddle = (segMiddle - runStart)%playInterval;
    
    if (segStart < 600 || segEnd < 600 || segMiddle < 600){
        return TRUE;
    }

    return FALSE;
}

/************************************************************************************
 *
 * The main program
 *
 ***********************************************************************************/

int main(int argc, char **argv)
{
    static LALStatus         stat;
    INT4                     retVal=0;
    INT4                     inarg=1;
    FILE                     *fpin=NULL;
    INT4                     fileCounter=0;
    BOOLEAN                  playground=TRUE;

    REAL4                    snrstar;
    BOOLEAN                  checkMasses;
    REAL4                    mass1,mass2,dm;

    CHAR                     inputfile[MAXSTR],line[MAXSTR];

    CHAR                     *tablename=NULL;
    CHAR                     *outfileName=NULL;
    const CHAR               *searchSummaryName="search_summary";
    SearchSummaryIndex        searchSummaryIndex;
    SearchSummaryTable        searchSummaryTable;
    SnglInspiralIndex         tableIndex;
    SnglInspiralTable        *currentEvent=NULL,*prevEvent=NULL;
    SnglInspiralTable         inspiralEvent,*inspiralEventList=NULL;
    MetadataTable             myTable;
    LIGOLwXMLStream           xmlStream;
    INT4                      i;

    struct MetaioParseEnvironment triggerParseEnv;
    const MetaioParseEnv triggerEnv = &triggerParseEnv;

    /*******************************************************************
     * BEGIN PARSE ARGUMENTS (inarg stores the current position)        *
     *******************************************************************/
    if (argc <= 1){
        LALPrintError( USAGE, *argv );
        return 0;
    }

    while ( inarg < argc ) {
        /* file containing a list of xml files to use */
        if ( !strcmp( argv[inarg], "--input" ) )
        {
            inarg++;
            sprintf(inputfile,"%s",argv[inarg++]);
        }
        /* table to read from xml file */
        else if ( !strcmp( argv[inarg], "--table" ) )
        {
            inarg++;
            tablename = argv[inarg++];
        }
        /* output file name  */
        else if ( !strcmp( argv[inarg], "--outfile" ) )
        {
            inarg++;
            outfileName = argv[inarg++];
        }
        /* return all events.  By default,  returns only playground */
        else if ( !strcmp( argv[inarg], "--noplayground" ) )
        {
            inarg++;
            playground=FALSE;;
        }
        /* select events with SNR > snrstar */
        else if ( !strcmp( argv[inarg], "--snrstar" ) ) {
            inarg++;
            snrstar = atof(argv[inarg++]);
        }
        /* specify template masses & error */
        else if ( !strcmp( argv[inarg], "--template" ) ) {
            checkMasses = TRUE;
            inarg++;
            if (argv[inarg] && (argv[inarg][0] != '-'))
            {
                mass1 = atof(argv[inarg++]);
                if (argv[inarg] && (argv[inarg][0] != '-'))
                {
                    mass2 = atof(argv[inarg++]);
                    if (argv[inarg] && (argv[inarg][0] != '-'))
                    {
                        dm = atof(argv[inarg++]);
                    }
                }
            }
        }
        /* print a help message */
        else if ( !strcmp( argv[inarg], "--help" ) ) {
            LALPrintError( USAGE, *argv );
            return SNGLINSPIRALREADER_EARG;
        }
        /* Check for unrecognized options. */
        else if ( argv[inarg][0] == '-' ) {
            LALPrintError( "Unrecognized options\n" );
            return SNGLINSPIRALREADER_EARG;
        }
        /* Default case for unknown arguments */
        else
        {
            LALPrintError( "Unknown arguments\n" );
            return SNGLINSPIRALREADER_EARG;
        }
    }
    if (inputfile == NULL){
        LALPrintError( "Must supply an xml file to parse\n" );
        return SNGLINSPIRALREADER_EARG;
    }

    /*******************************************************************
    * END PARSE ARGUMENTS                                              *
    *******************************************************************/


    /*******************************************************************
    * initialize things
    *******************************************************************/
    lal_errhandler = LAL_ERR_EXIT;
    set_debug_level( "3" );
    memset( &inspiralEvent, 0, sizeof(inspiralEvent) );
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );

    
    /*****************************************************************
     * OPEN FILE WITH LIST OF XML FILES (one file per line)
     ****************************************************************/
    if ( !(fpin = fopen(inputfile,"r")) ){
        LALPrintError("Could not open input file\n");
    }

    
    /*****************************************************************
     * open output xml file
     *****************************************************************/
    LAL_CALL( LALOpenLIGOLwXMLFile(&stat, &xmlStream, outfileName), &stat);
    LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, sngl_inspiral_table), &stat);

    
    /*****************************************************************
     * loop over the xml files
     *****************************************************************/
    while ( getline(line, MAXSTR, fpin) ){
        fileCounter++;

        /* open xml file at search summary table */
        if ( (retVal = MetaioOpenTable( triggerEnv, line, searchSummaryName)) !=0 ){
            fprintf(stderr, "Error opening injection file %s for %s\n", 
                    line, searchSummaryName );
            MetaioAbort( triggerEnv ); 
            exit(SNGLINSPIRALREADER_EFILE);
        }

        /* get the search summary table */
        LAL_CALL( buildSearchSummaryIndex(&stat, triggerEnv, &searchSummaryIndex), 
                &stat);
        LAL_CALL( getSearchSummaryTable(&stat, triggerEnv, &searchSummaryTable, 
                    &searchSummaryIndex), &stat);

        /* check for events and playground */
        if ( ( playground && 
                    !(isPlayground(searchSummaryTable.out_start_time.gpsSeconds,
                           searchSummaryTable.out_end_time.gpsSeconds ))) ){
            fprintf(stdout,"File %i %s: not in playground, continuing\n",
                    fileCounter,line);
            continue;
        }
        else if ( searchSummaryTable.nevents == 0 ){
            fprintf(stdout,"File %i %s: no events, continuing\n",
                    fileCounter,line);
            continue;
        }
        else {
            fprintf(stdout,"File %i %s: processing\n",
                    fileCounter,line);
        }

        /* close the stream */
        MetaioAbort( triggerEnv );
        
        /* open xml file at inspiral table */
        if ( (retVal = MetaioOpenTable( triggerEnv, line, tablename)) !=0 ){
            fprintf(stderr, "Error opening injection file %s\n", line );
            MetaioAbort( triggerEnv ); 
            exit(SNGLINSPIRALREADER_EFILE);
        }

        /* Locate the relevant columns in the inspiral table */
        buildSnglInspiralIndex(&stat, triggerEnv, &tableIndex);

        /* Loop over the triggers */
        while (1) {

            /* get the next row */
            retVal = MetaioGetRow(triggerEnv);
            if ( retVal == -1 ) {
                printf( "Error while getting row from file %s\n",line);
                MetaioAbort( triggerEnv ); 
                return SNGLINSPIRALREADER_EROW;
            } else if ( retVal == 0 ) {
                /*-- Reached end of file --*/
                break;
            }

            /* get the inspiral event */
            getSnglInspiralEvent(&stat, triggerEnv, &inspiralEvent, &tableIndex);

            /* check for events and playground */
            if ( ( playground ) && 
                        (inspiralEvent.end_time.gpsSeconds-729273613)%6370 >600 ){
                continue;
            }

            /* check that event satisfies threshold */
            if (inspiralEvent.snr < snrstar){
                continue;
            }
            
            /* allocate memory for the inspiral event */
            if ( (inspiralEventList) == NULL )
            {
                currentEvent=inspiralEventList=
                    (SnglInspiralTable *) LALMalloc( sizeof(SnglInspiralTable) );
            }
            else 
            {
                currentEvent = 
                    (SnglInspiralTable *) LALMalloc( sizeof(SnglInspiralTable) );
            }

            /* copy event into linked list element */
            memcpy(currentEvent, &inspiralEvent, sizeof(SnglInspiralTable) );
            
            /* point to the next event */
            currentEvent->next = NULL;
            if (prevEvent != NULL) prevEvent->next = currentEvent;
            prevEvent = currentEvent;
            currentEvent = currentEvent->next;
        }

        MetaioAbort(triggerEnv);

        /* Write the results to the burst table */
        myTable.snglInspiralTable = inspiralEventList;
        LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, myTable, 
                    sngl_inspiral_table), &stat);

        /* free the temporary memory containing the events */
        while (inspiralEventList)
        {
            SnglInspiralTable *thisEvent;
            thisEvent = inspiralEventList;
            inspiralEventList = inspiralEventList->next;
            LALFree( thisEvent );
        }

        /* clean up for next file */
        currentEvent=prevEvent=inspiralEventList=NULL;
        memset(&inspiralEvent, 0, sizeof(inspiralEvent));
    }

    /* close the output file */
    LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);
    LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);
    LALCheckMemoryLeaks();

    return 0;
}
