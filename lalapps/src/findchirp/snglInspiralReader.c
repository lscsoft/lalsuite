#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lalapps.h>
#include <lal/LIGOMetadataUtils.h>
#include "event_utils.h"
RCSID("$Id$");

#define MAXSTR 2048

/* Usage format string. */
#define USAGE "Usage: %s --input infile --table tablename --outfile filename \
    [--snrstar snrstar] [--noplayground] [--sort] [--cluster msec] \
    [--clusteralgorithm clusterchoice] [--help]\n"

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

/****************************************************************************
 *
 * The main program
 *
 *****************************************************************************/

int main(int argc, char **argv)
{
    static LALStatus         stat;
    INT4                     retVal=0;
    FILE                     *fpin=NULL;
    INT4                     fileCounter=0;
    BOOLEAN                  playground=TRUE;
    INT4                     sort=FALSE;
    INT4                     cluster=FALSE;   
    Clusterchoice	     clusterchoice;
	
    REAL4                    snrstar=0.0;
    INT4                     dtime=0;

    CHAR                     inputfile[MAXSTR],line[MAXSTR];

    CHAR                     *tablename=NULL;
    CHAR                     *outfileName=NULL;
    const CHAR               *searchSummaryName="search_summary";
    const CHAR               *snglInspiralName="sngl_inspiral";
    SearchSummaryIndex        searchSummaryIndex;
    SearchSummaryTable        searchSummaryTable;
    SnglInspiralIndex         tableIndex;
    SnglInspiralTable        *currentEvent=NULL,*prevEvent=NULL;
    SnglInspiralTable         inspiralEvent,*inspiralEventList=NULL;
    MetadataTable             myTable;
    LIGOLwXMLStream           xmlStream;
    INT4                      numEvents;

    static int		      verbose_flag = 0;

    struct MetaioParseEnvironment triggerParseEnv;
    const MetaioParseEnv triggerEnv = &triggerParseEnv;

    /*******************************************************************
     * BEGIN PARSE ARGUMENTS (inarg stores the current position)        *
     *******************************************************************/
   int c;
    while (1)
    {
	/* getopt arguments */
	static struct option long_options[] = 
	{
	    /* these options set a flag */
	    {"verbose",         no_argument,	&verbose_flag, 1 },
 	    /* parameters which determine the output xml file */
	    {"input",		required_argument,  0,	'a'},
	    {"table",		required_argument,  0,	'b'},
	    {"outfile",		required_argument,  0,	'c'},
	    {"snrstar",		required_argument,  0,	'd'},
	    {"cluster",		required_argument,  0,	'e'},
	    {"clusteralgorithm",required_argument,  0,	'f'},
	    {"noplayground",	required_argument,  0,	'g'},
	    {"help",		no_argument,	    0,	'h'}, 
	    {"sort",		no_argument,	    0,	'i'},
	    {0, 0, 0, 0}
	};
	/* getopt_long stores the option index here. */
	int option_index = 0;
 
    	c = getopt_long ( argc, argv, "a:b:c:d:e:f:g:h:i:", 
				long_options, &option_index );

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

	case 'a':
	    /* file containing list of xml files to use */
	    {
		sprintf(inputfile,"%s",optarg);
    	    }
	    break;	
	
	case 'b':
	    /* table to be read from xml file */
	    {
		tablename = optarg;
	    }
	    break;
	
	case 'c':
	    /* output file name */
	    {
		outfileName = optarg;
	    }
	    break;

	case 'd':
	    /* a threshold cut on signal to noise, only triggers with
	    * SNR > snrstar are kept */
	    {
		snrstar = atof( optarg );
	    }
	    break;

	case 'e':
	    /* cluster the events */ 
	    {
		dtime = atoi( optarg );
		cluster = TRUE;
	    }
	    break;
    
	case 'f':
	    /* choose the clustering algorithm, default is 0
	     * options are detailed in event_utils.c */
	    {
		if ( ! strcmp( "snr_and_chisq", optarg ) )
	        {
		    clusterchoice = snr_and_chisq;
		}
		else if ( ! strcmp( "snrsq_over_chisq", optarg) )
		{
		    clusterchoice = snrsq_over_chisq;
		}	    
		else
		{
		    fprintf( stderr, "invalid argument to  --%s:\n"
			"unknown clustering specified:\n "
			"%s (must be one of: snr_and_chisq, snrsq_over_chisq)\n",
			long_options[option_index].name, optarg);
		    exit( 1 );
		}
	    }
	    break;
   
	case 'g':
	    /* don't restrict to the playground data */
	    {
		playground = FALSE;
	    }
	    break;
    
	case 'h':
	    /* print help */
	    {
		LALPrintError( USAGE, *argv );
		return SNGLINSPIRALREADER_EARG;
	    }
	    break;

	    
	case 'i':
	    /* sort the events in time */
	    {
		sort = TRUE;
	    }
	    break;

	    
	default:
	    {
		return SNGLINSPIRALREADER_EARG;
	    }
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

    if (inputfile == NULL)
    {
        LALPrintError( "Must supply an xml file to parse\n" );
        return SNGLINSPIRALREADER_EARG;
    }

    if ( ! tablename )
    { 
	LALPrintError( "Table name must be specified\n" );
	return SNGLINSPIRALREADER_EARG;
    }

    if ( ! outfileName )
    {
    	LALPrintError( "Outfile name must be specified\n" );
	return SNGLINSPIRALREADER_EARG;
    }


    /*******************************************************************
    * END PARSE ARGUMENTS                                              *
    *******************************************************************/


    /*******************************************************************
    * initialize things
    *******************************************************************/
    lal_errhandler = LAL_ERR_EXIT;
    set_debug_level( "1" );
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
    LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, sngl_inspiral_table), 
		    &stat);

    
    /*****************************************************************
     * loop over the xml files
     *****************************************************************/
    while ( getline(line, MAXSTR, fpin) ){
        fileCounter++;

        /* open xml file at search summary table */
        if ( (retVal = MetaioOpenTable( triggerEnv, line, searchSummaryName))
		       	!=0 ){
            fprintf(stderr, "Error opening injection file %s for %s\n", 
                    line, searchSummaryName );
            if (fileCounter>1){
            MetaioAbort( triggerEnv ); 
            exit(SNGLINSPIRALREADER_EFILE);
            } else {
                fprintf(stderr,"Proceeding anyway\n");
            }
        }

        if ( !retVal ){
        /* get the search summary table */
        LAL_CALL( buildSearchSummaryIndex(&stat, triggerEnv, 
		    &searchSummaryIndex), &stat);
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
        }
        
        /* open xml file at inspiral table */
        if ( (retVal = MetaioOpenTable( triggerEnv, line, tablename)) !=0 ){
            fprintf(stderr, "Error opening injection file %s\n", line );
            MetaioAbort( triggerEnv ); 
            exit(SNGLINSPIRALREADER_EFILE);
        }

        /* Locate the relevant columns in the inspiral table */
        buildSnglInspiralIndex(&stat, triggerEnv, &tableIndex);

        numEvents=0;
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
            getSnglInspiralEvent(&stat, triggerEnv, &inspiralEvent, 
			    &tableIndex);

            /* check for events and playground */
            if ( ( playground ) && 
                        (inspiralEvent.end_time.gpsSeconds-729273613)%6370 
			    >600 ){
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

            numEvents++;
        }

        MetaioAbort(triggerEnv);

        /* sort the events */
        if (sort){
            LAL_CALL( LALSortSnglInspiral( &stat, &inspiralEventList,
                        *LALCompareSnglInspiralByTime), &stat);
        }

        /* cluster the events */
        if ( cluster && inspiralEventList){
            LAL_CALL( LALClusterSnglInspiralTable( &stat, inspiralEventList,
                        dtime, clusterchoice), &stat);
        }

        /* Write the results to the inspiral table */
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
