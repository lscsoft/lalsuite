#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOMetadataUtils.h>
#include <lalapps.h>
#include "event_utils.h"
RCSID("$Id$");

#define MAXSTR 2048

/* Usage format string. */
#define USAGE \
"Usage: %s --input infile [--inject injectfile] --outfile outfile\n" \
"[--snrstar snrstar] [--sort] [--noplayground] [--deltat dt]\n" \
"[--help]  [--cluster clust]  [--clusteralgorithm clusterchoice]\n" \
"[--missedinjections missedfile]\n"

#define INSPINJFIND_EARG   1
#define INSPINJFIND_EROW   2
#define INSPINJFIND_EFILE  3

#define INSPINJFIND_MSGEARG   "Error parsing arguments"
#define INSPINJFIND_MSGROW    "Error reading row from XML table"
#define INSPINJFIND_MSGEFILE  "Could not open file"

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
 

/*****************************************************************************
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
    INT4		     saveMissedInjections=FALSE;
    INT4		     coincidence = FALSE;
	
    INT8		     simTime;
    INT8		     inspiralTime;
    
    REAL4                    snrstar=0.0;
    INT4                     clust=0;
    INT4		     dt = 20;
    
    CHAR                     inputfile[MAXSTR],line[MAXSTR];

    CHAR                     *outfileName=NULL;
    CHAR		     *injectfile=NULL;
    CHAR		     *missedInjectionsFile=NULL;
    
    const CHAR               *searchSummaryName="search_summary";
    const CHAR		     *simInspiralName="sim_inspiral";
    const CHAR		     *snglInspiralName="sngl_inspiral";


    SearchSummaryIndex        searchSummaryIndex;
    SearchSummaryTable        searchSummaryTable;

    SnglInspiralIndex         tableIndex;
    SnglInspiralTable        *currentEvent=NULL,*prevEvent=NULL;
    SnglInspiralTable         inspiralEvent,*inspiralEventList=NULL;
    SnglInspiralTable	    *currentInspiralEvent=NULL,*prevInspiralEvent=NULL;

    SimInspiralIndex	      simIndex;
    SimInspiralTable	     *currentSimEvent=NULL,*prevSimEvent=NULL;
    SimInspiralTable	      simEvent,*simEventList=NULL;
    SimInspiralTable	     *simEventMissed=NULL;
    SimInspiralTable	     *prevSimEventMissed=NULL;

    MetadataTable             myTable;
    MetadataTable	      simFoundTable;
    MetadataTable	      simMissedTable;
    LIGOLwXMLStream           xmlStream;
    INT4		      numEvents = 0;
    INT4		      numInjects = 0; 
    INT4		      numKeptInjects = 0; 
    INT4		      numFound = 0;

    static int		      verbose_flag;
    
    struct MetaioParseEnvironment triggerParseEnv;
    const MetaioParseEnv triggerEnv = &triggerParseEnv;

    
    /********************************************************************
     * BEGIN PARSE ARGUMENTS						*
     ********************************************************************/
    
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
	    {"inject",		required_argument,  0,	'b'},
	    {"outfile",		required_argument,  0,	'c'},
	    {"snrstar",		required_argument,  0,	'd'},
	    {"sort",		no_argument,	    0,	'e'},
	    {"noplayground",	required_argument,  0,	'f'},
	    {"deltat",		required_argument,  0,  'g'},
	    {"help",		no_argument,	    0,	'h'},    
	    {"cluster",		required_argument,  0,	'j'},
	    {"clusteralgorithm",required_argument,  0,	'k'},
	    {"missedinjections",required_argument,  0,	'l'},
	    {0, 0, 0, 0}
	};
	/* getopt_long stores the option index here. */
	int option_index = 0;
 
    	c = getopt_long ( argc, argv, "a:b:c:d:e:f:g:h:j:k:l:", 
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
	    /* give name of file containing the injections */
	    {
	    	injectfile = optarg;
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
	    /* sort the events in time */
	    {
		sort = TRUE;
	    }
	    break;
    
	case 'f':
	    /* don't restrict to the playground data */
	    {
		playground = FALSE;
	    }
	    break;

	case 'g':
	    /* provide the time discrepancy allowed when checking for 
	     * coincidence of injections and inspirals 
	     * If not specified, then default is 20ms */ 
	    {
		dt = atoi( optarg );
	    }
	    break;

	    
    
	case 'h':
	    /* print help */
	    {
		LALPrintError( USAGE, *argv );
		return INSPINJFIND_EARG;
	    }
	    break;

	case 'j':
	    /* cluster the events */ 
	    {
		clust = atoi( optarg );
		cluster = TRUE;
	    }
	    break;
    
	case 'k':
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
 
	case 'l':
	    /* save a copy of those injected events which were missed, 
	     * and specify the filename where those events should be stored */
	    {
		missedInjectionsFile = optarg;
		saveMissedInjections = TRUE;
	    }
	    break;
	    
	default:
	    {
		return INSPINJFIND_EARG;
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



    /* check validity of required arguments */

    /* check validity of input file */
    if ( inputfile == NULL)
    {
    	LALPrintError( "Input file must be specified\n" );
    	exit( 1 );
    }
    
     /* check outfile has been given */
    if ( ! outfileName )
    {
    	LALPrintError( "Outfile name must be specified\n" );
	exit( 1 );
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
    * OPEN THE FILE WITH THE INJECTIONS (if appropriate)
    ******************************************************************/ 

    if ( ! (injectfile == NULL))
    {	
      
    /****************************************************************
     *  open xml file at injection table 
     ****************************************************************/
        
	if ( (retVal = MetaioOpenTable( triggerEnv, 
		    		injectfile, simInspiralName)) !=0 )
	{
            fprintf(stderr, "Error opening injection file %s\n", injectfile );
            MetaioAbort( triggerEnv ); 
            exit(INSPINJFIND_EFILE);
        }

	/* Locate the relevant columns in the injection table */
        buildSimInspiralIndex(&stat, triggerEnv, &simIndex);

        numInjects=0;
	/* Loop over the injections */
	while (1) 
	{
	    /* get the next row */
	    retVal = MetaioGetRow(triggerEnv);
	    if ( retVal == -1 ) 
	    {
		printf( "Error while getting row from file %s\n",line);
	        MetaioAbort( triggerEnv ); 
	        return INSPINJFIND_EROW;
	    } 
	    else if ( retVal == 0 ) 
	    {
		/*-- Reached end of file --*/
	        break;
	    }

	    /* get the injection details */
	    getSimInspiralVars(&stat, triggerEnv, &simEvent, &simIndex);
             
	    /* check that the injection occured during playground time */
            if ( ( playground ) && 
                (simEvent.geocent_end_time.gpsSeconds - 729273613)%6370 > 600)
	    {
                continue;
            }

	    /* allocate memory for the injection */
	    if ( (simEventList) == NULL )
	    {
		currentSimEvent = simEventList = (SimInspiralTable *) 
			LALCalloc( 1 , sizeof(SimInspiralTable) );
	    }
	    else 
	    {
		currentSimEvent = (SimInspiralTable *) 
			LALCalloc( 1 , sizeof(SimInspiralTable) );
	    }

	   	    
	    /* copy event into linked list element */
	    memcpy(currentSimEvent, &simEvent, sizeof(SimInspiralTable) );
            
	    /* point to the next event */
	    currentSimEvent->next = NULL;
	    if (prevSimEvent != NULL) prevSimEvent->next = currentSimEvent;
	    prevSimEvent = currentSimEvent;
	    currentSimEvent = currentSimEvent->next;

	    numInjects++;
	}
		
	/* set currentSimEvent back to the beginning of the list */
	currentSimEvent = simEventList;
	prevSimEvent = NULL;
	numKeptInjects = 0;
	
	MetaioAbort(triggerEnv);
    }

    

    /*****************************************************************
     * OPEN FILE WITH LIST OF XML FILES (one file per line)
     ****************************************************************/
    if ( !(fpin = fopen(inputfile,"r")) )
    {
        LALPrintError("Could not open input file\n");
	exit (1);
    }

 
    /*****************************************************************
     * loop over the xml files
     *****************************************************************/
    numEvents=0;
    
    while ( getline(line, MAXSTR, fpin) )
    {
        fileCounter++;

        /**************************************************************
	 *  open xml file at search summary table 
	 **************************************************************/
        if ( (retVal = MetaioOpenTable( triggerEnv, line, 
					searchSummaryName)) !=0 )
	{
            fprintf(stderr, "Error opening injection file %s for %s\n", 
                    line, searchSummaryName );
            if (fileCounter>1)
	    {
		MetaioAbort( triggerEnv ); 
		exit(INSPINJFIND_EFILE);
            } 
	    else 
	    {
                fprintf(stderr,"Proceeding anyway\n");
            }
        }

        if ( !retVal )
	{
        /* get the search summary table */
	    LAL_CALL( buildSearchSummaryIndex(&stat, triggerEnv, 
				&searchSummaryIndex), &stat);
	    LAL_CALL( getSearchSummaryTable(&stat, triggerEnv, 
			&searchSummaryTable, &searchSummaryIndex), &stat);

	    /* check for events and playground */
	    if ( ( playground && 
                !(isPlayground(searchSummaryTable.out_start_time.gpsSeconds,
                           searchSummaryTable.out_end_time.gpsSeconds ))) )
	    {
		fprintf(stdout,"File %i %s: not in playground, continuing\n",
                    fileCounter,line);
		continue;
	    }
	    else 
	    {
				
		/* if there are injections, keep only those which occur
		 * during the times searched in the input files */
		if( numInjects != 0)
		{
		    while (currentSimEvent->geocent_end_time.gpsSeconds <
			   searchSummaryTable.out_end_time.gpsSeconds)
		    {
			/*  check if injection is before file start time */
			if (currentSimEvent->geocent_end_time.gpsSeconds < 
			    searchSummaryTable.out_start_time.gpsSeconds)
			{
			/* discard the current injection */    
			    if (prevSimEvent != NULL)
			    {
				prevSimEvent->next = currentSimEvent->next;
				LALFree(currentSimEvent);
				currentSimEvent = prevSimEvent->next;
			    }
			    else
			    {
				simEventList = simEventList->next;
				LALFree(currentSimEvent);
			        currentSimEvent = simEventList;
			    }
			}
			else 
			{
			/* keep the current injection */    
			    numKeptInjects++;
			    prevSimEvent = currentSimEvent;
			    currentSimEvent = currentSimEvent->next;
			}
		    }		
		}    
		/* check to see whether there are any events in the file */
		if ( searchSummaryTable.nevents == 0 )
	        {
		    fprintf(stdout,"File %i %s: no events, continuing\n",
			    fileCounter,line);
	            continue;
	        }
	        else 
	        {
		    fprintf(stdout,"File %i %s: processing\n", 
	    	        fileCounter,line);
	        }

		
		
	    }
	    /* close the stream */
	    MetaioAbort( triggerEnv );
        }
        
        /**************************************************************
	 *  open xml file at inspiral table 
	 **************************************************************/
	
        if ( (retVal = MetaioOpenTable( triggerEnv, line, 
					snglInspiralName)) !=0 )
	{
            fprintf(stderr, "Error opening injection file %s\n", line );
            MetaioAbort( triggerEnv ); 
            exit(INSPINJFIND_EFILE);
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
                return INSPINJFIND_EROW;
            } else if ( retVal == 0 ) {
                /*-- Reached end of file --*/
                break;
            }

            /* get the inspiral event */
            getSnglInspiralEvent(&stat, triggerEnv, 
			    &inspiralEvent, &tableIndex);

            /* check for events and playground */
            if ( ( playground ) && 
                    (inspiralEvent.end_time.gpsSeconds-729273613)%6370 >600 )
	    {
                continue;
            }

            /* check that event satisfies threshold */
            if (inspiralEvent.snr < snrstar)
	    {
                continue;
            }
	    	   
            /* allocate memory for the inspiral event */
            if ( (inspiralEventList) == NULL )
            {
                currentEvent=inspiralEventList= (SnglInspiralTable *) 
			LALCalloc( 1, sizeof(SnglInspiralTable) );
            }
            else 
            {
                currentEvent = (SnglInspiralTable *) 
			LALCalloc( 1, sizeof(SnglInspiralTable) );
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
    }

    /* sort the events */
    if (sort)
    {
	LAL_CALL( LALSortSnglInspiral( &stat, &inspiralEventList,
                        *LALCompareSnglInspiralByTime), &stat);
    }

        
    /* discard the remaining injections which occured after the last file */
    while (currentSimEvent != NULL) 
    {
	if (prevSimEvent != NULL)
	{
	    prevSimEvent->next = currentSimEvent->next;
	    LALFree(currentSimEvent);
	    currentSimEvent = prevSimEvent->next;
	}
	else
	{
	    simEventList = simEventList->next;
	    LALFree(currentSimEvent);
	    currentSimEvent = simEventList;
	}
    }
    
    fprintf( stdout, "number of unclustered events %d \n", numEvents);
    fprintf( stdout, "number of injections in playground %d \n", numInjects);	
    fprintf( stdout, "number of relevant injections %d \n", numKeptInjects);
    /*****************************************************************
     * open output xml file
     *****************************************************************/
    LAL_CALL( LALOpenLIGOLwXMLFile(&stat, &xmlStream, outfileName), &stat);

    
     
	/*************************************************************
	 * CHECK FOR COINCIDENCE OF ARRIVAL TIME BETWEEN INJECTIONS
	 * AND FOUND INSPIRALS 
	 **************************************************************/
    if ( ! (injectfile == NULL))
    {
	/* allowed delay time given by dt */
	
	numFound = 0;
	coincidence = FALSE;
	/*  Note: we are assuming that both the inspiral and injection events
	*  are time sorted */
    
	currentSimEvent = simEventList;
	prevSimEvent = NULL;
	currentInspiralEvent = inspiralEventList;
	prevInspiralEvent = NULL;

	/* loop over all injections */
	while ( currentInspiralEvent != NULL && currentSimEvent != NULL )
	{
	
	    /* compute the time in nanosec for the injection */
	    LALGPStoINT8(&stat, &simTime, 
			&(currentSimEvent->geocent_end_time));

	    /* compute the time in nanosec for the inspiral */
	    LALGPStoINT8(&stat, &inspiralTime, 
			&(currentInspiralEvent->end_time));
	    /*  check if inspiral is before injection 
	     *  (up to allowed difference -- dt) */
	    if ( ((inspiralTime) < (simTime + 1000000LL * dt)) )
	    {
		/*  check if inspiral is coincident with injection
		*  (up to allowed difference -- dt) */
		if ( ((inspiralTime) > (simTime - 1000000LL * dt)) )
		{
		    /* record the fact that there's an inspiral coincident 
		     * with the injection */
		    if (coincidence == FALSE)
		    {
			++numFound;
			coincidence = TRUE;
    		    }
		    /* keep the inspiral event */
		    prevInspiralEvent = currentInspiralEvent;
		    currentInspiralEvent = currentInspiralEvent->next;
    		}
		else
		/* discard the inspiral event and check the next one,
		* do nothing to the injection event */
		{
		   if (prevInspiralEvent != NULL)
		    {
			prevInspiralEvent->next = currentInspiralEvent->next;
			LALFree(currentInspiralEvent);
			currentInspiralEvent = prevInspiralEvent->next;
		    }
		    
		   else
		    {
			inspiralEventList = inspiralEventList->next;
			LALFree(currentInspiralEvent);
			currentInspiralEvent = inspiralEventList;	
		    }
		}
	    }
	    else if (coincidence == TRUE)
	    /* keep the injection event */
	    {
		prevSimEvent = currentSimEvent;
		currentSimEvent = currentSimEvent->next;
		coincidence = FALSE;
	    }
	    else if (coincidence == FALSE)
	    /* Put the injection event into the simEventMissed list and 
	     * check the next one, do nothing to the inspiral event. */
	    {
		if (prevSimEvent != NULL)
		{
		    prevSimEvent->next = currentSimEvent->next;
    		    currentSimEvent->next = NULL;
		    if(prevSimEventMissed != NULL)
		    {    
			prevSimEventMissed->next = currentSimEvent;
			prevSimEventMissed = prevSimEventMissed->next;
		    }
		    else
		    {
			prevSimEventMissed = simEventMissed = currentSimEvent;
		    }	
		    currentSimEvent = prevSimEvent->next;		    
		}
		else
		{
		    simEventList = simEventList->next;
		    currentSimEvent->next = NULL;
		    if(prevSimEventMissed != NULL)
		    {    
			prevSimEventMissed->next = currentSimEvent;
			prevSimEventMissed = prevSimEventMissed->next;
		    }
		    else
		    {
			prevSimEventMissed = simEventMissed = currentSimEvent;
		    } 
		    currentSimEvent = simEventList;
		}
	    }
	}
    
	fprintf( stdout, "done sorting, number of coincidences found is %d\n", 
			numFound);

	if( currentInspiralEvent == NULL )
	{
	/* Remove all remaining injection events since they occur 
	 * after the last recorded inspiral -- put them in the missed table 
	 */
		
	    while ( currentSimEvent != NULL )
	    {
		if (prevSimEvent != NULL)
		{
		    prevSimEvent->next = currentSimEvent->next;
    		    currentSimEvent->next = NULL;
		    if(prevSimEventMissed != NULL)
		    {    
			prevSimEventMissed->next = currentSimEvent;
			prevSimEventMissed = prevSimEventMissed->next;
		    }
		    else
		    {
			prevSimEventMissed = simEventMissed = currentSimEvent;
		    }	
		    currentSimEvent = prevSimEvent->next;		    
		}
		else
		{
		    simEventList = simEventList->next;
		    currentSimEvent->next = NULL;
		    if(prevSimEventMissed != NULL)
		    {    
			prevSimEventMissed->next = currentSimEvent;
			prevSimEventMissed = prevSimEventMissed->next;
		    }
		    else
		    {
			prevSimEventMissed = simEventMissed = currentSimEvent;
		    } 
		    currentSimEvent = simEventList;
		}
	
	    }
	}
	
	if( currentSimEvent == NULL )
	{
	    /* Remove all remaining inspiral events since they occur after 
	     * the last recorded injection */
	    while (currentInspiralEvent != NULL )	
	    {
		if (prevInspiralEvent != NULL)
		{
		    prevInspiralEvent->next = currentInspiralEvent->next;
		    LALFree(currentInspiralEvent);
		    currentInspiralEvent = prevInspiralEvent->next;
		}
		else
		{
		    inspiralEventList = inspiralEventList->next;
		    LALFree(currentInspiralEvent);
		    currentInspiralEvent = inspiralEventList;	
		}
    	    }
	}
	
	
	/* Write the found injections to the sim table */
	LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, 
				 sim_inspiral_table),&stat);
        simFoundTable.simInspiralTable = simEventList;
        LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, simFoundTable, 
                    sim_inspiral_table), &stat);
	LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);


        /* free the temporary memory containing the events */
        while (simEventList)
        {
            SimInspiralTable *currentSimEvent;
            currentSimEvent = simEventList;
            simEventList = simEventList->next;
            LALFree( currentSimEvent );
        }

           
    }
/*********************************************************************/

    
    /******************************************************************
     * Now, cluster the events which were found to be coincident
     * (within the allowed window) with the injections.  Note, here it 
     * is natural to cluster with a time window  
     *		    clust = 2 * dt
     ******************************************************************/
	
    /* cluster the events */
    if ( cluster && inspiralEventList)
    {
        LAL_CALL( LALClusterSnglInspiralTable( &stat, inspiralEventList,
                        clust, clusterchoice), &stat);
    }
         
    /* Write the results to the inspiral table */
    LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, 
				 sngl_inspiral_table),&stat);
    myTable.snglInspiralTable = inspiralEventList;
    LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, myTable, 
                    sngl_inspiral_table), &stat);

    /* free the temporary memory containing the events */
    while (inspiralEventList)
    {
        SnglInspiralTable *currentInspiralEvent;
        currentInspiralEvent = inspiralEventList;
        inspiralEventList = inspiralEventList->next;
        LALFree( currentInspiralEvent );
    }
    

    /* close the output file */
    LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);
    LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);
   

    /* Open the missedInjectionsFile and copy the missed injections into it */
    if(saveMissedInjections == TRUE)
    {    
	/* Write the missed injections to the sim table */
	LAL_CALL( LALOpenLIGOLwXMLFile(&stat, &xmlStream, 
		    missedInjectionsFile), &stat);
	LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, 
				 sim_inspiral_table),&stat);
        simMissedTable.simInspiralTable = simEventMissed;
        LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, simMissedTable, 
                    sim_inspiral_table), &stat);
	LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);
	LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);
    }
    
    /* free the temporary memory containing the missed injections */
    while (simEventMissed)
        {
            SimInspiralTable *currentSimEventMissed;
            currentSimEventMissed = simEventMissed;
            simEventMissed = simEventMissed->next;
            LALFree( currentSimEventMissed );
        }

    LALCheckMemoryLeaks();
    return 0;
}
