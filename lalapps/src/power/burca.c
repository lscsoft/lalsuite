#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataUtils.h>
#include <lalapps.h>

RCSID("$Id$");

#define INCA_EARG   1
#define INCA_EROW   2
#define INCA_EFILE  3

#define INCA_MSGEARG   "Error parsing arguments"
#define INCA_MSGROW    "Error reading row from XML table"
#define INCA_MSGEFILE  "Could not open file"

#define INCA_TIMEWARNING "only triggers before 00:00 Dec 31, \
2010 will be used"

#define BURCA_SEARCHSUMWARNING "The Search summary info of \
of the two ifos do not match "

#define TRUE     1
#define FALSE    0 
#define MAXIFO   16
#define MAXFILES 128
#define MSEC   (1000000LL)

/* Usage format string. */
#define USAGE "Usage: %s --ifo-a trigfile.a --ifo-b trigfile.b \
--start-time startCoincidence --stop-time stopCoincidence \
--drho-plus dRhoPlus --drho-minus dRhoMinus --dt deltat \
--outfile outfilename --noplayground [--help]\n"

         
 
/************************************************************************************
 * The main program
 ***********************************************************************************/
int main(int argc, char **argv)
{
    static LALStatus       stat;
    INT4                   usePlayground= 1;
    INT4                   verbose = FALSE;

    INT4                   isPlay = 0;
    INT4                   deltaT=0;
    INT8                   triggerTime = 0;
    INT8                   startCoincidenceNS=0;
    INT8                   endCoincidenceNS=0;
    LIGOTimeGPS            startCoincidence={0,0};
    LIGOTimeGPS            endCoincidence={977788813,0};
    LIGOTimeGPS            slideData={0,0};
    INT8                   slideDataNS=0;

    CHAR                 **trigFile;
    INT4                   nTrigFile[MAXIFO];

    CHAR                   *outfileName=NULL;
    SnglBurstTable         *currentEvent=NULL,*tmpEvent=NULL;
    SnglBurstTable         *prevEvent=NULL,*outEvent=NULL;
    SnglBurstTable         *coincidentEvents=NULL;
    SnglBurstTable        **burstEventList=NULL;
    SnglBurstTable        **currentTrigger=NULL;
    SnglBurstAccuracy       accParams;
    SearchSummaryTable      *searchSummary=NULL;
    MetadataTable           myTable;
    MetadataTable           searchsumm;
    LIGOLwXMLStream         xmlStream;
    INT4                    i, j;

    REAL8                   tmpOutStart[MAXIFO];
    REAL8                   tmpOutEnd[MAXIFO];
    REAL8                   tmpInStart[MAXIFO];
    REAL8                   tmpInEnd[MAXIFO];

    /* getopt arguments */
    struct option long_options[] =
    {
        /* these options set a flag */
        {"verbose",                 no_argument,       &verbose,           1 },
        {"noplayground",            no_argument,       &usePlayground,     0 },
        /* parameters used to generate calibrated power spectrum */
        {"ifo-a",                   required_argument, 0,                'a'},
        {"ifo-b",                   required_argument, 0,                'b'},
        {"drhoplus",                required_argument, 0,                'c'},
        {"drhominus",               required_argument, 0,                'd'},
        {"dt",                      required_argument, 0,                't'},
        {"start-time",              required_argument, 0,                'r'},
        {"stop-time",               required_argument, 0,                's'},
        {"slide-time",              required_argument, 0,                'X'},
        {"slide-time-ns",           required_argument, 0,                'Y'},
        {"outfile",                 required_argument, 0,                'o'},
        {"help",                    no_argument,       0,                'h'}, 
        {0, 0, 0, 0}
    };
    int c;

    /*******************************************************************
     * initialize things
     *******************************************************************/
    lal_errhandler = LAL_ERR_EXIT;
    set_debug_level( "2" );
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );


    trigFile = (CHAR **) LALCalloc( MAXIFO * MAXFILES, sizeof(CHAR*) );
    burstEventList = (SnglBurstTable **) LALCalloc( 2, 
            sizeof(SnglBurstTable) );
    currentTrigger    = (SnglBurstTable **) LALCalloc( 2, 
            sizeof(SnglBurstTable) );
    memset( &accParams, 0, sizeof(SnglBurstAccuracy) );
    memset( &nTrigFile, 0, MAXIFO * sizeof(INT4) );

    /*******************************************************************
     * BEGIN PARSE ARGUMENTS (inarg stores the current position)        *
     *******************************************************************/
    /* parse the arguments */
    while ( 1 )
    {
        /* getopt_long stores long option here */
        int option_index = 0;

        c = getopt_long_only( argc, argv, 
                "a:b:c:d:m:t:o:r:s:h", long_options, &option_index );

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
                /* trigger files for interferometer A */
                {
                    size_t len = strlen(optarg) + 1;
                    trigFile[nTrigFile[0]] = 
                        (CHAR *) LALCalloc( len, sizeof(CHAR) );
                    memcpy( trigFile[nTrigFile[0]], optarg, len);
                    nTrigFile[0]++;
                }
                break;

            case 'b':
                /* trigger files for interferometer B */
                {
                    size_t len = strlen(optarg) + 1;
                    *(trigFile+MAXFILES+nTrigFile[1]) = 
                        (CHAR *) LALCalloc( len, sizeof(CHAR) );
                    memcpy( *(trigFile+MAXFILES+nTrigFile[1]), optarg, len);
                    nTrigFile[1]++;
                }
                break;

            case 'c':
                /* SNR error upward */
                accParams.dRhoPlus = atof(optarg);
                break;

            case 'd':
                /* SNR error downward */
                accParams.dRhoMinus = atof(optarg);
                break;

            case 't':
                /* time coincidence window */
                accParams.dtime = atof(optarg) * MSEC;
                deltaT = accParams.dtime;
                break;
           
            case 'r':
                /* time coincidence window */
                startCoincidence.gpsSeconds = atoi(optarg);
                break;

            case 's':
                /* time coincidence window */
                endCoincidence.gpsSeconds = atoi(optarg);
                break;

            case 'X':
                slideData.gpsSeconds = (INT4) atoi( optarg );
                break;

            case 'Y':
                slideData.gpsNanoSeconds = (INT4) atoi( optarg );
                break;

            case 'o':
                /* the name of the output file */
                {
                  size_t len = strlen(optarg) + 1;
                  outfileName = (CHAR *) LALCalloc( len, sizeof(CHAR) );
                  memcpy( outfileName, optarg, len);
                }
                break;

            case 'h':
                /* help message */
                fprintf( stderr, USAGE , argv[0]);
                exit( 1 );
                break;

            case '?':
                fprintf( stderr, USAGE , argv[0]);
                exit( 1 );
                break;

            default:
                fprintf( stderr, "unknown error while parsing options\n" );
                fprintf( stderr, USAGE, argv[0] );
                exit( 1 );
        }
    }

    if ( *trigFile == NULL || (*trigFile+MAXFILES) == NULL ||
            outfileName == NULL ){
        LALPrintError( "Must supply two input and one output file\n" );
        return INCA_EARG;
    }

    /*******************************************************************
     * END PARSE ARGUMENTS                                              *
     *******************************************************************/
    if (endCoincidence.gpsSeconds == 977788813)
        fprintf(stderr,"Warning: %s\n", INCA_TIMEWARNING);

    /****************************************************************
     * Create the search summary table [Info form ifo-a is only     *
     * recorded as of now ]                                         *
     ***************************************************************/

    for(j=0 ; j<2; j++)
      {	
	/* Read in the search summary info of ifo*/
	SearchSummaryTableFromLIGOLw( &searchSummary, *(trigFile+j*MAXFILES));

        LAL_CALL( LALGPStoFloat( &stat,
				 &tmpOutStart[j], &(searchSummary->out_start_time) ),&stat );
	LAL_CALL( LALGPStoFloat( &stat,
				 &tmpOutEnd[j], &(searchSummary->out_end_time) ),&stat );
	LAL_CALL( LALGPStoFloat( &stat,
				 &tmpInStart[j], &(searchSummary->in_start_time) ),&stat );
	LAL_CALL( LALGPStoFloat( &stat,
				 &tmpInEnd[j], &(searchSummary->in_end_time) ),&stat );


	while (searchSummary)
	  {
	    SearchSummaryTable *thisEvent;
	    thisEvent = searchSummary;
	    searchSummary = searchSummary->next;
	    LALFree( thisEvent );
	  }
	searchSummary = NULL; 
      }
    /* Compare to check if the two IFO's have the same search summary info*/
    if ((tmpOutStart[0] != tmpOutStart[1]) || (tmpInEnd[0] != tmpInEnd[1]))
      {       
	fprintf(stderr,"Warning: %s\n",  BURCA_SEARCHSUMWARNING);
      }

    searchsumm.searchSummaryTable = (SearchSummaryTable *)
      LALCalloc( 1, sizeof(SearchSummaryTable) );

      /* Store the info. for output with the output of burca*/
      LAL_CALL( LALFloatToGPS( &stat,
			       &(searchsumm.searchSummaryTable->out_start_time), &tmpOutStart[0] ),&stat );
      LAL_CALL( LALFloatToGPS( &stat,
			       &(searchsumm.searchSummaryTable->out_end_time), &tmpOutEnd[0] ),&stat );
      LAL_CALL( LALFloatToGPS( &stat,
			       &(searchsumm.searchSummaryTable->in_start_time), &tmpInStart[0] ),&stat );
      LAL_CALL( LALFloatToGPS( &stat,
			       &(searchsumm.searchSummaryTable->in_end_time), &tmpInEnd[0] ),&stat );


   /*****************************************************************
     * loop over input files for both ifos
     *****************************************************************/
    if (verbose) fprintf(stdout,"Looping over files\n");
    for(j=0 ; j<2 ; j++)
    {

        currentEvent = tmpEvent = burstEventList[j] = NULL;
        for(i=0; i<nTrigFile[j] ; i++)
        {
            LAL_CALL( LALSnglBurstTableFromLIGOLw (&stat, &tmpEvent, 
                        *(trigFile+j*MAXFILES+i)), &stat);

            /* connect results to linked list */
            if (currentEvent == NULL)
            {
                burstEventList[j] = currentEvent = tmpEvent;
            }
            else
            {
                currentEvent->next = tmpEvent;
            }

            /* move to the end of the linked list for next input file */
            while (currentEvent->next != NULL)
            {
                currentEvent = currentEvent->next;
            }
            tmpEvent = currentEvent->next;
        }

        /* sort the triggers in time */
        if ( burstEventList[j] )
        {
            LAL_CALL( LALSortSnglBurst(&stat, &(burstEventList[j]), 
                        XLALCompareSnglBurstByStartTime ), &stat);
        }
    }

    /*****************************************************************
     * slide the data if requested:  always slide burstEventList[0]
     * and put it on a ring with start and end time given by 
     * startCoincidence and endCoincidence ........
     *****************************************************************/
    LAL_CALL( LALGPStoINT8( &stat, &slideDataNS, &slideData ), &stat );
    if ( slideDataNS )
    {
      INT8 tStart=0;
      INT8 tPeak=0;

      currentEvent = burstEventList[0];

      while( currentEvent != NULL ){

        LAL_CALL( LALGPStoINT8( &stat, &tStart,
              &(currentEvent->start_time)), &stat );
        LAL_CALL( LALGPStoINT8( &stat, &tPeak,
              &(currentEvent->peak_time)), &stat );

       tStart += slideDataNS;
       tPeak += slideDataNS;

        if ( (tStart-startCoincidenceNS)*(tStart-endCoincidenceNS) > 0 )
        {
          if ( tStart < startCoincidenceNS )
          {
            tStart += endCoincidenceNS - startCoincidenceNS;
            tPeak += endCoincidenceNS - startCoincidenceNS;
          }
          if ( tStart > endCoincidenceNS )
          {
            tStart += startCoincidenceNS - endCoincidenceNS;
            tPeak += startCoincidenceNS - endCoincidenceNS;
          }
        }

        LAL_CALL( LALINT8toGPS( &stat, &(currentEvent->start_time),
            &tStart), &stat );
        LAL_CALL( LALINT8toGPS( &stat, &(currentEvent->peak_time),
            &tPeak), &stat );

        currentEvent = currentEvent->next;
      }
    }


    /*****************************************************************
     * find the first trigger after coincidence start time
     *****************************************************************/
    if (verbose) fprintf(stdout,"Moving to first trigger in window\n");
    for(j=0 ; j<2 ; j++)
    {
        currentTrigger[j] = burstEventList[j];

        while ( (currentTrigger[j] != NULL) && 
                (currentTrigger[j]->start_time.gpsSeconds < startCoincidence.gpsSeconds) )
        {
            currentTrigger[j] = currentTrigger[j]->next;
        }
    }


    /*****************************************************************
     * outer loop over triggers from interferometer A
     ****************************************************************/
    if (verbose) fprintf(stdout,"Start loop over ifo A\n");
    while ( (currentTrigger[0] != NULL) && 
            (currentTrigger[0]->start_time.gpsSeconds < endCoincidence.gpsSeconds) )
    {
        INT8 ta, tb;

        LAL_CALL( LALGPStoINT8(&stat, &ta, &(currentTrigger[0]->peak_time)), &stat);

        /*catch up triggers from ifo B */
        while (currentTrigger[1] != NULL)
        {
            LAL_CALL( LALGPStoINT8(&stat, &tb, &(currentTrigger[1]->peak_time)), &stat);
            if (tb > ta-deltaT)
            {
                break;
            }
            currentTrigger[1] = currentTrigger[1]->next;
        }

        if (verbose) fprintf(stdout,"Start loop over ifo B\n");

        /* look for coincident events in B within the time window */
        LAL_CALL( LALGPStoINT8( &stat, &triggerTime, 
              &(currentTrigger[0]->start_time) ), &stat );

        LAL_CALL( LALINT8NanoSecIsPlayground( &stat, &isPlay, 
              &triggerTime ), &stat );

        if ( verbose )
        {
          if ( isPlay )
          {
            fprintf( stdout, "trigger is playground\n" );
          } 
          else
          {
            fprintf( stdout, "trigger is not playground\n" );
          }
        }

        /* if we are playground only and the trigger is in playground or
         * we are not using playground and the trigger is not in the 
         * playground  ... */

        if ( ( usePlayground && isPlay ) || ( ! usePlayground && ! isPlay) )
        {
            tmpEvent = currentTrigger[1];
            while (tmpEvent != NULL)
            {
                LAL_CALL( LALGPStoINT8(&stat, &tb, &(tmpEvent->peak_time)), &stat);
                if (tb > ta+deltaT)
                {
                    break;
                }
                else
                {
                    /* this is a LAL function which compares events */
                    LAL_CALL( LALCompareSnglBurst(&stat, currentTrigger[0],
                                tmpEvent, &accParams), &stat);
                }

                if (accParams.match )
                {
                    if (coincidentEvents == NULL)
                    {
                        outEvent = coincidentEvents = (SnglBurstTable *)
                            LALCalloc(1, sizeof(SnglBurstTable) );
                        prevEvent = outEvent;
                    }
                    else 
                    {
                        outEvent = (SnglBurstTable *)
                            LALCalloc(1, sizeof(SnglBurstTable) );
                        prevEvent->next = outEvent;
                    }
                    memcpy( outEvent, currentTrigger[0], sizeof(SnglBurstTable));
                    prevEvent = outEvent;
                    outEvent = outEvent->next = NULL;
                }
                tmpEvent = tmpEvent->next;
            }
        }

        currentTrigger[0] = currentTrigger[0]->next;
    }


    /*****************************************************************
     * open output xml file
     *****************************************************************/
    LAL_CALL( LALOpenLIGOLwXMLFile(&stat, &xmlStream, outfileName), &stat);

    /* write the search cummary table */

    LAL_CALL( LALBeginLIGOLwXMLTable( &stat, &xmlStream, search_summary_table ), 
	      &stat );
    LAL_CALL( LALWriteLIGOLwXMLTable( &stat, &xmlStream, searchsumm, 
				      search_summary_table ), &stat );
    LAL_CALL( LALEndLIGOLwXMLTable ( &stat, &xmlStream ), &stat );
    LALFree( searchsumm.searchSummaryTable );

    /*write the triggers */

    LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, sngl_burst_table), &stat);
    myTable.snglBurstTable = coincidentEvents;
    LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, myTable,
                    sngl_burst_table), &stat);
    LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);
    LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);


    /*****************************************************************
     * clean up the memory that has been allocated 
     *****************************************************************/
    while ( coincidentEvents != NULL)
    {
        prevEvent = coincidentEvents;
        coincidentEvents = coincidentEvents->next;
        LALFree( prevEvent );
    }

    for(j=0; j<2 ; j++)
    {
        tmpEvent = burstEventList[j];
        while (tmpEvent)
        {
            prevEvent = tmpEvent;
            tmpEvent = tmpEvent->next;
            LALFree( prevEvent );
        }
        for (i=0; i<nTrigFile[j] ; i++){
            LALFree( trigFile[j*MAXFILES +i] );
        }
    }

    LALFree( outfileName );
    LALFree( currentTrigger );
    LALFree( burstEventList );
    LALFree( trigFile );
    
    LALCheckMemoryLeaks();

    return 0;
}
