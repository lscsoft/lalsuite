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

/* cluster options */
enum { undefined, clusterbypeaktimeandfreq, clusterbytimeandfreq } clusterchoice = undefined;

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
    INT4                   cluster = 0;

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

    CHAR                   *outburstfileName=NULL;
    CHAR                   *outinspiralfileName=NULL;

    SnglBurstTable         *currentBurstEvent=NULL,*tmpBurstEvent=NULL;
    SnglBurstTable         *prevBurstEvent=NULL,*outBurstEvent=NULL;
    SnglBurstTable         *coincidentBurstEvents=NULL;
    SnglBurstTable         *burstEventList=NULL;
    SnglBurstTable         *currentBurstTrigger=NULL;

    SnglInspiralTable      *currentInspiralEvent=NULL,*tmpInspiralEvent=NULL;
    SnglInspiralTable      *prevInspiralEvent=NULL,*outInspiralEvent=NULL;
    SnglInspiralTable      *coincidentInspiralEvents=NULL;
    SnglInspiralTable      *inspiralEventList=NULL;
    SnglInspiralTable      *currentInspiralTrigger=NULL;

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
        {"burst-trig-file",         required_argument, 0,                'a'},
        {"inspiral-trig-file",      required_argument, 0,                'b'},
        {"drhoplus",                required_argument, 0,                'c'},
        {"drhominus",               required_argument, 0,                'd'},
        {"dt",                      required_argument, 0,                't'},
        {"start-time",              required_argument, 0,                'r'},
        {"stop-time",               required_argument, 0,                's'},
        {"slide-time",              required_argument, 0,                'X'},
        {"slide-time-ns",           required_argument, 0,                'Y'},
        {"cluster",                 required_argument, 0,                'f'},
        {"outburstfile",            required_argument, 0,                'o'},
        {"outinspiralfile",         required_argument, 0,                'e'},
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
                "a:b:c:d:m:t:o:e:f:r:s:h", long_options, &option_index );

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
                  outburstfileName = (CHAR *) LALCalloc( len, sizeof(CHAR) );
                  memcpy( outburstfileName, optarg, len);
                }
                break;

            case 'e':
                /* the name of the output file */
                {
                  size_t len = strlen(optarg) + 1;
                  outinspiralfileName = (CHAR *) LALCalloc( len, sizeof(CHAR) );
                  memcpy( outinspiralfileName, optarg, len);
                }
                break;

            case 'f':
	        /*
		 * set the cluster option
		 */			
	        {
		  if ( ! strcmp( "clusterbypeaktimeandfreq", optarg ) )
		    {
		      cluster = 1;
		      clusterchoice = clusterbypeaktimeandfreq;
		    }
		  else if ( ! strcmp( "clusterbytimeandfreq", optarg ) )
		    {
		      cluster = 1;
		      clusterchoice = clusterbytimeandfreq;
		    }
		  else
		    {
		      fprintf( stderr, "invalid argument to --cluster\n"
			       "unknown clusterchoice specified;\n"
			       "(must be one of:clusterbypeaktimeandfreq ,clusterbytimeandfreq )\n");
		    }
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
            outburstfileName == NULL ){
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
     * loop over input files for bursts
     *****************************************************************/
    if (verbose) fprintf(stdout,"Looping over burst trig files\n");

    currentBurstEvent = tmpBurstEvent = burstEventList = NULL;
    for(i=0; i<nTrigFile[0] ; i++)
      {
	LAL_CALL( LALSnglBurstTableFromLIGOLw (&stat, &tmpBurstEvent, 
					       *(trigFile+0*MAXFILES+i)), &stat);
	
	/* connect results to linked list */
	if (currentBurstEvent == NULL)
	  {
	    burstEventList = currentBurstEvent = tmpBurstEvent;
	  }
	else
	  {
	    currentBurstEvent->next = tmpBurstEvent;
	  }
	
	/* move to the end of the linked list for next input file */
	while (currentBurstEvent->next != NULL)
	  {
	    currentBurstEvent = currentBurstEvent->next;
	  }
	tmpBurstEvent = currentBurstEvent->next;
      }
    
    /* sort the triggers in time */
    if ( burstEventList )
      {
	LAL_CALL( LALSortSnglBurst(&stat, &(burstEventList), 
				   XLALCompareSnglBurstByStartTime ), &stat);
      }

   /*****************************************************************
     * loop over input files for inspirals
     *****************************************************************/
    if (verbose) fprintf(stdout,"Looping over inspiral trig files\n");

    currentInspiralEvent = tmpInspiralEvent = inspiralEventList = NULL;
    for(i=0; i<nTrigFile[1] ; i++)
      {
	LALSnglInspiralTableFromLIGOLw (&tmpInspiralEvent, 
					       *(trigFile+1*MAXFILES+i), 0, -1);
	
	/* connect results to linked list */
	if (currentInspiralEvent == NULL)
	  {
	    inspiralEventList = currentInspiralEvent = tmpInspiralEvent;
	  }
	else
	  {
	    currentInspiralEvent->next = tmpInspiralEvent;
	  }
	
	/* move to the end of the linked list for next input file */
	while (currentInspiralEvent->next != NULL)
	  {
	    currentInspiralEvent = currentInspiralEvent->next;
	  }
	tmpInspiralEvent = currentInspiralEvent->next;
      }
    
    /* sort the inspiral triggers in time */
    if ( inspiralEventList )
      {
	LAL_CALL( LALSortSnglInspiral(&stat, &(inspiralEventList), 
				   LALCompareSnglInspiralByTime ), &stat);
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

      currentBurstEvent = burstEventList;

      while( currentBurstEvent != NULL ){

        LAL_CALL( LALGPStoINT8( &stat, &tStart,
              &(currentBurstEvent->start_time)), &stat );
        LAL_CALL( LALGPStoINT8( &stat, &tPeak,
              &(currentBurstEvent->peak_time)), &stat );

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

        LAL_CALL( LALINT8toGPS( &stat, &(currentBurstEvent->start_time),
            &tStart), &stat );
        LAL_CALL( LALINT8toGPS( &stat, &(currentBurstEvent->peak_time),
            &tPeak), &stat );

        currentBurstEvent = currentBurstEvent->next;
      }
    }


    /*****************************************************************
     * find the first burst trigger after coincidence start time
     *****************************************************************/
    if (verbose) fprintf(stdout,"Moving to first burst trigger in window\n");

    currentBurstTrigger = burstEventList;

    while ( (currentBurstTrigger != NULL) && 
	    (currentBurstTrigger->start_time.gpsSeconds < startCoincidence.gpsSeconds) )
      {
	currentBurstTrigger = currentBurstTrigger->next;
      }

    /*****************************************************************
     * find the first inspiral trigger after coincidence start time
     *****************************************************************/
    if (verbose) fprintf(stdout,"Moving to first inspiral trigger in window\n");

    currentInspiralTrigger = inspiralEventList;

    while ( (currentInspiralTrigger != NULL) && 
	    (currentInspiralTrigger->end_time.gpsSeconds < startCoincidence.gpsSeconds) )
      {
	currentInspiralTrigger = currentInspiralTrigger->next;
      }

    /*****************************************************************
     * outer loop over burst triggers
     ****************************************************************/
    if (verbose) fprintf(stdout,"Start loop over burst triggers\n");
    while ( (currentBurstTrigger != NULL) && 
            (currentBurstTrigger->start_time.gpsSeconds < endCoincidence.gpsSeconds) )
    {
        INT8 burststart, burstend, inspiralend;

        LAL_CALL( LALGPStoINT8(&stat, &burststart, &(currentBurstTrigger->start_time)), &stat);
	burstend = burststart + 1e9 * currentBurstTrigger->duration;

        /*catch up inspiral triggers*/
        while (currentInspiralTrigger != NULL)
        {
            LAL_CALL( LALGPStoINT8(&stat, &inspiralend, &(currentInspiralTrigger->end_time)), &stat);
            if (inspiralend > burststart-deltaT)
            {
                break;
            }
            currentInspiralTrigger = currentInspiralTrigger->next;
        }


        /* look for coincident events in B within the time window */
        LAL_CALL( LALGPStoINT8( &stat, &triggerTime, 
              &(currentBurstTrigger->start_time) ), &stat );

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

	if (verbose) fprintf(stdout,"Start loop over inspiral triggers\n");

        if ( ( usePlayground && isPlay ) || ( ! usePlayground && ! isPlay) )
        {
            tmpInspiralEvent = currentInspiralTrigger;
            while (tmpInspiralEvent != NULL)
            {
                LAL_CALL( LALGPStoINT8(&stat, &inspiralend, &(tmpInspiralEvent->end_time)), &stat);
                if (inspiralend > burstend)
                    break;

                /* this is a LAL function which compares events */
                LAL_CALL( LALCompareSnglBurstSnglInspiral(&stat, currentBurstTrigger,
                            tmpInspiralEvent, &accParams.difference, deltaT), &stat);

                if (!accParams.difference)
                {
                    if (coincidentBurstEvents == NULL)
                    {
                        outBurstEvent = coincidentBurstEvents = (SnglBurstTable *)
                            LALCalloc(1, sizeof(SnglBurstTable) );
                        prevBurstEvent = outBurstEvent;
                    }
                    else 
                    {
                        outBurstEvent = (SnglBurstTable *)
                            LALCalloc(1, sizeof(SnglBurstTable) );
                        prevBurstEvent->next = outBurstEvent;
                    }
                    memcpy( outBurstEvent, currentBurstTrigger, sizeof(SnglBurstTable));
                    prevBurstEvent = outBurstEvent;
                    outBurstEvent = outBurstEvent->next = NULL;

		    if (coincidentInspiralEvents == NULL)
                    {
                        outInspiralEvent = coincidentInspiralEvents = (SnglInspiralTable *)
                            LALCalloc(1, sizeof(SnglInspiralTable) );
                        prevInspiralEvent = outInspiralEvent;
                    }
                    else 
                    {
                        outInspiralEvent = (SnglInspiralTable *)
                            LALCalloc(1, sizeof(SnglInspiralTable) );
                        prevInspiralEvent->next = outInspiralEvent;
                    }
                    memcpy( outInspiralEvent,tmpInspiralEvent , sizeof(SnglInspiralTable));
                    prevInspiralEvent = outInspiralEvent;
                    outInspiralEvent = outInspiralEvent->next = NULL;
                }
                tmpInspiralEvent = tmpInspiralEvent->next;
            }
        }

        currentBurstTrigger = currentBurstTrigger->next;
    }


    /* cluster the burst triggers if asked for */
    if(cluster && clusterchoice == clusterbypeaktimeandfreq && coincidentBurstEvents)
      LAL_CALL(LALClusterSnglBurstTable(&stat, &coincidentBurstEvents, XLALCompareSnglBurstByPeakTime, XLALCompareSnglBurstByPeakTimeAndFreq), &stat);
    else if (cluster && clusterchoice == clusterbytimeandfreq && coincidentBurstEvents)
      LAL_CALL(LALClusterSnglBurstTable(&stat, &coincidentBurstEvents,  NULL, XLALCompareSnglBurst), &stat);


    if ( coincidentInspiralEvents )
      {
	LAL_CALL( LALSortSnglInspiral(&stat, &(coincidentInspiralEvents), 
				   LALCompareSnglInspiralByTime ), &stat);
      }

    /*****************************************************************
     * open output xml file
     *****************************************************************/
    LAL_CALL( LALOpenLIGOLwXMLFile(&stat, &xmlStream, outburstfileName), &stat);

    /* write the search cummary table */

    LAL_CALL( LALBeginLIGOLwXMLTable( &stat, &xmlStream, search_summary_table ), 
	      &stat );
    LAL_CALL( LALWriteLIGOLwXMLTable( &stat, &xmlStream, searchsumm, 
				      search_summary_table ), &stat );
    LAL_CALL( LALEndLIGOLwXMLTable ( &stat, &xmlStream ), &stat );
    LALFree( searchsumm.searchSummaryTable );

    /*write the coincident burst triggers */

    LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, sngl_burst_table), &stat);
    myTable.snglBurstTable = coincidentBurstEvents;
    LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, myTable,
                    sngl_burst_table), &stat);
    LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);
    LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);

    /*write the coincident inspiral triggers */
    if ( outinspiralfileName ){
      LAL_CALL( LALOpenLIGOLwXMLFile(&stat, &xmlStream, outinspiralfileName), &stat);
      
      LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, sngl_inspiral_table), &stat);
      myTable.snglInspiralTable = coincidentInspiralEvents;
      LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, myTable,
					sngl_inspiral_table), &stat);
      LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);
      LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);
      LALFree( outinspiralfileName );
    }
    /*****************************************************************
     * clean up the memory that has been allocated 
     *****************************************************************/
    while ( coincidentBurstEvents != NULL)
    {
        prevBurstEvent = coincidentBurstEvents;
        coincidentBurstEvents = coincidentBurstEvents->next;
        LALFree( prevBurstEvent );
    }

    while ( coincidentInspiralEvents != NULL)
    {
        prevInspiralEvent = coincidentInspiralEvents;
        coincidentInspiralEvents = coincidentInspiralEvents->next;
	LAL_CALL( LALFreeSnglInspiral( &stat, &prevInspiralEvent), &stat );
    }

    while ( burstEventList != NULL)
    {
        prevBurstEvent = burstEventList;
        burstEventList = burstEventList->next;
        LALFree( prevBurstEvent );
    }
 
    while ( inspiralEventList != NULL)
    {
        prevInspiralEvent = inspiralEventList;
        inspiralEventList = inspiralEventList->next;
	LAL_CALL( LALFreeSnglInspiral( &stat, &prevInspiralEvent), &stat );
    }
   
    for(j=0; j<2 ; j++)
    {
      for (i=0; i<nTrigFile[j] ; i++){
	LALFree( trigFile[j*MAXFILES +i] );
      }
    }
    
    LALFree( outburstfileName );
    LALFree( trigFile );
    
    LALCheckMemoryLeaks();

    return 0;
}
