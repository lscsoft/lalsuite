#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/BurstSearch.h>
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

         
/****************************************************************************
 * FUNCTION TESTS IF TRIGGER IS IN PLAYGROUND DATA
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

    if (segStart < playLength || segEnd < playLength || segMiddle < playLength){
        return TRUE;
    }

    return FALSE;
}

 
/************************************************************************************
 * The main program
 ***********************************************************************************/
int main(int argc, char **argv)
{
    static LALStatus       stat;
    INT4                   playground=TRUE;
    INT4                   verbose = FALSE;

    INT4                   startCoincidence=0;
    INT4                   endCoincidence=977788813;
    INT4                   deltaT=0;
    INT4                   triggerTime = 0;

    CHAR                 **trigFile;
    INT4                   nTrigFile[MAXIFO];

    CHAR                   *outfileName=NULL;
    SnglBurstTable         *currentEvent=NULL,*tmpEvent=NULL;
    SnglBurstTable         *prevEvent=NULL,*outEvent=NULL;
    SnglBurstTable         *coincidentEvents=NULL;
    SnglBurstTable        **burstEventList=NULL;
    SnglBurstTable        **currentTrigger=NULL;
    SnglBurstAccuracy       accParams;
    MetadataTable           myTable;
    LIGOLwXMLStream         xmlStream;
    INT4                    i, j;

    /* getopt arguments */
    struct option long_options[] =
    {
        /* these options set a flag */
        {"verbose",                 no_argument,       &verbose,           TRUE },
        {"noplayground",            no_argument,       &playground,        FALSE },
        /* parameters used to generate calibrated power spectrum */
        {"ifo-a",                   required_argument, 0,                'a'},
        {"ifo-b",                   required_argument, 0,                'b'},
        {"drhoplus",                required_argument, 0,                'c'},
        {"drhominus",               required_argument, 0,                'd'},
        {"dt",                      required_argument, 0,                't'},
        {"start-time",              required_argument, 0,                'r'},
        {"stop-time",               required_argument, 0,                's'},
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
                startCoincidence = atoi(optarg);
                break;

            case 's':
                /* time coincidence window */
                endCoincidence = atoi(optarg);
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
    if (endCoincidence == 977788813)
        fprintf(stderr,"Warning: %s\n", INCA_TIMEWARNING);


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
                        LALCompareSnglBurstByTime ), &stat);
        }
    }


    /*****************************************************************
     * find the first trigger after coincidence start time
     *****************************************************************/
    if (verbose) fprintf(stdout,"Moving to first trigger in window\n");
    for(j=0 ; j<2 ; j++)
    {
        currentTrigger[j] = burstEventList[j];

        while (currentTrigger[j] != NULL && 
                (currentTrigger[j]->start_time.gpsSeconds < startCoincidence) )
        {
            currentTrigger[j] = currentTrigger[j]->next;
        }
    }


    /*****************************************************************
     * outer loop over triggers from interferometer A
     ****************************************************************/
    if (verbose) fprintf(stdout,"Start loop over ifo A\n");
    while ( (currentTrigger[0] != NULL) && 
            (currentTrigger[0]->start_time.gpsSeconds < endCoincidence) )
    {
        INT8 ta, tb;

        LAL_CALL( LALGPStoINT8(&stat, &ta, &(currentTrigger[0]->start_time)), &stat);

        /*catch up triggers from ifo B */
        while (currentTrigger[1] != NULL)
        {
            LAL_CALL( LALGPStoINT8(&stat, &tb, &(currentTrigger[1]->start_time)), &stat);
            if (tb > ta-deltaT)
            {
                break;
            }
            currentTrigger[1] = currentTrigger[1]->next;
        }

        if (verbose) fprintf(stdout,"Start loop over ifo B\n");
        /* look for coincident events in B within the time window */
        triggerTime = currentTrigger[0]->start_time.gpsSeconds;
        if ( (!(isPlayground(triggerTime, triggerTime)) && playground) == 0 )
        {
            tmpEvent = currentTrigger[1];
            while (tmpEvent != NULL)
            {
                LAL_CALL( LALGPStoINT8(&stat, &tb, &(tmpEvent->start_time)), &stat);
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
