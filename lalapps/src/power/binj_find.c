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
#include <lal/BurstSearch.h>
#include <lalapps.h>

RCSID("$Id$");

#define MAXSTR 2048

/* Usage format string. */
#define USAGE "Usage: %s --input infile --outfile filename \
    [--snrstar snrstar] [--noplayground] [--sort] [--cluster msec] \
    [--clusterchoice choicenumber] [--help]\n"

#define BINJ_FIND_EARG   1
#define BINJ_FIND_EROW   2
#define BINJ_FIND_EFILE  3

#define BINJ_FIND_MSGEARG   "Error parsing arguments"
#define BINJ_FIND_MSGROW    "Error reading row from XML table"
#define BINJ_FIND_MSGEFILE  "Could not open file"

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
    
    if (segStart < playLength || segEnd < playLength || segMiddle < playLength){
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
    FILE                     *fpin=NULL;
    INT4                     fileCounter=0;
    BOOLEAN                  playground=TRUE;
    INT4                     sort=FALSE;
    INT4                     pass=TRUE;
    size_t                   len=0;
    INT4                     ndetected=0;

    /* filenames */
    CHAR                    *inputFile=NULL;
    CHAR                    *injectionFile=NULL;
    CHAR                    *outputFile=NULL;

    /* times of comparison */
    INT4                     gpsStartTime=0;
    INT4                     gpsEndTime=0;
    INT8                     injPeakTime=0;
    INT8                     burstStartTime=0;

    /* confidence threshold */
    INT4                     maxConfidenceFlag=0;
    REAL4                    maxConfidence=0.0;

    /* duration thresholds */
    INT4                     minDurationFlag=0;
    REAL4                    minDuration=0.0;
    INT4                     maxDurationFlag=0;
    REAL4                    maxDuration=0.0;

    /* central_freq threshold */
    INT4                     maxCentralfreqFlag=0;
    REAL4                    maxCentralfreq=0.0;
    INT4                     minCentralfreqFlag=0;
    REAL4                    minCentralfreq=0.0;

    /* bandwidth threshold */ 
    INT4                     maxBandwidthFlag=0;
    REAL4                    maxBandwidth=0.0;

    /* amplitude threshold */
    INT4                     maxAmplitudeFlag=0;
    REAL4                    maxAmplitude=0.0;
    INT4                     minAmplitudeFlag=0;
    REAL4                    minAmplitude=0.0;

    /* snr threshold */
    INT4                     maxSnrFlag=0;
    REAL4                    maxSnr=0.0;                           
    INT4                     minSnrFlag=0;
    REAL4                    minSnr=0.0;         

    CHAR                     line[MAXSTR];

    /* triggers */
    SnglBurstTable        *tmpEvent=NULL,*currentEvent=NULL,*prevEvent=NULL;
    SnglBurstTable         burstEvent,*burstEventList=NULL,*outEventList=NULL;

    /* injections */
    SimBurstTable            *simBurstList=NULL;
    SimBurstTable            *currentSimBurst=NULL;

    /* comparison parameters */
    SnglBurstAccuracy         accParams;
    
    /* Table outputs */
    MetadataTable             myTable;
    LIGOLwXMLStream           xmlStream;

    static int		      verbose_flag = 0;

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
	    {"trigfile",	required_argument,  0,	'a'},
	    {"injfile",		required_argument,  0,	'b'},
	    {"outfile",		required_argument,  0,	'c'},
	    {"max-confidence",	required_argument,  0,	'd'},
	    {"gps-start-time",	required_argument,  0,	'e'},
	    {"gps-end-time",	required_argument,  0,	'f'},
	    {"noplayground",	required_argument,  0,	'n'},
	    {"help",		no_argument,	    0,	'o'}, 
	    {"sort",		no_argument,	    0,	'p'},
	    {0, 0, 0, 0}
	};
	/* getopt_long stores the option index here. */
	int option_index = 0;
 
    	c = getopt_long ( argc, argv, "a:c:d:e:f:g:h:i:", 
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
            /* create storage for the input file name */
            len = strlen( optarg ) + 1;
            inputFile = (CHAR *) calloc( len, sizeof(CHAR));
            memcpy( inputFile, optarg, len );
	    break;	
	
	case 'b':
            /* create storage for the injection file name */
            len = strlen( optarg ) + 1;
            injectionFile = (CHAR *) calloc( len, sizeof(CHAR));
            memcpy( injectionFile, optarg, len );
	    break;	
	
	case 'c':
            /* create storage for the injection file name */
            len = strlen( optarg ) + 1;
            outputFile = (CHAR *) calloc( len, sizeof(CHAR));
            memcpy( outputFile, optarg, len );
	    break;	

	case 'd':
	    /* the confidence must be smaller than this number */
	    {
              maxConfidenceFlag = 1;
              maxConfidence = atof( optarg );
	    }
	    break;

	case 'e':
	    /* only events with duration greater than this are selected */
	    {
              gpsStartTime = atoi( optarg );
	    }
	    break;
    
	case 'f':
	    /* only events with duration less than this are selected */
	    {
              gpsEndTime = atoi( optarg );
	    }
	    break;

	case 'n':
	    /* don't restrict to the playground data */
	    {
		playground = FALSE;
	    }
	    break;
    
	case 'o':
	    /* print help */
	    {
		LALPrintError( USAGE, *argv );
		return BINJ_FIND_EARG;
	    }
	    break;

	case 'p':
	    /* sort the events in time */
	    {
		sort = TRUE;
	    }
	    break;
	    
	default:
	    {
		return BINJ_FIND_EARG;
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

    if ( ! outputFile || ! inputFile || ! injectionFile )
    {
    	LALPrintError( "Input file, injection file, and output file 
            name must be specified\n" );
	return BINJ_FIND_EARG;
    }


    /*******************************************************************
    * END PARSE ARGUMENTS                                              *
    *******************************************************************/


    /*******************************************************************
    * initialize things
    *******************************************************************/
    lal_errhandler = LAL_ERR_EXIT;
    set_debug_level( "1" );
    memset( &burstEvent, 0, sizeof(SnglBurstTable) );
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
    xmlStream.fp = NULL;

    
    /*****************************************************************
     * OPEN FILE WITH LIST OF XML FILES (one file per line)
     ****************************************************************/
    if ( !(fpin = fopen(inputFile,"r")) ){
      LALPrintError("Could not open input file\n");
    }


    /*****************************************************************
     * READ IN THE INJECTIONS
     ****************************************************************/
    if ( verbose_flag )
      fprintf(stdout, "Reading in SimBurst Table\n");

    LAL_CALL( LALSimBurstTableFromLIGOLw ( &stat, &simBurstList, injectionFile,
          gpsStartTime, gpsEndTime), &stat );


    /*****************************************************************
     * loop over the xml files
     *****************************************************************/
    currentEvent = tmpEvent = burstEventList = NULL;
    while ( getline(line, MAXSTR, fpin) ){
      fileCounter++;
      if (verbose_flag)
      {
        fprintf(stderr,"Working on file %s\n", line);
      }

      LAL_CALL( LALSnglBurstTableFromLIGOLw (&stat, &tmpEvent, 
            line), &stat);

      /* connect results to linked list */
      if (currentEvent == NULL)
      {
        burstEventList = currentEvent = tmpEvent;
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

    /****************************************************************
     * do any requested cuts
     ***************************************************************/
    tmpEvent = burstEventList;
    while ( tmpEvent ){

      /* check the confidence */
      if( maxConfidenceFlag && !(tmpEvent->confidence < maxConfidence) )
        pass = FALSE;

      /* check min duration */
      if( minDurationFlag && !(tmpEvent->duration > minDuration) )
        pass = FALSE;

      /* check max duration */
      if( maxDurationFlag && !(tmpEvent->duration < maxDuration) )
        pass = FALSE;

      /* check min centralfreq */
      if( minCentralfreqFlag && !(tmpEvent->central_freq > minCentralfreq) )
        pass = FALSE;

      /* check max centralfreq */
      if( maxCentralfreqFlag && !(tmpEvent->central_freq < maxCentralfreq) )
        pass = FALSE;

      /* check max bandwidth */
      if( maxBandwidthFlag && !(tmpEvent->bandwidth < maxBandwidth) )
        pass = FALSE;

      /* check min amplitude */
      if( minAmplitudeFlag && !(tmpEvent->amplitude > minAmplitude) )
        pass = FALSE;

      /* check max amplitude */
      if( maxAmplitudeFlag && !(tmpEvent->amplitude < maxAmplitude) )
        pass = FALSE;

      /* check min snr */
      if( minSnrFlag && !(tmpEvent->snr > minSnr) )
        pass = FALSE;

      /* check max snr */
      if( maxSnrFlag && !(tmpEvent->snr < maxSnr) )
        pass = FALSE;

      /* check if trigger starts in playground */
      if ( playground && !(isPlayground(tmpEvent->start_time.gpsSeconds,
              tmpEvent->start_time.gpsSeconds)) )
        pass = FALSE;

      /* set it for output if it passes */  
      if ( pass )
      {
        if (outEventList == NULL)
        {
          outEventList = currentEvent = (SnglBurstTable *)
            LALCalloc(1, sizeof(SnglBurstTable) );
          prevEvent = currentEvent;
        }
        else 
        {
          currentEvent = (SnglBurstTable *)
            LALCalloc(1, sizeof(SnglBurstTable) );
          prevEvent->next = currentEvent;
        }
        memcpy( currentEvent, tmpEvent, sizeof(SnglBurstTable));
        prevEvent = currentEvent;
        currentEvent = currentEvent->next = NULL;
      }
      tmpEvent = tmpEvent->next;
      pass = TRUE;
    }

    /*****************************************************************
     * sort the remaining triggers
     *****************************************************************/
    LAL_CALL( LALSortSnglBurst(&stat, &(outEventList), 
                        LALCompareSnglBurstByTime ), &stat);

    
    /*****************************************************************
     * first event in list
     *****************************************************************/
    currentSimBurst = simBurstList;
    currentEvent = outEventList;
    while ( currentSimBurst != NULL ){
      /* convert injection time to INT8 */
      LAL_CALL( LALGPStoINT8(&stat, &injPeakTime, 
            &(currentSimBurst->geocent_peak_time)), &stat);

      /* loop over the burst events */
      while( currentEvent != NULL )
      {

        /* convert injection time to INT8 */
        LAL_CALL( LALGPStoINT8(&stat, &burstStartTime,
              &(currentEvent->start_time)), &stat);

        if( injPeakTime < burstStartTime )
          break;

        LAL_CALL( LALCompareSimBurstAndSnglBurst(&stat,
              currentSimBurst, currentEvent, &accParams), &stat);

        if( accParams.match )
        {
          ndetected++;
          break;
        }

        currentEvent = currentEvent->next;
      } 

      currentSimBurst = currentSimBurst->next;
    }

    fprintf(stdout,"Detected %i injections\n",ndetected);


    /*****************************************************************
     * open output xml file
     *****************************************************************/
    LAL_CALL( LALOpenLIGOLwXMLFile(&stat, &xmlStream, outputFile), &stat);
    LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, sngl_burst_table), &stat);
    myTable.snglBurstTable = outEventList;
    LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, myTable,
                    sngl_burst_table), &stat);
    LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);
    LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);

    return 0;
}
