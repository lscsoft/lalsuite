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

#define SNGLBURSTREADER_EARG   1
#define SNGLBURSTREADER_EROW   2
#define SNGLBURSTREADER_EFILE  3

#define SNGLBURSTREADER_MSGEARG   "Error parsing arguments"
#define SNGLBURSTREADER_MSGROW    "Error reading row from XML table"
#define SNGLBURSTREADER_MSGEFILE  "Could not open file"

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

    /* confidence threshold */
    INT4                     maxConfidenceFlag=0;
    REAL4                    maxConfidence=0.0;

    /* duration thresholds */
    INT4                     minDurationFlag=0;
    REAL4                    minDuration=0.0;
    INT4                     maxDurationFlag=0;
    REAL4                    maxDuration=0.0;

    CHAR                     inputfile[MAXSTR],line[MAXSTR];

    CHAR                     *outfileName=NULL;
    SnglBurstTable        *tmpEvent=NULL,*currentEvent=NULL,*prevEvent=NULL;
    SnglBurstTable         burstEvent,*burstEventList=NULL,*outEventList=NULL;
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
	    {"input",		required_argument,  0,	'a'},
	    {"outfile",		required_argument,  0,	'c'},
	    {"max-confidence",	required_argument,  0,	'd'},
	    {"min-duration",	required_argument,  0,	'e'},
	    {"max-duration",	required_argument,  0,	'f'},
	    {"noplayground",	required_argument,  0,	'g'},
	    {"help",		no_argument,	    0,	'h'}, 
	    {"sort",		no_argument,	    0,	'i'},
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
	    /* file containing list of xml files to use */
	    {
		sprintf(inputfile,"%s",optarg);
    	    }
	    break;	
	
	case 'c':
	    /* output file name */
	    {
		outfileName = optarg;
	    }
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
              minDurationFlag = 1;
              minDuration = atof( optarg );
	    }
	    break;
    
	case 'f':
	    /* only events with duration less than this are selected */
	    {
              maxDurationFlag = 1;
              maxDuration = atof( optarg );
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
		return SNGLBURSTREADER_EARG;
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
		return SNGLBURSTREADER_EARG;
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
        return SNGLBURSTREADER_EARG;
    }

    if ( ! outfileName )
    {
    	LALPrintError( "Outfile name must be specified\n" );
	return SNGLBURSTREADER_EARG;
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
    if ( !(fpin = fopen(inputfile,"r")) ){
        LALPrintError("Could not open input file\n");
    }

    

    /*****************************************************************
     * loop over the xml files
     *****************************************************************/
    currentEvent = tmpEvent = burstEventList = NULL;
    while ( getline(line, MAXSTR, fpin) ){
      fileCounter++;
      if (verbose_flag)
      {
        fprintf(stderr,"WOrking on file %s\n", line);
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
     * open output xml file
     *****************************************************************/
    LAL_CALL( LALOpenLIGOLwXMLFile(&stat, &xmlStream, outfileName), &stat);
    LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, sngl_burst_table), &stat);
    myTable.snglBurstTable = outEventList;
    LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, myTable,
                    sngl_burst_table), &stat);
    LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);
    LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);

    return 0;
}
