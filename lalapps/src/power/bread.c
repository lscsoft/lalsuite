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
    [--max-confidence maximum conf] [--noplayground] [--sort] [--min-duration min dur] \
    [--max-duration max dur] [--min-centralfreq min central_freq] \
    [--max-centralfreq max central_freq] [--max-bandwidth max bw] \
    [--min-bandwidth min bw] [--min-amplitude min amp] [--max-amplitude max amp] \
    [--min-snr min snr] [--max-snr max snr] [--help]\n"

#define SNGLBURSTREADER_EARG   1
#define SNGLBURSTREADER_EROW   2
#define SNGLBURSTREADER_EFILE  3

#define SNGLBURSTREADER_MSGEARG   "Error parsing arguments"
#define SNGLBURSTREADER_MSGROW    "Error reading row from XML table"
#define SNGLBURSTREADER_MSGEFILE  "Could not open file"

#define TRUE  1
#define FALSE 0


/*
 * Read a line of text from a file, striping the newline character if read.
 * Return an empty string on any error.
 */

static int getline(char *line, int max, FILE *fpin)
{
	char *end;

	if(!fgets(line, max, fpin))
		line[0] = '\0';
	end = strchr(line, '\n');
	if(end)
		*end = '\0';
	return(strlen(line));
}


/*
 * Command-line options
 */

struct options_t {
	int verbose;
	int playground;
	int noplayground;
	int sort;

	/* confidence threshold */
	BOOLEAN maxConfidenceFlag;
	REAL4 maxConfidence;

	/* time window */
	BOOLEAN trigStartTimeFlag;
	INT4 trigStartTime;
	BOOLEAN trigStopTimeFlag;
	INT4 trigStopTime;

	/* duration thresholds */
	BOOLEAN minDurationFlag;
	REAL4 minDuration;
	BOOLEAN maxDurationFlag;
	REAL4 maxDuration;

	/* central_freq threshold */
	BOOLEAN maxCentralfreqFlag;
	REAL4 maxCentralfreq;
	BOOLEAN minCentralfreqFlag;
	REAL4 minCentralfreq;

	/* bandwidth threshold */
	BOOLEAN maxBandwidthFlag;
	REAL4 maxBandwidth;

	/* amplitude threshold */
	BOOLEAN maxAmplitudeFlag;
	REAL4 maxAmplitude;
	BOOLEAN minAmplitudeFlag;
	REAL4 minAmplitude;

	/* snr threshold */
	BOOLEAN maxSnrFlag;
	REAL4 maxSnr;
	BOOLEAN minSnrFlag;
	REAL4 minSnr;
};


/*
 * Set defaults for command-line options.
 * (ANSI C structure initialization sucks)
 */

static void set_option_defaults(struct options_t *options)
{
	options->verbose = FALSE;
	options->playground = FALSE;
	options->noplayground = FALSE;
	options->sort = FALSE;

	options->maxConfidenceFlag = FALSE;
	options->maxConfidence = 0.0;

	options->trigStartTimeFlag = FALSE;
	options->trigStartTime = 0;
	options->trigStopTimeFlag = FALSE;
	options->trigStopTime = 0;

	options->minDurationFlag = FALSE;
	options->minDuration = 0.0;
	options->maxDurationFlag = FALSE;
	options->maxDuration = 0.0;

	options->maxCentralfreqFlag = FALSE;
	options->maxCentralfreq = 0.0;
	options->minCentralfreqFlag = FALSE;
	options->minCentralfreq = 0.0;

	options->maxBandwidthFlag = FALSE;
	options->maxBandwidth = 0.0;

	options->maxAmplitudeFlag = FALSE;
	options->maxAmplitude = 0.0;
	options->minAmplitudeFlag = FALSE;
	options->minAmplitude = 0.0;

	options->maxSnrFlag = FALSE;
	options->maxSnr = 0.0;
	options->minSnrFlag = FALSE;
	options->minSnr = 0.0;
}


/*
 * Parse command line arguments.
 */

static void parse_command_line(int argc, char **argv, struct options_t *options, char **infile, char **outfile)
{
	struct option long_options[] = {
		/* these options set a flag */
		{"verbose",         no_argument,        &options->verbose, TRUE},
		/* parameters which determine the output xml file */
		{"input",           required_argument,  NULL,  'a'},
		{"outfile",         required_argument,  NULL,  'c'},
		{"max-confidence",  required_argument,  NULL,  'd'},
		{"min-duration",    required_argument,  NULL,  'e'},
		{"max-duration",    required_argument,  NULL,  'f'},
		{"min-centralfreq", required_argument,  NULL,  'g'},
		{"max-centralfreq", required_argument,  NULL,  'h'},
		{"max-bandwidth",   required_argument,  NULL,  'i'},
		{"min-amplitude",   required_argument,  NULL,  'j'},
		{"max-amplitude",   required_argument,  NULL,  'k'},
		{"min-snr",         required_argument,  NULL,  'l'},
		{"max-snr",         required_argument,  NULL,  'm'},
		{"trig-start-time", required_argument,  NULL,  'q'},
		{"trig-stop-time",  required_argument,  NULL,  'r'},
		{"playground",      no_argument,        &options->playground, TRUE},
		{"noplayground",    no_argument,        &options->noplayground, TRUE},
		{"help",            no_argument,        NULL,  'o'}, 
		{"sort",            no_argument,        &options->sort,  TRUE},
		{0, 0, 0, 0}
	};
	int c;
	int option_index;

	while(1) {
		c = getopt_long(argc, argv, "a:c:d:e:f:g:h:i:", long_options, &option_index);

		/* detect the end of the options */
		if(c == -1)
			break;

		switch(c) {
			case 0:
			/* if this option set a flag, do nothing else now */
			if ( long_options[option_index].flag != 0 )
				break;
			else {
				fprintf(stderr, "error parsing option %s with argument %s\n", long_options[option_index].name, optarg);
				exit( 1 );
			}
			break;

			case 'a':
			/*
			 * file containing list of xml files to use
			 */
			*infile = optarg;
			break;	

			case 'c':
			/*
			 * output file name
			 */
			*outfile = optarg;
			break;

			case 'd':
			/*
			 * the confidence must be smaller than this number
			 */
			options->maxConfidenceFlag = TRUE;
			options->maxConfidence = atof(optarg);
			break;

			case 'e':
			/*
			 * only events with duration greater than this are
			 * selected
			 */
			options->minDurationFlag = TRUE;
			options->minDuration = atof(optarg);
			break;

			case 'f':
			/*
			 * only events with duration less than this are
			 * selected
			 */
			options->maxDurationFlag = TRUE;
			options->maxDuration = atof(optarg);
			break;

			case 'g':
			/*
			 * only events with centralfreq greater than this
			 * are selected
			 */
			options->minCentralfreqFlag = TRUE;
			options->minCentralfreq = atof(optarg);
			break;

			case 'h':
			/*
			 * only events with centralfreq less than this are
			 * selected
			 */
			options->maxCentralfreqFlag = TRUE;
			options->maxCentralfreq = atof(optarg);
			break;

			case 'i': 
			/*
			 * only events with bandwidth less than this are
			 * selected
			 */
			options->maxBandwidthFlag = TRUE;
			options->maxBandwidth = atof(optarg);
			break;
    
			case 'j':
			/*
			 * only events with amp. more than this are
			 * selected
			 */
			options->minAmplitudeFlag = TRUE;
			options->minAmplitude = atof(optarg);
			break;

			case 'k':
			/*
			 * only events with amp. less than this are
			 * selected
			 */
			options->maxAmplitudeFlag = TRUE;
			options->maxAmplitude = atof(optarg);
			break;

			case 'l':
			/*
			 * only events with snr more than this are selected
			 */
			options->minSnrFlag = TRUE;
			options->minSnr = atof(optarg);
			break;

			case 'm':
			/*
			 * only events with snr less than this are selected
			 */
			options->maxSnrFlag = TRUE;
			options->maxSnr = atof(optarg);
			break;

			case 'r':
			/*
			 * only events with time before this are selected
			 */
			options->trigStopTimeFlag = TRUE;
			options->trigStopTime = atoi(optarg);
			break;

			case 'o':
			/*
			 * print help
			 */
			LALPrintError(USAGE, *argv);
			exit(SNGLBURSTREADER_EARG);

			default:
			exit(SNGLBURSTREADER_EARG);
		}   
	}

	if(optind < argc) {
		fprintf(stderr, "extraneous command line arguments:\n");
		while(optind < argc)
			fprintf(stderr, "%s\n", argv[optind++]);
		exit(SNGLBURSTREADER_EARG);
	}

	if(!*infile) {
		LALPrintError( "Must supply an xml file to parse\n" );
		exit(SNGLBURSTREADER_EARG);
	}

	if(!*outfile) {
		LALPrintError( "Outfile name must be specified\n" );
		exit(SNGLBURSTREADER_EARG);
	}
}


/*
 * Read and discard the search summary table from a LIGO LW burst trigger
 * file.  Return the time (seconds) encompassed by the file.
 */

static INT4 read_search_summary(char *filename, FILE *fpout)
{
	SearchSummaryTable *searchSummary = NULL;
	SearchSummaryTable *tmp;
	INT4 start;
	INT4 end;

	SearchSummaryTableFromLIGOLw(&searchSummary, filename);

	start = searchSummary->in_start_time.gpsSeconds;
	end = searchSummary->in_end_time.gpsSeconds;

	if(fpout)
		fprintf(fpout, "%d  %d  %d\n", start, end, end - start);

	while(searchSummary) {
		tmp = searchSummary;
		searchSummary = searchSummary->next;
		LALFree(tmp);
	}

	return(end - start);
}


/****************************************************************************
 * 
 * FUNCTION TESTS IF THE FILE CONTAINS ANY PLAYGROUND DATA
 * Remember to check if doing S2 or S3
 ***************************************************************************/
static int isPlayground(INT4 gpsStart, INT4 gpsEnd)
{
    INT4 runStart=729273613;
    INT4 playInterval=6370;
    INT4 playLength=600;
    INT4 segStart,segEnd,segMiddle;

    segStart = (gpsStart - runStart)%playInterval;
    segEnd   = (gpsEnd - runStart)%playInterval;
    segMiddle = gpsStart + (INT4) (0.5 * (gpsEnd - gpsStart));
    segMiddle = (segMiddle - runStart)%playInterval;
    
    return(segStart < playLength || segEnd < playLength || segMiddle < playLength);
}


/****************************************************************************
 *
 * The main program
 *
 *****************************************************************************/

int main(int argc, char **argv)
{
	static LALStatus  stat;
	FILE              *fpin=NULL;
	FILE              *fpout;

	/*number of events*/
	INT8              numEvents;

	/*searchsummary info */
	INT4              timeAnalyzed;

	CHAR              line[MAXSTR];

	CHAR              *infile=NULL,*outfile=NULL;
	SnglBurstTable    *tmpEvent=NULL,*currentEvent,*prevEvent=NULL;
	SnglBurstTable    burstEvent,*burstEventList,*outEventList=NULL;
	SnglBurstTable    **addpoint;
	MetadataTable     myTable;
	LIGOLwXMLStream   xmlStream;
	struct options_t  options;


	/*******************************************************************
	* initialize things
	*******************************************************************/

	set_option_defaults(&options);
	parse_command_line(argc, argv, &options, &infile, &outfile);

	lal_errhandler = LAL_ERR_EXIT;
	set_debug_level( "1" );
	memset( &burstEvent, 0, sizeof(SnglBurstTable) );
	memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
	xmlStream.fp = NULL;
	numEvents = 0;

	/*****************************************************************
	 * loop over the xml input files
	 *****************************************************************/

	if ( !(fpin = fopen(infile,"r")) )
		LALPrintError("Could not open list of input files\n");

	if(options.verbose) {
		fpout = fopen("./EPjobstartstop.dat","w");
		fprintf(fpout, "# This file contains the start & stop times of all jobs that succeded\n");
	} else
		fpout = NULL;

	burstEventList = NULL;
	addpoint = &burstEventList;
	timeAnalyzed = 0;

	while(getline(line, MAXSTR, fpin)) {
		if(options.verbose)
			fprintf(stderr, "Working on file %s\n", line);

		/*
		 * Read the search summary table
		 */

		timeAnalyzed += read_search_summary(line, fpout);

		/*
		 * Read the Sngl_Burst table
		 */

		LAL_CALL(LALSnglBurstTableFromLIGOLw(&stat, addpoint, line), &stat);

		/*
		 * move addpoint to the end of the linked list
		 */

		while(*addpoint)
			addpoint = &(*addpoint)->next;
	}

	/* print out the total time analysed */
	if(options.verbose) {
		fprintf(fpout, "# Total time analysed = %d\n", timeAnalyzed);
		fclose(fpout);
	}

	/****************************************************************
	* do any requested cuts
	***************************************************************/
	tmpEvent = burstEventList;
	while(tmpEvent) {
		BOOLEAN pass = TRUE;

		/* check if in after specified time window */
		if( options.trigStartTimeFlag && !(tmpEvent->start_time.gpsSeconds > options.trigStartTime) )
			pass = FALSE;

		/* check if in before specified end time */
		if( options.trigStopTimeFlag && !(tmpEvent->start_time.gpsSeconds < options.trigStopTime) )
			pass = FALSE;

		/* check the confidence */
		if( options.maxConfidenceFlag && !(tmpEvent->confidence < options.maxConfidence) )
			pass = FALSE;

		/* check min duration */
		if( options.minDurationFlag && !(tmpEvent->duration > options.minDuration) )
			pass = FALSE;

		/* check max duration */
		if( options.maxDurationFlag && !(tmpEvent->duration < options.maxDuration) )
			pass = FALSE;

		/* check min centralfreq */
		if( options.minCentralfreqFlag && !(tmpEvent->central_freq > options.minCentralfreq) )
			pass = FALSE;

		/* check max centralfreq */
		if( options.maxCentralfreqFlag && !(tmpEvent->central_freq < options.maxCentralfreq) )
			pass = FALSE;

		/* check max bandwidth */
		if( options.maxBandwidthFlag && !(tmpEvent->bandwidth < options.maxBandwidth) )
			pass = FALSE;

		/* check min amplitude */
		if( options.minAmplitudeFlag && !(tmpEvent->amplitude > options.minAmplitude) )
			pass = FALSE;

		/* check max amplitude */
		if( options.maxAmplitudeFlag && !(tmpEvent->amplitude < options.maxAmplitude) )
			pass = FALSE;

		/* check min snr */
		if( options.minSnrFlag && !(tmpEvent->snr > options.minSnr) )
			pass = FALSE;

		/* check max snr */
		if( options.maxSnrFlag && !(tmpEvent->snr < options.maxSnr) )
			pass = FALSE;

		/* check if trigger starts in playground */
		if ( options.playground && !(isPlayground(tmpEvent->peak_time.gpsSeconds, tmpEvent->peak_time.gpsSeconds)) )
			pass = FALSE;

		if ( options.noplayground && (isPlayground(tmpEvent->peak_time.gpsSeconds, tmpEvent->peak_time.gpsSeconds))  )
			pass = FALSE;
	
		/* set it for output if it passes */  
		if ( pass ) {
			if (outEventList == NULL) {
				outEventList = currentEvent = (SnglBurstTable *) LALCalloc(1, sizeof(SnglBurstTable) );
				prevEvent = currentEvent;
			} else {
				currentEvent = (SnglBurstTable *) LALCalloc(1, sizeof(SnglBurstTable) );
				prevEvent->next = currentEvent;
			}
			memcpy(currentEvent, tmpEvent, sizeof(SnglBurstTable));
			prevEvent = currentEvent;
			currentEvent = currentEvent->next = NULL;
		}
		tmpEvent = tmpEvent->next;
	}
  
	tmpEvent = outEventList;
	while( tmpEvent ) {
		tmpEvent = tmpEvent->next;
		numEvents++;
	}

	if(options.verbose)
		fprintf(stderr, "Total no. of triggers %ld\n", numEvents);


    /*****************************************************************
     * sort the triggers
     *****************************************************************/
    LAL_CALL(LALSortSnglBurst(&stat, &(outEventList), LALCompareSnglBurstByTime), &stat);



    /*****************************************************************
     * open output xml file
     *****************************************************************/
    LAL_CALL(LALOpenLIGOLwXMLFile(&stat, &xmlStream, outfile), &stat);
    LAL_CALL(LALBeginLIGOLwXMLTable (&stat, &xmlStream, sngl_burst_table), &stat);
    myTable.snglBurstTable = outEventList;
    LAL_CALL(LALWriteLIGOLwXMLTable (&stat, &xmlStream, myTable, sngl_burst_table), &stat);
    LAL_CALL(LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);
    LAL_CALL(LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);

	exit(0);
}
