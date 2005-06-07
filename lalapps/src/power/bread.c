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
#include <lal/FrameCache.h>
#include <lalapps.h>

int snprintf(char *str, size_t size, const char *format, ...);

RCSID("$Id$");

#include <unistd.h>
#ifndef HAVE_GETHOSTNAME_PROTOTYPE
int gethostname(char *name, int len);
#endif

#ifndef MAXHOSTNAMELEN
#define MAXHOSTNAMELEN 256
#endif

#include <errno.h>
#include <unistd.h>
#include <sys/param.h>
#define STR( x ) #x
#define XSTR( x ) STR( x )
#define MAXPROTOCOLLEN 16

#define MAXSTR 2048

/* Usage format string. */
#define USAGE "Usage: %s --input filename --output filename " \
	"--trig-start-time time --trig-stop-time time " \
	"[--outtxt txt filename] " \
	"[--max-confidence maximum conf] [--noplayground] [--sort] " \
        "[--cluster clusterchoice(clusterbypeaktimeandfreq/clusterbytimeandfreq/stringcluster)] " \
	"[--min-duration min dur] [--max-duration max dur] " \
	"[--min-centralfreq min central_freq] [--max-centralfreq max central_freq] " \
	"[--max-bandwidth max bw] [--min-bandwidth min bw] " \
	"[--min-amplitude min amp] [--max-amplitude max amp] " \
	"[--min-snr min snr] [--max-snr max snr] [--globtrig dirname]" \
	"[--help] [--printsum]\n"

#define SNGLBURSTREADER_EARG   1
#define SNGLBURSTREADER_EROW   2
#define SNGLBURSTREADER_EFILE  3

#define SNGLBURSTREADER_MSGEARG   "Error parsing arguments"
#define SNGLBURSTREADER_MSGROW    "Error reading row from XML table"
#define SNGLBURSTREADER_MSGEFILE  "Could not open file"

#define TRUE  1
#define FALSE 0

/* cluster options */
enum { undefined, clusterbypeaktimeandfreq, clusterbytimeandfreq , stringcluster} clusterchoice = undefined;

/*
 * Command-line options
 */

struct options_t {
	int verbose;
	int printsum;
	int cluster;
	int globtrig;
	char *globdir;
	int outtxt;
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
	options->printsum = FALSE;
	options->cluster = FALSE;
	options->globtrig = FALSE;
	options->globdir = NULL;
	options->outtxt = FALSE;
	options->playground = FALSE;
	options->noplayground = FALSE;
	options->sort = FALSE;

	options->maxConfidenceFlag = FALSE;
	options->maxConfidence = 0.0;

	options->trigStartTimeFlag = FALSE;
	options->trigStartTime = -1;
	options->trigStopTimeFlag = FALSE;
	options->trigStopTime = -1;

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

static void parse_command_line(int argc, char **argv, struct options_t *options, char **infile, char **outfile, char **outtxt)
{
	struct option long_options[] = {
		/* these options set a flag */
		{"verbose",         no_argument,        &options->verbose, TRUE},
		{"printsum",        no_argument,        &options->printsum, TRUE},
		/* parameters which determine the output xml file */
		{"input",           required_argument,  NULL,  'a'},
		{"outtxt",          required_argument,  NULL,  'b'},
		{"output",          required_argument,  NULL,  'c'},
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
		{"cluster",         required_argument,  NULL,  'n'},
		{"trig-start-time", required_argument,  NULL,  'q'},
		{"trig-stop-time",  required_argument,  NULL,  'r'},
		{"globtrig",        required_argument,  NULL,  's'},
		{"playground",      no_argument,        &options->playground, TRUE},
		{"noplayground",    no_argument,        &options->noplayground, TRUE},
		{"help",            no_argument,        NULL,  'o'}, 
		{"sort",            no_argument,        &options->sort,  TRUE},
		{NULL, 0, NULL, 0}
	};
	int c;
	int option_index;

	*infile = *outfile = *outtxt = NULL;

	do {
		switch(c = getopt_long(argc, argv, "a:c:d:e:f:g:h:i:", long_options, &option_index)) {
			case -1:
			case 0:
			break;

			case 'a':
			/*
			 * file containing list of xml files to use
			 */
			*infile = optarg;
			break;

			case 'b':
			/*
			 * output txt file name
			 */
			options->outtxt = TRUE;
			*outtxt = optarg;
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
			
		        case 'n':
			/*
			 * set the cluster option
			 */			
			{
			  if ( ! strcmp( "clusterbypeaktimeandfreq", optarg ) )
			    {
			      options->cluster = TRUE;
			      clusterchoice = clusterbypeaktimeandfreq;
			    }
			  else if ( ! strcmp( "clusterbytimeandfreq", optarg ) )
			    {
			      options->cluster = TRUE;
			      clusterchoice = clusterbytimeandfreq;
			    }
			  else if ( ! strcmp( "stringcluster", optarg ) )
			    {
			      options->cluster = TRUE;
			      clusterchoice = stringcluster;
			    }
			  else
			    {
			      fprintf( stderr, "invalid argument to --cluster\n"
				       "unknown clusterchoice specified;\n"
				       "(must be one of:clusterbypeaktimeandfreq ,clusterbytimeandfreq )\n");
			    }
			}
			break;

			case 'm':
			/*
			 * only events with snr less than this are selected
			 */
			options->maxSnrFlag = TRUE;
			options->maxSnr = atof(optarg);
			break;

			case 'q':
			/*
			 * only events with time after this are selected
			 */
			options->trigStartTimeFlag = TRUE;
			options->trigStartTime = atoi(optarg);
			break;

			case 'r':
			/*
			 * only events with time before this are selected
			 */
			options->trigStopTimeFlag = TRUE;
			options->trigStopTime = atoi(optarg);
			break;

			case 's':
			/*
			 * only events with time before this are selected
			 */
			options->globtrig = TRUE;
			options->globdir = optarg;
			break;

			case ':':
			case '?':
			case 'o':
			default:
			/*
			 * print usage
			 */
			LALPrintError(USAGE, *argv);
			exit(SNGLBURSTREADER_EARG);
		}
	} while(c != -1);

	if(optind < argc) {
		fprintf(stderr, "extraneous command line arguments:\n");
		while(optind < argc)
			fprintf(stderr, "%s\n", argv[optind++]);
		exit(SNGLBURSTREADER_EARG);
	}

	if(!*infile && (options->globtrig == FALSE)) {
		LALPrintError( "Must supply an xml file to parse\n" );
		exit(SNGLBURSTREADER_EARG);
	}

	if(!*outtxt && options->outtxt) {
		LALPrintError( "Output txt file name must be specified\n" );
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

static INT8 read_search_summary_start_end(LALStatus *stat, char *filename, INT8 *start, INT8 *end, FILE *fpout)
{
	SearchSummaryTable *searchSummary = NULL;
	SearchSummaryTable *tmp;
	INT8 local_start, local_end;

	/* allow for NULL pointers if the calling code doesn't care about these
	 * results */
	if(!start)
		start = &local_start;
	if(!end)
		end = &local_end;

	SearchSummaryTableFromLIGOLw(&searchSummary, filename);

	if(!searchSummary)
		return(0);

	LAL_CALL(LALGPStoINT8(stat, start, &searchSummary->out_start_time), stat);
	LAL_CALL(LALGPStoINT8(stat, end, &searchSummary->out_end_time), stat);

	if(fpout)
		fprintf(fpout, "%lld  %lld  %lld\n", *start, *end, *end - *start);

	while(searchSummary) {
		tmp = searchSummary;
		searchSummary = searchSummary->next;
		LALFree(tmp);
	}

	return(*end - *start);
}


static void populate_search_summary_ifo(char *filename, MetadataTable *searchsumm)
{
	SearchSummaryTable *searchSummary = NULL;
	SearchSummaryTable *tmp;

	SearchSummaryTableFromLIGOLw(&searchSummary, filename);

	if(!searchSummary){
	  fprintf(stderr,"No search summary table found in input file\n");
	  exit(1);
	}

	snprintf(searchsumm->searchSummaryTable->ifos, LIGOMETA_IFOS_MAX, "%s", searchSummary->ifos);

	while(searchSummary) {
		tmp = searchSummary;
		searchSummary = searchSummary->next;
		LALFree(tmp);
	}
}


/*
 * Function tests if the given interval contains any playground data
 * FIXME: check if doing S2 or S3
 */

static int contains_playground(INT4 gpsStart, INT4 gpsEnd)
{
	INT4 S2Start = 729273613;
	INT4 playInterval = 6370;
	INT4 playLength = 600;
	INT4 start_offset = gpsStart - S2Start;
	INT4 end_offset = gpsEnd - S2Start;

	return((start_offset % playInterval < playLength) || (end_offset / playInterval != start_offset / playInterval));
}


/*
 * Trim unwanted events from the event list.
 */

static BOOLEAN keep_this_event(SnglBurstTable *event, struct options_t options)
{
	/* check if in after specified time window */
	if( options.trigStartTimeFlag && !(event->start_time.gpsSeconds >= options.trigStartTime) )
		return(FALSE);

	/* check if in before specified end time */
	if( options.trigStopTimeFlag && !(event->start_time.gpsSeconds <= options.trigStopTime) )
		return(FALSE);

	/* check the confidence */
	if( options.maxConfidenceFlag && !(event->confidence <= options.maxConfidence) )
		return(FALSE);

	/* check min duration */
	if( options.minDurationFlag && !(event->duration >= options.minDuration) )
		return(FALSE);

	/* check max duration */
	if( options.maxDurationFlag && !(event->duration <= options.maxDuration) )
		return(FALSE);

	/* check min centralfreq */
	if( options.minCentralfreqFlag && !(event->central_freq >= options.minCentralfreq) )
		return(FALSE);

	/* check max centralfreq */
	if( options.maxCentralfreqFlag && !(event->central_freq <= options.maxCentralfreq) )
		return(FALSE);

	/* check max bandwidth */
	if( options.maxBandwidthFlag && !(event->bandwidth <= options.maxBandwidth) )
		return(FALSE);

	/* check min amplitude */
	if( options.minAmplitudeFlag && !(event->amplitude >= options.minAmplitude) )
		return(FALSE);

	/* check max amplitude */
	if( options.maxAmplitudeFlag && !(event->amplitude <= options.maxAmplitude) )
		return(FALSE);

	/* check min snr */
	if( options.minSnrFlag && !(event->snr >= options.minSnr) )
		return(FALSE);

	/* check max snr */
	if( options.maxSnrFlag && !(event->snr <= options.maxSnr) )
		return(FALSE);

	/* check if trigger starts in playground */
	if ( options.playground && !(contains_playground(event->peak_time.gpsSeconds, event->peak_time.gpsSeconds)) )
		return(FALSE);

	if ( options.noplayground && (contains_playground(event->peak_time.gpsSeconds, event->peak_time.gpsSeconds))  )
		return(FALSE);
	
	return(TRUE);
}


static SnglBurstTable *free_this_event(SnglBurstTable *event)
{
	SnglBurstTable *next = event ? event->next : NULL;
	LALFree(event);
	return(next);
}


static SnglBurstTable **trim_event_list(SnglBurstTable **list, struct options_t options)
{
	SnglBurstTable *event;

	while(*list && !keep_this_event(*list, options))
		*list = free_this_event(*list);

	if(!*list)
		return(list);

	for(event = *list; event->next; ) {
		if(keep_this_event(event->next, options))
			event = event->next;
		else
			event->next = free_this_event(event->next);
	}
	return(&event->next);
}


/*
 * Count the number of events in a burst event list.
 */

static long int count_events(SnglBurstTable *event)
{
	long int count;

	for(count = 0; event; count++)
		event = event->next;
	
	return(count);
}


static int output_txt_file(FILE *fpout, SnglBurstTable *snglBursts){
  SnglBurstTable *thisEvent=NULL;

  thisEvent = snglBursts;
  /*fprintf(fpout,"# %s\n",thisEvent->ifo);
  fprintf(fpout,"# start_time,start_time_ns,peak_time,peak_time_ns,duration,central_freq,bandwidth,snr,confidence\n");
  fflush(fpout);*/

  while ( thisEvent ){
    fprintf(fpout,"%d %d %d %d %f %f %f %f %e\n",
	    thisEvent->start_time.gpsSeconds,
	    thisEvent->start_time.gpsNanoSeconds,
	    thisEvent->peak_time.gpsSeconds,
	    thisEvent->peak_time.gpsNanoSeconds,
	    thisEvent->duration,
	    thisEvent->central_freq,
	    thisEvent->bandwidth,
	    thisEvent->snr,
	    thisEvent->confidence
	    );
    thisEvent = thisEvent->next;
  }

  return 0;
}


static int path_from_fileurl( CHAR *filepath, FrStat *file )
{
  char prot[MAXPROTOCOLLEN + 1] = "";
  char host[MAXHOSTNAMELEN + 1] = "";
  int n;

  /* get protocol, hostname, and path */
  if ( ! file || ! file->url )
  {
	  fprintf(stderr,"No file passed in\n");
	  return 1;
  }
  n = sscanf( file->url, "%" XSTR( MAXPROTOCOLLEN ) "[^:]://%"
      XSTR( MAXHOSTNAMELEN ) "[^/]%" XSTR( MAXPATHLEN ) "s", prot, host, filepath );
  if ( n != 3 ) /* perhaps the hostname has been omitted */
  {
    n = sscanf( file->url, "%" XSTR( MAXPROTOCOLLEN ) "[^:]://%"
        XSTR( MAXPATHLEN ) "s", prot, filepath );
    if ( n == 2 )
      strcpy( host, "localhost" );
    else
    {
      strncpy( filepath, file->url, MAXPATHLEN );
      strcpy( prot, "none" );
    }
  }

  if ( ! strcmp( prot, "file" ) )
  {
    if ( strcmp( host, "localhost" ) )
    { /* make sure the host *is* localhost */
      char localhost[MAXHOSTNAMELEN + 1];
      gethostname( localhost, MAXHOSTNAMELEN );
      if ( strcmp( host, localhost ) )
      { /* nope */
        fprintf( stderr, "Can not read files from remote hosts.\n" );
        return 1;
      }
    }
  }
  else
  {
    fprintf( stderr, "Unsupported protocol %s.", prot );
    return 1;
  }

  return 0;
}





/*
 * Entry Point
 */

int main(int argc, char **argv)
{
	static LALStatus  stat;
	FILE              *fpout;
	INT8              timeAnalyzed;
	CHAR              line[MAXSTR];
	CHAR              *infile, *outfile, *outtxt;
	SnglBurstTable    *burstEventList;
	SnglBurstTable    **addpoint;
	MetadataTable     myTable;
	MetadataTable     searchsumm;
	LIGOLwXMLStream   xmlStream;
	struct options_t  options;
	FrCache *fileCache = NULL;
	FrCacheSieve fileCacheSieve;
	UINT4			  j;



	/*
	 * Initialize things
	 */

	set_option_defaults(&options);
	parse_command_line(argc, argv, &options, &infile, &outfile, &outtxt);

	lal_errhandler = LAL_ERR_EXIT;
	set_debug_level("1");
	memset(&xmlStream, 0, sizeof(LIGOLwXMLStream));
	xmlStream.fp = NULL;


	/*
	 * Store the times between which the triggers are selected in 
	 * a search summary table
	 */
	/* create the search summary table */
	searchsumm.searchSummaryTable = LALCalloc(1, sizeof(SearchSummaryTable));
	if(!searchsumm.searchSummaryTable->out_start_time.gpsSeconds){
	  searchsumm.searchSummaryTable->in_start_time.gpsSeconds = options.trigStartTime;
	  searchsumm.searchSummaryTable->out_start_time.gpsSeconds = options.trigStartTime;
	}
	if(!searchsumm.searchSummaryTable->out_end_time.gpsSeconds){
	  searchsumm.searchSummaryTable->in_end_time.gpsSeconds = options.trigStopTime;
	  searchsumm.searchSummaryTable->out_end_time.gpsSeconds = options.trigStopTime;
	}


	/*
	 * Loop over the xml input files
	 */

	/* open the file cache and create an internal cache if the file cache 
	 * is provided, else
	 * use glob to get the xml files into an internal cache
	 */
	if (options.globtrig)
	  LAL_CALL( LALFrCacheGenerate(&stat, &fileCache, options.globdir, "*.xml" ), &stat);
	else
	  LAL_CALL(LALFrCacheImport(&stat, &fileCache, infile), &stat);
	if ( options.trigStartTimeFlag || options.trigStopTimeFlag ){
		fileCacheSieve.srcRegEx=NULL;
		fileCacheSieve.dscRegEx=NULL;
		fileCacheSieve.urlRegEx=NULL;
		fileCacheSieve.earliestTime=options.trigStartTime;
		fileCacheSieve.latestTime=options.trigStopTime;
		fileCache = XLALFrSieveCache( fileCache, &fileCacheSieve );
	}


	if(options.printsum) {
		fpout = fopen("./EPjobstartstop.dat","w");
		fprintf(fpout, "# This file contains the start & stop times of all jobs that succeded\n");
	} else
		fpout = NULL;

	burstEventList = NULL;
	addpoint = &burstEventList;
	timeAnalyzed = 0;

	for (j=0 ; j<fileCache->numFrameFiles ; j++) {

		if ( path_from_fileurl(line, &(fileCache->frameFiles[j])) )
			exit(1);
		
		if(options.verbose)
			fprintf(stderr, "Working on file %s\n", fileCache->frameFiles[j].url);

		/*
		 * Read the search summary table
		 */

		timeAnalyzed += read_search_summary_start_end(&stat, line, NULL, NULL, fpout);


		populate_search_summary_ifo(line, &searchsumm);
		
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


	/*
	 * Cluster the events
	 */

	if(options.cluster && clusterchoice == clusterbypeaktimeandfreq)
	  LAL_CALL(LALClusterSnglBurstTable(&stat, &burstEventList, XLALCompareSnglBurstByPeakTime, XLALCompareSnglBurstByPeakTimeAndFreq), &stat);
	else if (options.cluster && clusterchoice == clusterbytimeandfreq)
	  LAL_CALL(LALClusterSnglBurstTable(&stat, &burstEventList,  NULL, XLALCompareSnglBurst), &stat);
	else if(options.cluster && clusterchoice == stringcluster)
	  XLALClusterStringBurstTable(&burstEventList, NULL, XLALCompareStringBurstByTime);

	/*
	 * Do any requested cuts
	 */

	trim_event_list(&burstEventList, options);

	/*
	 * print out the total time analysed and the number of triggers found
	 */
	if(options.printsum){
	  fprintf(fpout, "# Total time analysed = %lld nanosec.s\n", timeAnalyzed);
	  fprintf(fpout, "# Total no. of triggers = %ld\n", count_events(burstEventList));
	  fclose(fpout);
	}


	/*
	 * Sort the triggers
	 */

	LAL_CALL(LALSortSnglBurst(&stat, &burstEventList, XLALCompareSnglBurstByStartTime), &stat);


	/*
	 * Write output txt file if asked for
	 */

	if( options.outtxt ){	
	  fpout = fopen(outtxt,"w");  
	  output_txt_file(fpout, burstEventList);
	  fclose(fpout);
	}

	/*
	 * Write output xml file
	 */

	LAL_CALL(LALOpenLIGOLwXMLFile(&stat, &xmlStream, outfile), &stat);

	/* search summary table */
	searchsumm.searchSummaryTable->nevents = XLALCountSnglBurst(burstEventList);
	LAL_CALL(LALBeginLIGOLwXMLTable(&stat, &xmlStream, search_summary_table), &stat);
	LAL_CALL(LALWriteLIGOLwXMLTable(&stat, &xmlStream, searchsumm, search_summary_table), &stat);
	LAL_CALL(LALEndLIGOLwXMLTable(&stat, &xmlStream), &stat);
 	LALFree(searchsumm.searchSummaryTable);

	/* sngl_burst_table */
	LAL_CALL(LALBeginLIGOLwXMLTable (&stat, &xmlStream, sngl_burst_table), &stat);
	myTable.snglBurstTable = burstEventList;
	LAL_CALL(LALWriteLIGOLwXMLTable (&stat, &xmlStream, myTable, sngl_burst_table), &stat);
	LAL_CALL(LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);
	LAL_CALL(LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);

	exit(0);
}
