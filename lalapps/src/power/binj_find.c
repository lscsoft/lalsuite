#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/FrameCache.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataUtils.h>
#include <lalapps.h>

long long int llabs(long long int j);	/* LAL's not written in ANSI C.  Get over it */

RCSID("$Id$");

#define MAXSTR 2048
#define TRUE  1
#define FALSE 0

enum { undefined, comparebytimeandfreq, comparebytime } comparechoice = comparebytimeandfreq;
 
/*
 * =============================================================================
 *                            Command Line Options
 * =============================================================================
 */

struct options_t {
	CHAR *inputFile;
	CHAR *burstinjectionFile;
	CHAR *inspinjectionFile;
	CHAR *injmadeFile;
	CHAR *injFoundFile;
	CHAR *outSnglFile;

	int verbose;
	int printresult;

	int playground;
	int noplayground;

	int best_confidence;
	int best_peaktime;

	/* times of comparison */
	INT4 gpsStartTime;
	INT4 gpsEndTime;

	/* central_freq threshold */
	INT4 maxCentralfreqFlag;
	REAL4 maxCentralfreq;
	INT4 minCentralfreqFlag;
	REAL4 minCentralfreq;

	/* confidence threshold */
	INT4 maxConfidenceFlag;
	REAL4 maxConfidence;

	/* duration thresholds */
	INT4 minDurationFlag;
	REAL4 minDuration;
	INT4 maxDurationFlag;
	REAL4 maxDuration;

	/* bandwidth threshold */
	INT4 maxBandwidthFlag;
	REAL4 maxBandwidth;

	/* amplitude threshold */
	INT4 maxAmplitudeFlag;
	REAL4 maxAmplitude;
	INT4 minAmplitudeFlag;
	REAL4 minAmplitude;

	/* snr threshold */
	INT4 maxSnrFlag;
	REAL4 maxSnr;
	INT4 minSnrFlag;
	REAL4 minSnr;
};


static void set_option_defaults(struct options_t *options)
{
	const long S2StartTime = 729273613;	/* Feb 14 2003 16:00:00 UTC */
	const long S2StopTime = 734367613;	/* Apr 14 2003 15:00:00 UTC */

	options->inputFile = NULL;
	options->burstinjectionFile = NULL;
	options->inspinjectionFile = NULL;
	options->injmadeFile = NULL;
	options->injFoundFile = NULL;
	options->outSnglFile = NULL;

	options->verbose = FALSE;
	options->printresult = FALSE;

	options->playground = FALSE;
	options->noplayground = FALSE;

	options->best_confidence = FALSE;
	options->best_peaktime = FALSE;

	options->gpsStartTime = S2StartTime;
	options->gpsEndTime = S2StopTime;

	options->maxCentralfreqFlag = FALSE;
	options->maxCentralfreq = 0.0;
	options->minCentralfreqFlag = FALSE;
	options->minCentralfreq = 0.0;

	options->maxConfidenceFlag = FALSE;
	options->maxConfidence = 0.0;

	options->minDurationFlag = FALSE;
	options->minDuration = 0.0;
	options->maxDurationFlag = FALSE;
	options->maxDuration = 0.0;

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


static void print_usage(FILE *handle, char *prog)
{
	fprintf(handle,
"Usage: %s --input-trig <filename> --input-burstinj <filename>\n" \
"       --input-inspinj <filename>\n" \
"	--output-trig <filename> --output-inj-made <filename>\n" \
"	--output-inj-found <filename> --best-confidence|--best-peaktime\n" \
"	[--compare-choice comparebytimeandfreq(default)/comparebytime][--noplayground] [--help]\n", prog);
}


static void parse_command_line(int argc, char *argv[], struct options_t *options)
{
	struct option long_options[] = {
	        {"verbose",          no_argument, &options->verbose, TRUE},
	        {"printresult",      no_argument, &options->printresult, TRUE},
		{"playground",       no_argument, &options->playground, TRUE},
		{"noplayground",     no_argument, &options->noplayground, TRUE},
		{"best-confidence",  no_argument, &options->best_confidence, TRUE},
		{"best-peaktime",    no_argument, &options->best_peaktime, TRUE},
		{"help",             no_argument, NULL, 'o'},
		{"input-trig",       required_argument, NULL, 'a'},
		{"input-burstinj",   required_argument, NULL, 'b'},
		{"input-inspinj",    required_argument, NULL, 'j'},
		{"output-inj-made",  required_argument, NULL, 'c'},
		{"max-confidence",   required_argument, NULL, 'd'},
		{"gps-start-time",   required_argument, NULL, 'e'},
		{"gps-end-time",     required_argument, NULL, 'f'},
		{"min-centralfreq",  required_argument, NULL, 'g'},
		{"max-centralfreq",  required_argument, NULL, 'h'},
		{"output-inj-found", required_argument, NULL, 'i'},
		{"compare-choice",   required_argument, NULL, 'k'},
		{"output-trig",      required_argument, NULL, 'q'},
		{NULL, 0, NULL, 0}
	};
	int c, index;

	do switch(c = getopt_long(argc, argv, "a:b:j:c:d:e:f:g:h:k:i:oq:", long_options, &index)) {
	case 'a':
		options->inputFile = optarg;
		break;

	case 'b':
		options->burstinjectionFile = optarg;
		break;

	case 'j':
		options->inspinjectionFile = optarg;
		break;

	case 'c':
		options->injmadeFile = optarg;
		break;

	case 'd':
		options->maxConfidenceFlag = TRUE;
		options->maxConfidence = atof(optarg);
		break;

	case 'e':
		options->gpsStartTime = atoi(optarg);
		break;

	case 'f':
		options->gpsEndTime = atoi(optarg);
		break;

	case 'g':
		options->minCentralfreqFlag = TRUE;
		options->minCentralfreq = atof(optarg);
		break;

	case 'h':
		options->maxCentralfreqFlag = TRUE;
		options->maxCentralfreq = atof(optarg);
		break;

	case 'i':
		options->injFoundFile = optarg;
		break;

	case 'k':
		/*
		 * set the cluster option
		 */			
		{
		  if ( ! strcmp( "comparebytimeandfreq", optarg ) )
		    {
		      comparechoice = comparebytimeandfreq;
		    }
		  else if ( ! strcmp( "comparebytime", optarg ) )
		    {
		      comparechoice = comparebytime ;
		    }
		  else
		    {
		      fprintf( stderr, "invalid argument to --compare-choice\n"
			       "unknown comparechoice specified;\n"
			       "(must be one of:comparebytimeandfreq,comparebytime )\n");
			    }
		}
		break;
		

	case 'o':
		print_usage(stderr, argv[0]);
		exit(1);

	case 'q':
		options->outSnglFile = optarg;
		break;

	/* option sets a flag */
	case 0:
		break;

	/* end of arguments */
	case -1:
		break;

	/* unrecognized argument */
	case '?':
		print_usage(stderr, argv[0]);
		exit(1);

	/* missing argument for an option */
	case ':':
		print_usage(stderr, argv[0]);
		exit(1);
	} while(c != -1);

	if(optind < argc) {
		fprintf(stderr, "%s: error: extraneous command line arguments:\n", argv[0]);
		while(optind < argc)
			fprintf(stderr, "\t%s\n", argv[optind++]);
		print_usage(stderr, argv[0]);
		exit(1);
	}

	if(!(options->best_confidence ^ options->best_peaktime)) {
		fprintf(stderr, "%s: error: must specify exactly one of --best-confidence or --best-peaktime\n", argv[0]);
		print_usage(stderr, argv[0]);
		exit(1);
	}

	if(!options->inputFile || (!options->burstinjectionFile &&  !options->inspinjectionFile) || !options->injmadeFile  || !options->injFoundFile || !options->outSnglFile) {
		fprintf(stderr, "%s: error: must specify all of --input-trig, --input-inj(burst/insp), --output-trig, --output-inj-made, and --output-inj-found\n", argv[0]);
		print_usage(stderr, argv[0]);
		exit(1);
	}
}


/*
 * =============================================================================
 *                        Miscellaneous Functions
 * =============================================================================
 */

/*
 * Read a line of text from a file, striping the newline character if read.
 * Return an empty string on any error.
 */

static int getline(char *line, int max, FILE *file)
{
	char *end;

	if(!fgets(line, max, file))
		line[0] = '\0';
	end = strchr(line, '\n');
	if(end)
		*end = '\0';
	return(strlen(line));
}


/*
 * Test if the given interval contains any playground data.
 * FIXME: check if doing S2 or S3
 */

static int contains_playground(INT4 gpsStart, INT4 gpsEnd)
{
	const INT4 S2Start = 729273613;
	const INT4 playInterval = 6370;
	const INT4 playLength = 600;
	INT4 start_offset = gpsStart - S2Start;
	INT4 end_offset = gpsEnd - S2Start;

	return((start_offset % playInterval < playLength) || (end_offset / playInterval != start_offset / playInterval));
}


/*
 * Pick the best of two triggers.
 */

static SnglBurstTable *select_burstinj_event(LALStatus *stat, SimBurstTable *injection, SnglBurstTable *a, SnglBurstTable *b, struct options_t options)
{
	if(!a)
		return(b);
	if(!b)
		return(a);

	if(options.best_confidence)
		return(a->confidence < b->confidence ? a : b);

	if(options.best_peaktime) {
		INT8 injtime, atime, btime;

		if(! strcmp("ZENITH",injection->coordinates))
		  LAL_CALL(LALGPStoINT8(stat, &injtime, &injection->geocent_peak_time), stat);
		else {
		  if(! strcmp("H1", a->ifo))
		    LAL_CALL(LALGPStoINT8(stat, &injtime, &injection->h_peak_time), stat);
		  else if(! strcmp("H2", a->ifo))
		    LAL_CALL(LALGPStoINT8(stat, &injtime, &injection->h_peak_time), stat);
		  else if(! strcmp("L1", a->ifo))
		    LAL_CALL(LALGPStoINT8(stat, &injtime, &injection->l_peak_time), stat);
		}	
		
		LAL_CALL(LALGPStoINT8(stat, &atime, &a->peak_time), stat);
		LAL_CALL(LALGPStoINT8(stat, &btime, &b->peak_time), stat);
		
		return(llabs(atime - injtime) < llabs(btime - injtime) ? a : b);
	}
	
	return(NULL);
}

static SnglBurstTable *select_inspinj_event(LALStatus *stat, SimInspiralTable *injection, SnglBurstTable *a, SnglBurstTable *b, struct options_t options)
{
	if(!a)
		return(b);
	if(!b)
		return(a);

	if(options.best_confidence)
		return(a->confidence < b->confidence ? a : b);

	if(options.best_peaktime) {
		INT8 injtime, atime, btime;
		
		if(! strcmp("H1", a->ifo))
		  LAL_CALL(LALGPStoINT8(stat, &injtime, &injection->h_end_time), stat);
	  	else if(! strcmp("H2", a->ifo))
		  LAL_CALL(LALGPStoINT8(stat, &injtime, &injection->h_end_time), stat);
	  	else if(! strcmp("L1", a->ifo))
		  LAL_CALL(LALGPStoINT8(stat, &injtime, &injection->l_end_time), stat);
		

		LAL_CALL(LALGPStoINT8(stat, &atime, &a->peak_time), stat);
		LAL_CALL(LALGPStoINT8(stat, &btime, &b->peak_time), stat);
		
		return(llabs(atime - injtime) < llabs(btime - injtime) ? a : b);
	}
	
	return(NULL);
}
/*
 * =============================================================================
 *                             Injection Handling
 * =============================================================================
 */

/*
 * Trim a SimBurstTable.  The trim function takes as input the address of the
 * pointer pointing to the head of the trigger list.  Upon exit, this pointer
 * will be set to point to the new head of the trigger list.  The return value
 * is the address of the "next" pointer in the last element of the list (which
 * might be same as the pointer to the head of the list if all elements were
 * deleted).
 */

static int keep_this_burstinjection(SimBurstTable *injection, struct options_t options)
{
	if(options.minCentralfreqFlag && !(injection->freq > options.minCentralfreq))
		return(FALSE);
	if(options.maxCentralfreqFlag && !(injection->freq < options.maxCentralfreq))
		return(FALSE);
	if(options.playground && !contains_playground(injection->geocent_peak_time.gpsSeconds, injection->geocent_peak_time.gpsSeconds))
		return(FALSE);

	return(TRUE);
}

static SimBurstTable *free_this_burstinjection(SimBurstTable *injection)
{
	SimBurstTable *next = injection ? injection->next : NULL;
	LALFree(injection);
	return(next);
}

/* not used anymore, but might come in handy, so don't delete just yet */
#if 0
static void free_burstinjections(SimBurstTable *list)
{
	while(list)
		list = free_this_injection(list);
}
#endif

static SimBurstTable **trim_burstinjection_list(SimBurstTable **list, struct options_t options)
{
	SimBurstTable *injection = *list;

	while(*list && !keep_this_burstinjection(*list, options))
		*list = free_this_burstinjection(*list);

	if(!*list)
		return(list);

	for(injection = *list; injection->next; ) {
		if(keep_this_burstinjection(injection->next, options))
			injection = injection->next;
		else
			injection->next = free_this_burstinjection(injection->next);
	}
	return(&injection->next);
}

/*
 * Trim a SimInspiralTable.  The trim function takes as input the address of the
 * pointer pointing to the head of the trigger list.  Upon exit, this pointer
 * will be set to point to the new head of the trigger list.  The return value
 * is the address of the "next" pointer in the last element of the list (which
 * might be same as the pointer to the head of the list if all elements were
 * deleted).
 */

static int keep_this_inspinjection(SimInspiralTable *injection, struct options_t options)
{
	if(options.playground && !contains_playground(injection->geocent_end_time.gpsSeconds, injection->geocent_end_time.gpsSeconds))
		return(FALSE);

	return(TRUE);
}

static SimInspiralTable *free_this_inspinjection(SimInspiralTable *injection)
{
	SimInspiralTable *next = injection ? injection->next : NULL;
	LALFree(injection);
	return(next);
}

/* not used anymore, but might come in handy, so don't delete just yet */
#if 0
static void free_inspinjections(SimInspiralTable *list)
{
	while(list)
		list = free_this_injection(list);
}
#endif

static SimInspiralTable **trim_inspinjection_list(SimInspiralTable **list, struct options_t options)
{
	SimInspiralTable *injection = *list;

	while(*list && !keep_this_inspinjection(*list, options))
		*list = free_this_inspinjection(*list);

	if(!*list)
		return(list);

	for(injection = *list; injection->next; ) {
		if(keep_this_inspinjection(injection->next, options))
			injection = injection->next;
		else
			injection->next = free_this_inspinjection(injection->next);
	}
	return(&injection->next);
}

/*
 * Read injection data.
 *
 * The input file should contain a list of injection XML files, one file name
 * per line.
 */

static SimBurstTable *read_burstinjection_list(LALStatus *stat, char *filename, INT4 start_time, INT4 end_time, struct options_t options)
{
	FILE *infile;
	char line[MAXSTR];
	SimBurstTable *list = NULL;
	SimBurstTable **addpoint = &list;

	if(options.verbose)
		fprintf(stderr, "Reading in SimBurst Table\n");

	if(!(infile = fopen(filename, "r")))
		LALPrintError("Could not open input file\n");

	while(getline(line, MAXSTR, infile)) {
		if (options.verbose)
			fprintf(stderr, "Working on file %s\n", line);

		LAL_CALL(LALSimBurstTableFromLIGOLw(stat, addpoint, line, start_time, end_time), stat);

		addpoint = trim_burstinjection_list(addpoint, options);
	}

	fclose(infile);

	return(list);
}

static SimInspiralTable *read_inspiralinjection_list(char *filename, INT4 start_time, INT4 end_time, struct options_t options)
{
	FILE *infile;
	char line[MAXSTR];
	SimInspiralTable *list = NULL;
	SimInspiralTable **addpoint = &list;

	if(options.verbose)
		fprintf(stderr, "Reading in SimInspiral Table\n");

	if(!(infile = fopen(filename, "r")))
		LALPrintError("Could not open input file\n");

	while(getline(line, MAXSTR, infile)) {
		if (options.verbose)
			fprintf(stderr, "Working on file %s\n", line);

		SimInspiralTableFromLIGOLw(addpoint, line, start_time, end_time);

		addpoint = trim_inspinjection_list(addpoint, options);

	}

	fclose(infile);

	return(list);
}



/*
 * Extract the burst injections that lie between the given start and end times.
 */

static SimBurstTable **extract_burstinjections(LALStatus *stat, SimBurstTable **addpoint, SimBurstTable *injection, INT8 start, INT8 end)
{
	INT8 peaktime;

	for(; injection; injection = injection->next) {
		LAL_CALL(LALGPStoINT8(stat, &peaktime, &injection->geocent_peak_time), stat);
		if((start < peaktime) && (peaktime < end)) {
			*addpoint = LALMalloc(sizeof(**addpoint));
			**addpoint = *injection;
			addpoint = &(*addpoint)->next;
			*addpoint = NULL;
		}
	}

	return(addpoint);
}

/*
 * Extract the inspiral injections that lie between the given start and end times.
 */

static SimInspiralTable **extract_inspiralinjections(LALStatus *stat, SimInspiralTable **addpoint, SimInspiralTable *injection, INT8 start, INT8 end)
{
	INT8 endtime;

	for(; injection; injection = injection->next) {
		LAL_CALL(LALGPStoINT8(stat, &endtime, &injection->geocent_end_time), stat);
		if((start < endtime) && (endtime < end)) {
			*addpoint = LALMalloc(sizeof(**addpoint));
			**addpoint = *injection;
			addpoint = &(*addpoint)->next;
			*addpoint = NULL;
		}
	}

	return(addpoint);
}
/*
 * =============================================================================
 *                             Injection Finding
 * =============================================================================
 */

/*
 * For each injection in the list point to by injection, the triggers in the
 * list pointed to by triglist are searched for the best match (if any).
 * 
 * When this exits, the pointer whose address was passed as detinj will be
 * pointing to the head of a list of injections that were found.  Likewise, the
 * pointer whose address was passed as dettrig will point to the list of
 * matching triggers.
 */

static void find_burstinjections(LALStatus *stat, SimBurstTable *injection, SnglBurstTable *triglist, SimBurstTable **detinj, SnglBurstTable **dettrig, int *ninjected, int *ndetected, struct options_t options)
{
	SnglBurstTable *event, *bestmatch;

	for(; injection; (*ninjected)++, injection = injection->next) {
		if(options.verbose)
			fprintf(stderr, "\tSearching for injection at time %d.%09d s\n", injection->geocent_peak_time.gpsSeconds, injection->geocent_peak_time.gpsNanoSeconds);

		bestmatch = NULL;
		for(event = triglist; event; event = event->next) {
			int match;

			/* if the injection's centre frequency and peak time
			 * don't lie within the trigger's time-frequency
			 * volume, move to next event */
			if (comparechoice == comparebytimeandfreq)
			  LAL_CALL(LALCompareSimBurstAndSnglBurst(stat, injection, event, XLALCompareSimBurstAndSnglBurstByTimeandFreq, &match), stat);
			else if (comparechoice == comparebytime)
			  LAL_CALL(LALCompareSimBurstAndSnglBurst(stat, injection, event, XLALCompareSimBurstAndSnglBurstByTime, &match), stat);

			if(!match)
				continue;

			/* compare this trigger to the best so far */
			bestmatch = select_burstinj_event(stat, injection, bestmatch, event, options);
		}

		/* if we didn't detect a matching event, continue to next
		 * injection */
		if(!bestmatch)
			continue;
		(*ndetected)++;

		/* record the detected injection */
		*detinj = LALMalloc(sizeof(**detinj));
		**detinj = *injection;
		detinj = &(*detinj)->next;
		*detinj = NULL;

		/* record the matching trigger */
		*dettrig = LALMalloc(sizeof(**dettrig));
		**dettrig = *bestmatch;
		dettrig = &(*dettrig)->next;
		*dettrig = NULL;
	}
}


static void find_inspinjections(LALStatus *stat, SimInspiralTable *injection, SnglBurstTable *triglist, SimInspiralTable **detinj, SnglBurstTable **dettrig, int *ninjected, int *ndetected, struct options_t options)
{
	SnglBurstTable *event, *bestmatch;

	for(; injection; (*ninjected)++, injection = injection->next) {
		if(options.verbose)
			fprintf(stderr, "\tSearching for injection at time %d.%09d s\n", injection->geocent_end_time.gpsSeconds, injection->geocent_end_time.gpsNanoSeconds);

		bestmatch = NULL;
		for(event = triglist; event; event = event->next) {
			int match;

			/* if the injection's centre frequency and peak time
			 * don't lie within the trigger's time-frequency
			 * volume, move to next event */
			LAL_CALL(LALCompareSimInspiralAndSnglBurst(stat, injection, event, &match), stat);
			if(!match)
				continue;

			/* compare this trigger to the best so far */
			bestmatch = select_inspinj_event(stat, injection, bestmatch, event, options);
		}

		/* if we didn't detect a matching event, continue to next
		 * injection */
		if(!bestmatch)
			continue;
		(*ndetected)++;

		/* record the detected injection */
		*detinj = LALMalloc(sizeof(**detinj));
		**detinj = *injection;
		detinj = &(*detinj)->next;
		*detinj = NULL;

		/* record the matching trigger */
		*dettrig = LALMalloc(sizeof(**dettrig));
		**dettrig = *bestmatch;
		dettrig = &(*dettrig)->next;
		*dettrig = NULL;
	}
}
/*
 * =============================================================================
 *                              Trigger Handling
 * =============================================================================
 */

/*
 * Determine the start and end times of an analysis from the search summary
 * table.
 */

static INT8 read_search_summary_start_end(LALStatus *stat, char *filename, INT8 *start, INT8 *end, FILE *fpout)
{
	SearchSummaryTable *searchSummary = NULL;
	INT8 local_start, local_end;

	/* allow for NULL pointers if the calling code doesn't care about these
	 * results */
	if(!start)
		start = &local_start;
	if(!end)
		end = &local_end;

	SearchSummaryTableFromLIGOLw(&searchSummary, filename);

	LAL_CALL(LALGPStoINT8(stat, start, &searchSummary->out_start_time), stat);
	LAL_CALL(LALGPStoINT8(stat, end, &searchSummary->out_end_time), stat);

	if(fpout)
		fprintf(fpout, "%lld  %lld  %lld\n", *start, *end, *end - *start);

	while(searchSummary) {
		SearchSummaryTable *tmp = searchSummary;
		searchSummary = searchSummary->next;
		LALFree(tmp);
	}

	return(*end - *start);
}


/*
 * Trim a SnglBurstTable.  The trim function takes as input the address of the
 * pointer pointing to the head of the trigger list.  Upon exit, this pointer
 * will be set to point to the new head of the trigger list.  The return value
 * is the address of the "next" pointer in the last element of the list (which
 * might be same as the pointer to the head of the list if all elements were
 * deleted).
 */

static int keep_this_event(SnglBurstTable *event, struct options_t options)
{
	if(options.maxConfidenceFlag && !(event->confidence < options.maxConfidence))
		return(FALSE);

	if(options.minDurationFlag && !(event->duration > options.minDuration))
		return(FALSE);

	if(options.maxDurationFlag && !(event->duration < options.maxDuration))
		return(FALSE);

	if(options.minCentralfreqFlag && !(event->central_freq > options.minCentralfreq))
		return(FALSE);

	if(options.maxCentralfreqFlag && !(event->central_freq < options.maxCentralfreq))
		return(FALSE);

	if(options.maxBandwidthFlag && !(event->bandwidth < options.maxBandwidth))
		return(FALSE);

	if(options.minAmplitudeFlag && !(event->amplitude > options.minAmplitude))
		return(FALSE);

	if(options.maxAmplitudeFlag && !(event->amplitude < options.maxAmplitude))
		return(FALSE);

	if(options.minSnrFlag && !(event->snr > options.minSnr))
		return(FALSE);

	if(options.maxSnrFlag && !(event->snr < options.maxSnr))
		return(FALSE);

	if(options.playground && !(contains_playground(event->start_time.gpsSeconds, event->start_time.gpsSeconds)))
		return(FALSE);

	return(TRUE);
}

static SnglBurstTable *free_this_event(SnglBurstTable *event)
{
	SnglBurstTable *next = event ? event->next : NULL;
	LALFree(event);
	return(next);
}

static void free_events(SnglBurstTable *list)
{
	while(list)
		list = free_this_event(list);
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
 * =============================================================================
 *                                Entry Point
 * =============================================================================
 */

int main(int argc, char **argv)
{
	/*
	 * Data.
	 */

	/* LAL */
	LALStatus stat;

	/* command line */
	struct options_t options;

	/* triggers */
	INT8 timeAnalyzed = 0;
	SnglBurstTable *trigger_list;
	SnglBurstTable *detectedTriggers = NULL;
	SnglBurstTable **dettrigaddpoint = &detectedTriggers;

	/* injections */
	INT4 ninjected = 0;
	INT4 ndetected = 0;

	/* burst injection tables */
	SimBurstTable *burstinjection_list = NULL;
	SimBurstTable *made_burstlist = NULL;
	SimBurstTable **made_burstaddpoint = &made_burstlist;
	SimBurstTable *detectedburstInjections = NULL;
	SimBurstTable **detburstinjaddpoint = &detectedburstInjections;

	/* inspiral injection tables */
	SimInspiralTable *inspinjection_list = NULL;
	SimInspiralTable *made_inspirallist = NULL;
	SimInspiralTable **made_inspiraladdpoint = &made_inspirallist;
	SimInspiralTable *detectedinspInjections = NULL;
	SimInspiralTable **detinspinjaddpoint = &detectedinspInjections;

	/* outputs */
	MetadataTable myTable;
	LIGOLwXMLStream xmlStream;

	/* input loop */
	FILE *infile;
	FILE *fp;
	char line[MAXSTR];
	INT8 SearchStart, SearchEnd;

	/*
	 * Initialize things.
	 */

	memset(&stat, 0, sizeof(stat));
	lal_errhandler = LAL_ERR_EXIT;
	set_debug_level("35");
	set_option_defaults(&options);
	parse_command_line(argc, argv, &options);

	/*
	 * Read and trim the burst/inspiral injection list.
	 */
	if(options.burstinjectionFile)
	  burstinjection_list = read_burstinjection_list(&stat, options.burstinjectionFile, options.gpsStartTime, options.gpsEndTime, options);
	else 
	  inspinjection_list = read_inspiralinjection_list(options.inspinjectionFile, options.gpsStartTime, options.gpsEndTime, options);

	/*
	 * Loop over trigger files, searching each for the appropriate
	 * injections.
	 */

	/* Open the list of trigger files */
	if(!(infile = fopen(options.inputFile, "r")))
		LALPrintError("Could not open input file\n");

	/* For each trigger file named in the input file... */
	while(getline(line, MAXSTR, infile)) {
		trigger_list = NULL;

		if(options.verbose)
			fprintf(stderr, "Working on file %s\n", line);

		/* Determine the times encompassed by this file */
		timeAnalyzed += read_search_summary_start_end(&stat, line, &SearchStart, &SearchEnd, NULL);

		/* Select the burst/inspiral injections made during these times */
		if( burstinjection_list )
		  extract_burstinjections(&stat, made_burstaddpoint, burstinjection_list, SearchStart, SearchEnd);
		else 
		  extract_inspiralinjections(&stat, made_inspiraladdpoint, inspinjection_list, SearchStart, SearchEnd);

		/* Read and trim the triggers from this file */
		LAL_CALL(LALSnglBurstTableFromLIGOLw(&stat, &trigger_list, line), &stat);
		trim_event_list(&trigger_list, options);

		/* Search the triggers for matches against the selected
		 * burst/inspiral injections 
		 */
		if( burstinjection_list )
		  find_burstinjections(&stat, *made_burstaddpoint, trigger_list, detburstinjaddpoint, dettrigaddpoint, &ninjected, &ndetected, options);
		else 
		  find_inspinjections(&stat, *made_inspiraladdpoint, trigger_list, detinspinjaddpoint, dettrigaddpoint, &ninjected, &ndetected, options);

		if( burstinjection_list ){
		  while(*detburstinjaddpoint)
		    detburstinjaddpoint = &(*detburstinjaddpoint)->next;
		  while(*made_burstaddpoint)
		    made_burstaddpoint = &(*made_burstaddpoint)->next;
		}
		else {
		  while(*detinspinjaddpoint)
		    detinspinjaddpoint = &(*detinspinjaddpoint)->next;
		  while(*made_inspiraladdpoint)
		    made_inspiraladdpoint = &(*made_inspiraladdpoint)->next;
		}

		while(*dettrigaddpoint)
			dettrigaddpoint = &(*dettrigaddpoint)->next;

		/* Clean up */
		free_events(trigger_list);
	}

	fclose(infile);

	/*
	 * Output some summary information.
	 */
	if(options.printresult){
	  fp = fopen("BinjFindResults.dat","w");
	  fprintf(fp,"%19.9f seconds = %.1f hours analyzed\n", timeAnalyzed / 1e9, timeAnalyzed / 3.6e12);
	  fprintf(fp, "Total injections: %d\n", ninjected);
	  fprintf(fp, "Total detected: %d\n", ndetected);
	  fprintf(fp, "Efficiency: %f\n", (double) ndetected / ninjected);
	  fclose(fp);
	}
	/*
	 * Write output XML files.
	 */

	memset(&xmlStream, 0, sizeof(LIGOLwXMLStream));
	xmlStream.fp = NULL;

	/* List of injections that were actually made */
	LAL_CALL(LALOpenLIGOLwXMLFile(&stat, &xmlStream, options.injmadeFile), &stat);
	if( made_burstlist ){
	  LAL_CALL(LALBeginLIGOLwXMLTable(&stat, &xmlStream, sim_burst_table), &stat);
	  myTable.simBurstTable = made_burstlist;
	  LAL_CALL(LALWriteLIGOLwXMLTable(&stat, &xmlStream, myTable, sim_burst_table), &stat);
	  LAL_CALL(LALEndLIGOLwXMLTable(&stat, &xmlStream), &stat);
	  LAL_CALL(LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);
	}
	else {
	  LAL_CALL(LALBeginLIGOLwXMLTable(&stat, &xmlStream, sim_inspiral_table), &stat);
	  myTable.simInspiralTable = made_inspirallist;
	  LAL_CALL(LALWriteLIGOLwXMLTable(&stat, &xmlStream, myTable, sim_inspiral_table), &stat);
	  LAL_CALL(LALEndLIGOLwXMLTable(&stat, &xmlStream), &stat);
	  LAL_CALL(LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);
	}

	/* List of injections which were detected */
	LAL_CALL(LALOpenLIGOLwXMLFile(&stat, &xmlStream, options.injFoundFile), &stat);
	if( detectedburstInjections ){
	  LAL_CALL(LALBeginLIGOLwXMLTable(&stat, &xmlStream, sim_burst_table), &stat);
	  myTable.simBurstTable = detectedburstInjections;
	  LAL_CALL(LALWriteLIGOLwXMLTable(&stat, &xmlStream, myTable, sim_burst_table), &stat);
	  LAL_CALL(LALEndLIGOLwXMLTable(&stat, &xmlStream), &stat);
	  LAL_CALL(LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);
	}
	else {
	  LAL_CALL(LALBeginLIGOLwXMLTable(&stat, &xmlStream, sim_inspiral_table), &stat);
	  myTable.simInspiralTable = detectedinspInjections;
	  LAL_CALL(LALWriteLIGOLwXMLTable(&stat, &xmlStream, myTable, sim_inspiral_table), &stat);
	  LAL_CALL(LALEndLIGOLwXMLTable(&stat, &xmlStream), &stat);
	  LAL_CALL(LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);
	}

	/* List of matching triggers */
	LAL_CALL(LALOpenLIGOLwXMLFile(&stat, &xmlStream, options.outSnglFile), &stat);
	LAL_CALL(LALBeginLIGOLwXMLTable(&stat, &xmlStream, sngl_burst_table), &stat);
	myTable.snglBurstTable = detectedTriggers;
	LAL_CALL(LALWriteLIGOLwXMLTable(&stat, &xmlStream, myTable, sngl_burst_table), &stat);
	LAL_CALL(LALEndLIGOLwXMLTable(&stat, &xmlStream), &stat);
	LAL_CALL(LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);

	exit(0);
}
