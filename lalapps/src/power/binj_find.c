#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataUtils.h>
#include <lalapps.h>

long long int llabs(long long int j);	/* LAL's not written in ANSI C.  Get over it */

RCSID("$Id$");

#define MAXSTR 2048
#define TRUE  1
#define FALSE 0

/*
 * =============================================================================
 *                            Command Line Options
 * =============================================================================
 */

struct options_t {
	CHAR *inputFile;
	CHAR *injectionFile;
	CHAR *injmadeFile;
	CHAR *injFoundFile;
	CHAR *outSnglFile;

	int verbose;

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
	options->injectionFile = NULL;
	options->injmadeFile = NULL;
	options->injFoundFile = NULL;
	options->outSnglFile = NULL;

	options->verbose = FALSE;

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
"Usage: %s --input-trig <filename> --input-inj <filename>\n" \
"	--output-trig <filename> --output-inj-made <filename>\n" \
"	--output-inj-found <filename> --best-confidence|--best-peaktime\n" \
"	[--noplayground] [--help]\n", prog);
}


static void parse_command_line(int argc, char *argv[], struct options_t *options)
{
	struct option long_options[] = {
		{"verbose", no_argument, &options->verbose, TRUE},
		{"input-trig", required_argument, NULL, 'a'},
		{"input-inj", required_argument, NULL, 'b'},
		{"output-inj-made", required_argument, NULL, 'c'},
		{"max-confidence", required_argument, NULL, 'd'},
		{"gps-start-time", required_argument, NULL, 'e'},
		{"gps-end-time", required_argument, NULL, 'f'},
		{"min-centralfreq", required_argument, NULL, 'g'},
		{"max-centralfreq", required_argument, NULL, 'h'},
		{"output-inj-found", required_argument, NULL, 'i'},
		{"playground", no_argument, &options->playground, TRUE},
		{"noplayground", no_argument, &options->noplayground, TRUE},
		{"help", no_argument, NULL, 'o'},
		{"output-trig", required_argument, NULL, 'q'},
		{"best-confidence", no_argument, &options->best_confidence, TRUE},
		{"best-peaktime", no_argument, &options->best_peaktime, TRUE},
		{NULL, 0, NULL, 0}
	};
	int c, index;

	do switch(c = getopt_long(argc, argv, "a:c:d:e:f:g:h:i:oq:", long_options, &index)) {
	case 'a':
		options->inputFile = optarg;
		break;

	case 'b':
		options->injectionFile = optarg;
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

	if(!options->inputFile || !options->injectionFile || !options->injmadeFile  || !options->injFoundFile || !options->outSnglFile) {
		fprintf(stderr, "%s: error: must specify all of --input-trig, --input-inj, --output-trig, --output-inj-made, and --output-inj-found\n", argv[0]);
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
 * Function tests if the given interval contains any playground data.
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
 * Pick the best of two events.
 */

static SnglBurstTable *select_event(LALStatus *stat, SimBurstTable *injection, SnglBurstTable *a, SnglBurstTable *b, struct options_t options)
{
	if(!a)
		return(b);
	if(!b)
		return(a);

	if(options.best_confidence)
		return(a->confidence < b->confidence ? a : b);

	if(options.best_peaktime) {
		INT8 injtime, atime, btime;

		LAL_CALL(LALGPStoINT8(stat, &injtime, &injection->l_peak_time), stat);
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
 * Trim a SimBurstTable.
 */

static int keep_this_injection(SimBurstTable *injection, struct options_t options)
{
	if (options.minCentralfreqFlag && !(injection->freq > options.minCentralfreq))
		return(FALSE);
	if (options.maxCentralfreqFlag && !(injection->freq < options.maxCentralfreq))
		return(FALSE);
	if (options.playground && !contains_playground(injection->l_peak_time.gpsSeconds, injection->l_peak_time.gpsSeconds))
		return(FALSE);

	return(TRUE);
}

static SimBurstTable *free_this_injection(SimBurstTable *injection)
{
	SimBurstTable *next = injection ? injection->next : NULL;
	LALFree(injection);
	return(next);
}

static SimBurstTable *trim_injection_list(SimBurstTable *injection, struct options_t options)
{
	SimBurstTable *head;

	while(injection && !keep_this_injection(injection, options))
		injection = free_this_injection(injection);
	head = injection;

	if(injection)
		while(injection->next) {
			if(keep_this_injection(injection->next, options))
				injection = injection->next;
			else
				injection->next = free_this_injection(injection->next);
		}

	return(head);
}


/*
 * Read injection data.
 *
 * The input file should contain a list of injection XML files, one file name
 * per line.
 */

static SimBurstTable *read_injection_list(LALStatus *stat, char *filename, INT4 start_time, INT4 end_time, struct options_t options)
{
	FILE *infile;
	char line[MAXSTR];
	SimBurstTable *list = NULL;
	SimBurstTable **addpoint = &list;

	if (options.verbose)
		fprintf(stderr, "Reading in SimBurst Table\n");

	if (!(infile = fopen(filename, "r")))
		LALPrintError("Could not open input file\n");

	while (getline(line, MAXSTR, infile)) {
		if (options.verbose)
			fprintf(stderr, "Working on file %s\n", line);

		LAL_CALL(LALSimBurstTableFromLIGOLw(stat, addpoint, line, start_time, end_time), stat);

		*addpoint = trim_injection_list(*addpoint, options);

		while (*addpoint)
			addpoint = &(*addpoint)->next;
	}

	fclose(infile);

	return(list);
}


/*
 * Extract the injections that lie between the given start and end times.
 */

static SimBurstTable **extract_injections(LALStatus *stat, SimBurstTable **addpoint, SimBurstTable *injection, INT8 start, INT8 end)
{
	INT8 peaktime;

	for(; injection; injection = injection->next) {
		LAL_CALL(LALGPStoINT8(stat, &peaktime, &injection->l_peak_time), stat);
		if ((start < peaktime) && (peaktime < end)) {
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

	LAL_CALL(LALGPStoINT8(stat, start, &searchSummary->in_start_time), stat);
	LAL_CALL(LALGPStoINT8(stat, end, &searchSummary->in_end_time), stat);

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
 * Trim a SnglBurstTable
 */

static int keep_this_event(SnglBurstTable *event, struct options_t options)
{
	if (options.maxConfidenceFlag && !(event->confidence < options.maxConfidence))
		return(FALSE);

	if (options.minDurationFlag && !(event->duration > options.minDuration))
		return(FALSE);

	if (options.maxDurationFlag && !(event->duration < options.maxDuration))
		return(FALSE);

	if (options.minCentralfreqFlag && !(event->central_freq > options.minCentralfreq))
		return(FALSE);

	if (options.maxCentralfreqFlag && !(event->central_freq < options.maxCentralfreq))
		return(FALSE);

	if (options.maxBandwidthFlag && !(event->bandwidth < options.maxBandwidth))
		return(FALSE);

	if (options.minAmplitudeFlag && !(event->amplitude > options.minAmplitude))
		return(FALSE);

	if (options.maxAmplitudeFlag && !(event->amplitude < options.maxAmplitude))
		return(FALSE);

	if (options.minSnrFlag && !(event->snr > options.minSnr))
		return(FALSE);

	if (options.maxSnrFlag && !(event->snr < options.maxSnr))
		return(FALSE);

	if (options.playground && !(contains_playground(event->start_time.gpsSeconds, event->start_time.gpsSeconds)))
		return(FALSE);

	return(TRUE);
}

static SnglBurstTable *free_this_event(SnglBurstTable * event)
{
	SnglBurstTable *next = event ? event->next : NULL;
	LALFree(event);
	return(next);
}

static SnglBurstTable *trim_event_list(SnglBurstTable *event, struct options_t options)
{
	SnglBurstTable *head;

	while(event && !keep_this_event(event, options))
		event = event->next;
	head = event;

	if(event)
		while(event->next) {
			if(keep_this_event(event->next, options))
				event = event->next;
			else
				event->next = free_this_event(event->next);
		}

	return(head);
}


/*
 * Read trigger list, and generate a new injection list from the injections
 * that occur in the data that was analyzed.  We also trim out unwanted
 * triggers as we go to try to control memory usage.
 */

static SnglBurstTable *read_trigger_list(LALStatus *stat, char *filename, INT8 *timeAnalyzed, SimBurstTable **injection, struct options_t options)
{
	FILE *infile;
	char line[MAXSTR];
	SnglBurstTable *list = NULL;
	SnglBurstTable **eventaddpoint = &list;
	SimBurstTable *newinjection = NULL;
	SimBurstTable **injaddpoint = &newinjection;
	INT8 start, end;

	if(!(infile = fopen(filename, "r")))
		LALPrintError("Could not open input file\n");

	*timeAnalyzed = 0;
	while(getline(line, MAXSTR, infile)) {
		if(options.verbose)
			fprintf(stderr, "Working on file %s\n", line);

		*timeAnalyzed += read_search_summary_start_end(stat, line, &start, &end, NULL);

		injaddpoint = extract_injections(stat, injaddpoint, *injection, start, end);

		LAL_CALL(LALSnglBurstTableFromLIGOLw(stat, eventaddpoint, line), stat);

		*eventaddpoint = trim_event_list(*eventaddpoint, options);

		while(*eventaddpoint)
			eventaddpoint = &(*eventaddpoint)->next;
	}

	/*
	 * FIXME: LAL's memory management blows.  On Saikat's S3 injections,
	 * this loop takes about .25 seconds PER ITERATION!  (KCC)
	while(*injection) {
		SimBurstTable *tmp = *injection;
		*injection = (*injection)->next;
		LALFree(tmp);
	}
	*/
	*injection = newinjection;

	fclose(infile);

	return(list);
}


/*
 * =============================================================================
 *                                Entry Point
 * =============================================================================
 */

int main(int argc, char **argv)
{
	static LALStatus stat;

	INT4 ndetected;
	INT4 ninjected;
	INT8 timeAnalyzed;

	static struct options_t options;

	/* triggers */
	SnglBurstTable *event, *bestmatch;
	SnglBurstTable *burstEventList = NULL;
	SnglBurstTable *detectedTriggers = NULL;
	SnglBurstTable **detTriggersAddPoint = &detectedTriggers;

	/* injections */
	SimBurstTable *injection;
	SimBurstTable *simBurstList = NULL;
	SimBurstTable *detectedInjections = NULL;
	SimBurstTable **detInjectionsAddPoint = &detectedInjections;

	/* Table outputs */
	MetadataTable myTable;
	LIGOLwXMLStream xmlStream;

	/*
	 * Initialize things.
	 */

	lal_errhandler = LAL_ERR_EXIT;
	set_debug_level("1");

	/*
	 * Parse arguments.
	 */

	set_option_defaults(&options);
	parse_command_line(argc, argv, &options);

	/*
	 * Read and trim the injection list.
	 */

	simBurstList = read_injection_list(&stat, options.injectionFile, options.gpsStartTime, options.gpsEndTime, options);

	/*
	 * Read and trim the trigger list;  remove injections from the
	 * injection list that lie outside the time intervals that were
	 * actually analyzed according to the search summary tables.  Sort the
	 * trigger list too.
	 */

	burstEventList = read_trigger_list(&stat, options.inputFile, &timeAnalyzed, &simBurstList, options);
	if(options.verbose)
		fprintf(stderr, "Sorting triggers...\n");
	LAL_CALL(LALSortSnglBurst(&stat, &burstEventList, LALCompareSnglBurstByTime), &stat);

	/*
	 * For each injection, search the entire trigger list for the best
	 * match (if any).
	 */

	ninjected = ndetected = 0;
	for (injection = simBurstList; injection; ninjected++, injection = injection->next) {
		if(options.verbose)
			fprintf(stderr, "Searching for injection at time %d.%09d s\n", injection->l_peak_time.gpsSeconds, injection->l_peak_time.gpsNanoSeconds);

		bestmatch = NULL;
		for (event = burstEventList; event; event = event->next) {
			SnglBurstAccuracy accParams;

			/* if the injection's centre frequency and peak time
			 * don't lie within the trigger's time-frequency
			 * volume, move to next event */
			LAL_CALL(LALCompareSimBurstAndSnglBurst(&stat, injection, event, &accParams), &stat);
			if (!accParams.match)
				continue;

			/* compare this trigger to the best so far */
			bestmatch = select_event(&stat, injection, bestmatch, event, options);
		}

		/* if we didn't detect a matching event, continue to next
		 * injection */
		if(!bestmatch)
			continue;
		ndetected++;

		/* record the detected trigger */
		*detTriggersAddPoint = LALMalloc(sizeof(**detTriggersAddPoint));
		**detTriggersAddPoint = *bestmatch;
		detTriggersAddPoint = &(*detTriggersAddPoint)->next;
		*detTriggersAddPoint = NULL;

		/* record the detected injection */
		*detInjectionsAddPoint = LALMalloc(sizeof(**detInjectionsAddPoint));
		**detInjectionsAddPoint = *injection;
		detInjectionsAddPoint = &(*detInjectionsAddPoint)->next;
		*detInjectionsAddPoint = NULL;
	}

	fprintf(stdout,"%19.9f seconds = %.1f hours analyzed\n", timeAnalyzed / 1e9, timeAnalyzed / 3.6e12);
	fprintf(stdout, "Total injections: %d\n", ninjected);
	fprintf(stdout, "Total detected: %d\n", ndetected);
	fprintf(stdout, "Efficiency: %f\n", (double) ndetected / ninjected);

	/*
	 * Write output XML files.
	 */

	memset(&xmlStream, 0, sizeof(LIGOLwXMLStream));
	xmlStream.fp = NULL;

	/* List of injections that were actually made */
	LAL_CALL(LALOpenLIGOLwXMLFile(&stat, &xmlStream, options.injmadeFile), &stat);
	LAL_CALL(LALBeginLIGOLwXMLTable(&stat, &xmlStream, sim_burst_table), &stat);
	myTable.simBurstTable = simBurstList;
	LAL_CALL(LALWriteLIGOLwXMLTable(&stat, &xmlStream, myTable, sim_burst_table), &stat);
	LAL_CALL(LALEndLIGOLwXMLTable(&stat, &xmlStream), &stat);
	LAL_CALL(LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);

	/* List of injections which were detected */
	LAL_CALL(LALOpenLIGOLwXMLFile(&stat, &xmlStream, options.injFoundFile), &stat);
	LAL_CALL(LALBeginLIGOLwXMLTable(&stat, &xmlStream, sim_burst_table), &stat);
	myTable.simBurstTable = detectedInjections;
	LAL_CALL(LALWriteLIGOLwXMLTable(&stat, &xmlStream, myTable, sim_burst_table), &stat);
	LAL_CALL(LALEndLIGOLwXMLTable(&stat, &xmlStream), &stat);
	LAL_CALL(LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);

	/* List of triggers corresponding to injection */
	LAL_CALL(LALOpenLIGOLwXMLFile(&stat, &xmlStream, options.outSnglFile), &stat);
	LAL_CALL(LALBeginLIGOLwXMLTable(&stat, &xmlStream, sngl_burst_table), &stat);
	myTable.snglBurstTable = detectedTriggers;
	LAL_CALL(LALWriteLIGOLwXMLTable(&stat, &xmlStream, myTable, sngl_burst_table), &stat);
	LAL_CALL(LALEndLIGOLwXMLTable(&stat, &xmlStream), &stat);
	LAL_CALL(LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);

	exit(0);
}
