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

RCSID("$Id$");

#define MAXSTR 2048
#define PLAYGROUND_START_TIME 729273613
#define PLAYGROUND_INTERVAL 6370
#define PLAYGROUND_LENGTH 600

/* Usage format string. */
#define USAGE "Usage: %s --input infile  --injfile injectionfile \
     --injmadefile filename --detsnglfile filename --injfoundfile filename \
     [--noplayground] [--help]\n"

#define BINJ_FIND_EARG   1
#define BINJ_FIND_EROW   2
#define BINJ_FIND_EFILE  3

#define BINJ_FIND_MSGEARG   "Error parsing arguments"
#define BINJ_FIND_MSGROW    "Error reading row from XML table"
#define BINJ_FIND_MSGEFILE  "Could not open file"

#define TRUE  1
#define FALSE 0

/*
 * =============================================================================
 *                            Command Line Options
 * =============================================================================
 */

struct options_t {
	int verbose;

	int playground;
	int noplayground;

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

/*
 * Set the options to their defaults.
 */

static void set_option_defaults(struct options_t *options)
{
	options->verbose = 0;

	options->playground = FALSE;
	options->noplayground = FALSE;

	options->maxCentralfreqFlag = 0;
	options->maxCentralfreq = 0.0;
	options->minCentralfreqFlag = 0;
	options->minCentralfreq = 0.0;

	options->maxConfidenceFlag = 0;
	options->maxConfidence = 0.0;

	options->minDurationFlag = 0;
	options->minDuration = 0.0;
	options->maxDurationFlag = 0;
	options->maxDuration = 0.0;

	options->maxBandwidthFlag = 0;
	options->maxBandwidth = 0.0;

	options->maxAmplitudeFlag = 0;
	options->maxAmplitude = 0.0;
	options->minAmplitudeFlag = 0;
	options->minAmplitude = 0.0;

	options->maxSnrFlag = 0;
	options->maxSnr = 0.0;
	options->minSnrFlag = 0;
	options->minSnr = 0.0;
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
 * Function tests if the file contains any playground data.
 * FIXME: check if doing S2 or S3
 */

static int isPlayground(INT4 gpsStart, INT4 gpsEnd)
{
	INT4 runStart = 729273613;
	INT4 playInterval = 6370;
	INT4 playLength = 600;
	INT4 segStart, segEnd, segMiddle;

	segStart = (gpsStart - runStart) % playInterval;
	segEnd = (gpsEnd - runStart) % playInterval;
	segMiddle = gpsStart + (INT4) (0.5 * (gpsEnd - gpsStart));
	segMiddle = (segMiddle - runStart) % playInterval;

	return(segStart < playLength || segEnd < playLength || segMiddle < playLength);
}


/*
 * =============================================================================
 *                             Injection Handling
 * =============================================================================
 */

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
		fprintf(stdout, "Reading in SimBurst Table\n");

	if (!(infile = fopen(filename, "r")))
		LALPrintError("Could not open input file\n");

	while (getline(line, MAXSTR, infile)) {
		if (options.verbose)
			fprintf(stderr, "Working on file %s\n", line);

		LAL_CALL(LALSimBurstTableFromLIGOLw(stat, addpoint, line, start_time, end_time), stat);

		while (*addpoint)
			addpoint = &(*addpoint)->next;
	}

	fclose(infile);

	return(list);
}


/*
 * Trim a SimBurstTable.
 */

static int keep_this_injection(SimBurstTable *injection, struct options_t options)
{
	if (options.minCentralfreqFlag && !(injection->freq > options.minCentralfreq))
		return(FALSE);
	if (options.maxCentralfreqFlag && !(injection->freq < options.maxCentralfreq))
		return(FALSE);
	if (options.playground && !isPlayground(injection->l_peak_time.gpsSeconds, injection->l_peak_time.gpsSeconds))
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
	SimBurstTable *head, *prev;

	while(injection && !keep_this_injection(injection, options))
		injection = free_this_injection(injection);
	head = injection;

	/* FIXME: don't check the first event again */
	for(prev = injection; injection; injection = prev->next) {
		if(keep_this_injection(injection, options))
			prev = injection;
		else
			prev->next = free_this_injection(injection);
	}

	return(head);
}


/*
 * Extract the injections that lie between the given start and end times.
 */

static SimBurstTable *extract_injections(LALStatus *stat, SimBurstTable *injection, INT8 start, INT8 end)
{
	SimBurstTable *head = NULL;
	SimBurstTable **addpoint = &head;
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

	return(head);
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

static INT4 read_search_summary_start_end(char *filename, INT4 *start, INT4 *end, FILE *fpout)
{
	SearchSummaryTable *searchSummary = NULL;
	SearchSummaryTable *tmp;
	INT4 local_start, local_end;

	/* allow for NULL pointers if the calling code doesn't care about these
	 * results */
	if(!start)
		start = &local_start;
	if(!end)
		end = &local_end;

	SearchSummaryTableFromLIGOLw(&searchSummary, filename);

	*start = searchSummary->in_start_time.gpsSeconds;
	*end = searchSummary->in_end_time.gpsSeconds;

	if(fpout)
		fprintf(fpout, "%d  %d  %d\n", *start, *end, *end - *start);

	while(searchSummary) {
		tmp = searchSummary;
		searchSummary = searchSummary->next;
		LALFree(tmp);
	}

	return(*end - *start);
}


/*
 * Read trigger list, and generate a new injection list from the injections
 * that occur in the data that was analyzed.
 */

static SnglBurstTable *read_trigger_list(LALStatus *stat, char *filename, INT4 *timeAnalyzed, SimBurstTable **injection, struct options_t options)
{
	FILE *infile;
	char line[MAXSTR];
	SnglBurstTable *list = NULL;
	SnglBurstTable **eventaddpoint = &list;
	SimBurstTable *newinjection = NULL;
	SimBurstTable **injectionaddpoint = &newinjection;
	INT4 start, end;

	if (!(infile = fopen(filename, "r")))
		LALPrintError("Could not open input file\n");

	*timeAnalyzed = 0;
	while(getline(line, MAXSTR, infile)) {
		if(options.verbose)
			fprintf(stderr, "Working on file %s\n", line);

		*timeAnalyzed += read_search_summary_start_end(line, &start, &end, NULL);

		*injectionaddpoint = extract_injections(stat, *injection, start + 1000000000LL, end * 1000000000LL);
		while (*injectionaddpoint)
			injectionaddpoint = &(*injectionaddpoint)->next;

		LAL_CALL(LALSnglBurstTableFromLIGOLw(stat, eventaddpoint, line), stat);
		while (*eventaddpoint)
			eventaddpoint = &(*eventaddpoint)->next;
	}

	while(*injection) {
		SimBurstTable *tmp = *injection;
		*injection = (*injection)->next;
		LALFree(tmp);
	}
	*injection = newinjection;

	fclose(infile);

	return(list);
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

	if (options.playground && !(isPlayground(event->start_time.gpsSeconds, event->start_time.gpsSeconds)))
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
	SnglBurstTable *head, *prev;

	while(event && !keep_this_event(event, options))
		event = event->next;
	head = event;

/* FIXME: don't check the first event again */
	for(prev = event; event; event = prev->next) {
		if(keep_this_event(event, options))
			prev = event;
		else
			prev->next = free_this_event(event);
	}

	return(head);
}


/*
 * =============================================================================
 *                                Entry Point
 * =============================================================================
 */

int main(int argc, char **argv)
{
	static LALStatus stat;

	INT4 sort = FALSE;
	INT4 ndetected;
	INT4 ninjected;
	INT4 timeAnalyzed;
	const long S2StartTime = 729273613;	/* Feb 14 2003 16:00:00 UTC */
	const long S2StopTime = 734367613;	/* Apr 14 2003 15:00:00 UTC */

	static struct options_t options;

	/* filenames */
	CHAR *inputFile = NULL;
	CHAR *injectionFile = NULL;
	CHAR *injmadeFile = NULL;
	CHAR *injFoundFile = NULL;
	CHAR *outSnglFile = NULL;

	/* times of comparison */
	INT4 gpsStartTime = S2StartTime;
	INT4 gpsEndTime = S2StopTime;
	INT8 injPeakTime = 0;
	INT8 burstStartTime = 0;

	/* triggers */
	SnglBurstTable *event;
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
	int c;
	while (1) {
		/* getopt arguments */
		static struct option long_options[] = {
			/* these options set a flag */
			{"verbose", no_argument, &options.verbose, 1},
			/* parameters which determine the output xml file */
			{"input", required_argument, NULL, 'a'},
			{"injfile", required_argument, NULL, 'b'},
			{"injmadefile", required_argument, NULL, 'c'},
			{"max-confidence", required_argument, NULL, 'd'},
			{"gps-start-time", required_argument, NULL, 'e'},
			{"gps-end-time", required_argument, NULL, 'f'},
			{"min-centralfreq", required_argument, NULL, 'g'},
			{"max-centralfreq", required_argument, NULL, 'h'},
			{"injfoundfile", required_argument, NULL, 'j'},
			{"playground", no_argument, &options.playground, 1},
			{"noplayground", no_argument, &options.noplayground, 1},
			{"help", no_argument, NULL, 'o'},
			{"sort", no_argument, NULL, 'p'},
			{"detsnglfile", required_argument, NULL, 'q'},
			{NULL, 0, NULL, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long(argc, argv, "a:c:d:e:f:g:h:i:", long_options, &option_index);

		/* detect the end of the options */
		if (c == -1)
			break;

		switch (c) {
		case 0:
			/* if this option set a flag, do nothing else now */
			if (long_options[option_index].flag != 0)
				break;
			else {
				fprintf(stderr, "error parsing option %s with argument %s\n", long_options[option_index].name, optarg);
				exit(1);
			}
			break;

		case 'a':
			inputFile = optarg;
			break;

		case 'b':
			injectionFile = optarg;
			break;

		case 'c':
			injmadeFile = optarg;
			break;

		case 'd':
			/* the confidence must be smaller than this number */
			options.maxConfidenceFlag = 1;
			options.maxConfidence = atof(optarg);
			break;

		case 'e':
			/* only events with duration greater than this are selected */
			gpsStartTime = atoi(optarg);
			break;

		case 'f':
			/* only events with duration less than this are selected */
			gpsEndTime = atoi(optarg);
			break;

		case 'g':
			/* only events with centralfreq greater than this are selected */
			options.minCentralfreqFlag = 1;
			options.minCentralfreq = atof(optarg);
			break;

		case 'h':
			/* only events with centralfreq less than this are selected */
			options.maxCentralfreqFlag = 1;
			options.maxCentralfreq = atof(optarg);
			break;

		case 'j':
			injFoundFile = optarg;
			break;

		case 'o':
			/* print help */
			LALPrintError(USAGE, *argv);
			exit(BINJ_FIND_EARG);

		case 'p':
			/* sort the events in time */
			sort = TRUE;
			break;


		case 'q':
			outSnglFile = optarg;
			break;

		default:
			exit(BINJ_FIND_EARG);
		}
	}

	if (optind < argc) {
		fprintf(stderr, "extraneous command line arguments:\n");
		while (optind < argc)
			fprintf(stderr, "%s\n", argv[optind++]);
		exit(1);
	}

	if (!injmadeFile || !inputFile || !injectionFile || !outSnglFile) {
		LALPrintError("Input file, injection file, output trig. file and output inj. file " "            names must be specified\n");
		exit(BINJ_FIND_EARG);
	}

	/*
	 * Read and trim the injection list.
	 */

	simBurstList = read_injection_list(&stat, injectionFile, gpsStartTime, gpsEndTime, options);
	simBurstList = trim_injection_list(simBurstList, options);

	/*
	 * Read the trigger list;  remove injections from the injection list
	 * that lie outside the time intervals that were actually analyzed
	 * according to the search summary tables.
	 */

	burstEventList = read_trigger_list(&stat, inputFile, &timeAnalyzed, &simBurstList, options);

	/*
	 * Do any requested cuts and sort the remaining triggers.
	 */

	burstEventList = trim_event_list(burstEventList, options);
	LAL_CALL(LALSortSnglBurst(&stat, &burstEventList, LALCompareSnglBurstByTime), &stat);

	/*
	 * For each injection, search the entire trigger list for the best
	 * match (if any).
	 */

	ninjected = ndetected = 0;
	for (injection = simBurstList; injection; injection = injection->next) {
		ninjected++;

		/* convert injection time to INT8 */
		LAL_CALL(LALGPStoINT8(&stat, &injPeakTime, &injection->l_peak_time), &stat);

		/* loop over the burst events */
		for (event = burstEventList; event; event = event->next) {
			SnglBurstAccuracy accParams;

			LAL_CALL(LALGPStoINT8(&stat, &burstStartTime, &event->start_time), &stat);

			if (injPeakTime < burstStartTime)
				break;

			LAL_CALL(LALCompareSimBurstAndSnglBurst(&stat, injection, event, &accParams), &stat);

			if (!accParams.match)
				continue;

			ndetected++;

			/* record the detected triggers */
			*detTriggersAddPoint = LALCalloc(1, sizeof(**detTriggersAddPoint));
			**detTriggersAddPoint = *event;
			detTriggersAddPoint = &(*detTriggersAddPoint)->next;
			*detTriggersAddPoint = NULL;

			/* record the detected injections */
			*detInjectionsAddPoint = LALCalloc(1, sizeof(**detInjectionsAddPoint));
			**detInjectionsAddPoint = *injection;
			detInjectionsAddPoint = &(*detInjectionsAddPoint)->next;
			*detInjectionsAddPoint = NULL;
			break;
		}
	}

	fprintf(stdout,"%d sec = %d hours analyzed\n", timeAnalyzed, timeAnalyzed/3600);
	fprintf(stdout, "Detected %i injections out of %i made\n", ndetected, ninjected);
	fprintf(stdout, "Efficiency is %f \n", ((REAL4) ndetected / (REAL4) ninjected));

	/*
	 * Write output xml file.
	 */

	memset(&xmlStream, 0, sizeof(LIGOLwXMLStream));
	xmlStream.fp = NULL;

	/* List of injections that were actually made */
	LAL_CALL(LALOpenLIGOLwXMLFile(&stat, &xmlStream, injmadeFile), &stat);
	LAL_CALL(LALBeginLIGOLwXMLTable(&stat, &xmlStream, sim_burst_table), &stat);
	myTable.simBurstTable = simBurstList;
	LAL_CALL(LALWriteLIGOLwXMLTable(&stat, &xmlStream, myTable, sim_burst_table), &stat);
	LAL_CALL(LALEndLIGOLwXMLTable(&stat, &xmlStream), &stat);
	LAL_CALL(LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);

	/* List of injections which were detected */
	LAL_CALL(LALOpenLIGOLwXMLFile(&stat, &xmlStream, injFoundFile), &stat);
	LAL_CALL(LALBeginLIGOLwXMLTable(&stat, &xmlStream, sim_burst_table), &stat);
	myTable.simBurstTable = detectedInjections;
	LAL_CALL(LALWriteLIGOLwXMLTable(&stat, &xmlStream, myTable, sim_burst_table), &stat);
	LAL_CALL(LALEndLIGOLwXMLTable(&stat, &xmlStream), &stat);
	LAL_CALL(LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);

	/* List of triggers corresponding to injection */
	LAL_CALL(LALOpenLIGOLwXMLFile(&stat, &xmlStream, outSnglFile), &stat);
	LAL_CALL(LALBeginLIGOLwXMLTable(&stat, &xmlStream, sngl_burst_table), &stat);
	myTable.snglBurstTable = detectedTriggers;
	LAL_CALL(LALWriteLIGOLwXMLTable(&stat, &xmlStream, myTable, sngl_burst_table), &stat);
	LAL_CALL(LALEndLIGOLwXMLTable(&stat, &xmlStream), &stat);
	LAL_CALL(LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);

	exit(0);
}
