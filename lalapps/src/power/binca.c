#include <math.h>
#include <processtable.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataUtils.h>
#include <lalapps.h>

int snprintf(char *str, size_t size, const char *format, ...);
long long int llabs(long long int j);

RCSID("$Id$");
#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source: /usr/local/cvs/lscsoft/lalapps/src/power/binca.c,v"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "binca"

#define BINCA_EARG   1
#define BINCA_EROW   2
#define BINCA_EFILE  3

#define BINCA_MSGEARG   "Error parsing arguments"
#define BINCA_MSGROW    "Error reading row from XML table"
#define BINCA_MSGEFILE  "Could not open file"
#define BINCA_TIMEWARNING "only triggers before 00:00 Dec 31, \
2010 will be used"

#define BURCA_SEARCHSUMWARNING "The Search summary info of \
of the two ifos do not match "

#define TRUE     1
#define FALSE    0 
#define MAXIFO   16
#define MAXFILES 128
#define MSEC   (1000000LL)
#define NSEC   (1000000000LL)
#define MAXSTR 2048

/* cluster options */
enum { undefined, clusterbypeaktimeandfreq, clusterbytimeandfreq } burstclusterchoice = undefined;

/* Usage format string. */
#define USAGE "Usage: %s --burst-trig-file input-bursttrigfile \
--inspiral-trig-file input-inspiraltrigfile \
--start-time startCoincidence --stop-time stopCoincidence \
--dt deltat --slide-time sec.s to slide \
--slide-time-ns nsec.s to slide --number-slides no. of slides \
--output-burst-dir outputburstdir --output-insp-dir outputinspdir \n"

/* ======================================================================
 *          Gets a line from the input file
 * ======================================================================
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

/* ======================================================================
 *  Adds the process param table
 * ======================================================================
 */

static ProcessParamsTable **add_process_param(ProcessParamsTable **proc_param, const char *type, const char *param, const char *value)
{
  *proc_param = LALCalloc(1, sizeof(**proc_param));
  (*proc_param)->next = NULL;
  snprintf((*proc_param)->program, LIGOMETA_PROGRAM_MAX, PROGRAM_NAME);
  snprintf((*proc_param)->type, LIGOMETA_TYPE_MAX, type);
  snprintf((*proc_param)->param, LIGOMETA_PARAM_MAX, "--%s", param);
  snprintf((*proc_param)->value, LIGOMETA_VALUE_MAX, value);

  return(&(*proc_param)->next);
}

#define ADD_PROCESS_PARAM(type) \
	do { paramaddpoint = add_process_param(paramaddpoint, type, long_options[option_index].name, optarg); } while(0)

/* 
 * ==========================================================================
 *         Command Line Options
 * ==========================================================================
 */

struct options_t {
  CHAR *comment;
  CHAR *outputburstdir;
  CHAR *outputinspdir;
  CHAR *inputBurstFiles;
  CHAR *inputInspiralFiles;
  CHAR *ifoCut;
  CHAR *inspiralclusterchoice;

  INT4  combineSnr;
  
  INT4  usePlayground;
  INT4  ignorePlayground;
  INT4  verbose;
  INT4  cluster_burst;
  INT4  cluster_inspiral;
  INT4  writeInspiral;
  INT4  freeinspmemory;

  INT4  cluster_inspiral_dt;
  INT4  deltaT;
  INT4  startCoincidence;
  INT4  endCoincidence;

  REAL4 slideData;
  INT8  numSlides;
};

/* Sets the default values of command line options */

static void set_option_defaults(struct options_t *options)
{
  const long StartTime = 729273613;	/* Feb 14 2003 16:00:00 UTC */
  const long StopTime = 888888888;	/* Mar 07 2008 01:34:34 UTC */

  options->comment = "";
  options->outputburstdir = "./";
  options->outputinspdir = "./";
  options->inputBurstFiles = NULL;
  options->inputInspiralFiles = NULL;
  options->ifoCut = NULL;
  options->inspiralclusterchoice = NULL;

  options->combineSnr = FALSE;

  options->usePlayground = TRUE;
  options->ignorePlayground = FALSE;
  options->verbose = FALSE;
  options->cluster_burst = FALSE;
  options->cluster_inspiral = FALSE;
  options->writeInspiral = FALSE;
  options->freeinspmemory = FALSE;

  options->cluster_inspiral_dt = 0;
  options->deltaT = 0;
  options->startCoincidence = StartTime;
  options->endCoincidence = StopTime;

  options->slideData = 0.0;
  options->numSlides = 0;
}

/* Parses the command line options */

static void parse_command_line(int argc, char *argv[], struct options_t *options,  MetadataTable *procparams )
{
  ProcessParamsTable **paramaddpoint = &procparams->processParamsTable;
  int option_index;
  struct option long_options[] = {
    /* these options set a flag */
    {"verbose",                 no_argument, &options->verbose,          TRUE},
    {"noplayground",            no_argument, &options->usePlayground,    FALSE},
    {"ignore-playground",       no_argument, &options->ignorePlayground, TRUE},
    {"write-inspiral",          no_argument, &options->writeInspiral,    TRUE},
    {"free-inspmemory",         no_argument, &options->freeinspmemory,   TRUE},
    {"combine-snr",             no_argument, &options->combineSnr,       TRUE},
    /* sets the values */
    {"burst-list-files",        required_argument, NULL, 'a'},
    {"inspiral-list-files",     required_argument, NULL, 'b'},
    {"deltat",                  required_argument, NULL, 'c'},
    {"start-time",              required_argument, NULL, 'd'},
    {"stop-time",               required_argument, NULL, 'e'},
    {"slide-time",              required_argument, NULL, 'f'},
    {"number-slides",           required_argument, NULL, 'g'},
    {"output-burst-dir",        required_argument, NULL, 'h'},
    {"output-insp-dir",         required_argument, NULL, 'i'},
    {"burst-clusterchoice",     required_argument, NULL, 'j'},
    {"user-tag",                required_argument, NULL, 'k'},
    {"ifo-cut",                 required_argument, NULL, 'l'},
    {"inspiral-clusterchoice",  required_argument, NULL, 'm'},
    {"inspiral-cluster-dt",     required_argument, NULL, 'n'},
    {"help",                    no_argument,       NULL, 'o'}, 
    {NULL, 0, NULL, 0}
  };
  int c;

  do switch(c = getopt_long( argc, argv, "a:b:c:d:e:f:g:h:i:j:k:l:m:n:o", long_options, &option_index )) {
  case 'a':
    /* input burst trigger file */
    options->inputBurstFiles = optarg;
    ADD_PROCESS_PARAM("string");
    break;

  case 'b':
    /* input inspiral trigger file */
    options->inputInspiralFiles = optarg;
    ADD_PROCESS_PARAM("string");
    break;

  case 'c':
    /* coincidence window: specified as msec.s and then gets 
     * converted to nanoseconds.
     */
    options->deltaT = atoi(optarg) * MSEC;
    ADD_PROCESS_PARAM("int");
    break;

  case 'd':
    /* start time to check coincidence */
    options->startCoincidence = atoi(optarg);
    ADD_PROCESS_PARAM("int");
    break;

  case 'e':
    /* end time to check coincidence */
    options->endCoincidence = atoi(optarg);
    ADD_PROCESS_PARAM("int");
    break;

  case 'f':
    /* no. of seconds the triggers to be slided*/
    options->slideData = atof(optarg);
    ADD_PROCESS_PARAM("double");
    break;

  case 'g':
    /* no. of slides */
    options->numSlides = atoi(optarg);
    ADD_PROCESS_PARAM("int");      
    break;

  case 'h':
    /* output directory for binca burst files */
    options->outputburstdir = optarg;
    ADD_PROCESS_PARAM("string");
    break;

  case 'i':
    /* output directory for binca insp files */
    options->outputinspdir = optarg;
    ADD_PROCESS_PARAM("string");
    break;
	  
  case 'j':
    /*
     * set the burst cluster option
     */			
    {
      if ( ! strcmp( "clusterbypeaktimeandfreq", optarg ) )
	{
	  options->cluster_burst = TRUE;
	  burstclusterchoice = clusterbypeaktimeandfreq;
	}
      else if ( ! strcmp( "clusterbytimeandfreq", optarg ) )
	{
	  options->cluster_burst = TRUE;
	  burstclusterchoice = clusterbytimeandfreq;
	}
      else
	{
	  fprintf( stderr, "invalid argument to --burst-clusterchoice\n"
		   "unknown clusterchoice specified;\n"
		   "(must be one of:clusterbypeaktimeandfreq ,clusterbytimeandfreq )\n");
	}
    }
    ADD_PROCESS_PARAM("string");
    break;

  case 'k':
    /* sets the user-tag */
    options->comment = optarg;
    ADD_PROCESS_PARAM("string");
    break;

  case 'l':
    /* ifo to be cut out from the inspiral table */
    options->ifoCut = optarg;
    ADD_PROCESS_PARAM("string");
    break;

  case 'm':
    /* sets the inspiral cluster choice */
    {
      if ( ! strcmp( "snr", optarg ) )
	{
	  options->cluster_inspiral = TRUE;
	  options->inspiralclusterchoice = optarg;
	}
      else
	{
	  fprintf( stderr, "invalid argument to --inspiral-clusterchoice\n"
		   "unknown clusterchoice specified;\n"
		   "(must be one of:snr )\n");
	}
    }
    ADD_PROCESS_PARAM("string");
    break;
    break;

  case 'n':
    /* inspiral cluster time: specified as msec.s and then gets 
     * converted to nanoseconds.
     */
    options->cluster_inspiral_dt = atoi(optarg) * MSEC;
    ADD_PROCESS_PARAM("int");
    break;

  case 'o':
    fprintf(stderr, USAGE, argv[0]);
    exit(1);

  case 0:
    break;

  case -1:
    break;

  case '?':
    fprintf( stderr, USAGE , argv[0]);
    exit(1);

  case ':':
    fprintf(stderr, USAGE, argv[0]);
    exit(1);
      
  default:
    fprintf( stderr, "unknown error while parsing options\n" );
    fprintf( stderr, USAGE, argv[0] );
    exit(1);
  } while(c!= -1);
}

/* 
 * ==========================================================================
 * Funtion to read in the inspiral triggers. All the triggers are read
 * in from the files listed in the input .txt file to a linked list of 
 * SnglInspirals
 * ==========================================================================
 */
 
static SnglInspiralTable *read_inspiraltriggers(LALStatus *stat, char *filename, struct options_t options)
{
  FILE *infile;
  char line[MAXSTR];
  SnglInspiralTable *list = NULL;
  SnglInspiralTable **addpoint = &list;
  SnglInspiralClusterChoice clusterchoice = snr;

  if(options.verbose)
    fprintf(stderr, "Reading in SnglInspiral Table\n");

  if(!(infile = fopen(filename, "r")))
    LALPrintError("Could not open input inspiral file\n");

  while(getline(line, MAXSTR, infile)) 
    {
      if (options.verbose)
	fprintf(stderr, "Working on file %s\n", line);
      
      LALSnglInspiralTableFromLIGOLw(addpoint, line, 0, -1);
      
      if(*addpoint)
	{
	  LAL_CALL( LALClusterSnglInspiralTable(stat, *addpoint, 0, clusterchoice ), stat);
	  if(options.ifoCut)
	    LAL_CALL( LALIfoCutSingleInspiral(stat, &(*addpoint), options.ifoCut ), stat);
	}

      while(*addpoint != NULL)
	addpoint = &(**addpoint).next;
    }

  fclose(infile);

  if(list)
    LAL_CALL( LALSortSnglInspiral( stat, &list, *LALCompareSnglInspiralByTime ), stat );

  return(list);
}

/*
 * =========================================================================
 * Reads in the search summary table from the input burst file as in 
 * *filename
 * =========================================================================
 */

static INT8 read_search_summary_start_end(LALStatus *stat, char *filename, INT8 *start, INT8 *end, FILE *fpout, char *ifo)
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
  snprintf(ifo, 16, "%s", searchSummary->ifos);
  LAL_CALL(LALGPStoINT8(stat, start, &searchSummary->out_start_time), stat);
  LAL_CALL(LALGPStoINT8(stat, end, &searchSummary->out_end_time), stat);

  if(fpout)
    fprintf(fpout, "%lld  %lld  %lld\n", *start, *end, *end - *start);

  while(searchSummary) 
    {
      SearchSummaryTable *tmp = searchSummary;
      searchSummary = searchSummary->next;
      LALFree(tmp);
    }

  return(*end - *start);
}

/* 
 * ======================================================================
 * Function to adjust search times 
 * ======================================================================
 */
static INT4 adjust_times( INT8 start, INT8 end, struct options_t options)
{ 
  INT8 optstart = options.startCoincidence * NSEC;
  INT8 optend = options.endCoincidence * NSEC;

 if( (end < optstart) || (start > optend) )
   return 0;
 else if ((start < optstart) && (end > optstart) && (end < optend) )
   {
     start = optstart;
     return 1;
   }
 else if ((end > optend) && (start > optstart) && (start < optend))
   { 
     end = optend;
     return 1;
   }
 else if ((start < optstart) && (end > optend))
   {
     start = optstart;
     end = optend;
     return 1;
   }
 else
   return 1;
}

/*
 * =========================================================================
 * Select the inspiral trigger with the highest snr
 * =========================================================================
 */

static SnglInspiralTable *select_inspiral(LALStatus *stat, SnglInspiralTable *a, SnglInspiralTable *b, struct options_t options)
{
  if(!a)
    return(b);
  if(!b)
    return(a);

  return(b->snr < a->snr ? a : b);
}

/*
 * ==========================================================================
 * Function to find the coincident burst and inspiral triggers 
 * ==========================================================================
 */

static void find_coincident_events(LALStatus *stat, SnglInspiralTable *inspiraltriggers, SnglBurstTable *bursttriggers, SnglInspiralTable **coincinsp, SnglBurstTable **coincburst, INT4 *ncoincident, INT4 *nbcoincident, INT8 start, INT8 end, struct options_t options)
{
  SnglInspiralTable *tmpInspiralEvent, *bestsnr;

  /*****************************************************************
   * First sort the burst and inspiral triggers
   *****************************************************************/
  LAL_CALL( LALSortSnglBurst(stat, &bursttriggers, XLALCompareSnglBurstByStartTime ), stat);

  LAL_CALL( LALSortSnglInspiral( stat, &inspiraltriggers, *LALCompareSnglInspiralByTime ), stat );

  /*****************************************************************
   * find the first burst trigger after search star
   *****************************************************************/

  if (options.verbose) 
    fprintf(stderr,"Moving to first burst trigger in window\n");
  
  while ( (bursttriggers != NULL) && (bursttriggers->start_time.gpsSeconds * NSEC < start) )
    {
      bursttriggers = bursttriggers->next;
    }
	
  /*****************************************************************
   * find the first inspiral trigger after search start
   *****************************************************************/

  if (options.verbose) 
    fprintf(stderr,"Moving to first inspiral trigger in window\n");
       
  while ( (inspiraltriggers != NULL) && (inspiraltriggers->end_time.gpsSeconds * NSEC < start) )
    {
      inspiraltriggers = inspiraltriggers->next;
    }

  /*****************************************************************
   * outer loop over burst triggers
   ****************************************************************/

  if (options.verbose) 
    fprintf(stderr,"Start loop over burst triggers\n");

  while ( (bursttriggers != NULL) && (bursttriggers->start_time.gpsSeconds * NSEC < end) )
    {
      INT8 burststart, burstend, burstpeak, inspiralend;
      INT4 isPlay = 0;
      LIGOTimeGPS burstendgps={0,0};
      INT4 match = 0;

      /* Convert the start time of the burst trigger to nanoseconds(ns.) */ 
      LAL_CALL( LALGPStoINT8(stat, &burststart, &(bursttriggers->start_time)), stat);
      /* Covert the peak time into nanoseconds */
      LAL_CALL( LALGPStoINT8(stat, &burstpeak, &(bursttriggers->peak_time)), stat);

      /* Find the gps end time of the burst trigger and convert to ns. */
      LAL_CALL(LALAddFloatToGPS(stat, &burstendgps, &(bursttriggers->start_time), bursttriggers->duration), stat);
      LAL_CALL( LALGPStoINT8(stat, &burstend, &(burstendgps)), stat);

      /* Catch up inspiral triggers: Find the first inspiral trigger
       * which has its end time within deltat ns. of the start time of
       * the burst trigger
       */
      while (inspiraltriggers != NULL)
	{
	  LAL_CALL( LALGPStoINT8(stat, &inspiralend, &(inspiraltriggers->end_time)), stat);
	  if (inspiralend > burstpeak - options.deltaT)
	    {
	      break;
	    }
	  inspiraltriggers = inspiraltriggers->next;
	}

      /* Check to see if the burst trigger is in playgrounnd 
       * (if asked to)
       */
      if (!options.ignorePlayground)
	{ 
	  LAL_CALL( LALINT8NanoSecIsPlayground( stat, &isPlay, &burststart ), stat );

	  if ( options.verbose )
	    {
	      if ( isPlay )
		{
		  fprintf( stderr, "trigger is playground\n" );
		} 
	      else
		{
		  fprintf( stderr, "trigger is not playground\n" );
		}
	    }
	}

      /* if we are playground only and the trigger is in playground or
       * we are not using playground and the trigger is not in the 
       * playground or if we want to ignore playground check ... */

      if ( (options.usePlayground && isPlay) || (!options.usePlayground && !isPlay) || options.ignorePlayground )
	{
	  tmpInspiralEvent = inspiraltriggers;

	  /* Looping over the inspiral triggers till the inspiral end time
	   * crosses the end time of the burst trigger
	   */
	  if (options.verbose) 
	    fprintf(stderr,"Start loop over inspiral triggers\n");

	  bestsnr = NULL;
	  while (tmpInspiralEvent != NULL)
	    {
	      LAL_CALL( LALGPStoINT8(stat, &inspiralend, &(tmpInspiralEvent->end_time)), stat);

	      /* As long as the inspiral end time is less than the
	       * burst peak time + deltat keep on considering them as 
	       * coincident triggers 
	       */
	      if (inspiralend > burstpeak + options.deltaT)
		break;


	      /* compare the trigger for the maximum snr */
	      bestsnr = select_inspiral(stat, bestsnr, tmpInspiralEvent, options);

	      /* Move on to the next inspiral trigger */  
	      tmpInspiralEvent = tmpInspiralEvent->next;
	    }

	      /* record the coincident inspiral trigger: 
	       * As long as the inspiral trigger end time
	       * is within deltaT ns. of the start time
	       * of the burst trigger to deltaT after the 
	       * end of the burst trigger, we consider that to 
	       * be coincident.  
	       */
	  if(bestsnr)
	    {
	      (*ncoincident)++;
	      match = 1;
	      *coincinsp = LALMalloc(sizeof(**coincinsp));
	      **coincinsp = *bestsnr;
	      coincinsp = &(*coincinsp)->next;
	      *coincinsp = NULL;
	    }
	}

      /* If there was a coincident inspiral trigger 
       * record the burst trigger as a coincident one 
       */

      if (match)
	{
	  REAL4 tmpsnr = 0.0;

	  (*nbcoincident)++;

	  if (options.combineSnr)
	    {
	      tmpsnr = bursttriggers->snr;
	      bursttriggers->snr = bursttriggers->snr + pow(bestsnr->snr,2);
	    }

	  *coincburst = LALMalloc(sizeof(**coincburst));
	  **coincburst = *bursttriggers;
	  coincburst = &(*coincburst)->next;
	  *coincburst = NULL;

	  if (options.combineSnr)
	    bursttriggers->snr = tmpsnr;
	}

      /* Move on to the next burst trigger */
      bursttriggers = bursttriggers->next;
    }
}

/*
 * ===================================================================
 * Functions to free the trigger lists 
 * ===================================================================
 */

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

static SnglInspiralTable *free_this_inspiral_event(SnglInspiralTable *event)
{
  SnglInspiralTable *next = event ? event->next : NULL;
  LALFree(event);
  return(next);
}

static void free_inspiral_events(SnglInspiralTable *list)
{
  while(list)
    list = free_this_inspiral_event(list);
}



/*
 * ==========================================================================
 *
 * The main program
 *
 * ==========================================================================
 */

int main(int argc, char **argv)
{
  static LALStatus       stat;
  struct options_t       options;
  LALLeapSecAccuracy     accuracy = LALLEAPSEC_LOOSE;
  SnglInspiralClusterChoice clusterchoice = snr;

  INT8                   timeAnalyzed = 0;
  INT8                   SearchStart;
  INT8                   SearchEnd;
  INT4                   ncoincident;
  INT4                   nbcoincident;

  CHAR                   outfileName[512];
  CHAR                   line[MAXSTR];
  CHAR                   ifo[3];

  SnglBurstTable         *burst_trigger_list;
  SnglBurstTable         *coincidentBurstEvents=NULL;
  SnglBurstTable         **coinc_burstaddpoint=&coincidentBurstEvents;

  SnglInspiralTable      *inspiral_trigger_list=NULL;
  SnglInspiralTable      *inspiral_list=NULL;
  SnglInspiralTable     **inspiral_trimmed_trigger_list=&inspiral_list;
  SnglInspiralTable      *coincidentInspiralEvents=NULL;
  SnglInspiralTable     **coinc_inspiraladdpoint=&coincidentInspiralEvents;

  MetadataTable           myTable;
  MetadataTable           procTable;
  MetadataTable           procparams;
  MetadataTable           searchsumm;
  LIGOLwXMLStream         xmlStream;

  FILE *infile;

  /*==================================================================
   * Initialize things
   *=================================================================*/
  memset(&stat, 0, sizeof(stat));
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "2" );
  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  set_option_defaults(&options);

  /* create the process and process params tables */
  procTable.processTable = LALCalloc(1, sizeof(ProcessTable));
  LAL_CALL( LALGPSTimeNow(&stat, &(procTable.processTable->start_time), &accuracy ), &stat);
  /*LAL_CALL( populate_process_table( &stat, procTable.processTable, PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &stat );*/
  procparams.processParamsTable = NULL;

  /* Parse the command line arguments */
  parse_command_line(argc, argv, &options, &procparams);

  /* fill the comment */
  snprintf(procTable.processTable->comment, LIGOMETA_COMMENT_MAX, "%s", options.comment);
 
  /* ==================================================================
   * End Initializing things
   * ==================================================================
   */
  
  /* =======================================================================
   * Read in the inspiral triggers (All the files in input .txt file) 
   * =======================================================================   
   */
 
  if (options.inputInspiralFiles)
    inspiral_trigger_list = read_inspiraltriggers(&stat, options.inputInspiralFiles, options);

  if(!inspiral_trigger_list && options.verbose)
    fprintf(stderr,"Warning: No inspiral triggers within the specified times\n");

  /* ======================================================================
   * Read in the burst triggers one file at a time 
   * ======================================================================
   */

  if(!(infile = fopen(options.inputBurstFiles, "r")))
    LALPrintError("Could not open input burst file\n");

  while(getline(line, MAXSTR, infile))
    {
      burst_trigger_list = NULL;
      ncoincident = 0;
      nbcoincident = 0;

      if (options.verbose)
	fprintf(stderr, "Working on burst file %s\n", line);

      /*
       * ================================================================= 
       * Determine the times encompassed by this burst file 
       * =================================================================
       */
      read_search_summary_start_end(&stat, line, &SearchStart, &SearchEnd, NULL, ifo);

      /* Check if the burst and inspiral ifo.s match */
      if(options.ifoCut && (strcmp(ifo, options.ifoCut) != 0 ))
	{
	  fprintf(stderr,"Error:burst trigger ifo and inspiral trigger ifo do not match\n");
	  exit(1);
	}

      /* ================================================================
       * Check the times with the requested startcoincidence 
       * stopcoincidence and adjust if required to
       * ================================================================
       */

      if(!adjust_times(SearchStart, SearchEnd, options))
	continue;

      timeAnalyzed += SearchEnd - SearchStart;
     
      /* =================================================================
       * Populate the search summary table, to be written in the 
       * output file. The searchsummary start and end times are  
       * same as in the input burst file
       * =================================================================
       */
      searchsumm.searchSummaryTable = LALCalloc(1, sizeof(SearchSummaryTable));
      snprintf(searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX, "%s", options.comment);
      searchsumm.searchSummaryTable->nnodes = 1;
      XLALINT8NSToGPS( &(searchsumm.searchSummaryTable->out_start_time), SearchStart );
      XLALINT8NSToGPS( &(searchsumm.searchSummaryTable->out_end_time), SearchEnd );
      XLALINT8NSToGPS( &(searchsumm.searchSummaryTable->in_start_time), SearchStart );
      XLALINT8NSToGPS( &(searchsumm.searchSummaryTable->in_end_time), SearchEnd );

      /* =================================================================
       * Read the burst triggers from this file 
       * =================================================================
       */
      LAL_CALL(LALSnglBurstTableFromLIGOLw(&stat, &burst_trigger_list, line), &stat); 

      if(!burst_trigger_list && options.verbose)
	fprintf(stderr," No burst trigger in this file, moving on to the next\n");

      /* 
       * =================================================================
       * If there are burst and inspiral triggers then check for 
       * coincidence
       * =================================================================
       */
      if(burst_trigger_list && inspiral_trigger_list)
	{
	  /*
	   * Find the coincident inspiral and burst triggers 
	   */
	  find_coincident_events(&stat, inspiral_trigger_list, burst_trigger_list, coinc_inspiraladdpoint, coinc_burstaddpoint, &ncoincident, &nbcoincident, SearchStart, SearchEnd, options);

	  if (options.verbose)
	    fprintf(stderr,"%d burst and %d insp triggers have been found in times %lld - %lld\n", nbcoincident, ncoincident, SearchStart, SearchEnd);

	  /* 
	   * Sort the triggers before output
	   */
	  if(coincidentInspiralEvents)
	    LAL_CALL( LALSortSnglInspiral( &stat, &coincidentInspiralEvents, *LALCompareSnglInspiralByTime ), &stat );

	  if(coincidentBurstEvents)
	    LAL_CALL( LALSortSnglBurst(&stat, &coincidentBurstEvents, XLALCompareSnglBurstByStartTime ), &stat);

	  /* 
	   * Cluster the burst triggers if asked for 
	   */
	  if(options.cluster_burst && burstclusterchoice == clusterbypeaktimeandfreq && coincidentBurstEvents)
	    XLALClusterSnglBurstTable(&coincidentBurstEvents, XLALCompareSnglBurstByPeakTime, XLALCompareSnglBurstByPeakTimeAndFreq, XLALSnglBurstCluster);
	  else if (options.cluster_burst && burstclusterchoice == clusterbytimeandfreq && coincidentBurstEvents)
	    XLALClusterSnglBurstTable(&coincidentBurstEvents,  NULL, XLALCompareSnglBurst, XLALSnglBurstCluster);
	}

      /* ================================================================
       *   Write the xml output 
       * ================================================================
       */

      snprintf(outfileName, sizeof(outfileName)-1, "%s/%s-BINCA_%s-%d-%d.xml", options.outputburstdir,ifo,options.comment,searchsumm.searchSummaryTable->in_start_time.gpsSeconds,searchsumm.searchSummaryTable->in_end_time.gpsSeconds - searchsumm.searchSummaryTable->in_start_time.gpsSeconds);
      outfileName[sizeof(outfileName)-1] = '\0';

      LAL_CALL( LALOpenLIGOLwXMLFile(&stat, &xmlStream, outfileName), &stat);
   
      /* write process table */
      LAL_CALL(LALGPSTimeNow(&stat, &(procTable.processTable->end_time), &accuracy), &stat);
      LAL_CALL(LALBeginLIGOLwXMLTable(&stat, &xmlStream, process_table), &stat);
      LAL_CALL(LALWriteLIGOLwXMLTable(&stat, &xmlStream, procTable, process_table), &stat);
      LAL_CALL(LALEndLIGOLwXMLTable(&stat, &xmlStream), &stat);
   
      /* write process params table */
      LAL_CALL(LALBeginLIGOLwXMLTable(&stat, &xmlStream, process_params_table), &stat);
      LAL_CALL(LALWriteLIGOLwXMLTable(&stat, &xmlStream, procparams, process_params_table), &stat);
      LAL_CALL(LALEndLIGOLwXMLTable(&stat, &xmlStream), &stat);
   
      /* write search summary */
      searchsumm.searchSummaryTable->nevents = ncoincident;
      LAL_CALL( LALBeginLIGOLwXMLTable( &stat, &xmlStream, search_summary_table ), &stat );
      LAL_CALL( LALWriteLIGOLwXMLTable( &stat, &xmlStream, searchsumm, search_summary_table ), &stat );
      LAL_CALL( LALEndLIGOLwXMLTable ( &stat, &xmlStream ), &stat );
  
      /*write the coincident burst triggers */
      LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, sngl_burst_table), &stat);
      myTable.snglBurstTable = coincidentBurstEvents;
      LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, myTable, sngl_burst_table), &stat);
      LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);
      

      /*write the coincident inspiral triggers*/      
      LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, sngl_inspiral_table), &stat);
      myTable.snglInspiralTable = coincidentInspiralEvents;
      LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, myTable, sngl_inspiral_table), &stat);
      LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);

      LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);
  
      /* ===============================================================
       * Free the allocated memory
       * ==============================================================
       */
      if (options.verbose)
	fprintf(stderr,"Freeing coincident lists \n");

      while(coincidentBurstEvents)
	{
	  SnglBurstTable *tmp;
	  tmp = coincidentBurstEvents->next;
	  LALFree(coincidentBurstEvents);
	  coincidentBurstEvents = tmp;
	}

      while(coincidentInspiralEvents)
	{
	  SnglInspiralTable *tmp;
	  tmp = coincidentInspiralEvents->next;
	  LALFree(coincidentInspiralEvents);
	  coincidentInspiralEvents = tmp;
	}
	
      /* Free the burst event list and the searchsummary table
       * before moving on to the next burst input file in the
       * .txt file
       */
      free_events(burst_trigger_list);
      LALFree(searchsumm.searchSummaryTable);
    }
    
  /* Close the burst input .txt file */
  fclose(infile);

  LALFree(procTable.processTable);
 
  while(procparams.processParamsTable) 
    {
      ProcessParamsTable *table = procparams.processParamsTable;
      procparams.processParamsTable = table->next;
      LALFree(table);
    }

  /* Free the memory allocated for the inspiral list if
   * asked to. This is not always done since LALFree is 
   * really slow and it takes a long time. 
   */
  if(options.freeinspmemory)
    {
      free_inspiral_events(inspiral_trigger_list);
      LALCheckMemoryLeaks();
    }

  return 0;
}
