#include "power.h"

#define PROGRAM_NAME "power"

typedef struct
{
   INT4    argc;
   CHAR**  argv;
}LALInitSearchParams;

/* declare the parsing function which is at the end of the file */
int snprintf(char *str, size_t size, const  char  *format, ...);
int initializeEPSearch( int argc, char *argv[], void **searchParams,
        MetadataTable  *procparams );

NRCSID( POWERC, "power $Id$");
RCSID( "power $Id$");

#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"

#include <config.h>
#ifndef HAVE_LIBLALFRAME
int main( void )
{
  fputs( "Disabled: LALApps compiled with non-frame-enabled LAL\n", stderr );
  return 77;
}
#else

#define POWERC_NARGS 19 

#define POWERC_ENORM    0
#define POWERC_ESUB     1
#define POWERC_EARG     2
#define POWERC_EVAL     3
#define POWERC_EFILE    4
#define POWERC_EINPUT   5
#define POWERC_EMEM     6
#define POWERC_ENOEPOCH 7

#define POWERC_MSGENORM         "Normal exit"
#define POWERC_MSGESUB          "Subroutine failed"
#define POWERC_MSGEARG          "Error parsing arguments"
#define POWERC_MSGEVAL          "Input argument out of valid range"
#define POWERC_MSGEFILE         "Could not open file"
#define POWERC_MSGEINPUT        "Error reading file"
#define POWERC_MSGEMEM          "Out of memory"
#define POWERC_MSGENOEPOCH      "Epoch must be provided"

#define TRUE       1
#define FALSE      0

/* Usage format string. */
#define USAGE "Usage: %s --npts npoints --nseg nsegments --olap overlap \
    --olapfctr olapfactor --minfbin nfbin --mintbin ntbin --flow flow \
    --delf df --lngth bandwidth --nsigma sigma --alphdef alpha \
    --segdcle nsegs --threshold threshold --etomstr nevents \
    --framedir framedir --channel channel --simtype simtype --srate srate\
    --spectype spectype --window window --epoch sec nsec --numpts npts \
    [--comment comment] [--framecache filename] [--printSpectrum] \
    [--verbose] [--dbglevel lalDebugLevel] [--help]\n"

/* Fill an array with white noise */
static void makeWhiteNoise(
        LALStatus             *status,
        REAL4TimeSeries       *series, 
        REAL4                  noiseAmpl,
        INT4                   seed
        )
{
    long   length  = series->data->length;
    long   index   = 0;
    static RandomParams *rparams = NULL;

    INITSTATUS (status, "makeWhiteNoise", POWERC);
    ATTATCHSTATUSPTR (status);

    /* generate Gaussian white noise with unit variance */
    LALCreateRandomParams(status->statusPtr, &rparams, seed);
    CHECKSTATUSPTR (status);
    LALNormalDeviates(status->statusPtr, series->data, rparams);
    CHECKSTATUSPTR (status);
    LALDestroyRandomParams(status->statusPtr, &rparams);
    CHECKSTATUSPTR (status);

    /* apply gain factor to noise in time series */
    for (index = 0; index < length; index++)
    {
        double sample = series->data->data[index];
        sample *= noiseAmpl;
        series->data->data[index] = sample;
    }

    DETATCHSTATUSPTR (status);
    RETURN (status);
}



/* some global output flags */
INT4               verbose       = FALSE;
INT4               cluster       = FALSE;
INT4               printSpectrum = FALSE;
INT4               printData     = FALSE;
INT4               whiteNoise    = FALSE;   /* insertion of Gaussian white noise */
INT4               sineBurst     = FALSE;   /* insertion of shaped sine burst  */

/* global variables */
FrChanIn   channelIn;               /* channnel information               */
CHAR       site[2];                 /* one character site                 */
CHAR       ifo[3];                  /* two character interferometer       */
CHAR       comment[LIGOMETA_COMMENT_MAX]; /* string in output file name   */
CHAR      *cachefile     = NULL;    /* name of file with frame cache info */
CHAR      *dirname       = NULL;    /* name of directory with frames      */
REAL4      noiseAmpl     = 1.0;     /* gain factor for white noise        */
INT4       seed          = 1;       /* set non-zero to generate noise     */
INT4       numPoints     = 4096;    /* number of samples from frames      */
INT8       gpsStartTimeNS = 0;      /* gps start time in nsec             */
LIGOTimeGPS   epoch;                /* gps start time                     */
INT4       sampleRate    = 2048;    /* sample rate in Hz                  */

/* data conditioning parameters */
CHAR      *calCacheFile  = NULL;    /* name of the calibration cache file */
CHAR      *injectionFile  = NULL;   /* file with list of injections       */

/* parameters for the sine-gaussian injection testing */
REAL4                 sineFreq      = 100.0;   /* nominal frequency of sine burst */
REAL4                 sineOffset    = 1.0;     /* sin-burst center in time series */
REAL4                 sineAmpl      = 1.0;     /* peak amplitude of sine burst    */
REAL4                 sineWidth     = 1.0;     /* width (in sigmas) of sine burst */


int main( int argc, char *argv[])
{
    static LALStatus      stat;
    LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;

    INT4                  inarg         = 1;
    void                 *searchParams  = NULL;
    LALInitSearchParams   initSearchParams;
    EPSearchParams       *params        = NULL;
    FrStream             *stream        = NULL;
    FrCache              *frameCache    = NULL;
    CHAR                  fname[256];
    BOOLEAN               epochSet      = FALSE;

    /* data storage */
    REAL4TimeSeries            series;
    COMPLEX8FrequencySeries    resp;

    /* Burst events */
    SnglBurstTable      *burstEvent    = NULL;
    SnglBurstTable      *nextEvent     = NULL;
    MetadataTable        myTable;
    MetadataTable        procTable;
    MetadataTable        procparams;
    MetadataTable        searchsumm;
    ProcessParamsTable   *this_proc_param;
    LIGOLwXMLStream      xmlStream;

    /*  used in injections */
    COMPLEX8FrequencySeries     *injRespPtr;    
 
    /* units and other things */
    const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};

    
    /* Must allocate memory for the argv in dso code */
    initSearchParams.argv = (CHAR **) LALMalloc( (size_t)(POWERC_NARGS * 
                sizeof(CHAR*)) );
    
    /* Which error handler to use */
    lal_errhandler = LAL_ERR_EXIT;
    set_debug_level( "3" );

    
    /*******************************************************************
    * INITIALIZE EVERYTHING                                            *
    *******************************************************************/

    /* create the process and process params tables */
    procTable.processTable = (ProcessTable *) 
        LALCalloc( 1, sizeof(ProcessTable) );
    LAL_CALL( LALGPSTimeNow ( &stat, &(procTable.processTable->start_time),
                &accuracy ), &stat );
    LAL_CALL( populate_process_table( &stat, procTable.processTable, 
                PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &stat );
    this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
        LALCalloc( 1, sizeof(ProcessParamsTable) );

    /* parse arguments and fill procparams table */
    initializeEPSearch( argc, argv, &searchParams, &procparams);

    params = (EPSearchParams *) searchParams;
    params->printSpectrum = printSpectrum;
    params->cluster = cluster;
    
    /* create the search summary table */
    searchsumm.searchSummaryTable = (SearchSummaryTable *)
        LALCalloc( 1, sizeof(SearchSummaryTable) );

    /* fill the comment, if a user has specified one, or leave it blank */
    if ( ! *comment )
    {
        snprintf( procTable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
    } else {
        snprintf( procTable.processTable->comment, LIGOMETA_COMMENT_MAX,
                "%s", comment );
    }

    /* create and initialize the time series vector */
    series.data = NULL;
    LAL_CALL( LALCreateVector( &stat, &series.data, numPoints), &stat);
    memset( series.data->data, 0, series.data->length*sizeof(REAL4) );
    series.epoch.gpsSeconds     = epoch.gpsSeconds;
    series.epoch.gpsNanoSeconds = epoch.gpsNanoSeconds;
    strcpy(series.name, params->channelName);
    series.deltaT = 1.0/((REAL8) sampleRate);
    series.f0 = 0.0;
    series.sampleUnits = lalADCCountUnit;

   
    /*******************************************************************
    * GET AND CONDITION THE DATA                                       *
    *******************************************************************/
    
    /* only try to load frame if name is specified */
    if (dirname || cachefile)
    {
        if(dirname){
            /* Open frame stream */
            LAL_CALL( LALFrOpen( &stat, &stream, dirname, "*.gwf" ), &stat);
        }
        else if (cachefile){
            /* Open frame cache */
            LAL_CALL( LALFrCacheImport( &stat, &frameCache, cachefile ), &stat);
            LAL_CALL( LALFrCacheOpen( &stat, &stream, frameCache ), &stat);
            LAL_CALL( LALDestroyFrCache( &stat, &frameCache ), &stat );
        }
        /*
         * Determine information about the channel and seek to the
         * right place in the fram files 
         */
        LAL_CALL( LALFrGetREAL4TimeSeries( &stat, &series, &channelIn, stream), &stat);
        series.epoch.gpsSeconds     = epoch.gpsSeconds;
        series.epoch.gpsNanoSeconds = epoch.gpsNanoSeconds;
        LAL_CALL( LALFrSeek(&stat, &(series.epoch), stream), &stat);

        /* get the data */
        LAL_CALL( LALFrGetREAL4TimeSeries( &stat, &series, &channelIn, stream), &stat);

        /* close the frame stream */
        LAL_CALL( LALFrClose( &stat, &stream ), &stat);
    }

    /* populate time series with white noise if specified */
    if (whiteNoise)
    {
        makeWhiteNoise(&stat, &series, noiseAmpl, seed);
    }

    /* write diagnostic info to disk */
    if ( printData ){
      LALPrintTimeSeries( &series, "./timeseries.dat" );
    }

    /* create storage for the response function */
    memset( &resp, 0, sizeof(COMPLEX8FrequencySeries) );
    if ( calCacheFile )
    {
      LAL_CALL( LALCCreateVector( &stat, &(resp.data), numPoints / 2 + 1 ), 
          &stat );

      /* set the parameters of the response to match the data */
      resp.epoch = epoch;
      resp.deltaF = (REAL8) sampleRate / (REAL8) numPoints;
      resp.sampleUnits = strainPerCount;
      strcpy( resp.name, channelIn.name );

      /* generate the response function for the current time */
      if ( verbose ) fprintf( stdout, "generating response at time %d sec %d ns\n",
          resp.epoch.gpsSeconds, resp.epoch.gpsNanoSeconds );
      LAL_CALL( LALExtractFrameResponse( &stat, &resp, calCacheFile, ifo ),
          &stat );
    } 

    /*****************************************************************
     * Add injections into the time series:  UNTESTED
     *****************************************************************/
    if( sineBurst )
    {
      SimBurstTable *injections = NULL;

       injections =  (SimBurstTable *)LALMalloc( sizeof(SimBurstTable) );
    
      /* Fill in the injection structure */
      injections->geocent_peak_time.gpsSeconds = 723345762;
      injections->geocent_peak_time.gpsNanoSeconds = 0;
      injections->longitude                      = 50.0;
      injections->latitude                       = 50.0;
      injections->hrss                           = sineAmpl;
      injections->freq                           = sineFreq;
      injections->tau                            = sineWidth;
      injections->next                           = NULL;
      /* Inject */

      injRespPtr = &resp;
      LAL_CALL( LALBurstInjectSignals( &stat, &series, injections, injRespPtr ), 
          &stat ); 
    }

    /* Finally call condition data */
    LAL_CALL( EPConditionData( &stat, &series, searchParams), &stat);

    /*******************************************************************
    * DO THE SEARCH                                                    *
    *******************************************************************/
    while ( params->currentSegment < params->initParams->numSegments )
    {
        UINT4                tmpDutyCycle=0;
        UINT4                dumCurrentSeg=params->currentSegment;
        SnglBurstTable      *tmpEvent     = NULL;

        /* count the segments to analyze */
        for ( tmpDutyCycle=0 ; tmpDutyCycle < params->initParams->segDutyCycle && 
                dumCurrentSeg < params->initParams->numSegments ; tmpDutyCycle++ )
        {
            dumCurrentSeg++;
        }
                
        /* tell operator how we are doing */
        if (verbose){
            fprintf(stdout,"Analyzing segments %i -- %i\n", params->currentSegment,
                    params->currentSegment + tmpDutyCycle - 1);
        }

        /* This is the main engine of the excess power method */ 
        LAL_CALL( EPSearch (&stat, params, &tmpEvent, tmpDutyCycle), &stat);

        if ( tmpEvent != NULL ){

            /* add events to event list */
            if (burstEvent == NULL){
                nextEvent = burstEvent = 
                    (SnglBurstTable *)LALMalloc( sizeof(SnglBurstTable) );
            } else {
                nextEvent->next = (SnglBurstTable *)LALMalloc( sizeof(SnglBurstTable) );
                nextEvent = nextEvent->next;
            }
            memcpy(nextEvent, tmpEvent, sizeof(SnglBurstTable));

            /* locate end of event list */
            while (nextEvent->next){
                nextEvent = nextEvent->next;
            }

            /* free the head of the temporary list */
            LALFree( tmpEvent );
            tmpEvent = NULL;

        }

        /* increment to the next segment number to be analyzed */
        params->currentSegment += tmpDutyCycle;
    }


    /*******************************************************************
    * OUTPUT THE RESULTS 
    *******************************************************************/
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
    snprintf( fname, sizeof(fname), "%s-%s-POWER-%d-%d.xml",
            site, comment, epoch.gpsSeconds, (INT4)(series.deltaT * numPoints));
    LAL_CALL( LALOpenLIGOLwXMLFile(&stat, &xmlStream, fname), &stat);


    /* write the process table */
    snprintf( procTable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s", ifo );
    LAL_CALL( LALGPSTimeNow ( &stat, &(procTable.processTable->end_time),
                &accuracy ), &stat );
    LAL_CALL( LALBeginLIGOLwXMLTable( &stat, &xmlStream, process_table ), 
            &stat );
    LAL_CALL( LALWriteLIGOLwXMLTable( &stat, &xmlStream, procTable, 
                process_table ), &stat );
    LAL_CALL( LALEndLIGOLwXMLTable ( &stat, &xmlStream ), &stat );
    LALFree( procTable.processTable );

    
    /* write the process params table */
    LAL_CALL( LALBeginLIGOLwXMLTable( &stat, &xmlStream, process_params_table ), 
            &stat );
    LAL_CALL( LALWriteLIGOLwXMLTable( &stat, &xmlStream, procparams, 
                process_params_table ), &stat );
    LAL_CALL( LALEndLIGOLwXMLTable ( &stat, &xmlStream ), &stat );
    while( procparams.processParamsTable )
    {
        this_proc_param = procparams.processParamsTable;
        procparams.processParamsTable = this_proc_param->next;
        LALFree(this_proc_param);
    }

    
    /* write the search summary table */
    LAL_CALL( LALBeginLIGOLwXMLTable( &stat, &xmlStream, search_summary_table ), 
            &stat );
    LAL_CALL( LALWriteLIGOLwXMLTable( &stat, &xmlStream, searchsumm, 
                search_summary_table ), &stat );
    LAL_CALL( LALEndLIGOLwXMLTable ( &stat, &xmlStream ), &stat );
    LALFree( searchsumm.searchSummaryTable );


    /* Write the results to the burst table */
    LAL_CALL( LALBeginLIGOLwXMLTable (&stat, &xmlStream, sngl_burst_table), &stat);
    myTable.snglBurstTable=burstEvent;
    LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, myTable, 
                sngl_burst_table), &stat);
    while (burstEvent)
    {
        SnglBurstTable *thisEvent;
        thisEvent = burstEvent;
        burstEvent = burstEvent->next;
        LALFree( thisEvent );
    }
    LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);

    /* close the xml stream */
    LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);


    /*******************************************************************
    * FINALIZE EVERYTHING                                            *
    *******************************************************************/
    LAL_CALL( LALSDestroyVector( &stat, &(series.data) ), &stat);
    LAL_CALL( EPFinalizeSearch( &stat, &searchParams), &stat);
    if ( calCacheFile )
    {
      LAL_CALL( LALCDestroyVector( &stat, &(resp.data) ), &stat);
    }
    if ( cachefile ) LALFree( cachefile ); 
    if ( dirname ) LALFree( dirname ); 

    LALCheckMemoryLeaks();

    return 0;
}



/*************************************************************************
 * FUNCTION TO INITIALIZE AND PARSE ARGUMENTS
 *************************************************************************/

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
this_proc_param->next = (ProcessParamsTable *) \
  LALCalloc( 1, sizeof(ProcessParamsTable) ); \
this_proc_param = this_proc_param->next ; \
  snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );\
  snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
      PROGRAM_NAME ); \
      snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
          long_options[option_index].name ); \
          snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); 

int initializeEPSearch( 
        int            argc, 
        char          *argv[], 
        void         **searchParams,
        MetadataTable  *procparams 
        )
{ 
    /* getopt arguments */
    struct option long_options[] =
    {
        /* these options set a flag */
        {"verbose",                 no_argument,       &verbose,           TRUE },
        {"cluster",                 no_argument,       &cluster,           TRUE },
        /* these options don't set a flag */
        {"alphdef",                 required_argument, 0,                 'a'}, 
        {"channel",                 required_argument, 0,                 'b'}, 
        {"comment",                 required_argument, 0,                 'c'}, 
        {"delf",                    required_argument, 0,                 'd'}, 
        {"etomstr",                 required_argument, 0,                 'e'}, 
        {"flow",                    required_argument, 0,                 'f'}, 
        {"framecache",              required_argument, 0,                 'g'}, 
        {"framedir",                required_argument, 0,                 'i'}, 
        {"help",                    no_argument,       0,                 'h'}, 
        {"lngth",                   required_argument, 0,                 'j'}, 
        {"minfbin",                 required_argument, 0,                 'k'}, 
        {"mintbin",                 required_argument, 0,                 'l'}, 
        {"npts",                    required_argument, 0,                 'm'},
        {"noiseamp",                required_argument, 0,                 'n'}, 
        {"numpts",                  required_argument, 0,                 'o'}, 
        {"nseg",                    required_argument, 0,                 'p'}, 
        {"nsigma",                  required_argument, 0,                 'q'}, 
        {"olap",                    required_argument, 0,                 'r'}, 
        {"olapfctr",                required_argument, 0,                 's'}, 
        {"segdcle",                 required_argument, 0,                 't'}, 
        {"simtype",                 required_argument, 0,                 'u'}, 
        {"spectype",                required_argument, 0,                 'v'}, 
        {"start_time",              required_argument, 0,                 'x'}, 
        {"start_time_ns",           required_argument, 0,                 'y'}, 
        {"srate",                   required_argument, 0,                 'z'}, 
        {"sinFreq",                 required_argument, 0,                 'A'}, 
        {"seed",                    required_argument, 0,                 'E'}, 
        {"threshold",               required_argument, 0,                 'B'}, 
        {"window",                  required_argument, 0,                 'C'}, 
        {"dbglevel",                required_argument, 0,                 'D'},
	{"sinOffset",               required_argument, 0,                 'F'},
	{"sinAmpl",                 required_argument, 0,                 'G'},
	{"sinWidth",                required_argument, 0,                 'H'},
        {"calcache",                required_argument, 0,                 'I'},
        /* output options */
        {"printData",               no_argument,       &printData,         1 },
        {"printSpectrum",           no_argument,       &printData,         1 },
        {0, 0, 0, 0}
    };
    int c;
    size_t len=0;
    ProcessParamsTable *this_proc_param = procparams->processParamsTable;
    EPSearchParams     *params;
    LALStatus           stat = blank_status;

    /*
     *
     * allocate memory 
     *
     */

    *searchParams = params = LALMalloc (sizeof( EPSearchParams )); 
    if ( !params )
    {
        fprintf(stderr, "Memory allocation failed for searchParams\n");
        exit(1);
    }

    params->tfTilingInput = 
        (CreateTFTilingIn *) LALMalloc (sizeof(CreateTFTilingIn));
    if ( !params->tfTilingInput )
    {
        LALFree (*searchParams); *searchParams = NULL;
        fprintf(stderr, "Memory allocation failed for tfTilingInput\n");
        exit(1);
    }

    params->initParams = (EPInitParams *) LALMalloc (sizeof(EPInitParams));
    if ( !params-> initParams)
    {
        LALFree (params->tfTilingInput); params->tfTilingInput = NULL;
        LALFree (*searchParams); *searchParams = NULL;
        fprintf(stderr, "Memory allocation failed for initParams\n");
        exit(1);
    }

    params->compEPInput = 
        (ComputeExcessPowerIn *) LALMalloc (sizeof(ComputeExcessPowerIn));
    if ( !params-> compEPInput)
    {
        LALFree (params->initParams); params->initParams = NULL;
        LALFree (params->tfTilingInput); params->tfTilingInput = NULL;
        LALFree (*searchParams); *searchParams = NULL;
        fprintf(stderr, "Memory allocation failed for compEPInput\n");
        exit(1);
    }

    while ( 1 )
    {
        /* getopt_long stores long option here */
        int option_index = 0;

        c = getopt_long_only( argc, argv, 
                "a:b:c:d:e:f:g:i:h:j:k:l:m:n:o:p:q:r:s:t:u:v:x:y:z:A:E:B:C:D:F:G:H:I:",
                long_options, &option_index );

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
                /* default value of alpha used in search code */
                {
                    REAL4 alphaDflt = atof( optarg );
                    if ( alphaDflt <= 0.0 || alphaDflt >= 1.0 ){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be in range [0,1]\n"
                                "(%f specified)\n",
                                long_options[option_index].name, alphaDflt );
                        exit( 1 );
                    }
                    /* default alpha value for tiles with sigma < numSigmaMin */
                    params->compEPInput->alphaDefault     = alphaDflt;                    
                    ADD_PROCESS_PARAM( "float", "%e", alphaDflt );
                }
                break;

            case 'b':
                /* channel to be used in the analysis */
                params->channelName = (CHAR *) LALMalloc(strlen( optarg )+1 );
                if (! params->channelName ){
                    fprintf(stderr,"Error allocating memory for channel name\n");
                }
                strcpy( params->channelName, optarg );
                channelIn.name = params->channelName;
                channelIn.type = ADCDataChannel;

                /* copy the first character to site and the first two to ifo */
                memset( site, 0, sizeof(site) );
                memset( ifo, 0, sizeof(ifo) );
                memcpy( site, channelIn.name, sizeof(site) - 1 );
                memcpy( ifo, channelIn.name, sizeof(ifo) - 1 );
                ADD_PROCESS_PARAM( "string", "%s", optarg );
                break;

            case 'c':
                /* comment string to be used in output file name */
                snprintf( comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
                ADD_PROCESS_PARAM( "string", "%s", optarg );
                break;

            case 'd':
                /* Frequency resolution of first TF plane */
                {
                    REAL4 delf = atof( optarg );
                    if ( delf <= 0.0 ){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be >= 0 (%f specified)\n",
                                long_options[option_index].name, delf );
                        exit( 1 );
                    }
                    /* default alpha value for tiles with sigma < numSigmaMin */
                    params->tfTilingInput->deltaF = delf;
                    ADD_PROCESS_PARAM( "float", "%e", delf );
                }
                break;

            case 'e':
                /* the number of communicated events (integer) */
                {
                    INT4 events = atoi( optarg );
                    if ( events < 1 || events > 999 ){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be in range [1:999] (%i specified)\n",
                                long_options[option_index].name, events );
                        exit( 1 );
                    }
                    /* EK - Max. number of events to communicate to master */
                    params->events2Master = events;
                    ADD_PROCESS_PARAM( "int", "%d", events );
                }
                break;

            case 'f':
                /* Lowest frequency in Hz to be searched */
                {
                    REAL4 flow = atof( optarg );
                    if ( flow < 0.0 ){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be >= 0.0 (%f specified)\n",
                                long_options[option_index].name, flow );
                        exit( 1 );
                    }
                    params->tfTilingInput->flow = flow;
                    ADD_PROCESS_PARAM( "float", "%e", flow );
                }
                break;

            case 'g':
                /* the frame cache file name */
                {
                    len = strlen(optarg) + 1;
                    cachefile = (CHAR *) LALCalloc( len , sizeof(CHAR));
                    memcpy( cachefile, optarg, len);
                    ADD_PROCESS_PARAM( "string", "%s", optarg);
                }
                break;

            case 'i':
                /* directory containing frames */
                {
                    len = strlen(optarg) + 1;
                    dirname =  (CHAR *) LALCalloc( len , sizeof(CHAR));
                    memcpy( dirname, optarg, len);
                    ADD_PROCESS_PARAM( "string", "%s", optarg);
                }
                break;

            case 'h':
                /* print out short help information */
                fprintf(stderr,  USAGE, *argv );
                exit( 0 );
                break;

            case 'j':
                /* frequency bandwidth */
                {
                    INT4 tmplength = atoi(optarg);
                    if (tmplength <= 0){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%i specified)\n",
                                long_options[option_index].name, tmplength);
                        exit( 1 );
                    }
                    params->tfTilingInput->length = tmplength;
                    ADD_PROCESS_PARAM( "int", "%d", tmplength );                    
                }
                break;

            case 'k':
                /* minimum number of freq bins allowed in a search */
                {
                    INT4 minfbin = atoi(optarg);
                    if (minfbin <= 0){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%i specified)\n",
                                long_options[option_index].name, minfbin);
                        exit( 1 );
                    }
                    params->tfTilingInput->minFreqBins = minfbin;
                    ADD_PROCESS_PARAM( "int", "%d", minfbin );                    
                }
                break;

            case 'l':
                /* minimum number of time bins allowed in a search */
                {
                    INT4 mintbin = atoi(optarg);
                    if (mintbin <= 0){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%i specified)\n",
                                long_options[option_index].name, mintbin);
                        exit( 1 );
                    }
                    params->tfTilingInput->minTimeBins = mintbin;
                    ADD_PROCESS_PARAM( "int", "%d", mintbin );                    
                }
                break;


            case 'm':
                /* minimum number of time bins allowed in a search */
                {
                    INT4 tmpm = atoi(optarg);
                    if (tmpm <= 0){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%i specified)\n",
                                long_options[option_index].name, tmpm);
                        exit( 1 );
                    }
                    params->initParams->numPoints = tmpm;
                    params->ntotT = params->initParams->numPoints;      
                    ADD_PROCESS_PARAM( "int", "%d", tmpm );                    
                }
                break;

            case 'n':
                /* turn on white noise simulation at this amplitude */
                {
                    REAL4 tmpamp = atof(optarg);
                    whiteNoise = TRUE;
                    if (tmpamp <= 0.0){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%f specified)\n",
                                long_options[option_index].name, tmpamp);
                        exit( 1 );
                    }
                    noiseAmpl = tmpamp;
                    ADD_PROCESS_PARAM( "float", "%e", tmpamp );                    
                }
                break;

            case 'o':
                /* number of data points to read in */
                {
                    INT4 tmpn = atoi(optarg);
                    if (tmpn <= 0){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%i specified)\n",
                                long_options[option_index].name, tmpn);
                        exit( 1 );
                    }
                    numPoints = tmpn;
                    ADD_PROCESS_PARAM( "int", "%d", tmpn );                    
                }
                break;

            case 'p':
                /* number of data points to read in */
                {
                    INT4 tmpseg = atoi(optarg);
                    if (tmpseg <= 0){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%i specified)\n",
                                long_options[option_index].name, tmpseg);
                        exit( 1 );
                    }
                    params->initParams->numSegments = tmpseg;
                    ADD_PROCESS_PARAM( "int", "%d", tmpseg );                    
                }
                break;

            case 'q':
                /* number of sigma below which we ignore tile */
                {
                    REAL4 tmpsigma = atof(optarg);
                    if (tmpsigma <= 1.0){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%f specified)\n",
                                long_options[option_index].name, tmpsigma);
                        exit( 1 );
                    }
                    params->compEPInput->numSigmaMin  = tmpsigma;
                    ADD_PROCESS_PARAM( "float", "%e", tmpsigma );  
                }
                break;

            case 'r':
                /* Overlap betweeen segments (# of points) */
                {
                    INT4 tmpolap = atoi(optarg);
                    if (tmpolap < 0){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%i specified)\n",
                                long_options[option_index].name, tmpolap);
                        exit( 1 );
                    }
                    params->ovrlap = tmpolap;
                    ADD_PROCESS_PARAM( "int", "%d", tmpolap ); 
                }
                break;

            case 's':
                /* Amount of overlap between neighboring TF tiles */
                {
                    INT4 tmpolap = atoi(optarg);
                    if (tmpolap < 0){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%i specified)\n",
                                long_options[option_index].name, tmpolap);
                        exit( 1 );
                    }
                    params->tfTilingInput->overlapFactor = tmpolap;
                    ADD_PROCESS_PARAM( "int", "%d", tmpolap ); 
                }
                break;


            case 't':
                /* Number of segments sent to slave */
                {
                    INT4 tmpseg = atoi(optarg);
                    if (tmpseg < 0){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%i specified)\n",
                                long_options[option_index].name, tmpseg);
                        exit( 1 );
                    }
                    params->initParams->segDutyCycle      = tmpseg; 
                    ADD_PROCESS_PARAM( "int", "%d", tmpseg ); 
                }
                break;


            case 'u':
                /* Simulation type:  currently ignored. */
                {
                    INT4 tmpsim = atoi(optarg);
                    if (tmpsim < 0 || tmpsim > 3){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%i specified)\n",
                                long_options[option_index].name, tmpsim);
                        exit( 1 );
                    }
                    tmpsim = 0;
                    params->simType =  tmpsim;
                    ADD_PROCESS_PARAM( "int", "%d", tmpsim ); 
                }
                break;

            case 'v':
                /* Spectrum method to use */
                if ( !strcmp( optarg, "useMean" ) ) {
                    params->initParams->method            = useMean;
                } 
                else if ( !strcmp( optarg, "useMedian" ) ) {
                    params->initParams->method            = useMedian;
                }
                else {

                    fprintf(stderr,"invalid argument to --%s:\n"
                            "Must be useMean/useMedian (%s specified)\n",
                            long_options[option_index].name, optarg);
                    exit( 1 );
                }
                ADD_PROCESS_PARAM( "string", "%s", optarg);
                break;


            case 'x':
                {
                    long int gstartt = atol( optarg );
                    if ( gstartt < 441417609 )
                    {
                        fprintf( stderr, "invalid argument to --%s:\n"
                                "GPS start time is prior to " 
                                "Jan 01, 1994  00:00:00 UTC:\n"
                                "(%ld specified)\n",
                                long_options[option_index].name, gstartt );
                        exit( 1 );
                    }
                    if ( gstartt > 999999999 )
                    {
                        fprintf( stderr, "invalid argument to --%s:\n"
                                "GPS start time is after " 
                                "Sep 14, 2011  01:46:26 UTC:\n"
                                "(%ld specified)\n", 
                                long_options[option_index].name, gstartt );
                        exit( 1 );
                    }
                    gpsStartTimeNS += (UINT8) gstartt * 1000000000LL;
                    ADD_PROCESS_PARAM( "int", "%ld", gstartt );
                }
                break;

            case 'y':
                {
                    long int gstarttns = atol( optarg );
                    if ( gstarttns < 0 )
                    {
                        fprintf( stderr, "invalid argument to --%s:\n"
                                "GPS start time nanoseconds is negative\n",
                                long_options[option_index].name );
                        exit( 1 );
                    }
                    if ( gstarttns > 999999999 )
                    {
                        fprintf( stderr, "invalid argument to --%s:\n"
                                "GPS start time nanoseconds is greater than unity:\n" 
                                "Must be <= 999999999 (%ld specified)\n", 
                                long_options[option_index].name, gstarttns );
                        exit( 1 );
                    }
                    gpsStartTimeNS += (UINT8) gstarttns;
                    ADD_PROCESS_PARAM( "int", "%ld", gstarttns );
                }
                break;

            case 'z':
                /* sample rate in Hz */
                {
                    INT4 tmpsrate = atoi(optarg);
                    if (tmpsrate <= 0){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%i specified)\n",
                                long_options[option_index].name, tmpsrate);
                        exit( 1 );
                    }
                    sampleRate = tmpsrate; 
                    ADD_PROCESS_PARAM( "int", "%d", tmpsrate ); 
                }
                break;

	    case 'A':
	        /* inject Sine-Gaussians: set the freq. */
	        {
                  REAL4 tmpsineFreq = atof(optarg);
		  sineBurst = TRUE;
		  if (tmpsineFreq <= 0){
		      fprintf(stderr,"invalid argument to --%s:\n",
                              long_options[option_index].name);
		      exit( 1 );
		  }
		  sineFreq = tmpsineFreq;
		  ADD_PROCESS_PARAM( "float", "%e", tmpsineFreq );
		}      
		break;

            case 'E':
                /* seed for noise simulations */
                {
                    INT4 tmpseed = atoi(optarg);
                    if (tmpseed <= 0){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%i specified)\n",
                                long_options[option_index].name, tmpseed);
                        exit( 1 );
                    }
                    seed = tmpseed; 
                    ADD_PROCESS_PARAM( "int", "%d", tmpseed ); 
                }
                break;

            case 'B':
                /* Identify events with alpha less that this value */
                {
                    REAL4 tmpth = atof(optarg);
                    if (tmpth < 0.0){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%f specified)\n",
                                long_options[option_index].name, tmpth);
                        exit( 1 );
                    }
                    params->alphaThreshold = tmpth;              
                    ADD_PROCESS_PARAM( "float", "%e", tmpth );  
                }
                break;

            case 'C':
                /* Window to use on the data */
                {
                    INT4 tmpwin = atoi(optarg);
                    if (tmpwin < 0 || tmpwin > 6){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%i specified)\n",
                                long_options[option_index].name, tmpwin);
                        exit( 1 );
                    }
                    params->winParams.type                = tmpwin;
                    ADD_PROCESS_PARAM( "int", "%d", tmpwin ); 
                }
                break;

            case 'D':
                /* set the lalDebugLevel to something useful */
                {
                    set_debug_level( optarg );
                    ADD_PROCESS_PARAM( "int", "%d", atoi(optarg) ); 
                }
                break;

            case 'F':
	        /* set the offset for the injections */
	        {
                  REAL4 tmpsineOffset = atof(optarg);
		  if (tmpsineOffset <= 0){
		      fprintf(stderr,"invalid argument to --%s:\n",
                              long_options[option_index].name);
		      exit( 1 );
		  }
		  sineOffset = tmpsineOffset;
		} 
		break;

            case 'G':
	        /* set the amplitude of injections. */
	        {
                  REAL4 tmpsineAmpl = atof(optarg);
		  if (tmpsineAmpl <= 0){
		      fprintf(stderr,"invalid argument to --%s:\n",
                              long_options[option_index].name);
		      exit( 1 );
		  }
		  sineAmpl = tmpsineAmpl;
		} 
		break;

            case 'H':
	        /* set the Width of injections. */
	        {
                  REAL4 tmpsineWidth = atof(optarg);
		  if (tmpsineWidth <= 0){
		      fprintf(stderr,"invalid argument to --%s:\n",
                              long_options[option_index].name);
		      exit( 1 );
		  }
		  sineWidth = tmpsineWidth;
                }
                break;

            case 'I':
                /* create storage for the calibration frame cache name */
                len = strlen( optarg ) + 1;
                calCacheFile = (CHAR *) calloc( len, sizeof(CHAR));
                memcpy( calCacheFile, optarg, len );
                ADD_PROCESS_PARAM( "string", "%s", optarg );
                break;

             case 'J':
                /* create storage for the calibration frame cache name */
                len = strlen( optarg ) + 1;
                injectionFile = (CHAR *) calloc( len, sizeof(CHAR));
                memcpy( injectionFile, optarg, len );
                ADD_PROCESS_PARAM( "string", "%s", optarg );
                break;

           case '?':
                exit( 1 );
                break;

            default:
                fprintf( stderr, "unknown error while parsing options\n" );
                exit( 1 );
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

    /* check that start time was specified & fill the epoch structure */
    if ( !gpsStartTimeNS ){
        fprintf(stderr,"--start-time must be specified");
        exit (1);
    }
    LAL_CALL( LALINT8toGPS( &stat, &epoch, &gpsStartTimeNS), &stat);

    /* initialize parameter structures */
    params->tfTiling     = NULL;
    params->epSegVec     = NULL;
    params->numSlaves    = NULL;

    /* allocate memory for the conditioned segments */
    LAL_CALL( LALCreateEPDataSegmentVector (&stat, &(params->epSegVec), 
                params->initParams), &stat);  

    /* initialize parameters */
    params->haveData        = 0;
    params->currentSegment  = 0;
    params->numEvents       = 0;
    params->searchMaster    = 0;
    params->tfTilingInput->maxTileBand = 64.0;

    return 0;
}

#undef ADD_PROCESS_PARAMS

#endif
