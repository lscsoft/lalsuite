#include "power.h"

/* declare the parsing function which is at the end of the file */
int snprintf(char *str, size_t size, const  char  *format, ...);
int initializeEPSearch( int argc, char *argv[], void **searchParams,
        MetadataTable  *procparams );

NRCSID( POWERC, "power $Id$");
RCSID( "power $Id$");

#define PROGRAM_NAME "power"
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
INT4               geodata       = FALSE;
INT4               printSpectrum = FALSE;
INT4               printData     = FALSE;
INT4               whiteNoise    = FALSE;   /* insertion of Gaussian white noise */
INT4               sineBurst     = FALSE;   /* insertion of shaped sine burst  */
INT4               injFlag       = FALSE;
INT4               calFlag       = FALSE;
INT4               mdcFlag       = FALSE;
 
/* global variables */
FrChanIn   channelIn;               /* channnel information               */
FrChanIn   mdcchannelIn;            /* mdc signal only channnel info      */
EPSearchParams  *mdcparams = NULL;  /* mdc search param                   */
CHAR       site[2];                 /* one character site                 */
CHAR       ifo[3];                  /* two character interferometer       */
CHAR       comment[LIGOMETA_COMMENT_MAX]; /* string in output file name   */
CHAR      *cachefile      = NULL;   /* name of file with frame cache info */
CHAR      *dirname        = NULL;   /* name of directory with frames      */
REAL4      noiseAmpl      = 1.0;    /* gain factor for white noise        */
INT4       seed           = 1;      /* set non-zero to generate noise     */
INT4       numPoints      = 4096;   /* number of samples from frames      */
INT4       totalNumPoints = 0;      /* total number of points to analyze  */
UINT4       totalNumSegs   = 0;      /* total number of points to analyze  */
INT8       gpsStartTimeNS = 0;      /* gps start time in nsec             */
INT8       gpsStopTimeNS  = 0;      /* gps start time in nsec             */
LIGOTimeGPS   startEpoch;           /* gps start time                     */
LIGOTimeGPS   stopEpoch;            /* gps stop time                      */
INT4       sampleRate     = 2048;   /* sample rate in Hz                  */

/* data conditioning parameters */
CHAR      *calCacheFile  = NULL;    /* name of the calibration cache file */
CHAR      *injectionFile = NULL;    /* file with list of injections       */
CHAR      *mdcCacheFile  = NULL;    /* name of mdc signal cache file */

/* GEO data high pass corner freq. */
REAL8	  fcorner	= 100.0;	/* corner frequency in Hz */

/* parameters for the sine-gaussian injection testing */
REAL4                 sineFreq      = 100.0;   /* nominal frequency of sine burst */
REAL4                 sineOffset    = 1.0;     /* sin-burst center in time series */
REAL4                 sineAmpl      = 1.0;     /* peak amplitude of sine burst    */
REAL4                 sineWidth     = 1.0;     /* width (in sigmas) of sine burst */


int main( int argc, char *argv[])
{
    static LALStatus      stat;
    LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;

    void                 *searchParams  = NULL;
    EPSearchParams       *params        = NULL;
    FrStream             *stream        = NULL;
    FrCache              *frameCache    = NULL;
    CHAR                  fname[256];
    PassBandParamStruc    highpassParam;
    REAL4                 fsafety=0;
    LIGOTimeGPS		  duration = { 0, 0};
    LIGOTimeGPS		  tmpEpoch = { 0, 0};
    LALTimeInterval tmpInterval;
    REAL8                 tmpOffset = 0.0;
    
    /* data storage */
    REAL4TimeSeries            series;
    REAL8TimeSeries            geoSeries;
    COMPLEX8FrequencySeries    resp;
    REAL4TimeSeries            mdcSeries;

    /* Burst events */
    SnglBurstTable      *burstEvent    = NULL;
    SnglBurstTable      *nextEvent     = NULL;
    MetadataTable        myTable;
    MetadataTable        procTable;
    MetadataTable        procparams;
    MetadataTable        searchsumm;
    ProcessParamsTable   *this_proc_param;
    LIGOLwXMLStream      xmlStream;

    /* units and other things */
    const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};

    
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
      snprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX, " " );    
    } 
    else 
    {
      snprintf( procTable.processTable->comment, LIGOMETA_COMMENT_MAX,
          "%s", comment );
      snprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX,
          "%s", comment );
    }

    /* the number of nodes for a standalone job is always 1 */
    searchsumm.searchSummaryTable->nnodes = 1;

    /* set the temporary time variable indicating start of chunk */
    tmpEpoch.gpsSeconds = startEpoch.gpsSeconds;
    tmpEpoch.gpsNanoSeconds = startEpoch.gpsNanoSeconds;

    /* loop over chunks of data small enough to fit into memory */
    while (totalNumSegs>0){

      /* make sure you don't expect too many segments */
      if (params->initParams->numSegments > totalNumSegs)
        params->initParams->numSegments = totalNumSegs;

      /* decrement the total number of segments */
      totalNumSegs -= params->initParams->numSegments;

      /* compute the number of points in a chunk */
      numPoints = params->initParams->numSegments * (
          params->initParams->numPoints - params->ovrlap )
        + 3 * params->ovrlap;
      
    /* create and initialize the time series vector */
    series.data = NULL;
    LAL_CALL( LALCreateVector( &stat, &series.data, numPoints), &stat);
    memset( series.data->data, 0, series.data->length*sizeof(REAL4) );
    series.epoch.gpsSeconds     = tmpEpoch.gpsSeconds;
    series.epoch.gpsNanoSeconds = tmpEpoch.gpsNanoSeconds;
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
        REAL8 tmpTime=0;

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
        if (geodata){
          INT4 i;

          /* create and initialize the time series vector */
          geoSeries.data = NULL;
          LAL_CALL( LALDCreateVector( &stat, &geoSeries.data, numPoints), &stat);
          memset( geoSeries.data->data, 0, geoSeries.data->length*sizeof(REAL8) );
          geoSeries.epoch.gpsSeconds     = tmpEpoch.gpsSeconds;
          geoSeries.epoch.gpsNanoSeconds = tmpEpoch.gpsNanoSeconds;
          strcpy(geoSeries.name, params->channelName);
          geoSeries.deltaT = 1.0/((REAL8) sampleRate);
          geoSeries.f0 = 0.0;
          geoSeries.sampleUnits = lalADCCountUnit;
          LAL_CALL( LALFrGetREAL8TimeSeries( &stat, &geoSeries, &channelIn, stream), &stat);
          geoSeries.epoch.gpsSeconds     = tmpEpoch.gpsSeconds;
          geoSeries.epoch.gpsNanoSeconds = tmpEpoch.gpsNanoSeconds;
          LAL_CALL( LALFrSeek(&stat, &(geoSeries.epoch), stream), &stat);

          /* get the data */
          LAL_CALL( LALFrGetREAL8TimeSeries( &stat, &geoSeries, &channelIn, stream), &stat);

	  /* high pass filter before casting REAL8 to REAL4 */

	  highpassParam.nMax = 4;
    	  fsafety = params->tfTilingInput->flow - 10.0;
    	  highpassParam.f2 = fsafety > fcorner ? fcorner : fsafety;
    	  highpassParam.f1 = -1.0;
    	  highpassParam.a2 = 0.1;
    	  highpassParam.a1 = -1.0;
    	  LAL_CALL( LALButterworthREAL8TimeSeries(&stat, &geoSeries, &highpassParam), &stat);

          for(i=0;i<numPoints;i++){
            series.data->data[i] = (REAL4) geoSeries.data->data[i];
          }
          series.epoch.gpsSeconds = geoSeries.epoch.gpsSeconds;
          series.epoch.gpsNanoSeconds = geoSeries.epoch.gpsNanoSeconds;
          strcpy(series.name, geoSeries.name);
          series.deltaT = geoSeries.deltaT;
          series.f0 = geoSeries.f0;
          series.sampleUnits = lalADCCountUnit;
	  LAL_CALL( LALDDestroyVector( &stat, &geoSeries.data), &stat);

        }
        else
        {
          LAL_CALL( LALFrGetREAL4TimeSeries( &stat, &series, &channelIn, stream), &stat);
          series.epoch.gpsSeconds     = tmpEpoch.gpsSeconds;
          series.epoch.gpsNanoSeconds = tmpEpoch.gpsNanoSeconds;
          LAL_CALL( LALFrSeek(&stat, &(series.epoch), stream), &stat);

          /* get the data */
          LAL_CALL( LALFrGetREAL4TimeSeries( &stat, &series, &channelIn, stream), &stat);
        }
        
        /* store the start time of the raw channel in the search summary */
        if ( !(searchsumm.searchSummaryTable->in_start_time.gpsSeconds) ) 
          searchsumm.searchSummaryTable->in_start_time = series.epoch;

        /* store the stop time of the raw channel in the search summary */
        if (totalNumSegs<=0){
          LAL_CALL( LALGPStoFloat( &stat, &tmpTime, &(series.epoch) ), 
              &stat );
          tmpTime += series.deltaT * (REAL8) series.data->length;
          LAL_CALL( LALFloatToGPS( &stat, 
                &(searchsumm.searchSummaryTable->in_end_time), &tmpTime ), &stat );
        }

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
      LALPrintTimeSeries( &series, "./timeseriesasq.dat" );
    }

    /* create storage for the response function */
    memset( &resp, 0, sizeof(COMPLEX8FrequencySeries) );
    if ( calCacheFile )
    {
      INT4 i;

      LAL_CALL( LALCCreateVector( &stat, &(resp.data), numPoints / 2 + 1 ), 
          &stat );

      /* set the parameters of the response to match the data */
      resp.epoch = tmpEpoch;
      resp.deltaF = (REAL8) sampleRate / (REAL8) numPoints;
      resp.sampleUnits = strainPerCount;
      strcpy( resp.name, channelIn.name );

      /* generate the response function for the current time */
      if ( verbose ) 
        fprintf( stdout, "generating response at time %d sec %d ns\n",
          resp.epoch.gpsSeconds, resp.epoch.gpsNanoSeconds );

      /* getting the response is handled differently for geo */
      if(geodata){
        for(i=0;i<numPoints / 2 + 1;i++){
          resp.data->data[i].re = 1.0;
          resp.data->data[i].im = 0.0;
        }
      }
      else{
        LAL_CALL( LALExtractFrameResponse( &stat, &resp, calCacheFile, ifo, 
              &duration ), &stat );
      }
    } 

    /*****************************************************************
     * Add injections into the time series:  
     *****************************************************************/
    if( injFlag )
    {
      INT4  startTime = series.epoch.gpsSeconds;
      INT4  stopTime = startTime + (INT4)( series.data->length * series.deltaT );
      SimBurstTable *injections = NULL;

      if ( !calFlag )
      {
        fprintf(stderr, "Must supply calibration information for injectoins\n");
        exit(1);
      }

      /* read in list from file and make the injections */
      if ( verbose )
        fprintf(stdout, "Reading in SimBurst Table\n");

      LAL_CALL( LALSimBurstTableFromLIGOLw ( &stat, &injections, injectionFile,
            startTime, stopTime), &stat );

      if ( verbose )
        fprintf(stdout, "Injecting signals into time series\n");

      LAL_CALL( LALBurstInjectSignals( &stat, &series, injections, &resp ), 
          &stat ); 

      while (injections)
      {
        SimBurstTable *thisEvent;
        thisEvent = injections;
        injections = injections->next;
        LALFree( thisEvent );
      }

      if ( verbose )
        fprintf(stdout, "Finished making the injections\n");

      /* write diagnostic info to disk */
      if ( printData ){
        LALPrintTimeSeries( &series, "./injections.dat" );
      }
    }

    /* if one wants to use mdc signals for injections */

    if (mdcFlag)
    {
      INT4 i;

      /*open mdc cache */
      LAL_CALL( LALFrCacheImport( &stat, &frameCache, mdcCacheFile ), &stat);
      LAL_CALL( LALFrCacheOpen( &stat, &stream, frameCache ), &stat);
      LAL_CALL( LALDestroyFrCache( &stat, &frameCache ), &stat );

      /* create and initialize the mdc time series vector */
      mdcSeries.data = NULL;
      LAL_CALL( LALCreateVector( &stat, &mdcSeries.data, numPoints), &stat);
      memset( mdcSeries.data->data, 0, mdcSeries.data->length*sizeof(REAL4) );
      mdcSeries.epoch.gpsSeconds     = tmpEpoch.gpsSeconds;
      mdcSeries.epoch.gpsNanoSeconds = tmpEpoch.gpsNanoSeconds;
      strcpy(mdcSeries.name, mdcparams->channelName);
      mdcSeries.deltaT = 1.0/((REAL8) sampleRate);
      mdcSeries.f0 = 0.0;
      mdcSeries.sampleUnits = lalADCCountUnit;
      LAL_CALL( LALFrGetREAL4TimeSeries( &stat, &mdcSeries, &mdcchannelIn, stream), &stat);
      mdcSeries.epoch.gpsSeconds     = tmpEpoch.gpsSeconds;
      mdcSeries.epoch.gpsNanoSeconds = tmpEpoch.gpsNanoSeconds;
      LAL_CALL( LALFrSeek(&stat, &(mdcSeries.epoch), stream), &stat);

      /* get the mdc signal data */
      LAL_CALL( LALFrGetREAL4TimeSeries( &stat, &mdcSeries, &mdcchannelIn, stream), &stat);

      /* write diagnostic info to disk */
      if ( printData ){
        LALPrintTimeSeries( &mdcSeries, "./timeseriesmdc.dat" );
      }

      /* add the signal to the As_Q data */

      for(i=0;i<numPoints;i++)
      {
        series.data->data[i] += mdcSeries.data->data[i];
      }

      /* write diagnostic info to disk */
      if ( printData ){
        LALPrintTimeSeries( &series, "./timeseriesasqmdc.dat" );
      }

      /* destroy the mdc data vector */
      LAL_CALL( LALDestroyVector( &stat, &mdcSeries.data), &stat);
      /* close the frame stream */
      LAL_CALL( LALFrClose( &stat, &stream ), &stat);
    }

    /* Finally call condition data */
    LAL_CALL( EPConditionData( &stat, &series, searchParams), &stat);

    /* add information about times to summary table */
    {
      REAL8 tmpTime=0;

      /* store the start and end time of the raw channel in the search summary */
      LAL_CALL( LALGPStoFloat( &stat, &tmpTime, &(series.epoch) ), 
          &stat );
      tmpTime += series.deltaT * (REAL8) params->ovrlap;
      LAL_CALL( LALFloatToGPS( &stat, 
            &(searchsumm.searchSummaryTable->out_start_time), &tmpTime ), &stat );
      tmpTime += series.deltaT * ((REAL8) series.data->length - (REAL8) params->ovrlap);
      LAL_CALL( LALFloatToGPS( &stat, 
            &(searchsumm.searchSummaryTable->out_end_time), &tmpTime ), &stat );
    }
    
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

    /* compute the start time for the next chunk */
    tmpOffset = (REAL8)(numPoints - 3 * params->ovrlap)/((REAL8) sampleRate);
    LAL_CALL( LALFloatToInterval(&stat, &tmpInterval, 
          &tmpOffset), &stat );
    LAL_CALL( LALIncrementGPS(&stat, &(tmpEpoch), &(tmpEpoch), 
          &tmpInterval), &stat );

    /* clean up memory from that run */
    LAL_CALL( LALSDestroyVector( &stat, &(series.data) ), &stat);
    if ( calCacheFile )
    {
      LAL_CALL( LALCDestroyVector( &stat, &(resp.data) ), &stat);
    }
    if ( cachefile ) LALFree( cachefile ); 
    if ( dirname ) LALFree( dirname ); 

    }


    /*******************************************************************
    * OUTPUT THE RESULTS 
    *******************************************************************/
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
    snprintf( fname, sizeof(fname), "%s-%s-POWER-%d-%d.xml",
            ifo, comment, startEpoch.gpsSeconds, 
            stopEpoch.gpsSeconds-startEpoch.gpsSeconds);
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
    LAL_CALL( EPFinalizeSearch( &stat, &searchParams), &stat);

    if ( calCacheFile )
    {
      LAL_CALL( LALCDestroyVector( &stat, &(resp.data) ), &stat);
    }
    if ( cachefile ) LALFree( cachefile ); 
    if ( dirname ) LALFree( dirname ); 
    if ( mdcFlag ) LALFree( mdcparams->channelName );
    if ( mdcFlag ) LALFree( mdcparams );

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
        {"stop_time",               required_argument, 0,                 'X'}, 
        {"stop_time_ns",            required_argument, 0,                 'Y'}, 
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
        {"injfile",                 required_argument, 0,                 'J'},
	/* geo data flag, argument is corner freq. of high pass filter */
	{"geodata",		    required_argument, 0,		  'K'},
	/* mdc data & channel information */
        {"mdccache",                required_argument, 0,                 'L'},
        {"mdcchannel",              required_argument, 0,                 'M'},
        /* output options */
        {"printData",               no_argument,       &printData,         TRUE },
        {"printSpectrum",           no_argument,       &printSpectrum,     TRUE },
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
                "a:b:c:d:e:f:g:i:h:j:k:l:m:n:o:p:q:r:s:t:u:v:x:y:z:A:E:B:C:D:F:G:H:I:J:K:L:M:",
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

            case 'X':
                {
                    long int gstopt = atol( optarg );
                    if ( gstopt < 441417609 )
                    {
                        fprintf( stderr, "invalid argument to --%s:\n"
                                "GPS start time is prior to " 
                                "Jan 01, 1994  00:00:00 UTC:\n"
                                "(%ld specified)\n",
                                long_options[option_index].name, gstopt );
                        exit( 1 );
                    }
                    if ( gstopt > 999999999 )
                    {
                        fprintf( stderr, "invalid argument to --%s:\n"
                                "GPS start time is after " 
                                "Sep 14, 2011  01:46:26 UTC:\n"
                                "(%ld specified)\n", 
                                long_options[option_index].name, gstopt );
                        exit( 1 );
                    }
                    gpsStopTimeNS += (UINT8) gstopt * 1000000000LL;
                    ADD_PROCESS_PARAM( "int", "%ld", gstopt );
                }
                break;

            case 'Y':
                {
                    long int gstoptns = atol( optarg );
                    if ( gstoptns < 0 )
                    {
                        fprintf( stderr, "invalid argument to --%s:\n"
                                "GPS start time nanoseconds is negative\n",
                                long_options[option_index].name );
                        exit( 1 );
                    }
                    if ( gstoptns > 999999999 )
                    {
                        fprintf( stderr, "invalid argument to --%s:\n"
                                "GPS start time nanoseconds is greater than unity:\n" 
                                "Must be <= 999999999 (%ld specified)\n", 
                                long_options[option_index].name, gstoptns );
                        exit( 1 );
                    }
                    gpsStopTimeNS += (UINT8) gstoptns;
                    ADD_PROCESS_PARAM( "int", "%ld", gstoptns );
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
                /* Identify events with alpha less than this value */
                {
                    REAL8 tmpth = atof(optarg);
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
                calFlag = TRUE;
                ADD_PROCESS_PARAM( "string", "%s", optarg );
                break;

             case 'J':
                /* create storage for the calibration frame cache name */
                len = strlen( optarg ) + 1;
                injectionFile = (CHAR *) calloc( len, sizeof(CHAR));
                memcpy( injectionFile, optarg, len );
                injFlag = TRUE;
                ADD_PROCESS_PARAM( "string", "%s", optarg );
                break;

	     case 'K':
	       {
		/* read flag parameter to determine high pass corner freq. */
		REAL8 tmpfcorner = atof(optarg);
                  if (tmpfcorner <= 0){
                      fprintf(stderr,"invalid argument to --%s:\n",
                              long_options[option_index].name);
                      exit( 1 );
                  }
                  fcorner = tmpfcorner;
		  geodata = TRUE;
                }
                break;

             case 'L':
                /* create storage for the mdc signal frame cache name */
                len = strlen( optarg ) + 1;
                mdcCacheFile = (CHAR *) calloc( len, sizeof(CHAR));
                memcpy( mdcCacheFile, optarg, len );
                mdcFlag = TRUE;
                ADD_PROCESS_PARAM( "string", "%s", optarg );
                break;

	     case 'M':
	        mdcparams = LALMalloc (sizeof( EPSearchParams ));
	        if ( !mdcparams )
		  {
		    fprintf(stderr, "Memory allocation failed for searchParams\n");
		    exit(1);
		  }
                /* mdc channel to be used in the analysis */
                mdcparams->channelName = (CHAR *) LALMalloc(strlen( optarg )+1 );
                if (! mdcparams->channelName ){
                    fprintf(stderr,"Error allocating memory for channel name\n");
                }
                strcpy( mdcparams->channelName, optarg );
                mdcchannelIn.name = mdcparams->channelName;
                mdcchannelIn.type = ADCDataChannel;
  
                /* copy the first character to site and the first two to ifo 
                memset( site, 0, sizeof(site) );
                memset( ifo, 0, sizeof(ifo) );
                memcpy( site, mdcchannelIn.name, sizeof(site) - 1 );
                memcpy( ifo, mdcchannelIn.name, sizeof(ifo) - 1 );*/
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
    LAL_CALL( LALINT8toGPS( &stat, &startEpoch, &gpsStartTimeNS), &stat);

    /* check that stop time was specified & fill the epoch structure */
    if ( !gpsStopTimeNS ){
        fprintf(stderr,"--stop-time must be specified");
        exit (1);
    }
    LAL_CALL( LALINT8toGPS( &stat, &stopEpoch, &gpsStartTimeNS), &stat);

    /* compute the total number of points to be analyzed */
    totalNumPoints = (stopEpoch.gpsSeconds - startEpoch.gpsSeconds ) 
      * sampleRate;

    /* compute the total number of segments to be analyzed */
    totalNumSegs = (INT4) (( totalNumPoints - 3 * params->ovrlap ) /
      ( params->initParams->numPoints - params->ovrlap ));

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
