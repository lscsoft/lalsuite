#include "power.h"
#include <lal/TimeSeries.h>
#include <lal/EPSearch.h>

/* declare the parsing function which is at the end of the file */
int snprintf(char *str, size_t size, const  char  *format, ...);
int initializeEPSearch( int argc, char *argv[], EPSearchParams **params,
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

#define TRUE       1
#define FALSE      0

/* Usage format string. */
#define USAGE "Usage: %s [options] \n"\
"\n"\
"  --sample-rate F              filter data at F Hz, downsampling if necessary\n"\
"  --resample-filter TYPE       set resample filter to TYPE (ldas|butterworth)\n"\
"\n"\
    "--npts             npoints \n"\
    "--nseg             nsegments \n" \
    "--olap             overlap \n"\
    "--olapfctr         olapfactor \n"\
    "--minfbin          nfbin \n"\
    "--mintbin          ntbin \n"\
    "--flow             flow \n"\
    "--delf             df \n"\
    "--lngth            bandwidth \n"\
    "--nsigma           sigma \n"\
    "--alphdef          alpha \n"\
    "--segdcle          nsegs \n"\
    "--threshold        threshold\n"\
    "--etomstr          nevents \n"\
    "--framecache       filename \n"\
    "--channel-name     channel \n"\
    "--simtype          simtype \n"\
    "--spectype         spectype \n"\
    "--window           window \n"\
    "--gps-start-time       start-sec \n"\
    "--gps-start-time-ns    start-nsec \n"\
    "--gps-stop-time        stop-sec \n"\
    "--gps-stop-time-ns     stop-nsec \n"\
    "[--injection-file      injection file] \n"\
    "[--calibration-cache   calibration file] \n"\
    "[--mdccache        mdccache filename] \n"\
    "[--mdcchannel      mdcchannel] \n"\
    "[--geodata         high pass corner freq] \n"\
    "[--user-tag        comment] \n"\
    "[--framedir        dirname] \n"\
    "[--printSpectrum] \n"\
    "[--verbose] \n"\
    "[--dbglevel        lalDebugLevel] \n"\
    "[--help]\n"

/* this is a temporary function to make some code below more readable.
 * Obviously the argument and return types are sub-optimal... */
static size_t min(size_t a, size_t b)
{
	return(a < b ? a : b);
}

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
UINT4       totalNumSegs   = 0;     /* total number of segments to analyze  */
INT8       gpsStartTimeNS = 0;      /* gps start time in nsec             */
INT8       gpsStopTimeNS  = 0;      /* gps start time in nsec             */
LIGOTimeGPS   startEpoch;           /* gps start time                     */
LIGOTimeGPS   stopEpoch;            /* gps stop time                      */
INT4       frameSampleRate  = -1;   /* sample rate of the frame data      */
INT4       targetSampleRate = -1;   /* sample rate after resampling  */
INT4       psdAveragePoints = 0;    /* no. of points used to estimate the spec*/ 
INT4       ramLimitPoints = 0;      /* no. of points to be read in at a time */

/* data conditioning parameters */
CHAR      *calCacheFile  = NULL;    /* name of the calibration cache file */
CHAR      *injectionFile = NULL;    /* file with list of injections       */
CHAR      *mdcCacheFile  = NULL;    /* name of mdc signal cache file */
ResampleTSFilter resampFiltType = -1;

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

    EPSearchParams       *params        = NULL;
    FrStream             *stream        = NULL;
    FrCache              *frameCache    = NULL;
    CHAR                  fname[256];
    PassBandParamStruc    highpassParam;
    REAL4                 fsafety=0;
    CalibrationUpdateParams calfacts;
    LIGOTimeGPS		  tmpEpoch = { 0, 0};
    LALTimeInterval tmpInterval;
    REAL8                 tmpOffset = 0.0;
    REAL4                 minFreq = 0.0;    
    int                   start_sample;
    int                   usedNumPoints = 0;

    /* data storage */
    REAL4TimeSeries            series;
    REAL8TimeSeries            geoSeries;
    COMPLEX8FrequencySeries    resp;
    REAL4TimeSeries            mdcSeries;
 
    /* Burst events */
    SnglBurstTable      *burstEvent = NULL;
    SnglBurstTable     **EventAddPoint = &burstEvent;
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
    initializeEPSearch( argc, argv, &params, &procparams);

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

    /*set the min. freq to be searched for */
    minFreq = params->tfTilingInput->flow;

    /* set the temporary time variable indicating start of chunk */
    tmpEpoch = startEpoch;

    /******************************************************************
     * OUTER LOOP over data small enough to fit into memory 
     ******************************************************************/
    while ((totalNumPoints-usedNumPoints)>(2*(int)params->ovrlap*(frameSampleRate/targetSampleRate))){

      /* tell operator how we are doing */
      if (verbose){
        fprintf(stdout,"%i points analysed && %i points left\n",usedNumPoints,totalNumPoints-usedNumPoints);
      }

      /* compute the number of points in a chunk */
      numPoints = min(ramLimitPoints,(totalNumPoints - usedNumPoints)) ;

      if (verbose){
        fprintf(stdout,"read in %i points\n",numPoints);
      }
      
      /* count the no. of points that are being used */
      usedNumPoints += numPoints;

      /* create and initialize the time series vector */
      series.data = NULL;
      LAL_CALL( LALCreateVector( &stat, &series.data, numPoints), &stat);
      memset( series.data->data, 0, series.data->length*sizeof(REAL4) );
      series.epoch = tmpEpoch;
      strcpy(series.name, params->channelName);
      series.f0 = 0.0;
      series.deltaT = 1.0/((REAL8)frameSampleRate);
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
          geoSeries.epoch = tmpEpoch;
          strcpy(geoSeries.name, params->channelName);
          geoSeries.deltaT = 1.0/((REAL8)frameSampleRate);
          geoSeries.f0 = 0.0;
          geoSeries.sampleUnits = lalADCCountUnit;
          LAL_CALL( LALFrGetREAL8TimeSeries( &stat, &geoSeries, &channelIn, stream), &stat);
          geoSeries.epoch = tmpEpoch;
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
          series.epoch = geoSeries.epoch;
          strcpy(series.name, geoSeries.name);
          series.deltaT = geoSeries.deltaT;
          series.f0 = geoSeries.f0;
          series.sampleUnits = lalADCCountUnit;
          LAL_CALL( LALDDestroyVector( &stat, &geoSeries.data), &stat);

        }
        else
        {
          LAL_CALL( LALFrGetREAL4TimeSeries( &stat, &series, &channelIn, stream), &stat);
          series.epoch = tmpEpoch;
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
        if (verbose )
        {
          REAL4 norm=0;
          UINT4 j ;
          /* PRB - The normalization constant */
          norm = 0.0;
          for( j=0 ; j<series.data->length ; j++)
          {
            REAL4 re = series.data->data[j];

            norm += (re*re);
          }
          norm = sqrt(norm/series.data->length);
          fprintf(stderr,"the norm is %e\n",norm);
        }
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
        resp.deltaF = (REAL8) frameSampleRate / (REAL8) numPoints;
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
          memset( &calfacts, 0, sizeof(CalibrationUpdateParams) );
          calfacts.ifo = ifo;

          LAL_CALL( LALExtractFrameResponse( &stat, &resp, calCacheFile,  
                &calfacts ), &stat );
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
          fprintf(stderr, "Must supply calibration information for injections\n");
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

        if ( verbose )
          fprintf(stdout, "Using MDC frames for injections\n");

        /*open mdc cache */
        LAL_CALL( LALFrCacheImport( &stat, &frameCache, mdcCacheFile ), &stat);
        LAL_CALL( LALFrCacheOpen( &stat, &stream, frameCache ), &stat);
        LAL_CALL( LALDestroyFrCache( &stat, &frameCache ), &stat );

        /* create and initialize the mdc time series vector */
        mdcSeries.data = NULL;
        LAL_CALL( LALCreateVector( &stat, &mdcSeries.data, numPoints), &stat);
        memset( mdcSeries.data->data, 0, mdcSeries.data->length*sizeof(REAL4) );
        mdcSeries.epoch = tmpEpoch;
        strcpy(mdcSeries.name, mdcparams->channelName);
        mdcSeries.deltaT = 1.0/((REAL8)frameSampleRate);
        mdcSeries.f0 = 0.0;
        mdcSeries.sampleUnits = lalADCCountUnit;
        LAL_CALL( LALFrGetREAL4TimeSeries( &stat, &mdcSeries, &mdcchannelIn, stream), &stat);
        mdcSeries.epoch = tmpEpoch;
        LAL_CALL( LALFrSeek(&stat, &(mdcSeries.epoch), stream), &stat);

        /* get the mdc signal data */
        LAL_CALL(LALFrGetREAL4TimeSeries( &stat, &mdcSeries, &mdcchannelIn, stream), &stat);


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
      LAL_CALL( EPConditionData( &stat, &series, minFreq, 1.0/targetSampleRate, resampFiltType, params), &stat);

      /* add information about times to summary table */
      {
        REAL8 tmpTime=0;

        /* store the 'actual' start and end time(accounting for the  
         * 0.5 sec's at the begining & end) in the search summary */

        LAL_CALL( LALGPStoFloat( &stat, &tmpTime, &(series.epoch) ), 
		  &stat );
        if ( !(searchsumm.searchSummaryTable->out_start_time.gpsSeconds) )
        {
          LAL_CALL( LALFloatToGPS( &stat, 
                &(searchsumm.searchSummaryTable->out_start_time), &tmpTime ), &stat );
        }
        tmpTime += series.deltaT * ((REAL8) series.data->length);
        if ( totalNumSegs <= 0 )
        {
          LAL_CALL( LALFloatToGPS( &stat, 
                &(searchsumm.searchSummaryTable->out_end_time), &tmpTime ), &stat );
        }
      }

      /*******************************************************************
       * DO THE SEARCH                                                    *
       *******************************************************************/
      if (verbose)
	fprintf(stdout,"Got %i points to analyse after conditioning\n", series.data->length);

      /* first check if the psdAveragePoints of data are available or not?
       * If not, then set the psdAveragePoints to be equal to the length 
       * of data available 
       */
      if ((int)series.data->length < psdAveragePoints)
	psdAveragePoints = series.data->length;

      for(start_sample = 0; start_sample < ((int)series.data->length - 3*params->ovrlap); start_sample += (psdAveragePoints - 3* params->ovrlap))
      {
        REAL4TimeSeries *interval;
	int shifted_start_sample = 0;

	/* if the no. of points left is less than the psdAveragePoints
	 * then move the epoch back so that it can cut out 
	 * psdAveragePoints of data 
	 */
	if ( ((int)series.data->length - start_sample) < psdAveragePoints ){

	 shifted_start_sample = series.data->length - psdAveragePoints;

	 LAL_CALL(LALCutREAL4TimeSeries(&stat, &interval, &series, shifted_start_sample,psdAveragePoints), &stat);
	}
	else
	  LAL_CALL(LALCutREAL4TimeSeries(&stat, &interval, &series, start_sample, psdAveragePoints), &stat);
	
        if (verbose)
          fprintf(stdout,"Analyzing samples %i -- %i\n", start_sample,min(shifted_start_sample,start_sample) + interval->data->length);

        LAL_CALL(EPSearch(&stat, interval, params, EventAddPoint), &stat);

        while(*EventAddPoint)
          EventAddPoint = &(*EventAddPoint)->next;

        LAL_CALL(LALDestroyREAL4TimeSeries(&stat, interval), &stat);
      } 

 
      /* compute the start time for the next chunk */
      tmpOffset = (REAL8)(numPoints - 2 * (params->ovrlap*(frameSampleRate/targetSampleRate)))/((REAL8)frameSampleRate);
      LAL_CALL( LALFloatToInterval(&stat, &tmpInterval, 
            &tmpOffset), &stat );
      LAL_CALL( LALIncrementGPS(&stat, &(tmpEpoch), &(tmpEpoch), 
            &tmpInterval), &stat );

      /*recalculate the used no. of points considering the 
       *fact that we moved back by 2*params->ovrlap
       *when we calculated the offset in the previous step
       */
      usedNumPoints -= 2*params->ovrlap*(frameSampleRate/targetSampleRate);
      
     /* clean up memory from that run */
      LAL_CALL( LALSDestroyVector( &stat, &(series.data) ), &stat);
      if ( calCacheFile )
      {
        LAL_CALL( LALCDestroyVector( &stat, &(resp.data) ), &stat);
      }
    } 
    /*******************************************************************
     * while (totalNumSegs>0){..  ends here
     *******************************************************************/

    if ( cachefile ) LALFree( cachefile ); 
    if ( dirname ) LALFree( dirname ); 
    if ( mdcFlag ) LALFree( mdcparams->channelName );
    if ( mdcFlag ) LALFree( mdcparams );

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
    myTable.snglBurstTable = burstEvent;
    LAL_CALL( LALWriteLIGOLwXMLTable (&stat, &xmlStream, myTable, 
                sngl_burst_table), &stat);

    /* clean up memory allocated to hold burst events */
    while(burstEvent) {
        SnglBurstTable *event = burstEvent;
        burstEvent = burstEvent->next;
        LALFree(event);
    }
    LAL_CALL( LALEndLIGOLwXMLTable (&stat, &xmlStream), &stat);
    
    /* close the xml stream */
    LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);


    /*******************************************************************
    * FINALIZE EVERYTHING                                            *
    *******************************************************************/
    LAL_CALL( EPFinalizeSearch( &stat, &params), &stat);
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
        EPSearchParams **params,
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
        {"default-alpha",           required_argument, 0,                 'a'},
        {"channel-name",            required_argument, 0,                 'b'}, 
        {"user-tag",                required_argument, 0,                 'c'}, 
	{"etomstr",                 required_argument, 0,                 'e'},
        {"event-limit",             required_argument, 0,                 'e'},
        {"flow",                    required_argument, 0,                 'f'}, 
        {"low-freq-cutoff",         required_argument, 0,                 'f'}, 
        {"frame-cache",             required_argument, 0,                 'g'}, 
        {"framedir",                required_argument, 0,                 'i'}, 
        {"help",                    no_argument,       0,                 'h'}, 
        {"lngth",                   required_argument, 0,                 'j'},
	{"bandwidth",               required_argument, 0,                 'j'},
        {"minfbin",                 required_argument, 0,                 'k'}, 
        {"min-freq-bin",            required_argument, 0,                 'k'}, 
        {"mintbin",                 required_argument, 0,                 'l'}, 
        {"min-time-bin",            required_argument, 0,                 'l'},
        {"npts",                    required_argument, 0,                 'm'},
        {"noiseamp",                required_argument, 0,                 'n'},
	{"noise-amplitude",         required_argument, 0,                 'n'}, 
        {"nseg",                    required_argument, 0,                 'p'}, 
        {"ram-limit-points",        required_argument, 0,                 'p'},
        {"nsigma",                  required_argument, 0,                 'q'}, 
        {"olap",                    required_argument, 0,                 'r'},
	{"shift-points",            required_argument, 0,                 'r'},
        {"olapfctr",                required_argument, 0,                 's'}, 
        {"tileoverlap-factor",      required_argument, 0,                 's'},
        {"segdcle",                 required_argument, 0,                 't'},
        {"psd-average-points",      required_argument, 0,                 't'},  
        {"spectype",                required_argument, 0,                 'v'},
	{"psd-average-method",      required_argument, 0,                 'v'},
        {"gps-start-time",          required_argument, 0,                 'x'}, 
        {"gps-start-time-ns",       required_argument, 0,                 'y'}, 
        {"gps-end-time",            required_argument, 0,                 'X'}, 
        {"gps-end-time-ns",         required_argument, 0,                 'Y'}, 
        {"sample-rate",             required_argument, 0,                 'z'},
	{"frame-sample-rate",       required_argument, 0,                 'z'},
        {"target-sample-rate",      required_argument, 0,                 'A'},
	{"resample-filter",         required_argument, 0,                 'N'}, 
        {"seed",                    required_argument, 0,                 'E'}, 
        {"threshold",               required_argument, 0,                 'B'}, 
        {"window",                  required_argument, 0,                 'C'}, 
        {"dbglevel",                required_argument, 0,                 'D'},
        {"debug-level",             required_argument, 0,                 'D'},
        {"calibration-cache",       required_argument, 0,                 'I'},
        {"injection-file",          required_argument, 0,                 'J'},
	/* geo data flag, argument is corner freq. of high pass filter */
	{"geodata",		    required_argument, 0,		  'K'},
	/* mdc data & channel information */
        {"mdccache",                required_argument, 0,                 'L'},
        {"mdc-cache",               required_argument, 0,                 'L'},
        {"mdcchannel",              required_argument, 0,                 'M'},
        {"mdc-channel",             required_argument, 0,                 'M'},
        /* output options */
        {"printData",               no_argument,       &printData,         TRUE },
        {"printSpectrum",           no_argument,       &printSpectrum,     TRUE },
        {0, 0, 0, 0}
    };
    int c;
    size_t len=0;
    ProcessParamsTable *this_proc_param = procparams->processParamsTable;
    LALStatus           stat = blank_status;

    /*
     *
     * allocate memory 
     *
     */

    *params = LALMalloc (sizeof( EPSearchParams )); 
    if ( !*params )
    {
        fprintf(stderr, "Memory allocation failed for EPSearchParams\n");
        exit(1);
    }

    (*params)->tfTilingInput = LALMalloc (sizeof(CreateTFTilingIn));
    (*params)->initParams = LALMalloc (sizeof(EPInitParams));
    (*params)->compEPInput = LALMalloc (sizeof(ComputeExcessPowerIn));
    if ( !(*params)->tfTilingInput || !(*params)->initParams ||
         !(*params)->compEPInput )
    {
        LALFree ((*params)->tfTilingInput);
        LALFree ((*params)->initParams);
        LALFree ((*params)->compEPInput);
        LALFree (*params); *params = NULL;
        fprintf(stderr, "Memory allocation failed for tfTilingInput\n");
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
                    fprintf( stderr, "Type %s --help for options\n", *argv );
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
                    (*params)->compEPInput->alphaDefault     = alphaDflt;                    
                    ADD_PROCESS_PARAM( "float", "%e", alphaDflt );
                }
                break;

            case 'b':
                /* channel to be used in the analysis */
                (*params)->channelName = (CHAR *) LALMalloc(strlen( optarg )+1 );
                if (! (*params)->channelName ){
                    fprintf(stderr,"Error allocating memory for channel name\n");
                }
                strcpy( (*params)->channelName, optarg );
                channelIn.name = (*params)->channelName;
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
                    (*params)->events2Master = events;
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
                    (*params)->tfTilingInput->flow = flow;
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
                    (*params)->tfTilingInput->length = tmplength;
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
                    (*params)->tfTilingInput->minFreqBins = minfbin;
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
                    (*params)->tfTilingInput->minTimeBins = mintbin;
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
                    (*params)->initParams->numPoints = tmpm;
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



            case 'p':
                /* number of data points to read in at a time*/
                {
                    ramLimitPoints = atoi(optarg);
                    if (ramLimitPoints <= 0){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%i specified)\n",
                                long_options[option_index].name, ramLimitPoints);
                        exit( 1 );
                    }
                    (*params)->initParams->numSegments = ramLimitPoints;
                    ADD_PROCESS_PARAM( "int", "%d", ramLimitPoints );                    
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
                    (*params)->compEPInput->numSigmaMin  = tmpsigma;
                    ADD_PROCESS_PARAM( "float", "%e", tmpsigma );  
                }
                break;

            case 'r':
                /* Overlap betweeen segments (# of points)/no. of pts
		 *to be shifted 
		 */
                {
                    INT4 tmpolap = atoi(optarg);
                    if (tmpolap < 0){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%i specified)\n",
                                long_options[option_index].name, tmpolap);
                        exit( 1 );
                    }
                    (*params)->ovrlap = tmpolap;
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
                    (*params)->tfTilingInput->overlapFactor = tmpolap;
                    ADD_PROCESS_PARAM( "int", "%d", tmpolap ); 
                }
                break;


            case 't':
                /* Number of segments sent to slave */
                {
                    psdAveragePoints = atoi(optarg);
                    if (psdAveragePoints < 0){
                        fprintf(stderr,"invalid argument to --%s:\n"
                                "Must be > 0 (%i specified)\n",
                                long_options[option_index].name, psdAveragePoints);
                        exit( 1 );
                    }
                    (*params)->initParams->segDutyCycle = psdAveragePoints;

                    ADD_PROCESS_PARAM( "int", "%d", psdAveragePoints ); 
                }
                break;

            case 'v':
                /* Spectrum method to use */
                if ( !strcmp( optarg, "useMean" ) ) {
                    (*params)->initParams->method            = useMean;
                } 
                else if ( !strcmp( optarg, "useMedian" ) ) {
                    (*params)->initParams->method            = useMedian;
                }
                else if ( !strcmp( optarg, "useUnity" ) ) {
                    (*params)->initParams->method            = useUnity;
                }
                else {

                    fprintf(stderr,"invalid argument to --%s:\n"
                            "Must be useMean/useMedian/useUnity (%s specified)\n",
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
                frameSampleRate = (INT4) atoi( optarg );
                if ( frameSampleRate < 2 || frameSampleRate > 16384 || frameSampleRate % 2 )
                {
                  fprintf( stderr, "invalid argument to --%s:\n"
                      "rate must be power of 2 between 2 and 16384 inclusive: "
                      "(%d specified)\n", 
                      long_options[option_index].name, frameSampleRate );
                  exit( 1 );
                }
                ADD_PROCESS_PARAM( "int", "%d", frameSampleRate );
                break;

            case 'A':
                targetSampleRate = (INT4) atoi( optarg );
                if ( targetSampleRate < 2 || targetSampleRate > 16384 || targetSampleRate % 2 )
                {
                  fprintf( stderr, "invalid argument to --%s:\n"
                      "rate must be power of 2 between 2 and 16384 inclusive: "
                      "(%d specified)\n", 
                      long_options[option_index].name, targetSampleRate );
                  exit( 1 );
                }
                ADD_PROCESS_PARAM( "int", "%d", targetSampleRate );
                break;

	    case 'N':
	      if ( ! strcmp( "ldas", optarg ) )
		{
		  resampFiltType = 0;
		}
	      else if ( ! strcmp( "butterworth", optarg ) )
		{
		  resampFiltType = 1;
		}
	      else
		{
		  fprintf( stderr, "invalid argument to --%s:\n"
			   "unknown resampling filter type: "
			   "%s (must be ldas or butterworth)\n", 
			   long_options[option_index].name, optarg );
		  exit( 1 );
		}
	      ADD_PROCESS_PARAM( "string", "%s", optarg );
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
                    (*params)->alphaThreshold = tmpth;              
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
                    (*params)->winParams.type                = tmpwin;
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

            case 'I':
                /* create storage for the calibration frame cache name */
                len = strlen( optarg ) + 1;
                calCacheFile = (CHAR *) calloc( len, sizeof(CHAR));
                memcpy( calCacheFile, optarg, len );
                calFlag = TRUE;
                ADD_PROCESS_PARAM( "string", "%s", optarg );
                break;

             case 'J':
                /* create storage for the injetion file name */
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
		    fprintf(stderr, "Memory allocation failed for EPSearchParams\n");
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
                ADD_PROCESS_PARAM( "string", "%s", optarg );
                break;

           case '?':
                fprintf( stderr, "Type %s --help for options\n", *argv );
                exit( 1 );
                break;

            default:
                fprintf( stderr, "unknown error while parsing options\n" );
                exit( 1 );
        }
    }

    if ( optind < argc )
    {
        fprintf( stderr, "Error: extraneous command line arguments:\n" );
        while ( optind < argc )
        {
            fprintf ( stderr, "%s\n", argv[optind++] );
        }
        fprintf( stderr, "Type %s --help for options\n", *argv );
        exit( 1 );
    }

    /* check that start time was specified & fill the epoch structure */
    if ( !gpsStartTimeNS ){
        fprintf(stderr,"Error: --start-time must be specified\n");
        exit (1);
    }
    LAL_CALL( LALINT8toGPS( &stat, &startEpoch, &gpsStartTimeNS), &stat);

    /* check that stop time was specified & fill the epoch structure */
    if ( !gpsStopTimeNS ){
        fprintf(stderr,"Error: --stop-time must be specified\n");
        exit (1);
    }
    LAL_CALL( LALINT8toGPS( &stat, &stopEpoch, &gpsStopTimeNS), &stat);

    /* check that the sample rate has been supplied */
    if ( frameSampleRate < 0 )
    {
      fprintf( stderr, "--frame-sample-rate must be specified. This is the expected sampling rate of data in frames\n" );
      exit( 1 );
    }

   /* check that the target sample rate has been supplied */
    if ( targetSampleRate < 0 )
    {
      fprintf( stderr, "--target-sample-rate must be specified. This will be the sampling rate of the resampled time series \n" );
      exit( 1 );
    }

    /* compute the total number of points to be analyzed */
    totalNumPoints = ( stopEpoch.gpsSeconds  - startEpoch.gpsSeconds ) 
      * frameSampleRate;

    /* initialize parameter structures */
    (*params)->numSlaves    = NULL;

    /* initialize parameters */
    (*params)->haveData        = 0;
    (*params)->numEvents       = 0;
    (*params)->searchMaster    = 0;
    (*params)->tfTilingInput->maxTileBand = 64.0;

    return 0;
}

#undef ADD_PROCESS_PARAMS

#endif
