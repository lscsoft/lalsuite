#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lal/AVFactories.h>
#include <lal/BurstSearch.h>
#include <lal/Date.h>
#include <lal/EPSearch.h>
#include <lal/ExcessPower.h>
#include <lal/FrameCalibration.h>
#include <lal/FrameStream.h>
#include <lal/GenerateBurst.h>
#include <lal/IIRFilter.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALError.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/PrintFTSeries.h>
#include <lal/Random.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#include <lalapps.h>
#include <processtable.h>

/* declare the parsing function which is at the end of the file */
int snprintf(char *str, size_t size, const  char  *format, ...);
void initializeEPSearch(int argc, char *argv[], EPSearchParams **params, MetadataTable *procparams);

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

/* this is a temporary function to make some code below more readable.
 * Obviously the argument and return types are sub-optimal... */
static size_t min(size_t a, size_t b)
{
	return(a < b ? a : b);
}

/* Parameters from command line */
static struct {
	LIGOTimeGPS startEpoch;     /* gps start time                      */
	LIGOTimeGPS stopEpoch;      /* gps stop time                       */
	INT4 maxSeriesLength;       /* RAM-limited input length            */
	size_t psdAverageLength;    /* number of samples to use for PSD    */
} options;

/* some global output flags */
INT4 verbose;
INT4 geodata;
INT4 printData;
INT4 whiteNoise;                    /* insertion of Gaussian white noise   */
INT4 sineBurst;                     /* insertion of shaped sine burst      */
 
/* global variables */
FrChanIn channelIn;                 /* channnel information                */
FrChanIn mdcchannelIn;              /* mdc signal only channnel info       */
EPSearchParams *mdcparams;          /* mdc search param                    */
CHAR ifo[3];                        /* two character interferometer        */
CHAR comment[LIGOMETA_COMMENT_MAX]; /* string in output file name          */
CHAR *cachefile;                    /* name of file with frame cache info  */
CHAR *dirname;                      /* name of directory with frames       */
REAL4 noiseAmpl;                    /* gain factor for white noise         */
INT4 seed;                          /* set non-zero to generate noise      */
INT4 totalNumPoints;                /* total number of points to analyze   */
UINT4 totalNumSegs;                 /* total number of segments to analyze */
INT4 frameSampleRate;               /* sample rate of the frame data       */
INT4 targetSampleRate;              /* sample rate after resampling        */

/* data conditioning parameters */
CHAR *calCacheFile;                 /* name of the calibration cache file  */
CHAR *injectionFile;                /* file with list of injections        */
CHAR *mdcCacheFile;                 /* name of mdc signal cache file       */
ResampleTSFilter resampFiltType;

/* GEO data high pass corner freq. */
REAL8 fcorner;                      /* corner frequency in Hz              */


/* Fill an array with white noise */
static void makeWhiteNoise(
        LALStatus             *status,
        REAL4TimeSeries       *series
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


/*
 * round a psd length down to a value that is commensurate with the window
 * length and window shift.
 */

static size_t window_commensurate(size_t psdlength, size_t windowlength, size_t windowshift)
{
	if(psdlength < windowlength)
		return(0);
	return(((psdlength - windowlength) / windowshift) * windowshift + windowlength);
}


/*
 * Entry point
 */

int main( int argc, char *argv[])
{
    LALStatus               stat;
    LALLeapSecAccuracy      accuracy = LALLEAPSEC_LOOSE;

    EPSearchParams         *params = NULL;
    FrStream               *stream = NULL;
    FrCache                *frameCache = NULL;
    CHAR                    fname[256];
    PassBandParamStruc      highpassParam;
    REAL4                   fsafety = 0;
    CalibrationUpdateParams calfacts;
    LIGOTimeGPS		    tmpEpoch = {0, 0};
    LALTimeInterval         tmpInterval;
    REAL8                   tmpOffset = 0.0;
    REAL4                   minFreq = 0.0;    
    size_t                  start_sample;
    int                     usedNumPoints;
    INT4                    numPoints;  /* number of samples from frames */

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

    memset(&stat, 0, sizeof(stat));

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
    tmpEpoch = options.startEpoch;

    /******************************************************************
     * OUTER LOOP over data small enough to fit into memory 
     ******************************************************************/
    for(usedNumPoints = 0; (totalNumPoints - usedNumPoints) > (int) (2 * params->windowShift * (frameSampleRate / targetSampleRate)); usedNumPoints += numPoints - 2 * params->windowShift * (frameSampleRate / targetSampleRate)) {

      /* tell operator how we are doing */
      if (verbose){
        fprintf(stdout,"%i points analysed && %i points left\n", usedNumPoints, totalNumPoints - usedNumPoints);
      }

      /* compute the number of points to use in this run */
      numPoints = min(options.maxSeriesLength, (totalNumPoints - usedNumPoints));

      if(verbose) {
        fprintf(stdout,"reading %i points...\n", numPoints);
      }

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
        makeWhiteNoise(&stat, &series);
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
      if( injectionFile )
      {
        INT4  startTime = series.epoch.gpsSeconds;
        INT4  stopTime = startTime + (INT4)( series.data->length * series.deltaT );
        SimBurstTable *injections = NULL;

        if ( !calCacheFile )
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

      if(mdcCacheFile) {
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
      LAL_CALL( EPConditionData( &stat, &series, minFreq, 1.0/targetSampleRate, resampFiltType, params->windowShift, params), &stat);

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

		/***********************************************
		 * DO THE SEARCH                               *
		 ***********************************************/

		if (verbose)
			fprintf(stdout,"Got %i points to analyse after conditioning\n", series.data->length);

		if(options.psdAverageLength > series.data->length)
			options.psdAverageLength = window_commensurate(series.data->length, params->windowLength, params->windowShift);

		for(start_sample = 0; start_sample + (params->windowLength - params->windowShift) < series.data->length; start_sample += options.psdAverageLength - (params->windowLength - params->windowShift)) {
			REAL4TimeSeries *interval;

			if(start_sample + options.psdAverageLength > series.data->length)
				start_sample = series.data->length - options.psdAverageLength;

			LAL_CALL(LALCutREAL4TimeSeries(&stat, &interval, &series, start_sample, options.psdAverageLength), &stat);

			if(verbose)
				fprintf(stdout,"Analyzing samples %i -- %i\n", start_sample, start_sample + interval->data->length);

			LAL_CALL(EPSearch(&stat, interval, params, EventAddPoint), &stat);
			while(*EventAddPoint)
				EventAddPoint = &(*EventAddPoint)->next;

			LAL_CALL(LALDestroyREAL4TimeSeries(&stat, interval), &stat);
		}
 

		/* compute the start time for the next chunk */
		tmpOffset = (numPoints - 2.0 * params->windowShift * frameSampleRate / targetSampleRate) / (REAL8) frameSampleRate;
		LAL_CALL(LALFloatToInterval(&stat, &tmpInterval, &tmpOffset), &stat);
		LAL_CALL(LALIncrementGPS(&stat, &tmpEpoch, &tmpEpoch, &tmpInterval), &stat);

		/* clean up memory from that run */
		LAL_CALL(LALSDestroyVector(&stat, &series.data), &stat);
		if(calCacheFile)
			LAL_CALL(LALCDestroyVector(&stat, &resp.data), &stat);
	}
	/*******************************************************************
	 * outer-most loop ends here
	 *******************************************************************/

	LALFree(cachefile); 
	LALFree(calCacheFile);
	LALFree(dirname); 
	if(mdcparams) {
		LALFree(mdcparams->channelName);
		LALFree(mdcparams);
	}

    /*******************************************************************
    * OUTPUT THE RESULTS 
    *******************************************************************/
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
    snprintf( fname, sizeof(fname), "%s-%s-POWER-%d-%d.xml",
            ifo, comment, options.startEpoch.gpsSeconds, 
            options.stopEpoch.gpsSeconds-options.startEpoch.gpsSeconds);
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


	/*****************************************************************
	 * FINAL CLEANUP                                                 *
	 *****************************************************************/
	if(params) {
		LALFree(params->channelName);
		LALFree(params->compEPInput);
		LALFree(params->tfTilingInput);
		LALFree(params);
	}
	LALCheckMemoryLeaks();
	exit(0);
}


/*************************************************************************
 * INITIALIZE AND PARSE ARGUMENTS
 *************************************************************************/

static void print_usage(char *program)
{
	fprintf(stderr,
"Usage:  %s <options...>\n" \
"	[--cluster]\n" \
"	 --bandwidth <bandwidth>\n" \
"	[--calibration-cache <cache file>]\n" \
"	 --channel-name <string>\n" \
"	[--debug-level <level>]\n" \
"	 --default-alpha <alpha>\n" \
"	 --event-limit <count>\n" \
"	 --frame-cache <cache file>\n" \
"	 --frame-dir <directory>\n" \
"	 --frame-sample-rate <Hz>\n" \
"	[--geodata <high pass corner frequency>]\n" \
"	 --gps-end-time <seconds>\n" \
"	 --gps-end-time-ns <nanoseconds>\n" \
"	 --gps-start-time <seconds>\n" \
"	 --gps-start-time-ns <nanoseconds>\n" \
"	[--help]\n" \
"	[--injection-file <file name>]\n" \
"	 --low-freq-cutoff <flow>\n" \
"	[--mdc-cache <cache file>]\n" \
"	[--mdc-channel <channel name>]\n" \
"	 --min-freq-bin <nfbin>\n" \
"	 --min-time-bin <ntbin>\n" \
"	 --noise-amplitude <amplitude>\n" \
"	 --nsigma <sigma>\n" \
"	[--printData]\n" \
"	[--printSpectrum]\n" \
"	 --psd-average-method <method>\n" \
"	 --psd-average-points <samples>\n" \
"	 --ram-limit <MebiBytes>\n" \
"	 --resample-filter <filter type>\n" \
"	 --seed <seed>\n" \
"	 --target-sample-rate <Hz>\n" \
"	 --tile-overlap-factor <factor>\n" \
"	 --threshold <threshold>\n" \
"	[--user-tag <comment>]\n" \
"	[--verbose]\n" \
"	 --window <window>\n" \
"	 --window-length <samples>\n" \
"	 --window-shift <samples>\n", program);
}

static void print_bad_argument(const char *prog, const char *arg, const char *msg)
{
	fprintf(stderr, "%s: error: invalid argument for --%s: %s\n", prog, arg, msg);
}

static void print_missing_argument(const char *prog, const char *arg)
{
	fprintf(stderr, "%s: error: --%s not specified\n", prog, arg);
}

static void print_alloc_fail(const char *prog, const char *msg)
{
	fprintf(stderr, "%s: error: memory allocation failure %s\n", prog, msg);
}


#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) { \
	proc_param->next = LALCalloc(1, sizeof(ProcessParamsTable)); \
	proc_param = proc_param->next; \
	snprintf( proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue ); \
	snprintf( proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME ); \
	snprintf( proc_param->param, LIGOMETA_PARAM_MAX, "--%s", long_options[option_index].name ); \
	snprintf( proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
}

static int check_for_missing_parameters(LALStatus *stat, char *prog, struct option *long_options, EPSearchParams *params)
{
	int index;
	int got_all_arguments = TRUE;
	int arg_is_missing;
	INT8 gpstmp;

	for(index = 0; long_options[index].name; index++) {
		switch(long_options[index].val) {
			case 'I':
			arg_is_missing = !frameSampleRate;
			break;

			case 'K':
			LAL_CALL(LALGPStoINT8(stat, &gpstmp, &options.stopEpoch), stat);
			arg_is_missing = !gpstmp;
			break;

			case 'M':
			LAL_CALL(LALGPStoINT8(stat, &gpstmp, &options.startEpoch), stat);
			arg_is_missing = !gpstmp;
			break;

			case 'Z':
			arg_is_missing = !options.psdAverageLength;
			break;

			case 'a':
			arg_is_missing = !options.maxSeriesLength;
			break;

			case 'd':
			arg_is_missing = !params->windowShift;
			break;

			case 'e':
			arg_is_missing = !targetSampleRate;
			break;

			default:
			arg_is_missing = FALSE;
			break;
		}
		if(arg_is_missing) {
			print_missing_argument(prog, long_options[index].name);
			got_all_arguments = FALSE;
		}
	}
	return(got_all_arguments);
}


static INT8 DeltaGPStoINT8(LALStatus *stat, LIGOTimeGPS *stop, LIGOTimeGPS *start)
{
	LALTimeInterval deltainterval;
	LIGOTimeGPS deltagps;
	INT8 deltaNS;

	LAL_CALL(LALDeltaGPS(stat, &deltainterval, stop, start), stat);
	deltagps.gpsSeconds = deltainterval.seconds;
	deltagps.gpsNanoSeconds = deltainterval.nanoSeconds;
	LAL_CALL(LALGPStoINT8(stat, &deltaNS, &deltagps), stat);

	return(deltaNS);
}


void initializeEPSearch( 
	int argc, 
	char *argv[], 
	EPSearchParams **params,
	MetadataTable *procparams 
)
{ 
	char msg[240];
	int cluster;
	int printSpectrum;
	int c;
	int option_index;
	int ram = 0;
	INT8 gpstmp;
	INT8 gpsStartTimeNS = 0;
	INT8 gpsStopTimeNS = 0;
	ProcessParamsTable *proc_param = procparams->processParamsTable;
	LALStatus stat = blank_status;
	struct option long_options[] = {
		{"cluster",             no_argument,       &cluster,       TRUE},
		{"bandwidth",           required_argument, NULL,           'A'},
		{"calibration-cache",   required_argument, NULL,           'B'},
		{"channel-name",        required_argument, NULL,           'C'},
		{"debug-level",         required_argument, NULL,           'D'},
		{"default-alpha",       required_argument, NULL,           'E'},
		{"event-limit",         required_argument, NULL,           'F'},
		{"frame-cache",         required_argument, NULL,           'G'},
		{"frame-dir",           required_argument, NULL,           'H'},
		{"frame-sample-rate",   required_argument, NULL,           'I'},
		{"geodata",		required_argument, NULL,           'J'},
		{"gps-end-time",        required_argument, NULL,           'K'},
		{"gps-end-time-ns",     required_argument, NULL,           'L'},
		{"gps-start-time",      required_argument, NULL,           'M'},
		{"gps-start-time-ns",   required_argument, NULL,           'N'},
		{"help",                no_argument,       NULL,           'O'},
		{"injection-file",      required_argument, NULL,           'P'},
		{"low-freq-cutoff",     required_argument, NULL,           'Q'},
		{"mdc-cache",           required_argument, NULL,           'R'},
		{"mdc-channel",         required_argument, NULL,           'S'},
		{"min-freq-bin",        required_argument, NULL,           'T'},
		{"min-time-bin",        required_argument, NULL,           'U'},
		{"noise-amplitude",     required_argument, NULL,           'V'},
		{"nsigma",              required_argument, NULL,           'X'},
		{"printData",           no_argument,       &printData,     TRUE},
		{"printSpectrum",       no_argument,       &printSpectrum, TRUE},
		{"psd-average-method",  required_argument, NULL,           'Y'},
		{"psd-average-points",  required_argument, NULL,           'Z'},
		{"ram-limit",           required_argument, NULL,           'a'},
		{"resample-filter",     required_argument, NULL,           'b'},
		{"seed",                required_argument, NULL,           'c'},
		{"target-sample-rate",  required_argument, NULL,           'e'},
		{"tile-overlap-factor", required_argument, NULL,           'f'},
		{"threshold",           required_argument, NULL,           'g'},
		{"user-tag",            required_argument, NULL,           'h'},
		{"verbose",             no_argument,       &verbose,       TRUE},
		{"window",              required_argument, NULL,           'i'},
		{"window-length",       required_argument, NULL,           'W'},
		{"window-shift",        required_argument, NULL,           'd'},
		{NULL, 0, NULL, 0}
	};

	/*
	 * Allocate memory.
	 */

	*params = LALMalloc(sizeof(EPSearchParams));
	if(!*params) {
		print_alloc_fail(argv[0], "for *params");
		exit(1);
	}

	(*params)->tfTilingInput = LALMalloc(sizeof(*(*params)->tfTilingInput));
	(*params)->compEPInput = LALMalloc(sizeof(*(*params)->compEPInput));
	if(!(*params)->tfTilingInput || !(*params)->compEPInput) {
		LALFree((*params)->tfTilingInput);
		LALFree((*params)->compEPInput);
		LALFree(*params);
		*params = NULL;
		print_alloc_fail(argv[0], "");
		exit(1);
	}

	/*
	 * Set parameter defaults.
	 */

	(*params)->channelName = NULL;
	(*params)->eventLimit = 999;
	(*params)->tfTilingInput->maxTileBand = 64.0;
	(*params)->windowShift = 0;

	cachefile = NULL;
	calCacheFile = NULL;
	cluster = FALSE;
	dirname = NULL;
	fcorner = 100.0;
	frameSampleRate = 0;
	geodata = FALSE;
	memset(ifo, 0, sizeof(ifo));
	injectionFile = NULL;
	mdcCacheFile = NULL;
	mdcparams = NULL;
	noiseAmpl = 1.0;
	printData = FALSE;
	printSpectrum = FALSE;
	options.psdAverageLength = 0;
	resampFiltType = -1;
	seed = 1;
	sineBurst = FALSE;
	targetSampleRate = 0;
	totalNumPoints = 0;
	totalNumSegs = 0;
	verbose = FALSE;
	whiteNoise = FALSE;

	/*
	 * Parse command line.
	 */

	do switch(c = getopt_long(argc, argv, "", long_options, &option_index)) {
		case 'A':
		(*params)->tfTilingInput->length = atoi(optarg);
		if((*params)->tfTilingInput->length <= 0) {
			sprintf(msg, "must be > 0 (%i specified)", (*params)->tfTilingInput->length);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		ADD_PROCESS_PARAM("int", "%d", (*params)->tfTilingInput->length);
		break;

		case 'B':
		calCacheFile = LALCalloc(strlen(optarg) + 1, sizeof(*calCacheFile));
		if(calCacheFile) {
			strcpy(calCacheFile, optarg);
			ADD_PROCESS_PARAM("string", "%s", optarg);
		}
		break;

		case 'C':
		(*params)->channelName = (CHAR *) LALMalloc(strlen(optarg) + 1);
		if(!(*params)->channelName) {
			print_alloc_fail(argv[0], "for (*params)->channelName");
			exit(1);
		}
		strcpy((*params)->channelName, optarg);
		channelIn.name = (*params)->channelName;
		channelIn.type = ADCDataChannel;
		memcpy(ifo, channelIn.name, sizeof(ifo) - 1);
		ADD_PROCESS_PARAM("string", "%s", optarg);
		break;

		case 'D':
		set_debug_level(optarg);
		ADD_PROCESS_PARAM("int", "%d", atoi(optarg));
		break;

		case 'E':
		(*params)->compEPInput->alphaDefault = atof(optarg);
		if((*params)->compEPInput->alphaDefault <= 0.0 || (*params)->compEPInput->alphaDefault >= 1.0) {
			sprintf(msg, "must be in range [0,1] (%f specified)", (*params)->compEPInput->alphaDefault);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		ADD_PROCESS_PARAM("float", "%e", (*params)->compEPInput->alphaDefault);
		break;

		case 'F':
		(*params)->eventLimit = atoi(optarg);
		if((*params)->eventLimit < 1 || (*params)->eventLimit > 999) {
			sprintf(msg, "must be in range [1,999] (%i specified)", (*params)->eventLimit);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit( 1 );
		}
		ADD_PROCESS_PARAM("int", "%d", (*params)->eventLimit);
		break;

		case 'G':
		cachefile = LALCalloc(strlen(optarg) + 1, sizeof(*cachefile));
		if(cachefile) {
			strcpy(cachefile, optarg);
			ADD_PROCESS_PARAM("string", "%s", optarg);
		}
		break;

		case 'H':
		dirname =  LALCalloc(strlen(optarg) + 1, sizeof(*dirname));
		if(dirname) {
			strcpy(dirname, optarg);
			ADD_PROCESS_PARAM("string", "%s", optarg);
		}
		break;

		case 'I':
		frameSampleRate = (INT4) atoi(optarg);
		if(frameSampleRate < 2 || frameSampleRate > 16384 || frameSampleRate % 2) {
			sprintf(msg, "must be a power of 2 in the range [2,16384] (%d specified)", frameSampleRate);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		ADD_PROCESS_PARAM("int", "%d", frameSampleRate);
		break;

		case 'J':
		fcorner = atof(optarg);
		if(fcorner <= 0) {
			print_bad_argument(argv[0], long_options[option_index].name, "must be > 0");
			exit(1);
		}
		geodata = TRUE;
		break;

		case 'K':
		gpstmp = atol(optarg);
		if(gpstmp < 441417609 || gpstmp > 999999999) {
			sprintf(msg, "must be in the range [Jan 01 1994 00:00:00 UTC, Sep 14 2011 01:46:26 UTC] (%lld specified)", gpstmp);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		gpsStopTimeNS += gpstmp * 1000000000LL;
		ADD_PROCESS_PARAM("int", "%lld", gpstmp);
		break;

		case 'L':
		gpstmp = atol(optarg);
		if(gpstmp < 0 || gpstmp > 999999999) {
			sprintf(msg, "must be in the range [0,999999999] (%lld specified)", gpstmp);
			print_bad_argument(argv[0], long_options[option_index].name, "GPS start time nanoseconds is negative");
			exit(1);
		}
		gpsStopTimeNS += gpstmp;
		ADD_PROCESS_PARAM("int", "%lld", gpstmp);
		break;

		case 'M':
		gpstmp = atol(optarg);
		if(gpstmp < 441417609 || gpstmp > 999999999) {
			sprintf(msg, "must be in the range [Jan 01 1994 00:00:00 UTC, Sep 14 2011 01:46:26 UTC] (%lld specified)", gpstmp);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		gpsStartTimeNS += gpstmp * 1000000000LL;
		ADD_PROCESS_PARAM("int", "%lld", gpstmp);
		break;

		case 'N':
		gpstmp = atol(optarg);
		if(gpstmp < 0 || gpstmp > 999999999) {
			sprintf(msg, "must be in the range [0,999999999] (%lld specified)", gpstmp);
			print_bad_argument(argv[0], long_options[option_index].name, "GPS start time nanoseconds is negative");
			exit(1);
		}
		gpsStartTimeNS += gpstmp;
		ADD_PROCESS_PARAM("int", "%lld", gpstmp);
		break;

		case 'O':
		print_usage(argv[0]);
		exit(0);
		break;

		case 'P':
		injectionFile = LALCalloc(strlen(optarg) + 1, sizeof(*injectionFile));
		if(injectionFile) {
			strcpy(injectionFile, optarg);
			ADD_PROCESS_PARAM("string", "%s", optarg);
		}
		break;

		case 'Q':
		(*params)->tfTilingInput->flow = atof(optarg);
		if((*params)->tfTilingInput->flow < 0.0) {
			sprintf(msg,"must be >= 0.0 (%f specified)", (*params)->tfTilingInput->flow);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		ADD_PROCESS_PARAM("float", "%e", (*params)->tfTilingInput->flow);
		break;

		case 'R':
		mdcCacheFile = LALCalloc(strlen(optarg) + 1, sizeof(*mdcCacheFile));
		if(mdcCacheFile) {
			strcpy(mdcCacheFile, optarg);
			ADD_PROCESS_PARAM("string", "%s", optarg);
		}
		break;

		case 'S':
		mdcparams = LALMalloc(sizeof(*mdcparams));
		if(!mdcparams) {
			print_alloc_fail(argv[0], "for mdcparams");
			exit(1);
		}
		mdcparams->channelName = LALMalloc(strlen(optarg) + 1);
		if(!mdcparams->channelName) {
			print_alloc_fail(argv[0], "for mdcparams->channelName");
			exit(1);
		}
		strcpy(mdcparams->channelName, optarg);
		mdcchannelIn.name = mdcparams->channelName;
		mdcchannelIn.type = ADCDataChannel;
		ADD_PROCESS_PARAM("string", "%s", optarg);
		break;

		case 'T':
		(*params)->tfTilingInput->minFreqBins = atoi(optarg);
		if((*params)->tfTilingInput->minFreqBins <= 0) {
			sprintf(msg,"must be > 0 (%i specified)", (*params)->tfTilingInput->minFreqBins);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		ADD_PROCESS_PARAM("int", "%d", (*params)->tfTilingInput->minFreqBins);
		break;

		case 'U':
		(*params)->tfTilingInput->minTimeBins = atoi(optarg);
		if((*params)->tfTilingInput->minTimeBins <= 0) {
			sprintf(msg,"must be > 0 (%i specified)", (*params)->tfTilingInput->minTimeBins);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		ADD_PROCESS_PARAM("int", "%d", (*params)->tfTilingInput->minTimeBins);
		break;

		case 'V':
		noiseAmpl = atof(optarg);
		whiteNoise = TRUE;
		if(noiseAmpl <= 0.0) {
			sprintf(msg, "must be > 0 (%f specified)", noiseAmpl);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		ADD_PROCESS_PARAM("float", "%e", noiseAmpl);
		break;

		case 'W':
		(*params)->windowLength = atoi(optarg);
		if((*params)->windowLength <= 0) {
			sprintf(msg, "must be > 0 (%i specified)", (*params)->windowLength);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		ADD_PROCESS_PARAM("int", "%d", (*params)->windowLength);
		break;

		case 'X':
		(*params)->compEPInput->numSigmaMin = atof(optarg);
		if((*params)->compEPInput->numSigmaMin <= 1.0) {
			sprintf(msg, "must be > 0 (%f specified)", (*params)->compEPInput->numSigmaMin);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		ADD_PROCESS_PARAM("float", "%e", (*params)->compEPInput->numSigmaMin);
		break;

		case 'Y':
		if(!strcmp(optarg, "useMean"))
			(*params)->method = useMean;
		else if(!strcmp(optarg, "useMedian"))
			(*params)->method = useMedian;
		else if (!strcmp(optarg, "useUnity"))
			(*params)->method = useUnity;
		else {
			sprintf(msg, "must be \"useMean\", \"useMedian\", or \"useUnity\"");
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		ADD_PROCESS_PARAM("string", "%s", optarg);
		break;

		case 'Z':
		options.psdAverageLength = atoi(optarg);
		ADD_PROCESS_PARAM("int", "%d", options.psdAverageLength);
		break;

		case 'a':
		ram = atoi(optarg);
		if(ram <= 0) {
			sprintf(msg, "must be > 0 (%i specified)", ram);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		ADD_PROCESS_PARAM("int", "%d", ram);
		break;

		case 'b':
		if(!strcmp("ldas", optarg))
			resampFiltType = LDASfirLP;
		else if(!strcmp("butterworth", optarg))
			resampFiltType = defaultButterworth;
		else {
			sprintf(msg, "must be \"ldas\", or \"butterworth\"");
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		ADD_PROCESS_PARAM("string", "%s", optarg);
		break;

		case 'c':
		seed = atoi(optarg);
		if(seed <= 0) {
			sprintf(msg, "must be > 0 (%i specified)", seed);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		ADD_PROCESS_PARAM("int", "%d", seed);
		break;

		case 'd':
		(*params)->windowShift = atoi(optarg);
		ADD_PROCESS_PARAM("int", "%d", (*params)->windowShift);
		break;

		case 'e':
		targetSampleRate = (INT4) atoi(optarg);
		if(targetSampleRate < 2 || targetSampleRate > 16384 || targetSampleRate % 2) {
			sprintf(msg, "must be a power of 2 in the rage [2,16384] (%d specified)", targetSampleRate);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		ADD_PROCESS_PARAM("int", "%d", targetSampleRate);
		break;

		case 'f':
		(*params)->tfTilingInput->overlapFactor = atoi(optarg);
		if((*params)->tfTilingInput->overlapFactor < 0) {
			sprintf(msg, "must be > 0 (%i specified)", (*params)->tfTilingInput->overlapFactor);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		ADD_PROCESS_PARAM("int", "%d", (*params)->tfTilingInput->overlapFactor);
		break;

		case 'g':
		(*params)->alphaThreshold = atof(optarg);
		if((*params)->alphaThreshold < 0.0) {
			sprintf(msg, "must be > 0 (%f specified)", (*params)->alphaThreshold);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		ADD_PROCESS_PARAM("float", "%e", (*params)->alphaThreshold);
		break;

		case 'h':
		snprintf(comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
		ADD_PROCESS_PARAM("string", "%s", optarg);
		break;

		case 'i':
		(*params)->windowType = atoi(optarg);
		if((*params)->windowType >= NumberWindowTypes) {
			sprintf(msg, "must be <= %d (%i specified)", NumberWindowTypes, (*params)->windowType);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			exit(1);
		}
		ADD_PROCESS_PARAM("int", "%d", (*params)->windowType);
		break;

		/* option sets a flag */
		case 0:
		break;

		/* end of arguments */
		case -1:
		break;

		/* unrecognized option */
		case '?':
		print_usage(argv[0]);
		exit(1);
		break;

		/* missing argument for an option */
		case ':':
		print_usage(argv[0]);
		exit(1);
		break;
	} while(c != -1);

	/*
	 * Convert the start and stop epochs, along with the sample rate, to a
	 * total number of samples to analyze.
	 */

	LAL_CALL(LALINT8toGPS(&stat, &options.startEpoch, &gpsStartTimeNS), &stat);
	LAL_CALL(LALINT8toGPS(&stat, &options.stopEpoch, &gpsStopTimeNS), &stat);
	totalNumPoints = DeltaGPStoINT8(&stat, &options.stopEpoch, &options.startEpoch) * frameSampleRate / 1000000000L;

	/*
	 * Convert the amount of available RAM to a limit on the length of a
	 * time series to read in.
	 */

	options.maxSeriesLength = (ram * 1024 * 1024) / (4 * sizeof(REAL4));

	/*
	 * Check for missing parameters
	 */

	if(!check_for_missing_parameters(&stat, argv[0], long_options, *params))
		exit(1);

	/*
	 * Ensure psdAverageLength is comensurate with the analysis window
	 * length and shift.
	 */

	options.psdAverageLength = window_commensurate(options.psdAverageLength, (*params)->windowLength, (*params)->windowShift);

	/*
	 * Miscellaneous chores.
	 */

	(*params)->cluster = cluster;
	(*params)->printSpectrum = printSpectrum;

	if(verbose) {
		fprintf(stderr, "%s: available RAM limits analysis to %d samples at a time\n", argv[0], options.maxSeriesLength);
		fprintf(stderr, "%s: using --psd-average-points %d\n", argv[0], options.psdAverageLength);
	}
}

#undef ADD_PROCESS_PARAMS

#endif
