/*----------------------------------------------------------------------- 
 * 
 * File Name: tmpltbank.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <time.h>

#include <FrameL.h>
#include <lalapps.h>
#include <series.h>
#include <processtable.h>
#include <lalappsfrutils.h>

#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/PrintFTSeries.h>
#include <lal/FrameStream.h>
#include <lal/FrameCalibration.h>
#include <lal/Window.h>
#include <lal/TimeFreqFFT.h>
#include <lal/IIRFilter.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>
#include "inspiral.h"

RCSID( "$Id$" );
#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "tmpltbank"

int arg_parse_check( int argc, char *argv[], MetadataTable procparams );

/* type of data to analyze */
enum
{
  undefined,
  real_4,
  real_8
} calData = undefined;

/*
 *
 * variables that control program behaviour
 *
 */


/* debugging */
extern int vrbflg;                      /* verbocity of lal function    */

/* parameters used to generate calibrated power spectrum */
LIGOTimeGPS gpsStartTime = { 0, 0 };    /* input data GPS start time    */
LIGOTimeGPS gpsEndTime = { 0, 0 };      /* input data GPS end time      */
INT4  padData = 0;                      /* saftety margin on input data */
CHAR  *fqChanName       = NULL;         /* name of data channel         */
INT4  globFrameData     = 0;            /* glob *.gwf to get frame data */
CHAR  *frInCacheName    = NULL;         /* cache file containing frames */
CHAR  *frInType         = NULL;         /* type of data frames          */
INT4  numPoints         = -1;           /* points in a segment          */
INT4  numSegments       = -1;           /* number of segments           */
CHAR  ifo[3];                           /* two character ifo code       */
CHAR *channelName = NULL;               /* channel string               */
INT4  inputDataLength = 0;              /* number of points in input    */
INT4   resampFiltType   = -1;           /* low pass filter used for res */
INT4   sampleRate       = -1;           /* sample rate of filter data   */
INT4   highPass         = -1;           /* enable high pass on raw data */
REAL4  highPassFreq     = 0;            /* high pass frequency          */
INT4   highPassOrder    = -1;           /* order of the td iir filter   */
REAL4  highPassAtten    = -1;           /* attenuation of the td filter */
REAL4  fLow             = -1;           /* low frequency cutoff         */
INT4   specType         = -1;           /* use median or mean psd       */
CHAR  *calCacheName     = NULL;         /* location of calibration data */
INT4   globCalData      = 0;            /* glob for calibration frames  */
INT4   pointCal         = 0;            /* don't average cal over chunk */
REAL4  dynRangeExponent = 0;            /* exponent of dynamic range    */
REAL4 geoHighPassFreq = -1;             /* GEO high pass frequency      */
INT4  geoHighPassOrder = -1;            /* GEO high pass filter order   */
REAL4 geoHighPassAtten = -1;            /* GEO high pass attenuation    */

/* template bank generation parameters */
REAL4   minMass         = -1;           /* minimum component mass       */
REAL4   maxMass         = -1;           /* maximum component mass       */
REAL4   psi0Min         = 0;            /* minimum value of psi0        */
REAL4   psi0Max         = 0;            /* maximum value of psi0        */
REAL4   psi3Min         = 0;            /* minimum value of psi3        */
REAL4   psi3Max         = 0;            /* maximum value of psi3        */
REAL4   alpha           = 0;            /* BCV amplitude correction     */
INT4    maxFcutTmplts   = -1;           /* num tmplts in fcut direction */
REAL4   minMatch        = -1;           /* minimum requested match      */
REAL4   fUpper          = -1;           /* upper frequency cutoff       */
Order   order;                          /* post-Newtonian order         */
Approximant approximant;                /* approximation method         */
CoordinateSpace space;                  /* coordinate space used        */
GridSpacing gridSpacing = SquareNotOriented; /* grid spacing (square or hexa)*/
INT4   	isMaxTotMass    = 0;            /* Use a maximum total mass?	*/

/* standard candle parameters */
INT4    computeCandle = 0;              /* should we compute a candle?  */
REAL4   candleSnr = -1;                 /* candle signal to noise ratio */
REAL4   candleMass1 = -1;               /* standard candle mass (solar) */
REAL4   candleMass2 = -1;               /* standard candle mass (solar) */

/* output parameters */
CHAR  *userTag          = NULL;
int    writeRawData     = 0;            /* write the raw data to a file */
int    writeResponse    = 0;            /* write response function used */
int    writeSpectrum    = 0;            /* write computed psd to file   */
int    writeStrainSpec  = 0;            /* write computed stain spec    */

/* other command line args */
CHAR comment[LIGOMETA_COMMENT_MAX];     /* process param comment        */

int main ( int argc, char *argv[] )
{
  /* lal function variables */
  LALStatus             status = blank_status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;

  /* frame input data */
  FrCache      *frInCache = NULL;
  FrCache      *frGlobCache = NULL;
  FrCache      *calCache = NULL;
  FrStream     *frStream = NULL;
  FrChanIn      frChan;
  FrCacheSieve  sieve;
  const size_t  calGlobLen = 12;
  CHAR         *calGlobPattern;

  /* frame output data */
  struct FrFile *frOutFile = NULL;
  struct FrameH *outFrame  = NULL;

  /* raw input data storage */
  REAL4TimeSeries               chan;
  REAL8TimeSeries               geoChan;
  REAL4FrequencySeries          spec;
  COMPLEX8FrequencySeries       resp;

  /* structures for preconditioning */
  ResampleTSParams              resampleParams;
  LALWindowParams               wpars;
  AverageSpectrumParams         avgSpecParams;

  /* templates */
  InspiralCoarseBankIn          bankIn;
  SnglInspiralTable            *tmplt  = NULL;
  INT4                          numCoarse = 0;

  /* output data */
  MetadataTable         templateBank;
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         searchsumm;
  MetadataTable         searchsummvars;
  MetadataTable         candle;
  SummValueTable      **this_summvalue = &(candle.summValueTable);
  SearchSummvarsTable  *this_search_summvar;
  ProcessParamsTable   *this_proc_param;
  LIGOLwXMLStream       results;

  /* counters and other variables */
  UINT4 cut, i, j, k;
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
  CHAR  fname[256];
  LALUnitPair pair;
  REAL8 respRe, respIm;
  REAL8 shf;
  REAL8 inputLengthNS;
  UINT4 numInputPoints;
  const REAL8 epsilon = 1.0e-8;
  UINT4 resampleChan = 0;
  REAL8 tsLength;
  CalibrationUpdateParams calfacts;
  REAL8 dynRange = 0;


  /*
   *
   * initialization
   *
   */


  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) 
    calloc( 1, sizeof(ProcessTable) );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->start_time),
        &accuracy ), &status );
  LAL_CALL( populate_process_table( &status, proctable.processTable, 
        PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &status );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );

  /* create the search summary and zero out the summvars table */
  searchsumm.searchSummaryTable = (SearchSummaryTable *)
    calloc( 1, sizeof(SearchSummaryTable) );
  searchsummvars.searchSummvarsTable = NULL;

  /* call the argument parse and check function */
  arg_parse_check( argc, argv, procparams );

  /* can use LALMalloc() / LALCalloc() from here */

  /* fill the comment, if a user has specified on, or leave it blank */
  if ( ! *comment )
  {
    LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
    LALSnprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX, 
        " " );
  }
  else
  {
    LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
    LALSnprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
  }

  /* make sure the pointer to the first template is null */
  templateBank.snglInspiralTable = NULL;

  /* the number of nodes for a standalone job is always 1 */
  searchsumm.searchSummaryTable->nnodes = 1;

  if ( dynRangeExponent )
  {
    /* compute the dynamic range scaling for the psd computation */
    dynRange = (REAL8) pow( 2.0, dynRangeExponent );
  }
  else
  {
    dynRange = 1.0;
  }
  if ( vrbflg )
    fprintf( stdout, "using dynamic range scaling %le\n", dynRange );


  /*
   *
   * read in the input data channel
   *
   */


  /* set the time series parameters of the input data and resample params */
  memset( &resampleParams, 0, sizeof(ResampleTSParams) );
  resampleParams.deltaT = 1.0 / (REAL8) sampleRate;

  /* set the params of the input data time series */
  memset( &chan, 0, sizeof(REAL4TimeSeries) );
  memset( &geoChan, 0, sizeof(REAL8TimeSeries) );
  chan.epoch = gpsStartTime;
  chan.epoch.gpsSeconds -= padData; /* subtract pad seconds from start */

  /* copy the start time into the GEO time series */
  geoChan.epoch = chan.epoch;

  if ( globFrameData )
  {
    CHAR ifoRegExPattern[6];

    if ( vrbflg ) fprintf( stdout, "globbing for *.gwf frame files from %c "
        "of type %s in current directory\n", fqChanName[0], frInType );

    frGlobCache = NULL;

    /* create a frame cache by globbing all *.gwf files in the pwd */
    LAL_CALL( LALFrCacheGenerate( &status, &frGlobCache, NULL, NULL ), 
        &status );

    /* check we globbed at least one frame file */
    if ( ! frGlobCache->numFrameFiles )
    {
      fprintf( stderr, "error: no frame file files of type %s found\n",
          frInType );
      exit( 1 );
    }

    /* sieve out the requested data type */
    memset( &sieve, 0, sizeof(FrCacheSieve) );
    LALSnprintf( ifoRegExPattern, 
        sizeof(ifoRegExPattern) / sizeof(*ifoRegExPattern), ".*%c.*", 
        fqChanName[0] );
    sieve.srcRegEx = ifoRegExPattern;
    sieve.dscRegEx = frInType;
    LAL_CALL( LALFrCacheSieve( &status, &frInCache, frGlobCache, &sieve ), 
        &status );

    /* check we got at least one frame file back after the sieve */
    if ( ! frInCache->numFrameFiles )
    {
      fprintf( stderr, "error: no frame files of type %s globbed as input\n",
          frInType );
      exit( 1 );
    }

    LAL_CALL( LALDestroyFrCache( &status, &frGlobCache ), &status );
  }
  else
  {
    if ( vrbflg ) fprintf( stdout, 
        "reading frame file locations from cache file: %s\n", frInCacheName );

    /* read a frame cache from the specified file */
    LAL_CALL( LALFrCacheImport( &status, &frInCache, frInCacheName), &status );
  }

  /* open the input data frame stream from the frame cache */
  LAL_CALL( LALFrCacheOpen( &status, &frStream, frInCache ), &status );

  /* set the mode of the frame stream to fail on gaps or time errors */
  frStream->mode = LAL_FR_VERBOSE_MODE;

  /* seek to required epoch and set chan name */
  LAL_CALL( LALFrSeek( &status, &(chan.epoch), frStream ), &status );
  frChan.name = fqChanName;

  if ( calData == real_8 )
  {
    /* determine the sample rate of the raw data */
    LAL_CALL( LALFrGetREAL8TimeSeries( &status, &geoChan, &frChan, frStream ),
        &status );

    /* copy the data paramaters from the GEO channel to input data channel */
    LALSnprintf( chan.name, LALNameLength * sizeof(CHAR), "%s", geoChan.name );
    chan.epoch          = geoChan.epoch;
    chan.deltaT         = geoChan.deltaT;
    chan.f0             = geoChan.f0;
    chan.sampleUnits    = geoChan.sampleUnits;
  }
  else
  {
    /* determine the sample rate of the raw data and allocate enough memory */
    LAL_CALL( LALFrGetREAL4TimeSeries( &status, &chan, &frChan, frStream ),
        &status );
  }

  /* store the input sample rate */
  this_search_summvar = searchsummvars.searchSummvarsTable = 
    (SearchSummvarsTable *) LALCalloc( 1, sizeof(SearchSummvarsTable) );
  LALSnprintf( this_search_summvar->name, LIGOMETA_NAME_MAX * sizeof(CHAR),
      "raw data sample rate" );
  this_search_summvar->value = chan.deltaT;

  /* determine if we need to resample the channel */
  if ( vrbflg )
  {
    fprintf( stdout, "resampleParams.deltaT = %e\n", resampleParams.deltaT );
    fprintf( stdout, "chan.deltaT = %e\n", chan.deltaT );
  }
  if ( ! ( fabs( resampleParams.deltaT - chan.deltaT ) < epsilon ) )
  {
    resampleChan = 1;
    if ( vrbflg )
      fprintf( stdout, "input channel will be resampled\n" );

    if ( resampFiltType == 0 )
    {
      resampleParams.filterType = LDASfirLP;
    }
    else if ( resampFiltType == 1 )
    {
      resampleParams.filterType = defaultButterworth;
    }
  }

  /* determine the number of points to get and create storage forr the data */
  inputLengthNS = (REAL8) ( 1000000000LL * 
      ( gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds + 2 * padData ) );
  chan.deltaT *= 1.0e9;
  numInputPoints = (UINT4) floor( inputLengthNS / chan.deltaT + 0.5 );
  if ( calData == real_8 )
  {
    /* create storage for the GEO input data */
    LAL_CALL( LALDCreateVector( &status, &(geoChan.data), numInputPoints ), 
        &status );
  }
  LAL_CALL( LALSCreateVector( &status, &(chan.data), numInputPoints ), 
      &status );

  if ( vrbflg ) fprintf( stdout, "input channel %s has sample interval "
      "(deltaT) = %e\nreading %d points from frame stream\n", fqChanName, 
      chan.deltaT / 1.0e9, numInputPoints );

  if ( calData == real_8 )
  {
    /* read in the GEO data here */
    PassBandParamStruc geoHighpassParam;

    /* read the GEO data from the time series into geoChan      */
    /* which already has the correct amount of memory allocated */
    if ( vrbflg ) fprintf( stdout, "reading GEO data from frames... " );

    LAL_CALL( LALFrGetREAL8TimeSeries( &status, &geoChan, &frChan, frStream ),
        &status);

    if ( vrbflg ) fprintf( stdout, "done\n" );

    /* high pass the GEO data using the parameters specified on the cmd line */
    geoHighpassParam.nMax = geoHighPassOrder;
    geoHighpassParam.f1 = -1.0;
    geoHighpassParam.f2 = (REAL8) geoHighPassFreq;
    geoHighpassParam.a1 = -1.0;
    geoHighpassParam.a2 = (REAL8)(1.0 - geoHighPassAtten);
    if ( vrbflg ) fprintf( stdout, "applying %d order high pass to GEO data: "
        "%3.2f of signal passes at %4.2f Hz\n", 
        geoHighpassParam.nMax, geoHighpassParam.a2, geoHighpassParam.f2 );

    LAL_CALL( LALButterworthREAL8TimeSeries( &status, &geoChan, 
          &geoHighpassParam ), &status );

    /* cast the GEO data to REAL4 in the chan time series       */
    /* which already has the correct amount of memory allocated */
    for ( j = 0 ; j < numInputPoints ; ++j )
    {
      chan.data->data[j] = (REAL4) ( geoChan.data->data[j] * dynRange );
    }

    /* re-copy the data paramaters from the GEO channel to input data channel */
    LALSnprintf( chan.name, LALNameLength * sizeof(CHAR), "%s", geoChan.name );
    chan.epoch          = geoChan.epoch;
    chan.deltaT         = geoChan.deltaT;
    chan.f0             = geoChan.f0;
    chan.sampleUnits    = geoChan.sampleUnits;

    /* free the REAL8 GEO input data */
    LAL_CALL( LALDDestroyVector( &status, &(geoChan.data) ), &status );
    geoChan.data = NULL;
  }
  else
  {
    /* read the data channel time series from frames */
    LAL_CALL( LALFrGetREAL4TimeSeries( &status, &chan, &frChan, frStream ),
        &status );

    if ( calData == real_4 )
    {
      /* multiply the input data by dynRange */
      for ( j = 0 ; j < numInputPoints ; ++j )
      {
        chan.data->data[j] *= dynRange;
      }
    } 
  }
  memcpy( &(chan.sampleUnits), &lalADCCountUnit, sizeof(LALUnit) );

  /* store the start and end time of the raw channel in the search summary */
  searchsumm.searchSummaryTable->in_start_time = chan.epoch;
  LAL_CALL( LALGPStoFloat( &status, &tsLength, &(chan.epoch) ), 
      &status );
  tsLength += chan.deltaT * (REAL8) chan.data->length;
  LAL_CALL( LALFloatToGPS( &status, 
        &(searchsumm.searchSummaryTable->in_end_time), &tsLength ), &status );

  /* close the frame file stream and destroy the cache */
  LAL_CALL( LALFrClose( &status, &frStream ), &status );
  LAL_CALL( LALDestroyFrCache( &status, &frInCache ), &status );

  /* write the raw channel data as read in from the frame files */
  if ( writeRawData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
      &chan, "ct", "RAW" );

  if ( vrbflg ) fprintf( stdout, "read channel %s from frame stream\n"
      "got %d points with deltaT %e\nstarting at GPS time %d sec %d ns\n", 
      chan.name, chan.data->length, chan.deltaT, 
      chan.epoch.gpsSeconds, chan.epoch.gpsNanoSeconds );

  /* resample the input data */
  if ( resampleChan )
  {
    if (vrbflg) fprintf( stdout, "resampling input data from %e to %e\n",
        chan.deltaT, resampleParams.deltaT );

    LAL_CALL( LALResampleREAL4TimeSeries( &status, &chan, &resampleParams ),
        &status );

    if ( vrbflg ) fprintf( stdout, "channel %s resampled:\n"
        "%d points with deltaT %e\nstarting at GPS time %d sec %d ns\n", 
        chan.name, chan.data->length, chan.deltaT, 
        chan.epoch.gpsSeconds, chan.epoch.gpsNanoSeconds );

    /* write the resampled channel data as read in from the frame files */
    if ( writeRawData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
        &chan, "ct", "RAW_RESAMP" );
  }

  /* store the filter data sample rate */
  this_search_summvar = this_search_summvar->next = 
    (SearchSummvarsTable *) LALCalloc( 1, sizeof(SearchSummvarsTable) );
  LALSnprintf( this_search_summvar->name, LIGOMETA_NAME_MAX * sizeof(CHAR),
      "filter data sample rate" );
  this_search_summvar->value = chan.deltaT;

  /*
   * 
   * compute a calibrated strain spectrum
   *
   */

  /* create storage for the response and spectrum */
  memset( &spec, 0, sizeof(REAL4FrequencySeries) );
  LAL_CALL( LALSCreateVector( &status, &(spec.data), numPoints / 2 + 1 ), 
      &status );
  memset( &resp, 0, sizeof(COMPLEX8FrequencySeries) );
  LAL_CALL( LALCCreateVector( &status, &(resp.data), numPoints / 2 + 1 ), 
      &status );
  resp.epoch = spec.epoch = gpsStartTime;

  /* iir filter to remove low frequencies from data channel */
  if ( highPass )
  {
    PassBandParamStruc highpassParam;
    highpassParam.nMax = highPassOrder;
    highpassParam.f1 = -1.0;
    highpassParam.f2 = (REAL8) highPassFreq;
    highpassParam.a1 = -1.0;
    highpassParam.a2 = (REAL8)(1.0 - highPassAtten); /* a2 is not attenuation */

    if ( vrbflg ) fprintf( stdout, "applying %d order high pass: "
        "%3.2f of signal passes at %4.2f Hz\n", 
        highpassParam.nMax, highpassParam.a2, highpassParam.f2 );

    LAL_CALL( LALDButterworthREAL4TimeSeries( &status, &chan, &highpassParam ),
        &status );
  }

  /* remove pad from requested data from start and end of time series */
  memmove( chan.data->data, chan.data->data + padData * sampleRate, 
      (chan.data->length - 2 * padData * sampleRate) * sizeof(REAL4) );
  LALRealloc( chan.data->data, 
      (chan.data->length - 2 * padData * sampleRate) * sizeof(REAL4) );
  chan.data->length -= 2 * padData * sampleRate;
  chan.epoch.gpsSeconds += padData;

  if ( vrbflg ) fprintf( stdout, "after removal of %d second padding at "
      "start and end:\ndata channel sample interval (deltaT) = %e\n"
      "data channel length = %d\nstarting at %d sec %d ns\n", 
      padData , chan.deltaT , chan.data->length, 
      chan.epoch.gpsSeconds, chan.epoch.gpsNanoSeconds );

  /* store the start and end time of the filter channel in the search summ */
  searchsumm.searchSummaryTable->out_start_time = chan.epoch;
  LAL_CALL( LALGPStoFloat( &status, &tsLength, &(chan.epoch) ), 
      &status );
  tsLength += chan.deltaT * (REAL8) chan.data->length;
  LAL_CALL( LALFloatToGPS( &status, 
        &(searchsumm.searchSummaryTable->out_end_time), &tsLength ), 
      &status );

  /* compute the windowed power spectrum for the data channel */
  avgSpecParams.window = NULL;
  avgSpecParams.plan = NULL;
  LAL_CALL( LALCreateForwardRealFFTPlan( &status, 
        &(avgSpecParams.plan), numPoints, 0 ), &status );
  switch ( specType )
  {
    case 0:
      avgSpecParams.method = useMean;
      if ( vrbflg ) fprintf( stdout, "computing mean psd" );
      break;
    case 1:
      avgSpecParams.method = useMedian;
      if ( vrbflg ) fprintf( stdout, "computing median psd" );
      break;
  }

  wpars.type = Hann;
  wpars.length = numPoints;
  avgSpecParams.overlap = numPoints / 2;
  if ( vrbflg ) 
    fprintf( stdout, " with overlap %d\n", avgSpecParams.overlap );

  LAL_CALL( LALCreateREAL4Window( &status, &(avgSpecParams.window),
        &wpars ), &status );
  LAL_CALL( LALREAL4AverageSpectrum( &status, &spec, &chan, &avgSpecParams ),
      &status );
  LAL_CALL( LALDestroyREAL4Window( &status, &(avgSpecParams.window) ), 
      &status );
  LAL_CALL( LALDestroyRealFFTPlan( &status, &(avgSpecParams.plan) ), &status );

  /* write the spectrum data to a file */
  if ( writeSpectrum )
  {
    strcpy( spec.name, chan.name );
    outFrame = fr_add_proc_REAL4FrequencySeries( outFrame, 
        &spec, "ct/sqrtHz", "PSD" );
  }

  /* set the parameters of the response to match the data and spectrum */
  resp.deltaF = spec.deltaF;
  resp.f0 = spec.f0;
  resp.sampleUnits = strainPerCount;

  if ( calData )
  {
    /* if we are using calibrated data set the response to unity */
    if ( vrbflg ) fprintf( stdout, "generating unity response function\n" );
    for( k = 0; k < resp.data->length; ++k )
    {
      resp.data->data[k].re = (REAL4) (1.0 / dynRange);
      resp.data->data[k].im = 0.0;
    }
  }
  else
  {
    /* initialize the calfacts */
    memset( &calfacts, 0, sizeof(CalibrationUpdateParams) );

    if ( pointCal )
    {
      calfacts.duration.gpsSeconds = 1; 
      calfacts.duration.gpsNanoSeconds = 0;
    }
    else
    { 
      calfacts.duration.gpsSeconds = gpsEndTime.gpsSeconds 
        - gpsStartTime.gpsSeconds;
    }
    calfacts.ifo = ifo;

    /* create the lal calibration frame cache */
    if ( globCalData )
    {
      calGlobPattern = (CHAR *) LALCalloc( calGlobLen, sizeof(CHAR) );
      LALSnprintf( calGlobPattern, calGlobLen * sizeof(CHAR), 
          "*%c*CAL*.gwf", fqChanName[0] );
      if ( vrbflg ) fprintf( stdout, "globbing for %s calibration frame files "
          "in current directory\n", calGlobPattern );
    }
    else
    {
      calGlobPattern = NULL;
      if ( vrbflg ) fprintf( stdout, 
          "reading calibration data from cache: %s\n", calCacheName );
    }

    LAL_CALL( LALCreateCalibFrCache( &status, &calCache, calCacheName, 
          NULL, calGlobPattern ), &status );

    if ( calGlobPattern ) LALFree( calGlobPattern );

    /* store the name of the calibration files used */
    for ( i = 0; i < calCache->numFrameFiles; ++i )
    {
      this_search_summvar = this_search_summvar->next = 
        (SearchSummvarsTable *) LALCalloc( 1, sizeof(SearchSummvarsTable) );
      LALSnprintf( this_search_summvar->name, LIGOMETA_NAME_MAX * sizeof(CHAR),
          "calibration frame %d", i );
      LALSnprintf( this_search_summvar->string, 
          LIGOMETA_STRING_MAX * sizeof(CHAR), "%s", 
          calCache->frameFiles[i].url );
    }

    /* generate the response function for the current time */
    if ( vrbflg ) fprintf( stdout, "generating response at time %d sec %d ns\n"
        "response parameters f0 = %e, deltaF = %e, length = %d\n",
        resp.epoch.gpsSeconds, resp.epoch.gpsNanoSeconds,
        resp.f0, resp.deltaF, resp.data->length );
    LAL_CALL( LALExtractFrameResponse( &status, &resp, calCache, 
          &calfacts ), &status );

    /* descroy the frame cache for the calibrated data */
    LAL_CALL( LALDestroyFrCache( &status, &calCache ), &status );

    if ( vrbflg ) fprintf( stdout, "Values of calibration coefficients \n"
        "alpha = %f, alpha_beta = %f\n",
        calfacts.alpha.re, calfacts.alphabeta.re );
  }

  /* write the calibration data to a file */
  if ( writeResponse )
  {
    strcpy( resp.name, chan.name );
    outFrame = fr_add_proc_COMPLEX8FrequencySeries( outFrame, 
        &resp, "strain/ct", "RESPONSE" );
  }

  /* set low frequency cutoff of power spectrum */
  cut = fLow / spec.deltaF > 1 ?  fLow / spec.deltaF : 1;

  /* compute a calibrated strain power spectrum */
  bankIn.shf.epoch = spec.epoch;
  memcpy( bankIn.shf.name, spec.name, LALNameLength * sizeof(CHAR) );
  bankIn.shf.deltaF = spec.deltaF;
  bankIn.shf.f0 = spec.f0;
  bankIn.shf.data = NULL;
  pair.unitOne = &(spec.sampleUnits);
  pair.unitTwo = &(resp.sampleUnits);
  LAL_CALL( LALUnitMultiply( &status, &(bankIn.shf.sampleUnits), &pair ), 
      &status );
  LAL_CALL( LALDCreateVector( &status, &(bankIn.shf.data), spec.data->length ),
      &status );
  memset( bankIn.shf.data->data, 0, 
      bankIn.shf.data->length * sizeof(COMPLEX8) ); 

  shf = spec.data->data[cut] * 
    ( resp.data->data[cut].re * resp.data->data[cut].re +
      resp.data->data[cut].im * resp.data->data[cut].im );
  for ( k = 1; k < cut ; ++k )
  {
    bankIn.shf.data->data[k] = shf;
  }
  for ( k = cut; k < bankIn.shf.data->length; ++k )
  {
    respRe = (REAL8) resp.data->data[k].re;
    respIm = (REAL8) resp.data->data[k].im;
    bankIn.shf.data->data[k] = (REAL8) spec.data->data[k] *
      ( respRe * respRe + respIm * respIm );
  }

  /* write the scaled strain spectrum data to a file */
  if ( writeStrainSpec )
  {
#if 0
    strcpy( spec.name, chan.name );
    outFrame = fr_add_proc_REAL8FrequencySeries( outFrame, 
        &(bankIn.shf), "strain/sqrtHz", "STRAIN_PSD" );
#endif
    LALSnprintf( fname, sizeof(fname), "%s-TMPLTBANK-%d-%d.strainspec.txt",
        ifo, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
    LALDPrintFrequencySeries( &(bankIn.shf), fname );
  }


  /*
   *
   * compute the standard candle distance
   *
   */


  if ( computeCandle )
  {
    CHAR  candleComment[LIGOMETA_SUMMVALUE_COMM_MAX];
    REAL8 distance = 0;

    while ( candleMass1 < 50.0 )
    {
      /* experimental code to ease the computation of the standard candle */
      distance = compute_candle_distance(candleMass1, candleMass2,
          candleSnr, chan.deltaT, numPoints, &(bankIn.shf), cut);

      if ( vrbflg ) fprintf( stdout, "maximum distance for (%3.2f,%3.2f) "
          "at signal-to-noise %3.2f = ", candleMass1, candleMass2, candleSnr );

      /* experimental code to populate the summValue table */
      LALSnprintf( candleComment, LIGOMETA_SUMMVALUE_COMM_MAX,
          "%3.2f_%3.2f_%3.2f", candleMass1, candleMass2, candleSnr );

      this_summvalue =
        add_summvalue_table(this_summvalue, gpsStartTime, gpsEndTime,
            PROGRAM_NAME, ifo, "inspiral_effective_distance", 
            candleComment, distance);

      if ( vrbflg ) fprintf( stdout, "%e Mpc\n", (*this_summvalue)->value );

      candleMass2 = candleMass1 = candleMass1 + 1.0;
      this_summvalue = &(*this_summvalue)->next;
    }
  }


  /*
   *
   * compute the template bank
   *
   */


  /* bank generation parameters */
  if (isMaxTotMass)
    bankIn.massRange   = MinComponentMassMaxTotalMass;
  else
    bankIn.massRange     = MinMaxComponentMass;
  bankIn.mMin          = (REAL8) minMass;
  bankIn.mMax          = (REAL8) maxMass;
  if (isMaxTotMass)
  {
    bankIn.MMax	       = (REAL8) maxMass;
    bankIn.mMax	       = bankIn.MMax - bankIn.mMin;
  }
  else
    bankIn.MMax          = bankIn.mMax * 2.0;
  bankIn.psi0Min       = (REAL8) psi0Min;
  bankIn.psi0Max       = (REAL8) psi0Max;
  bankIn.psi3Min       = (REAL8) psi3Min;
  bankIn.psi3Max       = (REAL8) psi3Max;
  bankIn.numFcutTemplates = (UINT4) maxFcutTmplts;
  bankIn.alpha         = (REAL8) alpha;
  bankIn.mmCoarse      = (REAL8) minMatch;
  bankIn.mmFine        = 0.99; /* doesn't matter since no fine bank yet */
  bankIn.fLower        = (REAL8) fLow;
  bankIn.fUpper        = (REAL8) fUpper;
  bankIn.iflso         = 0; /* currently not implemented */
  bankIn.tSampling     = (REAL8) sampleRate;
  bankIn.order         = order;
  bankIn.approximant   = approximant;
  bankIn.gridSpacing   = gridSpacing;
  bankIn.space         = space;
  bankIn.etamin        = bankIn.mMin * ( bankIn.MMax - bankIn.mMin) /
    ( bankIn.MMax * bankIn.MMax );
  bankIn.LowGM  = -2.;
  bankIn.HighGM = 6.;
  /* generate the template bank */
  if ( vrbflg )
  {
    fprintf( stdout, "generating template bank parameters... " );
    fflush( stdout );
  }
  LAL_CALL( LALInspiralBankGeneration( &status, &bankIn, &tmplt, &numCoarse),
      &status );
  if ( vrbflg )
  {
    fprintf( stdout, "done. Got %d templates\n", numCoarse );
    fflush( stdout );
  }

  if ( numCoarse )
  {
    templateBank.snglInspiralTable = tmplt;
    LALSnprintf( tmplt->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR), ifo );
    LALSnprintf( tmplt->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
        "tmpltbank" );
    LALSnprintf( tmplt->channel, LIGOMETA_CHANNEL_MAX * sizeof(CHAR),
        channelName );
    while( (tmplt = tmplt->next) )
    {
      LALSnprintf( tmplt->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR), ifo );
      LALSnprintf( tmplt->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
          "tmpltbank" );
      LALSnprintf( tmplt->channel, LIGOMETA_CHANNEL_MAX * sizeof(CHAR),
          channelName );
    }
  }

  /* save the number of templates in the search summary table */
  searchsumm.searchSummaryTable->nevents = numCoarse;


  /*
   *
   * free the data storage
   *
   */


  LAL_CALL( LALDDestroyVector( &status, &(bankIn.shf.data) ), &status );
  LAL_CALL( LALSDestroyVector( &status, &(chan.data) ), &status );
  LAL_CALL( LALSDestroyVector( &status, &(spec.data) ), &status );
  LAL_CALL( LALCDestroyVector( &status, &(resp.data) ), &status );


  /*
   *
   * write the result results to disk
   *
   */


  /* write the output frame */
  if ( writeRawData || writeResponse || writeSpectrum )
  {
    LALSnprintf( fname, sizeof(fname), "%s-TMPLTBANK-%d-%d.gwf",
        ifo, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
    frOutFile = FrFileONew( fname, 0 );
    FrameWrite( outFrame, frOutFile );
    FrFileOEnd( frOutFile );
  }

  /* open the output xml file */
  memset( &results, 0, sizeof(LIGOLwXMLStream) );
  if ( userTag )
  {
    LALSnprintf( fname, sizeof(fname), "%s-TMPLTBANK_%s-%d-%d.xml",
        ifo, userTag, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else
  {
    LALSnprintf( fname, sizeof(fname), "%s-TMPLTBANK-%d-%d.xml",
        ifo, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &results, fname), &status );

  /* write the process table */
  LALSnprintf( proctable.processTable->ifos, LIGOMETA_IFO_MAX, "%s", ifo );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->end_time),
        &accuracy ), &status );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, proctable, 
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  free( proctable.processTable );

  /* erase the first empty process params entry */
  {
    ProcessParamsTable *emptyPPtable = procparams.processParamsTable;
    procparams.processParamsTable = procparams.processParamsTable->next;
    free( emptyPPtable );
  }

  /* write the process params table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_params_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, procparams, 
        process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  while( procparams.processParamsTable )
  {
    this_proc_param = procparams.processParamsTable;
    procparams.processParamsTable = this_proc_param->next;
    free( this_proc_param );
  }

  /* write the search summary table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, 
        search_summary_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, searchsumm, 
        search_summary_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );

  /* write the search summvars table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, 
        search_summvars_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, searchsummvars, 
        search_summvars_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  while( searchsummvars.searchSummvarsTable )
  {
    this_search_summvar = searchsummvars.searchSummvarsTable;
    searchsummvars.searchSummvarsTable = this_search_summvar->next;
    LALFree( this_search_summvar );
  }

  /* write the standard candle to the file */
  if ( computeCandle )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, summ_value_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, candle, 
          summ_value_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
    LALFree( candle.summValueTable );
    candle.summValueTable = NULL;
  }

  /* write the template bank to the file */
  if ( templateBank.snglInspiralTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, sngl_inspiral_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, templateBank, 
          sngl_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  }
  while ( templateBank.snglInspiralTable )
  {
    tmplt = templateBank.snglInspiralTable;
    templateBank.snglInspiralTable = templateBank.snglInspiralTable->next;
    LALFree( tmplt );
  }

  /* close the output xml file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &results ), &status );

  /* free the rest of the memory, check for memory leaks and exit */
  if ( calCacheName ) free( calCacheName );
  if ( frInCacheName ) free( frInCacheName );
  if ( frInType ) free( frInType );
  if ( channelName ) free( channelName );
  if ( fqChanName ) free( fqChanName );
  LALCheckMemoryLeaks();

  /* print a success message to stdout for parsing by exitcode */
  fprintf( stdout, "%s: EXITCODE0\n", argv[0] );
  exit( 0 );
}

/* ------------------------------------------------------------------------- */

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  calloc( 1, sizeof(ProcessParamsTable) ); \
  LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
   PROGRAM_NAME ); \
   LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
     long_options[option_index].name ); \
     LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
     LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

#define USAGE \
"  --help                       display this message\n"\
"  --verbose                    print progress information\n"\
"  --version                    print version information and exit\n"\
"  --debug-level LEVEL          set the LAL debug level to LEVEL\n"\
"  --user-tag STRING            set the process_params usertag to STRING\n"\
"  --comment STRING             set the process table comment to STRING\n"\
"\n"\
"  --gps-start-time SEC         GPS second of data start time\n"\
"  --gps-end-time SEC           GPS second of data end time\n"\
"  --pad-data T                 pad the data start and end time by T seconds\n"\
"\n"\
"  --glob-frame-data            glob *.gwf files in the pwd to obtain frame data\n"\
"  --frame-type TAG             input data is contained in frames of type TAG\n"\
"  --frame-cache                obtain frame data from LAL frame cache FILE\n"\
"  --calibration-cache FILE     obtain calibration from LAL frame cache FILE\n"\
"  --glob-calibration-data      obtain calibration by globbing in working dir\n"\
"\n"\
"  --channel-name CHAN          read data from interferometer channel CHAN\n"\
"  --calibrated-data TYPE       calibrated data of TYPE real_4 or real_8\n"\
"  --geo-high-pass-freq F       high pass GEO data above F Hz using an IIR filter\n"\
"  --geo-high-pass-order O      set the order of the GEO high pass filter to O\n"\
"  --geo-high-pass-atten A      set the attenuation of the high pass filter to A\n"\
"  --point-calibration          use the first point in the chunk to calibrate\n"\
"\n"\
"  --sample-rate F              filter data at F Hz, downsampling if necessary\n"\
"  --resample-filter TYPE       set resample filter to TYPE [ldas|butterworth]\n"\
"\n"\
"  --disable-high-pass          turn off the IIR highpass filter\n"\
"  --enable-high-pass F         high pass data above F Hz using an IIR filter\n"\
"  --high-pass-order O          set the order of the high pass filter to O\n"\
"  --high-pass-attenuation A    set the attenuation of the high pass filter to A\n"\
"  --spectrum-type TYPE         use PSD estimator TYPE [mean|median]\n"\
"  --dynamic-range-exponent X   set dynamic range scaling to 2^X\n"\
"\n"\
"  --segment-length N           set data segment length to N points\n"\
"  --number-of-segments N       set number of data segments to N\n"\
"\n"\
"  --standard-candle            compute a standard candle from the PSD\n"\
"  --candle-snr SNR             signal-to-noise ration of standard candle\n"\
"  --candle-mass1 M             mass of first component in candle binary\n"\
"  --candle-mass2 M             mass of second component in candle binary\n"\
"\n"\
"  --low-frequency-cutoff F     do not filter below F Hz\n"\
"  --high-frequency-cutoff F    upper frequency cutoff in Hz\n"\
"\n"\
"  --minimum-mass MASS          set minimum component mass of bank to MASS\n"\
"  --maximum-mass MASS          set maximum component mass of bank to MASS\n"\
"  --maximum-total-mass         make --maximum-mass refer to total mass\n"\
"\n"\
"  --minimum-psi0 PSI0          set minimum range of BCV parameter psi0 to PSI0\n"\
"  --maximum-psi0 PSI0          set maximum range of BCV parameter psi0 to PSI0\n"\
"  --minimum-psi3 PSI3          set minimum range of BCV parameter psi3 to PSI3\n"\
"  --maximum-psi3 PSI3          set maximum range of BCV parameter psi3 to PSI3\n"\
"  --maximum-fcut-tmplts N      maximum number of tmplts in fcut direction is N\n"\
"  --alpha ALPHA                set alpha for the BCV bank generation\n"\
"\n"\
"  --minimal-match M            generate bank with minimal match M\n"\
"\n"\
"  --order ORDER                set post-Newtonian order of the waveform to ORDER\n"\
"                                 (newtonian|oneHalfPN|onePN|onePointFivePN|\n"\
"                                 twoPN|twoPointFive|threePN|threePointFivePN)\n"\
"  --approximant APPROX         set approximant of the waveform to APPROX\n"\
"                                 (TaylorT1|TaylorT2|TaylorT3|TaylorF1|TaylorF2|\n"\
"                                 PadeT1|PadeT2|EOB|BCV|SpinTaylorT3|BCVSpin)\n"\
"  --space SPACE                grid up template bank with mass parameters SPACE\n"\
"                                 (Tau0Tau2|Tau0Tau3|Psi0Psi3)\n"\
"  --grid-spacing GRIDSPACING   grid up template bank with GRIDSPACING\n"\
"                                 (Square|Hexagonal|SquareNotOriented|\n"\
"                                 HexagonalNotOriented)\n"\
"\n"\
"  --write-raw-data             write raw data to a frame file\n"\
"  --write-response             write the computed response function to a frame\n"\
"  --write-spectrum             write the uncalibrated psd to a frame\n"\
"  --write-strain-spectrum      write the calibrated strain psd to a text file\n"\
"\n"

int arg_parse_check( int argc, char *argv[], MetadataTable procparams )
{
  /* getopt arguments */
  struct option long_options[] =
  {
    /* these options set a flag */
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"disable-high-pass",       no_argument,       &highPass,         0 },
    {"standard-candle",         no_argument,       &computeCandle,    1 },
    {"glob-frame-data",         no_argument,       &globFrameData,    1 },
    {"glob-calibration-data",   no_argument,       &globCalData,      1 },
    {"point-calibration",       no_argument,       &pointCal,         1 },
    /* parameters used to generate calibrated power spectrum */
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"channel-name",            required_argument, 0,                'c'},
    {"segment-length",          required_argument, 0,                'd'},
    {"number-of-segments",      required_argument, 0,                'e'},
    {"sample-rate",             required_argument, 0,                'g'},
    {"calibrated-data",         required_argument, 0,                'M'},
    {"geo-high-pass-freq",      required_argument, 0,                'J'},
    {"geo-high-pass-order",     required_argument, 0,                'K'},
    {"geo-high-pass-atten",     required_argument, 0,                'L'},
    {"help",                    no_argument,       0,                'h'},
    {"low-frequency-cutoff",    required_argument, 0,                'i'},
    {"spectrum-type",           required_argument, 0,                'j'},
    {"dynamic-range-exponent",  required_argument, 0,                'f'},
    {"calibration-cache",       required_argument, 0,                'p'},
    {"comment",                 required_argument, 0,                's'},
    {"enable-high-pass",        required_argument, 0,                't'},
    {"high-pass-order",         required_argument, 0,                'H'},
    {"high-pass-attenuation",   required_argument, 0,                'I'},
    {"frame-cache",             required_argument, 0,                'u'},
    {"frame-type",              required_argument, 0,                'n'},
    {"pad-data",                required_argument, 0,                'x'},
    {"debug-level",             required_argument, 0,                'z'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {"version",                 no_argument,       0,                'V'},    
    {"resample-filter",         required_argument, 0,                'r'},
    /* template bank generation parameters */
    {"minimum-mass",            required_argument, 0,                'A'},
    {"maximum-mass",            required_argument, 0,                'B'},
    {"minimum-psi0",            required_argument, 0,                'P'},
    {"maximum-psi0",            required_argument, 0,                'Q'},
    {"minimum-psi3",            required_argument, 0,                'R'},
    {"maximum-psi3",            required_argument, 0,                'S'},
    {"maximum-fcut-tmplts",     required_argument, 0,                'U'},
    {"alpha",                   required_argument, 0,                'T'},
    {"minimal-match",           required_argument, 0,                'C'},
    {"high-frequency-cutoff",   required_argument, 0,                'D'},
    {"order",                   required_argument, 0,                'E'},
    {"approximant",             required_argument, 0,                'F'},
    {"space",                   required_argument, 0,                'G'},
    {"grid-spacing",            required_argument, 0,                'v'},
    {"maximum-total-mass", 	no_argument, 	   0,		     'y'},
    /* standard candle parameters */
    {"candle-snr",              required_argument, 0,                'k'},
    {"candle-mass1",            required_argument, 0,                'l'},
    {"candle-mass2",            required_argument, 0,                'm'},
    /* frame writing options */
    {"write-raw-data",          no_argument,       &writeRawData,     1 },
    {"write-response",          no_argument,       &writeResponse,    1 },
    {"write-spectrum",          no_argument,       &writeSpectrum,    1 },
    {"write-strain-spectrum",   no_argument,       &writeStrainSpec,  1 },
    {0, 0, 0, 0}
  };
  int c;
  ProcessParamsTable *this_proc_param = procparams.processParamsTable;
  UINT4   haveOrder       = 0;
  UINT4   haveApprox      = 0;
  UINT4   haveSpace       = 0;
  UINT4   havePsi0Min     = 0;
  UINT4   havePsi0Max     = 0;
  UINT4   havePsi3Min     = 0;
  UINT4   havePsi3Max     = 0;
  UINT4   haveAlpha       = 0;


  /*
   *
   * parse command line arguments
   *
   */


  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "a:b:c:d:e:f:g:hi:j:k:l:m:n:p:r:s:t:u:v:x:yz:"
        "A:B:C:D:E:F:G:H:I:J:K:L:M:P:Q:R:S:T:U:VZ:",
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
          gpsStartTime.gpsSeconds = (INT4) gstartt;
          gpsStartTime.gpsNanoSeconds = 0;
          ADD_PROCESS_PARAM( "int", "%ld", gstartt );
        }
        break;

      case 'b':
        {
          long int gendt = atol( optarg );
          if ( gendt > 999999999 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS end time is after " 
                "Sep 14, 2011  01:46:26 UTC:\n"
                "(%ld specified)\n", 
                long_options[option_index].name, gendt );
            exit( 1 );
          }
          else if ( gendt < 441417609 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS end time is prior to " 
                "Jan 01, 1994  00:00:00 UTC:\n"
                "(%ld specified)\n", 
                long_options[option_index].name, gendt );
            exit( 1 );
          }            
          gpsEndTime.gpsSeconds = (INT4) gendt;
          gpsEndTime.gpsNanoSeconds = 0;
          ADD_PROCESS_PARAM( "int", "%ld", gendt );
        }
        break;

      case 'c':
        {
          /* create storage for the channel name and copy it */
          char *channamptr = NULL;
          optarg_len = strlen( optarg ) + 1;
          fqChanName = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
          memcpy( fqChanName, optarg, optarg_len );
          ADD_PROCESS_PARAM( "string", "%s", optarg );

          /* check that we have a proper channel name */
          if ( ! (channamptr = strstr( fqChanName, ":" ) ) )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "channel name must be a full LIGO channel name "
                "e.g. L1:LSC-AS_Q\n(%s specified)\n",
                long_options[option_index].name, optarg );
            exit( 1 );
          }
          optarg_len = strlen( ++channamptr ) + 1;
          channelName = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
          memcpy( channelName, channamptr, optarg_len );

          /* copy the first two characters to ifo */
          memset( ifo, 0, sizeof(ifo) );
          memcpy( ifo, optarg, sizeof(ifo) - 1 );
        }
        break;

      case 'd':
        numPoints = (INT4) atoi( optarg );
        if ( numPoints < 2 || numPoints % 2 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "number of points must be a non-zero power of 2: "
              "(%d specified) \n", 
              long_options[option_index].name, numPoints );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", numPoints );
        break;

      case 'e':
        numSegments = (INT4) atoi( optarg );
        if ( numSegments < 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "number of data segment must be greater than 0: "
              "(%d specified)\n", 
              long_options[option_index].name, numSegments );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", numSegments );
        break;

      case 'g':
        sampleRate = (INT4) atoi( optarg );
        if ( sampleRate < 2 || sampleRate > 16384 || sampleRate % 2 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "rate must be power of 2 between 2 and 16384 inclusive: "
              "(%d specified)\n", 
              long_options[option_index].name, sampleRate );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", sampleRate );
        break;

      case 'M':
        /* specify which type of calibrated data */
        {
          if ( ! strcmp( "real_4", optarg ) )
          {
            calData = real_4;
          }
          else if ( ! strcmp( "real_8", optarg ) )
          {
            calData = real_8;
          }
          else
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "unknown data type specified;\n"
                "%s (must be one of: real_4, real_8)\n",
                long_options[option_index].name, optarg);
          }
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'J':
        geoHighPassFreq = (REAL4) atof( optarg );
        if ( geoHighPassFreq <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GEO high pass filter frequency must be greater than 0 Hz: "
              "(%f Hz specified)\n",
              long_options[option_index].name, geoHighPassFreq );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", geoHighPassFreq );
        break;

      case 'K':
        geoHighPassOrder = (INT4) atoi( optarg );
        if ( geoHighPassOrder <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GEO high pass filter order must be greater than 0: "
              "(%d specified)\n",
              long_options[option_index].name, geoHighPassOrder );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", geoHighPassOrder );
        break;

      case 'L':
        geoHighPassAtten = (REAL4) atof( optarg );
        if ( geoHighPassAtten < 0.0 || geoHighPassAtten > 1.0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GEO high pass attenuation must be in the range [0:1]: "
              "(%f specified)\n",
              long_options[option_index].name, geoHighPassAtten );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", geoHighPassAtten );
        break;

      case 'h':
        fprintf( stdout, USAGE );
        exit( 0 );
        break;

      case 'i':
        fLow = (REAL4) atof( optarg );
        if ( fLow < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "low frequency cutoff is less than 0 Hz: "
              "(%f Hz specified)\n",
              long_options[option_index].name, fLow );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", fLow );
        break;

      case 'j':
        if ( ! strcmp( "mean", optarg ) )
        {
          specType = 0;
        }
        else if ( ! strcmp( "median", optarg ) )
        {
          specType = 1;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown power spectrum type: "
              "%s (must be mean or median)\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'f':
        dynRangeExponent = (REAL4) atof( optarg );
        ADD_PROCESS_PARAM( "float", "%e", dynRangeExponent );
        break;

      case 'p':
        /* create storage for the calibration frame cache name */
        optarg_len = strlen( optarg ) + 1;
        calCacheName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( calCacheName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'r':
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

      case 's':
        if ( strlen( optarg ) > LIGOMETA_COMMENT_MAX - 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "comment must be less than %d characters\n",
              long_options[option_index].name, LIGOMETA_COMMENT_MAX );
          exit( 1 );
        }
        else
        {
          LALSnprintf( comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
        }
        break;

      case 't':
        highPass = 1;
        highPassFreq = (REAL4) atof( optarg );
        if ( highPassFreq < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "low frequency cutoff is less than 0 Hz: "
              "(%f Hz specified)\n",
              long_options[option_index].name, highPassFreq );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", highPassFreq );
        break;

      case 'H':
        highPassOrder = (INT4) atoi( optarg );
        if ( highPassOrder <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "high pass filter order must be greater than 0: "
              "(%d specified)\n",
              long_options[option_index].name, highPassOrder );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", highPassOrder );
        break;

      case 'I':
        highPassAtten = (REAL4) atof( optarg );
        if ( highPassAtten < 0.0 || highPassAtten > 1.0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "high pass attenuation must be in the range [0:1]: "
              "(%f specified)\n",
              long_options[option_index].name, highPassAtten );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", highPassAtten );
        break;

      case 'u':
        optarg_len = strlen( optarg ) + 1;
        frInCacheName = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( frInCacheName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'n':
        optarg_len = strlen( optarg ) + 1;
        frInType = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( frInType, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'x':
        padData = (UINT4) atoi( optarg );
        if ( padData < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "number of seconds to pad from input data"
              "must be greater than 0: (%d specified)\n", 
              long_options[option_index].name, padData );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", padData );
        break;

      case 'z':
        set_debug_level( optarg );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'Z':
        /* create storage for the usertag */
        optarg_len = strlen( optarg ) + 1;
        userTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( userTag, optarg, optarg_len );

        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
            PROGRAM_NAME );
        LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "-userTag" );
        LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            optarg );
        break;

      case 'A':
        minMass = (REAL4) atof( optarg );
        if ( minMass <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "miniumum component mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, minMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", minMass );
        break;

      case 'B':
        maxMass = (REAL4) atof( optarg );
        if ( maxMass <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "maxiumum component mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, maxMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", maxMass );
        break;

      case 'P':
        psi0Min = (REAL4) atof( optarg );
        if ( psi0Min <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "miniumum value of psi0 must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, psi0Min );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", psi0Min );
        havePsi0Min = 1;
        break;

      case 'Q':
        psi0Max = (REAL4) atof( optarg );
        if ( psi0Max <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "maximum value of psi0 must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, psi0Max );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", psi0Max );
        havePsi0Max = 1;
        break;

      case 'R':
        psi3Min = (REAL4) atof( optarg );
        if ( psi0Min <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "miniumum value of psi3 must be < 0: "
              "(%f specified)\n",
              long_options[option_index].name, psi3Min );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", psi3Min );
        havePsi3Min = 1;
        break;

      case 'S':
        psi3Max = (REAL4) atof( optarg );
        ADD_PROCESS_PARAM( "float", "%e", psi3Max );
        havePsi3Max = 1;
        break;

      case 'U':
        maxFcutTmplts = (INT4) atof( optarg );
        if ( maxFcutTmplts < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "number of templates in f_final direction must be >= 0"
              "(%d specified)\n",
              long_options[option_index].name, maxFcutTmplts );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", maxFcutTmplts );
        break;

      case 'T':
        alpha = (REAL4) atof( optarg );
        if ( alpha < -1 || alpha > 1 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "value of alpha must be the range [0:1]"
              "(%f specified)\n",
              long_options[option_index].name, alpha );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", alpha );
        haveAlpha = 1;
        break;

      case 'C':
        minMatch = (REAL4) atof( optarg );
        if ( minMatch <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "minimum match of bank must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, minMatch );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", minMatch );
        break;

      case 'D':
        fUpper = (REAL4) atof( optarg );
        if ( fUpper <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "miniumu component mass must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, fUpper );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", fUpper );
        break;

      case 'E':
        if ( ! strcmp( "newtonian", optarg ) )
        {
          order = newtonian;
        }
        else if ( ! strcmp( "oneHalfPN", optarg ) )
        {
          order = oneHalfPN;
        }
        else if ( ! strcmp( "onePN", optarg ) )
        {
          order = onePN;
        }
        else if ( ! strcmp( "onePointFivePN", optarg ) )
        {
          order = onePointFivePN;
        }
        else if ( ! strcmp( "twoPN", optarg ) )
        {
          order = twoPN;
        }
        else if ( ! strcmp( "twoPointFive", optarg ) )
        {
          order = twoPointFivePN;
        }
        else if ( ! strcmp( "threePN", optarg ) )
        {
          order = threePN;
        }
        else if ( ! strcmp( "threePointFivePN", optarg ) )
        {
          order = threePointFivePN;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown order specified: "
              "%s (must be one of: newtonian, oneHalfPN, onePN,\n"
              "onePointFivePN, twoPN, twoPointFivePN, threePN or\n"
              "threePointFivePN)\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        haveOrder = 1;
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'F':
        if ( ! strcmp( "TaylorT1", optarg ) )
        {
          approximant = TaylorT1;
        }
        else if ( ! strcmp( "TaylorT2", optarg ) )
        {
          approximant = TaylorT2;
        }
        else if ( ! strcmp( "TaylorT3", optarg ) )
        {
          approximant = TaylorT3;
        }
        else if ( ! strcmp( "TaylorF1", optarg ) )
        {
          approximant = TaylorF1;
        }
        else if ( ! strcmp( "TaylorF2", optarg ) )
        {
          approximant = TaylorF2;
        }
        else if ( ! strcmp( "PadeT1", optarg ) )
        {
          approximant = PadeT1;
        }
        else if ( ! strcmp( "PadeF1", optarg ) )
        {
          approximant = PadeF1;
        }
        else if ( ! strcmp( "EOB", optarg ) )
        {
          approximant = EOB;
        }
        else if ( ! strcmp( "BCV", optarg ) )
        {
          approximant = BCV;
        }
        else if ( ! strcmp( "SpinTaylorT3", optarg ) )
        {
          approximant = SpinTaylorT3;
        }
        else if ( ! strcmp( "BCVSpin", optarg ) )
        {
          approximant = BCVSpin;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown order specified: "
              "%s (must be one of: TaylorT1, TaylorT2, TaylorT3, TaylorF1,\n"
              "TaylorF2, PadeT1, PadeF1, EOB, BCV, SpinTaylorT3, or BCVSpin)\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        haveApprox = 1;
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'G':
        if ( ! strcmp( "Tau0Tau2", optarg ) )
        {
          space = Tau0Tau2;
        }
        else if ( ! strcmp( "Tau0Tau3", optarg ) )
        {
          space = Tau0Tau3;
        }
        else if ( ! strcmp( "Psi0Psi3", optarg ) )
        {
          space = Psi0Psi3;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown space specified: "
              "%s (must be one of: Tau0Tau2, Tau0Tau3 or Psi0Psi3)\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        haveSpace = 1;
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'v':
        if ( ! strcmp( "Hexagonal", optarg) )
        {
          gridSpacing = Hexagonal;
        }
        else if ( ! strcmp( "SquareNotOriented", optarg) )
        {
          gridSpacing = SquareNotOriented;
        }
        else
        {
          fprintf(stderr, "invalid argument to --%s:\n"
              "unknown grid spacing specified: "
              "%s (must be one of  Hexagonal, SquareNotOriented )\n", 
              long_options[option_index].name, optarg );
          exit(1);
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'y':
        isMaxTotMass = 1;
        ADD_PROCESS_PARAM( "int", "%d", 1 );
        break;
         
      case 'k':
        candleSnr = (REAL4) atof( optarg );
        if ( candleSnr <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "standard candle signal-to-noise ratio must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, candleSnr );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", candleSnr );
        break;

      case 'l':
        candleMass1 = (REAL4) atof( optarg );
        if ( candleMass1 <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "standard candle first component mass must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, candleMass1 );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", candleMass1 );
        break;

      case 'm':
        candleMass2 = (REAL4) atof( optarg );
        if ( candleMass2 <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "standard candle second component mass must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, candleMass2 );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", candleMass2 );
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "LIGO/LSC Standalone Inspiral Template Bank Code\n" 
            "Duncan Brown <duncan@gravity.phys.uwm.edu>\n"
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n" );
        exit( 0 );
        break;

      case '?':
        fprintf( stderr, USAGE );
        exit( 1 );
        break;

      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        fprintf( stderr, USAGE );
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

  /*
   *
   * check validity of arguments
   *
   */


  /* check validity of input data time */
  if ( ! gpsStartTime.gpsSeconds )
  {
    fprintf( stderr, "--gps-start-time must be specified\n" );
    exit( 1 );
  }
  if ( ! gpsEndTime.gpsSeconds )
  {
    fprintf( stderr, "--gps-end-time must be specified\n" );
    exit( 1 );
  }
  if ( gpsEndTime.gpsSeconds <= gpsStartTime.gpsSeconds )
  {
    fprintf( stderr, "invalid gps time range: "
        "start time: %d, end time %d\n",
        gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds );
    exit( 1 );
  }

  /* check validity of data length parameters */
  if ( numPoints < 0 )
  {
    fprintf( stderr, "--segment-length must be specified\n" );
    exit( 1 );
  }
  if ( numSegments < 0 )
  {
    fprintf( stderr, "--number-of-segments must be specified\n" );
    exit( 1 );
  }

  /* check sample rate has been given */
  if ( sampleRate < 0 )
  {
    fprintf( stderr, "--sample-rate must be specified\n" );
    exit( 1 );
  }

  /* check high pass option has been given */
  if ( highPass < 0 )
  {
    fprintf( stderr, "--disable-high-pass or --enable-high-pass (freq)"
        " must be specified\n" );
    exit( 1 );
  }
  else if ( ! highPass )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--disable-high-pass" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }
  else
  {
    /* check that all the high pass parameters have been specified */
    if ( highPassOrder < 0 )
    {
      fprintf( stderr, "--high-pass-order must be specified\n" );
      exit( 1 );
    }
    if ( highPassAtten < 0 )
    {
      fprintf( stderr, "--high-pass-attenuation must be specified\n" );
      exit( 1 );
    }
  }

  if ( calData == real_8 )
  {
    /* check that geo high pass parameters have been specified */
    if ( geoHighPassFreq < 0 )
    {
      fprintf( stderr, 
          "--geo-high-pass-freq must be specified for GEO data\n" );
      exit( 1 );
    }
    if ( geoHighPassOrder < 0 )
    {
      fprintf( stderr, 
          "--geo-high-pass-order must be specified for GEO data\n" );
      exit( 1 );
    }
    if ( geoHighPassAtten < 0 )
    {
      fprintf( stderr, 
          "--geo-high-pass-atten must be specified for GEO data\n" );
      exit( 1 );
    }
  }

  /* check validity of input data length */
  inputDataLength = numPoints * numSegments - ( numSegments - 1 ) * 
    (numPoints / 2);
  {
    UINT8 gpsChanIntervalNS = gpsEndTime.gpsSeconds * 1000000000LL - 
      gpsStartTime.gpsSeconds * 1000000000LL;
    UINT8 inputDataLengthNS = (UINT8) inputDataLength * 1000000000LL / 
      (UINT8) sampleRate;

    if ( inputDataLengthNS != gpsChanIntervalNS )
    {
      fprintf( stderr, "length of input data and data chunk do not match\n" );
      fprintf( stderr, "start time: %d, end time %d\n",
          gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds );
      fprintf( stderr, "gps channel time interval: %lld ns\n"
          "computed input data length: %lld ns\n", 
          gpsChanIntervalNS, inputDataLengthNS );
      exit( 1 );
    }
  }

  /* check standard candle arguments */
  if ( computeCandle )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--standard-candle" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );

    if ( candleSnr < 0 )
    {
      fprintf( stderr, 
          "--candle-snr must be specified if --standard-candle is given\n" );
      exit( 1 );
    }
    if ( candleMass1 < 0 )
    {
      fprintf( stderr, 
          "--candle-mass1 must be specified if --standard-candle is given\n" );
      exit( 1 );
    }
    if ( candleMass2 < 0 )
    {
      fprintf( stderr, 
          "--candle-mass2 must be specified if --standard-candle is given\n" );
      exit( 1 );
    }
  }

  /* check that the spectrum generation parameters have been given */
  if ( fLow < 0 )
  {
    fprintf( stderr, "--low-frequency-cutoff must be specified\n" );
    exit( 1 );
  }
  if ( resampFiltType < 0 )
  {
    fprintf( stderr, "--resample-filter must be specified\n" );
    exit( 1 );
  }
  if ( specType < 0 )
  {
    fprintf( stderr, "--spectrum-type must be specified\n" );
    exit( 1 );
  }

  /* check that a channel has been requested and fill the ifo */
  if ( ! fqChanName )
  {
    fprintf( stderr, "--channel-name must be specified\n" );
    exit( 1 );
  }

  /* check that we can correctly obtain the input frame data */
  if ( globFrameData )
  {
    if ( frInCacheName )
    {
      fprintf( stderr, 
          "--frame-cache must not be specified when globbing frame data\n" );
      exit( 1 );
    }

    if ( ! frInType )
    {
      fprintf( stderr, 
          "--frame-type must be specified when globbing frame data\n" );
      exit( 1 );
    }
  }
  else
  {
    if ( ! frInCacheName )
    {
      fprintf( stderr, 
          "--frame-cache must be specified when not globbing frame data\n" );
      exit( 1 );
    }

    if ( frInType )
    {
      fprintf( stderr, "--frame-type must not be specified when obtaining "
          "frame data from a cache file\n" );
      exit( 1 );
    }
  }

  /* record the glob frame data option in the process params */
  if ( globFrameData )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--glob-frame-data" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  /* store point calibration option */
  if ( pointCal )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--point-calibration" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  /* check we can calibrate the data if it's not h(t) */
  if ( ! calData )
  {
    if ( ! ( calCacheName || globCalData ) )
    {
      fprintf( stderr, "either --calibration-cache or "
          "--glob-calibration-data must be specified\n" );
      exit( 1 );
    }
    else if ( calCacheName && globCalData )
    {
      fprintf( stderr, "only one of --calibration-cache or "
          "--glob-calibration-data can be specified\n" );
      exit( 1 );
    }
  }
  else
  {
    if ( calCacheName || globCalData )
    {
      fprintf( stderr, "neither --calibration-cache nor "
          "--glob-calibration-data\nshould be given for calibrated data\n" );
      exit( 1 );
    }
  }

  /* record the glob calibration data option in the process params */
  if ( globCalData )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--glob-calibration-data" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  /* check that the bank type has been specified */
  if ( ! haveOrder )
  {
    fprintf( stderr, "--order must be specified\n" );
    exit( 1 );
  }
  if ( ! haveApprox )
  {
    fprintf( stderr, "--approximant must be specified\n" );
    exit( 1 );
  }
  if ( ! haveSpace )
  {
    fprintf( stderr, "--space must be specified\n" );
    exit( 1 );
  }

  /* check validity of grid spacing with respect to approximant */
  if (gridSpacing != SquareNotOriented && gridSpacing != Hexagonal)
  {
    fprintf( stderr, "--grid-spacing  must be either SquareNotOriented or Hexagonal\n" );
    exit( 1 );
  }
  
  /* check that the correct range parameters have been given for the bank */
  if ( approximant == BCV )
  {
    if ( ! havePsi0Min )
    {
      fprintf( stderr, "--minimum-psi0 must be specified\n" );
      exit( 1 );
    }
    if ( ! havePsi0Max )
    {
      fprintf( stderr, "--maximum-psi0 must be specified\n" );
      exit( 1 );
    }
    if ( ! havePsi3Min )
    {
      fprintf( stderr, "--minimum-psi3 must be specified\n" );
      exit( 1 );
    }
    if ( ! havePsi3Max )
    {
      fprintf( stderr, "--maximum-psi3 must be specified\n" );
      exit( 1 );
    }
    if ( ! haveAlpha )
    {
      fprintf( stderr, "--alpha must be specified\n" );
      exit( 1 );
    }
    if ( maxFcutTmplts < 0 )
    {
      fprintf( stderr, "--maximum-fcut-tmplts must be specified\n" );
      exit( 1 );
    }

    if ( psi3Max <= psi3Min )
    {
      fprintf( stdout, "invalid argument to --maximum-psi3:\n"
          "maximum value of psi3 must be greater than minimum value of psi3: "
          "(%f specified)\n",
          psi3Max );
      exit( 1 );
    }

    minMass = maxMass = 0;
  }
  else
  {
    if ( minMass < 0 )
    {
      fprintf( stderr, "--minimum-mass must be specified\n" );
      exit( 1 );
    }
    if ( maxMass < 0 )
    {
      fprintf( stderr, "--maximum-mass must be specified\n" );
      exit( 1 );
    }
  }

  /* check that the bank parameters have been specified */
  if ( minMatch < 0 )
  {
    fprintf( stderr, "--minimal-match must be specified\n" );
    exit( 1 );
  }
  if ( fUpper < 0 )
  {
    fprintf( stderr, "--high-frequency-cutoff must be specified\n" );
    exit( 1 );
  }

  return 0;
}

#undef ADD_PROCESS_PARAM

