/*----------------------------------------------------------------------- 
 * 
 * File Name: inspiral.c
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
#include <math.h>

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
#include <lal/FrameStream.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/Calibration.h>
#include <lal/FrameCalibration.h>
#include <lal/Window.h>
#include <lal/TimeFreqFFT.h>
#include <lal/IIRFilter.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpTD.h>
#include <lal/FindChirpBCV.h>
#include <lal/FindChirpBCVSpin.h>
#include <lal/FindChirpChisq.h>

RCSID( "$Id$" );

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "inspiral"

/* define the parameters for a 1.4,1.4 sloar mass standard candle with snr 8 */
#define CANDLE_MASS1 1.4
#define CANDLE_MASS2 1.4
#define CANDLE_RHOSQ 64.0


#define ADD_SUMM_VALUE( sv_name, sv_comment, val, intval ) \
if ( this_summ_value ) \
{ \
  this_summ_value = this_summ_value->next = (SummValueTable *) \
    LALCalloc( 1, sizeof(SummValueTable) ); \
} \
else \
{ \
  summvalue.summValueTable = this_summ_value = (SummValueTable *) \
    LALCalloc( 1, sizeof(SummValueTable) ); \
} \
this_summ_value->version = 0; \
this_summ_value->start_time = searchsumm.searchSummaryTable->out_start_time; \
this_summ_value->end_time = searchsumm.searchSummaryTable->out_end_time; \
this_summ_value->value = (REAL4) val; \
this_summ_value->intvalue = (INT4) intval; \
LALSnprintf( this_summ_value->name, LIGOMETA_SUMMVALUE_NAME_MAX, "%s", \
    sv_name ); \
LALSnprintf( this_summ_value->ifo, LIGOMETA_IFO_MAX, "%s", ifo ); \
LALSnprintf( this_summ_value->comment, LIGOMETA_SUMMVALUE_COMM_MAX, \
    "%s", sv_comment ); \

double rint(double x);
int arg_parse_check( int argc, char *argv[], MetadataTable procparams );

#ifdef LALAPPS_CONDOR
extern int condor_compress_ckpt;
void init_image_with_file_name( char *ckpt_file_name );
void ckpt_and_exit( void );
#endif


/*
 *
 * variables that control program behaviour
 *
 */


/* debugging */
extern int vrbflg;                      /* verbocity of lal function    */

/* checkpointing */
INT4  dataCheckpoint = 0;               /* condor checkpoint after data */
CHAR  ckptPath[FILENAME_MAX];           /* input and ckpt file path     */
CHAR  outputPath[FILENAME_MAX];         /* output data file path        */

/* input data parameters */
INT8  gpsStartTimeNS   = 0;             /* input data GPS start time ns */
LIGOTimeGPS gpsStartTime;               /* input data GPS start time    */
INT8  gpsEndTimeNS     = 0;             /* input data GPS end time ns   */
LIGOTimeGPS gpsEndTime;                 /* input data GPS end time      */
INT4  padData = 0;                      /* saftety margin on input data */
CHAR  *fqChanName       = NULL;         /* name of data channel         */
INT4  globFrameData     = 0;            /* glob *.gwf to get frame data */
CHAR  *frInCacheName    = NULL;         /* cache file containing frames */
CHAR  *frInType         = NULL;         /* type of data frames          */
INT4  numPoints         = -1;           /* points in a segment          */
INT4  numSegments       = -1;           /* number of segments           */
INT4  ovrlap            = -1;           /* overlap between segments     */
CHAR  ifo[3];                           /* two character ifo code       */
CHAR *channelName = NULL;               /* channel string               */
UINT4 inputDataLength = 0;              /* number of points in input    */
REAL4 minimalMatch = -1;                /* override bank minimal match  */
REAL4 geoHighPassFreq = -1;             /* GEO high pass frequency      */
INT4  geoHighPassOrder = -1;            /* GEO high pass filter order   */
REAL4 geoHighPassAtten = -1;            /* GEO high pass attenuation    */
enum { undefined, real_4, real_8 } calData = undefined; /* cal data type */

/* data conditioning parameters */
LIGOTimeGPS slideData   = {0,0};        /* slide data for time shifting */
INT4   resampFiltType   = -1;           /* low pass filter used for res */
INT4   sampleRate       = -1;           /* sample rate of filter data   */
INT4   highPass         = -1;           /* enable high pass on raw data */
REAL4  highPassFreq     = 0;            /* high pass frequency          */
INT4   highPassOrder    = -1;           /* order of the td iir filter   */
REAL4  highPassAtten    = -1;           /* attenuation of the td filter */

REAL4  fLow             = -1;           /* low frequency cutoff         */
INT4   specType         = -1;           /* use median or mean psd       */
INT4   badMeanPsd       = 0;            /* use a mean with no overlap   */
INT4   invSpecTrunc     = -1;           /* length of inverse spec (s)   */
REAL4  dynRangeExponent = -1;           /* exponent of dynamic range    */
CHAR  *calCacheName     = NULL;         /* location of calibration data */
INT4   globCalData      = 0;            /* glob for calibration frames  */
INT4   pointCal         = 0;            /* don't average cal over chunk */
CHAR  *injectionFile    = NULL;         /* name of file containing injs */
int   injectOverhead    = 0;            /* inject h+ into detector      */

/* matched filter parameters */
CHAR *bankFileName      = NULL;         /* name of input template bank  */
INT4  startTemplate     = -1;           /* index of first template      */
INT4  stopTemplate      = -1;           /* index of last template       */
INT4  numChisqBins      = -1;           /* number of chisq bins         */
REAL4 snrThresh         = -1;           /* signal to noise thresholds   */
REAL4 chisqThresh       = -1;           /* chisq veto thresholds        */
Clustering clusterMethod;               /* chosen clustering algorithm  */  
REAL4 clusterWindow     = -1;           /* cluster over time window     */  
Approximant approximant;                /* waveform approximant         */
INT4 bcvConstraint      = 0;            /* constraint BCV filter        */

/* generic simulation parameters */
enum { unset, urandom, user } randSeedType = unset;    /* sim seed type */
INT4  randomSeed        = 0;            /* value of sim rand seed       */
REAL4 gaussVar          = 64.0;         /* variance of gaussian noise   */
INT4  gaussianNoise     = 0;            /* make input data gaussian     */

/* template bank simulation params */
INT4  bankSim           = 0;            /* number of template bank sims */
Approximant bankSimApproximant;         /* waveform to test tmplt bank  */
REAL4 bankMinMass       = -1;           /* minimum mass of injection    */
REAL4 bankMaxMass       = -1;           /* maximum mass of injection    */

/* output parameters */
CHAR  *userTag          = NULL;         /* string the user can tag with */
CHAR  *ifoTag           = NULL;         /* string to tag parent IFOs    */
CHAR   fileName[FILENAME_MAX];          /* name of output files         */
INT8   trigStartTimeNS  = 0;            /* write triggers only after    */
INT8   trigEndTimeNS    = 0;            /* write triggers only before   */
INT8   outTimeNS        = 0;            /* search summ out time         */
int    enableOutput     = -1;           /* write out inspiral events    */
int    writeRawData     = 0;            /* write the raw data to a file */
int    writeFilterData  = 0;            /* write post injection data    */
int    writeResponse    = 0;            /* write response function used */
int    writeSpectrum    = 0;            /* write computed psd to file   */
int    writeRhosq       = 0;            /* write rhosq time series      */
int    writeChisq       = 0;            /* write chisq time series      */
int    writeCData       = 0;            /* write complex time series c  */

/* other command line args */
CHAR comment[LIGOMETA_COMMENT_MAX];     /* process param comment        */

int main( int argc, char *argv[] )
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
  struct FrFile *frOutFile  = NULL;
  struct FrameH *outFrame   = NULL;
  UINT4          nRhosqFr = 0;
  UINT4          nChisqFr = 0;
  UINT4          nCDataFr = 0;

  /* raw input data storage */
  REAL4TimeSeries               chan;
  REAL8TimeSeries               geoChan;
  REAL4FrequencySeries          spec;
  COMPLEX8FrequencySeries       resp;
  DataSegmentVector            *dataSegVec = NULL;
  COMPLEX8TimeSeries           *coherentInputData = NULL;

  /* structures for preconditioning */
  COMPLEX8FrequencySeries       injResp;
  COMPLEX8FrequencySeries      *injRespPtr;
  ResampleTSParams              resampleParams;
  LALWindowParams               wpars;
  AverageSpectrumParams         avgSpecParams;

  /* findchirp data structures */
  FindChirpInitParams          *fcInitParams   = NULL;
  FindChirpSegmentVector       *fcSegVec       = NULL;
  FindChirpDataParams          *fcDataParams   = NULL;
  FindChirpTmpltParams         *fcTmpltParams  = NULL;
  FindChirpFilterParams        *fcFilterParams = NULL;
  FindChirpFilterInput         *fcFilterInput  = NULL;
  FindChirpStandardCandle       candle;

  /* inspiral template structures */
  INT4                          numTmplts    = 0;
  InspiralTemplate             *bankHead     = NULL;
  InspiralTemplate             *bankCurrent  = NULL;
  InspiralTemplateNode         *tmpltHead    = NULL;
  InspiralTemplateNode         *tmpltCurrent = NULL;
  InspiralTemplateNode         *tmpltInsert  = NULL;

  /* inspiral events */
  INT4                          numEvents   = 0;
  SnglInspiralTable            *event       = NULL;
  SnglInspiralTable            *eventList   = NULL;
  SnglInspiralTable            *tempTmplt   = NULL;
  MetadataTable                 savedEvents;

  /* output data */
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         searchsumm;
  MetadataTable         searchsummvars;
  SearchSummvarsTable  *this_search_summvar;
  MetadataTable         summvalue;
  SummValueTable       *this_summ_value = NULL;
  ProcessParamsTable   *this_proc_param;
  LIGOLwXMLStream       results;

  /* counters and other variables */
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
  UINT4 i, j, k;
  INT4  inserted;
  INT4  currentLevel;
  INT4  cDataForFrame;
  CHAR  fname[FILENAME_MAX];
  CHAR  cdataStr[LALNameLength];
  REAL8 inputLengthNS;
  UINT4 numInputPoints;
  const REAL8 epsilon = 1.0e-8;
  UINT4 resampleChan = 0;
  REAL8 tsLength;
  INT8  durationNS = 0;
  CalibrationUpdateParams calfacts, inj_calfacts;
  REAL4 alpha = 0;
  REAL4 alphabeta = 0;
  REAL4 inj_alpha = 0;
  REAL4 inj_alphabeta = 0;
  CHAR tmpChName[LALNameLength];
  REAL8 inputDeltaT;
  REAL8 dynRange = 1.0;
  REAL8 trigTime = 0.0;
  REAL8 lowerBound = 0.0;
  REAL8 upperBound = 0.0;

  /* random number generator parameters */
  RandomParams *randParams = NULL;

  /* template bank simulation variables */
  UINT4 cut = 0;
  REAL4 psdMin = 0;
  INT4  bankSimCount = 0;
  SimInspiralTable bankInjection;
  REAL4 matchNorm = 0;
  const REAL8 psdScaleFac = 1.0e-40;
  SnglInspiralTable *loudestEvent = NULL;
  SnglInspiralTable *prevLoudestEvent = NULL;
  SnglInspiralTable *loudestEventHead = NULL;
  SimInstParamsTable *thisSimInstParams = NULL;
  SimInstParamsTable *prevSimInstParams = NULL;
  MetadataTable simResults;


  /*
   *
   * initialization
   *
   */


  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
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

  /* zero out the checkpoint and output paths */
  memset( ckptPath, 0, FILENAME_MAX * sizeof(CHAR) );
  memset( outputPath, 0, FILENAME_MAX * sizeof(CHAR) );

  /* call the argument parse and check function */
  arg_parse_check( argc, argv, procparams );

  /* wind to the end of the process params table */
  for ( this_proc_param = procparams.processParamsTable; this_proc_param->next;
      this_proc_param = this_proc_param->next );

  /* can use LALMalloc() and LALCalloc() from here onwards */

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

  /* set the name of the output file */
  if ( userTag && ifoTag )
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-INSPIRAL_%s_%s-%d-%d", ifo, 
        ifoTag, userTag, gpsStartTime.gpsSeconds, 
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds ); 
  }
  else if ( userTag && !ifoTag )
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-INSPIRAL_%s-%d-%d", ifo, 
        userTag,  gpsStartTime.gpsSeconds, 
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else if ( !userTag && ifoTag )
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-INSPIRAL_%s-%d-%d", ifo, 
        ifoTag,  gpsStartTime.gpsSeconds, 
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds ); 
  }
  else
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-INSPIRAL-%d-%d", ifo,
        gpsStartTime.gpsSeconds, 
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );     
  }

  /* the number of nodes for a standalone job is always 1 */
  searchsumm.searchSummaryTable->nnodes = 1;

  /* fill the ifos field of the search summary table */
  LALSnprintf( searchsumm.searchSummaryTable->ifos, LIGOMETA_IFOS_MAX, ifo );

  /* make sure the pointer to the first event is null */
  savedEvents.snglInspiralTable = NULL;

  /* create the standard candle and database table */
  memset( &candle, 0, sizeof(FindChirpStandardCandle) );
  strncpy( candle.ifo, ifo, 2 * sizeof(CHAR) );
  candle.tmplt.mass1 = CANDLE_MASS1;
  candle.tmplt.mass2 = CANDLE_MASS2;
  candle.rhosq       = CANDLE_RHOSQ;
  candle.tmplt.totalMass = candle.tmplt.mass1 + candle.tmplt.mass2;
  candle.tmplt.mu = candle.tmplt.mass1 * candle.tmplt.mass2 / 
    candle.tmplt.totalMass;
  candle.tmplt.eta = candle.tmplt.mu / candle.tmplt.totalMass;

  /* create the dynamic range exponent */
  dynRange = pow( 2.0, dynRangeExponent );


  /*
   *
   * create a (possibly heirarcical) template bank
   *
   */


  /* read in the template bank from a ligo lw xml file */
  numTmplts = InspiralTmpltBankFromLIGOLw( &bankHead, bankFileName,
      startTemplate, stopTemplate );
  if ( numTmplts < 0 )
  {
    fprintf( stderr, "error: unable to read templates from %s\n", 
        bankFileName );
    exit( 1 );
  }
  else if ( numTmplts == 0 )
  {
    /* if there are no tmplts, store the time we would have analyzed and exit */
    fprintf( stdout, "no templates found in template bank file: %s\n"
        "exiting without searching for events.\n", bankFileName );

    searchsumm.searchSummaryTable->out_start_time.gpsSeconds = 
        gpsStartTime.gpsSeconds + (numPoints / (4 * sampleRate));

    LAL_CALL( LALGPStoINT8( &status, &outTimeNS,
            &(searchsumm.searchSummaryTable->out_start_time) ), &status );

    if ( trigStartTimeNS && (trigStartTimeNS > outTimeNS) )
    {
      LAL_CALL( LALINT8toGPS( &status, 
            &(searchsumm.searchSummaryTable->out_start_time), 
            &trigStartTimeNS ), &status );
    }
    
    searchsumm.searchSummaryTable->out_end_time.gpsSeconds = 
        gpsEndTime.gpsSeconds - (numPoints / (4 * sampleRate));

    LAL_CALL( LALGPStoINT8( &status, &outTimeNS,
            &(searchsumm.searchSummaryTable->out_end_time) ), &status );

    if ( trigEndTimeNS && (trigEndTimeNS < outTimeNS) )
    {
      LAL_CALL( LALINT8toGPS( &status, 
            &(searchsumm.searchSummaryTable->out_end_time), 
            &trigEndTimeNS ), &status );
    }
  }

  if ( vrbflg ) fprintf( stdout, "parsed %d templates from %s\n", 
      numTmplts, bankFileName );

  /* override the minimal match of the bank if specified on the command line */
  if ( minimalMatch >= 0 )
  {
    if ( vrbflg && bankHead )
    {
      fprintf( stdout, "Overriding bank minimal match:\n   value in bank = %e,"
          " new value = %e\n", bankHead->minMatch, minimalMatch );
    }
    for ( bankCurrent = bankHead; bankCurrent; bankCurrent = bankCurrent->next )
    {
      bankCurrent->minMatch = minimalMatch;
    }
  }

  if ( numTmplts )
  {
    /* save the minimal match of the bank in the process params: create */
    /* the table entry with calloc() since it will be freed with free() */
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) 
      calloc( 1, sizeof(ProcessParamsTable) ); 
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
        PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--minimal-match" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "float" ); 
    LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%e", 
        bankHead->minMatch );
  }

  /* create the linked list of template nodes for coarse templates */
  for ( bankCurrent = bankHead; bankCurrent; bankCurrent = bankCurrent->next )
  {
    LAL_CALL( LALFindChirpCreateTmpltNode( &status, 
          bankCurrent, &tmpltCurrent ), &status );
    if ( !tmpltHead ) tmpltHead = tmpltCurrent;
  }


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
  chan.epoch.gpsSeconds -= padData;   /* subtract pad seconds from start */
  /* subtract slide from start */
  chan.epoch.gpsSeconds -= slideData.gpsSeconds; 
  chan.epoch.gpsNanoSeconds -= slideData.gpsNanoSeconds;
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
    /* determine the sample rate of the raw data */
    LAL_CALL( LALFrGetREAL4TimeSeries( &status, &chan, &frChan, frStream ),
        &status );
  }

  /* store the input sample rate */
  this_search_summvar = searchsummvars.searchSummvarsTable = 
    (SearchSummvarsTable *) LALCalloc( 1, sizeof(SearchSummvarsTable) );
  LALSnprintf( this_search_summvar->name, LIGOMETA_NAME_MAX * sizeof(CHAR),
      "raw data sample rate" );
  this_search_summvar->value = inputDeltaT = chan.deltaT;

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

  /* determine the number of points to get and create storage for the data */
  inputLengthNS = 
    (REAL8) ( gpsEndTimeNS - gpsStartTimeNS + 2000000000LL * padData );
  numInputPoints = (UINT4) floor( inputLengthNS / (chan.deltaT * 1.0e9) + 0.5 );
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
      chan.deltaT, numInputPoints );

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


  /*
   *
   * create the random seed if it will be needed
   *
   */


  if ( randSeedType != unset )
  {
    /* store the seed in the search summvars table */
    this_search_summvar = this_search_summvar->next = 
      (SearchSummvarsTable *) LALCalloc( 1, sizeof(SearchSummvarsTable) );
    LALSnprintf( this_search_summvar->name, LIGOMETA_NAME_MAX * sizeof(CHAR),
        "template bank simulation seed" );

    if ( randSeedType == urandom )
    {
      FILE   *fpRand = NULL;
      INT4    randByte;

      if ( vrbflg ) 
        fprintf( stdout, "obtaining random seed from /dev/urandom: " );

      randomSeed = 0;
      fpRand = fopen( "/dev/urandom", "r" );
      if ( fpRand )
      {
        for ( randByte = 0; randByte < 4 ; ++randByte )
        {
          INT4 tmpSeed = (INT4) fgetc( fpRand );
          randomSeed += tmpSeed << ( randByte * 8 );
        }
        fclose( fpRand );
      }
      else
      {
        perror( "error obtaining random seed from /dev/urandom" );
        exit( 1 );
      }
    }
    else if ( randSeedType == user )
    {
      if ( vrbflg ) 
        fprintf( stdout, "using user specified random seed: " );
    }
    else
    {
      /* should never get here */
      fprintf( stderr, "error obtaining random seed\n" );
      exit( 1 );
    }

    this_search_summvar->value = randomSeed;
    if ( vrbflg ) fprintf( stdout, "%d\n", randomSeed );

    /* create the tmplt bank random parameter structure */
    LAL_CALL( LALCreateRandomParams( &status, &randParams, randomSeed ),
        &status );
  }

  /* replace the input data with gaussian noise if necessary */
  if ( gaussianNoise )
  {
    if ( vrbflg ) fprintf( stdout, 
        "setting input data to gaussian noise with variance %e... ", gaussVar );
    memset( chan.data->data, 0, chan.data->length * sizeof(REAL4) );
    LAL_CALL( LALNormalDeviates( &status, chan.data, randParams ), &status );
    for ( j = 0; j < chan.data->length; ++j )
    {
      chan.data->data[j] *= gaussVar;
    }
    if ( vrbflg ) fprintf( stdout, "done\n" );

    /* write the raw channel data as read in from the frame files */
    if ( writeRawData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
        &chan, "ct", "RAW_GAUSSIAN" );
  }


  /*
   *
   * generate the response function for the requested time
   *
   */


  /* create storage for the response function */
  memset( &resp, 0, sizeof(COMPLEX8FrequencySeries) );
  LAL_CALL( LALCCreateVector( &status, &(resp.data), numPoints / 2 + 1 ), 
      &status );

  /* set the parameters of the response to match the data */
  resp.epoch.gpsSeconds = chan.epoch.gpsSeconds + padData;
  resp.epoch.gpsNanoSeconds = chan.epoch.gpsNanoSeconds;
  resp.deltaF = (REAL8) sampleRate / (REAL8) numPoints;
  resp.sampleUnits = strainPerCount;
  strcpy( resp.name, chan.name );

  /* generate the response function for the current time */
  if ( vrbflg ) fprintf( stdout, "generating response at time %d sec %d ns\n",
      resp.epoch.gpsSeconds, resp.epoch.gpsNanoSeconds );

  /* initialize the calfacts */
  memset( &calfacts, 0, sizeof(CalibrationUpdateParams) );
  calfacts.ifo = ifo;

  /* determine length of chunk */
  if ( pointCal )
  {
    calfacts.duration.gpsSeconds = 1;
    calfacts.duration.gpsNanoSeconds = 0;
  }
  else
  {
    durationNS = gpsEndTimeNS - gpsStartTimeNS;
    LAL_CALL( LALINT8toGPS( &status, &(calfacts.duration), 
          &durationNS ), &status );
  }

  if ( calData )
  {
    /* if we are using calibrated data set the response to unity */
    for( k = 0; k < resp.data->length; ++k )
    {
      resp.data->data[k].re = (REAL4) (1.0 / dynRange);
      resp.data->data[k].im = 0.0;
    }
    if ( writeResponse ) outFrame = fr_add_proc_COMPLEX8FrequencySeries( 
        outFrame, &resp, "strain/ct", "RESPONSE_GEO" );
  }
  else
  {
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
    
    /* get the response from the frame data */
    LAL_CALL( LALExtractFrameResponse( &status, &resp, calCache, 
          &calfacts), &status );
    LAL_CALL( LALDestroyFrCache( &status, &calCache), &status );
    alpha = (REAL4) calfacts.alpha.re;
    alphabeta = (REAL4) calfacts.alphabeta.re;
    if ( vrbflg ) fprintf( stdout, 
        "for calibration of data, alpha = %f and alphabeta = %f\n",
        alpha, alphabeta);

    if ( writeResponse ) 
      outFrame = fr_add_proc_COMPLEX8FrequencySeries( outFrame, &resp, 
          "strain/ct", "RESPONSE" );
  }

  if ( gaussianNoise )
  {
    /* replace the response function with unity if */
    /* we are filtering gaussian noise             */
    if ( vrbflg ) fprintf( stdout, "setting response to unity... " );
    for ( k = 0; k < resp.data->length; ++k )
    {
      resp.data->data[k].re = 1.0;
      resp.data->data[k].im = 0;
    }
    if ( vrbflg ) fprintf( stdout, "done\n" );

    if ( writeResponse ) outFrame = fr_add_proc_COMPLEX8FrequencySeries( 
        outFrame, &resp, "strain/ct", "RESPONSE_GAUSSIAN" );
  }

  /* slide the channel back to the fake time for background studies */
  chan.epoch.gpsSeconds += slideData.gpsSeconds;
  chan.epoch.gpsNanoSeconds += slideData.gpsNanoSeconds;


  /*
   *
   * inject signals into the raw, unresampled data
   *
   */


  if ( injectionFile )
  {
    /* get injections within 500 seconds of either end of the segment.   */
    /* a 0.4,0.4 MACHO starting at 30.0 Hz has length 435.374683 seconds */
    /* so this should be plenty of safety. better to waste cpu than miss */
    /* injected signals...                                               */
    INT4 injSafety = 500;
    int  numInjections = 0;
    SimInspiralTable    *injections = NULL;
    SimInspiralTable    *thisInj = NULL;

    /* read in the injection data from XML */
    numInjections = SimInspiralTableFromLIGOLw( &injections, injectionFile,
        gpsStartTime.gpsSeconds - injSafety, 
        gpsEndTime.gpsSeconds + injSafety );

    if ( numInjections < 0 )
    {
      fprintf( stderr, "error: cannot read injection file" );
      exit( 1 );
    }
    else if ( numInjections )
    {
      /* see if we need a higher resolution response to do the injections */
      if ( resampleChan )
      {
        /* we need a different resolution of response function for injections */
        UINT4 rateRatio = floor( resampleParams.deltaT / chan.deltaT + 0.5 );
        UINT4 rawNumPoints = rateRatio * numPoints;

        if ( vrbflg ) fprintf( stdout, "rateRatio = %d\nrawNumPoints = %d\n"
            "chan.deltaT = %e\n", rateRatio, rawNumPoints, chan.deltaT );

        memset( &injResp, 0, sizeof(COMPLEX8FrequencySeries) );
        LAL_CALL( LALCCreateVector( &status, &(injResp.data), 
              rawNumPoints / 2 + 1 ), &status );
        injResp.epoch = resp.epoch;
        injResp.deltaF = 1.0 / ( rawNumPoints * chan.deltaT );
        injResp.sampleUnits = strainPerCount;
        strcpy( injResp.name, chan.name );

        if ( calData )
        {
          /* if we are using calibrated data set the response to unity */
          if ( vrbflg ) fprintf( stdout, 
              "setting injection response to inverse dynRange... " );
          for ( k = 0; k < injResp.data->length; ++k )
          {
            injResp.data->data[k].re = (REAL4)(1.0/dynRange);
            injResp.data->data[k].im = 0.0;
          }
          injRespPtr = &injResp;
          if ( writeResponse ) 
            outFrame = fr_add_proc_COMPLEX8FrequencySeries( outFrame, 
                &injResp, "strain/ct", "RESPONSE_INJ_CAL" );
        }
        else
        {
          /* generate the response function for the current time */
          if ( vrbflg ) fprintf( stdout, 
              "generating high resolution response at time %d sec %d ns\n"
              "length = %d points, deltaF = %e Hz\n",
              resp.epoch.gpsSeconds, resp.epoch.gpsNanoSeconds,
              injResp.data->length, injResp.deltaF );

          /* initialize the inj_calfacts */
          memset( &inj_calfacts, 0, sizeof(CalibrationUpdateParams) );
          inj_calfacts.ifo = ifo;
          LAL_CALL( LALINT8toGPS( &status, &(inj_calfacts.duration), 
                &durationNS ), &status );

          /* create the lal calibration frame cache */
          if ( globCalData )
          {
            calGlobPattern = (CHAR *) LALCalloc( calGlobLen, sizeof(CHAR) );
            LALSnprintf( calGlobPattern, calGlobLen * sizeof(CHAR), 
                "*%c*CAL*.gwf", fqChanName[0] );
            if ( vrbflg ) fprintf( stdout, 
                "globbing for %s calibration frame files "
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
              (SearchSummvarsTable *) 
              LALCalloc( 1, sizeof(SearchSummvarsTable) );
            LALSnprintf( this_search_summvar->name, 
                LIGOMETA_NAME_MAX * sizeof(CHAR), 
                "injection calibration frame %d", i );
            LALSnprintf( this_search_summvar->string, 
                LIGOMETA_STRING_MAX * sizeof(CHAR), "%s", 
                calCache->frameFiles[i].url );
          }
          
          /* extract the calibration from frames */
          LAL_CALL( LALExtractFrameResponse( &status, &injResp, calCache, 
                &inj_calfacts ), &status );
          LAL_CALL( LALDestroyFrCache( &status, &calCache), &status );

          inj_alpha = (REAL4) calfacts.alpha.re;
          inj_alphabeta = (REAL4) calfacts.alphabeta.re;
          if ( vrbflg ) fprintf( stdout, 
              "for injections, alpha = %f and alphabeta = %f\n",
              inj_alpha, inj_alphabeta);

          injRespPtr = &injResp;

          if ( writeResponse ) 
            outFrame = fr_add_proc_COMPLEX8FrequencySeries( 
                outFrame, &injResp, "strain/ct", "RESPONSE_INJ" );
        }

        if ( gaussianNoise )
        {
          /* replace the response function with unity if */
          /* we are filtering gaussian noise             */
          if ( vrbflg ) fprintf( stdout, "setting response to unity... " );
          for ( k = 0; k < injResp.data->length; ++k )
          {
            injResp.data->data[k].re = 1.0;
            injResp.data->data[k].im = 0;
          }
          if ( vrbflg ) fprintf( stdout, "done\n" );

          if ( writeResponse ) outFrame = fr_add_proc_COMPLEX8FrequencySeries( 
              outFrame, &injResp, "strain/ct", "RESPONSE_INJ_GAUSSIAN" );
        }
      }
      else
      {
        /* the data is already at the correct sample rate, just do injections */
        injRespPtr = &resp;
        memset( &injResp, 0, sizeof(COMPLEX8FrequencySeries) );
      }

      /* inject the signals, preserving the channel name (Tev mangles it) */
      LALSnprintf( tmpChName, LALNameLength * sizeof(CHAR), "%s", chan.name );

      /* if injectOverhead option, then set chan.name to "ZENITH".  
       * This causes no detector site to be found in the injection code so
       * that the injection is done directly overhead (i.e. with a response 
       * function of F+ = 1; Fx = 0) */
      if ( injectOverhead )
      {
        LALSnprintf( chan.name, LALNameLength * sizeof(CHAR), "ZENITH" );
      }

      LAL_CALL( LALFindChirpInjectSignals( &status, &chan, injections, 
            injRespPtr ), &status );
      LALSnprintf( chan.name,  LALNameLength * sizeof(CHAR), "%s", tmpChName );

      if ( vrbflg ) fprintf( stdout, "injected %d signals from %s into %s\n", 
          numInjections, injectionFile, chan.name );

      while ( injections )
      {
        thisInj = injections;
        injections = injections->next;
        LALFree( thisInj );
      }

      /* write the raw channel data plus injections to the output frame file */
      if ( writeRawData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
          &chan, "ct", "RAW_INJ" );

      if ( injResp.data )
        LAL_CALL( LALCDestroyVector( &status, &(injResp.data) ), &status );
    }
    else
    {
      if ( vrbflg ) fprintf( stdout, "no injections in this chunk\n" );
    }
  }


  /*
   *
   * resample the data to the requested rate
   *
   */


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
   * we have all the input data that we need, so checkpoint if requested
   *
   */


  if ( dataCheckpoint )
  {
#ifdef LALAPPS_CONDOR
    condor_compress_ckpt = 1;
    if ( ckptPath[0] )
    {
      LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), "%s/%s.ckpt", 
          ckptPath, fileName );
    }
    else
    {
      LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), "%s.ckpt", fileName );
    }
    if ( vrbflg ) fprintf( stdout, "checkpointing to file %s\n", fname );
    init_image_with_file_name( fname );
    ckpt_and_exit();
#else
    fprintf( stderr, "--data-checkpoint cannot be used unless "
        "lalapps is condor compiled\n" );
    exit( 1 );
#endif
  }


  /* 
   *
   * high pass the data, removed pad from time series and check length of data
   *
   */


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

    LAL_CALL( LALButterworthREAL4TimeSeries( &status, &chan, &highpassParam ),
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

  if ( writeFilterData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
      &chan, "ct", "FILTER" );

  /* check data length */
  if ( chan.data->length != inputDataLength )
  {
    fprintf( stderr, "error: computed channel length and requested\n"
        "input data length do not match:\nchan.data->length = %d\n"
        "inputDataLength = %d\nyou have found a bug in the code.\n"
        "please report this to <duncan@gravity.phys.uwm.edu>\n",
        chan.data->length, inputDataLength );
    exit( 1 );
  }

  /* store the start and end time of the filter channel in the search summ */
  /* noting that we don't look for events in the first and last quarter    */
  /* of each findchirp segment of the input data                           */
  LAL_CALL( LALGPStoFloat( &status, &tsLength, &(chan.epoch) ), 
      &status );
  tsLength += (REAL8) (numPoints / 4) * chan.deltaT;
  LAL_CALL( LALFloatToGPS( &status, 
      &(searchsumm.searchSummaryTable->out_start_time), &tsLength ), 
      &status );


  LAL_CALL( LALGPStoINT8( &status, &outTimeNS,
      &(searchsumm.searchSummaryTable->out_start_time) ), &status );

  if ( trigStartTimeNS && (trigStartTimeNS > outTimeNS) )
  {
    /* override with trigger start time */
    LAL_CALL( LALINT8toGPS( &status,
        &(searchsumm.searchSummaryTable->out_start_time), 
        &trigStartTimeNS ), &status );
  }

  LAL_CALL( LALGPStoFloat( &status, &tsLength, &(chan.epoch) ), 
      &status );
  tsLength += chan.deltaT * ((REAL8) chan.data->length - (REAL8) (numPoints/4));
  LAL_CALL( LALFloatToGPS( &status, 
      &(searchsumm.searchSummaryTable->out_end_time), &tsLength ), 
      &status );

  
  LAL_CALL( LALGPStoINT8( &status, &outTimeNS,
        &(searchsumm.searchSummaryTable->out_end_time) ), &status );

  if ( trigEndTimeNS && (trigEndTimeNS < outTimeNS) )
  {
    /* override with trigger end time */
    LAL_CALL( LALINT8toGPS( &status, 
        &(searchsumm.searchSummaryTable->out_end_time), 
        &trigEndTimeNS ), &status );
  }

  /* 
   *
   * create and populate findchip initialization structure 
   *
   */


  if ( ! ( fcInitParams = (FindChirpInitParams *) 
        LALCalloc( 1, sizeof(FindChirpInitParams) ) ) )
  {
    fprintf( stderr, "could not allocate memory for findchirp init params\n" );
    exit( 1 );
  }

  fcInitParams->numPoints      = numPoints;
  fcInitParams->numSegments    = numSegments;
  fcInitParams->numChisqBins   = numChisqBins;
  fcInitParams->createRhosqVec = writeRhosq;
  fcInitParams->ovrlap         = ovrlap;
  fcInitParams->approximant    = approximant;
  fcInitParams->createCVec     = writeCData;


  /*
   *
   * power spectrum estimation
   *
   */


  /* initialize findchirp data conditioning routine */
  LAL_CALL( LALFindChirpDataInit( &status, &fcDataParams, fcInitParams ), 
      &status );
  fcDataParams->invSpecTrunc = invSpecTrunc * sampleRate;
  fcDataParams->fLow = fLow;

  /* create storage for the power spectral estimate */
  memset( &spec, 0, sizeof(REAL4FrequencySeries) );
  LAL_CALL( LALSCreateVector( &status, &(spec.data), numPoints / 2 + 1 ), 
      &status );

  /* compute the windowed power spectrum for the data channel */
  avgSpecParams.window = NULL;
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
    case 2:
      avgSpecParams.method = useUnity;
      if ( vrbflg ) fprintf( stdout, "simulation gaussian noise psd" );
      break;
  }

  /* use the fft plan created by findchirp */
  avgSpecParams.plan = fcDataParams->fwdPlan;

  wpars.type = Hann;
  wpars.length = numPoints;
  if ( badMeanPsd )
  {
    avgSpecParams.overlap = 0;
    if ( vrbflg ) fprintf( stdout, " without overlap\n" );
  }
  else
  {
    avgSpecParams.overlap = numPoints / 2;
    if ( vrbflg ) 
      fprintf( stdout, " with overlap %d\n", avgSpecParams.overlap );
  }

  LAL_CALL( LALCreateREAL4Window( &status, &(avgSpecParams.window),
        &wpars ), &status );
  LAL_CALL( LALREAL4AverageSpectrum( &status, &spec, &chan, &avgSpecParams ),
      &status );
  LAL_CALL( LALDestroyREAL4Window( &status, &(avgSpecParams.window) ), 
      &status );
  strcpy( spec.name, chan.name );

  if ( specType == 2 )
  {
    /* multiply the unit power spectrum to get a gaussian psd */
    REAL4 gaussVarSq = gaussVar * gaussVar;
    if ( resampleChan )
    {
      /* reduce the variance as we have resampled the data */
      gaussVarSq *= inputDeltaT / chan.deltaT;
    }
    for ( k = 0; k < spec.data->length; ++k )
    {
      spec.data->data[k] *= 2.0 * gaussVarSq * (REAL4) chan.deltaT;
    }

    if ( vrbflg ) 
      fprintf( stdout, "set psd to constant value = %e\n", spec.data->data[0] );
  }

  /* write the spectrum data to a file */
  if ( writeSpectrum ) outFrame = fr_add_proc_REAL4FrequencySeries( outFrame, 
      &spec, "ct^2/Hz", "PSD" );


  /*
   *
   * initialize the template bank simulation
   *
   */


  if ( bankSim )
  {
    if ( vrbflg ) fprintf( stdout, "initializing template bank simulation\n" );

    /* clear the output table */
    simResults.simInstParamsTable = NULL;

    /* override the findchirp initialization parameters */
    fcInitParams->numSegments  = 1;
    fcInitParams->ovrlap       = 0;

    /* set the thresholds for the bank sim */
    snrThresh    = 0;
    chisqThresh  = LAL_REAL4_MAX;

    /* set low frequency cutoff index of the psd and    */
    /* the value the psd at cut                         */
    cut = fLow / spec.deltaF > 1 ?  fLow / spec.deltaF : 1;
    if ( vrbflg ) 
      fprintf( stdout, "psd low frequency cutoff index = %d\n", cut );

    psdMin = spec.data->data[cut] * 
      ( ( resp.data->data[cut].re * resp.data->data[cut].re +
          resp.data->data[cut].im * resp.data->data[cut].im ) / psdScaleFac );

    /* calibrate the input power spectrum, scale to the */
    /* range of REAL4 and store the as S_v(f)           */
    for ( k = 0; k < cut; ++k )
    {
      spec.data->data[k] = psdMin;
    }
    for ( k = cut; k < spec.data->length; ++k )
    {
      REAL4 respRe = resp.data->data[k].re;
      REAL4 respIm = resp.data->data[k].im;
      spec.data->data[k] = spec.data->data[k] *
        ( ( respRe * respRe + respIm * respIm ) / psdScaleFac );
    }

    if ( writeSpectrum ) outFrame = fr_add_proc_REAL4FrequencySeries( outFrame, 
        &spec, "strain^2/Hz", "PSD_SIM" );

    /* set the response function to the sqrt of the inverse */ 
    /* of the psd scale factor since S_h = |R|^2 S_v        */
    for ( k = 0; k < resp.data->length; ++k )
    {
      resp.data->data[k].re = sqrt( psdScaleFac );
      resp.data->data[k].im = 0;
    }

    if ( writeResponse ) outFrame = fr_add_proc_COMPLEX8FrequencySeries( 
        outFrame, &resp, "strain/ct", "RESPONSE_SIM" );
  }


  /*
   *
   * create the data structures needed for findchirp
   *
   */


  if ( vrbflg ) fprintf( stdout, "initializing findchirp... " );

  /* create the data segment vector from the input data */
  LAL_CALL( LALInitializeDataSegmentVector( &status, &dataSegVec,
        &chan, &spec, &resp, fcInitParams ), &status );

  /* create the findchirp data storage */
  LAL_CALL( LALCreateFindChirpSegmentVector( &status, &fcSegVec, 
        fcInitParams ), &status );

  /* initialize the template functions */
  LAL_CALL( LALFindChirpTemplateInit( &status, &fcTmpltParams, 
        fcInitParams ), &status );

  fcDataParams->dynRange = fcTmpltParams->dynRange = dynRange;
  fcTmpltParams->deltaT = chan.deltaT;
  fcTmpltParams->fLow = fLow;

  /* initialize findchirp filter functions */
  LAL_CALL( LALFindChirpFilterInit( &status, &fcFilterParams, fcInitParams ), 
      &status );
  fcFilterParams->deltaT = chan.deltaT;
  fcFilterParams->chisqParams->approximant = approximant;

  /* set up parameters for the filter output veto */
#if 0
  fcFilterParams->filterOutputVetoParams = (FindChirpFilterOutputVetoParams *)
    LALCalloc( 1, sizeof(FindChirpFilterOutputVetoParams) );
#endif

  LAL_CALL( LALCreateFindChirpInput( &status, &fcFilterInput, fcInitParams ), 
      &status );

  LAL_CALL( LALFindChirpChisqVetoInit( &status, fcFilterParams->chisqParams, 
        fcInitParams->numChisqBins, fcInitParams->numPoints ), 
      &status );

  /* initialize findchirp clustering params */     /*XXX*/
  fcFilterParams->clusterMethod = clusterMethod;   /*XXX*/
  fcFilterParams->clusterWindow = clusterWindow;   /*XXX*/

  /* parse the thresholds */
  fcFilterParams->rhosqThresh = snrThresh * snrThresh;
  fcFilterParams->chisqThresh = chisqThresh;

  if ( vrbflg ) fprintf( stdout, "done\n" );


  /*
   *
   * matched filtering engine
   *
   */


  /* begin loop over number of template bank simulation */
  /* if we are not doing a template bank simulation,    */
  /* this exectues the main part of the filtering code  */
  /* just one (which is what we want to do)             */
  bankSimCount = 0;
  do
  {
    if ( bankSim )
    {

      if ( vrbflg ) 
        fprintf( stdout, "bank simulation %d/%d\n", bankSimCount, bankSim );


      /*
       *
       * inject a random signal if we are doing a template bank simulation
       *
       */


      if ( vrbflg ) fprintf( stdout, 
          "zeroing data stream and adding random injection for bank sim... " );

      /* zero out the input data segment and the injection params */
      memset( chan.data->data, 0, chan.data->length * sizeof(REAL4) );
      memset( &bankInjection, 0, sizeof(SimInspiralTable) );

      /* generate random parameters for the injection */
      LAL_CALL( LALUniformDeviate( &status, &(bankInjection.mass1), 
            randParams ), &status );
      bankInjection.mass1 *= (bankMaxMass - bankMinMass);
      bankInjection.mass1 += bankMinMass;
      LAL_CALL( LALUniformDeviate( &status, &(bankInjection.mass2), 
            randParams ), &status );
      bankInjection.mass2 *= (bankMaxMass - bankMinMass);
      bankInjection.mass2 += bankMinMass;
      bankInjection.eta = bankInjection.mass1 * bankInjection.mass2 /
        ( ( bankInjection.mass1 + bankInjection.mass2 ) *
          ( bankInjection.mass1 + bankInjection.mass2 ) );

      if ( bankSimApproximant == TaylorT2 )
      {
        LALSnprintf( bankInjection.waveform, LIGOMETA_WAVEFORM_MAX,
            "GeneratePPNtwoPN" );
      }
      else
      {
        fprintf( stderr, 
            "error: unknown waveform for bank simulation injection\n" );
        exit( 1 );
      }

      /* set the injection distance to 1 Mpc */
      bankInjection.distance = 1.0;

      /* inject the signals, preserving the channel name (Tev mangles it) */
      LALSnprintf( tmpChName, LALNameLength * sizeof(CHAR), "%s", 
          dataSegVec->data->chan->name );
      /* make sure the injection is hplus with no time delays */
      dataSegVec->data->chan->name[0] = 'P';
      LAL_CALL( LALFindChirpInjectSignals( &status, dataSegVec->data->chan, 
            &bankInjection, &resp ), &status );
      /* restore the saved channel name */
      LALSnprintf( dataSegVec->data->chan->name,  
          LALNameLength * sizeof(CHAR), "%s", tmpChName );

      /* write the channel data plus injection to the output frame file */
      if ( writeRawData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
          &chan, "ct", "SIM" );

      if ( vrbflg ) fprintf( stdout, "done\n" );
    }


    /*
     *
     * condition data segments for filtering
     *
     */


    switch ( approximant )
    {
      case TaylorT1:
      case TaylorT2:
      case TaylorT3:
      case GeneratePPN:
      case PadeT1:
      case EOB:
        if ( vrbflg ) 
          fprintf( stdout, "findchirp conditioning data for TD\n" );
        LAL_CALL( LALFindChirpTDData( &status, fcSegVec, dataSegVec, 
              fcDataParams ), &status );
        break;

      case TaylorF2:
        if ( vrbflg ) 
          fprintf( stdout, "findchirp conditioning data for SP\n" );
        LAL_CALL( LALFindChirpSPData( &status, fcSegVec, dataSegVec, 
              fcDataParams ), &status );
        break;

      case BCV:
        if ( vrbflg ) 
          fprintf( stdout, "findchirp conditioning data for BCV\n" );
        LAL_CALL( LALFindChirpBCVData( &status, fcSegVec, dataSegVec, 
              fcDataParams ), &status );
        break;

      case BCVSpin:
        if ( vrbflg ) 
          fprintf( stdout, "findchirp conditioning data for BCVSpin\n" );
        LAL_CALL( LALFindChirpBCVSpinData( &status, fcSegVec, dataSegVec, 
              fcDataParams ), &status );
        break;
      default:
        fprintf( stderr, "error: unknown waveform approximant for data\n" );
        exit( 1 );
        break;
    }

    if ( bankSim )
    {
      /* compute the minimal match normalization */
      REAL4 *tmpltPower = fcDataParams->tmpltPowerVec->data;
      COMPLEX8 *fcData = fcSegVec->data->data->data->data;

      if ( vrbflg ) fprintf( stdout,
          "computing minimal match normalization... " );

      matchNorm = 0;

      for ( k = 0; k < fcDataParams->tmpltPowerVec->length; ++k )
      {
        if ( tmpltPower[k] ) matchNorm += ( fcData[k].re * fcData[k].re +
            fcData[k].im * fcData[k].im ) / tmpltPower[k];
      }

      matchNorm *= ( 4.0 * (REAL4) fcSegVec->data->deltaT ) / 
        (REAL4) dataSegVec->data->chan->data->length;
      matchNorm = sqrt( matchNorm );

      if ( vrbflg ) fprintf( stdout, "%e\n", matchNorm );
    }
    else
    {
      if ( approximant == TaylorF2 ) 
      {
        /* compute the standard candle */
        REAL4 cannonDist = 1.0; /* Mpc */
        REAL4 m  = (REAL4) candle.tmplt.totalMass;
        REAL4 mu = (REAL4) candle.tmplt.mu;
        REAL4 distNorm = 2.0 * LAL_MRSUN_SI / (cannonDist * 1e6 * LAL_PC_SI);
        REAL4 candleTmpltNorm = sqrt( (5.0*mu) / 96.0 ) *
          pow( m / (LAL_PI*LAL_PI) , 1.0/3.0 ) *
          pow( LAL_MTSUN_SI / (REAL4) chan.deltaT, -1.0/6.0 );

        distNorm *= fcTmpltParams->dynRange;
        candleTmpltNorm *= candleTmpltNorm;
        candleTmpltNorm *= distNorm * distNorm;

        candle.sigmasq = 4.0 * ( (REAL4) chan.deltaT / (REAL4) numPoints );
        candle.sigmasq *= candleTmpltNorm * 
          fcSegVec->data->segNorm->data[fcSegVec->data->segNorm->length-1];

        candle.distance = sqrt( candle.sigmasq / candle.rhosq );

        if ( vrbflg ) 
        {
          fprintf( stdout, "candle m = %e\ncandle mu = %e\n"
              "candle.rhosq = %e\nchan.deltaT = %e\nnumPoints = %d\n"
              "fcSegVec->data->segNorm->data[fcSegVec->data->segNorm->length-1]"
              " = %e\ncandleTmpltNorm = %e\ncandle.distance = %e Mpc\n"
              "candle.sigmasq = %e\n",
              m, mu, candle.rhosq, chan.deltaT, numPoints, 
              fcSegVec->data->segNorm->data[fcSegVec->data->segNorm->length-1], 
              candleTmpltNorm, candle.distance, candle.sigmasq );
          fflush( stdout );
        }
      }
      else if ( approximant == BCV )
      {
        /* compute the standard candle for a  5-5 Msun inspiral*/
        REAL4 cannonDist = 1.0; /* Mpc */
        REAL4 m  = 10.0;
        REAL4 mu = 2.5;
        REAL4 k1 = numPoints * chan.deltaT ;
        UINT4 kmax = 432 * k1 ; /* 432 = fISCO for 5-5 */
        REAL4 distNorm = 2.0 * LAL_MRSUN_SI / (cannonDist * 1e6 * LAL_PC_SI);
        REAL4 candleTmpltNorm = sqrt( (5.0*mu) / 96.0 ) *
          pow( m/(LAL_PI*LAL_PI) , 1.0/3.0 ) * 
          pow(LAL_MTSUN_SI / (REAL4) chan.deltaT, -1.0/6.0);
        distNorm *= fcTmpltParams->dynRange;
        candleTmpltNorm *= candleTmpltNorm;
        candleTmpltNorm *= distNorm * distNorm;
        candle.sigmasq = 4.0 * ( (REAL4) chan.deltaT / (REAL4) numPoints );
        candle.sigmasq *= candleTmpltNorm *
          fcSegVec->data->segNorm->data[kmax];
        candle.distance = sqrt( candle.sigmasq / candle.rhosq );

        /* for 5-5 Msun... */
        candle.tmplt.mass1 = 5.0;
        candle.tmplt.mass2 = 5.0;
        candle.tmplt.totalMass = 10.0;
        candle.tmplt.mu = 2.5;
        candle.tmplt.eta = 0.25;

        if ( vrbflg ) 
        {
          fprintf( stdout, "candle m = %e\ncandle mu = %e\n"
              "candle.rhosq = %e\nchan.deltaT = %e\nnumPoints = %d\n"
              "fcSegVec->data->segNorm->data[kmax] = %e\n"
              "kmax = %d\ncandleTmpltNorm = %e\ncandle.distance = %e Mpc \n"
              "candle.sigmasq=%e\n",
              m,mu,candle.rhosq,chan.deltaT,
              numPoints,fcSegVec->data->segNorm->data[kmax],kmax,
              candleTmpltNorm,candle.distance,candle.sigmasq);
          fflush(stdout);
        }
      }
      else 
      {
        if ( vrbflg )
        {
          fprintf( stdout, "standard candle not calculated; \n"
              "chan.deltaT = %e\nnumPoints = %d\n",
              chan.deltaT, numPoints );
          fflush( stdout );
        }
      }  
    }


    /*
     *
     * hierarchial search engine
     *
     */


    for ( tmpltCurrent = tmpltHead, inserted = 0; tmpltCurrent; 
        tmpltCurrent = tmpltCurrent->next, inserted = 0 )
    {
      /*  generate template */
      switch ( approximant )
      {
        case TaylorT1:
        case TaylorT2:
        case TaylorT3:
        case GeneratePPN:
        case PadeT1:
        case EOB:
          LAL_CALL( LALFindChirpTDTemplate( &status, fcFilterInput->fcTmplt, 
                tmpltCurrent->tmpltPtr, fcTmpltParams ), &status );
          break;

        case TaylorF2:
          LAL_CALL( LALFindChirpSPTemplate( &status, fcFilterInput->fcTmplt, 
                tmpltCurrent->tmpltPtr, fcTmpltParams ), &status );
          break;

        case BCV:
          LAL_CALL( LALFindChirpBCVTemplate( &status, fcFilterInput->fcTmplt, 
                tmpltCurrent->tmpltPtr, fcTmpltParams ), &status );
          break;

        case BCVSpin:
          LAL_CALL( LALFindChirpBCVSpinTemplate( &status, 
                fcFilterInput->fcTmplt, tmpltCurrent->tmpltPtr, 
                fcTmpltParams, fcDataParams ), &status );
          break;

        default:
          fprintf( stderr, "error: unknown waveform template approximant \n" );
          exit( 1 );
          break;
      }

      /* loop over data segments */
      for ( i = 0; i < fcSegVec->length ; ++i )
      {
        INT8 fcSegStartTimeNS;
        INT8 fcSegEndTimeNS;

        LAL_CALL( LALGPStoINT8( &status, &fcSegStartTimeNS, 
              &(fcSegVec->data[i].data->epoch) ), &status );
        fcSegEndTimeNS = fcSegStartTimeNS + (INT8)
          ( (REAL8) numPoints * 1e9 * fcSegVec->data[i].deltaT );

        /* skip segment if it is not contained in the trig start or end times */
        if ( (trigStartTimeNS && (trigStartTimeNS > fcSegEndTimeNS)) || 
            (trigEndTimeNS && (trigEndTimeNS < fcSegStartTimeNS)) )
        { 
          if ( vrbflg ) fprintf( stdout, 
              "skipping segment %d/%d [%lld-%lld] (outside trig time)\n", 
              fcSegVec->data[i].number, fcSegVec->length, 
              fcSegStartTimeNS, fcSegEndTimeNS );

          continue;
        }

        /* filter data segment */ 
        if ( fcSegVec->data[i].level == tmpltCurrent->tmpltPtr->level )
        {
          if ( vrbflg ) fprintf( stdout, 
              "filtering segment %d/%d [%lld-%lld] "
              "against template %d/%d (%e,%e)\n", 
              fcSegVec->data[i].number,  fcSegVec->length,
              fcSegStartTimeNS, fcSegEndTimeNS,
              tmpltCurrent->tmpltPtr->number, numTmplts,
              fcFilterInput->fcTmplt->tmplt.mass1, 
              fcFilterInput->fcTmplt->tmplt.mass2 );

          fcFilterInput->segment = fcSegVec->data + i;

          /* decide which filtering routine to use */
          switch ( approximant )
          {
            case TaylorT1:
            case TaylorT2:
            case TaylorT3:
            case GeneratePPN:
            case PadeT1:
            case EOB:
              /* construct normalization for time domain templates... */
              LAL_CALL( LALFindChirpTDNormalize( &status, 
                    fcFilterInput->fcTmplt, fcFilterInput->segment, 
                    fcDataParams ), &status );
              /* ...and fall through to FindChirpFilterSegment() */
            case TaylorF2:
              LAL_CALL( LALFindChirpFilterSegment( &status, 
                    &eventList, fcFilterInput, fcFilterParams ), &status ); 
              break;
              
            case BCV:
              if (!bcvConstraint){
		      LAL_CALL( LALFindChirpBCVFilterSegment( &status,
                    &eventList, fcFilterInput, fcFilterParams ), &status );
	      }
	      else
	      { 
		      LAL_CALL( LALFindChirpBCVCFilterSegment( &status,
                    &eventList, fcFilterInput, fcFilterParams ), &status );
	      }
              break;
              
            case BCVSpin:
              LAL_CALL( LALFindChirpBCVSpinFilterSegment( &status,
                    &eventList, fcFilterInput, fcFilterParams, fcDataParams 
                    ), &status );
              break;

            default:
              fprintf( stderr, 
                  "error: unknown waveform approximant for filter\n" );
              exit( 1 );
              break;
          }

          if ( writeRhosq )
          {
            CHAR snrsqStr[LALNameLength];
            LALSnprintf( snrsqStr, LALNameLength*sizeof(CHAR), 
                "SNRSQ_%d", nRhosqFr++ );
            strcpy( fcFilterParams->rhosqVec->name, chan.name );
            outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
                fcFilterParams->rhosqVec, "none", snrsqStr );
          }

          if ( writeCData && ! strcmp(ifo,tmpltCurrent->tmpltPtr->ifo) )
          {
	    cDataForFrame = 0;
	    trigTime = tmpltCurrent->tmpltPtr->end_time.gpsSeconds + 1e-9 * tmpltCurrent->tmpltPtr->end_time.gpsNanoSeconds;
	    lowerBound = gpsStartTime.gpsSeconds + numPoints/(4 * sampleRate );
	    upperBound = gpsEndTime.gpsSeconds - numPoints/(4 * sampleRate );

	    if( trigTime < lowerBound || trigTime > upperBound )
	      {
		fprintf(stderr,"The trigger time is outside of the segment\n");
		fprintf(stderr,"Not writing C-data for this segment\n");
		goto noCdataLoopExitPoint;
	      }
 
	    tempTmplt = (SnglInspiralTable *) 
                 LALCalloc(1, sizeof(SnglInspiralTable) );
	    tempTmplt->event_id = (EventIDColumn *) 
                 LALCalloc(1, sizeof(EventIDColumn) );
	    tempTmplt->mass1 = tmpltCurrent->tmpltPtr->mass1;
	    tempTmplt->end_time.gpsSeconds = tmpltCurrent->tmpltPtr->end_time.gpsSeconds;
	    tempTmplt->end_time.gpsNanoSeconds = tmpltCurrent->tmpltPtr->end_time.gpsNanoSeconds;
	    tempTmplt->event_id->id = tmpltCurrent->tmpltPtr->event_id->id;
 
           LALFindChirpCreateCoherentInput( &status,
                  &coherentInputData, fcFilterParams->cVec, 
                  tempTmplt, 2.0, numPoints / 4 );

	    if( coherentInputData )
	      {
		cDataForFrame = 1;
		LALSnprintf( cdataStr, LALNameLength*sizeof(CHAR),
		       "CData_%d", nCDataFr++ );
		strcpy( coherentInputData->name, chan.name );
		outFrame = fr_add_proc_COMPLEX8TimeSeries( outFrame,
		       coherentInputData, "none", cdataStr );
		LAL_CALL( LALCDestroyVector( &status, 
                       &(coherentInputData->data) ), &status );
		coherentInputData = NULL;
	      }
          }

	noCdataLoopExitPoint:

          if ( writeChisq )
          {
            CHAR chisqStr[LALNameLength];
            REAL4TimeSeries chisqts;
            LALSnprintf( chisqStr, LALNameLength*sizeof(CHAR), 
                "CHISQ_%d", nChisqFr++ );
            chisqts.epoch = fcFilterInput->segment->data->epoch;
            memcpy( &(chisqts.name), fcFilterInput->segment->data->name,
                LALNameLength * sizeof(CHAR) );
            chisqts.deltaT = fcFilterInput->segment->deltaT;
            chisqts.data = fcFilterParams->chisqVec;
            outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
                &chisqts, "none", chisqStr );
          }
if ( vrbflg ) fprintf (stdout, "epoch = %d\n",fcFilterInput->segment->data->epoch );
        }
        else
        {
          if ( vrbflg ) fprintf( stdout, "skipping segment %d/%d [%lld-%lld] "
              "(segment level %d, template level %d)\n", 
              fcSegVec->data[i].number, fcSegVec->length, 
              fcSegStartTimeNS, fcSegEndTimeNS,
              fcSegVec->data[i].level, tmpltCurrent->tmpltPtr->level );
        }

        /*  test if filter returned any events */
        if ( eventList )
        {
          if ( vrbflg ) fprintf( stdout, 
              "segment %d rang template [m (%e,%e)] [psi (%e,%e)]\n",
              fcSegVec->data[i].number,
              fcFilterInput->fcTmplt->tmplt.mass1, 
              fcFilterInput->fcTmplt->tmplt.mass2,
              fcFilterInput->fcTmplt->tmplt.psi0, 
              fcFilterInput->fcTmplt->tmplt.psi3 );

          if ( tmpltCurrent->tmpltPtr->fine != NULL && inserted == 0 )
          {
            if ( vrbflg ) fprintf( stdout, 
                "inserting fine templates into list\n" );

            tmpltInsert = tmpltCurrent;
            inserted = 1;
            fcSegVec->data[i].level += 1;

            for ( bankCurrent = tmpltCurrent->tmpltPtr->fine ; 
                bankCurrent; bankCurrent = bankCurrent->next )
            {
              LAL_CALL( LALFindChirpCreateTmpltNode( &status, 
                    bankCurrent, &tmpltInsert ), &status );
            }

          }
          else if ( ! tmpltCurrent->tmpltPtr->fine )
          {
            if ( vrbflg ) fprintf( stdout, "***>  dumping events  <***\n" );

            if ( ! savedEvents.snglInspiralTable )
            {
              savedEvents.snglInspiralTable = eventList;
            }
            else
            {
              event->next = eventList;
            }
          }
          else
          {
            if ( vrbflg ) fprintf( stdout, 
                "already inserted fine templates, skipping\n" ); 

            fcSegVec->data[i].level += 1;
          } 

          /* save a ptr to the last event in the list and count the events */
          ++numEvents;
          while ( eventList->next )
          {
            eventList = eventList->next;
            ++numEvents;
          }
          event = eventList;
          eventList = NULL;
        } /* end if ( events ) */

        /* if going up a level, remove inserted nodes, reset segment levels */ 
        if ( tmpltCurrent->next && (tmpltCurrent->next->tmpltPtr->level < 
              tmpltCurrent->tmpltPtr->level) )
        {
          /* record the current number */
          currentLevel = tmpltCurrent->tmpltPtr->level;

          /* decrease segment filter levels if the have been increased */
          for ( i = 0 ; i < fcSegVec->length; i++ )
          {
            if ( fcSegVec->data[i].level == currentLevel )
            {
              fcSegVec->data[i].level -= 1;
            }
          }

          if ( vrbflg ) fprintf( stdout, "removing inserted fine templates\n" );

          while ( tmpltCurrent->tmpltPtr->level == currentLevel )
          {
            LAL_CALL( LALFindChirpDestroyTmpltNode( &status, &tmpltCurrent ),
                &status );
          }          
        } /* end if up a level */

      } /* end loop over data segments */

      /* delete the chisq bins for time domain templates */
      switch ( approximant )
      {
        case TaylorT1:
        case TaylorT2:
        case TaylorT3:
        case GeneratePPN:
        case PadeT1:
        case EOB:
          /* the chisq bins need to be re-computed for the next template */
          for ( i = 0; i < fcSegVec->length ; ++i )
          {
            if ( fcSegVec->data[i].chisqBinVec->data )
            {
              LALFree( fcSegVec->data[i].chisqBinVec->data );
              fcSegVec->data[i].chisqBinVec->data = NULL;
            }
          }
        default:
          break;
      }

    } /* end loop over linked list */

    if ( bankSim )
    {
      /* allocate memory for the loudest event over the template bank */
      loudestEvent = (SnglInspiralTable *) 
        LALCalloc( 1, sizeof(SnglInspiralTable) );

      /* find the loudest snr over the template bank */
      event = savedEvents.snglInspiralTable; 
      while ( event )
      {
        SnglInspiralTable *tmpEvent = event;

        if ( event->snr > loudestEvent->snr )
        {
          memcpy( loudestEvent, event, sizeof(SnglInspiralTable) );
          loudestEvent->next = NULL;
        }

        event = event->next;
        LALFree( tmpEvent );
      }
      savedEvents.snglInspiralTable = NULL;

      /* link the list of loudest events */
      if ( ! loudestEventHead ) loudestEventHead = loudestEvent;
      if ( prevLoudestEvent ) prevLoudestEvent->next = loudestEvent;
      prevLoudestEvent = loudestEvent;

      /* create sim_inst_params structure for mass1 */
      thisSimInstParams = (SimInstParamsTable *) 
        LALCalloc( 1, sizeof(SimInstParamsTable) );
      LALSnprintf( thisSimInstParams->name, LIGOMETA_SIMINSTPARAMS_NAME_MAX,
          "mass1" );
      thisSimInstParams->value = bankInjection.mass1;

      /* link the linked list */
      if ( ! simResults.simInstParamsTable )
        simResults.simInstParamsTable = thisSimInstParams;
      if ( prevSimInstParams ) prevSimInstParams->next = thisSimInstParams;

      /* create sim_inst_params structure for mass2 */
      thisSimInstParams = thisSimInstParams->next = (SimInstParamsTable *) 
        LALCalloc( 1, sizeof(SimInstParamsTable) );
      LALSnprintf( thisSimInstParams->name, LIGOMETA_SIMINSTPARAMS_NAME_MAX,
          "mass2" );
      thisSimInstParams->value = bankInjection.mass2;

      /* create sim_inst_params structure for minimal_match */
      thisSimInstParams = thisSimInstParams->next = (SimInstParamsTable *) 
        LALCalloc( 1, sizeof(SimInstParamsTable) );
      LALSnprintf( thisSimInstParams->name, LIGOMETA_SIMINSTPARAMS_NAME_MAX,
          "minimal_match" );
      thisSimInstParams->value = loudestEvent->snr / matchNorm;

      /* store the last created sim_inst_params table */
      prevSimInstParams = thisSimInstParams;
    }

    ++bankSimCount;

  } while ( bankSimCount < bankSim ); /* end loop over bank simulations */

  if ( bankSim )
  {
    /* save the number of bank simulation in the search summary table */
    searchsumm.searchSummaryTable->nevents = bankSim;

    /* point the saved events to the linked list of loudest events */
    savedEvents.snglInspiralTable = loudestEventHead;
  }
  else
  {
    /* save the number of events in the search summary table */
    searchsumm.searchSummaryTable->nevents = numEvents;
  }


  /*
   *
   * free memory used by filtering code
   *
   */


  if ( vrbflg ) fprintf( stdout, "freeing memory\n" );

  /* free memory used by findchirp */
  LAL_CALL( LALFindChirpChisqVetoFinalize( &status, 
        fcFilterParams->chisqParams, fcInitParams->numChisqBins ), 
      &status );
  LAL_CALL( LALDestroyFindChirpInput( &status, &fcFilterInput ), 
      &status );
  LAL_CALL( LALFindChirpFilterFinalize( &status, &fcFilterParams ), 
      &status );
  LAL_CALL( LALFindChirpTemplateFinalize( &status, &fcTmpltParams ), 
      &status );
  LAL_CALL( LALFindChirpDataFinalize( &status, &fcDataParams ),
      &status );
  LAL_CALL( LALDestroyFindChirpSegmentVector( &status, &fcSegVec ),
      &status );
  LALFree( fcInitParams );

  /* free the template bank */
  while ( bankHead )
  {
    bankCurrent = bankHead;
    bankHead = bankHead->next;
    LALFree( bankCurrent );
    bankCurrent = NULL;
  }
  /* destroy linked list of template nodes */
  while ( tmpltHead )
  {
    LAL_CALL( LALFindChirpDestroyTmpltNode( &status, &tmpltHead ), &status );
  }

  /* free the data storage */
  LAL_CALL( LALFinalizeDataSegmentVector( &status, &dataSegVec ), &status );
  LAL_CALL( LALSDestroyVector( &status, &(chan.data) ), &status );
  LAL_CALL( LALSDestroyVector( &status, &(spec.data) ), &status );
  LAL_CALL( LALCDestroyVector( &status, &(resp.data) ), &status );

  /* free the random parameters structure */
  if ( randSeedType != unset )
  {
    LAL_CALL( LALDestroyRandomParams( &status, &randParams ), &status );
  }


  /*
   *
   * write the result results to disk
   *
   */


  /* write the output frame */
  if ( writeRawData || writeFilterData || writeResponse || writeSpectrum ||
      writeRhosq || writeChisq || (writeCData && cDataForFrame) )
  {
    if ( outputPath[0] )
    {
      LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), "%s/%s.gwf", 
          outputPath, fileName );
    }
    else
    {
      LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), "%s.gwf", fileName );
    }
    if ( vrbflg ) fprintf( stdout, "writing frame data to %s... ", fname );
    frOutFile = FrFileONew( fname, 0 );
    FrameWrite( outFrame, frOutFile );
    FrFileOEnd( frOutFile );
    if ( vrbflg ) fprintf( stdout, "done\n" );
  }

  /* open the output xml file */
  memset( &results, 0, sizeof(LIGOLwXMLStream) );
  if ( outputPath[0] )
  {
    LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), "%s/%s.xml", 
        outputPath, fileName );
  }
  else
  {
    LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), "%s.xml", fileName );
  }
  if ( vrbflg ) fprintf( stdout, "writing XML data to %s...\n", fname );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &results, fname), &status );

  /* write the process table */
  if ( vrbflg ) fprintf( stdout, "  process table...\n" );
  LALSnprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s", ifo );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->end_time),
        &accuracy ), &status );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, proctable, 
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  free( proctable.processTable );

  /* write the process params table */
  if ( vrbflg ) fprintf( stdout, "  process_params table...\n" );
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
  if ( vrbflg ) fprintf( stdout, "  search_summary table...\n" );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, 
        search_summary_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, searchsumm, 
        search_summary_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );

  /* write the search summvars table */
  if ( numTmplts )
  {
    if ( vrbflg ) fprintf( stdout, "  search_summvars table...\n" );
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, 
          search_summvars_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, searchsummvars, 
          search_summvars_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  }
  while( searchsummvars.searchSummvarsTable )
  {
    this_search_summvar = searchsummvars.searchSummvarsTable;
    searchsummvars.searchSummvarsTable = this_search_summvar->next;
    LALFree( this_search_summvar );
  }

  /* write the summ_value table with the standard candle distance */
  if ( ! bankSim )
  {
    if ( approximant == TaylorF2 )
    {
      if ( vrbflg ) fprintf( stdout, "  summ_value table...\n" );
      ADD_SUMM_VALUE( "inspiral_effective_distance", "1.4_1.4_8", 
          candle.distance, 0);
    }
    else if ( approximant == BCV )
    {
      if ( vrbflg ) fprintf( stdout, "  summ_value table...\n" );
      ADD_SUMM_VALUE( "inspiral_effective_distance", "5.0_5.0_8",
          candle.distance, 0);
    }
  }

  /* store calibration information */
  ADD_SUMM_VALUE( "calibration alpha", "analysis", alpha, 0 );
  ADD_SUMM_VALUE( "calibration alphabeta", "analysis", alphabeta, 0 );
  if (injectionFile) 
  {
    ADD_SUMM_VALUE( "calibration alpha", "injection", inj_alpha, 0 );
    ADD_SUMM_VALUE( "calibration alphabeta", "injection", inj_alphabeta, 0 );
  }
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, summ_value_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, summvalue, 
        summ_value_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );


  while ( summvalue.summValueTable )
  {
    this_summ_value = summvalue.summValueTable;
    summvalue.summValueTable = summvalue.summValueTable->next;
    LALFree( this_summ_value );
  }

  /* free the search summary table */
  free( searchsumm.searchSummaryTable );

  /* write the sngl_inspiral triggers to the output xml */
  if ( savedEvents.snglInspiralTable )
  {
    SnglInspiralTable *tmpEventHead = NULL;
    SnglInspiralTable *lastEvent = NULL;

    /* sort the inspiral events by time */
    if ( vrbflg ) fprintf( stdout, "  sorting events by time... " );
    LAL_CALL( LALSortSnglInspiral( &status, &(savedEvents.snglInspiralTable),
          LALCompareSnglInspiralByTime), &status );
    if ( vrbflg ) fprintf( stdout, "done\n" );

    /* discard any triggers outside the trig start/end time window */
    event = savedEvents.snglInspiralTable;
    if ( trigStartTimeNS || trigEndTimeNS )
    {
      if ( vrbflg ) fprintf( stdout, 
          "  discarding triggers outside trig start/end time... " );

      while ( event )
      {
        INT8 trigTimeNS;
        LAL_CALL( LALGPStoINT8( &status, &trigTimeNS, &(event->end_time) ), 
            &status );

        if ( trigTimeNS &&
            ((trigStartTimeNS && (trigTimeNS < trigStartTimeNS)) ||
             (trigEndTimeNS && (trigTimeNS >= trigEndTimeNS))) )
        {
          /* throw this trigger away */
          SnglInspiralTable *tmpEvent = event;

          if ( lastEvent )
          {
            lastEvent->next = event->next;
          }

          /* increment the linked list by one and free the event */
          event = event->next;
          LALFree( tmpEvent );
        }
        else 
        {
          /* store the first event as the head of the new linked list */
          if ( ! tmpEventHead ) tmpEventHead = event;

          /* save the last event and increment the linked list by one */
          lastEvent = event;
          event = event->next;
        }
      }

      savedEvents.snglInspiralTable = tmpEventHead;

      if ( vrbflg ) fprintf( stdout, "done\n" );
    }

    /* if we haven't thrown all the triggers away, write sngl_inspiral table */
    if ( savedEvents.snglInspiralTable )
    {
      if ( vrbflg ) fprintf( stdout, "  sngl_inspiral table...\n" );
      LAL_CALL( LALBeginLIGOLwXMLTable( &status, 
            &results, sngl_inspiral_table ), &status );
      LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, savedEvents, 
            sngl_inspiral_table ), &status );
      LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
    }
  }
  while ( savedEvents.snglInspiralTable )
  {
    event = savedEvents.snglInspiralTable;
    savedEvents.snglInspiralTable = savedEvents.snglInspiralTable->next;
    LALFree( event );
  }

  /* write the template bank simulation results to the xml */
  if ( bankSim && simResults.simInstParamsTable )
  {
    if ( vrbflg ) fprintf( stdout, "  sim_inst table...\n" );
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, 
          &results, sim_inst_params_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, simResults, 
          sim_inst_params_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );

    while ( simResults.simInstParamsTable )
    {
      thisSimInstParams = simResults.simInstParamsTable;
      simResults.simInstParamsTable = simResults.simInstParamsTable->next;
      LALFree( thisSimInstParams );
    }
  }

  /* close the output xml file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &results ), &status );
  if ( vrbflg ) fprintf( stdout, "done. XML file closed\n" );

  /* free the rest of the memory, check for memory leaks and exit */
  if ( injectionFile ) free ( injectionFile ); 
  if ( calCacheName ) free( calCacheName );
  if ( frInCacheName ) free( frInCacheName );
  if ( frInType ) free( frInType );
  if ( bankFileName ) free( bankFileName );
  if ( channelName ) free( channelName );
  if ( fqChanName ) free( fqChanName );

  if ( vrbflg ) fprintf( stdout, "checking memory leaks and exiting\n" );
  LALCheckMemoryLeaks();
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
"lalapps_inspiral [options]\n\n"\
"  --help                       display this message\n"\
"  --verbose                    print progress information\n"\
"  --version                    print version information and exit\n"\
"  --debug-level LEVEL          set the LAL debug level to LEVEL\n"\
"  --user-tag STRING            set the process_params usertag to STRING\n"\
"  --ifo-tag STRING             set the ifotag to STRING - for file naming\n"\
"  --comment STRING             set the process table comment to STRING\n"\
"\n"\
"  --gps-start-time SEC         GPS second of data start time\n"\
"  --gps-start-time-ns NS       GPS nanosecond of data start time\n"\
"  --gps-end-time SEC           GPS second of data end time\n"\
"  --gps-end-time-ns NS         GPS nanosecond of data end time\n"\
"  --pad-data T                 pad the data start and end time by T seconds\n"\
"  --slide-time T               slide data start epoch by T seconds\n"\
"  --slide-time-ns T            slide data start epoch by T nanoseconds\n"\
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
"  --injection-file FILE        inject simulated inspiral signals from FILE\n"\
"\n"\
"  --inject-overhead            inject signals from overhead detector\n"\
"  --bank-file FILE             read template bank parameters from FILE\n"\
"  --minimal-match M            override bank minimal match with M (sets delta)\n"\
"  --start-template N           start filtering at template number N in bank\n"\
"  --stop-templateN             stop filtering at template number N in bank\n"\
"\n"\
"  --sample-rate F              filter data at F Hz, downsampling if necessary\n"\
"  --resample-filter TYPE       set resample filter to TYPE (ldas|butterworth)\n"\
"\n"\
"  --disable-high-pass          turn off the IIR highpass filter\n"\
"  --enable-high-pass F         high pass data above F Hz using an IIR filter\n"\
"  --high-pass-order O          set the order of the high pass filter to O\n"\
"  --high-pass-attenuation A    set the attenuation of the high pass filter to A\n"\
"  --spectrum-type TYPE         use PSD estimator TYPE (mean|median|gaussian)\n"\
"\n"\
"  --segment-length N           set data segment length to N points\n"\
"  --number-of-segments N       set number of data segments to N\n"\
"  --segment-overlap N          overlap data segments by N points\n"\
"\n"\
"  --low-frequency-cutoff F     do not filter below F Hz\n"\
"  --inverse-spec-length T      set length of inverse spectrum to T seconds\n"\
"  --dynamic-range-exponent X   set dynamic range scaling to 2^X\n"\
"\n"\
"  --approximant APPROX         set approximant of the waveform to APPROX\n"\
"                                 (TaylorF2|BCV|BCVSpin)\n"\
"  --chisq-bins P               set number of chisq veto bins to P\n"\
"  --snr-threshold RHO          set signal-to-noise threshold to RHO\n"\
"  --chisq-threshold X          threshold on chi^2 < X * ( p + rho^2 * delta^2 )\n"\
"  --cluster-method MTHD        set maximize over chirp MTHD (tmplt|window|noClustering)\n"\
"  --cluster-window SEC         set length of clustering time window if required\n"\
"\n"\
"  --enable-output              write the results to a LIGO LW XML file\n"\
"  --disable-output             do not write LIGO LW XML output file\n"\
"  --trig-start-time SEC        only output triggers after GPS time SEC\n"\
"  --trig-end-time SEC          only output triggers before GPS time SEC\n"\
"\n"\
"  --gaussian-noise VAR         replace data with gaussian noise of variance VAR\n"\
"  --random-seed SEED           set random number seed for injections to SEED\n"\
"                                 (urandom|integer)\n"\
"\n"\
"  --bank-simulation N          perform N injections to test the template bank\n"\
"  --sim-approximant APX        set approximant of the injected waveform to APX\n"\
"  --sim-minimum-mass M         set minimum mass of bank injected signal to M\n"\
"  --sim-maximum-mass M         set maximum mass of bank injected signal to M\n"\
"\n"\
"  --data-checkpoint            checkpoint and exit after data is read in\n"\
"  --checkpoint-path PATH       write checkpoint file under PATH\n"\
"  --output-path PATH           write output data to PATH\n"\
"\n"\
"  --write-raw-data             write raw data to a frame file\n"\
"  --write-filter-data          write data that is passed to filter to a frame\n"\
"  --write-response             write the computed response function to a frame\n"\
"  --write-spectrum             write the uncalibrated psd to a frame\n"\
"  --write-snrsq                write the snr time series for each data segment\n"\
"  --write-chisq                write the r^2 time series for each data segment\n"\
"  --write-cdata                write the complex filter output\n"\
"\n"

int arg_parse_check( int argc, char *argv[], MetadataTable procparams )
{
  /* getopt arguments */
  struct option long_options[] =
  {
    /* these options set a flag */
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"enable-output",           no_argument,       &enableOutput,     1 },
    {"disable-output",          no_argument,       &enableOutput,     0 },
    {"disable-high-pass",       no_argument,       &highPass,         0 },
    {"inject-overhead",         no_argument,       &injectOverhead,   1 },
    {"data-checkpoint",         no_argument,       &dataCheckpoint,   1 },
    {"glob-frame-data",         no_argument,       &globFrameData,    1 },
    {"glob-calibration-data",   no_argument,       &globCalData,      1 },
    {"point-calibration",       no_argument,       &pointCal,         1 },
    /* these options don't set a flag */
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-start-time-ns",       required_argument, 0,                'A'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"gps-end-time-ns",         required_argument, 0,                'B'},
    {"channel-name",            required_argument, 0,                'c'},
    {"segment-length",          required_argument, 0,                'd'},
    {"trig-start-time",         required_argument, 0,                'C'},
    {"trig-end-time",           required_argument, 0,                'D'},
    {"number-of-segments",      required_argument, 0,                'e'},
    {"segment-overlap",         required_argument, 0,                'f'},
    {"sample-rate",             required_argument, 0,                'g'},
    {"calibrated-data",         required_argument, 0,                'y'},
    {"geo-high-pass-freq",      required_argument, 0,                'E'},
    {"geo-high-pass-order",     required_argument, 0,                'P'},
    {"geo-high-pass-atten",     required_argument, 0,                'Q'},
    {"help",                    no_argument,       0,                'h'},
    {"low-frequency-cutoff",    required_argument, 0,                'i'},
    {"spectrum-type",           required_argument, 0,                'j'},
    {"inverse-spec-length",     required_argument, 0,                'k'},
    {"dynamic-range-exponent",  required_argument, 0,                'l'},
    {"start-template",          required_argument, 0,                'm'},
    {"minimal-match",           required_argument, 0,                'M'},
    {"stop-template",           required_argument, 0,                'n'},
    {"chisq-bins",              required_argument, 0,                'o'},
    {"calibration-cache",       required_argument, 0,                'p'},
    {"approximant",             required_argument, 0,                'F'},
    {"snr-threshold",           required_argument, 0,                'q'},
    {"chisq-threshold",         required_argument, 0,                'r'},
    {"resample-filter",         required_argument, 0,                'R'},
    {"comment",                 required_argument, 0,                's'},
    {"enable-high-pass",        required_argument, 0,                't'},
    {"high-pass-order",         required_argument, 0,                'H'},
    {"high-pass-attenuation",   required_argument, 0,                'T'},
    {"frame-cache",             required_argument, 0,                'u'},
    {"frame-type",              required_argument, 0,                'S'},
    {"bank-file",               required_argument, 0,                'v'},
    {"injection-file",          required_argument, 0,                'w'},
    {"pad-data",                required_argument, 0,                'x'},
    {"slide-time",              required_argument, 0,                'X'},
    {"slide-time-ns",           required_argument, 0,                'Y'},
    {"bank-simulation",         required_argument, 0,                'K'},
    {"sim-approximant",         required_argument, 0,                'L'},
    {"random-seed",             required_argument, 0,                'J'},
    {"sim-minimum-mass",        required_argument, 0,                'U'},
    {"sim-maximum-mass",        required_argument, 0,                'W'},
    {"gaussian-noise",          required_argument, 0,                'G'},
    {"checkpoint-path",         required_argument, 0,                'N'},
    {"output-path",             required_argument, 0,                'O'},
    {"debug-level",             required_argument, 0,                'z'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {"ifo-tag",                 required_argument, 0,                'I'},
    {"version",                 no_argument,       0,                'V'},
    {"cluster-method",          required_argument, 0,                '*'},  
    {"cluster-window",          required_argument, 0,                '#'},  
    /* frame writing options */
    {"write-raw-data",          no_argument,       &writeRawData,     1 },
    {"write-filter-data",       no_argument,       &writeFilterData,  1 },
    {"write-response",          no_argument,       &writeResponse,    1 },
    {"write-spectrum",          no_argument,       &writeSpectrum,    1 },
    {"write-snrsq",             no_argument,       &writeRhosq,       1 },
    {"write-chisq",             no_argument,       &writeChisq,       1 },
    {"write-cdata",             no_argument,       &writeCData,       1 },
    {0, 0, 0, 0}
  };
  int c;
  INT4 haveDynRange = 0;
  INT4 haveApprox = 0;
  INT4 haveClusterMethod = 0;
  INT4 haveBankSimApprox = 0;
  ProcessParamsTable *this_proc_param = procparams.processParamsTable;
  LALStatus             status = blank_status;


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
        "A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:V:W:X:Y:Z:"
        "a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:",
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
          gpsStartTimeNS += (INT8) gstartt * 1000000000LL;
          ADD_PROCESS_PARAM( "int", "%ld", gstartt );
        }
        break;

      case 'A':
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
          gpsStartTimeNS += (INT8) gstarttns;
          ADD_PROCESS_PARAM( "int", "%ld", gstarttns );
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
          gpsEndTimeNS += (INT8) gendt * 1000000000LL;
          ADD_PROCESS_PARAM( "int", "%ld", gendt );
        }
        break;

      case 'B':
        {
          long int gendtns = atol( optarg );
          if ( gendtns < 0 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS end time nanoseconds is negative\n",
                long_options[option_index].name );
            exit( 1 );
          }
          else if ( gendtns > 999999999 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS end time nanoseconds is greater than unity:\n" 
                "Must be <= 999999999:\n"
                "(%ld specified)\n", 
                long_options[option_index].name, gendtns );
            exit( 1 );
          }            
          gpsEndTimeNS += (INT8) gendtns;
          ADD_PROCESS_PARAM( "int", "%ld", gendtns );
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

          /* copy the first two characters to the ifo name */
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

      case 'C':
        {
          long int gstartt = atol( optarg );
          /* ignore a value of zero */
          if ( gstartt )
          {
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
          trigStartTimeNS = (INT8) gstartt * 1000000000LL;
          }
          ADD_PROCESS_PARAM( "int", "%ld", gstartt );
        }
        break;

      case 'D':
        {
          long int gendt = atol( optarg );
          /* ignore a value of zero */
          if ( gendt )
          {
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
            trigEndTimeNS = (INT8) gendt * 1000000000LL;
          }
          ADD_PROCESS_PARAM( "int", "%ld", gendt );
        }
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

      case 'f':
        ovrlap = (INT4) atoi( optarg );
        if ( ovrlap < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "data segment overlap must be positive: "
              "(%d specified)\n", 
              long_options[option_index].name, ovrlap );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", ovrlap );
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

      case 'y':
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

      case 'E':
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

      case 'P':
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

      case 'Q':
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
          fprintf( stderr, "invalid argument to --%s:\n"
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
        else if ( ! strcmp( "bad-mean", optarg ) )
        {
          specType = 0;
          badMeanPsd = 1;
        }
        else if ( ! strcmp( "gaussian", optarg ) )
        {
          specType = 2;
          fprintf( stderr,
              "WARNING: replacing psd with white gaussian spectrum\n" );
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown power spectrum type: "
              "%s (must be mean, median or gaussian)\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'k':
        invSpecTrunc = (INT4) atoi( optarg );
        if ( invSpecTrunc < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "inverse spectrum length must be positive or zero: "
              "(%d specified)\n", 
              long_options[option_index].name, invSpecTrunc );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", invSpecTrunc );
        break;

      case 'l':
        dynRangeExponent = (REAL4) atof( optarg );
        haveDynRange = 1;
        ADD_PROCESS_PARAM( "float", "%e", dynRangeExponent );
        break;

      case 'm':
        startTemplate = (INT4) atoi( optarg );
        if ( startTemplate < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "template bank start index must be positive: "
              "(%d specified)\n", 
              long_options[option_index].name, startTemplate );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", startTemplate );
        break;

      case 'M':
        minimalMatch = (REAL4) atof( optarg );
        if ( minimalMatch < 0 || minimalMatch > 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "minimal match must be in the range [0,1]: "          
              "(%e specified)\n", 
              long_options[option_index].name, minimalMatch );
        }
        /* process param added after bank is generated so that a */
        /* value in the bank looks like a command line option.   */
        break;

      case 'n':
        stopTemplate = (INT4) atoi( optarg );
        if ( stopTemplate < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "template bank stop index must be positive: "
              "(%d specified)\n", 
              long_options[option_index].name, stopTemplate );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", stopTemplate );
        break;

      case 'o':
        numChisqBins = (INT4) atoi( optarg );
        if ( numChisqBins < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "number of chisq veto bins must be positive: "
              "(%d specified)\n", 
              long_options[option_index].name, numChisqBins );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", numChisqBins );
        break;

      case 'p':
        /* create storage for the calibration frame cache name */
        optarg_len = strlen( optarg ) + 1;
        calCacheName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( calCacheName, optarg, optarg_len );
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
        else if ( ! strcmp( "GeneratePPN", optarg ) )
        {
          approximant = GeneratePPN;
        }
        else if ( ! strcmp( "PadeT1", optarg ) )
        {
          approximant = PadeT1;
        }
        else if ( ! strcmp( "EOB", optarg ) )
        {
          approximant = EOB;
        }
        else if ( ! strcmp( "TaylorF2", optarg ) )
        {
          approximant = TaylorF2;
        }
        else if ( ! strcmp( "BCV", optarg ) )
        {
          approximant = BCV;
        }
	else if ( ! strcmp( "BCVC", optarg ) )
        {
          approximant = BCV;
	  bcvConstraint = 1;
        }
        else if ( ! strcmp( "BCVSpin", optarg ) )
        {
          approximant = BCVSpin;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown order specified: "
              "%s (must be either TaylorF2 or BCV or BCVC or BCVSpin)\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        haveApprox = 1;
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'q':
        snrThresh = atof( optarg );
        if ( snrThresh < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "signal to noise threshold must be positive: "
              "(%f specified)\n", 
              long_options[option_index].name, snrThresh );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'r':
        chisqThresh = atof( optarg );
        if ( chisqThresh < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "chi squared threshold must be positive: "
              "(%f specified)\n", 
              long_options[option_index].name, chisqThresh );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'R':
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
        if ( highPassFreq <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "high pass filter frequency must be greater than 0 Hz: "
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
          fprintf( stderr, "invalid argument to --%s:\n"
              "high pass filter order must be greater than 0: "
              "(%d specified)\n",
              long_options[option_index].name, highPassOrder );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", highPassOrder );
        break;

      case 'T':
        highPassAtten = (REAL4) atof( optarg );
        if ( highPassAtten < 0.0 || highPassAtten > 1.0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "high pass attenuation must be in the range [0:1]: "
              "(%f specified)\n",
              long_options[option_index].name, highPassAtten );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", highPassAtten );
        break;

      case 'u':
        /* create storage for the input frame cache name */
        optarg_len = strlen( optarg ) + 1;
        frInCacheName = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( frInCacheName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'S':
        optarg_len = strlen( optarg ) + 1;
        frInType = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( frInType, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'v':
        /* create storage for the calibration frame cache name */
        optarg_len = strlen( optarg ) + 1;
        bankFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( bankFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'w':
        /* create storage for the injection file name */
        optarg_len = strlen( optarg ) + 1;
        injectionFile = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( injectionFile, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'x':
        padData = (INT4) atoi( optarg );
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

      case 'X':
        slideData.gpsSeconds = (INT4) atoi( optarg );
        ADD_PROCESS_PARAM( "int", "%d", slideData.gpsSeconds );
        break;

      case 'Y':
        slideData.gpsNanoSeconds = (INT4) atoi( optarg );
        ADD_PROCESS_PARAM( "int", "%d", slideData.gpsNanoSeconds );
        break;

      case 'K':
        bankSim = (INT4) atoi( optarg );
        if ( bankSim < 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "number of template bank simulations"
              "must be greater than 1: (%d specified)\n", 
              long_options[option_index].name, bankSim );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", bankSim );
        break;

      case 'L':
        if ( ! strcmp( "TaylorT2", optarg ) )
        {
          bankSimApproximant = TaylorT2;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown order specified: "
              "%s (must be TaylorT2)\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        haveBankSimApprox = 1;
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'J':
        if ( ! strcmp( "urandom", optarg ) )
        {
          randSeedType = urandom;
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        else
        {
          randSeedType = user;
          randomSeed = (INT4) atoi( optarg );
          ADD_PROCESS_PARAM( "int", "%d", randomSeed );
        }
        break;

      case 'U':
        bankMinMass = (REAL4) atof( optarg );
        if ( bankMinMass <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "miniumum component mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, bankMinMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", bankMinMass );
        break;

      case 'W':
        bankMaxMass = (REAL4) atof( optarg );
        if ( bankMaxMass <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "maxiumum component mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, bankMaxMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", bankMaxMass );
        break;

      case 'G':
        gaussVar = (REAL4) atof( optarg );
        if ( gaussVar < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "variance of gaussian noise must be >= 0: "
              "(%f specified)\n",
              long_options[option_index].name, gaussVar );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", gaussVar );
        gaussianNoise = 1;
        fprintf( stderr,
            "WARNING: replacing input data with white gaussian noise\n"
            "WARNING: replacing response function with unity\n" );
        break;

      case 'N':
        if ( LALSnprintf( ckptPath, FILENAME_MAX * sizeof(CHAR), 
              "%s", optarg ) < 0 )
        {
          fprintf( stderr, "invalid argument to --%s\n"
              "local path %s too long: string truncated\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );

      case 'O':
        if ( LALSnprintf( outputPath, FILENAME_MAX * sizeof(CHAR), 
              "%s", optarg ) < 0 )
        {
          fprintf( stderr, "invalid argument to --%s\n"
              "output path %s too long: string truncated\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );

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

      case 'I':
        /* create storaged for the ifo-tag */
        optarg_len = strlen( optarg ) + 1;
        ifoTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( ifoTag, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "LIGO/LSC Standalone Inspiral Search Engine\n" 
            "Duncan Brown <duncan@gravity.phys.uwm.edu>\n"
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n" );
        exit( 0 );
        break;

      case '*':								/*XXX*/
        if ( ! strcmp( "none", optarg ) )
        {
          clusterMethod = noClustering;
        }
	else if ( ! strcmp( "template", optarg ) )
        {
 	  clusterMethod = tmplt;
	} 
        else if ( ! strcmp( "window", optarg ) )
        {
 	  clusterMethod = window;
	} 
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown clustering method: "
              "%s (must be 'none', 'template' or 'window')\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }
	haveClusterMethod = 1;
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case '#':								/*XXX*/
        clusterWindow = (REAL4) atof( optarg );
        if ( clusterWindow <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "clustering time window is less than or equal to 0 secs: "
              "(%f secs specified)\n",
              long_options[option_index].name, clusterWindow );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", clusterWindow );
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

  /* enable output is stored in the first process param row */
  if ( enableOutput == 1 )
  {
    LALSnprintf( procparams.processParamsTable->program, 
        LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
    LALSnprintf( procparams.processParamsTable->param,
        LIGOMETA_PARAM_MAX, "--enable-output" );
    LALSnprintf( procparams.processParamsTable->type, 
        LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( procparams.processParamsTable->value, 
        LIGOMETA_TYPE_MAX, " " );
  }
  else if ( enableOutput == 0 )
  {
    LALSnprintf( procparams.processParamsTable->program, 
        LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
    LALSnprintf( procparams.processParamsTable->param,
        LIGOMETA_PARAM_MAX, "--disable-output" );
    LALSnprintf( procparams.processParamsTable->type, 
        LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( procparams.processParamsTable->value, 
        LIGOMETA_TYPE_MAX, " " );
  }
  else
  {
    fprintf( stderr, "--enable-output or --disable-output "
        "argument must be specified\n" );
    exit( 1 );
  } 

  /* check inject-overhead option */
  if ( injectOverhead )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--inject-overhead" );
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


  /*
   *
   * check validity of arguments
   *
   *
   */


  /* check validity of input data time */
  if ( ! gpsStartTimeNS )
  {
    fprintf( stderr, "--gps-start-time must be specified\n" );
    exit( 1 );
  }
  LAL_CALL( LALINT8toGPS( &status, &gpsStartTime, &gpsStartTimeNS ), 
      &status );
  if ( ! gpsEndTimeNS )
  {
    fprintf( stderr, "--gps-end-time must be specified\n" );
    exit( 1 );
  }
  LAL_CALL( LALINT8toGPS( &status, &gpsEndTime, &gpsEndTimeNS ), 
      &status );
  if ( gpsEndTimeNS <= gpsStartTimeNS )
  {
    fprintf( stderr, "invalid gps time range: "
        "start time: %d, end time %d\n",
        gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds );
    exit( 1 );
  }

  /* check trigger generation time is within input time */
  if ( trigStartTimeNS )
  {
    if ( trigStartTimeNS < gpsStartTimeNS )
    {
      fprintf( stderr, 
          "trigStartTimeNS = %lld\nis less than gpsStartTimeNS = %lld", 
          trigStartTimeNS, gpsStartTimeNS );
    }
  }
  if ( trigEndTimeNS )
  {
    if ( trigEndTimeNS > gpsEndTimeNS )
    {
      fprintf( stderr, 
          "trigEndTimeNS = %lld\nis greater than gpsEndTimeNS = %lld", 
          trigEndTimeNS, gpsEndTimeNS );
    }
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
  if ( ovrlap < 0 )
  {
    fprintf( stderr, "--segment-overlap must be specified\n" );
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
  inputDataLength = numPoints * numSegments - ( numSegments - 1 ) * ovrlap;
  {
    INT8 gpsChanIntervalNS = gpsEndTimeNS - gpsStartTimeNS;
    INT8 inputDataLengthNS = (INT8) inputDataLength * 1000000000LL / 
      (INT8) sampleRate;

    if ( inputDataLengthNS != gpsChanIntervalNS )
    {
      fprintf( stderr, "length of input data and data chunk do not match\n" );
      fprintf( stderr, "start time: %lld, end time %lld\n",
          gpsStartTimeNS / 1000000000LL, gpsEndTimeNS / 1000000000LL );
      fprintf( stderr, "gps channel time interval: %lld ns\n"
          "computed input data length: %lld ns\n", 
          gpsChanIntervalNS, inputDataLengthNS );
      exit( 1 );
    }
  }

  /* check filter parameters have been specified */
  if ( numChisqBins < 0 )
  {
    fprintf( stderr, "--chisq-bins must be specified\n" );
    exit( 1 );
  }
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
  if ( clusterWindow * sampleRate > numPoints )				
  {
    fprintf( stderr, "--cluster-window must be less than "
	"--segment-length\n" );
    exit( 1 );
  }
  if ( ! haveClusterMethod )
  {
    fprintf( stderr, "--cluster-method must be specified\n" );
    exit( 1 );
  }
  if ( clusterMethod == window && clusterWindow == -1 )
  {
    fprintf( stderr, "--cluster-window must be specified "
 	"if --clustering method 'window' chosen\n" );
    exit( 1 );
  }
  if ( clusterMethod != window && clusterWindow != -1 )
  {
    fprintf( stderr, "--cluster-window specified "
 	"but --clustering method 'window' not chosen\n" );
    exit( 1 );
  }
  if ( invSpecTrunc < 0 )
  {
    fprintf( stderr, "--inverse-spec-length must be specified\n" );
    exit( 1 );
  }
  else if ( invSpecTrunc * sampleRate > numPoints )
  {
    fprintf( stderr, "--inverse-spec-length must be less than "
        "--segment-length\n" );
    exit( 1 );
  }

  if ( ! haveDynRange )
  {
    fprintf( stderr, "--dynamic-range-exponent must be specified\n" );
    exit( 1 );
  }
  if ( ! haveApprox )
  {
    fprintf( stderr, "--approximant must be specified\n" );
    exit( 1 );
  }

  /* check that a channel has been requested and fill the ifo */
  if ( ! fqChanName )
  {
    fprintf( stderr, "--channel-name must be specified\n" );
    exit( 1 );
  }

  if ( ! bankSim )
  {
    /* check that the thresholds have been specified */
    if ( snrThresh < 0 )
    {
      fprintf( stderr, "--snr-threshold must be specified\n" );
      exit( 1 );
    }
    if ( chisqThresh < 0 )
    {
      fprintf( stderr, "--chisq-threshold must be specified\n" );
      exit( 1 );
    }
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

  /* check that a template bank has been specified */
  if ( ! bankFileName )
  {
    fprintf( stderr, "--bank-file must be specified\n" );
    exit( 1 );
  }

  /* check that a random seed for gaussian noise generation has been given */
  if ( gaussianNoise && randSeedType == unset )
  {
    fprintf( stderr, "--random-seed must be specified if "
        "--gaussian-noise is given\n" );
    exit( 1 );
  }

  /* check that if we are doing a bank sim that an approx has been given */
  if ( bankSim )
  {
    if ( ! haveBankSimApprox )
    {
      fprintf( stderr, "--sim-approximant must be specified if "
          "--bank-simulation is given\n" );
      exit( 1 );
    }
    if ( randSeedType == unset )
    {
      fprintf( stderr, "--random-seed must be specified if "
          "--bank-simulation is given\n" );
      exit( 1 );
    }
    if ( bankMinMass < 0 )
    {
      fprintf( stderr, "--sim-minimum-mass must be specified if "
          "--bank-simulation is given\n" );
      exit( 1 );
    }
    if ( bankMaxMass < 0 )
    {
      fprintf( stderr, "--sim-maximum-mass must be specified "
          "--bank-simulation is given\n" );
      exit( 1 );
    }
  }

  return 0;
}

#undef ADD_PROCESS_PARAM
