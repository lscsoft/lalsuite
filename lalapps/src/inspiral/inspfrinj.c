/*----------------------------------------------------------------------- 
 * 
 * File Name: inspfinj.c
 *
 * Author: Fairhurst, S. (based on inspiral.c by Brown, D.A.)
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
#include <lal/Calibration.h>
#include <lal/FrameCalibration.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/FindChirpSP.h>
#include <lal/Inject.h>

RCSID( "$Id$" );

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "inspiral"


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
this_summ_value->start_time = searchsumm.searchSummaryTable->in_start_time; \
this_summ_value->end_time = searchsumm.searchSummaryTable->in_end_time; \
this_summ_value->value = (REAL4) val; \
this_summ_value->intvalue = (INT4) intval; \
LALSnprintf( this_summ_value->name, LIGOMETA_SUMMVALUE_NAME_MAX, "%s", \
    sv_name ); \
LALSnprintf( this_summ_value->comment, LIGOMETA_SUMMVALUE_COMM_MAX, \
    "%s", sv_comment ); \

int arg_parse_check( int argc, char *argv[], MetadataTable procparams );


/*
 *
 * variables that control program behaviour
 *
 */


/* debugging */
extern int vrbflg;                      /* verbocity of lal function    */

/* input data parameters */
INT8  gpsStartTimeNS    = 0;            /* input data GPS start time ns */
LIGOTimeGPS gpsStartTime;               /* input data GPS start time    */
INT8  gpsEndTimeNS      = 0;            /* input data GPS end time ns   */
LIGOTimeGPS gpsEndTime;                 /* input data GPS end time      */
INT8  inputLengthNS     = 0;            /* input data length ns         */
INT4  numRespPoints     = -1;            /* num points for calc response */

CHAR  *fqChanName       = NULL;         /* name of data channel         */
CHAR  *frInCacheName    = NULL;         /* cache file containing frames */
CHAR   ifo[3];                          /* two character ifo code       */
CHAR  *channelName      = NULL;         /* channel string               */
CHAR   outfileName[FILENAME_MAX];       /* output file name             */

enum { undefined, real_4, real_8 } calData = undefined; /* cal data type */

/* data conditioning parameters */
INT4   sampleRate       = -1;           /* sample rate of filter data   */
INT4   frameLength      = -1;           /* length of output frames      */
INT4   injectSafety     = 0;            /* safety length in injections  */
UINT4  numFiles         = 0;            /* number of output files needed*/

CHAR  *calCacheName     = NULL;         /* location of calibration data */
CHAR  *injectionFile    = NULL;         /* name of file containing injs */

int   injectOverhead	= 0;		/* inject h+ into detector	*/
int   numInjections     = 0;
SimInspiralTable *injections = NULL;
SimInspiralTable    *thisInj = NULL;


/* output parameters */
CHAR  *userTag          = NULL;         /* string the user can tag with */
int    writeRawData     = 0;            /* write the raw data to frame  */
int    writeInjOnly     = 0;            /* write the inj data to frame  */
int    writeRawPlusInj  = 0;            /* write raw plus inj to frame  */
/* other command line args */
CHAR comment[LIGOMETA_COMMENT_MAX];     /* process param comment        */

int main( int argc, char *argv[] )
{
  /* lal function variables */
  LALStatus             status = blank_status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;

  /* frame input data */
  FrCache      *frInCache = NULL;
  FrCache      *calCache = NULL;
  FrStream     *frStream = NULL;
  FrChanIn      frChan;

  /* raw input data storage */
  REAL4TimeSeries               chan;
  REAL4TimeSeries               inj;

  /* structures for preconditioning */
  COMPLEX8FrequencySeries       injResp;

  /* output data */
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         searchsumm;
  MetadataTable         searchsummvars;
  MetadataTable         siminspiral;
  SearchSummvarsTable  *this_search_summvar;
  MetadataTable         summvalue;
  SummValueTable       *this_summ_value = NULL;
  ProcessParamsTable   *this_proc_param;
  LIGOLwXMLStream       results;

  /* counters and other variables */
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
  UINT4 k,n;
  CHAR  fname[FILENAME_MAX];
  UINT4 numPoints = 0;
  REAL8 tsLength;
  INT8  durationNS	= 0;
  CalibrationUpdateParams inj_calfacts;
  REAL4 inj_alpha = 0;
  REAL4 inj_alphabeta = 0;
  CHAR tmpChName[LALNameLength];
  REAL8 inputDeltaT;

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

  /* zero out the outfileName */
  memset( outfileName, 0, FILENAME_MAX * sizeof(CHAR) );

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


  /* the number of nodes for a standalone job is always 1 */
  searchsumm.searchSummaryTable->nnodes = 1;

  /*
   *
   * read in the input data channel
   *
   */

  memset( &chan, 0, sizeof(REAL4TimeSeries) );

  if ( frInCacheName )
  {
    /* set the params of the input data time series */
    chan.epoch = gpsStartTime;

    /* open a frame cache */
    LAL_CALL( LALFrCacheImport( &status, &frInCache, frInCacheName), 
	&status );
    LAL_CALL( LALFrCacheOpen( &status, &frStream, frInCache ), &status );

    /* set the mode of the frame stream to fail on gaps or time errors */
    frStream->mode = LAL_FR_VERBOSE_MODE;

    /* seek to required epoch and set chan name */
    LAL_CALL( LALFrSeek( &status, &(chan.epoch), frStream ), &status );
    frChan.name = fqChanName;

    /* determine the sample rate of the raw data */
    LAL_CALL( LALFrGetREAL4TimeSeries( &status, &chan, &frChan, frStream ),
	&status );

    /* store the input sample rate */
    this_search_summvar = searchsummvars.searchSummvarsTable = 
      (SearchSummvarsTable *) LALCalloc( 1, sizeof(SearchSummvarsTable) );
    LALSnprintf( this_search_summvar->name, LIGOMETA_NAME_MAX * sizeof(CHAR),
	"raw data sample rate" );
    this_search_summvar->value = inputDeltaT = chan.deltaT;

    /* determine the number of points to get and create storage for the data */
    numPoints = (UINT4) floor( ((REAL8) inputLengthNS) / (chan.deltaT * 1.0e9) 
	+ 0.5 );
    LAL_CALL( LALSCreateVector( &status, &(chan.data), numPoints ), 
	&status );

    if ( vrbflg ) fprintf( stdout, "input channel %s has sample interval "
	"(deltaT) = %e\nreading %d points from frame stream\n", fqChanName, 
	chan.deltaT, numPoints );

    /* read the data channel time series from frames */
    LAL_CALL( LALFrGetREAL4TimeSeries( &status, &chan, &frChan, frStream ),
	&status );
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
    if ( frInCacheName ) LAL_CALL( LALDestroyFrCache( &status, &frInCache ), 
	&status );

    if ( vrbflg ) fprintf( stdout, "read channel %s from frame stream\n"
	"got %d points with deltaT %e\nstarting at GPS time %d sec %d ns\n", 
	chan.name, chan.data->length, chan.deltaT, 
	chan.epoch.gpsSeconds, chan.epoch.gpsNanoSeconds );
  }


  /* 
   *
   * Create zeros on top of which to do the injections 
   *
   */

  memset( &inj, 0, sizeof(REAL4TimeSeries) );
  if ( frInCacheName )
  {
    inj.deltaT = chan.deltaT;
    inj.epoch = chan.epoch;
    inj.sampleUnits = chan.sampleUnits;
    strcpy( inj.name, chan.name );
  }  
  else
  {
    inj.deltaT = 1.0/ sampleRate;
    inj.epoch = gpsStartTime;
    inj.sampleUnits = lalADCCountUnit;
    LALSnprintf( inj.name, LIGOMETA_CHANNEL_MAX * sizeof(CHAR), "%s:STRAIN", 
	ifo );
    searchsumm.searchSummaryTable->in_start_time = gpsStartTime;
    searchsumm.searchSummaryTable->in_end_time = gpsEndTime;
    numPoints = (UINT4) floor( ((REAL8) inputLengthNS) / (inj.deltaT * 1.0e9) 
	+ 0.5 );
  }
  LAL_CALL( LALSCreateVector( &status, &(inj.data), numPoints ), 
      &status );
  memset( inj.data->data, 0, numPoints * sizeof(REAL4) );


  /*
   *
   * inject signals into the raw, unresampled data
   *
   */

  if ( injectionFile )
  {
    /* read in the injection data from XML */
    numInjections = SimInspiralTableFromLIGOLw( &injections, injectionFile,
	gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds + injectSafety );

    if ( numInjections < 0 )
    {
      fprintf( stderr, "error: cannot read injection file" );
      exit( 1 );
    }
    else if ( numInjections )
    {
      searchsumm.searchSummaryTable->nevents = numInjections;
      memset( &injResp, 0, sizeof(COMPLEX8FrequencySeries) );
      LAL_CALL( LALCCreateVector( &status, &(injResp.data), 
	    numRespPoints / 2 + 1 ), &status );
      injResp.epoch = inj.epoch ;
      injResp.deltaF = 1.0 / ( numRespPoints * inj.deltaT );
      strcpy( injResp.name, inj.name );

      if ( calCache )
      {
	/* generate the response function for the current time */
	if ( vrbflg ) fprintf( stdout, 
	    "generating response function at time %d sec %d ns\n"
	    "length = %d points, deltaF = %e Hz\n",
	    injResp.epoch.gpsSeconds, injResp.epoch.gpsNanoSeconds,
	    injResp.data->length, injResp.deltaF );
	injResp.sampleUnits = strainPerCount;

	/* initialize the inj_calfacts */
	memset( &inj_calfacts, 0, sizeof(CalibrationUpdateParams) );
	inj_calfacts.ifo = ifo;
	durationNS = gpsEndTimeNS - gpsStartTimeNS;
	LAL_CALL( LALINT8toGPS( &status, &(inj_calfacts.duration), 
	      &durationNS ), &status );

        LAL_CALL( LALFrCacheImport( &status, &calCache, calCacheName ), 
            &status );
	LAL_CALL( LALExtractFrameResponse( &status, &injResp, calCache, 
	      &inj_calfacts ), &status );
        LAL_CALL( LALDestroyFrCache( &status, &calCache ), &status );
	inj_alpha = (REAL4) inj_calfacts.alpha.re;
	inj_alphabeta = (REAL4) inj_calfacts.alphabeta.re;
	if ( vrbflg ) fprintf( stdout, 
	    "for injections, alpha = %f and alphabeta = %f\n",
	    inj_alpha, inj_alphabeta);
      }
      else 
      {
	/* generate a unity response function for h(t) */
	if ( vrbflg ) fprintf( stdout, "setting response to unity... " );
	injResp.sampleUnits = strainPerCount;
	for ( k = 0; k < injResp.data->length; ++k )
	{
	  injResp.data->data[k].re = 1.0;
	  injResp.data->data[k].im = 0;
	}
	if ( vrbflg ) fprintf( stdout, "done.\n" );

      }


      /* inject the signals, preserving the channel name (Tev mangles it) */
      LALSnprintf( tmpChName, LALNameLength * sizeof(CHAR), "%s", inj.name );

      /* if injectOverhead option, then set inj.name to "ZENITH".  
       * This causes no detector site to be found in the injection code so
       * that the injection is done directly overhead (i.e. with a response 
       * function of F+ = 1; Fx = 0) */
      if ( injectOverhead )
      {
	LALSnprintf( inj.name, LALNameLength * sizeof(CHAR), "ZENITH" );
      }

      LAL_CALL( LALFindChirpInjectSignals( &status, &inj, injections, 
	    &injResp ), &status );
      LALSnprintf( inj.name,  LALNameLength * sizeof(CHAR), "%s", tmpChName );

      if ( vrbflg ) fprintf( stdout, "injected %d signals from %s into %s\n", 
	  numInjections, injectionFile, inj.name );

      LAL_CALL( LALCDestroyVector( &status, &(injResp.data) ), &status );
    }
    else
    {
      if ( vrbflg ) fprintf( stdout, "no injections in this data\n" );
    }
  }

  /* 
   *
   * Create data segments of the desired length and save them as frames
   *
   */

  if ( writeRawData || writeRawPlusInj || writeInjOnly )
  {
    /* frame output data */
    struct FrFile *frOutFile  = NULL;
    struct FrameH *outFrame   = NULL;
    REAL4TimeSeries output;
    UINT4 length;


    memset( &output, 0, sizeof(REAL4TimeSeries) );
    output.deltaT = inj.deltaT;
    output.sampleUnits = inj.sampleUnits;
    output.data = (REAL4Vector *) LALCalloc( 1, sizeof(REAL4Vector) );
    length = numPoints / numFiles;
    output.data->length = length;

    for ( n = 0; n < numFiles; ++n )
    {
      outFrame = NULL;
      output.epoch.gpsSeconds  = gpsStartTime.gpsSeconds + n * frameLength;

      /* write the injection channel to frame */
      if ( writeInjOnly  ) 
      {
	strcpy( output.name, inj.name );
	output.data->data = inj.data->data + n * length;
	outFrame = fr_add_proc_REAL4TimeSeries( outFrame, &output, "ct", 
	    "INSP_INJ_ONLY" );
      }
      /* write the raw data to frame */
      if ( writeRawData || writeRawPlusInj ) 
      {
	strcpy( output.name, chan.name );
	output.data->data = chan.data->data + n * length;

	if ( writeRawData )
	{
	  outFrame = fr_add_proc_REAL4TimeSeries( outFrame, &output, "ct", 
	      NULL );
	}
	/* perform injections into this file's data only, preserve name*/
	LAL_CALL( LALSSInjectTimeSeries( &status, &output, &inj ), &status );
	strcpy( output.name, inj.name );

	if ( writeRawPlusInj )
	{
	  outFrame = fr_add_proc_REAL4TimeSeries( outFrame, &output, "ct", 
	      "RAW_PLUS_INSP_INJ" );
	}
      }
      if( !outfileName[0] )
      {
	/* output name not specified, set to IFO-INSPFRINJ-EPOCH-LENGTH.gwf */
	LALSnprintf( outfileName, FILENAME_MAX * sizeof(CHAR), 
	    "%s-INSPFRINJ", ifo );
      }

      if( userTag )
      {
	LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), 
	    "%s_%s-%d-%d.gwf", outfileName, userTag, output.epoch.gpsSeconds, 
	    frameLength );
      }
      else
      {
	LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), 
	    "%s-%d-%d.gwf", outfileName, output.epoch.gpsSeconds, 
	    frameLength );

      }

      if ( vrbflg ) fprintf( stdout, "writing frame data to %s... ", fname );
      frOutFile = FrFileONew( fname, 3 );
      FrameWrite( outFrame, frOutFile );
      FrFileOEnd( frOutFile );
      if ( vrbflg ) fprintf( stdout, "done\n" );
    }

    LALFree ( output.data );
  }

  /* free the data storage */
  if ( vrbflg ) fprintf( stdout, "freeing memory\n" );
  if ( chan.data )
  {
    LAL_CALL( LALSDestroyVector( &status, &(chan.data) ), &status );
  }
  if ( inj.data )
  {
    LAL_CALL( LALSDestroyVector( &status, &(inj.data) ), &status );
  }


  /* open the output xml file */
  memset( &results, 0, sizeof(LIGOLwXMLStream) );
  if( userTag )
  {
    LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), 
	"%s_%s-%d-%d.xml", outfileName, userTag, gpsStartTime.gpsSeconds, 
	gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else
  {
    LALSnprintf( fname, FILENAME_MAX * sizeof(CHAR), 
	"%s-%d-%d.xml", outfileName, gpsStartTime.gpsSeconds, 
	gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );

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
  if ( searchsummvars.searchSummvarsTable )
  {
    if ( vrbflg ) fprintf( stdout, "  search_summvars table...\n" );
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
  }
  /* write the summ_value table with the calibration data used for injections */
  if ( frInCacheName && injectionFile)
  {
    ADD_SUMM_VALUE( "calibration alpha", "injection", inj_alpha, 0 );
    ADD_SUMM_VALUE( "calibration alphabeta", "injection", inj_alphabeta, 0 );

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
  }

  /* free the search summary table */
  free( searchsumm.searchSummaryTable );

  /* write the sim_inspiral table */
  if ( injections )
  {
    if ( vrbflg ) fprintf( stdout, "sim_inspiral... " );
    siminspiral.simInspiralTable = injections;
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, 
	  sim_inspiral_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, siminspiral, 
	  sim_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &results ), &status );

    /* free the temporary memory containing the events */
    while ( injections )
    {
      thisInj = injections;
      injections = injections->next;
      LALFree( thisInj );
    }
  }

  /* close the output xml file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &results ), &status );
  if ( vrbflg ) fprintf( stdout, "done. XML file closed\n" );

  /* free the rest of the memory, check for memory leaks and exit */
  if ( injectionFile ) free ( injectionFile ); 
  if ( calCacheName ) free( calCacheName );
  if ( frInCacheName ) free( frInCacheName );
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
  "lalapps_inspfrinj [options]\n\n"\
"  --help                    display this message\n"\
"  --verbose                 print progress information\n"\
"  --version                 print version information and exit\n"\
"  --debug-level LEVEL       set the LAL debug level to LEVEL\n"\
"  --user-tag STRING         set the process_params usertag to STRING\n"\
"  --comment STRING          set the process table comment to STRING\n"\
"\n"\
"  --gps-start-time SEC      GPS second of data start time\n"\
"  --gps-start-time-ns NS    GPS nanosecond of data start time\n"\
"  --gps-end-time SEC        GPS second of data end time\n"\
"  --gps-end-time-ns NS      GPS nanosecond of data end time\n"\
"\n"\
"  --frame-cache             obtain frame data from LAL frame cache FILE\n"\
"  --calibration-cache FILE  obtain calibration from LAL frame cache FILE\n"\
"  --calibrated-data TYPE    calibrated data of TYPE real_4 or real_8\n"\
"  --num-resp-points N       num points to determine response function (4194304)\n"\
"  --channel-name CHAN       read data from interferometer channel CHAN\n"\
"\n"\
"  --injection-file FILE     inject simulated inspiral signals from FILE\n"\
"  --inject-overhead         inject signals from overhead detector\n"\
"  --inject-safety SEC       inject signals ending up to SEC after gps end time\n"\
"\n"\
"  --write-raw-data          write out the raw frame files\n"\
"  --write-inj-only          write out frames containing only injections\n"\
"  --write-raw-plus-inj      write out frames containing raw data with inj\n"\
"\n"\
"  --output-frame-length SEC write out data in frames of length SEC\n"\
"  --output-file-name OUTPUT set output file names to OUTPUT-GPSTIME-LENGTH.gwf\n"\
"                   if not set, default to IFO-INSPFRINJ-GPSTIME-LENGTH.gwf\n"\
"\n"\
"  --ifo  IFO                specify the IFO (only if not reading frames)\n"\
"  --sample-rate             data sample rate (only if not reading frames)\n"\
"\n"

int arg_parse_check( int argc, char *argv[], MetadataTable procparams )
{
  /* getopt arguments */
  struct option long_options[] =
  {
    /* these options set a flag */
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"write-raw-data",          no_argument,       &writeRawData,     1 },
    {"write-inj-only",          no_argument,       &writeInjOnly,     1 },
    {"write-raw-plus-inj",      no_argument,       &writeRawPlusInj,  1 },
    {"inject-overhead",		no_argument,	   &injectOverhead,   1 },
    /* these options don't set a flag */
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-start-time-ns",       required_argument, 0,                'A'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"gps-end-time-ns",         required_argument, 0,                'B'},
    {"channel-name",            required_argument, 0,                'c'},
    {"output-frame-length",     required_argument, 0,                'd'},
    {"output-file-name",        required_argument, 0,                'f'},
    {"help",                    no_argument,       0,                'h'},
    {"dynamic-range-exponent",  required_argument, 0,                'l'},
    {"calibrated-data",         required_argument, 0,                'y'},
    {"calibration-cache",       required_argument, 0,                'p'},
    {"num-resp-points",         required_argument, 0,                'N'},
    {"sample-rate",             required_argument, 0,                'r'},
    {"ifo",                     required_argument, 0,                'i'},
    {"comment",                 required_argument, 0,                's'},
    {"frame-cache",             required_argument, 0,                'u'},
    {"injection-file",          required_argument, 0,                'w'},
    {"inject-safety",           required_argument, 0,                'S'},
    {"debug-level",             required_argument, 0,                'z'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {"version",                 no_argument,       0,                'V'},
    {0, 0, 0, 0}
  };
  int c;
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
	"A:B:N:S:V:Z:"
	"a:b:c:d:f:hi:l:p:r:s:u:w:y:z:",
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
	/* set gps start seconds */
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
	/* set gps start nanoseconds */
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
	/* set gps end seconds */
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
	/* set gps end nanoseconds */
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
	frameLength = (INT4) atoi( optarg );
	if ( frameLength < 1 )
	{
	  fprintf( stderr, "invalid argument to --%s:\n"
	      "length of frame must be a positive integer: "
	      "(%d specified) \n", 
	      long_options[option_index].name, frameLength );
	  exit( 1 );
	}
	ADD_PROCESS_PARAM( "int", "%d", frameLength );
	break;

      case 'f':
	if ( LALSnprintf( outfileName, FILENAME_MAX * sizeof(CHAR), 
	      "%s", optarg ) < 0 )
	{
	  fprintf( stderr, "invalid argument to --%s\n"
	      "outfile name %s too long: string truncated\n",
	      long_options[option_index].name, optarg );
	  exit( 1 );
	}
	ADD_PROCESS_PARAM( "string", "%s", optarg );
	break;

      case 'h':
	fprintf( stdout, USAGE );
	exit( 0 );
	break;

      case 'p':
	/* create storage for the calibration frame cache name */
	optarg_len = strlen( optarg ) + 1;
	calCacheName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
	memcpy( calCacheName, optarg, optarg_len );
	ADD_PROCESS_PARAM( "string", "%s", optarg );
	break;


      case 'N':
	/* store the number of points used in computing the response */
	numRespPoints = (INT4) atoi( optarg );
	if ( numRespPoints < 2 || numRespPoints % 2 )
	{
	  fprintf( stderr, "invalid argument to --%s:\n"
	      "number of points must be an even positive integer,\n"
	      "(%d specified) \n", 
	      long_options[option_index].name, numRespPoints );
	  exit( 1 );
	}
	ADD_PROCESS_PARAM( "int", "%d", numRespPoints );
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
	    fprintf( stderr, "Sorry, code not currently set up to\n"
		"run on real_8 data\n" );
	    exit( 1 );
	  }
	  else
	  {
	    fprintf( stderr, "invalid argument to --%s:\n"
		"unknown data type specified;\n"
		"%s (must be one of: real_4, real_8)\n",
		long_options[option_index].name, optarg);
	    exit( 1 );
	  }
	  ADD_PROCESS_PARAM( "string", "%s", optarg );
	}
	break;

      case 'r':
	sampleRate = (INT4) atoi( optarg );
	if ( sampleRate < 1 )
	{
	  fprintf( stderr, "invalid argument to --%s:\n"
	      "sample rate must be a positive integer: "
	      "(%d specified) \n", 
	      long_options[option_index].name, sampleRate );
	  exit( 1 );
	}
	ADD_PROCESS_PARAM( "int", "%d", sampleRate );
	break;

      case 'i':
	{
	  /* create storage for the ifo name and copy it */
	  memset( ifo, 0, sizeof(ifo) );
	  memcpy( ifo, optarg, sizeof(ifo) - 1 );
	  ADD_PROCESS_PARAM( "string", "%s", optarg );
	}
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

      case 'S':
	injectSafety = (INT4) atoi( optarg );
	if ( injectSafety < 1 )
	{
	  fprintf( stderr, "invalid argument to --%s:\n"
	      "injection safety must be a positive integer: "
	      "(%d specified) \n", 
	      long_options[option_index].name, frameLength );
	  exit( 1 );
	}
	ADD_PROCESS_PARAM( "int", "%d", injectSafety );
	break;

      case 'u':
	/* create storage for the input frame cache name */
	optarg_len = strlen( optarg ) + 1;
	frInCacheName = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
	memcpy( frInCacheName, optarg, optarg_len );
	ADD_PROCESS_PARAM( "string", "%s", optarg );
	break;

      case 'w':
	/* create storage for the injection file name */
	optarg_len = strlen( optarg ) + 1;
	injectionFile = (CHAR *) calloc( optarg_len, sizeof(CHAR));
	memcpy( injectionFile, optarg, optarg_len );
	ADD_PROCESS_PARAM( "string", "%s", optarg );
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

      case 'V':
	/* print version information and exit */
	fprintf( stdout, "LIGO/LSC Inspiral Injection Program\n" 
	    "Steve Fairhurst <sfairhur@gravity.phys.uwm.edu>\n"
	    "CVS Version: " CVS_ID_STRING "\n"
	    "CVS Tag: " CVS_NAME_STRING "\n" );
	exit( 0 );
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

  /* check flags and store in process_params */

  if ( writeRawData )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
	"%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
	"--write-raw-data" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  if ( writeInjOnly )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
	"%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
	"--write-inj-only" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }
  if ( writeRawPlusInj )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
	"%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
	"--write-raw-plus-inj" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
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

  /*
   *
   * check validity of arguments
   *
   */


  /* check validity of input data time */
  
  /* start time specified */
  if ( ! gpsStartTimeNS )
  {
    fprintf( stderr, "--gps-start-time must be specified\n" );
    exit( 1 );
  }
  LAL_CALL( LALINT8toGPS( &status, &gpsStartTime, &gpsStartTimeNS ), 
      &status );

  /* end time specified */
  if ( ! gpsEndTimeNS )
  {
    fprintf( stderr, "--gps-end-time must be specified\n" );
    exit( 1 );
  }
  LAL_CALL( LALINT8toGPS( &status, &gpsEndTime, &gpsEndTimeNS ), 
      &status );

  /* end after start */
  if ( gpsEndTimeNS <= gpsStartTimeNS )
  {
    fprintf( stderr, "invalid gps time range: "
	"start time: %d, end time %d\n",
	gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds );
    exit( 1 );
  }
  
  /* calculate the length of the data in NS */
  inputLengthNS = ( gpsEndTimeNS - gpsStartTimeNS );

  /* check that the output frame length has been specified */
  if ( frameLength == -1 )
  {
    fprintf( stderr, "--output-frame-length must be specified\n" );
    exit( 1 );
  }
  /* and it divides the total time exactly */
  if ( inputLengthNS % ( (INT8) frameLength * 1000000000LL ) )
  {
    fprintf(stderr, "data length %d must be a multiple of frame length %d",
	(INT4) (inputLengthNS / 1000000000LL), frameLength);
    exit( 1 );
  }
  else
  {
    numFiles = (INT4) (inputLengthNS / 1000000000LL) / frameLength;
  }

  /* if a frame cache has been specified, check we have everything else
     which is necessary */
  if ( frInCacheName )
  {
    /* check that a channel has been requested */
    if (! fqChanName )
    {
      fprintf( stderr, "--channel-name must be specified\n" );
      exit( 1 );
    }
  }
  else
  {
    /* check sample rate has been given */
    if ( sampleRate < 0 )
    {
      fprintf( stderr, "If --frame-cache not specified,\n"
	  "--sample-rate must be specified\n" );
      exit( 1 );
    }
    /* check that neither write-raw-data or write-raw-plus-inj given */
    if ( writeRawData || writeRawPlusInj )
    {
      fprintf( stderr, "Neither --write-raw-data nor --write-raw-plus-inj\n"
	  "can be specified when --frame-cache not given\n");
    }
  }

  /* check that have one of: calibrated-data or a calibration-cache  */
  if ( ! calCacheName && ! calData )
  {
    fprintf( stderr, "Either --calibration-cache must be specified,\n" 
	"or must run on --calibrated-data.\n");
    exit( 1 );
  }
  else if ( calCacheName && calData )
  {
    fprintf( stderr, 
	"Only one of --calibration-cache and --calibrated-data\n"
	"should be specified\n.");
    exit( 1 );
  }
    

  /* check that we have then number of points to determine response fn */
  if ( numRespPoints < 0 )
  {
    fprintf( stderr, "--num-resp-points not specified\n"
	"This gives the number of points used to obtain response,\n"
	"Response has numRespPoints/2 + 1 points,\n"
	"Frequency resolution of 1/(numRespPoints * delta T).\n"
	"Setting it to a default value of 256 * 16384 = 4194304\n");
    numRespPoints = 4194304;
  }


  return 0;
}

#undef ADD_PROCESS_PARAM
