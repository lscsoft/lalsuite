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

#include "tmpltbank.h"

RCSID( "$Id$" );
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "tmpltbank"


/*
 *
 * variables that control program behaviour
 *
 */


/* debugging */
extern int vrbflg;                      /* verbocity of lal function    */

/* parameters used to generate calibrated power spectrum */
LIGOTimeGPS gpsStartTime = { 0, 0 };    /* input data GPS start time    */
LIGOTimeGPS gpsEndTime = { 0, 0};       /* input data GPS end time      */
INT4  padData = 0;                      /* saftety margin on input data */
CHAR  *fqChanName       = NULL;         /* name of data channel         */
CHAR  *frInCacheName    = NULL;         /* cache file containing frames */
INT4  numPoints         = -1;           /* points in a segment          */
INT4  numSegments       = -1;           /* number of segments           */
CHAR  ifo[3];                           /* two character ifo code       */
CHAR *channelName = NULL;               /* channel string               */
INT4  inputDataLength = 0;              /* number of points in input    */
INT4   resampFiltType   = -1;           /* low pass filter used for res */
INT4   sampleRate       = -1;           /* sample rate of filter data   */
INT4   highPass         = -1;           /* enable high pass on raw data */
REAL4  highPassFreq     = 0;            /* high pass frequency          */
REAL4  fLow             = -1;           /* low frequency cutoff         */
INT4   specType         = -1;           /* use median or mean psd       */
CHAR  *calCacheName     = NULL;         /* location of calibration data */

/* template bank generation parameters */
REAL4   minMass         = -1;           /* minimum component mass       */
REAL4   maxMass         = -1;           /* maximum component mass       */
REAL4   minMatch        = -1;           /* minimun requested match      */
REAL4   fUpper          = -1;           /* upper frequency cutoff       */
Order   order;                          /* post-Newtonian order         */
Approximant approximant;                /* approximation method         */
CoordinateSpace space;                  /* coordinate space used        */

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
  FrStream     *frStream = NULL;
  FrChanIn      frChan;
  
  /* frame output data */
  struct FrFile *frOutFile = NULL;
  struct FrameH *outFrame  = NULL;
  
  /* raw input data storage */
  REAL4TimeSeries               chan;
  REAL4FrequencySeries          spec;
  COMPLEX8FrequencySeries       resp;

  /* structures for preconditioning */
  ResampleTSParams              resampleParams;
  LALWindowParams               wpars;
  AverageSpectrumParams         avgSpecParams;

  /* templates */
  InspiralCoarseBankIn          bankIn;
  InspiralTemplateList         *coarseList = NULL;
  SnglInspiralTable            *tmplt  = NULL;
  INT4                          numCoarse = 0;

  /* output data */
  MetadataTable         templateBank;
  MetadataTable         proctable;
  MetadataTable         procparams;
  ProcessParamsTable   *this_proc_param;
  LIGOLwXMLStream       results;

  /* counters and other variables */
  INT4 i;
  UINT4 cut, k;
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
  CHAR  fname[256];
  LALUnitPair pair;
  REAL8 respRe, respIm;
  REAL8 shf;
  REAL8 inputLengthNS;
  UINT4 numInputPoints;
  const REAL8 epsilon = 1.0e-8;
  UINT4 resampleChan = 0;


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

  /* call the argument parse and check function */
  arg_parse_check( argc, argv, procparams );

  /* can use LALMalloc() / LALCalloc() from here */

  /* fill the comment, if a user has specified on, or leave it blank */
  if ( ! *comment )
  {
    LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
  } else {
    LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
  }

  /* make sure the pointer to the first template is null */
  templateBank.snglInspiralTable = NULL;


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
  chan.epoch = gpsStartTime;
  chan.epoch.gpsSeconds -= padData; /* subtract pad seconds from start */

  /* open a frame cache or files, seek to required epoch and set chan name */
  if ( frInCacheName )
  {
    LAL_CALL( LALFrCacheImport( &status, &frInCache, frInCacheName), &status );
    LAL_CALL( LALFrCacheOpen( &status, &frStream, frInCache ), &status );
  }
  else
  {
    LAL_CALL( LALFrOpen( &status, &frStream, NULL, "*.gwf" ), &status );
  }
  LAL_CALL( LALFrSeek( &status, &(chan.epoch), frStream ), &status );
  frChan.name = fqChanName;

  /* determine the sample rate of the raw data and allocate enough memory */
  LAL_CALL( LALFrGetREAL4TimeSeries( &status, &chan, &frChan, frStream ),
      &status );

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
  LAL_CALL( LALSCreateVector( &status, &(chan.data), numInputPoints ), 
      &status );

  if ( vrbflg ) fprintf( stdout, "input channel %s has sample interval "
      "(deltaT) = %e\nreading %d points from frame stream\n", fqChanName, 
      chan.deltaT / 1.0e9, numInputPoints );

  /* read the data channel time series from frames */
  LAL_CALL( LALFrGetREAL4TimeSeries( &status, &chan, &frChan, frStream ),
      &status );
  memcpy( &(chan.sampleUnits), &lalADCCountUnit, sizeof(LALUnit) );

  /* close the frame file stream and destroy the cache */
  LAL_CALL( LALFrClose( &status, &frStream ), &status );
  if ( frInCacheName ) LAL_CALL( LALDestroyFrCache( &status, &frInCache ), 
      &status );

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

  /* high pass the data to get a prevent bleeding of low frequences in psd */
  if ( highPass )
  {
    PassBandParamStruc highpassParam;
    highpassParam.nMax = 4;
    highpassParam.f1 = highPassFreq;
    highpassParam.f2 = -1.0;
    highpassParam.a1 = 0.1;
    highpassParam.a2 = -1.0;

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

  /* generate the response function for the current time */
  if ( vrbflg ) fprintf( stdout, "generating response at time %d sec %d ns\n"
      "response parameters f0 = %e, deltaF = %e, length = %d\n",
      resp.epoch.gpsSeconds, resp.epoch.gpsNanoSeconds,
      resp.f0, resp.deltaF, resp.data->length );
  LAL_CALL( LALExtractFrameResponse( &status, &resp, calCacheName, ifo ),
      &status );

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

  /* bank generation parameters */
  bankIn.massRange     = MinMaxComponentMass;
  bankIn.mMin          = (REAL8) minMass;
  bankIn.mMax          = (REAL8) maxMass;
  bankIn.MMax          = bankIn.mMax * 2.0;
  bankIn.mmCoarse      = (REAL8) minMatch;
  bankIn.mmFine        = 0.99; /* doesn't matter since no fine bank yet */
  bankIn.fLower        = (REAL8) fLow;
  bankIn.fUpper        = (REAL8) fUpper;
  bankIn.iflso         = 0; /* currently not implemented */
  bankIn.tSampling     = (REAL8) sampleRate;
  bankIn.order         = order;
  bankIn.approximant   = approximant;
  bankIn.space         = space;
  bankIn.etamin        = bankIn.mMin * ( bankIn.MMax - bankIn.mMin) /
    ( bankIn.MMax * bankIn.MMax );

  /* generate the template bank */
  if ( vrbflg )
  {
    fprintf( stdout, "generating template bank parameters... " );
    fflush( stdout );
  }
  LAL_CALL( LALInspiralCreateCoarseBank( &status, &coarseList, &numCoarse, 
        bankIn ), &status );
  if ( vrbflg )
  {
    fprintf( stdout, "done\n" );
    fflush( stdout );
  }

  /* convert the templates to sngl_inspiral structures for writing to XML */
  if ( numCoarse )
  {
    tmplt = templateBank.snglInspiralTable = (SnglInspiralTable *)
      LALCalloc( 1, sizeof(SnglInspiralTable) );
    LALSnprintf( tmplt->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR), ifo );
    LALSnprintf( tmplt->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR), 
        "tmpltbank" );
    LALSnprintf( tmplt->channel, LIGOMETA_CHANNEL_MAX * sizeof(CHAR),
        channelName );
    tmplt->mass1  = (REAL4) coarseList[0].params.mass1;
    tmplt->mass2  = (REAL4) coarseList[0].params.mass2;
    tmplt->mchirp = (REAL4) coarseList[0].params.chirpMass;
    tmplt->eta    = (REAL4) coarseList[0].params.eta;
    tmplt->tau0   = (REAL4) coarseList[0].params.t0;
    tmplt->tau2   = (REAL4) coarseList[0].params.t2;
    tmplt->tau3   = (REAL4) coarseList[0].params.t3;
    tmplt->tau4   = (REAL4) coarseList[0].params.t4;
    tmplt->tau5   = (REAL4) coarseList[0].params.t5;
    tmplt->ttotal = (REAL4) coarseList[0].params.tC;
    for ( i = 1; i < numCoarse; ++i )
    {
      tmplt = tmplt->next = (SnglInspiralTable *)
        LALCalloc( 1, sizeof(SnglInspiralTable) );
      LALSnprintf( tmplt->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR), ifo );
      LALSnprintf( tmplt->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR), 
          "tmpltbank" );
      LALSnprintf( tmplt->channel, LIGOMETA_CHANNEL_MAX * sizeof(CHAR),
          channelName );
      tmplt->mass1  = (REAL4) coarseList[i].params.mass1;
      tmplt->mass2  = (REAL4) coarseList[i].params.mass2;
      tmplt->mchirp = (REAL4) coarseList[i].params.chirpMass;
      tmplt->eta    = (REAL4) coarseList[i].params.eta;
      tmplt->tau0   = (REAL4) coarseList[i].params.t0;
      tmplt->tau2   = (REAL4) coarseList[i].params.t2;
      tmplt->tau3   = (REAL4) coarseList[i].params.t3;
      tmplt->tau4   = (REAL4) coarseList[i].params.t4;
      tmplt->tau5   = (REAL4) coarseList[i].params.t5;
      tmplt->ttotal = (REAL4) coarseList[i].params.tC;
    }
  }


  /*
   *
   * free the data storage
   *
   */


  LALFree( coarseList );
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
  free( calCacheName );
  free( frInCacheName );
  free( channelName );
  free( fqChanName );
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
"  --help                       display this message\n"\
"  --verbose                    print progress information\n"\
"  --debug-level LEVEL          set the LAL debug level to LEVEL\n"\
"  --user-tag STRING            set the process_params usertag to STRING\n"\
"  --comment STRING             set the process table comment to STRING\n"\
"\n"\
"  --gps-start-time SEC         GPS second of data start time\n"\
"  --gps-end-time SEC           GPS second of data end time\n"\
"  --pad-data T                 pad the data start and end time by T seconds\n"\
"\n"\
"  --frame-cache                obtain frame data from LAL frame cache FILE\n"\
"  --calibration-cache FILE     obtain calibration from LAL frame cache FILE\n"\
"  --channel-name CHAN          read data from interferometer channel CHAN\n"\
"\n"\
"  --sample-rate F              filter data at F Hz, downsampling if necessary\n"\
"  --resample-filter TYPE       set resample filter to TYPE [ldas|butterworth]\n"\
"\n"\
"  --disable-high-pass          turn off the IIR highpass filter\n"\
"  --enable-high-pass F         high pass data above F Hz using an IIR filter\n"\
"  --spectrum-type TYPE         use PSD estimator TYPE [mean|median]\n"\
"\n"\
"  --segment-length N           set data segment length to N points\n"\
"  --number-of-segments N       set number of data segments to N\n"\
"\n"\
"  --low-frequency-cutoff F     do not filter below F Hz\n"\
"  --high-frequency-cutoff F    upper frequency cutoff in Hz\n"\
"\n"\
"  --minimum-mass MASS          set minimum component mass of bank to MASS\n"\
"  --maximum-mass MASS          set maximum component mass of bank to MASS\n"\
"  --minimal-match M            generate bank with minimal match M\n"\
"\n"\
"  --order ORDER                set post-Newtonian order of the waveform to ORDER\n"\
"                                 (newtonian|oneHalfPN|onePN|onePointFivePN|\n"\
"                                 twoPN|twoPointFive|threePN|threePointFivePN)\n"\
"  --approximant APPROX         set approximant of the waveform to APPROX\n"\
"                                 (TaylorT1|TaylorT2|TaylorT3|TaylorF1|TaylorF2|\n"\
"                                 PadeT1|PadeT2|EOB|BCV|SpinTaylorT3)\n"\
"  --space SPACE                grid up template bank with mass parameters SPACE\n"\
"                                 (Tau0Tau2|Tau0Tau3)\n"\
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
    /* parameters used to generate calibrated power spectrum */
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"channel-name",            required_argument, 0,                'c'},
    {"segment-length",          required_argument, 0,                'd'},
    {"number-of-segments",      required_argument, 0,                'e'},
    {"sample-rate",             required_argument, 0,                'g'},
    {"help",                    no_argument,       0,                'h'},
    {"low-frequency-cutoff",    required_argument, 0,                'i'},
    {"spectrum-type",           required_argument, 0,                'j'},
    {"calibration-cache",       required_argument, 0,                'p'},
    {"comment",                 required_argument, 0,                's'},
    {"enable-high-pass",        required_argument, 0,                't'},
    {"frame-cache",             required_argument, 0,                'u'},
    {"pad-data",                required_argument, 0,                'x'},
    {"debug-level",             required_argument, 0,                'z'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {"resample-filter",         required_argument, 0,                'R'},
    /* template bank generation parameters */
    {"minimum-mass",            required_argument, 0,                'A'},
    {"maximum-mass",            required_argument, 0,                'B'},
    {"minimal-match",           required_argument, 0,                'C'},
    {"high-frequency-cutoff",   required_argument, 0,                'D'},
    {"order",                   required_argument, 0,                'E'},
    {"approximant",             required_argument, 0,                'F'},
    {"space",                   required_argument, 0,                'G'},
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
        "a:b:c:d:e:g:h:i:j:p:s:t:u:x:z:A:B:C:D:E:F:G:R:Z:", 
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

      case 'h':
        fprintf( stdout, USAGE );
        exit( 0 );
        break;

      case 'i':
        fLow = (REAL4) atof( optarg );
        if ( fLow < 40 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "low frequency cutoff is less than 40 Hz: "
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

      case 'p':
        /* create storage for the calibration frame cache name */
        optarg_len = strlen( optarg ) + 1;
        calCacheName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( calCacheName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
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

      case 'u':
        optarg_len = strlen( optarg ) + 1;
        frInCacheName = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( frInCacheName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'x':
        padData = (UINT4) atoi( optarg );
        if ( padData < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "number of seconds to pad from input data"
              "must be greater than 0: (%d specified)\n", 
              long_options[option_index].name, numSegments );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", numSegments );
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
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown order specified: "
              "%s (must be one of: TaylorT1, TaylorT2, TaylorT3, TaylorF1,\n"
              "TaylorF2, PadeT1, PadeF1, EOB, BCV or SpinTaylorT3)\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        haveApprox = 1;
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'G':
        if ( ! strcmp( "Tau0Tau2", optarg ) )
        {
          approximant = Tau0Tau2;
        }
        else if ( ! strcmp( "Tau0Tau3", optarg ) )
        {
          approximant = Tau0Tau3;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown order specified: "
              "%s (must be one of: Tau0Tau2 or Tau0Tau3)\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        haveSpace = 1;
        ADD_PROCESS_PARAM( "string", "%s", optarg );
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

  /* check that the frame caches have been specified */
  if ( ! frInCacheName )
  {
    fprintf( stderr, "--frame-cache must be specified\n" );
    exit( 1 );
  }
  if ( ! calCacheName )
  {
    fprintf( stderr, "--calibration-cache must be specified\n" );
    exit( 1 );
  }

  /* check that the bank parameters have been specified */
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
  
  return 0;
}

#undef ADD_PROCESS_PARAM

