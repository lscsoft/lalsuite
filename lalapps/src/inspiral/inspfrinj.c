/*
*  Copyright (C) 2007 Duncan Brown, Stephen Fairhurst
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/*----------------------------------------------------------------------- 
 * 
 * File Name: inspfrinj.c
 *
 * Author: Fairhurst, S. (based on inspiral.c by Brown, D.A.)
 * 
 * 
 *-----------------------------------------------------------------------
 */

/**
 * \file
 * \ingroup lalapps_inspiral
 *
 * <dl>
 * <dt>Name</dt><dd>
 * \c lalapps_inspfrinj --- performs inspiral injections into frame data.</dd>
 *
 * <dt>Synopsis</dt><dd>
 * \code
 * lalapps_inspfrinj [options]
 *
 *  [--help ]
 *  [--verbose ]
 *  [--version ]
 *  [--user-tag STRING ]
 *  [--comment STRING ]
 *
 *   --gps-start-time SEC
 *  [--gps-start-time-ns NS ]
 *   --gps-end-time SEC
 *  [--gps-end-time-ns NS ]
 *
 *  [--frame-cache FRAME_CACHE
 *   --channel-name CHAN ]
 *
 *  [--ifo  IFO
 *   --sample-rate SAMPLE_RATE ]
 *
 *   --calibration-cache CAL_CACHE
 *   --calibrated-data TYPE
 *
 *  [--num-resp-points N ]
 *
 *   --injection-file FILE
 *   --injection-channel INJ
 *  [--inject-overhead ]
 *  [--inject-safety SEC ]
 *
 *  [--write-raw-data ]
 *  [--write-inj-only ]
 *  [--write-raw-plus-inj ]
 *
 *   --output-frame-length OUTPUT_LENGTH
 *  [--output-file-name OUTPUT_NAME ]
 *
 * \endcode</dd>
 *
 * <dt>Description</dt><dd>
 *
 * \c lalapps_frinspinj generates injections and stores them as
 * frame files.  This code is essentially a copy of those parts of
 * \c lalapps_inspiral responsible for reading in the data and
 * performing the injections.  The injection details are either read in
 * from frame files containin an injection channel or from the
 * \c sim_inspiral table of a LIGO lightweight xml file or Is a
 * \c sim_inspiral table is given, the injections are generated
 * using LALFindChirpInjectSignals.  The injection data is output to frame
 * files.  The length of the output frames is specified by
 * <tt>output-frame-length</tt>.
 *
 * The code is capable of generating both calibrated and uncalibrated
 * injections.  To obtain calibrated h(t) injections the option
 * <tt>calibrated-data</tt> must be specified (at present the code only
 * works for real_4 data).  Alternatively, if a \c CAL_CACHE is
 * given, then uncalibrated injections will be produced using the a
 * response function generated from the calibration data in the cache.  The
 * calibration coefficients will by default be averaged over the duration
 * of the data.  When working with uncalibrated data, the number of points
 * used to determine the response function will have an effect, albeit
 * minor, on the injection signals.  By default, the value
 * <tt>num-resp-points</tt> is set to 4194304.  This matches the value used
 * by \c lalapps_inspiral when working with 256 second segments and
 * a channel sampled at 16384 Hz.
 *
 * If input data is given by a \c FRAME_CACHE, the code will read in
 * this data.  In this case, the sample rate and ifo will be determined
 * from the data and channel name respectively.  The output frames can
 * contain any of the following: raw data, injection only and raw plus
 * injection.
 *
 * If no input data is provided, then the \c IFO and
 * \c SAMPLE_RATE must be given on the command line.  Since no input
 * data is read in, the output data can only contain the injection only
 * channel.
 *
 * The output of \c lalapps_frinspinj is an xml file and one or several
 * frame files.  The xml file contains
 * \c process, \c process_params, \c search_summary and
 * \c sim_inspiral tables and is named (unless set on the command line):
 *
 * <tt>IFO-INSPFRINJ_USERTAG-GPSSTARTTIME-DURATION.xml</tt>\\
 *
 * where \c GPSSTARTTIME and \c DURATION are the start time and
 * duration passed to the code.  If <tt>--output-file-name</tt> is set to
 * \c OUTPUT on the command line, the xml will be named
 *
 * <tt>OUTPUT_USERTAG-GPSSTARTTIME-DURATION.xml</tt>.\\
 *
 * The \c sim_inspiral table contains a list of all injections which
 * were performed into the specified data.  The length of output frame
 * files is specified by <tt>output-frame-length</tt>.  Unless otherwise set
 * on the command line, output frames are named:
 *
 * <tt>IFO-INSPFRINJ_USERTAG-FRAMESTARTTIME-FRAMELENGTH.gwf</tt>\\
 *
 * The channels stored in the output frame are determined by which of
 * <tt>--write-raw-data</tt>, <tt>--write-inj-only</tt> and
 * <tt>--write-raw-plus-inj</tt> are selected.  The channel names are
 * \c CHAN, \c CHAN_INSP_INJ_ONLY and
 * \c CHAN_RAW_PLUS_INSP_INJ.  In the case where no input data is
 * given, the single output channel name is \c IFO:STRAIN_INSP_INJ_ONLY.</dd>
 *
 * <dt>Options</dt><dd>
 * <ul>
 *
 * <li><tt>--help</tt>:  Optional.  Print a help message and exit.</li>
 *
 * <li><tt>--verbose</tt>: Optional.  Print out the author, CVS version
 * and tag information and exit.</li>
 *
 * <li><tt>--version</tt>: Optional.  Print out the author, CVS version and
 * tag information and exit.
 *   </li>
 * <li><tt>--user-tag</tt> \c USERTAG: Optional. Set the user tag
 * for this job to be \c USERTAG. May also be specified on the command
 * line as <tt>-userTag</tt> for LIGO database compatibility.  This will
 * affect the naming of the output file.</li>
 *
 * <li><tt>--comment</tt> \c string: Optional. Add \c string
 * to the comment field in the process table. If not specified, no comment
 * is added. </li>
 *
 * <li><tt>--gps-start-time</tt> <tt>GPS seconds</tt>: Required.  Start
 * time, GPS seconds.</li>
 *
 * <li><tt>--gps-start-time-ns</tt> <tt>GPS seconds</tt>: Optional.
 * Start time, GPS nanoseconds.  If not specified, then set to zero.</li>
 *
 * <li><tt>--gps-end-time</tt> <tt>GPS seconds</tt>: Required.  End
 * time, GPS seconds.</li>
 *
 * <li><tt>--gps-end-time-ns</tt> <tt>GPS seconds</tt>: Optional.
 * Start time, GPS nanoseconds.  If not specified, then set to zero.</li>
 *
 * <li><tt>--frame-cache</tt> \c FRAME_CACHE: Optional.  Name of
 * the frame cache file.  If this is not specified then the code will not
 * read in any input data.  Hence it will be unable to produce the raw or
 * raw_plus_inj frames.</li>
 *
 * <li><tt>--channel-name</tt> \c CHAN: Optional.  Input data
 * channel name.  This must be specified if the <tt>frame-cache</tt> file
 * is specified.
 *  </li>
 * <li><tt>--ifo</tt> \c IFO: Optional.  Set input interferometer
 * name to \c IFO.  This information is required unless a
 * <tt>channel-name</tt> is specified (in this case the \c IFO is
 * obtained from the channel name).</li>
 *
 * <li><tt>--sample-rate</tt> \c SAMPLE_RATE: Optional.  Set the
 * sample rate (in Hz) at which the injections should be produced.  This is
 * required unless a frame cache is given (in this case, the sample rate of
 * the injections matches the input data).</li>
 *
 * <li><tt>--calibration-cache</tt> \c CAL_CACHE: Optional.
 * Specify the \c CAL_CACHE to obtain the calibration information.</li>
 *
 * <li><tt>--calibrated-data </tt>\c TYPE:  Optional.  Use
 * calibrated data as input/output.  This will produce strain data
 * injections.  This option must be specified if no calibration cache is
 * given.
 *   </li>
 * <li><tt>--num-resp-points</tt> \c N: Optional.  The number of
 * points used to generate the response function.  If generating
 * uncalibrated injections, this can actually change slightly the injection
 * data.  If not specified, then set to default value of 4194304.
 *  </li>
 * <li><tt>--injection-channel</tt> \c INJ: Optional.  If this
 * option is specified, the code will read in injection data from the INJ
 * channel and add it to the raw data.  Either this option or
 * <tt>--injection-file</tt> must be specified.
 *  </li>
 * <li><tt>--injection-file</tt> \c FILE: Optional. Read in the
 * injection data from \c FILE.  This should contain a
 * \c sim_inspiral table with a list of injections.  Either this
 * option or <tt>--injection-channel</tt> must be specified.</li>
 *
 * <li><tt>--inject-overhead</tt>: Optional.  Perform the injections
 * directly overhead the interferometer and optimally oriented.</li>
 *
 * <li><tt>--inject-safety</tt> \c SEC: Optional.  Inject signals
 * whose end time is up to \c SEC after the gps end time.
 *  </li>
 * <li><tt>--write-raw-data</tt>:  Optional.  Write out a frame
 * containing the raw data.  This can only be done if the raw data has been
 * read in from frames.</li>
 *
 * <li><tt>--write-inj-only</tt>:  Optional.  Write out a frame
 * containing only the injection data.</li>
 *
 * <li><tt>--write-raw-plus-inj</tt>:  Optional.  Write out a frame
 * containing the raw data with the injections added.
 *   </li>
 * <li><tt>--output-frame-length</tt> \c OUTPUT_LENGTH:  Required.
 * Specify the length of the output frames to be \c output_length
 * seconds.</li>
 *
 * <li><tt>--output-file-name</tt> \c OUTPUT_NAME : Optional.
 * Set the output file name.  This overrides the default file naming
 * convention given above.
 * </li>
 * </ul></dd>
 *
 * <dt>Example</dt><dd>
 * \code
 * lalapps_frinspinj --gps-start-time 732758030  --gps-end-time 732760078 \
 *   --frame-cache cache/L-732758022-732763688.cache --channel-name L1:LSC-AS_Q\
 *   --calibration-cache cache_files/L1-CAL-V03-729273600-734367600.cache \
 *   --injection-file HL-INJECTIONS_4096-729273613-5094000.xml \
 *   --inject-safety 50 --output-frame-length 16 \
 *   --write-raw-data   --write-inj-only   --write-raw-plus-inj
 * \endcode</dd>
 *
 * <dt>Author</dt><dd>
 * Steve Fairhurst</dd>
 * </dl>
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <time.h>
#include <math.h>

#include <series.h>
#include <lalappsfrutils.h>

#include <lal/LALConfig.h>
#include <lal/LALgetopt.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/LALFrStream.h>
#include <lal/Calibration.h>
#include <lal/LALFrameIO.h>
#include <lal/FindChirp.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOMetadataInspiralUtils.h>
#include <lal/LIGOLwXMLlegacy.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>

#include <LALAppsVCSInfo.h>

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "inspiral"

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
INT4  numRespPoints     = -1;           /* num points for calc response */
CHAR  *fqChanName       = NULL;         /* name of data channel         */
CHAR  *frInCacheName    = NULL;         /* cache file containing frames */
CHAR  *injCacheName     = NULL;         /* inj cache file for inj frames*/
CHAR   ifo[3];                          /* two character ifo code       */
CHAR   outfileName[FILENAME_MAX];       /* output file name             */

enum { undefined, real_4, real_8 } calData = undefined; /* cal data type*/
/* data conditioning parameters */
INT4   sampleRate       = -1;           /* sample rate of filter data   */
INT4   frameLength      = -1;           /* length of output frames      */
INT4   injectSafety     = 0;            /* safety length in injections  */
UINT4  numFiles         = 0;            /* number of output files needed*/

CHAR  *calCacheName     = NULL;         /* location of calibration data */
CHAR  *calFileName      = NULL;         /* location of calibration file */
CHAR  *injectionFile    = NULL;         /* name of file containing injs */
CHAR  *injChanName      = NULL;         /* the injection channel name   */

REAL4 injFlow           = 0;            /* injection start frequency    */
int   injectOverhead    = 0;            /* inject h+ into detector      */
int   numInjections     = 0;
SimInspiralTable *injections = NULL;
SimInspiralTable    *thisInj = NULL;

/* output parameters */
CHAR  *userTag          = NULL;         /* string the user can tag with */
int    writeRawData     = 0;            /* write the raw data to frame  */
int    writeInjOnly     = 0;            /* write the inj data to frame  */
int    writeRawPlusInj  = 0;            /* write raw plus inj to frame  */
int    writeReal8Frame  = 0;            /* write frames as real 8       */
int    outCompress = 0;
/* other command line args */
CHAR comment[LIGOMETA_COMMENT_MAX];     /* process param comment        */

int main( int argc, char *argv[] )
{
  /* lal function variables */
  LALStatus             XLAL_INIT_DECL(status);

  /* frame input data */
  LALCache     *frInCache = NULL;
  LALFrStream     *frStream = NULL;
  FrChanIn      frChan;
  FrChanIn      injChan;

  /* raw input data storage */
  REAL4TimeSeries               chan;
  REAL4TimeSeries               inj;

  /* structures for preconditioning */
  COMPLEX8FrequencySeries       injResp;

  /* output data */
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         searchsumm;
  MetadataTable         siminspiral;
  ProcessParamsTable   *this_proc_param;
  LIGOLwXMLStream       results;

  /* counters and other variables */
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
  UINT4 k,n,j;
  CHAR  fname[FILENAME_MAX];
  UINT4 numPoints = 0;
  REAL8 tsLength;
  CHAR tmpChName[LALNameLength];

  /*
   *
   * initialization
   *
   */


  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  XLALSetErrorHandler( XLALAbortErrorHandler );


  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  XLALGPSTimeNow(&(proctable.processTable->start_time));
  XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME, lalAppsVCSIdentInfo.vcsId,
      lalAppsVCSIdentInfo.vcsStatus, lalAppsVCSIdentInfo.vcsDate, 0);
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );

  /* create the search summary table */
  searchsumm.searchSummaryTable = (SearchSummaryTable *)
    calloc( 1, sizeof(SearchSummaryTable) );

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
    snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
    snprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX, 
        " " );
  } 
  else 
  {
    snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
    snprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
  }


  /* the number of nodes for a standalone job is always 1 */
  searchsumm.searchSummaryTable->nnodes = 1;

  /* initialize the raw and injection data */
  memset( &chan, 0, sizeof(REAL4TimeSeries) );
  memset( &inj, 0, sizeof(REAL4TimeSeries) );


  /*
   *
   * read in the input data channels
   *
   */


  if ( frInCacheName )
  {
    /* set the params of the input data time series */
    chan.epoch = gpsStartTime;

    /* open a frame cache */
    frInCache = XLALCacheImport(frInCacheName);
    LAL_CALL( LALFrCacheOpen( &status, &frStream, frInCache ), &status );

    /* set the mode of the frame stream to fail on gaps or time errors */
    frStream->mode = LAL_FR_STREAM_VERBOSE_MODE;


    /*
     *
     *  Read in the raw data
     *
     */

    /* seek to required epoch and set chan name */
    LAL_CALL( LALFrSeek( &status, &(chan.epoch), frStream ), &status );
    frChan.name = fqChanName;

    /* determine the sample rate of the raw data */
    LAL_CALL( LALFrGetREAL4TimeSeries( &status, &chan, &frChan, frStream ),
        &status );

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
    /* FIXME:  loss of precision;  consider
    searchsumm.searchSummaryTable->in_start_time = searchsumm.searchSummaryTable->in_end_time = chan.epoch;
    XLALGPSAdd(&searchsumm.searchSummaryTable->in_end_time, chan.deltaT * (REAL8) chan.data->length);
    */
    searchsumm.searchSummaryTable->in_start_time = chan.epoch;
    tsLength = XLALGPSGetREAL8( &(chan.epoch) );
    tsLength += chan.deltaT * (REAL8) chan.data->length;
    XLALGPSSetREAL8( &(searchsumm.searchSummaryTable->in_end_time), tsLength );

    if ( vrbflg ) fprintf( stdout, "read channel %s from frame stream\n"
        "got %d points with deltaT %e\nstarting at GPS time %d sec %d ns\n", 
        chan.name, chan.data->length, chan.deltaT, 
        chan.epoch.gpsSeconds, chan.epoch.gpsNanoSeconds );

    /* 
     *
     * Read in the injection data from frames
     *
     */

    if ( injChanName )
    {
      /* set the params of the input data time series */
      inj.epoch = gpsStartTime;

      if ( injCacheName )
      {
        /* close current frame cache */
        LAL_CALL( LALFrClose( &status, &frStream ), &status );
        if ( frInCacheName ) XLALDestroyCache( frInCache );

        /* open injection frame cache */
        frInCache = XLALCacheImport(injCacheName);
        LAL_CALL( LALFrCacheOpen( &status, &frStream, frInCache ), &status );
      }

      /* seek to required epoch and set inj name */
      LAL_CALL( LALFrSeek( &status, &(inj.epoch), frStream ), &status );
      injChan.name = injChanName;

      /* determine the sample rate of the inj data */
      LAL_CALL( LALFrGetREAL4TimeSeries( &status, &inj, &injChan, frStream ),
          &status );

      /* determine the number of points to get, create storage for the data */
      numPoints = (UINT4) floor( ((REAL8) inputLengthNS) / 
          (inj.deltaT * 1.0e9) + 0.5 );
      LAL_CALL( LALSCreateVector( &status, &(inj.data), numPoints ), 
          &status );

      if ( vrbflg ) fprintf( stdout, "input inj Channel %s has sample interval"
          " (deltaT) = %e\nreading %d points from frame stream\n", fqChanName, 
          inj.deltaT, numPoints );

      /* read the data inj Channel time series from frames */
      LAL_CALL( LALFrGetREAL4TimeSeries( &status, &inj, &injChan, frStream ),
          &status );
      memcpy( &(inj.sampleUnits), &lalADCCountUnit, sizeof(LALUnit) );


      if ( vrbflg ) fprintf( stdout, "read inj Channel %s from frame stream\n"
          "got %d points with deltaT %e\nstarting at GPS time %d sec %d ns\n", 
          inj.name, inj.data->length, inj.deltaT, 
          inj.epoch.gpsSeconds, inj.epoch.gpsNanoSeconds );
    }

    /* close the frame file stream and destroy the frame cache */
    LAL_CALL( LALFrClose( &status, &frStream ), &status );
    if ( frInCacheName ) XLALDestroyCache( frInCache );
  }

  /* 
   *
   * Injections from file
   *
   */

  if ( injectionFile )
  {

    /* Create zeros on top of which to do the injections */

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
      snprintf( inj.name, LALNameLength, "%s:STRAIN", ifo );
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
     * inject signals into the zero data
     *
     */

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
      /* store number of injections in search summary */
      searchsumm.searchSummaryTable->nevents = numInjections;

      /* set the injection start frequency */
      if ( injFlow )
      {
        for( thisInj = injections; thisInj; thisInj = thisInj->next )
        {
          thisInj->f_lower = injFlow;
        }
      }

      /* create the response function */
      memset( &injResp, 0, sizeof(COMPLEX8FrequencySeries) );
      LAL_CALL( LALCCreateVector( &status, &(injResp.data), 
            numRespPoints / 2 + 1 ), &status );
      injResp.epoch = inj.epoch ;
      injResp.deltaF = 1.0 / ( numRespPoints * inj.deltaT );
      strcpy( injResp.name, inj.name );

      if ( calFileName )
      {
        XLAL_ERROR(XLAL_EERR, "Calibration frames no longer supported");
      }
      else if ( calCacheName )
      {
        XLAL_ERROR(XLAL_EERR, "Calibration cache no longer supported");
      }
      else 
      {
        /* generate a unity response function for h(t) */
        if ( vrbflg ) fprintf( stdout, "setting response to unity... " );
        injResp.sampleUnits = strainPerCount;
        for ( k = 0; k < injResp.data->length; ++k )
        {
          injResp.data->data[k] = 1.0;
        }
        if ( vrbflg ) fprintf( stdout, "done.\n" );

      }


      /* inject the signals, preserving the channel name (Tev mangles it) */
      snprintf( tmpChName, LALNameLength, "%s", inj.name );

      /* if injectOverhead option, then set inj.name to "ZENITH".  
       * This causes no detector site to be found in the injection code so
       * that the injection is done directly overhead (i.e. with a response 
       * function of F+ = 1; Fx = 0) */
      if ( injectOverhead )
      {
        snprintf( inj.name, LALNameLength, "ZENITH" );
      }

      LAL_CALL( LALFindChirpInjectSignals( &status, &inj, injections, 
            &injResp ), &status );
      snprintf( inj.name,  LALNameLength, "%s", tmpChName );

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
    REAL8TimeSeries real8Output;
    UINT4 length;

    memset( &output, 0, sizeof(REAL4TimeSeries) );
    output.deltaT = inj.deltaT;
    output.sampleUnits = inj.sampleUnits;
    output.data = (REAL4Vector *) LALCalloc( 1, sizeof(REAL4Vector) );
    length = numPoints / numFiles;
    output.data->length = length;

    memset( &real8Output, 0, sizeof(REAL8TimeSeries) );
    real8Output.deltaT         = output.deltaT;
    real8Output.sampleUnits    = output.sampleUnits;
    LAL_CALL( LALDCreateVector( &status, &(real8Output.data), 
         output.data->length ), &status );
    real8Output.data->length = output.data->length;
    
    for ( n = 0; n < numFiles; ++n )
    {
      outFrame = NULL;
      output.epoch.gpsSeconds  = gpsStartTime.gpsSeconds + n * frameLength;
      real8Output.epoch        = output.epoch;
      
      /* write the injection channel to frame */
      if ( writeInjOnly  ) 
      {
        strcpy( output.name, inj.name );
        output.data->data = inj.data->data + n * length;
        if ( injectionFile )
        {
          /* write out injection channel to INSP_INJ_ONLY */
          outFrame = fr_add_proc_REAL4TimeSeries( outFrame, &output, "ct", 
              "INSP_INJ_ONLY" );
        }
        else if ( injChanName )
        {
          /* write out injections, preserving input frame name */
          outFrame = fr_add_proc_REAL4TimeSeries( outFrame, &output, "ct", 
              NULL );
        }
      }

      /* write the raw/raw plus inj data to frame */
      if ( writeRawData || writeRawPlusInj ) 
      {

        strcpy( output.name, chan.name );
        output.data->data = chan.data->data + n * length;

        if ( !writeReal8Frame )
        {
          if ( writeRawData )
          {
            outFrame = fr_add_proc_REAL4TimeSeries( outFrame, &output, "ct", 
                NULL );
          }
          /* add injections into this file's data */
	  XLALAddREAL4TimeSeries( &output, &inj );

          if ( writeRawPlusInj )
          {
            strcpy( output.name, chan.name );
            outFrame = fr_add_proc_REAL4TimeSeries( outFrame, &output, "ct", 
                "PLUS_INSP_INJ" );
          }
        }
        else
        {
          strcpy( real8Output.name, output.name );
          
          if ( vrbflg ) fprintf( stdout, 
              "Casting data to real 8 before writing frame\n" );

          for ( j = 0 ; j < output.data->length ; ++j )
          {
            real8Output.data->data[j] = (REAL8) ( output.data->data[j] );
          }

          if ( writeRawData )
          {
            outFrame = fr_add_proc_REAL8TimeSeries( outFrame, &real8Output, 
                "ct", NULL );
          }

          for ( j = 0 ; j < output.data->length ; ++j )
          {
            real8Output.data->data[j] += (REAL8) 
              ( inj.data->data[j + n * length] );
          }

          if ( writeRawPlusInj )
          {
             outFrame = fr_add_proc_REAL8TimeSeries( outFrame, &real8Output, 
                 "ct", "PLUS_INSP_INJ" );
          }

        }

      }

      /* set the output file name */
      if( !outfileName[0] )
      {
        /* output name not specified, set to IFO-INSPFRINJ-EPOCH-LENGTH.gwf */
        snprintf( outfileName, FILENAME_MAX, "%s-INSPFRINJ", ifo );
      }

      if( userTag )
      {
        if(snprintf( fname, FILENAME_MAX, "%s_%s-%d-%d.gwf", outfileName,
            userTag, output.epoch.gpsSeconds, frameLength ) >= FILENAME_MAX)
          abort();
      }
      else
      {
        if(snprintf( fname, FILENAME_MAX, "%s-%d-%d.gwf", outfileName,
            output.epoch.gpsSeconds, frameLength ) >= FILENAME_MAX)
          abort();
      }

      if ( vrbflg ) fprintf( stdout, "writing frame data to %s... ", fname );
      frOutFile = FrFileONew( fname, 3 );
      FrameWrite( outFrame, frOutFile );
      FrFileOEnd( frOutFile );
      if ( vrbflg ) fprintf( stdout, "done\n" );
    }

    LALFree ( output.data );

    if( real8Output.data )
    { 
      LAL_CALL( LALDDestroyVector( &status, &(real8Output.data) ), &status );
    }
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
  if( userTag && outCompress )
  {
    if(snprintf( fname, FILENAME_MAX, "%s_%s-%d-%d.xml.gz", outfileName,
        userTag, gpsStartTime.gpsSeconds, 
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds ) >= FILENAME_MAX)
      abort();
  }
  else if( userTag && !outCompress )
  {
    if(snprintf( fname, FILENAME_MAX, "%s_%s-%d-%d.xml", outfileName,
        userTag, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds ) >= FILENAME_MAX)
      abort();
  }
  else if( !userTag && outCompress )
  {
    if(snprintf( fname, FILENAME_MAX, "%s-%d-%d.xml.gz", outfileName,
        gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds ) >= FILENAME_MAX)
      abort();
  }
  else
  {
    if(snprintf( fname, FILENAME_MAX, "%s-%d-%d.xml", outfileName,
        gpsStartTime.gpsSeconds, 
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds ) >= FILENAME_MAX)
      abort();
  }

  if ( vrbflg ) fprintf( stdout, "writing XML data to %s...\n", fname );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &results, fname ), &status );

  /* write the process table */
  if ( vrbflg ) fprintf( stdout, "  process table...\n" );
  snprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s", ifo );
  XLALGPSTimeNow(&(proctable.processTable->end_time));
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
      XLALFreeSimInspiral( &thisInj );
    }
  }

  /* close the output xml file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &results ), &status );
  if ( vrbflg ) fprintf( stdout, "done. XML file closed\n" );

  /* free the rest of the memory, check for memory leaks and exit */
  if ( injectionFile ) free ( injectionFile ); 
  if ( calCacheName ) free( calCacheName );
  if ( calFileName ) free( calFileName );
  if ( frInCacheName ) free( frInCacheName );
  if ( fqChanName ) free( fqChanName );

  if ( vrbflg ) fprintf( stdout, "checking memory leaks and exiting\n" );
  LALCheckMemoryLeaks();
  exit( 0 );
}

/* ------------------------------------------------------------------------- */

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
calloc( 1, sizeof(ProcessParamsTable) ); \
snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
    PROGRAM_NAME ); \
snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
    long_options[option_index].name ); \
snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

/*
 * 
 * USAGE
 *
 */
static void print_usage(char *program)
{
  fprintf(stderr,
      "Usage:  %s [options] [LIGOLW XML input files]\n" \
      "The following options are recognized.  Options not surrounded in [] are\n" \
      "required.\n" \
      " [--help]                           display this message\n"\
      " [--verbose]                        print progress information\n"\
      " [--version]                        print version information and exit\n"\
      " [--user-tag]             usertag   set the process_params usertag to usertag\n"\
      " [--comment]              string    set the process table comment to string\n"\
      "\n"\
      "  --gps-start-time        start_time GPS second of data start time\n"\
      " [--gps-start-time-ns]    start_ns   GPS nanosecond of data start time\n"\
      "  --gps-end-time          end_time   GPS second of data end time\n"\
      " [--gps-end-time-ns]      end_ns     GPS nanosecond of data end time\n"\
      "\n"\
      " [--frame-cache]          cache      frame cache with locations of data\n"\
      " [--calibration-file]     cal_file   frame file containing calibration data\n"\
      " [--calibration-cache]    cal_cache  file with location of calibration data\n"\
      " [--calibrated-data]      type       calibrated data of type (real_4 | real_8)\n"\
      " [--num-resp-points]      N          num points to determine response function (4194304)\n"\
      " [--channel-name]         chan       channel from which to read data\n"\
      "\n"\
      " [--injection-channel]    inj_chan   channel from which to read inj data\n"\
      " [--injection-cache]      inj_cache  cache with location of injection data\n"\
      " [--injection-file]       inj_file   xml file with injection details\n"\
      " [--inject-overhead]                 inject signals from overhead detector\n"\
      " [--inject-safety]        safety     inject signals ending up to safety\n" 
      "                                       seconds after end_time\n"\
      " [--injection-start-freq] flow       inject signals starting at flow (40Hz)\n"\
      "\n"\
      " [--write-raw-data]                  write out raw-data channel\n"\
      " [--write-inj-only]                  write out inj-only channel\n"\
      " [--write-raw-plus-inj]              write out raw plus inj channel\n"\
      " [--write-real8-frame]               write out real 8 frames\n"\
      "\n"\
      " [--output-frame-length]  len        length of output frames\n"\
      " [--output-file-name]     out        set file names to out-gpstime-length.gwf\n"\
      "                      if not set, default to ifo-inspfrinj-gpstime-length.gwf\n"\
      " [--write-compress]                  write compressed xml files\n"\
      "\n"\
      " [--ifo]                  ifo        specify the ifo (if not reading frames)\n"\
      " [--sample-rate]          rate       data sample rate (if not reading frames)\n"\
      "\n", program );
}

int arg_parse_check( int argc, char *argv[], MetadataTable procparams )
{
  /* LALgetopt arguments */
  struct LALoption long_options[] =
  {
    /* these options set a flag */
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"write-compress",          no_argument,       &outCompress,      1 },
    {"write-raw-data",          no_argument,       &writeRawData,     1 },
    {"write-inj-only",          no_argument,       &writeInjOnly,     1 },
    {"write-raw-plus-inj",      no_argument,       &writeRawPlusInj,  1 },
    {"inject-overhead",         no_argument,       &injectOverhead,   1 },
    {"write-real8-frame",       no_argument,       &writeReal8Frame,  1 },
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
    {"calibration-file",        required_argument, 0,                'q'},
    {"num-resp-points",         required_argument, 0,                'N'},
    {"sample-rate",             required_argument, 0,                'r'},
    {"ifo",                     required_argument, 0,                'i'},
    {"comment",                 required_argument, 0,                's'},
    {"frame-cache",             required_argument, 0,                'u'},
    {"injection-cache",         required_argument, 0,                'C'},
    {"injection-file",          required_argument, 0,                'w'},
    {"inject-safety",           required_argument, 0,                'S'},
    {"injection-channel",       required_argument, 0,                'I'},
    {"injection-start-freq",    required_argument, 0,                'L'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {"version",                 no_argument,       0,                'V'},
    {0, 0, 0, 0}
  };
  int c;
  ProcessParamsTable *this_proc_param = procparams.processParamsTable;


  /*
   *
   * parse command line arguments
   *
   */


  while ( 1 )
  {
    /* LALgetopt_long stores long option here */
    int option_index = 0;
    size_t LALoptarg_len;

    c = LALgetopt_long_only( argc, argv,
        "A:B:C:I:L:N:S:V:Z:"
        "a:b:c:d:f:hi:l:p:q:r:s:u:w:y:",
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
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        break;

      case 'a':
        /* set gps start seconds */
        {
          long int gstartt = atol( LALoptarg );
          if ( gstartt < 441417609 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS start time is prior to " 
                "Jan 01, 1994  00:00:00 UTC:\n"
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
          long int gstarttns = atol( LALoptarg );
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
          long int gendt = atol( LALoptarg );
          if ( gendt < 441417609 )
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
          long int gendtns = atol( LALoptarg );
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
          LALoptarg_len = strlen( LALoptarg ) + 1;
          fqChanName = (CHAR *) calloc( LALoptarg_len, sizeof(CHAR) );
          memcpy( fqChanName, LALoptarg, LALoptarg_len );
          ADD_PROCESS_PARAM( "string", "%s", LALoptarg );

          /* check that we have a proper channel name */
          if ( ! (channamptr = strstr( fqChanName, ":" ) ) )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "channel name must be a full LIGO channel name "
                "e.g. L1:LSC-AS_Q\n(%s specified)\n",
                long_options[option_index].name, LALoptarg );
            exit( 1 );
          }

          /* copy the first two characters to the ifo name */
          memset( ifo, 0, sizeof(ifo) );
          memcpy( ifo, LALoptarg, sizeof(ifo) - 1 );
        }
        break;

      case 'd':
        /* set length of output frames */
        frameLength = (INT4) atoi( LALoptarg );
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
        /* set output file name */
        if ( snprintf( outfileName, FILENAME_MAX, "%s", LALoptarg ) < 0 )
        {
          fprintf( stderr, "invalid argument to --%s\n"
              "outfile name %s too long: string truncated\n",
              long_options[option_index].name, LALoptarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        break;

      case 'h':
        /* help message */
        print_usage(argv[0]);
        exit( 0 );
        break;

      case 'p':
        /* create storage for the calibration frame cache name */
        LALoptarg_len = strlen( LALoptarg ) + 1;
        calCacheName = (CHAR *) calloc( LALoptarg_len, sizeof(CHAR));
        memcpy( calCacheName, LALoptarg, LALoptarg_len );
        ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        break;

      case 'q':
        /* create storage for the calibration frame file name */
        LALoptarg_len = strlen( LALoptarg ) + 1;
        calFileName = (CHAR *) calloc( LALoptarg_len, sizeof(CHAR));
        memcpy( calFileName, LALoptarg, LALoptarg_len );
        ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        break;


      case 'N':
        /* store the number of points used in computing the response */
        numRespPoints = (INT4) atoi( LALoptarg );
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
          if ( ! strcmp( "real_4", LALoptarg ) )
          {
            calData = real_4;
          }
          else if ( ! strcmp( "real_8", LALoptarg ) )
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
                long_options[option_index].name, LALoptarg);
            exit( 1 );
          }
          ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        }
        break;

      case 'r':
        /* set the sample rate */
        sampleRate = (INT4) atoi( LALoptarg );
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

      case 'L':
        /* set the injection start frequency */
        injFlow = (REAL4) atof ( LALoptarg );
        if ( injFlow <= 0 )
        {
          fprintf( stderr, "invalide argument to --%s:\n"
              "injections must start at a positive frequency."
              "(%f specified) \n",
              long_options[option_index].name, injFlow );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%f", injFlow );
        break;

      case 'i':
        {
          /* create storage for the ifo name and copy it */
          memset( ifo, 0, sizeof(ifo) );
          memcpy( ifo, LALoptarg, sizeof(ifo) - 1 );
          ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        }
        break;


      case 's':
        if ( strlen( LALoptarg ) > LIGOMETA_COMMENT_MAX - 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "comment must be less than %d characters\n",
              long_options[option_index].name, LIGOMETA_COMMENT_MAX );
          exit( 1 );
        }
        else
        {
          snprintf( comment, LIGOMETA_COMMENT_MAX, "%s", LALoptarg);
        }
        break;

      case 'S':
        injectSafety = (INT4) atoi( LALoptarg );
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

      case 'I':
        {
          /* create storage for the injection channel name and copy it */
          LALoptarg_len = strlen( LALoptarg ) + 1;
          injChanName = (CHAR *) calloc( LALoptarg_len, sizeof(CHAR) );
          memcpy( injChanName, LALoptarg, LALoptarg_len );
          ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        }
        break;

      case 'u':
        /* create storage for the input frame cache name */
        LALoptarg_len = strlen( LALoptarg ) + 1;
        frInCacheName = (CHAR *) calloc( LALoptarg_len, sizeof(CHAR) );
        memcpy( frInCacheName, LALoptarg, LALoptarg_len );
        ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        break;

      case 'C':
        /* create storage for the input frame cache name */
        LALoptarg_len = strlen( LALoptarg ) + 1;
        injCacheName = (CHAR *) calloc( LALoptarg_len, sizeof(CHAR) );
        memcpy( injCacheName, LALoptarg, LALoptarg_len );
        ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        break;

      case 'w':
        /* create storage for the injection file name */
        LALoptarg_len = strlen( LALoptarg ) + 1;
        injectionFile = (CHAR *) calloc( LALoptarg_len, sizeof(CHAR));
        memcpy( injectionFile, LALoptarg, LALoptarg_len );
        ADD_PROCESS_PARAM( "string", "%s", LALoptarg );
        break;

      case 'Z':
        /* create storage for the usertag */
        LALoptarg_len = strlen( LALoptarg ) + 1;
        userTag = (CHAR *) calloc( LALoptarg_len, sizeof(CHAR) );
        memcpy( userTag, LALoptarg, LALoptarg_len );

        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
            PROGRAM_NAME );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "-userTag" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            LALoptarg );
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "LIGO/LSC Inspiral Injection Program\n" 
            "Steve Fairhurst <sfairhur@gravity.phys.uwm.edu>\n");
        XLALOutputVCSInfo(stderr, lalAppsVCSInfoList, 0, "%% ");
        exit( 0 );
        break;

      case '?':
        print_usage(argv[0]);
        exit( 1 );
        break;

      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        print_usage(argv[0]);
        exit( 1 );
    }
  }

  if ( LALoptind < argc )
  {
    fprintf( stderr, "extraneous command line arguments:\n" );
    while ( LALoptind < argc )
    {
      fprintf ( stderr, "%s\n", argv[LALoptind++] );
    }
    exit( 1 );
  }

  /* check flags and store in process_params */

  if ( writeRawData )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--write-raw-data" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  if ( writeInjOnly )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--write-inj-only" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }
  if ( writeRawPlusInj )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--write-raw-plus-inj" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  if ( writeReal8Frame )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--write-real8-frame" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }
  
  /* check inject-overhead option */
  if ( injectOverhead )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--inject-overhead" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
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
  XLALINT8NSToGPS( &gpsStartTime, gpsStartTimeNS );

  /* end time specified */
  if ( ! gpsEndTimeNS )
  {
    fprintf( stderr, "--gps-end-time must be specified\n" );
    exit( 1 );
  }
  XLALINT8NSToGPS( &gpsEndTime, gpsEndTimeNS );

  /* end after start */
  if ( gpsEndTimeNS <= gpsStartTimeNS )
  {
    fprintf( stderr, "invalid gps time range: "
        "start time: %d, end time %d\n",
        gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds );
    exit( 1 );
  }

  /* check that we have injections, either in xml or in the frames already */
  if ( !injectionFile && !injChanName )
  {
    fprintf( stderr, 
        "either --injection-file or --injection-channel must be specified\n");
    exit( 1 );
  }
  if ( injectionFile && injChanName )
  {
    fprintf( stderr, 
        "Only one of --injection-file and --injection-channel may be given\n");
    exit( 1 );
  }

  /* if an injection channel has been specified, need a frame cache */
  if ( injChanName && !frInCacheName )
  {
    fprintf( stderr, 
        "If --injection-channel specified, also require --frame-cache\n");
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

  /* check that have one of: calibrated-data or a calibration-cache
   * or a calibration frame file */
  if ( ! (calCacheName || calFileName) && ! calData )
  {
    fprintf( stderr, "Either --calibration-cache must be specified,\n" 
        "or must run on --calibrated-data.\n");
    exit( 1 );
  }
  else if ( (calCacheName && calData) || (calFileName && calData) || (calCacheName && calFileName) )
  {
    fprintf( stderr, 
        "Only one of --calibration-cache and --calibration and --calibrated-data\n"
        "should be specified\n.");
    exit( 1 );
  }

  /* if a calibration frame file, need to know which readout point is to
   * be used (e.g., DARM_ERR or AS_Q) -- this is taken from a channel name */
  if ( calFileName && ! fqChanName )
  {
    fprintf( stderr, "--channel-name must be specified\n" );
    exit( 1 );
  }

  /* check that we have then number of points to determine response fn */
  if ( numRespPoints < 0 )
  {
    if ( vrbflg ) fprintf( stdout, "--num-resp-points not specified\n"    
        "This gives the number of points used to obtain response,\n"
        "Response has numRespPoints/2 + 1 points,\n"
        "Frequency resolution of 1/(numRespPoints * delta T).\n"
        "Setting it to a default value of 256 * 16384 = 4194304\n");
    numRespPoints = 4194304;
  }


  return 0;
}

#undef ADD_PROCESS_PARAM
