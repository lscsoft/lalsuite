/*
 *  Copyright (C) 2007 Gareth Jones
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
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

/*-----------------------------------------------------------------------
 *
 * File Name: sned.c
 *
 * Author: Keppel, D. G., Jones, G. W.
 *
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
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpTD.h>
#include <lal/FindChirpBCV.h>
#include <lal/FindChirpBCVSpin.h>
#include <lal/FindChirpChisq.h>
#include <lal/LALTrigScanCluster.h>
#include <lal/PrintFTSeries.h>
#include <lal/ReadFTSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/GenerateInspiral.h>
#include <lal/TimeSeries.h>
#include <lal/VectorOps.h>
#include <lal/LALFrameL.h>

#include <LALAppsVCSInfo.h>

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "sned"

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
calloc( 1, sizeof(ProcessParamsTable) ); \
snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
    PROGRAM_NAME ); \
snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
    long_options[option_index].name ); \
snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

#define MAX_PATH 4096

#define USAGE \
  "lalapps_sned [options]\n"\
"\nDefaults are shown in brackets\n\n" \
"  --help                   display this message\n"\
"  --verbose                be verbose\n"\
"  --version                version info\n"\
"  --debug-level LEVEL      set the LAL debug level to LEVEL\n"\
"  --spinning-search        use the normalization for a spinning search\n"\
"                           instead of for a non-spinning search\n"\
"  --inject-overhead        inject signals from overhead detector\n"\
"  --write-chan             write out time series showing inspiral waveform\n"\
"  --inj-file    FILE       xml FILE contains injections\n"\
"  --coire-flag             use this if inj file is a coire file\n"\
"  --output-file FILE       FILE for output\n"\
"  --f-lower     FREQ       freq at which to begin integration\n"\
"  --ligo-only              only normalize the eff_dist columns for\n"\
"                           LIGO detectors\n"\
"  --snr-threshold SNR      for simulating full search, will print triggers with\n\
SNRs greater than this to a Found xml file. Disables --output-file.\n"\
"\n"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

static void destroyCoherentGW( CoherentGW *waveform );

static void destroyCoherentGW( CoherentGW *waveform )
{
  if ( waveform->h )
  {
    XLALDestroyREAL4VectorSequence( waveform->h->data );
    LALFree( waveform->a );
  }
  if ( waveform->a )
  {
    XLALDestroyREAL4VectorSequence( waveform->a->data );
    LALFree( waveform->a );
  }
  if ( waveform->phi )
  {
    XLALDestroyREAL8Vector( waveform->phi->data );
    LALFree( waveform->phi );
  }
  if ( waveform->f )
  {
    XLALDestroyREAL4Vector( waveform->f->data );
    LALFree( waveform->f );
  }
  if ( waveform->shift )
  {
    XLALDestroyREAL4Vector( waveform->shift->data );
    LALFree( waveform->shift );
  }

  return;
}

extern int vrbflg;           /* verbocity of lal function */
int writechan;               /* whether to write chan txt files */
int injoverhead;             /* perform inj overhead if this option is set */
int coireflg;                /* is input file coire (1) or inj (null) */
int nonSpinningSearch = 1;   /* normalize for a non-spinning search */
int ligoOnly = 0;            /* only normalize the LIGO eff_dist columns */
int snrflg;                  /* use sned for toy MC (1) or normal (null) */

int main( int argc, char *argv[] )
{
  LALStatus                     status = blank_status;

  UINT4                         k;
  UINT4                         kLow;
  UINT4                         kHi;
  REAL4                         fSampling;
  REAL4                         fReSampling;
  REAL4                         fLow = -1.0;
  INT4                          numPoints;
  INT4                          numRawPoints;
  REAL8                         deltaT;
  REAL8                         deltaTReSample;
  REAL8                         deltaF;

  ResampleTSParams              resampleParams;

  REAL4                         UNUSED statValue;

  /* vars required to make freq series */
  LIGOTimeGPS                   epoch = { 0, 0 };
  LIGOTimeGPS                   gpsStartTime = {0, 0}; 
  REAL8                         f0 = 0.;
  REAL8                         offset = 0.;
  INT8                          waveformStartTime = 0;

  /* files contain PSD info */
  CHAR                         *injectionFile = NULL;         
  CHAR                         *outputFile    = NULL;         

  COMPLEX8Vector               *unity = NULL;
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};

  int                           numInjections = 0;
  int                           numTriggers = 0;
  UINT8                         eventNumber = 0;

  /* template bank simulation variables */
  INT4                         injSimCount = 0;
  SimInspiralTable            *injectionHead  = NULL;
  SimInspiralTable            *thisInjection  = NULL;
  SimInspiralTable            thisNonSpinningInjection;
  SnglInspiralTable           *snglHead       = NULL;
  SearchSummaryTable          *searchSummHead = NULL;
  /* SummValueTable              *summValueHead  = NULL; */

  /* raw input data storage */
  COMPLEX8FrequencySeries       *resp          = NULL;
  COMPLEX8FrequencySeries       *detTransDummy = NULL;
  REAL4TimeSeries               *chan          = NULL;
  REAL4TimeSeries               *chanDummy     = NULL;
  RealFFTPlan                   *pfwd          = NULL;
  COMPLEX8FrequencySeries       *fftData       = NULL;
  COMPLEX8FrequencySeries       *fftStandardData = NULL;
  REAL8                          thisSigmasq     = 0;
  REAL8                          spinningSigmasqVec[6] = {1, 1, 1, 1, 1, 1};
  REAL8                          thisStandardSigmasq = 0;
  REAL8                          standardSigmasqVec[6] = {1, 1, 1, 1, 1, 1};
  REAL8                          thisMixedSigmasq = 0;
  REAL8                          mixedSigmasqVec[6] = {1, 1, 1, 1, 1, 1};
  REAL8                          dynRange      = 1./(3.0e-23);

  /* needed for inj */
  CoherentGW                 waveform;
  CoherentGW                 nonSpinningWaveform;
  PPNParamStruc              ppnParams;
  DetectorResponse           detector;
  INT4                       ifos[6] = {1, 1, 1, 1, 1, 1};
  InterferometerNumber       ifoNumber = LAL_UNKNOWN_IFO;

  /* output data */
  LIGOLwXMLStream       xmlStream;
  MetadataTable         proctable;
  MetadataTable         outputTable;
  MetadataTable         procparams;
  CHAR                  fname[256];         
  CHAR                  comment[LIGOMETA_COMMENT_MAX];
  ProcessParamsTable   *this_proc_param = NULL;

  CHAR   chanfilename[FILENAME_MAX];

  /* vars needed for toy Monte Carlo */ 
  REAL8                 SNRsq_threshold = 0;
  REAL8                 Store_ifoSNRsq[6];
  SimInspiralTable      *thisMissedSim = NULL;
  SimInspiralTable      *thisFoundSim = NULL;
  SimInspiralTable      *headMissedSim = NULL;
  SimInspiralTable      *headFoundSim = NULL;
  SnglInspiralTable     *H1FoundSngl = NULL;
  SnglInspiralTable     *headFoundSngl = NULL;
  SnglInspiralTable     *L1FoundSngl = NULL;

  /* set initial debug level */
  set_debug_level("1");

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  XLALGPSTimeNow(&(proctable.processTable->start_time));
  XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME, LALAPPS_VCS_IDENT_ID,
      LALAPPS_VCS_IDENT_STATUS, LALAPPS_VCS_IDENT_DATE, 0);
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );

  /* look at input args, write process params where required */
  while ( 1 )
  {
    /* getopt arguments */
    static struct option long_options[] =
    {
      /* these options set a flag */
      /* these options do not set a flag */
      {"help",                    no_argument,       0,                 'h'},
      {"verbose",                 no_argument,       &vrbflg,            1 },
      {"version",                 no_argument,       0,                 'V'},
      {"inj-file",                required_argument, 0,                 'd'},
      {"comment",                 required_argument, 0,                 'e'},
      {"output-file",             required_argument, 0,                 'f'},
      {"coire-flag",              no_argument,       &coireflg,          1 },
      {"spinning-search",         no_argument,       &nonSpinningSearch, 0 },
      {"write-chan",              no_argument,       &writechan,         1 },
      {"inject-overhead",         no_argument,       &injoverhead,       1 },
      {"f-lower",                 required_argument, 0,                 'g'},
      {"ligo-only",               no_argument,       &ligoOnly,          1 },
      {"debug-level",             required_argument, 0,                 'z'}, 
      {"snr-threshold",           required_argument, 0,                 's'},
      {0, 0, 0, 0}
    };
    int c;

    /*
     *
     * parse command line arguments
     *
     */

    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, "a:b:c:d:e:f:g:z:s:hV", long_options,
        &option_index );

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

      case 'h':
        fprintf( stderr, USAGE );
        exit( 0 );
        break;

      case 'd':
        /* create storage for the injection file name */
        optarg_len = strlen( optarg ) + 1;
        injectionFile = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( injectionFile, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'f':
        /* create storage for the output file name */
        optarg_len = strlen( optarg ) + 1;
        outputFile = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( outputFile, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'g':
        fLow = (REAL4) atof( optarg );
        ADD_PROCESS_PARAM( "float", "%e", fLow );
        break;

      case 'e':
        if ( strlen( optarg ) > LIGOMETA_COMMENT_MAX - 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "comment must be less than %d characters\n",
              long_options[option_index].name, LIGOMETA_COMMENT_MAX );
          exit( 1 );
        }
        else
        {
          snprintf( comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
        }
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout,
            "spin-normalize the effective distance for spinning injections\n"
            "Drew Keppel and Gareth Jones\n");
        XLALOutputVersionString(stderr, 0);
        exit( 0 );
        break;

      case 'z':
        set_debug_level( optarg );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 's':
        SNRsq_threshold = (REAL8) pow( atof( optarg ), 2);
        snrflg = 1;
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

  if (fLow <= 0 )
  {
    fprintf( stderr, "Error: --f-lower is a required option.\n" );
    exit( 1 );
  }

  if ( snrflg )
  {
    fSampling = 4096.;
  }
  else
  {
    fSampling = 16384.;
  }

  fReSampling     = 4096.;
  numPoints       = 1048576; 
  numRawPoints    = floor( numPoints * fSampling / fReSampling + 0.5 );
  deltaT          = 1./fSampling;
  deltaTReSample  = 1./fReSampling;
  deltaF          = fReSampling / numPoints;

  /* check the input arguments */
  if ( injectionFile == NULL )
  {
    fprintf( stderr, "Must specify the --injection-file\n" );
    exit( 1 );
  }

  if ( !snrflg && outputFile == NULL )
  {
    fprintf( stderr, "Must specify the --output-file\n" );
    exit( 1 );
  }

  if ( vrbflg )
  {
    fprintf( stdout, "injection file is %s\n", injectionFile );
    if ( snrflg )
    {
      fprintf( stdout, "output files are Found_Injections.xml and Missed_Injections.xml\n" );
      fprintf(stdout, "SNRsq_threshold is %f\n", SNRsq_threshold );
    }
    else
    {
      fprintf( stdout, "output file is %s\n", outputFile );
    }
  }

  /* read in injections from injection file */
  /* set endtime to 0 so that we read in all events */
  if ( vrbflg ) fprintf( stdout, "Reading sim_inspiral table of %s\n",
      injectionFile );
  LAL_CALL(numInjections = SimInspiralTableFromLIGOLw( &injectionHead,
        injectionFile, 0, 0), &status);
  if ( vrbflg ) fprintf( stdout,
      "Read %d injections from sim_inspiral table of %s\n", numInjections,
      injectionFile );

  if (coireflg)
  {
    if ( vrbflg ) fprintf( stdout, "Reading sngl_inspiral table of %s\n",
        injectionFile );
    LAL_CALL(numTriggers = LALSnglInspiralTableFromLIGOLw(&snglHead,
          injectionFile, 0, -1), &status);
    if ( vrbflg ) fprintf( stdout,
        "Read %d triggers from sngl_inspiral table of %s\n", numTriggers,
        injectionFile );
    if ( vrbflg )
    {
      fprintf( stdout, "Reading search_summary table of %s ...",
          injectionFile );
      fflush( stdout );
    }
    searchSummHead = XLALSearchSummaryTableFromLIGOLw (injectionFile);
    if ( vrbflg ) fprintf( stdout, " done\n");
  }

  if ( ligoOnly )
  {
    ifos[0] = 0;
    ifos[2] = 0;
    ifos[4] = 0;
    ifos[5] = 0;
  }

  /* make sure we start at head of linked list */
  thisInjection = injectionHead;

  chanDummy = XLALCreateREAL4TimeSeries( "", &epoch, f0,
        deltaT, &lalADCCountUnit, numRawPoints );
  if ( !chanDummy ){
    XLALPrintError("failure allocating chanDummy");
    exit(1);
  }

  /*
   *
   * set up the response function
   *
   */
  resp = XLALCreateCOMPLEX8FrequencySeries( chanDummy->name, &chanDummy->epoch,
         f0, deltaF, &strainPerCount, (numRawPoints / 2 + 1) );
  if ( !resp ){
    XLALPrintError("failure allocating response function");
    exit(1);
  }

  /* create vector that will contain detector.transfer info, since this 
   * is constant I calculate it once outside of all the loops and pass it 
   * in to detector.transfer when required 
   */
  detTransDummy = XLALCreateCOMPLEX8FrequencySeries( chanDummy->name,
                  &chanDummy->epoch, f0, deltaF, &strainPerCount,
                  (numRawPoints / 2 + 1) );
  if ( !detTransDummy ){
    XLALPrintError("failure allocating detector.transfer");
    exit(1);
  }

  /* invert the response function to get the transfer function */
  unity = XLALCreateCOMPLEX8Vector( resp->data->length );
  for ( k = 0; k < unity->length; ++k )
  {
    unity->data[k] = 1.0;
  }

  /* set response */
  for ( k = 0; k < resp->data->length; ++k )
  {
    resp->data->data[k] = 1.0;
  }

  XLALCCVectorDivide( detTransDummy->data, unity, resp->data );
  XLALDestroyCOMPLEX8Vector( unity );

  /* setting fixed waveform injection parameters */
  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  ppnParams.deltaT   = deltaT;
  ppnParams.lengthIn = 0;
  ppnParams.ppn      = NULL;

  /* set up resampling parameters */   
  memset( &resampleParams, 0, sizeof(ResampleTSParams) );
  resampleParams.deltaT = deltaTReSample;
  resampleParams.filterType = LDASfirLP;

  /* loop over injections */
  injSimCount = 0;

  do
  {
    fprintf( stdout, "injection %d/%d\n", injSimCount+1, numInjections );

    /* reset waveform structure */
    memset( &waveform, 0, sizeof(CoherentGW) );
    memset( &nonSpinningWaveform, 0, sizeof(CoherentGW) );

    if (thisInjection->f_lower >= fLow)
    {
      fprintf( stderr, "Error: f_lower in sim_inspiral is %e, which is less than fLow.\n", thisInjection->f_lower );
      fprintf( stderr, "Try setting fLow higher with --f-lower option\n");
      exit( 1 );
    }

    /* create the waveform, amp, freq phase etc */
    LAL_CALL( LALGenerateInspiral(&status, &waveform, thisInjection,
          &ppnParams), &status );

    if ( ! snrflg )
    {
      /* create the non-spinning waveform, amp, freq phase etc */
      memcpy( &thisNonSpinningInjection, thisInjection, sizeof(SimInspiralTable) );
      strcpy(thisNonSpinningInjection.waveform, "TaylorT1threePointFivePN\0");
      LAL_CALL( LALGenerateInspiral(&status, &nonSpinningWaveform,
            &thisNonSpinningInjection, &ppnParams), &status);
      if (vrbflg) fprintf( stdout, "ppnParams.tc %e\n ", ppnParams.tc);
    }
    statValue = 0.;

    /* calc lower index for integration */
    kLow = ceil(fLow / deltaF);
    if ( vrbflg )
    {
      fprintf( stdout, "starting integration to find sigmasq at frequency %e ",
          fLow);
      fprintf( stdout, "at index %d \n", kLow);
    }
    /* calc upper index for integration */
    kHi = floor(fReSampling / (2. * deltaF));
    if ( vrbflg )
    {
      fprintf( stdout, "ending integration to find sigmasq at frequency %e ",
          fReSampling / 2.);
      fprintf( stdout, "at index %d \n", kHi);
    }

    /* loop over ifo */
    for ( ifoNumber = 0; ifoNumber < 6; ifoNumber++ )
    {
      if ( ifos[ifoNumber] )
      {
        /* allocate memory and copy the parameters describing the freq series */
        memset( &detector, 0, sizeof( DetectorResponse ) );
        detector.site = (LALDetector *) LALMalloc( sizeof(LALDetector) );

        if (injoverhead)
        { 
          if ( vrbflg ) fprintf( stdout,
              "WARNING: perform overhead injections\n");
          /* setting detector.site to NULL causes SimulateCoherentGW to
           * perform overhead injections */  
          detector.site = NULL; 
        }
        else
        {
          /* if not overhead, set detector.site using ifonumber */  
          XLALReturnDetector( detector.site, ifoNumber );
        } 

        if (vrbflg) fprintf(stdout,
            "generating chan to put waveform in\n" );

        /* get the gps start time of the signal to inject */
        waveformStartTime = XLALGPSToINT8NS( &(thisInjection->geocent_end_time) );
        waveformStartTime -= (INT8) ( 1000000000.0 * ppnParams.tc );

        offset = (3.0 * numRawPoints / 4.0) * deltaT;
        /* test whether the time series is long enough to encompass the
         * injection */
        fprintf(stdout, "offset is %f, length of injection is %f\n", offset, (thisInjection->geocent_end_time.gpsSeconds - ppnParams.tc) );

        if ( offset < ppnParams.tc )
        {
          fprintf(stderr, "The time series is not large enough for this injection.\n" );
          exit (1);
        }
        gpsStartTime.gpsSeconds =
          thisInjection->geocent_end_time.gpsSeconds - (INT4) floor(offset);
        gpsStartTime.gpsNanoSeconds =
          thisInjection->geocent_end_time.gpsNanoSeconds;
        if (vrbflg) fprintf(stdout,
            "offset start time of injection by %f seconds \n", offset ); 

        /* is this okay? copying in detector transfer which so far only
         * contains response info  */
        detector.transfer = detTransDummy;

        XLALUnitInvert( &(detector.transfer->sampleUnits),
            &(resp->sampleUnits) );

        /* set the start times for spinning injection */
        XLALINT8NSToGPS( &(waveform.a->epoch), waveformStartTime );
        memcpy(&(waveform.f->epoch), &(waveform.a->epoch),
            sizeof(LIGOTimeGPS) );
        memcpy(&(waveform.phi->epoch), &(waveform.a->epoch),
            sizeof(LIGOTimeGPS) );

        if ( ! snrflg )
        {
          /* set the start time for the non-spinning injection */
          memcpy(&(nonSpinningWaveform.a->epoch), &(waveform.a->epoch),
              sizeof(LIGOTimeGPS) );
          memcpy(&(nonSpinningWaveform.f->epoch), &(waveform.a->epoch),
              sizeof(LIGOTimeGPS) );
          memcpy(&(nonSpinningWaveform.phi->epoch), &(waveform.a->epoch),
              sizeof(LIGOTimeGPS) );

          chan = XLALCreateREAL4TimeSeries( "", &epoch, f0,
                deltaT, &lalADCCountUnit, numRawPoints );
          if ( !chan ) {
            XLALPrintError("failure allocating chan");
            exit(1);
          }

          /* reset chan structure */
          memcpy( &chan->epoch, &gpsStartTime, sizeof(LIGOTimeGPS) );
          memset( chan->data->data, 0, chan->data->length * sizeof(REAL4) );

          /* perform the non-spinning injection */
          LAL_CALL( LALSimulateCoherentGW(&status, chan, &nonSpinningWaveform,
                &detector ), &status);

          LAL_CALL( LALResampleREAL4TimeSeries( &status, chan,
                &resampleParams ), &status );

          if (writechan)
          {
            /* write out channel data */
            if (vrbflg) fprintf(stdout, "writing channel data to file... \n" );
            switch ( ifoNumber )
            {
              case LAL_IFO_G1:
                snprintf(chanfilename, FILENAME_MAX,
                    "nonspinning_G1_inj%d.dat", injSimCount+1);
                if (vrbflg) fprintf( stdout,
                    "writing G1 channel time series out to %s\n", chanfilename );
                LALSPrintTimeSeries(chan, chanfilename );
                break;
              case LAL_IFO_H1:
                snprintf(chanfilename, FILENAME_MAX,
                    "nonspinning_H1_inj%d.dat", injSimCount+1);
                if (vrbflg) fprintf( stdout,
                    "writing H1 channel time series out to %s\n", chanfilename );
                LALSPrintTimeSeries(chan, chanfilename );
                break;
              case LAL_IFO_H2:
                snprintf(chanfilename, FILENAME_MAX,
                    "nonspinning_H2_inj%d.dat", injSimCount+1);
                if (vrbflg) fprintf( stdout,
                    "writing H2 channel time series out to %s\n", chanfilename );
                LALSPrintTimeSeries(chan, chanfilename );
                break;
              case LAL_IFO_L1:
                snprintf(chanfilename, FILENAME_MAX,
                    "nonspinning_L1_inj%d.dat", injSimCount+1);
                if (vrbflg) fprintf( stdout,
                    "writing L1 channel time series out to %s\n", chanfilename );
                LALSPrintTimeSeries(chan, chanfilename );
                break;
              case LAL_IFO_T1:
                snprintf(chanfilename, FILENAME_MAX,
                    "nonspinning_T1_inj%d.dat", injSimCount+1);
                if (vrbflg) fprintf( stdout,
                    "writing T1 channel time series out to %s\n", chanfilename );
                LALSPrintTimeSeries(chan, chanfilename );
                break;
              case LAL_IFO_V1:
                snprintf(chanfilename, FILENAME_MAX,
                    "nonspinning_V1_inj%d.dat", injSimCount+1);
                if (vrbflg) fprintf( stdout,
                    "writing V1 channel time series out to %s\n", chanfilename );
                LALSPrintTimeSeries(chan, chanfilename );
                break;
              default:
                fprintf( stderr,
                    "Error: ifoNumber %d does not correspond to a known IFO: \n",
                    ifoNumber );
                exit( 1 );
            }
          }

          LAL_CALL( LALCreateForwardRealFFTPlan( &status, &pfwd,
                chan->data->length, 0), &status);
          fftStandardData = XLALCreateCOMPLEX8FrequencySeries( 
                chan->name, &chan->epoch, f0, deltaF, &lalDimensionlessUnit,
                (numPoints / 2 + 1) );
          if ( !fftStandardData ){
            XLALPrintError("failure allocating fftStandardData");
            exit(1);
          }
          LAL_CALL( LALTimeFreqRealFFT( &status, fftStandardData, chan, pfwd ),
              &status);
          LAL_CALL( LALDestroyRealFFTPlan( &status, &pfwd ), &status);
          pfwd = NULL;

          /* reset chan structure */
          XLALDestroyREAL4TimeSeries( chan );
        } /* end if ( ! snrflg ) */

        chan = XLALCreateREAL4TimeSeries( "", &epoch, f0,
              deltaT, &lalADCCountUnit, numRawPoints );
        if ( !chan ){
          XLALPrintError("failure allocating chan");
          exit(1);
        }
        memcpy( &chan->epoch, &gpsStartTime, sizeof(LIGOTimeGPS) );

        /* reset chan structure */
        memset( chan->data->data, 0, chan->data->length * sizeof(REAL4) );

        /* get the gps start time of the signal to inject */
        waveformStartTime = XLALGPSToINT8NS( &(thisInjection->geocent_end_time) );
        waveformStartTime -= (INT8) ( 1000000000.0 * ppnParams.tc );

        offset = (chan->data->length / 2.0) * chan->deltaT;
        gpsStartTime.gpsSeconds =
          thisInjection->geocent_end_time.gpsSeconds - offset;
        gpsStartTime.gpsNanoSeconds =
          thisInjection->geocent_end_time.gpsNanoSeconds;
        chan->epoch = gpsStartTime;

        /* perform the spinning injection */
        LAL_CALL( LALSimulateCoherentGW(&status, chan, &waveform, &detector ),
            &status);
        LAL_CALL( LALResampleREAL4TimeSeries( &status, chan,
              &resampleParams ), &status );

        if (writechan)
        {
          /* write out channel data */
          if (vrbflg) fprintf(stdout, "writing channel data to file... \n" );
          switch ( ifoNumber )
          {
            case LAL_IFO_G1:
              snprintf( chanfilename, FILENAME_MAX, "spinning_G1_inj%d.dat",
                  injSimCount+1);
              if (vrbflg) fprintf( stdout,
                  "writing G1 channel time series out to %s\n", chanfilename );
              LALSPrintTimeSeries(chan, chanfilename );
              break;
            case LAL_IFO_H1:
              snprintf( chanfilename, FILENAME_MAX, "spinning_H1_inj%d.dat",
                  injSimCount+1);
              if (vrbflg) fprintf( stdout,
                  "writing H1 channel time series out to %s\n", chanfilename );
              LALSPrintTimeSeries(chan, chanfilename );
              break;
            case LAL_IFO_H2:
              snprintf( chanfilename, FILENAME_MAX, "spinning_H2_inj%d.dat",
                  injSimCount+1);
              if (vrbflg) fprintf( stdout,
                  "writing H2 channel time series out to %s\n", chanfilename );
              LALSPrintTimeSeries(chan, chanfilename );
              break;
            case LAL_IFO_L1:
              snprintf( chanfilename, FILENAME_MAX, "spinning_L1_inj%d.dat",
                  injSimCount+1);
              if (vrbflg) fprintf( stdout,
                  "writing L1 channel time series out to %s\n", chanfilename );
              LALSPrintTimeSeries(chan, chanfilename );
              break;
            case LAL_IFO_T1:
              snprintf( chanfilename, FILENAME_MAX, "spinning_T1_inj%d.dat",
                  injSimCount+1);
              if (vrbflg) fprintf( stdout,
                  "writing T1 channel time series out to %s\n", chanfilename );
              LALSPrintTimeSeries(chan, chanfilename );
              break;
            case LAL_IFO_V1:
              snprintf( chanfilename, FILENAME_MAX, "spinning_V1_inj%d.dat",
                  injSimCount+1);
              if (vrbflg) fprintf( stdout,
                  "writing V1 channel time series out to %s\n", chanfilename );
              LALSPrintTimeSeries(chan, chanfilename );
              break;
            default:
              fprintf( stderr,
                  "Error: ifoNumber %d does not correspond to a known IFO: \n",
                  ifoNumber );
              exit( 1 );
          }
        }

        LAL_CALL( LALCreateForwardRealFFTPlan( &status, &pfwd,
              chan->data->length, 0), &status);
        fftData = XLALCreateCOMPLEX8FrequencySeries(
              chan->name, &chan->epoch, f0, deltaF, &lalDimensionlessUnit,
              (numPoints / 2 + 1) );
        if ( !fftData ){
          XLALPrintError("failure allocating fftData");
          exit(1);
        }
        LAL_CALL( LALTimeFreqRealFFT( &status, fftData, chan, pfwd ), &status);
        LAL_CALL( LALDestroyRealFFTPlan( &status, &pfwd ), &status);
        pfwd = NULL;

        if ( !snrflg )
        {
          /* compute the Standard Sigmasq */
          thisStandardSigmasq = 0;

          /* avoid f=0 part of psd */
          {
            if (vrbflg)
            {
              switch ( ifoNumber )
              {
                case LAL_IFO_G1: fprintf( stdout, "using GEO PSD \n"); break;
                case LAL_IFO_H1:
                                 fprintf( stdout, "using LIGOI PSD with Hanford Location \n");
                                 break;
                case LAL_IFO_H2:
                                 fprintf( stdout, "using LIGOI PSD with Hanford Location \n");
                                 break;
                case LAL_IFO_L1:
                                 fprintf( stdout, "using LIGOI PSD with Livingston Location \n");
                                 break;
                case LAL_IFO_T1: fprintf( stdout, "using TAMA PSD \n"); break;
                case LAL_IFO_V1: fprintf( stdout, "using VIRGO PSD \n"); break;
                default:
                                 fprintf( stderr,
                                     "Error: ifoNumber %d does not correspond to a known IFO: \n",
                                     ifoNumber );
                                 exit( 1 );

              }
            }
            for ( k = kLow; k < kHi; k++ )
            {
              REAL8 freq;
              REAL8 sim_psd_value;
              freq = fftStandardData->deltaF * k;
              switch( ifoNumber )
              {
                case LAL_IFO_G1: LALGEOPsd( NULL, &sim_psd_value, freq ); break;
                case LAL_IFO_H1: LALLIGOIPsd( NULL, &sim_psd_value, freq ); break;
                case LAL_IFO_H2: LALLIGOIPsd( NULL, &sim_psd_value, freq ); break;
                case LAL_IFO_L1: LALLIGOIPsd( NULL, &sim_psd_value, freq ); break;
                case LAL_IFO_T1: LALTAMAPsd( NULL, &sim_psd_value, freq ); break;
                case LAL_IFO_V1: LALVIRGOPsd( NULL, &sim_psd_value, freq ); break;
                default:
                                 fprintf( stderr,
                                     "Error: ifoNumber %d does not correspond to a known IFO: \n",
                                     ifoNumber );
                                 exit( 1 );
              }

              thisStandardSigmasq +=
                ((crealf(fftStandardData->data->data[k]) * dynRange) *
                 (crealf(fftStandardData->data->data[k]) * dynRange)) /
                sim_psd_value;
              thisStandardSigmasq +=
                ((cimagf(fftStandardData->data->data[k]) * dynRange) *
                 (cimagf(fftStandardData->data->data[k]) * dynRange)) /
                sim_psd_value;
            }
          }

          thisStandardSigmasq *= 4*fftStandardData->deltaF;
          standardSigmasqVec[ifoNumber] = thisStandardSigmasq;
          if ( vrbflg )
          {
            fprintf( stdout, "thisStandardSigmasq %e\n", thisStandardSigmasq );
            fprintf( stdout, "standardSigmasqVec  %e\n",
                standardSigmasqVec[ifoNumber] );
            fflush( stdout );
          }


          /* compute the Mixed Sigmasq */
          thisMixedSigmasq = 0;

          /* avoid f=0 part of psd */
          {
            if (vrbflg)
            {
              switch ( ifoNumber )
              {
                case LAL_IFO_G1: fprintf( stdout, "using GEO PSD \n"); break;
                case LAL_IFO_H1:
                                 fprintf( stdout, "using LIGOI PSD with Hanford Location \n");
                                 break;
                case LAL_IFO_H2:
                                 fprintf( stdout, "using LIGOI PSD with Hanford Location \n");
                                 break;
                case LAL_IFO_L1:
                                 fprintf( stdout, "using LIGOI PSD with Livingston Location \n");
                                 break;
                case LAL_IFO_T1: fprintf( stdout, "using TAMA PSD \n"); break;
                case LAL_IFO_V1: fprintf( stdout, "using VIRGO PSD \n"); break;
                default:
                                 fprintf( stderr,
                                     "Error: ifoNumber %d does not correspond to a known IFO: \n",
                                     ifoNumber );
                                 exit( 1 );

              }
            }
            for ( k = kLow; k < kHi; k++ )
            {
              REAL8 numerator = 0.0;
              REAL8 freq;
              REAL8 sim_psd_value;
              freq = fftData->deltaF * k;
              switch( ifoNumber )
              { 
                case LAL_IFO_G1: LALGEOPsd( NULL, &sim_psd_value, freq ); break;
                case LAL_IFO_H1: LALLIGOIPsd( NULL, &sim_psd_value, freq ); break;
                case LAL_IFO_H2: LALLIGOIPsd( NULL, &sim_psd_value, freq ); break;
                case LAL_IFO_L1: LALLIGOIPsd( NULL, &sim_psd_value, freq ); break;
                case LAL_IFO_T1: LALTAMAPsd( NULL, &sim_psd_value, freq ); break;
                case LAL_IFO_V1: LALVIRGOPsd( NULL, &sim_psd_value, freq ); break;
                default:
                                 fprintf( stderr,
                                     "Error: ifoNumber %d does not correspond to a known IFO: \n",
                                     ifoNumber );
                                 exit( 1 );
              }

              numerator += pow((crealf(fftStandardData->data->data[k]) * dynRange) *
                  (crealf(fftData->data->data[k]) * dynRange) +
                  (cimagf(fftStandardData->data->data[k]) * dynRange) *
                  (cimagf(fftData->data->data[k]) * dynRange),2.0);
              numerator += pow((cimagf(fftStandardData->data->data[k]) * dynRange) *
                  (crealf(fftData->data->data[k]) * dynRange) -
                  (crealf(fftStandardData->data->data[k]) * dynRange) *
                  (cimagf(fftData->data->data[k]) * dynRange),2.0);

              thisMixedSigmasq += pow(numerator,0.5) / sim_psd_value;
            }
          }

          thisMixedSigmasq *= 4*fftData->deltaF;
          mixedSigmasqVec[ifoNumber] = thisMixedSigmasq;

          if ( vrbflg )
          {
            fprintf( stdout, "thisMixedSigmasq %e\n", thisMixedSigmasq );
            fprintf( stdout, "mixedSigmasqVec  %e\n",
                mixedSigmasqVec[ifoNumber] );
            fflush( stdout );
          }
        } /* end if (! snrflg ) */

        /* compute the Spinning sigmasq */
        thisSigmasq = 0;

        /* avoid f=0 part of psd */  
        {
          if (vrbflg)
          {
            switch ( ifoNumber )
            {
              case LAL_IFO_G1: fprintf( stdout, "using GEO PSD \n"); break;
              case LAL_IFO_H1:
                               fprintf( stdout, "using LIGOI PSD with Hanford Location \n");
                               break;
              case LAL_IFO_H2:
                               fprintf( stdout, "using LIGOI PSD with Hanford Location \n");
                               break;
              case LAL_IFO_L1:
                               fprintf( stdout, "using LIGOI PSD with Livingston Location \n");
                               break;
              case LAL_IFO_T1: fprintf( stdout, "using TAMA PSD \n"); break;
              case LAL_IFO_V1: fprintf( stdout, "using VIRGO PSD \n"); break;
              default:
                               fprintf( stderr,
                                   "Error: ifoNumber %d does not correspond to a known IFO: \n",
                                   ifoNumber );
                               exit( 1 );

            }
          }
          for ( k = kLow; k < kHi; k++ )
          {
            REAL8 freq;
            REAL8 sim_psd_value;
            freq = fftData->deltaF * k;
            switch( ifoNumber )
            { 
              case LAL_IFO_G1: LALGEOPsd( NULL, &sim_psd_value, freq ); break;
              case LAL_IFO_H1: LALLIGOIPsd( NULL, &sim_psd_value, freq ); break;
              case LAL_IFO_H2: LALLIGOIPsd( NULL, &sim_psd_value, freq ); break;
              case LAL_IFO_L1: LALLIGOIPsd( NULL, &sim_psd_value, freq ); break;
              case LAL_IFO_T1: LALTAMAPsd( NULL, &sim_psd_value, freq ); break;
              case LAL_IFO_V1: LALVIRGOPsd( NULL, &sim_psd_value, freq ); break;
              default:
                               fprintf( stderr,
                                   "Error: ifoNumber %d does not correspond to a known IFO: \n",
                                   ifoNumber );
                               exit( 1 );
            }

            thisSigmasq +=
              ((crealf(fftData->data->data[k]) * dynRange) * 
               (crealf(fftData->data->data[k]) * dynRange)) /
              sim_psd_value;
            thisSigmasq +=
              ((cimagf(fftData->data->data[k]) * dynRange) * 
               (cimagf(fftData->data->data[k]) * dynRange)) /
              sim_psd_value;
          }
        }

        thisSigmasq *= 4*fftData->deltaF;
        spinningSigmasqVec[ifoNumber] = thisSigmasq; 
        XLALDestroyCOMPLEX8FrequencySeries(fftData);

        if ( ! snrflg )
          XLALDestroyCOMPLEX8FrequencySeries(fftStandardData);

        if ( vrbflg )
        {
          fprintf( stdout, "thisSigmasq        %e\n", thisSigmasq );
          fprintf( stdout, "spinningSigmasqVec %e\n",
              spinningSigmasqVec[ifoNumber] );
          fflush( stdout );
        }

        /* Store StandardSigmasq in Store_ifoSNRsq[] for later test */
        Store_ifoSNRsq[ifoNumber] = thisSigmasq;

        /* free some memory */
        if (detector.transfer) detector.transfer = NULL;
        if ( detector.site )
        {
          LALFree( detector.site);
          detector.site = NULL;
        }
        XLALDestroyREAL4TimeSeries( chan );
      }
    }
    /* end loop over ifo */

    if ( snrflg )
    {
      /* For toy MC, compare StandardSigmasq to threshold value to see if it is 
         found or missed; if found, save parameters to SnglFound and thisFoundSim */
      if ( Store_ifoSNRsq[LAL_IFO_H1] >= SNRsq_threshold && Store_ifoSNRsq[LAL_IFO_L1] >= SNRsq_threshold )
      {
        /* create SimTable */
        if ( ! headFoundSim ) /*first found injection; store as head */
        {
          headFoundSim = thisFoundSim = (SimInspiralTable *) LALCalloc(1, sizeof(SimInspiralTable));
        }
        else
        {
          thisFoundSim = thisFoundSim->next = (SimInspiralTable *) LALCalloc(1, sizeof(SimInspiralTable));
        }
        /* copy SimTable from thisInjection */
        memcpy(thisFoundSim, thisInjection, sizeof(SimInspiralTable));
        thisFoundSim->next = NULL;
        /* create and populate SnglTable */
        if ( ! headFoundSngl )
        {
          headFoundSngl = H1FoundSngl = (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
        }
        else
        {
          H1FoundSngl = L1FoundSngl->next = (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
        }
        /* H1 entry */
        snprintf( H1FoundSngl->ifo, LIGOMETA_IFO_MAX, "H1" );
        snprintf( H1FoundSngl->channel, LIGOMETA_CHANNEL_MAX, "LSC-STRAIN" );
        snprintf( H1FoundSngl->search, LIGOMETA_SEARCH_MAX, "sned" );
        H1FoundSngl->eff_distance = thisInjection->eff_dist_h; 
        H1FoundSngl->end_time = thisInjection->h_end_time;
        H1FoundSngl->mass1 = thisInjection->mass1;
        H1FoundSngl->mass2 = thisInjection->mass2;
        H1FoundSngl->mtotal = thisInjection->mass1 + thisInjection->mass2;
        H1FoundSngl->eta = thisInjection->eta;
        H1FoundSngl->snr = sqrt(thisSigmasq);
        H1FoundSngl->sigmasq = thisSigmasq;
        H1FoundSngl->event_id = (EventIDColumn *)LALCalloc(1, sizeof(EventIDColumn) );
        H1FoundSngl->event_id->id = thisFoundSim->geocent_end_time.gpsSeconds * 1000000000LL + eventNumber++;
        /*storing injection number to alpha column (this means INT4 goes to REAL4)*/
        H1FoundSngl->alpha = injSimCount + 1; 
        /* L1 entry */
        L1FoundSngl = H1FoundSngl->next = (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
        memcpy( L1FoundSngl, H1FoundSngl, sizeof(SnglInspiralTable) );
        L1FoundSngl->next = NULL;
        snprintf( L1FoundSngl->ifo, LIGOMETA_IFO_MAX, "L1" );
        L1FoundSngl->eff_distance = thisInjection->eff_dist_l; 
        L1FoundSngl->end_time = thisInjection->l_end_time;
      }
      else /* Missed; save to thisMissedSim */
      {
        if ( ! thisMissedSim)
        {
          headMissedSim = thisMissedSim = (SimInspiralTable *) LALCalloc(1, sizeof(SimInspiralTable));
        }
        else
        {
          thisMissedSim = thisMissedSim->next = (SimInspiralTable *) LALCalloc(1, sizeof(SimInspiralTable));
        }
        memcpy(thisMissedSim, thisInjection, sizeof(SimInspiralTable));
        thisMissedSim->next = NULL;
      }
    }

    destroyCoherentGW( &waveform );
    destroyCoherentGW( &nonSpinningWaveform );
    if ( !snrflg )
    {
      /* normalize the eff_dist columns */
      if ( nonSpinningSearch )
      {
        thisInjection->eff_dist_g *=
          standardSigmasqVec[LAL_IFO_G1]/mixedSigmasqVec[LAL_IFO_G1];
        thisInjection->eff_dist_h *=
          standardSigmasqVec[LAL_IFO_H1]/mixedSigmasqVec[LAL_IFO_H1];
        thisInjection->eff_dist_l *=
          standardSigmasqVec[LAL_IFO_L1]/mixedSigmasqVec[LAL_IFO_L1];
        thisInjection->eff_dist_t *=
          standardSigmasqVec[LAL_IFO_T1]/mixedSigmasqVec[LAL_IFO_T1];
        thisInjection->eff_dist_v *=
          standardSigmasqVec[LAL_IFO_V1]/mixedSigmasqVec[LAL_IFO_V1];
      }
      else
      {
        thisInjection->eff_dist_g *=
          pow( standardSigmasqVec[LAL_IFO_G1]/spinningSigmasqVec[LAL_IFO_G1],
              0.5 );
        thisInjection->eff_dist_h *=
          pow( standardSigmasqVec[LAL_IFO_H1]/spinningSigmasqVec[LAL_IFO_H1],
              0.5 );
        thisInjection->eff_dist_l *=
          pow( standardSigmasqVec[LAL_IFO_L1]/spinningSigmasqVec[LAL_IFO_L1],
              0.5 );
        thisInjection->eff_dist_t *=
          pow( standardSigmasqVec[LAL_IFO_T1]/spinningSigmasqVec[LAL_IFO_T1],
              0.5 );
        thisInjection->eff_dist_v *=
          pow( standardSigmasqVec[LAL_IFO_V1]/spinningSigmasqVec[LAL_IFO_V1],
              0.5 );
      }
    } /* endif */
    /* increment the bank sim sim_inspiral table if necessary */
    if ( injectionHead )
    {
      thisInjection = thisInjection->next;
    }
  } while ( ++injSimCount < numInjections ); 
  /* end loop over injections */

  /* try opening, writing and closing the xml files */
  if ( snrflg )
  {
    /* open and write to the Found xml file */
    if ( vrbflg ) fprintf( stdout, "Opening Found xml file... " );
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
    snprintf( fname, sizeof(fname), "Found_Injections.xml" );
    LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlStream, fname ), &status );

    /* write out the process and process params tables */
    XLALGPSTimeNow(&(proctable.processTable->end_time));
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ),
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, proctable,
          process_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status );

    /* write the process params table */
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream,
          process_params_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, procparams,
          process_params_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status );

    /* write the Found SnglTable */
    if ( vrbflg ) fprintf( stdout, "writing SnglTable... ");
    outputTable.snglInspiralTable = headFoundSngl;
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream,
          sngl_inspiral_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable,
          sngl_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status );

    /* write the Found SimTable */
    if ( vrbflg ) fprintf( stdout, "writing SimTable... ");
    outputTable.simInspiralTable = headFoundSim;
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, sim_inspiral_table ),
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable,
          sim_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status );

    /* close the Found xml file */ 
    if ( vrbflg ) fprintf( stdout, " Closing\n");
    LAL_CALL( LALCloseLIGOLwXMLFile( &status, &xmlStream ), &status );

    /* open and write to the Missed xml file */
    if ( vrbflg ) fprintf( stdout, "Opening Missed xml file... ");
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
    snprintf( fname, sizeof(fname), "Missed_Injections.xml" );
    LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlStream, fname ), &status );

    /* write out the process and process params tables */
    XLALGPSTimeNow(&(proctable.processTable->end_time));
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ),
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, proctable,
          process_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status );

    /* write the process params table */
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream,
          process_params_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, procparams,
          process_params_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status );

    /* write the Missed SimTable */
    if ( vrbflg ) fprintf( stdout, "writing SimTable... " );
    outputTable.simInspiralTable = headMissedSim;
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, sim_inspiral_table ),
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable,
          sim_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status );

    /* close the Missed xml file */ 
    if ( vrbflg ) fprintf( stdout, "Closing\n");
    LAL_CALL( LALCloseLIGOLwXMLFile( &status, &xmlStream ), &status );
  }
  else
  {
    /* open the output xml file */
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
    snprintf( fname, sizeof(fname), "%s", outputFile );
    LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlStream, fname ), &status );

    /* write out the process and process params tables */
    if ( vrbflg ) fprintf( stdout, "process... " );
    XLALGPSTimeNow(&(proctable.processTable->end_time));
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ),
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, proctable,
          process_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status );
    free( proctable.processTable );
    /* Just being pedantic here ... */
    proctable.processTable = NULL;

    /* free the unused process param entry */
    this_proc_param = procparams.processParamsTable;
    procparams.processParamsTable = procparams.processParamsTable->next;
    free( this_proc_param );
    this_proc_param = NULL;

    /* write the process params table */
    if ( vrbflg ) fprintf( stdout, "process_params... " );
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream,
          process_params_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, procparams,
          process_params_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status );

    /* write the search summary table */
    if ( coireflg )
    {
      if ( vrbflg ) fprintf( stdout, "search_summary... " );
      outputTable.searchSummaryTable = searchSummHead;
      LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream,
            search_summary_table ), &status );
      LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable,
            search_summary_table ), &status );
      LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status );
    }

    /* write the sim inspiral table */
    if ( vrbflg ) fprintf( stdout, "sim_inspiral... " );
    outputTable.simInspiralTable = injectionHead;
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, sim_inspiral_table ),
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable,
          sim_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status );

    /* write the sngl inspiral table */
    if ( coireflg )
    {
      if ( vrbflg ) fprintf( stdout, "sngl_inspiral... " );
      outputTable.snglInspiralTable = snglHead;
      LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream,
            sngl_inspiral_table ), &status );
      LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable,
            sngl_inspiral_table ), &status );
      LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status );
    } 

    /* close the xml file */ 
    LAL_CALL( LALCloseLIGOLwXMLFile( &status, &xmlStream ), &status );
  } /* endif */
  free( injectionFile ); 
  injectionFile = NULL;



  /* free the process params */
  while( procparams.processParamsTable )
  {
    this_proc_param = procparams.processParamsTable;
    procparams.processParamsTable = this_proc_param->next;
    free( this_proc_param );
    this_proc_param = NULL;
  }

  /* free the sim inspiral tables */
  while ( injectionHead )
  {
    thisInjection = injectionHead;
    injectionHead = injectionHead->next;
    LALFree( thisInjection );
  }
  if ( snrflg )
  {
    while ( headFoundSim )
    {
      thisFoundSim = headFoundSim;
      headFoundSim = headFoundSim->next;
      LALFree( thisFoundSim );
    }

    while ( headMissedSim )
    {
      thisMissedSim = headMissedSim;
      headMissedSim = headMissedSim->next;
      LALFree( thisMissedSim );
    } 

    for ( H1FoundSngl = headFoundSngl; H1FoundSngl; H1FoundSngl = H1FoundSngl->next )
    {
      H1FoundSngl = H1FoundSngl->next;
      LALFree( H1FoundSngl->event_id );
    }

    while ( headFoundSngl )
    {
      H1FoundSngl = headFoundSngl;
      headFoundSngl = headFoundSngl->next;
      LALFree( H1FoundSngl );
    }
  }
  /* Freeing memory */
  XLALDestroyREAL4TimeSeries( chanDummy );
  XLALDestroyCOMPLEX8FrequencySeries( resp );
  XLALDestroyCOMPLEX8FrequencySeries( detTransDummy );


  /*check for memory leaks */
  LALCheckMemoryLeaks(); 

  exit( 0 ); 
}
