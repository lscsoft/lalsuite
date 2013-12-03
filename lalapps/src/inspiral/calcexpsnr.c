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
 * File Name: calexpsnr.c
 *
 * Author: Jones, G. W.
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
#include <lal/LALFrStream.h>
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
#define PROGRAM_NAME "calcexpsnr"

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
  "lalapps_calcexpsnr [options]\n"\
"\nDefaults are shown in brackets\n\n" \
"  --help                   display this message\n"\
"  --verbose                be verbose\n"\
"  --version                version info\n"\
"  --spectrum-H1 FILE       FILE contains PSD info for H1\n"\
"  --spectrum-H2 FILE       FILE contains PSD info for H2\n"\
"  --spectrum-L1 FILE       FILE contains PSD info for L1\n"\
"  --ligo-srd               use LIGOI SRD curve\n"\
"  --inject-overhead        inject signals from overhead detector\n"\
"  --write-chan             write out time series showing inspiral waveform\n"\
"  --inj-file    FILE       xml FILE contains injections \n"\
"  --coire-flag             use this if inj file is a coire file \n"\
"  --output-file FILE       FILE for output \n"\
"  --f-lower     FREQ       freq at which to begin integration \n"\
"\n"

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

extern int vrbflg;           /* verbocity of lal function    */
int ligosrd;                 /* whether to use ligo srd      */
int writechan;               /* whether to write chan txt files  */
int injoverhead;             /* perform inj overhead if this option is set  */
int coireflg;                /* is input file coire (1) or inj (null) */

int main( int argc, char *argv[] )
{
  LALStatus                     status = blank_status;

  UINT4                         k;
  UINT4                         kLow;
  UINT4                         kHi;
  INT4                          numPoints       = 524288;
  REAL4                         fSampling       = 2048.;
  REAL4                         fLow            = 70.;
  REAL4                         fLowInj         = 40.;
  REAL8                         deltaT          = 1./fSampling;
  REAL8                         deltaF          = fSampling / numPoints;

  REAL4                          statValue;
 
  /* vars required to make freq series */
  LIGOTimeGPS                   epoch = { 0, 0 };
  LIGOTimeGPS                   gpsStartTime = {0, 0}; 
  REAL8                         f0 = 0.;
  REAL8                         offset = 0.;
  INT8                          waveformStartTime = 0;

  /* files contain PSD info */
  CHAR                         *injectionFile = NULL;         
  CHAR                         *outputFile    = NULL;         
  CHAR                         *specFileH1    = NULL;         
  CHAR                         *specFileH2    = NULL;         
  CHAR                         *specFileL1    = NULL;         

  COMPLEX8Vector               *unity = NULL;
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};

  int                           numInjections = 0;
  int                           numTriggers = 0;

  /* template bank simulation variables */
  INT4                         injSimCount = 0;
  SimInspiralTable            *injectionHead  = NULL;
  SimInspiralTable            *thisInjection  = NULL;
  SnglInspiralTable           *snglHead       = NULL;
  SearchSummaryTable          *searchSummHead = NULL;
  /*SummValueTable              *summValueHead  = NULL;    */

  /* raw input data storage */
  REAL8FrequencySeries          *specH1        = NULL;
  REAL8FrequencySeries          *specH2        = NULL;
  REAL8FrequencySeries          *specL1        = NULL;
  REAL8FrequencySeries          *thisSpec      = NULL;
  COMPLEX8FrequencySeries       *resp          = NULL;
  COMPLEX8FrequencySeries       *detTransDummy = NULL;
  REAL4TimeSeries               *chan          = NULL;
  RealFFTPlan                   *pfwd          = NULL;
  COMPLEX8FrequencySeries       *fftData       = NULL;
  REAL8                          thisSnrsq     = 0;
  REAL8                          thisSnr       = 0;
  REAL8                          thisCombSnr   = 0;
  REAL8                          snrVec[3];
  REAL8                          dynRange      = 1./(3.0e-23);

  /* needed for inj */
  CoherentGW                 waveform;
  PPNParamStruc              ppnParams;
  DetectorResponse           detector;
  InterferometerNumber       ifoNumber   = LAL_UNKNOWN_IFO;

  /* output data */
  LIGOLwXMLStream       xmlStream;
  MetadataTable         proctable;
  MetadataTable         outputTable;
  MetadataTable         procparams;
  CHAR                  fname[256];         
  CHAR                  comment[LIGOMETA_COMMENT_MAX];
  ProcessParamsTable   *this_proc_param = NULL;

  CHAR   chanfilename[FILENAME_MAX];

  REAL4 sum = 0;
  REAL4 bitten_H1 = 0;
  REAL4 bitten_H2 = 0;
  REAL4 thisCombSnr_H1H2 = 0;

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
    {"help",                    no_argument,       0,                'h'},
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"version",                 no_argument,       0,                'V'},
    {"spectrum-H1",             required_argument, 0,                'a'},
    {"spectrum-H2",             required_argument, 0,                'b'},
    {"spectrum-L1",             required_argument, 0,                'c'},
    {"inj-file",                required_argument, 0,                'd'},
    {"comment",                 required_argument, 0,                'e'},
    {"output-file",             required_argument, 0,                'f'},
    {"coire-flag",              no_argument,       &coireflg,         1 },
    {"ligo-srd",                no_argument,       &ligosrd,          1 },
    {"write-chan",              no_argument,       &writechan,        1 },
    {"inject-overhead",         no_argument,       &injoverhead,      1 },
    {"f-lower",                 required_argument, 0,                'g'},
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

    c = getopt_long_only( argc, argv, "a:b:c:d:e:f:g:hV", long_options, &option_index );

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

      case 'a':
        /* create storage for the spectrum file name */
        optarg_len = strlen( optarg ) + 1;
        specFileH1 = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( specFileH1, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'b':
        /* create storage for the spectrum file name */
        optarg_len = strlen( optarg ) + 1;
        specFileH2 = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( specFileH2, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'c':
        /* create storage for the spectrum file name */
        optarg_len = strlen( optarg ) + 1;
        specFileL1 = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( specFileL1, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
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
        fLow = (INT4) atof( optarg );
        if ( fLow < 40 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "f-lower must be > 40Hz (%e specified)\n",
              long_options[option_index].name, fLow );
          exit( 1 );
        }
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
        fprintf( stdout, "calculation of expected SNR of injections\n"
            "Gareth Jones\n");
        XLALOutputVersionString(stderr, 0);
        exit( 0 );
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

  /* check the input arguments */
  if ( injectionFile == NULL )
  {
    fprintf( stderr, "Must specify the --injection-file\n" );
    exit( 1 );
  }

  if ( outputFile == NULL )
  {
    fprintf( stderr, "Must specify the --output-file\n" );
    exit( 1 );
  }

  if ( !ligosrd && specFileH1 == NULL )
  {
    fprintf( stderr, "Must specify the --spectrum-H1\n" );
    exit( 1 );
  }

  if ( !ligosrd && specFileH2 == NULL )
  {
    fprintf( stderr, "Must specify the --spectrum-H2\n" );
    exit( 1 );
  }

  if ( !ligosrd && specFileL1 == NULL )
  {
    fprintf( stderr, "Must specify the --spectrum-L1\n" );
    exit( 1 );
  }

  if ( ligosrd && (specFileH1 || specFileH2 || specFileL1 ))
  {
    fprintf( stdout, "WARNING: using LIGOI SRD power spectral density \n" );
  } 
 
  if ( vrbflg ){
    fprintf( stdout, "injection file is %s\n", injectionFile );
    fprintf( stdout, "output file is %s\n", outputFile );
    fprintf( stdout, "H1 spec file is   %s\n", specFileH1 );
    fprintf( stdout, "H2 spec file is   %s\n", specFileH2 );
    fprintf( stdout, "L1 spec file is   %s\n", specFileL1 );
  }

  /* create vector for H1, H2 and L1 spectrums */
  specH1 = XLALCreateREAL8FrequencySeries ( "",&epoch, f0, deltaF, &lalADCCountUnit, (numPoints / 2 + 1) );
  specH2 = XLALCreateREAL8FrequencySeries ( "",&epoch, f0, deltaF, &lalADCCountUnit, (numPoints / 2 + 1) );
  specL1 = XLALCreateREAL8FrequencySeries ( "",&epoch, f0, deltaF, &lalADCCountUnit, (numPoints / 2 + 1) );
  if (!specH1 || !specH2 || !specL1){
    XLALDestroyREAL8FrequencySeries ( specH1 );
    XLALDestroyREAL8FrequencySeries ( specH2 );
    XLALDestroyREAL8FrequencySeries ( specL1 );
    XLALPrintError("failure allocating H1, H2 and L1 spectra");
    exit(1);
  }

  if (!ligosrd){
    /* read in H1 spectrum */ 
    LAL_CALL( LALDReadFrequencySeries(&status, specH1, specFileH1), &status );
    if ( vrbflg ){
       fprintf( stdout, "read in H1 spec file\n" );
       fflush( stdout );
    } 

    /* read in H2 spectrum */ 
    LAL_CALL( LALDReadFrequencySeries(&status, specH2, specFileH2), &status );
    if ( vrbflg ){
       fprintf( stdout, "read in H2 spec file\n" );
       fflush( stdout );
    }

    /* read in L1 spectrum */ 
    LAL_CALL( LALDReadFrequencySeries(&status, specL1, specFileL1), &status );
    if ( vrbflg ){
       fprintf( stdout, "read in L1 spec file\n" );
       fflush( stdout );
     }
  }

  chan = XLALCreateREAL4TimeSeries( "", &epoch, f0, deltaT, 
                                     &lalADCCountUnit, numPoints );
  if ( !chan ){
    XLALPrintError("failure allocating chan");
    exit(1);
  }

  /*
   *
   * set up the response function
   *
   */
  resp = XLALCreateCOMPLEX8FrequencySeries( chan->name, 
     &chan->epoch, f0, deltaF, &strainPerCount, (numPoints / 2 + 1) );
  if ( !resp ){
    XLALPrintError("failure allocating response function");
    exit(1);
  }

  /* create vector that will contain detector.transfer info, since this 
   * is constant I calculate it once outside of all the loops and pass it 
   * in to detector.transfer when required 
   */
  detTransDummy = XLALCreateCOMPLEX8FrequencySeries( chan->name, &chan->epoch,
                  f0, deltaF, &strainPerCount, (numPoints / 2 + 1) );
  if ( !detTransDummy ){
    XLALPrintError("failure allocating detector.transfer info");
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

  /* read in injections from injection file */
  /* set endtime to 0 so that we read in all events */
  if ( vrbflg ) fprintf( stdout, "Reading sim_inspiral table of %s\n", injectionFile );
  LAL_CALL(numInjections = SimInspiralTableFromLIGOLw( &injectionHead, injectionFile, 0, 0), &status);
  if ( vrbflg ) fprintf( stdout, "Read %d injections from sim_inspiral table of %s\n", 
                                    numInjections, injectionFile );

  if (coireflg){
     if ( vrbflg ) fprintf( stdout, "Reading sngl_inspiral table of %s\n", injectionFile );
     LAL_CALL(numTriggers = LALSnglInspiralTableFromLIGOLw(&snglHead, injectionFile, 0, -1), &status);
     if ( vrbflg ) fprintf( stdout, "Read %d triggers from sngl_inspiral table of %s\n", 
                                    numTriggers, injectionFile );
     if ( vrbflg ) {
           fprintf( stdout, "Reading search_summary table of %s ...", injectionFile );
           fflush( stdout );
           }
     searchSummHead = XLALSearchSummaryTableFromLIGOLw (injectionFile);
     if ( vrbflg ) fprintf( stdout, " done\n");
  }

 /* make sure we start at head of linked list */
 thisInjection = injectionHead;

  /* setting fixed waveform injection parameters */
  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  ppnParams.deltaT   = deltaT;
  ppnParams.lengthIn = 0;
  ppnParams.ppn      = NULL;

  /* loop over injections */
  injSimCount = 0;
    
        
  do
  {
     fprintf( stdout, "injection %d/%d\n", injSimCount+1, numInjections );

     /* reset waveform structure */
     memset( &waveform, 0, sizeof(CoherentGW) );

     /* reset chan structure */
     memset( chan->data->data, 0, chan->data->length * sizeof(REAL4) );

     if (thisInjection->f_lower == 0){
        fprintf( stdout, "WARNING: f_lower in sim_inpiral = 0, ");
        fprintf( stdout, "changing this to %e\n ", fLowInj);
        thisInjection->f_lower = fLowInj;
     }

     /* create the waveform, amp, freq phase etc */
     LAL_CALL( LALGenerateInspiral(&status, &waveform, thisInjection, &ppnParams), &status);
     if (vrbflg) fprintf( stdout, "ppnParams.tc %e\n ", ppnParams.tc);

    statValue = 0.;
  
    /* calc lower index for integration */
    kLow = ceil(fLow / deltaF);
    if ( vrbflg ) {
        fprintf( stdout, "starting integration to find SNR at frequency %e ", fLow);
        fprintf( stdout, "at index %d \n", kLow);
    }
    /* calc upper index for integration */
    kHi = floor(fSampling / (2. * deltaF));
    if ( vrbflg ) {
        fprintf( stdout, "ending integration to find SNR at frequency %e ", fSampling / 2.);
        fprintf( stdout, "at index %d \n", kHi);
    }

    /* loop over ifo */
    for ( ifoNumber = 1; ifoNumber < 4; ifoNumber++ )
    {
        /* allocate memory and copy the parameters describing the freq series */
        memset( &detector, 0, sizeof( DetectorResponse ) );
        detector.site = (LALDetector *) LALMalloc( sizeof(LALDetector) );

        if (injoverhead){ 
           if ( vrbflg ) fprintf( stdout, "WARNING: perform overhead injections\n");
           /* setting detector.site to NULL causes SimulateCoherentGW to
            * perform overhead injections */  
           detector.site = NULL; 
        }
        else {
           /* if not overhead, set detector.site using ifonumber */  
           XLALReturnDetector( detector.site, ifoNumber );
        } 

        switch ( ifoNumber )
        {
        case 1:
           if ( vrbflg ) fprintf( stdout, "looking at H1 \n");
           thisSpec = specH1;
           break;
        case 2:
           if ( vrbflg ) fprintf( stdout, "looking at H2 \n");
           thisSpec = specH2;
           break;
        case 3:
           if ( vrbflg ) fprintf( stdout, "looking at L1 \n");
           thisSpec = specL1;
           break;
        default:
           fprintf( stderr, "Error: ifoNumber %d does not correspond to H1, H2 or L1: \n", ifoNumber );
           exit( 1 );
        }

        /* get the gps start time of the signal to inject */
        waveformStartTime = XLALGPSToINT8NS( &(thisInjection->geocent_end_time) );
        waveformStartTime -= (INT8) ( 1000000000.0 * ppnParams.tc );

        offset = (chan->data->length / 2.0) * chan->deltaT;
        gpsStartTime.gpsSeconds     = thisInjection->geocent_end_time.gpsSeconds - offset;
        gpsStartTime.gpsNanoSeconds = thisInjection->geocent_end_time.gpsNanoSeconds;
        chan->epoch = gpsStartTime;


       if (vrbflg) fprintf(stdout, "offset start time of injection by %f seconds \n", offset ); 
       
       /* is this okay? copying in detector transfer which so far only contains response info  */
       detector.transfer = detTransDummy;

       XLALUnitInvert( &(detector.transfer->sampleUnits), &(resp->sampleUnits) );

       /* set the start times for injection */
       XLALINT8NSToGPS( &(waveform.a->epoch), waveformStartTime );
       memcpy(&(waveform.f->epoch), &(waveform.a->epoch), sizeof(LIGOTimeGPS) );
       memcpy(&(waveform.phi->epoch), &(waveform.a->epoch), sizeof(LIGOTimeGPS) );
 
       /* perform the injection */
       LAL_CALL( LALSimulateCoherentGW(&status, chan, &waveform, &detector ), &status); 

       if (writechan){ 
          /* write out channel data */
          if (vrbflg) fprintf(stdout, "writing channel data to file... \n" ); 
          switch ( ifoNumber )
          {
          case 1:
             snprintf( chanfilename, FILENAME_MAX, "chanTest_H1_inj%d.dat", injSimCount+1);
             if (vrbflg) fprintf( stdout, "writing H1 channel time series out to %s\n", chanfilename );
             LALSPrintTimeSeries(chan, chanfilename );
             break;
          case 2:
             snprintf( chanfilename, FILENAME_MAX, "chanTest_H2_inj%d.dat", injSimCount+1);
             if (vrbflg) fprintf( stdout, "writing H2 channel time series out to %s\n", chanfilename );
             LALSPrintTimeSeries(chan, chanfilename );
             break;
          case 3:
             snprintf( chanfilename, FILENAME_MAX, "chanTest_L1_inj%d.dat", injSimCount+1);
             if (vrbflg) fprintf( stdout, "writing L1 channel time series out to %s\n", chanfilename );
             LALSPrintTimeSeries(chan, chanfilename );
             break;
         default:
             fprintf( stderr, "Error: ifoNumber %d does not correspond to H1, H2 or L1: \n", ifoNumber );
             exit( 1 );
         }  
      } 

      LAL_CALL( LALCreateForwardRealFFTPlan( &status, &pfwd, chan->data->length, 0), &status);

      fftData = XLALCreateCOMPLEX8FrequencySeries( chan->name, &chan->epoch, f0, deltaF, 
                                                   &lalDimensionlessUnit, (numPoints / 2 + 1) );
      if ( !fftData ){
        XLALPrintError("failure allocating fftData");
        exit(1);
      }
   
      LAL_CALL( LALTimeFreqRealFFT( &status, fftData, chan, pfwd ), &status);
   
      LAL_CALL( LALDestroyRealFFTPlan( &status, &pfwd ), &status);
      pfwd = NULL;

       /* compute the SNR */
       thisSnrsq = 0;
       /* avoid f=0 part of psd */  

       if (ligosrd){
          if (vrbflg) fprintf( stdout, "using LIGOI PSD \n");
          for ( k = kLow; k < kHi; k++ )
          {
           REAL8 freq;
           REAL8 sim_psd_value;
           freq = fftData->deltaF * k;
           LALLIGOIPsd( NULL, &sim_psd_value, freq ); 

           thisSnrsq += ((crealf(fftData->data->data[k]) * dynRange) * 
                      (crealf(fftData->data->data[k]) * dynRange)) / sim_psd_value;
           thisSnrsq += ((cimagf(fftData->data->data[k]) * dynRange) * 
                      (cimagf(fftData->data->data[k]) * dynRange)) / sim_psd_value;
           }
       }
       else {
          if (vrbflg) fprintf( stdout, "using input spectra \n");
          for ( k = kLow; k < kHi; k++ )
          {
           thisSnrsq += ((crealf(fftData->data->data[k]) * dynRange) * 
              (crealf(fftData->data->data[k]) * dynRange))  /
              (thisSpec->data->data[k] * dynRange * dynRange);
           thisSnrsq += ((cimagf(fftData->data->data[k]) * dynRange) * 
              (cimagf(fftData->data->data[k]) * dynRange)) /
              (thisSpec->data->data[k] * dynRange * dynRange);
        } 
      }

       thisSnrsq *= 4*fftData->deltaF;
       thisSnr    = pow(thisSnrsq, 0.5);
       /* Note indexing on snrVec, ifoNumber runs from 1..3 to get source correct,
        * we must index snrVec 0..2 
        */ 
       snrVec[ifoNumber-1] = thisSnr; 
       XLALDestroyCOMPLEX8FrequencySeries(fftData);

       if ( vrbflg ){
          fprintf( stdout, "thisSnrsq %e\n", thisSnrsq );
          fprintf( stdout, "snrVec    %e\n", snrVec[ifoNumber-1] );
          fflush( stdout );
       }

       /* sum thisSnrsq to eventually get combined snr*/
       statValue += thisSnrsq; 

       /* free some memory */
       if (detector.transfer) detector.transfer = NULL;
       if ( detector.site ) {LALFree( detector.site); detector.site = NULL;}
     }
     /* end loop over ifo */
  
    destroyCoherentGW( &waveform );

    /* store inverse eff snrs in eff_dist columns */
    thisInjection->eff_dist_h = 1./snrVec[0];
    thisInjection->eff_dist_g = 1./snrVec[1];
    thisInjection->eff_dist_l = 1./snrVec[2];

    /* store inverse sum of squares snr in eff_dist_t */
    thisCombSnr = pow(statValue, 0.5);
    if ( vrbflg ) fprintf( stdout, "thisCombSnr %e\n", thisCombSnr);
    thisInjection->eff_dist_t = 1./thisCombSnr;

    /* calc inverse bittenL snr for H1H2 and store in eff_dist_v */
    thisCombSnr_H1H2 = 0.;
    sum = snrVec[0] * snrVec[0] + snrVec[1] * snrVec[1];
    bitten_H1 = 3 * snrVec[0] -3;
    bitten_H2 = 3 * snrVec[1] -3;

    if (sum < bitten_H1){
       thisCombSnr_H1H2 = sum;
    }
    else
    {
       thisCombSnr_H1H2 = bitten_H1;
    }

    if (bitten_H2 < thisCombSnr_H1H2){
       thisCombSnr_H1H2 = bitten_H2;
    }
    thisInjection->eff_dist_v = 1./thisCombSnr_H1H2;


    /* increment the bank sim sim_inspiral table if necessary */
    if ( injectionHead )
    {
      thisInjection = thisInjection->next;
    }

  } while ( ++injSimCount < numInjections ); 
  /* end loop over injections */

  /* try opening, writing and closing an xml file */

  /* open the output xml file */
  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  snprintf( fname, sizeof(fname), "%s", outputFile);
  LAL_CALL( LALOpenLIGOLwXMLFile  ( &status, &xmlStream, fname), &status);

  /* write out the process and process params tables */
  if ( vrbflg ) fprintf( stdout, "process... " );
  XLALGPSTimeNow(&(proctable.processTable->end_time));
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, proctable, process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
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
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_params_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, procparams, process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

  /* write the search summary table */
  if ( coireflg ){
     if ( vrbflg ) fprintf( stdout, "search_summary... " );
     outputTable.searchSummaryTable = searchSummHead;
     LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, search_summary_table), &status);
     LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable, search_summary_table), &status);
     LAL_CALL( LALEndLIGOLwXMLTable  ( &status, &xmlStream), &status);
   }

  /* write the sim inspiral table */
  if ( vrbflg ) fprintf( stdout, "sim_inspiral... " );
  outputTable.simInspiralTable = injectionHead;
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, sim_inspiral_table), &status);
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable, sim_inspiral_table), &status);
  LAL_CALL( LALEndLIGOLwXMLTable  ( &status, &xmlStream), &status);

  /* write the sngl inspiral table */
  if ( coireflg ){
     if ( vrbflg ) fprintf( stdout, "sngl_inspiral... " );
     outputTable.snglInspiralTable = snglHead;
     LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, sngl_inspiral_table), &status);
     LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable, sngl_inspiral_table), &status);
     LAL_CALL( LALEndLIGOLwXMLTable  ( &status, &xmlStream), &status);
  } 

  /* close the xml file */ 
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlStream), &status);

  /* Freeing memory */
  XLALDestroyREAL4TimeSeries(chan);
  XLALDestroyCOMPLEX8FrequencySeries(resp);
  XLALDestroyCOMPLEX8FrequencySeries(detTransDummy);
  XLALDestroyREAL8FrequencySeries ( specH1 );
  XLALDestroyREAL8FrequencySeries ( specH2 );
  XLALDestroyREAL8FrequencySeries ( specL1 );


  free( specFileH1 );
  specFileH1 = NULL;
  free( specFileH2 );
  specFileH2 = NULL;
  free( specFileL1 );
  specFileL1 = NULL;
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

  /*check for memory leaks */
  LALCheckMemoryLeaks(); 

  exit( 0 ); 
}
