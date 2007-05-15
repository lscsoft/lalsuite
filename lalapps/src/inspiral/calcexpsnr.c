/*-----------------------------------------------------------------------
 *
 * File Name: calexpsnr.c
 *
 * Author: Jones, G. W.
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
#include <lal/LALTrigScanCluster.h>
#include <lal/PrintFTSeries.h>
#include <lal/ReadFTSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/GenerateInspiral.h>
#include <lal/TimeSeries.h>
#include <lal/VectorOps.h>

RCSID( "$Id$" );

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "calcexpsnr"

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
calloc( 1, sizeof(ProcessParamsTable) ); \
LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
    PROGRAM_NAME ); \
LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
    long_options[option_index].name ); \
LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

#define MAX_PATH 4096


#define USAGE \
  "lalapps_calcexpsnr [options]\n"\
"\nDefaults are shown in brackets\n\n" \
"  --help                   display this message\n"\
"  --verbose                be verbose\n"\
"  --version                version info\n"\
"  --debug-level LEVEL      set the LAL debug level to LEVEL\n"\
"  --spectrum-H1 FILE       FILE contains PSD info for H1\n"\
"  --spectrum-H2 FILE       FILE contains PSD info for H2\n"\
"  --spectrum-L1 FILE       FILE contains PSD info for L1\n"\
"  --inj-file    FILE       xml FILE contains injections\n"\
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

int main( int argc, char *argv[] )
{
  LALStatus                     status = blank_status;
  REAL8                         dynRange;
  REAL4                         dynRangeExponent = 65.;

  UINT4                         k;
  UINT4                         kLow;
  UINT4                         kHi;
  INT4                          numPoints       = 524288;
  REAL4                         fSampling       = 2048.;
  REAL4                         fLow            = 70.;
  REAL8                         deltaT          = 1./fSampling;
  REAL8                         deltaF          = fSampling / numPoints;

  REAL4                          statValue;
 
  /* vars required to make freq series */
  LIGOTimeGPS                   epoch = { 0, 0 }; 
  REAL8                         f0 = 0.;

  /* files contain PSD info */
  CHAR                         *injectionFile = NULL;         
  CHAR                         *specFileH1    = NULL;         
  CHAR                         *specFileH2    = NULL;         
  CHAR                         *specFileL1    = NULL;         

  COMPLEX8Vector               *unity = NULL;
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
  INT8                          waveformStartTime;

  int                           numInjections = 0;

  /* template bank simulation variables */
  INT4                         injSimCount = 0;
  SimInspiralTable            *injectionHead = NULL;
  SimInspiralTable            *thisInjection = NULL;

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
  REAL4                          thisSnrsq     = 0;
  REAL4                          thisSnr       = 0;
  REAL4                          thisCombSnr   = 0;
  REAL4                          snrVec[3];

  /* needed for inj */
  CoherentGW                 waveform;
  PPNParamStruc              ppnParams;
  DetectorResponse           detector;
  InterferometerNumber       ifoNumber   = LAL_UNKNOWN_IFO;

  /* output data */
  LIGOLwXMLStream       xmlStream;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;
  MetadataTable         proctable;
  MetadataTable         outputTable;
  MetadataTable         procparams;
  CHAR                  fname[256];         
  CHAR                  comment[LIGOMETA_COMMENT_MAX];
  ProcessParamsTable   *this_proc_param = NULL;

  /* set initial debug level */
  set_debug_level("1");

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  LAL_CALL(LALGPSTimeNow(&status, &(proctable.processTable->start_time), 
                                                         &accuracy), &status);
  LAL_CALL( populate_process_table( &status, proctable.processTable, 
                PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &status );
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
    {"debug-level",             required_argument, 0,                'z'}, 
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

    c = getopt_long_only( argc, argv, "a:b:c:d:e:z:hV", long_options, &option_index );

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
          LALSnprintf( comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
        }
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "calculation of expected SNR of injections\n"
            "Gareth Jones\n"
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n" );
        exit( 0 );
        break;

      case 'z':
        set_debug_level( optarg );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
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

  if ( specFileH1 == NULL )
  {
    fprintf( stderr, "Must specify the --spectrum-H1\n" );
    exit( 1 );
  }

  if ( specFileH2 == NULL )
  {
    fprintf( stderr, "Must specify the --spectrum-H2\n" );
    exit( 1 );
  }

  if ( specFileL1 == NULL )
  {
    fprintf( stderr, "Must specify the --spectrum-L1\n" );
    exit( 1 );
  }

  
  if ( vrbflg ){
    fprintf( stdout, "injection file is %s\n", injectionFile );
    fprintf( stdout, "H1 spec file is   %s\n", specFileH1 );
    fprintf( stdout, "H2 spec file is   %s\n", specFileH2 );
    fprintf( stdout, "L1 spec file is   %s\n", specFileL1 );
  }

  /* create vector for H1, H2 and L1 spectrums */
  LAL_CALL( LALCreateREAL8FrequencySeries ( &status, &specH1, "",epoch, f0, deltaF, lalADCCountUnit, (numPoints / 2 + 1) ), &status);
  LAL_CALL( LALCreateREAL8FrequencySeries ( &status, &specH2, "",epoch, f0, deltaF, lalADCCountUnit, (numPoints / 2 + 1) ), &status);
  LAL_CALL( LALCreateREAL8FrequencySeries ( &status, &specL1, "",epoch, f0, deltaF, lalADCCountUnit, (numPoints / 2 + 1) ), &status);

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
  
  /* write out spectrum test */
  /*LALDPrintFrequencySeries(specH1, "Test.dat" );
  fprintf( stdout, "written H1 spec file\n" );
  fflush( stdout );
  */

  /* create the dynamic range exponent */
  dynRange = pow( 2.0, dynRangeExponent );
 
  LAL_CALL( LALCreateREAL4TimeSeries( &status, &chan, "", epoch, f0, deltaT, 
                                     lalADCCountUnit, numPoints ), &status );

  /*
   *
   * set up the response function
   *
   */
  LAL_CALL( LALCreateCOMPLEX8FrequencySeries( &status, &resp, chan->name, 
     chan->epoch, f0, deltaF, strainPerCount, (numPoints / 2 + 1) ), &status );

  /* create vector that will contain detector.transfer info, since this 
   * is constant I calculate it once outside of all the loops and pass it 
   * in to detector.transfer when required 
   */
  LAL_CALL( LALCreateCOMPLEX8FrequencySeries( &status, &detTransDummy, 
                  chan->name, chan->epoch, f0, deltaF, strainPerCount, 
                  (numPoints / 2 + 1) ), &status );

  /* invert the response function to get the transfer function */
  unity = XLALCreateCOMPLEX8Vector( resp->data->length );
  for ( k = 0; k < unity->length; ++k )
     {
        unity->data[k].re = 1.0;
        unity->data[k].im = 0.0;
     }

  /* set response */
  for ( k = 0; k < resp->data->length; ++k )
  {
      resp->data->data[k].re = 1.0;
      resp->data->data[k].im = 0.0;
  }

  XLALCCVectorDivide( detTransDummy->data, unity, resp->data );
  XLALDestroyCOMPLEX8Vector( unity );

  LALCPrintFrequencySeries ( detTransDummy, "transTest.dat" );
  LALCPrintFrequencySeries ( resp, "respTest.dat" );
  
  /* read in injections from injection file */
  /* set endtime to 0 so that we read in all events */
  LAL_CALL(numInjections = SimInspiralTableFromLIGOLw( &injectionHead, injectionFile, 0, 0), &status);

  for ( thisInjection = injectionHead; thisInjection; thisInjection = thisInjection->next )
  {
     /* set the time of the injection to zero so it is injected in  */
     /* the middle of the data segment: the other times are ignored */
    thisInjection->geocent_end_time.gpsSeconds = 0;
    thisInjection->geocent_end_time.gpsNanoSeconds = 0;
  }
  thisInjection = injectionHead;

  if ( vrbflg ) fprintf( stdout, "Read %d inj parameters from %s\n", 
                                    numInjections, injectionFile );

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

     /* create the waveform, amp, freq phase etc */
     LAL_CALL( LALGenerateInspiral(&status, &waveform, thisInjection, &ppnParams), &status);

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
        XLALReturnDetector( detector.site, ifoNumber );
 
        fprintf( stdout, "ifoNumber %d\n", ifoNumber );
        fflush( stdout );

        switch ( ifoNumber )
        {
        case 1:
           fprintf( stdout, "looking at H1 \n");
           thisSpec = specH1;
           break;
        case 2:
           fprintf( stdout, "looking at H2 \n");
           thisSpec = specH2;
           break;
        case 3:
           fprintf( stdout, "looking at L1 \n");
           thisSpec = specL1;
           break;
        default:
           fprintf( stderr, "Error: ifoNumber %d does not correspond to H1, H2 or L1: \n", ifoNumber );
           exit( 1 );
        }

        /*LAL_CALL( LALCreateCOMPLEX8FrequencySeries( &status, &detector.transfer, chan->name, chan->epoch, f0, 
         *                                           deltaF, strainPerCount, (numPoints / 2 + 1) ), &status );
         */

        /* is this okay? copying in detector transfer which so far only contains response info  */
        detector.transfer = detTransDummy;

        XLALUnitInvert( &(detector.transfer->sampleUnits), &(resp->sampleUnits) );

        /* from bank sim routines */
        thisInjection->geocent_end_time.gpsSeconds = 0;
        thisInjection->geocent_end_time.gpsNanoSeconds = 0;

        waveformStartTime = XLALGPStoINT8( &(thisInjection->geocent_end_time));
        waveformStartTime -= (INT8) ( 1000000000.0 * ppnParams.tc );
        waveformStartTime += (INT8) ( 1000000000.0 *
          ((REAL8) (chan->data->length - ppnParams.length) / 2.0) * chan->deltaT);
    
        XLALINT8toGPS( &(waveform.a->epoch), waveformStartTime );
        memcpy(&(waveform.f->epoch), &(waveform.a->epoch), sizeof(LIGOTimeGPS) );
        memcpy(&(waveform.phi->epoch), &(waveform.a->epoch), sizeof(LIGOTimeGPS) );

  
       /* perform the injection */
       LAL_CALL( LALSimulateCoherentGW(&status, chan, &waveform, &detector ), &status); 

       /* will this work? No! */
       /*LAL_CALL( LALDestroyCOMPLEX8FrequencySeries( &status, detector.transfer), &status );
        *LALFree( detector.site);
        */
 
      /* write out channel  */
      LALSPrintTimeSeries(chan, "chanTest.dat" );
      fprintf( stdout, "written out chan\n" );
      fflush( stdout );

      LAL_CALL( LALCreateForwardRealFFTPlan( &status, &pfwd, chan->data->length, 0), &status);

      LAL_CALL( LALCreateCOMPLEX8FrequencySeries( &status, &fftData, chan->name, chan->epoch, f0, deltaF, 
                                                   lalDimensionlessUnit, (numPoints / 2 + 1) ), &status );
   
      LAL_CALL( LALTimeFreqRealFFT( &status, fftData, chan, pfwd ), &status);
   
      LAL_CALL( LALDestroyRealFFTPlan( &status, &pfwd ), &status);
      pfwd = NULL;

       /* compute the SNR for initial LIGO at design */
       thisSnrsq = 0;
       /* avoid f=0 part of psd */  
       for ( k = kLow; k < kHi; k++ )
       {
         REAL8 freq;
         /* use correct psd !!!!! */
         /*REAL8 sim_psd_value*/;
         freq = fftData->deltaF * k;
         /* LALLIGOIPsd( NULL, &psd_value, freq ); */

         thisSnrsq += fftData->data->data[k].re * fftData->data->data[k].re /
           thisSpec->data->data[k];
         thisSnrsq += fftData->data->data[k].im * fftData->data->data[k].im /
           thisSpec->data->data[k];
       }
       thisSnrsq *= 4*fftData->deltaF;
       thisSnr    = pow(thisSnrsq, 0.5);
       snrVec[ifoNumber] = thisSnr; 
       LAL_CALL( LALDestroyCOMPLEX8FrequencySeries( &status, fftData), &status );

       fprintf( stdout, "thisSnrsq %e\n", thisSnrsq );
       fprintf( stdout, "thisSnr   %e\n", thisSnr );
       fprintf( stdout, "snrVec    %e\n", snrVec[ifoNumber] );
       fflush( stdout );

       /* sum thisSnrsq to eventually get combined snr*/
       statValue += thisSnrsq; 
     }
     /* end loop over ifo */

    thisCombSnr = pow(statValue, 0.5);
    fprintf( stdout, "thisCombSnr %e\n", thisCombSnr);

    /* increment the bank sim sim_inspiral table if necessary */
    if ( injectionHead )
    {
      thisInjection = thisInjection->next;
    }

  } while ( ++injSimCount < numInjections ); 
  /* end loop over injections */

  fprintf( stdout, "out of loop\n" );
  fflush( stdout );
 
  /* try opening, writing and closing an xml file */

  /* open the output xml file */
  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  LALSnprintf( fname, sizeof(fname), "outFile.xml");
  LAL_CALL( LALOpenLIGOLwXMLFile  ( &status, &xmlStream, fname), &status);

  /* write out the process and process params tables */
  if ( vrbflg ) fprintf( stdout, "process... " );
  LAL_CALL(LALGPSTimeNow(&status, &(proctable.processTable->end_time), &accuracy), &status);
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, proctable, process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
  free( proctable.processTable );
 
  /* free the unused process param entry */
  this_proc_param = procparams.processParamsTable;
  procparams.processParamsTable = procparams.processParamsTable->next;
  free( this_proc_param );

  /* write the process params table */
  if ( vrbflg ) fprintf( stdout, "process_params... " );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_params_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, procparams, process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

  /* write the sim inspiral table */
  if ( vrbflg ) fprintf( stdout, "sim_inspiral... " );
  outputTable.simInspiralTable = injectionHead;
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, sim_inspiral_table), &status);
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable, sim_inspiral_table), &status);
  LAL_CALL( LALEndLIGOLwXMLTable  ( &status, &xmlStream), &status);
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlStream), &status);

  /* Freeing memory */
  LAL_CALL( LALDestroyREAL4TimeSeries( &status, chan), &status );
  LAL_CALL( LALDestroyCOMPLEX8FrequencySeries( &status, resp), &status );
  LAL_CALL( LALDestroyCOMPLEX8FrequencySeries( &status, detTransDummy), &status );
  LAL_CALL( LALDestroyREAL8FrequencySeries ( &status, specH1 ), &status);
  LAL_CALL( LALDestroyREAL8FrequencySeries ( &status, specH2 ), &status);
  LAL_CALL( LALDestroyREAL8FrequencySeries ( &status, specL1 ), &status);

  free( specFileH1 );
  free( specFileH2 );
  free( specFileL1 );
  free( injectionFile ); 

  fprintf( stdout, "111111\n" );
  fflush( stdout );

  /* try these lines here */
  /*LAL_CALL( LALDestroyCOMPLEX8FrequencySeries( &status, detector.transfer), &status );*/
 
  destroyCoherentGW( &waveform );

  fprintf( stdout, "222222\n" );
  fflush( stdout );

  if ( detector.site ) LALFree( detector.site);

  fprintf( stdout, "333333\n" );
  fflush( stdout );

  /* free the process params */
  while( procparams.processParamsTable )
  {
    this_proc_param = procparams.processParamsTable;
    procparams.processParamsTable = this_proc_param->next;
    free( this_proc_param );
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

  /*print a success message to stdout for parsing by exitcode */
  fprintf( stdout, "%s: EXITCODE0\n", argv[0] );
  fflush( stdout );

  exit( 0 ); 
}
