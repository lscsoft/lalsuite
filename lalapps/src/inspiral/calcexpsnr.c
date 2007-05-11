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

RCSID( "$Id$" );

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "calcexpsnr"

#define USAGE \
  "lalapps_calcexpsnr [options]\n"\
"\nDefaults are shown in brackets\n\n" \
"  --help                   display this message\n"\
"  --verbose                be verbose\n"\
"  --version                version info\n"\
"  --spectrum-H1 FILE       FILE contains PSD info for H1\n"\
"  --spectrum-H2 FILE       FILE contains PSD info for H2\n"\
"  --spectrum-L1 FILE       FILE contains PSD info for L1\n"\
"  --inj-file    FILE       xml FILE contains injections\n"\
"\n"

/* debugging */
extern int vrbflg;                      /* verbocity of lal function    */


int main( int argc, char *argv[] )
{
  LALStatus                     status = blank_status;
  REAL8                         dynRange = 1.0;
  REAL4                         dynRangeExponent = 65.0;

  UINT4                          k;
  INT4                          numPoints       = 524288;
  REAL4                         fSampling       = 2048.;
  REAL8                         deltaT          = 1./fSampling;
  REAL8                         deltaF          = fSampling / numPoints;

  /* vars required to make freq series */
  LIGOTimeGPS                   epoch = { 0, 0 }; 
  REAL8                         f0 = 0.;
  LALUnit                       sampleUnits; 

  /* files contain PSD info */
  CHAR                         *injectionFile = NULL;         
  CHAR                         *specFileH1    = NULL;         
  CHAR                         *specFileH2    = NULL;         
  CHAR                         *specFileL1    = NULL;         

  /* from blindinj */
  COMPLEX8Vector               *unity;
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
  INT8                       waveformStartTime;

  /* injection information */
  int                  numInjections = 0;

  /* template bank simulation variables */
  INT4  injSimCount = 0;
  SimInspiralTable   *injectionHead = NULL;
  SimInspiralTable   *thisInjection = NULL;
  /*FindChirpBankSimParams *bankSimParamsPtr = NULL;*/

  /* raw input data storage */
  REAL8FrequencySeries          *specH1 = NULL;
  COMPLEX8FrequencySeries       *resp;
  REAL4TimeSeries               *chan;
  REAL4FFTPlan                  *pfwd;
  COMPLEX8FrequencySeries       *fftData;
  REAL4                          thisSnrsq = 0;

  /* needed for inj */
  CoherentGW                 waveform;
  PPNParamStruc              ppnParams;
  DetectorResponse           detector;
  InterferometerNumber       ifoNumber   = LAL_UNKNOWN_IFO;

  /* getopt arguments */
  struct option long_options[] =
  {
    /* these options set a flag */
    {"verbose",                 no_argument,       &vrbflg,           1 },
    /* these options do not set a flag */
    {"help",                    no_argument,       0,                'h'},
    {"version",                 no_argument,       0,                'V'},
    {"spectrum-H1",             required_argument, 0,                'a'},
    {"spectrum-H2",             required_argument, 0,                'b'},
    {"spectrum-L1",             required_argument, 0,                'c'},
    {"inj-file",                required_argument, 0,                'd'},
    {0, 0, 0, 0}
  };
  int c;

  set_debug_level("39");
  
  
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

    c = getopt_long_only( argc, argv, "a:b:c:d:hV", long_options, &option_index );

    /* detect the end of the options */
    if ( c == - 1 )
    {
      break;
    }

    switch ( c )
    {
      case 'h':
        fprintf( stderr, USAGE );
        exit( 0 );
        break;

      case 'a':
        /* create storage for the spectrum file name */
        optarg_len = strlen( optarg ) + 1;
        specFileH1 = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( specFileH1, optarg, optarg_len );
        /*ADD_PROCESS_PARAM( "string", "%s", optarg );*/
        break;

      case 'b':
        /* create storage for the spectrum file name */
        optarg_len = strlen( optarg ) + 1;
        specFileH2 = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( specFileH2, optarg, optarg_len );
        /*ADD_PROCESS_PARAM( "string", "%s", optarg );*/
        break;

      case 'c':
        /* create storage for the spectrum file name */
        optarg_len = strlen( optarg ) + 1;
        specFileL1 = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( specFileL1, optarg, optarg_len );
        /*ADD_PROCESS_PARAM( "string", "%s", optarg );*/
        break;

      case 'd':
        /* create storage for the injection file name */
        optarg_len = strlen( optarg ) + 1;
        injectionFile = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( injectionFile, optarg, optarg_len );
        /*ADD_PROCESS_PARAM( "string", "%s", optarg );*/
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "calculation of expected SNR of injections\n"
            "Gareth Jones\n"
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n" );
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

  fprintf( stdout, "injection file is %s\n", injectionFile );
  fprintf( stdout, "H1 spec file is   %s\n", specFileH1 );
  fprintf( stdout, "H2 spec file is   %s\n", specFileH2 );
  fprintf( stdout, "L1 spec file is   %s\n", specFileL1 );

  /* create vector for H1 spectrum */
  memcpy( &sampleUnits, &lalADCCountUnit, sizeof(LALUnit) );
  LAL_CALL( LALCreateREAL8FrequencySeries ( &status, &specH1, "",epoch, f0, deltaF, sampleUnits, (numPoints / 2 + 1) ), &status);

  /* read in H1 spectrum */ 
  LAL_CALL( LALDReadFrequencySeries(&status, specH1, specFileH1), &status );
  fprintf( stdout, "read in H1 spec file\n" );
  fflush( stdout );

  /*
  *fprintf( stdout, "PSD value is %e \n", specH1->data->data[1] );
  *fflush( stdout );
  */
  
  /* write out spectrum test */
  LALDPrintFrequencySeries(specH1, "Test.dat" );
  fprintf( stdout, "written H1 spec file\n" );
  fflush( stdout );

  /* create the dynamic range exponent */
/*  dynRange = pow( 2.0, dynRangeExponent );*/

  /* create storage for the response function */
/*  memset( &resp, 0, sizeof(COMPLEX8FrequencySeries) );
  *LAL_CALL( LALCCreateVector( &status, &(resp.data), numPoints / 2 + 1 ), &status );
*/
 /* if we are using calibrated data set the response to unity */
/* for( k = 0; k < resp.data->length; ++k )
 *{
 *   resp.data->data[k].re = (REAL4) (1.0 / dynRange);
 *   resp.data->data[k].im = 0.0;
 *}
*/

  /* set up the channel to which we add the injection */
/*  XLALReturnIFO( ifo, ifoNumber );
 * LALSnprintf( name, LALNameLength * sizeof(CHAR), "%s:INJECT", ifo );
 */ 

 LAL_CALL( LALCreateREAL4TimeSeries( &status, &chan, "", epoch, f0, deltaT, sampleUnits, numPoints ), &status );
  if ( ! chan )
  {
    exit( 1 );
  }

  /*
   *
   * set up the response function
   *
   */
 LAL_CALL( LALCreateCOMPLEX8FrequencySeries( &status, &resp, chan->name, chan->epoch, f0, deltaF, strainPerCount, (numPoints / 2 + 1) ), &status );

  /* read in the simuation parameters from a sim_inspiral_table */
  /* bankSim = SimInspiralTableFromLIGOLw( &bankSimHead, injectionFile, 0, 0 ); */

  /* set endtime to 0 so that we read in all events */
  LAL_CALL(numInjections = SimInspiralTableFromLIGOLw( &injectionHead, injectionFile, 0, 0), &status);

  for ( thisInjection = injectionHead; thisInjection; thisInjection = thisInjection->next )
  /*for ( thisBankSim = bankSimHead; thisBankSim; thisBankSim = thisBankSim->next )*/
  {
     /* set the time of the injection to zero so it is injected in  */
     /* the middle of the data segment: the other times are ignored */
    thisInjection->geocent_end_time.gpsSeconds = 0;
    thisInjection->geocent_end_time.gpsNanoSeconds = 0;
  }
  thisInjection = injectionHead;

  /* if ( vrbflg ) */
  fprintf( stdout, "Read %d inj parameters from %s\n", numInjections, injectionFile );


  
  /* setting fixed waveform injection parameters */
  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  ppnParams.deltaT   = deltaT;
  ppnParams.lengthIn = 0;
  ppnParams.ppn      = NULL;


  injSimCount = 0;
  do
  {
     fprintf( stdout, "injection %d/%d\n", injSimCount, numInjections );

     /* reset waveform structure */
     memset( &waveform, 0, sizeof(CoherentGW) );

     /* reset chan structure */
     memset( chan->data->data, 0, chan->data->length * sizeof(REAL4) );

     /* create the waveform, amp, freq phase etc */
     LAL_CALL( LALGenerateInspiral(&status, &waveform, thisInjection, &ppnParams), &status);

    /* loop over ifo */
    for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
    {
        /* allocate memory and copy the parameters describing the freq series */
        memset( &detector, 0, sizeof( DetectorResponse ) );
        detector.site = (LALDetector *) LALMalloc( sizeof(LALDetector) );
        XLALReturnDetector( detector.site, ifoNumber );
 
        fprintf( stdout, "ifoNumber %d\n", ifoNumber );
        fflush( stdout );


        /************** wholesale copy from blindinj  ******************************/
        detector.transfer = XLALCreateCOMPLEX8FrequencySeries( chan->name,
           &(chan->epoch), f0, deltaF, &strainPerCount, ( numPoints / 2 + 1 ) );

        XLALUnitInvert( &(detector.transfer->sampleUnits), &(resp->sampleUnits) );

        /* invert the response function to get the transfer function */
        unity = XLALCreateCOMPLEX8Vector( resp->data->length );
        for ( k = 0; k < unity->length; ++k )
        {
           unity->data[k].re = 1.0;
           unity->data[k].im = 0.0;
        }

        XLALCCVectorDivide( detector.transfer->data, unity, resp->data );
        XLALDestroyCOMPLEX8Vector( unity );
         
        /* from bank sim routines */
        thisInjection->geocent_end_time.gpsSeconds = 0;
        thisInjection->geocent_end_time.gpsNanoSeconds = 0;

        waveformStartTime = XLALGPStoINT8( &(thisInjection->geocent_end_time));
        waveformStartTime -= (INT8) ( 1000000000.0 * ppnParams.tc );
        waveformStartTime += (INT8) ( 1000000000.0 *
          ((REAL8) (chan->data->length - ppnParams.length) / 2.0) * chan->deltaT
          );
    
        XLALINT8toGPS( &(waveform.a->epoch), waveformStartTime );
        memcpy(&(waveform.f->epoch), &(waveform.a->epoch), sizeof(LIGOTimeGPS) );
        memcpy(&(waveform.phi->epoch), &(waveform.a->epoch), sizeof(LIGOTimeGPS) );
        /***************** ends here *****************************************/
  
       /* perform the injection */
       LAL_CALL( LALSimulateCoherentGW(&status, chan, &waveform, &detector ), &status); 


       /* more theft */

       /* fft the output */
       pfwd = XLALCreateForwardREAL4FFTPlan( chan->data->length, 0 );
       fftData = XLALCreateCOMPLEX8FrequencySeries( chan->name,
           &(chan->epoch), 0, 1.0/chan->deltaT, &lalDimensionlessUnit,
           chan->data->length/2 + 1 );
       XLALREAL4TimeFreqFFT( fftData, chan, pfwd );
       XLALDestroyREAL4FFTPlan( pfwd );

       /* compute the SNR for initial LIGO at design */
       thisSnrsq = 0;
       for ( k = 0; k < fftData->data->length; k++ )
       {
         REAL8 freq;
         /* use correct psd !!!!! */
         REAL8 sim_psd_value;
         freq = fftData->deltaF * k;
         LALLIGOIPsd( NULL, &sim_psd_value, freq );

         fprintf( stdout, "k= %d  freq = %e sim_psd_value = %e  \n", k, freq, sim_psd_value );
         fprintf( stdout, fftData->data->data[k].re           

         
         thisSnrsq += fftData->data->data[k].re * fftData->data->data[k].re /
           sim_psd_value;
         thisSnrsq += fftData->data->data[k].re * fftData->data->data[k].re /
           sim_psd_value;
       }
       fprintf( stdout, "thisSnrsq %f\n", thisSnrsq );
       fprintf( stdout, "thisSnrsq %e\n", thisSnrsq );
       thisSnrsq *= 4*fftData->deltaF;
       XLALDestroyCOMPLEX8FrequencySeries( fftData );

       fprintf( stdout, "thisSnrsq %f\n", thisSnrsq );
       fprintf( stdout, "thisSnrsq %e\n", thisSnrsq );
       fflush( stdout );
       /* end of theft */ 





     }
     /* end loop over ifo */




    /* increment the bank sim sim_inspiral table if necessary */
    if ( injectionHead )
    {
      thisInjection = thisInjection->next;
    }

  } while ( ++injSimCount < numInjections ); /* end loop over injections */

  fprintf( stdout, "out of loop\n" );
  fflush( stdout );
  

  /* I need to remember to destroy all the vectors i have created */
  LAL_CALL( LALDestroyREAL4TimeSeries( &status, chan), &status );
  LAL_CALL( LALDestroyCOMPLEX8FrequencySeries( &status, resp), &status );
  LAL_CALL( LALDestroyREAL8FrequencySeries ( &status, specH1 ), &status);

  /* print a success message to stdout for parsing by exitcode */
  fprintf( stdout, "%s: EXITCODE0\n", argv[0] );
  exit( 0 );
}

