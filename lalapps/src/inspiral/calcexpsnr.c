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

/* template bank simulation params */
INT4  bankSim           = 0;            /* number of template bank sims */
CHAR *bankSimFileName   = NULL;         /* file contining sim_inspiral  */
FindChirpBankSimParams bankSimParams = { 0, 0, -1, -1, NULL, -1, NULL, NULL };
                                        /* template bank sim params     */

/* input filenames */
CHAR  *injectionFile    = NULL;         /* name of file containing injs */
CHAR  *specFileH1       = NULL;         /* name of file containing spec of H1 */
CHAR  *specFileH2       = NULL;         /* name of file containing spec of H2 */
CHAR  *specFileL1       = NULL;         /* name of file containing spec of L1 */

static INT4             numPoints       = 64;
static REAL4            fLow            = 40.0;

int main( int argc, char *argv[] )
{
  LALStatus                     status = blank_status;
  INT4                          flag = 0;
  INT4                          i;
  INT4                          j;
  INT4                          k;
  REAL8                         dynRange = 1.0;
  REAL4                         dynRangeExponent = 65.0;
  UINT4                         numInputPoints = numPoints;

  /* files contain PSD info */
  FILE                         *fpSpecH1;
  FILE                         *fpSpecH2;
  FILE                         *fpSpecL1;

  /* injection information */
  int                  numInjections = 0;
  SimInspiralTable    *injections = NULL;
  SimInspiralTable    *thisInj = NULL;


  /* template bank simulation variables */
  UINT4 bankSimCutLowIndex = 0;
  UINT4 bankSimCutHighIndex = 0;
  INT4  bankSimCount = 0;
  REAL4 matchNorm = 0;
  SnglInspiralTable  *loudestEventHead = NULL;
  SnglInspiralTable  *thisLoudestEvent = NULL;
  SimInspiralTable   *thisSimInspiral = NULL;
  SimInstParamsTable *thisSimInstParams = NULL;
  SimInspiralTable   *bankSimHead = NULL;
  SimInspiralTable   *thisBankSim = NULL;
  FindChirpBankSimParams *bankSimParamsPtr = NULL;

  /* raw input data storage */
  REAL4TimeSeries               chan;
  REAL8TimeSeries               strainChan;
  REAL4FrequencySeries          specH1;
  COMPLEX8FrequencySeries       resp;
  DataSegmentVector            *dataSegVec = NULL;
  COMPLEX8TimeSeries           *coherentInputData = NULL;

  /* output data */
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         searchsumm;
  MetadataTable         searchsummvars;
  MetadataTable         siminspiral;
  MetadataTable         siminstparams;
  MetadataTable         summvalue;
  MetadataTable         filtertable;
  SearchSummvarsTable  *this_search_summvar = NULL;
  SummValueTable       *this_summ_value = NULL;
  ProcessParamsTable   *this_proc_param = NULL;
  LIGOLwXMLStream       results;

  /* make sure all the output table pointers are null */
  siminspiral.simInspiralTable = NULL;
  siminstparams.simInstParamsTable = NULL;


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

  /* Try and read in power spectrum from files */
  /* H1 file */
  if ( !(fpSpecH1 = fopen( specFileH1, "r" )) )
  {
     fprintf( stdout, "unable to open the H1 spectrum file for reading\n" );
     fflush( stdout );
     exit (1);
  }

  /* H2 file */
  if ( !(fpSpecH2 = fopen( specFileH2, "r" )) )
  {
     fprintf( stdout, "unable to open the H2 spectrum file for reading\n" );
     fflush( stdout );
     close( fpSpecH1 );
     exit (1);
  }

  /* L1 file */
  if ( !(fpSpecL1 = fopen( specFileL1, "r" )) )
  {
     fprintf( stdout, "unable to open the L1 spectrum file for reading\n" );
     fflush( stdout );
     close( fpSpecH1 );
     close( fpSpecH2 );
     exit (1);
  }

  /* create storage for the power spectral estimate */
  memset( &specH1, 0, sizeof(REAL4FrequencySeries) );
  LAL_CALL( LALSCreateVector( &status, &(specH1.data), numPoints / 2 + 1 ),
      &status );

  /* read in spec and resp */
  for ( k = 0; k < numPoints/2 + 1; ++k )
  {
    if (( (flag = fscanf( fpSpecH1, "%f\n",
            &(specH1.data->data[k]) )) != 1 || flag == EOF )
        && k < numPoints/2 + 1 )
    {
       fprintf( stdout, "error reading input spectrum\n" );
       fflush( stdout );
       fclose( fpSpecH1 );
       fclose( fpSpecH2 );
       fclose( fpSpecL1 );
       exit(1);
    }
  }


  /* close all these files */
  close( fpSpecH1 );
  close( fpSpecH2 );
  close( fpSpecL1 );
  fprintf( stdout, "closed all spec files\n" );
  fflush( stdout );

  /* create the dynamic range exponent */
  dynRange = pow( 2.0, dynRangeExponent );

  /* create storage for the response function */
  memset( &resp, 0, sizeof(COMPLEX8FrequencySeries) );
  LAL_CALL( LALCCreateVector( &status, &(resp.data), numPoints / 2 + 1 ),
      &status );

 /* if we are using calibrated data set the response to unity */
 for( k = 0; k < resp.data->length; ++k )
 {
    resp.data->data[k].re = (REAL4) (1.0 / dynRange);
    resp.data->data[k].im = 0.0;
 }

  /* read in the simuation parameters from a sim_inspiral_table */
  bankSim = SimInspiralTableFromLIGOLw( &bankSimHead, injectionFile, 0, 0 );
  for ( thisBankSim = bankSimHead; thisBankSim; thisBankSim = thisBankSim->next )
  {
     /* set the time of the injection to zero so it is injected in  */
     /* the middle of the data segment: the other times are ignored */
    thisBankSim->geocent_end_time.gpsSeconds = 0;
    thisBankSim->geocent_end_time.gpsNanoSeconds = 0;
  }
  thisBankSim = bankSimHead;
  /* if ( vrbflg ) */
  fprintf( stdout, "Read %d bank sim parameters from %s\n", bankSim, injectionFile );

  bankSimCutLowIndex = XLALFindChirpBankSimInitialize( &specH1, &resp, fLow );
  /* if ( vrbflg ) */
  fprintf( stdout, "psd low frequency cutoff index = %d\n", bankSimCutLowIndex );


  /* stuff to do with chan */
  /* set the params of the input data time series */
  memset( &chan, 0, sizeof(REAL4TimeSeries) );
  LAL_CALL( LALSCreateVector( &status, &(chan.data), numInputPoints ), &status );

  bankSimCount = 0;
  do
  {
     fprintf( stdout, "bank simulation %d/%d\n", bankSimCount, bankSim );

     /* zero out the input data segment and the injection params */
     fprintf( stdout,
          "zeroing data stream and adding random injection for bank sim... " );
     memset( chan.data->data, 0, chan.data->length * sizeof(REAL4) );

     /* inject a random signal if we are doing a template bank simulation */
     if ( ! siminspiral.simInspiralTable )
        thisSimInspiral = siminspiral.simInspiralTable =
        XLALFindChirpBankSimInjectSignal( dataSegVec, &resp, thisBankSim,
        bankSimParamsPtr );
     else
        thisSimInspiral = thisSimInspiral->next =
        XLALFindChirpBankSimInjectSignal( dataSegVec, &resp, thisBankSim,
        bankSimParamsPtr );

  } while ( ++bankSimCount < bankSim ); /* end loop over bank simulations */



  /* print a success message to stdout for parsing by exitcode */
  fprintf( stdout, "%s: EXITCODE0\n", argv[0] );
  exit( 0 );
}

