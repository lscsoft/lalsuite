/*
*  Copyright (C) 2007 Patrick Brady, Stephen Fairhurst
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
 * File Name: blindinj.c
 *
 * Author: Fairhurst, S
 * 
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <config.h>

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <lalapps.h>
#include <processtable.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataInspiralUtils.h>
#include <lal/LIGOMetadataRingdownUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/Inject.h>
#include <lal/InspiralInjectionParams.h>
#include <lal/GenerateInspiral.h>
#include <lal/GenerateInspRing.h>
#include <lal/FindChirp.h>
#include <lal/GenerateRing.h>
#include <lal/RingUtils.h>
#include <lal/LALNoiseModels.h>
#include <lal/RealFFT.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/VectorOps.h>
#include <LALAppsVCSInfo.h>

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "blindinj"

#define USAGE \
  "lalapps_blindinj [options]\n"\
"\nDefaults are shown in brackets\n\n" \
"  --help                   display this message\n"\
"  --version                print version information and exit\n"\
"  --verbose                be verbose\n"\
"  --gps-start-time TIME    start time of injection\n"\
"  --injection-type TYPE    type of injection, must be one of \n"\
"                           (strain, etmx, etmy)\n"\
"  --seed           SEED    seed random number generator with SEED (1)\n"\
"  --debug-level    LEVEL   set the LAL debug level to LEVEL\n"\
"\n"

/* global definitions */

typedef enum 
{ 
  noResponse, 
  unityResponse, 
  LIGOdesign, 
  actuationX, 
  actuationY 
} 
ResponseFunction;

typedef struct actuationparameters {
  REAL4           ETMXcal;
  REAL4           ETMYcal;
  REAL4           pendFX;
  REAL4           pendFY;
  REAL4           pendQX;
  REAL4           pendQY;
  REAL4           length;
} ActuationParameters;

static ProcessParamsTable *next_process_param( 
    const char *name, 
    const char *type,
    const char *fmt, ... );

static void destroyCoherentGW( CoherentGW *waveform );

static REAL4TimeSeries *injectWaveform( 
    LALStatus            *status,
    SimInspiralTable     *inspInj,
    SimRingdownTable     *ringdownevents,
    ResponseFunction      responseType,
    InterferometerNumber  ifoNumber,
    LIGOTimeGPS           epoch);

extern int vrbflg;

/* Actuation function parameters */
ActuationParameters actuationParams[LAL_NUM_IFO];
REAL8         dynRange = 1.0/3.0e-23;/* dynamic range rescaling       */
INT4          duration = 100.0;      /* length of output data stream  */
INT4          sampleRate = 16384;    /* sample rate of output channel */

/* default values of various things */
REAL4         fLower = 30;         /* start frequency of injection */
REAL4         mergerLength  = 10;  /* length in ms of the merger */
REAL8         longestSignal = 95.0;/* length of 1.0 - 1.0 from 30 Hz */
REAL8         timeWindow = 4;      /* injection can be delayed by up to this
                                      number of seconds */
UINT4         numInjections = 0;   /* number of injections we have generated */
INT4          randSeed = 1;
REAL4         minNSMass = 1.0;     /* minimum NS component mass */
REAL4         maxNSMass = 2.0;     /* maximum NS component mass */
REAL4         minBHMass = 2.0;     /* minimum BH component mass */
REAL4         maxBHMass = 30.0;    /* maximum BH component mass */
REAL4         minTotalMass = 0.0;  /* minimum total mass */
REAL4         maxTotalMass = 35.0; /* maximum total mass */

REAL4         minNSSpin = 0.0;     /* minimum NS component spin */
REAL4         maxNSSpin = 0.2;     /* maximum NS component spin */
REAL4         minBHSpin = 0.0;     /* minimum BH component spin */
REAL4         maxBHSpin = 1.0;     /* maximum BH component spin */

REAL4         BNSfrac = 0.35;      /* fraction of injections which are BNS */
REAL4         BBHfrac = 0.35;      /* fraction of injections which are BBH */
/* 1 - BNSfrac - BBHfrac are NS-BH inj  */

REAL4         bnsSnrMean = 9.0;  /* mean single ifo snr of bns injection */
REAL4         bnsSnrStd = 0.5;   /* std of single ifo snr of bns injection */
REAL4         snrMean = 12.0;    /* mean single ifo snr of injection */
REAL4         snrStd  = 1.0;     /* std of single ifo snr of injection */
/* snrs assume detectors at design     */
REAL4Vector  *normalDev;         /* vector to store normally distributed vars*/

/* functions */

ProcessParamsTable *next_process_param( 
    const char *name, 
    const char *type,
    const char *fmt, ... )
{
  ProcessParamsTable *pp;
  va_list ap;
  pp = calloc( 1, sizeof( *pp ) );

  if ( ! pp )
  {
    perror( "next_process_param" );
    exit( 1 );
  }
  strncpy( pp->program, PROGRAM_NAME, LIGOMETA_PROGRAM_MAX );
  if ( ! strcmp( name, "userTag" ) || ! strcmp( name, "user-tag" ) )
    snprintf( pp->param, LIGOMETA_PARAM_MAX, "-userTag" );
  else
    snprintf( pp->param, LIGOMETA_PARAM_MAX, "--%s", name );
  strncpy( pp->type, type, LIGOMETA_TYPE_MAX );
  va_start( ap, fmt );
  vsnprintf( pp->value, LIGOMETA_VALUE_MAX, fmt, ap );
  va_end( ap );

  return pp;
}

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

static REAL4TimeSeries *injectWaveform( 
    LALStatus            *status,
    SimInspiralTable     *inspInj,
    SimRingdownTable     *ringdownevents,
    ResponseFunction      responseType,
    InterferometerNumber  ifoNumber,
    LIGOTimeGPS           epoch)
{
  REAL4TimeSeries           *chan;
  COMPLEX8FrequencySeries   *resp;
  COMPLEX8Vector            *unity;
  CHAR                       name[LALNameLength];
  CHAR                       ifo[LIGOMETA_IFO_MAX];
  PPNParamStruc              ppnParams;
  SimRingdownTable          *thisRingdownEvent = NULL;

  INT8                       waveformStartTime;
  DetectorResponse           detector;
  CoherentGW                 waveform, *wfm;
  ActuationParameters        actData = actuationParams[ifoNumber];
  UINT4 i,k;
  int injectSignalType = LALRINGDOWN_IMR_INJECT;
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
  FILE  *fp = NULL;
  char  fileName[FILENAME_MAX];

  /* set up the channel to which we add the injection */
  XLALReturnIFO( ifo, ifoNumber );
  snprintf( name, LALNameLength, "%s:INJECT", ifo );
  chan = XLALCreateREAL4TimeSeries( name, &epoch, 0, 1./sampleRate, 
      &lalADCCountUnit, sampleRate * duration );
  if ( ! chan )
  {
    exit( 1 );
  }

  memset( chan->data->data, 0, chan->data->length * sizeof(REAL4) );
 
  thisRingdownEvent = ringdownevents;
  
  /*
   *
   * Generate the Waveforms
   *
   */

  /* fixed waveform injection parameters */
  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  ppnParams.deltaT   = chan->deltaT;
  ppnParams.lengthIn = 0;
  ppnParams.ppn      = NULL;

  memset( &waveform, 0, sizeof(CoherentGW) );

  LAL_CALL( LALGenerateInspiral(status, &waveform, inspInj, &ppnParams), 
      status);

  /* add the ringdown */
  wfm = XLALGenerateInspRing( &waveform, inspInj, thisRingdownEvent, 
      injectSignalType );

  if ( !wfm )
  {
    fprintf( stderr, "Failed to generate the waveform \n" );
    if (xlalErrno == XLAL_EFAILED)
    {
      fprintf( stderr, "Too much merger\n");
      XLALDestroyREAL4TimeSeries( chan );     
      xlalErrno = XLAL_SUCCESS;
      return ( NULL );
    }
    else exit ( 1 );
  }

  waveform = *wfm;

  /* write out the waveform */
  if ( ifoNumber == LAL_IFO_H1 && vrbflg )
  {
    fprintf(  stdout, 
        "Writing out A+, Ax, f, phi, shift for waveform "
        "to file named INSPIRAL_WAVEFORM.dat\n");
    snprintf( fileName, FILENAME_MAX, "INSPIRAL_WAVEFORM.dat");
    fp = fopen(fileName, "w");
    for ( i = 0; i < waveform.phi->data->length; i++ )
    {
      if ( waveform.shift ) fprintf( fp, "%e\t %e\t %f\t %f\t %f\n", 
          waveform.a->data->data[2*i],
          waveform.a->data->data[2*i+1], waveform.f->data->data[i], 
          waveform.phi->data->data[i] , waveform.shift->data->data[i] );

      else fprintf( fp, "%e\t %e\t %f\t %f\n", 
          waveform.a->data->data[2*i],
          waveform.a->data->data[2*i+1], waveform.f->data->data[i], 
          waveform.phi->data->data[i] );
    }
    fclose( fp );
  }

  /* 
   *
   * set up the response function
   *
   */

  resp = XLALCreateCOMPLEX8FrequencySeries( chan->name, &(chan->epoch), 0,
      1.0 / ( duration ), &strainPerCount, ( sampleRate * duration / 2 + 1 ) );

  switch ( responseType )
  {
    case noResponse:
      fprintf( stderr, "Must specify the response function\n" );
      exit( 1 );
      break;

    case unityResponse:
      /* set the response function to unity */
      for ( k = 0; k < resp->data->length; ++k )
      {
        resp->data->data[k].re = (REAL4) (1.0 /dynRange);
        resp->data->data[k].im = 0.0;
      }
      break;

    case LIGOdesign:
      /* set the response function to LIGO design */
      for ( k = 0; k < resp->data->length; ++k )
      {
        REAL8 sim_psd_freq = (REAL8) k * resp->deltaF;
        REAL8 sim_psd_value;
        LALLIGOIPsd( NULL, &sim_psd_value, sim_psd_freq );
        resp->data->data[k].re = (REAL4) pow( sim_psd_value, 0.5 ) / dynRange;
        resp->data->data[k].im = 0.0;
      }
      break;

    case actuationX:
      /* actuation units are m/count so we must divide by arm length */
      resp = generateActuation( resp, actData.ETMXcal / actData.length,
          actData.pendFX, actData.pendQX);
      break;

    case actuationY:
      /* actuation units are m/count so we must divide by arm length */
      resp = generateActuation( resp, actData.ETMYcal / actData.length,
          actData.pendFY, actData.pendQY);
      break;
  }

  /*
   *
   * set up the detector structure
   *
   */

  /* allocate memory and copy the parameters describing the freq series */
  memset( &detector, 0, sizeof( DetectorResponse ) );
  detector.site = (LALDetector *) LALMalloc( sizeof(LALDetector) );
  XLALReturnDetector( detector.site, ifoNumber );
  detector.transfer = XLALCreateCOMPLEX8FrequencySeries( chan->name, 
      &(chan->epoch), 0, 1.0 / ( duration ), &strainPerCount, 
      ( sampleRate * duration / 2 + 1 ) );

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

  /* set the injection time */
  waveformStartTime = XLALGPSToINT8NS( &(inspInj->geocent_end_time));
  waveformStartTime -= (INT8) ( 1000000000.0 * ppnParams.tc );

  XLALINT8NSToGPS( &(waveform.a->epoch), waveformStartTime );
  memcpy(&(waveform.f->epoch), &(waveform.a->epoch), 
      sizeof(LIGOTimeGPS) );
  memcpy(&(waveform.phi->epoch), &(waveform.a->epoch),
      sizeof(LIGOTimeGPS) );

  /* perform the injection */
  LAL_CALL( LALSimulateCoherentGW( status, chan, &waveform, &detector ),
      status);

  destroyCoherentGW( &waveform );
  XLALDestroyCOMPLEX8FrequencySeries( detector.transfer );
  if ( detector.site ) LALFree( detector.site );

  XLALDestroyCOMPLEX8FrequencySeries( resp ); 
  return( chan );
}

/*
 *
 * MAIN CODE
 *
 */

int main( int argc, char *argv[] )
{
  LALStatus             status = blank_status;
  LIGOTimeGPS           gpsStartTime = {0, 0};
  LIGOTimeGPS           earliestEndTime = {0, 0};
  ResponseFunction      injectionResponse = noResponse;

  /* program variables */
  RandomParams         *randParams  = NULL;
  REAL4                 massPar     = 0;
  MassDistribution      mDist       = logComponentMass;
  InterferometerNumber  ifoNumber   = LAL_UNKNOWN_IFO;
  REAL4                 desiredSnr  = 0;           
  REAL4                 snrsqAt1Mpc = 0;           
  REAL4TimeSeries      *chan        = NULL;
  REAL4FFTPlan         *pfwd;
  COMPLEX8FrequencySeries *fftData;

  /* waveform */
  CHAR waveform[LIGOMETA_WAVEFORM_MAX];

  /* xml output data */
  CHAR                  fname[FILENAME_MAX];
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         inspInjections;
  MetadataTable         ringInjections;
  ProcessParamsTable   *this_proc_param;
  SimInspiralTable     *inj  = NULL;
  SimRingdownTable     *ringList = NULL;
  LIGOLwXMLStream       xmlfp;
  FILE                 *fp = NULL;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"help",                    no_argument,       0,                'h'},
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"version",                 no_argument,       0,                'V'},
    {"gps-start-time",          required_argument, 0,                'a'},
    {"injection-type",          required_argument, 0,                't'},
    {"seed",                    required_argument, 0,                's'},
    {"debug-level",             required_argument, 0,                'z'},
    {0, 0, 0, 0}
  };
  int c;


  /*taken from Calibration CVS file: 
   * calibration/frequencydomain/runs/S5/H1/model/V3/H1DARMparams_849677446.m */
  actuationParams[LAL_IFO_H1].ETMXcal = -0.795e-9;
  actuationParams[LAL_IFO_H1].pendFX  = 0.767;
  actuationParams[LAL_IFO_H1].pendQX  = 10.0;
  actuationParams[LAL_IFO_H1].ETMYcal = -0.827e-9;
  actuationParams[LAL_IFO_H1].pendFY  = 0.761;
  actuationParams[LAL_IFO_H1].pendQY  = 10.0;
  actuationParams[LAL_IFO_H1].length  = 4000.0;

  /*taken from Calibration CVS file: 
   * calibration/frequencydomain/runs/S5/H2/model/V3/H2DARMparams_849678155.m */
  actuationParams[LAL_IFO_H2].ETMXcal = -0.876e-9;
  actuationParams[LAL_IFO_H2].pendFX  = 0.749;
  actuationParams[LAL_IFO_H2].pendQX  = 10.0;
  actuationParams[LAL_IFO_H2].ETMYcal = -0.912e-9;
  actuationParams[LAL_IFO_H2].pendFY  = 0.764;
  actuationParams[LAL_IFO_H2].pendQY  = 10.0;
  actuationParams[LAL_IFO_H2].length  = 2000.0;

  /*taken from Calibration CVS file: 
   * calibration/frequencydomain/runs/S5/L1/model/V3/L1DARMparams_841930071.m */
  actuationParams[LAL_IFO_L1].ETMXcal = -0.447e-9;
  actuationParams[LAL_IFO_L1].pendFX  = 0.766;
  actuationParams[LAL_IFO_L1].pendQX  = 100.0;
  actuationParams[LAL_IFO_L1].ETMYcal = -0.438e-9;
  actuationParams[LAL_IFO_L1].pendFY  = 0.756;
  actuationParams[LAL_IFO_L1].pendQY  = 100.0;
  actuationParams[LAL_IFO_L1].length  = 4000.0;

  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "33" );


  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) 
    calloc( 1, sizeof(ProcessTable) );
  XLALGPSTimeNow(&(proctable.processTable->start_time));
  XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME, LALAPPS_VCS_IDENT_ID,
      LALAPPS_VCS_IDENT_STATUS, LALAPPS_VCS_IDENT_DATE, 0);
  snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );

  /* clear the waveform field */
  memset( waveform, 0, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR) );
  snprintf( waveform, LIGOMETA_WAVEFORM_MAX, "SpinTaylorthreePN");

  /*
   *
   * parse command line arguments
   *
   */


  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    long int gpsinput;

    c = getopt_long_only( argc, argv, "a:hs:t:V", long_options, &option_index );

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
        gpsinput = atol( optarg );
        if ( gpsinput < 441417609 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is prior to " 
              "Jan 01, 1994  00:00:00 UTC:\n"
              "(%ld specified)\n",
              long_options[option_index].name, gpsinput );
          exit( 1 );
        }
        gpsStartTime.gpsSeconds = gpsinput;

        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "int", 
              "%ld", gpsinput );
        break;

      case 's':
        randSeed = atoi( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "int", 
              "%d", randSeed );
        break;

      case 't':
        /* set the injection type */
        if ( ! strcmp( "strain", optarg ) )
        {
          injectionResponse = unityResponse;
        }
        else if ( ! strcmp( "etmx", optarg ) )
        {
          injectionResponse = actuationX;
        }
        else if ( ! strcmp( "etmy", optarg ) )
        {
          injectionResponse = actuationY;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown injection type specified: "
              "%s (must be strain, etmx or etmy)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        next_process_param( long_options[option_index].name, "string", "%s", 
            optarg );
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "blind hardware injection generation routine\n" 
            "Stephen Fairhurst\n");
        XLALOutputVersionString(stderr, 0);
        exit( 0 );
        break;

      case 'z':
        set_debug_level( optarg );
        next_process_param( long_options[option_index].name, "int", "%d",
            optarg );
        break;

      case 'h':
      case '?':
        fprintf( stderr, USAGE );
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
  if ( injectionResponse == noResponse )
  {
    fprintf( stderr, "Must specify the --injection-type\n" );
    exit( 1 );
  }

  if ( gpsStartTime.gpsSeconds <= 0 )
  {
    fprintf( stderr, "Must specify the --gps-start-time\n" );
    exit( 1 );
  }

  /*
   *
   * initialization
   *
   */


  /* initialize the random number generator */
  randParams = XLALCreateRandomParams( randSeed );

  /* null out the head of the linked list */
  memset( &inspInjections, 0, sizeof(MetadataTable) );
  memset( &ringInjections, 0, sizeof(MetadataTable) );
  

  /*
   *
   * create a random injection
   *
   */

  /* create the sim_inspiral table */
  inspInjections.simInspiralTable = inj = (SimInspiralTable *)
    LALCalloc( 1, sizeof(SimInspiralTable) );
  
  ringInjections.simRingdownTable = ringList = (SimRingdownTable *) 
    LALCalloc( 1, sizeof(SimRingdownTable) );
  

  /* set the geocentric end time of the injection */
  earliestEndTime = gpsStartTime;
  earliestEndTime = *XLALGPSAdd( &earliestEndTime, longestSignal );
  inj = XLALRandomInspiralTime( inj, randParams, earliestEndTime, 
      timeWindow );

  /* set the distance of the injection to 1 Mpc */
  inj->distance = 1.0;

  /* set the masses */
  massPar = XLALUniformDeviate( randParams );

  if ( vrbflg) fprintf(  stdout, 
      "Random variable to determine inj type = %f\n", 
      massPar);
  
  while ( numInjections == 0 )
  {

    normalDev = XLALCreateREAL4Vector( 1 );
    XLALNormalDeviates(normalDev, randParams);

    if ( massPar < BNSfrac )
    {
      if ( vrbflg ) fprintf( stdout, "Generating a BNS injection\n" );
      inj = XLALRandomInspiralMasses( inj, randParams, mDist,
          minNSMass, maxNSMass, minNSMass, maxNSMass, minTotalMass, maxTotalMass );
      inj = XLALRandomInspiralSpins( inj, randParams, minNSSpin,
          maxNSSpin, minNSSpin, maxNSSpin, -1.0, 1.0, 0.0, 0.1, 0, uniformSpinDist, 0.0, 0.0, 0.0, 0.0);
      desiredSnr = bnsSnrMean + bnsSnrStd * normalDev->data[0]; 
    }
    else if ( massPar < (BNSfrac + BBHfrac) )
    {
      if ( vrbflg ) fprintf( stdout, "Generating a BBH injection\n" );
      inj = XLALRandomInspiralMasses( inj, randParams, mDist,
          minBHMass, maxBHMass, minBHMass, maxBHMass, minTotalMass, maxTotalMass );
      inj = XLALRandomInspiralSpins( inj, randParams, minBHSpin,
          maxBHSpin, minBHSpin, maxBHSpin , -1.0, 1.0, 0.0, 0.1, 0, uniformSpinDist, 0.0, 0.0, 0.0, 0.0);
      desiredSnr = snrMean + snrStd * normalDev->data[0]; 
    }
    else
    {
      if ( vrbflg ) fprintf( stdout, "Generating an NS - BH injection\n" );
      inj = XLALRandomInspiralMasses( inj, randParams, mDist,
          minNSMass, maxNSMass, minBHMass, maxBHMass, minTotalMass, maxTotalMass );
      inj = XLALRandomInspiralSpins( inj, randParams, minNSSpin,
          maxNSSpin, minBHSpin, maxBHSpin , -1.0, 1.0, 0.0, 0.1, 0, uniformSpinDist, 0.0, 0.0, 0.0, 0.0);
      desiredSnr = snrMean + snrStd * normalDev->data[0]; 
    }
    XLALDestroyVector( normalDev );


    /* set the sky location */
    InclDistribution iDist = uniformInclDist ;
    inj = XLALRandomInspiralSkyLocation( inj, randParams );
    inj = XLALRandomInspiralOrientation( inj, randParams, iDist, 0);

    /* set the source and waveform fields */
    snprintf( inj->source, LIGOMETA_SOURCE_MAX, "???" );
    memcpy( inj->waveform, waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR));

    /* populate the site specific information */
    inj = XLALPopulateSimInspiralSiteInfo( inj );

    /* finally populate the flower */
    inj->f_lower = fLower;


    /*
     *
     * perform the injection
     *
     */

    if ( vrbflg ) fprintf( stdout, "Injection details\n"
        "mass 1 = %.2f;\t spin 1 = (%.2f, %.2f, %.2f)\n"
        "mass 2 = %.2f;\t spin 2 = (%.2f, %.2f, %.2f)\n"
        "Hanford effective distance = %.2f Mpc\n"
        "Livingston effective distance = %.2f Mpc\n",
        inj->mass1, inj->spin1x, inj->spin1y, inj->spin1z,
        inj->mass2, inj->spin2x, inj->spin2y, inj->spin2z,
        inj->eff_dist_h, 
        inj->eff_dist_l );

    for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
    {
      if ( ifoNumber == LAL_IFO_H1 || ifoNumber == LAL_IFO_L1 )
      {
        ResponseFunction  responseType = unityResponse;
        REAL4             thisSnrsq = 0;
        UINT4             k;
        
        chan = injectWaveform( &status, inj, ringList, responseType, ifoNumber, 
            gpsStartTime);

        if ( ! chan )
        {
          if ( vrbflg ) fprintf( stdout, 
              "Unable to generate a realistic merger for this inspiral\n" 
              "Trying again with different parameters\n" );
          break;
        }
        else if ( ifoNumber == LAL_IFO_H1 ) numInjections++;


        /*
         *
         * calculate the snr for the injection
         *
         */

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
          REAL8 sim_psd_value;
          freq = fftData->deltaF * k;
          LALLIGOIPsd( NULL, &sim_psd_value, freq );
          thisSnrsq += fftData->data->data[k].re * fftData->data->data[k].re /
            sim_psd_value;
          thisSnrsq += fftData->data->data[k].im * fftData->data->data[k].im /
            sim_psd_value;
        }
        thisSnrsq *= 4*fftData->deltaF;
        XLALDestroyCOMPLEX8FrequencySeries( fftData );

        if ( ifoNumber == LAL_IFO_H1 )
        {
          /* add in H2, assuming half the snr */
          snrsqAt1Mpc += 5.0/4.0 * thisSnrsq;
        }
        else
        {
          snrsqAt1Mpc += thisSnrsq;
        }
        if ( vrbflg ) 
        {
          CHAR  ifo[LIGOMETA_IFO_MAX];
          XLALReturnIFO( ifo, ifoNumber );
          fprintf(stdout, 
              "For %s, the SNR at distance of 1 Mpc is %.2f\n", ifo, 
              pow(thisSnrsq, 0.5));
        }
        XLALDestroyREAL4TimeSeries( chan );
      }
    }

    if ( numInjections )
    {
      /* scale the distance so the combined snr is equal to desired value */
      desiredSnr *= 1.5;
      inj->distance = 1.0 * pow( snrsqAt1Mpc, 0.5 ) / desiredSnr;
      inj = XLALPopulateSimInspiralSiteInfo( inj );
      if ( vrbflg ) fprintf( stdout, 
          "Rescaling the distance to %.2f to give a combined snr of %.2f\n", 
          inj->distance, desiredSnr);

      /*
       *
       * compute the waveforms for injection
       *
       */

      for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
      {
        if ( ifoNumber == LAL_IFO_H1 || ifoNumber == LAL_IFO_H2 || 
            ifoNumber == LAL_IFO_L1 )
        {
          UINT4             k;
          CHAR              ifo[LIGOMETA_IFO_MAX];
          CHAR              type[LIGOMETA_COMMENT_MAX];

          XLALReturnIFO( ifo, ifoNumber );

          if ( injectionResponse == unityResponse )
          {
            snprintf( type, LIGOMETA_COMMENT_MAX, "STRAIN");
            if ( vrbflg ) fprintf( stdout, 
                "Generating h(t) injection for %s\n", ifo );
          }
          if ( injectionResponse == actuationX ) 
          {
            snprintf( type, LIGOMETA_COMMENT_MAX, "ETMX");
            if ( vrbflg ) fprintf( stdout, 
                "Generating ETMX hardware injection for %s\n", ifo );
          }
          if ( injectionResponse == actuationY ) 
          {
            snprintf( type, LIGOMETA_COMMENT_MAX, "ETMY");
            if (vrbflg ) fprintf( stdout, 
                "Generating ETMY hardware injection for %s\n", ifo );
          }
          
          chan = injectWaveform( &status, inj, ringList, injectionResponse, 
              ifoNumber, gpsStartTime);

          snprintf( fname, FILENAME_MAX, 
              "%s-HARDWARE-INJECTION_%d_%s-%d-%d.txt",
              ifo, randSeed, type, gpsStartTime.gpsSeconds, duration );
          if ( vrbflg ) fprintf( stdout, "Writing waveform to %s\n", fname);

          fp = fopen( fname, "w" ); 
          for ( k = 0; k < chan->data->length; k++ )
          {
            if ( injectionResponse == unityResponse )
            {
              /* we have to fix the dynamic range scaling */
              fprintf( fp, "%e\n", chan->data->data[k] / dynRange );
            }
            else
            {
              fprintf( fp, "%e\n", chan->data->data[k] );
            }
          }

          fclose( fp );

          XLALDestroyREAL4TimeSeries( chan );
        }
      }
    }
  }

  /*
   *
   * write output to LIGO_LW XML file
   *
   */

  /* create the output file name */
  snprintf( fname, sizeof(fname), "HL-INJECTIONS_%d-%d-%d.xml", 
      randSeed, gpsStartTime.gpsSeconds,  duration );

  if ( vrbflg ) fprintf( stdout, "Writing the injection details to %s\n",
      fname);

  /* open the xml file */
  memset( &xmlfp, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlfp, fname), &status );

  /* write the process table */
  snprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "H1H2L1" );
  XLALGPSTimeNow(&(proctable.processTable->end_time));
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, process_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, proctable, 
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
  free( proctable.processTable );

  /* free the unused process param entry */
  this_proc_param = procparams.processParamsTable;
  procparams.processParamsTable = procparams.processParamsTable->next;
  free( this_proc_param );

  /* write the process params table */
  if ( procparams.processParamsTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, process_params_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, procparams, 
          process_params_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
    while( procparams.processParamsTable )
    {
      this_proc_param = procparams.processParamsTable;
      procparams.processParamsTable = this_proc_param->next;
      free( this_proc_param );
    }
  }

  /* write the sim_inspiral table */
  if ( inspInjections.simInspiralTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, sim_inspiral_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, inspInjections, 
          sim_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
  }

  while ( inspInjections.simInspiralTable )
  {
    inj = inspInjections.simInspiralTable;
    inspInjections.simInspiralTable = inspInjections.simInspiralTable->next;
    LALFree( inj );
  }

   
  /* write the sim_ringdown table */
  if ( ringInjections.simRingdownTable )
  {
    LALBeginLIGOLwXMLTable( &status, &xmlfp, sim_ringdown_table );
    LALWriteLIGOLwXMLTable( &status, &xmlfp, ringInjections,
        sim_ringdown_table );
    LALEndLIGOLwXMLTable ( &status, &xmlfp );
  }
  
  while ( ringInjections.simRingdownTable )
  {
    ringList=ringInjections.simRingdownTable;
    ringInjections.simRingdownTable = ringInjections.simRingdownTable->next;
    LALFree( ringList );
  }


  /* close the injection file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlfp ), &status );

  /* destroy random parameters */
  XLALDestroyRandomParams( randParams );

  /* check for memory leaks and exit */
  LALCheckMemoryLeaks();
  return 0;
}
