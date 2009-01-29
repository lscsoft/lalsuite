/*
*  Copyright (C) 2007 Matt Pitkin
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
 * File Name: generic_ring.c
 * 
 * Author: Goggin, L. M., and Brown, D. A. (hacks by Pitkin, M. D.)
 *
 * Hacked version of rinj.c to not be specific for black-hole ring-downs.
 * It just takes in a range of amplitudes, frequencies and quailty 
 * factors (instead of black-hole parameters) and produces a set of 
 * injections based on these. In the future this should probably just
 * be merged with ring.c as another set of input options.
 *
 * $Id$
 *
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <config.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <glob.h>
#include <lalapps.h>
#include <processtable.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/Ring.h>
#include <lal/lalGitID.h>
#include <lalappsGitID.h>

RCSID( "$Id$" );
#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"

#define PROGRAM_NAME "generic_rinj"


#define USAGE \
"lalapps_generic_rinj [options]\n"\
"\nDefaults are shown in brackets\n\n" \
"  --help                   display this message\n"\
"  --verbose                turn verbose flag on\n"\
"  --playground             consider only playground data\n"\
"  --gps-start-time TIME    start injections at GPS time TIME (793130413)\n"\
"  --gps-end-time TIME      end injections at GPS time TIME (795679213)\n"\
"  --time-step STEP         space injections by STEP / pi seconds apart (2630)\n"\
"  --time-interval TIME     distribute injections in interval TIME (0)\n"\
"  --seed SEED              seed random number generator with SEED (1)\n"\
"  --user-tag STRING        set the usertag to STRING\n"\
"  --minimum-quality MIN    set minimum quality factor to MIN (1000)\n"\
"  --maximum-quality MAX    set maximum quality factor to MAX (10000)\n"\
"  --minimum-frequency MIN  set minimum frequency to MIN (1000)\n"\
"  --maximum-frequency MAX  set maximum frequency to MAX (4000)\n"\
"  --minimum-amplitude MIN  set minimum amplitude to MIN (1e-26)\n"\
"  --maximum-amplitude MAX  set maximum amplitude to MAX (1e-20)\n"\
"  --waveform WVF           set the injection waveform to WVF\n"\
"\n"

extern int vrbflg;
int plygnd=0;
ProcessParamsTable *next_process_param( const char *name, const char *type,
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
    LALSnprintf( pp->param, LIGOMETA_PARAM_MAX, "-userTag" );
  else
    LALSnprintf( pp->param, LIGOMETA_PARAM_MAX, "--%s", name );
  strncpy( pp->type, type, LIGOMETA_TYPE_MAX );
  va_start( ap, fmt );
  LALVsnprintf( pp->value, LIGOMETA_VALUE_MAX, fmt, ap );
  va_end( ap );
  return pp;
}

static int is_outside_playground( LIGOTimeGPS gpsStartTime ); 
static REAL8 step( LIGOTimeGPS gpsStartTime );

int main( int argc, char *argv[] )
{
  LALStatus             status = blank_status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;
  const INT4            S4StartTime = 793130413; /* Tues Feb 22, 2005 at 12:00 CST (10:00 PST) */
  const INT4            S4StopTime  = 795679213; /* Wed March 23, 2005 at 24:00 CST (22:00 PST) */

  /* command line options */
  LIGOTimeGPS   gpsStartTime;
  LIGOTimeGPS   gpsEndTime;
  REAL8         meanTimeStep = 2630 / LAL_PI;
  REAL8         timeInterval = 0;
  REAL8         tstep = 0;
  UINT4         randSeed = 1;
  CHAR         *userTag = NULL;
  REAL4         minFreq = 1000.;
  REAL4         maxFreq = 4000.;
  REAL4         minQuality = 1000.;
  REAL4         maxQuality = 10000.;
  REAL4         minAmp = 1.e-26;
  REAL4         maxAmp = 1.e-20;
   
  /* program variables */
  RandomParams *randParams = NULL;
  REAL4  u;
  LALMSTUnitsAndAcc     gmstUnits = { MST_HRS, LALLEAPSEC_STRICT };
  LALGPSandAcc          gpsAndAcc;
  SkyPosition           skyPos;
  LALSource             source;
  LALPlaceAndGPS        placeAndGPS;
  DetTimeAndASource     detTimeAndSource;
  LALDetector           lho = lalCachedDetectors[LALDetectorIndexLHODIFF];
  LALDetector           llo = lalCachedDetectors[LALDetectorIndexLLODIFF];
  LALDetAndSource       detAndSource;
  LALDetAMResponse      resp;
  REAL8                 time_diff_ns;
  REAL4                 splus, scross, cosiota;
  
  /* waveform */
  CHAR waveform[LIGOMETA_WAVEFORM_MAX];
  CHAR coordinates[LIGOMETA_COORDINATES_MAX];  
  
  LALGPSCompareResult        compareGPS;

  /*  xml output data */
  CHAR                  fname[256];
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         injections;
  ProcessParamsTable   *this_proc_param;
  SimRingdownTable     *this_inj = NULL; 
  LIGOLwXMLStream       xmlfp;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"help",                    no_argument,       0,                'h'},
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"playground",              no_argument,       0,                'x'},
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"time-step",               required_argument, 0,                't'},
    {"time-interval",           required_argument, 0,                'i'},
    {"seed",                    required_argument, 0,                's'},
    {"minimum-frequency",       required_argument, 0,                'A'},
    {"maximum-frequency",       required_argument, 0,                'B'},
    {"minimum-quality",         required_argument, 0,                'P'},
    {"maximum-quality",         required_argument, 0,                'Q'},
    {"minimum-amplitude",       required_argument, 0,                'V'},
    {"maximum-amplitude",       required_argument, 0,                'W'},
    {"debug-level",             required_argument, 0,                'z'},
    {"waveform",                required_argument, 0,                'w'},
    {"coordinates",             required_argument, 0,                'c'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {0, 0, 0, 0}
  };
  int c;
  
  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );

  gpsStartTime.gpsSeconds = S4StartTime;
  gpsEndTime.gpsSeconds   = S4StopTime;
  gpsStartTime.gpsNanoSeconds = 0;
  gpsEndTime.gpsNanoSeconds   = 0;

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) 
    calloc( 1, sizeof(ProcessTable) );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->start_time),
        &accuracy ), &status );
  if (strcmp(CVS_REVISION, "$Revi" "sion$"))
  {
    XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME,
        CVS_REVISION, CVS_SOURCE, CVS_DATE, 0);
  }
  else
  {
    XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME,
        lalappsGitCommitID, lalappsGitGitStatus, lalappsGitCommitDate, 0);
  }
  LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );
     
  /* clear the waveform field */
  memset( waveform, 0, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR) );
  memset( coordinates, 0, LIGOMETA_COORDINATES_MAX * sizeof(CHAR) );
    
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
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "a:A:b:B:hP:Q:s:t:V:W:vz:Z:", long_options, &option_index );

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
        if ( gpsinput > 999999999 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is after " 
              "Sep 14, 2011  01:46:26 UTC:\n"
              "(%ld specified)\n", 
              long_options[option_index].name, gpsinput );
          exit( 1 );
        }
        gpsStartTime.gpsSeconds = gpsinput;

        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "int", 
              "%ld", gpsinput );
        break;

      case 'b':
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
        if ( gpsinput > 999999999 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is after " 
              "Sep 14, 2011  01:46:26 UTC:\n"
              "(%ld specified)\n", 
              long_options[option_index].name, gpsinput );
          exit( 1 );
        }
        gpsEndTime.gpsSeconds = gpsinput;
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
        tstep = (REAL8) atof( optarg );
        meanTimeStep = tstep / LAL_PI;
        if ( meanTimeStep <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "time step must be > 0: (%le seconds specified)\n",
              long_options[option_index].name, meanTimeStep );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "float", 
              "%le", meanTimeStep );
        break;
      
      case 'i':
        timeInterval = atof( optarg );
        if ( timeInterval < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "time interval must be >= 0: (%le seconds specified)\n",
              long_options[option_index].name, meanTimeStep );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%le", timeInterval );
        break;

      case 'A':
        minFreq = (REAL4) atof( optarg );
        if ( minFreq <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "miniumum frequency must be > 0: "
              "(%f Hz specified)\n",
              long_options[option_index].name, minFreq );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", minFreq );
        break;

      case 'B':
        maxFreq = (REAL4) atof( optarg );
        if ( maxFreq <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "maxiumum frequency must be > 0: "
              "(%f Hz specified)\n",
              long_options[option_index].name, maxFreq );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", maxFreq );
        break;

      case 'P':
        /* minimum quality factor */
        minQuality = (REAL4) atof( optarg );
        if ( minQuality < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "the minimum quality factor must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, minQuality );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", minQuality );
        break;

      case 'Q':
        /* maximum quality factor */
        maxQuality = (REAL4) atof( optarg );
        if ( maxQuality < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "the maximum quality factor must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, maxQuality );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", maxQuality );
        break;
       
      case 'V':
        /* minimum amplitude */
        minAmp = (REAL4) atof( optarg );
        if ( minAmp < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "minimum amplitude must be >= 0: "
              "(%f specified)\n",
              long_options[option_index].name, minAmp );
              exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", minAmp );
        break;

      case 'W':
        /* max amplitude */
        maxAmp = (REAL4) atof( optarg );
        if ( maxAmp <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "maximum amplitude must be greater than 0: "
              "(%f specified)\n",
              long_options[option_index].name, maxAmp );
              exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
        next_process_param( long_options[option_index].name, 
            "float", "%e", maxAmp );
        break;

      case 'Z':
       /* create storage for the usertag */
        optarg_len = strlen( optarg ) + 1;
        userTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( userTag, optarg, optarg_len );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "string", "%s", optarg );
        break;
      
      case 'w':
        LALSnprintf( waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "%s",
            optarg);
        this_proc_param = this_proc_param->next =
           next_process_param( long_options[option_index].name, "string",
              "%s", optarg);
        break;
      
      case 'c':
        LALSnprintf( coordinates, LIGOMETA_COORDINATES_MAX * sizeof(CHAR), "%s",
            optarg);
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "%s", optarg);
        break;
                                
      case 'v':
        vrbflg = 1;
        break;
      
      case 'x':
        plygnd = 1;
        break;
        
      case 'z':
        set_debug_level( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
            "string", "%s", optarg );
        break;

      case 'h':
        fprintf( stderr, USAGE );
        exit( 0 );
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

  if ( !*waveform )
    {
      /* use Ringdown as the default waveform */
      LALSnprintf( waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR),
          "Ringdown");
      }
  
  if ( !*coordinates )
        {
          /* use equatorial as the default system */
          LALSnprintf( coordinates, LIGOMETA_COORDINATES_MAX * sizeof(CHAR),
                              "EQUATORIAL");
                      }
  
                      
  /* check freq, Q and Amp ranges are valid */
  if( maxFreq < minFreq ){
    fprintf(stderr, "maximum frequency is less than minimum frequency.\n");
    exit( 1 );
  }
  
  if( maxQuality < minQuality ){
    fprintf(stderr, "maximum quality factor is less than minimum quality.\n");
    exit( 1 );
  }
  
  if( maxAmp < minAmp ){
    fprintf(stderr, "maximum amplitude is less than minimum amplitude.\n");
    exit( 1 );
  }
  
  /*
   *
   * initialization
   *
   */


  /* initialize the random number generator */
  LAL_CALL( LALCreateRandomParams( &status, &randParams, randSeed ), &status );
  
  /* null out the head of the linked list */
  injections.simRingdownTable = NULL;

  /* create the output file name */
  if ( userTag )
  {
    LALSnprintf( fname, sizeof(fname), "HL-INJECTIONS_%d_%s-%d-%d.xml", 
        randSeed, userTag, gpsStartTime.gpsSeconds, 
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else
  {
    LALSnprintf( fname, sizeof(fname), "HL-INJECTIONS_%d-%d-%d.xml", 
        randSeed, gpsStartTime.gpsSeconds, 
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }

  /* check that the start time is before the end time */
  LAL_CALL( LALCompareGPS( &status, &compareGPS, &gpsStartTime, &gpsEndTime ),
      &status );
  
  /* check if gps_start_time is in playground */
  /* if it is not then move it to a random place in the next playground
   * interval */
  if ( plygnd )
  {
    if ( vrbflg )
             fprintf( stdout, "injecting into playground only\n");
    if ( is_outside_playground( gpsStartTime ) )
    {
      LALStatus status = blank_status;
      REAL8 randstep;
      if ( vrbflg )
        fprintf( stdout, "gps-start-time outside of playground, ... shifting to next playground interval\n");
/*      REAL8 nextplay;*/
/*      nextplay = step(gpsStartTime);*/  /*find start of next playground interval*/
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      randstep = step(gpsStartTime) + u * 600;  /* find a random time within this 600s */
      LAL_CALL( LALAddFloatToGPS( &status, &gpsStartTime, &gpsStartTime,
            randstep ), &status );  /* add this to the old gpsStartTime */
          
    }
  }

  
  /*
   *
   * loop over duration of desired output times
   *
   */

  while ( compareGPS == LALGPS_EARLIER )
  {
    /* create the sim_ringdown table */
    if ( injections.simRingdownTable )
      {
        this_inj = this_inj->next = (SimRingdownTable *)
          LALCalloc( 1, sizeof(SimRingdownTable) );
        }
    else
      {
        injections.simRingdownTable = this_inj = (SimRingdownTable *)
          LALCalloc( 1, sizeof(SimRingdownTable) );
        }

    /* set the waveform and coordinates fields */
    memcpy( this_inj->waveform, waveform, LIGOMETA_WAVEFORM_MAX *
        sizeof(CHAR));
    memcpy( this_inj->coordinates, coordinates, LIGOMETA_COORDINATES_MAX *
                sizeof(CHAR));
    
    this_inj->epsilon = 0.;
    /* set the geocentric start time of the injection */
    /* XXX CHECK XXX */
    
     this_inj->geocent_start_time = gpsStartTime;
    if ( timeInterval )
    {
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      LAL_CALL( LALAddFloatToGPS( &status, &(this_inj->geocent_start_time),
          &(this_inj->geocent_start_time), u * timeInterval ), &status );
    }    

    /* set all black-hole parameters to zero as not needed */
    this_inj->mass = 0.;
    this_inj->spin = 0.;
    
    /* set frequency, f0, and quality factor Q */
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->frequency = minFreq + u * (maxFreq - minFreq);
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->quality = minQuality + u * (maxQuality - minQuality);
            
    /* spatial distribution */
    /* compute random longitude and latitude */ 
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->longitude = LAL_TWOPI * u ;
    
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->latitude = asin( 2.0 * u - 1.0 ) ;
   
    /* initial phase */
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->phase = u * LAL_TWOPI;
    
    /* set distance to zero */
    this_inj->distance = 0.;

    /* compute random inclination, polarization */
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->inclination = acos( 2.0 * u - 1.0 );
    
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->polarization = LAL_TWOPI * u ;
    
    gpsAndAcc.gps = this_inj->geocent_start_time;

    /* set gmst */
    LAL_CALL( LALGPStoGMST1( &status, &(this_inj->start_time_gmst),
          &(this_inj->geocent_start_time), &gmstUnits ), &status);
    
    memset( &skyPos, 0, sizeof(SkyPosition) );
    memset( &source, 0, sizeof(LALSource) );
    memset( &placeAndGPS, 0, sizeof(LALPlaceAndGPS) );
    memset( &detTimeAndSource, 0, sizeof(DetTimeAndASource) );
    memset( &detAndSource, 0, sizeof(LALDetAndSource) );

    skyPos.longitude = this_inj->longitude;
    skyPos.latitude  = this_inj->latitude;
    skyPos.system    = COORDINATESYSTEM_EQUATORIAL;

    source.equatorialCoords = skyPos;
    source.orientation      = this_inj->polarization;
    
    placeAndGPS.p_gps = &(this_inj->geocent_start_time);
    
    detTimeAndSource.p_det_and_time = &placeAndGPS;
    detTimeAndSource.p_source = &skyPos;
    detAndSource.pSource = &source;
                    
    gpsAndAcc.accuracy = LALLEAPSEC_STRICT;
    gpsAndAcc.gps = this_inj->geocent_start_time;
    
    /* calculate h0 */
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->amplitude = minAmp + u * (maxAmp - minAmp);
      
    /* calculate hrss */
    this_inj->hrss = this_inj->amplitude * sqrt( 2 / LAL_PI / this_inj->frequency ) * 
      pow( ( 2.0 * pow( this_inj->quality, 3.0 ) + this_inj->quality ) / 
          ( 1.0 + 4.0 * pow ( this_inj->quality, 2 ) ) , 0.5);
      
    /* initialize end times with geocentric value */
    this_inj->h_start_time = this_inj->l_start_time = this_inj->geocent_start_time;
    
    /* initialize distances with real distance and compute splus and scross*/
    this_inj->eff_dist_h = this_inj->eff_dist_l = 2.0 * this_inj->distance;
    cosiota = cos( this_inj->inclination );
    splus = -( 1.0 + cosiota * cosiota );
    scross = -2.0 * cosiota; 
      
    /* lho */
    placeAndGPS.p_detector = &lho;
    LAL_CALL( LALTimeDelayFromEarthCenter( &status, &time_diff_ns,
          &detTimeAndSource ), &status );
    LAL_CALL( LALAddFloatToGPS( &status, &(this_inj->h_start_time),
          &(this_inj->h_start_time), time_diff_ns ), &status );

    /* compute the response of the LHO detectors */
    detAndSource.pDetector = &lho;
    LAL_CALL( LALComputeDetAMResponse( &status, &resp, &detAndSource,
          &gpsAndAcc ), &status );
    
    /* compute the effective distance for LHO */
    this_inj->eff_dist_h /= sqrt( splus*splus*resp.plus*resp.plus +
        scross*scross*resp.cross*resp.cross );

    /* fprintf( stdout, "antenna factor = %e\n",sqrt( splus*splus*resp.plus*resp.plus +
        scross*scross*resp.cross*resp.cross ) );    */

    /* compute hrss at LHO */ 
    this_inj->hrss_h = this_inj->amplitude * pow ( ( 
          (2*pow(this_inj->quality,3)+this_inj->quality ) * splus*splus*resp.plus*resp.plus +
          2*pow(this_inj->quality,2) * splus*scross*resp.plus*resp.cross +
          2*pow(this_inj->quality,3) * scross*scross*resp.cross*resp.cross )
        /  2.0 / LAL_PI / this_inj->frequency / ( 1.0 + 4.0 * pow ( this_inj->quality, 2 ) ) , 0.5 );
      
    /* llo */
    placeAndGPS.p_detector = &llo;
    LAL_CALL( LALTimeDelayFromEarthCenter( &status,  &time_diff_ns,
          &detTimeAndSource ), &status);
    LAL_CALL( LALAddFloatToGPS( &status,  &(this_inj->l_start_time),
          &(this_inj->l_start_time), time_diff_ns ), &status);

    /* compute the response of the LLO detector */
    detAndSource.pDetector = &llo;
    LAL_CALL( LALComputeDetAMResponse( &status, &resp, &detAndSource,
          &gpsAndAcc ), &status);
    
    /* compute the effective distance for LLO */
    this_inj->eff_dist_l /= sqrt( splus*splus*resp.plus*resp.plus 
        + scross*scross*resp.cross*resp.cross );
    
    /* compute hrss at LLO */
    this_inj->hrss_l = this_inj->amplitude * pow ( (
          (2*pow(this_inj->quality,3)+this_inj->quality ) * splus*splus*resp.plus*resp.plus +
          2*pow(this_inj->quality,2) * splus*scross*resp.plus*resp.cross +
          2*pow(this_inj->quality,3) * scross*scross*resp.cross*resp.cross )
          /  2.0 / LAL_PI / this_inj->frequency / ( 1.0 + 4.0 * pow ( this_inj->quality, 2 ) ) , 0.5 );
        
    /* increment the injection time */
    if ( plygnd )
    {
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      meanTimeStep = step(gpsStartTime) + u * 600;
    }
    LAL_CALL( LALAddFloatToGPS( &status, &gpsStartTime, &gpsStartTime, 
          meanTimeStep ), &status );
    LAL_CALL( LALCompareGPS( &status, &compareGPS, &gpsStartTime, 
          &gpsEndTime ), &status );


  } /* end loop over injection times */

  
  /* destroy random parameters */
  LAL_CALL( LALDestroyRandomParams( &status, &randParams ), &status );
   
  /*
   *
   * write output to LIGO_LW XML file
   *
   */

  /* open the xml file */
  memset( &xmlfp, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlfp, fname), &status );

  /* write the process table */
  LALSnprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "H1H2L1" );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->end_time),
        &accuracy ), &status );
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

  /* write the sim_ringdown table */
  if ( injections.simRingdownTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, sim_ringdown_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, injections, 
          sim_ringdown_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
  }
  while ( injections.simRingdownTable )
  {
    this_inj = injections.simRingdownTable;
    injections.simRingdownTable = injections.simRingdownTable->next;
    LALFree( this_inj );
  }
  /* close the injection file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlfp ), &status );

  /* check for memory leaks and exit */
  LALCheckMemoryLeaks();
  return 0;

}

/* function to check if time outside playground */
  static int is_outside_playground(  LIGOTimeGPS gpsStartTime )
{
  int ans = 0;
  INT8 startTime = 0;
  startTime = XLALGPStoINT8 ( &gpsStartTime ) ;
  if ( ( ( startTime - 729273613LL * 1000000000LL ) % ( 6370 * 1000000000LL ) ) > ( 600 * 1000000000LL ) )
    ans = 1;
  return ans;
}


/* function to calculate the time interval to the start of next 600s of playground data */
static REAL8 step( LIGOTimeGPS gpsStartTime )
{
  INT8 istep;
  REAL8 fstep;
  INT8 startTime = 0;
  startTime = XLALGPStoINT8 ( &gpsStartTime ) ;
  istep = ( 6370 * 1000000000LL ) - ( ( startTime - 729273613LL * 1000000000LL  ) % ( 6370 * 1000000000LL ) );
  fstep = (REAL8) istep; /*convert integer to float*/
  fstep = fstep / 1000000000 ; /*convert ns to s */
  return fstep;
}

