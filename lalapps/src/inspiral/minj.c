/*
*  Copyright (C) 2007 Duncan Brown
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
 * File Name: minj.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <config.h>

#if !defined HAVE_GSL_GSL_FFT_REAL_H || !defined HAVE_LIBGSL
int main( void ) { fprintf( stderr, "no gsl: disabled\n" ); return 77; }
#else

#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <lalapps.h>
#include <processtable.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
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
#include <lal/lalGitID.h>
#include <lalappsGitID.h>

RCSID( "$Id$" );
#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "minj"

#define USAGE \
"lalapps_minj [options]\n"\
"\nDefaults are shown in brackets\n\n" \
"  --help                   display this message\n"\
"  --verbose                print mass and galactocentic cartesian coordinates\n"\
"  --gps-start-time TIME    start injections at GPS time TIME (729273613)\n"\
"  --gps-end-time TIME      end injections at GPS time TIME (734367613)\n"\
"  --time-step STEP         space injections by ave of STEP sec (2630/PI)\n"\
"  --time-interval TIME     distribute injections in interval TIME (0)\n"\
"  --seed SEED              seed random number generator with SEED (1)\n"\
"  --user-tag STRING        set the usertag to STRING\n"\
"  --minimum-mass MIN       set the minimum componenet mass to MIN (0.1)\n"\
"  --maximum-mass MAX       set the maximum componenet mass to MAX (1.0)\n"\
"  --core-radius A          set the galactic core radius to A kpc (8.5)\n"\
"  --flatten-halo Q         set the halo flattening parameter to Q (1)\n"\
"  --halo-radius RMAX       set the maximum halo radius to RMAX kpc (50)\n"\
"  --write-compress         write a compressed xml file\n"\
"\n"

/* all units are in kpc since this is what GalacticInspiralParamStruc expects */

extern int vrbflg;
struct halo_pdf_params
{
  double u;
  double a;
};


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
    snprintf( pp->param, LIGOMETA_PARAM_MAX, "-userTag" );
  else
    snprintf( pp->param, LIGOMETA_PARAM_MAX, "--%s", name );
  strncpy( pp->type, type, LIGOMETA_TYPE_MAX );
  va_start( ap, fmt );
  vsnprintf( pp->value, LIGOMETA_VALUE_MAX, fmt, ap );
  va_end( ap );
  return pp;
}


double halo_pdf( 
    double      r, 
    void       *p
    )
{
  struct halo_pdf_params *params = (struct halo_pdf_params *)p;
  double u = (params->u);
  double a = (params->a);

  return r - a * atan2(r,a) - u;
}


int main( int argc, char *argv[] )
{
  LALStatus             status = blank_status;
  const INT4            S2StartTime = 729273613; /* Feb 14 2003 16:00:00 UTC */
  const INT4            S2StopTime  = 734367613; /* Apr 14 2003 15:00:00 UTC */

  /* command line options */
  LIGOTimeGPS   gpsStartTime = {S2StartTime, 0};
  LIGOTimeGPS   gpsEndTime   = {S2StopTime, 0};
  REAL8         meanTimeStep = 2630 / LAL_PI;
  REAL8         timeInterval = 0;
  UINT4         randSeed = 1;
  CHAR         *userTag = NULL;
  REAL4         minMass = 0.1;
  REAL4         maxMass = 1.0;
  REAL4         r_core = 5.0;     /* kpc core radius */
  REAL4         r_max = 50.0;     /* kpc halo radius */
  REAL4         q = 1.0;          /* flatten halo */

  /* program variables */
  RandomParams *randParams = NULL;
  REAL4  u;
  REAL4  deltaM;

  int i, stat;
  const int maxIter = 1000;
  const double tol = 1e-6;
  const gsl_root_fsolver_type *solver_type;
  gsl_root_fsolver *solver;
  gsl_function pdf;

  struct halo_pdf_params pdf_params;

  double r_lo, r_hi, r;
  double cosphi, sinphi;
  double pdf_norm;

  GalacticInspiralParamStruc galacticPar;
  LALGPSCompareResult        compareGPS;

  /* xml output data */
  CHAR                  fname[256];
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         injections;
  ProcessParamsTable   *this_proc_param;
  SimInspiralTable     *this_inj = NULL;
  LIGOLwXMLStream       xmlfp;
  UINT4                 outCompress = 0;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"help",                    no_argument,       0,                'h'},
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"write-compress",          no_argument,       &outCompress,      1 },
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"time-step",               required_argument, 0,                't'},
    {"time-interval",                required_argument, 0,                     'i'},
    {"seed",                    required_argument, 0,                's'},
    {"minimum-mass",            required_argument, 0,                'A'},
    {"maximum-mass",            required_argument, 0,                'B'},
    {"core-radius",             required_argument, 0,                'p'},
    {"flatten-halo",            required_argument, 0,                'q'},
    {"halo-radius",             required_argument, 0,                'r'},
    {"debug-level",             required_argument, 0,                'z'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {0, 0, 0, 0}
  };
  int c;

  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) 
    calloc( 1, sizeof(ProcessTable) );
  XLALGPSTimeNow(&(proctable.processTable->start_time));
  if (strcmp(CVS_REVISION,"$Revi" "sion$"))
    {
      LAL_CALL( populate_process_table( &status, proctable.processTable, 
                                        PROGRAM_NAME, CVS_REVISION,
                                        CVS_SOURCE, CVS_DATE ), &status );
    }
  else
    {
      LAL_CALL( populate_process_table( &status, proctable.processTable, 
                                        PROGRAM_NAME, lalappsGitCommitID,
                                        lalappsGitGitStatus,
                                        lalappsGitCommitDate ), &status );
    }
  snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );
  

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
        "a:A:b:B:hi:p:q:r:s:t:vz:Z:", long_options, &option_index );

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
        meanTimeStep = (REAL8) atof( optarg );
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
        minMass = (REAL4) atof( optarg );
        if ( minMass <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "miniumum component mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, minMass );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", minMass );
        break;

      case 'B':
        maxMass = (REAL4) atof( optarg );
        if ( maxMass <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "maxiumum component mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, maxMass );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", maxMass );
        break;

      case 'p':
        /* core-radius */
        r_core = (REAL4) atof( optarg );
        if ( r_core <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "galactic core radius must be > 0: "
              "(%f kpc specified)\n",
              long_options[option_index].name, r_core );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", r_core );
        break;

      case 'q':
        /* flatten-halo */
        q = (REAL4) atof( optarg );
        if ( q <= 0 || q > 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "halo flattening parameter must be in range (0,1]: "
              "(%f specified)\n",
              long_options[option_index].name, q );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", q );
        break;

      case 'r':
        /* max halo radius */
        r_max = (REAL4) atof( optarg );
        if ( r_max <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "halo radius must be greater than 0: "
              "(%f kpc specified)\n",
              long_options[option_index].name, r_max );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", r_max );
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

      case 'v':
        vrbflg = 1;
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


  /*
   *
   * initialization
   *
   */


  /* initialize the random number generator */
  LAL_CALL( LALCreateRandomParams( &status, &randParams, randSeed ), &status );

  /* initialize the gsl solver with the spatial pdf */
  pdf.function = &halo_pdf;
  pdf.params   = &pdf_params;
  solver_type = gsl_root_fsolver_bisection;
  solver = gsl_root_fsolver_alloc( solver_type );
  pdf_params.a = r_core;

  /* normalization for the spatial pdf */
  pdf_norm = r_max - r_core * atan2( r_max, r_core );

  /* mass range */
  deltaM = maxMass - minMass;

  /* null out the head of the linked list */
  injections.simInspiralTable = NULL;

  /* create the output file name */
  if ( userTag && outCompress )
  {
    snprintf( fname, sizeof(fname), "HL-INJECTIONS_%d_%s-%d-%d.xml.gz",
        randSeed, userTag, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else if ( userTag && !outCompress )
  {
    snprintf( fname, sizeof(fname), "HL-INJECTIONS_%d_%s-%d-%d.xml", 
        randSeed, userTag, gpsStartTime.gpsSeconds, 
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else if ( !userTag && outCompress )
  {
    snprintf( fname, sizeof(fname), "HL-INJECTIONS_%d-%d-%d.xml.gz",
        randSeed, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else
  {
    snprintf( fname, sizeof(fname), "HL-INJECTIONS_%d-%d-%d.xml", 
        randSeed, gpsStartTime.gpsSeconds, 
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }

  /* check that the start time is before the end time */
  LAL_CALL( LALCompareGPS( &status, &compareGPS, &gpsStartTime, &gpsEndTime ),
      &status );


  /*
   *
   * loop over duration of desired output times
   *
   */


  while ( compareGPS == LALGPS_EARLIER )
  {
    /* uniformly distributed masses */
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    galacticPar.m1 = minMass + u * deltaM;
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    galacticPar.m2 = minMass + u * deltaM;
    
    /* spatial distribution */
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    pdf_params.u = u * pdf_norm;

    r_lo = 0.0;
    r_hi = r_max;
    gsl_root_fsolver_set( solver, &pdf, r_lo, r_hi );

    for ( i = 0; i < maxIter; ++i )
    {
      gsl_root_fsolver_iterate( solver );
      r = gsl_root_fsolver_root( solver );
      r_lo = gsl_root_fsolver_x_lower( solver );
      r_hi = gsl_root_fsolver_x_upper( solver );
      stat = gsl_root_test_interval( r_lo, r_hi, 0, tol );
      if ( stat == GSL_SUCCESS )
        break;
    }

    if ( stat != GSL_SUCCESS )
    {
      fprintf( stderr, "could not find root after %d iterations\n", maxIter );
      exit( 1 );
    }

    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    sinphi = 2.0 * u - 1.0;
    cosphi = sqrt( 1.0 - sinphi*sinphi );
    
    galacticPar.rho = r * cosphi;
    galacticPar.z   = q * r * sinphi;

    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    galacticPar.lGal = LAL_TWOPI * u;    

    if ( vrbflg ) fprintf( stdout, "%e %e %e %e %e\n", 
        galacticPar.m1, galacticPar.m2,
        galacticPar.rho * cos( galacticPar.lGal ),
        galacticPar.rho * sin( galacticPar.lGal ),
        galacticPar.z );

    /* create the sim_inspiral table */
    if ( injections.simInspiralTable )
    {
      this_inj = this_inj->next = (SimInspiralTable *)
        LALCalloc( 1, sizeof(SimInspiralTable) );
    }
    else
    {
      injections.simInspiralTable = this_inj = (SimInspiralTable *)
        LALCalloc( 1, sizeof(SimInspiralTable) );
    }

    /* set the geocentric end time of the injection */
    galacticPar.geocentEndTime = gpsStartTime;
    if ( timeInterval )
    {
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      XLALGPSAdd( &(galacticPar.geocentEndTime), u * timeInterval );
    }

    /* populate the sim_inspiral table */
    LAL_CALL( LALGalacticInspiralParamsToSimInspiralTable( &status,
          this_inj, &galacticPar, randParams ), &status );

    /* set the source and waveform fields */
    snprintf( this_inj->source, LIGOMETA_SOURCE_MAX * sizeof(CHAR), "MW" );
    snprintf( this_inj->waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), 
        "GeneratePPNtwoPN" );

    /* increment the injection time */
    XLALGPSAdd( &gpsStartTime, meanTimeStep );
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
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlfp, fname ), &status );

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
  if ( injections.simInspiralTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, sim_inspiral_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, injections, 
          sim_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
  }
  while ( injections.simInspiralTable )
  {
    this_inj = injections.simInspiralTable;
    injections.simInspiralTable = injections.simInspiralTable->next;
    LALFree( this_inj );
  }

  /* close the injection file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlfp ), &status );

  /* check for memory leaks and exit */
  LALCheckMemoryLeaks();
  return 0;
}
#endif
