/*
*  Copyright (C) 2007 Alexander Dietz, Drew Keppel, Duncan Brown, Eirini Messaritaki, Patrick Brady, Anand Sengupta, Stephen Fairhurst, Thomas Cokelaer
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
 * File Name: bbhinj.c
 *
 * Author: Brown, D. A., Messaritaki E.
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
 * \c lalapps_bbhinj --- produces inspiral injection data files.</dd>
 *
 * <dt>Synopsis</dt><dd>
 * \c lalapps_bbhinj
 * [<tt>--help</tt>]
 * [<tt>--verbose</tt>]
 * [<tt>--gps-start-time</tt> \c tstart]
 * [<tt>--gps-end-time</tt> \c tend]
 * [<tt>--time-step</tt> \c tstep]
 * [<tt>--time-interval</tt> \c tinterval]
 * [<tt>--seed</tt> \c seed]
 * [<tt>--usertag</tt> \c tag]
 * [<tt>--min-mass</tt> \c mmin]
 * [<tt>--max-mass</tt> \c mmax]
 * [<tt>--max-total-mass</tt> \c totalmass]
 * [<tt>--min-distance</tt> \c dmin]
 * [<tt>--max-distance</tt> \c dmax]
 * [<tt>--m-distr</tt> \c mdistr]
 * [<tt>--d-distr</tt> \c ddistr]
 * [<tt>--waveform</tt> \c wvf] </dd>
 *
 * <dt>Description</dt><dd>
 * \c lalapps_bbhinj
 * generates a number of inspiral  parameters suitable  for using in a Monte
 * Carlo injection to test the efficiency of an inspiral search.  The  various
 * parameters (detailed  below)  are randomly chosen and are appropriate for
 * a population of binaries that extends over all space between
 * the minimum and maximum distances specified.
 * Despite its name, it can be used for BNS and for BBH parameter generation.
 *
 * The output of this program  is  a  list  of  the  injected events,  starting
 * at  the specified start time and ending at the specified end time.  One
 * injection with random inspiral parameters will be made every specified time
 * step, and will be randomly placed within the specified time interval.
 * The output is written to a file name in the standard inspiral pipeline format:
 *
 * \code
 * HL-INJECTIONS_USERTAG_SEED-GPSSTART-DURATION.xml
 * \endcode
 *
 * where \c USERTAG is \c tag as specfied on the command line,
 * \c SEED is the  value  of  the random number seed chosen and
 * \c GPSSTART and \c DURATION describes the GPS time interval that
 * the file covers. The file is in the standard LIGO lightweight XML format
 * containing a \c sim_inspiral table that describes the injections.
 * In addition, an ASCII log file called <tt>injlog.txt</tt> is also written.
 * If a <tt>--user-tag</tt> is not specified on the command line, the
 * \c _USERTAG part of the filename will be omitted.
 *
 * The GEO, TAMA and VIRGO end times and effective distances are currently
 * not being populated.</dd>
 *
 * <dt>Options</dt><dd>
 * <ul>
 * <li><tt>--help</tt>: Print a help message.</li>
 *
 * <li><tt>--gps-start-time</tt> \c tstart:
 * Optional.  Start time of the injection data to be created. Defaults to the
 * start of S2, Feb 14 2003 16:00:00 UTC (GPS time 729273613)</li>
 *
 * <li><tt>--gps-end-time</tt> \c tend:
 * Optional. End time of the injection data to be created. Defaults to the end of
 * S2, Apr 14 2003 15:00:00 UTC (GPS time 734367613).</li>
 *
 * <li><tt>--time-step</tt> \c tstep:
 * Optional. Sets the time step interval between injections. The injections will
 * occur with an average spacing of \c tstep seconds. Defaults to
 * \f$2630/\pi\f$.</li>
 *
 * <li><tt>--time-interval</tt> \c tinterval:
 * Optional. Sets the time interval during which an injection can occur.
 * Injections are uniformly distributed over the interval.  Setting \c tstep
 * to \f$6370\f$ and \c tinterval to 600 guarantees there will be one injection
 * into each playground segment and they will be randomly distributed within the
 * playground times.</li>
 *
 * <li><tt>--seed</tt> \c seed:
 * Optional. Seed the random number generator with the integer \c seed.
 * Defaults to \f$1\f$.</li>
 *
 * <li><tt>--user-tag</tt> \c string: Optional. Set the user tag for this
 * job to be \c string. May also be specified on the command line as
 * <tt>-userTag</tt> for LIGO database compatibility.</li>
 *
 * <li><tt>--min-mass</tt> \c mmin: Optional. The minimum component
 * mass of the binary. If not specified, it is automatically set to 3
 * \f$M_\odot\f$.</li>
 *
 * <li><tt>--max-mass</tt> \c mmax: Optional. The maximum component
 * mass of the binary. If not specified, it is automatically set to 20
 * \f$M_\odot\f$.</li>
 *
 * <li><tt>-max-total-mass</tt> \c totalmass: Optional. The maximum of the total
 * masses of the two components. If not specified, no restriction is set to
 * the total mass of the two components. </li>
 *
 * <li><tt>--min-distance</tt> \c dmin: Optional. The minimum distance
 * (in kpc) from the earth of the signals. If not specified, it is automatically
 * set to 1 kpc.</li>
 *
 * <li><tt>--max-distance</tt> \c dmax: Optional. The maximum distance
 * (in kpc) from the earth of the signals. If not specified, it is automatically
 * set to 20 Mpc (20000 kpc).</li>
 *
 * <li><tt>--m-distr</tt> \c mdistr: Optional. It lets the user specify
 * the mass distribution of the population. The choices are:
 *  <ul>
 *  <li> \c mdistr = 0 the distribution is uniform over total mass
 *  (DEFAULT VALUE)
 *  </li><li> \c mdistr = 1 the distribution is uniform over mass1 and
 *     over mass2
 *  </li></ul></li>
 *
 * <li><tt>--d-distr</tt> \c ddistr: Optional. It lets the used
 * specify the distance distribution of the population. The choices are:
 *  <ul>
 *  <li> \c ddistr = 0 uniform-distance distribution
 *  (DEFAULT VALUE)
 *  </li><li> \c ddistr = 1 uniform \f$\log\f$(distance) distribution
 *  </li><li> \c ddistr = 2 uniform volume distribution
 *  </li></ul></li>
 *
 * <li><tt>--waveform</tt> \c wvf: Optional. Sets the waveforms to
 * be injected to \c wvf. The \c wvf consists of two parts,
 * which are specified as one word. The first part is one of:
 *   <ul>
 *   <li> EOB
 *   </li><li> PadeT1
 *   </li><li> TaylorT1
 *   </li><li> TaylorT3
 *   </li><li> GeneratePPN
 *   </li></ul>
 * The second part is one of:
 *   <ul>
 *   <li> newtonian
 *   </li><li> onePN
 *   </li><li> onePointFivePN
 *   </li><li> twoPN (ONLY CHOICE if GeneratePPN is used before!)
 *   </li><li> twoPointFivePN
 *   </li><li> threePN
 *   </li></ul>
 * If nothing is specified, the default value is EOBtwoPN. It should be noted
 * that if GeneratePPNtwoPN is used as the waveform, the code used to
 * perform the injections is different than otherwise. For GeneratePPNtwoPN, the
 * older injection code (that does only standard post-Newtonian injections) is
 * used. That is the recommended approximant for the case of binary neutron
 * star injections.
 * </li>
 * </ul></dd>
 *
 * <dt>Example:</dt><dd>
 * The following command will produce an injection file for a population
 * of binary black holes of total mass between 6 and 20 \f$M_\odot\f$ (component
 * mass between 3 and 10 \f$M_\odot\f$), with uniform-distance distribution
 * between 100 kpc and 3 Mpc. Since no mass-distribution is specified, the
 * mass-distribution will be uniform over total mass (default value).
 * \code
 * lalapps_bbhinj --min-mass 3.0 --max-mass 10.0 --min-distance 100
 * --max-distance 3000 --d-distr 0 --waveform PadeT1twoPN
 * \endcode</dd>
 *
 * <dt>Author</dt><dd>
 * Jolien Creighton, Patrick Brady, Duncan Brown and Eirini Messaritaki</dd>
 * </dl>
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
#include <time.h>
#include <lalapps.h>
#include <processtable.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <lal/LALStdio.h>
#include <lal/LALgetopt.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateInspiral.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/lalGitID.h>
#include <lalappsGitID.h>

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "bbhinj"

#define USAGE \
"lalapps_bbhinj [options]\n"\
"\nDefaults are shown in brackets\n\n" \
"  --help                   display this message\n"\
"  --version                print version information and exit\n"\
"  --verbose                print mass and galactocentic cartesian coordinates\n"\
"  --write-compress         write a compressed xml file\n"\
"  --f-lower FREQUENCY      lower cut-off frequency.\n"\
"                           If so, LAL code will use 40Hz by default.\n"\
"  --gps-start-time TIME    start injections at GPS time TIME (729273613)\n"\
"  --gps-end-time TIME      end injections at GPS time TIME (734367613)\n"\
"  --time-step STEP         space injections by ave of STEP sec (2630/PI)\n"\
"  --time-interval TIME     distribute injections in interval TIME (0)\n"\
"  --seed SEED              seed random number generator with SEED (1)\n"\
"  --user-tag STRING        set the usertag to STRING\n"\
"  --min-mass MIN           set the minimum component mass to MIN (3.0)\n"\
"  --max-mass MAX           set the maximum component mass to MAX (20.0)\n"\
"  --max-total-mass TOTAL   set the total mass of the two components\n"\
"  --min-distance DMIN      set the minimum distance to DMIN kpc (1)\n"\
"  --max-distance DMAX      set the maximum distance to DMAX kpc (20000)\n"\
"  --d-distr DDISTR         distribute injections uniformly over\n"\
"                           d (DDISTR = 0), or over log10(d) (DDISTR = 1)\n"\
"                           or over volume (DDISTR = 2)\n"\
"                           (default: DDISTR = 0)\n"\
"  --m-distr MDISTR         distribute injections uniformly over\n"\
"                           total mass (MDISTR = 0), or over mass1 and\n"\
"                           over mass2 (MDISTR = 1) (default: MDISTR=0)\n"\
"  --waveform WVF           set the injection waveform to WVF\n"\
"                           (EOB, GeneratePPN, TaylorT1, TaylorT3,PadeT1);\n"\
"                           followed by the order (newtonian, onePN,\n"\
"                           onePointFivePN, twoPN, twoPointFivePN, threePN)\n"\
"                           (default: EOBtwoPN)\n"\
"\n"


#define S2StartTime  729273613 /* Feb 14 2003 16:00:00 UTC */
#define S2StopTime   734367613 /* Apr 14 2003 15:00:00 UTC */
/* all units are in kpc since this is what GalacticInspiralParamStruc expects */
static ProcessParamsTable *next_process_param( 
        const char *name, 
        const char *type,
            const char *fmt, ... );

extern int vrbflg;


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


int main( int argc, char *argv[] )
{
  LALStatus             status = blank_status;
#if 0
  const INT4            S2StartTime = 729273613; /* Feb 14 2003 16:00:00 UTC */
  const INT4            S2StopTime  = 734367613; /* Apr 14 2003 15:00:00 UTC */
#endif
  /* command line options */
  LIGOTimeGPS   gpsStartTime = {S2StartTime, 0};
  LIGOTimeGPS   gpsEndTime   = {S2StopTime, 0};
  REAL8         meanTimeStep = 2630 / LAL_PI;
  REAL8         timeInterval = 0;
  UINT4         randSeed = 1;
  CHAR         *userTag = NULL;
  REAL4         minMass = 3.0;       /* minimum component mass */
  REAL4         maxMass = 20.0;      /* maximum component mass */
  REAL4         sumMaxMass = 0.0;    /* maximum total mass sum */
  UINT4         sumMaxMassUse=0;     /* flag indicating to use the sumMaxMass */
  REAL4         dmin = 1.0;          /* minimum distance from earth (kpc) */
  REAL4         dmax = 20000.0 ;     /* maximum distance from earth (kpc) */
  REAL4         fLower = 0;          /* default value for th lower cut off frequency */
/* REAL4         Rcore = 0.0; */
  UINT4         ddistr = 0, mdistr=0;

  /* program variables */
  RandomParams *randParams = NULL;
  REAL4  u, exponent, d2;
  REAL4  deltaM, mtotal;


  /* waveform */
  CHAR waveform[LIGOMETA_WAVEFORM_MAX];

#if 0
  int i, stat;

  double d, cosphi, sinphi;
#endif

/*  GalacticInspiralParamStruc galacticPar; */

  /* xml output data */
  CHAR                  fname[256];
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         injections;
  ProcessParamsTable   *this_proc_param;
  SimInspiralTable     *this_inj = NULL;
  LIGOLwXMLStream       xmlfp;
  UINT4                 outCompress = 0;

  /* LALgetopt arguments */
  struct LALoption long_options[] =
  {
    {"help",                    no_argument,       0,                'h'},
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"write-compress",          no_argument,       &outCompress,      1 },
    {"version",                 no_argument,       0,                'V'},
    {"f-lower",                        required_argument, 0,                'f'},
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"time-step",               required_argument, 0,                't'},
    {"time-interval",           required_argument, 0,                'i'},
    {"seed",                    required_argument, 0,                's'},
    {"min-mass",                required_argument, 0,                'A'},
    {"max-mass",                required_argument, 0,                'B'},
    {"max-total-mass",          required_argument, 0,                'x'},
    {"min-distance",            required_argument, 0,                'p'},
    {"max-distance",            required_argument, 0,                'r'},
    {"d-distr",                 required_argument, 0,                'd'},
    {"m-distr",                 required_argument, 0,                'm'},
    {"waveform",                required_argument, 0,                'w'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {0, 0, 0, 0}
  };
  int c;

  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;


  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) 
    calloc( 1, sizeof(ProcessTable) );
  XLALGPSTimeNow(&(proctable.processTable->start_time));
  if (strcmp(CVS_REVISION,"$Revi" "sion$"))
    {
      XLALPopulateProcessTable(proctable.processTable,
          PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE, 0);
    }
  else
    {
      XLALPopulateProcessTable(proctable.processTable,
          PROGRAM_NAME, lalappsGitCommitID,
          lalappsGitGitStatus,
          lalappsGitCommitDate, 0);
    }
  snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );

  /* clear the waveform field */
  memset( waveform, 0, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR) );
  

  /*
   *
   * parse command line arguments
   *
   */

     
  while ( 1 )
  {
    /* LALgetopt_long stores long option here */
    int option_index = 0;
    long int gpsinput;
    size_t LALoptarg_len;

    c = LALgetopt_long_only( argc, argv,
        "a:A:b:B:d:f:hi:m:p:q:r:s:t:vZ:w:", long_options, &option_index );

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
        gpsinput = atol( LALoptarg );
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

      case 'b':
        gpsinput = atol( LALoptarg );
        if ( gpsinput < 441417609 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is prior to " 
              "Jan 01, 1994  00:00:00 UTC:\n"
              "(%ld specified)\n",
              long_options[option_index].name, gpsinput );
          exit( 1 );
        }
        gpsEndTime.gpsSeconds = gpsinput;
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "int", 
              "%ld", gpsinput );
        break;

      case 'f':
        fLower = atof( LALoptarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "float", 
              "%f", fLower );
        break;

      case 's':
        randSeed = atoi( LALoptarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "int", 
              "%d", randSeed );
        break;

      case 't':
        meanTimeStep = (REAL8) atof( LALoptarg );
        if ( meanTimeStep <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "time step must be > 0: (%e seconds specified)\n",
              long_options[option_index].name, meanTimeStep );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "float", 
              "%e", meanTimeStep );
        break;

      case 'i':
        timeInterval = atof( LALoptarg );
        if ( timeInterval < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "time interval must be >= 0: (%e seconds specified)\n",
              long_options[option_index].name, meanTimeStep );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", timeInterval );
        break;

      case 'A':
        /* minimum component mass */
        minMass = (REAL4) atof( LALoptarg );
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
        /* maximum component mass */
        maxMass = (REAL4) atof( LALoptarg );
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

    case 'x':
      /* maximum sum of components */
      sumMaxMass = (REAL4) atof( LALoptarg );
      if ( sumMaxMass <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "sum of two component masses must be > 0: "
                   "(%f solar masses specified)\n",
                   long_options[option_index].name, sumMaxMass );
          exit( 1 );
        }
      sumMaxMassUse=1;
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
                              "float", "%e", sumMaxMass );
        break;

      case 'p':
        /* minimum distance from earth */
        dmin = (REAL4) atof( LALoptarg );
        if ( dmin <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "minimum distance must be > 0: "
              "(%f kpc specified)\n",
              long_options[option_index].name, dmin );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", dmin );
        break;

      case 'r':
        /* max distance from earth */
        dmax = (REAL4) atof( LALoptarg );
        if ( dmax <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "maximum distance must be greater than 0: "
              "(%f kpc specified)\n",
              long_options[option_index].name, dmax );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", dmax );
        break;

      case 'd':
        ddistr = (UINT4) atoi( LALoptarg );
        if ( ddistr != 0 && ddistr != 1 && ddistr != 2)
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "DDISTR must be either 0 or 1 or 2\n",
              long_options[option_index].name);
          exit(1);
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "int", "%d", ddistr );

        break;

      case 'm':
        mdistr = (UINT4) atoi( LALoptarg );
        if ( mdistr != 0 && mdistr != 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "MDISTR must be either 0 or 1\n",
              long_options[option_index].name);
          exit(1);
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "int", "%d", mdistr );

  
        break;

      case 'Z':
        /* create storage for the usertag */
        LALoptarg_len = strlen( LALoptarg ) + 1;
        userTag = (CHAR *) calloc( LALoptarg_len, sizeof(CHAR) );
        memcpy( userTag, LALoptarg, LALoptarg_len );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "string", "%s", LALoptarg );
        break;

      case 'w':
        snprintf( waveform, LIGOMETA_WAVEFORM_MAX, "%s", LALoptarg);
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "%s", LALoptarg);
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "Binary Black Hole INJection generation routine\n" 
            "Duncan A Brown and Eirini Messaritaki\n"
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n" );
        fprintf( stdout, lalappsGitID );
        exit( 0 );
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

  if ( LALoptind < argc )
  {
    fprintf( stderr, "extraneous command line arguments:\n" );
    while ( LALoptind < argc )
    {
      fprintf ( stderr, "%s\n", argv[LALoptind++] );
    }
    exit( 1 );
  }


  if ( !*waveform )
  {
    /* use EOBtwoPN as the default waveform */
    snprintf( waveform, LIGOMETA_WAVEFORM_MAX, "EOBtwoPN");
  }

  if ( !fLower )
  {
    fprintf( stderr, "--f-lower must be specified and non-zero\n" );
    exit( 1 );
  }

  /*
   *
   * initialization
   *
   */


  /* initialize the random number generator */
  LAL_CALL( LALCreateRandomParams( &status, &randParams, randSeed ), &status );

  /* mass range, per component */
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


  /*
   *
   * loop over duration of desired output times
   *
   */


  while ( XLALGPSCmp( &gpsStartTime, &gpsEndTime ) < 0 )
  {

    /* rho, z and lGal are the galactocentric galactic axial coordinates */
    /* r and phi are the geocentric galactic spherical coordinates       */
#if 0
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    galacticPar.lGal = LAL_TWOPI * u;
    
    galacticPar.z = r * sinphi ;
    galacticPar.rho = 0.0 - Rcore * cos(galacticPar.lGal) +
      sqrt( r * r - galacticPar.z * galacticPar.z - 
      Rcore * Rcore * sin(galacticPar.lGal) * sin(galacticPar.lGal) );
#endif

#if 0
    if ( vrbflg ) fprintf( stdout, "%e %e %e %e %e\n", 
        galacticPar.m1, galacticPar.m2,
        galacticPar.rho * cos( galacticPar.lGal ),
        galacticPar.rho * sin( galacticPar.lGal ),
        galacticPar.z );
#endif

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
    /* XXX CHECK XXX */
    this_inj->geocent_end_time = gpsStartTime;
    if ( timeInterval )
    {
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      XLALGPSAdd( &(this_inj->geocent_end_time), u * timeInterval );
    }

    /* set gmst */
    this_inj->end_time_gmst = fmod(XLALGreenwichMeanSiderealTime(
        &this_inj->geocent_end_time), LAL_TWOPI) * 24.0 / LAL_TWOPI; /* hours */
    if( XLAL_IS_REAL8_FAIL_NAN(this_inj->end_time_gmst) )
    {
      fprintf(stderr, "XLALGreenwichMeanSiderealTime() failed\n");
      exit(1);
    }
    /* XXX END CHECK XXX */

    /* populate the sim_inspiral table */

    if (mdistr == 1)
    /* uniformly distributed mass1 and uniformly distributed mass2 */
    {
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      this_inj->mass1 = minMass + u * deltaM;
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      this_inj->mass2 = minMass + u * deltaM;
      mtotal = this_inj->mass1 + this_inj->mass2 ;
      this_inj->eta = this_inj->mass1 * this_inj->mass2 / ( mtotal * mtotal );
      this_inj->mchirp = (this_inj->mass1 + this_inj->mass2) * 
        pow(this_inj->eta, 0.6);
    }
    else if (mdistr == 0)
    /*uniformly distributed total mass */
    {
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status);
      mtotal = 2.0 * minMass + u * 2.0 *deltaM ;

      if (sumMaxMassUse==1) {
        while (mtotal > sumMaxMass) {
          LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status);
          mtotal = 2.0 * minMass + u * 2.0 *deltaM ;          
        }
      }

      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      this_inj->mass1 = minMass + u * deltaM;
      this_inj->mass2 = mtotal - this_inj->mass1;

      while (this_inj->mass1 >= mtotal || 
          this_inj->mass2 >= maxMass || this_inj->mass2 <= minMass )
      {
        LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
        this_inj->mass1 = minMass + u * deltaM;
        this_inj->mass2 = mtotal - this_inj->mass1;
      }
      this_inj->eta = this_inj->mass1 * this_inj->mass2 / ( mtotal * mtotal );
      this_inj->mchirp = (this_inj->mass1 + this_inj->mass2) * 
        pow(this_inj->eta, 0.6);

    }

     /* spatial distribution */

#if 0
     LAL_CALL( LALUniformDeviate( &status, &u, randParams ),
            &status );
     sinphi = 2.0 * u - 1.0;
     cosphi = sqrt( 1.0 - sinphi*sinphi );
#endif

     if (ddistr == 0)
     /* uniform distribution in distance */
     {
       REAL4 deltaD = dmax - dmin ;
       LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
       this_inj->distance = dmin + deltaD * u ;
      }
      else if (ddistr == 1)
      /* uniform distribution in log(distance) */
      {
        REAL4 lmin = log10(dmin);
        REAL4 lmax = log10(dmax);
        REAL4 deltaL = lmax - lmin;
        LAL_CALL(  LALUniformDeviate(&status,&u,randParams),&status );
        exponent = lmin + deltaL * u;
        this_inj->distance = pow(10.0,(REAL4) exponent);
      }
     else if (ddistr == 2)
     /* uniform volume distribution */
     {
       REAL4 d2min = pow(dmin,3.0) ;
       REAL4 d2max = pow(dmax,3.0) ;
       REAL4 deltad2 = d2max - d2min ;
       LAL_CALL(  LALUniformDeviate(&status,&u,randParams),&status );
       d2 = d2min + u * deltad2 ;
       this_inj->distance = pow(d2,1.0/3.0);
     }

     this_inj->distance = this_inj->distance / 1000.0; /*convert to Mpc */
       

      /* compute random longitude and latitude */
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      this_inj->latitude = asin( 2.0 * u - 1.0 ) ;
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      this_inj->longitude = LAL_TWOPI * u ;
     

#if 0
    LAL_CALL( LALGalacticInspiralParamsToSimInspiralTable( &status,
          this_inj, &galacticPar, randParams ), &status );
    if (vrbflg)
    { 
      fprintf( stdout, "%e\n",
      sqrt(galacticPar.z*galacticPar.z+galacticPar.rho*galacticPar.rho
          + Rcore*Rcore + 2.0*Rcore*galacticPar.rho*cos(galacticPar.lGal))-
      this_inj->distance*1000.0);
    }
#endif


    /* set the source and waveform fields */
    snprintf( this_inj->source, LIGOMETA_SOURCE_MAX, "???" );
    memcpy( this_inj->waveform, waveform, LIGOMETA_WAVEFORM_MAX *
        sizeof(CHAR));

    /* XXX CHECK XXX */
    /* compute random inclination, polarization and coalescence phase */
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->inclination = acos( 2.0 * u - 1.0 );
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->polarization = LAL_TWOPI * u ;
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->coa_phase = LAL_TWOPI * u ;
    
    /* populate the site specific information */
    LAL_CALL(LALPopulateSimInspiralSiteInfo( &status, this_inj ), 
        &status);
    
    /* increment the injection time */
    XLALGPSAdd( &gpsStartTime, meanTimeStep );

    /* finally populate the flower */
    if (fLower > 0)         
    {
        this_inj->f_lower = fLower;
    }
    else
    {
        this_inj->f_lower = 0;
    }
    /* XXX END CHECK XXX */

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
  XLALSimInspiralAssignIDs ( injections.simInspiralTable, 0, 0 );
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
    XLALFreeSimInspiral( &this_inj );
  }

  /* close the injection file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlfp ), &status );

  /* check for memory leaks and exit */
  LALCheckMemoryLeaks();
  return 0;
}
#endif
