/*
 * Copyright (C) 2006, 2008, 2009 Reinhard Prix
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

/**
 * \file
 * \ingroup pulsarApps
 * \author Reinhard Prix
 * \brief
 * This module deals with calculating various F-statistic metric approximations,
 * Fisher matrices and mismatches. Contrary to previous implementations
 * this consistently uses gsl-integration, and allows for very generic
 * handling of different coordinate systems in the Doppler parameter space.
 *
 * In particular, it allows for easy extensions to new coordinates and Doppler parameters.
 */

/* ---------- includes ---------- */
#include <math.h>
#include <errno.h>
#include <string.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


#include <lal/UserInput.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/LALInitBarycenter.h>
#include <lal/AVFactories.h>
#include <lal/SkyCoordinates.h>
#include <lal/ComputeFstat.h>
#include <lal/PulsarTimes.h>
#include <lal/SFTutils.h>

#include <lal/ComputeFstat.h>
#include <lal/LogPrintf.h>
#include <lal/StringVector.h>

#include <lal/FlatPulsarMetric.h>
#include <lal/UniversalDopplerMetric.h>

#include <lalapps.h>

/* ---------- Error codes and messages ---------- */
#define FSTATMETRIC_EMEM 	1
#define FSTATMETRIC_EINPUT	2
#define FSTATMETRIC_ENULL	3
#define FSTATMETRIC_ENONULL	4
#define FSTATMETRIC_EFILE	5
#define FSTATMETRIC_EXLAL	6


#define FSTATMETRIC_MSGEMEM	"Out of memory."
#define FSTATMETRIC_MSGEINPUT	"Invalid input."
#define FSTATMETRIC_MSGENULL	"Illegal NULL input"
#define FSTATMETRIC_MSGENONULL	"Output field not initialized to NULL"
#define FSTATMETRIC_MSGEFILE	"File IO error"
#define FSTATMETRIC_MSGEXLAL	"XLAL function call failed"

/* ----- compile switches ----- */
/* uncomment the following to turn off range-checking in GSL vector-functions */
/* #define GSL_RANGE_CHECK_OFF 1*/

#define METRIC_FORMAT	"% .16e"		/* fprintf-format used for printing metric components */

/*---------- local defines ---------- */
#define TRUE 		(1==1)
#define FALSE 		(1==0)

#define OneBillion 	1.0e9
#define ORB_V0		1e-4

/* ----- Macros ----- */
/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

/** Simple Euklidean scalar product for two 3-dim vectors in cartesian coords */
#define SCALAR(u,v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])

/** copy 3 components of Euklidean vector */
#define COPY_VECT(dst,src) do { (dst)[0] = (src)[0]; (dst)[1] = (src)[1]; (dst)[2] = (src)[2]; } while(0)

#define SQ(x) ((x) * (x))

/* ---------- local types ---------- */

/** Rudimentary first sketch of a history type, to encode all
 * the history-trail leading to a certain result from primal inputs.
 *
 * This will be extended in the future and moved into LAL.
 */
typedef struct
{
  CHAR *app_name;			/**< name (and path) of this application */
  CHAR *cmdline;			/**< commandline used to produce this result */
  CHAR *VCSInfoString;			/**< lal+lalapps source-version + configure string */
} ResultHistory_t;


typedef struct
{
  EphemerisData *edat;			/**< ephemeris data (from XLALInitBarycenter()) */
  LIGOTimeGPS startTime;		/**< start time of observation */
  PulsarParams signalParams;		/**< GW signal parameters: Amplitudes + doppler */
  MultiDetectorInfo detInfo;		/**< (multi-)detector info */
  DopplerCoordinateSystem coordSys; 	/**< array of enums describing Doppler-coordinates to compute metric in */
  ResultHistory_t *history;		/**< history trail leading up to and including this application */
} ConfigVariables;


typedef struct
{
  BOOLEAN help;

  LALStringVector* IFOs;	/**< list of detector-names "H1,H2,L1,.." or single detector*/
  LALStringVector* IFOweights; /**< list of relative detector-weights "w1, w2, w3, .." */

  REAL8 Freq;		/**< target-frequency */
  REAL8 Alpha;		/**< skyposition Alpha: radians, equatorial coords. */
  REAL8 Delta;		/**< skyposition Delta: radians, equatorial coords. */
  REAL8 f1dot;		/**< target 1. spindown-value df/dt */
  REAL8 f2dot;		/**< target 2. spindown-value d2f/dt2 */
  REAL8 f3dot;		/**< target 3. spindown-value d3f/dt3 */

  CHAR *ephemDir;	/**< directory to look for ephemeris files */
  CHAR *ephemYear;	/**< date-range string on ephemeris-files to use */

  REAL8 startTime;	/**< GPS start time of observation */
  REAL8 refTime;	/**< GPS reference time of Doppler parameters */
  REAL8 duration;	/**< length of observation in seconds */
  INT4 Nseg;		/**< number of segments to split duration into */

  REAL8 h0;		/**< GW amplitude h_0 */
  REAL8 cosi;		/**< cos(iota) */
  REAL8 psi;		/**< polarization-angle psi */
  REAL8 phi0;           /**< initial GW phase phi_0 */

  CHAR* outputMetric;	/**< filename to write metrics into */

  INT4 detMotionType;	/**< enum-value DetectorMotionType specifying type of detector-motion to use */

  INT4 metricType;	/**< type of metric to compute: 0=phase-metric, 1=F-metric(s), 2=both */

  INT4 projection;     /**< project metric onto subspace orthogonal to this coordinate-axis (0=none, 1=1st-coordinate axis, ...) */

  LALStringVector* coords; /**< list of Doppler-coordinates to compute metric in, see --coordsHelp for possible values */
  BOOLEAN coordsHelp;	/**< output help-string explaining all the possible Doppler-coordinate names for --cords */

  BOOLEAN approxPhase;	/**< use an approximate phase-model, neglecting Roemer delay in spindown coordinates */

  BOOLEAN version;	/**< output code versions */

} UserVariables_t;


/*---------- empty structs for initializations ----------*/
ConfigVariables empty_ConfigVariables;
PulsarTimesParamStruc empty_PulsarTimesParamStruc;
UserVariables_t empty_UserVariables;
/* ---------- global variables ----------*/
extern int vrbflg;

/* ---------- local prototypes ---------- */
int initUserVars (UserVariables_t *uvar);
int XLALInitCode ( ConfigVariables *cfg, const UserVariables_t *uvar, const char *app_name);

EphemerisData *InitEphemeris (const CHAR *ephemDir, const CHAR *ephemYear );

int XLALOutputDopplerMetric ( FILE *fp, const DopplerMetric *metric, const ResultHistory_t *history );

int XLALDestroyConfig ( ConfigVariables *cfg );
void XLALDestroyResultHistory ( ResultHistory_t * history );

/*============================================================
 * FUNCTION definitions
 *============================================================*/

int
main(int argc, char *argv[])
{
  ConfigVariables config = empty_ConfigVariables;
  UserVariables_t uvar = empty_UserVariables;
  DopplerMetricParams metricParams = empty_DopplerMetricParams;

  lalDebugLevel = 0;
  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register user-variables */
  if ( XLALGetDebugLevel (argc, argv, 'v') != XLAL_SUCCESS ) {
    XLALPrintError( "%s(): XLALGetDebugLevel() failed\n", __func__ );
    return EXIT_FAILURE;
  }

  /* set log-level */
  LogSetLevel ( lalDebugLevel );

  if ( initUserVars(&uvar) != XLAL_SUCCESS ) {
    XLALPrintError( "%s(): initUserVars() failed\n", __func__ );
    return EXIT_FAILURE;
  }

  /* read cmdline & cfgfile  */
  if ( XLALUserVarReadAllInput(argc,argv) != XLAL_SUCCESS ) {
    XLALPrintError( "%s(): XLALUserVarReadAllInput() failed\n", __func__ );
    return EXIT_FAILURE;
  }

  if (uvar.help) 	/* help requested: we're done */
    return 0;

  CHAR *VCSInfoString;
  if ( (VCSInfoString = XLALGetVersionString(0)) == NULL ) {
    XLALPrintError("XLALGetVersionString(0) failed.\n");
    exit(1);
  }

  if ( uvar.version ) {
    printf ( "%s\n", VCSInfoString );
    return 0;
  }


  if ( uvar.coordsHelp )
    {
      CHAR *helpstr;
      if ( (helpstr = XLALDopplerCoordinateHelpAll()) == NULL )
	{
	  LogPrintf ( LOG_CRITICAL, "XLALDopplerCoordinateHelpAll() failed!\n\n");
	  return -1;
	}
      printf ( "\n%s\n", helpstr );
      XLALFree ( helpstr );
      return 0;
    } /* if coordsHelp */

  /* basic setup and initializations */
  XLAL_CHECK ( XLALInitCode( &config, &uvar, argv[0] ) == XLAL_SUCCESS, XLAL_EFUNC, "XLALInitCode() failed with xlalErrno = %d\n\n", xlalErrno );
  config.history->VCSInfoString = VCSInfoString;


  metricParams.coordSys      = config.coordSys;
  metricParams.detMotionType = uvar.detMotionType;
  metricParams.metricType    = uvar.metricType;
  metricParams.detInfo       = config.detInfo;
  metricParams.signalParams  = config.signalParams;
  metricParams.projectCoord  = uvar.projection - 1;	/* user-input counts from 1, but interally we count 0=1st coord. (-1==no projection) */
  metricParams.approxPhase   = uvar.approxPhase;

  XLAL_CHECK ( XLALSegListInitSimpleSegments ( &metricParams.segmentList, config.startTime, uvar.Nseg, uvar.duration / uvar.Nseg ) == XLAL_SUCCESS,
               XLAL_EFUNC, "XLALSegListInitSimpleSegments() failed with xlalErrno = %d\n", xlalErrno );

  /* ----- compute metric full metric + Fisher matrix ---------- */
  DopplerMetric *metric;
  if ( (metric = XLALDopplerFstatMetric ( &metricParams, config.edat )) == NULL ) {
    LogPrintf (LOG_CRITICAL, "Something failed in XLALDopplerFstatMetric(). xlalErrno = %d\n\n", xlalErrno);
    return -1;
  }

  /* ---------- output results ---------- */
  if ( uvar.outputMetric )
    {
      FILE *fpMetric;
      if ( (fpMetric = fopen ( uvar.outputMetric, "wb" )) == NULL ) {
	LogPrintf (LOG_CRITICAL, "%s: failed to open '%s' for writing. error = '%s'\n",
		   __func__, uvar.outputMetric, strerror(errno));
	return FSTATMETRIC_EFILE;
      }

      if ( XLALOutputDopplerMetric ( fpMetric, metric, config.history ) != XLAL_SUCCESS  ) {
	LogPrintf (LOG_CRITICAL, "%s: failed to write Doppler metric into output-file '%s'. xlalErrno = %d\n\n",
		   __func__, uvar.outputMetric, xlalErrno );
	return FSTATMETRIC_EFILE;
      }

      fclose ( fpMetric );

    } /* if outputMetric */

  /* ----- done: free all memory */
  XLALDestroyDopplerMetric ( metric );
  if ( XLALDestroyConfig( &config ) != XLAL_SUCCESS ) {
    LogPrintf (LOG_CRITICAL, "%s: XLADestroyConfig() failed, xlalErrno = %d.\n\n", __func__, xlalErrno );
    return FSTATMETRIC_EXLAL;
  }

  LALCheckMemoryLeaks();

  return 0;
} /* main */


/** register all "user-variables" */
int
initUserVars (UserVariables_t *uvar)
{

  /* set a few defaults */
  uvar->help = FALSE;

#define EPHEM_YEAR  "00-04"
  uvar->ephemYear = XLALCalloc (1, strlen(EPHEM_YEAR)+1);
  strcpy (uvar->ephemYear, EPHEM_YEAR);

  uvar->Freq = 100;
  uvar->f1dot = 0.0;
  uvar->f2dot = 0.0;
  uvar->f3dot = 0.0;
  uvar->h0 = 1;
  uvar->phi0 = 0;

  uvar->startTime = 714180733;
  uvar->duration = 10 * 3600;
  uvar->Nseg = 1;
  uvar->refTime = -1;	/* default: use mid-time */

  uvar->projection = 0;
  if ( (uvar->IFOs = XLALCreateStringVector ( "H1", NULL )) == NULL ) {
    LogPrintf (LOG_CRITICAL, "Call to XLALCreateStringVector() failed with xlalErrno = %d\n", xlalErrno );
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  uvar->IFOweights = NULL;

  uvar->detMotionType = DETMOTION_SPIN_ORBIT;
  uvar->metricType = 0;	/* by default: compute only phase metric */

  if ( (uvar->coords = XLALCreateStringVector ( "freq", "alpha", "delta", "f1dot", NULL )) == NULL ) {
    LogPrintf (LOG_CRITICAL, "Call to XLALCreateStringVector() failed with xlalErrno = %d\n", xlalErrno );
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  uvar->approxPhase = FALSE;

  /* register all our user-variables */

  XLALregBOOLUserStruct(help,		'h', UVAR_HELP,		"Print this help/usage message");
  XLALregLISTUserStruct(IFOs,		'I', UVAR_OPTIONAL, 	"Comma-separated list of detectors, eg. \"H1,H2,L1,G1, ...\" ");
  XLALregLISTUserStruct(IFOweights,	 0,  UVAR_OPTIONAL, 	"Comma-separated list of relative noise-weights, eg. \"w1,w2,w3,..\" ");
  XLALregREALUserStruct(Alpha,		'a', UVAR_OPTIONAL,	"skyposition Alpha in radians, equatorial coords.");
  XLALregREALUserStruct(Delta, 		'd', UVAR_OPTIONAL,	"skyposition Delta in radians, equatorial coords.");
  XLALregREALUserStruct(Freq, 		'f', UVAR_OPTIONAL, 	"target frequency");
  XLALregREALUserStruct(f1dot, 		's', UVAR_OPTIONAL, 	"first spindown-value df/dt");
  XLALregREALUserStruct(f2dot, 		 0 , UVAR_OPTIONAL, 	"second spindown-value d2f/dt2");
  XLALregREALUserStruct(f3dot, 		 0 , UVAR_OPTIONAL, 	"third spindown-value d3f/dt3");
  XLALregREALUserStruct(startTime,      't', UVAR_OPTIONAL, 	"GPS start time of observation");
  XLALregREALUserStruct(refTime,         0,  UVAR_OPTIONAL, 	"GPS reference time of Doppler parameters. Special values: 0=startTime, -1=mid-time");
  XLALregREALUserStruct(duration,	'T', UVAR_OPTIONAL,	"Duration of observation in seconds");
  XLALregINTUserStruct(Nseg,		'N', UVAR_OPTIONAL, 	"Compute semi-coherent metric for this number of segments within 'duration'" );
  XLALregSTRINGUserStruct(ephemDir, 	'E', UVAR_OPTIONAL,     "Directory where Ephemeris files are located");
  XLALregSTRINGUserStruct(ephemYear, 	'y', UVAR_OPTIONAL,     "Year (or range of years) of ephemeris files to be used");

  XLALregREALUserStruct(h0,	 	 0, UVAR_OPTIONAL,	"GW amplitude h0" );
  XLALregREALUserStruct(cosi,	 	 0, UVAR_OPTIONAL,	"Pulsar orientation-angle cos(iota) [-1,1]" );
  XLALregREALUserStruct(psi,		 0, UVAR_OPTIONAL,	"Wave polarization-angle psi [0, pi]" );
  XLALregREALUserStruct(phi0,		 0, UVAR_OPTIONAL,	"GW initial phase phi_0 [0, 2pi]" );

  XLALregINTUserStruct(metricType,	 0,  UVAR_OPTIONAL,	"type of metric to compute: 0=phase-metric, 1=F-metric(s), 2=both" );
  XLALregSTRINGUserStruct(outputMetric,	'o', UVAR_OPTIONAL,	"Output the metric components (in octave format) into this file.");
  XLALregINTUserStruct(projection,      0,  UVAR_OPTIONAL,     "Project onto subspace orthogonal to this axis: 0=none, 1=1st-coord, 2=2nd-coord etc");

  XLALregLISTUserStruct(coords,		'c', UVAR_OPTIONAL, 	"Doppler-coordinates to compute metric in (see --coordsHelp)");
  XLALregBOOLUserStruct(coordsHelp,      0,  UVAR_OPTIONAL,     "output help-string explaining all the possible Doppler-coordinate names for --coords");

  XLALregINTUserStruct(detMotionType,	 0,  UVAR_DEVELOPER,	"Detector-motion: 0=spin+orbit, 1=orbit, 2=spin, 3=spin+ptoleorbit, 4=ptoleorbit, 5=orbit+spin_z, 6=orbit+spin_xy");
  XLALregBOOLUserStruct(approxPhase,     0,  UVAR_DEVELOPER,	"Use an approximate phase-model, neglecting Roemer delay in spindown coordinates (or orders >= 1)");

  XLALregBOOLUserStruct(version,        'V', UVAR_SPECIAL,      "Output code version");

  return XLAL_SUCCESS;

} /* initUserVars() */


/** basic initializations: set-up 'ConfigVariables'
 */
int
XLALInitCode ( ConfigVariables *cfg, const UserVariables_t *uvar, const char *app_name)
{
  if ( !cfg || !uvar || !app_name ) {
    LogPrintf (LOG_CRITICAL, "%s: illegal NULL pointer input.\n\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL );
  }

  /* ----- check user-input consistency ---------- */
  XLAL_CHECK ( uvar->Nseg >= 1, XLAL_EDOM, "Invalid input --Nseg=%d: number of segments must be >= 1\n", uvar->Nseg );
  XLAL_CHECK ( uvar->duration >= 1, XLAL_EDOM, "Invalid input --duration=%f: duration must be >= 1 s\n", uvar->duration );

  /* ----- determine start-time from user-input */
  XLALGPSSetREAL8( &(cfg->startTime), uvar->startTime );

  if ( (cfg->edat = InitEphemeris ( uvar->ephemDir, uvar->ephemYear)) == NULL ) {
    LogPrintf (LOG_CRITICAL, "%s: InitEphemeris() Failed to initialize ephemeris data!\n\n", __func__);
    XLAL_ERROR ( XLAL_EFUNC );
  }

  /* ----- figure out reference time */
  REAL8 refTime;
  LIGOTimeGPS refTimeGPS;
  /* treat special values first */
  if ( uvar->refTime == 0 )		/* 0 = use startTime */
    refTime = uvar->startTime;
  else if ( uvar->refTime == -1 )	/* -1 = use mid-time of observation */
    refTime = uvar->startTime + 0.5 * uvar->duration;
  else
    refTime = uvar->refTime;

  XLALGPSSetREAL8( &refTimeGPS, refTime );

  /* ----- get parameter-space point from user-input) */
  cfg->signalParams.Amp.h0 = uvar->h0;
  cfg->signalParams.Amp.cosi = uvar->cosi;
  cfg->signalParams.Amp.psi = uvar->psi;
  cfg->signalParams.Amp.phi0 = uvar->phi0;

  {
    PulsarDopplerParams *dop = &(cfg->signalParams.Doppler);
    (*dop) = empty_PulsarDopplerParams;
    dop->refTime = refTimeGPS;
    dop->Alpha = uvar->Alpha;
    dop->Delta = uvar->Delta;
    dop->fkdot[0] = uvar->Freq;
    dop->fkdot[1] = uvar->f1dot;
    dop->fkdot[2] = uvar->f2dot;
    dop->fkdot[3] = uvar->f3dot;
    dop->orbit = NULL;
  }

  /* ----- initialize IFOs and (Multi-)DetectorStateSeries  ----- */
  {
    UINT4 numDet = uvar->IFOs->length;

    if ( numDet > DOPPLERMETRIC_MAX_DETECTORS ) {
      LogPrintf (LOG_CRITICAL, "%s: More detectors (%d) specified than allowed (%d)\n", __func__, numDet, DOPPLERMETRIC_MAX_DETECTORS );
      XLAL_ERROR ( XLAL_EINVAL );
    }
    if ( uvar->IFOweights && (uvar->IFOweights->length != numDet ) )
      {
	LogPrintf (LOG_CRITICAL, "%s: number of IFOweights (%d) must agree with the number of IFOs (%d)!\n\n",
		   __func__, uvar->IFOweights->length, numDet );
	XLAL_ERROR ( XLAL_EINVAL );
      }

    cfg->detInfo = empty_MultiDetectorInfo;
    if ( XLALParseMultiDetectorInfo ( &cfg->detInfo, uvar->IFOs, uvar->IFOweights ) != XLAL_SUCCESS ) {
      LogPrintf (LOG_CRITICAL, "%s: XLALParseMultiDetectorInfo() failed to parse detector names and/or weights. errno = %d.\n\n", __func__, xlalErrno);
      XLAL_ERROR ( XLAL_EFUNC );
    }

  } /* handle detector input */


  /* ---------- translate coordinate system into internal representation ---------- */
  if ( XLALDopplerCoordinateNames2System ( &cfg->coordSys, uvar->coords ) ) {
    LogPrintf (LOG_CRITICAL, "%s: Call to XLALDopplerCoordinateNames2System() failed. errno = %d\n\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  /* ---------- record full 'history' up to and including this application ---------- */
  {
    CHAR *cmdline = NULL;
    CHAR *tmp;
    size_t len = strlen ( app_name ) + 1;

    if ( (cfg->history = XLALCalloc ( 1, sizeof(*cfg->history))) == NULL ) {
      LogPrintf (LOG_CRITICAL, "%s: XLALCalloc(1,%d) failed.\n\n", __func__, sizeof(*cfg->history));
      XLAL_ERROR ( XLAL_ENOMEM );
    }

    if ( (tmp = XLALMalloc ( len )) == NULL ) {
      LogPrintf (LOG_CRITICAL, "%s: XLALMalloc (%s) failed.\n\n", __func__, len );
      XLAL_ERROR ( XLAL_ENOMEM );
    }
    strcpy ( tmp, app_name );
    cfg->history->app_name = tmp;

    /* get commandline describing search*/
    if ( (cmdline = XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE )) == NULL ) {
      LogPrintf (LOG_CRITICAL, "%s: XLALUserVarGetLog() failed with xlalErrno = %d.\n\n", __func__, xlalErrno );
      XLAL_ERROR ( XLAL_EFUNC );
    }
    cfg->history->cmdline = cmdline;
  } /* record history */


  return XLAL_SUCCESS;

} /* XLALInitCode() */


/** Destructor for internal configuration struct */
int
XLALDestroyConfig ( ConfigVariables *cfg )
{
  if ( !cfg ) {
    LogPrintf (LOG_CRITICAL, "%s: invalid NULL input!\n\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL );
  }

  XLALDestroyUserVars ();

  XLALDestroyResultHistory ( cfg->history );

  XLALDestroyEphemerisData ( cfg->edat );

  return XLAL_SUCCESS;

} /* XLALDestroyConfig() */



/** Load Ephemeris from ephemeris data-files  */
EphemerisData *
InitEphemeris (const CHAR *ephemDir,	/**< directory containing ephems */
	       const CHAR *ephemYear	/**< which years do we need? */
	       )
{
#define FNAME_LENGTH 1024
  EphemerisData *edat;
  CHAR EphemEarth[FNAME_LENGTH];	/* filename of earth-ephemeris data */
  CHAR EphemSun[FNAME_LENGTH];	/* filename of sun-ephemeris data */

  if ( !ephemYear ) {
    XLALPrintError ("\n%s: NULL pointer passed as ephemeris year range!\n", __func__);
    return NULL;
  }

  if ( ephemDir )
    {
      snprintf(EphemEarth, FNAME_LENGTH, "%s/earth%s.dat", ephemDir, ephemYear);
      snprintf(EphemSun, FNAME_LENGTH, "%s/sun%s.dat", ephemDir, ephemYear);
    }
  else
    {
      snprintf(EphemEarth, FNAME_LENGTH, "earth%s.dat", ephemYear);
      snprintf(EphemSun, FNAME_LENGTH, "sun%s.dat",  ephemYear);
    }

  EphemEarth[FNAME_LENGTH-1] = 0;
  EphemSun[FNAME_LENGTH-1] = 0;

  if ( (edat = XLALInitBarycenter(EphemEarth, EphemSun)) == NULL ) {
    XLALPrintError ( "%s: XLALInitBarycenter() failed with xlalErrno = %d.\n\n", __func__, xlalErrno );
    return NULL;
  }

  return edat;

} /* InitEphemeris() */

int
XLALOutputDopplerMetric ( FILE *fp, const DopplerMetric *metric, const ResultHistory_t *history )
{
  UINT4 i;
  REAL8 A, B, C, D;
  const DopplerMetricParams *meta;
  const PulsarDopplerParams *doppler;
  const PulsarAmplitudeParams *Amp;

  // ----- input sanity checks
  if ( !fp || !metric ) {
    LogPrintf (LOG_CRITICAL, "%s: illegal NULL input.\n\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  XLAL_CHECK ( XLALSegListIsInitialized ( &(metric->meta.segmentList) ), XLAL_EINVAL, "Got un-initialized segment list in 'metric->meta.segmentList'\n" );
  UINT4 Nseg = metric->meta.segmentList.length;
  XLAL_CHECK ( Nseg >= 1, XLAL_EDOM, "Got invalid zero-length segment list 'metric->meta.segmentList'\n" );

  /* useful shortcuts */
  meta = &(metric->meta);
  doppler = &(meta->signalParams.Doppler);
  Amp = &(meta->signalParams.Amp);

  /* output history info */
  if ( history )
    {
      if ( history->app_name ) fprintf (fp, "%%%% app_name: %s\n", history->app_name );
      if ( history->cmdline) fprintf (fp, "%%%% commandline: %s\n", history->cmdline );
      if ( history->VCSInfoString ) fprintf (fp, "%%%% Code Version: %s\n", history->VCSInfoString );
    }

  fprintf ( fp, "%%%% DopplerCoordinates = [ " );
  for ( i=0; i < meta->coordSys.dim; i ++ )
    {
      if ( i > 0 ) fprintf ( fp, ", " );
      fprintf ( fp, "%s", XLALDopplerCoordinateName(meta->coordSys.coordIDs[i]));
    }
  fprintf ( fp, "];\n");

  { /* output projection info */
    const char *pname;
    if ( meta->projectCoord < 0 )
      pname = "None";
    else
      pname = XLALDopplerCoordinateName ( meta->coordSys.coordIDs[meta->projectCoord] );

    fprintf ( fp, "%%%% Projection onto subspace orthogonal to coordinate: '%s'\n", pname);
  }

  fprintf ( fp, "%%%% DetectorMotionType = '%s'\n", XLALDetectorMotionName(meta->detMotionType) );
  fprintf ( fp, "%%%% h0 = %g; cosi = %g; psi = %g; phi0 = %g;\n", Amp->h0, Amp->cosi, Amp->psi, Amp->phi0 );
  fprintf ( fp, "%%%% DopplerPoint = {\n");
  fprintf ( fp, "%%%% 	refTime = {%d, %d}\n",
	    doppler->refTime.gpsSeconds, doppler->refTime.gpsNanoSeconds );
  fprintf ( fp, "%%%% 	Alpha = %f rad; Delta = %f rad\n", doppler->Alpha, doppler->Delta );
  fprintf ( fp, "%%%% 	fkdot = [%f, %g, %g, %g ]\n",
	    doppler->fkdot[0], doppler->fkdot[1], doppler->fkdot[2], doppler->fkdot[3] );
  if ( doppler->orbit )
    {
      const BinaryOrbitParams *orbit = doppler->orbit;
      fprintf ( fp, "%%%% 	   orbit = { \n");
      fprintf ( fp, "%%%% 		tp = {%d, %d}\n", orbit->tp.gpsSeconds, orbit->tp.gpsNanoSeconds );
      fprintf ( fp, "%%%% 		argp  = %g\n", orbit->argp );
      fprintf ( fp, "%%%% 		asini = %g\n", orbit->asini );
      fprintf ( fp, "%%%% 		ecc = %g\n", orbit->ecc );
      fprintf ( fp, "%%%% 		period = %g\n", orbit->period );
      fprintf ( fp, "%%%% 	   }\n");
    } /* if doppler->orbit */
  fprintf ( fp, "%%%% }\n");

  LIGOTimeGPS *tStart = &(meta->segmentList.segs[0].start);
  LIGOTimeGPS *tEnd   = &(meta->segmentList.segs[Nseg-1].end);
  REAL8 Tspan = XLALGPSDiff ( tEnd, tStart );
  fprintf ( fp, "%%%% startTime = {%d, %d}\n", tStart->gpsSeconds, tStart->gpsNanoSeconds );
  fprintf ( fp, "%%%% Tspan     = %f\n", Tspan );
  fprintf ( fp, "%%%% Nseg      = %d\n", Nseg );
  // FIXME: Output full segment list here:

  fprintf ( fp, "%%%% detectors = [");
  for ( i=0; i < meta->detInfo.length; i ++ )
    {
      if ( i > 0 ) fprintf ( fp, ", ");
      fprintf ( fp, "%s", meta->detInfo.sites[i].frDetector.name );
    }
  fprintf ( fp, "];\n");
  fprintf ( fp, "%%%% detectorWeights = [");
  for ( i=0; i < meta->detInfo.length; i ++ )
    {
      if ( i > 0 ) fprintf ( fp, ", ");
      fprintf ( fp, "%f", meta->detInfo.detWeights[i] );
    }
  fprintf ( fp, "];\n");

  /* ----- output phase metric ---------- */
  if ( metric->g_ij )
    {
      fprintf ( fp, "\ng_ij = \\\n" ); XLALfprintfGSLmatrix ( fp, METRIC_FORMAT,  metric->g_ij );
      fprintf ( fp, "maxrelerr_gPh = %.2e;\n", metric->maxrelerr_gPh );

      gsl_matrix *gN_ij;
      if ( (gN_ij = XLALNaturalizeMetric ( metric->g_ij, meta )) == NULL ) {
        XLALPrintError ("%s: something failed Naturalizing phase metric g_ij!\n", __func__ );
        XLAL_ERROR ( XLAL_EFUNC );
      }
      fprintf ( fp, "\ngN_ij = \\\n" ); XLALfprintfGSLmatrix ( fp, METRIC_FORMAT,  gN_ij );
      gsl_matrix_free ( gN_ij );

      gsl_matrix *gDN_ij;
      if ( (gDN_ij = XLALDiagNormalizeMetric ( metric->g_ij )) == NULL ) {
        XLALPrintError ("%s: something failed NormDiagonalizing phase metric g_ij!\n", __func__ );
        XLAL_ERROR ( XLAL_EFUNC );
      }
      fprintf ( fp, "\ngDN_ij = \\\n" ); XLALfprintfGSLmatrix ( fp, METRIC_FORMAT,  gDN_ij );
      gsl_matrix_free ( gDN_ij );
    }

  /* ----- output F-metric (and related matrices ---------- */
  if ( metric->gF_ij )
    {
      fprintf ( fp, "\ngF_ij = \\\n" );   XLALfprintfGSLmatrix ( fp, METRIC_FORMAT,  metric->gF_ij );
      fprintf ( fp, "\ngFav_ij = \\\n" ); XLALfprintfGSLmatrix ( fp, METRIC_FORMAT,  metric->gFav_ij );
      fprintf ( fp, "\nm1_ij = \\\n" );   XLALfprintfGSLmatrix ( fp, METRIC_FORMAT,  metric->m1_ij );
      fprintf ( fp, "\nm2_ij = \\\n" );   XLALfprintfGSLmatrix ( fp, METRIC_FORMAT,  metric->m2_ij );
      fprintf ( fp, "\nm3_ij = \\\n" );   XLALfprintfGSLmatrix ( fp, METRIC_FORMAT,  metric->m3_ij );
      fprintf ( fp, "maxrelerr_gF = %.2e;\n", metric->maxrelerr_gF );
    }

  /*  ----- output Fisher matrix ---------- */
  if ( metric->Fisher_ab )
    {
      A = gsl_matrix_get ( metric->Fisher_ab, 0, 0 );
      B = gsl_matrix_get ( metric->Fisher_ab, 1, 1 );
      C = gsl_matrix_get ( metric->Fisher_ab, 0, 1 );

      D = A * B - C * C;

      fprintf ( fp, "\nA = %.16g;\nB = %.16g;\nC = %.16g;\nD = %.16g;\n", A, B, C, D );
      fprintf ( fp, "\nrho2 = %.16g;\n", metric->rho2 );

      fprintf (fp, "\nFisher_ab = \\\n" ); XLALfprintfGSLmatrix ( fp, METRIC_FORMAT,  metric->Fisher_ab );
    }

  return XLAL_SUCCESS;

} /* XLALOutputDopplerMetric() */


/** Destructor for ResultHistory type
 */
void
XLALDestroyResultHistory ( ResultHistory_t * history )
{

  if ( !history )
    return;

  if ( history->app_name )
    XLALFree ( history->app_name );

  if ( history->cmdline )
    XLALFree ( history->cmdline );

  if ( history->VCSInfoString )
    XLALFree ( history->VCSInfoString );

  XLALFree ( history );

  return;

} /* XLALDestroyResultHistory() */
