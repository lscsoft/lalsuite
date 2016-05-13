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
 * \ingroup lalapps_pulsar_Fstatistic
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
#include <lal/LALString.h>

#include <lal/ComputeFstat.h>
#include <lal/LogPrintf.h>
#include <lal/StringVector.h>

#include <lal/FlatPulsarMetric.h>
#include <lal/UniversalDopplerMetric.h>
#include <lal/MetricUtils.h>

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

/**
 * Rudimentary first sketch of a history type, to encode all
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
  LALSegList segmentList;		/**< segment list contains start- and end-times of all segments */
  PulsarParams signalParams;		/**< GW signal parameters: Amplitudes + doppler */
  MultiLALDetector multiIFO;		/**< (multi-)detector info */
  MultiNoiseFloor multiNoiseFloor;	/**< ... and corresponding noise-floors to be used as weights */
  DopplerCoordinateSystem coordSys; 	/**< array of enums describing Doppler-coordinates to compute metric in */
  ResultHistory_t *history;		/**< history trail leading up to and including this application */
} ConfigVariables;


typedef struct
{
  LALStringVector *IFOs;	/**< list of detector-names "H1,H2,L1,.." or single detector*/
  LALStringVector *sqrtSX; 	/**< (string-) list of per-IFO sqrt{Sn} values, \f$\sqrt{S_X}\f$ */

  REAL8 Freq;		/**< target-frequency */
  REAL8 Alpha;		/**< skyposition Alpha: radians, equatorial coords. */
  REAL8 Delta;		/**< skyposition Delta: radians, equatorial coords. */
  REAL8 f1dot;		/**< target 1. spindown-value df/dt */
  REAL8 f2dot;		/**< target 2. spindown-value d2f/dt2 */
  REAL8 f3dot;		/**< target 3. spindown-value d3f/dt3 */
  REAL8 orbitasini;	/**< target projected semimajor axis of binary orbit [Units: light seconds] */
  REAL8 orbitPeriod;	/**< target period of binary orbit [Units: s]. */
  LIGOTimeGPS orbitTp;	/**< target time of periapse passage of the CW source in a binary orbit [Units: GPS seconds] */
  REAL8 orbitArgp;	/**< target argument of periapse of binary orbit [Units: rad] */
  REAL8 orbitEcc;	/**< target eccentricity of binary orbit [Units: none] */

  CHAR *ephemEarth;	/**< Earth ephemeris file to use */
  CHAR *ephemSun;	/**< Sun ephemeris file to use */

  LIGOTimeGPS refTime;	/**< GPS reference time of Doppler parameters */

  LIGOTimeGPS startTime;	/**< GPS start time of observation */
  REAL8 duration;	/**< length of observation in seconds */
  INT4 Nseg;		/**< number of segments to split duration into */
  CHAR *segmentList;	/**< ALTERNATIVE: specify segment file with format: repeated lines <startGPS endGPS duration[h] NumSFTs>" */

  REAL8 h0;		/**< GW amplitude h_0 */
  REAL8 cosi;		/**< cos(iota) */
  REAL8 psi;		/**< polarization-angle psi */
  REAL8 phi0;           /**< initial GW phase phi_0 */

  CHAR* outputMetric;	/**< filename to write metrics into */

  CHAR *detMotionStr;	/**< string specifying type of detector-motion to use */

  INT4 metricType;	/**< type of metric to compute: 0=phase-metric, 1=F-metric(s), 2=both */

  INT4 projection;     /**< project metric onto subspace orthogonal to this coordinate-axis (0=none, 1=1st-coordinate axis, ...) */

  LALStringVector* coords; /**< list of Doppler-coordinates to compute metric in, see --coordsHelp for possible values */
  BOOLEAN coordsHelp;	/**< output help-string explaining all the possible Doppler-coordinate names for --cords */

  BOOLEAN approxPhase;	/**< use an approximate phase-model, neglecting Roemer delay in spindown coordinates */

  BOOLEAN version;	/**< output code versions */

} UserVariables_t;

/* ---------- global variables ----------*/
extern int vrbflg;

/* ---------- local prototypes ---------- */
int initUserVars (UserVariables_t *uvar);
int XLALInitCode ( ConfigVariables *cfg, const UserVariables_t *uvar, const char *app_name);

int XLALOutputDopplerMetric ( FILE *fp, const DopplerPhaseMetric *Pmetric, const DopplerFstatMetric *Fmetric, const ResultHistory_t *history );

int XLALDestroyConfig ( ConfigVariables *cfg );
void XLALDestroyResultHistory ( ResultHistory_t * history );

/*============================================================
 * FUNCTION definitions
 *============================================================*/

int
main(int argc, char *argv[])
{
  ConfigVariables XLAL_INIT_DECL(config);
  UserVariables_t XLAL_INIT_DECL(uvar);
  DopplerMetricParams XLAL_INIT_DECL(metricParams);

  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register user-variables */
  if ( initUserVars(&uvar) != XLAL_SUCCESS ) {
    XLALPrintError( "%s(): initUserVars() failed\n", __func__ );
    return EXIT_FAILURE;
  }

  /* read cmdline & cfgfile  */
  BOOLEAN should_exit = 0;
  XLAL_CHECK( XLALUserVarReadAllInput( &should_exit, argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit )
    return EXIT_FAILURE;

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

  /* parse detector motion string */
  int detMotionType = XLALParseDetectorMotionString( uvar.detMotionStr );
  XLAL_CHECK ( detMotionType != XLAL_FAILURE, XLAL_EFUNC, "Failed to pass detector motion string '%s'", uvar.detMotionStr );
  metricParams.detMotionType = detMotionType;

  metricParams.segmentList   = config.segmentList;
  metricParams.coordSys      = config.coordSys;
  metricParams.multiIFO      = config.multiIFO;
  metricParams.multiNoiseFloor = config.multiNoiseFloor;
  metricParams.signalParams  = config.signalParams;
  metricParams.projectCoord  = uvar.projection - 1;	/* user-input counts from 1, but interally we count 0=1st coord. (-1==no projection) */
  metricParams.approxPhase   = uvar.approxPhase;


  /* ----- compute metric full metric + Fisher matrix ---------- */
  DopplerPhaseMetric *Pmetric = NULL;
  if ( uvar.metricType == 0 || uvar.metricType == 2 ) {
    if ( (Pmetric = XLALComputeDopplerPhaseMetric ( &metricParams, config.edat )) == NULL ) {
      LogPrintf (LOG_CRITICAL, "Something failed in XLALComputeDopplerPhaseMetric(). xlalErrno = %d\n\n", xlalErrno);
      return -1;
    }
  }
  DopplerFstatMetric *Fmetric = NULL;
  if ( uvar.metricType == 1 || uvar.metricType == 2 ) {
    if ( (Fmetric = XLALComputeDopplerFstatMetric ( &metricParams, config.edat )) == NULL ) {
      LogPrintf (LOG_CRITICAL, "Something failed in XLALComputeDopplerFstatMetric(). xlalErrno = %d\n\n", xlalErrno);
      return -1;
    }
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

      if ( XLALOutputDopplerMetric ( fpMetric, Pmetric, Fmetric, config.history ) != XLAL_SUCCESS  ) {
	LogPrintf (LOG_CRITICAL, "%s: failed to write Doppler metric into output-file '%s'. xlalErrno = %d\n\n",
		   __func__, uvar.outputMetric, xlalErrno );
	return FSTATMETRIC_EFILE;
      }

      fclose ( fpMetric );

    } /* if outputMetric */

  /* ----- done: free all memory */
  XLALDestroyDopplerPhaseMetric ( Pmetric );
  XLALDestroyDopplerFstatMetric ( Fmetric );
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
  uvar->ephemEarth = XLALStringDuplicate("earth00-19-DE405.dat.gz");
  uvar->ephemSun = XLALStringDuplicate("sun00-19-DE405.dat.gz");

  uvar->Freq = 100;
  uvar->f1dot = 0.0;
  uvar->f2dot = 0.0;
  uvar->f3dot = 0.0;
  uvar->h0 = 1;
  uvar->phi0 = 0;

  uvar->startTime.gpsSeconds = 714180733;
  uvar->duration = 10 * 3600;
  uvar->Nseg = 1;
  uvar->segmentList = NULL;

  uvar->refTime.gpsSeconds = -1;	/* default: use mid-time */

  uvar->projection = 0;
  if ( (uvar->IFOs = XLALCreateStringVector ( "H1", NULL )) == NULL ) {
    LogPrintf (LOG_CRITICAL, "Call to XLALCreateStringVector() failed with xlalErrno = %d\n", xlalErrno );
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  uvar->sqrtSX = NULL;

  uvar->detMotionStr = XLALStringDuplicate(XLALDetectorMotionName(DETMOTION_SPIN | DETMOTION_ORBIT));
  uvar->metricType = 0;	/* by default: compute only phase metric */

  if ( (uvar->coords = XLALCreateStringVector ( "freq", "alpha", "delta", "f1dot", NULL )) == NULL ) {
    LogPrintf (LOG_CRITICAL, "Call to XLALCreateStringVector() failed with xlalErrno = %d\n", xlalErrno );
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  uvar->approxPhase = FALSE;

  /* register all our user-variables */

  XLALRegisterUvarMember(IFOs,		STRINGVector, 'I', OPTIONAL, 	"CSV list of detectors, eg. \"H1,H2,L1,G1, ...\" ");
  XLALRegisterUvarMember(sqrtSX,	 	 STRINGVector, 0,  OPTIONAL, 	"[for F-metric weights] CSV list of detectors' noise-floors sqrt{Sn}");
  XLALRegisterUvarMember(Alpha,		RAJ, 'a', OPTIONAL,	"Sky: equatorial J2000 right ascension (in radians or hours:minutes:seconds)");
  XLALRegisterUvarMember(Delta, 		DECJ, 'd', OPTIONAL,	"Sky: equatorial J2000 declination (in radians or degrees:minutes:seconds)");
  XLALRegisterUvarMember(Freq, 		REAL8, 'f', OPTIONAL, 	"Target frequency");
  XLALRegisterUvarMember(f1dot, 		REAL8, 's', OPTIONAL, 	"First spindown-value df/dt");
  XLALRegisterUvarMember(f2dot, 		 REAL8, 0 , OPTIONAL, 	"Second spindown-value d2f/dt2");
  XLALRegisterUvarMember(f3dot, 		 REAL8, 0 , OPTIONAL, 	"Third spindown-value d3f/dt3");

  XLALRegisterUvarMember ( orbitasini,	REAL8, 0, OPTIONAL, 	"Target projected semimajor axis of binary orbit (Units: light seconds)");
  XLALRegisterUvarMember ( orbitPeriod,	REAL8, 0, OPTIONAL, 	"Target period of binary orbit (Units: s).");
  XLALRegisterUvarMember ( orbitTp,	EPOCH, 0, OPTIONAL, 	"Target time of periapse passage of the CW source in a binary orbit (Units: GPS seconds)");
  XLALRegisterUvarMember ( orbitArgp,	REAL8, 0, OPTIONAL, 	"Target argument of periapse of binary orbit (Units: rad)");
  XLALRegisterUvarMember ( orbitEcc,	REAL8, 0, OPTIONAL, 	"Target eccentricity of binary orbit (Units: none)");

  XLALRegisterUvarMember(refTime,         EPOCH, 0,  OPTIONAL, 	"Reference epoch for phase-evolution parameters (format 'xx.yy[GPS|MJD]'). [0=startTime, default=mid-time]");
  XLALRegisterUvarMember(startTime,      EPOCH, 't', OPTIONAL, 	"Start time of observation (format 'xx.yy[GPS|MJD]')");

  XLALRegisterUvarMember(duration,	REAL8, 'T', OPTIONAL,	"Duration of observation in seconds");
  XLALRegisterUvarMember(Nseg,		INT4, 'N', OPTIONAL, 	"Compute semi-coherent metric for this number of segments within 'duration'" );
  XLALRegisterUvarMember(segmentList,   STRING, 0,  OPTIONAL,     "ALTERNATIVE: specify segment file with format: repeated lines <startGPS endGPS duration[h] NumSFTs>");
  XLALRegisterUvarMember( ephemEarth,   STRING, 0,  OPTIONAL,     "Earth ephemeris file to use");
  XLALRegisterUvarMember( ephemSun,     STRING, 0,  OPTIONAL,     "Sun ephemeris file to use");

  XLALRegisterUvarMember(h0,	 	 REAL8, 0, OPTIONAL,	"GW amplitude h0" );
  XLALRegisterUvarMember(cosi,	 	 REAL8, 0, OPTIONAL,	"Pulsar orientation-angle cos(iota) [-1,1]" );
  XLALRegisterUvarMember(psi,		 REAL8, 0, OPTIONAL,	"Wave polarization-angle psi [0, pi]" );
  XLALRegisterUvarMember(phi0,		 REAL8, 0, OPTIONAL,	"GW initial phase phi_0 [0, 2pi]" );

  XLALRegisterUvarMember(metricType,	 INT4, 0,  OPTIONAL,	"type of metric to compute: 0=phase-metric, 1=F-metric(s), 2=both" );
  XLALRegisterUvarMember(outputMetric,	STRING, 'o', OPTIONAL,	"Output the metric components (in octave format) into this file.");
  XLALRegisterUvarMember(projection,      INT4, 0,  OPTIONAL,     "Project onto subspace orthogonal to this axis: 0=none, 1=1st-coord, 2=2nd-coord etc");

  XLALRegisterUvarMember(coords,		STRINGVector, 'c', OPTIONAL, 	"Doppler-coordinates to compute metric in (see --coordsHelp)");
  XLALRegisterUvarMember(coordsHelp,      BOOLEAN, 0,  OPTIONAL,     "output help-string explaining all the possible Doppler-coordinate names for --coords");

  XLALRegisterUvarMember(detMotionStr,  STRING, 0,  DEVELOPER,	"Detector-motion string: S|O|S+O where S=spin|spinz|spinxy and O=orbit|ptoleorbit");
  XLALRegisterUvarMember(approxPhase,     BOOLEAN, 0,  DEVELOPER,	"Use an approximate phase-model, neglecting Roemer delay in spindown coordinates (or orders >= 1)");

  XLALRegisterUvarMember(version,        BOOLEAN, 'V', SPECIAL,      "Output code version");

  return XLAL_SUCCESS;

} /* initUserVars() */


/**
 * basic initializations: set-up 'ConfigVariables'
 */
int
XLALInitCode ( ConfigVariables *cfg, const UserVariables_t *uvar, const char *app_name)
{
  if ( !cfg || !uvar || !app_name ) {
    LogPrintf (LOG_CRITICAL, "%s: illegal NULL pointer input.\n\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL );
  }

  /* Init ephemerides */
  XLAL_CHECK ( (cfg->edat = XLALInitBarycenter ( uvar->ephemEarth, uvar->ephemSun )) != NULL, XLAL_EFUNC );

  // ----- figure out which segments to use
  BOOLEAN manualSegments = XLALUserVarWasSet(&uvar->duration) || XLALUserVarWasSet(&uvar->startTime) || XLALUserVarWasSet(&uvar->Nseg);
  if ( manualSegments && uvar->segmentList ) {
    XLAL_ERROR ( XLAL_EDOM, "Can specify EITHER {--startTime, --duration, --Nseg} OR --segmentList\n");
  }
  LIGOTimeGPS startTimeGPS;
  REAL8 duration;
  if ( uvar->segmentList == NULL )
    {
      XLAL_CHECK ( uvar->Nseg >= 1, XLAL_EDOM, "Invalid input --Nseg=%d: number of segments must be >= 1\n", uvar->Nseg );
      XLAL_CHECK ( uvar->duration >= 1, XLAL_EDOM, "Invalid input --duration=%f: duration must be >= 1 s\n", uvar->duration );
      startTimeGPS = uvar->startTime;
      int ret = XLALSegListInitSimpleSegments ( &cfg->segmentList, startTimeGPS, uvar->Nseg, uvar->duration / uvar->Nseg );
      XLAL_CHECK ( ret == XLAL_SUCCESS, XLAL_EFUNC, "XLALSegListInitSimpleSegments() failed with xlalErrno = %d\n", xlalErrno );
      duration = uvar->duration;
    }
  else
    {
      LALSegList *segList = XLALReadSegmentsFromFile ( uvar->segmentList );
      XLAL_CHECK ( segList != NULL, XLAL_EIO, "XLALReadSegmentsFromFile() failed to load segment list from file '%s', xlalErrno = %d\n", uvar->segmentList, xlalErrno );
      cfg->segmentList = (*segList);	// copy *contents*
      XLALFree ( segList );
      startTimeGPS = cfg->segmentList.segs[0].start;
      UINT4 Nseg = cfg->segmentList.length;
      LIGOTimeGPS endTimeGPS = cfg->segmentList.segs[Nseg-1].end;
      duration = XLALGPSDiff( &endTimeGPS, &startTimeGPS );
    }

  /* ----- figure out reference time */
  LIGOTimeGPS refTimeGPS;

  /* treat special values first */
  if ( uvar->refTime.gpsSeconds == 0 )		/* 0 = use startTime */
    {
      refTimeGPS = uvar->startTime;
    }
  else if ( !XLALUserVarWasSet ( &uvar->refTime ) )	/* default = use mid-time of observation */
    {
      refTimeGPS = startTimeGPS;
      XLALGPSAdd( &refTimeGPS, duration / 2.0 );
    }
  else
    {
      refTimeGPS = uvar->refTime;
    }

  /* ----- get parameter-space point from user-input) */
  cfg->signalParams.Amp.h0 = uvar->h0;
  cfg->signalParams.Amp.cosi = uvar->cosi;
  cfg->signalParams.Amp.psi = uvar->psi;
  cfg->signalParams.Amp.phi0 = uvar->phi0;

  {
    PulsarDopplerParams *dop = &(cfg->signalParams.Doppler);
    XLAL_INIT_MEM((*dop));
    dop->refTime = refTimeGPS;
    dop->Alpha    = uvar->Alpha;
    dop->Delta    = uvar->Delta;
    dop->fkdot[0] = uvar->Freq;
    dop->fkdot[1] = uvar->f1dot;
    dop->fkdot[2] = uvar->f2dot;
    dop->fkdot[3] = uvar->f3dot;
    dop->asini    = uvar->orbitasini;
    dop->period   = uvar->orbitPeriod;
    dop->tp       = uvar->orbitTp;
    dop->ecc      = uvar->orbitEcc;
    dop->argp     = uvar->orbitArgp;
  }

  /* ----- initialize IFOs and (Multi-)DetectorStateSeries  ----- */
  XLAL_CHECK ( XLALParseMultiLALDetector ( &cfg->multiIFO, uvar->IFOs ) == XLAL_SUCCESS, XLAL_EFUNC );
  UINT4 numDet = cfg->multiIFO.length;
  XLAL_CHECK ( numDet >= 1, XLAL_EINVAL );

  if ( uvar->sqrtSX ) {
    XLAL_CHECK ( XLALParseMultiNoiseFloor ( &cfg->multiNoiseFloor, uvar->sqrtSX, numDet ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

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
      LogPrintf (LOG_CRITICAL, "%s: XLALCalloc(1,%zu) failed.\n\n", __func__, sizeof(*cfg->history));
      XLAL_ERROR ( XLAL_ENOMEM );
    }

    if ( (tmp = XLALMalloc ( len )) == NULL ) {
      LogPrintf (LOG_CRITICAL, "%s: XLALMalloc (%zu) failed.\n\n", __func__, len );
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

  XLALSegListClear ( &(cfg->segmentList) );
  return XLAL_SUCCESS;

} /* XLALDestroyConfig() */


int
XLALOutputDopplerMetric ( FILE *fp, const DopplerPhaseMetric *Pmetric, const DopplerFstatMetric *Fmetric, const ResultHistory_t *history )
{
  UINT4 i;
  REAL8 A, B, C, D;

  // ----- input sanity checks
  XLAL_CHECK ( fp != NULL, XLAL_EFAULT );
  XLAL_CHECK ( Pmetric != NULL || Fmetric != NULL, XLAL_EFAULT );
  const DopplerMetricParams *meta = (Pmetric != NULL) ? &(Pmetric->meta) : &(Fmetric->meta);
  XLAL_CHECK ( XLALSegListIsInitialized ( &(meta->segmentList) ), XLAL_EINVAL, "Got un-initialized segment list in 'metric->meta.segmentList'\n" );
  UINT4 Nseg = meta->segmentList.length;
  XLAL_CHECK ( Nseg >= 1, XLAL_EDOM, "Got invalid zero-length segment list 'metric->meta.segmentList'\n" );

  /* useful shortcuts */
  const PulsarDopplerParams *doppler = &(meta->signalParams.Doppler);
  const PulsarAmplitudeParams *Amp = &(meta->signalParams.Amp);

  /* output history info */
  if ( history )
    {
      if ( history->app_name ) fprintf (fp, "%%%% app_name: %s\n", history->app_name );
      if ( history->cmdline) fprintf (fp, "%%%% commandline: %s\n", history->cmdline );
      if ( history->VCSInfoString ) fprintf (fp, "%%%% Code Version: %s\n", history->VCSInfoString );
    }

  fprintf ( fp, "DopplerCoordinates = { " );
  for ( i=0; i < meta->coordSys.dim; i ++ )
    {
      if ( i > 0 ) fprintf ( fp, ", " );
      fprintf ( fp, "\"%s\"", XLALDopplerCoordinateName(meta->coordSys.coordIDs[i]));
    }
  fprintf ( fp, "};\n");

  { /* output projection info */
    const char *pname;
    if ( meta->projectCoord < 0 )
      pname = "None";
    else
      pname = XLALDopplerCoordinateName ( meta->coordSys.coordIDs[meta->projectCoord] );

    fprintf ( fp, "%%%% Projection onto subspace orthogonal to coordinate: '%s'\n", pname);
  }

  fprintf ( fp, "%%%% DetectorMotionType = '%s'\n", XLALDetectorMotionName(meta->detMotionType) );
  fprintf ( fp, "h0 = %g;\ncosi = %g;\npsi = %g;\nphi0 = %g;\n", Amp->h0, Amp->cosi, Amp->psi, Amp->phi0 );
  fprintf ( fp, "%%%% DopplerPoint = {\n");
  fprintf ( fp, "refTime = %.1f;\n", XLALGPSGetREAL8 ( &doppler->refTime ) );
  fprintf ( fp, "Alpha   = %f;\nDelta = %f;\n", doppler->Alpha, doppler->Delta );
  fprintf ( fp, "fkdot   = [%f, %g, %g, %g ];\n", doppler->fkdot[0], doppler->fkdot[1], doppler->fkdot[2], doppler->fkdot[3] );
  if ( doppler->asini > 0 )
    {
      fprintf ( fp, "%%%% 	   orbit = { \n");
      fprintf ( fp, "%%%% 		tp = {%d, %d}\n", doppler->tp.gpsSeconds, doppler->tp.gpsNanoSeconds );
      fprintf ( fp, "%%%% 		argp  = %g\n", doppler->argp );
      fprintf ( fp, "%%%% 		asini = %g\n", doppler->asini );
      fprintf ( fp, "%%%% 		ecc = %g\n", doppler->ecc );
      fprintf ( fp, "%%%% 		period = %g\n", doppler->period );
      fprintf ( fp, "%%%% 	   }\n");
    } /* if doppler->orbit */
  fprintf ( fp, "%%%% }\n");

  LIGOTimeGPS *tStart = &(meta->segmentList.segs[0].start);
  LIGOTimeGPS *tEnd   = &(meta->segmentList.segs[Nseg-1].end);
  REAL8 Tspan = XLALGPSDiff ( tEnd, tStart );
  fprintf ( fp, "startTime = %.1f;\n", XLALGPSGetREAL8 ( tStart ) );
  fprintf ( fp, "Tspan     = %.1f;\n", Tspan );
  fprintf ( fp, "Nseg      = %d;\n", Nseg );
  fprintf ( fp, "detectors = {");
  for ( i=0; i < meta->multiIFO.length; i ++ )
    {
      if ( i > 0 ) fprintf ( fp, ", ");
      fprintf ( fp, "\"%s\"", meta->multiIFO.sites[i].frDetector.name );
    }
  fprintf ( fp, "};\n");
  fprintf ( fp, "detectorWeights = [");
  for ( i=0; i < meta->multiNoiseFloor.length; i ++ )
    {
      if ( i > 0 ) fprintf ( fp, ", ");
      fprintf ( fp, "%f", meta->multiNoiseFloor.sqrtSn[i] );
    }
  fprintf ( fp, "];\n");

  /* ----- output phase metric ---------- */
  if ( Pmetric != NULL )
    {
      fprintf ( fp, "\ng_ij = \\\n" ); XLALfprintfGSLmatrix ( fp, METRIC_FORMAT,  Pmetric->g_ij );
      fprintf ( fp, "maxrelerr_gPh = %.2e;\n", Pmetric->maxrelerr );

      gsl_matrix *gN_ij = NULL;
      if ( XLALNaturalizeMetric ( &gN_ij, NULL, Pmetric->g_ij, meta ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: something failed Naturalizing phase metric g_ij!\n", __func__ );
        XLAL_ERROR ( XLAL_EFUNC );
      }
      fprintf ( fp, "\ngN_ij = \\\n" ); XLALfprintfGSLmatrix ( fp, METRIC_FORMAT,  gN_ij );
      gsl_matrix_free ( gN_ij );

      gsl_matrix *gDN_ij = NULL;
      if ( XLALDiagNormalizeMetric ( &gDN_ij, NULL, Pmetric->g_ij ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: something failed NormDiagonalizing phase metric g_ij!\n", __func__ );
        XLAL_ERROR ( XLAL_EFUNC );
      }
      fprintf ( fp, "\ngDN_ij = \\\n" ); XLALfprintfGSLmatrix ( fp, METRIC_FORMAT,  gDN_ij );
      gsl_matrix_free ( gDN_ij );
    }

  /* ----- output F-metric (and related matrices ---------- */
  if ( Fmetric != NULL )
    {
      fprintf ( fp, "\ngF_ij = \\\n" );   XLALfprintfGSLmatrix ( fp, METRIC_FORMAT,  Fmetric->gF_ij );
      fprintf ( fp, "\ngFav_ij = \\\n" ); XLALfprintfGSLmatrix ( fp, METRIC_FORMAT,  Fmetric->gFav_ij );
      fprintf ( fp, "\nm1_ij = \\\n" );   XLALfprintfGSLmatrix ( fp, METRIC_FORMAT,  Fmetric->m1_ij );
      fprintf ( fp, "\nm2_ij = \\\n" );   XLALfprintfGSLmatrix ( fp, METRIC_FORMAT,  Fmetric->m2_ij );
      fprintf ( fp, "\nm3_ij = \\\n" );   XLALfprintfGSLmatrix ( fp, METRIC_FORMAT,  Fmetric->m3_ij );
      fprintf ( fp, "maxrelerr_gF = %.2e;\n", Fmetric->maxrelerr );
    }

  /*  ----- output Fisher matrix ---------- */
  if ( Fmetric != NULL && Fmetric->Fisher_ab != NULL )
    {
      A = gsl_matrix_get ( Fmetric->Fisher_ab, 0, 0 );
      B = gsl_matrix_get ( Fmetric->Fisher_ab, 1, 1 );
      C = gsl_matrix_get ( Fmetric->Fisher_ab, 0, 1 );

      D = A * B - C * C;

      fprintf ( fp, "\nA = %.16g;\nB = %.16g;\nC = %.16g;\nD = %.16g;\n", A, B, C, D );
      fprintf ( fp, "\nrho2 = %.16g;\n", Fmetric->rho2 );

      fprintf (fp, "\nFisher_ab = \\\n" ); XLALfprintfGSLmatrix ( fp, METRIC_FORMAT,  Fmetric->Fisher_ab );
    }

  // ---------- output segment list at the end, as this can potentially become quite long and distracting
  char *seglist_octave;
  XLAL_CHECK ( (seglist_octave = XLALSegList2String ( &(meta->segmentList) )) != NULL, XLAL_EFUNC, "XLALSegList2String() with xlalErrno = %d\n", xlalErrno );
  fprintf ( fp, "\n\nsegmentList = %s;\n", seglist_octave );
  XLALFree ( seglist_octave );


  return XLAL_SUCCESS;

} /* XLALOutputDopplerMetric() */


/**
 * Destructor for ResultHistory type
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
