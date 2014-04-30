/*
 * Copyright (C) 2012 David Keitel, Santiago Caride, Reinhard Prix
 * Copyright (C) 2008, 2013 Karl Wette
 * Copyright (C) 2007 Chris Messenger
 * Copyright (C) 2004, 2007, 2010 Reinhard Prix
 * Copyright (C) 2005, 2006 Reinhard Prix, Iraj Gholami
 * Copyright (C) 2002, 2003, 2004 M.A. Papa, X. Siemens, Y. Itoh
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

/*********************************************************************************/
/**
 * \author R. Prix, I. Gholami, Y. Ioth, Papa, X. Siemens, C. Messenger, K. Wette
 * \file
 * \ingroup pulsarApps
 * \brief
 * Calculate the F-statistic for a given parameter-space of pulsar GW signals.
 * Implements the so-called "F-statistic" as introduced in \cite JKS98.
 *
 * This code is a descendant of an earlier implementation 'ComputeFStatistic.[ch]'
 * by Bruce Allen, Bernd Machenschalk, David Hammer, Jolien Creighton, Maria Alessandra Papa,
 * Reinhard Prix, Xavier Siemens, Scott Koranda, Yousuke Itoh
 *
 */
#include "config.h"

/* System includes */
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <strings.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <sys/time.h>

int finite(double);

/* LAL-includes */
#include <lal/LALString.h>
#include <lal/AVFactories.h>
#include <lal/GSLSupport.h>
#include <lal/LALInitBarycenter.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/ExtrapolatePulsarSpins.h>

#include <lal/NormalizeSFTRngMed.h>
#include <lal/ComputeFstat.h>
#include <lal/LALHough.h>

#include <lal/LogPrintf.h>
#include <lal/DopplerFullScan.h>
#include <lal/BinaryPulsarTiming.h>

#include <lal/TransientCW_utils.h>
#include <lal/LineRobustStats.h>

#include <lalapps.h>

/* local includes */
#include "HeapToplist.h"

/*---------- DEFINES ----------*/

#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename */

#define TRUE (1==1)
#define FALSE (1==0)

/*----- SWITCHES -----*/
#define NUM_SPINS 4		/* number of spin-values to consider: {f, fdot, f2dot, ... } */

/*----- Error-codes -----*/
#define COMPUTEFSTATISTIC_ENULL 	1
#define COMPUTEFSTATISTIC_ESYS     	2
#define COMPUTEFSTATISTIC_EINPUT   	3
#define COMPUTEFSTATISTIC_EMEM   	4
#define COMPUTEFSTATISTIC_ENONULL 	5
#define COMPUTEFSTATISTIC_EXLAL		6

#define COMPUTEFSTATISTIC_MSGENULL 	"Arguments contained an unexpected null pointer"
#define COMPUTEFSTATISTIC_MSGESYS	"System call failed (probably file IO)"
#define COMPUTEFSTATISTIC_MSGEINPUT   	"Invalid input"
#define COMPUTEFSTATISTIC_MSGEMEM   	"Out of memory. Bad."
#define COMPUTEFSTATISTIC_MSGENONULL 	"Output pointer is non-NULL"
#define COMPUTEFSTATISTIC_MSGEXLAL	"XLALFunction-call failed"

/*----- Macros -----*/

/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )
#define SQ(x) ( (x) * (x) )

#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

#define LAL_INT4_MAX 2147483647

#define INIT_MEM(x) memset(&(x), 0, sizeof((x)))

/*---------- internal types ----------*/

/** What info do we want to store in our toplist? */
typedef struct {
  PulsarDopplerParams doppler;		/**< Doppler params of this 'candidate' */
  REAL4 twoF;				/**< F-statistic value */
  UINT4 numDetectors;			/**< number of detectors = effective vector length. numDetectors=0 should make all code ignore the FX field. */
  REAL4 twoFX[PULSAR_MAX_DETECTORS];	/**< vector of single-detector F-statistic values (array of fixed size) */
  LIGOTimeGPS FaFb_refTime;		/**< reference time of Fa and Fb */
  COMPLEX16 Fa;				/**< complex amplitude Fa */
  COMPLEX16 Fb;				/**< complex amplitude Fb */
  AntennaPatternMatrix Mmunu;		/**< antenna-pattern matrix Mmunu = (h_mu|h_nu) */
  REAL4 LVstat;				/**< Line Veto statistic */
} FstatCandidate;

/**
 * moving 'Scanline window' of candidates on the scan-line,
 * which is used to find local 1D maxima.
 */
typedef struct
{
  UINT4 length;
  FstatCandidate *window; 		/**< array holding candidates */
  FstatCandidate *center;		/**< pointer to middle candidate in window */
} scanlineWindow_t;

/**
 * Struct holding various timing measurements and relevant search parameters.
 * This is used to fit timing-models with measured times to predict search run-times
 */
typedef struct
{
  UINT4 NSFTs;			/**< total number of SFTs */
  REAL8 tauFstat;		/**< time to compute one Fstatistic over full data-duration (NSFT atoms) [in seconds]*/
  REAL8 tauTemplate;		/**< total loop time per template, includes candidate-handling (transient stats, toplist etc) */

  /* transient-specific timings */
  UINT4 tauMin;			/**< shortest transient timescale [s] */
  UINT4 tauMax;			/**< longest transient timescale [s] */
  UINT4 NStart;			/**< number of transient start-time steps in FstatMap matrix */
  UINT4 NTau;			/**< number of transient timescale steps in FstatMap matrix */

  REAL8 tauTransFstatMap;	/**< time to compute transient-search Fstatistic-map over {t0, tau} [s]     */
  REAL8 tauTransMarg;		/**< time to marginalize the Fstat-map to compute transient-search Bayes [s] */
} timingInfo_t;

/**
 * Enum for which statistic is used to "rank" significance of candidates
 */
typedef enum {
  RANKBY_2F  = 0, 	/**< rank candidates by F-statistic */
  RANKBY_NC = 1,	/**< HierarchSearchGCT also has RANKBY_NC = 1, not applicable here */
  RANKBY_LV = 2  	/**< rank candidates by LV-statistic */
} RankingStat_t;


/**
 * Configuration settings required for and defining a coherent pulsar search.
 * These are 'pre-processed' settings, which have been derived from the user-input.
 */
typedef struct {
  LIGOTimeGPS startTime;		    /**< starting timestamp of SFTs */
  REAL8 Tspan;				    /**< time spanned by all SFTs */
  REAL8 Alpha;                              /**< sky position alpha in radians */
  REAL8 Delta;                              /**< sky position delta in radians */
  REAL8 Tsft;                               /**< length of one SFT in seconds */
  LIGOTimeGPS internalRefTime;	            /**< internal reference time used purely for F-stat computation (defaults to midTime) */
  DopplerRegion searchRegion;		    /**< parameter-space region to search over */
  DopplerFullScanState *scanState;          /**< current state of the Doppler-scan */
  PulsarDopplerParams stepSizes;	    /**< user-preferences on Doppler-param step-sizes */
  EphemerisData *ephemeris;		    /**< ephemeris data (from LALInitBarycenter()) */
  UINT4 NSFTs;				    /**< total number of all SFTs used */
  UINT4Vector *numSFTsPerDet;		    /**< number of SFTs per detector, for log strings, etc. */
  LALStringVector *detectorIDs;		    /**< detector ID names, for column headings string */
  FstatInput *Fstat_in;		    /**< Fstat input data struct */
  FstatQuantities Fstat_what;		    /**< Fstat quantities to compute */
  toplist_t* FstatToplist;		    /**< sorted 'toplist' of the NumCandidatesToKeep loudest candidates */
  scanlineWindow_t *scanlineWindow;         /**< moving window of candidates on scanline to find local maxima */
  CHAR *VCSInfoString;                      /**< LAL + LALapps Git version string */
  CHAR *logstring;                          /**< log containing max-info on the whole search setup */
  transientWindowRange_t transientWindowRange; /**< search range parameters for transient window */
  REAL8 LVlogRhoTerm;                       /**< log(rho^4/70) of LV line-prior amplitude 'rho' */
  REAL8Vector *LVloglX;                     /**< vector of line-prior ratios per detector {l1, l2, ... } */
  RankingStat_t RankingStatistic;           /**< rank candidates according to F or LV */
} ConfigVariables;


/* ----- User-variables: can be set from config-file or command-line */
typedef struct {
  INT4 Dterms;			/**< number of terms in LALDemod Dirichlet kernel is Dterms+1 */
  CHAR *IFO;			/**< IFO name: only required if using v1 SFTs */
  BOOLEAN SignalOnly;		/**< FALSE: estimate noise-floor from data, TRUE: assume Sh=1 */
  BOOLEAN UseNoiseWeights;	/**< use SFT-specific noise-weights for each segment in Fstat-computation */

  REAL8 Freq;			/**< start-frequency of search */
  REAL8 FreqBand;		/**< Frequency-band to search over */
  REAL8 dFreq;			/**< user-specifyable Frequency stepsize */

  REAL8 Alpha;			/**< equatorial right-ascension in rad */
  CHAR *RA;
  REAL8 dAlpha;
  REAL8 AlphaBand;

  REAL8 Delta;			/**< equatorial declination in rad */
  CHAR *Dec;
  REAL8 dDelta;
  REAL8 DeltaBand;

  REAL8 f1dot;			/**< 1st spindown dFreq/dt */
  REAL8 df1dot;
  REAL8 f1dotBand;

  REAL8 f2dot;			/**< 2nd spindown d^2Freq/dt^2 */
  REAL8 df2dot;
  REAL8 f2dotBand;

  REAL8 f3dot;			/**< 3rd spindown d^3Freq/dt^3 */
  REAL8 df3dot;
  REAL8 f3dotBand;

  /* orbital parameters */
  REAL8 orbitPeriod;		/**< binary-system orbital period in s */
  REAL8 orbitasini;		/**< amplitude of radial motion */
  INT4 orbitTpSSBsec;		/**< time of periapse passage */
  INT4 orbitTpSSBnan;
  REAL8 orbitTpSSBMJD;		/**< in MJD format */
  REAL8 orbitArgp;		/**< angle of periapse */
  REAL8 orbitEcc;		/**< orbital eccentricity */

  /* extra parameters for --gridType=GRID_SPINDOWN_AGEBRK parameter space */
  REAL8 spindownAge;            /**< spindown age of the object */
  REAL8 minBraking;             /**< minimum braking index */
  REAL8 maxBraking;             /**< maximum braking index */

  /* --- */
  REAL8 TwoFthreshold;		/**< output threshold on 2F */
  CHAR *ephemEarth;		/**< Earth ephemeris file to use */
  CHAR *ephemSun;		/**< Sun ephemeris file to use */

  INT4  gridType;		/**< type of template grid in parameter space */
  INT4  metricType;		/**< type of metric to use in metric template grids */
  BOOLEAN projectMetric;	/**< should we use frequency-projected metric? */
  REAL8 metricMismatch;		/**< maximal *nominal* metric mismatch *per* dimension */
  CHAR *skyRegion;		/**< list of skypositions defining a search-polygon */
  CHAR *DataFiles;		/**< glob-pattern for SFT data-files to use */

  BOOLEAN help;			/**< output help-string */
  CHAR *outputLogfile;		/**< write a log-file */
  CHAR *outputFstat;		/**< filename to output Fstatistic in */
  CHAR *outputLoudest;		/**< filename for loudest F-candidate plus parameter estimation */
  CHAR *outputLogPrintf;        /**< send output from LogPrintf statements to this file */

  CHAR *outputFstatHist;        /**< output discrete histogram of all Fstatistic values */
  REAL8 FstatHistBin;           /**< width of an Fstatistic histogram bin */

  BOOLEAN countTemplates;       /**< just count templates (if supported) instead of search */

  INT4 NumCandidatesToKeep;	/**< maximal number of toplist candidates to output */
  REAL8 FracCandidatesToKeep;	/**< fractional number of candidates to output in toplist */
  INT4 clusterOnScanline;	/**< number of points on "scanline" to use for 1-D local maxima finding */

  CHAR *gridFile;		/**< read template grid from this file */
  REAL8 dopplermax;		/**< maximal possible doppler-effect */

  INT4 RngMedWindow;		/**< running-median window for noise floor estimation */
  REAL8 refTime;		/**< reference-time for definition of pulsar-parameters [GPS] */
  REAL8 refTimeMJD;		/**< the same in MJD */

  REAL8 internalRefTime;	/**< which reference time to use internally for template-grid */
  INT4 SSBprecision;		/**< full relativistic timing or Newtonian */

  INT4 minStartTime;		/**< earliest start-time to use data from */
  INT4 maxEndTime;		/**< latest end-time to use data from */
  CHAR *workingDir;		/**< directory to use for output files */
  REAL8 timerCount;		/**< output progress-meter every timerCount seconds */

  BOOLEAN version;		/**< output version information */

  CHAR *outputFstatAtoms;	/**< output per-SFT, per-IFO 'atoms', ie quantities required to compute F-stat */

  BOOLEAN outputSingleFstats;	/**< in multi-detector case, also output single-detector F-stats */
  BOOLEAN computeLV;		/**< get single-IFO F-stats and compute Line Veto stat */
  CHAR *RankingStatistic;	/**< rank candidates according to F or LV */
  BOOLEAN LVuseAllTerms;	/**< Use only leading term or all terms in Line Veto computation */
  REAL8   LVrho;		/**< Prior parameter rho_max_line for LineVeto statistic */
  LALStringVector *LVlX;	/**< Line-to-gauss prior ratios lX for LineVeto statistic */
  REAL8 LVthreshold;		/**< output threshold on LV */

  CHAR *outputTransientStats;	/**< output file for transient B-stat values */
  CHAR *transient_WindowType;	/**< name of transient window ('none', 'rect', 'exp',...) */
  REAL8 transient_t0Days;	/**< earliest GPS start-time for transient window search, as offset in days from dataStartGPS */
  REAL8 transient_t0DaysBand;	/**< Range of GPS start-times to search in transient search, in days */
  INT4  transient_dt0;		/**< Step-size for search/marginalization over transient-window start-time, in seconds */
  REAL8 transient_tauDays;	/**< smallest transient window length for marginalization, in days */
  REAL8 transient_tauDaysBand;	/**< Range of transient-window timescales to search, in days */
  INT4  transient_dtau;		/**< Step-size for search/marginalization over transient-window timescale, in seconds */
  BOOLEAN transient_useFReg;  	/**< FALSE: use 'standard' e^F for marginalization, TRUE: use e^FReg = (1/D)*e^F */

  CHAR *outputTiming;		/**< output timing measurements and parameters into this file [append!]*/

  BOOLEAN useResamp;		/**< use FFT-resampling method instead of LALDemod() */
} UserInput_t;

/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

/* empty initializers */
static UserInput_t empty_UserInput;
static timingInfo_t empty_timingInfo;

/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);
void initUserVars (LALStatus *, UserInput_t *uvar);
void InitFStat ( LALStatus *, ConfigVariables *cfg, const UserInput_t *uvar );
void Freemem(LALStatus *,  ConfigVariables *cfg);

void checkUserInputConsistency (LALStatus *, const UserInput_t *uvar);
int outputBeamTS( const CHAR *fname, const AMCoeffs *amcoe, const DetectorStateSeries *detStates );
void getUnitWeights ( LALStatus *, MultiNoiseWeights **multiWeights, const MultiSFTVector *multiSFTs );

int write_FstatCandidate_to_fp ( FILE *fp, const FstatCandidate *thisFCand );
int write_PulsarCandidate_to_fp ( FILE *fp,  const PulsarCandidate *pulsarParams, const FstatCandidate *Fcand );

int compareFstatCandidates ( const void *candA, const void *candB );
int compareFstatCandidates_LV ( const void *candA, const void *candB );

void WriteFStatLog ( LALStatus *status, const CHAR *log_fname, const CHAR *logstr );
CHAR *XLALGetLogString ( const ConfigVariables *cfg );

int write_TimingInfo_to_fp ( FILE * fp, const timingInfo_t *ti );

/* ---------- scanline window functions ---------- */
scanlineWindow_t *XLALCreateScanlineWindow ( UINT4 windowWings );
void XLALDestroyScanlineWindow ( scanlineWindow_t *scanlineWindow );
int XLALAdvanceScanlineWindow ( const FstatCandidate *nextCand, scanlineWindow_t *scanWindow );
BOOLEAN XLALCenterIsLocalMax ( const scanlineWindow_t *scanWindow, const UINT4 sortingStatistic );

/*---------- empty initializers ---------- */
static const ConfigVariables empty_ConfigVariables;
static const FstatCandidate empty_FstatCandidate;

/* ----- which timing function to use ----- */
#ifdef HIGHRES_TIMING
REAL8 XLALGetUserCPUTime ( void );
#define GETTIME XLALGetUserCPUTime
#else
#define GETTIME XLALGetTimeOfDay
#endif

/*----------------------------------------------------------------------*/
/* Function definitions start here */
/*----------------------------------------------------------------------*/

/**
 * MAIN function of ComputeFStatistic code.
 * Calculate the F-statistic over a given portion of the parameter-space
 * and write a list of 'candidates' into a file(default: 'Fstats').
 */
int main(int argc,char *argv[])
{
  LALStatus status = blank_status;	/* initialize status */

  FILE *fpFstat = NULL, *fpTransientStats = NULL;
  REAL8 numTemplates, templateCounter;
  time_t clock0;
  PulsarDopplerParams dopplerpos = empty_PulsarDopplerParams;		/* current search-parameters */
  FstatCandidate loudestFCand = empty_FstatCandidate, thisFCand = empty_FstatCandidate;
  FILE *fpLogPrintf = NULL;
  gsl_vector_int *Fstat_histogram = NULL;

  UserInput_t uvar = empty_UserInput;
  ConfigVariables GV = empty_ConfigVariables;		/**< global container for various derived configuration settings */

  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register all user-variable */
  LAL_CALL (initUserVars(&status, &uvar), &status);

  if ( (GV.VCSInfoString = XLALGetVersionString(0)) == NULL ) {
    XLALPrintError("XLALGetVersionString(0) failed.\n");
    exit(1);
  }

  /* do ALL cmdline and cfgfile handling */
  LAL_CALL (LALUserVarReadAllInput(&status, argc, argv), &status);

  if (uvar.help)	/* if help was requested, we're done here */
    exit (0);

  if ( uvar.version )
    {
      printf ( "%s\n", GV.VCSInfoString );
      exit (0);
    }

  /* set log-level and open log-file */
  LogSetLevel ( lalDebugLevel );
  if (LALUserVarWasSet(&uvar.outputLogPrintf)) {
    if ((fpLogPrintf = fopen(uvar.outputLogPrintf, "wb")) == NULL) {
      XLALPrintError ("\nError opening file '%s' for writing..\n\n", uvar.outputLogPrintf);
      return (COMPUTEFSTATISTIC_ESYS);
    }
    LogSetFile(fpLogPrintf);
  }

  /* do some sanity checks on the user-input before we proceed */
  LAL_CALL ( checkUserInputConsistency(&status, &uvar), &status);

  /* Initialization the common variables of the code, */
  /* like ephemeries data and template grids: */
  LAL_CALL ( InitFStat(&status, &GV, &uvar), &status);

  /* ----- produce a log-string describing the specific run setup ----- */
  if ( (GV.logstring = XLALGetLogString ( &GV )) == NULL ) {
    XLALPrintError ( "XLALGetLogString() failed!\n");
    return COMPUTEFSTATISTIC_EXLAL;
  }
  LogPrintfVerbatim( LOG_DEBUG, "%s", GV.logstring );

  /* keep a log-file recording all relevant parameters of this search-run */
  if ( uvar.outputLogfile ) {
    LAL_CALL (WriteFStatLog ( &status, uvar.outputLogfile, GV.logstring ), &status );
  }

  /* if a complete output of the F-statistic file was requested,
   * we open and prepare the output-file here */
  if (uvar.outputFstat)
    {
      if ( (fpFstat = fopen (uvar.outputFstat, "wb")) == NULL)
	{
	  XLALPrintError ("\nError opening file '%s' for writing..\n\n", uvar.outputFstat);
	  return (COMPUTEFSTATISTIC_ESYS);
	}

      fprintf (fpFstat, "%s", GV.logstring );

      /* assemble column headings string */
      char colum_headings_string_base[] = "freq alpha delta f1dot f2dot f3dot 2F";
      UINT4 column_headings_string_length = sizeof(colum_headings_string_base);
      char column_headings_string[column_headings_string_length];
      INIT_MEM( column_headings_string );
      strcat ( column_headings_string, colum_headings_string_base );
      if ( uvar.computeLV )
        {
          const UINT4 numDetectors = GV.detectorIDs->length;
          column_headings_string_length += numDetectors*6; /* 6 per detector for " 2F_XY" */
          column_headings_string_length += 3; /* 3 for " LV"*/
          strcat ( column_headings_string, " LV" );
          for ( UINT4 X = 0; X < numDetectors ; X ++ )
            {
              char headingX[7];
              snprintf ( headingX, sizeof(headingX), " 2F_%s", GV.detectorIDs->data[X] );
              strcat ( column_headings_string, headingX );
            } /* for X < numDet */
        }
      fprintf (fpFstat, "%%%% columns:\n%%%% %s\n", column_headings_string );

    } /* if outputFstat */

  if ( uvar.outputTransientStats )
    {
      if ( (fpTransientStats = fopen (uvar.outputTransientStats, "wb")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar.outputTransientStats );
	  return (COMPUTEFSTATISTIC_ESYS);
	}

      fprintf (fpTransientStats, "%s", GV.logstring );			/* write search log comment */
      write_transientCandidate_to_fp ( fpTransientStats, NULL );	/* write header-line comment */
    }

  /* start Fstatistic histogram with a single empty bin */
  if (uvar.outputFstatHist) {
    if ((Fstat_histogram = gsl_vector_int_alloc(1)) == NULL) {
      XLALPrintError("\nCouldn't allocate 'Fstat_histogram'\n");
      return COMPUTEFSTATISTIC_EMEM;
    }
    gsl_vector_int_set_zero(Fstat_histogram);
  }

  /* setup binary parameters */
  REAL8 orbit_asini = 0 /* isolated pulsar */;
  REAL8 orbit_period = 0;
  REAL8 orbit_ecc = 0;
  LIGOTimeGPS orbit_tp = LIGOTIMEGPSZERO;
  REAL8 orbit_argp = 0;
  if ( LALUserVarWasSet(&uvar.orbitasini) && (uvar.orbitasini > 0) )
    {
      orbit_tp.gpsSeconds = uvar.orbitTpSSBsec;
      orbit_tp.gpsNanoSeconds = uvar.orbitTpSSBnan;
      orbit_argp = uvar.orbitArgp;
      orbit_asini = uvar.orbitasini;
      orbit_ecc = uvar.orbitEcc;
      orbit_period = uvar.orbitPeriod;
      if (LALUserVarWasSet(&uvar.orbitTpSSBMJD))
	{
	  /* convert MJD peripase to GPS using Matt Pitkins code found at lal/packages/pulsar/src/BinaryPulsarTimeing.c */
	  REAL8 GPSfloat;
	  GPSfloat = XLALTTMJDtoGPS(uvar.orbitTpSSBMJD);
	  XLALGPSSetREAL8(&(orbit_tp),GPSfloat);
	}
      else
	{
	  orbit_tp.gpsSeconds = uvar.orbitTpSSBsec;
	  orbit_tp.gpsNanoSeconds = uvar.orbitTpSSBnan;
	}
    }

  /* count number of templates */
  numTemplates = XLALNumDopplerTemplates ( GV.scanState );
  if (uvar.countTemplates)
    printf("%%%% Number of templates: %0.0f\n", numTemplates);

  /*----------------------------------------------------------------------
   * main loop: demodulate data for each point in the sky-position grid
   * and for each value of the frequency-spindown
   */
  templateCounter = 0.0;
  clock0 = time(NULL);

  REAL8 tic0, tic, toc, timeOfLastProgressUpdate = 0;	// high-precision timing counters
  timingInfo_t timing = empty_timingInfo;	// timings of Fstatistic computation, transient Fstat-map, transient Bayes factor
  timing.NSFTs = GV.NSFTs;

  /* fixed time-offset between internalRefTime and refTime */
  REAL8 DeltaTRefInt = XLALGPSDiff ( &(GV.internalRefTime), &(GV.searchRegion.refTime) ); // tRefInt - tRef

  UINT4 numFreqBins_FBand = 1;	// number of frequency-bins in the frequency-band used for resampling (1 if not --useResamp)
  REAL8 dFreqResamp = 0; // frequency resolution used to allocate vector of F-stat values for resampling
  if ( uvar.useResamp )	// handle special resampling case, where we deal with a vector of F-stat values instead of one
    {
	if ( LALUserVarWasSet(&uvar.dFreq) ) {
		dFreqResamp = uvar.dFreq;
	} else {
		dFreqResamp = 1.0/(2*GV.Tspan);
	}
      numFreqBins_FBand = (UINT4) ( 1 + floor ( GV.searchRegion.fkdotBand[0] / dFreqResamp ) );
    }

  // pointer to Fstat results structure, will be allocated by XLALComputeFstat()
  FstatResults* Fstat_res = NULL;

  /* skip search if user supplied --countTemplates */
  while ( !uvar.countTemplates && (XLALNextDopplerPos( &dopplerpos, GV.scanState ) == 0) )
    {
      /* temporary solution until binary-gridding exists */
      dopplerpos.asini  = orbit_asini;
      dopplerpos.period = orbit_period;
      dopplerpos.ecc    = orbit_ecc;
      dopplerpos.tp     = orbit_tp;
      dopplerpos.argp   = orbit_argp;

      tic0 = tic = GETTIME();

      /* use internalRefTime in order to safely computing F-statistic (avoid large |t - tRef|^s) */
      PulsarDopplerParams internalDopplerpos = dopplerpos;
      XLALExtrapolatePulsarSpins ( internalDopplerpos.fkdot, dopplerpos.fkdot, DeltaTRefInt );	// can't fail
      internalDopplerpos.refTime = GV.internalRefTime;

      /* main function call: compute F-statistic for this template */
      const int retn = XLALComputeFstat(&Fstat_res, GV.Fstat_in, &internalDopplerpos, dFreqResamp, numFreqBins_FBand, GV.Fstat_what);
      if ( retn != XLAL_SUCCESS ) {
        XLALPrintError ("%s: XLALComputeFstat() failed with errno=%d\n", __func__, xlalErrno );
        return xlalErrno;
      }
      /* if single-only flag is given, add +4 to F-statistic */
      if ( uvar.SignalOnly ) {
        if (XLALAdd4ToFstatResults(Fstat_res) != XLAL_SUCCESS) {
          XLALPrintError ("%s: XLALAdd4ToFstatResults() failed with errno=%d\n", __func__, xlalErrno );
          return xlalErrno;
        }
      }

      toc = GETTIME();
      timing.tauFstat += (toc - tic);   // pure Fstat-calculation time

      /* Progress meter */
      templateCounter += 1.0;
      if ( lalDebugLevel && ( (toc - timeOfLastProgressUpdate) > uvar.timerCount) )
        {
          REAL8 diffSec = time(NULL) - clock0 ;  /* seconds since start of loop*/
          REAL8 taup = diffSec / templateCounter ;
          REAL8 timeLeft = (numTemplates - templateCounter) *  taup;
          LogPrintf (LOG_DEBUG, "Progress: %g/%g = %.2f %% done, Estimated time left: %.0f s\n",
                     templateCounter, numTemplates, templateCounter/numTemplates * 100.0, timeLeft);
          timeOfLastProgressUpdate = toc;
        }

      // here we use Santiago's trick to hack the resampling Fstat(f) into the single-F rest of the
      // main-loop: we simply loop the remaining body over all frequency-bins in the Fstat-vector,
      // this way nothing needs to be changed!  in the non-resampling case, this loop iterates only
      // once, so nothing is changed ...
      for ( UINT4 iFreq = 0; iFreq < numFreqBins_FBand; iFreq ++ )
      {

        /* collect data on current 'Fstat-candidate' */
        thisFCand.doppler = dopplerpos;	// use 'original' dopplerpos @ refTime !
        thisFCand.doppler.fkdot[0] += iFreq * dFreqResamp; // this only does something for the resampling post-loop over frequency-bins, 0 otherwise ...
        thisFCand.twoF = Fstat_res->twoF[iFreq];
        if (GV.Fstat_what & FSTATQ_2F_PER_DET) {
          thisFCand.numDetectors = Fstat_res->numDetectors;
          for (UINT4 X = 0; X < thisFCand.numDetectors; ++X ) {
            thisFCand.twoFX[X] = Fstat_res->twoFPerDet[X][iFreq];
          }
        } else {
          thisFCand.numDetectors = 0;
        }
        if (GV.Fstat_what & FSTATQ_FAFB) {
          thisFCand.FaFb_refTime = Fstat_res->doppler.refTime; // this will be 'internal' reference time, used only for parameter estimation
          thisFCand.Fa = Fstat_res->FaFb[iFreq].Fa;
          thisFCand.Fb = Fstat_res->FaFb[iFreq].Fb;
        } else {
          thisFCand.Fa = thisFCand.Fb = crect(NAN,NAN);
        }
        thisFCand.Mmunu = Fstat_res->Mmunu;
        MultiFstatAtomVector* thisFAtoms = NULL;
        if (GV.Fstat_what & FSTATQ_ATOMS_PER_DET) {
          thisFAtoms = Fstat_res->multiFatoms[iFreq];
        }

        /* sanity check on the result */
        if ( !finite(thisFCand.twoF) )
	{
	  LogPrintf(LOG_CRITICAL, "non-finite 2F = %.16g, Fa=(%.16g,%.16g), Fb=(%.16g,%.16g)\n",
		    thisFCand.twoF, creal(thisFCand.Fa), cimag(thisFCand.Fa), creal(thisFCand.Fb), cimag(thisFCand.Fb) );
	  LogPrintf (LOG_CRITICAL, "[Alpha,Delta] = [%.16g,%.16g],\nfkdot=[%.16g,%.16g,%.16g,%16.g]\n",
		     thisFCand.doppler.Alpha, thisFCand.doppler.Delta,
		     thisFCand.doppler.fkdot[0], thisFCand.doppler.fkdot[1], thisFCand.doppler.fkdot[2], thisFCand.doppler.fkdot[3] );
	  if (thisFCand.doppler.asini > 0)
          {
            LogPrintf (LOG_CRITICAL, "tp = {%d s, %d ns}, argp = %f, asini = %f, ecc = %f, period = %f\n",
                       thisFCand.doppler.tp.gpsSeconds, thisFCand.doppler.tp.gpsNanoSeconds,
                       thisFCand.doppler.argp, thisFCand.doppler.asini,
                       thisFCand.doppler.ecc, thisFCand.doppler.period);
          }
	  return -1;
	}

      if ( uvar.computeLV )
        {
          REAL8* LVlogLX = NULL;
          if ( GV.LVloglX ) {
            LVlogLX = GV.LVloglX->data;
          }
          thisFCand.LVstat = XLALComputeLineVetoArray ( thisFCand.twoF, thisFCand.numDetectors, thisFCand.twoFX, GV.LVlogRhoTerm, LVlogLX, uvar.LVuseAllTerms );
          if ( xlalErrno ) {
            XLALPrintError ("%s: XLALComputeLineVetoArray() failed with errno=%d\n", __func__, xlalErrno );
            return xlalErrno;
          }
        }
      else
        {
          thisFCand.LVstat = NAN; /* in non-LV case, block field with NAN, needed for output checking in write_PulsarCandidate_to_fp() */
        }

      /* push new value onto scan-line buffer */
      XLALAdvanceScanlineWindow ( &thisFCand, GV.scanlineWindow );

      /* two types of threshold: fixed (TwoF- and/or LV-threshold) and dynamic (NumCandidatesToKeep) */
      BOOLEAN is1DlocalMax = FALSE;
      if ( XLALCenterIsLocalMax ( GV.scanlineWindow, GV.RankingStatistic ) ) /* must be 1D local maximum */
        is1DlocalMax = TRUE;
      BOOLEAN isOver2FThreshold = FALSE; /* will always be checked, so start at 'FALSE' */
      if ( GV.scanlineWindow->center->twoF >= uvar.TwoFthreshold ) /* fixed 2F threshold */
        isOver2FThreshold = TRUE;
      BOOLEAN isOverLVThreshold = TRUE;  /* will not be checked in non-LV case, so start at 'TRUE' */
      if ( uvar.computeLV && ( GV.scanlineWindow->center->LVstat < uvar.LVthreshold ) ) /* fixed LV threshold */
        isOverLVThreshold = FALSE;
      if ( is1DlocalMax && isOver2FThreshold && isOverLVThreshold )
        {
	  FstatCandidate *writeCand = GV.scanlineWindow->center;

	  /* insert this into toplist if requested */
	  if ( GV.FstatToplist  )			/* dynamic threshold */
	    {
	      if ( insert_into_toplist(GV.FstatToplist, (void*)writeCand ) ) {
		LogPrintf ( LOG_DETAIL, "Added new candidate into toplist: 2F = %f", writeCand->twoF );
		if ( uvar.computeLV )
		  LogPrintfVerbatim ( LOG_DETAIL, ", 2F_H1 = %f, 2F_L1 = %f, LV = %f", writeCand->twoFX[0], writeCand->twoFX[1], writeCand->LVstat );
	      }
	      else {
		LogPrintf ( LOG_DETAIL, "NOT added the candidate into toplist: 2F = %f", writeCand->twoF );
		if ( uvar.computeLV )
		  LogPrintfVerbatim ( LOG_DETAIL, ", 2F_H1 = %f, 2F_L1 = %f, LV = %f", writeCand->twoFX[0], writeCand->twoFX[1], writeCand->LVstat );
	      }
	      LogPrintfVerbatim ( LOG_DETAIL, "\n" );
	    }
	  else if ( fpFstat ) 				/* no toplist :write out immediately */
	    {
	      if ( write_FstatCandidate_to_fp ( fpFstat, writeCand ) != 0 )
		{
		  LogPrintf (LOG_CRITICAL, "Failed to write candidate to file.\n");
		  return -1;
		}
	    } /* if outputFstat */

	} /* if candidate is local maximum and over threshold */

      /* separately keep track of loudest candidate (for --outputLoudest) */
      switch ( GV.RankingStatistic )
        {
        case RANKBY_2F:
          if ( thisFCand.twoF > loudestFCand.twoF )
            loudestFCand = thisFCand;
          break;
        case RANKBY_LV:
          if ( thisFCand.LVstat > loudestFCand.LVstat )
            loudestFCand = thisFCand;
          break;
        default:
          XLAL_ERROR ( XLAL_EINVAL, "Invalid ranking statistic '%d', supported are 'F=0', and 'LV=2'\n", GV.RankingStatistic );
          break;
        }

      /* add Fstatistic to histogram if needed */
      if (uvar.outputFstatHist) {

	/* compute bin */
	const size_t bin = thisFCand.twoF / uvar.FstatHistBin;

	/* resize histogram vector if needed */
	if (!Fstat_histogram || bin >= Fstat_histogram->size)
	  if (NULL == (Fstat_histogram = XLALResizeGSLVectorInt(Fstat_histogram, bin + 1, 0))) {
	    XLALPrintError("\nCouldn't (re)allocate 'Fstat_histogram'\n");
	    return COMPUTEFSTATISTIC_EMEM;
	  }

	/* add to bin */
	gsl_vector_int_set(Fstat_histogram, bin,
			   gsl_vector_int_get(Fstat_histogram, bin) + 1);

      }


      /* ----- output F-stat atoms into files: one file per Doppler-position ---------- */
      /* F-stat 'atoms' = per-SFT components {a,b,Fa,Fb}_alpha) */
      if (uvar.outputFstatAtoms)
	{
	  FILE *fpFstatAtoms = NULL;
	  CHAR *fnameAtoms = NULL;
	  CHAR *dopplerName;
	  UINT4 len;

	  if ( (dopplerName = XLALPulsarDopplerParams2String ( &dopplerpos )) == NULL )
	    {
	      return COMPUTEFSTATISTIC_EXLAL;
	    }
	  len = strlen(uvar.outputFstatAtoms) + strlen(dopplerName) + 10;
	  if ( (fnameAtoms = LALMalloc (len)) == NULL )
	    {
	      LogPrintf( LOG_CRITICAL, "Failed to LALMalloc(%d)\n", len );
	      return COMPUTEFSTATISTIC_EMEM;
	    }
	  sprintf (fnameAtoms, "%s_%s.dat", uvar.outputFstatAtoms, dopplerName );

	  if ( (fpFstatAtoms = fopen (fnameAtoms, "wb")) == NULL)
	    {
	      XLALPrintError ("\n%s: Error opening file '%s' for writing..\n\n", __func__, fnameAtoms );
	      return COMPUTEFSTATISTIC_ESYS;
	    }
	  LALFree ( fnameAtoms );
	  LALFree ( dopplerName );

	  fprintf (fpFstatAtoms, "%s", GV.logstring );

	  if ( write_MultiFstatAtoms_to_fp ( fpFstatAtoms, thisFAtoms ) != XLAL_SUCCESS ) {
            XLALPrintError ("%s: failed to write atoms to output file. xlalErrno = %d\n", __func__, xlalErrno );
            return COMPUTEFSTATISTIC_ESYS;
          }

	  fclose (fpFstatAtoms);

	} /* if outputFstatAtoms */

      /* ----- compute transient-CW statistics if their output was requested  ----- */
      if ( fpTransientStats )
        {
          transientCandidate_t transientCand = empty_transientCandidate;

          /* compute Fstat map F_mn over {t0, tau} */
          tic = GETTIME();
          if ( (transientCand.FstatMap = XLALComputeTransientFstatMap ( thisFAtoms, GV.transientWindowRange, uvar.transient_useFReg)) == NULL ) {
            XLALPrintError ("%s: XLALComputeTransientFstatMap() failed with xlalErrno = %d.\n", __func__, xlalErrno );
            return COMPUTEFSTATISTIC_EXLAL;
          }
          toc = GETTIME();
          timing.tauTransFstatMap += (toc - tic); // time to compute transient Fstat-map

          /* compute marginalized Bayes factor */
          tic = GETTIME();
          transientCand.logBstat = XLALComputeTransientBstat ( GV.transientWindowRange, transientCand.FstatMap );
          UINT4 err = xlalErrno;
          if ( err ) {
            XLALPrintError ("%s: XLALComputeTransientBstat() failed with xlalErrno = %d\n", __func__, err );
            return COMPUTEFSTATISTIC_EXLAL;
          }
          toc = GETTIME();
          timing.tauTransMarg += (toc - tic);

          /* ----- compute parameter posteriors for {t0, tau} */
          pdf1D_t *pdf_t0  = NULL;
          pdf1D_t *pdf_tau = NULL;

          if ( (pdf_t0 = XLALComputeTransientPosterior_t0 ( GV.transientWindowRange, transientCand.FstatMap )) == NULL ) {
              XLALPrintError ("%s: failed to compute t0-posterior\n", __func__ );
              XLAL_ERROR ( xlalErrno );
          }
          if ( (pdf_tau = XLALComputeTransientPosterior_tau ( GV.transientWindowRange, transientCand.FstatMap )) == NULL ) {
              XLALPrintError ("%s: failed to compute tau-posterior\n", __func__ );
              XLAL_ERROR ( xlalErrno );
          }
          /* get maximum-posterior estimate (MP) from the modes of these pdfs */
          transientCand.t0_MP = XLALFindModeOfPDF1D ( pdf_t0 );
          if ( xlalErrno ) {
              XLALPrintError ("%s: mode-estimation failed for pdf_t0. xlalErrno = %d\n", __func__, xlalErrno );
              XLAL_ERROR ( xlalErrno );
          }
          transientCand.tau_MP =  XLALFindModeOfPDF1D ( pdf_tau );
          if ( xlalErrno ) {
              XLALPrintError ("%s: mode-estimation failed for pdf_tau. xlalErrno = %d\n", __func__, xlalErrno );
              XLAL_ERROR ( xlalErrno );
          }

          /* record timing-relevant transient search params */
          timing.tauMin  = GV.transientWindowRange.tau;
          timing.tauMax  = timing.tauMin + GV.transientWindowRange.tauBand;
          timing.NStart  = transientCand.FstatMap->F_mn->size1;
          timing.NTau    = transientCand.FstatMap->F_mn->size2;

          /* add meta-info on current transient-CW candidate */
          transientCand.doppler = dopplerpos;
          transientCand.windowRange = GV.transientWindowRange;

          /* correct for missing bias in signal-only case */
          if ( uvar.SignalOnly )
            transientCand.FstatMap->maxF += 2;

          /* output everything into stats-file (one line per candidate) */
          if ( write_transientCandidate_to_fp ( fpTransientStats, &transientCand ) != XLAL_SUCCESS ) {
            XLALPrintError ("%s: write_transientCandidate_to_fp() failed.\n", __func__ );
            return COMPUTEFSTATISTIC_EXLAL;
          }

          /* free dynamically allocated F-stat map */
          XLALDestroyTransientFstatMap ( transientCand.FstatMap );
          XLALDestroyPDF1D ( pdf_t0 );
          XLALDestroyPDF1D ( pdf_tau );

        } /* if fpTransientStats */

      } // for ( iFreq < numFreqBins_FBand )

      /* now measure total loop time per template */
      toc = GETTIME();
      timing.tauTemplate += (toc - tic0);

    } /* while more Doppler positions to scan */


  /* if requested: output timings into timing-file */
  if ( uvar.outputTiming )
    {
      FILE *fpTiming = NULL;
      REAL8 num_templates = numTemplates * numFreqBins_FBand;	// 'templates' now refers to number of 'frequency-bands' in resampling case

      // compute averages:
      timing.tauFstat         /= num_templates;
      timing.tauTemplate      /= num_templates;
      timing.tauTransFstatMap /= num_templates;
      timing.tauTransMarg     /= num_templates;

      if ( ( fpTiming = fopen ( uvar.outputTiming, "ab" )) == NULL ) {
        XLALPrintError ("%s: failed to open timing file '%s' for writing \n", __func__, uvar.outputTiming );
        return COMPUTEFSTATISTIC_ESYS;
      }
      /* write header comment line */
      if ( write_TimingInfo_to_fp ( fpTiming, NULL ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: write_TimingInfo_to_fp() failed to write header-comment into file '%s'\n", __func__, uvar.outputTiming );
        return COMPUTEFSTATISTIC_EXLAL;
      }

      if ( write_TimingInfo_to_fp ( fpTiming, &timing ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: write_TimingInfo_to_fp() failed.\n", __func__ );
        return COMPUTEFSTATISTIC_EXLAL;
      }

      fclose ( fpTiming );
    } /* if timing output requested */

  /* ----- if using toplist: sort and write it out to file now ----- */
  if ( fpFstat && GV.FstatToplist )
    {
      UINT4 el;

      /* sort toplist */
      LogPrintf ( LOG_DEBUG, "Sorting toplist ... ");
      if ( GV.RankingStatistic == RANKBY_2F )
        qsort_toplist ( GV.FstatToplist, compareFstatCandidates );
      else if ( GV.RankingStatistic == RANKBY_LV )
        qsort_toplist ( GV.FstatToplist, compareFstatCandidates_LV );
      else
        XLAL_ERROR ( XLAL_EINVAL, "Ranking statistic '%d' undefined here, allowed are 'F=0' and 'LV=2'\n", GV.RankingStatistic );
      LogPrintfVerbatim ( LOG_DEBUG, "done.\n");

      for ( el=0; el < GV.FstatToplist->elems; el ++ )
	{
	  const FstatCandidate *candi;
	  if ( ( candi = (const FstatCandidate *) toplist_elem ( GV.FstatToplist, el )) == NULL ) {
	    LogPrintf ( LOG_CRITICAL, "Internal consistency problems with toplist: contains fewer elements than expected!\n");
	    return -1;
	  }
	  if ( write_FstatCandidate_to_fp ( fpFstat, candi ) != 0 )
	    {
	      LogPrintf (LOG_CRITICAL, "Failed to write candidate to file.\n");
	      return -1;
	    }
	} /* for el < elems in toplist */

    } /* if fpFstat && toplist */

  if ( fpFstat )
    {
      fprintf (fpFstat, "%%DONE\n");
      fclose (fpFstat);
      fpFstat = NULL;
    }
  if ( fpTransientStats ) fclose (fpTransientStats);

  /* ----- estimate amplitude-parameters for the loudest canidate and output into separate file ----- */
  if ( uvar.outputLoudest )
    {
      FILE *fpLoudest;
      PulsarCandidate pulsarParams = empty_PulsarCandidate;
      pulsarParams.Doppler = loudestFCand.doppler;

      if ( XLALEstimatePulsarAmplitudeParams ( &pulsarParams, &loudestFCand.FaFb_refTime, loudestFCand.Fa, loudestFCand.Fb, &loudestFCand.Mmunu )
           != XLAL_SUCCESS )
      {
        XLALPrintError ("%s: XLALEstimatePulsarAmplitudeParams() failed with errno=%d\n", __func__, xlalErrno );
        return COMPUTEFSTATISTIC_ESYS;
      }

      if ( (fpLoudest = fopen (uvar.outputLoudest, "wb")) == NULL)
	{
	  XLALPrintError ("\nError opening file '%s' for writing..\n\n", uvar.outputLoudest);
	  return COMPUTEFSTATISTIC_ESYS;
	}

      /* write header with run-info */
      fprintf (fpLoudest, "%s", GV.logstring );

      /* write this 'candidate' to disc */
      if ( write_PulsarCandidate_to_fp ( fpLoudest,  &pulsarParams, &loudestFCand) != XLAL_SUCCESS )
	{
	  LogPrintf(LOG_CRITICAL, "call to write_PulsarCandidate_to_fp() failed!\n");
	  return COMPUTEFSTATISTIC_ESYS;
	}

      fclose (fpLoudest);

      gsl_matrix_free ( pulsarParams.AmpFisherMatrix );

    } /* write loudest candidate to file */

  LogPrintf (LOG_DEBUG, "Search finished.\n");

  /* write out the Fstatistic histogram */
  if (uvar.outputFstatHist) {

    size_t i = 0;
    FILE *fpFstatHist = fopen(uvar.outputFstatHist, "wb");

    if (fpFstatHist == NULL) {
      XLALPrintError ("\nError opening file '%s' for writing..\n\n", uvar.outputFstat);
      return (COMPUTEFSTATISTIC_ESYS);
    }
    fprintf(fpFstatHist, "%s", GV.logstring);

    for (i = 0; i < Fstat_histogram->size; ++i)
      fprintf(fpFstatHist, "%0.3f %0.3f %i\n",
	      uvar.FstatHistBin * i,
	      uvar.FstatHistBin * (i + 1),
	      gsl_vector_int_get(Fstat_histogram, i));

    fprintf(fpFstatHist, "%%DONE\n");
    fclose(fpFstatHist);

  }

  /* Free memory */
  LogPrintf (LOG_DEBUG, "Freeing Doppler grid ... ");
  LAL_CALL ( FreeDopplerFullScan(&status, &GV.scanState), &status);
  LogPrintfVerbatim ( LOG_DEBUG, "done.\n");

  XLALDestroyFstatResults ( Fstat_res );

  LAL_CALL ( Freemem(&status, &GV), &status);

  if (Fstat_histogram)
    gsl_vector_int_free(Fstat_histogram);

  /* close log-file */
  if (fpLogPrintf) {
    fclose(fpLogPrintf);
    LogSetFile(fpLogPrintf = NULL);
  }

  /* did we forget anything ? */
  LALCheckMemoryLeaks();

  return 0;

} /* main() */


/**
 * Register all our "user-variables" that can be specified from cmd-line and/or config-file.
 * Here we set defaults for some user-variables and register them with the UserInput module.
 */
void
initUserVars (LALStatus *status, UserInput_t *uvar)
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* set a few defaults */
  uvar->FreqBand = 0.0;
  uvar->Alpha 	= 0.0;
  uvar->Delta 	= 0.0;
  uvar->RA       = NULL;
  uvar->Dec      = NULL;
  uvar->AlphaBand = 0;
  uvar->DeltaBand = 0;
  uvar->skyRegion = NULL;
  // Dterms-default used to be 16, but has to be 8 for SSE version
  uvar->Dterms 	= OptimisedHotloopDterms;

  uvar->ephemEarth = XLALStringDuplicate("earth00-19-DE405.dat.gz");
  uvar->ephemSun = XLALStringDuplicate("sun00-19-DE405.dat.gz");

  uvar->SignalOnly = FALSE;
  uvar->UseNoiseWeights = TRUE;

  /* default step-sizes for GRID_FLAT */
  uvar->dAlpha 	= 0.001;
  uvar->dDelta 	= 0.001;
  uvar->dFreq 	 = 0.0;		/* zero indicates 'not set by user==> i.e. use canonical default */
  uvar->df1dot    = 0.0;
  uvar->df2dot    = 0.0;
  uvar->df3dot    = 0.0;

  /* define default orbital semi-major axis */
  uvar->orbitasini = 0 /* isolated pulsar */;

  uvar->TwoFthreshold = 0.0;
  uvar->NumCandidatesToKeep = 0;
  uvar->FracCandidatesToKeep = 0.0;
  uvar->clusterOnScanline = 0;

  uvar->metricType =  LAL_PMETRIC_NONE;
  uvar->projectMetric = TRUE;
  uvar->gridType = GRID_FLAT;

  uvar->metricMismatch = 0.02;

  uvar->help = FALSE;
  uvar->version = FALSE;

  uvar->outputLogfile = NULL;
  uvar->outputFstat = NULL;
  uvar->outputLoudest = NULL;
  uvar->outputLogPrintf = NULL;

  uvar->outputFstatHist = NULL;
  uvar->FstatHistBin = 0.1;

  uvar->countTemplates = FALSE;

  uvar->gridFile = NULL;

  uvar->dopplermax =  1.05e-4;
  uvar->RngMedWindow = 50;	/* for running-median */

  uvar->SSBprecision = SSBPREC_RELATIVISTIC;

  uvar->minStartTime = 0;
  uvar->maxEndTime = LAL_INT4_MAX;

  uvar->workingDir = (CHAR*)LALMalloc(512);
  strcpy(uvar->workingDir, ".");

  uvar->timerCount = 10;	/* output a timer/progress count every N seconds */

  uvar->spindownAge = 0.0;
  uvar->minBraking = 0.0;
  uvar->maxBraking = 0.0;

  uvar->useResamp = FALSE;

  uvar->outputSingleFstats = FALSE;
  uvar->computeLV = FALSE;
  #define DEFAULT_RANKINGSTATISTIC "F"
  uvar->RankingStatistic = LALCalloc (1, strlen(DEFAULT_RANKINGSTATISTIC)+1);
  strcpy (uvar->RankingStatistic, DEFAULT_RANKINGSTATISTIC);
  uvar->LVuseAllTerms = TRUE;
  uvar->LVrho = 0.0;
  uvar->LVlX = NULL;       /* NULL is intepreted as LVlX[X] = 1.0 for all X by XLALComputeLineVeto(Array) */
  uvar->LVthreshold = - LAL_REAL8_MAX;

#define DEFAULT_TRANSIENT "none"
  uvar->transient_WindowType = LALMalloc(strlen(DEFAULT_TRANSIENT)+1);
  strcpy ( uvar->transient_WindowType, DEFAULT_TRANSIENT );
  uvar->transient_useFReg = 0;

  /* ---------- register all user-variables ---------- */
  LALregBOOLUserStruct(status, 	help, 		'h', UVAR_HELP,     "Print this message");

  LALregREALUserStruct(status, 	Alpha, 		'a', UVAR_OPTIONAL, "Sky position alpha (equatorial coordinates) in radians");
  LALregREALUserStruct(status, 	Delta, 		'd', UVAR_OPTIONAL, "Sky position delta (equatorial coordinates) in radians");
  LALregSTRINGUserStruct(status,RA, 		 0 , UVAR_OPTIONAL, "Sky position alpha (equatorial coordinates) in format hh:mm:ss.sss");
  LALregSTRINGUserStruct(status,Dec, 		 0 , UVAR_OPTIONAL, "Sky position delta (equatorial coordinates) in format dd:mm:ss.sss");
  LALregSTRINGUserStruct(status,skyRegion, 	'R', UVAR_OPTIONAL, "ALTERNATIVE: Sky-region by polygon of form '(ra1,dec1),(ra2,dec2),(ra3,dec3),...' or 'allsky'");

  LALregREALUserStruct(status, 	Freq, 		'f', UVAR_OPTIONAL, "Starting search frequency Freq in Hz");
  LALregREALUserStruct(status, 	f1dot, 		's', UVAR_OPTIONAL, "First spindown parameter  f1dot = dFreq/dt");
  LALregREALUserStruct(status, 	f2dot, 		 0 , UVAR_OPTIONAL, "Second spindown parameter f2dot = d^2Freq/dt^2");
  LALregREALUserStruct(status, 	f3dot, 		 0 , UVAR_OPTIONAL, "Third spindown parameter  f3dot = d^3Freq/dt^3");

  LALregREALUserStruct(status, 	AlphaBand, 	'z', UVAR_OPTIONAL, "Search band in alpha (equatorial coordinates) in radians");
  LALregREALUserStruct(status, 	DeltaBand, 	'c', UVAR_OPTIONAL, "Search band in delta (equatorial coordinates) in radians");
  LALregREALUserStruct(status, 	FreqBand, 	'b', UVAR_OPTIONAL, "Search band in frequency in Hz");
  LALregREALUserStruct(status, 	f1dotBand, 	'm', UVAR_OPTIONAL, "Search band in f1dot in Hz/s");
  LALregREALUserStruct(status, 	f2dotBand, 	 0 , UVAR_OPTIONAL, "Search band in f2dot in Hz/s^2");
  LALregREALUserStruct(status, 	f3dotBand, 	 0 , UVAR_OPTIONAL, "Search band in f3dot in Hz/s^3");

  LALregREALUserStruct(status, 	dAlpha, 	'l', UVAR_OPTIONAL, "Stepsize in alpha (equatorial coordinates) in radians");
  LALregREALUserStruct(status, 	dDelta, 	'g', UVAR_OPTIONAL, "Stepsize in delta (equatorial coordinates) in radians");
  LALregREALUserStruct(status,  dFreq,          'r', UVAR_OPTIONAL, "Stepsize for frequency in Hz");
  LALregREALUserStruct(status, 	df1dot, 	'e', UVAR_OPTIONAL, "Stepsize for f1dot in Hz/s");
  LALregREALUserStruct(status, 	df2dot, 	 0 , UVAR_OPTIONAL, "Stepsize for f2dot in Hz/s^2");
  LALregREALUserStruct(status, 	df3dot, 	 0 , UVAR_OPTIONAL, "Stepsize for f3dot in Hz/s^3");

  LALregREALUserStruct(status, 	orbitasini, 	 0,  UVAR_OPTIONAL, "Binary Orbit: Projected semi-major axis in light-seconds [Default: 0.0]");
  LALregREALUserStruct(status, 	orbitPeriod, 	 0,  UVAR_OPTIONAL, "Binary Orbit: Period in seconds");
  LALregINTUserStruct(status, 	orbitTpSSBsec, 	 0,  UVAR_OPTIONAL, "Binary Orbit: (true) time of periapsis in SSB frame, GPS seconds");
  LALregINTUserStruct(status, 	orbitTpSSBnan, 	 0,  UVAR_OPTIONAL, "Binary Orbit: (true) time of periapsis in SSB frame, GPS nanoseconds part");
  LALregREALUserStruct(status, 	orbitTpSSBMJD, 	 0,  UVAR_OPTIONAL, "ALTERNATIVE: (true) time of periapsis in the SSB frame in MJD");
  LALregREALUserStruct(status, 	orbitArgp, 	 0,  UVAR_OPTIONAL, "Binary Orbit: Orbital argument of periapse in radians");
  LALregREALUserStruct(status, 	orbitEcc, 	 0,  UVAR_OPTIONAL, "Binary Orbit: Orbital eccentricity");

  LALregSTRINGUserStruct(status,DataFiles, 	'D', UVAR_REQUIRED, "File-pattern specifying (also multi-IFO) input SFT-files");
  LALregSTRINGUserStruct(status,IFO, 		'I', UVAR_OPTIONAL, "Detector: 'G1', 'L1', 'H1', 'H2' ...(useful for single-IFO v1-SFTs only!)");

  LALregBOOLUserStruct(status, 	SignalOnly, 	'S', UVAR_OPTIONAL, "Signal only flag");

  LALregREALUserStruct(status, 	TwoFthreshold,	'F', UVAR_OPTIONAL, "Set the threshold for selection of 2F");
  LALregINTUserStruct(status, 	gridType,	 0 , UVAR_OPTIONAL, "Grid: 0=flat, 1=isotropic, 2=metric, 3=skygrid-file, 6=grid-file, 7=An*lattice, 8=spin-square, 9=spin-age-brk");
  LALregINTUserStruct(status, 	metricType,	'M', UVAR_OPTIONAL, "Metric: 0=none,1=Ptole-analytic,2=Ptole-numeric, 3=exact");
  LALregREALUserStruct(status, 	metricMismatch,	'X', UVAR_OPTIONAL, "Maximal allowed mismatch for metric tiling");
  LALregSTRINGUserStruct(status,outputLogfile,	 0,  UVAR_OPTIONAL, "Name of log-file identifying the code + search performed");
  LALregSTRINGUserStruct(status,gridFile,	 0,  UVAR_OPTIONAL, "Load grid from this file: sky-grid or full-grid depending on --gridType.");
  LALregREALUserStruct(status,	refTime,	 0,  UVAR_OPTIONAL, "SSB reference time for pulsar-parameters [Default: startTime]");
  LALregREALUserStruct(status,	refTimeMJD,	 0,  UVAR_OPTIONAL, "ALTERNATIVE: SSB reference time for pulsar-parameters in MJD [Default: startTime]");

  LALregSTRINGUserStruct(status,outputFstat,	 0,  UVAR_OPTIONAL, "Output-file for F-statistic field over the parameter-space");
  LALregSTRINGUserStruct(status,outputLoudest,	 0,  UVAR_OPTIONAL, "Loudest F-statistic candidate + estimated MLE amplitudes");

  LALregSTRINGUserStruct(status,outputFstatHist, 0,  UVAR_OPTIONAL, "Output-file for a discrete histogram of all Fstatistic values");
  LALregREALUserStruct(status,  FstatHistBin,    0,  UVAR_OPTIONAL, "Width of an Fstatistic histogram bin");

  LALregINTUserStruct(status,  NumCandidatesToKeep,0, UVAR_OPTIONAL, "Number of Fstat 'candidates' to keep. (0 = All)");
  LALregREALUserStruct(status,FracCandidatesToKeep,0, UVAR_OPTIONAL, "Fraction of Fstat 'candidates' to keep.");
  LALregINTUserStruct(status,   clusterOnScanline, 0, UVAR_OPTIONAL, "Neighbors on each side for finding 1D local maxima on scanline");

  LALregINTUserStruct ( status, minStartTime, 	 0,  UVAR_OPTIONAL, "Earliest SFT-timestamp to include");
  LALregINTUserStruct ( status, maxEndTime, 	 0,  UVAR_OPTIONAL, "Latest SFT-timestamps to include");

  LALregSTRINGUserStruct(status,outputFstatAtoms,0,  UVAR_OPTIONAL, "Output filename *base* for F-statistic 'atoms' {a,b,Fa,Fb}_alpha. One file per doppler-point.");
  LALregBOOLUserStruct(status,  outputSingleFstats,0,  UVAR_OPTIONAL, "In multi-detector case, also output single-detector F-stats?");
  LALregBOOLUserStruct(status,  computeLV,	0,  UVAR_OPTIONAL, "Get single-detector F-stats and compute Line Veto statistic.");
  LALregSTRINGUserStruct(status,RankingStatistic,0,  UVAR_DEVELOPER, "Rank toplist candidates according to 'F' or 'LV' statistic");
  LALregREALUserStruct(status,  LVrho,		0,  UVAR_OPTIONAL, "LineVeto: Prior rho_max_line, must be >=0");
  LALregLISTUserStruct(status,  LVlX,		0,  UVAR_OPTIONAL, "LineVeto: line-to-gauss prior ratios lX for different detectors X, length must be numDetectors. Defaults to lX=1,1,..");
  LALregREALUserStruct(status, 	LVthreshold,	0,  UVAR_OPTIONAL, "Set the threshold for selection of LV");

  LALregSTRINGUserStruct(status,outputTransientStats,0,  UVAR_OPTIONAL, "TransientCW: Output filename for transient-CW statistics.");
  LALregSTRINGUserStruct(status, transient_WindowType,0,UVAR_OPTIONAL,  "TransientCW: Type of transient signal window to use. ('none', 'rect', 'exp').");
  LALregREALUserStruct (status, transient_t0Days, 0,  UVAR_OPTIONAL,    "TransientCW: Earliest GPS start-time for transient window search, as offset in days from dataStartGPS");
  LALregREALUserStruct (status, transient_t0DaysBand,0,UVAR_OPTIONAL,   "TransientCW: Range of GPS start-times to search in transient search, in days");
  LALregINTUserStruct (status, transient_dt0,    0,  UVAR_OPTIONAL,     "TransientCW: Step-size in transient-CW start-time in seconds [Default:Tsft]");
  LALregREALUserStruct(status, transient_tauDays,0,  UVAR_OPTIONAL,     "TransientCW: Minimal transient-CW duration timescale, in days");
  LALregREALUserStruct(status, transient_tauDaysBand,0,  UVAR_OPTIONAL, "TransientCW: Range of transient-CW duration timescales to search, in days");
  LALregINTUserStruct (status, transient_dtau,   0,  UVAR_OPTIONAL,     "TransientCW: Step-size in transient-CW duration timescale, in seconds [Default:Tsft]");

  LALregBOOLUserStruct( status, version,	'V', UVAR_SPECIAL,  "Output version information");

  LALregBOOLUserStruct(status,  useResamp,       0,  UVAR_OPTIONAL,  "Use FFT-resampling method instead of LALDemod()");

  /* ----- more experimental/expert options ----- */
  LALregREALUserStruct(status, 	dopplermax, 	'q', UVAR_DEVELOPER, "Maximum doppler shift expected");
  LALregBOOLUserStruct(status, 	UseNoiseWeights,'W', UVAR_DEVELOPER, "Use per-SFT noise weights");
  LALregSTRINGUserStruct(status,ephemEarth, 	 0,  UVAR_DEVELOPER, "Earth ephemeris file to use");
  LALregSTRINGUserStruct(status,ephemSun, 	 0,  UVAR_DEVELOPER, "Sun ephemeris file to use");

  LALregINTUserStruct (status, 	SSBprecision,	 0,  UVAR_DEVELOPER, "Precision to use for time-transformation to SSB: 0=Newtonian 1=relativistic");

  LALregINTUserStruct(status, 	RngMedWindow,	'k', UVAR_DEVELOPER, "Running-Median window size");
  LALregINTUserStruct(status,	Dterms,		't', UVAR_DEVELOPER, "Number of terms to keep in Dirichlet kernel sum");

  LALregSTRINGUserStruct(status,workingDir,     'w', UVAR_DEVELOPER, "Directory to use as work directory.");
  LALregREALUserStruct(status, 	timerCount, 	 0,  UVAR_DEVELOPER, "N: Output progress/timer info every N seconds");
  LALregREALUserStruct(status,	internalRefTime, 0,  UVAR_DEVELOPER, "internal reference time to use for Fstat-computation [Default: midTime]");

  LALregBOOLUserStruct(status, 	projectMetric, 	 0,  UVAR_DEVELOPER, "Use projected metric on Freq=const subspact");

  LALregSTRINGUserStruct(status,outputLogPrintf, 0,  UVAR_DEVELOPER, "Send all output from LogPrintf statements to this file");

  LALregBOOLUserStruct(status, 	countTemplates,  0,  UVAR_DEVELOPER, "Count number of templates (if supported) instead of search");

  LALregREALUserStruct(status,  spindownAge,     0,  UVAR_DEVELOPER, "Spindown age for --gridType=9");
  LALregREALUserStruct(status,  minBraking,      0,  UVAR_DEVELOPER, "Minimum braking index for --gridType=9");
  LALregREALUserStruct(status,  maxBraking,      0,  UVAR_DEVELOPER, "Maximum braking index for --gridType=9");

  XLALregBOOLUserStruct ( transient_useFReg,   	 0,  UVAR_DEVELOPER, "FALSE: use 'standard' e^F for marginalization, if TRUE: use e^FReg = (1/D)*e^F (BAD)");

  XLALregSTRINGUserStruct( outputTiming,         0,  UVAR_DEVELOPER, "Append timing measurements and parameters into this file");

  LALregBOOLUserStruct(status,LVuseAllTerms	,0,  UVAR_DEVELOPER, "Use only leading term or all terms in Line Veto computation?");

  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* initUserVars() */

/** Initialized Fstat-code: handle user-input and set everything up.
 * NOTE: the logical *order* of things in here is very important, so be careful
 */
void
InitFStat ( LALStatus *status, ConfigVariables *cfg, const UserInput_t *uvar )
{
  REAL8 fCoverMin, fCoverMax;	/* covering frequency-band to read from SFTs */
  SFTCatalog *catalog = NULL;
  SFTConstraints constraints = empty_SFTConstraints;
  LIGOTimeGPS minStartTimeGPS = empty_LIGOTimeGPS;
  LIGOTimeGPS maxEndTimeGPS = empty_LIGOTimeGPS;

  LIGOTimeGPS endTime;
  size_t toplist_length = uvar->NumCandidatesToKeep;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* set the current working directory */
  if(chdir(uvar->workingDir) != 0)
    {
      LogPrintf (LOG_CRITICAL,  "Unable to change directory to workinDir '%s'\n", uvar->workingDir);
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }

  /* use IFO-contraint if one given by the user */
  if ( LALUserVarWasSet ( &uvar->IFO ) )
    if ( (constraints.detector = XLALGetChannelPrefix ( uvar->IFO )) == NULL ) {
      ABORT ( status,  COMPUTEFSTATISTIC_EINPUT,  COMPUTEFSTATISTIC_MSGEINPUT);
    }
  minStartTimeGPS.gpsSeconds = uvar->minStartTime;
  maxEndTimeGPS.gpsSeconds = uvar->maxEndTime;
  constraints.startTime = &minStartTimeGPS;
  constraints.endTime = &maxEndTimeGPS;

  /* get full SFT-catalog of all matching (multi-IFO) SFTs */
  LogPrintf (LOG_DEBUG, "Finding all SFTs to load ... ");
  TRY ( LALSFTdataFind ( status->statusPtr, &catalog, uvar->DataFiles, &constraints ), status);
  LogPrintfVerbatim (LOG_DEBUG, "done. (found %d SFTs)\n", catalog->length);

  if ( constraints.detector )
    LALFree ( constraints.detector );

  if ( !catalog || catalog->length == 0 )
    {
      XLALPrintError ("\nSorry, didn't find any matching SFTs with pattern '%s'!\n\n", uvar->DataFiles );
      ABORT ( status,  COMPUTEFSTATISTIC_EINPUT,  COMPUTEFSTATISTIC_MSGEINPUT);
    }

  /* deduce start- and end-time of the observation spanned by the data */
  UINT4 numSFTfiles = catalog->length;
  cfg->Tsft = 1.0 / catalog->data[0].header.deltaF;
  cfg->startTime = catalog->data[0].header.epoch;
  endTime   = catalog->data[numSFTfiles - 1].header.epoch;
  XLALGPSAdd(&endTime, cfg->Tsft);	/* add on Tsft to last SFT start-time */

  // time spanned by the SFTs
  cfg->Tspan = XLALGPSDiff ( &endTime, &cfg->startTime );

  { /* ----- load ephemeris-data ----- */
    cfg->ephemeris = XLALInitBarycenter( uvar->ephemEarth, uvar->ephemSun );
    if ( !cfg->ephemeris ) {
      XLALPrintError("XLALInitBarycenter failed: could not load Earth ephemeris '%s' and Sun ephemeris '%s'\n", uvar->ephemEarth, uvar->ephemSun);
      ABORT ( status,  COMPUTEFSTATISTIC_EINPUT,  COMPUTEFSTATISTIC_MSGEINPUT);
    }
  }

  /* ----- get reference-times (from user if given, use startTime otherwise): ----- */
  LIGOTimeGPS refTime;
  if ( LALUserVarWasSet(&uvar->refTime) )
    {
      XLALGPSSetREAL8 ( &refTime, uvar->refTime );
    }
  else if (LALUserVarWasSet(&uvar->refTimeMJD))
    {
      /* convert MJD peripase to GPS using Matt Pitkins code found at lal/packages/pulsar/src/BinaryPulsarTimeing.c */
      REAL8 GPSfloat;
      GPSfloat = XLALTDBMJDtoGPS(uvar->refTimeMJD);
      XLALGPSSetREAL8 ( &refTime, GPSfloat );
    }
  else
    refTime = cfg->startTime;

  /* define sky position variables from user input */
  if (LALUserVarWasSet(&uvar->RA))
    {
      /* use Matt Pitkins conversion code found in lal/packages/pulsar/src/BinaryPulsarTiming.c */
      cfg->Alpha = XLALhmsToRads(uvar->RA);
    }
  else cfg->Alpha = uvar->Alpha;
  if (LALUserVarWasSet(&uvar->Dec))
    {
      /* use Matt Pitkins conversion code found in lal/packages/pulsar/src/BinaryPulsarTiming.c */
      cfg->Delta = XLALdmsToRads(uvar->Dec);
    }
  else cfg->Delta = uvar->Delta;


  REAL8 fMin, fMax, f1dotMin, f1dotMax, f2dotMin, f2dotMax, f3dotMin, f3dotMax;
  { /* ----- prepare spin-range at refTime (in *canonical format*, ie all Bands >= 0) ----- */
    fMin = MYMIN ( uvar->Freq, uvar->Freq + uvar->FreqBand );
    fMax = MYMAX ( uvar->Freq, uvar->Freq + uvar->FreqBand );

    /* Specific to --gridType=GRID_SPINDOWN_AGEBRK parameter space */
    if (uvar->gridType == GRID_SPINDOWN_AGEBRK) {

      /* Set the first and second spindown ranges to the maximum that will be
       * encountered by the age-braking index parameter space. These will *ONLY*
       * be used by older code to, e.g. load the correct band of SFTs, and
       * will *NOT* be used by the tiling code itself.
       * The formulas used here are, with age=a, min braking=n, max braking=N:
       *
       * -f0/((n-1)*a) <= f1 <= -f0/((N-1)*a)
       * n*(f1^2)/f0 <= f2 <= N*(f1^2)/f0
       *
       * f0/f1 are taken as the maximum/minimum value appropriate for getting the
       * maximum/minimum values of f1/f2 (note that f1 is strictly negative here).
       */

      f1dotMin = -1.0 * fMax / ((uvar->minBraking - 1.0) * uvar->spindownAge);
      f1dotMax = -1.0 * fMin / ((uvar->maxBraking - 1.0) * uvar->spindownAge);

      f2dotMin = uvar->minBraking * (f1dotMax * f1dotMax) / fMax;
      f2dotMax = uvar->maxBraking * (f1dotMin * f1dotMin) / fMin;

    }
    else {     /* Used for all other --gridTypes */

      f1dotMin = MYMIN ( uvar->f1dot, uvar->f1dot + uvar->f1dotBand );
      f1dotMax = MYMAX ( uvar->f1dot, uvar->f1dot + uvar->f1dotBand );

      f2dotMin = MYMIN ( uvar->f2dot, uvar->f2dot + uvar->f2dotBand );
      f2dotMax = MYMAX ( uvar->f2dot, uvar->f2dot + uvar->f2dotBand );

    }

    f3dotMin = MYMIN ( uvar->f3dot, uvar->f3dot + uvar->f3dotBand );
    f3dotMax = MYMAX ( uvar->f3dot, uvar->f3dot + uvar->f3dotBand );

  } /* spin-range at refTime */


  { /* ----- set up Doppler region to scan ----- */
    BOOLEAN haveAlphaDelta = (LALUserVarWasSet(&uvar->Alpha) && LALUserVarWasSet(&uvar->Delta)) || (LALUserVarWasSet(&uvar->RA) && LALUserVarWasSet(&uvar->Dec));

    if (uvar->skyRegion)
      {
	cfg->searchRegion.skyRegionString = (CHAR*)LALCalloc(1, strlen(uvar->skyRegion)+1);
	if ( cfg->searchRegion.skyRegionString == NULL ) {
	  ABORT (status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM);
	}
	strcpy (cfg->searchRegion.skyRegionString, uvar->skyRegion);
      }
    else if (haveAlphaDelta)    /* parse this into a sky-region */
      {
	TRY ( SkySquare2String( status->statusPtr, &(cfg->searchRegion.skyRegionString),
				cfg->Alpha, cfg->Delta,	uvar->AlphaBand, uvar->DeltaBand), status);
      }

    /* spin searchRegion defined by spin-range at reference-time */
    cfg->searchRegion.refTime = refTime;

    cfg->searchRegion.fkdot[0] = fMin;
    cfg->searchRegion.fkdot[1] = f1dotMin;
    cfg->searchRegion.fkdot[2] = f2dotMin;
    cfg->searchRegion.fkdot[3] = f3dotMin;

    cfg->searchRegion.fkdotBand[0] = fMax - fMin;
    cfg->searchRegion.fkdotBand[1] = f1dotMax - f1dotMin;
    cfg->searchRegion.fkdotBand[2] = f2dotMax - f2dotMin;
    cfg->searchRegion.fkdotBand[3] = f3dotMax - f3dotMin;

  } /* get DopplerRegion */

  /* ----- set fixed grid step-sizes from user-input: only used for GRID_FLAT ----- */
  cfg->stepSizes.Alpha = uvar->dAlpha;
  cfg->stepSizes.Delta = uvar->dDelta;
  cfg->stepSizes.fkdot[0] = uvar->dFreq;
  cfg->stepSizes.fkdot[1] = uvar->df1dot;
  cfg->stepSizes.fkdot[2] = uvar->df2dot;
  cfg->stepSizes.fkdot[3] = uvar->df3dot;


  REAL8 tmpFreqBandRef = cfg->searchRegion.fkdotBand[0];

  /* initialize full multi-dimensional Doppler-scanner */
  {
    DopplerFullScanInit scanInit;			/* init-structure for DopperScanner */

    scanInit.searchRegion = cfg->searchRegion;
    if ( uvar->useResamp )	// in the resampling-case, temporarily take out frequency-dimension of the Doppler template bank
      scanInit.searchRegion.fkdotBand[0] = 0;

    scanInit.gridType = uvar->gridType;
    scanInit.gridFile = uvar->gridFile;
    scanInit.metricType = uvar->metricType;
    scanInit.projectMetric = uvar->projectMetric;
    scanInit.metricMismatch = uvar->metricMismatch;
    scanInit.stepSizes = cfg->stepSizes;
    scanInit.ephemeris = cfg->ephemeris;		/* used by Ephemeris-based metric */
    scanInit.startTime = cfg->startTime;
    scanInit.Tspan     = cfg->Tspan;

    // just use first SFTs' IFO for metric (should be irrelevant)
    LALDetector *detector;
    if ( ( detector = XLALGetSiteInfo ( catalog->data[0].header.name ) ) == NULL ) {
      LogPrintf ( LOG_CRITICAL, "\nXLALGetSiteInfo() failed for detector '%s'\n", catalog->data[0].header.name );
      ABORT ( status, COMPUTEFSTATISTIC_EXLAL, COMPUTEFSTATISTIC_MSGEXLAL );
    }
    scanInit.Detector  = detector;

    /* Specific to --gridType=GRID_SPINDOWN_AGEBRK parameter space */
    if (uvar->gridType == GRID_SPINDOWN_AGEBRK) {
      scanInit.extraArgs[0] = uvar->spindownAge;
      scanInit.extraArgs[1] = uvar->minBraking;
      scanInit.extraArgs[2] = uvar->maxBraking;
    }

    LogPrintf (LOG_DEBUG, "Setting up template grid ... ");
    TRY ( InitDopplerFullScan ( status->statusPtr, &cfg->scanState, &scanInit), status);
    LogPrintfVerbatim (LOG_DEBUG, "template grid ready.\n");
    XLALNumDopplerTemplates ( cfg->scanState );
    XLALFree ( detector );
  }


  { /* ----- What frequency-band do we need to read from the SFTs?
     * propagate spin-range from refTime to startTime and endTime of observation
     */
    PulsarSpinRange spinRangeRef, spinRangeStart, spinRangeEnd;	/* temporary only */
    REAL8 fmaxStart, fmaxEnd, fminStart, fminEnd;

    // extract spanned spin-range at reference-time from the template-bank
    if ( XLALGetDopplerSpinRange ( &spinRangeRef, cfg->scanState ) != XLAL_SUCCESS ) {
      LogPrintf ( LOG_CRITICAL, "\nXLALGetDopplerSpinRange() failed\n" );
      ABORT ( status, COMPUTEFSTATISTIC_EXLAL, COMPUTEFSTATISTIC_MSGEXLAL );
    }

    // in the resampling case, we need to restore the frequency-band info now, which we set to 0
    // before calling the DopplerInit template bank construction
    if ( uvar->useResamp )
      spinRangeRef.fkdotBand[0] = tmpFreqBandRef;

    // write this search spin-range@refTime back into 'cfg' struct for users' reference
    cfg->searchRegion.refTime = spinRangeRef.refTime;
    memcpy ( cfg->searchRegion.fkdot, spinRangeRef.fkdot, sizeof(cfg->searchRegion.fkdot) );
    memcpy ( cfg->searchRegion.fkdotBand, spinRangeRef.fkdotBand, sizeof(cfg->searchRegion.fkdotBand) );

    /* compute spin-range at startTime of observation */
    TRY ( LALExtrapolatePulsarSpinRange (status->statusPtr, &spinRangeStart, cfg->startTime, &spinRangeRef ), status );
    /* compute spin-range at endTime of these SFTs */
    TRY ( LALExtrapolatePulsarSpinRange (status->statusPtr, &spinRangeEnd, endTime, &spinRangeStart ), status );

    fminStart = spinRangeStart.fkdot[0];
    /* ranges are in canonical format! */
    fmaxStart = fminStart + spinRangeStart.fkdotBand[0];
    fminEnd   = spinRangeEnd.fkdot[0];
    fmaxEnd   = fminEnd + spinRangeEnd.fkdotBand[0];

    /*  get covering frequency-band  */
    fCoverMax = MYMAX ( fmaxStart, fmaxEnd );
    fCoverMin = MYMIN ( fminStart, fminEnd );

    /* correct for doppler-shift */
    fCoverMax *= 1.0 + uvar->dopplermax;
    fCoverMin *= 1.0 - uvar->dopplermax;

  } /* extrapolate spin-range */

  /* ----- set computational parameters for F-statistic from User-input ----- */
  if ( uvar->useResamp ) {	// use resampling

    cfg->Fstat_in = XLALCreateFstatInput_Resamp();
    if ( cfg->Fstat_in == NULL ) {
      XLALPrintError("%s: XLALCreateFstatInput_Resamp() failed with errno=%d", __func__, xlalErrno);
      ABORT ( status, COMPUTEFSTATISTIC_EXLAL, COMPUTEFSTATISTIC_MSGEXLAL );
    }

  } else {			// use demodulation

    cfg->Fstat_in = XLALCreateFstatInput_Demod( uvar->Dterms, DEMODHL_BEST );
    if ( cfg->Fstat_in == NULL ) {
      XLALPrintError("%s: XLALCreateFstatInput_Demod() failed with errno=%d", __func__, xlalErrno);
      ABORT ( status, COMPUTEFSTATISTIC_EXLAL, COMPUTEFSTATISTIC_MSGEXLAL );
    }

  }

  /* if single-only flag is given, assume a PSD with sqrt(S) = 1.0 */
  MultiNoiseFloor assumeSqrtSX, *p_assumeSqrtSX;
  if ( uvar->SignalOnly ) {
    assumeSqrtSX.length = XLALCountIFOsInCatalog(catalog);
    for (UINT4 X = 0; X < assumeSqrtSX.length; ++X) {
      assumeSqrtSX.sqrtSn[X] = 1.0;
    }
    p_assumeSqrtSX = &assumeSqrtSX;
  } else {
    p_assumeSqrtSX = NULL;
  }

  PulsarParamsVector *injectSources = NULL;
  MultiNoiseFloor *p_injectSqrtSX = NULL;
  if ( XLALSetupFstatInput( cfg->Fstat_in, catalog, fCoverMin, fCoverMax,
                                injectSources, p_injectSqrtSX, p_assumeSqrtSX,
                                uvar->RngMedWindow, cfg->ephemeris, uvar->SSBprecision, 0 ) != XLAL_SUCCESS ) {
    XLALPrintError("%s: XLALSetupFstatInput() failed with errno=%d", __func__, xlalErrno);
    ABORT ( status, COMPUTEFSTATISTIC_EXLAL, COMPUTEFSTATISTIC_MSGEXLAL );
  }

  XLALDestroySFTCatalog(catalog);

  cfg->Fstat_what = FSTATQ_2F;   // always calculate multi-detector 2F
  if ( LALUserVarWasSet( &uvar->outputLoudest ) ) {
    cfg->Fstat_what |= FSTATQ_FAFB;   // also calculate Fa,b parts for parameter estimation
  }

  /* get SFT detectors and timestamps */
  const MultiLALDetector *multiIFO = XLALGetFstatInputDetectors( cfg->Fstat_in );
  if ( multiIFO == NULL ) {
    XLALPrintError("%s: XLALGetFstatInputDetectors() failed with errno=%d", __func__, xlalErrno);
    ABORT ( status, COMPUTEFSTATISTIC_EXLAL, COMPUTEFSTATISTIC_MSGEXLAL );
  }
  const MultiLIGOTimeGPSVector *multiTS = XLALGetFstatInputTimestamps( cfg->Fstat_in );
  if ( multiTS == NULL ) {
    XLALPrintError("%s: XLALGetFstatInputTimestamps() failed with errno=%d", __func__, xlalErrno);
    ABORT ( status, COMPUTEFSTATISTIC_EXLAL, COMPUTEFSTATISTIC_MSGEXLAL );
  }

  /* count total number of SFTs loaded */
  cfg->NSFTs = 0;
  for ( UINT4 X = 0; X < multiTS->length; X++ ) {
    cfg->NSFTs += multiTS->data[X]->length;
  }

  /* for column headings string, get number of detectors, detector name vector, and SFTs per detector vector */
  {
    const UINT4 numDetectors = multiIFO->length;
    cfg->numSFTsPerDet = XLALCreateUINT4Vector( numDetectors );
    if ( cfg->numSFTsPerDet == NULL ) {
      XLALPrintError ("%s: XLALCreateUINT4Vector( %u ) failed with errno=%d\n", __func__, numDetectors, xlalErrno );
      ABORT ( status, COMPUTEFSTATISTIC_EXLAL, COMPUTEFSTATISTIC_MSGEXLAL );
    }
    cfg->detectorIDs = NULL;
    for (UINT4 X = 0; X < numDetectors; X++) {
      cfg->numSFTsPerDet->data[X] = multiTS->data[X]->length;
      if ( (cfg->detectorIDs = XLALAppendString2Vector ( cfg->detectorIDs, multiIFO->sites[X].frDetector.prefix )) == NULL ) {
        XLALPrintError ("%s: XLALAppendString2Vector() failed with errno=%d\n", __func__, xlalErrno );
        ABORT ( status, COMPUTEFSTATISTIC_EXLAL, COMPUTEFSTATISTIC_MSGEXLAL );
      }
    } /* for X < numDetectors */
  }

  /* internal refTime is used for computing the F-statistic at, to avoid large (t - tRef)^2 values */
  if ( LALUserVarWasSet ( &uvar->internalRefTime ) ) {
    XLALGPSSetREAL8 ( &(cfg->internalRefTime), uvar->internalRefTime);
  }
  else
    {
      LIGOTimeGPS midTime = cfg->startTime;
      XLALGPSAdd ( &midTime, 0.5 * XLALGPSDiff( &endTime, &cfg->startTime ) );	// mid-time of observation
      cfg->internalRefTime = midTime;
    }

  /* ----- set up scanline-window if requested for 1D local-maximum clustering on scanline ----- */
  if ( (cfg->scanlineWindow = XLALCreateScanlineWindow ( uvar->clusterOnScanline )) == NULL ) {
    ABORT (status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM );
  }

  /* set number of toplist candidates from fraction if asked to */
  if (0.0 < uvar->FracCandidatesToKeep && uvar->FracCandidatesToKeep <= 1.0) {
    if (XLALNumDopplerTemplates(cfg->scanState) <= 0.0) {
      LogPrintf(LOG_CRITICAL, "Cannot use FracCandidatesToKeep because number of templates was counted to be zero!\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
    toplist_length = ceil(XLALNumDopplerTemplates(cfg->scanState) * uvar->FracCandidatesToKeep);
  }

  /* ----- set up toplist if requested ----- */
  if ( toplist_length > 0 ) {
    if ( strcmp(uvar->RankingStatistic, "F") == 0 )
     cfg->RankingStatistic = RANKBY_2F;
    else if ( strcmp(uvar->RankingStatistic, "LV") == 0 )
      {
        if ( !uvar->computeLV ) {
          XLALPrintError ("\nERROR: Ranking by LV-stat only possible if --computeLV given.\n\n");
          ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT );
        }
        cfg->RankingStatistic = RANKBY_LV;
      }
    else
      {
        XLALPrintError ("\nERROR: Invalid value specified for candidate ranking - supported are 'F' and 'LV'.\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT );
      }

    if ( cfg->RankingStatistic == RANKBY_LV )
      {
        if ( create_toplist( &(cfg->FstatToplist), toplist_length, sizeof(FstatCandidate), compareFstatCandidates_LV) != 0 )
          ABORT (status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM );
      }
    else // rank by F-stat
      {
        if ( create_toplist( &(cfg->FstatToplist), toplist_length, sizeof(FstatCandidate), compareFstatCandidates) != 0 )
          ABORT (status, COMPUTEFSTATISTIC_EMEM, COMPUTEFSTATISTIC_MSGEMEM );
      }
  } /* if toplist_length > 0 */


  /* ----- transient-window related parameters ----- */
  int twtype;
  if ( (twtype = XLALParseTransientWindowName ( uvar->transient_WindowType )) < 0 ) {
    ABORT (status, COMPUTEFSTATISTIC_EXLAL, COMPUTEFSTATISTIC_MSGEXLAL );
  }
  cfg->transientWindowRange.type = twtype;

  /* make sure user doesn't set window=none but sets window-parameters => indicates she didn't mean 'none' */
  if ( cfg->transientWindowRange.type == TRANSIENT_NONE )
    if ( XLALUserVarWasSet ( &uvar->transient_t0Days ) || XLALUserVarWasSet ( &uvar->transient_t0DaysBand ) || XLALUserVarWasSet ( &uvar->transient_dt0 ) ||
         XLALUserVarWasSet ( &uvar->transient_tauDays ) || XLALUserVarWasSet ( &uvar->transient_tauDaysBand ) || XLALUserVarWasSet ( &uvar->transient_dtau ) ) {
      XLALPrintError ("%s: ERROR: transientWindow->type == NONE, but window-parameters were set! Use a different window-type!\n", __func__ );
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }

  if (   uvar->transient_t0DaysBand < 0 || uvar->transient_tauDaysBand < 0 ) {
    XLALPrintError ("%s: only positive t0/tau bands allowed (%f, %f)\n", __func__, uvar->transient_t0DaysBand, uvar->transient_tauDaysBand );
    ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
  }

  cfg->transientWindowRange.t0      = cfg->startTime.gpsSeconds + uvar->transient_t0Days * DAY24;
  cfg->transientWindowRange.t0Band  = uvar->transient_t0DaysBand * DAY24;


  if ( XLALUserVarWasSet ( &uvar->transient_dt0 ) )
    cfg->transientWindowRange.dt0 = uvar->transient_dt0;
  else
    cfg->transientWindowRange.dt0 = cfg->Tsft;

  cfg->transientWindowRange.tau     = (UINT4) ( uvar->transient_tauDays * DAY24 );
  cfg->transientWindowRange.tauBand = (UINT4) ( uvar->transient_tauDaysBand * DAY24 );

  if ( XLALUserVarWasSet ( &uvar->transient_dtau ) )
    cfg->transientWindowRange.dtau = uvar->transient_dtau;
  else
    cfg->transientWindowRange.dtau = cfg->Tsft;


  /* get atoms back from Fstat-computing, either if atoms-output or transient-Bstat output was requested */
  if ( ( uvar->outputFstatAtoms != NULL ) || ( uvar->outputTransientStats != NULL ) ) {
    cfg->Fstat_what |= FSTATQ_ATOMS_PER_DET;
  }

  /* return single-IFO Fstat values for Line-veto statistic */
  if ( uvar->outputSingleFstats || uvar->computeLV ) {
    cfg->Fstat_what |= FSTATQ_2F_PER_DET;
  }

  /* ---------- prepare Line Veto statistics parameters ---------- */
  if ( uvar->LVrho < 0.0 ) {
    XLALPrintError("Invalid LV prior rho (given rho=%f, need rho>=0)!\n", uvar->LVrho);
    ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
  }
  else if ( uvar->LVrho > 0.0 )
    cfg->LVlogRhoTerm = 4.0 * log(uvar->LVrho) - log(70.0);
  else /* if uvar.LVrho == 0.0, logRhoTerm should become irrelevant in summation */
    cfg->LVlogRhoTerm = - LAL_REAL8_MAX;

  if ( uvar->computeLV && uvar->LVlX )
    {
      const UINT4 numDetectors = multiIFO->length;
      if (  uvar->LVlX->length != numDetectors ) {
        XLALPrintError( "Length of LV prior ratio vector does not match number of detectors! (%d != %d)\n", uvar->LVlX->length, numDetectors);
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }
      if ( (cfg->LVloglX = XLALCreateREAL8Vector ( numDetectors )) == NULL ) {
        XLALPrintError ("Failed to XLALCreateREAL8Vector ( %d )\n", numDetectors );
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }
      for (UINT4 X = 0; X < numDetectors; X++)
        {
          REAL4 LVlX;
          if ( 1 != sscanf ( uvar->LVlX->data[X], "%" LAL_REAL4_FORMAT, &LVlX ) ) {
            XLALPrintError ( "Illegal REAL4 commandline argument to --LVlX[%d]: '%s'\n", X, uvar->LVlX->data[X]);
            ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
          }
          if ( LVlX < 0.0 ) {
            XLALPrintError ( "Negative input prior-ratio for detector X=%d lX[X]=%f\n", X, LVlX );
            ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
          }
          else if ( LVlX > 0.0 )
            cfg->LVloglX->data[X] = log ( LVlX );
          else /* if zero prior ratio, approximate log(0)=-inf by -LAL_REA4_MAX to avoid raising underflow exceptions */
            cfg->LVloglX->data[X] = - LAL_REAL8_MAX;
        } /* for X < numDetectors */
    } /* if ( uvar.computeLV && uvar.LVlX ) */

  // ----- check that resampling option was used sensibly ...
  if ( uvar->useResamp )
    {
      // FIXME: probably should check a few more things, can't think of any right now ...
      // let's hope users are sensible
    }

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* InitFStat() */

/**
 * Produce a log-string describing the present run-setup
 */
CHAR *
XLALGetLogString ( const ConfigVariables *cfg )
{
  XLAL_CHECK_NULL ( cfg != NULL, XLAL_EINVAL );

  CHAR *logstr = NULL;
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, "%% cmdline: " )) != NULL, XLAL_EFUNC );
  CHAR *cmdline;
  XLAL_CHECK_NULL ( (cmdline = XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, cmdline )) != NULL, XLAL_EFUNC );
  XLALFree ( cmdline );
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, "\n" )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, cfg->VCSInfoString )) != NULL, XLAL_EFUNC );

  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, "%% Started search: " )) != NULL, XLAL_EFUNC );
  time_t tp = time(NULL);
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, asctime( gmtime( &tp ) ))) != NULL, XLAL_EFUNC );
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, "%% Loaded SFTs: [ " )) != NULL, XLAL_EFUNC );

#define BUFLEN 1024
  CHAR buf[BUFLEN];

  UINT4 numDet = cfg->detectorIDs->length;
  for ( UINT4 X=0; X < numDet; X ++ )
    {
      XLAL_CHECK_NULL ( snprintf ( buf, BUFLEN, "%s:%d%s",  cfg->detectorIDs->data[X], cfg->numSFTsPerDet->data[X], (X < numDet - 1) ? ", " : " ]\n" ) < BUFLEN, XLAL_EBADLEN );
      XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, buf )) != NULL, XLAL_EFUNC );
    }
  INT4 startTimeSeconds = cfg->startTime.gpsSeconds;
  struct tm startTimeUTC = *XLALGPSToUTC ( &startTimeUTC, startTimeSeconds );
  {
    CHAR *startTimeUTCString = XLALStringDuplicate ( asctime(&startTimeUTC) );
    startTimeUTCString[strlen(startTimeUTCString)-2] = 0;	// kill trailing newline
    XLAL_CHECK_NULL ( snprintf ( buf, BUFLEN, "%%%% GPS starttime         = %d (%s GMT)\n", startTimeSeconds, startTimeUTCString ) < BUFLEN, XLAL_EBADLEN );
    XLALFree ( startTimeUTCString );
  }
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, buf )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_NULL ( snprintf ( buf, BUFLEN, "%%%% Total time spanned    = %.0f s (%.2f hours)\n", cfg->Tspan, cfg->Tspan/3600.0 ) < BUFLEN, XLAL_EBADLEN );
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, buf )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_NULL ( snprintf (buf, BUFLEN, "%%%% InternalRefTime       = %.16g \n", XLALGPSGetREAL8 ( &(cfg->internalRefTime)) ) < BUFLEN, XLAL_EBADLEN );
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, buf )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_NULL ( snprintf (buf, BUFLEN, "%%%% Pulsar-params refTime = %.16g \n", XLALGPSGetREAL8 ( &(cfg->searchRegion.refTime) )) < BUFLEN, XLAL_EBADLEN );
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, buf )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, "%% Spin-range at refTime: fkdot = [ " )) != NULL, XLAL_EFUNC );
  for (UINT4 k=0; k < PULSAR_MAX_SPINS; k ++ )
    {
      XLAL_CHECK_NULL ( snprintf(buf, BUFLEN, "%.16g:%.16g%s", cfg->searchRegion.fkdot[k], cfg->searchRegion.fkdot[k] + cfg->searchRegion.fkdotBand[k],
                                 (k < PULSAR_MAX_SPINS - 1)?", ":" ]\n") < BUFLEN, XLAL_EBADLEN );
      XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, buf )) != NULL, XLAL_EFUNC );
    }

  /* return result */
  return logstr;

} // XLALGetLogString()



/***********************************************************************/
/**
 * Log the all relevant parameters of the present search-run to a log-file.
 * The name of the log-file is log_fname
 * <em>NOTE:</em> Currently this function only logs the user-input and code-versions.
 */
void
WriteFStatLog ( LALStatus *status, const CHAR *log_fname, const CHAR *log_string )
{
  FILE *fplog;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  if ( !log_fname || !log_string ) {	/* no logfile given */
    ABORT (status, COMPUTEFSTATISTIC_ENULL, COMPUTEFSTATISTIC_MSGENULL);
  }

  /* prepare log-file for writing */
  if ( (fplog = fopen(log_fname, "wb" )) == NULL) {
    LogPrintf ( LOG_CRITICAL , "Failed to open log-file '%s' for writing.\n\n", log_fname );
    ABORT (status, COMPUTEFSTATISTIC_ESYS, COMPUTEFSTATISTIC_MSGESYS);
  }

  fprintf (fplog, "%%%% LOG-FILE for ComputeFStatistic run\n\n");
  fprintf (fplog, "%s", log_string);
  fclose (fplog);


  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* WriteFStatLog() */


/** Free all globally allocated memory. */
void
Freemem(LALStatus *status,  ConfigVariables *cfg)
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  XLALDestroyUINT4Vector ( cfg->numSFTsPerDet );
  XLALDestroyStringVector ( cfg->detectorIDs );

  XLALDestroyFstatInput ( cfg->Fstat_in );

  /* destroy FstatToplist if any */
  if ( cfg->FstatToplist )
    free_toplist( &(cfg->FstatToplist) );

  if ( cfg->scanlineWindow )
    XLALDestroyScanlineWindow ( cfg->scanlineWindow );

  /* Free config-Variables and userInput stuff */
  TRY (LALDestroyUserVars (status->statusPtr), status);

  if ( cfg->searchRegion.skyRegionString )
    LALFree ( cfg->searchRegion.skyRegionString );

  /* Free ephemeris data */
  XLALDestroyEphemerisData ( cfg->ephemeris );

  if ( cfg->VCSInfoString )
    XLALFree ( cfg->VCSInfoString );
  if ( cfg->logstring )
    LALFree ( cfg->logstring );

  if ( cfg->LVloglX )
    XLALDestroyREAL8Vector ( cfg->LVloglX );

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* Freemem() */


/*----------------------------------------------------------------------*/
/**
 * Some general consistency-checks on user-input.
 * Throws an error plus prints error-message if problems are found.
 */
void
checkUserInputConsistency (LALStatus *status, const UserInput_t *uvar)
{
  INITSTATUS(status);

  /* check that only alpha OR RA has been set */
  if ( LALUserVarWasSet(&uvar->Alpha) && (LALUserVarWasSet(&uvar->RA)) )
    {
      XLALPrintError ("\nInput either Alpha OR RA, not both!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
  /* check that only delta OR Dec has been set */
  if ( LALUserVarWasSet(&uvar->Delta) && (LALUserVarWasSet(&uvar->Dec)) )
    {
      XLALPrintError ("\nInput either Delta OR Dec, not both!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }

  /* check for negative stepsizes in Freq, Alpha, Delta */
  if ( LALUserVarWasSet(&uvar->dAlpha) && (uvar->dAlpha < 0) )
    {
      XLALPrintError ("\nNegative value of stepsize dAlpha not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar->dDelta) && (uvar->dDelta < 0) )
    {
      XLALPrintError ("\nNegative value of stepsize dDelta not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar->dFreq) && (uvar->dFreq < 0) )
    {
      XLALPrintError ("\nNegative value of stepsize dFreq not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }

  /* check that reference time has not been set twice */
  if ( LALUserVarWasSet(&uvar->refTime) && LALUserVarWasSet(&uvar->refTimeMJD) )
    {
      XLALPrintError ("\nSet only uvar->refTime OR uvar->refTimeMJD OR leave empty to use SSB start time as Tref!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }

  /* binary parameter checks */
  if ( LALUserVarWasSet(&uvar->orbitPeriod) && (uvar->orbitPeriod <= 0) )
    {
      XLALPrintError ("\nNegative or zero value of orbital period not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar->orbitasini) && (uvar->orbitasini < 0) )
    {
      XLALPrintError ("\nNegative value of projected orbital semi-major axis not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
   if ( LALUserVarWasSet(&uvar->orbitTpSSBMJD) && (LALUserVarWasSet(&uvar->orbitTpSSBsec) || LALUserVarWasSet(&uvar->orbitTpSSBnan)))
    {
      XLALPrintError ("\nSet only uvar->orbitTpSSBMJD OR uvar->orbitTpSSBsec/nan to specify periapse passage time!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
   if ( LALUserVarWasSet(&uvar->orbitTpSSBMJD) && (uvar->orbitTpSSBMJD < 0) )
    {
      XLALPrintError ("\nNegative value of the true time of orbital periapsis not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar->orbitTpSSBsec) && (uvar->orbitTpSSBsec < 0) )
    {
      XLALPrintError ("\nNegative value of seconds part of the true time of orbital periapsis not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar->orbitTpSSBnan) && ((uvar->orbitTpSSBnan < 0) || (uvar->orbitTpSSBnan >= 1e9)) )
    {
      XLALPrintError ("\nTime of nanoseconds part the true time of orbital periapsis must lie in range (0, 1e9]!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar->orbitArgp) && ((uvar->orbitArgp < 0) || (uvar->orbitArgp >= LAL_TWOPI)) )
    {
      XLALPrintError ("\nOrbital argument of periapse must lie in range [0 2*PI)!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar->orbitEcc) && (uvar->orbitEcc < 0) )
    {
      XLALPrintError ("\nNegative value of orbital eccentricity not allowed!\n\n");
      ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
    }

  /* grid-related checks */
  {
    BOOLEAN haveAlphaBand = LALUserVarWasSet( &uvar->AlphaBand );
    BOOLEAN haveDeltaBand = LALUserVarWasSet( &uvar->DeltaBand );
    BOOLEAN haveSkyRegion, haveAlphaDelta, haveGridFile;
    BOOLEAN useSkyGridFile, useFullGridFile, haveMetric, useMetric;

    haveSkyRegion  	= (uvar->skyRegion != NULL);
    haveAlphaDelta 	= (LALUserVarWasSet(&uvar->Alpha) && LALUserVarWasSet(&uvar->Delta) ) || (LALUserVarWasSet(&uvar->RA) && LALUserVarWasSet(&uvar->Dec) );
    haveGridFile      	= (uvar->gridFile != NULL);
    useSkyGridFile   	= (uvar->gridType == GRID_FILE_SKYGRID);
    useFullGridFile	= (uvar->gridType == GRID_FILE_FULLGRID);
    haveMetric     	= (uvar->metricType > LAL_PMETRIC_NONE);
    useMetric     	= (uvar->gridType == GRID_METRIC);

    if ( !useFullGridFile && !useSkyGridFile && haveGridFile )
      {
        XLALPrintError ("\nERROR: gridFile was specified but not needed for gridType=%d\n\n", uvar->gridType );
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }
    if ( useSkyGridFile && !haveGridFile )
      {
        XLALPrintError ("\nERROR: gridType=SKY-FILE, but no --gridFile specified!\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }
    if ( useFullGridFile && !haveGridFile )
      {
	XLALPrintError ("\nERROR: gridType=GRID-FILE, but no --gridFile specified!\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }

    if ( (haveAlphaBand && !haveDeltaBand) || (haveDeltaBand && !haveAlphaBand) )
      {
	XLALPrintError ("\nERROR: Need either BOTH (AlphaBand, DeltaBand) or NONE.\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }

    if ( haveSkyRegion && haveAlphaDelta )
      {
        XLALPrintError ("\nOverdetermined sky-region: only use EITHER (Alpha,Delta) OR skyRegion!\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }
    if ( !useMetric && haveMetric)
      {
        LALWarning (status, "\nWARNING: Metric was specified for non-metric grid... will be ignored!\n");
      }
    if ( useMetric && !haveMetric)
      {
        XLALPrintError ("\nERROR: metric grid-type selected, but no metricType selected\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }

    /* Specific checks for --gridType=GRID_SPINDOWN_{SQUARE,AGEBRK} parameter spaces */
    if (uvar->gridType == GRID_SPINDOWN_SQUARE || uvar->gridType == GRID_SPINDOWN_AGEBRK) {

      /* Check that no third spindown range were given */
      if (uvar->f3dot != 0.0 || uvar->f3dotBand != 0.0) {
        XLALPrintError ("\nERROR: f3dot and f3dotBand cannot be used with gridType={8,9}\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }

      /* Check that no grid spacings were given */
      if (uvar->df1dot != 0.0 || uvar->df2dot != 0.0 || uvar->df3dot != 0.0) {
        XLALPrintError ("\nERROR: df{1,2,3}dot cannot be used with gridType={8,9}\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }

    }

    /* Specific checks for --gridType=GRID_SPINDOWN_AGEBRK parameter space */
    if (uvar->gridType == GRID_SPINDOWN_AGEBRK) {

      /* Check age and braking indices */
      if (uvar->spindownAge <= 0.0) {
        XLALPrintError ("\nERROR: spindownAge must be strictly positive with gridType=9\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }
      if (uvar->minBraking <= 0.0) {
        XLALPrintError ("\nERROR: minBraking must be strictly positive with gridType=9\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }
      if (uvar->maxBraking <= 0.0) {
        XLALPrintError ("\nERROR: minBraking must be strictly positive with gridType=9\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }
      if (uvar->minBraking >= uvar->maxBraking) {
        XLALPrintError ("\nERROR: minBraking must be strictly less than maxBraking with gridType=9\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }

      /* Check that no first and second spindown ranges were given */
      if (uvar->f1dot != 0.0 || uvar->f1dotBand != 0.0) {
        XLALPrintError ("\nERROR: f1dot and f1dotBand cannot be used with gridType=9\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }
      if (uvar->f2dot != 0.0 || uvar->f2dotBand != 0.0) {
        XLALPrintError ("\nERROR: f2dot and f2dotBand cannot be used with gridType=9\n\n");
        ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
      }

    }

  } /* Grid-related checks */

  /* check NumCandidatesToKeep and FracCandidatesToKeep */
  if (LALUserVarWasSet(&uvar->NumCandidatesToKeep) && LALUserVarWasSet(&uvar->FracCandidatesToKeep)) {
    XLALPrintError ("\nERROR: NumCandidatesToKeep and FracCandidatesToKeep are mutually exclusive\n\n");
    ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
  }
  if (LALUserVarWasSet(&uvar->FracCandidatesToKeep) && (uvar->FracCandidatesToKeep <= 0.0 || 1.0 < uvar->FracCandidatesToKeep)) {
    XLALPrintError ("\nERROR: FracCandidatesToKeep must be greater than 0.0 and less than or equal to 1.0\n\n");
    ABORT (status, COMPUTEFSTATISTIC_EINPUT, COMPUTEFSTATISTIC_MSGEINPUT);
  }

  RETURN (status);
} /* checkUserInputConsistency() */

/* debug-output a(t) and b(t) into given file.
 * return 0 = OK, -1 on error
 */
int
outputBeamTS( const CHAR *fname, const AMCoeffs *amcoe, const DetectorStateSeries *detStates )
{
  FILE *fp;
  UINT4 i, len;

  if ( !fname || !amcoe || !amcoe->a || !amcoe->b || !detStates)
    return -1;

  len = amcoe->a->length;
  if ( (len != amcoe->b->length) || ( len != detStates->length ) )
    return -1;

  if ( (fp = fopen(fname, "wb")) == NULL )
    return -1;

  for (i=0; i < len; i ++ )
    {
      INT4 ret;
      ret = fprintf (fp, "%9d %f %f %f \n",
		     detStates->data[i].tGPS.gpsSeconds, detStates->data[i].LMST, amcoe->a->data[i], amcoe->b->data[i] );
      if ( ret < 0 )
	{
	  fprintf (fp, "ERROR\n");
	  fclose(fp);
	  return -1;
	}
    }

  fclose(fp);
  return 0;
} /* outputBeamTS() */

/**
 * write full 'PulsarCandidate' (i.e. Doppler params + Amplitude params + error-bars + Fa,Fb, F, + A,B,C,D
 * RETURN 0 = OK, -1 = ERROR
 */
int
write_PulsarCandidate_to_fp ( FILE *fp,  const PulsarCandidate *pulsarParams, const FstatCandidate *Fcand )
{
  if ( !fp || !pulsarParams || !Fcand  )
    return -1;

  fprintf (fp, "\n");

  fprintf (fp, "refTime  = % 9d;\n", pulsarParams->Doppler.refTime.gpsSeconds );   /* forget about ns... */

  fprintf (fp, "\n");

  /* Amplitude parameters with error-estimates */
  fprintf (fp, "h0       = % .6g;\n", pulsarParams->Amp.h0 );
  fprintf (fp, "dh0      = % .6g;\n", pulsarParams->dAmp.h0 );
  fprintf (fp, "cosi     = % .6g;\n", pulsarParams->Amp.cosi );
  fprintf (fp, "dcosi    = % .6g;\n", pulsarParams->dAmp.cosi );
  fprintf (fp, "phi0     = % .6g;\n", pulsarParams->Amp.phi0 );
  fprintf (fp, "dphi0    = % .6g;\n", pulsarParams->dAmp.phi0 );
  fprintf (fp, "psi      = % .6g;\n", pulsarParams->Amp.psi );
  fprintf (fp, "dpsi     = % .6g;\n", pulsarParams->dAmp.psi );

  fprintf (fp, "\n");

  /* Doppler parameters */
  fprintf (fp, "Alpha    = % .16g;\n", pulsarParams->Doppler.Alpha );
  fprintf (fp, "Delta    = % .16g;\n", pulsarParams->Doppler.Delta );
  fprintf (fp, "Freq     = % .16g;\n", pulsarParams->Doppler.fkdot[0] );
  fprintf (fp, "f1dot    = % .16g;\n", pulsarParams->Doppler.fkdot[1] );
  fprintf (fp, "f2dot    = % .16g;\n", pulsarParams->Doppler.fkdot[2] );
  fprintf (fp, "f3dot    = % .16g;\n", pulsarParams->Doppler.fkdot[3] );

  fprintf (fp, "\n");

  /* Binary parameters */
  if (pulsarParams->Doppler.asini > 0)
    {
      fprintf (fp, "orbitPeriod       = % .16g;\n", pulsarParams->Doppler.period );
      fprintf (fp, "orbitasini        = % .16g;\n", pulsarParams->Doppler.asini );
      fprintf (fp, "orbitTpSSBsec     = % .8d;\n", pulsarParams->Doppler.tp.gpsSeconds );
      fprintf (fp, "orbitTpSSBnan     = % .8d;\n", pulsarParams->Doppler.tp.gpsNanoSeconds );
      fprintf (fp, "orbitArgp         = % .16g;\n", pulsarParams->Doppler.argp );
      fprintf (fp, "orbitEcc          = % .16g;\n", pulsarParams->Doppler.ecc );
    }

  /* Amplitude Modulation Coefficients */
  fprintf (fp, "Ad       = % .6g;\n", Fcand->Mmunu.Ad );
  fprintf (fp, "Bd       = % .6g;\n", Fcand->Mmunu.Bd );
  fprintf (fp, "Cd       = % .6g;\n", Fcand->Mmunu.Cd );
  fprintf (fp, "Ed       = % .6g;\n", Fcand->Mmunu.Ed );
  fprintf (fp, "Sinv_Tsft= % .6g;\n", Fcand->Mmunu.Sinv_Tsft );
  fprintf (fp, "\n");

  /* Fstat-values */
  fprintf (fp, "Fa       = % .6g  %+.6gi;\n", creal(Fcand->Fa), cimag(Fcand->Fa) );
  fprintf (fp, "Fb       = % .6g  %+.6gi;\n", creal(Fcand->Fb), cimag(Fcand->Fb) );
  fprintf (fp, "twoF     = % .6g;\n", Fcand->twoF );
  /* single-IFO Fstat-values, if present */
  UINT4 X, numDet = Fcand->numDetectors;
  for ( X = 0; X < numDet ; X ++ )
    fprintf (fp, "twoF%d    = % .6g;\n", X, Fcand->twoFX[X] );
  /* LVstat */
  if ( !isnan(Fcand->LVstat) ) /* if --computeLV=FALSE, the LV field was initialised to NAN - do not output LV */
    fprintf (fp, "LV       = % .6g;\n", Fcand->LVstat );

  fprintf (fp, "\nAmpFisher = \\\n" );
  XLALfprintfGSLmatrix ( fp, "%.9g",pulsarParams->AmpFisherMatrix );

  return 0;

} /* write_PulsarCandidate_to_fp() */

/** comparison function for our candidates toplist */
int
compareFstatCandidates ( const void *candA, const void *candB )
{
  REAL8 twoF1 = ((const FstatCandidate *)candA)->twoF;
  REAL8 twoF2 = ((const FstatCandidate *)candB)->twoF;
  if ( twoF1 < twoF2 )
    return 1;
  else if ( twoF1 > twoF2 )
    return -1;
  else
    return 0;

} /* compareFstatCandidates() */

/** comparison function for our candidates toplist with alternate LV sorting statistic */
int
compareFstatCandidates_LV ( const void *candA, const void *candB )
{
  REAL8 LV1 = ((const FstatCandidate *)candA)->LVstat;
  REAL8 LV2 = ((const FstatCandidate *)candB)->LVstat;
  if ( LV1 < LV2 )
    return 1;
  else if ( LV1 > LV2 )
    return -1;
  else
    return 0;

} /* compareFstatCandidates_LV() */

/**
 * write one 'FstatCandidate' (i.e. only Doppler-params + Fstat) into file 'fp'.
 * Return: 0 = OK, -1 = ERROR
 */
int
write_FstatCandidate_to_fp ( FILE *fp, const FstatCandidate *thisFCand )
{

  if ( !fp || !thisFCand )
    return -1;

  /* add extra output-field containing per-detector FX if non-NULL */
  char extraStatsStr[256] = "";     /* defaults to empty */
  char buf0[256];
  /* LVstat */
  if ( !isnan(thisFCand->LVstat) ) /* if --computeLV=FALSE, the LV field was initialised to NAN - do not output LV */
      snprintf ( extraStatsStr, sizeof(extraStatsStr), " %.9g", thisFCand->LVstat );
  if ( thisFCand->numDetectors > 0 )
    {
      for ( UINT4 X = 0; X < thisFCand->numDetectors; X ++ )
        {
          snprintf ( buf0, sizeof(buf0), " %.9g", thisFCand->twoFX[X] );
          UINT4 len1 = strlen ( extraStatsStr ) + strlen ( buf0 ) + 1;
          if ( len1 > sizeof ( extraStatsStr ) ) {
            XLALPrintError ("%s: assembled output string too long! (%d > %d)\n", __func__, len1, sizeof(extraStatsStr ));
            break;      /* we can't really terminate with error in this function, but at least we avoid crashing */
          }
          strcat ( extraStatsStr, buf0 );
        } /* for X < numDet */
    } /* if FX */

  fprintf (fp, "%.16g %.16g %.16g %.16g %.16g %.16g %.9g%s\n",
	   thisFCand->doppler.fkdot[0], thisFCand->doppler.Alpha, thisFCand->doppler.Delta,
	   thisFCand->doppler.fkdot[1], thisFCand->doppler.fkdot[2], thisFCand->doppler.fkdot[3],
	   thisFCand->twoF, extraStatsStr );

  return 0;

} /* write_FstatCandidate_to_fp */

/* --------------------------------------------------------------------------------
 * Scanline window functions
 * FIXME: should go into a separate file once implementation is settled down ...
 *
 * --------------------------------------------------------------------------------*/

/**
 * Create a scanline window, with given windowWings >= 0.
 * Note: the actual window-size is 1 + 2 * windowWings
 */
scanlineWindow_t *
XLALCreateScanlineWindow ( UINT4 windowWings ) /**< number of neighbors on each side in scanlineWindow */
{
  scanlineWindow_t *ret = NULL;
  UINT4 windowLen = 1 + 2 * windowWings;

  if ( ( ret = LALCalloc ( 1, sizeof(*ret)) ) == NULL ) {
    XLAL_ERROR_NULL( COMPUTEFSTATISTIC_EMEM );
  }

  ret->length = windowLen;

  if ( (ret->window = LALCalloc ( windowLen, sizeof( ret->window[0] ) )) == NULL ) {
    LALFree ( ret );
    XLAL_ERROR_NULL( COMPUTEFSTATISTIC_EMEM );
  }

  ret->center = &(ret->window[ windowWings ]);	/* points to central bin */

  return ret;

} /* XLALCreateScanlineWindow() */

void
XLALDestroyScanlineWindow ( scanlineWindow_t *scanlineWindow )
{
  if ( !scanlineWindow )
    return;

  if ( scanlineWindow->window )
    LALFree ( scanlineWindow->window );

  LALFree ( scanlineWindow );

  return;

} /* XLALDestroyScanlineWindow() */

/**
 * Advance by pushing a new candidate into the scanline-window
 */
int
XLALAdvanceScanlineWindow ( const FstatCandidate *nextCand, scanlineWindow_t *scanWindow )
{
  UINT4 i;

  if ( !nextCand || !scanWindow || !scanWindow->window ) {
    XLAL_ERROR ( XLAL_EINVAL );
  }

  for ( i=1; i < scanWindow->length; i ++ )
    scanWindow->window[i - 1] = scanWindow->window[i];

  scanWindow->window[ scanWindow->length - 1 ] = *nextCand;	/* struct-copy */

  return XLAL_SUCCESS;

} /* XLALAdvanceScanlineWindow() */

/**
 * check wether central candidate in Scanline-window is a local maximum
 */
BOOLEAN
XLALCenterIsLocalMax ( const scanlineWindow_t *scanWindow, const UINT4 rankingStatistic )
{

  if ( !scanWindow || !scanWindow->center )
    return FALSE;

  if ( rankingStatistic == RANKBY_2F ) /* F statistic */
    {

      REAL8 twoF0 = scanWindow->center->twoF;

      for ( UINT4 i=0; i < scanWindow->length; i ++ )
        if ( scanWindow->window[i].twoF > twoF0 )
         return FALSE;

    }

  else if ( rankingStatistic == RANKBY_LV ) /* LV statistic */
    {

      REAL8 LV0 = scanWindow->center->LVstat;

      for ( UINT4 i=0; i < scanWindow->length; i ++ )
        if ( scanWindow->window[i].LVstat > LV0 )
         return FALSE;

    }
  else
    {
      XLALPrintError ("Unsupported ranking statistic '%d' ! Supported: 'F=0' and 'LV=2'.\n", rankingStatistic );
      return FALSE;
    }

  return TRUE;

} /* XLALCenterIsLocalMax() */

/**
 * Function to append one timing-info line to open output file.
 *
 * NOTE: called with NULL timing pointer writes header-comment line.
 */
int
write_TimingInfo_to_fp ( FILE * fp, const timingInfo_t *ti )
{
  /* input sanity */
  if ( !fp ) {
    XLALPrintError ("%s: invalid NULL input 'fp'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  /* if timingInfo == NULL ==> write header comment line */
  if ( ti == NULL )
    {
      fprintf ( fp, "%%%%NSFTs  costFstat[s]   tauMin[s]  tauMax[s]  NStart    NTau    costTransFstatMap[s]  costTransMarg[s] costTemplate[s]  tauF0 [s]\n");
      return XLAL_SUCCESS;
    } /* if ti == NULL */

  // compute fundamental timing constant 'tauF0' = F-stat time per template per SFT
  REAL8 tauF0 = ti->tauTemplate / ti->NSFTs;
  fprintf ( fp, "% 5d    %10.6e      %6d     %6d    %5d   %5d           %10.6e       %10.6e	%10.6e   %10.6e\n",
            ti->NSFTs, ti->tauFstat, ti->tauMin, ti->tauMax, ti->NStart, ti->NTau, ti->tauTransFstatMap, ti->tauTransMarg, ti->tauTemplate, tauF0 );

  return XLAL_SUCCESS;

} /* write_TimingInfo_to_fp() */

#ifdef HIGHRES_TIMING
/**
 * Return process User CPU time used.
 */
REAL8
XLALGetUserCPUTime ( void )
{
  struct timespec res;
  struct timespec ut;
  clockid_t clk_id = CLOCK_PROCESS_CPUTIME_ID;

  if ( clock_getres ( clk_id, &res ) != 0 ) {
    XLALPrintError ("%s: failed to call clock_getres(), errno = %d\n", __func__, errno );
    XLAL_ERROR_REAL8 ( XLAL_ESYS );
  }
  XLALPrintError ("%s: Clock-precision: {%ld s, %ld ns}\n", __func__, res.tv_sec, res.tv_nsec );

  if ( clock_gettime ( clk_id, &ut) != 0 ) {
    XLALPrintError ("%s: failed to call clock_gettime(), errno = %d\n", __func__, errno );
    XLAL_ERROR_REAL8 ( XLAL_ESYS );
  }

  return ut.tv_sec + (ut.tv_nsec/1.e9);

} /* XLALGetUserCPUTime() */
#endif
