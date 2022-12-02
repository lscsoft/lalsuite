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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

/*********************************************************************************/
/**
 * \defgroup lalpulsar_bin_Fstatistic Fstatistic Search Applications
 * \ingroup lalpulsar_bin_Apps
 */

/**
 * \author R. Prix, I. Gholami, Y. Ioth, Papa, X. Siemens, C. Messenger, K. Wette
 * \file
 * \ingroup lalpulsar_bin_Fstatistic
 * \brief
 * Calculate the F-statistic for a given parameter-space of pulsar GW signals.
 * Implements the so-called "F-statistic" as introduced in \cite JKS98 .
 *
 * This code is a descendant of an earlier implementation 'ComputeFStatistic.[ch]'
 * by Bruce Allen, Bernd Machenschalk, David Hammer, Jolien Creighton, Maria Alessandra Papa,
 * Reinhard Prix, Xavier Siemens, Scott Koranda, Yousuke Itoh
 *
 */
#include "config.h"

#include <math.h>
#include <stdio.h>
#include <strings.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <gsl/gsl_math.h>

#include <lal/LALString.h>
#include <lal/AVFactories.h>
#include <lal/LALInitBarycenter.h>
#include <lal/UserInput.h>
#include <lal/TranslateAngles.h>
#include <lal/TranslateMJD.h>
#include <lal/SFTfileIO.h>
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/FstatisticTools.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/ComputeFstat.h>
#include <lal/LALHough.h>
#include <lal/LogPrintf.h>
#include <lal/DopplerFullScan.h>
#include <lal/BinaryPulsarTiming.h>
#include <lal/TransientCW_utils.h>
#include <lal/LineRobustStats.h>
#include <lal/HeapToplist.h>
#include <lal/LALPulsarVCSInfo.h>

/*---------- DEFINES ----------*/

#define TRUE (1==1)
#define FALSE (1==0)
#define SQ(x) ( (x) * (x) )

/*----- SWITCHES -----*/
/*----- Macros -----*/

/** convert GPS-time to REAL8 */
#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

// NOTE: LAL's nan is more portable than either of libc or gsl !
#define LAL_NAN XLALREAL4FailNaN()
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
  REAL4 log10BSGL;			/**< Line-robust statistic \f$\log_{10}B_{\mathrm{SGL}}\f$ */
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
  UINT4 NFreq;			/**< number of frequency bins */
  REAL8 tauFstat;		/**< time to compute one Fstatistic over full data-duration (NSFT atoms) [in seconds]*/
  REAL8 tauTemplate;		/**< total loop time per template, includes candidate-handling (transient stats, toplist etc) */
  REAL8 tauF0;			/**< Demod timing constant = time per template per SFT */

  /* ----- transient-specific timings */
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
  RANKBY_BSGL = 2  	/**< rank candidates by BSGListic */
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
  DopplerRegion searchRegion;		    /**< parameter-space region to search over */
  DopplerFullScanState *scanState;          /**< current state of the Doppler-scan */
  PulsarDopplerParams stepSizes;	    /**< user-preferences on Doppler-param step-sizes */
  EphemerisData *ephemeris;		    /**< ephemeris data (from XLALInitBarycenter()) */
  UINT4 NSFTs;				    /**< total number of all SFTs used */
  UINT4Vector *numSFTsPerDet;		    /**< number of SFTs per detector, for log strings, etc. */
  LALStringVector *detectorIDs;		    /**< detector ID names, for column headings string */
  FstatInput *Fstat_in;		    /**< Fstat input data struct */
  FstatQuantities Fstat_what;		    /**< Fstat quantities to compute */
  toplist_t* FstatToplist;		    /**< sorted 'toplist' of the NumCandidatesToKeep loudest candidates */
  scanlineWindow_t *scanlineWindow;         /**< moving window of candidates on scanline to find local maxima */
  CHAR *VCSInfoString;                      /**< Git version string */
  CHAR *logstring;                          /**< log containing max-info on the whole search setup */
  transientWindowRange_t transientWindowRange; /**< search range parameters for transient window */
  BSGLSetup *BSGLsetup;                    /**< pre-computed setup for line-robust statistic */
  RankingStat_t RankingStatistic;           /**< rank candidates according to F or BSGL */
  BOOLEAN useResamp;
  UINT4 numFreqBins_FBand;
  REAL8 dFreq;
  PulsarParamsVector *injectionSources;    /**< Source parameters to inject: comma-separated list of file-patterns and/or direct config-strings ('{...}') */
  MultiLIGOTimeGPSVector *multiTimestamps; /**< a vector of timestamps (only set if provided from dedicated time stamp files) */
  MultiLALDetector multiIFO;		   /**< detectors to generate data for (if provided by user and not via noise files) */
  BOOLEAN runSearch;		   /**< whether to actually perform the search, or just generate a grid and/or count templates */
} ConfigVariables;


/* ----- User-variables: can be set from config-file or command-line */
typedef struct {
  INT4 Dterms;			/**< number of terms in LALDemod Dirichlet kernel is Dterms+1 */

  LALStringVector* assumeSqrtSX;/**< Assume stationary Gaussian noise with detector noise-floors sqrt{SX}" */
  BOOLEAN UseNoiseWeights;	/**< use SFT-specific noise-weights for each segment in Fstat-computation */

  REAL8 Freq;			/**< start-frequency of search */
  REAL8 FreqBand;		/**< Frequency-band to search over */
  REAL8 dFreq;			/**< user-specifyable Frequency stepsize */

  REAL8 Alpha;			/**< equatorial right-ascension in rad */
  REAL8 dAlpha;
  REAL8 AlphaBand;

  REAL8 Delta;			/**< equatorial declination in rad */
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
  REAL8 dorbitPeriod;
  REAL8 orbitPeriodBand;
  REAL8 orbitasini;		/**< amplitude of radial motion */
  REAL8 dorbitasini;
  REAL8 orbitasiniBand;
  LIGOTimeGPS orbitTp;		/**< epoch of periapse passage */
  REAL8 dorbitTp;
  REAL8 orbitTpBand;
  REAL8 orbitArgp;		/**< angle of periapse */
  REAL8 dorbitArgp;
  REAL8 orbitArgpBand;
  REAL8 orbitEcc;		/**< orbital eccentricity */
  REAL8 dorbitEcc;
  REAL8 orbitEccBand;

  /* extra parameters for --gridType==GRID_SPINDOWN_SQUARE */
  BOOLEAN strictSpindownBounds; /**< suppress spindown grid points outside the [fkdot,fkdotBand] ranges? */

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

  CHAR *outputLogfile;		/**< write a log-file */
  CHAR *outputFstat;		/**< filename to output Fstatistic in */
  CHAR *outputLoudest;		/**< filename for loudest F-candidate plus parameter estimation */

  CHAR *outputFstatHist;        /**< output discrete histogram of all Fstatistic values */
  REAL8 FstatHistBin;           /**< width of an Fstatistic histogram bin */

  BOOLEAN countTemplates;       /**< just count templates (if supported) instead of search */
  CHAR *outputGrid;		/**< filename to output grid points in */

  INT4 NumCandidatesToKeep;	/**< maximal number of toplist candidates to output */
  REAL8 FracCandidatesToKeep;	/**< fractional number of candidates to output in toplist */
  INT4 clusterOnScanline;	/**< number of points on "scanline" to use for 1-D local maxima finding */

  CHAR *gridFile;		/**< read template grid from this file */

  INT4 RngMedWindow;		/**< running-median window for noise floor estimation */
  LIGOTimeGPS refTime;		/**< reference-time for definition of pulsar-parameters [GPS] */

  int SSBprecision;		/**< full relativistic timing or Newtonian */

  LIGOTimeGPS minStartTime;	/**< Only use SFTs with timestamps starting from (including) this epoch (format 'xx.yy[GPS]' or 'xx.yyMJD') */
  LIGOTimeGPS maxStartTime;	/**< Only use SFTs with timestamps up to (excluding) this epoch (format 'xx.yy[GPS]' or 'xx.yyMJD') */
  CHAR *workingDir;		/**< directory to use for output files */
  REAL8 timerCount;		/**< output progress-meter every timerCount seconds */

  CHAR *outputFstatAtoms;	/**< output per-SFT, per-IFO 'atoms', ie quantities required to compute F-stat */

  BOOLEAN outputSingleFstats;	/**< in multi-detector case, also output single-detector F-stats */
  CHAR *RankingStatistic;	/**< rank candidates according to F-stat or BSGL */

  BOOLEAN computeBSGL;		/**< get single-IFO F-stats and compute line-robust statistic */
  BOOLEAN BSGLlogcorr;		/**< BSGL: compute log-correction (slower) or not (faster) */
  REAL8   Fstar0;		/**< BSGL: transition-scale parameter 'Fstar0', see documentation for XLALCreateBSGLSetup() for details */
  LALStringVector *oLGX;	/**< BSGL: prior per-detector line-vs-Gauss odds 'oLGX', see XLALCreateBSGLSetup() for details */
  REAL8 BSGLthreshold;		/**< output threshold on BSGL */

  CHAR *outputTransientStats;	/**< output file for transient B-stat values */
  CHAR *outputTransientStatsAll;	/**< output file for transient Fstat map F(t0,tau) for each candidate */
  CHAR *transient_WindowType;	/**< name of transient window ('none', 'rect', 'exp',...) */
  LIGOTimeGPS transient_t0Epoch;	/**< earliest GPS start-time for transient window search, in seconds */
  UINT4 transient_t0Offset;	/**< earliest start-time for transient window search, as offset in seconds from dataStartGPS */
  UINT4 transient_t0Band;	/**< Range of GPS start-times to search in transient search, in seconds */
  INT4  transient_dt0;		/**< Step-size for search/marginalization over transient-window start-time, in seconds */
  UINT4 transient_tau;	/**< smallest transient window length for marginalization, in seconds */
  UINT4 transient_tauBand;	/**<  Range of transient-window timescales to search, in seconds */
  INT4  transient_dtau;		/**< Step-size for search/marginalization over transient-window timescale, in seconds */
  BOOLEAN transient_useFReg;  	/**< FALSE: use 'standard' e^F for marginalization, TRUE: use e^FReg = (1/D)*e^F */

  CHAR *outputTiming;		/**< output timing measurements and parameters into this file [append!]*/
  CHAR *outputFstatTiming;	/**< output F-statistic timing measurements and parameters into this file [append!]*/

  int FstatMethod;		//!< select which method/algorithm to use to compute the F-statistic

  BOOLEAN resampFFTPowerOf2;	//!< in Resamp: enforce FFT length to be a power of two (by rounding up)
  REAL8 allowedMismatchFromSFTLength; /**< maximum allowed mismatch from SFTs being too long */

  LALStringVector *injectionSources;    /**< Source parameters to inject: comma-separated list of file-patterns and/or direct config-strings ('{...}') */
  LALStringVector *injectSqrtSX; 	/**< Add Gaussian noise: list of respective detectors' noise-floors sqrt{Sn}" */
  LALStringVector *IFOs;	/**< list of detector-names "H1,H2,L1,.." or single detector */
  LALStringVector *timestampsFiles;     /**< Names of numDet timestamps files */
  REAL8 Tsft;                           /**< length of one SFT in seconds, used in combination with timestamps files (otherwise taken from SFT files) */
  INT4 randSeed;		/**< allow user to specify random-number seed for reproducible noise-realizations */

} UserInput_t;

/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lal/lib/std/LALError.c */

/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);
int initUserVars ( UserInput_t *uvar);
int InitFstat ( ConfigVariables *cfg, const UserInput_t *uvar );
void Freemem(ConfigVariables *cfg);

int checkUserInputConsistency (const UserInput_t *uvar);
int outputBeamTS( const CHAR *fname, const AMCoeffs *amcoe, const DetectorStateSeries *detStates );
MultiNoiseWeights *getUnitWeights ( const MultiSFTVector *multiSFTs );

int write_FstatCandidate_to_fp ( FILE *fp, const FstatCandidate *thisFCand, const BOOLEAN output_stats, const BOOLEAN output_orbit );
int write_PulsarCandidate_to_fp ( FILE *fp,  const PulsarCandidate *pulsarParams, const FstatCandidate *Fcand );

int compareFstatCandidates ( const void *candA, const void *candB );
int compareFstatCandidates_BSGL ( const void *candA, const void *candB );

int WriteFstatLog ( const CHAR *log_fname, const CHAR *logstr );
CHAR *XLALGetLogString ( const ConfigVariables *cfg, const BOOLEAN verbose );

int write_TimingInfo ( const CHAR *timingFile, const timingInfo_t *ti, const ConfigVariables *cfg );

gsl_vector_int *resize_histogram(gsl_vector_int *old_hist, size_t size);

/* ---------- scanline window functions ---------- */
scanlineWindow_t *XLALCreateScanlineWindow ( UINT4 windowWings );
void XLALDestroyScanlineWindow ( scanlineWindow_t *scanlineWindow );
int XLALAdvanceScanlineWindow ( const FstatCandidate *nextCand, scanlineWindow_t *scanWindow );
BOOLEAN XLALCenterIsLocalMax ( const scanlineWindow_t *scanWindow, const UINT4 sortingStatistic );

/* ----- which timing function to use ----- */
#define GETTIME XLALGetCPUTime

/*----------------------------------------------------------------------*/
/* Function definitions start here */
/*----------------------------------------------------------------------*/

/**
 * MAIN function of ComputeFstatistic code.
 * Calculate the F-statistic over a given portion of the parameter-space
 * and write a list of 'candidates' into a file(default: 'Fstats').
 */
int main(int argc,char *argv[])
{
  FILE *fpFstat = NULL;
  LALFILE *fpTransientStats = NULL, *fpTransientStatsAll = NULL;
  REAL8 numTemplates, templateCounter;
  time_t clock0;
  PulsarDopplerParams XLAL_INIT_DECL(dopplerpos);
  FstatCandidate XLAL_INIT_DECL(loudestFCand);
  FstatCandidate XLAL_INIT_DECL(thisFCand);
  gsl_vector_int *Fstat_histogram = NULL;

  UserInput_t XLAL_INIT_DECL(uvar);
  ConfigVariables XLAL_INIT_DECL(GV);

  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register all user-variable */
  XLAL_CHECK_MAIN ( initUserVars ( &uvar ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK_MAIN ( (GV.VCSInfoString = XLALVCSInfoString(lalPulsarVCSInfoList, 0, "%% ")) != NULL, XLAL_EFUNC );

  /* do ALL cmdline and cfgfile handling */
  BOOLEAN should_exit = 0;
  XLAL_CHECK_MAIN( XLALUserVarReadAllInput( &should_exit, argc, argv, lalPulsarVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit ) {
    exit (1);
  }

  /* do some sanity checks on the user-input before we proceed */
  XLAL_CHECK_MAIN ( checkUserInputConsistency ( &uvar ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* Initialization the common variables of the code, */
  /* like ephemeries data and template grids: */
  XLAL_CHECK_MAIN ( InitFstat( &GV, &uvar) == XLAL_SUCCESS, XLAL_EFUNC );

  /* ----- produce a log-string describing the specific run setup ----- */
  BOOLEAN logstrVerbose = TRUE;
  XLAL_CHECK_MAIN ( (GV.logstring = XLALGetLogString ( &GV, logstrVerbose )) != NULL, XLAL_EFUNC );
  LogPrintfVerbatim( LOG_NORMAL, "%s", GV.logstring );

  /* keep a log-file recording all relevant parameters of this search-run */
  if ( uvar.outputLogfile ) {
    XLAL_CHECK_MAIN ( WriteFstatLog ( uvar.outputLogfile, GV.logstring ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  /* count number of binary orbit parameters (done early to be able to prepare output file headers) */
  const UINT4 n_orbitasini = 1 + ( XLALUserVarWasSet(&uvar.dorbitasini) ? (UINT4) floor ( uvar.orbitasiniBand / uvar.dorbitasini ) : 0 );
  const UINT4 n_orbitPeriod = 1 + ( XLALUserVarWasSet(&uvar.dorbitPeriod) ? (UINT4) floor ( uvar.orbitPeriodBand / uvar.dorbitPeriod ) : 0 );
  const UINT4 n_orbitTp = 1 + ( XLALUserVarWasSet(&uvar.dorbitTp) ? (UINT4) floor ( uvar.orbitTpBand / uvar.dorbitTp ) : 0 );
  const UINT4 n_orbitArgp = 1 + ( XLALUserVarWasSet(&uvar.dorbitArgp) ? (UINT4) floor ( uvar.orbitArgpBand / uvar.dorbitArgp ) : 0 );
  const UINT4 n_orbitEcc = 1 + ( XLALUserVarWasSet(&uvar.dorbitEcc) ? (UINT4) floor ( uvar.orbitEccBand / uvar.dorbitEcc ) : 0 );
  const UINT4 n_orbit = n_orbitasini * n_orbitPeriod * n_orbitTp * n_orbitArgp * n_orbitEcc;

  /* if a complete output of the F-statistic results, or a grid output file, was requested,
   * we open and prepare the output-file here */
  if ( ( uvar.outputFstat && GV.runSearch ) || uvar.outputGrid )
    {
      if ( uvar.outputGrid )
        {
          XLAL_CHECK_MAIN ( (fpFstat = fopen (uvar.outputGrid, "wb")) != NULL, XLAL_ESYS, "\nError opening file '%s' for writing..\n\n", uvar.outputGrid );
        } else {
          XLAL_CHECK_MAIN ( (fpFstat = fopen (uvar.outputFstat, "wb")) != NULL, XLAL_ESYS, "\nError opening file '%s' for writing..\n\n", uvar.outputFstat );
        }
      fprintf (fpFstat, "%s", GV.logstring );

      /* assemble column headings string */
      char column_headings_string[1024];
      XLAL_INIT_MEM( column_headings_string );
      strcat ( column_headings_string, "freq alpha delta f1dot f2dot f3dot" );
      if ( !uvar.outputGrid )
        {
          strcat ( column_headings_string, " 2F");
          if ( uvar.computeBSGL )
            {
              strcat ( column_headings_string, " log10BSGL" );
            }
          if ( uvar.outputSingleFstats || uvar.computeBSGL )
            {
              const UINT4 numDetectors = GV.detectorIDs->length;
              for ( UINT4 X = 0; X < numDetectors ; X ++ )
                {
                  char headingX[7];
                  snprintf ( headingX, sizeof(headingX), " 2F_%s", GV.detectorIDs->data[X] );
                  strcat ( column_headings_string, headingX );
                } /* for X < numDet */
            }
        }
      if ( n_orbit > 1 )
        {
          strcat ( column_headings_string, " asini period tp argp ecc" );
      }
      fprintf (fpFstat, "%%%% columns:\n%%%% %s\n", column_headings_string );

    } /* if ( ( uvar.outputFstat && GV.runSearch ) || uvar.outputGrid ) */

  /* count number of templates */
  numTemplates = XLALNumDopplerTemplates ( GV.scanState );
  numTemplates *= n_orbit;
  if (uvar.countTemplates) {
    printf("%%%% Number of templates: %0.0f\n", numTemplates);
  }

  if ( uvar.outputGrid ) { /* alternative to main search loop: loop through the same grid but only write out the parameters, no F-stats */
    while ( XLALNextDopplerPos( &dopplerpos, GV.scanState ) == 0 ) {
    for (UINT4 i_orbitasini = 0; i_orbitasini < n_orbitasini; ++i_orbitasini) {
    for (UINT4 i_orbitPeriod = 0; i_orbitPeriod < n_orbitPeriod; ++i_orbitPeriod) {
    for (UINT4 i_orbitTp = 0; i_orbitTp < n_orbitTp; ++i_orbitTp) {
    for (UINT4 i_orbitArgp = 0; i_orbitArgp < n_orbitArgp; ++i_orbitArgp) {
    for (UINT4 i_orbitEcc = 0; i_orbitEcc < n_orbitEcc; ++i_orbitEcc) {
      dopplerpos.asini = uvar.orbitasini + i_orbitasini * uvar.dorbitasini;
      dopplerpos.period = uvar.orbitPeriod + i_orbitPeriod * uvar.dorbitPeriod;
      dopplerpos.tp = uvar.orbitTp; XLALGPSAdd( &dopplerpos.tp, i_orbitTp * uvar.dorbitTp );
      dopplerpos.argp = uvar.orbitArgp + i_orbitArgp * uvar.dorbitArgp;
      dopplerpos.ecc = uvar.orbitEcc + i_orbitEcc * uvar.dorbitEcc;
      for ( UINT4 iFreq = 0; iFreq < GV.numFreqBins_FBand; iFreq ++ )
      {
        /* collect data on current 'Fstat-candidate' */
        thisFCand.doppler = dopplerpos;	// use 'original' dopplerpos @ refTime !
        thisFCand.doppler.fkdot[0] += iFreq * GV.dFreq; // this only does something for the resampling post-loop over frequency-bins, 0 otherwise ...
        if ( fpFstat ) { /* no search, no toplist: write out grid point immediately */
          if ( write_FstatCandidate_to_fp ( fpFstat, &thisFCand, 0, n_orbit > 1 ) != 0 ) {
            LogPrintf (LOG_CRITICAL, "Failed to write candidate to file.\n");
            return -1;
          }
        }
      }
    }
    }
    }
    }
    }
    } /* while more Doppler positions to scan */
  } /* if ( uvar.outputGrid ) */

  if ( uvar.outputTransientStats && GV.runSearch )
    {
      XLAL_CHECK_MAIN ( (fpTransientStats = XLALFileOpen (uvar.outputTransientStats, "wb")) != NULL, XLAL_ESYS, "\nError opening file '%s' for writing..\n\n", uvar.outputTransientStats );
      XLALFilePrintf (fpTransientStats, "%s", GV.logstring );			/* write search log comment */
      XLAL_CHECK_MAIN ( write_transientCandidate_to_fp ( fpTransientStats, NULL, 's' ) == XLAL_SUCCESS, XLAL_EFUNC );	/* write header-line comment */
    }

  if ( uvar.outputTransientStatsAll && GV.runSearch )
    {
      XLAL_CHECK_MAIN ( (fpTransientStatsAll = XLALFileOpen (uvar.outputTransientStatsAll, "wb")) != NULL, XLAL_ESYS, "\nError opening file '%s' for writing..\n\n", uvar.outputTransientStatsAll );
      XLALFilePrintf (fpTransientStatsAll, "%s", GV.logstring );			/* write search log comment */
      XLAL_CHECK_MAIN ( write_transientCandidateAll_to_fp ( fpTransientStatsAll, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );	/* write header-line comment */
    }

  /* start Fstatistic histogram with a single empty bin */
  if ( uvar.outputFstatHist && GV.runSearch ) {
    XLAL_CHECK_MAIN ((Fstat_histogram = gsl_vector_int_alloc(1)) != NULL, XLAL_ENOMEM );
    gsl_vector_int_set_zero(Fstat_histogram);
  }

  /*----------------------------------------------------------------------
   * main loop: demodulate data for each point in the sky-position grid
   * and for each value of the frequency-spindown
   */
  templateCounter = 0.0;
  clock0 = GETTIME();

  // ----- prepare timing info
  REAL8 tic0, tic, toc, timeOfLastProgressUpdate = 0;	// high-precision timing counters
  timingInfo_t XLAL_INIT_DECL(timing);			// timings of Fstatistic computation, transient Fstat-map, transient Bayes factor

  // pointer to Fstat results structure, will be allocated by XLALComputeFstat()
  FstatResults* Fstat_res = NULL;

  while ( GV.runSearch && ( XLALNextDopplerPos( &dopplerpos, GV.scanState ) == 0 ) ) {

    for (UINT4 i_orbitasini = 0; i_orbitasini < n_orbitasini; ++i_orbitasini) {
    for (UINT4 i_orbitPeriod = 0; i_orbitPeriod < n_orbitPeriod; ++i_orbitPeriod) {
    for (UINT4 i_orbitTp = 0; i_orbitTp < n_orbitTp; ++i_orbitTp) {
    for (UINT4 i_orbitArgp = 0; i_orbitArgp < n_orbitArgp; ++i_orbitArgp) {
    for (UINT4 i_orbitEcc = 0; i_orbitEcc < n_orbitEcc; ++i_orbitEcc) {
    
      dopplerpos.asini = uvar.orbitasini + i_orbitasini * uvar.dorbitasini;
      dopplerpos.period = uvar.orbitPeriod + i_orbitPeriod * uvar.dorbitPeriod;
      dopplerpos.tp = uvar.orbitTp; XLALGPSAdd( &dopplerpos.tp, i_orbitTp * uvar.dorbitTp );
      dopplerpos.argp = uvar.orbitArgp + i_orbitArgp * uvar.dorbitArgp;
      dopplerpos.ecc = uvar.orbitEcc + i_orbitEcc * uvar.dorbitEcc;

      tic0 = tic = GETTIME();

      /* main function call: compute F-statistic for this template */
      XLAL_CHECK_MAIN ( XLALComputeFstat ( &Fstat_res, GV.Fstat_in, &dopplerpos, GV.numFreqBins_FBand, GV.Fstat_what) == XLAL_SUCCESS, XLAL_EFUNC );

      toc = GETTIME();
      timing.tauFstat += (toc - tic);   // pure Fstat-calculation time

      /* Progress meter */
      templateCounter += 1.0;
      if ( LogLevel() >= LOG_NORMAL && ( (toc - timeOfLastProgressUpdate) > uvar.timerCount) )
        {
          REAL8 diffSec = GETTIME() - clock0 ;  /* seconds since start of loop*/
          REAL8 taup = diffSec / templateCounter ;
          REAL8 timeLeft = (numTemplates - templateCounter) *  taup;
          LogPrintf (LOG_NORMAL, "Progress: %g/%g = %.2f %% done, Estimated time left: %.0f s\n",
                     templateCounter, numTemplates, templateCounter/numTemplates * 100.0, timeLeft);
          timeOfLastProgressUpdate = toc;
        }

      // here we use Santiago's trick to hack the resampling Fstat(f) into the single-F rest of the
      // main-loop: we simply loop the remaining body over all frequency-bins in the Fstat-vector,
      // this way nothing needs to be changed!  in the non-resampling case, this loop iterates only
      // once, so nothing is changed ...
      for ( UINT4 iFreq = 0; iFreq < GV.numFreqBins_FBand; iFreq ++ )
      {

        /* collect data on current 'Fstat-candidate' */
        thisFCand.doppler = dopplerpos;	// use 'original' dopplerpos @ refTime !
        thisFCand.doppler.fkdot[0] += iFreq * GV.dFreq; // this only does something for the resampling post-loop over frequency-bins, 0 otherwise ...
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
          thisFCand.FaFb_refTime = Fstat_res->refTimePhase; // 'internal' reference time, used only for global phase estimate
          thisFCand.Fa = Fstat_res->Fa[iFreq];
          thisFCand.Fb = Fstat_res->Fb[iFreq];
        } else {
          thisFCand.Fa = thisFCand.Fb = crect(LAL_NAN,LAL_NAN);
        }
        thisFCand.Mmunu = Fstat_res->Mmunu;
        MultiFstatAtomVector* thisFAtoms = NULL;
        if (GV.Fstat_what & FSTATQ_ATOMS_PER_DET) {
          thisFAtoms = Fstat_res->multiFatoms[iFreq];
        }

        /* sanity check on the result */
        if ( !isfinite(thisFCand.twoF) )
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

      if ( uvar.computeBSGL )
        {
          thisFCand.log10BSGL = XLALComputeBSGL ( thisFCand.twoF, thisFCand.twoFX, GV.BSGLsetup );
          XLAL_CHECK_MAIN ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeBSGL() failed with errno=%d\n", xlalErrno );
        }
      else
        {
          thisFCand.log10BSGL = LAL_NAN; /* in non-BSGL case, block field with NAN, needed for output checking in write_PulsarCandidate_to_fp() */
        }

      /* push new value onto scan-line buffer */
      XLALAdvanceScanlineWindow ( &thisFCand, GV.scanlineWindow );

      /* two types of threshold: fixed (TwoF- and/or BSGL-threshold) and dynamic (NumCandidatesToKeep) */
      BOOLEAN is1DlocalMax = FALSE;
      if ( XLALCenterIsLocalMax ( GV.scanlineWindow, GV.RankingStatistic ) ) /* must be 1D local maximum */
        is1DlocalMax = TRUE;
      BOOLEAN isOver2FThreshold = FALSE; /* will always be checked, so start at 'FALSE' */
      if ( GV.scanlineWindow->center->twoF >= uvar.TwoFthreshold ) /* fixed 2F threshold */
        isOver2FThreshold = TRUE;
      BOOLEAN isOverBSGLthreshold = TRUE;  /* will not be checked in non-BSGL case, so start at 'TRUE' */
      if ( uvar.computeBSGL && ( GV.scanlineWindow->center->log10BSGL < uvar.BSGLthreshold ) ) /* fixed threshold on log10BSGL */
        isOverBSGLthreshold = FALSE;
      if ( is1DlocalMax && isOver2FThreshold && isOverBSGLthreshold )
        {
	  FstatCandidate *writeCand = GV.scanlineWindow->center;

	  /* insert this into toplist if requested */
	  if ( GV.FstatToplist  )			/* dynamic threshold */
	    {
	      if ( insert_into_toplist(GV.FstatToplist, (void*)writeCand ) ) {
		LogPrintf ( LOG_DEBUG, "Added new candidate into toplist: 2F = %f", writeCand->twoF );
		if ( uvar.computeBSGL ) {
		  LogPrintfVerbatim ( LOG_DEBUG, ", 2F_H1 = %f, 2F_L1 = %f, log10BSGL = %f", writeCand->twoFX[0], writeCand->twoFX[1], writeCand->log10BSGL );
                }
	      }
	      else {
		LogPrintf ( LOG_DEBUG, "NOT added the candidate into toplist: 2F = %f", writeCand->twoF );
		if ( uvar.computeBSGL ) {
		  LogPrintfVerbatim ( LOG_DEBUG, ", 2F_H1 = %f, 2F_L1 = %f, log10BSGL = %f", writeCand->twoFX[0], writeCand->twoFX[1], writeCand->log10BSGL );
                }
	      }
	      LogPrintfVerbatim ( LOG_DEBUG, "\n" );
	    }
	  else if ( fpFstat ) 				/* no toplist :write out immediately */
	    {
	      if ( write_FstatCandidate_to_fp ( fpFstat, writeCand, 1, n_orbit > 1 ) != 0 )
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
          if ( thisFCand.twoF > loudestFCand.twoF ) {
            loudestFCand = thisFCand;
          }
          break;
        case RANKBY_BSGL:
          if ( thisFCand.log10BSGL > loudestFCand.log10BSGL ) {
            loudestFCand = thisFCand;
          }
          break;
        default:
          XLAL_ERROR_MAIN ( XLAL_EINVAL, "Invalid ranking statistic '%d', supported are 'F=0', and 'BSGL=2'\n", GV.RankingStatistic );
          break;
        }

      /* add Fstatistic to histogram if needed */
      if (uvar.outputFstatHist) {

	/* compute bin */
	const size_t bin = thisFCand.twoF / uvar.FstatHistBin;

	/* resize histogram vector if needed */
	if (!Fstat_histogram || bin >= Fstat_histogram->size) {
	  XLAL_CHECK_MAIN ( (Fstat_histogram = resize_histogram(Fstat_histogram, bin + 1)) != NULL, XLAL_EFUNC, "\nCouldn't (re)allocate 'Fstat_histogram'\n" );
        }

	/* add to bin */
	gsl_vector_int_set(Fstat_histogram, bin,
			   gsl_vector_int_get(Fstat_histogram, bin) + 1);

      }


      /* ----- output F-stat atoms into files: one file per Doppler-position ---------- */
      /* F-stat 'atoms' = per-SFT components {a,b,Fa,Fb}_alpha) */
      if (uvar.outputFstatAtoms)
	{
	  LALFILE *fpFstatAtoms = NULL;
	  CHAR *fnameAtoms = NULL;
	  CHAR *dopplerName;
	  UINT4 len;

	  XLAL_CHECK_MAIN ( (dopplerName = XLALPulsarDopplerParams2String ( &dopplerpos )) != NULL, XLAL_EFUNC );
	  len = strlen(uvar.outputFstatAtoms) + strlen(dopplerName) + 10;
	  XLAL_CHECK_MAIN ( (fnameAtoms = LALMalloc (len)) != NULL, XLAL_ENOMEM, "Failed to LALMalloc(%d)\n", len );
	  sprintf (fnameAtoms, "%s_%s.dat", uvar.outputFstatAtoms, dopplerName );

	  XLAL_CHECK_MAIN ( (fpFstatAtoms = XLALFileOpen (fnameAtoms, "wb")) != NULL, XLAL_ESYS, "Error opening file '%s' for writing..\n\n", fnameAtoms );
	  XLALFree ( fnameAtoms );
	  XLALFree ( dopplerName );

	  XLALFilePrintf (fpFstatAtoms, "%s", GV.logstring );

	  XLAL_CHECK_MAIN ( write_MultiFstatAtoms_to_fp ( fpFstatAtoms, thisFAtoms ) == XLAL_SUCCESS, XLAL_EFUNC );
	  XLALFileClose (fpFstatAtoms);

	} /* if outputFstatAtoms */

      /* ----- compute transient-CW statistics if their output was requested  ----- */
      if ( fpTransientStats || fpTransientStatsAll )
        {
          transientCandidate_t XLAL_INIT_DECL(transientCand);

          /* compute Fstat map F_mn over {t0, tau} */
          tic = GETTIME();
          XLAL_CHECK_MAIN ( (transientCand.FstatMap = XLALComputeTransientFstatMap ( thisFAtoms, GV.transientWindowRange, uvar.transient_useFReg)) != NULL, XLAL_EFUNC );
          toc = GETTIME();
          timing.tauTransFstatMap += (toc - tic); // time to compute transient Fstat-map

          /* compute marginalized Bayes factor */
          tic = GETTIME();
          transientCand.logBstat = XLALComputeTransientBstat ( GV.transientWindowRange, transientCand.FstatMap );
          XLAL_CHECK_MAIN ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeTransientBstat() failed with xlalErrno = %d\n", xlalErrno );
          toc = GETTIME();
          timing.tauTransMarg += (toc - tic);

          /* ----- compute parameter posteriors for {t0, tau} */
          pdf1D_t *pdf_t0  = NULL;
          pdf1D_t *pdf_tau = NULL;

          XLAL_CHECK_MAIN ( (pdf_t0 = XLALComputeTransientPosterior_t0 ( GV.transientWindowRange, transientCand.FstatMap )) != NULL, XLAL_EFUNC );
          XLAL_CHECK_MAIN ( (pdf_tau = XLALComputeTransientPosterior_tau ( GV.transientWindowRange, transientCand.FstatMap )) != NULL, XLAL_EFUNC );

          if ( fpTransientStats )
            {
              /* get maximum-posterior estimate (MP) from the modes of these pdfs */
              transientCand.t0_MP = XLALFindModeOfPDF1D ( pdf_t0 );
              XLAL_CHECK_MAIN ( xlalErrno == 0, XLAL_EFUNC, "mode-estimation failed for pdf_t0. xlalErrno = %d\n", xlalErrno );
              transientCand.tau_MP =  XLALFindModeOfPDF1D ( pdf_tau );
              XLAL_CHECK_MAIN ( xlalErrno == 0, XLAL_EFUNC, "mode-estimation failed for pdf_tau. xlalErrno = %d\n", xlalErrno );
            }

          /* record timing-relevant transient search params */
          timing.tauMin  = GV.transientWindowRange.tau;
          timing.tauMax  = timing.tauMin + GV.transientWindowRange.tauBand;
          timing.NStart  = transientCand.FstatMap->F_mn->size1;
          timing.NTau    = transientCand.FstatMap->F_mn->size2;

          /* add meta-info on current transient-CW candidate */
          transientCand.doppler = dopplerpos;
          transientCand.windowRange = GV.transientWindowRange;

          if ( fpTransientStats )
            {
              /* output everything into stats-file (one line per candidate) */
              XLAL_CHECK_MAIN ( write_transientCandidate_to_fp ( fpTransientStats, &transientCand, 's' ) == XLAL_SUCCESS, XLAL_EFUNC );
            }

          if ( fpTransientStatsAll )
            {
              /* output everything into stats-file (one block for the whole (t0,tau) grid per candidate) */
              XLAL_CHECK_MAIN ( write_transientCandidateAll_to_fp ( fpTransientStatsAll, &transientCand ) == XLAL_SUCCESS, XLAL_EFUNC );
            }

          /* free dynamically allocated F-stat map */
          XLALDestroyTransientFstatMap ( transientCand.FstatMap );
          XLALDestroyPDF1D ( pdf_t0 );
          XLALDestroyPDF1D ( pdf_tau );

        } /* if ( fpTransientStats || fpTransientStatsAll ) */

      } // for ( iFreq < numFreqBins_FBand )

      /* now measure total loop time per template */
      toc = GETTIME();
      timing.tauTemplate += (toc - tic0);

    }
    }
    }
    }
    }

  } /* while more Doppler positions to scan */


  /* if requested: output timings into timing-file */
  if ( uvar.outputTiming && GV.runSearch )
    {
      REAL8 num_templates = numTemplates * GV.numFreqBins_FBand;	// 'templates' now refers to number of 'frequency-bands' in resampling case

      timing.NSFTs = GV.NSFTs;
      timing.NFreq = (UINT4) ( 1 + floor ( GV.searchRegion.fkdotBand[0] / GV.dFreq ) );

      // compute averages:
      timing.tauFstat    /= num_templates;
      timing.tauTemplate /= num_templates;
      timing.tauF0       =  timing.tauFstat / timing.NSFTs;
      timing.tauTransFstatMap /= num_templates;
      timing.tauTransMarg     /= num_templates;

      XLAL_CHECK_MAIN ( write_TimingInfo ( uvar.outputTiming, &timing, &GV ) == XLAL_SUCCESS, XLAL_EFUNC );

    } /* if timing output requested */

  /* if requested: output F-statistic timings into F-statistic-timing-file */
  if ( uvar.outputFstatTiming && GV.runSearch )
    {

      FILE *fp;
      if ( (fp = fopen( uvar.outputFstatTiming, "rb" )) != NULL )
        {
          fclose(fp);
          XLAL_CHECK ( (fp = fopen( uvar.outputFstatTiming, "ab" ) ), XLAL_ESYS, "Failed to open existing timing-file '%s' for appending\n", uvar.outputFstatTiming );
          XLAL_CHECK_MAIN ( XLALAppendFstatTiming2File ( GV.Fstat_in, fp, 0 ) == XLAL_SUCCESS, XLAL_EFUNC );
        }
      else
        {
          XLAL_CHECK ( (fp = fopen( uvar.outputFstatTiming, "wb" ) ), XLAL_ESYS, "Failed to open new timing-file '%s' for writing\n", uvar.outputFstatTiming );
          XLAL_CHECK_MAIN ( XLALAppendFstatTiming2File ( GV.Fstat_in, fp, 1 ) == XLAL_SUCCESS, XLAL_EFUNC );
        }
      fclose(fp);

    } /* if timing output requested */

  /* ----- if using toplist: sort and write it out to file now ----- */
  if ( fpFstat && GV.FstatToplist && GV.runSearch )
    {
      UINT4 el;

      /* sort toplist */
      LogPrintf ( LOG_NORMAL, "Sorting toplist ... ");
      if ( GV.RankingStatistic == RANKBY_2F )
        qsort_toplist ( GV.FstatToplist, compareFstatCandidates );
      else if ( GV.RankingStatistic == RANKBY_BSGL )
        qsort_toplist ( GV.FstatToplist, compareFstatCandidates_BSGL );
      else
        XLAL_ERROR ( XLAL_EINVAL, "Ranking statistic '%d' undefined here, allowed are 'F=0' and 'BSGL=2'\n", GV.RankingStatistic );
      LogPrintfVerbatim ( LOG_NORMAL, "done.\n");

      for ( el=0; el < GV.FstatToplist->elems; el ++ )
	{
	  const FstatCandidate *candi;
	  if ( ( candi = (const FstatCandidate *) toplist_elem ( GV.FstatToplist, el )) == NULL ) {
	    LogPrintf ( LOG_CRITICAL, "Internal consistency problems with toplist: contains fewer elements than expected!\n");
	    return -1;
	  }
	  if ( write_FstatCandidate_to_fp ( fpFstat, candi, 1, n_orbit > 1 ) != 0 )
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
  XLALFileClose (fpTransientStats);
  XLALFileClose (fpTransientStatsAll);

  /* ----- estimate amplitude-parameters for the loudest canidate and output into separate file ----- */
  if ( uvar.outputLoudest && GV.runSearch )
    {
      FILE *fpLoudest;
      PulsarCandidate XLAL_INIT_DECL(pulsarParams);
      pulsarParams.Doppler = loudestFCand.doppler;

      XLAL_CHECK_MAIN ( XLALEstimatePulsarAmplitudeParams ( &pulsarParams, &loudestFCand.FaFb_refTime, loudestFCand.Fa, loudestFCand.Fb, &loudestFCand.Mmunu ) == XLAL_SUCCESS, XLAL_EFUNC );

      XLAL_CHECK_MAIN ( (fpLoudest = fopen (uvar.outputLoudest, "wb")) != NULL, XLAL_ESYS, "Error opening file '%s' for writing..\n\n", uvar.outputLoudest );

      /* write header with run-info */
      fprintf (fpLoudest, "%s", GV.logstring );

      /* write this 'candidate' to disc */
      XLAL_CHECK_MAIN ( write_PulsarCandidate_to_fp ( fpLoudest,  &pulsarParams, &loudestFCand) == XLAL_SUCCESS, XLAL_EFUNC );
      fclose (fpLoudest);

      gsl_matrix_free ( pulsarParams.AmpFisherMatrix );

    } /* write loudest candidate to file */

  LogPrintf (LOG_NORMAL, "Search finished.\n");

  /* write out the Fstatistic histogram */
  if  (uvar.outputFstatHist && GV.runSearch ) {

    size_t i = 0;
    FILE *fpFstatHist = fopen(uvar.outputFstatHist, "wb");

    XLAL_CHECK_MAIN (fpFstatHist != NULL, XLAL_ESYS, "\nError opening file '%s' for writing..\n\n", uvar.outputFstat );
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
  XLALDestroyDopplerFullScan ( GV.scanState);
  XLALDestroyFstatResults ( Fstat_res );

  Freemem ( &GV );

  if (Fstat_histogram) {
    gsl_vector_int_free(Fstat_histogram);
  }

  /* did we forget anything ? */
  LALCheckMemoryLeaks();

  return 0;

} /* main() */


/**
 * Register all our "user-variables" that can be specified from cmd-line and/or config-file.
 * Here we set defaults for some user-variables and register them with the UserInput module.
 */
int
initUserVars ( UserInput_t *uvar )
{
  XLAL_CHECK ( uvar != NULL, XLAL_EINVAL );

  /* set a few defaults */
  uvar->FreqBand = 0.0;
  uvar->Alpha 	= 0.0;
  uvar->Delta 	= 0.0;
  uvar->AlphaBand = 0;
  uvar->DeltaBand = 0;
  uvar->skyRegion = NULL;
  uvar->Dterms = FstatOptionalArgsDefaults.Dterms;

  uvar->ephemEarth = XLALStringDuplicate("earth00-40-DE405.dat.gz");
  uvar->ephemSun = XLALStringDuplicate("sun00-40-DE405.dat.gz");

  uvar->assumeSqrtSX = NULL;
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

  uvar->outputLogfile = NULL;
  uvar->outputFstat = NULL;
  uvar->outputLoudest = NULL;

  uvar->outputFstatHist = NULL;
  uvar->FstatHistBin = 0.1;

  uvar->countTemplates = FALSE;

  uvar->gridFile = NULL;

  uvar->RngMedWindow = FstatOptionalArgsDefaults.runningMedianWindow;

  uvar->SSBprecision = FstatOptionalArgsDefaults.SSBprec;

  uvar->minStartTime.gpsSeconds = 0;
  uvar->maxStartTime.gpsSeconds = LAL_INT4_MAX;

  uvar->workingDir = XLALStringDuplicate ( "." );

  uvar->timerCount = 10;	/* output a timer/progress count every N seconds */

  uvar->strictSpindownBounds = FALSE;

  uvar->spindownAge = 0.0;
  uvar->minBraking = 0.0;
  uvar->maxBraking = 0.0;

  uvar->FstatMethod = FstatOptionalArgsDefaults.FstatMethod;

  uvar->outputSingleFstats = FALSE;
  uvar->RankingStatistic = XLALStringDuplicate ( "F" );

  uvar->computeBSGL = FALSE;
  uvar->BSGLlogcorr = TRUE;
  uvar->Fstar0 = 0.0;
  uvar->oLGX = NULL;       /* NULL is intepreted as oLGX[X] = 1.0/Ndet for all X */
  uvar->BSGLthreshold = - LAL_REAL8_MAX;

  uvar->transient_WindowType = XLALStringDuplicate ( "none" );
  uvar->transient_useFReg = 0;
  uvar->resampFFTPowerOf2 = FstatOptionalArgsDefaults.resampFFTPowerOf2;
  uvar->allowedMismatchFromSFTLength = 0;
  uvar->injectionSources = NULL;
  uvar->injectSqrtSX = NULL;
  uvar->IFOs = NULL;
  uvar->timestampsFiles = NULL;

  uvar->Tsft=1800.0;

  /* ---------- register all user-variables ---------- */
  XLALRegisterUvarMember( 	Alpha, 		RAJ, 'a', OPTIONAL, "Sky: equatorial J2000 right ascension (in radians or hours:minutes:seconds)");
  XLALRegisterUvarMember( 	Delta, 		DECJ, 'd', OPTIONAL, "Sky: equatorial J2000 declination (in radians or degrees:minutes:seconds)");
  XLALRegisterUvarMember( 	skyRegion,      STRING, 'R', OPTIONAL, "ALTERNATIVE: Sky-region polygon '(Alpha1,Delta1),(Alpha2,Delta2),...' or 'allsky'");

  XLALRegisterUvarMember( 	Freq, 		REAL8, 'f', OPTIONAL, "Starting search frequency Freq in Hz");
  XLALRegisterUvarMember( 	f1dot, 		REAL8, 's', OPTIONAL, "First spindown parameter  f1dot = dFreq/dt");
  XLALRegisterUvarMember( 	f2dot, 		 REAL8, 0 , OPTIONAL, "Second spindown parameter f2dot = d^2Freq/dt^2");
  XLALRegisterUvarMember( 	f3dot, 		 REAL8, 0 , OPTIONAL, "Third spindown parameter  f3dot = d^3Freq/dt^3");

  XLALRegisterUvarMember( 	AlphaBand, 	RAJ, 'z', OPTIONAL, "Sky: search band from Alpha to Alpha+AlphaBand (in radians or h:m:s)");
  XLALRegisterUvarMember( 	DeltaBand, 	DECJ, 'c', OPTIONAL, "Sky: search band from Delta to Delta+DeltaBand (in radians or d:m:s)");
  XLALRegisterUvarMember( 	FreqBand, 	REAL8, 'b', OPTIONAL, "Search band in frequency in Hz");
  XLALRegisterUvarMember( 	f1dotBand, 	REAL8, 'm', OPTIONAL, "Search band in f1dot in Hz/s");
  XLALRegisterUvarMember( 	f2dotBand, 	 REAL8, 0 , OPTIONAL, "Search band in f2dot in Hz/s^2");
  XLALRegisterUvarMember( 	f3dotBand, 	 REAL8, 0 , OPTIONAL, "Search band in f3dot in Hz/s^3");

  XLALRegisterUvarMember( 	dAlpha, 	RAJ, 'l', OPTIONAL, "Sky: stepsize in Alpha (in radians or h:m:s)");
  XLALRegisterUvarMember( 	dDelta, 	DECJ, 'g', OPTIONAL, "Sky: stepsize in Delta (in radians or d:m:s)");
  XLALRegisterUvarMember(  	dFreq,          REAL8, 'r', OPTIONAL, "Stepsize for frequency in Hz");
  XLALRegisterUvarMember( 	df1dot, 	REAL8, 'e', OPTIONAL, "Stepsize for f1dot in Hz/s");
  XLALRegisterUvarMember( 	df2dot, 	 REAL8, 0 , OPTIONAL, "Stepsize for f2dot in Hz/s^2");
  XLALRegisterUvarMember( 	df3dot, 	 REAL8, 0 , OPTIONAL, "Stepsize for f3dot in Hz/s^3");

  XLALRegisterUvarMember( 	orbitasini, 	 REAL8, 0,  OPTIONAL, "Binary Orbit: Projected semi-major axis in light-seconds [Default: 0.0]");
  XLALRegisterUvarMember( 	orbitPeriod, 	 REAL8, 0,  OPTIONAL, "Binary Orbit: Period in seconds");
  XLALRegisterUvarMember(	orbitTp, 	 EPOCH, 0,  OPTIONAL, "Binary Orbit: (true) epoch of periapsis: use 'xx.yy[GPS|MJD]' format.");
  XLALRegisterUvarMember( 	orbitArgp, 	 REAL8, 0,  OPTIONAL, "Binary Orbit: Orbital argument of periapse in radians");
  XLALRegisterUvarMember( 	orbitEcc, 	 REAL8, 0,  OPTIONAL, "Binary Orbit: Orbital eccentricity");

  XLALRegisterUvarMember( 	orbitasiniBand,	 REAL8, 0,  OPTIONAL, "Binary Orbit: Band in Projected semi-major axis in light-seconds [Default: 0.0]");
  XLALRegisterUvarMember( 	orbitPeriodBand, REAL8, 0,  OPTIONAL, "Binary Orbit: Band in Period in seconds");
  XLALRegisterUvarMember(	orbitTpBand, 	 REAL8, 0,  OPTIONAL, "Binary Orbit: Band in (true) epoch of periapsis: use 'xx.yy[GPS|MJD]' format.");
  XLALRegisterUvarMember( 	orbitArgpBand, 	 REAL8, 0,  OPTIONAL, "Binary Orbit: Band in Orbital argument of periapse in radians");
  XLALRegisterUvarMember( 	orbitEccBand, 	 REAL8, 0,  OPTIONAL, "Binary Orbit: Band in Orbital eccentricity");

  XLALRegisterUvarMember( 	dorbitasini, 	 REAL8, 0,  OPTIONAL, "Binary Orbit: Spacing in Projected semi-major axis in light-seconds [Default: 0.0]");
  XLALRegisterUvarMember( 	dorbitPeriod, 	 REAL8, 0,  OPTIONAL, "Binary Orbit: Spacing in Period in seconds");
  XLALRegisterUvarMember(	dorbitTp, 	 REAL8, 0,  OPTIONAL, "Binary Orbit: Spacing in (true) epoch of periapsis: use 'xx.yy[GPS|MJD]' format.");
  XLALRegisterUvarMember( 	dorbitArgp, 	 REAL8, 0,  OPTIONAL, "Binary Orbit: Spacing in Orbital argument of periapse in radians");
  XLALRegisterUvarMember( 	dorbitEcc, 	 REAL8, 0,  OPTIONAL, "Binary Orbit: Spacing in Orbital eccentricity");

  XLALRegisterUvarMember(DataFiles, 	STRING, 'D', OPTIONAL, "File-pattern specifying (also multi-IFO) input SFT-files. Possibilities are:\n"
                         " - '<SFT file>;<SFT file>;...', where <SFT file> may contain wildcards\n - 'list:<file containing list of SFT files>'");

  XLALRegisterUvarMember( assumeSqrtSX,	 STRINGVector, 0,  OPTIONAL, "Don't estimate noise-floors but assume (stationary) per-IFO sqrt{SX} (if single value: use for all IFOs).\nNote that, unlike the historic --SignalOnly flag, this option will not lead to explicitly adding a +4 'correction' for noiseless SFTs to the output F-statistic.");

  XLALRegisterUvarMember( 	TwoFthreshold,	REAL8, 'F', OPTIONAL, "Set the threshold for selection of 2F");
  XLALRegisterUvarMember( 	gridType,	 INT4, 0 , OPTIONAL, "Grid: 0=flat, 1=isotropic, 2=metric, 3=skygrid-file, 6=grid-file, 8=spin-square, 9=spin-age-brk");
  XLALRegisterUvarMember( 	metricType,	INT4, 'M', OPTIONAL, "Metric: 0=none,1=Ptole-analytic,2=Ptole-numeric, 3=exact");
  XLALRegisterUvarMember( 	metricMismatch,	REAL8, 'X', OPTIONAL, "Maximal allowed mismatch for metric tiling");
  XLALRegisterUvarMember(outputLogfile,	 STRING, 0,  OPTIONAL, "Name of log-file identifying the code + search performed");
  XLALRegisterUvarMember(gridFile,	 STRING, 0,  OPTIONAL, "Load grid from this file: sky-grid or full-grid depending on --gridType.");
  XLALRegisterUvarMember(refTime,	 	 EPOCH, 0,  OPTIONAL, "Reference SSB epoch for pulsar-parameters: use 'xx.yy[GPS|MJD]' format [Default: startTime]");

  XLALRegisterUvarMember(outputFstat,	 STRING, 0,  OPTIONAL, "Output-file for F-statistic field over the parameter-space");
  XLALRegisterUvarMember(outputLoudest,	 STRING, 0,  OPTIONAL, "Loudest F-statistic candidate + estimated MLE amplitudes");

  XLALRegisterUvarMember(outputFstatHist, STRING, 0,  OPTIONAL, "Output-file for a discrete histogram of all Fstatistic values");
  XLALRegisterUvarMember(  FstatHistBin,    REAL8, 0,  OPTIONAL, "Width of an Fstatistic histogram bin");

  XLALRegisterUvarMember(  NumCandidatesToKeep,INT4, 0, OPTIONAL, "Number of Fstat 'candidates' to keep. (0 = All)");
  XLALRegisterUvarMember(FracCandidatesToKeep,REAL8, 0, OPTIONAL, "Fraction of Fstat 'candidates' to keep.");
  XLALRegisterUvarMember(   clusterOnScanline, INT4, 0, OPTIONAL, "Neighbors on each side for finding 1D local maxima on scanline");

  XLALRegisterUvarMember(minStartTime, 	 EPOCH, 0,  OPTIONAL, "Only use SFTs with timestamps starting from (including) this epoch (format 'xx.yy[GPS|MJD]') ");
  XLALRegisterUvarMember(maxStartTime, 	 EPOCH, 0,  OPTIONAL, "Only use SFTs with timestamps up to (excluding) this epoch (format 'xx.yy[GPS|MJD]')");

  XLALRegisterUvarMember(outputFstatAtoms,STRING, 0,  OPTIONAL, "Output filename *base* for F-statistic 'atoms' {a,b,Fa,Fb}_alpha. One file per doppler-point.");
  XLALRegisterUvarMember(  outputSingleFstats,BOOLEAN, 0,  OPTIONAL, "In multi-detector case, also output single-detector F-stats?");
  XLALRegisterUvarMember(RankingStatistic,STRING, 0,  DEVELOPER, "Rank toplist candidates according to 'F' or 'BSGL' statistic");

  // ----- Line robust stats parameters ----------
  XLALRegisterUvarMember(  computeBSGL,	BOOLEAN, 0,  OPTIONAL, "Compute and output line-robust statistic BSGL ");
  XLALRegisterUvarMember(  Fstar0,		REAL8, 0,  OPTIONAL, "BSGL: transition-scale parameter 'Fstar0'");
  XLALRegisterUvarMember(  oLGX,		STRINGVector, 0,  OPTIONAL, "BSGL: prior per-detector line-vs-Gauss odds 'oLGX' (Defaults to oLGX=1/Ndet)");
  XLALRegisterUvarMember(  BSGLlogcorr,	BOOLEAN, 0,  DEVELOPER,"BSGL: include log-correction terms (slower) or not (faster)");
  XLALRegisterUvarMember( 	BSGLthreshold,	REAL8, 0,  OPTIONAL, "BSGL threshold for candidate output");
  // --------------------------------------------

  XLALRegisterUvarMember(outputTransientStats,STRING, 0,  OPTIONAL, "TransientCW: Output filename for transient-CW statistics.");
  XLALRegisterUvarMember(outputTransientStatsAll,STRING, 0,  DEVELOPER, "TransientCW: Output filename for transient-CW statistics -- including the whole (t0,tau) grid for each candidate -- WARNING: CAN BE HUGE!.");
  XLALRegisterUvarMember( transient_WindowType,STRING, 0,OPTIONAL,  "TransientCW: Type of transient signal window to use. ('none', 'rect', 'exp').");
  XLALRegisterUvarMember( transient_t0Epoch, 	  	 EPOCH, 0, OPTIONAL, 	 "TransientCW: Earliest start-time for transient window search, in seconds (format 'xx.yy[GPS|MJD]')");
  XLALRegisterUvarMember( transient_t0Offset, 	  	 UINT4, 0, OPTIONAL, 	 "TransientCW: Earliest start-time for transient window search, as offset in seconds from dataStartGPS");
  XLALRegisterUvarMember( transient_t0Band, 	  	 UINT4, 0, OPTIONAL, 	 "TransientCW: Range of GPS start-times to search in transient search, in seconds");
  XLALRegisterUvarMember( transient_dt0, 	  	  	 INT4,  0, OPTIONAL, 	 "TransientCW: Step-size in transient-CW start-time in seconds [Default:Tsft]");
  XLALRegisterUvarMember( transient_tau, 	  	  	 UINT4, 0, OPTIONAL, 	 "TransientCW: Minimal transient-CW duration timescale, in seconds");
  XLALRegisterUvarMember( transient_tauBand, 	  	 UINT4, 0, OPTIONAL, 	 "TransientCW: Range of transient-CW duration timescales to search, in seconds");
  XLALRegisterUvarMember( transient_dtau, 	  	  	 INT4,  0, OPTIONAL, 	 "TransientCW: Step-size in transient-CW duration timescale, in seconds [Default:Tsft]");

  XLALRegisterUvarAuxDataMember( FstatMethod, UserEnum, XLALFstatMethodChoices(), 0, OPTIONAL,  "F-statistic method to use" );

  XLALRegisterUvarMember( 	countTemplates,  BOOLEAN, 0,  OPTIONAL, "Count number of templates (if supported) instead of search");
  XLALRegisterUvarMember( 	outputGrid,	 STRING, 0,  OPTIONAL, "Output-file for parameter-space grid (without running a search!)");

  /* ----- more experimental/expert options ----- */
  XLALRegisterUvarMember( 	UseNoiseWeights,BOOLEAN, 'W', DEVELOPER, "Use per-SFT noise weights");
  XLALRegisterUvarMember(ephemEarth, 	 STRING, 0,  DEVELOPER, "Earth ephemeris file to use");
  XLALRegisterUvarMember(ephemSun, 	 STRING, 0,  DEVELOPER, "Sun ephemeris file to use");

  XLALRegisterUvarAuxDataMember( SSBprecision, UserEnum, &SSBprecisionChoices, 0, DEVELOPER, "Precision to use for time-transformation to SSB" );

  XLALRegisterUvarMember( 	RngMedWindow,	INT4, 'k', DEVELOPER, "Running-Median window size");
  XLALRegisterUvarMember(	Dterms,		INT4, 't', DEVELOPER, "Number of kernel terms (single-sided) to use in\na) Dirichlet kernel if FstatMethod=Demod*\nb) sinc-interpolation kernel if FstatMethod=Resamp*" );

  XLALRegisterUvarMember(workingDir,     STRING, 'w', DEVELOPER, "Directory to use as work directory.");
  XLALRegisterUvarMember( 	timerCount, 	 REAL8, 0,  DEVELOPER, "N: Output progress/timer info every N seconds");

  XLALRegisterUvarMember( 	projectMetric, 	 BOOLEAN, 0,  DEVELOPER, "Use projected metric on Freq=const subspact");

  XLALRegisterUvarMember( strictSpindownBounds, BOOLEAN, 0, DEVELOPER, "suppress spindown grid points outside the [fkdot,fkdotBand] ranges? (only supported for --gridType=8)");

  XLALRegisterUvarMember(  spindownAge,     REAL8, 0,  DEVELOPER, "Spindown age for --gridType=9");
  XLALRegisterUvarMember(  minBraking,      REAL8, 0,  DEVELOPER, "Minimum braking index for --gridType=9");
  XLALRegisterUvarMember(  maxBraking,      REAL8, 0,  DEVELOPER, "Maximum braking index for --gridType=9");

  XLALRegisterUvarMember(transient_useFReg,   	 BOOLEAN, 0,  DEVELOPER, "FALSE: use 'standard' e^F for marginalization, if TRUE: use e^FReg = (1/D)*e^F (BAD)");

  XLALRegisterUvarMember(outputTiming,         STRING, 0,  DEVELOPER, "Append timing measurements and parameters into this file");
  XLALRegisterUvarMember(outputFstatTiming,    STRING, 0,  DEVELOPER, "Append F-statistic timing measurements and parameters into this file");

  XLALRegisterUvarMember(resampFFTPowerOf2,  BOOLEAN, 0,  DEVELOPER, "For Resampling methods: enforce FFT length to be a power of two (by rounding up)" );

  XLALRegisterUvarMember(allowedMismatchFromSFTLength, REAL8, 0, DEVELOPER, "Maximum allowed mismatch from SFTs being too long [Default: what's hardcoded in XLALFstatMaximumSFTLength]" );

  /* inject signals into the data being analyzed */
  XLALRegisterUvarMember(injectionSources,  STRINGVector, 0, DEVELOPER, "%s", InjectionSourcesHelpString );
  XLALRegisterUvarMember(injectSqrtSX,	    STRINGVector, 0, DEVELOPER, "Generate Gaussian Noise SFTs on-the-fly: CSV list of detectors' noise-floors sqrt{Sn}");
  XLALRegisterUvarMember(IFOs,	            STRINGVector, 0, DEVELOPER, "CSV list of detectors, eg. \"H1,H2,L1,G1, ...\", when no SFT files are specified");
  XLALRegisterUvarMember(timestampsFiles,   STRINGVector, 0, DEVELOPER, 
    "Files containing timestamps for the generated SFTs when using "UVAR_STR( injectSqrtSX ) ". "
    "Arguments correspond to the detectors given by " UVAR_STR( IFOs )
    "; for example, if " UVAR_STR( IFOs ) " is set to 'H1,L1', then an argument of "
    "'t1.txt,t2.txt' to this option will read H1 timestamps from the file 't1.txt', and L1 timestamps from the file 't2.txt'. "
    "The timebase of the generated SFTs is specified by " UVAR_STR( Tsft ) ". "
  );
  XLALRegisterUvarMember(Tsft,              REAL8, 0, DEVELOPER, "Generate SFTs with this timebase (in seconds) instead of loading from files. Requires --injectSqrtSX, --IFOs, --timestampsFiles");
  XLALRegisterUvarMember(randSeed,          INT4, 0, DEVELOPER, "Specify random-number seed for reproducible noise (0 means use /dev/urandom for seeding).");

  return XLAL_SUCCESS;

} /* initUserVars() */

/** Initialized Fstat-code: handle user-input and set everything up.
 * NOTE: the logical *order* of things in here is very important, so be careful
 */
int
InitFstat ( ConfigVariables *cfg, const UserInput_t *uvar )
{
  XLAL_CHECK ( (cfg != NULL) && (uvar != NULL), XLAL_EINVAL );

  REAL8 fCoverMin, fCoverMax;	/* covering frequency-band to read from SFTs */
  SFTCatalog *catalog = NULL;
  SFTConstraints XLAL_INIT_DECL(constraints);
  LIGOTimeGPS endTime;
  size_t toplist_length = uvar->NumCandidatesToKeep;

  /* check compatibility of some output options */
  cfg->runSearch = !( uvar->countTemplates || uvar->outputGrid );
  XLAL_CHECK ( cfg->runSearch || !( uvar->outputFstat || uvar->outputFstatAtoms || uvar->outputFstatHist || uvar->outputFstatTiming || uvar->outputLoudest || uvar->outputTransientStats || uvar->outputTransientStatsAll ), XLAL_EINVAL, "Invalid input: --countTemplates or --outputGrid means that no search is actually run, so we can't do any of --outputFstat --outputFstatAtoms --outputFstatHist --outputFstatTiming --outputLoudest --outputTransientStats --outputTransientStatsAll" );

  /* set the current working directory */
  XLAL_CHECK ( chdir ( uvar->workingDir ) == 0, XLAL_EINVAL, "Unable to change directory to workinDir '%s'\n", uvar->workingDir );

  /* ----- set computational parameters for F-statistic from User-input ----- */
  cfg->useResamp = ( uvar->FstatMethod >= FMETHOD_RESAMP_GENERIC ); // use resampling;

  /* check that resampling is compatible with gridType */
  if ( cfg->useResamp && uvar->gridType > GRID_SKY_LAST /* end-marker for factored grid types */ ) {
    XLAL_ERROR ( XLAL_EINVAL, "\nUse of resampling FstatMethod is incompatible with non-factored gridType=%i\n\n", uvar->gridType );
  }

  /* if IFO string vector was passed by user, parse it for later use */
  if ( uvar->IFOs != NULL ) {
    XLAL_CHECK ( XLALParseMultiLALDetector ( &(cfg->multiIFO), uvar->IFOs ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  /* read timestamps if requested by the user */
  if (uvar->timestampsFiles != NULL) {
      XLAL_CHECK ( (cfg->multiTimestamps = XLALReadMultiTimestampsFiles ( uvar->timestampsFiles )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( (cfg->multiTimestamps->length > 0) && (cfg->multiTimestamps->data != NULL), XLAL_EINVAL, "Got empty timestamps-list from XLALReadMultiTimestampsFiles()\n" );

      for ( UINT4 X=0; X < cfg->multiTimestamps->length; X ++ ) {
        cfg->multiTimestamps->data[X]->deltaT = uvar->Tsft;	// Tsft information not given by timestamps-file
      }
  }

  LIGOTimeGPS minStartTime = uvar->minStartTime;
  LIGOTimeGPS maxStartTime = uvar->maxStartTime;
  constraints.minStartTime = &minStartTime;
  constraints.maxStartTime = &maxStartTime;

  /* get full SFT-catalog of all matching (multi-IFO) SFTs */
  /* DataFiles optional because of injectSqrtSX option, don't try to load if no files given **/

  if( uvar->DataFiles != NULL ) {
    LogPrintf (LOG_NORMAL, "Finding all SFTs to load ... ");
    XLAL_CHECK ( (catalog = XLALSFTdataFind ( uvar->DataFiles, &constraints )) != NULL, XLAL_EFUNC );
    LogPrintfVerbatim (LOG_NORMAL, "done. (found %d SFTs)\n", catalog->length);
  } else {
    /* Build a fake catalog with timestamps and IFOs given on the commandline instead of noise data files */
    /* the data missing in the locators then signal to the Fstat code that fake noise needs to be generated, */
    /* the noise level is passed from the sqrtSX option */

    XLAL_CHECK ( (catalog = XLALMultiAddToFakeSFTCatalog ( NULL, uvar->IFOs,cfg-> multiTimestamps)) != NULL, XLAL_EFUNC );
  }

  XLAL_CHECK ( catalog->length > 0, XLAL_EINVAL, "\nSorry, didn't find any matching SFTs with pattern '%s'!\n\n", uvar->DataFiles );

  /* deduce start- and end-time of the observation spanned by the data */
  UINT4 numSFTfiles = catalog->length;

  cfg->Tsft = 1.0 / catalog->data[0].header.deltaF;
  cfg->startTime = catalog->data[0].header.epoch;
  endTime   = catalog->data[numSFTfiles - 1].header.epoch;
  XLALGPSAdd ( &endTime, cfg->Tsft );	/* add on Tsft to last SFT start-time */

  // time spanned by the SFTs
  cfg->Tspan = XLALGPSDiff ( &endTime, &cfg->startTime );

  /* ----- load ephemeris-data ----- */
  XLAL_CHECK ( (cfg->ephemeris = XLALInitBarycenter( uvar->ephemEarth, uvar->ephemSun )) != NULL, XLAL_EFUNC );

  /* ----- get reference-times (from user if given, use startTime otherwise): ----- */
  LIGOTimeGPS refTime;
  if ( XLALUserVarWasSet(&uvar->refTime) )
    {
      refTime = uvar->refTime;
    }
  else
    {
      refTime = cfg->startTime;
    }

  /* define sky position variables from user input */
  cfg->Alpha = uvar->Alpha;
  cfg->Delta = uvar->Delta;


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
    BOOLEAN haveAlphaDelta = (XLALUserVarWasSet(&uvar->Alpha) && XLALUserVarWasSet(&uvar->Delta));

    if (uvar->skyRegion)
      {
	XLAL_CHECK ( (cfg->searchRegion.skyRegionString = XLALCalloc(1, strlen(uvar->skyRegion)+1)) != NULL, XLAL_ENOMEM );
	strcpy (cfg->searchRegion.skyRegionString, uvar->skyRegion);
      }
    else if (haveAlphaDelta)    /* parse this into a sky-region */
      {
	XLAL_CHECK ( (cfg->searchRegion.skyRegionString = XLALSkySquare2String ( cfg->Alpha, cfg->Delta, uvar->AlphaBand, uvar->DeltaBand)) != NULL, XLAL_EFUNC );
      }
    else if (!uvar->gridFile)
      {
	XLAL_ERROR ( XLAL_EINVAL, "\nCould not setup searchRegion, have neither skyRegion nor (Alpha, Delta) nor a gridFile!\n\n" );
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

  /*read signal parameters to be injected, if requested by the user*/
  if ( uvar->injectionSources != NULL ) {
    XLAL_CHECK_MAIN ( (cfg->injectionSources = XLALPulsarParamsFromUserInput ( uvar->injectionSources, NULL ) ) != NULL, XLAL_EFUNC );
  }

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
    if ( cfg->useResamp ) {	// in the resampling-case, temporarily take out frequency-dimension of the Doppler template bank
      scanInit.searchRegion.fkdotBand[0] = 0;
    }

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
    const LALDetector *detector;
    XLAL_CHECK ( (detector = XLALGetSiteInfo ( catalog->data[0].header.name ) ) != NULL, XLAL_EFUNC );
    scanInit.Detector  = detector;

    /* Specific options for some gridTypes */
    if (uvar->gridType == GRID_SPINDOWN_SQUARE) {
      scanInit.extraArgs[0] = !uvar->strictSpindownBounds;
    }
    else if (uvar->gridType == GRID_SPINDOWN_AGEBRK) {
      scanInit.extraArgs[0] = uvar->spindownAge;
      scanInit.extraArgs[1] = uvar->minBraking;
      scanInit.extraArgs[2] = uvar->maxBraking;
    }

    LogPrintf (LOG_NORMAL, "Setting up template grid ... ");
    XLAL_CHECK ( (cfg->scanState = XLALInitDopplerFullScan ( &scanInit)) != NULL, XLAL_EFUNC );
    LogPrintfVerbatim (LOG_NORMAL, "template grid ready.\n");
    XLALNumDopplerTemplates ( cfg->scanState );
  }

  /* maximum ranges of binary orbit parameters */
  const REAL8 binaryMaxAsini = MYMAX( uvar->orbitasini, uvar->orbitasini + uvar->orbitasiniBand );
  const REAL8 binaryMinPeriod = MYMIN( uvar->orbitPeriod, uvar->orbitPeriod + uvar->orbitPeriodBand );
  const REAL8 binaryMaxEcc = MYMAX( uvar->orbitEcc, uvar->orbitEcc + uvar->orbitEccBand );

  { /* ----- What frequency-band do we need to read from the SFTs?
     * propagate spin-range from refTime to startTime and endTime of observation
     */
    PulsarSpinRange spinRangeRef;	/* temporary only */

    // extract spanned spin-range at reference-time from the template-bank
    XLAL_CHECK ( XLALGetDopplerSpinRange ( &spinRangeRef, cfg->scanState ) == XLAL_SUCCESS, XLAL_EFUNC );

    // in the resampling case, we need to restore the frequency-band info now, which we set to 0
    // before calling the DopplerInit template bank construction
    if ( cfg->useResamp ) {
      spinRangeRef.fkdotBand[0] = tmpFreqBandRef;
    }

    // write this search spin-range@refTime back into 'cfg' struct for users' reference
    cfg->searchRegion.refTime = spinRangeRef.refTime;
    memcpy ( cfg->searchRegion.fkdot, spinRangeRef.fkdot, sizeof(cfg->searchRegion.fkdot) );
    memcpy ( cfg->searchRegion.fkdotBand, spinRangeRef.fkdotBand, sizeof(cfg->searchRegion.fkdotBand) );

    // Compute covering frequency range, accounting for Doppler modulation due to the Earth and any binary companion */
    XLALCWSignalCoveringBand( &fCoverMin, &fCoverMax, &cfg->startTime, &endTime, &spinRangeRef, binaryMaxAsini, binaryMinPeriod, binaryMaxEcc );

  } /* extrapolate spin-range */

  /* check that SFT length is within allowed maximum */
  /* use fCoverMax here to work with loading grid from file (--gridType=6) */
  XLAL_CHECK ( XLALFstatCheckSFTLengthMismatch ( cfg->Tsft, fCoverMax, binaryMaxAsini, binaryMinPeriod, uvar->allowedMismatchFromSFTLength ) == XLAL_SUCCESS, XLAL_EFUNC, "Excessive mismatch would be incurred due to SFTs being too long for the current search setup. Please double-check your parameter ranges or provide shorter input SFTs. If you really know what you're doing, you could also consider using the --allowedMismatchFromSFTLength override." );

  /* If requested, assume a certain PSD instead of estimating it from data.
   * NOTE that, unlike for the now removed --SignalOnly flag,
   * no extra +4 is added to the F-statistic output in this case.
   */
  MultiNoiseFloor s_assumeSqrtSX, *assumeSqrtSX;
  if ( uvar->assumeSqrtSX != NULL ) {
    XLAL_CHECK( XLALParseMultiNoiseFloor( &s_assumeSqrtSX, uvar->assumeSqrtSX, XLALCountIFOsInCatalog(catalog) ) == XLAL_SUCCESS, XLAL_EFUNC );
    assumeSqrtSX = &s_assumeSqrtSX;
  } else {
    assumeSqrtSX = NULL;
  }

  if ( XLALUserVarWasSet ( &uvar->dFreq) ) {
    cfg->dFreq = uvar->dFreq;
  } else {
    cfg->dFreq = 1.0/(2*cfg->Tspan);
  }
  if ( cfg->useResamp )	{ // handle special resampling case, where we deal with a vector of F-stat values instead of one
    cfg->numFreqBins_FBand = (UINT4) ( 1 + floor ( cfg->searchRegion.fkdotBand[0] / cfg->dFreq ) );
  } else {
    cfg->numFreqBins_FBand = 1;	// number of frequency-bins in the frequency-band used for resampling (1 if not using Resampling)
  }

  MultiNoiseFloor *injectSqrtSX = NULL;
  MultiNoiseFloor s_injectSqrtSX;

  if(uvar->injectSqrtSX != NULL) {
    XLAL_CHECK( XLALParseMultiNoiseFloor( &s_injectSqrtSX, uvar->injectSqrtSX, XLALCountIFOsInCatalog(catalog) ) == XLAL_SUCCESS, XLAL_EFUNC );
    injectSqrtSX = &s_injectSqrtSX;
  }

  FstatOptionalArgs optionalArgs = FstatOptionalArgsDefaults;
  optionalArgs.Dterms  = uvar->Dterms;
  optionalArgs.SSBprec = uvar->SSBprecision;
  optionalArgs.runningMedianWindow = uvar->RngMedWindow;
  optionalArgs.injectSources = cfg->injectionSources;
  optionalArgs.injectSqrtSX = injectSqrtSX;
  optionalArgs.randSeed=uvar->randSeed;
  optionalArgs.assumeSqrtSX = assumeSqrtSX;
  optionalArgs.FstatMethod = uvar->FstatMethod;
  optionalArgs.resampFFTPowerOf2 = uvar->resampFFTPowerOf2;
  optionalArgs.collectTiming = XLALUserVarWasSet ( &uvar->outputFstatTiming );
  optionalArgs.allowedMismatchFromSFTLength = uvar->allowedMismatchFromSFTLength;


  XLAL_CHECK ( (cfg->Fstat_in = XLALCreateFstatInput( catalog, fCoverMin, fCoverMax, cfg->dFreq, cfg->ephemeris, &optionalArgs )) != NULL, XLAL_EFUNC );
  XLALDestroySFTCatalog(catalog);

  cfg->Fstat_what = FSTATQ_2F;   // always calculate multi-detector 2F
  if ( XLALUserVarWasSet( &uvar->outputLoudest ) ) {
    cfg->Fstat_what |= FSTATQ_FAFB;   // also calculate Fa,b parts for parameter estimation
  }

  /* get SFT detectors and timestamps */
  const MultiLALDetector *multiIFO;
  XLAL_CHECK ( (multiIFO = XLALGetFstatInputDetectors( cfg->Fstat_in )) != NULL, XLAL_EFUNC );
  const MultiLIGOTimeGPSVector *multiTS;
  XLAL_CHECK ( (multiTS = XLALGetFstatInputTimestamps( cfg->Fstat_in )) != NULL, XLAL_EFUNC );

  /* count total number of SFTs loaded */
  cfg->NSFTs = 0;
  for ( UINT4 X = 0; X < multiTS->length; X++ ) {
    cfg->NSFTs += multiTS->data[X]->length;
  }

  /* for column headings string, get number of detectors, detector name vector, and SFTs per detector vector */
  {
    const UINT4 numDetectors = multiIFO->length;
    XLAL_CHECK ( (cfg->numSFTsPerDet = XLALCreateUINT4Vector( numDetectors )) != NULL, XLAL_EFUNC );
    cfg->detectorIDs = NULL;
    for (UINT4 X = 0; X < numDetectors; X++) {
      cfg->numSFTsPerDet->data[X] = multiTS->data[X]->length;
      XLAL_CHECK ( (cfg->detectorIDs = XLALAppendString2Vector ( cfg->detectorIDs, multiIFO->sites[X].frDetector.prefix )) != NULL, XLAL_EFUNC );
    } /* for X < numDetectors */
  }

  /* ----- set up scanline-window if requested for 1D local-maximum clustering on scanline ----- */
  XLAL_CHECK ( (cfg->scanlineWindow = XLALCreateScanlineWindow ( uvar->clusterOnScanline )) != NULL, XLAL_EFUNC );

  /* set number of toplist candidates from fraction if asked to */
  if (0.0 < uvar->FracCandidatesToKeep && uvar->FracCandidatesToKeep <= 1.0) {
    if (XLALNumDopplerTemplates(cfg->scanState) <= 0.0) {
      XLAL_ERROR ( XLAL_EINVAL , "Cannot use FracCandidatesToKeep because number of templates was counted to be zero!\n");
    }
    toplist_length = ceil(XLALNumDopplerTemplates(cfg->scanState) * uvar->FracCandidatesToKeep);
  }

  /* ----- set up toplist if requested ----- */
  if ( toplist_length > 0 ) {
    if ( strcmp(uvar->RankingStatistic, "F") == 0 )
     cfg->RankingStatistic = RANKBY_2F;
    else if ( strcmp(uvar->RankingStatistic, "BSGL") == 0 )
      {
        if ( !uvar->computeBSGL ) {
          XLAL_ERROR ( XLAL_EINVAL, "\nERROR: Ranking by BSGL only possible if --computeBSGL given.\n\n");
        }
        cfg->RankingStatistic = RANKBY_BSGL;
      }
    else
      {
        XLAL_ERROR ( XLAL_EINVAL, "\nERROR: Invalid value specified for candidate ranking - supported are 'F' and 'BSGL'.\n\n");
      }

    if ( cfg->RankingStatistic == RANKBY_BSGL )
      {
        XLAL_CHECK ( create_toplist( &(cfg->FstatToplist), toplist_length, sizeof(FstatCandidate), compareFstatCandidates_BSGL) == 0, XLAL_EFUNC );
      }
    else // rank by F-stat
      {
        XLAL_CHECK ( create_toplist( &(cfg->FstatToplist), toplist_length, sizeof(FstatCandidate), compareFstatCandidates) == 0, XLAL_EFUNC );
      }
  } /* if toplist_length > 0 */


  /* ----- transient-window related parameters ----- */
  int twtype;
  XLAL_CHECK ( (twtype = XLALParseTransientWindowName ( uvar->transient_WindowType )) >= 0, XLAL_EFUNC );
  cfg->transientWindowRange.type = twtype;

  /* make sure user doesn't set window=none but sets window-parameters => indicates she didn't mean 'none' */
  XLAL_CHECK ( !( ( cfg->transientWindowRange.type == TRANSIENT_NONE )
                  && ( XLALUserVarWasSet ( &uvar->transient_dt0 )
                       || XLALUserVarWasSet ( &uvar->transient_t0Epoch )
                       || XLALUserVarWasSet ( &uvar->transient_t0Offset )
                       || XLALUserVarWasSet ( &uvar->transient_t0Band )
                       || XLALUserVarWasSet ( &uvar->transient_tau )
                       || XLALUserVarWasSet ( &uvar->transient_tauBand )
                       || XLALUserVarWasSet ( &uvar->transient_dtau ) )
                ), XLAL_EINVAL, "ERROR: transientWindow->type == NONE, but window-parameters were set! Use a different window-type!" );

  if ( XLALUserVarWasSet ( &uvar->transient_t0Offset ) ) {
    XLAL_CHECK ( !XLALUserVarWasSet ( &uvar->transient_t0Epoch ), XLAL_EINVAL, "ERROR: only one of transient_t0Epoch, transient_t0Offset may be used!" );
    cfg->transientWindowRange.t0      = cfg->startTime.gpsSeconds + uvar->transient_t0Offset;
  }
  else if ( XLALUserVarWasSet ( &uvar->transient_t0Epoch ) ) {
    cfg->transientWindowRange.t0      = uvar->transient_t0Epoch.gpsSeconds; /* just dropping ns part here */
  }

  if ( XLALUserVarWasSet ( &uvar->transient_t0Band ) ) {
    cfg->transientWindowRange.t0Band  = uvar->transient_t0Band;
  }

  if ( XLALUserVarWasSet ( &uvar->transient_dt0 ) ) {
    cfg->transientWindowRange.dt0 = uvar->transient_dt0;
  }
  else {
    cfg->transientWindowRange.dt0 = cfg->Tsft;
  }

  if ( XLALUserVarWasSet ( &uvar->transient_tau ) ) {
    cfg->transientWindowRange.tau  = uvar->transient_tau;
  }

  if ( XLALUserVarWasSet ( &uvar->transient_tauBand ) ) {
    cfg->transientWindowRange.tauBand  = uvar->transient_tauBand;
  }

  if ( XLALUserVarWasSet ( &uvar->transient_dtau ) ) {
    cfg->transientWindowRange.dtau = uvar->transient_dtau;
  }
  else {
    cfg->transientWindowRange.dtau = cfg->Tsft;
  }

  /* get atoms back from Fstat-computing, either if atoms-output or transient-Bstat output was requested */
  if ( ( uvar->outputFstatAtoms != NULL ) || ( uvar->outputTransientStats != NULL ) || ( uvar->outputTransientStatsAll != NULL ) ) {
    cfg->Fstat_what |= FSTATQ_ATOMS_PER_DET;
  }

  /* return single-IFO Fstat values for Line-robust statistic */
  if ( uvar->outputSingleFstats || uvar->computeBSGL ) {
    cfg->Fstat_what |= FSTATQ_2F_PER_DET;
  }

  /* ---------- prepare Line-robust statistics parameters ---------- */
  if ( uvar->computeBSGL )
    {
      UINT4 numDetectors = multiIFO->length;
      /* take BSGL user input and pre-compute the corresponding BSGLsetup */
      REAL4 *oLGX_p = NULL;
      REAL4 oLGX[PULSAR_MAX_DETECTORS];
      if ( uvar->oLGX != NULL )
        {
          if ( uvar->oLGX->length != numDetectors ) {
            XLAL_ERROR ( XLAL_EINVAL, "Invalid input: length(oLGX) = %d differs from number of detectors (%d)'\n", uvar->oLGX->length, numDetectors );
          }
          XLAL_CHECK ( XLALParseLinePriors ( &oLGX[0], uvar->oLGX ) == XLAL_SUCCESS, XLAL_EFUNC );
          oLGX_p = &oLGX[0];
        } // if uvar->oLGX != NULL

      XLAL_CHECK ( (cfg->BSGLsetup = XLALCreateBSGLSetup ( numDetectors, uvar->Fstar0, oLGX_p, uvar->BSGLlogcorr, 1 )) != NULL, XLAL_EFUNC ); // coherent F-stat: NSeg=1
    } // if uvar_computeBSGL

  return XLAL_SUCCESS;

} /* InitFstat() */

/**
 * Produce a log-string describing the present run-setup
 */
CHAR *
XLALGetLogString ( const ConfigVariables *cfg, const BOOLEAN verbose )
{
  XLAL_CHECK_NULL ( cfg != NULL, XLAL_EINVAL );

#define BUFLEN 1024
  CHAR buf[BUFLEN];

  CHAR *logstr = NULL;
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, "%% cmdline: " )) != NULL, XLAL_EFUNC );
  CHAR *cmdline;
  XLAL_CHECK_NULL ( (cmdline = XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, cmdline )) != NULL, XLAL_EFUNC );
  XLALFree ( cmdline );
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, "\n" )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, cfg->VCSInfoString )) != NULL, XLAL_EFUNC );

  if ( !verbose ) {
      /* don't include search details */
      return logstr;
  }

  XLAL_CHECK_NULL ( snprintf ( buf, BUFLEN, "%%%% FstatMethod used: '%s'\n", XLALGetFstatInputMethodName ( cfg->Fstat_in ) ) < BUFLEN, XLAL_EBADLEN );
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, buf )) != NULL, XLAL_EFUNC );

  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, "%% Started search: " )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, LogGetTimestamp() )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, "\n%% Loaded SFTs: [ " )) != NULL, XLAL_EFUNC );

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
    startTimeUTCString[strlen(startTimeUTCString)-1] = 0;	// kill trailing newline
    XLAL_CHECK_NULL ( snprintf ( buf, BUFLEN, "%%%% GPS starttime         = %d (%s GMT)\n", startTimeSeconds, startTimeUTCString ) < BUFLEN, XLAL_EBADLEN );
    XLALFree ( startTimeUTCString );
  }
  char bufGPS[32];
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, buf )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_NULL ( snprintf ( buf, BUFLEN, "%%%% Total time spanned    = %.0f s (%.2f hours)\n", cfg->Tspan, cfg->Tspan/3600.0 ) < BUFLEN, XLAL_EBADLEN );
  XLAL_CHECK_NULL ( (logstr = XLALStringAppend ( logstr, buf )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_NULL ( snprintf (buf, BUFLEN, "%%%% Pulsar-params refTime = %s\n", XLALGPSToStr ( bufGPS, &(cfg->searchRegion.refTime) )) < BUFLEN, XLAL_EBADLEN );
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
int
WriteFstatLog ( const CHAR *log_fname, const CHAR *log_string )
{
  XLAL_CHECK ( (log_fname != NULL) && (log_string != NULL), XLAL_EINVAL );

  /* prepare log-file for writing */
  FILE *fplog;
  XLAL_CHECK ( (fplog = fopen(log_fname, "wb" )) != NULL, XLAL_ESYS, "Failed to open log-file '%s' for writing.\n\n", log_fname );

  fprintf (fplog, "%%%% LOG-FILE for ComputeFstatistic run\n\n");
  fprintf (fplog, "%s", log_string);
  fclose (fplog);

  return XLAL_SUCCESS;

} /* WriteFstatLog() */


/** Free all globally allocated memory. */
void
Freemem( ConfigVariables *cfg )
{
  if ( !cfg ) {
    return;
  }
  XLALDestroyUINT4Vector ( cfg->numSFTsPerDet );
  XLALDestroyStringVector ( cfg->detectorIDs );

  XLALDestroyFstatInput ( cfg->Fstat_in );

  /* destroy FstatToplist if any */
  if ( cfg->FstatToplist ) {
    free_toplist( &(cfg->FstatToplist) );
  }

  if ( cfg->scanlineWindow ) {
    XLALDestroyScanlineWindow ( cfg->scanlineWindow );
  }

  /* Free config-Variables and userInput stuff */
  XLALDestroyUserVars();

  if ( cfg->searchRegion.skyRegionString ) {
    XLALFree ( cfg->searchRegion.skyRegionString );
  }

  /* Free ephemeris data */
  XLALDestroyEphemerisData ( cfg->ephemeris );

  if ( cfg->VCSInfoString ) {
    XLALFree ( cfg->VCSInfoString );
  }
  if ( cfg->logstring ) {
    XLALFree ( cfg->logstring );
  }

  XLALFree ( cfg->BSGLsetup );

  /* free source injection if any */
  if ( cfg->injectionSources ) {
    XLALDestroyPulsarParamsVector ( cfg->injectionSources );
  }

  if( cfg->multiTimestamps ) {
    XLALDestroyMultiTimestamps( cfg->multiTimestamps );
  }

  return;

} /* Freemem() */


/*----------------------------------------------------------------------*/
/**
 * Some general consistency-checks on user-input.
 * Throws an error plus prints error-message if problems are found.
 */
int
checkUserInputConsistency ( const UserInput_t *uvar )
{
  XLAL_CHECK ( uvar != NULL, XLAL_EINVAL );

  /* check for negative stepsizes in Freq, Alpha, Delta */
  if ( XLALUserVarWasSet(&uvar->dAlpha) && (uvar->dAlpha < 0) )
    {
      XLALPrintError ("\nNegative value of stepsize dAlpha not allowed!\n\n");
      XLAL_ERROR ( XLAL_EINVAL );
    }
  if ( XLALUserVarWasSet(&uvar->dDelta) && (uvar->dDelta < 0) )
    {
      XLALPrintError ("\nNegative value of stepsize dDelta not allowed!\n\n");
      XLAL_ERROR ( XLAL_EINVAL );
    }
  if ( XLALUserVarWasSet(&uvar->dFreq) && (uvar->dFreq < 0) )
    {
      XLALPrintError ("\nNegative value of stepsize dFreq not allowed!\n\n");
      XLAL_ERROR ( XLAL_EINVAL );
    }

  /* binary parameter checks */
  if ( XLALUserVarWasSet(&uvar->orbitPeriod) && (uvar->orbitPeriod <= 0) )
    {
      XLALPrintError ("\nNegative or zero value of orbital period not allowed!\n\n");
      XLAL_ERROR ( XLAL_EINVAL );
    }
  if ( XLALUserVarWasSet(&uvar->orbitasini) && (uvar->orbitasini < 0) )
    {
      XLALPrintError ("\nNegative value of projected orbital semi-major axis not allowed!\n\n");
      XLAL_ERROR ( XLAL_EINVAL );
    }
  if ( XLALUserVarWasSet(&uvar->orbitArgp) && ((uvar->orbitArgp < 0) || (uvar->orbitArgp >= LAL_TWOPI)) )
    {
      XLALPrintError ("\nOrbital argument of periapse must lie in range [0 2*PI)!\n\n");
      XLAL_ERROR ( XLAL_EINVAL );
    }
  if ( XLALUserVarWasSet(&uvar->orbitEcc) && (uvar->orbitEcc < 0) )
    {
      XLALPrintError ("\nNegative value of orbital eccentricity not allowed!\n\n");
      XLAL_ERROR ( XLAL_EINVAL );
    }

  /* grid-related checks */
  {
    BOOLEAN haveAlphaBand = XLALUserVarWasSet( &uvar->AlphaBand );
    BOOLEAN haveDeltaBand = XLALUserVarWasSet( &uvar->DeltaBand );
    BOOLEAN haveSkyRegion, haveAlphaDelta, haveGridFile;
    BOOLEAN useSkyGridFile, useFullGridFile, haveMetric, useMetric;

    haveSkyRegion  	= (uvar->skyRegion != NULL);
    haveAlphaDelta 	= (XLALUserVarWasSet(&uvar->Alpha) && XLALUserVarWasSet(&uvar->Delta));
    haveGridFile      	= (uvar->gridFile != NULL);
    useSkyGridFile   	= (uvar->gridType == GRID_FILE_SKYGRID);
    useFullGridFile	= (uvar->gridType == GRID_FILE_FULLGRID);
    haveMetric     	= (uvar->metricType > LAL_PMETRIC_NONE);
    useMetric     	= (uvar->gridType == GRID_METRIC);

    if ( !useFullGridFile && !useSkyGridFile && haveGridFile )
      {
        XLALPrintError ("\nERROR: gridFile was specified but not needed for gridType=%d\n\n", uvar->gridType );
        XLAL_ERROR ( XLAL_EINVAL );
      }
    if ( useSkyGridFile && !haveGridFile )
      {
        XLALPrintError ("\nERROR: gridType=SKY-FILE, but no --gridFile specified!\n\n");
        XLAL_ERROR ( XLAL_EINVAL );
      }
    if ( useFullGridFile && !haveGridFile )
      {
	XLALPrintError ("\nERROR: gridType=GRID-FILE, but no --gridFile specified!\n\n");
        XLAL_ERROR ( XLAL_EINVAL );
      }

    if ( (haveAlphaBand && !haveDeltaBand) || (haveDeltaBand && !haveAlphaBand) )
      {
	XLALPrintError ("\nERROR: Need either BOTH (AlphaBand, DeltaBand) or NONE.\n\n");
        XLAL_ERROR ( XLAL_EINVAL );
      }

    if ( haveSkyRegion && haveAlphaDelta )
      {
        XLALPrintError ("\nOverdetermined sky-region: only use EITHER (Alpha,Delta) OR skyRegion!\n\n");
        XLAL_ERROR ( XLAL_EINVAL );
      }
    if ( !haveSkyRegion && !haveAlphaDelta && !useSkyGridFile && !useFullGridFile )
      {
        XLALPrintError ("\nUnderdetermined sky-region: use one of (Alpha,Delta), skyRegion or a gridFile!\n\n");
        XLAL_ERROR ( XLAL_EINVAL );
      }

    if ( !useMetric && haveMetric)
      {
        XLALPrintWarning ("\nWARNING: Metric was specified for non-metric grid... will be ignored!\n");
      }
    if ( useMetric && !haveMetric)
      {
        XLALPrintError ("\nERROR: metric grid-type selected, but no metricType selected\n\n");
        XLAL_ERROR ( XLAL_EINVAL );
      }

    /* Specific checks for --gridType=GRID_SPINDOWN_{SQUARE,AGEBRK} parameter spaces */
    if (uvar->gridType == GRID_SPINDOWN_SQUARE || uvar->gridType == GRID_SPINDOWN_AGEBRK) {

      /* Check that no grid spacings were given */
      if (uvar->df1dot != 0.0 || uvar->df2dot != 0.0 || uvar->df3dot != 0.0) {
        XLALPrintError ("\nERROR: df{1,2,3}dot cannot be used with gridType={8,9}\n\n");
        XLAL_ERROR ( XLAL_EINVAL );
      }

    }

    if (uvar->strictSpindownBounds && ! (uvar->gridType == GRID_SPINDOWN_SQUARE) ) {
      XLALPrintError("\nERROR: strictSpindownBounds can only be used with gridType=8\n\n");
      XLAL_ERROR ( XLAL_EINVAL );
    }

    /* Specific checks for --gridType=GRID_SPINDOWN_AGEBRK parameter space */
    if (uvar->gridType == GRID_SPINDOWN_AGEBRK) {

      /* Check that no third spindown range were given */
      if (uvar->f3dot != 0.0 || uvar->f3dotBand != 0.0) {
        XLALPrintError ("\nERROR: f3dot and f3dotBand cannot be used with gridType=9\n\n");
        XLAL_ERROR ( XLAL_EINVAL );
      }

      /* Check age and braking indices */
      if (uvar->spindownAge <= 0.0) {
        XLALPrintError ("\nERROR: spindownAge must be strictly positive with gridType=9\n\n");
        XLAL_ERROR ( XLAL_EINVAL );
      }
      if (uvar->minBraking <= 0.0) {
        XLALPrintError ("\nERROR: minBraking must be strictly positive with gridType=9\n\n");
        XLAL_ERROR ( XLAL_EINVAL );
      }
      if (uvar->maxBraking <= 0.0) {
        XLALPrintError ("\nERROR: minBraking must be strictly positive with gridType=9\n\n");
        XLAL_ERROR ( XLAL_EINVAL );
      }
      if (uvar->minBraking >= uvar->maxBraking) {
        XLALPrintError ("\nERROR: minBraking must be strictly less than maxBraking with gridType=9\n\n");
        XLAL_ERROR ( XLAL_EINVAL );
      }

      /* Check that no first and second spindown ranges were given */
      if (uvar->f1dot != 0.0 || uvar->f1dotBand != 0.0) {
        XLALPrintError ("\nERROR: f1dot and f1dotBand cannot be used with gridType=9\n\n");
        XLAL_ERROR ( XLAL_EINVAL );
      }
      if (uvar->f2dot != 0.0 || uvar->f2dotBand != 0.0) {
        XLALPrintError ("\nERROR: f2dot and f2dotBand cannot be used with gridType=9\n\n");
        XLAL_ERROR ( XLAL_EINVAL );
      }

    }

  } /* Grid-related checks */

  /* check NumCandidatesToKeep and FracCandidatesToKeep */
  if (XLALUserVarWasSet(&uvar->NumCandidatesToKeep) && XLALUserVarWasSet(&uvar->FracCandidatesToKeep)) {
    XLALPrintError ("\nERROR: NumCandidatesToKeep and FracCandidatesToKeep are mutually exclusive\n\n");
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if (XLALUserVarWasSet(&uvar->FracCandidatesToKeep) && (uvar->FracCandidatesToKeep <= 0.0 || 1.0 < uvar->FracCandidatesToKeep)) {
    XLALPrintError ("\nERROR: FracCandidatesToKeep must be greater than 0.0 and less than or equal to 1.0\n\n");
    XLAL_ERROR ( XLAL_EINVAL );
  }

  /* check options for input data and injections */
  XLAL_CHECK ( (uvar->DataFiles ==NULL) ^ (uvar->injectSqrtSX == NULL), XLAL_EINVAL,  "Must pass exactly one out of --DataFiles or --injectSqrtSX \n");
  if(uvar->DataFiles !=NULL) {
     XLAL_CHECK ( !XLALUserVarWasSet(&uvar->Tsft) , XLAL_EINVAL, UVAR_STR(Tsft) " can only be used for data generation with " UVAR_STR(injectSqrtSX)", not when loading existing "  UVAR_STR(DataFiles) "\n");
     XLAL_CHECK ( uvar->IFOs == NULL , XLAL_EINVAL, UVAR_STR(IFOs) " can only be used for data generation with " UVAR_STR(injectSqrtSX) ", not when loading existing "  UVAR_STR(DataFiles) "\n");
     XLAL_CHECK ( uvar->timestampsFiles == NULL , XLAL_EINVAL,UVAR_STR(timestampsFiles) " can only be used for data generation with " UVAR_STR(injectSqrtSX) ", not when loading existing "  UVAR_STR(DataFiles) "\n");
  }
  if(uvar->injectSqrtSX !=NULL) {
     XLAL_CHECK ( uvar->timestampsFiles != NULL &&  uvar->IFOs != NULL , XLAL_EINVAL,"--injectSqrtSX requires --IFOs, --timestampsFiles \n");
  }

  return XLAL_SUCCESS;

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
  char bufGPS[32];
  fprintf (fp, "refTime  = %s;\n", XLALGPSToStr ( bufGPS, &pulsarParams->Doppler.refTime ) );
  fprintf (fp, "\n");

  /* Amplitude parameters with error-estimates */
  fprintf (fp, "aPlus    = % .6g;\n", pulsarParams->Amp.aPlus );
  fprintf (fp, "daPlus   = % .6g;\n", pulsarParams->dAmp.aPlus );
  fprintf (fp, "aCross   = % .6g;\n", pulsarParams->Amp.aCross );
  fprintf (fp, "daCross  = % .6g;\n", pulsarParams->dAmp.aCross );
  fprintf (fp, "phi0     = % .6g;\n", pulsarParams->Amp.phi0 );
  fprintf (fp, "dphi0    = % .6g;\n", pulsarParams->dAmp.phi0 );
  fprintf (fp, "psi      = % .6g;\n", pulsarParams->Amp.psi );
  fprintf (fp, "dpsi     = % .6g;\n", pulsarParams->dAmp.psi );

  /* To allow arbitrary aPlus/aCross amplitudes, we have switched to output-variables A^\nu = (aPlus, aCross, phi0, psi)
   * and hence using a different Jacobian matrix. The dh0 and dcosi recovered below are only valid when GW is emitted at
   * twice the spin frequency, and are derived from dA+, dAx. Since both the current and original calculations ignore 
   * off-diagonal items in computing the errors, the dh0 and dcosi below are expected to be slightly different from the
   * original estimate when the output-variables are B^\nu = (h0, cosi, phi0, psi). Since they are only used for sanity 
   * checks, it should not impact any usage. */
  if ( fabs(pulsarParams->Amp.aPlus) >= fabs(pulsarParams->Amp.aCross) ){ /* Assume GW at twice the spin frequency only */
    REAL8 h0, dh0, cosi, dcosi;
    h0 = pulsarParams->Amp.aPlus + sqrt( SQ(pulsarParams->Amp.aPlus) - SQ(pulsarParams->Amp.aCross) );
    if (h0 > 0)
      cosi = pulsarParams->Amp.aCross/h0;
    else
      cosi = 0;
    dh0 = (pulsarParams->Amp.aPlus + pulsarParams->dAmp.aPlus) + sqrt( SQ(pulsarParams->Amp.aPlus + pulsarParams->dAmp.aPlus) - SQ(fabs(pulsarParams->Amp.aCross) + pulsarParams->dAmp.aCross) ) - h0;
    if ( (h0+dh0) > 0 )
      dcosi = fabs( (fabs(pulsarParams->Amp.aCross) + pulsarParams->dAmp.aCross) / (h0 + dh0) - fabs(cosi) );
    else
      dcosi = 0;
    
    fprintf (fp, "h0    = % .6g;\n", h0 );
    fprintf (fp, "dh0   = % .6g;\n", dh0 );
    fprintf (fp, "cosi  = % .6g;\n", cosi );
    fprintf (fp, "dcosi = % .6g;\n", dcosi );
  }

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
      fprintf (fp, "orbitTp           = %s;\n", XLALGPSToStr ( bufGPS, &(pulsarParams->Doppler.tp) ));
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

  /* F-stat-values */
  fprintf (fp, "Fa       = % .6g  %+.6gi;\n", creal(Fcand->Fa), cimag(Fcand->Fa) );
  fprintf (fp, "Fb       = % .6g  %+.6gi;\n", creal(Fcand->Fb), cimag(Fcand->Fb) );
  fprintf (fp, "twoF     = % .6g;\n", Fcand->twoF );
  /* single-IFO F-stat values, if present */
  UINT4 X, numDet = Fcand->numDetectors;
  for ( X = 0; X < numDet ; X ++ )
    fprintf (fp, "twoF%d    = % .6g;\n", X, Fcand->twoFX[X] );
  /* BSGL */
  if ( !XLALIsREAL4FailNaN(Fcand->log10BSGL) ) /* if --computeBSGL=FALSE, the log10BSGL field was initialised to NAN - do not output it */
    fprintf (fp, "log10BSGL = % .6g;\n", Fcand->log10BSGL );

  fprintf (fp, "\nAmpFisher = " );
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

/** comparison function for our candidates toplist with alternate sorting statistic BSGL=log10BSGL */
int
compareFstatCandidates_BSGL ( const void *candA, const void *candB )
{
  REAL8 BSGL1 = ((const FstatCandidate *)candA)->log10BSGL;
  REAL8 BSGL2 = ((const FstatCandidate *)candB)->log10BSGL;
  if ( BSGL1 < BSGL2 )
    return 1;
  else if ( BSGL1 > BSGL2 )
    return -1;
  else
    return 0;

} /* compareFstatCandidates_BSGL() */

/**
 * write one 'FstatCandidate' (i.e. only Doppler-params + Fstat) into file 'fp'.
 * Return: 0 = OK, -1 = ERROR
 */
int
write_FstatCandidate_to_fp ( FILE *fp, const FstatCandidate *thisFCand, const BOOLEAN output_stats, const BOOLEAN output_orbit )
{

  if ( !fp || !thisFCand )
    return -1;

  /* write main Doppler parameters */
  fprintf (fp, "%.16g %.16g %.16g %.16g %.16g %.16g",
	   thisFCand->doppler.fkdot[0], thisFCand->doppler.Alpha, thisFCand->doppler.Delta,
	   thisFCand->doppler.fkdot[1], thisFCand->doppler.fkdot[2], thisFCand->doppler.fkdot[3] );

  if (output_stats) {
    /* add extra output-field containing per-detector FX if non-NULL */
    char extraStatsStr[256] = "";     /* defaults to empty */
    char buf0[256];
    /* BSGL */
    if ( !XLALIsREAL4FailNaN(thisFCand->log10BSGL) ) { /* if --computeBSGL=FALSE, the log10BSGL field was initialised to NAN - do not output it */
      snprintf ( extraStatsStr, sizeof(extraStatsStr), " %.9g", thisFCand->log10BSGL );
    }
    if ( thisFCand->numDetectors > 0 ) {
      for ( UINT4 X = 0; X < thisFCand->numDetectors; X ++ ) {
        snprintf ( buf0, sizeof(buf0), " %.9g", thisFCand->twoFX[X] );
        UINT4 len1 = strlen ( extraStatsStr ) + strlen ( buf0 ) + 1;
        if ( len1 > sizeof ( extraStatsStr ) ) {
          XLAL_ERROR ( XLAL_EINVAL, "assembled output string too long! (%d > %zu)\n", len1, sizeof(extraStatsStr ));
          break;      /* we can't really terminate with error in this function, but at least we avoid crashing */
        }
        strcat ( extraStatsStr, buf0 );
      } /* for X < numDet */
    } /* if FX */
    fprintf (fp, " %.9g%s", thisFCand->twoF, extraStatsStr );
  }

  if (output_orbit) {
    fprintf( fp, " %.16g %.16g %.16g %.16g %.16g",
             thisFCand->doppler.asini, thisFCand->doppler.period,
             XLALGPSGetREAL8(&thisFCand->doppler.tp),
             thisFCand->doppler.argp, thisFCand->doppler.ecc );
  }

  fprintf( fp, "\n" );

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

  XLAL_CHECK_NULL ( ( ret = LALCalloc ( 1, sizeof(*ret)) ) != NULL, XLAL_ENOMEM );
  ret->length = windowLen;

  XLAL_CHECK_NULL ( (ret->window = LALCalloc ( windowLen, sizeof( ret->window[0] ) )) != NULL, XLAL_ENOMEM );

  ret->center = &(ret->window[ windowWings ]);	/* points to central bin */

  return ret;

} /* XLALCreateScanlineWindow() */

void
XLALDestroyScanlineWindow ( scanlineWindow_t *scanlineWindow )
{
  if ( !scanlineWindow )
    return;

  if ( scanlineWindow->window )
    XLALFree ( scanlineWindow->window );

  XLALFree ( scanlineWindow );

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

  if ( rankingStatistic == RANKBY_2F ) /* F-statistic */
    {

      REAL8 twoF0 = scanWindow->center->twoF;

      for ( UINT4 i=0; i < scanWindow->length; i ++ )
        if ( scanWindow->window[i].twoF > twoF0 )
         return FALSE;

    }

  else if ( rankingStatistic == RANKBY_BSGL ) /* line-robust statistic log10BSGL */
    {

      REAL8 BSGL0 = scanWindow->center->log10BSGL;

      for ( UINT4 i=0; i < scanWindow->length; i ++ )
        if ( scanWindow->window[i].log10BSGL > BSGL0 )
         return FALSE;

    }
  else
    {
      XLALPrintError ("Unsupported ranking statistic '%d' ! Supported: 'F=0' and 'BSGL=2'.\n", rankingStatistic );
      return FALSE;
    }

  return TRUE;

} /* XLALCenterIsLocalMax() */

/**
 * Function to append one timing-info line to output file.
 *
 */
int
write_TimingInfo ( const CHAR *fname, const timingInfo_t *ti, const ConfigVariables *cfg )
{
  XLAL_CHECK ( (fname != NULL) && (ti != NULL), XLAL_EINVAL );

  CHAR *logstr_short;
  BOOLEAN logstrVerbose = FALSE;
  XLAL_CHECK ( (logstr_short = XLALGetLogString ( cfg, logstrVerbose )) != NULL, XLAL_EFUNC );

  FILE *fp;
  if ( (fp = fopen(fname,"rb" )) != NULL )
    {
      fclose(fp);
      XLAL_CHECK ( (fp = fopen( fname, "ab" ) ), XLAL_ESYS, "Failed to open existing timing-file '%s' for appending\n", fname );
    }
  else
    {
      XLAL_CHECK ( (fp = fopen( fname, "wb" ) ), XLAL_ESYS, "Failed to open new timing-file '%s' for writing\n", fname );
      fprintf ( fp, "%s", logstr_short );
      fprintf ( fp, "%2s%6s %10s %10s %10s %10s %10s %10s %10s %10s %10s %16s %12s\n",
                "%%", "NSFTs", "NFreq", "tauF[s]", "tauFEff[s]", "tauF0[s]", "FstatMethod",
                "tauMin", "tauMax", "NStart", "NTau", "tauTransFstatMap", "tauTransMarg"
              );
    }

  fprintf ( fp, "%8d %10d %10.1e %10.1e %10.1e %10s %10d %10d %10d %10d %16.1e %12.1e\n",
            ti->NSFTs, ti->NFreq, ti->tauFstat, ti->tauTemplate, ti->tauF0, XLALGetFstatInputMethodName(cfg->Fstat_in),
            ti->tauMin, ti->tauMax, ti->NStart, ti->NTau, ti->tauTransFstatMap, ti->tauTransMarg
          );

  fclose ( fp );
  XLALFree ( logstr_short );
  return XLAL_SUCCESS;

} /* write_TimingInfo() */

/* Resize histogram */
gsl_vector_int *resize_histogram(gsl_vector_int *old_hist, size_t size) {
  gsl_vector_int *new_hist = gsl_vector_int_alloc(size);
  XLAL_CHECK_NULL(new_hist != NULL, XLAL_ENOMEM);
  gsl_vector_int_set_zero(new_hist);
  if (old_hist != NULL) {
    gsl_vector_int_view old = gsl_vector_int_subvector(old_hist, 0, GSL_MIN(old_hist->size, size));
    gsl_vector_int_view new = gsl_vector_int_subvector(new_hist, 0, GSL_MIN(old_hist->size, size));
    gsl_vector_int_memcpy(&new.vector, &old.vector);
    gsl_vector_int_free(old_hist);
  }
  return new_hist;
}
