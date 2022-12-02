/*
 *  Copyright (C) 2011-2013 Karl Wette.
 *  Copyright (C) 2009-2010 Holger Pletsch.
 *
 *  Based on HierarchicalSearch.c by
 *  Copyright (C) 2005-2008 Badri Krishnan, Alicia Sintes, Bernd Machenschalk.
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
 *
 */

/*********************************************************************************/
/**
 * \defgroup lalpulsar_bin_GCT GCT Search Application
 * \ingroup lalpulsar_bin_Apps
 */

/**
 * \author Holger Pletsch
 * \file
 * \ingroup lalpulsar_bin_GCT
 * \brief Hierarchical semicoherent CW search code based on F-Statistic,
 * exploiting global-correlation coordinates (Phys.Rev.Lett. 103, 181102, 2009)
 *
 */

/* ---------- Includes -------------------- */
#include <lal/Segments.h>
#include <lal/LALString.h>
#include <lal/LineRobustStats.h>
#include <RecalcToplistStats.h>

#include "HierarchSearchGCT.h"

#ifdef GC_SSE2_OPT
#include <gc_hotloop_sse2.h>
#else
#define ALRealloc LALRealloc
#define ALFree LALFree
#endif

/* ---------- Defines -------------------- */
/* #define DIAGNOSISMODE 1 */
#define NUDGE	10*LAL_REAL8_EPS

#define TRUE (1==1)
#define FALSE (1==0)

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define GET_GCT_CHECKPOINT read_gct_checkpoint // (cptname, semiCohToplist, NULL, &count)
#define SET_GCT_CHECKPOINT write_gct_checkpoint
char**global_argv;
int global_argc;

#define FSTART          100.0	/**< Default Start search frequency */
#define FBAND           0.0  /**< Default search band */
#define FDOT            0.0       /**< Default value of first spindown */
#define DFDOT           0.0       /**< Default range of first spindown parameter */
#define F2DOT           0.0       /**< Default value of second spindown */
#define F3DOT           0.0       /**< Default value of third spindown */
#define DF2DOT          0.0       /**< Default range of second spindown parameter */
#define DF3DOT          0.0       /**< Default range of third spindown parameter */
#define SKYREGION       "allsky" /**< default sky region to search over -- just a single point*/
#define MISMATCH        0.3       /**< Default for metric grid maximal mismatch value */
#define DALPHA          0.001   /**< Default resolution for isotropic or flat grids */
#define DDELTA          0.001   /**< Default resolution for isotropic or flat grids */
#define FSTATTHRESHOLD  2.6	/**< Default threshold on Fstatistic for peak selection */
#define NCAND1          10      /**< Default number of candidates to be followed up from first stage */
#define FNAMEOUT        "./HS_GCT.out"  /**< Default output file basename */
#ifndef LAL_INT4_MAX
#define LAL_INT4_MAX    2147483647
#endif

#define BLOCKSIZE_REALLOC 50

#define Vorb_GCT   = 2.9785e04;
#define Vspin_GCT  = 465.10;
#define REARTH_GCT = 6.378140e06;
#define C_GCT      = 299792458;

/* ---------- Macros -------------------- */
#define HSMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define HSMIN(x,y) ( (x) < (y) ? (x) : (y) )

#define GETTIME() (uvar_outputTiming ? XLALGetCPUTime() : 0)
//#define GETTIME() (uvar_outputTiming ? XLALGetTimeOfDay() : 0)

/* ---------- Exported types ---------- */
/** useful variables for each hierarchical stage */
typedef struct {
  CHAR  *sftbasename;    /**< filename pattern for sfts */
  LIGOTimeGPS tStartGPS; /**< start and end time of stack */
  REAL8 tObs;            /**< tEndGPS - tStartGPS */
  REAL8 refTime;         /**< reference time for pulsar params */
  PulsarSpinRange spinRange_startTime; /**< freq and fdot range at start-time of observation */
  PulsarSpinRange spinRange_endTime;   /**< freq and fdot range at end-time of observation */
  PulsarSpinRange spinRange_refTime;   /**< freq and fdot range at the reference time */
  PulsarSpinRange spinRange_midTime;   /**< freq and fdot range at mid-time of observation */
  EphemerisData *edat;             /**< ephemeris data for XLALBarycenter */
  LIGOTimeGPSVector *midTstack;    /**< timestamps vector for mid time of each stack */
  LIGOTimeGPSVector *startTstack;  /**< timestamps vector for start time of each stack */
  LIGOTimeGPSVector *endTstack;    /**< timestamps vector for end time of each stack */
  LIGOTimeGPS minStartTimeGPS;     /**< all sft data must be after this time */
  LIGOTimeGPS maxStartTimeGPS;       /**< all sft timestamps must be before this GPS time */
  UINT4 blocksRngMed;              /**< blocksize for running median noise floor estimation */
  UINT4 Dterms;                    /**< size of Dirichlet kernel for Fstat calculation */
  UINT4 DtermsRecalc;              /**< Recalc: size of Dirichlet kernel for Fstat calculation */
  LALStringVector* assumeSqrtSX;   /**< Assume stationary Gaussian noise with detector noise-floors sqrt{SX}" */
  /* parameters describing the coherent data-segments */
  REAL8 tStack;                    /**< duration of stacks */
  UINT4 nStacks;                   /**< number of stacks */
  LALSegList *segmentList;         /**< parsed segment list read from user-specified input file --segmentList */
  BSGLSetup *BSGLsetup;           /**< pre-computed setup for line-robust statistic BSGL */
  REAL8 dFreqStack;                /**< frequency resolution of Fstat calculation */
  REAL8 df1dot;                    /**< coarse grid resolution in spindown */
  REAL8 df2dot;                    /**< coarse grid resolution in 2nd spindown */
  REAL8 df3dot;                    /**< coarse grid resolution in 3rd spindown */
  UINT4 extraBinsFstat;            /**< Extra Fstat frequency bins required to cover residual spindowns */
  UINT4 binsFstatSearch;	   /**< nominal number of Fstat frequency bins in search band */
  UINT4 nf1dot;			/**< number of 1st spindown Fstat bins */
  UINT4 nf2dot;			/**< number of 2nd spindown Fstat bins */
  UINT4 nf3dot;			/**< number of 3rd spindown Fstat bins */
  int SSBprec;                     /**< SSB transform precision */
  FstatMethodType Fmethod;         //!< which Fstat-method/algorithm to use
  BOOLEAN recalcToplistStats;	   //!< do additional analysis for all toplist candidates, output F, FXvector for postprocessing */
  FstatMethodType FmethodRecalc;   //!< which Fstat-method/algorithm to use for the recalc step
  REAL8 mismatch1;                 /**< 'mismatch1' user-input needed here internally ... */
  UINT4 nSFTs;                     /**< total number of SFTs */
  LALStringVector *detectorIDs;    /**< vector of detector IDs */
  REAL4 NSegmentsInvX[PULSAR_MAX_DETECTORS]; /**< effective inverse number of segments per detector (needed for correct averaging in single-IFO F calculation) */
  FstatInputVector* Fstat_in_vec;	/**< Original wide-parameter search: vector of Fstat input data structures for XLALComputeFstat(), one per stack */
  FstatInputVector* Fstat_in_vec_recalc; /**< Recalculate the toplist: Vector of Fstat input data structures for XLALComputeFstat(), one per stack */
  PulsarParamsVector *injectionSources; ///< Source parameters to inject: comma-separated list of file-patterns and/or direct config-strings ('{...}')
  BOOLEAN collectFstatTiming;		///< flag whether to collect and output F-stat timing info
} UsefulStageVariables;


/**
 * Struct holding various timing measurements and relevant search parameters.
 * This is used to fit timing-models with measured times to predict search run-times
 */
typedef struct
{
  UINT4 Nseg;			///< number of semi-coherent segments
  UINT4 Ndet;			///< number of detectors
  UINT4 Tcoh;			///< length of coherent segments in seconds
  UINT4 Nsft;			///< total number of SFTs
  UINT4 Ncand;			///< length of toplists

  UINT4 NFreqCo;		///< total number of frequency bins computed in coarse grid (including sidebands!)
  REAL8 Ncoh;			///< number of coarse-grid Fstat templates ('coherent')
  REAL8 Ninc;			///< number of fine-grid templates ('incoherent')

  const char* FstatMethodStr;	///< Fstat-method used
  const char* RecalcMethodStr;	///< Fstat-method used

  // ----------
  // extended timing model:
  // runtime = Nseg * Ndet * Ncoh * tauF + Nseg * Ninc * tauSumFine + Ninc * tauExtraStats + Ntop * tauRecalc + time_Other
  REAL8 tau_Fstat;		//< time to compute F-stat for one coarse-grid template, one detector, one segment
  REAL8 tau_SumF;		//< time to sum F-stat values of one segment for one fine-grid point
  REAL8 tau_Bayes;		//< time to compute final Bayes-factor statistics for one fine-grid point
  REAL8 tau_Recalc;		//< time to recalc one toplist-entry
  REAL8 time_Other;		//< all non-scaling leftovers in overall timeing
  // ----------

} timingInfo_t;


/* ------------------------ Functions -------------------------------- */
void SetUpSFTs( LALStatus *status, UsefulStageVariables *in );
void PrintFstatVec( LALStatus *status, FstatResults *in, FILE *fp, PulsarDopplerParams *thisPoint,
                    LIGOTimeGPS refTime, INT4 stackIndex);
void PrintCatalogInfo( LALStatus *status, const SFTCatalog *catalog, FILE *fp );
void PrintStackInfo( LALStatus *status, const SFTCatalogSequence *catalogSeq, FILE *fp );
void UpdateSemiCohToplists ( LALStatus *status, toplist_t *list1, toplist_t *list2, toplist_t *list3, FineGrid *in, REAL8 f1dot_fg, REAL8 f2dot_fg, REAL8 f3dot_fg, UsefulStageVariables *usefulparams, REAL4 NSegmentsInv, REAL4 *NSegmentsInvX, BOOLEAN have_f3dot );

void UpdateSemiCohToplistsOptimTriple ( LALStatus *status,
                             SortBy_t toplist_sortby,
                             toplist_t *list1,
                             toplist_t *list2,
                             toplist_t *list3,
                             FineGrid *in,
                             REAL8 f1dot_fg,
                             REAL8 f2dot_fg,
                             REAL8 f3dot_fg,
                             UsefulStageVariables *usefulparams,
                             REAL4 NSegmentsInv,
                             REAL4 *NSegmentsInvX,
                             BOOLEAN have_f3dot
                             );

void GetSegsPosVelAccEarthOrb( LALStatus *status, REAL8VectorSequence **posSeg,
                               REAL8VectorSequence **velSeg, REAL8VectorSequence **accSeg,
                               UsefulStageVariables *usefulparams );
static inline INT4 ComputeU1idx( REAL8 freq_event, REAL8 f1dot_event, REAL8 A1, REAL8 B1, REAL8 U1start, REAL8 U1winInv );
void ComputeU2idx( REAL8 freq_event, REAL8 f1dot_event, REAL8 A2, REAL8 B2, REAL8 U2start, REAL8 U2winInv,
                   INT4 *U2idx);
int compareCoarseGridUindex( const void *a, const void *b );
int compareFineGridNC( const void *a,const void *b );
int compareFineGridsumTwoF( const void *a,const void *b );

SFTCatalogSequence *XLALSetUpStacksFromSegmentList ( const SFTCatalog *catalog, const LALSegList *segList );

int XLALExtrapolateToplistPulsarSpins ( toplist_t *list,
					const LIGOTimeGPS usefulParamsRefTime,
					const LIGOTimeGPS finegridRefTime);

static int write_TimingInfo ( const CHAR *fname, const timingInfo_t *ti );
static inline REAL4 findLoudestTwoF ( const FstatResults *in );

/* ---------- Global variables -------------------- */
LALStatus *global_status; /* a global pointer to main()s head of the LALStatus structure */
char *global_column_headings_stringp;

// XLALReadSegmentsFromFile(): applications which still must support
// the deprecated 4-column format should set this variable to non-zero
extern int XLALReadSegmentsFromFile_support_4column_format;

/* ###################################  MAIN  ################################### */

int main( int argc, char *argv[]) {
  LALStatus XLAL_INIT_DECL(status);

  /* temp loop variables: generally k loops over segments and j over SFTs in a stack */
  UINT4 k;
  UINT4 skyGridCounter; /* coarse sky position counter */
  UINT4 f1dotGridCounter; /* coarse f1dot position counter */

  /* GPS timestamp vectors */
  LIGOTimeGPSVector *midTstack = NULL;
  LIGOTimeGPSVector *startTstack = NULL;
  LIGOTimeGPSVector *endTstack = NULL;

  /* General GPS times */
  LIGOTimeGPS XLAL_INIT_DECL(refTimeGPS);
  LIGOTimeGPS XLAL_INIT_DECL(tMidGPS);

  /* GPS time used for each segment's midpoint */
  LIGOTimeGPS XLAL_INIT_DECL(midTstackGPS);
  REAL8 timeDiffSeg; /* Difference to tMidGPS (midpoint) of Tobs */

  /* pos, vel, acc at midpoint of segments */
  REAL8VectorSequence *posStack = NULL;
  REAL8VectorSequence *velStack = NULL;
  REAL8VectorSequence *accStack = NULL;

  /* duration of each segment */
  REAL8 tStack;

  /* number of segments */
  UINT4 nStacks;

  /* Total observation time */
  REAL8 tObs;

  /* SFT related stuff */
  static LIGOTimeGPS minStartTimeGPS, maxStartTimeGPS;

  /* some useful variables for each stage */
  UsefulStageVariables XLAL_INIT_DECL(usefulParams);

  /* F-statistic computation related stuff */
  FstatResults* Fstat_res = NULL;			// Pointer to Fstat results structure, will be allocated by XLALComputeFstat()
  FstatQuantities Fstat_what = FSTATQ_2F;		// Quantities to be computed by XLALComputeFstat()
  UINT4 binsFstat1;

  /* Semicoherent variables */
  static SemiCoherentParams semiCohPar;

  /* coarse grid */
  CoarseGrid XLAL_INIT_DECL(coarsegrid);
  REAL8 dFreqStack; /* frequency resolution of Fstat calculation */
  REAL8 df1dot;  /* coarse grid resolution in spindown */
  UINT4 ifdot;  /* counter for coarse-grid spindown values */
  REAL8 df2dot;  /* coarse grid resolution in 2nd spindown */
  UINT4 if2dot;  /* counter for coarse-grid 2nd spindown values */
  REAL8 df3dot;  /* coarse grid resolution in 3rd spindown */
  UINT4 if3dot;  /* counter for coarse-grid 3rd spindown values */

  /* fine grid */
  FineGrid finegrid;
  UINT4 nf1dots_fg = 1; /* number of frequency and spindown values */
  REAL8 gammaRefine, sigmasq;  /* refinement factor and variance */
  UINT4 nf2dots_fg=1;          /* number of second spindown values */
  REAL8 gamma2Refine, sigma4;  /* 2nd spindown refinement factor and 4th moment */
  UINT4 nf3dots_fg=1;          /* number of third spindown values */
  REAL8 gamma3Refine=1;  /* 3rd spindown refinement */
  
  /* GCT helper variables */
  UINT4 if1dot_fg, if2dot_fg, if3dot_fg;
  UINT4 ifreq;
  INT4  U1idx;
  REAL8 myf0, freq_event, f1dot_event;
  REAL8 dfreq_fg, df1dot_fg, freqmin_fg, f1dotmin_fg, freqband_fg;
  REAL8 df2dot_fg, f2dotmin_fg;
  REAL8 df3dot_fg, f3dotmin_fg;
  REAL8 u1start, u1win, u1winInv;
  REAL8 freq_fg, f1dot_fg, f2dot_fg, f3dot_fg, f1dot_event_fg;
  REAL8 A1, B1;
  // currently unused: REAL8 A2;
  REAL8 B2; /* GCT helper variables for faster calculation of u1 or u2 */
  REAL8 pos[3];
  REAL8 vel[3];
  // currently unused: REAL8 acc[3];
  REAL8 cosAlpha, sinAlpha, cosDelta, sinDelta;
  REAL8 nvec[3]; /* unit vector pointing to sky position */

  /* These vars are currently not used, but eventually in the future.
     INT4  U2idx, NumU2idx;
     REAL8 myf0max, u2start, u2end;
     REAL8 u2win, u2winInv;
  */

  /* fstat candidate structure for candidate toplist*/
  toplist_t *semiCohToplist=NULL;
  toplist_t *semiCohToplist2=NULL;	// only used for SORTBY_DUAL_F_BSGL, SORTBY_TRIPLE_BStSGLtL or SORTBY_F_BSGLtL_BtSGLtL
  toplist_t *semiCohToplist3=NULL;	// only used for SORTBY_TRIPLE_BStSGLtL and SORTBY_F_BSGLtL_BtSGLtL

  /* template and grid variables */
  static DopplerSkyScanInit scanInit;   /* init-structure for DopperScanner */
  DopplerSkyScanState XLAL_INIT_DECL(thisScan); /* current state of the Doppler-scan */
  static PulsarDopplerParams dopplerpos;               /* current search-parameters */
  static PulsarDopplerParams thisPoint;
  UINT4 oldcg=0, oldfg=0;

  /* temporary storage for spinrange vector */
  static PulsarSpinRange spinRange_Temp;

  /* variables for logging */
  CHAR *fnamelog=NULL;
  FILE *fpLog=NULL;
  CHAR *logstr=NULL;

  /* output candidate files and file pointers */
  CHAR *fnameSemiCohCand=NULL;
  CHAR *fnameFstatVec1=NULL;
  FILE *fpFstat1=NULL;

  /* checkpoint filename */
  CHAR *uvar_fnameChkPoint = NULL;

  /* user variables */
  BOOLEAN uvar_log = FALSE;     /* logging done if true */

  BOOLEAN uvar_printCand1 = FALSE;      /* if 1st stage candidates are to be printed */
  BOOLEAN uvar_printFstat1 = FALSE;
  BOOLEAN uvar_loudestTwoFPerSeg = FALSE;	// output loudest per-segment Fstat candidates
  BOOLEAN uvar_semiCohToplist = TRUE; /* if overall first stage candidates are to be output */

  LALStringVector* uvar_assumeSqrtSX = NULL;    /* Assume stationary Gaussian noise with detector noise-floors sqrt{SX}" */

  BOOLEAN uvar_recalcToplistStats = FALSE; 	/* Do additional analysis for all toplist candidates, output F, FXvector for postprocessing */
  BOOLEAN uvar_loudestSegOutput = FALSE; 	/* output extra info about loudest segment; requires recalcToplistStats */

  // ----- Line robust stats parameters ----------
  BOOLEAN uvar_computeBSGL = FALSE;          	/* In Fstat loop, compute line-robust statistic (BSGL=log10BSGL) using single-IFO F-stats */
  BOOLEAN uvar_BSGLlogcorr = FALSE;		/* compute log-correction in line-robust statistic BSGL (slower) or not (faster) */
  REAL8   uvar_Fstar0sc = 0.0;			/* (semi-coherent) BSGL transition-scale parameter 'Fstar0sc=Nseg*Fstar0coh', see documentation for XLALCreateBSGLSetup() for details */
  LALStringVector *uvar_oLGX = NULL;       	/* prior per-detector line-vs-Gauss odds ratios 'oLGX', see XLALCreateBSGLSetup() for details */
  BOOLEAN uvar_getMaxFperSeg = FALSE;          	/* In Fstat loop, compute maximum F and FX over segments */
  // --------------------------------------------

  REAL8 uvar_dAlpha = DALPHA;   /* resolution for flat or isotropic grids -- coarse grid*/
  REAL8 uvar_dDelta = DDELTA;
  REAL8 uvar_f1dot = FDOT;      /* first spindown value */
  REAL8 uvar_f1dotBand = DFDOT; /* range of first spindown parameter */
  REAL8 uvar_f2dot = F2DOT;     /* second spindown value */
  REAL8 uvar_f2dotBand = DF2DOT; /* range of second spindown parameter */
  REAL8 uvar_f3dot = F3DOT;     /* second spindown value */
  REAL8 uvar_f3dotBand = DF3DOT; /* range of second spindown parameter */
  REAL8 uvar_Freq = FSTART;
  REAL8 uvar_FreqBand = FBAND;

  REAL8 uvar_dFreq = 0;
  REAL8 uvar_df1dot = 0; /* coarse grid frequency and spindown resolution */
  REAL8 uvar_df2dot = 0; /* coarse grid second spindown resolution */
  REAL8 uvar_df3dot = 0; /* coarse grid third spindown resolution */

  REAL8 uvar_ThrF = FSTATTHRESHOLD; /* threshold of Fstat to select peaks */
  REAL8 uvar_mismatch1 = MISMATCH; /* metric mismatch for first stage coarse grid */

  REAL8 uvar_minStartTime1 = 0;
  REAL8 uvar_maxStartTime1 = LAL_INT4_MAX;

  REAL8 uvar_refTime = 0;
  INT4 uvar_nCand1 = NCAND1; /* number of candidates to be followed up from first stage */

  INT4 uvar_blocksRngMed = FstatOptionalArgsDefaults.runningMedianWindow;

  REAL8 uvar_tStack = 0;
  INT4  uvar_nStacksMax = 1;
  CHAR *uvar_segmentList = NULL;	/**< ALTERNATIVE: file containing a pre-computed segment list of tuples (startGPS endGPS duration[h] NumSFTs) */

  INT4 uvar_Dterms = FstatOptionalArgsDefaults.Dterms;
  INT4 uvar_DtermsRecalc = FstatOptionalArgsDefaults.Dterms;
  INT4 uvar_SSBprecision = FstatOptionalArgsDefaults.SSBprec;
  INT4 uvar_gammaRefine = 1;
  INT4 uvar_gamma2Refine = 1;
  INT4 uvar_metricType1 = LAL_PMETRIC_COH_PTOLE_ANALYTIC;
  INT4 uvar_gridType1 = GRID_METRIC;
  INT4 uvar_skyPointIndex = -1;

  CHAR *uvar_ephemEarth;	/**< Earth ephemeris file to use */
  CHAR *uvar_ephemSun;		/**< Sun ephemeris file to use */

  CHAR *uvar_skyRegion = NULL;
  CHAR *uvar_fnameout = NULL;
  CHAR *uvar_DataFiles1 = NULL;
  CHAR *uvar_skyGridFile=NULL;
  INT4 uvar_numSkyPartitions = 0;
  INT4 uvar_partitionIndex = 0;
  INT4 uvar_SortToplist = 0;

  CHAR *uvar_outputTiming = NULL;
  CHAR *uvar_outputTimingDetails = NULL;

  int uvar_FstatMethod = FstatOptionalArgsDefaults.FstatMethod;
  int uvar_FstatMethodRecalc = FstatOptionalArgsDefaults.FstatMethod;

  timingInfo_t XLAL_INIT_DECL(timing);

  LALStringVector *uvar_injectionSources = NULL;

  // timing values
  REAL8 tic_RecalcToplist, time_RecalcToplist = 0;
  REAL8 tic_Fstat, time_Fstat = 0;
  REAL8 tic_SumFine, time_SumFine = 0;
  REAL8 tic_ExtraStats, time_ExtraStats = 0;
  REAL8 tic_Start, time_Total = 0;

  global_status = &status;

  global_argv = argv;
  global_argc = argc;

  // XLALReadSegmentsFromFile(): continue to support deprecated 4-column format (startGPS endGPS duration NumSFTs, duration is ignored)
  XLALReadSegmentsFromFile_support_4column_format = 1;

  uvar_ephemEarth = XLALStringDuplicate("earth00-40-DE405.dat.gz");
  uvar_ephemSun = XLALStringDuplicate("sun00-40-DE405.dat.gz");

  uvar_skyRegion = LALCalloc( strlen(SKYREGION) + 1, sizeof(CHAR) );
  strcpy(uvar_skyRegion, SKYREGION);

  uvar_fnameout = LALCalloc( strlen(FNAMEOUT) + 1, sizeof(CHAR) );
  strcpy(uvar_fnameout, FNAMEOUT);

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register user input variables */
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_log,                 "log",                 BOOLEAN,      0,   OPTIONAL,   "Write log file") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_semiCohToplist,      "semiCohToplist",      BOOLEAN,      0,   OPTIONAL,   "Print toplist of semicoherent candidates" ) == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_DataFiles1,          "DataFiles1",          STRING,       0,   REQUIRED,   "1st SFT file pattern. Possibilities are:\n"
                                          " - '<SFT file>;<SFT file>;...', where <SFT file> may contain wildcards\n - 'list:<file containing list of SFT files>'") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_skyRegion,           "skyRegion",           STRING,       0,   OPTIONAL,   "sky-region polygon (or 'allsky')") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_numSkyPartitions,    "numSkyPartitions",    INT4,         0,   OPTIONAL,   "No. of (equi-)partitions to split skygrid into") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_partitionIndex,      "partitionIndex",      INT4,         0,   OPTIONAL,   "Index [0,numSkyPartitions-1] of sky-partition to generate") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_skyGridFile,         "skyGridFile",         STRING,       0,   OPTIONAL,   "sky-grid file") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_dAlpha,              "dAlpha",              REAL8,        0,   OPTIONAL,   "Resolution for flat or isotropic coarse grid") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_dDelta,              "dDelta",              REAL8,        0,   OPTIONAL,   "Resolution for flat or isotropic coarse grid") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_Freq,                "Freq",                REAL8,        'f', OPTIONAL,   "Start search frequency") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_dFreq,               "dFreq",               REAL8,        0,   OPTIONAL,   "Frequency resolution (required if nonzero FreqBand)") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_FreqBand,            "FreqBand",            REAL8,        'b', OPTIONAL,   "Search frequency band") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_f1dot,               "f1dot",               REAL8,        0,   OPTIONAL,   "Spindown parameter") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_df1dot,              "df1dot",              REAL8,        0,   OPTIONAL,   "Spindown resolution (required if nonzero f1dotBand)") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_f1dotBand,           "f1dotBand",           REAL8,        0,   OPTIONAL,   "Spindown Range") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_f2dot,               "f2dot",               REAL8,        0,   OPTIONAL,   "2nd spindown parameter") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_df2dot,              "df2dot",              REAL8,        0,   OPTIONAL,   "2nd spindown resolution (required if nonzero f2dotBand)") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_f2dotBand,           "f2dotBand",           REAL8,        0,   OPTIONAL,   "2nd spindown Range") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_f3dot,               "f3dot",               REAL8,        0,   OPTIONAL,   "3rd spindown parameter") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_df3dot,              "df3dot",              REAL8,        0,   OPTIONAL,   "3rd spindown resolution (required if nonzero f3dotBand)") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_f3dotBand,           "f3dotBand",           REAL8,        0,   OPTIONAL,   "3rd spindown Range") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_ThrF,                "peakThrF",            REAL8,        0,   OPTIONAL,   "Fstat Threshold") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_mismatch1,           "mismatch1",           REAL8,        'm', OPTIONAL,   "1st stage mismatch") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_gridType1,           "gridType1",           INT4,         0,   OPTIONAL,   "0=flat, 1=isotropic, 2=metric, 3=file") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_metricType1,         "metricType1",         INT4,         0,   OPTIONAL,   "0=none, 1=Ptole-analytic, 2=Ptole-numeric, 3=exact") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_gammaRefine,         "gammaRefine",         INT4,         'g', OPTIONAL,   "Refinement of fine grid (default: use segment times)") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_gamma2Refine,        "gamma2Refine",        INT4,         'G', OPTIONAL,   "Refinement of f2dot fine grid (default: use segment times, -1=use gammaRefine)") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_fnameout,            "fnameout",            STRING,       'o', OPTIONAL,   "Output filename") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_fnameChkPoint,       "fnameChkPoint",       STRING,       0,   OPTIONAL,   "Checkpoint filename") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_nCand1,              "nCand1",              INT4,         'n', OPTIONAL,   "No. of candidates to output") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_printCand1,          "printCand1",          BOOLEAN,      0,   OPTIONAL,   "Print 1st stage candidates") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_refTime,             "refTime",             REAL8,        0,   OPTIONAL,   "Ref. time for pulsar pars [Default: mid-time]") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_ephemEarth,          "ephemEarth",          STRING,       0,   OPTIONAL,   "Location of Earth ephemeris file") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_ephemSun,            "ephemSun",            STRING,       0,   OPTIONAL,   "Location of Sun ephemeris file") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_minStartTime1,       "minStartTime1",       REAL8,        0,   OPTIONAL,   "1st stage: Only use SFTs with timestamps starting from (including) this GPS time") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_maxStartTime1,       "maxStartTime1",       REAL8,        0,   OPTIONAL,   "1st stage: Only use SFTs with timestamps up to (excluding) this GPS time") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_printFstat1,         "printFstat1",         BOOLEAN,      0,   OPTIONAL,   "Print 1st stage Fstat vectors") == XLAL_SUCCESS, XLAL_EFUNC);

  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_assumeSqrtSX,        "assumeSqrtSX",        STRINGVector, 0,   OPTIONAL,   "Don't estimate noise-floors but assume (stationary) per-IFO sqrt{SX} (if single value: use for all IFOs)") == XLAL_SUCCESS, XLAL_EFUNC);

  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_nStacksMax,          "nStacksMax",          INT4,         0,   OPTIONAL,   "Maximum No. of segments" ) == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_tStack,              "tStack",              REAL8,        'T', OPTIONAL,   "Duration of segments (sec)" ) == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_segmentList,         "segmentList",         STRING,       0,   OPTIONAL,   "ALTERNATIVE: file containing a segment list: lines of form <startGPS endGPS duration[h] NumSFTs>") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_recalcToplistStats,  "recalcToplistStats",  BOOLEAN,      0,   OPTIONAL,   "Additional analysis for toplist candidates, recalculate 2F, 2FX at finegrid") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_loudestSegOutput,    "loudestSegOutput",    BOOLEAN,      0,   OPTIONAL,   "Output extra info about loudest segment; (requires --recalcToplistStats)") == XLAL_SUCCESS, XLAL_EFUNC);

  // ----- Line robust stats parameters ----------
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_computeBSGL,         "computeBSGL",         BOOLEAN,      0,   OPTIONAL,   "Compute and output line-robust statistic (BSGL)") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_Fstar0sc,            "Fstar0sc",            REAL8,        0,   OPTIONAL,   "BSGL: semi-coh transition-scale parameter 'Fstar0sc=Nseg*Fstar0coh'" ) == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_oLGX,                "oLGX",                STRINGVector, 0,   OPTIONAL,   "BSGL: prior per-detector line-vs-Gauss odds 'oLGX' (Defaults to oLGX=1/Ndet)") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_BSGLlogcorr,         "BSGLlogcorr",         BOOLEAN,      0,   DEVELOPER,  "BSGL: include log-correction terms (slower) or not (faster)") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_getMaxFperSeg,       "getMaxFperSeg",       BOOLEAN,      0,   OPTIONAL,   "Compute and output maximum F and FX over segments") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_SortToplist,         "SortToplist",         INT4,         0,   OPTIONAL,   "Sort toplist by: 0=Fstat, 1=nc, 2=B_S/GL, 3='Fstat + B_S/GL', 4=B_S/GLtL, 5=B_tS/GLtL, 6='B_S/GL + B_S/GLtL + B_tS/GLtL', 7='Fstat + B_S/GLtL + B_tS/GLtL' ") == XLAL_SUCCESS, XLAL_EFUNC);
  // --------------------------------------------

  XLAL_CHECK_MAIN( XLALRegisterNamedUvarAuxData( &uvar_FstatMethod, "FstatMethod", UserEnum, XLALFstatMethodChoices(), 0, OPTIONAL, "F-statistic method to use" ) == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvarAuxData( &uvar_FstatMethodRecalc, "FstatMethodRecalc", UserEnum, XLALFstatMethodChoices(), 0, OPTIONAL, "F-statistic method to use for recalc" ) == XLAL_SUCCESS, XLAL_EFUNC);

  /* developer user variables */
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_blocksRngMed,        "blocksRngMed",        INT4,         0,   DEVELOPER,  "RngMed block size") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvarAuxData( &uvar_SSBprecision, "SSBprecision", UserEnum, &SSBprecisionChoices, 0, DEVELOPER, "Precision for SSB transform") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_Dterms,              "Dterms",              INT4,         0,   DEVELOPER,  "Number of kernel terms (single-sided) to use in\na) Dirichlet kernel if FstatMethod=Demod*\nb) sinc-interpolation kernel if FstatMethod=Resamp*" ) == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_DtermsRecalc,        "DtermsRecalc",        INT4,         0,   DEVELOPER,  "Same as 'Dterms', applies to 'Recalc' step" ) == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_skyPointIndex,       "skyPointIndex",       INT4,         0,   DEVELOPER,  "Only analyze this skypoint in grid" ) == XLAL_SUCCESS, XLAL_EFUNC);

  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_outputTiming,        "outputTiming",        STRING,       0,   DEVELOPER,  "Append timing information into this file") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_outputTimingDetails, "outputTimingDetails", STRING,       0,   DEVELOPER,  "Append detailed averaged F-stat timing information to this file") == XLAL_SUCCESS, XLAL_EFUNC);

  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_loudestTwoFPerSeg,   "loudestTwoFPerSeg",   BOOLEAN,      0, DEVELOPER, "Output loudest per-segment Fstat values into file '_loudestTwoFPerSeg'" ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* inject signals into the data being analyzed */
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar ( &uvar_injectionSources, "injectionSources",      STRINGVector, 0, DEVELOPER, "%s", InjectionSourcesHelpString) == XLAL_SUCCESS, XLAL_EFUNC );

  /* read all command line variables */
  BOOLEAN should_exit = 0;
  XLAL_CHECK_MAIN( XLALUserVarReadAllInput(&should_exit, argc, argv, lalPulsarVCSInfoList) == XLAL_SUCCESS, XLAL_EFUNC);
  if (should_exit)
    return(1);

  /* assemble version string */
  CHAR *VCSInfoString;
  if ( (VCSInfoString = XLALVCSInfoString(lalPulsarVCSInfoList, 0, "%% ")) == NULL ) {
    XLALPrintError("XLALVCSInfoString failed.\n");
    return HIERARCHICALSEARCH_ESUB;
  }

  LogPrintfVerbatim( LOG_DEBUG, "Code-version: %s\n", VCSInfoString );
  // LogPrintfVerbatim( LOG_DEBUG, "CFS Hotloop variant: %s\n", OptimisedHotloopSource );

  /* some basic sanity checks on user vars */
  if ( uvar_nStacksMax < 1) {
    fprintf(stderr, "Invalid number of segments!\n");
    return( HIERARCHICALSEARCH_EBAD );
  }

#ifndef EXP_NO_NUM_COUNT
  /* check that the numbercount can't exceed the data type */
  {
    UINT8 maxseg = 1;
    maxseg = maxseg << (8*sizeof(FINEGRID_NC_T));
    maxseg -= 1;
    if ( (UINT8)uvar_nStacksMax > maxseg) {
      fprintf(stderr,
              "Number of segments exceeds %" LAL_UINT8_FORMAT "!\n"
              "Compile without GC_SSE2_OPT to extend the available segment range\n",
              maxseg);
      return( HIERARCHICALSEARCH_EBAD );
    }
  }
#endif

  if ( uvar_blocksRngMed < 1 ) {
    fprintf(stderr, "Invalid Running Median block size\n");
    return( HIERARCHICALSEARCH_EBAD );
  }

  if ( uvar_ThrF < 0 ) {
    fprintf(stderr, "Invalid value of F-statistic threshold\n");
    return( HIERARCHICALSEARCH_EBAD );
  }

  if ( uvar_f3dotBand != 0 && ( !XLALUserVarWasSet(&uvar_gammaRefine) || !XLALUserVarWasSet(&uvar_gammaRefine) || uvar_gammaRefine != 1 || uvar_gamma2Refine != 1 )){
	fprintf(stderr, "Search over 3rd spindown is available only with gammaRefine AND gamma2Refine manually set to 1!\n");
	return( HIERARCHICALSEARCH_EVAL );
  }

  /* 2F threshold for semicoherent stage */
#ifndef EXP_NO_NUM_COUNT
  REAL4 TwoFthreshold = 2.0 * uvar_ThrF;
#endif

  if ( (uvar_SortToplist < 0) || (uvar_SortToplist >= SORTBY_LAST) ) {
    XLALPrintError ( "Invalid value %d specified for toplist sorting, must be within [0, %d]\n", uvar_SortToplist, SORTBY_LAST - 1 );
    return( HIERARCHICALSEARCH_EBAD );
  }
  if ( (uvar_SortToplist == SORTBY_BSGL || uvar_SortToplist == SORTBY_DUAL_F_BSGL ||
        uvar_SortToplist == SORTBY_BSGLtL || uvar_SortToplist == SORTBY_BtSGLtL ||
        uvar_SortToplist == SORTBY_TRIPLE_BStSGLtL ||
        uvar_SortToplist == SORTBY_F_BSGLtL_BtSGLtL) && !uvar_computeBSGL ) {
    fprintf(stderr, "Toplist sorting by BSGL[tL] only possible if --computeBSGL given.\n");
    return( HIERARCHICALSEARCH_EBAD );
  }
  if ( ( uvar_SortToplist == SORTBY_BSGLtL || uvar_SortToplist == SORTBY_BtSGLtL ||
        uvar_SortToplist == SORTBY_TRIPLE_BStSGLtL ||
        uvar_SortToplist == SORTBY_F_BSGLtL_BtSGLtL ) && !uvar_getMaxFperSeg ) {
    fprintf(stderr, "Toplist sorting by B[t]SGLtL only possible if --getMaxFperSeg given.\n");
    return( HIERARCHICALSEARCH_EBAD );
  }

  /* create toplist -- semiCohToplist has the same structure
     as a fstat candidate, so treat it as a fstat candidate */
  if ( uvar_SortToplist == SORTBY_DUAL_F_BSGL )	// special treatement of 'dual' toplists: 1st one sorted by 'F', 2nd one by 'BSGL'
    {
      XLAL_CHECK ( 0 == create_gctFstat_toplist ( &semiCohToplist, uvar_nCand1, SORTBY_F ),
                   XLAL_EFUNC, "create_gctFstat_toplist() failed for nCand=%d and sortBy=%d\n", uvar_nCand1, SORTBY_F );
      XLAL_CHECK ( 0 == create_gctFstat_toplist ( &semiCohToplist2, uvar_nCand1, SORTBY_BSGL ),
                   XLAL_EFUNC, "create_gctFstat_toplist() failed for nCand=%d and sortBy=%d\n", uvar_nCand1, SORTBY_BSGL );
    }
  else if ( uvar_SortToplist == SORTBY_TRIPLE_BStSGLtL )// special treatement of 'triple' toplists: 1st one sorted by 'B_S/GL', 2nd one by 'B_S/GLtL', 3rd by 'B_tS/GLtL'
    {
      XLAL_CHECK ( 0 == create_gctFstat_toplist ( &semiCohToplist, uvar_nCand1, SORTBY_BSGL ),
                   XLAL_EFUNC, "create_gctFstat_toplist() failed for nCand=%d and sortBy=%d\n", uvar_nCand1, SORTBY_BSGL );
      XLAL_CHECK ( 0 == create_gctFstat_toplist ( &semiCohToplist2, uvar_nCand1, SORTBY_BSGLtL ),
                   XLAL_EFUNC, "create_gctFstat_toplist() failed for nCand=%d and sortBy=%d\n", uvar_nCand1, SORTBY_BSGLtL );
      XLAL_CHECK ( 0 == create_gctFstat_toplist ( &semiCohToplist3, uvar_nCand1, SORTBY_BtSGLtL ),
                   XLAL_EFUNC, "create_gctFstat_toplist() failed for nCand=%d and sortBy=%d\n", uvar_nCand1, SORTBY_BtSGLtL );
    }
  else if ( uvar_SortToplist == SORTBY_F_BSGLtL_BtSGLtL )// special treatement of 'triple' toplists: 1st one sorted by 'F', 2nd one by 'B_S/GLtL', 3rd by 'B_tS/GLtL'
    {
      XLAL_CHECK ( 0 == create_gctFstat_toplist ( &semiCohToplist, uvar_nCand1, SORTBY_F ),
                   XLAL_EFUNC, "create_gctFstat_toplist() failed for nCand=%d and sortBy=%d\n", uvar_nCand1, SORTBY_F );
      XLAL_CHECK ( 0 == create_gctFstat_toplist ( &semiCohToplist2, uvar_nCand1, SORTBY_BSGLtL ),
                   XLAL_EFUNC, "create_gctFstat_toplist() failed for nCand=%d and sortBy=%d\n", uvar_nCand1, SORTBY_BSGLtL );
      XLAL_CHECK ( 0 == create_gctFstat_toplist ( &semiCohToplist3, uvar_nCand1, SORTBY_BtSGLtL ),
                   XLAL_EFUNC, "create_gctFstat_toplist() failed for nCand=%d and sortBy=%d\n", uvar_nCand1, SORTBY_BtSGLtL );
    }
  else	// 'normal' single-sorting toplist cases (sortby 'F', 'nc' or 'BSGL')
    {
      XLAL_CHECK ( 0 == create_gctFstat_toplist ( &semiCohToplist, uvar_nCand1, uvar_SortToplist),
                   XLAL_EFUNC, "create_gctFstat_toplist() failed for nCand=%d and sortBy=%d\n", uvar_nCand1, uvar_SortToplist );
    }

  /* write the log file */
  if ( uvar_log )
    {
      fnamelog = LALCalloc( strlen(uvar_fnameout) + 1 + 4, sizeof(CHAR) );
      strcpy(fnamelog, uvar_fnameout);
      strcat(fnamelog, ".log");
      /* open the log file for writing */
      if ((fpLog = fopen(fnamelog, "wb")) == NULL) {
        fprintf(stderr, "Unable to open file %s for writing\n", fnamelog);
        LALFree(fnamelog);
        /*exit*/
        return(HIERARCHICALSEARCH_EFILE);
      }

      /* get the log string */
      XLAL_CHECK_MAIN( ( logstr = XLALUserVarGetLog(UVAR_LOGFMT_CFGFILE) ) != NULL, XLAL_EFUNC);

      fprintf( fpLog, "# Log file for HierarchSearchGCT.c\n\n");
      fprintf( fpLog, "# User Input:\n");
      fprintf( fpLog, "#-------------------------------------------\n");
      fprintf( fpLog, "# cmdline: %s\n", logstr );
      LALFree(logstr);

      /* add code version ID (only useful for git-derived versions) */
      fprintf ( fpLog, "# version: %s\n", VCSInfoString );

      fclose (fpLog);
      LALFree(fnamelog);

    } /* end of logging */

  tic_Start = GETTIME();
  usefulParams.collectFstatTiming = ( uvar_outputTimingDetails != NULL );

  /* initializations of coarse and fine grids */
  coarsegrid.TwoF=NULL;
  coarsegrid.TwoFX=NULL;
  coarsegrid.Uindex=NULL;
  finegrid.nc= NULL;
  finegrid.sumTwoF=NULL;
  finegrid.sumTwoFX=NULL;
  finegrid.maxTwoFl=NULL;
  finegrid.maxTwoFXl=NULL;
  finegrid.maxTwoFlIdx=NULL;
  finegrid.maxTwoFXlIdx=NULL;

  /* initialize ephemeris info */
  EphemerisData *edat;
  XLAL_CHECK ( (edat = XLALInitBarycenter ( uvar_ephemEarth, uvar_ephemSun )) != NULL, XLAL_EFUNC );

  XLALGPSSetREAL8(&minStartTimeGPS, uvar_minStartTime1);
  XLALGPSSetREAL8(&maxStartTimeGPS, uvar_maxStartTime1);

  /* create output files for writing if requested by user */
  if ( uvar_printCand1 )
    {
      fnameSemiCohCand = LALCalloc( strlen(uvar_fnameout) + 1, sizeof(CHAR) );
      if ( fnameSemiCohCand == NULL) {
        fprintf(stderr, "error allocating memory [HierarchSearchGCT.c %d]\n" , __LINE__);
        return(HIERARCHICALSEARCH_EMEM);
      }
      strcpy(fnameSemiCohCand, uvar_fnameout);
    }

  if ( uvar_printFstat1 )
    {
      const CHAR *append = "_fstatVec1.dat";
      fnameFstatVec1 = LALCalloc( strlen(uvar_fnameout) + strlen(append) + 1, sizeof(CHAR) );
      strcpy(fnameFstatVec1, uvar_fnameout);
      strcat(fnameFstatVec1, append);
      if ( !(fpFstat1 = fopen( fnameFstatVec1, "wb")))  {
        fprintf ( stderr, "Unable to open Fstat file fstatvec1.out for writing.\n");
        return (HIERARCHICALSEARCH_EFILE);
      }
    }

  /*------------ Set up stacks, detector states etc. */
  /* initialize spin range vectors */
  XLAL_INIT_MEM( spinRange_Temp );

  /* some useful first stage params */
  usefulParams.sftbasename = uvar_DataFiles1;

  /* ----- prepare generation of coherent segment list */
  if ( XLALUserVarWasSet ( &uvar_segmentList ) )
    {
      if ( XLALUserVarWasSet ( &uvar_nStacksMax ) || XLALUserVarWasSet ( &uvar_tStack ) ) {
        XLALPrintError ( "Use EITHER (--nStacksMax and --tStack) OR --segmentList to define the coherent segments!\n\n" );
        return HIERARCHICALSEARCH_EBAD;
      }
      if ( (usefulParams.segmentList = XLALReadSegmentsFromFile ( uvar_segmentList )) == NULL ) {
        XLALPrintError ("Failed to parse segment-list file '%s'. xlalErrno = %d\n", uvar_segmentList, xlalErrno );
        return HIERARCHICALSEARCH_ESUB;
      }
    }
  else /* set up maximally nStacksMax fixed-size segments of length tStack */
    {
      if ( !XLALUserVarWasSet ( &uvar_tStack ) ) {
        XLALPrintError ( "Need to set --tStack or --segmentList to define the coherent segments!\n\n" );
        return HIERARCHICALSEARCH_EBAD;
      }

      usefulParams.nStacks = uvar_nStacksMax;
      usefulParams.tStack = uvar_tStack;
      usefulParams.segmentList = NULL;
    }
  /* ----- */


  XLAL_INIT_MEM ( usefulParams.spinRange_startTime );
  XLAL_INIT_MEM ( usefulParams.spinRange_endTime );
  XLAL_INIT_MEM ( usefulParams.spinRange_refTime );
  XLAL_INIT_MEM ( usefulParams.spinRange_midTime );

  /* copy user specified spin variables at reftime  */
  /* the reference time value in spinRange_refTime will be set in SetUpSFTs() */
  usefulParams.spinRange_refTime.fkdot[0] = uvar_Freq; /* frequency */
  usefulParams.spinRange_refTime.fkdot[1] = uvar_f1dot;  /* 1st spindown */
  usefulParams.spinRange_refTime.fkdot[2] = uvar_f2dot;  /* 2nd spindown */
  usefulParams.spinRange_refTime.fkdot[3] = uvar_f3dot;  /* 3rd spindown */  
  usefulParams.spinRange_refTime.fkdotBand[0] = uvar_FreqBand; /* frequency range */
  usefulParams.spinRange_refTime.fkdotBand[1] = uvar_f1dotBand; /* spindown range */
  usefulParams.spinRange_refTime.fkdotBand[2] = uvar_f2dotBand; /* spindown range */
  usefulParams.spinRange_refTime.fkdotBand[3] = uvar_f3dotBand; /* spindown range */

  usefulParams.edat = edat;
  usefulParams.minStartTimeGPS = minStartTimeGPS;
  usefulParams.maxStartTimeGPS = maxStartTimeGPS;
  usefulParams.blocksRngMed = uvar_blocksRngMed;
  usefulParams.Dterms = uvar_Dterms;
  usefulParams.DtermsRecalc = uvar_DtermsRecalc;
  usefulParams.assumeSqrtSX = uvar_assumeSqrtSX;
  usefulParams.SSBprec = uvar_SSBprecision;
  usefulParams.Fmethod = uvar_FstatMethod;
  usefulParams.FmethodRecalc = uvar_FstatMethodRecalc;
  usefulParams.recalcToplistStats = uvar_recalcToplistStats;

  usefulParams.mismatch1 = uvar_mismatch1;

  /* set reference time for pulsar parameters */
  if ( XLALUserVarWasSet(&uvar_refTime))
    usefulParams.refTime = uvar_refTime;
  else {
    XLALPrintWarning("Reference time will be set to mid-time of observation time\n");
    usefulParams.refTime = -1;
  }

  /* set Fstat calculation frequency resolution (coarse grid) */
  if ( XLALUserVarWasSet(&uvar_FreqBand) ) {
    if ( XLALUserVarWasSet(&uvar_dFreq) ) {
      usefulParams.dFreqStack = uvar_dFreq;
    } else {
      XLALPrintError("--dFreq is required if --FreqBand is given\n");
      return( HIERARCHICALSEARCH_EBAD );
    }
  } else {
    usefulParams.dFreqStack = 0;
  }

  /* set Fstat spindown resolution (coarse grid) */
  if ( XLALUserVarWasSet(&uvar_f1dotBand) ) {
    if ( XLALUserVarWasSet(&uvar_df1dot) ) {
      usefulParams.df1dot = uvar_df1dot;
    } else {
      XLALPrintError("--df1dot is required if --f1dotBand is given\n");
      return( HIERARCHICALSEARCH_EBAD );
    }
  } else {
    usefulParams.df1dot = 0;
  }

  /* set Fstat 2nd spindown resolution (coarse grid) */
  if ( XLALUserVarWasSet(&uvar_f2dotBand) ) {
    if ( XLALUserVarWasSet(&uvar_df2dot) ) {
      usefulParams.df2dot = uvar_df2dot;
    }
    else {
      XLALPrintError("--df2dot is required if --f2dotBand is given\n");
      return( HIERARCHICALSEARCH_EBAD );
    }
  }
  else {
    usefulParams.df2dot = 0;
  }
  
  /* set Fstat 3rd spindown resolution (coarse grid) */
  if ( XLALUserVarWasSet(&uvar_f3dotBand) ) {
    if ( XLALUserVarWasSet(&uvar_df3dot) ) {
      usefulParams.df3dot = uvar_df3dot;
    }
    else {
      XLALPrintError("--df3dot is required if --f3dotBand is given\n");
      return( HIERARCHICALSEARCH_EBAD );
    }
  }
  else {
    usefulParams.df3dot = 0;
  }

  // read signal parameters to be injected, if requested by the user

  if ( uvar_injectionSources != NULL ) {
    XLAL_CHECK_MAIN ( (usefulParams.injectionSources = XLALPulsarParamsFromUserInput ( uvar_injectionSources, NULL ) ) != NULL, XLAL_EFUNC );
  }

  /* for 1st stage: read sfts, calculate detector states */
  LogPrintf( LOG_NORMAL,"Reading input data ... ");
  LAL_CALL( SetUpSFTs( &status, &usefulParams ), &status);
  LogPrintfVerbatim ( LOG_NORMAL, " done.\n");

  /* some useful params computed by SetUpSFTs */
  tStack = usefulParams.tStack;
  tObs = usefulParams.tObs;
  nStacks = usefulParams.nStacks;
  midTstack = usefulParams.midTstack;
  startTstack = usefulParams.startTstack;
  endTstack = usefulParams.endTstack;
  tMidGPS = usefulParams.spinRange_midTime.refTime;
  refTimeGPS = usefulParams.spinRange_refTime.refTime;
  fprintf(stderr, "%% --- GPS reference time = %.4f ,  GPS data mid time = %.4f\n",
          XLALGPSGetREAL8(&refTimeGPS), XLALGPSGetREAL8(&tMidGPS) );

  REAL4 *loudestTwoFPerSeg = NULL;
  if ( uvar_loudestTwoFPerSeg )
    {
      loudestTwoFPerSeg = XLALCalloc ( nStacks, sizeof(REAL4) );
    } // if loudestTwoFPerSeg

  /* free segment list */
  if ( usefulParams.segmentList )
    if ( XLALSegListClear( usefulParams.segmentList ) != XLAL_SUCCESS )
      XLAL_ERROR ( XLAL_EFUNC );
  XLALFree ( usefulParams.segmentList );
  usefulParams.segmentList = NULL;

  /*------- set frequency and spindown resolutions and ranges for Fstat and semicoherent steps -----*/

  dFreqStack = usefulParams.dFreqStack;
  df1dot = usefulParams.df1dot;
  df2dot = usefulParams.df2dot;
  df3dot = usefulParams.df3dot;
  LogPrintf(LOG_NORMAL, "dFreqStack = %e, df1dot = %e, df2dot = %e, df3dot = %e\n", dFreqStack, df1dot, df2dot, df3dot);

  /* set number of fine-grid spindowns */
  if ( XLALUserVarWasSet(&uvar_gammaRefine) ) {
    gammaRefine = uvar_gammaRefine;
  }
  else {
    sigmasq = 0.0; /* second moment of segments' midpoints */
    for (k = 0; k < nStacks; k++) {
      midTstackGPS = midTstack->data[k];
      timeDiffSeg = XLALGPSDiff( &midTstackGPS, &tMidGPS );
      sigmasq = sigmasq + (timeDiffSeg * timeDiffSeg);
    }
    sigmasq = sigmasq / (nStacks * tStack * tStack);
    /* Refinement factor (approximate) */
    gammaRefine = sqrt(1.0 + 60 * sigmasq);   /* Eq. from PRL, page 3 */
  }

  /* set number of fine-grid 2nd spindowns */
  if ( XLALUserVarWasSet(&uvar_gamma2Refine) ) {
    /* use 1st spindown refinement if user value < 0 */
    if ( uvar_gamma2Refine < 0 ) {
      gamma2Refine = gammaRefine;
    }
    else {
      gamma2Refine = uvar_gamma2Refine;
    }
  }
  else {
    sigmasq = sigma4 = 0.0; /* second and 4th moment of segments' midpoints */
    for (k = 0; k < nStacks; k++) {
      midTstackGPS = midTstack->data[k];
      timeDiffSeg = XLALGPSDiff( &midTstackGPS, &tMidGPS );
      sigmasq = sigmasq + (timeDiffSeg * timeDiffSeg);
      sigma4  = sigma4  + (timeDiffSeg * timeDiffSeg * timeDiffSeg * timeDiffSeg);
    }
    sigmasq = sigmasq / (nStacks * tStack * tStack);
    sigma4  = sigma4  / (nStacks * tStack * tStack * tStack * tStack);
    /* 2nd spindown refinement factor gamma2/gamma1 (approximate)
       see Pletsch, PRD 82 042002, 2010 */
    gamma2Refine = sqrt( 2100.0 * (sigma4 - sigmasq*sigmasq) );
  }

  /**** debugging information ******/
  /* print some debug info about spinrange */
  LogPrintf(LOG_DETAIL, "Frequency and spindown range at refTime (%d): [%f,%f], [%e,%e], [%e,%e], [%e,%e]\n",
            usefulParams.spinRange_refTime.refTime.gpsSeconds,
            usefulParams.spinRange_refTime.fkdot[0],
            usefulParams.spinRange_refTime.fkdot[0] + usefulParams.spinRange_refTime.fkdotBand[0],
            usefulParams.spinRange_refTime.fkdot[1],
            usefulParams.spinRange_refTime.fkdot[1] + usefulParams.spinRange_refTime.fkdotBand[1],
            usefulParams.spinRange_refTime.fkdot[2],
            usefulParams.spinRange_refTime.fkdot[2] + usefulParams.spinRange_refTime.fkdotBand[2],
			usefulParams.spinRange_refTime.fkdot[3],
            usefulParams.spinRange_refTime.fkdot[3] + usefulParams.spinRange_refTime.fkdotBand[3]);
  
  LogPrintf(LOG_DETAIL, "Frequency and spindown range at startTime (%d): [%f,%f], [%e,%e], [%e,%e], [%e,%e]\n",
            usefulParams.spinRange_startTime.refTime.gpsSeconds,
            usefulParams.spinRange_startTime.fkdot[0],
            usefulParams.spinRange_startTime.fkdot[0] + usefulParams.spinRange_startTime.fkdotBand[0],
            usefulParams.spinRange_startTime.fkdot[1],
            usefulParams.spinRange_startTime.fkdot[1] + usefulParams.spinRange_startTime.fkdotBand[1],
            usefulParams.spinRange_startTime.fkdot[2],
            usefulParams.spinRange_startTime.fkdot[2] + usefulParams.spinRange_startTime.fkdotBand[2],
            usefulParams.spinRange_startTime.fkdot[3],
            usefulParams.spinRange_startTime.fkdot[3] + usefulParams.spinRange_startTime.fkdotBand[3]);
  

  LogPrintf(LOG_DETAIL, "Frequency and spindown range at midTime (%d): [%f,%f], [%e,%e], [%e,%e], [%e,%e]\n",
            usefulParams.spinRange_midTime.refTime.gpsSeconds,
            usefulParams.spinRange_midTime.fkdot[0],
            usefulParams.spinRange_midTime.fkdot[0] + usefulParams.spinRange_midTime.fkdotBand[0],
            usefulParams.spinRange_midTime.fkdot[1],
            usefulParams.spinRange_midTime.fkdot[1] + usefulParams.spinRange_midTime.fkdotBand[1],
            usefulParams.spinRange_midTime.fkdot[2],
            usefulParams.spinRange_midTime.fkdot[2] + usefulParams.spinRange_midTime.fkdotBand[2],
			usefulParams.spinRange_midTime.fkdot[3],
            usefulParams.spinRange_midTime.fkdot[3] + usefulParams.spinRange_midTime.fkdotBand[3]);

  LogPrintf(LOG_DETAIL, "Frequency and spindown range at endTime (%d): [%f,%f], [%e,%e], [%e,%e], [%e,%e]\n",
            usefulParams.spinRange_endTime.refTime.gpsSeconds,
            usefulParams.spinRange_endTime.fkdot[0],
            usefulParams.spinRange_endTime.fkdot[0] + usefulParams.spinRange_endTime.fkdotBand[0],
            usefulParams.spinRange_endTime.fkdot[1],
            usefulParams.spinRange_endTime.fkdot[1] + usefulParams.spinRange_endTime.fkdotBand[1],
            usefulParams.spinRange_endTime.fkdot[2],
            usefulParams.spinRange_endTime.fkdot[2] + usefulParams.spinRange_endTime.fkdotBand[2],
            usefulParams.spinRange_endTime.fkdot[3],
            usefulParams.spinRange_endTime.fkdot[3] + usefulParams.spinRange_endTime.fkdotBand[3]);

  /* print debug info about stacks */
  fprintf(stderr, "%% --- Setup, N = %d, T = %.0f s, Tobs = %.0f s, gammaRefine = %.0f, gamma2Refine = %.0f, gamma3Refine = %.0f\n",
          nStacks, tStack, tObs, gammaRefine, gamma2Refine, gamma3Refine);


  /*---------- set up F-statistic calculation stuff ---------*/
  /* set reference time for calculating Fstatistic */
  thisPoint.refTime = tMidGPS; /* midpoint of data spanned */

  /* binary orbit and higher spindowns not considered */
  thisPoint.asini = 0 /* isolated pulsar */;
  XLAL_INIT_MEM ( thisPoint.fkdot );

  /*---------- set up stuff for semi-coherent part ---------*/
  /* set up some semiCoherent parameters */
  semiCohPar.tsMid = midTstack;
  semiCohPar.refTime = tMidGPS;	// unused??

  /* reference time for finegrid is midtime */
  finegrid.refTime = tMidGPS;

  /* allocate memory for pos vel acc vectors */
  posStack = XLALCreateREAL8VectorSequence( nStacks, 3 );
  velStack = XLALCreateREAL8VectorSequence( nStacks, 3 );
  accStack = XLALCreateREAL8VectorSequence( nStacks, 3 );

  /* calculate Earth orbital positions, velocities and accelerations */
  LAL_CALL( GetSegsPosVelAccEarthOrb( &status, &posStack, &velStack, &accStack, &usefulParams), &status);

  /* semicoherent configuration parameters */
  semiCohPar.pos = posStack;
  semiCohPar.vel = velStack;
  semiCohPar.acc = accStack;
  semiCohPar.outBaseName = uvar_fnameout;

  /* get number of detectors and detector name vector */
  LALStringVector *detectorIDs = usefulParams.detectorIDs;
  const UINT4 numDetectors = detectorIDs->length;

  // compute single-IFO F-statistics for line-robust stats
  if ( uvar_computeBSGL )
    {
      Fstat_what |= FSTATQ_2F_PER_DET;

      /* take BSGL user vars and pre-compute the corresponding BSGLsetup */
      REAL4 *oLGX_p = NULL;
      REAL4 oLGX[PULSAR_MAX_DETECTORS];
      if ( uvar_oLGX != NULL )
        {
          if ( uvar_oLGX->length != numDetectors ) {
            fprintf ( stderr, "Invalid input: length(oLGX) = %d differs from number of detectors (%d)'\n", uvar_oLGX->length, numDetectors );
            return( HIERARCHICALSEARCH_EBAD );
          }
          if ( XLALParseLinePriors ( &oLGX[0], uvar_oLGX ) != XLAL_SUCCESS ) {
            fprintf(stderr, "Invalid input oLGX'\n" );
            return( HIERARCHICALSEARCH_EBAD );
          }
          oLGX_p = &oLGX[0];
        } // if uvar_oLGX != NULL

      usefulParams.BSGLsetup = XLALCreateBSGLSetup ( numDetectors, uvar_Fstar0sc, oLGX_p, uvar_BSGLlogcorr, nStacks );
      if ( usefulParams.BSGLsetup == NULL ) {
        fprintf(stderr, "XLALCreateBSGLSetup() failed\n");
        return( HIERARCHICALSEARCH_EBAD );
      }
    } // if uvar_computeBSGL

  /* assemble column headings string for output file */
  CHAR column_headings_string_base[256];
  if (XLALUserVarWasSet(&uvar_f3dot)) {
		sprintf(column_headings_string_base,"freq alpha delta f1dot f2dot f3dot nc <2F>");
  }
  else {
		sprintf(column_headings_string_base,"freq alpha delta f1dot f2dot nc <2F>");
  }

  UINT4 column_headings_string_length = sizeof(column_headings_string_base);
  if ( uvar_computeBSGL ) {
    column_headings_string_length += 10 + numDetectors*8; /* 10 for " log10BSGL" and 8 per detector for " <2F_XY>" */
  }
  if ( uvar_getMaxFperSeg ) {
    column_headings_string_length += 12 + 13; /* 12 for " log10BSGLtL" and 13 for " log10BtSGLtL" */
    column_headings_string_length += 6 +9 + numDetectors*(9+12); /* 6 for " max2F", 9 for " max2Fseg" and 9+12 per detector for " max2F_XY" and " max2F_XYseg"*/
  }
  if ( uvar_recalcToplistStats ) {
    column_headings_string_length += 6 + 11 + numDetectors*9; /* 6 for " <2Fr>" and 9 per detector for " <2Fr_XY>" */
    if ( uvar_computeBSGL) {
      column_headings_string_length += 11; /* for " log10BSGLr" */
      if ( uvar_getMaxFperSeg ) {
        column_headings_string_length += 13; /* for " log10BSGLtLr" */
      }
    }
    if (XLALUserVarWasSet(&uvar_f3dot)) {
      column_headings_string_length += 1;
    }
    if ( uvar_loudestSegOutput ) {
      column_headings_string_length += 9 + numDetectors*7; /* for " lseg 2Fl" and 7 per detector for " 2Fl_XY" */
    }
  }
  char column_headings_string[column_headings_string_length];
  XLAL_INIT_MEM( column_headings_string );
  strcat ( column_headings_string, column_headings_string_base );
  if ( uvar_computeBSGL ) {
    strcat ( column_headings_string, " log10BSGL" );
    for ( UINT4 X = 0; X < numDetectors ; X ++ ) {
      char headingX[9];
      snprintf ( headingX, sizeof(headingX), " <2F_%s>", detectorIDs->data[X] );
      strcat ( column_headings_string, headingX );
    } /* for X < numDet */
  }
  if ( uvar_getMaxFperSeg ) {
    strcat ( column_headings_string, " log10BSGLtL log10BtSGLtL" );
    strcat ( column_headings_string, " max2F" );
    strcat ( column_headings_string, " max2Fseg" );
    for ( UINT4 X = 0; X < numDetectors ; X ++ ) {
      char headingX[16];
      snprintf ( headingX, sizeof(headingX), " max2F_%s", detectorIDs->data[X] );
      strcat ( column_headings_string, headingX );
    }
    for ( UINT4 X = 0; X < numDetectors ; X ++ ) {
      char headingX[16];
      snprintf ( headingX, sizeof(headingX), " max2F_%sseg", detectorIDs->data[X] );
      strcat ( column_headings_string, headingX );
    } /* for X < numDet */
  }
  if ( uvar_recalcToplistStats ) {
    strcat ( column_headings_string, " <2Fr>" );
    if ( uvar_computeBSGL) {
      strcat ( column_headings_string, " log10BSGLr" );
    }
    for ( UINT4 X = 0; X < numDetectors ; X ++ ) {
      char headingX[10];
      snprintf ( headingX, sizeof(headingX), " <2Fr_%s>", detectorIDs->data[X] );
      strcat ( column_headings_string, headingX );
    } /* for X < numDet */
    if ( uvar_loudestSegOutput ) {
      strcat ( column_headings_string, " lseg 2Fl" );
      for ( UINT4 X = 0; X < numDetectors ; X ++ ) {
        char headingX[10];
        snprintf ( headingX, sizeof(headingX), " 2Fl_%s", detectorIDs->data[X] );
        strcat ( column_headings_string, headingX );
      }
    }
    if ( uvar_computeBSGL && uvar_getMaxFperSeg ) {
      strcat ( column_headings_string, " log10BSGLtLr" );
    }
  }
  global_column_headings_stringp = column_headings_string;

  /* get effective inverse number of segments per detector (needed for correct averaging in single-IFO F calculation) */
  REAL4 NSegmentsInv = 1.0 / nStacks; /* also need this for multi-detector F-stat averaging later on */

  /*-----------Create template grid for first stage ---------------*/
  /* prepare initialization of DopplerSkyScanner to step through paramter space */
  scanInit.dAlpha = uvar_dAlpha;
  scanInit.dDelta = uvar_dDelta;
  scanInit.gridType = uvar_gridType1;
  scanInit.metricType = uvar_metricType1;
  scanInit.metricMismatch = uvar_mismatch1;
  scanInit.projectMetric = TRUE;
  scanInit.obsDuration = tStack;
  scanInit.obsBegin = tMidGPS;
  scanInit.ephemeris = edat;
  scanInit.skyGridFile = uvar_skyGridFile;
  scanInit.skyRegionString = (CHAR*)LALCalloc(1, strlen(uvar_skyRegion)+1);
  if ( scanInit.skyRegionString == NULL) {
    fprintf(stderr, "error allocating memory [HierarchSearchGCT.c.c %d]\n" , __LINE__);
    return(HIERARCHICALSEARCH_EMEM);
  }
  strcpy (scanInit.skyRegionString, uvar_skyRegion);

  // just use first SFTs' IFO for metric (should be irrelevant)
  const LALDetector* firstDetector = XLALGetSiteInfo( detectorIDs->data[0] );
  if ( firstDetector == NULL ) {
    LogPrintf ( LOG_CRITICAL, "\nXLALGetSiteInfo() failed for detector '%s'\n", detectorIDs->data[0] );
    return HIERARCHICALSEARCH_EXLAL;
  }
  scanInit.Detector = firstDetector;

  scanInit.numSkyPartitions = uvar_numSkyPartitions;
  scanInit.partitionIndex = uvar_partitionIndex;

  scanInit.Freq = usefulParams.spinRange_midTime.fkdot[0] + usefulParams.spinRange_midTime.fkdotBand[0]; /* Only used if metric is employed */

  /* initialize skygrid  */
  LogPrintf(LOG_DETAIL, "Setting up coarse sky grid... \n");
  LAL_CALL ( InitDopplerSkyScan ( &status, &thisScan, &scanInit), &status);


  /* ----- start main calculations by going over coarse grid points --------*/
  skyGridCounter = 0;
  f1dotGridCounter = 0;

  XLALNextDopplerSkyPos(&dopplerpos, &thisScan);

  /* "spool forward" if we found a checkpoint */
  {
    UINT4 count = 0; /* The first checkpoint should have value 1 */
    UINT4 skycount = 0;

    if( 0 > GET_GCT_CHECKPOINT (uvar_fnameChkPoint, semiCohToplist, semiCohToplist2, semiCohToplist3, &count)) {
      XLALPrintError ("%s : '%s' \n", HIERARCHICALSEARCH_MSGECHECKPT,uvar_fnameChkPoint);
      return (HIERARCHICALSEARCH_ECHECKPT); 
    }

    if (count) {
      f1dotGridCounter = (UINT4) (count % usefulParams.nf1dot);  /* Checkpointing counter = i_sky * nf1dot + i_f1dot */
      skycount = (UINT4) ((count - f1dotGridCounter) / usefulParams.nf1dot);
    }
   fprintf (stderr, "%% --- Cpt:%d,  total:%d,  sky:%d/%d,  f1dot:%d/%d\n",
             count, thisScan.numSkyGridPoints*usefulParams.nf1dot, skycount+1, thisScan.numSkyGridPoints, f1dotGridCounter+1, usefulParams.nf1dot);

    for(skyGridCounter = 0; skyGridCounter < skycount; skyGridCounter++)
      XLALNextDopplerSkyPos(&dopplerpos, &thisScan);

    if ( count == thisScan.numSkyGridPoints * usefulParams.nf1dot )
      thisScan.state = STATE_FINISHED;

  }

  /* spool forward if uvar_skyPointIndex is set
     ---This probably doesn't make sense when checkpointing is turned on */
  if ( XLALUserVarWasSet(&uvar_skyPointIndex)) {
    UINT4 count = uvar_skyPointIndex;
    for(skyGridCounter = 0; (skyGridCounter < count)&&(thisScan.state != STATE_FINISHED) ; skyGridCounter++)
      XLALNextDopplerSkyPos(&dopplerpos, &thisScan);
  }

  /* ################## loop over SKY coarse-grid points ################## */
  while(thisScan.state != STATE_FINISHED)
    {
      /*------------- calculate F-Statistic for each segment --------------*/

      /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
      thisPoint.Alpha = dopplerpos.Alpha;
      thisPoint.Delta = dopplerpos.Delta;

      /* Calculate unit vector n */
      cosAlpha = cos(thisPoint.Alpha);
      sinAlpha = sin(thisPoint.Alpha);
      cosDelta = cos(thisPoint.Delta);
      sinDelta = sin(thisPoint.Delta);
      nvec[0] = cosAlpha * cosDelta;
      nvec[1] = sinAlpha * cosDelta;
      nvec[2] = sinDelta;

        /* calculate number of bins for Fstat overhead due to residual spin-down */
        semiCohPar.extraBinsFstat = usefulParams.extraBinsFstat;
        binsFstat1 = usefulParams.binsFstatSearch + 2 * semiCohPar.extraBinsFstat;

      /* ################## loop over coarse-grid F1DOT values ################## */
      ifdot = 0;

      while ( ifdot < usefulParams.nf1dot ) {

        /* if checkpoint read, spool forward */
        if (f1dotGridCounter > 0) {
          ifdot = f1dotGridCounter;
          f1dotGridCounter = 0;
        }

        /* ################## loop over coarse-grid F2DOT values ################## */
        if2dot = 0;

        while ( if2dot < usefulParams.nf2dot ) {

        /* ################## loop over coarse-grid F3DOT values ################## */
        if3dot = 0;

        while ( if3dot < usefulParams.nf3dot ) {

          /* show progress */
          LogPrintf( LOG_NORMAL, "Coarse grid sky:%d/%d f1dot:%d/%d f2dot:%d/%d f3dot:%d/%d\n",
                     skyGridCounter+1, thisScan.numSkyGridPoints, ifdot+1, usefulParams.nf1dot, if2dot+1, usefulParams.nf2dot, if3dot+1, usefulParams.nf3dot );


          /* ------------- Set up coarse grid --------------------------------------*/
          coarsegrid.freqlength = (UINT4) (binsFstat1);
          coarsegrid.nStacks = nStacks;
          coarsegrid.length = coarsegrid.freqlength * coarsegrid.nStacks;
          coarsegrid.numDetectors = numDetectors;

          /* allocate memory for coarsegrid */
          coarsegrid.TwoF = (REAL4 *)LALRealloc( coarsegrid.TwoF, coarsegrid.length * sizeof(REAL4));
          if ( uvar_computeBSGL ) {
            coarsegrid.TwoFX = (REAL4 *)LALRealloc( coarsegrid.TwoFX, coarsegrid.numDetectors * coarsegrid.length * sizeof(REAL4));
          }
          coarsegrid.Uindex = (UINT4 *)LALRealloc( coarsegrid.Uindex, coarsegrid.length * sizeof(UINT4));

          if ( coarsegrid.TwoF == NULL || coarsegrid.Uindex == NULL) {
            fprintf(stderr, "ERROR: Memory allocation  [HierarchSearchGCT.c %d]\n" , __LINE__);
            return(HIERARCHICALSEARCH_EMEM);
          }

          /* ------------- Set up fine grid --------------------------------------*/

          /* frequency fine-grid borders */
          freqmin_fg = usefulParams.spinRange_midTime.fkdot[0];
          freqband_fg = usefulParams.spinRange_midTime.fkdotBand[0];

          /* fine-grid frequency resolution */
          dfreq_fg = dFreqStack;
          UINT4 nfreqs_fg = ceil(freqband_fg / dfreq_fg);  /* number of points in frequency */
          if (nfreqs_fg == 0) {
            nfreqs_fg++;
          }

          /* copy frequency setup parameters to fine-grid struct */
          finegrid.freqmin_fg = freqmin_fg;
          finegrid.dfreq_fg = dfreq_fg;
          finegrid.freqlength = nfreqs_fg ;
#define ALIGN_REAL4 4  /* 16 bytes / sizeof(REAL4) = 4 */
          finegrid.freqlengthAL = ALIGN_REAL4 * ((UINT4)ceil ( 1.0 * finegrid.freqlength / ALIGN_REAL4 ));

          /* fine-grid f1dot resolution */
          if (usefulParams.nf1dot == 1) {
            nf1dots_fg = 1;
          }
          else {
            nf1dots_fg = ceil(gammaRefine);        /* number of spindown fine-grid points */
            if ( (nf1dots_fg % 2) == 0 ) {    /* if even, add one (to refine symmetrically) */
              nf1dots_fg++;
            }
          }
          df1dot_fg = df1dot / nf1dots_fg;  /* spindown fine-grid  stepsize */

          /* adjust f1dotmin_fg, so that f1dot finegrid is centered around coarse-grid f1dot point */
          f1dotmin_fg = (usefulParams.spinRange_midTime.fkdot[1] + ifdot * df1dot) - df1dot_fg * floor(nf1dots_fg / 2.0);

          /* fine-grid f2dot resolution */
          if ( uvar_f2dotBand == 0 ) {
            nf2dots_fg = 1;
          }
          else {
            nf2dots_fg = ceil(gamma2Refine);        /* number of 2nd spindown fine-grid points */
            if ( (nf2dots_fg % 2) == 0 ) {    /* if even, add one (to refine symmetrically) */
              nf2dots_fg++;
            }
          }
          df2dot_fg = df2dot / nf2dots_fg;  /* 2nd spindown fine-grid  stepsize */

          /* adjust f2dotmin_fg, so that f2dot finegrid is centered around coarse-grid f2dot point */
          f2dotmin_fg = (usefulParams.spinRange_midTime.fkdot[2] + if2dot * df2dot) - df2dot_fg * floor(nf2dots_fg / 2.0);

          /* fine-grid f3dot resolution */
          nf3dots_fg = 1;
          df3dot_fg = 1;  /* 3rd spindown fine-grid  stepsize */

          /* adjust f3dotmin_fg, so that f3dot finegrid is centered around coarse-grid f3dot point */
          f3dotmin_fg = (usefulParams.spinRange_midTime.fkdot[3] + if3dot * df3dot) - df3dot_fg * floor(nf3dots_fg / 2.0);

          /* total number of fine-grid points */
          finegrid.length = finegrid.freqlength;

          if(!oldcg) {
            oldcg = coarsegrid.length;
            LogPrintfVerbatim(LOG_NORMAL, "%% --- CG:%d ",coarsegrid.length);
          }
          if(!oldfg) {
            oldfg = finegrid.length;
            LogPrintfVerbatim(LOG_NORMAL, "FG:%d  f1dotmin_fg:%.13g df1dot_fg:%.13g f2dotmin_fg:%.13g df2dot_fg:%.13g f3dotmin_fg:%.13g df3dot_fg:%.13g\n",
                              finegrid.length,f1dotmin_fg,df1dot_fg,f2dotmin_fg,df2dot_fg,f3dotmin_fg,df3dot_fg);
          }
          if((coarsegrid.length != oldcg) || (finegrid.length != oldfg)) {
            LogPrintfVerbatim(LOG_CRITICAL, "ERROR: Grid-sizes disagree!\nPrevious CG:%d FG:%d, currently CG:%d FG:%d\n",
                              oldcg,oldfg,coarsegrid.length,finegrid.length);
            return(HIERARCHICALSEARCH_EVAL);
          }

          /* number of detectors, needed for sumTwoFX array */
          finegrid.numDetectors = coarsegrid.numDetectors;

          /* allocate memory for finegrid points */
          /* FIXME: The SSE2 optimized code relies on an identical alignment modulo 16 bytes of the
             arrays finegrid.nc and finegrid.sumTwoF !!
             MacOS enforces 16 byte alignment of memory blocks with size >= 16 bytes allocated
             with malloc and realloc, but e.g. for 32 bit Linux this will NOT hold.
             Alternatives might be using (posix_)memalign under Linux.
             Windows ==>???
          */

          finegrid.nc = (FINEGRID_NC_T *)ALRealloc( finegrid.nc, finegrid.length * sizeof(FINEGRID_NC_T));
          finegrid.sumTwoF = (REAL4 *)ALRealloc( finegrid.sumTwoF, finegrid.length * sizeof(REAL4));
          if ( uvar_getMaxFperSeg ) {
            finegrid.maxTwoFl = (REAL4 *)ALRealloc( finegrid.maxTwoFl, finegrid.length * sizeof(REAL4));
            finegrid.maxTwoFlIdx = (UINT4 *)ALRealloc( finegrid.maxTwoFlIdx, finegrid.length * sizeof(UINT4));

          }
          if ( uvar_computeBSGL ) {
            finegrid.sumTwoFX = (REAL4 *)ALRealloc( finegrid.sumTwoFX, finegrid.numDetectors * finegrid.freqlengthAL * sizeof(REAL4));
            if ( uvar_getMaxFperSeg ) {
              finegrid.maxTwoFXl = (REAL4 *)ALRealloc( finegrid.maxTwoFXl, finegrid.numDetectors * finegrid.freqlengthAL * sizeof(REAL4));
              finegrid.maxTwoFXlIdx = (UINT4 *)ALRealloc( finegrid.maxTwoFXlIdx, finegrid.numDetectors * finegrid.freqlengthAL * sizeof(UINT4));
            }
          }

          if ( finegrid.nc == NULL || finegrid.sumTwoF == NULL) {
            fprintf(stderr, "ERROR: Memory allocation [HierarchSearchGCT.c %d]\n" , __LINE__);
            return(HIERARCHICALSEARCH_EMEM);
          }

          /* copy sky coarse-grid point to finegrid, because sky is not refined */
          finegrid.alpha = thisPoint.Alpha;
          finegrid.delta = thisPoint.Delta;

          /* ---------- Walk through fine grid f1dot --------------- */
          for( if1dot_fg = 0; if1dot_fg < nf1dots_fg; if1dot_fg++ ) {

            /* get the 1st spindown of this fine-grid point */
            f1dot_fg = f1dotmin_fg + if1dot_fg * df1dot_fg;

            /* ---------- Walk through fine grid f2dot --------------- */
            for( if2dot_fg = 0; if2dot_fg < nf2dots_fg; if2dot_fg++ ) {

              /* get the 2nd spindown of this fine-grid point */
              f2dot_fg = f2dotmin_fg + if2dot_fg * df2dot_fg;

            /* ---------- Walk through fine grid f3dot --------------- */
            for( if3dot_fg = 0; if3dot_fg < nf3dots_fg; if3dot_fg++ ) {

              /* get the 3rd spindown of this fine-grid point */
              f3dot_fg = f3dotmin_fg + if3dot_fg * df3dot_fg;

              /* initialize the entire finegrid ( 2F-sum and number count set to 0 ) */
              memset( finegrid.nc, 0, finegrid.length * sizeof(FINEGRID_NC_T) );
              memset( finegrid.sumTwoF, 0, finegrid.length * sizeof(REAL4) );
              if ( uvar_getMaxFperSeg ) {
                memset( finegrid.maxTwoFl, 0, finegrid.length * sizeof(REAL4) );
		memset( finegrid.maxTwoFlIdx,0,finegrid.length * sizeof(UINT4) );
              }
              if ( uvar_computeBSGL ) {
                memset( finegrid.sumTwoFX, 0, finegrid.numDetectors * finegrid.freqlengthAL * sizeof(REAL4) );
                if ( uvar_getMaxFperSeg ) {
                  memset( finegrid.maxTwoFXl, 0, finegrid.numDetectors * finegrid.freqlengthAL * sizeof(REAL4) );
                  memset( finegrid.maxTwoFXlIdx, 0, finegrid.numDetectors * finegrid.freqlengthAL * sizeof(UINT4) );
                }
              }

              /* compute F-statistic values for coarse grid the first time through fine grid fdots loop */
              const BOOLEAN doComputeFstats = ( (if1dot_fg == 0) && (if2dot_fg == 0) && (if3dot_fg == 0));

              /* #########################################################################*/
              /* ------------- MAIN LOOP over Segments for F-statistic -------------------*/

              for (k = 0; k < nStacks; k++) {

                /* Get pos vel acc for this segment */
                pos[0] = semiCohPar.pos->data[3*k];
                pos[1] = semiCohPar.pos->data[3*k + 1];
                pos[2] = semiCohPar.pos->data[3*k + 2];

                vel[0] = semiCohPar.vel->data[3*k];
                vel[1] = semiCohPar.vel->data[3*k + 1];
                vel[2] = semiCohPar.vel->data[3*k + 2];

                /* currently unused:
                REAL8 acc[3];
                acc[0] = semiCohPar.acc->data[3*k];
                acc[1] = semiCohPar.acc->data[3*k + 1];
                acc[2] = semiCohPar.acc->data[3*k + 2];
                */

                /* Midpoint in time of current segment */
                midTstackGPS = midTstack->data[k];

                /* Difference in time between this segment's midpoint and Fstat reftime */
                timeDiffSeg = XLALGPSDiff( &midTstackGPS, &thisPoint.refTime );

                /* ---------------------------------------------------------------------------------------- */

                /* Compute sky position associated dot products
                   (for global-correlation coordinates, from Eq. (1) in PRL) */

                B1 = ( pos[0] * nvec[0] \
                       + pos[1] * nvec[1] \
                       + pos[2] * nvec[2] ); // This is \vec r \dot \vec n

                /* current unused:
                   REAL8 A2 = ( acc[0] * nvec[0]        \
                       + acc[1] * nvec[1] \
                       + acc[2] * nvec[2] ); // This is \vec a \dot \vec n
                */

                B2 = ( vel[0] * nvec[0] \
                       + vel[1] * nvec[1] \
                       + vel[2] * nvec[2] ); // This is \vec v \dot \vec n

                A1 = 1.0 + B2;

                /* Setup of the cell size (windows) in u1 */
                u1win = dFreqStack * A1;
                u1winInv = 1.0/u1win; /* Precomputing */

                /* Setup of the cell size (windows) in u2 */
                /* --- currently only u1 is needed.
                   u2win = df1dot;
                   u2winInv = 1.0/u2win; */

                /* Set starting frequency for Fstat calculation */
                thisPoint.fkdot[0] = usefulParams.spinRange_midTime.fkdot[0] - semiCohPar.extraBinsFstat * dFreqStack;

                /* Set spindown value for Fstat calculation */
                thisPoint.fkdot[1] = usefulParams.spinRange_midTime.fkdot[1] + ifdot * df1dot;

                /* Set spindown value for Fstat calculation */
                thisPoint.fkdot[2] = usefulParams.spinRange_midTime.fkdot[2] + if2dot * df2dot;

                /* Set spindown value for Fstat calculation */
                thisPoint.fkdot[3] = usefulParams.spinRange_midTime.fkdot[3] + if3dot * df3dot;

                /* Frequency at the segment's midpoint for later use */
                f1dot_event = thisPoint.fkdot[1] + thisPoint.fkdot[2] * timeDiffSeg;
                myf0 = thisPoint.fkdot[0] + thisPoint.fkdot[1] * timeDiffSeg +
                  + 0.5 * thisPoint.fkdot[2] * timeDiffSeg * timeDiffSeg;

                /* Smallest values of u1 and u2 (to be subtracted) */
                u1start = myf0 * A1 + f1dot_event * B1; /* Eq. (1a) */
                U1idx = 0;

                /* Holger: current code structure of loops (processing f1dot by f1dot) needs only U1 calculation.
                   u2start = f1dot_event + myf0 * A2 + 2.0 * f1dot_event * B2;
                   myf0max = myf0 + (Fstat_res->numFreqBins - 1) * dFreqStack;
                   u2end = f1dot_event + myf0max * A2 + 2.0 * f1dot_event * B2;
                   NumU2idx = ceil(fabs(u2start - u2end) * u2winInv);
                   U2idx = 0;
                */

                /* ----------------------------------------------------------------- */
                /************************ Compute F-Statistic ************************/
                if (doComputeFstats) { /* if first time through fine grid fdots loop */

                  tic_Fstat = GETTIME();
                  const int retn = XLALComputeFstat(&Fstat_res, usefulParams.Fstat_in_vec->data[k], &thisPoint, binsFstat1, Fstat_what);
                  if ( retn != XLAL_SUCCESS ) {
                    XLALPrintError ("%s: XLALComputeFstat() failed with errno=%d\n", __func__, xlalErrno );
                    return xlalErrno;
                  }

                  /* Loop over coarse-grid frequency bins */
                  for (ifreq = 0; ifreq < Fstat_res->numFreqBins; ifreq++) {

                    /* go to next frequency coarse-grid point */
                    freq_event = myf0 + ifreq * dFreqStack;

                    /* compute the global-correlation coordinate indices */
                    U1idx = ComputeU1idx ( freq_event, f1dot_event, A1, B1, u1start, u1winInv );

                    /* Holger: current code structure of loops (processing f1dot by f1dot) needs only U1 calculation.
                       ComputeU2idx ( freq_event, f1dot_event, A2, B2, u2start, u2winInv, &U2idx);
                    */

                    /* Check U1 index value */
                    if ( (INT4)ifreq != U1idx ) {
                      fprintf(stderr, "ERROR:  Incorrect Frequency-Index!\n ----> Seg: %03d  ifreq: %d   cg U1: %d \n",
                              k, ifreq, U1idx);
                      return(HIERARCHICALSEARCH_ECG);
                    }
                    else {
                      /* Holger: current code structure of loops (processing f1dot by f1dot) needs only U1 calculation.
                         coarsegrid.Uindex = U1idx * NumU2idx + U2idx; */
                      coarsegrid.Uindex[CG_INDEX(coarsegrid, k, ifreq)] = U1idx;
                    }

                    /* ============ Copy the *2F* value ============ */
                    coarsegrid.TwoF[CG_INDEX(coarsegrid, k, ifreq)] = Fstat_res->twoF[ifreq];
                    if ( uvar_computeBSGL ) {
                      for (UINT4 X = 0; X < coarsegrid.numDetectors; X++) {
                        INT4 detid = -1;
                        for (UINT4 Y = 0; Y < Fstat_res->numDetectors; Y++) { /* look for matching detector ID in this segment */
                          if ( strcmp( Fstat_res->detectorNames[Y], detectorIDs->data[X] ) == 0 ) {
                            detid = Y;
                          }
                        }
                        if ( detid == -1 ) { /* if no match found, detector X was not present in this segment, so use 2FX=0.0 */
                          coarsegrid.TwoFX[CG_FX_INDEX(coarsegrid, X, k, ifreq)] = 0.0;
                        } else { /* if a match was found, get the corresponding 2F value */
                          coarsegrid.TwoFX[CG_FX_INDEX(coarsegrid, X, k, ifreq)] = Fstat_res->twoFPerDet[detid][ifreq];
                        }
                      } /* for X < numDetectors */
                    } /* if ( uvar_computeBSGL ) */

                  } /* END: Loop over coarse-grid frequency bins (ifreq) */

                  /* print fstat vector if required -- mostly for debugging */
                  if ( uvar_printFstat1 )
                    {
                      LAL_CALL( PrintFstatVec ( &status, Fstat_res, fpFstat1, &thisPoint, refTimeGPS, k+1), &status);
                    }

                  // if requested, keep track of loudest candidate from each segment
                  if ( loudestTwoFPerSeg )
                    {
                      REAL4 loudestFstat_k = findLoudestTwoF ( Fstat_res );
                      loudestTwoFPerSeg[k] = fmaxf ( loudestTwoFPerSeg[k], loudestFstat_k );
                    }
                  /* --- Holger: This is not needed in U1-only case. Sort the coarse grid in Uindex --- */
                  /* qsort(coarsegrid.list, (size_t)coarsegrid.length, sizeof(CoarseGridPoint), compareCoarseGridUindex); */

                  time_Fstat += (GETTIME() - tic_Fstat);

                } // if (doComputeFstats)
                /* -------------------- END Compute F-Statistic -------------------- */

                /* -------------------- Map fine grid to coarse grid -------------------- */
                tic_SumFine = GETTIME();

                /* get the frequency of this fine-grid point at mid point of segment */
                /* OLD: ifreq_fg = 0; freq_tmp = finegrid.freqmin_fg + ifreq_fg * finegrid.dfreq_fg + f1dot_tmp * timeDiffSeg; */
                f1dot_event_fg = f1dot_fg + f2dot_fg * timeDiffSeg;
                freq_fg = finegrid.freqmin_fg + f1dot_fg * timeDiffSeg +
                  0.5 * f2dot_fg * timeDiffSeg * timeDiffSeg; /* first fine-grid frequency */

                /* compute the global-correlation coordinate indices */
                U1idx = ComputeU1idx ( freq_fg, f1dot_event_fg, A1, B1, u1start, u1winInv );

                if (U1idx < 0) {
                  fprintf(stderr,"ERROR: Stepped outside the coarse grid (%d)! \n", U1idx);
                  return(HIERARCHICALSEARCH_ECG);
                }

                if (U1idx + finegrid.freqlength - 1 >= Fstat_res->numFreqBins) {
                  fprintf(stderr,"ERROR: Stepped outside the coarse grid (%d:%d:%d:%d)! \n",
                          U1idx, finegrid.freqlength, U1idx + finegrid.freqlength - 1, Fstat_res->numFreqBins);
                  return(HIERARCHICALSEARCH_ECG);
                }

                /* coarse grid over frequency for this stack */
                REAL4 * cgrid2F = coarsegrid.TwoF + CG_INDEX(coarsegrid, k, U1idx);

                /* fine grid over frequency */
                REAL4 * fgrid2F = finegrid.sumTwoF + FG_INDEX(finegrid, 0);
#ifndef EXP_NO_NUM_COUNT
                FINEGRID_NC_T * fgridnc = finegrid.nc + FG_INDEX(finegrid, 0);
#endif

#ifdef GC_SSE2_OPT
                if ( uvar_getMaxFperSeg ) {
                  /* disables number count keeping */
                  REAL4 * fgrid2Fmax = finegrid.maxTwoFl + FG_INDEX(finegrid, 0);
                  UINT4 * fgrid2FmaxIdx = finegrid.maxTwoFlIdx + FG_INDEX(finegrid, 0);

                  gc_hotloop_2Fmax_tracking (fgrid2F, fgrid2Fmax, fgrid2FmaxIdx, cgrid2F, k, finegrid.freqlength);
                } else {
#ifndef EXP_NO_NUM_COUNT
                  gc_hotloop( fgrid2F, cgrid2F, fgridnc, TwoFthreshold, finegrid.freqlength );
#else
                  gc_hotloop_no_nc ( fgrid2F, cgrid2F, finegrid.freqlength );
#endif
		}
                if ( uvar_computeBSGL ) {
                  for (UINT4 X = 0; X < finegrid.numDetectors; X++) {
                    REAL4 * cgrid2FX = coarsegrid.TwoFX + CG_FX_INDEX(coarsegrid, X, k, U1idx);
                    REAL4 * fgrid2FX = finegrid.sumTwoFX + FG_FX_INDEX(finegrid, X, 0);

                    if ( uvar_getMaxFperSeg ) {
                      REAL4 * fgrid2FXmax = finegrid.maxTwoFXl + FG_FX_INDEX(finegrid,X, 0);
                      UINT4 * fgrid2FXmaxIdx = finegrid.maxTwoFXlIdx + FG_FX_INDEX(finegrid,X, 0);

                      gc_hotloop_2Fmax_tracking (fgrid2FX, fgrid2FXmax, fgrid2FXmaxIdx, cgrid2FX, k, finegrid.freqlength  );
                    } else {
                      gc_hotloop_no_nc( fgrid2FX, cgrid2FX, finegrid.freqlength );
                    }
                  } /* for  X  */
                }
#else // GC_SSE2_OPT
                for(UINT4 ifreq_fg=0; ifreq_fg < finegrid.freqlength; ifreq_fg++) {
                  fgrid2F[0] += cgrid2F[0];
#ifndef EXP_NO_NUM_COUNT
                  fgridnc[0] += (TwoFthreshold < cgrid2F[0]);
                  fgridnc++;
#endif // EXP_NO_NUM_COUNT
                  fgrid2F++;
                  cgrid2F++;
                }
                if ( uvar_computeBSGL ) {
                  for (UINT4 X = 0; X < finegrid.numDetectors; X++) {
                    REAL4 * cgrid2FX = coarsegrid.TwoFX + CG_FX_INDEX(coarsegrid, X, k, U1idx);
                    REAL4 * fgrid2FX = finegrid.sumTwoFX + FG_FX_INDEX(finegrid, X, 0);
                    for(UINT4 ifreq_fg=0; ifreq_fg < finegrid.freqlength; ifreq_fg++) {
                      fgrid2FX[0] += cgrid2FX[0];
                      fgrid2FX++;
                      cgrid2FX++;
                    }
                  }
                }


                if ( uvar_getMaxFperSeg ) {
                  cgrid2F = coarsegrid.TwoF + CG_INDEX(coarsegrid, k, U1idx);
                  REAL4 * fgridMax2Fl = finegrid.maxTwoFl + FG_INDEX(finegrid, 0);
                  UINT4 * fgrid2FmaxIdx = finegrid.maxTwoFlIdx + FG_INDEX(finegrid, 0);
                  int isLouder;
                  for (UINT4 ifreq_fg=0; ifreq_fg < finegrid.freqlength; ifreq_fg++) {
                    isLouder=(fgridMax2Fl[0] <= cgrid2F[0]);
                    fgridMax2Fl[0] = fmaxf ( fgridMax2Fl[0], cgrid2F[0] );
                    fgrid2FmaxIdx[0]= isLouder*k + (1-isLouder)*fgrid2FmaxIdx[0];
                    fgridMax2Fl++;
                    fgrid2FmaxIdx++;
                    cgrid2F++;
                  }
                  for (UINT4 X = 0; X < finegrid.numDetectors; X++) {
                    REAL4 * cgrid2FX = coarsegrid.TwoFX + CG_FX_INDEX(coarsegrid, X, k, U1idx);
                    REAL4 * fgridMax2FXl = finegrid.maxTwoFXl + FG_FX_INDEX(finegrid, X, 0);
                    UINT4 * fgrid2FXmaxIdx = finegrid.maxTwoFXlIdx + FG_FX_INDEX(finegrid,X, 0);

                    for(UINT4 ifreq_fg=0; ifreq_fg < finegrid.freqlength; ifreq_fg++) {
                      isLouder=(fgridMax2FXl[0] <= cgrid2FX[0]);
                      fgridMax2FXl[0] = fmaxf ( fgridMax2FXl[0], cgrid2FX[0] );
                      fgrid2FXmaxIdx[0] = isLouder*k + (1-isLouder)*fgrid2FXmaxIdx[0];
                      fgrid2FXmaxIdx++;
                      fgridMax2FXl++;
                      cgrid2FX++;
                    }
                  }
                }
#endif // GC_SSE2_OPT
                time_SumFine += ( GETTIME() - tic_SumFine );

              } /* end: ------------- MAIN LOOP over Segments --------------------*/

              /* ############################################################### */

              tic_ExtraStats = GETTIME();

              if( uvar_semiCohToplist ) {
                /* this is necessary here, because UpdateSemiCohToplists() might set
                   a checkpoint that needs some information from here */

                if(( uvar_SortToplist == SORTBY_TRIPLE_BStSGLtL || uvar_SortToplist == SORTBY_F_BSGLtL_BtSGLtL ) &&
                   finegrid.sumTwoFX && finegrid.maxTwoFXl ) {
                  LAL_CALL( UpdateSemiCohToplistsOptimTriple (&status, uvar_SortToplist, semiCohToplist, semiCohToplist2, semiCohToplist3, &finegrid, f1dot_fg, f2dot_fg, f3dot_fg, &usefulParams, NSegmentsInv, usefulParams.NSegmentsInvX, XLALUserVarWasSet(&uvar_f3dot) ), &status);
                } else {
                  LAL_CALL( UpdateSemiCohToplists (&status, semiCohToplist, semiCohToplist2, semiCohToplist3, &finegrid, f1dot_fg, f2dot_fg, f3dot_fg, &usefulParams, NSegmentsInv, usefulParams.NSegmentsInvX, XLALUserVarWasSet(&uvar_f3dot) ), &status);
                }
              } // if semiCohToplist
              time_ExtraStats += ( GETTIME() - tic_ExtraStats );

            } /* for( if1dot_fg = 0; if1dot_fg < nf1dots_fg; if1dot_fg++ ) */
            } /* for( if2dot_fg = 0; if2dot_fg < nf2dots_fg; if2dot_fg++ ) */
          } /* for( if3dot_fg = 0; if3dot_fg < nf3dots_fg; if3dot_fg++ ) */
          /* ---------- END walk through fine grid fdots --------------- */
          if3dot++;  /* Increment if3dot counter */

        } /* ########## End of loop over coarse-grid f3dot values (if3dot) ########## */
          if2dot++;  /* Increment if2dot counter */

        } /* ########## End of loop over coarse-grid f2dot values (if2dot) ########## */
        ifdot++;  /* Increment ifdot counter BEFORE SET_GCT_CHECKPOINT */

        if ( !uvar_outputTiming ) {
          SET_GCT_CHECKPOINT (uvar_fnameChkPoint, semiCohToplist, semiCohToplist2, semiCohToplist3,skyGridCounter*usefulParams.nf1dot+ifdot, TRUE);
        }

      } /* ########## End of loop over coarse-grid f1dot values (ifdot) ########## */

      /* continue forward till the end if uvar_skyPointIndex is set
         ---This probably doesn't make sense when checkpointing is turned on */
      if ( XLALUserVarWasSet(&uvar_skyPointIndex) ) {
        while(thisScan.state != STATE_FINISHED) {
          skyGridCounter++;
          XLALNextDopplerSkyPos(&dopplerpos, &thisScan);
        }
      }
      else {
        skyGridCounter++;

        XLALNextDopplerSkyPos( &dopplerpos, &thisScan );
      }

    } /* ######## End of while loop over 1st stage SKY coarse-grid points ############ */
  /*---------------------------------------------------------------------------------*/

  /* now that we have the final toplist, translate all pulsar parameters to correct reftime */
  xlalErrno = 0;
  XLALExtrapolateToplistPulsarSpins ( semiCohToplist, usefulParams.spinRange_refTime.refTime, finegrid.refTime );
  if ( semiCohToplist2 ) {	// handle (optional) second toplist
    XLALExtrapolateToplistPulsarSpins ( semiCohToplist2, usefulParams.spinRange_refTime.refTime, finegrid.refTime );
  }
  if ( semiCohToplist3 ) {	// handle (optional) second toplist
    XLALExtrapolateToplistPulsarSpins ( semiCohToplist3, usefulParams.spinRange_refTime.refTime, finegrid.refTime );
  }
  if ( xlalErrno != 0 ) {
    XLALPrintError ("%s line %d : XLALExtrapolateToplistPulsarSpins() failed with xlalErrno = %d.\n\n", __func__, __LINE__, xlalErrno );
    return(HIERARCHICALSEARCH_EXLAL);
  }
  LogPrintf( LOG_NORMAL, "Finished main analysis.\n");


  // output averaged F-stat timing model (--outputTimingDetails)
  if ( uvar_outputTimingDetails != NULL ) {
    FILE *fp;
    XLAL_CHECK ( (fp = fopen ( uvar_outputTimingDetails, "ab" )) != NULL, XLAL_ESYS, "Failed to open '%s' for writing\n", uvar_outputTimingDetails );
    for (size_t l=0; l < nStacks; l++) {
      XLAL_CHECK ( XLALAppendFstatTiming2File ( usefulParams.Fstat_in_vec->data[l], fp, (l==0) ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
    fclose ( fp );
  }

  /* Also compute F, FX (for line-robust statistics) for all candidates in final toplist */
  tic_RecalcToplist = GETTIME();
  if ( uvar_recalcToplistStats ) {

    LogPrintf( LOG_NORMAL, "Recalculating statistics for the final toplist...\n");
    RecalcStatsParams XLAL_INIT_DECL(recalcParams);
    recalcParams.listEntryTypeName = "GCTtop";
    if ( usefulParams.Fstat_in_vec_recalc != NULL ) {
      recalcParams.Fstat_in_vec		= usefulParams.Fstat_in_vec_recalc;
    } else {
      recalcParams.Fstat_in_vec		= usefulParams.Fstat_in_vec;
    }
    timing.RecalcMethodStr = XLALGetFstatInputMethodName ( recalcParams.Fstat_in_vec->data[0] );
    recalcParams.detectorIDs		= usefulParams.detectorIDs;
    recalcParams.startTstack		= usefulParams.startTstack;
    recalcParams.refTimeGPS         = refTimeGPS;
    recalcParams.BSGLsetup		    = usefulParams.BSGLsetup;
    recalcParams.loudestSegOutput	= uvar_loudestSegOutput;
    recalcParams.computeBSGLtL		= uvar_getMaxFperSeg;
    XLAL_CHECK ( XLAL_SUCCESS == XLALComputeExtraStatsForToplist ( semiCohToplist, &recalcParams ),
                 HIERARCHICALSEARCH_EXLAL, "XLALComputeExtraStatsForToplist() failed with xlalErrno = %d.\n\n", xlalErrno
                 );
    // also recalc optional 2nd toplist if present
    if ( semiCohToplist2 ) {
      XLAL_CHECK ( XLAL_SUCCESS == XLALComputeExtraStatsForToplist ( semiCohToplist2, &recalcParams ),
                   HIERARCHICALSEARCH_EXLAL, "XLALComputeExtraStatsForToplist() failed for 2nd toplist with xlalErrno = %d.\n\n", xlalErrno
                   );
    }
    // also recalc optional 3rd toplist if present
    if ( semiCohToplist3 ) {
      XLAL_CHECK ( XLAL_SUCCESS == XLALComputeExtraStatsForToplist ( semiCohToplist3, &recalcParams ),
                   HIERARCHICALSEARCH_EXLAL, "XLALComputeExtraStatsForToplist() failed for 3rd toplist with xlalErrno = %d.\n\n", xlalErrno
                   );
    }

    LogPrintf( LOG_NORMAL, "Finished recalculating toplist statistics.\n");
  } // if recalcToplist
  time_RecalcToplist = GETTIME() - tic_RecalcToplist;

  if ( uvar_outputTiming )
    {
      timing.FstatMethodStr  = XLALGetFstatInputMethodName ( usefulParams.Fstat_in_vec->data[0] );
      timing.RecalcMethodStr = (timing.RecalcMethodStr == NULL) ? "NONE" : timing.RecalcMethodStr;
      timing.Nseg = coarsegrid.nStacks;
      timing.Ndet = coarsegrid.numDetectors;
      timing.Tcoh = usefulParams.tStack;
      timing.Nsft = usefulParams.nSFTs;
      timing.Ncand = uvar_nCand1;

      timing.NFreqCo = coarsegrid.freqlength;		// includes Fstat sideband bins
      timing.Ncoh    = thisScan.numSkyGridPoints * timing.NFreqCo * usefulParams.nf1dot * usefulParams.nf2dot;

      REAL8 nf1dot_fine = usefulParams.nf1dot * nf1dots_fg;	// 'nf1dots_fg' is the number of fine-grid points *per coarse-grid point*!
      REAL8 nf2dot_fine = usefulParams.nf2dot * nf2dots_fg;	// 'nf1dots_fg' is the number of fine-grid points *per coarse-grid point*!
      timing.Ninc = thisScan.numSkyGridPoints * usefulParams.binsFstatSearch * nf1dot_fine * nf2dot_fine;	// excludes F-stat sideband bins

      timing.tau_Fstat 		= time_Fstat / ( timing.Nseg * timing.Ncoh * timing.Ndet );
      timing.tau_SumF	 	= time_SumFine  / ( timing.Nseg * timing.Ninc );
      timing.tau_Bayes	 	= time_ExtraStats / timing.Ninc;
      timing.tau_Recalc     	= time_RecalcToplist / timing.Ncand;
      time_Total = GETTIME() - tic_Start;
      timing.time_Other	= time_Total - (time_Fstat + time_SumFine + time_ExtraStats + time_RecalcToplist);	// what's left

      if ( uvar_outputTiming ) {
        XLAL_CHECK ( write_TimingInfo ( uvar_outputTiming, &timing ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
    } // if uvar_outputTiming

  // output loudest per-segment candidates
  if ( loudestTwoFPerSeg )
    {
      CHAR *fname = XLALStringDuplicate ( uvar_fnameout );
      fname = XLALStringAppend ( fname, "_loudestTwoFPerSeg");
      FILE *fp;
      if ( (fp = fopen ( fname, "wb" )) == NULL ) {
        XLALPrintError ( "Unable to open '%s' for writing\n", fname );
        return(HIERARCHICALSEARCH_EFILE);
      }
      for ( UINT4 l = 0; l < nStacks; l ++ ) {
        fprintf ( fp, "%.6" LAL_REAL4_FORMAT "\n", loudestTwoFPerSeg[l] );
      }
      fclose ( fp );
      XLALFree ( fname );
    } // if loudestTwoFPerSeg

  char t1_suffix[16];
  char t2_suffix[16];
  char t3_suffix[16];

  XLAL_INIT_MEM( t1_suffix );
  XLAL_INIT_MEM( t2_suffix );
  XLAL_INIT_MEM( t3_suffix );

  /* additional toplist(s) suffix selection */

  switch ( uvar_SortToplist ) {
    case SORTBY_DUAL_F_BSGL : {
      strcpy(t1_suffix,"");
      strcpy(t2_suffix,"-BSGL");
      strcpy(t3_suffix,"");
      break;
    };
    case SORTBY_TRIPLE_BStSGLtL:
    case SORTBY_F_BSGLtL_BtSGLtL : {
      strcpy(t1_suffix,"");
      strcpy(t2_suffix,"-BSGLtL");
      strcpy(t3_suffix,"-BtSGLtL");
      break;
    };
    default : {
      strcpy(t1_suffix,"");
      strcpy(t2_suffix,"");
      strcpy(t3_suffix,"");
    };
  }

  LogPrintf ( LOG_DEBUG, "Writing output ... ");
  {
    UINT4 newlen = strlen(uvar_fnameout) +  strlen(t1_suffix) +1;;
    CHAR *fname1;
    XLAL_CHECK ( (fname1 = XLALCalloc ( 1, newlen )) != NULL, XLAL_ENOMEM, "Failed to XLALCalloc(1, %d)\n\n", newlen );
    sprintf ( fname1, "%s%s", uvar_fnameout,t1_suffix );
    XLAL_CHECK ( write_hfs_oputput ( fname1, semiCohToplist) != -1, XLAL_EFAILED, "write_hfs_oputput('%s', toplist) failed for 1st toplist!\n", fname1 );
    XLALFree ( fname1 );
  }



  /* output optional additional toplists if any */

  if ( semiCohToplist2 ) {

    LogPrintfVerbatim ( LOG_DEBUG, "toplist2 ... ");
    UINT4 newlen = strlen(uvar_fnameout) + strlen(t2_suffix) +1;
    CHAR *fname2;
    XLAL_CHECK ( (fname2 = XLALCalloc ( 1, newlen )) != NULL, XLAL_ENOMEM, "Failed to XLALCalloc(1, %d)\n\n", newlen );
    sprintf ( fname2, "%s%s", uvar_fnameout,t2_suffix );
    XLAL_CHECK ( write_hfs_oputput ( fname2, semiCohToplist2) != -1, XLAL_EFAILED, "write_hfs_oputput('%s', toplist2) failed for 2nd toplist!\n", fname2 );
    XLALFree ( fname2 );

    if ( semiCohToplist3 )
      {
        LogPrintfVerbatim ( LOG_DEBUG, "toplist3 ... ");
        newlen = strlen(uvar_fnameout) +  strlen(t3_suffix) +1;
        CHAR *fname3;
        XLAL_CHECK ( (fname3 = XLALCalloc ( 1, newlen )) != NULL, XLAL_ENOMEM, "Failed to XLALCalloc(1, %d)\n\n", newlen );
        sprintf ( fname3, "%s%s", uvar_fnameout,t3_suffix );
        XLAL_CHECK ( write_hfs_oputput ( fname3, semiCohToplist3) != -1, XLAL_EFAILED, "write_hfs_oputput('%s', toplist3) failed for 3rd toplist!\n", fname3 );
        XLALFree ( fname3 );
      }
  }
  LogPrintfVerbatim ( LOG_DEBUG, "done.\n");

  clear_gct_checkpoint (uvar_fnameChkPoint);

  /*------------ free all remaining memory -----------*/

  if ( uvar_printCand1 ) {
    LALFree( fnameSemiCohCand );
  }

  if ( uvar_printFstat1 ) {
    fclose(fpFstat1);
    LALFree( fnameFstatVec1 );
  }

  if ( usefulParams.injectionSources ) {
    XLALDestroyPulsarParamsVector ( usefulParams.injectionSources );
  }

  XLALDestroyFstatInputVector(usefulParams.Fstat_in_vec);
  XLALDestroyFstatInputVector(usefulParams.Fstat_in_vec_recalc);

  XLALDestroyFstatResults(Fstat_res);

  XLALDestroyTimestampVector(startTstack);
  XLALDestroyTimestampVector(midTstack);
  XLALDestroyTimestampVector(endTstack);

  /* free Vel/Pos/Acc vectors and ephemeris */
  XLALDestroyREAL8VectorSequence( posStack );
  XLALDestroyREAL8VectorSequence( velStack );
  XLALDestroyREAL8VectorSequence( accStack );
  XLALDestroyEphemerisData ( edat );

  /* free dopplerscan stuff */
  LAL_CALL ( FreeDopplerSkyScan(&status, &thisScan), &status);
  if ( scanInit.skyRegionString )
    LALFree ( scanInit.skyRegionString );

  XLALDestroyStringVector ( detectorIDs );

  /* free fine grid and coarse grid */
  if (finegrid.nc) {
    ALFree(finegrid.nc);
  }
  if (finegrid.sumTwoF) {
    ALFree(finegrid.sumTwoF);
  }
  if (finegrid.sumTwoFX) {
    ALFree(finegrid.sumTwoFX);
  }
  if (finegrid.maxTwoFl) {
    ALFree(finegrid.maxTwoFl);
  }
  if (finegrid.maxTwoFlIdx) {
    ALFree(finegrid.maxTwoFlIdx);
  }

  if (finegrid.maxTwoFXl) {
    ALFree(finegrid.maxTwoFXl);
  }

  if (finegrid.maxTwoFXlIdx) {
    ALFree(finegrid.maxTwoFXlIdx);
  }

  if (coarsegrid.TwoF) {
    LALFree(coarsegrid.TwoF);
  }
  if (coarsegrid.TwoFX) {
    LALFree(coarsegrid.TwoFX);
  }
  if (coarsegrid.Uindex) {
    LALFree(coarsegrid.Uindex);
  }

  free_gctFstat_toplist ( &semiCohToplist );
  if ( semiCohToplist2 ) {
    free_gctFstat_toplist ( &semiCohToplist2 );
  }
  if ( semiCohToplist3 ) {
    free_gctFstat_toplist ( &semiCohToplist3 );
  }

  XLALDestroyBSGLSetup ( usefulParams.BSGLsetup );

  XLALDestroyUserVars();

  XLALFree ( VCSInfoString );

  XLALFree ( loudestTwoFPerSeg );

  LALCheckMemoryLeaks();

  return HIERARCHICALSEARCH_ENORM;
} /* main */







/**
 * Set up stacks, read SFTs, calculate SFT noise weights and calculate
 * detector-state
 */
void SetUpSFTs( LALStatus *status,			/**< pointer to LALStatus structure */
                UsefulStageVariables *in		/**< input params */
                )
{

  SFTCatalog *catalog = NULL;
  static SFTConstraints constraints;
  REAL8 timebase, tObs, deltaFsft;
  UINT4 k,numSFT;
  LIGOTimeGPS tStartGPS, tEndGPS, refTimeGPS, tMidGPS, midTstackGPS, startTstackGPS, endTstackGPS;
  SFTCatalogSequence catalogSeq;
  REAL8 midTseg,startTseg,endTseg;

  BOOLEAN crc_check;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* get sft catalog */
  constraints.minStartTime = &(in->minStartTimeGPS);
  constraints.maxStartTime = &(in->maxStartTimeGPS);
  XLAL_CHECK_LAL( status, ( catalog = XLALSFTdataFind( in->sftbasename, &constraints) ) != NULL, XLAL_EFUNC);

  /* check CRC sums of SFTs */
  XLAL_CHECK_LAL ( status, XLALCheckCRCSFTCatalog ( &crc_check, catalog ) == XLAL_SUCCESS, XLAL_EFUNC );
  if (!crc_check) {
    LogPrintf(LOG_CRITICAL,"SFT validity check failed\n");
    ABORT ( status, HIERARCHICALSEARCH_ESFT, HIERARCHICALSEARCH_MSGESFT );
  }

  /* set some sft parameters */
  deltaFsft = catalog->data[0].header.deltaF;
  timebase = 1.0/deltaFsft;

  /* get sft catalogs for each stack */
  if ( in->segmentList )	/* if segment list was given by user */
    {
      SFTCatalogSequence *catalogSeq_p;
      if ( (catalogSeq_p = XLALSetUpStacksFromSegmentList ( catalog, in->segmentList )) == NULL ) {
        XLALPrintError ( "%s: XLALSetUpStacksFromSegmentList() failed to set up segments from given list.\n", __func__ );
        ABORT ( status, HIERARCHICALSEARCH_ESUB, HIERARCHICALSEARCH_MSGESUB );
      }
      catalogSeq = (*catalogSeq_p);/* copy top-level struct */
      XLALFree ( catalogSeq_p );   /* free alloc'ed top-level struct after copying (contents survive in catalogSeq!) */

      /* we need to set tStack here:
       * this will be used for setting Freq,f1dot resolution on segments, therefore we use the longest segment duration
       */
      UINT4 iSeg;
      REAL8 maxT = 0;
      for ( iSeg=0; iSeg < in->segmentList->length; iSeg++)
        {
          REAL8 T = XLALGPSDiff ( &(in->segmentList->segs[iSeg].end), &(in->segmentList->segs[iSeg].start) );
          maxT = HSMAX ( maxT, T );
        }
      in->tStack = maxT;
    }
  else	/* set up nStacks segments of fixed span tStack */
    {
      TRY( SetUpStacks( status->statusPtr, &catalogSeq, in->tStack, catalog, in->nStacks), status);
    }

  /* reset number of stacks */
  UINT4 numSegments = catalogSeq.length;
  in->nStacks = numSegments;

  /* calculate start and end times and tobs from segmented catalog*/
  tStartGPS = catalogSeq.data[0].data[0].header.epoch;
  in->tStartGPS = tStartGPS;
  SFTCatalog *LastSegmentCat = &(catalogSeq.data[numSegments - 1]);
  UINT4 numSFTsInLastSeg = LastSegmentCat->length;
  tEndGPS = LastSegmentCat->data[numSFTsInLastSeg-1].header.epoch;
  XLALGPSAdd(&tEndGPS, timebase);
  tObs = XLALGPSDiff(&tEndGPS, &tStartGPS);
  in->tObs = tObs;

  /* get timestamps of start, mid and end times of each stack */
  /* set up vector containing mid times of stacks */
  in->midTstack =  XLALCreateTimestampVector ( in->nStacks );

  /* set up vector containing start times of stacks */
  in->startTstack =  XLALCreateTimestampVector ( in->nStacks );

  /* set up vector containing end times of stacks */
  in->endTstack =  XLALCreateTimestampVector ( in->nStacks );

  /* now loop over stacks and get time stamps */
  for (k = 0; k < in->nStacks; k++) {

    if ( catalogSeq.data[k].length == 0 ) {
      /* something is wrong */
      ABORT ( status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
    }

    /* start time of stack = time of first sft in stack */
    in->startTstack->data[k] = catalogSeq.data[k].data[0].header.epoch;

    /* end time of stack = time of last sft in stack */
    numSFT = catalogSeq.data[k].length;
    in->endTstack->data[k] = catalogSeq.data[k].data[numSFT - 1].header.epoch;

    /* reference time for Fstat has to be midpoint of segment */
    startTstackGPS = in->startTstack->data[k];
    endTstackGPS = in->endTstack->data[k];

    startTseg = XLALGPSGetREAL8( &startTstackGPS );
    endTseg = XLALGPSGetREAL8( &endTstackGPS );
    /*
      TRY ( LALGPStoFloat( status->statusPtr, &startTseg, &startTstackGPS ), status);
      TRY ( LALGPStoFloat( status->statusPtr, &endTseg, &endTstackGPS ), status);
    */
    midTseg = startTseg + ((endTseg - startTseg + timebase)*0.5);

    XLALGPSSetREAL8( &midTstackGPS, midTseg );
    /*
      TRY ( LALFloatToGPS( status->statusPtr, &midTstackGPS, &midTseg), status);
    */
    in->midTstack->data[k] = midTstackGPS;

  } /* loop over k */


  /* set reference time for pulsar parameters */
  /* first calculate the mid time of observation time span*/
  {
    REAL8 tStart8, tEnd8, tMid8;

    tStart8 = XLALGPSGetREAL8( &tStartGPS );
    tEnd8   = XLALGPSGetREAL8( &tEndGPS );
    tMid8 = 0.5 * (tStart8 + tEnd8);
    XLALGPSSetREAL8( &tMidGPS, tMid8 );
  }

  if ( in->refTime > 0 )  {
    REAL8 refTime = in->refTime;
    XLALGPSSetREAL8(&refTimeGPS, refTime);
  }
  else {  /* set refTime to exact midtime of the total observation-time spanned */
    refTimeGPS = tMidGPS;
  }

  /* get frequency and fdot bands at start time of sfts by extrapolating from reftime */
  in->spinRange_refTime.refTime = refTimeGPS;
  XLAL_CHECK_LAL( status, XLALExtrapolatePulsarSpinRange( &in->spinRange_startTime, &in->spinRange_refTime, XLALGPSDiff( &tStartGPS, &(&in->spinRange_refTime)->refTime ) ) == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_LAL( status, XLALExtrapolatePulsarSpinRange( &in->spinRange_endTime, &in->spinRange_refTime, XLALGPSDiff( &tEndGPS, &(&in->spinRange_refTime)->refTime ) ) == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_LAL( status, XLALExtrapolatePulsarSpinRange( &in->spinRange_midTime, &in->spinRange_refTime, XLALGPSDiff( &tMidGPS, &(&in->spinRange_refTime)->refTime ) ) == XLAL_SUCCESS, XLAL_EFUNC);

  /* set Fstat spindown resolution (coarse grid) */
  in->df1dot = HSMIN(in->df1dot, in->spinRange_midTime.fkdotBand[1]);

  /* calculate number of bins for Fstat overhead due to residual spin-down */
  in->extraBinsFstat = (UINT4)( 0.25*(in->tObs*in->df1dot + in->tObs*in->tObs*in->df2dot)/in->dFreqStack + 1e-6) + 1;

  /* calculate total number of bins for Fstat */
  if ( in->dFreqStack == 0 ) {
    in->binsFstatSearch = 1;
  } else {
    in->binsFstatSearch = (UINT4)(in->spinRange_midTime.fkdotBand[0]/in->dFreqStack + 1e-6) + 1;
  }
  /* number of coarse grid spindown values */
  if ( in->df1dot == 0 ) {
    in->nf1dot = 1;
  } else {
    in->nf1dot = (UINT4) ceil( in->spinRange_midTime.fkdotBand[1] / in->df1dot) + 1;
  }
  /* number of coarse grid 2nd spindown values */
  if ( in->df2dot == 0 ) {
    in->nf2dot = 1;
  } else {
    in->nf2dot = (UINT4) floor( in->spinRange_midTime.fkdotBand[2] / in->df2dot + NUDGE) + 1;
  }
    /* number of coarse grid 3rd spindown values */
  if ( in->df3dot == 0 ) {
    in->nf3dot = 1;
  } else {
    in->nf3dot = (UINT4) floor( in->spinRange_midTime.fkdotBand[3] / in->df3dot + NUDGE) + 1;
  }

  /* set wings of sfts to be read */
  REAL8 minCoverFreq, maxCoverFreq;
  REAL8 asiniMax = 0, PeriodMin = 0, maxEcc = 0;
  // NOTE: *must* use spin-range at *mid-time* (not reftime), which is where the GCT code sets up its
  // template bank. This is potentially 'wider' than the physically-requested template bank, and
  // can therefore also require more SFT frequency bins!
  // NOTE2: a second trap here is that the GCT code does not strictly respect the given bands, but can exceed them on the upside
  // due to 'conversative' discretization on bins, seen above. (binsFstatSearch,nf1dot,nf2dot,nf3dot)
  // therefore we use this 'effective' spinrange in order to be able to estimate the required number of SFT bins correctly:
  PulsarSpinRange spinRangeEff = in->spinRange_midTime;
  spinRangeEff.fkdotBand[0] = (in->binsFstatSearch-1) * in->dFreqStack;
  spinRangeEff.fkdotBand[1] = (in->nf1dot-1) * in->df1dot;
  spinRangeEff.fkdotBand[2] = (in->nf2dot-1) * in->df2dot;
  spinRangeEff.fkdotBand[3] = (in->nf3dot-1) * in->df3dot;
  XLALCWSignalCoveringBand ( &minCoverFreq, &maxCoverFreq, &tStartGPS, &tEndGPS, &(spinRangeEff), asiniMax, PeriodMin, maxEcc );

  REAL8 freqmin = minCoverFreq - in->extraBinsFstat * in->dFreqStack;
  REAL8 freqmax = maxCoverFreq + in->extraBinsFstat * in->dFreqStack;

  /* fill detector name vector with all detectors present in any data sements */
  if ( ( in->detectorIDs = XLALListIFOsInCatalog( catalog ) ) == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }
  const UINT4 numDetectors = in->detectorIDs->length;

  /* set up vector of Fstat input data structs */
  in->Fstat_in_vec = XLALCreateFstatInputVector( in->nStacks );
  if ( in->Fstat_in_vec == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }
  // if using different Fstat-method for main search and recalc: set-up separate Fstat input
  if ( in->recalcToplistStats && (in->Fmethod != in->FmethodRecalc) )
    {
      in->Fstat_in_vec_recalc = XLALCreateFstatInputVector( in->nStacks );
      if ( in->Fstat_in_vec == NULL ) {
        ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
      }
    }

  in->nSFTs = 0;
  for (UINT4 X = 0; X < numDetectors; X++) {
    in->NSegmentsInvX[X] = 0;
  }

  FstatOptionalArgs optionalArgs = FstatOptionalArgsDefaults;
  optionalArgs.SSBprec = in->SSBprec;
  optionalArgs.Dterms = in->Dterms;
  optionalArgs.runningMedianWindow = in->blocksRngMed;
  optionalArgs.FstatMethod = in->Fmethod;
  optionalArgs.collectTiming = in->collectFstatTiming;
  optionalArgs.injectSources = in->injectionSources;

  FstatOptionalArgs XLAL_INIT_DECL(optionalArgsRecalc);

  /* loop over segments and read sfts */
  for (k = 0; k < in->nStacks; k++) {

    /* if flag is given, assume a PSD with sqrt(S) = 1.0 */
    MultiNoiseFloor s_assumeSqrtSX;
    if ( in->assumeSqrtSX != NULL ) {
      const SFTCatalog *catalog_k = &(catalogSeq.data[k]);
      LALStringVector *detectorIDs_k = NULL;
      if ( ( detectorIDs_k = XLALListIFOsInCatalog( catalog_k ) ) == NULL ) {
        ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
      }
      if ( XLALParseMultiNoiseFloorMapped( &s_assumeSqrtSX, detectorIDs_k, in->assumeSqrtSX, in->detectorIDs ) != XLAL_SUCCESS ) {
        XLALPrintError("%s: XLALParseMultiNoiseFloorMapped() failed with errno=%d", __func__, xlalErrno);
        ABORT ( status, HIERARCHICALSEARCH_EXLAL, HIERARCHICALSEARCH_MSGEXLAL );
      }
      optionalArgs.assumeSqrtSX = &s_assumeSqrtSX;
      XLALDestroyStringVector( detectorIDs_k );
    } else {
      optionalArgs.assumeSqrtSX = NULL;
    }

    /* ----- create Fstat input data struct ----- */
    if ( k == 0 ) {
      optionalArgs.prevInput = NULL;
    } else {
      optionalArgs.prevInput = in->Fstat_in_vec->data[0];     // re-use shared workspace from first segment for all subsequent segments
    }

    in->Fstat_in_vec->data[k] = XLALCreateFstatInput ( &catalogSeq.data[k], freqmin, freqmax, in->dFreqStack, in->edat, &optionalArgs );
    if ( in->Fstat_in_vec->data[k] == NULL ) {
      XLALPrintError("%s: XLALCreateFstatInput() failed with errno=%d", __func__, xlalErrno);
      ABORT ( status, HIERARCHICALSEARCH_EXLAL, HIERARCHICALSEARCH_MSGEXLAL );
    }
    // ----- if recalc uses a different Fstat-method from main search, we'll setup its own Fstat setup struct
    if ( in->recalcToplistStats && (in->Fstat_in_vec_recalc != NULL) )
      {
        optionalArgsRecalc = optionalArgs;
        optionalArgsRecalc.FstatMethod = in->FmethodRecalc;
        optionalArgsRecalc.Dterms = in->DtermsRecalc;
        if ( k == 0 ) {
          optionalArgsRecalc.prevInput = NULL;
        } else {
          optionalArgsRecalc.prevInput = in->Fstat_in_vec_recalc->data[0];     // re-use shared workspace from first segment for all subsequent segments
        }

        in->Fstat_in_vec_recalc->data[k] = XLALCreateFstatInput ( &catalogSeq.data[k], freqmin, freqmax, in->dFreqStack, in->edat, &optionalArgsRecalc );
        if ( in->Fstat_in_vec->data[k] == NULL ) {
          XLALPrintError("%s: XLALCreateFstatInput() failed with errno=%d", __func__, xlalErrno);
          ABORT ( status, HIERARCHICALSEARCH_EXLAL, HIERARCHICALSEARCH_MSGEXLAL );
        }
      } // if Fstat_in_vec_recalc

    if ( k == 0 )
      {
        LogPrintf (LOG_NORMAL, "Search FstatMethod used: '%s'\n", XLALGetFstatInputMethodName( in->Fstat_in_vec->data[0] ) );
        LogPrintf (LOG_NORMAL, "Recalc FstatMethod used: '%s'\n", XLALGetFstatInputMethodName( in->Fstat_in_vec_recalc ? in->Fstat_in_vec_recalc->data[0] : in->Fstat_in_vec->data[0] ) );
      }

    // --------------------------------------------------
    /* get SFT detectors and timestamps */
    const MultiLALDetector *multiIFO = XLALGetFstatInputDetectors( in->Fstat_in_vec->data[k] );
    if ( multiIFO == NULL ) {
      XLALPrintError("%s: XLALGetFstatInputDetectors() failed with errno=%d", __func__, xlalErrno);
      ABORT ( status, HIERARCHICALSEARCH_EXLAL, HIERARCHICALSEARCH_MSGEXLAL );
    }
    const MultiLIGOTimeGPSVector *multiTS = XLALGetFstatInputTimestamps( in->Fstat_in_vec->data[k] );
    if ( multiTS == NULL ) {
      XLALPrintError("%s: XLALGetFstatInputTimestamps() failed with errno=%d", __func__, xlalErrno);
      ABORT ( status, HIERARCHICALSEARCH_EXLAL, HIERARCHICALSEARCH_MSGEXLAL );
    }

    /* ----- get effective inverse number of segments per detector (needed for correct averaging in single-IFO F calculation) ----- */
    for (UINT4 X = 0; X < numDetectors; X++) {
      /* for each detector, check if present in each segment, and save the number of segments where it is */
      for (UINT4 Y = 0; Y < multiTS->length; Y++) {
        if ( strcmp( multiIFO->sites[Y].frDetector.prefix, in->detectorIDs->data[X] ) == 0 )
          in->NSegmentsInvX[X] += 1;
      } /* for Y < numDetectors */
    } /* for X < numDetectors */

    /* ----- print debug info about SFTs in this stack ----- */
    LogPrintf(LOG_DETAIL, "Segment %d ", k+1);
    for ( UINT4 j = 0; j < multiIFO->length; j++) {
      LogPrintfVerbatim(LOG_DETAIL, "%s: %d  ", multiIFO->sites[j].frDetector.prefix, multiTS->data[j]->length);
    }
    LogPrintfVerbatim(LOG_DETAIL, "\n");

    /* ----- count the total and per-segment number of SFTs used ----- */
    UINT4 nSFTsInSeg = 0;
    for ( UINT4 X = 0; X < multiTS->length; ++X ) {
      nSFTsInSeg += multiTS->data[X]->length;
    }
    in->nSFTs += nSFTsInSeg;

    /* ----- if we have a segment-list: double-check number of SFTs ----- */
    if ( in->segmentList ) {
      /* check the number of SFTs we found in this segment against the nominal value, stored in the segment list field 'id' */
      UINT4 nSFTsExpected = in->segmentList->segs[k].id;
      if ( (nSFTsExpected > 0) && (nSFTsInSeg != nSFTsExpected) ) {
        XLALPrintError ("%s: Segment list seems inconsistent with data read: segment %d contains %d SFTs, should hold %d SFTs\n", __func__, k, nSFTsInSeg, nSFTsExpected );
        ABORT ( status, HIERARCHICALSEARCH_EBAD, HIERARCHICALSEARCH_MSGEBAD );
      }
    } /* if have segmentList */

  } /* loop over k */
  for (UINT4 X = 0; X < numDetectors; X++) {
    in->NSegmentsInvX[X] = 1.0 / in->NSegmentsInvX[X]; /* now it is the inverse number */
  }
  LogPrintf( LOG_NORMAL, "Number of segments: %d, total number of SFTs in segments: %d\n", in->nStacks, in->nSFTs );

  /* realloc if nStacks != in->nStacks */
  /*   if ( in->nStacks > nStacks ) { */

  /*     in->midTstack->length = nStacks; */
  /*     in->midTstack->data = (LIGOTimeGPS *)LALRealloc( in->midTstack->data, nStacks * sizeof(LIGOTimeGPS)); */

  /*     in->startTstack->length = nStacks; */
  /*     in->startTstack->data = (LIGOTimeGPS *)LALRealloc( in->startTstack->data, nStacks * sizeof(LIGOTimeGPS)); */

  /*     stackMultiSFT->length = nStacks; */
  /*     stackMultiSFT->data = (MultiSFTVector **)LALRealloc( stackMultiSFT->data, nStacks * sizeof(MultiSFTVector *)); */

  /*   }  */

  /* we don't need the original catalog anymore*/
  XLALDestroySFTCatalog(catalog );

  /* free catalog sequence */
  for (k = 0; k < in->nStacks; k++)
    {
      if ( catalogSeq.data[k].length > 0 ) {
        LALFree(catalogSeq.data[k].data);
      } /* end if */
    } /* loop over stacks */
  LALFree( catalogSeq.data);

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* SetUpSFTs */






/**
 * \brief Breaks up input sft catalog into specified number of stacks
 * Loops over elements of the catalog, assigns a bin index and
 * allocates memory to the output catalog sequence appropriately.  If
 * there are long gaps in the data, then some of the catalogs in the
 * output catalog sequence may be of zero length.
 */
void SetUpStacks(LALStatus *status,        /**< pointer to LALStatus structure */
                 SFTCatalogSequence  *out, /**< Output catalog of sfts -- one for each stack */
                 REAL8 tStack,             /**< Output duration of each stack */
                 SFTCatalog  *in,          /**< Input sft catalog to be broken up into stacks (ordered in increasing time)*/
                 UINT4 nStacksMax )        /**< User specified number of stacks */
{
  UINT4 j, stackCounter, length;
  REAL8 tStart, thisTime;
  REAL8 Tsft;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* check input parameters */
  ASSERT ( in != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( in->length > 0, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  ASSERT ( nStacksMax > 0, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  ASSERT ( in != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( tStack > 0, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( out != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

  /* set memory of output catalog sequence to maximum possible length */
  out->length = nStacksMax;
  out->data = (SFTCatalog *)LALCalloc( 1, nStacksMax * sizeof(SFTCatalog));
  if ( out->data == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }


  Tsft = 1.0 / in->data[0].header.deltaF;

  /* get first sft timestamp */
  /* tStart will be start time of a given stack.
     This initializes tStart to the first sft time stamp as this will
     be the start time of the first stack */
  tStart = XLALGPSGetREAL8(&(in->data[0].header.epoch));

  /* loop over the sfts */
  stackCounter = 0;
  for( j = 0; j < in->length; j++)
    {
      /* thisTime is current sft timestamp */
      thisTime = XLALGPSGetREAL8(&(in->data[j].header.epoch));

      /* if sft lies in stack duration then add
         this sft to the stack. Otherwise move
         on to the next stack */
      if ( (thisTime - tStart + Tsft <= tStack) )
        {
          out->data[stackCounter].length += 1;

          length = out->data[stackCounter].length;

          /* realloc to increase length of catalog */
          out->data[stackCounter].data = (SFTDescriptor *)LALRealloc( out->data[stackCounter].data, length * sizeof(SFTDescriptor));
          if ( out->data[stackCounter].data == NULL ) {
            ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
          }

          out->data[stackCounter].data[length - 1] = in->data[j];
        }
      else /* move onto the next stack */
        {
          if ( stackCounter + 1 == nStacksMax )
            break;

          stackCounter++;

          /* reset start time of stack */
          tStart = XLALGPSGetREAL8(&(in->data[j].header.epoch));

          /* realloc to increase length of catalog and copy data */
          out->data[stackCounter].length = 1;    /* first entry in new stack */
          out->data[stackCounter].data = (SFTDescriptor *)LALRealloc( out->data[stackCounter].data, sizeof(SFTDescriptor));
          if ( out->data[stackCounter].data == NULL ) {
            ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
          }

          out->data[stackCounter].data[0] = in->data[j];
        } /* if new stack */

    } /* loop over sfts */

  /* realloc catalog sequence length to actual number of stacks */
  out->length = stackCounter + 1;
  out->data = (SFTCatalog *)LALRealloc( out->data, (stackCounter+1) * sizeof(SFTCatalog) );
  if ( out->data == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* SetUpStacks() */







/** Print some sft catalog info */
void PrintCatalogInfo( LALStatus  *status,
                       const SFTCatalog *catalog,
                       FILE *fp)
{

  INT4 nSFT;
  LIGOTimeGPS start, end;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( fp != NULL, status, HIERARCHICALSEARCH_EFILE, HIERARCHICALSEARCH_MSGEFILE );
  ASSERT ( catalog != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

  nSFT = catalog->length;
  start = catalog->data[0].header.epoch;
  end = catalog->data[nSFT-1].header.epoch;

  fprintf(fp, "## Number of SFTs: %d\n", nSFT);
  fprintf(fp, "## First SFT timestamp: %d %d\n", start.gpsSeconds, start.gpsNanoSeconds);
  fprintf(fp, "## Last SFT timestamp: %d %d\n", end.gpsSeconds, end.gpsNanoSeconds);

  DETATCHSTATUSPTR (status);
  RETURN(status);

}






/** Print some stack info from sft catalog sequence*/
void PrintStackInfo( LALStatus  *status,
                     const SFTCatalogSequence *catalogSeq,
                     FILE *fp)
{

  INT4 nStacks, k;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( fp != NULL, status, HIERARCHICALSEARCH_EFILE, HIERARCHICALSEARCH_MSGEFILE );
  ASSERT ( catalogSeq != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( catalogSeq->length > 0, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  ASSERT ( catalogSeq->data != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

  nStacks = catalogSeq->length;
  fprintf(fp, "## Number of stacks: %d\n", nStacks);

  for ( k = 0; k < nStacks; k++) {
    fprintf(fp, "## Stack No. %d : \n", k+1);
    TRY ( PrintCatalogInfo( status->statusPtr, catalogSeq->data + k, fp), status);
  }

  fprintf(fp, "\n\n");

  DETATCHSTATUSPTR (status);
  RETURN(status);

}

#ifdef __GNUC__
#define likely(x)      __builtin_expect(!!(x), 1)
#define unlikely(x)    __builtin_expect(!!(x), 0)
#else
#define likely(x)      (x)
#define unlikely(x)    (x)
#endif

/**
 * Get SemiCoh candidates into toplist(s)
 * This function allows for inserting candidates into up to 3 toplists at once, which might be sorted differently!
 */
void UpdateSemiCohToplistsOptimTriple ( LALStatus *status,
                             SortBy_t toplists_sortby,
                             toplist_t *list1,
                             toplist_t *list2,  //< optional (can be NULL): insert candidate into this 2nd toplist as well
                             toplist_t *list3,  //< optional (can be NULL): insert candidate into this 3rd toplist as well
                             FineGrid *in,
                             REAL8 f1dot_fg,
                             REAL8 f2dot_fg,
                             REAL8 f3dot_fg,
                             UsefulStageVariables *usefulparams,
                             REAL4 NSegmentsInv,
                             REAL4 *NSegmentsInvX,
                             BOOLEAN have_f3dot
                             )
{

  REAL8 freq_fg;
  UINT4 ifreq_fg;
  GCTtopOutputEntry line;
  BOOLEAN delay_compute_BSGL = (toplists_sortby == SORTBY_F_BSGLtL_BtSGLtL);

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( list1 != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( in != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( usefulparams != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

 /* Optimized version for triple toplist cases: BSGL, BSGLtL,BtSGLtL or
                                                  2F, BSGLtL,BtSGLtL

    First compute all detection metrics by which any of the toplists is sorted (even as seconary
    comparison criterion).

    Next test if the candidate is loud enough to make it into at least one of the toplists.
    If so, build the candidate structure.
    If not, just try the next fine grid entry.

 */


  /* ---------- Walk through fine-grid and insert candidates into toplist--------------- */
  for( ifreq_fg = 0; ifreq_fg < in->freqlength; ifreq_fg++ ) {

    freq_fg = in->freqmin_fg + ifreq_fg * in->dfreq_fg;


    /* local placeholders for summed 2F value over segments, not averages yet */
    REAL4 sumTwoF = in->sumTwoF[ifreq_fg];
    REAL4 sumTwoFX[PULSAR_MAX_DETECTORS];

    /* compute BSGL */

    line.maxTwoFl = in->maxTwoFl[ifreq_fg];
    for (UINT4 X = 0; X < in->numDetectors; X++) {
      int fg_FX_idx= FG_FX_INDEX(*in, X, ifreq_fg);
      sumTwoFX[X] = in->sumTwoFX[fg_FX_idx]; /* here it's still the summed 2F value over segments, not the average */
      line.maxTwoFXl[X] = in->maxTwoFXl[fg_FX_idx];
    }
    xlalErrno = 0;

    /* Note: sorting a toplist by 2F takes BSGL as a secondary sorting criterion, so even if the 1st toplist is
       sorted by 2F, we still provide an overestimate of BSGL for a test if it can make it into the toplist in 
       the first place */

    if(delay_compute_BSGL) {
      line.log10BSGL = LAL_REAL4_MAX;
    } else {
      line.log10BSGL = XLALComputeBSGL ( sumTwoF, sumTwoFX, usefulparams->BSGLsetup );
    }
    if ( xlalErrno != 0 ) {
      XLALPrintError ("%s line %d : XLALComputeBSGL() failed with xlalErrno = %d.\n\n", __func__, __LINE__, xlalErrno );
      ABORT ( status, HIERARCHICALSEARCH_EXLAL, HIERARCHICALSEARCH_MSGEXLAL );
    }
    if ( unlikely(line.log10BSGL < -LAL_REAL4_MAX*0.1) ) {
      line.log10BSGL = -LAL_REAL4_MAX*0.1; /* avoid minimum value, needed for output checking in print_gctFstatline_to_str() */
    }

    line.log10BSGLtL  = XLALComputeBSGLtL ( sumTwoF, sumTwoFX, line.maxTwoFXl, usefulparams->BSGLsetup );
    if ( unlikely(xlalErrno != 0) ) {
      XLALPrintError ("%s line %d : XLALComputeBSGLtL() failed with xlalErrno = %d.\n\n", __func__, __LINE__, xlalErrno );
      ABORT ( status, HIERARCHICALSEARCH_EXLAL, HIERARCHICALSEARCH_MSGEXLAL );
    }
    if ( unlikely(line.log10BSGLtL < -LAL_REAL4_MAX*0.1) ) {
      line.log10BSGLtL = -LAL_REAL4_MAX*0.1; /* avoid minimum value, needed for output checking in print_gctFstatline_to_str() */
    }

    line.log10BtSGLtL = XLALComputeBtSGLtL ( line.maxTwoFl, sumTwoFX, line.maxTwoFXl, usefulparams->BSGLsetup );
    if ( unlikely(xlalErrno != 0) ) {
      XLALPrintError ("%s line %d : XLALComputeBSGLtL() failed with xlalErrno = %d.\n\n", __func__, __LINE__, xlalErrno );
      ABORT ( status, HIERARCHICALSEARCH_EXLAL, HIERARCHICALSEARCH_MSGEXLAL );
    }
    if ( unlikely(line.log10BtSGLtL < -LAL_REAL4_MAX*0.1) ) {
      line.log10BtSGLtL = -LAL_REAL4_MAX*0.1; /* avoid minimum value, needed for output checking in print_gctFstatline_to_str() */
    }

    line.Freq = freq_fg; /* NOTE: this is not the final output frequency! For performance reasons, it will only later get correctly extrapolated for the final toplist */
    line.Alpha = in->alpha;
    line.Delta = in->delta;
    line.F1dot = f1dot_fg;
    line.F2dot = f2dot_fg;
    line.F3dot = f3dot_fg;
    line.have_f3dot = have_f3dot;
    line.nc = in->nc[ifreq_fg];

    /* take F-stat averages over segments */
    line.avTwoF = sumTwoF*NSegmentsInv; /* average multi-2F by full number of segments */

    /* now test if this candidate makes it into any of the toplists, if not we don't need
       to do any more copying and computing of data from finegrid to toplist entry structure */


    int isIncludedToplists = TEST_FSTAT_TOPLIST_INCLUSION( list1, &line) ||
                             TEST_FSTAT_TOPLIST_INCLUSION( list2, &line) ||
                             TEST_FSTAT_TOPLIST_INCLUSION( list3, &line);


    if(likely(! isIncludedToplists)) continue;


    /* if we delayed the computation of BSGL for a 2F sorted first toplist, now is the time to actually compute it */
    if(delay_compute_BSGL) {
      line.log10BSGL = XLALComputeBSGL ( sumTwoF, sumTwoFX, usefulparams->BSGLsetup );
    }

    line.numDetectors = in->numDetectors;
    line.avTwoFrecalc = -1.0; /* initialise this to -1.0, so that it only gets written out by print_gctFstatline_to_str if later overwritten in recalcToplistStats step */
    line.log10BSGLrecalc = -LAL_REAL4_MAX; /* for now, block field with minimal value, needed for output checking in print_gctFstatline_to_str() */
    line.log10BSGLtLrecalc = -LAL_REAL4_MAX; /* for now, block field with minimal value, needed for output checking in print_gctFstatline_to_str() */
    line.loudestSeg = -1;
    line.twoFloudestSeg = -1.0;
    for (UINT4 X = 0; X < PULSAR_MAX_DETECTORS; X++) { /* initialise single-IFO F-stat arrays to zero */
      line.twoFXloudestSeg[X] = -1.0;
      line.avTwoFXrecalc[X] = 0.0;
    }


    line.maxTwoFlSeg = in->maxTwoFlIdx[ifreq_fg];
    for (UINT4 X = 0; X < in->numDetectors; X++) {
      line.avTwoFX[X] = sumTwoFX[X]*NSegmentsInvX[X]; /* average single-2F by per-IFO number of segments */
      line.maxTwoFXlSeg[X] = in->maxTwoFXlIdx[FG_FX_INDEX(*in, X, ifreq_fg)];
    }

    insert_into_gctFstat_toplist( list1, &line);
    insert_into_gctFstat_toplist( list2, &line);
    insert_into_gctFstat_toplist( list3, &line);

  } // for ifreq_fg

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* UpdateSemiCohToplistsOptimTriple() */




/**
 * Get SemiCoh candidates into toplist(s)
 * This function allows for inserting candidates into up to 3 toplists at once, which might be sorted differently!
 */
void UpdateSemiCohToplists ( LALStatus *status,
                             toplist_t *list1,
                             toplist_t *list2,	//< optional (can be NULL): insert candidate into this 2nd toplist as well
                             toplist_t *list3,	//< optional (can be NULL): insert candidate into this 3rd toplist as well
                             FineGrid *in,
                             REAL8 f1dot_fg,
                             REAL8 f2dot_fg,
                             REAL8 f3dot_fg,
                             UsefulStageVariables *usefulparams,
                             REAL4 NSegmentsInv,
                             REAL4 *NSegmentsInvX,
                             BOOLEAN have_f3dot
                             )
{

  REAL8 freq_fg;
  UINT4 ifreq_fg;
  GCTtopOutputEntry line;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( list1 != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( in != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( usefulparams != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

  /* ---------- Walk through fine-grid and insert candidates into toplist--------------- */
  for( ifreq_fg = 0; ifreq_fg < in->freqlength; ifreq_fg++ ) {

    freq_fg = in->freqmin_fg + ifreq_fg * in->dfreq_fg;

    line.Freq = freq_fg; /* NOTE: this is not the final output frequency! For performance reasons, it will only later get correctly extrapolated for the final toplist */
    line.Alpha = in->alpha;
    line.Delta = in->delta;
    line.F1dot = f1dot_fg;
    line.F2dot = f2dot_fg;
    line.F3dot = f3dot_fg;
    line.nc = in->nc[ifreq_fg];
    line.avTwoF = 0.0; /* will be set to average over segments later */
    line.maxTwoFl = -1.0; /* initialise this to -1.0, so that it only gets written out by print_gctFstatline_to_str if actually computed */
    line.maxTwoFlSeg = -1;
    line.log10BSGL    = -LAL_REAL4_MAX; /* for now, block field with minimal value, needed for output checking in print_gctFstatline_to_str() */
    line.log10BSGLtL  = -LAL_REAL4_MAX;
    line.log10BtSGLtL = -LAL_REAL4_MAX;

    line.numDetectors = in->numDetectors;
    for (UINT4 X = 0; X < PULSAR_MAX_DETECTORS; X++) { /* initialise single-IFO F-stat arrays to zero */
      line.avTwoFX[X] = 0.0;
      line.maxTwoFXl[X] = 0.0;
      line.maxTwoFXlSeg[X] = -1;
      line.avTwoFXrecalc[X] = 0.0;
    }
    line.avTwoFrecalc = -1.0; /* initialise this to -1.0, so that it only gets written out by print_gctFstatline_to_str if later overwritten in recalcToplistStats step */
    line.log10BSGLrecalc = -LAL_REAL4_MAX; /* for now, block field with minimal value, needed for output checking in print_gctFstatline_to_str() */
    line.log10BSGLtLrecalc = -LAL_REAL4_MAX; /* for now, block field with minimal value, needed for output checking in print_gctFstatline_to_str() */
    line.have_f3dot = have_f3dot;
    line.loudestSeg = -1;
    line.twoFloudestSeg = -1.0;
    for (UINT4 X = 0; X < PULSAR_MAX_DETECTORS; X++) {
      line.twoFXloudestSeg[X] = -1.0;
    }

    /* local placeholders for summed 2F value over segments, not averages yet */
    REAL4 sumTwoF = in->sumTwoF[ifreq_fg];
    REAL4 sumTwoFX[PULSAR_MAX_DETECTORS];
    if ( in->sumTwoFX ) { /* if we already have FX values from the main loop, insert these, and calculate BSGL here */
      for (UINT4 X = 0; X < in->numDetectors; X++) {
        sumTwoFX[X] = in->sumTwoFX[FG_FX_INDEX(*in, X, ifreq_fg)]; /* here it's still the summed 2F value over segments, not the average */
      }
      xlalErrno = 0;

      line.log10BSGL = XLALComputeBSGL ( sumTwoF, sumTwoFX, usefulparams->BSGLsetup );
      if ( xlalErrno != 0 ) {
        XLALPrintError ("%s line %d : XLALComputeBSGL() failed with xlalErrno = %d.\n\n", __func__, __LINE__, xlalErrno );
        ABORT ( status, HIERARCHICALSEARCH_EXLAL, HIERARCHICALSEARCH_MSGEXLAL );
      }
      if ( line.log10BSGL < -LAL_REAL4_MAX*0.1 ) {
        line.log10BSGL = -LAL_REAL4_MAX*0.1; /* avoid minimum value, needed for output checking in print_gctFstatline_to_str() */
      }
    }
    else {
      line.log10BSGL = -LAL_REAL4_MAX; /* in non-BSGL case, block field with minimal value, needed for output checking in print_gctFstatline_to_str() */
    }

    /* take F-stat averages over segments */
    line.avTwoF = sumTwoF*NSegmentsInv; /* average multi-2F by full number of segments */
    if ( in->sumTwoFX ) {
      for (UINT4 X = 0; X < in->numDetectors; X++) {
        line.avTwoFX[X] = sumTwoFX[X]*NSegmentsInvX[X]; /* average single-2F by per-IFO number of segments */
      }
    }

    if ( in->maxTwoFXl ) { /* if we already have max-per-segment values from the main loop, insert these too */
      line.maxTwoFl = in->maxTwoFl[ifreq_fg];
      line.maxTwoFlSeg = in->maxTwoFlIdx[ifreq_fg];
      for (UINT4 X = 0; X < in->numDetectors; X++) {
       line.maxTwoFXl[X] = in->maxTwoFXl[FG_FX_INDEX(*in, X, ifreq_fg)];
       line.maxTwoFXlSeg[X] = in->maxTwoFXlIdx[FG_FX_INDEX(*in, X, ifreq_fg)];
      }

      line.log10BSGLtL  = XLALComputeBSGLtL ( sumTwoF, sumTwoFX, line.maxTwoFXl, usefulparams->BSGLsetup );
      if ( xlalErrno != 0 ) {
        XLALPrintError ("%s line %d : XLALComputeBSGLtL() failed with xlalErrno = %d.\n\n", __func__, __LINE__, xlalErrno );
        ABORT ( status, HIERARCHICALSEARCH_EXLAL, HIERARCHICALSEARCH_MSGEXLAL );
      }
      if ( line.log10BSGLtL < -LAL_REAL4_MAX*0.1 ) {
        line.log10BSGLtL = -LAL_REAL4_MAX*0.1; /* avoid minimum value, needed for output checking in print_gctFstatline_to_str() */
      }

      line.log10BtSGLtL = XLALComputeBtSGLtL ( line.maxTwoFl, sumTwoFX, line.maxTwoFXl, usefulparams->BSGLsetup );
      if ( xlalErrno != 0 ) {
        XLALPrintError ("%s line %d : XLALComputeBtSGLtL() failed with xlalErrno = %d.\n\n", __func__, __LINE__, xlalErrno );
        ABORT ( status, HIERARCHICALSEARCH_EXLAL, HIERARCHICALSEARCH_MSGEXLAL );
      }
      if ( line.log10BtSGLtL < -LAL_REAL4_MAX*0.1 ) {
        line.log10BtSGLtL = -LAL_REAL4_MAX*0.1; /* avoid minimum value, needed for output checking in print_gctFstatline_to_str() */
      }

    }

    insert_into_gctFstat_toplist( list1, &line);
    if ( list2 ){	// also insert candidate into (optional) second toplist
      insert_into_gctFstat_toplist( list2, &line);
    }
    if ( list3 ){	// also insert candidate into (optional) 3rd toplist
      insert_into_gctFstat_toplist( list3, &line);
    }

  } // for ifreq_fg

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* UpdateSemiCohToplists() */


// simply return maximal (multi-IFO) Fstat-value over frequency-bins in in->twoF
static inline REAL4
findLoudestTwoF ( const FstatResults *in )
{
  REAL4 maxTwoF = 0;
  for ( UINT4 k = 0; k < in->numFreqBins; k ++ )
    {
      maxTwoF = fmaxf ( maxTwoF, in->twoF[k] );
    } // k < numFreqBins

  return maxTwoF;
} // findLoudestTwoF()


/** Print Fstat vectors */
void PrintFstatVec (LALStatus *status,
                    FstatResults         *in,
                    FILE                 *fp,
                    PulsarDopplerParams  *thisPoint,
                    LIGOTimeGPS          refTime,
                    INT4                 stackIndex)
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  fprintf(fp, "%% Fstat values from stack %d (reftime -- %d %d)\n", stackIndex, refTime.gpsSeconds, refTime.gpsNanoSeconds);

  REAL8 alpha = thisPoint->Alpha;
  REAL8 delta = thisPoint->Delta;

  PulsarSpins fkdot;
  memcpy ( fkdot, thisPoint->fkdot, sizeof(fkdot) );

  UINT4 length = in->numFreqBins;
  REAL8 deltaF = in->dFreq;

  REAL8 f0 = fkdot[0];
  for (UINT4 k=0; k<length; k++)
    {
      fkdot[0] = f0 + k*deltaF;

      /* propagate fkdot back to reference-time  */
      XLAL_CHECK_LAL ( status, XLALExtrapolatePulsarSpins( fkdot, fkdot, XLALGPSDiff( &refTime, &thisPoint->refTime  ) ) == XLAL_SUCCESS, XLAL_EFUNC );

      fprintf(fp, "%d %.13g %.12g %.12g %.13g %.13g %.6g\n",
              stackIndex, fkdot[0], alpha, delta, fkdot[1], fkdot[2], in->twoF[k]);
    }

  fprintf(fp, "\n");

  DETATCHSTATUSPTR (status);
  RETURN(status);

}















/**
 * Calculate Earth orbital position, velocity and acceleration
 * at midpoint of each segment
 */
void GetSegsPosVelAccEarthOrb( LALStatus *status,
                               REAL8VectorSequence **posSeg,
                               REAL8VectorSequence **velSeg,
                               REAL8VectorSequence **accSeg,
                               UsefulStageVariables *usefulparams)
{

  UINT4 k, nStacks;
  LIGOTimeGPSVector *tsMid;
  vect3Dlist_t *pvaUR = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( usefulparams != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( usefulparams->nStacks > 0, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  ASSERT ( usefulparams->midTstack != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( posSeg != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( velSeg != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( accSeg != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

  /* local copies */
  nStacks = usefulparams->nStacks;
  tsMid =  usefulparams->midTstack;

  /* get pos,vel,acc at midpoint of each segment*/
  for (k = 0; k < nStacks; k++)
    {
      /* initialize velocities and positions */
      posSeg[0]->data[3*k]   = 0.0;
      posSeg[0]->data[3*k+1] = 0.0;
      posSeg[0]->data[3*k+2] = 0.0;

      velSeg[0]->data[3*k]   = 0.0;
      velSeg[0]->data[3*k+1] = 0.0;
      velSeg[0]->data[3*k+2] = 0.0;

      accSeg[0]->data[3*k]   = 0.0;
      accSeg[0]->data[3*k+1] = 0.0;
      accSeg[0]->data[3*k+2] = 0.0;

      /* get Earth's orbital pos vel acc  */
      if ( (pvaUR = XLALComputeOrbitalDerivatives(3, &tsMid->data[k], usefulparams->edat) ) == NULL ) {
        LogPrintf(LOG_CRITICAL,"GetSegsPosVelAccEarthOrb(): XLALComputeOrbitalDerivatives() failed.\n");
        ABORT ( status, HIERARCHICALSEARCH_ESFT, HIERARCHICALSEARCH_MSGESFT );
      }

      posSeg[0]->data[3*k]   = pvaUR->data[0][0];
      posSeg[0]->data[3*k+1] = pvaUR->data[0][1];
      posSeg[0]->data[3*k+2] = pvaUR->data[0][2];

      velSeg[0]->data[3*k]   = pvaUR->data[1][0];
      velSeg[0]->data[3*k+1] = pvaUR->data[1][1];
      velSeg[0]->data[3*k+2] = pvaUR->data[1][2];

      accSeg[0]->data[3*k]   = pvaUR->data[2][0];
      accSeg[0]->data[3*k+1] = pvaUR->data[2][1];
      accSeg[0]->data[3*k+2] = pvaUR->data[2][2];

      XLALDestroyVect3Dlist ( pvaUR );

    } /* loop over segment -- end pos vel acc of Earth's orbital motion */


  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* GetSegsPosVelAccEarthOrb() */







/** Calculate the U1 index for a given point in parameter space */
static inline INT4 ComputeU1idx( REAL8 freq_event,
                                 REAL8 f1dot_event,
                                 REAL8 A1,
                                 REAL8 B1,
                                 REAL8 U1start,
                                 REAL8 U1winInv)
{
  /* compute the index of global-correlation coordinate U1, Eq. (1) */
  return (((freq_event * A1 + f1dot_event * B1) - U1start) * U1winInv) + 0.5;

} /* ComputeU1idx */





/** Calculate the U2 index for a given point in parameter space */
void ComputeU2idx( REAL8 freq_event,
                   REAL8 f1dot_event,
                   REAL8 A2,
                   REAL8 B2,
                   REAL8 U2start,
                   REAL8 U2winInv,
                   INT4 *U2idx)
{

  /* compute the index of global-correlation coordinate U2 */
  *U2idx = (INT4) ((((f1dot_event + freq_event * A2 + 2.0 * f1dot_event * B2) - U2start) * U2winInv) + 0.5);

  return;

} /* ComputeU2idx */

/**
 * Set up 'segmented' SFT-catalogs for given list of segments and a total SFT-catalog.
 *
 * Note: this function does not allow 'empty' segments to be returned, i.e. if there is any
 * segment that would contain no SFTs from the given SFT-catalog, an error is returned.
 * These segment-lists are 'precomputed' and therefore one can assume that empty segments
 * are not intended.
 * However, the function will not complain if some SFTs from the catalog are 'unused', i.e.
 * they didn't fit into any of the given segments.
 *
 * \note the input segment list must be sorted, otherwise an error is returned
 *
 */
SFTCatalogSequence *
XLALSetUpStacksFromSegmentList ( const SFTCatalog *catalog,	/**< complete list of SFTs read in */
                                 const LALSegList *segList	/**< pre-computed list of segments to split SFTs into */
                                 )
{
  SFTCatalogSequence *stacks;	/* output: segmented SFT-catalogs */

  /* check input consistency */
  if ( !catalog || !segList ) {
    XLALPrintError ("%s: invalid NULL input\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }
  /* check that segment list is sorted */
  if ( ! segList->sorted ) {
    XLALPrintError ("%s: input segment list must be sorted! -> Use XLALSegListSort()\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EDOM );
  }

  UINT4 numSegments = segList->length;
  UINT4 numSFTs = catalog->length;

  /* set memory of output catalog sequence to maximum possible length */
  if ( (stacks = XLALCalloc ( 1, sizeof(*stacks) ) ) == NULL ) {
    XLALPrintError ("%s: XLALCalloc(%zu) failed.\n", __func__, sizeof(*stacks) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  stacks->length = numSegments;
  if ( (stacks->data = XLALCalloc( stacks->length, sizeof(*stacks->data) )) == NULL ) {
    XLALPrintError ("%s: failed to allocate segmented SFT-catalog\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  /* Step through segment list:
   * for every segment:
   *  - find earliest and last SFT *starting* within given segment
   *     this ensures we use all SFTs and dont lose some that fall on segment boundaries
   *  - copy this range of SFT-headers into the segmented SFT-catalog 'stacks'
   */
  UINT4 iSeg;
  INT4 iSFT0 = 0, iSFT1 = 0;	/* indices of earliest and last SFT fitting into segment iSeg */
  for ( iSeg = 0; iSeg < numSegments; iSeg ++ )
    {
      LALSeg *thisSeg = &(segList->segs[iSeg]);

      /* ----- find earliest SFT fitting into this segment */
      iSFT0 = iSFT1;	/* start from previous segment's last SFT */
      while ( 1 )
        {
          int cmp = XLALCWGPSinRange( catalog->data[iSFT0].header.epoch, &thisSeg->start, &thisSeg->end );

          if ( cmp < 0 )	/* iSFT0 lies *before* current segment => advance */
            iSFT0 ++;
          if ( cmp == 0 )	/* iSFT0 lies *inside* current segment ==> stop */
            break;

          /* if no more SFTs or iSFT0 lies *past* current segment => ERROR: empty segment! */
          if ( cmp > 0 || iSFT0 == (INT4)numSFTs )
            {
              XLALPrintError ("%s: Empty segment! No SFTs fit into segment iSeg=%d\n", __func__, iSeg );
              XLAL_ERROR_NULL ( XLAL_EDOM );
            }
        } /* while true */

      /* ----- find last SFT still starting within segment */
      iSFT1 = iSFT0;
      while ( 1 )
        {
          int cmp = XLALCWGPSinRange( catalog->data[iSFT1].header.epoch, &thisSeg->start, &thisSeg->end );

          if ( cmp < 0 ) {      /* start of iSFT1 lies *before* current segment ==> something is screwed up! */
            XLALPrintError ("%s: start of current SFT %d lies before current segment %d ==> code seems inconsistent!\n", __func__, iSFT1, iSeg );
            XLAL_ERROR_NULL ( XLAL_EFAILED );
          }
          if ( cmp == 0 )	/* start of iSFT1 lies *inside* current segment ==> advance */
            iSFT1 ++;

          if ( cmp > 0 || iSFT1 == (INT4)numSFTs ) {	/* last SFT reached or start of iSFT1 lies *past* current segment => step back once and stop */
            iSFT1 --;
            break;
          }

        } /* while true */

      INT4 numSFTsInSeg = iSFT1 - iSFT0 + 1;

      /* ----- allocate and copy this range of SFTs into the segmented catalog */
      stacks->data[iSeg].length = (UINT4)numSFTsInSeg;
      UINT4 size = sizeof(*stacks->data[iSeg].data);
      if ( (stacks->data[iSeg].data = XLALCalloc ( numSFTsInSeg, size)) == NULL ) {
        XLALPrintError ("%s: failed to XLALCalloc(%d, %d)\n", __func__, numSFTsInSeg, size );
        XLAL_ERROR_NULL ( XLAL_ENOMEM );
      }

      INT4 iSFT;
      for ( iSFT = iSFT0; iSFT <= iSFT1; iSFT++ )
        stacks->data[iSeg].data[iSFT - iSFT0] = catalog->data[iSFT];

    } /* for iSeg < numSegments */

  return stacks;

} /* XLALSetUpStacksFromSegmentList() */




/**
 * XLAL function to extrapolate the pulsar spin parameters of all toplist candidates
 * from reftime of the input toplist ('inRefTime') to a user-specified output reftime 'outRefTime'
 */
int XLALExtrapolateToplistPulsarSpins ( toplist_t *list,              /**< [out/in] toplist with GCTtopOutputEntry items, 'Freq,F1dot,F2dot,F3dot' fields will be overwritten  */
					const LIGOTimeGPS outRefTime, /**< reference time as requested for the final candidate output */
					const LIGOTimeGPS inRefTime   /**< reference time of the input toplist */
				        )
{

  /* check input parameters */
  if ( !list )
    XLAL_ERROR ( XLAL_EFAULT, "\nNULL pointer given instead of toplist.\n" );

  /* convert LIGOTimeGPS into real number difference for XLALExtrapolatePulsarSpins */
  REAL8 deltaTau = XLALGPSDiff( &outRefTime, &inRefTime );

  /* check if translation to reference time of toplist is necessary */
  if  ( deltaTau == 0 )
    return XLAL_SUCCESS; /* can skip this step if reftimes are equal */

  PulsarSpins XLAL_INIT_DECL(fkdot);

  UINT4 numElements = list->elems;
  for (UINT4 j = 0; j < numElements; j++ ) /* loop over toplist */
    {
      /* get fkdot of each candidate */
      GCTtopOutputEntry *elem = toplist_elem ( list, j );
      fkdot[0] = elem->Freq;
      fkdot[1] = elem->F1dot;
      fkdot[2] = elem->F2dot;
	  fkdot[3] = elem->F3dot;
      /* propagate fkdot to reference-time  */
      if ( XLALExtrapolatePulsarSpins( fkdot, fkdot, deltaTau ) != XLAL_SUCCESS )
        {
          XLALPrintError ("\n%s, line %d : XLALExtrapolatePulsarSpins() failed.\n\n", __func__, __LINE__);
          XLAL_ERROR ( XLAL_EFUNC );
        }
      /* write back propagated frequency to toplist */
      elem->Freq  = fkdot[0];
      elem->F1dot = fkdot[1];
      elem->F2dot = fkdot[2];
	  elem->F3dot = fkdot[3];
    }

  return XLAL_SUCCESS;

} /* XLALExtrapolateToplistPulsarSpins() */


/**
 * Function to append one timing-info line to output file.
 *
 */
static int
write_TimingInfo ( const CHAR *fname, const timingInfo_t *ti )
{
  /* input sanity */
  if ( !fname || !ti ) {
    XLALPrintError ("%s: invalid NULL input 'fp' | 'ti'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  FILE *fp;
  if ( (fp = fopen ( fname,"rb" )) == NULL )
    {
      XLAL_CHECK ( (fp = fopen( fname, "wb" )) != NULL, XLAL_ESYS, "Failed to open new timing-file '%s' for writing\n", fname );
      fprintf ( fp, "%%%%--------------------------------------------------------------------------------\n");
      fprintf ( fp, "%%%% GCT Timing model:\n");
      fprintf ( fp, "%%%% runtime = Nseg * Ndet * Ncoh * tau_Fstat + Nseg * Ninc * tau_SumF + Ninc * tau_Bayes + Ncan * tau_Recalc + time_Other\n" );
      fprintf ( fp, "%%%%--------------------------------------------------------------------------------\n");
      fprintf ( fp, "%%%%\n");
      fprintf ( fp, "%2s%10s %10s %10s %10s %10s | %6s %6s %10s %6s %10s %10s %10s %10s %%%15s %15s\n",
                "%%", "tau_Fstat", "tau_SumF", "tau_Bayes", "tau_Recalc", "time_Other",
                "Nseg", "Ndet", "Tcoh[s]", "Nsft", "NFreqCo", "Ncoh", "Ninc", "Ncand", "FstatMethod", "RecalcMethod" );
    }
  else
    {
      fclose(fp);
      XLAL_CHECK ( (fp = fopen( fname, "ab" )) != NULL, XLAL_ESYS, "Failed to open existing timing-file '%s' for appending\n", fname );
    }

  fprintf ( fp, "%12.1e %10.1e %10.1e %10.1e %10.1e   %6d %6d %10d %6d %10d %10.1e %10.1e %10d %%%15s %15s\n",
            ti->tau_Fstat, ti->tau_SumF, ti->tau_Bayes, ti->tau_Recalc, ti->time_Other,
            ti->Nseg, ti->Ndet, ti->Tcoh, ti->Nsft, ti->NFreqCo, ti->Ncoh, ti->Ninc, ti->Ncand, ti->FstatMethodStr, ti->RecalcMethodStr  );

  fclose ( fp );
  return XLAL_SUCCESS;

} // write_TimingInfo()
