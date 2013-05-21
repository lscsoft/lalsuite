/*
 *  Copyright (C) 2011 Karl Wette.
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
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 *
 */

/*********************************************************************************/
/** \author Holger Pletsch
 * \file
 * \ingroup pulsarApps
 * \brief Hierarchical semicoherent CW search code based on F-Statistic,
 *  exploiting global-correlation coordinates (Phys.Rev.Lett. 103, 181102, 2009)
 *
 *********************************************************************************/

/* ---------- Includes -------------------- */
#include <lal/Segments.h>

#include "HierarchSearchGCT.h"
#include <lal/TransientCW_utils.h> /* for XLALFastNegExp */

#include "LineVeto.h"

#ifdef GC_SSE2_OPT
#include <gc_hotloop_sse2.h>
#else
#define ALRealloc LALRealloc
#define ALFree LALFree
#endif

/* ---------- Defines -------------------- */
/* #define OUTPUT_TIMING 1 */
/* #define DIAGNOSISMODE 1 */
#define NUDGE	10*LAL_REAL8_EPS

#define TRUE (1==1)
#define FALSE (1==0)

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* Hooks for Einstein@Home / BOINC
   These are defined to do nothing special in the standalone case
   and will be set in boinc_extras.h if EAH_BOINC is set
*/
#ifdef EAH_BOINC
#include "hs_boinc_extras.h"
#define COMPUTEFSTATFREQBAND_RS ComputeFStatFreqBand_RS
#else
#define GET_GCT_CHECKPOINT read_gct_checkpoint // (cptname, semiCohToplist, NULL, &count)
#define SET_GCT_CHECKPOINT write_gct_checkpoint
#define SHOW_PROGRESS(rac,dec,skyGridCounter,tpl_total,freq,fband)
#define MAIN  main
#ifdef HS_OPTIMIZATION
extern void
LocalComputeFStatFreqBand ( LALStatus *status,
                            REAL4FrequencySeries *FstatVector,
                            const PulsarDopplerParams *doppler,
                            const MultiSFTVector *multiSFTs,
                            const MultiNoiseWeights *multiWeights,
                            const MultiDetectorStateSeries *multiDetStates,
                            const ComputeFParams *params);
#define COMPUTEFSTATFREQBAND LocalComputeFStatFreqBand
extern void
LocalComputeFStat ( LALStatus *status,
		    Fcomponents *Fstat,
		    const PulsarDopplerParams *doppler,
		    const MultiSFTVector *multiSFTs,
		    const MultiNoiseWeights *multiWeights,
		    const MultiDetectorStateSeries *multiDetStates,
		    const ComputeFParams *params,
		    ComputeFBuffer *cfBuffer);
#define COMPUTEFSTAT LocalComputeFStat
#else
#define COMPUTEFSTATFREQBAND ComputeFStatFreqBand
#define COMPUTEFSTAT ComputeFStat
#endif
#define COMPUTEFSTATFREQBAND_RS ComputeFStatFreqBand_RS
char**global_argv;
int global_argc;
#endif /* EAH_BOINC */

#define EARTHEPHEMERIS  "earth00-19-DE405.dat"
#define SUNEPHEMERIS    "sun00-19-DE405.dat"
#define BLOCKSRNGMED    101     /**< Default running median window size */
#define FSTART          100.0	/**< Default Start search frequency */
#define FBAND           0.0  /**< Default search band */
#define FDOT            0.0       /**< Default value of first spindown */
#define DFDOT           0.0       /**< Default range of first spindown parameter */
#define F2DOT           0.0       /**< Default value of second spindown */
#define DF2DOT          0.0       /**< Default range of second spindown parameter */
#define SKYREGION       "allsky" /**< default sky region to search over -- just a single point*/
#define DTERMS          8    /**< Default number of dirichlet kernel terms for calculating Fstat */

/**< Default number of dirichlet kernel terms for calculating Fstat */
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
#define INIT_MEM(x) memset(&(x), 0, sizeof((x)))


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
  EphemerisData *edat;             /**< ephemeris data for LALBarycenter */
  LIGOTimeGPSVector *midTstack;    /**< timestamps vector for mid time of each stack */
  LIGOTimeGPSVector *startTstack;  /**< timestamps vector for start time of each stack */
  LIGOTimeGPSVector *endTstack;    /**< timestamps vector for end time of each stack */
  LIGOTimeGPS minStartTimeGPS;     /**< all sft data must be after this time */
  LIGOTimeGPS maxEndTimeGPS;       /**< all sft data must be before this time */
  UINT4 blocksRngMed;              /**< blocksize for running median noise floor estimation */
  UINT4 Dterms;                    /**< size of Dirichlet kernel for Fstat calculation */
  BOOLEAN SignalOnly;              /**< FALSE: estimate noise-floor from data, TRUE: assume Sh=1 */
  REAL8 dopplerMax;                /**< extra sft wings for doppler motion */
  /* parameters describing the coherent data-segments */
  REAL8 tStack;                    /**< duration of stacks */
  UINT4 nStacks;                   /**< number of stacks */
  LALSegList *segmentList;         /**< parsed segment list read from user-specified input file --segmentList */
  BOOLEAN LVuseAllTerms;           /**< which terms to use in LineVeto computation - FALSE: only leading term, TRUE: all terms */
  REAL8 LVlogRhoTerm;              /**< For LineVeto statistic: extra term coming from prior normalization: log(rho_max_line^4/70) */
  REAL8Vector *LVloglX;            /**< For LineVeto statistic: vector of logs of line prior ratios lX per detector */
  REAL8 dFreqStack;                /**< frequency resolution of Fstat calculation */
  REAL8 df1dot;                    /**< coarse grid resolution in spindown */
  REAL8 df2dot;                    /**< coarse grid resolution in 2nd spindown */
  UINT4 extraBinsFstat;            /**< Extra bins required for Fstat calculation */
} UsefulStageVariables;


/* ------------------------ Functions -------------------------------- */
void SetUpSFTs( LALStatus *status, MultiSFTVectorSequence *stackMultiSFT,
                MultiNoiseWeightsSequence *stackMultiNoiseWeights,
                MultiDetectorStateSeriesSequence *stackMultiDetStates, UsefulStageVariables *in, BOOLEAN useWholeSFTs, REAL8 mismatch1);
void PrintFstatVec( LALStatus *status, REAL4FrequencySeries *in, FILE *fp, PulsarDopplerParams *thisPoint,
                    LIGOTimeGPS refTime, INT4 stackIndex);
void PrintCatalogInfo( LALStatus *status, const SFTCatalog *catalog, FILE *fp );
void PrintStackInfo( LALStatus *status, const SFTCatalogSequence *catalogSeq, FILE *fp );
void UpdateSemiCohToplists ( LALStatus *status, toplist_t *list1, toplist_t *list2, FineGrid *in, REAL8 f1dot_fg, REAL8 f2dot_fg, UsefulStageVariables *usefulparams, REAL4 NSegmentsInv, REAL4 *NSegmentsInvX );
void GetSegsPosVelAccEarthOrb( LALStatus *status, REAL8VectorSequence **posSeg,
                               REAL8VectorSequence **velSeg, REAL8VectorSequence **accSeg,
                               UsefulStageVariables *usefulparams );
static inline INT4 ComputeU1idx( REAL8 freq_event, REAL8 f1dot_eventB1, REAL8 A1, REAL8 U1start, REAL8 U1winInv );
void ComputeU2idx( REAL8 freq_event, REAL8 f1dot_event, REAL8 A2, REAL8 B2, REAL8 U2start, REAL8 U2winInv,
                   INT4 *U2idx);
int compareCoarseGridUindex( const void *a, const void *b );
int compareFineGridNC( const void *a,const void *b );
int compareFineGridsumTwoF( const void *a,const void *b );

SFTCatalogSequence *XLALSetUpStacksFromSegmentList ( const SFTCatalog *SFTCatalog, const LALSegList *segList );

int XLALComputeFStatFreqBand (  MultiFstatFrequencySeries **fstatSeries,
				const PulsarDopplerParams *doppler,
				const MultiSFTVector *multiSFTs,
				const MultiNoiseWeights *multiWeights,
				const MultiDetectorStateSeries *multiDetStates,
				const ComputeFParams *params);

int XLALExtrapolateToplistPulsarSpins ( toplist_t *list,
					const LIGOTimeGPS usefulParamsRefTime,
					const LIGOTimeGPS finegridRefTime);

/* ---------- Global variables -------------------- */
LALStatus *global_status; /* a global pointer to MAIN()s head of the LALStatus structure */
extern int lalDebugLevel;
char *global_column_headings_stringp;

/* ###################################  MAIN  ################################### */

int MAIN( int argc, char *argv[]) {
  LALStatus status = blank_status;

  /* temp loop variables: generally k loops over segments and j over SFTs in a stack */
  UINT4 j;
  UINT4 k;
  UINT4 skyGridCounter; /* coarse sky position counter */
  UINT4 f1dotGridCounter; /* coarse f1dot position counter */

  /* ephemeris */
  EphemerisData *edat = NULL;

  /* GPS timestamp vectors */
  LIGOTimeGPSVector *midTstack = NULL;
  LIGOTimeGPSVector *startTstack = NULL;
  LIGOTimeGPSVector *endTstack = NULL;

  /* General GPS times */
  LIGOTimeGPS refTimeGPS = empty_LIGOTimeGPS;
  LIGOTimeGPS tMidGPS = empty_LIGOTimeGPS;

  /* GPS time used for each segment's midpoint */
  LIGOTimeGPS midTstackGPS = empty_LIGOTimeGPS;
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
  static MultiSFTVectorSequence stackMultiSFT;
  static MultiNoiseWeightsSequence stackMultiNoiseWeights;
  static MultiDetectorStateSeriesSequence stackMultiDetStates;
  static LIGOTimeGPS minStartTimeGPS, maxEndTimeGPS;
  SFTtype *firstSFT;
  REAL8 Tsft;

  /* some useful variables for each stage */
  UsefulStageVariables usefulParams;

  /* F-statistic computation related stuff */
  REAL4FrequencySeriesVector fstatVector; /* F-statistic vectors for each segment */
  UINT4 binsFstat1, binsFstatSearch=0;
  static ComputeFParams CFparams;
  ComputeFBufferVector_RS resampbuffers;  /* used to store the buffered quantities used in repeated calls to ComputeFstatFreqBand_RS */

  /* Semicoherent variables */
  static SemiCoherentParams semiCohPar;

  /* coarse grid */
  CoarseGrid coarsegrid;
  REAL8 dFreqStack; /* frequency resolution of Fstat calculation */
  REAL8 df1dot;  /* coarse grid resolution in spindown */
  UINT4 nf1dot;  /* number of coarse-grid spindown values */
  UINT4 ifdot;  /* counter for coarse-grid spindown values */
  REAL8 df2dot;  /* coarse grid resolution in 2nd spindown */
  UINT4 nf2dot;  /* number of coarse-grid 2nd spindown values */
  UINT4 if2dot;  /* counter for coarse-grid 2nd spindown values */

  /* fine grid */
  FineGrid finegrid;
  UINT4 nf1dots_fg = 1; /* number of frequency and spindown values */
  REAL8 gammaRefine, sigmasq;  /* refinement factor and variance */
  UINT4 nf2dots_fg=1;          /* number of second spindown values */
  REAL8 gamma2Refine, sigma4;  /* 2nd spindown refinement factor and 4th moment */

  /* GCT helper variables */
  UINT4 if1dot_fg, if2dot_fg;
  UINT4 fveclength, ifreq;
  INT4  U1idx;
  REAL8 myf0, freq_event, f1dot_event, deltaF, f1dot_eventB1;
  REAL8 dfreq_fg, df1dot_fg, freqmin_fg, f1dotmin_fg, freqband_fg;
  REAL8 df2dot_fg, f2dotmin_fg;
  REAL8 u1start, u1win, u1winInv;
  REAL8 freq_fg, f1dot_fg, f2dot_fg;
  REAL4 Fstat;
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
  toplist_t *semiCohToplist2=NULL;	// only used for SORTBY_DUAL_F_LV: 1st toplist sorted by 'F', 2nd one by 'LV'

  /* template and grid variables */
  static DopplerSkyScanInit scanInit;   /* init-structure for DopperScanner */
  DopplerSkyScanState thisScan = empty_DopplerSkyScanState; /* current state of the Doppler-scan */
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
  BOOLEAN uvar_help = FALSE;    /* true if -h option is given */
  BOOLEAN uvar_log = FALSE;     /* logging done if true */
  INT4 uvar_loglevel = 0;

  BOOLEAN uvar_printCand1 = FALSE;      /* if 1st stage candidates are to be printed */
  BOOLEAN uvar_printFstat1 = FALSE;
  BOOLEAN uvar_semiCohToplist = TRUE; /* if overall first stage candidates are to be output */
  BOOLEAN uvar_useResamp = FALSE;      /* use resampling to compute F-statistic instead of SFT method */
  BOOLEAN uvar_SignalOnly = FALSE;     /* if Signal-only case (for SFT normalization) */
  BOOLEAN uvar_recalcToplistStats = FALSE; /* Do additional analysis for all toplist candidates, output F, FXvector for postprocessing */
  BOOLEAN uvar_computeLV = FALSE;          /* In Fstat loop, get single-IFO F-stats [and, in future, compute Line Veto stat] */
  BOOLEAN uvar_LVuseAllTerms = TRUE;       /* Use only leading term or all terms in Line Veto computation */
  REAL8 uvar_LVrho = 0.0;                  /* Prior parameter rho_max_line for LineVeto statistic */
  LALStringVector *uvar_LVlX = NULL;       /* Line-to-gauss prior ratios lX for LineVeto statistic */
  CHAR *uvar_outputSingleSegStats = NULL; /* Additionally output single-segment Fstats for each final toplist candidate */

  REAL8 uvar_dAlpha = DALPHA;   /* resolution for flat or isotropic grids -- coarse grid*/
  REAL8 uvar_dDelta = DDELTA;
  REAL8 uvar_f1dot = FDOT;      /* first spindown value */
  REAL8 uvar_f1dotBand = DFDOT; /* range of first spindown parameter */
  REAL8 uvar_f2dot = F2DOT;     /* second spindown value */
  REAL8 uvar_f2dotBand = DF2DOT; /* range of second spindown parameter */
  REAL8 uvar_Freq = FSTART;
  REAL8 uvar_FreqBand = FBAND;

  REAL8 uvar_dFreq = 0;
  REAL8 uvar_df1dot = 0; /* coarse grid frequency and spindown resolution */
  REAL8 uvar_df2dot = 0; /* coarse grid second spindown resolution */

  REAL8 uvar_ThrF = FSTATTHRESHOLD; /* threshold of Fstat to select peaks */
  REAL8 uvar_mismatch1 = MISMATCH; /* metric mismatch for first stage coarse grid */

  REAL8 uvar_minStartTime1 = 0;
  REAL8 uvar_maxEndTime1 = LAL_INT4_MAX;
  REAL8 uvar_dopplerMax = 1.05e-4;

  REAL8 uvar_refTime = 0;
  INT4 uvar_nCand1 = NCAND1; /* number of candidates to be followed up from first stage */

  INT4 uvar_blocksRngMed = BLOCKSRNGMED;

  REAL8 uvar_tStack = 0;
  INT4  uvar_nStacksMax = 1;
  CHAR *uvar_segmentList = NULL;	/**< ALTERNATIVE: file containing a pre-computed segment list of tuples (startGPS endGPS duration[h] NumSFTs) */

  INT4 uvar_Dterms = DTERMS;
  INT4 uvar_SSBprecision = SSBPREC_RELATIVISTIC;
  INT4 uvar_gammaRefine = 1;
  INT4 uvar_gamma2Refine = 1;
  INT4 uvar_metricType1 = LAL_PMETRIC_COH_PTOLE_ANALYTIC;
  INT4 uvar_gridType1 = GRID_METRIC;
  INT4 uvar_sftUpsampling = 1;
  INT4 uvar_skyPointIndex = -1;

  CHAR *uvar_ephemE = NULL;
  CHAR *uvar_ephemS = NULL;

  CHAR *uvar_skyRegion = NULL;
  CHAR *uvar_fnameout = NULL;
  CHAR *uvar_DataFiles1 = NULL;
  CHAR *uvar_skyGridFile=NULL;
  INT4 uvar_numSkyPartitions = 0;
  INT4 uvar_partitionIndex = 0;
  INT4 uvar_SortToplist = 0;
  BOOLEAN uvar_version = 0;

  CHAR *uvar_outputTiming = NULL;

  BOOLEAN uvar_useWholeSFTs = 0;

  global_status = &status;

#ifndef EAH_BOINC
  global_argv = argv;
  global_argc = argc;
#endif

  /* LALDebugLevel must be called before any LALMallocs have been used */
  lalDebugLevel = 0;
#ifdef EAH_LALDEBUGLEVEL
  lalDebugLevel = EAH_LALDEBUGLEVEL;
#endif
  LAL_CALL( LALGetDebugLevel( &status, argc, argv, 'd'), &status);

  /* set log-level */
#ifdef EAH_LOGLEVEL
  uvar_loglevel = EAH_LOGLEVEL;
#else
  uvar_loglevel = lalDebugLevel & LALALLDBG;
#endif

  uvar_ephemE = LALCalloc( strlen( EARTHEPHEMERIS ) + 1, sizeof(CHAR) );
  strcpy(uvar_ephemE, EARTHEPHEMERIS);

  uvar_ephemS = LALCalloc( strlen(SUNEPHEMERIS) + 1, sizeof(CHAR) );
  strcpy(uvar_ephemS, SUNEPHEMERIS);

  uvar_skyRegion = LALCalloc( strlen(SKYREGION) + 1, sizeof(CHAR) );
  strcpy(uvar_skyRegion, SKYREGION);

  uvar_fnameout = LALCalloc( strlen(FNAMEOUT) + 1, sizeof(CHAR) );
  strcpy(uvar_fnameout, FNAMEOUT);

  /* set LAL error-handler */
#ifdef EAH_BOINC
  lal_errhandler = BOINC_LAL_ErrHand;
#else
  lal_errhandler = LAL_ERR_EXIT;
#endif

  /* register user input variables */
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "help",        'h', UVAR_HELP,     "Print this message", &uvar_help), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "logLevel",     0,  UVAR_OPTIONAL, "Set logLevel", &uvar_loglevel), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "log",          0,  UVAR_OPTIONAL, "Write log file", &uvar_log), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "semiCohToplist",0, UVAR_OPTIONAL, "Print toplist of semicoherent candidates", &uvar_semiCohToplist ), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "DataFiles1",   0,  UVAR_REQUIRED, "1st SFT file pattern", &uvar_DataFiles1), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "skyRegion",    0,  UVAR_OPTIONAL, "sky-region polygon (or 'allsky')", &uvar_skyRegion), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "numSkyPartitions",0,UVAR_OPTIONAL, "No. of (equi-)partitions to split skygrid into", &uvar_numSkyPartitions), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "partitionIndex",0,UVAR_OPTIONAL, "Index [0,numSkyPartitions-1] of sky-partition to generate", &uvar_partitionIndex), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "skyGridFile",  0,  UVAR_OPTIONAL, "sky-grid file", &uvar_skyGridFile), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dAlpha",       0,  UVAR_OPTIONAL, "Resolution for flat or isotropic coarse grid", &uvar_dAlpha), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dDelta",       0,  UVAR_OPTIONAL, "Resolution for flat or isotropic coarse grid", &uvar_dDelta), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "Freq",        'f', UVAR_OPTIONAL, "Start search frequency", &uvar_Freq), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dFreq",        0,  UVAR_OPTIONAL, "Frequency resolution (default \\propto 1/Tstack)", &uvar_dFreq), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "FreqBand",    'b', UVAR_OPTIONAL, "Search frequency band", &uvar_FreqBand), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "f1dot",        0,  UVAR_OPTIONAL, "Spindown parameter", &uvar_f1dot), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "df1dot",       0,  UVAR_OPTIONAL, "Spindown resolution (default \\propto 1/Tstack^2)", &uvar_df1dot), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "f1dotBand",    0,  UVAR_OPTIONAL, "Spindown Range", &uvar_f1dotBand), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "f2dot",        0,  UVAR_OPTIONAL, "2nd spindown parameter", &uvar_f2dot), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "df2dot",       0,  UVAR_OPTIONAL, "2nd spindown resolution (default \\propto 1/Tstack^3)", &uvar_df2dot), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "f2dotBand",    0,  UVAR_OPTIONAL, "2nd spindown Range", &uvar_f2dotBand), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "peakThrF",     0,  UVAR_OPTIONAL, "Fstat Threshold", &uvar_ThrF), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "mismatch1",   'm', UVAR_OPTIONAL, "1st stage mismatch", &uvar_mismatch1), &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "gridType1",    0,  UVAR_OPTIONAL, "0=flat,1=isotropic,2=metric,3=file", &uvar_gridType1),  &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "metricType1",  0,  UVAR_OPTIONAL, "0=none,1=Ptole-analytic,2=Ptole-numeric,3=exact", &uvar_metricType1), &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "gammaRefine", 'g', UVAR_OPTIONAL, "Refinement of fine grid (default: use segment times)", &uvar_gammaRefine), &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "gamma2Refine",'G', UVAR_OPTIONAL, "Refinement of f2dot fine grid (default: use segment times, -1=use gammaRefine)", &uvar_gamma2Refine), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "fnameout",    'o', UVAR_OPTIONAL, "Output filename", &uvar_fnameout), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "fnameChkPoint",0,  UVAR_OPTIONAL, "Checkpoint filename", &uvar_fnameChkPoint), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "nCand1",      'n', UVAR_OPTIONAL, "No. of candidates to output", &uvar_nCand1), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printCand1",   0,  UVAR_OPTIONAL, "Print 1st stage candidates", &uvar_printCand1), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "refTime",      0,  UVAR_OPTIONAL, "Ref. time for pulsar pars [Default: mid-time]", &uvar_refTime), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "ephemE",       0,  UVAR_OPTIONAL, "Location of Earth ephemeris file", &uvar_ephemE),  &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "ephemS",       0,  UVAR_OPTIONAL, "Location of Sun ephemeris file", &uvar_ephemS),  &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "minStartTime1",0,  UVAR_OPTIONAL, "1st stage min start time of observation", &uvar_minStartTime1), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "maxEndTime1",  0,  UVAR_OPTIONAL, "1st stage max end time of observation",   &uvar_maxEndTime1),   &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printFstat1",  0,  UVAR_OPTIONAL, "Print 1st stage Fstat vectors", &uvar_printFstat1), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "useResamp",    0,  UVAR_OPTIONAL, "Use resampling to compute F-statistic", &uvar_useResamp), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "SignalOnly",  'S', UVAR_OPTIONAL, "Signal only flag", &uvar_SignalOnly), &status);

  LAL_CALL( LALRegisterINTUserVar(    &status, "nStacksMax",   0,  UVAR_OPTIONAL, "Maximum No. of segments", &uvar_nStacksMax ),&status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "tStack",      'T', UVAR_OPTIONAL, "Duration of segments (sec)", &uvar_tStack ),&status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "segmentList",  0, UVAR_OPTIONAL, "ALTERNATIVE: file containing a segment list: lines of form <startGPS endGPS duration[h] NumSFTs>", &uvar_segmentList),  &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "recalcToplistStats", 0, UVAR_OPTIONAL, "Additional analysis for toplist candidates, recalculate 2F, 2FX at finegrid", &uvar_recalcToplistStats), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "computeLV",    0, UVAR_OPTIONAL,  "Compute LineVeto stat for all candidates from single- and multi-IFO F-stats, can be used as main toplist statistic", &uvar_computeLV), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "LVrho",        0, UVAR_OPTIONAL,  "LineVeto: Prior rho_max_line, must be >=0", &uvar_LVrho), &status);
  LAL_CALL( LALRegisterLISTUserVar(   &status, "LVlX",         0, UVAR_OPTIONAL,  "LineVeto: line-to-gauss prior ratios lX for different detectors X, length must be numDetectors. Defaults to lX=1,1,..", &uvar_LVlX), &status);

  /* developer user variables */
  LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed", 0, UVAR_DEVELOPER, "RngMed block size", &uvar_blocksRngMed), &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "SSBprecision", 0, UVAR_DEVELOPER, "Precision for SSB transform.", &uvar_SSBprecision),    &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "Dterms",       0, UVAR_DEVELOPER, "No. of terms to keep in Dirichlet Kernel", &uvar_Dterms ), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "skyPointIndex",0, UVAR_DEVELOPER, "Only analyze this skypoint in grid", &uvar_skyPointIndex ), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dopplerMax",   0, UVAR_DEVELOPER, "Max Doppler shift",  &uvar_dopplerMax), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "sftUpsampling",0, UVAR_DEVELOPER, "Upsampling factor for fast LALDemod",  &uvar_sftUpsampling), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "SortToplist",  0, UVAR_DEVELOPER, "Sort toplist by: 0=avg2F, 1=numbercount, 2=LV-stat, 3=dual-toplists 'avg2F+LV'",  &uvar_SortToplist), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "LVuseAllTerms",0, UVAR_DEVELOPER, "LineVeto: which terms to include - FALSE: only leading term, TRUE: all terms", &uvar_LVuseAllTerms), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "outputSingleSegStats", 0,  UVAR_DEVELOPER, "Base filename for single-segment Fstat output (1 file per final toplist candidate!)", &uvar_outputSingleSegStats),  &status);

  LAL_CALL( LALRegisterSTRINGUserVar( &status, "outputTiming", 0, UVAR_DEVELOPER, "Append timing information into this file", &uvar_outputTiming), &status);

  LAL_CALL( LALRegisterBOOLUserVar( &status, "useWholeSFTs", 0, UVAR_DEVELOPER, "Read in all SFTs bins (workaround for code searching outside input band)", &uvar_useWholeSFTs), &status);

  LAL_CALL ( LALRegisterBOOLUserVar(  &status, "version",     'V', UVAR_SPECIAL,  "Output version information", &uvar_version), &status);

  /* read all command line variables */
  LAL_CALL( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* set log level */
  LogSetLevel(uvar_loglevel);

  /* assemble version string */
  CHAR *VCSInfoString;
  if ( (VCSInfoString = XLALGetVersionString(0)) == NULL ) {
    XLALPrintError("XLALGetVersionString(0) failed.\n");
    return HIERARCHICALSEARCH_ESUB;
  }

  LogPrintfVerbatim( LOG_DEBUG, "Code-version: %s\n", VCSInfoString );

  if ( uvar_version )
    {
      printf ("%s\n", VCSInfoString );
      return (0);
    }

  /* exit if help was required */
  if (uvar_help)
    return(0);


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

  /* 2F threshold for semicoherent stage */
#ifndef EXP_NO_NUM_COUNT
  REAL4 TwoFthreshold = 2.0 * uvar_ThrF;
#endif

  if ( (uvar_SortToplist < 0) || (uvar_SortToplist >= SORTBY_LAST) ) {
    XLALPrintError ( "Invalid value %d specified for toplist sorting, must be within [0, %d]\n", uvar_SortToplist, SORTBY_LAST - 1 );
    return( HIERARCHICALSEARCH_EBAD );
  }
  if ( (uvar_SortToplist == SORTBY_LV || uvar_SortToplist == SORTBY_DUAL_F_LV) && !uvar_computeLV ) {
    fprintf(stderr, "Toplist sorting by LV-stat only possible if --computeLV given.\n");
    return( HIERARCHICALSEARCH_EBAD );
  }

  /* take LV user vars and save them in usefulParams */
  usefulParams.LVuseAllTerms = uvar_LVuseAllTerms;
  usefulParams.LVloglX = NULL;
  if ( uvar_LVrho < 0.0 ) {
    fprintf(stderr, "Invalid LV prior rho (given rho=%f, need rho>=0)!\n", uvar_LVrho);
    return( HIERARCHICALSEARCH_EBAD );
  }
  else if ( uvar_LVrho > 0.0 )
    usefulParams.LVlogRhoTerm = 4.0 * log(uvar_LVrho) - log(70.0);
  else /* if uvar_LVrho == 0.0, logRhoTerm should become irrelevant in summation */
    usefulParams.LVlogRhoTerm = - LAL_REAL4_MAX;

  /* create toplist -- semiCohToplist has the same structure
     as a fstat candidate, so treat it as a fstat candidate */
  if ( uvar_SortToplist == SORTBY_DUAL_F_LV )	// special treatement of 'dual' toplists: 1st one sorted by 'F', 2nd one by 'LV'
    {
      XLAL_CHECK ( 0 == create_gctFStat_toplist ( &semiCohToplist, uvar_nCand1, SORTBY_F ),
                   XLAL_EFUNC, "create_gctFStat_toplist() failed for nCand=%d and sortBy=%d\n", uvar_nCand1, SORTBY_F );
      XLAL_CHECK ( 0 == create_gctFStat_toplist ( &semiCohToplist2, uvar_nCand1, SORTBY_LV ),
                   XLAL_EFUNC, "create_gctFStat_toplist() failed for nCand=%d and sortBy=%d\n", uvar_nCand1, SORTBY_LV );
    }
  else	// 'normal' single-sorting toplist cases (sortby 'F', 'nc' or 'LV')
    {
      XLAL_CHECK ( 0 == create_gctFStat_toplist ( &semiCohToplist, uvar_nCand1, uvar_SortToplist),
                   XLAL_EFUNC, "create_gctFStat_toplist() failed for nCand=%d and sortBy=%d\n", uvar_nCand1, uvar_SortToplist );
    }

#ifdef EAH_BOINC
  // BOINC Apps always checkpoint, so set a default filename here
  if (uvar_fnameChkPoint == NULL) {
    CHAR*fname = "checkpoint.cpt";
    uvar_fnameChkPoint = XLALMalloc(strlen(fname)+1);
    if (uvar_fnameChkPoint == NULL) {
      fprintf(stderr, "error allocating memory [HierarchSearchGCT.c %d]\n" , __LINE__);
      return(HIERARCHICALSEARCH_EMEM);
    }
    strcpy(uvar_fnameChkPoint, fname);
  }
#endif

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
      LAL_CALL( LALUserVarGetLog(&status, &logstr, UVAR_LOGFMT_CFGFILE), &status);

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

  /* prepare timing-file: write header-line with columns headings */
  if ( uvar_outputTiming )
    {
      FILE *timing_fp;
      if ( (timing_fp = fopen ( uvar_outputTiming, "ab" )) == NULL ) {
        XLALPrintError ("Failed to open timing file '%s' for writing/appending.\n", uvar_outputTiming );
        return HIERARCHICALSEARCH_EFILE;
      }
      /* write column headings */
      fprintf ( timing_fp, "%10s %10s %10s %7s %7s    %9s %9s %9s %9s    %10s %10s\n",
                "%% Ncoarse", "NSB", "Nfine", "Nsft", "Nseg", "tauWU[s]", "tauCo[s]", "tauIc[s]", "tauLV[s]", "c0co[s]", "c0ic[s]" );
      fclose ( timing_fp );
    } /* if outputTiming */

  /* initializations of coarse and fine grids */
  coarsegrid.TwoF=NULL;
  coarsegrid.TwoFX=NULL;
  coarsegrid.Uindex=NULL;
  finegrid.nc= NULL;
  finegrid.sumTwoF=NULL;
  finegrid.sumTwoFX=NULL;
  Fstat = 0;

  /* initialize ephemeris info */
  edat = (EphemerisData *)LALCalloc(1, sizeof(EphemerisData));
  if ( edat == NULL) {
    fprintf(stderr, "error allocating memory [HierarchSearchGCT.c %d]\n" , __LINE__);
    return(HIERARCHICALSEARCH_EMEM);
  }
  (*edat).ephiles.earthEphemeris = uvar_ephemE;
  (*edat).ephiles.sunEphemeris = uvar_ephemS;

  /* read in ephemeris data */
  LAL_CALL( LALInitBarycenter( &status, edat), &status);

  XLALGPSSetREAL8(&minStartTimeGPS, uvar_minStartTime1);
  XLALGPSSetREAL8(&maxEndTimeGPS, uvar_maxEndTime1);

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
  INIT_MEM( spinRange_Temp );

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


  INIT_MEM ( usefulParams.spinRange_startTime );
  INIT_MEM ( usefulParams.spinRange_endTime );
  INIT_MEM ( usefulParams.spinRange_refTime );
  INIT_MEM ( usefulParams.spinRange_midTime );

  /* copy user specified spin variables at reftime  */
  /* the reference time value in spinRange_refTime will be set in SetUpSFTs() */
  usefulParams.spinRange_refTime.fkdot[0] = uvar_Freq; /* frequency */
  usefulParams.spinRange_refTime.fkdot[1] = uvar_f1dot;  /* 1st spindown */
  usefulParams.spinRange_refTime.fkdot[2] = uvar_f2dot;  /* 2nd spindown */
  usefulParams.spinRange_refTime.fkdotBand[0] = uvar_FreqBand; /* frequency range */
  usefulParams.spinRange_refTime.fkdotBand[1] = uvar_f1dotBand; /* spindown range */
  usefulParams.spinRange_refTime.fkdotBand[2] = uvar_f2dotBand; /* spindown range */

  usefulParams.edat = edat;
  usefulParams.minStartTimeGPS = minStartTimeGPS;
  usefulParams.maxEndTimeGPS = maxEndTimeGPS;
  usefulParams.blocksRngMed = uvar_blocksRngMed;
  usefulParams.Dterms = uvar_Dterms;
  usefulParams.SignalOnly = uvar_SignalOnly;
  usefulParams.dopplerMax = uvar_dopplerMax;

  /* set reference time for pulsar parameters */
  if ( LALUserVarWasSet(&uvar_refTime))
    usefulParams.refTime = uvar_refTime;
  else {
    XLALPrintWarning("Reference time will be set to mid-time of observation time\n");
    usefulParams.refTime = -1;
  }

  /* set Fstat calculation frequency resolution (coarse grid) */
  if ( LALUserVarWasSet(&uvar_FreqBand) ) {
    if ( LALUserVarWasSet(&uvar_dFreq) ) {
      usefulParams.dFreqStack = uvar_dFreq;
    } else {
      XLALPrintError("--dFreq is required if --FreqBand is given\n");
      return( HIERARCHICALSEARCH_EBAD );
    }
  } else {
    usefulParams.dFreqStack = 0;
  }

  /* set Fstat spindown resolution (coarse grid) */
  if ( LALUserVarWasSet(&uvar_f1dotBand) ) {
    if ( LALUserVarWasSet(&uvar_df1dot) ) {
      usefulParams.df1dot = uvar_df1dot;
    } else {
      XLALPrintError("--df1dot is required if --f1dotBand is given\n");
      return( HIERARCHICALSEARCH_EBAD );
    }
  } else {
    usefulParams.df1dot = 0;
  }

  /* set Fstat 2nd spindown resolution (coarse grid) */
  if ( LALUserVarWasSet(&uvar_f2dotBand) ) {
    if ( LALUserVarWasSet(&uvar_df2dot) ) {
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

  /* for 1st stage: read sfts, calculate detector states */
  LogPrintf( LOG_NORMAL,"Reading input data ... ");
  LAL_CALL( SetUpSFTs( &status, &stackMultiSFT, &stackMultiNoiseWeights, &stackMultiDetStates, &usefulParams, uvar_useWholeSFTs, uvar_mismatch1), &status);
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
  firstSFT = &(stackMultiSFT.data[0]->data[0]->data[0]); /* use  first SFT from  first detector */
  Tsft = 1.0 / firstSFT->deltaF; /* define the length of an SFT (assuming 1/Tsft resolution) */

  /* count the total and per-segment number of SFTs used */
  UINT4 iTS, nSFTs = 0;
  for ( iTS = 0; iTS < nStacks; iTS ++ )
    {
      UINT4 nSFTsInSeg = 0;
      for ( UINT4 X=0; X < stackMultiSFT.data[iTS]->length; X ++ )
        nSFTsInSeg += stackMultiSFT.data[iTS]->data[X]->length;
      nSFTs += nSFTsInSeg;
      /* if we have a segment-list: double-check number of SFTs */
      if ( usefulParams.segmentList )
        {
          /* check the number of SFTs we found in this segment against the nominal value,
           * stored in the segment list field 'id' */
          UINT4 nSFTsExpected = usefulParams.segmentList->segs[iTS].id;
          if ( nSFTsInSeg != nSFTsExpected ) {
            XLALPrintError ("%s: Segment list seems inconsistent with data read: segment %d contains %d SFTs, should hold %d SFTs\n", __func__, iTS, nSFTsInSeg, nSFTsExpected );
            XLAL_ERROR ( XLAL_EDOM );
          }

        } /* if have segmentList */

    } /* for iTS < nStacks */
  XLALPrintWarning ("Number of segments: %d, total number of SFTs in segments: %d\n", nStacks, nSFTs );

  /* free segment list */
  if ( usefulParams.segmentList )
    if ( XLALSegListClear( usefulParams.segmentList ) != XLAL_SUCCESS )
      XLAL_ERROR ( XLAL_EFUNC );
  XLALFree ( usefulParams.segmentList );
  usefulParams.segmentList = NULL;


  /* special treatment of SFTs if upsampling is used */
  if ( uvar_sftUpsampling > 1 )
    {
      LogPrintf (LOG_DEBUG, "Upsampling SFTs by factor %d ... ", uvar_sftUpsampling );
      for (k = 0; k < nStacks; k++) {
        LAL_CALL ( upsampleMultiSFTVector ( &status, stackMultiSFT.data[k], uvar_sftUpsampling, 16 ), &status );
      }
      LogPrintfVerbatim (LOG_DEBUG, "done.\n");
    }

  /*------- set frequency and spindown resolutions and ranges for Fstat and semicoherent steps -----*/

  dFreqStack = usefulParams.dFreqStack;
  df1dot = usefulParams.df1dot;
  df2dot = usefulParams.df2dot;
  LogPrintf(LOG_NORMAL, "dFreqStack = %e, df1dot = %e, df2dot = %e\n", dFreqStack, df1dot, df2dot);

  /* number of coarse grid spindown values */
  if ( df1dot == 0 ) {
    nf1dot = 1;
  } else {
    nf1dot = (UINT4) ceil( usefulParams.spinRange_midTime.fkdotBand[1] / df1dot) + 1;
  }

  /* set number of fine-grid spindowns */
  if ( LALUserVarWasSet(&uvar_gammaRefine) ) {
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

  /* number of coarse grid 2nd spindown values */
  if ( df2dot == 0 ) {
    nf2dot = 1;
  } else {
    nf2dot = (UINT4) floor( usefulParams.spinRange_midTime.fkdotBand[2] / df2dot - NUDGE) + 1;
  }

  /* set number of fine-grid 2nd spindowns */
  if ( LALUserVarWasSet(&uvar_gamma2Refine) ) {
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
  LogPrintf(LOG_DETAIL, "Frequency and spindown range at refTime (%d): [%f,%f], [%e,%e], [%e,%e]\n",
            usefulParams.spinRange_refTime.refTime.gpsSeconds,
            usefulParams.spinRange_refTime.fkdot[0],
            usefulParams.spinRange_refTime.fkdot[0] + usefulParams.spinRange_refTime.fkdotBand[0],
            usefulParams.spinRange_refTime.fkdot[1],
            usefulParams.spinRange_refTime.fkdot[1] + usefulParams.spinRange_refTime.fkdotBand[1],
            usefulParams.spinRange_refTime.fkdot[2],
            usefulParams.spinRange_refTime.fkdot[2] + usefulParams.spinRange_refTime.fkdotBand[2]);

  LogPrintf(LOG_DETAIL, "Frequency and spindown range at startTime (%d): [%f,%f], [%e,%e], [%e,%e]\n",
            usefulParams.spinRange_startTime.refTime.gpsSeconds,
            usefulParams.spinRange_startTime.fkdot[0],
            usefulParams.spinRange_startTime.fkdot[0] + usefulParams.spinRange_startTime.fkdotBand[0],
            usefulParams.spinRange_startTime.fkdot[1],
            usefulParams.spinRange_startTime.fkdot[1] + usefulParams.spinRange_startTime.fkdotBand[1],
            usefulParams.spinRange_startTime.fkdot[2],
            usefulParams.spinRange_startTime.fkdot[2] + usefulParams.spinRange_startTime.fkdotBand[2]);

  LogPrintf(LOG_DETAIL, "Frequency and spindown range at midTime (%d): [%f,%f], [%e,%e], [%e,%e]\n",
            usefulParams.spinRange_midTime.refTime.gpsSeconds,
            usefulParams.spinRange_midTime.fkdot[0],
            usefulParams.spinRange_midTime.fkdot[0] + usefulParams.spinRange_midTime.fkdotBand[0],
            usefulParams.spinRange_midTime.fkdot[1],
            usefulParams.spinRange_midTime.fkdot[1] + usefulParams.spinRange_midTime.fkdotBand[1],
            usefulParams.spinRange_midTime.fkdot[2],
            usefulParams.spinRange_midTime.fkdot[2] + usefulParams.spinRange_midTime.fkdotBand[2]);

  LogPrintf(LOG_DETAIL, "Frequency and spindown range at endTime (%d): [%f,%f], [%e,%e], [%e,%e]\n",
            usefulParams.spinRange_endTime.refTime.gpsSeconds,
            usefulParams.spinRange_endTime.fkdot[0],
            usefulParams.spinRange_endTime.fkdot[0] + usefulParams.spinRange_endTime.fkdotBand[0],
            usefulParams.spinRange_endTime.fkdot[1],
            usefulParams.spinRange_endTime.fkdot[1] + usefulParams.spinRange_endTime.fkdotBand[1],
            usefulParams.spinRange_endTime.fkdot[2],
            usefulParams.spinRange_endTime.fkdot[2] + usefulParams.spinRange_endTime.fkdotBand[2]);

  /* print debug info about stacks */
  fprintf(stderr, "%% --- Setup, N = %d, T = %.0f s, Tobs = %.0f s, gammaRefine = %.0f, gamma2Refine = %.0f\n",
          nStacks, tStack, tObs, gammaRefine, gamma2Refine);

  for (k = 0; k < nStacks; k++) {

    LogPrintf(LOG_DETAIL, "Segment %d ", k+1);
    for ( j = 0; j < stackMultiSFT.data[k]->length; j++) {

      INT4 tmpVar = stackMultiSFT.data[k]->data[j]->length;
      LogPrintfVerbatim(LOG_DETAIL, "%s: %d  ", stackMultiSFT.data[k]->data[j]->data[0].name, tmpVar);
    } /* loop over ifos */

    LogPrintfVerbatim(LOG_DETAIL, "\n");

  } /* loop over segments */



  /*---------- set up F-statistic calculation stuff ---------*/
  /* set reference time for calculating Fstatistic */
  thisPoint.refTime = tMidGPS; /* midpoint of data spanned */

  /* binary orbit and higher spindowns not considered */
  thisPoint.orbit = NULL;
  INIT_MEM ( thisPoint.fkdot );

  /* some compute F-Stat params */
  CFparams.Dterms = uvar_Dterms;
  CFparams.SSBprec = uvar_SSBprecision;
  CFparams.upsampling = uvar_sftUpsampling;
  CFparams.edat = edat;

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

  /* allocate some fstat memory */
  fstatVector.length = nStacks; /* for EACH segment generate a fstat-vector */
  fstatVector.data = NULL;
  fstatVector.data = (REAL4FrequencySeries *)LALCalloc( 1, nStacks * sizeof(REAL4FrequencySeries));
  if ( fstatVector.data == NULL) {
    fprintf(stderr, "error allocating memory [HierarchSearchGCT.c %d]\n" , __LINE__);
    return(HIERARCHICALSEARCH_EMEM);
  }

  /* allocate buffer memory for resampling - initialise first */
  resampbuffers.length = nStacks;
  resampbuffers.data = NULL;
  if (uvar_useResamp) {
    if ( (resampbuffers.data = (ComputeFBuffer_RS **)XLALCalloc(nStacks,sizeof(ComputeFBuffer_RS *))) == NULL ) {
      fprintf(stderr, "error allocating memory [HierarchSearchGCT.c %d]\n" , __LINE__);
      return(HIERARCHICALSEARCH_EMEM);
    }
  }

  /* get numer of detectors and, in line veto case, detector name vector */
  UINT4 numDetectors = 0;
  LALStringVector *detectorIDs = NULL;

  /* fill detector name vector with all detectors present in any data sements */
  if ( ( detectorIDs = XLALGetDetectorIDs ( &stackMultiSFT )) == NULL )
    XLAL_ERROR ( XLAL_EFUNC );
  numDetectors = detectorIDs->length;

  /* assemble column headings string for output file */
  char colum_headings_string_base[] = "freq alpha delta f1dot f2dot nc <2F>";
  UINT4 column_headings_string_length = sizeof(colum_headings_string_base);
  if ( uvar_computeLV ) {
    column_headings_string_length += 3 + numDetectors*8; /* 3 for " LV" and 8 per detector for " <2F_XY>" */
  }
  if ( uvar_recalcToplistStats ) {
    column_headings_string_length += 6 + numDetectors*9; /* 6 for " <2Fr>" and 9 per detector for " <2Fr_XY>" */
  }
  char column_headings_string[column_headings_string_length];
  INIT_MEM( column_headings_string );
  strcat ( column_headings_string, colum_headings_string_base );
  if ( uvar_computeLV ) {
    strcat ( column_headings_string, " LV" );
    for ( UINT4 X = 0; X < numDetectors ; X ++ ) {
      char headingX[9];
      snprintf ( headingX, sizeof(headingX), " <2F_%s>", detectorIDs->data[X] );
      strcat ( column_headings_string, headingX );
    } /* for X < numDet */
  }
  if ( uvar_recalcToplistStats ) {
    strcat ( column_headings_string, " <2Fr>" );
    for ( UINT4 X = 0; X < numDetectors ; X ++ ) {
      char headingX[10];
      snprintf ( headingX, sizeof(headingX), " <2Fr_%s>", detectorIDs->data[X] );
      strcat ( column_headings_string, headingX );
    } /* for X < numDet */
  }
  global_column_headings_stringp = column_headings_string;

  /* get effective inverse number of segments per detector (needed for correct averaging in single-IFO F calculation) */
  REAL4 * NSegmentsInvX = NULL;
  if ( uvar_computeLV )
    {
      NSegmentsInvX = ALRealloc( NSegmentsInvX, numDetectors * sizeof(*NSegmentsInvX));
      for (UINT4 X = 0; X < numDetectors; X++)
        {
          NSegmentsInvX[X] = 0;
          for (k = 0; k < nStacks; k++)
            { /* for each detector, check if present in each segment, and save the number of segments where it is */
              for (UINT4 Y = 0; Y < stackMultiSFT.data[k]->length; Y++)
                {
                  if ( strcmp( stackMultiSFT.data[k]->data[Y]->data[0].name, detectorIDs->data[X] ) == 0 )
                    NSegmentsInvX[X] += 1;
                } /* for Y < numDetectors */
            } /* for k < nStacks */
          NSegmentsInvX[X] = 1.0 / NSegmentsInvX[X]; /* now it is the inverse number */
        } /* for X < numDetectors */
    } /* if ( uvar_computeLV ) */
  REAL4 NSegmentsInv = 1.0 / nStacks; /* also need this for multi-detector F-stat averaging later on */

  /* set up line prior ratios: either given by user, then convert from string to REAL4 vector; else, pass NULL, which is interpreted as lX=1.0 for all X */
  usefulParams.LVloglX = NULL;
  if ( uvar_computeLV && uvar_LVlX ) {
    if (  uvar_LVlX->length != numDetectors ) {
      fprintf(stderr, "Length of LV prior ratio vector does not match number of detectors! (%d != %d)\n", uvar_LVlX->length, numDetectors);
      return( HIERARCHICALSEARCH_EBAD );
    }
    if ( (usefulParams.LVloglX = XLALCreateREAL8Vector ( numDetectors )) == NULL ) {
      fprintf(stderr, "Failed call to XLALCreateREAL8Vector( %d )\n", numDetectors );
      return( HIERARCHICALSEARCH_EXLAL );
    }
    for (UINT4 X = 0; X < numDetectors; X++) {
      if ( 1 != sscanf ( uvar_LVlX->data[X], "%" LAL_REAL8_FORMAT, &usefulParams.LVloglX->data[X] ) ) {
        fprintf(stderr, "Illegal REAL8 commandline argument to --LVlX[%d]: '%s'\n", X, uvar_LVlX->data[X]);
        return ( HIERARCHICALSEARCH_EBAD );
      }
      if ( usefulParams.LVloglX->data[X] < 0.0 ) {
        fprintf(stderr, "Negative input prior-ratio for detector X=%d lX[X]=%f\n", X, usefulParams.LVloglX->data[X] );
        return( HIERARCHICALSEARCH_EBAD );
      }
      else if ( usefulParams.LVloglX->data[X] > 0.0 )
        usefulParams.LVloglX->data[X] = log(usefulParams.LVloglX->data[X]);
      else /* if zero prior ratio, approximate log(0)=-inf by -LAL_REA4_MAX to avoid raising underflow exceptions */
        usefulParams.LVloglX->data[X] = - LAL_REAL8_MAX;
    } /* for X < numDetectors */
  } /* if ( computeLV && uvar_LVlX ) */

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
  scanInit.Detector = &(stackMultiDetStates.data[0]->data[0]->detector); /* Only used if metric is employed */
  scanInit.ephemeris = edat;
  scanInit.skyGridFile = uvar_skyGridFile;
  scanInit.skyRegionString = (CHAR*)LALCalloc(1, strlen(uvar_skyRegion)+1);
  if ( scanInit.skyRegionString == NULL) {
    fprintf(stderr, "error allocating memory [HierarchSearchGCT.c.c %d]\n" , __LINE__);
    return(HIERARCHICALSEARCH_EMEM);
  }
  strcpy (scanInit.skyRegionString, uvar_skyRegion);

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

    GET_GCT_CHECKPOINT (uvar_fnameChkPoint, semiCohToplist, semiCohToplist2, &count);

    if (count) {
      f1dotGridCounter = (UINT4) (count % nf1dot);  /* Checkpointing counter = i_sky * nf1dot + i_f1dot */
      skycount = (UINT4) ((count - f1dotGridCounter) / nf1dot);
    }
   fprintf (stderr, "%% --- Cpt:%d,  total:%d,  sky:%d/%d,  f1dot:%d/%d\n",
             count, thisScan.numSkyGridPoints*nf1dot, skycount+1, thisScan.numSkyGridPoints, f1dotGridCounter+1, nf1dot);

    for(skyGridCounter = 0; skyGridCounter < skycount; skyGridCounter++)
      XLALNextDopplerSkyPos(&dopplerpos, &thisScan);

    if ( count == thisScan.numSkyGridPoints*nf1dot )
      thisScan.state = STATE_FINISHED;

  }

  /* spool forward if uvar_skyPointIndex is set
     ---This probably doesn't make sense when checkpointing is turned on */
  if ( LALUserVarWasSet(&uvar_skyPointIndex)) {
    UINT4 count = uvar_skyPointIndex;
    for(skyGridCounter = 0; (skyGridCounter < count)&&(thisScan.state != STATE_FINISHED) ; skyGridCounter++)
      XLALNextDopplerSkyPos(&dopplerpos, &thisScan);
  }


  /* timing */
  REAL8 timeStart = 0.0, timeEnd = 0.0;
  REAL8 coherentTime = 0.0, incoherentTime = 0.0, vetoTime = 0.0;
  REAL8 timeStamp1 = 0.0, timeStamp2 = 0.0;
  if ( uvar_outputTiming )
    timeStart = XLALGetTimeOfDay();

  /* ################## loop over SKY coarse-grid points ################## */
  while(thisScan.state != STATE_FINISHED)
    {

      SHOW_PROGRESS(dopplerpos.Alpha, dopplerpos.Delta,
                    skyGridCounter * nf1dot + f1dotGridCounter,
                    thisScan.numSkyGridPoints * nf1dot, uvar_Freq, uvar_FreqBand);

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

      {  /********Allocate fstat vector memory *****************/

        /* calculate number of bins for Fstat overhead due to residual spin-down */
        semiCohPar.extraBinsFstat = usefulParams.extraBinsFstat;

        /* calculate total number of bins for Fstat */
        if ( dFreqStack == 0 ) {
          binsFstatSearch = 1;
        } else {
          binsFstatSearch = (UINT4)(usefulParams.spinRange_midTime.fkdotBand[0]/dFreqStack + 1e-6) + 1;
        }
        binsFstat1 = binsFstatSearch + 2 * semiCohPar.extraBinsFstat;

        /* loop over segments for memory allocation */
        for (k = 0; k < nStacks; k++) {

          /* watch out: the epoch here is not the reference time for f0! */
          fstatVector.data[k].epoch = startTstack->data[k];
          fstatVector.data[k].deltaF = dFreqStack;
          fstatVector.data[k].f0 = usefulParams.spinRange_midTime.fkdot[0] - semiCohPar.extraBinsFstat * dFreqStack;

          if (fstatVector.data[k].data == NULL) {
            fstatVector.data[k].data = (REAL4Sequence *)LALCalloc( 1, sizeof(REAL4Sequence));
            if ( fstatVector.data[k].data == NULL) {
              fprintf(stderr, "ERROR: Memory allocation  [HierarchSearchGCT.c %d]\n" , __LINE__);
              return(HIERARCHICALSEARCH_EMEM);
            }

            fstatVector.data[k].data->length = binsFstat1;
            fstatVector.data[k].data->data = (REAL4 *)LALCalloc( 1, binsFstat1 * sizeof(REAL4));
            if ( fstatVector.data[k].data->data == NULL) {
              fprintf(stderr, "ERROR: Memory allocation  [HierarchSearchGCT.c %d]\n" , __LINE__);
              return(HIERARCHICALSEARCH_EMEM);
            }

          }
          else {
            fstatVector.data[k].data = (REAL4Sequence *)LALRealloc( fstatVector.data[k].data, sizeof(REAL4Sequence));
            if ( fstatVector.data[k].data == NULL) {
              fprintf(stderr, "ERROR: Memory allocation  [HierarchSearchGCT.c %d]\n" , __LINE__);
              return(HIERARCHICALSEARCH_EMEM);
            }

            fstatVector.data[k].data->length = binsFstat1;
            fstatVector.data[k].data->data = (REAL4 *)LALRealloc( fstatVector.data[k].data->data, binsFstat1 * sizeof(REAL4));
            if ( fstatVector.data[k].data->data == NULL) {
              fprintf(stderr, "ERROR: Memory allocation  [HierarchSearchGCT.c %d]\n" , __LINE__);
              return(HIERARCHICALSEARCH_EMEM);
            }
          }
        } /* loop over segments */
      } /* fstat memory allocation block */



      /* ################## loop over coarse-grid F1DOT values ################## */
      ifdot = 0;

      while ( ifdot < nf1dot ) {

        /* if checkpoint read, spool forward */
        if (f1dotGridCounter > 0) {
          ifdot = f1dotGridCounter;
          f1dotGridCounter = 0;
        }

        /* ################## loop over coarse-grid F2DOT values ################## */
        if2dot = 0;

        while ( if2dot < nf2dot ) {

          /* show progress */
          LogPrintf( LOG_NORMAL, "Coarse grid sky:%d/%d f1dot:%d/%d f2dot:%d/%d\n", skyGridCounter+1, thisScan.numSkyGridPoints, ifdot+1, nf1dot, if2dot+1, nf2dot );

          /* ------------- Set up coarse grid --------------------------------------*/
          coarsegrid.freqlength = (UINT4) (binsFstat1);
          coarsegrid.nStacks = nStacks;
          coarsegrid.length = coarsegrid.freqlength * coarsegrid.nStacks;
          coarsegrid.numDetectors = numDetectors;

          /* allocate memory for coarsegrid */
          coarsegrid.TwoF = (REAL4 *)LALRealloc( coarsegrid.TwoF, coarsegrid.length * sizeof(REAL4));
	  if ( uvar_computeLV )
            coarsegrid.TwoFX = (REAL4 *)LALRealloc( coarsegrid.TwoFX, coarsegrid.numDetectors * coarsegrid.length * sizeof(REAL4));
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
          if (nf1dot == 1) {
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

          /* total number of fine-grid points */
          finegrid.length = finegrid.freqlength;

          if(!oldcg) {
            oldcg = coarsegrid.length;
            LogPrintfVerbatim(LOG_NORMAL, "%% --- CG:%d ",coarsegrid.length);
          }
          if(!oldfg) {
            oldfg = finegrid.length;
            LogPrintfVerbatim(LOG_NORMAL, "FG:%ld  f1dotmin_fg:%.13g df1dot_fg:%.13g f2dotmin_fg:%.13g df2dot_fg:%.13g\n",
                              finegrid.length,f1dotmin_fg,df1dot_fg,f2dotmin_fg,df2dot_fg);
          }
          if((coarsegrid.length != oldcg) || (finegrid.length != oldfg)) {
            LogPrintfVerbatim(LOG_CRITICAL, "ERROR: Grid-sizes disagree!\nPrevious CG:%d FG:%ld, currently CG:%d FG:%ld\n",
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
	  if ( uvar_computeLV )
            finegrid.sumTwoFX = (REAL4 *)ALRealloc( finegrid.sumTwoFX, finegrid.numDetectors * finegrid.freqlengthAL * sizeof(REAL4));

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

              /* initialize the entire finegrid ( 2F-sum and number count set to 0 ) */
              memset( finegrid.nc, 0, finegrid.length * sizeof(FINEGRID_NC_T) );
              memset( finegrid.sumTwoF, 0, finegrid.length * sizeof(REAL4) );
	      if ( uvar_computeLV )
                memset( finegrid.sumTwoFX, 0, finegrid.numDetectors * finegrid.freqlengthAL * sizeof(REAL4) );

              /* compute F-statistic values for coarse grid the first time through fine grid fdots loop */
              const BOOLEAN doComputeFstats = ( (if1dot_fg == 0) && (if2dot_fg == 0) );

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
                       + pos[2] * nvec[2] ); /* This is \vec r \dot \vec n */

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
                thisPoint.fkdot[0] = fstatVector.data[k].f0;

                /* Length and spacing of the Fstat vector in frequency */
                fveclength = fstatVector.data[k].data->length;
                deltaF = fstatVector.data[k].deltaF;

                /* Set spindown value for Fstat calculation */
                thisPoint.fkdot[1] = usefulParams.spinRange_midTime.fkdot[1] + ifdot * df1dot;

                /* Set spindown value for Fstat calculation */
                thisPoint.fkdot[2] = usefulParams.spinRange_midTime.fkdot[2] + if2dot * df2dot;

                /* Frequency at the segment's midpoint for later use */
                f1dot_event = thisPoint.fkdot[1] + thisPoint.fkdot[2] * timeDiffSeg;
                myf0 = thisPoint.fkdot[0] + thisPoint.fkdot[1] * timeDiffSeg +
                  + 0.5 * thisPoint.fkdot[2] * timeDiffSeg * timeDiffSeg;

                /* Smallest values of u1 and u2 (to be subtracted) */
                u1start = myf0 * A1 + f1dot_event * B1; /* Eq. (1a) */
                U1idx = 0;

                /* Holger: current code structure of loops (processing f1dot by f1dot) needs only U1 calculation.
                   u2start = f1dot_event + myf0 * A2 + 2.0 * f1dot_event * B2;
                   myf0max = myf0 + (fveclength - 1) * deltaF;
                   u2end = f1dot_event + myf0max * A2 + 2.0 * f1dot_event * B2;
                   NumU2idx = ceil(fabs(u2start - u2end) * u2winInv);
                   U2idx = 0;
                */

                /* pre-compute product */
                f1dot_eventB1 = f1dot_fg * B1;

                /* timing */
                if ( uvar_outputTiming )
                  timeStamp1 = XLALGetTimeOfDay();

                /* ----------------------------------------------------------------- */
                /************************ Compute F-Statistic ************************/
                if (doComputeFstats) { /* if first time through fine grid fdots loop */

                  /* prepare different Fstat structure for uvar_computeLV case */
                  MultiFstatFrequencySeries *multiFstatVector = NULL;

                  if (uvar_useResamp) {

                    /* point the params buffer to the current segment buffer */
                    CFparams.buffer = resampbuffers.data[k];
                    /* Resampling method implementation to compute the F-statistic */
                    LAL_CALL( COMPUTEFSTATFREQBAND_RS ( &status, &fstatVector.data[k], &thisPoint,
                                                        stackMultiSFT.data[k], stackMultiNoiseWeights.data[k],
                                                        &CFparams), &status);

                    /* repoint the buffer vector element to the potentially modified buffer */
                    resampbuffers.data[k] = CFparams.buffer;

                  }
                  else if (uvar_computeLV) {
                    thisPoint.dFreq = dFreqStack;
                    thisPoint.numFreqBins = binsFstat1;
                    CFparams.returnSingleF = TRUE;
                    xlalErrno = 0;
                    XLALComputeFStatFreqBand ( &multiFstatVector, &thisPoint, stackMultiSFT.data[k], stackMultiNoiseWeights.data[k], stackMultiDetStates.data[k], &CFparams );
                    if ( xlalErrno != 0 ) {
                      XLALPrintError ("%s line %d : XLALComputeFStatFreqBand() failed with xlalErrno = %d.\n\n", __func__, __LINE__, xlalErrno );
                      return(HIERARCHICALSEARCH_EXLAL);
                    }
                  }
                  else {

                    /* LALDemod method implementation to compute the F-statistic */
                    LAL_CALL( COMPUTEFSTATFREQBAND ( &status, &fstatVector.data[k], &thisPoint,
                                                     stackMultiSFT.data[k], stackMultiNoiseWeights.data[k],
                                                     stackMultiDetStates.data[k], &CFparams), &status);
                  }

                  /* Loop over coarse-grid frequency bins */
                  for (ifreq = 0; ifreq < fveclength; ifreq++) {

                    /* Get the F-statistic value ( Recall here it's *F*, not yet 2F ) */
                    if (uvar_computeLV) {
                      Fstat = multiFstatVector->F->data[ifreq];
                      fstatVector.data[k].data->data[ifreq] = Fstat;
                    }
                    else
                      Fstat = fstatVector.data[k].data->data[ifreq];

                    if ( uvar_SignalOnly )
                      {
                        /* Correct normalization in --SignalOnly case:
                         * we didn't normalize data by 1/sqrt(Tsft * 0.5 * Sh) in terms of
                         * the single-sided PSD Sh: the SignalOnly case is characterized by
                         * setting Sh->1, so we need to divide F by (0.5*Tsft)
                         */
                        Fstat *= 2.0 / Tsft;
                        Fstat += 2;		/* HERE it's *F*, but recall E[2F]:= 4 + SNR^2, so just add 2 here. */
                        fstatVector.data[k].data->data[ifreq] = Fstat; /* Reinhard: check if used later */
                        if (uvar_computeLV) {
                          for (UINT4 X = 0; X < multiFstatVector->FX->length; X++) {
                            multiFstatVector->FX->data[FX_INDEX(multiFstatVector->FX, X, ifreq)] *= 2.0 / Tsft;
                            multiFstatVector->FX->data[FX_INDEX(multiFstatVector->FX, X, ifreq)] += 2;
                          }
                        }
                      }

                    /* go to next frequency coarse-grid point */
                    freq_event = myf0 + ifreq * deltaF;

                    /* compute the global-correlation coordinate indices */
                    U1idx = ComputeU1idx ( freq_event, f1dot_eventB1, A1, u1start, u1winInv );

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
                    coarsegrid.TwoF[CG_INDEX(coarsegrid, k, ifreq)] = 2.0 * Fstat;
                    if ( uvar_computeLV ) {
                      for (UINT4 X = 0; X < coarsegrid.numDetectors; X++) {
                        INT4 detid = -1;
                        for (UINT4 Y = 0; Y < multiFstatVector->FX->length; Y++) { /* look for matching detector ID in this segment */
                          if ( strcmp( stackMultiSFT.data[k]->data[Y]->data[0].name, detectorIDs->data[X] ) == 0 )
                            detid = Y;
                        }
                        if ( detid == -1 ) /* if no match found, detector X was not present in this segment, so use 2FX=0.0 */
                          coarsegrid.TwoFX[CG_FX_INDEX(coarsegrid, X, k, ifreq)] = 0.0;
                        else /* if a match was found, get the corresponding F value and multiply by 2 */
                          coarsegrid.TwoFX[CG_FX_INDEX(coarsegrid, X, k, ifreq)] = 2.0 * multiFstatVector->FX->data[FX_INDEX(multiFstatVector->FX, detid, ifreq)];
                      } /* for X < numDetectors */
                    } /* if ( uvar_computeLV ) */

                  } /* END: Loop over coarse-grid frequency bins (ifreq) */

                  if ( multiFstatVector ) { /* free struct from XLALComputeFStatFreqBand() */
                    XLALDestroyREAL4Vector( multiFstatVector->F );
                    XLALDestroyREAL4VectorSequence( multiFstatVector->FX );
                    LALFree(multiFstatVector);
                  }

                  /* print fstat vector if required -- mostly for debugging */
                  if ( uvar_printFstat1 )
                    {
                      LAL_CALL( PrintFstatVec ( &status, &fstatVector.data[k], fpFstat1, &thisPoint, refTimeGPS, k+1), &status);
                    }

                  /* --- Holger: This is not needed in U1-only case. Sort the coarse grid in Uindex --- */
                  /* qsort(coarsegrid.list, (size_t)coarsegrid.length, sizeof(CoarseGridPoint), compareCoarseGridUindex); */

                }
                /* -------------------- END Compute F-Statistic -------------------- */

                /* timing */
                if ( uvar_outputTiming ) {
                  timeStamp2 = XLALGetTimeOfDay();
                  if (doComputeFstats) {
                    coherentTime += timeStamp2 - timeStamp1;
                  }
                  timeStamp1 = timeStamp2;
                }

                /* -------------------- Map fine grid to coarse grid -------------------- */

                /* get the frequency of this fine-grid point at mid point of segment */
                /* OLD: ifreq_fg = 0; freq_tmp = finegrid.freqmin_fg + ifreq_fg * finegrid.dfreq_fg + f1dot_tmp * timeDiffSeg; */
                freq_fg = finegrid.freqmin_fg + f1dot_fg * timeDiffSeg +
                  0.5 * f2dot_fg * timeDiffSeg * timeDiffSeg; /* first fine-grid frequency */

                /* compute the global-correlation coordinate indices */
                U1idx = ComputeU1idx ( freq_fg, f1dot_eventB1, A1, u1start, u1winInv );

                if (U1idx < 0) {
                  fprintf(stderr,"ERROR: Stepped outside the coarse grid (%d)! \n", U1idx);
                  return(HIERARCHICALSEARCH_ECG);
                }

                if (U1idx + finegrid.freqlength >= fveclength) {
                  fprintf(stderr,"ERROR: Stepped outside the coarse grid (%d:%d:%d:%d)! \n",
                          U1idx, finegrid.freqlength, U1idx + finegrid.freqlength, fveclength);
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
#ifndef EXP_NO_NUM_COUNT
                gc_hotloop( fgrid2F, cgrid2F, fgridnc, TwoFthreshold, finegrid.freqlength );
#else
                gc_hotloop_no_nc ( fgrid2F, cgrid2F, finegrid.freqlength );
#endif
                if ( uvar_computeLV ) {
                  for (UINT4 X = 0; X < finegrid.numDetectors; X++) {
                    REAL4 * cgrid2FX = coarsegrid.TwoFX + CG_FX_INDEX(coarsegrid, X, k, U1idx);
                    REAL4 * fgrid2FX = finegrid.sumTwoFX + FG_FX_INDEX(finegrid, X, 0);
                    gc_hotloop_no_nc( fgrid2FX, cgrid2FX, finegrid.freqlength );
                  }
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
                if ( uvar_computeLV ) {
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
#endif // GC_SSE2_OPT

                /* timing */
                if ( uvar_outputTiming ) {
                  timeStamp2 = XLALGetTimeOfDay();
                  incoherentTime += timeStamp2 - timeStamp1;
                }

              } /* end: ------------- MAIN LOOP over Segments --------------------*/

              /* ############################################################### */

              if( uvar_semiCohToplist ) {
                /* this is necessary here, because UpdateSemiCohToplists() might set
                   a checkpoint that needs some information from here */
                LAL_CALL( UpdateSemiCohToplists (&status, semiCohToplist, semiCohToplist2, &finegrid, f1dot_fg, f2dot_fg, &usefulParams, NSegmentsInv, NSegmentsInvX ), &status);
              }

            } /* for( if1dot_fg = 0; if1dot_fg < nf1dots_fg; if1dot_fg++ ) */
          } /* for( if2dot_fg = 0; if2dot_fg < nf2dots_fg; if2dot_fg++ ) */
          /* ---------- END walk through fine grid fdots --------------- */

          if2dot++;  /* Increment if2dot counter */

        } /* ########## End of loop over coarse-grid f2dot values (if2dot) ########## */
        ifdot++;  /* Increment ifdot counter BEFORE SET_GCT_CHECKPOINT */

        SHOW_PROGRESS(dopplerpos.Alpha, dopplerpos.Delta,
                      skyGridCounter * nf1dot + ifdot,
                      thisScan.numSkyGridPoints * nf1dot, uvar_Freq, uvar_FreqBand);

        SET_GCT_CHECKPOINT (uvar_fnameChkPoint, semiCohToplist, semiCohToplist2, skyGridCounter*nf1dot+ifdot, TRUE);

      } /* ########## End of loop over coarse-grid f1dot values (ifdot) ########## */

      /* continue forward till the end if uvar_skyPointIndex is set
         ---This probably doesn't make sense when checkpointing is turned on */
      if ( LALUserVarWasSet(&uvar_skyPointIndex) ) {
        while(thisScan.state != STATE_FINISHED) {
          skyGridCounter++;
          XLALNextDopplerSkyPos(&dopplerpos, &thisScan);
        }
      }
      else {
        skyGridCounter++;

        /* this is necessary here, because the checkpoint needs some information from here */
        SHOW_PROGRESS(dopplerpos.Alpha, dopplerpos.Delta,
                      skyGridCounter * nf1dot,
                      thisScan.numSkyGridPoints * nf1dot, uvar_Freq, uvar_FreqBand);

        XLALNextDopplerSkyPos( &dopplerpos, &thisScan );
      }

    } /* ######## End of while loop over 1st stage SKY coarse-grid points ############ */
  /*---------------------------------------------------------------------------------*/

  /* now that we have the final toplist, translate all pulsar parameters to correct reftime */
  xlalErrno = 0;
  XLALExtrapolateToplistPulsarSpins ( semiCohToplist, usefulParams.spinRange_refTime.refTime, finegrid.refTime );
  if ( semiCohToplist2 )	// handle (optional) second toplist
    XLALExtrapolateToplistPulsarSpins ( semiCohToplist2, usefulParams.spinRange_refTime.refTime, finegrid.refTime );
  if ( xlalErrno != 0 ) {
    XLALPrintError ("%s line %d : XLALExtrapolateToplistPulsarSpins() failed with xlalErrno = %d.\n\n", __func__, __LINE__, xlalErrno );
    return(HIERARCHICALSEARCH_EXLAL);
  }

  /* timing */
  if ( uvar_outputTiming )
    timeEnd = XLALGetTimeOfDay();

  LogPrintf( LOG_NORMAL, "Finished main analysis.\n");

  /* Also compute F, FX (for line veto statistics) for all candidates in final toplist */
  if ( uvar_recalcToplistStats ) {

    LogPrintf( LOG_NORMAL, "Recalculating statistics for the final toplist...\n");

    /* timing */
    if ( uvar_outputTiming )
      timeStamp1 = XLALGetTimeOfDay();

    /* need pre-sorted toplist to have right segment-Fstat file numbers (do not rely on this feature, could be messed up by precision issues in sorting!) */
    if ( uvar_outputSingleSegStats )
      {
        sort_gctFStat_toplist(semiCohToplist);
        if ( semiCohToplist2 )
          sort_gctFStat_toplist(semiCohToplist2);
      }

    XLAL_CHECK ( XLAL_SUCCESS == XLALComputeExtraStatsForToplist ( semiCohToplist, "GCTtop", &stackMultiSFT, &stackMultiNoiseWeights, &stackMultiDetStates, &CFparams, refTimeGPS, uvar_SignalOnly, uvar_outputSingleSegStats ),
                 HIERARCHICALSEARCH_EXLAL, "XLALComputeExtraStatsForToplist() failed with xlalErrno = %d.\n\n", xlalErrno
                 );
    // also recalc optional 2nd toplist if present
    if ( semiCohToplist2 )
      XLAL_CHECK ( XLAL_SUCCESS == XLALComputeExtraStatsForToplist ( semiCohToplist2, "GCTtop", &stackMultiSFT, &stackMultiNoiseWeights, &stackMultiDetStates, &CFparams, refTimeGPS, uvar_SignalOnly, uvar_outputSingleSegStats ),
                   HIERARCHICALSEARCH_EXLAL, "XLALComputeExtraStatsForToplist() failed for 2nd toplist with xlalErrno = %d.\n\n", xlalErrno
                   );

    /* timing */
    if ( uvar_outputTiming ) {
      timeStamp2 = XLALGetTimeOfDay();
      vetoTime = timeStamp2 - timeStamp1;
    }

    LogPrintf( LOG_NORMAL, "Finished recalculating toplist statistics.\n");

  }

  if ( uvar_outputTiming )
    {
      FILE *timing_fp;
      if ( ( timing_fp = fopen ( uvar_outputTiming, "ab" )) == NULL ) {
        XLALPrintError ("%s: failed to open timing-file '%s' for appending.\n", __func__, uvar_outputTiming );
        return HIERARCHICALSEARCH_EFILE;
      }
      REAL8 tauWU = timeEnd - timeStart;

      /* compute fundamental timing-model constants:
       * 'tauF0' = Fstat time per template per SFT,
       * 'tauS0' = time to add one per-segment F-stat value per fine-grid point
       */
      REAL8 Ncoarse = thisScan.numSkyGridPoints * nf1dot * nf2dot * ( binsFstatSearch + 2 * semiCohPar.extraBinsFstat);	// includes GCSideband bins
      REAL8 NSB     = thisScan.numSkyGridPoints * nf1dot * nf2dot * 2 * semiCohPar.extraBinsFstat;	// pure coarse GCSideband template count
      REAL8 Nfine   = thisScan.numSkyGridPoints * binsFstatSearch * nf1dots_fg * nf2dots_fg;	// doesn't include F-stat sideband bins
      REAL8 c0co    = coherentTime / ( Ncoarse * nSFTs );
      // Note: we use (total-FstatTime) instead of incoherentTime to ensure accurate prediction power
      // whatever extra time isn't captured by incoherentTime+coherentTime also needs to be accounted as 'incoherent time'
      // in our model...
      REAL8 c0ic   = (tauWU - coherentTime) / ( Nfine * nStacks );

      fprintf ( timing_fp, "%10.3g %10.3g %10.3g %7d %7d    %9.3g %9.3g %9.3g %9.3g    %10.3g %10.3g\n",
                Ncoarse, NSB, Nfine, nSFTs, nStacks, tauWU, coherentTime, incoherentTime, vetoTime, c0co, c0ic );
      fclose ( timing_fp );
    } // if uvar_outputTiming

  LogPrintf ( LOG_DEBUG, "Writing output ... ");
  XLAL_CHECK ( write_hfs_oputput(uvar_fnameout, semiCohToplist) != -1, XLAL_EFAILED, "write_hfs_oputput('%s', toplist) failed.!\n", uvar_fnameout );
  // output optional second toplist, if it exists, into "<uvar_fnameout>-LV"
  if ( semiCohToplist2 )
    {
      LogPrintf ( LOG_DEBUG, "toplist2 ... ");
      UINT4 newlen = strlen(uvar_fnameout) + 10;
      CHAR *fname2;
      XLAL_CHECK ( (fname2 = XLALCalloc ( 1, newlen )) != NULL, XLAL_ENOMEM, "Failed to XLALCalloc(1, %d)\n\n", newlen );
      sprintf ( fname2, "%s-LV", uvar_fnameout );
      XLAL_CHECK ( write_hfs_oputput ( fname2, semiCohToplist2) != -1, XLAL_EFAILED, "write_hfs_oputput('%s', toplist2) failed for 2nd toplist!\n", fname2 );
      XLALFree ( fname2 );
    }

  LogPrintfVerbatim ( LOG_DEBUG, "done.\n");

  // in BOINC App the checkpoint is left behind to be cleaned up by the Core Client
#ifndef EAH_BOINC
  clear_gct_checkpoint (uvar_fnameChkPoint);
#endif

  /*------------ free all remaining memory -----------*/

  if ( uvar_printCand1 ) {
    LALFree( fnameSemiCohCand );
  }

  if ( uvar_printFstat1 ) {
    fclose(fpFstat1);
    LALFree( fnameFstatVec1 );
  }

  /* free first stage memory */
  for ( k = 0; k < nStacks; k++) {
    LAL_CALL( LALDestroyMultiSFTVector ( &status, stackMultiSFT.data + k), &status);
    LAL_CALL( LALDestroyMultiNoiseWeights ( &status, stackMultiNoiseWeights.data + k), &status);
    XLALDestroyMultiDetectorStateSeries ( stackMultiDetStates.data[k] );
  }

  if (stackMultiSFT.data)
    LALFree(stackMultiSFT.data);
  if (stackMultiNoiseWeights.data)
    LALFree(stackMultiNoiseWeights.data);
  if (stackMultiDetStates.data)
    LALFree(stackMultiDetStates.data);

  XLALDestroyTimestampVector(startTstack);
  XLALDestroyTimestampVector(midTstack);
  XLALDestroyTimestampVector(endTstack);

  /* free Fstat vectors  */
  for(k = 0; k < nStacks; k++)
    if (fstatVector.data[k].data) {
      if (fstatVector.data[k].data->data)
        LALFree(fstatVector.data[k].data->data);
      LALFree(fstatVector.data[k].data);
    }
  LALFree(fstatVector.data);

  /* if resampling is used then free buffer */
  if ( uvar_useResamp ) {
    for (k=0;k<resampbuffers.length;k++) {
      XLALEmptyComputeFBuffer_RS( resampbuffers.data[k] );
      XLALFree(resampbuffers.data[k]);
    }
    XLALFree(resampbuffers.data);
  }

  /* free Vel/Pos/Acc vectors and ephemeris */
  XLALDestroyREAL8VectorSequence( posStack );
  XLALDestroyREAL8VectorSequence( velStack );
  XLALDestroyREAL8VectorSequence( accStack );
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);

  /* free dopplerscan stuff */
  LAL_CALL ( FreeDopplerSkyScan(&status, &thisScan), &status);
  if ( scanInit.skyRegionString )
    LALFree ( scanInit.skyRegionString );

  XLALDestroyStringVector ( detectorIDs );
  if (NSegmentsInvX)
    ALFree(NSegmentsInvX);

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

  if (coarsegrid.TwoF) {
    LALFree(coarsegrid.TwoF);
  }
  if (coarsegrid.TwoFX) {
    LALFree(coarsegrid.TwoFX);
  }
  if (coarsegrid.Uindex) {
    LALFree(coarsegrid.Uindex);
  }

  free_gctFStat_toplist ( &semiCohToplist );
  if ( semiCohToplist2 ) free_gctFStat_toplist ( &semiCohToplist2 );

  XLALDestroyREAL8Vector ( usefulParams.LVloglX );

  XLALDestroyExpLUT(); /* lookup table for fast exponential function, used in computeLV case */
  XLALDestroyLogLUT(); /* lookup table for fast logarithm function, used in computeLV case */

  LAL_CALL (LALDestroyUserVars(&status), &status);

  XLALFree ( VCSInfoString );

  LALCheckMemoryLeaks();

  return HIERARCHICALSEARCH_ENORM;
} /* main */







/** Set up stacks, read SFTs, calculate SFT noise weights and calculate
    detector-state */
void SetUpSFTs( LALStatus *status,			/**< pointer to LALStatus structure */
                MultiSFTVectorSequence *stackMultiSFT, /**< output multi sft vector for each stack */
                MultiNoiseWeightsSequence *stackMultiNoiseWeights, /**< output multi noise weights for each stack */
                MultiDetectorStateSeriesSequence *stackMultiDetStates, /**< output multi detector states for each stack */
                UsefulStageVariables *in, /**< input params */
                BOOLEAN useWholeSFTs,	/**< special switch: load all given frequency bins from SFTs */
                REAL8 UNUSED mismatch1		/**< 'mismatch1' user-input needed here internally ... */
                )
{
  SFTCatalog *catalog = NULL;
  static SFTConstraints constraints;
  REAL8 timebase, tObs, deltaFsft;
  UINT4 k,numSFT;
  LIGOTimeGPS tStartGPS, tEndGPS, refTimeGPS, tMidGPS, midTstackGPS, startTstackGPS, endTstackGPS;
  SFTCatalogSequence catalogSeq;
  REAL8 midTseg,startTseg,endTseg;

  REAL8 doppWings, freqmin, freqmax;
  REAL8 startTime_freqLo, startTime_freqHi;
  REAL8 endTime_freqLo, endTime_freqHi;
  REAL8 freqLo, freqHi;
  INT4 extraBins;

  INT4 sft_check_result = 0;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* get sft catalog */
  constraints.startTime = &(in->minStartTimeGPS);
  constraints.endTime = &(in->maxEndTimeGPS);
  TRY( LALSFTdataFind( status->statusPtr, &catalog, in->sftbasename, &constraints), status);

  /* check CRC sums of SFTs */
  TRY ( LALCheckSFTCatalog ( status->statusPtr, &sft_check_result, catalog ), status );
  if (sft_check_result) {
    LogPrintf(LOG_CRITICAL,"SFT validity check failed (%d)\n", sft_check_result);
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
  TRY( LALExtrapolatePulsarSpinRange( status->statusPtr, &in->spinRange_startTime, tStartGPS, &in->spinRange_refTime), status);
  TRY( LALExtrapolatePulsarSpinRange( status->statusPtr, &in->spinRange_endTime, tEndGPS, &in->spinRange_refTime), status);
  TRY( LALExtrapolatePulsarSpinRange( status->statusPtr, &in->spinRange_midTime, tMidGPS, &in->spinRange_refTime), status);

  /* set Fstat spindown resolution (coarse grid) */
  in->df1dot = HSMIN(in->df1dot, in->spinRange_midTime.fkdotBand[1]);

  /* set Fstat 2nd spindown resolution (coarse grid) */
  in->df2dot = HSMIN(in->df2dot, in->spinRange_midTime.fkdotBand[2]);

  /* calculate number of bins for Fstat overhead due to residual spin-down */
  in->extraBinsFstat = (UINT4)( 0.25*(in->tObs*in->df1dot + in->tObs*in->tObs*in->df2dot)/in->dFreqStack + 1e-6) + 1;

  /* set wings of sfts to be read */
  /* the wings must be enough for the Doppler shift and extra bins
     for the running median block size and Dterms for Fstat calculation.
     In addition, it must also include wings for the spindown correcting
     for the reference time  */
  /* calculate Doppler wings at the highest frequency */
  startTime_freqLo = in->spinRange_startTime.fkdot[0]; /* lowest search freq at start time */
  startTime_freqHi = startTime_freqLo + in->spinRange_startTime.fkdotBand[0]; /* highest search freq. at start time*/
  endTime_freqLo = in->spinRange_endTime.fkdot[0];
  endTime_freqHi = endTime_freqLo + in->spinRange_endTime.fkdotBand[0];

  freqLo = HSMIN ( startTime_freqLo, endTime_freqLo );
  freqHi = HSMAX ( startTime_freqHi, endTime_freqHi );
  doppWings = freqHi * in->dopplerMax;    /* maximum Doppler wing -- probably larger than it has to be */
  extraBins = HSMAX ( in->blocksRngMed/2 + 1, in->Dterms );

  if (useWholeSFTs) {
    freqmin = freqmax = -1;
  }
  else {
    freqmin = freqLo - doppWings - extraBins * deltaFsft - in->extraBinsFstat * in->dFreqStack;
    freqmax = freqHi + doppWings + extraBins * deltaFsft + in->extraBinsFstat * in->dFreqStack;
  }

  /* ----- finally memory for segments of multi sfts ----- */
  stackMultiSFT->length = in->nStacks;
  stackMultiSFT->data = (MultiSFTVector **)LALCalloc(1, in->nStacks * sizeof(MultiSFTVector *));
  if ( stackMultiSFT->data == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }

  stackMultiDetStates->length = in->nStacks;
  stackMultiDetStates->data = (MultiDetectorStateSeries **)LALCalloc(1, in->nStacks * sizeof(MultiDetectorStateSeries *));
  if ( stackMultiDetStates->data == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }

  stackMultiNoiseWeights->length = in->nStacks;
  if ( in->SignalOnly )  {
    stackMultiNoiseWeights->data = (MultiNoiseWeights **)LALMalloc(in->nStacks * sizeof(MultiNoiseWeights *));
    if ( stackMultiNoiseWeights->data == NULL ) {
      ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
    }
  }
  else {
    stackMultiNoiseWeights->data = (MultiNoiseWeights **)LALCalloc(1, in->nStacks * sizeof(MultiNoiseWeights *));
    if ( stackMultiNoiseWeights->data == NULL ) {
      ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
    }
  }


  /* loop over segments and read sfts */
  for (k = 0; k < in->nStacks; k++) {

    /* ----- load the multi-IFO SFT-vectors ----- */
    stackMultiSFT->data[k] = XLALLoadMultiSFTs( catalogSeq.data + k, freqmin, freqmax );
    if ( stackMultiSFT->data[k] == NULL ) {
      ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
    }

    /* ----- obtain the (multi-IFO) 'detector-state series' for all SFTs ----- */
    TRY ( LALGetMultiDetectorStates ( status->statusPtr, stackMultiDetStates->data + k,
                                      stackMultiSFT->data[k], in->edat ), status );

    /* ----- normalize sfts and compute noise weights ----- */
    if ( in->SignalOnly )  {
      stackMultiNoiseWeights->data[k] = NULL;
    }
    else {
      MultiPSDVector *psd = NULL;
      TRY( LALNormalizeMultiSFTVect ( status->statusPtr, &psd, stackMultiSFT->data[k],
                                      in->blocksRngMed ), status );
      TRY( LALComputeMultiNoiseWeights  ( status->statusPtr, stackMultiNoiseWeights->data + k,
                                          psd, in->blocksRngMed, 0 ), status );
      TRY ( LALDestroyMultiPSDVector ( status->statusPtr, &psd ), status );
    } /* if ( in->SignalOnly )  */


  } /* loop over k */



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
  TRY( LALDestroySFTCatalog( status->statusPtr, &catalog ), status);

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






/** \brief Breaks up input sft catalog into specified number of stacks

    Loops over elements of the catalog, assigns a bin index and
    allocates memory to the output catalog sequence appropriately.  If
    there are long gaps in the data, then some of the catalogs in the
    output catalog sequence may be of zero length.
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






/** Get SemiCoh candidates into toplist(s)
 * This function allows for inserting candidates into up to 2 toplists at once, which might be sorted differently!
 */
void UpdateSemiCohToplists ( LALStatus *status,
                             toplist_t *list1,
                             toplist_t *list2,	//< optional (can be NULL): insert candidate into this 2nd toplist as well
                             FineGrid *in,
                             REAL8 f1dot_fg,
                             REAL8 f2dot_fg,
                             UsefulStageVariables *usefulparams,
                             REAL4 NSegmentsInv,
                             REAL4 *NSegmentsInvX
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
    line.nc = in->nc[ifreq_fg];
    line.sumTwoF = in->sumTwoF[ifreq_fg]; /* here it's still the summed 2F value over segments, not the average */
    line.numDetectors = in->numDetectors;
    for (UINT4 X = 0; X < GCTTOP_MAX_IFOS; X++) { /* initialise single-IFO F-stat arrays to zero */
      line.sumTwoFX[X] = 0.0;
      line.sumTwoFXrecalc[X] = 0.0;
    }
    line.sumTwoFrecalc = -1.0; /* initialise this to -1.0, so that it only gets written out by print_gctFStatline_to_str if later overwritten in recalcToplistStats step */

    if ( in->sumTwoFX ) { /* if we already have FX values from the main loop, insert these, and calculate LV-stat here */
      for (UINT4 X = 0; X < in->numDetectors; X++)
        line.sumTwoFX[X] = in->sumTwoFX[FG_FX_INDEX(*in, X, ifreq_fg)]; /* here it's still the summed 2F value over segments, not the average */
      xlalErrno = 0;

      REAL8 *loglX = NULL;
      if ( usefulparams->LVloglX )
        loglX = usefulparams->LVloglX->data;

      line.LV = XLALComputeLineVetoArray ( line.sumTwoF, line.numDetectors, line.sumTwoFX, usefulparams->LVlogRhoTerm, loglX, usefulparams->LVuseAllTerms );
      line.LV *= NSegmentsInv; /* normalize by number of segments */

      if ( xlalErrno != 0 ) {
        XLALPrintError ("%s line %d : XLALComputeLineVeto() failed with xlalErrno = %d.\n\n", __func__, __LINE__, xlalErrno );
        ABORT ( status, HIERARCHICALSEARCH_EXLAL, HIERARCHICALSEARCH_MSGEXLAL );
      }
      if ( line.LV < -LAL_REAL4_MAX*0.1 )
        line.LV = -LAL_REAL4_MAX*0.1; /* avoid minimum value, needed for output checking in print_gctFStatline_to_str() */
    }
    else
      line.LV = -LAL_REAL4_MAX; /* in non-LV case, block field with minimal value, needed for output checking in print_gctFStatline_to_str() */

    /* take F-stat averages over segments */
    line.sumTwoF *= NSegmentsInv; /* average multi-2F by full number of segments */
    if ( in->sumTwoFX ) {
      for (UINT4 X = 0; X < in->numDetectors; X++) {
        line.sumTwoFX[X] *= NSegmentsInvX[X]; /* average single-2F by per-IFO number of segments */
      }
    }

    insert_into_gctFStat_toplist( list1, line);
    if ( list2 )	// also insert candidate into (optional) second toplist
      insert_into_gctFStat_toplist( list2, line);

  }

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* UpdateSemiCohToplists() */






/** Print Fstat vectors */
void PrintFstatVec (LALStatus *status,
                    REAL4FrequencySeries *in,
                    FILE                 *fp,
                    PulsarDopplerParams  *thisPoint,
                    LIGOTimeGPS          refTime,
                    INT4                 stackIndex)
{
  INT4 length, k;
  REAL8 f0, deltaF, alpha, delta;
  PulsarSpins fkdot;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  INIT_MEM(fkdot);

  fprintf(fp, "%% Fstat values from stack %d (reftime -- %d %d)\n", stackIndex, refTime.gpsSeconds, refTime.gpsNanoSeconds);

  alpha = thisPoint->Alpha;
  delta = thisPoint->Delta;
  fkdot[1] = thisPoint->fkdot[1];
  f0 = thisPoint->fkdot[0];

  length = in->data->length;
  deltaF = in->deltaF;

  for (k=0; k<length; k++)
    {
      fkdot[0] = f0 + k*deltaF;

      /* propagate fkdot back to reference-time  */
      TRY ( LALExtrapolatePulsarSpins (status->statusPtr, fkdot, refTime, fkdot, thisPoint->refTime ), status );

      fprintf(fp, "%d %.13g %.12g %.12g %.13g %.13g %.6g\n",
              stackIndex, fkdot[0], alpha, delta, fkdot[1], fkdot[2], 2*in->data->data[k]);
    }

  fprintf(fp, "\n");

  DETATCHSTATUSPTR (status);
  RETURN(status);

}















/** Calculate Earth orbital position, velocity and acceleration
    at midpoint of each segment */
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
                                 REAL8 f1dot_eventB1,
                                 REAL8 A1,
                                 REAL8 U1start,
                                 REAL8 U1winInv)
{
  /* compute the index of global-correlation coordinate U1, Eq. (1) */
  return (((freq_event * A1 + f1dot_eventB1) - U1start) * U1winInv) + 0.5;

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

/** Set up 'segmented' SFT-catalogs for given list of segments and a total SFT-catalog.
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
    XLALPrintError ("%s: XLALCalloc(%d) failed.\n", __func__, sizeof(*stacks) );
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
          LIGOTimeGPS gpsStart = catalog->data[iSFT0].header.epoch;
          int cmp = XLALGPSInSeg ( &gpsStart, thisSeg );

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
          LIGOTimeGPS gpsEnd = catalog->data[iSFT1].header.epoch;
          int cmp = XLALGPSInSeg ( &gpsEnd, thisSeg );

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




/** XLAL function to (multi-IFO) F-statistic over a vector of number of frequency bins,
    returned in (*fstatSeries)->F.

    if params->returnSingleF==true: also returns per-IFO F-stat values over frequency bins in (*fstatSeries)->FX.

    Note: Contrary to ComputeFStatFreqBand(), the output (*fstatSeries) can be allocated
    before this function is called, or passed as a NULL pointer, in which case it will be allocated.

    This allows one to re-use the output structure vectors without unneccessary alloc/free's,
    and simplifies usage of this function for the called.

    Note2: the start frequency, step size in frequency and the number of frequency bins to be
    computed are *always* read from input 'doppler'.

    Note3: This function is currently simply a wrapper for ComputeFstat(), while future implementations
    will also include resampling.
*/
int XLALComputeFStatFreqBand ( MultiFstatFrequencySeries **fstatSeries,	/**< [out] Combined vectors of multi- and single-IFO Fstat values */
                               const PulsarDopplerParams *doppler,		/**< parameter-space point to compute F for (and freq band info) */
                               const MultiSFTVector *multiSFTs, 		/**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
                               const MultiNoiseWeights *multiWeights,		/**< noise-weights of all SFTs */
                               const MultiDetectorStateSeries *multiDetStates,/**< 'trajectories' of the different IFOs */
                               const ComputeFParams *params			/**< addition computational params */
                               )
{

  /* check input parameters */
  if ( !fstatSeries )
    XLAL_ERROR ( XLAL_EFAULT, "\nNULL input pointer 'fstatSeries'\n" );
  if ( !doppler )
    XLAL_ERROR ( XLAL_EFAULT, "\nInput dopplerParams pointer is NULL !\n\n");
  if ( doppler->orbit )
    XLAL_ERROR ( XLAL_EDOM, "\ndoppler->orbit != NULL, but binary parameters currently not supported by this function!\n\n");
  if ( !multiSFTs )
    XLAL_ERROR ( XLAL_EFAULT, "\nInput multiSFTs pointer is NULL !\n\n");
  if ( !multiDetStates )
    XLAL_ERROR ( XLAL_EFAULT, "\nInput multiDetStates pointer is NULL !\n\n");
  if ( multiDetStates->length != multiSFTs->length )
    XLAL_ERROR ( XLAL_EBADLEN, "\nInput vector lengths do not match (len(multiDetStates)=%d and len(multiSFTs)=%d) !\n\n", multiDetStates->length, multiSFTs->length);
  if ( !params )
    XLAL_ERROR ( XLAL_EFAULT, "\nInput CFParams pointer is NULL !\n\n");
  if ( params->returnAtoms )
    XLAL_ERROR ( XLAL_EINVAL, "\nUsing the option 'returnAtoms' is not supported in this function!\n\n");

  /* some useful shortcuts */
  UINT4 numBins      = doppler->numFreqBins;
  UINT4 numDetectors = multiSFTs->length;
  MultiFstatFrequencySeries * retFstatSeries = (*fstatSeries);       /* build up new return structure either from scratch or by reallocating input */

  /* ---------- check if output structure retFstatSeries exists and (re-)alloc if necessary ---------- */
  if ( retFstatSeries == NULL  && (retFstatSeries = XLALCalloc( 1, sizeof(*retFstatSeries))) == NULL )
    XLAL_ERROR ( XLAL_ENOMEM, "\nFailed to allocate memory for MultiFstatFrequencySeries !\n\n" );

  /* (re)set output structure meta info (search parameters and frequency band) from input values */
  retFstatSeries->doppler = (*doppler);		// struct-copy  FIXME: this would break for binary-NS searches

  /* check and (re)alloc F vector */
  if ( retFstatSeries->F == NULL && (retFstatSeries->F = XLALCalloc ( 1, sizeof(*retFstatSeries->F))) == NULL )
    XLAL_ERROR ( XLAL_ENOMEM, "\nFailed to allocate memory for MultiFstatFrequencySeries->F !\n\n" );
  retFstatSeries->F->length = numBins;
  if ( (retFstatSeries->F->data = XLALRealloc(retFstatSeries->F->data, numBins * sizeof(*retFstatSeries->F->data))) == NULL )
    XLAL_ERROR ( XLAL_ENOMEM, "\nFailed to re-allocate %d elements for MultiFstatFrequencySeries->F->data !\n\n", numBins );

  /* if FX requested, check and (re)alloc FX vector sequence */
  if ( params->returnSingleF )
    {
      if ( retFstatSeries->FX == NULL && (retFstatSeries->FX = XLALCalloc ( 1, sizeof(*retFstatSeries->FX))) == NULL )
        XLAL_ERROR ( XLAL_ENOMEM, "\nFailed to allocate memory for MultiFstatFrequencySeries->FX !\n\n" );
      retFstatSeries->FX->length = numDetectors;
      retFstatSeries->FX->vectorLength = numBins;
      if ( (retFstatSeries->FX->data = XLALRealloc( retFstatSeries->FX->data, numDetectors * numBins * sizeof(*retFstatSeries->FX->data))) == NULL )
        XLAL_ERROR ( XLAL_ENOMEM, "\nFailed to re-allocate %d elements for MultiFstatFrequencySeries->FX->data !\n\n", numBins * numDetectors );
    }
  else
    { /* if no FX return requested, destroy FX field if it exists */
      if ( retFstatSeries->FX != NULL )
        XLALDestroyREAL4VectorSequence ( retFstatSeries->FX );
    }
  /* ---------- END: memory-handling of output structure retFstatSeries ----------*/

  /* copy values from 'doppler' to local variable 'thisPoint' */
  PulsarDopplerParams thisPoint = (*doppler);	// struct copy
  REAL8 dFreq   = thisPoint.dFreq;
  REAL8 fStart  = thisPoint.fkdot[0];

  ComputeFBuffer cfBuffer = empty_ComputeFBuffer;
  LALStatus fakeStatus = blank_status;	    /* fake LAL status structure, needed as long as ComputeFStat is LAL function and not XLAL */

  /* loop over frequency values and fill up values in fstatSeries */
  for ( UINT4 k = 0; k < numBins; k++) {
    Fcomponents Fstat;

    thisPoint.fkdot[0] = fStart + k*dFreq;

    COMPUTEFSTAT ( &fakeStatus, &Fstat, &thisPoint, multiSFTs, multiWeights, multiDetStates, params, &cfBuffer );
    if ( fakeStatus.statusCode )
      XLAL_ERROR (XLAL_EFUNC, "\nFailure in LAL function ComputeFStat(). statusCode=%d\n\n", fakeStatus.statusCode);

    retFstatSeries->F->data[k] = Fstat.F;

    if ( params->returnSingleF )
      for ( UINT4 X=0; X < numDetectors; X ++)
        retFstatSeries->FX->data[FX_INDEX(retFstatSeries->FX, X, k) ] = Fstat.FX[X]; /* fstatSeries->FX->data is ordered as (det1bin1,det1bin2,..,det1binN,det2bin1,...detMbinN) */

  } /* for k < numBins */

  XLALEmptyComputeFBuffer ( &cfBuffer );

  /* return result */
  (*fstatSeries) = retFstatSeries;

  return XLAL_SUCCESS;

} /* XLALComputeFStatFreqBand() */



/** XLAL function to extrapolate the pulsar spin parameters of all toplist candidates
 * from reftime of the input toplist ('inRefTime') to a user-specified output reftime 'outRefTime'
 */
int XLALExtrapolateToplistPulsarSpins ( toplist_t *list,              /**< [out/in] toplist with GCTtopOutputEntry items, 'Freq,F1dot,F2dot' fields will be overwritten  */
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

  PulsarSpins fkdot;
  INIT_MEM ( fkdot );

  UINT4 numElements = list->elems;
  for (UINT4 j = 0; j < numElements; j++ ) /* loop over toplist */
    {
      /* get fkdot of each candidate */
      GCTtopOutputEntry *elem = toplist_elem ( list, j );
      fkdot[0] = elem->Freq;
      fkdot[1] = elem->F1dot;
      fkdot[2] = elem->F2dot;
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
    }

  return XLAL_SUCCESS;

} /* XLALExtrapolateToplistPulsarSpins() */
