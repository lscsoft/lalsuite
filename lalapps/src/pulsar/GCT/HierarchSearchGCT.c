/*
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
 * \brief Hierarchical semicoherent CW search code based on F-Statistic,
 *  exploiting global-correlation coordinates (Phys.Rev.Lett. 103, 181102, 2009)
 *
 *********************************************************************************/

/* ---------- Includes -------------------- */
#include "HierarchSearchGCT.h"

RCSID( "$Id$");

/* ---------- Defines -------------------- */
/* #define OUTPUT_TIMING 1 */
/* #define DIAGNOSISMODE 1 */

#define TRUE (1==1)
#define FALSE (1==0)

/* Hooks for Einstein@Home / BOINC
   These are defined to do nothing special in the standalone case
   and will be set in boinc_extras.h if EAH_BOINC is set
 */
#ifdef EAH_BOINC
#include "hs_boinc_extras.h"
#define COMPUTEFSTATFREQBAND_RS ComputeFStatFreqBand_RS
#else
#define GET_CHECKPOINT(toplist,total,countp,outputname,cptname) if(read_hfs_checkpoint("checkpoint.cpt", semiCohToplist, &count)) count=0
#define SET_CHECKPOINT write_hfs_checkpoint("checkpoint.cpt",semiCohToplist,skyGridCounter*nf1dot+ifdot,1)
#define SHOW_PROGRESS(rac,dec,skyGridCounter,tpl_total,freq,fband)
#define MAIN  main
#define FOPEN fopen
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
#else
#define COMPUTEFSTATFREQBAND ComputeFStatFreqBand
#endif
#define COMPUTEFSTATFREQBAND_RS ComputeFStatFreqBand_RS
char**global_argv;
int global_argc;
#endif /* EAH_BOINC */

#define EARTHEPHEMERIS  "earth05-09.dat"
#define SUNEPHEMERIS 	"sun05-09.dat"
#define BLOCKSRNGMED 	101 	/**< Default running median window size */
#define FSTART 		100.0	/**< Default Start search frequency */
#define FBAND           0.01  /**< Default search band */
#define FDOT            0.0	  /**< Default value of first spindown */
#define DFDOT           0.0	  /**< Default range of first spindown parameter */
#define SKYREGION       "allsky" /**< default sky region to search over -- just a single point*/
#define DTERMS          16    /**< Default number of dirichlet kernel terms for calculating Fstat */

/**< Default number of dirichlet kernel terms for calculating Fstat */
#define MISMATCH        0.3 	  /**< Default for metric grid maximal mismatch value */
#define DALPHA          0.001 	/**< Default resolution for isotropic or flat grids */
#define DDELTA          0.001 	/**< Default resolution for isotropic or flat grids */
#define FSTATTHRESHOLD 	2.6	/**< Default threshold on Fstatistic for peak selection */
#define NCAND1          10 	/**< Default number of candidates to be followed up from first stage */
#define FNAMEOUT        "./HS_GCT.out"  /**< Default output file basename */
#ifndef LAL_INT4_MAX
#define LAL_INT4_MAX 	2147483647
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
  REAL8 tStack;          /**< duration of stacks */
  UINT4 nStacks;         /**< number of stacks */
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
} UsefulStageVariables;


/* ------------------------ Functions -------------------------------- */
void SetUpSFTs( LALStatus *status, MultiSFTVectorSequence *stackMultiSFT,
               MultiNoiseWeightsSequence *stackMultiNoiseWeights,
               MultiDetectorStateSeriesSequence *stackMultiDetStates, UsefulStageVariables *in );
void PrintFstatVec( LALStatus *status, REAL4FrequencySeries *in, FILE *fp, PulsarDopplerParams *thisPoint,
                   LIGOTimeGPS refTime, INT4 stackIndex);
void PrintCatalogInfo( LALStatus *status, const SFTCatalog *catalog, FILE *fp );
void PrintStackInfo( LALStatus *status, const SFTCatalogSequence *catalogSeq, FILE *fp );
void UpdateSemiCohToplist( LALStatus *status, toplist_t *list, FineGrid *in, UsefulStageVariables *usefulparams );
void GetSegsPosVelAccEarthOrb( LALStatus *status, REAL8VectorSequence **posSeg,
                              REAL8VectorSequence **velSeg, REAL8VectorSequence **accSeg,
                              UsefulStageVariables *usefulparams );
static inline INT4 ComputeU1idx( REAL8 freq_event, REAL8 f1dot_eventB1, REAL8 A1, REAL8 U1start, REAL8 U1winInv );
void ComputeU2idx( REAL8 freq_event, REAL8 f1dot_event, REAL8 A2, REAL8 B2, REAL8 U2start, REAL8 U2winInv,
                  INT4 *U2idx);
int compareCoarseGridUindex( const void *a, const void *b );
int compareFineGridNC( const void *a,const void *b );
int compareFineGridsumTwoF( const void *a,const void *b );

/* ---------- Global variables -------------------- */
LALStatus *global_status; /* a global pointer to MAIN()s head of the LALStatus structure */
extern int lalDebugLevel;

#ifdef OUTPUT_TIMING
time_t clock0;
UINT4 nSFTs;
#endif


/* ###################################  MAIN  ################################### */

int MAIN( int argc, char *argv[]) {

  LALStatus status = blank_status;

  /* temp loop variables: generally k loops over segments and j over SFTs in a stack */
  INT4 j;
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
  LIGOTimeGPS tStartGPS = empty_LIGOTimeGPS;
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
  UINT4 binsFstat1, binsFstatSearch;
  static ComputeFParams CFparams;
  ComputeFBufferVector_RS resampbuffers;  /* used to store the buffered quantities used in repeated calls to ComputeFstatFreqBand_RS */
  
  /* Semicoherent variables */
  static SemiCoherentParams semiCohPar;

  /* coarse grid */
  CoarseGrid coarsegrid;
  CoarseGridPoint thisCgPoint; /* a single coarse-grid point */
  REAL8 dFreqStack; /* frequency resolution of Fstat calculation */
  REAL8 df1dot;  /* coarse grid resolution in spindown */
  UINT4 nf1dot;  /* number of coarse-grid spindown values */
  UINT4 ifdot;  /* counter for coarse-grid spindown values */

  /* fine grid */
  FineGrid finegrid;
  FineGridPoint thisFgPoint; /* a single fine-grid point */
  UINT4 nfreqs_fg, nf1dots_fg=1; /* number of frequency and spindown values */
  REAL8 gammaRefine, sigmasq;  /* refinement factor and variance */

  /* GCT helper variables */
  UINT4 ic, ic2, ic3, ifine, ifreq_fg, if1dot_fg;
  UINT4 fveclength, ifreq;
  INT4  U1idx;
  REAL8 myf0, freq_event, f1dot_event, deltaF, f1dot_eventB1;
  REAL8 dfreq_fg, df1dot_fg, freqmin_fg, f1dotmin_fg, freqband_fg;
  REAL8 u1start, u1win, u1winInv;
  REAL8 freq_tmp, f1dot_tmp;
  REAL4 Fstat, TwoFthreshold, sumTwoF_tmp, TwoF_tmp, sumTwoFmax; /* REAL4 precision of Fstat values */
  UINT4 nc_max;
  REAL8 A1, B1, A2, B2; /* GCT helper variables for faster calculation of u1 or u2 */
  REAL8 pos[3];
  REAL8 vel[3];
  REAL8 acc[3];
  REAL8 cosAlpha, sinAlpha, cosDelta, sinDelta;
  REAL8 nvec[3]; /* unit vector pointing to sky position */

  /* These vars are currently not used, but eventually in the future.
  INT4  U2idx, NumU2idx;
  REAL8 myf0max, u2start, u2end;
  REAL8 u2win, u2winInv;
  */

  /* fstat candidate structure for candidate toplist*/
  toplist_t *semiCohToplist=NULL;

  /* template and grid variables */
  static DopplerSkyScanInit scanInit;   /* init-structure for DopperScanner */
  DopplerSkyScanState thisScan = empty_DopplerSkyScanState; /* current state of the Doppler-scan */
  static PulsarDopplerParams dopplerpos;	       /* current search-parameters */
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

  /* checkpoint filename and index of loop over skypoints */
  /* const CHAR *fnameChkPoint="checkpoint.cpt"; */
  /* FILE *fpChkPoint=NULL; */
  /* UINT4 loopindex, loopcounter; */

  /* user variables */
  BOOLEAN uvar_help = FALSE; 	/* true if -h option is given */
  BOOLEAN uvar_log = FALSE; 	/* logging done if true */
  INT4 uvar_loglevel = 0;

  BOOLEAN uvar_printCand1 = FALSE; 	/* if 1st stage candidates are to be printed */
  BOOLEAN uvar_printFstat1 = FALSE;
  BOOLEAN uvar_semiCohToplist = TRUE; /* if overall first stage candidates are to be output */
  BOOLEAN uvar_useResamp = FALSE;      /* use resampling to compute F-statistic instead of SFT method */
  BOOLEAN uvar_SignalOnly = FALSE;     /* if Signal-only case (for SFT normalization) */
  BOOLEAN uvar_SepDetVeto = FALSE;     /* Do separate detector analysis for the top candidate */
  
  REAL8 uvar_dAlpha = DALPHA; 	/* resolution for flat or isotropic grids -- coarse grid*/
  REAL8 uvar_dDelta = DDELTA;
  REAL8 uvar_f1dot = FDOT; 	/* first spindown value */
  REAL8 uvar_f1dotBand = 0.0; /* range of first spindown parameter */
  REAL8 uvar_Freq = FSTART;
  REAL8 uvar_FreqBand = FBAND;

  REAL8 uvar_dFreq = 0;
  REAL8 uvar_df1dot = 0; /* coarse grid frequency and spindown resolution */

  REAL8 uvar_ThrF = FSTATTHRESHOLD; /* threshold of Fstat to select peaks */
  REAL8 uvar_mismatch1 = MISMATCH; /* metric mismatch for first stage coarse grid */

  REAL8 uvar_minStartTime1 = 0;
  REAL8 uvar_maxEndTime1 = LAL_INT4_MAX;
  REAL8 uvar_dopplerMax = 1.05e-4;

  REAL8 uvar_refTime = 0;
  REAL8 uvar_tStack = 0;
  INT4 uvar_nCand1 = NCAND1; /* number of candidates to be followed up from first stage */

  INT4 uvar_blocksRngMed = BLOCKSRNGMED;
  INT4 uvar_nStacksMax = 1;
  INT4 uvar_Dterms = DTERMS;
  INT4 uvar_SSBprecision = SSBPREC_RELATIVISTIC;
  INT4 uvar_gammaRefine = 1;
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
  uvar_loglevel = lalDebugLevel;
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
  LAL_CALL( LALRegisterREALUserVar(   &status, "dFreq",        0,  UVAR_OPTIONAL, "Frequency resolution (default=1/Tstack)", &uvar_dFreq), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "FreqBand",    'b', UVAR_OPTIONAL, "Search frequency band", &uvar_FreqBand), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "f1dot",        0,  UVAR_OPTIONAL, "Spindown parameter", &uvar_f1dot), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "df1dot",       0,  UVAR_OPTIONAL, "Spindown resolution (default=1/Tstack^2)", &uvar_df1dot), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "f1dotBand",    0,  UVAR_OPTIONAL, "Spindown Range", &uvar_f1dotBand), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "nStacksMax",   0,  UVAR_OPTIONAL, "Maximum No. of 1st stage segments", &uvar_nStacksMax ),&status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "tStack",      'T', UVAR_REQUIRED, "Duration of 1st stage segments (sec)", &uvar_tStack ),&status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "peakThrF",     0,  UVAR_OPTIONAL, "Fstat Threshold", &uvar_ThrF), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "mismatch1",   'm', UVAR_OPTIONAL, "1st stage mismatch", &uvar_mismatch1), &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "gridType1",    0,  UVAR_OPTIONAL, "0=flat,1=isotropic,2=metric,3=file", &uvar_gridType1),  &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "metricType1",  0,  UVAR_OPTIONAL, "0=none,1=Ptole-analytic,2=Ptole-numeric,3=exact", &uvar_metricType1), &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "gammaRefine", 'g', UVAR_OPTIONAL, "Refinement of fine grid (default: use segment times)", &uvar_gammaRefine), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "fnameout",    'o', UVAR_OPTIONAL, "Output filename", &uvar_fnameout), &status);
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

  /* developer user variables */
  LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed", 0, UVAR_DEVELOPER, "RngMed block size", &uvar_blocksRngMed), &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "SSBprecision", 0, UVAR_DEVELOPER, "Precision for SSB transform.", &uvar_SSBprecision),    &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "Dterms",       0, UVAR_DEVELOPER, "No. of terms to keep in Dirichlet Kernel", &uvar_Dterms ), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "skyPointIndex",0, UVAR_DEVELOPER, "Only analyze this skypoint in grid", &uvar_skyPointIndex ), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dopplerMax",   0, UVAR_DEVELOPER, "Max Doppler shift",  &uvar_dopplerMax), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "sftUpsampling",0, UVAR_DEVELOPER, "Upsampling factor for fast LALDemod",  &uvar_sftUpsampling), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "SortToplist",  0, UVAR_DEVELOPER, "Sort toplist by: 0=average2F, 1=numbercount",  &uvar_SortToplist), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "SepDetVeto",   0, UVAR_OPTIONAL,  "Separate detector veto with top candidate", &uvar_SepDetVeto), &status);
  
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

  if ( uvar_blocksRngMed < 1 ) {
    fprintf(stderr, "Invalid Running Median block size\n");
    return( HIERARCHICALSEARCH_EBAD );
  }

  if ( uvar_ThrF < 0 ) {
    fprintf(stderr, "Invalid value of F-statistic threshold\n");
    return( HIERARCHICALSEARCH_EBAD );
  }

  /* 2F threshold for semicoherent stage */
  TwoFthreshold = 2.0 * uvar_ThrF;

  if ( (uvar_SortToplist != 0) && (uvar_SortToplist != 1) ) {
    fprintf(stderr, "Invalid value specified for toplist sorting\n");
    return( HIERARCHICALSEARCH_EBAD );
  }

  /* create toplist -- semiCohToplist has the same structure
     as a fstat candidate, so treat it as a fstat candidate */
  create_gctFStat_toplist(&semiCohToplist, uvar_nCand1, uvar_SortToplist);

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

  /* initializations of coarse and fine grids */
  coarsegrid.list = NULL;
  finegrid.list = NULL;
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
      if ( !(fpFstat1 = fopen( fnameFstatVec1, "wb"))) 	{
        fprintf ( stderr, "Unable to open Fstat file fstatvec1.out for writing.\n");
        return (HIERARCHICALSEARCH_EFILE);
      }
    }

  /*------------ Set up stacks, detector states etc. */
  /* initialize spin range vectors */
  INIT_MEM( spinRange_Temp );

  /* some useful first stage params */
  usefulParams.sftbasename = uvar_DataFiles1;
  usefulParams.nStacks = uvar_nStacksMax;
  usefulParams.tStack = uvar_tStack;

  INIT_MEM ( usefulParams.spinRange_startTime );
  INIT_MEM ( usefulParams.spinRange_endTime );
  INIT_MEM ( usefulParams.spinRange_refTime );
  INIT_MEM ( usefulParams.spinRange_midTime );

  /* copy user specified spin variables at reftime  */
  /* the reference time value in spinRange_refTime will be set in SetUpSFTs() */
  usefulParams.spinRange_refTime.fkdot[0] = uvar_Freq; /* frequency */
  usefulParams.spinRange_refTime.fkdot[1] = uvar_f1dot;  /* 1st spindown */
  usefulParams.spinRange_refTime.fkdotBand[0] = uvar_FreqBand; /* frequency range */
  usefulParams.spinRange_refTime.fkdotBand[1] = uvar_f1dotBand; /* spindown range */

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
    LogPrintf(LOG_DETAIL, "Reference time will be set to mid-time of observation time\n");
    usefulParams.refTime = -1;
  }

  /* for 1st stage: read sfts, calculate detector states */
  LogPrintf( LOG_NORMAL,"Reading input data ... ");
  LAL_CALL( SetUpSFTs( &status, &stackMultiSFT, &stackMultiNoiseWeights, &stackMultiDetStates, &usefulParams), &status);
  fprintf(stderr," done.\n");

  /* some useful params computed by SetUpSFTs */
  tStack = usefulParams.tStack;
  tObs = usefulParams.tObs;
  nStacks = usefulParams.nStacks;
  tStartGPS = usefulParams.tStartGPS;
  midTstack = usefulParams.midTstack;
  startTstack = usefulParams.startTstack;
  endTstack = usefulParams.endTstack;
  tMidGPS = usefulParams.spinRange_midTime.refTime;
  refTimeGPS = usefulParams.spinRange_refTime.refTime;
  fprintf(stderr, "%% --- GPS reference time = %.4f ,  GPS data mid time = %.4f\n",
          XLALGPSGetREAL8(&refTimeGPS), XLALGPSGetREAL8(&tMidGPS) );
  firstSFT = &(stackMultiSFT.data[0]->data[0]->data[0]); /* use  first SFT from  first detector */
  Tsft = 1.0 / firstSFT->deltaF; /* define the length of an SFT (assuming 1/Tsft resolution) */

  if ( uvar_sftUpsampling > 1 )
    {
      LogPrintf (LOG_DEBUG, "Upsampling SFTs by factor %d ... ", uvar_sftUpsampling );
      for (k = 0; k < nStacks; k++) {
        LAL_CALL ( upsampleMultiSFTVector ( &status, stackMultiSFT.data[k], uvar_sftUpsampling, 16 ), &status );
      }
      LogPrintfVerbatim (LOG_DEBUG, "done.\n");
    }

  /*------- set frequency and spindown resolutions and ranges for Fstat and semicoherent steps -----*/

  /* set Fstat calculation frequency resolution (coarse grid) */
  if ( LALUserVarWasSet(&uvar_dFreq) ) {
    dFreqStack = uvar_dFreq;
  }
  else {
    dFreqStack = 1.0/tStack;
  }

  /* set Fstat spindown resolution (coarse grid) */
  if ( LALUserVarWasSet(&uvar_df1dot) ) {
    df1dot = uvar_df1dot;
  }
  else {
    df1dot = 1.0/(tStack*tStack);
  }

  /* number of coarse grid spindown values */
  nf1dot = (UINT4)( usefulParams.spinRange_midTime.fkdotBand[1] / df1dot + 1e-6) + 1;

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



  /**** debugging information ******/
  /* print some debug info about spinrange */
  LogPrintf(LOG_DETAIL, "Frequency and spindown range at refTime (%d): [%f-%f], [%e-%e]\n",
	    usefulParams.spinRange_refTime.refTime.gpsSeconds,
	    usefulParams.spinRange_refTime.fkdot[0],
	    usefulParams.spinRange_refTime.fkdot[0] + usefulParams.spinRange_refTime.fkdotBand[0],
	    usefulParams.spinRange_refTime.fkdot[1],
	    usefulParams.spinRange_refTime.fkdot[1] + usefulParams.spinRange_refTime.fkdotBand[1]);

  LogPrintf(LOG_DETAIL, "Frequency and spindown range at startTime (%d): [%f-%f], [%e-%e]\n",
	    usefulParams.spinRange_startTime.refTime.gpsSeconds,
	    usefulParams.spinRange_startTime.fkdot[0],
	    usefulParams.spinRange_startTime.fkdot[0] + usefulParams.spinRange_startTime.fkdotBand[0],
	    usefulParams.spinRange_startTime.fkdot[1],
	    usefulParams.spinRange_startTime.fkdot[1] + usefulParams.spinRange_startTime.fkdotBand[1]);

  LogPrintf(LOG_DETAIL, "Frequency and spindown range at midTime (%d): [%f-%f], [%e-%e]\n",
	    usefulParams.spinRange_midTime.refTime.gpsSeconds,
	    usefulParams.spinRange_midTime.fkdot[0],
	    usefulParams.spinRange_midTime.fkdot[0] + usefulParams.spinRange_midTime.fkdotBand[0],
	    usefulParams.spinRange_midTime.fkdot[1],
	    usefulParams.spinRange_midTime.fkdot[1] + usefulParams.spinRange_midTime.fkdotBand[1]);

  LogPrintf(LOG_DETAIL, "Frequency and spindown range at endTime (%d): [%f-%f], [%e-%e]\n",
	    usefulParams.spinRange_endTime.refTime.gpsSeconds,
	    usefulParams.spinRange_endTime.fkdot[0],
	    usefulParams.spinRange_endTime.fkdot[0] + usefulParams.spinRange_endTime.fkdotBand[0],
	    usefulParams.spinRange_endTime.fkdot[1],
	    usefulParams.spinRange_endTime.fkdot[1] + usefulParams.spinRange_endTime.fkdotBand[1]);

  /* print debug info about stacks */
  fprintf(stderr, "%% --- Setup, N = %d, T = %.0fs, Tobs = %.0fs, gammaRefine = %f\n",  
          nStacks, tStack, tObs, gammaRefine);

  for (k = 0; k < nStacks; k++) {

    LogPrintf(LOG_DETAIL, "Segment %d ", k+1);
    for ( j = 0; j < (INT4)stackMultiSFT.data[k]->length; j++) {

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
  semiCohPar.refTime = tMidGPS;

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
    
    GET_CHECKPOINT(semiCohToplist, &count, thisScan.numSkyGridPoints * nf1dot, fnameSemiCohCand, NULL);
        
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

#ifdef OUTPUT_TIMING
    clock0 = time(NULL);
#endif

 /* ################## loop over SKY coarse-grid points ################## */
  while(thisScan.state != STATE_FINISHED)
  {

    SkyPosition skypos;

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

    /* get amplitude modulation weights */
    skypos.longitude = thisPoint.Alpha;
    skypos.latitude = thisPoint.Delta;
    skypos.system = COORDINATESYSTEM_EQUATORIAL;

      {  /********Allocate fstat vector memory *****************/

        /* calculate number of bins for Fstat overhead due to residual spin-down */
        semiCohPar.extraBinsFstat = (UINT4)( (0.25 * tObs * df1dot)/dFreqStack + 1e-6) + 1;

        /* calculate total number of bins for Fstat */
        binsFstatSearch = (UINT4)(usefulParams.spinRange_midTime.fkdotBand[0]/dFreqStack + 1e-6) + 1;
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
        
        /* show progress */
#ifdef EAH_BOINC
        LogPrintf( LOG_NORMAL, "%d/%d\n", skyGridCounter+1, ifdot+1 );
#else
        LogPrintf( LOG_NORMAL, "sky:%d f1dot:%d\n", skyGridCounter+1, ifdot+1 );
#endif

        /* ------------- Set up coarse grid --------------------------------------*/
        coarsegrid.length = (UINT4) (binsFstat1);

        /* allocate memory for coarsegrid */
        coarsegrid.list = (CoarseGridPoint *)LALRealloc( coarsegrid.list, coarsegrid.length * sizeof(CoarseGridPoint));
        if ( coarsegrid.list == NULL) {
          fprintf(stderr, "ERROR: Memory allocation  [HierarchSearchGCT.c %d]\n" , __LINE__);
          return(HIERARCHICALSEARCH_EMEM);
        }

        /* Initialize first coarsegrid point */
        thisCgPoint.Uindex = 0;
        thisCgPoint.TwoF = 0.0;

        /* ------------- Set up fine grid --------------------------------------*/

        /* frequency fine-grid borders */
        freqmin_fg = usefulParams.spinRange_midTime.fkdot[0];
        freqband_fg = usefulParams.spinRange_midTime.fkdotBand[0];

        /* fine-grid frequency resolution */
        dfreq_fg = dFreqStack;
        nfreqs_fg = ceil(freqband_fg / dfreq_fg);  /* number of points in frequency */

        /* copy frequency setup parameters to fine-grid struct */
        finegrid.freqmin_fg = freqmin_fg;
        finegrid.dfreq_fg = dfreq_fg;
        finegrid.freqlength = nfreqs_fg ;

        /* fine-grid f1dot resolution */
        nf1dots_fg = ceil(gammaRefine);        /* number of spindown fine-grid points */
        if ( (nf1dots_fg % 2) == 0 ) {    /* if even, add one (to refine symmetrically) */
          nf1dots_fg++;
        }
        df1dot_fg = df1dot / nf1dots_fg;  /* spindown fine-grid  stepsize */

        /* adjust f1dotmin_fg, so that f1dot finegrid is centered around coarse-grid f1dot point */
        f1dotmin_fg = (usefulParams.spinRange_midTime.fkdot[1] + ifdot * df1dot) - df1dot_fg * floor(nf1dots_fg / 2.0);

        /* copy 1st spindown setup parameters to fine-grid struct */
        finegrid.f1dotmin_fg = f1dotmin_fg;
        finegrid.df1dot_fg = df1dot_fg;
        finegrid.f1dotlength = nf1dots_fg;

        /* total number of fine-grid points */
        finegrid.length = nf1dots_fg * nfreqs_fg;

        if(!oldcg) {
          oldcg = coarsegrid.length;
          LogPrintfVerbatim(LOG_NORMAL, "%% --- CG:%d ",coarsegrid.length);
        }
        if(!oldfg) {
          oldfg = finegrid.length;
          LogPrintfVerbatim(LOG_NORMAL, "FG:%ld  f1dotmin_fg:%.13g df1dot_fg:%.13g\n",finegrid.length,f1dotmin_fg,df1dot_fg);
        }
        if((coarsegrid.length != oldcg) || (finegrid.length != oldfg)) {
          LogPrintfVerbatim(LOG_CRITICAL, "ERROR: Grid-sizes disagree!\nPrevious CG:%d FG:%ld, currently CG:%d FG:%ld\n",
                            oldcg,oldfg,coarsegrid.length,finegrid.length);
          return(HIERARCHICALSEARCH_EVAL);
        }
        
        /* reference time for finegrid is midtime */
        finegrid.refTime = tMidGPS;

        /* allocate memory for finegrid points */
        finegrid.list = (FineGridPoint *)LALRealloc( finegrid.list, finegrid.length * sizeof(FineGridPoint));
        if ( finegrid.list == NULL) {
          fprintf(stderr, "ERROR: Memory allocation [HierarchSearchGCT.c %d]\n" , __LINE__);
          return(HIERARCHICALSEARCH_EMEM);
        }

        /* copy sky coarse-grid point to finegrid, because sky is not refined */
        finegrid.alpha = thisPoint.Alpha;
        finegrid.delta = thisPoint.Delta;

        /* initialize finegrid point */
        thisFgPoint.nc=0;
        thisFgPoint.sumTwoF=0.0;

        /* initialize the entire finegrid ( 2F-sum and number count set to 0 ) */
        ic = 0;
        for( ic3 = 0; ic3 < nf1dots_fg; ic3++ ) {
          /*f1dot_tmp = f1dotmin_fg + ic3 * df1dot_fg;*/
          for( ic2 = 0; ic2 < nfreqs_fg; ic2++ ) {
            /*freq_tmp = freqmin_fg + ic2 * dfreq_fg;*/
            finegrid.list[ic] = thisFgPoint;
            ic++;
          }
        }

        /* Keeping track of maximum number count in DIAGNOSE mode */
        nc_max = 0;    /* initialize */
        sumTwoFmax = 0.0;


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

          acc[0] = semiCohPar.acc->data[3*k];
          acc[1] = semiCohPar.acc->data[3*k + 1];
          acc[2] = semiCohPar.acc->data[3*k + 2];


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

          A2 = ( acc[0] * nvec[0] \
               + acc[1] * nvec[1] \
               + acc[2] * nvec[2] ); /* This is \vec a \dot \vec n */

          B2 = ( vel[0] * nvec[0] \
               + vel[1] * nvec[1] \
               + vel[2] * nvec[2] ); /* This is \vec v \dot \vec n */

          A1 = 1.0 + B2;

          /* Setup of the cell size (windows) in u1 */
          u1win = dFreqStack * A1;
          u1winInv = 1.0/u1win; /* Precomputing */

          /* Setup of the cell size (windows) in u2 */
          /* --- currently only u1 is needed.
          u2win = df1dot;
          u2winInv = 1.0/u2win; */

          /* ----------------------------------------------------------------- */
          /************************ Compute F-Statistic ************************/

          /* Set starting frequency for Fstat calculation */
          thisPoint.fkdot[0] = fstatVector.data[k].f0;

          /* Length and spacing of the Fstat vector in frequency */
          fveclength = fstatVector.data[k].data->length;
          deltaF = fstatVector.data[k].deltaF;

          /* Set spindown value for Fstat calculation */
          thisPoint.fkdot[1] = usefulParams.spinRange_midTime.fkdot[1] + ifdot * df1dot;

          /* Frequency at the segment's midpoint for later use */
          f1dot_event = thisPoint.fkdot[1];
          myf0 = thisPoint.fkdot[0] + thisPoint.fkdot[1] * timeDiffSeg;

          if (uvar_useResamp) {
	   
	    /* point the params buffer to the current segment buffer */
	    CFparams.buffer = resampbuffers.data[k];
	    printf("k = %d\n",k);
            /* Resampling method implementation to compute the F-statistic */
            LAL_CALL( COMPUTEFSTATFREQBAND_RS ( &status, &fstatVector.data[k], &thisPoint,
                                               stackMultiSFT.data[k], stackMultiNoiseWeights.data[k],
                                               stackMultiDetStates.data[k], &CFparams), &status);
	    
	    /* repoint the buffer vector element to the potentially modified buffer */
	    resampbuffers.data[k] = CFparams.buffer;

          }
          else {

            /* LALDemod method implementation to compute the F-statistic */
            LAL_CALL( COMPUTEFSTATFREQBAND ( &status, &fstatVector.data[k], &thisPoint,
                                            stackMultiSFT.data[k], stackMultiNoiseWeights.data[k],
                                            stackMultiDetStates.data[k], &CFparams), &status);
          }

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

          f1dot_eventB1 = f1dot_event * B1;

          /* Loop over coarse-grid frequency bins */
          for (ifreq = 0; ifreq < fveclength; ifreq++) {

            /* Get the F-statistic value */
            Fstat = fstatVector.data[k].data->data[ifreq]; /* Recall here it's *F*, not yet 2F */

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
                 thisCgPoint.Uindex = U1idx * NumU2idx + U2idx; */
              thisCgPoint.Uindex = U1idx;
            }

            /* ============ Copy the *2F* value ============ */
            thisCgPoint.TwoF = 2.0 * Fstat;

            /* Add this point to the coarse grid */
            coarsegrid.list[ifreq] = thisCgPoint;

          } /* END: Loop over coarse-grid frequency bins (ifreq) */

          
          /* print fstat vector if required -- mostly for debugging */
          if ( uvar_printFstat1 )
          {
            LAL_CALL( PrintFstatVec ( &status, &fstatVector.data[k], fpFstat1, &thisPoint, refTimeGPS, k+1), &status);
          }

          /* --- Holger: This is not needed in U1-only case. Sort the coarse grid in Uindex --- */
          /* qsort(coarsegrid.list, (size_t)coarsegrid.length, sizeof(CoarseGridPoint), compareCoarseGridUindex); */

          /* ---------- Walk through fine grid and map to coarse grid --------------- */
          ifine = 0;

          for( if1dot_fg = 0; if1dot_fg < finegrid.f1dotlength; if1dot_fg++ ) {

            /* get the 1st spindown of this fine-grid point */
            f1dot_tmp = finegrid.f1dotmin_fg + if1dot_fg * finegrid.df1dot_fg;

            /* pre-compute prouduct */
            f1dot_eventB1 = f1dot_tmp * B1;
                   
            /* get the frequency of this fine-grid point at mid point of segment */
            /* OLD: ifreq_fg = 0; freq_tmp = finegrid.freqmin_fg + ifreq_fg * finegrid.dfreq_fg + f1dot_tmp * timeDiffSeg; */
            freq_tmp = finegrid.freqmin_fg + f1dot_tmp * timeDiffSeg; /* first fine-grid frequency */
             
            /* compute the global-correlation coordinate indices */
            U1idx = ComputeU1idx ( freq_tmp, f1dot_eventB1, A1, u1start, u1winInv );

            if (U1idx < 0) {
              fprintf(stderr,"ERROR: Stepped outside the coarse grid (%d)! \n", U1idx);
              return(HIERARCHICALSEARCH_ECG);
            }

            if (U1idx + finegrid.freqlength >= fveclength) {
              fprintf(stderr,"ERROR: Stepped outside the coarse grid (%d:%d:%d:%d)! \n",
		      U1idx, finegrid.freqlength, U1idx + finegrid.freqlength, fveclength);
              return(HIERARCHICALSEARCH_ECG);
            }

            for( ifreq_fg = 0; ifreq_fg < finegrid.freqlength; ifreq_fg++ ) {

              /* Add the 2F value to the 2F sum */
              TwoF_tmp = coarsegrid.list[U1idx].TwoF;
              sumTwoF_tmp = finegrid.list[ifine].sumTwoF + TwoF_tmp;
              finegrid.list[ifine].sumTwoF = sumTwoF_tmp;

              /* Increase the number count */
              if (TwoF_tmp > TwoFthreshold) {
                finegrid.list[ifine].nc++;
              }
                
#ifdef DIAGNOSISMODE
              /* Keep track of strongest candidate (maximum 2F-sum and maximum number count) */
              if (finegrid.list[ifine].nc > nc_max) {
                nc_max = finegrid.list[ifine].nc;
              }
              if (sumTwoF_tmp > sumTwoFmax) {
                sumTwoFmax = sumTwoF_tmp;
              }
#endif

              /* -------------- Single-trial check ------------- */
              /*
               if ( ifine == 850642 && (k+1) == nStacks ) {
               fprintf(stderr, "MyFineGridPoint,%d f: %.13f fdot: %g  NC: %d  2F: %f\n",
               k+1, finegrid.freqmin_fg + ifreq_fg * finegrid.dfreq_fg,
               finegrid.f1dotmin_fg + if1dot_fg * finegrid.df1dot_fg,
               finegrid.list[ifine].nc, (finegrid.list[ifine].sumTwoF / nStacks)
               );
               }
               */

              U1idx++; /* increment U1 index */
              ifine++; /* increment fine-grid index */

            } /* for( ifreq_fg = 0; ifreq_fg < finegrid.freqlength; ifreq_fg++ ) { */

          } /* for( if1dot_fg = 0; if1dot_fg < finegrid.f1dotlength; if1dot_fg++ ) { */


#ifdef DIAGNOSISMODE
          fprintf(stderr, "  --- Seg: %03d  nc_max: %03d  avesumTwoFmax: %f \n", k, nc_max, sumTwoFmax/(k+1));
#endif

        } /* end: ------------- MAIN LOOP over Segments --------------------*/

        /* ############################################################### */

        if( uvar_semiCohToplist ) {
          /* this is necessary here, because UpdateSemiCohToplist() might set
           a checkpoint that needs some information from here */
          LogPrintf(LOG_DETAIL, "Updating toplist with semicoherent candidates\n");
          LAL_CALL( UpdateSemiCohToplist(&status, semiCohToplist, &finegrid, &usefulParams), &status);
        }

        ifdot++;  /* Increment ifdot counter BEFORE SET_CHECKPOINT */
        
        SHOW_PROGRESS(dopplerpos.Alpha, dopplerpos.Delta,
                      skyGridCounter * nf1dot + ifdot,
                      thisScan.numSkyGridPoints * nf1dot, uvar_Freq, uvar_FreqBand);
#ifdef EAH_BOINC
        SET_CHECKPOINT;
#endif
        
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

#ifdef OUTPUT_TIMING
  {
    time_t tau = time(NULL) - clock0;
    UINT4 Nrefine = nf1dots_fg;
    FILE *timing_fp = fopen ( "HS_timing.dat", "ab" );
    fprintf ( timing_fp, "%d 	%d 	%d 	%d 	%d 	%d 	%d 	%d\n",
	      thisScan.numSkyGridPoints, nf1dot, binsFstatSearch, 2 * semiCohPar.extraBinsFstat, nSFTs, nStacks, Nrefine, tau );
    fclose ( timing_fp );
  }
#endif

  LogPrintf( LOG_NORMAL, "Finished analysis.\n");

  LogPrintf ( LOG_DEBUG, "Writing output ...");
  
  write_hfs_oputput(uvar_fnameout, semiCohToplist);

  /* --- Further analysis with the top candidate if desired --- 
         This veto computes the average F-statistic for each detector
         data stream and compares it to the multi-IFO F-statistic */
  if ( uvar_SepDetVeto ) {

    UINT8 icand, icandMax;
    UINT4 numDetectors, X, topNC;
    REAL4 topTwoF=0.0, maxTopTwoF=0.0;
    Fcomponents FstatSeg;
    ComputeFBuffer cfBuffer2 = empty_ComputeFBuffer;
    PulsarSpins fkdotTMP;
    MultiSFTVector *SFTsSingleDet=NULL;
    MultiNoiseWeights *NoiseSingleDet=NULL;
    MultiDetectorStateSeries *DetStatesSingleDet=NULL;
    
    INIT_MEM( fkdotTMP );
    numDetectors = stackMultiSFT.data[0]->length;
    
    REAL4 aveTwoFstat[numDetectors+1];
    
    if ( (SFTsSingleDet = (MultiSFTVector *)LALCalloc(1, sizeof(MultiSFTVector))) == NULL ){
      fprintf(stderr,"SFTsSingleDet Calloc failed\n");
      return(HIERARCHICALSEARCH_EMEM);
    }
    if ( (NoiseSingleDet = (MultiNoiseWeights *)LALCalloc(1, sizeof(MultiNoiseWeights))) == NULL ){
      fprintf(stderr,"NoiseSingleDet Calloc failed\n");
      return(HIERARCHICALSEARCH_EMEM);
    }
    if ( (DetStatesSingleDet = (MultiDetectorStateSeries *)LALCalloc(1, sizeof(MultiDetectorStateSeries))) == NULL ){
      fprintf(stderr,"DetStatesSingleDet Calloc failed\n");
      return(HIERARCHICALSEARCH_EMEM);
    }
   
    fprintf(stderr, "%% --- Starting separate detector analysis of the top candidate, No. of IFOs: %d\n",numDetectors);
    
    /* Sort the toplist by average 2F to get the strongest candidates */
    sort_gctFStat_toplist_strongest(semiCohToplist);
        
    icand=0; /* At the moment, just the top candidate is analyzed */
    icandMax = icand;
    
    /* find loudest candidate */
    while ( !((*(GCTtopOutputEntry*)semiCohToplist->heap[0]).sumTwoF 
                > (*(GCTtopOutputEntry*)semiCohToplist->heap[icand]).sumTwoF) ) {

      /* Initialize */
      for (X=0; X < (numDetectors+1); X++) 
        aveTwoFstat[X] = 0.0;
      
      thisPoint.Alpha = (*(GCTtopOutputEntry*)semiCohToplist->heap[icand]).Alpha;
      thisPoint.Delta = (*(GCTtopOutputEntry*)semiCohToplist->heap[icand]).Delta;
      fkdotTMP[0] = (*(GCTtopOutputEntry*)semiCohToplist->heap[icand]).Freq;
      fkdotTMP[1] = (*(GCTtopOutputEntry*)semiCohToplist->heap[icand]).F1dot;
      topNC = (*(GCTtopOutputEntry*)semiCohToplist->heap[icand]).nc;
      topTwoF = (*(GCTtopOutputEntry*)semiCohToplist->heap[icand]).sumTwoF;
      /* 
      fprintf(stderr, "  At GPS time %.4f, %.14g %.13g %.13g %.14g %d %.6f  %d\n",
              XLALGPSGetREAL8( &usefulParams.spinRange_refTime.refTime ),
              fkdotTMP[0], thisPoint.Alpha, thisPoint.Delta, fkdotTMP[1], topNC, topTwoF, icand );
      */
      LAL_CALL ( LALExtrapolatePulsarSpins (&status,
                                            thisPoint.fkdot, thisPoint.refTime,
                                            fkdotTMP, refTimeGPS), &status );
      /*
      fprintf(stderr, "  At GPS time %.4f, %.14g %.13g %.13g %.14g %d %.6f\n",
              XLALGPSGetREAL8( &thisPoint.refTime ),
              thisPoint.fkdot[0], thisPoint.Alpha, thisPoint.Delta, thisPoint.fkdot[1],
              topNC, topTwoF );
      */
      for (k = 0; k < nStacks; k++) {
        
        /* --- Compute multi-IFO F-statistic --- */
        LAL_CALL( ComputeFStat ( &status, &FstatSeg, &thisPoint, stackMultiSFT.data[k], 
                                stackMultiNoiseWeights.data[k], stackMultiDetStates.data[k], 
                                &CFparams, &cfBuffer2 ), &status);
        
        if ( uvar_SignalOnly ) {
          FstatSeg.F *= 2.0 / Tsft;
          FstatSeg.F += 2;		
        }
        aveTwoFstat[0] += 2.0 * FstatSeg.F / nStacks;
        
      }
      if (aveTwoFstat[0] > maxTopTwoF) {
        maxTopTwoF = aveTwoFstat[0];
        icandMax = icand;
      }
      fprintf(stderr,"  icand: %" LAL_UINT8_FORMAT "  aveTwoFstat: %f \n",icand,aveTwoFstat[0]);
      icand++;
      
    } /* end while ( !((*(GCTtopOutputEntry*)semiCohToplist->heap[0]).sumTwoF ... */

    
    thisPoint.Alpha = (*(GCTtopOutputEntry*)semiCohToplist->heap[icandMax]).Alpha;
    thisPoint.Delta = (*(GCTtopOutputEntry*)semiCohToplist->heap[icandMax]).Delta;
    fkdotTMP[0] = (*(GCTtopOutputEntry*)semiCohToplist->heap[icandMax]).Freq;
    fkdotTMP[1] = (*(GCTtopOutputEntry*)semiCohToplist->heap[icandMax]).F1dot;
    topNC = (*(GCTtopOutputEntry*)semiCohToplist->heap[icandMax]).nc;
    topTwoF = (*(GCTtopOutputEntry*)semiCohToplist->heap[icandMax]).sumTwoF;
    aveTwoFstat[0] = topTwoF;
    
    fprintf(stderr, "  At GPS time %.4f, %.14g %.13g %.13g %.14g %d %.6f  %" LAL_UINT8_FORMAT " (%" LAL_UINT8_FORMAT ")\n",
            XLALGPSGetREAL8( &usefulParams.spinRange_refTime.refTime ),
            fkdotTMP[0], thisPoint.Alpha, thisPoint.Delta, fkdotTMP[1], topNC, topTwoF, icandMax, icand );
            
    LAL_CALL ( LALExtrapolatePulsarSpins (&status,
                                     thisPoint.fkdot, thisPoint.refTime,
                                     fkdotTMP, refTimeGPS), &status );
    
    fprintf(stderr, "  At GPS time %.4f, %.14g %.13g %.13g %.14g %d %.6f\n",
            XLALGPSGetREAL8( &thisPoint.refTime ),
            thisPoint.fkdot[0], thisPoint.Alpha, thisPoint.Delta, thisPoint.fkdot[1],
            topNC, topTwoF );
    
    /* --- Compute separate-IFO F-statistic for each data segment --- */
    for (k = 0; k < nStacks; k++) {
            
        for (X=0; X < numDetectors; X++) {

        cfBuffer2 = empty_ComputeFBuffer;
        SFTsSingleDet->length = 1;
        SFTsSingleDet->data = &(stackMultiSFT.data[k]->data[X]);
        
        if ( uvar_SignalOnly ) {      
          NoiseSingleDet = NULL;
        }
        else {
          NoiseSingleDet->length = 1;
          NoiseSingleDet->data = &(stackMultiNoiseWeights.data[k]->data[X]);
          NoiseSingleDet->Sinv_Tsft = stackMultiNoiseWeights.data[k]->Sinv_Tsft;
        }
      
        DetStatesSingleDet->length = 1;
        DetStatesSingleDet->data = &(stackMultiDetStates.data[k]->data[X]);
        DetStatesSingleDet->startTime = stackMultiDetStates.data[k]->startTime;
        DetStatesSingleDet->Tspan = stackMultiDetStates.data[k]->Tspan;
        
        LAL_CALL( ComputeFStat ( &status, &FstatSeg, &thisPoint, SFTsSingleDet, 
                              NoiseSingleDet, DetStatesSingleDet, 
                              &CFparams, &cfBuffer2 ), &status);
      
        if ( uvar_SignalOnly ) {
          FstatSeg.F *= 2.0 / Tsft;
          FstatSeg.F += 2;		
        }
        aveTwoFstat[X+1] += 2.0 * FstatSeg.F / nStacks;
        
      }
    }
    
       
    for (X=0; X < (numDetectors+1); X++) {
      if (X>0) {
        fprintf(stderr, "%% --- average2F[%o]= %.6f\t (%s)\t Z= %.4f \n", 
                X, aveTwoFstat[X], 
		(CHAR*) &(stackMultiDetStates.data[0]->data[X-1]->detector.frDetector.name),
		aveTwoFstat[0]/aveTwoFstat[X] );
      }
      else {
        fprintf(stderr, "%% --- average2F[%o]= %.6f\n", X, aveTwoFstat[X]);
      }
    }
    
    XLALEmptyComputeFBuffer ( &cfBuffer2 );
    LALFree(SFTsSingleDet);
    LALFree(NoiseSingleDet);
    LALFree(DetStatesSingleDet);
    
  }
  
  
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

  /* free fine grid and coarse grid */
  if (finegrid.list) {
    LALFree(finegrid.list);
  }
  if (coarsegrid.list) {
    LALFree(coarsegrid.list);
  }
  
  /* free candidate toplist */
  free_gctFStat_toplist(&semiCohToplist);

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
		UsefulStageVariables *in /**< input params */)
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

  INITSTATUS( status, "SetUpSFTs", rcsid );
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

  /* calculate start and end times and tobs from catalog*/
  tStartGPS = catalog->data[0].header.epoch;
  in->tStartGPS = tStartGPS;
  tEndGPS = catalog->data[catalog->length - 1].header.epoch;
  XLALGPSAdd(&tEndGPS, timebase);
  tObs = XLALGPSDiff(&tEndGPS, &tStartGPS);
  in->tObs = tObs;

  /* get sft catalogs for each stack */
  TRY( SetUpStacks( status->statusPtr, &catalogSeq, in->tStack, catalog, in->nStacks), status);

  /* reset number of stacks */
  in->nStacks = catalogSeq.length;

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

  freqmin = freqLo - doppWings - extraBins * deltaFsft;
  freqmax = freqHi + doppWings + extraBins * deltaFsft;

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
    TRY( LALLoadMultiSFTs ( status->statusPtr, stackMultiSFT->data + k,  catalogSeq.data + k,
                           freqmin, freqmax ), status);

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


#ifdef OUTPUT_TIMING
  /* need to count the total number of SFTs */
  nSFTs = 0;
  for ( k = 0; k < in->nStacks; k ++ )
    {
      UINT4 X;
      for ( X=0; X < stackMultiSFT->data[k]->length; X ++ )
        nSFTs += stackMultiSFT->data[k]->data[X]->length;
    } /* for k < stacks */
#endif

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* SetUpSFTs */






/** \brief Breaks up input sft catalog into specified number of stacks

    Loops over elements of the catalog, assigns a bin index and
    allocates memory to the output catalog sequence appropriately.  If
    there are long gaps in the data, then some of the catalogs in the
    output catalog sequence may be of zero length.
*/
void SetUpStacks(LALStatus *status,	   /**< pointer to LALStatus structure */
		 SFTCatalogSequence  *out, /**< Output catalog of sfts -- one for each stack */
		 REAL8 tStack,             /**< Output duration of each stack */
		 SFTCatalog  *in,          /**< Input sft catalog to be broken up into stacks (ordered in increasing time)*/
		 UINT4 nStacksMax )        /**< User specified number of stacks */
{
  UINT4 j, stackCounter, length;
  REAL8 tStart, thisTime;
  REAL8 Tsft;

  INITSTATUS( status, "SetUpStacks", rcsid );
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

  INITSTATUS( status, "PrintCatalogInfo", rcsid );
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

  INITSTATUS( status, "PrintStackInfo", rcsid );
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






/** Read checkpointing file
    This does not (yet) check any consistency of
    the existing results file */
void GetChkPointIndex( LALStatus *status,
		       INT4 *loopindex,
		       const CHAR *fnameChkPoint)
{

  FILE  *fp=NULL;
  UINT4 tmpIndex;
  CHAR lastnewline='\0';

  INITSTATUS( status, "GetChkPointIndex", rcsid );
  ATTATCHSTATUSPTR (status);

  /* if something goes wrong later then lopindex will be 0 */
  *loopindex = 0;

  /* try to open checkpoint file */
  if (!(fp = fopen(fnameChkPoint, "rb")))
    {
      if ( lalDebugLevel )
	fprintf (stdout, "Checkpoint-file '%s' not found.\n", fnameChkPoint);

      DETATCHSTATUSPTR (status);
      RETURN(status);
    }

  /* if we are here then checkpoint file has been found */
  if ( lalDebugLevel )
    fprintf ( stdout, "Found checkpoint-file '%s' \n", fnameChkPoint);

  /* check the checkpointfile -- it should just have one integer
     and a DONE on the next line */
  if ( ( 2 != fscanf (fp, "%" LAL_UINT4_FORMAT "\nDONE%c", &tmpIndex, &lastnewline) ) || ( lastnewline!='\n' ) )
    {
      fprintf ( stdout, "Failed to read checkpoint index from '%s'!\n", fnameChkPoint);
      fclose(fp);

      DETATCHSTATUSPTR (status);
      RETURN(status);
    }

  /* everything seems ok -- set loop index */
  *loopindex = tmpIndex;

  fclose( fp );

  DETATCHSTATUSPTR (status);
  RETURN(status);

}







/** Get SemiCoh candidates toplist */
void UpdateSemiCohToplist(LALStatus *status,
                       toplist_t *list,
                       FineGrid *in,
                       UsefulStageVariables *usefulparams)
{

  BOOLEAN translateSpins = FALSE;
  PulsarSpins fkdot;
  REAL8 freq_tmp, f1dot_tmp;
  UINT4 ifine, if1dot_fg, ifreq_fg, Nsegments;
  INT4 debug;
  GCTtopOutputEntry line;

  INIT_MEM(fkdot);

  INITSTATUS( status, "UpdateSemiCohToplist", rcsid );
  ATTATCHSTATUSPTR (status);

  ASSERT ( list != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( in != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( usefulparams != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

  Nsegments = usefulparams->nStacks;

  /* check if translation to reference time of fine-grid is necessary */
  if  ( XLALGPSDiff( &in->refTime, &usefulparams->spinRange_refTime.refTime) != 0 ) {
    translateSpins = TRUE;
    /*fprintf(stderr,"translateSpins = TRUE\n");*/
  }

  /* ---------- Walk through fine-grid and insert candidates into toplist--------------- */
  ifine = 0;
  for( if1dot_fg = 0; if1dot_fg < in->f1dotlength; if1dot_fg++ ) {

    f1dot_tmp = in->f1dotmin_fg + if1dot_fg * in->df1dot_fg;

    for( ifreq_fg = 0; ifreq_fg < in->freqlength; ifreq_fg++ ) {

      freq_tmp = in->freqmin_fg + ifreq_fg * in->dfreq_fg;

      if ( translateSpins ) { /* propagate fkdot to reference-time  */
        fkdot[0] = freq_tmp;
        fkdot[1] = f1dot_tmp;

        TRY ( LALExtrapolatePulsarSpins (status->statusPtr,
              fkdot, usefulparams->spinRange_refTime.refTime, fkdot, in->refTime), status );

        freq_tmp = fkdot[0];
      }

      line.Freq = freq_tmp;
      line.Alpha = in->alpha;
      line.Delta = in->delta;
      line.F1dot = f1dot_tmp;
      line.nc = in->list[ifine].nc;
      line.sumTwoF = (in->list[ifine].sumTwoF) / Nsegments; /* save the average 2F value */

      debug = insert_into_gctFStat_toplist( list, line);

      ifine++;
    }
  }

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* UpdateSemiCohToplist() */








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

  INITSTATUS( status, "PrintFstatVec", rcsid );
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

    fprintf(fp, "%d %.13g %.12g %.12g %.13g %.6g\n",
            stackIndex, fkdot[0], alpha, delta, fkdot[1], 2*in->data->data[k]);
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

  INITSTATUS( status, "GetSegsPosVelAccEarthOrb", rcsid );
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






/** Comparison function for sorting the coarse grid in u1 and u2 */
int compareCoarseGridUindex(const void *a,const void *b) {
  CoarseGridPoint a1, b1;
  a1 = *((const CoarseGridPoint *)a);
  b1 = *((const CoarseGridPoint *)b);

  if( a1.Uindex < b1.Uindex )
    return(-1);
  else if( a1.Uindex > b1.Uindex)
    return(1);
  else
    return(0);
}







/** Comparison function for sorting the fine grid in number count */
int compareFineGridNC(const void *a,const void *b) {
  FineGridPoint a1, b1;
  a1 = *((const FineGridPoint *)a);
  b1 = *((const FineGridPoint *)b);

  if( a1.nc < b1.nc )
    return(1);
  else if( a1.nc > b1.nc)
    return(-1);
  else
    return(0);
}





/** Comparison function for sorting the fine grid in summed 2F */
int compareFineGridsumTwoF(const void *a,const void *b) {
  FineGridPoint a1, b1;
  a1 = *((const FineGridPoint *)a);
  b1 = *((const FineGridPoint *)b);

  if( a1.sumTwoF < b1.sumTwoF )
    return(1);
  else if( a1.sumTwoF > b1.sumTwoF)
    return(-1);
  else
    return(0);
}



