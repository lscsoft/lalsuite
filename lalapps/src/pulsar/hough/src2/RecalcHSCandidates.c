/*
 *  Copyright (C) 2013 Karl Wette
 *  Copyright (C) 2005 Badri Krishnan, Alicia Sintes, Reinhard Prix, Bernd Machenschalk
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; withoulistt even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */


/**
 * \author Badri Krishnan, Alicia Sintes, Reinhard Prix, Bernd Machenschalk
 * \file
 * \ingroup pulsarApps
 * \brief Program for calculating F-stat values for different time segments
 * and combining them semi-coherently using the Hough transform, and following
 * up candidates using a longer coherent integration.
 *
 * \par  Description
 *
 * This code implements a hierarchical strategy to look for unknown gravitational
 * wave pulsars. It scans through the parameter space using a less sensitive but
 * computationally inexpensive search and follows up the candidates using
 * more sensitive methods.
 *
 * \par Algorithm
 *
 * Currently the code does a single stage hierarchical search using the Hough
 * algorithm and follows up the candidates using a full coherent integration.
 *
 * - The user specifies a directory containing SFTs, and the number of
 * \e stacks that this must be broken up into.  Stacks are chosen by
 * breaking up the total time spanned by the sfts into equal
 * portions, and then choosing the sfts which lie in each stack
 * portion.  The SFTs can be from multiple IFOs.
 *
 * - The user specifies a region in parameter space to search over.
 * The code sets up a grid (the "coarse" grid) in this region and
 * calculates the F-statistic for each stack at each point of the
 * coarse grid, for the entire frequency band.  Alternatively, the
 * user can specify a sky-grid file.
 *
 * - The different Fstat vactors for each stack are combined using a
 * semi-coherent method such as the Hough transform or stack slide
 * algorithms.
 *
 * - For Hough, a threshold is set on the F-statistic to convert the
 * F-statistic vector into a vector of 0s and 1s known as a \e
 * peakgram -- there is one peakgram for each stack and for each
 * grid point.  These peakgrams are combined using the Hough algorithm.
 *
 * - For stack-slide, we add the different Fstatistic values to get
 * the summed Fstatistic power.
 *
 * - The semi-coherent part of the search constructs a grid (the
 * "fine" grid) in a small patch around every coarse grid point, and
 * combines the different stacks following the \e master equation
 * \f[ f(t) - F_0(t) = \xi(t).(\hat{n} - \hat{n}_0) \f]
 * where
 * \f[ F_0 = f_0 + \sum \Delta f_k \frac{(\Delta t)^k}{k!}  \f]
 * Here \f$ \hat{n}_0 \f$ is the sky-point at which the F-statistic is
 * calculated and \f$ \Delta f_k \f$ is the \e residual spindown
 * parameter.  For details see Phys.Rev.D 70, 082001 (2004).  The
 * size of the patch depends on the validity of the above master
 * equation.
 *
 * - The output of the Hough search is a \e number \e count at each point
 * of the grid.  A threshold is set on the number count, leading to
 * candidates in parameter space.  For stack slide, instead of the
 * number count, we get the summed Fstatistic power.  Alternatively,
 * the user can specify exactly how many candidates should be
 * followed up for each coarse grid point.
 *
 * - These candidates are followed up using a second set of SFTs (also
 * specified by the user).  The follow up consists of a full
 * coherent integration, i.e. the F-statistic is calculated for the
 * whole set of SFTs without breaking them up into stacks.  The user
 * can choose to output the N highest significance candidates.
 *
 * \par Questions/To-do
 *
 * - Should we over-resolve the Fstat calculation to reduce loss in signal power?  We would
 * still calculate the peakgrams at the 1/T resolution, but the peak selection would
 * take into account Fstat values over several over-resolved bins.
 *
 * - What is the best grid for calculating the F-statistic?  At first glance, the
 * Hough patch and the metric F-statistic patch do not seem to be compatible.  If we
 * use the Hough patches to break up the sky, it could be far from optimal.  What is
 * the exact connection between the two grids?
 *
 * - Implement multiple semi-coherent stages
 *
 * - Get timings and optimize the pipeline parameters
 *
 * - Checkpointing for running on Einstein\@Home
 *
 * - Incorporate stack slide as an alternative to Hough in the semi-coherent stages
 *
 * - ....
 *
 */

#include <lal/LALString.h>
#include "HierarchicalSearch.h"
#include <../../GCT/RecalcToplistStats.h>

#define TRUE (1==1)
#define FALSE (1==0)

/* Hooks for Einstein\@Home / BOINC
   These are defined to do nothing special in the standalone case
   and will be set in boinc_extras.h if EAH_BOINC is set
 */
#ifdef EAH_BOINC
#include "EinsteinAtHome/hs_boinc_extras.h"
#else /* EAH_BOINC */
/* checkpointing */
#define HS_CHECKPOINTING 0 /* no checkpointing in the non-BOINC case (yet) */
#define GET_CHECKPOINT(toplist,total,count,outputname,cptname) *total=0;
#define INSERT_INTO_HOUGHFSTAT_TOPLIST insert_into_houghFStat_toplist
#define SHOW_PROGRESS(rac,dec,tpl_count,tpl_total,freq,fband)
#define SET_CHECKPOINT
/* BOINC */
#define MAIN main
#endif /* EAH_BOINC */

/* These might have been set differently in hs_boinc_extras.h or ComputeFStatREAL4.h */
#ifndef GPUREADY_DEFAULT
#define GPUREADY_DEFAULT 0
#endif
#ifndef REARRANGE_SFT_DATA
#define REARRANGE_SFT_DATA
#endif
#ifndef INITIALIZE_COPROCESSOR_DEVICE
#define INITIALIZE_COPROCESSOR_DEVICE
#endif
#ifndef UNINITIALIZE_COPROCESSOR_DEVICE
#define UNINITIALIZE_COPROCESSOR_DEVICE
#endif


BOOLEAN uvar_printMaps = FALSE; /**< global variable for printing Hough maps */
BOOLEAN uvar_printGrid = FALSE; /**< global variable for printing Hough grid */
BOOLEAN uvar_printStats = FALSE;/**< global variable for calculating Hough map stats */
BOOLEAN uvar_dumpLUT = FALSE;  	/**< global variable for printing Hough look-up-tables for debugging */
BOOLEAN uvar_validateLUT = FALSE;

#define HSMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define HSMIN(x,y) ( (x) < (y) ? (x) : (y) )

#define BLOCKSIZE_REALLOC 50

/** Useful stuff for a single stage of the Hierarchical search */
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
  EphemerisData *edat;   /**< ephemeris data for LALBarycenter */
  LIGOTimeGPSVector *midTstack;    /**< timestamps vector for mid time of each stack */
  LIGOTimeGPSVector *startTstack;  /**< timestamps vector for start time of each stack */
  LIGOTimeGPS minStartTimeGPS;     /**< all sft data must be after this time */
  LIGOTimeGPS maxStartTimeGPS;       /**< all sft data must be before this GPS time */
  UINT4 blocksRngMed;              /**< blocksize for running median noise floor estimation */
  UINT4 Dterms;                    /**< size of Dirichlet kernel for Fstat calculation */
  REAL8 dopplerMax;                /**< extra sft wings for doppler motion */
  SSBprecision SSBprec;            /**< SSB transform precision */
  LALStringVector *detectorIDs;    /**< vector of detector IDs */
} UsefulStageVariables;


/* functions for printing various stuff */
void ComputeStackNoiseWeights( LALStatus *status, REAL8Vector **out, FstatInputVector* Fstat_in_vec );

void ComputeStackNoiseAndAMWeights( LALStatus *status, REAL8Vector *out, FstatInputVector* Fstat_in_vec, SkyPosition skypos);

void GetStackVelPos( LALStatus *status, REAL8VectorSequence **velStack, REAL8VectorSequence **posStack,
		     FstatInputVector* Fstat_in_vec );


void SetUpSFTs( LALStatus *status, FstatInputVector** p_Fstat_in_vec, UsefulStageVariables *in );

void PrintFstatVec (LALStatus *status, REAL4FrequencySeries *in, FILE *fp, PulsarDopplerParams *thisPoint,
		    LIGOTimeGPS  refTime, INT4 stackIndex);

void PrintHoughGrid(LALStatus *status, HOUGHPatchGrid *patch, HOUGHDemodPar  *parDem, CHAR *fnameOut, INT4 iHmap);

void PrintSemiCohCandidates(LALStatus *status, SemiCohCandidateList *in, FILE *fp, LIGOTimeGPS refTime);

void PrintHoughHistogram(LALStatus *status, UINT8Vector *hist, CHAR *fnameOut);

void PrintCatalogInfo( LALStatus  *status, const SFTCatalog *catalog, FILE *fp);

void PrintStackInfo( LALStatus  *status, const SFTCatalogSequence *catalogSeq, FILE *fp);

void GetSemiCohToplist(LALStatus *status, toplist_t *list, SemiCohCandidateList *in, REAL8 meanN, REAL8 sigmaN);

void ComputeNumExtraBins(LALStatus *status, SemiCoherentParams *par, REAL8 fdot, REAL8 f0, REAL8 deltaF);

void DumpLUT2file(LALStatus *status, HOUGHptfLUT *lut, HOUGHPatchGrid *patch, CHAR *basename, INT4 ind);

void ValidateHoughLUT(LALStatus *status, HOUGHptfLUT *lut, HOUGHPatchGrid  *patch, CHAR *basename, INT4 ind, REAL4 alpha, REAL4 delta, REAL4 weight);

void GetXiInSingleStack (LALStatus         *status,
			 HOUGHSizePar      *size,
			 HOUGHDemodPar     *par);

void GetHoughPatchTopCandidate (LALStatus *status, SemiCohCandidate *topCand, HOUGHMapTotal *ht, HOUGHPatchGrid *patch, HOUGHDemodPar *parDem );

void RCComputeFstatHoughMap (LALStatus *status,
			     SemiCohCandidateList *out,
			     HOUGHPeakGramVector *pgV,
			     SemiCoherentParams *params,
                             INT8 fBin0
                             );


/* default values for input variables */
#define BLOCKSRNGMED 		101 	/**< Default running median window size */

#define NFDOT  			10    	/**< Default size of hough cylinder of look up tables */
#define DTERMS 			8     	/**< Default number of dirichlet kernel terms for calculating Fstat */
#define FSTATTHRESHOLD 		2.6	/**< Default threshold on Fstatistic for peak selection */
#define FNAMEOUT 		"./out/HS.dat"  /**< Default output file basename */
#define PIXELFACTOR 		2.0
#ifndef LAL_INT4_MAX
#define LAL_INT4_MAX 		2147483647
#endif

/* a global pointer to MAIN()s head of the LALStatus structure,
   made global so a signal handler can read it */
LALStatus *global_status;

#ifdef OUTPUT_TIMING
time_t clock0;
UINT4 nSFTs;
UINT4 nStacks;
UINT4 nSkyRefine;
#endif

/* ==================== NEW code to 'follow-up' Hough candidates for topcand, F1, F2, multi-F ==================== */

typedef struct
{
  REAL8 Freq;
  REAL8 Alpha;
  REAL8 Delta;
  REAL8 f1dot;
  REAL8 sig;
} HoughCandidate;

typedef struct
{
  UINT4 length;
  HoughCandidate *data;	/**< 'length' array of HoughCandidate entries */
  REAL8 FreqMin;	/**< smallest candidate frequency */
  REAL8 FreqBand;	/**< candidate frequency band */
  REAL8 f1dotMin;	/** smallest candidate f1dot */
  REAL8 f1dotBand;	/**< candidate f1dot Band */
  REAL8 AlphaMin;
  REAL8 AlphaBand;
  REAL8 DeltaMin;
  REAL8 DeltaBand;
} HoughCandidateList;

HoughCandidateList *XLALLoadHoughCandidateList ( const char *fname, REAL8 FreqShift );
HoughCandidateList *XLALCreateHoughCandidateList ( UINT4 length );
void XLALDestroyHoughCandidateList ( HoughCandidateList * list );

/* ==================== ==================== */

int MAIN( int argc, char *argv[]) {
  LALStatus status = blank_status;

  /* temp loop variables: generally k loops over stacks */
  UINT4 k;

  /* in general any variable ending with 1 is for the
     first stage, 2 for the second and so on */

  /* timestamp vectors */
  LIGOTimeGPSVector *midTstack=NULL;
  LIGOTimeGPSVector *startTstack=NULL;


  LIGOTimeGPS XLAL_INIT_DECL(refTimeGPS);
  LIGOTimeGPS XLAL_INIT_DECL(tMidGPS);

  /* velocities and positions at midTstack */
  REAL8VectorSequence *velStack=NULL;
  REAL8VectorSequence *posStack=NULL;

  /* weights for each stack*/
  REAL8Vector *weightsV=NULL;
  REAL8Vector *weightsNoise=NULL;

  /* duration of each stack */
  REAL8 tStack;

  /* sft related stuff */
  static LIGOTimeGPS minStartTimeGPS, maxStartTimeGPS;


  /* some useful variables for each stage */
  UsefulStageVariables usefulParams;

  /* number of stacks -- not necessarily same as uvar_nStacks! */
  UINT4 nStacks;

  /* F-statistic computation related stuff */
  static REAL4FrequencySeriesVector fstatVector;	/* Fstatistic vectors for each stack */
  FstatInputVector* Fstat_in_vec = NULL;		// Vector of Fstat input data structures for XLALComputeFstat(), one per stack
  FstatResults* Fstat_res = NULL;			// Pointer to Fstat results structure, will be allocated by XLALComputeFstat()
  FstatQuantities Fstat_what = FSTATQ_2F;		// Quantities to be computed by XLALComputeFstat()
  UINT4 binsFstat1, binsFstatSearch;

  /* hough variables */
  static HOUGHPeakGramVector pgV;
  static SemiCoherentParams semiCohPar;
  static SemiCohCandidateList semiCohCandList;
  REAL8 alphaPeak, sumWeightSquare, meanN=0, sigmaN=0;

  /* fstat candidate structure */
  toplist_t *semiCohToplist=NULL;

  /* template and grid variables */
  static PulsarDopplerParams thisPoint;

  /* temporary storage for spinrange vector */
  static PulsarSpinRange spinRange_Temp;

  /* variables for logging */
  CHAR *fnamelog=NULL;
  FILE *fpLog=NULL;
  CHAR *logstr=NULL;

  /* output candidate files and file pointers */
  CHAR *fnameSemiCohCand=NULL;
  CHAR *fnameFstatVec1=NULL;
  FILE *fpSemiCoh=NULL;
  FILE *fpFstat1=NULL;

  /* checkpoint filename and index of loop over skypoints */
  /* const CHAR *fnameChkPoint="checkpoint.cpt"; */
  /*   FILE *fpChkPoint=NULL; */
  /*   UINT4 loopindex, loopcounter; */

  /* user variables */
  BOOLEAN uvar_help = FALSE; 	/* true if -h option is given */
  BOOLEAN uvar_log = FALSE; 	/* logging done if true */

  BOOLEAN uvar_printFstat1 = FALSE;
  BOOLEAN uvar_useWeights  = FALSE;
  BOOLEAN uvar_outputFX    = TRUE; /* Do additional analysis for all toplist candidates, output F and FXvector for postprocessing */
  /* BOOLEAN uvar_useFstatWeights = TRUE; /\* Use noise weights in final toplist Fstat computation? *\/ */

  REAL8 uvar_peakThrF = FSTATTHRESHOLD; /* threshold of Fstat to select peaks */

  REAL8 uvar_pixelFactor = PIXELFACTOR;
  REAL8 uvar_minStartTime1 = 0;
  REAL8 uvar_maxStartTime1 = LAL_INT4_MAX;
  REAL8 uvar_dopplerMax = 1.05e-4;

  REAL8 uvar_refTime = 0;
  REAL8 uvar_semiCohPatchX = 0;
  REAL8 uvar_semiCohPatchY = 0;
  REAL8 uvar_tStack = 0;

  INT4 uvar_blocksRngMed = BLOCKSRNGMED;
  INT4 uvar_nStacksMax = 1;
  INT4 uvar_Dterms = 8;
  INT4 uvar_SSBprecision = SSBPREC_RELATIVISTIC;

  CHAR *uvar_ephemEarth = NULL;
  CHAR *uvar_ephemSun = NULL;

  CHAR *uvar_fnameout = NULL;
  CHAR *uvar_DataFiles1 = NULL;
  BOOLEAN uvar_version = 0;
  CHAR *uvar_outputSingleSegStats = NULL; /* Additionally output single-segment Fstats for each final toplist candidate */

  CHAR *uvar_followupList = NULL;	/* Hough candidate list to be 'followed up': compute Hough top-cand, F1, F2, multi-F */
  REAL8 uvar_WU_Freq = 0;	/* if given, use to compute frequency-correction of input-candidates coming from HS-code */
  REAL8 uvar_WU_dFreq;		/* Frequency-spacing of original HS search that produced the input candidates */
  REAL8 uvar_WU_FreqBand;	/* if given, used to load original SFT band, in order to obtain identical noise-weights */
  REAL8 uvar_WU_f1dot;		/* if given, used to load original SFT band, in order to obtain identical noise-weights */
  REAL8 uvar_WU_f1dotBand;	/* if given, used to load original SFT band, in order to obtain identical noise-weights */

#ifndef GPUREADY_DEFAULT
#define GPUREADY_DEFAULT 0
#endif
  BOOLEAN uvar_GPUready = GPUREADY_DEFAULT;
  global_status = &status;

  BOOLEAN uvar_correctFreqs = TRUE;

  uvar_ephemEarth = XLALStringDuplicate("earth00-19-DE405.dat.gz");
  uvar_ephemSun   = XLALStringDuplicate("sun00-19-DE405.dat.gz");

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
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "log",          0,  UVAR_OPTIONAL, "Write log file", &uvar_log), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "useWeights",   0,  UVAR_OPTIONAL, "Weight each stack using noise and AM?", &uvar_useWeights ), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "outputFX",     0,  UVAR_OPTIONAL, "Additional analysis for toplist candidates, output 2FX?", &uvar_outputFX ), &status);
  /* LAL_CALL( LALRegisterBOOLUserVar(   &status, "useFstatWeights", 0,  UVAR_DEVELOPER, "Use noise weights in toplist Fstat computation?", &uvar_useFstatWeights ), &status); */
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "DataFiles1",   0,  UVAR_REQUIRED, "1st SFT file pattern", &uvar_DataFiles1), &status);

  LAL_CALL( LALRegisterINTUserVar(    &status, "nStacksMax",   0,  UVAR_OPTIONAL, "Maximum No. of 1st stage stacks", &uvar_nStacksMax ),&status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "tStack",       0,  UVAR_REQUIRED, "Duration of 1st stage stacks (sec)", &uvar_tStack ),&status);

  /* ---------- */
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "followupList", 0,  UVAR_REQUIRED, "Hough candidate list to target (Freq Alpha Delta f1dot sig)", &uvar_followupList),  &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "WU_Freq",      0,  UVAR_OPTIONAL, "WU start-frequency: use to correct HS frequency-offset bug", &uvar_WU_Freq ),&status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "WU_dFreq",     0,  UVAR_REQUIRED, "Frequency-spacing of original HS search that produced the input candidates", &uvar_WU_dFreq ),&status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "WU_FreqBand",  0,  UVAR_OPTIONAL, "Used to load original SFT band, in order to obtain identical noise-weights", &uvar_WU_FreqBand ),&status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "WU_f1dot",     0,  UVAR_OPTIONAL, "Used to load original SFT band, in order to obtain identical noise-weights", &uvar_WU_f1dot ),&status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "WU_f1dotBand", 0,  UVAR_OPTIONAL, "Used to load original SFT band, in order to obtain identical noise-weights", &uvar_WU_f1dotBand ),&status);

  /* ---------- */

  LAL_CALL( LALRegisterREALUserVar(   &status, "pixelFactor",  0,  UVAR_OPTIONAL, "Semi coh. sky resolution = 1/v*pixelFactor*f*Tcoh", &uvar_pixelFactor), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "semiCohPatchX",0,  UVAR_OPTIONAL, "Semi coh. sky grid size (default = 1/f*Tcoh*Vepi)", &uvar_semiCohPatchX), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "semiCohPatchY",0,  UVAR_OPTIONAL, "Semi coh. sky grid size (default = 1/f*Tcoh*Vepi)", &uvar_semiCohPatchY), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "fnameout",    'o', UVAR_REQUIRED, "Output fileneme", &uvar_fnameout), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "peakThrF",     0,  UVAR_OPTIONAL, "Fstat Threshold", &uvar_peakThrF), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "refTime",      0,  UVAR_OPTIONAL, "Ref. time for pulsar pars [Default: mid-time]", &uvar_refTime), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "ephemEarth",   0,  UVAR_OPTIONAL, "Location of Earth ephemeris file", &uvar_ephemEarth),  &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "ephemSun",     0,  UVAR_OPTIONAL, "Location of Sun ephemeris file", &uvar_ephemSun),  &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "minStartTime1",0,  UVAR_OPTIONAL, "1st stage: Only use SFTs with timestamps starting from (including) this GPS time", &uvar_minStartTime1), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "maxStartTime1",0,  UVAR_OPTIONAL, "1st stage: Only use SFTs with timestamps up to (excluding) this GPS time",   &uvar_maxStartTime1),   &status);


  /* developer user variables */
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printFstat1",  0,  UVAR_DEVELOPER,  "Print 1st stage Fstat vectors", &uvar_printFstat1), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed", 0, UVAR_DEVELOPER, "RngMed block size", &uvar_blocksRngMed), &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "SSBprecision", 0, UVAR_DEVELOPER, "Precision for SSB transform.", &uvar_SSBprecision),    &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printMaps",    0, UVAR_DEVELOPER, "Print Hough maps -- for debugging", &uvar_printMaps), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printGrid",    0, UVAR_DEVELOPER, "Print Hough fine grid -- for debugging", &uvar_printGrid), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "dumpLUT",      0, UVAR_DEVELOPER, "Print Hough look-up-tables -- for debugging", &uvar_dumpLUT), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "validateLUT",  0, UVAR_DEVELOPER, "Validate Hough look-up-tables -- for debugging", &uvar_validateLUT), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printStats",   0, UVAR_DEVELOPER, "Print Hough map statistics", &uvar_printStats), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "Dterms",       0, UVAR_DEVELOPER, "No.of terms to keep in Dirichlet Kernel", &uvar_Dterms ), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dopplerMax",   0, UVAR_DEVELOPER, "Max Doppler shift",  &uvar_dopplerMax), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "GPUready",     0, UVAR_DEVELOPER, "Use single-precision 'GPU-ready' core routines", &uvar_GPUready), &status);
  LAL_CALL ( LALRegisterBOOLUserVar(  &status, "version",     'V', UVAR_SPECIAL,  "Output version information", &uvar_version), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "outputSingleSegStats", 0,  UVAR_OPTIONAL, "Base filename for single-segment Fstat output (1 file per final toplist candidate!)", &uvar_outputSingleSegStats),  &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "correctFreqs", 0, UVAR_DEVELOPER, "Correct candidate output frequencies (ie fix bug #147). Allows reproducing 'historical results'", &uvar_correctFreqs), &status);

  /* read all command line variables */
  LAL_CALL( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    return(0);


  /* assemble version string */
  CHAR *VCSInfoString;
  if ( (VCSInfoString = XLALGetVersionString(0)) == NULL ) {
    XLALPrintError("XLALGetVersionString(0) failed.\n");
    return( HIERARCHICALSEARCH_EBAD );
  }
  LogPrintfVerbatim( LOG_DEBUG, "Code-version: %s", VCSInfoString );

  if ( uvar_version )
    {
      printf ("%s\n", VCSInfoString );
      return (0);
    }

  /* set log-level */
#ifdef EAH_LOGLEVEL
  LogSetLevel ( EAH_LOGLEVEL );
#else
  LogSetLevel ( lalDebugLevel );
#endif

  if ( uvar_nStacksMax < 1) {
    fprintf(stderr, "Invalid number of stacks\n");
    return( HIERARCHICALSEARCH_EBAD );
  }


  if ( uvar_blocksRngMed < 1 ) {
    fprintf(stderr, "Invalid Running Median block size\n");
    return( HIERARCHICALSEARCH_EBAD );
  }

  if ( uvar_peakThrF < 0 ) {
    fprintf(stderr, "Invalid value of Fstatistic threshold\n");
    return( HIERARCHICALSEARCH_EBAD );
  }

  /* probability of peak selection */
  alphaPeak = (1+uvar_peakThrF)*exp(-uvar_peakThrF);

  /* followup-list given: load candidate list to follow up */
  if ( !uvar_followupList ) {
    XLALPrintError ("No candidate followup-list obtained ... ERROR!\n");
    return( HIERARCHICALSEARCH_EBAD );
  }

  REAL8 FreqShift = 0;	/* frequency-shift correction for frequency-quantization offset bug */
  if ( uvar_WU_Freq > 0 )
    FreqShift = uvar_WU_Freq - uvar_WU_dFreq * floor ( uvar_WU_Freq / uvar_WU_dFreq + 0.5 );
  LogPrintf (LOG_DEBUG, "Applying frequency-correction shift of %.9g Hz \n", FreqShift );

  HoughCandidateList *InputCandList;
  if ( ( InputCandList = XLALLoadHoughCandidateList ( uvar_followupList, FreqShift )) == NULL ) {
    XLALPrintError ("Failed to load followup-list '%s'\n", uvar_followupList );
    return( HIERARCHICALSEARCH_EBAD );
  }

  /* create toplist with exactly the length of the input candidate-list */
  create_houghFStat_toplist(&semiCohToplist, InputCandList->length );

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
	/*exit*/ return(HIERARCHICALSEARCH_EFILE);
      }

      /* get the log string */
      LAL_CALL( LALUserVarGetLog(&status, &logstr, UVAR_LOGFMT_CFGFILE), &status);

      fprintf( fpLog, "## Log file for HierarchicalSearch.c\n\n");
      fprintf( fpLog, "# User Input:\n");
      fprintf( fpLog, "#-------------------------------------------\n");
      fprintf( fpLog, "%s", logstr);
      LALFree(logstr);

      /* add code version ID */
      fprintf ( fpLog, "%s", VCSInfoString );

      fclose (fpLog);

      LALFree(fnamelog);

    } /* end of logging */


  /*--------- Some initializations ----------*/

  /* read in ephemeris data */
  EphemerisData * edat;
  XLAL_CHECK ( (edat = XLALInitBarycenter ( uvar_ephemEarth, uvar_ephemSun )) != NULL, XLAL_EFUNC );

  XLALGPSSetREAL8(&minStartTimeGPS, uvar_minStartTime1);
  XLALGPSSetREAL8(&maxStartTimeGPS, uvar_maxStartTime1);

  /* create output Hough file */
  fnameSemiCohCand = LALCalloc( strlen(uvar_fnameout) + 1, sizeof(CHAR) );
  if ( fnameSemiCohCand == NULL) {
    fprintf(stderr, "error allocating memory [HierarchicalSearch.c %d]\n" , __LINE__);
    return(HIERARCHICALSEARCH_EMEM);
  }

  strcpy(fnameSemiCohCand, uvar_fnameout);

  if ( uvar_printFstat1 )
    {
      const CHAR *append = "_fstatVec1.dat";
      fnameFstatVec1 = LALCalloc( strlen(uvar_fnameout) + strlen(append) + 1, sizeof(CHAR) );
      strcpy(fnameFstatVec1, uvar_fnameout);
      strcat(fnameFstatVec1, append);
      if ( !(fpFstat1 = fopen( fnameFstatVec1, "wb")))
	{
	  fprintf ( stderr, "Unable to open Fstat file fstatvec1.out for writing.\n");
	  return (HIERARCHICALSEARCH_EFILE);
	}
    }

  /*------------ Set up stacks, noise weights, detector states etc. */
  /* initialize spin range vectors */
  XLAL_INIT_MEM(spinRange_Temp);

  /* some useful first stage params */
  usefulParams.sftbasename = uvar_DataFiles1;
  usefulParams.nStacks = uvar_nStacksMax;
  usefulParams.tStack = uvar_tStack;
  usefulParams.SSBprec = uvar_SSBprecision;

  XLAL_INIT_MEM ( usefulParams.spinRange_startTime );
  XLAL_INIT_MEM ( usefulParams.spinRange_endTime );
  XLAL_INIT_MEM ( usefulParams.spinRange_refTime );
  XLAL_INIT_MEM ( usefulParams.spinRange_midTime );

  /* either use WU-original band parameters if given, otherwise adapt to candidate-list range */
  REAL8 Freq, FreqBand, f1dot, f1dotBand;
  /* frequency */
  if ( XLALUserVarWasSet(&uvar_WU_Freq) )
    Freq = uvar_WU_Freq;
  else
    Freq = InputCandList->FreqMin;
  /* frequency range */
  if ( XLALUserVarWasSet(&uvar_WU_FreqBand) )
    FreqBand = uvar_WU_FreqBand;
  else
    FreqBand = InputCandList->FreqBand;
  /* 1st spindown */
  if ( XLALUserVarWasSet(&uvar_WU_f1dot) )
    f1dot = uvar_WU_f1dot;
  else
    f1dot = InputCandList->f1dotMin;
  /* spindown range */
  if ( XLALUserVarWasSet(&uvar_WU_f1dotBand) )
    f1dotBand = uvar_WU_f1dotBand;
  else
    f1dotBand = InputCandList->f1dotBand;

  usefulParams.spinRange_refTime.fkdot[0] = Freq;
  usefulParams.spinRange_refTime.fkdotBand[0] = FreqBand;
  usefulParams.spinRange_refTime.fkdot[1] = f1dot;
  usefulParams.spinRange_refTime.fkdotBand[1] = f1dotBand;

  usefulParams.edat = edat;
  usefulParams.minStartTimeGPS = minStartTimeGPS;
  usefulParams.maxStartTimeGPS = maxStartTimeGPS;
  usefulParams.blocksRngMed = uvar_blocksRngMed;
  usefulParams.Dterms = uvar_Dterms;
  usefulParams.dopplerMax = uvar_dopplerMax;

  /* set reference time for pular parameters */
  if ( LALUserVarWasSet(&uvar_refTime))
    usefulParams.refTime = uvar_refTime;
  else {
    LogPrintf(LOG_DETAIL, "Reference time will be set to mid-time of observation time\n");
    usefulParams.refTime = -1;
  }

  /* for 1st stage: read sfts, calculate multi-noise weights and detector states */
  LogPrintf (LOG_DEBUG, "Reading SFTs and setting up stacks ... ");
  LAL_CALL( SetUpSFTs( &status, &Fstat_in_vec, &usefulParams ), &status);
  LogPrintfVerbatim (LOG_DEBUG, "done\n");

  /* some useful params computed by SetUpSFTs */
  tStack = usefulParams.tStack;
  nStacks = usefulParams.nStacks;
  /* currently unused: LIGOTimeGPS tStartGPS = usefulParams.tStartGPS; */
  midTstack = usefulParams.midTstack;
  startTstack = usefulParams.startTstack;
  tMidGPS = usefulParams.spinRange_midTime.refTime;
  refTimeGPS = usefulParams.spinRange_refTime.refTime;
  LogPrintf(LOG_DETAIL, "GPS Reference Time = %d\n", refTimeGPS.gpsSeconds);

  /*------- set frequency and spindown resolutions and ranges for Fstat and semicoherent steps -----*/

  /*---------- compute noise weight for each stack and initialize total weights vector
     -- for debugging purposes only -- we will calculate noise and AM weights later ----------*/
  if (lalDebugLevel) {
    LAL_CALL( ComputeStackNoiseWeights( &status, &weightsNoise, Fstat_in_vec), &status);
  }

  /* weightsV is the actual weights vector used */
  LAL_CALL( LALDCreateVector( &status, &weightsV, nStacks), &status);


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
  if ( weightsNoise ) {
    for (k = 0; k < nStacks; k++) {
      LogPrintf(LOG_DETAIL, "Stack %d (GPS start time = %d, Noise weight = %f )\n ",
                k, startTstack->data[k].gpsSeconds, weightsNoise->data[k]);
    } /* loop over stacks */
  }



  /*---------- set up F-statistic calculation stuff ---------*/

  /* set reference time for calculating Fstatistic */
  /* thisPoint.refTime = tStartGPS; */
  thisPoint.refTime = tMidGPS;
  /* binary orbit and higher spindowns not considered */
  thisPoint.asini = 0 /* isolated pulsar */;
  XLAL_INIT_MEM ( thisPoint.fkdot );

  /* set up some semiCoherent parameters */
  semiCohPar.useToplist = FALSE;
  semiCohPar.tsMid = midTstack;
  /* semiCohPar.refTime = tStartGPS; */
  semiCohPar.refTime = tMidGPS;
  /* calculate detector velocity and positions */
  LAL_CALL( GetStackVelPos( &status, &velStack, &posStack, Fstat_in_vec), &status);
  semiCohPar.vel = velStack;
  semiCohPar.pos = posStack;

  semiCohPar.outBaseName = uvar_fnameout;
  semiCohPar.pixelFactor = uvar_pixelFactor;
  semiCohPar.nfdot = 1;
  semiCohPar.dfdot = 1;

  /* allocate memory for Hough candidates */
  semiCohCandList.length = 1;	/* hardcoded: only keep top-candidate per Hough map */
  semiCohCandList.refTime = tMidGPS;
  semiCohCandList.nCandidates = 0; /* initialization */
  semiCohCandList.list = LALCalloc( semiCohCandList.length, sizeof(SemiCohCandidate));
  if ( semiCohCandList.list == NULL) {
    XLALPrintError( "Error allocating memory \n");
    return(HIERARCHICALSEARCH_EMEM);
  }

  /* set semicoherent patch size */
  if ( LALUserVarWasSet(&uvar_semiCohPatchX)) {
    semiCohPar.patchSizeX = uvar_semiCohPatchX;
  }
  else {
    semiCohPar.patchSizeX = 1.0 / (  tStack * usefulParams.spinRange_midTime.fkdot[0] * VEPI );
  }

  if ( LALUserVarWasSet(&uvar_semiCohPatchY)) {
    semiCohPar.patchSizeY = uvar_semiCohPatchY;
  }
  else {
    semiCohPar.patchSizeY = 1.0 / ( tStack * usefulParams.spinRange_midTime.fkdot[0] * VEPI );
  }

  LogPrintf(LOG_DETAIL,"Hough patchsize is %frad x %frad\n",
	    semiCohPar.patchSizeX, semiCohPar.patchSizeY);

  /* allocate some fstat memory */
  fstatVector.length = nStacks;
  fstatVector.data = NULL;
  fstatVector.data = (REAL4FrequencySeries *)LALCalloc( 1, nStacks * sizeof(REAL4FrequencySeries));
  if ( fstatVector.data == NULL) {
    fprintf(stderr, "error allocating memory [HierarchicalSearch.c %d]\n" , __LINE__);
    return(HIERARCHICALSEARCH_EMEM);
  }

  INT8 fBin0 = floor(usefulParams.spinRange_midTime.fkdot[0]/uvar_WU_dFreq + 0.5);
  //printf ("Freq0 = %.16g, fBin0 = %d\n", usefulParams.spinRange_midTime.fkdot[0], (INT4)fBin0 );

  /* ==================== loop over candidates ==================== */
  LogPrintf(LOG_DEBUG, "Total candidates = %d. Progress: ", InputCandList->length);

  UINT4 iCand;
  for ( iCand = 0; iCand < InputCandList->length; iCand ++ )
    {
      HoughCandidate *thisCand = &InputCandList->data[iCand];
      LogPrintfVerbatim(LOG_DEBUG, "%d, ", iCand );

      /* set coarse-grid point from followup-candidate */
      thisPoint.Alpha = thisCand->Alpha;
      thisPoint.Delta = thisCand->Delta;
      thisPoint.fkdot[0] = thisCand->Freq;
      thisPoint.fkdot[1] = thisCand->f1dot;

      LogPrintf (LOG_DETAIL, "Computing candidate %d/%d: f=%f, Alpha=%f, Delta=%f, f1dot=%g, refTime = %d ...",
                 iCand + 1, InputCandList->length, thisPoint.fkdot[0], thisPoint.Alpha, thisPoint.Delta, thisPoint.fkdot[1], thisPoint.refTime.gpsSeconds );

      /*------------- calculate F statistic for each stack --------------*/
      SkyPosition skypos;
      /* get amplitude modulation weights */
      skypos.longitude = thisPoint.Alpha;
      skypos.latitude = thisPoint.Delta;
      skypos.system = COORDINATESYSTEM_EQUATORIAL;

      /* initialize weights to unity */
      LAL_CALL( LALHOUGHInitializeWeights( &status, weightsV), &status);

      if (uvar_useWeights) {
	LAL_CALL( ComputeStackNoiseAndAMWeights( &status, weightsV, Fstat_in_vec, skypos), &status);
      }

      semiCohPar.weightsV = weightsV;
      semiCohPar.alpha = thisPoint.Alpha;
      semiCohPar.delta = thisPoint.Delta;

      LogPrintf(LOG_DETAIL, "Stack weights for alpha = %f, delta = %f are:\n", skypos.longitude, skypos.latitude);
      for (k = 0; k < nStacks; k++) {
	LogPrintf(LOG_DETAIL, "%f\n", weightsV->data[k]);
      }

      { /********Allocate fstat vector memory *****************/

	/* extra bins for fstat due to skypatch and spindowns */
	UINT4 extraBinsfdot;
	REAL8 freqHighest;

	/* calculate number of bins for fstat overhead */
	freqHighest = usefulParams.spinRange_midTime.fkdot[0] + usefulParams.spinRange_midTime.fkdotBand[0];
	ComputeNumExtraBins(&status, &semiCohPar, 0, freqHighest, uvar_WU_dFreq);

	extraBinsfdot = 0;

	semiCohPar.extraBinsFstat += extraBinsfdot;

	/* allocate fstat memory */
	binsFstatSearch = 1;
	//binsFstatSearch = (UINT4)(usefulParams.spinRange_midTime.fkdotBand[0]/uvar_WU_dFreq + 1e-6) + 1;
	binsFstat1 = binsFstatSearch + 2*semiCohPar.extraBinsFstat;

	for (k = 0; k < nStacks; k++) {
	  /* careful--the epoch here is not the reference time for f0! */
	  fstatVector.data[k].epoch = startTstack->data[k];
	  fstatVector.data[k].deltaF = uvar_WU_dFreq;
	  fstatVector.data[k].f0 = thisPoint.fkdot[0] - semiCohPar.extraBinsFstat * uvar_WU_dFreq;

	  if (fstatVector.data[k].data == NULL) {
	    fstatVector.data[k].data = (REAL4Sequence *)LALCalloc( 1, sizeof(REAL4Sequence));
	    if ( fstatVector.data[k].data == NULL) {
	      fprintf(stderr, "error allocating memory [HierarchicalSearch.c %d]\n" , __LINE__);
	      return(HIERARCHICALSEARCH_EMEM);
	    }

	    fstatVector.data[k].data->length = binsFstat1;
	    fstatVector.data[k].data->data = (REAL4 *)LALCalloc( 1, binsFstat1 * sizeof(REAL4));
	    if ( fstatVector.data[k].data->data == NULL) {
	      fprintf(stderr, "error allocating memory [HierarchicalSearch.c %d]\n" , __LINE__);
	      return(HIERARCHICALSEARCH_EMEM);
	    }
          } else {
          fstatVector.data[k].data = (REAL4Sequence *)LALRealloc( fstatVector.data[k].data, sizeof(REAL4Sequence));
          if ( fstatVector.data[k].data == NULL) {
            fprintf(stderr, "error allocating memory [HierarchicalSearch.c %d]\n" , __LINE__);
            return(HIERARCHICALSEARCH_EMEM);
          }

          fstatVector.data[k].data->length = binsFstat1;
          fstatVector.data[k].data->data = (REAL4 *)LALRealloc( fstatVector.data[k].data->data, binsFstat1 * sizeof(REAL4));
          if ( fstatVector.data[k].data->data == NULL) {
            fprintf(stderr, "error allocating memory [HierarchicalSearch.c %d]\n" , __LINE__);
            return(HIERARCHICALSEARCH_EMEM);
          }
          }
	} /* loop over stacks */


      } /* fstat memory allocation block */

      /* calculate the Fstatistic for each stack*/
      LogPrintf(LOG_DETAIL, "Starting Fstat calculation for each stack...");

      PulsarDopplerParams dummyPoint = thisPoint;
      dummyPoint.fkdot[0] = fstatVector.data[0].f0;	// unfortunately, doppler-point also needs to refer to lowest frequency in band

      for ( k = 0; k < nStacks; k++)
        {
          const int retn = XLALComputeFstat(&Fstat_res, Fstat_in_vec->data[k], &dummyPoint, fstatVector.data[0].deltaF, binsFstat1, Fstat_what);
          if ( retn != XLAL_SUCCESS ) {
            XLALPrintError ("%s: XLALComputeFstat() failed with errno=%d\n", __func__, xlalErrno );
            return xlalErrno;
          }
          for (UINT4 iFreq = 0; iFreq < binsFstat1; ++iFreq) {
            fstatVector.data[k].data->data[iFreq] = 0.5 * Fstat_res->twoF[iFreq];    // *** copy value of *1*F ***
          }
        } /* for k < nStacks */

      LogPrintfVerbatim(LOG_DETAIL, "done\n");

      /* print fstat vector if required -- mostly for debugging */
      if ( uvar_printFstat1 )
        {
          for (k = 0; k < nStacks; k++)
            LAL_CALL( PrintFstatVec ( &status, fstatVector.data + k, fpFstat1, &dummyPoint, refTimeGPS, k+1), &status);
        }


      /*--------------- get candidates from a semicoherent search ---------------*/

      /* the input to this section is the set of fstat vectors fstatVector and the
         parameters semiCohPar. The output is the list of candidates in semiCohCandList */

      /* set spindown for Hough grid -- same as for Fstat calculation */
      semiCohPar.fdot = thisPoint.fkdot[1];

      /* the hough option */
      /* select peaks */
      LogPrintf(LOG_DETAIL, "Starting Hough calculation...\n");
      sumWeightSquare = 0.0;
      for ( k = 0; k < nStacks; k++)
        sumWeightSquare += weightsV->data[k] * weightsV->data[k];

      /* set number count threshold based on significance threshold */
      meanN = nStacks * alphaPeak;
      sigmaN = sqrt(sumWeightSquare * alphaPeak * (1.0 - alphaPeak));
      LogPrintf(LOG_DETAIL, "Expected mean number count=%f, std=%f\n", meanN, sigmaN );

      /* convert fstat vector to peakgrams using the Fstat threshold */
      LAL_CALL( FstatVectToPeakGram( &status, &pgV, &fstatVector, (REAL4)uvar_peakThrF), &status);

     /* get candidates */
      /* this is the second most costly function. We here allow for using architecture-specific
         optimized functions from a local file instead of the standard RCComputeFstatHoughMap()
         below that refers to the LALHOUGH functions in LAL */
      LAL_CALL ( RCComputeFstatHoughMap ( &status, &semiCohCandList, &pgV, &semiCohPar, fBin0 ), &status);

      /* now correct the candidate frequency, which has suffered from a frequency-bin 'quantization' error in
       * the peak-gram step (which only stores fBinIni = round[f0/dFreq]). See bug #147.
       */
      if ( uvar_correctFreqs )
        {
          LogPrintf (LOG_DETAIL, "Correcting output candidate frequency: Hough f0 = %.9f, but actually f0 = %.9f (offset = %.9g)\n",
                     semiCohCandList.list[0].freq, thisPoint.fkdot[0], semiCohCandList.list[0].freq - thisPoint.fkdot[0] );
          semiCohCandList.list[0].freq = thisPoint.fkdot[0];
        }

      /* free peakgrams -- we don't need them now because we have the Hough maps */
      for (k=0; k<nStacks; k++)
        LALFree(pgV.pg[k].peak);
      LALFree(pgV.pg);
      LogPrintf(LOG_DETAIL, "...finished Hough calculation\n");

      /* end hough */

      LAL_CALL( GetSemiCohToplist(&status, semiCohToplist, &semiCohCandList, meanN, sigmaN), &status);

    } /* for InputCandList */

  LogPrintfVerbatim ( LOG_DEBUG, " done.\n");

  /* Also compute F, FX (for line veto statistics) for all candidates in final toplist */
  if ( uvar_outputFX ) {
    LogPrintfVerbatim ( LOG_DEBUG, "Computing FX ...");

    /* MultiNoiseWeightsSequence *multiNoiseWeightsPointer; */
    /* if ( uvar_useFstatWeights ) */
    /*   multiNoiseWeightsPointer = &stackMultiNoiseWeights; */
    /* else */
    /*   multiNoiseWeightsPointer = NULL; */

    xlalErrno = 0;
    XLALComputeExtraStatsForToplist ( semiCohToplist, "HoughFStat", Fstat_in_vec, usefulParams.detectorIDs,
                                      usefulParams.startTstack, refTimeGPS,  uvar_outputSingleSegStats );
    if ( xlalErrno != 0 ) {
      XLALPrintError ("%s line %d : XLALComputeExtraStatsForToplist() failed with xlalErrno = %d.\n\n", __func__, __LINE__, xlalErrno );
      return(HIERARCHICALSEARCH_EBAD);
    }
    LogPrintfVerbatim ( LOG_DEBUG, " done.\n");
  }

  LogPrintf(LOG_DETAIL, "Finished analysis and now printing results and cleaning up...");

  LogPrintfVerbatim ( LOG_DEBUG, "Writing output ...\n");

  /* print candidates */
  if (!(fpSemiCoh = fopen(fnameSemiCohCand, "wb")))
    {
      LogPrintf ( LOG_CRITICAL, "Unable to open output-file '%s' for writing.\n", fnameSemiCohCand);
      return HIERARCHICALSEARCH_EFILE;
    }
  /* write header-line comment explaining columns */
  fprintf ( fpSemiCoh, "%%%%  Freq            Alpha              Delta              f1dot                 HoughFStat        AlphaBest          DeltaBest          MeanSig   VarSig    <multiF>   <F1>    <F2> ...\n");

  sort_houghFStat_toplist(semiCohToplist);
  if ( write_houghFStat_toplist_to_fp( semiCohToplist, fpSemiCoh, NULL) < 0)
    fprintf( stderr, "Error in writing toplist to file\n");
  /*     LAL_CALL( AppendFstatCandidates( &status, &fStatCand, fpFstat), &status); */
  if (fprintf(fpSemiCoh,"%%DONE\n") < 0)
    fprintf(stderr, "Error writing end marker\n");
  fclose(fpSemiCoh);

  LogPrintfVerbatim ( LOG_DEBUG, " done.\n");

  /*------------ free all remaining memory -----------*/

  /* free memory */

 LALFree(fnameSemiCohCand);


  if ( VCSInfoString ) XLALFree ( VCSInfoString );

  if ( uvar_printFstat1 )
    {
      fclose(fpFstat1);
      LALFree( fnameFstatVec1 );
    }

  /* free first stage memory */
  XLALDestroyFstatInputVector( Fstat_in_vec );
  XLALDestroyFstatResults( Fstat_res );

  XLALDestroyTimestampVector(midTstack);
  XLALDestroyTimestampVector(startTstack);

  /* free Fstat vectors  */
  for(k = 0; k < nStacks; k++)
    if (fstatVector.data[k].data) {
      if (fstatVector.data[k].data->data)
	LALFree(fstatVector.data[k].data->data);
      LALFree(fstatVector.data[k].data);
    }
  LALFree(fstatVector.data);

  XLALDestroyStringVector ( usefulParams.detectorIDs );

  /* free Vel/Pos vectors and ephemeris */
  XLALDestroyEphemerisData(edat);
  LAL_CALL( LALDDestroyVectorSequence (&status,  &velStack), &status);
  LAL_CALL( LALDDestroyVectorSequence (&status,  &posStack), &status);

  if (weightsNoise) {
    LAL_CALL( LALDDestroyVector ( &status, &weightsNoise ), &status);
  }
  LAL_CALL( LALDDestroyVector ( &status, &weightsV ), &status);

  /* free input candidate list */
  XLALDestroyHoughCandidateList ( InputCandList );

  /* free candidates */
  LALFree(semiCohCandList.list);
  free_houghFStat_toplist(&semiCohToplist);

  LAL_CALL (LALDestroyUserVars(&status), &status);

  LALCheckMemoryLeaks();

  LogPrintfVerbatim(LOG_DETAIL, "done\n");

  return HIERARCHICALSEARCH_ENORM;
} /* main */




/**
 * Set up stacks, read SFTs, calculate SFT noise weights and calculate
 * detector-state
 */
void SetUpSFTs( LALStatus *status,			/**< pointer to LALStatus structure */
		FstatInputVector** p_Fstat_in_vec,	/**< pointer to vector of Fstat input data structures for XLALComputeFstat(), one per stack */
		UsefulStageVariables *in /**< input params */)
{

  SFTCatalog *catalog = NULL;
  static SFTConstraints constraints;
  REAL8 timebase, tObs, deltaFsft;
  UINT4 k,numSFTby2;
  LIGOTimeGPS tStartGPS, tEndGPS, refTimeGPS, tMidGPS;
  SFTCatalogSequence catalogSeq;

  REAL8 doppWings, fMin, fMax;
  REAL8 startTime_freqLo, startTime_freqHi;
  REAL8 endTime_freqLo, endTime_freqHi;
  REAL8 freqLo, freqHi;

  INT4 sft_check_result = 0;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* get sft catalog */
  constraints.minStartTime = &(in->minStartTimeGPS);
  constraints.maxStartTime = &(in->maxStartTimeGPS);
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
  TRY( SetUpStacks( status->statusPtr, &catalogSeq, in->tStack, catalog, in->nStacks), status);

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

  /* fill detector name vector with all detectors present in any data sements */
  in->detectorIDs = NULL;
  for (k = 0; k < in->nStacks; k++) {
    if ( ( in->detectorIDs = XLALGetDetectorIDsFromSFTCatalog ( in->detectorIDs, catalogSeq.data + k ) ) == NULL ) {
      ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
    }
  }

  /* get timestamps of start and mid of each stack */
  /* set up vector containing mid times of stacks */
  in->midTstack =  XLALCreateTimestampVector ( in->nStacks );

  /* set up vector containing start times of stacks */
  in->startTstack =  XLALCreateTimestampVector ( in->nStacks );

  /* now loop over stacks and get time stamps */
  for (k = 0; k < in->nStacks; k++) {

    if ( catalogSeq.data[k].length == 0 ) {
      /* something is wrong */
      ABORT ( status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
    }

    /* start time of stack = time of first sft in stack */
    in->startTstack->data[k] = catalogSeq.data[k].data[0].header.epoch;

    /* mid time of stack */
    numSFTby2 = catalogSeq.data[k].length/2;
    in->midTstack->data[k] = catalogSeq.data[k].data[numSFTby2].header.epoch;

  } /* loop over k */


  /* set reference time for pular parameters */
  /* first calculate the mid time of observation time span*/
  {
    REAL8 tStart8, tEnd8, tMid8;

    tStart8 = XLALGPSGetREAL8( &tStartGPS );
    tEnd8   = XLALGPSGetREAL8( &tEndGPS );
    tMid8 = 0.5 * (tStart8 + tEnd8);
    XLALGPSSetREAL8( &tMidGPS, tMid8 );
  }

  if ( in->refTime > 0)  {
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

  fMin = freqLo - doppWings;
  fMax = freqHi + doppWings;

  /* set up vector of Fstat input data structs */
  (*p_Fstat_in_vec) = XLALCreateFstatInputVector( in->nStacks );
  if ( (*p_Fstat_in_vec) == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

#ifdef OUTPUT_TIMING
  /* need to count the total number of SFTs */
  nStacks = in->nStacks;
  nSFTs = 0;
#endif

  /* loop over stacks and read sfts */
  for (k = 0; k < in->nStacks; k++) {

    /* create Fstat input data struct for demodulation */
    (*p_Fstat_in_vec)->data[k] = XLALCreateFstatInput_Demod( in->Dterms, FMETHOD_DEMOD_BEST );
    if ( (*p_Fstat_in_vec)->data[k] == NULL ) {
      XLALPrintError("%s: XLALCreateFstatInput_Demod() failed with errno=%d", __func__, xlalErrno);
      ABORT ( status, HIERARCHICALSEARCH_EXLAL, HIERARCHICALSEARCH_MSGEXLAL );
    }
    if ( XLALSetupFstatInput( (*p_Fstat_in_vec)->data[k], &catalogSeq.data[k], fMin, fMax, NULL,
                                  NULL, NULL, in->blocksRngMed, in->edat, in->SSBprec, 0 ) != XLAL_SUCCESS ) {
      XLALPrintError("%s: XLALSetupFstatInput() failed with errno=%d", __func__, xlalErrno);
      ABORT ( status, HIERARCHICALSEARCH_EXLAL, HIERARCHICALSEARCH_MSGEXLAL );
    }

    /* get SFT detectors and timestamps */
    const MultiLALDetector *multiIFO = XLALGetFstatInputDetectors( (*p_Fstat_in_vec)->data[k] );
    if ( multiIFO == NULL ) {
      XLALPrintError("%s: XLALGetFstatInputDetectors() failed with errno=%d", __func__, xlalErrno);
      ABORT ( status, HIERARCHICALSEARCH_EXLAL, HIERARCHICALSEARCH_MSGEXLAL );
    }
    const MultiLIGOTimeGPSVector *multiTS = XLALGetFstatInputTimestamps( (*p_Fstat_in_vec)->data[k] );
    if ( multiTS == NULL ) {
      XLALPrintError("%s: XLALGetFstatInputTimestamps() failed with errno=%d", __func__, xlalErrno);
      ABORT ( status, HIERARCHICALSEARCH_EXLAL, HIERARCHICALSEARCH_MSGEXLAL );
    }

    /* print debug info about this stack */
    LogPrintf(LOG_DETAIL, "Stack %d ", k);
    for ( UINT4 j = 0; j < multiTS->length; j++) {
      LogPrintfVerbatim(LOG_DETAIL, "%s: %d  ", multiIFO->sites[j].frDetector.prefix, multiTS->data[j]->length);
    } /* loop over ifos */
    LogPrintfVerbatim(LOG_DETAIL, "\n");

#ifdef OUTPUT_TIMING
    /* need to count the total number of SFTs */
    for ( UINT4 X = 0; X < multiTS->length; X++ ) {
      nSFTs += multiTS->data[X]->length;
    }
#endif

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



/**
 * Function for calculating Hough Maps and candidates.
 * This function takes a peakgram as input. This peakgram was constructed
 * by setting a threshold on a sequence of Fstatistic vectors.  The function
 * produces a Hough map in the sky for each value of the frequency and spindown.
 * The Hough nummber counts are then used to select candidates in
 * parameter space to be followed up in a more refined search.
 * This uses DriveHough_v3.c as a prototype suitably modified to work
 * on demodulated data instead of SFTs.
 */
void
RCComputeFstatHoughMap(LALStatus *status,		/**< pointer to LALStatus structure */
                       SemiCohCandidateList  *out,   /**< Candidates from thresholding Hough number counts */
                       HOUGHPeakGramVector *pgV, 	/**< HOUGHPeakGramVector obtained after thresholding Fstatistic vectors */
                       SemiCoherentParams *params,	/**< pointer to HoughParams -- parameters for calculating Hough maps */
                       INT8 fBin0			/**< UNDOCUMENTED */
                     )
{

  /* hough structures */
  HOUGHMapTotal ht;
  HOUGHptfLUTVector   lutV; /* the Look Up Table vector*/
  PHMDVectorSequence  phmdVS;  /* the partial Hough map derivatives */
  UINT8FrequencyIndexVector freqInd; /* for trajectory in time-freq plane */
  HOUGHResolutionPar parRes;   /* patch grid information */
  HOUGHPatchGrid  patch;   /* Patch description */
  HOUGHParamPLUT  parLut;  /* parameters needed to build lut  */
  HOUGHDemodPar   parDem;  /* demodulation parameters */
  HOUGHSizePar    parSize;

  UINT2  xSide, ySide, maxNBins, maxNBorders;
  INT8  fBinIni, fBinFin, fBin;
  INT4  iHmap, nfdot;
  UINT4 k, nStacks ;
  REAL8 deltaF, dfdot, alpha, delta;
  REAL8 patchSizeX, patchSizeY;
  REAL8VectorSequence *vel, *pos;
  REAL8 fdot, refTime;
  LIGOTimeGPS refTimeGPS;
  LIGOTimeGPSVector   *tsMid;
  REAL8Vector *timeDiffV=NULL;
  UINT8Vector hist; /* histogram vector */
  UINT8Vector histTotal; /* total histogram vector */
  HoughStats stats; /* statistics struct */
  CHAR *fileStats = NULL;
  FILE *fpStats = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* check input is not null */
  if ( out == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }
  if ( out->length == 0 ) {
    ABORT ( status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  }
  if ( out->list == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  }
  if ( pgV == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }
  if ( pgV->length == 0 ) {
    ABORT ( status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  }
  if ( pgV->pg == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }
  if ( params == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }



  /* copy some parameters from peakgram vector */
  deltaF = pgV->pg->deltaF;

  nStacks = pgV->length;
  fBinIni = pgV->pg[0].fBinIni;
  fBinFin = pgV->pg[0].fBinFin;

  /* copy some params to local variables */
  nfdot = params->nfdot;
  dfdot = params->dfdot;
  alpha = params->alpha;
  delta = params->delta;
  vel = params->vel;
  pos = params->pos;
  fdot = params->fdot;
  tsMid = params->tsMid;
  refTimeGPS = params->refTime;
  refTime = XLALGPSGetREAL8(&refTimeGPS);

  /* set patch size */
  /* this is supposed to be the "educated guess"
     delta theta = 1.0 / (Tcoh * f0 * Vepi )
     where Tcoh is coherent time baseline,
     f0 is frequency and Vepi is rotational velocity
     of detector */
  patchSizeX = params->patchSizeX;
  patchSizeY = params->patchSizeY;

  /* calculate time differences from start of observation time for each stack*/
  TRY( LALDCreateVector( status->statusPtr, &timeDiffV, nStacks), status);

  for (k=0; k<nStacks; k++) {
    REAL8 tMidStack;
    tMidStack = XLALGPSGetREAL8(tsMid->data + k);
    timeDiffV->data[k] = tMidStack - refTime;
  }



  /*--------------- first memory allocation --------------*/
  /* look up table vector */
  lutV.length = nStacks;
  lutV.lut = NULL;
  lutV.lut = (HOUGHptfLUT *)LALCalloc(1,nStacks*sizeof(HOUGHptfLUT));
  if ( lutV.lut == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }


  /* partial hough map derivative vector */
  phmdVS.length  = nStacks;

  {
    REAL8 maxTimeDiff, startTimeDiff, endTimeDiff;

    startTimeDiff = fabs(timeDiffV->data[0]);
    endTimeDiff = fabs(timeDiffV->data[timeDiffV->length - 1]);
    maxTimeDiff = HSMAX( startTimeDiff, endTimeDiff);

    /* set number of freq. bins for which LUTs will be calculated */
    /* this sets the range of residual spindowns values */
    /* phmdVS.nfSize  = 2*nfdotBy2 + 1; */
    phmdVS.nfSize  = 2 * floor((nfdot-1) * (REAL4)(dfdot * maxTimeDiff / deltaF) + 0.5f) + 1;
  }

  phmdVS.deltaF  = deltaF;
  phmdVS.phmd = NULL;
  phmdVS.phmd=(HOUGHphmd *)LALCalloc( 1,phmdVS.length * phmdVS.nfSize *sizeof(HOUGHphmd));
  if ( phmdVS.phmd == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }

  /* residual spindown trajectory */
  freqInd.deltaF = deltaF;
  freqInd.length = nStacks;
  freqInd.data = NULL;
  freqInd.data =  ( UINT8 *)LALCalloc(1,nStacks*sizeof(UINT8));
  if ( freqInd.data == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }

  /* resolution in space of residual spindowns */
  ht.dFdot.length = 1;
  ht.dFdot.data = NULL;
  ht.dFdot.data = (REAL8 *)LALCalloc( 1, ht.dFdot.length * sizeof(REAL8));
  if ( ht.dFdot.data == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }

  /* the residual spindowns */
  ht.spinRes.length = 1;
  ht.spinRes.data = NULL;
  ht.spinRes.data = (REAL8 *)LALCalloc( 1, ht.spinRes.length*sizeof(REAL8));
  if ( ht.spinRes.data == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }

  /* the residual spindowns */
  ht.spinDem.length = 1;
  ht.spinDem.data = NULL;
  ht.spinDem.data = (REAL8 *)LALCalloc( 1, ht.spinRes.length*sizeof(REAL8));
  if ( ht.spinDem.data == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }

  /* the demodulation params */
  parDem.deltaF = deltaF;
  parDem.skyPatch.alpha = alpha;
  parDem.skyPatch.delta = delta;
  parDem.spin.length = 1;
  parDem.spin.data = NULL;
  parDem.spin.data = (REAL8 *)LALCalloc(1, sizeof(REAL8));
  if ( parDem.spin.data == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }
  parDem.spin.data[0] = fdot;

  /* the skygrid resolution params */
  parRes.deltaF = deltaF;
  parRes.patchSkySizeX  = patchSizeX;
  parRes.patchSkySizeY  = patchSizeY;
  parRes.pixelFactor = params->pixelFactor;
  parRes.pixErr = PIXERR;
  parRes.linErr = LINERR;
  parRes.vTotC = VTOT;

  /* memory allocation for histogram and opening stats file*/
  if ( uvar_printStats ) {
    hist.length = nStacks+1;
    histTotal.length = nStacks+1;
    hist.data = NULL;
    histTotal.data = NULL;
    hist.data = (UINT8 *)LALCalloc(1, (nStacks+1)*sizeof(UINT8));
    if ( hist.data == NULL ) {
      ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
    }

    histTotal.data = (UINT8 *)LALCalloc(1, (nStacks+1)*sizeof(UINT8));
    if ( histTotal.data == NULL ) {
      ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
    }

    {
      UINT4   j;
      for(j=0; j< histTotal.length; ++j)
	histTotal.data[j]=0;
    }
    {
      const CHAR *append = "stats";
      fileStats = LALCalloc(strlen(params->outBaseName) + strlen(append) + 1, sizeof(CHAR) );
      if ( fileStats == NULL ) {
	ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
      }

      strcpy( fileStats, params->outBaseName);
      strcat( fileStats, append);
    }
    if ( !(fpStats = fopen(fileStats, "wb")))
      {
	fprintf(stderr, "Unable to open file '%s' for writing...continuing\n", fileStats);
      }
  }



  /* adjust fBinIni and fBinFin to take maxNBins into account */
  /* and make sure that we have fstat values for sufficient number of bins */
  parRes.f0Bin =  fBinIni;

  fBinIni += params->extraBinsFstat;
  fBinFin -= params->extraBinsFstat;

  LogPrintf (LOG_DETAIL, "fBinIni = %d,   fBinFin = %d, extraBins = %d\n", (INT4)fBinIni, (INT4)fBinFin, params->extraBinsFstat );

  /* this is not very clean -- the Fstat calculation has to know how many extra bins are needed */

  LogPrintf(LOG_DETAIL, "Freq. range analyzed by Hough = [%fHz - %fHz] (%d bins)\n", fBinIni*deltaF, fBinFin*deltaF, fBinFin - fBinIni + 1);
  ASSERT ( fBinIni <= fBinFin, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );

  /* initialise number of candidates -- this means that any previous candidates
     stored in the list will be lost for all practical purposes*/
  out->nCandidates = 0;

  INT4 numHmaps = (fBinFin - fBinIni + 1)*phmdVS.nfSize;
  LogPrintf (LOG_DETAIL, "numHmap = %d\n", numHmaps );
  if (out->length != numHmaps) {
    out->length = numHmaps;
    out->list = (SemiCohCandidate *)LALRealloc( out->list, out->length * sizeof(SemiCohCandidate));
    if ( out->list == NULL ) {
      ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
    }
  }

  /*------------------ start main Hough calculation ---------------------*/

  /* initialization */
  fBin= fBinIni; /* initial search bin */
  iHmap = 0; /* hough map index */

  while( fBin <= fBinFin ){
    INT8 fBinSearch, fBinSearchMax;
    UINT4 i,j;
    REAL8UnitPolarCoor sourceLocation;

    parRes.f0Bin =  fBin0;	// *fix* this value to first bin if at WU start search-frequency
    TRY( LALHOUGHComputeSizePar( status->statusPtr, &parSize, &parRes ),  status );
    xSide = parSize.xSide;
    ySide = parSize.ySide;

    maxNBins = parSize.maxNBins;
    maxNBorders = parSize.maxNBorders;

    /*------------------ create patch grid at fBin ----------------------*/
    patch.xSide = xSide;
    patch.ySide = ySide;
    patch.xCoor = NULL;
    patch.yCoor = NULL;
    patch.xCoor = (REAL8 *)LALCalloc(1,xSide*sizeof(REAL8));
    if ( patch.xCoor == NULL ) {
      ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
    }

    patch.yCoor = (REAL8 *)LALCalloc(1,ySide*sizeof(REAL8));
    if ( patch.yCoor == NULL ) {
      ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
    }
    TRY( LALHOUGHFillPatchGrid( status->statusPtr, &patch, &parSize ), status );

    /*------------- other memory allocation and settings----------------- */
    for(j=0; j<lutV.length; ++j){
      lutV.lut[j].maxNBins = maxNBins;
      lutV.lut[j].maxNBorders = maxNBorders;
      lutV.lut[j].border = (HOUGHBorder *)LALCalloc(1,maxNBorders*sizeof(HOUGHBorder));
      if ( lutV.lut[j].border == NULL ) {
	ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
      }

      lutV.lut[j].bin =	(HOUGHBin2Border *)LALCalloc(1,maxNBins*sizeof(HOUGHBin2Border));
      if ( lutV.lut[j].bin == NULL ) {
	ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
      }

      for (i=0; i<maxNBorders; ++i){
	lutV.lut[j].border[i].ySide = ySide;
	lutV.lut[j].border[i].xPixel = (COORType *)LALCalloc(1,ySide*sizeof(COORType));
	if ( lutV.lut[j].border[i].xPixel == NULL ) {
	  ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
	}
      }
    }

    for(j = 0; j < phmdVS.length * phmdVS.nfSize; ++j){
      phmdVS.phmd[j].maxNBorders = maxNBorders;
      phmdVS.phmd[j].leftBorderP = (HOUGHBorder **)LALCalloc(1,maxNBorders*sizeof(HOUGHBorder *));
      if ( phmdVS.phmd[j].leftBorderP == NULL ) {
	ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
      }

      phmdVS.phmd[j].rightBorderP = (HOUGHBorder **)LALCalloc(1,maxNBorders*sizeof(HOUGHBorder *));
      if ( phmdVS.phmd[j].rightBorderP == NULL ) {
	ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
      }

      phmdVS.phmd[j].ySide = ySide;
      phmdVS.phmd[j].firstColumn = NULL;
      phmdVS.phmd[j].firstColumn = (UCHAR *)LALCalloc(1,ySide*sizeof(UCHAR));
      if ( phmdVS.phmd[j].firstColumn == NULL ) {
	ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
      }
    }

    /*------------------- create all the LUTs at fBin ---------------------*/
    for (j=0; j < (UINT4)nStacks; j++){  /* create all the LUTs */
      parDem.veloC.x = vel->data[3*j];
      parDem.veloC.y = vel->data[3*j + 1];
      parDem.veloC.z = vel->data[3*j + 2];
      parDem.positC.x = pos->data[3*j];
      parDem.positC.y = pos->data[3*j + 1];
      parDem.positC.z = pos->data[3*j + 2];
      parDem.timeDiff = timeDiffV->data[j];

      /* calculate parameters needed for buiding the LUT */
      TRY( LALHOUGHCalcParamPLUT( status->statusPtr, &parLut, &parSize, &parDem), status);

      /* build the LUT */
      TRY( LALHOUGHConstructPLUT( status->statusPtr, &(lutV.lut[j]), &patch, &parLut ), status);


      /* for debugging
	 fprintf(stdout,"%d\n", lutV.lut[j].nBin);
      */

      /* for debugging */
      if ( uvar_validateLUT) {
	TRY( ValidateHoughLUT( status->statusPtr, &(lutV.lut[j]), &patch, params->outBaseName, j, alpha, delta, params->weightsV->data[j]), status);
      }

      /* for debugging */
      if ( uvar_dumpLUT) {
	TRY( DumpLUT2file( status->statusPtr, &(lutV.lut[j]), &patch, params->outBaseName, j), status);
      }
    }

    /*--------- build the set of  PHMD centered around fBin -------------*/
    phmdVS.fBinMin = fBin - phmdVS.nfSize/2;
    TRY( LALHOUGHConstructSpacePHMD(status->statusPtr, &phmdVS, pgV, &lutV), status );
    TRY( LALHOUGHWeighSpacePHMD(status->statusPtr, &phmdVS, params->weightsV), status);

    /*-------------- initializing the Total Hough map space ------------*/
    ht.xSide = xSide;
    ht.ySide = ySide;
    ht.skyPatch.alpha = alpha;
    ht.skyPatch.delta = delta;
    ht.mObsCoh = nStacks;
    ht.deltaF = deltaF;
    ht.spinDem.data[0] = fdot;
    ht.patchSizeX = patchSizeX;
    ht.patchSizeY = patchSizeY;
    ht.dFdot.data[0] = dfdot;
    ht.map   = NULL;
    ht.map   = (HoughTT *)LALCalloc(1,xSide*ySide*sizeof(HoughTT));
    if ( ht.map == NULL ) {
      ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
    }

    TRY( LALHOUGHInitializeHT( status->statusPtr, &ht, &patch), status); /*not needed */

    /*  Search frequency interval possible using the same LUTs */
    fBinSearch = fBin;
    fBinSearchMax = fBin + parSize.nFreqValid - 1;

    /* Study all possible frequencies with one set of LUT */
    while ( (fBinSearch <= fBinFin) && (fBinSearch < fBinSearchMax) )  {

      /* finally we can construct the hough maps and select candidates */
      {
	INT4   n, nfdotBy2;

	nfdotBy2 = nfdot/2;
	ht.f0Bin = fBinSearch;

	/*loop over all values of residual spindown */
	/* check limits of loop */
	for( n = -nfdotBy2; n <= nfdotBy2 ; n++ ){

	  ht.spinRes.data[0] =  n*dfdot;

	  for (j=0; j < (UINT4)nStacks; j++) {
	    freqInd.data[j] = fBinSearch + floor( (REAL4)(timeDiffV->data[j]*n*dfdot/deltaF) + 0.5f);
	  }

	  TRY( LALHOUGHConstructHMT_W(status->statusPtr, &ht, &freqInd, &phmdVS),status );

	  /* get top candidate from Hough skypatch */
          SemiCohCandidate XLAL_INIT_DECL(topCand);
          TRY( GetHoughPatchTopCandidate ( status->statusPtr, &topCand, &ht, &patch, &parDem ), status);

          /* append this to candidate list */
          INT4 numCandidates = out->nCandidates;

          if (numCandidates >= out->length)           /* extend list if necessary */
            {
              out->length += BLOCKSIZE_REALLOC;
              out->list = LALRealloc( out->list, out->length * sizeof(SemiCohCandidate));
              LogPrintf(LOG_DETAIL, "Need to realloc Hough candidate list to %d entries\n", out->length);
            }
          out->list[numCandidates] = topCand;
          numCandidates++;
          out->nCandidates = numCandidates;

	  /* calculate statistics and histogram */
	  if ( uvar_printStats && (fpStats != NULL) ) {
	    TRY( LALHoughStatistics ( status->statusPtr, &stats, &ht), status );
	    TRY( LALStereo2SkyLocation ( status->statusPtr, &sourceLocation,
					stats.maxIndex[0], stats.maxIndex[1],
					&patch, &parDem), status);

	    fprintf(fpStats, "%d %f %f %f %f %f %f %f %g \n", iHmap, sourceLocation.alpha, sourceLocation.delta,
		    (REAL4)stats.maxCount, (REAL4)stats.minCount, (REAL4)stats.avgCount, (REAL4)stats.stdDev,
		    fBinSearch*deltaF,  ht.spinRes.data[0] );

	    TRY( LALHoughHistogram ( status->statusPtr, &hist, &ht), status);
	    for(j=0; j< histTotal.length; ++j)
	      histTotal.data[j]+=hist.data[j];
	  }

	  /* print hough map */
	  if ( uvar_printMaps ) {
	    TRY( PrintHmap2file( status->statusPtr, &ht, params->outBaseName, iHmap), status);
	  }

	  if ( uvar_printGrid ) {
	    /* just print one grid */
	    if ( iHmap == 0 ) {
	      TRY( PrintHoughGrid( status->statusPtr, &patch, &parDem, params->outBaseName, iHmap), status);
	    }
	  }

	  /* increment hough map index */
	  ++iHmap;

	} /* end loop over spindown trajectories */

      } /* end of block for calculating total hough maps */


      /*------ shift the search freq. & PHMD structure 1 freq.bin -------*/
      ++fBinSearch;
      TRY( LALHOUGHupdateSpacePHMDup(status->statusPtr, &phmdVS, pgV, &lutV), status );
      TRY( LALHOUGHWeighSpacePHMD(status->statusPtr, &phmdVS, params->weightsV), status);

    }   /* closing while loop over fBinSearch */

#ifdef OUTPUT_TIMING
    /* printf ("xside x yside = %d x %d = %d\n", parSize.xSide, parSize.ySide, parSize.xSide * parSize.ySide ); */
    nSkyRefine = parSize.xSide * parSize.ySide;
#endif

    fBin = fBinSearch;

    /*--------------  Free partial memory -----------------*/
    LALFree(patch.xCoor);
    LALFree(patch.yCoor);
    LALFree(ht.map);

    for (j=0; j<lutV.length ; ++j){
      for (i=0; i<maxNBorders; ++i){
	LALFree( lutV.lut[j].border[i].xPixel);
      }
      LALFree( lutV.lut[j].border);
      LALFree( lutV.lut[j].bin);
    }
    for(j=0; j<phmdVS.length * phmdVS.nfSize; ++j){
      LALFree( phmdVS.phmd[j].leftBorderP);
      LALFree( phmdVS.phmd[j].rightBorderP);
      LALFree( phmdVS.phmd[j].firstColumn);
    }

  } /* closing first while */


  /* free remaining memory */
  LALFree(ht.spinRes.data);
  LALFree(ht.spinDem.data);
  LALFree(ht.dFdot.data);
  LALFree(lutV.lut);
  LALFree(phmdVS.phmd);
  LALFree(freqInd.data);
  LALFree(parDem.spin.data);

  TRY( LALDDestroyVector( status->statusPtr, &timeDiffV), status);

  /* copy toplist candidates to output structure if necessary */
  if (uvar_printStats ) {

    /* print the histogram */
    TRY( PrintHoughHistogram(status->statusPtr, &histTotal, params->outBaseName), status);

    /* close stats file */
    LALFree(fileStats);
    fclose(fpStats);

    /* free histograms */
    LALFree(hist.data);
    LALFree(histTotal.data);
  }

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* RCComputeFstatHoughMap() */


/**
 * Function for selecting frequency bins from a set of Fstatistic vectors.
 *
 * Input is a vector of Fstatistic vectors.  It allocates memory
 * for the peakgrams based on the frequency span of the Fstatistic vectors
 * and fills tyem up by setting a threshold on the Fstatistic.  Peakgram must be
 * deallocated outside the function.
 */
void FstatVectToPeakGram (LALStatus *status,			/**< pointer to LALStatus structure */
			  HOUGHPeakGramVector *pgV,		/**< a vector of peakgrams  */
                          REAL4FrequencySeriesVector *FstatVect,/**< sequence of Fstatistic vectors */
			  REAL4  thr				/**< REAL8 threshold for selecting frequency bins */
                          )
{
  INT4 j, k;
  INT4 nStacks, nSearchBins, nPeaks;
  UCHAR *upg;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  if ( FstatVect == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }
  if ( FstatVect->length == 0 ) {
    ABORT ( status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  }
  if ( FstatVect->data == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }
  if ( pgV == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }

  nStacks = FstatVect->length;
  nSearchBins = FstatVect->data->data->length;


  /* first memory allocation */
  pgV->length = nStacks;
  pgV->pg = (HOUGHPeakGram *)LALCalloc( 1, nStacks * sizeof(HOUGHPeakGram));
  if ( pgV->pg == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }

  upg = (UCHAR *)LALCalloc( 1, nSearchBins * sizeof(UCHAR));
  if ( upg == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }

  /* loop over each stack and set peakgram */
  for (k=0; k<nStacks; k++) {
    INT4 *pInt; /* temporary pointer */
    REAL4 *pV;  /* temporary pointer */
    REAL8 f0, deltaF;
    pV = FstatVect->data[k].data->data;

    /* loop over Fstat vector, count peaks, and set upg values */
    nPeaks = 0;
    for(j=0; j<nSearchBins; j++) {
      if (pV[j] > thr ) {
	nPeaks++;
	upg[j] = 1;
      }
      else
	upg[j] = 0;
    }

    /* fix length of peakgram and allocate memory appropriately */
    pgV->pg[k].length = nPeaks;
    pgV->pg[k].peak = (INT4 *)LALCalloc( 1, nPeaks * sizeof(INT4));
    if ( pgV->pg[k].peak == NULL ) {
      ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
    }


    /* fill up other peakgram parameters */
    f0     = FstatVect->data[k].f0;
    deltaF = FstatVect->data[k].deltaF;

    pgV->pg[k].deltaF = deltaF;

    pgV->pg[k].fBinIni = (UINT4)( f0/deltaF + 0.5);
    pgV->pg[k].fBinFin = pgV->pg[k].fBinIni + nSearchBins - 1;

    /* do loop again and fill peakgram vector */
    pInt = pgV->pg[k].peak;
    for (j=0; j<nSearchBins; j++) {
      if ( upg[j] == 1) {
	*pInt = j;
	pInt++;
      }
    }
  }

  /* free the UCHAR peakgram */
  LALFree(upg);

  DETATCHSTATUSPTR (status);
  RETURN(status);
}


/**
 * \brief Breaks up input sft catalog into specified number of stacks
 *
 * Loops over elements of the catalog, assigns a bin index and
 * allocates memory to the output catalog sequence appropriately.  If
 * there are long gaps in the data, then some of the catalogs in the
 * output catalog sequence may be of zero length.
 */
void SetUpStacks(LALStatus *status, 	   /**< pointer to LALStatus structure */
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


/** Print single Hough map to a specified output file */
void PrintHmap2file(LALStatus *status,
		    HOUGHMapTotal *ht,
		    CHAR *fnameOut,
		    INT4 iHmap)
{
  FILE  *fp=NULL;   /* Output file */
  CHAR filename[256], filenumber[16];
  INT4  k, i ;
  UINT2 xSide, ySide;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  strcpy(  filename, fnameOut);
  sprintf( filenumber, ".%06d",iHmap);
  strcat(  filename, filenumber);

  fp=fopen(filename,"wb");
  ASSERT ( fp != NULL, status, HIERARCHICALSEARCH_EFILE, HIERARCHICALSEARCH_MSGEFILE );

  ySide= ht->ySide;
  xSide= ht->xSide;

  for(k=ySide-1; k>=0; --k){
    for(i=0;i<xSide;++i){
      fprintf( fp ," %f", ht->map[k*xSide +i]);
      fflush( fp );
    }
    fprintf( fp ," \n");
    fflush( fp );
  }

  fclose( fp );

  DETATCHSTATUSPTR (status);
  RETURN(status);
}


void PrintHoughGrid(LALStatus *status,
		    HOUGHPatchGrid *patch,
		    HOUGHDemodPar  *parDem,
		    CHAR *fnameOut,
		    INT4 iHmap)
{

  UINT2 xSide, ySide;
  FILE  *fp1=NULL;   /* Output file */
  FILE  *fp2=NULL;   /* Output file */
  CHAR filename1[256], filename2[256], filenumber[16];
  INT4  k, i ;
  REAL8UnitPolarCoor sourceLocation;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  sprintf( filenumber, ".%06d",iHmap);

  strcpy( filename1, fnameOut);
  strcat( filename1, "_GridAlpha");
  strcat( filename1, filenumber);

  strcpy( filename2, fnameOut);
  strcat( filename2, "_GridDelta");
  strcat( filename2, filenumber);

  fp1=fopen(filename1,"wb");
  ASSERT ( fp1 != NULL, status, HIERARCHICALSEARCH_EFILE, HIERARCHICALSEARCH_MSGEFILE );

  fp2=fopen(filename2,"wb");
  ASSERT ( fp2 != NULL, status, HIERARCHICALSEARCH_EFILE, HIERARCHICALSEARCH_MSGEFILE );

  xSide = patch->xSide;
  ySide = patch->ySide;

  for(k=ySide-1; k>=0; --k){
    for(i=0;i<xSide;++i){

      TRY( LALStereo2SkyLocation ( status->statusPtr, &sourceLocation,
				   k, i, patch, parDem), status);

      fprintf( fp1 ," %f", sourceLocation.alpha);
      fflush( fp1 );

      fprintf( fp2 ," %f", sourceLocation.delta);
      fflush( fp2 );
    }

    fprintf( fp1 ," \n");
    fflush( fp1 );

    fprintf( fp2 ," \n");
    fflush( fp2 );
  }

  fclose( fp1 );
  fclose( fp2 );

  DETATCHSTATUSPTR (status);
  RETURN(status);

}

void ValidateHoughLUT(LALStatus       *status,
		      HOUGHptfLUT     *lut,
		      HOUGHPatchGrid  *patch,
		      CHAR            *basename,
		      INT4            ind,
		      REAL4           alpha,
		      REAL4           delta,
		      REAL4           weight)
{
  FILE  *fp=NULL;
  CHAR filename[256];
  INT4  j, i;
  UINT4 k;
  INT8  f0Bin;
  UINT2 xSide, ySide, maxNBins, maxNBorders;

  HOUGHPeakGram  pg;
  HOUGHphmd      phmd;
  HOUGHMapDeriv  hd;
  HOUGHMapTotal  ht;

  BOOLEAN validateFlag = FALSE;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  strcpy(  filename, basename);
  strcat(  filename, ".validate");

  maxNBins = lut->maxNBins;
  maxNBorders = lut->maxNBorders;
  f0Bin = lut->f0Bin;

  xSide = patch->xSide;
  ySide = patch->ySide;

  pg.deltaF = lut->deltaF;
  pg.fBinIni = f0Bin - maxNBins;
  pg.fBinFin = f0Bin + 5*maxNBins;
  pg.length = maxNBins;
  pg.peak = NULL;
  pg.peak = (INT4 *)LALCalloc(1, pg.length*sizeof(INT4));

  phmd.fBin = f0Bin;
  phmd.maxNBorders = maxNBorders;
  phmd.leftBorderP = (HOUGHBorder **)LALMalloc(maxNBorders*sizeof(HOUGHBorder *));
  phmd.rightBorderP =  (HOUGHBorder **)LALMalloc(maxNBorders*sizeof(HOUGHBorder *));

  phmd.ySide = ySide;
  phmd.firstColumn = NULL;
  phmd.firstColumn = (UCHAR *)LALMalloc(ySide*sizeof(UCHAR));

  ht.xSide = xSide;
  ht.ySide = ySide;
  ht.map   = (HoughTT *)LALMalloc(xSide*ySide*sizeof(HoughTT));

  hd.xSide = xSide;
  hd.ySide = ySide;
  hd.map   = (HoughDT *)LALMalloc((xSide+1)*ySide*sizeof(HoughDT));

  /* construct a fake peakgram to print all borders
     -- peakgram should be 1,0,1,0,1,0.... */
  for (k = 0; k < pg.length; ++k){
    pg.peak[k] = 2*k;
  }

  TRY( LALHOUGHPeak2PHMD(status->statusPtr, &phmd, lut, &pg ), status );

  TRY( LALHOUGHInitializeHT(status->statusPtr, &ht, patch ), status );

  TRY( LALHOUGHInitializeHD(status->statusPtr, &hd), status );

  TRY( LALHOUGHAddPHMD2HD(status->statusPtr, &hd, &phmd ), status );

  TRY( LALHOUGHIntegrHD2HT(status->statusPtr, &ht, &hd ), status );

  for(j = ySide-1; j >= 0;  --j){
    for(i = 0; i < xSide; ++i){

      if (( ht.map[j*xSide+i] > 1.0) || (ht.map[j*xSide+i] < 0.0 ))
	validateFlag = TRUE;
    }
  }

  LALFree(pg.peak);

  LALFree(phmd.leftBorderP);
  LALFree(phmd.rightBorderP);
  LALFree(phmd.firstColumn);

  LALFree(ht.map);
  LALFree(hd.map);

  if (validateFlag) {
    fp=fopen(filename,"a");
    ASSERT ( fp != NULL, status, HIERARCHICALSEARCH_EFILE, HIERARCHICALSEARCH_MSGEFILE );
    fprintf(fp ," %d  %f  %f  %f\n", ind, alpha, delta, weight);
    fflush(fp);
    fclose(fp);
  }

  DETATCHSTATUSPTR (status);
  RETURN(status);
}


/** Print single Hough map to a specified output file */
void DumpLUT2file(LALStatus       *status,
		  HOUGHptfLUT     *lut,
		  HOUGHPatchGrid  *patch,
		  CHAR            *basename,
		  INT4            ind)
{

  FILE  *fp=NULL;
  CHAR filename[256], filenumber[16];
  INT4  j, i;
  UINT4 k;
  INT8  f0Bin;
  UINT2 xSide, ySide, maxNBins, maxNBorders;

  HOUGHPeakGram  pg;
  HOUGHphmd      phmd;
  HOUGHMapDeriv  hd;
  HOUGHMapTotal  ht;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  strcpy(  filename, basename);
  strcat(  filename, ".lut");
  sprintf( filenumber, ".%06d",ind);
  strcat(  filename, filenumber);

  fp=fopen(filename,"w");
  ASSERT ( fp != NULL, status, HIERARCHICALSEARCH_EFILE, HIERARCHICALSEARCH_MSGEFILE );

  maxNBins = lut->maxNBins;
  maxNBorders = lut->maxNBorders;
  f0Bin = lut->f0Bin;

  xSide = patch->xSide;
  ySide = patch->ySide;

  pg.deltaF = lut->deltaF;
  pg.fBinIni = f0Bin - maxNBins;
  pg.fBinFin = f0Bin + 5*maxNBins;
  pg.length = maxNBins;
  pg.peak = NULL;
  pg.peak = (INT4 *)LALCalloc(1, pg.length*sizeof(INT4));

  phmd.fBin = f0Bin;
  phmd.maxNBorders = maxNBorders;
  phmd.leftBorderP = (HOUGHBorder **)LALMalloc(maxNBorders*sizeof(HOUGHBorder *));
  phmd.rightBorderP =  (HOUGHBorder **)LALMalloc(maxNBorders*sizeof(HOUGHBorder *));

  phmd.ySide = ySide;
  phmd.firstColumn = NULL;
  phmd.firstColumn = (UCHAR *)LALMalloc(ySide*sizeof(UCHAR));

  ht.xSide = xSide;
  ht.ySide = ySide;
  ht.map   = (HoughTT *)LALMalloc(xSide*ySide*sizeof(HoughTT));

  hd.xSide = xSide;
  hd.ySide = ySide;
  hd.map   = (HoughDT *)LALMalloc((xSide+1)*ySide*sizeof(HoughDT));

  /* construct a fake peakgram to print all borders
     -- peakgram should be 1,0,1,0,1,0.... */
  for (k = 0; k < pg.length; ++k){
   /*  pg.peak[k] = 4*k+3; */
    pg.peak[k] = 2*k;
  }

  TRY( LALHOUGHPeak2PHMD(status->statusPtr, &phmd, lut, &pg ), status );

  TRY( LALHOUGHInitializeHT(status->statusPtr, &ht, patch ), status );

  TRY( LALHOUGHInitializeHD(status->statusPtr, &hd), status );

  TRY( LALHOUGHAddPHMD2HD(status->statusPtr, &hd, &phmd ), status );

  TRY( LALHOUGHIntegrHD2HT(status->statusPtr, &ht, &hd ), status );

  for(j = ySide-1; j >= 0;  --j){
    for(i = 0; i < xSide; ++i){
      fprintf(fp ," %f", ht.map[j*xSide +i]);
      fflush(fp);
    }
    fprintf(fp," \n");
    fflush(fp);
  }

  LALFree(pg.peak);

  LALFree(phmd.leftBorderP);
  LALFree(phmd.rightBorderP);
  LALFree(phmd.firstColumn);

  LALFree(ht.map);
  LALFree(hd.map);

  fclose(fp);

  DETATCHSTATUSPTR (status);
  RETURN(status);
}



/** Get Hough candidates as a toplist */
void GetHoughCandidates_toplist(LALStatus *status,
				toplist_t *list,
				HOUGHMapTotal *ht,
				HOUGHPatchGrid  *patch,
				HOUGHDemodPar   *parDem)
{
  REAL8UnitPolarCoor sourceLocation;
  REAL8 deltaF, f0, fdot, dFdot, patchSizeX, patchSizeY;
  INT8 f0Bin;
  INT4 i,j, xSide, ySide;
  SemiCohCandidate thisCandidate;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  deltaF = ht->deltaF;
  f0Bin = ht->f0Bin;
  f0 = f0Bin * deltaF;

  fdot = ht->spinDem.data[0] + ht->spinRes.data[0];
  dFdot = ht->dFdot.data[0];

  xSide = ht->xSide;
  ySide = ht->ySide;
  patchSizeX = ht->patchSizeX;
  patchSizeY = ht->patchSizeY;

  thisCandidate.freq =  f0;
  thisCandidate.dFreq = deltaF;
  thisCandidate.fdot = fdot;
  thisCandidate.dFdot = dFdot;
  thisCandidate.dAlpha = patchSizeX / ((REAL8)xSide);
  thisCandidate.dDelta = patchSizeY / ((REAL8)ySide);

  for (i = 0; i < ySide; i++)
    {
      for (j = 0; j < xSide; j++)
	{

	  /* get sky location of pixel */
	  TRY( LALStereo2SkyLocation (status->statusPtr, &sourceLocation,
				      j, i, patch, parDem), status);

	  thisCandidate.alpha = sourceLocation.alpha;
	  thisCandidate.delta = sourceLocation.delta;
	  thisCandidate.significance =  ht->map[i*xSide + j];

	  insert_into_toplist(list, &thisCandidate);

	} /* end loop over xSide */

    } /* end loop over ySide */

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* end hough toplist selection */



/** Get top Hough candidates over hough sky-patch */
void
GetHoughPatchTopCandidate (LALStatus            *status,
                      SemiCohCandidate     *topCand,
                      HOUGHMapTotal        *ht,
                      HOUGHPatchGrid       *patch,
                      HOUGHDemodPar        *parDem )
{
  REAL8UnitPolarCoor sourceLocationBest;
  REAL8 deltaF, f0, fdot, dFdot, patchSizeX, patchSizeY;
  INT8 f0Bin;
  INT4 i,j, xSide, ySide;
  SemiCohCandidate XLAL_INIT_DECL(thisCandidate);

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( topCand != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

  ASSERT ( ht != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( patch != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( parDem != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

  deltaF = ht->deltaF;
  f0Bin = ht->f0Bin;
  f0 = f0Bin * deltaF;	// this is actually wrong by bin-offset quantization of round(freq/dFreq)

  fdot = ht->spinDem.data[0] + ht->spinRes.data[0];
  dFdot = ht->dFdot.data[0];

  xSide = ht->xSide;
  ySide = ht->ySide;
  patchSizeX = ht->patchSizeX;
  patchSizeY = ht->patchSizeY;

  REAL8 currentMax;
  INT4 jMax, iMax;
  REAL8 meanSig, varianceSig;

  currentMax = 0.0;
  jMax = iMax = 0;

  /* loop over hough map to get location of max */
  for (i = 0; i < ySide; i++)
    {
      for (j = 0; j < xSide; j++)
        {

          if ( ht->map[i*xSide + j] > currentMax ) {
            currentMax = ht->map[i*xSide + j];
            jMax = j;
            iMax = i;
          }

        } /* end loop over xSide */

    } /* end loop over ySide */


  thisCandidate.significance =  ht->map[iMax*xSide + jMax];

  thisCandidate.freq =  f0;
  thisCandidate.dFreq = deltaF;
  thisCandidate.fdot = fdot;
  thisCandidate.dFdot = dFdot;
  thisCandidate.dAlpha = 3.0 * patchSizeX / ((REAL8)xSide);
  thisCandidate.dDelta = 3.0 * patchSizeY / ((REAL8)ySide);

  TRY( LALStereo2SkyLocation (status->statusPtr, &sourceLocationBest, jMax, iMax, patch, parDem), status);

  thisCandidate.alphaBest = sourceLocationBest.alpha;
  thisCandidate.deltaBest = sourceLocationBest.delta;

  TRY( LALHoughmapMeanVariance( status->statusPtr, &meanSig, &varianceSig, ht), status);
  thisCandidate.meanSig = meanSig;
  thisCandidate.varianceSig = varianceSig;

  thisCandidate.alpha = parDem->skyPatch.alpha;
  thisCandidate.delta = parDem->skyPatch.delta;

  (*topCand) = thisCandidate;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* GetHoughPatchTopCandidate() */



/** Print Hough candidates */
void PrintSemiCohCandidates(LALStatus *status,
			    SemiCohCandidateList *in,
			    FILE *fp,
			    LIGOTimeGPS refTime)
{
  INT4 k;
  PulsarSpins  fkdotIn, fkdotOut;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  XLAL_INIT_MEM ( fkdotIn );
  XLAL_INIT_MEM ( fkdotOut );

  for (k=0; k < in->nCandidates; k++) {
    /*     fprintf(fp, "%e %e %e %g %g %g %g %e %e\n", in->list[k].significance, in->list[k].freq,  */
    /* 	    in->list[k].dFreq, in->list[k].alpha, in->list[k].dAlpha, in->list[k].delta, in->list[k].dDelta, */
    /* 	    in->list[k].fdot, in->list[k].dFdot); */

    fkdotIn[0] = in->list[k].freq;
    fkdotIn[1] = in->list[k].fdot;

    TRY( LALExtrapolatePulsarSpins ( status->statusPtr, fkdotOut, refTime, fkdotIn, in->refTime), status);

    fprintf(fp, "%f %f %f %e %f\n", fkdotOut[0], in->list[k].alpha, in->list[k].delta,
	    fkdotOut[1], in->list[k].significance);
  }


  DETATCHSTATUSPTR (status);
  RETURN(status);
}


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

  XLAL_INIT_MEM(fkdot);

  fprintf(fp, "## Fstat values from stack %d (reftime -- %d %d)\n", stackIndex, refTime.gpsSeconds, refTime.gpsNanoSeconds);

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

      fprintf(fp, "%.13g %.7g %.7g %.5g %.6g\n", fkdot[0], alpha, delta, fkdot[1], 2*in->data->data[k]);
    }

  fprintf(fp, "\n");

  DETATCHSTATUSPTR (status);
  RETURN(status);

}

/** Print hough histogram to a file */
void PrintHoughHistogram( LALStatus *status,
			  UINT8Vector *hist,
			  CHAR *fnameOut)
{

  FILE  *fp=NULL;   /* Output file */
  char filename[256];
  UINT4  i ;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  strcpy(  filename, fnameOut);
  strcat(  filename, "histo");
  if ( !(fp=fopen(filename,"w")))
    {
      fprintf(stderr, "Unable to open file '%s' for writing\n", filename);
      exit(1);
    }

  for (i=0; i < hist->length; i++)
    fprintf(fp,"%d  %" LAL_UINT8_FORMAT "\n", i, hist->data[i]);

  fclose( fp );

  DETATCHSTATUSPTR (status);
  RETURN(status);

}


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


/**
 * Read checkpointing file
 * This does not (yet) check any consistency of
 * the existing results file
 */
void GetChkPointIndex( LALStatus *status,
		       INT4 *loopindex,
		       const CHAR *fnameChkPoint)
{

  FILE  *fp=NULL;
  UINT4 tmpIndex;
  CHAR lastnewline='\0';

  INITSTATUS(status);
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

/**
 * Calculate average velocity and position of detector network during each
 * stack
 */
void GetStackVelPos( LALStatus *status,
		     REAL8VectorSequence **velStack,
		     REAL8VectorSequence **posStack,
		     FstatInputVector* Fstat_in_vec)
{
  UINT4 k, j, m, nStacks;
  INT4 counter, numifo;
  CreateVectorSequenceIn createPar;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);


  ASSERT ( Fstat_in_vec != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( Fstat_in_vec->length > 0, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  ASSERT ( Fstat_in_vec->data != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( velStack != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( posStack != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

  nStacks = Fstat_in_vec->length;

  /* create velocity and position vectors */
  createPar.length = nStacks; /* number of vectors */
  createPar.vectorLength = 3; /* length of each vector */
  TRY( LALDCreateVectorSequence ( status->statusPtr,  velStack, &createPar), status);
  TRY( LALDCreateVectorSequence ( status->statusPtr,  posStack, &createPar), status);


  /* calculate detector velocity and position for each stack*/
  /* which detector?  -- think carefully about this! */
  /* Since only the sfts are associated with a detector, the cleanest
     solution (unless something better comes up) seems to be to
     average the velocities and positions of the ifo within the sft
     time intervals in each stack.  Thus, for the purposes of the
     velocity and position calculation, we can view the sfts as coming
     from the a single hypothetical ifo which is moving in a strange
     way */

  for (k = 0; k < nStacks; k++)
    {
      const MultiDetectorStateSeries *multiDetStates = XLALGetFstatInputDetectorStates(Fstat_in_vec->data[k]);
      counter=0;
      numifo = multiDetStates->length;

      /* initialize velocities and positions */
      velStack[0]->data[3*k] = 0;
      velStack[0]->data[3*k+1] = 0;
      velStack[0]->data[3*k+2] = 0;

      posStack[0]->data[3*k] = 0;
      posStack[0]->data[3*k+1] = 0;
      posStack[0]->data[3*k+2] = 0;

      for ( j = 0; (INT4)j < numifo; j++)
	{
	  INT4 numsft = multiDetStates->data[j]->length;
	  for ( m = 0; (INT4)m < numsft; m++)
	    {
	      /* sum velocity components */
	      velStack[0]->data[3*k] += multiDetStates->data[j]->data[m].vDetector[0];
	      velStack[0]->data[3*k+1] += multiDetStates->data[j]->data[m].vDetector[1];
	      velStack[0]->data[3*k+2] += multiDetStates->data[j]->data[m].vDetector[2];
	      /* sum position components */
	      posStack[0]->data[3*k] += multiDetStates->data[j]->data[m].rDetector[0];
	      posStack[0]->data[3*k+1] += multiDetStates->data[j]->data[m].rDetector[1];
	      posStack[0]->data[3*k+2] += multiDetStates->data[j]->data[m].rDetector[2];

	      counter++;

	    } /* loop over sfts for this ifo */

	} /* loop over ifos in stack */

      /* divide by 'counter' to get average */
      velStack[0]->data[3*k] /= counter;
      velStack[0]->data[3*k+1] /= counter;
      velStack[0]->data[3*k+2] /= counter;

      posStack[0]->data[3*k] /= counter;
      posStack[0]->data[3*k+1] /= counter;
      posStack[0]->data[3*k+2] /= counter;

    } /* loop over stacks -- end velocity and position calculation */

  DETATCHSTATUSPTR (status);
  RETURN(status);

}


/** Calculate noise weight for each stack*/
void ComputeStackNoiseWeights( LALStatus *status,
			       REAL8Vector **out,
			       FstatInputVector* Fstat_in_vec )
{

  INT4 nStacks, k, j, i, numifo, numsft;
  REAL8Vector *weightsVec=NULL;


  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( Fstat_in_vec != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( Fstat_in_vec->length > 0, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  ASSERT ( Fstat_in_vec->data != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

  ASSERT ( out != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( *out == NULL, status, HIERARCHICALSEARCH_ENONULL, HIERARCHICALSEARCH_MSGENONULL );

  nStacks = Fstat_in_vec->length;

  TRY( LALDCreateVector( status->statusPtr, &weightsVec, nStacks), status);


  for (k=0; k<nStacks; k++) {

    const MultiNoiseWeights *multNoiseWts = XLALGetFstatInputNoiseWeights(Fstat_in_vec->data[k]);
    ASSERT ( multNoiseWts != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

    numifo = multNoiseWts->length;
    ASSERT ( numifo > 0, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );

    /* initialize */
    weightsVec->data[k] = 0;

    for ( j = 0; j < numifo; j++) {

      numsft = multNoiseWts->data[j]->length;

      for ( i = 0; i < numsft; i++) {

	weightsVec->data[k] += multNoiseWts->data[j]->data[i];

      } /* loop over sfts */

    }/* loop over ifos in stack */

  } /* loop over stacks*/


  LAL_CALL( LALHOUGHNormalizeWeights( status->statusPtr, weightsVec), status);

  *out = weightsVec;

  DETATCHSTATUSPTR (status);
  RETURN(status);

}




/** Calculate noise and AM weight for each stack for a given sky position*/
void ComputeStackNoiseAndAMWeights( LALStatus *status,
				    REAL8Vector *out,
				    FstatInputVector* Fstat_in_vec,
				    SkyPosition skypos)
{

  UINT4 nStacks, iStack, iIFO, iSFT, numifo, numsft;
  REAL8 a, b, n;
  MultiAMCoeffs *multiAMcoef = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( Fstat_in_vec != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( Fstat_in_vec->length > 0, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  ASSERT ( Fstat_in_vec->data != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

  nStacks = Fstat_in_vec->length;

  ASSERT ( out != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( out->length == nStacks, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  ASSERT ( out->data != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );


  for (iStack=0; iStack<nStacks; iStack++) {

    const MultiNoiseWeights *multNoiseWts = XLALGetFstatInputNoiseWeights(Fstat_in_vec->data[iStack]);
    ASSERT ( multNoiseWts != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

    const MultiDetectorStateSeries *multDetStates = XLALGetFstatInputDetectorStates(Fstat_in_vec->data[iStack]);
    ASSERT ( multDetStates != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

    numifo = multNoiseWts->length;
    ASSERT ( numifo > 0, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
    ASSERT ( multDetStates->length == numifo, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );

    /* initialize */
    out->data[iStack] = 0;

    multiAMcoef = NULL;
    TRY ( LALGetMultiAMCoeffs ( status->statusPtr, &multiAMcoef, multDetStates, skypos), status);


    for ( iIFO = 0; iIFO < numifo; iIFO++) {

      numsft = multNoiseWts->data[iIFO]->length;
      ASSERT ( multDetStates->data[iIFO]->length == numsft, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );

      for ( iSFT = 0; iSFT < numsft; iSFT++) {

	a = multiAMcoef->data[iIFO]->a->data[iSFT];
	b = multiAMcoef->data[iIFO]->b->data[iSFT];
	n =  multNoiseWts->data[iIFO]->data[iSFT];

	out->data[iStack] += (a*a + b*b)*n;

      } /* loop over sfts */

    }/* loop over ifos in stack */

    XLALDestroyMultiAMCoeffs ( multiAMcoef );

  } /* loop over stacks*/


  TRY( LALHOUGHNormalizeWeights( status->statusPtr, out), status);

  DETATCHSTATUSPTR (status);
  RETURN(status);

}


/** Get SemiCoh candidates toplist */
void GetSemiCohToplist(LALStatus            *status,
		       toplist_t            *list,
		       SemiCohCandidateList *in,
		       REAL8                meanN,
		       REAL8                sigmaN)
{

  INT4 k;
  HoughFStatOutputEntry line;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( list != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( in != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( in->length >= in->nCandidates, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );

  /* go through candidates and insert into toplist if necessary */
  for ( k = 0; k < in->nCandidates; k++) {

    line.Freq = in->list[k].freq;
    line.Alpha = in->list[k].alpha;
    line.Delta = in->list[k].delta;
    line.f1dot = in->list[k].fdot;
    line.HoughFStat = (in->list[k].significance - meanN)/sigmaN; 
    /* for debugging */
    /* line.HoughFStat = in->list[k].significance; */
    /* if (line.HoughFStat > 121) */
    /*   fprintf(stdout, "number count exceeded"); */

    line.AlphaBest = in->list[k].alphaBest;

    if ( line.AlphaBest < 0 ) {
      line.AlphaBest += LAL_TWOPI;
    }
    if ( line.AlphaBest > LAL_TWOPI ) {
      line.AlphaBest -= LAL_TWOPI;
    }

    line.DeltaBest = in->list[k].deltaBest;
    line.MeanSig = (in->list[k].meanSig - meanN) / sigmaN;
    line.VarianceSig = in->list[k].varianceSig / (sigmaN * sigmaN);

    /* initialize LV postprocessing entries to zero.
     * This will be filled later, if user requested it.
     */
    line.sumTwoF = -1;    	/* sum of 2F-values for LV postprocessing */
    line.sumTwoFX = NULL; 	/* sum of 2F-values per detector for LV postprocessing */

    INSERT_INTO_HOUGHFSTAT_TOPLIST( list, line);

  }

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* GetSemiCohToplist() */



/** Optimized calculation of Fstat overhead */
void ComputeNumExtraBins(LALStatus            *status,
			 SemiCoherentParams   *par,
			 REAL8                fdot,
			 REAL8                f0,
			 REAL8                deltaF)
{

  UINT4 i, j, nStacks, extraBins2, tmpExtraBins2;
  HOUGHptfLUT  lut;
  HOUGHDemodPar parDem;
  HOUGHParamPLUT  parLut;
  HOUGHSizePar   parSize;
  HOUGHResolutionPar parRes;
  REAL8 tMidStack, tStart, tEnd;
  REAL8Vector *timeDiffV=NULL;
  HOUGHPatchGrid  patch;
  REAL8 pixelFactor, alpha, delta, patchSizeX, patchSizeY, refTime;
  LIGOTimeGPS refTimeGPS;

  REAL8VectorSequence  *vel;
  REAL8VectorSequence  *pos;
  LIGOTimeGPSVector    *tsMid;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  vel = par->vel;
  pos = par->pos;
  pixelFactor = par->pixelFactor;
  fdot = par->fdot;
  alpha = par->alpha;
  delta = par->delta;
  patchSizeX = par->patchSizeX;
  patchSizeY = par->patchSizeY;
  tsMid = par->tsMid;

  refTimeGPS = par->refTime;
  refTime = XLALGPSGetREAL8(&refTimeGPS);

  nStacks = tsMid->length;;

  tStart = XLALGPSGetREAL8(tsMid->data);
  tEnd = XLALGPSGetREAL8(tsMid->data + nStacks - 1);
  refTime = 0.5 * (tStart + tEnd);

  /* the skygrid resolution params */
  parRes.deltaF = deltaF;
  parRes.patchSkySizeX  = patchSizeX;
  parRes.patchSkySizeY  = patchSizeY;
  parRes.pixelFactor = pixelFactor;
  parRes.pixErr = PIXERR;
  parRes.linErr = LINERR;
  parRes.vTotC = VTOT;
  parRes.f0Bin =  (INT8)(f0/deltaF + 0.5);

  TRY( LALHOUGHComputeSizePar( status->statusPtr, &parSize, &parRes ),  status );

  patch.xSide = parSize.xSide;
  patch.ySide = parSize.ySide;
  patch.xCoor = NULL;
  patch.yCoor = NULL;
  patch.xCoor = (REAL8 *)LALCalloc(1,patch.xSide*sizeof(REAL8));
  patch.yCoor = (REAL8 *)LALCalloc(1,patch.ySide*sizeof(REAL8));
  TRY( LALHOUGHFillPatchGrid( status->statusPtr, &patch, &parSize ), status );

  /* the demodulation params */
  parDem.deltaF = deltaF;
  parDem.skyPatch.alpha = alpha;
  parDem.skyPatch.delta = delta;
  parDem.spin.length = 1;
  parDem.spin.data = NULL;
  parDem.spin.data = (REAL8 *)LALCalloc(1, sizeof(REAL8));
  parDem.spin.data[0] = fdot;

  TRY( LALDCreateVector( status->statusPtr, &timeDiffV, nStacks), status);

  for (j=0; j<nStacks; j++) {
    tMidStack = XLALGPSGetREAL8(tsMid->data + j);
    timeDiffV->data[j] = tMidStack - refTime;
  }


  lut.maxNBins = parSize.maxNBins;
  lut.maxNBorders = parSize.maxNBorders;
  lut.border = (HOUGHBorder *)LALCalloc(1, parSize.maxNBorders*sizeof(HOUGHBorder));
  lut.bin = (HOUGHBin2Border *)LALCalloc(1, parSize.maxNBins*sizeof(HOUGHBin2Border));
  for (i = 0; i < parSize.maxNBorders; i++){
    lut.border[i].ySide = parSize.ySide;
    lut.border[i].xPixel = (COORType *)LALCalloc(1,parSize.ySide*sizeof(COORType));
  }

  /* loop over stacks and create LUTs */
  extraBins2 = 0;
  for (j=0; j < (UINT4)nStacks; j++){
    parDem.veloC.x = vel->data[3*j];
    parDem.veloC.y = vel->data[3*j + 1];
    parDem.veloC.z = vel->data[3*j + 2];
    parDem.positC.x = pos->data[3*j];
    parDem.positC.y = pos->data[3*j + 1];
    parDem.positC.z = pos->data[3*j + 2];
    parDem.timeDiff = timeDiffV->data[j];

    /* calculate parameters needed for buiding the LUT */
    TRY( LALHOUGHCalcParamPLUT( status->statusPtr, &parLut, &parSize, &parDem),status );
    /* build the LUT */
    TRY( LALHOUGHConstructPLUT( status->statusPtr, &lut, &patch, &parLut ),
	 status );

    tmpExtraBins2 = HSMAX( abs(lut.iniBin), abs(lut.nBin + lut.iniBin) );

    if ( tmpExtraBins2 > extraBins2)
      extraBins2 = tmpExtraBins2;
  } /* LUTs are created */


  par->extraBinsFstat = extraBins2 + 2;

  /* now free memory and exit */
  TRY( LALDDestroyVector( status->statusPtr, &timeDiffV), status);

  LALFree(patch.xCoor);
  LALFree(patch.yCoor);

  LALFree(parDem.spin.data);

  for (i = 0; i < lut.maxNBorders; i++){
    LALFree( lut.border[i].xPixel);
  }
  LALFree( lut.border);
  LALFree( lut.bin);


  DETATCHSTATUSPTR (status);
  RETURN(status);


  /*   REAL8 skypos[3], xi[3], xiProj[3], xhat[3], yhat[3]; */
  /*   REAL8 xiDotn, xiProjX, xiProjY, patchSizeNorm, xiProjNorm; */


  /*   nStacks = vel->vectorLength; */
  /*   extraBins = 0; */

  /*   skypos[0] = cos(delta) * cos(alpha); */
  /*   skypos[1] = cos(delta) * sin(alpha); */
  /*   skypos[2] = sin(delta); */

  /*   xhat[0] = -sin(alpha); */
  /*   xhat[1] = cos(alpha); */
  /*   xhat[3] = 0; */

  /*   yhat[0] = sin(delta)*cos(alpha); */
  /*   yhat[1] = sin(delta)*sin(alpha); */
  /*   yhat[2] = -cos(delta); */

  /*   /\* loop over stacks *\/ */
  /*   for ( k = 0; k < nStacks; k++) { */

  /*     UINT4 minBinsPar, minBinsPerp; */


  /*     xi[0] = vel->data[3*k]; */
  /*     xi[1] = vel->data[3*k+1]; */
  /*     xi[2] = vel->data[3*k+2]; */

  /*     xiDotn = xi[0]*skypos[0] + xi[1]*skypos[1] + xi[2]*skypos[2]; */

  /*     xiProj[0] = xi[0] - xiDotn * skypos[0]; */
  /*     xiProj[1] = xi[1] - xiDotn * skypos[1]; */
  /*     xiProj[2] = xi[2] - xiDotn * skypos[2]; */

  /*     xiProjX = xiProj[0]*xhat[0] + xiProj[1]*xhat[1] + xiProj[2]*xhat[2]; */
  /*     xiProjY = xiProj[0]*yhat[0] + xiProj[1]*yhat[1] + xiProj[2]*yhat[2]; */

  /*     xiProjNorm = sqrt(xiProj[0]*xiProj[0] + xiProj[0]*xiProj[0] + xiProj[0]*xiProj[0]); */
  /*     patchSizeNorm = patchSizeX*patchSizeX + patchSizeY*patchSizeY; */

  /*     minBinsPar = (UINT4)(f0 * fabs(xiDotn) * patchSizeNorm /deltaF + 0.5); */

  /*     minBinsPerp = (UINT4)(f0 * xiProjNorm * patchSizeNorm/ deltaF + 0.5); */

  /*     extraBinsInThisStack = minBinsPar + minBinsPerp; */

  /*     if (extraBins < extraBinsInThisStack) */
  /*       extraBins = extraBinsInThisStack; */

  /*   } /\* loop over stacks *\/ */



}   /*ComputeNumExtraBins()*/





void GetXiInSingleStack (LALStatus         *status,
			 HOUGHSizePar      *size,
			 HOUGHDemodPar     *par)  /* demodulation parameters */
{

  /* --------------------------------------------- */

  REAL8   f0;  /* frequency corresponding to f0Bin */
  INT8    f0Bin;
  REAL8   deltaF;  /*  df=1/TCOH  */
  REAL8   vFactor, xFactor;
  REAL8   xiX, xiY, xiZ;
  REAL8   modXi,invModXi;
  REAL8UnitPolarCoor   xiInit, xiFinal;
  UINT4   spinOrder, i;
  REAL8   *spinF;
  REAL8   timeDiff;    /*  T(t)-T(t0) */
  REAL8   timeDiffProd;
  /* --------------------------------------------- */

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (par, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL);
  ASSERT (size, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL);

  /*   Make sure f0Bin  is not zero: */
  f0Bin = size->f0Bin;
  ASSERT (f0Bin, status, HIERARCHICALSEARCH_EBAD, HIERARCHICALSEARCH_MSGEBAD);


  deltaF = size->deltaF;

  f0 = f0Bin * deltaF;


  /* *********** xi calculation *****************  */

  vFactor = f0;
  xFactor = 0.0;

  spinOrder = par->spin.length;

  if(spinOrder){
    ASSERT (par->spin.data , status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL);
    timeDiff = par->timeDiff;
    timeDiffProd = 1.0;
    spinF = par->spin.data;

    for (i=0; i<spinOrder; ++i ){
      xFactor += spinF[i] * timeDiffProd * (i+1.0);
      timeDiffProd *= timeDiff;
      vFactor += spinF[i] * timeDiffProd;
    }
  }

  xiX = vFactor * (par->veloC.x) + xFactor * (par->positC.x);
  xiY = vFactor * (par->veloC.y) + xFactor * (par->positC.y);
  xiZ = vFactor * (par->veloC.z) + xFactor * (par->positC.z);

  /* -------------------------------------------   */
  /* ***** convert xi into Polar coordinates ***** */

  modXi = sqrt(xiX*xiX + xiY*xiY + xiZ*xiZ);
  /* for testing we used:  modXi = F0* 1.06e-4; */
  invModXi = 1./modXi;

  xiInit.delta = asin( xiZ*invModXi);
  /* the arc sine is in the interval [-pi/2,pi/2] */

  if( xiX || xiY ){
    xiInit.alpha = atan2(xiY, xiX);
  }else{
    xiInit.alpha = 0.0;
  }

  /*   if( (xiX == 0.0 ) && (xiY == 0.0 ) ){ */
  /*     xiInit.alpha = 0.0; */
  /*   }else{  xiInit.alpha = atan2(xiY, xiX);  } */

  /* -------------------------------------------   */
  /* **** Rotate Patch, so that its center becomes */
  /* **** the south pole {x,y,z} = {0,0,-1}  ***** */
  /*  Calculate xi in the new coordinate system.   */

  TRY(LALRotatePolarU(status->statusPtr,
		      &xiFinal ,&xiInit, &(*par).skyPatch), status);


  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}


/**
 * Load hough-candidate list from given file
 */
HoughCandidateList *
XLALLoadHoughCandidateList ( const char *fname,	/**< input candidate-list file 'Freq Alpha Delta f1dot sig' */
                             REAL8 FreqShift	/**< apply this shift to input frequencies to correct HS offset-bug */
                             )
{
  if ( ! fname )
    XLAL_ERROR_NULL ( XLAL_EINVAL );

  LALParsedDataFile *data = NULL;
  if ( XLALParseDataFile (&data, fname) != XLAL_SUCCESS )
    XLAL_ERROR_NULL ( XLAL_EFUNC );

  UINT4 nCands = data->lines->nTokens;

  HoughCandidateList *out;
  if ( ( out = XLALCreateHoughCandidateList ( nCands )) == NULL )
    XLAL_ERROR_NULL ( XLAL_EFUNC );

  HoughCandidate candMax = { 0,   -1, -2, -1, 0 };	/* initialize with impossibly small values */
  HoughCandidate candMin = { 1e9, 10, 10,  1, 0 }; /* initialize with impossibly large values */

  /* parse the data-file as a list of lines */
  UINT4 i;
  for (i=0; i < nCands; i++)
    {
      HoughCandidate *cand = &out->data[i];
      if ( 5 != sscanf( data->lines->tokens[i], "%lg %lg %lg %lg %lg", &cand->Freq, &cand->Alpha, &cand->Delta, &cand->f1dot, &cand->sig ))
	{
	  XLALPrintError ( "%s: could not parse 5 numbers from line %d in candidate-file '%s':\n", __func__, i, fname);
          XLALPrintError ("'%s'\n", data->lines->tokens[i] );
          XLALDestroyHoughCandidateList ( out );
          XLAL_ERROR_NULL (   XLAL_EDATA );
	}
      /* apply frequency correction */
      cand->Freq += FreqShift;

      candMax.Freq  = HSMAX ( candMax.Freq,  cand->Freq );
      candMax.Alpha = HSMAX ( candMax.Alpha, cand->Alpha );
      candMax.Delta = HSMAX ( candMax.Delta, cand->Delta );
      candMax.f1dot = HSMAX ( candMax.f1dot, cand->f1dot );

      candMin.Freq  = HSMIN ( candMin.Freq,  cand->Freq );
      candMin.Alpha = HSMIN ( candMin.Alpha, cand->Alpha );
      candMin.Delta = HSMIN ( candMin.Delta, cand->Delta );
      candMin.f1dot = HSMIN ( candMin.f1dot, cand->f1dot );

    } /* for i < nLines */

  XLALDestroyParsedDataFile ( data );

  /* set frequency/f1dot search-space boundaries */
  out->FreqMin  = candMin.Freq;
  out->FreqBand = candMax.Freq - candMin.Freq;

  out->f1dotMin  = candMin.f1dot;
  out->f1dotBand = candMax.f1dot - candMin.f1dot;
  /* set sky (Alpha,Delta) 'box'-boundaries */
  out->AlphaMin  = candMin.Alpha;
  out->AlphaBand = candMax.Alpha - candMin.Alpha;

  out->DeltaMin  = candMin.Delta;
  out->DeltaBand = candMax.Delta - candMin.Delta;

  return out;

} /* XLALLoadHoughCandidateList */



HoughCandidateList *
XLALCreateHoughCandidateList ( UINT4 length )
{
  HoughCandidateList *out;
  if ( (out = XLALCalloc ( 1, sizeof(HoughCandidateList) )) == NULL )
    XLAL_ERROR_NULL ( XLAL_ENOMEM );

  if ( (out->data = XLALCalloc ( length, sizeof(*out->data) )) == NULL )
    XLAL_ERROR_NULL ( XLAL_ENOMEM );

  out->length = length;

  return out;

} /* XLALCreateHoughCandidateList() */


void
XLALDestroyHoughCandidateList ( HoughCandidateList * list )
{
  if ( !list )
    return;

  if ( list->data )
    XLALFree ( list->data );

  XLALFree ( list );
  return;
} /* XLALDestroyHoughCandidateList () */
