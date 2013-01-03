/*
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
   and combining them semi-coherently using the Hough transform, and following
   up candidates using a longer coherent integration.

   \par  Description

   This code implements a hierarchical strategy to look for unknown gravitational
   wave pulsars. It scans through the parameter space using a less sensitive but
   computationally inexpensive search and follows up the candidates using
   more sensitive methods.

   \par Algorithm

   Currently the code does a single stage hierarchical search using the Hough
   algorithm and follows up the candidates using a full coherent integration.

   - The user specifies a directory containing SFTs, and the number of
     \e stacks that this must be broken up into.  Stacks are chosen by
     breaking up the total time spanned by the sfts into equal
     portions, and then choosing the sfts which lie in each stack
     portion.  The SFTs can be from multiple IFOs.

   - The user specifies a region in parameter space to search over.
     The code sets up a grid (the "coarse" grid) in this region and
     calculates the F-statistic for each stack at each point of the
     coarse grid, for the entire frequency band.  Alternatively, the
     user can specify a sky-grid file.

   - The different Fstat vactors for each stack are combined using a
     semi-coherent method such as the Hough transform or stack slide
     algorithms.

   - For Hough, a threshold is set on the F-statistic to convert the
     F-statistic vector into a vector of 0s and 1s known as a \e
     peakgram -- there is one peakgram for each stack and for each
     grid point.  These peakgrams are combined using the Hough algorithm.

   - For stack-slide, we add the different Fstatistic values to get
     the summed Fstatistic power.

   - The semi-coherent part of the search constructs a grid (the
     "fine" grid) in a small patch around every coarse grid point, and
     combines the different stacks following the \e master equation
     \f[ f(t) - F_0(t) = \xi(t).(\hat{n} - \hat{n}_0) \f]
     where
     \f[ F_0 = f_0 + \sum \Delta f_k \frac{(\Delta t)^k}{k!}  \f]
     Here \f$ \hat{n}_0 \f$ is the sky-point at which the F-statistic is
     calculated and \f$ \Delta f_k \f$ is the \e residual spindown
     parameter.  For details see Phys.Rev.D 70, 082001 (2004).  The
     size of the patch depends on the validity of the above master
     equation.

   - The output of the Hough search is a \e number \e count at each point
     of the grid.  A threshold is set on the number count, leading to
     candidates in parameter space.  For stack slide, instead of the
     number count, we get the summed Fstatistic power.  Alternatively,
     the user can specify exactly how many candidates should be
     followed up for each coarse grid point.

   - These candidates are followed up using a second set of SFTs (also
     specified by the user).  The follow up consists of a full
     coherent integration, i.e. the F-statistic is calculated for the
     whole set of SFTs without breaking them up into stacks.  The user
     can choose to output the N highest significance candidates.


   \par Questions/To-do

   - Should we over-resolve the Fstat calculation to reduce loss in signal power?  We would
     still calculate the peakgrams at the 1/T resolution, but the peak selection would
     take into account Fstat values over several over-resolved bins.

   - What is the best grid for calculating the F-statistic?  At first glance, the
     Hough patch and the metric F-statistic patch do not seem to be compatible.  If we
     use the Hough patches to break up the sky, it could be far from optimal.  What is
     the exact connection between the two grids?

   - Implement multiple semi-coherent stages

   - Get timings and optimize the pipeline parameters

   - Checkpointing for running on Einstein\@Home

   - Incorporate stack slide as an alternative to Hough in the semi-coherent stages

   - ....

 */

#include "OptimizedCFS/ComputeFstatREAL4.h"
#include "HierarchicalSearch.h"


/* This need to go here after HierarchicalSearch.h is included;
   StackSlideFstat.h needs SemiCohCandidateList defined. */
#include "StackSlideFstat.h"

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
#ifndef COMPUTEFSTATHOUGHMAP
#define COMPUTEFSTATHOUGHMAP ComputeFstatHoughMap
#endif
#ifndef COMPUTEFSTATFREQBAND
#define COMPUTEFSTATFREQBAND ComputeFStatFreqBand
#endif
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

#if USE_OPENCL_KERNEL || defined(USE_CUDA)
extern int gpu_device_id;
#else
int gpu_device_id;
#endif

extern int lalDebugLevel;

BOOLEAN uvar_printMaps = FALSE; /**< global variable for printing Hough maps */
BOOLEAN uvar_printGrid = FALSE; /**< global variable for printing Hough grid */
BOOLEAN uvar_printStats = FALSE;/**< global variable for calculating Hough map stats */
BOOLEAN uvar_dumpLUT = FALSE;  	/**< global variable for printing Hough look-up-tables for debugging */
BOOLEAN uvar_validateLUT = FALSE;

#define HSMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define HSMIN(x,y) ( (x) < (y) ? (x) : (y) )

#define INIT_MEM(x) memset(&(x), 0, sizeof((x)))

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
  LIGOTimeGPS maxEndTimeGPS;       /**< all sft data must be before this time */
  UINT4 blocksRngMed;              /**< blocksize for running median noise floor estimation */
  UINT4 Dterms;                    /**< size of Dirichlet kernel for Fstat calculation */
  REAL8 dopplerMax;                /**< extra sft wings for doppler motion */
} UsefulStageVariables;


static int smallerHough(const void *a,const void *b) {
  SemiCohCandidate a1, b1;
  a1 = *((const SemiCohCandidate *)a);
  b1 = *((const SemiCohCandidate *)b);

  if( a1.significance < b1.significance )
    return(1);
  else if( a1.significance > b1.significance)
    return(-1);
  else
    return(0);
}

/* functions for printing various stuff */
void ComputeStackNoiseWeights( LALStatus *status, REAL8Vector **out, MultiNoiseWeightsSequence *in );

void ComputeStackNoiseAndAMWeights( LALStatus *status, REAL8Vector *out, MultiNoiseWeightsSequence *inNoise,
				    MultiDetectorStateSeriesSequence *inDetStates, SkyPosition skypos);

void GetStackVelPos( LALStatus *status, REAL8VectorSequence **velStack, REAL8VectorSequence **posStack,
		     MultiDetectorStateSeriesSequence *stackMultiDetStates);


void SetUpSFTs( LALStatus *status, MultiSFTVectorSequence *stackMultiSFT, MultiNoiseWeightsSequence *stackMultiNoiseWeights,
		MultiDetectorStateSeriesSequence *stackMultiDetStates, UsefulStageVariables *in);

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

/* default values for input variables */
#define EARTHEPHEMERIS 		"earth05-09.dat"
#define SUNEPHEMERIS 		"sun05-09.dat"

#define BLOCKSRNGMED 		101 	/**< Default running median window size */
#define FSTART 			310.0	/**< Default Start search frequency */

#define FBAND 			0.01	/**< Default search band */
#define FDOT 			0.0	/**< Default value of first spindown */
#define DFDOT 			0.0	/**< Default range of first spindown parameter */
#define SKYREGION 		"allsky" /**< default sky region to search over -- just a single point*/
#define NFDOT  			10    	/**< Default size of hough cylinder of look up tables */
#define DTERMS 			8     	/**< Default number of dirichlet kernel terms for calculating Fstat */
#define MISMATCH 		0.2 	/**< Default for metric grid maximal mismatch value */
#define DALPHA 			0.001 	/**< Default resolution for isotropic or flat grids */
#define DDELTA 			0.001 	/**< Default resolution for isotropic or flat grids */
#define FSTATTHRESHOLD 		2.6	/**< Default threshold on Fstatistic for peak selection */
#define NCAND1 			5 	/**< Default number of candidates to be followed up from first stage */
#define FNAMEOUT 		"./out/HS.dat"  /**< Default output file basename */
#define PIXELFACTOR 		2.0
#ifndef LAL_INT4_MAX
#define LAL_INT4_MAX 		2147483647
#endif

/* a global pointer to MAIN()s head of the LALStatus structure,
   made global so a signal handler can read it */
LALStatus *global_status;

int alloc_len;

#ifdef OUTPUT_TIMING
time_t clock0;
UINT4 nSFTs;
UINT4 nStacks;
UINT4 nSkyRefine;
#endif


int MAIN( int argc, char *argv[]) {

  LALStatus status = blank_status;

  /* temp loop variables: generally k loops over stacks and j over SFTs in a stack*/
  INT4 j;
  UINT4 k;
  UINT4 skyGridCounter;

  /* in general any variable ending with 1 is for the
     first stage, 2 for the second and so on */

  /* ephemeris */
  EphemerisData *edat = NULL;

  /* timestamp vectors */
  LIGOTimeGPSVector *midTstack=NULL;
  LIGOTimeGPSVector *startTstack=NULL;


  LIGOTimeGPS refTimeGPS = empty_LIGOTimeGPS;
  LIGOTimeGPS tMidGPS = empty_LIGOTimeGPS;
  REAL8 tObs;

  /* velocities and positions at midTstack */
  REAL8VectorSequence *velStack=NULL;
  REAL8VectorSequence *posStack=NULL;

  /* weights for each stack*/
  REAL8Vector *weightsV=NULL;
  REAL8Vector *weightsNoise=NULL;

  /* duration of each stack */
  REAL8 tStack;

  /* sft related stuff */
  static MultiSFTVectorSequence stackMultiSFT;
  static MultiNoiseWeightsSequence stackMultiNoiseWeights;
  static MultiDetectorStateSeriesSequence stackMultiDetStates;
  static LIGOTimeGPS minStartTimeGPS, maxEndTimeGPS;


  /* some useful variables for each stage */
  UsefulStageVariables usefulParams;

  /* number of stacks -- not necessarily same as uvar_nStacks! */
  UINT4 nStacks;
  REAL8 dFreqStack; /* frequency resolution of Fstat calculation */
  REAL8 df1dot, df1dotRes;  /* coarse grid resolution in spindown */
  UINT4 nf1dot, nf1dotRes; /* coarse and fine grid number of spindown values */

  /* LALdemod related stuff */
  static REAL4FrequencySeriesVector fstatVector; /* Fstatistic vectors for each stack */
  UINT4 binsFstat1, binsFstatSearch;
  static ComputeFParams CFparams;

  /* hough variables */
  static HOUGHPeakGramVector pgV;
  static SemiCoherentParams semiCohPar;
  static SemiCohCandidateList semiCohCandList;
  REAL8 alphaPeak, sumWeightSquare, meanN=0, sigmaN=0;

  /* fstat candidate structure */
  toplist_t *semiCohToplist=NULL;

  /* template and grid variables */
  static DopplerSkyScanInit scanInit;   /* init-structure for DopperScanner */
  DopplerSkyScanState thisScan = empty_DopplerSkyScanState; /* current state of the Doppler-scan */
  static PulsarDopplerParams dopplerpos;	       /* current search-parameters */
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

  BOOLEAN uvar_printCand1 = FALSE; 	/* if 1st stage candidates are to be printed */
  BOOLEAN uvar_printFstat1 = FALSE;
  BOOLEAN uvar_useToplist1 = FALSE;
  BOOLEAN uvar_useWeights  = FALSE;
  BOOLEAN uvar_semiCohToplist = FALSE; /* if overall first stage candidates are to be output */

  REAL8 uvar_dAlpha = DALPHA; 	/* resolution for flat or isotropic grids -- coarse grid*/
  REAL8 uvar_dDelta = DDELTA;
  REAL8 uvar_f1dot = FDOT; 	/* first spindown value */
  REAL8 uvar_f1dotBand = DFDOT; /* range of first spindown parameter */
  REAL8 uvar_Freq = FSTART;
  REAL8 uvar_FreqBand = FBAND;

  REAL8 uvar_dFreq = 0;
  REAL8 uvar_df1dot = 0; /* coarse grid frequency and spindown resolution */

  REAL8 uvar_peakThrF = FSTATTHRESHOLD; /* threshold of Fstat to select peaks */
  REAL8 uvar_mismatch1 = MISMATCH; /* metric mismatch for first stage coarse grid */

  REAL8 uvar_pixelFactor = PIXELFACTOR;
  REAL8 uvar_df1dotRes = 0;
  REAL8 uvar_threshold1 = 0;
  REAL8 uvar_minStartTime1 = 0;
  REAL8 uvar_maxEndTime1 = LAL_INT4_MAX;
  REAL8 uvar_dopplerMax = 1.05e-4;

  REAL8 uvar_refTime = 0;
  REAL8 uvar_semiCohPatchX = 0;
  REAL8 uvar_semiCohPatchY = 0;
  REAL8 uvar_tStack = 0;

  INT4 uvar_method = -1; 	/* hough = 0, stackslide = 1, -1 = pure fstat*/
  INT4 uvar_nCand1 = NCAND1; /* number of candidates to be followed up from first stage */

  INT4 uvar_blocksRngMed = BLOCKSRNGMED;
  INT4 uvar_nStacksMax = 1;
  INT4 uvar_Dterms = DTERMS;
  INT4 uvar_SSBprecision = SSBPREC_RELATIVISTIC;
  INT4 uvar_nf1dotRes = 1;
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
  BOOLEAN uvar_version = 0;
  INT4 uvar_partitionIndex = 0;

  BOOLEAN uvar_correctFreqs = TRUE;

#ifndef GPUREADY_DEFAULT
#define GPUREADY_DEFAULT 0
#endif
  BOOLEAN uvar_GPUready = GPUREADY_DEFAULT;
  INT4 uvar_gpu_device = -1;
  global_status = &status;

  /* LALDebugLevel must be called before any LALMallocs have been used */
  lalDebugLevel = 0;
  LAL_CALL( LALGetDebugLevel( &status, argc, argv, 'd'), &status);
#ifdef EAH_LALDEBUGLEVEL
  lalDebugLevel = EAH_LALDEBUGLEVEL;
#endif
  uvar_ephemE = LALCalloc( alloc_len = strlen( EARTHEPHEMERIS ) + 1, sizeof(CHAR) );
  XLAL_CHECK ( uvar_ephemE != NULL, XLAL_ENOMEM, "Failed to allocated memory LALCalloc(1, %d)\n", alloc_len );
  strcpy(uvar_ephemE, EARTHEPHEMERIS);

  uvar_ephemS = LALCalloc( alloc_len = strlen(SUNEPHEMERIS) + 1, sizeof(CHAR) );
  XLAL_CHECK ( uvar_ephemS != NULL, XLAL_ENOMEM, "Failed to allocated memory LALCalloc(1, %d)\n", alloc_len );
  strcpy(uvar_ephemS, SUNEPHEMERIS);

  uvar_skyRegion = LALCalloc( alloc_len = strlen(SKYREGION) + 1, sizeof(CHAR) );
  XLAL_CHECK ( uvar_skyRegion != NULL, XLAL_ENOMEM, "Failed to allocated memory LALCalloc(1, %d)\n", alloc_len );
  strcpy(uvar_skyRegion, SKYREGION);

  uvar_fnameout = LALCalloc( alloc_len = strlen(FNAMEOUT) + 1, sizeof(CHAR) );
  XLAL_CHECK ( uvar_fnameout != NULL, XLAL_ENOMEM, "Failed to allocated memory LALCalloc(1, %d)\n", alloc_len );
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
  LAL_CALL( LALRegisterINTUserVar(    &status, "method",       0,  UVAR_OPTIONAL, "0=Hough,1=stackslide,-1=fstat", &uvar_method ), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "semiCohToplist",0, UVAR_OPTIONAL, "Print semicoh toplist?", &uvar_semiCohToplist ), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "useWeights",   0,  UVAR_OPTIONAL, "Weight each stack using noise and AM?", &uvar_useWeights ), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "DataFiles1",   0,  UVAR_REQUIRED, "1st SFT file pattern", &uvar_DataFiles1), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "skyRegion",    0,  UVAR_OPTIONAL, "Sky-region by polygon of form '(ra1,dec1),(ra2,dec2),(ra3,dec3),...' or 'allsky'", &uvar_skyRegion), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "numSkyPartitions",0,UVAR_OPTIONAL, "Number of (equi-)partitions to split skygrid into", &uvar_numSkyPartitions), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "partitionIndex",0, UVAR_OPTIONAL, "Index [0,numSkyPartitions-1] of sky-partition to generate", &uvar_partitionIndex), &status);

  LAL_CALL( LALRegisterREALUserVar(   &status, "Freq",        'f', UVAR_OPTIONAL, "Start search frequency", &uvar_Freq), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dFreq",        0,  UVAR_OPTIONAL, "Frequency resolution (default=1/Tstack)", &uvar_dFreq), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "FreqBand",    'b', UVAR_OPTIONAL, "Search frequency band", &uvar_FreqBand), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "f1dot",        0,  UVAR_OPTIONAL, "Spindown parameter", &uvar_f1dot), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "df1dot",       0,  UVAR_OPTIONAL, "Spindown resolution (default=1/Tstack^2)", &uvar_df1dot), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "f1dotBand",    0,  UVAR_OPTIONAL, "Spindown Range", &uvar_f1dotBand), &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "nf1dotRes",    0,  UVAR_OPTIONAL, "No.of residual fdot values (default=nStacks)", &uvar_nf1dotRes), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "nStacksMax",   0,  UVAR_OPTIONAL, "Maximum No. of 1st stage stacks", &uvar_nStacksMax ),&status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "tStack",       0,  UVAR_REQUIRED, "Duration of 1st stage stacks (sec)", &uvar_tStack ),&status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "mismatch1",    0,  UVAR_OPTIONAL, "1st stage mismatch", &uvar_mismatch1), &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "gridType1",    0,  UVAR_OPTIONAL, "0=flat,1=isotropic,2=metric,3=file", &uvar_gridType1),  &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "metricType1",  0,  UVAR_OPTIONAL, "0=none,1=Ptole-analytic,2=Ptole-numeric,3=exact", &uvar_metricType1), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "skyGridFile",  0,  UVAR_OPTIONAL, "sky-grid file", &uvar_skyGridFile), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dAlpha",       0,  UVAR_OPTIONAL, "Resolution for flat or isotropic coarse grid", &uvar_dAlpha), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dDelta",       0,  UVAR_OPTIONAL, "Resolution for flat or isotropic coarse grid", &uvar_dDelta), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "pixelFactor",  0,  UVAR_OPTIONAL, "Semi coh. sky resolution = 1/v*pixelFactor*f*Tcoh", &uvar_pixelFactor), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "semiCohPatchX",0,  UVAR_OPTIONAL, "Semi coh. sky grid size (default = 1/f*Tcoh*Vepi)", &uvar_semiCohPatchX), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "semiCohPatchY",0,  UVAR_OPTIONAL, "Semi coh. sky grid size (default = 1/f*Tcoh*Vepi)", &uvar_semiCohPatchY), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "fnameout",    'o', UVAR_OPTIONAL, "Output fileneme", &uvar_fnameout), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "peakThrF",     0,  UVAR_OPTIONAL, "Fstat Threshold", &uvar_peakThrF), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "nCand1",       0,  UVAR_OPTIONAL, "No.of 1st stage candidates to be followed up", &uvar_nCand1), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "threshold1",   0,  UVAR_OPTIONAL, "Threshold on significance for 1st stage (if no toplist)", &uvar_threshold1), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printCand1",   0,  UVAR_OPTIONAL, "Print 1st stage candidates", &uvar_printCand1), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "refTime",      0,  UVAR_OPTIONAL, "Ref. time for pulsar pars [Default: mid-time]", &uvar_refTime), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "ephemE",       0,  UVAR_OPTIONAL, "Location of Earth ephemeris file", &uvar_ephemE),  &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "ephemS",       0,  UVAR_OPTIONAL, "Location of Sun ephemeris file", &uvar_ephemS),  &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "minStartTime1",0,  UVAR_OPTIONAL, "1st stage min start time of observation", &uvar_minStartTime1), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "maxEndTime1",  0,  UVAR_OPTIONAL, "1st stage max end time of observation",   &uvar_maxEndTime1),   &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printFstat1",  0,  UVAR_OPTIONAL, "Print 1st stage Fstat vectors", &uvar_printFstat1), &status);

  /* developer user variables */
  LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed", 0, UVAR_DEVELOPER, "RngMed block size", &uvar_blocksRngMed), &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "SSBprecision", 0, UVAR_DEVELOPER, "Precision for SSB transform.", &uvar_SSBprecision),    &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printMaps",    0, UVAR_DEVELOPER, "Print Hough maps -- for debugging", &uvar_printMaps), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printGrid",    0, UVAR_DEVELOPER, "Print Hough fine grid -- for debugging", &uvar_printGrid), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "dumpLUT",      0, UVAR_DEVELOPER, "Print Hough look-up-tables -- for debugging", &uvar_dumpLUT), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "validateLUT",  0, UVAR_DEVELOPER, "Validate Hough look-up-tables -- for debugging", &uvar_validateLUT), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printStats",   0, UVAR_DEVELOPER, "Print Hough map statistics", &uvar_printStats), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "Dterms",       0, UVAR_DEVELOPER, "No.of terms to keep in Dirichlet Kernel", &uvar_Dterms ), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "skyPointIndex",0, UVAR_DEVELOPER, "Only analyze this skypoint in grid", &uvar_skyPointIndex ), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dopplerMax",   0, UVAR_DEVELOPER, "Max Doppler shift",  &uvar_dopplerMax), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "sftUpsampling",0, UVAR_DEVELOPER, "Upsampling factor for fast LALDemod",  &uvar_sftUpsampling), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "useToplist1",  0, UVAR_DEVELOPER, "Use toplist for 1st stage candidates?", &uvar_useToplist1 ), &status);
  LAL_CALL( LALRegisterREALUserVar (  &status, "df1dotRes",    0, UVAR_DEVELOPER, "Resolution in residual fdot values (default=df1dot/nf1dotRes)", &uvar_df1dotRes), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "GPUready",     0, UVAR_OPTIONAL,  "Use single-precision 'GPU-ready' core routines", &uvar_GPUready), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "correctFreqs", 0, UVAR_DEVELOPER, "Correct candidate output frequencies (ie fix bug #147). Allows reproducing 'historical results'", &uvar_correctFreqs), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "device",       0, UVAR_DEVELOPER, "GPU device id", &uvar_gpu_device ), &status);
  LAL_CALL ( LALRegisterBOOLUserVar(  &status, "version",     'V', UVAR_SPECIAL,  "Output version information", &uvar_version), &status);

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

  if(uvar_gpu_device >= 0)
    gpu_device_id = uvar_gpu_device;

  /* some basic sanity checks on user vars */
  if ( (uvar_method != 0) && (uvar_method != 1) && (uvar_method != -1)) {
    fprintf(stderr, "Invalid method....must be 0, 1 or -1\n");
    return( HIERARCHICALSEARCH_EBAD );
  }

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

  /* create toplist -- semiCohToplist has the same structure
     as a fstat candidate, so treat it as a fstat candidate */
  create_houghFStat_toplist(&semiCohToplist, uvar_nCand1);

  /* write the log file */
  if ( uvar_log )
    {
      fnamelog = LALCalloc( alloc_len = strlen(uvar_fnameout) + 1 + 4, sizeof(CHAR) );
      XLAL_CHECK ( fnamelog != NULL, XLAL_ENOMEM, "Failed to allocated memory LALCalloc(1, %d)\n", alloc_len );
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

  /* initialize ephemeris info */

  edat = LALCalloc(1, sizeof(EphemerisData));
  XLAL_CHECK ( edat != NULL, XLAL_ENOMEM, "Error allocating memory LALCalloc ( 1, %d )\n", sizeof(EphemerisData) );

  (*edat).ephiles.earthEphemeris = uvar_ephemE;
  (*edat).ephiles.sunEphemeris = uvar_ephemS;

  /* read in ephemeris data */
  LAL_CALL( LALInitBarycenter( &status, edat), &status);

  XLALGPSSetREAL8(&minStartTimeGPS, uvar_minStartTime1);
  XLALGPSSetREAL8(&maxEndTimeGPS, uvar_maxEndTime1);

  /* create output Hough file for writing if requested by user */
  if ( uvar_printCand1 )
    {
      fnameSemiCohCand = LALCalloc( alloc_len = strlen(uvar_fnameout) + 1, sizeof(CHAR) );
      XLAL_CHECK ( fnameSemiCohCand != NULL, XLAL_ENOMEM, "Failed to allocate memory LALCalloc ( 1, %d )\n", alloc_len );

      strcpy(fnameSemiCohCand, uvar_fnameout);
    }


  if ( uvar_printFstat1 )
    {
      const CHAR *append = "_fstatVec1.dat";
      fnameFstatVec1 = LALCalloc( alloc_len = strlen(uvar_fnameout) + strlen(append) + 1, sizeof(CHAR) );
      XLAL_CHECK ( fnameFstatVec1 != NULL, XLAL_ENOMEM, "Failed to allocate memory LALCalloc ( 1, %d )\n", alloc_len );
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
  INIT_MEM(spinRange_Temp);

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
  LAL_CALL( SetUpSFTs( &status, &stackMultiSFT, &stackMultiNoiseWeights, &stackMultiDetStates, &usefulParams), &status);
  LogPrintfVerbatim (LOG_DEBUG, "done\n");

  /* some useful params computed by SetUpSFTs */
  tStack = usefulParams.tStack;
  tObs = usefulParams.tObs;
  nStacks = usefulParams.nStacks;
  /* currently unused: LIGOTimeGPS tStartGPS = usefulParams.tStartGPS; */
  midTstack = usefulParams.midTstack;
  startTstack = usefulParams.startTstack;
  tMidGPS = usefulParams.spinRange_midTime.refTime;
  refTimeGPS = usefulParams.spinRange_refTime.refTime;
  LogPrintf(LOG_DETAIL, "GPS Reference Time = %d\n", refTimeGPS.gpsSeconds);


  if ( uvar_sftUpsampling > 1 )
    {
      LogPrintf (LOG_DEBUG, "Upsampling SFTs by factor %d ... ", uvar_sftUpsampling );
      for (k = 0; k < nStacks; k++) {
	LAL_CALL ( upsampleMultiSFTVector ( &status, stackMultiSFT.data[k], uvar_sftUpsampling, 16 ), &status );
      }
      LogPrintfVerbatim (LOG_DEBUG, "done.\n");
    }

  /*------- set frequency and spindown resolutions and ranges for Fstat and semicoherent steps -----*/

  /* set Fstat calculation frequency resolution
     -- default is 1/tstack */
  if ( LALUserVarWasSet(&uvar_dFreq) ) {
    dFreqStack = uvar_dFreq;
  }
  else {
    dFreqStack = 1.0/tStack;
  }

  /* set Fstat spindown resolution
     -- default is 1/Tstack^2 */
  if ( LALUserVarWasSet(&uvar_df1dot) ) {
    df1dot = uvar_df1dot;
  }
  else {
    df1dot = 1.0/(tStack*tStack);
  }

  /* number of coarse grid spindown values */
  nf1dot = (UINT4)( usefulParams.spinRange_midTime.fkdotBand[1]/ df1dot + 1e-6) + 1;

  /* set number of residual spindowns for semi-coherent step
     --default = nStacks */
  if ( LALUserVarWasSet(&uvar_nf1dotRes) ) {
    nf1dotRes = uvar_nf1dotRes;
  }
  else {
    nf1dotRes = nStacks;
  }

  /* resolution of residual spindowns
     -- default = df1dot/nf1dotRes */
  if ( LALUserVarWasSet(&uvar_df1dotRes) ) {
    df1dotRes = uvar_df1dotRes;
  }
  else {
    df1dotRes = df1dot/nf1dotRes;
  }


  /*---------- compute noise weight for each stack and initialize total weights vector
     -- for debugging purposes only -- we will calculate noise and AM weights later ----------*/
  if (lalDebugLevel) {
    LAL_CALL( ComputeStackNoiseWeights( &status, &weightsNoise, &stackMultiNoiseWeights), &status);
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
  LogPrintf(LOG_DETAIL, "1st stage params: Nstacks = %d,  Tstack = %.0fsec, dFreq = %eHz, Tobs = %.0fsec\n",
	    nStacks, tStack, dFreqStack, tObs);
  for (k = 0; k < nStacks; k++) {

    LogPrintf(LOG_DETAIL, "Stack %d ", k);
    if ( weightsNoise )
      LogPrintfVerbatim(LOG_DETAIL, "(GPS start time = %d, Noise weight = %f ) ", startTstack->data[k].gpsSeconds,
			weightsNoise->data[k]);

    for ( j = 0; j < (INT4)stackMultiSFT.data[k]->length; j++) {

      INT4 tmpVar = stackMultiSFT.data[k]->data[j]->length;
      LogPrintfVerbatim(LOG_DETAIL, "%s: %d  ", stackMultiSFT.data[k]->data[j]->data[0].name, tmpVar);

    } /* loop over ifos */
    LogPrintfVerbatim(LOG_DETAIL, "\n");
  } /* loop over stacks */

  /** fix frequency-quantization offset bug #147, by applying a fixed frequency-correction to all candidates */
  REAL8 FreqBugCorrection = 0;
  if ( uvar_correctFreqs )
    {
      FreqBugCorrection = uvar_Freq - uvar_dFreq * floor ( uvar_Freq / uvar_dFreq + 0.5 );
      printf ("Applying frequency-correction shift of %.9g Hz \n", FreqBugCorrection );
    }
  else
    {
      printf ("WARNING: turned off frequency-shift bug correction! I hope you know what you're doing!\n");
    }

  /*---------- set up F-statistic calculation stuff ---------*/

  /* set reference time for calculating Fstatistic */
  /* thisPoint.refTime = tStartGPS; */
  thisPoint.refTime = tMidGPS;
  /* binary orbit and higher spindowns not considered */
  thisPoint.orbit = NULL;
  INIT_MEM ( thisPoint.fkdot );

  /* some compute F params */
  CFparams.Dterms = uvar_Dterms;
  CFparams.SSBprec = uvar_SSBprecision;
  CFparams.upsampling = uvar_sftUpsampling;


  /* set up some semiCoherent parameters */
  semiCohPar.useToplist = uvar_useToplist1;
  semiCohPar.tsMid = midTstack;
  /* semiCohPar.refTime = tStartGPS; */
  semiCohPar.refTime = tMidGPS;
  /* calculate detector velocity and positions */
  LAL_CALL( GetStackVelPos( &status, &velStack, &posStack, &stackMultiDetStates), &status);
  semiCohPar.vel = velStack;
  semiCohPar.pos = posStack;

  semiCohPar.outBaseName = uvar_fnameout;
  semiCohPar.pixelFactor = uvar_pixelFactor;
  semiCohPar.nfdot = nf1dotRes;
  semiCohPar.dfdot = df1dotRes;

  /* allocate memory for Hough candidates */
  semiCohCandList.length = uvar_nCand1;
  /* reference time for candidates */
  /* semiCohCandList.refTime = tStartGPS; */
  semiCohCandList.refTime = tMidGPS;
  semiCohCandList.nCandidates = 0; /* initialization */
  semiCohCandList.list = LALCalloc( 1, alloc_len = semiCohCandList.length * sizeof(SemiCohCandidate));
  XLAL_CHECK ( semiCohCandList.list != NULL, XLAL_ENOMEM, "Error allocating memory LALCalloc ( 1, %d )\n", alloc_len );

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

  LogPrintf(LOG_DETAIL,"Hough/Stackslide patchsize is %frad x %frad\n",
	    semiCohPar.patchSizeX, semiCohPar.patchSizeY);

  /* allocate some fstat memory */
  fstatVector.length = nStacks;
  fstatVector.data = NULL;
  fstatVector.data = LALCalloc( 1, alloc_len = nStacks * sizeof(REAL4FrequencySeries));
  XLAL_CHECK ( fstatVector.data != NULL, XLAL_ENOMEM, "Error allocating memory LALCalloc ( 1, %d )\n", alloc_len );

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
  scanInit.Detector = &(stackMultiDetStates.data[0]->data[0]->detector);
  scanInit.ephemeris = edat;
  scanInit.skyGridFile = uvar_skyGridFile;
  scanInit.skyRegionString = LALCalloc(1, alloc_len = strlen(uvar_skyRegion)+1);
  XLAL_CHECK ( scanInit.skyRegionString != NULL, XLAL_ENOMEM, "Failed to allocate memory LALCalloc ( 1, %d )\n", alloc_len );
  strcpy (scanInit.skyRegionString, uvar_skyRegion);

  scanInit.numSkyPartitions = uvar_numSkyPartitions;
  scanInit.partitionIndex = uvar_partitionIndex;

  scanInit.Freq = usefulParams.spinRange_midTime.fkdot[0] +  usefulParams.spinRange_midTime.fkdotBand[0];

  /* initialize skygrid  */
  LogPrintf(LOG_DETAIL, "Setting up coarse sky grid...");
  LAL_CALL ( InitDopplerSkyScan ( &status, &thisScan, &scanInit), &status);
  LogPrintfVerbatim(LOG_DETAIL, "done\n");

  /*----- start main calculations by going over coarse grid points and
          selecting candidates and following up --------*/

  /* loop over skygrid points */
  skyGridCounter = 0;

  XLALNextDopplerSkyPos(&dopplerpos, &thisScan);

  /* "spool forward" if we found a checkpoint */
  {
    UINT4 count = 0;
    GET_CHECKPOINT(semiCohToplist, &count, thisScan.numSkyGridPoints, fnameSemiCohCand, NULL);
    for(skyGridCounter = 0; skyGridCounter < count; skyGridCounter++)
      XLALNextDopplerSkyPos(&dopplerpos, &thisScan);
  }

  /* spool forward if uvar_skyPointIndex is set
     ---This probably doesn't make sense when checkpointing is turned on */
  if ( LALUserVarWasSet(&uvar_skyPointIndex)) {
    UINT4 count = uvar_skyPointIndex;
    for(skyGridCounter = 0; (skyGridCounter < count)&&(thisScan.state != STATE_FINISHED) ; skyGridCounter++)
      XLALNextDopplerSkyPos(&dopplerpos, &thisScan);
  }


#ifdef SKYPOS_PRECISION
  LogPrintf(LOG_DEBUG, "SKYPOS_PRECISION: %15f (0x%x)\n", (REAL4)SKYPOS_PRECISION,(INT8)SKYPOS_PRECISION);
#endif

  LogPrintf(LOG_DEBUG, "Total skypoints = %d. Progress: ", thisScan.numSkyGridPoints);

  INITIALIZE_COPROCESSOR_DEVICE

#ifdef OUTPUT_TIMING
    clock0 = time(NULL);
#endif


  while(thisScan.state != STATE_FINISHED)
    {
      UINT4 ifdot;  /* counter for spindown values */
      SkyPosition skypos;
      ComputeFBufferREAL4V cfvBuffer = empty_ComputeFBufferREAL4V;

      /* if (skyGridCounter == 25965) { */


#ifdef SKYPOS_PRECISION
      /* reduce precision of sky position */
      dopplerpos.Alpha = (INT4)(dopplerpos.Alpha * SKYPOS_PRECISION) / (REAL4)SKYPOS_PRECISION;
      dopplerpos.Delta = (INT4)(dopplerpos.Delta * SKYPOS_PRECISION) / (REAL4)SKYPOS_PRECISION;
#endif

      SHOW_PROGRESS(dopplerpos.Alpha,dopplerpos.Delta,
		    skyGridCounter,thisScan.numSkyGridPoints,
		    uvar_Freq, uvar_FreqBand);

      LogPrintfVerbatim(LOG_DEBUG, "%d, ", skyGridCounter );

      /*------------- calculate F statistic for each stack --------------*/

      /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
      thisPoint.Alpha = dopplerpos.Alpha;
      thisPoint.Delta = dopplerpos.Delta;

      /* get amplitude modulation weights */
      skypos.longitude = thisPoint.Alpha;
      skypos.latitude = thisPoint.Delta;
      skypos.system = COORDINATESYSTEM_EQUATORIAL;

      semiCohPar.alpha = thisPoint.Alpha;
      semiCohPar.delta = thisPoint.Delta;

      /* initialize weights to unity */
      LAL_CALL( LALHOUGHInitializeWeights( &status, weightsV), &status);

      if (uvar_useWeights) {
	LAL_CALL( ComputeStackNoiseAndAMWeights( &status, weightsV, &stackMultiNoiseWeights,
						 &stackMultiDetStates, skypos), &status);
      }
      semiCohPar.weightsV = weightsV;

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
	ComputeNumExtraBins(&status, &semiCohPar, 0, freqHighest, dFreqStack);

	/* conservative estimate */
	/* extraBinsSky = (UINT4)(0.5 * LAL_SQRT2 * VTOT   */
	/*      * (usefulParams.spinRange_midTime.fkdot[0] + usefulParams.spinRange_midTime.fkdotBand[0])  */
	/*      * HSMAX(semiCohPar.patchSizeX, semiCohPar.patchSizeY) / dFreqStack) + 10; */

	/* extra bins due to fdot is maximum number of frequency bins drift that can be
	   caused by the residual spindown.  The reference time for the spindown is the midtime,
	   so relevant interval is Tobs/2 and largest possible value of residual spindown is
	   (number of residual spindowns -1)*resolution in residual spindowns */
	extraBinsfdot = (UINT4)(tObs * (nf1dotRes - 1) * df1dotRes / dFreqStack + 0.5);

	semiCohPar.extraBinsFstat += extraBinsfdot;

	/* allocate fstat memory */
	binsFstatSearch = (UINT4)(usefulParams.spinRange_midTime.fkdotBand[0]/dFreqStack + 1e-6) + 1;
	binsFstat1 = binsFstatSearch + 2*semiCohPar.extraBinsFstat;

	for (k = 0; k < nStacks; k++) {
	  /* careful--the epoch here is not the reference time for f0! */
	  fstatVector.data[k].epoch = startTstack->data[k];
	  fstatVector.data[k].deltaF = dFreqStack;
	  fstatVector.data[k].f0 = usefulParams.spinRange_midTime.fkdot[0] - semiCohPar.extraBinsFstat * dFreqStack;
	  if (fstatVector.data[k].data == NULL) {
	    fstatVector.data[k].data = LALCalloc( 1, alloc_len = sizeof(REAL4Sequence));
	    XLAL_CHECK( fstatVector.data[k].data != NULL, XLAL_ENOMEM, "Failed to allocate memory LALCalloc ( 1, %d )\n", alloc_len );

	    fstatVector.data[k].data->length = binsFstat1;
	    fstatVector.data[k].data->data = LALCalloc( 1, alloc_len = binsFstat1 * sizeof(REAL4));
	    XLAL_CHECK( fstatVector.data[k].data->data != NULL, XLAL_ENOMEM, "Failed to allocate memory LALCalloc ( 1, %d )\n", alloc_len );
	  }
	  else {
	    fstatVector.data[k].data = LALRealloc( fstatVector.data[k].data, sizeof(REAL4Sequence));
	    XLAL_CHECK( fstatVector.data[k].data != NULL, XLAL_ENOMEM, "Failed to re-allocate memory LALRealloc ( %d )\n", sizeof(REAL4Sequence) );

	    fstatVector.data[k].data->length = binsFstat1;
	    fstatVector.data[k].data->data = LALRealloc( fstatVector.data[k].data->data, alloc_len = binsFstat1 * sizeof(REAL4));
	    XLAL_CHECK( fstatVector.data[k].data->data != NULL, XLAL_ENOMEM, "Failed to re-allocate memory LALRealloc ( %d )\n", alloc_len );
	  }
	} /* loop over stacks */

	REARRANGE_SFT_DATA

      } /* fstat memory allocation block */

      cfvBuffer = empty_ComputeFBufferREAL4V;

	/* loop over fdot values
	   -- spindown range and resolutionhas been set earlier */
      for ( ifdot = 0; ifdot < nf1dot; ifdot++)
	{

	  LogPrintf(LOG_DETAIL, "Analyzing %d/%d Coarse sky grid points and %d/%d spindown values\n",
		    skyGridCounter, thisScan.numSkyGridPoints, ifdot+1, nf1dot);

	  SHOW_PROGRESS(dopplerpos.Alpha,dopplerpos.Delta, skyGridCounter + (REAL4)ifdot / (REAL4)nf1dot,
			thisScan.numSkyGridPoints, uvar_Freq, uvar_FreqBand);


	  /* calculate the Fstatistic for each stack*/
	  LogPrintf(LOG_DETAIL, "Starting Fstat calculation for each stack...");
          if ( !uvar_GPUready )
            {
              for ( k = 0; k < nStacks; k++)
                {
                  /* set spindown value for Fstat calculation */
                  thisPoint.fkdot[1] = usefulParams.spinRange_midTime.fkdot[1] + ifdot * df1dot;
                  thisPoint.fkdot[0] = fstatVector.data[k].f0;
                  /* thisPoint.fkdot[0] = usefulParams.spinRange_midTime.fkdot[0]; */

                  /* this is the most costly function. We here allow for using an architecture-specific optimized
                     function from e.g. a local file instead of the standard ComputeFStatFreqBand() from LAL */
                  LAL_CALL( COMPUTEFSTATFREQBAND ( &status, fstatVector.data + k, &thisPoint,
                                                   stackMultiSFT.data[k], stackMultiNoiseWeights.data[k],
                                                   stackMultiDetStates.data[k], &CFparams), &status);
                } /* for k < nStacks */
            }
          else	/* run "GPU-ready" REAL4 version aiming at maximum parallelism for GPU optimization */
            {
              /* set spindown value for Fstat calculation */
              thisPoint.fkdot[1] = usefulParams.spinRange_midTime.fkdot[1] + ifdot * df1dot;
              thisPoint.fkdot[0] = fstatVector.data[0].f0;
              /* thisPoint.fkdot[0] = usefulParams.spinRange_midTime.fkdot[0]; */

              /* this is the most costly function. We here allow for using an architecture-specific optimized
               * function from e.g. a local file instead of the standard ComputeFStatFreqBand() from LAL */
              if ( XLALComputeFStatFreqBandVector ( &fstatVector, &thisPoint, &stackMultiSFT, &stackMultiNoiseWeights,
                                                    &stackMultiDetStates, CFparams.Dterms, &cfvBuffer) != XLAL_SUCCESS )
                {
                  LogPrintf (LOG_CRITICAL, "main(): XLALComputeFStatFreqBandVector() failed with errno = %d\n", xlalErrno );
                  return XLAL_EFUNC;
                }

            } /* if GPUready */

          LogPrintfVerbatim(LOG_DETAIL, "done\n");

	  /* print fstat vector if required -- mostly for debugging */
	  if ( uvar_printFstat1 )
	    {
	      for (k = 0; k < nStacks; k++)
		LAL_CALL( PrintFstatVec ( &status, fstatVector.data + k, fpFstat1, &thisPoint, refTimeGPS, k+1), &status);
	    }


	  /*--------------- get candidates from a semicoherent search ---------------*/

	  /* the input to this section is the set of fstat vectors fstatVector and the
	     parameters semiCohPar. The output is the list of candidates in semiCohCandList */


	  /* set spindown for Hough grid -- same as for Fstat calculation */
	  semiCohPar.fdot = thisPoint.fkdot[1];

	  /* the hough option */
	  /* select peaks */
	  if (uvar_method == 0) {

	    LogPrintf(LOG_DETAIL, "Starting Hough calculation...\n");
	    sumWeightSquare = 0.0;
	    for ( k = 0; k < nStacks; k++)
	      sumWeightSquare += weightsV->data[k] * weightsV->data[k];

	    /* set number count threshold based on significance threshold */
	    meanN = nStacks * alphaPeak;
	    sigmaN = sqrt(sumWeightSquare * alphaPeak * (1.0 - alphaPeak));
	    semiCohPar.threshold = uvar_threshold1*sigmaN + meanN;
	    LogPrintf(LOG_DETAIL, "Expected mean number count=%f, std=%f, threshold=%f\n",
		      meanN, sigmaN, semiCohPar.threshold);

	    /* convert fstat vector to peakgrams using the Fstat threshold */
	    LAL_CALL( FstatVectToPeakGram( &status, &pgV, &fstatVector, (REAL4)uvar_peakThrF), &status);

	    /* get candidates */
	    /* this is the second most costly function. We here allow for using architecture-specific
	       optimized functions from a local file instead of the standard ComputeFstatHoughMap()
	       below that refers to the LALHOUGH functions in LAL */
	    LAL_CALL ( COMPUTEFSTATHOUGHMAP ( &status, &semiCohCandList, &pgV, &semiCohPar), &status);

            /* ----- now apply frequency-shift correction to fix bug #147 to all candidates returned */
            INT4 iCand;
            for ( iCand = 0; iCand < semiCohCandList.nCandidates; iCand ++ )
              semiCohCandList.list[iCand].freq += FreqBugCorrection;
            /* ------------------------------------------------------------ */

	    /* free peakgrams -- we don't need them now because we have the Hough maps */
	    for (k=0; k<nStacks; k++)
	      LALFree(pgV.pg[k].peak);
	    LALFree(pgV.pg);
	    LogPrintf(LOG_DETAIL, "...finished Hough calculation\n");

	      /* end hough */
	  } else if ( uvar_method == 1 ) {
	    /* --- stackslide option --------*/
	    LogPrintf(LOG_DETAIL, "Starting StackSlide calculation...\n");
	    /* 12/18/06 gm; use threshold from command line as threshold on stackslide sum of F-stat values */
	    semiCohPar.threshold = uvar_threshold1;
	    LAL_CALL( StackSlideVecF( &status, &semiCohCandList, &fstatVector, &semiCohPar), &status);
	    LogPrintf(LOG_DETAIL, "...finished StackSlide calculation\n");
	  }

	  /* print candidates if desired */
	  /* 	  if ( uvar_printCand1 ) { */
	  /* 	    LAL_CALL ( PrintSemiCohCandidates ( &status, &semiCohCandList, fpSemiCoh, refTimeGPS), &status); */
	  /* 	  } */


	  if( uvar_semiCohToplist) {
	    /* this is necessary here, because GetSemiCohToplist() might set
	       a checkpoint that needs some information from here */
	    SHOW_PROGRESS(dopplerpos.Alpha,dopplerpos.Delta,
			  skyGridCounter,thisScan.numSkyGridPoints,
			  uvar_Freq, uvar_FreqBand);

	    LogPrintf(LOG_DETAIL, "Selecting toplist from semicoherent candidates\n");
	    LAL_CALL( GetSemiCohToplist(&status, semiCohToplist, &semiCohCandList, meanN, sigmaN), &status);
	  }

	} /* end loop over coarse grid fdot values */

      XLALEmptyComputeFBufferREAL4V ( &cfvBuffer );

      /* continue forward till the end if uvar_skyPointIndex is set
	 ---This probably doesn't make sense when checkpointing is turned on */
      if ( LALUserVarWasSet(&uvar_skyPointIndex) )
	{
	  while(thisScan.state != STATE_FINISHED) {
	    skyGridCounter++;
	    XLALNextDopplerSkyPos(&dopplerpos, &thisScan);
	  }
	}
      else
	{
	  skyGridCounter++;

	  /* this is necessary here, because the checkpoint needs some information from here */
	  SHOW_PROGRESS(dopplerpos.Alpha,dopplerpos.Delta,
			skyGridCounter,thisScan.numSkyGridPoints,
			uvar_Freq, uvar_FreqBand);

	  SET_CHECKPOINT;

	  XLALNextDopplerSkyPos( &dopplerpos, &thisScan );
	}

    } /* end while loop over 1st stage coarse skygrid */

#ifdef OUTPUT_TIMING
  {
    time_t tau = time(NULL) - clock0;
    UINT4 Nrefine = nSkyRefine * nf1dotRes;
    FILE *timing_fp = fopen ( "HS_timing.dat", "ab" );
    fprintf ( timing_fp, "%d 	%d 	%d 	%d 	%d 	%d 	%d 	%d\n",
	      thisScan.numSkyGridPoints, nf1dot, binsFstatSearch, 2 * semiCohPar.extraBinsFstat, nSFTs, nStacks, Nrefine, tau );
    fclose ( timing_fp );
  }
#endif

  UNINITIALIZE_COPROCESSOR_DEVICE

  LogPrintfVerbatim ( LOG_DEBUG, " done.\n");

  LogPrintf(LOG_DETAIL, "Finished analysis and now printing results and cleaning up...");

  LogPrintfVerbatim ( LOG_DEBUG, "Writing output ...\n");

#if (!HS_CHECKPOINTING)
  /* print 1st stage candidates */
  {
    if (!(fpSemiCoh = fopen(fnameSemiCohCand, "wb")))
      {
	LogPrintf ( LOG_CRITICAL, "Unable to open output-file '%s' for writing.\n", fnameSemiCohCand);
	return HIERARCHICALSEARCH_EFILE;
      }
    if ( uvar_printCand1 && uvar_semiCohToplist ) {
      fprintf ( fpSemiCoh, "%%%%  Freq            Alpha              Delta              f1dot                 HoughFStat        AlphaBest          DeltaBest          MeanSig    VarSig\n");
      sort_houghFStat_toplist(semiCohToplist);
      if ( write_houghFStat_toplist_to_fp( semiCohToplist, fpSemiCoh, NULL) < 0)
	fprintf( stderr, "Error in writing toplist to file\n");
      /*     LAL_CALL( AppendFstatCandidates( &status, &fStatCand, fpFstat), &status); */
      if (fprintf(fpSemiCoh,"%%DONE\n") < 0)
	fprintf(stderr, "Error writing end marker\n");
      fclose(fpSemiCoh);
    }
  }
#else
  write_and_close_checkpointed_file();
#endif

  LogPrintfVerbatim ( LOG_DEBUG, " done.\n");

  /*------------ free all remaining memory -----------*/

  /* free memory */

  if ( uvar_printCand1 )
    {
      LALFree(fnameSemiCohCand);
    }

  if ( VCSInfoString ) XLALFree ( VCSInfoString );

  if ( uvar_printFstat1 )
    {
      fclose(fpFstat1);
      LALFree( fnameFstatVec1 );
    }

  /* free first stage memory */
  for ( k = 0; k < nStacks; k++) {
    LAL_CALL( LALDestroyMultiSFTVector ( &status, stackMultiSFT.data + k), &status);
    LAL_CALL( LALDestroyMultiNoiseWeights ( &status, stackMultiNoiseWeights.data + k), &status);
    XLALDestroyMultiDetectorStateSeries ( stackMultiDetStates.data[k] );
  }
  LALFree(stackMultiSFT.data);
  LALFree(stackMultiNoiseWeights.data);
  LALFree(stackMultiDetStates.data);


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


  /* free Vel/Pos vectors and ephemeris */
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);
  LAL_CALL( LALDDestroyVectorSequence (&status,  &velStack), &status);
  LAL_CALL( LALDDestroyVectorSequence (&status,  &posStack), &status);

  if (weightsNoise) {
    LAL_CALL( LALDDestroyVector ( &status, &weightsNoise ), &status);
  }
  LAL_CALL( LALDDestroyVector ( &status, &weightsV ), &status);

  /* free dopplerscan stuff */
  LAL_CALL ( FreeDopplerSkyScan(&status, &thisScan), &status);
  if ( scanInit.skyRegionString )
    LALFree ( scanInit.skyRegionString );


  /* free candidates */
  LALFree(semiCohCandList.list);
  free_houghFStat_toplist(&semiCohToplist);

  LAL_CALL (LALDestroyUserVars(&status), &status);

  LALCheckMemoryLeaks();

  LogPrintfVerbatim(LOG_DETAIL, "done\n");

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
  UINT4 k,numSFTby2;
  LIGOTimeGPS tStartGPS, tEndGPS, refTimeGPS, tMidGPS;
  SFTCatalogSequence catalogSeq;

  REAL8 doppWings, fMin, fMax;
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
  extraBins = HSMAX ( in->blocksRngMed/2 + 1, in->Dterms );

  fMin = freqLo - doppWings - extraBins * deltaFsft;
  fMax = freqHi + doppWings + extraBins * deltaFsft;

  /* finally memory for stack of multi sfts */
  stackMultiSFT->length = in->nStacks;
  stackMultiSFT->data = LALCalloc(1, alloc_len = in->nStacks * sizeof(MultiSFTVector *));
  if ( stackMultiSFT->data == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  stackMultiNoiseWeights->length = in->nStacks;
  stackMultiNoiseWeights->data = LALCalloc(1, alloc_len = in->nStacks * sizeof(MultiNoiseWeights *));
  if ( stackMultiNoiseWeights->data == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  stackMultiDetStates->length = in->nStacks;
  stackMultiDetStates->data = LALCalloc(1, alloc_len = in->nStacks * sizeof(MultiDetectorStateSeries *));
  if ( stackMultiDetStates->data == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  /* loop over stacks and read sfts */
  for (k = 0; k < in->nStacks; k++) {

    MultiPSDVector *psd = NULL;

    /* load the sfts */
    TRY( LALLoadMultiSFTs ( status->statusPtr, stackMultiSFT->data + k,  catalogSeq.data + k,
			    fMin, fMax ), status);

    /* normalize sfts and compute noise weights and detector state */
    TRY( LALNormalizeMultiSFTVect ( status->statusPtr, &psd, stackMultiSFT->data[k],
				    in->blocksRngMed ), status );

    TRY( LALComputeMultiNoiseWeights  ( status->statusPtr, stackMultiNoiseWeights->data + k,
					psd, in->blocksRngMed, 0 ), status );

    TRY ( LALGetMultiDetectorStates ( status->statusPtr, stackMultiDetStates->data + k,
				      stackMultiSFT->data[k], in->edat ), status );

    TRY ( LALDestroyMultiPSDVector ( status->statusPtr, &psd ), status );

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
  nStacks = stackMultiSFT->length;
  nSFTs = 0;
  for ( k = 0; k < nStacks; k ++ )
    {
      UINT4 X;
      for ( X=0; X < stackMultiSFT->data[k]->length; X ++ )
	nSFTs += stackMultiSFT->data[k]->data[X]->length;
    } /* for k < stacks */
#endif




  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* SetUpSFTs */



/** Function for calculating Hough Maps and candidates.
 *
    This function takes a peakgram as input. This peakgram was constructed
    by setting a threshold on a sequence of Fstatistic vectors.  The function
    produces a Hough map in the sky for each value of the frequency and spindown.
    The Hough nummber counts are then used to select candidates in
    parameter space to be followed up in a more refined search.
    This uses DriveHough_v3.c as a prototype suitably modified to work
    on demodulated data instead of SFTs.
*/
void ComputeFstatHoughMap(LALStatus *status,		/**< pointer to LALStatus structure */
			  SemiCohCandidateList  *out,   /**< Candidates from thresholding Hough number counts */
			  HOUGHPeakGramVector *pgV, 	/**< HOUGHPeakGramVector obtained after thresholding Fstatistic vectors */
			  SemiCoherentParams *params	/**< pointer to HoughParams -- parameters for calculating Hough maps */
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

  toplist_t *houghToplist;

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
  lutV.lut = LALCalloc(1, alloc_len = nStacks*sizeof(HOUGHptfLUT));
  if ( lutV.lut == NULL ) {
    XLALPrintError ("Failed to LALCalloc(1,%d)\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
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
  phmdVS.phmd = LALCalloc( 1, alloc_len = phmdVS.length * phmdVS.nfSize *sizeof(HOUGHphmd));
  if ( phmdVS.phmd == NULL ) {
    XLALPrintError ("Failed to LALCalloc(1,%d)\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  /* residual spindown trajectory */
  freqInd.deltaF = deltaF;
  freqInd.length = nStacks;
  freqInd.data = LALCalloc(1, alloc_len = nStacks*sizeof(UINT8));
  if ( freqInd.data == NULL ) {
    XLALPrintError ("Failed to LALCalloc(1,%d)\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  /* resolution in space of residual spindowns */
  ht.dFdot.length = 1;
  ht.dFdot.data = LALCalloc( 1, alloc_len = ht.dFdot.length * sizeof(REAL8));
  if ( ht.dFdot.data == NULL ) {
    XLALPrintError ("Failed to LALCalloc(1,%d)\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  /* the residual spindowns */
  ht.spinRes.length = 1;
  ht.spinRes.data = LALCalloc( 1, alloc_len = ht.spinRes.length*sizeof(REAL8));
  if ( ht.spinRes.data == NULL ) {
    XLALPrintError ("Failed to LALCalloc(1,%d)\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  /* the residual spindowns */
  ht.spinDem.length = 1;
  ht.spinDem.data = LALCalloc( 1, alloc_len = ht.spinRes.length*sizeof(REAL8));
  if ( ht.spinDem.data == NULL ) {
    XLALPrintError ("Failed to LALCalloc(1,%d)\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  /* the demodulation params */
  parDem.deltaF = deltaF;
  parDem.skyPatch.alpha = alpha;
  parDem.skyPatch.delta = delta;
  parDem.spin.length = 1;
  parDem.spin.data = LALCalloc(1, alloc_len = sizeof(REAL8));
  if ( parDem.spin.data == NULL ) {
    XLALPrintError ("Failed to LALCalloc(1,%d)\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
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
    hist.data = LALCalloc(1, alloc_len = (nStacks+1)*sizeof(UINT8));
    if ( hist.data == NULL ) {
      XLALPrintError ("Failed to LALCalloc(1,%d)\n", alloc_len );
      ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
    }

    histTotal.data = LALCalloc(1, alloc_len = (nStacks+1)*sizeof(UINT8));
    if ( histTotal.data == NULL ) {
      XLALPrintError ("Failed to LALCalloc(1,%d)\n", alloc_len );
      ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
    }

    {
      UINT4   j;
      for(j=0; j< histTotal.length; ++j)
	histTotal.data[j]=0;
    }
    {
      const CHAR *append = "stats";
      fileStats = LALCalloc( alloc_len = strlen(params->outBaseName) + strlen(append) + 1, sizeof(CHAR) );
      if ( fileStats == NULL ) {
        XLALPrintError ("Failed to LALCalloc(1,%d)\n", alloc_len );
	ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
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
  /* this is not very clean -- the Fstat calculation has to know how many extra bins are needed */

  LogPrintf(LOG_DETAIL, "Freq. range analyzed by Hough = [%fHz - %fHz] (%d bins)\n",
	    fBinIni*deltaF, fBinFin*deltaF, fBinFin - fBinIni + 1);
  ASSERT ( fBinIni <= fBinFin, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );

  /* initialise number of candidates -- this means that any previous candidates
     stored in the list will be lost for all practical purposes*/
  out->nCandidates = 0;

  /* create toplist of candidates */
  if (params->useToplist) {
    create_toplist(&houghToplist, out->length, sizeof(SemiCohCandidate), smallerHough);
  }
  else {
    /* if no toplist then use number of hough maps */
    INT4 numHmaps = (fBinFin - fBinIni + 1)*phmdVS.nfSize;
    if (out->length != numHmaps) {
      out->length = numHmaps;
      out->list = LALRealloc( out->list, alloc_len = out->length * sizeof(SemiCohCandidate));
      if ( out->list == NULL ) {
        XLALPrintError ("Failed to LALRealloc( *, %d)\n", alloc_len );
	ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
      }
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

    parRes.f0Bin =  fBin;
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
    patch.xCoor = LALCalloc(1, alloc_len = xSide*sizeof(REAL8));
    if ( patch.xCoor == NULL ) {
      XLALPrintError ("Failed to LALCalloc( 1, %d)\n", alloc_len );
      ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
    }

    patch.yCoor = LALCalloc(1, alloc_len = ySide*sizeof(REAL8));
    if ( patch.yCoor == NULL ) {
      XLALPrintError ("Failed to LALCalloc( 1, %d)\n", alloc_len );
      ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
    }
    TRY( LALHOUGHFillPatchGrid( status->statusPtr, &patch, &parSize ), status );

    /*------------- other memory allocation and settings----------------- */
    for(j=0; j<lutV.length; ++j){
      lutV.lut[j].maxNBins = maxNBins;
      lutV.lut[j].maxNBorders = maxNBorders;
      lutV.lut[j].border = LALCalloc(1, alloc_len = maxNBorders*sizeof(HOUGHBorder));
      if ( lutV.lut[j].border == NULL ) {
        XLALPrintError ("Failed to LALCalloc( 1, %d)\n", alloc_len );
	ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
      }

      lutV.lut[j].bin =	LALCalloc(1, alloc_len = maxNBins*sizeof(HOUGHBin2Border));
      if ( lutV.lut[j].bin == NULL ) {
        XLALPrintError ("Failed to LALCalloc( 1, %d)\n", alloc_len );
	ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
      }

      for (i=0; i<maxNBorders; ++i){
	lutV.lut[j].border[i].ySide = ySide;
	lutV.lut[j].border[i].xPixel = LALCalloc(1, alloc_len = ySide*sizeof(COORType));
	if ( lutV.lut[j].border[i].xPixel == NULL ) {
          XLALPrintError ("Failed to LALCalloc( 1, %d)\n", alloc_len );
	  ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
	}
      }
    }

    for(j = 0; j < phmdVS.length * phmdVS.nfSize; ++j){
      phmdVS.phmd[j].maxNBorders = maxNBorders;
      phmdVS.phmd[j].leftBorderP = LALCalloc(1, alloc_len = maxNBorders*sizeof(HOUGHBorder *));
      if ( phmdVS.phmd[j].leftBorderP == NULL ) {
        XLALPrintError ("Failed to LALCalloc( 1, %d)\n", alloc_len );
	ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
      }

      phmdVS.phmd[j].rightBorderP = LALCalloc(1, alloc_len = maxNBorders*sizeof(HOUGHBorder *));
      if ( phmdVS.phmd[j].rightBorderP == NULL ) {
        XLALPrintError ("Failed to LALCalloc( 1, %d)\n", alloc_len );
	ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
      }

      phmdVS.phmd[j].ySide = ySide;
      phmdVS.phmd[j].firstColumn = NULL;
      phmdVS.phmd[j].firstColumn = LALCalloc(1, alloc_len = ySide*sizeof(UCHAR));
      if ( phmdVS.phmd[j].firstColumn == NULL ) {
        XLALPrintError ("Failed to LALCalloc( 1, %d)\n", alloc_len );
	ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
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
    ht.map   = LALCalloc(1, alloc_len = xSide*ySide*sizeof(HoughTT));
    if ( ht.map == NULL ) {
      XLALPrintError ("Failed to LALCalloc( 1, %d)\n", alloc_len );
      ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
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

	  /* get candidates */
	  if ( params->useToplist ) {
	    TRY(GetHoughCandidates_toplist( status->statusPtr, houghToplist, &ht, &patch, &parDem), status);
	  }
	  else {
	    TRY(GetHoughCandidates_threshold( status->statusPtr, out, &ht, &patch, &parDem, params->threshold), status);
	  }

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
  if ( params->useToplist ) {
    for ( k=0; k<houghToplist->elems; k++) {
      out->list[k] = *((SemiCohCandidate *)(toplist_elem(houghToplist, k)));
    }
    out->nCandidates = houghToplist->elems;
    free_toplist(&houghToplist);
  }

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

}

/** Function for selecting frequency bins from a set of Fstatistic vectors.

    Input is a vector of Fstatistic vectors.  It allocates memory
    for the peakgrams based on the frequency span of the Fstatistic vectors
    and fills tyem up by setting a threshold on the Fstatistic.  Peakgram must be
    deallocated outside the function.
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
  pgV->pg = LALCalloc( 1, alloc_len = nStacks * sizeof(HOUGHPeakGram));
  if ( pgV->pg == NULL ) {
    XLALPrintError ("Failed to LALCalloc( 1, %d)\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  upg = LALCalloc( 1, alloc_len = nSearchBins * sizeof(UCHAR));
  if ( upg == NULL ) {
    XLALPrintError ("Failed to LALCalloc( 1, %d)\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
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
    pgV->pg[k].peak = LALCalloc( 1, alloc_len = nPeaks * sizeof(INT4));
    if ( pgV->pg[k].peak == NULL ) {
      XLALPrintError ("Failed to LALCalloc( 1, %d)\n", alloc_len );
      ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
    }


    /* fill up other peakgram parameters */
    pgV->pg[k].deltaF = FstatVect->data[k].deltaF;
    f0 = FstatVect->data[k].f0;
    deltaF = FstatVect->data[k].deltaF;
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


/** \brief Breaks up input sft catalog into specified number of stacks

    Loops over elements of the catalog, assigns a bin index and
    allocates memory to the output catalog sequence appropriately.  If
    there are long gaps in the data, then some of the catalogs in the
    output catalog sequence may be of zero length.
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
  out->data = LALCalloc( 1, alloc_len = nStacksMax * sizeof(SFTCatalog));
  if ( out->data == NULL ) {
    XLALPrintError ("Failed to LALCalloc( 1, %d)\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
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
	  out->data[stackCounter].data = LALRealloc( out->data[stackCounter].data, alloc_len = length * sizeof(SFTDescriptor));
	  if ( out->data[stackCounter].data == NULL ) {
            XLALPrintError ("Failed to LALRealloc( *, %d)\n", alloc_len );
	    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
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
	  out->data[stackCounter].data = LALRealloc( out->data[stackCounter].data, alloc_len = sizeof(SFTDescriptor));
	  if ( out->data[stackCounter].data == NULL ) {
            XLALPrintError ("Failed to LALRealloc( *, %d)\n", alloc_len );
	    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
	  }

	  out->data[stackCounter].data[0] = in->data[j];
	} /* if new stack */

    } /* loop over sfts */

  /* realloc catalog sequence length to actual number of stacks */
  out->length = stackCounter + 1;
  out->data = LALRealloc( out->data, alloc_len = (stackCounter+1) * sizeof(SFTCatalog) );
  if ( out->data == NULL ) {
    XLALPrintError ("Failed to LALRealloc( *, %d)\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
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
  pg.peak = LALCalloc(1, alloc_len = pg.length*sizeof(INT4) );
  if ( pg.peak == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  phmd.fBin = f0Bin;
  phmd.maxNBorders = maxNBorders;
  phmd.leftBorderP = LALMalloc( alloc_len = maxNBorders*sizeof(HOUGHBorder *));
  if ( phmd.leftBorderP == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  phmd.rightBorderP = LALMalloc( alloc_len = maxNBorders*sizeof(HOUGHBorder *));
  if ( phmd.rightBorderP == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  phmd.ySide = ySide;
  phmd.firstColumn = LALMalloc( alloc_len = ySide*sizeof(UCHAR));
  if ( phmd.firstColumn == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  ht.xSide = xSide;
  ht.ySide = ySide;
  ht.map   = LALMalloc( alloc_len = xSide*ySide*sizeof(HoughTT));
  if ( ht.map == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  hd.xSide = xSide;
  hd.ySide = ySide;
  hd.map   = LALMalloc( alloc_len = (xSide+1)*ySide*sizeof(HoughDT));
  if ( hd.map == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

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
  pg.peak = LALCalloc(1, alloc_len = pg.length*sizeof(INT4));
  if ( pg.peak == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  phmd.fBin = f0Bin;
  phmd.maxNBorders = maxNBorders;
  phmd.leftBorderP  = LALMalloc( alloc_len = maxNBorders*sizeof(HOUGHBorder *));
  if ( phmd.leftBorderP == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  phmd.rightBorderP = LALMalloc( alloc_len = maxNBorders*sizeof(HOUGHBorder *));
  if ( phmd.rightBorderP == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  phmd.ySide = ySide;
  phmd.firstColumn = LALMalloc( alloc_len = ySide*sizeof(UCHAR));
  if ( phmd.firstColumn == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  ht.xSide = xSide;
  ht.ySide = ySide;
  ht.map   = LALMalloc( alloc_len = xSide*ySide*sizeof(HoughTT));
  if ( ht.map == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  hd.xSide = xSide;
  hd.ySide = ySide;
  hd.map   = LALMalloc( alloc_len = (xSide+1)*ySide*sizeof(HoughDT));
  if ( hd.map == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

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



/** Get Hough candidates as a toplist using a fixed threshold*/
void GetHoughCandidates_threshold(LALStatus            *status,
				  SemiCohCandidateList *out,
				  HOUGHMapTotal        *ht,
				  HOUGHPatchGrid       *patch,
				  HOUGHDemodPar        *parDem,
				  REAL8                threshold)
{
  REAL8UnitPolarCoor sourceLocation, sourceLocationBest;
  REAL8 deltaF, f0, fdot, dFdot, patchSizeX, patchSizeY;
  INT8 f0Bin;
  INT4 i,j, xSide, ySide, numCandidates;
  SemiCohCandidate thisCandidate;
  UINT2 criteria=1;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);


  ASSERT ( out != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( out->length > 0, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  ASSERT ( out->list != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

  ASSERT ( ht != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( patch != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( parDem != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( threshold > 0, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );

  deltaF = ht->deltaF;
  f0Bin = ht->f0Bin;
  f0 = f0Bin * deltaF;

  fdot = ht->spinDem.data[0] + ht->spinRes.data[0];
  dFdot = ht->dFdot.data[0];

  xSide = ht->xSide;
  ySide = ht->ySide;
  patchSizeX = ht->patchSizeX;
  patchSizeY = ht->patchSizeY;

  numCandidates = out->nCandidates;


  /* choose local max and threshold crossing */
  if (criteria == 0) {

    thisCandidate.freq =  f0;
    thisCandidate.dFreq = deltaF;
    thisCandidate.fdot = fdot;
    thisCandidate.dFdot = dFdot;
    thisCandidate.dAlpha = 3.0 * patchSizeX / ((REAL8)xSide);
    thisCandidate.dDelta = 3.0 * patchSizeY / ((REAL8)ySide);


    for (i = 1; i < ySide-1; i++)
      {
	for (j = 1; j < xSide-1; j++)
	  {

	    BOOLEAN isLocalMax = TRUE;

	    thisCandidate.significance =  ht->map[i*xSide + j];

	    /* check if this is a local maximum */
	    isLocalMax = (thisCandidate.significance > ht->map[i*xSide + j - 1])
	      && (thisCandidate.significance > ht->map[i*xSide + j + 1])
	      && (thisCandidate.significance > ht->map[(i+1)*xSide + j])
	      && (thisCandidate.significance > ht->map[(i-1)*xSide + j])
	      && (thisCandidate.significance > ht->map[(i-1)*xSide + j - 1])
	      && (thisCandidate.significance > ht->map[(i-1)*xSide + j + 1])
	      && (thisCandidate.significance > ht->map[(i+1)*xSide + j - 1])
	      && (thisCandidate.significance > ht->map[(i+1)*xSide + j + 1]);

	    /* realloc list if necessary */
	    if (numCandidates >= out->length) {
	      out->length += BLOCKSIZE_REALLOC;
	      out->list = LALRealloc( out->list, alloc_len = out->length * sizeof(SemiCohCandidate));
              if ( out->list == NULL ) {
                XLALPrintError ("Failed to LALRealloc ( 1, %d )\n", alloc_len );
                ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
              }
	      LogPrintf(LOG_DETAIL, "Need to realloc Hough candidate list to %d entries\n", out->length);
	    } /* need a safeguard to ensure that the reallocs don't happen too often */

	    /* add to list if candidate exceeds threshold and there is enough space in list */
	    if( ((REAL8)thisCandidate.significance > threshold) && (numCandidates < out->length) && isLocalMax) {
	      /* get sky location of pixel */

	      TRY( LALStereo2SkyLocation (status->statusPtr, &sourceLocation,
					  j, i, patch, parDem), status);

	      /* the above function uses atan2() which returns alpha in
		 the range [-pi,pi] => add pi to it */
	      thisCandidate.alpha = sourceLocation.alpha;
	      thisCandidate.delta = sourceLocation.delta;

	      out->list[numCandidates] = thisCandidate;
	      numCandidates++;
	      out->nCandidates = numCandidates;
	    }


	  } /* end loop over xSide */

      } /* end loop over ySide */

  } /* local maximum criteria for candidate selection */


  /* choose global maximum */
  if ( criteria == 1) {

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

    TRY( LALStereo2SkyLocation (status->statusPtr, &sourceLocationBest,
				jMax, iMax, patch, parDem), status);

    thisCandidate.alphaBest = sourceLocationBest.alpha;
    thisCandidate.deltaBest = sourceLocationBest.delta;

    TRY( LALHoughmapMeanVariance( status->statusPtr, &meanSig, &varianceSig, ht), status);
    thisCandidate.meanSig = meanSig;
    thisCandidate.varianceSig = varianceSig;

    /* realloc list if necessary */
    if (numCandidates >= out->length) {
      out->length += BLOCKSIZE_REALLOC;
      out->list = LALRealloc( out->list, alloc_len = out->length * sizeof(SemiCohCandidate));
      if ( out->list == NULL ) {
        XLALPrintError ("Failed to LALRealloc ( 1, %d )\n", alloc_len );
        ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
      }
      LogPrintf(LOG_DETAIL, "Need to realloc Hough candidate list to %d entries\n", out->length);
    } /* need a safeguard to ensure that the reallocs don't happen too often */


    /* get sky location of pixel */

    /*     TRY( LALStereo2SkyLocation (status->statusPtr, &sourceLocation,  */
    /* 				jMax, iMax, patch, parDem), status); */

   /* the above function uses atan2() which returns alpha in
       the range [-pi,pi] => add pi to it */
    /* thisCandidate.alpha = sourceLocation.alpha + LAL_PI; */
    /* thisCandidate.delta = sourceLocation.delta; */

    /* return coordinates of hough map center only */
    thisCandidate.alpha = parDem->skyPatch.alpha;
    thisCandidate.delta = parDem->skyPatch.delta;

    out->list[numCandidates] = thisCandidate;
    numCandidates++;
    out->nCandidates = numCandidates;
  }

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* end hough threshold selection */





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

  INIT_MEM ( fkdotIn );
  INIT_MEM ( fkdotOut );

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

  INIT_MEM(fkdot);

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



void GetFstatCandidates_toplist( LALStatus *status,
				 toplist_t *list,
				 REAL8FrequencySeries *in,
				 REAL8 alpha,
				 REAL8 delta,
				 REAL8 fdot)
{
  INT4 k, length;
  REAL8 deltaF, f0;
  HoughFStatOutputEntry line;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* fill up alpha, delta and fdot in fstatline */
  line.f1dot = fdot;
  line.Alpha = alpha;
  line.Delta = delta;

  length = in->data->length;
  deltaF = in->deltaF;
  f0 = in->f0;

  /* loop over Fstat vector and fill up toplist */
  for ( k=0; k<length; k++)
    {

      line.HoughFStat = in->data->data[k];
      line.Freq = f0 + k*deltaF;

      insert_into_houghFStat_toplist( list, line);

    }

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

/** Calculate average velocity and position of detector network during each
    stack */
void GetStackVelPos( LALStatus *status,
		     REAL8VectorSequence **velStack,
		     REAL8VectorSequence **posStack,
		     MultiDetectorStateSeriesSequence *stackMultiDetStates)
{
  UINT4 k, j, m, nStacks;
  INT4 counter, numifo;
  CreateVectorSequenceIn createPar;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);


  ASSERT ( stackMultiDetStates != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( stackMultiDetStates->length > 0, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  ASSERT ( stackMultiDetStates->data != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( velStack != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( posStack != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

  nStacks = stackMultiDetStates->length;

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
      counter=0;
      numifo = stackMultiDetStates->data[k]->length;

      /* initialize velocities and positions */
      velStack[0]->data[3*k] = 0;
      velStack[0]->data[3*k+1] = 0;
      velStack[0]->data[3*k+2] = 0;

      posStack[0]->data[3*k] = 0;
      posStack[0]->data[3*k+1] = 0;
      posStack[0]->data[3*k+2] = 0;

      for ( j = 0; (INT4)j < numifo; j++)
	{
	  INT4 numsft = stackMultiDetStates->data[k]->data[j]->length;
	  for ( m = 0; (INT4)m < numsft; m++)
	    {
	      /* sum velocity components */
	      velStack[0]->data[3*k] += stackMultiDetStates->data[k]->data[j]->data[m].vDetector[0];
	      velStack[0]->data[3*k+1] += stackMultiDetStates->data[k]->data[j]->data[m].vDetector[1];
	      velStack[0]->data[3*k+2] += stackMultiDetStates->data[k]->data[j]->data[m].vDetector[2];
	      /* sum position components */
	      posStack[0]->data[3*k] += stackMultiDetStates->data[k]->data[j]->data[m].rDetector[0];
	      posStack[0]->data[3*k+1] += stackMultiDetStates->data[k]->data[j]->data[m].rDetector[1];
	      posStack[0]->data[3*k+2] += stackMultiDetStates->data[k]->data[j]->data[m].rDetector[2];

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
			       MultiNoiseWeightsSequence *in )
{

  INT4 nStacks, k, j, i, numifo, numsft;
  REAL8Vector *weightsVec=NULL;
  MultiNoiseWeights *multNoiseWts;


  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( in != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( in->length > 0, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  ASSERT ( in->data != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

  ASSERT ( out != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( *out == NULL, status, HIERARCHICALSEARCH_ENONULL, HIERARCHICALSEARCH_MSGENONULL );

  nStacks = in->length;

  TRY( LALDCreateVector( status->statusPtr, &weightsVec, nStacks), status);


  for (k=0; k<nStacks; k++) {

    multNoiseWts = in->data[k];
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
				    MultiNoiseWeightsSequence *inNoise,
				    MultiDetectorStateSeriesSequence *inDetStates,
				    SkyPosition skypos)
{

  UINT4 nStacks, iStack, iIFO, iSFT, numifo, numsft;
  REAL8 a, b, n;
  MultiNoiseWeights *multNoiseWts;
  MultiDetectorStateSeries *multDetStates;
  MultiAMCoeffs *multiAMcoef = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( inNoise != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( inNoise->length > 0, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  ASSERT ( inNoise->data != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

  nStacks = inNoise->length;

  ASSERT ( inDetStates != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( inDetStates->length == nStacks, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  ASSERT ( inDetStates->data != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

  ASSERT ( out != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  ASSERT ( out->length == nStacks, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  ASSERT ( out->data != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );


  for (iStack=0; iStack<nStacks; iStack++) {

    multNoiseWts = inNoise->data[iStack];
    ASSERT ( multNoiseWts != NULL, status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );

    multDetStates = inDetStates->data[iStack];
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

  INIT_MEM ( line );

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
  patch.xCoor = LALCalloc(1, alloc_len = patch.xSide*sizeof(REAL8));
  if ( patch.xCoor == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  patch.yCoor = LALCalloc(1, alloc_len = patch.ySide*sizeof(REAL8));
  if ( patch.yCoor == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  TRY( LALHOUGHFillPatchGrid( status->statusPtr, &patch, &parSize ), status );

  /* the demodulation params */
  parDem.deltaF = deltaF;
  parDem.skyPatch.alpha = alpha;
  parDem.skyPatch.delta = delta;
  parDem.spin.length = 1;
  parDem.spin.data = LALCalloc(1, alloc_len = sizeof(REAL8));
  if ( parDem.spin.data == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  parDem.spin.data[0] = fdot;

  TRY( LALDCreateVector( status->statusPtr, &timeDiffV, nStacks), status);

  for (j=0; j<nStacks; j++) {
    tMidStack = XLALGPSGetREAL8(tsMid->data + j);
    timeDiffV->data[j] = tMidStack - refTime;
  }


  lut.maxNBins = parSize.maxNBins;
  lut.maxNBorders = parSize.maxNBorders;
  lut.border = LALCalloc(1, alloc_len = parSize.maxNBorders*sizeof(HOUGHBorder));
  if ( lut.border == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  lut.bin = LALCalloc(1, alloc_len = parSize.maxNBins*sizeof(HOUGHBin2Border));
  if ( lut.bin == NULL ) {
    XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
    ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
  }

  for (i = 0; i < parSize.maxNBorders; i++){
    lut.border[i].ySide = parSize.ySide;
    lut.border[i].xPixel = LALCalloc(1, alloc_len = parSize.ySide*sizeof(COORType));
    if ( lut.border[i].xPixel == NULL ) {
      XLALPrintError ("Failed to LALCalloc ( 1, %d )\n", alloc_len );
      ABORT ( status, HIERARCHICALSEARCH_EMEM, HIERARCHICALSEARCH_MSGEMEM );
    }
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
