/*
 *  Copyright (C) 2005 Badri Krishnan, Alicia Sintes  
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
 * \author Badri Krishnan, Alicia Sintes
 * \file HierarchicalSearch.c
 * \brief Program for calculating F-stat values for different time segments 
   and combining them semi-coherently using the Hough transform, and following 
   up candidates using a more sensitive search.

   \par  Description
   
   This code implements a hierarchical strategy to look for unknown gravitational 
   wave pulsars. It scans through the parameter space using a less sensitive but 
   computationally inexpensive search and follows up the candidates using 
   more sensitive methods.  

   \par Algorithm
   
   Currently the code does a single stage hierarchical search using the Hough
   algorithm and follows up the candidates using a full coherent integration.  

   - The user specifies a directory containing SFTs, and the number of \e stacks 
     that this must be broken up into.  
     At present two ways of breaking up the SFTs into stacks are supported. 
     These are equivalent if there are no gaps in the data. Either the SFTs are 
     divided up equally among the stacks, or the total time spanned by the data is 
     broken up equally

   - The user specifies a region in parameter space to search over.  The code 
     sets up a grid (the "coarse" grid) in this region and calculates the F-statistic 
     for each stack at each point of the coarse grid.  

   - A threshold is set on the F-statistic to convert the 
     F-statistic vector into a vector of 0s and 1s known as a \e peakgram -- there 
     is one peakgram for each stack and for each grid point.

   - The peakgrams are combined using the Hough transform.  

   - The Hough part of the search constructs a grid (the "fine" grid) in a small patch 
     around every coarse grid point, and combines the different stacks
     following the \e master equation 
     \f[
        f(t) - F_0(t) = \xi(t).(\hat{n} - \hat{n}_0)
     \f]	
     where 
     \f[
        F_0 = f_0 + \sum \Delta  f_k {(\Delta t)^k \over k!}
     \f]
     Here \f$ \hat{n}_0 \f$ is the sky-point at which the F-statistic is calculated
     and \f$ \Delta f_k \f$ is the \e residual spindown parameter.  For details see
     Phys.Rev.D 70, 082001 (2004).  The size of the patch depends on the validity of
     the above master equation.  

   - The output of the Hough search is a \e number \e count at point of the grid. 
     A threshold is set on the number count, leading to candidates in parameter space.

   - These candidates are followed up using a second set of SFTs (also specified by
     the user).  The follow up consists of a full coherent integration, i.e. the F-statistic
     is calculated for the whole set of SFTs without breaking them up into stacks. 
     A threshold can be set on the F-statistic to get the final list of candidates.  


   \par Immediate to-do list

   - The reference time is not yet handled correctly -- ok if it the default, i.e. the start 
     time of the first SFT but not in general

   - Use average velocity instead of mid-time

   \par Longer term

   - Should we over-resolve the Fstat calculation to reduce loss in signal power?  We would
     still calculate the peakgrams at the 1/T resolution, but the peak selection would 
     take into account Fstat values over several over-resolved bins.  
   
   - What is the best grid for calculating the F-statistic?  At first glance, the 
     Hough patch and the metric F-statistic patch do not seem to be compatible.  If we
     use the Hough patches to break up the sky, it could be far from optimal.  What is 
     the exact connection between the two grids?  

   - Implement multiple semi-coherent stages

   - Get timings and optimize the pipeline parameters

   - Checkpointing for running on Einstein@Home

   - Incorporate stack slide as an alternative to Hough in the semi-coherent stages

   - ....

 */



#include"./HierarchicalSearch.h"

RCSID( "$Id$");

#define TRUE (1==1)
#define FALSE (1==0)

extern int lalDebugLevel;

BOOLEAN uvar_printMaps; /**< global variable for printing Hough maps */
BOOLEAN uvar_printStats; /**< global variable for calculating Hough map stats */


#define HSMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define HSMIN(x,y) ( (x) < (y) ? (x) : (y) )


typedef struct {
  CHAR  *sftbasename;
  REAL8 deltaF;          /**< SFT frequency resolution */
  REAL8 deltaFstack;     /**< frequency resolution of stacks */
  UINT4 nStacks;         /**< number of stacks */
  LIGOTimeGPS tStartGPS; /**< start and end time of stack */
  REAL8 tObs;            /**< tEndGPS - tStartGPS */
  REAL8 refTime;
  LALPulsarSpinRange *spinRange_startTime;
  LALPulsarSpinRange *spinRange_endTime;
  LALPulsarSpinRange *spinRange_refTime;
  EphemerisData *edat;
  LIGOTimeGPSVector *midTstack; 
  LIGOTimeGPSVector *startTstack;   
  LIGOTimeGPS minStartTimeGPS;
  LIGOTimeGPS maxEndTimeGPS;
  UINT4 blocksRngMed;
  UINT4 Dterms;
  REAL8 dopplerMax;
} UsefulStageVariables;


/* functions for printing various stuff */
void GetStackVelPos( LALStatus *status, REAL8VectorSequence **velStack, REAL8VectorSequence **posStack,
		     MultiDetectorStateSeriesSequence *stackMultiDetStates);


void SetUpSFTs( LALStatus *status, MultiSFTVectorSequence *stackMultiSFT, MultiNoiseWeightsSequence *stackMultiNoiseWeights,
		MultiDetectorStateSeriesSequence *stackMultiDetStates, UsefulStageVariables *in);

void PrintFstatVec (LALStatus *status, REAL8FrequencySeries *in, FILE *fp, PulsarDopplerParams *thisPoint, 
		    LIGOTimeGPS  refTime, INT4 stackIndex);

void PrintSemiCohCandidates(LALStatus *status, SemiCohCandidateList *in, FILE *fp);

void PrintHoughHistogram(LALStatus *status, UINT8Vector *hist, CHAR *fnameOut);

void PrintCatalogInfo( LALStatus  *status, const SFTCatalog *catalog, FILE *fp);

void PrintStackInfo( LALStatus  *status, const SFTCatalogSequence *catalogSeq, FILE *fp);

/* default values for input variables */
#define EPHEMERISDIR "/local_data/badkri/lscsoft/share/lal/" /**< Default location of ephemeris files*/
#define EPHEMERISYEAR "05-09"   /**< Default location of sun ephemeris */
#define NSTACKS 10    /**< Default number of stacks */
#define BLOCKSRNGMED 101 /**< Default running median window size */
#define FSTART 310.0   /**< Default Start search frequency */

#define FBAND 0.01    /**< Default search band */
#define FDOT 0.0      /**< Default value of first spindown */
#define DFDOT 0.0   /**< Default range of first spindown parameter */
#define SKYREGION "(1,1),(1,1.3),(1.3,1.3),(1.3,1)" /**< default sky region to search over -- just a single point*/
#define NFDOT  10    /**< Default size of hough cylinder of look up tables */
#define DTERMS 8     /**< Default number of dirichlet kernel terms for calculating Fstat */
#define MISMATCH 0.2 /**< Default for metric grid maximal mismatch value */
#define DALPHA 0.001 /**< Default resolution for isotropic or flat grids */
#define DDELTA 0.001 /**< Default resolution for isotropic or flat grids */
#define FSTATTHRESHOLD 2.6  /**< Default threshold on Fstatistic for peak selection */
#define NCAND1 5 /**< Default number of candidates to be followed up from first stage */
#define SFTDIRECTORY "/home/badkri/fakesfts2/H-1_H1*.sft"  /**< Default directory containing sfts */
#define FNAMEOUT "./out/HS"  /**< Default output file basename */
#define PIXELFACTOR 2
#define LAL_INT4_MAX 2147483647

int main( int argc, char *argv[]) {

  /* initialize status */
  LALStatus status = blank_status;
  
  /* temp loop variables: generally k loops over stacks and j over SFTs in a stack*/
  INT4 j, k;

  /* in general any variable ending with 1 is for the 
     first stage, 2 for the second and so on */

  /* ephemeris */
  EphemerisData *edat = NULL;

  /* timestamp vectors */
  LIGOTimeGPSVector *midTstack1=NULL; 
  LIGOTimeGPSVector *midTstack2=NULL; 
  LIGOTimeGPSVector *startTstack1=NULL; 
  LIGOTimeGPSVector *startTstack2=NULL; 
  
  LIGOTimeGPSVector *sftTimeStamps1=NULL; 
  LIGOTimeGPSVector *sftTimeStamps2=NULL; 
  LIGOTimeGPS refTimeGPS, tStart1GPS, tStart2GPS;
  REAL8 tObs1, tObs2;

  /* velocities and positions at midTstack */
  REAL8VectorSequence *velStack1=NULL;
  REAL8VectorSequence *posStack1=NULL;

  /* duration of each stack */
  REAL8 tStack1, tStack2;

  /* sft related stuff */
  static MultiSFTVectorSequence stackMultiSFT1;
  static MultiSFTVectorSequence stackMultiSFT2;
  static MultiNoiseWeightsSequence stackMultiNoiseWeights1;
  static MultiNoiseWeightsSequence stackMultiNoiseWeights2;
  static MultiDetectorStateSeriesSequence stackMultiDetStates1;
  static MultiDetectorStateSeriesSequence stackMultiDetStates2;
  LIGOTimeGPS minStartTimeGPS1, maxEndTimeGPS1;
  LIGOTimeGPS minStartTimeGPS2, maxEndTimeGPS2;

  /* some useful variables for each stage */
  UsefulStageVariables usefulParams1, usefulParams2;

  /* number of stacks -- not necessarily same as uvar_nStacks! */
  INT4 nStacks1, nStacks2;
  REAL8 deltaF1, deltaF2, timebase1, timebase2; /* frequency resolution of SFTs */
  REAL8 deltaFstack1, deltaFstack2; /* frequency resolution of Fstat calculation */

  /* LALdemod related stuff */
  REAL8FrequencySeriesVector fstatVector1, fstatVector2; /* Fstatistic vectors for each stack */
  UINT4 binsFstat1, binsFstat2;
  ComputeFParams CFparams;		   

  /* hough variables */
  HOUGHPeakGramVector pgV;
  SemiCoherentParams semiCohPar;
  SemiCohCandidateList semiCohCandList1;

  /* fstat candidate structure */
  toplist_t *fstatToplist=NULL;

  /* template and grid variables */
  DopplerScanInit scanInit1, scanInit2;   /* init-structure for DopperScanner */
  DopplerScanState thisScan1 = empty_DopplerScanState; /* current state of the Doppler-scan */
  DopplerScanState thisScan2 = empty_DopplerScanState; /* current state of the Doppler-scan */
  DopplerPosition dopplerpos1, dopplerpos2;	       /* current search-parameters */
  PulsarDopplerParams thisPoint1, thisPoint2; 

  /* various spinrange vectors */
  LALPulsarSpinRange *spinRange_refTime, *spinRange_startTime1, *spinRange_startTime2;
  LALPulsarSpinRange *spinRange_endTime1, *spinRange_endTime2;
  LALPulsarSpinRange *spinRange_Temp=NULL;

  /* variables for logging */
  CHAR *fnamelog=NULL;
  FILE *fpLog=NULL;
  CHAR *logstr=NULL; 

  /* output candidate files and file pointers */
  CHAR *fnameFstatCand=NULL;
  CHAR *fnameSemiCohCand=NULL;
  CHAR *fnameFstatVec1=NULL;
  FILE *fpFstat=NULL;
  FILE *fpSemiCoh=NULL;
  FILE *fpFstat1=NULL;
  
  /* checkpoint filename and index of loop over skypoints */
  /*   CHAR *fnameChkPoint=NULL; */
  /*   FILE *fpChkPoint=NULL; */
  /*   UINT4 loopindex, loopcounter; */
  
  /* user variables */
  BOOLEAN uvar_help; /* true if -h option is given */
  BOOLEAN uvar_log; /* logging done if true */
  BOOLEAN uvar_printCand1; /* if 1st stage candidates are to be printed */
  BOOLEAN uvar_chkPoint;
  BOOLEAN uvar_followUp;
  BOOLEAN uvar_printFstat1;

  REAL8 uvar_dAlpha, uvar_dDelta; /* resolution for flat or isotropic grids */
  REAL8 uvar_f1dot; /* first spindown value */
  REAL8 uvar_f1dotBand; /* range of first spindown parameter */
  REAL8 uvar_Freq, uvar_FreqBand;
  REAL8 uvar_peakThrF; /* threshold of Fstat to select peaks */
  REAL8 uvar_houghThr; /* threshold on hough number count to select candidates */
  REAL8 uvar_fstatThr; /* threshold for selecting candidates from Fstat vector */
  REAL8 uvar_mismatch1; /* metric mismatch for first stage coarse grid */
  REAL8 uvar_mismatch2; /* metric mismatch for second stage coarse grid */
  REAL8 uvar_refTime;
  REAL8 uvar_minStartTime1, uvar_maxEndTime1;
  REAL8 uvar_minStartTime2, uvar_maxEndTime2;
  REAL8 uvar_dopplerMax;
  REAL8 uvar_pixelFactor;

  INT4 uvar_method; /* hough = 0, stackslide = 1*/
  INT4 uvar_nCand1; /* number of candidates to be followed up from first stage */
  INT4 uvar_nCand2; /* number of candidates from second stage */
  INT4 uvar_blocksRngMed;
  INT4 uvar_nStacks1, uvar_nStacks2;
  INT4 uvar_Dterms;
  INT4 uvar_SSBprecision;
  INT4 uvar_nfdot;
  INT4 uvar_gridType;
  INT4 uvar_metricType;
  INT4 uvar_reallocBlock;

  CHAR *uvar_ephemDir=NULL;
  CHAR *uvar_ephemYear=NULL;
  CHAR *uvar_sftData1=NULL;
  CHAR *uvar_sftData2=NULL;
  CHAR *uvar_fnameout=NULL;
  CHAR *uvar_skyGridFile=NULL;
  CHAR *uvar_skyRegion=NULL;


  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /*---------------------------------------------------------------*/
  /* set defaults, read user variables, log user variables and log cvs tags */
  /* LALDebugLevel must be called before anything else */
  lalDebugLevel = 0;
  LAL_CALL( LALGetDebugLevel( &status, argc, argv, 'd'), &status);

  /* now set the other defaults */
  uvar_help = FALSE;
  uvar_log = FALSE;
  uvar_method = 0;
  uvar_followUp = TRUE;
  uvar_printMaps = FALSE;
  uvar_printStats = FALSE;
  uvar_printCand1 = FALSE;
  uvar_printFstat1 = FALSE;
  uvar_chkPoint = FALSE;
  uvar_nStacks1 = NSTACKS;
  uvar_nStacks2 = 1;
  uvar_Dterms = DTERMS;
  uvar_dAlpha = DALPHA;
  uvar_dDelta = DDELTA;
  uvar_f1dot = FDOT;
  uvar_f1dotBand = DFDOT;
  uvar_Freq = FSTART;
  uvar_FreqBand = FBAND;
  uvar_blocksRngMed = BLOCKSRNGMED;
  uvar_nfdot = NFDOT;
  uvar_peakThrF = FSTATTHRESHOLD;
  uvar_nCand1 = uvar_nCand2 = NCAND1;
  uvar_houghThr = 0;
  uvar_fstatThr = FSTATTHRESHOLD;
  uvar_SSBprecision = SSBPREC_RELATIVISTIC;
  uvar_metricType = LAL_PMETRIC_COH_PTOLE_ANALYTIC;
  uvar_gridType = GRID_METRIC;
  uvar_mismatch1 = uvar_mismatch2 = MISMATCH;
  uvar_reallocBlock = 5000;
  uvar_pixelFactor = PIXELFACTOR;

  uvar_minStartTime1 = 0;
  uvar_maxEndTime1 = LAL_INT4_MAX;

  uvar_minStartTime2 = 0;
  uvar_maxEndTime2 = LAL_INT4_MAX;

  uvar_dopplerMax = 1.05e-4;

  uvar_skyGridFile = NULL;

  uvar_skyRegion = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_skyRegion, SKYREGION);

  uvar_ephemDir = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_ephemDir, EPHEMERISDIR);

  uvar_ephemYear = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_ephemYear, EPHEMERISYEAR);

  uvar_sftData1 = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_sftData1, SFTDIRECTORY);

  /* do not set default for sftData2 -- use only if user specifies */
  /*   uvar_sftData2 = (CHAR *)LALMalloc(512*sizeof(CHAR)); */

  uvar_fnameout = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_fnameout, FNAMEOUT);

  /* register user input variables */
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "help",        'h', UVAR_HELP,     "Print this message", &uvar_help), &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "log",          0,  UVAR_OPTIONAL, "Write log file", &uvar_log), &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "chkPoint",     0,  UVAR_OPTIONAL, "For checkpointing", &uvar_chkPoint), &status);  
  LAL_CALL( LALRegisterINTUserVar(    &status, "method",       0,  UVAR_OPTIONAL, "Hough=0, stackslide=1", &uvar_method ), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "followUp",     0,  UVAR_OPTIONAL, "Follow up stage?", &uvar_followUp), &status);  
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftData1",     0,  UVAR_OPTIONAL, "1st SFT file pattern", &uvar_sftData1), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftData2",     0,  UVAR_OPTIONAL, "2nd SFT file pattern", &uvar_sftData2), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "skyRegion",    0,  UVAR_OPTIONAL, "sky-region polygon (or 'allsky')", &uvar_skyRegion), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "Freq",        'f', UVAR_OPTIONAL, "Start search frequency", &uvar_Freq), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "FreqBand",    'b', UVAR_OPTIONAL, "Search frequency band", &uvar_FreqBand), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "f1dot",        0,  UVAR_OPTIONAL, "Spindown parameter", &uvar_f1dot), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "f1dotBand",    0,  UVAR_OPTIONAL, "Range of fdot", &uvar_f1dotBand), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "nStacks1",    'N', UVAR_OPTIONAL, "No.of 1st stage stacks", &uvar_nStacks1 ), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "nStacks2",     0,  UVAR_OPTIONAL, "No.of 2nd stage stacks", &uvar_nStacks2 ), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "mismatch1",    0,  UVAR_OPTIONAL, "1st stage mismatch", &uvar_mismatch1), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "mismatch2",    0,  UVAR_OPTIONAL, "2nd stage mismatch", &uvar_mismatch2), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "pixelFactor",  0,  UVAR_OPTIONAL, "Semi coherent sky grid resolution", &uvar_pixelFactor), &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "gridType",     0,  UVAR_OPTIONAL, "0=flat,1=isotropic,2=metric,3=file", &uvar_gridType),  &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "metricType",   0,  UVAR_OPTIONAL, "0=none,1=Ptole-analytic,2=Ptole-numeric,3=exact", &uvar_metricType), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "skyGridFile",  0,  UVAR_OPTIONAL, "sky-grid file", &uvar_skyGridFile), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dAlpha",       0,  UVAR_OPTIONAL, "Resolution for flat or isotropic grids", &uvar_dAlpha), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dDelta",       0,  UVAR_OPTIONAL, "Resolution for flat or isotropic grids", &uvar_dDelta), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "fnameout",    'o', UVAR_OPTIONAL, "Output basefileneme", &uvar_fnameout), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "peakThrF",     0,  UVAR_OPTIONAL, "Fstat Threshold", &uvar_peakThrF), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "nCand1",       0,  UVAR_OPTIONAL, "No.of 1st stage candidates to be followed up", &uvar_nCand1), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "nCand2",       0,  UVAR_OPTIONAL, "No.of 2nd stage candidates to be followed up",&uvar_nCand2), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "houghThr",     0,  UVAR_OPTIONAL, "Hough number count threshold (default --nCand1)",   &uvar_houghThr), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printCand1",   0,  UVAR_OPTIONAL, "Print 1st stage candidates", &uvar_printCand1), &status);  
  LAL_CALL( LALRegisterREALUserVar(   &status, "refTime",      0,  UVAR_OPTIONAL, "Ref. time for pulsar pars [start time]", &uvar_refTime), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "ephemDir",     0,  UVAR_OPTIONAL, "Location of ephemeris files", &uvar_ephemDir),  &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "ephemYear",    0,  UVAR_OPTIONAL, "Ephemeris year", &uvar_ephemYear),  &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "minStartTime1",0,  UVAR_OPTIONAL, "1st Min start time of observation", &uvar_minStartTime1), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "maxEndTime1",  0,  UVAR_OPTIONAL, "1st Max end time of observation",   &uvar_maxEndTime1),   &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "minStartTime2",0,  UVAR_OPTIONAL, "2nd Min start time of observation", &uvar_minStartTime2), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "maxEndTime2",   0,  UVAR_OPTIONAL, "2nd Max end time of observation",  &uvar_maxEndTime2),   &status);

  /* developer user variables */
  LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed", 0, UVAR_DEVELOPER, "RngMed block size", &uvar_blocksRngMed), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fstatThr",     0, UVAR_DEVELOPER, "Fstat threshold (default --nCand2)", &uvar_fstatThr), &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "reallocBlock", 0, UVAR_DEVELOPER,"Blocks to realloc for Fstat output if necessary",   &uvar_reallocBlock),    &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "SSBprecision", 0, UVAR_DEVELOPER,"Precision for SSB transform.", &uvar_SSBprecision),    &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "nfdot",        0, UVAR_DEVELOPER,"No.of residual fdot values to be searched", &uvar_nfdot), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printMaps",    0, UVAR_DEVELOPER,"Print Hough maps", &uvar_printMaps), &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printStats",   0, UVAR_DEVELOPER,"Print Hough map statistics", &uvar_printStats), &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printFstat1",  0, UVAR_DEVELOPER,"Print 1st stage Fstat vectors", &uvar_printFstat1), &status);  
  LAL_CALL( LALRegisterINTUserVar(    &status, "Dterms",       0, UVAR_DEVELOPER,"No.of terms to keep in Dirichlet Kernel", &uvar_Dterms ), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dopplerMax",   0, UVAR_DEVELOPER,"Max Doppler shift",  &uvar_dopplerMax), &status);


  /* read all command line variables */
  LAL_CALL( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 

  /* some basic sanity checks on user vars */
  if ( (uvar_method != 0) && (uvar_method != 1)) {
    fprintf(stderr, "Invalid method....must be 0 or 1\n");
    exit( HIERARCHICALSEARCH_EBAD );
  }

  if ( uvar_nStacks1 < 1) {
    fprintf(stderr, "Invalid number of stacks\n");
    exit( HIERARCHICALSEARCH_EBAD );
  }

  if ( uvar_nStacks2 < 1) {
    fprintf(stderr, "Invalid number of stacks\n");
    exit( HIERARCHICALSEARCH_EBAD );
  }

  if ( uvar_blocksRngMed < 1 ) {
    fprintf(stderr, "Invalid Running Median block size\n");
    exit( HIERARCHICALSEARCH_EBAD );
  }

  if ( uvar_peakThrF < 0 ) {
    fprintf(stderr, "Invalid value of Fstatistic threshold\n");
    exit( HIERARCHICALSEARCH_EBAD );
  }

  if ( uvar_followUp && (!LALUserVarWasSet(&uvar_sftData2))) {
    fprintf( stderr, "Must specify SFTs for second stage!\n");
    exit( HIERARCHICALSEARCH_EBAD );
  }

  /* no need to follow up zero 1st stage candidates */
  if ( uvar_nCand1 <= 0 )
    uvar_followUp = FALSE;

  /* no need to follow up if no candidates are to be saved */
  if ( uvar_nCand2 <= 0 )
    uvar_followUp = FALSE;


  /* write the log file */
  if ( uvar_log ) 
    {
      fnamelog = (CHAR *)LALMalloc( 512*sizeof(CHAR));
      strcpy(fnamelog, uvar_fnameout);
      strcat(fnamelog, ".log");
      /* open the log file for writing */
      if ((fpLog = fopen(fnamelog, "w")) == NULL) {
	fprintf(stderr, "Unable to open file %s for writing\n", fnamelog);
	LALFree(fnamelog);
	exit(1);
      }

      /* get the log string */
      LAL_CALL( LALUserVarGetLog(&status, &logstr, UVAR_LOGFMT_CFGFILE), &status);  
      
      fprintf( fpLog, "## Log file for DriveHoughFStat\n\n");
      fprintf( fpLog, "# User Input:\n");
      fprintf( fpLog, "#-------------------------------------------\n");
      fprintf( fpLog, logstr);
      LALFree(logstr);
      
      /*get the cvs tags */
      {
	CHAR command[1024] = "";
	fprintf (fpLog, "\n\n# CVS-versions of executable:\n");
	fprintf (fpLog, "# -----------------------------------------\n");
	fclose (fpLog);
	
	sprintf (command, "ident %s | sort -u >> %s", argv[0], fnamelog);
	system (command);	/* we don't check this. If it fails, we assume that */
	/* one of the system-commands was not available, and */
	/* therefore the CVS-versions will not be logged */
	
	LALFree(fnamelog); 
	
      } /* end of cvs tags block */
      
    } /* end of logging */
  

  /*--------- Some initializations ----------*/

  /* initialize ephemeris info */ 
  {
    CHAR ephemE[512], ephemS[512];
    
    sprintf(ephemE, "%s/earth%s.dat", uvar_ephemDir, uvar_ephemYear);
    sprintf(ephemS, "%s/sun%s.dat", uvar_ephemDir, uvar_ephemYear);
    
    edat = (EphemerisData *)LALMalloc(sizeof(EphemerisData));
    (*edat).ephiles.earthEphemeris = ephemE;
    (*edat).ephiles.sunEphemeris = ephemS;
    
    /* read in ephemeris data */
    LAL_CALL( LALInitBarycenter( &status, edat), &status);        
  }
  
  /* create spin range vectors */
  spinRange_refTime = XLALCreatePulsarSpinRange(2);
  spinRange_startTime1 = XLALCreatePulsarSpinRange(2);
  spinRange_startTime2 = XLALCreatePulsarSpinRange(2);
  spinRange_endTime1 = XLALCreatePulsarSpinRange(2);
  spinRange_endTime2 = XLALCreatePulsarSpinRange(2);

  /* utility temp variable */
  spinRange_Temp = XLALCreatePulsarSpinRange(2);

  /* copy user specified spin variables at reftime  */
  for (j = 0; j < 2; j++)  
    {
      spinRange_refTime->fkdot->data[0] = uvar_Freq; /* frequency */
      spinRange_refTime->fkdot->data[1] = uvar_f1dot;  /* 1st spindown */
      spinRange_refTime->fkdotBand->data[0] = uvar_FreqBand; /* frequency range */
      spinRange_refTime->fkdotBand->data[1] = uvar_f1dotBand; /* spindown range */
    }

  /* some compute F params */
  CFparams.Dterms = uvar_Dterms;
  CFparams.SSBprec = uvar_SSBprecision;

  LAL_CALL ( LALFloatToGPS( &status, &minStartTimeGPS1, &uvar_minStartTime1), &status);
  LAL_CALL ( LALFloatToGPS( &status, &maxEndTimeGPS1, &uvar_maxEndTime1), &status);
  LAL_CALL ( LALFloatToGPS( &status, &minStartTimeGPS2, &uvar_minStartTime2), &status);
  LAL_CALL ( LALFloatToGPS( &status, &maxEndTimeGPS2, &uvar_maxEndTime2), &status);
    
  /* create output Hough file for writing if requested by user */
  if ( uvar_printCand1 )
    {
      fnameSemiCohCand = (CHAR *)LALMalloc( 512*sizeof(CHAR));
      strcpy(fnameSemiCohCand, uvar_fnameout);
      strcat(fnameSemiCohCand, "_semicoh.txt");
      if (!(fpSemiCoh = fopen(fnameSemiCohCand, "w"))) 
	{
	  fprintf ( stderr, "Unable to open file '%s' for writing.\n", fnameSemiCohCand);
	  return HIERARCHICALSEARCH_EFILE;
	}
    }
  
  /* open output fstat file for writing */
  /* create output Fstat file name and open it for writing */
  if ( uvar_followUp ) 
    {
      
      fnameFstatCand = (CHAR *)LALMalloc( 512*sizeof(CHAR));
      strcpy(fnameFstatCand, uvar_fnameout);
      strcat(fnameFstatCand, "_fstat.txt");

      if  (!(fpFstat = fopen(fnameFstatCand, "w")))
	{
	  fprintf ( stderr, "Unable to open Fstat file '%s' for writing.\n", fnameFstatCand);
	  return HIERARCHICALSEARCH_EFILE;
	}
    } 

  if ( uvar_printFstat1 )
    {

      fnameFstatVec1 = (CHAR *)LALMalloc( 512*sizeof(CHAR));
      strcpy(fnameFstatVec1, uvar_fnameout);
      strcat(fnameFstatVec1, "_fstatVec1.txt");
      if ( !(fpFstat1 = fopen( fnameFstatVec1, "w")))
	{
	  fprintf ( stderr, "Unable to open Fstat file fstatvec1.out for writing.\n");
	  return HIERARCHICALSEARCH_EFILE;
	}
    }


  /*------------ Set up stacks, noise weights, detector states etc. */


  /* some useful first stage params */
  usefulParams1.sftbasename = uvar_sftData1;
  usefulParams1.nStacks = uvar_nStacks1;
  usefulParams1.spinRange_startTime = spinRange_startTime1;
  usefulParams1.spinRange_endTime = spinRange_endTime1;
  usefulParams1.spinRange_refTime = spinRange_refTime;
  usefulParams1.edat = edat;
  usefulParams1.minStartTimeGPS = minStartTimeGPS1;
  usefulParams1.maxEndTimeGPS = maxEndTimeGPS1;
  usefulParams1.blocksRngMed = uvar_blocksRngMed;
  usefulParams1.Dterms = uvar_Dterms;
  usefulParams1.dopplerMax = uvar_dopplerMax;

  /* set reference time for pular parameters */
  if ( LALUserVarWasSet(&uvar_refTime)) 
    usefulParams1.refTime = uvar_refTime;
  else 
    usefulParams1.refTime = -1;
  
  LAL_CALL( SetUpSFTs( &status, &stackMultiSFT1, &stackMultiNoiseWeights1, 
		       &stackMultiDetStates1, &usefulParams1), &status);


  /* some useful params computed by SetUpSFTs */
  deltaF1 = usefulParams1.deltaF;
  timebase1 = 1.0/deltaF1;
  deltaFstack1 = usefulParams1.deltaFstack;
  tStack1 = 1.0/deltaFstack1;
  tObs1 = usefulParams1.tObs;
  nStacks1 = usefulParams1.nStacks;
  tStart1GPS = usefulParams1.tStartGPS;
  midTstack1 = usefulParams1.midTstack;
  startTstack1 = usefulParams1.startTstack;
  refTimeGPS = spinRange_refTime->epoch;

  /* set reference time for calculating Fstatistic */
  thisPoint1.refTime = tStart1GPS;
  /* binary orbit and higher spindowns not considered so far */
  thisPoint1.orbit = NULL;
  thisPoint1.fkdot[2] = 0.0;
  thisPoint1.fkdot[3] = 0.0;
      
  /* allocate fstat memory */
  fstatVector1.length = nStacks1;
  fstatVector1.data = NULL;
  fstatVector1.data = (REAL8FrequencySeries *)LALCalloc( 1, nStacks1 * sizeof(REAL8FrequencySeries));
  binsFstat1 = (UINT4)(spinRange_startTime1->fkdotBand->data[0]/deltaFstack1 + 0.5) + 1;
  for (k = 0; k < nStacks1; k++) 
    { 
      /* careful--the epoch here is not the reference time for f0! */
      fstatVector1.data[k].epoch = startTstack1->data[k];
      fstatVector1.data[k].deltaF = deltaFstack1;
      fstatVector1.data[k].f0 = spinRange_startTime1->fkdot->data[0];
      fstatVector1.data[k].data = (REAL8Sequence *)LALCalloc( 1, sizeof(REAL8Sequence));
      fstatVector1.data[k].data->length = binsFstat1;
      fstatVector1.data[k].data->data = (REAL8 *)LALCalloc( 1, binsFstat1 * sizeof(REAL8));
    } 

  
  /*------------ calculate velocity and position for each 1st stage stack ------------*/

  LAL_CALL( GetStackVelPos( &status, &velStack1, &posStack1, &stackMultiDetStates1), &status);


  /*------------------ read sfts and set up stacks for follow up stage -----------------------*/
  /* check if user requested a follow up stage*/
  if ( uvar_followUp ) {
    
    /* some useful first stage params */
    usefulParams2.sftbasename = uvar_sftData2;
    usefulParams2.nStacks = uvar_nStacks2;
    usefulParams2.spinRange_startTime = spinRange_startTime2;
    usefulParams2.spinRange_endTime = spinRange_endTime2;
    usefulParams2.spinRange_refTime = spinRange_refTime;
    usefulParams2.edat = edat;
    usefulParams2.minStartTimeGPS = minStartTimeGPS2;
    usefulParams2.maxEndTimeGPS = maxEndTimeGPS2;
    usefulParams2.blocksRngMed = uvar_blocksRngMed;
    usefulParams2.Dterms = uvar_Dterms;
    usefulParams2.dopplerMax = uvar_dopplerMax;
    
    /* set reference time for pular parameters */
    if ( LALUserVarWasSet(&uvar_refTime)) 
      usefulParams2.refTime = uvar_refTime;
    else 
      usefulParams2.refTime = -1;
    
    LAL_CALL( SetUpSFTs( &status, &stackMultiSFT2, &stackMultiNoiseWeights2, 
			 &stackMultiDetStates2, &usefulParams2), &status);
    
    /* some useful params computed by SetUpSFTs */
    deltaF2 = usefulParams2.deltaF;
    timebase2 = 1.0/deltaF2;
    deltaFstack2 = usefulParams2.deltaFstack;
    tStack2 = 1.0/deltaFstack2;
    tObs2 = usefulParams2.tObs;
    nStacks2 = usefulParams2.nStacks;
    tStart2GPS = usefulParams2.tStartGPS;
    midTstack2 = usefulParams2.midTstack;
    startTstack2 = usefulParams2.startTstack;
    
    /* set reference time for calculating fstat */
    thisPoint2.refTime = tStart2GPS;
    /* binary and higher spindowns not implemented so far */
    thisPoint2.orbit = NULL;
    thisPoint2.fkdot[2] = 0.0;
    thisPoint2.fkdot[3] = 0.0;
              
  } /* end if(uvar_followup) */
  

  /* start setting up the semicoherent (Hough or stack slide) part of the search */
  
  /* set up some semiCoherent parameters */
  semiCohPar.tsMid = midTstack1;
  semiCohPar.refTime = tStart1GPS;
  semiCohPar.vel = velStack1;
  semiCohPar.pos = posStack1;
  semiCohPar.outBaseName = uvar_fnameout;
  semiCohPar.pixelFactor = uvar_pixelFactor;

  /* allocate memory for Hough candidates */
  semiCohCandList1.length = uvar_nCand1;
  semiCohCandList1.refTime = tStart1GPS;
  semiCohCandList1.nCandidates = 0; /* initialization */
  semiCohCandList1.minSigIndex = 0;
  semiCohCandList1.list = (SemiCohCandidate *)LALCalloc( 1, semiCohCandList1.length * sizeof(SemiCohCandidate));


  /*-----------Create template grid for first stage ---------------*/
  /* prepare initialization of DopplerScanner to step through paramter space */
  scanInit1.dAlpha = uvar_dAlpha;
  scanInit1.dDelta = uvar_dDelta;
  scanInit1.gridType = uvar_gridType;
  scanInit1.metricType = uvar_metricType;
  scanInit1.metricMismatch = uvar_mismatch1;
  scanInit1.projectMetric = TRUE;
  scanInit1.obsDuration = tStack1;
  scanInit1.obsBegin = startTstack1->data[ nStacks1/2 ];
  scanInit1.Detector = &(stackMultiDetStates1.data[0]->data[0]->detector);
  scanInit1.ephemeris = edat;
  scanInit1.skyGridFile = uvar_skyGridFile;
  scanInit1.searchRegion.Freq = uvar_Freq;
  scanInit1.searchRegion.FreqBand = uvar_FreqBand;
  scanInit1.searchRegion.f1dot = uvar_f1dot;
  scanInit1.searchRegion.f1dotBand = spinRange_startTime1->fkdotBand->data[1];;
  scanInit1.searchRegion.skyRegionString = (CHAR*)LALCalloc(1, strlen(uvar_skyRegion)+1);
  strcpy (scanInit1.searchRegion.skyRegionString, uvar_skyRegion);

  /* initialize skygrid  */  
  LAL_CALL ( InitDopplerScan ( &status, &thisScan1, &scanInit1), &status); 
  
  /* parameters for 2nd stage */
  /* set up parameters for second stage Fstat calculation */
  /* we don't set frequency and frequency band 
     -- these depend on candidates to be followed up */
  if ( uvar_followUp ) 
    {
      /* allocate memory for Fstat candidates */
      if ( LALUserVarWasSet(&uvar_fstatThr))
      	create_fstat_toplist(&fstatToplist, uvar_reallocBlock);
      else
      	create_fstat_toplist(&fstatToplist, uvar_nCand2);
      
     /* prepare initialization of DopplerScanner to step through paramter space */
      scanInit2.dAlpha = uvar_dAlpha;
      scanInit2.dDelta = uvar_dDelta;
      scanInit2.gridType = uvar_gridType;
      scanInit2.metricType = uvar_metricType;
      scanInit2.metricMismatch = uvar_mismatch2;
      scanInit2.projectMetric = TRUE;
      scanInit2.obsDuration = tStack2;
      scanInit2.obsBegin = startTstack2->data[ nStacks1/2 ];
      scanInit2.Detector = &(stackMultiDetStates2.data[0]->data[0]->detector);
      scanInit2.ephemeris = edat;  /* used by Ephemeris-based metric */
      scanInit2.skyGridFile = uvar_skyGridFile;
      /* the search region for scanInit2 will be set later 
	 -- it depends on the candidates from the first stage */
    } /* end if( uvar_followUp) */





  /*----- start main calculations by going over coarse grid points and 
          selecting candidates and following up --------*/
 
  /* loop over skygrid points */
  while(1)
    {
      UINT4 ifdot, nfdot;  /* counter and number of spindown values */
      REAL8 dfDot;  /* resolution in spindown */

      LAL_CALL (NextDopplerPos( &status, &dopplerpos1, &thisScan1 ), &status);
      if (thisScan1.state == STATE_FINISHED) /* scanned all DopplerPositions yet? */
	break;

      if (   lalDebugLevel )
	fprintf( stdout, "Coarse grid has %d points\n", thisScan1.numGridPoints);
      
      /*------------- calculate F statistic for each stack --------------*/
      
      /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
      thisPoint1.Alpha = dopplerpos1.Alpha;
      thisPoint1.Delta = dopplerpos1.Delta;
      
      /* number of fdot values */
      dfDot = thisScan1.df1dot;
      
      nfdot = (UINT4)( spinRange_startTime1->fkdotBand->data[1]/ dfDot + 0.5) + 1; 
      
      /* loop over fdot values */
      for ( ifdot=0; ifdot<nfdot; ifdot++)
	{
	  /* set spindown value for Fstat calculation */
  	  thisPoint1.fkdot[1] = spinRange_startTime1->fkdot->data[1] + ifdot * dfDot;
	  thisPoint1.fkdot[0] = spinRange_startTime1->fkdot->data[0];
	  	  
	  /* calculate the Fstatistic for each stack*/
	  for ( k = 0; k < nStacks1; k++) {
	    LAL_CALL( ComputeFStatFreqBand ( &status, fstatVector1.data + k, &thisPoint1, 
					     stackMultiSFT1.data[k], stackMultiNoiseWeights1.data[k], 
					     stackMultiDetStates1.data[k], &CFparams), &status);
	  }
	  
	  /* print fstat vector if required -- mostly for debugging */
	  if ( uvar_printFstat1 )
	    {
	      for (k=0; k<nStacks1; k++)
		LAL_CALL( PrintFstatVec ( &status, fstatVector1.data + k, fpFstat1, &thisPoint1, refTimeGPS, k+1), &status); 
	    }

	  
	  /*--------------- get candidates from a semicoherent search ---------------*/

	  /* the input to this section is the set of fstat vectors fstatVector1 and the 
	     parameters semiCohPar. The output is the list of candidates in semiCohCandList1 */

	  /* set sky location and spindown for Hough grid -- same as for Fstat calculation */	  
	  semiCohPar.alpha = thisPoint1.Alpha;
	  semiCohPar.delta = thisPoint1.Delta;
	  semiCohPar.fdot = thisPoint1.fkdot[1];
	  semiCohPar.nfdot = uvar_nfdot; /* look into this more carefully */
	  
	  /* the hough option */
	  /* select peaks */ 	      
	  if ( uvar_method == 0 ) {
	    LAL_CALL( FstatVectToPeakGram( &status, &pgV, &fstatVector1, uvar_peakThrF), &status);
	    	    
	    /* get candidates */
	    LAL_CALL ( ComputeFstatHoughMap ( &status, &semiCohCandList1, &pgV, &semiCohPar), &status);
	    
	    /* free peakgrams -- we don't need them now because we have the Hough maps */
	    for (k=0; k<nStacks1; k++) 
	      LALFree(pgV.pg[k].peak);
	    LALFree(pgV.pg);
	  }


	  /* --- stackslide option should come here --------*/

	  /* if ( uvar_method == 1 ) {.....} */	  


	  /* print candidates if desired */
	  if ( uvar_printCand1 ) {
	    LAL_CALL ( PrintSemiCohCandidates ( &status, &semiCohCandList1, fpSemiCoh), &status);
	  }
	  
	  /*------------- Follow up candidates --------------*/
	  
	  /* this part is more general than it has to be.
	     it is meant to be generalized to the case when 
	     nStacks2 is not necessarily 1 */
	  
	  /* check if user requested a follow up stage*/
	  if ( uvar_followUp ) 
	    {
	      
	      /* loop over candidates surviving 1st stage  */
	      for ( j=0; j < semiCohCandList1.nCandidates; j++) 
		{		  
		  /* 2nd stage frequency and spindown variables */
		  /* these are the parameters passed down from the 
		     first stage which is why they have the subscript 1 */
		  REAL8 fStart1, freqBand1, fdot1, fdotBand1;
		  REAL8 alpha1, delta1, alphaBand1, deltaBand1;

		  /* get spin range of candidate */
		  spinRange_Temp->epoch = tStart1GPS;
		  spinRange_Temp->fkdot->data[0] = semiCohCandList1.list[j].freq - 0.5 * semiCohCandList1.list[j].dFreq;
		  spinRange_Temp->fkdotBand->data[0] = semiCohCandList1.list[j].dFreq;
		  spinRange_Temp->fkdot->data[1] = semiCohCandList1.list[j].fdot - 0.5 * semiCohCandList1.list[j].dFdot;
		  spinRange_Temp->fkdotBand->data[1] = semiCohCandList1.list[j].dFdot;

		  /* extrapulate spin range to start time of second stack */
		  LAL_CALL( LALExtrapolatePulsarSpinRange(  &status, spinRange_Temp, tStart2GPS, spinRange_Temp ), &status);

		  /* set frequency and spindown ranges */
		  fStart1 = spinRange_Temp->fkdot->data[0];
		  freqBand1 = spinRange_Temp->fkdotBand->data[0];
		  fdot1 = spinRange_Temp->fkdot->data[1];
		  fdotBand1 = spinRange_Temp->fkdotBand->data[1];
		  
		  /* set sky region to be refined */
		  alpha1 = semiCohCandList1.list[j].alpha - 0.5 * semiCohCandList1.list[j].dAlpha;
		  delta1 = semiCohCandList1.list[j].delta - 0.5 * semiCohCandList1.list[j].dDelta;
		  alphaBand1 = semiCohCandList1.list[j].dAlpha;
		  deltaBand1 = semiCohCandList1.list[j].dDelta;
		  LAL_CALL (SkySquare2String( &status, &(scanInit2.searchRegion.skyRegionString),
					      alpha1, delta1, alphaBand1, deltaBand1), &status);
		  
		  /* set second doppler scan variables */
		  scanInit2.searchRegion.Freq = fStart1;
		  scanInit2.searchRegion.FreqBand = freqBand1;
		  scanInit2.searchRegion.f1dot = fdot1;
		  scanInit2.searchRegion.f1dotBand = fdotBand1;
		  
		  /* allocate fstat memory */
		  fstatVector2.length = nStacks2;
		  fstatVector2.data = NULL;
		  fstatVector2.data = (REAL8FrequencySeries *)LALCalloc( 1, nStacks2 * sizeof(REAL8FrequencySeries));
		  binsFstat2 = (UINT4)(freqBand1/deltaFstack2 + 0.5) + 1;
		  for (k = 0; k < nStacks2; k++) 
		    { 
		      fstatVector2.data[k].epoch = startTstack2->data[k];
		      fstatVector2.data[k].deltaF = deltaFstack2;
		      fstatVector2.data[k].f0 = fStart1;
		      fstatVector2.data[k].data = (REAL8Sequence *)LALCalloc( 1, sizeof(REAL8Sequence));
		      fstatVector2.data[k].data->length = binsFstat2;
		      fstatVector2.data[k].data->data = (REAL8 *)LALCalloc( 1, binsFstat2 * sizeof(REAL8));
		    } 


		  /* initialize skygrid  */  
		  LAL_CALL ( InitDopplerScan ( &status, &thisScan2, &scanInit2), &status); 
		  
		  
		  /* loop over fine skygrid points */
		  while(1)
		    {
		      UINT4 ifdot2, nfdot2;  /* counter and number of spindown values */
		      REAL8 dfDot2;  /* resolution in spindown */

		      LAL_CALL (NextDopplerPos( &status, &dopplerpos2, &thisScan2 ), &status);
		      if (thisScan2.state == STATE_FINISHED) /* scanned all DopplerPositions yet? */
			break;

		      if ( lalDebugLevel )
			fprintf( stdout, "Fine grid has %d points\n", thisScan2.numGridPoints);

		      
		      /*------------- calculate F statistic for each stack --------------*/
		      
		      /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
		      thisPoint2.Alpha = dopplerpos2.Alpha;
		      thisPoint2.Delta = dopplerpos2.Delta;
		      		      
		      /* number of fdot values */
		      dfDot2 = thisScan2.df1dot;
		      nfdot2 = (UINT4)( fdotBand1 / dfDot2 + 0.5) + 1; 
		      
		      /* loop over fdot values */
		      for ( ifdot2 = 0; ifdot2 < nfdot2; ifdot2++)
			{
			  
			  /* set spindown value for Fstat calculation */
			  thisPoint2.fkdot[1] = fdot1 + ifdot2 * dfDot2;
			  
			  /* calculate the Fstatistic */
			  for ( k = 0; k < nStacks2; k++) {
			    LAL_CALL( ComputeFStatFreqBand ( &status, fstatVector2.data + k, &thisPoint2, 
							     stackMultiSFT2.data[k], stackMultiNoiseWeights2.data[k], 
							     stackMultiDetStates2.data[k], &CFparams), &status);
			  }
			  
			  /* select candidates from 2nd stage */
			  for (k = 0; k < nStacks2; k++) 
			    {
			      if ( LALUserVarWasSet(&uvar_fstatThr))
				LAL_CALL( GetFstatCandidates( &status, fstatToplist, fstatVector2.data + k, uvar_fstatThr,
							      thisPoint2.Alpha, thisPoint2.Delta, 
							      thisPoint2.fkdot[1], uvar_reallocBlock ), &status);
			      else
				LAL_CALL( GetFstatCandidates_toplist( &status, fstatToplist, fstatVector2.data + k, 
								      thisPoint2.Alpha, thisPoint2.Delta, 
								      thisPoint2.fkdot[1] ), &status);
			      

			    } /* end loop over nstacks2 for selecting candidates */
			  
			  
			} /* end loop over refined spindown parameters */
		      
		    } /* end while loop over second stage refined sky grid */
		  

		  /* free Fstat vectors  */
		  for(k = 0; k < nStacks2; k++) {
		    LALFree(fstatVector2.data[k].data->data);
		    LALFree(fstatVector2.data[k].data);
		  }
		  LALFree(fstatVector2.data);


		  /* destroy dopplerscan2 variables */ 
		  LAL_CALL ( FreeDopplerScan(&status, &thisScan2), &status);
		  if ( scanInit2.searchRegion.skyRegionString )
		    LALFree ( scanInit2.searchRegion.skyRegionString );
		  
		  
		} /* end loop over candidates from 1st stage */
	      
	    } /* end block for follow-up stage */ 
	  
	} /* end loop over coarse grid fdot values */
      
    } /* end while loop over 1st stage coarse skygrid */
      

  /* print fstat candidates */  
  {
    UINT4 checksum;
    if ( uvar_followUp )
      if ( write_fstat_toplist_to_fp( fstatToplist, fpFstat, &checksum) < 0)
	fprintf( stderr, "Error in writing toplist to file\n");
    /*     LAL_CALL( AppendFstatCandidates( &status, &fStatCand, fpFstat), &status); */
  }
	 
  /*------------ free all remaining memory -----------*/
  
  /* free memory */

  if ( uvar_followUp )
    {
      LALFree(midTstack2->data);
      XLALDestroyTimestampVector(startTstack2);
      /* free sfts */
      for ( k = 0; k < nStacks2; k++) {
	LAL_CALL( LALDestroyMultiSFTVector ( &status, stackMultiSFT2.data + k), &status);
	LAL_CALL( LALDestroyMultiNoiseWeights ( &status, stackMultiNoiseWeights2.data + k), &status);
	XLALDestroyMultiDetectorStateSeries ( stackMultiDetStates2.data[k] );
      }
      LALFree(stackMultiSFT2.data);
      LALFree(stackMultiNoiseWeights2.data);
      LALFree(stackMultiDetStates2.data);
           
      if ( fstatToplist ) 
	free_fstat_toplist(&fstatToplist);
      
      fclose(fpFstat);
      LALFree(fnameFstatCand);

      LAL_CALL(LALDestroyTimestampVector ( &status, &sftTimeStamps2), &status);  	
    }

  if ( uvar_printCand1 )
    {
      LALFree(fnameSemiCohCand);
      fclose(fpSemiCoh);
    }
  
  if ( uvar_printFstat1 )
    {
      fclose(fpFstat1);
      LALFree( fnameFstatVec1 );
    }
  
  /* free first stage memory */
  for ( k = 0; k < nStacks1; k++) {
    LAL_CALL( LALDestroyMultiSFTVector ( &status, stackMultiSFT1.data + k), &status);
    LAL_CALL( LALDestroyMultiNoiseWeights ( &status, stackMultiNoiseWeights1.data + k), &status);
    XLALDestroyMultiDetectorStateSeries ( stackMultiDetStates1.data[k] );
  }
  LALFree(stackMultiSFT1.data);
  LALFree(stackMultiNoiseWeights1.data);
  LALFree(stackMultiDetStates1.data);
  
  LAL_CALL(LALDestroyTimestampVector ( &status, &sftTimeStamps1), &status);  	  

  XLALDestroyTimestampVector(midTstack1);
  XLALDestroyTimestampVector(startTstack1);

  /* free Fstat vectors  */
  for(k = 0; k < nStacks1; k++) {
    LALFree(fstatVector1.data[k].data->data);
    LALFree(fstatVector1.data[k].data);
  }
  LALFree(fstatVector1.data);
  

  /* free Vel/Pos vectors and ephemeris */
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);
  LAL_CALL( LALDDestroyVectorSequence (&status,  &velStack1), &status);
  LAL_CALL( LALDDestroyVectorSequence (&status,  &posStack1), &status);
  
  /* free dopplerscan stuff */
  LAL_CALL ( FreeDopplerScan(&status, &thisScan1), &status);
  if ( scanInit1.searchRegion.skyRegionString )
    LALFree ( scanInit1.searchRegion.skyRegionString );

 
  /* free candidates */
  LALFree(semiCohCandList1.list);

  XLALDestroyPulsarSpinRange(spinRange_endTime1);
  XLALDestroyPulsarSpinRange(spinRange_endTime2);
  XLALDestroyPulsarSpinRange(spinRange_startTime1);
  XLALDestroyPulsarSpinRange(spinRange_startTime2);
  XLALDestroyPulsarSpinRange(spinRange_Temp);

  LAL_CALL (LALDestroyUserVars(&status), &status);  

  LALCheckMemoryLeaks();

  return HIERARCHICALSEARCH_ENORM;
} /* main */





/** Set up stacks, read SFTs, calculate SFT noise weights and calculate 
    detector-state */
void SetUpSFTs( LALStatus *status,
		MultiSFTVectorSequence *stackMultiSFT,
		MultiNoiseWeightsSequence *stackMultiNoiseWeights,
		MultiDetectorStateSeriesSequence *stackMultiDetStates,
		UsefulStageVariables *in)
{

  SFTCatalog *catalog = NULL;
  static SFTConstraints constraints;
  REAL8 timebase, tStack, tObs, refTime, deltaF, S_hat;
  UINT4 nStacks, k;
  LIGOTimeGPS tStartGPS, tEndGPS, refTimeGPS;
  LIGOTimeGPSVector *sftTimeStamps=NULL; 
  SFTCatalogSequence catalogSeq;
  
  REAL8 doppWings, fMin, fMax;
  REAL8 startTime_freqLo, startTime_freqHi;
  REAL8 endTime_freqLo, endTime_freqHi;
  REAL8 freqLo, freqHi;
  INT4 extraBins, tmpLeap;
  
  /* leap second for LALBarycenter */
  LALLeapSecFormatAndAcc lsfas = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
    
  INITSTATUS( status, "SetUpSFTs", rcsid );
  ATTATCHSTATUSPTR (status);

  /* get sft catalog */
  constraints.startTime = &(in->minStartTimeGPS);
  constraints.endTime = &(in->maxEndTimeGPS);
  TRY( LALSFTdataFind( status->statusPtr, &catalog, in->sftbasename, &constraints), status);
  
  /*   if ( uvar_followUp ) { */
  /*     TRY ( PrintCatalogInfo ( status->statusPtr, catalog, fpFstat), status); */
  /*   } */
  
  /*   if ( uvar_printFstat1 ) { */
  /*     TRY ( PrintCatalogInfo ( status->statusPtr, catalog, fpFstat1), status); */
  /*   } */
  
  /*   if ( uvar_printCand1 ) { */
  /*     TRY ( PrintCatalogInfo ( status->statusPtr, catalog, fpSemiCoh), status); */
  /*   } */
  
  /* set some sft parameters */
  in->deltaF = deltaF = catalog->data->header.deltaF;
  timebase = 1.0/deltaF;
  
  /* get sft catalogs for each stack */
  TRY( SetUpStacks( status->statusPtr, &catalogSeq, &tStack, &sftTimeStamps, catalog, in->nStacks), status);
  
  /*   if ( uvar_followUp ) { */
  /*     TRY ( PrintStackInfo ( status->statusPtr, &catalogSeq, fpFstat), status); */
  /*   } */
  
  /*   if ( uvar_printFstat1 ) { */
  /*     TRY ( PrintStackInfo ( status->statusPtr, &catalogSeq, fpFstat1), status); */
  /*   } */
  
  /*   if ( uvar_printCand1 ) { */
  /*     TRY ( PrintStackInfo ( status->statusPtr, &catalogSeq, fpSemiCoh), status); */
  /*   } */
  
  /* calculate start and end times and tobs */
  tStartGPS = sftTimeStamps->data[0];
  tEndGPS = sftTimeStamps->data[sftTimeStamps->length - 1];
  TRY( LALAddFloatToGPS( status->statusPtr, &tEndGPS, &tEndGPS, timebase ), status);
  TRY ( LALDeltaFloatGPS ( status->statusPtr, &tObs, &tEndGPS, &tStartGPS), status);

  /* Leap seconds for the first timestamp */   
  TRY( LALLeapSecs( status->statusPtr, &tmpLeap, &tStartGPS, &lsfas), status);
  in->edat->leap = (INT2)tmpLeap;

  in->tStartGPS = tStartGPS;
  in->tObs = tObs;
  
  /* use stacks info to calculate search frequency resolution */
  in->deltaFstack = 1.0/tStack;
  
  /* set reference time for pular parameters */
  refTime = in->refTime;
  if ( refTime > 0)  {
      TRY ( LALFloatToGPS( status->statusPtr, &refTimeGPS, &refTime), status);
  }
  else {
    refTimeGPS = tStartGPS;
  }
  
  /* get frequency and fdot bands at start time of sfts by extrapolating from reftime */
  in->spinRange_refTime->epoch = refTimeGPS;
  TRY( LALExtrapolatePulsarSpinRange( status->statusPtr, in->spinRange_startTime, tStartGPS, in->spinRange_refTime), status); 
  TRY( LALExtrapolatePulsarSpinRange( status->statusPtr, in->spinRange_endTime, tEndGPS, in->spinRange_refTime), status); 
       
  /* set wings of sfts to be read */
  /* the wings must be enough for the Doppler shift and extra bins
     for the running median block size and Dterms for Fstat calculation.
     In addition, it must also include wings for the spindown correcting 
     for the reference time  */
  /* calculate Doppler wings at the highest frequency */
  startTime_freqLo = in->spinRange_startTime->fkdot->data[0]; /* lowest search freq at start time */
  startTime_freqHi = startTime_freqLo + in->spinRange_startTime->fkdotBand->data[0]; /* highest search freq. at start time*/
  endTime_freqLo = in->spinRange_endTime->fkdot->data[0];
  endTime_freqHi = endTime_freqLo + in->spinRange_endTime->fkdotBand->data[0];
  
  freqLo = HSMIN ( startTime_freqLo, endTime_freqLo );
  freqHi = HSMAX ( startTime_freqHi, endTime_freqHi );
  doppWings = freqHi * in->dopplerMax;    /* maximum Doppler wing -- probably larger than it has to be */
  extraBins = HSMAX ( in->blocksRngMed/2 + 1, in->Dterms );
  
  fMin = freqLo - doppWings - extraBins * deltaF; 
  fMax = freqHi + doppWings + extraBins * deltaF; 
  
  
  /* set up vector containing mid times of stacks */    
  in->midTstack =  XLALCreateTimestampVector ( in->nStacks );
  
  /* set up vector containing start times of stacks */    
  in->startTstack =  XLALCreateTimestampVector ( in->nStacks );
  
  /* finally memory for stack of multi sfts */
  stackMultiSFT->length = in->nStacks;
  stackMultiSFT->data = (MultiSFTVector **)LALCalloc(1, in->nStacks * sizeof(MultiSFTVector *));
  
  stackMultiNoiseWeights->length = in->nStacks;
  stackMultiNoiseWeights->data = (MultiNoiseWeights **)LALCalloc(1, in->nStacks * sizeof(MultiNoiseWeights *));
  
  stackMultiDetStates->length = in->nStacks;
  stackMultiDetStates->data = (MultiDetectorStateSeries **)LALCalloc(1, in->nStacks * sizeof(MultiDetectorStateSeries *));
  
  for (nStacks = 0, k = 0; k < in->nStacks; k++) {
    /* nStacks is number of non-empty stacks while
       in->nStacks1 is what the user specified */
    
    if ( catalogSeq.data[k].length > 0 ){
      /* if the stack is non-empty */
      
      MultiPSDVector *psd = NULL;	
      LIGOTimeGPSVector *startTsftStack = NULL; /* sft time stamps for sfts in current stack */
      
      /* get SFT timestamps in stack */
      TRY( LALSFTtimestampsFromCatalog(  status->statusPtr, &startTsftStack, catalogSeq.data + k ), status);  	
      
      /* start time of stack */
      in->startTstack->data[nStacks] = startTsftStack->data[0];
      
      /* mid time of stack */          
      TRY( LALAddFloatToGPS( status->statusPtr, in->midTstack->data + nStacks, in->startTstack->data + nStacks,  
				    0.5 * tStack ), status);
      
      /* load the sfts */
      TRY( LALLoadMultiSFTs ( status->statusPtr, stackMultiSFT->data + nStacks,  catalogSeq.data + k, 
				   fMin, fMax ), status);
      
      /* normalize sfts and compute noise weights and detector state */
      TRY( LALNormalizeMultiSFTVect ( status->statusPtr, &psd, stackMultiSFT->data[nStacks], 
					   in->blocksRngMed ), status );
      
      TRY( LALComputeMultiNoiseWeights  ( status->statusPtr, stackMultiNoiseWeights->data + nStacks, &S_hat, 
					       psd, in->blocksRngMed, 0 ), status );
	
      TRY ( LALGetMultiDetectorStates ( status->statusPtr, stackMultiDetStates->data + nStacks, 
					     stackMultiSFT->data[nStacks], in->edat ), status );

      TRY ( LALDestroyMultiPSDVector ( status->statusPtr, &psd ), status );
      
      TRY ( LALDestroyTimestampVector ( status->statusPtr, &startTsftStack), status);
      
      nStacks++;
      
      } /* if ( catalogSeq->data[k].length > 0 ) */
  } /* loop over k */
  
  
  /* realloc if nStacks != in->nStacks */
  if ( in->nStacks > nStacks ) {
    
    in->midTstack->length = nStacks;
    in->midTstack->data = (LIGOTimeGPS *)LALRealloc( in->midTstack->data, nStacks * sizeof(LIGOTimeGPS));
    
    in->startTstack->length = nStacks;
    in->startTstack->data = (LIGOTimeGPS *)LALRealloc( in->startTstack->data, nStacks * sizeof(LIGOTimeGPS));
    
    stackMultiSFT->length = nStacks;
    stackMultiSFT->data = (MultiSFTVector **)LALRealloc( stackMultiSFT->data, nStacks * sizeof(MultiSFTVector *));
    
  }  /* finish realloc */
    
  /* we don't need the original catalog */
  TRY( LALDestroySFTCatalog( status->statusPtr, &catalog ), status);  	
  
  /* free catalog sequence */
  for (k = 0; k < in->nStacks; k++)
    {
      if ( catalogSeq.data[k].length > 0 ) {
	LALFree(catalogSeq.data[k].data);
      } /* end if */
    } /* loop over stacks */
  LALFree( catalogSeq.data);  

  in->nStacks = nStacks;  
  
  DETATCHSTATUSPTR (status);
  RETURN(status);
  
} /* SetUpSFTs */



/** \brief Function for calculating Hough Maps and candidates 
    \param pgV is a HOUGHPeakGramVector obtained after thresholding Fstatistic vectors
    \param params is a pointer to HoughParams -- parameters for calculating Hough maps
    \out houghCand Candidates from thresholding Hough number counts

    This function takes a peakgram as input. This peakgram was constructed
    by setting a threshold on a sequence of Fstatistic vectors.  The function 
    produces a Hough map in the sky for each value of the frequency and spindown.
    The Hough nummber counts are then used to select candidates in 
    parameter space to be followed up in a more refined search.
    This uses DriveHough_v3.c as a prototype suitably modified to work 
    on demodulated data instead of SFTs.  
*/
void ComputeFstatHoughMap(LALStatus *status,
			  SemiCohCandidateList  *out,   /* output candidates */
			  HOUGHPeakGramVector *pgV, /* peakgram vector */
			  SemiCoherentParams *params)
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
  INT4  k, iHmap, nSpin1Max, nSpin1Min, nStacks, nfdot;
  REAL8 deltaF, alpha, delta;
  REAL8 patchSizeX, patchSizeY, f1jump;
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

  INITSTATUS( status, "ComputeFstatHoughMap", rcsid );
  ATTATCHSTATUSPTR (status);

  /* copy some parameters from peakgram vector */
  deltaF = pgV->pg->deltaF;
  nStacks = pgV->length;

  /* copy some params to local variables */
  nfdot = params->nfdot;
  fBinIni = pgV->pg[0].fBinIni;
  fBinFin = pgV->pg[0].fBinFin;
  alpha = params->alpha;
  delta = params->delta;
  vel = params->vel;
  pos = params->pos;
  fdot = params->fdot;
  tsMid = params->tsMid;
  refTimeGPS = params->refTime;
  TRY ( LALGPStoFloat( status->statusPtr, &refTime, &refTimeGPS), status);

  /* set patch size */
  /* this is supposed to be the "educated guess" 
     delta theta = 1.0 / (Tcoh * f0 * Vepi )
     where Tcoh is coherent time baseline, 
     f0 is frequency and Vepi is rotational velocity 
     of detector */
  patchSizeX = 0.5 / ( fBinIni * VEPI ); 
  patchSizeY = 0.5 / ( fBinIni * VEPI ); 

  /*--------------- first memory allocation --------------*/
  /* look up table vector */
  lutV.length = nStacks;
  lutV.lut = NULL;
  lutV.lut = (HOUGHptfLUT *)LALMalloc(nStacks*sizeof(HOUGHptfLUT));
  
  /* partial hough map derivative vector */
  phmdVS.length  = nStacks;
  phmdVS.nfSize  = 2*nfdot + 1;
  phmdVS.deltaF  = deltaF;
  phmdVS.phmd = NULL;
  phmdVS.phmd=(HOUGHphmd *)LALMalloc( phmdVS.length * phmdVS.nfSize *sizeof(HOUGHphmd));
  
  /* residual spindown trajectory */
  freqInd.deltaF = deltaF;
  freqInd.length = nStacks;
  freqInd.data = NULL;
  freqInd.data =  ( UINT8 *)LALMalloc(nStacks*sizeof(UINT8));
   
  /* resolution in space of residual spindowns */
  ht.dFdot.length = 1;
  ht.dFdot.data = NULL;
  ht.dFdot.data = (REAL8 *)LALMalloc( ht.dFdot.length * sizeof(REAL8));

  /* the residual spindowns */
  ht.spinRes.length = 1;
  ht.spinRes.data = NULL;
  ht.spinRes.data = (REAL8 *)LALMalloc(ht.spinRes.length*sizeof(REAL8));

  /* the residual spindowns */
  ht.spinDem.length = 1;
  ht.spinDem.data = NULL;
  ht.spinDem.data = (REAL8 *)LALMalloc(ht.spinRes.length*sizeof(REAL8));

  /* the demodulation params */
  parDem.deltaF = deltaF;
  parDem.skyPatch.alpha = alpha;
  parDem.skyPatch.delta = delta;
  parDem.spin.length = 1;
  parDem.spin.data = NULL;
  parDem.spin.data = (REAL8 *)LALCalloc(1, sizeof(REAL8));
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
    hist.data = (UINT8 *)LALMalloc((nStacks+1)*sizeof(UINT8));
    histTotal.data = (UINT8 *)LALMalloc((nStacks+1)*sizeof(UINT8));
    { 
      UINT4   j;
      for(j=0; j< histTotal.length; ++j) 
	histTotal.data[j]=0; 
    }
    fileStats = NULL;
    fileStats = (CHAR *)LALCalloc(1, 256 * sizeof(CHAR));
    strcpy( fileStats, params->outBaseName);
    strcat( fileStats, "stats");
    if ( !(fpStats = fopen(fileStats, "w")))
      {
	fprintf(stderr, "Unable to open file '%s' for writing...continuing\n", fileStats);
      }
  }


  /* calculate time differences from start of observation time */
  TRY( LALDCreateVector( status->statusPtr, &timeDiffV, nStacks), status);

  for (k=0; k<nStacks; k++) {
    REAL8 tMidStack;

    TRY ( LALGPStoFloat ( status->statusPtr, &tMidStack, tsMid->data + k), status);
    timeDiffV->data[k] = tMidStack - refTime;
  }

  /* if there are residual spindowns */
  f1jump = 1.0 / timeDiffV->data[nStacks - 1]; /* resolution in residual fdot */
 
  /* adjust fBinIni and fBinFin to take maxNBins into account */
  /* and make sure that we have fstat values for sufficient number of bins */
  parRes.f0Bin =  fBinIni;      
    /* need a function to more accurately calculate maxNBins apriori */
  TRY( LALHOUGHComputeSizePar( status->statusPtr, &parSize, &parRes ),  status );
  fBinIni += parSize.maxNBins;
  fBinFin -= parSize.maxNBins;

  ASSERT ( fBinIni < fBinFin, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );

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
    patch.xCoor = (REAL8 *)LALMalloc(xSide*sizeof(REAL8));
    patch.yCoor = (REAL8 *)LALMalloc(ySide*sizeof(REAL8));
    TRY( LALHOUGHFillPatchGrid( status->statusPtr, &patch, &parSize ), status );
    
    /*------------- other memory allocation and settings----------------- */
    for(j=0; j<lutV.length; ++j){
      lutV.lut[j].maxNBins = maxNBins;
      lutV.lut[j].maxNBorders = maxNBorders;
      lutV.lut[j].border =
	(HOUGHBorder *)LALMalloc(maxNBorders*sizeof(HOUGHBorder));
      lutV.lut[j].bin =
	(HOUGHBin2Border *)LALMalloc(maxNBins*sizeof(HOUGHBin2Border));
      for (i=0; i<maxNBorders; ++i){
	lutV.lut[j].border[i].ySide = ySide;
	lutV.lut[j].border[i].xPixel =
	  (COORType *)LALMalloc(ySide*sizeof(COORType));
      }
    }

    for(j=0; j<phmdVS.length * phmdVS.nfSize; ++j){
      phmdVS.phmd[j].maxNBorders = maxNBorders;
      phmdVS.phmd[j].leftBorderP =
	(HOUGHBorder **)LALMalloc(maxNBorders*sizeof(HOUGHBorder *));
      phmdVS.phmd[j].rightBorderP =
	(HOUGHBorder **)LALMalloc(maxNBorders*sizeof(HOUGHBorder *));
      phmdVS.phmd[j].ySide = ySide;
      phmdVS.phmd[j].firstColumn = NULL;
      phmdVS.phmd[j].firstColumn = (UCHAR *)LALMalloc(ySide*sizeof(UCHAR));
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
      TRY( LALHOUGHParamPLUT( status->statusPtr, &parLut, &parSize, &parDem),status );
      /* build the LUT */
      TRY( LALHOUGHConstructPLUT( status->statusPtr, &(lutV.lut[j]), &patch, &parLut ),
	   status );
    }
    
    /*--------- build the set of  PHMD centered around fBin -------------*/     
    phmdVS.fBinMin = fBin - nfdot;
    TRY( LALHOUGHConstructSpacePHMD(status->statusPtr, &phmdVS, pgV, &lutV), status );
    
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
    ht.dFdot.data[0] = deltaF * f1jump;
    ht.map   = NULL;
    ht.map   = (HoughTT *)LALMalloc(xSide*ySide*sizeof(HoughTT));
    TRY( LALHOUGHInitializeHT( status->statusPtr, &ht, &patch), status); /*not needed */
    
    /*  Search frequency interval possible using the same LUTs */
    fBinSearch = fBin;
    fBinSearchMax = fBin + parSize.nFreqValid - 1 - nfdot; /* do we need the -nfdot here? */
     
    /* Study all possible frequencies with one set of LUT */    
    while ( (fBinSearch <= fBinFin) && (fBinSearch < fBinSearchMax) )  { 

      /* finally we can construct the hough maps and select candidates */
      {
	INT4   n;
	REAL8  f1dis;

	ht.f0Bin = fBinSearch;

	nSpin1Max = (nfdot < fBinFin - fBinSearch)? nfdot : fBinFin - fBinSearch;
	nSpin1Min = (nfdot < fBinSearch - fBinIni)? nfdot : fBinSearch - fBinIni;

	/*loop over all values of residual spindown */
	for( n = -nSpin1Min; n <= nSpin1Max; n++ ){ 
	  f1dis =  n*f1jump;

	  ht.spinRes.data[0] =  f1dis*deltaF;
	  
	  for (j=0; j < (UINT4)nStacks; j++)
	    freqInd.data[j] = fBinSearch + floor(timeDiffV->data[j]*f1dis + 0.5);
	  
	  TRY( LALHOUGHConstructHMT(status->statusPtr, &ht, &freqInd, &phmdVS),status );

	  /* get candidates */
	  /* 	  if ( selectCandThr ) { */
	  /* 	    TRY(GetHoughCandidates( status->statusPtr, out, &ht, &patch, &parDem, houghThr), status); */
	  /* 	  } */
	  /* 	  else { */
	  /* 	    TRY(GetHoughCandidates_toplist( status->statusPtr, out, &ht, &patch, &parDem), status);  */
	  /* 	  } */

	  TRY(GetHoughCandidates_toplist( status->statusPtr, out, &ht, &patch, &parDem), status);

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
	  
	  /* increment hough map index */ 	  
	  ++iHmap;
	  
	} /* end loop over spindown trajectories */

      } /* end of block for calculating total hough maps */
      

      /*------ shift the search freq. & PHMD structure 1 freq.bin -------*/
      ++fBinSearch;
      TRY( LALHOUGHupdateSpacePHMDup(status->statusPtr, &phmdVS, pgV, &lutV), status );
      
    }   /* closing while loop over fBinSearch */
    
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

/** \brief Function for selecting frequency bins from a set of Fstatistic vectors
    \param FstatVect : sequence of Fstatistic vectors
    \param thr is a REAL8 threshold for selecting frequency bins
    \return pgV : a vector of peakgrams 

    Input is a vector of Fstatistic vectors.  It allocates memory 
    for the peakgrams based on the frequency span of the Fstatistic vectors
    and fills tyem up by setting a threshold on the Fstatistic.  Peakgram must be 
    deallocated outside the function.
*/
void FstatVectToPeakGram (LALStatus *status,
			  HOUGHPeakGramVector *pgV,
			  REAL8FrequencySeriesVector *FstatVect,
			  REAL8  thr)
{
  INT4 j, k;
  INT4 nStacks, nSearchBins, nPeaks;
  UCHAR *upg;  

  INITSTATUS( status, "FstatVectToPeakGram", rcsid );
  ATTATCHSTATUSPTR (status);

  nStacks = FstatVect->length;
  nSearchBins = FstatVect->data->data->length;


  /* first memory allocation */  
  pgV->length = nStacks;
  pgV->pg = (HOUGHPeakGram *)LALMalloc( nStacks * sizeof(HOUGHPeakGram));


  upg = (UCHAR *)LALMalloc( nSearchBins * sizeof(UCHAR));

  /* loop over each stack and set peakgram */
  for (k=0; k<nStacks; k++) {
    INT4 *pInt; /* temporary pointer */
    REAL8 *pV;  /* temporary pointer */
    REAL8 f0, deltaF;
    pV = FstatVect->data[k].data->data;

    /* loop over Fstat vector, count peaks, and set upg values */
    nPeaks = 0;
    for(j=0; j<nSearchBins; j++) {
      if ( pV[j] > thr ) {
	nPeaks++;	
	upg[j] = 1; 
      }
      else
	upg[j] = 0;
    }

    /* fix length of peakgram and allocate memory appropriately */
    pgV->pg[k].length = nPeaks; 
    pgV->pg[k].peak = (INT4 *)LALMalloc( nPeaks * sizeof(INT4)); 

    /* fill up other peakgram parameters */
    pgV->pg[k].deltaF = FstatVect->data[k].deltaF;
    f0 = FstatVect->data[k].f0;
    deltaF = FstatVect->data[k].deltaF;
    pgV->pg[k].fBinIni = floor( f0/deltaF + 0.5);
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



void SetUpStacks(LALStatus *status, 
		 SFTCatalogSequence  *out,
		 REAL8 *tStack,
		 LIGOTimeGPSVector **ts,
		 SFTCatalog  *in,
		 UINT4 nStacks)
{

  UINT4 nSFTs, length, j, k;
  REAL8 timeBase, deltaF, tStart, tEnd, thisTime;


  INITSTATUS( status, "SetUpStacks", rcsid );
  ATTATCHSTATUSPTR (status);

  out->length = nStacks;
  out->data = (SFTCatalog *)LALCalloc( 1, nStacks * sizeof(SFTCatalog));
  
  nSFTs = in->length;

  ASSERT ( nStacks < nSFTs, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );

  TRY( LALSFTtimestampsFromCatalog( status->statusPtr, ts, in), status);  	

  deltaF = in->data[0].header.deltaF;
  timeBase = 1.0/deltaF;

  TRY ( LALGPStoFloat ( status->statusPtr, &tStart, ts[0]->data), status);
  TRY ( LALGPStoFloat ( status->statusPtr, &tEnd, ts[0]->data + nSFTs - 1), status);
  *tStack = (tEnd + timeBase - tStart) / nStacks;

  /* loop over the sfts and find out if it belongs to the k^th stack */
  for( j = 0; j < nSFTs; j++) {

    /* calculate which stack the j^th SFT belongs to */
    TRY ( LALGPStoFloat ( status->statusPtr, &thisTime, ts[0]->data + j), status);
    k = (UINT4)((thisTime - tStart)/(*tStack));

    out->data[k].length += 1;    

    length = out->data[k].length;
    
    out->data[k].data = (SFTDescriptor *)LALRealloc( out->data[k].data, length * sizeof(SFTDescriptor));

    out->data[k].data[length - 1] = in->data[j];   

  } /* loop over sfts */


  DETATCHSTATUSPTR (status);
  RETURN(status);
}






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
  
  INITSTATUS( status, "PrintHmap2file", rcsid );
  ATTATCHSTATUSPTR (status);

  strcpy(  filename, fnameOut);
  sprintf( filenumber, ".%06d",iHmap); 
  strcat(  filename, filenumber);

  fp=fopen(filename,"w");  
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


/** Get Hough candidates */
void GetHoughCandidates(LALStatus *status,
			SemiCohCandidateList *semiCohCand,
			HOUGHMapTotal *ht,
			HOUGHPatchGrid  *patch,
			HOUGHDemodPar   *parDem,
			REAL8 houghThreshold)
{
  REAL8UnitPolarCoor sourceLocation;
  REAL8 deltaF, f0, fdot, dFdot, patchSizeX, patchSizeY;
  INT8 f0Bin;  
  INT4 nCandidates, i,j, xSide, ySide;
  
  INITSTATUS( status, "GetHoughCandidates", rcsid );
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

  for (i = 0; i < ySide; i++){
    for (j = 0; j < xSide; j++){ 
      /* if threshold is exceeded then add candidate */
      if ( ht->map[i*xSide + j] > houghThreshold ) {

	nCandidates = semiCohCand->nCandidates;

	/* if there isn't enough memory then realloc */
	if ( nCandidates >= semiCohCand->length ) {
	  semiCohCand->length += 5000;
	  semiCohCand->list = (SemiCohCandidate *)LALRealloc( semiCohCand->list, semiCohCand->length * sizeof(SemiCohCandidate));
	} /* end of reallocs */

	semiCohCand->list[nCandidates].freq = f0;
	semiCohCand->list[nCandidates].dFreq = deltaF;
	
	semiCohCand->list[nCandidates].fdot = fdot;
	semiCohCand->list[nCandidates].dFdot = dFdot;

	/* get sky location of pixel */
	TRY( LALStereo2SkyLocation (status->statusPtr, &sourceLocation, 
				    j, i, patch, parDem), status);

	semiCohCand->list[nCandidates].alpha = sourceLocation.alpha;
	semiCohCand->list[nCandidates].delta = sourceLocation.delta;

	semiCohCand->list[nCandidates].dAlpha = patchSizeX / ((REAL8)xSide);
	semiCohCand->list[nCandidates].dDelta = patchSizeY / ((REAL8)ySide);

	/* increment candidate count */
	semiCohCand->nCandidates += 1;

      } /* end if statement */
    } /* end loop over xSide */
  } /* end loop over ySide */
  DETATCHSTATUSPTR (status);
  RETURN(status);
}


/** Get Hough candidates */
void GetHoughCandidates_toplist(LALStatus *status,
				SemiCohCandidateList *semiCohCand,
				HOUGHMapTotal *ht,
				HOUGHPatchGrid  *patch,
				HOUGHDemodPar   *parDem)
{
  REAL8UnitPolarCoor sourceLocation;
  REAL8 deltaF, f0, fdot, dFdot, patchSizeX, patchSizeY;
  INT8 f0Bin;  
  INT4 nCandidates, maxCandidates; 
  INT4 i,j, xSide, ySide, minSigIndex;

  INITSTATUS( status, "GetHoughCandidates_toplist", rcsid );
  ATTATCHSTATUSPTR (status);

  deltaF = ht->deltaF;
  f0Bin = ht->f0Bin;
  f0 = f0Bin * deltaF;

  nCandidates = semiCohCand->nCandidates;
  maxCandidates = semiCohCand->length;

  fdot = ht->spinDem.data[0] + ht->spinRes.data[0];
  dFdot = ht->dFdot.data[0];

  xSide = ht->xSide;
  ySide = ht->ySide;

  patchSizeX = ht->patchSizeX;
  patchSizeY = ht->patchSizeY;

  for (i = 0; i < ySide; i++)
    {
      for (j = 0; j < xSide; j++)
	{ 
	  REAL8 tempSig = ht->map[i*xSide + j];

	  if ( nCandidates < maxCandidates )
	    {
	      semiCohCand->list[nCandidates].freq = f0;
	      semiCohCand->list[nCandidates].dFreq = deltaF;
	      
	      semiCohCand->list[nCandidates].fdot = fdot;
	      semiCohCand->list[nCandidates].dFdot = dFdot;
	      
	      /* get sky location of pixel */
	      TRY( LALStereo2SkyLocation (status->statusPtr, &sourceLocation, 
					  j, i, patch, parDem), status);

	      semiCohCand->list[nCandidates].alpha = sourceLocation.alpha;
	      semiCohCand->list[nCandidates].delta = sourceLocation.delta;
	      semiCohCand->list[nCandidates].dAlpha = patchSizeX / ((REAL8)xSide);
	      semiCohCand->list[nCandidates].dDelta = patchSizeY / ((REAL8)ySide);
	      semiCohCand->list[nCandidates].significance = tempSig;

	      TRY ( GetMinSigIndex_toplist ( status->statusPtr, &minSigIndex, semiCohCand), status);		
	      semiCohCand->minSigIndex = minSigIndex;

	      /* increment candidate count */
	      semiCohCand->nCandidates += 1;
	      nCandidates = semiCohCand->nCandidates;

	    }
	  else 
	    {
	      /* if event is more significant than least significant 
		 event stored in toplist then replace it */
	      minSigIndex = semiCohCand->minSigIndex;
	      if ( tempSig > semiCohCand->list[minSigIndex].significance ) 
		{
		  semiCohCand->list[minSigIndex].freq = f0;
		  semiCohCand->list[minSigIndex].dFreq = deltaF;
		  semiCohCand->list[minSigIndex].fdot = fdot;
		  semiCohCand->list[minSigIndex].dFdot = dFdot;
	  
		  /* get sky location of pixel */
		  TRY( LALStereo2SkyLocation (status->statusPtr, &sourceLocation, 
					      j, i, patch, parDem), status);
		  
		  semiCohCand->list[minSigIndex].alpha = sourceLocation.alpha;
		  semiCohCand->list[minSigIndex].delta = sourceLocation.delta;  
		  semiCohCand->list[minSigIndex].dAlpha = patchSizeX / ((REAL8)xSide);
		  semiCohCand->list[minSigIndex].dDelta = patchSizeY / ((REAL8)ySide);
		  semiCohCand->list[minSigIndex].significance = tempSig;

		  /* find index of new least significant event */
		  TRY ( GetMinSigIndex_toplist ( status->statusPtr, &minSigIndex, semiCohCand), status);
		  semiCohCand->minSigIndex = minSigIndex;

		} /* end if statement for replacing least significant candidate */

	    } /* end else part of main if statement */

	} /* end loop over xSide */

    } /* end loop over ySide */

  DETATCHSTATUSPTR (status);
  RETURN(status);
}






/** Get least significant Hough candidate index */
void GetMinSigIndex_toplist(LALStatus *status,
			    INT4 *minSigIndex,
			    SemiCohCandidateList *semiCohCand)
{

  REAL8 minSig;
  INT4 k, nCandidates = semiCohCand->nCandidates;

  INITSTATUS( status, "GetMinSigIndex_toplist", rcsid );
  ATTATCHSTATUSPTR (status);


  *minSigIndex = 0;
  minSig = semiCohCand->list[0].significance;
  for (k = 1; k < nCandidates; k++)
    {
      REAL8 tempSig = semiCohCand->list[k].significance;
      if ( tempSig < minSig )
	{
	  *minSigIndex = k;
	  minSig = tempSig;
	}
    } /* end loop over candidates */
  

  DETATCHSTATUSPTR (status);
  RETURN(status);
}






/** Print Hough candidates */
void PrintSemiCohCandidates(LALStatus *status,
			  SemiCohCandidateList *in,
			  FILE *fp)
{
  INT4 k;

  INITSTATUS( status, "PrintSemiCohCandidates", rcsid );
  ATTATCHSTATUSPTR (status);

 
  for (k=0; k < in->nCandidates; k++)
    fprintf(fp, "%e   %e   %g   %g   %g   %g   %e   %e\n", in->list[k].freq, 
	    in->list[k].dFreq, in->list[k].alpha, in->list[k].dAlpha, in->list[k].delta, in->list[k].dDelta,
	    in->list[k].fdot, in->list[k].dFdot);


  DETATCHSTATUSPTR (status);
  RETURN(status);
}


/** Print Fstat vectors */
  void PrintFstatVec (LALStatus *status,
		      REAL8FrequencySeries *in,
		      FILE                 *fp,
		      PulsarDopplerParams  *thisPoint,
		      LIGOTimeGPS          refTime,
		      INT4                 stackIndex)
{
  INT4 length, k;
  REAL8 f0, deltaF, alpha, delta;
  REAL8Vector *fkdot=NULL;

  INITSTATUS( status, "PrintFstatVec", rcsid );
  ATTATCHSTATUSPTR (status);

  TRY( LALDCreateVector( status->statusPtr, &fkdot, 2), status);

  fprintf(fp, "## Fstat values from stack %d (reftime -- %d %d)\n", stackIndex, refTime.gpsSeconds, refTime.gpsNanoSeconds);

  alpha = thisPoint->Alpha;
  delta = thisPoint->Delta;
  fkdot->data[1] = thisPoint->fkdot[1];
  f0 = thisPoint->fkdot[0];

  length = in->data->length;
  deltaF = in->deltaF;

  for (k=0; k<length; k++)
    {
      fkdot->data[0] = f0 + k*deltaF;
		      
      /* propagate fkdot back to reference-time  */
      TRY ( LALExtrapolatePulsarSpins(status->statusPtr, fkdot, refTime, fkdot, thisPoint->refTime ), status );

      fprintf(fp, "%.13g %.7g %.7g %.5g %.6g\n", fkdot->data[0], alpha, delta, fkdot->data[1], in->data->data[k]);
    }

  fprintf(fp, "\n");

  TRY( LALDDestroyVector( status->statusPtr, &fkdot), status);

  DETATCHSTATUSPTR (status);
  RETURN(status);

}

/** Sets threshold on Fstat vector and fills up frequency and deltaF values in candidate list. 
    Important -- Currently this function does not get the sky-loctaion and spindown
    values for the candidates because the input does not have this information.  This 
    is to be fixed as soon as possible!
*/
void GetFstatCandidates( LALStatus *status,
			 toplist_t *list,
			 REAL8FrequencySeries *in,
			 REAL8 FstatThr,
			 REAL8 alpha,
			 REAL8 delta,
			 REAL8 fdot,
			 INT4 blockRealloc)
{
  INT4 k, length, nCandidates, maxCandidates;
  REAL8 deltaF, f0;

  INITSTATUS( status, "GetFstatCandidates", rcsid );
  ATTATCHSTATUSPTR (status);

  if (list == NULL)
    create_fstat_toplist( &list, blockRealloc);

  length = in->data->length;
  deltaF = in->deltaF;
  f0 = in->f0;


  for ( k=0; k<length; k++) {
    
	nCandidates = list->elems;
	maxCandidates = list->length;

	/* if there isn't enough memory then realloc */
	if ( nCandidates >= maxCandidates ) 
	  {
	    list->length += blockRealloc;
	    list->data = (TOPLISTLINE *)LALRealloc( list->data, 
						    maxCandidates*sizeof(TOPLISTLINE));

	  } /* end of realloc */


	if ( in->data->data[k] > FstatThr ) {

	  list->data[nCandidates].Freq = f0 + k*deltaF;
	  list->data[nCandidates].Alpha = alpha;
	  list->data[nCandidates].Delta = delta;
	  list->data[nCandidates].f1dot = fdot;
	  list->data[nCandidates].Fstat = in->data->data[k];
	  list->elems += 1;
	}
  }
  
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
  INT4 k, length, debug;
  REAL8 deltaF, f0;
  TOPLISTLINE line;

  INITSTATUS( status, "GetFstatCandidates_toplist", rcsid );
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
      
      line.Fstat = in->data->data[k];
      line.Freq = f0 + k*deltaF;
      
      debug = insert_into_fstat_toplist( list, line);

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
 
  INITSTATUS( status, "PrintHoughHistogram", rcsid );
  ATTATCHSTATUSPTR (status);

  strcpy(  filename, fnameOut);
  strcat(  filename, "histo");
  if ( !(fp=fopen(filename,"w")))
    {
      fprintf(stderr, "Unable to open file '%s' for writing\n", filename);
      exit(1);
    }

  for (i=0; i < hist->length; i++)
    fprintf(fp,"%d  %llu\n", i, hist->data[i]);
  
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
 
  INITSTATUS( status, "PrintCatalogInfo", rcsid );
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
		       CHAR *fnameChkPoint)
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
  
  INITSTATUS( status, "GetStackVelPos", rcsid );
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
  /* Since only the sfts are associated with a detector, the cleanest solution 
     (unless something better comes up) seems to be to average the velocities 
     and positions of the sfts in each stack.  Thus, for the purposes of the velocity
     and position calculation, we can view the sfts as coming from the a single 
     hypothetical ifo which is moving in a strange way */

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
