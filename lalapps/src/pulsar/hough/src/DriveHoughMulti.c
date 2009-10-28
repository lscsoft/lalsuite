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
 *  MA  02111-1307  USA2
 */

/**
 * \file DriveHough_v3.c
 * \author Badri Krishnan, Alicia Sintes, Llucia Sancho 
 * \brief Driver code for performing Hough transform search on non-demodulated
   data using SFTs from possible multiple IFOs

   Revision: $Id$
 
   History:   Created by Sintes and Krishnan July 04, 2003
              Modifications for S4 January 2006
	      Modifications for S5 November 2007

   \par Description
   
   This is the main driver for the Hough transform routines. It takes as input 
   a set of SFTs from possibly more than one IFO and outputs the number counts 
   using the Hough transform.  For a single IFO, this should be essentially equivalent 
   to DriveHough_v3.  validatehoughmultichi2

   This code just does spin-downs values. 

   \par User input

   The user inputs the following parameters:

   - Search frequency range

   - A file containing list of skypatches to search over.  For each skypatch, 
      the information is:
      - RA and dec of skypatch center.
      - Size in RA and dec.

   - Location of Directory containing the SFTs (must be v2 SFTs).

   - Interferometer (optional)

   - List of linefiles containing information about known spectral disturbanclves

   - Location of output directory and basename of output files.

   - Block size of running median for estimating psd.

   - The parameter nfSizeCylinder which determines the range of spindown parameters
      to be searched over.

   - Boolean variable for deciding if the SFTs should be inverse noise weighed.

   - Boolean variable for deciding whether amplitude modulation weights should be used.

   - Boolean variables for deciding whether the Hough maps, the statistics, list of 
      events above a threshold, and logfile should be written

   /par Output

   The output is written in several sub-directories of the specified output directory.  The
   first two items are default while the rest are written according to the user input:

   - A directory called logfiles records the user input, contents of the skypatch file 
      and cvs tags contained in the executable (if the user has required logging)

   - A directory called nstarfiles containing the loudest event for each search frequency 
      bin maximised over all sky-locations and spindown parameters.  An event is said to be
      the loudest if it has the maximum significance defined as: (number count - mean)/sigma.

   - A directory for each skypatch containing the number count statistics, the histogram, 
      the list of events, and the Hough maps 
*/

/* lalapps/hough includes */
#include "./DriveHoughColor.h"
#include "./MCInjectHoughMulti.h"

/* lalapps includes */
#include <lalapps.h>
#include <FstatToplist.h>

/* lal includes */
#include <lal/DopplerScan.h>
#include <lal/LogPrintf.h>

/* gsl includes */
#include <gsl/gsl_permutation.h>



RCSID( "$Id$");



/* globals, constants and defaults */


extern int lalDebugLevel;

/* boolean global variables for controlling output */
BOOLEAN uvar_EnableExtraInfo, uvar_EnableChi2;

/* #define EARTHEPHEMERIS "./earth05-09.dat" */
/* #define SUNEPHEMERIS "./sun05-09.dat"    */

#define EARTHEPHEMERIS "/home/badkri/lscsoft/share/lal/earth05-09.dat"
#define SUNEPHEMERIS "/home/badkri/lscsoft/share/lal/sun05-09.dat"

#define HOUGHMAXFILENAMELENGTH 512 /* maximum # of characters  of a filename */

#define DIROUT "./outMulti"   /* output directory */
#define BASENAMEOUT "HM"    /* prefix file output */

#define THRESHOLD 1.6 /* thresold for peak selection, with respect to the
                              the averaged power in the search band */
#define SKYFILE "./skypatchfile"      
#define F0 310.0   /*  frequency to build the LUT and start search */
#define FBAND 0.05   /* search frequency band  */
#define NFSIZE  21   /* n-freq. span of the cylinder, to account for spin-down search */
#define BLOCKSRNGMED 101 /* Running median window size */

#define SKYREGION "allsky" /**< default sky region to search over -- just a single point*/

#define TRUE (1==1)
#define FALSE (1==0)

#define NBLOCKSTEST 8 /* number of data blocks to do Chi2 test */

#define SPINDOWNJUMP 1 /* "Jump" to the next spin-down being analyzed (to avoid doing them all) */

/* ****************************************
 * Structure, HoughParamsTest, typedef
 */

typedef struct tagHoughParamsTest{
    UINT4  length;            /* number p of blocks to split the data into */
    UINT4  *numberSFTp;       /* Ni SFTs in each data block */
    REAL8  *sumWeight;        /* Pointer to the sumWeight of each block of data */
    REAL8  *sumWeightSquare;  /* Pointer to the sumWeightSquare of each block of data */
}HoughParamsTest; 


 typedef struct tagUCHARPeakGramVector{
   UINT4             length; /**< number of elements */
   UCHARPeakGram     *upg;    /**< expanded Peakgrams */
 } UCHARPeakGramVector;

/******************************************/  

/* local function prototype */

void SplitSFTs(LALStatus *status, REAL8Vector *weightsV, HoughParamsTest *chi2Params);

void ComputeFoft_NM(LALStatus *status, REAL8Vector *foft, HoughTemplate *pulsarTemplate, REAL8Vector *timeDiffV, REAL8Cart3CoorVector *velV);

void ComputeandPrintChi2 ( LALStatus *status, toplist_t *tl, REAL8Vector *timeDiffV, REAL8Cart3CoorVector *velV, INT4 p, REAL8 alphaPeak, MultiDetectorStateSeries *mdetStates, REAL8Vector *weightsNoise, UCHARPeakGramVector *upgV);

void GetPeakGramFromMultSFTVector_NondestroyPg1(LALStatus *status, HOUGHPeakGramVector *out, UCHARPeakGramVector *upgV, MultiSFTVector *in, REAL8 thr) ; 

void PrintLogFile (LALStatus *status, CHAR *dir, CHAR *basename, CHAR *skyfile, LALStringVector *linefiles, CHAR *executable );

int CreateSkypatchDirs(CHAR *filestats, CHAR *base, INT4 index);

int PrintHistogram(UINT8Vector *hist, CHAR *fnameOut, REAL8 minSignificance, REAL8 maxSignificance);

int PrintHmap2m_file(HOUGHMapTotal *ht, CHAR *fnameOut, INT4 iHmap);

int PrintHmap2file(HOUGHMapTotal *ht, CHAR *fnameOut, INT4 iHmap);

int OpenExtraInfoFiles(CHAR *fileMaps, FILE **fp1_ptr, CHAR *filehisto, CHAR *dirname, CHAR *basename, INT4 index);

int PrintExtraInfo(CHAR *fileMaps, FILE **fp1_ptr, INT4 iHmap, HOUGHMapTotal *ht, REAL8UnitPolarCoor *sourceLocation, HoughStats *stats, INT8 fBinSearch, REAL8 deltaF);

void ReadTimeStampsFile (LALStatus *status, LIGOTimeGPSVector *ts, CHAR *filename);

void LALHoughHistogramSignificance(LALStatus *status, UINT8Vector *out, HOUGHMapTotal *in, 
				   REAL8 mean, REAL8 sigma, REAL8 minSignificance, REAL8 maxSignificance);

void GetSFTVelTime(LALStatus *status, REAL8Cart3CoorVector *velV, LIGOTimeGPSVector *timeV, MultiDetectorStateSeries *in);

void GetSFTNoiseWeights(LALStatus *status, REAL8Vector *out, MultiNoiseWeights  *in);

void GetPeakGramFromMultSFTVector(LALStatus *status, HOUGHPeakGramVector *out, MultiSFTVector *in, REAL8 thr);

void SetUpSkyPatches(LALStatus *status, HoughSkyPatchesInfo *out, CHAR *skyFileName, CHAR *skyRegion, REAL8 dAlpha, REAL8 dDelta, INT4 numSkyPartitions, INT4 partitionIndex);

void GetAMWeights(LALStatus *status, REAL8Vector *out, MultiDetectorStateSeries *mdetStates, REAL8 alpha, REAL8 delta);

void SelectBestStuff(LALStatus *status, BestVariables *out, BestVariables  *in,	UINT4 mObsCohBest);

void DuplicateBestStuff(LALStatus *status, BestVariables *out, BestVariables *in);

void GetToplistFromHoughmap(LALStatus *status, toplist_t *list, HOUGHMapTotal *ht, HOUGHPatchGrid *patch, HOUGHDemodPar *parDem, REAL8 mean, REAL8 sigma);


void LALHOUGHCreateLUTVector(LALStatus           *status,
			     HOUGHptfLUTVector   *lutV,
			     UINT4               length);

void LALHOUGHCreateLUTs(LALStatus           *status,
			HOUGHptfLUTVector   *lutV,
			UINT2               maxNBins, 
			UINT2               maxNBorders,
			UINT2               ySide);

void LALHOUGHCreatePHMDVS(LALStatus           *status,
			  PHMDVectorSequence  *phmdVS,
			  UINT4               length,
			  UINT4               nfSize);

void LALHOUGHCreatePHMDs(LALStatus           *status,
			 PHMDVectorSequence  *phmdVS,
			 UINT2               maxNBins, 
			 UINT2               maxNBorders,
			 UINT2               ySide);

void LALHOUGHDestroyPHMDs(LALStatus           *status,
			  PHMDVectorSequence  *phmdVS);

void LALHOUGHCreateHT(LALStatus             *status,
		      HOUGHMapTotal         *ht,
		      UINT2                 xSide,
		      UINT2                 ySide);

void LALHOUGHCreateFreqIndVector(LALStatus                 *status,
				 UINT8FrequencyIndexVector *freqInd,
				 UINT4                     length,
				 REAL8                     deltaF);


void LALHOUGHDestroyLUTs(LALStatus           *status,
			 HOUGHptfLUTVector   *lut);

/******************************************/

int main(int argc, char *argv[]){

  /* LALStatus pointer */
  static LALStatus  status;  
  
  /* time and velocity  */
  LIGOTimeGPSVector    *timeV=NULL;
  REAL8Vector  *timeDiffV=NULL;
  static REAL8Cart3CoorVector velV;


  LIGOTimeGPS firstTimeStamp, lastTimeStamp;
  REAL8 tObs;

  /* standard pulsar sft types */ 
  MultiSFTVector *inputSFTs = NULL;
  UINT4 binsSFT, sftFminBin;
  UINT4 numSearchBins;
  
  /* information about all the ifos */
  MultiDetectorStateSeries *mdetStates = NULL;
  UINT4 numifo;

  /* vector of weights */
  REAL8Vector *weightsV=NULL, *weightsNoise=NULL;

  /* ephemeris */
  EphemerisData    *edat=NULL;

  /* struct containing subset of weights, timediff, velocity and peakgrams */
  static BestVariables best;

  /* hough structures */
  static HOUGHptfLUTVector   lutV; /* the Look Up Table vector*/
  static HOUGHPeakGramVector pgV;  /* vector of peakgrams */
  static UCHARPeakGramVector upgV;  /* vector of expanded peakgrams */
  static PHMDVectorSequence  phmdVS;  /* the partial Hough map derivatives */
  static UINT8FrequencyIndexVector freqInd; /* for trajectory in time-freq plane */
  static HOUGHResolutionPar parRes;   /* patch grid information */
  static HOUGHPatchGrid  patch;   /* Patch description */ 
  static HOUGHParamPLUT  parLut;  /* parameters needed to build lut  */
  static HOUGHDemodPar   parDem;  /* demodulation parameters or  */
  static HOUGHSizePar    parSize; 
  static HOUGHMapTotal   ht;   /* the total Hough map */
  static UINT8Vector     *hist; /* histogram of number counts for a single map */
  static UINT8Vector     *histTotal; /* number count histogram for all maps */
  static HoughStats      stats;  /* statistical information about a Hough map */

  /* skypatch info */
  REAL8  *skyAlpha=NULL, *skyDelta=NULL, *skySizeAlpha=NULL, *skySizeDelta=NULL; 
  INT4   nSkyPatches, skyCounter=0; 
  static HoughSkyPatchesInfo skyInfo;

  /* output filenames and filepointers */
  CHAR   filehisto[ HOUGHMAXFILENAMELENGTH ]; 
  /* CHAR   filestats[ HOUGHMAXFILENAMELENGTH ];*/ 
  CHAR   fileMaps[ HOUGHMAXFILENAMELENGTH ];
  CHAR   fileSigma[ HOUGHMAXFILENAMELENGTH ];
  FILE   *fp1 = NULL;
  FILE   *fpSigma = NULL;

  /* the maximum number count */
  static HoughSignificantEventVector nStarEventVec;

  /* miscellaneous */
  INT4   iHmap, nSpin1Max ;
  UINT4  mObsCoh, mObsCohBest;
  INT8   f0Bin, fLastBin, fBin;
  REAL8  alpha, delta, timeBase, deltaF, f1jump;
  REAL8  patchSizeX, patchSizeY;
  REAL8  alphaPeak, minSignificance, maxSignificance;
  REAL8  meanN, sigmaN;
  UINT2  xSide, ySide;
  UINT2  maxNBins, maxNBorders;

  /* output toplist candidate structure */
  toplist_t *toplist=NULL;

  /* sft constraint variables */
  LIGOTimeGPS startTimeGPS, endTimeGPS;
  LIGOTimeGPSVector inputTimeStampsVector;

  /* user input variables */
  BOOLEAN  uvar_help, uvar_weighAM, uvar_weighNoise, uvar_printLog;
  INT4     uvar_blocksRngMed, uvar_nfSizeCylinder, uvar_maxBinsClean, uvar_binsHisto;  
  INT4     uvar_keepBestSFTs=1;
  INT4     uvar_numCand;
  REAL8    uvar_startTime, uvar_endTime;
  REAL8    uvar_pixelFactor;
  REAL8    uvar_f0, uvar_peakThreshold, uvar_freqBand;
  REAL8    uvar_dAlpha, uvar_dDelta; /* resolution for isotropic sky-grid */
  CHAR     *uvar_earthEphemeris=NULL;
  CHAR     *uvar_sunEphemeris=NULL;
  CHAR     *uvar_sftData=NULL;
  CHAR     *uvar_dirnameOut=NULL;
  CHAR     *uvar_fbasenameOut=NULL;
  CHAR     *uvar_skyfile=NULL;
  CHAR     *uvar_timeStampsFile=NULL;
  CHAR     *uvar_skyRegion=NULL;
  LALStringVector *uvar_linefiles=NULL;

  INT4     uvar_chiSqBins;

  INT4     uvar_spindownJump;

  INT4 uvar_numSkyPartitions = 0;
  INT4 uvar_partitionIndex = 0;

  LIGOTimeGPS refTimeGPS; /* reference time */
  REAL8    uvar_refTime;

  REAL8 uvar_deltaF1dot; 

  /* Set up the default parameters */

  /* LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;
  
  lalDebugLevel = 0;  /* LALDebugLevel must be called before anything else */
  LAL_CALL( LALGetDebugLevel( &status, argc, argv, 'd'), &status);
  
  uvar_help = FALSE;
  uvar_weighAM = TRUE;
  uvar_weighNoise = TRUE;
  uvar_printLog = FALSE;
  uvar_blocksRngMed = BLOCKSRNGMED;
  uvar_nfSizeCylinder = NFSIZE;
  uvar_f0 = F0;
  uvar_freqBand = FBAND;
  uvar_peakThreshold = THRESHOLD;
  uvar_maxBinsClean = 100;
  uvar_binsHisto = 1000;
  uvar_pixelFactor = PIXELFACTOR;
  uvar_dAlpha = 0.2;
  uvar_dDelta = 0.2;
  uvar_numCand=1;
  uvar_EnableExtraInfo=FALSE;
  uvar_EnableChi2=FALSE;
  uvar_chiSqBins = NBLOCKSTEST;
  uvar_spindownJump = SPINDOWNJUMP;


  uvar_earthEphemeris = (CHAR *)LALCalloc( HOUGHMAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_earthEphemeris,EARTHEPHEMERIS);

  uvar_sunEphemeris = (CHAR *)LALCalloc( HOUGHMAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_sunEphemeris,SUNEPHEMERIS);

  uvar_dirnameOut = (CHAR *)LALCalloc( HOUGHMAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_dirnameOut,DIROUT);

  uvar_fbasenameOut = (CHAR *)LALCalloc( HOUGHMAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_fbasenameOut,BASENAMEOUT);

  uvar_skyfile = (CHAR *)LALCalloc( HOUGHMAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_skyfile,SKYFILE);

  /* register user input variables */
  LAL_CALL( LALRegisterBOOLUserVar( &status, "help",             'h',  UVAR_HELP,     "Print this message", &uvar_help), &status);  
  LAL_CALL( LALRegisterREALUserVar( &status, "f0",               'f',  UVAR_OPTIONAL, "Start search frequency", &uvar_f0), &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "freqBand",         'b',  UVAR_OPTIONAL, "Search frequency band", &uvar_freqBand), &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "startTime",         0,  UVAR_OPTIONAL, "GPS start time of observation", &uvar_startTime), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "endTime",         0,  UVAR_OPTIONAL, "GPS end time of observation", &uvar_endTime), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "timeStampsFile",  0,  UVAR_OPTIONAL, "Input time-stamps file", &uvar_timeStampsFile),   &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "skyRegion",       0,  UVAR_OPTIONAL, "sky-region polygon (or 'allsky')", &uvar_skyRegion), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dAlpha",          0,  UVAR_OPTIONAL, "Resolution for flat or isotropic coarse grid (rad)", &uvar_dAlpha), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dDelta",          0,  UVAR_OPTIONAL, "Resolution for flat or isotropic coarse grid (rad)", &uvar_dDelta), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "skyfile",         0,  UVAR_OPTIONAL, "Alternative: input skypatch file", &uvar_skyfile),  &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "peakThreshold",   0,  UVAR_OPTIONAL, "Peak selection threshold", &uvar_peakThreshold),   &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "weighAM",         0,  UVAR_OPTIONAL, "Use amplitude modulation weights", &uvar_weighAM),  &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "weighNoise",      0,  UVAR_OPTIONAL, "Use SFT noise weights", &uvar_weighNoise), &status);  
  LAL_CALL( LALRegisterINTUserVar(    &status, "keepBestSFTs",    0,  UVAR_OPTIONAL, "Number of best SFTs to use (default--keep all)", &uvar_keepBestSFTs),  &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printLog",        0,  UVAR_OPTIONAL, "Print Log file", &uvar_printLog), &status);  
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "earthEphemeris", 'E', UVAR_OPTIONAL, "Earth Ephemeris file",  &uvar_earthEphemeris),  &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sunEphemeris",   'S', UVAR_OPTIONAL, "Sun Ephemeris file", &uvar_sunEphemeris), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftData",        'D', UVAR_REQUIRED, "SFT filename pattern", &uvar_sftData), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "dirnameOut",     'o', UVAR_OPTIONAL, "Output directory", &uvar_dirnameOut), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "fbasenameOut",    0,  UVAR_OPTIONAL, "Output file basename", &uvar_fbasenameOut), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "binsHisto",       0,  UVAR_OPTIONAL, "No. of bins for histogram", &uvar_binsHisto),  &status);
  LAL_CALL( LALRegisterLISTUserVar(   &status, "linefiles",       0,  UVAR_OPTIONAL, "Comma separated List of linefiles (filenames must contain IFO name)", &uvar_linefiles), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "nfSizeCylinder",  0,  UVAR_OPTIONAL, "Size of cylinder of PHMDs", &uvar_nfSizeCylinder),  &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "pixelFactor",    'p', UVAR_OPTIONAL, "sky resolution=1/v*pixelFactor*f*Tcoh", &uvar_pixelFactor), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "numCand",         0,  UVAR_OPTIONAL, "No. of toplist candidates", &uvar_numCand), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printExtraInfo",  0,  UVAR_OPTIONAL, "Print HoughMaps, HoughStatistics, expected number count stdev", &uvar_EnableExtraInfo), &status); 
  LAL_CALL( LALRegisterINTUserVar(    &status, "chiSqBins",       0,  UVAR_OPTIONAL, "Number of chi-square bins for veto tests",  &uvar_chiSqBins), &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "enableChi2",      0,  UVAR_OPTIONAL, "Print Chi2 value for each element in the Toplist", &uvar_EnableChi2), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "spindownJump",       0,  UVAR_OPTIONAL, "Jump to the next spin-down being analyzed (to avoid doing them all)",  &uvar_spindownJump), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "numSkyPartitions",0,UVAR_OPTIONAL, "Number of (equi-)partitions to split skygrid into", &uvar_numSkyPartitions), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "partitionIndex",0,UVAR_OPTIONAL, "Index [0,numSkyPartitions-1] of sky-partition to generate", &uvar_partitionIndex), &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "refTime",         0,  UVAR_OPTIONAL, "GPS reference time of observation", &uvar_refTime), &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "deltaF1dot",         0,  UVAR_OPTIONAL, "(Step size for f1dot)*Tcoh [Default: 1/Tobs]", &uvar_deltaF1dot), &status);

  /* developer input variables */
  LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed",    0, UVAR_DEVELOPER, "Running Median block size", &uvar_blocksRngMed), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "maxBinsClean",    0, UVAR_DEVELOPER, "Maximum number of bins in cleaning", &uvar_maxBinsClean), &status);


  /* read all command line variables */
  LAL_CALL( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 

  /* very basic consistency checks on user input */
  if ( uvar_f0 < 0 ) {
    LogPrintf(LOG_CRITICAL, "start frequency must be positive\n");
    exit(1);
  }

  if ( uvar_freqBand < 0 ) {
    LogPrintf(LOG_CRITICAL, "search frequency band must be positive\n");
    exit(1);
  }
 
  if ( uvar_peakThreshold < 0 ) {
    LogPrintf(LOG_CRITICAL, "peak selection threshold must be positive\n");
    exit(1);
  }

  /* probability of peak selection */
  alphaPeak = exp( -uvar_peakThreshold);


  if ( uvar_binsHisto < 1 ) {
    LogPrintf(LOG_CRITICAL, "binsHisto must be at least 1\n");
    exit(1);
  }

  if ( uvar_keepBestSFTs < 1 ) {
    LogPrintf(LOG_CRITICAL, "must keep at least 1 SFT\n");
    exit(1);
  }

  /* write log file with command line arguments, cvs tags, and contents of skypatch file */
  if ( uvar_printLog ) {
    LAL_CALL( PrintLogFile( &status, uvar_dirnameOut, uvar_fbasenameOut, uvar_skyfile, uvar_linefiles, argv[0]), &status);
  }


  /***** start main calculations *****/

  LogPrintf (LOG_NORMAL, "Setting up sky-patches...");
  /* set up skypatches */
  LAL_CALL( SetUpSkyPatches( &status, &skyInfo, uvar_skyfile, uvar_skyRegion, uvar_dAlpha, uvar_dDelta, uvar_numSkyPartitions, uvar_partitionIndex), &status);
  nSkyPatches = skyInfo.numSkyPatches;
  skyAlpha = skyInfo.alpha;
  skyDelta = skyInfo.delta;
  skySizeAlpha = skyInfo.alphaSize;
  skySizeDelta = skyInfo.deltaSize;
  LogPrintfVerbatim (LOG_NORMAL, "done\n");

  /* set up toplist */
  /* create toplist -- semiCohToplist has the same structure 
     as a fstat candidate, so treat it as a fstat candidate */
  if ( create_fstat_toplist(&toplist, uvar_numCand) != 0) {
    LogPrintf(LOG_CRITICAL,"Unable to create toplist\n");
  }


  LogPrintf (LOG_NORMAL, "Reading SFTs...");
  /* read sft Files and set up weights */
  {
    /* new SFT I/O data types */
    SFTCatalog *catalog = NULL;
    static SFTConstraints constraints;

    REAL8 doppWings, fmin, fmax;
    INT4 k;

    /* set detector constraint */
    constraints.detector = NULL;

    if ( LALUserVarWasSet( &uvar_startTime ) ) {
      XLALGPSSetREAL8(&startTimeGPS, uvar_startTime);
      constraints.startTime = &startTimeGPS;
    }

    if ( LALUserVarWasSet( &uvar_endTime ) ) {
      XLALGPSSetREAL8(&endTimeGPS, uvar_endTime);
      constraints.endTime = &endTimeGPS;
    }

    if ( LALUserVarWasSet( &uvar_timeStampsFile ) ) {
      LAL_CALL ( ReadTimeStampsFile ( &status, &inputTimeStampsVector, uvar_timeStampsFile), &status);
      constraints.timestamps = &inputTimeStampsVector;
    }
    
    /* get sft catalog */
    LAL_CALL( LALSFTdataFind( &status, &catalog, uvar_sftData, &constraints), &status);
    if ( (catalog == NULL) || (catalog->length == 0) ) {
      LogPrintf (LOG_CRITICAL,"Unable to match any SFTs with pattern '%s'\n", uvar_sftData );
      exit(1);
    }

    /* now we can free the inputTimeStampsVector */
    if ( LALUserVarWasSet( &uvar_timeStampsFile ) ) {
      LALFree( inputTimeStampsVector.data );
    }



    /* first some sft parameters */
    deltaF = catalog->data[0].header.deltaF;  /* frequency resolution */
    timeBase= 1.0/deltaF; /* coherent integration time */
    f0Bin = floor( uvar_f0 * timeBase + 0.5); /* initial search frequency */
    numSearchBins =  uvar_freqBand * timeBase; /* total number of search bins - 1 */
    fLastBin = f0Bin + numSearchBins;   /* final frequency bin to be analyzed */
    

    /* read sft files making sure to add extra bins for running median */
    /* add wings for Doppler modulation and running median block size*/
    doppWings = (uvar_f0 + uvar_freqBand) * VTOT;    
    fmin = uvar_f0 - doppWings - (uvar_blocksRngMed + uvar_nfSizeCylinder) * deltaF;
    fmax = uvar_f0 + uvar_freqBand + doppWings + (uvar_blocksRngMed + uvar_nfSizeCylinder) * deltaF;

    /* read the sfts */
    LAL_CALL( LALLoadMultiSFTs ( &status, &inputSFTs, catalog, fmin, fmax), &status);
    numifo = inputSFTs->length;    

    /* find number of sfts */
    /* loop over ifos and calculate number of sfts */
    /* note that we can't use the catalog to determine the number of SFTs
       because SFTs might be segmented in frequency */
    mObsCoh = 0; /* initialization */
    for (k = 0; k < (INT4)numifo; k++ ) {
      mObsCoh += inputSFTs->data[k]->length;
    } 
    
    /* set number of SFTs to be kept */
    /* currently mobscohbest is set equal to mobscoh if no weights are used 
       -- this probably will be changed in the future */
    if ( LALUserVarWasSet( &uvar_keepBestSFTs ) && (uvar_weighNoise||uvar_weighAM)) {
      mObsCohBest = uvar_keepBestSFTs;

      /* set to mobscoh if it is more than number of sfts */
      if ( mObsCohBest > mObsCoh )
	mObsCohBest = mObsCoh;
      }  
    else 
      mObsCohBest = mObsCoh;  


    /* catalog is ordered in time so we can get start, end time and tObs*/
    firstTimeStamp = catalog->data[0].header.epoch;
    lastTimeStamp = catalog->data[catalog->length - 1].header.epoch;

    if ( LALUserVarWasSet( &uvar_refTime ) ) 
    {
      XLALGPSSetREAL8(&refTimeGPS, uvar_refTime);
      tObs = XLALGPSDiff( &lastTimeStamp, &refTimeGPS ) + timeBase;
    }
    else
    {
      tObs = XLALGPSDiff( &lastTimeStamp, &firstTimeStamp ) + timeBase;  
    }

   


    /* clean sfts if required */
    if ( LALUserVarWasSet( &uvar_linefiles ) )
      {

	RandomParams *randPar=NULL;
	FILE *fpRand=NULL;
	INT4 seed, ranCount;  

	LogPrintfVerbatim (LOG_NORMAL, "...cleaning SFTs...");
	if ( (fpRand = fopen("/dev/urandom", "r")) == NULL ) {
	  LogPrintf(LOG_CRITICAL,"error in opening /dev/urandom" ); 
	  exit(1);
	} 

	if ( (ranCount = fread(&seed, sizeof(seed), 1, fpRand)) != 1 ) {
	  LogPrintf(LOG_CRITICAL,"error in getting random seed" );
	  exit(1);
	}

	LAL_CALL ( LALCreateRandomParams (&status, &randPar, seed), &status );

	LAL_CALL( LALRemoveKnownLinesInMultiSFTVector ( &status, inputSFTs, uvar_maxBinsClean, 
							uvar_blocksRngMed, uvar_linefiles, randPar), &status);

	LAL_CALL ( LALDestroyRandomParams (&status, &randPar), &status);
	fclose(fpRand);
      } /* end cleaning */


    /* SFT info -- assume all SFTs have same length */
    binsSFT = inputSFTs->data[0]->data->data->length;
    sftFminBin = (INT4) floor(inputSFTs->data[0]->data[0].f0 * timeBase + 0.5);


    LAL_CALL( LALDestroySFTCatalog( &status, &catalog ), &status);  	

  } /* end of sft reading block */
  LogPrintfVerbatim (LOG_NORMAL, "done\n");


  /** some memory allocations */

  /* allocate memory for velocity vector */
  velV.length = mObsCoh;
  velV.data = NULL;
  velV.data = (REAL8Cart3Coor *)LALCalloc(1, mObsCoh*sizeof(REAL8Cart3Coor));
  
  /* allocate memory for timestamps and timediff vectors */
  timeV = XLALCreateTimestampVector (mObsCoh);
  timeDiffV = XLALCreateREAL8Vector( mObsCoh);
  
  /* allocate and initialize noise and AMweights vectors */
  weightsV = XLALCreateREAL8Vector( mObsCoh);
  weightsNoise = XLALCreateREAL8Vector( mObsCoh);    
  LAL_CALL( LALHOUGHInitializeWeights( &status, weightsNoise), &status);
  LAL_CALL( LALHOUGHInitializeWeights( &status, weightsV), &status);
  
 
  LogPrintf (LOG_NORMAL, "Setting up weights...");  
  /* get detector velocities weights vector, and timestamps */
  { 
    MultiNoiseWeights *multweight = NULL;    
    MultiPSDVector *multPSD = NULL;  
    UINT4 j;

    /*  get ephemeris  */
    edat = (EphemerisData *)LALCalloc(1, sizeof(EphemerisData));
    (*edat).ephiles.earthEphemeris = uvar_earthEphemeris;
    (*edat).ephiles.sunEphemeris = uvar_sunEphemeris;
    LAL_CALL( LALInitBarycenter( &status, edat), &status);


    /* normalize sfts */
    LAL_CALL( LALNormalizeMultiSFTVect (&status, &multPSD, inputSFTs, uvar_blocksRngMed), &status);
   
    /* compute multi noise weights */
    if ( uvar_weighNoise ) {
      LAL_CALL ( LALComputeMultiNoiseWeights ( &status, &multweight, multPSD, uvar_blocksRngMed, 0), &status);
    }
    
    /* we are now done with the psd */
    LAL_CALL ( LALDestroyMultiPSDVector  ( &status, &multPSD), &status);

    /* get information about all detectors including velocity and timestamps */
    /* note that this function returns the velocity at the 
       mid-time of the SFTs -- should not make any difference */
    LAL_CALL ( LALGetMultiDetectorStates ( &status, &mdetStates, inputSFTs, edat), &status);

    LAL_CALL ( GetSFTVelTime( &status, &velV, timeV, mdetStates), &status);

    /* copy the noise-weights vector if required*/
    if ( uvar_weighNoise ) {

      LAL_CALL ( GetSFTNoiseWeights( &status, weightsNoise, multweight), &status);

      LAL_CALL ( LALDestroyMultiNoiseWeights ( &status, &multweight), &status);
    }

    /* compute the time difference relative to startTime for all SFTs */
    if ( LALUserVarWasSet( &uvar_refTime ) ) 
    { 
      for(j = 0; j < mObsCoh; j++)
      timeDiffV->data[j] = XLALGPSDiff( timeV->data + j, &refTimeGPS );
    }
    else
    {
      for(j = 0; j < mObsCoh; j++)
      timeDiffV->data[j] = XLALGPSDiff( timeV->data + j, &firstTimeStamp );
    }

  } /* end block for weights, velocity and time */
  LogPrintfVerbatim (LOG_NORMAL, "done\n");    


  LogPrintf (LOG_NORMAL, "Generating peakgrams...");     
  /* generating peakgrams  */  
  pgV.length = mObsCoh;
  pgV.pg = NULL;
  pgV.pg = (HOUGHPeakGram *)LALCalloc(1,mObsCoh*sizeof(HOUGHPeakGram));
  
  if (uvar_EnableChi2)
  {   
    upgV.length = mObsCoh;
    upgV.upg = NULL;
    upgV.upg = (UCHARPeakGram *)LALCalloc(1,mObsCoh*sizeof(UCHARPeakGram));
    
    LAL_CALL( GetPeakGramFromMultSFTVector_NondestroyPg1(&status, &pgV, &upgV, inputSFTs, uvar_peakThreshold),&status);
  }  
  else
  {
    LAL_CALL( GetPeakGramFromMultSFTVector( &status, &pgV, inputSFTs, uvar_peakThreshold), &status);
  }


  /* we are done with the sfts and ucharpeakgram now */
  LAL_CALL (LALDestroyMultiSFTVector(&status, &inputSFTs), &status );
  LogPrintfVerbatim (LOG_NORMAL, "done\n");

  /* if we want to print expected sigma for each skypatch */
  if ( uvar_EnableExtraInfo ) 
    {
      strcpy ( fileSigma, uvar_dirnameOut);
      strcat ( fileSigma, "/");
      strcat ( fileSigma, uvar_fbasenameOut);
      strcat ( fileSigma, "sigma");
      
      if ( (fpSigma = fopen(fileSigma,"w")) == NULL)
	{
	  LogPrintf(LOG_CRITICAL,"Unable to find file %s for writing\n", fileSigma);
	  return DRIVEHOUGHCOLOR_EFILE;
	}
    } /* end if( uvar_EnableExtraInfo) */


  /* min and max values significance that are possible */
  minSignificance = -sqrt(mObsCohBest * alphaPeak/(1-alphaPeak));
  maxSignificance = sqrt(mObsCohBest * (1-alphaPeak)/alphaPeak);
      

  LogPrintf (LOG_NORMAL, "Starting loop over skypatches...");
  /* loop over sky patches -- main Hough calculations */
  for (skyCounter = 0; skyCounter < nSkyPatches; skyCounter++)
    {
      UINT4 k;
      REAL8 sumWeightSquare;
      /*     REAL8  meanN, sigmaN;*/
      BestVariables temp;

      LogPrintfVerbatim (LOG_NORMAL, "%d/%d,",skyCounter, nSkyPatches);

      /* set sky positions and skypatch sizes */
      alpha = skyAlpha[skyCounter];
      delta = skyDelta[skyCounter];
      patchSizeX = skySizeDelta[skyCounter];
      patchSizeY = skySizeAlpha[skyCounter];

      /* copy noise weights if required */
      if ( uvar_weighNoise )
	memcpy(weightsV->data, weightsNoise->data, mObsCoh * sizeof(REAL8));

      /* calculate amplitude modulation weights if required */
      if (uvar_weighAM) {
	LAL_CALL( GetAMWeights( &status, weightsV, mdetStates, alpha, delta), &status);
      }
      
      /* sort weights vector to get the best sfts */
      temp.length = mObsCoh;
      temp.weightsV = weightsV;
      temp.timeDiffV = timeDiffV;
      temp.velV = &velV;
      temp.pgV = &pgV;

      if ( uvar_weighAM || uvar_weighNoise ) {
	LAL_CALL( SelectBestStuff( &status, &best, &temp, mObsCohBest), &status);	
      }
      else {
	LAL_CALL( DuplicateBestStuff( &status, &best, &temp), &status);	
      }

      /* Normalize the Best SFTs weights */
      LAL_CALL( LALHOUGHNormalizeWeights( &status, best.weightsV), &status);
      
      /* calculate the sum of the weights squared */
      sumWeightSquare = 0.0;
      for ( k = 0; k < mObsCohBest; k++)
	sumWeightSquare += best.weightsV->data[k] * best.weightsV->data[k];

      /* probability of selecting a peak expected mean and standard deviation for noise only */
      meanN = mObsCohBest * alphaPeak; 
      sigmaN = sqrt(sumWeightSquare * alphaPeak * (1.0 - alphaPeak));
      

      if ( uvar_EnableExtraInfo )
      {
	  fprintf(fpSigma, "%f\n", sigmaN);
	  if ( OpenExtraInfoFiles( fileMaps, &fp1, filehisto, uvar_dirnameOut, uvar_fbasenameOut, skyCounter ))
	       return DRIVEHOUGHCOLOR_EFILE;
      }

      /****  general parameter settings and 1st memory allocation ****/      
      
      LAL_CALL ( LALHOUGHCreateLUTVector( &status, &lutV, mObsCohBest), &status);
      
      LAL_CALL( LALHOUGHCreatePHMDVS( &status, &phmdVS, mObsCohBest, uvar_nfSizeCylinder), &status);
      phmdVS.deltaF  = deltaF;
      
      LAL_CALL( LALHOUGHCreateFreqIndVector( &status, &freqInd, mObsCohBest, deltaF), &status);

      /* allocating histogram of the number-counts in the Hough maps */
      if ( uvar_EnableExtraInfo ) {
	UINT4 k;	
	hist = XLALCreateUINT8Vector (uvar_binsHisto);
	histTotal = XLALCreateUINT8Vector (uvar_binsHisto);	

	/* initialize to 0 */
	/*  memset(histTotal->data, histTotal->length*sizeof(histTotal->data[0]), 0); */
	for (k = 0; k < histTotal->length; k++) {
	  histTotal->data[k] = 0;
	  hist->data[k] = 0;
	}
      }

      /* set demodulation pars for non-demodulated data (SFT input)*/
      parDem.deltaF = deltaF;
      parDem.skyPatch.alpha = alpha;
      parDem.skyPatch.delta = delta;
      parDem.timeDiff = 0.0;
      parDem.spin.length = 0;
      parDem.spin.data = NULL;
      parDem.positC.x = 0.0;
      parDem.positC.y = 0.0;
      parDem.positC.z = 0.0;
      
      /* sky-resolution parameters **/
      parRes.deltaF = deltaF;
      parRes.patchSkySizeX  = patchSizeX;
      parRes.patchSkySizeY  = patchSizeY;
      parRes.pixelFactor = uvar_pixelFactor;
      parRes.pixErr = PIXERR;
      parRes.linErr = LINERR;
      parRes.vTotC = VTOT;


      
      fBin= f0Bin;
      iHmap = 0;
      
      /* ***** for spin-down case ****/
      nSpin1Max = uvar_nfSizeCylinder - 1 ;
      /* nSpin1Max = floor(uvar_nfSizeCylinder/2.0) ;*/


      if ( LALUserVarWasSet( &uvar_deltaF1dot ) ) 
      {
	f1jump = uvar_deltaF1dot * uvar_spindownJump;
      }
      else
      {
	f1jump = 1./tObs * uvar_spindownJump;	
      }

      
      /* start of main loop over search frequency bins */
      /********** starting the search from f0Bin to fLastBin.
		  Note one set LUT might not cover all the interval.
		  This is taken into account *******************/

      while( fBin <= fLastBin){
	INT8 fBinSearch, fBinSearchMax;
	UINT4 j; 
	REAL8UnitPolarCoor sourceLocation;
	
	
	parRes.f0Bin =  fBin;      
	LAL_CALL( LALHOUGHComputeNDSizePar( &status, &parSize, &parRes ),  &status );
	xSide = parSize.xSide;
	ySide = parSize.ySide;
	maxNBins = parSize.maxNBins;
	maxNBorders = parSize.maxNBorders;
	
	/* *******************create patch grid at fBin ****************  */
	patch.xSide = xSide;
	patch.ySide = ySide;
	patch.xCoor = NULL;
	patch.yCoor = NULL;
	patch.xCoor = (REAL8 *)LALCalloc(1,xSide*sizeof(REAL8));
	patch.yCoor = (REAL8 *)LALCalloc(1,ySide*sizeof(REAL8));
	LAL_CALL( LALHOUGHFillPatchGrid( &status, &patch, &parSize ), &status );
	
	/*************** other memory allocation and settings************ */

	LAL_CALL( LALHOUGHCreateLUTs( &status, &lutV, maxNBins, maxNBorders, ySide), &status);
	
	LAL_CALL( LALHOUGHCreatePHMDs( &status, &phmdVS, maxNBins, maxNBorders, ySide), &status);


	/* ************* create all the LUTs at fBin ********************  */  
	for (j = 0; j < mObsCohBest; ++j){  /* create all the LUTs */
	  parDem.veloC.x = best.velV->data[j].x;
	  parDem.veloC.y = best.velV->data[j].y;
	  parDem.veloC.z = best.velV->data[j].z;      
	  /* calculate parameters needed for buiding the LUT */
	  LAL_CALL( LALNDHOUGHParamPLUT( &status, &parLut, &parSize, &parDem),&status );
	  /* build the LUT */
	  LAL_CALL( LALHOUGHConstructPLUT( &status, &(lutV.lut[j]), &patch, &parLut ),
	       &status );
	}
        
	/************* build the set of  PHMD centered around fBin***********/     
	phmdVS.fBinMin = fBin - uvar_nfSizeCylinder +1 ;
	/*phmdVS.fBinMin = fBin - floor( uvar_nfSizeCylinder/2.) ;*/
	
	LAL_CALL( LALHOUGHConstructSpacePHMD(&status, &phmdVS, best.pgV, &lutV), &status );
	if (uvar_weighAM || uvar_weighNoise) {
	  LAL_CALL( LALHOUGHWeighSpacePHMD(&status, &phmdVS, best.weightsV), &status);
	}
	
	/* ************ initializing the Total Hough map space *********** */   
	
	LAL_CALL( LALHOUGHCreateHT( &status, &ht, xSide, ySide), &status);
	ht.mObsCoh = mObsCohBest;
	ht.deltaF = deltaF;


	/*  Search frequency interval possible using the same LUTs */
	fBinSearch = fBin;
	fBinSearchMax = fBin + parSize.nFreqValid - 1 ;
	/*fBinSearchMax = fBin + parSize.nFreqValid - 1 - floor((uvar_nfSizeCylinder/2. ) ;*/
	

	/* Study all possible frequencies with one set of LUT */	

	while ( (fBinSearch <= fLastBin) && (fBinSearch < fBinSearchMax) ) 
	  {
	    
	    /**** study 1 spin-down. at  fBinSearch ****/

	    INT4   n;
	    REAL8  f1dis;
	    REAL8 significance;

	    ht.f0Bin = fBinSearch;
	    ht.spinRes.length = 1;
	    ht.spinRes.data = NULL;
	    ht.spinRes.data = (REAL8 *)LALCalloc(ht.spinRes.length, sizeof(REAL8));
	    
	    for ( n = 0; n <= floor(nSpin1Max/uvar_spindownJump); ++n) { 
	      /*loop over all spindown values */

	      f1dis = - n * f1jump;
	      ht.spinRes.data[0] =  f1dis * deltaF;

	      /* construct path in time-freq plane */	      
	      for (j = 0 ; j < mObsCohBest; ++j){
		freqInd.data[j] = fBinSearch + floor(best.timeDiffV->data[j]*f1dis + 0.5);
	      }
	      
	      if (uvar_weighAM || uvar_weighNoise) {
		LAL_CALL( LALHOUGHConstructHMT_W( &status, &ht, &freqInd, &phmdVS ), &status );
	      }
	      else {
		LAL_CALL( LALHOUGHConstructHMT( &status, &ht, &freqInd, &phmdVS ), &status );	  
	      }


	      /* ********************* perfom stat. analysis on the maps ****************** */
	      
	      if ( uvar_EnableExtraInfo ) {

		LAL_CALL( LALHoughStatistics ( &status, &stats, &ht), &status );
	        LAL_CALL( LALStereo2SkyLocation (&status, &sourceLocation, 
				       stats.maxIndex[0], stats.maxIndex[1], &patch, &parDem), &status);
		
		/*LAL_CALL( LALHoughHistogram ( &status, &hist, &ht), &status);*/
		LAL_CALL( LALHoughHistogramSignificance ( &status, hist, &ht, meanN, sigmaN, 
							  minSignificance, maxSignificance), &status);

		for(j = 0; j < histTotal->length; j++){ 
		  histTotal->data[j] += hist->data[j]; 
		}
		significance =  (stats.maxCount - meanN)/sigmaN;   
	      }	      	      

	      /* select candidates from hough maps */
	      LAL_CALL( GetToplistFromHoughmap( &status, toplist, &ht, &patch, &parDem, meanN, sigmaN), &status);


	      /* ***** print results *********************** */

	      if( uvar_EnableExtraInfo )
	      {
		  if( PrintExtraInfo( fileMaps, &fp1, iHmap, &ht, &sourceLocation, &stats, fBinSearch, deltaF))
		      return DRIVEHOUGHCOLOR_EFILE;
	      }

	      ++iHmap;
	    } /* end loop over spindown values */ 
	    
	    LALFree(ht.spinRes.data);
	    
	    
	    /***** shift the search freq. & PHMD structure 1 freq.bin ****** */
	    ++fBinSearch;
	    
	    LAL_CALL( LALHOUGHupdateSpacePHMDup(&status, &phmdVS, best.pgV, &lutV), &status );
	    
	    if (uvar_weighAM || uvar_weighNoise) {
	      LAL_CALL( LALHOUGHWeighSpacePHMD( &status, &phmdVS, best.weightsV), &status);	    
	    }

	  }   /*closing second while */
	
	fBin = fBinSearch;
	
	/* ********************  Free partial memory ******************* */
	LALFree(patch.xCoor);
	LALFree(patch.yCoor);
	LALFree(ht.map);
	
	LALHOUGHDestroyLUTs( &status, &lutV);
	
	LALHOUGHDestroyPHMDs( &status, &phmdVS);
	
	
      } /* closing while */
      

      /* printing total histogram */
      if ( uvar_EnableExtraInfo ) 
	{
	  if( PrintHistogram( histTotal, filehisto, minSignificance, maxSignificance) ) return 7;
	}

      /* --------------------------------------------------*/
      /* Closing files with statistics results and events*/
      if (uvar_EnableExtraInfo) fclose(fp1);
      
      /* Free memory allocated inside skypatches loop */
      LALFree(lutV.lut);  
      lutV.lut = NULL;
      
      LALFree(phmdVS.phmd);
      phmdVS.phmd = NULL;
      
      LALFree(freqInd.data);
      freqInd.data = NULL;
      
      if ( uvar_EnableExtraInfo ) {
	XLALDestroyUINT8Vector (hist);
	XLALDestroyUINT8Vector (histTotal);
      }
      
    } /* finish loop over skypatches */
    LogPrintfVerbatim (LOG_NORMAL, "...done\n");


  /* close sigma file */
  if ( uvar_EnableExtraInfo )
    fclose(fpSigma);

  /*********************************************************/
  /* print toplist */
  /********************************************************/
  
  /* If we want to print Chi2 value */
  
  if (uvar_EnableChi2){
    LogPrintf (LOG_NORMAL, "Starting chi-square follow-up of top candidates...");
    LAL_CALL(ComputeandPrintChi2(&status, toplist, timeDiffV, &velV, uvar_chiSqBins, alphaPeak, mdetStates, weightsNoise, &upgV), &status);    
    LogPrintfVerbatim (LOG_NORMAL, "done\n");
  }
  else {
    
    FILE   *fpToplist = NULL;

    LogPrintf (LOG_NORMAL, "Sort and print toplist...");
    fpToplist = fopen("hough_top.dat","w");
    
    sort_fstat_toplist(toplist);
    
    if ( write_fstat_toplist_to_fp( toplist, fpToplist, NULL) < 0)
      LogPrintf( LOG_CRITICAL, "error in writing toplist to file...\n");
    
    if (fprintf(fpToplist,"%%DONE\n") < 0)
      LogPrintf(LOG_CRITICAL, "error writing end marker...\n");
    
    fclose(fpToplist);
    LogPrintfVerbatim (LOG_NORMAL, "done\n");
  }
  /********************************************************/
  
  LogPrintf (LOG_NORMAL, "Free memory and exit...");
  {
    UINT4 j;
    for (j = 0; j < mObsCoh; ++j) LALFree( pgV.pg[j].peak); 
  }
  
  LALFree(pgV.pg);
  
  if (uvar_EnableChi2) 
    {
      {UINT4 j;
      for (j = 0; j < mObsCoh; ++j) LALFree( upgV.upg[j].data);
      }
      LALFree(upgV.upg);
      
    }
  
  XLALDestroyTimestampVector ( timeV);
  XLALDestroyREAL8Vector( timeDiffV);
  

  LALFree(velV.data);

  XLALDestroyREAL8Vector(weightsV);
  XLALDestroyREAL8Vector(weightsNoise);

  XLALDestroyMultiDetectorStateSeries ( mdetStates );

  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);

  LALFree(skyAlpha);
  LALFree(skyDelta);
  LALFree(skySizeAlpha);
  LALFree(skySizeDelta);

  if (nStarEventVec.event) {
    LALFree(nStarEventVec.event);
  }

  LALFree(best.weightsV->data);  
  LALFree(best.weightsV);
  LALFree(best.timeDiffV->data);
  LALFree(best.timeDiffV);
  LALFree(best.velV->data);
  LALFree(best.velV);
  LALFree(best.pgV->pg);
  LALFree(best.pgV);

  free_fstat_toplist(&toplist);

  LAL_CALL (LALDestroyUserVars(&status), &status);

  LALCheckMemoryLeaks();

  LogPrintfVerbatim (LOG_NORMAL, "bye\n");

  /*   if ( lalDebugLevel ) */
  /*     REPORTSTATUS ( &status); */

  return status.statusCode;
}



  
/******************************************************************/
/* printing the Histogram of all maps into a file                    */
/******************************************************************/

int CreateSkypatchDirs(CHAR *filestats,
		       CHAR *base,
		       INT4 index)
{
  CHAR tempstr[16];
  INT4 mkdir_result;

  /* create the directory name uvar_dirnameOut/skypatch_$j */
  strcpy(  filestats, base);

  strcat( filestats, "/skypatch_");
  sprintf( tempstr, "%d", index);
  strcat( filestats, tempstr);  
  strcat( filestats, "/");
  
  /* check whether file can be created or if it exists already 
     if not then exit */
  errno = 0;  
  mkdir_result = mkdir(filestats, S_IRWXU | S_IRWXG | S_IRWXO);

  if ( (mkdir_result == -1) && (errno != EEXIST) )
    {
      fprintf(stderr, "unable to create skypatch directory %d\n", index);
      return 1;
    }

  return 0;
  
}

  
int PrintHistogram(UINT8Vector *hist, CHAR *fnameOut, REAL8 minSignificance, REAL8 maxSignificance)
{

  INT4 binsHisto;
  REAL8 dSig;
  FILE  *fp=NULL;   /* Output file */
  char filename[ HOUGHMAXFILENAMELENGTH ];
  INT4  i ;
  
  strcpy(  filename, fnameOut);
  strcat(  filename, "histo");

  binsHisto = hist->length;
  dSig = (maxSignificance - minSignificance)/binsHisto;
  
  if ( (fp = fopen(filename,"w")) == NULL)
    {  
      fprintf(stderr,"Unable to find file %s\n",filename);
      return DRIVEHOUGHCOLOR_EFILE; 
    }

  for (i = 0; i < binsHisto; i++){
    fprintf(fp,"%g  %llu\n", minSignificance + i*dSig, hist->data[i]);
  }
  
  fclose( fp );  
  return 0;

}



/******************************************************************/
/* printing the HM into a file                    */
/******************************************************************/

int PrintHmap2file(HOUGHMapTotal *ht, CHAR *fnameOut, INT4 iHmap){

  FILE  *fp=NULL;   /* Output file */
  char filename[ HOUGHMAXFILENAMELENGTH ], filenumber[16]; 
  INT4  k, i ;
  UINT2 xSide, ySide;
   
  strcpy(  filename, fnameOut);
  sprintf( filenumber, ".%06d",iHmap); 
  strcat(  filename, filenumber);

  if ( (fp = fopen(filename,"w")) == NULL)
    {  
      fprintf(stderr,"Unable to find file %s\n",filename);
      return DRIVEHOUGHCOLOR_EFILE; 
    }

  ySide= ht->ySide;
  xSide= ht->xSide;

  for(k=ySide-1; k>=0; --k){
    for(i=0;i<xSide;++i){
      fprintf( fp ," %f", (REAL4)ht->map[k*xSide +i]);
      fflush( fp );
    }
    fprintf( fp ," \n");
    fflush( fp );
  }

  fclose( fp );  
  return 0;
}


/******************************************************************/
/* printing the HM into a m_file                    */
/******************************************************************/

int PrintHmap2m_file(HOUGHMapTotal *ht, CHAR *fnameOut, INT4 iHmap){

  FILE  *fp=NULL;   /* Output file */
  char filename[ HOUGHMAXFILENAMELENGTH ], filenumber[16]; 
  INT4  k, i ;
  UINT2 xSide, ySide;
  INT4 mObsCoh;
  REAL8 f0,f1;
   
  strcpy(  filename, fnameOut);
  sprintf( filenumber, "%06d.m",iHmap); 
  strcat(  filename, filenumber);
  fp=fopen(filename,"w");

  if ( !fp ){  
    fprintf(stderr,"Unable to find file %s\n",filename);
    return DRIVEHOUGHCOLOR_EFILE; 
  }

  ySide= ht->ySide;
  xSide= ht->xSide;
  f0=ht->f0Bin* ht->deltaF;
  mObsCoh = ht->mObsCoh;
  f1=0.0;
  if( ht->spinRes.length ){ f1=ht->spinRes.data[0]; }
  
  /* printing into matlab format */
  
  fprintf( fp ,"f0= %f ; \n", f0);
  fprintf( fp ,"f1= %g ; \n", f1);
  fprintf( fp ,"map= [ \n");

  for(k=ySide-1; k>=0; --k){
    for(i=0;i<xSide;++i){
      fprintf( fp ," %f", (REAL4)ht->map[k*xSide +i]);
      fflush( fp );
    }
    fprintf( fp ," \n");
    fflush( fp );
  }
  
  fprintf( fp ,"    ]; \n");
  fclose( fp );  
  return 0;
}

/******************************************************************************************/

void PrintLogFile (LALStatus       *status,
		   CHAR            *dir,
		   CHAR            *basename,
		   CHAR            *skyfile,
		   LALStringVector *linefiles,
		   CHAR            *executable )
{
  CHAR *fnameLog=NULL; 
  FILE *fpLog=NULL;
  CHAR *logstr=NULL; 
  UINT4 k;

  INITSTATUS (status, "PrintLogFile", rcsid);
  ATTATCHSTATUSPTR (status);
  
  /* open log file for writing */
  fnameLog = (CHAR *)LALCalloc( HOUGHMAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(fnameLog,dir);
  strcat(fnameLog, "/logfiles/");
  /* now create directory fdirOut/logfiles using mkdir */
  errno = 0;
  {
    /* check whether file can be created or if it exists already 
       if not then exit */
    INT4 mkdir_result;
    mkdir_result = mkdir(fnameLog, S_IRWXU | S_IRWXG | S_IRWXO);
    if ( (mkdir_result == -1) && (errno != EEXIST) )
      {
	fprintf(stderr, "unable to create logfiles directory %s\n", fnameLog);
        LALFree(fnameLog);
	exit(1);  /* stop the program */
      }
  }

  /* create the logfilename in the logdirectory */
  strcat(fnameLog, basename);
  strcat(fnameLog,".log");
  /* open the log file for writing */
  if ((fpLog = fopen(fnameLog, "w")) == NULL) {
    fprintf(stderr, "Unable to open file %s for writing\n", fnameLog);
    LALFree(fnameLog);
    exit(1);
  }
  
  /* get the log string */
  TRY( LALUserVarGetLog(status->statusPtr, &logstr, UVAR_LOGFMT_CFGFILE), status);  

  fprintf( fpLog, "## LOG FILE FOR Hough Driver\n\n");
  fprintf( fpLog, "# User Input:\n");
  fprintf( fpLog, "#-------------------------------------------\n");
  fprintf( fpLog, logstr);
  LALFree(logstr);

  /* copy contents of skypatch file into logfile */
  fprintf(fpLog, "\n\n# Contents of skypatch file:\n");
  fclose(fpLog);
  {
    CHAR command[1024] = "";
    sprintf(command, "cat %s >> %s", skyfile, fnameLog);
    system(command);
  }


  /* copy contents of linefile if necessary */
  if ( linefiles ) {

    for ( k = 0; k < linefiles->length; k++) {
      
      if ((fpLog = fopen(fnameLog, "a")) != NULL) {
	CHAR command[1024] = "";
	fprintf (fpLog, "\n\n# Contents of linefile %s :\n", linefiles->data[k]);
	fprintf (fpLog, "# -----------------------------------------\n");
	fclose (fpLog);
	sprintf(command, "cat %s >> %s", linefiles->data[k], fnameLog);      
	system (command);	 
      } 
    } 
  }

  /* append an ident-string defining the exact CVS-version of the code used */
  if ((fpLog = fopen(fnameLog, "a")) != NULL) 
    {
      CHAR command[1024] = "";
      fprintf (fpLog, "\n\n# CVS-versions of executable:\n");
      fprintf (fpLog, "# -----------------------------------------\n");
      fclose (fpLog);
      
      sprintf (command, "ident %s | sort -u >> %s", executable, fnameLog);
      system (command);	/* we don't check this. If it fails, we assume that */
    			/* one of the system-commands was not available, and */
    			/* therefore the CVS-versions will not be logged */ 
    }

  LALFree(fnameLog); 
  	 
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}    


/* read timestamps file */
void ReadTimeStampsFile (LALStatus          *status,
			 LIGOTimeGPSVector  *ts,
			 CHAR               *filename)
{

  FILE  *fp = NULL;
  INT4  numTimeStamps, r;
  UINT4 j;
  REAL8 temp1, temp2;

  INITSTATUS (status, "ReadTimeStampsFile", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT(ts, status, DRIVEHOUGHCOLOR_ENULL,DRIVEHOUGHCOLOR_MSGENULL); 
  ASSERT(ts->data == NULL, status, DRIVEHOUGHCOLOR_ENULL,DRIVEHOUGHCOLOR_MSGENULL); 
  ASSERT(ts->length == 0, status, DRIVEHOUGHCOLOR_ENULL,DRIVEHOUGHCOLOR_MSGENULL); 
  ASSERT(filename, status, DRIVEHOUGHCOLOR_ENULL,DRIVEHOUGHCOLOR_MSGENULL); 

  if ( (fp = fopen(filename, "r")) == NULL) {
    ABORT( status, DRIVEHOUGHCOLOR_EFILE, DRIVEHOUGHCOLOR_MSGEFILE);
  }

  /* count number of timestamps */
  numTimeStamps = 0;     

  do {
    r = fscanf(fp,"%lf%lf\n", &temp1, &temp2);
    /* make sure the line has the right number of entries or is EOF */
    if (r==2) numTimeStamps++;
  } while ( r != EOF);
  rewind(fp);

  ts->length = numTimeStamps;
  ts->data = LALCalloc (1, numTimeStamps * sizeof(LIGOTimeGPS));;
  if ( ts->data == NULL ) {
    fclose(fp);
    ABORT( status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  }
  
  for (j = 0; j < ts->length; j++)
    {
      r = fscanf(fp,"%lf%lf\n", &temp1, &temp2);
      ts->data[j].gpsSeconds = (INT4)temp1;
      ts->data[j].gpsNanoSeconds = (INT4)temp2;
    }
  
  fclose(fp);
  	 
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}    



/** given a total hough map, this function produces a histogram of
   the number count significance */
void LALHoughHistogramSignificance(LALStatus      *status,
				   UINT8Vector    *out,
				   HOUGHMapTotal  *in,
				   REAL8          mean,
				   REAL8          sigma,
				   REAL8          minSignificance,
				   REAL8          maxSignificance)
{

  INT4   i, j, binsHisto, xSide, ySide, binIndex;
  REAL8  temp;

  INITSTATUS (status, "LALHoughHistogramSignificance", rcsid);
  ATTATCHSTATUSPTR (status);

  /* make sure arguments are not null */
  ASSERT (in, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (in->map, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (out, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);

  /* make sure input hough map is ok*/
  ASSERT (in->xSide > 0, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  ASSERT (in->ySide > 0, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  ASSERT (in->mObsCoh > 0, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);

  ASSERT (out->length > 0, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);

  binsHisto = out->length;
  xSide = in->xSide;
  ySide = in->ySide;

  /* initialize histogram vector*/
  for (i=0; i < binsHisto; i++) 
    out->data[i] = 0;

  /* loop over hough map and find histogram */
  for (i = 0; i < ySide; i++){
    for (j = 0; j < xSide; j++){

      /* calculate significance of number count */
      temp = (in->map[i*xSide + j] - mean)/sigma;

      /* make sure temp is in proper range */
      ASSERT (temp > minSignificance, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
      ASSERT (temp < maxSignificance, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);

      binIndex = (INT4) binsHisto * (temp - minSignificance)/(maxSignificance - minSignificance);

      /* make sure binIndex is in proper range */
      ASSERT (binIndex >= 0, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
      ASSERT (binIndex <= binsHisto-1, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);

      /* add to relevant entry in histogram */
      out->data[binIndex] += 1;
    }
  }

  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);
}



void GetSFTVelTime(LALStatus                *status,
		   REAL8Cart3CoorVector     *velV,
		   LIGOTimeGPSVector        *timeV, 
		   MultiDetectorStateSeries *in)
{

  UINT4 numifo, numsft, iIFO, iSFT, j;  
  
  INITSTATUS (status, "GetSFTVelTime", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT (in, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (in->length > 0, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);

  ASSERT (velV, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (velV->length > 0, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (velV->data, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);

  ASSERT (timeV, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (timeV->length > 0, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (timeV->data, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);

  ASSERT (velV->length == timeV->length, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);

  numifo = in->length;
  
  /* copy the timestamps, weights, and velocity vector */
  for (j = 0, iIFO = 0; iIFO < numifo; iIFO++ ) {

    ASSERT (in->data[iIFO], status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);    
    
    numsft = in->data[iIFO]->length;    
    ASSERT (numsft > 0, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);    

    for ( iSFT = 0; iSFT < numsft; iSFT++, j++) {
      
      velV->data[j].x = in->data[iIFO]->data[iSFT].vDetector[0];
      velV->data[j].y = in->data[iIFO]->data[iSFT].vDetector[1];
      velV->data[j].z = in->data[iIFO]->data[iSFT].vDetector[2];
      
      /* mid time of sfts */
      timeV->data[j] = in->data[iIFO]->data[iSFT].tGPS;
      
    } /* loop over SFTs */
    
  } /* loop over IFOs */
  
  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);
}




void GetSFTNoiseWeights(LALStatus          *status,
			REAL8Vector        *out,
			MultiNoiseWeights  *in)
{

  UINT4 numifo, numsft, iIFO, iSFT, j;  
  
  INITSTATUS (status, "GetSFTNoiseWeights", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT (in, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (in->length > 0, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);

  ASSERT (out, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (out->length > 0, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (out->data, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);

  numifo = in->length;
  
  /* copy the timestamps, weights, and velocity vector */
  for (j = 0, iIFO = 0; iIFO < numifo; iIFO++ ) {

    ASSERT (in->data[iIFO], status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);    
    
    numsft = in->data[iIFO]->length;    
    ASSERT (numsft > 0, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);    

    for ( iSFT = 0; iSFT < numsft; iSFT++, j++) {

      out->data[j] = in->data[iIFO]->data[iSFT];      
      
    } /* loop over SFTs */
    
  } /* loop over IFOs */

  TRY( LALHOUGHNormalizeWeights( status->statusPtr, out), status);
  
  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);
}


/** Loop over SFTs and set a threshold to get peakgrams.  SFTs must be normalized.  */
void GetPeakGramFromMultSFTVector(LALStatus           *status,
				  HOUGHPeakGramVector *out, /**< Output peakgrams */
				  MultiSFTVector      *in,  /**< Input SFTs */
				  REAL8               thr   /**< Threshold on SFT power */)  
{

  SFTtype  *sft;
  UCHARPeakGram  pg1;
  INT4   nPeaks;
  UINT4  iIFO, iSFT, numsft, numifo, j, binsSFT; 
  
  INITSTATUS (status, "GetSFTNoiseWeights", rcsid);
  ATTATCHSTATUSPTR (status);

  numifo = in->length;
  binsSFT = in->data[0]->data->data->length;
  
  pg1.length = binsSFT;
  pg1.data = NULL;
  pg1.data = (UCHAR *)LALCalloc( 1, binsSFT*sizeof(UCHAR));
  
  /* loop over sfts and select peaks */
  for ( j = 0, iIFO = 0; iIFO < numifo; iIFO++){
    
    numsft = in->data[iIFO]->length;
    
    for ( iSFT = 0; iSFT < numsft; iSFT++, j++) {
      
      sft = in->data[iIFO]->data + iSFT;
      
      TRY (SFTtoUCHARPeakGram( status->statusPtr, &pg1, sft, thr), status);
      
      nPeaks = pg1.nPeaks;

      /* compress peakgram */      
      out->pg[j].length = nPeaks;
      out->pg[j].peak = NULL;
      out->pg[j].peak = (INT4 *)LALCalloc( 1, nPeaks*sizeof(INT4));
      
      TRY( LALUCHAR2HOUGHPeak( status->statusPtr, &(out->pg[j]), &pg1), status );
      
    } /* loop over SFTs */

  } /* loop over IFOs */

  LALFree(pg1.data);

  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}


/** Set up location of skypatch centers and sizes 
    If user specified skyRegion then use DopplerScan function
    to construct an isotropic grid. Otherwise use skypatch file. */
void SetUpSkyPatches(LALStatus           *status,
		     HoughSkyPatchesInfo *out,   /**< output skypatches info */
		     CHAR                *skyFileName, /**< name of skypatch file */
		     CHAR                *skyRegion,  /**< skyregion (if isotropic grid is to be constructed) */
		     REAL8               dAlpha,      /**< alpha resolution (if isotropic grid is to be constructed) */
		     REAL8               dDelta,  /**< delta resolution (if isotropic grid is to be constructed) */
		     INT4                numSkyPartitions,/**<Number of (equi-)partitions to split skygrid into */                   
                     INT4                partitionIndex)/**< Index [0,numSkyPartitions-1] of sky-partition to generate */
{

  DopplerSkyScanInit scanInit = empty_DopplerSkyScanInit;   /* init-structure for DopperScanner */
  DopplerSkyScanState thisScan = empty_DopplerSkyScanState; /* current state of the Doppler-scan */
  UINT4 nSkyPatches, skyCounter;
  PulsarDopplerParams dopplerpos;	  
  
  INITSTATUS (status, "SetUpSkyPatches", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT (out, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (dAlpha > 0, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  ASSERT (dDelta > 0, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);

  if (skyRegion ) {
    
    scanInit.dAlpha = dAlpha;
    scanInit.dDelta = dDelta;
    scanInit.gridType = GRID_ISOTROPIC;
    scanInit.metricType =  LAL_PMETRIC_NONE;
    /* scanInit.metricMismatch = 0; */
    /* scanInit.projectMetric = TRUE; */
    /* scanInit.obsDuration = tStack; */
    /* scanInit.obsBegin = tMidGPS; */
    /* scanInit.Detector = &(stackMultiDetStates.data[0]->data[0]->detector); */
    /* scanInit.ephemeris = edat; */
    /* scanInit.skyGridFile = uvar_skyGridFile; */
    scanInit.skyRegionString = (CHAR*)LALCalloc(1, strlen(skyRegion)+1);
    strcpy (scanInit.skyRegionString, skyRegion);
    /*   scanInit.Freq = usefulParams.spinRange_midTime.fkdot[0] +  usefulParams.spinRange_midTime.fkdotBand[0]; */
    
    scanInit.numSkyPartitions = numSkyPartitions;
    scanInit.partitionIndex = partitionIndex;

    /* set up the grid */
    TRY ( InitDopplerSkyScan ( status->statusPtr, &thisScan, &scanInit), status); 
    
    nSkyPatches = out->numSkyPatches = thisScan.numSkyGridPoints;
    
    out->alpha = (REAL8 *)LALCalloc(1, nSkyPatches*sizeof(REAL8));
    out->delta = (REAL8 *)LALCalloc(1, nSkyPatches*sizeof(REAL8));     
    out->alphaSize = (REAL8 *)LALCalloc(1, nSkyPatches*sizeof(REAL8));
    out->deltaSize = (REAL8 *)LALCalloc(1, nSkyPatches*sizeof(REAL8));     
        
    /* loop over skygrid points */  
    XLALNextDopplerSkyPos(&dopplerpos, &thisScan);
    
    skyCounter = 0; 
    while(thisScan.state != STATE_FINISHED) {
      
      out->alpha[skyCounter] = dopplerpos.Alpha;
      out->delta[skyCounter] = dopplerpos.Delta;
      out->alphaSize[skyCounter] = dAlpha;
      out->deltaSize[skyCounter] = dDelta;
      
      /*  if ((dopplerpos.Delta>0) && (dopplerpos.Delta < atan(4*LAL_PI/dAlpha/dDelta) ))
        out->alphaSize[skyCounter] = dAlpha*cos(dopplerpos.Delta -0.5*dDelta)/cos(dopplerpos.Delta);

      if ((dopplerpos.Delta<0) && (dopplerpos.Delta > -atan(4*LAL_PI/dAlpha/dDelta) ))
        out->alphaSize[skyCounter] = dAlpha*cos(dopplerpos.Delta +0.5*dDelta)/cos(dopplerpos.Delta);
      */      

      XLALNextDopplerSkyPos(&dopplerpos, &thisScan);
      skyCounter++;
      
    } /* end while loop over skygrid */
      
    /* free dopplerscan stuff */
    TRY ( FreeDopplerSkyScan( status->statusPtr, &thisScan), status);
    if ( scanInit.skyRegionString )
      LALFree ( scanInit.skyRegionString );

  } else {

    /* read skypatch info */
    {
      FILE   *fpsky = NULL; 
      INT4   r;
      REAL8  temp1, temp2, temp3, temp4;

      ASSERT (skyFileName, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
      
      if ( (fpsky = fopen(skyFileName, "r")) == NULL)
	{
	  ABORT ( status, DRIVEHOUGHCOLOR_EFILE, DRIVEHOUGHCOLOR_MSGEFILE );
	}
      
      nSkyPatches = 0;
      do 
	{
	  r = fscanf(fpsky,"%lf%lf%lf%lf\n", &temp1, &temp2, &temp3, &temp4);
	  /* make sure the line has the right number of entries or is EOF */
	  if (r==4) nSkyPatches++;
	} while ( r != EOF);
      rewind(fpsky);

      out->numSkyPatches = nSkyPatches;      
      out->alpha = (REAL8 *)LALCalloc(nSkyPatches, sizeof(REAL8));
      out->delta = (REAL8 *)LALCalloc(nSkyPatches, sizeof(REAL8));     
      out->alphaSize = (REAL8 *)LALCalloc(nSkyPatches, sizeof(REAL8));
      out->deltaSize = (REAL8 *)LALCalloc(nSkyPatches, sizeof(REAL8));     
      
      for (skyCounter = 0; skyCounter < nSkyPatches; skyCounter++)
	{
	  r = fscanf(fpsky,"%lf%lf%lf%lf\n", out->alpha + skyCounter, out->delta + skyCounter, 
		     out->alphaSize + skyCounter,  out->deltaSize + skyCounter);
	}
      
      fclose(fpsky);     
    } /* end skypatchfile reading block */    
  } /* end setting up of skypatches */


  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}


void GetAMWeights(LALStatus                *status,
		  REAL8Vector              *out,
		  MultiDetectorStateSeries *mdetStates,		  
		  REAL8                    alpha,
		  REAL8                    delta)
{

  MultiAMCoeffs *multiAMcoef = NULL;
  UINT4 iIFO, iSFT, k, numsft, numifo;
  REAL8 a, b;
  SkyPosition skypos;
  
  INITSTATUS (status, "GetAMWeights", rcsid);
  ATTATCHSTATUSPTR (status);
  
  /* get the amplitude modulation coefficients */
  skypos.longitude = alpha;
  skypos.latitude = delta;
  skypos.system = COORDINATESYSTEM_EQUATORIAL;
  TRY ( LALGetMultiAMCoeffs ( status->statusPtr, &multiAMcoef, mdetStates, skypos), status);

  numifo = mdetStates->length;
  
  /* loop over the weights and multiply them by the appropriate
     AM coefficients */
  for ( k = 0, iIFO = 0; iIFO < numifo; iIFO++) {
    
    numsft = mdetStates->data[iIFO]->length;
    
    for ( iSFT = 0; iSFT < numsft; iSFT++, k++) {	  
            
      a = multiAMcoef->data[iIFO]->a->data[iSFT];
      b = multiAMcoef->data[iIFO]->b->data[iSFT];    
      out->data[k] *= (a*a + b*b);

    } /* loop over SFTs */

  } /* loop over IFOs */
  
  TRY( LALHOUGHNormalizeWeights( status->statusPtr, out), status);
  
  XLALDestroyMultiAMCoeffs ( multiAMcoef );
  

  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}



void SelectBestStuff(LALStatus      *status,
		     BestVariables  *out,
		     BestVariables  *in,		  
		     UINT4          mObsCohBest)
{

  size_t *index=NULL;
  UINT4 k, mObsCoh;

  INITSTATUS (status, "SelectBestStuff", rcsid);
  ATTATCHSTATUSPTR (status);

  /* check consistency of input */
  ASSERT (out, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);

  ASSERT (in, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (in->length > 0, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  mObsCoh = in->length;

  ASSERT (in->weightsV, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (in->weightsV->length == mObsCoh, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);

  ASSERT (in->timeDiffV, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (in->timeDiffV->length == mObsCoh, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);

  ASSERT (in->velV, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (in->velV->length == mObsCoh, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);

  ASSERT (mObsCohBest > 0, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  ASSERT (mObsCohBest <= in->length, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);

  /* memory allocation for output -- use reallocs because this 
     function may be called within the loop over sky positions, so memory might
     already have been allocated previously */
  out->length = mObsCohBest;


  if (out->weightsV == NULL) {
    out->weightsV = (REAL8Vector *)LALCalloc(1, sizeof(REAL8Vector));
    out->weightsV->length = mObsCohBest;
    out->weightsV->data = (REAL8 *)LALCalloc(1, mObsCohBest*sizeof(REAL8));	
  }

  if (out->timeDiffV == NULL) {
    out->timeDiffV = (REAL8Vector *)LALCalloc(1, sizeof(REAL8Vector));
    out->timeDiffV->length = mObsCohBest;
    out->timeDiffV->data = (REAL8 *)LALCalloc(1, mObsCohBest*sizeof(REAL8));	
  }

  if (out->velV == NULL) {
    out->velV =  (REAL8Cart3CoorVector *)LALCalloc(1, sizeof(REAL8Cart3CoorVector));
    out->velV->length = mObsCohBest;
    out->velV->data = (REAL8Cart3Coor *)LALCalloc(1, mObsCohBest*sizeof(REAL8Cart3Coor));
  }

  if (out->pgV == NULL) {
    out->pgV = (HOUGHPeakGramVector *)LALCalloc(1, sizeof(HOUGHPeakGramVector));
    out->pgV->length = mObsCohBest;
    out->pgV->pg = (HOUGHPeakGram *)LALCalloc(1, mObsCohBest*sizeof(HOUGHPeakGram));
  }

  index = LALCalloc(1, mObsCohBest*sizeof(size_t));  
  gsl_sort_largest_index( index, mObsCohBest, in->weightsV->data, 1, mObsCoh);	

  for ( k = 0; k < mObsCohBest; k++) {

    out->weightsV->data[k] = in->weightsV->data[index[k]];    
    out->timeDiffV->data[k] = in->timeDiffV->data[index[k]];    

    out->velV->data[k].x = in->velV->data[index[k]].x;
    out->velV->data[k].y = in->velV->data[index[k]].y;
    out->velV->data[k].z = in->velV->data[index[k]].z;

    /* this copies the pointers to the peakgram data from the input */
    out->pgV->pg[k] = in->pgV->pg[index[k]];

  }

  /*   gsl_sort_index( index, timeDiffV.data, 1, mObsCohBest);	 */
  
  LALFree(index);

  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}


/* copy data in BestVariables struct */
void DuplicateBestStuff(LALStatus      *status,
			BestVariables  *out,
			BestVariables  *in)
{

  UINT4 mObsCoh;

  INITSTATUS (status, "SelectBestStuff", rcsid);
  ATTATCHSTATUSPTR (status);

  /* check consistency of input */
  ASSERT (out, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);

  ASSERT (in, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (in->length > 0, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  mObsCoh = in->length;

  ASSERT (in->weightsV, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (in->weightsV->length == mObsCoh, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);

  ASSERT (in->timeDiffV, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (in->timeDiffV->length == mObsCoh, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);

  ASSERT (in->velV, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (in->velV->length == mObsCoh, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);

  /* memory allocation for output -- check output is null because
     function may be called within the loop over sky positions, so memory might
     already have been allocated previously */
  out->length = mObsCoh;

  if (out->weightsV == NULL) {
    out->weightsV = (REAL8Vector *)LALCalloc(1, sizeof(REAL8Vector));
    out->weightsV->length = mObsCoh;
    out->weightsV->data = (REAL8 *)LALCalloc(1, mObsCoh*sizeof(REAL8));	
  }

  if (out->timeDiffV == NULL) {
    out->timeDiffV = (REAL8Vector *)LALCalloc(1, sizeof(REAL8Vector));
    out->timeDiffV->length = mObsCoh;
    out->timeDiffV->data = (REAL8 *)LALCalloc(1, mObsCoh*sizeof(REAL8));	
  }

  if (out->velV == NULL) {
    out->velV =  (REAL8Cart3CoorVector *)LALCalloc(1, sizeof(REAL8Cart3CoorVector));
    out->velV->length = mObsCoh;
    out->velV->data = (REAL8Cart3Coor *)LALCalloc(1, mObsCoh*sizeof(REAL8Cart3Coor));
  }

  if (out->pgV == NULL) {
    out->pgV = (HOUGHPeakGramVector *)LALCalloc(1, sizeof(HOUGHPeakGramVector));
    out->pgV->length = mObsCoh;
    out->pgV->pg = (HOUGHPeakGram *)LALCalloc(1, mObsCoh*sizeof(HOUGHPeakGram));
  }

  memcpy(out->weightsV->data, in->weightsV->data, mObsCoh * sizeof(REAL8));
  memcpy(out->timeDiffV->data, in->timeDiffV->data, mObsCoh * sizeof(REAL8));
  memcpy(out->velV->data, in->velV->data, mObsCoh * sizeof(REAL8Cart3Coor));
  memcpy(out->pgV->pg, in->pgV->pg, mObsCoh * sizeof(HOUGHPeakGram));
  
  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}





/** helper function for creating LUT vector */
void LALHOUGHCreateLUTVector(LALStatus           *status,
			     HOUGHptfLUTVector   *lutV,
			     UINT4               length)
{

  INITSTATUS (status, "LALHOUGHCreateLUTVector", rcsid);
  ATTATCHSTATUSPTR (status);

  /* check input pars are ok */
  if (lutV == NULL) {
    ABORT (status, DRIVEHOUGHCOLOR_ENONULL, DRIVEHOUGHCOLOR_MSGENONULL);
  }

  if (lutV->lut != NULL) {
    ABORT (status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  }

  if (length <= 0) {
    ABORT (status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  }

  /* allocate memory */
  lutV->length = length;
  lutV->lut = (HOUGHptfLUT *)LALCalloc(lutV->length, sizeof(HOUGHptfLUT));

  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}


void LALHOUGHCreateLUTs(LALStatus           *status,
			HOUGHptfLUTVector   *lutV,
			UINT2               maxNBins, 
			UINT2               maxNBorders,
			UINT2               ySide)
{

  UINT4 j,i;

  INITSTATUS (status, "LALHOUGHCreateLUTs", rcsid);
  ATTATCHSTATUSPTR (status);

  /* check input pars are ok */
  if (lutV == NULL) {
    ABORT (status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  }

  if (maxNBins <= 0) {
    ABORT (status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  }

  if (maxNBorders <= 0) {
    ABORT (status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  }

  if (ySide <= 0) {
    ABORT (status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  }


  /* loop over luts and allocate memory */
  for(j = 0; j < lutV->length; j++){

    if (lutV->lut + j == NULL) {
      ABORT (status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
    }

    lutV->lut[j].maxNBins = maxNBins;
    lutV->lut[j].maxNBorders = maxNBorders;
    lutV->lut[j].border = (HOUGHBorder *)LALCalloc(maxNBorders, sizeof(HOUGHBorder));
    lutV->lut[j].bin = (HOUGHBin2Border *)LALCalloc(maxNBins, sizeof(HOUGHBin2Border));

    for (i = 0; i < maxNBorders; i++){
      lutV->lut[j].border[i].ySide = ySide;
      lutV->lut[j].border[i].xPixel = (COORType *)LALCalloc(ySide, sizeof(COORType));
    }

  }
  
  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}


void LALHOUGHDestroyLUTs(LALStatus           *status,
			 HOUGHptfLUTVector   *lutV)
{

  UINT4 j,i;

  INITSTATUS (status, "LALHOUGHDestroyLUTs", rcsid);
  ATTATCHSTATUSPTR (status);

  for (j = 0; j < lutV->length ; ++j){
    for (i = 0; i < lutV->lut[j].maxNBorders; ++i){

      if (lutV->lut[j].border[i].xPixel) {
	LALFree( lutV->lut[j].border[i].xPixel);
      }
    }

    if (lutV->lut[j].border) {
      LALFree( lutV->lut[j].border);
    }

    if (lutV->lut[j].bin) {
      LALFree( lutV->lut[j].bin);
    }
  }

  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}


void LALHOUGHCreatePHMDVS(LALStatus           *status,
			  PHMDVectorSequence  *phmdVS,
			  UINT4               length,
			  UINT4               nfSize)
{

  INITSTATUS (status, "LALHOUGHCreatePHMDVS", rcsid);
  ATTATCHSTATUSPTR (status);

  /* check input pars are ok */
  if (phmdVS == NULL) {
    ABORT (status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  }

  if (phmdVS->phmd != NULL) {
    ABORT (status, DRIVEHOUGHCOLOR_ENONULL, DRIVEHOUGHCOLOR_MSGENONULL);
  }

  if (length <= 0) {
    ABORT (status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  }

  if (nfSize <= 0) {
    ABORT (status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  }

  phmdVS->length  = length;
  phmdVS->nfSize  = nfSize;
  phmdVS->deltaF  = 0; /* initialization */

  phmdVS->phmd=(HOUGHphmd *)LALCalloc(1, length*nfSize*sizeof(HOUGHphmd));
  
  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}



void LALHOUGHCreatePHMDs(LALStatus           *status,
			 PHMDVectorSequence  *phmdVS,
			 UINT2               maxNBins, 
			 UINT2               maxNBorders,
			 UINT2               ySide)
{

  UINT4 j;

  INITSTATUS (status, "LALHOUGHCreatePHMDs", rcsid);
  ATTATCHSTATUSPTR (status);

  /* check input pars are ok */
  if (phmdVS == NULL) {
    ABORT (status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  }

  if (maxNBins <= 0) {
    ABORT (status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  }

  if (maxNBorders <= 0) {
    ABORT (status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  }

  if (ySide <= 0) {
    ABORT (status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  }

  for(j = 0; j < phmdVS->length * phmdVS->nfSize; j++){

    if ( phmdVS->phmd + j == NULL) {
      ABORT (status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
    }

    phmdVS->phmd[j].maxNBorders = maxNBorders;
    phmdVS->phmd[j].leftBorderP = (HOUGHBorder **)LALCalloc(maxNBorders, sizeof(HOUGHBorder *));
    phmdVS->phmd[j].rightBorderP = (HOUGHBorder **)LALCalloc(maxNBorders, sizeof(HOUGHBorder *));

    phmdVS->phmd[j].ySide = ySide;
    phmdVS->phmd[j].firstColumn = (UCHAR *)LALCalloc(ySide, sizeof(UCHAR));
  }

  
  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}


void LALHOUGHDestroyPHMDs(LALStatus           *status,
			  PHMDVectorSequence  *phmdVS)
{
  UINT4 j;

  INITSTATUS (status, "LALHOUGHDestroyPHMDs", rcsid);
  ATTATCHSTATUSPTR (status);

  for(j = 0; j < phmdVS->length * phmdVS->nfSize; j++){

    if (phmdVS->phmd + j) {
      LALFree( phmdVS->phmd[j].leftBorderP);
      LALFree( phmdVS->phmd[j].rightBorderP);
      LALFree( phmdVS->phmd[j].firstColumn);
    }

  }
  
  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}


void LALHOUGHCreateHT(LALStatus             *status,
		      HOUGHMapTotal         *ht,
		      UINT2                 xSide,
		      UINT2                 ySide)
{

  INITSTATUS (status, "LALHOUGHCreateHT", rcsid);
  ATTATCHSTATUSPTR (status);

  /* check input pars are ok */
  if (ht == NULL) {
    ABORT (status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  }

  if (xSide <= 0) {
    ABORT (status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  }

  if (ySide <= 0) {
    ABORT (status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  }


  ht->xSide = xSide;
  ht->ySide = ySide;

  ht->map = NULL;
  ht->map = (HoughTT *)LALCalloc(xSide*ySide, sizeof(HoughTT));

  ht->deltaF = 0; /* initialization */
  ht->mObsCoh = 0; /* initialization */
	
  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}



void LALHOUGHCreateFreqIndVector(LALStatus                 *status,
				 UINT8FrequencyIndexVector *freqInd,
				 UINT4                     length,
				 REAL8                     deltaF)
{

  INITSTATUS (status, "LALHOUGHCreateFreqIndVector", rcsid);
  ATTATCHSTATUSPTR (status);

  /* check input pars are ok */
  if (freqInd == NULL) {
    ABORT (status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  }

  if (length <= 0) {
    ABORT (status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  }

  if (deltaF <= 0) {
    ABORT (status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  }

  freqInd->deltaF = deltaF;
  freqInd->length = length;
  freqInd->data = NULL;
  freqInd->data =  ( UINT8 *)LALCalloc(1, freqInd->length*sizeof(UINT8));

  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}



/** Get Hough candidates as a toplist */
void GetToplistFromHoughmap(LALStatus *status,
			    toplist_t *list,
			    HOUGHMapTotal *ht,
			    HOUGHPatchGrid  *patch,
			    HOUGHDemodPar   *parDem,
			    REAL8 mean,
			    REAL8 sigma)
{
  REAL8UnitPolarCoor sourceLocation;
  REAL8 deltaF, f0, fdot;
  INT8 f0Bin;  
  INT4 i,j, xSide, ySide;
  FstatOutputEntry candidate;

  INITSTATUS( status, "GetToplistFromHoughMap", rcsid );
  ATTATCHSTATUSPTR (status);

  deltaF = ht->deltaF;
  f0Bin = ht->f0Bin;
  f0 = f0Bin * deltaF;

  fdot = ht->spinRes.data[0];

  xSide = ht->xSide;
  ySide = ht->ySide;

  candidate.Freq =  f0;
  candidate.f1dot = fdot;

  for (i = 0; i < ySide; i++)  {
      for (j = 0; j < xSide; j++) { 
	/* get sky location of pixel */
	TRY( LALStereo2SkyLocation (status->statusPtr, &sourceLocation, 
				    j, i, patch, parDem), status);
	
	candidate.Alpha = sourceLocation.alpha;
	candidate.Delta = sourceLocation.delta;
	candidate.Fstat =  (ht->map[i*xSide + j] - mean)/sigma;
	
	insert_into_fstat_toplist(list, candidate);
	
      } /* end loop over xSide */  
  } /* end loop over ySide */
  
  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* end hough toplist selection */


/* *********************************************************************************************/
/* Print some extra information: HoughMaps, HoughStatistics, sigma */
int OpenExtraInfoFiles(  CHAR          *fileMaps,
                         FILE          **fp1_ptr,
                         CHAR          *filehisto,
			 CHAR	       *dirname,
			 CHAR  	       *basename,
			 INT4          index)								
{
CHAR   filestats[ MAXFILENAMELENGTH ];
/*CHAR   fileMaps[ MAXFILENAMELENGTH ];*/
FILE    *fp1 = NULL;


/* --------------------------------------------------*/
/* Create directory fnameout/skypatch_$j using mkdir if required */
	if(CreateSkypatchDirs(filestats, dirname,index))
		return DRIVEHOUGHCOLOR_EFILE;
	
/* --------------------------------------------------*/	
/* Create the base filenames for the stats, histo and template files.*/
	strcat( filestats, basename);
	strcpy( filehisto, filestats);
	strcpy( fileMaps, filestats);
	
/* --------------------------------------------------*/
/* Create and open the stats file for writing */
	strcat( filestats, "stats");
	if( (fp1 = fopen(filestats,"w")) == NULL )
		{
			fprintf(stderr, "Unable to find file %s for writing\n", filestats);
			return DRIVEHOUGHCOLOR_EFILE;
		}
	setvbuf(fp1, (char * )NULL, _IOLBF, 0);
	*fp1_ptr = fp1;	
	return(0);

} /* end OpenExtraInfoFiles */
			

/* *********************************************************************************************/
/** Print some extra information: HoughMaps, HoughStatistics, sigma */
int PrintExtraInfo(      CHAR               *fileMaps,
                         FILE               **fp1_ptr,
                         INT4               iHmap,
                         HOUGHMapTotal      *ht,
		         REAL8UnitPolarCoor *sourceLocation,
		         HoughStats         *stats,
			 INT8               fBinSearch,
                         REAL8              deltaF)								
{
FILE *fp1 = *fp1_ptr;
		
/* --------------------------------------------------*/
/* Print results */

	/* Printing Hough Maps */
	if( PrintHmap2m_file( ht, fileMaps, iHmap ) ) return 5;
 
	/* Printing Statistics */
	fprintf(fp1, "%d %f %f %f %f %f %f %f %g\n",
			iHmap, sourceLocation->alpha, sourceLocation->delta,
			(REAL4)stats->maxCount, (REAL4)stats->minCount, stats->avgCount,stats->stdDev,
			(fBinSearch*deltaF), ht->spinRes.data[0]);
	
	return(0);

} /* end pirntExtraInfo */
					



/********************************************************************************/
/*                Computing the frequency path with no mismatch                 */
/********************************************************************************/
void ComputeFoft_NM(LALStatus   *status,
		 REAL8Vector          *foft,
                 HoughTemplate        *pulsarTemplate,
		 REAL8Vector          *timeDiffV,
		 REAL8Cart3CoorVector *velV){
  
  INT4   mObsCoh;
  REAL8   f0new, vcProdn, timeDiffN;
  REAL8   sourceDelta, sourceAlpha, cosDelta;
  INT4    j,i, nspin, factorialN; 
  REAL8Cart3Coor  sourceLocation;

  /* --------------------------------------------- */
  INITSTATUS (status, "ComputeFoft", rcsid);
  ATTATCHSTATUSPTR (status);
  
  /*   Make sure the arguments are not NULL: */
  ASSERT (foft,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (pulsarTemplate,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (timeDiffV,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (velV,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  
  ASSERT (foft->data,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (timeDiffV->data,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (velV->data,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);

  sourceDelta = pulsarTemplate->latitude;
  sourceAlpha = pulsarTemplate->longitude;
  cosDelta = cos(sourceDelta);  

  sourceLocation.x = cosDelta* cos(sourceAlpha);
  sourceLocation.y = cosDelta* sin(sourceAlpha);
  sourceLocation.z = sin(sourceDelta);
  
  mObsCoh = foft->length;    
  nspin = pulsarTemplate->spindown.length;
  
  for (j=0; j<mObsCoh; ++j){  /* loop for all different time stamps */
    vcProdn = velV->data[j].x * sourceLocation.x
      + velV->data[j].y * sourceLocation.y
      + velV->data[j].z * sourceLocation.z;
    f0new = pulsarTemplate->f0;
    factorialN = 1;
    timeDiffN = timeDiffV->data[j];
    
    for (i=0; i<nspin;++i){ /* loop for spin-down values */
      factorialN *=(i+1);
      f0new += pulsarTemplate->spindown.data[i]* timeDiffN / factorialN;
      timeDiffN *= timeDiffN;
    }
    foft->data[j] = f0new * (1.0 +vcProdn);
    /*foft->data[j] = floor(f0new*1800+0.5) * (1.0 +vcProdn)/1800;*/
  }    
  
  
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);

} /* end Computefoft */




/* *********************************************************************/

void SplitSFTs(LALStatus         *status,
	       REAL8Vector       *weightsV,
	       HoughParamsTest   *chi2Params){
  
  UINT4    j=0;           /* index of each block. It runs betwen 0 and p */ 
  UINT4   iSFT;       
  REAL8   *weights_ptr;  /* pointer to weightsV.data */
  REAL8   sumWeightpMax; /* Value of sumWeight we want to fix in each set of SFTs */
  UINT4   numberSFT;     /* Counter with the # of SFTs in each block */        
  UINT4   mObsCoh, p;
  REAL8   partialsumWeightp, partialsumWeightSquarep;
  
  /* --------------------------------------------- */
  INITSTATUS (status, "SplitSFTs", rcsid);
  ATTATCHSTATUSPTR (status);
  
  /*   Make sure the arguments are not NULL: */
  ASSERT (weightsV,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (chi2Params,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  
  ASSERT (weightsV->data,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (chi2Params->length,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (chi2Params->numberSFTp,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (chi2Params->sumWeight,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (chi2Params->sumWeightSquare,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (chi2Params->length < weightsV->length,  status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  
  mObsCoh = weightsV->length;    
  p = chi2Params->length;
  
  sumWeightpMax = (REAL8)(mObsCoh)/p;       /* Compute the value of the sumWeight we want to fix in each set of SFT's */
  weights_ptr = weightsV->data;    /* Make the pointer to point to the first position of the vector weightsV.data */
  
  iSFT = 0;
  for (j = 0; j < p; j++){

      partialsumWeightSquarep = 0;
      partialsumWeightp = 0;
    
      for(numberSFT = 0;(partialsumWeightp<sumWeightpMax)&&(iSFT<mObsCoh); numberSFT++, iSFT++){
 
	  partialsumWeightp += *weights_ptr;
	  partialsumWeightSquarep += (*weights_ptr)*(*weights_ptr);
	  weights_ptr++; 

      } /* loop over SFTs */
    
      ASSERT ( (UINT4)j < p, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);
  
      chi2Params->numberSFTp[j] = numberSFT;
      chi2Params->sumWeight[j] = partialsumWeightp;
      chi2Params->sumWeightSquare[j] = partialsumWeightSquarep;
    
  } /* loop over the p blocks of data */
       
  ASSERT ( iSFT == mObsCoh, status, DRIVEHOUGHCOLOR_EBAD, DRIVEHOUGHCOLOR_MSGEBAD);  

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);

}/* SplitSFTs */




/**********************************************************************/

void ComputeandPrintChi2 ( LALStatus                *status,
		           toplist_t                *tl,
		           REAL8Vector              *timeDiffV,
		           REAL8Cart3CoorVector     *velV,
			   INT4                     p,
			   REAL8                    alphaPeak,
			   MultiDetectorStateSeries *mdetStates,
			   REAL8Vector		    *weightsNoise,
			   UCHARPeakGramVector      *upgV /**< Expanded (UCHAR) peakgrams */)
{
    HoughTemplate    pulsarTemplate;
    UINT4            mObsCoh; 
    UINT4    k, i, j, ii, numberSFTp ;
    INT4    index;
    REAL8  sumWeightSquare, meanN, sigmaN;       
    REAL8   eta, numberCount;                 
    REAL8   nj, sumWeightj, sumWeightSquarej;
    FstatOutputEntry   readTopList;
    REAL8Vector  foft;
    REAL8Vector  weightsV;
    REAL8        timeBase;
    UINT4        sftFminBin ;
    FILE *fpChi2=NULL;
    REAL8 oldSig;

    /* Chi2Test parameters */
    HoughParamsTest chi2Params;
    REAL8  numberCountTotal;   /* Sum over all the numberCounts */
    REAL8  chi2;
    REAL8Vector numberCountV;  /* Vector with the number count of each block inside */

    /* --------------------------------------------- */
    INITSTATUS (status, "ComputeandPrintChi2", rcsid);
    ATTATCHSTATUSPTR (status);
  
    /*   Make sure the arguments are not NULL: */
    ASSERT (timeDiffV,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
    ASSERT (velV,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
    ASSERT (mdetStates,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
    ASSERT (weightsNoise,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
    ASSERT (upgV,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);

    ASSERT (weightsNoise->data,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
    ASSERT (timeDiffV->data,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
    ASSERT (velV->data,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);

    mObsCoh = velV->length;
    timeBase = upgV->upg[0].timeBase; 
    sftFminBin = upgV->upg[0].fminBinIndex;
    
    /* memory for weightsV */
    weightsV.length=mObsCoh;
    weightsV.data=(REAL8 *)LALCalloc(1, mObsCoh*sizeof(REAL8));
    
    /* memory for chi2Params */
    chi2Params.length = p;
    chi2Params.numberSFTp = NULL;
    chi2Params.sumWeight = NULL;
    chi2Params.sumWeightSquare = NULL;
    chi2Params.numberSFTp = (UINT4 *)LALMalloc( p*sizeof(UINT4));
    chi2Params.sumWeight = (REAL8 *)LALMalloc( p*sizeof(REAL8));
    chi2Params.sumWeightSquare = (REAL8 *)LALMalloc( p*sizeof(REAL8));

    /* memory for number Count Vector */
    numberCountV.length = p;
    numberCountV.data = NULL;
    numberCountV.data = (REAL8 *)LALMalloc( p*sizeof(REAL8));

    /* Memory for one spindown */
    pulsarTemplate.spindown.length = 1 ;
    pulsarTemplate.spindown.data = NULL;
    pulsarTemplate.spindown.data = (REAL8 *)LALMalloc(sizeof(REAL8));

    /* Memory for f(t) vector */
    foft.length = mObsCoh;
    foft.data = NULL;
    foft.data = (REAL8 *)LALMalloc( mObsCoh*sizeof(REAL8));

    /* Open file to write the toplist with 2 new columns: significance and chi2 */

    fpChi2 = fopen("hough_top.dat", "w");

    /* ----------------------------------------------------------------------------------*/
    /* Loop over all the elements in the TopList */

    for (i=0; i < tl->elems; i++){
    
	readTopList = *((FstatOutputEntry*)(tl->heap[i]));    

        /* Copy template parameters from the TopList */
	pulsarTemplate.f0= readTopList.Freq  ;
	pulsarTemplate.spindown.data[0] = readTopList.f1dot ;
	pulsarTemplate.latitude= readTopList.Delta ;
	pulsarTemplate.longitude= readTopList.Alpha ;
	oldSig = readTopList.Fstat;

	/* copy noise weights */
	memcpy(weightsV.data, weightsNoise->data, mObsCoh * sizeof(REAL8));

	/* calculate amplitude modulation weights */
	TRY(GetAMWeights( status->statusPtr, &weightsV, mdetStates, pulsarTemplate.longitude, pulsarTemplate.latitude), status);


	/**********************************************************************************/
	/* Split the SFTs into p blocks and calculate the number count in each block */  
	/**********************************************************************************/  

	/* compute mean and sigma for noise only */    
	/* first calculate the sum of the weights squared */
	sumWeightSquare = 0.0;
     
	for ( k = 0; k < (UINT4)mObsCoh; k++)
	    sumWeightSquare += weightsV.data[k] * weightsV.data[k];

        meanN = mObsCoh * alphaPeak;
	sigmaN = sqrt (sumWeightSquare * alphaPeak * (1.0 - alphaPeak));
	
	/* the received frequency as a function of time  */
	TRY( ComputeFoft_NM(status->statusPtr, &foft, &pulsarTemplate, timeDiffV, velV), status);   

	
	/* Split the SFTs into p blocs */
	TRY(SplitSFTs(status->statusPtr, &weightsV, &chi2Params), status);
	
            
	/* loop over SFT, generate peakgram and get number count */
	
	j=0;
   
	for (k=0 ; k<(UINT4)p ; k++ ){
	    
	    numberSFTp=chi2Params.numberSFTp[k];
	    numberCount = 0;
	    
	    for (ii=0 ; (ii < numberSFTp) ; ii++) {
		
		index = floor( foft.data[j]*timeBase - sftFminBin + 0.5); 
		
		numberCount += upgV->upg[j].data[index]*weightsV.data[j];
		
		j++;
		
	    } /* loop over SFTs */
	    
	    numberCountV.data[k]=numberCount;
	    
	}/* loop over blocks */
	
	
	/*-----------------------------*/  
	/* Chi2 Test */
	
	numberCountTotal=0;
	chi2=0;
	
	for(k=0; k<(UINT4)p ; k++){
	    numberCountTotal += numberCountV.data[k];
	}
	
	eta=numberCountTotal/mObsCoh;
	
	for(j=0 ; j<((UINT4)p) ; j++){
	    
	    nj=numberCountV.data[j];
	    sumWeightj=chi2Params.sumWeight[j];
	    sumWeightSquarej=chi2Params.sumWeightSquare[j];
	    
	    chi2 += (nj-sumWeightj*eta)*(nj-sumWeightj*eta)/(sumWeightSquarej*eta*(1-eta));
	}
	
	setvbuf(fpChi2, (char *)NULL, _IOLBF, 0);
	fprintf(fpChi2, "%g  %g %g  %g  %g %g  %g \n", pulsarTemplate.f0, pulsarTemplate.longitude, pulsarTemplate.latitude, pulsarTemplate.spindown.data[0], (numberCountTotal - meanN)/sigmaN, oldSig, chi2);
	
	/*-----------------------------*/
		
    } /* End of loop over top list elements */
    
    LALFree(pulsarTemplate.spindown.data);
    LALFree(foft.data);
    LALFree(numberCountV.data);
    LALFree(chi2Params.numberSFTp);
    LALFree(chi2Params.sumWeight);
    LALFree(chi2Params.sumWeightSquare);
    LALFree(weightsV.data);
    weightsV.data=NULL;
        
    fclose(fpChi2);

    
    DETATCHSTATUSPTR (status);
    
    /* normal exit */	
    RETURN (status);
    
} /* End of ComputeChi2 */



/************************************************************************************/
/** Loop over SFTs and set a threshold to get peakgrams.  SFTs must be normalized.  */
/** This function will create a vector with the uncompressed PeakGrams in it        */
/** (this is necesary if we want to compute the chi2 later )                        */
void GetPeakGramFromMultSFTVector_NondestroyPg1(LALStatus                   *status,
						HOUGHPeakGramVector         *out, /**< Output compressed peakgrams */
						UCHARPeakGramVector         *upgV, /**< Output uncompressed peakgrams */
						MultiSFTVector              *in,  /**< Input SFTs */
						REAL8                       thr   /**< Threshold on SFT power */)  
{
  SFTtype  *sft;
  INT4   nPeaks;
  UINT4  iIFO, iSFT, numsft, numifo, j, binsSFT; 
  
  INITSTATUS (status, "GetSFTNoiseWeights", rcsid);
  ATTATCHSTATUSPTR (status);

  numifo = in->length;
  binsSFT = in->data[0]->data->data->length;
 
  /* loop over sfts and select peaks */
  for ( j = 0, iIFO = 0; iIFO < numifo; iIFO++){
    
      numsft = in->data[iIFO]->length;
    
      for ( iSFT = 0; iSFT < numsft; iSFT++, j++) {
      
	  sft = in->data[iIFO]->data + iSFT;

          /* Store the expanded PeakGrams */ 
	  upgV->upg[j].length = binsSFT;
	  upgV->upg[j].data = NULL;
	  upgV->upg[j].data= (UCHAR *)LALCalloc( 1, binsSFT*sizeof(UCHAR));
     
	  TRY (SFTtoUCHARPeakGram( status->statusPtr, &(upgV->upg[j]), sft, thr), status);
      
	  nPeaks = upgV->upg[j].nPeaks;

          /* compress peakgram */      
	  out->pg[j].length = nPeaks;
	  out->pg[j].peak = NULL;
	  out->pg[j].peak = (INT4 *)LALCalloc( 1, nPeaks*sizeof(INT4));
      
	  TRY( LALUCHAR2HOUGHPeak( status->statusPtr, &(out->pg[j]), &(upgV->upg[j])), status );
      
      } /* loop over SFTs */
  } /* loop over IFOs */

 

  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}
