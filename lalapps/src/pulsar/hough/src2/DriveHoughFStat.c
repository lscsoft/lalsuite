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
 * \file DriveHoughFStat.c
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



#include"./DriveHoughFStat.h"

RCSID( "$Id$");

#define TRUE (1==1)
#define FALSE (1==0)

extern int lalDebugLevel;

BOOLEAN uvar_printMaps; /**< global variable for printing Hough maps */
BOOLEAN uvar_printStats; /**< global variable for calculating Hough map stats */


/* default values for input variables */
#define EARTHEPHEMERIS "../src/earth00-04.dat" /**< Default location of earth ephemeris */
#define SUNEPHEMERIS "../src/sun00-04.dat"   /**< Default location of sun ephemeris */
#define NSTACKS 10    /**< Default number of stacks */
#define IFO 2         /**<  Default detector, 1:GEO, 2:LLO, 3:LHO  */
#define BLOCKSRNGMED 101 /**< Default running median window size */
#define FSTART 255.0   /**< Default Start search frequency */
#define FBAND 0.001    /**< Default search band */
#define FDOT 0.0      /**< Default value of first spindown */
#define DFDOT 0.0   /**< Default range of first spindown parameter */
#define SKYREGION "(1,1),(1,1.5),(1.5,1.5),(1.5,1)" /**< default sky region to search over -- just a single point*/
#define NFSIZE  21    /**< Default size of hough cylinder of look up tables */
#define DTERMS 8     /**< Default number of dirichlet kernel terms for calculating Fstat */
#define MISMATCH 0.2 /**< Default for metric grid maximal mismatch value */
#define DALPHA 0.001 /**< Default resolution for isotropic or flat grids */
#define DDELTA 0.001 /**< Default resolution for isotropic or flat grids */
#define FSTATTHRESHOLD 2.6  /**< Default threshold on Fstatistic for peak selection */
#define NCAND1 5 /**< Default number of candidates to be followed up from first stage */
#define SFTDIRECTORY "/local_data/badkri/fakesfts/"  /**< Default directory containing sfts */
#define FNAMEOUT "./candidates"  /**< Default output file basename */

int main( int argc, char *argv[]) {

  /* initialize status */
  LALStatus status = blank_status;
  
  /* temp loop variables: generally k loops over stacks and j over SFTs in a stack*/
  INT4 j, k;

  /* in general any variable ending with 1 is for the 
     first stage, 2 for the second and so on */

  /* detectors and ephemeris */
  LALDetector detector1, detector2;
  EphemerisData *edat = NULL;

  /* timestamp vectors */
  LIGOTimeGPSVector midTstack1; 
  LIGOTimeGPSVector midTstack2; 
  LIGOTimeGPSVector *startTsft1=NULL; 
  LIGOTimeGPSVector *startTsft2=NULL; 
  REAL8 refTime;

  /* velocities and positions at midTstack */
  REAL8VectorSequence *velStack=NULL;
  REAL8VectorSequence *posStack=NULL;

  /* duration of each stack */
  REAL8 tObs1, tStart1, tEnd1;
  REAL8 tObs2, tStart2, tEnd2;
  REAL8 *tStack1, *tStack2, tStackAvg1, tStackAvg2;

  /* leap second for LALBarycenter */
  LALLeapSecFormatAndAcc lsfas = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
  INT4 tmpLeap; 

  /* sft related stuff */
  /* input sft vectors */
  SFTVector *inputSFTVec1=NULL;
  SFTVector *inputSFTVec2=NULL;

  /* sequence of sft vectors -- one for each stack */
  SFTVectorSequence stackSFTs1;
  SFTVectorSequence stackSFTs2;

  /* number of SFTs in each stack and total number of SFTs */
  INT4 *mCohSft1, *mCohSft2, nSFTs1, nSFTs2;

  /* number of stacks -- not necessarily same as uvar_nStacks! */
  INT4 nStacks1, nStacks2;
  INT4 sftlength1, sftlength2; /* number of bins in each sft */
  REAL8 deltaF1, deltaF2, timebase1, timebase2; /* frequency resolution of SFTs */
  INT8 sftFminBin1, sftFminBin2; /* first sft bin index */
  INT8 fHoughBinIni, fHoughBinFin, binsHough; /* frequency bins of start and end search frequencies */
  REAL8 deltaFstack1, deltaFstack2; /* frequency resolution of Fstat calculation */

  /* LALdemod related stuff */
  REAL8FrequencySeriesVector FstatVect1, FstatVect2; /* Fstatistic vectors for each stack */
  FstatStackParams FstatPar1, FstatPar2;

  /* hough variables */
  HOUGHPeakGramVector pgV;
  HoughParams houghPar;
  HoughCandidates houghCand1;

  /* fstat candidate structure */
  toplist_t *fstatToplist=NULL;

  /* template and grid variables */
  DopplerScanInit scanInit1, scanInit2;   /* init-structure for DopperScanner */
  DopplerScanState thisScan1 = empty_DopplerScanState; /* current state of the Doppler-scan */
  DopplerScanState thisScan2 = empty_DopplerScanState; /* current state of the Doppler-scan */
  DopplerPosition dopplerpos1, dopplerpos2;		/* current search-parameters */
  SkyPosition thisPoint1, thisPoint2;

  /* variables for logging */
  CHAR *fnamelog=NULL;
  FILE *fpLog=NULL;
  CHAR *logstr=NULL; 

  /* output candidate files and file pointers */
  CHAR *fnameFstatCand=NULL;
  CHAR *fnameHoughCand=NULL;
  FILE *fpFstat=NULL;
  FILE *fpHough=NULL;
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
  REAL8 uvar_fdot; /* first spindown value */
  REAL8 uvar_fdotBand; /* range of first spindown parameter */
  REAL8 uvar_fStart, uvar_fBand;
  REAL8 uvar_FstatPeakSelect; /* threshold of Fstat to select peaks */
  REAL8 uvar_houghThr; /* threshold on hough number count to select candidates */
  REAL8 uvar_fstatThr; /* threshold for selecting candidates from Fstat vector */
  REAL8 uvar_metricMismatch1;
  REAL8 uvar_metricMismatch2;
  REAL8 uvar_refTime;
  INT4 uvar_nCand1; /* number of candidates to be followed up from first stage */
  INT4 uvar_nCand2; /* number of candidates from second stage */
  INT4 uvar_ifo1, uvar_ifo2; 
  INT4 uvar_blocksRngMed;
  INT4 uvar_nStacks1, uvar_nStacks2;
  INT4 uvar_Dterms;
  INT4 uvar_SSBprecision;
  INT4 uvar_nfSize;
  INT4 uvar_gridType;
  INT4 uvar_metricType;
  INT4 uvar_reallocBlock;
  CHAR *uvar_earthEphemeris=NULL;
  CHAR *uvar_sunEphemeris=NULL;
  CHAR *uvar_sftDir1=NULL;
  CHAR *uvar_sftDir2=NULL;
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

  /* now set the defaults */
  uvar_help = FALSE;
  uvar_log = FALSE;
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
  uvar_fdot = FDOT;
  uvar_fdotBand = DFDOT;
  uvar_fStart = FSTART;
  uvar_fBand = FBAND;
  uvar_ifo1 = IFO;
  uvar_ifo2 = IFO;
  uvar_blocksRngMed = BLOCKSRNGMED;
  uvar_nfSize = NFSIZE;
  uvar_FstatPeakSelect = FSTATTHRESHOLD;
  uvar_nCand1 = uvar_nCand2 = NCAND1;
  uvar_houghThr = 0;
  uvar_fstatThr = FSTATTHRESHOLD;
  uvar_SSBprecision = SSBPREC_RELATIVISTIC;
  uvar_metricType = LAL_PMETRIC_COH_PTOLE_ANALYTIC;
  uvar_gridType = GRID_METRIC;
  uvar_metricMismatch1 = uvar_metricMismatch2 = MISMATCH;
  uvar_reallocBlock = 5000;

  uvar_skyGridFile = NULL;

  uvar_skyRegion = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_skyRegion, SKYREGION);

  uvar_earthEphemeris = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_earthEphemeris, EARTHEPHEMERIS);

  uvar_sunEphemeris = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_sunEphemeris, SUNEPHEMERIS);

  uvar_sftDir1 = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_sftDir1, SFTDIRECTORY);

  /* do not set default for sftDir2 -- use only if user specifies */
  /*   uvar_sftDir2 = (CHAR *)LALMalloc(512*sizeof(CHAR)); */

  uvar_fnameout = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_fnameout, FNAMEOUT);

  /* register user input variables */
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "help",            'h', UVAR_HELP,     "Print this message",                                &uvar_help),            &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "log",              0,  UVAR_OPTIONAL, "Write log file",                                    &uvar_log),             &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "chkPoint",         0,  UVAR_OPTIONAL, "For checkpointing (disabled)",                      &uvar_chkPoint),        &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "followUp",         0,  UVAR_OPTIONAL, "Follow up 1st stage candidates?",                   &uvar_followUp),        &status);  
  LAL_CALL( LALRegisterINTUserVar(    &status, "ifo1",            'i', UVAR_OPTIONAL, "Detector GEO(1) LLO(2) LHO(3)",                     &uvar_ifo1 ),           &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "ifo2",             0,  UVAR_OPTIONAL, "Detector for follow up stage",                      &uvar_ifo2 ),           &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftDir1",          0,  UVAR_OPTIONAL, "SFT Directory for 1st stage",                       &uvar_sftDir1),         &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftDir2",          0,  UVAR_OPTIONAL, "SFT Directory for 2nd stage",                       &uvar_sftDir2),         &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "skyRegion",	   0,  UVAR_OPTIONAL, "Specify sky-region by polygon (or use 'allsky')",   &uvar_skyRegion),       &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fStart",          'f', UVAR_OPTIONAL, "Start search frequency",                            &uvar_fStart),          &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fBand",           'b', UVAR_OPTIONAL, "Search frequency band",                             &uvar_fBand),           &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fdot",             0,  UVAR_OPTIONAL, "Spindown parameter",                                &uvar_fdot),            &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fdotBand",         0,  UVAR_OPTIONAL, "Range of spindown parameter",                       &uvar_fdotBand),        &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "nStacks1",        'N', UVAR_OPTIONAL, "Number of 1st stage stacks",                        &uvar_nStacks1 ),       &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "nStacks2",         0,  UVAR_OPTIONAL, "Number of 2nd stage stacks",                        &uvar_nStacks2 ),       &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "metricMismatch1",  0,  UVAR_OPTIONAL, "Maximal allowed 1st stage metric mismatch",         &uvar_metricMismatch1), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "metricMismatch2",  0,  UVAR_OPTIONAL, "Maximal allowed 2nd stage metric mismatch",         &uvar_metricMismatch2), &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "gridType",         0,  UVAR_OPTIONAL, "0=flat,1=isotropic,2=metric,3=file",                &uvar_gridType),        &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "metricType",       0,  UVAR_OPTIONAL, "0=none,1=Ptole-analytic,2=Ptole-numeric,3=exact",   &uvar_metricType),      &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "skyGridFile",	   0,  UVAR_OPTIONAL, "Load sky-grid from this file",                      &uvar_skyGridFile),     &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dAlpha",           0,  UVAR_OPTIONAL, "Resolution for flat or isotropic grids",            &uvar_dAlpha),          &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dDelta",           0,  UVAR_OPTIONAL, "Resolution for flat or isotropic grids",            &uvar_dDelta),          &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "fnameout",        'o', UVAR_OPTIONAL, "Output basefileneme",                               &uvar_fnameout),        &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "FstatPeakSelect",  0,  UVAR_OPTIONAL, "Threshold on Fstatistic for peak selection",        &uvar_FstatPeakSelect), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "nCand1",           0,  UVAR_OPTIONAL, "No. of 1st stage candidates to be followed up",     &uvar_nCand1),          &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "nCand2",           0,  UVAR_OPTIONAL, "No. of 2nd stage candidates to be followed up",     &uvar_nCand2),          &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "houghThr",         0,  UVAR_OPTIONAL, "Hough number count threshold (default --nCand1)",   &uvar_houghThr),        &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fstatThr",         0,  UVAR_OPTIONAL, "Fstat threshold (default --nCand2)",                &uvar_fstatThr),        &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printCand1",       0,  UVAR_OPTIONAL, "Print 1st stage candidates",                        &uvar_printCand1),      &status);  
  LAL_CALL( LALRegisterREALUserVar(   &status, "refTime",          0,  UVAR_OPTIONAL, "Ref. time for pulsar pars (default = start time)",  &uvar_refTime),         &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed",    'w', UVAR_OPTIONAL, "RngMed block size",                                 &uvar_blocksRngMed),    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "earthEphemeris",  'E', UVAR_OPTIONAL, "Earth Ephemeris file",                              &uvar_earthEphemeris),  &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sunEphemeris",    'S', UVAR_OPTIONAL, "Sun Ephemeris file",                                &uvar_sunEphemeris),    &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "reallocBlock",     0,  UVAR_DEVELOPER,"Blocks to realloc for Fstat output if necessary",   &uvar_reallocBlock),    &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "SSBprecision",	   0,  UVAR_DEVELOPER,"Precision for SSB transform.",                      &uvar_SSBprecision),    &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "nfSize",           0,  UVAR_DEVELOPER,"No.of LUTs to keep in memory",                      &uvar_nfSize),          &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printMaps",        0,  UVAR_DEVELOPER,"Print Hough maps",                                  &uvar_printMaps),       &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printStats",       0,  UVAR_DEVELOPER,"Print Hough map statistics",                        &uvar_printStats),      &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printFstat1",      0,  UVAR_DEVELOPER,"Print 1st stage Fstat vectors ",                    &uvar_printFstat1),     &status);  
  LAL_CALL( LALRegisterINTUserVar(    &status, "Dterms",           0,  UVAR_DEVELOPER,"For Dirichlet Kernel approx.",                      &uvar_Dterms ),         &status);

  /* read all command line variables */
  LAL_CALL( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 

  /* some basic sanity checks on user vars */
  if ( uvar_nStacks1 < 1) {
    fprintf(stderr, "Invalid number of stacks\n");
    exit( DRIVEHOUGHFSTAT_EBAD );
  }

  if ( uvar_nStacks2 < 1) {
    fprintf(stderr, "Invalid number of stacks\n");
    exit( DRIVEHOUGHFSTAT_EBAD );
  }

  if ( uvar_blocksRngMed < 1 ) {
    fprintf(stderr, "Invalid Running Median block size\n");
    exit( DRIVEHOUGHFSTAT_EBAD );
  }

  if ( uvar_FstatPeakSelect < 0 ) {
    fprintf(stderr, "Invalid value of Fstatistic threshold\n");
    exit( DRIVEHOUGHFSTAT_EBAD );
  }


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
      
      fprintf( fpLog, "## Log file for MCInjectHoughS2\n\n");
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
  

  /* set detector */
  if (uvar_ifo1 == 1) detector1 = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if (uvar_ifo1 == 2) detector1 = lalCachedDetectors[LALDetectorIndexLLODIFF];
  if (uvar_ifo1 == 3) detector1 = lalCachedDetectors[LALDetectorIndexLHODIFF];

  if (uvar_ifo2 == 1) detector2 = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if (uvar_ifo2 == 2) detector2 = lalCachedDetectors[LALDetectorIndexLLODIFF];
  if (uvar_ifo2 == 3) detector2 = lalCachedDetectors[LALDetectorIndexLHODIFF];



  /* create output Hough file for writing if requested by user */
  if ( uvar_printCand1 )
    {
      fnameHoughCand = (CHAR *)LALMalloc( 512*sizeof(CHAR));
      strcpy(fnameHoughCand, uvar_fnameout);
      strcat(fnameHoughCand, "_hough.txt");
      if (!(fpHough = fopen(fnameHoughCand, "w"))) 
	{
	  fprintf ( stderr, "Unable to open Hough file '%s' for writing.\n", fnameHoughCand);
	  return DRIVEHOUGHFSTAT_EFILE;
	}
    }


  /*------------- read sfts and set up sft timestamp vector ----------*/
  /* for first stage */
  {

    /* new SFT I/O data types */
    SFTCatalog *catalog1=NULL;
    SFTConstraints *constraints1=NULL;

    CHAR *tempDir;
    REAL8 doppWings, fMin, fMax;

    /* get sft catalog */
    tempDir = (CHAR *)LALMalloc(512*sizeof(CHAR));
    strcpy(tempDir, uvar_sftDir1);
    strcat(tempDir, "/*SFT*.*");
    LAL_CALL( LALSFTdataFind( &status, &catalog1, tempDir, constraints1), &status);

    /* set some sft parameters */
    nSFTs1 = catalog1->length;
    deltaF1 = catalog1->data->header.deltaF;
    timebase1 = 1.0/deltaF1;

    /* set wings of sfts to be read */
    /* the wings must be enough for the Doppler shift and extra bins
       for the running median block size and Dterms for Fstat calculation.
       In addition, it must also include wings for the spindown correcting 
       for the reference time -- this is not done yet */
    /* calculate Doppler wings at the highest frequency */
    doppWings = (uvar_fStart + uvar_fBand) * VTOT;    
    fMin = uvar_fStart - doppWings - (uvar_blocksRngMed + uvar_Dterms)*deltaF1; 
    fMax = uvar_fStart + uvar_fBand + (doppWings + uvar_blocksRngMed + uvar_Dterms)*deltaF1; 

    /* get SFT timestamps */
    LAL_CALL( LALSFTtimestampsFromCatalog(  &status, &startTsft1, catalog1 ), &status);  	

    /* calculate start and end times and tobs */
    LAL_CALL( LALGPStoFloat ( &status, &tStart1, startTsft1->data), &status);
    LAL_CALL( LALGPStoFloat ( &status, &tEnd1, startTsft1->data + startTsft1->length - 1), &status);
    tEnd1 += timebase1;
    tObs1 = tEnd1 - tStart1;

    /* set reference time for pular parameters */
    if ( LALUserVarWasSet(&uvar_refTime)) 
      refTime = uvar_refTime;
    else
      refTime = tStart1;

    /* read the sfts */
    LAL_CALL( LALLoadSFTs ( &status, &inputSFTVec1, catalog1, fMin, fMax), &status);

    /* set other sft parameters */
    sftlength1 = inputSFTVec1->data->data->length;
    sftFminBin1 = floor( timebase1 * inputSFTVec1->data->f0 + 0.5);	      

    /* normalize sfts */
    LAL_CALL( LALNormalizeSFTVect (&status, inputSFTVec1, uvar_blocksRngMed, 0), &status);

    /* free memory */
    LALFree(tempDir);
    LAL_CALL( LALDestroySFTCatalog( &status, &catalog1 ), &status);  	

  } /* end of 1st stage sft reading block */


  /*------------- set up 1st stage stacks -----------------*/

  if (uvar_nStacks1 > nSFTs1) {
    fprintf(stderr, "invalid number of stacks...exiting\n");
    exit(1);
  }

  /* set up the stacks for the first stage */
  /* if sfts are to be split evenly among stacks */
  LAL_CALL( SetUpStacks1( &status, &stackSFTs1, inputSFTVec1, uvar_nStacks1), &status);
  /* if time is to be split evenly between stacks */
  /* LAL_CALL( SetUpStacks2( &status, &stackSFTs1, inputSFTVec1, &startTsft1, uvar_nStacks1), &status); */

  /* set number of stacks -- may be different from uvar_nStacks1! */
  nStacks1 = stackSFTs1.length;

  /* set up vector of stack durations */
  tStack1 = NULL;
  tStack1 = (REAL8 *)LALMalloc( nStacks1 * sizeof(REAL8));

  /* set up vector of number of sfts in each stack */
  mCohSft1 = NULL;
  mCohSft1 = (INT4 *)LALMalloc( nStacks1 * sizeof(INT4));

  /* set up vector containing mid times of stacks */    
  midTstack1.length = nStacks1;
  midTstack1.data = (LIGOTimeGPS *)LALMalloc( nStacks1 * sizeof(LIGOTimeGPS));

  for (k=0; k<nStacks1; k++) {
    LIGOTimeGPS tempT1, tempT2;
    INT4 tempInt;
    REAL8 tempF1, tempF2, tempMid;

    /* number of sfts in stack */
    tempInt = stackSFTs1.data[k].length;
    mCohSft1[k] = tempInt;

    /* duration of each stack */
    tempT1 = stackSFTs1.data[k].data[0].epoch;
    tempT2 = stackSFTs1.data[k].data[tempInt-1].epoch;
    LAL_CALL ( LALGPStoFloat ( &status, &tempF1, &tempT1), &status);
    LAL_CALL ( LALGPStoFloat ( &status, &tempF2, &tempT2), &status);
    tStack1[k] = tempF2 + timebase1 - tempF1;
    
    /* mid timestamp of each stack */
    tempMid = tempF1 + 0.5 * tStack1[k];
    LAL_CALL ( LALFloatToGPS ( & status, midTstack1.data + k, &tempMid), &status);
  }

  /* use stacks info to calculate search frequencies */
  /* Fstat is calculated at the frequency resolutions of the stacks
     here the stacks may be of different durations so we take the average
     This is valid if the stack durations are not very different which 
     will, hopefully, be true */
  tStackAvg1 = 0.0;
  for (k=0; k<nStacks1; k++)
    tStackAvg1 += tStack1[k];
  tStackAvg1 /= nStacks1;
  deltaFstack1 = 1.0/tStackAvg1;


  
  /*------------ calculate velocity and position for each 1st stage stack ------------*/
  /* setting of ephemeris info */ 
  edat = (EphemerisData *)LALMalloc(sizeof(EphemerisData));
  (*edat).ephiles.earthEphemeris = uvar_earthEphemeris;
  (*edat).ephiles.sunEphemeris = uvar_sunEphemeris;

  /* create velocity and position vectors */
  {
    CreateVectorSequenceIn createPar;
    createPar.length = nStacks1; /* number of vectors */
    createPar.vectorLength = 3; /* length of each vector */
    LAL_CALL( LALDCreateVectorSequence ( &status,  &velStack, &createPar), &status);
    LAL_CALL( LALDCreateVectorSequence ( &status,  &posStack, &createPar), &status);
  }
    
  /* Leap seconds for the first timestamp */   
  LAL_CALL( LALLeapSecs(&status, &tmpLeap, midTstack1.data, &lsfas), &status);
  (*edat).leap = (INT2)tmpLeap;
  
  /* read in ephemeris data */
  LAL_CALL( LALInitBarycenter( &status, edat), &status);
  
  /* calculate detector velocity and position at mid time of stacks*/
  /* maybe calculate average over stack as well? */
  for (k=0; k<nStacks1; k++){
    LAL_CALL (LALDetectorVel ( &status, velStack->data + 3*k, midTstack1.data + k, detector1, edat), &status);
    LAL_CALL (LALDetectorPos ( &status, posStack->data + 3*k, midTstack1.data + k, detector1, edat), &status);
  }



  /*------------------ read sfts and set up stacks for follow up stage -----------------------*/
  /* check if user requested a follow up stage*/
  if ( uvar_followUp ) {
    
    /* a friendly warning */
    if (lalDebugLevel)
      fprintf(stdout, "Make sure IFOs are set correctly!\n");

    /* read sfts for second stage
       currently reads entire frequency band because it is probably inefficient to
       read sfts separately for each candidate, and range of frequencies
       to be followed up is probably very close to the full range.  However, this
       might have to be changed if there are a very few candidates to be followed up.
    */
    if ( LALUserVarWasSet(&uvar_sftDir2))     
      {
	/* new SFT I/O data types */
	SFTCatalog *catalog2=NULL;
	SFTConstraints *constraints2=NULL;
	
	CHAR *tempDir;
	REAL8 doppWings, fMin, fMax;
	
	/* get sft catalog */
	tempDir = (CHAR *)LALMalloc(512*sizeof(CHAR));
	strcpy(tempDir, uvar_sftDir2);
	strcat(tempDir, "/*SFT*.*");
	LAL_CALL( LALSFTdataFind( &status, &catalog2, tempDir, constraints2), &status);
	
	/* set some sft parameters */
	nSFTs2 = catalog2->length;
	deltaF2 = catalog2->data->header.deltaF;
	timebase2 = 1.0/deltaF2;
	
	/* set wings of sfts to be read */
	/* the wings must be enough for the Doppler shift and extra bins
	   for the running median block size and Dterms for Fstat calculation.
	   In addition, it must also include wings for the spindown correcting 
	   for the reference time -- this is not done yet */
	/* calculate Doppler wings at the highest frequency */
	/* this was already done for the first stage -- probably shouldn't repeat it */ 
	doppWings = (uvar_fStart + uvar_fBand) * VTOT;    
	fMin = uvar_fStart - doppWings - (uvar_blocksRngMed + uvar_Dterms)*deltaF2; 
	fMax = uvar_fStart + uvar_fBand + (doppWings + uvar_blocksRngMed + uvar_Dterms)*deltaF2; 
	
	/* get SFT timestamps */
	LAL_CALL( LALSFTtimestampsFromCatalog(  &status, &startTsft2, catalog2 ), &status);  	
	
	/* calculate start and end times and tobs */
	LAL_CALL( LALGPStoFloat ( &status, &tStart2, startTsft2->data), &status);
	LAL_CALL( LALGPStoFloat ( &status, &tEnd2, startTsft2->data + startTsft2->length - 1), &status);
	tEnd2 += timebase2;
	tObs2 = tEnd2 - tStart2;
	
	/* read the sfts */
	LAL_CALL( LALLoadSFTs ( &status, &inputSFTVec2, catalog2, fMin, fMax), &status);
	
	/* set other sft parameters */
	sftlength2 = inputSFTVec2->data->data->length;
	sftFminBin2 = floor( timebase2 * inputSFTVec2->data->f0 + 0.5);	      
	
	/* normalize sfts */
	LAL_CALL( LALNormalizeSFTVect (&status, inputSFTVec2, uvar_blocksRngMed, 0), &status);
      
	/* free memory */
	LALFree(tempDir);
	LAL_CALL( LALDestroySFTCatalog( &status, &catalog2 ), &status);  	
	
      } /* end sftdir2 if statement */
    else 
      {
	/* reuse 1st stage sfts */
	
	/* set some sft parameters */
	nSFTs2 = nSFTs1;
	deltaF2 = deltaF1;
	timebase2 = 1.0/deltaF2;

	/* time stamps are the same */
	startTsft2 = startTsft1;	
		
	/* calculate start and end times and tobs */
	LAL_CALL( LALGPStoFloat ( &status, &tStart2, startTsft2->data), &status);
	LAL_CALL( LALGPStoFloat ( &status, &tEnd2, startTsft2->data + startTsft2->length - 1), &status);
	tEnd2 += timebase2;
	tObs2 = tEnd2 - tStart2;
	
	/* set sftvec2 equal to 1st stage sft vector -- these are already normalized*/
	inputSFTVec2 = inputSFTVec1;
	
	/* set other sft parameters */
	sftlength2 = inputSFTVec2->data->data->length;
	sftFminBin2 = floor( timebase2 * inputSFTVec2->data->f0 + 0.5);	      
	

      } /* end 2nd stage sft reading */ 

    /* set up stacks for follow up stage */
    /* there is just one stack now */

    LAL_CALL( SetUpStacks1( &status, &stackSFTs2, inputSFTVec2, uvar_nStacks2), &status);
    nStacks2 = stackSFTs2.length;    

    /* set up vector of stack durations */
    tStack2 = NULL;
    tStack2 = (REAL8 *)LALMalloc( nStacks2 * sizeof(REAL8));
    
    /* set up vector of number of sfts in each stack */
    mCohSft2 = NULL;
    mCohSft2 = (INT4 *)LALMalloc( nStacks2 * sizeof(INT4));

    /* set up vector containing mid times of stacks */    
    midTstack2.length = nStacks2;
    midTstack2.data = (LIGOTimeGPS *)LALMalloc( nStacks2 * sizeof(LIGOTimeGPS));

    for (k=0; k<nStacks2; k++) {

      LIGOTimeGPS tempT1, tempT2;
      INT4 tempInt;
      REAL8 tempF1, tempF2, tempMid;

      /* number of sfts in stack */
      tempInt = stackSFTs2.data[k].length;
      mCohSft2[k] = tempInt;
      
      /* duration of each stack */
      tempT1 = stackSFTs2.data[k].data[0].epoch;
      tempT2 = stackSFTs2.data[k].data[tempInt-1].epoch;
      LAL_CALL ( LALGPStoFloat ( &status, &tempF1, &tempT1), &status);
      LAL_CALL ( LALGPStoFloat ( &status, &tempF2, &tempT2), &status);
      tStack2[k] = tempF2 + timebase2 - tempF1;
      
      /* mid timestamp of each stack */
      tempMid = tempF1 + 0.5 * tStack2[k];
      LAL_CALL ( LALFloatToGPS ( & status, midTstack2.data + k, &tempMid), &status);
    } /* end of loop for setting midpoints of 2nd stage stacks */
    
    /* use stacks info to calculate search frequencies */
    /* Fstat is calculated at the frequency resolutions of the stacks
       here the stacks may be of different durations so we take the average
       This is valid if the stack durations are not very different which 
       will, hopefully, be true */
    tStackAvg2 = 0.0;
    for (k=0; k<nStacks2; k++)
      tStackAvg2 += tStack2[k];
    tStackAvg2 /= nStacks2;
    deltaFstack2 = 1.0/tStackAvg2;


  } /* end of sft reading and setting up stacks for follow up stage */


  /* set up parameters for first stage Fstat calculation */
  FstatPar1.nStacks = nStacks1;
  FstatPar1.tsStack = &midTstack1;
  FstatPar1.tStackAvg = tStackAvg1;
  FstatPar1.fBand = uvar_fBand;
  FstatPar1.fStart = uvar_fStart;
  FstatPar1.nfSizeCylinder = uvar_nfSize;
  FstatPar1.mCohSft = mCohSft1;
  FstatPar1.refTime = refTime;
  FstatPar1.SSBprecision = uvar_SSBprecision;
  FstatPar1.Dterms = uvar_Dterms;
  FstatPar1.detector = detector1;
  FstatPar1.edat = edat;
  FstatPar1.ts = startTsft1;
  FstatPar1.fdot = NULL;
  LAL_CALL ( LALDCreateVector( &status, &(FstatPar1.fdot), 1), &status);

  
  /* set up the Hough parameters */
  /* start and end bin for calculating hough map */
  /* these are just what the user specified */
  binsHough = floor( tStackAvg1 * uvar_fBand );
  fHoughBinIni = floor( tStackAvg1 * uvar_fStart + 0.5);
  fHoughBinFin = fHoughBinIni + binsHough - 1;
  if ( ! LALUserVarWasSet(&uvar_houghThr))
    uvar_houghThr = 0.65 * nStacks1;
  houghPar.outBaseName = uvar_fnameout;
  houghPar.houghThr = uvar_houghThr;
  houghPar.tStart = tStart1;
  houghPar.fBinIni = fHoughBinIni;
  houghPar.fBinFin = fHoughBinFin;
  houghPar.nfSizeCylinder = uvar_nfSize;
  houghPar.detector = detector1;
  houghPar.ts = &midTstack1;
  houghPar.vel = velStack;
  houghPar.pos = posStack;
  houghPar.fdot = NULL;
  LAL_CALL ( LALDCreateVector( &status, &(houghPar.fdot), 1), &status);


  /* allocate memory for Hough candidates */
  houghCand1.length = uvar_nCand1;
  houghCand1.nCandidates = 0; /* initialization */
  houghCand1.minSigIndex = 0;
  houghCand1.list = (HoughList *)LALMalloc( houghCand1.length * sizeof(HoughList));


  /*-----------Create template grid for first stage ---------------*/
  /* prepare initialization of DopplerScanner to step through paramter space */
  scanInit1.dAlpha = uvar_dAlpha;
  scanInit1.dDelta = uvar_dDelta;
  scanInit1.gridType = uvar_gridType;
  scanInit1.metricType = uvar_metricType;
  scanInit1.metricMismatch = uvar_metricMismatch1;
  scanInit1.projectMetric = TRUE;
  scanInit1.obsDuration = tStackAvg1;
  scanInit1.obsBegin = midTstack1.data[ nStacks1/2 ];
  scanInit1.Detector = &detector1;
  scanInit1.ephemeris = edat;		/* used by Ephemeris-based metric */
  scanInit1.skyGridFile = uvar_skyGridFile;
  scanInit1.searchRegion.Freq = uvar_fStart;
  scanInit1.searchRegion.FreqBand = uvar_fBand;
  scanInit1.searchRegion.f1dot = uvar_fdot;
  scanInit1.searchRegion.f1dotBand = uvar_fdotBand;
  scanInit1.searchRegion.skyRegionString = (CHAR*)LALCalloc(1, strlen(uvar_skyRegion)+1);
  strcpy (scanInit1.searchRegion.skyRegionString, uvar_skyRegion);


  /* parameters for 2nd stage */
  /* set up parameters for second stage Fstat calculation */
  /* same as for FstatPar1 except we don't set frequency
     and frequency band -- these depend on candidates to be followed up */
  if ( uvar_followUp ) 
    {
      /* fill up parameters for Fstat calculation */
      FstatPar2.nStacks = nStacks2;
      FstatPar2.tsStack = &midTstack2;
      FstatPar2.tStackAvg = tStackAvg2;
      FstatPar2.nfSizeCylinder = uvar_nfSize;
      FstatPar2.mCohSft = mCohSft2;
      FstatPar2.refTime = refTime;
      FstatPar2.SSBprecision = uvar_SSBprecision;
      FstatPar2.Dterms = uvar_Dterms;
      FstatPar2.detector = detector2;
      FstatPar2.edat = edat;
      FstatPar2.ts = startTsft2;
      FstatPar2.fdot = NULL;
      LAL_CALL ( LALDCreateVector( &status, &(FstatPar2.fdot), 1), &status);


      /* allocate memory for Fstat candidates */
      if ( LALUserVarWasSet(&uvar_fstatThr))
	  create_toplist(&fstatToplist, uvar_reallocBlock); 
      else if ( uvar_nCand2 > 0 )
	create_toplist(&fstatToplist, uvar_nCand2); 


     /* prepare initialization of DopplerScanner to step through paramter space */
      scanInit2.dAlpha = uvar_dAlpha;
      scanInit2.dDelta = uvar_dDelta;
      scanInit2.gridType = uvar_gridType;
      scanInit2.metricType = uvar_metricType;
      scanInit2.metricMismatch = uvar_metricMismatch2;
      scanInit2.projectMetric = TRUE;
      scanInit2.obsDuration = tStackAvg2;
      scanInit2.obsBegin = midTstack2.data[ nStacks1/2 ];
      scanInit2.Detector = &detector2;
      scanInit2.ephemeris = edat;		/* used by Ephemeris-based metric */
      scanInit2.skyGridFile = uvar_skyGridFile;
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
	  return DRIVEHOUGHFSTAT_EFILE;
	}
    } 

  if ( uvar_printFstat1 )
    {
      if ( !(fpFstat1 = fopen( "./fstatvec1.out", "w")))
	{
	  fprintf ( stderr, "Unable to open Fstat file fstatvec1.out for writing.\n");
	  return DRIVEHOUGHFSTAT_EFILE;
	}
    }

 
  /* initialize skygrid  */  
  LAL_CALL ( InitDopplerScan ( &status, &thisScan1, &scanInit1), &status); 


  /* loop over skygrid points */
  while(1)
    {
      UINT4 ifdot, nfdot;  /* counter and number of spindown values */
      REAL8 dfDot;  /* resolution in spindown */
      
      LAL_CALL (NextDopplerPos( &status, &dopplerpos1, &thisScan1 ), &status);
      if (thisScan1.state == STATE_FINISHED) /* scanned all DopplerPositions yet? */
	break;
      
      /*------------- calculate F statistic for each stack --------------*/
      
      /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
      thisPoint1.longitude = dopplerpos1.Alpha;
      thisPoint1.latitude = dopplerpos1.Delta;
      thisPoint1.system = COORDINATESYSTEM_EQUATORIAL;
      LAL_CALL (LALNormalizeSkyPosition(&status, &thisPoint1, &thisPoint1), &status);
      
      /* set sky template points */
      FstatPar1.alpha = thisPoint1.longitude;
      FstatPar1.delta = thisPoint1.latitude;
      
      /* number of fdot values */
      dfDot = thisScan1.df1dot;
      nfdot = (UINT4)(uvar_fdotBand / dfDot + 0.5) + 1; 
      
      /* loop over fdot values */
      for ( ifdot=0; ifdot<nfdot; ifdot++)
	{
	  /* set spindown value for Fstat calculation */
	  FstatPar1.fdot->data[0] = uvar_fdot + ifdot * dfDot;
	  
	  /* set up memory for Fstat vectors */
	  LAL_CALL(SetUpFstatStack( &status, &FstatVect1, &FstatPar1), &status);
	  
	  /* calculate the Fstatistic */
	  LAL_CALL(ComputeFstatStack( &status, &FstatVect1, &stackSFTs1, &FstatPar1), &status);
	  
	  /* print fstat vector if required -- mostly for debugging */
	  if ( uvar_printFstat1 )
	    {
	      for (k=0; k<nStacks1; k++)
		LAL_CALL( PrintFstatVec_fp ( &status, FstatVect1.data + k, fpFstat1, FstatPar1.alpha, 
					     FstatPar1.delta, FstatPar1.fdot->data[0]), &status); 
	    }


	  
	  /*------------ select peaks ------------*/ 
	      
	  LAL_CALL( FstatVectToPeakGram( &status, &pgV, &FstatVect1, uvar_FstatPeakSelect), &status);
	  
	  /* free Fstat vectors -- we don't need them now because we have the peakgrams */
	  for(k=0; k<nStacks1; k++) {
	    LALFree(FstatVect1.data[k].data->data);
	    LALFree(FstatVect1.data[k].data);
	  }
	  LALFree(FstatVect1.data);
	  
	  
	  /*--------------- calculate Hough map and get candidates ---------------*/
	  
	  /* set sky location and spindown for Hough grid -- same as for Fstat calculation */	  
	  houghPar.alpha = FstatPar1.alpha;
	  houghPar.delta = FstatPar1.delta;
	  houghPar.fdot->data[0] = FstatPar1.fdot->data[0];
	  
	  /* get candidates */
	  LAL_CALL ( ComputeFstatHoughMap ( &status, &houghCand1, &pgV, &houghPar, LALUserVarWasSet(&uvar_houghThr) ), &status);
	  
	  /* free peakgrams -- we don't need them now because we have the Hough maps */
	  for (k=0; k<nStacks1; k++) 
	    LALFree(pgV.pg[k].peak);
	  LALFree(pgV.pg);
	  
	  /* print candidates */
	  if ( uvar_printCand1 )
	    LAL_CALL ( PrintHoughCandidates ( &status, &houghCand1, fpHough), &status);
	  
	  
	  /*------------- Follow up candidates --------------*/
	  
	  /* this part is more general than it has to be
	     it is meant to be generalized to the case when 
	     the number of follow-up stacks is not necessarily 1 */
	  
	  /* check if user requested a follow up stage*/
	  if ( uvar_followUp ) 
	    {
	      
	      /* loop over candidates surviving 1st stage  */
	      for ( j=0; j < houghCand1.nCandidates; j++) 
		{		  
		  /* 2nd stage frequency and spindown variables */
		  /* these are the parameters passed down from the 
		     first stage which is why they have the subscript 1 */
		  REAL8 fStart1, freqBand1, fdot1, fdotBand1;
		  REAL8 alpha1, delta1, alphaBand1, deltaBand1;
		  
		  /* set frequency and spindown ranges */
		  fStart1 = houghCand1.list[j].freq - 0.5 * houghCand1.list[j].dFreq;
		  freqBand1 = houghCand1.list[j].dFreq;
		  fdot1 = houghCand1.list[j].fdot - 0.5 * houghCand1.list[j].dFdot;
		  fdotBand1 = houghCand1.list[j].dFdot;
		  
		  /* set frequency range */
		  FstatPar2.fStart = fStart1;
		  FstatPar2.fBand = freqBand1;
		  
		  /* set sky region to be refined */
		  alpha1 = houghCand1.list[j].alpha - 0.5 * houghCand1.list[j].dAlpha;
		  delta1 = houghCand1.list[j].delta - 0.5 * houghCand1.list[j].dDelta;
		  alphaBand1 = houghCand1.list[j].dAlpha;
		  deltaBand1 = houghCand1.list[j].dDelta;
		  LAL_CALL (SkySquare2String( &status, &(scanInit2.searchRegion.skyRegionString),
					      alpha1, delta1, alphaBand1, deltaBand1), &status);
		  
		  /* set second doppler scan variables */
		  scanInit2.searchRegion.Freq = fStart1;
		  scanInit2.searchRegion.FreqBand = freqBand1;
		  scanInit2.searchRegion.f1dot = fdot1;
		  scanInit2.searchRegion.f1dotBand = fdotBand1;
		  
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
		      
		      /*------------- calculate F statistic for each stack --------------*/
		      
		      /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
		      thisPoint2.longitude = dopplerpos2.Alpha;
		      thisPoint2.latitude = dopplerpos2.Delta;
		      thisPoint2.system = COORDINATESYSTEM_EQUATORIAL;
		      LAL_CALL (LALNormalizeSkyPosition(&status, &thisPoint2, &thisPoint2), &status);
		      
		      /* set sky template points */
		      FstatPar2.alpha = thisPoint2.longitude;
		      FstatPar2.delta = thisPoint2.latitude;
		      
		      /* number of fdot values */
		      dfDot2 = thisScan2.df1dot;
		      nfdot2 = (UINT4)( fdotBand1 / dfDot2 + 0.5) + 1; 
		      
		      /* loop over fdot values */
		      for ( ifdot2=0; ifdot2<nfdot2; ifdot2++)
			{
			  
			  /* set spindown value for Fstat calculation */
			  FstatPar2.fdot->data[0] = fdot1 + ifdot2 * dfDot2;
			  
			  /* set up memory for Fstat vectors */
			  LAL_CALL(SetUpFstatStack( &status, &FstatVect2, &FstatPar2), &status);
			  
			  /* calculate the Fstatistic */
			  LAL_CALL(ComputeFstatStack( &status, &FstatVect2, &stackSFTs2, &FstatPar2), &status);
			  
			  /* select candidates from 2nd stage */
			  for (k=0; k<nStacks2; k++) 
			    {
			      if ( LALUserVarWasSet(&uvar_fstatThr))
				LAL_CALL( GetFstatCandidates( &status, fstatToplist, FstatVect2.data + k, uvar_fstatThr,
							      FstatPar2.alpha, FstatPar2.delta, 
							      FstatPar2.fdot->data[0], uvar_reallocBlock ), &status);
			      else
				LAL_CALL( GetFstatCandidates_toplist( &status, fstatToplist, FstatVect2.data + k, 
								      FstatPar2.alpha, FstatPar2.delta, 
								      FstatPar2.fdot->data[0] ), &status);
			      
			      /* free Fstat vector */
			      LALFree(FstatVect2.data[k].data->data);
			      LALFree(FstatVect2.data[k].data);
			    } /* end loop over nstacks2 for selecting candidates */
			  
			  LALFree(FstatVect2.data);
			  
			} /* end loop over refined spindown parameters */
		      
		    } /* end while loop over second stage refined sky grid */
		  
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
      if ( write_toplist_to_fp( fstatToplist, fpFstat, &checksum) < 0)
	 fprintf( stderr, "Error in writing toplist to file\n");
    /*    LAL_CALL( AppendFstatCandidates( &status, &fStatCand, fpFstat), &status); */
  }
	 
  /*------------ free all remaining memory -----------*/
  
  /* free memory */

  if ( uvar_followUp )
    {
      LALFree(tStack2);
      LALFree(mCohSft2);
      LAL_CALL (LALDDestroyVector (&status, &(FstatPar2.fdot)), &status); 
      LALFree(midTstack2.data);
      LALFree(stackSFTs2.data);
     
      if ( fstatToplist ) 
	free_toplist(&fstatToplist);
      
      fclose(fpFstat);
      LALFree(fnameFstatCand);
    }

  if ( uvar_printCand1 )
    {
      LALFree(fnameHoughCand);
      fclose(fpHough);
    }

  if ( uvar_printFstat1 )
    fclose(fpFstat1);

  /* free first stage memory */
  LAL_CALL(LALDestroyTimestampVector ( &status, &startTsft1), &status);  	  
  LALFree(stackSFTs1.data);
  LAL_CALL (LALDestroySFTVector(&status, &inputSFTVec1),&status );
  LALFree(mCohSft1);
  LALFree(midTstack1.data);
  LALFree(tStack1);
  LAL_CALL (LALDDestroyVector (&status, &(houghPar.fdot)), &status);
  LAL_CALL (LALDDestroyVector (&status, &(FstatPar1.fdot)), &status); 

  /* free Vel/Pos vectors and ephemeris */
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);
  LAL_CALL( LALDDestroyVectorSequence (&status,  &velStack), &status);
  LAL_CALL( LALDDestroyVectorSequence (&status,  &posStack), &status);
  
  /* free dopplerscan stuff */
  LAL_CALL ( FreeDopplerScan(&status, &thisScan1), &status);
  if ( scanInit1.searchRegion.skyRegionString )
    LALFree ( scanInit1.searchRegion.skyRegionString );

  /* free second stage memory if required */
  if ( ( LALUserVarWasSet(&uvar_sftDir2)) )
    {     
      LAL_CALL (LALDestroySFTVector(&status, &inputSFTVec2),&status );
      LAL_CALL(LALDestroyTimestampVector ( &status, &startTsft2), &status);  	
    }
 
  /* free candidates */
  LALFree(houghCand1.list);

  LAL_CALL (LALDestroyUserVars(&status), &status);  

  LALCheckMemoryLeaks();

  return DRIVEHOUGHFSTAT_ENORM;
}




/** \brief Function for calculating the Fstatistic for a stack of SFT data
    \param *stackSFTs : pointer to SFTVectorSequence -- a sequence of SFT vectors
    \param *params : structure of type FstatStackParams -- parameters for calculating Fstatistic 
    \return *out : pointer to REAL8FrequencySeriesVector -- sequence of Fstatistic vectors

    This function takes a set of SFTs broken up into stacks appropriately 
    and calculates the Fstatistic for each stack over a frequency range for 
    a single sky-location and a single value of the first spindown (the demodulation parameters). 
    It allocates memory for the Fstatistic vectors appropriately which must 
    be freed later outside the function.  It uses the ComputeFstatistic_v2 code
    as a prototype.  The output set of Fstatistic vectors might be 
    combined semicoherently in a small parameter space patch around the 
    demodulation parameters.  Thus, the Fstatistic is calculated for a somewhat 
    larger frequency range than what the user specified if the number of stacks exceeds unity.
*/
void ComputeFstatStack (LALStatus *status, 
			REAL8FrequencySeriesVector *out, 
			SFTVectorSequence *stackSFTs, 
			FstatStackParams *params)
{
  /* stuff copied from params */
  INT4 nStacks = params->nStacks;
  REAL8 refTime = params->refTime;
  INT4 *mCohSft = params->mCohSft;
  REAL8Vector *fdot = params->fdot;
  REAL8 fStart = params->fStart;
  REAL8 tStackAvg = params->tStackAvg;
  REAL8 deltaF = 1.0/tStackAvg;
  REAL8 fBand = params->fBand;

  /* copy timeBase from SFT vector */
  REAL8 timeBase = 1.0 / ( stackSFTs->data->data->deltaF);

  /* other variables */
  SSBtimes tSSB;
  AMCoeffs amcoe;
  SkyPosition skyPoint;
  DetectorStateSeries *DetectorStates=NULL;
  LIGOTimeGPSVector timeStack;
  INT4 k, j, indexSft;
  REAL8Vector *fkdot=NULL;
  LIGOTimeGPS tempRef;
  INT8 binsFstat;

  INITSTATUS( status, "ComputeFstatStack", rcsid );
  ATTATCHSTATUSPTR (status);
   
  /* other stuff copied from params */
  skyPoint.longitude = params->alpha;
  skyPoint.latitude = params->delta;
  skyPoint.system = COORDINATESYSTEM_EQUATORIAL;
  TRY (LALNormalizeSkyPosition( status->statusPtr, &skyPoint, &skyPoint), status);

  /* set reference time in GPS struct */
  TRY (LALFloatToGPS ( status->statusPtr, &tempRef, &refTime), status);

  /* copy spindown */
  fkdot = NULL;
  TRY ( LALDCreateVector( status->statusPtr, &(fkdot), 2), status);
  fkdot->data[1] = fdot->data[0];

  /*---------------- start loop over stacks ------------------*/
  indexSft = 0;
  for(k=0; k<nStacks; k++) {

    /* set timestamps vector for sfts in stack */
    timeStack.length = mCohSft[k];
    timeStack.data = params->ts->data + indexSft;

    /* obtain detector positions and velocities, together with LMSTs for the SFT midpoints (i.e. shifted by tSFT/2) */
    TRY ( LALGetDetectorStates ( status->statusPtr, &DetectorStates, &timeStack, 
				      &(params->detector), params->edat, 
				      timeBase / 2.0 ), status);


    /* allocate memory for am coeffs */
    amcoe.a = NULL;
    amcoe.b = NULL;
    TRY (LALSCreateVector(status->statusPtr, &(amcoe.a), mCohSft[k]), status);
    TRY (LALSCreateVector(status->statusPtr, &(amcoe.b), mCohSft[k]), status);

    /* allocate memory for tssb */
    tSSB.DeltaT = NULL;
    tSSB.Tdot = NULL;
    TRY ( LALDCreateVector(status->statusPtr, &(tSSB.DeltaT), mCohSft[k]), status );
    TRY ( LALDCreateVector(status->statusPtr, &(tSSB.Tdot), mCohSft[k]), status );

    /* loop over frequency bins and get Fstatistic */
    binsFstat = floor( tStackAvg * fBand );
    for(j=0; j<binsFstat; j++) {

      /* increase frequency value */
      fkdot->data[0] = fStart + j*deltaF;

      /* transform to SSB frame */ 
      TRY ( LALGetSSBtimes ( status->statusPtr, &tSSB, DetectorStates, skyPoint, 
			     tempRef , params->SSBprecision), status);
     

      /* calculate amplitude modulation coefficients */
      TRY ( LALGetAMCoeffs( status->statusPtr, &amcoe, DetectorStates, skyPoint), status);      

      {
	/* get the F statistic */
	Fcomponents FaFb;
	REAL4 fact;
	REAL4 At, Bt, Ct;
	REAL4 FaRe, FaIm, FbRe, FbIm;
	
	XLALComputeFaFb ( &FaFb, stackSFTs->data + k, fkdot, &tSSB, 
			  &amcoe, params->Dterms);
	At = amcoe.A;
	Bt = amcoe.B;
	Ct = amcoe.C;
	
	FaRe = FaFb.Fa.re;
	FaIm = FaFb.Fa.im;	    
	FbRe = FaFb.Fb.re;
	FbIm = FaFb.Fb.im;
	
	fact = 2.0f / (1.0f * stackSFTs->data[k].length * amcoe.D);    

	/* fill up output vector */	
	out->data[k].data->data[j] = fact * (Bt * (FaRe*FaRe + FaIm*FaIm) 
					     + At * (FbRe*FbRe + FbIm*FbIm) 
					     - 2.0f * Ct *(FaRe*FbRe + FaIm*FbIm) );
      } /* end fstat calculation block */
    
    } /* end loop over frequencies */

    /* increment over correct number of sfts */
    indexSft += mCohSft[k];

    /*---------- clear memory -----------*/
    /* destroy DetectorStateSeries */
    TRY ( LALDestroyDetectorStateSeries (status->statusPtr, &DetectorStates), status);

    /* Free AM-coefficients */
    TRY (LALSDestroyVector(status->statusPtr, &(amcoe.a)), status);
    TRY (LALSDestroyVector(status->statusPtr, &(amcoe.b)), status);
    /* Free SSB-times */
    TRY (LALDDestroyVector(status->statusPtr, &(tSSB.DeltaT)), status);
    TRY (LALDDestroyVector(status->statusPtr, &(tSSB.Tdot)), status);

  } /* end loop over stacks */

  TRY (LALDDestroyVector ( status->statusPtr, &(fkdot)), status);
  
  DETATCHSTATUSPTR (status);
  RETURN(status); 
}



/** Function for allocating memory for Fstat vectors */
void SetUpFstatStack (LALStatus *status, 
		      REAL8FrequencySeriesVector *out, 
		      FstatStackParams *params)
{
  /* stuff copied from params */
  INT4 nStacks = params->nStacks;
  INT4 nfSizeCylinder = params->nfSizeCylinder;
  REAL8 fStart = params->fStart;
  REAL8 fBand = params->fBand;
  REAL8 tStackAvg = params->tStackAvg;
  REAL8 deltaF = 1.0/tStackAvg;


  /* other variables */
  INT4 k;
  INT8 binsFstat;

  INITSTATUS( status, "ComputeFstatStack", rcsid );
  ATTATCHSTATUSPTR (status);


  /*---------- memory for Fstatistic Vector -----------*/

  /* number of bins for calculating Fstat */
  /* add extraBins on either side*/
  {
    /* extraBins = nfSizeCylinder/2 + maxNBins/2 
       nfSizeCylinder is parameter, but maxNBins must 
       be calculated from Hough routines.  It is the 
       largest number of bins affected by the skypatch */

    INT4 extraBins; /* the extra number of bins for which the Fstat must be calculated */
    HOUGHResolutionPar resPar;
    HOUGHSizePar sizePar;
    INT8 tempBin, fBinIni, fBinFin;
    REAL8 patchSize;

    /* extraBins required only if there is a hough follow up 
       which happens only if nStacks is more than 1 */
    if (nStacks > 1) {
      
      /* calculate sizePar for fStart */
      tempBin = (INT8) ( tStackAvg * fStart );    
      patchSize = 0.5 / ( tempBin * VEPI ); 

      resPar.f0Bin = tempBin;
      resPar.deltaF = deltaF;
      resPar.patchSkySizeX = patchSize;
      resPar.patchSkySizeY = patchSize;
      resPar.pixelFactor = PIXELFACTOR;
      resPar.pixErr = PIXERR;
      resPar.linErr = LINERR;
      resPar.vTotC = VTOT;
      
      TRY ( LALHOUGHComputeSizePar ( status->statusPtr, &sizePar, &resPar ), status);
      extraBins = nfSizeCylinder/2 + sizePar.maxNBins/2;
    }
    else
      extraBins = 0;

    /* now we can calculate required span of Fstat vector */
    binsFstat = floor( tStackAvg * fBand ) + 2*extraBins;
    fBinIni = floor( tStackAvg * fStart + 0.5) - extraBins;
    fBinFin = fBinIni + binsFstat - 1;

    out->length = nStacks;
    out->data = NULL;
    out->data = (REAL8FrequencySeries *)LALMalloc(nStacks * sizeof(REAL8FrequencySeries));
    for (k=0; k<nStacks; k++) {
      out->data[k].epoch = params->tsStack->data[k];
      out->data[k].deltaF = deltaF;
      out->data[k].f0 = deltaF * fBinIni;
      out->data[k].data = (REAL8Sequence *)LALMalloc(sizeof(REAL8Sequence));
      out->data[k].data->length = binsFstat;
      out->data[k].data->data = (REAL8 *)LALMalloc( binsFstat * sizeof(REAL8));
    }
  } /* end of Fstat memory allocation */
   
  DETATCHSTATUSPTR (status);
  RETURN(status); 
}




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
			  HoughCandidates  *out,   /* output candidates */
			  HOUGHPeakGramVector *pgV, /* peakgram vector */
			  HoughParams *params,
			  BOOLEAN selectCandThr)
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
  INT4  k, iHmap, nSpin1Max, nStacks, nfSizeCylinder;
  REAL8 deltaF, alpha, delta;
  REAL8 patchSizeX, patchSizeY, f1jump, tStart;
  REAL8VectorSequence *vel, *pos;
  REAL8Vector *fdot;
  LIGOTimeGPSVector   *ts;
  REAL8Vector timeDiffV;
  REAL8 houghThr;
  UINT4Vector hist; /* histogram vector */ 
  UINT4Vector histTotal; /* total histogram vector */
  HoughStats stats; /* statistics struct */
  CHAR *fileStats = NULL;
  FILE *fpStats = NULL;

  INITSTATUS( status, "ComputeFstatHoughMap", rcsid );
  ATTATCHSTATUSPTR (status);

  /* copy some params to local variables */
  houghThr = params->houghThr;
  nfSizeCylinder = params->nfSizeCylinder;
  fBinIni = params->fBinIni;
  fBinFin = params->fBinFin;
  alpha = params->alpha;
  delta = params->delta;
  vel = params->vel;
  pos = params->pos;
  fdot = params->fdot;
  ts = params->ts;
  tStart = params->tStart;

  /* copy some parameters from peakgram vector */
  deltaF = pgV->pg->deltaF;
  nStacks = pgV->length;

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
  phmdVS.nfSize  = nfSizeCylinder;
  phmdVS.deltaF  = deltaF;
  phmdVS.phmd = NULL;
  phmdVS.phmd=(HOUGHphmd *)LALMalloc(nStacks*nfSizeCylinder*sizeof(HOUGHphmd));
  
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

  /* Case: no spindown */
  parDem.deltaF = deltaF;
  parDem.skyPatch.alpha = alpha;
  parDem.skyPatch.delta = delta;
  parDem.spin.length = fdot->length;
  parDem.spin.data = fdot->data;
  
  parRes.deltaF = deltaF;
  parRes.patchSkySizeX  = patchSizeX;
  parRes.patchSkySizeY  = patchSizeY;
  parRes.pixelFactor = PIXELFACTOR;
  parRes.pixErr = PIXERR;
  parRes.linErr = LINERR;
  parRes.vTotC = VTOT;


  /* memory allocation for histogram and opening stats file*/
  if ( uvar_printStats ) {
    hist.length = nStacks+1;
    histTotal.length = nStacks+1;
    hist.data = NULL;
    histTotal.data = NULL;
    hist.data = (UINT4 *)LALMalloc((nStacks+1)*sizeof(UINT4));
    histTotal.data = (UINT4 *)LALMalloc((nStacks+1)*sizeof(UINT4));
    { 
      UINT4   j;
      for(j=0; j< histTotal.length; ++j) 
	histTotal.data[j]=0; 
    }
    fileStats = NULL;
    fileStats = (CHAR *)LALMalloc( 256 * sizeof(CHAR));
    strcpy( fileStats, params->outBaseName);
    strcat( fileStats, "stats");
    if ( !(fpStats = fopen(fileStats, "w")))
      {
	fprintf(stderr, "Unable to open file '%s' for writing\n", fileStats);
	exit(1);
      }
  }



  /* calculate time differences from start of observation time */
  timeDiffV.length = nStacks;
  timeDiffV.data = (REAL8 *)LALMalloc( nStacks * sizeof(REAL8));
  for (k=0; k<nStacks; k++) {
    REAL8 tMidStack;

    TRY ( LALGPStoFloat ( status->statusPtr, &tMidStack, ts->data + k), status);
    timeDiffV.data[k] = tMidStack - tStart;
  }

  /* if there are residual spindowns */
  nSpin1Max = nfSizeCylinder/2; /* integer division -- maximum number of spindowns */
  f1jump = 1.0 / timeDiffV.data[nStacks - 1]; /* resolution in residual fdot */
 
  /*------------------ start main Hough calculation ---------------------*/

  /* initialization */  
  fBin= fBinIni; /* initial search bin */
  iHmap = 0; /* hough map index */

  while( fBin <= fBinFin){
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
      parDem.timeDiff = timeDiffV.data[j];

      /* calculate parameters needed for buiding the LUT */
      TRY( LALHOUGHParamPLUT( status->statusPtr, &parLut, &parSize, &parDem),status );
      /* build the LUT */
      TRY( LALHOUGHConstructPLUT( status->statusPtr, &(lutV.lut[j]), &patch, &parLut ),
	   status );
    }
    
    /*--------- build the set of  PHMD centered around fBin -------------*/     
    phmdVS.fBinMin = fBin - floor(nfSizeCylinder/2.);
    TRY( LALHOUGHConstructSpacePHMD(status->statusPtr, &phmdVS, pgV, &lutV), status );
    
    /*-------------- initializing the Total Hough map space ------------*/   
    ht.xSide = xSide;
    ht.ySide = ySide;
    ht.skyPatch.alpha = alpha;
    ht.skyPatch.delta = delta;
    ht.mObsCoh = nStacks;
    ht.deltaF = deltaF;
    ht.spinDem.length = fdot->length;
    ht.spinDem.data = fdot->data;
    ht.patchSizeX = patchSizeX;
    ht.patchSizeY = patchSizeY;
    ht.dFdot.data[0] = deltaF * f1jump;
    ht.map   = NULL;
    ht.map   = (HoughTT *)LALMalloc(xSide*ySide*sizeof(HoughTT));
    TRY( LALHOUGHInitializeHT( status->statusPtr, &ht, &patch), status); /*not needed */
    
    /*  Search frequency interval possible using the same LUTs */
    fBinSearch = fBin;
    fBinSearchMax = fBin + parSize.nFreqValid - 1 - (nfSizeCylinder - 1 )/2;
     
    /* Study all possible frequencies with one set of LUT */    
    while ( (fBinSearch <= fBinFin) && (fBinSearch < fBinSearchMax) )  { 

      /* finally we can construct the hough maps and select candidates */
      {
	INT4   n;
	REAL8  f1dis;

	ht.f0Bin = fBinSearch;
	    
	/*loop over all values of residual spindown */
	for( n=-nSpin1Max; n<= nSpin1Max; n++ ){ 
	  f1dis =  n*f1jump;

	  ht.spinRes.data[0] =  f1dis*deltaF;
	  
	  for (j=0; j < (UINT4)nStacks; j++)
	    freqInd.data[j] = fBinSearch + floor(timeDiffV.data[j]*f1dis + 0.5);
	  
	  TRY( LALHOUGHConstructHMT(status->statusPtr, &ht, &freqInd, &phmdVS),status );

	  /* get candidates */
	  if ( selectCandThr )
	    TRY(GetHoughCandidates( status->statusPtr, out, &ht, &patch, &parDem, houghThr), status);
	  else
	    TRY(GetHoughCandidates_toplist( status->statusPtr, out, &ht, &patch, &parDem), status);

	  /* calculate statistics and histogram */
	  if ( uvar_printStats ) {
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
  LALFree(ht.dFdot.data);
  LALFree(lutV.lut);
  LALFree(phmdVS.phmd);
  LALFree(freqInd.data);
  LALFree(timeDiffV.data);

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


/** Given a sftVector, this function splits it up into several different sft vectors 
    there are basically two ways of doing this: either each stack contains the same number
    of SFTs, or each stack spans the same duration. These two methods are equivalent
    only if there are no gaps between the SFTs 
    
    The function SetUpStacks1 distributes the SFT equally among the stacks while 
    SetUpStacks2 makes each stack span the same time duration. 
*/
void SetUpStacks1(LALStatus *status, 
		 SFTVectorSequence  *out,  
		 SFTVector  *sftVect,
		 INT4 nStacks)
{
  INT4 k, mCohSft, nSFTs;

  INITSTATUS( status, "SetUpStacks1", rcsid );
  ATTATCHSTATUSPTR (status);

  out->length = nStacks;
  out->data = (SFTVector *)LALMalloc( nStacks * sizeof(SFTVector));

  nSFTs = sftVect->length;
  mCohSft = nSFTs/nStacks; /* integer division -- some sfts will be discarded */ 

  for (k=0; k<nStacks; k++) {
    SFTVector *tempVect;
    tempVect = out->data + k;
    tempVect->length = mCohSft;
    /* point the output to the right elements of sftVect */
    tempVect->data = sftVect->data + mCohSft*k;
  }  

  DETATCHSTATUSPTR (status);
  RETURN(status);
}


void SetUpStacks2(LALStatus *status, 
		  SFTVectorSequence  *out,  
		  SFTVector  *sftVect,
		  LIGOTimeGPSVector *ts,
		  INT4 nStacks)
{
  REAL8 tStart, tEnd, tStack, timeBase;
  INT4 k, j;
  REAL8 thisTime;
  INT4 nSFTs, sftCount;
  SFTVector *tempVect;

  INITSTATUS( status, "SetUpStacks2", rcsid );
  ATTATCHSTATUSPTR (status);

  out->length = nStacks;
  out->data = (SFTVector *)LALMalloc( nStacks * sizeof(SFTVector));

  nSFTs = ts->length;
  timeBase = 1.0/sftVect->data->deltaF;
  TRY ( LALGPStoFloat ( status->statusPtr, &tStart, ts->data), status);
  TRY ( LALGPStoFloat ( status->statusPtr, &tEnd, ts->data + nSFTs), status);
  tEnd += timeBase;
  tStack = (tEnd - tStart) / nStacks;

  /* count number of sfts between tStart + k*tStack and tStart + (k+1)*tStack */
  k = 0; /* initialization -- stack label */
  sftCount = 0; /* initialization -- number of sfts in a stack */

  /* loop over the sfts and find out if it belongs to the k^th stack */
  for( j=0; j<nSFTs; j++) {

    /* if sft time stamp is within the k^th stack, then go on to next sft 
       otherwise set up the k^th output vector */
    TRY ( LALGPStoFloat ( status->statusPtr, &thisTime, ts->data + j), status);
    if ( thisTime < tStart + (k+1)*tStack ) {
      sftCount++;
    }
    else {
      /* set up the output vector only if there are sfts */
      if (sftCount) {
	tempVect = out->data + k;
	tempVect->length = sftCount;
	tempVect->data = sftVect->data + j;

	/* increment stack counter and reset sft counter */
	k++;
	sftCount = 0;
      }
    }
  }
  
  /* if some stacks were empty then realloc output vector */
  if ( k<nStacks ) {
    out->length = k;
    out->data = (SFTVector *)LALRealloc( out->data, k*sizeof(SFTVector));
  }

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
  if ( !(fp=fopen(filename,"w"))) 
    { 
      fprintf(stderr, "Unable to open file '%s' for writing\n", filename);
      exit(1);
    }

  /* replace this by an assert */
  /*   if ( !fp ){   */
  /*     fprintf(stderr,"Unable to find file %s\n",filename); */
  /*     return DRIVEHOUGHFSTAT_EFILE;  */
  /*   } */

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
			HoughCandidates *houghCand,
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

	nCandidates = houghCand->nCandidates;

	/* if there isn't enough memory then realloc */
	if ( nCandidates >= houghCand->length ) {
	  houghCand->length += 5000;
	  houghCand->list = (HoughList *)LALRealloc( houghCand->list,
						     houghCand->length * sizeof(HoughList));
	} /* end of reallocs */

	houghCand->list[nCandidates].freq = f0;
	houghCand->list[nCandidates].dFreq = deltaF;
	
	houghCand->list[nCandidates].fdot = fdot;
	houghCand->list[nCandidates].dFdot = dFdot;

	/* get sky location of pixel */
	TRY( LALStereo2SkyLocation (status->statusPtr, &sourceLocation, 
				    j, i, patch, parDem), status);

	houghCand->list[nCandidates].alpha = sourceLocation.alpha;
	houghCand->list[nCandidates].delta = sourceLocation.delta;

	houghCand->list[nCandidates].dAlpha = patchSizeX / ((REAL8)xSide);
	houghCand->list[nCandidates].dDelta = patchSizeY / ((REAL8)ySide);

	/* increment candidate count */
	houghCand->nCandidates += 1;

      } /* end if statement */
    } /* end loop over xSide */
  } /* end loop over ySide */
  DETATCHSTATUSPTR (status);
  RETURN(status);
}


/** Get Hough candidates */
void GetHoughCandidates_toplist(LALStatus *status,
				HoughCandidates *houghCand,
				HOUGHMapTotal *ht,
				HOUGHPatchGrid  *patch,
				HOUGHDemodPar   *parDem)
{
  REAL8UnitPolarCoor sourceLocation;
  REAL8 deltaF, f0, fdot, dFdot, patchSizeX, patchSizeY;
  REAL8 mean, std;
  INT8 f0Bin;  
  INT4 nCandidates, maxCandidates; 
  INT4 i,j, xSide, ySide, minSigIndex;
  HoughStats stats;

  INITSTATUS( status, "GetHoughCandidates_toplist", rcsid );
  ATTATCHSTATUSPTR (status);

  deltaF = ht->deltaF;
  f0Bin = ht->f0Bin;
  f0 = f0Bin * deltaF;

  nCandidates = houghCand->nCandidates;
  maxCandidates = houghCand->length;

  fdot = ht->spinDem.data[0] + ht->spinRes.data[0];
  dFdot = ht->dFdot.data[0];

  xSide = ht->xSide;
  ySide = ht->ySide;

  patchSizeX = ht->patchSizeX;
  patchSizeY = ht->patchSizeY;

  TRY( LALHoughStatistics ( status->statusPtr, &stats, ht), status );
  mean = stats.avgCount;
  std = stats.stdDev;

  for (i = 0; i < ySide; i++)
    {
      for (j = 0; j < xSide; j++)
	{ 
	  REAL8 tempSig = (ht->map[i*xSide + j] - mean)/std;

	  if ( nCandidates < maxCandidates )
	    {
	      houghCand->list[nCandidates].freq = f0;
	      houghCand->list[nCandidates].dFreq = deltaF;
	      
	      houghCand->list[nCandidates].fdot = fdot;
	      houghCand->list[nCandidates].dFdot = dFdot;
	      
	      /* get sky location of pixel */
	      TRY( LALStereo2SkyLocation (status->statusPtr, &sourceLocation, 
					  j, i, patch, parDem), status);

	      houghCand->list[nCandidates].alpha = sourceLocation.alpha;
	      houghCand->list[nCandidates].delta = sourceLocation.delta;
	      houghCand->list[nCandidates].dAlpha = patchSizeX / ((REAL8)xSide);
	      houghCand->list[nCandidates].dDelta = patchSizeY / ((REAL8)ySide);
	      houghCand->list[nCandidates].significance = tempSig;

	      TRY ( GetMinSigIndex_toplist ( status->statusPtr, &minSigIndex, houghCand), status);		
	      houghCand->minSigIndex = minSigIndex;

	      /* increment candidate count */
	      houghCand->nCandidates += 1;
	      nCandidates = houghCand->nCandidates;

	    }
	  else 
	    {
	      /* if event is more significant than least significant 
		 event stored in toplist then replace it */
	      minSigIndex = houghCand->minSigIndex;
	      if ( tempSig > houghCand->list[minSigIndex].significance ) 
		{
		  houghCand->list[minSigIndex].freq = f0;
		  houghCand->list[minSigIndex].dFreq = deltaF;
		  houghCand->list[minSigIndex].fdot = fdot;
		  houghCand->list[minSigIndex].dFdot = dFdot;
	  
		  /* get sky location of pixel */
		  TRY( LALStereo2SkyLocation (status->statusPtr, &sourceLocation, 
					      j, i, patch, parDem), status);
		  
		  houghCand->list[minSigIndex].alpha = sourceLocation.alpha;
		  houghCand->list[minSigIndex].delta = sourceLocation.delta;  
		  houghCand->list[minSigIndex].dAlpha = patchSizeX / ((REAL8)xSide);
		  houghCand->list[minSigIndex].dDelta = patchSizeY / ((REAL8)ySide);
		  houghCand->list[minSigIndex].significance = tempSig;

		  /* find index of new least significant event */
		  TRY ( GetMinSigIndex_toplist ( status->statusPtr, &minSigIndex, houghCand), status);
		  houghCand->minSigIndex = minSigIndex;

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
			    HoughCandidates *houghCand)
{

  REAL8 minSig;
  INT4 k, nCandidates = houghCand->nCandidates;

  INITSTATUS( status, "GetHoughCandidates_toplist", rcsid );
  ATTATCHSTATUSPTR (status);


  *minSigIndex = 0;
  minSig = houghCand->list[0].significance;
  for (k = 1; k < nCandidates; k++)
    {
      REAL8 tempSig = houghCand->list[k].significance;
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
void PrintHoughCandidates(LALStatus *status,
			  HoughCandidates *in,
			  FILE *fp)
{
  INT4 k;

  INITSTATUS( status, "GetHoughCandidates", rcsid );
  ATTATCHSTATUSPTR (status);

 
  for (k=0; k < in->nCandidates; k++)
    fprintf(fp, "%e   %e   %g   %g   %g   %g   %e   %e\n", in->list[k].freq, 
	    in->list[k].dFreq, in->list[k].alpha, in->list[k].dAlpha, in->list[k].delta, in->list[k].dDelta,
	    in->list[k].fdot, in->list[k].dFdot);


  DETATCHSTATUSPTR (status);
  RETURN(status);
}


/** Print Fstat vectors -- mostly for debugging purposes */
  void PrintFstatVec_fp (LALStatus *status,
			 REAL8FrequencySeries *in,
			 FILE *fp,
			 REAL8 alpha,
			 REAL8 delta,
			 REAL8 fdot)
{
  INT4 length, k;
  REAL8 f0, deltaF, freq;

  INITSTATUS( status, "PrintFstatVec_fp", rcsid );
  ATTATCHSTATUSPTR (status);


  length = in->data->length;
  deltaF = in->deltaF;
  f0 = in->f0;

  for (k=0; k<length; k++)
    {
      freq = f0 + k*deltaF;
      fprintf(fp, "%.13g %.7g %.7g %.5g %.6g\n", freq, alpha, delta, fdot, in->data->data[k]);
    }

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
    create_toplist( &list, blockRealloc); 

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
      
      debug = insert_into_toplist( list, line);

    }

  DETATCHSTATUSPTR (status);
  RETURN(status);
}


/** Print hough histogram to a file */
void PrintHoughHistogram( LALStatus *status,
		     UINT4Vector *hist, 
		     CHAR *fnameOut)
{

  FILE  *fp=NULL;   /* Output file */
  char filename[256];
  UINT4  i ;
 
  INITSTATUS( status, "GetHoughCandidates", rcsid );
  ATTATCHSTATUSPTR (status);

  strcpy(  filename, fnameOut);
  strcat(  filename, "histo");
  if ( !(fp=fopen(filename,"w")))
    {
      fprintf(stderr, "Unable to open file '%s' for writing\n", filename);
      exit(1);
    }

  for (i=0; i < hist->length; i++)
    fprintf(fp,"%d  %d\n", i, hist->data[i]);
  
  fclose( fp );  
  
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
