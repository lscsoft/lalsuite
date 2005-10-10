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

   - The user specifies a region in parameter space to search over.  At present, 
     only a single sky-location and spindown are allowed, though a frequency range 
     can be specified.  The F-statistic is calculated for each stack at the chosen
     sky-position, spindown, and frequency range.  

   - A threshold is set on the F-statistic to convert the 
     F-statistic vector into a vector of 0s and 1s known as a \e peakgram -- there 
     is one peakgram for each stack.

   - The peakgrams are combined using the Hough transform.

   - The Hough part of the search constructs a grid in a small patch around the 
     chosen sky-position and spindown and stacks up the different segments in frequency
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

   - Do Fstat memory allocation outside Fstat calculation function

   \par Longer term

   - Should we over-resolve the Fstat calculation to reduce loss in signal power?  We would
     still calculate the peakgrams at the 1/T resolution, but the peak selection would 
     take into account Fstat values over several over-resolved bins.  
   
   - What is the best grid for calculating the F-statistic grid?  At first glance, the 
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
#define ALPHA 1.57    /**< Default sky location -- right ascension */
#define DELTA  0.0    /**< Default sky location -- declination */
#define NFSIZE  21    /**< Default size of hough cylinder of look up tables */
#define DTERMS 8     /**< Default number of dirichlet kernel terms for calculating Fstat */
#define FSTATTHRESHOLD 2.6  /**< Default threshold on Fstatistic for peak selection */
#define SFTDIRECTORY "/home/badkri/fakesfts/"  /**< Default directory containing sfts */
#define FNAMEOUT "./temp/OutHoughFStat"  /**< Default output file basename */

int main( int argc, char *argv[]) {
  LALStatus status = blank_status;	/* initialize status */
  
  INT4 j, k; /* temp loop variables: k loops over stacks and j over SFTs in a stack*/

  /* detector, ephemeris and velocity/position vector */
  LALDetector detector, detector2;
  EphemerisData *edat = NULL;
  LIGOTimeGPSVector midTstack, startTsft; 
  REAL8VectorSequence *velStack=NULL, *posStack=NULL; /* velocities and positions at midTstack */
  LALLeapSecFormatAndAcc lsfas = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
  INT4 tmpLeap; 
  REAL8 tObs, tStart, tEnd;
  REAL8  *tStack, tStackAvg; /* duration of each stack */
  REAL8 refTime;

  /* sft related stuff */
  SFTVector *inputSFTs=NULL;  /* vector of SFTtypes and SFTtype is COMPLEX8FrequencySeries */
  SFTVectorSequence stackSFTs; /* sequence of sft vectors -- one for each stack */
  INT4 *mCohSft, nSFTs; /* number of SFTs in each stack and total number of SFTs */
  INT4 nStacks; /* number of stacks -- not necessarily same as uvar_nStacks! */
  INT4 sftlength; /* number of bins in each sft */
  REAL8 deltaF, timeBase; /* frequency resolution of SFTs */
  INT8 sftFminBin; /* first sft bin index */
  INT8 fHoughBinIni, fHoughBinFin, binsHough; /* frequency bins of start and end search frequencies */
  REAL8 deltaFstack; /* frequency resolution of Fstat calculation */

  /* LALdemod related stuff */
  REAL8FrequencySeriesVector FstatVect; /* Fstatistic vectors for each stack */
  FstatStackParams FstatPar;

  /* hough variables */
  HOUGHPeakGramVector pgV;
  HoughParams houghPar;
  HoughCandidates houghCand;

  /* variables for logging */
  CHAR *fnamelog=NULL;
  FILE *fpLog = NULL;
  CHAR *logstr=NULL; 

  REAL8 fStatMax = 0.0;

  /* user variables */
  BOOLEAN uvar_help; /* true if -h option is given */
  BOOLEAN uvar_log; /* logging done if true */
  REAL8 uvar_alpha, uvar_delta;  /* sky-location angles */
  REAL8 uvar_fdot; /* first spindown value */
  REAL8 uvar_fStart, uvar_fBand;
  REAL8 uvar_FstatThr; /* threshold of Fstat to select peaks */
  REAL8 uvar_houghThr; /* threshold on hough number count to select candidates */
  INT4 uvar_ifo, uvar_ifo2, uvar_blocksRngMed, uvar_nStacks, uvar_Dterms;
  REAL8 uvar_refTime;
  INT4 uvar_SSBprecision = SSBPREC_RELATIVISTIC;
  INT4 uvar_nfSize;
  CHAR *uvar_earthEphemeris=NULL;
  CHAR *uvar_sunEphemeris=NULL;
  CHAR *uvar_sftDir=NULL;
  CHAR *uvar_sftDir2=NULL;
  CHAR *uvar_fnameout=NULL;

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
  uvar_printMaps = FALSE;
  uvar_printStats = FALSE;
  uvar_nStacks = NSTACKS;
  uvar_Dterms = DTERMS;
  uvar_alpha = ALPHA;
  uvar_delta = DELTA;
  uvar_fdot = FDOT;
  uvar_fStart = FSTART;
  uvar_fBand = FBAND;
  uvar_ifo = IFO;
  uvar_ifo2 = IFO;
  uvar_blocksRngMed = BLOCKSRNGMED;
  uvar_nfSize = NFSIZE;
  uvar_FstatThr = FSTATTHRESHOLD;
  uvar_earthEphemeris = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_earthEphemeris,EARTHEPHEMERIS);

  uvar_sunEphemeris = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_sunEphemeris,SUNEPHEMERIS);

  uvar_sftDir = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_sftDir,SFTDIRECTORY);

  /* do not set default for sftDir2 -- use only if user specifies */
  uvar_sftDir2 = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_sftDir2,SFTDIRECTORY);

  uvar_fnameout = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_fnameout,FNAMEOUT);

  /* register user input variables */
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "help",            'h', UVAR_HELP,     "Print this message",            &uvar_help),            &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "log",              0,  UVAR_OPTIONAL, "Write log file",                &uvar_log),             &status);  
  LAL_CALL( LALRegisterINTUserVar(    &status, "ifo",             'i', UVAR_OPTIONAL, "Detector GEO(1) LLO(2) LHO(3)", &uvar_ifo ),            &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "ifo2",             0,  UVAR_OPTIONAL, "Detector for follow up stage",  &uvar_ifo2 ),           &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "nStacks",         'N', UVAR_OPTIONAL, "Number of stacks",              &uvar_nStacks ),        &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "Dterms",           0,  UVAR_OPTIONAL, "For Dirichlet Kernel approx.",  &uvar_Dterms ),         &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed",    'w', UVAR_OPTIONAL, "RngMed block size",             &uvar_blocksRngMed),    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "earthEphemeris",  'E', UVAR_OPTIONAL, "Earth Ephemeris file",          &uvar_earthEphemeris),  &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sunEphemeris",    'S', UVAR_OPTIONAL, "Sun Ephemeris file",            &uvar_sunEphemeris),    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftDir",           0,  UVAR_OPTIONAL, "SFT Directory",                 &uvar_sftDir),          &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftDir2",          0,  UVAR_OPTIONAL, "SFT Directory for follow up",   &uvar_sftDir2),         &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "fnameout",        'o', UVAR_OPTIONAL, "Output basefileneme",           &uvar_fnameout),        &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "alpha",            0,  UVAR_OPTIONAL, "Right ascension",               &uvar_alpha),           &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "delta",            0,  UVAR_OPTIONAL, "Declination",                   &uvar_delta),           &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fStart",          'f', UVAR_OPTIONAL, "Start search frequency",        &uvar_fStart),          &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fdot",             0,  UVAR_OPTIONAL, "Spindown parameter",            &uvar_fdot),            &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "FstatThr",         0,  UVAR_OPTIONAL, "Threshold on Fstatistic",       &uvar_FstatThr),        &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "houghThr",         0,  UVAR_OPTIONAL, "Hough number count threshold",  &uvar_houghThr),        &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "refTime",          0,  UVAR_OPTIONAL, "Reference time for pulsar par", &uvar_refTime),         &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "SSBprecision",	   0,  UVAR_DEVELOPER,"Precision for SSB transform.",  &uvar_SSBprecision),    &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "nfSize",           0,  UVAR_DEVELOPER,"No.of LUTs to keep in memory",  &uvar_nfSize),          &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printMaps",        0,  UVAR_DEVELOPER,"Print Hough maps",              &uvar_printMaps),       &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printStats",       0,  UVAR_DEVELOPER,"Print Hough map statistics",    &uvar_printStats),      &status);  

  /* read all command line variables */
  LAL_CALL( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 

  /* some very basic sanity checks on user vars */
  if ( uvar_nStacks < 1) {
    fprintf(stderr, "Invalid number of stacks\n");
    exit(1);
  }
  if ( uvar_blocksRngMed < 1 ) {
    fprintf(stderr, "Invalid Running Median block size\n");
    exit(1);
  }

  if ( uvar_FstatThr < 0 ) {
    fprintf(stderr, "Invalid value of Fstatistic threshold\n");
    exit(1);
  }

  /* set detector */
  if (uvar_ifo ==1) detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if (uvar_ifo ==2) detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  if (uvar_ifo ==3) detector = lalCachedDetectors[LALDetectorIndexLHODIFF];

  if (uvar_ifo2 ==1) detector2 = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if (uvar_ifo2 ==2) detector2 = lalCachedDetectors[LALDetectorIndexLLODIFF];
  if (uvar_ifo2 ==3) detector2 = lalCachedDetectors[LALDetectorIndexLHODIFF];

  /* write the log file */
  if (uvar_log) {

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
    } /* end of user var reading */
  } /* end of logging */

  /*------------- read sfts and set up sft timestamp vector ----------*/
  {
    CHAR *tempDir;
    REAL8 doppWings, fmin, fmax;

    tempDir = (CHAR *)LALMalloc(512*sizeof(CHAR));
    strcpy(tempDir, uvar_sftDir);
    strcat(tempDir, "/*SFT*.*");

    doppWings = (uvar_fStart + uvar_fBand) * VTOT;    
    fmin = uvar_fStart - doppWings;
    fmax = uvar_fStart + uvar_fBand + doppWings;
    LAL_CALL( LALReadSFTfiles ( &status, &inputSFTs, fmin, fmax, 
				uvar_blocksRngMed + uvar_Dterms, 
				tempDir), &status); 

    /* normalize sfts */
    LAL_CALL( LALNormalizeSFTVect (&status, inputSFTs, uvar_blocksRngMed, 0), &status);

    LALFree(tempDir);
    
    /* set other sft parameters */
    nSFTs = inputSFTs->length;
    sftlength = inputSFTs->data->data->length;
    deltaF = inputSFTs->data->deltaF;
    timeBase = 1.0/deltaF;
    sftFminBin = floor( timeBase * inputSFTs->data->f0 + 0.5);

    /* get start time of sfts */
    startTsft.length = nSFTs;
    startTsft.data = (LIGOTimeGPS *)LALMalloc( nSFTs * sizeof(LIGOTimeGPS));
    for (j=0; j<nSFTs; j++) 
      startTsft.data[j] = inputSFTs->data[j].epoch;

    /* calculate start and end times and tobs */
    LAL_CALL( LALGPStoFloat ( &status, &tStart, startTsft.data), &status);
    LAL_CALL( LALGPStoFloat ( &status, &tEnd, startTsft.data + nSFTs - 1), &status);
    tEnd += timeBase;
    tObs = tEnd - tStart;

    /* set reference time for pular parameters */
    if ( LALUserVarWasSet(&uvar_refTime)) 
      refTime = uvar_refTime;
    else
      refTime = tStart;

  } /* end of sft reading block */



  /*------------- set up stacks -----------------*/

  if (uvar_nStacks > nSFTs) {
    fprintf(stderr, "invalid number of stacks...exiting\n");
    exit(1);
  }

  /* set up the stacks */
  /* if sfts are to be split evenly among stacks */
  LAL_CALL( SetUpStacks1( &status, &stackSFTs, inputSFTs, uvar_nStacks), &status);
  /* if time is to be split evenly between stacks */
  /* LAL_CALL( SetUpStacks2( &status, &stackSFTs, inputSFTs, &startTsft, uvar_nStacks), &status); */

  /* set number of stacks -- may be different from uvar_nStacks! */
  nStacks = stackSFTs.length;

  /* set up vector of stack durations */
  tStack = NULL;
  tStack = (REAL8 *)LALMalloc( nStacks * sizeof(REAL8));

  /* set up vector of number of sfts in each stack */
  mCohSft = NULL;
  mCohSft = (INT4 *)LALMalloc( nStacks * sizeof(INT4));

  /* set up vector containing mid times of stacks */    
  midTstack.length = nStacks;
  midTstack.data = (LIGOTimeGPS *)LALMalloc( nStacks * sizeof(LIGOTimeGPS));

  for (k=0; k<nStacks; k++) {
    LIGOTimeGPS tempT1, tempT2;
    INT4 tempInt;
    REAL8 tempF1, tempF2, tempMid;

    /* number of sfts in stack */
    tempInt = stackSFTs.data[k].length;
    mCohSft[k] = tempInt;

    /* duration of each stack */
    tempT1 = stackSFTs.data[k].data[0].epoch;
    tempT2 = stackSFTs.data[k].data[tempInt-1].epoch;
    LAL_CALL ( LALGPStoFloat ( &status, &tempF1, &tempT1), &status);
    LAL_CALL ( LALGPStoFloat ( &status, &tempF2, &tempT2), &status);
    tStack[k] = tempF2 + timeBase - tempF1;
    
    /* mid timestamp of each stack */
    tempMid = tempF1 + 0.5 * tStack[k];
    LAL_CALL ( LALFloatToGPS ( & status, midTstack.data + k, &tempMid), &status);
  }

  /* use stacks info to calculate search frequencies */
  /* Fstat is calculated at the frequency resolutions of the stacks
     here the stacks may be of different durations so we take the average
     This is valid if the stack durations are not very different which 
     will, hopefully, be true */
  tStackAvg = 0.0;
  for (k=0; k<nStacks; k++)
    tStackAvg += tStack[k];
  tStackAvg /= nStacks;
  deltaFstack = 1.0/tStackAvg;

 
  /*------------- calculate velocity and position for each stack ------------*/
  /* setting of ephemeris info */ 
  edat = (EphemerisData *)LALMalloc(sizeof(EphemerisData));
  (*edat).ephiles.earthEphemeris = uvar_earthEphemeris;
  (*edat).ephiles.sunEphemeris = uvar_sunEphemeris;

  /* create velocity and position vectors */
  {
    CreateVectorSequenceIn createPar;
    createPar.length = nStacks; /* number of vectors */
    createPar.vectorLength = 3; /* length of each vector */
    LAL_CALL( LALDCreateVectorSequence ( &status,  &velStack, &createPar), &status);
    LAL_CALL( LALDCreateVectorSequence ( &status,  &posStack, &createPar), &status);
  }
    
  /* Leap seconds for the first timestamp */   
  LAL_CALL( LALLeapSecs(&status, &tmpLeap, midTstack.data, &lsfas), &status);
  (*edat).leap = (INT2)tmpLeap;
  
  /* read in ephemeris data */
  LAL_CALL( LALInitBarycenter( &status, edat), &status);
  
  /* calculate detector velocity and position at mid time of stacks*/
  /* maybe calculate average over stack as well? */
  for (k=0; k<nStacks; k++){
    LAL_CALL (LALDetectorVel ( &status, velStack->data + 3*k, midTstack.data + k, detector, edat), &status);
    LAL_CALL (LALDetectorPos ( &status, posStack->data + 3*k, midTstack.data + k, detector, edat), &status);
  }



  /*------------- calculate F statistic for each stack --------------*/

  /* set up parameters for Fstat calculation */
  FstatPar.nStacks = nStacks;
  FstatPar.tsStack = &midTstack;
  FstatPar.tStackAvg = tStackAvg;
  FstatPar.fBand = uvar_fBand;
  FstatPar.fStart = uvar_fStart;
  FstatPar.nfSizeCylinder = uvar_nfSize;
  FstatPar.mCohSft = mCohSft;
  FstatPar.refTime = refTime;
  FstatPar.SSBprecision = uvar_SSBprecision;
  FstatPar.Dterms = uvar_Dterms;
  FstatPar.detector = detector;
  FstatPar.edat = edat;
  FstatPar.ts = &startTsft;
  FstatPar.alpha = uvar_alpha;
  FstatPar.delta = uvar_delta;
  FstatPar.fdot = NULL;
  LAL_CALL ( LALDCreateVector( &status, &(FstatPar.fdot), 1), &status);
  FstatPar.fdot->data[0] = uvar_fdot;

  /* set up memory for Fstat vectors */
  LAL_CALL(SetUpFstatStack( &status, &FstatVect, &FstatPar), &status);

  /* calculate the Fstatistic */
  LAL_CALL(ComputeFstatStack( &status, &FstatVect, &stackSFTs, &FstatPar), &status);

  /* free sfts */
  LALFree(stackSFTs.data);
  LAL_CALL (LALDestroySFTVector(&status, &inputSFTs),&status );

  /* print fstat vector -- for debugging */
  /*   LAL_CALL( PrintFstat ( &status, FstatVect.data, uvar_fnameout, 0 ), &status); */



  /*------------ select peaks ------------*/ 

  LAL_CALL( FstatVectToPeakGram( &status, &pgV, &FstatVect, uvar_FstatThr), &status);

  /* free Fstat */
  for(k=0; k<nStacks; k++) {
    LALFree(FstatVect.data[k].data->data);
    LALFree(FstatVect.data[k].data);
  }
  LALFree(FstatVect.data);


  /*--------------- calculate Hough map and get candidates ---------------*/
  /* start and end bin for calculating hough map */
  /* these are just what the user specified */
  binsHough = floor( tStackAvg * uvar_fBand );
  fHoughBinIni = floor( tStackAvg * uvar_fStart + 0.5);
  fHoughBinFin = fHoughBinIni + binsHough - 1;

  /* set up the Hough parameters */
  if ( ! LALUserVarWasSet(&uvar_houghThr))
    uvar_houghThr = 0.65 * nStacks;
  houghPar.outBaseName = uvar_fnameout;
  houghPar.houghThr = uvar_houghThr;
  houghPar.tStart = tStart;
  houghPar.fBinIni = fHoughBinIni;
  houghPar.fBinFin = fHoughBinFin;
  houghPar.nfSizeCylinder = uvar_nfSize;
  houghPar.detector = detector;
  houghPar.ts = &midTstack;
  houghPar.vel = velStack;
  houghPar.pos = posStack;
  houghPar.alpha = uvar_alpha;
  houghPar.delta = uvar_delta;
  houghPar.fdot = NULL;
  LAL_CALL ( LALDCreateVector( &status, &(houghPar.fdot), 1), &status);
  houghPar.fdot->data[0] = uvar_fdot;

  /* allocate memory for candidates structure */
  houghCand.length = 5000; /* a starting value -- use realloc if this is insufficient */
  houghCand.nCandidates = 0; /* initialization */
  houghCand.freq = (REAL8 *)LALMalloc( houghCand.length * sizeof(REAL8));
  houghCand.alpha = (REAL8 *)LALMalloc( houghCand.length * sizeof(REAL8));
  houghCand.delta = (REAL8 *)LALMalloc( houghCand.length * sizeof(REAL8));
  houghCand.dFreq = (REAL8 *)LALMalloc( houghCand.length * sizeof(REAL8));
  houghCand.dAlpha = (REAL8 *)LALMalloc( houghCand.length * sizeof(REAL8));
  houghCand.dDelta = (REAL8 *)LALMalloc( houghCand.length * sizeof(REAL8));
  houghCand.fdot = (REAL8 *)LALMalloc( houghCand.length * sizeof(REAL8));
  houghCand.dFdot = (REAL8 *)LALMalloc( houghCand.length * sizeof(REAL8));

  /* get candidates */
  LAL_CALL ( ComputeFstatHoughMap ( &status, &houghCand, &pgV, &houghPar), &status);


  /* free memory */
  LALFree(startTsft.data);
  LALFree(mCohSft);

  LAL_CALL (LALDDestroyVector (&status, &(houghPar.fdot)), &status);

  /* free peakgrams */
  for (k=0; k<nStacks; k++) 
    LALFree(pgV.pg[k].peak);
  LALFree(pgV.pg);

  /* free timestamp and Vel/Pos vectors */
  LALFree(midTstack.data);
  LALFree(tStack);
  LAL_CALL( LALDDestroyVectorSequence (&status,  &velStack), &status);
  LAL_CALL( LALDDestroyVectorSequence (&status,  &posStack), &status);


  /* print candidates */
  LAL_CALL ( PrintHoughCandidates ( &status, &houghCand, uvar_fnameout), &status);


  /*------------- Follow up candidates --------------*/

  /* this part is more general than it has to be
     it is meant to be generalized to the case when 
     the number of follow-up stacks is not necessarily 1 */

  /* check if user requested a follow up stage*/
  if ( LALUserVarWasSet(&uvar_sftDir2)) {
    
    INT4 nStacks2;
    
    /* a friendly warning */
    if (lalDebugLevel)
      fprintf(stdout, "Make sure IFOs are set correctly!\n");

    /* read sfts 
       currently reads entire frequency band because it is probably inefficient to
       read sfts separately for each candidate, and range of frequencies
       to be followed up is probably very close to the full range.  However, this
       might have to be changed if there are a very few candidates to be followed up.
    */
    {
      CHAR *tempDir;
      REAL8 doppWings, fmin, fmax;
      
      tempDir = (CHAR *)LALMalloc(512*sizeof(CHAR));
      strcpy(tempDir, uvar_sftDir2);
      strcat(tempDir, "/*SFT*.*");
      
      doppWings = (uvar_fStart + uvar_fBand) * VTOT;    
      fmin = uvar_fStart - doppWings;
      fmax = uvar_fStart + uvar_fBand + doppWings;
      LAL_CALL( LALReadSFTfiles ( &status, &inputSFTs, fmin, fmax, 
				  uvar_blocksRngMed + uvar_Dterms, 
				  tempDir), &status); 
      
      /* normalize sfts */
      LAL_CALL( LALNormalizeSFTVect (&status, inputSFTs, uvar_blocksRngMed, 0), &status);

      LALFree(tempDir);    

      /* set other sft parameters */
      nSFTs = inputSFTs->length;
      sftlength = inputSFTs->data->data->length;
      deltaF = inputSFTs->data->deltaF;
      timeBase = 1.0/deltaF;
      sftFminBin = floor( timeBase * inputSFTs->data->f0 + 0.5);

      /* get start time of sfts */
      startTsft.length = nSFTs;
      startTsft.data = NULL;
      startTsft.data = (LIGOTimeGPS *)LALMalloc( nSFTs * sizeof(LIGOTimeGPS));
      for (j=0; j<nSFTs; j++) 
	startTsft.data[j] = inputSFTs->data[j].epoch;
      
      /* calculate start and end times and tobs */
      LAL_CALL( LALGPStoFloat ( &status, &tStart, startTsft.data), &status);
      LAL_CALL( LALGPStoFloat ( &status, &tEnd, startTsft.data + nSFTs - 1), &status);
      tEnd += timeBase;
      tObs = tEnd - tStart;

    } /* end sft reading block */

    /* there is just one stack now */
    nStacks2 = 1;
    LAL_CALL( SetUpStacks1( &status, &stackSFTs, inputSFTs, nStacks2), &status);
    
    /* set up vector of stack durations */
    tStack = NULL;
    tStack = (REAL8 *)LALMalloc( nStacks2 * sizeof(REAL8));
    
    /* set up vector of number of sfts in each stack */
    mCohSft = NULL;
    mCohSft = (INT4 *)LALMalloc( nStacks2 * sizeof(INT4));

    /* set up vector containing mid times of stacks */    
    midTstack.length = nStacks2;
    midTstack.data = (LIGOTimeGPS *)LALMalloc( nStacks2 * sizeof(LIGOTimeGPS));

    for (k=0; k<nStacks2; k++) {

      LIGOTimeGPS tempT1, tempT2;
      INT4 tempInt;
      REAL8 tempF1, tempF2, tempMid;

      /* number of sfts in stack */
      tempInt = stackSFTs.data[k].length;
      mCohSft[k] = tempInt;
      
      /* duration of each stack */
      tempT1 = stackSFTs.data[k].data[0].epoch;
      tempT2 = stackSFTs.data[k].data[tempInt-1].epoch;
      LAL_CALL ( LALGPStoFloat ( &status, &tempF1, &tempT1), &status);
      LAL_CALL ( LALGPStoFloat ( &status, &tempF2, &tempT2), &status);
      tStack[k] = tempF2 + timeBase - tempF1;
      
      /* mid timestamp of each stack */
      tempMid = tempF1 + 0.5 * tStack[k];
      LAL_CALL ( LALFloatToGPS ( & status, midTstack.data + k, &tempMid), &status);
    }
    
    /* use stacks info to calculate search frequencies */
    /* Fstat is calculated at the frequency resolutions of the stacks
       here the stacks may be of different durations so we take the average
       This is valid if the stack durations are not very different which 
       will, hopefully, be true */
    tStackAvg = 0.0;
    for (k=0; k<nStacks2; k++)
      tStackAvg += tStack[k];
    tStackAvg /= nStacks2;
    deltaFstack = 1.0/tStackAvg;

    /* set up parameters for Fstat calculation */
    FstatPar.nStacks = nStacks2;
    FstatPar.tsStack = &midTstack;
    FstatPar.tStackAvg = tStackAvg;
    FstatPar.nfSizeCylinder = uvar_nfSize;
    FstatPar.mCohSft = mCohSft;
    FstatPar.refTime = refTime;
    FstatPar.SSBprecision = uvar_SSBprecision;
    FstatPar.Dterms = uvar_Dterms;
    FstatPar.detector = detector2;
    FstatPar.edat = edat;
    FstatPar.ts = &startTsft;

    /*------------- calculate the Fstatistic ---------------*/
    /* loop over candidates */
    for ( j=0; j < houghCand.nCandidates; j++) {

      REAL8 tempMax;
      
      FstatPar.fStart = houghCand.freq[j] - 0.5 * houghCand.dFreq[j];
      FstatPar.fBand = houghCand.dFreq[j];
      FstatPar.alpha = houghCand.alpha[j];
      FstatPar.delta = houghCand.delta[j];
      FstatPar.fdot->data[0] = houghCand.fdot[j];

      LAL_CALL(SetUpFstatStack( &status, &FstatVect, &FstatPar), &status);

      LAL_CALL(ComputeFstatStack( &status, &FstatVect, &stackSFTs, &FstatPar), &status);

      LAL_CALL ( GetLoudestFstat ( &status, &tempMax, FstatVect.data), &status);

      if ( tempMax > fStatMax )
	fStatMax = tempMax;
      
      /*     LAL_CALL( PrintFstat ( &status, FstatVect.data, uvar_fnameout, 0 ), &status);  */
      
      /* free Fstat */
      for(k=0; k<nStacks2; k++) {
	LALFree(FstatVect.data[k].data->data);
	LALFree(FstatVect.data[k].data);
      }
      LALFree(FstatVect.data);
      
    } /* end loop over candidates */

    /* print loudest event */
    fprintf(stdout, "Loudest Candidate has F = %e\n", fStatMax);

    /* we are done with follow-up */
    /*free sfts */
    LALFree(startTsft.data);
    LALFree(stackSFTs.data);
    LAL_CALL (LALDestroySFTVector(&status, &inputSFTs),&status );

    LALFree(midTstack.data);
    LALFree(tStack);
    LALFree(mCohSft);
 } /* end block for follow-up stage */ 




  /*------------ free all remaining memory -----------*/
  /* free candidates */
  LALFree(houghCand.freq);
  LALFree(houghCand.dFreq);
  LALFree(houghCand.alpha);
  LALFree(houghCand.dAlpha);
  LALFree(houghCand.delta);
  LALFree(houghCand.dDelta);
  LALFree(houghCand.fdot);
  LALFree(houghCand.dFdot);
  
  LAL_CALL (LALDDestroyVector (&status, &(FstatPar.fdot)), &status); 

  /* free ephemeris */
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);

  LAL_CALL (LALDestroyUserVars(&status), &status);  

  LALCheckMemoryLeaks();

  return 0;
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
			const SFTVectorSequence *stackSFTs, 
			FstatStackParams *params)
{
  /* stuff copied from params */
  INT4 nStacks = params->nStacks;
  REAL8 refTime = params->refTime;
  INT4 *mCohSft = params->mCohSft;
  REAL8Vector *fdot = params->fdot;
  INT4 nfSizeCylinder = params->nfSizeCylinder;
  REAL8 fStart = params->fStart;
  REAL8 fBand = params->fBand;
  REAL8 tStackAvg = params->tStackAvg;
  REAL8 deltaF = 1.0/tStackAvg;

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
			  const HOUGHPeakGramVector *pgV, /* peakgram vector */
			  HoughParams *params)
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
    fpStats = fopen(fileStats, "w");
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
	  TRY(GetHoughCandidates( status->statusPtr, out, &ht, &patch, 
				  &parDem, houghThr), status);


	  /* calculate statistics and histogram */
	  if ( uvar_printStats ) {
	    TRY( LALHoughStatistics ( status->statusPtr, &stats, &ht), status );
	    TRY( LALStereo2SkyLocation ( status->statusPtr, &sourceLocation, 
					stats.maxIndex[0], stats.maxIndex[1], 
					&patch, &parDem), status);

	    fprintf(fpStats, "%d %f %f %d %d %f %f %f %g \n",
		    iHmap, sourceLocation.alpha, sourceLocation.delta,
		    stats.maxCount, stats.minCount, stats.avgCount,stats.stdDev,
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
			  const REAL8FrequencySeriesVector *FstatVect,
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
		 const SFTVector  *sftVect,
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
		  const SFTVector  *sftVect,
		  const LIGOTimeGPSVector *ts,
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


/** Prints Fstatistic values as a function of frequency to a specified output file */ 
void PrintFstat( LALStatus *status,
		 REAL8FrequencySeries *Fstat, 
		 CHAR *fname, 
		 INT4 stackIndex)
{

  FILE *fp=NULL;
  INT4 k, length;
  REAL8 freq, deltaF;
  CHAR filename[256], filenumber[16]; 

  INITSTATUS( status, "PrintFstat", rcsid );
  ATTATCHSTATUSPTR (status);

  length = Fstat->data->length;
  freq = Fstat->f0;
  deltaF = Fstat->deltaF;

  strcpy ( filename, fname);
  sprintf ( filenumber, ".%d", stackIndex);
  strcat ( filename, filenumber);
  
  fp = fopen(filename, "w");

  for (k=0; k<length; k++) {
    fprintf(fp, "%e   %e\n", freq, Fstat->data->data[k]);
    freq += deltaF;
  }

  fclose(fp);

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

  /* replace this by an assert */
  /*   if ( !fp ){   */
  /*     fprintf(stderr,"Unable to find file %s\n",filename); */
  /*     return DRIVEHOUGHFSTAT_EFILE;  */
  /*   } */

  ySide= ht->ySide;
  xSide= ht->xSide;

  for(k=ySide-1; k>=0; --k){
    for(i=0;i<xSide;++i){
      fprintf( fp ," %d", ht->map[k*xSide +i]);
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
			const HOUGHMapTotal *ht,
			const HOUGHPatchGrid  *patch,
			const HOUGHDemodPar   *parDem,
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
	  houghCand->freq = (REAL8 *)LALRealloc( houghCand->freq, 
						 houghCand->length * sizeof(REAL8));
	  houghCand->alpha = (REAL8 *)LALRealloc( houghCand->alpha, 
						 houghCand->length * sizeof(REAL8));
	  houghCand->delta = (REAL8 *)LALRealloc( houghCand->delta, 
						 houghCand->length * sizeof(REAL8));
	  houghCand->dFreq = (REAL8 *)LALRealloc( houghCand->dFreq, 
						 houghCand->length * sizeof(REAL8));
	  houghCand->dAlpha = (REAL8 *)LALRealloc( houghCand->dAlpha, 
						 houghCand->length * sizeof(REAL8));
	  houghCand->dDelta = (REAL8 *)LALRealloc( houghCand->dDelta, 
						 houghCand->length * sizeof(REAL8));
	  houghCand->fdot = (REAL8 *)LALRealloc( houghCand->fdot, 
						 houghCand->length * sizeof(REAL8));
	  houghCand->dFdot = (REAL8 *)LALRealloc( houghCand->dFdot, 
						  houghCand->length * sizeof(REAL8));
	} /* end of reallocs */

	houghCand->freq[nCandidates] = f0;
	houghCand->dFreq[nCandidates] = deltaF;
	
	houghCand->fdot[nCandidates] = fdot;
	houghCand->dFdot[nCandidates] = dFdot;

	/* get sky location of pixel */
	TRY( LALStereo2SkyLocation (status->statusPtr, &sourceLocation, 
				    j, i, patch, parDem), status);

	houghCand->alpha[nCandidates] = sourceLocation.alpha;
	houghCand->delta[nCandidates] = sourceLocation.delta;

	houghCand->dAlpha[nCandidates] = patchSizeX / ((REAL8)xSide);
	houghCand->dDelta[nCandidates] = patchSizeY / ((REAL8)ySide);

	/* increment candidate count */
	houghCand->nCandidates += 1;

      } /* end if statement */
    } /* end loop over xSide */
  } /* end loop over ySide */
  DETATCHSTATUSPTR (status);
  RETURN(status);
}



/** Print Hough candidates */
void PrintHoughCandidates(LALStatus *status,
			  HoughCandidates *in,
			  CHAR *fname)
{
  FILE  *fp=NULL;   /* Output file */
  CHAR filename[256]; 
  INT4 k;

  INITSTATUS( status, "GetHoughCandidates", rcsid );
  ATTATCHSTATUSPTR (status);

  strcpy(filename, fname);
  strcat(filename, ".cand");
  
  fp = fopen(filename, "w");

  for (k=0; k < in->nCandidates; k++)
    fprintf(fp, "%e   %e   %g   %g   %g   %g   %e   %e\n", in->freq[k], 
	    in->dFreq[k], in->alpha[k], in->dAlpha[k], in->delta[k], in->dDelta[k],
	    in->fdot[k], in->dFdot[k]);
  
  fclose(fp);

  DETATCHSTATUSPTR (status);
  RETURN(status);
}



/** Utility function for getting largest value in a Fstat vector */
void GetLoudestFstat(LALStatus *status,
		     REAL8 *max,
		     REAL8FrequencySeries *Fstat)
{

  REAL8 temp = 0;
  INT4 k, length;

  INITSTATUS( status, "GetHoughCandidates", rcsid );
  ATTATCHSTATUSPTR (status);

  length = Fstat->data->length;

  for (k=0; k<length; k++) {
    if ( Fstat->data->data[k] > temp ) 
      temp = Fstat->data->data[k];
  }

  *max = temp;

  DETATCHSTATUSPTR (status);
  RETURN(status);
}

/** Sets threshold on Fstat vector and fills up frequency and deltaF values in candidate list. 
    Important -- Currently this function does not get the sky-loctaion and spindown
    values for the candidates because the input does not have this information.  This 
    is to be fixed as soon as possible!
*/
void GetFstatCandidates( LALStatus *status,
			 HoughCandidates *cand,
			 const REAL8FrequencySeries *in,
			 REAL8 FstatThr)
{
  INT4 k, length, nCandidates;
  REAL8 deltaF, f0;

  INITSTATUS( status, "GetFstatCandidates", rcsid );
  ATTATCHSTATUSPTR (status);

  length = in->data->length;
  deltaF = in->deltaF;
  f0 = in->f0;

  for ( k=0; k<length; k++) {
    
	nCandidates = cand->nCandidates;
	
	/* if there isn't enough memory then realloc */
	if ( nCandidates >= cand->length ) {
	  cand->length += 5000;
	  cand->freq = (REAL8 *)LALRealloc( cand->freq, 
					    cand->length * sizeof(REAL8));
	  cand->alpha = (REAL8 *)LALRealloc( cand->alpha, 
					     cand->length * sizeof(REAL8));
	  cand->delta = (REAL8 *)LALRealloc( cand->delta, 
					     cand->length * sizeof(REAL8));
	  cand->dFreq = (REAL8 *)LALRealloc( cand->dFreq, 
					     cand->length * sizeof(REAL8));
	  cand->dAlpha = (REAL8 *)LALRealloc( cand->dAlpha, 
					      cand->length * sizeof(REAL8));
	  cand->dDelta = (REAL8 *)LALRealloc( cand->dDelta, 
					      cand->length * sizeof(REAL8));
	  cand->fdot = (REAL8 *)LALRealloc( cand->fdot, 
					    cand->length * sizeof(REAL8));
	  cand->dFdot = (REAL8 *)LALRealloc( cand->dFdot, 
					     cand->length * sizeof(REAL8)); 
	} /* end of reallocs */


	if ( in->data->data[k] > FstatThr ) {

	  cand->freq[nCandidates] = f0 + k*deltaF;
	  cand->dFreq[nCandidates] = deltaF;

	  cand->nCandidates += 1;
	}
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
  fp=fopen(filename,"w");

  for (i=0; i < hist->length; i++)
    fprintf(fp,"%d  %d\n", i, hist->data[i]);
  
  fclose( fp );  
  
  DETATCHSTATUSPTR (status);
  RETURN(status);

}
