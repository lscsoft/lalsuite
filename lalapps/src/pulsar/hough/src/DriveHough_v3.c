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
 * \file
 * \ingroup pulsarApps
 * \author Badri Krishnan, Alicia Sintes
 * \brief Driver code for performing Hough transform search on non-demodulated
 * data.
 *
 * History:   Created by Sintes and Krishnan July 04, 2003
 * Modifications for S4 January 2006
 *
 * \par Description
 *
 * This is the main driver for the Hough transform routines.
 * It takes as input a set of SFTs and search parameters and
 * outputs the number counts using the Hough transform.
 *
 * \par User input
 *
 * The user inputs the following parameters:
 *
 * - Search frequency range
 *
 * - A file containing list of skypatches to search over.  For each skypatch,
 * the information is:
 * - RA and dec of skypatch center.
 * - Size in RA and dec.
 *
 * - Location of Directory containing the SFTs (can be either v1 or v2).
 *
 * - Interferometer .
 * -- Not required for v2 SFTs.
 *
 * - Location of output directory and basename of output files.
 *
 * - Block size of running median for estimating psd.
 *
 * - The parameter nfSizeCylinder which determines the range of spindown parameters
 * to be searched over.
 *
 * - Boolean variable for deciding if the SFTs should be inverse noise weighed.
 *
 * - Boolean variable for deciding whether amplitude modulation weights should be used.
 *
 * - Boolean variables for deciding whether the Hough maps, the statistics and list of
 * events above a threshold should be written
 *
 * /par Output
 *
 * The output is written in several sub-directories of the specified output directory.  The
 * first two items are default while the rest are written according to the user input:
 *
 * - A directory called logfiles records the user input, contents of the skypatch file
 * and cvs tags contained in the executable
 *
 * - A directory called nstarfiles containing the loudest event for each search frequency
 * bin maximised over all sky-locations and spindown parameters.  An event is said to be
 * the loudest if it has the maximum significance defined as: (number count - mean)/sigma.
 *
 * - A directory for each skypatch containing the number count statistics, the histogram,
 * the list of events, and the Hough maps
 */


#include "./DriveHoughColor.h"

#ifdef TIMING
#include "./timer/cycle_counter/Intel/GCC/cycle_counter.h"
#endif

/* globals, constants and defaults */


/* boolean global variables for controlling output */
BOOLEAN uvar_printEvents, uvar_printTemplates, uvar_printMaps, uvar_printStats, uvar_printSigma;

/* #define EARTHEPHEMERIS "./earth05-09.dat" */
/* #define SUNEPHEMERIS "./sun05-09.dat"    */

#define EARTHEPHEMERIS "./earth00-04.dat"
#define SUNEPHEMERIS "./sun00-04.dat"

#define ACCURACY 0.00000001 /* of the velocity calculation -- irrelevant */

#define MAXFILES 3000 /* maximum number of files to read in a directory */
#define MAXFILENAMELENGTH 256 /* maximum # of characters  of a SFT filename */

#define SFTDIRECTORY "/local_data/badkri/fakesfts-multi"
/* #define SFTDIRECTORY "/nfs/morbo/geo600/hannover/sft/S4-LIGO/sft_1800.20050512.S4/S4-L1.1800-sft" */

#define DIROUT "./outHM1/"      /* prefix file output */
#define BASENAMEOUT "HM1"

#define THRESHOLD 1.6 /* thresold for peak selection, with respect to the
                              the averaged power in the search band */
#define FALSEALARM 1.0e-9 /* Hough false alarm for candidate selection */
#define SKYFILE "./sky1"      
#define F0 205.0          /*  frequency to build the LUT and start search */
#define FBAND 0.2          /* search frequency band  (in Hz) */
#define NFSIZE  21 /* n-freq. span of the cylinder, to account for spin-down search */
#define BLOCKSRNGMED 101 /* Running median window size */

#define TRUE (1==1)
#define FALSE (1==0)


void PrintLogFile (LALStatus *status, CHAR *dir, CHAR *basename, CHAR *skyfile, CHAR *linefile, CHAR *executable );

int PrintHistogram(UINT8Vector *, CHAR *);

int PrintHmap2file(HOUGHMapTotal *ht, CHAR *fnameOut, INT4 iHmap);

int PrintHmap2m_file(HOUGHMapTotal *ht, CHAR *fnameOut, INT4 iHmap);

void PrintHoughEvents (LALStatus *status, FILE *fpEvents, INT4 houghThreshold, HOUGHMapTotal *ht, HOUGHPatchGrid *patch, HOUGHDemodPar *parDem);

void PrintnStarFile (LALStatus *status, HoughSignificantEventVector *eventVec, CHAR *dirname, CHAR *basename);

/***********************************************/

int main(int argc, char *argv[]){

  /* LALStatus pointer */
  static LALStatus  status;  

  /* detector */
  static LALDetector  *detector;

  /* time and velocity vectors */
  static LIGOTimeGPSVector   *timeV=NULL;
  static REAL8Cart3CoorVector velV;
  static REAL8Vector         timeDiffV;
  REAL8  tObs;

  /* standard pulsar sft types */ 
  SFTVector *inputSFTs = NULL;

  /* vector of weights */
  REAL8Vector weightsV, weightsNoise;

  /* ephemeris */
  EphemerisData    *edat=NULL;

  /* hough structures */
  static HOUGHptfLUTVector   lutV; /* the Look Up Table vector*/
  static HOUGHPeakGramVector pgV;  /* vector of peakgrams */
  static PHMDVectorSequence  phmdVS;  /* the partial Hough map derivatives */
  static UINT8FrequencyIndexVector freqInd; /* for trajectory in time-freq plane */
  static HOUGHResolutionPar parRes;   /* patch grid information */
  static HOUGHPatchGrid  patch;   /* Patch description */ 
  static HOUGHParamPLUT  parLut;  /* parameters needed to build lut  */
  static HOUGHDemodPar   parDem;  /* demodulation parameters or  */
  static HOUGHSizePar    parSize; 
  static HOUGHMapTotal   ht;   /* the total Hough map */
  static UINT8Vector     hist; /* histogram of number counts for a single map */
  static UINT8Vector     histTotal; /* number count histogram for all maps */
  static HoughStats      stats;

  /* skypatch info */
  REAL8  *skyAlpha, *skyDelta, *skySizeAlpha, *skySizeDelta; 
  INT4   nSkyPatches, skyCounter=0; 

  /* output filenames and filepointers */
  CHAR   filehisto[256]; 
  CHAR   filestats[256]; 
  CHAR   fileEvents[256];
  CHAR   fileTemplates[256];
  CHAR   fileMaps[256];
  CHAR   fileSigma[256];
  FILE   *fpTemplates = NULL;
  FILE   *fpEvents = NULL;
  FILE   *fpStats = NULL;
  FILE   *fpSigma = NULL;

  /* the maximum number count */
  static HoughSignificantEventVector nStarEventVec;

  /* miscellaneous */
  INT4   houghThreshold=0, iHmap, nSpin1Max;
  UINT4  mObsCoh;
  INT8   f0Bin, fLastBin, fBin;
  REAL8  alpha, delta, timeBase, deltaF, f1jump;
  REAL8  patchSizeX, patchSizeY;
  UINT2  xSide, ySide;
  UINT2  maxNBins, maxNBorders;

  /* user input variables */
  BOOLEAN  uvar_help, uvar_weighAM, uvar_weighNoise;
  INT4     uvar_blocksRngMed, uvar_nfSizeCylinder, uvar_maxBinsClean;
  REAL8    uvar_f0, uvar_peakThreshold, uvar_houghFalseAlarm, uvar_fSearchBand;
  CHAR     *uvar_earthEphemeris=NULL;
  CHAR     *uvar_sunEphemeris=NULL;
  CHAR     *uvar_sftDir=NULL;
  CHAR     *uvar_dirnameOut=NULL;
  CHAR     *uvar_fbasenameOut=NULL;
  CHAR     *uvar_skyfile=NULL;
  CHAR     *uvar_ifo=NULL;
  CHAR     *uvar_linefile=NULL;

#ifdef TIMING
  unsigned long long start, stop;
#endif
  
#ifdef TIMING
  start = realcc();
#endif

  

  /* Set up the default parameters */
  /* LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;
  
  uvar_help = FALSE;
  uvar_weighAM = FALSE;
  uvar_weighNoise = FALSE;
  uvar_blocksRngMed = BLOCKSRNGMED;
  uvar_nfSizeCylinder = NFSIZE;
  uvar_maxBinsClean = 100;
  uvar_f0 = F0;
  uvar_fSearchBand = FBAND;
  uvar_peakThreshold = THRESHOLD;
  uvar_houghFalseAlarm = FALSEALARM;
  uvar_printEvents = FALSE;
  uvar_printTemplates = FALSE;
  uvar_printMaps = FALSE;
  uvar_printStats = FALSE;
  uvar_printSigma = FALSE;
  uvar_earthEphemeris = (CHAR *)LALCalloc(512, sizeof(CHAR));
  strcpy(uvar_earthEphemeris,EARTHEPHEMERIS);

  uvar_sunEphemeris = (CHAR *)LALCalloc(512, sizeof(CHAR));
  strcpy(uvar_sunEphemeris,SUNEPHEMERIS);

  uvar_sftDir = (CHAR *)LALCalloc(512, sizeof(CHAR));
  strcpy(uvar_sftDir,SFTDIRECTORY);

  uvar_dirnameOut = (CHAR *)LALCalloc(512, sizeof(CHAR));
  strcpy(uvar_dirnameOut,DIROUT);

  uvar_fbasenameOut = (CHAR *)LALCalloc(512, sizeof(CHAR));
  strcpy(uvar_fbasenameOut,BASENAMEOUT);

  uvar_skyfile = (CHAR *)LALCalloc(512, sizeof(CHAR));
  strcpy(uvar_skyfile,SKYFILE);

  /* register user input variables */
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "help",             'h', UVAR_HELP,     "Print this message",                  &uvar_help),            &status);  
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "ifo",              'i', UVAR_OPTIONAL, "Detector L1, H1, H2, G1",             &uvar_ifo ),            &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "f0",               'f', UVAR_OPTIONAL, "Start search frequency",              &uvar_f0),              &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fSearchBand",      'b', UVAR_OPTIONAL, "Search frequency band",               &uvar_fSearchBand),     &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "skyfile",           0,  UVAR_OPTIONAL, "Input skypatch file",                 &uvar_skyfile),         &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "peakThreshold",     0,  UVAR_OPTIONAL, "Peak selection threshold",            &uvar_peakThreshold),   &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "weighAM",           0,  UVAR_OPTIONAL, "Use amplitude modulation weights",    &uvar_weighAM),         &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "weighNoise",        0,  UVAR_OPTIONAL, "Use SFT noise weights",               &uvar_weighNoise),      &status);  
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "earthEphemeris",   'E', UVAR_OPTIONAL, "Earth Ephemeris file",                &uvar_earthEphemeris),  &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sunEphemeris",     'S', UVAR_OPTIONAL, "Sun Ephemeris file",                  &uvar_sunEphemeris),    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftDir",           'D', UVAR_OPTIONAL, "SFT Directory",                       &uvar_sftDir),          &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "dirnameOut",       'o', UVAR_OPTIONAL, "Output directory",                    &uvar_dirnameOut),      &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "fbasenameOut",      0,  UVAR_OPTIONAL, "Output file basename",                &uvar_fbasenameOut),    &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printMaps",         0,  UVAR_OPTIONAL, "Print Hough maps",                    &uvar_printMaps),       &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printTemplates",    0,  UVAR_OPTIONAL, "Print templates file",                &uvar_printTemplates),  &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "houghFalseAlarm",   0,  UVAR_OPTIONAL, "Hough false alarm to set threshold",  &uvar_houghFalseAlarm), &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printEvents",       0,  UVAR_OPTIONAL, "Print events above threshold",        &uvar_printEvents),     &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printStats",        0,  UVAR_OPTIONAL, "Print Hough statistics",              &uvar_printStats),      &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printSigma",        0,  UVAR_OPTIONAL, "Print expected number count stdev.",  &uvar_printSigma),      &status);  
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "linefile",          0,  UVAR_OPTIONAL, "list of known lines to clean SFTs",   &uvar_linefile),        &status);  

  /* developer input variables */
  LAL_CALL( LALRegisterINTUserVar(    &status, "nfSizeCylinder",    0, UVAR_DEVELOPER, "Size of cylinder of PHMDs",           &uvar_nfSizeCylinder),  &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed",      0, UVAR_DEVELOPER, "Running Median block size",           &uvar_blocksRngMed),    &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "maxBinsClean",      0, UVAR_DEVELOPER, "Maximum number of bins in cleaning",  &uvar_maxBinsClean),    &status);

  /* read all command line variables */
  LAL_CALL( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 

  /* very basic consistency checks on user input */
  if ( uvar_f0 < 0 ) {
    fprintf(stderr, "start frequency must be positive\n");
    exit(1);
  }

  if ( uvar_fSearchBand < 0 ) {
    fprintf(stderr, "search frequency band must be positive\n");
    exit(1);
  }
 
  if ( uvar_peakThreshold < 0 ) {
    fprintf(stderr, "peak selection threshold must be positive\n");
    exit(1);
  }

  if ( uvar_printEvents ) {
    if ((uvar_houghFalseAlarm > 1.0 ) || (uvar_houghFalseAlarm < 0.0) )  {
      fprintf(stderr, "hough false alarm must be between 0 and 1\n");
      exit(1);
    }
  }


  /* write log file with command line arguments, cvs tags, contents of skypatch file and, if required, contents of linefile*/
  if ( LALUserVarWasSet( &uvar_linefile ) )    
    LAL_CALL( PrintLogFile( &status, uvar_dirnameOut, uvar_fbasenameOut, uvar_skyfile, uvar_linefile, argv[0]), &status);
  else
    LAL_CALL( PrintLogFile( &status, uvar_dirnameOut, uvar_fbasenameOut, uvar_skyfile, NULL, argv[0]), &status);



  /***** start main calculations *****/

  /* read skypatch info */
  {
    FILE   *fpsky = NULL; 
    INT4   r;
    REAL8  temp1, temp2, temp3, temp4;
    
    if ( (fpsky = fopen(uvar_skyfile, "r")) == NULL)
      {
	fprintf(stderr, "Unable to find skyfile %s\n", uvar_skyfile);
	return DRIVEHOUGHCOLOR_EFILE;
      }
        
    nSkyPatches = 0;
    do 
      {
	r = fscanf(fpsky,"%lf%lf%lf%lf\n", &temp1, &temp2, &temp3, &temp4);
	/* make sure the line has the right number of entries or is EOF */
	if (r==4) nSkyPatches++;
      } while ( r != EOF);
    rewind(fpsky);
    
    skyAlpha = (REAL8 *)LALCalloc(nSkyPatches, sizeof(REAL8));
    skyDelta = (REAL8 *)LALCalloc(nSkyPatches, sizeof(REAL8));     
    skySizeAlpha = (REAL8 *)LALCalloc(nSkyPatches, sizeof(REAL8));
    skySizeDelta = (REAL8 *)LALCalloc(nSkyPatches, sizeof(REAL8));     
    
    for (skyCounter = 0; skyCounter < nSkyPatches; skyCounter++)
      {
	r = fscanf(fpsky,"%lf%lf%lf%lf\n", skyAlpha + skyCounter, skyDelta + skyCounter, 
		   skySizeAlpha + skyCounter,  skySizeDelta + skyCounter);
      }
    
    fclose(fpsky);     
  } /* end skyfile reading block */


  /* read sft files and set up weights and nstar vector */
  {
    /* new SFT I/O data types */
    SFTCatalog *catalog = NULL;
    static SFTConstraints constraints;

    CHAR *tempDir;
    REAL8 doppWings, f_min, f_max;
    INT4 length;

    /* set detector constraint */
    constraints.detector = NULL;
    if ( LALUserVarWasSet( &uvar_ifo ) )    
      constraints.detector = XLALGetChannelPrefix ( uvar_ifo );

    /* get sft catalog */
    tempDir = (CHAR *)LALCalloc(512, sizeof(CHAR));
    strcpy(tempDir, uvar_sftDir);
    strcat(tempDir, "/*SFT*.*");
    LAL_CALL( LALSFTdataFind( &status, &catalog, tempDir, &constraints), &status);

    /* exit if catalog has multi ifos and no detector constraint has been set */
    if ( ( XLALCountIFOsInCatalog(catalog) > 1) && !(LALUserVarWasSet( &uvar_ifo ))) {
      fprintf(stderr, "sft catalog contains more than one ifo and no ifo has been specified...exiting\n");
      exit(1);
    }

    /* set detector */
    detector = XLALGetSiteInfo ( catalog->data[0].header.name );

    /* get some sft parameters */
    mObsCoh = catalog->length; /* number of sfts */
    deltaF = catalog->data->header.deltaF;  /* frequency resolution */
    timeBase= 1.0/deltaF; /* coherent integration time */
    f0Bin = floor( uvar_f0 * timeBase + 0.5); /* initial search frequency */
    length =  uvar_fSearchBand * timeBase; /* total number of search bins - 1 */
    fLastBin = f0Bin + length;   /* final frequency bin to be analyzed */

    /* using value of length, allocate memory for most significant event nstar, fstar etc. */
    nStarEventVec.length = length + 1;
    nStarEventVec.event = NULL;
    nStarEventVec.event = (HoughSignificantEvent *)LALCalloc( length+1, sizeof(HoughSignificantEvent));
    /* initialize nstar values -- really unnecessary */
    memset( nStarEventVec.event, 0, length+1);
    
    /* get SFT timestamps */
    LAL_CALL( LALSFTtimestampsFromCatalog(  &status, &timeV, catalog ), &status);  	

    /* total observation time = last timestamp - first timestamp + timeBase*/
    tObs = XLALGPSDiff( timeV->data + mObsCoh - 1, timeV->data ) + timeBase;

    /* add wings for Doppler modulation and running median block size*/
    doppWings = (uvar_f0 + uvar_fSearchBand) * VTOT;    
    f_min = uvar_f0 - doppWings - (uvar_blocksRngMed + uvar_nfSizeCylinder) * deltaF;
    f_max = uvar_f0 + uvar_fSearchBand + doppWings + (uvar_blocksRngMed + uvar_nfSizeCylinder) * deltaF;

    /* read sft files making sure to add extra bins for running median */
    /* read the sfts */
    LAL_CALL( LALLoadSFTs ( &status, &inputSFTs, catalog, f_min, f_max), &status);

    /* clean sfts if required */
    if ( LALUserVarWasSet( &uvar_linefile ) )
      {
	RandomParams *randPar=NULL;
	FILE *fpRand=NULL;
	INT4 seed, ranCount;  

	if ( (fpRand = fopen("/dev/urandom", "r")) == NULL ) {
	  fprintf(stderr,"Error in opening /dev/urandom" ); 
	  exit(1);
	} 

	if ( (ranCount = fread(&seed, sizeof(seed), 1, fpRand)) != 1 ) {
	  fprintf(stderr,"Error in getting random seed" );
	  exit(1);
	}

	LAL_CALL ( LALCreateRandomParams (&status, &randPar, seed), &status );

	/* maxBinsClean is set to 100 by default -- should be large enough */
	LAL_CALL( LALRemoveKnownLinesInSFTVect ( &status, inputSFTs, uvar_maxBinsClean, uvar_blocksRngMed/2, uvar_linefile, randPar), &status);

	LAL_CALL ( LALDestroyRandomParams (&status, &randPar), &status);
	fclose(fpRand);
      } /* end cleaning */

    /* free memory */
    if ( LALUserVarWasSet( &uvar_ifo ) )    
      LALFree( constraints.detector );
    LALFree( tempDir);
    LAL_CALL( LALDestroySFTCatalog( &status, &catalog ), &status);  	

    /* set up weights -- this should be done before normalizing the sfts */
    weightsV.length = mObsCoh;
    weightsV.data = (REAL8 *)LALCalloc(mObsCoh, sizeof(REAL8));

    weightsNoise.length = mObsCoh;
    weightsNoise.data = (REAL8 *)LALCalloc(mObsCoh, sizeof(REAL8));

    /* initialize all weights to unity */
    LAL_CALL( LALHOUGHInitializeWeights( &status, &weightsNoise), &status);
    LAL_CALL( LALHOUGHInitializeWeights( &status, &weightsV), &status);

    /* calculate sft noise weights if required by user */
    if (uvar_weighNoise ) {
      LAL_CALL( LALComputeNoiseWeights( &status, &weightsNoise, inputSFTs, uvar_blocksRngMed, 0), &status); 
      LAL_CALL( LALHOUGHNormalizeWeights( &status, &weightsNoise), &status);
    }

    /* normalize sfts */
    LAL_CALL( LALNormalizeSFTVect (&status, inputSFTs, uvar_blocksRngMed), &status);

  } /* end of sft reading block */
  
    
  /* generating peakgrams  */  
  pgV.length = mObsCoh;
  pgV.pg = NULL;
  pgV.pg = (HOUGHPeakGram *)LALCalloc(mObsCoh, sizeof(HOUGHPeakGram));

  { 
    SFTtype  sft;
    UCHARPeakGram     pg1;
    INT4   nPeaks, length;
    UINT4  j; 
  
    length = inputSFTs->data->data->length;
    pg1.length = length;
    pg1.data = NULL;
    pg1.data = (UCHAR *)LALCalloc(length, sizeof(UCHAR));

    /* loop over sfts and select peaks */
    for (j=0; j < mObsCoh; j++){

      sft = inputSFTs->data[j];

      LAL_CALL (SFTtoUCHARPeakGram( &status, &pg1, &sft, uvar_peakThreshold), &status);
      
      nPeaks = pg1.nPeaks;

      /* compress peakgram */      
      pgV.pg[j].length = nPeaks;
      pgV.pg[j].peak = NULL;
      pgV.pg[j].peak = (INT4 *)LALCalloc(nPeaks, sizeof(INT4));

      LAL_CALL( LALUCHAR2HOUGHPeak( &status, &(pgV.pg[j]), &pg1), &status );
    } /* end loop over sfts */

    /* we are done with the sfts and ucharpeakgram now */
    LAL_CALL (LALDestroySFTVector(&status, &inputSFTs), &status );
    LALFree(pg1.data);

  }/* end block for selecting peaks */



  /* compute detector velocity for SFT timestamps  */
  velV.length = mObsCoh;
  velV.data = NULL;
  velV.data = (REAL8Cart3Coor *)LALCalloc(mObsCoh, sizeof(REAL8Cart3Coor));

  {  
    VelocityPar   velPar;
    REAL8     vel[3];
    UINT4     j;

    velPar.detector = *detector;
    velPar.tBase = timeBase;
    velPar.vTol = ACCURACY; /* irrelevant */
    velPar.edat = NULL; 

    /*  ephemeris info */
    edat = (EphemerisData *)LALCalloc(1, sizeof(EphemerisData));
   (*edat).ephiles.earthEphemeris = uvar_earthEphemeris;
   (*edat).ephiles.sunEphemeris = uvar_sunEphemeris;

    /* read in ephemeris data */
    LAL_CALL( LALInitBarycenter( &status, edat), &status);
    velPar.edat = edat;

    /* calculate average velocity for each SFT duration */    
    for(j=0; j< mObsCoh; ++j){
      velPar.startTime.gpsSeconds     = timeV->data[j].gpsSeconds;
      velPar.startTime.gpsNanoSeconds = timeV->data[j].gpsNanoSeconds;
      
      LAL_CALL( LALAvgDetectorVel ( &status, vel, &velPar), &status ); 
      
      velV.data[j].x= vel[0];
      velV.data[j].y= vel[1];
      velV.data[j].z= vel[2];  
    }
  } /* end velocity calculation block */


  /* compute the time difference relative to startTime for all SFT */
  timeDiffV.length = mObsCoh;
  timeDiffV.data = NULL; 
  timeDiffV.data = (REAL8 *)LALCalloc(mObsCoh, sizeof(REAL8));

  {
    UINT4   j; 

    for(j=0; j < mObsCoh; ++j)
      timeDiffV.data[j] = XLALGPSDiff( timeV->data + j, timeV->data ) + 0.5*timeBase;
  }


  /* if we want to print expected sigma for each skypatch */
  if ( uvar_printSigma ) 
    {
      strcpy ( fileSigma, uvar_dirnameOut);
      strcat ( fileSigma, "/");
      strcat ( fileSigma, uvar_fbasenameOut);
      strcat ( fileSigma, "sigma_values");
      
      if ( (fpSigma = fopen(fileSigma,"w")) == NULL)
	{
	  fprintf(stderr,"Unable to find file %s for writing\n", fileSigma);
	  return DRIVEHOUGHCOLOR_EFILE;
	}
    }
  
  
  /* loop over sky patches */
  for (skyCounter = 0; skyCounter < nSkyPatches; skyCounter++)
    {
      UINT4 k;
      REAL8 sumWeightSquare;
      REAL8  alphaPeak, meanN, sigmaN, erfcInv;

      /* set sky positions and skypatch sizes */
      alpha = skyAlpha[skyCounter];
      delta = skyDelta[skyCounter];
      patchSizeX = skySizeDelta[skyCounter];
      patchSizeY = skySizeAlpha[skyCounter];

      fprintf(stdout, "starting skypatch %d of %d\n", skyCounter, nSkyPatches);
      
      /* calculate amplitude modulation weights */
      if (uvar_weighAM) {

	memcpy(weightsV.data, weightsNoise.data, mObsCoh * sizeof(REAL8));
	
	/*for ( k = 0; k < mObsCoh; k++)
	  weightsV.data[k] = weightsNoise.data[k]; */
	
	LAL_CALL( LALHOUGHComputeAMWeights( &status, &weightsV, timeV, detector, 
					    edat, alpha, delta), &status);
	LAL_CALL( LALHOUGHNormalizeWeights( &status, &weightsV), &status);
      }

      /* calculate the sum of the weights squared */
      sumWeightSquare = 0.0;
      for ( k = 0; k < mObsCoh; k++)
	sumWeightSquare += weightsV.data[k] * weightsV.data[k];


      /* computing the Hough threshold for a given false alarm  */
      /* HoughThreshold = N*alpha +sqrt(2 ||w||^2 * alpha *(1-alpha))*erfcinv(2 alpha_h) */      

      /* probability of selecting a peak expected mean and standard deviation for noise only */
      alphaPeak = exp( - uvar_peakThreshold);
      meanN = mObsCoh* alphaPeak; 
      sigmaN = sqrt(sumWeightSquare * alphaPeak * (1.0 - alphaPeak));

      if ( uvar_printSigma )
	fprintf(fpSigma, "%f \n", sigmaN);

      /* this should be  erfcInv =erfcinv(2.0 *uvar_houghFalseAlarm) */
      /* the function used is basically the inverse of the CDF for the 
	 Gaussian distribution with unit variance and
	 erfcinv(x) = gsl_cdf_ugaussian_Qinv (0.5*x)/sqrt(2) */
      /* First check that false alarm is within bounds 
	 and set it to something reasonable if not */
      if ( uvar_printEvents ) {
	erfcInv = gsl_cdf_ugaussian_Qinv (uvar_houghFalseAlarm)/sqrt(2);    
	houghThreshold = meanN + sigmaN*sqrt(2.0)*erfcInv;    
      }
      

      /* opening the output statistics and event files */

      /*  create directory fnameout/skypatch_$j using mkdir if required */
      if ( uvar_printStats || uvar_printEvents || uvar_printTemplates || uvar_printMaps )
	{

	  /* create the directory name uvar_dirnameOut/skypatch_$j */
	  strcpy(  filestats, uvar_dirnameOut);
	  strcat( filestats, "/skypatch_");
	  {
	    CHAR tempstr[16];
	    sprintf(tempstr, "%d", skyCounter+1);
	    strcat( filestats, tempstr);
	  }
	  strcat( filestats, "/");
	  
	  errno = 0;
	  {
	    /* check whether file can be created or if it exists already 
	       if not then exit */
	    INT4 mkdir_result;
	    mkdir_result = mkdir(filestats, S_IRWXU | S_IRWXG | S_IRWXO);
	    if ( (mkdir_result == -1) && (errno != EEXIST) )
	      {
		fprintf(stderr, "unable to create skypatch directory %d\n", skyCounter);
		return 1;  /* stop the program */
	      }
	  } /* dir created */
            
	  /* create the base filenames for the stats, histo and event files and template files*/
	  strcat( filestats, uvar_fbasenameOut);
	  strcpy( filehisto, filestats);

	} /* if ( uvar_printStats || uvar_printEvents || uvar_printTemplates || uvar_printMaps ) */


      if ( uvar_printEvents )
	strcpy( fileEvents, filestats);

      if ( uvar_printTemplates )
	strcpy(fileTemplates, filestats);

      if ( uvar_printMaps )
	strcpy(fileMaps, filestats);

      /* create and open the stats file for writing */
      if ( uvar_printStats )
	{
	  strcat(  filestats, "stats");
	  if ( (fpStats = fopen(filestats,"w")) == NULL)
	    {
	      fprintf(stderr,"Unable to find file %s for writing\n", filestats);
	      return DRIVEHOUGHCOLOR_EFILE;
	    }
	  /*setlinebuf(fpStats);*/ /*line buffered on */  
	  setvbuf(fpStats, (char *)NULL, _IOLBF, 0);      
	}

      if ( uvar_printEvents )
	{
	  /* create and open the events list file */
	  strcat(  fileEvents, "events");
	  if ( (fpEvents=fopen(fileEvents,"w")) == NULL ) 
	    {
	      fprintf(stderr,"Unable to find file %s\n", fileEvents);
	      return DRIVEHOUGHCOLOR_EFILE;
	    }
	  setvbuf(fpEvents, (char *)NULL, _IOLBF, 0);      
	  /*setlinebuf(fpEvents);*/ /*line buffered on */  
	}

      if ( uvar_printTemplates )
	{
	  /* create and open templates file */
	  strcat( fileTemplates, "templates");
	  if ( (fpTemplates = fopen(fileTemplates, "w")) == NULL)
	    {
	      fprintf(stderr, "Unable to create file %s\n", fileTemplates);
	      return DRIVEHOUGHCOLOR_EFILE;
	    }
	  setvbuf(fpTemplates, (char *)NULL, _IOLBF, 0);      
	  /*setlinebuf(fpTemplates);*/ /*line buffered on */   
	}


      /****  general parameter settings and 1st memory allocation ****/      
      lutV.length    = mObsCoh;
      lutV.lut = NULL;
      lutV.lut = (HOUGHptfLUT *)LALCalloc(mObsCoh, sizeof(HOUGHptfLUT));
      
      phmdVS.length  = mObsCoh;
      phmdVS.nfSize  = uvar_nfSizeCylinder;
      phmdVS.deltaF  = deltaF;
      phmdVS.phmd = NULL;
      phmdVS.phmd=(HOUGHphmd *)LALCalloc(mObsCoh*uvar_nfSizeCylinder, sizeof(HOUGHphmd));
      
      freqInd.deltaF = deltaF;
      freqInd.length = mObsCoh;
      freqInd.data = NULL;
      freqInd.data =  ( UINT8 *)LALCalloc(mObsCoh, sizeof(UINT8));
      

      /* for non-demodulated data (SFT input)*/
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
      parRes.pixelFactor = PIXELFACTOR;
      parRes.pixErr = PIXERR;
      parRes.linErr = LINERR;
      parRes.vTotC = VTOT;

      /* allocating histogram of the number-counts in the Hough maps */
      if ( uvar_printStats ) {
	UINT4   j;
	
	hist.length = mObsCoh+1;
	histTotal.length = mObsCoh+1;
	hist.data = NULL;
	histTotal.data = NULL;
	hist.data = (UINT8 *)LALCalloc((mObsCoh+1), sizeof(UINT8));
	histTotal.data = (UINT8 *)LALCalloc((mObsCoh+1), sizeof(UINT8));
	for(j = 0; j < histTotal.length; ++j)
	  histTotal.data[j]=0;
      }

      
      fBin= f0Bin;
      iHmap = 0;
      
      /* ***** for spin-down case ****/
      nSpin1Max = floor(uvar_nfSizeCylinder/2.0);
      f1jump = 1.0 / tObs;
      

      /* start of main loop over search frequency bins */
      /********** starting the search from f0Bin to fLastBin.
		  Note one set LUT might not cover all the interval.
		  This is taken into account *******************/

      while( fBin <= fLastBin){
	INT8 fBinSearch, fBinSearchMax;
	UINT4 i,j; 
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
	patch.xCoor = (REAL8 *)LALCalloc(xSide, sizeof(REAL8));
	patch.yCoor = (REAL8 *)LALCalloc(ySide, sizeof(REAL8));
	LAL_CALL( LALHOUGHFillPatchGrid( &status, &patch, &parSize ), &status );
	
	/*************** other memory allocation and settings************ */
	for(j=0; j<lutV.length; ++j){
	  lutV.lut[j].maxNBins = maxNBins;
	  lutV.lut[j].maxNBorders = maxNBorders;
	  lutV.lut[j].border =
	    (HOUGHBorder *)LALCalloc(maxNBorders, sizeof(HOUGHBorder));
	  lutV.lut[j].bin =
	    (HOUGHBin2Border *)LALCalloc(maxNBins, sizeof(HOUGHBin2Border));
	  for (i=0; i<maxNBorders; ++i){
	    lutV.lut[j].border[i].ySide = ySide;
	    lutV.lut[j].border[i].xPixel =
	      (COORType *)LALCalloc(ySide, sizeof(COORType));
	  }
	}
	for(j=0; j<phmdVS.length * phmdVS.nfSize; ++j){
	  phmdVS.phmd[j].maxNBorders = maxNBorders;
	  phmdVS.phmd[j].leftBorderP =
	    (HOUGHBorder **)LALCalloc(maxNBorders, sizeof(HOUGHBorder *));
	  phmdVS.phmd[j].rightBorderP =
	    (HOUGHBorder **)LALCalloc(maxNBorders, sizeof(HOUGHBorder *));
	  phmdVS.phmd[j].ySide = ySide;
	  phmdVS.phmd[j].firstColumn = NULL;
	  phmdVS.phmd[j].firstColumn = (UCHAR *)LALCalloc(ySide, sizeof(UCHAR));
	}
	
	/* ************* create all the LUTs at fBin ********************  */  
	for (j=0;j< mObsCoh;++j){  /* create all the LUTs */
	  parDem.veloC.x = velV.data[j].x;
	  parDem.veloC.y = velV.data[j].y;
	  parDem.veloC.z = velV.data[j].z;      
	  /* calculate parameters needed for buiding the LUT */
	  LAL_CALL( LALNDHOUGHParamPLUT( &status, &parLut, &parSize, &parDem),&status );
	  /* build the LUT */
	  LAL_CALL( LALHOUGHConstructPLUT( &status, &(lutV.lut[j]), &patch, &parLut ),
	       &status );
	}
        
	/************* build the set of  PHMD centered around fBin***********/     
	phmdVS.fBinMin = fBin - floor( uvar_nfSizeCylinder/2. );
	LAL_CALL( LALHOUGHConstructSpacePHMD(&status, &phmdVS, &pgV, &lutV), &status );
	if (uvar_weighAM || uvar_weighNoise) {
	  LAL_CALL( LALHOUGHWeighSpacePHMD(&status, &phmdVS, &weightsV), &status);
	}
	
	/* ************ initializing the Total Hough map space *********** */   
	ht.xSide = xSide;
	ht.ySide = ySide;
	ht.mObsCoh = mObsCoh;
	ht.deltaF = deltaF;
	ht.map   = NULL;
	ht.map   = (HoughTT *)LALCalloc(xSide*ySide, sizeof(HoughTT));
	LAL_CALL( LALHOUGHInitializeHT( &status, &ht, &patch), &status); /*not needed */
	

	/*  Search frequency interval possible using the same LUTs */
	fBinSearch = fBin;
	fBinSearchMax= fBin + parSize.nFreqValid - 1 - floor( (uvar_nfSizeCylinder - 1)/2.);
	

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
	    
	    for( n=0; n<= nSpin1Max; ++n){ 
	      /*loop over all spindown values */

	      f1dis = - n*f1jump;
	      ht.spinRes.data[0] =  f1dis*deltaF;
	      
	      for (j=0;j< mObsCoh;++j){
		freqInd.data[j] = fBinSearch + floor(timeDiffV.data[j]*f1dis+0.5);
	      }
	      
	      if (uvar_weighAM || uvar_weighNoise) {
		LAL_CALL( LALHOUGHConstructHMT_W( &status, &ht, &freqInd, &phmdVS ), &status );
	      }
	      else {
		LAL_CALL( LALHOUGHConstructHMT( &status, &ht, &freqInd, &phmdVS ), &status );	  
	      }

	      /* ********************* perfom stat. analysis on the maps ****************** */
	      LAL_CALL( LALHoughStatistics ( &status, &stats, &ht), &status );
	      LAL_CALL( LALStereo2SkyLocation (&status, &sourceLocation, 
				       stats.maxIndex[0], stats.maxIndex[1], &patch, &parDem), &status);

	      /* calculate histogram if necessary */
	      if ( uvar_printStats ) {
		LAL_CALL( LALHoughHistogram ( &status, &hist, &ht), &status);

		for(j=0; j< histTotal.length; ++j){ 
		  histTotal.data[j]+=hist.data[j]; 
		}	      
	      }

	      

	      significance =  (stats.maxCount - meanN)/sigmaN;	      
	      if ( significance > nStarEventVec.event[fBinSearch-f0Bin].nStarSignificance )
		{
		  nStarEventVec.event[fBinSearch-f0Bin].nStar = stats.maxCount;
		  nStarEventVec.event[fBinSearch-f0Bin].nStarSignificance = significance;
		  nStarEventVec.event[fBinSearch-f0Bin].freqStar = fBinSearch * deltaF;
		  nStarEventVec.event[fBinSearch-f0Bin].alphaStar = sourceLocation.alpha;
		  nStarEventVec.event[fBinSearch-f0Bin].deltaStar = sourceLocation.delta;
		  nStarEventVec.event[fBinSearch-f0Bin].fdotStar = ht.spinRes.data[0];
		}

	      /* ***** print results *********************** */

	      if( uvar_printMaps )
		  if( PrintHmap2m_file( &ht, fileMaps, iHmap ) ) return 5;

	      if ( uvar_printStats )
		fprintf(fpStats, "%d %f %f %f %f %f %f %f %g\n",
			iHmap, sourceLocation.alpha, sourceLocation.delta,
			(REAL4)stats.maxCount, (REAL4)stats.minCount, stats.avgCount,stats.stdDev,
			(fBinSearch*deltaF), ht.spinRes.data[0]);

	      if ( uvar_printEvents )
		LAL_CALL( PrintHoughEvents (&status, fpEvents, houghThreshold, &ht,
					    &patch, &parDem), &status );

	      if ( uvar_printTemplates )
		LAL_CALL( PrintHoughEvents (&status, fpTemplates, 0.0, &ht, &patch, &parDem), &status);

	      ++iHmap;
	    } /* end loop over spindown values */ 
	    
	    LALFree(ht.spinRes.data);
	    
	    
	    /***** shift the search freq. & PHMD structure 1 freq.bin ****** */
	    ++fBinSearch;
	    
	    LAL_CALL( LALHOUGHupdateSpacePHMDup(&status, &phmdVS, &pgV, &lutV), &status );
	    
	    if (uvar_weighAM || uvar_weighNoise) {
	      LAL_CALL( LALHOUGHWeighSpacePHMD( &status, &phmdVS, &weightsV), &status);	    
	    }

	  }   /* ********>>>>>>  closing second while  <<<<<<<<**********<  */
	
	fBin = fBinSearch;
	
	/* ********************  Free partial memory ******************* */
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
	
      } /* closing while */
      
      /******************************************************************/
      /* printing total histogram */
      /******************************************************************/
      if ( uvar_printStats )
	if( PrintHistogram( &histTotal, filehisto) ) return 7;
      
      /******************************************************************/
      /* closing files with statistics results and events */
      /******************************************************************/  
      if ( uvar_printStats )
	fclose(fpStats);

      if ( uvar_printEvents )
	fclose(fpEvents);

      if ( uvar_printTemplates )
	fclose(fpTemplates);
     

      /******************************************************************/
      /* Free memory allocated inside skypatches loop */
      /******************************************************************/
      
      LALFree(lutV.lut);  
      
      LALFree(phmdVS.phmd);
      LALFree(freqInd.data);
      if ( uvar_printStats ) {
	LALFree(hist.data);
	LALFree(histTotal.data);
      }
      
#ifdef TIMING
      stop = realcc();
      printf(" All: %llu\n", stop-start);
#endif
  
    } /* finish loop over skypatches */

  /* close sigma file */
  if ( uvar_printSigma )
    fclose(fpSigma);

  /* print most significant events */
  LAL_CALL( PrintnStarFile( &status, &nStarEventVec, uvar_dirnameOut, 
			    uvar_fbasenameOut), &status);


  /* free memory allocated outside skypatches loop */ 
  LALFree( detector);

  {
    UINT4 j;
    for (j = 0; j < mObsCoh; ++j) LALFree( pgV.pg[j].peak); 
  }
  LALFree(pgV.pg);

  LAL_CALL(LALDestroyTimestampVector ( &status, &timeV), &status); 
  
  LALFree(timeDiffV.data);

  LALFree(velV.data);

  LALFree(weightsV.data);
  LALFree(weightsNoise.data);  

  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);

  LALFree(skyAlpha);
  LALFree(skyDelta);
  LALFree(skySizeAlpha);
  LALFree(skySizeDelta);

  LALFree( nStarEventVec.event );

  LAL_CALL (LALDestroyUserVars(&status), &status);

  LALCheckMemoryLeaks();

  if ( lalDebugLevel )
    REPORTSTATUS ( &status);

  return status.statusCode;
}



  
/******************************************************************/
/* printing the Histogram of all maps into a file                    */
/******************************************************************/
  
int PrintHistogram(UINT8Vector *hist, CHAR *fnameOut){

  FILE  *fp=NULL;   /* Output file */
  char filename[256];
  UINT4  i ;
 
  strcpy(  filename, fnameOut);
  strcat(  filename, "histo");
  
  if ( (fp = fopen(filename,"w")) == NULL)
    {  
      fprintf(stderr,"Unable to find file %s\n",filename);
      return DRIVEHOUGHCOLOR_EFILE; 
    }

  for (i=0; i < hist->length; i++){
    fprintf(fp,"%d  %" LAL_UINT8_FORMAT "\n", i, hist->data[i]);
  }
  
  fclose( fp );  
  return 0;
}



/******************************************************************/
/* printing the HM into a file                    */
/******************************************************************/

int PrintHmap2file(HOUGHMapTotal *ht, CHAR *fnameOut, INT4 iHmap){

  FILE  *fp=NULL;   /* Output file */
  char filename[256], filenumber[16]; 
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
  char filename[256], filenumber[16]; 
  INT4  k, i ;
  UINT2 xSide, ySide;
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


/******************************************************************/
/*  Find and print events to a given open file */
/******************************************************************/
void PrintHoughEvents (LALStatus       *status,
        	      FILE            *fpEvents,
	  	      INT4            houghThreshold,
		      HOUGHMapTotal   *ht,
	    	      HOUGHPatchGrid  *patch,
	 	      HOUGHDemodPar   *parDem)
{

  REAL8UnitPolarCoor sourceLocation;
  UINT2    xPos, yPos, xSide, ySide;
  REAL8    temp;
  REAL8    f0;
  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  
 /* make sure arguments are not null */
  ASSERT (patch , status, DRIVEHOUGHCOLOR_ENULL,DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (parDem, status, DRIVEHOUGHCOLOR_ENULL,DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (ht, status, DRIVEHOUGHCOLOR_ENULL,DRIVEHOUGHCOLOR_MSGENULL);

 /* make sure input hough map is ok*/
  ASSERT (ht->xSide > 0, status, DRIVEHOUGHCOLOR_EBAD,DRIVEHOUGHCOLOR_MSGEBAD);
  ASSERT (ht->ySide > 0, status, DRIVEHOUGHCOLOR_EBAD,DRIVEHOUGHCOLOR_MSGEBAD);
  
  /* read input parameters */
  xSide = ht->xSide;
  ySide = ht->ySide; 
  
  f0=(ht->f0Bin)*(ht->deltaF);
  
  for(yPos =0; yPos<ySide; yPos++){
    for(xPos =0; xPos<xSide; xPos++){
      /* read the current number count */
      temp = ht->map[yPos*xSide + xPos];
      if(temp > houghThreshold){
        TRY( LALStereo2SkyLocation(status->statusPtr, 
				&sourceLocation,xPos,yPos,patch, parDem), status);
	if (ht->spinRes.length) {
	  fprintf(fpEvents, "%g %g %g %g %g \n", 
		  temp, sourceLocation.alpha, sourceLocation.delta, 
		  f0, ht->spinRes.data[0]);
	}
	else {
	  fprintf(fpEvents, "%g %g %g %g %g \n", 
		  temp, sourceLocation.alpha, sourceLocation.delta, 
		  f0,0.00);
	}
      }      
    }
  }
  	 
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}    
/* >>>>>>>>>>>>>>>>>>>>>*************************<<<<<<<<<<<<<<<<<<<< */





void PrintLogFile (LALStatus       *status,
		   CHAR            *dir,
		   CHAR            *basename,
		   CHAR            *skyfile,
		   CHAR            *linefile,
		   CHAR            *executable )
{
  CHAR *fnameLog=NULL; 
  FILE *fpLog=NULL;
  CHAR *logstr=NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  
  /* open log file for writing */
  fnameLog = (CHAR *)LALCalloc( 512, sizeof(CHAR));
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
  fprintf( fpLog, "%s", logstr);
  LALFree(logstr);

  /* copy contents of skypatch file into logfile */
  fprintf(fpLog, "\n\n# Contents of skypatch file:\n");
  fclose(fpLog);
  {
    CHAR command[1024] = "";
    sprintf(command, "cat %s >> %s", skyfile, fnameLog);
    if ( system(command) ) fprintf (stderr, "\nsystem('%s') returned non-zero status!\n\n", command );

  }

  /* copy contents of linefile if necessary */
  if ( linefile ) {
    if ((fpLog = fopen(fnameLog, "a")) != NULL) 
      {
	CHAR command[1024] = "";
	fprintf (fpLog, "\n\n# Contents of linefile:\n");
	fprintf (fpLog, "# -----------------------------------------\n");
	fclose (fpLog);
	sprintf(command, "cat %s >> %s", linefile, fnameLog);
        if ( system(command) ) fprintf (stderr, "\nsystem('%s') returned non-zero status!\n\n", command );

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
      /* we don't check this. If it fails, we assume that */
      /* one of the system-commands was not available, and */
      /* therefore the CVS-versions will not be logged */
      if ( system(command) ) fprintf (stderr, "\nsystem('%s') returned non-zero status!\n\n", command );
    }

  LALFree(fnameLog); 
  	 
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}    


/* print most significant events */
void PrintnStarFile (LALStatus                   *status,
		     HoughSignificantEventVector *eventVec,
		     CHAR                        *dirname,
		     CHAR                        *basename)
{
  CHAR *filestar = NULL; 
  FILE *fpStar = NULL; 
  INT4 length, starIndex;
  HoughSignificantEvent *event;
  INT4 mkdir_result;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT(eventVec, status, DRIVEHOUGHCOLOR_ENULL,DRIVEHOUGHCOLOR_MSGENULL); 
  ASSERT(eventVec->event, status, DRIVEHOUGHCOLOR_ENULL,DRIVEHOUGHCOLOR_MSGENULL); 
  ASSERT(dirname, status, DRIVEHOUGHCOLOR_ENULL,DRIVEHOUGHCOLOR_MSGENULL); 
  ASSERT(basename, status, DRIVEHOUGHCOLOR_ENULL,DRIVEHOUGHCOLOR_MSGENULL); 

  /* create the directory for writing nstar */
  filestar = (CHAR *)LALCalloc( 512, sizeof(CHAR));
  strcpy( filestar, dirname);
  strcat( filestar, "/nstarfiles/");
  errno = 0;
  /* check whether file can be created or if it exists already 
     if not then exit */
  mkdir_result = mkdir(filestar, S_IRWXU | S_IRWXG | S_IRWXO);
  if ( (mkdir_result == -1) && (errno != EEXIST) ) {
    ABORT( status, DRIVEHOUGHCOLOR_EDIR, DRIVEHOUGHCOLOR_MSGEDIR);
  }
  
  strcat( filestar, basename );
  strcat( filestar, "nstar");

  /* open the nstar file for writing */
  if ( (fpStar = fopen(filestar,"w")) == NULL) {
    ABORT( status, DRIVEHOUGHCOLOR_EFILE, DRIVEHOUGHCOLOR_MSGEFILE);
  }
  
  /*line buffering */  
  setvbuf(fpStar, (char *)NULL, _IOLBF, 0);      

  /* write the nstar results */
  length = eventVec->length;
  event = eventVec->event;
  for(starIndex = 0; starIndex < length; starIndex++)
    {
      fprintf(fpStar, "%f %f %f %f %f %g \n", event->nStar, event->nStarSignificance, event->freqStar, 
	      event->alphaStar, event->deltaStar, event->fdotStar );
      event++;
    }

  /* close nstar file */
  fclose(fpStar);
  LALFree(filestar); 
  	 
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}    
