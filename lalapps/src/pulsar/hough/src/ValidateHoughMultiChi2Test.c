/*
 *  Copyright (C) 2006 Llucia Sancho de la Jordana,
 *  Badri Krishnan, Alicia M. Sintes  
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
 * data using SFTs from possible multiple IFOs
 *
 * History: Created by Sancho de la Jordana, Sintes and Krishnan December 15, 2006
 */


#include "./DriveHoughColor.h"
#include "./MCInjectHoughMulti.h"

/* globals, constants and defaults */



/* boolean global variables for controlling output */
BOOLEAN uvar_printEvents, uvar_printTemplates, uvar_printMaps, uvar_printStats, uvar_printSigma;

/* #define EARTHEPHEMERIS "./earth05-09.dat" */
/* #define SUNEPHEMERIS "./sun05-09.dat"    */

#define EARTHEPHEMERIS "/home/llucia/chi2/earth05-09.dat"
#define SUNEPHEMERIS "/home/llucia/chi2/sun05-09.dat"

#define MAXFILENAMELENGTH 512 /* maximum # of characters  of a filename */

#define DIROUT "./outMultiChi2Test"   /* output directory */
#define BASENAMEOUT "HM"    /* prefix file output */

#define THRESHOLD 1.6 /* thresold for peak selection, with respect to the
                              the averaged power in the search band */
#define FALSEALARM 1.0e-9 /* Hough false alarm for candidate selection */
#define SKYFILE "./skypatchfile"      
#define F0 310.0   /*  frequency to build the LUT and start search */
#define FBAND 0.05   /* search frequency band  */
#define NFSIZE  21   /* n-freq. span of the cylinder, to account for spin-down search */
#define BLOCKSRNGMED 101 /* Running median window size */

#define TRUE (1==1)
#define FALSE (1==0)

#define LAL_INT4_MAX 2147483647

#define NBLOCKSTEST 8 /* number of data blocks to do Chi2 test */


#define SFTDIRECTORY "/local_data/sintes/SFT-S5-120-130/*SFT*.*"

/* local function prototype */



/* ****************************************
 * Structure, HoughParamsTest, typedef
 */

typedef struct tagHoughParamsTest{
    UINT4  length;            /* number p of blocks to split the data into */
    UINT4  *numberSFTp;       /* Ni SFTs in each data block */
    REAL8  *sumWeight;        /* Pointer to the sumWeight of each block of data */
    REAL8  *sumWeightSquare;  /* Pointer to the sumWeightSquare of each block of data */
}HoughParamsTest; 

/******************************************/

/* function to Split the SFTs into p blocks */

void SplitSFTs(LALStatus         *status, 
	       REAL8Vector       *weightsV, 
	       HoughParamsTest   *chi2Params);


/******************************************/

int main(int argc, char *argv[]){

  /* LALStatus pointer */
  static LALStatus  status;  
  
  /* time and velocity  */
  static LIGOTimeGPSVector    timeV;
  static REAL8Cart3CoorVector velV;
  static REAL8Vector          timeDiffV;
  LIGOTimeGPS firstTimeStamp;

  /* standard pulsar sft types */ 
  MultiSFTVector *inputSFTs = NULL;
  UINT4 binsSFT;
  UINT4 sftFminBin;
  UINT4 numsft;

  INT4 k;
  FILE *fp=NULL;

  /* information about all the ifos */
  MultiDetectorStateSeries *mdetStates = NULL;
  UINT4 numifo;

  /* vector of weights */
  REAL8Vector weightsV;

  /* ephemeris */
  EphemerisData    *edat=NULL;

  static UCHARPeakGram       pg1;
  static HoughTemplate  pulsarTemplate;
  static REAL8Vector  foft; 

  /* miscellaneous */
  UINT4  mObsCoh;
  REAL8  timeBase, deltaF;
  REAL8  numberCount;

  /* Chi2Test parameters */
  HoughParamsTest chi2Params;
  REAL8Vector numberCountV;  /* Vector with the number count of each block inside */
  REAL8  numberCountTotal;   /* Sum over all the numberCounts */
  REAL8  chi2;

  /* sft constraint variables */
  LIGOTimeGPS startTimeGPS, endTimeGPS;
  LIGOTimeGPSVector *inputTimeStampsVector=NULL;


  REAL8  alphaPeak, meanN, sigmaN;

  /* user input variables */
  BOOLEAN  uvar_help, uvar_weighAM, uvar_weighNoise;
  INT4     uvar_blocksRngMed, uvar_nfSizeCylinder, uvar_maxBinsClean;
  REAL8    uvar_startTime, uvar_endTime;
  REAL8    uvar_fStart, uvar_peakThreshold, uvar_fSearchBand;
  REAL8    uvar_Alpha, uvar_Delta, uvar_Freq, uvar_fdot;
  REAL8    uvar_AlphaWeight, uvar_DeltaWeight;
  CHAR     *uvar_earthEphemeris=NULL;
  CHAR     *uvar_sunEphemeris=NULL;
  CHAR     *uvar_sftDir=NULL;
  CHAR     *uvar_timeStampsFile=NULL;
  CHAR     *uvar_outfile=NULL;
  LALStringVector *uvar_linefiles=NULL;
  INT4     uvar_p;


  /* Set up the default parameters */  

  /* LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;
    
  uvar_help = FALSE;
  uvar_weighAM = TRUE;
  uvar_weighNoise = TRUE;
  uvar_blocksRngMed = BLOCKSRNGMED;
  uvar_nfSizeCylinder = NFSIZE;
  uvar_fStart = F0;
  uvar_fSearchBand = FBAND;
  uvar_peakThreshold = THRESHOLD;
  uvar_maxBinsClean = 100;
  uvar_startTime= 0;
  uvar_endTime = LAL_INT4_MAX;

  uvar_Alpha = 1.0;
  uvar_Delta = 1.0;
  uvar_Freq = 310.0;
  uvar_fdot = 0.0;

  uvar_AlphaWeight = uvar_Alpha;
  uvar_DeltaWeight = uvar_Delta;

  uvar_p = NBLOCKSTEST;
  chi2Params.length=uvar_p;
  chi2Params.numberSFTp=NULL;
  chi2Params.sumWeight=NULL;
  chi2Params.sumWeightSquare=NULL;

  uvar_outfile = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_outfile, "./tempout");

  uvar_earthEphemeris = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_earthEphemeris,EARTHEPHEMERIS);

  uvar_sunEphemeris = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_sunEphemeris,SUNEPHEMERIS);


  uvar_sftDir = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_sftDir,SFTDIRECTORY);


  /* register user input variables */
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "help",            'h', UVAR_HELP,     "Print this message",                    &uvar_help),            &status);  
  LAL_CALL( LALRegisterREALUserVar(   &status, "fStart",          'f', UVAR_OPTIONAL, "Start search frequency",                &uvar_fStart),              &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fSearchBand",     'b', UVAR_OPTIONAL, "Search frequency band",                 &uvar_fSearchBand),     &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "startTime",        0,  UVAR_OPTIONAL, "GPS start time of observation",         &uvar_startTime),        &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "endTime",          0,  UVAR_OPTIONAL, "GPS end time of observation",           &uvar_endTime),          &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "timeStampsFile",   0,  UVAR_OPTIONAL, "Input time-stamps file",                &uvar_timeStampsFile),   &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "peakThreshold",    0,  UVAR_OPTIONAL, "Peak selection threshold",              &uvar_peakThreshold),   &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "weighAM",          0,  UVAR_OPTIONAL, "Use amplitude modulation weights",      &uvar_weighAM),         &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "weighNoise",       0,  UVAR_OPTIONAL, "Use SFT noise weights",                 &uvar_weighNoise),      &status);  
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "earthEphemeris",  'E', UVAR_OPTIONAL, "Earth Ephemeris file",                  &uvar_earthEphemeris),  &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sunEphemeris",    'S', UVAR_OPTIONAL, "Sun Ephemeris file",                    &uvar_sunEphemeris),    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftDir",          'D', UVAR_REQUIRED, "SFT filename pattern",                  &uvar_sftDir),          &status);
  LAL_CALL( LALRegisterLISTUserVar(   &status, "linefiles",        0,  UVAR_OPTIONAL, "Comma separated List of linefiles (filenames must contain IFO name)",  
	      &uvar_linefiles),       &status);

  LAL_CALL( LALRegisterREALUserVar(   &status, "Alpha",            0,  UVAR_OPTIONAL, "Sky location (longitude)",               &uvar_Alpha),              &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "Delta",            0,  UVAR_OPTIONAL, "Sky location (latitude)",                &uvar_Delta),              &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "Freq",             0,  UVAR_OPTIONAL, "Template frequency",                     &uvar_Freq),               &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fdot",             0,  UVAR_OPTIONAL, "First spindown",                         &uvar_fdot),               &status);

  LAL_CALL( LALRegisterREALUserVar(   &status, "AlphaWeight",      0,  UVAR_OPTIONAL, "sky Alpha for weight calculation",       &uvar_AlphaWeight),        &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "DeltaWeight",      0,  UVAR_OPTIONAL, "sky Delta for weight calculation",       &uvar_DeltaWeight),        &status);

  LAL_CALL( LALRegisterINTUserVar(    &status, "nfSizeCylinder",   0,  UVAR_OPTIONAL, "Size of cylinder of PHMDs",             &uvar_nfSizeCylinder),  &status);

  LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed",     0,  UVAR_OPTIONAL, "Running Median block size",             &uvar_blocksRngMed),    &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "maxBinsClean",     0,  UVAR_OPTIONAL, "Maximum number of bins in cleaning",    &uvar_maxBinsClean),    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "outfile",          0,  UVAR_OPTIONAL, "output file name",                      &uvar_outfile),         &status);

  LAL_CALL( LALRegisterINTUserVar(    &status, "pdatablock",     'p',  UVAR_OPTIONAL, "Number of data blocks for veto tests",  &uvar_p),               &status);


  /* read all command line variables */
  LAL_CALL( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 

  /* very basic consistency checks on user input */
  if ( uvar_fStart < 0 ) {
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


  /***** start main calculations *****/

  chi2Params.length=uvar_p;
  chi2Params.numberSFTp = (UINT4 *)LALMalloc( uvar_p*sizeof(UINT4));
  chi2Params.sumWeight = (REAL8 *)LALMalloc( uvar_p*sizeof(REAL8));
  chi2Params.sumWeightSquare = (REAL8 *)LALMalloc( uvar_p*sizeof(REAL8));


  /* read sft Files and set up weights */
  {
    /* new SFT I/O data types */
    SFTCatalog *catalog = NULL;
    static SFTConstraints constraints;

    REAL8 doppWings, f_min, f_max;

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
      LAL_CALL ( LALReadTimestampsFile ( &status, &inputTimeStampsVector, uvar_timeStampsFile), &status);
      constraints.timestamps = inputTimeStampsVector;
    }
    
    /* get sft catalog */
    LAL_CALL( LALSFTdataFind( &status, &catalog, uvar_sftDir, &constraints), &status);
    if ( (catalog == NULL) || (catalog->length == 0) ) {
      fprintf (stderr,"Unable to match any SFTs with pattern '%s'\n", uvar_sftDir );
      exit(1);
    }

    /* now we can free the inputTimeStampsVector */
    if ( LALUserVarWasSet( &uvar_timeStampsFile ) ) {
      LALDestroyTimestampVector ( &status, &inputTimeStampsVector);
    }

    /* get some sft parameters */
    mObsCoh = catalog->length; /* number of sfts */
    deltaF = catalog->data->header.deltaF;  /* frequency resolution */
    timeBase= 1.0/deltaF; /* coherent integration time */
    // unused: UINT8 f0Bin = floor( uvar_fStart * timeBase + 0.5); /* initial search frequency */
    // unused: INT4 length =  uvar_fSearchBand * timeBase; /* total number of search bins - 1 */
    
    /* catalog is ordered in time so we can get start, end time and tObs*/
    firstTimeStamp = catalog->data[0].header.epoch;
    // unused: LIGOTimeGPS lastTimeStamp = catalog->data[mObsCoh - 1].header.epoch;

    /* allocate memory for velocity vector */
    velV.length = mObsCoh;
    velV.data = NULL;
    velV.data = (REAL8Cart3Coor *)LALCalloc(mObsCoh, sizeof(REAL8Cart3Coor));

    /* allocate memory for timestamps vector */
    timeV.length = mObsCoh;
    timeV.data = NULL;
    timeV.data = (LIGOTimeGPS *)LALCalloc( mObsCoh, sizeof(LIGOTimeGPS));

    /* allocate memory for vector of time differences from start */
    timeDiffV.length = mObsCoh;
    timeDiffV.data = NULL; 
    timeDiffV.data = (REAL8 *)LALCalloc(mObsCoh, sizeof(REAL8));
  
    /* add wings for Doppler modulation and running median block size*/
    doppWings = (uvar_fStart + uvar_fSearchBand) * VTOT;    
    f_min = uvar_fStart - doppWings - (uvar_blocksRngMed + uvar_nfSizeCylinder) * deltaF;
    f_max = uvar_fStart + uvar_fSearchBand + doppWings + (uvar_blocksRngMed + uvar_nfSizeCylinder) * deltaF;

    /* read sft files making sure to add extra bins for running median */
    /* read the sfts */
    LAL_CALL( LALLoadMultiSFTs ( &status, &inputSFTs, catalog, f_min, f_max), &status);


    /* clean sfts if required */
    if ( LALUserVarWasSet( &uvar_linefiles ) )
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

	LAL_CALL( LALRemoveKnownLinesInMultiSFTVector ( &status, inputSFTs, uvar_maxBinsClean, uvar_blocksRngMed, uvar_linefiles, randPar), &status);

	LAL_CALL ( LALDestroyRandomParams (&status, &randPar), &status);
	fclose(fpRand);
      } /* end cleaning */


    /* SFT info -- assume all SFTs have same length */
    numifo = inputSFTs->length;
    binsSFT = inputSFTs->data[0]->data->data->length;
    sftFminBin = (INT4) floor(inputSFTs->data[0]->data[0].f0 * timeBase + 0.5);    


    LAL_CALL( LALDestroySFTCatalog( &status, &catalog ), &status);  	

  } /* end of sft reading block */



  /* get detector velocities weights vector, and timestamps */
  { 
    MultiNoiseWeights *multweight = NULL;    
    MultiPSDVector *multPSD = NULL;  
    UINT4 iIFO, iSFT, j;

    /*  get ephemeris  */
    edat = (EphemerisData *)LALCalloc(1, sizeof(EphemerisData));
    (*edat).ephiles.earthEphemeris = uvar_earthEphemeris;
    (*edat).ephiles.sunEphemeris = uvar_sunEphemeris;
    LAL_CALL( LALInitBarycenter( &status, edat), &status);


    /* normalize sfts */
    LAL_CALL( LALNormalizeMultiSFTVect (&status, &multPSD, inputSFTs, uvar_blocksRngMed), &status);

    /* set up weights */
    weightsV.length = mObsCoh;
    weightsV.data = (REAL8 *)LALCalloc(1, mObsCoh * sizeof(REAL8));

    /* initialize all weights to unity */
    LAL_CALL( LALHOUGHInitializeWeights( &status, &weightsV), &status);
   
    /* compute multi noise weights if required */
    if ( uvar_weighNoise ) {
      LAL_CALL ( LALComputeMultiNoiseWeights ( &status, &multweight, multPSD, uvar_blocksRngMed, 0), &status);
    }
    
    /* we are now done with the psd */
    LAL_CALL ( LALDestroyMultiPSDVector  ( &status, &multPSD), &status);

    /* get information about all detectors including velocity and timestamps */
    /* note that this function returns the velocity at the 
       mid-time of the SFTs -- should not make any difference */
    LAL_CALL ( LALGetMultiDetectorStates ( &status, &mdetStates, inputSFTs, edat), &status);


    /* copy the timestamps, weights, and velocity vector */
    for (j = 0, iIFO = 0; iIFO < numifo; iIFO++ ) {

      numsft = mdetStates->data[iIFO]->length;
      
      for ( iSFT = 0; iSFT < numsft; iSFT++, j++) {

	velV.data[j].x = mdetStates->data[iIFO]->data[iSFT].vDetector[0];
	velV.data[j].y = mdetStates->data[iIFO]->data[iSFT].vDetector[1];
	velV.data[j].z = mdetStates->data[iIFO]->data[iSFT].vDetector[2];

	if ( uvar_weighNoise )
	  weightsV.data[j] = multweight->data[iIFO]->data[iSFT];

	/* mid time of sfts */
	timeV.data[j] = mdetStates->data[iIFO]->data[iSFT].tGPS;

      } /* loop over SFTs */

    } /* loop over IFOs */

    if ( uvar_weighNoise ) {
      LAL_CALL( LALHOUGHNormalizeWeights( &status, &weightsV), &status);
    }

    /* compute the time difference relative to startTime for all SFTs */
    for(j = 0; j < mObsCoh; j++)
      timeDiffV.data[j] = XLALGPSDiff( timeV.data + j, &firstTimeStamp );

    if ( uvar_weighNoise ) {    
      LAL_CALL ( LALDestroyMultiNoiseWeights ( &status, &multweight), &status);
    }

  } /* end block for noise weights, velocity and time */
  

  
  /* calculate amplitude modulation weights if required */
  if (uvar_weighAM) 
    {
      MultiAMCoeffs *multiAMcoef = NULL;
      UINT4 iIFO, iSFT;
      SkyPosition skypos;      

      /* get the amplitude modulation coefficients */
      skypos.longitude = uvar_AlphaWeight;
      skypos.latitude = uvar_DeltaWeight;
      skypos.system = COORDINATESYSTEM_EQUATORIAL;
      LAL_CALL ( LALGetMultiAMCoeffs ( &status, &multiAMcoef, mdetStates, skypos), &status);
      
      /* loop over the weights and multiply them by the appropriate
	 AM coefficients */
      for ( k = 0, iIFO = 0; iIFO < numifo; iIFO++) {
	
	numsft = mdetStates->data[iIFO]->length;
	
	for ( iSFT = 0; iSFT < numsft; iSFT++, k++) {	  
	  
	  REAL8 a, b;
	  
	  a = multiAMcoef->data[iIFO]->a->data[iSFT];
	  b = multiAMcoef->data[iIFO]->b->data[iSFT];    
	  weightsV.data[k] *= (a*a + b*b);
	} /* loop over SFTs */
      } /* loop over IFOs */
      
      LAL_CALL( LALHOUGHNormalizeWeights( &status, &weightsV), &status);
      
      XLALDestroyMultiAMCoeffs ( multiAMcoef );
    } /* end AM weights calculation */


  /* misc. memory allocations */

  /* memory for one spindown */  
  pulsarTemplate.spindown.length = 1; 
  pulsarTemplate.spindown.data = NULL; 
  pulsarTemplate.spindown.data = (REAL8 *)LALMalloc(sizeof(REAL8));

  /* copy template parameters */
  pulsarTemplate.spindown.data[0] = uvar_fdot;
  pulsarTemplate.f0 = uvar_Freq;
  pulsarTemplate.latitude = uvar_Delta;
  pulsarTemplate.longitude = uvar_Alpha;

  /* memory for f(t) vector */
  foft.length = mObsCoh;
  foft.data = NULL;
  foft.data = (REAL8 *)LALMalloc(mObsCoh*sizeof(REAL8));

  /* memory for peakgram */
  pg1.length = binsSFT;
  pg1.data = NULL;
  pg1.data = (UCHAR *)LALCalloc( binsSFT, sizeof(UCHAR));
 
  /* memory for number Count Vector */
  numberCountV.length = uvar_p;
  numberCountV.data = NULL;
  numberCountV.data = (REAL8 *)LALMalloc( uvar_p*sizeof(REAL8));

  /* block for calculating peakgram and number count */  
  {
    UINT4 iIFO, iSFT, ii, numberSFTp;
    INT4 ind;
    REAL8 sumWeightSquare;
    SFTtype  *sft;        
  
    /* compute mean and sigma for noise only */    
    /* first calculate the sum of the weights squared */
    sumWeightSquare = 0.0;
    for ( k = 0; k < (INT4)mObsCoh; k++)
      sumWeightSquare += weightsV.data[k] * weightsV.data[k];
    
    /* probability of selecting a peak expected mean and standard deviation for noise only */
    alphaPeak = exp( - uvar_peakThreshold);
    meanN = mObsCoh* alphaPeak; 
    sigmaN = sqrt(sumWeightSquare * alphaPeak * (1.0 - alphaPeak));
        
 
   /* the received frequency as a function of time  */
   LAL_CALL( ComputeFoft(&status, &foft, &pulsarTemplate, &timeDiffV, &velV, timeBase), &status);   
   
   LAL_CALL(SplitSFTs(&status, &weightsV, &chi2Params), &status);
   
   /* loop over SFT, generate peakgram and get number count */
   UINT4 j;
   j=0;
   iIFO=0;
   iSFT=0;
   
   numsft = mdetStates->data[iIFO]->length;
  
     
      for (k=0 ; k<uvar_p ; k++ ){
       
         numberSFTp=chi2Params.numberSFTp[k];
         numberCount = 0;

         for (ii=0 ; (ii < numberSFTp)&&(iIFO<numifo) ; ii++) {
	   
	    sft = inputSFTs->data[iIFO]->data + iSFT;

            LAL_CALL (SFTtoUCHARPeakGram( &status, &pg1, sft, uvar_peakThreshold), &status);	    

            ind = floor( foft.data[j]*timeBase - sftFminBin + 0.5); 
            
	    numberCount += pg1.data[ind]*weightsV.data[j];
	    
	    j++;

	    iSFT++;

	    if (iSFT >= numsft){
		
		iIFO++;
		iSFT=0;
		if (iIFO<numifo){
		numsft = mdetStates->data[iIFO]->length;
		}
	    }

         } /* loop over SFTs */
       
         numberCountV.data[k]=numberCount;
       	  
      } /* loop over blocks */
  
  }
  
/* Chi2 Test */
  {
      REAL8   eta;                /* Auxiliar variable */ 
      REAL8   nj, sumWeightj, sumWeightSquarej;
                     

      numberCountTotal=0;
      chi2=0;
 
      for(k=0; k<uvar_p ; k++){
	  numberCountTotal += numberCountV.data[k];
      }
      
      eta=numberCountTotal/mObsCoh;
      INT4 j;
      for(j=0 ; j<(uvar_p) ; j++){
	  
	  nj=numberCountV.data[j];
	  sumWeightj=chi2Params.sumWeight[j];
	  sumWeightSquarej=chi2Params.sumWeightSquare[j];
	  
	  chi2 += (nj-sumWeightj*eta)*(nj-sumWeightj*eta)/(sumWeightSquarej*eta*(1-eta));
      }
   }


  
   
    fp = fopen(uvar_outfile , "w");
    setvbuf(fp, (char *)NULL, _IOLBF, 0);
    fprintf(fp, "%g  %g  %g  %g %g  %g  %g  %g \n", (numberCountTotal - meanN)/sigmaN, meanN ,sigmaN, chi2, uvar_Freq, uvar_Alpha, uvar_Delta, uvar_fdot);
/*    fprintf(stdout, "%g  %g  %g  %g %g  %g  %g  %g \n", (numberCountTotal - meanN)/sigmaN, meanN ,sigmaN, chi2, uvar_Freq, uvar_Alpha, uvar_Delta, uvar_fdot);*/
    fclose(fp);
  


  /* free memory */
  LALFree(pulsarTemplate.spindown.data);  
  LALFree(timeV.data);
  LALFree(timeDiffV.data);
  LALFree(foft.data);
  LALFree(velV.data);

  LALFree(weightsV.data);

  XLALDestroyMultiDetectorStateSeries ( mdetStates );

  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);

  LAL_CALL (LALDestroyMultiSFTVector(&status, &inputSFTs), &status );
  LALFree(pg1.data);
  
  LALFree(numberCountV.data);

  LALFree(chi2Params.numberSFTp);
  LALFree(chi2Params.sumWeight);
  LALFree(chi2Params.sumWeightSquare);

  LAL_CALL (LALDestroyUserVars(&status), &status);

  LALCheckMemoryLeaks();

  if ( lalDebugLevel )
    REPORTSTATUS ( &status);

  return status.statusCode;
}



  


void ComputeFoft(LALStatus   *status,
		 REAL8Vector          *foft,
                 HoughTemplate        *pulsarTemplate,
		 REAL8Vector          *timeDiffV,
		 REAL8Cart3CoorVector *velV,
                 REAL8                 timeBase){
  
  INT4   mObsCoh;
  REAL8   f0new, vcProdn, timeDiffN;
  INT4    f0newBin;
  REAL8   sourceDelta, sourceAlpha, cosDelta;
  INT4    j,i, nspin, factorialN; 
  REAL8Cart3Coor  sourceLocation;
  
  /* --------------------------------------------- */
  INITSTATUS(status);
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
    f0newBin = floor( f0new * timeBase + 0.5);
    foft->data[j] = f0newBin * (1.0 +vcProdn) / timeBase;
  }    
    
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}







void SplitSFTs(LALStatus         *status,
	       REAL8Vector       *weightsV,
	       HoughParamsTest   *chi2Params){
  
    UINT4    j=0;           /* index of each block. It runs betwen 0 and p */ 
    UINT4   iSFT=0;       
    REAL8   *weights_ptr;  /* pointer to weightsV.data */
    REAL8   sumWeightpMax; /* Value of sumWeight we want to fix in each set of SFTs */
    UINT4   numberSFT;     /* Counter with the # of SFTs in each block */        
    UINT4   mObsCoh, p;
    REAL8   partialsumWeightp, partialsumWeightSquarep;
  
  /* --------------------------------------------- */
  INITSTATUS(status);
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

  sumWeightpMax= (REAL8)(mObsCoh)/p;       /* Compute the value of the sumWeight we want to fix in each set of SFT's */
  weights_ptr=weightsV->data;    /* Make the pointer to point to the first position of the vector weightsV.data */

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
    
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}



   
