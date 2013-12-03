/*
 *  Copyright (C) 2005  Alicia Sintes, Badri Krishnan, 
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
 * \author Alicia Sintes, Badri Krishnan
 * \brief
 * Monte Carlo signal injections for several h_0 values and
 * compute the Hough transform for a small number of point in parameter space each time
 */

/* 
 The idea is that we would like to analize a 300 Hz band on a cluster of
 machines. Each process should analyze 1 Hz band  (or whatever).
 
 	- Read the  band to be analized and the wings needed to read the originals SFTs. 
	-Read the h_0 values to be analyzed in one go
	-loop over the MC injections:
		+ Generate random parameters (f, f', alpha, delata, i...)
		+ generate h(t), produce its FFT
		+ Add h(f) to SFT for a given h_o value (and all of them)
		+ get number count
	
  Input shoud be from
             SFT files 
	     band,  nh_0, h_01, h_02....
	     ephemeris info
         	     
   This code will output files containing the MC results and info about injected
   signals. 
*/

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include "./MCInjectHoughMulti.h" /* proper path*/


#define EARTHEPHEMERIS "./earth05-09.dat" 
#define SUNEPHEMERIS "./sun05-09.dat"    

#define MAXFILENAMELENGTH 512 /* maximum # of characters  of a filename */

#define ACCURACY 0.00000001 /* of the velocity calculation */
#define MAXFILES 3000 /* maximum number of files to read in a directory */
 
#define THRESHOLD 1.6 /* thresold for peak selection, with respect to the
                              the averaged power in the search band */
#define F0 250.0          /*  frequency to build the LUT and start search */
#define FBAND 2.0          /* search frequency band  (in Hz) */
#define ALPHA 0.0		/* center of the sky patch (in radians) */
#define DELTA  (-LAL_PI_2)
#define PATCHSIZEX (LAL_PI*0.99) /* patch size */
#define PATCHSIZEY (LAL_PI*0.99)
#define NFSIZE  21 /* n-freq. span of the cylinder, to account for spin-down
                          search */
#define BLOCKSRNGMED 101 /* Running median window size */
#define NH0 2 /* number of h0 values to be analyzed */
#define H0MIN 1.0e-23
#define H0MAX 1.0e-22
#define NMCLOOP 2 /* number of Monte-Carlos */
#define NTEMPLATES 16 /* number templates for each Monte-Carlo */
#define DTERMS 5

#define SFTDIRECTORY "/local_data/sintes/SFT-S5-120-130/*SFT*.*"
/* */
#define DIROUT "./"   /* output directory */
#define FILEOUT "/HoughMC"      /* prefix file output */ 

#define TRUE (1==1)
#define FALSE (1==0)


static REAL8Vector empty_REAL8Vector;

/******************************************/
void ComputeFoft_NM(LALStatus   *status,
                 REAL8Vector          *foft,
                 HoughTemplate        *pulsarTemplate,
                 REAL8Vector          *timeDiffV,
                 REAL8Cart3CoorVector *velV);


void PrintLogFile2(LALStatus       *status, 
                  CHAR            *dir, 
		  CHAR            *basename,  
		  LALStringVector *linefiles,
		  CHAR            *executable );

/******************************************/


/* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv------------------------------------ */
int main(int argc, char *argv[]){

  static LALStatus            status;  
   
  EphemerisData              *edat = NULL;
  static REAL8Cart3CoorVector  velV;
  static LIGOTimeGPSVector     timeV;
  MultiLIGOTimeGPSVector     *multiIniTimeV = NULL;
  static REAL8Vector          timeDiffV;
  LIGOTimeGPS                firstTimeStamp, lastTimeStamp;
  REAL8                      tObs;

  static REAL8Vector          foft;
  static REAL8Vector          foftV[NTEMPLATES];
  static REAL8Vector          h0V;
 
  static HoughInjectParams    injectPar;
  static PulsarData           pulsarInject;
  static HoughTemplate        pulsarTemplate;
  static HoughNearTemplates   closeTemplates;
  SkyConstAndZeroPsiAMResponse      *pSkyConstAndZeroPsiAMResponse = NULL;
  SFTandSignalParams                *pSFTandSignalParams=NULL;
  

  /* standard pulsar sft types */ 
  MultiSFTVector *inputSFTs  = NULL;
  MultiSFTVector *sumSFTs    = NULL;
  MultiSFTVector *signalSFTs = NULL;
  
  /* information about all the ifos */
  MultiDetectorStateSeries *mdetStates = NULL;
  UINT4     numifo;
  UINT4     binsSFT;

  /* vector of weights */
  REAL8Vector      weightsV, weightsNoise,  weightsAM = empty_REAL8Vector;
  REAL8      alphaPeak, meanN, sigmaN, significance;
  
  REAL4TimeSeries   *signalTseries = NULL;
  static PulsarSignalParams  params;
  static SFTParams           sftParams;
  static UCHARPeakGram       pg1;
    
  UINT4  msp = 1; /*number of spin-down parameters */ 
  REAL8  numberCount,maxNumberCount;
  REAL8  numberCountV[NTEMPLATES];
  UINT4  nTemplates;
   
  UINT4  mObsCoh;
  REAL8  normalizeThr;
  REAL8  timeBase, deltaF;
  REAL8  h0scale;
 
  INT4   sftFminBin;
  REAL8  fHeterodyne;
  REAL8  tSamplingRate;      
 
  INT4 MCloopId;
  INT4 h0loop;
  
  FILE  *fpPar = NULL;
  FILE  *fpH0 = NULL;
  FILE  *fpNc = NULL;

 /* sft constraint variables */
  LIGOTimeGPS startTimeGPS, endTimeGPS;
  LIGOTimeGPSVector *inputTimeStampsVector=NULL;

  /******************************************************************/ 
  /*    user input variables   */
  /******************************************************************/ 
  BOOLEAN uvar_help, uvar_weighAM, uvar_weighNoise, uvar_printLog, uvar_fast;
  INT4    uvar_blocksRngMed, uvar_nh0, uvar_nMCloop, uvar_AllSkyFlag;
  INT4    uvar_nfSizeCylinder, uvar_maxBinsClean, uvar_Dterms;
  REAL8   uvar_f0, uvar_fSearchBand, uvar_peakThreshold, uvar_h0Min, uvar_h0Max;
  REAL8   uvar_alpha, uvar_delta, uvar_patchSizeAlpha, uvar_patchSizeDelta;
  CHAR   *uvar_earthEphemeris=NULL;
  CHAR   *uvar_sunEphemeris=NULL;
  CHAR   *uvar_sftDir=NULL;
  CHAR   *uvar_dirnameOut=NULL;
  CHAR   *uvar_fnameOut=NULL;
  LALStringVector *uvar_linefiles=NULL;
  REAL8    uvar_startTime, uvar_endTime;
  CHAR     *uvar_timeStampsFile=NULL;
  /******************************************************************/ 
  /*  set up the default parameters  */
  /******************************************************************/ 
  /* LAL error-handler */
  lal_errhandler   =   LAL_ERR_EXIT;
  
  uvar_help = FALSE;
  uvar_AllSkyFlag = 1;
  
  uvar_weighAM = TRUE;
  uvar_weighNoise = TRUE;
  uvar_printLog = FALSE;
  uvar_fast = TRUE;
  
  uvar_nh0 = NH0;
  uvar_h0Min = H0MIN;
  uvar_h0Max = H0MAX;

  uvar_nMCloop = NMCLOOP;
  nTemplates = NTEMPLATES;  
  uvar_alpha = ALPHA;
  uvar_delta = DELTA;
  uvar_f0 =  F0;
  uvar_fSearchBand = FBAND;
  uvar_peakThreshold = THRESHOLD;
  uvar_nfSizeCylinder = NFSIZE;
  uvar_blocksRngMed = BLOCKSRNGMED;
  uvar_maxBinsClean = 100;
  uvar_Dterms = DTERMS;

  uvar_patchSizeAlpha = PATCHSIZEX;
  uvar_patchSizeDelta = PATCHSIZEY; 
  
  uvar_earthEphemeris = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_earthEphemeris,EARTHEPHEMERIS);

  uvar_sunEphemeris = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_sunEphemeris,SUNEPHEMERIS);

  uvar_sftDir = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_sftDir,SFTDIRECTORY);

  uvar_dirnameOut = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_dirnameOut,DIROUT);

  uvar_fnameOut = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_fnameOut, FILEOUT);
 

  /******************************************************************/ 
  /*      register user input variables    */
  /******************************************************************/ 
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "help",            'h', UVAR_HELP,     "Print this message",            &uvar_help),            &status);  
  LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed",    'w', UVAR_OPTIONAL, "RngMed block size",             &uvar_blocksRngMed),    &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "f0",              'f', UVAR_OPTIONAL, "Start search frequency",        &uvar_f0),              &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fSearchBand",     'b', UVAR_OPTIONAL, "Search frequency band",         &uvar_fSearchBand),     &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "peakThreshold",   't', UVAR_OPTIONAL, "Peak selection threshold",      &uvar_peakThreshold),   &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "earthEphemeris",  'E', UVAR_OPTIONAL, "Earth Ephemeris file",          &uvar_earthEphemeris),  &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sunEphemeris",    'S', UVAR_OPTIONAL, "Sun Ephemeris file",            &uvar_sunEphemeris),    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftDir",          'D', UVAR_OPTIONAL, "SFT Directory",                 &uvar_sftDir),          &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "dirnameOut",      'o', UVAR_OPTIONAL, "Output directory",                      &uvar_dirnameOut),      &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "fnameout",        '0', UVAR_OPTIONAL, "Output file prefix",            &uvar_fnameOut),        &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "alpha",           'r', UVAR_OPTIONAL, "Right ascension",               &uvar_alpha),           &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "delta",           'l', UVAR_OPTIONAL, "Declination",                   &uvar_delta),           &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "patchSizeAlpha",  'R', UVAR_OPTIONAL, "Patch size in right ascension", &uvar_patchSizeAlpha),  &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "patchSizeDelta",  'L', UVAR_OPTIONAL, "Patch size in declination",     &uvar_patchSizeDelta),  &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "patch",           'P', UVAR_OPTIONAL, "Inject in patch if 0",          &uvar_AllSkyFlag),      &status);  
  LAL_CALL( LALRegisterINTUserVar(    &status, "nMCloop",         'N', UVAR_OPTIONAL, "Number of MC injections",       &uvar_nMCloop),         &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "h0Min",           'm', UVAR_OPTIONAL, "Smallest h0 to inject",         &uvar_h0Min),           &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "h0Max",           'M', UVAR_OPTIONAL, "Largest h0 to inject",          &uvar_h0Max),           &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "nh0",             'n', UVAR_OPTIONAL, "Number of h0 values to inject", &uvar_nh0),             &status);  
  LAL_CALL( LALRegisterLISTUserVar(   &status, "linefiles",        0,  UVAR_OPTIONAL, "list of linefiles separated by commas", &uvar_linefiles),       &status);
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "weighAM",          0,  UVAR_OPTIONAL, "Use amplitude modulation weights",      &uvar_weighAM),         &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "weighNoise",       0,  UVAR_OPTIONAL, "Use SFT noise weights",                 &uvar_weighNoise),      &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "printLog",         0,  UVAR_OPTIONAL, "Print Log file",                        &uvar_printLog),        &status);  
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "fast",             0,  UVAR_OPTIONAL, "Use fast frequency domain SFT injections",    &uvar_fast),      &status);  
  LAL_CALL( LALRegisterREALUserVar(   &status, "startTime",        0,  UVAR_OPTIONAL, "GPS start time of observation",         &uvar_startTime),        &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "endTime",          0,  UVAR_OPTIONAL, "GPS end time of observation",           &uvar_endTime),          &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "timeStampsFile",   0,  UVAR_OPTIONAL, "Input time-stamps file",                &uvar_timeStampsFile),   &status);
  /* developer input variables */
  LAL_CALL( LALRegisterINTUserVar(    &status, "nfSizeCylinder",   0, UVAR_DEVELOPER, "Size of cylinder of PHMDs",             &uvar_nfSizeCylinder),  &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed",     0, UVAR_DEVELOPER, "Running Median block size",             &uvar_blocksRngMed),    &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "maxBinsClean",     0, UVAR_DEVELOPER, "Maximum number of bins in cleaning",    &uvar_maxBinsClean),    &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "Dterms",           0,  UVAR_DEVELOPER, "Number of f-bins in MC injection",     &uvar_Dterms),    &status);

  /******************************************************************/ 
  /* read all command line variables */
  /******************************************************************/ 
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
  
  LAL_CALL( LALRngMedBias( &status, &normalizeThr, uvar_blocksRngMed ), &status ); 

  /******************************************************************/ 
  /* write log file with command line arguments, cvs tags, and contents of skypatch file */
  /******************************************************************/ 
  if ( uvar_printLog ) {
    LAL_CALL( PrintLogFile2( &status, uvar_dirnameOut, uvar_fnameOut, uvar_linefiles, argv[0]), &status);
  }
   
  /******************************************************************/ 
  /* set fullsky flag */
  /******************************************************************/ 

  injectPar.fullSky = 1;
  if (uvar_AllSkyFlag == 0)
    injectPar.fullSky= 0;  /* patch case */  

  /******************************************************************/ 
  /* computing h0 values  */
  /******************************************************************/ 
  h0V.length=uvar_nh0;
  h0V.data = NULL;
  h0V.data = (REAL8 *)LALMalloc(uvar_nh0*sizeof(REAL8));
  h0V.data[0] = uvar_h0Min;
  
  if(uvar_nh0 >1){
    INT4 k;
    REAL8 steph0;   
    steph0 = (uvar_h0Max-uvar_h0Min)/(uvar_nh0-1.);
    for(k=1; k<uvar_nh0; ++k) h0V.data[k]= h0V.data[k-1]+steph0;
  }
  
  /******************************************************************/ 
  /*  preparing  output files */
  /******************************************************************/ 
  {
    INT4 k;
    CHAR filename[MAXFILENAMELENGTH];
    
    /* the paramerter file */
    strcpy( filename, uvar_dirnameOut);
     
    strcat( filename, uvar_fnameOut);
    strcat( filename, "_par");
    fpPar= fopen(filename, "w"); /* where to write the parameters */
    /*setlinebuf(fpPar);*/  /* line buffered on */
    setvbuf(fpPar, (char *)NULL, _IOLBF, 0);
    
    /* the  file  with the h0 values */
    strcpy( filename, uvar_dirnameOut);
    strcat( filename, uvar_fnameOut);
    strcat( filename, "_h0");
    fpH0= fopen(filename, "w"); /* where to write the parameters */
    /*setlinebuf(fpH0); */ /* line buffered on */
    setvbuf(fpH0, (char *)NULL, _IOLBF, 0); 
   
    /* the  file  with the the number-counts for different h0 values */
    strcpy( filename, uvar_dirnameOut);
    strcat( filename, uvar_fnameOut);
    strcat( filename, "_nc");
    fpNc= fopen(filename, "w"); /* where to write the parameters */
    /*setlinebuf(fpNc);*/  /* line buffered on */
    setvbuf(fpNc, (char *)NULL, _IOLBF, 0);

    for (k=0; k<uvar_nh0; ++k){ fprintf(fpH0, "%g \n",  h0V.data[k] ); }  
    fclose(fpH0);
    
  }

  /******************************************************************/ 
  /* sft reading */
  /******************************************************************/ 
 
  {
    /* new SFT I/O data types */
    SFTCatalog *catalog = NULL;
    static SFTConstraints constraints;

    REAL8   doppWings, f_min, f_max;
 
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

    /* catalog is ordered in time so we can get start, end time and tObs*/
    firstTimeStamp = catalog->data[0].header.epoch;
    lastTimeStamp = catalog->data[mObsCoh - 1].header.epoch;
    tObs = XLALGPSDiff( &lastTimeStamp, &firstTimeStamp ) + timeBase;

    /* add wings for Doppler modulation and running median block size*/
    doppWings = (uvar_f0 + uvar_fSearchBand) * VTOT;    
    f_min = uvar_f0 - doppWings - (uvar_blocksRngMed + uvar_nfSizeCylinder) * deltaF;
    f_max = uvar_f0 + uvar_fSearchBand + doppWings + (uvar_blocksRngMed + uvar_nfSizeCylinder) * deltaF;

    /* read sft files making sure to add extra bins for running median */
    LAL_CALL( LALLoadMultiSFTs ( &status, &inputSFTs, catalog, f_min, f_max), &status);
 
    /* SFT info -- assume all SFTs have same length */
    numifo = inputSFTs->length;
    binsSFT = inputSFTs->data[0]->data->data->length;
     
    /* some more sft parameetrs */  
    fHeterodyne =inputSFTs->data[0]->data[0].f0;   
    sftFminBin = (INT4) floor(fHeterodyne *timeBase +0.5);    
    tSamplingRate = 2.0*deltaF*(binsSFT -1.);
         
    /* free memory */    
    LAL_CALL( LALDestroySFTCatalog( &status, &catalog ), &status);  	 
  } 
  
  /******************************************************************/  
  /* allocate memory for sumSFTs of the same size of inputSFTs and place for
  signal only sfts. */
  /******************************************************************/ 
 
  {
    UINT4Vector  numsft;
    UINT4        iIFO;
    
    numsft.length =  numifo;
    numsft.data = NULL;
    numsft.data =(UINT4 *)LALCalloc(numifo, sizeof(UINT4));
       
    for ( iIFO = 0; iIFO < numifo; iIFO++) {
      numsft.data[iIFO] = inputSFTs->data[iIFO]->length;     
    }
    
    LAL_CALL( LALCreateMultiSFTVector(&status, &sumSFTs, binsSFT, &numsft), &status);
    LAL_CALL( LALCreateMultiSFTVector(&status, &signalSFTs, binsSFT, &numsft), &status);
    LALFree( numsft.data);
     
  }
   
  /******************************************************************/  
  /* allocate memory for velocity vector and timestamps */
  /******************************************************************/ 
  
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
  
   /* start time of the sfts sorted for the different IFOs */
  multiIniTimeV = (MultiLIGOTimeGPSVector *)LALMalloc(sizeof(MultiLIGOTimeGPSVector));
  multiIniTimeV->length = numifo;
  multiIniTimeV->data = (LIGOTimeGPSVector **)LALCalloc(numifo, sizeof(LIGOTimeGPSVector *));

  /******************************************************************/ 
  /* get detector velocities and timestamps */
  /******************************************************************/ 
  
    /*  setting of ephemeris info */ 
  edat = (EphemerisData *)LALMalloc(sizeof(EphemerisData));
  (*edat).ephiles.earthEphemeris = uvar_earthEphemeris;
  (*edat).ephiles.sunEphemeris = uvar_sunEphemeris;

  {
    UINT4   iIFO, iSFT, numsft, j;

    LAL_CALL( LALInitBarycenter( &status, edat), &status);
    
    /* get information about all detectors including velocity and timestamps */
    /* note that this function returns the velocity at the 
       mid-time of the SFTs --CAREFULL later on with the time stamps!!! velocity
       is ok */
    
    LAL_CALL ( LALGetMultiDetectorStates ( &status, &mdetStates, inputSFTs, edat), &status);
    
    /* copy the timestamps and velocity vector */
    for (j = 0, iIFO = 0; iIFO < numifo; iIFO++ ) {
      numsft = mdetStates->data[iIFO]->length;      
      for ( iSFT = 0; iSFT < numsft; iSFT++, j++) {
	velV.data[j].x = mdetStates->data[iIFO]->data[iSFT].vDetector[0];
	velV.data[j].y = mdetStates->data[iIFO]->data[iSFT].vDetector[1];
	velV.data[j].z = mdetStates->data[iIFO]->data[iSFT].vDetector[2];
	/* mid time of sfts */
        timeV.data[j] = mdetStates->data[iIFO]->data[iSFT].tGPS;
      } /* loop over SFTs */
      
      multiIniTimeV->data[iIFO]=NULL;
      LAL_CALL( LALGetSFTtimestamps(&status, &multiIniTimeV->data[iIFO],
                                     inputSFTs->data[iIFO] ), &status);
    
    } /* loop over IFOs */
    
    /* compute the time difference relative to startTime for all SFT */
    for(j = 0; j < mObsCoh; j++)
      timeDiffV.data[j] = XLALGPSDiff( timeV.data + j, &firstTimeStamp );
    
     /* removing mid time-stamps, no longer needed now */
    LALFree(timeV.data); 
   
  }
 

  /******************************************************************/ 
  /* probability of selecting a peak and expected mean for noise only */
  /******************************************************************/ 
  alphaPeak = exp( - uvar_peakThreshold);
  meanN = mObsCoh* alphaPeak;

  /******************************************************************/ 
  /*  initialization of fast injections parameters  TO BE FIXED for multiIFO
     and memory should be dealocated at the end */ 
  /******************************************************************/ 

  pSkyConstAndZeroPsiAMResponse = (SkyConstAndZeroPsiAMResponse *)
                   LALMalloc(sizeof(SkyConstAndZeroPsiAMResponse)*numifo);
  {	
    UINT4 numsft,iIFO;
     
     for (iIFO=0; iIFO<numifo; iIFO++){
       
       numsft = mdetStates->data[iIFO]->length;   
         
       pSkyConstAndZeroPsiAMResponse[iIFO].skyConst = (REAL8 *)LALMalloc((2*msp*(numsft+1)+2*numsft+3)*sizeof(REAL8));
       pSkyConstAndZeroPsiAMResponse[iIFO].fPlusZeroPsi = (REAL4 *)LALMalloc(numsft*sizeof(REAL4));
       pSkyConstAndZeroPsiAMResponse[iIFO].fCrossZeroPsi = (REAL4 *)LALMalloc(numsft*sizeof(REAL4));
      }
    }

  {    
    INT4 k;
    
    pSFTandSignalParams = (SFTandSignalParams *)LALMalloc(sizeof(SFTandSignalParams));
        /* create lookup table (LUT) values for doing trig */
        /* pSFTandSignalParams->resTrig = 64; */ /* length sinVal and cosVal; resolution of trig functions = 2pi/resTrig */
    pSFTandSignalParams->resTrig = 0; /* 08/02/04 gam; avoid serious bug when using LUT for trig calls */
    pSFTandSignalParams->trigArg = (REAL8 *)LALMalloc((pSFTandSignalParams->resTrig+1)*sizeof(REAL8));
    pSFTandSignalParams->sinVal  = (REAL8 *)LALMalloc((pSFTandSignalParams->resTrig+1)*sizeof(REAL8));
    pSFTandSignalParams->cosVal  = (REAL8 *)LALMalloc((pSFTandSignalParams->resTrig+1)*sizeof(REAL8));

    for (k=0; k<=pSFTandSignalParams->resTrig; k++) {
       pSFTandSignalParams->trigArg[k]= ((REAL8)LAL_TWOPI) * ((REAL8)k) / ((REAL8)pSFTandSignalParams->resTrig);
       pSFTandSignalParams->sinVal[k]=sin( pSFTandSignalParams->trigArg[k] );
       pSFTandSignalParams->cosVal[k]=cos( pSFTandSignalParams->trigArg[k] );
    }

    pSFTandSignalParams->pSigParams = &params;    /* as defined in Hough*/
    pSFTandSignalParams->pSFTParams = &sftParams; /* as defined in Hough*/
    pSFTandSignalParams->nSamples = (INT4)(0.5*timeBase);  /* nsample to get version 2 sfts */
    pSFTandSignalParams->Dterms = uvar_Dterms; 

  }
  
  /******************************************************************/ 
  /*   setting of parameters */ 
  /******************************************************************/ 
  injectPar.h0   = uvar_h0Min;
  injectPar.fmin = uvar_f0;
  injectPar.fSearchBand = uvar_fSearchBand;
  injectPar.deltaF = deltaF;
  injectPar.alpha = uvar_alpha;  /* patch center if not full sky */
  injectPar.delta = uvar_delta;
  injectPar.patchSizeAlpha = uvar_patchSizeAlpha; /* patch size if not full sky */
  injectPar.patchSizeDelta = uvar_patchSizeDelta; 
  injectPar.pixelFactor = PIXELFACTOR;
  injectPar.vTotC = VTOT;
  injectPar.timeObs =tObs;  
  injectPar.spnFmax.data = NULL; 
  injectPar.spnFmax.length=msp;   /*only 1 spin */
  injectPar.spnFmax.data = (REAL8 *)LALMalloc(msp*sizeof(REAL8));
  injectPar.spnFmax.data[0] = -(uvar_nfSizeCylinder/2) *deltaF/tObs;
  
  pulsarInject.spindown.length = msp;
  pulsarTemplate.spindown.length = msp; 
  pulsarInject.spindown.data = NULL;
  pulsarTemplate.spindown.data = NULL; 
  pulsarInject.spindown.data = (REAL8 *)LALMalloc(msp*sizeof(REAL8));
  pulsarTemplate.spindown.data = (REAL8 *)LALMalloc(msp*sizeof(REAL8));
 
  sftParams.Tsft = timeBase;
  sftParams.noiseSFTs = NULL;       
  
  params.orbit = NULL;
  /* params.transferFunction = NULL; */
  params.ephemerides = edat;
  params.startTimeGPS.gpsSeconds = firstTimeStamp.gpsSeconds;   /* start time of output time series */
  params.startTimeGPS.gpsNanoSeconds = firstTimeStamp.gpsNanoSeconds;   /* start time of output time series */
  params.duration = injectPar.timeObs; /* length of time series in seconds */
  params.samplingRate = tSamplingRate;
  params.fHeterodyne = fHeterodyne;  
  params.pulsar.refTime.gpsSeconds = firstTimeStamp.gpsSeconds; 
  params.pulsar.refTime.gpsNanoSeconds = firstTimeStamp.gpsNanoSeconds; 
  /* ****************************************************************/
  
  /* WE SHOULD LOOP OVER MC SIGNAL INJECTION HERE
     BEFORE THAT :
     -for each different h0 value create a file containing the h0 value
     LOOP over xxx Monte-Carlo signal Injections:
	- Generate signal injections parameters and the corresponding template
		parameters allowing some mismatch
	- Compute the frequency path for the template parameters
	-Generate the time series for injected signals and the
		corresponding SFTs with no added noise (for all times).
		 
	LOOP over the different h0 values
		     -Add SFT with the signal normalized to the SFT original noise
		     -clean lines, compute weights, select peaks
		     -get the significance
   */
  
  /* ****************************************************************/
  
  pg1.length = binsSFT; /*equal to binsSFT */
  pg1.data = NULL;
  pg1.data = (UCHAR *)LALMalloc(binsSFT* sizeof(UCHAR));
  
  /* ****************************************************************/
  foft.length = mObsCoh;
  foft.data = NULL;
  foft.data = (REAL8 *)LALMalloc(mObsCoh*sizeof(REAL8));
  {
    UINT4 j;
    for (j=0;j<nTemplates;++j) {
      foftV[j].length = mObsCoh;
      foftV[j].data = NULL;
      foftV[j].data = (REAL8 *)LALMalloc(mObsCoh*sizeof(REAL8));
    }
  } 
  
  /* ****************************************************************/
  /*  HERE SHOULD START THE MONTE-CARLO */
  
  for(MCloopId=0; MCloopId < uvar_nMCloop; ++MCloopId){
     
    LAL_CALL( GenerateInjectParamsNoVeto(&status, &pulsarInject, &pulsarTemplate,
					 &closeTemplates, &injectPar), &status );
    	     
    /* writing the parameters into fpPar, following the format
       MCloopId  I.f0 H.f0 I.f1 H.f1 I.alpha H.alpha I.delta H.delta I.phi0  I.psi
       (not cos iota)  */

    fprintf(fpPar," %d %f %f %g %g %f %f %f %f %f %f %f", 
	    MCloopId, pulsarInject.f0, pulsarTemplate.f0,
	    pulsarInject.spindown.data[0], pulsarTemplate.spindown.data[0],
	    pulsarInject.longitude, pulsarTemplate.longitude,
	    pulsarInject.latitude, pulsarTemplate.latitude,
	    pulsarInject.phi0, pulsarInject.psi, (pulsarInject.aCross)/injectPar.h0
	    );

	     
   /* ****************************************************************/
   /* Computing the frequency path f(t) = f0(t)* (1+v/c.n)  for all the different templates */
	     
   /* the geometrically nearest template */ 
   LAL_CALL( ComputeFoft_NM(&status, &foft,&pulsarTemplate,&timeDiffV,&velV), &status);

   /* for all the 16 near templates */
   {
     UINT4 j,i,k, itemplate;
     
     itemplate =0;
     for(j=0;j<2;++j){
       pulsarTemplate.f0 = closeTemplates.f0[j];
       for(i=0;i<2;++i){
	 pulsarTemplate.spindown.data[0] = closeTemplates.f1[i];
	 for(k=0;k<4;++k){
	   pulsarTemplate.latitude = closeTemplates.skytemp[k].delta;
	   pulsarTemplate.longitude = closeTemplates.skytemp[k].alpha;
	   LAL_CALL( ComputeFoft_NM(&status, &(foftV[itemplate]),
				 &pulsarTemplate,&timeDiffV,&velV), &status);
	   ++itemplate;
	 }
       }
     }
   }
   
   /* ****************************************************************/
		  	     
   /*  params.pulsar.TRefSSB=  ? ; */
   params.pulsar.position.longitude = pulsarInject.longitude;
   params.pulsar.position.latitude =  pulsarInject.latitude ;
   params.pulsar.position.system= COORDINATESYSTEM_EQUATORIAL; 
   params.pulsar.psi=    pulsarInject.psi;
   params.pulsar.aPlus=  pulsarInject.aPlus;
   params.pulsar.aCross= pulsarInject.aCross;
   params.pulsar.phi0=   pulsarInject.phi0;
   params.pulsar.f0=     pulsarInject.f0;
   params.pulsar.spindown=  &pulsarInject.spindown ;   
    
   {
     UINT4 iIFO, numsft, iSFT, j;    
           
     if(uvar_fast){
     
       for (iIFO=0; iIFO<numifo; iIFO++){       
         params.site = &(mdetStates->data[iIFO]->detector);
         sftParams.timestamps = multiIniTimeV->data[iIFO];
	 numsft = mdetStates->data[iIFO]->length; 
	 
	 /* initialize data to zero */
         for ( iSFT = 0; iSFT < numsft; iSFT++){	   
	   for (j=0; j < binsSFT; j++) {
	     signalSFTs->data[iIFO]->data[iSFT].data->data[j].realf_FIXME = 0.0;
	     signalSFTs->data[iIFO]->data[iSFT].data->data[j].imagf_FIXME = 0.0;	    
	   }	 
         }
     	  	 
	 LAL_CALL( LALComputeSkyAndZeroPsiAMResponse (&status,
	  &pSkyConstAndZeroPsiAMResponse[iIFO], pSFTandSignalParams), &status);
         LAL_CALL( LALFastGeneratePulsarSFTs (&status, &signalSFTs->data[iIFO],
	  &pSkyConstAndZeroPsiAMResponse[iIFO], pSFTandSignalParams), &status);	 
       }
     }
     else{
     
       for (iIFO=0; iIFO<numifo; iIFO++){
         params.site = &(mdetStates->data[iIFO]->detector);
         sftParams.timestamps = multiIniTimeV->data[iIFO];
	 
         LAL_CALL(LALDestroySFTVector(&status, &signalSFTs->data[iIFO]),&status );
         signalSFTs->data[iIFO] = NULL;
	 
         LAL_CALL( LALGeneratePulsarSignal(&status, &signalTseries, &params ), &status);
         LAL_CALL( LALSignalToSFTs(&status, &signalSFTs->data[iIFO], signalTseries, &sftParams), &status);
           
         LALFree(signalTseries->data->data);
         LALFree(signalTseries->data);
         LALFree(signalTseries);
         signalTseries =NULL;      
       }
     }   
   }
   
 
  /******************************************************************/ 
  /* initialize all weights to unity */
  /******************************************************************/ 
  
  /* set up weights -- this should be done before normalizing the sfts */
  weightsV.length = mObsCoh;
  weightsV.data = (REAL8 *)LALCalloc(mObsCoh, sizeof(REAL8));
  
  weightsNoise.length = mObsCoh;
  weightsNoise.data = (REAL8 *)LALCalloc(mObsCoh, sizeof(REAL8));
  
  /* initialize all weights to unity */
  LAL_CALL( LALHOUGHInitializeWeights( &status, &weightsNoise), &status);
  LAL_CALL( LALHOUGHInitializeWeights( &status, &weightsV), &status);
  
  /******************************************************************/ 
  /*   setting the weights considering only the AM coefficients to be only
       computed just where we need*/ 
  /******************************************************************/ 
  if (uvar_weighAM){
    SkyPosition      skypos;
    UINT4            iIFO, iSFT;
    UINT4 	      k, numsft;
    MultiAMCoeffs   *multiAMcoef = NULL;
    
    weightsAM.length = mObsCoh;
    weightsAM.data = (REAL8 *)LALCalloc(mObsCoh, sizeof(REAL8));
    skypos.system = COORDINATESYSTEM_EQUATORIAL;
    
    skypos.longitude = pulsarInject.longitude;
    skypos.latitude  = pulsarInject.latitude;
    LAL_CALL ( LALGetMultiAMCoeffs ( &status, &multiAMcoef, mdetStates, skypos), &status);
      
    /* loop over the weights and set them by the appropriate AM coefficients */
    for ( k = 0, iIFO = 0; iIFO < numifo; iIFO++) {	  
      numsft = mdetStates->data[iIFO]->length;	
      for ( iSFT = 0; iSFT < numsft; iSFT++, k++) {	  
	REAL8 a, b;
	  
	a = multiAMcoef->data[iIFO]->a->data[iSFT];
	b = multiAMcoef->data[iIFO]->b->data[iSFT];    
	weightsAM.data[k] = (a*a + b*b);
      } /* loop over SFTs */
    } /* loop over IFOs */
      
    XLALDestroyMultiAMCoeffs ( multiAMcoef );
      
  }

	     
   /* ****************************************************************/
   /*  HERE THE LOOP FOR DIFFERENT h0 VALUES */
	     
    fprintf(fpNc, " %d ",  MCloopId);
	     
    for(h0loop=0; h0loop <uvar_nh0; ++h0loop){
      
      UINT4       j, ind, itemplate; 
      UINT4       iIFO, numsft, iSFT;
      COMPLEX8   *noiseSFT;
      COMPLEX8   *signalSFT;
      COMPLEX8   *sumSFT;
      
      numberCount=0.0;
      for(itemplate=0; itemplate<nTemplates; ++itemplate){
        numberCountV[itemplate]=0.0;
      }
      
      h0scale =h0V.data[h0loop]/h0V.data[0]; /* different for different h0 values */
      
      /* ****************************************************************/
      /* adding signal+ noise SFT, TO BE CHECKED */ 
           
      for (iIFO=0; iIFO<numifo; iIFO++){
        numsft =  inputSFTs->data[iIFO]->length;
	for ( iSFT = 0; iSFT < numsft; iSFT++){
	  
	  noiseSFT  = inputSFTs->data[iIFO]->data[iSFT].data->data;
	  signalSFT = signalSFTs->data[iIFO]->data[iSFT].data->data;
	  sumSFT    = sumSFTs->data[iIFO]->data[iSFT].data->data;
	  	  
	  for (j=0; j < binsSFT; j++) {
	    sumSFT->realf_FIXME = crealf(*noiseSFT) + h0scale *crealf(*signalSFT);
	    sumSFT->imagf_FIXME = cimagf(*noiseSFT) + h0scale *cimagf(*signalSFT);
	    ++noiseSFT;
	    ++signalSFT;
	    ++sumSFT;
	  }	 
	}
      }
      
      /* ****************************************************************/
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
	  
	  LAL_CALL( LALRemoveKnownLinesInMultiSFTVector ( &status, sumSFTs, uvar_maxBinsClean, uvar_blocksRngMed, uvar_linefiles, randPar), &status);
	  
	  LAL_CALL ( LALDestroyRandomParams (&status, &randPar), &status);
	  fclose(fpRand);
	} /* end cleaning */
      
      
      /* ****************************************************************/
      /* normalize sfts compute weights */
      {   
	MultiNoiseWeights *multweight = NULL;    
	MultiPSDVector *multPSD = NULL;  
	REAL8 sumWeightSquare;
	
	
	/* initialize all weights to unity  each time*/
        LAL_CALL( LALHOUGHInitializeWeights( &status, &weightsNoise), &status);
        LAL_CALL( LALHOUGHInitializeWeights( &status, &weightsV), &status);
  
	
	/* normalize sfts */
	LAL_CALL( LALNormalizeMultiSFTVect (&status, &multPSD, sumSFTs, uvar_blocksRngMed), &status);
	
	/* compute multi noise weights */ 
	if ( uvar_weighNoise ) {
 	  LAL_CALL ( LALComputeMultiNoiseWeights ( &status, &multweight, multPSD, uvar_blocksRngMed, 0), &status);
	}
	
	/* we are now done with the psd */
	LAL_CALL ( LALDestroyMultiPSDVector  ( &status, &multPSD), &status);
	
	/* copy  weights */
        if ( uvar_weighNoise ) {
	  for (j = 0, iIFO = 0; iIFO < numifo; iIFO++ ) {
	    numsft = mdetStates->data[iIFO]->length;
	    for ( iSFT = 0; iSFT < numsft; iSFT++, j++) {
	      weightsNoise.data[j] = multweight->data[iIFO]->data[iSFT];
	    } /* loop over SFTs */
	  } /* loop over IFOs */
      
	  LAL_CALL ( LALDestroyMultiNoiseWeights ( &status, &multweight), &status);
	  memcpy(weightsV.data, weightsNoise.data, mObsCoh * sizeof(REAL8));
        }
	
	if (uvar_weighAM && weightsAM.data ) {
	  for (j=0; j<mObsCoh; j++){
	    weightsV.data[j] = weightsV.data[j]*weightsAM.data[j];
	  }
	}
	
	LAL_CALL( LALHOUGHNormalizeWeights( &status, &weightsV), &status);
	
        /* calculate the sum of the weights squared */
        sumWeightSquare = 0.0;
        for ( j = 0; j < mObsCoh; j++)
          sumWeightSquare += weightsV.data[j] * weightsV.data[j];

        /* standard deviation for noise only */
        sigmaN = sqrt(sumWeightSquare * alphaPeak * (1.0 - alphaPeak));	
      }
      
      /* ****************************************************************/
      /* loop over SFT, generate peakgram and get number count */
      {
        SFTtype  *sft;

	for ( j = 0, iIFO = 0; iIFO < numifo; iIFO++){
	  numsft = mdetStates->data[iIFO]->length;
	  
	  for ( iSFT = 0; iSFT < numsft; iSFT++, j++) {
	    
	    sft = sumSFTs->data[iIFO]->data + iSFT;
	    LAL_CALL (SFTtoUCHARPeakGram( &status, &pg1, sft, uvar_peakThreshold), &status);	    
	    ind = floor( foft.data[j]*timeBase -sftFminBin+0.5); 
	    numberCount+=pg1.data[ind]*weightsV.data[j]; /* adds 0 or 1 to the counter*/
	    
	    for (itemplate=0; itemplate<nTemplates; ++itemplate) {
	      ind = floor( foftV[itemplate].data[j]*timeBase -sftFminBin+0.5); 
	      numberCountV[itemplate]+=pg1.data[ind]*weightsV.data[j];
	    }	    
	  } /* loop over SFTs */	  
	} /* loop over IFOs */

      }
      
      /* ****************************************************************/
      /*check the max number count */
      maxNumberCount = numberCount;
      for (itemplate=0; itemplate<nTemplates; ++itemplate) {
	if( numberCountV[itemplate] > maxNumberCount ) {
	  maxNumberCount = numberCountV[itemplate];
	}
      }

      /******************************************************************/
      /* printing the significance in the proper file */
      /******************************************************************/
      significance = (maxNumberCount - meanN)/sigmaN;
      fprintf(fpNc, " %f ", significance);
      
    } /* closing loop for different h0 values */
    fprintf(fpNc, " \n");


  LALFree(weightsV.data);
  LALFree(weightsNoise.data);
  if (uvar_weighAM && weightsAM.data){ 
      LALFree(weightsAM.data); 
  }

	        
  } /* Closing MC loop */
  
  /******************************************************************/
  /* Closing files */
  /******************************************************************/  
  fclose(fpPar); 
  fclose(fpNc); 

  
  /******************************************************************/
  /* Free memory and exit */
  /******************************************************************/
  
  /* LALFree(fp); */
  LALFree(pg1.data);
 
  LALFree(velV.data);
  LALFree(timeDiffV.data);
  
  LALFree(pSFTandSignalParams->cosVal);
  LALFree(pSFTandSignalParams->sinVal);
  LALFree(pSFTandSignalParams->trigArg);
  LALFree(pSFTandSignalParams);

  
  {
     UINT4 iIFO;
      
     for(iIFO = 0; iIFO<numifo; iIFO++){
        LAL_CALL(LALDestroyTimestampVector(&status, &multiIniTimeV->data[iIFO]),&status );
	LALFree(pSkyConstAndZeroPsiAMResponse[iIFO].skyConst);
	LALFree(pSkyConstAndZeroPsiAMResponse[iIFO].fPlusZeroPsi);
	LALFree(pSkyConstAndZeroPsiAMResponse[iIFO].fCrossZeroPsi);
      }
  }
  
  LALFree(multiIniTimeV->data);
  LALFree(multiIniTimeV);
  LALFree(pSkyConstAndZeroPsiAMResponse);
   
  XLALDestroyMultiDetectorStateSeries ( mdetStates );

  LALFree(foft.data);
  LALFree(h0V.data);
  
  {
     UINT4 j;
     for (j=0;j<nTemplates;++j) {
        LALFree(foftV[j].data);
     }
  }

  
  LALFree(injectPar.spnFmax.data);
  LALFree(pulsarInject.spindown.data);
  LALFree(pulsarTemplate.spindown.data);
   
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);
      
  LAL_CALL(LALDestroyMultiSFTVector(&status, &inputSFTs),&status );
  LAL_CALL(LALDestroyMultiSFTVector(&status, &sumSFTs),&status );
  LAL_CALL(LALDestroyMultiSFTVector(&status, &signalSFTs),&status );
 
  LAL_CALL (LALDestroyUserVars(&status), &status);  

  LALCheckMemoryLeaks();
  
  if ( lalDebugLevel )
    REPORTSTATUS ( &status);

  return status.statusCode;
}



 
/***************************************************************************/
void GenerateInjectParams(LALStatus   *status,
                        PulsarData           *injectPulsar,
                        HoughTemplate        *templatePulsar,
			HoughNearTemplates   *closeTemplates,
                        HoughInjectParams    *params,
			LineNoiseInfo        *lines  ){
			
  INT4          seed=0; /* seed generated using current time */
  REAL4         randval;
  RandomParams  *randPar=NULL;
  FILE     *fpRandom;
  INT4     count;
  
  REAL4    cosiota, h0;
  REAL8    f0, deltaF, deltaX;
  REAL8    latitude, longitude;  /* of the source in radians */
  INT8    f0bin;
  UINT4    msp;
  
  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  
  /*   Make sure the arguments are not NULL: */
  ASSERT (injectPulsar,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (templatePulsar, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (params, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (lines, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  
  /*  ++++++++++++++++++from makefakedata
   * Modified so as to not create random number parameters with seed
   * drawn from clock.  Seconds don't change fast enough and sft's
   * look alike.  We open /dev/urandom and read a 4 byte integer from
   * it and use that as our seed.  Note: /dev/random is slow after the
   * first, few accesses.
   */

  fpRandom = fopen("/dev/urandom","r");
  ASSERT (fpRandom, status, DRIVEHOUGHCOLOR_EFILE,  DRIVEHOUGHCOLOR_MSGEFILE); 
  
  count = fread(&seed, sizeof(INT4),1, fpRandom);
  if ( count != 0 ) ABORT( status, DRIVEHOUGHCOLOR_EARG,  DRIVEHOUGHCOLOR_MSGEARG);
  
  fclose(fpRandom);
  
  TRY( LALCreateRandomParams(status->statusPtr, &randPar, seed), status);
  
 /*
  *   to create a single random deviate distributed uniforly between zero and unity		     
  *   TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
  */
  
  
  /* get random value phi0 [0, 2 pi] */ 
  TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
  injectPulsar->phi0 = randval * LAL_TWOPI;
  
  /* get random value cos iota [-1,1] */ 
  TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
  cosiota = 2.0* randval -1.0;
  
  h0=params->h0;
  injectPulsar->aCross = h0*cosiota;
  injectPulsar->aPlus  = 0.5*h0*(1.0 + cosiota*cosiota);
  
  /* get random value psi [0, 2 pi] */ 
  TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
  injectPulsar->psi = randval * LAL_TWOPI;

  /* getting random number for the frequency (and mismatch)*/
  TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
  f0 = params->fmin + (params->fSearchBand) * randval;
  
  /* veto the frequency if it is affected by a line */
  {
    INT4 veto=1;
    while( veto > 0 ){

      TRY( LALCheckLines (status->statusPtr, &veto, lines, f0 ), status); 
      if ( veto > 0 )
	{
	  TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
	  f0 = params->fmin + (params->fSearchBand) * randval;
	}
    } /* end of while loop */
  }
   
  injectPulsar->f0 = f0;
  deltaF = params->deltaF;
  f0bin  = floor(f0/deltaF +0.5);
  templatePulsar->f0 = f0bin*deltaF;
  closeTemplates->f0[0] = floor(f0/deltaF)*deltaF;
  closeTemplates->f0[1] = ceil(f0/deltaF)*deltaF;
 
  /* sky location, depending if  full sky or small patch is analyzed */
/*
 *   deltaX = deltaF/(params->vTotC * params->pixelFactor *
 *  	           (params->fmin + params->fSearchBand) );
 */
  deltaX = deltaF/(params->vTotC * params->pixelFactor * f0 );
  
  
  if (params->fullSky){ /*full sky*/   
    REAL8 kkcos;
    
    TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
    longitude = randval * LAL_TWOPI;
    TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
    kkcos = 2.0* randval -1.0;
    latitude = acos(kkcos) -LAL_PI_2;
  }
  else {  /*small patch */  
    TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
    longitude = params->alpha + (params->patchSizeAlpha) *(randval-0.5);
    TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
    latitude = params->delta + (params->patchSizeDelta) *(randval-0.5);    
  }
  
  injectPulsar->longitude = longitude;
  injectPulsar->latitude  = latitude;   
  
  {
    REAL8UnitPolarCoor    template, par; 
    REAL8UnitPolarCoor    templRotated;
    REAL8Cart2Coor        templProjected;
    REAL8      dX1[2], dX2[2];
    INT4      ii,jj,kk;
    
    par.alpha = injectPulsar->longitude;
    par.delta = injectPulsar->latitude; 

    /* mismatch with the template in stereographic plane */
    TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
    templProjected.x = dX1[0] = deltaX*(randval-0.5);
    TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
    templProjected.y = dX2[0] = deltaX*(randval-0.5);

    if (dX1[0]<0.0) { 
      dX1[1]= dX1[0]+deltaX;
    } else {
      dX1[1]= dX1[0]-deltaX;
    }
    
    if (dX2[0]<0.0) { 
      dX2[1]= dX2[0]+deltaX;
    } else {
      dX2[1]= dX2[0]-deltaX;
    }
    
    /* invert the stereographic projection for a point on the projected plane */
    TRY( LALStereoInvProjectCart( status->statusPtr,
                                &templRotated, &templProjected ), status );
    /* inverse rotate the mismatch from the south pole to desired location */
    TRY( LALInvRotatePolarU( status->statusPtr, &template, &templRotated, &par), status);
    templatePulsar->longitude = template.alpha; 
    templatePulsar->latitude = template.delta; 
     
    kk=0;
    for (ii=0; ii<2; ii++){
      for (jj=0; jj<2; jj++) {
      templProjected.x = dX1[ii];
      templProjected.y = dX2[jj];
      TRY( LALStereoInvProjectCart( status->statusPtr,
                                &templRotated, &templProjected ), status );
      TRY( LALInvRotatePolarU( status->statusPtr, &(closeTemplates->skytemp[kk]), &templRotated, 
                               &par), status);
      ++kk;
      }
    }
    
  }

  /* now the spindown if any */
  msp = params->spnFmax.length ;
  closeTemplates->f1[0] = 0.0;
  closeTemplates->f1[1] = 0.0;

  ASSERT (templatePulsar->spindown.length == msp, status, DRIVEHOUGHCOLOR_EBAD,
	  DRIVEHOUGHCOLOR_MSGEBAD);
  ASSERT (injectPulsar->spindown.length == msp, status, DRIVEHOUGHCOLOR_EBAD,
	  DRIVEHOUGHCOLOR_MSGEBAD);
  
  if(msp){ /*if there are spin-down values */
    REAL8 deltaFk, spink;
    REAL8 timeObsInv;
    UINT4   i;
    ASSERT (injectPulsar->spindown.data,  status, DRIVEHOUGHCOLOR_ENULL, 
	    DRIVEHOUGHCOLOR_MSGENULL);
    ASSERT (templatePulsar->spindown.data,  status, DRIVEHOUGHCOLOR_ENULL, 
	    DRIVEHOUGHCOLOR_MSGENULL);
    ASSERT (params->spnFmax.data,  status, DRIVEHOUGHCOLOR_ENULL, 
	    DRIVEHOUGHCOLOR_MSGENULL);
    
    /* delta f_k = k! deltaF/ [T_Obs}^k  spd grid resolution*/
    timeObsInv= 1.0/params->timeObs;
    deltaFk= deltaF*timeObsInv;
    
    /* first spin-down parameter, (only spin-down) */	    
    TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
    spink=params->spnFmax.data[0]* randval;
    
    injectPulsar->spindown.data[0]= spink;
    templatePulsar->spindown.data[0] = floor(spink/deltaFk +0.5)*deltaFk;
    
    closeTemplates->f1[0] = floor(spink/deltaFk)*deltaFk;
    closeTemplates->f1[1] = ceil( spink/deltaFk)*deltaFk;

    /* the rest of the spin orders */
    for (i=1; i< msp; ++i) {
      TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
      spink=params->spnFmax.data[i]* (2.0* randval-1.0);
      injectPulsar->spindown.data[i]= spink;   
      deltaFk= deltaFk*timeObsInv*(i+1.0);
      templatePulsar->spindown.data[i] = floor(spink/deltaFk +0.5)*deltaFk;
    }
  }
  /* free memory */
  TRY( LALDestroyRandomParams(status->statusPtr, &randPar), status);
  
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}

 
/***************************************************************************/
void GenerateInjectParamsNoVeto(LALStatus   *status,
                        PulsarData           *injectPulsar,
                        HoughTemplate        *templatePulsar,
			HoughNearTemplates   *closeTemplates,
                        HoughInjectParams    *params ){
			
  INT4          seed=0; /* seed generated using current time */
  REAL4         randval;
  RandomParams  *randPar=NULL;
  FILE     *fpRandom;
  INT4     count;
  
  REAL4    cosiota, h0;
  REAL8    f0, deltaF, deltaX;
  REAL8    latitude, longitude;  /* of the source in radians */
  INT8    f0bin;
  UINT4    msp;
  
  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  
  /*   Make sure the arguments are not NULL: */
  ASSERT (injectPulsar,  status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (templatePulsar, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  ASSERT (params, status, DRIVEHOUGHCOLOR_ENULL, DRIVEHOUGHCOLOR_MSGENULL);
  
  /*  ++++++++++++++++++from makefakedata
   * Modified so as to not create random number parameters with seed
   * drawn from clock.  Seconds don't change fast enough and sft's
   * look alike.  We open /dev/urandom and read a 4 byte integer from
   * it and use that as our seed.  Note: /dev/random is slow after the
   * first, few accesses.
   */

  fpRandom = fopen("/dev/urandom","r");
  ASSERT (fpRandom, status, DRIVEHOUGHCOLOR_EFILE,  DRIVEHOUGHCOLOR_MSGEFILE); 
  
  count = fread(&seed, sizeof(INT4),1, fpRandom);
  if (count != 0) ABORT ( status, DRIVEHOUGHCOLOR_EARG,  DRIVEHOUGHCOLOR_MSGEARG);
  
  fclose(fpRandom);
  
  TRY( LALCreateRandomParams(status->statusPtr, &randPar, seed), status);
  
 /*
  *   to create a single random deviate distributed uniforly between zero and unity		     
  *   TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
  */
  
  
  /* get random value phi0 [0, 2 pi] */ 
  TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
  injectPulsar->phi0 = randval * LAL_TWOPI;
  
  /* get random value cos iota [-1,1] */ 
  TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
  cosiota = 2.0* randval -1.0;
  
  h0=params->h0;
  injectPulsar->aCross = h0*cosiota;
  injectPulsar->aPlus  = 0.5*h0*(1.0 + cosiota*cosiota);
  
  /* get random value psi [0, 2 pi] */ 
  TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
  injectPulsar->psi = randval * LAL_TWOPI;

  /* getting random number for the frequency (and mismatch)*/
  TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
  f0 = params->fmin + (params->fSearchBand) * randval;
   
  injectPulsar->f0 = f0;
  deltaF = params->deltaF;
  f0bin  = floor(f0/deltaF +0.5);
  templatePulsar->f0 = f0bin*deltaF;
  closeTemplates->f0[0] = floor(f0/deltaF)*deltaF;
  closeTemplates->f0[1] = ceil(f0/deltaF)*deltaF;
 
  /* sky location, depending if  full sky or small patch is analyzed */
  deltaX = deltaF/(params->vTotC * params->pixelFactor *
 	           (params->fmin + params->fSearchBand) );
  
  
  if (params->fullSky){ /*full sky*/   
    REAL8 kkcos;
    
    TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
    longitude = randval * LAL_TWOPI;
    TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
    kkcos = 2.0* randval -1.0;
    latitude = acos(kkcos) -LAL_PI_2;
  }
  else {  /*small patch */  
    TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
    longitude = params->alpha + (params->patchSizeAlpha) *(randval-0.5);
    TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
    latitude = params->delta + (params->patchSizeDelta) *(randval-0.5);    
  }
  
  injectPulsar->longitude = longitude;
  injectPulsar->latitude  = latitude;   
  
  {
    REAL8UnitPolarCoor    template, par; 
    REAL8UnitPolarCoor    templRotated;
    REAL8Cart2Coor        templProjected;
    REAL8      dX1[2], dX2[2];
    INT4      ii,jj,kk;
    
    par.alpha = injectPulsar->longitude;
    par.delta = injectPulsar->latitude; 

    /* mismatch with the template in stereographic plane */
    TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
    templProjected.x = dX1[0] = deltaX*(randval-0.5);
    TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
    templProjected.y = dX2[0] = deltaX*(randval-0.5);

    if (dX1[0]<0.0) { 
      dX1[1]= dX1[0]+deltaX;
    } else {
      dX1[1]= dX1[0]-deltaX;
    }
    
    if (dX2[0]<0.0) { 
      dX2[1]= dX2[0]+deltaX;
    } else {
      dX2[1]= dX2[0]-deltaX;
    }
    
    /* invert the stereographic projection for a point on the projected plane */
    TRY( LALStereoInvProjectCart( status->statusPtr,
                                &templRotated, &templProjected ), status );
    /* inverse rotate the mismatch from the south pole to desired location */
    TRY( LALInvRotatePolarU( status->statusPtr, &template, &templRotated, &par), status);
    templatePulsar->longitude = template.alpha; 
    templatePulsar->latitude = template.delta; 
     
    kk=0;
    for (ii=0; ii<2; ii++){
      for (jj=0; jj<2; jj++) {
      templProjected.x = dX1[ii];
      templProjected.y = dX2[jj];
      TRY( LALStereoInvProjectCart( status->statusPtr,
                                &templRotated, &templProjected ), status );
      TRY( LALInvRotatePolarU( status->statusPtr, &(closeTemplates->skytemp[kk]), &templRotated, 
                               &par), status);
      ++kk;
      }
    }
    
  }

  /* now the spindown if any */
  msp = params->spnFmax.length ;
  closeTemplates->f1[0] = 0.0;
  closeTemplates->f1[1] = 0.0;

  ASSERT (templatePulsar->spindown.length == msp, status, DRIVEHOUGHCOLOR_EBAD,
	  DRIVEHOUGHCOLOR_MSGEBAD);
  ASSERT (injectPulsar->spindown.length == msp, status, DRIVEHOUGHCOLOR_EBAD,
	  DRIVEHOUGHCOLOR_MSGEBAD);
  
  if(msp){ /*if there are spin-down values */
    REAL8 deltaFk, spink;
    REAL8 timeObsInv;
    UINT4   i;
    ASSERT (injectPulsar->spindown.data,  status, DRIVEHOUGHCOLOR_ENULL, 
	    DRIVEHOUGHCOLOR_MSGENULL);
    ASSERT (templatePulsar->spindown.data,  status, DRIVEHOUGHCOLOR_ENULL, 
	    DRIVEHOUGHCOLOR_MSGENULL);
    ASSERT (params->spnFmax.data,  status, DRIVEHOUGHCOLOR_ENULL, 
	    DRIVEHOUGHCOLOR_MSGENULL);
    
    /* delta f_k = k! deltaF/ [T_Obs}^k  spd grid resolution*/
    timeObsInv= 1.0/params->timeObs;
    deltaFk= deltaF*timeObsInv;
    
    /* first spin-down parameter, (only spin-down) */	    
    TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
    spink=params->spnFmax.data[0]* randval;
    
    injectPulsar->spindown.data[0]= spink;
    templatePulsar->spindown.data[0] = floor(spink/deltaFk +0.5)*deltaFk;
    
    closeTemplates->f1[0] = floor(spink/deltaFk)*deltaFk;
    closeTemplates->f1[1] = ceil( spink/deltaFk)*deltaFk;

    /* the rest of the spin orders */
    for (i=1; i< msp; ++i) {
      TRY( LALUniformDeviate(status->statusPtr, &randval, randPar), status);
      spink=params->spnFmax.data[i]* (2.0* randval-1.0);
      injectPulsar->spindown.data[i]= spink;   
      deltaFk= deltaFk*timeObsInv*(i+1.0);
      templatePulsar->spindown.data[i] = floor(spink/deltaFk +0.5)*deltaFk;
    }
  }
  /* free memory */
  TRY( LALDestroyRandomParams(status->statusPtr, &randPar), status);
  
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}


/* ****************************************************************/
/* Computing the frequency path f(t) = f0(t)* (1+v/c.n)   */
/* without mismatch                                       */
/* ****************************************************************/   
/******************************************************************/
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
    foft->data[j] = f0new * (1.0 +vcProdn);
  }    
    
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}

/******************************************************************/

			
/* ****************************************************************/
/*    PrintLogFile2 (like in the Driver, but this one doesn't     */
/*    copy the contents of skypatch file)                               */
/* ****************************************************************/   
/******************************************************************/


void PrintLogFile2 (LALStatus       *status,
		   CHAR            *dir,
		   CHAR            *basename,
		   LALStringVector *linefiles,
		   CHAR            *executable )
{
  CHAR *fnameLog=NULL; 
  FILE *fpLog=NULL;
  CHAR *logstr=NULL; 
  UINT4 k;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  
  /* open log file for writing */
  fnameLog = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
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

  fprintf( fpLog, "## LOG FILE FOR MC Inject Hough\n\n");
  fprintf( fpLog, "# User Input:\n");
  fprintf( fpLog, "#-------------------------------------------\n");
  fprintf( fpLog, "%s", logstr);
  LALFree(logstr);

  /* copy contents of linefile if necessary */
  if ( linefiles ) {

    for ( k = 0; k < linefiles->length; k++) {
      
      if ((fpLog = fopen(fnameLog, "a")) != NULL) {
	CHAR command[1024] = "";
	fprintf (fpLog, "\n\n# Contents of linefile %s :\n", linefiles->data[k]);
	fprintf (fpLog, "# -----------------------------------------\n");
	fclose (fpLog);
	sprintf(command, "cat %s >> %s", linefiles->data[k], fnameLog);      
        if ( system(command) ) fprintf (stderr, "\nsystem('%s') returned non-zero status!\n\n", command );
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

