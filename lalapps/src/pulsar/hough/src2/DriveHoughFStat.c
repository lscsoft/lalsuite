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
 * \brief
 * Puts together the F-statistic and Hough routines.  It calculates F-stat values
 * for different time segments and combines them semi-coherently using the
 * Hough transform
 *                                                                          
 ****/



#include"./DriveHoughFStat.h"

RCSID( "$Id$");

#define TRUE (1==1)
#define FALSE (1==0)


extern int lalDebugLevel;

/* default values for input variables */
#define EARTHEPHEMERIS "../src/earth00-04.dat"
#define SUNEPHEMERIS "../src/sun00-04.dat"
#define NSTACKS 10
#define IFO 2         /*  detector, 1:GEO, 2:LLO, 3:LHO */
#define BLOCKSRNGMED 101 /* Running median window size */
#define FSTART 255.0       
#define FBAND 0.001
#define FDOT 0.0
#define ALPHA 1.57
#define DELTA  0.0
#define NFSIZE  21
#define DTERMS 8
#define FSTATTHRESHOLD 2.6
#define SFTDIRECTORY "/local_data/badkri/fakesfts/"
#define FNAMEOUT "./OutHoughFStat"

int main( int argc, char *argv[]) {
  LALStatus status = blank_status;	/* initialize status */
  
  INT4 j,k; /* temp loop variables: k loops over stacks and j over SFTs in a stack*/

  /* detector, ephemeris and velocity/position vector */
  LALDetector detector;
  EphemerisData *edat = NULL;
  LIGOTimeGPSVector midTstack, startTsft; 
  REAL8VectorSequence *velStack=NULL, *posStack=NULL; /* velocities and positions at midTstack */
  LALLeapSecFormatAndAcc lsfas = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
  INT4 tmpLeap; 
  REAL8 tObs, tStart, tEnd;
  REAL8  *tStack, tStackAvg; /* duration of each stack */
  REAL8 refTime;
  REAL8 extra, fStartExtra, fBandExtra; /* extra wings for Fstat calculation */

  /* sft related stuff */
  SFTVector *inputSFTs=NULL;  /* vector of SFTtypes and SFTtype is COMPLEX8FrequencySeries */
  SFTVectorSequence stackSFTs; /* sequence of sft vectors -- one for each stack */
  INT4 *mCohSft, nSFTs; /* number of SFTs in each stack and total number of SFTs */
  INT4 nStacks; /* number of stacks -- not necessarily same as uvar_nStacks! */
  INT4 sftlength; /* number of bins in each sft */
  REAL8 deltaF, timeBase; /* frequency resolution of SFTs */
  INT8 sftFminBin; /* first sft bin index */
  INT8 fBinIni, fBinFin, binsHough; /* frequency bins of start and end search frequencies */
  INT8 fSearchBinIni, fSearchBinFin, binsFstat; /* same as above for Fstat */
  REAL8 deltaFstack; /* frequency resolution of Fstat calculation */

  /* LALdemod related stuff */
  REAL8FrequencySeriesVector FstatVect; /* Fstatistic vectors for each stack */
  FstatStackParams FstatPar;

  /* hough variables */
  INT4  nfSizeCylinder=NFSIZE;
  HOUGHPeakGramVector pgV;
  HoughParams houghPar;
  HOUGHMapTotal ht;

  /* variables for logging */
  CHAR *fnamelog=NULL;
  FILE *fpLog = NULL;
  CHAR *logstr=NULL; 

  /* user variables */
  BOOLEAN uvar_help; /* true if -h option is given */
  REAL8 uvar_alpha, uvar_delta;  /* sky-location angles */
  REAL8 uvar_fdot; /* first spindown value */
  REAL8 uvar_fStart, uvar_fBand;
  REAL8 uvar_FstatThr; /* threshold of Fstat to select peaks */
  INT4 uvar_ifo, uvar_blocksRngMed, uvar_nStacks, uvar_Dterms;
  REAL8 uvar_refTime;
  INT4 uvar_SSBprecision = SSBPREC_RELATIVISTIC;
  CHAR *uvar_earthEphemeris=NULL;
  CHAR *uvar_sunEphemeris=NULL;
  CHAR *uvar_sftDir=NULL;
  CHAR *uvar_fnameout=NULL;

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /************/
  /* set defaults, read user variables, log user variables and log cvs tags */
  /* LALDebugLevel must be called before anything else */
  lalDebugLevel = 1;
  LAL_CALL( LALGetDebugLevel( &status, argc, argv, 'd'), &status);

  /* now set the defaults */
  uvar_help = FALSE;
  uvar_nStacks = NSTACKS;
  uvar_Dterms = DTERMS;
  uvar_alpha = ALPHA;
  uvar_delta = DELTA;
  uvar_fdot = FDOT;
  uvar_fStart = FSTART;
  uvar_fBand = FBAND;
  uvar_ifo = IFO;
  uvar_blocksRngMed = BLOCKSRNGMED;
  uvar_FstatThr = FSTATTHRESHOLD;
  uvar_earthEphemeris = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_earthEphemeris,EARTHEPHEMERIS);

  uvar_sunEphemeris = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_sunEphemeris,SUNEPHEMERIS);

  uvar_sftDir = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_sftDir,SFTDIRECTORY);

  uvar_fnameout = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_fnameout,FNAMEOUT);

  /* register user input variables */
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "help",            'h', UVAR_HELP,     "Print this message",            &uvar_help),            &status);  
  LAL_CALL( LALRegisterINTUserVar(    &status, "ifo",             'i', UVAR_OPTIONAL, "Detector GEO(1) LLO(2) LHO(3)", &uvar_ifo ),            &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "nStacks",         'N', UVAR_OPTIONAL, "Number of stacks",              &uvar_nStacks ),        &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "Dterms",          'n', UVAR_OPTIONAL, "For Dirichlet Kernel approx.",  &uvar_Dterms ),         &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed",    'w', UVAR_OPTIONAL, "RngMed block size",             &uvar_blocksRngMed),    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "earthEphemeris",  'E', UVAR_OPTIONAL, "Earth Ephemeris file",          &uvar_earthEphemeris),  &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sunEphemeris",    'S', UVAR_OPTIONAL, "Sun Ephemeris file",            &uvar_sunEphemeris),    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftDir",          'D', UVAR_OPTIONAL, "SFT Directory",                 &uvar_sftDir),          &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "fnameout",        'o', UVAR_OPTIONAL, "Output basefileneme",           &uvar_fnameout),        &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "alpha",           'r', UVAR_OPTIONAL, "Right ascension",               &uvar_alpha),           &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "delta",           'l', UVAR_OPTIONAL, "Declination",                   &uvar_delta),           &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fStart",          'f', UVAR_OPTIONAL, "Start search frequency",        &uvar_fStart),          &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fdot",            'p', UVAR_OPTIONAL, "Spindown parameter",            &uvar_fdot),            &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "FstatThr",        't', UVAR_OPTIONAL, "Threshold on Fstatistic",       &uvar_FstatThr),        &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "refTime",          0,  UVAR_OPTIONAL, "Reference time for pulsar par", &uvar_refTime),         &status);
  LAL_CALL( LALRegisterINTUserVar (   &status, "SSBprecision",	   0,  UVAR_DEVELOPER,"Precision for SSB transform.",  &uvar_SSBprecision),    &status);

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
  if (uvar_ifo ==1) detector=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if (uvar_ifo ==2) detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
  if (uvar_ifo ==3) detector=lalCachedDetectors[LALDetectorIndexLHODIFF];

  /* write the log file */
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


  /*------------- read sfts and set up sft timestamp vector ----------*/
  {
    CHAR *tempDir;
    REAL8 doppWings, fmin, fmax;

    tempDir = (CHAR *)LALMalloc(512*sizeof(CHAR));
    strcpy(tempDir, uvar_sftDir);
    strcat(tempDir, "/*SFT*.*");

    doppWings = (uvar_fStart + uvar_fBand) * VTOT;    
    extra = doppWings;
    fStartExtra = uvar_fStart - extra; /* read in a little bit extra ?*/
    fBandExtra = uvar_fBand + 2.0*extra;

    fmin = fStartExtra - doppWings;
    fmax = fStartExtra + fBandExtra + doppWings;
    LAL_CALL( LALReadSFTfiles ( &status, &inputSFTs, fmin, fmax, 
				nfSizeCylinder + uvar_blocksRngMed + uvar_Dterms, 
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

  /* number of bins for calculating Fstat */
  binsFstat = floor( fBandExtra * tStackAvg + 0.5);
  fSearchBinIni = floor( tStackAvg * fStartExtra + 0.5);
  fSearchBinFin = fSearchBinIni + binsFstat;

  /* start and end bin for calculating hough map */
  fBinIni = floor( tStackAvg * uvar_fStart + 0.5);
  binsHough = floor( tStackAvg * uvar_fBand + 0.5);
  fBinFin = fBinIni + binsHough - 1;
    

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
  /* memory for Fstatistic Vector */
  FstatVect.length = nStacks;
  FstatVect.data = NULL;
  FstatVect.data = (REAL8FrequencySeries *)LALMalloc(nStacks * sizeof(REAL8FrequencySeries));
  for (k=0; k<nStacks; k++) {
    FstatVect.data[k].epoch = midTstack.data[k];
    FstatVect.data[k].deltaF = deltaFstack;
    FstatVect.data[k].f0 = fStartExtra;
    FstatVect.data[k].data = (REAL8Sequence *)LALMalloc(sizeof(REAL8Sequence));
    FstatVect.data[k].data->length = binsFstat;
    FstatVect.data[k].data->data = (REAL8 *)LALMalloc( binsFstat * sizeof(REAL8));
  }
  

  /* set up parameters for Fstat calculation */
  FstatPar.mCohSft = mCohSft;
  FstatPar.timeBase = timeBase;
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

  /* calculate the Fstatistic */
  LAL_CALL(ComputeFstatStack( &status, &FstatVect, &stackSFTs, &FstatPar), &status);



  /*------------ select peaks ------------*/ 
  /* first allocate memory for peakgrams */
  pgV.length = nStacks;
  pgV.pg = (HOUGHPeakGram *)LALMalloc( nStacks * sizeof(HOUGHPeakGram));
  for (k=0; k<nStacks; k++) {
    pgV.pg[k].deltaF = deltaF;
    pgV.pg[k].fBinIni = fSearchBinIni;
    pgV.pg[k].fBinFin = fSearchBinFin;
  }

  /* compute the peakgrams */
  LAL_CALL( FstatVectToPeakGram( &status, &pgV, &FstatVect, uvar_FstatThr), &status);

  LAL_CALL( PrintFstat ( &status, FstatVect.data, uvar_fnameout), &status);


  /*--------------- calculate Hough map ---------------*/
  /* set up the Hough parameters */
  houghPar.tStart = tStart;
  houghPar.fBinIni = fBinIni;
  houghPar.fBinFin = fBinFin;
  houghPar.nfSizeCylinder = nfSizeCylinder;
  houghPar.detector = detector;
  houghPar.ts = &midTstack;
  houghPar.vel = velStack;
  houghPar.pos = posStack;
  houghPar.alpha = uvar_alpha;
  houghPar.delta = uvar_delta;
  houghPar.fdot = NULL;
  LAL_CALL ( LALDCreateVector( &status, &(houghPar.fdot), 1), &status);
  houghPar.fdot->data[0] = uvar_fdot;

  LAL_CALL ( ComputeFstatHoughMap ( &status, &ht, &pgV, &houghPar), &status);

  /*------------ free all remaining memory -----------*/


  /* these can be moved to after Fstat calculation */
  for(k=0; k<nStacks; k++) {
    LALFree(FstatVect.data[k].data->data);
    LALFree(FstatVect.data[k].data);
  }
  LALFree(FstatVect.data);
  LAL_CALL (LALDestroySFTVector(&status, &inputSFTs),&status );
  LAL_CALL (LALDDestroyVector (&status, &(FstatPar.fdot)), &status);
  LAL_CALL (LALDDestroyVector (&status, &(houghPar.fdot)), &status);
  LALFree(stackSFTs.data);
  LALFree(startTsft.data);
  LALFree(mCohSft);

  /* these can be moved to after hough map calculation */
  /* free peakgrams */
  for (k=0; k<nStacks; k++) 
    LALFree(pgV.pg[k].peak);
  LALFree(pgV.pg);
  /* free timestamp and Vel/Pos vectors */
  LALFree(midTstack.data);
  LALFree(tStack);
  LAL_CALL( LALDDestroyVectorSequence (&status,  &velStack), &status);
  LAL_CALL( LALDDestroyVectorSequence (&status,  &posStack), &status);

  /* free ephemeris */
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);

  LAL_CALL (LALDestroyUserVars(&status), &status);  

  LALCheckMemoryLeaks();

  return 0;
}


void ComputeFstatStack (LALStatus *status, 
			REAL8FrequencySeriesVector *out, 
			const SFTVectorSequence *stackSFTs, 
			FstatStackParams *params)
{
  /* stuff copied from params */
  REAL8 timeBase = params->timeBase; /* can also be copied from SFTs */
  REAL8 refTime = params->refTime;
  INT4 *mCohSft = params->mCohSft;
  REAL8Vector *fdot = params->fdot;

  /* stuff copied from output Fstat vector */
  INT4 binsFstat = out->data->data->length;
  INT4 nStacks = out->length;
  REAL8 deltaF = out->data->deltaF;
  REAL8 fStart = out->data->f0;

  /* other variables */
  SSBtimes tSSB;
  AMCoeffs amcoe;
  SkyPosition skyPoint;
  DetectorStateSeries *DetectorStates=NULL;
  LIGOTimeGPSVector timeStack;
  INT4 k, j, indexSft;
  REAL8Vector *fkdot=NULL;
  LIGOTimeGPS tempRef;

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

  /* start loop over stacks */
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




void ComputeFstatHoughMap(LALStatus *status,
			  HOUGHMapTotal   *ht,   /* the total Hough map */
			  const HOUGHPeakGramVector *pgV, /* peakgram vector */
			  HoughParams *params)
{

  /* hough structures */
  static HOUGHptfLUTVector   lutV; /* the Look Up Table vector*/
  static PHMDVectorSequence  phmdVS;  /* the partial Hough map derivatives */
  static UINT8FrequencyIndexVector freqInd; /* for trajectory in time-freq plane */
  static HOUGHResolutionPar parRes;   /* patch grid information */
  static HOUGHPatchGrid  patch;   /* Patch description */ 
  static HOUGHParamPLUT  parLut;  /* parameters needed to build lut  */
  static HOUGHDemodPar   parDem;  /* demodulation parameters */
  static HOUGHSizePar    parSize; 

  UINT2  xSide, ySide, maxNBins, maxNBorders;
  INT8  fBinIni, fBinFin, fBin;
  INT4  k, iHmap, nSpin1Max, nStacks, nfSizeCylinder;
  REAL8 deltaF, alpha, delta;
  REAL8 patchSizeX, patchSizeY, f1jump, tStart;
  REAL8VectorSequence *vel, *pos;
  REAL8Vector *fdot;
  LIGOTimeGPSVector   *ts;
  REAL8Vector timeDiffV;

  INITSTATUS( status, "ComputeFstatHoughMap", rcsid );
  ATTATCHSTATUSPTR (status);

  /* copy some params to local variables */
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

  /* first memory allocation */
  lutV.length = nStacks;
  lutV.lut = NULL;
  lutV.lut = (HOUGHptfLUT *)LALMalloc(nStacks*sizeof(HOUGHptfLUT));
  
  phmdVS.length  = nStacks;
  phmdVS.nfSize  = nfSizeCylinder;
  phmdVS.deltaF  = deltaF;
  phmdVS.phmd = NULL;
  phmdVS.phmd=(HOUGHphmd *)LALMalloc(nStacks*nfSizeCylinder*sizeof(HOUGHphmd));
  
  freqInd.deltaF = deltaF;
  freqInd.length = nStacks;
  freqInd.data = NULL;
  freqInd.data =  ( UINT8 *)LALMalloc(nStacks*sizeof(UINT8));
   
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

  /* initialization */  
  fBin= fBinIni;
  iHmap = 0;

  /* calculate time differences from start of observation time */
  timeDiffV.length = nStacks;
  timeDiffV.data = (REAL8 *)LALMalloc( nStacks * sizeof(REAL8));
  for (k=0; k<nStacks; k++) {
    REAL8 tMidStack;

    TRY ( LALGPStoFloat ( status->statusPtr, &tMidStack, ts->data + k), status);
    timeDiffV.data[k] = tMidStack - tStart;
  }

  /* if there are spindowns */
  nSpin1Max = floor(nfSizeCylinder/2.0); /* max number of spindowns */
  f1jump = 1.0 / timeDiffV.data[nStacks - 1]; /* resolution in fdot */
  
  /* start main Hough calculation */
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
    ht->xSide = xSide;
    ht->ySide = ySide;
    ht->mObsCoh = nStacks;
    ht->deltaF = deltaF;
    ht->map   = NULL;
    ht->map   = (HoughTT *)LALMalloc(xSide*ySide*sizeof(HoughTT));
    TRY( LALHOUGHInitializeHT( status->statusPtr, ht, &patch), status); /*not needed */
    

    /*  Search frequency interval possible using the same LUTs */
    fBinSearch = fBin;
    fBinSearchMax = fBin + parSize.nFreqValid - 1 - floor((nfSizeCylinder - 1 )/2.0);
     
    /* Study all possible frequencies with one set of LUT */    
    while ( (fBinSearch <= fBinFin) && (fBinSearch < fBinSearchMax) )  { 
      

      /* Case: No spin-down. Study the fBinSearch */
      ht->f0Bin = fBinSearch;
      ht->spinRes.length =0;
      ht->spinRes.data = NULL;
      for (j=0; j < (UINT4)nStacks; j++)
	freqInd.data[j]= fBinSearch; 
      TRY( LALHOUGHConstructHMT( status->statusPtr, ht, &freqInd, &phmdVS ), status );
            
      ++iHmap;
            

      /* Case: 1 spin-down. at  fBinSearch */
      {
	INT4   n;
	REAL8  f1dis;
	
	ht->spinRes.length = 1;
	ht->spinRes.data = NULL;
	ht->spinRes.data = (REAL8 *)LALMalloc(ht->spinRes.length*sizeof(REAL8));
	    
	for( n=1; n<= nSpin1Max; ++n){ /*loop over all values of f1 */
	  f1dis = - n*f1jump;
	  ht->spinRes.data[0] =  f1dis*deltaF;
	  
	  for (j=0; j < (UINT4)nStacks; j++)
	    freqInd.data[j] = fBinSearch +floor(timeDiffV.data[j]*f1dis+0.5);
	  
	  TRY( LALHOUGHConstructHMT(status->statusPtr, ht, &freqInd, &phmdVS),status );
	  	  
	  
	  ++iHmap;
	  
	} /* end loop over nSpin1Max */
	LALFree(ht->spinRes.data);
      }
      

      /*------ shift the search freq. & PHMD structure 1 freq.bin -------*/
      ++fBinSearch;
      TRY( LALHOUGHupdateSpacePHMDup(status->statusPtr, &phmdVS, pgV, &lutV), status );
      
    }   /* closing second while  */
    
    fBin = fBinSearch;
    
    /*--------------  Free partial memory -----------------*/
    LALFree(patch.xCoor);
    LALFree(patch.yCoor);
    LALFree(ht->map);

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
  LALFree(lutV.lut);
  LALFree(phmdVS.phmd);
  LALFree(freqInd.data);
  LALFree(timeDiffV.data);
  DETATCHSTATUSPTR (status);
  RETURN(status);

}


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

  nStacks = pgV->length;
  nSearchBins = FstatVect->data->data->length;

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


/* given a sftVector, this function splits it up into several different sft vectors */
/* there are basically two ways of doing this: either each stack contains the same number
   of SFTs, or each stack spans the same duration. These two methods are equivalent
   only if there are no gaps between the SFTs */
/* The function SetUpStacks1 distributes the SFT equally among the stacks while 
   SetUpStacks2 makes each stack span the same time duration. */
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



void PrintFstat( LALStatus *status,
		 REAL8FrequencySeries *Fstat, 
		 CHAR *fname)
{

  FILE *fp=NULL;
  INT4 k, length;
  REAL8 freq, deltaF;

  INITSTATUS( status, "PrintFstat", rcsid );
  ATTATCHSTATUSPTR (status);

  length = Fstat->data->length;
  freq = Fstat->f0;
  deltaF = Fstat->deltaF;
  
  fp = fopen(fname, "w");

  for (k=0; k<length; k++) {
    fprintf(fp, "%g   %g\n", freq, Fstat->data->data[k]);
    freq += deltaF;
  }

  fclose(fp);

  DETATCHSTATUSPTR (status);
  RETURN(status);
}
