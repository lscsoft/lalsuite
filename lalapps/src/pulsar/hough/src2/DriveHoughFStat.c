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
#define DELTAF 2.778e-05 /* resolution corresponding to 10h */
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
  LIGOTimeGPSVector **midTs=NULL; /* time stamps of mid points of sfts */
  LIGOTimeGPSVector midTstack, startTsft; 
  REAL8VectorSequence *velStack, *posStack; /* velocities and positions at midTstack */
  LALLeapSecFormatAndAcc lsfas = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
  INT4 tmpLeap; 
  REAL8 tObs, tStart, tEnd;
  REAL8  *tStack, tStackAvg; /* duration of each stack */

  /* sft related stuff */
  SFTVector *inputSFTs=NULL;  /* vector of SFTtypes and SFTtype is COMPLEX8FrequencySeries */
  SFTVectorSequence stackSFTs; /* sequence of sft vectors -- one for each stack */
  INT4 *mCohSft, nSFTs; /* number of SFTs in each stack and total number of SFTs */
  INT4 nStacks; /* number of stacks -- not necessarily same as uvar_nStacks! */
  INT4 sftlength; /* number of bins in each sft */
  REAL8 deltaF, timeBase; /* frequency resolution of SFTs */
  INT8 sftFminBin; /* first sft bin index */
  INT4 fSearchBinIni, fSearchBinFin; /* frequency bins of start and end search frequencies */
  REAL8 deltaFstack; /* frequency resolution of Fstat calculation */

  /* LALdemod related stuff */
  REAL8FrequencySeriesVector FstatVect; /* Fstatistic vectors for each stack */
  FstatStackParams FstatPar;
  INT4 binsFstat; /* number of frequency values where Fstat is calculated */

  /* hough variables */
  INT4  nfSizeCylinder=NFSIZE;
  HOUGHPeakGramVector pgV;

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

    tempDir = (CHAR *)LALMalloc(512*sizeof(CHAR));
    strcpy(tempDir, uvar_sftDir);
    strcat(tempDir, "/*SFT*.*"); 
    LAL_CALL( LALReadSFTfiles ( &status, &inputSFTs, uvar_fStart, uvar_fStart + uvar_fBand, nfSizeCylinder + uvar_blocksRngMed , tempDir), &status);

    /* normalize sfts */
    LAL_CALL( LALNormalizeSFTVect (&status, inputSFTs, uvar_blocksRngMed, 0), &status);

    LALFree(tempDir);
  
    /* get start time of sfts */
    startTsft.data = (LIGOTimeGPS *)LALMalloc( nSFTs * sizeof(LIGOTimeGPS));
    for (j=0; j<nSFTs; j++) 
      startTsft.data[j] = inputSFTs->data[j].epoch;
    
    /* set other sft parameters */
    nSFTs = inputSFTs->length;
    sftlength = inputSFTs->data->data->length;
    deltaF = inputSFTs->data->deltaF;
    timeBase = 1.0/deltaF;
    sftFminBin = floor( timeBase * inputSFTs->data->f0 + 0.5);

    /* calculate start and end times and tobs */
    LAL_CALL( LALGPStoFloat ( &status, &tStart, startTsft.data), &status);
    LAL_CALL( LALGPStoFloat ( &status, &tEnd, startTsft.data + nSFTs), &status);
    tEnd += timeBase;
    tObs = tEnd - tStart;

  } /* end of sft reading block */



  /*------------- set up stacks -----------------*/

  if (uvar_nStacks > nSFTs) {
    fprintf(stderr, "invalid number of stacks...exiting\n");
    exit(1);
  }

  /* set up the stacks */
  LAL_CALL( SetUpStacks1( &status, &stackSFTs, inputSFTs, uvar_nStacks), &status);

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
    tempT2 = stackSFTs.data[k].data[tempInt].epoch;
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
  /* starting bin, number of bins and final bin for calculating Fstat */
  fSearchBinIni = floor( tStackAvg * uvar_fStart + 0.5);
  binsFstat = floor( uvar_fBand * tStackAvg + 0.5);
  fSearchBinFin = fSearchBinIni + binsFstat - 1;
    

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
    LAL_CALL (LALDetectorVel ( &status, velStack->data + k*nStacks, midTstack.data + k, detector, edat), &status);
    LAL_CALL (LALDetectorPos ( &status, posStack->data + k*nStacks, midTstack.data + k, detector, edat), &status);
  }



  /*------------- calculate F statistic for each stack --------------*/
  /* memory for Fstatistic Vector */
  FstatVect.length = nStacks;
  FstatVect.data = NULL;
  FstatVect.data = (REAL8FrequencySeries *)LALMalloc(nStacks * sizeof(REAL8FrequencySeries));
  for (k=0; k<nStacks; k++) {
    FstatVect.data[k].epoch = midTstack.data[k];
    FstatVect.data[k].deltaF = 1.0/tStackAvg;
    FstatVect.data[k].f0 = uvar_fStart;
    FstatVect.data[k].data->length = binsFstat;
    FstatVect.data[k].data->data = (REAL8 *)LALMalloc( binsFstat * sizeof(REAL8));
  }
  

  /* now calculate F-statistic */


  /*------------ select peaks ------------*/ 
  /* first allocate some memory */
  pgV.length = nStacks;
  pgV.pg = (HOUGHPeakGram *)LALMalloc( nStacks * sizeof(HOUGHPeakGram));
  for (k=0; k<nStacks; k++) {
    pgV.pg[k].deltaF = deltaF;
    pgV.pg[k].fBinIni = fSearchBinIni;
    pgV.pg[k].fBinFin = fSearchBinFin;
  }

  /* compute the peakgrams */
  LAL_CALL( FstatVectToPeakGram( &status, &pgV, &FstatVect, uvar_FstatThr), &status);


  /*************** calculate Hough map *************/
  /* set up the Hough parameters */


  /*********** free all remaining memory **********/
  /* free timestamp and Vel/Pos vectors */
  for (k=0; k<nStacks; k++) {
      LALFree(midTs[k]->data);
      LALFree(midTs[k]);
    }  
  LALFree(midTs);
  LALFree(startTsft.data);
  LALFree(midTstack.data);
  LALFree(mCohSft);
  LALFree(tStack);

  LAL_CALL( LALDDestroyVectorSequence (&status,  &velStack), &status);
  LAL_CALL( LALDDestroyVectorSequence (&status,  &posStack), &status);

  for(k=0; k<nStacks; k++)
    LALFree(FstatVect.data[k].data->data);
  LALFree(FstatVect.data);

  LALFree(stackSFTs.data);

  /* free ephemeris */
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);

  /* free peakgrams */
  for (k=0; k<nStacks; k++) 
    LALFree(pgV.pg[k].peak);
  LALFree(pgV.pg);

  LAL_CALL (LALDestroySFTVector(&status, &inputSFTs),&status );
  LAL_CALL (LALDDestroyVector (&status, &(FstatPar.spindown)), &status);
  LAL_CALL (LALDestroyUserVars(&status), &status);  

  LALCheckMemoryLeaks();

  return 0;
}


void ComputeFstatStack (LALStatus *status, 
			REAL8FrequencySeriesVector *Fstat, 
			SFTVectorSequence *inputSFTs, 
			FstatStackParams *params)
{
  Fcomponents FaFb;
  SSBtimes *tSSB;
  AMCoeffs *amcoe;
  SkyPosition skyPoint;

  INT4 k, nStacks;

  INITSTATUS( status, "ComputeFstatStack", rcsid );
  ATTATCHSTATUSPTR (status);

  skyPoint.longitude = params->alpha;
  skyPoint.latitude = params->delta;
  skyPoint.system = COORDINATESYSTEM_EQUATORIAL;
  LAL_CALL (LALNormalizeSkyPosition( status, &skyPoint, &skyPoint), status);

  nStacks = params->nStacks;
  for(k=0; k<nStacks; k++) {
    
  }
  
  DETATCHSTATUSPTR (status);
  RETURN(status); 
}





void ComputeFstatHoughMap(LALStatus *status,
			  HOUGHMapTotal   *ht,   /* the total Hough map */
			  HOUGHPeakGramVector *pgV,
			  HoughParams *params)
{

  INITSTATUS( status, "ComputeFstatHoughMap", rcsid );
  ATTATCHSTATUSPTR (status);



  DETATCHSTATUSPTR (status);
  RETURN(status);

}


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

  nStacks = pgV->length;
  nSearchBins = FstatVect->data->data->length;

  upg = (UCHAR *)LALMalloc( nSearchBins * sizeof(UCHAR));

  /* loop over each stack and set peakgram */
  for (k=0; k<nStacks; k++) {
    INT4 *pInt; /* temporary pointer */
    REAL8 *pV;  /* temporary pointer */

    pV = FstatVect->data[k].data->data;

    /* loop over stack, count peaks, and set upg values */
    nPeaks = 0;
    for(j=0; j<nSearchBins; j++) {
      if ( pV[j] > thr ) {
	nPeaks++;	
	upg[j] = 1; 
      }
      else
	upg[j] = 0;
    }

    pgV->pg[k].length = nPeaks; 
    pgV->pg[k].peak = (INT4 *)LALMalloc( nPeaks * sizeof(INT4)); 
    pInt = pgV->pg[k].peak;
    /* do loop again and fill peakgram vector */
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
		 INT4 nStacks)
{

  INITSTATUS( status, "SetUpStacks2", rcsid );
  ATTATCHSTATUSPTR (status);

  

  DETATCHSTATUSPTR (status);
  RETURN(status);
}
