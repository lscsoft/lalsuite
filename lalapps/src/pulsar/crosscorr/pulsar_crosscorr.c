/*
 *  Copyright (C) 2007 Badri Krishnan
 *  Copyright (C) 2008 Christine Chung, Badri Krishnan and John Whelan
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
 * \author Christine Chung, Badri Krishnan, John Whelan
 * \date 2008
 * \file
 * \brief Perform CW cross-correlation search
 *
 * $Id: pulsar_crosscorr.c,v 1.23 2009/03/13 00:43:04 cchung Exp $
 *


/ lalapps includes */
#include <lalapps.h>
#include <pulsar_crosscorr.h>
#include <lal/PulsarCrossCorr.h>
#include <lal/DopplerScan.h>
#include <gsl/gsl_permutation.h>

RCSID( "$Id: pulsar_crosscorr.c,v 1.23 2009/03/13 00:43:04 cchung Exp $");

/* globals, constants and defaults */

extern int lalDebugLevel;

/* user input variables */
BOOLEAN  uvar_help;
BOOLEAN  uvar_averagePsi;
BOOLEAN  uvar_averageIota;
BOOLEAN  uvar_autoCorrelate;
BOOLEAN  uvar_spindownParams;
INT4     uvar_blocksRngMed; 
INT4     uvar_detChoice;
REAL8    uvar_startTime, uvar_endTime;
REAL8    uvar_f0, uvar_fdot, uvar_fBand, uvar_fdotBand;
REAL8    uvar_fddot, uvar_fddotBand;
REAL8    uvar_fResolution, uvar_fdotResolution, uvar_fddotResolution;
REAL8    uvar_dAlpha, uvar_dDelta; /* resolution for isotropic sky-grid */
REAL8    uvar_maxlag;
REAL8    uvar_psi;
REAL8    uvar_refTime;
REAL8    uvar_cosi;
REAL8    uvar_epsilon;
REAL8    uvar_magfield;
REAL8    uvar_brakingindex;
REAL8    uvar_epsilonBand;
REAL8    uvar_epsilonResolution;
REAL8    uvar_magfieldBand;
REAL8    uvar_magfieldResolution;
REAL8    uvar_brakingindexBand;
REAL8    uvar_brakingindexResolution;


CHAR     *uvar_ephemDir=NULL;
CHAR     *uvar_ephemYear=NULL;
CHAR     *uvar_sftDir=NULL;
CHAR     *uvar_dirnameOut=NULL;
CHAR     *uvar_skyfile=NULL;
CHAR     *uvar_skyRegion=NULL;
CHAR     *uvar_filenameOut=NULL;

#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
#define EPHEM_YEARS "05-09"

#define F0 100
#define FBAND 1
#define LAL_INT4_MAX 2147483647


#define BLOCKSRNGMED 51
#define MAXFILENAMELENGTH 512 /* maximum # of characters  of a filename */

#define DIROUT "./output/"   /* output directory */
#define FILEOUT "CrossCorr_out.dat"
#define BASENAMEOUT "radio"    /* prefix file output */

#define SKYFILE "./skypatchfile"
#define SKYREGION "allsky"

#define TRUE (1==1)
#define FALSE (1==0)

#define SQUARE(x) (x*x)
#define CUBE(x) (x*x*x)

#define N_SPINDOWN_DERIVS 6

void initUserVars (LALStatus *status);


int main(int argc, char *argv[]){
  /* LALStatus pointer */
  static LALStatus status;  

  /* sft related variables */ 

  SFTVector *inputSFTs = NULL;
  SFTtype *tmpSFT = NULL; 
  REAL8FrequencySeries *tmpPSD = NULL;

  SFTListElement *sftList, *sftHead = NULL, *sftTail = NULL;
  PSDListElement *psdList, *psdHead = NULL, *psdTail = NULL;
  REAL8ListElement *freqList, *freqHead = NULL, *freqTail = NULL;
  REAL8ListElement *phaseList, *phaseHead = NULL, *phaseTail = NULL;
  REAL8 deltaF_SFT, timeBase;
  REAL8  *psi = NULL; 
  UINT4 counter = 0, i=0; 
  COMPLEX8FrequencySeries *sft1 = NULL, *sft2 = NULL;
  INT4 sameDet;
  REAL8FrequencySeries *psd1, *psd2;
  COMPLEX16Vector *yalpha = NULL, *ualpha = NULL;

  CrossCorrAmps amplitudes;
  CrossCorrBeamFn *beamfns1, *beamfns2;
  CrossCorrBeamFnListElement *beamList, *beamHead = NULL, *beamTail = NULL;


  /* information about all the ifos */
  DetChoice detChoice;
  LIGOTimeGPS firstTimeStamp, lastTimeStamp;
  REAL8 tObs;
  REAL8 patchSizeX, patchSizeY;

  /* ephemeris */
  EphemerisData    *edat=NULL;
  INT4 tmpLeap;
  LALLeapSecFormatAndAcc  lsfas = {LALLEAPSEC_GPSUTC, LALLEAPSEC_LOOSE};
  CHAR EphemEarth[MAXFILENAMELENGTH];
  CHAR EphemSun[MAXFILENAMELENGTH];

  /* skypatch info */
  REAL8  *skyAlpha=NULL, *skyDelta=NULL,
    *skySizeAlpha=NULL, *skySizeDelta=NULL; 
  INT4   nSkyPatches, skyCounter=0; 
  SkyPatchesInfo skyInfo;

  /* frequency loop info */
  INT4 nfreqLoops=1, freqCounter = 0;
  INT4 nParams = 0;
  REAL8 f_current = 0.0;
  INT4 ualphacounter=0.0;

  /* frequency derivative loop info. we can go up to f_doubledot */
  INT4 nfdotLoops = 1, fdotCounter = 0;
  INT4 nfddotLoops = 1, fddotCounter = 0;
  REAL8 fdot_current = 0.0, delta_fdot = 0.0;
  REAL8 fddot_current = 0.0, delta_fddot = 0.0; 

  /* frequency derivative array, if we search over q1, q2, n */
  INT4 nq1Loops = 1, nq2Loops = 1, nnLoops = 1;
  INT4 q1Counter = 0, q2Counter = 0, nCounter = 0;
  REAL8 eps_current = 0.0, mag_current = 0.0, n_current = 0.0;
  REAL8 delta_eps = 0.0, delta_mag = 0.0, delta_n = 0.0;
  REAL8Vector *fdots = NULL;

  static INT4VectorSequence  *sftPairIndexList=NULL;
  REAL8Vector  *sigmasq;
  PulsarDopplerParams thisPoint;
  static REAL8Vector *rho, *stddev;
  REAL8 tmpstat, freq1, phase1, freq2, phase2;


  REAL8 doppWings, fMin, fMax;


  /* output file */
  FILE *fp = NULL;
  CHAR filename[MAXFILENAMELENGTH];

  /* sft constraint variables */
  LIGOTimeGPS startTimeGPS, endTimeGPS, refTime;

  FILE	   *skytest=NULL;

  SkyPosition skypos; 

  /* new SFT I/O data types */
  SFTCatalog *catalog = NULL;
  SFTCatalog *slidingcat = NULL;
  static SFTConstraints constraints;
  INT4 sftcounter = 0, slidingcounter = 1, listLength = 0;

  /*  MultiSFTVector *multiSFTs = NULL;   */

  /* LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  lalDebugLevel = 0;  /* LALDebugLevel must be called before anything else */
  LAL_CALL( LALGetDebugLevel( &status, argc, argv, 'd'), &status);
   
  LAL_CALL (initUserVars(&status), &status);

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

  if ( uvar_fBand < 0 ) {
    fprintf(stderr, "search frequency band must be positive\n");
    exit(1);
  }

  if ( LALUserVarWasSet(&uvar_fResolution) && (uvar_fResolution <= 0)) {
    fprintf(stderr, "search frequency resolution must be > 0\n");
    exit(1);
  }

  if (uvar_fdotBand < 0) {
    fprintf(stderr, "frequency derivative band must be positive\n");
    exit(1);
  }

  if ( LALUserVarWasSet(&uvar_fdotResolution) && (uvar_fdotResolution <= 0)) {
    fprintf(stderr, "frequency derivative resolution must be > 0\n");
    exit(1);
  }

  if (uvar_fddotBand < 0) {
    fprintf(stderr, "frequency double derivative band must be positive\n");
    exit(1);
  }

  if ( LALUserVarWasSet(&uvar_fddotResolution) && (uvar_fddotResolution <= 0)) {
    fprintf(stderr, "frequency double derivative resolution must be > 0\n");
    exit(1);
  }

  if (uvar_averagePsi && LALUserVarWasSet (&uvar_psi)) {
    fprintf(stderr, "if uvar_averagePsi = TRUE, psi will not be used\n");
    exit(1);
  } 

  if (uvar_averageIota && LALUserVarWasSet (&uvar_cosi)) {
    fprintf(stderr, "if uvar_averageIota = TRUE, cosi will not be used\n");
    exit(1);
  } 

  if (uvar_autoCorrelate) {
    fprintf(stderr, "autoCorrelate option is not currently usable\n");
    exit(1);
  }


  /* open output file */
  strcpy (filename, uvar_dirnameOut);
  strcat(filename, uvar_filenameOut);

  if ((fp = fopen(filename, "w")) == NULL) {
    fprintf(stderr, "error opening output file\n");
    exit(1);
  }

  fprintf(fp, "##Alpha\tDelta\tFrequency\t Fdot \t Fddot \t Normalised Power\n");


  /* set sft catalog constraints */
  constraints.detector = NULL;
  constraints.timestamps = NULL;
  constraints.startTime = NULL;
  constraints.endTime = NULL;

  XLALGPSSet(&startTimeGPS, 0, 0);
  XLALGPSSet(&endTimeGPS, LAL_INT4_MAX, 0); 

  if ( LALUserVarWasSet( &uvar_startTime ) ) {
    XLALGPSSetREAL8(&startTimeGPS, uvar_startTime);

    constraints.startTime = &startTimeGPS;
  }

  if ( LALUserVarWasSet( &uvar_endTime ) ) {
    XLALGPSSetREAL8(&endTimeGPS, uvar_endTime);
    constraints.endTime = &endTimeGPS;
  }


  /* get sft catalog */
  /* note that this code depends very heavily on the fact that the catalog
     returned by LALSFTdataFind is time sorted */
  LAL_CALL( LALSFTdataFind( &status, &catalog, uvar_sftDir, &constraints),
	    &status);
  if ( (catalog == NULL) || (catalog->length == 0) ) {
    fprintf (stderr, "Unable to match any SFTs with pattern '%s'\n",
	     uvar_sftDir );
    exit(1);
  }


  /* get SFT parameters so that we can initialise search frequency resolutions */
  /* calculate deltaF_SFT */
  deltaF_SFT = catalog->data[0].header.deltaF;  /* frequency resolution */
  timeBase= 1.0/deltaF_SFT; /* sft baseline */

  /* catalog is ordered in time so we can get start, end time and tObs */
  firstTimeStamp = catalog->data[0].header.epoch;
  lastTimeStamp = catalog->data[catalog->length - 1].header.epoch;
  tObs = XLALGPSDiff( &lastTimeStamp, &firstTimeStamp ) + timeBase;

  /*set pulsar reference time */
  if (LALUserVarWasSet ( &uvar_refTime )) {
    XLALGPSSetREAL8(&refTime, uvar_refTime);
  } 
  else {	/*if refTime is not set, set it to midpoint of sfts*/
    XLALGPSSetREAL8(&refTime, (0.5*tObs) + XLALGPSGetREAL8(&firstTimeStamp)); 
  }

  /* set frequency resolution defaults if not set by user */
  if (!(LALUserVarWasSet (&uvar_fResolution))) {
    uvar_fResolution = 1/tObs;
  }

  /*get number of frequency loops*/
  nfreqLoops = ceil(uvar_fBand/uvar_fResolution);

  /* if we are using spindown parameters, initialise the fdots array */
  if (uvar_spindownParams) {

    fdots = (REAL8Vector *)LALCalloc(1, sizeof(REAL8Vector));
    fdots->length = N_SPINDOWN_DERIVS;
    fdots->data = (REAL8 *)LALCalloc(fdots->length, sizeof(REAL8));

    nq1Loops = 1 + ceil(uvar_epsilonBand/uvar_epsilonResolution);
    
    nq2Loops = 1 + ceil(uvar_magfieldBand/uvar_magfieldResolution);

    nnLoops = 1 + ceil(uvar_brakingindexBand/uvar_brakingindexResolution);

    delta_eps = uvar_epsilonResolution;
  
    delta_mag = uvar_magfieldResolution;

    delta_n = uvar_brakingindexResolution; 
  }
  /*otherwise just search over f0, fdot, fddot */
  else {

    if (!(LALUserVarWasSet (&uvar_fdotResolution))) {
      uvar_fdotResolution = SQUARE(1/tObs);
    }

    if (!(LALUserVarWasSet (&uvar_fddotResolution))) {
      uvar_fddotResolution = CUBE(1/tObs);
    }

    nfdotLoops = 1 + ceil(uvar_fdotBand/uvar_fdotResolution);

    nfddotLoops = 1 + ceil(uvar_fddotBand/uvar_fddotResolution);

    delta_fdot = uvar_fdotResolution;
 
    delta_fddot = uvar_fddotResolution;
  }

  /*  set up ephemeris  */
  if(uvar_ephemDir) {
    LALSnprintf(EphemEarth, MAXFILENAMELENGTH, "%s/earth%s.dat",
		uvar_ephemDir, uvar_ephemYear);
    LALSnprintf(EphemSun, MAXFILENAMELENGTH, "%s/sun%s.dat",
		uvar_ephemDir, uvar_ephemYear);
  } else {
    LALSnprintf(EphemEarth, MAXFILENAMELENGTH, "earth%s.dat", uvar_ephemYear);
    LALSnprintf(EphemSun, MAXFILENAMELENGTH, "sun%s.dat", uvar_ephemYear);
  }

  EphemEarth[MAXFILENAMELENGTH-1] = 0;
  EphemSun[MAXFILENAMELENGTH-1] = 0;

  edat = (EphemerisData *)LALCalloc(1, sizeof(EphemerisData));
  (*edat).ephiles.earthEphemeris = EphemEarth;
  (*edat).ephiles.sunEphemeris = EphemSun;

  LAL_CALL( LALLeapSecs(&status, &tmpLeap, &firstTimeStamp, &lsfas), &status);
  (*edat).leap = (INT2)tmpLeap;

  LAL_CALL( LALInitBarycenter( &status, edat), &status);

  /* set up skypatches */
  if ((skytest = fopen(uvar_skyfile, "r")) == NULL) {
    fprintf(stderr, "skyfile doesn't exist\n");
  }

  LAL_CALL( SetUpRadiometerSkyPatches( &status, &skyInfo,
				       uvar_skyfile, uvar_skyRegion,
				       uvar_dAlpha, uvar_dDelta),
	    &status);
  nSkyPatches = skyInfo.numSkyPatches;
  skyAlpha = skyInfo.alpha;
  skyDelta = skyInfo.delta;
  skySizeAlpha = skyInfo.alphaSize;
  skySizeDelta = skyInfo.deltaSize;


  /* curly As */
  /* because we have the option of either averaging over i or not, we
     need to calculate A_{+,x}^2 and A_xA_+ rather than the individual
     values because <A_x> = 0 */
  if (uvar_averageIota) {
    amplitudes.Aplussq = 7.0/15.0;
    amplitudes.Acrosssq = 1.0/3.0;
    amplitudes.AplusAcross = 0;
  } else {
    amplitudes.Aplussq = SQUARE((1.0 + uvar_cosi*uvar_cosi)/2.0);
    amplitudes.Acrosssq = SQUARE(uvar_cosi);
    amplitudes.AplusAcross = (uvar_cosi/2) + (CUBE(uvar_cosi)/2);
  }

  /* polarisation angle */
  if (LALUserVarWasSet(&uvar_psi)) { 
    psi = (REAL8 *) LALCalloc(1, sizeof(REAL8));
    *psi = uvar_psi;
  }
  else {
    psi = NULL;
  }

  /* initialise output arrays */
  nParams = nSkyPatches * nfreqLoops *nfdotLoops * nfddotLoops;

  rho = (REAL8Vector *)LALCalloc(1, sizeof(REAL8Vector));
  stddev = (REAL8Vector *)LALCalloc(1, sizeof(REAL8Vector));
  rho->data = (REAL8 *)LALCalloc(nParams, sizeof(REAL8));
  stddev->data = (REAL8 *)LALCalloc(nParams, sizeof(REAL8));

  /*initialise detector choice*/
  detChoice = uvar_detChoice;


  /* add wings for Doppler modulation and running median block size */
  /* remove fBand from doppWings because we are going bin-by-bin (?) */
  doppWings = (uvar_f0 + (freqCounter*deltaF_SFT)) * VTOT;
  fMin = uvar_f0 - doppWings - uvar_blocksRngMed * deltaF_SFT;
  fMax = uvar_f0 + uvar_fBand + doppWings + uvar_blocksRngMed * deltaF_SFT;


  slidingcounter = 0;
 	   
  /***********start main calculations**************/

  /*outer loop over all sfts in catalog, so that we load only the relevant sfts each time*/
  for(sftcounter=0; sftcounter < (INT4)catalog->length -1; sftcounter++) {
    tmpSFT = NULL;
    tmpPSD = NULL;
    slidingcat = NULL;
    sftPairIndexList = NULL;
    yalpha = NULL;
    ualpha = NULL;
    sigmasq = NULL;

    counter = 0;
    /* throw away first sft from inputSFTs, and first psdvector, frequency, phase vectors, beamfns */
    if (sftcounter > 0) {

      LAL_CALL(DeleteSFTHead(&status, &sftHead), &status);

      LAL_CALL (DeletePSDHead (&status, &psdHead), &status);

      LAL_CALL (DeleteREAL8Head(&status, &freqHead), &status);
	
      LAL_CALL (DeleteREAL8Head (&status, &phaseHead), &status);

      LAL_CALL (DeleteBeamFnHead (&status, &beamHead), &status);


    }
    /* make a second sft catalog with only the ones within maxlag of the current sft*/
    /* do all sfts with same time together */
    while((slidingcounter < (INT4)catalog->length) &&
	  (XLALGPSDiff(&catalog->data[slidingcounter].header.epoch, &catalog->data[sftcounter].header.epoch) <= uvar_maxlag)) {

      inputSFTs = NULL;
	  
      LAL_CALL( CopySFTFromCatalog(&status, catalog, &inputSFTs, fMin, fMax, slidingcounter), &status);

      LAL_CALL( AddSFTtoList(&status, &sftHead, &sftTail, &(inputSFTs->data[0])), &status);

      LAL_CALL( AddPSDtoList(&status, &psdHead, &psdTail, sftTail->sft.data->length),&status);

      tmpSFT = &(sftTail->sft);
      tmpPSD = &(psdTail->psd);
    
      LAL_CALL( LALNormalizeSFT (&status, tmpPSD,
				 tmpSFT, uvar_blocksRngMed),	&status);

      LAL_CALL( AddREAL8toList(&status, &freqHead, &freqTail), &status);

      LAL_CALL( AddREAL8toList(&status, &phaseHead, &phaseTail), &status);

      LAL_CALL( AddBeamFntoList(&status, &beamHead, &beamTail), &status);

      slidingcounter++;
      XLALDestroySFTVector(inputSFTs);

    }

    listLength = slidingcounter - sftcounter;
   
    if (listLength > 1) {

      /* start frequency loop */
      for (freqCounter = 0; freqCounter < nfreqLoops; freqCounter++) {

        f_current = uvar_f0 + (uvar_fResolution*freqCounter);

        if (uvar_spindownParams) { /*if searching over q1, q2, n*/
 
  	  /* Q1 loop */
	  for (q1Counter = 0; q1Counter < nq1Loops; q1Counter++) {

	    eps_current = uvar_epsilon + (delta_eps*q1Counter);

	    /* Q2 loop */
	    for (q2Counter = 0; q2Counter < nq2Loops; q2Counter++) {

	      mag_current = uvar_magfield + (delta_mag*q2Counter);

  	      /* n loop */
	      for (nCounter = 0; nCounter < nnLoops; nCounter++) {

	        n_current = uvar_brakingindex + (delta_n*nCounter);

                CalculateFdots(&status, fdots, f_current, eps_current, mag_current, n_current);
	   
	        /* loop over sky patches -- main calculations  */
	        for (skyCounter = 0; skyCounter < nSkyPatches; skyCounter++) {

	         /* initialize Doppler parameters of the potential source */
	         thisPoint.Alpha = skyAlpha[skyCounter]; 
	         thisPoint.Delta = skyDelta[skyCounter]; 
	         thisPoint.fkdot[0] = f_current;
 		 for (i=1; i < PULSAR_MAX_SPINS; i++) {
	           thisPoint.fkdot[i] = fdots->data[i-1]; 
                 }
	         thisPoint.refTime = refTime;
   
	         /* set sky positions and skypatch sizes  */
	         patchSizeX = skySizeDelta[skyCounter]; 
	         patchSizeY = skySizeAlpha[skyCounter]; 

  
   	         /* get the amplitude modulation coefficients */
	         skypos.longitude = thisPoint.Alpha; 
	         skypos.latitude = thisPoint.Delta; 
	         skypos.system = COORDINATESYSTEM_EQUATORIAL; 
  	

 	         LAL_CALL( GetBeamInfo( &status, beamHead, sftHead, freqHead, phaseHead, skypos, 
				     edat, &thisPoint), &status);
 
  	         /* loop over SFT mini-list to get pairs */
	         ualphacounter = 0;

	         /*  correlate sft pairs  */
 	         sftList = sftHead;
	         psdList = psdHead;
	         freqList = freqHead;
	         phaseList = phaseHead;
	         beamList = beamHead;

	         sft1 = &(sftList->sft);
	         psd1 = &(psdList->psd);
	         freq1 = freqList->val;
	         phase1 = phaseList->val;
	         beamfns1 = &(beamList->beamfn);

	         /*while there are elements in the sft minilist, keep
	          going and check whether it should be paired with SFT1. 
	          there is no need to check the lag as the sfts must satisfy this condition 
 	          already to be in the mini-list*/
	         while (sftList->nextSFT) {
	  	   /*if we are autocorrelating, we want the head to be paired with itself first*/
   	    	   if ((sftList == sftHead) && uvar_autoCorrelate) { 
		     sft2 = &(sftList->sft);
		     psd2 = &(psdList->psd);
	  	     freq2 = freqList->val;
		     phase2 = phaseList->val;
		     beamfns2 = &(beamList->beamfn);
 
 		   sftList = (SFTListElement *)sftList->nextSFT;
		   psdList = (PSDListElement *)psdList->nextPSD;
		   freqList = (REAL8ListElement *)freqList->nextVal;
		   phaseList = (REAL8ListElement *)phaseList->nextVal;
		   beamList = (CrossCorrBeamFnListElement *)beamList->nextBeamfn;
 
		 } else { /*otherwise just step to the next sft*/

		   sftList = (SFTListElement *)sftList->nextSFT;
		   psdList = (PSDListElement *)psdList->nextPSD;
		   freqList = (REAL8ListElement *)freqList->nextVal;
		   phaseList = (REAL8ListElement *)phaseList->nextVal;
		   beamList = (CrossCorrBeamFnListElement *)beamList->nextBeamfn;

	  	   sft2 = &(sftList->sft);
		   psd2 = &(psdList->psd);
	  	   freq2 = freqList->val;
		   phase2 = phaseList->val;
		   beamfns2 = &(beamList->beamfn);
 		 }
       
                 /*strcmp returns 0 if strings are equal, >0 if strings are different*/
                 sameDet = strcmp(sft1->name, sft2->name);

    	         /* if they are different, set sameDet to 1 so that it will match if
	 	   detChoice == DIFFERENT */
      	 	 if (sameDet != 0) { sameDet = 1; }
      
      		 /* however, if detChoice = ALL, then we want sameDet to match it */
      		 if (detChoice == ALL) { sameDet = detChoice; }
	  	  
	         /* decide whether to add this pair or not */
     		 if ((sameDet == (INT4)detChoice)) {

	          /* increment the size of  Y, u, sigmasq vectors by 1  */
 	          yalpha = XLALResizeCOMPLEX16Vector(yalpha, 1 + ualphacounter);
      		  ualpha = XLALResizeCOMPLEX16Vector(ualpha, 1 + ualphacounter);
      		  sigmasq = XLALResizeREAL8Vector(sigmasq, 1 + ualphacounter);

 		  LAL_CALL( LALCorrelateSingleSFTPair( &status, &(yalpha->data[ualphacounter]),
						     sft1, sft2, psd1, psd2, freq1, freq2),
		  	    &status);

	  	  LAL_CALL( LALCalculateSigmaAlphaSq( &status, &sigmasq->data[ualphacounter],
						    freq1, freq2, psd1, psd2),
			    &status);

		  /*if we are averaging over psi and cos(iota), call the simplified 
 	    	    Ualpha function*/
	  	  if (uvar_averagePsi && uvar_averageIota) {
		    LAL_CALL( LALCalculateAveUalpha ( &status, &ualpha->data[ualphacounter], 
						    phase1, phase2, *beamfns1, *beamfns2, 
						    sigmasq->data[ualphacounter]),
			       &status);

		  } else {
		    LAL_CALL( LALCalculateUalpha ( &status, &ualpha->data[ualphacounter], amplitudes,
						 phase1, phase2, *beamfns1, *beamfns2,
						 sigmasq->data[ualphacounter], psi),
			      &status);
		
		  }
		  ualphacounter++;
                }
	      } /*finish loop over sft pairs*/

	      /*printf("finish loop over pairs\n");*/

	      /* calculate rho from Yalpha and Ualpha, if there were pairs */
 	      if (ualphacounter > 0) {
	        tmpstat = 0;
	        LAL_CALL( LALCalculateCrossCorrPower( &status, &tmpstat, yalpha, ualpha),
			&status);

	        rho->data[counter] += tmpstat;

	        /* calculate standard deviation of rho (Eq 4.6) */
	        LAL_CALL( LALNormaliseCrossCorrPower( &status, &tmpstat, ualpha, sigmasq),
			&status); 

	        stddev->data[counter] += tmpstat;
	      }

	      counter++;

 	     } /*finish loop over sky patches*/
	    } /* finish loop over n*/ 


	  } /*finish loop over q2*/

	} /*finish loop over q1*/
 
       
      } /*endif*/

      else { /* if searching through f, fdots instead */
 
	/* frequency derivative loop */
	for (fdotCounter = 0; fdotCounter < nfdotLoops; fdotCounter++) {

	  fdot_current = uvar_fdot + (delta_fdot*fdotCounter);

	  /* frequency double derivative loop */
	  for (fddotCounter = 0; fddotCounter < nfddotLoops; fddotCounter++) {

	    fddot_current = uvar_fddot + (delta_fddot*fddotCounter);

	   
	    /* loop over sky patches -- main calculations  */
	    for (skyCounter = 0; skyCounter < nSkyPatches; skyCounter++) {

	      /* initialize Doppler parameters of the potential source */
	      thisPoint.Alpha = skyAlpha[skyCounter]; 
	      thisPoint.Delta = skyDelta[skyCounter]; 
	      thisPoint.fkdot[0] = f_current;
	      thisPoint.fkdot[1] = fdot_current; 
	      thisPoint.fkdot[2] = fddot_current;
	      thisPoint.fkdot[3] = 0.0;
	      thisPoint.refTime = refTime;
   
	      /* set sky positions and skypatch sizes  */
	      patchSizeX = skySizeDelta[skyCounter]; 
	      patchSizeY = skySizeAlpha[skyCounter]; 

  
	      /* get the amplitude modulation coefficients */
	      skypos.longitude = thisPoint.Alpha; 
	      skypos.latitude = thisPoint.Delta; 
	      skypos.system = COORDINATESYSTEM_EQUATORIAL; 
  	

	      LAL_CALL( GetBeamInfo( &status, beamHead, sftHead, freqHead, phaseHead, skypos, 
				     edat, &thisPoint), &status);

	      /* loop over SFT mini-list to get pairs */
	      ualphacounter = 0;

	      /*  correlate sft pairs  */
 	      sftList = sftHead;
	      psdList = psdHead;
	      freqList = freqHead;
	      phaseList = phaseHead;
	      beamList = beamHead;

	      sft1 = &(sftList->sft);
	      psd1 = &(psdList->psd);
	      freq1 = freqList->val;
	      phase1 = phaseList->val;
	      beamfns1 = &(beamList->beamfn);

	      /*while there are elements in the sft minilist, keep
	        going and check whether it should be paired with SFT1. 
	        there is no need to check the lag as the sfts must satisfy this condition 
 	        already to be in the mini-list*/
	      while (sftList->nextSFT) {
		/*if we are autocorrelating, we want the head to be paired with itself first*/
 		if ((sftList == sftHead) && uvar_autoCorrelate) { 
		  sft2 = &(sftList->sft);
		  psd2 = &(psdList->psd);
	  	  freq2 = freqList->val;
		  phase2 = phaseList->val;
		  beamfns2 = &(beamList->beamfn);

 		  sftList = (SFTListElement *)sftList->nextSFT;
		  psdList = (PSDListElement *)psdList->nextPSD;
		  freqList = (REAL8ListElement *)freqList->nextVal;
		  phaseList = (REAL8ListElement *)phaseList->nextVal;
		  beamList = (CrossCorrBeamFnListElement *)beamList->nextBeamfn;

		} else { /*otherwise just step to the next sft*/

		  sftList = (SFTListElement *)sftList->nextSFT;
		  psdList = (PSDListElement *)psdList->nextPSD;
		  freqList = (REAL8ListElement *)freqList->nextVal;
		  phaseList = (REAL8ListElement *)phaseList->nextVal;
		  beamList = (CrossCorrBeamFnListElement *)beamList->nextBeamfn;

	  	  sft2 = &(sftList->sft);
		  psd2 = &(psdList->psd);
	  	  freq2 = freqList->val;
		  phase2 = phaseList->val;
		  beamfns2 = &(beamList->beamfn);
 		}
      
                /*strcmp returns 0 if strings are equal, >0 if strings are different*/
                sameDet = strcmp(sft1->name, sft2->name);

    	        /* if they are different, set sameDet to 1 so that it will match if
	 	   detChoice == DIFFERENT */
      		if (sameDet != 0) { sameDet = 1; }
      
      		/* however, if detChoice = ALL, then we want sameDet to match it */
      		if (detChoice == ALL) { sameDet = detChoice; }
	  	  
	        /* decide whether to add this pair or not */
     		if ((sameDet == (INT4)detChoice)) {

	          /* increment the size of  Y, u, sigmasq vectors by 1  */
 	          yalpha = XLALResizeCOMPLEX16Vector(yalpha, 1 + ualphacounter);
      		  ualpha = XLALResizeCOMPLEX16Vector(ualpha, 1 + ualphacounter);
      		  sigmasq = XLALResizeREAL8Vector(sigmasq, 1 + ualphacounter);

 		  LAL_CALL( LALCorrelateSingleSFTPair( &status, &(yalpha->data[ualphacounter]),
						     sft1, sft2, psd1, psd2, freq1, freq2),
		  	    &status);

	  	  LAL_CALL( LALCalculateSigmaAlphaSq( &status, &sigmasq->data[ualphacounter],
						    freq1, freq2, psd1, psd2),
			    &status);

		  /*if we are averaging over psi and cos(iota), call the simplified 
 	    	    Ualpha function*/
	  	  if (uvar_averagePsi && uvar_averageIota) {
		    LAL_CALL( LALCalculateAveUalpha ( &status, &ualpha->data[ualphacounter], 
						    phase1, phase2, *beamfns1, *beamfns2, 
						    sigmasq->data[ualphacounter]),
			       &status);

		  } else {
		    LAL_CALL( LALCalculateUalpha ( &status, &ualpha->data[ualphacounter], amplitudes,
						 phase1, phase2, *beamfns1, *beamfns2,
						 sigmasq->data[ualphacounter], psi),
			      &status);
		
		  }
		  ualphacounter++;
                }
	      } /*finish loop over sft pairs*/

	      /*printf("finish loop over pairs\n");*/

	      /* calculate rho from Yalpha and Ualpha, if there were pairs */
 	      if (ualphacounter > 0) {
	        tmpstat = 0;
	        LAL_CALL( LALCalculateCrossCorrPower( &status, &tmpstat, yalpha, ualpha),
			&status);

	        rho->data[counter] += tmpstat;

	        /* calculate standard deviation of rho (Eq 4.6) */
	        LAL_CALL( LALNormaliseCrossCorrPower( &status, &tmpstat, ualpha, sigmasq),
			&status); 

	        stddev->data[counter] += tmpstat;
	      }

	      counter++;

	    } /* finish loop over skypatches */ 


	  } /*finish loop over frequency double derivatives*/

	} /*finish loop over frequency derivatives*/
 
      } /*endelse*/  
    } /*finish loop over frequencies */

      XLALDestroyCOMPLEX16Vector(yalpha);
      XLALDestroyCOMPLEX16Vector(ualpha);

      XLALDestroyREAL8Vector(sigmasq);

    } /*end if listLength > 1 */

  } /* finish loop over all sfts */
  printf("finish loop over all sfts\n");

  counter = 0;

  /* print all variables to file */
  for (freqCounter = 0; freqCounter < nfreqLoops; freqCounter++) {

    f_current = uvar_f0 + (uvar_fResolution*freqCounter);

    /* frequency derivative loop */
    for (fdotCounter = 0; fdotCounter < nfdotLoops; fdotCounter++) {

      fdot_current = uvar_fdot + (delta_fdot*fdotCounter);

      /* frequency double derivative loop */
      for (fddotCounter = 0; fddotCounter < nfddotLoops; fddotCounter++) {
        for (skyCounter = 0; skyCounter < nSkyPatches; skyCounter++) { 
   	  /* initialize Doppler parameters of the potential source */
	  thisPoint.Alpha = skyAlpha[skyCounter]; 
	  thisPoint.Delta = skyDelta[skyCounter]; 

	  /*normalise rho*/
	  rho->data[counter] = rho->data[counter]/sqrt(stddev->data[counter]);
	  fprintf(fp, "%1.5f\t %1.5f\t %1.5f\t %e\t %e\t %1.10f\n", thisPoint.Alpha,
		  thisPoint.Delta, f_current,
		  fdot_current, fddot_current, rho->data[counter]);
	  counter++;
	}
      }
    }
  }

  /* select candidates  */



  fclose (fp);

  LAL_CALL( LALDestroySFTCatalog( &status, &catalog ), &status);

  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);
  LALFree(skyAlpha);
  LALFree(skyDelta);
  LALFree(skySizeAlpha);
  LALFree(skySizeDelta);
  XLALDestroyREAL8Vector(stddev);
  XLALDestroyREAL8Vector(rho);


  /*free the last few elements (if they're not already free). */
  if (beamHead) {
    if (beamHead != beamTail) { XLALFree(beamTail); }
    XLALFree(beamHead);
  }

  if(sftHead) {
    if (sftHead != sftTail) {
      tmpSFT = &(sftTail->sft);
      LAL_CALL(LALDestroySFTtype(&status, &tmpSFT),&status);
    }
    tmpSFT = &(sftHead->sft);
    LAL_CALL(LALDestroySFTtype(&status, &tmpSFT),&status);
  } 
  if (psdHead) {
    if (psdHead != psdTail) {
      XLALDestroyREAL8FrequencySeries(&(psdTail->psd));
    }

    XLALDestroyREAL8FrequencySeries (&(psdHead->psd));
  }

  if (phaseHead) {
    if (phaseHead != phaseTail) { XLALFree(phaseTail);}
    XLALFree(phaseHead);
  }

  if (freqHead) {
    if (freqHead != freqTail) {XLALFree(freqTail);}
    XLALFree(freqHead);
  }

  LAL_CALL (LALDestroyUserVars(&status), &status);

  LALCheckMemoryLeaks();

  if ( lalDebugLevel )
    REPORTSTATUS ( &status);

  return status.statusCode;

} /* main */


/** Set up location of skypatch centers and sizes
    If user specified skyRegion then use DopplerScan function
    to construct an isotropic grid. Otherwise use skypatch file. */
void SetUpRadiometerSkyPatches(LALStatus           *status,
			       SkyPatchesInfo      *out,   /**< output skypatches info */
			       CHAR                *skyFileName, /**< name of skypatch file */
			       CHAR                *skyRegion,  /**< skyregion (if isotropic grid is to be constructed) */
			       REAL8               dAlpha,      /**< alpha resolution (if isotropic grid is to be constructed) */
			       REAL8               dDelta)  /**< delta resolution (if isotropic grid is to be constructed) */
{

  DopplerSkyScanInit scanInit = empty_DopplerSkyScanInit;   /* init-structure for DopperScanner */
  DopplerSkyScanState thisScan = empty_DopplerSkyScanState; /* current state of the Doppler-scan */
  UINT4 nSkyPatches, skyCounter;
  PulsarDopplerParams dopplerpos;

  INITSTATUS (status, "SetUpRadiometerSkyPatches", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT (out, status, PULSAR_CROSSCORR_ENULL, PULSAR_CROSSCORR_MSGENULL);
  ASSERT (dAlpha > 0, status, PULSAR_CROSSCORR_EBAD, PULSAR_CROSSCORR_MSGEBAD);
  ASSERT (dDelta > 0, status, PULSAR_CROSSCORR_EBAD, PULSAR_CROSSCORR_MSGEBAD);

  if (skyRegion ) {

    scanInit.dAlpha = dAlpha;
    scanInit.dDelta = dDelta;
    scanInit.gridType = GRID_ISOTROPIC;
    scanInit.metricType =  LAL_PMETRIC_NONE;
    scanInit.skyRegionString = (CHAR*)LALCalloc(1, strlen(skyRegion)+1);
    strcpy (scanInit.skyRegionString, skyRegion);
    /*   scanInit.Freq = usefulParams.spinRange_midTime.fkdot[0] +  usefulParams.spinRange_midTime.fkdotBand[0]; */

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

      if ((dopplerpos.Delta>0) && (dopplerpos.Delta < atan(4*LAL_PI/dAlpha/dDelta) ))
        out->alphaSize[skyCounter] = dAlpha*cos(dopplerpos.Delta -0.5*dDelta)/cos(dopplerpos.Delta);

      if ((dopplerpos.Delta<0) && (dopplerpos.Delta > -atan(4*LAL_PI/dAlpha/dDelta) ))
        out->alphaSize[skyCounter] = dAlpha*cos(dopplerpos.Delta +0.5*dDelta)/cos(dopplerpos.Delta);

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

      ASSERT (skyFileName, status, PULSAR_CROSSCORR_ENULL, PULSAR_CROSSCORR_MSGENULL);

      if ( (fpsky = fopen(skyFileName, "r")) == NULL) {
	ABORT ( status, PULSAR_CROSSCORR_EFILE, PULSAR_CROSSCORR_MSGEFILE );
      }

      nSkyPatches = 0;
      do {
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

      for (skyCounter = 0; skyCounter < nSkyPatches; skyCounter++) {
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


void GetBeamInfo(LALStatus *status, 
		 CrossCorrBeamFnListElement *beamHead, 
		 SFTListElement *sftHead, 
		 REAL8ListElement *freqHead,
		 REAL8ListElement *phaseHead, 
		 SkyPosition skypos, 
		 EphemerisData *edat, 
		 PulsarDopplerParams *thisPoint){

  REAL8 freq1;
  REAL8 phase1;
  REAL8Vector thisVel, thisPos;
  LIGOTimeGPSVector *ts=NULL;
  LALDetector *det;
  REAL8 tOffs;
  AMCoeffs *AMcoef = NULL; 
  SFTListElement *sft = NULL;
  REAL8ListElement *freqtmp, *phasetmp;
  CrossCorrBeamFnListElement *beamtmp;
  LIGOTimeGPS *epoch = NULL;

  INITSTATUS (status, "GetBeamInfo", rcsid);
  ATTATCHSTATUSPTR (status);

  freq1 = 0;
  phase1 = 0;


  sft = sftHead;
  freqtmp = freqHead;
  phasetmp = phaseHead;
  beamtmp = beamHead;

  ts = XLALCreateTimestampVector(1);


  thisVel.length = 3;
  thisPos.length = 3;
  tOffs = 0.5/sft->sft.deltaF;


  /* get information about all detectors including velocity and
     timestamps */
  /*only have 1 element in detectorStateSeries and AMCoeffs because
    the detector has to be updated for every SFT */
  while (sft) {

    DetectorStateSeries *detState = NULL;
    AMcoef = (AMCoeffs *)LALCalloc(1, sizeof(AMCoeffs));
    AMcoef->a = XLALCreateREAL4Vector(1);
    AMcoef->b = XLALCreateREAL4Vector(1);

    epoch = &(sft->sft.epoch);
    det = XLALGetSiteInfo (sft->sft.name); 
    ts->data[0] = sft->sft.epoch;

    /* note that this function returns the velocity at the
       mid-time of the SFTs -- should not make any difference */

    LALGetDetectorStates ( status->statusPtr, &detState, ts, det,
			   edat, tOffs);

    detState->detector = *det;

    LALNewGetAMCoeffs ( status->statusPtr, AMcoef, detState, skypos);
		     
    thisVel.data = detState->data[0].vDetector;
    thisPos.data = detState->data[0].rDetector;


    LALGetSignalFrequencyInSFT( status->statusPtr, &freq1, epoch, thisPoint,
				&thisVel);

    LALGetSignalPhaseInSFT( status->statusPtr, &phase1, epoch, thisPoint,
			    &thisPos);

    freqtmp->val = freq1; 
    phasetmp->val = phase1;

    /* store a and b in the CrossCorrBeamFn */

      beamtmp->beamfn.a = (AMcoef->a->data[0]);
      beamtmp->beamfn.b = (AMcoef->b->data[0]);
    
		
    /* clean up AMcoefs */
    XLALDestroyAMCoeffs(AMcoef);
    XLALDestroyDetectorStateSeries(detState);
    XLALFree(det);
    sft = (SFTListElement *)sft->nextSFT;
    freqtmp = (REAL8ListElement *)freqtmp->nextVal;
    phasetmp = (REAL8ListElement *)phasetmp->nextVal;
    beamtmp = (CrossCorrBeamFnListElement *)beamtmp->nextBeamfn;
  }

  XLALDestroyTimestampVector(ts);
  sft = NULL;
  freqtmp = NULL;
  phasetmp = NULL;
  beamtmp = NULL;  

  DETATCHSTATUSPTR (status);

  RETURN (status);

}


void CopySFTFromCatalog(LALStatus *status,
		   	SFTCatalog *catalog,
			SFTVector **sft,
			REAL8 fMin,
			REAL8 fMax,
			INT4 sftindex) 
{
  SFTCatalog *slidingcat;
  SFTDescriptor *desc;

  INITSTATUS (status, "CopySFTFromCatalog", rcsid);
  ATTATCHSTATUSPTR (status);


  slidingcat = NULL;
  *sft = NULL;

  ASSERT ( catalog, status, PULSAR_CROSSCORR_ENULL, PULSAR_CROSSCORR_MSGENULL );

  /*check that we are loading an sensible frequency range*/
  if (fMin < catalog->data[sftindex].header.f0 || fMax > (catalog->data[sftindex].header.f0 + catalog->data[sftindex].numBins*catalog->data[sftindex].header.deltaF)) {
    ABORT(status, PULSAR_CROSSCORR_EVAL, PULSAR_CROSSCORR_MSGEVAL);
  }

  if ( (slidingcat = LALCalloc ( 1, sizeof ( SFTCatalog ))) == NULL ) {
    ABORT(status, PULSAR_CROSSCORR_EMEM, PULSAR_CROSSCORR_MSGEMEM);
  }

  slidingcat->length = 1;
  slidingcat->data = LALRealloc(slidingcat->data, 1*sizeof(*(slidingcat->data)));
  /* memset(&(slidingcat->data[0]), 0, sizeof(ret->data[0]));*/
  desc = &(slidingcat->data[0]);
  /*  desc->locator = LALCalloc(1, sizeof(*(desc->locator)));*/

  desc->locator = catalog->data[sftindex].locator;
  desc->header = catalog->data[sftindex].header;
  desc->comment = catalog->data[sftindex].comment;
  desc->numBins = catalog->data[sftindex].numBins;
  desc->version = catalog->data[sftindex].version;
  desc->crc64 = catalog->data[sftindex].crc64;

  LALLoadSFTs(status->statusPtr, sft, slidingcat, fMin, fMax); 
  /* LALDestroySFTCatalog ( status->statusPtr, &slidingcat );*/
  XLALFree(slidingcat->data);
  LALFree(slidingcat);

  DETATCHSTATUSPTR (status);
	
  RETURN (status);

}

void AddSFTtoList(LALStatus *status,
		  SFTListElement **sftHead,
		  SFTListElement **sftTail,
		  SFTtype *sft)
{

  SFTListElement *sftList;

  INITSTATUS (status, "AddSFTtoList", rcsid);
  ATTATCHSTATUSPTR (status);

  sftList = (SFTListElement *)LALCalloc(1, sizeof(SFTListElement));

  LALCopySFT(status->statusPtr, &(sftList->sft), sft);

  sftList->nextSFT = NULL;
  if (!(*sftHead)) { 
    *sftHead = sftList; 
  }
  else {
    (*sftTail)->nextSFT = (struct SFTListElement *)sftList;
  }
	
  *sftTail = sftList;

  DETATCHSTATUSPTR (status);
	
  RETURN (status);


}

void AddPSDtoList(LALStatus *status,
		  PSDListElement **psdHead,
		  PSDListElement **psdTail,
		  INT4 length)
{
  PSDListElement *psdList;

  INITSTATUS (status, "AddPSDtoList", rcsid);
  ATTATCHSTATUSPTR (status);

  psdList = (PSDListElement *)LALCalloc(1, sizeof(PSDListElement));
  psdList->psd.data = XLALCreateREAL8Sequence(length);
  psdList->nextPSD = NULL;
  if (!(*psdHead)) {
    *psdHead = psdList;
  } else {
    (*psdTail)->nextPSD = (struct PSDListElement *)psdList;
  }
  *psdTail = psdList;

  DETATCHSTATUSPTR (status);
	
  RETURN (status);


}

void AddREAL8toList(LALStatus *status,
		    REAL8ListElement **head,
		    REAL8ListElement **tail)
{
  REAL8ListElement *List;

  INITSTATUS (status, "AddREAL8toList", rcsid);
  ATTATCHSTATUSPTR (status);


  List = (REAL8ListElement *)LALCalloc(1, sizeof(REAL8ListElement));
  List->val = 0;
  List->nextVal = NULL;
  if (!(*head)) {
    *head = List;
  } else {
    (*tail)->nextVal = (struct REAL8ListElement *)List;
  }
  *tail = List;

  DETATCHSTATUSPTR (status);
	
  RETURN (status);


}

void AddBeamFntoList(LALStatus *status,
		     CrossCorrBeamFnListElement **beamHead,
		     CrossCorrBeamFnListElement **beamTail)
{

  CrossCorrBeamFnListElement *beamList;

  INITSTATUS (status, "AddBeamFntoList", rcsid);
  ATTATCHSTATUSPTR (status);


  beamList = (CrossCorrBeamFnListElement *)LALCalloc(1, sizeof(CrossCorrBeamFnListElement));
  beamList->beamfn.a = 0;
  beamList->beamfn.b = 0;
  beamList->nextBeamfn = NULL;
  if (!(*beamHead)) {
    *beamHead = beamList;	
  } else {
    (*beamTail)->nextBeamfn = (struct CrossCorrBeamFnListElement *)beamList;
  }
  *beamTail = beamList;

  DETATCHSTATUSPTR (status);
	
  RETURN (status);

}

void DeleteSFTHead (LALStatus *status, 
		    SFTListElement **sftHead)
{
  SFTListElement *sftList;
  SFTtype *tmpSFT;

  INITSTATUS (status, "DeleteSFTHead", rcsid);
  ATTATCHSTATUSPTR (status);


  sftList = *sftHead;
  *sftHead = (SFTListElement *)(*sftHead)->nextSFT;

  tmpSFT = &(sftList->sft);
  LALDestroySFTtype(status->statusPtr, &tmpSFT);
  sftList = NULL;

  DETATCHSTATUSPTR (status);
	
  RETURN (status);

}

void DeletePSDHead (LALStatus *status, 
		    PSDListElement **psdHead)
{
  PSDListElement *psdList;
  REAL8FrequencySeries *tmpPSD;

  INITSTATUS (status, "DeletePSDHead", rcsid);
  ATTATCHSTATUSPTR (status);

  psdList = *psdHead;
  *psdHead = (PSDListElement *)(*psdHead)->nextPSD;
  tmpPSD = &(psdList->psd);
  XLALDestroyREAL8FrequencySeries (tmpPSD);
  psdList = NULL;

  DETATCHSTATUSPTR (status);
	
  RETURN (status);

}

void DeleteREAL8Head (LALStatus *status,
		      REAL8ListElement **head)
{
  REAL8ListElement *List;
  REAL8 *tmpVal;

  INITSTATUS (status, "DeleteREAL8Head", rcsid);
  ATTATCHSTATUSPTR (status);


  List = *head;
  *head = (REAL8ListElement *)(*head)->nextVal;
  tmpVal = &(List->val);
  LALFree(tmpVal);
  List = NULL;

  DETATCHSTATUSPTR (status);
	
  RETURN (status);

}

void DeleteBeamFnHead (LALStatus *status,
		       CrossCorrBeamFnListElement **beamHead)
{
  CrossCorrBeamFnListElement *beamList;
  CrossCorrBeamFn *beamfns;

  INITSTATUS (status, "DeleteBeamFnHead", rcsid);
  ATTATCHSTATUSPTR (status);

  beamList = *beamHead;
  *beamHead = (CrossCorrBeamFnListElement *)(*beamHead)->nextBeamfn;
  beamfns = &(beamList->beamfn);
  LALFree(beamfns);
  beamList = NULL;

  DETATCHSTATUSPTR (status);
	
  RETURN (status);

}

void CalculateFdots (LALStatus *status,
		     REAL8Vector *fdots,
		     REAL8 f0,
		     REAL8 epsilon,
		     REAL8 magfield,
		     REAL8 n)
{
  REAL8 grav, rstar, inertia, c, mu0;
  REAL8 q1, q2;


  INITSTATUS (status, "CalculateFdots", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT(fdots->length >= 6, status, PULSAR_CROSSCORR_ENULL, PULSAR_CROSSCORR_MSGENULL);


  grav = 6.67e-11;
  rstar = 1e4;
  c = 2.998e8;
  inertia = 0.4*1.4*1.99e30*rstar*rstar;
  mu0 = 1.2566e-6;

  q1 = 32*grav*inertia*SQUARE(LAL_PI)*SQUARE(LAL_PI)*SQUARE(epsilon)/(5*CUBE(c)*SQUARE(c));

  q2 = 2*CUBE(LAL_PI)*CUBE(rstar)*CUBE(rstar)*SQUARE(magfield)/(3*mu0*inertia*CUBE(c));

  /* multiply q2 by the (pi R/c)^n-3 term*/
  q2 = q2*pow(LAL_PI*rstar/c, n-3);

  /* hard code each derivative. symbolic differentiation too hard? */
  fdots->data[0] = -(q1 * pow(f0, 5)) - (q2 * pow(f0, n));

  fdots->data[1] = -(5.0 * q1 * pow(f0, 4) * fdots->data[0])
		   -(n * q2 * pow(f0, n-1) * fdots->data[0]);

  fdots->data[2] = -q1 * (20.0 * CUBE(f0) * SQUARE(fdots->data[0]) + 5.0 * pow(f0, 4) * fdots->data[1])
		   -q2 * ((n-1) * n * pow(f0, n-2) * SQUARE(fdots->data[0]) + n * pow(f0,n-1) * fdots->data[1]);

  fdots->data[3] = -q1 * (60.0*SQUARE(f0)*CUBE(fdots->data[0]) + 60.0*CUBE(f0)*fdots->data[0]*fdots->data[1] 
			  + 5.0*pow(f0,4)*fdots->data[2])
		   -q2 * ((n-2)*(n-1)*n*pow(f0,n-3)*CUBE(fdots->data[0])
			   + 3*n*(n-1)*pow(f0,n-2)*fdots->data[0]*fdots->data[1] + n*pow(f0, n-1)*fdots->data[2]);

  fdots->data[4] = -q1 * (120.0*f0*pow(fdots->data[0],4) + 360.0*SQUARE(f0)*SQUARE(fdots->data[0])*fdots->data[1] 
		   	+ 60.0*CUBE(f0)*SQUARE(fdots->data[1]) + 80.0*CUBE(f0)*fdots->data[0]*fdots->data[2] 
			+ 5.0*SQUARE(f0)*SQUARE(f0)*fdots->data[3])
		   -q2 * ((n-3)*(n-2)*(n-1)*n*pow(f0,n-4)*pow(fdots->data[0],4)
			  + 6.0*(n-2)*(n-1)*n*pow(f0,n-3)*SQUARE(fdots->data[0])*fdots->data[1] 
			  + 3.0*(n-1)*n*pow(f0,n-2)*SQUARE(fdots->data[1]) 
			  + 4.0*(n-1)*n*pow(f0,n-2)*fdots->data[0]*fdots->data[2] + n*pow(f0, n-1)*fdots->data[3]);

  fdots->data[5] = -q1 * (120.0*pow(fdots->data[0],5) + 1200.0*f0*CUBE(fdots->data[0])*fdots->data[1] 
			  + 900.0*SQUARE(f0)*fdots->data[0]*SQUARE(fdots->data[1]) + 600.0*SQUARE(f0)*SQUARE(fdots->data[0])*fdots->data[2]
			  + 200.0*CUBE(f0)*fdots->data[1]*fdots->data[2] + 100.0*CUBE(f0)*fdots->data[0]*fdots->data[3]
			  + 5.0*pow(f0,4)*fdots->data[4])
		   -q2 * ((n-4)*(n-3)*(n-2)*(n-1)*n*pow(f0,n-5)*pow(fdots->data[0],5) 
			 + 10.0*(n-3)*(n-2)*(n-1)*n*pow(f0,n-4)*CUBE(fdots->data[0])*fdots->data[1]
			 + 15.0*(n-2)*(n-1)*n*pow(f0,n-3)*fdots->data[0]*SQUARE(fdots->data[1])
			 + 10.0*(n-2)*(n-1)*n*pow(f0,n-3)*SQUARE(fdots->data[1])*fdots->data[2]
			 + 10.0*(n-1)*n*pow(f0,n-2)*fdots->data[1]*fdots->data[2]
 			 + 5.0*(n-1)*n*pow(f0,n-2)*fdots->data[0]*fdots->data[3]
			 + n*pow(f0,n-1)*fdots->data[4]);

/*  for (i=1; i < fdots->length; i++) {
    counter = 1;
    gwcounter = gwBrakingIndex;
    emcounter = emBrakingIndex;
    gwfactorial = 1;
    emfactorial = 1;
    fdotfactorial = 1;

    while (counter <= i) {
 
      gwfactorial *= gwcounter;
      emfactorial *= emcounter;
      fdotfactorial *=

      fdots->data[i] += -(gwfactorial * q1 * pow(f0, gwcounter-1) * pow(fdots->data[i-counter], counter) )
 		        -(emfactorial * q2 * pow(fdots->data[i-counter], counter) *pow(f0, emcounter-1))
      gwcounter--;
      emcounter--;
      counter++;
    }
  }
*/

  DETATCHSTATUSPTR (status);

  RETURN (status);
}

void initUserVars (LALStatus *status)
{
  INITSTATUS( status, "initUserVars", rcsid );
  ATTATCHSTATUSPTR (status);


  uvar_maxlag = 0;

  uvar_help = FALSE;
  uvar_averagePsi = TRUE;
  uvar_averageIota = TRUE;
  uvar_autoCorrelate = FALSE;
  uvar_spindownParams = FALSE;
  uvar_blocksRngMed = BLOCKSRNGMED;
  uvar_detChoice = 2;
  uvar_f0 = F0;
  uvar_fResolution = 0.0;
  uvar_startTime = 0.0;
  uvar_endTime = LAL_INT4_MAX;
  uvar_fdot = 0.0;
  uvar_fdotBand = 0.0;
  uvar_fdotResolution = 0.0;
  uvar_fddot = 0.0;
  uvar_fddotBand = 0.0;
  uvar_fddotResolution = 0.0;
  uvar_fBand = FBAND;
  uvar_dAlpha = 0.2;
  uvar_dDelta = 0.2;
  uvar_psi = 0.0;
  uvar_refTime = 0.0;
  uvar_cosi = 0.0;
  uvar_epsilon = 1e-5;
  uvar_epsilonBand = 0.0;
  uvar_epsilonResolution = 0.0;
  uvar_magfield = 1e9;
  uvar_magfieldBand = 0.0;
  uvar_magfieldResolution = 0.0;
  uvar_brakingindex = 3;
  uvar_brakingindexBand = 0.0;
  uvar_brakingindexResolution = 0.0;

  uvar_ephemDir = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_ephemDir,DEFAULT_EPHEMDIR);

  uvar_ephemYear = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_ephemYear,EPHEM_YEARS);

  uvar_dirnameOut = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_dirnameOut,DIROUT);

  uvar_filenameOut = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_filenameOut,FILEOUT);

  uvar_skyfile = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_skyfile,SKYFILE);

  /* register user input variables */
  LALRegisterBOOLUserVar( status->statusPtr, "help",
			  'h', UVAR_HELP,
			  "Print this message",
			  &uvar_help);
  LALRegisterBOOLUserVar( status->statusPtr, "averagePsi",
			  0, UVAR_OPTIONAL,
			  "Use average over psi",
			  &uvar_averagePsi);
  LALRegisterBOOLUserVar( status->statusPtr, "averageIota",
			  0, UVAR_OPTIONAL,
			  "Use average over iota",
			  &uvar_averageIota);
  LALRegisterBOOLUserVar( status->statusPtr, "autoCorrelate",
			  0, UVAR_OPTIONAL,
			  "Include autocorrelations",
			  &uvar_autoCorrelate);
  LALRegisterREALUserVar( status->statusPtr, "f0",
			  'f', UVAR_OPTIONAL,
			  "Start search frequency",
			  &uvar_f0);
  LALRegisterREALUserVar( status->statusPtr, "fdot",
			  0, UVAR_OPTIONAL,
			  "Start search frequency derivative",
			  &uvar_fdot);
  LALRegisterREALUserVar( status->statusPtr, "fBand",
			  'b', UVAR_OPTIONAL,
			  "Search frequency band",
			  &uvar_fBand);
  LALRegisterREALUserVar( status->statusPtr, "fRes",
			  0, UVAR_OPTIONAL,
			  "Search frequency resolution. Default: 1/T",
			  &uvar_fResolution);
  LALRegisterREALUserVar( status->statusPtr, "fdotBand",
			  0, UVAR_OPTIONAL,
			  "Search frequency derivative band",
			  &uvar_fdotBand);
  LALRegisterREALUserVar( status->statusPtr, "fdotRes",
			  'r', UVAR_OPTIONAL,
			  "Search frequency derivative resolution. Default: 1/T^2",
			  &uvar_fdotResolution);
  LALRegisterREALUserVar( status->statusPtr, "fddot",
			  0, UVAR_OPTIONAL,
			  "Start frequency double derivative",
			  &uvar_fddot);
  LALRegisterREALUserVar( status->statusPtr, "fddotBand",
			  0, UVAR_OPTIONAL,
			  "Search frequency double derivative band",
			  &uvar_fddotBand);
  LALRegisterREALUserVar( status->statusPtr, "fddotRes",
			  0, UVAR_OPTIONAL,
			  "Search frequency double derivative resolution. Default: 1/T^3",
			  &uvar_fddotResolution);
  LALRegisterREALUserVar( status->statusPtr, "startTime",
			  0, UVAR_OPTIONAL,
			  "GPS start time of observation",
			  &uvar_startTime);
  LALRegisterREALUserVar( status->statusPtr, "endTime",
			  0, UVAR_OPTIONAL,
			  "GPS end time of observation",
			  &uvar_endTime);
  LALRegisterSTRINGUserVar( status->statusPtr, "skyRegion",
			    0, UVAR_OPTIONAL,
			    "sky-region polygon (or 'allsky')",
			    &uvar_skyRegion);
  LALRegisterREALUserVar( status->statusPtr, "dAlpha",
			  0, UVAR_OPTIONAL,
			  "Sky resolution (flat/isotropic) (rad)",
			  &uvar_dAlpha);
  LALRegisterREALUserVar( status->statusPtr, "dDelta",
			  0, UVAR_OPTIONAL,
			  "Sky resolution (flat/isotropic) (rad)",
			  &uvar_dDelta);
  LALRegisterREALUserVar( status->statusPtr, "psi",
			  0, UVAR_OPTIONAL,
			  "Polarisation angle (rad)",
			  &uvar_psi);
  LALRegisterSTRINGUserVar( status->statusPtr, "skyfile",
			    0, UVAR_OPTIONAL,
			    "Alternative: input skypatch file",
			    &uvar_skyfile);
  LALRegisterSTRINGUserVar( status->statusPtr, "ephemDir",
			    'E', UVAR_OPTIONAL,
			    "Directory where ephemeris files are located",
			    &uvar_ephemDir);
  LALRegisterSTRINGUserVar( status->statusPtr, "ephemYear",
			    'y', UVAR_OPTIONAL,
			    "Year (or range of years) of ephemeris files to be used",
			    &uvar_ephemYear);
  LALRegisterSTRINGUserVar( status->statusPtr, "sftDir",
			    'D', UVAR_REQUIRED,
			    "SFT filename pattern",
			    &uvar_sftDir);
  LALRegisterREALUserVar( status->statusPtr, "maxlag",
			  0,  UVAR_OPTIONAL,
			  "Maximum time lag for correlating sfts",
			  &uvar_maxlag);
  LALRegisterSTRINGUserVar( status->statusPtr, "dirnameOut",
			    'o', UVAR_OPTIONAL,
			    "Output directory",
			    &uvar_dirnameOut);
  LALRegisterSTRINGUserVar( status->statusPtr, "filenameOut",
			    'o', UVAR_OPTIONAL,
			    "Output filename",
			    &uvar_filenameOut);
  LALRegisterINTUserVar( status->statusPtr, "blocksRngMed",
			 0, UVAR_OPTIONAL,
			 "Running Median block size",
			 &uvar_blocksRngMed);
  LALRegisterINTUserVar( status->statusPtr, "detChoice",
			 0, UVAR_OPTIONAL,
			 "0: Correlate SFTs from same IFOs only; 1: different IFOs only; 2: all IFOs",
			 &uvar_detChoice);
  LALRegisterREALUserVar( status->statusPtr, "refTime",
			  0, UVAR_OPTIONAL,
			  "Pulsar reference time (SSB)",
			  &uvar_refTime);
  LALRegisterREALUserVar( status->statusPtr, "cosi",
			  0, UVAR_OPTIONAL,
			  "cos(iota) inclination angle",
			  &uvar_cosi); 
  LALRegisterREALUserVar( status->statusPtr, "epsilon",
			  0, UVAR_OPTIONAL,
			  "Start pulsar ellipticity",
			  &uvar_epsilon); 
  LALRegisterREALUserVar( status->statusPtr, "epsilonBand",
			  0, UVAR_OPTIONAL,
			  "Pulsar ellipticity search band",
			  &uvar_epsilonBand); 
  LALRegisterREALUserVar( status->statusPtr, "epsilonRes",
			  0, UVAR_OPTIONAL,
			  "Pulsar ellipticity search resolution",
			  &uvar_epsilonResolution); 
  LALRegisterREALUserVar( status->statusPtr, "magfield",
			  0, UVAR_OPTIONAL,
			  "Start pulsar magnetic field (T)",
			  &uvar_magfield); 
  LALRegisterREALUserVar( status->statusPtr, "magfieldBand",
			  0, UVAR_OPTIONAL,
			  "Pulsar magnetic field search band (T)",
			  &uvar_magfieldBand); 
  LALRegisterREALUserVar( status->statusPtr, "magfieldRes",
			  0, UVAR_OPTIONAL,
			  "Pulsar magnetic field search resolution (T)",
			  &uvar_magfieldResolution); 
  LALRegisterREALUserVar( status->statusPtr, "braking",
			  0, UVAR_OPTIONAL,
			  "Pulsar electromagnetic braking index",
			  &uvar_brakingindex); 
  LALRegisterREALUserVar( status->statusPtr, "brakingBand",
			  0, UVAR_OPTIONAL,
			  "Pulsar electromagnetic braking index search band",
			  &uvar_brakingindexBand); 
  LALRegisterREALUserVar( status->statusPtr, "brakingRes",
			  0, UVAR_OPTIONAL,
			  "Pulsar electromagnetic braking index search resolution",
			  &uvar_brakingindexResolution); 
  LALRegisterBOOLUserVar( status->statusPtr, "spindownParams",
			  0, UVAR_OPTIONAL,
			  "Search over pulsar spindown parameters instead of frequency derivatives",
			  &uvar_spindownParams); 


  DETATCHSTATUSPTR (status);
  RETURN (status);
}


