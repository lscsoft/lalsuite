/*
 *  Copyright (C) 2007 Badri Krishnan  
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



#include "./radiometer.h"
/* lalapps includes */
#include <lalapps.h>
#include <lal/DopplerScan.h>
#include <gsl/gsl_permutation.h>

RCSID( "$Id$");



/* globals, constants and defaults */


extern int lalDebugLevel;


#define EARTHEPHEMERIS "/Users/artax/work/opt/lscsoft/lal/share/lal/earth05-09.dat"
#define SUNEPHEMERIS "/Users/artax/work/opt/lscsoft/lal/share/lal/sun05-09.dat"

#define F0 100
#define FBAND 1
#define LAL_INT4_MAX 2147483647


#define BLOCKSRNGMED 51
#define MAXFILENAMELENGTH 512 /* maximum # of characters  of a filename */

#define DIROUT "./out"   /* output directory */
#define BASENAMEOUT "radio"    /* prefix file output */

#define SKYFILE "./skypatchfile"      
#define SKYREGION "allsky" 

#define TRUE (1==1)
#define FALSE (1==0)


int main(int argc, char *argv[]){
  /* LALStatus pointer */
  static LALStatus  status;  
  
  /* sft related variables */ 

  SFTVector *inputSFTs = NULL;
  MultiSFTVector *multiSFTs = NULL;
  REAL8 deltaF, timeBase, freq = 0, *freq1 = &freq;
  REAL8 phase = 0, *phase1 = &phase, psi = 0.0; 
  INT8  fLastBin, f0Bin;
  UINT4 binsSFT, numSearchBins, numsft, counter, paircounter;
  INT4 index1, index2;
  COMPLEX8FrequencySeries *sft, *sft1, *sft2;
  REAL8FrequencySeries *psd1, *psd2;
  COMPLEX16Vector *yalpha, *ualpha;
  REAL8Vector *Fplus, *Fcross;
  REAL8 Aplus = sqrt(7.0/15.0), Across = sqrt(1.0/3.0), rho0 = 0.0; 
  /*A's are averaged over cos i... need to change this later!*/ 

  /*SFTPairVec sftPairs;*/
  SFTPairParams pairParams;
  
  /* information about all the ifos */
  DetectorStateSeries *mdetStates = NULL, *tmpDetState = NULL;
  PSDVector *psdVec = NULL;  
  UINT4 numifo;
  LALDetector *det;

  LIGOTimeGPS firstTimeStamp, lastTimeStamp;
  REAL8 tObs, tOffs;
  REAL8 patchSizeX, patchSizeY;

  /* ephemeris */
  EphemerisData    *edat=NULL;
  INT4 tmpLeap;
  LALLeapSecFormatAndAcc  lsfas = {LALLEAPSEC_GPSUTC, LALLEAPSEC_LOOSE};

  /* skypatch info */
  REAL8  *skyAlpha=NULL, *skyDelta=NULL, *skySizeAlpha=NULL, *skySizeDelta=NULL; 
  INT4   nSkyPatches, skyCounter=0; 
  SkyPatchesInfo skyInfo;

  static INT4VectorSequence  *sftPairIndexList;
  static REAL8Vector *frequencyShiftList, *signalPhaseList;
  INT4  j;
  static PulsarDopplerParams  thisPoint;
  static REAL8Vector thisVel, thisPos, *weights;

    REAL8 doppWings, fMin, fMax;

  AMCoeffs *AMcoef = NULL; 



  /* sft constraint variables */
  LIGOTimeGPS startTimeGPS, endTimeGPS;
  LIGOTimeGPSVector *ts=NULL;

  /* user input variables */
  BOOLEAN  uvar_help;
  INT4     uvar_blocksRngMed; 
  REAL8    uvar_startTime, uvar_endTime;
  REAL8    uvar_f0, uvar_fdot, uvar_fBand;
  REAL8    uvar_dAlpha, uvar_dDelta; /* resolution for isotropic sky-grid */
  REAL8    uvar_maxlag;
  CHAR     *uvar_earthEphemeris=NULL;
  CHAR     *uvar_sunEphemeris=NULL;
  CHAR     *uvar_sftDir=NULL;
  CHAR     *uvar_dirnameOut=NULL;
  CHAR     *uvar_fbasenameOut=NULL;
  CHAR     *uvar_skyfile=NULL;
  CHAR     *uvar_skyRegion=NULL;
  FILE	   *skytest=NULL;


  SkyPosition skypos; 


  /* LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;
  
  lalDebugLevel = 0;  /* LALDebugLevel must be called before anything else */
  LAL_CALL( LALGetDebugLevel( &status, argc, argv, 'd'), &status);
  
  uvar_maxlag = 100000;

  uvar_help = FALSE;
  uvar_blocksRngMed = BLOCKSRNGMED;
  uvar_f0 = F0;
  uvar_startTime = 0.0;
  uvar_endTime = LAL_INT4_MAX;
  uvar_fdot = 0.0;
  uvar_fBand = FBAND;
  uvar_dAlpha = 0.2;
  uvar_dDelta = 0.2;

  uvar_earthEphemeris = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_earthEphemeris,EARTHEPHEMERIS);

  uvar_sunEphemeris = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_sunEphemeris,SUNEPHEMERIS);

  uvar_dirnameOut = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_dirnameOut,DIROUT);

  uvar_fbasenameOut = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_fbasenameOut,BASENAMEOUT);

  uvar_skyfile = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_skyfile,SKYFILE);

  /* register user input variables */
  LAL_CALL( LALRegisterBOOLUserVar( &status, "help",             'h', UVAR_HELP,     "Print this message", &uvar_help), &status);  
  LAL_CALL( LALRegisterREALUserVar( &status, "f0",               'f', UVAR_OPTIONAL, "Start search frequency", &uvar_f0), &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "fdot",               0, UVAR_OPTIONAL, "Start search frequency derivative", &uvar_fdot), &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "fBand",            'b', UVAR_OPTIONAL, "Search frequency band", &uvar_fBand), &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "startTime",         0,  UVAR_OPTIONAL, "GPS start time of observation", &uvar_startTime), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "endTime",         0,  UVAR_OPTIONAL, "GPS end time of observation", &uvar_endTime), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "skyRegion",       0,  UVAR_OPTIONAL, "sky-region polygon (or 'allsky')", &uvar_skyRegion), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dAlpha",          0,  UVAR_OPTIONAL, "Sky resolution (flat/isotropic) (rad)", &uvar_dAlpha), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "dDelta",          0,  UVAR_OPTIONAL, "Sky resolution (flat/isotropic) (rad)", &uvar_dDelta), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "skyfile",         0,  UVAR_OPTIONAL, "Alternative: input skypatch file", &uvar_skyfile),     &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "earthEphemeris", 'E', UVAR_OPTIONAL, "Earth Ephemeris file",  &uvar_earthEphemeris),  &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sunEphemeris",   'S', UVAR_OPTIONAL, "Sun Ephemeris file", &uvar_sunEphemeris), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftDir",         'D', UVAR_REQUIRED, "SFT filename pattern", &uvar_sftDir), &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "maxlag",          0,  UVAR_OPTIONAL, "Maximum time lag for correlating sfts", &uvar_maxlag), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "dirnameOut",     'o', UVAR_OPTIONAL, "Output directory", &uvar_dirnameOut), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "fbasenameOut",    0,  UVAR_OPTIONAL, "Output file basename", &uvar_fbasenameOut), &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed",    0,  UVAR_OPTIONAL, "Running Median block size", &uvar_blocksRngMed), &status);


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
 
    
  /* read sfts */
  {
    /* new SFT I/O data types */
    SFTCatalog *catalog = NULL;
    static SFTConstraints constraints;


    /* set detector constraint */
    constraints.detector = NULL;
    constraints.timestamps = NULL;

    if ( LALUserVarWasSet( &uvar_startTime ) ) {
      LAL_CALL ( LALFloatToGPS( &status, &startTimeGPS, &uvar_startTime), &status);
      constraints.startTime = &startTimeGPS;
    }

    if ( LALUserVarWasSet( &uvar_endTime ) ) {
      LAL_CALL ( LALFloatToGPS( &status, &endTimeGPS, &uvar_endTime), &status);
      constraints.endTime = &endTimeGPS;
    }

    
    /* get sft catalog */
    LAL_CALL( LALSFTdataFind( &status, &catalog, uvar_sftDir, &constraints), &status);
    if ( (catalog == NULL) || (catalog->length == 0) ) {
      fprintf (stderr,"Unable to match any SFTs with pattern '%s'\n", uvar_sftDir );
      exit(1);
    }
  

    /* first some sft parameters */
    deltaF = catalog->data[0].header.deltaF;  /* frequency resolution */
    timeBase= 1.0/deltaF; /* coherent integration time */
    f0Bin = floor( uvar_f0 * timeBase + 0.5); /* initial search frequency */
    numSearchBins =  uvar_fBand * timeBase; /* total number of search bins - 1 */
    fLastBin = f0Bin + numSearchBins;   /* final frequency bin to be analyzed */

    /* read sft files making sure to add extra bins for running median */
    /* add wings for Doppler modulation and running median block size*/
    doppWings = (uvar_f0 + uvar_fBand) * VTOT;    
    fMin = uvar_f0 - doppWings - uvar_blocksRngMed * deltaF;
    fMax = uvar_f0 + uvar_fBand + doppWings + uvar_blocksRngMed * deltaF;
    
    /* read the sfts */
    /* first load them into a MultiSFTVector, then concatenate the various vectors into one*/
    LAL_CALL( LALLoadMultiSFTs ( &status, &multiSFTs, catalog, fMin, fMax), &status);
    LAL_CALL (LoadAllSFTs (&status, &inputSFTs, multiSFTs, catalog->length), &status);
    numifo = inputSFTs->length;    

printf("there are %i sfts\n",catalog->length);
    for (j=0; j < catalog->length; j++){
printf("sft epochs %i\n", inputSFTs->data[j].epoch.gpsSeconds);
    } 


    /* find number of sfts */
    /* loop over ifos and calculate number of sfts */
    /* note that we can't use the catalog to determine the number of SFTs
       because SFTs might be segmented in frequency */
    numsft = inputSFTs->length; /* initialization */


    /* catalog is ordered in time so we can get start, end time and tObs*/
    firstTimeStamp = catalog->data[0].header.epoch;
    lastTimeStamp = catalog->data[catalog->length - 1].header.epoch;
    tObs = XLALGPSDiff( &lastTimeStamp, &firstTimeStamp ) + timeBase;

printf("observation time is %f\n",tObs);

    /* SFT info -- assume all SFTs have same length */
    binsSFT = inputSFTs->data->data->length;

    LAL_CALL(LALDestroyMultiSFTVector(&status, &multiSFTs), &status);

    LAL_CALL( LALDestroySFTCatalog( &status, &catalog ), &status);  	

  } /* end of sft reading block */


    
  /*  set up ephemeris  */
  edat = (EphemerisData *)LALCalloc(1, sizeof(EphemerisData));
  (*edat).ephiles.earthEphemeris = uvar_earthEphemeris;
  (*edat).ephiles.sunEphemeris = uvar_sunEphemeris;

  LAL_CALL( LALLeapSecs(&status, &tmpLeap, &firstTimeStamp, &lsfas), &status);
  (*edat).leap = (INT2)tmpLeap;

  LAL_CALL( LALInitBarycenter( &status, edat), &status);
 

  /* set up skypatches */
  if ((skytest = fopen(uvar_skyfile, "r")) == NULL) {
  fprintf(stderr, "skyfile doesn't exist\n");
  }


  LAL_CALL( SetUpRadiometerSkyPatches( &status, &skyInfo, uvar_skyfile, uvar_skyRegion, uvar_dAlpha, uvar_dDelta), &status);
  nSkyPatches = skyInfo.numSkyPatches;
  skyAlpha = skyInfo.alpha;
  skyDelta = skyInfo.delta;
  skySizeAlpha = skyInfo.alphaSize;
  skySizeDelta = skyInfo.deltaSize;


  /* normalize sfts */
  /* get information about all detectors including velocity and timestamps */
  /* note that this function returns the velocity at the 
     mid-time of the SFTs -- should not make any difference */
  psdVec = (PSDVector *) LALCalloc(1, sizeof(PSDVector)); 
  psdVec->length = inputSFTs->length;

  psdVec->data = (REAL8FrequencySeries *)LALCalloc (psdVec->length, sizeof(REAL8FrequencySeries));


  mdetStates = (DetectorStateSeries *) LALCalloc(1, sizeof(DetectorStateSeries));
  mdetStates->length = inputSFTs->length;
  mdetStates->data = (DetectorState *) LALCalloc (mdetStates->length, sizeof(DetectorState));

  ts = XLALCreateTimestampVector (1);
  tOffs = 0.5/deltaF;

  /*loop over all sfts and get the PSDs and detector states */
  for (j=0; j < (INT4)inputSFTs->length; j++) {
	tmpDetState = NULL;
  	psdVec->data[j].data = NULL;
	psdVec->data[j].data = (REAL8Sequence *)LALCalloc (1,sizeof(REAL8Sequence));
	psdVec->data[j].data->length = inputSFTs->data[j].data->length;
	psdVec->data[j].data->data = (REAL8 *)LALCalloc (inputSFTs->data[j].data->length, sizeof(REAL8));

	LAL_CALL( LALNormalizeSFT (&status, psdVec->data + j, inputSFTs->data + j, uvar_blocksRngMed), &status);

	det = XLALGetSiteInfo (inputSFTs->data[j].name); 
	ts->data[0] = inputSFTs->data[j].epoch;

	LAL_CALL ( LALGetDetectorStates ( &status, &tmpDetState, ts, det, edat, tOffs), &status);

printf("epoch %i\n",tmpDetState->data[0].tGPS.gpsSeconds);
  	mdetStates->data[j] = tmpDetState->data[0];

	
  }


  pairParams.lag = uvar_maxlag;  

  /* create sft pair indices */

  LAL_CALL ( CreateSFTPairsIndicesFrom2SFTvectors( &status, &sftPairIndexList, inputSFTs, &pairParams), &status);
printf("number of pairs %i\n",sftPairIndexList->vectorLength);

  /* initialise F_+, F_x vectors */

  Fplus = XLALCreateREAL8Vector(numsft);
  Fcross = XLALCreateREAL8Vector(numsft);
  frequencyShiftList = XLALCreateREAL8Vector(numsft);
  signalPhaseList = XLALCreateREAL8Vector(numsft);
  weights = XLALCreateREAL8Vector(numsft);

/*
 *initialise frequency and phase vectors *
  initialise Y, u, weight vectors */
    yalpha = (COMPLEX16Vector *) LALCalloc(1, sizeof(COMPLEX16));
    yalpha->length = sftPairIndexList->vectorLength;
    yalpha->data = LALCalloc(yalpha->length, sizeof(COMPLEX16));
    ualpha = (COMPLEX16Vector *) LALCalloc(1, sizeof(COMPLEX16));
    ualpha->length = sftPairIndexList->vectorLength;
    ualpha->data = LALCalloc(ualpha->length, sizeof(COMPLEX16));

   AMcoef = (AMCoeffs *)LALCalloc(1, sizeof(AMCoeffs));
   AMcoef->a = XLALCreateREAL4Vector(numsft);
   AMcoef->b = XLALCreateREAL4Vector(numsft);


  /*    loop over sky patches -- main calculations  */
    for (skyCounter = 0; skyCounter < nSkyPatches; skyCounter++) 
       { 

  /* initialize Doppler parameters of the potential source */
     thisPoint.Alpha = skyAlpha[skyCounter]; 
     thisPoint.Delta = skyDelta[skyCounter]; 
     thisPoint.fkdot[0] = uvar_f0;
     thisPoint.fkdot[1] = uvar_fdot; 
     thisPoint.refTime = inputSFTs->data->epoch; /*must be changed!*/


printf("coordinate is %f, %f\n",thisPoint.Alpha, thisPoint.Delta);
 
  thisVel.length = 3;
  thisPos.length = 3;
  
  /*     set sky positions and skypatch sizes  */
         patchSizeX = skySizeDelta[skyCounter]; 
         patchSizeY = skySizeAlpha[skyCounter]; 
  
  
  /*        get the amplitude modulation coefficients */
         skypos.longitude = thisPoint.Alpha; 
         skypos.latitude = thisPoint.Delta; 
         skypos.system = COORDINATESYSTEM_EQUATORIAL; 
         LAL_CALL ( LALGetAMCoeffs ( &status, AMcoef, mdetStates, skypos), &status); 

printf("A, B %f %f\n",AMcoef->A, AMcoef->B);


  /* loop over each SFT to get frequency, then store in frequencyShiftList */


  counter = 0;
  for (j=0; j < (INT4)numsft; j++) {
      
	 thisVel.data = mdetStates->data[j].vDetector;
 	 thisPos.data = mdetStates->data[j].rDetector;
	 sft = &(inputSFTs->data[j]);

         LAL_CALL( GetSignalFrequencyInSFT( &status, freq1, sft, &thisPoint, &thisVel), &status);
         LAL_CALL( GetSignalPhaseInSFT( &status, phase1, sft, &thisPoint, &thisPos), &status);

	 frequencyShiftList->data[j] = *freq1; 
	 signalPhaseList->data[j] = *phase1;
/*printf("shifted frequencies %f\n",frequencyShiftList.data[counter]);
printf("signal phases %f\n", signalPhaseList.data[counter]);*/

	Fplus->data[j] = (AMcoef->a->data[j] * cos(2.0*psi))
		 	     + (AMcoef->b->data[j] * sin(2.0*psi));

	Fcross->data[j] = (AMcoef->b->data[j] * cos(2.0*psi))
			     - (AMcoef->a->data[j] * sin(2.0*psi));

/*printf("fplus, fcross %f %f\n", Fplus->data[j], Fcross->data[j]);*/

      } 
    
     
  /* loop over SFT pairs */
    counter = 0;

    for (j=0; j < (INT4)sftPairIndexList->vectorLength; j++) {
  
     /*  correlate sft pairs  */
	index1 = sftPairIndexList->data[j];
	index2 = sftPairIndexList->data[j + sftPairIndexList->vectorLength];

	sft1 = &(inputSFTs->data[index1]);
	sft2 = &(inputSFTs->data[index2]);
	psd1 = &(psdVec->data[index1]);
	psd2 = &(psdVec->data[index2]);
 /*printf("f0_1 %f f0_2 %f\n",sft1->f0, sft2->f0);*/

  	LAL_CALL( CorrelateSingleSFTPair( &status, &(yalpha->data[j]), sft1, sft2, &(frequencyShiftList->data[index1]), &(frequencyShiftList->data[index2])), &status);


 	LAL_CALL( CalculateUalpha (&status, &ualpha->data[j], &Aplus, &Across, &signalPhaseList->data[index1], &signalPhaseList->data[index2], &Fplus->data[index1], &Fplus->data[index2], &Fcross->data[index1], &Fcross->data[index2], &frequencyShiftList->data[index1], &frequencyShiftList->data[index2], psd1, psd2), &status);

/*printf("Y %i real %1.11f imaginary %f\n", j, yalpha->data[j].re, yalpha->data[j].im);  */
/*printf("U %i real %f imaginary %f\n", j, ualpha->data[j].re, ualpha->data[j].im);  */

    }

  /* calculate rho (weights) from Yalpha and Ualpha */
	LAL_CALL( CalculateWeights( &status, weights, yalpha, ualpha), &status);

  /*   select candidates  */

    counter = 0;
    paircounter = 0;
printf("Final array of candidate pairs:\n");
    for (j=0; j < (INT4)weights->length; j++) {

     if (weights->data[j] >= rho0) {
printf("\tPair index: %i %i, weight: %f\n", sftPairIndexList->data[paircounter], sftPairIndexList->data[paircounter + sftPairIndexList->vectorLength], weights->data[j]);

     }
    paircounter++;
    }





  } /* finish loop over skypatches */ 
  
  
  /* free memory */
  XLALDestroyTimestampVector(ts);

  XLALDestroyAMCoeffs ( AMcoef ); 
  XLALDestroyCOMPLEX16Vector(yalpha);
  XLALDestroyCOMPLEX16Vector(ualpha);

  LAL_CALL (LALDestroySFTVector(&status, &inputSFTs), &status );
/*  LALFree(sftPairs.data);*/
  LALFree(sftPairIndexList->data);

  XLALDestroyDetectorStateSeries ( mdetStates );
  LAL_CALL ( LALDestroyPSDVector  ( &status, &psdVec), &status);

  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);

  LALFree(skyAlpha);
  LALFree(skyDelta);
  LALFree(skySizeAlpha);
  LALFree(skySizeDelta);

  LALFree(frequencyShiftList->data);
  LALFree(Fcross->data);
  LALFree(Fplus->data);
  LALFree(weights->data);


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

  ASSERT (out, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (dAlpha > 0, status, RADIOMETER_EBAD, RADIOMETER_MSGEBAD);
  ASSERT (dDelta > 0, status, RADIOMETER_EBAD, RADIOMETER_MSGEBAD);

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

      ASSERT (skyFileName, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
      
      if ( (fpsky = fopen(skyFileName, "r")) == NULL)
	{
	  ABORT ( status, RADIOMETER_EFILE, RADIOMETER_MSGEFILE );
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

