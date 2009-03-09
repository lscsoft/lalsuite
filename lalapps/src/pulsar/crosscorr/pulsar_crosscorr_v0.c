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
 * $Id$
 *


/ lalapps includes */
#include <lalapps.h>
#include <pulsar_crosscorr_v0.h>
#include <lal/PulsarCrossCorr_v0.h>
#include <lal/DopplerScan.h>
#include <gsl/gsl_permutation.h>

RCSID( "$Id$");

/* globals, constants and defaults */

extern int lalDebugLevel;

#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
#define EPHEM_YEARS "00-04"

#define F0 100
#define FBAND 1
#define LAL_INT4_MAX 2147483647


#define BLOCKSRNGMED 51
#define MAXFILENAMELENGTH 512 /* maximum # of characters  of a filename */

#define DIROUT "./output/"   /* output directory */
#define BASENAMEOUT "radio"    /* prefix file output */

#define SKYFILE "./skypatchfile"
#define SKYREGION "allsky"

#define TRUE (1==1)
#define FALSE (1==0)


int main(int argc, char *argv[]){
  /* LALStatus pointer */
  static LALStatus status;  
  
  /* sft related variables */ 

  SFTVector *inputSFTs = NULL;
  REAL8 deltaF, timeBase, freq = 0, *freq1 = &freq;
  REAL8 phase = 0, *phase1 = &phase, psi; 
  UINT4 numsft, counter = 0; 
  COMPLEX8FrequencySeries *sft = NULL, *sft1 = NULL, *sft2 = NULL;
  INT4 index1, index2;
  REAL8FrequencySeries *psd1, *psd2;
  COMPLEX16Vector *yalpha = NULL, *ualpha = NULL;

  CrossCorrAmps_v0 amplitudes;
  CrossCorrBeamFn_v0 *beamfns;
  SFTPairParams_v0 pairParams;

  /* information about all the ifos */
  PSDVector *psdVec = NULL;  
  LALDetector *det;
 
  LIGOTimeGPS firstTimeStamp, lastTimeStamp;
  REAL8 tObs, tOffs;
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
  SkyPatchesInfo_v0 skyInfo;
  INT4  j;

  /* frequency loop info */
  INT4 nfreqLoops, freqCounter = 0;
  INT4 nParams = 0;
  REAL8 f_current = 0.0;
  REAL8 *stddev, ualphacounter=0.0;

  /* frequency derivative loop info. we can go up to f_doubledot */
  INT4 nfdotLoops = 1, fdotCounter = 0;
  INT4 nfddotLoops = 1, fddotCounter = 0;
  REAL8 fdot_current = 0.0, delta_fdot = 0.0;
  REAL8 fddot_current = 0.0, delta_fddot = 0.0; 

  static INT4VectorSequence  *sftPairIndexList;
  REAL8Vector *frequencyShiftList, *signalPhaseList, *sigmasq;
  PulsarDopplerParams thisPoint;
  static REAL8Vector thisVel, thisPos, *weights;
 

  REAL8 doppWings, fMin, fMax;

  AMCoeffs *AMcoef = NULL; 

  /* output file */
  FILE *fp = NULL;
  CHAR filename[MAXFILENAMELENGTH];

  /* sft constraint variables */
  LIGOTimeGPS startTimeGPS, endTimeGPS, refTime;
  LIGOTimeGPSVector *ts=NULL;

  /* user input variables */
  BOOLEAN  uvar_help;
  BOOLEAN  uvar_averagePsi;
  BOOLEAN  uvar_averageIota;
  BOOLEAN  uvar_searchfddot;
  INT4     uvar_blocksRngMed; 
  INT4     uvar_detChoice;
  REAL8    uvar_startTime, uvar_endTime;
  REAL8    uvar_f0, uvar_fdot, uvar_fBand, uvar_fdotBand;
  REAL8    uvar_fddot, uvar_fddotBand;
  REAL8    uvar_fdotResolution, uvar_fddotResolution;
  REAL8    uvar_dAlpha, uvar_dDelta; /* resolution for isotropic sky-grid */
  REAL8    uvar_maxlag;
  REAL8    uvar_psi;
  REAL8    uvar_refTime;
  REAL8    uvar_cosi;
  CHAR     *uvar_ephemDir=NULL;
  CHAR     *uvar_ephemYear=NULL;
  CHAR     *uvar_sftDir=NULL;
  CHAR     *uvar_dirnameOut=NULL;
  CHAR     *uvar_fbasenameOut=NULL;
  CHAR     *uvar_skyfile=NULL;
  CHAR     *uvar_skyRegion=NULL;
  FILE	   *skytest=NULL;
 
 
  SkyPosition skypos; 

  /* new SFT I/O data types */
  SFTCatalog *catalog = NULL;
  static SFTConstraints constraints;

  MultiSFTVector *multiSFTs = NULL;   

  /* LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;
  
  lalDebugLevel = 0;  /* LALDebugLevel must be called before anything else */
  LAL_CALL( LALGetDebugLevel( &status, argc, argv, 'd'), &status);
   
  uvar_maxlag = 0;
 
  uvar_help = FALSE;
  uvar_averagePsi = TRUE;
  uvar_averageIota = TRUE;
  uvar_blocksRngMed = BLOCKSRNGMED;
  uvar_searchfddot = FALSE;
  uvar_detChoice = 2;
  uvar_f0 = F0;
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
 
  uvar_ephemDir = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_ephemDir,DEFAULT_EPHEMDIR);
 
  uvar_ephemYear = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_ephemYear,EPHEM_YEARS);
 
  uvar_dirnameOut = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_dirnameOut,DIROUT);
 
  uvar_fbasenameOut = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_fbasenameOut,BASENAMEOUT);
 
  uvar_skyfile = (CHAR *)LALCalloc( MAXFILENAMELENGTH , sizeof(CHAR));
  strcpy(uvar_skyfile,SKYFILE);
 
  /* register user input variables */
  LAL_CALL( LALRegisterBOOLUserVar( &status, "help",
				    'h', UVAR_HELP,
				    "Print this message",
				    &uvar_help),
	    &status);  
  LAL_CALL( LALRegisterBOOLUserVar( &status, "averagePsi",
				    0, UVAR_OPTIONAL,
				    "Use average over psi",
				    &uvar_averagePsi),
	    &status); 
  LAL_CALL( LALRegisterBOOLUserVar( &status, "averageIota",
				    0, UVAR_OPTIONAL,
				    "Use average over iota",
				    &uvar_averageIota),
	    &status); 
  LAL_CALL( LALRegisterREALUserVar( &status, "f0",
				    'f', UVAR_OPTIONAL,
				    "Start search frequency",
				    &uvar_f0),
	    &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "fdot",
				    0, UVAR_OPTIONAL,
				    "Start search frequency derivative",
				    &uvar_fdot),
	    &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "fBand",
				    'b', UVAR_OPTIONAL,
				    "Search frequency band",
				    &uvar_fBand),
	    &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "fdotBand",
				    0, UVAR_OPTIONAL,
				    "Search frequency derivative band",
				    &uvar_fdotBand),
	    &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "fdotRes",
				    'r', UVAR_OPTIONAL,
				    "Search frequency derivative resolution",
				    &uvar_fdotResolution),
	    &status);
  LAL_CALL( LALRegisterBOOLUserVar( &status, "searchfddot",
				    0, UVAR_OPTIONAL,
				    "Search a range of frequency double derivative",
				    &uvar_searchfddot),
	    &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "fddot",
				    0, UVAR_OPTIONAL,
				    "Start frequency double derivative",
				    &uvar_fddot),
	    &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "fddotBand",
				    0, UVAR_OPTIONAL,
				    "Search frequency double derivative band",
				    &uvar_fddotBand),
	    &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "fddotRes",
				    0, UVAR_OPTIONAL,
				    "Search frequency double derivative resolution",
				    &uvar_fddotResolution),
	    &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "startTime",
				    0, UVAR_OPTIONAL,
				    "GPS start time of observation",
				    &uvar_startTime),
	    &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "endTime",
				    0, UVAR_OPTIONAL,
				    "GPS end time of observation",
				    &uvar_endTime),
	    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "skyRegion",
				      0, UVAR_OPTIONAL,
				      "sky-region polygon (or 'allsky')",
				      &uvar_skyRegion),
	    &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "dAlpha",
				    0, UVAR_OPTIONAL,
				    "Sky resolution (flat/isotropic) (rad)",
				    &uvar_dAlpha),
	    &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "dDelta",
				    0, UVAR_OPTIONAL,
				    "Sky resolution (flat/isotropic) (rad)",
				    &uvar_dDelta),
	    &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "psi",
				    0, UVAR_OPTIONAL,
				    "Polarisation angle (rad)",
				    &uvar_psi),
	    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "skyfile",
				      0, UVAR_OPTIONAL,
				      "Alternative: input skypatch file",
				      &uvar_skyfile),
	    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "ephemDir",
				      'E', UVAR_OPTIONAL,
				  "Directory where ephemeris files are located",
				      &uvar_ephemDir),
	    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "ephemYear",
				      'y', UVAR_OPTIONAL,
		       "Year (or range of years) of ephemeris files to be used",
				      &uvar_ephemYear),
	    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftDir",
				      'D', UVAR_REQUIRED,
				      "SFT filename pattern",
				      &uvar_sftDir),
	    &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "maxlag",
				    0,  UVAR_OPTIONAL,
				    "Maximum time lag for correlating sfts",
				    &uvar_maxlag),
	    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "dirnameOut",
				      'o', UVAR_OPTIONAL,
				      "Output directory",
				      &uvar_dirnameOut),
	    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "fbasenameOut",
				      0, UVAR_OPTIONAL,
				      "Output file basename",
				      &uvar_fbasenameOut),
	    &status);
  LAL_CALL( LALRegisterINTUserVar( &status, "blocksRngMed",
				   0, UVAR_OPTIONAL,
				   "Running Median block size",
				   &uvar_blocksRngMed),
	    &status);
  LAL_CALL( LALRegisterINTUserVar( &status, "detChoice",
				   0, UVAR_OPTIONAL,
   "0: Correlate SFTs from same IFOs only; 1: different IFOs only; 2: all IFOs",
				   &uvar_detChoice),
	    &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "refTime",
				    0, UVAR_OPTIONAL,
				    "Pulsar reference time (SSB)",
				    &uvar_refTime),
	    &status); 
  LAL_CALL( LALRegisterREALUserVar( &status, "cosi",
				    0, UVAR_OPTIONAL,
				    "cos(iota) inclination angle",
				    &uvar_cosi),
	    &status); 
 
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

  if (uvar_fdotResolution < 0) {
    fprintf(stderr, "frequency derivative resolution must be positive\n");
    exit(1);
  }

  if (uvar_fddotResolution < 0) {
    fprintf(stderr, "frequency double derivative resolution must be positive\n");
    exit(1);
  }


  /* open output file */
  strcpy (filename, uvar_dirnameOut);
  strcat(filename, "Radiometer_out_v0.txt");

  if ((fp = fopen(filename, "w")) == NULL) {
    fprintf(stderr, "error opening output file\n");
    exit(1);
  }

  fprintf(fp, "##Alpha\tDelta\tFrequency\t Fdot \t Fddot \t Normalised Power\n");

  /* read sfts */

  /* set detector constraint */
  constraints.detector = NULL;
  constraints.timestamps = NULL;

  if ( LALUserVarWasSet( &uvar_startTime ) ) {
    LAL_CALL ( LALFloatToGPS( &status, &startTimeGPS, &uvar_startTime),
	       &status);
    constraints.startTime = &startTimeGPS;
  }

  if ( LALUserVarWasSet( &uvar_endTime ) ) {
    LAL_CALL ( LALFloatToGPS( &status, &endTimeGPS, &uvar_endTime), &status);
    constraints.endTime = &endTimeGPS;
  }

  if (LALUserVarWasSet ( &uvar_refTime )) {
    LAL_CALL(LALFloatToGPS(&status, &refTime, &uvar_refTime), &status);
  }

  /* get sft catalog */
  LAL_CALL( LALSFTdataFind( &status, &catalog, uvar_sftDir, &constraints),
	    &status);
  if ( (catalog == NULL) || (catalog->length == 0) ) {
    fprintf (stderr, "Unable to match any SFTs with pattern '%s'\n",
	     uvar_sftDir );
    exit(1);
  }

  /* first some sft parameters */
  deltaF = catalog->data[0].header.deltaF;  /* frequency resolution */
  timeBase= 1.0/deltaF; /* coherent integration time */

  /* catalog is ordered in time so we can get start, end time and tObs */
  firstTimeStamp = catalog->data[0].header.epoch;
  lastTimeStamp = catalog->data[catalog->length - 1].header.epoch;
  tObs = XLALGPSDiff( &lastTimeStamp, &firstTimeStamp ) + timeBase;

  nfreqLoops = ceil(uvar_fBand/deltaF);

  /* only loop through more than 1 fdot value if fdotBand and fdotRes
     are non-zero */
  if (uvar_fdotBand != 0 && uvar_fdotResolution > 0) {
    nfdotLoops = 1 + abs(ceil(uvar_fdotBand/uvar_fdotResolution));

    delta_fdot = (uvar_fdotBand > 0)
      ? uvar_fdotResolution : -(uvar_fdotResolution);
  }

  /* only loop through more than 1 fddot value if fddotBand and fddotRes
     are non-zero, and if we want to */
  if (uvar_searchfddot && uvar_fddotBand != 0 && uvar_fddotResolution > 0) {
    nfddotLoops = 1 + abs(ceil(uvar_fddotBand/uvar_fddotResolution));

    delta_fddot = (uvar_fddotBand > 0)
      ? uvar_fddotResolution : -(uvar_fddotResolution);
  }

  /* polarisation angle */
  psi = uvar_psi;

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


  /* curly As */
  /* because we have the option of either averaging over i or not, we
     need to calculate A_{+,x}^2 and A_xA_+ rather than the individual
     values because <A_x> = 0 */
  if (uvar_averageIota) {
    amplitudes.Aplussq = 7.0/15.0;
    amplitudes.Acrosssq = 1.0/3.0;
    amplitudes.AplusAcross = 0;
  } else {
    amplitudes.Aplussq = pow(((1.0 + uvar_cosi*uvar_cosi)/2.0),2);
    amplitudes.Acrosssq = pow(uvar_cosi,2);
    amplitudes.AplusAcross = (uvar_cosi/2) + (pow(uvar_cosi,3)/2);
  }


  /* set up skypatches */
  if ((skytest = fopen(uvar_skyfile, "r")) == NULL) {
    fprintf(stderr, "skyfile doesn't exist\n");
  }

  LAL_CALL( SetUpRadiometerSkyPatches_v0( &status, &skyInfo,
				       uvar_skyfile, uvar_skyRegion,
				       uvar_dAlpha, uvar_dDelta),
	    &status);
  nSkyPatches = skyInfo.numSkyPatches;
  skyAlpha = skyInfo.alpha;
  skyDelta = skyInfo.delta;
  skySizeAlpha = skyInfo.alphaSize;
  skySizeDelta = skyInfo.deltaSize;

  /* initialise output arrays */
  nParams = nSkyPatches * nfreqLoops *nfdotLoops;

  weights = XLALCreateREAL8Vector(nParams);

  stddev = LALCalloc(1, sizeof(REAL8));

  /* set the max. allowable time lag between 2 sfts */
  pairParams.lag = uvar_maxlag;

  /* read the sfts */
  /* first load them into a MultiSFTVector, then concatenate the
     various vectors into one */
  /* read sft files making sure to add extra bins for running median */
  /* add wings for Doppler modulation and running median block size */
  /* remove fBand from doppWings because we are going bin-by-bin (?) */
  doppWings = (uvar_f0 + (freqCounter*deltaF)) * VTOT;
  fMin = uvar_f0 - doppWings - uvar_blocksRngMed * deltaF;
  fMax = uvar_f0 + uvar_fBand + doppWings + uvar_blocksRngMed * deltaF;

  /* combine all sfts into a single sftvector */
  LAL_CALL( LALLoadMultiSFTs ( &status, &multiSFTs, catalog, fMin, fMax), &status);
  LAL_CALL( LALCombineAllSFTs_v0 (&status, &inputSFTs, multiSFTs, catalog->length), &status);
 
  /* clean up sft catalog and multisftvector */
  LAL_CALL(LALDestroyMultiSFTVector(&status, &multiSFTs), &status);
  LAL_CALL( LALDestroySFTCatalog( &status, &catalog ), &status);

  /* find number of sfts */
  numsft = inputSFTs->length;

  /* initialise PSD vector */
  psdVec = (PSDVector *) LALCalloc(1, sizeof(PSDVector));
  psdVec->length = numsft;

  psdVec->data = (REAL8FrequencySeries *)
  		  LALCalloc (psdVec->length, sizeof(REAL8FrequencySeries));

  /*loop over all sfts, get their PSDs and normalise them with LALNormalizeSFT*/
  for (j=0; j < (INT4)numsft; j++) {

    psdVec->data[j].data = NULL;
    psdVec->data[j].data = (REAL8Sequence *)
    			    LALCalloc (1,sizeof(REAL8Sequence));
    psdVec->data[j].data->length = inputSFTs->data[j].data->length;
    psdVec->data[j].data->data = (REAL8 *)
    				LALCalloc (inputSFTs->data[j].data->length, sizeof(REAL8));
	
    LAL_CALL( LALNormalizeSFT (&status, psdVec->data + j,
     			       inputSFTs->data + j, uvar_blocksRngMed),
	 	&status);
	
  }


  /* create sft pair indices */
  LAL_CALL ( LALCreateSFTPairsIndicesFrom2SFTvectors_v0( &status,
 						      &sftPairIndexList,
						      inputSFTs,
						      &pairParams,
						      uvar_detChoice),
	     &status);

  /* start frequency loop */

  for (freqCounter = 0; freqCounter < nfreqLoops; freqCounter++) {
    f_current = uvar_f0 + (deltaF*freqCounter);

    /* frequency derivative loop */
    for (fdotCounter = 0; fdotCounter < nfdotLoops; fdotCounter++) {

      fdot_current = uvar_fdot + (delta_fdot*fdotCounter);

      /* frequency double derivative loop */
      for (fddotCounter = 0; fddotCounter < nfddotLoops; fddotCounter++) {

	fddot_current = uvar_fddot + (delta_fddot*fddotCounter++);
 
        ts = XLALCreateTimestampVector (1);
        tOffs = 0.5/deltaF;

        beamfns = (CrossCorrBeamFn_v0 *) LALCalloc(numsft, sizeof(CrossCorrBeamFn_v0));

        frequencyShiftList = XLALCreateREAL8Vector(numsft);
        signalPhaseList = XLALCreateREAL8Vector(numsft);


        /* initialise frequency and phase vectors *
         * initialise Y, u, weight vectors */
        yalpha = XLALCreateCOMPLEX16Vector(sftPairIndexList->vectorLength);
        ualpha = XLALCreateCOMPLEX16Vector(sftPairIndexList->vectorLength);

        sigmasq = XLALCreateREAL8Vector(sftPairIndexList->vectorLength);

        /*printf("starting loops over sky patches\n");*/

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
 
	  thisVel.length = 3;
	  thisPos.length = 3;
  
    	  /* set sky positions and skypatch sizes  */
	  patchSizeX = skySizeDelta[skyCounter]; 
	  patchSizeY = skySizeAlpha[skyCounter]; 
  
  
    	  /* get the amplitude modulation coefficients */
	  skypos.longitude = thisPoint.Alpha; 
	  skypos.latitude = thisPoint.Delta; 
	  skypos.system = COORDINATESYSTEM_EQUATORIAL; 

	  /* loop over each SFT to get frequency, then store in
	     frequencyShiftList */
	  for (j=0; j < (INT4)numsft; j++) {
     
            /* get information about all detectors including velocity and
             timestamps */
	    /*only have 1 element in detectorStateSeries and AMCoeffs because
	      the detector has to be updated for every SFT */
	    DetectorStateSeries *detState = NULL;
	    AMcoef = (AMCoeffs *)LALCalloc(1, sizeof(AMCoeffs));
	    AMcoef->a = XLALCreateREAL4Vector(1);
	    AMcoef->b = XLALCreateREAL4Vector(1);

	    det = XLALGetSiteInfo (inputSFTs->data[j].name); 
	    ts->data[0] = inputSFTs->data[j].epoch;

            /* note that this function returns the velocity at the
               mid-time of the SFTs -- should not make any difference */

	    LAL_CALL ( LALGetDetectorStates ( &status, &detState, ts, det,
					    edat, tOffs),
		     &status);

	    detState->detector = *det;

	    LAL_CALL ( LALGetAMCoeffs ( &status, AMcoef, detState, skypos),
		     &status); 

  	    thisVel.data = detState->data[0].vDetector;
	    thisPos.data = detState->data[0].rDetector;
	    sft = &(inputSFTs->data[j]);


	    LAL_CALL( LALGetSignalFrequencyInSFT_v0( &status, freq1, sft, &thisPoint,
						&thisVel),
		    &status);
	    LAL_CALL( LALGetSignalPhaseInSFT_v0( &status, phase1, sft, &thisPoint,
	 				    &thisPos),
	  	    &status);

	    frequencyShiftList->data[j] = *freq1; 
	    signalPhaseList->data[j] = *phase1;

	    /* There is some ambiguity here. If uvar_averagePsi = true,
	       then we can't directly calculate the <F F> products
	       (because <F_+> = <F_x> = 0) so we store a and b in the
	       variables.

	       The <F F> products will be calculated later in the
	       calculateUalpha function.

	       If uvar_averagePsi = false, then there is no problem and
	       we calculate Fplus_or_a and Fcross_or_b here */

	    if(uvar_averagePsi) {
	      beamfns[j].Fplus_or_a = (AMcoef->a->data[0]);
	      beamfns[j].Fcross_or_b = (AMcoef->b->data[0]);
	    } else {
	      beamfns[j].Fplus_or_a = (AMcoef->a->data[0] * cos(2.0*psi))
	      + (AMcoef->b->data[0] * sin(2.0*psi));

	      beamfns[j].Fcross_or_b = (AMcoef->b->data[0] * cos(2.0*psi))
	      - (AMcoef->a->data[0] * sin(2.0*psi));
	    }

	    sft = NULL;
	
	    /* clean up AMcoefs */
	    XLALDestroyAMCoeffs(AMcoef);
	    XLALDestroyDetectorStateSeries(detState);
	    XLALFree(det);

	  } /*finish loop over individual sfts */

	  /* loop over SFT pairs */

	  ualphacounter = 0.0;
	  for (j=0; j < (INT4)sftPairIndexList->vectorLength; j++) {
  
	  /*  correlate sft pairs  */
	  index1 = sftPairIndexList->data[j];
	  index2 = sftPairIndexList->data[j + sftPairIndexList->vectorLength];

	  sft1 = &(inputSFTs->data[index1]);
	  sft2 = &(inputSFTs->data[index2]);

	  psd1 = &(psdVec->data[index1]);
	  psd2 = &(psdVec->data[index2]);
  	
	  LAL_CALL( LALCorrelateSingleSFTPair_v0( &status, &(yalpha->data[j]),
					       sft1, sft2, psd1, psd2,
					    &(frequencyShiftList->data[index1]),
					   &(frequencyShiftList->data[index2])),
		    &status);

	  LAL_CALL( LALCalculateSigmaAlphaSq_v0( &status, &sigmasq->data[j],
					      frequencyShiftList->data[index1],
					      frequencyShiftList->data[index2],
					      psd1, psd2),
		    &status);

	  /*if we are averaging over psi and cos(iota), call the simplified 
 	    Ualpha function*/
	  if (uvar_averagePsi && uvar_averageIota) {
	    LAL_CALL( LALCalculateAveUalpha_v0 ( &status, &ualpha->data[j], 
		 			   &signalPhaseList->data[index1],
					   &signalPhaseList->data[index2],
					   beamfns[index1],
					   beamfns[index2], &sigmasq->data[j]),
		      &status);
	  }
	  else {
	    LAL_CALL( LALCalculateUalpha_v0 ( &status, &ualpha->data[j], amplitudes,
		 			   &signalPhaseList->data[index1],
					   &signalPhaseList->data[index2],
					   beamfns[index1],
					   beamfns[index2], &sigmasq->data[j]),
		      &status);
	  }

	  ualphacounter = ualphacounter + 1;
  	  } /*finish loop over sft pairs*/


	  /* calculate rho (weights) from Yalpha and Ualpha */
	  LAL_CALL( LALCalculateCrossCorrPower_v0( &status,
					      &(weights->data[counter]),
					      yalpha, ualpha),
		  &status);

	  /* calculate standard deviation of rho (Eq 4.6) */
	  LAL_CALL( LALNormaliseCrossCorrPower_v0( &status, stddev, ualpha, sigmasq),
		  &status); 

	  weights->data[counter] = weights->data[counter]/(*stddev);

	  /* select candidates  */
	  /* print all variables to file */

	  fprintf(fp, "%1.5f\t %1.5f\t %1.5f\t %e\t %e\t %1.10f\n", thisPoint.Alpha,
	  	thisPoint.Delta, f_current,
		fdot_current, fddot_current, weights->data[counter]);

 	  counter++;

        } /* finish loop over skypatches */ 
      printf("Frequency %f\n", f_current);

      /* free memory */


      XLALDestroyTimestampVector(ts);

      XLALDestroyCOMPLEX16Vector(yalpha);
      XLALDestroyCOMPLEX16Vector(ualpha);
      XLALDestroyREAL8Vector(frequencyShiftList);
      XLALDestroyREAL8Vector(signalPhaseList);

      XLALDestroyREAL8Vector(sigmasq);
      XLALFree(beamfns); 

      /*  XLALDestroyDetectorStateSeries ( detStates );*/

      } /*finish loop over frequency double derivatives*/

    } /*finish loop over frequency derivatives*/

  } /*finish loop over frequencies */

  fclose (fp);

  LAL_CALL (LALDestroySFTVector(&status, &inputSFTs), &status );
  LAL_CALL ( LALDestroyPSDVector  ( &status, &psdVec), &status);
  XLALFree(sftPairIndexList->data);
  XLALFree(sftPairIndexList);

  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);
  LALFree(skyAlpha);
  LALFree(skyDelta);
  LALFree(skySizeAlpha);
  LALFree(skySizeDelta);
  LALFree(stddev);

  XLALDestroyREAL8Vector(weights);

  LAL_CALL (LALDestroyUserVars(&status), &status);

  LALCheckMemoryLeaks();

  if ( lalDebugLevel )
    REPORTSTATUS ( &status);

  return status.statusCode;

} /* main */


/** Set up location of skypatch centers and sizes
    If user specified skyRegion then use DopplerScan function
    to construct an isotropic grid. Otherwise use skypatch file. */
void SetUpRadiometerSkyPatches_v0(LALStatus           *status,
			       SkyPatchesInfo_v0      *out,   /**< output skypatches info */
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

  ASSERT (out, status, PULSAR_CROSSCORR_V0_ENULL, PULSAR_CROSSCORR_V0_MSGENULL);
  ASSERT (dAlpha > 0, status, PULSAR_CROSSCORR_V0_EBAD, PULSAR_CROSSCORR_V0_MSGEBAD);
  ASSERT (dDelta > 0, status, PULSAR_CROSSCORR_V0_EBAD, PULSAR_CROSSCORR_V0_MSGEBAD);

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

      ASSERT (skyFileName, status, PULSAR_CROSSCORR_V0_ENULL, PULSAR_CROSSCORR_V0_MSGENULL);

      if ( (fpsky = fopen(skyFileName, "r")) == NULL) {
	ABORT ( status, PULSAR_CROSSCORR_V0_EFILE, PULSAR_CROSSCORR_V0_MSGEFILE );
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
