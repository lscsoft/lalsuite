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
#include <pulsar_crosscorr.h>
#include <lal/PulsarCrossCorr.h>
#include <lal/DopplerScan.h>
#include <gsl/gsl_permutation.h>

RCSID( "$Id$");

/* globals, constants and defaults */

extern int lalDebugLevel;

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
  REAL8 deltaF, timeBase;
  REAL8  psi; 
  UINT4 counter = 0; 
  COMPLEX8FrequencySeries *sft1 = NULL, *sft2 = NULL;
  INT4 index1, index2;
  REAL8FrequencySeries *psd1, *psd2;
  COMPLEX16Vector *yalpha = NULL, *ualpha = NULL;

  CrossCorrAmps amplitudes;
  CrossCorrBeamFn *beamfns1, *beamfns2;
  CrossCorrBeamFnListElement *beamList, *beamHead = NULL, *beamTail = NULL;

  SFTPairParams pairParams;

  /* information about all the ifos */



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
  SkyPatchesInfo skyInfo;
  INT4  j,i;

  /* frequency loop info */
  INT4 nfreqLoops, freqCounter = 0;
  INT4 nParams = 0;
  REAL8 f_current = 0.0;
  INT4 ualphacounter=0.0;

  /* frequency derivative loop info. we can go up to f_doubledot */
  INT4 nfdotLoops = 1, fdotCounter = 0;
  INT4 nfddotLoops = 1, fddotCounter = 0;
  REAL8 fdot_current = 0.0, delta_fdot = 0.0;
  REAL8 fddot_current = 0.0, delta_fddot = 0.0; 

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
  strcat(filename, "Radiometer_out.txt");

  if ((fp = fopen(filename, "w")) == NULL) {
    fprintf(stderr, "error opening output file\n");
    exit(1);
  }

  fprintf(fp, "##Alpha\tDelta\tFrequency\t Fdot \t Fddot \t Normalised Power\n");

  /* read sfts */

  /* set detector constraint */
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

  if (LALUserVarWasSet ( &uvar_refTime )) {
    XLALGPSSetREAL8(&refTime, uvar_refTime);
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
  tOffs = 0.5/deltaF;


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

  /* XLALGPSLeapSeconds (tmpLeap);*/
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

  LAL_CALL( SetUpRadiometerSkyPatches( &status, &skyInfo,
				       uvar_skyfile, uvar_skyRegion,
				       uvar_dAlpha, uvar_dDelta),
	    &status);
  nSkyPatches = skyInfo.numSkyPatches;
  skyAlpha = skyInfo.alpha;
  skyDelta = skyInfo.delta;
  skySizeAlpha = skyInfo.alphaSize;
  skySizeDelta = skyInfo.deltaSize;


  /* initialise output arrays */
  nParams = nSkyPatches * nfreqLoops *nfdotLoops * nfddotLoops;

  rho = XLALCreateREAL8Vector(nParams);
  stddev = XLALCreateREAL8Vector(nParams);

  for (j=0; j < nParams; j++) {
    rho->data[j] = 0.0;
    stddev->data[j] = 0.0;
  }


  /* set the max. allowable time lag between 2 sfts */
  pairParams.lag = uvar_maxlag;

  /* add wings for Doppler modulation and running median block size */
  /* remove fBand from doppWings because we are going bin-by-bin (?) */
  doppWings = (uvar_f0 + (freqCounter*deltaF)) * VTOT;
  fMin = uvar_f0 - doppWings - uvar_blocksRngMed * deltaF;
  fMax = uvar_f0 + uvar_fBand + doppWings + uvar_blocksRngMed * deltaF;


  slidingcounter = 0;
 	   
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
	  (XLALGPSDiff(&catalog->data[slidingcounter].header.epoch, &catalog->data[sftcounter].header.epoch) <= pairParams.lag)) {

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
      /* create sft pair indices */
      LAL_CALL ( LALCreateSFTPairsIndicesFrom2SFTvectors( &status,
							  &sftPairIndexList,
							  sftHead,
							  &pairParams,
							  listLength,
							  uvar_detChoice),
		 &status);

      /* initialise Y, u, sigmasq vectors  */
      yalpha = XLALCreateCOMPLEX16Vector(sftPairIndexList->vectorLength);
      ualpha = XLALCreateCOMPLEX16Vector(sftPairIndexList->vectorLength);
      sigmasq = XLALCreateREAL8Vector(sftPairIndexList->vectorLength);

 
      /* start frequency loop */
      for (freqCounter = 0; freqCounter < nfreqLoops; freqCounter++) {

	f_current = uvar_f0 + (deltaF*freqCounter);

	/* frequency derivative loop */
	for (fdotCounter = 0; fdotCounter < nfdotLoops; fdotCounter++) {

	  fdot_current = uvar_fdot + (delta_fdot*fdotCounter);

	  /* frequency double derivative loop */
	  for (fddotCounter = 0; fddotCounter < nfddotLoops; fddotCounter++) {

	    fddot_current = uvar_fddot + (delta_fddot*fddotCounter++);


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
				     edat, &thisPoint, psi), &status);

	      /* loop over SFT pairs */
	      ualphacounter = 0;

	      for (j=0; j < (INT4)sftPairIndexList->vectorLength; j++) {


		/*  correlate sft pairs  */
		index1 = sftPairIndexList->data[j]; /*this is always 0?*/
		index2 = sftPairIndexList->data[j + sftPairIndexList->vectorLength];

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

		for (i = 1; i <= index2; i++) {
		  /*use sftlist, psdlist as tmps */
		  sftList = (SFTListElement *)sftList->nextSFT;
		  psdList = (PSDListElement *)psdList->nextPSD;
		  freqList = (REAL8ListElement *)freqList->nextVal;
		  phaseList = (REAL8ListElement *)phaseList->nextVal;
		  beamList = (CrossCorrBeamFnListElement *)beamList->nextBeamfn;
		}

		sft2 = &(sftList->sft);
		psd2 = &(psdList->psd);
		freq2 = freqList->val;
		phase2 = phaseList->val;
		beamfns2 = &(beamList->beamfn);
		
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
						    &sigmasq->data[ualphacounter]),
			    &status);

		} else {
		  LAL_CALL( LALCalculateUalpha ( &status, &ualpha->data[ualphacounter], amplitudes,
						 phase1, phase2, *beamfns1, *beamfns2,
						 &sigmasq->data[ualphacounter]),
			    &status);
		}
		/*printf("yalpha, sigmaalpha, ualpha %e\t %e\t %e\t %e\t %e\n", yalpha->data[ualphacounter].re, yalpha->data[ualphacounter].im, sigmasq->data[ualphacounter], ualpha->data[ualphacounter].re, ualpha->data[ualphacounter].im);*/

		ualphacounter++;
	      } /*finish loop over sft pairs*/

	      /*printf("finish loop over pairs\n");*/


	      /* calculate rho from Yalpha and Ualpha */
	      tmpstat = 0;
	      LAL_CALL( LALCalculateCrossCorrPower( &status, &tmpstat, yalpha, ualpha),
			&status);

	      rho->data[counter] += tmpstat;

	      /* calculate standard deviation of rho (Eq 4.6) */
	      LAL_CALL( LALNormaliseCrossCorrPower( &status, &tmpstat, ualpha, sigmasq),
			&status); 

	      stddev->data[counter] += tmpstat;

	      /*if (counter == 0) {printf("%e \n", stddev->data[counter]);}*/

	      counter++;

	    } /* finish loop over skypatches */ 


	  } /*finish loop over frequency double derivatives*/

	} /*finish loop over frequency derivatives*/
	/*   printf("Frequency %f\n", f_current);*/
      } /*finish loop over frequencies */

      XLALDestroyCOMPLEX16Vector(yalpha);
      XLALDestroyCOMPLEX16Vector(ualpha);

      XLALDestroyREAL8Vector(sigmasq);



      XLALFree(sftPairIndexList->data);
      XLALFree(sftPairIndexList);

    } /*end if listlength > 1*/
    /*printf("finish loop over all frequencies, sftcounter %d\n", sftcounter);*/

  } /* finish loop over all sfts */
  printf("finish loop over all sfts\n");

  counter = 0;

  /* print all variables to file */
  for (freqCounter = 0; freqCounter < nfreqLoops; freqCounter++) {

    f_current = uvar_f0 + (deltaF*freqCounter);

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
		 PulsarDopplerParams *thisPoint,
		 REAL8 psi){

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

    LALGetAMCoeffs ( status->statusPtr, AMcoef, detState, skypos);
		     
    thisVel.data = detState->data[0].vDetector;
    thisPos.data = detState->data[0].rDetector;


    LALGetSignalFrequencyInSFT( status->statusPtr, &freq1, epoch, thisPoint,
				&thisVel);

    LALGetSignalPhaseInSFT( status->statusPtr, &phase1, epoch, thisPoint,
			    &thisPos);

    freqtmp->val = freq1; 
    phasetmp->val = phase1;

    /* There is some ambiguity here. If uvar_averagePsi = true,
       then we can't directly calculate the <F F> products
       (because <F_+> = <F_x> = 0) so we store a and b in the
       variables.

       The <F F> products will be calculated later in the
       calculateUalpha function.

       If uvar_averagePsi = false, then there is no problem and
       we calculate Fplus_or_a and Fcross_or_b here */

    if(uvar_averagePsi) {
      beamtmp->beamfn.Fplus_or_a = (AMcoef->a->data[0]);
      beamtmp->beamfn.Fcross_or_b = (AMcoef->b->data[0]);
    } else {
      beamtmp->beamfn.Fplus_or_a = (AMcoef->a->data[0] * cos(2.0*psi))
	+ (AMcoef->b->data[0] * sin(2.0*psi));
      beamtmp->beamfn.Fcross_or_b = (AMcoef->b->data[0] * cos(2.0*psi))
	- (AMcoef->a->data[0] * sin(2.0*psi));
    }
		
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
			INT4 index) 
{
  SFTCatalog *slidingcat;
  SFTDescriptor *desc;

  INITSTATUS (status, "CopySFTFromCatalog", rcsid);
  ATTATCHSTATUSPTR (status);


  slidingcat = NULL;
  *sft = NULL;

  ASSERT ( catalog, status, PULSAR_CROSSCORR_ENULL, PULSAR_CROSSCORR_MSGENULL );


  if ( (slidingcat = LALCalloc ( 1, sizeof ( SFTCatalog ))) == NULL ) {
    ABORT(status, PULSAR_CROSSCORR_EMEM, PULSAR_CROSSCORR_MSGEMEM);
  }

  slidingcat->length = 1;
  slidingcat->data = LALRealloc(slidingcat->data, 1*sizeof(*(slidingcat->data)));
  /* memset(&(slidingcat->data[0]), 0, sizeof(ret->data[0]));*/
  desc = &(slidingcat->data[0]);
  /*  desc->locator = LALCalloc(1, sizeof(*(desc->locator)));*/

  desc->locator = catalog->data[index].locator;
  desc->header = catalog->data[index].header;
  desc->comment = catalog->data[index].comment;
  desc->numBins = catalog->data[index].numBins;
  desc->version = catalog->data[index].version;
  desc->crc64 = catalog->data[index].crc64;

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
  beamList->beamfn.Fplus_or_a = 0;
  beamList->beamfn.Fcross_or_b = 0;
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

void initUserVars (LALStatus *status)
{
  INITSTATUS( status, "initUserVars", rcsid );
  ATTATCHSTATUSPTR (status);


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
  LALRegisterREALUserVar( status->statusPtr, "fdotBand",
			  0, UVAR_OPTIONAL,
			  "Search frequency derivative band",
			  &uvar_fdotBand);
  LALRegisterREALUserVar( status->statusPtr, "fdotRes",
			  'r', UVAR_OPTIONAL,
			  "Search frequency derivative resolution",
			  &uvar_fdotResolution);
  LALRegisterBOOLUserVar( status->statusPtr, "searchfddot",
			  0, UVAR_OPTIONAL,
			  "Search a range of frequency double derivative",
			  &uvar_searchfddot);
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
			  "Search frequency double derivative resolution",
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
  LALRegisterSTRINGUserVar( status->statusPtr, "fbasenameOut",
			    0, UVAR_OPTIONAL,
			    "Output file basename",
			    &uvar_fbasenameOut);
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

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


