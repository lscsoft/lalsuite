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
   REAL8 deltaF, timeBase, freq = 0, *freq1 = &freq;
   REAL8 phase = 0, *phase1 = &phase, psi; 
   UINT4 numsft, counter = 0; 
   COMPLEX8FrequencySeries *sft = NULL, *sft1 = NULL, *sft2 = NULL;
   INT4 index1, index2;
   REAL8FrequencySeries *psd1, *psd2;
   COMPLEX16Vector *yalpha = NULL, *ualpha = NULL;
   REAL8Vector *Fplus, *Fcross;
   REAL8 Aplus, Across; 
   /*A's are averaged over cos i... need to change this later!*/ 

   SFTPairParams pairParams;
   
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

   /* skypatch info */
   REAL8  *skyAlpha=NULL, *skyDelta=NULL, *skySizeAlpha=NULL, *skySizeDelta=NULL; 
   INT4   nSkyPatches, skyCounter=0; 
   SkyPatchesInfo skyInfo;
   INT4  j;

   /* frequency loop info */
   INT4 nfreqLoops, freqCounter = 0;
   INT4 nParams = 0;
   REAL8 *stddev, ualphacounter=0.0;

   static INT4VectorSequence  *sftPairIndexList;
   REAL8Vector *frequencyShiftList, *signalPhaseList, *sigmasq;
   PulsarDopplerParams  thisPoint;
   static REAL8Vector thisVel, thisPos, *weights;
 

   REAL8 doppWings, fMin, fMax;

   AMCoeffs *AMcoef = NULL; 

   /* output file */
   FILE *fp = NULL;

   /* sft constraint variables */
   LIGOTimeGPS startTimeGPS, endTimeGPS, refTime;
   LIGOTimeGPSVector *ts=NULL;

   /* user input variables */
   BOOLEAN  uvar_help;
   INT4     uvar_blocksRngMed; 
   INT4     uvar_detChoice;
   REAL8    uvar_startTime, uvar_endTime;
   REAL8    uvar_f0, uvar_fdot, uvar_fBand;
   REAL8    uvar_dAlpha, uvar_dDelta; /* resolution for isotropic sky-grid */
   REAL8    uvar_maxlag;
   REAL8    uvar_psi;
   REAL8    uvar_refTime;
   REAL8    uvar_cosi;
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
   uvar_detChoice = 2;
   uvar_f0 = F0;
   uvar_startTime = 0.0;
   uvar_endTime = LAL_INT4_MAX;
   uvar_fdot = 0.0;
   uvar_fBand = FBAND;
   uvar_dAlpha = 0.2;
   uvar_dDelta = 0.2;
   uvar_psi = 0.0;
   uvar_refTime = 0.0;
   uvar_cosi = 0.0;
 
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
   LAL_CALL( LALRegisterREALUserVar(   &status, "psi",             0,  UVAR_OPTIONAL, "Polarisation angle (rad)", &uvar_psi), &status);
   LAL_CALL( LALRegisterSTRINGUserVar( &status, "skyfile",         0,  UVAR_OPTIONAL, "Alternative: input skypatch file", &uvar_skyfile),     &status);
   LAL_CALL( LALRegisterSTRINGUserVar( &status, "earthEphemeris", 'E', UVAR_OPTIONAL, "Earth Ephemeris file",  &uvar_earthEphemeris),  &status);
   LAL_CALL( LALRegisterSTRINGUserVar( &status, "sunEphemeris",   'S', UVAR_OPTIONAL, "Sun Ephemeris file", &uvar_sunEphemeris), &status);
   LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftDir",         'D', UVAR_REQUIRED, "SFT filename pattern", &uvar_sftDir), &status);
   LAL_CALL( LALRegisterREALUserVar(   &status, "maxlag",          0,  UVAR_OPTIONAL, "Maximum time lag for correlating sfts", &uvar_maxlag), &status);
   LAL_CALL( LALRegisterSTRINGUserVar( &status, "dirnameOut",     'o', UVAR_OPTIONAL, "Output directory", &uvar_dirnameOut), &status);
   LAL_CALL( LALRegisterSTRINGUserVar( &status, "fbasenameOut",    0,  UVAR_OPTIONAL, "Output file basename", &uvar_fbasenameOut), &status);
   LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed",    0,  UVAR_OPTIONAL, "Running Median block size", &uvar_blocksRngMed), &status);
   LAL_CALL( LALRegisterINTUserVar(    &status, "detChoice",       0,  UVAR_OPTIONAL, "0: Correlate SFTs from same IFOs only; 1: different IFOs only; 2: all IFOs", &uvar_detChoice), &status);
   LAL_CALL( LALRegisterREALUserVar(    &status, "refTime",    	   0,  UVAR_OPTIONAL, "Pulsar reference time (SSB)", &uvar_refTime), &status); 
   LAL_CALL( LALRegisterREALUserVar(    &status, "cosi",    	   0,  UVAR_OPTIONAL, "cos(iota) inclination angle", &uvar_cosi), &status); 
 
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

   /* open output file */
   if ((fp = fopen("Radiometer_out.txt", "w")) == NULL) {
    fprintf(stderr, "error opening output file\n");
    exit(1);
   }
 
   fprintf(fp, "##Alpha\tDelta\tFrequency\t\tRaw Power \t\t Sigma \t\t Normalised Power\n");
    
   /* read sfts */
   
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

   if (LALUserVarWasSet ( &uvar_refTime )) {
	LAL_CALL(LALFloatToGPS(&status, &refTime, &uvar_refTime), &status);
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



   /* catalog is ordered in time so we can get start, end time and tObs*/
   firstTimeStamp = catalog->data[0].header.epoch;
   lastTimeStamp = catalog->data[catalog->length - 1].header.epoch;
   tObs = XLALGPSDiff( &lastTimeStamp, &firstTimeStamp ) + timeBase;

   nfreqLoops = ceil(uvar_fBand/deltaF);


   /*  set up ephemeris  */
   edat = (EphemerisData *)LALCalloc(1, sizeof(EphemerisData));
   (*edat).ephiles.earthEphemeris = uvar_earthEphemeris;
   (*edat).ephiles.sunEphemeris = uvar_sunEphemeris;

   LAL_CALL( LALLeapSecs(&status, &tmpLeap, &firstTimeStamp, &lsfas), &status);
   (*edat).leap = (INT2)tmpLeap;

   LAL_CALL( LALInitBarycenter( &status, edat), &status);

   /* polarisation angle */
   psi = uvar_psi;

   /* curly As */
   Aplus = (1.0 + (uvar_cosi*uvar_cosi))/2.0;
   Across = uvar_cosi;

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


   /* initialise output arrays */
   nParams = nSkyPatches * nfreqLoops;

   weights = XLALCreateREAL8Vector(nParams);
  
   stddev = LALCalloc(1, sizeof(REAL8));

   /* start frequency loop */

   for (freqCounter = 0; freqCounter < nfreqLoops; freqCounter++) {
/*printf("start of loop, nparams %i \n", nParams);*/

   MultiSFTVector *multiSFTs = NULL;


    /* read the sfts */
    /* first load them into a MultiSFTVector, then concatenate the various vectors into one*/
    /* read sft files making sure to add extra bins for running median */
    /* add wings for Doppler modulation and running median block size*/
    doppWings = (uvar_f0 + (freqCounter*deltaF) + uvar_fBand) * VTOT;    
    fMin = uvar_f0 + (freqCounter*deltaF) - doppWings - uvar_blocksRngMed * deltaF;
    fMax = uvar_f0 + (freqCounter*deltaF) + uvar_fBand + doppWings + uvar_blocksRngMed * deltaF;
 
/*printf("fmin, fmax %f %f\n", fMin, fMax);*/

    LAL_CALL( LALLoadMultiSFTs ( &status, &multiSFTs, catalog, fMin, fMax), &status);
    LAL_CALL (CombineAllSFTs (&status, &inputSFTs, multiSFTs, catalog->length), &status);


    /* find number of sfts */
    /* loop over ifos and calculate number of sfts */
    /* note that we can't use the catalog to determine the number of SFTs
       because SFTs might be segmented in frequency */
    numsft = inputSFTs->length; 

    /* normalize sfts */
    /* get information about all detectors including velocity and timestamps */
    /* note that this function returns the velocity at the 
     mid-time of the SFTs -- should not make any difference */
    psdVec = (PSDVector *) LALCalloc(1, sizeof(PSDVector)); 
    psdVec->length = numsft;

    psdVec->data = (REAL8FrequencySeries *)LALCalloc (psdVec->length, sizeof(REAL8FrequencySeries));


    ts = XLALCreateTimestampVector (1);
    tOffs = 0.5/deltaF;
    /*loop over all sfts and get the PSDs and detector states */
    for (j=0; j < (INT4)numsft; j++) {
  	psdVec->data[j].data = NULL;
	psdVec->data[j].data = (REAL8Sequence *)LALCalloc (1,sizeof(REAL8Sequence));
	psdVec->data[j].data->length = inputSFTs->data[j].data->length;
	psdVec->data[j].data->data = (REAL8 *)LALCalloc (inputSFTs->data[j].data->length, sizeof(REAL8));

	LAL_CALL( LALNormalizeSFT (&status, psdVec->data + j, inputSFTs->data + j, uvar_blocksRngMed), &status);
	
    }

    pairParams.lag = uvar_maxlag;  

    /* create sft pair indices */
    LAL_CALL ( CreateSFTPairsIndicesFrom2SFTvectors( &status, &sftPairIndexList, inputSFTs, &pairParams, uvar_detChoice), &status);

    /* initialise F_+, F_x vectors */

    Fplus = XLALCreateREAL8Vector(numsft);
    Fcross = XLALCreateREAL8Vector(numsft);
    frequencyShiftList = XLALCreateREAL8Vector(numsft);
    signalPhaseList = XLALCreateREAL8Vector(numsft);

    /* initialise frequency and phase vectors *
     * initialise Y, u, weight vectors */
    yalpha = (COMPLEX16Vector *) LALCalloc(1, sizeof(COMPLEX16));
    yalpha->length = sftPairIndexList->vectorLength;
    yalpha->data = LALCalloc(yalpha->length, sizeof(COMPLEX16));
    ualpha = (COMPLEX16Vector *) LALCalloc(1, sizeof(COMPLEX16));
    ualpha->length = sftPairIndexList->vectorLength;
    ualpha->data = LALCalloc(ualpha->length, sizeof(COMPLEX16));

    sigmasq = XLALCreateREAL8Vector(sftPairIndexList->vectorLength);


  /*  AMcoef = (AMCoeffs *)LALCalloc(1, sizeof(AMCoeffs));
    AMcoef->a = XLALCreateREAL4Vector(numsft);
    AMcoef->b = XLALCreateREAL4Vector(numsft);
*/
/*printf("starting loops over sky patches\n");*/

    /* loop over sky patches -- main calculations  */
    for (skyCounter = 0; skyCounter < nSkyPatches; skyCounter++) 
       { 

   	/* initialize Doppler parameters of the potential source */
	thisPoint.Alpha = skyAlpha[skyCounter]; 
	thisPoint.Delta = skyDelta[skyCounter]; 
	thisPoint.fkdot[0] = uvar_f0 + (freqCounter*deltaF);
	thisPoint.fkdot[1] = uvar_fdot; 
	thisPoint.fkdot[2] = 0.0;
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

    /* loop over each SFT to get frequency, then store in frequencyShiftList */
    for (j=0; j < (INT4)numsft; j++) {

	/*only have 1 element in detectorStateSeries and AMCoeffs because
	  the detector has to be updated for every SFT */
	DetectorStateSeries *detState = NULL;
    	AMcoef = (AMCoeffs *)LALCalloc(1, sizeof(AMCoeffs));
    	AMcoef->a = XLALCreateREAL4Vector(1);
    	AMcoef->b = XLALCreateREAL4Vector(1);



	det = XLALGetSiteInfo (inputSFTs->data[j].name); 
	ts->data[0] = inputSFTs->data[j].epoch;

	LAL_CALL ( LALGetDetectorStates ( &status, &detState, ts, det, edat, tOffs), &status);

        detState->detector = *det;

	LAL_CALL ( LALGetAMCoeffs ( &status, AMcoef, detState, skypos), &status); 

	thisVel.data = detState->data[0].vDetector;
 	thisPos.data = detState->data[0].rDetector;
	sft = &(inputSFTs->data[j]);


        LAL_CALL( GetSignalFrequencyInSFT( &status, freq1, sft, &thisPoint, &thisVel, &firstTimeStamp), &status);
        LAL_CALL( GetSignalPhaseInSFT( &status, phase1, sft, &thisPoint, &thisPos), &status);

	frequencyShiftList->data[j] = *freq1; 
	signalPhaseList->data[j] = *phase1;

	Fplus->data[j] = (AMcoef->a->data[0] * cos(2.0*psi))
		 	     + (AMcoef->b->data[0] * sin(2.0*psi));

	Fcross->data[j] = (AMcoef->b->data[0] * cos(2.0*psi))
			     - (AMcoef->a->data[0] * sin(2.0*psi));
	sft = NULL;
	
	XLALDestroyDetectorStateSeries(detState);
	LALFree(det);
 	XLALDestroyAMCoeffs ( AMcoef ); 

    } 
    
     
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
  	
	LAL_CALL( CorrelateSingleSFTPair( &status, &(yalpha->data[j]), sft1, sft2, psd1, psd2,  &(frequencyShiftList->data[index1]), &(frequencyShiftList->data[index2])), &status);
/*printf("freqs index%f %i\n", frequencyShiftList->data[index1], index1);*/


 	LAL_CALL( CalculateSigmaAlphaSq( &status, &sigmasq->data[j], frequencyShiftList->data[index1], frequencyShiftList->data[index2], psd1, psd2), &status);
/*printf("freqs end %f\n", frequencyShiftList->data[index1]);*/
 
	LAL_CALL( CalculateUalpha (&status, &ualpha->data[j], &Aplus, &Across, &signalPhaseList->data[index1], &signalPhaseList->data[index2], &Fplus->data[index1], &Fplus->data[index2], &Fcross->data[index1], &Fcross->data[index2], &sigmasq->data[j]), &status);
     	ualphacounter = ualphacounter + 1;

    }


    /* calculate rho (weights) from Yalpha and Ualpha */
    LAL_CALL( CalculateCrossCorrPower( &status, &(weights->data[counter]), yalpha, ualpha), &status);

    /* calculate standard deviation of rho (Eq 4.6) */
    LAL_CALL( NormaliseCrossCorrPower( &status, stddev, ualpha, sigmasq), &status); 

    /* select candidates  */
    /* print all interesting variables to file */

    fprintf(fp, "%1.5f\t %1.5f\t %1.5f\t%1.10f\t%1.10f\t%1.10f\n", thisPoint.Alpha, thisPoint.Delta, uvar_f0 + (freqCounter*deltaF), weights->data[counter], *stddev, weights->data[counter]/(*stddev));

/*printf("%1.5f\t%1.10f\t%1.10f\t%1.10f\n", uvar_f0 + (counter*deltaF), weights->data[counter], *stddev, weights->data[counter]/(*stddev));*/

    weights->data[counter] = weights->data[counter++]/(*stddev);


   } /* finish loop over skypatches */ 
printf("Frequency %f\n", uvar_f0 + (freqCounter*deltaF));

   /* free memory */

   LAL_CALL(LALDestroyMultiSFTVector(&status, &multiSFTs), &status);

/*printf("end of freq loop\n");*/

   XLALDestroyTimestampVector(ts);

   XLALDestroyCOMPLEX16Vector(yalpha);
   XLALDestroyCOMPLEX16Vector(ualpha);
   XLALDestroyREAL8Vector(frequencyShiftList);
   XLALDestroyREAL8Vector(signalPhaseList);
   XLALDestroyREAL8Vector(Fcross);
   XLALDestroyREAL8Vector(Fplus);
   LAL_CALL ( LALDestroyPSDVector  ( &status, &psdVec), &status);
   LAL_CALL (LALDestroySFTVector(&status, &inputSFTs), &status );
   LALFree(sftPairIndexList->data);
   LALFree(sftPairIndexList);
   XLALDestroyREAL8Vector(sigmasq);
 /*  XLALDestroyDetectorStateSeries ( detStates );*/
   } /*finish loop over frequencies */

  
   fclose (fp);


   LAL_CALL( LALDestroySFTCatalog( &status, &catalog ), &status);  



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

