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


INT4 lalDebugLevel;

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
#define SFTDIRECTORY "/home/badkri/fakesfts/"
#define FNAMEOUT "./OutHoughFStat"

int main( int argc, char *argv[]) {
  LALStatus status = blank_status;	/* initialize status */
  
  INT4 j,k; /* temp loop variables: k loops over stacks and j over SFTs in a stack*/

  /* detector, ephemeris and velocity/position vector */
  LALDetector detector;
  EphemerisData *edat = NULL;
  TimeVelPosVector *timeVelPos=NULL;   
  LIGOTimeGPSVector *midTs=NULL; /* time stamps of mid points of SFTs */

  /* sft related stuff */
  SFTVector *inputSFTs=NULL;  /* vector of SFTtypes and SFTtype is COMPLEX8FrequencySeries */
  INT4 mCohSft, nSFTs; /* number of SFTs in each stack and total number of SFTs */
  INT4 Nstacks; /* number of stacks -- not necessarily same as uvar_Nstacks! */
  INT4 sftlength; /* number of bins in each sft */
  REAL8 deltaF, timeBase; /* frequency resolution */
  INT8 sftFminBin; /* first sft bin index */
  INT4 medianBias; /* running median bias */

  /* LALdemod related stuff */
  LALFstat *Fstat;	
  /* FFT is a structure containing *SFTtype 
     --- for some reason COMPLEX8FrequencySeriesVector is not used */
  FFT **SFTData=NULL;
  DemodPar *DemodParams;
  CSParams *csParams  = NULL;  /* ComputeSky parameters */
  AMCoeffsParams *amParams;
  AMCoeffs amc;   
  EarthState earth;
  EmissionTime emit;
  BarycenterInput baryinput; /* Stores detector location and other barycentering data */
  
  /* hough variables */
  INT4  nfSizeCylinder=NFSIZE;

  /* variables for logging */
  CHAR *fnamelog=NULL;
  FILE *fpLog = NULL;
  CHAR *logstr=NULL; 

  /* user variables */
  BOOLEAN uvar_help; /* true if -h option is given */
  REAL8 uvar_alpha, uvar_delta;  /* sky-location angles */
  REAL8 uvar_fdot; /* first spindown value */
  REAL8 uvar_fStart, uvar_fBand;
  INT4 uvar_ifo, uvar_blocksRngMed, uvar_Nstacks, uvar_Dterms;
  CHAR *uvar_earthEphemeris=NULL;
  CHAR *uvar_sunEphemeris=NULL;
  CHAR *uvar_sftDir=NULL;
  CHAR *uvar_fnameout=NULL;

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /************/
  /* set defaults, read user variables, log user variables and log cvs tags */
  /* LALDebugLevel must be called before anything else */
  lalDebugLevel = 0;
  LAL_CALL( LALGetDebugLevel( &status, argc, argv, 'd'), &status);

  /* now set the defaults */
  uvar_help = FALSE;
  uvar_Nstacks = NSTACKS;
  uvar_Dterms = DTERMS;
  uvar_alpha = ALPHA;
  uvar_delta = DELTA;
  uvar_fdot = FDOT;
  uvar_fStart = FSTART;
  uvar_fBand = FBAND;
  uvar_ifo = IFO;
  uvar_blocksRngMed = BLOCKSRNGMED;
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
  LAL_CALL( LALRegisterINTUserVar(    &status, "Nstacks",         'N', UVAR_OPTIONAL, "Number of stacks",              &uvar_Nstacks ),        &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "Dterms",          'N', UVAR_OPTIONAL, "For Dirichlet Kernel approx.",  &uvar_Dterms ),         &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed",    'w', UVAR_OPTIONAL, "RngMed block size",             &uvar_blocksRngMed),    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "earthEphemeris",  'E', UVAR_OPTIONAL, "Earth Ephemeris file",          &uvar_earthEphemeris),  &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sunEphemeris",    'S', UVAR_OPTIONAL, "Sun Ephemeris file",            &uvar_sunEphemeris),    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftDir",          'D', UVAR_OPTIONAL, "SFT Directory",                 &uvar_sftDir),          &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "fnameout",        'o', UVAR_OPTIONAL, "Output basefileneme",           &uvar_fnameout),        &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "alpha",           'r', UVAR_OPTIONAL, "Right ascension",               &uvar_alpha),           &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "delta",           'l', UVAR_OPTIONAL, "Declination",                   &uvar_delta),           &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fStart",          'f', UVAR_OPTIONAL, "Start search frequency",        &uvar_fStart),          &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fdot",            'd', UVAR_OPTIONAL, "Spindown parameter",            &uvar_fdot),            &status);

  /* read all command line variables */
  LAL_CALL( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 

  /* some sanity checks on user vars */
  if ( uvar_Nstacks < 1) {
    fprintf(stderr, "Invalid number of stacks\n");
    exit(1);
  }
  if ( uvar_blocksRngMed < 1 ) {
    fprintf(stderr, "Invalid Running Median block size\n");
    exit(1);
  }

  /* set detector */
  if (uvar_ifo ==1) detector=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if (uvar_ifo ==2) detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
  if (uvar_ifo ==3) detector=lalCachedDetectors[LALDetectorIndexLHODIFF];

  LAL_CALL ( LALRngMedBias( &status, &medianBias, uvar_blocksRngMed ), &status);

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
  }
  /********* end of user var reading etc. ******/


  /** read sfts and set up timestamp vectors **/
  {
    CHAR *tempDir;
    tempDir = (CHAR *)LALMalloc(512*sizeof(CHAR));
    strcpy(tempDir, uvar_sftDir);
    strcat(tempDir, "/*SFT*.*"); 
    LAL_CALL( LALReadSFTfiles ( &status, &inputSFTs, uvar_fStart, uvar_fStart + uvar_fBand, nfSizeCylinder + uvar_blocksRngMed , tempDir), &status);
    LALFree(tempDir);
  }

  /* set number of stacks */
  nSFTs = inputSFTs->length;
  Nstacks = uvar_Nstacks;
  mCohSft = nSFTs/Nstacks;
  
  /* if too many stacks were chosen by user */  
  if ( mCohSft == 0 ) {
    Nstacks = nSFTs;
    mCohSft = 1;
    if ( lalDebugLevel > 0)
      fprintf(stderr, "Warning: requested number of stacks exceeds number of SFTs...setting Nstacks = nSFTs\n"); 
  }
  /* notify user about any discarded sfts */
  if ( (lalDebugLevel > 0) && (nSFTs - mCohSft*Nstacks) )
    fprintf(stderr, "Warning: last %d SFTs are discarded due to rounding off\n", nSFTs - mCohSft*Nstacks);
  
  /* set other sft parameters */
  sftlength = inputSFTs->data->data->length;
  deltaF = inputSFTs->data->deltaF;
  timeBase = 1.0/deltaF;
  sftFminBin = floor( timeBase * inputSFTs->data->f0 + 0.5);
  
  /* create and fill timestamp vectors */
  timeVelPos = (TimeVelPosVector *)LALMalloc( Nstacks * sizeof(TimeVelPosVector));
  for (k=0; k<Nstacks; k++) {
    timeVelPos[k].length = mCohSft;
    timeVelPos[k].ts = (LIGOTimeGPS *)LALMalloc( mCohSft * sizeof(LIGOTimeGPS));
    timeVelPos[k].velx = (REAL8 *)LALMalloc( mCohSft * sizeof(REAL8));
    timeVelPos[k].vely = (REAL8 *)LALMalloc( mCohSft * sizeof(REAL8));
    timeVelPos[k].velz = (REAL8 *)LALMalloc( mCohSft * sizeof(REAL8));
    timeVelPos[k].posx = (REAL8 *)LALMalloc( mCohSft * sizeof(REAL8));
    timeVelPos[k].posy = (REAL8 *)LALMalloc( mCohSft * sizeof(REAL8));
    timeVelPos[k].posz = (REAL8 *)LALMalloc( mCohSft * sizeof(REAL8));

    /* assign timestamp values from SFTs */  
    for (j = 0; j<mCohSft; j++) {
      timeVelPos[k].ts[j].gpsSeconds = inputSFTs->data[j+mCohSft*k].epoch.gpsSeconds;
      timeVelPos[k].ts[j].gpsNanoSeconds = inputSFTs->data[j+mCohSft*k].epoch.gpsNanoSeconds;	
    }
  }  
  /* calculate mid time stamps of sfts */
  midTs = (LIGOTimeGPSVector *)LALMalloc( mCohSft * sizeof(LIGOTimeGPSVector));
  for (k=0; k<Nstacks; k++) {
    midTs[k].length = mCohSft;
    midTs[k].data  = NULL;
    midTs[k].data = (LIGOTimeGPS *)LALMalloc( mCohSft * sizeof(LIGOTimeGPS));
    for (j=0; j<mCohSft; j++) {
      REAL8 tempTime;

      LAL_CALL( LALGPStoFloat ( &status, &tempTime, &(timeVelPos[k].ts[j])), &status);
      tempTime += 0.5 * timeBase;
      LAL_CALL( LALFloatToGPS ( &status, &(midTs[k].data[j]), &tempTime), &status);
    }
  }

  /* calculate velocity and position */
  /* setting of ephemeris info */ 
  edat = (EphemerisData *)LALMalloc(sizeof(EphemerisData));
  (*edat).ephiles.earthEphemeris = uvar_earthEphemeris;
  (*edat).ephiles.sunEphemeris = uvar_sunEphemeris;
  {  
    VelocityPar   velPar;
    REAL8     vel[3], pos[3]; 
    
    LALLeapSecFormatAndAcc lsfas = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
    INT4 tmpLeap; /* need this because Date pkg defines leap seconds as
		     INT4, while EphemerisData defines it to be INT2. This won't
                   cause problems before, oh, I don't know, the Earth has been 
                   destroyed in nuclear holocaust. -- dwchin 2004-02-29 */
    
    velPar.detector = detector;
    velPar.tBase = timeBase;
    velPar.vTol = 0.01; /* accuracy for position calculation */
    velPar.edat = NULL;
    
    /* Leap seconds for the start time of the run */   
    LAL_CALL( LALLeapSecs(&status, &tmpLeap, &(timeVelPos[0].ts[0]), &lsfas), &status);
    (*edat).leap = (INT2)tmpLeap;
    
    /* read in ephemeris data */
    LAL_CALL( LALInitBarycenter( &status, edat), &status);
    velPar.edat = edat;

    /* calculate detector velocity and position */
    for (k=0; k<Nstacks; k++){
      for (j=0; j<mCohSft; ++j){
	velPar.startTime.gpsSeconds     = timeVelPos[k].ts[j].gpsSeconds;
	velPar.startTime.gpsNanoSeconds = timeVelPos[k].ts[j].gpsNanoSeconds;
	
	LAL_CALL( LALAvgDetectorVel ( &status, vel, &velPar), &status );
	LAL_CALL( LALAvgDetectorPos ( &status, pos, &velPar), &status );
	timeVelPos[k].velx[j] = vel[0];
	timeVelPos[k].vely[j] = vel[1];
	timeVelPos[k].velz[j] = vel[2];   
	timeVelPos[k].posx[j] = pos[0];
	timeVelPos[k].posy[j] = pos[1];  
	timeVelPos[k].posz[j] = pos[2];    	
      }  
    }
  }

  /************ setup parameters for calling  LALDemod**********/

  /* Detector location */
  baryinput.site.location[0] = detector.location[0]/LAL_C_SI;
  baryinput.site.location[1] = detector.location[1]/LAL_C_SI;
  baryinput.site.location[2] = detector.location[2]/LAL_C_SI;
  baryinput.alpha = uvar_alpha;
  baryinput.delta = uvar_delta;
  baryinput.dInv = 0.e0;

  /* amParams stuff */
  amParams = (AMCoeffsParams *)LALMalloc(sizeof(AMCoeffsParams));
  amParams->das = (LALDetAndSource *)LALMalloc(sizeof(LALDetAndSource));
  amParams->das->pSource = (LALSource *)LALMalloc(sizeof(LALSource));
  amParams->baryinput = &baryinput;
  amParams->earth = &earth; 
  amParams->edat = edat;
  amParams->das->pDetector = &detector; 
  amParams->das->pSource->equatorialCoords.latitude = uvar_delta;
  amParams->das->pSource->equatorialCoords.longitude = uvar_alpha;
  amParams->das->pSource->orientation = 0.0;
  amParams->das->pSource->equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  amParams->polAngle = amParams->das->pSource->orientation ;
  amParams->leapAcc = LALLEAPSEC_STRICT;
  /* memory for output of LALComputeAM */
  amc.a = NULL;
  amc.b = NULL;
  LAL_CALL (LALSCreateVector( &status, &(amc.a), mCohSft), &status);
  LAL_CALL (LALSCreateVector( &status, &(amc.b), mCohSft), &status);

  /* ComputeSky stuff*/
  csParams = (CSParams *)LALMalloc(sizeof(CSParams));
  csParams->skyPos = (REAL8 *)LALMalloc(2*sizeof(REAL8));
  csParams->mObsSFT = mCohSft;   
  csParams->tSFT = timeBase;
  csParams->edat = edat;
  csParams->baryinput = &baryinput;
  csParams->spinDwnOrder = 1;
  csParams->skyPos[0] = uvar_alpha;
  csParams->skyPos[1] = uvar_delta;
  csParams->earth = &earth;
  csParams->emit = &emit;

  /* DemodParams stuff */
  /* Allocate DemodParams structure */
  DemodParams = (DemodPar *)LALCalloc(1, sizeof(DemodPar));
  /* space for sky constants */
  /* Based on maximum index for array of as and bs sky constants as from ComputeSky.c */
  DemodParams->skyConst = (REAL8 *)LALMalloc(4 * mCohSft * sizeof(REAL8));
  /* space for spin down params */
  DemodParams->spinDwnOrder = 1;
  DemodParams->spinDwn = (REAL8 *)LALMalloc(sizeof(REAL8));
  DemodParams->spinDwn[0] = uvar_fdot; /* needs correction for reference time */
  DemodParams->SFTno = mCohSft;
  DemodParams->f0 = uvar_fStart;
  DemodParams->imax = (INT4)(uvar_fBand/deltaF + 1e-6) + 1;
  DemodParams->Dterms = uvar_Dterms;
  DemodParams->df   = deltaF;
  DemodParams->ifmin = (INT4) floor( (1.0 - VTOT)* uvar_fStart * timeBase) - uvar_Dterms;;
  DemodParams->returnFaFb = FALSE; /* don't require Fa and Fb */

  /* allocate memory for Fstat structure for each stack -- output from LALDemod*/
  Fstat = (LALFstat *)LALMalloc( Nstacks * sizeof(LALFstat));
  /* allocate memory for SFTData */
  SFTData = (FFT **)LALMalloc( mCohSft * sizeof(FFT *));  
  for (j=0; j<mCohSft; j++) {
    SFTData[j]=NULL;
    SFTData[j] = (FFT *)LALMalloc(sizeof(FFT));
  }

  /* calculate Fstat for each stack */
  for (k=0; k<Nstacks; k++) {
    
    Fstat[k].F = (REAL8 *)LALMalloc( DemodParams->imax * sizeof(REAL8));

    /* compute a(t) and b(t) and A,B,C,D */
    LAL_CALL ( LALComputeAM( &status, &amc, midTs[k].data, amParams), &status); 
    DemodParams->amcoe = &amc;

    /* compute the "sky-constants" */
    csParams->tGPS = timeVelPos[k].ts;
    LAL_CALL ( LALComputeSky( &status, DemodParams->skyConst, 0, csParams), &status);  

    for (j=0; j<mCohSft; j++) {
      SFTData[j]->fft = &(inputSFTs->data[k*mCohSft + j]);
    }

    LAL_CALL ( LALDemod(&status, Fstat + k, SFTData, DemodParams), &status);    
  }

  /* free timestamp and Vel/Pos vectors */
  for (k=0; k<Nstacks; k++) {
      LALFree(timeVelPos[k].ts);
      LALFree(timeVelPos[k].posx);
      LALFree(timeVelPos[k].posy);
      LALFree(timeVelPos[k].posz);
      LALFree(timeVelPos[k].velx);
      LALFree(timeVelPos[k].vely);
      LALFree(timeVelPos[k].velz);
      LALFree(midTs[k].data);
    }  
  LALFree(timeVelPos);
  LALFree(midTs);

  /* free ephemeris */
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);

  /* free amParams */
  LALFree(amParams->das->pSource);
  LALFree(amParams->das);
  LALFree(amParams);
  LAL_CALL ( LALSDestroyVector( &status, &(amc.a)), &status);
  LAL_CALL ( LALSDestroyVector( &status, &(amc.b)), &status);

  /* free ComputeSky Params */
  LALFree(csParams->skyPos);
  LALFree(csParams);

  /* free Fstat structure */
  for (k=0; k<Nstacks; k++)
    LALFree(Fstat[k].F);
  LALFree(Fstat);

  for (j=0; j<mCohSft; j++)
    LALFree(SFTData[j]);
  LALFree(SFTData);


  /* free Demod params */
  LALFree(DemodParams->skyConst);
  LALFree(DemodParams->spinDwn);
  LALFree(DemodParams);

  LAL_CALL (LALDestroySFTVector(&status, &inputSFTs),&status );

  LAL_CALL (LALDestroyUserVars(&status), &status);  

  LALCheckMemoryLeaks();

  return 0;
}
