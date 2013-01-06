/*
*  Copyright (C) 2007 Badri Krishnan, Reinhard Prix, Alicia Sintes Olives
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
 * \author Badri Krishnan and Alicia Sintes
 * \brief  Checking improvements in Hough if amplitude modulation is considered.
 */


#include "../src/MCInjectHoughS2.h"
#include <lal/LALComputeAM.h>
#include <lal/ComputeSky.h>
#include <gsl/gsl_cdf.h>

extern INT4 lalDebugLevel;

/* defaults */
#define EARTHEPHEMERIS "../src/earth00-04.dat"
#define SUNEPHEMERIS "../src/sun00-04.dat"
#define MAXFILES 3000 /* maximum number of files to read in a directory */
#define MAXFILENAMELENGTH 256 /* maximum # of characters  of a SFT filename */
#define IFO 2         /*  detector, 1:GEO, 2:LLO, 3:LHO */
#define THRESHOLD 1.6 /* thresold for peak selection, with respect to the */
#define NFSIZE  21 /* n-freq. span of the cylinder, to account for spin-down */
#define BLOCKSRNGMED 101 /* Running median window size */

/* default injected pulsar parameters */
#define F0 255.0          /*  frequency to build the LUT and start search */
#define FDOT 0.0 /* default spindown parameter */
#define ALPHA 1.57
#define DELTA  0.0
#define COSIOTA 0.5
#define PHI0 0.0
#define PSI 0.0
#define H0 (1.0e-23)

/* default file and directory names */
#define SFTDIRECTORY "/local_data/badkri/fakesfts/"
#define FILEOUT "./ValidateAMOut"   
#define TRUE (1==1)
#define FALSE (1==0)

int main( int argc, char *argv[]){

  static LALStatus            status;  
  static LALDetector          detector;
  static LIGOTimeGPSVector    timeV;
  static REAL8Cart3CoorVector velV;
  static REAL8Vector          timeDiffV;
  static REAL8Vector          foft;
  static PulsarSignalParams  params;
  static SFTParams sftParams;

  static UCHARPeakGram     pg1; 
  /*   static UCHARPeakGram     pg2; */
  static COMPLEX8SFTData1  sft1;
  static REAL8PeriodoPSD   periPSD;

  REAL4TimeSeries   *signalTseries = NULL;
  SFTVector    *inputSFTs  = NULL;  
  SFTVector    *outputSFTs = NULL;
  /* data about injected signal */
  static PulsarData           pulsarInject;

  /* the template */
  static HoughTemplate  pulsarTemplate, pulsarTemplate1;

  /* FILE  *fpOUT = NULL; */ /* output file pointer */
  FILE  *fpLog = NULL; /* log file pointer */
  CHAR  *logstr=NULL; /* log string containing user input variables */
  CHAR  *fnamelog=NULL; /* name of log file */
  INT4  nfSizeCylinder;
  
  EphemerisData   *edat = NULL;

  INT4   mObsCoh;
  REAL8  numberCount2, numberCount1;
  REAL8  nth1, nth2, erfcInv;
  /*   REAl8  meanAM, sigmaAM, varAM, mean, sigma; */
  REAL8  mean1, mean2, var1, var2, sumsq, peakprob;
  REAL8  sftBand;  
  REAL8  timeBase, deltaF, normalizeThr, threshold;
  /*   REAL8  thresholdAM; */
  UINT4  sftlength; 
  INT4   sftFminBin;
  REAL8  fHeterodyne;
  REAL8  tSamplingRate;

  /* grid spacings */
  REAL8 deltaTheta;
  INT4 mmP, mmT; /* for loop over mismatched templates */

  /* user input variables */
  BOOLEAN uvar_help;
  INT4 uvar_ifo, uvar_blocksRngMed;
  REAL8 uvar_houghFalseAlarm;
  REAL8 uvar_peakThreshold;
  REAL8 uvar_alpha, uvar_delta, uvar_h0, uvar_f0;
  REAL8 uvar_psi, uvar_phi0, uvar_fdot, uvar_cosiota;
  REAL8 uvar_mismatchW;
  CHAR *uvar_earthEphemeris=NULL;
  CHAR *uvar_sunEphemeris=NULL;
  CHAR *uvar_sftDir=NULL;
  CHAR *uvar_fnameout=NULL;

  /* vector of weights */
  REAL8Vector *weight;

  /*  set up the default parameters  */
  lalDebugLevel = 0;

  nfSizeCylinder = NFSIZE;
  
  /* LALDebugLevel must be called before anything else */
  LAL_CALL( LALGetDebugLevel( &status, argc, argv, 'd'), &status);

  /* *********************************************************************** */
  /* set other user input variables */
  uvar_help = FALSE;
  uvar_peakThreshold = THRESHOLD;
  uvar_ifo = IFO;
  uvar_blocksRngMed = BLOCKSRNGMED;
  uvar_mismatchW = 0.0;
  uvar_houghFalseAlarm = 1.0e-10;

  /* set default pulsar parameters */
  uvar_h0 = H0;
  uvar_alpha = ALPHA;
  uvar_delta = DELTA;
  uvar_f0 =  F0;
  uvar_fdot = FDOT;
  uvar_psi = PSI;
  uvar_cosiota = COSIOTA;
  uvar_phi0 = PHI0;
  
  /* now set the default filenames */
  uvar_earthEphemeris = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_earthEphemeris,EARTHEPHEMERIS);

  uvar_sunEphemeris = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_sunEphemeris,SUNEPHEMERIS);

  uvar_sftDir = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_sftDir,SFTDIRECTORY);

  uvar_fnameout = (CHAR *)LALMalloc(512*sizeof(CHAR));
  strcpy(uvar_fnameout, FILEOUT);

  /* *********************************************************************** */
  /* register user input variables */
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "help",            'h', UVAR_HELP,     "Print this message",             &uvar_help),            &status);  
  LAL_CALL( LALRegisterINTUserVar(    &status, "ifo",             'i', UVAR_OPTIONAL, "Detector GEO(1) LLO(2) LHO(3)",  &uvar_ifo ),            &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "blocksRngMed",    'w', UVAR_OPTIONAL, "RngMed block size",              &uvar_blocksRngMed),    &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "peakThreshold",   't', UVAR_OPTIONAL, "Peak selection threshold",       &uvar_peakThreshold),   &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "houghFalseAlarm",  0,  UVAR_OPTIONAL, "Overall Hough False alarm",      &uvar_houghFalseAlarm), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "earthEphemeris",  'E', UVAR_OPTIONAL, "Earth Ephemeris file",           &uvar_earthEphemeris),  &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sunEphemeris",    'S', UVAR_OPTIONAL, "Sun Ephemeris file",             &uvar_sunEphemeris),    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftDir",          'D', UVAR_OPTIONAL, "SFT Directory",                  &uvar_sftDir),          &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "fnameout",        'o', UVAR_OPTIONAL, "Output file prefix",             &uvar_fnameout),        &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "alpha",           'r', UVAR_OPTIONAL, "Right ascension",                &uvar_alpha),           &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "delta",           'l', UVAR_OPTIONAL, "Declination",                    &uvar_delta),           &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "h0",              'm', UVAR_OPTIONAL, "h0 to inject",                   &uvar_h0),              &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "f0",              'f', UVAR_OPTIONAL, "Signal frequency",               &uvar_f0),              &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "psi",             'p', UVAR_OPTIONAL, "Polarization angle",             &uvar_psi),             &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "phi0",            'P', UVAR_OPTIONAL, "Initial phase",                  &uvar_phi0),            &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "cosiota",         'c', UVAR_OPTIONAL, "Cosine of iota",                 &uvar_cosiota),         &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fdot",            'd', UVAR_OPTIONAL, "Spindown parameter",             &uvar_fdot),            &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "mismatch",        'M', UVAR_OPTIONAL, "Mismatch in weight calculation", &uvar_mismatchW),       &status);

  /* read all command line variables */
  LAL_CALL( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 
  
  /* *********************************************************************** */
  /* write the log file */
  fnamelog = (CHAR *)LALMalloc( 512*sizeof(CHAR));
  strcpy(fnamelog, uvar_fnameout);
  strcat(fnamelog, "_log");
  /* open the log file for writing */
  if ((fpLog = fopen(fnamelog, "w")) == NULL) {
    fprintf(stderr, "Unable to open file %s for writing\n", fnamelog);
    LALFree(fnamelog);
    exit(1);
  }

  /* get the log string */
  LAL_CALL( LALUserVarGetLog(&status, &logstr, UVAR_LOGFMT_CFGFILE), &status);  

  fprintf( fpLog, "## Log file for HoughValidateAM \n\n");
  fprintf( fpLog, "# User Input:\n");
  fprintf( fpLog, "#-------------------------------------------\n");
  fprintf( fpLog, "%s", logstr);
  LALFree(logstr);

  /* append an ident-string defining the exact CVS-version of the code used */
  {
    CHAR command[1024] = "";
    fprintf (fpLog, "\n\n# CVS-versions of executable:\n");
    fprintf (fpLog, "# -----------------------------------------\n");
    fclose (fpLog);
    
    sprintf (command, "ident %s | sort -u >> %s", argv[0], fnamelog);
    /* we don't check this. If it fails, we assume that */
    /* one of the system-commands was not available, and */
    /* therefore the CVS-versions will not be logged */
    if ( system(command) ) fprintf (stderr, "\nsystem('%s') returned non-zero status!\n\n", command );

    LALFree(fnamelog); 
  }

  /* *********************************************************************** */
  /* set peak selection threshold */
  LAL_CALL( LALRngMedBias( &status, &normalizeThr, uvar_blocksRngMed ), &status ); 
  threshold = uvar_peakThreshold/normalizeThr; 

  /* *********************************************************************** */
  /* set detector */
  if (uvar_ifo ==1) detector=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if (uvar_ifo ==2) detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
  if (uvar_ifo ==3) detector=lalCachedDetectors[LALDetectorIndexLHODIFF];


  /* *********************************************************************** */
  /* copy user input values */
  pulsarInject.f0 = uvar_f0;
  pulsarInject.latitude = uvar_delta;
  pulsarInject.longitude = uvar_alpha;
  pulsarInject.aPlus = 0.5 * uvar_h0 * ( 1.0 + uvar_cosiota * uvar_cosiota );
  pulsarInject.aCross = uvar_h0 * uvar_cosiota;
  pulsarInject.psi = uvar_psi;
  pulsarInject.phi0 = uvar_phi0;
  pulsarInject.spindown.length = 1;
  pulsarInject.spindown.data = NULL;
  pulsarInject.spindown.data = (REAL8 *)LALMalloc(sizeof(REAL8));
  pulsarInject.spindown.data[0] = uvar_fdot;

  /* *********************************************************************** */
  /* copy these values also to the pulsar template */
  /* template is complately matched at this point */
  pulsarTemplate.f0 = uvar_f0;
  pulsarTemplate.latitude = uvar_delta;
  pulsarTemplate.longitude = uvar_alpha;
  pulsarTemplate.spindown.length = 1;
  pulsarTemplate.spindown.data = NULL;
  pulsarTemplate.spindown.data = (REAL8 *)LALMalloc(sizeof(REAL8));
  pulsarTemplate.spindown.data[0] = uvar_fdot;

  /* *********************************************************************** */
  /* allocate memory for mismatched spindown template */
  pulsarTemplate1.spindown.length = 1;
  pulsarTemplate1.spindown.data = NULL;
  pulsarTemplate1.spindown.data = (REAL8 *)LALMalloc(sizeof(REAL8));

  /* *********************************************************************** */
  /* read sfts */
  /* note that sft must be long enough -- at least 1Hz */
  {
    CHAR *tempDir;
    tempDir = (CHAR *)LALMalloc(512*sizeof(CHAR));
    strcpy(tempDir, uvar_sftDir);
    strcat(tempDir, "/*SFT*.*"); 
    sftBand = 0.5; 
    LAL_CALL( LALReadSFTfiles ( &status, &inputSFTs, uvar_f0 - sftBand, uvar_f0 + sftBand, nfSizeCylinder + uvar_blocksRngMed , tempDir), &status);
    LALFree(tempDir);
  }


  /* *********************************************************************** */
  /* get sft parameters */
  mObsCoh = inputSFTs->length;
  sftlength = inputSFTs->data->data->length;
  deltaF = inputSFTs->data->deltaF;
  timeBase = 1.0/deltaF;
  sftFminBin = floor( timeBase * inputSFTs->data->f0 + 0.5);

  /* signal generation parameters */
  fHeterodyne = sftFminBin*deltaF;
  tSamplingRate = 2.0*deltaF*(sftlength -1.);

  /* *********************************************************************** */
  /* create timestamp vector */
  timeV.length = mObsCoh;
  timeV.data = NULL;  
  timeV.data = (LIGOTimeGPS *)LALMalloc(mObsCoh*sizeof(LIGOTimeGPS));

  /* *********************************************************************** */
  /* read timestamps */
  { 
    INT4    i; 
    SFTtype  *sft= NULL; 
    
    sft = inputSFTs->data;
    for (i=0; i < mObsCoh; i++){
      timeV.data[i].gpsSeconds = sft->epoch.gpsSeconds;
      timeV.data[i].gpsNanoSeconds = sft->epoch.gpsNanoSeconds;
      ++sft;
    }    
  }

  /* *********************************************************************** */
  /* compute the time difference relative to startTime for all SFT */
  timeDiffV.length = mObsCoh;
  timeDiffV.data = NULL; 
  timeDiffV.data = (REAL8 *)LALMalloc(mObsCoh*sizeof(REAL8));

  {   
    REAL8   t0, ts, tn, midTimeBase;
    INT4   j; 

    midTimeBase=0.5*timeBase;
    ts = timeV.data[0].gpsSeconds;
    tn = timeV.data[0].gpsNanoSeconds * 1.00E-9;
    t0=ts+tn;
    timeDiffV.data[0] = midTimeBase;

    for (j=1; j< mObsCoh; ++j){
      ts = timeV.data[j].gpsSeconds;
      tn = timeV.data[j].gpsNanoSeconds * 1.00E-9;  
      timeDiffV.data[j] = ts+tn -t0+midTimeBase; 
    }  
  }

  /* *********************************************************************** */
  /* setting of ephemeris info */ 
  edat = (EphemerisData *)LALMalloc(sizeof(EphemerisData));
  (*edat).ephiles.earthEphemeris = uvar_earthEphemeris;
  (*edat).ephiles.sunEphemeris = uvar_sunEphemeris;
  
  /* *********************************************************************** */
  /* compute detector velocity for those time stamps  */ 
  velV.length = mObsCoh; 
  velV.data = NULL;
  velV.data = (REAL8Cart3Coor *)LALMalloc(mObsCoh*sizeof(REAL8Cart3Coor));
  
  {  
    VelocityPar   velPar;
    REAL8     vel[3]; 
    UINT4     j; 

    velPar.detector = detector;
    velPar.tBase = timeBase;
    velPar.vTol = 0.0; /* irrelevant */
    velPar.edat = NULL;

    /* read in ephemeris data */
    LAL_CALL( LALInitBarycenter( &status, edat), &status);
    velPar.edat = edat;

    /* calculate detector velocity */    
    for(j=0; j< velV.length; ++j){
      velPar.startTime.gpsSeconds     = timeV.data[j].gpsSeconds;
      velPar.startTime.gpsNanoSeconds = timeV.data[j].gpsNanoSeconds;
      
      LAL_CALL( LALAvgDetectorVel ( &status, vel, &velPar), &status );
      velV.data[j].x= vel[0];
      velV.data[j].y= vel[1];
      velV.data[j].z= vel[2];   
    }  
  }

  /* calculate weights based on amplitude modulation functions */

  /* allocate memory for weight */
  weight = NULL;
  LAL_CALL( LALDCreateVector( &status, &weight, mObsCoh), &status);

  /* calculate amplitude modulation weights */
  LAL_CALL( LALHOUGHInitializeWeights( &status, weight), &status);
  LAL_CALL( LALHOUGHComputeAMWeights( &status, weight, &timeV, &detector, edat, 
				      uvar_alpha + uvar_mismatchW, 
				      uvar_delta + uvar_mismatchW), &status);
  LAL_CALL( LALHOUGHNormalizeWeights( &status, weight), &status);
  

  /* calculate mean and variances */
  peakprob = exp(-uvar_peakThreshold); /* probability of peak selection */
  erfcInv = gsl_cdf_ugaussian_Qinv (uvar_houghFalseAlarm)/sqrt(2); /* erfcinv(2*uvar_houghFalseAlarm */

  /* usual number count threshold */
  /* 1 is for the usual hough transform 
     2 is for weighted hough transform */ 
  mean1 = mObsCoh * peakprob;
  var1 = mObsCoh * peakprob * (1 - peakprob);
  nth1 = mean1 + sqrt(2.0 * var1) * erfcInv;

  /* threshold for weighted hough number count */
  mean2 = mean1; /* this is due to normalization of weights */
  /* find sum of weights squared */
  sumsq = 0.0;
  INT4 j;
  for (j=0; j<mObsCoh; j++) {
    REAL8 tempW;
    tempW = weight->data[j];
    sumsq += tempW * tempW;
  }
  var2 = sumsq * peakprob * (1 - peakprob);
  nth2 = mean2 + sqrt( 2.0 * var2) * erfcInv;

  /*   { /\* check sum weight = 1 *\/ */
  /*     REAL8 tempSum=0.0; */
  /*     for (j=0; j<mObsCoh; j++) { */
  /*       tempSum += weight->data[j]; */
  /*     } */
  /*     fprintf(stdout, "%f\n", tempSum); */
  /*   } */
  /* *********************************************************************** */
  /* computing varAM */
  /*   { */
  /*     INT4 j; */
  /*     REAL8  x; */
  /*     /\* REAL8  y; *\/ */
  
  /*     varAM = 0.0; */
  
  /*     /\* y=0.0; *\/ */
  /*     for(j=0; j<mObsCoh; j++){ */
  /*       a = aVec->data[j]; */
  /*       b = bVec->data[j]; */
  /*       x = 2.0*(a*a + b*b)/(A+B) -1.0; */
  /*       varAM += x*x; */
  /*       /\* y+=x; *\/ */
  /*     } */
  /*     varAM = varAM/mObsCoh; */
  /*     /\* fprintf(stdout, "zero ?= %g\n", y); *\/ */
  
  /*   } */
  
  
  
  /* *********************************************************************** */
  /* set grid spacings */
  {
    deltaTheta = 1.0 / ( VTOT * uvar_f0 * timeBase );
    /* currently unused: REAL8 deltaFdot = deltaF / timeBase; */
  }

  /* *********************************************************************** */
  /* allocate memory for f(t) pattern */
  foft.length = mObsCoh;
  foft.data = NULL;
  foft.data = (REAL8 *)LALMalloc(mObsCoh*sizeof(REAL8));

  /* allocate memory for Hough peripsd structure */
  periPSD.periodogram.length = sftlength;
  periPSD.periodogram.data = NULL;
  periPSD.periodogram.data = (REAL8 *)LALMalloc(sftlength* sizeof(REAL8));
  periPSD.psd.length = sftlength;
  periPSD.psd.data = NULL;
  periPSD.psd.data = (REAL8 *)LALMalloc(sftlength* sizeof(REAL8));

  /* allocate memory for peakgrams */
  pg1.length = sftlength;
  pg1.data = NULL;
  pg1.data = (UCHAR *)LALMalloc(sftlength* sizeof(UCHAR));

  /*   pg2.length = sftlength; */
  /*   pg2.data = NULL; */
  /*   pg2.data = (UCHAR *)LALMalloc(sftlength* sizeof(UCHAR)); */

  /* *********************************************************************** */
  /* generate signal and add to input sfts */  
  /* parameters for output sfts */
  sftParams.Tsft = timeBase;
  sftParams.timestamps = &(timeV);
  sftParams.noiseSFTs = inputSFTs;  

  /* signal generation parameters */
  params.orbit = NULL;
  /* params.transferFunction = NULL; */
  params.site = &(detector);
  params.ephemerides = edat;
  params.startTimeGPS.gpsSeconds = timeV.data[0].gpsSeconds;   /* start time of output time series */
  params.startTimeGPS.gpsNanoSeconds = timeV.data[0].gpsNanoSeconds;   /* start time of output time series */
  params.duration = timeDiffV.data[mObsCoh-1] + 0.5 * timeBase; /* length of time series in seconds */
  params.samplingRate = tSamplingRate;
  params.fHeterodyne = fHeterodyne;
  /* reference time for frequency and spindown is first timestamp */
  params.pulsar.refTime.gpsSeconds = timeV.data[0].gpsSeconds; 
  params.pulsar.refTime.gpsNanoSeconds = timeV.data[0].gpsNanoSeconds;

  params.pulsar.position.longitude = pulsarInject.longitude;
  params.pulsar.position.latitude = pulsarInject.latitude ;
  params.pulsar.position.system = COORDINATESYSTEM_EQUATORIAL; 
  params.pulsar.psi = pulsarInject.psi;
  params.pulsar.aPlus = pulsarInject.aPlus;
  params.pulsar.aCross = pulsarInject.aCross;
  params.pulsar.phi0 = pulsarInject.phi0;
  params.pulsar.f0 = pulsarInject.f0;
  params.pulsar.spindown = &pulsarInject.spindown ;

  LAL_CALL( LALGeneratePulsarSignal(&status, &signalTseries, &params ), &status);
  LAL_CALL( LALSignalToSFTs(&status, &outputSFTs, signalTseries, &sftParams), &status); 
  
  /* fill in elements of sft structure sft1 used in peak selection */
  sft1.length = sftlength;
  sft1.fminBinIndex = sftFminBin;
  sft1.timeBase = timeBase;


  /* *********************************************************************** */
  /* loop over mismatched templates */
  for (mmT = 0; mmT <= 0; mmT++)
    { 
      for (mmP = 0; mmP <= 0; mmP++)
	{
	  INT4 mmFactor;
	  
	  /* displace the template */
	  mmFactor = 0;
	  pulsarTemplate1.f0 = pulsarTemplate.f0 /*+ mmFactor * mm * deltaF*/; 
	  pulsarTemplate1.latitude = pulsarTemplate.latitude + mmFactor * mmT * deltaTheta;
	  pulsarTemplate1.longitude = pulsarTemplate.longitude + mmFactor * mmP * deltaTheta;
	  pulsarTemplate1.spindown.data[0] = pulsarTemplate.spindown.data[0] /*+ mmFactor * mm * deltaFdot*/;
	 
	  /* initialize number counts */ 
	  numberCount1 = 0.0; /* usual number count */
	  numberCount2 = 0.0; /* number count with amplitude modulation */
	  /* meanAM = 0.0; */  /* mean and variance for new distribution */
	  /* sigmaAM = 0.0; */
	  /* now calculate the number count for the template */
	  for (j=0; j < mObsCoh; j++)  
	    {
	      INT4 ind;
	      REAL8 tempW;
	      /* REAL8 realThrAM; */
	      
	      sft1.epoch.gpsSeconds = timeV.data[j].gpsSeconds;
	      sft1.epoch.gpsNanoSeconds = timeV.data[j].gpsNanoSeconds;
	      sft1.data = outputSFTs->data[j].data->data;
	      
	      LAL_CALL( COMPLEX8SFT2Periodogram1(&status, &periPSD.periodogram, &sft1), &status );	
	      
	      LAL_CALL( LALPeriodo2PSDrng( &status, 
				      &periPSD.psd, &periPSD.periodogram, &uvar_blocksRngMed), &status );	
	      
	   
	      /* construct peakgram with usual threshold */
	      LAL_CALL( LALSelectPeakColorNoise(&status,&pg1,&threshold,&periPSD), &status); 	 
	      /*LAL_CALL( LALSelectPeakColorNoise(&status,&pg2,&thresholdAM,&periPSD), &status); */


	      /*adjust threshold for amplitude modulation */
	      /* 	      a = aVec->data[j]; */
	      /* 	      b = bVec->data[j]; */
	      /* 	      thresholdAM = threshold * 0.5 * (A+B) / ( a*a + b*b); */

	   
	      LAL_CALL( ComputeFoft(&status, &foft, &pulsarTemplate1, &timeDiffV, &velV, timeBase), &status);
	      
	      ind = floor( foft.data[j]*timeBase - sftFminBin + 0.5); 
	      
	      tempW = weight->data[j];
	      numberCount1 += pg1.data[ind]; 
	      numberCount2 += tempW * pg1.data[ind];
	      
	      /* 	      realThrAM = thresholdAM * normalizeThr; */
	      
	      /* 	      meanAM += exp(-realThrAM);	       */
	      /* 	      sigmaAM += exp(-realThrAM) * (1.0 - exp(-realThrAM)); */
	    } /* end loop over sfts */

	  fprintf(stdout, "%g   %g   %g  %g    %g   %g    %g    %g    %g   %g    %g\n", 
	          uvar_alpha, uvar_delta, uvar_cosiota, mean1, var1, nth1, numberCount1, mean2, var2, nth2, numberCount2);

	  
	  /* calculate the number count thresholds */
	  /* 	  mean = mObsCoh * exp(-uvar_peakThreshold); */
	  /* 	  sigma = mean * ( 1.0 - exp(-uvar_peakThreshold)); */
	  
	  /* 	  uvar_houghFalseAlarm = 1.0e-10; */
	  /* 	  erfcInv = gsl_cdf_ugaussian_Qinv (uvar_houghFalseAlarm)/sqrt(2);     */
	  /* 	  nth1 = mean + sqrt( 2 * sigma ) * erfcInv;  */
	  /* 	  /\* use mean and variance to get approximate threshold in Gaussian approximation *\/ */
	  /* 	  nth2 = meanAM + sqrt( 2 * sigmaAM ) * erfcInv;  */
	  /* 	  /\* 	  fprintf(stdout, "%g    %g    %d    %d    %g    %g    %d    %g    %g    %g    %g\n",  *\/ */
	  /* 	  /\* 	          uvar_alpha, uvar_delta, numberCount1, numberCount2,  *\/ */
	  /* 	  /\* 	          nth1, nth2, mObsCoh, uvar_peakThreshold,  *\/ */
	  /* 	  /\* 	          meanAM, sqrt(sigmaAM), varAM ); *\/ */
	  /* 	  fprintf(stdout, "%g    %g    %d    %d    %g    %g    %g    %g \n",  */
	  /* 	          uvar_alpha, uvar_delta, numberCount1, numberCount2,  */
	  /* 	          nth1, nth2, (numberCount1 -mean)/sqrt(sigma), (numberCount2 - meanAM)/sqrt(sigmaAM) ); */

	  

	}/* end loop over mmP */
    }/* end loop over mmT */
  


  /* *********************************************************************** */
  /* free structures created by signal generation routines */
  LALFree(signalTseries->data->data);
  LALFree(signalTseries->data);
  LALFree(signalTseries);
  signalTseries =NULL;
  LAL_CALL(LALDestroySFTVector(&status, &outputSFTs),&status );

  /* destroy input sfts */
  LAL_CALL(LALDestroySFTVector(&status, &inputSFTs),&status );

  /* free other structures */
  LALFree(foft.data);  
  LALFree(pulsarInject.spindown.data);
  LALFree(pulsarTemplate.spindown.data);
  LALFree(pulsarTemplate1.spindown.data);
  LALFree(timeV.data);
  LALFree(timeDiffV.data);
  LALFree(velV.data);
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);
  LALFree(periPSD.periodogram.data);
  LALFree(periPSD.psd.data);

  LALFree(pg1.data);
  /*   LALFree(pg2.data); */

  /* free amParams */
  LAL_CALL( LALDDestroyVector( &status, &weight), &status);

  /*   LAL(amc.a->data); */
  /*   LALFree(amc.b->data); */

  LAL_CALL (LALDestroyUserVars(&status), &status);  
  LALCheckMemoryLeaks();
  
  INFO( DRIVEHOUGHCOLOR_MSGENORM );
  return DRIVEHOUGHCOLOR_ENORM;
}


/******************************************************************/
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
