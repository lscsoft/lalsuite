/****************************************************
 *
 * Filename: HoughMismatch.c
 *
 * Revision: $Id$
 * 
 * Estimating mismatch of grid used in Hough search
 * 
 ****************************************************/


#include "../src/MCInjectHoughS2.h"
#include <lal/LALComputeAM.h>
#include <lal/ComputeSky.h>
#include <gsl/gsl_cdf.h>

INT4 lalDebugLevel;

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
#define SFTDIRECTORY "/home/badkri/fakesfts/"
#define FILEOUT "./MismatchOut"   
#define TRUE (1==1)
#define FALSE (1==0)

NRCSID (HOUGHMISMATCHC, "$Id$");

int main( int argc, char *argv[]){

  static LALStatus            status;  
  static LALDetector          detector;
  static LIGOTimeGPSVector    timeV;
  static REAL8Cart3CoorVector velV;
  static REAL8Vector          timeDiffV;
  static REAL8Vector          foft;
  static PulsarSignalParams  params;
  static SFTParams sftParams;

  static UCHARPeakGram     pg1, pg2;
  static COMPLEX8SFTData1  sft1;
  static REAL8PeriodoPSD   periPSD;

  REAL4TimeSeries   *signalTseries = NULL;
  SFTVector    *inputSFTs  = NULL;  
  SFTVector    *outputSFTs = NULL;
  /* data about injected signal */
  static PulsarData           pulsarInject;

  /* the template */
  static HoughTemplate  pulsarTemplate, pulsarTemplate1;

  FILE  *fpOUT = NULL; /* output file pointer */
  FILE  *fpLog = NULL; /* log file pointer */
  CHAR  *logstr=NULL; /* log string containing user input variables */
  CHAR  *fnamelog=NULL; /* name of log file */
  INT4  nfSizeCylinder;
  
  EphemerisData   *edat = NULL;

  INT4   mObsCoh, j, numberCount1, numberCount2;
  REAL8  nth1, nth2, erfcInv, houghFalseAlarm, meanAM, sigmaAM;
  REAL8  sftBand;  
  REAL8  timeBase, deltaF, normalizeThr, threshold, thresholdAM;
  UINT4  sftlength; 
  INT4   sftFminBin;
  REAL8  fHeterodyne;
  REAL8  tSamplingRate;

  /* grid spacings */
  REAL8 deltaTheta, deltaFdot;
  INT4 mmP, mmT; /* for loop over mismatched templates */

  /* user input variables */
  BOOLEAN uvar_help;
  INT4 uvar_ifo, uvar_blocksRngMed;
  REAL8 uvar_peakThreshold;
  REAL8 uvar_alpha, uvar_delta, uvar_h0, uvar_f0;
  REAL8 uvar_psi, uvar_phi0, uvar_fdot, uvar_cosiota;
  CHAR *uvar_earthEphemeris=NULL;
  CHAR *uvar_sunEphemeris=NULL;
  CHAR *uvar_sftDir=NULL;
  CHAR *uvar_fnameout=NULL;

  /* amplitude modulation stuff */
  REAL4Vector *aVec, *bVec;
  REAL8 A, B, C, D, a, b;
  AMCoeffs amc; 
  AMCoeffsParams *amParams;
  EarthState earth;
  BarycenterInput baryinput;

  /*  set up the default parameters  */
  lalDebugLevel = 0;

  nfSizeCylinder = NFSIZE;
  /* LALDebugLevel must be called before anything else */
  SUB( LALGetDebugLevel( &status, argc, argv, 'd'), &status);

  /* set other user input variables */
  uvar_help = FALSE;
  uvar_peakThreshold = THRESHOLD;
  uvar_ifo = IFO;
  uvar_blocksRngMed = BLOCKSRNGMED;

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

  /* register user input variables */
  SUB( LALRegisterBOOLUserVar(   &status, "help",            'h', UVAR_HELP,     "Print this message",            &uvar_help),            &status);  
  SUB( LALRegisterINTUserVar(    &status, "ifo",             'i', UVAR_OPTIONAL, "Detector GEO(1) LLO(2) LHO(3)", &uvar_ifo ),            &status);
  SUB( LALRegisterINTUserVar(    &status, "blocksRngMed",    'w', UVAR_OPTIONAL, "RngMed block size",             &uvar_blocksRngMed),    &status);
  SUB( LALRegisterREALUserVar(   &status, "peakThreshold",   't', UVAR_OPTIONAL, "Peak selection threshold",      &uvar_peakThreshold),   &status);
  SUB( LALRegisterSTRINGUserVar( &status, "earthEphemeris",  'E', UVAR_OPTIONAL, "Earth Ephemeris file",          &uvar_earthEphemeris),  &status);
  SUB( LALRegisterSTRINGUserVar( &status, "sunEphemeris",    'S', UVAR_OPTIONAL, "Sun Ephemeris file",            &uvar_sunEphemeris),    &status);
  SUB( LALRegisterSTRINGUserVar( &status, "sftDir",          'D', UVAR_OPTIONAL, "SFT Directory",                 &uvar_sftDir),          &status);
  SUB( LALRegisterSTRINGUserVar( &status, "fnameout",        'o', UVAR_OPTIONAL, "Output file prefix",            &uvar_fnameout),        &status);
  SUB( LALRegisterREALUserVar(   &status, "alpha",           'r', UVAR_OPTIONAL, "Right ascension",               &uvar_alpha),           &status);
  SUB( LALRegisterREALUserVar(   &status, "delta",           'l', UVAR_OPTIONAL, "Declination",                   &uvar_delta),           &status);
  SUB( LALRegisterREALUserVar(   &status, "h0",              'm', UVAR_OPTIONAL, "h0 to inject",                  &uvar_h0),              &status);
  SUB( LALRegisterREALUserVar(   &status, "f0",              'f', UVAR_OPTIONAL, "Start search frequency",        &uvar_f0),              &status);
  SUB( LALRegisterREALUserVar(   &status, "psi",             'p', UVAR_OPTIONAL, "Polarization angle",            &uvar_psi),             &status);
  SUB( LALRegisterREALUserVar(   &status, "phi0",            'P', UVAR_OPTIONAL, "Initial phase",                 &uvar_phi0),            &status);
  SUB( LALRegisterREALUserVar(   &status, "cosiota",         'c', UVAR_OPTIONAL, "Cosine of iota",                &uvar_cosiota),         &status);
  SUB( LALRegisterREALUserVar(   &status, "fdot",            'd', UVAR_OPTIONAL, "Spindown parameter",            &uvar_fdot),            &status);

  /* read all command line variables */
  SUB( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 
  
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
  SUB( LALUserVarGetLog(&status, &logstr, UVAR_LOGFMT_CFGFILE), &status);  

  fprintf( fpLog, "## Log file for HoughMismatch\n\n");
  fprintf( fpLog, "# User Input:\n");
  fprintf( fpLog, "#-------------------------------------------\n");
  fprintf( fpLog, logstr);
  LALFree(logstr);

  /* append an ident-string defining the exact CVS-version of the code used */
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

  /* set peak selection threshold */
  SUB( LALRngMedBias( &status, &normalizeThr, uvar_blocksRngMed ), &status ); 
  threshold = uvar_peakThreshold/normalizeThr; 

  /* set detector */
  if (uvar_ifo ==1) detector=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if (uvar_ifo ==2) detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
  if (uvar_ifo ==3) detector=lalCachedDetectors[LALDetectorIndexLHODIFF];


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

  /* copy these values also to the pulsar template */
  /* template is complately matched at this point */
  pulsarTemplate.f0 = uvar_f0;
  pulsarTemplate.latitude = uvar_delta;
  pulsarTemplate.longitude = uvar_alpha;
  pulsarTemplate.spindown.length = 1;
  pulsarTemplate.spindown.data = NULL;
  pulsarTemplate.spindown.data = (REAL8 *)LALMalloc(sizeof(REAL8));
  pulsarTemplate.spindown.data[0] = uvar_fdot;

  /* allocate memory for mismatched spindown template */
  pulsarTemplate1.spindown.length = 1;
  pulsarTemplate1.spindown.data = NULL;
  pulsarTemplate1.spindown.data = (REAL8 *)LALMalloc(sizeof(REAL8));

  /* read sfts */
  {
    CHAR *tempDir;
    tempDir = (CHAR *)LALMalloc(512*sizeof(CHAR));
    strcpy(tempDir, uvar_sftDir);
    strcat(tempDir, "/*SFT*.*"); 
    sftBand = 0.5; 
    SUB( LALReadSFTfiles ( &status, &inputSFTs, uvar_f0 - sftBand, uvar_f0 + sftBand, nfSizeCylinder + uvar_blocksRngMed , tempDir), &status);
    LALFree(tempDir);
  }


  /* get sft parameters */
  mObsCoh = inputSFTs->length;
  sftlength = inputSFTs->data->data->length;
  deltaF = inputSFTs->data->deltaF;
  timeBase = 1.0/deltaF;
  sftFminBin = floor( timeBase * inputSFTs->data->f0 + 0.5);
  fHeterodyne = sftFminBin*deltaF;
  tSamplingRate = 2.0*deltaF*(sftlength -1.);

  /* create timestamp vector */
  timeV.length = mObsCoh;
  timeV.data = NULL;  
  timeV.data = (LIGOTimeGPS *)LALMalloc(mObsCoh*sizeof(LIGOTimeGPS));

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

  /* setting of ephemeris info */ 
  edat = (EphemerisData *)LALMalloc(sizeof(EphemerisData));
  (*edat).ephiles.earthEphemeris = uvar_earthEphemeris;
  (*edat).ephiles.sunEphemeris = uvar_sunEphemeris;
  
  /* compute detector velocity for those time stamps  */ 
  velV.length = mObsCoh; 
  velV.data = NULL;
  velV.data = (REAL8Cart3Coor *)LALMalloc(mObsCoh*sizeof(REAL8Cart3Coor));
  
  {  
    VelocityPar   velPar;
    REAL8     vel[3]; 
    UINT4     j; 
    
    LALLeapSecFormatAndAcc lsfas = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
    INT4 tmpLeap; /* need this because Date pkg defines leap seconds as
		     INT4, while EphemerisData defines it to be INT2. This won't
                   cause problems before, oh, I don't know, the Earth has been 
                   destroyed in nuclear holocaust. -- dwchin 2004-02-29 */
    
    velPar.detector = detector;
    velPar.tBase = timeBase;
    velPar.vTol = 0.0; /* irrelevant */
    velPar.edat = NULL;
    
    /* Leap seconds for the start time of the run */   
    SUB( LALLeapSecs(&status, &tmpLeap, &(timeV.data[0]), &lsfas), &status);
    (*edat).leap = (INT2)tmpLeap;
    
    /* read in ephemeris data */
    SUB( LALInitBarycenter( &status, edat), &status);
    velPar.edat = edat;

    /* calculate detector velocity */    
    for(j=0; j< velV.length; ++j){
      velPar.startTime.gpsSeconds     = timeV.data[j].gpsSeconds;
      velPar.startTime.gpsNanoSeconds = timeV.data[j].gpsNanoSeconds;
      
      SUB( LALAvgDetectorVel ( &status, vel, &velPar), &status );
      velV.data[j].x= vel[0];
      velV.data[j].y= vel[1];
      velV.data[j].z= vel[2];   
    }  
  }


  /* amplitude modulation stuff */


  /* detector location */
  baryinput.site.location[0] = detector.location[0]/LAL_C_SI;
  baryinput.site.location[1] = detector.location[1]/LAL_C_SI;
  baryinput.site.location[2] = detector.location[2]/LAL_C_SI;
  baryinput.dInv = 0.e0;
  /* alpha and delta must come from the skypatch */
  /* for now set it to something arbitrary */
  baryinput.alpha = 0.0;
  baryinput.delta = 0.0;
  
  /* Allocate space for amParams stucture */
  /* Here, amParams->das is the Detector and Source info */
  amParams = (AMCoeffsParams *)LALMalloc(sizeof(AMCoeffsParams));
  amParams->das = (LALDetAndSource *)LALMalloc(sizeof(LALDetAndSource));
  amParams->das->pSource = (LALSource *)LALMalloc(sizeof(LALSource));
  /* Fill up AMCoeffsParams structure */
  amParams->baryinput = &baryinput;
  amParams->earth = &earth; 
  amParams->edat = edat;
  amParams->das->pDetector = &detector; 
  /* make sure alpha and delta are correct */
  amParams->das->pSource->equatorialCoords.latitude = baryinput.delta;
  amParams->das->pSource->equatorialCoords.longitude = baryinput.alpha;
  amParams->das->pSource->orientation = 0.0;
  amParams->das->pSource->equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  amParams->polAngle = amParams->das->pSource->orientation ; /* These two have to be the same!!*/
  amParams->leapAcc = LALLEAPSEC_STRICT;

  /* allocate memory for a[i] and b[i] */
  amc.a = NULL;
  amc.b = NULL;
  SUB( LALSCreateVector( &status, &(amc.a), mObsCoh), &status);
  SUB( LALSCreateVector( &status, &(amc.b), mObsCoh), &status);
  /*   amc.a->length = mObsCoh; */
  /*   amc.a->data = (REAL4 *)LALMalloc(mObsCoh*sizeof(REAL4)); */
  /*   amc.b->length = mObsCoh; */
  /*   amc.b->data = (REAL4 *)LALMalloc(mObsCoh*sizeof(REAL4)); */

  /* timeV is start time ---> change to mid time */    
  SUB (LALComputeAM(&status, &amc, timeV.data, amParams), &status); 
  aVec = amc.a; /* a and b are as defined in JKS */
  bVec = amc.b;
  A = amc.A; /* note A is twice average of a[i]^2 */ 
  B = amc.B; /* note B is twice average of b[i]^2 */ 
  C = amc.C;
  D = amc.D;
  

  /* set grid spacings */
  {
    REAL8 vel2;

    /* use start velocity to set delta Theta */    
    vel2 = velV.data[0].x * velV.data[0].x + velV.data[0].y * velV.data[0].y + velV.data[0].z * velV.data[0].z;
    
    deltaTheta = 1.0 / ( VTOT * uvar_f0 * timeBase );
    deltaFdot = deltaF / timeBase;
  }

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

  pg2.length = sftlength;
  pg2.data = NULL;
  pg2.data = (UCHAR *)LALMalloc(sftlength* sizeof(UCHAR));

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
  params.pulsar.tRef.gpsSeconds = timeV.data[0].gpsSeconds; 
  params.pulsar.tRef.gpsNanoSeconds = timeV.data[0].gpsNanoSeconds;

  params.pulsar.position.longitude = pulsarInject.longitude;
  params.pulsar.position.latitude = pulsarInject.latitude ;
  params.pulsar.position.system = COORDINATESYSTEM_EQUATORIAL; 
  params.pulsar.psi = pulsarInject.psi;
  params.pulsar.aPlus = pulsarInject.aPlus;
  params.pulsar.aCross = pulsarInject.aCross;
  params.pulsar.phi0 = pulsarInject.phi0;
  params.pulsar.f0 = pulsarInject.f0;
  params.pulsar.spindown = &pulsarInject.spindown ;

  SUB( LALGeneratePulsarSignal(&status, &signalTseries, &params ), &status);
  SUB( LALSignalToSFTs(&status, &outputSFTs, signalTseries, &sftParams), &status); 
  
  /* fill in elements of sft structure sft1 used in peak selection */
  sft1.length = sftlength;
  sft1.fminBinIndex = sftFminBin;
  sft1.timeBase = timeBase;


  /* loop over mismatched templates */
  for (mmT = 0; mmT <= 0; mmT++)
    { 
      for (mmP = 0; mmP <= 0; mmP++)
	{
	  INT4 mmFactor;
	  
	  /* displace the template */
	  mmFactor = 1.0;
	  pulsarTemplate1.f0 = pulsarTemplate.f0 /*+ mmFactor * mm * deltaF*/; 
	  pulsarTemplate1.latitude = pulsarTemplate.latitude + mmFactor * mmT * deltaTheta;
	  pulsarTemplate1.longitude = pulsarTemplate.longitude + mmFactor * mmP * deltaTheta;
	  pulsarTemplate1.spindown.data[0] = pulsarTemplate.spindown.data[0] /*+ mmFactor * mm * deltaFdot*/;
	  
	  numberCount1 = 0; /* usual number count */
	  numberCount2 = 0; /* number count with amplitude modulation */
	  meanAM = 0.0;  /* mean and variance for new distribution */
	  sigmaAM = 0.0;
	  /* now calculate the number count for the template */
	  for (j=0; j < mObsCoh; j++)  
	    {
	      INT4 index;
	      
	      sft1.epoch.gpsSeconds = timeV.data[j].gpsSeconds;
	      sft1.epoch.gpsNanoSeconds = timeV.data[j].gpsNanoSeconds;
	      sft1.data = outputSFTs->data[j].data->data;
	      
	      SUB( COMPLEX8SFT2Periodogram1(&status, &periPSD.periodogram, &sft1), &status );	
	      
	      SUB( LALPeriodo2PSDrng( &status, 
				      &periPSD.psd, &periPSD.periodogram, &uvar_blocksRngMed), &status );	
	      
	   
	      /* construct peakgram with usual threshold */
	      SUB( LALSelectPeakColorNoise(&status,&pg1,&threshold,&periPSD), &status); 	 

	      /*adjust threshold for amplitude modulation */
	      a = aVec->data[j];
	      b = bVec->data[j];
	      thresholdAM = threshold * 0.5 * (A+B) / ( a*a + b*b);
	      SUB( LALSelectPeakColorNoise(&status,&pg2,&thresholdAM,&periPSD), &status); 	 
	   
	      SUB( ComputeFoft(&status, &foft, &pulsarTemplate1, &timeDiffV, &velV, timeBase), &status);
	      
	      index = floor( foft.data[j]*timeBase - sftFminBin + 0.5); 
	      
	      numberCount1 += pg1.data[index]; 
	      numberCount2 += pg2.data[index]; 

	      meanAM += exp(-thresholdAM);	      
	      sigmaAM += exp(-thresholdAM) * (1.0 - exp(-thresholdAM));
	    } 
	  /* print the number count */

	  /* calculate the number count thresholds */
	  houghFalseAlarm = 1.0e-10;
	  erfcInv = gsl_cdf_ugaussian_Qinv (houghFalseAlarm)/sqrt(2);    
	  nth1 = mObsCoh* exp(-threshold) + sqrt(2*mObsCoh*exp(-threshold)*(1-exp(-threshold)))*erfcInv; 
	  /* use mean and variance to get approximate threshold in Gaussian approximatio */
	  nth2 = meanAM + sqrt( 2 * sigmaAM )* erfcInv; 
	  fprintf(stdout, "%d    %d    %g    %g\n", numberCount1, numberCount2, nth1, nth2);
	}
    }
  
  /* free structures created by signal generation routines */
  LALFree(signalTseries->data->data);
  LALFree(signalTseries->data);
  LALFree(signalTseries);
  signalTseries =NULL;
  SUB(LALDestroySFTVector(&status, &outputSFTs),&status );

  /* destroy input sfts */
  SUB(LALDestroySFTVector(&status, &inputSFTs),&status );

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
  LALFree(pg2.data);

  /* free amParams */
  LALFree(amParams->das->pSource);
  LALFree(amParams->das);
  LALFree(amParams);
  SUB( LALSDestroyVector( &status, &(amc.a)), &status);
  SUB( LALSDestroyVector( &status, &(amc.b)), &status);

  /*   LAL(amc.a->data); */
  /*   LALFree(amc.b->data); */

  SUB (LALDestroyUserVars(&status), &status);  
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
  INITSTATUS (status, "ComputeFoft", HOUGHMISMATCHC);
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
