 /*-----------------------------------------------------------------------
 *
 * File Name: makefakedata_v4.c
 *
 * Authors: R. Prix, M.A. Papa, X. Siemens, B. Allen, C. Messenger
 *
 * Revision: $Id$
 *
 *           
 *-----------------------------------------------------------------------
 */

/* ---------- includes ---------- */
#include <lalapps.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALInitBarycenter.h>
#include <lal/Random.h>

#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/GeneratePulsarSignal.h>


RCSID ("$Id$");

/* Error codes and messages */
/* <lalErrTable file="MAKEFAKEDATACErrorTable"> */
#define MAKEFAKEDATAC_ENORM 	0
#define MAKEFAKEDATAC_ESUB  	1
#define MAKEFAKEDATAC_EARG  	2
#define MAKEFAKEDATAC_EBAD  	3
#define MAKEFAKEDATAC_EFILE 	4
#define MAKEFAKEDATAC_ENOARG 	5

#define MAKEFAKEDATAC_MSGENORM "Normal exit"
#define MAKEFAKEDATAC_MSGESUB  "Subroutine failed"
#define MAKEFAKEDATAC_MSGEARG  "Error parsing arguments"
#define MAKEFAKEDATAC_MSGEBAD  "Bad argument values"
#define MAKEFAKEDATAC_MSGEFILE "File IO error"
#define MAKEFAKEDATAC_MSGENOARG "Missing argument"

/* </lalErrTable> */

/***************************************************/
#define TRUE (1==1)
#define FALSE (1==0)

/*----------------------------------------------------------------------*/
/* config-variables derived from user-variables */
typedef struct 
{
  EphemerisData edat;
  LALDetector Detector;  
  LIGOTimeGPS startTime;
  LIGOTimeGPS refTime;
  LIGOTimeGPSVector *timestamps;
  REAL8 duration;
  REAL8 fmin_eff;	/* 'effective' fmin: round down such that fmin*Tsft = integer! */
  REAL8Vector *spindown;
  SFTVector *noiseSFTs;
  BOOLEAN binaryPulsar;
  LIGOTimeGPS orbitPeriTime;
} ConfigVars_t;

/* Locations of the earth and sun ephemeris data */
#define EPHEM_YEARS  "00-04"

/* local prototypes */
/* Prototypes for the functions defined in this file */
void FreeMem (LALStatus* stat, ConfigVars_t *cfg);
void InitUserVars (LALStatus *stat);
void ReadTimestamps (LALStatus* stat, LIGOTimeGPSVector **timestamps, const CHAR *fname);
void InitMakefakedata (LALStatus *stat, ConfigVars_t *cfg, int argc, char *argv[]);
void AddGaussianNoise (LALStatus* status, REAL4TimeSeries *outSeries, REAL4TimeSeries *inSeries, REAL4 sigma);
void GetOrbitalParams (LALStatus *stat, BinaryOrbitParams *orbit);

extern void write_timeSeriesR4 (FILE *fp, const REAL4TimeSeries *series);
extern void write_timeSeriesR8 (FILE *fp, const REAL8TimeSeries *series);

/*----------------------------------------------------------------------*/
static const PulsarSignalParams empty_params;
static const SFTParams empty_sftParams;
static const EphemerisData empty_edat;
static const ConfigVars_t empty_GV;
/*----------------------------------------------------------------------*/
/* User variables */
/*----------------------------------------------------------------------*/
CHAR *uvar_outSFTbname;
CHAR *uvar_outTDDFile;
CHAR *uvar_timestampsFile;
CHAR *uvar_detector;
REAL8 uvar_startTime;
REAL8 uvar_refTime;
CHAR *uvar_ephemDir;
CHAR *uvar_ephemYear;
CHAR *uvar_noiseSFTs;
BOOLEAN uvar_doWindowing;
BOOLEAN uvar_binaryoutput;
BOOLEAN uvar_nomagic;
BOOLEAN uvar_help;
REAL8 uvar_Tsft;
INT4 uvar_nTsft;
REAL8 uvar_fmin;
REAL8 uvar_Band;
REAL8 uvar_aPlus;
REAL8 uvar_aCross;
REAL8 uvar_psi;
REAL8 uvar_phi0;
REAL8 uvar_f0;
REAL8 uvar_longitude;
REAL8 uvar_latitude;
REAL8 uvar_f1dot;
REAL8 uvar_f2dot;
REAL8 uvar_f3dot;

/* orbital parameters */
REAL8 uvar_orbitSemiMajorAxis;
REAL8 uvar_orbitEccentricity;
INT4  uvar_orbitTperiSSBsec;
INT4  uvar_orbitTperiSSBns;
REAL8 uvar_orbitPeriod;
REAL8 uvar_orbitArgPeriapse;

/*----------------------------------------------------------------------*/


/* Gaussian noise */
REAL8 uvar_noiseSigma;

extern int vrbflg;

/*----------------------------------------------------------------------
 * main function 
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[]) 
{
  LALStatus status = blank_status;	/* initialize status */
  ConfigVars_t GV = empty_GV;
  PulsarSignalParams params = empty_params;
  SFTParams sftParams = empty_sftParams;
  SFTVector *SFTs = NULL;
  REAL4TimeSeries *Tseries = NULL;
  BinaryOrbitParams orbit;

  CHAR *fname = NULL;
  UINT4 i;

  lalDebugLevel = 0;
  vrbflg = 1;		/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;	/* exit with returned status-code on error */

  /* ------------------------------ 
   * read user-input and set up shop
   *------------------------------*/
  LAL_CALL (LALGetDebugLevel (&status, argc, argv, 'v'), &status);

  LAL_CALL (InitMakefakedata (&status, &GV, argc, argv), &status);

  if (uvar_help)	/* help was called, do nothing here.. */
    return (0);

  /*----------------------------------------
   * fill the PulsarSignalParams struct 
   *----------------------------------------*/
  /* pulsar params */
  params.pulsar.tRef = GV.refTime;
  params.pulsar.position.longitude = uvar_longitude;
  params.pulsar.position.latitude  = uvar_latitude;
  params.pulsar.position.system    = COORDINATESYSTEM_EQUATORIAL;
  params.pulsar.psi 		   = uvar_psi;
  params.pulsar.aPlus		   = uvar_aPlus;
  params.pulsar.aCross		   = uvar_aCross;
  params.pulsar.phi0		   = uvar_phi0;
  params.pulsar.f0		   = uvar_f0;

  /* orbital params */
  if (GV.binaryPulsar)
    {
      LAL_CALL ( GetOrbitalParams (&status, &orbit), &status);
      params.orbit = &orbit;
    }
  else
    params.orbit = NULL;

  /* detector params */
  params.transfer = NULL;	/* detector transfer function (NULL if not used) */	
  params.site = &(GV.Detector);	
  params.ephemerides = &(GV.edat);

  /* characterize the output time-series */
  params.startTimeGPS   = GV.startTime;
  params.duration     	= (UINT4) ceil(GV.duration); /* length of time-series in seconds */
  params.samplingRate 	= 2.0 * uvar_Band;	/* sampling rate of time-series (= 2 * frequency-Band) */
  params.fHeterodyne  	= GV.fmin_eff;		/* heterodyning frequency for output time-series */

  /*================================================================================*/
  /* test new function SimulatePulsarSignal() */
#if 0
  {
    REAL8 LMST, GMST;	/* local mean sidereal time */
    LALPlaceAndGPS place_and_gps;
    LALMSTUnitsAndAcc units_and_acc;
    LALDate date;
    CHARVector *dateString = NULL;
    LIGOTimeGPS gps;
    LALFrDetector det = params.site->frDetector;
    REAL8 zeta;
    REAL8TimeSeries *Tseries8 = NULL;
  
    LAL_CALL ( LALCHARCreateVector (&status, &dateString, 64), &status);	/* ARGHHHH ..... */

    units_and_acc.accuracy = LALLEAPSEC_STRICT;
    units_and_acc.units =   MST_DEG;	/* return-format of LMST */
 
    place_and_gps.p_detector = params.site;
    place_and_gps.p_gps = &(params.startTimeGPS);

    /*
    LAL_CALL ( LALDateString (&status, dateString, &date), &status);
    LAL_CALL (LALUTCtoGPS (&status, &gps, &date, &(units_and_acc.accuracy)), &status);
    */
                                                    
    LAL_CALL (LALGPStoLMST1 (&status, &LMST, &place_and_gps, &units_and_acc), &status);

    LAL_CALL (LALGPStoGMST1 (&status, &GMST, &(params.startTimeGPS), &units_and_acc), &status);

    LAL_CALL (LALGPStoUTC (&status, &date, &(params.startTimeGPS), &(units_and_acc.accuracy) ), &status);
    LAL_CALL ( LALDateString (&status, dateString, &date), &status);
    printf ("DEBUG: GPS= %d, UTC = %s\n", params.startTimeGPS.gpsSeconds, dateString->data);
    LAL_CALL ( LALDateString (&status, dateString, &date), &status);

    printf ("DEBUG: Start-time: GPS= %d, UTC= %s; long=%f deg (x=%f, y=%f), LMST = %f, GMST = %f\n", 
	    params.startTimeGPS.gpsSeconds, dateString->data, 
	    (REAL8) (params.site->frDetector.vertexLongitudeRadians * (REAL8)LAL_180_PI), 
	    params.site->location[0], params.site->location[1],
	    LMST, GMST );

    printf ("atan2 (y, x) = %f deg\n", atan2 (params.site->location[1], params.site->location[0]) * (REAL8)LAL_180_PI);


    zeta = 1.0/(sin(det.xArmAzimuthRadians - det.yArmAzimuthRadians));
    if(params.site->type == LALDETECTORTYPE_CYLBAR) zeta=1.0;
    printf ("\nDEBUG: xArmAzimuth - yArmAzimut = %.10f, sin zeta = %f\n",
	  det.xArmAzimuthRadians -  det.yArmAzimuthRadians, zeta );

    LAL_CALL (LALSimulatePulsarSignal (&status, &Tseries8, &params), &status);

  }
#endif
  /*================================================================================*/


  /*----------------------------------------
   * generate the heterodyned time-series 
   *----------------------------------------*/
  LAL_CALL (LALGeneratePulsarSignal (&status, &Tseries, &params), &status );

  if (lalDebugLevel >= 3)
    {  
      FILE *fp;
      fp = fopen ("Tseries_v4.dat", "w");
      write_timeSeriesR4 (fp, Tseries);
      fclose (fp);
    }

  /* add Gaussian noise if requested */
  if ( (REAL4)uvar_noiseSigma > 0)
    {
      LAL_CALL ( AddGaussianNoise (&status, Tseries, Tseries, (REAL4)uvar_noiseSigma), &status);
    }

  /*----------------------------------------
   * last step: turn this timeseries into SFTs
   *----------------------------------------*/
  sftParams.Tsft = uvar_Tsft;
  sftParams.timestamps = GV.timestamps;
  sftParams.noiseSFTs = GV.noiseSFTs;

  LAL_CALL ( LALSignalToSFTs (&status, &SFTs, Tseries, &sftParams), &status);


  /* write SFTs to disk */
  if (uvar_outSFTbname)
    {
      fname = LALCalloc (1, strlen (uvar_outSFTbname) + 10);
      for (i=0; i < SFTs->length; i++)
	{
	  sprintf (fname, "%s.%05d", uvar_outSFTbname, i);
	  LAL_CALL (LALWriteSFTfile (&status, &(SFTs->data[i]), fname), &status);
	} /* for i < nSFTs */
      LALFree (fname);
    } /* if SFTbname */

  /* free memory */
  if (Tseries) 
    {
      LAL_CALL (LALDestroyVector (&status, &(Tseries->data)), &status);
      LALFree (Tseries);
      Tseries = NULL;
    }
  if (SFTs) 
    LAL_CALL (LALDestroySFTVector(&status, &SFTs), &status);

  LAL_CALL (FreeMem (&status, &GV), &status);	/* free the rest */

  LALCheckMemoryLeaks(); 

  return 0;
} /* main */

/*---------------------------------------------------------------------- 
 *  handle user-input and set up shop accordingly 
 *
 * we do all consistency-checks on user-input here
 *----------------------------------------------------------------------*/
void
InitMakefakedata (LALStatus *stat, 
		  ConfigVars_t *cfg,
		  int argc, 
		  char *argv[])
{
  LALLeapSecFormatAndAcc leapParams = {  LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT };
  INT4 leapSecs;
  UINT4 msp = 0;	/* number of spindown-parameters */
  CHAR *earthdata;
  CHAR *sundata;
  EphemerisData edat = empty_edat;
  REAL8 duration;
  LIGOTimeGPSVector *timestamps = NULL;
  UINT4 imin;
  
  INITSTATUS( stat, "InitMakefakedata", rcsid );
  ATTATCHSTATUSPTR (stat);

  /* register all user-variables */
  TRY (InitUserVars (stat->statusPtr), stat);	  

  /* read cmdline & cfgfile  */	
  TRY (LALUserVarReadAllInput (stat->statusPtr, argc,argv), stat);  

  if (uvar_help) 	/* if help was requested, we're done */
    exit (0);

  /* check more complex input-dependencies */
  if ( (!LALUserVarWasSet(&uvar_timestampsFile) && !LALUserVarWasSet(&uvar_startTime))
       || (LALUserVarWasSet(&uvar_timestampsFile) && LALUserVarWasSet(&uvar_startTime)) )
    {
      LALPrintError ("\nPlease specify either timestampsFile or startTime !\n\n");
      ABORT (stat, MAKEFAKEDATAC_ENOARG, MAKEFAKEDATAC_MSGENOARG);
    }

  /* prepare vector of spindown parameters */
  if (uvar_f3dot != 0) 		msp = 3;	/* counter number of spindowns */
  else if (uvar_f2dot != 0)	msp = 2;
  else if (uvar_f1dot != 0)	msp = 1;
  else 				msp = 0;
  if (msp) {
    TRY (LALDCreateVector (stat->statusPtr, &(cfg->spindown), msp), stat);
  }
  switch (msp) 
    {
    case 3:
      cfg->spindown->data[2] = uvar_f3dot;
    case 2:
      cfg->spindown->data[1] = uvar_f2dot;
    case 1:
      cfg->spindown->data[0] = uvar_f1dot;
      break;
    case 0:
      break;
    default:
      LALPrintError ("\nmsp = %d makes no sense to me...\n\n", msp);
      ABORT (stat,  MAKEFAKEDATAC_EBAD,  MAKEFAKEDATAC_MSGEBAD);
    } /* switch(msp) */

  
  /* prepare detector */
  if      (!strcmp(uvar_detector,"LHO"))   cfg->Detector = lalCachedDetectors[LALDetectorIndexLHODIFF];
  else if (!strcmp(uvar_detector,"LLO"))   cfg->Detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  else if (!strcmp(uvar_detector,"VIRGO")) cfg->Detector = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  else if (!strcmp(uvar_detector,"GEO"))   cfg->Detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  else if (!strcmp(uvar_detector,"TAMA"))  cfg->Detector = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
  else if (!strcmp(uvar_detector,"CIT"))   cfg->Detector = lalCachedDetectors[LALDetectorIndexCIT40DIFF];
  else {
    LALPrintError ("\nUnknown detector specified: `%s\n\n", uvar_detector);
    ABORT (stat,  MAKEFAKEDATAC_EBAD,  MAKEFAKEDATAC_MSGEBAD);
  }

  /* read timestamps and set signal-duration */
  cfg->timestamps = NULL;
  if (uvar_timestampsFile) 
    {
      LIGOTimeGPS t1, t0;
      TRY (ReadTimestamps (stat->statusPtr, &timestamps, uvar_timestampsFile), stat);
      if ((UINT4)uvar_nTsft > timestamps->length) {
	LALPrintError ("\nTimestamps-file contains less than nTsft=%d entries!\n\n", uvar_nTsft);
	ABORT (stat,  MAKEFAKEDATAC_EBAD,  MAKEFAKEDATAC_MSGEBAD);
      }
      t1 = timestamps->data[uvar_nTsft-1];
      t0 = timestamps->data[0];
      TRY (LALDeltaFloatGPS(stat->statusPtr, &duration, &t1, &t0), stat);
      duration += uvar_Tsft;
    } 
  else
    duration = uvar_nTsft * uvar_Tsft;

  cfg->duration = duration;
  cfg->timestamps = timestamps;


  /* get observation start-time */
  if (LALUserVarWasSet (&uvar_startTime)) {
    TRY ( LALFloatToGPS (stat->statusPtr, &(cfg->startTime), &uvar_startTime), stat);
  }
  else
    cfg->startTime = cfg->timestamps->data[0];
      
  /* get leap-seconds since start of GPS-time */
  TRY ( LALLeapSecs (stat->statusPtr, &leapSecs,  &(cfg->startTime), &leapParams), stat);

  /* Prepare quantities for barycentering */
  earthdata = LALCalloc(1, strlen(uvar_ephemDir) + strlen(uvar_ephemYear) + 20);
  sundata = LALCalloc(1, strlen(uvar_ephemDir) + strlen(uvar_ephemYear) + 20);

  sprintf (earthdata, "%s/%s%s.dat", uvar_ephemDir, "/earth", uvar_ephemYear);
  sprintf (sundata, "%s/%s%s.dat", uvar_ephemDir, "/sun", uvar_ephemYear);
  edat.ephiles.earthEphemeris = earthdata;
  edat.ephiles.sunEphemeris   = sundata;
  edat.leap = (INT2) leapSecs;

  /* Init ephemerides */  
  TRY( LALInitBarycenter (stat->statusPtr, &edat), stat);   
  LALFree (earthdata);
  LALFree (sundata);

  cfg->edat = edat;

  /* if reference-time was given for pulsar-parameters */
  if (LALUserVarWasSet(&uvar_refTime)) {
      TRY ( LALFloatToGPS (stat->statusPtr, &(cfg->refTime), &uvar_refTime), stat);
  } else	/* otherwise set to 0, so startTime in SSB will be used */
    cfg->refTime.gpsSeconds = cfg->refTime.gpsNanoSeconds = 0;
  

  /* calculate "effective" fmin from uvar_fmin: following makefakedata_v2, we
   * make sure that fmin_eff * Tsft = integer, such that freqBinIndex corresponds
   * to a frequency-index of the non-heterodyned signal
   */
  imin = (UINT4) floor( uvar_fmin * uvar_Tsft);
  cfg->fmin_eff = (REAL8)imin / uvar_Tsft;

  /* if any orbital parameters specified, we need all of them (except for nano-seconds)! */
  {
    BOOLEAN set1 = LALUserVarWasSet(&uvar_orbitSemiMajorAxis);
    BOOLEAN set2 = LALUserVarWasSet(&uvar_orbitEccentricity);
    BOOLEAN set3 = LALUserVarWasSet(&uvar_orbitPeriod);
    BOOLEAN set4 = LALUserVarWasSet(&uvar_orbitArgPeriapse);
    BOOLEAN set5 = LALUserVarWasSet(&uvar_orbitTperiSSBsec);
    BOOLEAN set6 = LALUserVarWasSet(&uvar_orbitTperiSSBns);
    if (set1 || set2 || set3 || set4 || set5 || set6)
    {
      if ( (uvar_orbitSemiMajorAxis > 0) && !(set1 && set2 && set3 && set4 && set5) ) 
	{
	  LALPrintError ("\nPlease either specify  ALL orbital parameters or NONE!\n\n");
	  ABORT (stat, MAKEFAKEDATAC_EBAD, MAKEFAKEDATAC_MSGEBAD);
	}
      cfg->binaryPulsar = TRUE;

      if ( (uvar_orbitEccentricity < 0) || (uvar_orbitEccentricity > 1) )
	{
	  LALPrintError ("\nEccentricity has to lie within [0, 1]\n\n");
	  ABORT (stat, MAKEFAKEDATAC_EBAD, MAKEFAKEDATAC_MSGEBAD);
	}

    } /* if one or more orbital parameters were set */
  } 


  /* EITHER add Gaussian noise OR real noise-sft's */
  if ( LALUserVarWasSet(&uvar_noiseSFTs) && LALUserVarWasSet(&uvar_noiseSigma) )
    {
      LALPrintError ("\nERROR: only one of 'noiseSFTs' or 'noiseSigma' can be specified!\n\n");
      ABORT (stat, MAKEFAKEDATAC_EBAD, MAKEFAKEDATAC_MSGEBAD);
    }

  /* if real noise-SFTs: load them in now */
  if ( uvar_noiseSFTs && LALUserVarWasSet (&uvar_noiseSFTs))
    {
      REAL8 fMin, fMax;
      fMin = cfg->fmin_eff;
      fMax = fMin + uvar_Band;
      TRY ( LALReadSFTfiles (stat->statusPtr, &(cfg->noiseSFTs), fMin, fMax, 0, uvar_noiseSFTs), stat);
    } /* if uvar_noiseSFTs */


  /* catch (yet) unsupported options */
  if ( LALUserVarWasSet (&uvar_outTDDFile) )
    {
      LALPrintError ("\nSorry, time-domain output is currently not yet supported\n\n");
      ABORT (stat, MAKEFAKEDATAC_EBAD, MAKEFAKEDATAC_MSGEBAD);
    }
  
  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* InitMakefakedata() */



/*----------------------------------------------------------------------*/
/* register all our "user-variables" */
void
InitUserVars (LALStatus *stat)
{
  INITSTATUS( stat, "InitUserVars", rcsid );
  ATTATCHSTATUSPTR (stat);

  /* set a few defaults first */
  uvar_ephemYear = LALCalloc (1, strlen(EPHEM_YEARS)+1);
  strcpy (uvar_ephemYear, EPHEM_YEARS);

  uvar_doWindowing = FALSE;
  uvar_binaryoutput = FALSE;
  uvar_nomagic = FALSE;
  uvar_help = FALSE;
  uvar_Tsft = 1800;
  uvar_f1dot = 0.0;
  uvar_f2dot = 0.0;
  uvar_f3dot = 0.0;

  uvar_noiseSigma = 0;
  uvar_noiseSFTs = NULL;

  /* now register all our user-variable */

  LALregSTRINGUserVar(stat, outSFTbname,'n', UVAR_OPTIONAL, "Path and basefilename of output SFT files");
  LALregSTRINGUserVar(stat, outTDDFile,	't', UVAR_OPTIONAL, "Filename for output time-series");
  LALregSTRINGUserVar(stat, detector,  	'I', UVAR_REQUIRED, "Detector: LHO, LLO, VIRGO, GEO, TAMA, CIT, ROME");
  LALregREALUserVar(stat,   startTime,	'G', UVAR_OPTIONAL, "Detector GPS time to start data");
  LALregREALUserVar(stat,   refTime, 	'S', UVAR_OPTIONAL, "Reference time tRef (in SSB) at which pulsar is defined");
  LALregSTRINGUserVar(stat, ephemDir,	'E', UVAR_REQUIRED, "Directory path for ephemeris files");
  LALregSTRINGUserVar(stat, ephemYear, 	'y', UVAR_OPTIONAL, "Year (or range of years) of ephemeris files to be used");
  LALregSTRINGUserVar(stat,timestampsFile,0, UVAR_OPTIONAL, "Timestamps file");
  LALregREALUserVar(stat,   Tsft, 	 0 , UVAR_OPTIONAL, "SFT time baseline Tsft");
  LALregINTUserVar(stat,    nTsft,	 0 , UVAR_REQUIRED, "Number of SFTs nTsft");
  LALregREALUserVar(stat,   fmin,	 0 , UVAR_REQUIRED, "lowest frequency in output SFT");
  LALregREALUserVar(stat,   Band,	 0 , UVAR_REQUIRED, "bandwidth of output SFT");
  LALregREALUserVar(stat,   longitude,	 0 , UVAR_REQUIRED, "Right ascension [radians] alpha of pulsar");
  LALregREALUserVar(stat,   latitude, 	 0 , UVAR_REQUIRED, "Declination [radians] delta of pulsar");
  LALregREALUserVar(stat,   aPlus,	 0 , UVAR_REQUIRED, "Plus polarization amplitude aPlus");
  LALregREALUserVar(stat,   aCross, 	 0 , UVAR_REQUIRED, "Cross polarization amplitude aCross");
  LALregREALUserVar(stat,   psi,  	 0 , UVAR_REQUIRED, "Polarization angle psi");
  LALregREALUserVar(stat,   phi0,	 0 , UVAR_REQUIRED, "Initial phase phi");
  LALregREALUserVar(stat,   f0,  	 0 , UVAR_REQUIRED, "Gravitational wave-frequency f0 at tRef");
  LALregREALUserVar(stat,   f1dot,  	 0 , UVAR_OPTIONAL, "First spindown parameter f'");
  LALregREALUserVar(stat,   f2dot,  	 0 , UVAR_OPTIONAL, "Second spindown parameter f''");
  LALregREALUserVar(stat,   f3dot,  	 0 , UVAR_OPTIONAL, "Third spindown parameter f'''");
  LALregREALUserVar(stat,   noiseSigma,	 0 , UVAR_OPTIONAL, "Gaussian noise variance sigma");
  LALregSTRINGUserVar(stat, noiseSFTs,	'D', UVAR_OPTIONAL, "Glob-like pattern specifying noise-SFTs to be added to signal");  

  /* the orbital parameters */
  LALregREALUserVar(stat,   orbitSemiMajorAxis, 0, UVAR_OPTIONAL, "Projected orbital semi-major axis in seconds (a/c)");
  LALregREALUserVar(stat,   orbitEccentricity,  0, UVAR_OPTIONAL, "Orbital eccentricity");
  LALregINTUserVar(stat,    orbitTperiSSBsec,   0, UVAR_OPTIONAL, "'observed' (SSB) time of periapsis passage. Seconds.");
  LALregINTUserVar(stat,    orbitTperiSSBns,    0, UVAR_OPTIONAL, "'observed' (SSB) time of periapsis passage. Nanoseconds.");
  LALregREALUserVar(stat,   orbitPeriod,        0, UVAR_OPTIONAL, "Orbital period (seconds)");
  LALregREALUserVar(stat,   orbitArgPeriapse,   0, UVAR_OPTIONAL, "Argument of periapsis (radians)");                            

  LALregBOOLUserVar(stat,   help,	'h', UVAR_HELP    , "Print this help/usage message");
  
  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* InitUserVars() */


/*----------------------------------------------------------------------
 * This routine frees up all the memory 
 *----------------------------------------------------------------------*/
void FreeMem (LALStatus* stat, ConfigVars_t *cfg)
{

  INITSTATUS( stat, "FreeMem", rcsid );
  ATTATCHSTATUSPTR (stat);

  
  /* Free config-Variables and userInput stuff */
  TRY (LALDestroyUserVars (stat->statusPtr), stat);

  /* free timestamps if any */
  if (cfg->timestamps){
    TRY (LALDestroyTimestampVector (stat->statusPtr, &(cfg->timestamps)), stat);
  }

  /* free spindown-vector (REAL8) */
  if (cfg->spindown) {
    TRY (LALDDestroyVector (stat->statusPtr, &(cfg->spindown)), stat);
  }

  /* free noise-SFTs */
  if (cfg->noiseSFTs) {
    TRY (LALDestroySFTVector (stat->statusPtr, &(cfg->noiseSFTs)), stat);
  }

  /* Clean up earth/sun Ephemeris tables */
  LALFree(cfg->edat.ephemE);
  LALFree(cfg->edat.ephemS);

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* FreeMem() */


/*reads timestamps file and fills-in timestamps vector*/
void
ReadTimestamps (LALStatus* stat, LIGOTimeGPSVector **timestamps, const CHAR *fname)
{  
  FILE *fp;
  UINT4 i;
  INT4 secs, ns;

  INITSTATUS( stat, "ReadTimestamps", rcsid );
  ATTATCHSTATUSPTR (stat);

  ASSERT (fname, stat, MAKEFAKEDATAC_EBAD, MAKEFAKEDATAC_MSGEBAD );
  ASSERT (timestamps, stat, MAKEFAKEDATAC_EBAD,  MAKEFAKEDATAC_MSGEBAD);
  ASSERT (*timestamps == NULL, stat, MAKEFAKEDATAC_EBAD,  MAKEFAKEDATAC_MSGEBAD);

  if ( (fp = fopen( fname, "r")) == NULL) {
    LALPrintError("\nUnable to open timestampsname file %s\n\n", fname);
    ABORT (stat, MAKEFAKEDATAC_EFILE, MAKEFAKEDATAC_MSGEFILE);
  }

  TRY ( LALCreateTimestampVector (stat->statusPtr, timestamps, uvar_nTsft), stat);    

  for ( i=0; i < (UINT4)uvar_nTsft; i++)
    {
      if (fscanf ( fp, "%d  %d\n", &secs, &ns ) != 2)
	{
	  LALDestroyTimestampVector (stat->statusPtr, timestamps);
	  LALPrintError("\nUnable to read datum from line # %d from file %s\n\n", i+1, fname);
	  ABORT (stat, MAKEFAKEDATAC_EFILE, MAKEFAKEDATAC_MSGEFILE);
	} /* if read-error */
      (*timestamps)->data[i].gpsSeconds = secs;
      (*timestamps)->data[i].gpsNanoSeconds = ns;
    }  /* for i < nTsft */
  
  fclose(fp);

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* ReadTimestamps() */

/*----------------------------------------------------------------------
 * generate Gaussian noise with variance sigma, add it to inSeries
 * returns outSeries
 *
 * NOTE: inSeries is allowed to be identical to outSeries!
 *
 *----------------------------------------------------------------------*/
void
AddGaussianNoise (LALStatus* status, REAL4TimeSeries *outSeries, REAL4TimeSeries *inSeries, REAL4 sigma)
{

  REAL4Vector    *v1 = NULL;
  RandomParams   *randpar = NULL;
  UINT4          numPoints, i;
  INT4 seed;
  FILE *devrandom;
  REAL4Vector *bak;

  INITSTATUS( status, "AddGaussianNoise", rcsid );
  ATTATCHSTATUSPTR (status);
  

  ASSERT ( outSeries, status, MAKEFAKEDATAC_EBAD, MAKEFAKEDATAC_MSGEBAD);
  ASSERT ( inSeries, status, MAKEFAKEDATAC_EBAD, MAKEFAKEDATAC_MSGEBAD);
  ASSERT ( outSeries->data->length == inSeries->data->length, status, MAKEFAKEDATAC_EBAD, MAKEFAKEDATAC_MSGEBAD);
  
  numPoints = inSeries->data->length;

  TRY (LALCreateVector (status->statusPtr, &v1, numPoints), status);
  
  /*
   * Modified so as to not create random number parameters with seed
   * drawn from clock.  Seconds don't change fast enough and sft's
   * look alike.  We open /dev/urandom and read a 4 byte integer from
   * it and use that as our seed.  Note: /dev/random is slow after the
   * first, few accesses.
   */
  
  if ( (devrandom = fopen("/dev/urandom","r")) == NULL )
    {
      LALPrintError ("\nCould not open '/dev/urandom'\n\n");
      ABORT (status, MAKEFAKEDATAC_EFILE, MAKEFAKEDATAC_MSGEFILE);
    }

  if ( fread( (void*)&seed, sizeof(INT4), 1, devrandom) != 1)
    {
      LALPrintError ("\nCould not read from '/dev/urandom'\n\n");
      fclose(devrandom);
      ABORT (status, MAKEFAKEDATAC_EFILE, MAKEFAKEDATAC_MSGEFILE);
    }

  fclose(devrandom);
  
  TRY (LALCreateRandomParams (status->statusPtr, &randpar, seed), status);

  TRY (LALNormalDeviates (status->statusPtr, v1, randpar), status);

  for (i = 0; i < numPoints; i++)
    outSeries->data->data[i] = inSeries->data->data[i] + sigma * v1->data[i];

  /* copy the rest of the time-series structure */
  bak = outSeries->data;
  outSeries = inSeries;		/* copy all struct-entries */
  outSeries->data = bak;	/* restore data-pointer */


  /* destroy randpar*/
  TRY (LALDestroyRandomParams (status->statusPtr, &randpar), status);
  
  /*   destroy v1 */
  TRY (LALDestroyVector (status->statusPtr, &v1), status);


  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* AddGaussianNoise() */


/*----------------------------------------------------------------------
 * 
 * for signals from NS in binary orbits: "Translate" our input parameters 
 * into those required by LALGenerateSpinOrbitCW() (see LAL doc for equations)
 * 
 *----------------------------------------------------------------------*/
void
GetOrbitalParams (LALStatus *stat, BinaryOrbitParams *orbit)
{
  REAL8 OneMEcc;
  REAL8 OnePEcc;
  REAL8 correction;
  LIGOTimeGPS TperiTrue;
  LIGOTimeGPS TperiSSB;
  
  INITSTATUS( stat, "GetOrbitalParams", rcsid );
  ATTATCHSTATUSPTR (stat);
  
  OneMEcc = 1.0 - uvar_orbitEccentricity;
  OnePEcc = 1.0 + uvar_orbitEccentricity;

  /* we need to convert the observed time of periapse passage to the "true" time of periapse passage */
  /* this correction is due to the light travel time from binary barycenter to the source at periapse */
  /* it is what Teviets codes require */
  TperiSSB.gpsSeconds = uvar_orbitTperiSSBsec;
  TperiSSB.gpsNanoSeconds = uvar_orbitTperiSSBns;
  correction = uvar_orbitSemiMajorAxis * OneMEcc * sin(uvar_orbitArgPeriapse);

  TRY (LALAddFloatToGPS (stat->statusPtr, &TperiTrue, &TperiSSB, -correction), stat);

  orbit->orbitEpoch = TperiTrue;
  orbit->omega = uvar_orbitArgPeriapse;
  orbit->rPeriNorm = uvar_orbitSemiMajorAxis * OneMEcc;
  orbit->oneMinusEcc = OneMEcc;
  orbit->angularSpeed = (LAL_TWOPI/uvar_orbitPeriod) * sqrt(OnePEcc/(OneMEcc*OneMEcc*OneMEcc));

  DETATCHSTATUSPTR (stat);
  RETURN(stat);

} /* GetOrbitalParams() */

