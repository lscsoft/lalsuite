 /*-----------------------------------------------------------------------
 *
 * File Name: makefakedata_v4.c
 *
 * Authors: Prix, R., Papa, M.A., 
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

#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/GeneratePulsarSignal.h>


RCSID ("$Id");

/* Error codes and messages */
/* <lalErrTable file="MAKEFAKEDATACErrorTable"> */
#define MAKEFAKEDATAC_ENORM 	0
#define MAKEFAKEDATAC_ESUB  	1
#define MAKEFAKEDATAC_EARG  	2
#define MAKEFAKEDATAC_EBAD  	3
#define MAKEFAKEDATAC_EFILE 	4

#define MAKEFAKEDATAC_MSGENORM "Normal exit"
#define MAKEFAKEDATAC_MSGESUB  "Subroutine failed"
#define MAKEFAKEDATAC_MSGEARG  "Error parsing arguments"
#define MAKEFAKEDATAC_MSGEBAD  "Bad argument values"
#define MAKEFAKEDATAC_MSGEFILE "File IO error"

/* </lalErrTable> */


/***************************************************/

/* Locations of the earth and sun ephemeris data */
#define EARTHDATA "earth00-04.dat"
#define SUNDATA "sun00-04.dat"

CHAR timestampsname[128];

/* Prototypes for the functions defined in this file */
void freemem(LALStatus* stat);
int parseR4(FILE *fp, char* vname, REAL4 *data, const CHAR *fname);
int parseR8(FILE *fp, char* vname, REAL8 *data, const CHAR *fname);
int parseI4(FILE *fp, char* vname, INT4 *data, const CHAR *fname);


/* make it a bit easier for us to register all the user-variables in a constistent way */
#define regUserVar(name,type,option,help) LALRegisterUserVar(stat, #name, type, option, help, &(uvar_ ## name)) 

void initUserVars (LALStatus *stat);
void read_timestamps (LALStatus* status, LIGOTimeGPSVector **timestamps);


/*----------------------------------------------------------------------*/
/* User variables */
CHAR *uvar_inDataFilename;
CHAR *uvar_freqbasefilename;
CHAR *uvar_timebasefilename;
CHAR *uvar_timestampsname;
CHAR *uvar_detector;
REAL8 uvar_startTime;
REAL8 uvar_refTime;
CHAR *uvar_ephemDir;
CHAR *uvar_noiseDir;
BOOLEAN uvar_doWindowing;
BOOLEAN uvar_binaryoutput;
BOOLEAN uvar_nomagic;
INT4 uvar_debug;
BOOLEAN uvar_help;
REAL8 uvar_Tsft;
INT4 uvar_nTsft;
REAL8 uvar_fmin;
REAL8 uvar_Band;
REAL8 uvar_sigma;
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
/*----------------------------------------------------------------------*/

static const PulsarSignalParams empty_params;
static const SFTParams empty_sftParams;
static const LALStatus empty_status;
static const EphemerisData empty_edat;
/*----------------------------------------------------------------------
 * main function 
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[]) 
{
  LALStatus status = empty_status;	/* initialize status */
  PulsarSignalParams params = empty_params;
  SFTParams sftParams = empty_sftParams;;
  SFTVector *SFTs = NULL;
  REAL4TimeSeries *Tseries = NULL;

  EphemerisData edat = empty_edat;
  LALLeapSecFormatAndAcc leapParams = {  LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT };
  INT4 leapSecs;
  LALDetector Detector;
  LIGOTimeGPS startTime = {0, 0};
  LIGOTimeGPS refTime = {0, 0};
  LIGOTimeGPSVector *timestamps = NULL;
  REAL8 duration;
  UINT4 msp = 0;	/* number of spindown-parameters */
  REAL8Vector *spindown = NULL;
  CHAR earthdata[] = EARTHDATA;
  CHAR sundata[]   = SUNDATA;
  SFTVector *noiseSFTs = NULL;

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;	/* exit with returned status-code on error */

  /* ------------------------------ 
   * read user-input 
   *------------------------------*/

  LAL_CALL (initUserVars (&status), &status);	  /* register all user-variables */

  LAL_CALL (LALUserVarReadAllInput (&status, argc,argv), &status);  /* read cmdline & cfgfile  */	

  if (uvar_help) 
    {
      CHAR *helpstring = NULL;
      LAL_CALL (LALUserVarHelpString (&status, &helpstring), &status);
      printf ("\n%s\n", helpstring);
    }
  
  /* ------------------------------
   * do some pre-processing of the user-data 
   * in preparation for the call to LALGeneratePulsarSignal() 
   *------------------------------*/

  /* prepare vector of spindown parameters */
  if (uvar_f3dot != 0) 		msp = 3;	/* counter number of spindowns */
  else if (uvar_f2dot != 0)	msp = 2;
  else if (uvar_f1dot != 0)	msp = 1;
  else 				msp = 0;
  if (msp) {
    LAL_CALL (LALDCreateVector (&status, &spindown, msp), &status);
  }

  /* prepare detector */
  if (uvar_detector == NULL) {
    LALPrintError ("No detector specified!\n");
    return (-1);
  }
  if (!strcmp(uvar_detector,"LHO"))   Detector = lalCachedDetectors[LALDetectorIndexLHODIFF];
  else if (!strcmp(uvar_detector,"LLO"))   Detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  else if (!strcmp(uvar_detector,"VIRGO")) Detector = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  else if (!strcmp(uvar_detector,"GEO"))   Detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  else if (!strcmp(uvar_detector,"TAMA"))  Detector = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
  else if (!strcmp(uvar_detector,"CIT"))   Detector = lalCachedDetectors[LALDetectorIndexCIT40DIFF];
  else {
    LALPrintError ("Unknown detector specified: `%s\n`", uvar_detector);
    return (-1);
  }

  /* get leap-seconds since start of GPS-time */
  LAL_CALL ( LALFloatToGPS (&status, &startTime, &uvar_startTime), &status);

  LAL_CALL ( LALLeapSecs (&status, &leapSecs,  &startTime, &leapParams), &status);

  /* Prepare quantities for barycentering */
  edat.ephiles.earthEphemeris = earthdata;
  edat.ephiles.sunEphemeris   = sundata;
  edat.leap = (INT2) leapSecs;
  LAL_CALL( LALInitBarycenter (&status, &edat), &status);   /* Read in ephemerides */  

  /* read timestamps and set signal-duration */
  timestamps = NULL;
  if (uvar_timestampsname) 
    {
      LIGOTimeGPS t1, t0;
      LAL_CALL (read_timestamps (&status, &timestamps), &status);
      t1 = timestamps->data[uvar_nTsft-1];
      t0 = timestamps->data[0];
      LAL_CALL (LALDeltaFloatGPS(&status, &duration, &t1, &t0), &status);
      duration += uvar_Tsft;
    } 
  else
    duration = uvar_nTsft * uvar_Tsft;

  /* read noise-sft's if given */
  if (uvar_noiseDir)
    {
      CHAR *fglob;
      if( (fglob = LALCalloc (1, strlen(uvar_noiseDir) + 10)) == NULL) {
	LALPrintError ("Out of memory, .. arghhh\n");
	return (-1);
      }
      strcpy (fglob, uvar_noiseDir);
      strcat (fglob, "/SFT");		/* follow convention of makefakedata_v2 */
      LAL_CALL ( LALReadSFTfiles (&status, &noiseSFTs, uvar_fmin, uvar_fmin + uvar_Band, fglob), &status);
      LALFree (fglob);

    } /* if uvar_noisedir */

  /* ------------------------------ 
   * fill the PulsarSignalParams struct 
   *------------------------------*/
  LAL_CALL ( LALFloatToGPS (&status, &refTime, &uvar_refTime), &status);

  /* pulsar params */
  params.pulsar.tRef = refTime;
  params.pulsar.position.longitude = uvar_longitude;
  params.pulsar.position.latitude  = uvar_latitude;
  params.pulsar.position.system    = COORDINATESYSTEM_EQUATORIAL;
  params.pulsar.psi 		   = uvar_psi;
  params.pulsar.aPlus		   = uvar_aPlus;
  params.pulsar.aCross		   = uvar_aCross;
  params.pulsar.phi0		   = uvar_phi0;
  params.pulsar.f0		   = uvar_f0;

  /* orbital params (currently null) */
  params.orbit = NULL;

  /* detector params */
  params.transfer = NULL;	/* detector transfer function (NULL if not used) */	
  params.site = &Detector;	
  params.ephemerides = &edat;

  /* characterize the output time-series */
  params.startTimeGPS   = startTime;
  params.duration     	= (UINT4) ceil(duration); /* length of time-series in seconds */
  params.samplingRate 	= 2.0 * uvar_Band;	/* sampling rate of time-series (= 2 * frequency-Band) */
  params.fHeterodyne  	= uvar_fmin;		/* heterodyning frequency for output time-series */

  /* generate the heterodyned time-series */
  LAL_CALL (LALGeneratePulsarSignal (&status, &Tseries, &params), &status );

  if (lalDebugLevel >= 3) {
    LAL_CALL (PrintR4TimeSeries (&status, Tseries, "debug_tseries.agr"), &status);
  }

  /*----------------------------------------
   * last step: turn it timeseries into SFTs
   *----------------------------------------*/
  sftParams.Tsft = uvar_Tsft;
  sftParams.timestamps = timestamps;
  sftParams.noiseSFTs = noiseSFTs;

  LAL_CALL ( LALSignalToSFTs (&status, &SFTs, Tseries, &sftParams), &status);


  /* write SFTs to disk */


  /* free memory */

  LALCheckMemoryLeaks(); 

  return 0;
} /* main */


/* This routine frees up all the memory */
void freemem(LALStatus* stat)
{

  INITSTATUS( stat, "freemem", rcsid );
  ATTATCHSTATUSPTR (stat);

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* freemem() */


/*reads timestamps file and fills-in timestamps vector*/
void
read_timestamps (LALStatus* stat, LIGOTimeGPSVector **timestamps)
{  
  FILE *fp;
  UINT4 i;
  INT4 secs, ns;

  INITSTATUS( stat, "read_timestamps", rcsid );
  ATTATCHSTATUSPTR (stat);

  ASSERT (uvar_timestampsname, stat, MAKEFAKEDATAC_EBAD, MAKEFAKEDATAC_MSGEBAD );

  if ( (fp = fopen( uvar_timestampsname, "r")) == NULL) {
    LALPrintError("Unable to open timestampsname file %s\n", uvar_timestampsname);
    ABORT (stat, MAKEFAKEDATAC_EFILE, MAKEFAKEDATAC_MSGEFILE);
  }

  TRY ( LALCreateTimestampVector (stat->statusPtr, timestamps, uvar_nTsft), stat);    

  for ( i=0; i < (UINT4)uvar_nTsft; i++)
    {
      if (fscanf ( fp, "%d  %d\n", &secs, &ns ) != 2)
	{
	  LALDestroyTimestampVector (stat->statusPtr, timestamps);
	  LALPrintError("Unable to read datum from line # %d from file %s\n", i+1, uvar_timestampsname);
	  ABORT (stat, MAKEFAKEDATAC_EFILE, MAKEFAKEDATAC_MSGEFILE);
	} /* if read-error */
    }  /* for i < nTsft */
  
  fclose(fp);

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* read_timestamps() */


/* for backwards compatibility */
int parseR4(FILE *fp, char* vname, REAL4 *data, const CHAR *fname)
{
  char junk[1024], junk2[1024];
  char command[1024];
  int r;
  
  memset(junk, '\0', 1024);
  memset(junk2,'\0', 1024);
  
  r=fscanf(fp, "%f%[\t ]%[^\012]", data, junk, junk2);
  if (r!=3)  {
    LALPrintError ("Unable to assign %s from file: %s\n"
	  "The entry must be of the form:\n"
	  "NUMBER TEXT\n"
	  "with white space in between. TEXT is NOT optional!\n",
	  vname, fname);
    sprintf(command, "cat %s 1>&2\n", fname);
    system(command);
    return 1;
  }
      return 0;
} /* parseR4() */


int parseR8(FILE *fp, char* vname, REAL8 *data, const CHAR *fname)
{
  char junk[1024], junk2[1024];
  char command[1024];
  int r;
  
  memset(junk, '\0', 1024);
  memset(junk2,'\0', 1024);
  
  r=fscanf(fp, "%lf%[\t ]%[^\n]", data, junk, junk2);
  if (r!=3)  {
    LALPrintError ("Unable to assign %s from file: %s\n"
	  "The entry must be of the form:\n"
	  "NUMBER TEXT\n"
	  "with white space in between. TEXT is NOT optional!\n",
	  vname, fname);
    sprintf(command, "cat %s 1>&2\n", fname);
    system(command);
    return 1;
  }
      return 0;
} /* parseR8() */

int parseI4(FILE *fp, char* vname, INT4 *data, const CHAR *fname)
{
  char junk[1024], junk2[1024];
  char command[1024];
  int r;
  
  memset(junk, '\0', 1024);
  memset(junk2,'\0', 1024);
  
  r=fscanf(fp, "%d%[\t ]%[^\n]", data, junk, junk2);
  if (r!=3)  {
    LALPrintError ("Unable to assign %s from file: %s\n"
	  "The entry must be of the form:\n"
	  "NUMBER TEXT\n"
	  "with white space in between. TEXT is NOT optional!\n",
	  vname, fname);
    sprintf(command, "cat %s 1>&2\n", fname);
    system(command);
    return 1;
  }
      return 0;
} /* parseI4() */

/*----------------------------------------------------------------------*/
/* register all our "user-variables" */
void
initUserVars (LALStatus *stat)
{
  INITSTATUS( stat, "initUserVars", rcsid );
	      
  regUserVar(inDataFilename,	UVAR_STRING, 'i', "Name of input parameter file");
  regUserVar(freqbasefilename, 	UVAR_STRING, 'n', "Basefilename of output SFT files");
  regUserVar(timebasefilename, 	UVAR_STRING, 't', "Basefilename of output STRAIN files");
  regUserVar(detector, 		UVAR_STRING, 'I', "Detector: LHO, LLO, VIRGO, GEO, TAMA, CIT, ROME");
  regUserVar(startTime, 	UVAR_REAL8,  'G', "Detector GPS time to start data");
  regUserVar(refTime, 		UVAR_REAL8,  'S', "Reference time tRef at which pulsar is defined");
  regUserVar(ephemDir, 		UVAR_STRING, 'E', "Directory path for ephemeris files");
  regUserVar(noiseDir, 		UVAR_STRING, 'D', "Directory with noise-SFTs");
  regUserVar(doWindowing, 	UVAR_BOOL,   'w', "Window data in time domain before doing FFT");
  regUserVar(binaryoutput, 	UVAR_BOOL,   'b', "Output time-domain data in IEEE754 binary format");
  regUserVar(nomagic, 		UVAR_BOOL,   'm', "DON'T output 1234.5 before time-domain binary samples");
  regUserVar(debug, 		UVAR_INT4,   'v', "set debug-level");
  regUserVar(help, 		UVAR_BOOL,   'h', "Print this help/usage message");
  regUserVar(Tsft, 		UVAR_REAL8,  'T', "SFT time baseline Tsft");
  regUserVar(nTsft, 		UVAR_INT4,   'N', "Number of SFTs nTsft");
  regUserVar(fmin, 		UVAR_REAL8,  'F', "minimum frequency fmin of output SFT");
  regUserVar(Band, 		UVAR_REAL8,  'B', "bandwidth of output SFT");
  regUserVar(sigma, 		UVAR_REAL8,  's', "noise variance sigma");
  regUserVar(aPlus, 		UVAR_REAL8,  'a', "Plus polarization amplitude aPlus");
  regUserVar(aCross, 		UVAR_REAL8,  'x', "Cross polarization amplitude aCross");
  regUserVar(psi, 		UVAR_REAL8,  'P', "Polarization angle psi");
  regUserVar(phi0, 		UVAR_REAL8,  'p', "Initial phase phi");
  regUserVar(f0, 		UVAR_REAL8,  'f', "Pulsar frequency f0 at tRef");
  regUserVar(latitude, 		UVAR_REAL8,  'a', "Declination [radians] delta of pulsar");
  regUserVar(longitude, 	UVAR_REAL8,  'd', "Right ascension [radians] alpha of pulsar");
  regUserVar(f1dot, 		UVAR_REAL8,  '1', "First spindown parameter f'");
  regUserVar(f2dot, 		UVAR_REAL8,  '2', "Second spindown parameter f''");
  regUserVar(f3dot, 		UVAR_REAL8,  '3', "Third spindown parameter f'''");

  RETURN (stat);

} /* initUserVars() */
