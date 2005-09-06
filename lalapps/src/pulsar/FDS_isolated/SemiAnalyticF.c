/*********************************************************************************/
/*                 Semi-Analytic calculation of the F-statistic                  */
/*                                                                               */
/*			           X. Siemens                                    */
/*                                                                               */
/*                             UWM - February 2003                               */
/*********************************************************************************/

#include <errno.h>

#include <lal/AVFactories.h>
#include <lal/UserInput.h>
#include <lal/LALStdio.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LALComputeAM.h>


#include <lalapps.h>


RCSID( "$Id$");

/*---------- error-codes ---------- */
#define SEMIANALYTIC_ENORM 		0
#define SEMIANALYTIC_ESUB  		1
#define SEMIANALYTIC_EINPUT  		2
#define SEMIANALYTIC_EBAD  		3
#define SEMIANALYTIC_EFILE 		4
#define SEMIANALYTIC_ENOARG 		5
#define SEMIANALYTIC_EMEM 		6
#define SEMIANALYTIC_EREADFILE 		8

#define SEMIANALYTIC_MSGENORM 		"Normal exit"
#define SEMIANALYTIC_MSGESUB  		"Subroutine failed"
#define SEMIANALYTIC_MSGEINPUT 		"Invalid input"
#define SEMIANALYTIC_MSGEBAD  		"Bad argument values"
#define SEMIANALYTIC_MSGEFILE 		"File IO error"
#define SEMIANALYTIC_MSGENOARG 		"Missing argument"
#define SEMIANALYTIC_MSGEMEM 		"Out of memory..."
#define SEMIANALYTIC_MSGEREADFILE 	"Error reading in file"

/*---------- defines ---------- */
#define TRUE (1==1)
#define FALSE (1==0)

#define SQ(x) ((x)*(x))

/*---------- local types ---------- */
struct CommandLineArgsTag {
  REAL8 skyalpha;
  REAL8 skydelta;
  REAL8 tsft;
  UINT4 nTsft;
  CHAR *detector;
  CHAR *timestamps;
  UINT4 gpsStart;
  CHAR *efiles;
  REAL8 phi;
  REAL8 psi;
  REAL8 h0;
  REAL8 cosiota;
  REAL8 sqrtSh;
  REAL8 duration;
  CHAR *ephemYear;
  REAL8 aPlus;
  REAL8 aCross;
  BOOLEAN help;
} CommandLineArgs;

/*---------- global variables ---------- */
LIGOTimeGPS *timestamps = NULL;       /* Time stamps from SFT data */
AMCoeffs amc;

extern int vrbflg;

/*---------- local prototypes ---------- */
void InitUserVars (LALStatus *status, struct CommandLineArgsTag *CLA);
void ReadUserInput (LALStatus *, struct CommandLineArgsTag *CLA, int argc,char *argv[]);
void Freemem( LALStatus *);
void CreateNautilusDetector (LALStatus *, LALDetector *detector);
void Initialize (LALStatus *status, struct CommandLineArgsTag *CLA);
void ComputeF(LALStatus *, struct CommandLineArgsTag CLA);

int ReadTimeStamps(struct CommandLineArgsTag CLA);
int MakeTimeStamps(struct CommandLineArgsTag CLA);
void CheckUserInput (LALStatus *,  struct CommandLineArgsTag *CLA );

/*---------- function definitions ---------- */
int main(int argc,char *argv[]) 
{
  LALStatus status = blank_status;	/* initialize status */

  lalDebugLevel = 0;	/* default value */
  vrbflg = 1;		/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;	/* exit with returned status-code on error */

  /*----------  read user-input and set up shop ----------*/
  LAL_CALL (LALGetDebugLevel(&status, argc, argv, 'v'), &status);

  /* register all user-variables */
  LAL_CALL (InitUserVars(&status, &CommandLineArgs), &status);	  

  /* read cmdline & cfgfile  */	
  LAL_CALL (LALUserVarReadAllInput(&status, argc, argv), &status);
  if (CommandLineArgs.help)	/* help was called, do nothing here.. */
    return (0);
  LAL_CALL ( CheckUserInput (&status, &CommandLineArgs), &status);

  LAL_CALL ( Initialize (&status, &CommandLineArgs), &status);

  /*---------- central function: compute F-statistic ---------- */
  LAL_CALL ( ComputeF(&status, CommandLineArgs), &status); 

  
  /* Free remaining memory */
  LAL_CALL ( LALSDestroyVector(&status, &(amc.a) ), &status);
  LAL_CALL ( LALSDestroyVector(&status, &(amc.b) ), &status);
  LAL_CALL ( LALDestroyUserVars(&status), &status);

  LALCheckMemoryLeaks();

  return 0;

} /* main() */

/*******************************************************************************/

int MakeTimeStamps(struct CommandLineArgsTag CLA) 
{
  UINT4 i;
 
  /* allocate memory for timestamps */
  timestamps=(LIGOTimeGPS *)LALMalloc(CLA.nTsft*sizeof(LIGOTimeGPS)); 
      
  /* generate timetamps */
  for (i=0;i<CLA.nTsft;i++){
    timestamps[i].gpsSeconds=CLA.gpsStart+(int)(i*CLA.tsft);
    timestamps[i].gpsNanoSeconds=0;
  } 
  
  return 0;
  
}
/*******************************************************************************/

int ReadTimeStamps(struct CommandLineArgsTag CLA) 
{
  FILE *fp;
  UINT4 i;
  int r;
 
 
  /*   %strcpy(filename,inDataFilename); */
  fp=fopen(CLA.timestamps,"r");
  if (fp==NULL) {
    fprintf(stderr,"Unable to find file %s\n",CLA.timestamps);
    return 1;
  }
  timestamps=(LIGOTimeGPS *)LALMalloc(CLA.nTsft*sizeof(LIGOTimeGPS)); 
      
  
  for (i=0;i<CLA.nTsft;i++){
    r=fscanf(fp,"%d  %d\n", &timestamps[i].gpsSeconds, &timestamps[i].gpsNanoSeconds);
    if ( r !=2 ) {
      fprintf(stderr,"Unable to read datum # %d\n",i);
      fprintf(stderr,"from file %s\n",CLA.timestamps);
      return 1; 
    } 
  } 
  
  fclose(fp);
  return 0;
  
}
/*******************************************************************************/

void 
ComputeF( LALStatus *status, struct CommandLineArgsTag CLA)
{

  REAL8 A,B,C,D,alpha,beta,delta,kappa,A1,A2,A3,A4,h0,cosi, To,Sh,F;
  REAL8 aPlus, aCross;
  REAL8 twopsi, twophi;
  REAL8 lambda;

  INITSTATUS (status, "ComputeF", rcsid );
  ATTATCHSTATUSPTR ( status);
  
  A = amc.A;
  B = amc.B;
  C = amc.C;
  D = amc.D; 

  twophi = 2.0 * CLA.phi;
  twopsi = 2.0 * CLA.psi;

  h0 = CLA.h0;
  cosi = CLA.cosiota;

  if ( h0 != 0 ) 
    {
      aPlus = h0 * 0.5 * (1.0 + cosi*cosi );
      aCross= h0 * cosi;
    } 
  else   /* alternative way to specify amplitude (compatible with mfd_v4) */
    {
      aPlus = CLA.aPlus;
      aCross = CLA.aCross;
    }
  
  A1 = aPlus * cos(twopsi) * cos(twophi) - aCross * sin(twopsi) * sin(twophi);
  A2 = aPlus * sin(twopsi) * cos(twophi) + aCross * cos(twopsi) * sin(twophi);
  A3 =-aPlus * cos(twopsi) * sin(twophi) - aCross * sin(twopsi) * cos(twophi);
  A4 =-aPlus * sin(twopsi) * sin(twophi) + aCross * cos(twopsi) * cos(twophi);
  
  alpha = 0.5 * A * A1 + 0.5 * C * A2;
  beta  = 0.5 * B * A2 + 0.5 * C * A1;
  delta = 0.5 * A * A3 + 0.5 * C * A4;
  kappa = 0.5 * B * A4 + 0.5 * C * A3;
  
  To = CLA.nTsft * CLA.tsft;
  
  Sh=pow(CLA.sqrtSh,2);
  
  F = (B*pow(alpha,2) + A*pow(beta,2) - 2*C*alpha*beta) 
    + (B*pow(delta,2) + A*pow(kappa,2)- 2*C*delta*kappa);
  F *= To / (D * Sh);

  /* Note: the expectation-value of 2F is 4 + lambda ==> add 2 to Fstat*/
  F += 2.0;

  lambda = SQ(A1) * A + 2.0 * C * A1 * A2 + SQ(A2) * B + SQ(A3) * A + 2.0 * C * A3 * A4 + SQ(A4) * B;
  lambda *= 0.5 * To / Sh;
  
  fprintf( stdout, "%g\n", F );
  /*    fprintf(stdout,"\nSNR = %e\n",sqrt(F)); */


  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* ComputeF() */


/** 
 * register all our "user-variables" 
 */
void
InitUserVars (LALStatus *status, struct CommandLineArgsTag *CLA)
{

  INITSTATUS( status, "InitUserVars", rcsid );
  ATTATCHSTATUSPTR (status);

  /* Initialize default values */
  CLA->skyalpha=0.0;
  CLA->skydelta=0.0;
  CLA->tsft=1800;
  CLA->detector=0;               
  CLA->nTsft=0;            
  CLA->timestamps=NULL;
  CLA->gpsStart=-1;
  CLA->efiles=NULL;
  CLA->phi=0.0;
  CLA->psi=0.0;
  CLA->cosiota=0.0;
  CLA->h0 = 0;
  CLA->sqrtSh=1.0;
  
  /** Default year-span of ephemeris-files to be used */
#define EPHEM_YEARS  "00-04"
  CLA->ephemYear = LALCalloc(1, strlen(EPHEM_YEARS)+1);
  strcpy (CLA->ephemYear, EPHEM_YEARS);
  
#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
  CLA->efiles = LALCalloc(1, strlen(DEFAULT_EPHEMDIR)+1);
  strcpy (CLA->efiles, DEFAULT_EPHEMDIR);
  
  CLA->help = FALSE;
  
  /* ---------- register all our user-variable ---------- */
  TRY (LALRegisterBOOLUserVar(status->statusPtr, "help", 'h', UVAR_HELP, "Print this message",
			      &(CLA->help)), status); 
  TRY( LALRegisterREALUserVar(status->statusPtr, "latitude", 'd', UVAR_REQUIRED, 
			      "Sky position delta (equatorial coordinates) in radians", 
			      &(CLA->skydelta)), status);
  
  TRY( LALRegisterREALUserVar(status->statusPtr, "longitude", 'a', UVAR_REQUIRED, 
			      "Sky position alpha (equatorial coordinates) in radians", 
			      &(CLA->skyalpha)), status);
  
  TRY( LALRegisterREALUserVar(status->statusPtr, "phi0", 'Q', UVAR_OPTIONAL, 
			     "Phi_0: Initial phase in radians", &(CLA->phi)), status);

  TRY( LALRegisterREALUserVar(status->statusPtr, "psi", 'Y', UVAR_OPTIONAL, 
			     "Polarisation in radians", &(CLA->psi)), status);

  TRY( LALRegisterREALUserVar(status->statusPtr, "cosiota", 'i', UVAR_OPTIONAL, 
			      "Cos(iota)", &(CLA->cosiota)), status);
  TRY( LALRegisterREALUserVar(status->statusPtr, "h0", 's', UVAR_OPTIONAL, 
			      "Strain amplitude h_0", &(CLA->h0)), status);
  TRY( LALRegisterREALUserVar(status->statusPtr, "sqrtSh", 'N', UVAR_OPTIONAL, 
			      "Noise floor: one-sided sqrt(Sh) in 1/sqrt(Hz)", &(CLA->sqrtSh)), status);
  
  TRY( LALRegisterSTRINGUserVar(status->statusPtr, "timestampsFile", 'T', UVAR_OPTIONAL, 
				"Name of timestamps file", &(CLA->timestamps)), status);
  
  TRY( LALRegisterINTUserVar(status->statusPtr, "startTime", 'S', UVAR_OPTIONAL, 
			     "GPS start time of continuous observation", &(CLA->gpsStart)), status);
  
  TRY( LALRegisterREALUserVar(status->statusPtr, "Tsft", 't', UVAR_OPTIONAL, 
			      "Length of an SFT in seconds", &(CLA->tsft)), status);
  
  TRY( LALRegisterINTUserVar(status->statusPtr, "nTsft", 'n', UVAR_OPTIONAL, 
			     "Number of SFTs", &(CLA->nTsft)), status);
  
  TRY( LALRegisterSTRINGUserVar(status->statusPtr, "ephemDir", 'E', UVAR_OPTIONAL, 
				"Directory where Ephemeris files are located", 
				&(CLA->efiles)), status);
  
  TRY( LALRegisterSTRINGUserVar(status->statusPtr, "detector", 'D', UVAR_REQUIRED, 
				"Detector: GEO(0),LLO(1),LHO(2),NAUTILUS(3),VIRGO(4),TAMA(5),CIT(6)",
				&(CLA->detector)), status);
  
  /* ----- added for mfd_v4 compatibility ---------- */
  TRY ( LALRegisterREALUserVar(status->statusPtr, "duration", 0, UVAR_OPTIONAL,
			       "Duration of requested signal in seconds", 
			       &(CLA->duration)), status); 
  
  TRY ( LALRegisterSTRINGUserVar(status->statusPtr, "ephemYear", 0, UVAR_OPTIONAL,
				 "Year (or range of years) of ephemeris files to be used",
				 &(CLA->ephemYear)), status);
  
  TRY ( LALRegisterREALUserVar(status->statusPtr, "aPlus", 0, UVAR_OPTIONAL, 
			       "Plus polarization amplitude aPlus", 
			       &(CLA->aPlus)), status);
  TRY ( LALRegisterREALUserVar(status->statusPtr, "aCross", 0, UVAR_OPTIONAL, 
			       "Cross polarization amplitude aCross", 
			       &(CLA->aCross)), status);


  DETATCHSTATUSPTR (status);
  RETURN(status);
} /* InitUserVars() */
  
/** 
 * Handle user-input and check its validity. 
 * Load ephemeris and calculate AM-coefficients (stored globally) 
 */
void
Initialize (LALStatus *status, struct CommandLineArgsTag *CLA)
{
  EphemerisData *edat=NULL;          /* Stores earth/sun ephemeris data for barycentering */
  BarycenterInput baryinput;         /* Stores detector location and other barycentering data */
  EarthState earth;
  AMCoeffsParams *amParams;
  LIGOTimeGPS *midTS=NULL;           /* Time stamps for amplitude modulation coefficients */
  LALLeapSecFormatAndAcc formatAndAcc = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
  LALDetector Detector;              /* Our detector*/
  UINT4 k;

  INITSTATUS (status, "Initialize", rcsid);
  ATTATCHSTATUSPTR (status);


  if ( LALUserVarWasSet (&(CLA->duration) ) )
    CLA->nTsft = (UINT4) (CLA->duration / CLA->tsft + 0.5);	  /* we're cheating here */

  /* read or generate SFT timestamps */
  if ( LALUserVarWasSet(&(CLA->timestamps)) ) {
    if (ReadTimeStamps(*CLA)) {
      ABORT ( status,  SEMIANALYTIC_ESUB,  SEMIANALYTIC_MSGESUB);
    }
  } else {
    if (MakeTimeStamps(*CLA)) {
      ABORT ( status,  SEMIANALYTIC_ESUB,  SEMIANALYTIC_MSGESUB);
    }
  }


  /*---------- initialize detector ---------- */
  if ( !strcmp (CLA->detector, "GEO") || !strcmp (CLA->detector, "0") ) 
    Detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  else if ( !strcmp (CLA->detector, "LLO") || ! strcmp (CLA->detector, "1") ) 
    Detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  else if ( !strcmp (CLA->detector, "LHO") || !strcmp (CLA->detector, "2") )
    Detector = lalCachedDetectors[LALDetectorIndexLHODIFF];
  else if ( !strcmp (CLA->detector, "NAUTILUS") || !strcmp (CLA->detector, "3"))
    {
      TRY (CreateNautilusDetector (status->statusPtr, &(Detector)), status);
    }
  else if ( !strcmp (CLA->detector, "VIRGO") || !strcmp (CLA->detector, "4") )
    Detector = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  else if ( !strcmp (CLA->detector, "TAMA") || !strcmp (CLA->detector, "5") )
    Detector = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
  else if ( !strcmp (CLA->detector, "CIT") || !strcmp (CLA->detector, "6") )
    Detector = lalCachedDetectors[LALDetectorIndexCIT40DIFF];
  else
    {
      LALPrintError ("\nUnknown detector. Currently allowed are 'GEO', 'LLO', 'LHO',"
		     " 'NAUTILUS', 'VIRGO', 'TAMA', 'CIT' or '0'-'6'\n\n");
      ABORT (status, SEMIANALYTIC_EINPUT, SEMIANALYTIC_MSGEINPUT);
    }

  /* ---------- load ephemeris-files ---------- */
#define MAXFILENAME 256
  {
    CHAR filenameE[MAXFILENAME], filenameS[MAXFILENAME];
    INT4 leap;
    
    /* don't use absolute path if none was given, this
     * allows LAL to find the ephemeris in LAL_DATA_PATH */
    if ( LALUserVarWasSet (&(CLA->efiles)) ) 
      {
	LALSnprintf (filenameE, MAXFILENAME, "%s/earth%s.dat", CLA->efiles, CLA->ephemYear );
	LALSnprintf (filenameS, MAXFILENAME, "%s/sun%s.dat", CLA->efiles, CLA->ephemYear );
      }
    else
      {
	LALSnprintf (filenameE, MAXFILENAME, "earth%s.dat", CLA->ephemYear );
	LALSnprintf (filenameS, MAXFILENAME, "sun%s.dat", CLA->ephemYear );
      }
    filenameE[MAXFILENAME-1] = 0;
    filenameS[MAXFILENAME-1] = 0;

    edat = (EphemerisData *)LALMalloc(sizeof(EphemerisData));
    (*edat).ephiles.earthEphemeris = filenameE;     
    (*edat).ephiles.sunEphemeris = filenameS;         

    TRY ( LALLeapSecs(status->statusPtr,&leap,&timestamps[0],&formatAndAcc), status);
    (*edat).leap=leap; 

    /* Reads in ephemeris files */
    TRY( LALInitBarycenter (status->statusPtr, edat), status );

  } /* ephemeris-reading */


  /* ---------- calculate AM-coefficients ---------- */

  /* prepare call to barycentering routing */
  baryinput.site.location[0] = Detector.location[0]/LAL_C_SI;
  baryinput.site.location[1] = Detector.location[1]/LAL_C_SI;
  baryinput.site.location[2] = Detector.location[2]/LAL_C_SI;
  baryinput.alpha = CLA->skyalpha;
  baryinput.delta = CLA->skydelta;
  baryinput.dInv = 0.e0;

  /* amParams structure to compute a(t) and b(t) */

  /* Allocate space for amParams stucture */
  /* Here, amParams->das is the Detector and Source info */
  amParams = (AMCoeffsParams *)LALMalloc(sizeof(AMCoeffsParams));
  amParams->das = (LALDetAndSource *)LALMalloc(sizeof(LALDetAndSource));
  amParams->das->pSource = (LALSource *)LALMalloc(sizeof(LALSource));
  /* Fill up AMCoeffsParams structure */
  amParams->baryinput = &baryinput;
  amParams->earth = &earth; 
  amParams->edat = edat;
  amParams->das->pDetector = &Detector; 
  amParams->das->pSource->equatorialCoords.latitude = CLA->skydelta;
  amParams->das->pSource->equatorialCoords.longitude = CLA->skyalpha;
  amParams->das->pSource->orientation = 0.0;
  amParams->das->pSource->equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  amParams->polAngle = amParams->das->pSource->orientation ; /* These two have to be the same!!!!!!!!!*/
  amParams->leapAcc = formatAndAcc.accuracy;
  
  /* Allocate space for AMCoeffs */
  amc.a = NULL;
  amc.b = NULL;
  TRY ( LALSCreateVector(status->statusPtr, &(amc.a), (UINT4)  CLA->nTsft), status);
  TRY ( LALSCreateVector(status->statusPtr, &(amc.b), (UINT4)  CLA->nTsft), status);
  
  /* Mid point of each SFT */
  midTS = (LIGOTimeGPS *)LALCalloc(CLA->nTsft,sizeof(LIGOTimeGPS));
  for(k=0; k<CLA->nTsft; k++)
    {
      REAL8 teemp=0.0;
      TRY ( LALGPStoFloat(status->statusPtr, &teemp, &(timestamps[k])), status);
      teemp += 0.5*CLA->tsft;
      TRY ( LALFloatToGPS(status->statusPtr, &(midTS[k]), &teemp), status);
    }
  
  TRY ( LALComputeAM(status->statusPtr, &amc, midTS, amParams), status);

  /* Free memory */
  LALFree(timestamps);
  LALFree(midTS);

  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);

  LALFree(amParams->das->pSource);
  LALFree(amParams->das);
  LALFree(amParams);


  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* ParseUserInput() */


/*******************************************************************************/
/** Set up the \em LALDetector struct representing the NAUTILUS detector */
void
CreateNautilusDetector (LALStatus *status, LALDetector *detector)
{
  /*   LALDetector Detector;  */
  LALFrDetector detector_params;
  LALDetectorType bar;
  LALDetector Detector1;

  INITSTATUS (status, "CreateNautilusDetector", rcsid);
  ATTATCHSTATUSPTR (status);

  bar=LALDETECTORTYPE_CYLBAR;
  strcpy(detector_params.name, "NAUTILUS");
  detector_params.vertexLongitudeRadians=12.67*LAL_PI/180.0;
  detector_params.vertexLatitudeRadians=41.82*LAL_PI/180.0;
  detector_params.vertexElevation=300.0;
  detector_params.xArmAltitudeRadians=0.0;
  detector_params.xArmAzimuthRadians=44.0*LAL_PI/180.0;

  TRY (LALCreateDetector(status->statusPtr, &Detector1, &detector_params, bar), status);
  
  *detector=Detector1;

  DETATCHSTATUSPTR (status);
  RETURN (status);
  
} /* CreateNautilusDetector() */

/** 
 * Check validity of user-input
 */
void
CheckUserInput (LALStatus *status,  struct CommandLineArgsTag *CLA )
{
  INITSTATUS (status, "CheckUserInput", rcsid);

  /* set a few abbreviations */
  BOOLEAN have_timestamps= LALUserVarWasSet (&(CLA->timestamps));
  BOOLEAN have_gpsStart = LALUserVarWasSet  (&(CLA->gpsStart));
  BOOLEAN have_duration  = LALUserVarWasSet (&(CLA->duration));
  BOOLEAN have_nTsft     = LALUserVarWasSet (&(CLA->nTsft));
  
  if( have_timestamps && (have_gpsStart||have_duration) )
    {
      fprintf(stderr,"\nBoth start time/duration and timestamps file specified - just need one !!\n");
      fprintf(stderr,"Try ./lalapps_SemiAnalyticF -h \n\n");
      ABORT (status, SEMIANALYTIC_EINPUT, SEMIANALYTIC_MSGEINPUT);
    }   
  
  if( !have_timestamps && !have_gpsStart )
    {
      fprintf(stderr,"\nNeed to specify gpsStart time or a timestamps file !!\n");
      fprintf(stderr,"Try ./lalapps_SemiAnalyticF -h \n\n");
      ABORT (status, SEMIANALYTIC_EINPUT, SEMIANALYTIC_MSGEINPUT);
    }
  
  if ( have_duration && have_nTsft )
    {
      fprintf (stderr, "\nSpecify only one of {duration, nTsft}!\n\n");
      ABORT (status, SEMIANALYTIC_EINPUT, SEMIANALYTIC_MSGEINPUT);
    }

  /* now one can either specify {h0, cosiota} OR {aPlus, aCross} */
  {
    BOOLEAN have_h0 = LALUserVarWasSet (&(CLA->h0));
    BOOLEAN have_cosiota = LALUserVarWasSet (&(CLA->cosiota));
    BOOLEAN have_aPlus = LALUserVarWasSet (&(CLA->aPlus));
    BOOLEAN have_aCross = LALUserVarWasSet (&(CLA->aCross));
    
    if ( (have_h0 || have_cosiota) && (have_aPlus || have_aCross) ) 
      {
	fprintf (stderr, "\nSpecify only one set of {h0/cosiota} or {aPlus/aCross}\n\n");
	ABORT (status, SEMIANALYTIC_EINPUT, SEMIANALYTIC_MSGEINPUT);
      }

    if ( !have_h0 && !(have_aPlus) )
      {
	fprintf (stderr, "\nYou need to specify either h0 or aPlus/aCross\n\n");
	ABORT (status, SEMIANALYTIC_EINPUT, SEMIANALYTIC_MSGEINPUT);
      }

  }

  RETURN (status);

} /* CheckUserInput() */
