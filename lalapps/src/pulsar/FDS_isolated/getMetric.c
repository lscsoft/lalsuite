/**
 * \file getMetric.c
 *
 * \author{Reinhard Prix}
 *
 * Standalone code to calculated the metric in a given parameter-space point
 * using the MetricWrapper.
 *
 * Revision: $Id$
 *           
 *-----------------------------------------------------------------------*/

/* ---------- includes ---------- */
#include <math.h>

#include <lalapps.h>

#include <lal/UserInput.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/LALInitBarycenter.h>
#include <lal/AVFactories.h>

#include "DopplerScan.h"

RCSID ("$Id$");

/* Error codes and messages */
#define GETMETRIC_ENORM 	0
#define GETMETRIC_ESUB  	1
#define GETMETRIC_EARG  	2
#define GETMETRIC_EBAD  	3
#define GETMETRIC_EFILE 	4
#define GETMETRIC_ENOARG 	5
#define GETMETRIC_EINPUT 	6


#define GETMETRIC_MSGENORM 	"Normal exit"
#define GETMETRIC_MSGESUB  	"Subroutine failed"
#define GETMETRIC_MSGEARG  	"Error parsing arguments"
#define GETMETRIC_MSGEBAD  	"Bad argument values"
#define GETMETRIC_MSGEFILE 	"File IO error"
#define GETMETRIC_MSGENOARG 	"Missing argument"
#define GETMETRIC_MSGEINPUT 	"Invalid user input"

/*---------- local defines ---------- */
#define TRUE (1==1)
#define FALSE (1==0)

/* ---------- some local types ---------- */
typedef struct {
  CHAR EphemEarth[512];		/**< filename of earth-ephemeris data */
  CHAR EphemSun[512];		/**< filename of sun-ephemeris data */
  LALDetector Detector;         /**< detector of data to be searched */
  EphemerisData *ephemeris;	/**< ephemeris data (from LALInitBarycenter()) */
  LIGOTimeGPS startTimeGPS;	/**< starttime of observation */
} ConfigVariables;

/*---------- empty structs for initializations ----------*/
static const ConfigVariables empty_ConfigVariables;
static const PtoleMetricIn empty_metricpar;

/* ---------- local prototypes ---------- */
void initUserVars (LALStatus *stat);
void initGeneral (LALStatus *lstat, ConfigVariables *cfg);
void CreateNautilusDetector (LALStatus *lstat, LALDetector *Detector);

/* ---------- User variables ---------- */
BOOLEAN uvar_help;

CHAR* uvar_IFO;			/**< name of detector (LHO, LLO, GEO, VIRGO,..)*/
REAL8 uvar_Alpha;		/**< skyposition Alpha: radians, equatorial coords. */
REAL8 uvar_Delta;		/**< skyposition Delta: radians, equatorial coords. */
REAL8 uvar_Freq;		/**< target-frequency */
REAL8 uvar_f1dot;		/**< target 1. spindown-value df/dt */
CHAR *uvar_ephemDir;		/**< directory of ephemeris-files */
CHAR *uvar_ephemYear;		/**< year-range of ephemeris-file to use */
INT4  uvar_metricType;		/**< Metric function to use: ptole_analytic, ptole-numeric, ... */
REAL8 uvar_startTime;		/**< GPS start time of observation */
REAL8 uvar_duration;		/**< length of observation in seconds */
BOOLEAN uvar_projectMetric;	/**< project out frequency-component of metric */


extern int vrbflg;

/*----------------------------------------------------------------------
 * main function 
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[]) 
{
  LALStatus status = blank_status;
  ConfigVariables config = empty_ConfigVariables;

  lalDebugLevel = 0;  
  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register user-variables */
  LAL_CALL (LALGetDebugLevel (&status, argc, argv, 'v'), &status);
  LAL_CALL (initUserVars (&status), &status);	  

  /* read cmdline & cfgfile  */	
  LAL_CALL (LALUserVarReadAllInput (&status, argc,argv), &status);  

  if (uvar_help) 	/* help requested: we're done */
    exit (0);

  /* some general setup of the code */
  LAL_CALL (initGeneral (&status, &config), &status);


  /* ---------- main-task: get metric for given parameters ---------- */
  {
    REAL8Vector *metric = NULL;
    PtoleMetricIn metricpar = empty_metricpar;

    metricpar.position.system = COORDINATESYSTEM_EQUATORIAL;
    metricpar.position.longitude = uvar_Alpha;
    metricpar.position.latitude = uvar_Delta;

    LAL_CALL ( LALSCreateVector (&status, &(metricpar.spindown), 1), &status);
    metricpar.spindown->data[0] = uvar_f1dot / uvar_Freq; /* f1 = df/dt / f0 !!*/
      
    metricpar.epoch = config.startTimeGPS;
    metricpar.duration = (REAL4) uvar_duration;
    metricpar.maxFreq = uvar_Freq;
      
    metricpar.site = &(config.Detector);
    metricpar.ephemeris = config.ephemeris;	/* needed for ephemeris-metrics */
    metricpar.metricType = uvar_metricType;

    LAL_CALL ( LALMetricWrapper(&status, &metric, &metricpar), &status);
    LAL_CALL ( LALSDestroyVector(&status, &(metricpar.spindown)), &status);

    if (uvar_projectMetric) {
      LAL_CALL ( LALProjectMetric( &status, metric, 0 ), &status);
    }

    printf ("\n %s Metric: (f, alpha, delta, f1) \n", 
	    uvar_projectMetric ? "Projected":"Unprojected");
    printf (" %g \n", metric->data[INDEX_f0_f0]);
    printf (" %g  %g\n", metric->data[INDEX_f0_A], metric->data[INDEX_A_A]);
    printf (" %g  %g  %g\n", 
	    metric->data[INDEX_f0_D], metric->data[INDEX_A_D], metric->data[INDEX_D_D]);
    printf (" %g  %g  %g  %g\n\n",
	    metric->data[INDEX_f0_f1], metric->data[INDEX_A_f1], metric->data[INDEX_D_f1],
	    metric->data[INDEX_f1_f1]);

    LAL_CALL( LALDDestroyVector (&status, &metric), &status);
  } /* get metric */

  LAL_CALL ( LALDestroyUserVars (&status), &status);

  LALCheckMemoryLeaks(); 

  return 0;
} /* main */


/*----------------------------------------------------------------------*/
/** register all our "user-variables" */
void
initUserVars (LALStatus *stat)
{
  INITSTATUS( stat, "initUserVars", rcsid );
  ATTATCHSTATUSPTR (stat);

  /* set a few defaults */
#define EPHEM_YEARS  "00-04"
  uvar_ephemYear = (CHAR*)LALCalloc (1, strlen(EPHEM_YEARS)+1);
  strcpy (uvar_ephemYear, EPHEM_YEARS);

#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
  uvar_ephemDir = (CHAR*)LALCalloc (1, strlen(DEFAULT_EPHEMDIR)+1);
  strcpy (uvar_ephemDir, DEFAULT_EPHEMDIR);

  uvar_f1dot = 0.0;
  uvar_metricType = LAL_PMETRIC_COH_PTOLE_ANALYTIC;
  uvar_help = FALSE;

  uvar_startTime = 714180733;

  uvar_projectMetric = FALSE;

  /* now register all our user-variable */

  LALregBOOLUserVar(stat,	help,		'h', UVAR_HELP,     
		    "Print this help/usage message");
  LALregSTRINGUserVar(stat,	IFO,		'I', UVAR_REQUIRED, 
		      "Detector: GEO(0),LLO(1),LHO(2),NAUTILUS(3),VIRGO(4),TAMA(5),CIT(6)");

  LALregREALUserVar(stat,	Alpha,		'a', UVAR_REQUIRED,
		    "skyposition Alpha in radians, equatorial coords.");
  LALregREALUserVar(stat,	Delta, 		'd', UVAR_REQUIRED,
		    "skyposition Delta in radians, equatorial coords.");
  LALregREALUserVar(stat,	Freq, 		'f', UVAR_REQUIRED, 
		    "target frequency");
  LALregREALUserVar(stat,	f1dot, 		's', UVAR_OPTIONAL, 
		    "first spindown-value df/dt");
  LALregINTUserVar(stat,        metricType,     'M', UVAR_OPTIONAL, 
		   "Metric: 0=none,1=Ptole-analytic,2=Ptole-numeric, 3=exact");
  LALregBOOLUserVar(stat,	projectMetric,	 0,  UVAR_OPTIONAL,
		    "Project metric onto frequency-surface");
  LALregREALUserVar(stat,       startTime,      't', UVAR_OPTIONAL, 
		    "GPS start time of observation");
  LALregREALUserVar(stat,	duration,	'T', UVAR_REQUIRED, 
		    "Duration of observation in seconds");

  LALregSTRINGUserVar(stat,     ephemDir,       'E', UVAR_OPTIONAL, 
		      "Directory where Ephemeris files are located");
  LALregSTRINGUserVar(stat,     ephemYear,      'y', UVAR_OPTIONAL, 
		      "Year (or range of years) of ephemeris files to be used");


  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* initUserVars() */

/*----------------------------------------------------------------------*/
/** do some general initializations, 
 * e.g. load ephemeris-files (if required), setup detector etc
 */
void
initGeneral (LALStatus *lstat, ConfigVariables *cfg)
{

  INITSTATUS( lstat, "initGeneral", rcsid );
  ATTATCHSTATUSPTR (lstat);

  TRY ( LALFloatToGPS (lstat->statusPtr, &(cfg->startTimeGPS), &uvar_startTime), lstat);

  /* ---------- init ephemeris if needed ---------- */
  if ( uvar_metricType ==  LAL_PMETRIC_COH_EPHEM )
    {
      LALLeapSecFormatAndAcc formatAndAcc = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
      INT4 leap;

      if (LALUserVarWasSet (&uvar_ephemDir) )
	{
	  sprintf(cfg->EphemEarth, "%s/earth%s.dat", uvar_ephemDir, uvar_ephemYear);
	  sprintf(cfg->EphemSun, "%s/sun%s.dat", uvar_ephemDir, uvar_ephemYear);
	}
      else
	{
	  sprintf(cfg->EphemEarth, "earth%s.dat", uvar_ephemYear);
	  sprintf(cfg->EphemSun, "sun%s.dat", uvar_ephemYear);
	}

      cfg->ephemeris = (EphemerisData*) LALCalloc( 1, sizeof(EphemerisData) );
      cfg->ephemeris->ephiles.earthEphemeris = cfg->EphemEarth;
      cfg->ephemeris->ephiles.sunEphemeris = cfg->EphemSun;

      TRY (LALLeapSecs(lstat->statusPtr, &leap, &(cfg->startTimeGPS), &formatAndAcc), lstat);
      cfg->ephemeris->leap = leap;

      TRY (LALInitBarycenter (lstat->statusPtr, cfg->ephemeris), lstat);

  } /* end: init ephemeris data */


  /* ---------- initialize detector ---------- */
  if ( !strcmp (uvar_IFO, "GEO") || !strcmp (uvar_IFO, "0") ) 
    cfg->Detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  else if ( !strcmp (uvar_IFO, "LLO") || ! strcmp (uvar_IFO, "1") ) 
    cfg->Detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  else if ( !strcmp (uvar_IFO, "LHO") || !strcmp (uvar_IFO, "2") )
    cfg->Detector = lalCachedDetectors[LALDetectorIndexLHODIFF];
  else if ( !strcmp (uvar_IFO, "NAUTILUS") || !strcmp (uvar_IFO, "3")) {
    TRY (CreateNautilusDetector (lstat->statusPtr, &(cfg->Detector)), lstat);
  }
  else if ( !strcmp (uvar_IFO, "VIRGO") || !strcmp (uvar_IFO, "4") )
    cfg->Detector = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  else if ( !strcmp (uvar_IFO, "TAMA") || !strcmp (uvar_IFO, "5") )
    cfg->Detector = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
  else if ( !strcmp (uvar_IFO, "CIT") || !strcmp (uvar_IFO, "6") )
    cfg->Detector = lalCachedDetectors[LALDetectorIndexCIT40DIFF];
  else
    {
      LALPrintError ("\nUnknown detector. Currently allowed are 'GEO', 'LLO', 'LHO',"
		     " 'NAUTILUS', 'VIRGO', 'TAMA', 'CIT' or '0'-'6'\n\n");
      ABORT (lstat, GETMETRIC_EINPUT, GETMETRIC_MSGEINPUT);
    }

  DETATCHSTATUSPTR(lstat);
  RETURN(lstat);

} /* initGeneral() */

/*----------------------------------------------------------------------*/
/** Set up the \em LALDetector struct representing the NAUTILUS detector */
void
CreateNautilusDetector (LALStatus *lstat, LALDetector *Detector)
{
  /*   LALDetector Detector;  */
  LALFrDetector detector_params;
  LALDetectorType bar;
  LALDetector Detector1;

  INITSTATUS (lstat, "CreateNautilusDetector", rcsid);
  ATTATCHSTATUSPTR (lstat);

  bar=LALDETECTORTYPE_CYLBAR;
  strcpy(detector_params.name, "NAUTILUS");
  detector_params.vertexLongitudeRadians=12.67*LAL_PI/180.0;
  detector_params.vertexLatitudeRadians=41.82*LAL_PI/180.0;
  detector_params.vertexElevation=300.0;
  detector_params.xArmAltitudeRadians=0.0;
  detector_params.xArmAzimuthRadians=44.0*LAL_PI/180.0;

  TRY (LALCreateDetector(lstat->statusPtr, &Detector1, &detector_params, bar), lstat);
  
  *Detector=Detector1;

  DETATCHSTATUSPTR (lstat);
  RETURN (lstat);
  
} /* CreateNautilusDetector() */
