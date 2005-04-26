/**
 * \file getMesh.c
 *
 * \author{Reinhard Prix}
 *
 * Standalone code to produce a 'mesh' in the parameter-space 
 * Currently this is just some code to use InitDopplerScan()
 * directly outside of ComputeFStatistic, and can't do much more
 * than generating sky-grids. 
 *
 * UserInput-parameters a compatible with ComputeFStatistic.
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
#define GETMESH_ENORM 	0
#define GETMESH_ESUB  	1
#define GETMESH_EARG  	2
#define GETMESH_EBAD  	3
#define GETMESH_EFILE 	4
#define GETMESH_ENOARG 	5
#define GETMESH_EINPUT 	6
#define GETMESH_EMEM 	7


#define GETMESH_MSGENORM 	"Normal exit"
#define GETMESH_MSGESUB  	"Subroutine failed"
#define GETMESH_MSGEARG  	"Error parsing arguments"
#define GETMESH_MSGEBAD  	"Bad argument values"
#define GETMESH_MSGEFILE 	"File IO error"
#define GETMESH_MSGENOARG 	"Missing argument"
#define GETMESH_MSGEINPUT 	"Invalid user input"
#define GETMESH_MSGEMEM 	"Out of memory"

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
  DopplerRegion searchRegion;	/**< search-region to search over (here: just skyRegion) */
} ConfigVariables;

/*---------- empty structs for initializations ----------*/
static const ConfigVariables empty_ConfigVariables;
static const PtoleMetricIn empty_metricpar;

/* ---------- local prototypes ---------- */
void initUserVars (LALStatus *stat);
void initGeneral (LALStatus *lstat, ConfigVariables *cfg);
void CreateNautilusDetector (LALStatus *lstat, LALDetector *Detector);
void checkUserInputConsistency (LALStatus *lstat);
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
REAL8 uvar_metricMismatch;	/**< metric mismatch */
REAL8 uvar_startTime;		/**< GPS start time of observation */
REAL8 uvar_duration;		/**< length of observation in seconds */
BOOLEAN uvar_projectMetric;		/**< project out frequency-component of metric */

/* this are in addition to the above from getMetric.c */
REAL8 uvar_dAlpha;		/**< Alpha-resolution if GRID_FLAT or GRID_ISOTROPIC */
REAL8 uvar_dDelta;		/**< Delta-resolution if ... */
REAL8 uvar_AlphaBand;		/**< Sky-region interval in Alpha */
REAL8 uvar_DeltaBand;		/**< Sky-region interval in Delta */
INT4  uvar_gridType;		/**< GRID_FLAT, GRID_ISOTROPIC, GRID_METRIC or GRID_FILE */
INT4  uvar_metricType;		/**< if GRID_METRIC: what type of metric to use? */
REAL8 uvar_metricMismatch;	/**< what's the maximal mismatch of the grid? */
CHAR *uvar_skyRegion;		/**< alternative input: polygon describing skyregion */

CHAR *uvar_skyGridFile;		/**< read in sky-grid if GRID_FILE */
CHAR *uvar_outputSkyGrid;	/**< write sky-grid into this file */

extern int vrbflg;

/*----------------------------------------------------------------------
 * main function 
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[]) 
{
  LALStatus lstat = blank_status;
  ConfigVariables config = empty_ConfigVariables;
  DopplerScanInit scanInit = empty_DopplerScanInit; /* init-structure for DopperScanner */
  DopplerScanState thisScan = empty_DopplerScanState; /* current state of the Doppler-scan */

  lalDebugLevel = 0;  
  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register user-variables */
  LAL_CALL (LALGetDebugLevel (&lstat, argc, argv, 'v'), &lstat);
  LAL_CALL (initUserVars (&lstat), &lstat);	  

  /* read cmdline & cfgfile  */	
  LAL_CALL (LALUserVarReadAllInput (&lstat, argc,argv), &lstat);  

  if (uvar_help) 	/* help requested: we're done */
    exit (0);

  /* check if user-input acutally made sense... */
  LAL_CALL ( checkUserInputConsistency(&lstat), &lstat);

  /* do some general setup of the code */
  LAL_CALL (initGeneral (&lstat, &config), &lstat);


  /* ---------- main-task: run InitDopplerScan() ---------- */
  {
    if (lalDebugLevel) printf ("\nSetting up template grid ...");
    scanInit.metricType = (LALPulsarMetricType) uvar_metricType;
    scanInit.dAlpha = uvar_dAlpha;
    scanInit.dDelta = uvar_dDelta;
    scanInit.gridType = (DopplerGridType) uvar_gridType;
    scanInit.metricMismatch = uvar_metricMismatch;
    scanInit.obsBegin = config.startTimeGPS;
    scanInit.obsDuration = uvar_duration;
    scanInit.projectMetric = uvar_projectMetric;
    scanInit.fmax  = uvar_Freq;
    scanInit.Detector = &(config.Detector);
    scanInit.ephemeris = config.ephemeris;       /* used by Ephemeris-based metric */
    scanInit.skyGridFile = uvar_skyGridFile;      /* if applicable */

    scanInit.searchRegion = config.searchRegion;   /* initialize DopplerScan with search-region */
    
    /* the following call generates a skygrid plus determines dFreq, df1dot,..
     * Currently the complete template-grid is more like a 'foliation' of
     * the parameter-space by skygrids, stacked along f and f1dot, 
     * not a full 4D metric template-grid
     */
    LAL_CALL ( InitDopplerScan( &lstat, &thisScan, &scanInit), &lstat); 
    
    /* we write the sky-grid to disk? */
    if ( uvar_outputSkyGrid ) 
      {
	printf ("\nNow writing sky-grid into file '%s' ...", uvar_outputSkyGrid);
	LAL_CALL(writeSkyGridFile(&lstat, thisScan.grid, uvar_outputSkyGrid, &scanInit), &lstat);
	printf (" done.\n\n");
      }
    
  } /* ---------- run InitDopplerScan() ---------- */

  /* consider this finished */
  thisScan.state = STATE_FINISHED;

  /* Free DopplerScan-stuff (grid) */
  LAL_CALL ( FreeDopplerScan(&lstat, &thisScan), &lstat);

  /* Free User-variables and contents */
  LAL_CALL ( LALDestroyUserVars (&lstat), &lstat);

  if ( config.searchRegion.skyRegionString )
    LALFree ( config.searchRegion.skyRegionString );

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
  uvar_metricMismatch = 0.02;
  uvar_help = FALSE;

  uvar_startTime = 714180733;

  uvar_projectMetric = TRUE;

  uvar_dAlpha = 0.001;
  uvar_dDelta = 0.001;
  uvar_AlphaBand = 0;
  uvar_DeltaBand = 0;

  uvar_gridType = GRID_FLAT;
  uvar_metricType = LAL_PMETRIC_COH_PTOLE_ANALYTIC;

  uvar_metricMismatch = 0.02;
  uvar_skyRegion = NULL;
  uvar_skyGridFile = NULL;
  uvar_outputSkyGrid = NULL;

  /* now register all our user-variable */

  LALregBOOLUserVar(stat,	help,		'h', UVAR_HELP,     
		    "Print this help/usage message");
  LALregSTRINGUserVar(stat,	IFO,		'I', UVAR_REQUIRED, 
		      "Detector: GEO(0),LLO(1),LHO(2),NAUTILUS(3),VIRGO(4),TAMA(5),CIT(6)");

  LALregREALUserVar(stat,	Alpha,		'a', UVAR_OPTIONAL,
		    "skyposition Alpha in radians, equatorial coords.");
  LALregREALUserVar(stat,	Delta, 		'd', UVAR_OPTIONAL,
		    "skyposition Delta in radians, equatorial coords.");

  LALregREALUserVar(stat,       AlphaBand,      'z', UVAR_OPTIONAL, "Band in alpha (equatorial coordinates) in radians");
  LALregREALUserVar(stat,       DeltaBand,      'c', UVAR_OPTIONAL, "Band in delta (equatorial coordinates) in radians");
  LALregSTRINGUserVar(stat,     skyRegion,      'R', UVAR_OPTIONAL, "ALTERNATIVE: specify sky-region by polygon. Format: \"(a1,d1), (a2,d2),..(aN,dN)\"");
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

  /* ---------- */
  LALregINTUserVar(stat,        gridType,        0 , UVAR_OPTIONAL, "Template SKY-grid: 0=flat, 1=isotropic, 2=metric, 3=file");
  LALregINTUserVar(stat,        metricType,     'M', UVAR_OPTIONAL, "Metric: 0=none,1=Ptole-analytic,2=Ptole-numeric, 3=exact");
  LALregREALUserVar(stat,       metricMismatch, 'X', UVAR_OPTIONAL, "Maximal mismatch for SKY-grid (adjust value for more dimensions)");
  LALregREALUserVar(stat,       dAlpha,         'l', UVAR_OPTIONAL, "Resolution in alpha (equatorial coordinates) in radians");
  LALregREALUserVar(stat,       dDelta,         'g', UVAR_OPTIONAL, "Resolution in delta (equatorial coordinates) in radians");
  LALregSTRINGUserVar(stat,     skyGridFile,     0,  UVAR_OPTIONAL, "Load sky-grid from this file.");
  LALregSTRINGUserVar(stat,     outputSkyGrid,   0,  UVAR_OPTIONAL, "Write sky-grid into this file.");

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
      ABORT (lstat, GETMESH_EINPUT, GETMESH_MSGEINPUT);
    }

  /* determine search-region in parameter-space (in this case: just sky) */
  {
    BOOLEAN haveAlphaDelta = LALUserVarWasSet(&uvar_Alpha) && LALUserVarWasSet(&uvar_Delta);
    DopplerRegion ret = empty_DopplerRegion;

    if (uvar_skyRegion)
      {

	ret.skyRegionString = (CHAR*) LALCalloc(1, strlen( uvar_skyRegion ) + 1 );
	if ( ret.skyRegionString == NULL ) {
	  ABORT (lstat, GETMESH_EMEM, GETMESH_MSGEMEM);
	}
	strcpy (ret.skyRegionString, uvar_skyRegion);
      }
    else if (haveAlphaDelta)    /* parse this into a sky-region */
      {
	TRY ( SkySquare2String( lstat->statusPtr, &(ret.skyRegionString),
				uvar_Alpha, uvar_Delta, 
				uvar_AlphaBand, uvar_DeltaBand), lstat);
      }
    ret.Freq = uvar_Freq;
    ret.FreqBand = 0;
    ret.f1dot = uvar_f1dot;
    ret.f1dotBand = 0;

    cfg->searchRegion = ret;
  } /* ---------- get search-region ---------- */

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


/*----------------------------------------------------------------------*/
/** Some general consistency-checks on user-input.
 * Throws an error plus prints error-message if problems are found.
 */
void
checkUserInputConsistency (LALStatus *lstat)
{

  INITSTATUS (lstat, "checkUserInputConsistency", rcsid);  

  if (uvar_ephemYear == NULL)
    {
      LALPrintError ("\nNo ephemeris year specified (option 'ephemYear')\n\n");
      ABORT (lstat, GETMESH_EINPUT, GETMESH_MSGEINPUT);
    }      
  /* don't allow negative bands (for safty in griding-routines) */
  if ( (uvar_AlphaBand < 0) ||  (uvar_DeltaBand < 0) )
    {
      LALPrintError ("\nNegative value of sky-bands not allowed (alpha or delta)!\n\n");
      ABORT (lstat, GETMESH_EINPUT, GETMESH_MSGEINPUT);
    }
  /* check for negative stepsizes in Freq, Alpha, Delta */
  if ( LALUserVarWasSet(&uvar_dAlpha) && (uvar_dAlpha < 0) )
    {
      LALPrintError ("\nNegative value of stepsize dAlpha not allowed!\n\n");
      ABORT (lstat, GETMESH_EINPUT, GETMESH_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar_dDelta) && (uvar_dDelta < 0) )
    {
      LALPrintError ("\nNegative value of stepsize dDelta not allowed!\n\n");
      ABORT (lstat, GETMESH_EINPUT, GETMESH_MSGEINPUT);
    }

  /* grid-related checks */
  {
    BOOLEAN haveAlphaBand = LALUserVarWasSet( &uvar_AlphaBand );
    BOOLEAN haveDeltaBand = LALUserVarWasSet( &uvar_DeltaBand );

    BOOLEAN haveSkyRegion, haveAlphaDelta, haveGridFile, useGridFile, haveMetric, useMetric;

    haveSkyRegion  = (uvar_skyRegion != NULL);
    haveAlphaDelta = (LALUserVarWasSet(&uvar_Alpha) && LALUserVarWasSet(&uvar_Delta) );
    haveGridFile   = (uvar_skyGridFile != NULL);
    useGridFile   = (uvar_gridType == GRID_FILE);
    haveMetric     = (uvar_metricType > LAL_PMETRIC_NONE);
    useMetric     = (uvar_gridType == GRID_METRIC);

    if ( (haveAlphaBand && !haveDeltaBand) || (haveDeltaBand && !haveAlphaBand) )
      {
	LALPrintError ("\nERROR: Need either BOTH (AlphaBand, DeltaBand) or NONE.\n\n"); 
        ABORT (lstat, GETMESH_EINPUT, GETMESH_MSGEINPUT);
      }

    if ( !useGridFile && !(haveSkyRegion || haveAlphaDelta) )
      {
        LALPrintError ("\nNeed sky-region: either use (Alpha,Delta) or skyRegion!\n\n");
        ABORT (lstat, GETMESH_EINPUT, GETMESH_MSGEINPUT);
      }
    if ( haveSkyRegion && haveAlphaDelta )
      {
        LALPrintError ("\nOverdetermined sky-region: only use EITHER (Alpha,Delta)"
		       " OR skyRegion!\n\n");
        ABORT (lstat, GETMESH_EINPUT, GETMESH_MSGEINPUT);
      }
    if ( useGridFile && !haveGridFile )
      {
        LALPrintError ("\nERROR: gridType=FILE, but no skyGridFile specified!\n\n");
        ABORT (lstat, GETMESH_EINPUT, GETMESH_MSGEINPUT);  
      }
    if ( !useGridFile && haveGridFile )
      {
        LALWarning (lstat, "\nWARNING: skyGridFile was specified but not needed ..."
		    " will be ignored\n");
      }
    if ( useGridFile && (haveSkyRegion || haveAlphaDelta) )
      {
        LALWarning (lstat, "\nWARNING: We are using skyGridFile, but sky-region was"
		    " also specified ... will be ignored!\n");
      }
    if ( !useMetric && haveMetric) 
      {
        LALWarning (lstat, "\nWARNING: Metric was specified for non-metric grid..."
		    " will be ignored!\n");
      }
    if ( useMetric && !haveMetric) 
      {
        LALPrintError ("\nERROR: metric grid-type selected, but no metricType selected\n\n");
        ABORT (lstat, GETMESH_EINPUT, GETMESH_MSGEINPUT);      
      }

  } /* grid-related checks */

  RETURN (lstat);
} /* checkUserInputConsistency() */
