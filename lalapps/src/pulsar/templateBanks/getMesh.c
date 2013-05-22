/*
*  Copyright (C) 2005, 2006, 2007 Reinhard Prix
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
 * \author Reinhard Prix
 *
 * Standalone code to produce a 'mesh' in the parameter-space
 * Currently this is just some code to use InitDopplerScan()
 * directly outside of ComputeFStatistic, and can't do much more
 * than generating sky-grids.
 *
 * UserInput-parameters a compatible with ComputeFStatistic.
 */

/* ---------- includes ---------- */
#include <math.h>

#include <lalapps.h>

#include <lal/UserInput.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/LALInitBarycenter.h>
#include <lal/AVFactories.h>
#include <lal/SFTutils.h>
#include <lal/LogPrintf.h>

#include <lal/DopplerScan.h>

/* Error codes and messages */
#define GETMESH_ENORM 	0
#define GETMESH_ESUB  	1
#define GETMESH_EARG  	2
#define GETMESH_EBAD  	3
#define GETMESH_EFILE 	4
#define GETMESH_ENOARG 	5
#define GETMESH_EINPUT 	6
#define GETMESH_EMEM 	7
#define GETMESH_ENULL 	8

#define GETMESH_MSGENORM 	"Normal exit"
#define GETMESH_MSGESUB  	"Subroutine failed"
#define GETMESH_MSGEARG  	"Error parsing arguments"
#define GETMESH_MSGEBAD  	"Bad argument values"
#define GETMESH_MSGEFILE 	"File IO error"
#define GETMESH_MSGENOARG 	"Missing argument"
#define GETMESH_MSGEINPUT 	"Invalid user input"
#define GETMESH_MSGEMEM 	"Out of memory"
#define GETMESH_MSGENULL 	"Illegal NULL pointer in input"

/*---------- local defines ---------- */
#define TRUE (1==1)
#define FALSE (1==0)

/* ---------- some local types ---------- */
/* User variables */
typedef struct
{
  BOOLEAN help;

  CHAR* IFO;			/**< name of detector (LHO, LLO, GEO, VIRGO,..)*/
  REAL8 Alpha;			/**< skyposition Alpha: radians, equatorial coords. */
  REAL8 Delta;			/**< skyposition Delta: radians, equatorial coords. */
  REAL8 Freq;			/**< target-frequency */
  REAL8 f1dot;			/**< target 1. spindown-value df/dt */
  CHAR *ephemDir;		/**< directory of ephemeris-files */
  CHAR *ephemYear;		/**< year-range of ephemeris-file to use */

  REAL8 startTime;		/**< GPS start time of observation */
  REAL8 endTime;		/**< GPS end time of observation */
  REAL8 duration;		/**< OR alternatively: length of observation in seconds */
  BOOLEAN projectMetric;		/**< project out frequency-component of metric */

  REAL8 dAlpha;			/**< Alpha-resolution if GRID_FLAT or GRID_ISOTROPIC */
  REAL8 dDelta;			/**< Delta-resolution if ... */
  REAL8 AlphaBand;		/**< Sky-region interval in Alpha */
  REAL8 DeltaBand;		/**< Sky-region interval in Delta */
  REAL8 FreqBand;		/**< Frequency-band */
  REAL8 f1dotBand;		/**< spindown-band for f1dot */
  REAL8 dFreq;
  REAL8 df1dot;
  INT4  gridType;		/**< GRID_FLAT, GRID_ISOTROPIC, GRID_METRIC or GRID_FILE */
  INT4  metricType;		/**< if GRID_METRIC: what type of metric to use? */
  REAL8 metricMismatch;		/**< what's the maximal mismatch of the grid? */
  CHAR *skyRegion;		/**< alternative input: polygon describing skyregion */

  CHAR *skyGridFile;		/**< read in sky-grid if GRID_FILE */
  CHAR *outputSkyGrid;		/**< write sky-grid into this file */

  INT4 searchNeighbors;		/**< number of desired gridpoints/dimension around central point*/
  INT4 randomSeed;		/**< random-seed for searchNeighbors grid randomization */

  CHAR *mergedSFTFile;		/**< dummy for CFS-compatibility */

  INT4 numSkyPartitions;	/**< number of (roughly)equal partitions to split sky-grid into */
  INT4 partitionIndex;		/**< 0<= index < numSkyPartitions: index of sky-partition to generate */
} UserVariables_t;

typedef struct
{
  CHAR EphemEarth[512];		/**< filename of earth-ephemeris data */
  CHAR EphemSun[512];		/**< filename of sun-ephemeris data */
  LALDetector *Detector;        /**< detector of data to be searched */
  EphemerisData *ephemeris;	/**< ephemeris data (from LALInitBarycenter()) */
  LIGOTimeGPS startTimeGPS;	/**< starttime of observation */
  REAL8 duration;		/**< total observation-time spanned */
  DopplerRegion searchRegion;	/**< Doppler parameter-space to search */
} ConfigVariables;


/*---------- empty structs for initializations ----------*/
static const UserVariables_t empty_UserVariables;
static const ConfigVariables empty_ConfigVariables;
static const PtoleMetricIn empty_metricpar;

/* ---------- local prototypes ---------- */
void initUserVars (LALStatus *, UserVariables_t *uvar);

void initGeneral (LALStatus *, ConfigVariables *cfg, const UserVariables_t *uvar);
void checkUserInputConsistency (LALStatus *, const UserVariables_t *uvar);
void getSearchRegion (LALStatus *, DopplerRegion *searchRegion, const DopplerSkyScanInit *params, const UserVariables_t *uvar);
void setTrueRandomSeed(void);


extern int vrbflg;

/*----------------------------------------------------------------------
 * main function
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[])
{
  LALStatus status = blank_status;
  ConfigVariables config = empty_ConfigVariables;
  UserVariables_t uvar = empty_UserVariables;
  DopplerSkyScanInit scanInit = empty_DopplerSkyScanInit; /* init-structure for DopperScanner */
  DopplerSkyScanState thisScan = empty_DopplerSkyScanState; /* current state of the Doppler-scan */
  UINT4 nFreq, nf1dot;

  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register user-variables */
  LAL_CALL (initUserVars (&status, &uvar), &status);

  /* read cmdline & cfgfile  */
  LAL_CALL (LALUserVarReadAllInput (&status, argc,argv), &status);

  if (uvar.help) 	/* help requested: we're done */
    exit (0);

  /* check if user-input acutally made sense... */
  LAL_CALL ( checkUserInputConsistency(&status, &uvar), &status);

  /* do some general setup of the code */
  LAL_CALL (initGeneral (&status, &config, &uvar), &status);


  /* ---------- main-task: run InitDopplerSkyScan() ---------- */
  if (lalDebugLevel) printf ("\nSetting up template grid ...");
  scanInit.metricType = (LALPulsarMetricType) uvar.metricType;
  scanInit.dAlpha = uvar.dAlpha;
  scanInit.dDelta = uvar.dDelta;
  scanInit.gridType = (DopplerGridType) uvar.gridType;
  scanInit.metricMismatch = uvar.metricMismatch;
  scanInit.obsBegin = config.startTimeGPS;
  scanInit.obsDuration = config.duration;
  scanInit.projectMetric = uvar.projectMetric;
  scanInit.Detector = config.Detector;
  scanInit.ephemeris = config.ephemeris;       /* used by Ephemeris-based metric */
  scanInit.skyGridFile = uvar.skyGridFile;      /* if applicable */

  scanInit.numSkyPartitions = uvar.numSkyPartitions;
  scanInit.partitionIndex = uvar.partitionIndex;

  /* figure out searchRegion from UserInput and possibly --searchNeighbors */
  config.searchRegion = empty_DopplerRegion;	/* set to empty first */
  LAL_CALL ( getSearchRegion(&status, &(config.searchRegion), &scanInit, &uvar ), &status);
  scanInit.skyRegionString = config.searchRegion.skyRegionString;
  scanInit.Freq = config.searchRegion.fkdot[0] + config.searchRegion.fkdotBand[0];

  /* the following call generates a skygrid plus determines dFreq, df1dot,..
   * Currently the complete template-grid is more like a 'foliation' of
   * the parameter-space by skygrids, stacked along f and f1dot,
   * not a full 4D metric template-grid
   */
  LAL_CALL ( InitDopplerSkyScan( &status, &thisScan, &scanInit), &status);

  /* ---------- overload Frequency- and spindown-resolution if input by user ----------*/
  if ( LALUserVarWasSet( &uvar.dFreq ) )
    thisScan.dfkdot[0] = uvar.dFreq;

  if( LALUserVarWasSet( &uvar.df1dot) )
    thisScan.dfkdot[1] = uvar.df1dot;

  /* we write the sky-grid to disk? */
  if ( uvar.outputSkyGrid )
    {
      printf ("\nNow writing sky-grid into file '%s' ...", uvar.outputSkyGrid);
      LAL_CALL(writeSkyGridFile(&status, thisScan.skyGrid, uvar.outputSkyGrid ), &status);
      printf (" done.\n\n");
    }


  /* ----- output grid-info ----- */
  nFreq = (UINT4)(config.searchRegion.fkdotBand[0] / thisScan.dfkdot[0] + 1e-6) + 1;
  nf1dot =(UINT4)(config.searchRegion.fkdotBand[1]/ thisScan.dfkdot[1] + 1e-6) + 1;

  /* debug output */
  if ( lalDebugLevel )
    {
      printf ("DEBUG: Search-region:\n");
      printf ("       skyRegion = \"%s\"\n", scanInit.skyRegionString);
      printf ("       Freq in  = [%.16g, %.16g]\n",
	      config.searchRegion.fkdot[0],
	      config.searchRegion.fkdot[0] + config.searchRegion.fkdotBand[0]);
      printf ("       f1dot in = [%.16g, %.16g]\n",
	      config.searchRegion.fkdot[1],
	      config.searchRegion.fkdot[1] + config.searchRegion.fkdotBand[1]);

      printf ("\nDEBUG: actual grid-spacings: dFreq = %g, df1dot = %g\n\n",
	      thisScan.dfkdot[0], thisScan.dfkdot[1]);

      printf ("Templates: sky x Freq x f1dot = %d x %d x %d\n\n",
	      thisScan.numSkyGridPoints, nFreq, nf1dot );
    } /* debug-output */


  /* "official output": if NOT --searchNeighbors, just output the \
   * Total Number of templates:
   */
  if ( ! LALUserVarWasSet (&uvar.searchNeighbors) )
    {
      printf ("\n%g\n\n", (REAL8)thisScan.numSkyGridPoints * nFreq *  nf1dot );
    }
  /* output octave-compatible 'search-range' as a matrix:
   * [ Alpha, AlphaBand; Delta, DeltaBand; Freq, FreqBand; f1dot, f1dotBand ]
   */
  else
    {
      SkyRegion skyregion;
      DopplerRegion *searchRegion = &(config.searchRegion);
      REAL8 Alpha, AlphaBand, Delta, DeltaBand;
      REAL8 Freq, FreqBand, f1dot, f1dotBand;

      /* we need to parse the skyRegion string into a SkyRegion, which is easier to handle */
      LAL_CALL (ParseSkyRegionString(&status, &skyregion, searchRegion->skyRegionString), &status);
      if ( skyregion.vertices)
	LALFree ( skyregion.vertices );

      Alpha = skyregion.lowerLeft.longitude;
      AlphaBand = skyregion.upperRight.longitude - Alpha;
      Delta = skyregion.lowerLeft.latitude;
      DeltaBand = skyregion.upperRight.latitude - Delta;

      Freq = searchRegion->fkdot[0];
      FreqBand = searchRegion->fkdotBand[0];
      f1dot = searchRegion->fkdot[1];
      f1dotBand = searchRegion->fkdotBand[1];

      /* now write octave-compatible matrix containing the search-region */
      printf ("[%.12g, %.12g; %.12g, %.12g; %.12g, %.12g; %.12g, %.12g ]\n",
	      Alpha, AlphaBand,
	      Delta, DeltaBand,
	      Freq, FreqBand,
	      f1dot, f1dotBand );

    } /* if uvar.searchNeighbors */

  /* ----- clean up and exit ----- */
  /* Free DopplerSkyScan-stuff (grid) */
  thisScan.state = STATE_FINISHED;
  LAL_CALL ( FreeDopplerSkyScan(&status, &thisScan), &status);

  /* Free User-variables and contents */
  LAL_CALL ( LALDestroyUserVars (&status), &status);

  if ( config.searchRegion.skyRegionString )
    LALFree ( config.searchRegion.skyRegionString );

  LALFree ( config.Detector );

  LALCheckMemoryLeaks();

  return 0;
} /* main */


/*----------------------------------------------------------------------*/
/** register all our "user-variables" */
void
initUserVars (LALStatus *status, UserVariables_t *uvar)
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* set a few defaults */
#define EPHEM_YEARS  "00-19-DE405"
  uvar->ephemYear = (CHAR*)LALCalloc (1, strlen(EPHEM_YEARS)+1);
  strcpy (uvar->ephemYear, EPHEM_YEARS);

#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
  uvar->ephemDir = (CHAR*)LALCalloc (1, strlen(DEFAULT_EPHEMDIR)+1);
  strcpy (uvar->ephemDir, DEFAULT_EPHEMDIR);

  uvar->f1dot = 0.0;
  uvar->metricType = LAL_PMETRIC_COH_PTOLE_ANALYTIC;
  uvar->metricMismatch = 0.02;
  uvar->help = FALSE;

  uvar->startTime = 714180733;
  uvar->endTime = 0;

  uvar->projectMetric = TRUE;

  uvar->dAlpha = 0.001;
  uvar->dDelta = 0.001;

  uvar->AlphaBand = 0;
  uvar->DeltaBand = 0;
  uvar->FreqBand = 0;
  uvar->f1dotBand = 0;

  uvar->gridType = GRID_FLAT;
  uvar->metricType = LAL_PMETRIC_COH_PTOLE_ANALYTIC;

  uvar->metricMismatch = 0.02;
  uvar->skyRegion = NULL;
  uvar->skyGridFile = NULL;
  uvar->outputSkyGrid = NULL;

  uvar->searchNeighbors = 0;

  /* now register all user-variable */
  LALregBOOLUserStruct(status,	help,		'h', UVAR_HELP,		"Print this help/usage message");
  LALregSTRINGUserStruct(status,IFO,		'I', UVAR_REQUIRED, 	"Detector: H1, H2, L1, G1, ... ");

  LALregREALUserStruct(status,	Alpha,		'a', UVAR_OPTIONAL,	"skyposition Alpha in radians, equatorial coords.");
  LALregREALUserStruct(status,	Delta, 		'd', UVAR_OPTIONAL,	"skyposition Delta in radians, equatorial coords.");

  LALregREALUserStruct(status,	AlphaBand,      'z', UVAR_OPTIONAL, 	"Band in alpha (equatorial coordinates) in radians");
  LALregREALUserStruct(status,	DeltaBand,      'c', UVAR_OPTIONAL, 	"Band in delta (equatorial coordinates) in radians");
  LALregSTRINGUserStruct(status,skyRegion,      'R', UVAR_OPTIONAL, 	"ALTERNATIVE: sky-region polygon \"(a1,d1),(a2,d2),..(aN,dN)\"");

  LALregINTUserStruct(status,	numSkyPartitions, 0, UVAR_OPTIONAL,	"Number of (equi-)partitions to split skygrid into");
  LALregINTUserStruct(status,	partitionIndex,   0, UVAR_OPTIONAL,	"Index [0,numSkyPartitions-1] of sky-partition to generate");

  LALregREALUserStruct(status,	Freq, 		'f', UVAR_REQUIRED, 	"target frequency");
  LALregREALUserStruct(status,	f1dot, 		's', UVAR_OPTIONAL, 	"first spindown-value df/dt");
  LALregREALUserStruct(status,	FreqBand,       'b', UVAR_OPTIONAL, 	"Search frequency band in Hz");
  LALregREALUserStruct(status,  f1dotBand,      'm', UVAR_OPTIONAL, 	"Search-band for f1dot");

  LALregREALUserStruct(status,	dFreq,         	'r', UVAR_OPTIONAL, 	"Frequency resolution in Hz (default: 1/(2*Tsft*Nsft)");
  LALregREALUserStruct(status,	df1dot,        	'e', UVAR_OPTIONAL, 	"Resolution for f1dot (default: use metric or 1/(2*T^2))");

  LALregREALUserStruct(status,	startTime,       0, UVAR_OPTIONAL, 	"GPS start time of observation");
  LALregREALUserStruct(status, 	endTime,      	 0, UVAR_OPTIONAL,	"GPS end time of observation");
  LALregREALUserStruct(status,	duration,	'T', UVAR_OPTIONAL,	"Alternative: Duration of observation in seconds");

  LALregSTRINGUserStruct(status,ephemDir,       'E', UVAR_OPTIONAL,	"Directory where Ephemeris files are located");
  LALregSTRINGUserStruct(status,ephemYear,      'y', UVAR_OPTIONAL, 	"Year (or range of years) of ephemeris files to be used");

  /* ---------- */
  LALregINTUserStruct(status,	gridType,        0 , UVAR_OPTIONAL, 	"0=flat, 1=isotropic, 2=metric, 3=file, 4=SkyGridFILE+metric");
  LALregINTUserStruct(status, 	metricType,     'M', UVAR_OPTIONAL, 	"Metric: 0=none,1=Ptole-analytic,2=Ptole-numeric, 3=exact");
  LALregBOOLUserStruct(status,	projectMetric,	 0,  UVAR_OPTIONAL,	"Project metric onto frequency-surface");
  LALregREALUserStruct(status,	metricMismatch, 'X', UVAR_OPTIONAL, 	"Maximal mismatch for SKY-grid (adjust value for more dimensions)");
  LALregREALUserStruct(status,	dAlpha,         'l', UVAR_OPTIONAL, 	"Resolution in alpha (equatorial coordinates) in radians");
  LALregREALUserStruct(status,	dDelta,         'g', UVAR_OPTIONAL, 	"Resolution in delta (equatorial coordinates) in radians");
  LALregSTRINGUserStruct(status,skyGridFile,     0,  UVAR_OPTIONAL, 	"Load sky-grid from this file.");
  LALregSTRINGUserStruct(status,outputSkyGrid,   0,  UVAR_OPTIONAL, 	"Write sky-grid into this file.");
  LALregINTUserStruct(status, 	searchNeighbors, 0,  UVAR_OPTIONAL, 	"Find search-grid with roughly this many points/dimension");
  LALregINTUserStruct(status, 	randomSeed, 	 0,  UVAR_OPTIONAL, 	"Random-seed to use for searchNeighbors grid-randomization");

  LALregSTRINGUserStruct(status,mergedSFTFile, 	 0,  UVAR_DEVELOPER, 	"Dummy for CFS-compatibility ");

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* initUserVars() */

/*----------------------------------------------------------------------*/
/** do some general initializations,
 * e.g. load ephemeris-files (if required), setup detector etc
 */
void
initGeneral (LALStatus *status, ConfigVariables *cfg, const UserVariables_t *uvar)
{

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* ----- set up Tspan */
  XLALGPSSetREAL8(&(cfg->startTimeGPS), uvar->startTime);

  if ( LALUserVarWasSet (&uvar->endTime) )
    cfg->duration = uvar->endTime - uvar->startTime;
  else
    cfg->duration = uvar->duration;


  /* ---------- init ephemeris if needed ---------- */
  if ( uvar->metricType ==  LAL_PMETRIC_COH_EPHEM )
    {
      if (LALUserVarWasSet (&uvar->ephemDir) )
	{
	  sprintf(cfg->EphemEarth, "%s/earth%s.dat", uvar->ephemDir, uvar->ephemYear);
	  sprintf(cfg->EphemSun, "%s/sun%s.dat", uvar->ephemDir, uvar->ephemYear);
	}
      else
	{
	  sprintf(cfg->EphemEarth, "earth%s.dat", uvar->ephemYear);
	  sprintf(cfg->EphemSun, "sun%s.dat", uvar->ephemYear);
	}

      cfg->ephemeris = (EphemerisData*) LALCalloc( 1, sizeof(EphemerisData) );
      cfg->ephemeris->ephiles.earthEphemeris = cfg->EphemEarth;
      cfg->ephemeris->ephiles.sunEphemeris = cfg->EphemSun;

      TRY (LALInitBarycenter (status->statusPtr, cfg->ephemeris), status);

  } /* end: init ephemeris data */


  /* ---------- initialize detector ---------- */
  if ( (cfg->Detector = XLALGetSiteInfo ( uvar->IFO )) == NULL ) {
    ABORT (status, GETMESH_EINPUT, GETMESH_MSGEINPUT);
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* initGeneral() */

/*----------------------------------------------------------------------*/
/** Some general consistency-checks on user-input.
 * Throws an error plus prints error-message if problems are found.
 */
void
checkUserInputConsistency (LALStatus *status, const UserVariables_t *uvar)
{

  INITSTATUS(status);

  if (uvar->ephemYear == NULL)
    {
      XLALPrintError ("\nNo ephemeris year specified (option 'ephemYear')\n\n");
      ABORT (status, GETMESH_EINPUT, GETMESH_MSGEINPUT);
    }
  /* don't allow negative bands (for safty in griding-routines) */
  if ( (uvar->AlphaBand < 0) ||  (uvar->DeltaBand < 0) )
    {
      XLALPrintError ("\nNegative value of sky-bands not allowed (alpha or delta)!\n\n");
      ABORT (status, GETMESH_EINPUT, GETMESH_MSGEINPUT);
    }
  /* check for negative stepsizes in Freq, Alpha, Delta */
  if ( LALUserVarWasSet(&uvar->dAlpha) && (uvar->dAlpha < 0) )
    {
      XLALPrintError ("\nNegative value of stepsize dAlpha not allowed!\n\n");
      ABORT (status, GETMESH_EINPUT, GETMESH_MSGEINPUT);
    }
  if ( LALUserVarWasSet(&uvar->dDelta) && (uvar->dDelta < 0) )
    {
      XLALPrintError ("\nNegative value of stepsize dDelta not allowed!\n\n");
      ABORT (status, GETMESH_EINPUT, GETMESH_MSGEINPUT);
    }

  /* grid-related checks */
  {
    BOOLEAN haveAlphaBand = LALUserVarWasSet( &uvar->AlphaBand );
    BOOLEAN haveDeltaBand = LALUserVarWasSet( &uvar->DeltaBand );

    BOOLEAN haveSkyRegion, haveAlphaDelta, haveSkyGridFile, useSkyGridFile, haveMetric, useMetric;

    haveSkyRegion  = (uvar->skyRegion != NULL);
    haveAlphaDelta = (LALUserVarWasSet(&uvar->Alpha) && LALUserVarWasSet(&uvar->Delta) );
    haveSkyGridFile   = (uvar->skyGridFile != NULL);
    useSkyGridFile   = (uvar->gridType == GRID_FILE_SKYGRID) || (uvar->gridType == GRID_METRIC_SKYFILE);
    haveMetric     = (uvar->metricType > LAL_PMETRIC_NONE);
    useMetric     = (uvar->gridType == GRID_METRIC);

    if ( (haveAlphaBand && !haveDeltaBand) || (haveDeltaBand && !haveAlphaBand) )
      {
	XLALPrintError ("\nERROR: Need either BOTH (AlphaBand, DeltaBand) or NONE.\n\n");
        ABORT (status, GETMESH_EINPUT, GETMESH_MSGEINPUT);
      }

    if ( !useSkyGridFile && !(haveSkyRegion || haveAlphaDelta) )
      {
        XLALPrintError ("\nNeed sky-region: either use (Alpha,Delta) or skyRegion!\n\n");
        ABORT (status, GETMESH_EINPUT, GETMESH_MSGEINPUT);
      }
    if ( haveSkyRegion && haveAlphaDelta )
      {
        XLALPrintError ("\nOverdetermined sky-region: only use EITHER (Alpha,Delta)"
		       " OR skyRegion!\n\n");
        ABORT (status, GETMESH_EINPUT, GETMESH_MSGEINPUT);
      }
    if ( useSkyGridFile && !haveSkyGridFile )
      {
        XLALPrintError ("\nERROR: gridType=FILE, but no skyGridFile specified!\n\n");
        ABORT (status, GETMESH_EINPUT, GETMESH_MSGEINPUT);
      }
    if ( !useSkyGridFile && haveSkyGridFile )
      {
        LALWarning (status, "\nWARNING: skyGridFile was specified but not needed ..."
		    " will be ignored\n");
      }
    if ( useSkyGridFile && (haveSkyRegion || haveAlphaDelta) )
      {
        LALWarning (status, "\nWARNING: We are using skyGridFile, but sky-region was"
		    " also specified ... will be ignored!\n");
      }
    if ( !useMetric && haveMetric)
      {
        LALWarning (status, "\nWARNING: Metric was specified for non-metric grid..."
		    " will be ignored!\n");
      }
    if ( useMetric && !haveMetric)
      {
        XLALPrintError ("\nERROR: metric grid-type selected, but no metricType selected\n\n");
        ABORT (status, GETMESH_EINPUT, GETMESH_MSGEINPUT);
      }

  } /* grid-related checks */


  /* check that observation start+duration were specified correctly */
  {
    BOOLEAN haveEndTime = LALUserVarWasSet(&uvar->endTime);
    BOOLEAN haveDuration =LALUserVarWasSet(&uvar->duration);
    if ( !haveEndTime && !haveDuration )
      {
	XLALPrintError ("\nERROR: need either --endTime or --duration!\n\n");
        ABORT (status, GETMESH_EINPUT, GETMESH_MSGEINPUT);
      }
    if ( haveEndTime && haveDuration )
      {
	XLALPrintError ("\nERROR: can specify only one of --endTime or --duration!\n\n");
        ABORT (status, GETMESH_EINPUT, GETMESH_MSGEINPUT);
      }

    if ( haveEndTime && (uvar->endTime < uvar->startTime) )
      {
	XLALPrintError ("\nERROR: endTime must be later than startTime!\n\n");
        ABORT (status, GETMESH_EINPUT, GETMESH_MSGEINPUT);
      }

  } /* check Tspan */

  RETURN (status);
} /* checkUserInputConsistency() */


/** Determine the DopplerRegion in parameter-space to search over.
 *
 * Normally this is just given directly by the user and therefore trivial.
 *
 * However, for Monte-Carlo-runs testing the metric-grid we want to search
 * a (randomized) 'small' region around a given parameter-space point,
 * specifying only the (approx.) number of grid-points desired per dimension.
 * This is done by getMCDopplerCube(), which determines the appropriate
 * Bands in each direction using the metric.
 * This behavior is triggered by the user-input  --searchNeighbors.
 *
 * This function allows the user additionally to manually override these Bands.
 *
 */
void
getSearchRegion (LALStatus *status,		/**< pointer to LALStatus structure */
		 DopplerRegion *searchRegion,	/**< OUT: the DopplerRegion to search over */
		 const DopplerSkyScanInit *params,	/**< IN: DopplerSkyScan params might be needed */
		 const UserVariables_t *uvar	/**< [in] set of all user-input */
		 )
{

  DopplerRegion ret = empty_DopplerRegion;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( searchRegion, status, GETMESH_ENULL, GETMESH_MSGENULL);
  ASSERT ( searchRegion->skyRegionString == NULL, status,
	   GETMESH_EINPUT, GETMESH_MSGEINPUT);

  /* set defaults */
  ret.fkdot[0] = uvar->Freq;
  ret.fkdotBand[0] = uvar->FreqBand;

  ret.fkdot[1] = uvar->f1dot;
  ret.fkdotBand[1] = uvar->f1dotBand;

  /* 'normalize' if negative bands where given */
  if ( ret.fkdotBand[0] < 0 )
    {
      ret.fkdotBand[0] *= -1.0;
      ret.fkdot[0] -= ret.fkdotBand[0];
    }
  if ( ret.fkdotBand[1] < 0 )
    {
      ret.fkdotBand[1] *= -1.0;
      ret.fkdot[1] -= ret.fkdotBand[1];
    }



  /* if user specified the option -searchNeighbors=N, we generate an
   * automatic search-region of N grid-steps in each dimension
   * around the given search-point
   */
  if ( LALUserVarWasSet(&uvar->searchNeighbors) )
    {
      DopplerRegion cube = empty_DopplerRegion;
      PulsarDopplerParams signal_params = empty_PulsarDopplerParams;

      signal_params.Alpha = uvar->Alpha;
      signal_params.Delta = uvar->Delta;
      signal_params.fkdot[0]  = uvar->Freq;
      signal_params.fkdot[1] = uvar->f1dot;

      /* set random-seed for MC grid-randomization */
      if ( LALUserVarWasSet(&uvar->randomSeed) )
	srand(uvar->randomSeed);
      else
	setTrueRandomSeed();

      /* construct MC doppler-cube around signal-location */
      TRY ( getMCDopplerCube(status->statusPtr,
			     &cube, signal_params, uvar->searchNeighbors, params), status);

      /* free previous skyRegionString */
      if ( ret.skyRegionString )
	LALFree (ret.skyRegionString);

      /* overload defaults with automatic search-region */
      ret = cube;

    } /* if searchNeighbors */

  /* ---------- finally, the user can override the 'neighborhood' search-Region
   * if he explicitly specified search-bands in some (or all) directions.
   *
   * The motivation is to allow more selective control over the search-region,
   * e.g. by setting certain bands to zero.
   */
  {
    BOOLEAN haveSkyBands = LALUserVarWasSet(&uvar->AlphaBand) && LALUserVarWasSet(&uvar->DeltaBand);

    if ( haveSkyBands )     /* manually set skyRegion */
      {
	CHAR *str = NULL;
	TRY ( SkySquare2String( status->statusPtr, &str,
				uvar->Alpha, uvar->Delta,
				uvar->AlphaBand, uvar->DeltaBand), status);
	if ( ret.skyRegionString)
	  LALFree ( ret.skyRegionString );

	ret.skyRegionString = str;
      }
    else if ( uvar->skyRegion )
      {
	if ( ret.skyRegionString)
	  LALFree ( ret.skyRegionString );
	ret.skyRegionString = LALCalloc(1, strlen(uvar->skyRegion) + 1);
	strcpy ( ret.skyRegionString, uvar->skyRegion);
      }

    if ( LALUserVarWasSet(&uvar->FreqBand) )	/* manually set Frequency-interval */
      {
	ret.fkdot[0] = uvar->Freq;
	ret.fkdotBand[0] = uvar->FreqBand;
      }
    if ( LALUserVarWasSet(&uvar->f1dotBand) )	/* manually set spindown-interval */
      {
	ret.fkdot[1] = uvar->f1dot;
	ret.fkdotBand[1] = uvar->f1dotBand;
      }

  } /* user-override of search-bands */


  /* 'normalize' all spin-bands to be positive */
  if ( ret.fkdotBand[0] < 0 )
    {
      ret.fkdotBand[0]  *= -1.0;
      ret.fkdot[0]  -= ret.fkdotBand[0];
    }
  if ( ret.fkdotBand[1] < 0 )
    {
      ret.fkdotBand[1] *= -1.0;
      ret.fkdot[1] -= ret.fkdot[1];
    }

  *searchRegion = ret;	/* return the result */

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* getSearchRegion() */

/** set random-seed from /dev/urandom if possible, otherwise
 * from uninitialized local-var ;)
 */
void
setTrueRandomSeed(void)
{
  FILE *fpRandom;
  UINT4 seed;

  fpRandom = fopen("/dev/urandom", "r");	/* read Linux random-pool for seed */
  if ( fpRandom == NULL )
    {
      seed = (UINT4) ( 1e6 * XLALGetTimeOfDay() );
      XLALPrintError ("\nCould not open /dev/urandom ... using clock microseconds to set seed to %d.\n\n", seed );
    }
  else
    {
      if ( fread(&seed, sizeof(UINT4),1, fpRandom) != 1 )
        XLALPrintError ("\nCould not read from /dev/urandom ... using default seed.\n\n");
      fclose(fpRandom);
    }

  srand(seed);

  return;
} /* setTrueRandomSeed() */

