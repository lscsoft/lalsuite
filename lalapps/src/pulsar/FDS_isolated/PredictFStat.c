/*
 * Copyright (C) 2006 Iraj Gholami, Reinhard Prix
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

/*********************************************************************************/
/**
 * \author I. Gholami, R. Prix
 * \file
 * \ingroup pulsarApps
 * \brief
 * Calculate the *expected* (multi-IFO) F-statistic for pulsar GW signals, without actually
 * performing a search. The "F-statistic" was introduced in \ref JKS98 and Cutler-Schutz 2005.
 * Contrary to SemiAnalyticF this code can use (multi-IFO) SFTs to specify the startTime,
 * duration, detectors and noise-floors to use in the estimation.
 *
 */
#include "config.h"

/* System includes */
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

/* LAL-includes */
#include <lal/AVFactories.h>
#include <lal/LALInitBarycenter.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/ComputeFstat.h>
#include <lal/LALHough.h>
#include <lal/LogPrintf.h>

#include <lal/TransientCW_utils.h>

#include <lalapps.h>

/* local includes */

/*---------- DEFINES ----------*/

#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename */

#define EPHEM_YEARS  "00-19-DE405"	/**< default range: override with --ephemYear */

#define TRUE (1==1)
#define FALSE (1==0)

/*----- Error-codes -----*/
#define PREDICTFSTAT_ENULL 	1
#define PREDICTFSTAT_ESYS     	2
#define PREDICTFSTAT_EINPUT   	3
#define PREDICTFSTAT_EMEM   	4
#define PREDICTFSTAT_ENONULL 	5
#define PREDICTFSTAT_EXLAL	6

#define PREDICTFSTAT_MSGENULL 	"Arguments contained an unexpected null pointer"
#define PREDICTFSTAT_MSGESYS	"System call failed (probably file IO)"
#define PREDICTFSTAT_MSGEINPUT  "Invalid input"
#define PREDICTFSTAT_MSGEMEM   	"Out of memory. Bad."
#define PREDICTFSTAT_MSGENONULL "Output pointer is non-NULL"
#define PREDICTFSTAT_MSGEXLAL	"XLALFunction-call failed"

/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

#define SQ(x) ((x)*(x))

#define LAL_INT4_MAX 2147483647

/**
 * Configuration settings required for and defining a coherent pulsar search.
 * These are 'pre-processed' settings, which have been derived from the user-input.
 */
typedef struct {
  CHAR *dataSummary;            /**< descriptive string describing the data */
  REAL8 aPlus, aCross;		/**< internally always use Aplus, Across */
  AntennaPatternMatrix Mmunu;	/**< antenna-pattern matrix and normalization */
  UINT4 numSFTs;		/**< number of SFTs = Tobs/Tsft */
} ConfigVariables;

/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

ConfigVariables GV;		/**< global container for various derived configuration settings */

/* ----- User-variables: can be set from config-file or command-line */
typedef struct {
  BOOLEAN help;		/**< trigger output of help string */

  INT4 RngMedWindow;	/**< running-median window to use for noise-floor estimation */

  REAL8 aPlus;		/**< '+' polarization amplitude: aPlus  [alternative to {h0, cosi}: aPlus = 0.5*h0*(1+cosi^2)] */
  REAL8 aCross;		/**< 'x' polarization amplitude: aCross [alternative to {h0, cosi}: aCross= h0 * cosi] */
  REAL8 h0;		/**< overall GW amplitude h0 [alternative to {aPlus, aCross}] */
  REAL8 cosi;		/**< cos(inclination angle)  [alternative to {aPlus, aCross}] */
  REAL8 psi;		/**< polarization angle psi */
  REAL8 phi0;		/**< initial GW phase phi_0 in radians */
  REAL8 Freq;		/**< GW signal frequency */
  REAL8 Alpha;		/**< sky-position angle 'alpha', which is right ascencion in equatorial coordinates */
  REAL8 Delta;		/**< sky-position angle 'delta', which is declination in equatorial coordinates */

  BOOLEAN SignalOnly;	/**< switch wheter to assume Sh=1 instead of estimating noise-floors from SFTs */

  CHAR *IFO;		/**< GW detector short-name, only useful if not using v2-SFTs as input */

  CHAR *ephemDir;	/**< directory where ephemeris-files are to be found */
  CHAR *ephemYear;	/**< year-range of ephemeris-files to use */

  CHAR *DataFiles;	/**< SFT input-files to use to determine startTime, duration, IFOs and for noise-floor estimation */
  CHAR *outputFstat;	/**< output file to write F-stat estimation results into */
  BOOLEAN printFstat;	/**< print F-stat estimation results to terminal? */
  INT4 minStartTime;	/**< limit start-time of input SFTs to use */
  INT4 maxEndTime;	/**< limit end-time of input SFTs to use */

  CHAR *transientWindowType;	/**< name of transient window ('rect', 'exp',...) */
  REAL8 transientStartTime;	/**< GPS start-time of transient window */
  REAL8 transientTauDays;	/**< time-scale in days of transient window */

  REAL8 cosiota;	/* DEPRECATED in favor of cosi */

  BOOLEAN version;	/**< output version-info */

} UserInput_t;

static UserInput_t empty_UserInput;

/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);

void initUserVars (LALStatus *status, UserInput_t *uvar );
void InitPFS ( LALStatus *, ConfigVariables *cfg, const UserInput_t *uvar );
void InitEphemeris (LALStatus *, EphemerisData *edat, const CHAR *ephemDir, const CHAR *ephemYear);

/*---------- empty initializers ---------- */

/*----------------------------------------------------------------------*/
/* Main Function starts here */
/*----------------------------------------------------------------------*/
/**
 * MAIN function of PredictFStat code.
 * Calculates the F-statistic for a given position in the sky and detector
 * semi-analytically and outputs the final 2F value.
 */
int main(int argc,char *argv[])
{
  LALStatus status = blank_status;	/* initialize status */
  REAL8 rho2;	/* SNR^2 */

  UserInput_t uvar = empty_UserInput;
  CHAR *VCSInfoString;          /**< LAL + LALapps Git version string */

  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register all user-variable */
  LAL_CALL (initUserVars(&status, &uvar), &status);

  /* do ALL cmdline and cfgfile handling */
  LAL_CALL (LALUserVarReadAllInput(&status, argc, argv), &status);

  if (uvar.help)	/* if help was requested, we're done here */
    exit (0);


  if ( (VCSInfoString = XLALGetVersionString(0)) == NULL ) {
    XLALPrintError("XLALGetVersionString(0) failed.\n");
    exit(1);
  }

  if ( uvar.version ) {
    printf ("%s\n", VCSInfoString );
    exit(0);
  }

  /* Initialize code-setup */
  LAL_CALL ( InitPFS(&status, &GV, &uvar ), &status);

  { /* Calculating the F-Statistic */
    REAL8 al1, al2, al3;
    REAL8 Ap2 = SQ(GV.aPlus);
    REAL8 Ac2 = SQ(GV.aCross);
    REAL8 cos2psi2 = SQ( cos(2*uvar.psi) );
    REAL8 sin2psi2 = SQ( sin(2*uvar.psi) );

    al1 = Ap2 * cos2psi2 + Ac2 * sin2psi2;	/* A1^2 + A3^2 */
    al2 = Ap2 * sin2psi2 + Ac2 * cos2psi2;	/* A2^2 + A4^2 */
    al3 = ( Ap2 - Ac2 ) * sin(2.0*uvar.psi) * cos(2.0*uvar.psi);	/* A1 A2 + A3 A4 */

    /* SNR^2 */
    rho2 = GV.Mmunu.Sinv_Tsft * (GV.Mmunu.Ad * al1 + GV.Mmunu.Bd * al2 + 2.0 * GV.Mmunu.Cd * al3 );
  }

  if (uvar.printFstat)
    fprintf(stdout, "\n%.1f\n", 4.0 + rho2);

  /* output predicted Fstat-value into file, if requested */
  if (uvar.outputFstat)
    {
      FILE *fpFstat = NULL;
      CHAR *logstr = NULL;

      if ( (fpFstat = fopen (uvar.outputFstat, "wb")) == NULL)
	{
	  XLALPrintError ("\nError opening file '%s' for writing..\n\n", uvar.outputFstat);
	  return (PREDICTFSTAT_ESYS);
	}

      /* log search-footprint at head of output-file */
      LAL_CALL( LALUserVarGetLog (&status, &logstr,  UVAR_LOGFMT_CMDLINE ), &status );

      fprintf(fpFstat, "%%%% cmdline: %s\n", logstr );
      LALFree ( logstr );

      fprintf ( fpFstat, "%s\n", VCSInfoString );

      /* append 'dataSummary' */
      fprintf (fpFstat, "%s", GV.dataSummary );
      /* output E[2F] and std[2F] */
      fprintf (fpFstat, "twoF_expected = %g;\n", 4.0 + rho2);
      fprintf (fpFstat, "twoF_sigma    = %g;\n", sqrt( 4.0 * ( 2.0 + rho2 ) ) );

      /* output antenna-pattern matrix MNat_mu_nu = matrix(A, B, C) */
      {
	/* compute A = <a^2>, B=<b^2>, C=<ab> from the 'discretized versions Ad, Bc, Cd */
	REAL8 A = GV.Mmunu.Ad / GV.numSFTs;
	REAL8 B = GV.Mmunu.Bd / GV.numSFTs;
	REAL8 C = GV.Mmunu.Cd / GV.numSFTs;
	REAL8 D = A * B - C * C;
	fprintf (fpFstat, "A = %f;\n", A );
	fprintf (fpFstat, "B = %f;\n", B );
	fprintf (fpFstat, "C = %f;\n", C );
	fprintf (fpFstat, "D = %f;\n", D );
      }
      fclose (fpFstat);
    } /* if outputFstat */

  /* Free config-Variables and userInput stuff */
  LAL_CALL (LALDestroyUserVars (&status), &status);
  LALFree ( GV.dataSummary );
  XLALFree ( VCSInfoString );

  /* did we forget anything ? */
  LALCheckMemoryLeaks();

  return 0;

} /* main() */

/**
 * Register all our "user-variables" that can be specified from cmd-line and/or config-file.
 * Here we set defaults for some user-variables and register them with the UserInput module.
 */
void
initUserVars (LALStatus *status, UserInput_t *uvar )
{

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* set a few defaults */
  uvar->RngMedWindow = 50;	/* for running-median */

  uvar->ephemYear = LALCalloc (1, strlen(EPHEM_YEARS)+1);
  strcpy (uvar->ephemYear, EPHEM_YEARS);

#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
  uvar->ephemDir = LALCalloc (1, strlen(DEFAULT_EPHEMDIR)+1);
  strcpy (uvar->ephemDir, DEFAULT_EPHEMDIR);

  uvar->help = FALSE;
  uvar->outputFstat = NULL;
  uvar->printFstat = TRUE;

  uvar->minStartTime = 0;
  uvar->maxEndTime = LAL_INT4_MAX;

  uvar->phi0 = 0;

#define DEFAULT_TRANSIENT "none"
  uvar->transientWindowType = LALMalloc(strlen(DEFAULT_TRANSIENT)+1);
  strcpy ( uvar->transientWindowType, DEFAULT_TRANSIENT );

  /* register all our user-variables */
  LALregBOOLUserStruct(status,	help, 		'h', UVAR_HELP,     "Print this message");

  LALregREALUserStruct(status, 	aPlus,	 	 0 , UVAR_OPTIONAL, "'Plus' polarization amplitude: aPlus  [alternative to {h0, cosi}");
  LALregREALUserStruct(status,	aCross,  	 0 , UVAR_OPTIONAL, "'Cross' polarization amplitude: aCross [alternative to {h0, cosi}");
  LALregREALUserStruct(status,	h0,		's', UVAR_OPTIONAL, "Overall GW amplitude h0 [alternative to {aPlus, aCross}]");
  LALregREALUserStruct(status,	cosi,		'i', UVAR_OPTIONAL, "Inclination angle of rotation axis cos(iota) [alternative to {aPlus, aCross}]");

  LALregREALUserStruct(status,	psi,		 0, UVAR_REQUIRED, "Polarisation angle in radians");
  LALregREALUserStruct(status,	phi0,		 0, UVAR_OPTIONAL, "Initial GW phase phi0 in radians");

  LALregREALUserStruct(status,	Alpha,		'a', UVAR_REQUIRED, "Sky position alpha (equatorial coordinates) in radians");
  LALregREALUserStruct(status,	Delta,		'd', UVAR_REQUIRED, "Sky position delta (equatorial coordinates) in radians");
  LALregREALUserStruct(status,	Freq,		'F', UVAR_REQUIRED, "GW signal frequency (only used for noise-estimation in SFTs)");

  LALregSTRINGUserStruct(status,DataFiles, 	'D', UVAR_REQUIRED, "File-pattern specifying (multi-IFO) input SFT-files");
  LALregSTRINGUserStruct(status,IFO, 		'I', UVAR_OPTIONAL, "Detector-constraint: 'G1', 'L1', 'H1', 'H2' ...(useful for single-IFO v1-SFTs only!)");
  LALregSTRINGUserStruct(status,ephemDir, 	'E', UVAR_OPTIONAL, "Directory where Ephemeris files are located");
  LALregSTRINGUserStruct(status,ephemYear, 	'y', UVAR_OPTIONAL, "Year (or range of years) of ephemeris files to be used");
  LALregSTRINGUserStruct(status,outputFstat,     0,  UVAR_OPTIONAL, "Output-file for predicted F-stat value" );
  LALregBOOLUserStruct(status,printFstat,	 0,  UVAR_OPTIONAL, "Print predicted F-stat value to terminal" );

  LALregINTUserStruct ( status,	minStartTime, 	 0,  UVAR_OPTIONAL, "Earliest SFT-timestamp to include");
  LALregINTUserStruct ( status,	maxEndTime, 	 0,  UVAR_OPTIONAL, "Latest SFT-timestamps to include");

  LALregBOOLUserStruct(status,	SignalOnly,	'S', UVAR_OPTIONAL, "Assume Sh=1 instead of estimating noise-floor from SFTs");

  LALregBOOLUserStruct(status,	version,        'V', UVAR_SPECIAL,   "Output code version");

  LALregINTUserStruct(status,	RngMedWindow,	'k', UVAR_DEVELOPER, "Running-Median window size");
  LALregREALUserStruct(status,	cosiota,	 0 , UVAR_DEVELOPER, "[DEPRECATED] Use --cosi instead!");

  /* transient signal window properties (name, start, duration) */
  LALregSTRINGUserStruct(status, transientWindowType,  0, UVAR_OPTIONAL, "Name of transient signal window to use. ('none', 'rect', 'exp').");
  LALregREALUserStruct(status,   transientStartTime,    0, UVAR_OPTIONAL, "GPS start-time 't0' of transient signal window.");
  LALregREALUserStruct(status,   transientTauDays,   0, UVAR_OPTIONAL, "Timescale 'tau' of transient signal window in seconds.");

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* initUserVars() */

/** Load Ephemeris from ephemeris data-files  */
void
InitEphemeris (LALStatus * status,	/**< pointer to LALStatus structure */
	       EphemerisData *edat,	/**< [out] the ephemeris-data */
	       const CHAR *ephemDir,	/**< directory containing ephems */
	       const CHAR *ephemYear	/**< which years do we need? */
	       )
{
#define FNAME_LENGTH 1024
  CHAR EphemEarth[FNAME_LENGTH];	/* filename of earth-ephemeris data */
  CHAR EphemSun[FNAME_LENGTH];	/* filename of sun-ephemeris data */

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( edat, status, PREDICTFSTAT_ENULL, PREDICTFSTAT_MSGENULL );
  ASSERT ( ephemYear, status, PREDICTFSTAT_ENULL, PREDICTFSTAT_MSGENULL );

  if ( ephemDir )
    {
      snprintf(EphemEarth, FNAME_LENGTH, "%s/earth%s.dat", ephemDir, ephemYear);
      snprintf(EphemSun, FNAME_LENGTH, "%s/sun%s.dat", ephemDir, ephemYear);
    }
  else
    {
      snprintf(EphemEarth, FNAME_LENGTH, "earth%s.dat", ephemYear);
      snprintf(EphemSun, FNAME_LENGTH, "sun%s.dat",  ephemYear);
    }
  EphemEarth[FNAME_LENGTH-1]=0;
  EphemSun[FNAME_LENGTH-1]=0;

  /* NOTE: the 'ephiles' are ONLY ever used in LALInitBarycenter, which is
   * why we can use local variables (EphemEarth, EphemSun) to initialize them.
   */
  edat->ephiles.earthEphemeris = EphemEarth;
  edat->ephiles.sunEphemeris = EphemSun;

  TRY (LALInitBarycenter(status->statusPtr, edat), status);

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* InitEphemeris() */

/** Initialized Fstat-code: handle user-input and set everything up. */
void
InitPFS ( LALStatus *status, ConfigVariables *cfg, const UserInput_t *uvar )
{
  static const char *fn = "InitPFS()";

  SFTCatalog *catalog = NULL;
  SFTConstraints constraints = empty_SFTConstraints;
  SkyPosition skypos;

  LIGOTimeGPS startTime, endTime;
  REAL8 duration, Tsft;
  LIGOTimeGPS minStartTimeGPS, maxEndTimeGPS;

  EphemerisData *edat = NULL;		    	/* ephemeris data */
  MultiAMCoeffs *multiAMcoef = NULL;
  MultiPSDVector *multiRngmed = NULL;
  MultiNoiseWeights *multiNoiseWeights = NULL;
  MultiSFTVector *multiSFTs = NULL;	    	/* multi-IFO SFT-vectors */
  MultiDetectorStateSeries *multiDetStates = NULL; /* pos, vel and LMSTs for detector at times t_i */
  UINT4 X, i;


  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  { /* Check user-input consistency */
    BOOLEAN have_h0, have_cosi, have_cosiota, have_Ap, have_Ac;
    REAL8 cosi = 0;

    have_h0 = LALUserVarWasSet ( &uvar->h0 );
    have_cosi = LALUserVarWasSet ( &uvar->cosi );
    have_cosiota = LALUserVarWasSet ( &uvar->cosiota );
    have_Ap = LALUserVarWasSet ( &uvar->aPlus );
    have_Ac = LALUserVarWasSet ( &uvar->aCross );

    /* ----- handle cosi/cosiota ambiguity */
    if ( (have_cosi && have_cosiota)  ) {
      LogPrintf (LOG_CRITICAL, "Need EITHER --cosi [preferred] OR --cosiota [deprecated]!\n");
      ABORT ( status, PREDICTFSTAT_EINPUT, PREDICTFSTAT_MSGEINPUT );
    }
    if ( have_cosiota ) {
      cosi = uvar->cosiota;
      have_cosi = TRUE;
    }
    else if ( have_cosi ) {
      cosi = uvar->cosi;
      have_cosi = TRUE;
    }

    /* ----- handle {h0,cosi} || {aPlus,aCross} freedom ----- */
    if ( ( have_h0 && !have_cosi ) || ( !have_h0 && have_cosi ) )
      {
	LogPrintf (LOG_CRITICAL, "Need both (h0, cosi) to specify signal!\n");
	ABORT ( status, PREDICTFSTAT_EINPUT, PREDICTFSTAT_MSGEINPUT );
      }
    if ( ( have_Ap && !have_Ac) || ( !have_Ap && have_Ac ) )
      {
	LogPrintf (LOG_CRITICAL, "Need both (aPlus, aCross) to specify signal!\n");
	ABORT ( status, PREDICTFSTAT_EINPUT, PREDICTFSTAT_MSGEINPUT );
      }
    if ( have_h0 && have_Ap )
      {
	LogPrintf (LOG_CRITICAL, "Overdetermined: specify EITHER (h0,cosi) OR (aPlus,aCross)!\n");
	ABORT ( status, PREDICTFSTAT_EINPUT, PREDICTFSTAT_MSGEINPUT );
      }
    /* ----- internally we always use Aplus, Across */
    if ( have_h0 )
      {
	cfg->aPlus = 0.5 * uvar->h0 * ( 1.0 + SQ( cosi) );
	cfg->aCross = uvar->h0 * uvar->cosi;
      }
    else
      {
	cfg->aPlus = uvar->aPlus;
	cfg->aCross = uvar->aCross;
      }
  }/* check user-input */


  /* ----- prepare SFT-reading ----- */
  if ( LALUserVarWasSet ( &uvar->IFO ) )
    if ( (constraints.detector = XLALGetChannelPrefix ( uvar->IFO )) == NULL ) {
      ABORT ( status,  PREDICTFSTAT_EINPUT,  PREDICTFSTAT_MSGEINPUT);
    }

  minStartTimeGPS.gpsSeconds = uvar->minStartTime;
  minStartTimeGPS.gpsNanoSeconds = 0;
  maxEndTimeGPS.gpsSeconds = uvar->maxEndTime;
  maxEndTimeGPS.gpsNanoSeconds = 0;
  constraints.startTime = &minStartTimeGPS;
  constraints.endTime = &maxEndTimeGPS;

  /* ----- get full SFT-catalog of all matching (multi-IFO) SFTs */
  LogPrintf (LOG_DEBUG, "Finding all SFTs to load ... ");
  TRY ( LALSFTdataFind ( status->statusPtr, &catalog, uvar->DataFiles, &constraints ), status);
  LogPrintfVerbatim (LOG_DEBUG, "done. (found %d SFTs)\n", catalog->length);
  if ( constraints.detector )
    LALFree ( constraints.detector );

  if ( catalog->length == 0 )
    {
      LogPrintf (LOG_CRITICAL, "No matching SFTs for pattern '%s'!\n", uvar->DataFiles );
      ABORT ( status,  PREDICTFSTAT_EINPUT,  PREDICTFSTAT_MSGEINPUT);
    }

  /* ----- deduce start- and end-time of the observation spanned by the data */
  {
    GV.numSFTs = catalog->length;	/* total number of SFTs */
    Tsft = 1.0 / catalog->data[0].header.deltaF;
    startTime = catalog->data[0].header.epoch;
    endTime   = catalog->data[GV.numSFTs-1].header.epoch;
    XLALGPSAdd(&endTime, Tsft);
    duration = GPS2REAL8(endTime) - GPS2REAL8 (startTime);
  }

  { /* ----- load ephemeris-data ----- */
    CHAR *ephemDir;

    edat = LALCalloc(1, sizeof(EphemerisData));
    if ( LALUserVarWasSet ( &uvar->ephemDir ) )
      ephemDir = uvar->ephemDir;
    else
      ephemDir = NULL;
    TRY( InitEphemeris (status->statusPtr, edat, ephemDir, uvar->ephemYear ), status);
  }

  {/* ----- load the multi-IFO SFT-vectors ----- */
    UINT4 wings = uvar->RngMedWindow/2 + 10;   /* extra frequency-bins needed for rngmed */
    REAL8 fMax = uvar->Freq + 1.0 * wings / Tsft;
    REAL8 fMin = uvar->Freq - 1.0 * wings / Tsft;

    LogPrintf (LOG_DEBUG, "Loading SFTs ... ");
    TRY ( LALLoadMultiSFTs ( status->statusPtr, &multiSFTs, catalog, fMin, fMax ), status );
    LogPrintfVerbatim (LOG_DEBUG, "done.\n");
    TRY ( LALDestroySFTCatalog ( status->statusPtr, &catalog ), status );
  }

  TRY ( LALNormalizeMultiSFTVect (status->statusPtr, &multiRngmed, multiSFTs, uvar->RngMedWindow ), status);
  TRY ( LALComputeMultiNoiseWeights (status->statusPtr, &multiNoiseWeights, multiRngmed, uvar->RngMedWindow, 0 ), status );

  /* correctly handle the --SignalOnly case:
   * set noise-weights to 1, and
   * set Sh->1 (single-sided)
   */
  if ( uvar->SignalOnly )
    {
      multiNoiseWeights->Sinv_Tsft = Tsft;
      for ( X=0; X < multiNoiseWeights->length; X ++ )
        for ( i=0; i < multiNoiseWeights->data[X]->length; i ++ )
          multiNoiseWeights->data[X]->data[i] = 1.0;
    }

  /* ----- handle transient-signal window if given ----- */
  if ( LALUserVarWasSet ( &uvar->transientWindowType ) && strcmp ( uvar->transientWindowType, "none") )
    {
      transientWindow_t transientWindow;	/**< properties of transient-signal window */
      MultiLIGOTimeGPSVector *mTS;

      if ( !strcmp ( uvar->transientWindowType, "rect" ) )
        transientWindow.type = TRANSIENT_RECTANGULAR;		/* rectangular window [t0, t0+tau] */
      else if ( !strcmp ( uvar->transientWindowType, "exp" ) )
        transientWindow.type = TRANSIENT_EXPONENTIAL;		/* exponential decay window e^[-(t-t0)/tau for t>t0, 0 otherwise */
      else
	{
	  XLALPrintError ("Illegal transient window '%s' specified: valid are 'none', 'rect' or 'exp'\n", uvar->transientWindowType);
          ABORT ( status, PREDICTFSTAT_EINPUT, PREDICTFSTAT_MSGEINPUT );
	}

      if ( LALUserVarWasSet ( &uvar->transientStartTime ) )
        transientWindow.t0 = uvar->transientStartTime;
      else
        transientWindow.t0 = XLALGPSGetREAL8( &startTime ); /* if not set, default window startTime == startTime here */

      transientWindow.tau  = uvar->transientTauDays;

      if ( (mTS = XLALExtractMultiTimestampsFromSFTs ( multiSFTs )) == NULL ) {
        XLALPrintError ("%s: failed to XLALExtractMultiTimestampsFromSFTs() from SFTs. xlalErrno = %d.\n", fn, xlalErrno );
        ABORT ( status, PREDICTFSTAT_EXLAL, PREDICTFSTAT_MSGEXLAL );
      }

      if ( XLALApplyTransientWindow2NoiseWeights ( multiNoiseWeights, mTS, transientWindow ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: XLALApplyTransientWindow2NoiseWeights() failed! xlalErrno = %d\n", fn, xlalErrno );
        ABORT ( status, PREDICTFSTAT_EXLAL, PREDICTFSTAT_MSGEXLAL );
      }

      XLALDestroyMultiTimestamps ( mTS );

    } /* apply transient window to noise-weights */


  /* ----- obtain the (multi-IFO) 'detector-state series' for all SFTs ----- */
  TRY (LALGetMultiDetectorStates( status->statusPtr, &multiDetStates, multiSFTs, edat), status );

  /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
  skypos.longitude = uvar->Alpha;
  skypos.latitude = uvar->Delta;
  skypos.system = COORDINATESYSTEM_EQUATORIAL;
  TRY (LALNormalizeSkyPosition ( status->statusPtr, &skypos, &skypos), status);

  TRY ( LALGetMultiAMCoeffs ( status->statusPtr, &multiAMcoef, multiDetStates, skypos ), status);
  /* noise-weighting of Antenna-patterns and compute A,B,C */
  if ( XLALWeightMultiAMCoeffs ( multiAMcoef, multiNoiseWeights ) != XLAL_SUCCESS ) {
    LogPrintf (LOG_CRITICAL, "XLALWeightMultiAMCoeffs() failed with error = %d\n\n", xlalErrno );
    ABORT ( status, PREDICTFSTAT_EXLAL, PREDICTFSTAT_MSGEXLAL );
  }

  /* OK: we only need the antenna-pattern matrix M_mu_nu */
  cfg->Mmunu = multiAMcoef->Mmunu;

  /* ----- produce a log-string describing the data-specific setup ----- */
  {
    struct tm utc;
    time_t tp;
    CHAR dateStr[512], line[512], summary[1024];
    UINT4 j, numDet;
    numDet = multiSFTs->length;
    tp = time(NULL);
    sprintf (summary, "%%%% Date: %s", asctime( gmtime( &tp ) ) );
    strcat (summary, "%% Loaded SFTs: [ " );
    for ( j=0; j < numDet; j ++ ) {
      sprintf (line, "%s:%d%s",  multiSFTs->data[j]->data->name, multiSFTs->data[j]->length,
	       (j < numDet - 1)?", ":" ]\n");
      strcat ( summary, line );
    }
    utc = *XLALGPSToUTC( &utc, (INT4)GPS2REAL8(startTime) );
    strcpy ( dateStr, asctime(&utc) );
    dateStr[ strlen(dateStr) - 1 ] = 0;
    sprintf (line, "%%%% Start GPS time tStart = %12.3f    (%s GMT)\n", GPS2REAL8(startTime), dateStr);
    strcat ( summary, line );
    sprintf (line, "%%%% Total time spanned    = %12.3f s  (%.1f hours)\n", duration, duration/3600 );
    strcat ( summary, line );

    if ( (cfg->dataSummary = LALCalloc(1, strlen(summary) + 1 )) == NULL ) {
      ABORT (status, PREDICTFSTAT_EMEM, PREDICTFSTAT_MSGEMEM);
    }
    strcpy ( cfg->dataSummary, summary );

    LogPrintfVerbatim( LOG_DEBUG, cfg->dataSummary );
  } /* write dataSummary string */

  /* free everything not needed any more */
  TRY ( LALDestroyMultiPSDVector (status->statusPtr, &multiRngmed ), status );
  TRY ( LALDestroyMultiNoiseWeights (status->statusPtr, &multiNoiseWeights ), status );
  TRY ( LALDestroyMultiSFTVector (status->statusPtr, &multiSFTs ), status );
  XLALDestroyMultiDetectorStateSeries ( multiDetStates );
  XLALDestroyMultiAMCoeffs ( multiAMcoef );

  /* Free ephemeris data */
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);


  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* InitPFS() */
