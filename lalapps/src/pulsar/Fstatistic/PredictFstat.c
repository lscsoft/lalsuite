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
 * \ingroup lalapps_pulsar_Fstatistic
 * \brief
 * Calculate the *expected* (multi-IFO) F-statistic for pulsar GW signals, without actually
 * performing a search. The "F-statistic" was introduced in \cite JKS98 and Cutler-Schutz 2005.
 * Contrary to SemiAnalyticF this code can use (multi-IFO) SFTs to specify the startTime,
 * duration, detectors and noise-floors to use in the estimation.
 *
 */
#include "config.h"

/* System includes */
#include <stdio.h>

/* LAL-includes */
#include <lal/LALString.h>
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

  LALStringVector* assumeSqrtSX;/**< Assume stationary Gaussian noise with detector noise-floors sqrt{SX}" */
  BOOLEAN SignalOnly;	/**< DEPRECATED: ALTERNATIVE switch to assume Sh=1 instead of estimating noise-floors from SFTs */

  CHAR *IFO;		/**< GW detector short-name, only useful if not using v2-SFTs as input */

  CHAR *ephemEarth;	/**< Earth ephemeris file to use */
  CHAR *ephemSun;	/**< Sun ephemeris file to use */

  CHAR *DataFiles;	/**< SFT input-files to use to determine startTime, duration, IFOs and for noise-floor estimation */
  CHAR *outputFstat;	/**< output file to write F-stat estimation results into */
  BOOLEAN printFstat;	/**< print F-stat estimation results to terminal? */
  INT4 minStartTime;	/**< Only use SFTs with timestamps starting from (including) this GPS time */
  INT4 maxStartTime;	/**< Only use SFTs with timestamps up to (excluding) this GPS time */

  CHAR *transientWindowType;	/**< name of transient window ('rect', 'exp',...) */
  REAL8 transientStartTime;	/**< GPS start-time of transient window */
  REAL8 transientTauDays;	/**< time-scale in days of transient window */

  BOOLEAN version;	/**< output version-info */

} UserInput_t;

/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);

int initUserVars ( UserInput_t *uvar );
int InitPFS ( ConfigVariables *cfg, const UserInput_t *uvar );

/*---------- empty initializers ---------- */

/*----------------------------------------------------------------------*/
/* Main Function starts here */
/*----------------------------------------------------------------------*/
/**
 * MAIN function of PredictFstat code.
 * Calculates the F-statistic for a given position in the sky and detector
 * semi-analytically and outputs the final 2F value.
 */
int main(int argc,char *argv[])
{
  REAL8 rho2;	/* SNR^2 */

  UserInput_t XLAL_INIT_DECL(uvar);
  CHAR *VCSInfoString;          /**< LAL + LALapps Git version string */

  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register all user-variable */
  XLAL_CHECK_MAIN ( initUserVars( &uvar) == XLAL_SUCCESS, XLAL_EFUNC );

  /* do ALL cmdline and cfgfile handling */
  BOOLEAN should_exit = 0;
  XLAL_CHECK( XLALUserVarReadAllInput( &should_exit, argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit ) {
    exit (1);
  }

  XLAL_CHECK_MAIN ( (VCSInfoString = XLALGetVersionString(0)) != NULL, XLAL_EFUNC );

  if ( uvar.version ) {
    printf ("%s\n", VCSInfoString );
    exit(0);
  }

  /* Initialize code-setup */
  XLAL_CHECK_MAIN ( InitPFS ( &GV, &uvar ) == XLAL_SUCCESS, XLAL_EFUNC );

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

  if (uvar.printFstat) {
    fprintf(stdout, "\n%.1f\n", 4.0 + rho2);
  }

  /* output predicted Fstat-value into file, if requested */
  if (uvar.outputFstat)
    {
      FILE *fpFstat = NULL;
      CHAR *logstr = NULL;

      XLAL_CHECK_MAIN ( (fpFstat = fopen (uvar.outputFstat, "wb")) != NULL, XLAL_ESYS, "\nError opening file '%s' for writing..\n\n", uvar.outputFstat );

      /* log search-footprint at head of output-file */
      XLAL_CHECK_MAIN ( (logstr = XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE )) != NULL, XLAL_EFUNC );

      fprintf(fpFstat, "%%%% cmdline: %s\n", logstr );
      XLALFree ( logstr );

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
  XLALDestroyUserVars();
  XLALFree ( GV.dataSummary );
  XLALFree ( VCSInfoString );

  /* did we forget anything ? */
  LALCheckMemoryLeaks();

  return 0;

} /* main() */

/**
 * Register all our "user-variables" that can be specified from cmd-line and/or config-file.
 * Here we set defaults for some user-variables and register them with the UserInput module.
 */
int
initUserVars ( UserInput_t *uvar )
{
  XLAL_CHECK ( uvar != NULL, XLAL_EINVAL );

  /* set a few defaults */
  uvar->RngMedWindow = 50;	/* for running-median */

  uvar->ephemEarth = XLALStringDuplicate("earth00-19-DE405.dat.gz");
  uvar->ephemSun = XLALStringDuplicate("sun00-19-DE405.dat.gz");

  uvar->outputFstat = NULL;
  uvar->printFstat = 1;

  uvar->minStartTime = 0;
  uvar->maxStartTime = LAL_INT4_MAX;

  uvar->assumeSqrtSX = NULL;
  uvar->SignalOnly = 0;

  uvar->phi0 = 0;
  uvar->transientWindowType = XLALStringDuplicate ( "none" );

  /* register all our user-variables */
  XLALRegisterUvarMember( aPlus, 	 REAL8, 0 , OPTIONAL, "'Plus' polarization amplitude: aPlus  [alternative to {h0, cosi}");
  XLALRegisterUvarMember( aCross,  	 REAL8, 0 , OPTIONAL, "'Cross' polarization amplitude: aCross [alternative to {h0, cosi}");
  XLALRegisterUvarMember( h0,		REAL8, 's', OPTIONAL, "Overall GW amplitude h0 [alternative to {aPlus, aCross}]");
  XLALRegisterUvarMember( cosi,		REAL8, 'i', OPTIONAL, "Inclination angle of rotation axis cos(iota) [alternative to {aPlus, aCross}]");

  XLALRegisterUvarMember( psi,		 REAL8, 0, REQUIRED, "Polarisation angle in radians");
  XLALRegisterUvarMember( phi0,		 REAL8, 0, OPTIONAL, "Initial GW phase phi0 in radians");

  XLALRegisterUvarMember( Alpha,	REAL8, 'a', REQUIRED, "Sky position alpha (equatorial coordinates) in radians");
  XLALRegisterUvarMember( Delta,	REAL8, 'd', REQUIRED, "Sky position delta (equatorial coordinates) in radians");
  XLALRegisterUvarMember( Freq,		REAL8, 'F', REQUIRED, "GW signal frequency (only used for noise-estimation in SFTs)");

  XLALRegisterUvarMember( DataFiles, 	STRING, 'D', REQUIRED, "File-pattern specifying (multi-IFO) input SFT-files");
  XLALRegisterUvarMember( IFO, 		STRING, 'I', OPTIONAL, "Detector-constraint: 'G1', 'L1', 'H1', 'H2' ...(useful for single-IFO v1-SFTs only!)");
  XLALRegisterUvarMember( ephemEarth, 	 STRING, 0,  OPTIONAL, "Earth ephemeris file to use");
  XLALRegisterUvarMember( ephemSun, 	 STRING, 0,  OPTIONAL, "Sun ephemeris file to use");
  XLALRegisterUvarMember( outputFstat,     STRING, 0,  OPTIONAL, "Output-file for predicted F-stat value" );
  XLALRegisterUvarMember( printFstat,	 BOOLEAN, 0,  OPTIONAL, "Print predicted F-stat value to terminal" );

  XLALRegisterUvarMember( minStartTime, 	 INT4, 0,  OPTIONAL, "Only use SFTs with timestamps starting from (including) this GPS time");
  XLALRegisterUvarMember( maxStartTime, 	 INT4, 0,  OPTIONAL, "Only use SFTs with timestamps up to (excluding) this GPS time");

  XLALRegisterUvarMember( assumeSqrtSX,	 STRINGVector, 0,  OPTIONAL, "Don't estimate noise-floors but assume (stationary) per-IFO sqrt{SX} (if single value: use for all IFOs)");
  XLALRegisterUvarMember( SignalOnly,	BOOLEAN, 'S', DEVELOPER,"DEPRECATED ALTERNATIVE: Don't estimate noise-floors but assume sqrtSX=1 instead");

  XLALRegisterUvarMember( version,        BOOLEAN, 'V', SPECIAL,  "Output code version");

  XLALRegisterUvarMember( RngMedWindow,	INT4, 'k', DEVELOPER, "Running-Median window size");

  /* transient signal window properties (name, start, duration) */
  XLALRegisterUvarMember( transientWindowType,  STRING, 0, OPTIONAL, "Name of transient signal window to use. ('none', 'rect', 'exp').");
  XLALRegisterUvarMember( transientStartTime,    REAL8, 0, OPTIONAL, "GPS start-time 't0' of transient signal window.");
  XLALRegisterUvarMember( transientTauDays,   REAL8, 0, OPTIONAL, "Timescale 'tau' of transient signal window in seconds.");

  return XLAL_SUCCESS;

} /* initUserVars() */

/** Initialized Fstat-code: handle user-input and set everything up. */
int
InitPFS ( ConfigVariables *cfg, const UserInput_t *uvar )
{
  XLAL_CHECK ( (cfg != NULL) && (uvar != NULL), XLAL_EINVAL );

  SFTCatalog *catalog = NULL;
  SFTConstraints XLAL_INIT_DECL(constraints);
  SkyPosition skypos;

  LIGOTimeGPS minStartTimeGPS, maxStartTimeGPS;

  EphemerisData *edat = NULL;		    	/* ephemeris data */
  MultiAMCoeffs *multiAMcoef = NULL;
  MultiDetectorStateSeries *multiDetStates = NULL; /* pos, vel and LMSTs for detector at times t_i */

  { /* Check user-input consistency */
    BOOLEAN have_h0   = XLALUserVarWasSet ( &uvar->h0 );
    BOOLEAN have_cosi = XLALUserVarWasSet ( &uvar->cosi );
    BOOLEAN have_Ap   = XLALUserVarWasSet ( &uvar->aPlus );
    BOOLEAN have_Ac   = XLALUserVarWasSet ( &uvar->aCross );

    /* ----- handle {h0,cosi} || {aPlus,aCross} freedom ----- */
    if ( ( have_h0 && !have_cosi ) || ( !have_h0 && have_cosi ) ) {
      XLAL_ERROR ( XLAL_EINVAL, "Need both (h0, cosi) to specify signal!\n");
    }
    if ( ( have_Ap && !have_Ac) || ( !have_Ap && have_Ac ) ) {
      XLAL_ERROR ( XLAL_EINVAL, "Need both (aPlus, aCross) to specify signal!\n");
    }
    if ( have_h0 && have_Ap ) {
      XLAL_ERROR ( XLAL_EINVAL, "Overdetermined: specify EITHER (h0,cosi) OR (aPlus,aCross)!\n");
    }
    /* ----- internally we always use Aplus, Across */
    if ( have_h0 )
      {
	cfg->aPlus = 0.5 * uvar->h0 * ( 1.0 + SQ( uvar->cosi) );
	cfg->aCross = uvar->h0 * uvar->cosi;
      }
    else
      {
	cfg->aPlus = uvar->aPlus;
	cfg->aCross = uvar->aCross;
      }
  }/* check user-input */

  /* ----- prepare SFT-reading ----- */
  if ( XLALUserVarWasSet ( &uvar->IFO ) ) {
    XLAL_CHECK ( (constraints.detector = XLALGetChannelPrefix ( uvar->IFO )) != NULL, XLAL_EFUNC );
  }

  minStartTimeGPS.gpsSeconds = uvar->minStartTime;
  minStartTimeGPS.gpsNanoSeconds = 0;
  maxStartTimeGPS.gpsSeconds = uvar->maxStartTime;
  maxStartTimeGPS.gpsNanoSeconds = 0;
  constraints.minStartTime = &minStartTimeGPS;
  constraints.maxStartTime = &maxStartTimeGPS;

  /* ----- get full SFT-catalog of all matching (multi-IFO) SFTs */
  XLALPrintInfo ( "Finding all SFTs to load ... ");
  XLAL_CHECK ( (catalog = XLALSFTdataFind ( uvar->DataFiles, &constraints )) != NULL, XLAL_EFUNC );
  XLALPrintInfo ( "done. (found %d SFTs)\n", catalog->length );
  if ( constraints.detector ) {
    XLALFree ( constraints.detector );
  }
  XLAL_CHECK ( catalog->length > 0, XLAL_EINVAL, "No matching SFTs for pattern '%s'!\n", uvar->DataFiles );

  /* ----- deduce start- and end-time of the observation spanned by the data */
  GV.numSFTs = catalog->length;	/* total number of SFTs */
  REAL8 Tsft = 1.0 / catalog->data[0].header.deltaF;
  LIGOTimeGPS startTime = catalog->data[0].header.epoch;
  LIGOTimeGPS endTime   = catalog->data[GV.numSFTs-1].header.epoch;
  XLALGPSAdd ( &endTime, Tsft );
  REAL8 duration = GPS2REAL8(endTime) - GPS2REAL8 (startTime);

  /* ----- load ephemeris-data ----- */
  XLAL_CHECK ( (edat = XLALInitBarycenter( uvar->ephemEarth, uvar->ephemSun )) != NULL, XLAL_EFUNC );

  MultiSFTCatalogView *multiCatalogView;
  XLAL_CHECK ( (multiCatalogView = XLALGetMultiSFTCatalogView ( catalog )) != NULL, XLAL_EFUNC );

  UINT4 numDetectors = multiCatalogView->length;
  // ----- get the (multi-IFO) 'detector-state series' for given catalog
  MultiLIGOTimeGPSVector *mTS;
  XLAL_CHECK ( (mTS = XLALTimestampsFromMultiSFTCatalogView ( multiCatalogView )) != NULL, XLAL_EFUNC );

  MultiLALDetector XLAL_INIT_DECL(multiIFO);
  XLAL_CHECK ( XLALMultiLALDetectorFromMultiSFTCatalogView ( &multiIFO, multiCatalogView ) == XLAL_SUCCESS, XLAL_EFUNC );

  REAL8 tOffset = 0.5 * Tsft;
  XLAL_CHECK ( ( multiDetStates = XLALGetMultiDetectorStates( mTS, &multiIFO, edat, tOffset )) != NULL, XLAL_EFUNC );

  // ----- compute or estimate multiNoiseWeights ----------
  MultiNoiseWeights *multiNoiseWeights = NULL;

  // noise-sqrtSX provided by user ==> don't need to load the SFTs at all
  if ( uvar->SignalOnly || (uvar->assumeSqrtSX != NULL) )
    {
      XLAL_CHECK ( !uvar->SignalOnly || (uvar->assumeSqrtSX == NULL), XLAL_EINVAL, "Cannot pass --SignalOnly AND --assumeSqrtSX at the same time!\n");
      LALStringVector *assumeSqrtSX_input;
      if ( uvar->SignalOnly ) {
        assumeSqrtSX_input = XLALCreateStringVector ( "1", NULL );
      } else {
        assumeSqrtSX_input = uvar->assumeSqrtSX;
      }

      MultiNoiseFloor XLAL_INIT_DECL(assumeSqrtSX);
      XLAL_CHECK ( XLALParseMultiNoiseFloor ( &assumeSqrtSX, assumeSqrtSX_input, numDetectors ) == XLAL_SUCCESS, XLAL_EFUNC );

      if ( uvar->SignalOnly ) {
        XLALDestroyStringVector ( assumeSqrtSX_input );
      }

      XLAL_CHECK ( (multiNoiseWeights = XLALCalloc(1,sizeof(*multiNoiseWeights))) != NULL, XLAL_ENOMEM );
      XLAL_CHECK ( (multiNoiseWeights->data = XLALCalloc ( numDetectors, sizeof(*multiNoiseWeights->data) )) != NULL, XLAL_ENOMEM );
      multiNoiseWeights->length = numDetectors;

      REAL8 Sinv = 0;
      UINT4 numSFTs = 0;
      for ( UINT4 X = 0; X < numDetectors; X ++ )
        {
          UINT4 numSFTsX = mTS->data[X]->length;
          numSFTs += numSFTsX;
          XLAL_CHECK ( (multiNoiseWeights->data[X] = XLALCreateREAL8Vector ( numSFTsX )) != NULL, XLAL_EFUNC );

          REAL8 SXinv = 1.0  / ( SQ(assumeSqrtSX.sqrtSn[X]) );
          Sinv += numSFTsX * SXinv;
          for ( UINT4 j = 0; j < numSFTsX; j ++ ) {
            multiNoiseWeights->data[X]->data[j] = SXinv;
          } // for j < numSFTsX
        }  // for X < numDetectors

      // ----- now properly normalize this
      Sinv /= 1.0 * numSFTs;
      for ( UINT4 X = 0; X < numDetectors; X ++ )
        {
          UINT4 numSFTsX = mTS->data[X]->length;
          for ( UINT4 j = 0; j < numSFTsX; j ++ ) {
            multiNoiseWeights->data[X]->data[j] /= Sinv;
          } // for j < numSFTsX
        } // for X < numDetectors

      multiNoiseWeights->Sinv_Tsft = Sinv * Tsft;

    } // if SignalOnly OR assumeSqrtSX given
  else
    {// load the multi-IFO SFT-vectors for noise estimation
      UINT4 wings = uvar->RngMedWindow/2 + 10;   /* extra frequency-bins needed for rngmed */
      REAL8 fMax = uvar->Freq + 1.0 * wings / Tsft;
      REAL8 fMin = uvar->Freq - 1.0 * wings / Tsft;

      MultiSFTVector *multiSFTs;
      XLAL_CHECK ( (multiSFTs = XLALLoadMultiSFTsFromView ( multiCatalogView, fMin, fMax )) != NULL, XLAL_EFUNC );

      MultiPSDVector *multiRngmed = NULL;
      XLAL_CHECK ( (multiRngmed = XLALNormalizeMultiSFTVect ( multiSFTs, uvar->RngMedWindow, NULL )) != NULL, XLAL_EFUNC );
      XLALDestroyMultiSFTVector ( multiSFTs );

      XLAL_CHECK ( (multiNoiseWeights = XLALComputeMultiNoiseWeights ( multiRngmed, uvar->RngMedWindow, 0 )) != NULL, XLAL_EFUNC );
      XLALDestroyMultiPSDVector ( multiRngmed );
    } // if noise-estimation done from SFTs

  /* ----- handle transient-signal window if given ----- */
  if ( XLALUserVarWasSet ( &uvar->transientWindowType ) && strcmp ( uvar->transientWindowType, "none") )
    {
      transientWindow_t transientWindow;	/**< properties of transient-signal window */

      if ( !strcmp ( uvar->transientWindowType, "rect" ) ) {
        transientWindow.type = TRANSIENT_RECTANGULAR;		/* rectangular window [t0, t0+tau] */
      }
      else if ( !strcmp ( uvar->transientWindowType, "exp" ) ) {
        transientWindow.type = TRANSIENT_EXPONENTIAL;		/* exponential decay window e^[-(t-t0)/tau for t>t0, 0 otherwise */
      }
      else {
        XLAL_ERROR ( XLAL_EINVAL, "Illegal transient window '%s' specified: valid are 'none', 'rect' or 'exp'\n", uvar->transientWindowType);
      }

      if ( XLALUserVarWasSet ( &uvar->transientStartTime ) ) {
        transientWindow.t0 = uvar->transientStartTime;
      }
      else {
        transientWindow.t0 = XLALGPSGetREAL8( &startTime ); /* if not set, default window startTime == startTime here */
      }

      transientWindow.tau  = uvar->transientTauDays;

      XLAL_CHECK ( XLALApplyTransientWindow2NoiseWeights ( multiNoiseWeights, mTS, transientWindow ) == XLAL_SUCCESS, XLAL_EFUNC );

    } /* apply transient window to noise-weights */


  /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
  skypos.longitude = uvar->Alpha;
  skypos.latitude = uvar->Delta;
  skypos.system = COORDINATESYSTEM_EQUATORIAL;
  XLALNormalizeSkyPosition ( &skypos.longitude, &skypos.latitude);

  XLAL_CHECK ( (multiAMcoef = XLALComputeMultiAMCoeffs ( multiDetStates, multiNoiseWeights, skypos )) != NULL, XLAL_EFUNC );

  /* OK: we only need the antenna-pattern matrix M_mu_nu */
  cfg->Mmunu = multiAMcoef->Mmunu;

  /* ----- produce a log-string describing the data-specific setup ----- */
  {
    struct tm utc;
    time_t tp;
    CHAR dateStr[512], line[512], summary[1024];
    tp = time(NULL);
    sprintf (summary, "%%%% Date: %s", asctime( gmtime( &tp ) ) );
    strcat (summary, "%% Loaded SFTs: [ " );
    for ( UINT4 X = 0; X < numDetectors; X ++ ) {
      sprintf (line, "%s:%d%s",  multiIFO.sites[X].frDetector.name, mTS->data[X]->length, (X < numDetectors - 1)?", ":" ]\n");
      strcat ( summary, line );
    }
    utc = *XLALGPSToUTC( &utc, (INT4)GPS2REAL8(startTime) );
    strcpy ( dateStr, asctime(&utc) );
    dateStr[ strlen(dateStr) - 1 ] = 0;
    sprintf (line, "%%%% Start GPS time tStart = %12.3f    (%s GMT)\n", GPS2REAL8(startTime), dateStr);
    strcat ( summary, line );
    sprintf (line, "%%%% Total time spanned    = %12.3f s  (%.1f hours)\n", duration, duration/3600 );
    strcat ( summary, line );

    XLAL_CHECK ( (cfg->dataSummary = LALCalloc(1, strlen(summary) + 1 )) != NULL, XLAL_ENOMEM );
    strcpy ( cfg->dataSummary, summary );

    LogPrintfVerbatim( LOG_DEBUG, "%s", cfg->dataSummary );
  } /* write dataSummary string */

  /* free everything not needed any more */
  XLALDestroyMultiTimestamps ( mTS );
  XLALDestroySFTCatalog ( catalog );
  XLALDestroyMultiSFTCatalogView ( multiCatalogView );
  XLALDestroyMultiNoiseWeights ( multiNoiseWeights );
  XLALDestroyMultiDetectorStateSeries ( multiDetStates );
  XLALDestroyMultiAMCoeffs ( multiAMcoef );

  /* Free ephemeris data */
  XLALDestroyEphemerisData (edat);

  return XLAL_SUCCESS;

} /* InitPFS() */
