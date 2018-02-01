/*
 * Copyright (C) 2017 Maximillian Bensch, Reinhard Prix
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
#include <lal/FstatisticTools.h>
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
  PulsarAmplitudeParams pap;    /**< PulsarAmplitudeParameter {h0, cosi, psi, phi0} */
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

  BOOLEAN PureSignal;   /**< If true, calculate 2F for pure signal, i.e. E[2F] = 2F = rho^2 */
  LALStringVector* assumeSqrtSX;/**< Assume stationary Gaussian noise with detector noise-floors sqrt{SX}" */

  CHAR *ephemEarth;	/**< Earth ephemeris file to use */
  CHAR *ephemSun;	/**< Sun ephemeris file to use */

  CHAR *DataFiles;	/**< SFT input-files to use to determine startTime, duration, IFOs and for noise-floor estimation */
  CHAR *outputFstat;	/**< output file to write F-stat estimation results into */
  BOOLEAN printFstat;	/**< print F-stat estimation results to terminal? */
  LIGOTimeGPS minStartTime;	/**< Limit duration to this earliest GPS SFT start-time */
  LIGOTimeGPS maxStartTime;	/**< Limit duration to this latest GPS SFT start-time */

  LALStringVector *timestampsFiles; /**< Names of timestamps files, one per detector */
  LIGOTimeGPS startTime;	/**< Start-time in detector-frame (GPS seconds) */
  INT4 duration;		/**< Duration in seconds */
  LALStringVector* IFOs;	/**< list of detector-names "H1,H2,L1,.." */
  REAL8 Tsft;		        /**< SFT time baseline Tsft */

  CHAR *transientWindowType;	/**< name of transient window ('rect', 'exp',...) */
  LIGOTimeGPS transientStartTime;	/**< GPS start-time of transient window */
  REAL8 transientTau;	        /**< time-scale in seconds of transient window */
  REAL8 transientTauDays;       /**< DEFUNCT */


  BOOLEAN SignalOnly;	/**< DEPRECATED: ALTERNATIVE --assumeSqrtSX and/or --PureSignal */
} UserInput_t;

/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);

int initUserVars ( UserInput_t *uvar );
int InitPFS ( ConfigVariables *cfg, UserInput_t *uvar );

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
  XLAL_CHECK( XLALUserVarReadAllInput( &should_exit, argc, argv, lalAppsVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit ) {
    exit (1);
  }

  XLAL_CHECK_MAIN ( (VCSInfoString = XLALGetVersionString(0)) != NULL, XLAL_EFUNC );

  /* Initialize code-setup */
  XLAL_CHECK_MAIN ( InitPFS ( &GV, &uvar ) == XLAL_SUCCESS, XLAL_EFUNC );

  rho2 = XLALComputeOptimalSNR2FromMmunu ( GV.pap, GV.Mmunu );
  XLAL_CHECK_MAIN ( xlalErrno == XLAL_SUCCESS, XLAL_EFUNC );

  /* F-statistic expected mean and standard deviation */
  const REAL8 twoF_expected = uvar.PureSignal ? ( rho2 ) : ( 4.0 + rho2 );
  const REAL8 twoF_sigma    = uvar.PureSignal ? (    0 ) : ( sqrt( 8.0 + 4.0 * rho2 ) );

  /* output predicted Fstat-value, if requested */
  if (uvar.printFstat) {
    fprintf(stdout, "\n%.1f\n", twoF_expected);
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
      fprintf (fpFstat, "twoF_expected = %g;\n", twoF_expected);
      fprintf (fpFstat, "twoF_sigma    = %g;\n", twoF_sigma);

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

  uvar->minStartTime.gpsSeconds = 0;
  uvar->maxStartTime.gpsSeconds = LAL_INT4_MAX;

  uvar->PureSignal = 0;

  uvar->assumeSqrtSX = NULL;
  uvar->SignalOnly = 0;

  uvar->phi0 = 0;
  uvar->transientWindowType = XLALStringDuplicate ( "none" );

  uvar->Tsft=1800;

  /* register all our user-variables */
  lalUserVarHelpOptionSubsection = "Output control";
  XLALRegisterUvarMember( outputFstat,   STRING,      0,  NODEFAULT,"Output-file (octave/matlab) for predicted F-stat value, variance and antenna-patterns" );
  XLALRegisterUvarMember( printFstat,    BOOLEAN,     0,  OPTIONAL, "Print predicted F-stat value to terminal" );
  XLALRegisterUvarMember( PureSignal,    BOOLEAN,    'P', OPTIONAL, "If true return 2F=SNR^2 ('pure signal without noise'). Otherwise return E[2F] = 4 + SNR^2.");

  lalUserVarHelpOptionSubsection = "Signal parameters";
  XLALRegisterUvarMember( h0,            REAL8,      's', NODEFAULT,"Overall GW amplitude h0 (alternatively use " UVAR_STR2AND(aPlus, aCross) ")." );
  XLALRegisterUvarMember( cosi,          REAL8,      'i', NODEFAULT,"Inclination angle of rotation axis cos(iota) (alternatively use " UVAR_STR2AND(aPlus, aCross) ").");
  XLALRegisterUvarMember( aPlus,         REAL8,       0,  NODEFAULT,"'Plus' polarization amplitude A_+ (alternative to " UVAR_STR2AND(h0, cosi) ").");
  XLALRegisterUvarMember( aCross,        REAL8,       0,  NODEFAULT,"'Cross' polarization amplitude: A_x (alternative to " UVAR_STR2AND(h0, cosi) ").");

  XLALRegisterUvarMember( psi,           REAL8,       0,  REQUIRED, "Polarisation angle in radians");
  XLALRegisterUvarMember( phi0,          REAL8,       0,  OPTIONAL, "Initial GW phase phi0 in radians");

  XLALRegisterUvarMember( Alpha,         RAJ,        'a', REQUIRED, "Sky position: equatorial J2000 right ascension");
  XLALRegisterUvarMember( Delta,         DECJ,       'd', REQUIRED, "Sky position: equatorial J2000 right declination");

  lalUserVarHelpOptionSubsection = "Data and noise properties";
  XLALRegisterUvarMember( DataFiles,     STRING,     'D', NODEFAULT,"Per-detector SFTs (for detectors, timestamps and noise-estimate)\n"
                          "(Alternatives: " UVAR_STR(assumeSqrtSX)", "UVAR_STR(timestampsFiles)", "UVAR_STR(IFOs) ", "UVAR_STR2AND(minStartTime,maxStartTime)").");
  XLALRegisterUvarMember( Freq,          REAL8,      'F', NODEFAULT,"Frequency for noise-floor estimation (required if not given " UVAR_STR(assumeSqrtSX) ").");
  XLALRegisterUvarMember( RngMedWindow,  INT4,       'k', OPTIONAL, "Running median size for noise-floor estimation (only used if not given " UVAR_STR(assumeSqrtSX) ").");

  XLALRegisterUvarMember( assumeSqrtSX,  STRINGVector,0,  OPTIONAL, "Assume stationary per-detector noise-floor sqrt(S[X]) instead of estimating "
                          "(required if not given " UVAR_STR(DataFiles)").");

  XLALRegisterUvarMember( IFOs,          STRINGVector,0,  NODEFAULT,"CSV list of detectors, eg. \"H1,L1,...\" (required if not given " UVAR_STR(DataFiles)").");
  XLALRegisterUvarMember(timestampsFiles,STRINGVector,0,  NODEFAULT,"CSV list of SFT timestamps files, one per detector (conflicts with " UVAR_STR(DataFiles) ").");
  XLALRegisterUvarMember( minStartTime,  EPOCH,       0,  OPTIONAL, "Limit duration to [" UVAR_STR(minStartTime) ", " UVAR_STR(maxStartTime) " + " UVAR_STR(Tsft) ").");
  XLALRegisterUvarMember( maxStartTime,  EPOCH,       0,  OPTIONAL, "Limit duration to [" UVAR_STR(minStartTime) ", " UVAR_STR(maxStartTime) " + " UVAR_STR(Tsft) ").");

  XLALRegisterUvarMember( Tsft,          REAL8,       0,  OPTIONAL, "Time baseline of SFTs in seconds (conflicts with " UVAR_STR(DataFiles) ")." );

  lalUserVarHelpOptionSubsection = "Transient signal properties";
  XLALRegisterUvarMember( transientWindowType,STRING, 0,  OPTIONAL, "Transient-signal window function to assume. ('none', 'rect', 'exp').");
  XLALRegisterUvarMember( transientStartTime, EPOCH,  0,  NODEFAULT,"GPS start-time 't0' of transient signal window.");
  XLALRegisterUvarMember( transientTau,       REAL8,  0,  NODEFAULT,"Timescale 'tau' of transient signal window in seconds.");

  // ---------- developer options ----------
  lalUserVarHelpOptionSubsection = "";
  XLALRegisterUvarMember( ephemEarth,    STRING,      0,  DEVELOPER, "Earth ephemeris file to use");
  XLALRegisterUvarMember( ephemSun,      STRING,      0,  DEVELOPER, "Sun ephemeris file to use");

  // ---------- deprecated options ---------
  XLALRegisterUvarMember( SignalOnly,    BOOLEAN,    'S', DEPRECATED, "use --assumeSqrtSX and/or --PureSignal instead");
  XLALRegisterUvarMember( startTime,     EPOCH,       0,  DEPRECATED, "Start-time of the signal in detector-frame");
  XLALRegisterUvarMember( duration,      INT4,        0,  DEPRECATED, "Duration of the signal in seconds");

  // ---------- defunct options ---------
  XLALRegisterUvarMember( transientTauDays,REAL8,  0,  DEFUNCT, "DEFUNCT: use " UVAR_STR(transientTau) " instead.");

  return XLAL_SUCCESS;

} /* initUserVars() */

/** Initialized Fstat-code: handle user-input and set everything up. */
int
InitPFS ( ConfigVariables *cfg, UserInput_t *uvar )
{
  XLAL_CHECK ( (cfg != NULL) && (uvar != NULL), XLAL_EINVAL );

  SFTCatalog *catalog = NULL;
  SkyPosition skypos;

  EphemerisData *edat = NULL;		    	/* ephemeris data */
  MultiAMCoeffs *multiAMcoef = NULL;
  MultiDetectorStateSeries *multiDetStates = NULL; /* pos, vel and LMSTs for detector at times t_i */

  { /* Check user-input consistency */
    BOOLEAN have_h0   = XLALUserVarWasSet ( &uvar->h0 );
    BOOLEAN have_cosi = XLALUserVarWasSet ( &uvar->cosi );
    BOOLEAN have_Ap   = XLALUserVarWasSet ( &uvar->aPlus );
    BOOLEAN have_Ac   = XLALUserVarWasSet ( &uvar->aCross );

    /* ----- handle {h0,cosi} || {aPlus,aCross} freedom ----- */
    XLAL_CHECK( (( have_h0 && !have_cosi ) || ( !have_h0 && have_cosi )) == 0, XLAL_EINVAL, "Need both (h0, cosi) to specify signal!\n");
    XLAL_CHECK(( ( have_Ap && !have_Ac) || ( !have_Ap && have_Ac ) ) == 0, XLAL_EINVAL, "Need both (aPlus, aCross) to specify signal!\n");
    XLAL_CHECK(( have_h0 && have_Ap ) == 0, XLAL_EINVAL, "Overdetermined: specify EITHER (h0,cosi) OR (aPlus,aCross)!\n");
    /* ----- internally we always use h0, cosi */
    if ( have_h0 )
      {
        cfg->pap.h0=uvar->h0;
        cfg->pap.cosi=uvar->cosi;
      }
    else
      {
        cfg->pap.h0 = uvar->aPlus + sqrt( SQ( uvar->aPlus ) - SQ( uvar->aCross ) );
        cfg->pap.cosi= uvar->aCross / cfg->pap.h0;
      }
    cfg->pap.psi=uvar->psi;
    cfg->pap.phi0=uvar->phi0;
  }/* check user-input */

  MultiLALDetector XLAL_INIT_DECL(multiIFO);
  // ----- IFOs : only from one of {--IFOs, --DataFiles }: mutually exclusive
  BOOLEAN have_IFOs      = (uvar->IFOs != NULL);
  BOOLEAN have_SFTs = (uvar->DataFiles != NULL);
  XLAL_CHECK ( have_IFOs || have_SFTs, XLAL_EINVAL, "Need one of --IFOs, --DataFiles to determine detectors\n");
  if ( have_IFOs )
    {
      XLAL_CHECK ( XLALParseMultiLALDetector ( &multiIFO, uvar->IFOs ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

  // ----- TIMESTAMPS: either from --timestampsFiles, --startTime+duration, or --SFTs
  BOOLEAN have_startTime = XLALUserVarWasSet ( &uvar->startTime );
  BOOLEAN have_duration = XLALUserVarWasSet ( &uvar->duration );
  BOOLEAN have_timestampsFiles = ( uvar->timestampsFiles != NULL );
  BOOLEAN have_assumeSqrtSX = ( uvar->assumeSqrtSX != NULL );
  // need BOTH startTime+duration or none
  XLAL_CHECK ( ( have_duration && have_startTime) || !( have_duration || have_startTime ), XLAL_EINVAL, "Need BOTH {--startTime,--duration} or NONE\n");
  // at least one of {startTime,timestamps,SFTs} required
  XLAL_CHECK ( have_timestampsFiles || have_startTime || have_SFTs , XLAL_EINVAL, "Need at least one of {--timestampsFiles, --startTime+duration, --DataFiles}\n" );
  // don't allow timestamps + {startTime+duration OR SFTs}
  XLAL_CHECK ( !have_timestampsFiles || !(have_startTime||have_SFTs), XLAL_EINVAL, "--timestampsFiles incompatible with {--DataFiles or --startTime+duration}\n");
  // if we have timestamps or startTime+duration we also need assumeSqrtSX
  XLAL_CHECK ( have_SFTs || ( (have_timestampsFiles || have_startTime) && have_assumeSqrtSX ), XLAL_EINVAL, "Need --assumeSqrtSX \n");

  MultiLIGOTimeGPSVector *mTS = NULL;
  REAL8 Tsft=0, duration=0;
  LIGOTimeGPS startTime = {0,0};
  UINT4 numDetectors=0;
  MultiSFTCatalogView *multiCatalogView = NULL;
  // ----- compute or estimate multiTimestamps ----------
  if ( have_SFTs )
    {
      SFTConstraints XLAL_INIT_DECL(constraints);
      /* ----- prepare SFT-reading ----- */
      constraints.minStartTime = &uvar->minStartTime;
      constraints.maxStartTime = &uvar->maxStartTime;

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
      Tsft = 1.0 / catalog->data[0].header.deltaF;
      startTime = catalog->data[0].header.epoch;
      LIGOTimeGPS endTime   = catalog->data[GV.numSFTs-1].header.epoch;
      XLALGPSAdd ( &endTime, Tsft );
      duration = GPS2REAL8(endTime) - GPS2REAL8 (startTime);

      XLAL_CHECK ( (multiCatalogView = XLALGetMultiSFTCatalogView ( catalog )) != NULL, XLAL_EFUNC );

      numDetectors = multiCatalogView->length;
      // ----- get the (multi-IFO) 'detector-state series' for given catalog
      XLAL_CHECK ( (mTS = XLALTimestampsFromMultiSFTCatalogView ( multiCatalogView )) != NULL, XLAL_EFUNC );

      XLAL_CHECK ( XLALMultiLALDetectorFromMultiSFTCatalogView ( &multiIFO, multiCatalogView ) == XLAL_SUCCESS, XLAL_EFUNC );
    } // endif have_SFTs
  else
    {
      Tsft=uvar->Tsft;
      numDetectors=multiIFO.length;
      if ( have_timestampsFiles )
        {
          XLAL_CHECK ( (mTS = XLALReadMultiTimestampsFiles ( uvar->timestampsFiles )) != NULL, XLAL_EFUNC );
          XLAL_CHECK ( (mTS->length > 0) && (mTS->data != NULL), XLAL_EINVAL, "Got empty timestamps-list from XLALReadMultiTimestampsFiles()\n" );
          for ( UINT4 X=0; X < mTS->length; X ++ )
            {
              mTS->data[X]->deltaT = uvar->Tsft;	// Tsft information not given by timestamps-file
            }
          duration = mTS->data[0]->length * mTS->data[0]->deltaT;
          startTime=mTS->data[0]->data[0];
        } // endif have_timestampsFiles
      else if ( have_startTime && have_duration )
        {
          XLAL_CHECK ( ( mTS = XLALMakeMultiTimestamps ( uvar->startTime, uvar->duration, uvar->Tsft, 0, multiIFO.length )) != NULL, XLAL_EFUNC );
          duration=uvar->duration;
          startTime=uvar->startTime;
        } // endif have_startTime
    }// calculate or estimate timestamps
  /* ----- load ephemeris-data ----- */
  XLAL_CHECK ( (edat = XLALInitBarycenter( uvar->ephemEarth, uvar->ephemSun )) != NULL, XLAL_EFUNC );

  // ----- compute or estimate multiDetectorStates ----------
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
      GV.numSFTs=numSFTs;

    } // if SignalOnly OR assumeSqrtSX given
  else if(have_SFTs)
    {// load the multi-IFO SFT-vectors for noise estimation
      XLAL_CHECK ( XLALUserVarWasSet ( &uvar->Freq ), XLAL_EINVAL, "Need a frequency for noise-estimation in SFTs!\n");
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

      int twtype;
      XLAL_CHECK ( (twtype = XLALParseTransientWindowName ( uvar->transientWindowType )) >= 0, XLAL_EFUNC );
      transientWindow.type = twtype;

      if ( XLALUserVarWasSet ( &uvar->transientStartTime ) ) {
        transientWindow.t0 = uvar->transientStartTime.gpsSeconds; // dropping ns part
      } else {
        XLAL_ERROR ( XLAL_EINVAL, "Required input " UVAR_STR(transientStartTime) " missing for transient window type '%s'", uvar->transientWindowType );
      }

      if ( XLALUserVarWasSet ( &uvar->transientTau ) ) {
        transientWindow.tau  = uvar->transientTau;
      } else {
        XLAL_ERROR ( XLAL_EINVAL, "Required input " UVAR_STR(transientStartTime) " missing for transient window type '%s'", uvar->transientWindowType );
      }

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
