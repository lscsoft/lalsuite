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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

/*********************************************************************************/
/**
 * \author I. Gholami, R. Prix
 * \file
 * \ingroup lalpulsar_bin_Fstatistic
 * \brief
 * Calculate the *expected* (multi-IFO) F-statistic for pulsar GW signals, without actually
 * performing a search. The "F-statistic" was introduced in \cite JKS98 and Cutler-Schutz 2005.
 * Contrary to SemiAnalyticF this code can use (multi-IFO) SFTs to specify the startTime,
 * duration, detectors and noise-floors to use in the estimation.
 *
 */
#include "config.h"

#include <stdio.h>

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
#include <lal/LALPulsarVCSInfo.h>

/* local includes */

/*---------- DEFINES ----------*/
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
extern int vrbflg;		/**< defined in lal/lib/std/LALError.c */

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
  LIGOTimeGPS minStartTime;	/**< SFT timestamps (from DataFiles or generated internally) must be >= this GPS timestamp */
  LIGOTimeGPS maxStartTime;	/**< SFT timestamps (from DataFiles) must be < this GPS timestamp */
  REAL8 duration;	/**< SFT timestamps (generated internally) will be < minStartTime+duration */

  LALStringVector *timestampsFiles; /**< Names of timestamps files, one per detector */
  LALStringVector* IFOs;	/**< list of detector-names "H1,H2,L1,.." */
  REAL8 Tsft;		        /**< SFT time baseline Tsft */

  CHAR *transientWindowType;	/**< name of transient window ('rect', 'exp',...) */
  LIGOTimeGPS transientStartTime;	/**< GPS start-time of transient window */
  REAL8 transientTau;	        /**< time-scale in seconds of transient window */
  REAL8 transientTauDays;       /**< DEFUNCT */


  BOOLEAN SignalOnly;	/**< DEFUNCT: use --assumeSqrtSX */
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
  CHAR *VCSInfoString;          /* Git version string */

  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register all user-variable */
  XLAL_CHECK_MAIN ( initUserVars( &uvar) == XLAL_SUCCESS, XLAL_EFUNC );

  /* do ALL cmdline and cfgfile handling */
  BOOLEAN should_exit = 0;
  XLAL_CHECK( XLALUserVarReadAllInput( &should_exit, argc, argv, lalPulsarVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit ) {
    exit (1);
  }

  XLAL_CHECK_MAIN ( (VCSInfoString = XLALVCSInfoString(lalPulsarVCSInfoList, 0, "%% ")) != NULL, XLAL_EFUNC );

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
  uvar->RngMedWindow = FstatOptionalArgsDefaults.runningMedianWindow;

  uvar->ephemEarth = XLALStringDuplicate("earth00-40-DE405.dat.gz");
  uvar->ephemSun = XLALStringDuplicate("sun00-40-DE405.dat.gz");

  uvar->outputFstat = NULL;
  uvar->printFstat = 1;

  uvar->minStartTime.gpsSeconds = 0;
  uvar->maxStartTime.gpsSeconds = LAL_INT4_MAX;
  uvar->duration = 0;

  uvar->PureSignal = 0;

  uvar->assumeSqrtSX = NULL;

  uvar->phi0 = 0;
  uvar->transientWindowType = XLALStringDuplicate ( "none" );

  uvar->Tsft=1800;

  /* register all our user-variables */
  lalUserVarHelpOptionSubsection = "Output control";
  XLALRegisterUvarMember( outputFstat,   STRING,      0,  NODEFAULT,"Output-file (octave/matlab) for predicted F-stat value, variance and antenna-patterns" );
  XLALRegisterUvarMember( printFstat,    BOOLEAN,     0,  OPTIONAL, "Print predicted F-stat value to terminal" );
  XLALRegisterUvarMember( PureSignal,    BOOLEAN,    'P', OPTIONAL, "If true return 2F=SNR^2 ('pure signal without noise'). Otherwise return E[2F] = 4 + SNR^2.");

  lalUserVarHelpOptionSubsection = "Signal parameters";
  XLALRegisterUvarMember( h0,            REAL8,      's', NODEFAULT, "Overall GW amplitude h0 (required if " UVAR_STR(cosi) " used, alternatively use " UVAR_STR2AND(aPlus, aCross) ")." );
  XLALRegisterUvarMember( cosi,          REAL8,      'i', NODEFAULT, "Inclination angle of rotation axis cos(iota) (required if " UVAR_STR(h0) " used, alternatively use " UVAR_STR2AND(aPlus, aCross) ").");
  XLALRegisterUvarMember( aPlus,         REAL8,       0,  NODEFAULT, "'Plus' polarization amplitude A_+  (required if " UVAR_STR(aCross) " used, alternative to " UVAR_STR2AND(h0, cosi) ").");
  XLALRegisterUvarMember( aCross,        REAL8,       0,  NODEFAULT, "'Cross' polarization amplitude: A_x  (required if " UVAR_STR(aPlus) " used, alternative to " UVAR_STR2AND(h0, cosi) ").");

  XLALRegisterUvarMember( psi,           REAL8,       0,  OPTIONAL, "Polarisation angle in radians (required if " UVAR_STR4OR(h0, cosi, aPlus, aCross) " used)");
  XLALRegisterUvarMember( phi0,          REAL8,       0,  OPTIONAL, "Initial GW phase phi0 in radians");

  XLALRegisterUvarMember( Alpha,         RAJ,        'a', OPTIONAL, "Sky position: equatorial J2000 right ascension (required if " UVAR_STR4OR(h0, cosi, aPlus, aCross) " used)");
  XLALRegisterUvarMember( Delta,         DECJ,       'd', OPTIONAL, "Sky position: equatorial J2000 right declination (required if " UVAR_STR4OR(h0, cosi, aPlus, aCross) " used)");

  lalUserVarHelpOptionSubsection = "Data and noise properties";
  XLALRegisterUvarMember( DataFiles,     STRING,     'D', NODEFAULT,"Per-detector SFTs (for detectors, timestamps and noise-estimate). Possibilities are:\n"
                          " - '<SFT file>;<SFT file>;...', where <SFT file> may contain wildcards\n - 'list:<file containing list of SFT files>'\n"
                          "(Alternatives: " UVAR_STR2AND(assumeSqrtSX,IFOs)" and one of "UVAR_STR(timestampsFiles)" or "UVAR_STR2AND(minStartTime,duration)").");
  XLALRegisterUvarMember( Freq,          REAL8,      'F', NODEFAULT,"Frequency for noise-floor estimation (required if not given " UVAR_STR(assumeSqrtSX) ").");
  XLALRegisterUvarMember( RngMedWindow,  INT4,       'k', OPTIONAL, "Running median size for noise-floor estimation (only used if not given " UVAR_STR(assumeSqrtSX) ").");

  XLALRegisterUvarMember( assumeSqrtSX,  STRINGVector,0,  OPTIONAL, "Assume stationary per-detector noise-floor sqrt(S[X]) instead of estimating "
                          "(required if not given " UVAR_STR(DataFiles)").");

  XLALRegisterUvarMember( IFOs,          STRINGVector,0,  NODEFAULT,"CSV list of detectors, eg. \"H1,L1,...\" (required if not given " UVAR_STR(DataFiles)").");
  XLALRegisterUvarMember(timestampsFiles,STRINGVector,0,  NODEFAULT,"CSV list of SFT timestamps files, one per detector (conflicts with " UVAR_STR(DataFiles) ").");
  XLALRegisterUvarMember( minStartTime,  EPOCH,       0,  OPTIONAL, "SFT timestamps must be >= this GPS timestamp.");
  XLALRegisterUvarMember( maxStartTime,  EPOCH,       0,  OPTIONAL, "SFT timestamps must be < this GPS timestamp (only valid with " UVAR_STR(DataFiles) ").");
  XLALRegisterUvarMember( duration,      REAL8,       0,  OPTIONAL, "SFT timestamps will be < " UVAR_STR(minStartTime) "+" UVAR_STR(duration) " (only valid without " UVAR_STR(DataFiles) ").");

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

  // ---------- defunct options ---------
  XLALRegisterUvarMember( transientTauDays,REAL8,     0,  DEFUNCT, "use " UVAR_STR(transientTau) " instead.");
  XLALRegisterUvarMember( SignalOnly,    BOOLEAN,    'S', DEFUNCT, "use --assumeSqrtSX instead");

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
    BOOLEAN have_psi   = XLALUserVarWasSet ( &uvar->psi );
    BOOLEAN have_Alpha = XLALUserVarWasSet ( &uvar->Alpha );
    BOOLEAN have_Delta = XLALUserVarWasSet ( &uvar->Delta );

    /* ----- handle {h0,cosi} || {aPlus,aCross} freedom ----- */
    XLAL_CHECK( (( have_h0 && !have_cosi ) || ( !have_h0 && have_cosi )) == 0, XLAL_EINVAL, "Need both (h0, cosi) to specify signal!\n");
    XLAL_CHECK(( ( have_Ap && !have_Ac) || ( !have_Ap && have_Ac ) ) == 0, XLAL_EINVAL, "Need both (aPlus, aCross) to specify signal!\n");
    XLAL_CHECK(( have_h0 && have_Ap ) == 0, XLAL_EINVAL, "Overdetermined: specify EITHER (h0,cosi) OR (aPlus,aCross)!\n");
    XLAL_CHECK(( ( have_h0 || have_Ap ) && !( have_psi && have_Alpha && have_Delta ) ) == 0, XLAL_EINVAL, "If " UVAR_STR4OR(h0, cosi, aPlus, aCross) ", also need a full set of " UVAR_STR3AND(psi, Alpha, Delta) ".");
    /* ----- internally we always use h0, cosi */
    if ( have_h0 )
      {
        cfg->pap.aPlus = 0.5 * uvar->h0 * (1.0 + SQ(uvar->cosi));
        cfg->pap.aCross = uvar->h0 * uvar->cosi;
      }
    else
      {
        cfg->pap.aPlus = uvar->aPlus;
        cfg->pap.aCross = uvar->aCross;
      }
    cfg->pap.psi=uvar->psi;
    cfg->pap.phi0=uvar->phi0;
  }/* check user-input */

  // ----- the following quantities need to specified via the given user inputs
  REAL8 Tsft;
  UINT4 numDetectors;

  MultiLALDetector XLAL_INIT_DECL(multiIFO);
  MultiLIGOTimeGPSVector *mTS = NULL;
  MultiNoiseWeights *multiNoiseWeights = NULL;

  // ----- IFOs : only from one of {--IFOs, --DataFiles }: mutually exclusive
  BOOLEAN have_IFOs = UVAR_SET(IFOs);
  BOOLEAN have_SFTs = UVAR_SET(DataFiles);
  BOOLEAN have_Tsft = UVAR_SET(Tsft);
  XLAL_CHECK ( (have_IFOs || have_SFTs) && !(have_IFOs && have_SFTs), XLAL_EINVAL, "Need exactly one of " UVAR_STR2OR(IFOs,DataFiles) " to determine detectors\n");
  XLAL_CHECK ( !(have_SFTs && have_Tsft), XLAL_EINVAL, UVAR_STR(Tsft) " cannot be specified with " UVAR_STR(DataFiles) ".");

  // ----- get timestamps from EITHER one of --timestampsFiles, (--minStartTime,--duration) or --SFTs
  BOOLEAN have_timeSpan     = UVAR_SET2(minStartTime,duration);
  BOOLEAN have_maxStartTime = UVAR_SET(maxStartTime);
  BOOLEAN have_duration     = UVAR_SET(duration);
  BOOLEAN have_timestamps   = UVAR_SET(timestampsFiles);
  // at least one of {timestamps,minStartTime+duration,SFTs} required
  XLAL_CHECK ( have_timestamps || have_timeSpan == 2 || have_SFTs, XLAL_EINVAL,
               "Need at least one of {" UVAR_STR(timestampsFiles)", "UVAR_STR2AND(minStartTime,duration)", or "UVAR_STR(DataFiles)"}." );
  // don't allow timestamps AND SFTs
  XLAL_CHECK ( !(have_timestamps && have_SFTs), XLAL_EINVAL, UVAR_STR(timestampsFiles) " is incompatible with " UVAR_STR(DataFiles) ".");
  // don't allow maxStartTime WITHOUT SFTs because the old behaviour in that case
  // was inconsistent and we now require duration instead to avoid ambiguity
  XLAL_CHECK ( !(have_maxStartTime && !have_SFTs), XLAL_EINVAL, "Using " UVAR_STR(maxStartTime) " is incompatible with NOT providing " UVAR_STR(DataFiles) ", please use " UVAR_STR2AND(minStartTime,duration) " instead." );
  // don't allow duration AND SFTs either
  XLAL_CHECK ( !(have_duration && have_SFTs), XLAL_EINVAL, "Using " UVAR_STR(duration) " is incompatible with " UVAR_STR(DataFiles) "; if you want to set constraints on the loaded SFTs, please use " UVAR_STR(maxStartTime) " instead." );

  // if we don't have SFTs, then we need assumeSqrtSX
  BOOLEAN have_assumeSqrtSX = UVAR_SET(assumeSqrtSX);
  BOOLEAN have_Freq         = UVAR_SET(Freq);
  XLAL_CHECK ( have_SFTs || have_assumeSqrtSX, XLAL_EINVAL, "Need at least one of " UVAR_STR2OR(assumeSqrtSX,DataFiles) " for noise-floor.");
  // need --Freq for noise-floor estimation if we don't have --assumeSqrtSX
  XLAL_CHECK ( have_assumeSqrtSX || have_Freq, XLAL_EINVAL, "Need at least one of " UVAR_STR2OR(assumeSqrtSX,Freq) " for noise-floor.");

  // ----- compute or estimate multiTimestamps ----------
  if ( have_SFTs )
    {
      SFTConstraints XLAL_INIT_DECL(constraints);
      MultiSFTCatalogView *multiCatalogView = NULL;
      /* ----- prepare SFT-reading ----- */
      constraints.minStartTime = &uvar->minStartTime;
      constraints.maxStartTime = &uvar->maxStartTime;

      /* ----- get full SFT-catalog of all matching (multi-IFO) SFTs */
      XLALPrintInfo ( "Finding all SFTs to load ... ");
      XLAL_CHECK ( (catalog = XLALSFTdataFind ( uvar->DataFiles, &constraints )) != NULL, XLAL_EFUNC );
      XLALPrintInfo ( "done. (found %d SFTs)\n", catalog->length );
      XLAL_CHECK ( catalog->length > 0, XLAL_EINVAL, "No matching SFTs for pattern '%s' and GPS constraints [%d,%d)!\n", uvar->DataFiles, uvar->minStartTime.gpsSeconds, uvar->maxStartTime.gpsSeconds );

      /* ----- deduce start- and end-time of the observation spanned by the data */
      Tsft = 1.0 / catalog->data[0].header.deltaF;

      XLAL_CHECK ( (multiCatalogView = XLALGetMultiSFTCatalogView ( catalog )) != NULL, XLAL_EFUNC );

      numDetectors = multiCatalogView->length;
      // ----- get the (multi-IFO) 'detector-state series' for given catalog
      XLAL_CHECK ( (mTS = XLALTimestampsFromMultiSFTCatalogView ( multiCatalogView )) != NULL, XLAL_EFUNC );
      XLAL_CHECK( mTS->length == numDetectors, XLAL_EINVAL, "Got %d detectors but %d sets of timestamps.", numDetectors, mTS->length );
      XLAL_CHECK ( XLALMultiLALDetectorFromMultiSFTCatalogView ( &multiIFO, multiCatalogView ) == XLAL_SUCCESS, XLAL_EFUNC );

      // ----- estimate noise-floor from SFTs if --assumeSqrtSX was not given:
      if ( !have_assumeSqrtSX )
        {
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
        }

      XLALDestroyMultiSFTCatalogView ( multiCatalogView );
      XLALDestroySFTCatalog ( catalog );

    } // if have_SFTs
  else
    {
      XLAL_CHECK ( XLALParseMultiLALDetector ( &multiIFO, uvar->IFOs ) == XLAL_SUCCESS, XLAL_EFUNC );
      numDetectors = multiIFO.length;
      Tsft = uvar->Tsft;

      if ( have_timestamps )
        {
          LIGOTimeGPS *endTimeGPS = NULL;
          if ( have_duration ) {
            LIGOTimeGPS tempGPS = uvar->minStartTime;
            XLALGPSAdd( &tempGPS, uvar->duration);
            endTimeGPS = &tempGPS;
          }
          XLAL_CHECK ( (mTS = XLALReadMultiTimestampsFilesConstrained ( uvar->timestampsFiles, &(uvar->minStartTime), endTimeGPS )) != NULL, XLAL_EFUNC );
          XLAL_CHECK ( (mTS->length > 0) && (mTS->data != NULL), XLAL_EINVAL, "Got empty timestamps-list from XLALReadMultiTimestampsFiles()\n" );
          XLAL_CHECK( mTS->length == numDetectors, XLAL_EINVAL, "Got %d detectors but %d sets of timestamps.", numDetectors, mTS->length );
          for ( UINT4 X=0; X < mTS->length; X ++ ) {
            mTS->data[X]->deltaT = Tsft;	// Tsft information not given by timestamps-file
          }
        } // if have_timestamps
      else if ( have_timeSpan == 2 ) // if timespan only
        {
          XLAL_CHECK ( ( mTS = XLALMakeMultiTimestamps ( uvar->minStartTime, uvar->duration, Tsft, 0, numDetectors )) != NULL, XLAL_EFUNC );
        } // have_timeSpan
      else {
        XLAL_ERROR (XLAL_EINVAL, "Something has gone wrong: couldn't deduce timestamps");
      }
    } // if !have_SFTs

  // ---------- determine start-time and total amount of data
  LIGOTimeGPS startTime = {LAL_INT4_MAX,0};
  cfg->numSFTs = 0;
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      if ( XLALGPSCmp ( &(mTS->data[X]->data[0]), &startTime ) < 0 ) {
        startTime = mTS->data[X]->data[0];
      }
      GV.numSFTs += mTS->data[X]->length;
    }
  REAL8 Tdata = GV.numSFTs * Tsft;

  // ---------- compute noise-weights from --assumeSqrtSX instead of from SFTs
  if ( have_assumeSqrtSX )
    {
      MultiNoiseFloor XLAL_INIT_DECL(assumeSqrtSX);
      XLAL_CHECK ( XLALParseMultiNoiseFloor ( &assumeSqrtSX, uvar->assumeSqrtSX, numDetectors ) == XLAL_SUCCESS, XLAL_EFUNC );
      // no need to check assumeSqrtSX.length - XLALParseMultiNoiseFloor() does it internally
      for ( UINT4 X = 0; X < numDetectors; X ++ )
        {
            XLAL_CHECK(assumeSqrtSX.sqrtSn[X] > 0, XLAL_EDOM, "all entries of assumeSqrtSX must be >0");
        }

      XLAL_CHECK ( (multiNoiseWeights = XLALCalloc(1,sizeof(*multiNoiseWeights))) != NULL, XLAL_ENOMEM );
      XLAL_CHECK ( (multiNoiseWeights->data = XLALCalloc ( numDetectors, sizeof(*multiNoiseWeights->data) )) != NULL, XLAL_ENOMEM );
      multiNoiseWeights->length = numDetectors;

      REAL8 Sinv = 0;
      for ( UINT4 X = 0; X < numDetectors; X ++ )
        {
          UINT4 numSFTsX = mTS->data[X]->length;
          XLAL_CHECK ( (multiNoiseWeights->data[X] = XLALCreateREAL8Vector ( numSFTsX )) != NULL, XLAL_EFUNC );

          REAL8 SXinv = 1.0  / ( SQ(assumeSqrtSX.sqrtSn[X]) );
          Sinv += numSFTsX * SXinv;
          for ( UINT4 j = 0; j < numSFTsX; j ++ ) {
            multiNoiseWeights->data[X]->data[j] = SXinv;
          } // for j < numSFTsX
        }  // for X < numDetectors

      // ----- now properly normalize this
      Sinv /= 1.0 * GV.numSFTs;
      for ( UINT4 X = 0; X < numDetectors; X ++ )
        {
          UINT4 numSFTsX = mTS->data[X]->length;
          for ( UINT4 j = 0; j < numSFTsX; j ++ ) {
            multiNoiseWeights->data[X]->data[j] /= Sinv;
          } // for j < numSFTsX
        } // for X < numDetectors

      multiNoiseWeights->Sinv_Tsft = Sinv * Tsft;

    } // if --assumeSqrtSX given

  /* ----- load ephemeris-data ----- */
  XLAL_CHECK ( (edat = XLALInitBarycenter( uvar->ephemEarth, uvar->ephemSun )) != NULL, XLAL_EFUNC );

  // ----- get multiDetectorStates ----------
  REAL8 tOffset = 0.5 * Tsft;
  XLAL_CHECK ( ( multiDetStates = XLALGetMultiDetectorStates( mTS, &multiIFO, edat, tOffset )) != NULL, XLAL_EFUNC );

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

  // ----- compute antenna-pattern matrix M_munu ----------
  skypos.longitude = uvar->Alpha;
  skypos.latitude = uvar->Delta;
  skypos.system = COORDINATESYSTEM_EQUATORIAL;
  XLALNormalizeSkyPosition ( &skypos.longitude, &skypos.latitude);

  XLAL_CHECK ( (multiAMcoef = XLALComputeMultiAMCoeffs ( multiDetStates, multiNoiseWeights, skypos )) != NULL, XLAL_EFUNC );

  cfg->Mmunu = multiAMcoef->Mmunu;
  XLALDestroyMultiAMCoeffs ( multiAMcoef );

  /* ----- produce a log-string describing the data-specific setup ----- */
  {
    struct tm utc;
    time_t tp;
    CHAR dateStr[512], line[1024], summary[2048];
    tp = time(NULL);
    sprintf (summary, "%%%% Date: %s", asctime( gmtime( &tp ) ) );
    strcat (summary, "%% Loaded SFTs: [ " );
    for ( UINT4 X = 0; X < numDetectors; X ++ ) {
      sprintf (line, "%s:%d%s",  multiIFO.sites[X].frDetector.name, mTS->data[X]->length, (X < numDetectors - 1)?", ":" ]\n");
      strcat ( summary, line );
    }
    utc = *XLALGPSToUTC( &utc, (INT4)XLALGPSGetREAL8(&startTime) );
    strcpy ( dateStr, asctime(&utc) );
    dateStr[ strlen(dateStr) - 1 ] = 0;
    snprintf (line, sizeof(line), "%%%% Start GPS time tStart = %12.3f    (%s GMT)\n", XLALGPSGetREAL8(&startTime), dateStr);
    strcat ( summary, line );
    sprintf (line, "%%%% Total amount of data: Tdata = %12.3f s  (%.2f days)\n", Tdata, Tdata/86400 );
    strcat ( summary, line );

    XLAL_CHECK ( (cfg->dataSummary = LALCalloc(1, strlen(summary) + 1 )) != NULL, XLAL_ENOMEM );
    strcpy ( cfg->dataSummary, summary );

    LogPrintfVerbatim( LOG_DEBUG, "%s", cfg->dataSummary );
  } /* write dataSummary string */

  /* free everything not needed any more */
  XLALDestroyMultiTimestamps ( mTS );
  XLALDestroyMultiNoiseWeights ( multiNoiseWeights );
  XLALDestroyMultiDetectorStateSeries ( multiDetStates );

  /* Free ephemeris data */
  XLALDestroyEphemerisData (edat);

  return XLAL_SUCCESS;

} /* InitPFS() */
