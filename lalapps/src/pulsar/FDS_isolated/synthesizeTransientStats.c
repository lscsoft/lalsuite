/*
 * Copyright (C) 2010 Reinhard Prix
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
/** \author R. Prix
 * \file
 * \brief
 * Generate N samples of various statistics (F-stat, B-stat,...) drawn from their
 * respective distributions, assuming Gaussian noise, and drawing signal params from
 * their (given) priors
 *
 * This is based on synthesizeBstat, and is mostly meant to be used for Monte-Carlos
 * studies of ROC curves
 *
 *********************************************************************************/

/*
 *
 * Some possible use-cases to consider
 * - transient search BF-stat (synthesize atoms)
 * - line-veto studies (generate line-realizations)
 * - different B-stats from different prior models (to avoid integration)
 *
 */

#include "config.h"

/* System includes */
#include <stdio.h>
#include <stdbool.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

/* GSL includes */
#include <gsl/gsl_rng.h>

/* LAL-includes */
#include <lal/SkyCoordinates.h>
#include <lal/AVFactories.h>
#include <lal/LALInitBarycenter.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lal/ComputeFstat.h>

#include <lal/TransientCW_utils.h>
#include <lal/ProbabilityDensity.h>

#include <lal/SynthesizeCWDraws.h>

#include <lalapps.h>

/*---------- DEFINES ----------*/
#define EPHEM_YEARS  "00-19-DE405"	/**< default range, covering S5: override with --ephemYear */
#define SQ(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUAD(x) ((x)*(x)*(x)*(x))

/* ---------- local types ---------- */

/** User-variables: can be set from config-file or command-line */
typedef struct {
  BOOLEAN help;		/**< trigger output of help string */

  /* amplitude parameters + ranges: 4 alternative ways to specify the h0-prior (set to < 0 to deactivate all but one!) */
  REAL8 fixedh0Nat;	/**< Alternative 1: if >=0 ==> fix the GW amplitude: h0/sqrt(Sn) */
  REAL8 fixedSNR;	/**< Alternative 2: if >=0 ==> fix the optimal SNR of the injected signals */
  REAL8 fixedh0NatMax;	/**< Alternative 3: if >=0 ==> draw GW amplitude h0 in [0, h0NatMax ]: <==> 'regularized' F-stat prior (obsolete) */
  REAL8 fixedRhohMax;	/**< Alternative 4: if >=0 ==> draw rhoh=h0*(detM)^(1/8) in [0, rhohMax]: <==> canonical F-stat prior */

  REAL8 cosi;		/**< cos(inclination angle). If not set: randomize within [-1,1] */
  REAL8 psi;		/**< polarization angle psi. If not set: randomize within [-pi/4,pi/4] */
  REAL8 phi0;		/**< initial GW phase phi_0. If not set: randomize within [0, 2pi] */
  INT4 AmpPriorType;	/**< enumeration of types of amplitude-priors: 0=physical, 1=canonical */

  /* Doppler parameters */
  REAL8 Alpha;		/**< skyposition Alpha (RA) in radians */
  REAL8 Delta;		/**< skyposition Delta (Dec) in radians */

  /* transient window ranges: for injection ... */
  CHAR *injectWindow_type; 	/**< Type of transient window to inject ('none', 'rect', 'exp'). */
  REAL8 injectWindow_t0Days;	/**< Earliest start-time of transient window to inject, as offset in days from dataStartGPS */
  REAL8 injectWindow_t0DaysBand;/**< Range of start-times of transient window to inject, in days */
  REAL8 injectWindow_tauDays;	/**< Shortest transient-window timescale to inject, in days */
  REAL8 injectWindow_tauDaysBand; /**< Range of transient-window timescale to inject, in days */
  /* ... and for search */
  CHAR *searchWindow_type; 	/**< Type of transient window to search with ('none', 'rect', 'exp'). */
  REAL8 searchWindow_t0Days;	/**< Earliest start-time of transient window to search, as offset in days from dataStartGPS */
  REAL8 searchWindow_t0DaysBand;/**< Range of start-times of transient window to inject, in days */
  REAL8 searchWindow_tauDays;	/**< Shortest transient-window timescale to search, in days */
  REAL8 searchWindow_tauDaysBand; /**< Range of transient-window timescale to search, in days */
  INT4  searchWindow_dt0;	/**< Step-size for search/marginalization over transient-window start-time, in seconds [Default:Tsft] */
  INT4  searchWindow_dtau;       /**< Step-size for search/marginalization over transient-window timescale, in seconds [Default:Tsft] */

  /* other parameters */
  CHAR *IFO;		/**< IFO name */
  INT4 dataStartGPS;	/**< data start-time in GPS seconds */
  INT4 dataDuration;	/**< data-span to generate */
  INT4 TAtom;		/**< Fstat atoms time baseline */

  BOOLEAN computeFtotal; /**< Also compute 'total' F-statistic over the full data-span */
  INT4 numDraws;	/**< number of random 'draws' to simulate for F-stat and B-stat */

  CHAR *outputStats;	/**< output file to write numDraw resulting statistics into */
  CHAR *outputAtoms;	/**< output F-statistic atoms into a file with this basename */
  CHAR *outputFstatMap;	/**< output F-statistic over 2D parameter space {t0, tau} into file with this basename */
  CHAR *outputInjParams;/**< output injection parameters into this file */
  CHAR *outputPosteriors;/**< output posterior pdfs on t0 and tau */
  BOOLEAN SignalOnly;	/**< dont generate noise-draws: will result in non-random 'signal only' values of F and B */

  BOOLEAN useFReg;	/**< use 'regularized' Fstat (1/D)*e^F for marginalization, or 'standard' e^F */

  CHAR *ephemYear;	/**< date-range string on ephemeris-files to use */

  BOOLEAN version;	/**< output version-info */
  INT4 randSeed;	/**< GSL random-number generator seed value to use */
} UserInput_t;

/** Configuration settings required for and defining a coherent pulsar search.
 * These are 'pre-processed' settings, which have been derived from the user-input.
 */
typedef struct {
  AmplitudePrior_t AmpPrior;	/**< amplitude-parameter priors to draw signals from */

  SkyPosition skypos;		/**< (Alpha,Delta,system). Use Alpha < 0 to signal 'allsky' */
  BOOLEAN SignalOnly;		/**< dont generate noise-draws: will result in non-random 'signal only' values of F and B */

  transientWindowRange_t transientSearchRange;	/**< transient-window range for the search (flat priors) */
  transientWindowRange_t transientInjectRange;	/**< transient-window range for injections (flat priors) */

  MultiDetectorStateSeries *multiDetStates;	/**< multi-detector state series covering observation time */

  gsl_rng *rng;			/**< gsl random-number generator */
  CHAR *logString;		/**< logstring for file-output, containing cmdline-options + code VCS version info */

} ConfigVariables;

/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);

int XLALInitUserVars ( UserInput_t *uvar );
int XLALInitCode ( ConfigVariables *cfg, const UserInput_t *uvar );
EphemerisData * XLALInitEphemeris (const CHAR *ephemYear );
int XLALInitAmplitudePrior ( AmplitudePrior_t *AmpPrior, const UserInput_t *uvar );

/* exportable API */

/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

/*---------- empty initializers ---------- */
ConfigVariables empty_ConfigVariables;
UserInput_t empty_UserInput;

/*----------------------------------------------------------------------*/
/* Main Function starts here */
/*----------------------------------------------------------------------*/
/**
 * MAIN function
 * Generates samples of B-stat and F-stat according to their pdfs for given signal-params.
 */
int main(int argc,char *argv[])
{
  UserInput_t uvar = empty_UserInput;
  ConfigVariables cfg = empty_ConfigVariables;		/**< various derived configuration settings */

  vrbflg = 1;	/* verbose error-messages */
  LogSetLevel(lalDebugLevel);

  /* turn off default GSL error handler */
  gsl_set_error_handler_off ();

  /* ----- register and read all user-variables ----- */
  LogSetLevel(lalDebugLevel);

  if ( XLALInitUserVars( &uvar ) != XLAL_SUCCESS ) {
    LogPrintf ( LOG_CRITICAL, "%s: XLALInitUserVars() failed with errno=%d\n", __func__, xlalErrno );
    return 1;
  }

  /* do ALL cmdline and cfgfile handling */
  if ( XLALUserVarReadAllInput ( argc, argv ) != XLAL_SUCCESS ) {
    LogPrintf ( LOG_CRITICAL, "%s: XLALUserVarReadAllInput() failed with errno=%d\n", __func__, xlalErrno );
    return 1;
  }

  if (uvar.help)	/* if help was requested, we're done here */
    return 0;

  if ( uvar.version ) {
    /* output verbose VCS version string if requested */
    CHAR *vcs;
    if ( (vcs = XLALGetVersionString (lalDebugLevel)) == NULL ) {
      LogPrintf ( LOG_CRITICAL, "%s:XLALGetVersionString(%d) failed with errno=%d.\n", __func__, lalDebugLevel, xlalErrno );
      return 1;
    }
    printf ( "%s\n", vcs );
    XLALFree ( vcs );
    return 0;
  }

  /* ---------- Initialize code-setup ---------- */
  if ( XLALInitCode( &cfg, &uvar ) != XLAL_SUCCESS ) {
    LogPrintf (LOG_CRITICAL, "%s: XLALInitCode() failed with error = %d\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  /* ----- prepare stats output ----- */
  FILE *fpTransientStats = NULL;
  if ( uvar.outputStats )
    {
      if ( (fpTransientStats = fopen (uvar.outputStats, "wb")) == NULL)
	{
	  LogPrintf (LOG_CRITICAL, "Error opening file '%s' for writing..\n\n", uvar.outputStats );
	  XLAL_ERROR ( XLAL_EIO );
	}
      fprintf (fpTransientStats, "%s", cfg.logString );		/* write search log comment */
      if ( write_transientCandidate_to_fp ( fpTransientStats, NULL ) != XLAL_SUCCESS ) { /* write header-line comment */
        XLAL_ERROR ( XLAL_EFUNC );
      }
    } /* if outputStats */

  /* ----- prepare injection params output ----- */
  FILE *fpInjParams = NULL;
  if ( uvar.outputInjParams )
    {
      if ( (fpInjParams = fopen (uvar.outputInjParams, "wb")) == NULL)
	{
	  LogPrintf (LOG_CRITICAL, "Error opening file '%s' for writing..\n\n", uvar.outputInjParams );
	  XLAL_ERROR ( XLAL_EIO );
	}
      fprintf (fpInjParams, "%s", cfg.logString );		/* write search log comment */
      if ( write_InjParams_to_fp ( fpInjParams, NULL, 0 ) != XLAL_SUCCESS ) { /* write header-line comment */
        XLAL_ERROR ( XLAL_EFUNC );
      }
    } /* if outputInjParams */

  /* ----- main MC loop over numDraws trials ---------- */
  multiAMBuffer_t multiAMBuffer = empty_multiAMBuffer;	  /* prepare AM-buffer */
  INT4 i;

  for ( i=0; i < uvar.numDraws; i ++ )
    {
      InjParams_t injParamsDrawn = empty_InjParams_t;

      /* ----- generate signal random draws from ranges and generate Fstat atoms */
      MultiFstatAtomVector *multiAtoms;
      multiAtoms = XLALSynthesizeTransientAtoms ( &injParamsDrawn, cfg.skypos, cfg.AmpPrior, cfg.transientInjectRange, cfg.multiDetStates, cfg.SignalOnly, &multiAMBuffer, cfg.rng, -1);
      if ( multiAtoms ==NULL ) {
        LogPrintf ( LOG_CRITICAL, "%s: XLALSynthesizeTransientAtoms() failed with xlalErrno = %d\n", __func__, xlalErrno );
        XLAL_ERROR ( XLAL_EFUNC );
      }

      /* ----- if requested, output signal injection parameters into file */
      if ( fpInjParams && (write_InjParams_to_fp ( fpInjParams, &injParamsDrawn, uvar.dataStartGPS ) != XLAL_SUCCESS ) ) {
        XLAL_ERROR ( XLAL_EFUNC );
      } /* if fpInjParams & failure*/


      /* ----- add meta-info on current transient-CW candidate */
      transientCandidate_t cand = empty_transientCandidate;
      cand.doppler.Alpha = multiAMBuffer.skypos.longitude;
      cand.doppler.Delta = multiAMBuffer.skypos.latitude;
      cand.windowRange   = cfg.transientSearchRange;

      /* ----- if needed: compute transient-Bstat search statistic on these atoms */
      if ( fpTransientStats || uvar.outputFstatMap || uvar.outputPosteriors )
        {
          /* compute Fstat map F_mn over {t0, tau} */
          if ( (cand.FstatMap = XLALComputeTransientFstatMap ( multiAtoms, cand.windowRange, uvar.useFReg)) == NULL ) {
            XLALPrintError ("%s: XLALComputeTransientFstatMap() failed with xlalErrno = %d.\n", __func__, xlalErrno );
            XLAL_ERROR ( XLAL_EFUNC );
          }
        } /* if we'll need the Fstat-map F_mn */

      /* ----- if requested compute marginalized Bayes factor */
      if ( fpTransientStats )
        {
          cand.logBstat = XLALComputeTransientBstat ( cand.windowRange, cand.FstatMap );
          UINT4 err = xlalErrno;
          if ( err ) {
            XLALPrintError ("%s: XLALComputeTransientBstat() failed with xlalErrno = %d\n", __func__, err );
            XLAL_ERROR ( XLAL_EFUNC );
          }

          if ( uvar.SignalOnly )
            {
              cand.FstatMap->maxF += 2;
              cand.logBstat += 2;
            }

        } /* if Bstat requested */

      /* ----- if requested, compute parameter posteriors for {t0, tau} */
      pdf1D_t *pdf_t0  = NULL;
      pdf1D_t *pdf_tau = NULL;
      if ( fpTransientStats || uvar.outputPosteriors )
        {
          if ( (pdf_t0 = XLALComputeTransientPosterior_t0 ( cand.windowRange, cand.FstatMap )) == NULL ) {
            XLALPrintError ("%s: failed to compute t0-posterior\n", __func__ );
            XLAL_ERROR ( XLAL_EFUNC );
          }
          if ( (pdf_tau = XLALComputeTransientPosterior_tau ( cand.windowRange, cand.FstatMap )) == NULL ) {
            XLALPrintError ("%s: failed to compute tau-posterior\n", __func__ );
            XLAL_ERROR ( XLAL_EFUNC );
          }
          /* get maximum-posterior estimate (MP) from the modes of these pdfs */
          cand.t0_MP = XLALFindModeOfPDF1D ( pdf_t0 );
          if ( xlalErrno ) {
            XLALPrintError ("%s: mode-estimation failed for pdf_t0. xlalErrno = %d\n", __func__, xlalErrno );
            XLAL_ERROR ( XLAL_EFUNC );
          }
          cand.tau_MP =  XLALFindModeOfPDF1D ( pdf_tau );
          if ( xlalErrno ) {
            XLALPrintError ("%s: mode-estimation failed for pdf_tau. xlalErrno = %d\n", __func__, xlalErrno );
            XLAL_ERROR ( XLAL_EFUNC );
          }

        } // if posteriors required

      /* ----- if requested, compute Ftotal over full data-span */
      if ( uvar.computeFtotal )
        {
          transientFstatMap_t *FtotalMap;
          /* prepare special window to cover all the data with one F-stat calculation == Ftotal */
          transientWindowRange_t winRangeAll = empty_transientWindowRange;
          winRangeAll.type = TRANSIENT_NONE;

          BOOLEAN useFReg = false;
          if ( (FtotalMap = XLALComputeTransientFstatMap ( multiAtoms, winRangeAll, useFReg)) == NULL ) {
            XLALPrintError ("%s: XLALComputeTransientFstatMap() failed with xlalErrno = %d.\n", __func__, xlalErrno );
            XLAL_ERROR ( XLAL_EFUNC );
          }

          /* we only use twoFtotal = 2 * maxF from this single-Fstat calculation */
          REAL8 twoFtotal = 2.0 * FtotalMap->maxF;
          if ( uvar.SignalOnly )
            twoFtotal += 4;

          /* ugly hack: lacking a good container for twoFtotal, we borrow fkdot[3] for this here ;) [only used for paper-MCs] */
          cand.doppler.fkdot[3] = twoFtotal;

          /* good riddance .. */
          XLALDestroyTransientFstatMap ( FtotalMap );

        } /* if computeFtotal */

      /* ----- if requested, output atoms-vector into file */
      if ( uvar.outputAtoms )
        {

          FILE *fpAtoms;
          char *fnameAtoms;
          UINT4 len = strlen ( uvar.outputAtoms ) + 20;
          if ( (fnameAtoms = XLALCalloc ( 1, len )) == NULL ) {
            XLALPrintError ("%s: failed to XLALCalloc ( 1, %d )\n", __func__, len );
            XLAL_ERROR ( XLAL_EFUNC );
          }
          sprintf ( fnameAtoms, "%s_%04d_of_%04d.dat", uvar.outputAtoms, i + 1, uvar.numDraws );

          if ( ( fpAtoms = fopen ( fnameAtoms, "wb" )) == NULL ) {
            XLALPrintError ("%s: failed to open atoms-output file '%s' for writing.\n", __func__, fnameAtoms );
            XLAL_ERROR ( XLAL_EFUNC );
          }
	  fprintf ( fpAtoms, "%s", cfg.logString );	/* output header info */

	  if ( write_MultiFstatAtoms_to_fp ( fpAtoms, multiAtoms ) != XLAL_SUCCESS ) {
            XLALPrintError ("%s: failed to write atoms to output file '%s'. xlalErrno = %d\n", __func__, fnameAtoms, xlalErrno );
            XLAL_ERROR ( XLAL_EFUNC );
          }

          XLALFree ( fnameAtoms );
	  fclose (fpAtoms);
        } /* if outputAtoms */

      /* ----- if requested, output Fstat-map over {t0, tau} */
      if ( uvar.outputFstatMap )
        {
          FILE *fpFstatMap;
          char *fnameFstatMap;
          UINT4 len = strlen ( uvar.outputFstatMap ) + 20;
          if ( (fnameFstatMap = XLALCalloc ( 1, len )) == NULL ) {
            XLALPrintError ("%s: failed to XLALCalloc ( 1, %d )\n", __func__, len );
            XLAL_ERROR ( XLAL_EFUNC );
          }
          sprintf ( fnameFstatMap, "%s_%04d_of_%04d.dat", uvar.outputFstatMap, i + 1, uvar.numDraws );

          if ( ( fpFstatMap = fopen ( fnameFstatMap, "wb" )) == NULL ) {
            XLALPrintError ("%s: failed to open Fstat-map output file '%s' for writing.\n", __func__, fnameFstatMap );
            XLAL_ERROR ( XLAL_EFUNC );
          }
	  fprintf ( fpFstatMap, "%s", cfg.logString );	/* output header info */

          fprintf (fpFstatMap, "\nFstat_mn = \\\n" );
          if ( XLALfprintfGSLmatrix ( fpFstatMap, "%.9g", cand.FstatMap->F_mn ) != XLAL_SUCCESS ) {
            XLALPrintError ("%s: XLALfprintfGSLmatrix() failed.\n", __func__ );
            XLAL_ERROR ( XLAL_EFUNC );
          }

          XLALFree ( fnameFstatMap );
	  fclose (fpFstatMap);

        } /* if outputFstatMap */

      /* ----- if requested, output posterior pdfs on transient params {t0, tau} into a file */
      if ( uvar.outputPosteriors )
        {
          FILE *fpPosteriors;
          char *fnamePosteriors;
          UINT4 len = strlen ( uvar.outputPosteriors ) + 20;
          if ( (fnamePosteriors = XLALCalloc ( 1, len )) == NULL ) {
            XLALPrintError ("%s: failed to XLALCalloc ( 1, %d )\n", __func__, len );
            XLAL_ERROR ( XLAL_EFUNC );
          }
          sprintf ( fnamePosteriors, "%s_%04d_of_%04d.dat", uvar.outputPosteriors, i + 1, uvar.numDraws );

          if ( ( fpPosteriors = fopen ( fnamePosteriors, "wb" )) == NULL ) {
            XLALPrintError ("%s: failed to open posteriors-output file '%s' for writing.\n", __func__, fnamePosteriors );
            XLAL_ERROR ( XLAL_EFUNC );
          }
	  fprintf ( fpPosteriors, "%s", cfg.logString );	/* output header info */

          /* write them to file, using pdf-method */
	  if ( XLALOutputPDF1D_to_fp ( fpPosteriors, pdf_t0, "pdf_t0" ) != XLAL_SUCCESS ) {
            XLALPrintError ("%s: failed to output t0-posterior to file '%s'.\n", __func__, fnamePosteriors );
            XLAL_ERROR ( XLAL_EFUNC );
          }
	  if ( XLALOutputPDF1D_to_fp ( fpPosteriors, pdf_tau, "pdf_tau" ) != XLAL_SUCCESS ) {
            XLALPrintError ("%s: failed to output tau-posterior to file '%s'.\n", __func__, fnamePosteriors );
            XLAL_ERROR ( XLAL_EFUNC );
          }

          /* free mem, close file */
          XLALFree ( fnamePosteriors );
	  fclose (fpPosteriors);

        } /* if outputPosteriors */


      /* ----- if requested, output transient-cand statistics */
      if ( fpTransientStats && write_transientCandidate_to_fp ( fpTransientStats, &cand ) != XLAL_SUCCESS ) {
        XLALPrintError ( "%s: write_transientCandidate_to_fp() failed.\n", __func__ );
        XLAL_ERROR ( XLAL_EFUNC );
      }

      /* ----- free Memory */
      XLALDestroyTransientFstatMap ( cand.FstatMap );
      XLALDestroyMultiFstatAtomVector ( multiAtoms );
      XLALDestroyPDF1D ( pdf_t0 );
      XLALDestroyPDF1D ( pdf_tau );

    } /* for i < numDraws */

  /* ----- close files ----- */
  if ( fpTransientStats) fclose ( fpTransientStats );
  if ( fpInjParams ) fclose ( fpInjParams );

  /* ----- free memory ---------- */
  XLALDestroyMultiDetectorStateSeries ( cfg.multiDetStates );
  XLALDestroyMultiAMCoeffs ( multiAMBuffer.multiAM );
  XLALDestroyExpLUT();
  /* ----- free amplitude prior pdfs ----- */
  XLALDestroyPDF1D ( cfg.AmpPrior.pdf_h0Nat );
  XLALDestroyPDF1D ( cfg.AmpPrior.pdf_cosi );
  XLALDestroyPDF1D ( cfg.AmpPrior.pdf_psi );
  XLALDestroyPDF1D ( cfg.AmpPrior.pdf_phi0 );

  if ( cfg.logString ) XLALFree ( cfg.logString );
  gsl_rng_free ( cfg.rng );

  XLALDestroyUserVars();

  /* did we forget anything ? (doesn't cover gsl-memory!) */
  LALCheckMemoryLeaks();

  return 0;

} /* main() */

/**
 * Register all our "user-variables" that can be specified from cmd-line and/or config-file.
 * Here we set defaults for some user-variables and register them with the UserInput module.
 */
int
XLALInitUserVars ( UserInput_t *uvar )
{
  /* set a few defaults */
  uvar->help = 0;
  uvar->outputStats = NULL;

  uvar->Alpha = -1;	/* Alpha < 0 indicates "allsky" */
  uvar->Delta = 0;

  uvar->phi0 = 0;
  uvar->psi = 0;

  uvar->dataStartGPS = 814838413;	/* 1 Nov 2005, ~ start of S5 */
  uvar->dataDuration = (INT4) round ( LAL_YRSID_SI );	/* 1 year of data */

  uvar->ephemYear = XLALCalloc (1, strlen(EPHEM_YEARS)+1);
  strcpy (uvar->ephemYear, EPHEM_YEARS);

  uvar->numDraws = 1;
  uvar->TAtom = 1800;

  uvar->computeFtotal = 0;
  uvar->useFReg = 0;

  uvar->fixedh0Nat = -1;
  uvar->fixedSNR = -1;
  uvar->fixedh0NatMax = -1;
  uvar->fixedRhohMax = -1;

#define DEFAULT_IFO "H1"
  uvar->IFO = XLALMalloc ( strlen(DEFAULT_IFO)+1 );
  strcpy ( uvar->IFO, DEFAULT_IFO );

  /* ---------- transient window defaults ---------- */
#define DEFAULT_TRANSIENT "rect"
  uvar->injectWindow_type = XLALMalloc(strlen(DEFAULT_TRANSIENT)+1);
  strcpy ( uvar->injectWindow_type, DEFAULT_TRANSIENT );
  uvar->searchWindow_type = XLALMalloc(strlen(DEFAULT_TRANSIENT)+1);
  strcpy ( uvar->searchWindow_type, DEFAULT_TRANSIENT );

  uvar->injectWindow_tauDays     = 1.0;
  uvar->injectWindow_tauDaysBand = 13.0;

  REAL8 tauMaxDays = ( uvar->injectWindow_tauDays +  uvar->injectWindow_tauDaysBand );
  /* default window-ranges are t0 in [dataStartTime, dataStartTime - 3 * tauMax] */
  uvar->injectWindow_t0Days     = 0; // offset in days from uvar->dataStartGPS
  uvar->injectWindow_t0DaysBand = fmax ( 0.0, 1.0*uvar->dataDuration/DAY24 - TRANSIENT_EXP_EFOLDING * tauMaxDays );	/* make sure it's >= 0 */

  /* search-windows by default identical to inject-windows */
  uvar->searchWindow_t0Days = uvar->injectWindow_t0Days;
  uvar->searchWindow_t0DaysBand = uvar->injectWindow_t0DaysBand;
  uvar->searchWindow_tauDays = uvar->injectWindow_tauDays;
  uvar->searchWindow_tauDaysBand = uvar->injectWindow_tauDaysBand;

  uvar->searchWindow_dt0  = uvar->TAtom;
  uvar->searchWindow_dtau = uvar->TAtom;

  /* register all our user-variables */
  XLALregBOOLUserStruct ( help, 		'h',     UVAR_HELP, "Print this message");

  /* signal Doppler parameters */
  XLALregREALUserStruct ( Alpha, 		'a', UVAR_OPTIONAL, "Sky position alpha (equatorial coordinates) in radians [Default:allsky]");
  XLALregREALUserStruct ( Delta, 		'd', UVAR_OPTIONAL, "Sky position delta (equatorial coordinates) in radians [Default:allsky]");

  /* signal amplitude parameters */
  XLALregREALUserStruct ( fixedh0Nat,		 0, UVAR_OPTIONAL, "Alternative 1: if >=0 fix the GW amplitude: h0/sqrt(Sn)");
  XLALregREALUserStruct ( fixedSNR, 		 0, UVAR_OPTIONAL, "Alternative 2: if >=0 fix the optimal SNR of the injected signals");
  XLALregREALUserStruct ( fixedh0NatMax,	 0, UVAR_OPTIONAL, "Alternative 3: if >=0 draw GW amplitude h0 in [0, h0NatMax ] (FReg prior)");
  XLALregREALUserStruct ( fixedRhohMax, 	 0, UVAR_OPTIONAL, "Alternative 4: if >=0 draw rhoh=h0*(detM)^(1/8) in [0, rhohMax] (canonical F-stat prior)");

  XLALregREALUserStruct ( cosi,			'i', UVAR_OPTIONAL, "cos(inclination angle). If not set: randomize within [-1,1].");
  XLALregREALUserStruct ( psi,			 0,  UVAR_OPTIONAL, "polarization angle psi. If not set: randomize within [-pi/4,pi/4].");
  XLALregREALUserStruct ( phi0,		 	 0,  UVAR_OPTIONAL, "initial GW phase phi_0. If not set: randomize within [0, 2pi]");

  XLALregINTUserStruct  ( AmpPriorType,	 	 0,  UVAR_OPTIONAL, "Enumeration of types of amplitude-priors: 0=physical, 1=canonical");

  XLALregSTRINGUserStruct ( IFO,	        'I', UVAR_OPTIONAL, "Detector: 'G1','L1','H1,'H2', 'V1', ... ");
  XLALregINTUserStruct ( dataStartGPS,	 	 0,  UVAR_OPTIONAL, "data start-time in GPS seconds");
  XLALregINTUserStruct ( dataDuration,	 	 0,  UVAR_OPTIONAL, "data-span to generate (in seconds)");

  /* transient window ranges: for injection ... */
  XLALregSTRINGUserStruct( injectWindow_type,    0, UVAR_OPTIONAL, "Type of transient window to inject ('none', 'rect', 'exp')");
  XLALregREALUserStruct  ( injectWindow_tauDays, 0, UVAR_OPTIONAL, "Shortest transient-window timescale to inject, in days");
  XLALregREALUserStruct  ( injectWindow_tauDaysBand,0,UVAR_OPTIONAL,"Range of transient-window timescale to inject, in days");
  XLALregREALUserStruct  ( injectWindow_t0Days,  0, UVAR_OPTIONAL, "Earliest start-time of transient window to inject, as offset in days from dataStartGPS");
  XLALregREALUserStruct  ( injectWindow_t0DaysBand,0,UVAR_OPTIONAL,"Range of GPS start-time of transient window to inject, in days [Default:dataDuration-3*tauMax]");
  /* ... and for search */
  XLALregSTRINGUserStruct( searchWindow_type,    0, UVAR_OPTIONAL, "Type of transient window to search with ('none', 'rect', 'exp') [Default:injectWindow]");
  XLALregREALUserStruct  ( searchWindow_tauDays, 0, UVAR_OPTIONAL, "Shortest transient-window timescale to search, in days [Default:injectWindow]");
  XLALregREALUserStruct  ( searchWindow_tauDaysBand,0,UVAR_OPTIONAL, "Range of transient-window timescale to search, in days [Default:injectWindow]");
  XLALregREALUserStruct  ( searchWindow_t0Days,  0, UVAR_OPTIONAL, "Earliest start-time of transient window to search, as offset in days from dataStartGPS [Default:injectWindow]");
  XLALregREALUserStruct  ( searchWindow_t0DaysBand,0,UVAR_OPTIONAL, "Range of GPS start-time of transient window to search, in days [Default:injectWindow]");

  XLALregINTUserStruct   ( searchWindow_dtau, 	 0, UVAR_OPTIONAL, "Step-size for search/marginalization over transient-window timescale, in seconds [Default:TAtom]");
  XLALregINTUserStruct   ( searchWindow_dt0, 	 0, UVAR_OPTIONAL, "Step-size for search/marginalization over transient-window start-time, in seconds [Default:TAtom]");

  /* misc params */
  XLALregBOOLUserStruct ( computeFtotal,	 0, UVAR_OPTIONAL, "Also compute 'total' F-statistic over the full data-span" );

  XLALregINTUserStruct  ( numDraws,		'N', UVAR_OPTIONAL,"Number of random 'draws' to simulate");
  XLALregINTUserStruct  ( randSeed,		 0, UVAR_OPTIONAL, "GSL random-number generator seed value to use");

  XLALregSTRINGUserStruct ( outputStats,	'o', UVAR_OPTIONAL, "Output file containing 'numDraws' random draws of stats");
  XLALregSTRINGUserStruct ( outputAtoms,	 0,  UVAR_OPTIONAL, "Output F-statistic atoms into a file with this basename");
  XLALregSTRINGUserStruct ( outputFstatMap,	 0,  UVAR_OPTIONAL, "Output F-statistic over 2D parameter space {t0, tau} into file with this basename");

  XLALregSTRINGUserStruct ( outputInjParams,	 0,  UVAR_OPTIONAL,  "Output injection parameters into this file");
  XLALregSTRINGUserStruct ( outputPosteriors,	 0,  UVAR_OPTIONAL,  "output posterior pdfs on t0 and tau (in octave format) into this file ");

  XLALregBOOLUserStruct ( SignalOnly,        	'S', UVAR_OPTIONAL, "Signal only: generate pure signal without noise");
  XLALregBOOLUserStruct ( useFReg,        	 0,  UVAR_OPTIONAL, "use 'regularized' Fstat (1/D)*e^F (if TRUE) for marginalization, or 'standard' e^F (if FALSE)");

  XLALregSTRINGUserStruct( ephemYear, 	        'y', UVAR_OPTIONAL, "Year (or range of years) of ephemeris files to be used");

  XLALregBOOLUserStruct ( version,        	'V', UVAR_SPECIAL,  "Output code version");

  /* 'hidden' stuff */
  XLALregINTUserStruct ( TAtom,		  	  0, UVAR_DEVELOPER, "Time baseline for Fstat-atoms (typically Tsft) in seconds." );


  if ( xlalErrno ) {
    XLALPrintError ("%s: something failed in initializing user variabels .. errno = %d.\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

} /* XLALInitUserVars() */


/** Initialize Fstat-code: handle user-input and set everything up. */
int
XLALInitCode ( ConfigVariables *cfg, const UserInput_t *uvar )
{
  /* generate log-string for file-output, containing cmdline-options + code VCS version info */
  char *vcs;
  if ( (vcs = XLALGetVersionString(0)) == NULL ) {	  /* short VCS version string */
    XLALPrintError ( "%s: XLALGetVersionString(0) failed with errno=%d.\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }
  char *cmdline;
  if ( (cmdline = XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE )) == NULL ) {
    XLALPrintError ( "%s: XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE ) failed with errno=%d.\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }
  const char fmt[] = "%%%% cmdline: %s\n%%%%\n%s%%%%\n";
  UINT4 len = strlen(vcs) + strlen(cmdline) + strlen(fmt) + 1;
  if ( ( cfg->logString = XLALMalloc ( len  )) == NULL ) {
    XLALPrintError ("%s: XLALMalloc ( %d ) failed.\n", len );
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  sprintf ( cfg->logString, fmt, cmdline, vcs );
  XLALFree ( cmdline );
  XLALFree ( vcs );

  /* trivial settings from user-input */
  cfg->SignalOnly = uvar->SignalOnly;

  /* ----- parse user-input on signal amplitude-paramters + ranges ----- */
  /* skypos */
  cfg->skypos.longitude = uvar->Alpha;	/* Alpha < 0 indicates 'allsky' */
  cfg->skypos.latitude  = uvar->Delta;
  cfg->skypos.system = COORDINATESYSTEM_EQUATORIAL;

  /* ----- amplitude-params: create prior pdfs reflecting the user-input */
  if ( XLALInitAmplitudePrior ( &cfg->AmpPrior, uvar ) != XLAL_SUCCESS )
    XLAL_ERROR ( XLAL_EFUNC );

  /* ----- initialize random-number generator ----- */
  /* read out environment variables GSL_RNG_xxx
   * GSL_RNG_SEED: use to set random seed: default = 0, override by using --randSeed on cmdline
   * GSL_RNG_TYPE: type of random-number generator to use: default = 'mt19937'
   */
  gsl_rng_env_setup ();
  /* allow overriding the random-seed per command-line */
  if ( XLALUserVarWasSet ( &uvar->randSeed ) )
    gsl_rng_default_seed = uvar->randSeed;
  cfg->rng = gsl_rng_alloc (gsl_rng_default);

  LogPrintf ( LOG_DEBUG, "random-number generator type: %s\n", gsl_rng_name (cfg->rng));
  LogPrintf ( LOG_DEBUG, "seed = %lu\n", gsl_rng_default_seed );

  /* init ephemeris-data */
  EphemerisData *edat;
  if ( (edat = XLALInitEphemeris ( uvar->ephemYear )) == NULL ) {
    LogPrintf ( LOG_CRITICAL, "%s: Failed to init ephemeris data for year-span '%s'\n", __func__, uvar->ephemYear );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  /* init detector info */
  LALDetector *site;
  if ( (site = XLALGetSiteInfo ( uvar->IFO )) == NULL ) {
    XLALPrintError ("%s: Failed to get site-info for detector '%s'\n", __func__, uvar->IFO );
    XLAL_ERROR ( XLAL_EFUNC );
  }
  MultiLALDetector *multiDet;
  if ( (multiDet = XLALCreateMultiLALDetector ( 1 )) == NULL ) {   /* currently only implemented for single-IFO case */
    XLALPrintError ("%s: XLALCreateMultiLALDetector(1) failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }
  multiDet->data[0] = (*site); 	/* copy! */
  XLALFree ( site );

  /* init timestamps vector covering observation time */
  UINT4 numSteps = (UINT4) ceil ( uvar->dataDuration / uvar->TAtom );
  MultiLIGOTimeGPSVector * multiTS;
  if ( (multiTS = XLALCalloc ( 1, sizeof(*multiTS))) == NULL ) {
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  multiTS->length = 1;
  if ( (multiTS->data = XLALCalloc (1, sizeof(*multiTS->data))) == NULL ) {
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  if ( (multiTS->data[0] = XLALCreateTimestampVector (numSteps)) == NULL ) {
    XLALPrintError ("%s: XLALCreateTimestampVector(%d) failed.\n", __func__, numSteps );
  }
  multiTS->data[0]->deltaT = uvar->TAtom;
  UINT4 i;
  for ( i=0; i < numSteps; i ++ )
    {
      UINT4 ti = uvar->dataStartGPS + i * uvar->TAtom;
      multiTS->data[0]->data[i].gpsSeconds = ti;
      multiTS->data[0]->data[i].gpsNanoSeconds = 0;
    }

  /* get detector states */
  if ( (cfg->multiDetStates = XLALGetMultiDetectorStates ( multiTS, multiDet, edat, 0.5 * uvar->TAtom )) == NULL ) {
    XLALPrintError ( "%s: XLALGetMultiDetectorStates() failed.\n", __func__ );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  /* get rid of all temporary memory allocated for this step */
  XLALDestroyMultiLALDetector ( multiDet );
  XLALFree(edat->ephemE);
  XLALFree(edat->ephemS);
  XLALFree ( edat );
  XLALDestroyMultiTimestamps ( multiTS );
  multiTS = NULL;


  /* ---------- initialize transient window ranges, for injection ... ---------- */
  transientWindowRange_t InjectRange = empty_transientWindowRange;
  if ( !uvar->injectWindow_type || !strcmp ( uvar->injectWindow_type, "none") )
    InjectRange.type = TRANSIENT_NONE;			/* default: no transient signal window */
  else if ( !strcmp ( uvar->injectWindow_type, "rect" ) )
    InjectRange.type = TRANSIENT_RECTANGULAR;		/* rectangular window [t0, t0+tau] */
  else if ( !strcmp ( uvar->injectWindow_type, "exp" ) )
    InjectRange.type = TRANSIENT_EXPONENTIAL;		/* exponential window starting at t0, charact. time tau */
  else
    {
      XLALPrintError ("%s: Illegal transient inject window '%s' specified: valid are 'none', 'rect' or 'exp'\n", __func__, uvar->injectWindow_type);
      XLAL_ERROR ( XLAL_EINVAL );
    }
  /* make sure user doesn't set window=none but sets window-parameters => indicates she didn't mean 'none' */
  if ( InjectRange.type == TRANSIENT_NONE )
    if ( XLALUserVarWasSet ( &uvar->injectWindow_t0Days ) || XLALUserVarWasSet ( &uvar->injectWindow_t0DaysBand ) ||
         XLALUserVarWasSet ( &uvar->injectWindow_tauDays ) || XLALUserVarWasSet ( &uvar->injectWindow_tauDaysBand ) ) {
      XLALPrintError ("%s: ERROR: injectWindow_type == NONE, but window-parameters were set! Use a different window-type!\n", __func__ );
      XLAL_ERROR ( XLAL_EINVAL );
    }

  if ( uvar->injectWindow_t0DaysBand < 0 || uvar->injectWindow_tauDaysBand < 0 ) {
    XLALPrintError ("%s: only positive t0/tau window injection bands allowed (%d, %f)\n", __func__, uvar->injectWindow_t0DaysBand, uvar->injectWindow_tauDaysBand );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  /* apply correct defaults if unset: t0=dataStart, t0Band=dataDuration-3*tauMax */
  InjectRange.t0 = uvar->dataStartGPS + uvar->injectWindow_t0Days * DAY24;

  REAL8 tauMax = ( uvar->injectWindow_tauDays +  uvar->injectWindow_tauDaysBand ) * DAY24;
  if ( XLALUserVarWasSet (&uvar->injectWindow_t0DaysBand ) )
    InjectRange.t0Band  = uvar->injectWindow_t0DaysBand * DAY24;
  else
    InjectRange.t0Band  = fmax ( 0.0, uvar->dataDuration - TRANSIENT_EXP_EFOLDING * tauMax - InjectRange.t0 ); 	/* make sure it's >= 0 */

  InjectRange.tau     = (UINT4) ( uvar->injectWindow_tauDays * DAY24 );
  InjectRange.tauBand = (UINT4) ( uvar->injectWindow_tauDaysBand * DAY24 );

  cfg->transientInjectRange = InjectRange;

  /* ---------- ... and for search -------------------- */
  transientWindowRange_t SearchRange = empty_transientWindowRange;
  if ( !uvar->searchWindow_type || !strcmp ( uvar->searchWindow_type, "none") )
    SearchRange.type = TRANSIENT_NONE;			/* default: no transient signal window */
  else if ( !strcmp ( uvar->searchWindow_type, "rect" ) )
    SearchRange.type = TRANSIENT_RECTANGULAR;		/* rectangular window [t0, t0+tau] */
  else if ( !strcmp ( uvar->searchWindow_type, "exp" ) )
    SearchRange.type = TRANSIENT_EXPONENTIAL;		/* exponential window starting at t0, charact. time tau */
  else
    {
      XLALPrintError ("%s: Illegal transient search window '%s' specified: valid are 'none', 'rect' or 'exp'\n", __func__, uvar->searchWindow_type);
      XLAL_ERROR ( XLAL_EINVAL );
    }

  /* apply correct defaults if unset: use inect window */
  if ( !XLALUserVarWasSet ( &uvar->searchWindow_type ) )
    SearchRange.type    = InjectRange.type;
  if ( !XLALUserVarWasSet ( &uvar->searchWindow_t0Days ) )
    SearchRange.t0      = InjectRange.t0;
  else
    SearchRange.t0      = uvar->dataStartGPS + uvar->searchWindow_t0Days * DAY24;
  if ( !XLALUserVarWasSet ( &uvar->searchWindow_t0DaysBand ) )
    SearchRange.t0Band = InjectRange.t0Band;
  else
    SearchRange.t0Band  = (UINT4) (uvar->searchWindow_t0DaysBand * DAY24);
  if ( !XLALUserVarWasSet ( &uvar->searchWindow_tauDays ) )
    SearchRange.tau = InjectRange.tau;
  else
    SearchRange.tau     = (UINT4) ( uvar->searchWindow_tauDays * DAY24 );
  if ( !XLALUserVarWasSet ( &uvar->searchWindow_tauDaysBand ) )
    SearchRange.tauBand = InjectRange.tauBand;
  else
    SearchRange.tauBand = (UINT4) ( uvar->searchWindow_tauDaysBand * DAY24 );

  if ( XLALUserVarWasSet ( &uvar->searchWindow_dt0 ) )
    SearchRange.dt0 = uvar->searchWindow_dt0;
  else
    SearchRange.dt0 = uvar->TAtom;

  if ( XLALUserVarWasSet ( &uvar->searchWindow_dtau ) )
    SearchRange.dtau = uvar->searchWindow_dtau;
  else
    SearchRange.dtau = uvar->TAtom;

  /* make sure user doesn't set window=none but sets window-parameters => indicates she didn't mean 'none' */
  if ( SearchRange.type == TRANSIENT_NONE )
    if ( XLALUserVarWasSet ( &uvar->searchWindow_t0Days ) || XLALUserVarWasSet ( &uvar->searchWindow_t0DaysBand ) ||
         XLALUserVarWasSet ( &uvar->searchWindow_tauDays ) || XLALUserVarWasSet ( &uvar->searchWindow_tauDaysBand ) ) {
      XLALPrintError ("%s: ERROR: searchWindow_type == NONE, but window-parameters were set! Use a different window-type!\n", __func__ );
      XLAL_ERROR ( XLAL_EINVAL );
    }

  if (   uvar->searchWindow_t0DaysBand < 0 || uvar->searchWindow_tauDaysBand < 0 ) {
    XLALPrintError ("%s: only positive t0/tau window injection bands allowed (%d, %f)\n", __func__, uvar->searchWindow_t0DaysBand, uvar->searchWindow_tauDaysBand );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  cfg->transientSearchRange = SearchRange;

  return XLAL_SUCCESS;

} /* XLALInitCode() */


/** Load Ephemeris from ephemeris data-files
 */
EphemerisData *
XLALInitEphemeris (const CHAR *ephemYear )	/**< which years do we need? */
{
#define FNAME_LENGTH 1024
  CHAR EphemEarth[FNAME_LENGTH];	/* filename of earth-ephemeris data */
  CHAR EphemSun[FNAME_LENGTH];	/* filename of sun-ephemeris data */

  /* check input consistency */
  if ( !ephemYear ) {
    XLALPrintError ("%s: invalid NULL input for 'ephemYear'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  snprintf(EphemEarth, FNAME_LENGTH, "earth%s.dat", ephemYear);
  snprintf(EphemSun, FNAME_LENGTH, "sun%s.dat",  ephemYear);

  EphemEarth[FNAME_LENGTH-1]=0;
  EphemSun[FNAME_LENGTH-1]=0;

  EphemerisData *edat;
  if ( (edat = XLALInitBarycenter ( EphemEarth, EphemSun)) == NULL ) {
    XLALPrintError ("%s: XLALInitBarycenter() failed.\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  /* return ephemeris */
  return edat;

} /* XLALInitEphemeris() */


/** Initialize amplitude-prior pdfs from the user-input
 */
int
XLALInitAmplitudePrior ( AmplitudePrior_t *AmpPrior, const UserInput_t *uvar )
{
  const UINT4 AmpPriorBins = 100;	// defines the binnning accuracy of our prior-pdfs

  /* consistency check */
  if ( !AmpPrior || !uvar ) {
    XLALPrintError ( "%s: invalid NULL input 'AmpPrior' or 'uvar'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( AmpPrior->pdf_h0Nat || AmpPrior->pdf_cosi || AmpPrior->pdf_psi || AmpPrior->pdf_phi0 ) {
    XLALPrintError ("%s: AmplitudePriors must be set to NULL before calling this function!\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  /* first check that user only provided *one* method of determining the amplitude-prior range */
  UINT4 numSets = 0;
  if ( uvar->fixedh0Nat >= 0 ) numSets ++;
  if ( uvar->fixedSNR >= 0 ) numSets ++;
  if ( uvar->fixedh0NatMax >= 0 ) numSets ++ ;
  if ( uvar->fixedRhohMax >= 0 ) numSets ++;
  if ( numSets != 1 ) {
    XLALPrintError ("%s: Specify (>=0) exactly *ONE* amplitude-prior range of {fixedh0Nat, fixedSNR, fixedh0NatMax, fixedRhohMax}\n", __func__);
    XLAL_ERROR ( XLAL_EINVAL );
  }

  /* ===== first pass: deal with all user-supplied fixed values ==> singular priors! ===== */

  /* ----- h0 ----- */
  if ( uvar->fixedh0Nat >= 0 )	/* fix h0Nat */
    if ( (AmpPrior->pdf_h0Nat = XLALCreateSingularPDF1D ( uvar->fixedh0Nat )) == NULL )
      XLAL_ERROR ( XLAL_EFUNC );

  if ( uvar->fixedSNR >= 0 )   /* dummy-pdf, as signal will be computed with h0Nat=1 then rescaled to fixedSNR */
    if ( (AmpPrior->pdf_h0Nat = XLALCreateSingularPDF1D ( 1.0 )) == NULL )
      XLAL_ERROR ( XLAL_EFUNC );

  AmpPrior->fixedSNR   = uvar->fixedSNR;
  AmpPrior->fixRhohMax = (uvar->fixedRhohMax >= 0);

  /* ----- cosi ----- */
  if ( XLALUserVarWasSet ( &uvar->cosi ) )
    if ( (AmpPrior->pdf_cosi = XLALCreateSingularPDF1D (  uvar->cosi )) == NULL )
      XLAL_ERROR ( XLAL_EFUNC );
  /* ----- psi ----- */
  if ( XLALUserVarWasSet ( &uvar->psi ) )
    if ( (AmpPrior->pdf_psi = XLALCreateSingularPDF1D (  uvar->psi )) == NULL )
      XLAL_ERROR ( XLAL_EFUNC );
  /* ----- phi0 ----- */
  if ( XLALUserVarWasSet ( &uvar->phi0 ) )
    if ( (AmpPrior->pdf_phi0 = XLALCreateSingularPDF1D (  uvar->phi0 )) == NULL )
      XLAL_ERROR ( XLAL_EFUNC );


  /* ===== second pass: deal with non-singular prior ranges, taking into account the type of priors to use */
  REAL8 h0NatMax = 0;
  if (  uvar->fixedh0NatMax >= 0 ) /* draw h0Nat from [0, h0NatMax] */
    h0NatMax = uvar->fixedh0NatMax;
  if ( uvar->fixedRhohMax >= 0 ) /* draw h0 from [0, rhohMax/(detM)^(1/8)] */
    h0NatMax = uvar->fixedRhohMax;	/* at first, will be rescaled by (detM)^(1/8) after the fact */

  switch ( uvar->AmpPriorType )
    {
    case AMP_PRIOR_TYPE_PHYSICAL:

      /* ----- h0 ----- */ // uniform in [0, h0NatMax] : not that 'physical', but simple
      if ( AmpPrior->pdf_h0Nat == NULL )
        if ( (AmpPrior->pdf_h0Nat = XLALCreateUniformPDF1D ( 0, h0NatMax )) == NULL )
          XLAL_ERROR ( XLAL_EFUNC );
      /* ----- cosi ----- */
      if ( AmpPrior->pdf_cosi == NULL )
        if ( (AmpPrior->pdf_cosi = XLALCreateUniformPDF1D ( -1.0, 1.0 )) == NULL )
          XLAL_ERROR ( XLAL_EFUNC );
      /* ----- psi ----- */
      if ( AmpPrior->pdf_psi == NULL )
        if ( (AmpPrior->pdf_psi = XLALCreateUniformPDF1D ( -LAL_PI_4, LAL_PI_4 )) == NULL )
          XLAL_ERROR ( XLAL_EFUNC );
      /* ----- phi0 ----- */
      if ( AmpPrior->pdf_phi0 == NULL )
        if ( (AmpPrior->pdf_phi0 = XLALCreateUniformPDF1D ( 0, LAL_TWOPI )) == NULL )
          XLAL_ERROR ( XLAL_EFUNC );

      break;

    case AMP_PRIOR_TYPE_CANONICAL:
      /* ----- pdf(h0) ~ h0^3 ----- */
      if ( AmpPrior->pdf_h0Nat == NULL )
        {
          UINT4 i;
          pdf1D_t *pdf;
          if ( ( pdf = XLALCreateDiscretePDF1D ( 0, h0NatMax, AmpPriorBins )) == NULL )
            XLAL_ERROR ( XLAL_EFUNC );

          for ( i=0; i < pdf->probDens->length; i ++ )
            {
              REAL8 xMid = 0.5 * ( pdf->xTics->data[i] + pdf->xTics->data[i+1] );
              pdf->probDens->data[i] = CUBE( xMid );	// pdf(h0) ~ h0^3
            }
          AmpPrior->pdf_h0Nat = pdf;
        }
      /* ----- pdf(cosi) ~ ( 1 - cosi^2)^3 ----- */
      if ( AmpPrior->pdf_cosi == NULL )
        {
          UINT4 i;
          pdf1D_t *pdf;
          if ( ( pdf = XLALCreateDiscretePDF1D ( -1.0, 1.0, AmpPriorBins )) == NULL )
            XLAL_ERROR ( XLAL_EFUNC );

          for ( i=0; i < pdf->probDens->length; i ++ )
            {
              REAL8 xMid = 0.5 * ( pdf->xTics->data[i] + pdf->xTics->data[i+1] );
              REAL8 y = 1.0 - SQ(xMid);
              pdf->probDens->data[i] = CUBE( y );
            }
          AmpPrior->pdf_cosi = pdf;
        }
      /* ----- psi ----- */
      if ( AmpPrior->pdf_psi == NULL )
        if ( (AmpPrior->pdf_psi = XLALCreateUniformPDF1D ( -LAL_PI_4, LAL_PI_4 )) == NULL )
          XLAL_ERROR ( XLAL_EFUNC );
      /* ----- phi0 ----- */
      if ( AmpPrior->pdf_phi0 == NULL )
        if ( (AmpPrior->pdf_phi0 = XLALCreateUniformPDF1D ( 0, LAL_TWOPI )) == NULL )
          XLAL_ERROR ( XLAL_EFUNC );

      break;

    default:
      XLALPrintError ("%s: something went wrong ... unknown priorType = %d\n", __func__, uvar->AmpPriorType );
      XLAL_ERROR ( XLAL_EINVAL );
      break;

    } // switch( uvar->AmpPriorType )

  return XLAL_SUCCESS;

} /* XLALInitAmplitudePrior() */
