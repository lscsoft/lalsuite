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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>


/* LAL-includes */
#include <lal/SkyCoordinates.h>
#include <lal/AVFactories.h>
#include <lal/LALInitBarycenter.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lal/ComputeFstat.h>


#include <lalapps.h>


#include "../transientCW_utils.h"
#include "../ProbabilityDensity.h"

/*---------- DEFINES ----------*/
#define EPHEM_YEARS  "05-09"	/**< default range, covering S5: override with --ephemYear */
#define SQ(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUAD(x) ((x)*(x)*(x)*(x))

/* ---------- local types ---------- */

/** Enumeration of allowed amplitude-prior types
 */
typedef enum {
  AMP_PRIOR_TYPE_PHYSICAL = 0,	/**< 'physical' priors: isotropic pdf{cosi,psi,phi0} AND flat pdf(h0) */
  AMP_PRIOR_TYPE_CANONICAL,	/**< 'canonical' priors: uniform in A^mu up to h_max */
  AMP_PRIOR_TYPE_LAST
} AmpPriorType_t;

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


/** Signal (amplitude) parameter ranges
 */
typedef struct {
  pdf1D_t *pdf_h0Nat;	/**< pdf for h0/sqrt{Sn} */
  REAL8 fixedSNR;	/**< alternative 1: adjust h0 to fix the optimal SNR of the signal */
  BOOLEAN fixRhohMax;	/**< alternative 2: draw h0 with fixed rhohMax = h0Max * (detM)^(1/8) <==> canonical Fstat prior */

  pdf1D_t *pdf_cosi;	/**< pdf(cosi) */
  pdf1D_t *pdf_psi;	/**< pdf(psi) */
  pdf1D_t *pdf_phi0;	/**< pdf(phi0) */
} AmplitudePrior_t;

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
  MultiLIGOTimeGPSVector *multiTS;		/**< corresponding timestamps vector for convenience */

  gsl_rng *rng;			/**< gsl random-number generator */
  CHAR *logString;		/**< logstring for file-output, containing cmdline-options + code VCS version info */

} ConfigVariables;

/** struct for buffering of AM-coeffs, if signal for same sky-position is injected
 */
typedef struct {
  SkyPosition skypos;		/**< sky-position for which we have AM-coeffs computed already */
  MultiAMCoeffs *multiAM;;	/**< pre-computed AM-coeffs for skypos */
} multiAMBuffer_t;


/** Hold all (generally) randomly drawn injection parameters: skypos, amplitude-params, M_mu_nu, transient-window, SNR
 */
typedef struct
{
  SkyPosition skypos;
  PulsarAmplitudeParams ampParams;
  PulsarAmplitudeVect ampVect;
  AntennaPatternMatrix M_mu_nu;
  transientWindow_t transientWindow;
  REAL8 SNR;
  REAL8 detM1o8;	// (detMp)^(1/8): rescale param between h0, and rhoh = h0 * (detMp)^(1/8)
} InjParams_t;


/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);

int XLALInitUserVars ( UserInput_t *uvar );
int XLALInitCode ( ConfigVariables *cfg, const UserInput_t *uvar );
EphemerisData * XLALInitEphemeris (const CHAR *ephemYear );

/* exportable API */
int XLALDrawCorrelatedNoise ( PulsarAmplitudeVect n_mu, const gsl_matrix *L, gsl_rng * rng );

FstatAtomVector* XLALGenerateFstatAtomVector ( const LIGOTimeGPSVector *TS, const AMCoeffs *amcoeffs );
MultiFstatAtomVector* XLALGenerateMultiFstatAtomVector ( const  MultiLIGOTimeGPSVector *multiTS, const MultiAMCoeffs *multiAM );

int XLALAddNoiseToFstatAtomVector ( FstatAtomVector *atoms, gsl_rng * rng );
int XLALAddNoiseToMultiFstatAtomVector ( MultiFstatAtomVector *multiAtoms, gsl_rng * rng );

REAL8 XLALAddSignalToFstatAtomVector ( FstatAtomVector* atoms, AntennaPatternMatrix *M_mu_nu, const PulsarAmplitudeVect A_Mu, transientWindow_t transientWindow );
REAL8 XLALAddSignalToMultiFstatAtomVector ( MultiFstatAtomVector* multiAtoms, AntennaPatternMatrix *M_mu_nu, const PulsarAmplitudeVect A_Mu, transientWindow_t transientWindow );

MultiFstatAtomVector *XLALSynthesizeTransientAtoms ( InjParams_t *injParams, const ConfigVariables *cfg, multiAMBuffer_t *multiAMBuffer );

int XLALRescaleMultiFstatAtomVector ( MultiFstatAtomVector* multiAtoms,	REAL8 rescale );
int XLALInitAmplitudePrior ( AmplitudePrior_t *AmpPrior, const UserInput_t *uvar );

int write_InjParams_to_fp ( FILE * fp, const InjParams_t *par, UINT4 dataStartGPS );


/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

/*---------- empty initializers ---------- */
ConfigVariables empty_ConfigVariables;
UserInput_t empty_UserInput;
multiAMBuffer_t empty_multiAMBuffer;
InjParams_t empty_InjParams_t;

/*----------------------------------------------------------------------*/
/* Main Function starts here */
/*----------------------------------------------------------------------*/
/**
 * MAIN function
 * Generates samples of B-stat and F-stat according to their pdfs for given signal-params.
 */
int main(int argc,char *argv[])
{
  const char *fn = __func__;

  UserInput_t uvar = empty_UserInput;
  ConfigVariables cfg = empty_ConfigVariables;		/**< various derived configuration settings */

  lalDebugLevel = 0;
  vrbflg = 1;	/* verbose error-messages */
  LogSetLevel(lalDebugLevel);

  /* turn off default GSL error handler */
  gsl_set_error_handler_off ();

  /* ----- register and read all user-variables ----- */
  if ( XLALGetDebugLevel ( argc, argv, 'v') != XLAL_SUCCESS ) {
    LogPrintf ( LOG_CRITICAL, "%s: XLALGetDebugLevel() failed with errno=%d\n", fn, xlalErrno );
    return 1;
  }
  LogSetLevel(lalDebugLevel);

  if ( XLALInitUserVars( &uvar ) != XLAL_SUCCESS ) {
    LogPrintf ( LOG_CRITICAL, "%s: XLALInitUserVars() failed with errno=%d\n", fn, xlalErrno );
    return 1;
  }

  /* do ALL cmdline and cfgfile handling */
  if ( XLALUserVarReadAllInput ( argc, argv ) != XLAL_SUCCESS ) {
    LogPrintf ( LOG_CRITICAL, "%s: XLALUserVarReadAllInput() failed with errno=%d\n", fn, xlalErrno );
    return 1;
  }

  if (uvar.help)	/* if help was requested, we're done here */
    return 0;

  if ( uvar.version ) {
    /* output verbose VCS version string if requested */
    CHAR *vcs;
    if ( (vcs = XLALGetVersionString (lalDebugLevel)) == NULL ) {
      LogPrintf ( LOG_CRITICAL, "%s:XLALGetVersionString(%d) failed with errno=%d.\n", fn, lalDebugLevel, xlalErrno );
      return 1;
    }
    printf ( "%s\n", vcs );
    XLALFree ( vcs );
    return 0;
  }

  /* ---------- Initialize code-setup ---------- */
  if ( XLALInitCode( &cfg, &uvar ) != XLAL_SUCCESS ) {
    LogPrintf (LOG_CRITICAL, "%s: XLALInitCode() failed with error = %d\n", fn, xlalErrno );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }

  /* ----- prepare stats output ----- */
  FILE *fpTransientStats = NULL;
  if ( uvar.outputStats )
    {
      if ( (fpTransientStats = fopen (uvar.outputStats, "wb")) == NULL)
	{
	  LogPrintf (LOG_CRITICAL, "Error opening file '%s' for writing..\n\n", uvar.outputStats );
	  XLAL_ERROR ( fn, XLAL_EIO );
	}
      fprintf (fpTransientStats, "%s", cfg.logString );		/* write search log comment */
      if ( write_transientCandidate_to_fp ( fpTransientStats, NULL ) != XLAL_SUCCESS ) { /* write header-line comment */
        XLAL_ERROR ( fn, XLAL_EFUNC );
      }
    } /* if outputStats */

  /* ----- prepare injection params output ----- */
  FILE *fpInjParams = NULL;
  if ( uvar.outputInjParams )
    {
      if ( (fpInjParams = fopen (uvar.outputInjParams, "wb")) == NULL)
	{
	  LogPrintf (LOG_CRITICAL, "Error opening file '%s' for writing..\n\n", uvar.outputInjParams );
	  XLAL_ERROR ( fn, XLAL_EIO );
	}
      fprintf (fpInjParams, "%s", cfg.logString );		/* write search log comment */
      if ( write_InjParams_to_fp ( fpInjParams, NULL, 0 ) != XLAL_SUCCESS ) { /* write header-line comment */
        XLAL_ERROR ( fn, XLAL_EFUNC );
      }
    } /* if outputInjParams */

  /* ----- main MC loop over numDraws trials ---------- */
  multiAMBuffer_t multiAMBuffer = empty_multiAMBuffer;	  /* prepare AM-buffer */
  INT4 i;

  for ( i=0; i < uvar.numDraws; i ++ )
    {
      InjParams_t injParams = empty_InjParams_t;

      /* ----- generate signal random draws from ranges and generate Fstat atoms */
      MultiFstatAtomVector *multiAtoms;
      if ( (multiAtoms = XLALSynthesizeTransientAtoms ( &injParams, &cfg, &multiAMBuffer )) == NULL ) {
        LogPrintf ( LOG_CRITICAL, "%s: XLALSynthesizeTransientAtoms() failed with xlalErrno = %d\n", fn, xlalErrno );
        XLAL_ERROR ( fn, XLAL_EFUNC );
      }

      /* ----- if requested, output signal injection parameters into file */
      if ( fpInjParams && (write_InjParams_to_fp ( fpInjParams, &injParams, uvar.dataStartGPS ) != XLAL_SUCCESS ) ) {
        XLAL_ERROR ( fn, XLAL_EFUNC );
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
            XLALPrintError ("%s: XLALComputeTransientFstatMap() failed with xlalErrno = %d.\n", fn, xlalErrno );
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }
        } /* if we'll need the Fstat-map F_mn */

      /* ----- if requested compute marginalized Bayes factor */
      if ( fpTransientStats )
        {
          cand.logBstat = XLALComputeTransientBstat ( cand.windowRange, cand.FstatMap );
          UINT4 err = xlalErrno;
          if ( err ) {
            XLALPrintError ("%s: XLALComputeTransientBstat() failed with xlalErrno = %d\n", fn, err );
            XLAL_ERROR ( fn, XLAL_EFUNC );
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
            XLALPrintError ("%s: failed to compute t0-posterior\n", fn );
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }
          if ( (pdf_tau = XLALComputeTransientPosterior_tau ( cand.windowRange, cand.FstatMap )) == NULL ) {
            XLALPrintError ("%s: failed to compute tau-posterior\n", fn );
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }
          /* get maximum-posterior estimate (MP) from the modes of these pdfs */
          cand.t0_MP = XLALFindModeOfPDF1D ( pdf_t0 );
          if ( xlalErrno ) {
            XLALPrintError ("%s: mode-estimation failed for pdf_t0. xlalErrno = %d\n", fn, xlalErrno );
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }
          cand.tau_MP =  XLALFindModeOfPDF1D ( pdf_tau );
          if ( xlalErrno ) {
            XLALPrintError ("%s: mode-estimation failed for pdf_tau. xlalErrno = %d\n", fn, xlalErrno );
            XLAL_ERROR ( fn, XLAL_EFUNC );
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
            XLALPrintError ("%s: XLALComputeTransientFstatMap() failed with xlalErrno = %d.\n", fn, xlalErrno );
            XLAL_ERROR ( fn, XLAL_EFUNC );
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
            XLALPrintError ("%s: failed to XLALCalloc ( 1, %d )\n", fn, len );
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }
          sprintf ( fnameAtoms, "%s_%04d_of_%04d.dat", uvar.outputAtoms, i + 1, uvar.numDraws );

          if ( ( fpAtoms = fopen ( fnameAtoms, "wb" )) == NULL ) {
            XLALPrintError ("%s: failed to open atoms-output file '%s' for writing.\n", fn, fnameAtoms );
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }
	  fprintf ( fpAtoms, "%s", cfg.logString );	/* output header info */

	  if ( write_MultiFstatAtoms_to_fp ( fpAtoms, multiAtoms ) != XLAL_SUCCESS ) {
            XLALPrintError ("%s: failed to write atoms to output file '%s'. xlalErrno = %d\n", fn, fnameAtoms, xlalErrno );
            XLAL_ERROR ( fn, XLAL_EFUNC );
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
            XLALPrintError ("%s: failed to XLALCalloc ( 1, %d )\n", fn, len );
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }
          sprintf ( fnameFstatMap, "%s_%04d_of_%04d.dat", uvar.outputFstatMap, i + 1, uvar.numDraws );

          if ( ( fpFstatMap = fopen ( fnameFstatMap, "wb" )) == NULL ) {
            XLALPrintError ("%s: failed to open Fstat-map output file '%s' for writing.\n", fn, fnameFstatMap );
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }
	  fprintf ( fpFstatMap, "%s", cfg.logString );	/* output header info */

          fprintf (fpFstatMap, "\nFstat_mn = \\\n" );
          if ( XLALfprintfGSLmatrix ( fpFstatMap, "%.9g", cand.FstatMap->F_mn ) != XLAL_SUCCESS ) {
            XLALPrintError ("%s: XLALfprintfGSLmatrix() failed.\n", fn );
            XLAL_ERROR ( fn, XLAL_EFUNC );
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
            XLALPrintError ("%s: failed to XLALCalloc ( 1, %d )\n", fn, len );
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }
          sprintf ( fnamePosteriors, "%s_%04d_of_%04d.dat", uvar.outputPosteriors, i + 1, uvar.numDraws );

          if ( ( fpPosteriors = fopen ( fnamePosteriors, "wb" )) == NULL ) {
            XLALPrintError ("%s: failed to open posteriors-output file '%s' for writing.\n", fn, fnamePosteriors );
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }
	  fprintf ( fpPosteriors, "%s", cfg.logString );	/* output header info */

          /* write them to file, using pdf-method */
	  if ( XLALOutputPDF1D_to_fp ( fpPosteriors, pdf_t0, "pdf_t0" ) != XLAL_SUCCESS ) {
            XLALPrintError ("%s: failed to output t0-posterior to file '%s'.\n", fn, fnamePosteriors );
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }
	  if ( XLALOutputPDF1D_to_fp ( fpPosteriors, pdf_tau, "pdf_tau" ) != XLAL_SUCCESS ) {
            XLALPrintError ("%s: failed to output tau-posterior to file '%s'.\n", fn, fnamePosteriors );
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }

          /* free mem, close file */
          XLALFree ( fnamePosteriors );
	  fclose (fpPosteriors);

        } /* if outputPosteriors */


      /* ----- if requested, output transient-cand statistics */
      if ( fpTransientStats && write_transientCandidate_to_fp ( fpTransientStats, &cand ) != XLAL_SUCCESS ) {
        XLALPrintError ( "%s: write_transientCandidate_to_fp() failed.\n", fn );
        XLAL_ERROR ( fn, XLAL_EFUNC );
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
  XLALDestroyMultiTimestamps ( cfg.multiTS );
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
  const char *fn = __func__;

  /* set a few defaults */
  uvar->help = 0;
  uvar->outputStats = NULL;

  uvar->Alpha = -1;	/* Alpha < 0 indicates "allsky" */
  uvar->Delta = 0;

  uvar->phi0 = 0;
  uvar->psi = 0;

  uvar->dataStartGPS = 814838413;	/* 1 Nov 2005, ~ start of S5 */
  uvar->dataDuration = LAL_YRSID_SI;	/* 1 year of data */

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
    XLALPrintError ("%s: something failed in initializing user variabels .. errno = %d.\n", fn, xlalErrno );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

} /* XLALInitUserVars() */


/** Initialize Fstat-code: handle user-input and set everything up. */
int
XLALInitCode ( ConfigVariables *cfg, const UserInput_t *uvar )
{
  const char *fn = __func__;

  /* generate log-string for file-output, containing cmdline-options + code VCS version info */
  char *vcs;
  if ( (vcs = XLALGetVersionString(0)) == NULL ) {	  /* short VCS version string */
    XLALPrintError ( "%s: XLALGetVersionString(0) failed with errno=%d.\n", fn, xlalErrno );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }
  char *cmdline;
  if ( (cmdline = XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE )) == NULL ) {
    XLALPrintError ( "%s: XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE ) failed with errno=%d.\n", fn, xlalErrno );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }
  const char fmt[] = "%%%% cmdline: %s\n%%%%\n%s%%%%\n";
  UINT4 len = strlen(vcs) + strlen(cmdline) + strlen(fmt) + 1;
  if ( ( cfg->logString = XLALMalloc ( len  )) == NULL ) {
    XLALPrintError ("%s: XLALMalloc ( %d ) failed.\n", len );
    XLAL_ERROR ( fn, XLAL_ENOMEM );
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
    XLAL_ERROR ( fn, XLAL_EFUNC );

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
    LogPrintf ( LOG_CRITICAL, "%s: Failed to init ephemeris data for year-span '%s'\n", fn, uvar->ephemYear );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }

  /* init detector info */
  LALDetector *site;
  if ( (site = XLALGetSiteInfo ( uvar->IFO )) == NULL ) {
    XLALPrintError ("%s: Failed to get site-info for detector '%s'\n", fn, uvar->IFO );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }
  MultiLALDetector *multiDet;
  if ( (multiDet = XLALCreateMultiLALDetector ( 1 )) == NULL ) {   /* currently only implemented for single-IFO case */
    XLALPrintError ("%s: XLALCreateMultiLALDetector(1) failed with errno=%d\n", fn, xlalErrno );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }
  multiDet->data[0] = (*site); 	/* copy! */
  XLALFree ( site );

  /* init timestamps vector covering observation time */
  UINT4 numSteps = (UINT4) ceil ( uvar->dataDuration / uvar->TAtom );
  MultiLIGOTimeGPSVector * multiTS;
  if ( (multiTS = XLALCalloc ( 1, sizeof(*multiTS))) == NULL ) {
    XLAL_ERROR ( fn, XLAL_ENOMEM );
  }
  multiTS->length = 1;
  if ( (multiTS->data = XLALCalloc (1, sizeof(*multiTS->data))) == NULL ) {
    XLAL_ERROR ( fn, XLAL_ENOMEM );
  }
  if ( (multiTS->data[0] = XLALCreateTimestampVector (numSteps)) == NULL ) {
    XLALPrintError ("%s: XLALCreateTimestampVector(%d) failed.\n", fn, numSteps );
  }
  multiTS->data[0]->deltaT = uvar->TAtom;
  UINT4 i;
  for ( i=0; i < numSteps; i ++ )
    {
      UINT4 ti = uvar->dataStartGPS + i * uvar->TAtom;
      multiTS->data[0]->data[i].gpsSeconds = ti;
      multiTS->data[0]->data[i].gpsNanoSeconds = 0;
    }
  cfg->multiTS = multiTS;   /* store timestamps vector */

  /* get detector states */
  if ( (cfg->multiDetStates = XLALGetMultiDetectorStates ( multiTS, multiDet, edat, 0.5 * uvar->TAtom )) == NULL ) {
    XLALPrintError ( "%s: XLALGetMultiDetectorStates() failed.\n", fn );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }
  XLALDestroyMultiLALDetector ( multiDet );
  XLALFree(edat->ephemE);
  XLALFree(edat->ephemS);
  XLALFree ( edat );


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
      XLALPrintError ("%s: Illegal transient inject window '%s' specified: valid are 'none', 'rect' or 'exp'\n", fn, uvar->injectWindow_type);
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }
  /* make sure user doesn't set window=none but sets window-parameters => indicates she didn't mean 'none' */
  if ( InjectRange.type == TRANSIENT_NONE )
    if ( XLALUserVarWasSet ( &uvar->injectWindow_t0Days ) || XLALUserVarWasSet ( &uvar->injectWindow_t0DaysBand ) ||
         XLALUserVarWasSet ( &uvar->injectWindow_tauDays ) || XLALUserVarWasSet ( &uvar->injectWindow_tauDaysBand ) ) {
      XLALPrintError ("%s: ERROR: injectWindow_type == NONE, but window-parameters were set! Use a different window-type!\n", fn );
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }

  if ( uvar->injectWindow_t0DaysBand < 0 || uvar->injectWindow_tauDaysBand < 0 ) {
    XLALPrintError ("%s: only positive t0/tau window injection bands allowed (%d, %f)\n", fn, uvar->injectWindow_t0DaysBand, uvar->injectWindow_tauDaysBand );
    XLAL_ERROR ( fn, XLAL_EINVAL );
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
      XLALPrintError ("%s: Illegal transient search window '%s' specified: valid are 'none', 'rect' or 'exp'\n", fn, uvar->searchWindow_type);
      XLAL_ERROR ( fn, XLAL_EINVAL );
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
      XLALPrintError ("%s: ERROR: searchWindow_type == NONE, but window-parameters were set! Use a different window-type!\n", fn );
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }

  if (   uvar->searchWindow_t0DaysBand < 0 || uvar->searchWindow_tauDaysBand < 0 ) {
    XLALPrintError ("%s: only positive t0/tau window injection bands allowed (%d, %f)\n", fn, uvar->searchWindow_t0DaysBand, uvar->searchWindow_tauDaysBand );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  cfg->transientSearchRange = SearchRange;

  return XLAL_SUCCESS;

} /* XLALInitCode() */


/** Generate 4 random-noise draws n_mu = {n_1, n_2, n_3, n_4} with correlations according to
 * the matrix M = L L^T, which is passed in as input.
 *
 * Note: you need to pass a pre-allocated 4-vector n_mu.
 * Note2: this function is meant as a lower-level noise-generation utility, called
 * from a higher-level function to translate the antenna-pattern functions into pre-factorized Lcor
 */
int
XLALDrawCorrelatedNoise ( PulsarAmplitudeVect n_mu,	/**< [out] 4d vector of noise-components {n_mu}, with correlation L * L^T */
                          const gsl_matrix *L,		/**< [in] correlator matrix to get n_mu = L_mu_nu * norm_nu from 4 uncorr. unit variates norm_nu */
                          gsl_rng * rng			/**< gsl random-number generator */
                          )
{
  const CHAR *fn = __func__;
  int gslstat;

  /* ----- check input arguments ----- */
  if ( !L || (L->size1 != 4) || (L->size2 != 4) ) {
    XLALPrintError ( "%s: Invalid correlator matrix, n_mu must be pre-allocate 4x4 matrix", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( !rng ) {
    XLALPrintError ("%s: invalid NULL input as gsl random-number generator!\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  /* ----- generate 4 normal-distributed, uncorrelated random numbers ----- */
  PulsarAmplitudeVect n;

  n[0] = gsl_ran_gaussian ( rng, 1.0 );
  n[1] = gsl_ran_gaussian ( rng, 1.0 );
  n[2] = gsl_ran_gaussian ( rng, 1.0 );
  n[3] = gsl_ran_gaussian ( rng, 1.0 );

  /* use four normal-variates {norm_nu} with correlator matrix L to get: n_mu = L_{mu nu} norm_nu
   * which gives {n_\mu} satisfying cov(n_mu,n_nu) = (L L^T)_{mu nu} = M_{mu nu}
   */
  gsl_vector_view n_view = gsl_vector_view_array ( n, 4 );
  gsl_vector_view n_mu_view = gsl_vector_view_array ( n_mu, 4 );

  /* int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)
   * compute the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H
   * for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
   */
  if ( (gslstat = gsl_blas_dgemv (CblasNoTrans, 1.0, L, &n_view.vector, 0.0, &n_mu_view.vector)) != 0 ) {
    XLALPrintError ( "%s: gsl_blas_dgemv(L * norm) failed: %s\n", fn, gsl_strerror (gslstat) );
    XLAL_ERROR ( fn, XLAL_EFAILED );
  }

  return XLAL_SUCCESS;

} /* XLALDrawCorrelatedNoise() */

/** Generate an FstatAtomVector for given antenna-pattern functions.
 * Simply creates FstatAtomVector and initializes with antenna-pattern function.
 */
FstatAtomVector*
XLALGenerateFstatAtomVector ( const LIGOTimeGPSVector *TS,	/**< input timestamps vector t_i */
                              const AMCoeffs *amcoeffs		/**< input antenna-pattern functions {a_i, b_i} */
                              )
{
  const char *fn = __func__;

  /* check input consistency */
  if ( !TS || !TS->data ) {
    XLALPrintError ("%s: invalid NULL input in TS=%p or TS->data=%p\n", fn, TS, TS->data );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }
  if ( !amcoeffs || !amcoeffs->a || !amcoeffs->b ) {
    XLALPrintError ("%s: invalid NULL input in amcoeffs=%p or amcoeffs->a=%p, amcoeffs->b=%p\n", fn, amcoeffs, amcoeffs->a, amcoeffs->b );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }
  UINT4 numAtoms = TS->length;
  if ( numAtoms != amcoeffs->a->length || numAtoms != amcoeffs->b->length ) {
    XLALPrintError ("%s: inconsistent lengths numTS=%d amcoeffs->a = %d, amecoeffs->b = %d\n", fn, numAtoms, amcoeffs->a->length, amcoeffs->b->length );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  /* prepare output vector */
  FstatAtomVector *atoms;
  if ( ( atoms = XLALCreateFstatAtomVector ( numAtoms ) ) == NULL ) {
    XLALPrintError ("%s: XLALCreateFstatAtomVector(%d) failed.\n", fn, numAtoms );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }
  atoms->TAtom = TS->deltaT;

  UINT4 alpha;
  for ( alpha=0; alpha < numAtoms; alpha ++ )
    {
      REAL8 a = amcoeffs->a->data[alpha];
      REAL8 b = amcoeffs->b->data[alpha];

      atoms->data[alpha].timestamp = TS->data[alpha].gpsSeconds;	/* don't care about nanoseconds for atoms */
      atoms->data[alpha].a2_alpha = a * a;
      atoms->data[alpha].b2_alpha = b * b;
      atoms->data[alpha].ab_alpha = a * b;

      /* Fa,Fb are zero-initialized from XLALCreateFstatAtomVector() */

    } /* for alpha < numAtoms */

  /* return result */
  return atoms;

} /* XLALGenerateFstatAtomVector() */


/** Generate a multi-FstatAtomVector for given antenna-pattern functions.
 * Simply creates MultiFstatAtomVector and initializes with antenna-pattern function.
 */
MultiFstatAtomVector*
XLALGenerateMultiFstatAtomVector ( const MultiLIGOTimeGPSVector *multiTS,	/**< input multi-timestamps vector t_i */
                                   const MultiAMCoeffs *multiAM			/**< input antenna-pattern functions {a_i, b_i} */
                                   )
{
  const char *fn = __func__;

  /* check input consistency */
  if ( !multiTS || !multiTS->data ) {
    XLALPrintError ("%s: invalid NULL input in 'multiTS'\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }
  if ( !multiAM || !multiAM->data || !multiAM->data[0] ) {
    XLALPrintError ("%s: invalid NULL input in 'mutiAM'\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  UINT4 numDet = multiTS->length;
  if ( numDet != multiAM->length ) {
    XLALPrintError ("%s: inconsistent number of detectors in multiTS (%d) and multiAM (%d)\n", fn, multiTS->length, multiAM->length );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  /* create multi-atoms vector */
  MultiFstatAtomVector *multiAtoms;
  if ( ( multiAtoms = XLALCalloc ( 1, sizeof(*multiAtoms) )) == NULL ) {
    XLALPrintError ("%s: XLALCalloc ( 1, %d) failed.\n", fn, sizeof(*multiAtoms) );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }
  multiAtoms->length = numDet;
  if ( ( multiAtoms->data = XLALCalloc ( numDet, sizeof(*multiAtoms->data) ) ) == NULL ) {
    XLALPrintError ("%s: XLALCalloc ( %d, %d) failed.\n", fn, numDet, sizeof(*multiAtoms->data) );
    XLALFree ( multiAtoms );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  /* loop over detectors and generate each atoms-vector individually */
  UINT4 X;
  for ( X=0; X < numDet; X ++ )
    {
      if ( ( multiAtoms->data[X] = XLALGenerateFstatAtomVector ( multiTS->data[X], multiAM->data[X] )) == NULL ) {
        XLALPrintError ("%s: XLALGenerateFstatAtomVector() failed.\n", fn );
        XLALDestroyMultiFstatAtomVector ( multiAtoms );
        XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
      }

    } /* for X < numDet */

  /* return result */
  return multiAtoms;

} /* XLALGenerateMultiFstatAtomVector() */

/** Add Gaussian-noise components to given FstatAtomVector
 */
int
XLALAddNoiseToFstatAtomVector ( FstatAtomVector *atoms,	/**< input atoms-vector, noise will be added to this */
                                gsl_rng * rng		/**< random-number generator */
                                )
{
  const char *fn = __func__;

  /* check input consistency */
  if ( !atoms || !rng ) {
    XLALPrintError ("%s: invalid NULL input for 'atoms'=%p or random-number generator 'rng'=%p\n", fn, atoms, rng );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  UINT4 numAtoms = atoms->length;

  /* prepare gsl-matrix for correlator L = 1/2 * [ a, a ; b , b ] */
  gsl_matrix *Lcor;
  if ( (Lcor = gsl_matrix_calloc ( 4, 4 )) == NULL ) {
    XLALPrintError ("%s: gsl_matrix_calloc ( 4, 4 ) failed.\n", fn );
    XLAL_ERROR ( fn, XLAL_ENOMEM );
  }
  /* prepare placeholder for 4 n_mu noise draws */
  PulsarAmplitudeVect n_mu;

  /* ----- step through atoms and synthesize noise ----- */
  UINT4 alpha;
  for ( alpha=0; alpha < numAtoms; alpha ++ )
    {
      /* unfortunately we need {a,b} here, but
       * the atoms only store {a^2, b^2, ab }
       * so we need to invert this [module arbitrary relative sign)
       */
      REAL8 a2 = atoms->data[alpha].a2_alpha;
      REAL8 b2 = atoms->data[alpha].b2_alpha;
      REAL8 ab = atoms->data[alpha].ab_alpha;

      REAL8 a_by_2 = 0.5 * sqrt(a2);
      REAL8 b_by_2 = 0.5 * sqrt(b2);
      /* convention: always set sign on b */
      if ( ab < 0 )
        b_by_2 = - b_by_2;

      /* upper-left block */
      gsl_matrix_set ( Lcor, 0, 0, a_by_2 );
      gsl_matrix_set ( Lcor, 0, 1, a_by_2 );
      gsl_matrix_set ( Lcor, 1, 0, b_by_2);
      gsl_matrix_set ( Lcor, 1, 1, b_by_2);
      /* lower-right block: +2 on all components */
      gsl_matrix_set ( Lcor, 2, 2, a_by_2 );
      gsl_matrix_set ( Lcor, 2, 3, a_by_2 );
      gsl_matrix_set ( Lcor, 3, 2, b_by_2);
      gsl_matrix_set ( Lcor, 3, 3, b_by_2);

      if ( XLALDrawCorrelatedNoise ( n_mu, Lcor, rng ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: failed to XLALDrawCorrelatedNoise().\n", fn );
        XLAL_ERROR ( fn, XLAL_EFUNC );
      }
      REAL8 x1,x2,x3,x4;
      x1 = n_mu[0];
      x2 = n_mu[1];
      x3 = n_mu[2];
      x4 = n_mu[3];

      /* add this to Fstat-atom */
      /* relation Fa,Fb <--> x_mu: see Eq.(72) in CFSv2-LIGO-T0900149-v2.pdf */
      atoms->data[alpha].Fa_alpha.re +=   x1;
      atoms->data[alpha].Fa_alpha.im += - x3;
      atoms->data[alpha].Fb_alpha.re +=   x2;
      atoms->data[alpha].Fb_alpha.im += - x4;

    } /* for i < numAtoms */

  /* free internal memory */
  gsl_matrix_free ( Lcor );

  return XLAL_SUCCESS;

} /* XLALAddNoiseToFstatAtomVector() */


/** Add Gaussian-noise components to given multi-FstatAtomVector
 */
int
XLALAddNoiseToMultiFstatAtomVector ( MultiFstatAtomVector *multiAtoms,	/**< input multi atoms-vector, noise will be added to this */
                                     gsl_rng * rng			/**< random-number generator */
                                     )
{
  const char *fn = __func__;

  /* check input consistency */
  if ( !multiAtoms || !multiAtoms->data ) {
    XLALPrintError ("%s: invalid NULL input in 'multiAtoms'\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( !rng ) {
    XLALPrintError ("%s: invalid NULL input for random-number generator 'rng'\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  UINT4 numDetectors = multiAtoms->length;

  UINT4 X;
  for ( X=0; X < numDetectors; X ++ )
    {
      if ( XLALAddNoiseToFstatAtomVector ( multiAtoms->data[X], rng ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: XLALAddNoiseToFstatAtomVector() failed.\n", fn );
        XLAL_ERROR ( fn, XLAL_EFUNC );
      }

    } /* for X < numDetectors */

  return XLAL_SUCCESS;

} /* XLALSynthesizeMultiFstatAtomVector4Noise() */


/** Add signal s_mu = M_mu_nu A^nu within the given transient-window
 *  to given atoms.
 *
 * RETURN: SNR^2 of the injected signal
 * and the effective AntennaPatternMatrix M_mu_nu for this signal.
 */
REAL8
XLALAddSignalToFstatAtomVector ( FstatAtomVector* atoms,	 /**< [in/out] atoms vectors containing antenna-functions and possibly noise {Fa,Fb} */
                                 AntennaPatternMatrix *M_mu_nu,	 /**< [out] effective antenna-pattern matrix for the injected signal */
                                 const PulsarAmplitudeVect A_Mu, /**< [in] input canonical amplitude vector A^mu = {A1,A2,A3,A4} */
                                 transientWindow_t transientWindow /**< transient signal window */
                                 )
{
  const char *fn = __func__;
  int gslstat;

  /* check input consistency */
  if ( !atoms || !atoms->data ) {
    XLALPrintError ( "%s: Invalid NULL input 'atoms'\n", fn );
    XLAL_ERROR_REAL8 ( fn, XLAL_EINVAL );
  }
  if ( !M_mu_nu ) {
    XLALPrintError ( "%s: Invalid NULL input 'M_mu_nu'\n", fn );
    XLAL_ERROR_REAL8 ( fn, XLAL_EINVAL );
  }

  /* prepare transient-window support */
  UINT4 t0, t1;
  if ( XLALGetTransientWindowTimespan ( &t0, &t1, transientWindow ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALGetTransientWindowTimespan() failed.\n", fn );
    XLAL_ERROR_REAL8 ( fn, XLAL_EFUNC );
  }

  /* prepare gsl-matrix for Mh_mu_nu = [ a^2, a*b ; a*b , b^2 ] */
  gsl_matrix *Mh_mu_nu;
  if ( (Mh_mu_nu = gsl_matrix_calloc ( 4, 4 )) == NULL ) {
    XLALPrintError ("%s: gsl_matrix_calloc ( 4, 4 ) failed.\n", fn );
    XLAL_ERROR_REAL8 ( fn, XLAL_ENOMEM );
  }

  gsl_vector_const_view A_Mu_view = gsl_vector_const_view_array ( A_Mu, 4 );

  REAL8 TAtom = atoms->TAtom;
  UINT4 numAtoms = atoms->length;
  UINT4 alpha;
  REAL8 Ad = 0, Bd = 0, Cd = 0;		// usual non-transient antenna-pattern functions
  REAL8 Ap = 0, Bp = 0, Cp = 0;		// "effective" transient-CW antenna-pattern matrix M'_mu_nu

  for ( alpha=0; alpha < numAtoms; alpha ++ )
    {
      UINT4 ti = atoms->data[alpha].timestamp;
      REAL8 win = XLALGetTransientWindowValue ( ti, t0, t1, transientWindow.tau, transientWindow.type );

      if ( win == 0 )
        continue;

      /* compute sh_mu = sqrt(gamma/2) * Mh_mu_nu A^nu * win, where Mh_mu_nu is now just
       * the per-atom block matrix [a^2,  ab; ab, b^2 ]
       * where Sn=1, so gamma = Sinv*TAtom = TAtom
       */
      // NOTE: for sh_mu: only LINEAR in window-function, NOT quadratic! -> see notes
      REAL8 a2 = win * atoms->data[alpha].a2_alpha;
      REAL8 b2 = win * atoms->data[alpha].b2_alpha;
      REAL8 ab = win * atoms->data[alpha].ab_alpha;

      Ad += a2;
      Bd += b2;
      Cd += ab;

      // we also compute M'_mu_nu, which will be used to estimate optimal SNR
      // NOTE: M'_mu_nu is QUADRATIC in window-function!, so we multiply by win again
      Ap += win * a2;
      Bp += win * b2;
      Cp += win * ab;

      /* upper-left block */
      gsl_matrix_set ( Mh_mu_nu, 0, 0, a2 );
      gsl_matrix_set ( Mh_mu_nu, 1, 1, b2 );
      gsl_matrix_set ( Mh_mu_nu, 0, 1, ab );
      gsl_matrix_set ( Mh_mu_nu, 1, 0, ab );
      /* lower-right block: +2 on all components */
      gsl_matrix_set ( Mh_mu_nu, 2, 2, a2 );
      gsl_matrix_set ( Mh_mu_nu, 3, 3, b2 );
      gsl_matrix_set ( Mh_mu_nu, 2, 3, ab );
      gsl_matrix_set ( Mh_mu_nu, 3, 2, ab );

      /* placeholder for 4-vector sh_mu */
      PulsarAmplitudeVect sh_mu = {0,0,0,0};
      gsl_vector_view sh_mu_view = gsl_vector_view_array ( sh_mu, 4 );

      /* int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)
       * compute the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H
       * for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
       *
       * sh_mu = sqrt(gamma/2) * Mh_mu_nu A^nu, where here gamma = TAtom (as Sinv=1)
       */
      REAL8 norm_s = sqrt(TAtom / 2.0);
      if ( (gslstat = gsl_blas_dgemv (CblasNoTrans, norm_s, Mh_mu_nu, &A_Mu_view.vector, 0.0, &sh_mu_view.vector)) != 0 ) {
        XLALPrintError ( "%s: gsl_blas_dgemv(L * norm) failed: %s\n", fn, gsl_strerror (gslstat) );
        XLAL_ERROR_REAL8 ( fn, XLAL_EFAILED );
      }

      /* add this signal to the atoms, using the relation Fa,Fb <--> x_mu: see Eq.(72) in CFSv2-LIGO-T0900149-v2.pdf */
      REAL8 s1,s2,s3,s4;
      s1 = sh_mu[0];
      s2 = sh_mu[1];
      s3 = sh_mu[2];
      s4 = sh_mu[3];

      atoms->data[alpha].Fa_alpha.re +=   s1;
      atoms->data[alpha].Fa_alpha.im += - s3;
      atoms->data[alpha].Fb_alpha.re +=   s2;
      atoms->data[alpha].Fb_alpha.im += - s4;

    } /* for alpha < numAtoms */

  /* compute optimal SNR^2 expected for this signal,
   * using rho2 = A^mu M'_mu_nu A^nu = T/Sn( A' [A1^2+A3^2] + 2C' [A1A2 +A3A4] + B' [A2^2+A4^2])
   * NOTE: here we use the "effective" transient-CW antenna-pattern matrix M'_mu_nu
   */
  REAL8 A1 = A_Mu[0];
  REAL8 A2 = A_Mu[1];
  REAL8 A3 = A_Mu[2];
  REAL8 A4 = A_Mu[3];

  REAL8 rho2 = TAtom  * ( Ap * ( SQ(A1) + SQ(A3) ) + 2.0*Cp * ( A1*A2 + A3*A4 ) + Bp * ( SQ(A2) + SQ(A4) ) );

  /* return "effective" transient antenna-pattern matrix */
  M_mu_nu->Sinv_Tsft = TAtom;	/* everything here in units of Sn, so effectively Sn=1 */
  M_mu_nu->Ad = Ap;
  M_mu_nu->Bd = Bp;
  M_mu_nu->Cd = Cp;
  M_mu_nu->Dd = Ap * Bp - Cp * Cp;

  /* free memory */
  gsl_matrix_free ( Mh_mu_nu );

  /* return SNR^2 */
  return rho2;

} /* XLALAddSignalToFstatAtomVector() */


/** Add given signal s_mu = M_mu_nu A^nu within the given transient-window
 * to multi-IFO noise-atoms.
 *
 * RETURN: SNR^2 of the injected signal
 * and the effective AntennaPatternMatrix M_mu_nu for this signal.
 */
REAL8
XLALAddSignalToMultiFstatAtomVector ( MultiFstatAtomVector* multiAtoms,	 /**< [in/out] multi atoms vectors containing antenna-functions and possibly noise {Fa,Fb} */
                                      AntennaPatternMatrix *M_mu_nu,	 /**< [out] effective antenna-pattern matrix for the injected signal */
                                      const PulsarAmplitudeVect A_Mu, 	/**< [in] input canonical amplitude vector A^mu = {A1,A2,A3,A4} */
                                      transientWindow_t transientWindow /**< transient signal window */
                                      )
{
  const char *fn = __func__;

  /* check input consistency */
  if ( !multiAtoms || !multiAtoms->data ) {
    XLALPrintError ( "%s: Invalid NULL input 'multiAtoms'\n", fn );
    XLAL_ERROR_REAL8 ( fn, XLAL_EINVAL );
  }
  if ( !M_mu_nu ) {
    XLALPrintError ( "%s: Invalid NULL input 'M_mu_nu'\n", fn );
    XLAL_ERROR_REAL8 ( fn, XLAL_EINVAL );
  }

  UINT4 numDet = multiAtoms->length;
  UINT4 X;

  REAL8 rho2 = 0;
  (*M_mu_nu) = empty_AntennaPatternMatrix;

  for ( X=0; X < numDet; X ++ )
    {
      REAL8 rho2X;
      AntennaPatternMatrix M_mu_nu_X;
      rho2X = XLALAddSignalToFstatAtomVector ( multiAtoms->data[X], &M_mu_nu_X, A_Mu, transientWindow );
      if ( xlalErrno ) {
        XLALPrintError ("%s: XLALAddSignalToFstatAtomVector() failed.\n", fn );
        XLAL_ERROR_REAL8 ( fn, XLAL_EFUNC );
      }

      rho2 += rho2X;			/* multi-IFO SNR^2 = sum_X SNR_X^2 */
      M_mu_nu->Ad += M_mu_nu_X.Ad;	/* multi-IFO M_mu_nu = sum_X M_mu_nu_X */
      M_mu_nu->Bd += M_mu_nu_X.Bd;
      M_mu_nu->Cd += M_mu_nu_X.Cd;

      M_mu_nu->Sinv_Tsft += M_mu_nu_X.Sinv_Tsft;	/* noise adds harmonically 1/S = sum_X (1/S_X) */

    } /* for X < numDet */

  /* update sub-determinant */
  M_mu_nu->Dd = M_mu_nu->Ad * M_mu_nu->Bd - SQ(M_mu_nu->Cd);

  /* return SNR^2 */
  return rho2;

} /* XLALAddSignalToMultiFstatAtomVector() */


/** Load Ephemeris from ephemeris data-files
 */
EphemerisData *
XLALInitEphemeris (const CHAR *ephemYear )	/**< which years do we need? */
{
  const char *fn = __func__;

#define FNAME_LENGTH 1024
  CHAR EphemEarth[FNAME_LENGTH];	/* filename of earth-ephemeris data */
  CHAR EphemSun[FNAME_LENGTH];	/* filename of sun-ephemeris data */

  /* check input consistency */
  if ( !ephemYear ) {
    XLALPrintError ("%s: invalid NULL input for 'ephemYear'\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  snprintf(EphemEarth, FNAME_LENGTH, "earth%s.dat", ephemYear);
  snprintf(EphemSun, FNAME_LENGTH, "sun%s.dat",  ephemYear);

  EphemEarth[FNAME_LENGTH-1]=0;
  EphemSun[FNAME_LENGTH-1]=0;

  EphemerisData *edat;
  if ( (edat = XLALInitBarycenter ( EphemEarth, EphemSun)) == NULL ) {
    XLALPrintError ("%s: XLALInitBarycenter() failed.\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  /* return ephemeris */
  return edat;

} /* XLALInitEphemeris() */

/** Generates a multi-Fstat atoms vector for given parameters, drawing random parameters wherever required.
 *
 * Input: detector states, signal sky-pos (or allsky), amplitudes (or range), transient window range
 *
 */
MultiFstatAtomVector *
XLALSynthesizeTransientAtoms ( InjParams_t *injParams,		/**< [out] return summary of injected signal parameters (can be NULL) */
                               const ConfigVariables *cfg,	/**< [in] input params for transient atoms synthesis */
                               multiAMBuffer_t *multiAMBuffer	/**< buffer for AM-coefficients if re-using same skyposition (must be !=NULL) */
                               )
{
  const char *fn = __func__;

  /* check input */
  if ( !cfg || !cfg->rng || !multiAMBuffer ) {
    XLALPrintError ("%s: invalid NULL input in cfg, random-number generator cfg->rng, or multiAMBuffer\n", fn );
    XLAL_ERROR_NULL (fn, XLAL_EINVAL );
  }

  SkyPosition skypos;

  /* -----  if Alpha < 0 ==> draw skyposition isotropically from all-sky */
  if ( cfg->skypos.longitude < 0 )
    {
      skypos.longitude = gsl_ran_flat ( cfg->rng, 0, LAL_TWOPI );	/* alpha uniform in [ 0, 2pi ] */
      skypos.latitude  = acos ( gsl_ran_flat ( cfg->rng, -1, 1 ) ) - LAL_PI_2;	/* cos(delta) uniform in [ -1, 1 ] */
      skypos.system    = COORDINATESYSTEM_EQUATORIAL;
      /* never re-using buffered AM-coeffs here, as we randomly draw new skypositions */
      if ( multiAMBuffer->multiAM ) XLALDestroyMultiAMCoeffs ( multiAMBuffer->multiAM );
      multiAMBuffer->multiAM = NULL;
    } /* if random skypos to draw */
  else /* otherwise: re-use AM-coeffs if already computed, or initialize them if for the first time */
    {
      if ( multiAMBuffer->multiAM )
        if ( multiAMBuffer->skypos.longitude != cfg->skypos.longitude ||
             multiAMBuffer->skypos.latitude  != cfg->skypos.latitude  ||
             multiAMBuffer->skypos.system    != cfg->skypos.system )
          {
            XLALDestroyMultiAMCoeffs ( multiAMBuffer->multiAM );
            multiAMBuffer->multiAM = NULL;
          }
      skypos = cfg->skypos;
    } /* if single skypos given */

  /* ----- generate antenna-pattern functions for this sky-position */
  const MultiNoiseWeights *weights = NULL;	/* NULL = unit weights */
  if ( !multiAMBuffer->multiAM && (multiAMBuffer->multiAM = XLALComputeMultiAMCoeffs ( cfg->multiDetStates, weights, skypos )) == NULL ) {
    XLALPrintError ( "%s: XLALComputeMultiAMCoeffs() failed with xlalErrno = %d\n", fn, xlalErrno );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }
  multiAMBuffer->skypos = skypos; /* store buffered skyposition */

  /* ----- generate a pre-initialized F-stat atom vector containing only the antenna-pattern coefficients */
  MultiFstatAtomVector *multiAtoms;
  if ( (multiAtoms = XLALGenerateMultiFstatAtomVector ( cfg->multiTS, multiAMBuffer->multiAM )) == NULL ) {
    XLALPrintError ( "%s: XLALGenerateMultiFstatAtomVector() failed with xlalErrno = %d\n", fn, xlalErrno );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  /* ----- draw amplitude vector A^mu from given ranges in {h0, cosi, psi, phi0} */
  PulsarAmplitudeParams Amp;
  if ( cfg->AmpPrior.fixedSNR > 0 )	/* special treatment of fixed-SNR: use h0=1, later rescale signal */
    Amp.h0 = 1.0;
  else if ( cfg->AmpPrior.fixedSNR == 0 )/* same as setting h0 = 0 */
    Amp.h0 = 0;
  else					/* otherwise, draw from h0-prior */
    Amp.h0 = XLALDrawFromPDF1D ( cfg->AmpPrior.pdf_h0Nat, cfg->rng );

  Amp.cosi = XLALDrawFromPDF1D ( cfg->AmpPrior.pdf_cosi, cfg->rng );
  Amp.psi  = XLALDrawFromPDF1D ( cfg->AmpPrior.pdf_psi,  cfg->rng );
  Amp.phi0 = XLALDrawFromPDF1D ( cfg->AmpPrior.pdf_phi0, cfg->rng );

  /* convert amplitude params to 'canonical' vector coordinates */
  PulsarAmplitudeVect A_Mu;
  if ( XLALAmplitudeParams2Vect ( A_Mu, Amp ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALAmplitudeParams2Vect() failed with xlalErrno = %d\n", fn, xlalErrno );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  /* ----- draw transient-window parameters from given ranges using flat priors */
  transientWindow_t injectWindow = empty_transientWindow;
  injectWindow.type = cfg->transientInjectRange.type;
  if ( injectWindow.type != TRANSIENT_NONE )	/* nothing to be done if no window */
    {
      injectWindow.t0  = (UINT4) gsl_ran_flat ( cfg->rng, cfg->transientInjectRange.t0, cfg->transientInjectRange.t0 + cfg->transientInjectRange.t0Band );
      injectWindow.tau = (UINT4) gsl_ran_flat ( cfg->rng, cfg->transientInjectRange.tau, cfg->transientInjectRange.tau + cfg->transientInjectRange.tauBand );
    }

  /* ----- add transient signal to the Fstat atoms */
  AntennaPatternMatrix M_mu_nu;
  REAL8 rho2 = XLALAddSignalToMultiFstatAtomVector ( multiAtoms, &M_mu_nu, A_Mu, injectWindow );
  if ( xlalErrno ) {
    XLALPrintError ( "%s: XLALAddSignalToMultiFstatAtomVector() failed with xlalErrno = %d\n", fn, xlalErrno );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  /* ----- special signal rescaling if 1) fixedSNR OR 2) fixed rhohMax */
  REAL8 detM1o8 = sqrt ( M_mu_nu.Sinv_Tsft ) * pow ( M_mu_nu.Dd, 0.25 );	// (detM)^(1/8) = sqrt(Tsft/Sn) * (Dp)^(1/4)

  /* 1) if fixedSNR signal is requested: rescale everything to the desired SNR now */
  if ( (cfg->AmpPrior.fixedSNR > 0) || cfg->AmpPrior.fixRhohMax )
    {
      if ( (cfg->AmpPrior.fixedSNR > 0) && cfg->AmpPrior.fixRhohMax ) { /* double-check consistency: only one is allowed */
        XLALPrintError ("%s: Something went wrong: both [cfg->fixedSNR = %f > 0] and [cfg->fixedRhohMax==true] are not allowed!\n", fn, cfg->AmpPrior.fixedSNR );
        XLAL_ERROR_NULL ( fn, XLAL_EDOM );
      }

      REAL8 rescale = 1.0;

      if ( cfg->AmpPrior.fixedSNR > 0 )
        rescale = cfg->AmpPrior.fixedSNR / sqrt(rho2);	// rescale atoms by this factor, such that SNR = cfg->fixedSNR
      if ( cfg->AmpPrior.fixRhohMax )
        rescale = 1.0 / detM1o8;	// we drew h0 in [0, rhohMax], so we now need to rescale as h0Max = rhohMax/(detM)^(1/8)

      if ( XLALRescaleMultiFstatAtomVector ( multiAtoms, rescale ) != XLAL_SUCCESS ) {	      /* rescale atoms */
        XLALPrintError ( "%s: XLALRescaleMultiFstatAtomVector() failed with xlalErrno = %d\n", fn, xlalErrno );
        XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
      }

      Amp.h0 *= rescale;	      /* rescale amplitude-params for consistency */
      UINT4 i; for (i=0; i < 4; i ++) A_Mu[i] *= rescale;

      rho2 *= SQ(rescale);	      /* rescale reported optimal SNR */

    } /* if fixedSNR > 0 OR fixedRhohMax */

  /* ----- add noise to the Fstat atoms, unless --SignalOnly was specified */
  if ( !cfg->SignalOnly )
    if ( XLALAddNoiseToMultiFstatAtomVector ( multiAtoms, cfg->rng ) != XLAL_SUCCESS ) {
      XLALPrintError ("%s: XLALAddNoiseToMultiFstatAtomVector() failed with xlalErrno = %d\n", fn, xlalErrno );
      XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
    }

  /* ----- if requested: return all inject signal parameters */
  if ( injParams )
    {
      injParams->skypos = skypos;
      injParams->ampParams = Amp;
      UINT4 i; for (i=0; i < 4; i ++) injParams->ampVect[i] = A_Mu[i];
      injParams->M_mu_nu = M_mu_nu;
      injParams->transientWindow = injectWindow;
      injParams->SNR = sqrt(rho2);
      injParams->detM1o8 = detM1o8;
    } /* if injParams */


  return multiAtoms;

} /* XLALSynthesizeTransientAtoms() */

/** Rescale a given multi-Fstat atoms vector {Fa,Fb} by given scalar factor.
 * This is used to rescale signal atoms to desired fixed SNR.
 */
int
XLALRescaleMultiFstatAtomVector ( MultiFstatAtomVector* multiAtoms,	/**< [in/out] multi atoms vectors containing a signal in {Fa,Fb} to be rescaled */
                                  REAL8 rescale				/**< rescale factor: Fa' = rescale * Fa, and Fb'= rescale * Fb */
                                  )
{
  const char *fn = __func__;

  /* check input */
  if ( !multiAtoms ) {
    XLALPrintError ("%s: invalid NULL input 'multiAtoms'\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  UINT4 numDet = multiAtoms->length;
  UINT4 X;

  for ( X=0; X < numDet; X ++ )
    {
      FstatAtomVector *atoms = multiAtoms->data[X];
      UINT4 numAtoms = atoms->length;
      UINT4 i;
      for ( i=0; i < numAtoms; i ++ )
        {
          FstatAtom *thisAtom = &atoms->data[i];

          thisAtom->Fa_alpha.re *= rescale;
          thisAtom->Fa_alpha.im *= rescale;

          thisAtom->Fb_alpha.re *= rescale;
          thisAtom->Fb_alpha.im *= rescale;

        } /* for i < numAtoms */

    } /* for X < numDet */

  return XLAL_SUCCESS;

} /* XLALRescaleMultiFstatAtomVector() */


/** Write an injection-parameters structure to the given file-pointer,
 * adding one line with the injection parameters
 */
int
write_InjParams_to_fp ( FILE * fp,		/**< [in] file-pointer to output file */
                        const InjParams_t *par,	/**< [in] injection params to write. NULL means write header-comment line */
                        UINT4 dataStartGPS	/**< [in] data start-time in GPS seconds (used to turn window 't0' into offset from dataStartGPS) */
                        )
{
  const char *fn = __func__;

  /* input consistency */
  if ( ! fp ) {
    XLALPrintError ("%s: invalid NULL input 'fp'\n", fn);
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  int ret;
  /* if requested, write header-comment line */
  if ( par == NULL ) {
    ret = fprintf(fp, "%%%%Alpha Delta      SNR       h0   cosi    psi   phi0          A1       A2       A3       A4         Ad       Bd       Cd       Dd       t0[d]     tau[d] type   (detMp)^(1/8)\n");
    if ( ret < 0 ) {
      XLALPrintError ("%s: failed to fprintf() to given file-pointer 'fp'.\n", fn );
      XLAL_ERROR ( fn, XLAL_EIO );
    }

    return XLAL_SUCCESS;	/* we're done here */

  } /* if par == NULL */

  REAL8 t0_d = 1.0 * ( par->transientWindow.t0 - dataStartGPS ) / DAY24;
  REAL8 tau_d = 1.0 * par->transientWindow.tau / DAY24;

  /* if injParams given, output them to the file */
  ret = fprintf ( fp, " %5.3f %6.3f   %6.3f  %7.3g %6.3f %6.3f %6.3f    %8.3g %8.3g %8.3g %8.3g   %8.3g %8.3g %8.3g %8.3g    %8.3f  %8.3f    %1d   %8.3g\n",
                  par->skypos.longitude, par->skypos.latitude,						/* skypos */
                  par->SNR,										/* SNR */
                  par->ampParams.h0, par->ampParams.cosi, par->ampParams.psi, par->ampParams.phi0,	/* amplitude params {h0,cosi,psi,phi0}*/
                  par->ampVect[0], par->ampVect[1], par->ampVect[2], par->ampVect[3],			/* ampltiude vector A^mu */
                  par->M_mu_nu.Ad, par->M_mu_nu.Bd, par->M_mu_nu.Cd, par->M_mu_nu.Dd,			/* antenna-pattern matrix components */
                  t0_d, tau_d, par->transientWindow.type,		/* transient-window params */
                  par->detM1o8										/* rescale parameter (detMp)^(1/8) */
                  );
  if ( ret < 0 ) {
    XLALPrintError ("%s: failed to fprintf() to given file-pointer 'fp'.\n", fn );
    XLAL_ERROR ( fn, XLAL_EIO );
  }

 return XLAL_SUCCESS;

} /* write_InjParams_to_fp() */


/** Initialize amplitude-prior pdfs from the user-input
 */
int
XLALInitAmplitudePrior ( AmplitudePrior_t *AmpPrior, const UserInput_t *uvar )
{
  const char *fn = __func__;

  const UINT4 AmpPriorBins = 100;	// defines the binnning accuracy of our prior-pdfs

  /* consistency check */
  if ( !AmpPrior || !uvar ) {
    XLALPrintError ( "%s: invalid NULL input 'AmpPrior' or 'uvar'\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( AmpPrior->pdf_h0Nat || AmpPrior->pdf_cosi || AmpPrior->pdf_psi || AmpPrior->pdf_phi0 ) {
    XLALPrintError ("%s: AmplitudePriors must be set to NULL before calling this function!\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  /* first check that user only provided *one* method of determining the amplitude-prior range */
  UINT4 numSets = 0;
  if ( uvar->fixedh0Nat >= 0 ) numSets ++;
  if ( uvar->fixedSNR >= 0 ) numSets ++;
  if ( uvar->fixedh0NatMax >= 0 ) numSets ++ ;
  if ( uvar->fixedRhohMax >= 0 ) numSets ++;
  if ( numSets != 1 ) {
    XLALPrintError ("%s: Specify (>=0) exactly *ONE* amplitude-prior range of {fixedh0Nat, fixedSNR, fixedh0NatMax, fixedRhohMax}\n", fn);
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  /* ===== first pass: deal with all user-supplied fixed values ==> singular priors! ===== */

  /* ----- h0 ----- */
  if ( uvar->fixedh0Nat >= 0 )	/* fix h0Nat */
    if ( (AmpPrior->pdf_h0Nat = XLALCreateSingularPDF1D ( uvar->fixedh0Nat )) == NULL )
      XLAL_ERROR ( fn, XLAL_EFUNC );

  if ( uvar->fixedSNR >= 0 )   /* dummy-pdf, as signal will be computed with h0Nat=1 then rescaled to fixedSNR */
    if ( (AmpPrior->pdf_h0Nat = XLALCreateSingularPDF1D ( 1.0 )) == NULL )
      XLAL_ERROR ( fn, XLAL_EFUNC );

  AmpPrior->fixedSNR   = uvar->fixedSNR;
  AmpPrior->fixRhohMax = (uvar->fixedRhohMax >= 0);

  /* ----- cosi ----- */
  if ( XLALUserVarWasSet ( &uvar->cosi ) )
    if ( (AmpPrior->pdf_cosi = XLALCreateSingularPDF1D (  uvar->cosi )) == NULL )
      XLAL_ERROR ( fn, XLAL_EFUNC );
  /* ----- psi ----- */
  if ( XLALUserVarWasSet ( &uvar->psi ) )
    if ( (AmpPrior->pdf_psi = XLALCreateSingularPDF1D (  uvar->psi )) == NULL )
      XLAL_ERROR ( fn, XLAL_EFUNC );
  /* ----- phi0 ----- */
  if ( XLALUserVarWasSet ( &uvar->phi0 ) )
    if ( (AmpPrior->pdf_phi0 = XLALCreateSingularPDF1D (  uvar->phi0 )) == NULL )
      XLAL_ERROR ( fn, XLAL_EFUNC );


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
          XLAL_ERROR ( fn, XLAL_EFUNC );
      /* ----- cosi ----- */
      if ( AmpPrior->pdf_cosi == NULL )
        if ( (AmpPrior->pdf_cosi = XLALCreateUniformPDF1D ( -1.0, 1.0 )) == NULL )
          XLAL_ERROR ( fn, XLAL_EFUNC );
      /* ----- psi ----- */
      if ( AmpPrior->pdf_psi == NULL )
        if ( (AmpPrior->pdf_psi = XLALCreateUniformPDF1D ( -LAL_PI_4, LAL_PI_4 )) == NULL )
          XLAL_ERROR ( fn, XLAL_EFUNC );
      /* ----- phi0 ----- */
      if ( AmpPrior->pdf_phi0 == NULL )
        if ( (AmpPrior->pdf_phi0 = XLALCreateUniformPDF1D ( 0, LAL_TWOPI )) == NULL )
          XLAL_ERROR ( fn, XLAL_EFUNC );

      break;

    case AMP_PRIOR_TYPE_CANONICAL:
      /* ----- pdf(h0) ~ h0^3 ----- */
      if ( AmpPrior->pdf_h0Nat == NULL )
        {
          UINT4 i;
          pdf1D_t *pdf;
          if ( ( pdf = XLALCreateDiscretePDF1D ( 0, h0NatMax, AmpPriorBins )) == NULL )
            XLAL_ERROR ( fn, XLAL_EFUNC );

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
            XLAL_ERROR ( fn, XLAL_EFUNC );

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
          XLAL_ERROR ( fn, XLAL_EFUNC );
      /* ----- phi0 ----- */
      if ( AmpPrior->pdf_phi0 == NULL )
        if ( (AmpPrior->pdf_phi0 = XLALCreateUniformPDF1D ( 0, LAL_TWOPI )) == NULL )
          XLAL_ERROR ( fn, XLAL_EFUNC );

      break;

    default:
      XLALPrintError ("%s: something went wrong ... unknown priorType = %d\n", fn, uvar->AmpPriorType );
      XLAL_ERROR ( fn, XLAL_EINVAL );
      break;

    } // switch( uvar->AmpPriorType )

  return XLAL_SUCCESS;

} /* XLALInitAmplitudePrior() */
