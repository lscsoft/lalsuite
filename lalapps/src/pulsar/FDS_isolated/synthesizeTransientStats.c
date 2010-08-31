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

#include "../transientCW_utils.h"

/* GSL includes */
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_bessel.h>


/* LAL-includes */
#include <lal/SkyCoordinates.h>
#include <lal/AVFactories.h>
#include <lal/LALInitBarycenter.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lal/ComputeFstat.h>


#include <lalapps.h>

/*---------- DEFINES ----------*/
#define EPHEM_YEARS  "05-09"	/**< default range, covering S5: override with --ephemYear */
#define SQ(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUAD(x) ((x)*(x)*(x)*(x))

/** Enumeration of allowed amplitude-prior types
 */
typedef enum
  {
    AMP_PRIOR_PHYSICAL = 0,	/**< 'physical' priors: isotropic pdf{cosi,psi,phi0} AND flat pdf(h0) */
    AMP_PRIOR_CANONICAL,	/**< 'canonical' priors: uniform in A^mu up to h_max */
    AMP_PRIOR_LAST
  } AmpPriorType_t;

/** Encode a probability distribution pdf(x) over an interval [xMin, xMin+xBand].
 * Set xBand=0 for the special case of one certain value with p(x_0)=1 (pdf is not used in this case).
 * The optional field 'sampling' allows to use gsl_ran_discrete() to draw samples from that distribution.
 */
typedef struct
{
  REAL8 xMin;			/**< lower bound on random-variable domain [xMin, xMin+xBand] */
  REAL8 xBand;			/**< size of random-variable domain, can be 0 */
  REAL8Vector *prob;		/**< vector of pdf(x_i) over x-domain. NULL means uniform pdf over [xMin,xMax] */
  BOOLEAN normalized;		/**< true if the pdf(x) is normalized: sum_i pdf(x_i) dx = 1 */
  gsl_ran_discrete_t *sampling;	/**< optional: sampling distribution input for drawing samples using gsl_ran_discrete() */
} pdf1D_t;


/** Signal (amplitude) parameter ranges
 */
typedef struct {
  pdf1D_t pdf_h0Nat;	/**< pdf for h0/sqrt{Sn} */
  REAL8 fixedSNR;	/**< alternative: if > 0 => fix SNR of injected signals */
  pdf1D_t pdf_cosi;	/**< pdf(cosi) */
  pdf1D_t pdf_psi;	/**< pdf(psi) */
  pdf1D_t pdf_phi0;	/**< pdf(phi0) */
} AmplitudePrior_t;

/** Configuration settings required for and defining a coherent pulsar search.
 * These are 'pre-processed' settings, which have been derived from the user-input.
 */
typedef struct {
  AmplitudePrior_t AmpPrior;	/**< amplitude-parameter priors to draw signals from */
  SkyPosition skypos;		/**< (Alpha,Delta,system). Use Alpha < 0 to signal 'allsky' */

  transientWindowRange_t transientSearchRange;	/**< transient-window range for the search (flat priors) */
  transientWindowRange_t transientInjectRange;	/**< transient-window range for injections (flat priors) */

  MultiDetectorStateSeries *multiDetStates;	/**< multi-detector state series covering observation time */
  MultiLIGOTimeGPSVector *multiTS;		/**< corresponding timestamps vector for convenience */

  gsl_rng *rng;			/**< gsl random-number generator */
  CHAR *logString;		/**< logstring for file-output, containing cmdline-options + code VCS version info */

} ConfigVariables;

/* struct for buffering of AM-coeffs, if signal for same sky-position is injected */
typedef struct {
  SkyPosition skypos;		/**< sky-position for which we have AM-coeffs computed already */
  MultiAMCoeffs *multiAM;;	/**< pre-computed AM-coeffs for skypos */
} multiAMBuffer_t;


/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

/* ----- User-variables: can be set from config-file or command-line */
typedef struct {
  BOOLEAN help;		/**< trigger output of help string */

  /* amplitude parameters + ranges */
  REAL8 h0Nat;		/**< GW amplitude: h0/sqrt(Sn) */
  REAL8 h0NatBand;	/**< Band to draw GW amplitude [h0Nat, h0Nat + h0NatBand ] */
  REAL8 fixedSNR;	/**< Alternative: fixe the optimal SNR of the injected signals */

  REAL8 cosi;		/**< cos(inclination angle). If not set: randomize within [-1,1] */
  REAL8 psi;		/**< polarization angle psi. If not set: randomize within [-pi/4,pi/4] */
  REAL8 phi0;		/**< initial GW phase phi_0. If not set: randomize within [0, 2pi] */
  INT4 AmpPriorType;	/**< enumeration of types of amplitude-priors: 0=physical, 1=canonical */

  /* Doppler parameters */
  REAL8 Alpha;		/**< skyposition Alpha (RA) in radians */
  REAL8 Delta;		/**< skyposition Delta (Dec) in radians */

  /* transient window ranges: for injection ... */
  CHAR *injectWindow_type; 	/**< Type of transient window to inject ('none', 'rect', 'exp'). */
  INT4  injectWindow_t0;	/**< Earliest GPS start-time of transient window to inject, in seconds */
  INT4  injectWindow_t0Band;	/**< Range of GPS start-time of transient window to inject, in seconds */
  REAL8 injectWindow_tauDays;	/**< Shortest transient-window timescale to inject, in days */
  REAL8 injectWindow_tauDaysBand; /**< Range of transient-window timescale to inject, in days */
  /* ... and for search */
  CHAR *searchWindow_type; 	/**< Type of transient window to search with ('none', 'rect', 'exp'). */
  INT4  searchWindow_t0;		/**< Earliest GPS start-time of transient window to search, in seconds */
  INT4  searchWindow_t0Band;	/**< Range of GPS start-time of transient window to search, in seconds */
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
  BOOLEAN SignalOnly;	/**< dont generate noise-draws: will result in non-random 'signal only' values of F and B */

  BOOLEAN useFReg;	/**< use 'regularized' Fstat (1/D)*e^F for marginalization, or 'standard' e^F */

  CHAR *ephemYear;	/**< date-range string on ephemeris-files to use */

  BOOLEAN version;	/**< output version-info */

} UserInput_t;




/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);

int XLALInitUserVars ( UserInput_t *uvar );
int XLALInitCode ( ConfigVariables *cfg, const UserInput_t *uvar );
EphemerisData * XLALInitEphemeris (const CHAR *ephemYear );

/* exportable API */
int XLALDrawCorrelatedNoise ( gsl_vector *n_mu, const gsl_matrix *L, gsl_rng * rng );

FstatAtomVector* XLALGenerateFstatAtomVector ( const LIGOTimeGPSVector *TS, const AMCoeffs *amcoeffs );
MultiFstatAtomVector* XLALGenerateMultiFstatAtomVector ( const  MultiLIGOTimeGPSVector *multiTS, const MultiAMCoeffs *multiAM );

int XLALAddNoiseToFstatAtomVector ( FstatAtomVector *atoms, gsl_rng * rng );
int XLALAddNoiseToMultiFstatAtomVector ( MultiFstatAtomVector *multiAtoms, gsl_rng * rng );

int XLALAddSignalToFstatAtomVector ( FstatAtomVector* atoms, REAL8 *rho2, const gsl_vector *A_Mu, transientWindow_t transientWindow );
int XLALAddSignalToMultiFstatAtomVector ( MultiFstatAtomVector* multiAtoms, REAL8 *rho2, const gsl_vector *A_Mu, transientWindow_t transientWindow );

MultiFstatAtomVector *XLALSynthesizeTransientAtoms ( ConfigVariables *cfg, BOOLEAN SignalOnly, multiAMBuffer_t *multiAMBuffer );

int XLALRescaleMultiFstatAtomVector ( MultiFstatAtomVector* multiAtoms,	REAL8 rescale );


int XLALInitAmplitudePrior ( AmplitudePrior_t *AmpPrior, UINT4 Npoints, AmpPriorType_t priorType );
REAL8 XLALDrawFromPDF ( const pdf1D_t *probDist, const gsl_rng *rng );


/*---------- empty initializers ---------- */
ConfigVariables empty_ConfigVariables;
UserInput_t empty_UserInput;
multiAMBuffer_t empty_multiAMBuffer;


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
    return 1;
  }

  FILE *fpTransientStats = NULL;
  /* ----- prepare stats output ----- */
  if ( uvar.outputStats )
    {
      if ( (fpTransientStats = fopen (uvar.outputStats, "wb")) == NULL)
	{
	  LogPrintf (LOG_CRITICAL, "Error opening file '%s' for writing..\n\n", uvar.outputStats );
	  return 1;
	}
      fprintf (fpTransientStats, "%s", cfg.logString );		/* write search log comment */
      write_TransientCandidate_to_fp ( fpTransientStats, NULL );	/* write header-line comment */
    }

  /* ----- main MC loop over numDraws trials ---------- */
  multiAMBuffer_t multiAMBuffer = empty_multiAMBuffer;	  /* prepare AM-buffer */
  INT4 i;

  for ( i=0; i < uvar.numDraws; i ++ )
    {
      /* generate signal random draws from ranges and generate Fstat atoms */
      MultiFstatAtomVector *multiAtoms;
      if ( (multiAtoms = XLALSynthesizeTransientAtoms ( &cfg, uvar.SignalOnly, &multiAMBuffer )) == NULL ) {
        LogPrintf ( LOG_CRITICAL, "%s: XLALSynthesizeTransientAtoms() failed with xlalErrno = %d\n", fn, xlalErrno );
        return 1;
      }

      /* compute transient-Bstat search statistic on these atoms */
      TransientCandidate_t cand = empty_TransientCandidate;
      if ( XLALComputeTransientBstat ( &cand, multiAtoms,  cfg.transientSearchRange, uvar.useFReg ) != XLAL_SUCCESS ) {
        LogPrintf ( LOG_CRITICAL, "%s: XLALComputeTransientBstat() failed with xlalErrno = %d\n", fn, xlalErrno );
        return 1;
      }

      /* if requested, also compute Ftotal over full data-span */
      if ( uvar.computeFtotal )
        {
          TransientCandidate_t candTotal;
          transientWindowRange_t winRangeAll = empty_transientWindowRange;
          winRangeAll.type = TRANSIENT_NONE;	/* window 'none' will simply cover all the data with 1 F-stat calculation */
          if ( XLALComputeTransientBstat ( &candTotal, multiAtoms,  winRangeAll, uvar.useFReg ) != XLAL_SUCCESS ) {
            LogPrintf ( LOG_CRITICAL, "%s: XLALComputeTransientBstat() failed for totalFstat (winRangeAll) with xlalErrno = %d\n", fn, xlalErrno );
            return 1;
          }
          /* we only carry over twoFtotal = maxTwoF from this single-Fstat calculation */
          cand.twoFtotal = candTotal.maxTwoF;
        } /* if computeFtotal */

      /* if requested, output atoms-vector into file */
      if ( uvar.outputAtoms )
        {
          FILE *fpAtoms;
          char *fnameAtoms;
          UINT4 len = strlen ( uvar.outputAtoms ) + 20;
          if ( (fnameAtoms = XLALCalloc ( 1, len )) == NULL ) {
            XLALPrintError ("%s: failed to XLALCalloc ( 1, %d )\n", fn, len );
            return 1;
          }
          sprintf ( fnameAtoms, "%s_%04d_of_%04d.dat", uvar.outputAtoms, i + 1, uvar.numDraws );

          if ( ( fpAtoms = fopen ( fnameAtoms, "wb" )) == NULL ) {
            XLALPrintError ("%s: failed to open atoms-output file '%s' for writing.\n", fn, fnameAtoms );
            return 1;
          }
	  fprintf ( fpAtoms, "%s", cfg.logString );	/* output header info */

	  if ( write_MultiFstatAtoms_to_fp ( fpAtoms, multiAtoms ) != XLAL_SUCCESS ) {
            XLALPrintError ("%s: failed to write atoms to output file '%s'. xlalErrno = %d\n", fn, fnameAtoms, xlalErrno );
            return 1;
          }

          XLALFree ( fnameAtoms );
	  fclose (fpAtoms);
        } /* if outputAtoms */


      /* free atoms */
      XLALDestroyMultiFstatAtomVector ( multiAtoms );

      /* add info on current transient-CW candidate */
      cand.doppler.Alpha = multiAMBuffer.skypos.longitude;
      cand.doppler.Delta = multiAMBuffer.skypos.latitude;

      if ( uvar.SignalOnly )
        cand.maxTwoF += 4;

      if ( fpTransientStats && write_TransientCandidate_to_fp ( fpTransientStats, &cand ) != XLAL_SUCCESS ) {
        LogPrintf ( LOG_CRITICAL, "%s: write_TransientCandidate_to_fp() failed.\n", fn );
        return 1;
      }

    } /* for i < numDraws */

  /* ----- free memory ---------- */
  if ( fpTransientStats) fclose ( fpTransientStats );
  XLALDestroyMultiDetectorStateSeries ( cfg.multiDetStates );
  XLALDestroyMultiTimestamps ( cfg.multiTS );
  XLALDestroyMultiAMCoeffs ( multiAMBuffer.multiAM );
  if ( cfg.transientSearchRange.exp_buffer )
    gsl_matrix_free ( cfg.transientSearchRange.exp_buffer );

  if ( cfg.logString ) XLALFree ( cfg.logString );
  gsl_rng_free ( cfg.rng );

  XLALDestroyUserVars();

  /* did we forget anything ? */
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

  uvar->h0Nat = 1;
  uvar->h0NatBand = 0;

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
  REAL8 tauMax = ( uvar->injectWindow_tauDays +  uvar->injectWindow_tauDaysBand ) * DAY24;
  /* default window-ranges are t0 in [dataStartTime, dataStartTime - 3 * tauMax] */
  uvar->injectWindow_t0     = uvar->dataStartGPS;
  uvar->injectWindow_t0Band = fmax ( 0.0, uvar->dataDuration - TRANSIENT_EXP_EFOLDING * tauMax );	/* make sure it's >= 0 */
  /* search-windows by default identical to inject-windows */
  uvar->searchWindow_t0 = uvar->injectWindow_t0;
  uvar->searchWindow_t0Band = uvar->injectWindow_t0Band;
  uvar->searchWindow_tauDays = uvar->injectWindow_tauDays;
  uvar->searchWindow_tauDaysBand = uvar->injectWindow_tauDaysBand;

  uvar->searchWindow_dt0  = uvar->TAtom;
  uvar->searchWindow_dtau = uvar->TAtom;

  /* register all our user-variables */
  XLALregBOOLUserStruct ( help, 		'h',     UVAR_HELP, "Print this message");

  /* signal Doppler parameters */
  XLALregREALUserStruct ( Alpha, 		'a', UVAR_OPTIONAL, "Sky position alpha (equatorial coordinates) in radians [Default: allsky]");
  XLALregREALUserStruct ( Delta, 		'd', UVAR_OPTIONAL, "Sky position delta (equatorial coordinates) in radians [Default: allsky]");

  /* signal amplitude parameters */
  XLALregREALUserStruct ( h0Nat,		 0, UVAR_OPTIONAL, "GW amplitude h0 measured in units of noise-PSD sqrt(Sn)");
  XLALregREALUserStruct ( h0NatBand,		 0, UVAR_OPTIONAL, "Draw GW 'amplitude' randomly from [h0Nat, h0Nat+h0NatBand]");
  XLALregREALUserStruct ( fixedSNR,		 0, UVAR_OPTIONAL, "Alternative: fixe the optimal SNR of the injected signals");

  XLALregREALUserStruct ( cosi,			'i', UVAR_OPTIONAL, "cos(inclination angle). If not set: randomize within [-1,1].");
  XLALregREALUserStruct ( psi,			 0,  UVAR_OPTIONAL, "polarization angle psi. If not set: randomize within [-pi/4,pi/4].");
  XLALregREALUserStruct ( phi0,		 	 0,  UVAR_OPTIONAL, "initial GW phase phi_0. If not set: randomize within [0, 2pi]");

  XLALregINTUserStruct ( AmpPriorType,	 	 0,  UVAR_OPTIONAL, "Enumeration of types of amplitude-priors: 0=physical, 1=canonical");

  XLALregSTRINGUserStruct ( IFO,	        'I', UVAR_OPTIONAL, "Detector: 'G1','L1','H1,'H2', 'V1', ... ");
  XLALregINTUserStruct ( dataStartGPS,	 	 0,  UVAR_OPTIONAL, "data start-time in GPS seconds");
  XLALregINTUserStruct ( dataDuration,	 	 0,  UVAR_OPTIONAL, "data-span to generate (in seconds)");

  /* transient window ranges: for injection ... */
  XLALregSTRINGUserStruct( injectWindow_type,    0, UVAR_OPTIONAL, "Type of transient window to inject ('none', 'rect', 'exp')");
  XLALregREALUserStruct  ( injectWindow_tauDays, 0, UVAR_OPTIONAL, "Shortest transient-window timescale to inject, in days");
  XLALregREALUserStruct  ( injectWindow_tauDaysBand,0,UVAR_OPTIONAL,"Range of transient-window timescale to inject, in days");
  XLALregINTUserStruct   ( injectWindow_t0, 	 0, UVAR_OPTIONAL, "Earliest GPS start-time of transient window to inject, in seconds [dataStartGPS]");
  XLALregINTUserStruct   ( injectWindow_t0Band,  0, UVAR_OPTIONAL, "Range of GPS start-time of transient window to inject, in seconds [dataDuration-3*tauMax]");
  /* ... and for search */
  XLALregSTRINGUserStruct( searchWindow_type,    0, UVAR_OPTIONAL, "Type of transient window to search with ('none', 'rect', 'exp')");
  XLALregREALUserStruct  ( searchWindow_tauDays, 0, UVAR_OPTIONAL, "Shortest transient-window timescale to search, in days");
  XLALregREALUserStruct  ( searchWindow_tauDaysBand,0,UVAR_OPTIONAL, "Range of transient-window timescale to search, in days");
  XLALregINTUserStruct   ( searchWindow_dtau, 	 0, UVAR_OPTIONAL, "Step-size for search/marginalization over transient-window timescale, in seconds [Default:TAtom]");
  XLALregINTUserStruct   ( searchWindow_t0,      0, UVAR_OPTIONAL, "Earliest GPS start-time of transient window to search, in seconds [dataStartGPS]");
  XLALregINTUserStruct   ( searchWindow_t0Band,  0, UVAR_OPTIONAL, "Range of GPS start-time of transient window to search, in seconds [dataDuration-3*tauMax]");
  XLALregINTUserStruct   ( searchWindow_dt0, 	 0, UVAR_OPTIONAL, "Step-size for search/marginalization over transient-window start-time, in seconds [Default:TAtom]");

  /* misc params */
  XLALregBOOLUserStruct ( computeFtotal,	 0, UVAR_OPTIONAL, "Also compute 'total' F-statistic over the full data-span" );

  XLALregINTUserStruct ( numDraws,		'N', UVAR_OPTIONAL, "Number of random 'draws' to simulate");

  XLALregSTRINGUserStruct ( outputStats,	'o', UVAR_OPTIONAL, "Output file containing 'numDraws' random draws of stats");
  XLALregSTRINGUserStruct ( outputAtoms,	 0,  UVAR_OPTIONAL,  "Output F-statistic atoms into a file with this basename");
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

  /* ----- parse user-input on signal amplitude-paramters + ranges ----- */
  /* skypos */
  cfg->skypos.longitude = uvar->Alpha;	/* Alpha < 0 indicates 'allsky' */
  cfg->skypos.latitude  = uvar->Delta;
  cfg->skypos.system = COORDINATESYSTEM_EQUATORIAL;

  /* amplitude-params: setup pdf ranges, then initialize the priors */
  cfg->AmpPrior.pdf_h0Nat.xMin = uvar->h0Nat;
  if ( uvar->h0NatBand < 0 ) {
    XLALPrintError ("%s: only non-negatives values >= 0 for 'h0NatBand' allowed, got '%g'\n", fn, uvar->h0NatBand );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  cfg->AmpPrior.pdf_h0Nat.xBand = uvar->h0NatBand;

  /* implict ranges on cosi, psi and phi0 if not specified by user */
  if ( XLALUserVarWasSet ( &uvar->cosi ) ) {
    cfg->AmpPrior.pdf_cosi.xMin  = uvar->cosi;
    cfg->AmpPrior.pdf_cosi.xBand = 0;
  } else {
    cfg->AmpPrior.pdf_cosi.xMin  = -1;
    cfg->AmpPrior.pdf_cosi.xBand = 2;
  }
  if ( XLALUserVarWasSet ( &uvar->psi ) ) {
    cfg->AmpPrior.pdf_psi.xMin = uvar->psi;
    cfg->AmpPrior.pdf_psi.xBand= 0;
  } else {
    cfg->AmpPrior.pdf_psi.xMin  = - LAL_PI_4;
    cfg->AmpPrior.pdf_psi.xBand =   LAL_PI_2;
  }
  if ( XLALUserVarWasSet ( &uvar->phi0 ) ) {
    cfg->AmpPrior.pdf_phi0.xMin  = uvar->phi0;
    cfg->AmpPrior.pdf_phi0.xBand = 0;
  } else {
    cfg->AmpPrior.pdf_phi0.xMin  = 0;
    cfg->AmpPrior.pdf_phi0.xBand = LAL_TWOPI;
  }
  UINT4 priorSampling = 100;
  if ( XLALInitAmplitudePrior ( &cfg->AmpPrior, priorSampling, uvar->AmpPriorType ) != XLAL_SUCCESS ) {
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }

  /* ----- initialize random-number generator ----- */
  /* read out environment variables GSL_RNG_xxx
   * GSL_RNG_SEED: use to set random seed: defult = 0
   * GSL_RNG_TYPE: type of random-number generator to use: default = 'mt19937'
   */
  gsl_rng_env_setup ();
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
    if ( XLALUserVarWasSet ( &uvar->injectWindow_t0 ) || XLALUserVarWasSet ( &uvar->injectWindow_t0Band ) ||
         XLALUserVarWasSet ( &uvar->injectWindow_tauDays ) || XLALUserVarWasSet ( &uvar->injectWindow_tauDaysBand ) ) {
      XLALPrintError ("%s: ERROR: injectWindow_type == NONE, but window-parameters were set! Use a different window-type!\n", fn );
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }

  if ( uvar->injectWindow_t0Band < 0 || uvar->injectWindow_tauDaysBand < 0 ) {
    XLALPrintError ("%s: only positive t0/tau window injection bands allowed (%d, %f)\n", fn, uvar->injectWindow_t0Band, uvar->injectWindow_tauDaysBand );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  /* apply correct defaults if unset: t0=dataStart, t0Band=dataDuration-3*tauMax */
  if ( !XLALUserVarWasSet ( &uvar->injectWindow_t0 ) )
    InjectRange.t0 = uvar->dataStartGPS;
  else
    InjectRange.t0      = uvar->injectWindow_t0;

  REAL8 tauMax = ( uvar->injectWindow_tauDays +  uvar->injectWindow_tauDaysBand ) * DAY24;
  if ( !XLALUserVarWasSet (&uvar->injectWindow_t0Band ) )
    InjectRange.t0Band  = fmax ( 0.0, uvar->dataDuration - TRANSIENT_EXP_EFOLDING * tauMax ); 	/* make sure it's >= 0 */
  else
    InjectRange.t0Band  = uvar->injectWindow_t0Band;

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
  /* make sure user doesn't set window=none but sets window-parameters => indicates she didn't mean 'none' */
  if ( SearchRange.type == TRANSIENT_NONE )
    if ( XLALUserVarWasSet ( &uvar->searchWindow_t0 ) || XLALUserVarWasSet ( &uvar->searchWindow_t0Band ) ||
         XLALUserVarWasSet ( &uvar->searchWindow_tauDays ) || XLALUserVarWasSet ( &uvar->searchWindow_tauDaysBand ) ) {
      XLALPrintError ("%s: ERROR: searchWindow_type == NONE, but window-parameters were set! Use a different window-type!\n", fn );
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }

  if (   uvar->searchWindow_t0Band < 0 || uvar->searchWindow_tauDaysBand < 0 ) {
    XLALPrintError ("%s: only positive t0/tau window injection bands allowed (%d, %f)\n", fn, uvar->searchWindow_t0Band, uvar->searchWindow_tauDaysBand );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  /* apply correct defaults if unset: use inect window */
  if ( !XLALUserVarWasSet ( &uvar->searchWindow_t0 ) )
    SearchRange.t0      = InjectRange.t0;
  else
    SearchRange.t0      = uvar->searchWindow_t0;
  if ( !XLALUserVarWasSet ( &uvar->searchWindow_t0Band ) )
    SearchRange.t0Band = InjectRange.t0Band;
  else
    SearchRange.t0Band  = uvar->searchWindow_t0Band;
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

  /* if exponential window --> prebuffer the window-values */
  if ( SearchRange.type == TRANSIENT_EXPONENTIAL )
    if ( XLALFillExpWindowBuffer ( &SearchRange ) != XLAL_SUCCESS ) {
      XLALPrintError ("%s: XLALFillExpWindowBuffer() failed.\n", fn );
      XLAL_ERROR ( fn, XLAL_EFUNC );
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
XLALDrawCorrelatedNoise ( gsl_vector *n_mu,		/**< [out] pre-allocated 4-vector of noise-components {n_mu}, with correlation L * L^T */
                          const gsl_matrix *L,		/**< [in] correlator matrix to get n_mu = L_mu_nu * norm_nu from 4 uncorr. unit variates norm_nu */
                          gsl_rng * rng			/**< gsl random-number generator */
                          )
{
  const CHAR *fn = __func__;

  static gsl_vector *normal = NULL;    /* keep a static copy to avoid repeated malloc/free cycles */
  int gslstat;

  /* ----- check input arguments ----- */
  if ( !n_mu || n_mu->size != 4 ) {
    XLALPrintError ( "%s: Invalid input, n_mu must be pre-allocate 4-vector", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( !L || (L->size1 != 4) || (L->size2 != 4) ) {
    XLALPrintError ( "%s: Invalid correlator matrix, n_mu must be pre-allocate 4x4 matrix", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( !rng ) {
    XLALPrintError ("%s: invalid NULL input as gsl random-number generator!\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }


  /* ----- generate 4 normal-distributed, uncorrelated random numbers ----- */
  if ( !normal && (normal = gsl_vector_alloc ( 4 ) ) == NULL) {
    XLALPrintError ( "%s: gsl_vector_calloc(4) failed\n", fn );
    XLAL_ERROR ( fn, XLAL_ENOMEM );
  }

  gsl_vector_set (normal, 0,  gsl_ran_gaussian ( rng, 1.0 ) );
  gsl_vector_set (normal, 1,  gsl_ran_gaussian ( rng, 1.0 ) );
  gsl_vector_set (normal, 2,  gsl_ran_gaussian ( rng, 1.0 ) );
  gsl_vector_set (normal, 3,  gsl_ran_gaussian ( rng, 1.0 ) );


  /* use four normal-variates {norm_nu} with correlator matrix L to get: n_mu = L_{mu nu} norm_nu
   * which gives {n_\mu} satisfying cov(n_mu,n_nu) = (L L^T)_{mu nu} = M_{mu nu}
   */

  /* int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)
   * compute the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H
   * for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
   */
  if ( (gslstat = gsl_blas_dgemv (CblasNoTrans, 1.0, L, normal, 0.0, n_mu)) != 0 ) {
    XLALPrintError ( "%s: gsl_blas_dgemv(L * norm) failed: %s\n", fn, gsl_strerror (gslstat) );
    XLAL_ERROR ( fn, XLAL_EFAILED );
  }

  /* ---------- free memory ---------- */
  // gsl_vector_free ( normal ); /* we keep this statically as it's always a 4D vector */

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
  gsl_vector *n_mu;
  if ( (n_mu = gsl_vector_calloc ( 4 ) ) == NULL ) {
    XLALPrintError ("%s: gsl_vector_calloc ( 4 ) failed.\n", fn );
    gsl_matrix_free ( Lcor );
    XLAL_ERROR ( fn, XLAL_ENOMEM );
  }

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
      x1 = gsl_vector_get ( n_mu, 0 );
      x2 = gsl_vector_get ( n_mu, 1 );
      x3 = gsl_vector_get ( n_mu, 2 );
      x4 = gsl_vector_get ( n_mu, 3 );

      /* add this to Fstat-atom */
      /* relation Fa,Fb <--> x_mu: see Eq.(72) in CFSv2-LIGO-T0900149-v2.pdf */
      atoms->data[alpha].Fa_alpha.re +=   x1;
      atoms->data[alpha].Fa_alpha.im += - x3;
      atoms->data[alpha].Fb_alpha.re +=   x2;
      atoms->data[alpha].Fb_alpha.im += - x4;

    } /* for i < numAtoms */

  /* free internal memory */
  gsl_vector_free ( n_mu );
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


/** Add given signal s_mu = M_mu_nu A^nu within the given transient-window
 * to noise-atoms
 */
int
XLALAddSignalToFstatAtomVector ( FstatAtomVector* atoms,	 /**< [in/out] atoms vectors containing antenna-functions and possibly noise {Fa,Fb} */
                                 REAL8 *rho2,			 /**< [out] expected optimal SNR^2 for this signal */
                                 const gsl_vector *A_Mu,	 /**< [in] input canonical amplitude vector A^mu = {A1,A2,A3,A4} */
                                 transientWindow_t transientWindow /**< transient signal window */
                                 )
{
  const char *fn = __func__;
  int gslstat;

  /* check input consistency */
  if ( !atoms || !atoms->data ) {
    XLALPrintError ( "%s: Invalid NULL input 'atoms'\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( !rho2 ) {
    XLALPrintError ( "%s: Invalid NULL input 'rho2'\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( !A_Mu || (A_Mu->size != 4) ) {
    XLALPrintError ( "%s: Invalid input vector A_Mu: must be allocated 4D\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  /* prepare transient-window support */
  UINT4 t0, t1;
  if ( XLALGetTransientWindowTimespan ( &t0, &t1, transientWindow ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALGetTransientWindowTimespan() failed.\n", fn );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }

  /* prepare gsl-matrix for Mh_mu_nu = [ a^2, a*b ; a*b , b^2 ] */
  gsl_matrix *Mh_mu_nu;
  if ( (Mh_mu_nu = gsl_matrix_calloc ( 4, 4 )) == NULL ) {
    XLALPrintError ("%s: gsl_matrix_calloc ( 4, 4 ) failed.\n", fn );
    XLAL_ERROR ( fn, XLAL_ENOMEM );
  }
  /* prepare placeholder for 4-vector sh_mu */
  gsl_vector *sh_mu;
  if ( (sh_mu = gsl_vector_calloc ( 4 ) ) == NULL ) {
    XLALPrintError ("%s: gsl_vector_calloc ( 4 ) failed.\n", fn );
    gsl_matrix_free ( Mh_mu_nu );
    XLAL_ERROR ( fn, XLAL_ENOMEM );
  }

  REAL8 TAtom = atoms->TAtom;
  UINT4 numAtoms = atoms->length;
  UINT4 alpha;
  REAL8 Ad = 0, Bd = 0, Cd = 0;
  for ( alpha=0; alpha < numAtoms; alpha ++ )
    {
      UINT4 ti = atoms->data[alpha].timestamp;
      REAL8 win = XLALGetTransientWindowValue ( ti, t0, t1, transientWindow.tau, transientWindow.type );

      if ( win == 0 )
        continue;

      /* compute sh_mu = sqrt(gamma/2) * Mh_mu_nu A^nu, where Mh_mu_nu is now just
       * the per-atom block matrix [a^2,  ab; ab, b^2 ]
       * where Sn=1, so gamma = Sinv*TAtom = TAtom
       */
      REAL8 norm = win * win;
      REAL8 a2 = norm * atoms->data[alpha].a2_alpha;
      REAL8 b2 = norm * atoms->data[alpha].b2_alpha;
      REAL8 ab = norm * atoms->data[alpha].ab_alpha;

      Ad += a2;
      Bd += b2;
      Cd += ab;

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

      /* int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)
       * compute the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H
       * for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
       *
       * sh_mu = sqrt(gamma/2) * Mh_mu_nu A^nu, where here gamma = TAtom (as Sinv=1)
       */
      REAL8 norm_s = sqrt(TAtom / 2.0);
      if ( (gslstat = gsl_blas_dgemv (CblasNoTrans, norm_s, Mh_mu_nu, A_Mu, 0.0, sh_mu)) != 0 ) {
        XLALPrintError ( "%s: gsl_blas_dgemv(L * norm) failed: %s\n", fn, gsl_strerror (gslstat) );
        XLAL_ERROR ( fn, XLAL_EFAILED );
      }

      /* add this signal to the atoms, using the relation Fa,Fb <--> x_mu: see Eq.(72) in CFSv2-LIGO-T0900149-v2.pdf */
      REAL8 s1,s2,s3,s4;
      s1 = gsl_vector_get ( sh_mu, 0 );
      s2 = gsl_vector_get ( sh_mu, 1 );
      s3 = gsl_vector_get ( sh_mu, 2 );
      s4 = gsl_vector_get ( sh_mu, 3 );

      atoms->data[alpha].Fa_alpha.re +=   s1;
      atoms->data[alpha].Fa_alpha.im += - s3;
      atoms->data[alpha].Fb_alpha.re +=   s2;
      atoms->data[alpha].Fb_alpha.im += - s4;

    } /* for alpha < numAtoms */

  /* compute optimal SNR^2 expected for this signal */
  REAL8 A1 = gsl_vector_get ( A_Mu, 0 );
  REAL8 A2 = gsl_vector_get ( A_Mu, 1 );
  REAL8 A3 = gsl_vector_get ( A_Mu, 2 );
  REAL8 A4 = gsl_vector_get ( A_Mu, 3 );

  (*rho2) = TAtom  * ( Ad * ( SQ(A1) + SQ(A3) ) + 2.0*Cd * ( A1*A2 + A3*A4 ) + Bd * ( SQ(A2) + SQ(A4) ) );

  /* free memory */
  gsl_vector_free ( sh_mu );
  gsl_matrix_free ( Mh_mu_nu );

  /* return status */
  return XLAL_SUCCESS;

} /* XLALAddSignalToFstatAtomVector() */


/** Add given signal s_mu = M_mu_nu A^nu within the given transient-window
 * to multi-IFO noise-atoms
 */
int
XLALAddSignalToMultiFstatAtomVector ( MultiFstatAtomVector* multiAtoms,	 /**< [in/out] multi atoms vectors containing antenna-functions and possibly noise {Fa,Fb} */
                                      REAL8 *rho2,			 /**< [out] expected optimal SNR^2 for this signal */
                                      const gsl_vector *A_Mu,	 	/**< [in] input canonical amplitude vector A^mu = {A1,A2,A3,A4} */
                                      transientWindow_t transientWindow /**< transient signal window */
                                      )
{
  const char *fn = __func__;

  /* check input consistency */
  if ( !multiAtoms || !multiAtoms->data ) {
    XLALPrintError ( "%s: Invalid NULL input 'multiAtoms'\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( !rho2 ) {
    XLALPrintError ( "%s: Invalid NULL input 'rho2'\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( !A_Mu || (A_Mu->size != 4) ) {
    XLALPrintError ( "%s: Invalid input vector A_Mu: must be allocated 4D\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  UINT4 numDet = multiAtoms->length;
  UINT4 X;

  (*rho2) = 0;
  for ( X=0; X < numDet; X ++ )
    {
      REAL8 rho2X;

      if ( XLALAddSignalToFstatAtomVector ( multiAtoms->data[X], &rho2X, A_Mu, transientWindow ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: XLALAddSignalToFstatAtomVector() failed.\n", fn );
        XLAL_ERROR ( fn, XLAL_EFUNC );
      }

      (*rho2) += rho2X;		/* multi-IFO SNR^2 = sum_X SNR_X^2 */

    } /* for X < numDet */


  return XLAL_SUCCESS;

} /* XLALAddSignalToMultiFstatAtomVector() */

/** Load Ephemeris from ephemeris data-files  */
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
XLALSynthesizeTransientAtoms ( ConfigVariables *cfg,		/**< [in] input params for transient atoms synthesis */
                               BOOLEAN SignalOnly,		/**< signal-only switch: add noise or not */
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

  /* if Alpha < 0 ==> draw skyposition isotropically from all-sky */
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

  /* generate antenna-pattern functions for this sky-position */
  const MultiNoiseWeights *weights = NULL;	/* NULL = unit weights */
  if ( !multiAMBuffer->multiAM && (multiAMBuffer->multiAM = XLALComputeMultiAMCoeffs ( cfg->multiDetStates, weights, skypos )) == NULL ) {
    XLALPrintError ( "%s: XLALComputeMultiAMCoeffs() failed with xlalErrno = %d\n", fn, xlalErrno );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }
  multiAMBuffer->skypos = skypos; /* store buffered skyposition */

  /* generate a pre-initialized F-stat atom vector containing only the antenna-pattern coefficients */
  MultiFstatAtomVector *multiAtoms;
  if ( (multiAtoms = XLALGenerateMultiFstatAtomVector ( cfg->multiTS, multiAMBuffer->multiAM )) == NULL ) {
    XLALPrintError ( "%s: XLALGenerateMultiFstatAtomVector() failed with xlalErrno = %d\n", fn, xlalErrno );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  /* draw amplitude vector A^mu from given ranges in {h0, cosi, psi, phi0} */
  PulsarAmplitudeParams Amp;
  if ( cfg->AmpPrior.fixedSNR > 0 )	/* special treatment of fixed-SNR: use h0=1, later rescale signal */
    Amp.h0 = 1.0;
  else
    Amp.h0   = XLALDrawFromPDF ( &cfg->AmpPrior.pdf_h0Nat, cfg->rng );

  Amp.cosi = XLALDrawFromPDF ( &cfg->AmpPrior.pdf_cosi, cfg->rng );
  Amp.psi  = XLALDrawFromPDF ( &cfg->AmpPrior.pdf_psi,  cfg->rng );
  Amp.phi0 = XLALDrawFromPDF ( &cfg->AmpPrior.pdf_phi0, cfg->rng );

  gsl_vector *A_Mu;	/* output vector */
  if ( (A_Mu = gsl_vector_alloc ( 4 )) == NULL ) {
    XLALPrintError ( "%s: gsl_vector_alloc (4) failed.\n", fn);
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  /* convert amplitude params to 'canonical' vector coordinates */
  if ( XLALAmplitudeParams2Vect ( A_Mu, &Amp ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALAmplitudeParams2Vect() failed with xlalErrno = %d\n", fn, xlalErrno );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  /* draw transient-window parameters from given ranges using flat priors */
  transientWindow_t injectWindow = empty_transientWindow;
  injectWindow.type = cfg->transientInjectRange.type;
  if ( injectWindow.type != TRANSIENT_NONE )	/* nothing to be done if no window */
    {
      injectWindow.t0  = (UINT4) gsl_ran_flat ( cfg->rng, cfg->transientInjectRange.t0, cfg->transientInjectRange.t0 + cfg->transientInjectRange.t0Band );
      injectWindow.tau = (UINT4) gsl_ran_flat ( cfg->rng, cfg->transientInjectRange.tau, cfg->transientInjectRange.tau + cfg->transientInjectRange.tauBand );
    }

  /* add transient signal to the Fstat atoms */
  REAL8 rho2 = 0;
  if ( XLALAddSignalToMultiFstatAtomVector ( multiAtoms, &rho2, A_Mu, injectWindow ) != XLAL_SUCCESS ) {
    XLALPrintError ( "%s: XLALAddSignalToMultiFstatAtomVectorn() failed with xlalErrno = %d\n", fn, xlalErrno );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  if ( cfg->AmpPrior.fixedSNR > 0 )
    {
      REAL8 rescale = cfg->AmpPrior.fixedSNR / sqrt(rho2);	/* rescale atoms by this factor, s.t. SNR = cfg->AmpPrior.SNR */
      if ( XLALRescaleMultiFstatAtomVector ( multiAtoms, rescale ) != XLAL_SUCCESS ) {
        XLALPrintError ( "%s: XLALRescaleMultiFstatAtomVector() failed with xlalErrno = %d\n", fn, xlalErrno );
        XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
      }
    } /* if fixed SNR given */

  /* add noise to the Fstat atoms, unless --SignalOnly was specified */
  if ( !SignalOnly )
    if ( XLALAddNoiseToMultiFstatAtomVector ( multiAtoms, cfg->rng ) != XLAL_SUCCESS ) {
      XLALPrintError ("%s: XLALAddNoiseToMultiFstatAtomVector() failed with xlalErrno = %d\n", fn, xlalErrno );
      XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
    }

  /* we're done, cleanup and return the final atoms-vector */
  gsl_vector_free ( A_Mu );

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


/** Function to setup prior pdfs and sampling-buffers for given ranges (in AmpPrior)
 * and prior types for amplitude-priors ('physical', 'canonical', ... )
 */
int
XLALInitAmplitudePrior ( AmplitudePrior_t *AmpPrior,	/**< [in/out] prior value ranges [xMin,xMax] and pdfs to be initialized */
                         UINT4 Npoints,			/**< [in] how many points to sample pdfs with */
                         AmpPriorType_t priorType	/**< [in] type of amplitude prior */
                         )
{
  const char *fn = __func__;

  /* check input consistency */
  if ( !AmpPrior || !Npoints ){
    XLALPrintError ("%s: invalid NULL input 'AmpPrior' or Npoints==0 \n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( AmpPrior->pdf_h0Nat.prob || AmpPrior->pdf_cosi.prob || AmpPrior->pdf_psi.prob || AmpPrior->pdf_phi0.prob ) {
    XLALPrintError ("%s: non-NULL 'prob' pointer found in AmpPrior struct. All need to be NULL initialized!\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( AmpPrior->pdf_h0Nat.sampling || AmpPrior->pdf_cosi.sampling || AmpPrior->pdf_psi.sampling || AmpPrior->pdf_phi0.sampling ) {
    XLALPrintError ("%s: non-NULL 'sampling' pointer found in AmpPrior struct. All need to be NULL initialized!\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( priorType >= AMP_PRIOR_LAST ) {
    XLALPrintError ("%s: Unknown amplitude-prior type '%d' specified, must be within [0,%d]\n", fn, priorType, AMP_PRIOR_LAST - 1 );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( AmpPrior->pdf_h0Nat.xBand < 0 || AmpPrior->pdf_cosi.xBand < 0 || AmpPrior->pdf_psi.xBand < 0 || AmpPrior->pdf_phi0.xBand < 0 ) {
    XLALPrintError ("%s: only non-negatives values >= 0 for 'xBand' allowed in AmpPrior\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  UINT4 i;
  pdf1D_t *thisPDF;

  switch ( priorType )
    {
    case AMP_PRIOR_PHYSICAL:

      /* flat prior for h0(or SNR) */
      AmpPrior->pdf_h0Nat.prob = NULL;
      /* flat prior for cosi */
      AmpPrior->pdf_cosi.prob = NULL;
      /* flat prior for psi */
      AmpPrior->pdf_psi.prob = NULL;
      /* flat prior for phi0 */
      AmpPrior->pdf_phi0.prob = NULL;

      break;

    case AMP_PRIOR_CANONICAL:

      /* ----- pdf(h0) ~ h0^3 amplitude prior ----- */
      thisPDF = &AmpPrior->pdf_h0Nat;
      if ( thisPDF->xBand > 0 )
        {
          if ( (thisPDF->prob = XLALCreateREAL8Vector ( Npoints )) == NULL ) {
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }
          REAL8 hMin = thisPDF->xMin;
          REAL8 dh   = thisPDF->xBand / Npoints;
          for ( i=0; i < Npoints; i ++ )
            {
              REAL8 hi = hMin + i * dh;
              thisPDF->prob->data[i] = CUBE(hi);
            } /* for i < Npoints */
          thisPDF->normalized = false;

          /* setup precomputed sampling data */
          if ( (thisPDF->sampling = gsl_ran_discrete_preproc ( thisPDF->prob->length, thisPDF->prob->data ) ) == NULL ) {
            XLALPrintError ("%s: gsl_ran_discrete_preproc() failed\n", fn );
            XLAL_ERROR ( fn, XLAL_EFAILED );
          }
        } /* if h0NatBand > 0 */


      /* ----- pdf(cosi) ~ ( 1 - cosi^2)^3 ----- */
      thisPDF = &AmpPrior->pdf_cosi;
      if ( thisPDF->xBand > 0 )
        {
          if ( (thisPDF->prob = XLALCreateREAL8Vector ( Npoints )) == NULL ) {
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }
          REAL8 cosiMin = thisPDF->xMin;
          REAL8 dcosi   = thisPDF->xBand / Npoints;
          for ( i=0; i < Npoints; i ++ )
            {
              REAL8 cosi = cosiMin + i * dcosi;
              REAL8 y = 1.0 - SQ(cosi);
              thisPDF->prob->data[i] = CUBE( y );
            } /* for i < Npoints */
        } /* if cosiBand > 0 */
      thisPDF->normalized = false;

      /* setup precomputed sampling data */
      if ( (thisPDF->sampling = gsl_ran_discrete_preproc ( thisPDF->prob->length, thisPDF->prob->data ) ) == NULL ) {
        XLALPrintError ("%s: gsl_ran_discrete_preproc() failed\n", fn );
        XLAL_ERROR ( fn, XLAL_EFAILED );
      }


      /* ----- flat prior for psi ----- */
      AmpPrior->pdf_psi.prob = NULL;
      /* ----- flat prior for phi0 ----- */
      AmpPrior->pdf_phi0.prob = NULL;

      break;

    default:
      XLALPrintError ("%s: something went wrong ... priorType = %d\n", fn, priorType );
      XLAL_ERROR ( fn, XLAL_EINVAL );
      break;
    } /*   switch ( priorType ) */


  return XLAL_SUCCESS;

} /* XLALInitAmplitudePrior() */


/** Function to generate random samples drawn from the given pdf(x)
 */
REAL8
XLALDrawFromPDF ( const pdf1D_t *probDist,	/**< probability distribution to sample from */
                  const gsl_rng *rng		/**< random-number generator */
                  )
{
  const char *fn = __func__;

  /* input consistency */
  if ( !probDist || !rng ) {
    XLALPrintError ("%s: NULL input 'probDist' or 'rng'\n", fn );
    XLAL_ERROR_REAL8 ( fn, XLAL_EINVAL );
  }
  if ( probDist->xBand < 0 ) {
    XLALPrintError ("%s: xBand '%g' needs to be non-negative >= 0 !\n", fn, probDist->xBand );
    XLAL_ERROR_REAL8 ( fn, XLAL_EINVAL );
  }

  /* ----- special case: single value with certainty */
  if ( probDist->xBand == 0 )
    return probDist->xMin;

  /* ----- special case: uniform pdf in [xMin, xMax] */
  REAL8 xMin = probDist->xMin;
  REAL8 xMax = xMin + probDist->xBand;
  if ( !probDist->prob )
    return gsl_ran_flat ( rng, xMin, xMax );

  /* ----- general case: draw from discretized pdf, buffer gsl_ran_discrete_t if not computed previously */
  REAL8 dx = probDist->xBand / probDist->prob->length;

  if ( !probDist->sampling ) {
    XLALPrintError ("%s: no pre-computed sampling data in pdf.\n", fn );
    XLAL_ERROR_REAL8 ( fn, XLAL_EINVAL );
  }

  UINT4 ind = gsl_ran_discrete (rng, probDist->sampling );
  REAL8 x = xMin + ind * dx;

  return x;

} /* XLALDrawFromPDF() */
