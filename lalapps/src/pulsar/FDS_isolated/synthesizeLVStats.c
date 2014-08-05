/*
 * Copyright (C) 2010 Reinhard Prix
 * Copyright (C) 2011, 2014 David Keitel
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
 * \author R. Prix, D. Keitel
 * \file
 * \brief
 * Generate N samples of various statistics (F-stat, LV-stat,...) drawn from their
 * respective distributions, assuming Gaussian noise, and drawing signal params from
 * their (given) priors
 *
 * This is based on synthesizeBstat and synthesizeTransientStats, and is mostly meant
 * to be used for Monte-Carlos studies of ROC curves
 *
 */

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
#include <lal/LALString.h>
#include <lal/SkyCoordinates.h>
#include <lal/AVFactories.h>
#include <lal/LALInitBarycenter.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lal/ComputeFstat.h>

#include <lal/TransientCW_utils.h>
#include <lal/ProbabilityDensity.h>

#include <lal/SynthesizeCWDraws.h>

#include <lal/LineRobustStats.h>

#include <lal/StringVector.h>

#include <lalapps.h>

/*---------- DEFINES ----------*/
#define SQ(x) ((x)*(x))
#define SQUARE(x) ( (x) * (x) )
#define CUBE(x) ((x)*(x)*(x))
#define QUAD(x) ((x)*(x)*(x)*(x))
#define TRUE (1==1)
#define FALSE (1==0)

/*----- Macros ----- */

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

  /* other parameters */
  LALStringVector* IFOs; /**< list of detector-names "H1,H2,L1,.." or single detector*/
  CHAR *lineIFO;         /**< name of IFO into which line gets inserted */
  INT4 dataStartGPS;	/**< data start-time in GPS seconds */
  INT4 dataDuration;	/**< data-span to generate */
  INT4 TAtom;		/**< Fstat atoms time baseline */

  BOOLEAN computeLV;	/**< Also compute LineVeto-statistic */
  REAL8 LVrho;		/**< prior rho_max_line for LineVeto-statistic */
  LALStringVector* LVlX; /**< Line-to-gauss prior ratios lX for LineVeto statistic */

  LALStringVector* sqrtSX; /**< per-detector noise PSD sqrt(SX) */

  INT4 numDraws;	/**< number of random 'draws' to simulate for F-stat and B-stat */

  CHAR *outputStats;	/**< output file to write numDraw resulting statistics into */
  CHAR *outputAtoms;	/**< output F-statistic atoms into a file with this basename */
  CHAR *outputInjParams;/**< output injection parameters into this file */
  BOOLEAN outputMmunuX;	/**< Whether to write the per-IFO antenna pattern matrices into the parameter file */
  BOOLEAN SignalOnly;	/**< dont generate noise-draws: will result in non-random 'signal only' values of F and B */

  BOOLEAN useFReg;	/**< use 'regularized' Fstat (1/D)*e^F for marginalization, or 'standard' e^F */

  CHAR *ephemEarth;	/**< Earth ephemeris file to use */
  CHAR *ephemSun;	/**< Sun ephemeris file to use */

  BOOLEAN version;	/**< output version-info */
  INT4 randSeed;	/**< GSL random-number generator seed value to use */
} UserInput_t;

/**
 * Configuration settings required for and defining a coherent pulsar search.
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

  REAL8 LVlogRhoTerm;		/**< For LineVeto statistic: extra term coming from prior normalization: log(rho_max_line^4/70) */
  REAL8Vector *LVloglX;		/**< For LineVeto statistic: vector of logs of line prior ratios lX per detector */

  MultiNoiseWeights *multiNoiseWeights;	/**< per-detector noise weights SX^-1/S^-1, no per-SFT variation (can be NULL for unit weights) */

} ConfigVariables;

/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);

int XLALInitUserVars ( UserInput_t *uvar );
int XLALInitCode ( ConfigVariables *cfg, const UserInput_t *uvar );
int XLALInitAmplitudePrior ( AmplitudePrior_t *AmpPrior, const UserInput_t *uvar );
MultiLIGOTimeGPSVector * XLALCreateMultiLIGOTimeGPSVector ( UINT4 numDetectors );
int write_LV_candidate_to_fp ( FILE *fp, const LVcomponents *LVstat, const LALStringVector *IFOs, const InjParams_t *injParams );
MultiNoiseWeights * XLALComputeConstantMultiNoiseWeightsFromNoiseFloor (const MultiNoiseFloor *multiNoiseFloor, const MultiLIGOTimeGPSVector *multiTS, const UINT4 Tsft );

/* exportable API */

/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

/*----------------------------------------------------------------------*/
/* Main Function starts here */
/*----------------------------------------------------------------------*/
/**
 * MAIN function
 * Generates samples of B-stat and F-stat according to their pdfs for given signal-params.
 */
int main(int argc,char *argv[])
{
  UserInput_t XLAL_INIT_DECL(uvar);
  ConfigVariables XLAL_INIT_DECL(cfg);

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

  REAL8 *loglX = NULL;
  if ( cfg.LVloglX ) {
    loglX = cfg.LVloglX->data;
  }

  /* compare IFO name for line injection with IFO list, find the corresponding index, or throw an error if not found */
  UINT4 numDetectors = cfg.multiDetStates->length;
  INT4 lineX = -1;
  if ( uvar.lineIFO ) {
    for ( UINT4 X=0; X < numDetectors; X++ ) {
      if ( strcmp( uvar.lineIFO, uvar.IFOs->data[X] ) == 0 )
        lineX = X;
    }
    if ( lineX == -1 ) {
      XLALPrintError ("\nError in function %s, line %d : Could not match detector ID \"%s\" for line injection to any detector.\n\n", __func__, __LINE__, uvar.lineIFO);
      XLAL_ERROR ( XLAL_EFAILED );
    }
  }

  /* ----- prepare stats output ----- */
  FILE *fpStats = NULL;
  if ( uvar.outputStats )
    {
      if ( (fpStats = fopen (uvar.outputStats, "wb")) == NULL)
	{
	  LogPrintf (LOG_CRITICAL, "Error opening file '%s' for writing..\n\n", uvar.outputStats );
	  XLAL_ERROR ( XLAL_EIO );
	}
      fprintf (fpStats, "%s", cfg.logString );		/* write search log comment */
      if ( write_LV_candidate_to_fp ( fpStats, NULL, uvar.IFOs, NULL ) != XLAL_SUCCESS ) { /* write header-line comment */
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
      if ( write_InjParams_to_fp ( fpInjParams, NULL, 0, uvar.outputMmunuX, numDetectors ) != XLAL_SUCCESS ) { /* write header-line comment */
        XLAL_ERROR ( XLAL_EFUNC );
      }
    } /* if outputInjParams */

  multiAMBuffer_t XLAL_INIT_DECL(multiAMBuffer);      /* prepare AM-buffer */

  /* ----- main MC loop over numDraws trials ---------- */
  INT4 i;
  for ( i=0; i < uvar.numDraws; i ++ )
    {
      InjParams_t XLAL_INIT_DECL(injParamsDrawn);

      /* ----- generate signal random draws from ranges and generate Fstat atoms */
      MultiFstatAtomVector *multiAtoms;

      multiAtoms = XLALSynthesizeTransientAtoms ( &injParamsDrawn, cfg.skypos, cfg.AmpPrior, cfg.transientInjectRange, cfg.multiDetStates, cfg.SignalOnly, &multiAMBuffer, cfg.rng, lineX, cfg.multiNoiseWeights );
      XLAL_CHECK ( multiAtoms != NULL, XLAL_EFUNC );

      /* ----- if requested, output signal injection parameters into file */
      if ( fpInjParams && (write_InjParams_to_fp ( fpInjParams, &injParamsDrawn, uvar.dataStartGPS, uvar.outputMmunuX, numDetectors ) != XLAL_SUCCESS ) ) {
        XLAL_ERROR ( XLAL_EFUNC );
      } /* if fpInjParams & failure*/

      /* initialise LVcomponents structure and allocate memory */
      LVcomponents   lvstats;      /* struct containing multi-detector Fstat, single-detector Fstats, Line Veto stat */
      if ( (lvstats.TwoFX = XLALCreateREAL4Vector ( numDetectors )) == NULL ) {
        XLALPrintError ("%s: failed to XLALCreateREAL4Vector( %d )\n", __func__, numDetectors );
        XLAL_ERROR ( XLAL_EFUNC );
      }

      /* compute F and LV statistics from atoms */
      UINT4 X;
      for ( X=0; X < numDetectors; X++ )    {
        lvstats.TwoFX->data[X] = XLALComputeFstatFromAtoms ( multiAtoms, X );
        if ( xlalErrno != 0 ) {
          XLALPrintError ("\nError in function %s, line %d : Failed call to XLALComputeFstatFromAtoms().\n\n", __func__, __LINE__);
          XLAL_ERROR ( XLAL_EFUNC );
        }
      }

      lvstats.TwoF = XLALComputeFstatFromAtoms ( multiAtoms, -1 );
      if ( xlalErrno != 0 ) {
        XLALPrintError ("\nError in function %s, line %d : Failed call to XLALComputeFstatFromAtoms().\n\n", __func__, __LINE__);
        XLAL_ERROR ( XLAL_EFUNC );
      }

      if ( uvar.computeLV ) {
        BOOLEAN useAllTerms = TRUE;
        lvstats.LV = XLALComputeLineVetoArray ( (REAL4)lvstats.TwoF, numDetectors, (REAL4*)lvstats.TwoFX->data, cfg.LVlogRhoTerm, loglX, useAllTerms );
        if ( xlalErrno != 0 ) {
          XLALPrintError ("\nError in function %s, line %d : Failed call to XLALComputeLineVeto().\n\n", __func__, __LINE__);
          XLAL_ERROR ( XLAL_EFUNC );
        }
      }
      else {
       lvstats.LV = 0.0;
      }

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


      /* ----- if requested, output transient-cand statistics */
      if ( fpStats && write_LV_candidate_to_fp ( fpStats, &lvstats, uvar.IFOs, &injParamsDrawn ) != XLAL_SUCCESS ) {
        XLALPrintError ( "%s: write_transientCandidate_to_fp() failed.\n", __func__ );
        XLAL_ERROR ( XLAL_EFUNC );
      }

      /* ----- free Memory */
      XLALDestroyREAL4Vector ( lvstats.TwoFX );
      XLALDestroyMultiFstatAtomVector ( multiAtoms );

    } /* for i < numDraws */

  /* ----- close files ----- */
  if ( fpStats ) fclose ( fpStats );
  if ( fpInjParams ) fclose ( fpInjParams );

  /* ----- free memory ---------- */
  XLALDestroyMultiDetectorStateSeries ( cfg.multiDetStates );
  XLALDestroyMultiNoiseWeights ( cfg.multiNoiseWeights );
  XLALDestroyExpLUT();
  XLALDestroyMultiAMCoeffs ( multiAMBuffer.multiAM );
  /* ----- free amplitude prior pdfs ----- */
  XLALDestroyPDF1D ( cfg.AmpPrior.pdf_h0Nat );
  XLALDestroyPDF1D ( cfg.AmpPrior.pdf_cosi );
  XLALDestroyPDF1D ( cfg.AmpPrior.pdf_psi );
  XLALDestroyPDF1D ( cfg.AmpPrior.pdf_phi0 );

  if ( cfg.logString ) XLALFree ( cfg.logString );
  gsl_rng_free ( cfg.rng );

  XLALDestroyREAL8Vector ( cfg.LVloglX );

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
  uvar->dataDuration = (INT4) round ( LAL_YRSID_SI ) ;	/* 1 year of data */

  uvar->ephemEarth = XLALStringDuplicate("earth00-19-DE405.dat.gz");
  uvar->ephemSun = XLALStringDuplicate("sun00-19-DE405.dat.gz");

  uvar->numDraws = 1;
  uvar->TAtom = 1800;

  uvar->computeLV = 0;
  uvar->LVrho = 0.0;
  uvar->LVlX = NULL;
  uvar->sqrtSX = NULL;
  uvar->useFReg = 0;

  uvar->fixedh0Nat = -1;
  uvar->fixedSNR = -1;
  uvar->fixedh0NatMax = -1;
  uvar->fixedRhohMax = -1;

  if ( (uvar->IFOs = XLALCreateStringVector ( "H1", NULL )) == NULL ) {
    LogPrintf (LOG_CRITICAL, "Call to XLALCreateStringVector() failed with xlalErrno = %d\n", xlalErrno );
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  uvar->lineIFO = NULL;

  /* ---------- transient window defaults ---------- */
#define DEFAULT_TRANSIENT "rect"

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

  XLALregLISTUserStruct( IFOs,                  'I', UVAR_OPTIONAL, "Comma-separated list of detectors, eg. \"H1,H2,L1,G1, ...\" ");
  XLALregSTRINGUserStruct ( lineIFO,             0,  UVAR_OPTIONAL, "Insert a line (signal in this one IFO, pure gaussian noise in others), e.g. \"H1\"");
  XLALregINTUserStruct ( dataStartGPS,	 	 0,  UVAR_OPTIONAL, "data start-time in GPS seconds");
  XLALregINTUserStruct ( dataDuration,	 	 0,  UVAR_OPTIONAL, "data-span to generate (in seconds)");

  /* misc params */
  XLALregBOOLUserStruct ( computeLV,		 0, UVAR_OPTIONAL, "Also compute LineVeto-statistic");
  XLALregREALUserStruct ( LVrho,		 0, UVAR_OPTIONAL, "LineVeto: prior rho_max_line for LineVeto-statistic");
  XLALregLISTUserStruct ( LVlX,			 0, UVAR_OPTIONAL, "LineVeto: line-to-gauss prior ratios lX for different detectors X, length must be numDetectors. Defaults to lX=1,1,..");

  XLALregLISTUserStruct ( sqrtSX,		 0, UVAR_OPTIONAL, "Per-detector noise PSD sqrt(SX). Only ratios relevant to compute noise weights. Defaults to 1,1,...");

  XLALregINTUserStruct  ( numDraws,		'N', UVAR_OPTIONAL,"Number of random 'draws' to simulate");
  XLALregINTUserStruct  ( randSeed,		 0, UVAR_OPTIONAL, "GSL random-number generator seed value to use");

  XLALregSTRINGUserStruct ( outputStats,	'o', UVAR_OPTIONAL, "Output file containing 'numDraws' random draws of stats");
  XLALregSTRINGUserStruct ( outputAtoms,	 0,  UVAR_OPTIONAL, "Output F-statistic atoms into a file with this basename");
  XLALregSTRINGUserStruct ( outputInjParams,	 0,  UVAR_OPTIONAL, "Output injection parameters into this file");
  XLALregBOOLUserStruct ( outputMmunuX,        	 0,  UVAR_OPTIONAL, "Write the per-IFO antenna pattern matrices into the parameter file");

  XLALregBOOLUserStruct ( SignalOnly,        	'S', UVAR_OPTIONAL, "Signal only: generate pure signal without noise");
  XLALregBOOLUserStruct ( useFReg,        	 0,  UVAR_OPTIONAL, "use 'regularized' Fstat (1/D)*e^F (if TRUE) for marginalization, or 'standard' e^F (if FALSE)");

  XLALregSTRINGUserStruct ( ephemEarth, 	 0,  UVAR_OPTIONAL, "Earth ephemeris file to use");
  XLALregSTRINGUserStruct ( ephemSun, 	 	 0,  UVAR_OPTIONAL, "Sun ephemeris file to use");

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
  EphemerisData *edat = XLALInitBarycenter( uvar->ephemEarth, uvar->ephemSun );
  if ( !edat ) {
    LogPrintf ( LOG_CRITICAL, "%s: XLALInitBarycenter failed: could not load Earth ephemeris '%s' and Sun ephemeris '%s'\n", __func__, uvar->ephemEarth, uvar->ephemSun);
    XLAL_ERROR ( XLAL_EFUNC );
  }

  UINT4 numDetectors = uvar->IFOs->length;
  MultiLALDetector multiDet;
  XLAL_CHECK ( XLALParseMultiLALDetector ( &multiDet, uvar->IFOs ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* init timestamps vector covering observation time */
  UINT4 numSteps = (UINT4) ceil ( uvar->dataDuration / uvar->TAtom );
  MultiLIGOTimeGPSVector * multiTS;
  if ( (multiTS = XLALCreateMultiLIGOTimeGPSVector (numDetectors)) == NULL ) {
     XLALPrintError ("%s: XLALCreateMultiLIGOTimeGPSVector(%d) failed.\n", __func__, numDetectors );
  }

  for ( UINT4 X=0; X < numDetectors; X++ )    {
    if ( (multiTS->data[X] = XLALCreateTimestampVector (numSteps)) == NULL ) {
      XLALPrintError ("%s: XLALCreateTimestampVector(%d) failed.\n", __func__, numSteps );
    }
    multiTS->data[X]->deltaT = uvar->TAtom;
    UINT4 i;
    for ( i=0; i < numSteps; i ++ )
      {
	UINT4 ti = uvar->dataStartGPS + i * uvar->TAtom;
	multiTS->data[X]->data[i].gpsSeconds = ti;
	multiTS->data[X]->data[i].gpsNanoSeconds = 0;
      }
  }

  /* get detector states */
  if ( (cfg->multiDetStates = XLALGetMultiDetectorStates ( multiTS, &multiDet, edat, 0.5 * uvar->TAtom )) == NULL ) {
    XLALPrintError ( "%s: XLALGetMultiDetectorStates() failed.\n", __func__ );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  /* process LV user vars */
  if ( uvar->LVrho < 0.0 ) {
    fprintf(stderr, "Invalid LV prior rho (given rho=%f, need rho>=0)!\n", uvar->LVrho);
    return( XLAL_EINVAL );
  }
  else if ( uvar->LVrho > 0.0 ) {
    cfg->LVlogRhoTerm = 4.0 * log(uvar->LVrho) - log(70.0);
  }
  else { /* if uvar_LVrho == 0.0, logRhoTerm should become irrelevant in summation */
    cfg->LVlogRhoTerm = - LAL_REAL4_MAX;
  }

  /* set up line prior ratios: either given by user, then convert from string to REAL4 vector; else, pass NULL, which is interpreted as lX=1.0 for all X */
  cfg->LVloglX = NULL;

  if ( uvar->computeLV && uvar->LVlX ) {

    if (  uvar->LVlX->length != numDetectors ) {
      fprintf(stderr, "Length of LV prior ratio vector does not match number of detectors! (%d != %d)\n", uvar->LVlX->length, numDetectors);
      XLAL_ERROR ( XLAL_EINVAL );
    }

    if ( (cfg->LVloglX = XLALCreateREAL8Vector ( numDetectors )) == NULL ) {
      fprintf(stderr, "Failed call to XLALCreateREAL8Vector( %d )\n", numDetectors );
      XLAL_ERROR ( XLAL_EFUNC );
    }

    for (UINT4 X = 0; X < numDetectors; X++) {

      if ( 1 != sscanf ( uvar->LVlX->data[X], "%" LAL_REAL8_FORMAT, &cfg->LVloglX->data[X] ) ) {
        fprintf(stderr, "Illegal REAL8 commandline argument to --LVlX[%d]: '%s'\n", X, uvar->LVlX->data[X]);
        XLAL_ERROR ( XLAL_EINVAL );
      }

      if ( cfg->LVloglX->data[X] < 0.0 ) {
        fprintf(stderr, "Negative input prior-ratio for detector X=%d lX[X]=%f\n", X, cfg->LVloglX->data[X] );
        XLAL_ERROR ( XLAL_EINVAL );
      }
      else if ( cfg->LVloglX->data[X] > 0.0 ) {
        cfg->LVloglX->data[X] = log(cfg->LVloglX->data[X]);
      }
      else { /* if zero prior ratio, approximate log(0)=-inf by -LAL_REA4_MAX to avoid raising underflow exceptions */
        cfg->LVloglX->data[X] = - LAL_REAL8_MAX;
      }

    } /* for X < numDetectors */

  } /* if ( uvar->computeLV && uvar->LVlX ) */

  if ( uvar->sqrtSX ) { /* translate user-input PSD sqrt(SX) to noise-weights (this actually does not care whether they were normalized or not) */

    /* parse input comma-separated list */
    MultiNoiseFloor multiNoiseFloor;
    XLAL_CHECK ( XLALParseMultiNoiseFloor ( &multiNoiseFloor, uvar->sqrtSX, numDetectors ) == XLAL_SUCCESS, XLAL_EFUNC );

    /* translate to noise weights */
    XLAL_CHECK ( ( cfg->multiNoiseWeights = XLALComputeConstantMultiNoiseWeightsFromNoiseFloor ( &multiNoiseFloor, multiTS, uvar->TAtom ) ) != NULL, XLAL_EFUNC );

  } /* if ( uvar->sqrtSX ) */

  /* get rid of all temporary memory allocated for this step */
  XLALDestroyEphemerisData ( edat );
  XLALDestroyMultiTimestamps ( multiTS );
  multiTS = NULL;

  /* ---------- initialize transient window ranges, for injection ... ---------- */
  cfg->transientInjectRange.type = TRANSIENT_NONE;			/* default: no transient signal window */
  /* apply correct defaults if unset: t0=dataStart, t0Band=dataDuration-3*tauMax */
//   cfg->transientInjectRange.t0 = uvar->dataStartGPS + uvar->injectWindow_t0Days * DAY24;

  cfg->transientSearchRange = cfg->transientInjectRange;
  return XLAL_SUCCESS;

} /* XLALInitCode() */


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


/**
 * Simple creator function for MultiLIGOTimeGPSVector with numDetectors entries
 */
MultiLIGOTimeGPSVector *
XLALCreateMultiLIGOTimeGPSVector ( UINT4 numDetectors )
{
  MultiLIGOTimeGPSVector *ret;

  if ( (ret = XLALMalloc ( sizeof(*ret) )) == NULL ) {
    XLALPrintError ("%s: XLALMalloc(%d) failed.\n", __func__, sizeof(*ret) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  ret->length = numDetectors;
  if ( (ret->data = XLALCalloc ( numDetectors, sizeof(*ret->data) )) == NULL ) {
    XLALPrintError ("%s: XLALCalloc(%d, %d) failed.\n", __func__, numDetectors, sizeof(*ret->data) );
    XLALFree ( ret );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  return ret;

} /* XLALCreateMultiLIGOTimeGPSVector() */


/**
 * Write one line for given LV candidate into output file.
 *
 * NOTE: input LVstat and injParams can be NULL pointers, then writes header comment instead
 *
 */
int
write_LV_candidate_to_fp ( FILE *fp, const LVcomponents *LVstat, const LALStringVector *IFOs, const InjParams_t *injParams )
{

  XLAL_CHECK ( fp, XLAL_EINVAL, "Invalid NULL filepointer input." );

  /* if requested, write header-comment line */
  if ( !LVstat && !injParams ) {

    char stat_header_string[256] = "";
    char buf0[256];
    snprintf ( stat_header_string, sizeof(stat_header_string), "2F      " );
    for ( UINT4 X = 0; X < IFOs->length ; X ++ ) {
      snprintf ( buf0, sizeof(buf0), " 2F_%s   ", IFOs->data[X] );
      UINT4 len1 = strlen ( stat_header_string ) + strlen ( buf0 ) + 1;
      XLAL_CHECK ( len1 <= sizeof ( stat_header_string ), XLAL_EBADLEN, "Assembled output string too long! (%d > %d)", len1, sizeof(stat_header_string));
      strcat ( stat_header_string, buf0 );
    }
    snprintf ( buf0, sizeof(buf0), " LV" );
    strcat ( stat_header_string, buf0 );
    fprintf(fp, "%%%% freq  alpha    delta    f1dot    %s\n", stat_header_string);

    return XLAL_SUCCESS;	/* we're done here */

  } /* if par == NULL */

  /* sanity checks */
  XLAL_CHECK ( LVstat && LVstat->TwoFX && LVstat->TwoFX->data, XLAL_EFAULT, "Invalid LVstat pointer as input parameter!" );
  XLAL_CHECK ( injParams, XLAL_EFAULT, "Invalid injParams pointer as input parameter!" );

  /* add output-field containing twoF and per-detector sumTwoFX */
  char statString[256] = "";	/* defaults to empty */
  char buf0[256];
  snprintf ( statString, sizeof(statString), "%.6f", LVstat->TwoF );
  for ( UINT4 X = 0; X < IFOs->length ; X ++ ) {
    snprintf ( buf0, sizeof(buf0), " %.6f", LVstat->TwoFX->data[X] );
    UINT4 len1 = strlen ( statString ) + strlen ( buf0 ) + 1;
    XLAL_CHECK ( len1 <= sizeof ( statString ), XLAL_EBADLEN, "Assembled output string too long! (%d > %d)", len1, sizeof(statString));
    strcat ( statString, buf0 );
  } /* for X < IFOs->length */
  snprintf ( buf0, sizeof(buf0), " %.6f", LVstat->LV );
  strcat ( statString, buf0 );
  fprintf (fp, "%.6f %.6f %.6f %.6f %s\n",
		0.0, /* legacy field from synthesizeTransientStats: freq, not needed here */
		injParams->skypos.longitude,
		injParams->skypos.latitude,
		0.0,
		statString /* legacy field from synthesizeTransientStats: f1dot, not needed here */
	    );

  return XLAL_SUCCESS;

} /* write_LV_candidate_to_fp() */


/**
 *
 */
MultiNoiseWeights *
XLALComputeConstantMultiNoiseWeightsFromNoiseFloor ( const MultiNoiseFloor *multiNoiseFloor,	/**< [in] noise floor values sqrt(S) for all detectors */
                                                     const MultiLIGOTimeGPSVector *multiTS,	/**< [in] timestamps vectors for all detectors, only needed for their lengths */
                                                     const UINT4 Tsft				/**< [in] length of SFTs in secons, needed for normalization factor Sinv_Tsft */
                                                     )
{

  /* check input parameters */
  XLAL_CHECK_NULL ( multiNoiseFloor != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( multiTS != NULL, XLAL_EINVAL );
  UINT4 numDet = multiNoiseFloor->length;
  XLAL_CHECK_NULL ( numDet == multiTS->length, XLAL_EINVAL, "Inconsistent length between multiNoiseFloor (%d) and multiTimeStamps (%d) structs.\n", numDet, multiTS->length );

  /* create multi noise weights for output */
  MultiNoiseWeights *multiWeights = NULL;
  XLAL_CHECK_NULL ( (multiWeights = XLALCalloc(1, sizeof(*multiWeights))) != NULL, XLAL_ENOMEM, "Failed call to XLALCalloc ( 1, %d )\n", sizeof(*multiWeights) );
  XLAL_CHECK_NULL ( (multiWeights->data = XLALCalloc(numDet, sizeof(*multiWeights->data))) != NULL, XLAL_ENOMEM, "Failed call to XLALCalloc ( %d, %d )\n", numDet, sizeof(*multiWeights->data) );
  multiWeights->length = numDet;

  REAL8 sqrtSnTotal = 0;
  UINT4 numSFTsTotal = 0;
  for (UINT4 X = 0; X < numDet; X++) { /* first loop over detectors: compute total noise floor normalization */
    sqrtSnTotal  += multiTS->data[X]->length / SQ ( multiNoiseFloor->sqrtSn[X] ); /* actually summing up 1/Sn, not sqrtSn yet */
    numSFTsTotal += multiTS->data[X]->length;
  }
  sqrtSnTotal = sqrt ( numSFTsTotal / sqrtSnTotal ); /* SnTotal = harmonicMean{S_X} assuming per-detector stationarity */

  for (UINT4 X = 0; X < numDet; X++) { /* second loop over detectors: compute the weights */

    /* compute per-IFO weights */
    REAL8 noise_weight_X = SQ(sqrtSnTotal/multiNoiseFloor->sqrtSn[X]); /* w_Xalpha = S_Xalpha^-1/S^-1 = S / S_Xalpha */

    /* create k^th weights vector */
    if( ( multiWeights->data[X] = XLALCreateREAL8Vector ( multiTS->data[X]->length ) ) == NULL ) {
      /* free weights vectors created previously in loop */
      XLALDestroyMultiNoiseWeights ( multiWeights );
      XLAL_ERROR_NULL ( XLAL_EFUNC, "Failed to allocate noiseweights for IFO X = %d\n", X );
    } /* if XLALCreateREAL8Vector() failed */

    /* loop over SFT timestamps and use same weights for all */
    for ( UINT4 alpha = 0; alpha < multiTS->data[X]->length; alpha++) {
      multiWeights->data[X]->data[alpha] = noise_weight_X;
    }

  } /* for X < numDet */

  multiWeights->Sinv_Tsft = Tsft / SQ ( sqrtSnTotal );

  return multiWeights;

} /* XLALComputeConstantMultiNoiseWeightsFromNoiseFloor() */
