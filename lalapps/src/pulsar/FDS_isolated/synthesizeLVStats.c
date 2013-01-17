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
/** \author R. Prix, D. Keitel
 * \file
 * \brief
 * Generate N samples of various statistics (F-stat, LV-stat,...) drawn from their
 * respective distributions, assuming Gaussian noise, and drawing signal params from
 * their (given) priors
 *
 * This is based on synthesizeBstat and synthesizeTransientStats, and is mostly meant
 * to be used for Monte-Carlos studies of ROC curves
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

#include <../GCT/LineVeto.h>

#include <lal/StringVector.h>

#include <lalapps.h>

/*---------- DEFINES ----------*/
#define EPHEM_YEARS  "00-19-DE405"	/**< default range, covering S5: override with --ephemYear */
#define SQ(x) ((x)*(x))
#define SQUARE(x) ( (x) * (x) )
#define CUBE(x) ((x)*(x)*(x))
#define QUAD(x) ((x)*(x)*(x)*(x))
#define TRUE (1==1)
#define FALSE (1==0)

/*----- Macros ----- */
#define INIT_MEM(x) memset(&(x), 0, sizeof((x)))


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
  INT4 numDraws;	/**< number of random 'draws' to simulate for F-stat and B-stat */

  CHAR *outputStats;	/**< output file to write numDraw resulting statistics into */
  CHAR *outputAtoms;	/**< output F-statistic atoms into a file with this basename */
  CHAR *outputInjParams;/**< output injection parameters into this file */
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

/** multi-detector array of InjParams_t types */
typedef struct {
   UINT4 length;            /**< number of detectors */
   InjParams_t *data;       /**< array of InjParams_t (pointers), one for each detector X */
 } MultiInjParams;

/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);

int XLALInitUserVars ( UserInput_t *uvar );
int XLALInitCode ( ConfigVariables *cfg, const UserInput_t *uvar );
EphemerisData * XLALInitEphemeris (const CHAR *ephemYear );
int XLALInitAmplitudePrior ( AmplitudePrior_t *AmpPrior, const UserInput_t *uvar );
MultiLIGOTimeGPSVector * XLALCreateMultiLIGOTimeGPSVector ( UINT4 numDetectors );
int write_LV_candidate_to_fp ( FILE *fp, const LVcomponents *LVstat, const PulsarDopplerParams *dopplerParams_in );
MultiInjParams * XLALCreateMultiInjParams ( UINT4 numDetectors );
void XLALDestroyMultiInjParams ( MultiInjParams *multipar );
InjParams_t * XLALCombineInjParamsForLine( const MultiInjParams *injParamsDrawnX, const UINT4 lineX );

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

  lalDebugLevel = 0;
  vrbflg = 1;	/* verbose error-messages */
  LogSetLevel(lalDebugLevel);

  /* turn off default GSL error handler */
  gsl_set_error_handler_off ();

  /* ----- register and read all user-variables ----- */
  if ( XLALGetDebugLevel ( argc, argv, 'v') != XLAL_SUCCESS ) {
    LogPrintf ( LOG_CRITICAL, "%s: XLALGetDebugLevel() failed with errno=%d\n", __func__, xlalErrno );
    return 1;
  }
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
  FILE *fpStats = NULL;
  if ( uvar.outputStats )
    {
      if ( (fpStats = fopen (uvar.outputStats, "wb")) == NULL)
	{
	  LogPrintf (LOG_CRITICAL, "Error opening file '%s' for writing..\n\n", uvar.outputStats );
	  XLAL_ERROR ( XLAL_EIO );
	}
      fprintf (fpStats, "%s", cfg.logString );		/* write search log comment */
//       if ( write_LV_candidate_to_fp ( fpStats, NULL, NULL ) != XLAL_SUCCESS ) { /* write header-line comment */
//         XLAL_ERROR ( XLAL_EFUNC );
//       }
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

  /* ----- main MC loop over numDraws trials ---------- */
  INT4 i;
  for ( i=0; i < uvar.numDraws; i ++ )
    {
      InjParams_t injParamsDrawn = empty_InjParams_t;

      /* ----- generate signal random draws from ranges and generate Fstat atoms */
      multiAMBuffer_t multiAMBuffer = empty_multiAMBuffer;      /* prepare AM-buffer */
      MultiFstatAtomVector *multiAtoms;

      if ( !uvar.lineIFO ) { /* signal injection in all detectors, 1 call to XLALSynthesizeTransientAtoms for multi-detector results */
        multiAtoms = XLALSynthesizeTransientAtoms ( &injParamsDrawn, cfg.skypos, cfg.AmpPrior, cfg.transientInjectRange, cfg.multiDetStates, cfg.SignalOnly, &multiAMBuffer, cfg.rng);
        if ( multiAtoms ==NULL ) {
          LogPrintf ( LOG_CRITICAL, "%s: XLALSynthesizeTransientAtoms() failed with xlalErrno = %d\n", __func__, xlalErrno );
          XLAL_ERROR ( XLAL_EFUNC );
        }
      }

      else { /* inject a signal in detector lineIFO only, pure gaussian noise in the others */

        /* prepare multiAtoms structure (will be filled by hand from single-IFO results of individual calls to XLALSynthesizeTransientAtoms ) */
        if ( ( multiAtoms = XLALCalloc ( 1, sizeof(*multiAtoms) )) == NULL ) {
          XLALPrintError ("%s: XLALCalloc ( 1, %d) failed.\n", __func__, sizeof(*multiAtoms) );
          XLAL_ERROR ( XLAL_ENOMEM );
        }
        multiAtoms->length = numDetectors;
        if ( ( multiAtoms->data = XLALCalloc ( numDetectors, sizeof(*multiAtoms->data) ) ) == NULL ) {
          XLALPrintError ("%s: XLALCalloc ( %d, %d) failed.\n", __func__, numDetectors, sizeof(*multiAtoms->data) );
          XLALFree ( multiAtoms );
          XLAL_ERROR ( XLAL_ENOMEM );
        }

        /* prepare array of injection parameters per detector, so that they can be combined for output afterwards */
        MultiInjParams *injParamsDrawnX = XLALCreateMultiInjParams ( numDetectors );

        for ( UINT4 X=0; X < numDetectors; X++ ) { /* loop through detectors */

          /* finish preparing multiAtoms structure for insertion of atoms for detector X */
          UINT4 numAtoms = cfg.multiDetStates->data[X]->length;
          if ( ( multiAtoms->data[X] = XLALCreateFstatAtomVector ( numAtoms ) ) == NULL ) {
            XLALPrintError ("%s: XLALCreateFstatAtomVector(%d) failed.\n", __func__, numAtoms );
            XLAL_ERROR ( XLAL_EFUNC );
          }

          MultiFstatAtomVector *multiAtomsX; /* temporary multiAtoms structure with only 1 detector entry */
          AmplitudePrior_t AmpPriorX = cfg.AmpPrior; /* for all detectors without line, use temporary AmpPrior struct to set signal strength to 0 */
          if ( X != (UINT4)(lineX) )
            AmpPriorX.fixedSNR = 0.0;
          injParamsDrawnX->data[X] = empty_InjParams_t; /* initialize injection parameter structure for 1 detector */

          /* temporary DetectorStateSeries structure so that XLALSynthesizeTransientAtoms will only synth for detector X */
          MultiDetectorStateSeries multiDetStatesX;
          if ( ( multiDetStatesX.data = LALCalloc ( 1, sizeof( *(multiDetStatesX.data) ) )) == NULL ) {
            XLALPrintError ("%s: LALCalloc ( %d, sizeof(%d)) failed\n", __func__, 1, sizeof(*(multiDetStatesX.data)) );
            XLAL_ERROR ( XLAL_ENOMEM );
          }
          multiDetStatesX.length = 1;
          multiDetStatesX.startTime = cfg.multiDetStates->startTime;
          multiDetStatesX.Tspan = cfg.multiDetStates->Tspan;
          multiDetStatesX.data[0] = cfg.multiDetStates->data[X];

          /* finally, the synth call for this detector X with temporary DetStates and AmpPrior */
          multiAtomsX = XLALSynthesizeTransientAtoms ( &injParamsDrawnX->data[X], cfg.skypos, AmpPriorX, cfg.transientInjectRange, &multiDetStatesX, cfg.SignalOnly, &multiAMBuffer, cfg.rng);
          if ( multiAtomsX == NULL ) {
            LogPrintf ( LOG_CRITICAL, "%s: XLALSynthesizeTransientAtoms() failed with xlalErrno = %d\n", __func__, xlalErrno );
            XLAL_ERROR ( XLAL_EFUNC );
          }

          /* copy single-IFO atoms into multiAtoms struct (manually to avoid LALFree errors - should be done better! */
          multiAtoms->data[X]->length = multiAtomsX->data[0]->length;
          multiAtoms->data[X]->TAtom = multiAtomsX->data[0]->TAtom;
          for ( UINT4 alpha=0; alpha < multiAtomsX->data[0]->length; alpha++ )
            multiAtoms->data[X]->data[alpha] = multiAtomsX->data[0]->data[alpha];

          /* free temporary structs for this detector */
          XLALFree ( multiDetStatesX.data );
          XLALDestroyMultiFstatAtomVector ( multiAtomsX );

        } /* for X < numDetectors */

        /* combine single-IFO injection parameters into one struct for output, free temp struct afterwards */
        injParamsDrawn = *XLALCombineInjParamsForLine( injParamsDrawnX, lineX );
        XLALDestroyMultiInjParams ( injParamsDrawnX );

      } /* finished if/else statement for line injection */

      /* ----- if requested, output signal injection parameters into file */
      if ( fpInjParams && (write_InjParams_to_fp ( fpInjParams, &injParamsDrawn, uvar.dataStartGPS ) != XLAL_SUCCESS ) ) {
        XLAL_ERROR ( XLAL_EFUNC );
      } /* if fpInjParams & failure*/

      /* initialise LVcomponents structure and allocate memory */
      LVcomponents   lvstats;      /* struct containing multi-detector Fstat, single-detector Fstats, Line Veto stat */
      if ( (lvstats.TwoFX = XLALCreateREAL4Vector ( numDetectors )) == NULL ) {
        XLALPrintError ("%s: failed to XLALCreateREAL4Vector( %d )\n", __func__, numDetectors );
        XLAL_ERROR ( XLAL_EFUNC );
      }

      REAL8Vector *linepriorX;
      if ( (linepriorX = XLALCreateREAL8Vector ( numDetectors )) == NULL ) {
        XLALPrintError ("%s: failed to XLALCreateREAL8Vector( %d )\n", __func__, numDetectors );
        XLAL_ERROR ( XLAL_EFUNC );
      }

      /* compute F and LV statistics from atoms */
      UINT4 X;
      for ( X=0; X < numDetectors; X++ )    {
        linepriorX->data[X] = 0.5;
        lvstats.TwoFX->data[X] = 2.0*XLALComputeFstatFromAtoms ( multiAtoms, X );
        if ( xlalErrno != 0 ) {
          XLALPrintError ("\nError in function %s, line %d : Failed call to XLALComputeFstatFromAtoms().\n\n", __func__, __LINE__);
          XLAL_ERROR ( XLAL_EFUNC );
        }
      }

      lvstats.TwoF = 2.0*XLALComputeFstatFromAtoms ( multiAtoms, -1 );
      if ( xlalErrno != 0 ) {
        XLALPrintError ("\nError in function %s, line %d : Failed call to XLALComputeFstatFromAtoms().\n\n", __func__, __LINE__);
        XLAL_ERROR ( XLAL_EFUNC );
      }

      if ( uvar.computeLV ) {
        BOOLEAN useAllTerms = TRUE;
        lvstats.LV = XLALComputeLineVeto ( (REAL4)lvstats.TwoF, (REAL4Vector*)lvstats.TwoFX, uvar.LVrho, linepriorX, useAllTerms );
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
      if ( fpStats && write_LV_candidate_to_fp ( fpStats, &lvstats, NULL ) != XLAL_SUCCESS ) {
        XLALPrintError ( "%s: write_transientCandidate_to_fp() failed.\n", __func__ );
        XLAL_ERROR ( XLAL_EFUNC );
      }

      /* ----- free Memory */
      XLALDestroyREAL4Vector ( lvstats.TwoFX );
      XLALDestroyREAL8Vector ( linepriorX );
      XLALDestroyMultiFstatAtomVector ( multiAtoms );
      XLALDestroyMultiAMCoeffs ( multiAMBuffer.multiAM );

    } /* for i < numDraws */

  /* ----- close files ----- */
  if ( fpStats ) fclose ( fpStats );
  if ( fpInjParams ) fclose ( fpInjParams );

  /* ----- free memory ---------- */
  XLALDestroyMultiDetectorStateSeries ( cfg.multiDetStates );
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
  uvar->dataDuration = (INT4) round ( LAL_YRSID_SI ) ;	/* 1 year of data */

  uvar->ephemYear = XLALCalloc (1, strlen(EPHEM_YEARS)+1);
  strcpy (uvar->ephemYear, EPHEM_YEARS);

  uvar->numDraws = 1;
  uvar->TAtom = 1800;

  uvar->computeLV = 0;
  uvar->LVrho = 0.0;
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
  XLALregREALUserStruct ( LVrho,		 0, UVAR_OPTIONAL, "prior rho_max_line for LineVeto-statistic");

  XLALregINTUserStruct  ( numDraws,		'N', UVAR_OPTIONAL,"Number of random 'draws' to simulate");
  XLALregINTUserStruct  ( randSeed,		 0, UVAR_OPTIONAL, "GSL random-number generator seed value to use");

  XLALregSTRINGUserStruct ( outputStats,	'o', UVAR_OPTIONAL, "Output file containing 'numDraws' random draws of stats");
  XLALregSTRINGUserStruct ( outputAtoms,	 0,  UVAR_OPTIONAL, "Output F-statistic atoms into a file with this basename");
  XLALregSTRINGUserStruct ( outputInjParams,	 0,  UVAR_OPTIONAL,  "Output injection parameters into this file");

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

  UINT4 numDetectors = uvar->IFOs->length;
  UINT4 X;
  MultiLALDetector *multiDet;
  if ( (multiDet = XLALCreateMultiLALDetector ( numDetectors )) == NULL ) {
    XLALPrintError ("%s: XLALCreateMultiLALDetector(1) failed with errno=%d\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }
  LALDetector *site = NULL;
  for ( X=0; X < numDetectors; X++ )    {
    if ( (site = XLALGetSiteInfo ( uvar->IFOs->data[X] )) == NULL ) {
      XLALPrintError ("%s: Failed to get site-info for detector '%s'\n", __func__, uvar->IFOs->data[X] );
      XLAL_ERROR ( XLAL_EFUNC );
    }
    multiDet->data[X] = (*site); 	/* copy! */
    XLALFree ( site );
  }

  /* init timestamps vector covering observation time */
  UINT4 numSteps = (UINT4) ceil ( uvar->dataDuration / uvar->TAtom );
  MultiLIGOTimeGPSVector * multiTS;
  if ( (multiTS = XLALCreateMultiLIGOTimeGPSVector (numDetectors)) == NULL ) {
     XLALPrintError ("%s: XLALCreateMultiLIGOTimeGPSVector(%d) failed.\n", __func__, numDetectors );
  }

  for ( X=0; X < numDetectors; X++ )    {
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
  cfg->transientInjectRange.type = TRANSIENT_NONE;			/* default: no transient signal window */
  /* apply correct defaults if unset: t0=dataStart, t0Band=dataDuration-3*tauMax */
//   cfg->transientInjectRange.t0 = uvar->dataStartGPS + uvar->injectWindow_t0Days * DAY24;

  cfg->transientSearchRange = cfg->transientInjectRange;
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


/** Simple creator function for MultiLIGOTimeGPSVector with numDetectors entries
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


/** Write one line for given LV candidate into output file.
 *
 * NOTE: input dopplerParams can be NULL pointer, then just writes 0 in doppler fields (useful for synthetic LV draws)
 *
 */
int
write_LV_candidate_to_fp ( FILE *fp, const LVcomponents *LVstat, const PulsarDopplerParams *dopplerParams_in )
{
  /* sanity checks */
  if ( !fp ) {
    XLALPrintError ( "%s: invalid NULL filepointer input.\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( !LVstat || !LVstat->TwoFX || !LVstat->TwoFX->data ) {
    XLALPrintError ("\nError in function %s, line %d : Empty LVstat pointer as input parameter!\n\n", __func__, __LINE__);
    XLAL_ERROR ( XLAL_EFAULT);
  }

  PulsarDopplerParams dopplerParams = empty_PulsarDopplerParams;
  INIT_MEM( dopplerParams.fkdot );
  if ( dopplerParams_in == NULL ) { /* just write zeros */
    dopplerParams.Alpha = 0.0;
    dopplerParams.Delta = 0.0;
  }
  else {
    dopplerParams = *dopplerParams_in;
  }

  /* add output-field containing twoF and per-detector sumTwoFX */
  char statString[256] = "";	/* defaults to empty */
  char buf0[256];
  snprintf ( statString, sizeof(statString), " %.6f", LVstat->TwoF );
  UINT4 numDet = LVstat->TwoFX->length;
  UINT4 X;
  for ( X = 0; X < numDet ; X ++ ) {
    snprintf ( buf0, sizeof(buf0), " %.6f", LVstat->TwoFX->data[X] );
    UINT4 len1 = strlen ( statString ) + strlen ( buf0 ) + 1;
    if ( len1 > sizeof ( statString ) ) {
      XLALPrintError ("%s: assembled output string too long! (%d > %d)\n", __func__, len1, sizeof(statString));
      break;	/* we can't really terminate with error in this function, but at least we avoid crashing */
    }
    strcat ( statString, buf0 );
  } /* for X < numDet */
  snprintf ( buf0, sizeof(buf0), " %.6f", LVstat->LV );
  strcat ( statString, buf0 );
  fprintf (fp, "%.16g %.13g %.13g %.13g%s\n",
		dopplerParams.fkdot[0],
		dopplerParams.Alpha,
		dopplerParams.Delta,
		dopplerParams.fkdot[1],
		statString
	    );

  return XLAL_SUCCESS;

} /* write_LV_candidate_to_fp() */


/** Simple creator function for MultiInjParams with numDetectors entries */
MultiInjParams *
XLALCreateMultiInjParams ( UINT4 numDetectors )
{
  MultiInjParams *ret;

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

} /* XLALCreateMultiInjParams() */


/** Corresponding destructor function for MultiInjParams.
  * As usual this allows NULL input.
  */
void
XLALDestroyMultiInjParams ( MultiInjParams *multipar )
{
  if ( !multipar )
    return;

  if ( multipar->data )
    XLALFree ( multipar->data );

  XLALFree ( multipar );

  return;

} /* XLALDestroyMultiInjParams() */


InjParams_t *
XLALCombineInjParamsForLine( const MultiInjParams *injParamsX,  /**< array of the single-IFO injection parameters */
                      const UINT4 lineX                      /**< detector number where line was injected */
                      )
{

  /* check input parameters and report errors */
  if ( !injParamsX ) {
    XLALPrintError ("\nError in function %s, line %d : received empty injParamsX pointer!\n\n", __func__, __LINE__);
    XLAL_ERROR_NULL ( XLAL_EFAULT);
  }

  UINT4 numDetectors = injParamsX->length;
//   for ( UINT4 X=0; X < numDetectors; X ++ ) {
//     if ( !injParamsX->data[X] ) {
//       XLALPrintError ("\nError in function %s, line %d : injParams[%d] is empty!\n\n", __func__, __LINE__, X);
//       XLAL_ERROR_NULL ( XLAL_EFAULT);
//     }
//   }

  InjParams_t *multiInjParams = &empty_InjParams_t;

  /* pulsar parameters (skypos, amplitudes, transient stuff) can be taken from the single IFO with a line injection */
  multiInjParams->skypos = injParamsX->data[lineX].skypos;
  multiInjParams->ampParams = injParamsX->data[lineX].ampParams;
  for ( UINT4 i=0; i < 4; i++ )
    multiInjParams->ampVect[i] = injParamsX->data[lineX].ampVect[i];
  multiInjParams->transientWindow = injParamsX->data[lineX].transientWindow;

  /* SNR and Mmunu have to be summed over detectors */
  REAL8 rho2 = 0.0;
  for ( UINT4 X=0; X < numDetectors; X ++ ) {
    multiInjParams->M_mu_nu.Ad += injParamsX->data[X].M_mu_nu.Ad; /* multi-IFO M_mu_nu = sum_X M_mu_nu_X */
    multiInjParams->M_mu_nu.Bd += injParamsX->data[X].M_mu_nu.Bd;
    multiInjParams->M_mu_nu.Cd += injParamsX->data[X].M_mu_nu.Cd;
    multiInjParams->M_mu_nu.Sinv_Tsft += injParamsX->data[X].M_mu_nu.Sinv_Tsft; /* noise adds harmonically 1/S = sum_X (1/S_X) */
    rho2 += SQ(injParamsX->data[X].SNR); /* multi-IFO SNR^2 = sum_X SNR_X^2 */
  } /* for X < numDetectors */

  /* compute multi-SNR and -determinants */
  multiInjParams->M_mu_nu.Dd = multiInjParams->M_mu_nu.Ad * multiInjParams->M_mu_nu.Bd - SQ(multiInjParams->M_mu_nu.Cd); /* update sub-determinant */
  multiInjParams->SNR = sqrt(rho2); /* multi-IFO SNR^2 = sum_X SNR_X^2 */
  multiInjParams->detM1o8 = sqrt ( multiInjParams->M_mu_nu.Sinv_Tsft ) * pow ( multiInjParams->M_mu_nu.Dd, 0.25 ); /* (detM)^(1/8) = sqrt(Tsft/Sn) * (Dp)^(1/4) */

  return multiInjParams;
} /* XLALCombineInjParamsForLine */
