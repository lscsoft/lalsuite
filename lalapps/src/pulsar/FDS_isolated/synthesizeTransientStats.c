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
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

/* GSL includes */
#include <lal/LALGSL.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
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


/** Signal (amplitude) parameter ranges
 */
typedef struct {
  REAL8 h0Nat;		/**< h0 in *natural units* ie h0Nat = h0/sqrt(Sn) */
  REAL8 h0NatBand;	/**< draw h0Sn from Band [h0Sn, h0Sn + Band ] */
  REAL8 SNR;		/**< if > 0: alternative to h0Nat/h0NatBand: fix optimal signal SNR */
  REAL8 cosi;
  REAL8 cosiBand;
  REAL8 psi;
  REAL8 psiBand;
  REAL8 phi0;
  REAL8 phi0Band;
} AmpParamsRange_t;

/** Complete signal ranges to be considered for random-drawing
 * of signals
 */
typedef struct {
  SkyPosition skypos;	/**< (Alpha,Delta,system). Use Alpha < 0 to signal 'allsky' */
  AmpParamsRange_t ampRange;
} SignalParamsRange_t;

/** Configuration settings required for and defining a coherent pulsar search.
 * These are 'pre-processed' settings, which have been derived from the user-input.
 */
typedef struct {
  AmpParamsRange_t AmpRange;	/**< signal parameter ranges: lower bounds + bands */
  gsl_rng * rng;		/**< gsl random-number generator */
  LALDetector *site;		/**< detector site */
  EphemerisData *edat;		/**< ephemeris data */
  CHAR *version_string;		/**< code VCS version info */
} ConfigVariables;

/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

/* ----- User-variables: can be set from config-file or command-line */
typedef struct {
  BOOLEAN help;		/**< trigger output of help string */

  /* amplitude parameters + ranges */
  REAL8 h0;		/**< instantaneous GW amplitude h0 measured in units of sqrt(Sn): (h0/sqrt(Sn)) */
  REAL8 h0Band;		/**< randomize signal within [h0, h0+Band] with uniform prior */
  REAL8 SNR;		/**< specify fixed SNR: adjust h0 to obtain signal of this optimal SNR */
  REAL8 cosi;		/**< cos(inclination angle). If not set: randomize within [-1,1] */
  REAL8 psi;		/**< polarization angle psi. If not set: randomize within [-pi/4,pi/4] */
  REAL8 phi0;		/**< initial GW phase phi_0. If not set: randomize within [0, 2pi] */

  REAL8 Alpha;		/**< skyposition Alpha (RA) in radians */
  REAL8 Delta;		/**< skyposition Delta (Dec) in radians */

  CHAR *IFO;		/**< IFO name */

  REAL8 numDraws;	/**< number of random 'draws' to simulate for F-stat and B-stat */

  CHAR *outputStats;	/**< output file to write numDraw resulting statistics into */

  BOOLEAN SignalOnly;	/**< dont generate noise-draws: will result in non-random 'signal only' values of F and B */

  BOOLEAN version;	/**< output version-info */

} UserInput_t;


/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);

int XLALInitUserVars ( UserInput_t *uvar );
int XLALInitCode ( ConfigVariables *cfg, const UserInput_t *uvar );

/* exportable API */
int XLALDrawCorrelatedNoise ( gsl_vector *n_mu, const gsl_matrix *L, gsl_rng * rng );
FstatAtomVector* XLALSynthesizeFstatAtomVector4Noise ( const AMCoeffs *amcoeffs, gsl_rng * rng );
MultiFstatAtomVector* XLALSynthesizeMultiFstatAtomVector4Noise ( const MultiAMCoeffs *multiAM, REAL8 TAtom, gsl_rng * rng);

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

  if ( (cfg.version_string = XLALGetVersionString(lalDebugLevel)) == NULL ) {
    LogPrintf ( LOG_CRITICAL, "%s:XLALGetVersionString(%d) failed with errno=%d.\n", fn, lalDebugLevel, xlalErrno );
    return 1;
  }

  if ( uvar.version ) {
    printf ( "%s\n", cfg.version_string );
    return 0;
  }

  /* ---------- Initialize code-setup ---------- */
  if ( XLALInitCode( &cfg, &uvar ) != XLAL_SUCCESS ) {
    LogPrintf (LOG_CRITICAL, "%s: XLALInitCode() failed with error = %d\n", fn, xlalErrno );
    return 1;
  }

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

  uvar->phi0 = 0;
  uvar->psi = 0;

  /* register all our user-variables */
  XLALregBOOLUserStruct ( help, 		'h',     UVAR_HELP, "Print this message");
  /* signal amplitude parameters */
  XLALregREALUserStruct ( h0,			's', UVAR_OPTIONAL, "Overall GW amplitude h0, in natural units of sqrt{Sn}");
  XLALregREALUserStruct ( h0Band,	 	 0,  UVAR_OPTIONAL, "Randomize amplitude within [h0, h0+h0Band] with uniform prior");
  XLALregREALUserStruct ( SNR,			 0,  UVAR_OPTIONAL, "Alternative: adjust h0 to obtain signal of exactly this optimal SNR");

  XLALregREALUserStruct ( cosi,			'i', UVAR_OPTIONAL, "cos(inclination angle). If not set: randomize within [-1,1].");
  XLALregREALUserStruct ( psi,			 0,  UVAR_OPTIONAL, "polarization angle psi. If not set: randomize within [-pi/4,pi/4].");
  XLALregREALUserStruct ( phi0,		 	 0,  UVAR_OPTIONAL, "initial GW phase phi_0. If not set: randomize within [0, 2pi]");


  XLALregREALUserStruct ( numDraws,		'N', UVAR_OPTIONAL, "Number of random 'draws' to simulate for statistics");

  XLALregSTRINGUserStruct ( outputStats,	'o', UVAR_OPTIONAL, "Output file containing 'numDraws' random draws of stats");

  XLALregBOOLUserStruct ( SignalOnly,        	'S', UVAR_OPTIONAL, "Signal only: generate pure signal without noise");

  XLALregBOOLUserStruct ( version,        	'V', UVAR_SPECIAL,  "Output code version");

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

  /* ----- parse user-input on signal amplitude-paramters + ranges ----- */

  if ( XLALUserVarWasSet ( &uvar->SNR ) && ( XLALUserVarWasSet ( &uvar->h0 ) || XLALUserVarWasSet (&uvar->h0Band) ) )
    {
      LogPrintf (LOG_CRITICAL, "%s: specify only one of either {--h0,--h0Band} or --SNR\n", fn);
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }

  cfg->AmpRange.h0Nat = uvar->h0;
  cfg->AmpRange.h0NatBand = uvar->h0Band;
  cfg->AmpRange.SNR = uvar->SNR;

  /* implict ranges on cosi, psi and phi0 if not specified by user */
  if ( XLALUserVarWasSet ( &uvar->cosi ) )
    {
      cfg->AmpRange.cosi = uvar->cosi;
      cfg->AmpRange.cosiBand = 0;
    }
  else
    {
      cfg->AmpRange.cosi = -1;
      cfg->AmpRange.cosiBand = 2;
    }
  if ( XLALUserVarWasSet ( &uvar->psi ) )
    {
      cfg->AmpRange.psi = uvar->psi;
      cfg->AmpRange.psiBand = 0;
    }
  else
    {
      cfg->AmpRange.psi = - LAL_PI_4;
      cfg->AmpRange.psiBand = LAL_PI_2;
    }
  if ( XLALUserVarWasSet ( &uvar->phi0 ) )
    {
      cfg->AmpRange.phi0 = uvar->phi0;
      cfg->AmpRange.phi0Band = 0;
    }
  else
    {
      cfg->AmpRange.phi0 = 0;
      cfg->AmpRange.phi0Band = LAL_TWOPI;
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

  gsl_vector *normal;
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
  if ( (normal = gsl_vector_calloc ( 4 ) ) == NULL) {
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
  gsl_vector_free ( normal );

  return XLAL_SUCCESS;

} /* XLALDrawCorrelatedNoise() */

/** Generate an FstatAtomVector of pure noise for given antenna-pattern functions
 */
FstatAtomVector*
XLALSynthesizeFstatAtomVector4Noise ( const AMCoeffs *amcoeffs,	/**< input antenna-pattern functions {a_i, b_i} */
                                      gsl_rng * rng		/**< random-number generator */
                                      )
{
  const char *fn = __func__;

  /* check input consistency */
  if ( !amcoeffs || !amcoeffs->a || !amcoeffs->b ) {
    XLALPrintError ("%s: invalid NULL input in amcoeffs=%p or amcoeffs->a=%p, amcoeffs->b=%p\n", fn, amcoeffs, amcoeffs->a, amcoeffs->b );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }
  if ( !rng ) {
    XLALPrintError ("%s: invalid NULL input for random-number generator 'rng'\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  UINT4 numSFTs = amcoeffs->a->length;
  if ( numSFTs != amcoeffs->b->length ) {
    XLALPrintError ("%s: inconsistent lengths amcoeffs->a = %d, amecoeffs->b = %d\n", fn, amcoeffs->a->length, amcoeffs->b->length );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  /* prepare output vector */
  FstatAtomVector *atoms;
  if ( ( atoms = XLALCreateFstatAtomVector ( numSFTs ) ) == NULL ) {
    XLALPrintError ("%s: XLALCreateFstatAtomVector(%d) failed.\n", fn, numSFTs );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  /* prepare gsl-matrix for correlator L = 1/4 * [ a, a ; b , b ] */
  gsl_matrix *Lcor;
  if ( (Lcor = gsl_matrix_calloc ( 4, 4 )) == NULL ) {
    XLALPrintError ("%s: gsl_matrix_calloc ( 4, 4 ) failed.\n", fn );
    XLALDestroyFstatAtomVector ( atoms );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }
  /* prepare placeholder for 4 n_mu draws */
  gsl_vector *n_mu;
  if ( (n_mu = gsl_vector_calloc ( 4 ) ) == NULL ) {
    XLALPrintError ("%s: gsl_vector_calloc ( 4 ) failed.\n", fn );
    gsl_matrix_free ( Lcor );
    XLALDestroyFstatAtomVector ( atoms );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  /* ----- step through atoms and synthesize them ----- */
  UINT4 alpha;
  for ( alpha=0; alpha < numSFTs; alpha ++ )
    {
      REAL8 a = amcoeffs->a->data[alpha];
      REAL8 b = amcoeffs->b->data[alpha];

      REAL8 a4 = 0.25 * a;
      REAL8 b4 = 0.25 * b;
      /* upper-left block */
      gsl_matrix_set ( Lcor, 0, 0, a4 );
      gsl_matrix_set ( Lcor, 0, 1, a4 );
      gsl_matrix_set ( Lcor, 1, 0, b4 );
      gsl_matrix_set ( Lcor, 1, 1, b4 );
      /* lower-right block: +2 on all components */
      gsl_matrix_set ( Lcor, 2, 2, a4 );
      gsl_matrix_set ( Lcor, 2, 3, a4 );
      gsl_matrix_set ( Lcor, 3, 2, b4 );
      gsl_matrix_set ( Lcor, 3, 3, b4 );

      if ( XLALDrawCorrelatedNoise ( n_mu, Lcor, rng ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: failed to XLALDrawCorrelatedNoise().\n", fn );
        XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
      }
      REAL8 x1,x2,x3,x4;
      x1 = gsl_vector_get ( n_mu, 0 );
      x2 = gsl_vector_get ( n_mu, 1 );
      x3 = gsl_vector_get ( n_mu, 2 );
      x4 = gsl_vector_get ( n_mu, 3 );

      /* store this in Fstat-atom */
      atoms->data[alpha].a2_alpha = a * a;
      atoms->data[alpha].b2_alpha = b * b;
      atoms->data[alpha].ab_alpha = a * b;

      /* relation Fa,Fb <--> x_mu: see Eq.(72) in CFSv2-LIGO-T0900149-v2.pdf */
      atoms->data[alpha].Fa_alpha.re =   x1;
      atoms->data[alpha].Fa_alpha.im = - x3;

      atoms->data[alpha].Fb_alpha.re =   x2;
      atoms->data[alpha].Fb_alpha.im = - x4;

    } /* for i < numSFTs */

  /* free internal memory */
  gsl_vector_free ( n_mu );
  gsl_matrix_free ( Lcor );

  /* return result */
  return atoms;

} /* XLALSynthesizeFstatAtomVector4Noise() */


/** Generate an FstatAtomVector of pure noise for given antenna-pattern functions
 */
MultiFstatAtomVector*
XLALSynthesizeMultiFstatAtomVector4Noise ( const MultiAMCoeffs *multiAM,/**< input antenna-pattern functions {a_i, b_i} */
                                           REAL8 TAtom,			/**< atom time-base (typically Tsft) */
                                           gsl_rng * rng		/**< random-number generator */
                                           )
{
  const char *fn = __func__;

  /* check input consistency */
  if ( !multiAM || !multiAM->data || !multiAM->data[0] ) {
    XLALPrintError ("%s: invalid NULL input in 'mutiAM'\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }
  if ( !rng ) {
    XLALPrintError ("%s: invalid NULL input for random-number generator 'rng'\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  UINT4 numDetectors = multiAM->length;

  /* create output vector */
  MultiFstatAtomVector *multiAtoms;
  if ( (multiAtoms = XLALCalloc ( 1, sizeof(*multiAtoms) ) ) == NULL ) {
    XLALPrintError ("%s: XLALCalloc ( 1, %d) failed.\n", fn, sizeof(*multiAtoms) );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }
  if ( (multiAtoms->data = XLALCalloc ( numDetectors, sizeof(*multiAtoms->data) )) == NULL ) {
    XLALPrintError ("%s: XLALCalloc ( %, %d) failed.\n", fn, numDetectors, sizeof(*multiAtoms->data) );
    XLALFree ( multiAtoms );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  UINT4 X;
  for ( X=0; X < numDetectors; X ++ )
    {
      if ( ( multiAtoms->data[X] = XLALSynthesizeFstatAtomVector4Noise ( multiAM->data[X], rng )) == NULL ) {
        XLALPrintError ("%s: XLALSynthesizeFstatAtomVector4Noise() failed.\n", fn );
        XLALDestroyMultiFstatAtomVector ( multiAtoms );
        XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
      }
      multiAtoms->data[X]->TAtom = TAtom;
    } /* for X < numDetectors */

  /* return result */
  return multiAtoms;

} /* XLALSynthesizeMultiFstatAtomVector4Noise() */
