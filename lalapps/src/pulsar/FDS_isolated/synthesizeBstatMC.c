/*
 * Copyright (C) 2008 Reinhard Prix
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
 * Generate N samples of B-statistic (and F-statistic) values drawn from their
 * respective distributions, assuming Gaussian noise, for given signal parameters.
 *
 * This is mostly meant to be used for Monte-Carlos studies of ROC curves
 *
 *********************************************************************************/
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
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>


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

#include <lal/lalGitID.h>
#include <lalappsGitID.h>

#include <lalapps.h>

/*---------- DEFINES ----------*/

#define MAXFILENAMELENGTH 256   /* Maximum # of characters of a SFT filename */

#define EPHEM_YEARS  "00-04"	/**< default range: override with --ephemYear */

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

/** Signal (amplitude) parameter ranges
 */
typedef struct {
  REAL8 h0;
  REAL8 cosi;
  REAL8 psi;
  REAL8 phi0;
} AmpParams_t;

/** Configuration settings required for and defining a coherent pulsar search.
 * These are 'pre-processed' settings, which have been derived from the user-input.
 */
typedef struct {
  gsl_matrix *M_mu_nu;		/**< antenna-pattern matrix and normalization */
  AmpParams_t AmpMin;		/**< signal amplitude-parameters: lower bound on value-ranges */
  AmpParams_t AmpBand;		/**< signal ampltiude-parameters: band-widths on value-ranges */
  gsl_rng * rng;		/**< gsl random-number generator */
  /* signal + data vectors */
  gsl_matrix *A_Mu_i;		/**< list of 'numDraws' signal amplitude vectors {A^mu} */
  gsl_matrix *s_mu_i;		/**< list of 'numDraws' (covariant) signal amplitude vectors {s_mu = M_mu_nu A^nu} */
  gsl_matrix *x_mu_i;		/**< list of 'numDraws' data vectors {x_mu = s_mu + n_mu} */
  /* detection statistics */
  gsl_vector *lnL;		/**< list of 'numDraws' log-likelihood statistics */
  gsl_vector *Fstat;		/**< list of 'numDraws' F-statistics */
  gsl_vector *Bstat;		/**< list of 'numDraws' B-statistics */
} ConfigVariables;

/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

/* ----- User-variables: can be set from config-file or command-line */
typedef struct {
  BOOLEAN help;		/**< trigger output of help string */

  /* amplitude parameters + ranges */
  REAL8 h0;		/**< overall GW amplitude h0 */
  REAL8 h0Band;		/**< randomize signal within [h0, h0+Band] with uniform prior */
  REAL8 cosi;		/**< cos(inclination angle). If not set: randomize within [-1,1] */
  REAL8 psi;		/**< polarization angle psi. If not set: randomize within [-pi/4,pi/4] */
  REAL8 phi0;		/**< initial GW phase phi_0. If not set: randomize within [0, 2pi] */

  /* Doppler parameters + ranges */
  REAL8 Freq;		/**< GW signal frequency */
  REAL8 Alpha;		/**< sky-position angle 'alpha', which is right ascencion in equatorial coordinates */
  REAL8 Delta;		/**< sky-position angle 'delta', which is declination in equatorial coordinates */

  REAL8 M11;		/**< componentent {1,1} of M_{mu,nu}: T Sinv A */
  REAL8 M22;		/**< componentent {2,2} of M_{mu,nu}: T Sinv B */
  REAL8 M12;		/**< componentent {1,2} of M_{mu,nu}: T Sinv C */

  REAL8 numDraws;	/**< number of random 'draws' to simulate for F-stat and B-stat */

  REAL8 numMCpoints;	/**< number of points to use for Monte-Carlo integration */

  CHAR *outputStats;	/**< output file to write F-stat estimation results into */

  INT4 integrationMethod; /**< 0 = 2D Vegas Monte-Carlo, 1 = 2D Gauss-Kronod */

  BOOLEAN version;	/**< output version-info */
} UserInput_t;


typedef struct {
  double M11;
  double M22;
  double M12;
  const gsl_vector *x_mu;
  double cosi;		/**< only used for non-Monte-Carlo gsl-integration: value of cosi at which to integrate over psi */
} integrationParams_t;


RCSID( "$Id$");

/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);

void initUserVars (LALStatus *status, UserInput_t *uvar );
int InitCode ( ConfigVariables *cfg, const UserInput_t *uvar );

int XLALsynthesizeSignals ( gsl_matrix *A_Mu_i, gsl_matrix *s_mu_i, const gsl_matrix *M_mu_nu, AmpParams_t AmpMin, AmpParams_t AmpBand, gsl_rng * rng );
int XLALsynthesizeData ( gsl_matrix *x_mu_i, const gsl_matrix *M_mu_nu, const gsl_matrix *s_mu, gsl_rng * rng );

int XLALcomputeLogLikelihood ( gsl_vector *lnL, const gsl_matrix *A_Mu_i, const gsl_matrix *s_mu_i, const gsl_matrix *x_mu_i);
int XLALcomputeFstatistic ( gsl_vector *Fstat, const gsl_matrix *M_mu_nu, const gsl_matrix *x_mu_i );

int XLALcomputeBstatisticMC ( gsl_vector *Bstat, const gsl_matrix *M_mu_nu, const gsl_matrix *x_mu_i, gsl_rng * rng, UINT4 numMCpoints );
int XLALcomputeBstatisticGauss ( gsl_vector *Bstat, const gsl_matrix *M_mu_nu, const gsl_matrix *x_mu_i );

double BstatIntegrandOuter ( double cosi, void *p );
double BstatIntegrandInner ( double psi, void *p );
double BstatIntegrand ( double A[], size_t dim, void *p );

int XLALAmplitudeParams2Vect ( gsl_vector *A_Mu, REAL8 h0, REAL8 cosi, REAL8 psi, REAL8 phi0 );


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
  LALStatus status = blank_status;
  UserInput_t uvar = empty_UserInput;
  ConfigVariables GV = empty_ConfigVariables;		/**< various derived configuration settings */
  UINT4 i;

  lalDebugLevel = 0;
  vrbflg = 1;	/* verbose error-messages */

  /* turn off default GSL error handler */
  gsl_set_error_handler_off ();

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register all user-variable */
  LAL_CALL (LALGetDebugLevel(&status, argc, argv, 'v'), &status);
  LogSetLevel(lalDebugLevel);
  LAL_CALL (initUserVars(&status, &uvar), &status);

  /* do ALL cmdline and cfgfile handling */
  LAL_CALL (LALUserVarReadAllInput(&status, argc, argv), &status);

  if (uvar.help)	/* if help was requested, we're done here */
    exit (0);

  if ( uvar.version ) {
    printf ( "%s\n", lalGitID );
    printf ( "%s\n", lalappsGitID );
    return 0;
  }

  /* ---------- Initialize code-setup ---------- */
  if ( InitCode( &GV, &uvar ) ) {
    LogPrintf (LOG_CRITICAL, "InitCode() failed with error = %d\n", xlalErrno );
    return 1;
  }

  /* ---------- generate numDraws random draws of signals (A^mu, s_mu) */
  if ( XLALsynthesizeSignals( GV.A_Mu_i, GV.s_mu_i, GV.M_mu_nu, GV.AmpMin, GV.AmpBand, GV.rng ) ) {
    LogPrintf (LOG_CRITICAL, "XLALsynthesizeData() failed with error = %d\n", xlalErrno );
    return 1;
  }

  /* ---------- generate numDraws random draws of signal + noise ==> x_mu_i */
  if ( XLALsynthesizeData( GV.x_mu_i, GV.M_mu_nu, GV.s_mu_i, GV.rng ) ) {
    LogPrintf (LOG_CRITICAL, "XLALsynthesizeData() failed with error = %d\n", xlalErrno );
    return 1;
  }

  /* ---------- compute log likelihood ratio lnL ---------- */
  if ( XLALcomputeLogLikelihood ( GV.lnL, GV.A_Mu_i, GV.s_mu_i, GV.x_mu_i) ) {
    LogPrintf (LOG_CRITICAL, "XLALcomputeLogLikelihood() failed with error = %d\n", xlalErrno );
    return 1;
  }

  /* ---------- compute F-statistic ---------- */
  if ( XLALcomputeFstatistic ( GV.Fstat, GV.M_mu_nu, GV.x_mu_i ) ) {
    LogPrintf (LOG_CRITICAL, "XLALcomputeFstatistic() failed with error = %d\n", xlalErrno );
    return 1;
  }

  /* ---------- compute B-statistic ---------- */
  switch ( uvar.integrationMethod )
    {
    case 0:
      if ( XLALcomputeBstatisticGauss ( GV.Bstat, GV.M_mu_nu, GV.x_mu_i ) ) {
	LogPrintf (LOG_CRITICAL, "XLALcomputeBstatisticGauss() failed with error = %d\n", xlalErrno );
	return 1;
      }
      break;

    case 1:
      if ( XLALcomputeBstatisticMC ( GV.Bstat, GV.M_mu_nu, GV.x_mu_i, GV.rng, (UINT4)uvar.numMCpoints ) ) {
	LogPrintf (LOG_CRITICAL, "XLALcomputeBstatisticMC() failed with error = %d\n", xlalErrno );
	return 1;
      }
      break;

    default:
      LogPrintf (LOG_CRITICAL, "Sorry, --integrationMethod = %d not implemented!\n", uvar.integrationMethod );
      return 1;
      break;
    } /* switch integrationMethod */

  /* ---------- output F-statistic and B-statistic samples into file, if requested */
  if (uvar.outputStats)
    {
      FILE *fpStat = NULL;
      CHAR *logstr = NULL;
      CHAR *id1, *id2;

      if ( (fpStat = fopen (uvar.outputStats, "wb")) == NULL)
	{
	  LALPrintError ("\nError opening file '%s' for writing..\n\n", uvar.outputStats);
	  return (PREDICTFSTAT_ESYS);
	}

      /* log search-footprint at head of output-file */
      LAL_CALL( LALUserVarGetLog (&status, &logstr,  UVAR_LOGFMT_CMDLINE ), &status );

      fprintf(fpStat, "%%%% cmdline: %s\n", logstr );
      LALFree ( logstr );
      id1 = XLALClearLinebreaks ( lalGitID );
      id2 = XLALClearLinebreaks ( lalappsGitID );
      fprintf ( fpStat, "%%%% %s\n%%%%%s\n", id1, id2 );
      LALFree ( id1 ); LALFree ( id2 );

      /* append 'dataSummary' */
      fprintf (fpStat, "%%%% x_1       x_2        x_3        x_4              lnL              2F             Bstat\n");
      for ( i=0; i < GV.Bstat->size; i ++ )
	fprintf ( fpStat, "%10f %10f %10f %10f     %12f     %12f     %12f\n",
		  gsl_matrix_get ( GV.x_mu_i, i, 0 ),
		  gsl_matrix_get ( GV.x_mu_i, i, 1 ),
		  gsl_matrix_get ( GV.x_mu_i, i, 2 ),
		  gsl_matrix_get ( GV.x_mu_i, i, 3 ),

		  gsl_vector_get ( GV.lnL, i ),
		  gsl_vector_get ( GV.Fstat, i ),
		  gsl_vector_get ( GV.Bstat, i )
		  );

      fclose (fpStat);
    } /* if outputStat */

  /* Free config-Variables and userInput stuff */
  LAL_CALL (LALDestroyUserVars (&status), &status);

  gsl_matrix_free ( GV.M_mu_nu );
  gsl_matrix_free ( GV.s_mu_i );
  gsl_matrix_free ( GV.A_Mu_i );
  gsl_vector_free ( GV.lnL );
  gsl_vector_free ( GV.Fstat );
  gsl_matrix_free ( GV.x_mu_i );
  gsl_vector_free ( GV.Bstat );

  gsl_rng_free (GV.rng);

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

  INITSTATUS( status, "initUserVars", rcsid );
  ATTATCHSTATUSPTR (status);

  /* set a few defaults */
  uvar->help = FALSE;
  uvar->outputStats = NULL;

  uvar->phi0 = 0;
  uvar->psi = 0;

  uvar->numDraws = 1;
  uvar->numMCpoints = 1e4;

  uvar->integrationMethod = 0;	/* 2D Vegas MonteCarlo */

  /* register all our user-variables */
  LALregBOOLUserStruct(status,	help, 		'h', UVAR_HELP,     "Print this message");

  LALregREALUserStruct(status,	h0,		's', UVAR_OPTIONAL, "Overall GW amplitude h0");
  LALregREALUserStruct(status,	cosi,		'i', UVAR_OPTIONAL, "cos(inclination angle). If not set: randomize within [-1,1].");
  LALregREALUserStruct(status,	psi,		'Y', UVAR_OPTIONAL, "polarization angle psi. If not set: randomize within [-pi/4,pi/4].");
  LALregREALUserStruct(status,	phi0,		'Y', UVAR_OPTIONAL, "initial GW phase phi_0. If not set: randomize within [0, 2pi]");

  LALregREALUserStruct(status,	h0Band,		 0,  UVAR_OPTIONAL, "Randomize amplitude within [h0, h0+h0Band] with uniform prior");

  LALregREALUserStruct(status,	M11,	  	 0,  UVAR_REQUIRED, "Antenna-pattern matrix M_mu_nu: component {1,1} = T Sinv A");
  LALregREALUserStruct(status,	M22,	  	 0,  UVAR_REQUIRED, "Antenna-pattern matrix M_mu_nu: component {2,2} = T Sinv B");
  LALregREALUserStruct(status,	M12,	  	 0,  UVAR_REQUIRED, "Antenna-pattern matrix M_mu_nu: component {1,2} = T Sinv C");

  LALregREALUserStruct(status,	numDraws,	'N', UVAR_OPTIONAL, "Number of random 'draws' to simulate for F-stat and B-stat");

  LALregSTRINGUserStruct(status, outputStats,	'o', UVAR_OPTIONAL, "Output file containing 'numDraws' random draws of lnL, 2F and B");

  LALregINTUserStruct(status, integrationMethod,'m', UVAR_OPTIONAL, "2D Integration-method: 0=Gauss-Kronod, 1=Monte-Carlo(Vegas)");
  LALregREALUserStruct(status,	numMCpoints,	'M', UVAR_OPTIONAL, "Number of points to use in Monte-Carlo integration");

  LALregBOOLUserStruct(status,	version,        'V', UVAR_SPECIAL,   "Output code version");

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* initUserVars() */


/** Initialized Fstat-code: handle user-input and set everything up. */
int
InitCode ( ConfigVariables *cfg, const UserInput_t *uvar )
{
  const char *fn = "InitCode()";
  UINT4 numDraws;

  /* ----- parse user-input on signal amplitude-paramters + ranges ----- */
  /* explicit range on h0 */
  cfg->AmpMin.h0 = uvar->h0;
  cfg->AmpBand.h0 = uvar->h0Band;
  /* implict ranges on cosi, psi and phi0 if not specified by user */
  if ( LALUserVarWasSet ( &uvar->cosi ) )
    {
      cfg->AmpMin.cosi = uvar->cosi;
      cfg->AmpBand.cosi = 0;
    }
  else
    {
      cfg->AmpMin.cosi = -1;
      cfg->AmpBand.cosi = 2;
    }
  if ( LALUserVarWasSet ( &uvar->psi ) )
    {
      cfg->AmpMin.psi = uvar->psi;
      cfg->AmpBand.psi = 0;
    }
  else
    {
      cfg->AmpMin.psi = - LAL_PI_4;
      cfg->AmpBand.psi = LAL_PI_2;
    }
  if ( LALUserVarWasSet ( &uvar->phi0 ) )
    {
      cfg->AmpMin.phi0 = uvar->phi0;
      cfg->AmpBand.phi0 = 0;
    }
  else
    {
      cfg->AmpMin.phi0 = 0;
      cfg->AmpBand.phi0 = LAL_TWOPI;
    }

  /* ---------- allocate input/output-vectors of numDraws entries ---------- */
  numDraws = (UINT4)uvar->numDraws;	/* cast from REAL8 into UINT4 */

  /* ----- allocate signal amplitude vectors ---------- */
  if ( (cfg->A_Mu_i = gsl_matrix_calloc ( numDraws, 4 )) == NULL ) {
    LogPrintf ( LOG_CRITICAL, "%s: Out of memory.\n", fn);
    return XLAL_ENOMEM;
  }
  /* ----- signal components s_mu ----- */
  if ( (cfg->s_mu_i = gsl_matrix_calloc ( numDraws, 4 )) == NULL ) {
    LogPrintf ( LOG_CRITICAL, "%s: Out of memory.\n", fn);
    return XLAL_ENOMEM;
  }
  /* ----- data-components x_mu */
  if ( ( cfg->x_mu_i = gsl_matrix_calloc ( numDraws, 4 ) ) == NULL ) {
    LogPrintf ( LOG_CRITICAL, "Out of memory.\n");
    return XLAL_ENOMEM;
  }
  /* ----- log-likelihood */
  if ( (cfg->lnL = gsl_vector_calloc ( numDraws ) ) == NULL) {
    LogPrintf ( LOG_CRITICAL, "Out of memory?\n");
    return 1;
  }
  /* ----- F-statistic */
  if ( (cfg->Fstat = gsl_vector_calloc ( numDraws ) ) == NULL) {
    LogPrintf ( LOG_CRITICAL, "Out of memory?\n");
    return 1;
  }
  /* ----- B-statistic */
  if ( (cfg->Bstat = gsl_vector_calloc ( numDraws ) ) == NULL) {
    LogPrintf ( LOG_CRITICAL, "Out of memory?\n");
    return 1;
  }

  /* ----- set up M_mu_nu matrix ----- */
  if ( ( cfg->M_mu_nu = gsl_matrix_calloc ( 4, 4 )) == NULL ) {
    LogPrintf (LOG_CRITICAL, "%s: Seem to be out of memory!\n", fn);
    return PREDICTFSTAT_EMEM;
  }
  {
    gsl_matrix_set (cfg->M_mu_nu, 0, 0,   uvar->M11 );
    gsl_matrix_set (cfg->M_mu_nu, 1, 1,   uvar->M22 );
    gsl_matrix_set (cfg->M_mu_nu, 0, 1,   uvar->M12 );
    gsl_matrix_set (cfg->M_mu_nu, 1, 0,   uvar->M12 );

    /*
      gsl_matrix_set (cfg->M_mu_nu, 3, 0,   Ed );
      gsl_matrix_set (cfg->M_mu_nu, 1, 2,  -Ed );
    */

    gsl_matrix_set (cfg->M_mu_nu, 2, 2,   uvar->M11 );
    gsl_matrix_set (cfg->M_mu_nu, 3, 3,   uvar->M22 );
    gsl_matrix_set (cfg->M_mu_nu, 2, 3,   uvar->M12 );
    gsl_matrix_set (cfg->M_mu_nu, 3, 2,   uvar->M12 );

    /*
      gsl_matrix_set (cfg->M_mu_nu, 2, 1,  -Ed );
      gsl_matrix_set (cfg->M_mu_nu, 3, 2,   Ed );
    */
  }

  /* ----- initialize random-number generator ----- */
  /* read out environment variables GSL_RNG_xxx */
  gsl_rng_env_setup ();
  cfg->rng = gsl_rng_alloc (gsl_rng_default);

  LogPrintf ( LOG_DEBUG, "random-number generator type: %s\n", gsl_rng_name (cfg->rng));
  LogPrintf ( LOG_DEBUG, "seed = %lu\n", gsl_rng_default_seed );

  return 0;

} /* InitCode() */


/** Generate random signal draws with uniform priors in given bands  [h0, cosi, psi, phi0], and
 *  return list of 'numDraws' {s_mu} vectors.
 */
int
XLALsynthesizeSignals ( gsl_matrix *A_Mu_i,		/**< [OUT] list of numDraws 4D line-vectors {A^nu} */
			gsl_matrix *s_mu_i,		/**< [OUT] list of numDraws 4D line-vectors {s_mu = M_mu_nu A^nu} */
			const gsl_matrix *M_mu_nu,	/**< antenna-pattern matrix M_mu_nu */
			AmpParams_t AmpMin,		/**< signal amplitude-parameters: lower bound on value-ranges */
			AmpParams_t AmpBand,		/**< signal ampltiude-parameters: band-widths on value-ranges */
			gsl_rng * rng			/**< gsl random-number generator */
			)
{
  const char *fn = "XLALsynthesizeSignals()";
  UINT4 row, numDraws;

  REAL8 h0Min   = AmpMin.h0;
  REAL8 h0Max   = h0Min + AmpBand.h0;
  REAL8 cosiMin = AmpMin.cosi;
  REAL8 cosiMax = cosiMin + AmpBand.cosi;
  REAL8 psiMin  = AmpMin.psi;
  REAL8 psiMax  = psiMin + AmpBand.psi;
  REAL8 phi0Min = AmpMin.phi0;
  REAL8 phi0Max = phi0Min + AmpBand.phi0;
  REAL8 h0, cosi, psi, phi0;
  gsl_vector *A_Mu, *s_mu;
  int gslstat;

  /* ----- check input arguments ----- */
  if ( !M_mu_nu || M_mu_nu->size1 != M_mu_nu->size2 || M_mu_nu->size1 != 4 ) {
    LogPrintf ( LOG_CRITICAL, "%s: Invalid input, M_mu_nu must be a 4x4 matrix.", fn );
    return XLAL_EINVAL;
  }

  if ( !A_Mu_i || A_Mu_i->size2 != 4 ) {
    LogPrintf ( LOG_CRITICAL, "%s: Invalid input, A_Mu_i must be a numDrawsx4 matrix.", fn );
    return XLAL_EINVAL;
  }

  if ( (A_Mu = gsl_vector_alloc ( 4 )) == NULL ) {
    LogPrintf ( LOG_CRITICAL, "%s: Out of memory!\n", fn );
    return XLAL_ENOMEM;
  }

  if ( (s_mu = gsl_vector_alloc ( 4 )) == NULL ) {
    LogPrintf ( LOG_CRITICAL, "%s: Out of memory!\n", fn );
    return XLAL_ENOMEM;
  }

  numDraws = A_Mu_i->size1;

  for ( row = 0; row < numDraws; row ++ )
    {
      h0   = gsl_ran_flat ( rng, h0Min, h0Max );
      cosi = gsl_ran_flat ( rng, cosiMin, cosiMax );
      psi  = gsl_ran_flat ( rng, psiMin, psiMax );
      phi0 = gsl_ran_flat ( rng, phi0Min, phi0Max );

      XLALAmplitudeParams2Vect ( A_Mu, h0, cosi, psi, phi0 );

      /* GSL-doc: int gsl_blas_dsymv (CBLAS_UPLO_t Uplo, double alpha, const gsl_matrix * A,
       *                              const gsl_vector * x, double beta, gsl_vector * y )
       *
       * compute the matrix-vector product and sum: y = alpha A x + beta y
       * for the symmetric matrix A. Since the matrix A is symmetric only its
       * upper half or lower half need to be stored. When Uplo is CblasUpper
       * then the upper triangle and diagonal of A are used, and when Uplo
       * is CblasLower then the lower triangle and diagonal of A are used.
       */
      if ( (gslstat = gsl_blas_dsymv (CblasUpper, 1.0, M_mu_nu, A_Mu, 0.0, s_mu)) ) {
	LogPrintf ( LOG_CRITICAL, "%s: gsl_blas_dsymv(M_mu_nu * A^mu failed): %s\n", fn, gsl_strerror (gslstat) );
	return XLAL_EDOM;
      }

      gsl_matrix_set_row ( A_Mu_i, row, A_Mu );
      gsl_matrix_set_row ( s_mu_i, row, s_mu );

    } /* row < numDraws */

  gsl_vector_free ( A_Mu );
  gsl_vector_free ( s_mu );

  return 0;

} /* XLALsynthesizeSignals() */


/** Generate random-noise draws and combine with (FIXME: single!) signal.
 *  Returns a list of numDraws vectors {x_mu}
 */
int
XLALsynthesizeData ( gsl_matrix *x_mu_i,		/**< [OUT] list of numDraws 4D line-vectors {x_mu = n_mu + s_mu} */
		     const gsl_matrix *M_mu_nu,		/**< 4x4 antenna-pattern matrix */
		     const gsl_matrix *s_mu_i,		/**< numDraws x 4D vector of 'signal components' s_mu = (s|h_mu) */
		     gsl_rng * rng			/**< gsl random-number generator */
		     )
{
  const CHAR *fn = "XLALsynthesizeData()";
  gsl_matrix *tmp, *M_chol;
  UINT4 row, col;
  gsl_matrix *normal;
  int gslstat;
  UINT4 numDraws;

  /* ----- check input arguments ----- */
  if ( !M_mu_nu || M_mu_nu->size1 != M_mu_nu->size2 || M_mu_nu->size1 != 4 ) {
    LogPrintf ( LOG_CRITICAL, "%s: Invalid input, M_mu_nu must be a 4x4 matrix.", fn );
    return XLAL_EINVAL;
  }

  if ( !x_mu_i || (x_mu_i->size2 != 4) ) {
    LogPrintf ( LOG_CRITICAL, "%s: Invalid input, x_mu_i must be a numDrawsx4 matrix.", fn );
    return XLAL_EINVAL;
  }

  numDraws = x_mu_i->size1;

  if ( !s_mu_i || (s_mu_i->size2 != 4)|| (s_mu_i->size1 != numDraws) ) {
    LogPrintf ( LOG_CRITICAL, "%s: Invalid input, s_mu_i must be a numDrawsx4 matrix.", fn );
    return XLAL_EINVAL;
  }

  /* ----- Cholesky decompose M_mu_nu = L^T * L ----- */
  if ( (M_chol = gsl_matrix_calloc ( 4, 4 ) ) == NULL) {
    LogPrintf ( LOG_CRITICAL, "Out of memory?\n");
    return XLAL_ENOMEM;
  }
  if ( (tmp = gsl_matrix_calloc ( 4, 4 ) ) == NULL) {
    LogPrintf ( LOG_CRITICAL, "Out of memory?\n");
    return XLAL_ENOMEM;
  }
  if ( (gslstat = gsl_matrix_memcpy ( tmp, M_mu_nu )) ) {
    LogPrintf ( LOG_CRITICAL, "gsl_matrix_memcpy() failed: %s\n", gsl_strerror (gslstat) );
    return XLAL_EDOM;
  }
  if ( (gslstat = gsl_linalg_cholesky_decomp ( tmp ) ) ) {
    LogPrintf ( LOG_CRITICAL, "gsl_linalg_cholesky_decomp(M_mu_nu) failed: %s\n", gsl_strerror (gslstat) );
    return XLAL_EDOM;
  }
  /* copy lower triangular matrix, which is L */
  for ( row = 0; row < 4; row ++ )
    for ( col = 0; col <= row; col ++ )
      gsl_matrix_set ( M_chol, row, col, gsl_matrix_get ( tmp, row, col ) );

  /* ----- generate 'numDraws' normal-distributed random numbers ----- */
  if ( (normal = gsl_matrix_calloc ( numDraws, 4 ) ) == NULL) {
    LogPrintf ( LOG_CRITICAL, "Out of memory?\n");
    return XLAL_ENOMEM;
  }

  for ( row = 0; row < numDraws; row ++ )
    {
      gsl_matrix_set (normal, row, 0,  gsl_ran_gaussian ( rng, 1.0 ) );
      gsl_matrix_set (normal, row, 1,  gsl_ran_gaussian ( rng, 1.0 ) );
      gsl_matrix_set (normal, row, 2,  gsl_ran_gaussian ( rng, 1.0 ) );
      gsl_matrix_set (normal, row, 3,  gsl_ran_gaussian ( rng, 1.0 ) );
    } /* for row < numDraws */

  /* ----- initialize x_mu = s_mu ----- */
  if ( (gslstat = gsl_matrix_memcpy (x_mu_i, s_mu_i)) ) {
    LogPrintf ( LOG_CRITICAL, "%s: gsl_matrix_memcpy() failed: %s\n", gsl_strerror (gslstat) );
    return XLAL_EDOM;
  }

  /* use normal-variates with Cholesky decomposed matrix to get n_mu with cov(n_mu,n_nu) = M_mu_nu */
  for ( row = 0; row < numDraws; row ++ )
    {
      gsl_vector_const_view normi = gsl_matrix_const_row ( normal, row );
      gsl_vector_view xi = gsl_matrix_row ( x_mu_i, row );

      /* int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)
       * compute the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H
       * for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
       */
      if ( (gslstat = gsl_blas_dgemv (CblasNoTrans, 1.0, M_chol, &(normi.vector), 1.0, &(xi.vector))) ) {
	LogPrintf ( LOG_CRITICAL, "gsl_blas_dgemv(M_chol * ni) failed: %s\n", gsl_strerror (gslstat) );
	return 1;
      }
    } /* for row < numDraws */

  /* ---------- free memory ---------- */
  gsl_matrix_free ( tmp );
  gsl_matrix_free ( M_chol );
  gsl_matrix_free ( normal );

  return XLAL_SUCCESS;

} /* XLALsynthesizeData() */



/** Compute log-likelihood function for given input data
 */
int
XLALcomputeLogLikelihood ( gsl_vector *lnL,		/**< [OUT] log-likelihood vector */
			   const gsl_matrix *A_Mu_i,	/**< 4D amplitude-vector (FIXME: numDraws) */
			   const gsl_matrix *s_mu_i,	/**< 4D signal-component vector s_mu = (s|h_mu) [FIXME] */
			   const gsl_matrix *x_mu_i	/**< numDraws x 4D data-vectors x_mu */
			   )
{
  const char *fn = "XLALcomputeLogLikelihood()";
  int gslstat;
  UINT4 row, numDraws;
  gsl_vector_view A_Mu, d_mu;
  gsl_matrix *tmp;
  REAL8 res;

  /* ----- check input arguments ----- */
  if ( !lnL || !A_Mu_i || !s_mu_i || !x_mu_i ) {
    LogPrintf ( LOG_CRITICAL, "%s: illegal NULL input vector passed.\n", fn);
    return XLAL_EINVAL;
  }
  numDraws = lnL->size;

  if ( (A_Mu_i->size1 != numDraws) || (A_Mu_i->size2 != 4) ) {
    LogPrintf ( LOG_CRITICAL, "%s: input Amplitude-vector A^mu must be numDraws x 4D.\n", fn);
    return XLAL_EINVAL;
  }
  if ( (s_mu_i->size1 != numDraws) || (s_mu_i->size2 != 4) ) {
    LogPrintf ( LOG_CRITICAL, "%s: input Amplitude-vector s_mu must be numDraws x 4D.\n", fn);
    return XLAL_EINVAL;
  }
  if ( (x_mu_i->size1 != numDraws) || (x_mu_i->size2 != 4) ) {
    LogPrintf ( LOG_CRITICAL, "%s: input vector-list x_mu_i must be numDraws x 4.\n", fn);
    return XLAL_EINVAL;
  }

  if ( (tmp = gsl_matrix_alloc ( numDraws, 4 ) ) == NULL ) {
    LogPrintf ( LOG_CRITICAL, "%s: Out of memory.\n", fn);
    return XLAL_ENOMEM;
  }

  /* STEP1: compute tmp_mu = x_mu - 0.5 s_mu */
  gsl_matrix_memcpy ( tmp, s_mu_i );
  gsl_matrix_scale ( tmp, - 0.5 );
  gsl_matrix_add ( tmp, x_mu_i );


  /* STEP2: compute A^mu tmp_mu */
  for ( row=0; row < numDraws; row ++ )
    {
      A_Mu = gsl_matrix_row (A_Mu_i, row);
      d_mu = gsl_matrix_row (tmp, row);

      /* Function: int gsl_blas_ddot (const gsl_vector * x, const gsl_vector * y, double * result)
       * These functions compute the scalar product x^T y for the vectors x and y, returning the result in result.
       */
      if ( (gslstat = gsl_blas_ddot (&A_Mu.vector, &d_mu.vector, &res)) ) {
	LogPrintf ( LOG_CRITICAL, "%s: lnL = gsl_blas_ddot(A^mu * (x_mu - 0.5 s_mu) failed: %s\n", fn, gsl_strerror (gslstat) );
	return 1;
      }
      gsl_vector_set ( lnL, row, res );

    } /* for row < numDraws */

  gsl_matrix_free ( tmp );

  return 0;

} /* XLALcomputeLogLikelihood() */


/** Compute F-statistic for given input data
 */
int
XLALcomputeFstatistic ( gsl_vector *Fstat,		/**< [OUT] F-statistic vector */
			const gsl_matrix *M_mu_nu,	/**< antenna-pattern matrix M_mu_nu */
			const gsl_matrix *x_mu_i	/**< data-vectors x_mu: numDraws x 4 */
			)
{
  const char *fn = "XLALcomputeFstatistic()";

  int sig;
  gsl_vector *x_Mu = gsl_vector_alloc ( 4 );
  gsl_permutation *perm = gsl_permutation_calloc ( 4 );
  gsl_matrix *Mmunu_LU = gsl_matrix_calloc ( 4, 4 );
  int gslstat;
  UINT4 row, numDraws;

  /* ----- check input arguments ----- */
  if ( !Fstat || !M_mu_nu || !x_mu_i ) {
    LogPrintf ( LOG_CRITICAL, "%s: illegal NULL input vector passed.\n", fn);
    return XLAL_EINVAL;
  }
  numDraws = Fstat->size;
  if ( (M_mu_nu->size1 != 4) || (M_mu_nu->size2 != 4) ) {
    LogPrintf ( LOG_CRITICAL, "%s: antenna-pattern matrix M_mu_nu must be 4x4.\n", fn);
    return XLAL_EINVAL;
  }
  if ( (x_mu_i->size1 != numDraws) || (x_mu_i->size2 != 4) ) {
    LogPrintf ( LOG_CRITICAL, "%s: input vector-list x_mu_i must be numDraws x 4.\n", fn);
    return XLAL_EINVAL;
  }

  gsl_matrix_memcpy (Mmunu_LU, M_mu_nu);

  /* Function: int gsl_linalg_LU_decomp (gsl_matrix * A, gsl_permutation * p, int * signum)
   *
   * These functions factorize the square matrix A into the LU decomposition PA = LU.
   * On output the diagonal and upper triangular part of the input matrix A contain the matrix U. The lower
   * triangular part of the input matrix (excluding the diagonal) contains L. The diagonal elements of L are
   * unity, and are not stored. The permutation matrix P is encoded in the permutation p. The j-th column of
   * the matrix P is given by the k-th column of the identity matrix, where k = p_j the j-th element of the
   * permutation vector. The sign of the permutation is given by signum. It has the value (-1)^n, where n is
   * the number of interchanges in the permutation.
   * The algorithm used in the decomposition is Gaussian Elimination with partial pivoting
   * (Golub & Van Loan, Matrix Computations, Algorithm 3.4.1).
   */
  if( (gslstat = gsl_linalg_LU_decomp (Mmunu_LU, perm, &sig)) ) {
    LogPrintf ( LOG_CRITICAL, "%s: gsl_linalg_LU_decomp (Mmunu) failed: %s\n", fn, gsl_strerror (gslstat) );
    return 1;
  }

  for ( row=0; row < numDraws; row ++ )
    {
      gsl_vector_const_view xi = gsl_matrix_const_row ( x_mu_i, row );
      double x2;

      /* STEP 1: compute x^mu = M^{mu,nu} x_nu */

      /* Function: int gsl_linalg_LU_solve (const gsl_matrix * LU, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x)
       *
       * These functions solve the square system A x = b using the LU decomposition of A into (LU, p) given by
       * gsl_linalg_LU_decomp or gsl_linalg_complex_LU_decomp.
       */
      if ( (gslstat = gsl_linalg_LU_solve (Mmunu_LU, perm, &(xi.vector), x_Mu)) ) {
	LogPrintf ( LOG_CRITICAL, "%s: gsl_linalg_LU_solve (x^Mu = M^{mu,nu} x_nu) failed: %s\n", fn, gsl_strerror (gslstat) );
	return 1;
      }

      /* STEP 2: compute scalar product x_mu x^mu */

      /* Function: int gsl_blas_ddot (const gsl_vector * x, const gsl_vector * y, double * result)
       *
       * These functions compute the scalar product x^T y for the vectors x and y, returning the result in result.
       */
      if ( (gslstat = gsl_blas_ddot (&(xi.vector), x_Mu, &x2)) ) {
	LogPrintf ( LOG_CRITICAL, "%s: row = %d: int gsl_blas_ddot (x_mu x^mu) failed: %s\n", fn, row, gsl_strerror (gslstat) );
	return 1;
      }

      /* write result into Fstat (=2F) vector */
      gsl_vector_set ( Fstat, row, x2 );

    } /* for row < numDraws */

  gsl_permutation_free ( perm );
  gsl_matrix_free ( Mmunu_LU );
  gsl_vector_free ( x_Mu );

  return 0;

} /* XLALcomputeFstatistic () */


/** Compute the B-statistic for given input data, using Monte-Carlo integration for
 * the marginalization over {cosi, psi}, while {h0, phi0} have been marginalized analytically.
 *
 * Currently uses the Vegas Monte-Carlo integrator, which samples more densely where the integrand is larger.
 */
int
XLALcomputeBstatisticMC ( gsl_vector *Bstat,		/**< [OUT] vector of numDraws B-statistic values */
			  const gsl_matrix *M_mu_nu,	/**< antenna-pattern matrix M_mu_nu */
			  const gsl_matrix *x_mu_i,	/**< data-vectors x_mu: numDraws x 4 */
			  gsl_rng * rng,		/**< gsl random-number generator */
			  UINT4 numMCpoints		/**< number of points to use in Monte-Carlo integration */
			  )
{
  const char *fn = "XLALcomputeBstatisticMC()";

  gsl_monte_vegas_state * MCS_vegas = gsl_monte_vegas_alloc ( 2 );
  gsl_monte_function F;
  integrationParams_t pars;
  double prefact = 2.0 * 7.87480497286121;	/* 2 * sqrt(2) * pi^(3/2) */
  UINT4 row, numDraws;
  int gslstat;

  /* ----- check input arguments ----- */
  if ( !Bstat || !M_mu_nu || !x_mu_i || !rng) {
    LogPrintf ( LOG_CRITICAL, "%s: illegal NULL input vector passed.\n", fn);
    return XLAL_EINVAL;
  }
  numDraws = Bstat->size;
  if ( (M_mu_nu->size1 != 4) || (M_mu_nu->size2 != 4) ) {
    LogPrintf ( LOG_CRITICAL, "%s: antenna-pattern matrix M_mu_nu must be 4x4.\n", fn);
    return XLAL_EINVAL;
  }
  if ( (x_mu_i->size1 != numDraws) || (x_mu_i->size2 != 4) ) {
    LogPrintf ( LOG_CRITICAL, "%s: input vector-list x_mu_i must be numDraws x 4.\n", fn);
    return XLAL_EINVAL;
  }

  /* ----- prepare Monte-Carlo integrator ----- */
  pars.M11 = gsl_matrix_get ( M_mu_nu, 0, 0 );
  pars.M22 = gsl_matrix_get ( M_mu_nu, 1, 1 );
  pars.M12 = gsl_matrix_get ( M_mu_nu, 0, 1 );

  F.f = &BstatIntegrand;
  F.dim = 2;
  F.params = &pars;

  for ( row=0; row < numDraws; row ++ )
    {
      gsl_vector_const_view xi = gsl_matrix_const_row ( x_mu_i, row );
      double Bb;
      double AmpLower[2], AmpUpper[2];
      double abserr;
      pars.x_mu = &(xi.vector);

      gsl_monte_vegas_init ( MCS_vegas );

      /* Integration boundaries */
      AmpLower[0] = -1;		/* cosi */
      AmpUpper[0] =  1;		/* cosi */

      AmpLower[1] = -LAL_PI_4;	/* psi */
      AmpUpper[1] =  LAL_PI_4;	/* psi */

      /* Function: int gsl_monte_vegas_integrate (gsl_monte_function * f, double * xl, double * xu, size_t dim, size_t calls,
       *                                          gsl_rng * r, gsl_monte_vegas_state * s, double * result, double * abserr)
       *
       * This routines uses the vegas Monte Carlo algorithm to integrate the function f over the dim-dimensional hypercubic
       * region defined by the lower and upper limits in the arrays xl and xu, each of size dim. The integration uses a
       * fixed number of function calls calls, and obtains random sampling points using the random number generator r.
       * A previously allocated workspace s must be supplied. The result of the integration is returned in result, with
       * an estimated absolute error abserr. The result and its error estimate are based on a weighted average of independent
       * samples. The chi-squared per degree of freedom for the weighted average is returned via the state struct component,
       * s->chisq, and must be consistent with 1 for the weighted average to be reliable.
       */
      if ( (gslstat = gsl_monte_vegas_integrate ( &F, AmpLower, AmpUpper, 2, numMCpoints, rng, MCS_vegas, &Bb, &abserr)) ) {
	LogPrintf ( LOG_CRITICAL, "%s: row = %d: gsl_monte_vegas_integrate() failed: %s\n", fn, row, gsl_strerror (gslstat) );
	return 1;
      }
      gsl_vector_set ( Bstat, row, prefact * Bb );

    } /* row < numDraws */

  gsl_monte_vegas_free ( MCS_vegas );

  return 0;

} /* XLALcomputeBstatisticMC() */


/** Compute the B-statistic for given input data, using standard Gauss-Kronod integration (gsl_integration_qng)
 * for the marginalization over {cosi, psi}, while {h0, phi0} have been marginalized analytically.
 *
 */
int
XLALcomputeBstatisticGauss ( gsl_vector *Bstat,		/**< [OUT] vector of numDraws B-statistic values */
			     const gsl_matrix *M_mu_nu,	/**< antenna-pattern matrix M_mu_nu */
			     const gsl_matrix *x_mu_i	/**< data-vectors x_mu: numDraws x 4 */
			     )
{
  const char *fn = "XLALcomputeBstatisticGauss()";

  double prefact = 2.0 * 7.87480497286121;	/* 2 * sqrt(2) * pi^(3/2) */
  UINT4 row, numDraws;
  int gslstat;
  double epsabs = 0;
  double epsrel = 1e-2;
  double abserr;
  gsl_function F;
  integrationParams_t pars;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  /* ----- check input arguments ----- */
  if ( !Bstat || !M_mu_nu || !x_mu_i ) {
    LogPrintf ( LOG_CRITICAL, "%s: illegal NULL input vector passed.\n", fn);
    return XLAL_EINVAL;
  }
  numDraws = Bstat->size;
  if ( (M_mu_nu->size1 != 4) || (M_mu_nu->size2 != 4) ) {
    LogPrintf ( LOG_CRITICAL, "%s: antenna-pattern matrix M_mu_nu must be 4x4.\n", fn);
    return XLAL_EINVAL;
  }
  if ( (x_mu_i->size1 != numDraws) || (x_mu_i->size2 != 4) ) {
    LogPrintf ( LOG_CRITICAL, "%s: input vector-list x_mu_i must be numDraws x 4.\n", fn);
    return XLAL_EINVAL;
  }

  /* ----- prepare Gauss-Kronod integrator ----- */
  pars.M11 = gsl_matrix_get ( M_mu_nu, 0, 0 );
  pars.M22 = gsl_matrix_get ( M_mu_nu, 1, 1 );
  pars.M12 = gsl_matrix_get ( M_mu_nu, 0, 1 );

  F.function = &BstatIntegrandOuter;
  F.params = &pars;

  for ( row=0; row < numDraws; row ++ )
    {
      gsl_vector_const_view xi = gsl_matrix_const_row ( x_mu_i, row );
      double Bb;
      double CosiLower, CosiUpper;

      pars.x_mu = &(xi.vector);

      /* Integration boundaries */
      CosiLower = -1;
      CosiUpper =  1;

      /* Function: int gsl_integration_qags (const gsl_function * f, double a, double b, double epsabs, double epsrel,
       *                                     size_t limit, gsl_integration_workspace * workspace, double * result, double * abserr)
       *
       * This function applies the Gauss-Kronrod 21-point integration rule adaptively until an estimate of the integral
       * of f over (a,b) is achieved within the desired absolute and relative error limits, epsabs and epsrel. The results
       * are extrapolated using the epsilon-algorithm, which accelerates the convergence of the integral in the presence of
       * discontinuities and integrable singularities. The function returns the final approximation from the extrapolation,
       * result, and an estimate of the absolute error, abserr. The subintervals and their results are stored in the memory
       * provided by workspace. The maximum number of subintervals is given by limit, which may not exceed the allocated size of the workspace.
       */
      if ( (gslstat = gsl_integration_qags ( &F, CosiLower, CosiUpper, epsabs, epsrel, 1000, w, &Bb, &abserr)) ) {
	LogPrintf ( LOG_CRITICAL, "%s: row = %d: gsl_integration_qag() failed: res=%f, abserr=%f, intervals=%d, %s\n",
		    fn, row, Bb, abserr, w->size, gsl_strerror (gslstat) );
	return 1;
      }

      gsl_vector_set ( Bstat, row, prefact * Bb );

    } /* row < numDraws */

  gsl_integration_workspace_free (w);

  return 0;

} /* XLALcomputeBstatisticGauss() */


/** log likelihood ratio lnL marginalized over {h0, phi0} (analytical) and integrated over psi in [-pi/4,pi/4], for
 * given cosi: BstatIntegrandOuter(cosi) = int lnL dh0 dphi0 dpsi
 *
 * This function is of type gsl_function for gsl-integration over cosi
 *
 */
double
BstatIntegrandOuter ( double cosi, void *p )
{
  const char *fn = "BstatIntegrandOuter()";
  integrationParams_t *par = (integrationParams_t *) p;
  gsl_function F;
  double epsabs = 0;
  double epsrel = 1e-3;
  double abserr;
  double ret;
  double PsiLower, PsiUpper;
  int gslstat;
  static gsl_integration_workspace * w = NULL;

  if ( !w )
    w = gsl_integration_workspace_alloc (1000);

  par->cosi = cosi;
  F.function = &BstatIntegrandInner;
  F.params = p;

  /* Integration boundaries */
  PsiLower = -LAL_PI_4;
  PsiUpper =  LAL_PI_4;

  /* Function: int gsl_integration_qags (const gsl_function * f, double a, double b, double epsabs, double epsrel,
   *                                     size_t limit, gsl_integration_workspace * workspace, double * result, double * abserr)
   *
   * This function applies the Gauss-Kronrod 21-point integration rule adaptively until an estimate of the integral
   * of f over (a,b) is achieved within the desired absolute and relative error limits, epsabs and epsrel. The results
   * are extrapolated using the epsilon-algorithm, which accelerates the convergence of the integral in the presence of
   * discontinuities and integrable singularities. The function returns the final approximation from the extrapolation,
   * result, and an estimate of the absolute error, abserr. The subintervals and their results are stored in the memory
   * provided by workspace. The maximum number of subintervals is given by limit, which may not exceed the allocated size of the workspace.
   */
  if ( (gslstat = gsl_integration_qags ( &F, PsiLower, PsiUpper, epsabs, epsrel, 1000, w, &ret, &abserr)) ) {
    LogPrintf ( LOG_CRITICAL, "%s: gsl_integration_qag() failed: res=%f, abserr=%f, intervals=%d, %s\n",
		fn, ret, abserr, w->size, gsl_strerror (gslstat) );
    return 1;
  }

  return ret;

} /* BstatIntegrandOuter() */


/** log likelihood ratio lnL marginalized over {h0, phi0} (analytical) for given psi and pars->cosi
 * BstatIntegrandInner(cosi,psi) = int lnL dh0 dphi0
 *
 * This function is of type gsl_function for gsl-integration over psi at fixed cosi,
 * and represents a simple wrapper around BstatIntegrand() for gsl-integration.
 *
 */
double
BstatIntegrandInner ( double psi, void *p )
{
  const char *fn = "BstatIntegrandInner()";
  integrationParams_t *par = (integrationParams_t *) p;
  double Amp[2], ret;

  Amp[0] = par->cosi;
  Amp[1] = psi;

  ret = BstatIntegrand ( Amp, 2, p );

  return ret;

} /* BstatIntegrandInner() */


/** compute log likelihood ratio lnL for given Amp = {h0, cosi, psi, phi0} and M_{mu,nu}.
 * computes lnL = A^mu x_mu - 1/2 A^mu M_mu_nu A^nu.
 *
 * This function is of type gsl_monte_function for gsl Monte-Carlo integration.
 *
 */
double
BstatIntegrand ( double Amp[], size_t dim, void *p )
{
  integrationParams_t *par = (integrationParams_t *) p;
  double x1, x2, x3, x4;
  double al1, al2, al3, al4;
  double eta, etaSQ, etaSQp1SQ;
  double psi, sin2psi, cos2psi, sin2psiSQ, cos2psiSQ;
  double AMA, qSQ, arg0;
  double integrand;

  if ( dim != 2 ) {
    LogPrintf (LOG_CRITICAL, "Error: BstatIntegrand() was called with illegal dim = %d != 2\n", dim );
    abort ();
  }

  /* introduce a few handy shortcuts */
  x1 = gsl_vector_get ( par->x_mu, 0 );
  x2 = gsl_vector_get ( par->x_mu, 1 );
  x3 = gsl_vector_get ( par->x_mu, 2 );
  x4 = gsl_vector_get ( par->x_mu, 3 );

  eta = Amp[0];
  etaSQ = SQ(eta);			/* eta^2 */
  etaSQp1SQ = SQ ( (1.0 + etaSQ) );	/* (1+eta^2)^2 */

  psi = Amp[1];
  sin2psi = sin ( 2.0 * psi );
  cos2psi = cos ( 2.0 * psi );
  sin2psiSQ = SQ(sin2psi);
  cos2psiSQ = SQ(cos2psi);

  /* compute amplitude-params alpha1, alpha2, alpha3 and alpha4 */
  al1 = 0.25 * etaSQp1SQ * cos2psiSQ + etaSQ * sin2psiSQ;
  al2 = 0.25 * etaSQp1SQ * sin2psiSQ + etaSQ * cos2psiSQ;
  al3 = 0.25 * SQ( (1.0 - etaSQ) ) * sin2psi * cos2psi;
  al4 = eta * ( 1.0 + etaSQ );

  /* STEP 1: compute AMA = A^mu M_mu_nu A^nu / h0^2 */
  AMA = al1 * par->M11 + al2 * par->M22 + 2.0 * al3 * par->M12;

  /* STEP2: compute q^2 */
  qSQ = al1 * ( SQ(x1) + SQ(x3) ) + al2 * ( SQ(x2) + SQ(x4) ) + 2.0 * al3 * ( x1 * x2 + x3 * x4 ) + al4 * ( x1 * x4 - x2 * x3 );

  /* STEP3 : put all the pieces together */
  arg0 = 0.25 * qSQ  / AMA;

  integrand = pow(AMA, -0.5) * exp(arg0) * gsl_sf_bessel_I0(arg0);

  if ( lalDebugLevel >= 2 )
    printf ("%f   %f    %f   %f %f\n", eta, psi, integrand, AMA, arg0 );

  return integrand;

} /* BstatIntegrand() */


/** Convert amplitude params {h0, cosi, psi, phi0} into amplitude vector A^mu.
 *
 * NOTE: We rely on Amu[] being allocated 4-dim vector!
 */
int
XLALAmplitudeParams2Vect ( gsl_vector *A_Mu, REAL8 h0, REAL8 cosi, REAL8 psi, REAL8 phi0 )
{
  REAL8 aPlus = 0.5 * h0 * ( 1.0 + SQ(cosi) );
  REAL8 aCross = h0 * cosi;
  REAL8 cos2psi = cos ( 2.0 * psi );
  REAL8 sin2psi = sin ( 2.0 * psi );
  REAL8 cosphi0 = cos ( phi0 );
  REAL8 sinphi0 = sin ( phi0 );

  if ( !A_Mu || (A_Mu->size != 4) )
    return -1;

  gsl_vector_set ( A_Mu, 0,  aPlus * cos2psi * cosphi0 - aCross * sin2psi * sinphi0 );
  gsl_vector_set ( A_Mu, 1,  aPlus * sin2psi * cosphi0 + aCross * cos2psi * sinphi0 );
  gsl_vector_set ( A_Mu, 2, -aPlus * cos2psi * sinphi0 - aCross * sin2psi * cosphi0 );
  gsl_vector_set ( A_Mu, 3, -aPlus * sin2psi * sinphi0 + aCross * cos2psi * cosphi0 );

  return 0;

} /* XLALAmplitudeParams2Vect() */

