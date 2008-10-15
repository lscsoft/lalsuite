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

#define GSL_RANGE_CHECK_OFF

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

/** Configuration settings required for and defining a coherent pulsar search.
 * These are 'pre-processed' settings, which have been derived from the user-input.
 */
typedef struct {
  gsl_matrix *M_mu_nu;		/**< antenna-pattern matrix and normalization */
  gsl_matrix *M_Mu_Nu;		/**< matrix invsers M^{mu,nu} of M_{mu,nu} */
  gsl_vector  *A_Mu;		/**< signal amplitude vector */
} ConfigVariables;

/*---------- Global variables ----------*/
extern int vrbflg;		/**< defined in lalapps.c */

ConfigVariables GV;		/**< global container for various derived configuration settings */

/* ----- User-variables: can be set from config-file or command-line */
typedef struct {
  BOOLEAN help;		/**< trigger output of help string */

  REAL8 h0;		/**< overall GW amplitude h0 */
  REAL8 cosi;		/**< cos(inclination angle) */
  REAL8 psi;		/**< polarization angle psi */
  REAL8 phi0;		/**< initial GW phase phi_0 in radians */
  REAL8 Freq;		/**< GW signal frequency */
  REAL8 Alpha;		/**< sky-position angle 'alpha', which is right ascencion in equatorial coordinates */
  REAL8 Delta;		/**< sky-position angle 'delta', which is declination in equatorial coordinates */

  REAL8 M11;		/**< componentent {1,1} of M_{mu,nu}: T Sinv A */
  REAL8 M22;		/**< componentent {2,2} of M_{mu,nu}: T Sinv B */
  REAL8 M12;		/**< componentent {1,2} of M_{mu,nu}: T Sinv C */

  INT4 numDraws;	/**< number of random 'draws' to simulate for F-stat and B-stat */

  REAL8 numMCpoints;	/**< number of points to use for Monte-Carlo integration */

  CHAR *outputStats;	/**< output file to write F-stat estimation results into */

  BOOLEAN version;	/**< output version-info */
} UserInput_t;

static UserInput_t empty_UserInput;


typedef struct {
  double M11;
  double M22;
  double M12;
  const gsl_vector *x_mu;
} MCparams;


RCSID( "$Id$");

/* ---------- local prototypes ---------- */
int main(int argc,char *argv[]);

void initUserVars (LALStatus *status, UserInput_t *uvar );
void InitCode ( LALStatus *, ConfigVariables *cfg, const UserInput_t *uvar );
double BstatIntegrand ( double A[], size_t dim, void *p );
int XLALAmplitudeParams2Vect ( REAL8 Amu[], REAL8 h0, REAL8 cosi, REAL8 psi, REAL8 phi0 );

/*---------- empty initializers ---------- */


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

  gsl_matrix *M_chol;

  gsl_vector *s_mu, *FStat, *BStat, *dBStat, *lnL;
  gsl_matrix *normal, *x_mu_i;

  const gsl_rng_type * T;
  gsl_rng * rng;  /* gsl random-number generator */
  UINT4 i, row, col;
  int gslstat;
  REAL8 rho2;	/* optimal SNR^2 */

  lalDebugLevel = 0;
  vrbflg = 1;	/* verbose error-messages */

  /* turn off default GSL error handler */
  gsl_set_error_handler_off ();

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register all user-variable */
  LAL_CALL (LALGetDebugLevel(&status, argc, argv, 'v'), &status);
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

  /* Initialize code-setup */
  LAL_CALL ( InitCode(&status, &GV, &uvar ), &status);

  /* ----- compute Cholesky decomposition of M_mu_nu */
  if ( (M_chol = gsl_matrix_calloc ( 4, 4 ) ) == NULL) {
    LogPrintf ( LOG_CRITICAL, "Out of memory?\n");
    return 1;
  }
  {
    gsl_matrix *tmp;
    if ( (tmp = gsl_matrix_calloc ( 4, 4 ) ) == NULL) {
      LogPrintf ( LOG_CRITICAL, "Out of memory?\n");
      return 1;
    }
    if ( (gslstat = gsl_matrix_memcpy ( tmp, GV.M_mu_nu )) ) {
      LogPrintf ( LOG_CRITICAL, "gsl_matrix_memcpy() failed: %s\n", gsl_strerror (gslstat) );
      return 1;
    }
    /* Cholesky decompose M_mu_nu = L^T * L */
    if ( (gslstat = gsl_linalg_cholesky_decomp ( tmp ) ) ) {
      LogPrintf ( LOG_CRITICAL, "gsl_linalg_cholesky_decomp(M_mu_nu) failed: %s\n", gsl_strerror (gslstat) );
      return 1;
    }
    /* copy lower triangular matrix, which is L */
    for ( row = 0; row < 4; row ++ )
      for ( col = 0; col <= row; col ++ )
	gsl_matrix_set ( M_chol, row, col, gsl_matrix_get ( tmp, row, col ) );
    gsl_matrix_free ( tmp );
  } /* Cholesky-decompose M_mu_nu */


  /* compute signal components s_mu (simply A^mu with index lowered via M_mu_nu) */
  if ( (s_mu = gsl_vector_calloc ( 4 )) == NULL ) {
    LogPrintf ( LOG_CRITICAL, "Out of memory?\n");
    return 1;
  }
  /* GSL-doc: int gsl_blas_dsymv (CBLAS_UPLO_t Uplo, double alpha, const gsl_matrix * A,
   *                              const gsl_vector * x, double beta, gsl_vector * y )
   *
   * compute the matrix-vector product and sum: y = alpha A x + beta y
   * for the symmetric matrix A. Since the matrix A is symmetric only its
   * upper half or lower half need to be stored. When Uplo is CblasUpper
   * then the upper triangle and diagonal of A are used, and when Uplo
   * is CblasLower then the lower triangle and diagonal of A are used.
   */
  if ( (gslstat = gsl_blas_dsymv (CblasUpper, 1.0, GV.M_mu_nu, GV.A_Mu, 0.0, s_mu)) ) {
    LogPrintf ( LOG_CRITICAL, "gsl_blas_dsymv(M_mu_nu * A^mu failed): %s\n", gsl_strerror (gslstat) );
    return 1;
  }

  /* ----- generate 'numDraws' normal-distributed random numbers ----- */
  /* initialize gsl random-number generator */
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc (T);

  /*
     printf ("generator type: %s\n", gsl_rng_name (rng));
     printf ("seed = %lu\n", gsl_rng_default_seed);
  */
  /* printf ("first value = %lu\n", gsl_rng_get (rng)); */

  /* generate 'numDraws' random number with normal distribution: mean=0, sigma=1 */
  if ( (normal = gsl_matrix_calloc ( uvar.numDraws, 4 ) ) == NULL) {
    LogPrintf ( LOG_CRITICAL, "Out of memory?\n");
    return 1;
  }
  for ( i = 0; i < (UINT4)uvar.numDraws; i ++ )
    {
      gsl_matrix_set (normal, i, 0,  gsl_ran_gaussian ( rng, 1.0 ) );
      gsl_matrix_set (normal, i, 1,  gsl_ran_gaussian ( rng, 1.0 ) );
      gsl_matrix_set (normal, i, 2,  gsl_ran_gaussian ( rng, 1.0 ) );
      gsl_matrix_set (normal, i, 3,  gsl_ran_gaussian ( rng, 1.0 ) );
    }

  /* use normal-variates with Cholesky decomposed matrix to get n_mu with cov(n_mu,n_nu) = M_mu_nu */
  if ( ( x_mu_i = gsl_matrix_calloc ( uvar.numDraws, 4 ) ) == NULL ) {
    LogPrintf ( LOG_CRITICAL, "Out of memory?\n");
    return 1;
  }
  /*
  FILE *fpTemp = fopen ( "__tmp.dat", "wb" );
  */
  for ( i = 0; i < (UINT4)uvar.numDraws; i ++ )
    {
      gsl_vector_const_view normi = gsl_matrix_const_row ( normal, i );
      gsl_vector_view xi = gsl_matrix_row ( x_mu_i, i );

      gsl_matrix_set_row (x_mu_i, i, s_mu);	/* initialize x_mu = s_mu */

      /* int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)
       * compute the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H
       * for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
       */
      if ( (gslstat = gsl_blas_dgemv (CblasNoTrans, 1.0, M_chol, &(normi.vector), 1.0, &(xi.vector))) ) {
	LogPrintf ( LOG_CRITICAL, "gsl_blas_dgemv(M_chol * ni) failed: %s\n", gsl_strerror (gslstat) );
	return 1;
      }
      /*
      fprintf ( fpTemp, "%f  %f  %f  %f\n",
	gsl_matrix_get ( x_mu_i, i, 0 ),
	gsl_matrix_get ( x_mu_i, i, 1 ),
	gsl_matrix_get ( x_mu_i, i, 2 ),
	gsl_matrix_get ( x_mu_i, i, 3 ) );
      */
    }
  /*
  fclose ( fpTemp );
  */

  /* ---------- compute log likelihood ratio lnL ---------- */
  /* this assumes *known* A^Mu of course
   * computes lnL = A^mu x_mu - 1/2 A^mu M_mu_nu A^nu
   *              = A^mu x_mu - 1/2 A^mu s_mu
   */
  {
    if ( (lnL = gsl_vector_calloc ( uvar.numDraws ) ) == NULL) {
      LogPrintf ( LOG_CRITICAL, "Out of memory?\n");
      return 1;
    }

    /* STEP 1: compute A^mu x_mu (--> vector of numDraws scalars)*/

    /* Function: int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)
     *
     * These functions compute the matrix-vector product and sum y = \alpha op(A) x + \beta y,
     * where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
     */
    if ( (gslstat = gsl_blas_dgemv (CblasNoTrans, 1.0, x_mu_i, GV.A_Mu, 0.0, lnL)) ) {
      LogPrintf ( LOG_CRITICAL, "gsl_blas_dgemv(x_mu_i * A^mu) failed: %s\n", gsl_strerror (gslstat) );
      return 1;
    }

    /* STEP 2: compute the single scalar "optimal SNR^2": rho2 = A^mu M_mu_nu A^nu = A^mu s_mu */

    /* Function: int gsl_blas_ddot (const gsl_vector * x, const gsl_vector * y, double * result)
     *
     * These functions compute the scalar product x^T y for the vectors x and y, returning the result in result.
     */
    if ( (gslstat = gsl_blas_ddot (GV.A_Mu, s_mu, &rho2)) ) {
      LogPrintf ( LOG_CRITICAL, "rho2 = gsl_blas_ddot(A^mu * s_mu) failed: %s\n", gsl_strerror (gslstat) );
      return 1;
    }

    /* STEP 3: combine results for lnL = A^mu x_mu - 1/2 A^mu M_mu_nu A^nu */
    /* Function: int gsl_vector_add_constant (gsl_vector * a, const double x)
     *
     * This function adds the constant value x to the elements of the vector a, a'_i = a_i + x.
     */
    if ( ( gslstat = gsl_vector_add_constant (lnL, - 0.5 * rho2)) ) {
      LogPrintf ( LOG_CRITICAL, "gsl_vector_add_constant(lnL - 1/2 A^mu s_mu) failed: %s\n", gsl_strerror (gslstat) );
      return 1;
    }

  } /* compute lnL */

  /* ---------- compute F-stat values ---------- */
  {
    int sig;
    gsl_vector *x_Mu = gsl_vector_alloc ( 4 );
    gsl_permutation *perm = gsl_permutation_calloc ( 4 );
    gsl_matrix *Mmunu_LU = gsl_matrix_calloc ( 4, 4 );
    gsl_matrix_memcpy (Mmunu_LU, GV.M_mu_nu);

    if ( (FStat = gsl_vector_calloc ( uvar.numDraws ) ) == NULL) {
      LogPrintf ( LOG_CRITICAL, "Out of memory?\n");
      return 1;
    }

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
      LogPrintf ( LOG_CRITICAL, "gsl_linalg_LU_decomp (Mmunu) failed: %s\n", gsl_strerror (gslstat) );
      return 1;
    }

    for ( i=0; i < (UINT4)uvar.numDraws; i ++ )
      {
	gsl_vector_const_view xi = gsl_matrix_const_row ( x_mu_i, i );
	double x2;

	/* STEP 1: compute x^mu = M^{mu,nu} x_nu */

	/* Function: int gsl_linalg_LU_solve (const gsl_matrix * LU, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x)
	 *
	 * These functions solve the square system A x = b using the LU decomposition of A into (LU, p) given by
	 * gsl_linalg_LU_decomp or gsl_linalg_complex_LU_decomp.
	 */
	if ( (gslstat = gsl_linalg_LU_solve (Mmunu_LU, perm, &(xi.vector), x_Mu)) ) {
	  LogPrintf ( LOG_CRITICAL, "gsl_linalg_LU_solve (x^Mu = M^{mu,nu} x_nu) failed: %s\n", gsl_strerror (gslstat) );
	  return 1;
	}

	/* STEP 2: compute scalar product x_mu x^mu */

	/* Function: int gsl_blas_ddot (const gsl_vector * x, const gsl_vector * y, double * result)
	 *
	 * These functions compute the scalar product x^T y for the vectors x and y, returning the result in result.
	 */
	if ( (gslstat = gsl_blas_ddot (&(xi.vector), x_Mu, &x2)) ) {
	  LogPrintf ( LOG_CRITICAL, "i = %d: int gsl_blas_ddot (x_mu x^mu) failed: %s\n", i, gsl_strerror (gslstat) );
	  return 1;
	}

	/* write result into FStat (=2F) vector */
	gsl_vector_set ( FStat, i, x2 );

      } /* for i < numDraws */

    gsl_permutation_free ( perm );
    gsl_matrix_free ( Mmunu_LU );
    gsl_vector_free ( x_Mu );
  } /* compute F-stat draws */

  /* ---------- compute B-stat values ---------- */
  if ( (BStat = gsl_vector_calloc ( uvar.numDraws ) ) == NULL) {
    LogPrintf ( LOG_CRITICAL, "Out of memory?\n");
    return 1;
  }
  if ( (dBStat = gsl_vector_calloc ( uvar.numDraws ) ) == NULL) {
    LogPrintf ( LOG_CRITICAL, "Out of memory?\n");
    return 1;
  }

  {
    gsl_monte_vegas_state * MCS_vegas = gsl_monte_vegas_alloc ( 2 );
    gsl_monte_function F;
    MCparams pars;
    double prefact = 2.0 * 7.87480497286121;	/* 2 * sqrt(2) * pi^(3/2) */

    pars.M11 = uvar.M11;
    pars.M22 = uvar.M22;
    pars.M12 = uvar.M12;

    F.f = &BstatIntegrand;
    F.dim = 2;
    F.params = &pars;

    for ( i=0; i < (UINT4)uvar.numDraws; i ++ )
      {
	gsl_vector_const_view xi = gsl_matrix_const_row ( x_mu_i, i );
	double Bb;
	double Amp[2] = {0,0};
	double xlower[2], xupper[2];
	double abserr;
	pars.x_mu = &(xi.vector);

	Amp[0] = uvar.cosi;
	Amp[1] = uvar.psi;

	gsl_monte_vegas_init ( MCS_vegas );

	/* Integration boundaries */
	xlower[0] = -1;		/* cosi */
	xlower[1] = -LAL_PI_4;	/* psi */

	xupper[0] = 1;		/* cosi */
	xupper[1] = LAL_PI_4;	/* psi */

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
	if ( (gslstat = gsl_monte_vegas_integrate ( &F, xlower, xupper, 2, (size_t)uvar.numMCpoints, rng, MCS_vegas, &Bb, &abserr)) ) {
	  LogPrintf ( LOG_CRITICAL, "i = %d: gsl_monte_vegas_integrate() failed: %s\n", i, gsl_strerror (gslstat) );
	  return 1;
	}
	gsl_vector_set ( BStat, i, prefact * Bb );
	gsl_vector_set ( dBStat, i, abserr );
	/*
	  printf ( "Vegas Monte-Carlo integration: Bb = %f, abserr = %f ==> relerr = %f\n", Bb, abserr, abserr / Bb );
	*/

      } /* i < numDraws */

    gsl_monte_vegas_free ( MCS_vegas );
  } /* compute B-stat */

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
      fprintf (fpStat, "%%%% rho2 = %f\n", rho2 );
      fprintf (fpStat, "%%%% s_mu = ");
      XLALfprintfGSLvector ( fpStat, "%f", s_mu );
      fprintf (fpStat, "%%%% A_Mu = ");
      XLALfprintfGSLvector ( fpStat, "%f", GV.A_Mu );
      fprintf (fpStat, "%%%% x_1       x_2        x_3        x_4              lnL              2F             B_vegas         dB\n");
      for ( i=0; i < (UINT4)uvar.numDraws; i ++ )
	fprintf ( fpStat, "%10f %10f %10f %10f     %12f     %12f     %12f  %12f\n",
		  gsl_matrix_get ( x_mu_i, i, 0 ),
		  gsl_matrix_get ( x_mu_i, i, 1 ),
		  gsl_matrix_get ( x_mu_i, i, 2 ),
		  gsl_matrix_get ( x_mu_i, i, 3 ),

		  gsl_vector_get ( lnL, i ),
		  gsl_vector_get ( FStat, i ),
		  gsl_vector_get ( BStat, i ),
		  gsl_vector_get ( dBStat, i )
		  );

      fclose (fpStat);
    } /* if outputStat */

  /* Free config-Variables and userInput stuff */
  LAL_CALL (LALDestroyUserVars (&status), &status);

  gsl_matrix_free ( GV.M_mu_nu );
  gsl_vector_free ( s_mu );
  gsl_vector_free ( GV.A_Mu );
  gsl_vector_free ( lnL );
  gsl_vector_free ( FStat );
  gsl_matrix_free ( normal );
  gsl_matrix_free ( x_mu_i );
  gsl_vector_free ( BStat );
  gsl_vector_free ( dBStat );

  gsl_rng_free (rng);

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

  /* register all our user-variables */
  LALregBOOLUserStruct(status,	help, 		'h', UVAR_HELP,     "Print this message");

  LALregREALUserStruct(status,	h0,		's', UVAR_OPTIONAL, "Overall GW amplitude h0]");
  LALregREALUserStruct(status,	cosi,		'i', UVAR_OPTIONAL, "Inclination angle of rotation axis cos(iota)");

  LALregREALUserStruct(status,	psi,		'Y', UVAR_OPTIONAL, "Polarisation angle in radians");
  LALregREALUserStruct(status,	phi0,		'Y', UVAR_OPTIONAL, "Initial GW phase phi0 in radians");

  LALregREALUserStruct(status,	M11,	  	 0,  UVAR_REQUIRED, "Antenna-pattern matrix M_mu_nu: component {1,1} = T Sinv A");
  LALregREALUserStruct(status,	M22,	  	 0,  UVAR_REQUIRED, "Antenna-pattern matrix M_mu_nu: component {2,2} = T Sinv B");
  LALregREALUserStruct(status,	M12,	  	 0,  UVAR_REQUIRED, "Antenna-pattern matrix M_mu_nu: component {1,2} = T Sinv C");

  LALregINTUserStruct(status,	numDraws,	'N', UVAR_OPTIONAL, "Number of random 'draws' to simulate for F-stat and B-stat");
  LALregREALUserStruct(status,	numMCpoints,	'M', UVAR_OPTIONAL, "Number of points to use in Monte-Carlo integration");

  LALregSTRINGUserStruct(status, outputStats,	'o', UVAR_OPTIONAL, "Output filename containing random draws of lnL, 2F and B");

  LALregBOOLUserStruct(status,	version,        'V', UVAR_SPECIAL,   "Output code version");

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* initUserVars() */


/** Initialized Fstat-code: handle user-input and set everything up. */
void
InitCode ( LALStatus *status, ConfigVariables *cfg, const UserInput_t *uvar )
{
  INITSTATUS (status, "InitCode", rcsid);
  ATTATCHSTATUSPTR (status);

  /* ----- set up M_mu_nu matrix ----- */
  if ( ( cfg->M_mu_nu = gsl_matrix_calloc ( 4, 4 )) == NULL ) {
    LogPrintf (LOG_CRITICAL, "Seem to be out of memory!\n");
    ABORT ( status, PREDICTFSTAT_EMEM, PREDICTFSTAT_MSGEMEM );
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

  /* ----- generate signal amplitude vectors ---------- */
  if ( (cfg->A_Mu = gsl_vector_calloc ( 4 )) == NULL ) {
    LogPrintf ( LOG_CRITICAL, "Out of memory?\n");
    ABORT ( status, PREDICTFSTAT_EMEM, PREDICTFSTAT_MSGEMEM );
  }

  {
    double Amu[4];

    XLALAmplitudeParams2Vect ( Amu, uvar->h0, uvar->cosi, uvar->psi, uvar->phi0 );
    gsl_vector_set ( cfg->A_Mu, 0, Amu[0] );
    gsl_vector_set ( cfg->A_Mu, 1, Amu[1] );
    gsl_vector_set ( cfg->A_Mu, 2, Amu[2] );
    gsl_vector_set ( cfg->A_Mu, 3, Amu[3] );
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* InitPFS() */


/** compute log likelihood ratio lnL for given Amp = {h0, cosi, psi, phi0} and M_{mu,nu}.
 * computes lnL = A^mu x_mu - 1/2 A^mu M_mu_nu A^nu.
 *
 * This function is of type gsl_monte_function for gsl Monte-Carlo integration.
 * The parameter-struct is simply gsl_matrix *M_mu_nu.
 *
 */
double
BstatIntegrand ( double Amp[], size_t dim, void *p )
{
  MCparams *par = (MCparams *) p;
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

  printf ("%f   %f    %f\n", eta, psi, integrand );

  return integrand;

} /* BstatIntegrand() */


/** Convert amplitude params {h0, cosi, psi, phi0} into amplitude vector A^mu.
 *
 * NOTE: We rely on Amu[] being allocated 4-dim vector!
 */
int
XLALAmplitudeParams2Vect ( REAL8 Amu[], REAL8 h0, REAL8 cosi, REAL8 psi, REAL8 phi0 )
{
  REAL8 aPlus = 0.5 * h0 * ( 1.0 + SQ(cosi) );
  REAL8 aCross = h0 * cosi;
  REAL8 cos2psi = cos ( 2.0 * psi );
  REAL8 sin2psi = sin ( 2.0 * psi );
  REAL8 cosphi0 = cos ( phi0 );
  REAL8 sinphi0 = sin ( phi0 );

  if ( !Amu )
    return -1;

  Amu[0] =  aPlus * cos2psi * cosphi0 - aCross * sin2psi * sinphi0;
  Amu[1] =  aPlus * sin2psi * cosphi0 + aCross * cos2psi * sinphi0;
  Amu[2] = -aPlus * cos2psi * sinphi0 - aCross * sin2psi * cosphi0;
  Amu[3] = -aPlus * sin2psi * sinphi0 + aCross * cos2psi * cosphi0;

  return 0;

} /* XLALAmplitudeParams2Vect() */
