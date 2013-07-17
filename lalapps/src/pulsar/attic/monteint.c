/*
*  Copyright (C) 2007 Yousuke Itoh
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

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>


/* gsl-include */
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

/* LAL-includes */
#include <lal/AVFactories.h>
#include <lal/UserInput.h>
#include <lal/LALDatatypes.h>
#include <lal/DetectorSite.h>
#include <lal/LALBarycenter.h>
#include <lal/LALComputeAM.h>
#include <lal/ComputeSky.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LALGSL.h>
#include <lal/LALError.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdlib.h>
#include <lal/StringInput.h>
#include <lal/ConfigFile.h>


/* LALApps-includes */
#include <lalapps/lalapps.h>

/* Definitions */
/*
#define USE_S2ALLSKY
#define TEST_MONTEINT
#define TEST_MARGINALIZEPDF
#define TEST_INTERPOLATE2D
*/
#define EPHEM_YEARS  "00-04"



#if ( defined(TEST_MONTEINT) || defined(TEST_MARGINALIZEPDF) || defined(TEST_INTERPOLATE2D) )
#define TEST_CODE
#endif


/* Pointer to function */
typedef REAL8 (*FindRootFunction)(REAL8 input, void *Params);
typedef REAL8 (*MonteCarloIntegrand)(REAL8 *inputs, size_t Dimension, void *Params);


/* Structures */
typedef struct {
  REAL8 Fstat;
  REAL8 tsft;			/**< length of an SFT in seconds */
  LALDetector Detector;         /**< Our detector*/
  CHAR *timestampsfile;
  CHAR *EphemEarth; 	/**< filename of earth-ephemeris data */
  CHAR *EphemSun;     /**< filename of sun-ephemeris data */
  REAL8 confidence;
} ConfigVariables;


typedef struct {
  INT4 verboseflag;
  size_t warmupcalls;
  size_t maincalls;
  INT4 maxVegasStages;
  const CHAR *method;
} MonteCarloIntParams;

typedef struct {
  MonteCarloIntegrand integrand;
  size_t dimension;
  void *params;
  REAL8 *domain;
} MonteCarloIntIn;

typedef struct {
  REAL8 result;
  REAL8 sigma;
} MonteCarloIntOut;


typedef struct {
  REAL8 result;
  REAL8 lower;
  REAL8 upper;
} findRootOut;


typedef struct {
  FindRootFunction function;
  void *params;
  REAL8 lower;
  REAL8 upper;
  REAL8 tolerance;
} findRootParams;


typedef struct {
  REAL8 gridsize;
  UINT4Vector *dimlen;
  REAL8Vector *alpha;
  REAL8Vector *delta;
  REAL8Array *As;
  REAL8Array *Bs;
  REAL8Array *Cs;
} ABCcoeffs;


typedef struct {
  REAL8 Fstat;
  REAL8 confidence;
  LIGOTimeGPS *midTS;
  BarycenterInput baryinput;
  AMCoeffsParams *amParams;
  AMCoeffs amc;
  ABCcoeffs *ABCs;
  REAL8 perr; /* error in p value*/
  REAL8 norm; /* norm of the marginalized pdf */
  REAL8 normsigma; /* error in norm */
} BayesFstatParams;



/* global uvar_variables */
BOOLEAN uvar_help;
BOOLEAN uvar_useLUT;
CHAR *uvar_outputLabel;
CHAR* uvar_IFO;
REAL8 uvar_Tsft;
REAL8 uvar_Fstat;
REAL8 uvar_confidence;
CHAR *uvar_ephemDir;
CHAR *uvar_ephemYear;
CHAR *uvar_timestampsfile;



/* global LALStatus for initialization */
LALStatus global_status;


/*----------------------------------------------------------------------*/
/* Error-codes */

#define BAYESFSTATC_ENULL 		1
#define BAYESFSTATC_ESYS     		2
#define BAYESFSTATC_EINPUT   		3
#define BAYESFSTATC_EMEM   		4


#define BAYESFSTATC_MSGENULL 		"Arguments contained an unexpected null pointer"
#define BAYESFSTATC_MSGESYS		"System call failed (probably file IO)"
#define BAYESFSTATC_MSGEINPUT   	"Invalid input"
#define BAYESFSTATC_MSGEMEM   	        "Memory allocation failure or Out of memory. Bad."



/*----------------------------------------------------------------------*/
/* Exit values */
#define BAYESFSTAT_EXIT_OK              0  /* normal exit */
#define BAYESFSTAT_EXIT_USAGE           7  /* user requested help */
#define BAYESFSTAT_EXIT_NOMEM          15  /* out of memory */




/*----------------------------------------------------------------------*/
/* Function declaration */
/*----------------------------------------------------------------------*/
void ABCLookUpTable(LALStatus *, ABCcoeffs *, BayesFstatParams *, REAL8 gridsize, UINT4 Ndelta, UINT4 Nalpha);
void AllocateABCcoeffs(LALStatus *, ABCcoeffs *,UINT4 Ndelta, UINT4 Nalpha);

void computeNorm(REAL8 *theNorm, REAL8 *theNormsigma, MonteCarloIntIn *, MonteCarloIntParams *);
REAL8 confminuscdf(REAL8 theStrain, void *theParams);
REAL8 DFunction(REAL8 *theVars, size_t theDimension, void *theParams);
REAL8 DInverseWithIntMeasure(REAL8 *theVars, size_t theDimension, void *theParams);
REAL8 DSquareWithIntMeasure(REAL8 *k, size_t dim, void *params);
void filelength(LALStatus *, INT8 *length, CHAR *filename);
void findRoot(LALStatus *, findRootOut *theOutputs, findRootParams *theFrparams);
void FreeABCcoeffs(LALStatus *, ABCcoeffs *);
void FreeBayesFstatParams(LALStatus *, BayesFstatParams *);
void FreeConfigVariables(LALStatus *, ConfigVariables *);
void FreeMem(LALStatus *, ConfigVariables *, BayesFstatParams *);
void InitAMParams(LALStatus *, BayesFstatParams *, ConfigVariables *);
void InitBayesFstatParams(LALStatus *, BayesFstatParams *, INT4 theArgc, CHAR *theArgv[]);
void InitConfigVariables(LALStatus *, ConfigVariables *);
void InitUserVars(LALStatus *);
void integrateDSquare(LALStatus *status, BayesFstatParams *bfparams);  
void Interpolate2DSphere(REAL8 *finterp, REAL8Vector *xvec, REAL8Vector *yvec, REAL8Array *func, REAL8 xval, REAL8 yval, REAL8 gridsize);
REAL8 marginalizedNCPDF(REAL8 theStrain, BayesFstatParams *);
void MonteCarloIntegrate(MonteCarloIntOut *, MonteCarloIntIn *, MonteCarloIntParams *);
REAL8 ncpdfWithIntMeasure(REAL8 *theVars, size_t theDimension, void *theParams);
REAL8 NormFunction(REAL8 theFstat);		    
void ReadTimeStamps(LALStatus *, LIGOTimeGPS *timestamps, INT8 ntimestamps, CHAR *timestampsfile); 
void testCode(LALStatus *theStatus, BayesFstatParams *theBfparams); 
void testInterPolate2DSphere(LALStatus *status, BayesFstatParams *bfparams);
void testInterPolate2DSphereCore(LALStatus *, BayesFstatParams *, REAL8 alpha, REAL8 delta);








/*----------------------------------------------------------------------*/
/* CODE starts here */
/*----------------------------------------------------------------------*/


/*----------------------------------------------------------------------
 *
 *
 *
 *
 *
 *
 *----------------------------------------------------------------------*/

extern INT4 vrbflg;

INT4 
main (INT4 argc, CHAR *argv[])
{
  LALStatus *status = &global_status;

  BayesFstatParams  *bfparams = NULL;
  findRootOut frout = {0, 0, 0};
  findRootParams frparams;
  REAL8 pvalue;

  lal_errhandler = LAL_ERR_EXIT;
  vrbflg = 1; /* verbose flag (see lalapps.h) */


  bfparams = (BayesFstatParams *) LALMalloc(sizeof(BayesFstatParams));
  LAL_CALL( InitBayesFstatParams(status, bfparams, argc, argv), status);


#ifdef TEST_CODE
  LAL_CALL( testCode(status, bfparams), status); 
#endif

#ifndef TEST_CODE
  frparams.lower = 0.8;
  frparams.upper = 1.5;
  frparams.tolerance = 0.01; /* 1 % error */
  frparams.function  = confminuscdf;
  frparams.params = (void *) bfparams;


  LAL_CALL( findRoot(status, &frout, &frparams ), status);

  pvalue = marginalizedNCPDF(frout.result,bfparams);

  if (lalDebugLevel>=1) {
    {REAL8 plower, pupper;
    plower = marginalizedNCPDF(frout.lower, bfparams);
    pupper = marginalizedNCPDF(frout.upper, bfparams);
        
    printf("%% Target p value, F stat, Obtained p value," 
	   "1 sgima lower bound in p, 1 sigma upper bound in p, estimated error in p,"
	   " normalized strain, 1 sigma lower bound in strain, 1 sigma upper bound in strain \n");
    
    fprintf(stdout,"%10.6g %10.6g %10.6g %10.6g %10.6g %10.6g %10.6g %10.6g %10.6g\n", 
	    bfparams->confidence, bfparams->Fstat, pvalue, plower, pupper, bfparams->perr, frout.result, frout.lower, frout.upper);
    }
  } else {
    {REAL8 magic = 6.0; /* Numeric normalization factor to make 95 upper limits on the normalized strain around 1. 
		     See REAL8 ncpdfWithIntMeasure() */
    fprintf(stdout,"%10.6g %10.6g %10.6g\n", /* F stat, p-value, normalized strain h0\sqrt{T_0/S_h}. */
	    bfparams->Fstat, pvalue, magic * sqrt( 2*bfparams->Fstat ) * frout.result );
    }
  }


#endif
  fflush(stdout);
 

  LAL_CALL( FreeBayesFstatParams(status, bfparams), status);
  LAL_CALL( LALDestroyUserVars(status), status);
  /* did we forget anything ? */
  LALCheckMemoryLeaks();

  return status->statusCode;
} /* INT4 main() */


/*----------------------------------------------------------------------
 *FUNCTION void testCode().
 *
 * + test the code.
 *
 *
 *
 *----------------------------------------------------------------------*/
void testCode(LALStatus *status, BayesFstatParams *bfparams) 
{

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  if ( bfparams == NULL ) {
    ABORT (status, BAYESFSTATC_ENULL, BAYESFSTATC_MSGENULL);
  }

#ifdef TEST_MONTEINT
  TRY( integrateDSquare(status->statusPtr, bfparams), status);
#endif

#ifdef TEST_MARGINALIZEPDF
  {REAL8 strain, pvalue = 0.0;
  strain = 1.0;
  pvalue = marginalizedNCPDF(strain, bfparams);
  printf("%% Obtained p value,  normalized strain, and F stat\n %g \t %g \t %g\n",
	  pvalue, strain, bfparams->Fstat);
  }
#endif
#ifdef TEST_INTERPOLATE2D
  TRY( testInterPolate2DSphere(status->statusPtr, bfparams), status);
#endif

  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* void testCode(). */




/*----------------------------------------------------------------------
 *FUNCTION void findRoot(). From the GSL example.
 *
 * + Find root.
 *
 *
 *
 *----------------------------------------------------------------------*/
void 
findRoot (LALStatus *status, findRootOut *frout, findRootParams *frparams)
{
  INT4 gsl_status = 0;
  INT4 iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T = NULL;
  gsl_root_fsolver *s = NULL;
  REAL8 r = 0;
  gsl_function F;
  REAL8 x_lo, x_hi;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  F.function = frparams->function;
  F.params = frparams->params;

  x_lo = frparams->lower;
  x_hi = frparams->upper;


  TRYGSL( T = gsl_root_fsolver_brent, status);
  TRYGSL( s = gsl_root_fsolver_alloc (T), status);
  TRYGSL( gsl_root_fsolver_set (s, &F, x_lo, x_hi), status);

  if(lalDebugLevel>=2) {
    printf ("%% using %s method\n", 
	    gsl_root_fsolver_name (s) );
    
    printf ("%% %5s [%9s, %9s] %9s %9s\n",
	    "iter", "lower", "upper", "root", 
	    "err(est)");
  }

  do
    {
      iter++;
      TRYGSL( gsl_status = gsl_root_fsolver_iterate (s), status);
      TRYGSL( r = gsl_root_fsolver_root (s), status);
      TRYGSL( x_lo = gsl_root_fsolver_x_lower (s), status);
      TRYGSL( x_hi = gsl_root_fsolver_x_upper (s), status);
      TRYGSL( gsl_status = gsl_root_test_interval (x_lo, x_hi,
                                       0, frparams->tolerance), status);

      if(lalDebugLevel>=2) {
	if (gsl_status == GSL_SUCCESS)
	  printf ("%%Converged:\n");
      }

      if(lalDebugLevel>=2) {
	printf ("%% %5d [%.7f, %.7f] %.7f %.7f\n",
		iter, x_lo, x_hi,
		r, 
		x_hi - x_lo);
      }
    }
  while (gsl_status == GSL_CONTINUE && iter < max_iter);

  if(lalDebugLevel>=2) {
    printf ("%% %5d  %12.7f  %12.7f  %12.7f  %12.7f\n",
	    iter, x_lo, x_hi,
	    r, 
	    x_hi - x_lo);
  }

  frout->result=r;
  frout->lower=x_lo;
  frout->upper=x_hi;


  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* void findRoot () */



/*----------------------------------------------------------------------
 *FUNCTION void confminuscdf(). 
 *
 * + function to be root-finded.
 *
 *
 *
 *----------------------------------------------------------------------*/
REAL8 confminuscdf(REAL8 strain, void *params)
{
  BayesFstatParams *bfparams;
  bfparams = (BayesFstatParams *) params;

  return (bfparams->confidence) -  marginalizedNCPDF(strain, bfparams);
}  /* void confminuscdf().  */ 



/*----------------------------------------------------------------------
 *FUNCTION void MonteCarloIntegrate(). From the GSL web-site example.
 *
 * 
 *
 *
 *
 *----------------------------------------------------------------------*/

void MonteCarloIntegrate(MonteCarloIntOut *mciout, MonteCarloIntIn *mciin, MonteCarloIntParams *mciparams)
{
  REAL8 res, err;
  REAL8 *xl,*xu;
  INT4 ic;
  const gsl_rng_type *T;

  size_t warmupcalls = 10000;
  size_t calls = 500000;
  INT4 verboseflag = 0;
  INT4 maxVegasStages = 100;

  gsl_rng *r;
  gsl_monte_function G;
  gsl_monte_vegas_state *s;

  gsl_rng_env_setup ();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);


  G.f = mciin->integrand;
  G.dim = mciin->dimension;
  G.params =  mciin->params;

  warmupcalls = mciparams->warmupcalls;
  calls = mciparams->maincalls;
  verboseflag = mciparams->verboseflag;
  maxVegasStages = mciparams->maxVegasStages;

  if(verboseflag >= 3) {
    fprintf(stderr,"%% \n%% We are using the %s method.\n",mciparams->method);
    fprintf(stderr,"%% maximum vegas stages %d.\n",maxVegasStages);
    fprintf(stderr,"%% number of warm-up calls %zu.\n",warmupcalls);
    fprintf(stderr,"%% number of main calls %zu.\n",calls);
  }

  xl = (REAL8 *) malloc(mciin->dimension*sizeof(REAL8)); 
  xu = (REAL8 *) malloc(mciin->dimension*sizeof(REAL8)); 
  for(ic=0;ic< (INT4) mciin->dimension;ic++) {
    xl[ic] = mciin->domain[2*ic];
    xu[ic] = mciin->domain[2*ic+1];
    if(verboseflag >= 2) {
      printf("%% domain of integration for the %d-th variable %f %f\n",ic+1, xl[ic],xu[ic]);
    }
  } 

  s = gsl_monte_vegas_alloc (mciin->dimension);

  gsl_monte_vegas_integrate (&G, xl, xu, mciin->dimension, warmupcalls, r, s,
			     &res, &err);
  if(verboseflag >= 2) {
    printf ("%% vegas warm-up: result = %g; sigma = %g\n", res, err);
  }

  ic = 0;
  do
    {
      gsl_monte_vegas_integrate (&G, xl, xu, mciin->dimension, calls/5, r, s,
				 &res, &err);
      if(verboseflag >= 2) {
	printf ("%% %d-th stage:  ",ic+1 );
	printf ("%% result = % .6f sigma = % .6f "
		"chisq/dof = %.1f\n", res, err, s->chisq);
      }
      ic++;
    }
  while ((fabs (s->chisq - 1.0) > 0.5) && (ic < maxVegasStages));

  if(ic==maxVegasStages) {
    fprintf(stderr,"%% Maximum stages reached in the vegas method.\n");
    fprintf(stderr,"%% The chi-square value is %10.6g\n",s->chisq);
  }

  if(verboseflag >= 2) {
    printf ("%% vegas final: result = %22.12g; sigma = %22.12g\n", res, err);
  }

  gsl_monte_vegas_free (s);

  gsl_rng_free(r);
  free(xl);
  free(xu);

  mciout->result = res;
  mciout->sigma = err;

  return;
} /* void MonteCarloIntegrate() */




/*----------------------------------------------------------------------
 * FUNCTION REAL8 marginalizedNCPDF()
 *  
 * + marginalizing Non-Central-Chi-Square pdf. 
 *  
 *
 *----------------------------------------------------------------------*/
REAL8 marginalizedNCPDF(REAL8 strain, BayesFstatParams *bfparams)
{
  REAL8 dinv,dinvsigma,normFfactor;
  REAL8 norm, normsigma;
  REAL8 res, ressigma;
  MonteCarloIntParams mciparams;
  MonteCarloIntIn mciin;
  MonteCarloIntOut mciout;


  mciparams.verboseflag = lalDebugLevel;
  mciparams.method = "vegas";
  mciparams.warmupcalls = 10000;
  mciparams.maincalls = 500000;
  mciparams.maxVegasStages = 100;

  mciin.params = (void *) bfparams;

  mciin.domain = (REAL8 *) calloc(10,sizeof(REAL8));
  mciin.dimension = 4;
  /* alpha  */
  mciin.domain[0] = 0.0;        mciin.domain[1] = LAL_TWOPI;
  /* delta */
  mciin.domain[2] = -LAL_PI_2;  mciin.domain[3] = LAL_PI_2;
  /* cosiota */
  mciin.domain[4] = 0.0;        mciin.domain[5] = 1.0;
  /* psi */
  mciin.domain[6] = - LAL_PI_4; mciin.domain[7] = LAL_PI_4;


#if 1
  if( bfparams->norm < 0.0 ) {
    computeNorm(&dinv,&dinvsigma,&mciin,&mciparams);
    normFfactor = NormFunction(bfparams->Fstat);		 
    norm = dinv * normFfactor;
    normsigma = dinvsigma * normFfactor;
    bfparams->norm = norm;
    bfparams->normsigma = normsigma;
  }
  norm = bfparams->norm;
  normsigma = bfparams->normsigma;
#endif

#if 0
    computeNorm(&dinv,&dinvsigma,&mciin,&mciparams);
    normFfactor = NormFunction(bfparams->Fstat);		 
    norm = dinv * normFfactor;
    normsigma = dinvsigma * normFfactor;
    bfparams->norm = norm;
    bfparams->normsigma = normsigma;
#endif

  res = 0.0;
  ressigma = pow(normsigma/norm,2);

  mciin.dimension = 5; 
  mciin.integrand = ncpdfWithIntMeasure;  
  mciin.domain[8] = 0.0;
  mciin.domain[9] = strain;

  MonteCarloIntegrate(&mciout, &mciin, &mciparams);
  ressigma += pow(mciout.sigma/mciout.result, 2);
  res = mciout.result;

  if(lalDebugLevel>=2) {
    printf("%% Fstat , normalized strain ,  confidence ,  err(confidence) \n");
    printf("%% %10.6g %10.6g  %10.6g %10.6g\n",bfparams->Fstat, strain, res/norm, res/norm*sqrt(ressigma));
  }
  bfparams->perr = res/norm*sqrt(ressigma);

  free(mciin.domain);

  return res/norm;
} /* REAL8 marginalizedNCPDF() */




/*----------------------------------------------------------------------
 * FUNCTION void computeNorm()
 * 
 * + Integrate 1/D over alpha, delta, cosiota, psi.
 * + Result is the norm apart from the F-stat depending part.
 *
 *----------------------------------------------------------------------*/
void computeNorm(REAL8 *norm, REAL8 *normsigma, MonteCarloIntIn *mciin, MonteCarloIntParams *mciparams)
{
  MonteCarloIntOut mciout;

  mciout.result = 0.0;
  mciout.sigma = 0.0;
  /* 4-dimensional integration. alpha, delta, cosiota, psi.  */
  mciin->dimension = 4;

#ifdef USE_S2ALLSKY /* to speed-up... */
  /* These values are for S2 10 hours all-sky. */
  if( !strcmp(uvar_IFO,"LLO") ) {
    mciout.result = 61.78979323250208;
    mciout.sigma = 0.0126945;
  } else if( !strcmp(uvar_IFO,"LHO") ) {
    mciout.result = 63.07430522601882;
    mciout.sigma = 0.0163082;
  } else {
    /* We need this norm only once for a given timestamps and an IFO. */
    mciin->integrand = DInverseWithIntMeasure;
    MonteCarloIntegrate(&mciout, mciin, mciparams);
  } 
#else  /* in general... */
  /* We need this norm only once for a given timestamps and an IFO. */
  mciin->integrand = DInverseWithIntMeasure;
  MonteCarloIntegrate(&mciout, mciin, mciparams);
#endif

  *normsigma = (mciout.sigma);
  *norm = (mciout.result);

  if(lalDebugLevel >= 3) {
    printf("%% \\integrate D^{-1} = %22.16g, sigma = %g\n", mciout.result,mciout.sigma);
  }

  return;
} /* void computeNorm() */


/*----------------------------------------------------------------------
 * FUNCTION REAL8 ncpdfWithInteMeasure()
 * 
 * + Non-central chi-square distribution with 4 degrees of freedom multiplied 
 * by the integral measure cos(delta) for the delta integration.
 *
 *
 *----------------------------------------------------------------------*/

REAL8 
ncpdfWithIntMeasure(REAL8 *k, size_t dim, void *params)
{
  /* 
     k[0]=alpha;
     k[1]=delta;
     k[2]=mu;
     k[3]=psi;
     k[4]=H = h_0 \sin\theta \sin\zeta \sqrt{T/S_h} / 6 / \sqrt{2F}; 
     The magic factor 6 makes the 95 % UL H_{0.95} around 1. 
  */
  REAL8 ncpdfvalue = 0.0;
  REAL8 k1[4];
  REAL8 F;
  REAL8 dvalH;
  INT4 dim1;
  INT4 ic;
  REAL8 magic = 6.0;
  BayesFstatParams *bfparams;

  if(dim != 5) {
    fprintf(stderr,"Error: the dimension should be 5 for ncpdfWithIntMeasure().\n");
  }
  dim1 = 4;
  for(ic=0;ic<dim1;ic++) {
    k1[ic]=k[ic];
  }

  bfparams = params;

  F = bfparams->Fstat;
  dvalH = magic * k[4] * DFunction(k1, dim1, params);

  ncpdfvalue = 2.0*F*dvalH - F - F*pow(dvalH,2);
  ncpdfvalue = exp(ncpdfvalue);
  ncpdfvalue *= sqrt(2*F);
  ncpdfvalue *= gsl_sf_bessel_In_scaled(3,2.0*F*dvalH);
  ncpdfvalue *= 1.0/pow(dvalH,3);
  ncpdfvalue *= cos(k[1]);

  ncpdfvalue *= magic;

  return ncpdfvalue;
} /* REAL8 ncpdfWithInteMeasure() */




/*----------------------------------------------------------------------
 * FUNCTION void InitBayesFstatParams(). 
 *
 *
 *----------------------------------------------------------------------*/
void 
InitBayesFstatParams(LALStatus *status, BayesFstatParams *bfparams, INT4 argc, CHAR *argv[])
{

  ConfigVariables *cfg = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT( bfparams != NULL, status, BAYESFSTATC_ENULL, BAYESFSTATC_MSGENULL);

  TRY( InitUserVars(status->statusPtr), status);
  TRY( LALUserVarReadAllInput(status->statusPtr,argc,argv), status);

  if (uvar_help)	/* if help was requested, we're done here */
    exit(BAYESFSTAT_EXIT_OK);

  cfg = (ConfigVariables *) LALMalloc(sizeof(ConfigVariables));

  TRY( InitConfigVariables(status->statusPtr, cfg), status);
  TRY( InitAMParams(status->statusPtr, bfparams, cfg), status);

  /* Store the user-specified F stat to our parameters structure */
  bfparams->Fstat = cfg->Fstat;
  bfparams->confidence = cfg->confidence;

  bfparams->perr = - 1.0;
  bfparams->norm = - 1.0;
  bfparams->normsigma = - 1.0;

  /* Look-Up-Table Jaranowski Krolak Schutz A,B,C coefficients. */
  {
  REAL8 gridsize;
  UINT4 Nalpha, Ndelta;

  bfparams->ABCs = NULL;
  bfparams->ABCs = (ABCcoeffs *) LALMalloc(sizeof(ABCcoeffs));
  bfparams->ABCs->gridsize = 0.01;
  gridsize = bfparams->ABCs->gridsize;
  Nalpha = (UINT4) (floor(LAL_TWOPI/gridsize) + 3);
  Ndelta = (UINT4) (floor(LAL_PI/gridsize) + 3);
  AllocateABCcoeffs(status->statusPtr, bfparams->ABCs,Ndelta,Nalpha);
  BEGINFAIL(status) {
    TRY( FreeMem(status->statusPtr, cfg, bfparams), status);
  } ENDFAIL(status);
  ABCLookUpTable(status->statusPtr, bfparams->ABCs, bfparams, gridsize, Ndelta, Nalpha);
  BEGINFAIL(status) {
    TRY( FreeMem(status->statusPtr, cfg, bfparams), status);
  } ENDFAIL(status);
  }


  FreeConfigVariables(status->statusPtr, cfg);
  BEGINFAIL(status) {TRY( FreeMem(status->statusPtr, cfg,bfparams),status);} ENDFAIL(status);

  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* void InitBayesFstatParams().  */



/*----------------------------------------------------------------------
 * FUNCTION void InitUserVars(). Mostly from ComputeFStatistic.
 *
 *
 *----------------------------------------------------------------------*/
void InitUserVars(LALStatus *status)
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* set-up: IFO, timestampsfilename, ephemeris stuffs, tsft. */

  uvar_ephemYear = (CHAR*)LALCalloc (1, strlen(EPHEM_YEARS)+1);
  if(uvar_ephemYear == NULL) {
    ABORT (status, BAYESFSTATC_EMEM, BAYESFSTATC_MSGEMEM);
  }
  strcpy (uvar_ephemYear, EPHEM_YEARS);
#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
  uvar_ephemDir = (CHAR*)LALCalloc (1, strlen(DEFAULT_EPHEMDIR)+1);
  if(uvar_ephemDir == NULL) {
    LALFree(uvar_ephemYear);
    ABORT (status, BAYESFSTATC_EMEM, BAYESFSTATC_MSGEMEM);
  }
  strcpy (uvar_ephemDir, DEFAULT_EPHEMDIR);

  uvar_useLUT = 1;
  uvar_Fstat = 10.0;
  uvar_Tsft = 1800.0;
  uvar_confidence = 0.95; /* targeting confidence level, given F stat */

  /* register all our user-variables */
  LALregBOOLUserVar(status, 	help, 		'h', UVAR_HELP,     "Print this message"); 
  LALregBOOLUserVar(status, 	useLUT, 		'L', UVAR_OPTIONAL,     "use Look-Up-Table for JKS ABC coefficients."); 
  LALregSTRINGUserVar(status, 	IFO, 		'I', UVAR_REQUIRED, "Detector: GEO(0),LLO(1),LHO(2)");
  LALregSTRINGUserVar(status, 	timestampsfile, 	'T', UVAR_REQUIRED, "timestamps file name");
  LALregREALUserVar(status, 	Tsft, 		't', UVAR_OPTIONAL, "time-length of SFT in seconds");
  LALregREALUserVar(status, 	Fstat, 		'F', UVAR_OPTIONAL, "F statistic");
  LALregREALUserVar(status, 	confidence, 		'C', UVAR_OPTIONAL, "targeting confidence level, for a given F stat");
  LALregSTRINGUserVar(status,	ephemDir, 	'E', UVAR_OPTIONAL, "Directory where Ephemeris files are located");
  LALregSTRINGUserVar(status,	ephemYear, 	'y', UVAR_OPTIONAL, "Year (or range of years) of ephemeris files to be used");
  LALregSTRINGUserVar(status,	outputLabel,	'o', UVAR_OPTIONAL, "Label to be appended to all output file-names");


  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* void InitUserVars(). */





/*----------------------------------------------------------------------
 * FUNCTION void InitConfigVariables(). Mostly from ComputeFStatistic.
 *
 * + Input: command line arguments 
 * + Output: ConfigureVariables *cfg
 * + Initialize a ConfigVariables *cfg defined at the header of this file. 
 *
 *----------------------------------------------------------------------*/

void 
InitConfigVariables(LALStatus *status, ConfigVariables *cfg)
{


  LIGOTimeGPS starttime;
  FILE *fp;
  INT4 nout;


  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  ASSERT( cfg != NULL, status, BAYESFSTATC_ENULL, BAYESFSTATC_MSGENULL);

  /*----------------------------------------------------------------------
   * timestamps file
   */
  
  if( !LALUserVarWasSet (&uvar_timestampsfile) ) {
    fprintf(stderr,"\ntimestamps file name should be specified.\n\n");
    ABORT (status, BAYESFSTATC_EINPUT, BAYESFSTATC_MSGEINPUT);
  }
   cfg->timestampsfile = (CHAR *) LALMalloc(strlen(uvar_timestampsfile)*sizeof(CHAR)+1);
  if( cfg->timestampsfile == NULL ) {
    fprintf(stderr,"\nCannot allocate memory for timestampsfile. \n\n");
    ABORT (status, BAYESFSTATC_EMEM, BAYESFSTATC_MSGEMEM);
  }
  strcpy(cfg->timestampsfile, uvar_timestampsfile);

  fp=fopen(cfg->timestampsfile,"r");
  if (fp==NULL) {
    fprintf(stderr,"Unable to find file %s\n",cfg->timestampsfile);
    LALFree(cfg->timestampsfile);
    ABORT (status, BAYESFSTATC_ESYS, BAYESFSTATC_MSGESYS);
  }
  nout=fscanf(fp,"%d  %d\n", &(starttime.gpsSeconds), &(starttime.gpsNanoSeconds));
  if ( nout !=2 ) {
    fprintf(stderr,"Unable to read datum # %d\n",1);
    fprintf(stderr,"from file %s\n",cfg->timestampsfile);
    LALFree(cfg->timestampsfile);
    ABORT (status, BAYESFSTATC_ESYS, BAYESFSTATC_MSGESYS);
  }  
  fclose(fp);


  /*----------------------------------------------------------------------
   * initilize Tsft 
   */

  if( !LALUserVarWasSet (&uvar_Tsft) ) {
    cfg->tsft = 1800.0;
    /*
    fprintf(stderr,"\n time length of SFT should be specified.\n\n");
    LALFree(cfg->timestampsfile);
    ABORT (status, BAYESFSTATC_EINPUT, BAYESFSTATC_MSGEINPUT);
    */
  } else {
    cfg->tsft = uvar_Tsft;
  }


  /*----------------------------------------------------------------------
   * initilize Tsft 
   */

  if( LALUserVarWasSet (&uvar_Fstat) ) {
    cfg->Fstat = uvar_Fstat;
  } else {
    cfg->Fstat = 10.0;
  }


  /*----------------------------------------------------------------------
   * initilize Tsft 
   */

  if( LALUserVarWasSet (&uvar_confidence) ) {
    cfg->confidence = uvar_confidence;
  } else {
    cfg->confidence = 0.95;
  }


  /*----------------------------------------------------------------------
   * initialize detector 
   */
  if( !LALUserVarWasSet (&uvar_IFO) ) {
    fprintf(stderr,"\n IFO should be specified.\n\n");
    LALFree(cfg->timestampsfile);
    ABORT (status, BAYESFSTATC_EINPUT, BAYESFSTATC_MSGEINPUT);
  }

  if ( !strcmp (uvar_IFO, "GEO") || !strcmp (uvar_IFO, "0") ) 
    cfg->Detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  else if ( !strcmp (uvar_IFO, "LLO") || ! strcmp (uvar_IFO, "1") ) 
    cfg->Detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  else if ( !strcmp (uvar_IFO, "LHO") || !strcmp (uvar_IFO, "2") )
    cfg->Detector = lalCachedDetectors[LALDetectorIndexLHODIFF];
  else
    {
      fprintf(stderr,"\nUnknown detector. Currently allowed are 'GEO', 'LLO', 'LHO', '0'-'2'\n\n");
      LALFree(cfg->timestampsfile);
      ABORT (status, BAYESFSTATC_EINPUT, BAYESFSTATC_MSGEINPUT);
    }




  /* ----------------------------------------------------------------------*/
  /*
   * ephemeris 
   */

  if (uvar_ephemYear == NULL)
    {
      fprintf(stderr,"\nNo ephemeris year specified (option 'ephemYear')\n\n");
      LALFree(cfg->timestampsfile);
      ABORT (status, BAYESFSTATC_EINPUT, BAYESFSTATC_MSGEINPUT);
    }     

#define EPHEM_EXT ".dat"
  {INT4 filenamelenE, filenamelenS;
  filenamelenE = strlen("/earth") + strlen(EPHEM_EXT) + strlen(uvar_ephemYear) + 1;
  filenamelenS = strlen("/sun") + strlen(EPHEM_EXT) + strlen(uvar_ephemYear) + 1;

  if (LALUserVarWasSet (&uvar_ephemDir) ) {
    filenamelenE += strlen(uvar_ephemDir);
    filenamelenS += strlen(uvar_ephemDir);
  }
  cfg->EphemEarth = (CHAR *) LALCalloc(1,filenamelenE);
  cfg->EphemSun = (CHAR *) LALCalloc(1,filenamelenS);
  if( cfg->EphemEarth == NULL || 
      cfg->EphemSun == NULL ) {
    LALFree(cfg->timestampsfile);
    fprintf(stderr,"\nCannot allocate memory for ephemeris files. \n\n");
    ABORT (status, BAYESFSTATC_EMEM, BAYESFSTATC_MSGEMEM);
  }

  if (LALUserVarWasSet (&uvar_ephemDir) ) 
    {
      sprintf(cfg->EphemEarth,"%s/earth%s" EPHEM_EXT, uvar_ephemDir, uvar_ephemYear);
      sprintf(cfg->EphemSun, "%s/sun%s" EPHEM_EXT, uvar_ephemDir, uvar_ephemYear);
    } 
  else 
    {
      sprintf(cfg->EphemEarth, "earth%s" EPHEM_EXT, uvar_ephemYear);
      sprintf(cfg->EphemSun, "sun%s" EPHEM_EXT,  uvar_ephemYear);
    }
  }
  

  if( (fp=fopen(cfg->EphemEarth,"r")) == NULL ) {
    fprintf(stderr,"\nCann not open the earth ephemeris file %s.\n\n",cfg->EphemEarth);
    LALFree(cfg->timestampsfile);
    ABORT (status, BAYESFSTATC_EINPUT, BAYESFSTATC_MSGEINPUT);
  } else {
    fclose(fp);
  }
  if( (fp=fopen(cfg->EphemSun,"r")) == NULL ) {
    fprintf(stderr,"\nCann not open the sun ephemeris file %s.\n\n",cfg->EphemSun);
    LALFree(cfg->timestampsfile);
    ABORT (status, BAYESFSTATC_EINPUT, BAYESFSTATC_MSGEINPUT);
  } else {
    fclose(fp);
  }


  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* void InitConfigVariables().  */


/*----------------------------------------------------------------------
 *FUNCTION void filelength().
 *
 * + Count a number of rows in a file given by a user as CHAR *filename.
 * + Count it up to MaxFileLen.
 * + Return the number of rows in the file.
 * + Return -1 if failed.
 *
 *----------------------------------------------------------------------*/


void  
filelength(LALStatus *status, INT8 *length, CHAR *filename)
{
  INT8 len = 0, counter=0;
  FILE *fp = NULL;
  INT8 MaxFileLen = 65536;
  INT4 buffsize = 256;
  CHAR line[256];
  CHAR *ptr;

  INITSTATUS(status);

  fp=fopen(filename,"r");
  if (fp==NULL) {
    fprintf(stderr,"Unable to find file %s\n",filename);
    ABORT (status, BAYESFSTATC_EINPUT, BAYESFSTATC_MSGEINPUT); 
  }
  
  while(counter < MaxFileLen ) 
    {
      ptr = fgets(line, buffsize, fp);
      if( ptr == NULL ) { 
	break;
      }
      counter++;
     }
  fclose(fp);

  len = counter;
  if( len == MaxFileLen ) {
    fprintf(stderr, "Maximum number of rows reached during reading the file %s\n", filename);
  }

  *length = len;

  RETURN (status);

} /* void filelength() */


/*----------------------------------------------------------------------
 * FUNCTION void InitAMParams(). Mostly from ComputeFStatistic.
 *
 * + Input: ConfigureVariables *cfg defined at the header of this file.
 * + Output: BayesFstatParams *bfparams defined at the header of this file.
 * + *bfparams contains all the necessary inputs for LALComputeAM().
 * 
 *----------------------------------------------------------------------*/


void
InitAMParams(LALStatus *status, BayesFstatParams *bfparams, ConfigVariables *cfg)
{


  /* ASK ME!!!! IS THIS CORRECT? CAN THIS DISAPPEAR AFTER EXITING THIS FUNCTION? */
  EarthState earth;

  LIGOTimeGPS *timestamps = NULL;
  REAL8 Alpha = 0.0, Delta = 0.0;

  INT8 ntimestamps = 0;
  INT4 ic;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  ASSERT( cfg != NULL, status, BAYESFSTATC_ENULL, BAYESFSTATC_MSGENULL);

  /* Detector location: MAKE INTO INPUT!!!!! */
  bfparams->baryinput.site.location[0]=cfg->Detector.location[0]/LAL_C_SI;
  bfparams->baryinput.site.location[1]=cfg->Detector.location[1]/LAL_C_SI;
  bfparams->baryinput.site.location[2]=cfg->Detector.location[2]/LAL_C_SI;
  bfparams->baryinput.dInv=0.e0;
  bfparams->baryinput.alpha=Alpha;
  bfparams->baryinput.delta=Delta;

/* amParams structure to compute a(t) and b(t) */
/* Allocate space for amParams stucture */
/* Here, amParams->das is the Detector and Source info */
  bfparams->amParams = (AMCoeffsParams *)LALMalloc(sizeof(AMCoeffsParams));
  bfparams->amParams->das = (LALDetAndSource *)LALMalloc(sizeof(LALDetAndSource));
  bfparams->amParams->das->pSource = (LALSource *)LALMalloc(sizeof(LALSource));
  bfparams->amParams->das->pDetector = (LALDetector *)LALMalloc(sizeof(LALDetector));
  bfparams->amParams->edat = (EphemerisData*)LALCalloc(1, sizeof(EphemerisData));
  if( bfparams->amParams == NULL || 
      bfparams->amParams->das == NULL || 
      bfparams->amParams->das->pSource == NULL || 
      bfparams->amParams->das->pDetector == NULL || 
      bfparams->amParams->edat == NULL
      ) {
    TRY( FreeMem(status->statusPtr, cfg,bfparams),status);
    ABORT (status, BAYESFSTATC_EMEM, BAYESFSTATC_MSGEMEM); 
  }


  /* ----------------------------------------------------------------------*/
  /*
   * initialize Ephemeris-data 
   */
  {
    bfparams->amParams->edat->ephiles.earthEphemeris = (CHAR *) LALMalloc(strlen(cfg->EphemEarth)+1);
    bfparams->amParams->edat->ephiles.sunEphemeris = (CHAR *) LALMalloc(strlen(cfg->EphemEarth)+1);
    if( bfparams->amParams->edat->ephiles.earthEphemeris == NULL || 
	bfparams->amParams->edat->ephiles.sunEphemeris == NULL 
	) {
      TRY( FreeMem(status->statusPtr, cfg,bfparams),status);
      ABORT (status, BAYESFSTATC_EMEM, BAYESFSTATC_MSGEMEM); 
    }
    strcpy(bfparams->amParams->edat->ephiles.earthEphemeris, cfg->EphemEarth);
    strcpy(bfparams->amParams->edat->ephiles.sunEphemeris, cfg->EphemSun);

    LALInitBarycenter(status->statusPtr, bfparams->amParams->edat);               
    BEGINFAIL(status) {
      TRY( FreeMem(status->statusPtr, cfg, bfparams),status);
    } ENDFAIL(status);

  } /* end: init ephemeris data */


/* Fill up AMCoeffsParams structure */
  bfparams->amParams->earth = &earth; 
  *(bfparams->amParams->das->pDetector) = cfg->Detector; 
  bfparams->amParams->das->pSource->orientation = 0.0;
  bfparams->amParams->das->pSource->equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  bfparams->amParams->polAngle = bfparams->amParams->das->pSource->orientation ; /* These two have to be the same!!!!!!!!!*/
  bfparams->amParams->das->pSource->equatorialCoords.latitude = Delta;
  bfparams->amParams->das->pSource->equatorialCoords.longitude = Alpha;
  bfparams->amParams->baryinput = &(bfparams->baryinput);

  filelength(status->statusPtr, &ntimestamps, cfg->timestampsfile);
  BEGINFAIL(status) { TRY( FreeMem(status->statusPtr, cfg,bfparams),status); } ENDFAIL(status);

  timestamps=(LIGOTimeGPS *)LALMalloc(ntimestamps*sizeof(LIGOTimeGPS)); 
  if( timestamps == NULL ) {
    TRY( FreeMem(status->statusPtr, cfg,bfparams),status);
    ABORT (status, BAYESFSTATC_EMEM, BAYESFSTATC_MSGEMEM); 
  }

  ReadTimeStamps(status->statusPtr, timestamps, ntimestamps, cfg->timestampsfile);
  BEGINFAIL(status) {
    fprintf(stderr,"timestamps file reading error %s\n",cfg->timestampsfile);
    TRY( FreeMem(status->statusPtr, cfg,bfparams),status); 
  } ENDFAIL(status);


 /* Mid point of each SFT */
  bfparams->midTS = (LIGOTimeGPS *)LALCalloc(ntimestamps, sizeof(LIGOTimeGPS));
  if( bfparams->midTS == NULL ) {
    TRY( FreeMem(status->statusPtr, cfg,bfparams),status);
    ABORT (status, BAYESFSTATC_EMEM, BAYESFSTATC_MSGEMEM); 
  }

  for(ic=0; ic< ntimestamps; ic++)
    { 
      /* FIXME:  loss of precision; consider
      bfparams->midTS[ic] = timestamps[ic];
      XLALGPSAdd(&bfparams->midTS[ic], 0.5*(cfg->tsft));
      */
      REAL8 teemp=0.0;
      teemp = XLALGPSGetREAL8(&(timestamps[ic]));
      teemp += 0.5*(cfg->tsft);
      XLALGPSSetREAL8(&(bfparams->midTS[ic]), teemp);
    }

  LALFree(timestamps);

  /* Allocate space for AMCoeffs */
  bfparams->amc.a = NULL;
  bfparams->amc.b = NULL;
  LALSCreateVector(status->statusPtr, &(bfparams->amc.a), (UINT4) ntimestamps);
  BEGINFAIL(status) {TRY( FreeMem(status->statusPtr, cfg,bfparams),status); } ENDFAIL(status);
  LALSCreateVector(status->statusPtr, &(bfparams->amc.b), (UINT4) ntimestamps);
  BEGINFAIL(status) {TRY( FreeMem(status->statusPtr, cfg,bfparams),status);} ENDFAIL(status);


  DETATCHSTATUSPTR (status);
  RETURN (status);
}  /* void InitAMParams(). */
 


/*----------------------------------------------------------------------
 * FUNCTION int ReadTimeStamps(). Mostly from SemiAnalyticF.
 *   
 * + Read timestamps file given by a user as CHAR *timestampsfile.
 * + Read timestamps file up to a number of rows give by a user as long ntimestamps.
 * + Store the timestamps into LIGOTimeGPS *timestamps, memory of which must be allocated .
 *   outside of this function.
 * + Return 0 if succeed, 1 otherwise. 
 *----------------------------------------------------------------------*/


void
ReadTimeStamps(LALStatus *status, LIGOTimeGPS *timestamps, INT8 ntimestamps, CHAR *timestampsfile) 
{
  FILE *fp = NULL;
  INT4 ic = 0, nout = -1;

  INITSTATUS(status);

  fp=fopen(timestampsfile,"r");
  if (fp==NULL) {
    fprintf(stderr,"Unable to find file %s\n",timestampsfile);
    ABORT (status, BAYESFSTATC_EINPUT, BAYESFSTATC_MSGEINPUT); 
  }

  for (ic=0;ic< ntimestamps;ic++){
    nout=fscanf(fp,"%d  %d\n", &timestamps[ic].gpsSeconds, &timestamps[ic].gpsNanoSeconds);
    if ( nout !=2 ) {
      fprintf(stderr,"Unable to read datum # %d\n",ic);
      fprintf(stderr,"from file %s\n",timestampsfile);
      ABORT (status, BAYESFSTATC_EINPUT, BAYESFSTATC_MSGEINPUT); 
    } 
  } 
  
  fclose(fp);

  RETURN (status);
} /* int ReadTimeStamps() */



/*----------------------------------------------------------------------
 * FUNCTION void FreeMem().
 *
 *
 *
 *
 *
 *----------------------------------------------------------------------*/


void FreeMem(LALStatus *status, ConfigVariables *cfg, BayesFstatParams *bfparams)
{

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  ASSERT( cfg != NULL, status, BAYESFSTATC_ENULL, BAYESFSTATC_MSGENULL);
  ASSERT( bfparams != NULL, status, BAYESFSTATC_ENULL, BAYESFSTATC_MSGENULL);

  FreeBayesFstatParams(status->statusPtr, bfparams);
  BEGINFAIL(status) {
      TRY( FreeConfigVariables(status->statusPtr, cfg), status);
  } ENDFAIL(status);
  TRY( FreeConfigVariables(status->statusPtr, cfg), status);

  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* void FreeMem(). */




/*----------------------------------------------------------------------
 * FUNCTION void FreeConfigVariables().
 *
 *
 *
 *
 *
 *----------------------------------------------------------------------*/


void FreeConfigVariables(LALStatus *status, ConfigVariables *cfg)
{

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  ASSERT( cfg != NULL, status, BAYESFSTATC_ENULL, BAYESFSTATC_MSGENULL);

  if(cfg->EphemSun != NULL) {
    LALFree(cfg->EphemSun);
  }

  if(cfg->EphemEarth != NULL) {
    LALFree(cfg->EphemEarth);
  }

  if(cfg->timestampsfile != NULL) {
    LALFree(cfg->timestampsfile);
  }
  if(cfg != NULL) 
    LALFree(cfg);

  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* void FreeConfigVariables(). */


/*----------------------------------------------------------------------
 * FUNCTION void FreeBayesFstatParams().
 *
 *
 *
 *
 *
 *----------------------------------------------------------------------*/


void FreeBayesFstatParams(LALStatus *status, BayesFstatParams *bfparams)
{

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  ASSERT( bfparams != NULL, status, BAYESFSTATC_ENULL, BAYESFSTATC_MSGENULL);


  if(bfparams->ABCs != NULL ) {
    TRY( FreeABCcoeffs(status->statusPtr, bfparams->ABCs), status);
  }

  if(bfparams->amc.a != NULL ) {
    TRY( LALSDestroyVector(status->statusPtr, &(bfparams->amc.a)), status);
  }
  if(bfparams->amc.b != NULL ) {
    TRY( LALSDestroyVector(status->statusPtr, &(bfparams->amc.b)), status);
  }

  if(bfparams->midTS != NULL ) {
    LALFree(bfparams->midTS);
  }


  if(bfparams->amParams->edat != NULL) {
    LALFree(bfparams->amParams->edat->ephemE);
    LALFree(bfparams->amParams->edat->ephemS);
    LALFree(bfparams->amParams->edat->ephiles.earthEphemeris);
    LALFree(bfparams->amParams->edat->ephiles.sunEphemeris);
    LALFree(bfparams->amParams->edat);
  }

  if(bfparams->amParams != NULL ) {
    LALFree(bfparams->amParams->das->pDetector);
    LALFree(bfparams->amParams->das->pSource);
    LALFree(bfparams->amParams->das);
    LALFree(bfparams->amParams);
  }
  
  if(bfparams != NULL ) {
    LALFree(bfparams);
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* void FreeBayesFstatParams(). */








/*----------------------------------------------------------------------
 * FUNCTION REAL8 NormFunction().
 *
 *
 *
 *
 *
 *----------------------------------------------------------------------*/

REAL8 
NormFunction(REAL8 F) 		    
{ 
  REAL8 result = 0.0;
  result = F * ( 1.0 + 2.0*F ) * gsl_sf_bessel_I0_scaled(F/2.0);
  result -= ( 4.0 + 3.0*F + 2.0*pow(F,2) ) * gsl_sf_bessel_I1_scaled(F/2.0);
  result *= sqrt(LAL_TWOPI)/15 * F;
  return result;
} /* REAL8 NormFunction(). */


/*----------------------------------------------------------------------
 * FUNCTION REAL8 DFunction().
 *
 *
 *
 *
 *
 *----------------------------------------------------------------------*/

REAL8 
DFunction(REAL8 *k,
	  size_t dim,
	  void *params)
{
  /* 
     k[0]=alpha;
     k[1]=delta;
     k[2]=mu;
     k[3]=psi;
   */
  /*  LALStatus *status = &global_status; */
  LALStatus status = global_status;
  BayesFstatParams *bfparams;


  REAL8 Dvalue = 0.0;
  REAL8 A, B, C;
  INT4 useLUT;

  bfparams = (BayesFstatParams *) params;

  if(dim != 4) {
    fprintf(stderr,"Error: the dimension should be 4 for DFunction().\n");
  }

  useLUT = uvar_useLUT;
  if( useLUT == 1) {
    REAL8 res = 0.0;
    Interpolate2DSphere(&res, bfparams->ABCs->alpha, bfparams->ABCs->delta, bfparams->ABCs->As, k[0], k[1],bfparams->ABCs->gridsize);
    A = res;

    Interpolate2DSphere(&res, bfparams->ABCs->alpha, bfparams->ABCs->delta, bfparams->ABCs->Bs, k[0], k[1],bfparams->ABCs->gridsize);
    B = res;

    Interpolate2DSphere(&res, bfparams->ABCs->alpha, bfparams->ABCs->delta, bfparams->ABCs->Cs, k[0], k[1],bfparams->ABCs->gridsize);
    C = res;

  } else {

    bfparams->baryinput.alpha = k[0];
    bfparams->baryinput.delta = k[1];
    bfparams->amParams->das->pSource->equatorialCoords.latitude = k[1];
    bfparams->amParams->das->pSource->equatorialCoords.longitude = k[0];
    
    LALComputeAM(&status, &(bfparams->amc), bfparams->midTS, bfparams->amParams); 
    
    A = bfparams->amc.A;
    B = bfparams->amc.B;
    C = bfparams->amc.C;
  }

  Dvalue = ( A + B ) * ( 1.0 + 6.0*pow(k[2],2) + pow(k[2],4) );
  Dvalue += ( A - B ) * pow( 1.0 - pow(k[2],2), 2 ) * cos(4.0*k[3]);
  Dvalue += 2.0*C * pow( 1.0 - pow(k[2],2), 2 ) * sin(4.0*k[3]);
  Dvalue = sqrt(Dvalue);
  Dvalue *= 1.0/4.0;

  return Dvalue;
} /* REAL8 DFunction(). */



/*----------------------------------------------------------------------
 * FUNCTION void testInterPolate2DSphere() 
 *
 *
 *
 *
 *
 *----------------------------------------------------------------------*/
void testInterPolate2DSphere(LALStatus *status, BayesFstatParams *bfparams)
{
  REAL8 alpha, delta;
  UINT4 ic,N;


  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT( bfparams != NULL, status, BAYESFSTATC_ENULL, BAYESFSTATC_MSGENULL);

  N=100000;
  for(ic=0;ic<N;ic++) {
    alpha = LAL_TWOPI*rand()/(RAND_MAX+1.0);
    delta = LAL_PI*rand()/(RAND_MAX+1.0) - LAL_PI_2;
    TRY( testInterPolate2DSphereCore(status->statusPtr, bfparams, alpha, delta), status);
  }
  /*  alpha = 3.98806; delta = 1.56974 */;
  alpha = 3.15986 ; delta = 1.56922;
  TRY( testInterPolate2DSphereCore(status->statusPtr, bfparams, alpha, delta), status);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}  /* void testInterPolate2DSphere() */


/*----------------------------------------------------------------------
 * FUNCTION void testInterPolate2DSphereCore() 
 *
 *
 *
 *
 *
 *----------------------------------------------------------------------*/

void testInterPolate2DSphereCore(LALStatus *status, BayesFstatParams *bfparams, REAL8 alpha, REAL8 delta) 
{
  REAL8 A, B, C;
  REAL8 res = 0.0;



  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  ASSERT( bfparams != NULL, status, BAYESFSTATC_ENULL, BAYESFSTATC_MSGENULL);

  Interpolate2DSphere(&res, bfparams->ABCs->alpha, bfparams->ABCs->delta, bfparams->ABCs->As, alpha, delta,bfparams->ABCs->gridsize);
  A = res;
  Interpolate2DSphere(&res, bfparams->ABCs->alpha, bfparams->ABCs->delta, bfparams->ABCs->Bs, alpha, delta,bfparams->ABCs->gridsize);
  B = res;
  Interpolate2DSphere(&res, bfparams->ABCs->alpha, bfparams->ABCs->delta, bfparams->ABCs->Cs, alpha, delta,bfparams->ABCs->gridsize);
  C = res;

  bfparams->baryinput.alpha = alpha;
  bfparams->baryinput.delta = delta;
  bfparams->amParams->das->pSource->equatorialCoords.latitude = delta;
  bfparams->amParams->das->pSource->equatorialCoords.longitude = alpha;
  
  TRY( LALComputeAM(status->statusPtr, &(bfparams->amc), bfparams->midTS, bfparams->amParams), status); 

#if 0
  printf("Interpolate results: %g %g %g\n",A ,B ,C);
  printf("LALComputeAM results: %g %g %g\n",bfparams->amc.A,bfparams->amc.B,bfparams->amc.C);
  fflush(stdout);
#endif 

  printf("%g \t %g \t ",alpha ,delta);
  printf("%g \t %g \t %g \t ",A ,B ,C);
  printf("%g \t %g \t %g\n",bfparams->amc.A,bfparams->amc.B,bfparams->amc.C);
  fflush(stdout);


  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* void testInterPolate2DSphereCore()  */

/*----------------------------------------------------------------------
 * FUNCTION void Interpolate2DSphere()
 *
 *
 *
 *
 *
 *----------------------------------------------------------------------*/

void Interpolate2DSphere(REAL8 *finterp, REAL8Vector *xvec, REAL8Vector *yvec, REAL8Array *func, REAL8 xval, REAL8 yval, REAL8 gridsize)
{
  UINT4 idx1, idx2, idx3, idx4, xIdx1, xIdx2, yIdx1, yIdx2, Nx, Ny;
  REAL8 sval,tval,f1,f2,f3,f4,fval;

  Nx = xvec->length;
  Ny = yvec->length;

  xIdx1 = (UINT4) floor(xval/gridsize);
  yIdx1 = (UINT4) floor((yval+LAL_PI_2)/gridsize);

  xIdx2 = xIdx1+1;
  yIdx2 = yIdx1+1;

  if(yIdx2==Ny) {
    xIdx2 = (UINT4) ( ((REAL8) Nx)/2.0 + ((REAL8) xIdx1) );
    yIdx2 = yIdx1;
  }
  if(xIdx2>=Nx) {
    xIdx2 = xIdx2-Nx;
  }

  sval = (xval - xvec->data[xIdx1])/gridsize;
  tval = (yval - yvec->data[yIdx1])/gridsize;

#ifdef TEST_INTERPOLATE2D
#if 0
  printf("%u %u %u %u\n",xIdx1,yIdx1,xIdx2,yIdx2);


  printf("(%g,%g),(%g,%g),(%g,%g),(%g,%g),(%g,%g)\n",
	 xval,yval,
	 xvec->data[xIdx1],yvec->data[yIdx1],
	 xvec->data[xIdx2],yvec->data[yIdx1],
	 xvec->data[xIdx2],yvec->data[yIdx2],
	 xvec->data[xIdx1],yvec->data[yIdx2]
	 );
#endif
#endif

  idx1 = xIdx1 + Nx * yIdx1;
  idx2 = xIdx2 + Nx * yIdx1;
  idx3 = xIdx2 + Nx * yIdx2;
  idx4 = xIdx1 + Nx * yIdx2;

  f1 = func->data[idx1];
  f2 = func->data[idx2];
  f3 = func->data[idx3];
  f4 = func->data[idx4];

#ifdef TEST_INTERPOLATE2D
#if 0
  printf("%u %u %u %u, %u %u, %u %u\n",
	 idx1,idx2,idx3,idx4, Nx, Ny, func->dimLength->data[0],func->dimLength->data[1]
	 );
  printf("%g %g\n",
	 sval,tval
	 );
  printf("%g %g %g %g\n",
	 f1,f2,f3,f4
	 );
  fflush(stdout);
#endif
#endif

  fval = (1.0-sval)*(1.0-tval)*f1;
  fval += sval*(1.0-tval)*f2;
  fval += sval*tval*f3;
  fval += (1.0-sval)*tval*f4; 

  *finterp = fval;

  return;
}  /* void Interpolate2DSphere() */




/*----------------------------------------------------------------------
 * FUNCTION ABCLookUpTable.
 *
 *
 *
 *
 *
 *----------------------------------------------------------------------*/
void
ABCLookUpTable(LALStatus *status, ABCcoeffs *ABCs, BayesFstatParams *bfparams, REAL8 gridsize, UINT4 Ndelta, UINT4 Nalpha)
{
  UINT4 ic,jc, idx, idx2;


  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  ASSERT( ABCs != NULL, status, BAYESFSTATC_ENULL, BAYESFSTATC_MSGENULL);
  ASSERT( bfparams != NULL, status, BAYESFSTATC_ENULL, BAYESFSTATC_MSGENULL);


  bfparams->baryinput.alpha = 0.0;
  bfparams->baryinput.delta = 0.0;

  for(ic=0;ic<Nalpha;ic++) { 
    ABCs->alpha->data[ic] = gridsize * ( (REAL8) ic) ;
  }
  for(ic=0;ic<Ndelta;ic++) { 
    ABCs->delta->data[ic] = gridsize * ( (REAL8) ic) - LAL_PI_2;
  }


/* Not use LALNormalizeSkyPosition */
#ifndef DEBUG_ABCLOOKUPTABLE
  for(jc=0;jc<Ndelta;jc++) { /* loop over delta */
    idx = Nalpha*jc;
    bfparams->baryinput.delta = ABCs->delta->data[jc];
    for(ic=0;ic<Nalpha;ic++) {     /* loop over alpha */
      bfparams->baryinput.alpha = ABCs->alpha->data[ic];
      bfparams->amParams->das->pSource->equatorialCoords.latitude =  bfparams->baryinput.delta;
      bfparams->amParams->das->pSource->equatorialCoords.longitude =  bfparams->baryinput.alpha;
      TRY( LALComputeAM(status->statusPtr, &(bfparams->amc), bfparams->midTS, bfparams->amParams), status); 
      idx2=idx+ic;
      ABCs->As->data[idx2] = bfparams->amc.A;
      ABCs->Bs->data[idx2] = bfparams->amc.B;
      ABCs->Cs->data[idx2] = bfparams->amc.C;
    }     /* loop over alpha */
  }     /* loop over delta */
#endif


/* Use LALNormalizeSkyPosition */
#ifdef DEBUG_ABCLOOKUPTABLE 
  {SkyPosition ThatPosition, ThisPosition;
  ThisPosition.system = COORDINATESYSTEM_EQUATORIAL; 
  for(jc=0;jc<Ndelta;jc++) { /* loop over delta */
    idx = Nalpha*jc;
    ThisPosition.latitude = ABCs->delta->data[jc];
    for(ic=0;ic<Nalpha;ic++) {     /* loop over alpha */
      ThisPosition.longitude = ABCs->alpha->data[ic];
      TRY( LALNormalizeSkyPosition(status->statusPtr, &ThatPosition, &ThisPosition), status);
      bfparams->baryinput.delta = ThatPosition.latitude;
      bfparams->baryinput.alpha = ThatPosition.longitude;
      bfparams->amParams->das->pSource->equatorialCoords.latitude =  bfparams->baryinput.delta;
      bfparams->amParams->das->pSource->equatorialCoords.longitude =  bfparams->baryinput.alpha;
      TRY( LALComputeAM(status->statusPtr, &(bfparams->amc), bfparams->midTS, bfparams->amParams), status); 
      idx2=idx+ic;
      ABCs->As->data[idx2] = bfparams->amc.A;
      ABCs->Bs->data[idx2] = bfparams->amc.B;
      ABCs->Cs->data[idx2] = bfparams->amc.C;
    }     /* loop over alpha */
  }     /* loop over delta */
  }
#endif


#ifdef DEBUG_ABCLOOKUPTABLE
#if 0
  {FILE *fp;
  fp=fopen("ABClut.dat","w");
  if(fp==NULL) fprintf(stderr,"File open error for ABC look up table.\n");
  for(jc=0;jc<Ndelta;jc++) { 
    for(ic=0;ic<Nalpha;ic++) { 
      fprintf(fp,"%g %g %g\n",
	     ABCs->As->data[ic+Nalpha*jc],
	     ABCs->Bs->data[ic+Nalpha*jc],
	     ABCs->Cs->data[ic+Nalpha*jc]);
    }
  }
  fclose(fp);
  }
#endif
#endif


  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* void ABCLookUpTable() */



/*----------------------------------------------------------------------
 * FUNCTION void AllocateABCcoeffs()
 *
 *
 *
 *
 *
 *----------------------------------------------------------------------*/

void 
AllocateABCcoeffs(LALStatus *status, ABCcoeffs *ABCs,UINT4 Ndelta, UINT4 Nalpha)
{

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  ASSERT( ABCs != NULL, status, BAYESFSTATC_ENULL, BAYESFSTATC_MSGENULL);

  ABCs->As = NULL;
  ABCs->Bs = NULL;
  ABCs->Cs = NULL;
  ABCs->dimlen = NULL;
  ABCs->alpha = NULL;
  ABCs->delta = NULL;

  TRY( LALU4CreateVector(status->statusPtr,&(ABCs->dimlen),2), status);
  ABCs->dimlen->data[0]=Ndelta;
  ABCs->dimlen->data[1]=Nalpha;

  LALDCreateVector(status->statusPtr,&(ABCs->alpha),Nalpha);
  BEGINFAIL(status) { TRY( FreeABCcoeffs(status->statusPtr, ABCs), status); } ENDFAIL(status);

  LALDCreateVector(status->statusPtr,&(ABCs->delta),Ndelta);
  BEGINFAIL(status) { TRY( FreeABCcoeffs(status->statusPtr, ABCs), status); } ENDFAIL(status);


  LALDCreateArray(status->statusPtr,&(ABCs->As),ABCs->dimlen);
  BEGINFAIL(status) { TRY( FreeABCcoeffs(status->statusPtr, ABCs), status); } ENDFAIL(status);


  LALDCreateArray(status->statusPtr,&(ABCs->Bs),ABCs->dimlen);
  BEGINFAIL(status) { TRY( FreeABCcoeffs(status->statusPtr, ABCs), status); } ENDFAIL(status);


  LALDCreateArray(status->statusPtr,&(ABCs->Cs),ABCs->dimlen);
  BEGINFAIL(status) { TRY( FreeABCcoeffs(status->statusPtr, ABCs), status); } ENDFAIL(status);


  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* void AllocateABCcoeffs() */



/*----------------------------------------------------------------------
 * FUNCTION void FreeABCcoeffs(LALStatus *status, ABCcoeffs *ABCs)
 *
 *
 *
 *
 *
 *----------------------------------------------------------------------*/


void
FreeABCcoeffs(LALStatus *status, ABCcoeffs *ABCs)
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  ASSERT( ABCs != NULL, status, BAYESFSTATC_ENULL, BAYESFSTATC_MSGENULL);

  if(ABCs->alpha != NULL) 
    LALDDestroyVector(status->statusPtr,&(ABCs->alpha));
  if(ABCs->delta != NULL) 
    LALDDestroyVector(status->statusPtr,&(ABCs->delta));

  if(ABCs->As != NULL) 
    LALDDestroyArray(status->statusPtr,&(ABCs->As));
  if(ABCs->Bs != NULL) 
    LALDDestroyArray(status->statusPtr,&(ABCs->Bs));
  if(ABCs->Cs != NULL) 
    LALDDestroyArray(status->statusPtr,&(ABCs->Cs));

  if(ABCs->dimlen != NULL) 
    LALU4DestroyVector(status->statusPtr,&(ABCs->dimlen));

  if(ABCs != NULL) 
    LALFree(ABCs);

  DETATCHSTATUSPTR (status);
  RETURN (status);
} /* void FreeABCcoeffs() */




/*----------------------------------------------------------------------
 * FUNCTION DInverseWithIntMeasure().
 *
 *
 *
 *
 *
 *----------------------------------------------------------------------*/

REAL8 
DInverseWithIntMeasure(REAL8 *k,
	 size_t dim,
	 void *params)
{
  REAL8 Dinv = 0.0;

  /* 
     k[0]=alpha;
     k[1]=delta;
     k[2]=mu;
     k[3]=psi;
  */

  Dinv = 1.0/DFunction(k, dim, params);
  /* integral measure */
  Dinv *= cos(k[1]);

  return Dinv;
} /* REAL8 DInverseWithIntMeasure() */


/*----------------------------------------------------------------------
 * FUNCTION DSquareWithIntMeasure().
 *
 *
 *
 *
 *
 *----------------------------------------------------------------------*/

REAL8 
DSquareWithIntMeasure(REAL8 *k,
	 size_t dim,
	 void *params)
{
  REAL8 D2 = 0.0;

  /* 
     k[0]=alpha;
     k[1]=delta;
     k[2]=mu;
     k[3]=psi;
  */



  D2 = pow(DFunction(k, dim, params),2);
  /* integral measure for the integration over delta. */
  D2 *= cos(k[1]);
  /*
  printf("%g %g %g %g %g\n",k[0],k[1],k[2],k[3],D2);
  */
  return D2;
} /* REAL8 DSquareWithIntMeasure() */


/*----------------------------------------------------------------------
 * FUNCTION REAL8 integrateDSquare
 * 
 * + Integrate D-squre for a checking and a demonstration purpose.
 *
 *
 *----------------------------------------------------------------------*/
void integrateDSquare(LALStatus *status, BayesFstatParams *bfparams)
{
  MonteCarloIntOut mciout;
  MonteCarloIntParams mciparams;
  MonteCarloIntIn mciin;

  REAL8 norm;

  INITSTATUS(status);
  ASSERT( bfparams != NULL, status, BAYESFSTATC_ENULL, BAYESFSTATC_MSGENULL);

  mciparams.verboseflag = lalDebugLevel;
  mciparams.method = "vegas";
  mciparams.warmupcalls = 10000;
  mciparams.maincalls = 500000;
  mciparams.maxVegasStages = 100;

  mciin.params = (void *) bfparams;

  mciin.domain = (REAL8 *) LALCalloc(10,sizeof(REAL8));
  if(mciin.domain == NULL) {
    ABORT (status, BAYESFSTATC_EMEM, BAYESFSTATC_MSGEMEM); 
  }
  mciin.dimension = 4;
  /* alpha  */
  mciin.domain[0] = 0.0;        mciin.domain[1] = LAL_TWOPI;
  /* delta */
  mciin.domain[2] = -LAL_PI_2;  mciin.domain[3] = LAL_PI_2;
  /* cosiota */
  mciin.domain[4] = 0.0;        mciin.domain[5] = 1.0;
  /* psi */
  mciin.domain[6] = - LAL_PI_4; mciin.domain[7] = LAL_PI_4;

  /*
  {REAL8 k[4], res;
  k[0]=3.15986; k[1]=1.56922;k[2]=0.997443; k[3]=-0.205399;
  res=DSquareWithIntMeasure(k,4,bfparams);
  printf("%g\n",res);
  }
  */

  /* just for check. This, after normalization, should give us 4/25 = 1.6, as the JKS paper shows. */ 
  mciin.integrand = DSquareWithIntMeasure;
  MonteCarloIntegrate(&mciout, &mciin, &mciparams);
  fprintf(stdout,"%% \n%% Testing the code ....");
  fprintf(stdout,"%% Integrating d^2 of JKS. \n");
  fprintf(stdout,"%% The correct answer = 4/25 = 0.16\n");
  norm = LAL_TWOPI * 2.0 * 1.0 * LAL_PI_2;
  fprintf(stdout,"%% numerical integration result = %22.16g\n",mciout.result/norm);
  fprintf(stdout,"%% estimated absolute error (1 standard deviation) is = %22.16g\n",mciout.sigma);
  fprintf(stdout,"%% estimated relative error (1 standard deviation) is = %22.16g\n",mciout.sigma/mciout.result);

  LALFree(mciin.domain);

  RETURN (status); 
} /* void integrateDSquare() */




