/*
*  Copyright (C) 2014 Matthew Pitkin
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

/**
 * \file
 * \ingroup lalapps_pulsar_HeterodyneSearch
 * \author Matthew Pitkin
 *
 * \brief Functions for use in testing the parameter estimation codes for
 * targeted pulsar searches.
 */

#include "config.h"
#include "ppe_testing.h"

/* *****************************************************************************/
/*                          TESTING FUNCTIONS                                 */
/* *****************************************************************************/

/* upper limit calculation distribution helper function prototypes */
static double ul_gauss_cdf_function( double x, void *params );
static double ul_gauss_CDFRoot( double mu, double sigma, double min, double max );

/**
 * \brief A test function to calculate a 1D posterior on a grid
 *
 * This function is only to be used as a check/test of the code and will be run
 * if the \c grid command line argument is present. It will calculate the
 * posterior for one parameter (given by \c gridpar), between the ranges given
 * by \c gridmin and \c gridmax (which default to 0 and 1) at a number of points
 * given by \c gridsteps (which default to 100).
 *
 * \param runState [in] The analysis information structure
 */
void gridOutput( LALInferenceRunState *runState ){
  LALInferenceThreadState *threadState=runState->threads[0];
  REAL8 h0min = 0.;
  REAL8 h0max = 0.;
  REAL8 h0range = 0, h0step = 0;
  INT4 h0steps = 0, i = 0;
  UINT4 verbose = LALInferenceCheckVariable(runState->algorithmParams, "verbose");

  ProcessParamsTable *ppt;
  REAL8 scale = 1., scalemin = 0., tmpscale = 0., tmpmin = 0., tmpgridval = 0.;

  ProcessParamsTable *commandLine = runState->commandLine;

  FILE *fp = NULL;
  REAL8 minL = LAL_REAL8_MAX;
  REAL8 sumPost = 0.;

  REAL8Vector *logL = NULL;

  CHAR *parname = NULL, parscale[256], parmin[256], outputgrid[256];

  /*------------------------------------------------------------*/
  /* test output on a h0 grid */
  ppt = LALInferenceGetProcParamVal( commandLine, "--grid" );
  if ( ppt ){
    ProcessParamsTable *ppt2;

    /* parameters over which to perform the grid search */
    ppt2 = LALInferenceGetProcParamVal( commandLine, "--gridpar" );

    if( ppt2 ){
      parname = XLALStringDuplicate( ppt2->value );

      if( !recognised_parameter( parname ) ){
        fprintf(stderr, "Error... parameter %s not recognised\n", parname );
        exit(0);
      }

      sprintf(parscale, "%s_scale", parname);
      sprintf(parmin, "%s_scale_min", parname);
    }
    else{
      fprintf(stderr, USAGEGRID, commandLine->program);
      exit(0);
    }

    ppt2 = LALInferenceGetProcParamVal( commandLine, "--gridmin" );

    if( ppt2 ) h0min = atof( ppt2->value );
    else h0min = 0.; /* default to zero */

    ppt2 = LALInferenceGetProcParamVal( commandLine, "--gridmax" );

    if( ppt2 ) h0max = atof( ppt2->value );
    else h0max = 1.; /* default to 1 */

    ppt2 = LALInferenceGetProcParamVal( commandLine, "--gridsteps" );

    if( ppt2 ) h0steps = atoi( ppt2->value );
    else h0steps = 100; /* default to 100 steps */
  }
  else{
    return;
  }

  if ( verbose ){
    fprintf(stderr, "Calculating posterior on %s over a grid from:\n", parname);
    fprintf(stderr, "\t%le --> %le in %d steps.\n", h0min, h0max, h0steps);
  }

  h0range = h0max - h0min;
  h0step = h0range / (REAL8)(h0steps-1.);

  logL = XLALCreateREAL8Vector( h0steps );

  /* reset rescale value for h0 */
  tmpscale = *(REAL8*)LALInferenceGetVariable( threadState->model->ifo->params,
                                               parscale );
  tmpmin = *(REAL8*)LALInferenceGetVariable( threadState->model->ifo->params,
                                             parmin );
  LALInferenceRemoveVariable( threadState->model->ifo->params, parscale );
  LALInferenceAddVariable( threadState->model->ifo->params, parscale, &scale,
                           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );
  LALInferenceRemoveVariable( threadState->model->ifo->params, parmin );
  LALInferenceAddVariable( threadState->model->ifo->params, parmin, &scalemin,
                           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );

  tmpgridval = *(REAL8*)LALInferenceGetVariable( threadState->currentParams,
                                                 parname );

  sprintf(outputgrid, "%s_grid_posterior.txt", parname);

  if ( (fp = fopen(outputgrid, "w")) == NULL ){
    fprintf(stderr, "Error... cannot open grid posterior file %s.\n",
            outputgrid);
    exit(0);
  }

  for( i = 0; i < h0steps; i++ ){
    REAL8 h0val = h0min + i*h0step;

    LALInferenceSetVariable( threadState->currentParams, parname, &h0val );

    logL->data[i] = runState->likelihood( threadState->currentParams,
                                          runState->data, threadState->model );

    if ( logL->data[i] < minL ) minL = logL->data[i];
  }

  /* integrate area under posterior - trapezium rule */
  for( i = 0; i < h0steps-1; i++ ){
    sumPost += ( exp(logL->data[i] - minL) + exp(logL->data[i+1] - minL) ) *
      h0step / 2.;
  }

  /* output posterior */
  for( i = 0; i < h0steps; i++ ){
    REAL8 h0val = h0min + i*h0step;
    fprintf(fp, "%le\t%le\n", h0val, exp( logL->data[i] - minL ) / sumPost);
  }

  fclose(fp);

  XLALDestroyREAL8Vector( logL );

  /* reset scale value and parameter value in currentParams */
  LALInferenceRemoveVariable( threadState->model->ifo->params, parscale );
  LALInferenceAddVariable( threadState->model->ifo->params, parscale, &tmpscale,
                           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );

  LALInferenceRemoveVariable( threadState->model->ifo->params, parmin );
  LALInferenceAddVariable( threadState->model->ifo->params, parmin, &tmpmin,
                           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );

  LALInferenceSetVariable( threadState->currentParams, parname, &tmpgridval );
}

/**
 * \brief Test the sampler using a Gaussian likelihood
 *
 * This is a testing function that can be substituted for the standard
 * likelihood function. It calculates only the \c h0 parameter posterior based
 * on a Gaussian likelihood with mean of 0.0 and standard deviation of 0.025 -
 * these values can be changed if required. It is just to be used to test the
 * sampling routine (e.g. Nested Sampling) with a well defined likelihood
 * function.
 *
 * \param vars [in] A set of pulsar parameters
 * \param data [in] A data structure
 * \param get_model UNDOCUMENTED
 *
 * \return Natural logarithm of the likelihood
 */
REAL8 test_gaussian_log_likelihood( LALInferenceVariables *vars,
                                    LALInferenceIFOData *data,
                                    LALInferenceModel *get_model ){
  REAL8 loglike = 0.; /* the log likelihood */

  REAL8 like_mean = LALInferenceGetREAL8Variable( vars, "H0MEAN" );
  REAL8 like_sigma = LALInferenceGetREAL8Variable( vars, "H0SIGMA" );

  if ( !LALInferenceCheckVariable( vars, "H0" ) ){
    fprintf(stderr, "Error... testing Gaussian likelihood required the \"H0\" parameter to be set in the prior file.\n");
    exit(1);
  }

  REAL8 h0 = LALInferenceGetREAL8Variable( vars, "H0" );

  /* search over a simple 1D Gaussian with x defined by the h0 variable */
  loglike = -log(sqrt(2.*LAL_PI)*like_sigma);
  loglike -= 0.5*(h0-like_mean)*(h0-like_mean) / (like_sigma*like_sigma);
  get_model->ifo_loglikelihoods[0] = loglike;

  data->likeli_counter += 1;

  return loglike;
}


/**
 * \brief Output the analytic evidence for the test Gaussian likelihood and a 95% upper limit
 *
 * This function calculated the analytical evidence for a given test Gaussian likelihood and
 * also a 95% upper limit on the likelihood. The values are output to the HDF5 file if given,
 * but otherwise are just output to screen.
 *
 * \param runState [in] The LALInference run state variable
 */
void test_gaussian_output( LALInferenceRunState *runState ){
  ProcessParamsTable *ppt = NULL;
  FILE *fp = NULL;

  /* get the minimum and maximum ranges for the Gaussian */
  REAL8 min = 0., max = INFINITY, Z;
  if( LALInferenceCheckMinMaxPrior( runState->priorArgs, "H0" ) ){
    LALInferenceGetMinMaxPrior( runState->priorArgs, "H0", &min, &max );
  }
  else{
    fprintf(stderr, "Error... no prior range set for the test Gaussian likelihood");
    exit(1);
  }

  /* get mean and standard deviation */
  REAL8 h0mean = 0., h0sigma = 0., h0sigma2 = 0., h095 = 0.;
  h0mean = LALInferenceGetREAL8Variable( runState->threads[0]->currentParams, "H0MEAN" );
  h0sigma = LALInferenceGetREAL8Variable( runState->threads[0]->currentParams, "H0SIGMA" );
  h0sigma2 = h0sigma*h0sigma;

  /* calculate the log evidence */
  Z = log(0.5*(erf(LAL_SQRT1_2*(h0mean - min)/h0sigma) - erf(LAL_SQRT1_2*(h0mean - max)/h0sigma)));
  Z -= log(max-min); /* multiply by prior */

  /* calculate the 95% upper limit */
  h095 = ul_gauss_CDFRoot(h0mean, h0sigma, min, max);

  /* calculate the KL divergence */
  REAL8 lnC = -log(max-min); /* log of prior */
  REAL8 p_Z = exp(lnC - Z); /* prior divided by the evidence */
  REAL8 L = Z + 0.5*log(LAL_TWOPI*h0sigma2);
  REAL8 D = (1. + 2.*L)*(erf((h0mean-min)*LAL_SQRT1_2/h0sigma) - erf((h0mean-max)*LAL_SQRT1_2/h0sigma));
  REAL8 G = (1./(sqrt(LAL_TWOPI)*h0sigma))*((min-h0mean)*exp(-0.5*(min-h0mean)*(min-h0mean)/h0sigma2) - (max-h0mean)*exp(-0.5*(max-h0mean)*(max-h0mean)/h0sigma2));
  REAL8 KLdiv = -0.25*p_Z*(D + 2.*G);

  /* Open output file (called test_gauss.txt) using the path of the --outfile value */
  ppt = LALInferenceGetProcParamVal(runState->commandLine, "--outfile");
  char *outfile = XLALStringDuplicate(ppt->value), *loc = NULL;
  /* find last '/' and replace with null termination character */
  loc = strrchr(outfile, '/');
  if ( loc ){ outfile[strlen(outfile)-strlen(loc)+1] = '\0'; }
  else{ outfile = NULL; }
  outfile = XLALStringAppend(outfile, "test_gauss.txt");

  if ( ( fp = fopen(outfile, "w") ) != NULL ){
    fprintf(fp, "%.12le\t%.12le\t%.12le\n", Z, h095, KLdiv);
    fclose(fp);
  }
  else{
    fprintf(stderr, "Warning... could not open test Gaussian output file '%s'.", outfile);
  }
}


typedef struct
tagul_params{
  double mu;
  double sigma;
  double min;
  double area;
} ul_params;

/* internal Gaussian CDF function for 95% upper limit finding */
double ul_gauss_cdf_function( double x, void *params ){
  ul_params *p = (ul_params*)params;

  double mu = p->mu;
  double sigma = p->sigma;
  double min = p->min;
  double area = p->area;
  double ul = 0.95; /* calculating 95% upper limit */
  double C = 1./(sqrt(2.)*sigma);

  return (0.5*(erf(C*(mu-min)) - erf(C*(mu-x)))/area) - ul;
}


/** \brief Find the root of the Gaussian CDF to give a 95% upper limit
 *
 * Use the Steffenson method to find the root of the function defined in \c ul_gauss_function.
 */
static double ul_gauss_CDFRoot( double mu, double sigma, double min, double max ){
  int gslstatus;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double x_lo, x_hi;
  double epsrel = 1e-4; /* relative error tolerance */

  double C = 1./(sqrt(2.)*sigma);
  double area = 0.5*(erf(C*(mu - min)) - erf(C*(mu - max)));

  gsl_function F;
  ul_params params = {mu, sigma, min, area};

  double x;
  /* initial bounds of value */
  x_lo = mu-10.*sigma;
  x_hi = mu+10.*sigma;

  F.function = &ul_gauss_cdf_function;
  F.params = &params;

  T = gsl_root_fsolver_brent; /* use Brent method */
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  do{
    iter++;
    gslstatus = gsl_root_fsolver_iterate( s );
    x = gsl_root_fsolver_root( s );
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);

    /* test relative error of bounds */
    gslstatus = gsl_root_test_interval (x_lo, x_hi, 0., epsrel);
  }
  while( gslstatus == GSL_CONTINUE && iter < max_iter );

  if ( gslstatus != GSL_SUCCESS ){
    XLALPrintError("%s: Failed to converge when drawing from Fermi-Dirac distribution.", __func__);
    XLAL_ERROR_REAL8(XLAL_EFAILED);
  }

  gsl_root_fsolver_free (s);

  return x;
}


/**
 * \brief Output a number of prior samples based on the initial live points
 *
 * This function will output prior samples for variable parameters (as create by
 * the LALInferenceSetupLivePointsArray function).
 *
 * \param runState [in]
 */
void outputPriorSamples( LALInferenceRunState *runState ){
  ProcessParamsTable *ppt = LALInferenceGetProcParamVal( runState->commandLine,
                                                         "--output-prior" );
  INT4 Nlive = *(INT4 *)LALInferenceGetVariable( runState->algorithmParams,
                                                 "Nlive" );

  if( ppt ){
    FILE *fp = NULL;

    /* loop over the live points */
    INT4 i = 0;

    if( ( fp = fopen("prior_samples.txt", "w") ) == NULL ){
      fprintf(stderr, "Error... could not open prior samples file\n");
      exit(1);
    }

    for ( i = 0; i < Nlive; i++ ){
      /* output variable parameters, rescaled to their proper ranges */
      LALInferenceVariableItem *item = runState->livePoints[i]->head;

      /* print out variables names */
      if ( i == 0 ){
        fprintf(fp, "%% ");

        while( item ){
          if( item->vary == LALINFERENCE_PARAM_LINEAR ||
              item->vary == LALINFERENCE_PARAM_CIRCULAR )
            fprintf(fp, "%s\t", item->name);

          item = item->next;
        }

        fprintf(fp, "\n");
      }

      item = runState->livePoints[i]->head; /* reset */

      while( item ){
        if( item->vary == LALINFERENCE_PARAM_LINEAR ||
            item->vary == LALINFERENCE_PARAM_CIRCULAR ){
          REAL8 var = 0.;
          var = *(REAL8 *)item->value;
          fprintf(fp, "%.16le\t", var);
        }

        item = item->next;
      }
      fprintf(fp, "\n");
    }

    fclose(fp);
  }
}


/**
 * \brief Read in an ascii text file of nested samples and compare the log likelihoods
 *
 * This function reads in a file containing nested samples, including their log likelihoods, recomputes
 * the full likelihood and compares it to the previously output value. This is useful for comparing
 * outputs using ROQ to outputs using the full likelihood.
 */
void compare_likelihoods( LALInferenceRunState *rs ){
  ProcessParamsTable *ppt = NULL;
  CHAR *sampfile = NULL, *namefile = NULL;
  CHAR *parambuf = NULL, *databuf = NULL;
  LALInferenceVariables *curparams = NULL;
  UINT4 j = 0, k = 0;
  FILE *fp = NULL;
  REAL8 logLnew = 0., logL = 0.;

  /* allocate memory for nested samples */
  curparams = XLALCalloc( 1, sizeof(LALInferenceVariables) );
  LALInferenceCopyVariables( rs->threads[0]->currentParams, curparams );

  /* read in parameter names from header file */
  ppt = LALInferenceGetProcParamVal( rs->commandLine, "--sample-file-header" );
  if ( ppt != NULL ){ namefile = ppt->value; }
  else{ XLAL_ERROR_VOID(XLAL_EINVAL, "Error... no nested samples header file given.\n"); }
  parambuf = XLALFileLoad( namefile );
  TokenList *paramNames = NULL;
  /* seperate parameter names */
  if ( XLALCreateTokenList( &paramNames, parambuf, " \t" ) != XLAL_SUCCESS ){
    XLAL_ERROR_VOID(XLAL_EINVAL, "Error... could not read in parameter names.\n");
  }

  /* get names of nested sample file columns */
  ppt = LALInferenceGetProcParamVal( rs->commandLine, "--sample-file" );
  if ( ppt != NULL ){ sampfile = ppt->value; }
  else{ XLAL_ERROR_VOID(XLAL_EINVAL, "Error... no nested samples file given.\n"); }
  databuf = XLALFileLoad( sampfile );
  TokenList *sampsList = NULL;
  /* seperate samples */
  if ( XLALCreateTokenList( &sampsList, databuf, "\n" ) != XLAL_SUCCESS ){
    XLAL_ERROR_VOID(XLAL_EINVAL, "Error... could not read in nested samples.\n");
  }

  fp = fopen("likelihoodComp.txt", "w");

  /* loop over samples and calculate likelihoods */
  for ( k = 0; k < sampsList->nTokens; k++){
    /* remove "fixed variable" logL from curparams */
    if ( LALInferenceCheckVariable(curparams, "logL" ) ){ LALInferenceRemoveVariable( curparams, "logL" ); }

    /* split line */
    TokenList *paramVals = NULL;
    if ( XLALCreateTokenList( &paramVals, sampsList->tokens[k], " \t" ) != XLAL_SUCCESS ){
      XLAL_ERROR_VOID(XLAL_EINVAL, "Error... could not separate parameter values.\n");
    }

    /* copy parameters into LALInferenceVariables structure */
    for ( j = 0; j < paramNames->nTokens; j++ ){
      REAL8 value = atof(paramVals->tokens[j]);
      LALInferenceAddVariable( curparams, paramNames->tokens[j], &value, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );
    }

    if ( !LALInferenceCheckVariable(curparams, "logL" ) ){
      XLAL_ERROR_VOID(XLAL_EINVAL, "Error... no log likelihood value present in nested samples file.\n");
    }

    logL = LALInferenceGetREAL8Variable( curparams, "logL" );
    logLnew = rs->likelihood( curparams, rs->data, rs->threads[0]->model );
    fprintf(stderr, "%.16le\t%.16le\t%.16le\n", logL, logLnew, logL-logLnew);
    fprintf(fp, "%.16le\t%.16le\t%.16le\n", logL, logLnew, logL-logLnew);

    XLALDestroyTokenList( paramVals );
  }

  fclose(fp);

  XLALFree( parambuf );
  XLALFree( databuf );
  XLALDestroyTokenList( sampsList );
  XLALDestroyTokenList( paramNames );
}

/*----------------------- END OF TESTING FUNCTIONS ---------------------------*/
