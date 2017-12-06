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
 * on a Gaussian likelihood with mean of 0.5 and standard deviation of 0.025 -
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
                                    UNUSED LALInferenceIFOData *data,
                                    LALInferenceModel *get_model ){
  REAL8 loglike = 0.; /* the log likelihood */

  REAL8 like_mean = 0.5;
  REAL8 like_sigma = 0.025;
  REAL8 h0 = *(REAL8 *)LALInferenceGetVariable( vars, "h0" );
  REAL8 h0scale = *(REAL8 *)LALInferenceGetVariable( get_model->ifo->params, "h0_scale" );
  REAL8 h0min = *(REAL8 *)LALInferenceGetVariable( get_model->ifo->params, "h0_scale_min" );

  get_model = NULL;

  h0 = h0*h0scale + h0min;

  /* search over a simple 1D Gaussian with x defined by the h0 variable */
  loglike = -log(sqrt(2.*LAL_PI)*like_sigma);
  loglike -= (h0-like_mean)*(h0-like_mean) / (2.*like_sigma*like_sigma);

  return loglike;
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
