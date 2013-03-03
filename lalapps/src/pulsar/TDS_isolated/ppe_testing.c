/**
 * \file
 * \ingroup pulsarApps
 * \author Matthew Pitkin
 *
 * \brief Functions for use in testing the parameter estimation codes for
 * targeted pulsar searches.
 */

#include "ppe_testing.h"

/* *****************************************************************************/
/*                          TESTING FUNCTIONS                                 */
/* *****************************************************************************/

/** \brief A test function to calculate a 1D posterior on a grid
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
  REAL8 h0min = 0.;
  REAL8 h0max = 0.;
  REAL8 h0range = 0, h0step = 0;
  INT4 h0steps = 0, i = 0;
  UINT4 verbose =
    LALInferenceCheckVariable(runState->algorithmParams,"verbose");
  
  ProcessParamsTable *ppt;
  REAL8 scaleval = 1., minval = 0., tmpscale = 0., tmpmin = 0., tmpgridval = 0.;
  
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
  tmpscale = *(REAL8*)LALInferenceGetVariable( runState->data->dataParams,
                                               parscale );
  tmpmin = *(REAL8*)LALInferenceGetVariable( runState->data->dataParams,
                                             parmin );
  LALInferenceRemoveVariable( runState->data->dataParams, parscale );
  LALInferenceAddVariable( runState->data->dataParams, parscale, &scaleval,
                           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );
  LALInferenceRemoveVariable( runState->data->dataParams, parmin );
  LALInferenceAddVariable( runState->data->dataParams, parmin, &minval,
                           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );
  
  tmpgridval = *(REAL8*)LALInferenceGetVariable( runState->currentParams,
                                                 parname );
  
  sprintf(outputgrid, "%s_grid_posterior.txt", parname);
  
  if ( (fp = fopen(outputgrid, "w")) == NULL ){
    fprintf(stderr, "Error... cannot open grid posterior file %s.\n",
            outputgrid);
    exit(0);
  }
  
  for( i = 0; i < h0steps; i++ ){
    REAL8 h0val = h0min + i*h0step;
    
    LALInferenceSetVariable( runState->currentParams, parname, &h0val );
    
    logL->data[i] = runState->likelihood( runState->currentParams,
                                          runState->data, runState->templt );
    
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
  LALInferenceRemoveVariable( runState->data->dataParams, parscale );
  LALInferenceAddVariable( runState->data->dataParams, parscale, &tmpscale,
                           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );
  
  LALInferenceRemoveVariable( runState->data->dataParams, parmin );
  LALInferenceAddVariable( runState->data->dataParams, parmin, &tmpmin,
                           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );
                           
  LALInferenceSetVariable( runState->currentParams, parname, &tmpgridval );
}

/** \brief Test the sampler using a Gaussian likelihood
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
                                    LALInferenceIFOData *data,
                                    LALInferenceTemplateFunction UNUSED
                                      get_model ){
  REAL8 loglike = 0.; /* the log likelihood */
  
  REAL8 like_mean = 0.5;
  REAL8 like_sigma = 0.025;
  REAL8 h0 = *(REAL8 *)LALInferenceGetVariable( vars, "h0" );
  REAL8 h0scale = *(REAL8 *)LALInferenceGetVariable( data->dataParams,
                                                     "h0_scale" );
  REAL8 h0min = *(REAL8 *)LALInferenceGetVariable( data->dataParams,
                                                   "h0_scale_min" );
  
  get_model = NULL;
                                                   
  h0 = h0*h0scale + h0min;

  /* search over a simple 1D Gaussian with x defined by the h0 variable */
  loglike = -log(sqrt(2.*LAL_PI)*like_sigma);
  loglike -= (h0-like_mean)*(h0-like_mean) / (2.*like_sigma*like_sigma);
  
  return loglike;
}


/** \brief Output a number of prior samples based on the initial live points
 * 
 * This function will output prior samples for variable parameters (as create by
 * the LALInferenceSetupLivePointsArray function) making sure to rescale the
 * values.
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
          /* rescale sample */
          CHAR scalePar[VARNAME_MAX] = "";
          CHAR scaleMinPar[VARNAME_MAX] = "";
          REAL8 scale = 0., scaleMin = 0., var = 0.;
          
          /* get scale factors */
          sprintf(scalePar, "%s_scale", item->name);
          scale = *(REAL8 *)LALInferenceGetVariable( runState->data->dataParams,
                                                     scalePar );
    
          sprintf(scaleMinPar, "%s_scale_min", item->name);
          scaleMin = *(REAL8 *)LALInferenceGetVariable(
            runState->data->dataParams, scaleMinPar );
          
          var = scaleMin + *(REAL8 *)item->value*scale;

          fprintf(fp, "%.8le\t", var);
        }
      
        item = item->next;
      }
      fprintf(fp, "\n");
    }
  
    fclose(fp);
  }
}

/*----------------------- END OF TESTING FUNCTIONS ---------------------------*/
