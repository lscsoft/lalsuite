/**
 * \file
 * \ingroup pulsarApps
 * \author Matthew Pitkin, John Veitch, Colin Gill
 *
 * \brief Pulsar likelihood and prior functions for use in parameter estimation
 * codes for targeted pulsar searches.
 */

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include "ppe_likelihood.h"

/******************************************************************************/
/*                     LIKELIHOOD AND PRIOR FUNCTIONS                         */
/******************************************************************************/
/** \brief The log likelihood function
 * 
 * This function calculates natural logarithm of the likelihood of a signal
 * model (specified by a given set of parameters) given the data from a set of 
 * detectors.
 * 
 * The likelihood is the joint likelihood of chunks of data over which the noise
 * is assumed stationary and Gaussian. For each chunk a Gaussian likelihood for
 * the noise and data has been marginalised over the unknown noise standard 
 * deviation using a Jeffreys prior on the standard deviation. Given the
 * data consisting of independent real and imaginary parts this gives a
 * Students-t distribution for each chunk (of length \f$m\f$) with \f$m/2\f$
 * degrees of freedom:
 * \f[
 * p(\mathbf{\theta}|\mathbf{B}) = \prod_{j=1}^M \frac{(m_j-1)!}{2\pi^{m_j}} 
 * \left( \sum_{k=k_0}^{k_0+(m_j-1)} |B_k - y(\mathbf{\theta})_k|^2 
 * \right)^{-m_j},
 * \f]
 * where \f$\mathbf{B}\f$ is a vector of the complex data, 
 * \f$y(\mathbf{\theta})\f$ is the model for a set of parameters
 * \f$\mathbf{\theta}\f$, \f$M\f$ is the total number of independent data chunks
 * with lengths \f$m_j\f$ and \f$k_0 = \sum_{i=1}^j 1 + m_{i-1}\f$ (with
 * \f$m_0 = 0\f$) is the index of the first data point in each chunk. The
 * product of this for each detector will give the full joint likelihood. In
 * the calculation here the unnecessary proportionality factors are left out
 * (this would effect the actual value of the marginal likelihood/evidence, but
 * since we are only interested in evidence ratios/Bayes factors these
 * factors would cancel out anyway. See [\ref DupuisWoan2005] for a more
 * detailed description.
 *
 * In this function data in chunks smaller than a certain minimum length 
 * \c chunkMin are ignored.
 * 
 * \param vars [in] The parameter values
 * \param data [in] The detector data and initial signal phase template
 * \param get_model [in] The signal template/model function
 * 
 * \return The natural logarithm of the likelihood function
 */
REAL8 pulsar_log_likelihood( LALInferenceVariables *vars, 
                             LALInferenceIFOData *data,
                             LALInferenceTemplateFunction *get_model ){
  REAL8 loglike = 0.; /* the log likelihood */
  UINT4 i = 0;
  CHAR *modeltype = NULL;/*need to check model type in this function*/
  
  modeltype = *(CHAR**)LALInferenceGetVariable( data->dataParams, "modeltype" );
  
  LALInferenceIFOData *datatemp1 = data, *datatemp2 = data;
 
  /* copy model parameters to data parameters */
  while( datatemp1 ){
    LALInferenceCopyVariables( vars, datatemp1->modelParams );
    datatemp1 = datatemp1->next;
  }
 
  /* get pulsar model */
  while( datatemp2 ){
    /*fprintf(bugtest,"getting model in log like func\n");*/
    get_model( datatemp2 );
    datatemp2 = datatemp2->next;
    /* If modeltype is pinsf need to advance data on to next, so this loop only
     runs once if there is only 1 det*/
    if ( !strcmp( modeltype, "pinsf" ) ) datatemp2 = datatemp2->next;
  }

  while ( data ){
    UINT4 j = 0, count = 0, cl = 0;
    UINT4 length = 0, chunkMin;
    REAL8 chunkLength = 0.;
    REAL8 logliketmp = 0.;

    REAL8 sumModel = 0., sumDataModel = 0.;
    REAL8 chiSquare = 0.;
    COMPLEX16 B, M;
  
    REAL8Vector *sumDat = NULL;
    UINT4Vector *chunkLengths = NULL;
    
    /*fprintf(bugtest,"calc log like for one set of data\n");*/
    
    sumDat = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams, 
                                                       "sumData" );
    chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( data->dataParams,
                                                             "chunkLength" );
    chunkMin = *(INT4*)LALInferenceGetVariable( data->dataParams, "chunkMin" );
  
    length = data->compTimeData->data->length;
  
    for( i = 0 ; i < length ; i += chunkLength ){
      chunkLength = (REAL8)chunkLengths->data[count];
    
      /* skip section of data if its length is less than the minimum allowed
        chunk length */
      if( chunkLength < chunkMin ){
        count++;
        continue;
      }

      sumModel = 0.;
      sumDataModel = 0.;

      cl = i + (INT4)chunkLength;
    
      for( j = i ; j < cl ; j++ ){
        B.re = data->compTimeData->data->data[j].re;
        B.im = data->compTimeData->data->data[j].im;

        M.re = data->compModelData->data->data[j].re;
        M.im = data->compModelData->data->data[j].im;
      
        /* sum over the model */
        sumModel += M.re*M.re + M.im*M.im;
        
        /* sum over that data and model */
        sumDataModel += B.re*M.re + B.im*M.im;
        /*fprintf(bugtest,"B.re= %e, B.im= %e, M.re: %e, M.im: %e\n",B.re, B.im,
M.re, M.im);*/
      }
 
      chiSquare = sumDat->data[count];
      chiSquare -= 2.*sumDataModel;
      chiSquare += sumModel;
      
      logliketmp -= chunkLength*log(chiSquare);
      
      count++;
    }
    loglike += logliketmp;
    data = data->next;
  }

  return loglike;
}


/** \brief The prior function
 *
 * This function calculates the natural logarithm of the prior for a set of
 * parameters. If the prior on a particular parameter is uniform over a given
 * range then \f$p(\theta) = 1/(\theta_{\rm max} - \theta_{\rm min})\f$. If the 
 * prior is Gaussian then the probability of that value given the mean and 
 * standard deviation of the Gaussian is calculated.
 * 
 * \param runState [in] A pointer to the LALInferenceRunState
 * \param params [in] The set of parameter values
 * 
 * \return The natural logarithm of the prior value for a set of parameters 
 */
REAL8 priorFunction( LALInferenceRunState *runState, 
                     LALInferenceVariables *params ){
  LALInferenceIFOData *data = runState->data;
  (void)runState;
  LALInferenceVariableItem *item = params->head;
  REAL8 min, max, mu, sigma, prior = 0, value = 0.;

  REAL8Vector *corVals = NULL;
  UINT4 cori = 0;

  /* if some correlated priors exist allocate corVals */
  if ( corlist ) corVals = XLALCreateREAL8Vector( corlist->length );
  
  for(; item; item = item->next ){
    /* get scale factor */
    CHAR scalePar[VARNAME_MAX] = "";
    CHAR scaleMinPar[VARNAME_MAX] = "";
    REAL8 scale = 0., scaleMin = 0.;
    
    if( item->vary == LALINFERENCE_PARAM_FIXED || 
      item->vary == LALINFERENCE_PARAM_OUTPUT ){ continue; }
    
    sprintf(scalePar, "%s_scale", item->name);
    scale = *(REAL8 *)LALInferenceGetVariable( data->dataParams, scalePar );
    
    sprintf(scaleMinPar, "%s_scale_min", item->name);
    scaleMin = *(REAL8 *)LALInferenceGetVariable( data->dataParams, 
                                                  scaleMinPar );
    
    if( item->vary == LALINFERENCE_PARAM_LINEAR || 
      item->vary == LALINFERENCE_PARAM_CIRCULAR ){
      /* Check for a gaussian */
      if ( LALInferenceCheckGaussianPrior(runState->priorArgs, item->name) ){
        LALInferenceGetGaussianPrior( runState->priorArgs, item->name, &mu,
                                      &sigma );
      
       value = (*(REAL8 *)item->value) * scale + scaleMin;
       mu += scaleMin;
       sigma *= scale;
       prior -= log(sqrt(2.*LAL_PI)*sigma);
       prior -= (value - mu)*(value - mu) / (2.*sigma*sigma);
      }
      /* check for a flat prior */
      else if( LALInferenceCheckMinMaxPrior(runState->priorArgs, item->name) ){
        LALInferenceGetMinMaxPrior( runState->priorArgs, item->name, 
                                    &min, &max );
      
        if( (*(REAL8 *) item->value) < min || (*(REAL8 *)item->value) > max ){
          return -INFINITY;
        }
        else prior -= log( (max - min) * scale );
      }
      else if( LALInferenceCheckCorrelatedPrior(runState->priorArgs,
        item->name) && corlist ){      
        /* set item in correct position given the order of the correlation
           matrix given by corlist */
        for( cori = 0; cori < corlist->length; cori++ ){
          if( !strcmp(item->name, corlist->data[cori]) ){
            corVals->data[cori] = *(REAL8 *)item->value;
            break;
          }
        }
      }
      else{
        XLALPrintError("Error... no prior specified!\n");
        XLAL_ERROR_REAL8( XLAL_EFUNC );
      }
    }
  }
  
  /* if there are values for which the priors are defined by a correlation
     coefficient matrix then get add the prior from that */
  if ( corlist ){
    gsl_matrix *cor = NULL;
    gsl_vector_view vals;
    gsl_vector *vm = gsl_vector_alloc( corVals->length );
    UINT4 idx = 0;
    REAL8 ptmp = 0;
    
    if ( LALInferenceCheckVariable( runState->priorArgs, "matrix_inverse" ) ){
      cor = *(gsl_matrix **)LALInferenceGetVariable( runState->priorArgs,
                                                     "matrix_inverse" );
    }
    else{
      LALInferenceGetCorrelatedPrior( runState->priorArgs, corlist->data[0],
                                      &cor, &idx );
    
      /* check for positive definiteness */
      if( !LALInferenceCheckPositiveDefinite( cor, cor->size1 ) ){
        XLALPrintError("Error... matrix is not positive definite!\n");
        XLAL_ERROR_REAL8(XLAL_EFUNC);
      }
   
      /* gsl_linalg_cholesky_invert is not supported in GSL versions < 1.9, so
         until this requirement is changed just use the LU decomposition method
         of calculating the matrix inverse. */
      /* XLAL_CALLGSL( gsl_linalg_cholesky_decomp( cor ) );
      XLAL_CALLGSL( gsl_linalg_cholesky_invert( cor ) );
      */
      gsl_permutation *p = gsl_permutation_alloc ( cor->size1 );
      gsl_matrix *invcor = gsl_matrix_alloc( cor->size1, cor->size2 );
      INT4 s;
      
      XLAL_CALLGSL( gsl_linalg_LU_decomp( cor, p, &s ) );
      XLAL_CALLGSL( gsl_linalg_LU_invert( cor, p, invcor ) );
      XLAL_CALLGSL( gsl_matrix_memcpy( cor, invcor ) );
      gsl_matrix_free( invcor );
      gsl_permutation_free( p );
      
      LALInferenceAddVariable( runState->priorArgs, "matrix_inverse",
                               &cor, LALINFERENCE_gslMatrix_t,
                               LALINFERENCE_PARAM_FIXED );
      
    }

    /* get the log prior (this only works properly if the parameter values have 
       been prescaled so as to be from a Gaussian of zero mean and unit
       variance, which should be the case in this code) */
    vals = gsl_vector_view_array( corVals->data, corVals->length );

    XLAL_CALLGSL( gsl_blas_dgemv(CblasNoTrans, 1., cor, &vals.vector, 0., vm) );
    XLAL_CALLGSL( gsl_blas_ddot(&vals.vector, vm, &ptmp) );

    /* divide by the 2 in the denominator of the Gaussian */
    ptmp /= 2.;
    
    prior -= ptmp;
  }
  
  return prior;
}

/*--------------- END OF LIKELIHOOD AND PRIOR FUNCTIONS ----------------------*/
