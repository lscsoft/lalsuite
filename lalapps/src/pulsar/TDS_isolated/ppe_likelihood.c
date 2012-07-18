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
                             LALInferenceTemplateFunction *get_model){
  REAL8 loglike = 0.; /* the log likelihood */
  UINT4 i = 0;
  CHAR *modeltype = NULL; /*need to check model type in this function*/
  
  modeltype = *(CHAR**)LALInferenceGetVariable( data->dataParams, "modeltype" );
  
  LALInferenceIFOData *datatemp1 = data, *datatemp2 = data, *datatemp3 = data;
 
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
  
  while ( datatemp3 ){
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
    
    sumDat = *(REAL8Vector **)LALInferenceGetVariable( datatemp3->dataParams, 
                                                       "sumData" );
    chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( datatemp3->dataParams,
                                                             "chunkLength" );
    chunkMin = *(INT4*)LALInferenceGetVariable( datatemp3->dataParams, "chunkMin" );
  
    length = datatemp3->compTimeData->data->length;
  
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
        B.re = datatemp3->compTimeData->data->data[j].re;
        B.im = datatemp3->compTimeData->data->data[j].im;

        M.re = datatemp3->compModelData->data->data[j].re;
        M.im = datatemp3->compModelData->data->data[j].im;
        
        /* sum over the model */
        sumModel += M.re*M.re + M.im*M.im;
        
        /* sum over that data and model */
        sumDataModel += B.re*M.re + B.im*M.im;
        /*fprintf(bugtest,"B.re= %e, B.im= %e, M.re: %e, M.im: %e\n",B.re, B.im,M.re, M.im);*/
      }
 
      chiSquare = sumDat->data[count];
      chiSquare -= 2.*sumDataModel;
      chiSquare += sumModel;
      
      logliketmp -= chunkLength*log(chiSquare) + LAL_LN2 * (chunkLength-1.) + gsl_sf_lnfact(chunkLength);
      
      count++;
    }
    loglike += logliketmp;
    datatemp3 = datatemp3->next;
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
  
  const CHAR *fn = __func__;
  
  /* check that parameters are with their prior ranges */
  if( !in_range( runState->priorArgs, params ) ) return -INFINITY;
  //LALInferenceCyclicReflectiveBound( params, runState->priorArgs );
  
  /* if a k-d tree prior exists ONLY use that */
  if( LALInferenceCheckVariable( runState->priorArgs, "kDTreePrior" ) &&
      LALInferenceCheckVariable( runState->priorArgs, "kDTreePriorTemplate" ) ){
    /* get tree */
    LALInferenceKDTree *tree =
      *(LALInferenceKDTree **)LALInferenceGetVariable(runState->priorArgs,
                                                      "kDTreePrior");
    
    /* get parameter template */
    LALInferenceVariables *template = 
      *(LALInferenceVariables **)LALInferenceGetVariable(runState->priorArgs,
                                                         "kDTreePriorTemplate");
    
    UINT4 Ncell = 16; /* number of points in a prior cell - i.e. controls
                         how fine or coarse the prior looks (default to 16) */ 
      
    if( LALInferenceCheckVariable( runState->priorArgs, "kDTreePriorNcell" ) ){
      Ncell = *(UINT4 *)LALInferenceGetVariable( runState->priorArgs,
                                                 "kDTreePriorNcell" );
    }
    
    if ( tree->npts == 0 ) {
      XLALPrintError("%s: no points in prior k-d tree.\n", fn);
      XLAL_ERROR_REAL8( XLAL_EFUNC );
    }

    REAL8 *pt = XLALCalloc(tree->dim, sizeof(REAL8));

    /* Get the coordinates of the current point - points in the tree are
       already scaled, so there's no need to rescale the current point */
    LALInferenceKDVariablesToREAL8(params, pt, template);
    
    /* find cell of current point */
    LALInferenceKDTree *currentCell = LALInferenceKDFindCell(tree, pt, Ncell);
    
    /* get log probability of current point - taken from the function
       LALInferenceKDLogProposalRatio() in LALInference.c */
    //REAL8 logVolume = LALInferenceKDLogCellEigenVolume(currentCell);
    REAL8 logVolume = LALInferenceKDLogCellVolume(currentCell);
    REAL8 logCellFactor = log((REAL8)currentCell->npts / (REAL8)tree->npts);
    
    prior = logVolume + logCellFactor;
    
    return prior;
  }
 
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
        prior -= log( (max - min) * scale );
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

/** \brief Convert an array of nested samples to posterior samples
 *
 * This function generates an array of posterior samples from the nested samples
 * by drawing points from the nested samples weighted by their prior weighting.
 * This assumes that the nested samples are in the array in ascending
 * likelihood order (which should be the case for the output of the
 * \c LALInferenceNestedSampler() function.
 * 
 * Within the input runstate->algorthimParams there needs to be: an array of
 * LALInferenceVariables called "nestedsamples" containing nested samples to be
 * converted in to the posterior; a value of "Nsamp" giving the number of
 * nested samples; and, a value of "numberlive" giving the number live points
 * used to generate the posterior samples.
 * 
 * The posterior samples will be output in runstate->algorthimParams as an
 * array of LALInferenceVariables called "posteriorsamples", and the number of
 * these in a value called "Nposterior".
 *
 * \param runState [in]
 */
void ns_to_posterior( LALInferenceRunState *runState ){
  UINT4 i = 0, count = 0, k = 0;
  UINT4Vector *Nsamp = NULL, *Nlive = NULL;
  REAL8 maxlogw = -INFINITY; /* maximum log weight */ 
  
  LALInferenceVariables **psamples = NULL;
  LALInferenceVariables ***nsamples = 
    *(LALInferenceVariables ****)LALInferenceGetVariable(
    runState->algorithmParams, "nestedsamples" );
    
  /* get number of samples and live points */ 
  Nsamp = *(UINT4Vector **)LALInferenceGetVariable( runState->algorithmParams,
                                                    "Nsamps" );
  Nlive = *(UINT4Vector **)LALInferenceGetVariable( runState->algorithmParams,
                                                    "numberlive" );

  if ( Nsamp->length != Nlive->length ){
    XLALPrintError("%s: Number of nested sample arrays not equal to number of \
live point for each array!", __func__);
    XLAL_ERROR_VOID( XLAL_EBADLEN );
  }
  
  for( k = 0; k < Nsamp->length; k++ ){
    /* vector of prior weights */
    REAL8Vector *logw = XLALCreateREAL8Vector( Nsamp->data[k] );
  
    maxlogw = -INFINITY;
    
    /* fill in sample weights */
    for ( i = 0; i < Nsamp->data[k]; i++ ){
      REAL8 logL = *(REAL8 *)LALInferenceGetVariable( nsamples[k][i], "logL" );
    
      if( i < Nsamp->data[k]-Nlive->data[k] ) logw->data[i] = (REAL8)(i+1);
      else logw->data[i] = (REAL8)(Nsamp->data[k] - Nlive->data[k]);
    
      logw->data[i] = -(logw->data[i]/(REAL8)Nlive->data[k]) + logL;
    
      if ( logw->data[i] > maxlogw ) maxlogw = logw->data[i];
    }

    /* get posterior samples */
    for ( i = 0; i < Nsamp->data[k]; i++ ){
      logw->data[i] -= maxlogw; /* normalise weights */
    
      /* if log weight is greater than a uniform random number then accept as
         a posterior sample */
      if ( logw->data[i] > log( gsl_rng_uniform( runState->GSLrandom ) ) ){
        psamples = XLALRealloc( psamples, 
                                (count+1)*sizeof(LALInferenceVariables*) );
        psamples[count] = XLALCalloc( 1, sizeof(LALInferenceVariables) );
        LALInferenceCopyVariables( nsamples[k][i], psamples[count] );
        count++;
      }
    }
    
    XLALDestroyREAL8Vector( logw );
  }
  
  LALInferenceAddVariable( runState->algorithmParams, "posteriorsamples",
                           &psamples, LALINFERENCE_void_ptr_t,
                           LALINFERENCE_PARAM_FIXED );
  LALInferenceAddVariable( runState->algorithmParams, "Nposterior",
                           &count, LALINFERENCE_UINT4_t,
                           LALINFERENCE_PARAM_FIXED );
}

/** \brief Create a k-d tree from prior samples
 *
 * This function creates a k-d tree from prior samples for use as a prior
 * distribution in the algorithm. The points in the output tree are scaled to
 * the prior ranges specified on the command line.
 * 
 * \param runState [in] A pointer to the LALInferenceRunState
 */
void create_kdtree_prior( LALInferenceRunState *runState ){
  LALInferenceKDTree *priortree = NULL;
  LALInferenceVariables **posterior = NULL; /* use these samples as prior */
  UINT4 nsamp = 0, i = 0, cnt = 0;
  
  LALInferenceVariableItem *samp = NULL; /* a single sample */
  REAL8 *low = NULL, *high = NULL; /* upper and lower bounds of tree */
  REAL8 *pt = NULL;
  size_t ndim = 0;
  
  LALInferenceVariables *template =
    XLALCalloc(1,sizeof(LALInferenceVariables));
  
  const CHAR *fn = __func__;
  
  /* get posterior samples to use as prior */
  if ( LALInferenceCheckVariable( runState->algorithmParams, 
    "posteriorsamples" ) ){
    posterior = *(LALInferenceVariables ***)LALInferenceGetVariable(
      runState->algorithmParams, "posteriorsamples" );
  }
  else{
    XLALPrintError("%s: No posterior samples set to use as prior.\n", fn);
    XLAL_ERROR_VOID( XLAL_EFUNC );
  }
  
  /* get the number of posterior samples */
  if ( LALInferenceCheckVariable( runState->algorithmParams, 
    "Nposterior" ) ){
    nsamp = *(UINT4 *)LALInferenceGetVariable( runState->algorithmParams,
                                               "Nposterior" );
  }
  else{
    XLALPrintError("%s: Number of posterior samples not set.\n", fn);
    XLAL_ERROR_VOID( XLAL_EFUNC );
  }
  
  /* get the upper and lower bounds for each variable parameter i.e. we won't
     be adding log likelihood, or log prior values */
  samp = posterior[0]->head;
  while ( samp ){
    if ( samp->vary != LALINFERENCE_PARAM_FIXED &&
         samp->vary != LALINFERENCE_PARAM_OUTPUT ) {
      if( LALInferenceCheckMinMaxPrior( runState->priorArgs,
                                        samp->name ) ){
        cnt++;
         
        low = XLALRealloc(low, sizeof(REAL8)*cnt);
        high = XLALRealloc(high, sizeof(REAL8)*cnt);
      
        LALInferenceGetMinMaxPrior( runState->priorArgs, samp->name,
                                    &(low[cnt-1]), &(high[cnt-1]) );
      }
      else if( LALInferenceCheckGaussianPrior( runState->priorArgs,
                                               samp->name ) ){
        /* REAL8 mn, stddiv; */
        REAL8 postlow, posthigh, difflh;
        
        cnt++;
        
        low = XLALRealloc(low, sizeof(REAL8)*cnt);
        high = XLALRealloc(high, sizeof(REAL8)*cnt);

        /* LALInferenceGetGaussianPrior( runState->priorArgs, currentItem->name,
                                      &mn, &stddiv ); */

        /* find the maximum and minimum posterior point values */
        REAL8 maxvaltmp = -INFINITY, minvaltmp = INFINITY;
  
        for ( UINT4 k = 0; k < nsamp; k++ ){
          REAL8 val = *(REAL8 *)LALInferenceGetVariable( posterior[k], 
                                                         samp->name );
    
          if ( val < minvaltmp ) minvaltmp = val;
          if ( val > maxvaltmp ) maxvaltmp = val;
        }
  
        posthigh = maxvaltmp;
        postlow = minvaltmp;
        
        difflh = posthigh - postlow;

        /* to add a bit of room at either side add on half the difference */
        low[cnt-1] = postlow - difflh/2.;
        high[cnt-1] = posthigh + difflh/2.;
      }
    }
    
    samp = samp->next;
  }

  ndim = (size_t)cnt;
  pt = XLALMalloc(cnt*sizeof(REAL8));

  /* set up tree */
  priortree = LALInferenceKDEmpty( low, high, ndim );
  
  /* get template */
  LALInferenceCopyVariables( posterior[0], template );
  
  /* add points to tree */
  for( i = 0; i < nsamp; i++ ){
    samp = posterior[i]->head;
    
    /* rescale sample */
    CHAR scalePar[VARNAME_MAX] = "";
    CHAR scaleMinPar[VARNAME_MAX] = "";
    REAL8 scale = 0., scaleMin = 0.;
    
    while( samp ){
      if ( samp->vary != LALINFERENCE_PARAM_FIXED &&
         samp->vary != LALINFERENCE_PARAM_OUTPUT ) {
        
        sprintf(scalePar, "%s_scale", samp->name);
        scale = *(REAL8 *)LALInferenceGetVariable( runState->data->dataParams,
                                                   scalePar );
    
        sprintf(scaleMinPar, "%s_scale_min", samp->name);
        scaleMin = *(REAL8 *)LALInferenceGetVariable(
          runState->data->dataParams, scaleMinPar );
        
        REAL8 val = *(REAL8 *)samp->value;
        val = (val - scaleMin)/scale;
       
        LALInferenceSetVariable( posterior[i], samp->name, &val );
      }
      
      samp = samp->next;
    }

    LALInferenceKDVariablesToREAL8( posterior[i], pt, template );
   
    LALInferenceKDAddPoint( priortree, pt );
  }
  
  /* add tree */
  LALInferenceAddVariable( runState->priorArgs, "kDTreePrior", &priortree,
                           LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED );

  LALInferenceAddVariable( runState->priorArgs, "kDTreePriorTemplate",
                           &template, LALINFERENCE_void_ptr_t,
                           LALINFERENCE_PARAM_FIXED );

  XLALFree( high );
  XLALFree( low );
  XLALFree( pt );
}


/** \brief Check that any parameters with minimum and maximum ranges are within
 * that range
 *
 * This function performs any cylcic/reflective transform and then makes sure
 * that all parameters in \c params, that have a defined minimum and maximum 
 * value, are within their allowed prior ranges.
 * 
 * \param priors [in] A pointer to the prior args LALInferenceVariables
 * \param params [in] The current parameters
 * 
 * \return 0 if out of range and 1 if in range
 */
UINT4 in_range( LALInferenceVariables *priors, LALInferenceVariables *params ){
  LALInferenceVariableItem *item = params->head;
  REAL8 min, max;
  
  /* loop over variables */
  for(; item; item = item->next ){
    if( item->vary == LALINFERENCE_PARAM_FIXED || 
      item->vary == LALINFERENCE_PARAM_OUTPUT ){ continue; }
    
    if( LALInferenceCheckMinMaxPrior( priors, item->name ) ){
      LALInferenceGetMinMaxPrior( priors, item->name, &min, &max );
      
      /* For cyclic boundaries, mod out by range. */
      if( item->vary == LALINFERENCE_PARAM_CIRCULAR ) {
        REAL8 val = *(REAL8 *)item->value;
        REAL8 delta = max - min;

        if (val > max) {
          REAL8 offset = val - min;

          *(REAL8 *)item->value = min + fmod(offset, delta);
        } 
        else {
          REAL8 offset = max - val;

          *(REAL8 *)item->value = max - fmod(offset, delta);
        }
        
        continue;
      }
      
      if( (*(REAL8 *) item->value) < min || (*(REAL8 *)item->value) > max )
        return 0;
    }
  }
  
  return 1;
}
  

/*--------------- END OF LIKELIHOOD AND PRIOR FUNCTIONS ----------------------*/
