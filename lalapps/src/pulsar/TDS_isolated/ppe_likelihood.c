/**
 * \file
 * \ingroup pulsarApps
 * \author Matthew Pitkin, John Veitch, Colin Gill
 *
 * \brief Pulsar likelihood and prior functions for use in parameter estimation
 * codes for targeted pulsar searches.
 */

#include "ppe_likelihood.h"

/******************************************************************************/
/*                     LIKELIHOOD AND PRIOR FUNCTIONS                         */
/******************************************************************************/
/**
 * \brief The log likelihood function
 *
 * This function calculates natural logarithm of the likelihood of a signal model (specified by a given set of
 * parameters) given the data from a set of detectors.
 *
 * The likelihood is the joint likelihood of chunks of data over which the noise is assumed stationary and Gaussian. For
 * each chunk a Gaussian likelihood for the noise and data has been marginalised over the unknown noise standard
 * deviation using a Jeffreys prior on the standard deviation. Given the data consisting of independent real and
 * imaginary parts this gives a Students-t distribution for each chunk (of length \f$m\f$) with \f$m/2\f$ degrees of
 * freedom:
 * \f[
 * p(\mathbf{\theta}|\mathbf{B}) = \prod_{j=1}^M \frac{(m_j-1)!}{2\pi^{m_j}}
 * \left( \sum_{k=k_0}^{k_0+(m_j-1)} |B_k - y(\mathbf{\theta})_k|^2
 * \right)^{-m_j},
 * \f]
 * where \f$\mathbf{B}\f$ is a vector of the complex data, \f$y(\mathbf{\theta})\f$ is the model for a set of parameters
 * \f$\mathbf{\theta}\f$, \f$M\f$ is the total number of independent data chunks with lengths \f$m_j\f$ and \f$k_0 =
 * \sum_{i=1}^j 1 + m_{i-1}\f$ (with \f$m_0 = 0\f$) is the index of the first data point in each chunk. The product of
 * this for each detector will give the full joint likelihood. In the calculation here the unnecessary proportionality
 * factors are left out (this would effect the actual value of the marginal likelihood/evidence, but since we are only
 * interested in evidence ratios/Bayes factors these factors would cancel out anyway. See \cite DupuisWoan2005 for a
 * more detailed description.
 *
 * In this function data in chunks smaller than a certain minimum length \c chunkMin are ignored.
 *
 * \param vars [in] The parameter values
 * \param data [in] The detector data and initial signal phase template
 * \param get_model [in] The signal model structure
 *
 * \return The natural logarithm of the likelihood function
 */
REAL8 pulsar_log_likelihood( LALInferenceVariables *vars, LALInferenceIFOData *data,
                             LALInferenceModel *get_model){
  REAL8 loglike = 0.; /* the log likelihood */
  UINT4 i = 0;
  REAL8Vector *freqFactors = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams, "freqfactors" );

  LALInferenceIFOData *datatemp1 = data, *datatemp2 = data, *datatemp3 = data;

  /* copy model parameters to data parameters */
  while( datatemp1 ){
    LALInferenceCopyVariables( vars, get_model->params );
    datatemp1 = datatemp1->next;
  }

  /* get pulsar model */
  while( datatemp2 ){
    get_model->templt( datatemp2 );

    for( i = 0; i < freqFactors->length; i++ ) { datatemp2 = datatemp2->next; }
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

    sumDat = *(REAL8Vector **)LALInferenceGetVariable( datatemp3->dataParams, "sumData" );
    chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( datatemp3->dataParams, "chunkLength" );
    chunkMin = *(INT4*)LALInferenceGetVariable( datatemp3->dataParams, "chunkMin" );

    length = datatemp3->compTimeData->data->length;

    for( i = 0 ; i < length ; i += chunkLength ){
      chunkLength = (REAL8)chunkLengths->data[count];

      /* skip section of data if its length is less than the minimum allowed chunk length */
      if( chunkLength < chunkMin ){
        count++;
        continue;
      }

      sumModel = 0.;
      sumDataModel = 0.;

      cl = i + (INT4)chunkLength;

      for( j = i ; j < cl ; j++ ){
        B = datatemp3->compTimeData->data->data[j];

        M = get_model->compSignal->data->data[j];

        /* sum over the model */
        sumModel += creal(M)*creal(M) + cimag(M)*cimag(M);

        /* sum over that data and model */
        sumDataModel += creal(B)*creal(M) + cimag(B)*cimag(M);
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

/**
 * \brief The prior function
 *
 * This function calculates the natural logarithm of the prior for a set of parameters. If the prior on a particular
 * parameter is uniform over a given range then \f$p(\theta) = 1/(\theta_{\rm max} - \theta_{\rm min})\f$. If the
 * prior is Gaussian then the probability of that value given the mean and standard deviation of the Gaussian is
 * calculated.
 *
 * \param runState [in] A pointer to the LALInferenceRunState
 * \param params [in] The set of parameter values
 *
 * \return The natural logarithm of the prior value for a set of parameters
 */
REAL8 priorFunction( LALInferenceRunState *runState, LALInferenceVariables *params, UNUSED LALInferenceModel *model ){
  LALInferenceIFOData *data = runState->data;
  (void)runState;
  LALInferenceVariableItem *item = params->head;
  REAL8 min, max, mu, sigma, prior = 0, value = 0.;

  REAL8Vector *corVals = NULL;
  UINT4 cori = 0;

  const CHAR *fn = __func__;

  /* check that parameters are within their prior ranges */
  if( !in_range( runState->priorArgs, params ) ) { return -INFINITY; }

  /* if a k-d tree prior exists ONLY use that */
  if( LALInferenceCheckVariable( runState->priorArgs, "kDTreePrior" ) &&
      LALInferenceCheckVariable( runState->priorArgs, "kDTreePriorTemplate" ) ){
    /* get tree */
    LALInferenceKDTree *tree = *(LALInferenceKDTree **)LALInferenceGetVariable(runState->priorArgs, "kDTreePrior");

    /* get parameter template */
    LALInferenceVariables *template =*(LALInferenceVariables **)LALInferenceGetVariable(runState->priorArgs,
                                                                                        "kDTreePriorTemplate");

    /* number of points in a prior cell - i.e. controls how fine or coarse the prior looks (default to 8) */
    UINT4 Ncell = 8;

    if( LALInferenceCheckVariable( runState->priorArgs, "kDTreePriorNcell" ) ){
      Ncell = *(UINT4 *)LALInferenceGetVariable( runState->priorArgs, "kDTreePriorNcell" );
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
    REAL8 logVolume = LALInferenceKDLogCellEigenVolume(currentCell);
    // REAL8 logVolume = LALInferenceKDLogCellVolume(currentCell);
    REAL8 logCellFactor = log((REAL8)currentCell->npts / (REAL8)tree->npts);

    /* probability is proportional to the inverse of the cell volume */
    prior = -(logVolume + logCellFactor);

    return prior;
  }

  /* if some correlated priors exist allocate corVals */
  if ( corlist ) { corVals = XLALCreateREAL8Vector( corlist->length ); }

  /* I31 and I21 values required to check/set I31 >= I21 */
  REAL8 I31 = -INFINITY;
  REAL8 I21 = -INFINITY;

  for(; item; item = item->next ){
    /* get scale factor */
    CHAR scalePar[VARNAME_MAX] = "";
    CHAR scaleMinPar[VARNAME_MAX] = "";
    REAL8 scale = 0., scaleMin = 0.;

    if( item->vary == LALINFERENCE_PARAM_FIXED || item->vary == LALINFERENCE_PARAM_OUTPUT ){ continue; }

    sprintf(scalePar, "%s_scale", item->name);
    scale = *(REAL8 *)LALInferenceGetVariable( data->dataParams, scalePar );

    sprintf(scaleMinPar, "%s_scale_min", item->name);
    scaleMin = *(REAL8 *)LALInferenceGetVariable( data->dataParams, scaleMinPar );

    if( item->vary == LALINFERENCE_PARAM_LINEAR || item->vary == LALINFERENCE_PARAM_CIRCULAR ){
      /* Check for a gaussian */
      if ( LALInferenceCheckGaussianPrior(runState->priorArgs, item->name) ){
        LALInferenceGetGaussianPrior( runState->priorArgs, item->name, &mu, &sigma );

        value = (*(REAL8 *)item->value) * scale + scaleMin;
        mu += scaleMin;
        sigma *= scale;
        prior -= log(sqrt(2.*LAL_PI)*sigma);
        prior -= (value - mu)*(value - mu) / (2.*sigma*sigma);

        if ( !strcmp(item->name, "I21") ){ I21 = value; }
        if ( !strcmp(item->name, "I31") ){ I31 = value; }
      }
      /* check for a flat prior */
      else if( LALInferenceCheckMinMaxPrior(runState->priorArgs, item->name) ){
        LALInferenceGetMinMaxPrior( runState->priorArgs, item->name, &min, &max );
        prior -= log( (max - min) * scale );

        if ( !strcmp(item->name, "I21") ){ I21 = (*(REAL8 *)item->value) * scale + scaleMin; }
        if ( !strcmp(item->name, "I31") ){ I31 = (*(REAL8 *)item->value) * scale + scaleMin; }
      }
      else if( LALInferenceCheckCorrelatedPrior(runState->priorArgs, item->name) && corlist ){
        /* set item in correct position given the order of the correlation matrix given by corlist */
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

    /* check if the I21 and I31 value have been set, and if so set the prior I31 >= I21 */
    if ( I21 != -INFINITY && I31 != -INFINITY ){
      if ( I31 < I21 ) { return -DBL_MAX; }
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
      cor = *(gsl_matrix **)LALInferenceGetVariable( runState->priorArgs, "matrix_inverse" );
    }
    else{
      LALInferenceGetCorrelatedPrior( runState->priorArgs, corlist->data[0], &cor, &idx );

      /* check for positive definiteness */
      if( !LALInferenceCheckPositiveDefinite( cor, cor->size1 ) ){
        XLALPrintError("Error... matrix is not positive definite!\n");
        XLAL_ERROR_REAL8(XLAL_EFUNC);
      }

      /* gsl_linalg_cholesky_invert is not supported in GSL versions < 1.9, so until this requirement is changed just
       * use the LU decomposition method of calculating the matrix inverse. */
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

      LALInferenceAddVariable( runState->priorArgs, "matrix_inverse", &cor, LALINFERENCE_gslMatrix_t,
                               LALINFERENCE_PARAM_FIXED );

    }

    /* get the log prior (this only works properly if the parameter values have been prescaled so as to be from a
     * Gaussian of zero mean and unit variance, which should be the case in this code) */
    vals = gsl_vector_view_array( corVals->data, corVals->length );

    XLAL_CALLGSL( gsl_blas_dgemv(CblasNoTrans, 1., cor, &vals.vector, 0., vm) );
    XLAL_CALLGSL( gsl_blas_ddot(&vals.vector, vm, &ptmp) );

    /* divide by the 2 in the denominator of the Gaussian */
    ptmp /= 2.;

    prior -= ptmp;
  }

  return prior;
}

/**
 * \brief Convert an array of nested samples to posterior samples
 *
 * This function generates an array of posterior samples from the nested samples by drawing points from the nested
 * samples weighted by their prior weighting. This assumes that the nested samples are in the array in ascending
 * likelihood order (which should be the case for the output of the \c LALInferenceNestedSampler() function. The
 * posterior sample generation is based on the method used in lalapps/src/inspiral/posterior/nest2pos.py
 *
 * Within the input runstate->algorthimParams there needs to be: an array of LALInferenceVariables called
 * "nestedsamples" containing nested samples to be converted in to the posterior; a value of "Nsamp" giving the number
 * of nested samples; and, a value of "numberlive" giving the number live points used to generate the posterior samples.
 *
 * The posterior samples will be output in runstate->algorthimParams as an array of LALInferenceVariables called
 * "posteriorsamples", and the number of these in a value called "Nposterior".
 *
 * \param runState [in]
 */
void ns_to_posterior( LALInferenceRunState *runState ){
  UINT4 i = 0, count = 0, k = 0;
  UINT4Vector *Nsamp = NULL, *Nlive = NULL;

  UINT4 Npost = 0; /* no. of posterior samples required from each NS run */
  ProcessParamsTable *ppt = NULL;

  LALInferenceVariables **psamples = NULL;
  LALInferenceVariables ***nsamples = *(LALInferenceVariables ****)LALInferenceGetVariable( runState->algorithmParams,
                                                                                            "nestedsamples" );

  /* get number of samples and live points */
  Nsamp = *(UINT4Vector **)LALInferenceGetVariable( runState->algorithmParams, "Nsamps" );
  Nlive = *(UINT4Vector **)LALInferenceGetVariable( runState->algorithmParams, "numberlive" );
  /* get (approximate) number of posterior samples to generate from each nested sample file */
  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--Npost" );
  if ( ppt != NULL ) { Npost = atoi(ppt->value); }
  else { Npost = 1000; } /* default to 1000 */

  if ( Nsamp->length != Nlive->length ){
    XLALPrintError("%s: Number of nested sample arrays not equal to number of live point for each array!", __func__);
    XLAL_ERROR_VOID( XLAL_EBADLEN );
  }

  REAL8Vector *log_evs = XLALCreateREAL8Vector( Nsamp->length );
  REAL8Vector **log_ws = XLALCalloc(Nsamp->length, (sizeof(REAL8Vector*)));
  REAL8 log_tot_ev = -INFINITY;

  for( k = 0; k < Nsamp->length; k++ ){
    /* vector of prior weights */
    log_ws[k] = XLALCreateREAL8Vector( Nsamp->data[k] );

    REAL8 log_vol_factor = log(1. - (1./(REAL8)Nlive->data[k]));
    REAL8 log_dvol = -1./(REAL8)Nlive->data[k];
    REAL8 log_vol = 0., log_ev = -INFINITY;
    REAL8 avg_log_like_end = -INFINITY;

    /* fill in sample weights */
    for ( i = 0; i < Nsamp->data[k]; i++ ){
      REAL8 logL = *(REAL8 *)LALInferenceGetVariable( nsamples[k][i], "logL" );

      if( i < Nsamp->data[k]-Nlive->data[k] ){
        log_ws[k]->data[i] = logL + log_vol + log_dvol;
        log_ev = LOGPLUS(log_ev, log_ws[k]->data[i]);
        log_vol += log_vol_factor;
      }
      else{
        avg_log_like_end = LOGPLUS(avg_log_like_end, logL);
        log_ws[k]->data[i] = log_vol + logL;
      }
    }

    avg_log_like_end -= log((REAL8)Nlive->data[k]);
    log_ev = LOGPLUS(log_ev, avg_log_like_end + log_vol);

    log_evs->data[k] = log_ev;
    log_tot_ev = LOGPLUS( log_tot_ev, log_ev );
  }

  for( k = 0; k < Nsamp->length; k++ ){
    /* round Npost based on the evidence for each data set */
    UINT4 Ns = (UINT4)ROUND((REAL8)Npost*exp(log_evs->data[k]-log_tot_ev));

    /* generate cumulative sum of weigths */
    REAL8Vector *log_cumsums = XLALCreateREAL8Vector( Nsamp->data[k] + 1 );
    log_cumsums->data[0] = -INFINITY;

    /* get posterior samples */
    for ( i = 0; i < Nsamp->data[k]; i++ ){
      log_ws[k]->data[i] -= log_evs->data[k]; /* normalise weights */

      log_cumsums->data[i+1] = LOGPLUS( log_cumsums->data[i], log_ws[k]->data[i] );
    }

    for( UINT4 j = 0; j < Ns; j++ ){
      REAL8 us = log( gsl_rng_uniform( runState->GSLrandom ) );

      /* draw the samples */
      for( i = 1; i < Nsamp->data[k]+1; i++ )
        if ( log_cumsums->data[i-1] < us && us <= log_cumsums->data[i] ) { break; }

      /* add the posterior sample */
      psamples = XLALRealloc( psamples, (count+1)*sizeof(LALInferenceVariables*) );
      psamples[count] = XLALCalloc( 1, sizeof(LALInferenceVariables) );
      LALInferenceCopyVariables( nsamples[k][i-1], psamples[count] );
      count++;
    }
  }

  /* free weights and evidence */
  for ( k = 0; k < Nsamp->length; k++ ) { XLALDestroyREAL8Vector( log_ws[k] ); }
  XLALFree( log_ws );
  XLALDestroyREAL8Vector( log_evs );

  LALInferenceAddVariable( runState->algorithmParams, "posteriorsamples", &psamples, LALINFERENCE_void_ptr_t,
                           LALINFERENCE_PARAM_FIXED );
  LALInferenceAddVariable( runState->algorithmParams, "Nposterior", &count, LALINFERENCE_UINT4_t,
                           LALINFERENCE_PARAM_FIXED );
}

/**
 * \brief Create a k-d tree from prior samples
 *
 * This function creates a k-d tree from prior samples for use as a prior distribution in the algorithm. The ranges of
 * the tree dimensions (parameters) are calculated from the maximum and minimum values of the parameter. If the the
 * lower and upper values are closer than half of there overall range to those specified in the prior file then the
 * prior file values are used. Otherwise the upper/lower value + half the overall range is used. In this case the prior
 * ranges set in priorArgs, and scale factors, will need to be replaced. [NOTE: This case will mean that the whole prior
 * range specified may not be searched, however the prior samples would suggest that there is very little probability of
 * information existing in the excluded areas]. The reason for not just using the full prior ranges is that if the prior
 * samples are tightly peaked in a small area of the prior space the k-D tree resolution is poor and can cause an
 * unwanted broadening of the prior e.g. if the h0 samples are peaked at zero with a distribution width of ~1e-24, but
 * the prior range is set from 0 to 1e-22 then new samples drawn from a k-D tree produced with the full range in h0
 * yields a broader distribution than required.
 *
 * [NOTE: it is not obvious how this would affect evidence comparisons!]
 *
 * \param runState [in] A pointer to the LALInferenceRunState
 */
void create_kdtree_prior( LALInferenceRunState *runState ){
  LALInferenceKDTree *priortree = NULL;
  LALInferenceVariables **posterior = NULL; /* use these samples as prior */
  UINT4 nsamp = 0, i = 0, cnt = 0;

  LALInferenceVariableItem *samp = NULL; /* a single sample */
  REAL8 *low = NULL, *high = NULL;
  REAL8 *lownew = NULL, *highnew = NULL; /* upper and lower bounds of tree */
  REAL8 *pt = NULL;
  size_t ndim = 0;

  LALInferenceVariables *template = XLALCalloc(1,sizeof(LALInferenceVariables));

  const CHAR *fn = __func__;

  /* get posterior samples to use as prior */
  if ( LALInferenceCheckVariable( runState->algorithmParams, "posteriorsamples" ) ){
    posterior = *(LALInferenceVariables ***)LALInferenceGetVariable( runState->algorithmParams, "posteriorsamples" );
  }
  else{
    XLALPrintError("%s: No posterior samples set to use as prior.\n", fn);
    XLAL_ERROR_VOID( XLAL_EFUNC );
  }

  /* get the number of posterior samples */
  if ( LALInferenceCheckVariable( runState->algorithmParams, "Nposterior" ) ){
    nsamp = *(UINT4 *)LALInferenceGetVariable( runState->algorithmParams, "Nposterior" );
  }
  else{
    XLALPrintError("%s: Number of posterior samples not set.\n", fn);
    XLAL_ERROR_VOID( XLAL_EFUNC );
  }

  /* get the upper and lower bounds for each variable parameter i.e. we won't be adding log likelihood, or log prior
   * values */
  samp = posterior[0]->head;
  while ( samp ){
    UINT4 change = 0; /* set it a range has changed */

    if ( samp->vary != LALINFERENCE_PARAM_FIXED && samp->vary != LALINFERENCE_PARAM_OUTPUT ) {
      /* get the minimum and maximum prior sample and the range */
      REAL8 maxvaltmp = -INFINITY, minvaltmp = INFINITY;
      REAL8 posthigh = 0., postlow = 0., difflh = 0.;
      REAL8 lowr = 0., highr = 0.;

      /* get the original scale factors */
      CHAR scalePar[VARNAME_MAX] = "";
      CHAR scaleMinPar[VARNAME_MAX] = "";
      REAL8 scale = 0., scaleMin = 0.;

      sprintf(scalePar, "%s_scale", samp->name);
      scale = *(REAL8 *)LALInferenceGetVariable( runState->data->dataParams, scalePar );
      sprintf(scaleMinPar, "%s_scale_min", samp->name);
      scaleMin = *(REAL8 *)LALInferenceGetVariable( runState->data->dataParams, scaleMinPar );

      for ( UINT4 k = 0; k < nsamp; k++ ){
        REAL8 val = *(REAL8 *)LALInferenceGetVariable( posterior[k], samp->name );

        if ( val < minvaltmp ) { minvaltmp = val; }
        if ( val > maxvaltmp ) { maxvaltmp = val; }
      }

      /* max and min of the prior samples */
      posthigh = maxvaltmp;
      postlow = minvaltmp;

      difflh = posthigh - postlow; /* range of the samples */

      /* dynamically allocate memory for the k-D tree ranges */
      low = XLALRealloc(low, sizeof(REAL8)*(cnt+1));
      high = XLALRealloc(high, sizeof(REAL8)*(cnt+1));
      lownew = XLALRealloc(lownew, sizeof(REAL8)*(cnt+1));
      highnew = XLALRealloc(highnew, sizeof(REAL8)*(cnt+1));

      /* if there's a minmax prior check the ranges */
      if( LALInferenceCheckMinMaxPrior( runState->priorArgs, samp->name ) ){
        LALInferenceGetMinMaxPrior( runState->priorArgs, samp->name, &lowr, &highr );

        highr = scaleMin + highr*scale;
        lowr = scaleMin;

        /* check how far min and max samples are from the prior ranges */
        if( posthigh > highr || posthigh < lowr ){
          fprintf(stderr, "Error... prior samples are out of range!\n");
          exit(1);
        }
        else if( posthigh + difflh/2. > highr ) { high[cnt] = highr; } /* just stick with prior range */
        else{
          high[cnt] = posthigh + difflh/2.; /* use max from samples+diff/2 */
          change++; /* we have changed the range */
        }

        if( postlow > highr || postlow < lowr ){
          fprintf(stderr, "Error... prior samples are out of range!\n");
          exit(1);
        }
        else if( postlow - difflh/2. < lowr ) { low[cnt] = lowr; } /* just stick with prior range */
        else{
          low[cnt] = postlow - difflh/2.; /* use min from samples-diff/2 */
          change++; /* we have changed the range */
        }
      }
      else if( LALInferenceCheckGaussianPrior( runState->priorArgs, samp->name ) ){
        /* to add a bit of room at either side add on half the difference */
        low[cnt] = postlow - difflh/2.;
        high[cnt] = posthigh + difflh/2.;

        /* remove the Gaussian prior and replace with min/max range */
        LALInferenceRemoveGaussianPrior( runState->priorArgs, samp->name );

        change++; /* the range will change */
      }
      else{
        fprintf(stderr, "Error... Prior type not specified!\n");
        exit(1);
      }

      /* change the prior ranges, the scale factors and the current params */
      if( change != 0 ){
        LALInferenceIFOData *data = runState->data;
        fprintf(stderr, "Here\n");
        REAL8 newscale = high[cnt] - low[cnt], newscaleMin = low[cnt];

        /* with the scaled parameters the k-D tree ranges will be between 0 and 1 */
        lownew[cnt] = 0.;
        highnew[cnt] = 1.;

        INT4 ii = 0;

        /* now change the current param to reflect new ranges */
        REAL8 var = *(REAL8 *)LALInferenceGetVariable( runState->currentParams, samp->name );
        while( data ){
          /* rescale current parameter value */
          if ( ii == 0 ){
            ii++;

            var = scaleMin + var*scale;
            var = (var-newscaleMin)/newscale;

            LALInferenceSetVariable( runState->currentParams, samp->name, &var );
          }

          /* change the scale factors */
          LALInferenceRemoveVariable( data->dataParams, scalePar );
          LALInferenceRemoveVariable( data->dataParams, scaleMinPar );

          LALInferenceAddVariable( data->dataParams, scalePar, &newscale, LALINFERENCE_REAL8_t,
                                   LALINFERENCE_PARAM_FIXED );
          LALInferenceAddVariable( data->dataParams, scaleMinPar, &newscaleMin, LALINFERENCE_REAL8_t,
                                   LALINFERENCE_PARAM_FIXED );

          data = data->next;
        }

        /* reset (or add) priorArgs */
        LALInferenceAddMinMaxPrior( runState->priorArgs, samp->name, &(lownew[cnt]), &(highnew[cnt]),
                                    LALINFERENCE_REAL8_t );
      }
      else{
        lownew[cnt] = (low[cnt] - scaleMin)/scale;
        highnew[cnt] = (high[cnt] - scaleMin)/scale;
      }

      cnt++;
    }

    samp = samp->next;
  }

  ndim = (size_t)cnt;
  pt = XLALMalloc(cnt*sizeof(REAL8));

  /* set up tree */
  priortree = LALInferenceKDEmpty( lownew, highnew, ndim );

  /* get template */
  LALInferenceCopyVariables( posterior[0], template );

  /* add points to tree */
  for( i = 0; i < nsamp; i++ ){
    samp = posterior[i]->head;

    cnt = 0;

    /* rescale samples with new scales */
    while( samp ){
      if ( samp->vary != LALINFERENCE_PARAM_FIXED && samp->vary != LALINFERENCE_PARAM_OUTPUT ) {
        REAL8 newscale = high[cnt] - low[cnt], newscaleMin = low[cnt];
        REAL8 val = *(REAL8 *)samp->value;
        val = (val - newscaleMin)/newscale;

        LALInferenceSetVariable( posterior[i], samp->name, &val );

        cnt++;
      }

      samp = samp->next;
    }

    LALInferenceKDVariablesToREAL8( posterior[i], pt, template );
    LALInferenceKDAddPoint( priortree, pt );
  }

  /* add tree */
  LALInferenceAddVariable( runState->priorArgs, "kDTreePrior", &priortree, LALINFERENCE_void_ptr_t,
                           LALINFERENCE_PARAM_FIXED );

  LALInferenceAddVariable( runState->priorArgs, "kDTreePriorTemplate", &template, LALINFERENCE_void_ptr_t,
                           LALINFERENCE_PARAM_FIXED );

  XLALFree( high );
  XLALFree( low );
  XLALFree( pt );
}


/**
 * \brief Check that any parameters with minimum and maximum ranges are within that range
 *
 * This function performs any cylcic/reflective transform and then makes sure that all parameters in \c params, that
 * have a defined minimum and maximum value, are within their allowed prior ranges.
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
    if( item->vary == LALINFERENCE_PARAM_FIXED || item->vary == LALINFERENCE_PARAM_OUTPUT ){ continue; }

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
