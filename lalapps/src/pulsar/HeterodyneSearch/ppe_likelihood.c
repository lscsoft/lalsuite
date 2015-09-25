/**
 * \file
 * \ingroup lalapps_pulsar_HeterodyneSearch
 * \author Matthew Pitkin, John Veitch, Colin Gill
 *
 * \brief Pulsar likelihood and prior functions for use in parameter estimation
 * codes for targeted pulsar searches.
 */

#include "ppe_likelihood.h"

#ifndef _OPENMP
#define omp ignore
#endif


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
 * this for each detector will give the full joint likelihood. See \cite DupuisWoan2005 for a
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
                             LALInferenceModel *get_model ){
  REAL8 loglike = 0.; /* the log likelihood */
  UINT4 i = 0, ifo = 0, nonGR = 0, roq = 0;
  INT4 gaussianLike = 0;
  LALInferenceIFOModel *ifomodeltemp = get_model->ifo;
  LALInferenceIFOData *tempdata = data;

  /* copy model parameters to data parameters */
  LALInferenceCopyVariables( vars, get_model->params );

  /* get pulsar model (this will internally loop over all data frequency streams and IFOs) */
  get_model->templt( get_model );

  if ( LALInferenceCheckVariable( ifomodeltemp->params, "nonGR" ) ){ nonGR = 1; }
  if ( LALInferenceCheckVariable( ifomodeltemp->params, "roq" ) ){ roq = 1; }

  while ( tempdata ){
    /* check if using a Gaussian likelihood - default is the students-t */
    if ( LALInferenceCheckVariable( ifomodeltemp->params, "gaussianLikelihood" ) ){  gaussianLike = 1; }

    get_model->ifo_loglikelihoods[ifo] = 0.0;

    REAL8Vector *sumDat = NULL;
    UINT4Vector *chunkLengths = NULL;
    UINT4 chunkMin, j = 0, count = 0, cl = 0;
    REAL8 chunkLength = 0., logliketmp = 0., chiSquare = 0.;

    sumDat = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumData" );
    chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( ifomodeltemp->params, "chunkLength" );
    chunkMin = *(INT4*)LALInferenceGetVariable( ifomodeltemp->params, "chunkMin" );

    if ( !roq ){ /* not using reduced order quadrature to calculate likelihood */
      UINT4 length = 0;
      REAL8 sumModel = 0., sumDataModel = 0.;
      COMPLEX16 B = 0., M = 0., Mp = 0., Mc = 0., vari = 0.;

      REAL8Vector *sumP = NULL, *sumC = NULL, *sumX = NULL, *sumY = NULL, *sumB = NULL, *sumL = NULL;
      REAL8Vector *sumPC = NULL, *sumPX = NULL, *sumPY = NULL, *sumPB = NULL, *sumPL = NULL;
      REAL8Vector *sumCX = NULL, *sumCY = NULL, *sumCB = NULL, *sumCL = NULL;
      REAL8Vector *sumXY = NULL, *sumXB = NULL, *sumXL = NULL;
      REAL8Vector *sumYB = NULL, *sumYL = NULL;
      REAL8Vector *sumBL = NULL;

      COMPLEX16Vector *sumDataP = NULL, *sumDataC = NULL, *sumDataX = NULL, *sumDataY = NULL, *sumDataB = NULL, *sumDataL = NULL;
      INT4 varyphase = 0;

      if ( LALInferenceCheckVariable( ifomodeltemp->params, "varyphase" ) ){ varyphase = 1; }

      length = tempdata->compTimeData->data->length;

      /* get pre-summed antenna pattern values if required */
      if ( !varyphase ){
        sumP = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumP" );
        sumC = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumC" );
        sumPC = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumPC" );
        sumDataP = *(COMPLEX16Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumDataP" );
        sumDataC = *(COMPLEX16Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumDataC" );

        /* get non-GR components */
        if ( nonGR ){
          sumX = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumX" );
          sumY = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumY" );
          sumB = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumB" );
          sumL = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumL" );

          sumPX = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumPX" );
          sumPY = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumPY" );
          sumPB = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumPB" );
          sumPL = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumPL" );
          sumCX = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumCX" );
          sumCY = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumCY" );
          sumCB = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumCB" );
          sumCL = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumCL" );
          sumXY = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumXY" );
          sumXB = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumXB" );
          sumXL = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumXL" );
          sumYB = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumYB" );
          sumYL = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumYL" );
          sumBL = *(REAL8Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumBL" );

          sumDataX = *(COMPLEX16Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumDataX" );
          sumDataY = *(COMPLEX16Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumDataY" );
          sumDataB = *(COMPLEX16Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumDataB" );
          sumDataL = *(COMPLEX16Vector **)LALInferenceGetVariable( ifomodeltemp->params, "sumDataL" );
        }
      }

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

        if ( varyphase ){ /* not using pre-summed values */
          #pragma omp parallel for default(shared) private(B, M) reduction(+:sumModel,sumDataModel)
          for( j = i ; j < cl ; j++ ){
            B = tempdata->compTimeData->data->data[j];
            M = ifomodeltemp->compTimeSignal->data->data[j];

            vari = 1.;
            if ( gaussianLike ) { vari = tempdata->varTimeData->data->data[j]; }

            /* sum over the model */
            sumModel += (creal(M)*creal(M) + cimag(M)*cimag(M))/vari;

            /* sum over that data and model */
            sumDataModel += (creal(B)*creal(M) + cimag(B)*cimag(M))/vari;
          }
        }
        else{ /* using pre-summed data */
          Mp = ifomodeltemp->compTimeSignal->data->data[0];
          Mc = ifomodeltemp->compTimeSignal->data->data[1];

          sumModel += sumP->data[count]*(creal(Mp)*creal(Mp) + cimag(Mp)*cimag(Mp)) +
                      sumC->data[count]*(creal(Mc)*creal(Mc) + cimag(Mc)*cimag(Mc)) +
                      2.*sumPC->data[count]*(creal(Mp)*creal(Mc) + cimag(Mp)*cimag(Mc));

          sumDataModel += creal(sumDataP->data[count])*creal(Mp) + cimag(sumDataP->data[count])*cimag(Mp) +
                          creal(sumDataC->data[count])*creal(Mc) + cimag(sumDataC->data[count])*cimag(Mc);

          if ( nonGR ){
            COMPLEX16 Mx = 0., My = 0., Mb = 0., Ml = 0.;

            Mx = ifomodeltemp->compTimeSignal->data->data[2];
            My = ifomodeltemp->compTimeSignal->data->data[3];
            Mb = ifomodeltemp->compTimeSignal->data->data[4];
            Ml = ifomodeltemp->compTimeSignal->data->data[5];

            sumModel += sumX->data[count]*(creal(Mx)*creal(Mx) + cimag(Mx)*cimag(Mx)) +
                        sumY->data[count]*(creal(My)*creal(My) + cimag(My)*cimag(My)) +
                        sumB->data[count]*(creal(Mb)*creal(Mb) + cimag(Mb)*cimag(Mb)) +
                        sumL->data[count]*(creal(Ml)*creal(Ml) + cimag(Ml)*cimag(Ml)) +
                        2.*(sumPX->data[count]*(creal(Mp)*creal(Mx) + cimag(Mp)*cimag(Mx)) +
                        sumPY->data[count]*(creal(Mp)*creal(My) + cimag(Mp)*cimag(My)) +
                        sumPB->data[count]*(creal(Mp)*creal(Mb) + cimag(Mp)*cimag(Mb)) +
                        sumPL->data[count]*(creal(Mp)*creal(Ml) + cimag(Mp)*cimag(Ml)) +
                        sumCX->data[count]*(creal(Mc)*creal(Mx) + cimag(Mc)*cimag(Mx)) +
                        sumCY->data[count]*(creal(Mc)*creal(My) + cimag(Mc)*cimag(My)) +
                        sumCB->data[count]*(creal(Mc)*creal(Mb) + cimag(Mc)*cimag(Mb)) +
                        sumCL->data[count]*(creal(Mc)*creal(Ml) + cimag(Mc)*cimag(Ml)) +
                        sumXY->data[count]*(creal(Mx)*creal(My) + cimag(Mx)*cimag(My)) +
                        sumXB->data[count]*(creal(Mx)*creal(Mb) + cimag(Mx)*cimag(Mb)) +
                        sumXL->data[count]*(creal(Mx)*creal(Ml) + cimag(Mx)*cimag(Ml)) +
                        sumYB->data[count]*(creal(My)*creal(Mb) + cimag(My)*cimag(Mb)) +
                        sumYL->data[count]*(creal(My)*creal(Ml) + cimag(My)*cimag(Ml)) +
                        sumBL->data[count]*(creal(Mb)*creal(Ml) + cimag(Mb)*cimag(Ml)));

            sumDataModel += creal(sumDataP->data[count])*creal(Mp) + cimag(sumDataP->data[count])*cimag(Mp) +
                            creal(sumDataC->data[count])*creal(Mc) + cimag(sumDataC->data[count])*cimag(Mc) +
                            creal(sumDataX->data[count])*creal(Mx) + cimag(sumDataX->data[count])*cimag(Mx) +
                            creal(sumDataY->data[count])*creal(My) + cimag(sumDataY->data[count])*cimag(My) +
                            creal(sumDataB->data[count])*creal(Mb) + cimag(sumDataB->data[count])*cimag(Mb) +
                            creal(sumDataL->data[count])*creal(Ml) + cimag(sumDataL->data[count])*cimag(Ml);
          }
        }

        chiSquare = sumDat->data[count] - 2.*sumDataModel + sumModel;

        if ( !gaussianLike ){ /* using Students-t likelihood */
          logliketmp += (gsl_sf_lnfact(chunkLength-1) - LAL_LN2 - chunkLength*LAL_LNPI);
          logliketmp -= chunkLength*log(chiSquare);
        }
        else{ /* using Gaussian likelihood */
          logliketmp -= 0.5*chiSquare;
        }

        count++;
      }
    }
    else{ /* using ROQ to get likelihoods */
      UINT4Vector *numbases = *(UINT4Vector **)LALInferenceGetVariable( ifomodeltemp->params, "numBases" );

      size_t vstart = 0, mstart = 0;

      gsl_vector_complex_view dmview, mmview;
      XLAL_CALLGSL( dmview = gsl_matrix_complex_row(tempdata->roq->weights, 0) );
      XLAL_CALLGSL( mmview = gsl_matrix_complex_row(tempdata->roq->mmweights, 0) );

      /* loop over chunks */
      for ( i = 0; i < chunkLengths->length; i++ ){
        chunkLength = (REAL8)chunkLengths->data[i];

        /* get the weights for this chunk */
        gsl_vector_complex_view dmweights, mmweightsub;
        gsl_matrix_complex_view mmweights;

        dmweights = gsl_vector_complex_subvector(&dmview.vector, vstart, numbases->data[i]);
        mmweightsub = gsl_vector_complex_subvector(&mmview.vector, mstart, numbases->data[i]*numbases->data[i]);
        mmweights = gsl_matrix_complex_view_vector(&mmweightsub.vector, numbases->data[i], numbases->data[i]);

        /* get data/model term */
        gsl_vector_complex_view cmodel, cmodelsub;
        cmodel = gsl_vector_complex_view_array((double*)ifomodeltemp->compTimeSignal->data->data, ifomodeltemp->compTimeSignal->data->length);
        cmodelsub = gsl_vector_complex_subvector(&cmodel.vector, vstart, numbases->data[i]);
        COMPLEX16 dm = LALInferenceROQCOMPLEX16DataDotModel(&dmweights.vector, &cmodelsub.vector);
        COMPLEX16 mm = LALInferenceROQCOMPLEX16ModelDotModel(&mmweights.matrix, &cmodelsub.vector);

        chiSquare = sumDat->data[count] - 2.*creal(dm) + creal(mm);

        if ( !gaussianLike ){ /* using Students-t likelihood */
          logliketmp += (gsl_sf_lnfact(chunkLength-1) - LAL_LN2 - chunkLength*LAL_LNPI);
          logliketmp -= chunkLength*log(chiSquare);
        }
        else{ /* using Gaussian likelihood */
          logliketmp -= 0.5*chiSquare;
        }

        /* shift start indices */
        vstart += numbases->data[i];
        mstart += numbases->data[i]*numbases->data[i];
      }
    }

    /* add normalisation constant for the Gaussian likelihoods */
    if ( gaussianLike ){
      REAL8 logGaussianNorm = *(REAL8 *)LALInferenceGetVariable( ifomodeltemp->params,  "logGaussianNorm" );
      logliketmp += logGaussianNorm;
    }

    get_model->ifo_loglikelihoods[ifo] = logliketmp;

    loglike += logliketmp;
    tempdata->likeli_counter += 1;
    tempdata = tempdata->next;
    ifomodeltemp = ifomodeltemp->next;
    ifo++;
  }

  return loglike;
}


/**
 * \brief Calculate the natural logarithm of the evidence that the data consists of only Gaussian noise
 * The function will calculate the natural logarithm of the evidence that the data (from one or more detectors) consists
 * of stationary segments/chunks described by a Gaussian with zero mean and unknown variance.
 *
 * The evidence is obtained from the joint likelihood given in \c pulsar_log_likelihood with the model term \f$y\f$ set
 * to zero.
 *
 * \param runState [in] The algorithm run state
 *
 * \return The natural logarithm of the noise only evidence
 */
REAL8 noise_only_likelihood( LALInferenceRunState *runState ){
  LALInferenceIFOData *data = runState->data;
  LALInferenceIFOModel *ifo = runState->model->ifo;

  REAL8 logL = 0.0;
  UINT4 i = 0;
  INT4 k = 0;

  REAL8Vector *freqFactors = NULL;
  FILE *fp = NULL;
  CHAR *Znoisefile = NULL;
  ProcessParamsTable *ppt;
  ProcessParamsTable *commandLine = runState->commandLine;
  INT4 gaussianLike = 0;
  /*-----------------------------*/
  /*get the outfile name*/
  ppt = LALInferenceGetProcParamVal( commandLine, "--outfile" );

  freqFactors = *(REAL8Vector **)LALInferenceGetVariable( ifo->params, "freqfactors" );

  /* open the file to output noise evidence (or null signal evidence) for each individual data stream */
  /* set the Znoise filename to the outfile name with "_Znoise" appended */
  Znoisefile = XLALStringDuplicate( ppt->value );
  Znoisefile = XLALStringAppend( Znoisefile, "_Znoise" );

  if( (fp = fopen(Znoisefile, "w")) == NULL ){
    fprintf(stderr, "Error... cannot open output Znoise file!\n");
    exit(0);
  }

  /* check whether using Gaussian or students-t likelihood - using Gaussian if --gaussian
   * has been given as a command line argument. */
  if ( LALInferenceGetProcParamVal( commandLine, "--gaussian-like" ) ){ gaussianLike = 1; }

  /*calculate the evidence */
  while ( ifo ){
    UINT4Vector *chunkLengths = NULL;
    REAL8Vector *sumDat = NULL;

    REAL8 chunkLength = 0.;
    REAL8 logLtmp = 0.;

    chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( ifo->params,  "chunkLength" );
    sumDat = *(REAL8Vector **)LALInferenceGetVariable( ifo->params, "sumData" );
    /*Sum the logL over the datachunks*/
    if ( !gaussianLike ){
      for (i=0; i<chunkLengths->length; i++){
        chunkLength = (REAL8)chunkLengths->data[i];
        logLtmp += (gsl_sf_lnfact(chunkLength-1) - LAL_LN2 - chunkLength*LAL_LNPI);
        logLtmp -= chunkLength * log(sumDat->data[i]);
      }
    }
    else{
      /* for Gaussian likelihood there's only one chunk */
      REAL8 logGaussianNorm = *(REAL8 *)LALInferenceGetVariable( ifo->params,  "logGaussianNorm" );
      logLtmp += logGaussianNorm;
      logLtmp -= 0.5*sumDat->data[0];
    }

    data->nullloglikelihood = logLtmp;

    logL += logLtmp;

    /* output the noise evidence for each data stream for each detector  */
    fprintf(fp, "%s\t%.3lf\t%.16le\n", data->name, freqFactors->data[k], logLtmp);

    k += 1; /* advance counter now, as freqfactors array index starts at zero.*/

    /* reset k, freqfactor counter once all datastreamns for a detector are done */
    if( k >= (INT4)freqFactors->length ) { k = 0; }

    data = data->next;
    ifo = ifo->next;
  }

  fclose(fp);

  return logL;
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
 * \param model UNDOCUMENTED
 *
 * \return The natural logarithm of the prior value for a set of parameters
 */
REAL8 priorFunction( LALInferenceRunState *runState, LALInferenceVariables *params, UNUSED LALInferenceModel *model ){
  LALInferenceIFOModel *ifo = runState->model->ifo;
  (void)runState;
  LALInferenceVariableItem *item = params->head;
  REAL8 prior = 0, value = 0.;

  REAL8Vector *corVals = NULL;
  UINT4 cori = 0;

  /* check that parameters are within their prior ranges */
  if( !in_range( runState->priorArgs, params ) ) { return -DBL_MAX; }

  /* if a k-d tree prior exists ONLY use that */
  if( LALInferenceCheckVariable( runState->priorArgs, "kDTreePrior" ) &&
      LALInferenceCheckVariable( runState->priorArgs, "kDTreePriorTemplate" ) ){
    /* get tree */
    LALInferenceKDTree *tree = *(LALInferenceKDTree **)LALInferenceGetVariable(runState->priorArgs, "kDTreePrior");

    /* get parameter template */
    LALInferenceVariables *template =*(LALInferenceVariables **)LALInferenceGetVariable(runState->priorArgs, "kDTreePriorTemplate");

    /* number of points in a prior cell - i.e. controls how fine or coarse the prior looks (default to 8) */
    UINT4 Ncell = 8;

    if( LALInferenceCheckVariable( runState->priorArgs, "kDTreePriorNcell" ) ){
      Ncell = *(UINT4 *)LALInferenceGetVariable( runState->priorArgs, "kDTreePriorNcell" );
    }

    if ( tree->npts == 0 ) {
      XLALPrintError("%s: no points in prior k-d tree.\n", __func__ );
      XLAL_ERROR_REAL8( XLAL_EFUNC );
    }

    REAL8 *pt = XLALCalloc(tree->dim, sizeof(REAL8));

    /* Get the coordinates of the current point */
    LALInferenceKDVariablesToREAL8(params, pt, template);

    /* find cell of current point */
    LALInferenceKDTree *currentCell = LALInferenceKDFindCell(tree, pt, Ncell);

    /* get log probability of current point - taken from the function LALInferenceKDLogProposalRatio() in LALInference.c */
    REAL8 logVolume = LALInferenceKDLogCellEigenVolume(currentCell);
    REAL8 logCellFactor = log((REAL8)currentCell->npts / (REAL8)tree->npts);

    /* probability is proportional to the inverse of the cell volume */
    prior = -(logVolume + logCellFactor);

    return prior;
  }

  /* if some correlated priors exist allocate corVals */
  if ( corlist ) { corVals = XLALCreateREAL8Vector( corlist->length ); }

  /* I31 and I21 values required to check/set I31 >= I21 */
  REAL8 I31 = -INFINITY, I21 = -INFINITY;

  /* C21 and C22 values required to check/set that they are either both positive or both
   * negative for the case of a biaxial star */
  REAL8 C21 = -INFINITY, C22 = -INFINITY;

  for(; item; item = item->next ){
    /* get scale factor */
    if( item->vary == LALINFERENCE_PARAM_FIXED || item->vary == LALINFERENCE_PARAM_OUTPUT ){ continue; }

    if( item->vary == LALINFERENCE_PARAM_LINEAR || item->vary == LALINFERENCE_PARAM_CIRCULAR ){
      /* Check for a gaussian (note that this is also set if using a correlation coefficient, so
       * also check that is not being used) */
      if ( LALInferenceCheckGaussianPrior(runState->priorArgs, item->name) && !LALInferenceCheckCorrelatedPrior(runState->priorArgs, item->name) ){
        REAL8 mu = 0., sigma = 0.;
        LALInferenceGetGaussianPrior(runState->priorArgs, item->name, &mu, &sigma);

        value = (*(REAL8 *)item->value);
        prior -= 0.5*((*(REAL8 *)item->value - mu)*(*(REAL8 *)item->value - mu))/sigma;

        if ( !strcmp(item->name, "I21") ){ I21 = value; }
        if ( !strcmp(item->name, "I31") ){ I31 = value; }
        if ( !strcmp(item->name, "C21") ){ C21 = value; }
        if ( !strcmp(item->name, "C22") ){ C22 = value; }
      }
      /* check for a flat prior */
      else if( LALInferenceCheckMinMaxPrior(runState->priorArgs, item->name) ){
        value = (*(REAL8 *)item->value);

        /* check if either using theta or iota rather than their cosines */
        if ( !strcmp(item->name, "IOTA") || !strcmp(item->name, "THETA") ){
          prior += theta_prior( value );
        }
        else { /* note: we don't need to update the prior as it is a constant */
          if ( !strcmp(item->name, "I21") ){ I21 = value; }
          if ( !strcmp(item->name, "I31") ){ I31 = value; }
          if ( !strcmp(item->name, "C21") ){ C21 = value; }
          if ( !strcmp(item->name, "C22") ){ C22 = value; }
        }
      }
      else if( LALInferenceCheckCorrelatedPrior(runState->priorArgs, item->name) && corlist ){
        /* set item in correct position given the order of the correlation matrix given by corlist */
        REAL8 mu = 0., sigma = 0.;
        LALInferenceGetGaussianPrior(runState->priorArgs, item->name, &mu, &sigma);

        for( cori = 0; cori < corlist->length; cori++ ){
          if( !strcmp(item->name, corlist->data[cori]) ){
            /* scale to a zero mean and unit variance Gaussian */
            corVals->data[cori] = (*(REAL8 *)item->value - mu)/sigma;
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

  /* if a biaxial star check that C21 and C22 are the same sign */
  if ( LALInferenceCheckVariable( ifo->params, "biaxial" ) ){
    /* in case one parameter is fixed check that */
    if ( C21 == -INFINITY ){ C21 = LALInferenceGetREAL8Variable( runState->currentParams, "C21" ); }
    if ( C22 == -INFINITY ){ C22 = LALInferenceGetREAL8Variable( runState->currentParams, "C22" ); }
    if ( (C21 < 0. && C22 > 0.) || (C21 > 0. && C22 < 0.) ) { return -DBL_MAX; } /* if same sign this will be positive */
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
        XLAL_ERROR_REAL8( XLAL_EFUNC );
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

      LALInferenceAddVariable( runState->priorArgs, "matrix_inverse", &cor, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED );
    }

    /* get the log prior (this only works properly if the parameter values have been prescaled so as to be from a
     * Gaussian of zero mean and unit variance, which happens on line 501) */
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
 * \brief Prior for angle that is equivalent to an latitude value to give a uniform prior on a sphere
 *
 * If you have two angles that define spherical polar coordinates and you want these to have a prior
 * that is uniform over the sphere then you can instead work with the cosine of the latitude-like angle
 * and have this to be uniform between -1 and 1. However, if you want to work in the actual angle then
 * you need to set the correct prior which will be \f$p(\theta) \propto \sin{\theta)\f$.
 *
 * \param theta [in] The angle in radians
 *
 * \return The log prior as defined above.
 */
REAL8 theta_prior( REAL8 theta ){
  return log(sin(theta));
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
        LALInferenceIFOModel *ifo = runState->model->ifo;

        /* with the scaled parameters the k-D tree ranges will be between 0 and 1 */
        lownew[cnt] = 0.;
        highnew[cnt] = 1.;

        INT4 ii = 0;

        /* now change the current param to reflect new ranges */
        REAL8 var = *(REAL8 *)LALInferenceGetVariable( runState->currentParams, samp->name );
        while( ifo ){
          /* rescale current parameter value */
          if ( ii == 0 ){
            ii++;
            LALInferenceSetVariable( runState->currentParams, samp->name, &var );
          }

          ifo = ifo->next;
        }

        /* reset (or add) priorArgs */
        LALInferenceAddMinMaxPrior( runState->priorArgs, samp->name, &(lownew[cnt]), &(highnew[cnt]),
                                    LALINFERENCE_REAL8_t );
      }
      else{
        lownew[cnt] = low[cnt];
        highnew[cnt] = high[cnt];
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

    while( samp ){
      if ( samp->vary != LALINFERENCE_PARAM_FIXED && samp->vary != LALINFERENCE_PARAM_OUTPUT ) {
        REAL8 val = *(REAL8 *)samp->value;
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
