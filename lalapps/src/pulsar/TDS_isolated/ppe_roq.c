/**
 * \file
 * \ingroup pulsarApps
 * \author Matthew Pitkin, John Veitch, Colin Gill
 *
 * \brief Reduced order quadrature generation functions for use in parameter estimation
 * codes for targeted pulsar searches.
 */

#include "ppe_roq.h"

/******************************************************************************/
/*                     REDUCED ORDER QUADRATURE FUNCTIONS                     */
/******************************************************************************/

/**
 * \brief Generate Chebyshev-Gauss-Lobatto nodes in frequency
 *
 * @param[in] fmin The minimum frequency
 * @param[in] fmax The maximum frequency
 * @param[in] nnodes The number of nodes
 *
 * @return An array with the node freqeuncy values
 */
REAL8 *chebyshev_gauss_lobatto_nodes( REAL8 fmin, REAL8 fmax, UINT4 nnodes ){
  UINT4 i = 0;
  REAL8 *fnodes = NULL;
  REAL8 df, fplus, n;

  XLAL_CHECK_NULL( fnodes = XLALMalloc(nnodes*sizeof(REAL8)), XLAL_EFUNC, "Couldn't allocate memory for frequency nodes." );

  df = (fmax-fmin)/2.;
  fplus = (fmax+fmin)/2.;
  n = (REAL8)nnodes - 1.;

  for ( i=0; i<nnodes; i++ ){ fnodes[i] = -cos(LAL_PI*(REAL8)i/n)*df + fplus; }

  return fnodes;
}


/**
 * \brief Generate an orthonormal basis set of waveforms from a training set
 *
 * This function will use the prior ranges on the parameters to generate a set of
 * training waveforms. From these training waveforms it will generate a set of
 * orthonormal basis functions, with the number generated controlled by a stopping
 * tolerance criterion. In general the values of the parameters for the training set
 * will be drawn randomly from across the prior ranges. However, ff frequency is one
 * of the parameters being used then values will be placed the the Chebyshev-Gauss-Lobatto
 * nodes.
 *
 * This basis set will be generated seperately for the real and imaginary parts of
 * the model.
 *
 * @param[in] runState The algorithm run state
 *
 * @return A structure holding the real and complex orthonormal bases
 */
void generate_interpolant( LALInferenceRunState *runState ){
  REAL8 tolerance = ROQTOLERANCE;
  UINT4 ntraining = 0, i = 0, j = 0, verbose = 0;
  INT4 gaussianLike = 0;
  ProcessParamsTable *ppt;

  LALInferenceIFOModel *ifo = NULL;
  LALInferenceIFOData *data = runState->data;

  /* timing values */
  struct timeval time1, time2;
  REAL8 tottime;

  if ( LALInferenceCheckVariable( runState->algorithmParams, "timefile" ) ){ gettimeofday(&time1, NULL); }

  /* check whether to use ROQ */
  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--roq" );
  if ( !ppt ){ return; }

  /* get tolerance stopping criterion */
  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--roq-tolerance" );
  if ( ppt ){ tolerance = atof( ppt->value ); }
  XLAL_CHECK_VOID( tolerance < 1. && tolerance > 0., XLAL_EFUNC, "ROQ tolerence (%le) is not within allowed range.", tolerance );

  /* get number of training sets */
  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--ntraining" );
  if ( ppt ){ ntraining = atoi( ppt->value ); }
  else { XLAL_ERROR_VOID( XLAL_EFUNC, "Number of training sets must be specifed if running with --roq" ); }
  XLAL_CHECK_VOID( ntraining > 1, XLAL_EFUNC, "Number of training sets (%d) is too small!", ntraining );

  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--verbose" );
  if ( ppt ){ verbose = 1; }

  /* check if using a Gaussian likelihoodm and therefore need variance weighted ROQ values */
  if ( LALInferenceGetProcParamVal( runState->commandLine, "--gaussian-like" ) ){ gaussianLike = 1; }

  ifo = runState->model->ifo;

  while ( ifo ){
    UINT4Vector *chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( ifo->params, "chunkLength" );
    REAL8 dt = *(REAL8*)LALInferenceGetVariable( ifo->params, "dt" );
    INT4 startidx = 0;

    UINT4 dlen = 0;
    UINT4Vector *nbases = XLALCreateUINT4Vector( chunkLengths->length );

    /* vectors to hold the vector and matrices with the interpolation weights */
    COMPLEX16 *dmweights = NULL;
    COMPLEX16 *mmweights = NULL;
    UINT4 dmlength = 0, mmlength = 0;

    /* create a copy of the model vector */
    LIGOTimeGPS gpstime = {0, 0};
    REAL8Vector *sidDayFrac = *(REAL8Vector**)LALInferenceGetVariable( ifo->params, "siderealDay" );
    REAL8Vector *ssbdelays = NULL, *bsbdelays = NULL;

    LIGOTimeGPSVector *timenodes = NULL;
    REAL8Vector *sidtimenodes = NULL;
    REAL8TimeSeries *timedatanodes = XLALCreateREAL8TimeSeries( "", &gpstime, 0., 1., &lalSecondUnit, 1 );
    REAL8Vector *ssbnodes = NULL, *bsbnodes = NULL;

    data->roq = XLALMalloc(sizeof(LALInferenceROQData));

    if ( LALInferenceCheckVariable( ifo->params, "ssb_delays") ){
      ssbdelays = *(REAL8Vector**)LALInferenceGetVariable( ifo->params, "ssb_delays" );
    }

    if ( LALInferenceCheckVariable( ifo->params, "bsb_delays") ){
      bsbdelays = *(REAL8Vector**)LALInferenceGetVariable( ifo->params, "bsb_delays" );
    }

    /* get chunk */
    for ( i=0; i<chunkLengths->length; i++ ){
      /* a temporary run state containing just the required data for a single detector */
      LALInferenceRunState *tmpRS = XLALMalloc(sizeof(LALInferenceRunState));
      LALInferenceIFOModel *ifotmp = XLALMalloc(sizeof(LALInferenceIFOModel));
      ifotmp->next = NULL;
      tmpRS->model = XLALMalloc(sizeof(LALInferenceModel));
      tmpRS->model->templt = runState->model->templt;
      tmpRS->currentParams = runState->currentParams;
      tmpRS->model->params = XLALCalloc(1, sizeof(LALInferenceVariables));
      tmpRS->model->ifo = ifotmp;
      tmpRS->GSLrandom = runState->GSLrandom;
      tmpRS->priorArgs = runState->priorArgs;

      gsl_matrix_complex *ts = NULL;

      UINT4 tlen = chunkLengths->data[i];

      ifotmp->times = XLALCreateTimestampVector( tlen );
      /* create a new model vector */
      ifotmp->compTimeSignal = XLALCreateCOMPLEX16TimeSeries( "", &gpstime, 0., 1., &lalSecondUnit, tlen );
      ifotmp->timeData = XLALCreateREAL8TimeSeries( "", &gpstime, 0., 1., &lalSecondUnit, tlen );

      ifotmp->params = ifo->params; /* just use a pointer to params in runState */
      ifotmp->ephem = ifo->ephem;
      ifotmp->detector = ifo->detector;
      ifotmp->tdat = ifo->tdat;
      ifotmp->ttype = ifo->ttype;

      /* get chunk times */
      for ( j=0; j<tlen; j++ ){
        ifotmp->times->data[j] = ifo->times->data[startidx+j];
        ifotmp->timeData->data->data[j] = ifo->timeData->data->data[startidx+j];
      }

      /* generate the training set */
      ts = generate_training_set( tmpRS, ntraining, 1 );

      /* generate reduced basis */
      gsl_vector_view weights = gsl_vector_view_array(&dt, 1);
      size_t nbases0 = 0;
      COMPLEX16 *RB = NULL;
      RB = LALInferenceGenerateCOMPLEX16OrthonormalBasis(&weights.vector, tolerance, ts, &nbases0);
      XLAL_CHECK_VOID( RB != NULL, XLAL_EFUNC, "Could not produce basis set");
      nbases->data[i] = nbases0;

      if ( verbose ){
       fprintf(stderr, "Number of reduced bases for ROQ generation is %zu\n", nbases0);
      }

      XLAL_CALLGSL( gsl_matrix_complex_free( ts ) );

      if ( nbases0 > (size_t)tlen-1 ){
        XLAL_ERROR_VOID( XLAL_EFUNC, "Number of bases is longer than date length, so no point using ROQ" );
      }

      /* if required test the basis */
      if ( LALInferenceGetProcParamVal( runState->commandLine, "--test-basis" ) ){
        /* generate another set of waveforms and test that the reduced basis can match them well enough */
        ts = generate_training_set( tmpRS, ntraining, 0 );

        /* check values are within tolerance (use 100 time the defined tolerance for the basis generation) */
        gsl_matrix_complex_view m = gsl_matrix_complex_view_array( (double*)RB, nbases0, tlen );
        if ( LALInferenceTestCOMPLEX16OrthonormalBasis( &weights.vector, tolerance*100., &m.matrix, ts ) != XLAL_SUCCESS ){
          XLAL_ERROR_VOID( XLAL_EFUNC, "Basis does not cover the space to the required tolerance" );
        }

        XLAL_CALLGSL( gsl_matrix_complex_free( ts ) );
      }

      XLALDestroyTimestampVector( ifotmp->times );
      XLALDestroyCOMPLEX16TimeSeries( ifotmp->compTimeSignal );
      XLALDestroyREAL8TimeSeries( ifotmp->timeData );
      LALInferenceClearVariables( tmpRS->model->params );
      XLALFree( tmpRS->model );
      XLALFree( tmpRS );

      /* generate the interpolants (pass the data noise variances for weighting the interpolants in the case of using a Gaussian likelihood) */
      LALInferenceCOMPLEXROQInterpolant *interp = NULL;
      gsl_matrix_complex_view RBview = gsl_matrix_complex_view_array((double*)RB, nbases0, tlen);
      interp = LALInferenceGenerateCOMPLEXROQInterpolant(&RBview.matrix);

      /* free the reduced basis set */
      XLALFree( RB );

      /* get pointer to data and variance */
      REAL8 varone = 1.;
      gsl_vector_view vari;
      gsl_vector_complex_view dataview, datasub;

      dataview = gsl_vector_complex_view_array((double*)data->compTimeData->data->data, data->compTimeData->data->length);
      datasub = gsl_vector_complex_subvector( &dataview.vector, startidx, tlen );

      if ( gaussianLike ){
        gsl_vector_view varview = gsl_vector_view_array(data->varTimeData->data->data, data->varTimeData->data->length);
        vari = gsl_vector_subvector( &varview.vector, startidx, tlen );
      }
      else{ vari = gsl_vector_view_array( &varone, 1 ); } /* using Students-t likelihood, so just set to 1 */

      /* create the data/model and model/model inner product weights */
      gsl_vector_complex *dmw = LALInferenceGenerateCOMPLEX16DataModelWeights(interp->B, &datasub.vector, &vari.vector);
      gsl_matrix_complex *mmw = LALInferenceGenerateCOMPLEXModelModelWeights(interp->B, &vari.vector);

      /* put weights into a vector */
      mmlength += (nbases0*nbases0);
      dmlength += nbases0;

      dmweights = XLALRealloc( dmweights, sizeof(COMPLEX16)*dmlength );
      mmweights = XLALRealloc( mmweights, sizeof(COMPLEX16)*mmlength );

      gsl_vector_complex_view mmptr, dmptr, mmptrsub, dmptrsub;
      gsl_matrix_complex_view cmmptr;

      dmptr = gsl_vector_complex_view_array((double*)dmweights, dmlength);
      dmptrsub = gsl_vector_complex_subvector(&dmptr.vector, dmlength-nbases0, nbases0);

      mmptr = gsl_vector_complex_view_array((double*)mmweights, mmlength);
      mmptrsub = gsl_vector_complex_subvector(&mmptr.vector, mmlength-(nbases0*nbases0), (nbases0*nbases0));

      cmmptr = gsl_matrix_complex_view_vector(&mmptrsub.vector, nbases0, nbases0);

      XLAL_CALLGSL( gsl_vector_complex_memcpy(&dmptrsub.vector, dmw) );
      XLAL_CALLGSL( gsl_matrix_complex_memcpy(&cmmptr.matrix, mmw) );

      /* get times at the interpolation nodes */
      timenodes = XLALResizeTimestampVector( timenodes, dmlength );
      sidtimenodes = XLALResizeREAL8Vector( sidtimenodes, dmlength );
      timedatanodes = XLALResizeREAL8TimeSeries( timedatanodes, 0, dmlength );

      if ( ssbdelays != NULL ){ ssbnodes = XLALResizeREAL8Vector( ssbnodes, dmlength ); }
      if ( bsbdelays != NULL ){ bsbnodes = XLALResizeREAL8Vector( bsbnodes, dmlength ); }

      for ( j=0; j<nbases0; j++ ){
        timenodes->data[dlen+j] = ifo->times->data[startidx+interp->nodes[j]];
        sidtimenodes->data[dlen+j] = sidDayFrac->data[startidx+interp->nodes[j]];
        timedatanodes->data->data[dlen+j] = ifo->timeData->data->data[startidx+interp->nodes[j]];

        if ( ssbdelays != NULL ){ ssbnodes->data[dlen+j] = ssbdelays->data[startidx+interp->nodes[j]]; }
        if ( bsbdelays != NULL ){ bsbnodes->data[dlen+j] = bsbdelays->data[startidx+interp->nodes[j]]; }
      }

      LALInferenceRemoveCOMPLEXROQInterpolant( interp );

      XLAL_CALLGSL( gsl_vector_complex_free( dmw ) );
      XLAL_CALLGSL( gsl_matrix_complex_free( mmw ) );

      startidx += tlen;
      dlen += nbases0;
    }

    /* resize time vectors and model vectors, so they just contain values at the interpolation nodes */
    LIGOTimeGPSVector *timescopy = XLALCreateTimestampVector( ifo->times->length );
    memcpy(timescopy->data, ifo->times->data, sizeof(LIGOTimeGPS)*ifo->times->length);
    LALInferenceAddVariable( ifo->params, "timeStampVectorFull", &timescopy, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED );
    XLALDestroyTimestampVector(ifo->times);
    ifo->times = timenodes;

    /* make a copy of the original sideral time vectors and time stamp vectors - this is needed when
     * calculating the SNR of the maximum likelihood point */
    REAL8Vector *siddaycopy = XLALCreateREAL8Vector( sidDayFrac->length );
    memcpy(siddaycopy->data, sidDayFrac->data, sizeof(REAL8)*sidDayFrac->length);
    LALInferenceRemoveVariable( ifo->params, "siderealDay" );
    LALInferenceAddVariable( ifo->params, "siderealDay", &sidtimenodes, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
    LALInferenceAddVariable( ifo->params, "siderealDayFull", &siddaycopy, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );

    REAL8TimeSeries *timedatacopy = XLALCreateREAL8TimeSeries( "", &gpstime, 0., 1., &lalSecondUnit, ifo->timeData->data->length );
    memcpy(timedatacopy->data->data, ifo->timeData->data->data, sizeof(REAL8)*ifo->timeData->data->length);
    LALInferenceAddVariable( ifo->params, "timeDataFull", &timedatacopy, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED );
    XLALDestroyREAL8TimeSeries( ifo->timeData );
    ifo->timeData = timedatanodes;

    if ( ssbdelays != NULL ){
      REAL8Vector *ssbcopy = XLALCreateREAL8Vector( ssbdelays->length );
      memcpy(ssbcopy->data, ssbdelays->data, sizeof(REAL8)*ssbdelays->length);
      LALInferenceRemoveVariable( ifo->params, "ssb_delays" );
      LALInferenceAddVariable( ifo->params, "ssb_delays", &ssbnodes, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
      LALInferenceAddVariable( ifo->params, "ssb_delays_full", &ssbcopy, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
    }

    if ( bsbdelays != NULL ){
      REAL8Vector *bsbcopy = XLALCreateREAL8Vector( bsbdelays->length );
      memcpy(bsbcopy->data, bsbdelays->data, sizeof(REAL8)*bsbdelays->length);
      LALInferenceRemoveVariable( ifo->params, "bsb_delays" );
      LALInferenceAddVariable( ifo->params, "bsb_delays", &bsbnodes, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
      LALInferenceAddVariable( ifo->params, "bsb_delays_full", &bsbcopy, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
    }

    ifo->compTimeSignal = XLALResizeCOMPLEX16TimeSeries( ifo->compTimeSignal, 0, dmlength );

    /* fill in data/model weights into roq->weights complex matrix */
    gsl_matrix_complex_view dmview;
    XLAL_CALLGSL( dmview = gsl_matrix_complex_view_array((double*)dmweights, 1, dmlength) );
    XLAL_CALLGSL( data->roq->weights = gsl_matrix_complex_alloc(1, dmlength) );
    XLAL_CALLGSL( gsl_matrix_complex_memcpy(data->roq->weights, &dmview.matrix) );
    XLALFree( dmweights );

    gsl_matrix_complex_view mmview;
    XLAL_CALLGSL( mmview = gsl_matrix_complex_view_array((double*)mmweights, 1, mmlength) );
    XLAL_CALLGSL( data->roq->mmweights = gsl_matrix_complex_alloc(1, mmlength) );
    XLAL_CALLGSL( gsl_matrix_complex_memcpy(data->roq->mmweights, &mmview.matrix) );
    XLALFree( mmweights );

    /* add interpolation weights and nodes to a variable in runState->model->ifo->params */
    LALInferenceAddVariable( ifo->params, "numBases", &nbases, LALINFERENCE_UINT4Vector_t, LALINFERENCE_PARAM_FIXED );

    ifo = ifo->next;
    data = data->next;
  }

  if ( LALInferenceCheckVariable( runState->algorithmParams, "timefile" ) ){
    gettimeofday(&time2, NULL);

    FILE *timefile = *(FILE **)LALInferenceGetVariable( runState->algorithmParams, "timefile" );
    UINT4 timenum = *(UINT4 *)LALInferenceGetVariable( runState->algorithmParams, "timenum" );
    tottime = (REAL8)((time2.tv_sec + time2.tv_usec*1.e-6) - (time1.tv_sec + time1.tv_usec*1.e-6));
    fprintf(timefile, "[%d] %s: %.9le secs\n", timenum, __func__, tottime);
    timenum++;
    check_and_add_fixed_variable( runState->algorithmParams, "timenum", &timenum, LALINFERENCE_UINT4_t );
  }
}


/**
 * \brief Generate a training set of waveforms for the basis generation
 *
 * This function will create a set of \a n waveforms randomly placed over the
 * prior parameter space. If \a freqnodes is 1 then if a frequency parameter
 * is required in the training set generation then it will use values calculated
 * at the Chebyshev-Gauss-Lobatto nodes.
 *
 * @param[in] rs A temporary run state
 * @param[in] n The number of training waveforms
 * @param[in] freqnodes A flag for whether to set frequencies at Chebyshev-Gauss-Lobatto nodes or not
 *
 * @return A complex matrix containing an array of training waveforms
 */
gsl_matrix_complex *generate_training_set( LALInferenceRunState *rs, UINT4 n, UINT4 freqnodes ){
  UINT4 j = 0;
  gsl_matrix_complex *ts = gsl_matrix_complex_alloc(n, rs->model->ifo->times->length);
  REAL8 *fnodes = NULL;

  for ( j=0; j<n; j++ ){
    /* choose random variables values and fill in runState->model->params */
    LALInferenceVariableItem *item = rs->currentParams->head;

    /* there's no need to re-scale values, so I can drawn from the scaled values */
    for(; item; item = item->next ){
      REAL8 value;

      if( item->vary == LALINFERENCE_PARAM_FIXED || item->vary == LALINFERENCE_PARAM_OUTPUT ){
        LALInferenceAddVariable( rs->model->params, item->name, item->value, item->type, item->vary );
        continue;
      }

      if( item->vary == LALINFERENCE_PARAM_LINEAR || item->vary == LALINFERENCE_PARAM_CIRCULAR ){
        /* Check for a gaussian (generate a point between +/-5 sigma) */
        if ( LALInferenceCheckGaussianPrior(rs->priorArgs, item->name) ){
          /* as we're using scaled values the mean is zero and sigma is 1 */
          value = -5. + 10.*gsl_rng_uniform(rs->GSLrandom);
        }
        /* check for a flat prior */
        else if( LALInferenceCheckMinMaxPrior(rs->priorArgs, item->name) ){
        /* as we're using scaled values they will be between zero and 1 */

          /* if variable is f0 (frequency) then set the Chebyshev-Gauss-Lobatto nodes */
          if ( ( !strcmp(item->name, "F0") || !strcmp(item->name, "f0") ) && freqnodes ){
            if ( j == 0 ){ fnodes = chebyshev_gauss_lobatto_nodes( 0., 1., n ); }
            value = fnodes[j];
          }
          else{ value = gsl_rng_uniform(rs->GSLrandom); }
          //value = gsl_rng_uniform(rs->GSLrandom);
        }
        else if( LALInferenceCheckCorrelatedPrior(rs->priorArgs, item->name) && corlist ){
          /* for this correlation matrix values are also scaled as for the Gaussian prior */
          value = -5. + 10.*gsl_rng_uniform(rs->GSLrandom);
        }
        else{ XLAL_ERROR_NULL( XLAL_EFUNC, "Error... no prior specified!\n" ); }

        LALInferenceAddVariable( rs->model->params, item->name, &value, item->type, item->vary );
      }
    }

    /* generate model */
    rs->model->templt( rs->model );

    /* place model into an array */
    gsl_vector_complex_view cview;
    cview = gsl_vector_complex_view_array((double*)rs->model->ifo->compTimeSignal->data->data, rs->model->ifo->times->length);
    gsl_matrix_complex_set_row(ts, j, &cview.vector);
  }

  if( fnodes != NULL ){ XLALFree( fnodes ); }

  return ts;
}


/*--------------- END OF ROQ FUNCTIONS ----------------------*/
