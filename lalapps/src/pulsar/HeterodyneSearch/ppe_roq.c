/**
 * \file
 * \ingroup lalapps_pulsar_HeterodyneSearch
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
 * @param[in] freqmin The minimum frequency
 * @param[in] freqmax The maximum frequency
 * @param[in] nnodes The number of nodes
 *
 * @return An array with the node freqeuncy values
 */
REAL8 *chebyshev_gauss_lobatto_nodes( REAL8 freqmin, REAL8 freqmax, UINT4 nnodes ){
  UINT4 i = 0;
  REAL8 *fnodes = NULL;
  REAL8 df, fplus, n;

  XLAL_CHECK_NULL( fnodes = XLALMalloc(nnodes*sizeof(REAL8)), XLAL_EFUNC, "Couldn't allocate memory for frequency nodes." );

  df = (freqmax-freqmin)/2.;
  fplus = (freqmax+freqmin)/2.;
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
  UINT4 ntraining = 0, i = 0, j = 0, verbose = 0, outputroq = 0, inputroq = 0, counter = 0;
  INT4 gaussianLike = 0;
  ProcessParamsTable *ppt;

  /* variables if reading in weights file */
  UINT4 nstreams = 0, nchunks = 0;
  FILE *fpin = NULL, *fpout = NULL;

  LALInferenceIFOModel *ifo = NULL;
  LALInferenceIFOData *data = runState->data;

  /* timing values */
  struct timeval time1, time2;
  REAL8 tottime;

  if ( LALInferenceCheckVariable( runState->algorithmParams, "timefile" ) ){ gettimeofday(&time1, NULL); }

  /* check whether to use ROQ */
  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--roq" );
  if ( !ppt ){ return; }

  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--input-weights" );
  if ( !ppt ){
    /* get tolerance stopping criterion */
    ppt = LALInferenceGetProcParamVal( runState->commandLine, "--roq-tolerance" );
    if ( ppt ){ tolerance = atof( ppt->value ); }
    XLAL_CHECK_VOID( tolerance < 1. && tolerance > 0., XLAL_EFUNC, "ROQ tolerence (%le) is not within allowed range.", tolerance );

    /* get number of training sets */
    ppt = LALInferenceGetProcParamVal( runState->commandLine, "--ntraining" );
    if ( ppt ){ ntraining = atoi( ppt->value ); }
    else { XLAL_ERROR_VOID( XLAL_EFUNC, "Number of training sets must be specifed if running with --roq" ); }
    XLAL_CHECK_VOID( ntraining > 1, XLAL_EFUNC, "Number of training sets (%d) is too small!", ntraining );

    /* check whether to output weights */
    ppt = LALInferenceGetProcParamVal( runState->commandLine, "--output-weights" );
    if ( ppt ){
      CHAR *outputweights = ppt->value;
      XLAL_CHECK_VOID( (fpout = fopen(outputweights, "wb")) != NULL, XLAL_EIO, "Could not open weights file for output." );
      outputroq = 1;

      /* write out the number of data streams */
      nstreams = *(UINT4*)LALInferenceGetVariable( runState->algorithmParams, "numstreams" );
      XLAL_CHECK_VOID( fwrite(&nstreams, sizeof(UINT4), 1, fpout) == 1, XLAL_EIO, "Could not write the number of data steams to file!" );
    }
  }
  else{
    /* read in weights from a file */
    CHAR *inputweights = ppt->value;
    UINT4 nstreamcheck = *(UINT4*)LALInferenceGetVariable( runState->algorithmParams, "numstreams" );

    XLAL_CHECK_VOID( (fpin = fopen(inputweights, "rb")) != NULL, XLAL_EIO, "Could not open weights file for reading.");

    /* read in the first bit of information, which is the number of datastreams as an UINT4 */
    XLAL_CHECK_VOID( fread((void*)&nstreams, sizeof(UINT4), 1, fpin) == 1, XLAL_EIO, "Could not read in first value from weights file\n");
    XLAL_CHECK_VOID( nstreams == nstreamcheck, XLAL_EFUNC, "Number of data streams inconsistent with ROQ weights file" );

    inputroq = 1;
  }

  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--verbose" );
  if ( ppt ){ verbose = 1; }

  /* check if using a Gaussian likelihoodm and therefore need variance weighted ROQ values */
  if ( LALInferenceGetProcParamVal( runState->commandLine, "--gaussian-like" ) ){ gaussianLike = 1; }

  ifo = runState->threads[0]->model->ifo;

  /* we need to get the frequency factors */
  REAL8Vector *freqFactors = *(REAL8Vector**)LALInferenceGetVariable( ifo->params, "freqfactors" );
  REAL8Vector *freqsCopy = XLALCreateREAL8Vector( freqFactors->length );
  memcpy(freqsCopy->data, freqFactors->data, sizeof(REAL8)*freqFactors->length);
  UINT4 nfreqs = freqsCopy->length;

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

    if ( inputroq ) {
      /* we assume that the data has been split into the same number and length chunks as in that
       * used to calculate the weights, but the second value(s) will be the number of chunks,
       * for the given stream, which can be checked */
      XLAL_CHECK_VOID( fread((void*)&nchunks, sizeof(UINT4), 1, fpin) == 1, XLAL_EIO, "Could not read chunk numbers from weights file\n");
      XLAL_CHECK_VOID( nchunks == chunkLengths->length, XLAL_EFUNC, "Number of chunks is not consistent!");
    }

    if ( outputroq ){
      XLAL_CHECK_VOID( fwrite(&chunkLengths->length, sizeof(UINT4), 1, fpout) == 1, XLAL_EIO, "Could not write number of chunks!" );
    }

    /* as we only use one datastream at a time (e.g. one frequency time series) we need to set the
     * freqFactor to just contain that one frequency */
    REAL8Vector *freqsTemp = XLALCreateREAL8Vector( 1 );
    freqsTemp->data[0] = freqsCopy->data[counter%nfreqs];
    check_and_add_fixed_variable( ifo->params, "freqfactors", &freqsTemp, LALINFERENCE_REAL8Vector_t );

    /* get chunk */
    for ( i=0; i<chunkLengths->length; i++ ){
      UINT4 tlen = chunkLengths->data[i];
      LALInferenceCOMPLEXROQInterpolant *interp = NULL;
      size_t nbases0 = 0;

      if ( !inputroq ){
        /* a temporary run state containing just the required data for a single detector */
        LALInferenceRunState *tmpRS = XLALMalloc(sizeof(LALInferenceRunState));
        LALInferenceIFOModel *ifotmp = XLALMalloc(sizeof(LALInferenceIFOModel));
        ifotmp->next = NULL;
        tmpRS->threads = LALInferenceInitThreads(1);
        tmpRS->threads[0]->model = XLALMalloc(sizeof(LALInferenceModel));
        tmpRS->threads[0]->model->templt = runState->threads[0]->model->templt;
        tmpRS->threads[0]->currentParams = runState->threads[0]->currentParams;
        tmpRS->threads[0]->model->params = XLALCalloc(1, sizeof(LALInferenceVariables));
        tmpRS->threads[0]->model->ifo = ifotmp;
        tmpRS->GSLrandom = runState->GSLrandom;
        tmpRS->priorArgs = runState->priorArgs;

        gsl_matrix_complex *ts = NULL;

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
        COMPLEX16 *RB = NULL;
        RB = LALInferenceGenerateCOMPLEX16OrthonormalBasis(&weights.vector, tolerance, ts, &nbases0);
        XLAL_CHECK_VOID( RB != NULL, XLAL_EFUNC, "Could not produce basis set");
        nbases->data[i] = nbases0;

        if ( verbose ){
          fprintf(stderr, "Number of reduced bases for ROQ generation is %zu\n", nbases0);
        }

        XLAL_CALLGSL( gsl_matrix_complex_free( ts ) );

        if ( nbases0 > (size_t)tlen-1 || nbases0 == ntraining ){
          fprintf( stderr, "Number of bases is longer than data length (or uses all training points), so not using ROQ for this segment\n" );
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
        LALInferenceClearVariables( tmpRS->threads[0]->model->params );
        XLALFree( tmpRS->threads[0]->model );
        XLALFree( tmpRS->threads );
        XLALFree( tmpRS );

        /* generate the interpolants (pass the data noise variances for weighting the interpolants in the case of using a Gaussian likelihood) */
        if ( nbases0 <= (size_t)tlen-1 && nbases0 < ntraining ){
          gsl_matrix_complex_view RBview = gsl_matrix_complex_view_array((double*)RB, nbases0, tlen);
          interp = LALInferenceGenerateCOMPLEXROQInterpolant(&RBview.matrix);
        }
        else{
          /* if the number of bases is greater than the length just use all the data points in a chunk */
          interp =  XLALMalloc(sizeof(LALInferenceCOMPLEXROQInterpolant));
          interp->B = NULL;
          interp->nodes = XLALMalloc(sizeof(UINT4)*tlen);

          for ( j=0; j<tlen; j++ ){ interp->nodes[j] = (UINT4)j; }
        }

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
        gsl_vector_complex *dmw;
        gsl_matrix_complex *mmw;
        if ( nbases0 <= (size_t)tlen-1 && nbases0 < ntraining ){
          dmw = LALInferenceGenerateCOMPLEX16DataModelWeights(interp->B, &datasub.vector, &vari.vector);
          mmw = LALInferenceGenerateCOMPLEXModelModelWeights(interp->B, &vari.vector);
        }
        else{
          /* if the number of bases is longer than the data then fill in the data-model weights with just the data
           * and the model-model weights with a diagonal matrix filled with the variances */
          dmw = gsl_vector_complex_alloc( tlen );
          mmw = gsl_matrix_complex_calloc( tlen, tlen );

          for ( j=0; j<tlen; j++ ){
            gsl_complex compval1;

            if ( gaussianLike ){ GSL_SET_COMPLEX(&compval1, 1./gsl_vector_get(&vari.vector, j), 0.); }
            else{ GSL_SET_COMPLEX(&compval1, 1./gsl_vector_get(&vari.vector, 0), 0.); }
            gsl_matrix_complex_set(mmw, j, j, compval1);
            if ( gaussianLike ){ compval1 = gsl_complex_div_real(gsl_vector_complex_get(&datasub.vector, j), gsl_vector_get(&vari.vector, j)); }
            else{ compval1 = gsl_complex_div_real(gsl_vector_complex_get(&datasub.vector, j), gsl_vector_get(&vari.vector, 0)); }
            gsl_vector_complex_set(dmw, j, compval1);
          }
        }

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

        /* output the node indices */
        if ( outputroq ){
          XLAL_CHECK_VOID( fwrite(&nbases0, sizeof(size_t), 1, fpout) == 1, XLAL_EIO, "Could not output number of nodes to file." );
          XLAL_CHECK_VOID( fwrite(&interp->nodes[0], sizeof(UINT4), nbases0, fpout) == nbases0, XLAL_EIO, "Could not output interpolation node indices." );
        }

        XLAL_CALLGSL( gsl_vector_complex_free( dmw ) );
        XLAL_CALLGSL( gsl_matrix_complex_free( mmw ) );
      }
      else{
        interp = XLALMalloc(sizeof(LALInferenceCOMPLEXROQInterpolant));
        interp->B = NULL;

        /* read in the number of interpolant nodes for the current chunk and then the node indices */
        XLAL_CHECK_VOID( fread((void*)&nbases0, sizeof(size_t), 1, fpin) == 1, XLAL_EIO, "Could not read number of nodes" );

        if ( verbose ){
          fprintf(stderr, "Number of reduced bases for ROQ generation is %zu\n", nbases0);
        }

        /* read in the number of nodes for this chunk */
        interp->nodes = XLALMalloc(nbases0*sizeof(UINT4));
        XLAL_CHECK_VOID( fread((void*)interp->nodes, sizeof(UINT4), nbases0, fpin) == nbases0, XLAL_EIO, "Could not read in interpolation indices" );

        nbases->data[i] = nbases0;
        mmlength += (nbases0*nbases0);
        dmlength += nbases0;
      }

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

    if ( inputroq ){
      /* now read in all the weights */
      dmweights = XLALMalloc( sizeof(COMPLEX16)*dmlength );
      XLAL_CHECK_VOID( fread((void*)dmweights, sizeof(COMPLEX16), dmlength, fpin) == dmlength, XLAL_EIO, "Could not read in data-model product weights" );

      mmweights = XLALMalloc( sizeof(COMPLEX16)*mmlength );
      XLAL_CHECK_VOID( fread((void*)mmweights, sizeof(COMPLEX16), mmlength, fpin) == mmlength, XLAL_EIO, "Could not read in model-model product weights" );
    }

    if ( outputroq ){
      XLAL_CHECK_VOID( fwrite(dmweights, sizeof(COMPLEX16), dmlength, fpout) == dmlength, XLAL_EIO, "Could not output data-model product weights" );
      XLAL_CHECK_VOID( fwrite(mmweights, sizeof(COMPLEX16), mmlength, fpout) == mmlength, XLAL_EIO, "Could not output model-model product weights" );
    }

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

    /* add interpolation weights and nodes to a variable in runState->threads[0]->model->ifo->params */
    LALInferenceAddVariable( ifo->params, "numBases", &nbases, LALINFERENCE_UINT4Vector_t, LALINFERENCE_PARAM_FIXED );

    /* reset freqfactors to the correct value */
    check_and_add_fixed_variable( ifo->params, "freqfactors", &freqsCopy, LALINFERENCE_REAL8Vector_t );

    ifo = ifo->next;
    data = data->next;
    counter++;
  }

  if ( inputroq ){ fclose(fpin); }

  if ( LALInferenceCheckVariable( runState->algorithmParams, "timefile" ) ){
    gettimeofday(&time2, NULL);

    FILE *timefile = *(FILE **)LALInferenceGetVariable( runState->algorithmParams, "timefile" );
    UINT4 timenum = *(UINT4 *)LALInferenceGetVariable( runState->algorithmParams, "timenum" );
    tottime = (REAL8)((time2.tv_sec + time2.tv_usec*1.e-6) - (time1.tv_sec + time1.tv_usec*1.e-6));
    fprintf(timefile, "[%d] %s: %.9le secs\n", timenum, __func__, tottime);
    timenum++;
    check_and_add_fixed_variable( runState->algorithmParams, "timenum", &timenum, LALINFERENCE_UINT4_t );
  }

  if ( outputroq ){
    fclose(*(FILE **)LALInferenceGetVariable( runState->algorithmParams, "timefile" ));
    fclose(fpout);
    if ( verbose ){ fprintf(stderr, "ROQ weights have been written to file.\nExiting program.\n"); }
    exit(0); /* exit the programme succussfully */
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
  gsl_matrix_complex *ts = gsl_matrix_complex_alloc(n, rs->threads[0]->model->ifo->times->length);
  REAL8 *fnodes = NULL;

  for ( j=0; j<n; j++ ){
    /* choose random variables values and fill in runState->threads[0]->model->params */
    LALInferenceVariableItem *item = rs->threads[0]->currentParams->head;

    /* there's no need to re-scale values, so I can drawn from the scaled values */
    for(; item; item = item->next ){
      REAL8 value;

      if( item->vary == LALINFERENCE_PARAM_FIXED || item->vary == LALINFERENCE_PARAM_OUTPUT ){
        LALInferenceAddVariable( rs->threads[0]->model->params, item->name, item->value, item->type, item->vary );
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

        LALInferenceAddVariable( rs->threads[0]->model->params, item->name, &value, item->type, item->vary );
      }
    }

    /* generate model */
    rs->threads[0]->model->templt( rs->threads[0]->model );

    /* place model into an array */
    gsl_vector_complex_view cview;
    cview = gsl_vector_complex_view_array((double*)rs->threads[0]->model->ifo->compTimeSignal->data->data, rs->threads[0]->model->ifo->times->length);
    gsl_matrix_complex_set_row(ts, j, &cview.vector);
  }

  if( fnodes != NULL ){ XLALFree( fnodes ); }

  return ts;
}


/*--------------- END OF ROQ FUNCTIONS ----------------------*/
