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

/******************************************************************************/
/*                       SOFTWARE INJECTION FUNCTIONS                         */
/******************************************************************************/

#include "ppe_inject.h"

/**
 * \brief Inject a simulated signal into the data
 *
 * This function will create an simulated signal (of the required model) to inject into the data from multiple
 * detectors. The parameters of the signal to be injected must be specified in a TEMPO-stype .par file given with the
 * \c inject-file command line argument. The parameters do not have to be the same as those in the .par file controlling
 * the analysis (although should ideally contain a signal within the bandwidth of the data).
 *
 * If a signal of a specific signal-to-noise ratio is required then the \c scale-snr command line argument can be used
 * to give the multi-detector SNR to which the signal needs to be scaled.
 *
 * The injected signal can be output if \c inject-output is set. Two files will be output: one containing the signal
 * only, and one containing the signal plus noise. These will both be in the format of a standard data input file. The
 * files will have names given by the \c inject-output value, with a prefix of the detector name, and a suffix of of \c
 * _signal_only, respectively.
 *
 * \param runState [in] the program information structure
 *
 * \sa calculate_time_domain_snr
 */
void inject_signal( LALInferenceRunState *runState ){
  LALInferenceIFOData *data = runState->data;
  LALInferenceIFOModel *ifo_model = runState->model->ifo;

  CHAR *injectfile = NULL, *snrfile = NULL;
  ProcessParamsTable *ppt;
  ProcessParamsTable *commandLine = runState->commandLine;

  PulsarParameters *injpars = XLALCalloc(sizeof(*injpars),1);;

  FILE *fpsnr = NULL; /* output file for SNRs */
  INT4 ndats = 0, j = 1;

  REAL8Vector *freqFactors = NULL;
  REAL8 snrmulti = 0.;
  REAL8 snrscale = 0;

  ppt = LALInferenceGetProcParamVal( commandLine, "--outfile" );
  if( !ppt ){
    fprintf(stderr, "Error... no output file specified!\n");
    exit(0);
  }

  snrfile = XLALStringDuplicate( ppt->value );
  snrfile = XLALStringAppend( snrfile, "_SNR" );

  if( (fpsnr = fopen(snrfile, "w")) == NULL ){
    fprintf(stderr, "Error... cannot open output SNR file!\n");
    exit(0);
  }

  ppt = LALInferenceGetProcParamVal( commandLine, "--inject-file" );
  if( ppt ){
    injectfile = XLALStringDuplicate( ppt->value );

    /* check that the file exists */
    if ( fopen(injectfile, "r") == NULL ){
      fprintf(stderr, "Error... Injection specified, but the injection parameter file %s is wrong.\n", injectfile);
      exit(3);
    }

    /* read in injection parameter file */
    injpars = XLALReadTEMPOParFileNew( injectfile );

    /* check RA and DEC are set (if only RAJ and DECJ are given in the par file) */
    if ( !PulsarCheckParam( injpars, "RA" ) ){
      if ( PulsarCheckParam( injpars, "RAJ" ) ){
        REAL8 ra = PulsarGetREAL8Param( injpars, "RAJ" );
        PulsarAddParam( injpars, "RA", &ra, PULSARTYPE_REAL8_t );
      }
      else {
        XLALPrintError ("%s: No source right ascension specified!", __func__ );
        XLAL_ERROR_VOID( XLAL_EINVAL );
      }
    }
    if ( !PulsarCheckParam( injpars, "DEC" ) ){
      if ( PulsarCheckParam( injpars, "DECJ" ) ) {
        REAL8 dec = PulsarGetREAL8Param( injpars, "DECJ" );
        PulsarAddParam( injpars, "DEC", &dec, PULSARTYPE_REAL8_t );
      }
      else {
        XLALPrintError ("%s: No source declination specified!", __func__ );
        XLAL_ERROR_VOID( XLAL_EINVAL );
      }
    }

    /* make sure that we have parameters in terms of amplitude and phase parameters */
    invert_source_params( injpars );
  }
  else{
    fclose( fpsnr );
    return;
  }

  freqFactors = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "freqfactors" );

  /* get the SNR scale factor if required */
  ppt = LALInferenceGetProcParamVal( commandLine, "--scale-snr" );
  if ( ppt ) { snrscale = atof( ppt->value ); }

  /* create signal to inject */
  /* for injection always attempt to include the signal phase model even if the search is not going to be over phase */
  INT4 varyphase = 1;

  if ( !LALInferenceCheckVariable( ifo_model->params, "varyphase" ) ){
    LALInferenceAddVariable( ifo_model->params, "varyphase", &varyphase, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );
    varyphase = 0; /* set to zero so varyphase is removed after the model is created */
  }

  /* check whether to inject a non-GR signal (the default will ALWAYS be GR) */
  UINT4 nonGR_search = LALInferenceCheckVariable( ifo_model->params, "nonGR" );

  char* injection_model = NULL;
  ProcessParamsTable *nonGR_injection;
  nonGR_injection = LALInferenceGetProcParamVal(commandLine, "--inject-nonGR");

  if ( nonGR_injection ){
    /* set injection model */
    injection_model = XLALStringDuplicate( nonGR_injection->value );
    if ( !nonGR_search ){
      /* add "nonGR" flag so that pulsar_model() runs in nonGR mode */
      UINT4 nonGRval = 1;
      while ( ifo_model ){
        LALInferenceAddVariable( ifo_model->params, "nonGR", &nonGRval, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );
        ifo_model = ifo_model->next;
      }
      ifo_model = runState->model->ifo;
      /* setup nonGR lookup tables */
      LALSource psr;
      psr.equatorialCoords.longitude = PulsarGetREAL8ParamOrZero( injpars, "RA" );
      psr.equatorialCoords.latitude = PulsarGetREAL8ParamOrZero( injpars, "DEC" );
      psr.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
      setup_lookup_tables( runState, &psr );
    }
    if ( *injection_model!='\0' ){
      /* set parameters corresponding to specific model */
      set_nonGR_model_parameters( injpars, injection_model);
    }
  }
  else {
    /* remove "nonGR" flag so that pulsar_model() runs in GR mode */
      while ( ifo_model ){
        LALInferenceRemoveVariable(ifo_model->params, "nonGR");
        ifo_model = ifo_model->next;
      }
      ifo_model = runState->model->ifo;
  }

  pulsar_model( injpars, ifo_model );

  /* get summed data for use in SNR calculation */
  sum_data( runState );

  fprintf(fpsnr, "# Injected SNR\n");

  /* reset model to head */
  ifo_model = runState->model->ifo;

  /* calculate SNRs */
  while ( data ){
    REAL8 snrval = 0.;

    if ( !LALInferenceCheckVariable( ifo_model->params, "varyphase" ) ){
      LALInferenceAddVariable( ifo_model->params, "varyphase", &varyphase, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );
    }

    snrval = calculate_time_domain_snr( data, ifo_model );

    /* reset varyphase to its original value */
    if ( !varyphase ){ LALInferenceRemoveVariable( ifo_model->params, "varyphase" ); }

    snrmulti += SQUARE(snrval);

    /* if not scaling print out individual detector/datastream SNRs */
    if( snrscale == 0. ){
      fprintf(fpsnr, "%s\t%.3lf\t%le\n", data->name, freqFactors->data[ndats%(INT4)freqFactors->length], snrval);
    }

    data = data->next;
    ifo_model = ifo_model->next;

    ndats++;
  }

  /* get overall multi-detector SNR */
  snrmulti = sqrt(snrmulti);

  /* only need to print out multi-detector snr if the were multiple detectors or data streams */
  if ( snrscale == 0. ){
    if ( ndats > 1 ) { fprintf(fpsnr, "Coherent\t%le\n", snrmulti); }
  }
  else{
    /* rescale the signal and calculate the SNRs */
    snrscale /= snrmulti;

    REAL8 C22 = PulsarGetREAL8ParamOrZero( injpars, "C22" )*snrscale;
    REAL8 C21 = PulsarGetREAL8ParamOrZero( injpars, "C21" )*snrscale;
    PulsarAddParam( injpars, "C22", &C22, PULSARTYPE_REAL8_t );
    PulsarAddParam( injpars, "C21", &C21, PULSARTYPE_REAL8_t );

    /* recreate the signal with scaled amplitude */
    varyphase = 1;

    /* reset to head */
    ifo_model = runState->model->ifo;
    data = runState->data;

    if ( !LALInferenceCheckVariable( ifo_model->params, "varyphase" ) ){
      LALInferenceAddVariable( ifo_model->params, "varyphase", &varyphase, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );
      varyphase = 0; /* set to zero so varyphase is removed after the model is created */
    }

    pulsar_model( injpars, ifo_model );

    /* get new summed data for use in SNR calculation */
    sum_data( runState );

    /* get new snrs */
    snrmulti = 0;
    ndats = 0;

    ifo_model = runState->model->ifo;

    while( data ){
      REAL8 snrval = 0.;

      if ( !LALInferenceCheckVariable( ifo_model->params, "varyphase" ) ){
        LALInferenceAddVariable( ifo_model->params, "varyphase", &varyphase, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );
      }

      /* recalculate the SNR */
      snrval = calculate_time_domain_snr( data, ifo_model );

      /* reset varyphase to its original value */
      if ( !varyphase ){ LALInferenceRemoveVariable( ifo_model->params, "varyphase" ); }

      snrmulti += SQUARE(snrval);

      fprintf(fpsnr, "%s\t%.3lf\t%le\n", data->name, freqFactors->data[ndats%(INT4)freqFactors->length], snrval);

      data = data->next;
      ifo_model = ifo_model->next;

      ndats++;
    }

    snrmulti = sqrt( snrmulti );
    //fprintf(stderr, "scaled multi data snr: %le\n", snrmulti);

    if( ndats > 1 ){
      fprintf(fpsnr, "Coherent\t%le\n", snrmulti);
    }
  }

  fclose( fpsnr );

  data = runState->data;
  ifo_model = runState->model->ifo;

  /* add signal to data */
  while( data ){
    FILE *fp = NULL, *fpso = NULL;
    ProcessParamsTable *ppt2 = LALInferenceGetProcParamVal( commandLine, "--inject-output" );
    INT4 i = 0, length = ifo_model->times->length;

    /* check whether to output the data */
    if ( ppt2 ){
      /* add the site prefix to the start of the output name */
      CHAR *outfile = NULL;
      CHAR *signalonly = NULL; /* file containing only signal and no noise */
      CHAR suffix[5];
      INT4 sf;

      outfile = XLALStringDuplicate( ppt2->value );

      /* append detector name to file */
      outfile = XLALStringAppend( outfile, "_" );
      outfile = XLALStringAppend( outfile, data->detector->frDetector.prefix );

      /* append the harmonic frequency of the signal in the file */
      sf = sprintf(suffix, "_%.1lf", freqFactors->data[(INT4)fmod(j,freqFactors->length)]);
      outfile = XLALStringAppend( outfile, suffix );

      if ( (fp = fopen(outfile, "w")) == NULL || !sf ){
        fprintf(stderr, "Non-fatal error... unable to open file %s to output injection\n", outfile);
      }

      signalonly = XLALStringDuplicate( outfile );
      signalonly = XLALStringAppend( signalonly, "_signal_only" );

      if ( (fpso = fopen(signalonly, "w")) == NULL ){
        fprintf(stderr, "Non-fatal error... unable to open file %s to output injection\n", signalonly);
      }
    }

    /* add the signal to the data */
    for ( i = 0; i < length; i++ ){
      data->compTimeData->data->data[i] += ifo_model->compTimeSignal->data->data[i];

      /* write out injection to file */
      if( fp != NULL && fpso != NULL ){
        /* print out data - time stamp, real and imaginary parts of data (injected signal + noise) */
        fprintf(fp, "%.5lf\t%le\t%le\n", XLALGPSGetREAL8( &ifo_model->times->data[i] ),
                creal(data->compTimeData->data->data[i]), cimag(data->compTimeData->data->data[i]) );

        /* print signal only data - time stamp, real and imaginary parts of signal */
        fprintf(fpso, "%.5lf\t%le\t%le\n", XLALGPSGetREAL8( &ifo_model->times->data[i] ),
                creal(ifo_model->compTimeSignal->data->data[i]), cimag(ifo_model->compTimeSignal->data->data[i]) );
      }
    }

    if ( fp != NULL ) { fclose( fp ); }
    if ( fpso != NULL ) { fclose( fpso ); }

    /* recalculate the data noise variance with the signal added (this mimics the variance that
     * would be calculated for a real signal.
     */
    compute_variance( data, ifo_model );

    data = data->next;
    ifo_model = ifo_model->next;
    j++;
  }

  /* reset nonGR status */
  ifo_model = runState->model->ifo;
  if ( nonGR_search ){
    if ( !nonGR_injection ){
      UINT4 nonGRval = 1;
      while ( ifo_model ){
        LALInferenceAddVariable( ifo_model->params, "nonGR", &nonGRval, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );
        ifo_model = ifo_model->next;
      }
      ifo_model = runState->model->ifo;
    }
  }
  else {
    while ( ifo_model ){
      LALInferenceRemoveVariable(ifo_model->params, "nonGR");
      ifo_model = ifo_model->next;
    }
    ifo_model = runState->model->ifo;
  }

  PulsarFreeParams( injpars ); /* free memory */
}

/*-------------------- END OF SOFTWARE INJECTION FUNCTIONS -------------------*/



/**
 * \brief Calculates the optimal matched filter signal-to-noise ratio for a given signal
 *
 * This function calculates the optimal matched filter signal-to-noise ratio (SNR) of a given signal model for a set of
 * detector data via:
 * \f[
 * \rho = \sqrt{\sum_{i=1}^N \frac{d_i^2}{\sigma^2}},
 * \f]
 * where \f$\{d\}\f$ is a time series of data, and \f$\sigma^2\f$ is its variance. As the data and model used here are
 * complex the real and imaginary SNRs are added in quadrature to give the total SNR.
 *
 * The data variance \f$\sigma^2\f$ is calculated on data that has had the running median subtracted in order to remove
 * any underlying trends (e.g. caused by a string signal). The variance is assumed constant over segments given in the
 * \c chunkLength vector and the SNR from each segment is added in quadrature.
 *
 * \param data [in] A data pointer containing detector data and the signal model
 * \param ifo_model [in] A model structure containing detector parameters and buffers
 *
 * \return The optimal matched filter signal-to-noise ratio
 */
REAL8 calculate_time_domain_snr( LALInferenceIFOData *data, LALInferenceIFOModel *ifo_model ){
  REAL8 snrval = 0., chunkLength = 0;

  INT4 i = 0, j = 0, length = 0, cl = 0;

  INT4 chunkMin = 0, count = 0, varyphase = 0, nonGR = 0, roq = 0;

  UINT4Vector *chunkLengths = NULL;
  chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( ifo_model->params, "chunkLength" );
  chunkMin = *(INT4*)LALInferenceGetVariable( ifo_model->params, "chunkMin" );

  REAL8Vector *sumP = NULL, *sumC = NULL, *sumX = NULL, *sumY = NULL, *sumB = NULL, *sumL = NULL;
  REAL8Vector *sumPC = NULL, *sumPX = NULL, *sumPY = NULL, *sumPB = NULL, *sumPL = NULL;
  REAL8Vector *sumCX = NULL, *sumCY = NULL, *sumCB = NULL, *sumCL = NULL;
  REAL8Vector *sumXY = NULL, *sumXB = NULL, *sumXL = NULL;
  REAL8Vector *sumYB = NULL, *sumYL = NULL;
  REAL8Vector *sumBL = NULL;

  if ( LALInferenceCheckVariable( ifo_model->params, "varyphase" ) ){ varyphase = 1; }
  if ( LALInferenceCheckVariable( ifo_model->params, "nonGR" ) ){ nonGR = 1; }
  if ( LALInferenceCheckVariable( ifo_model->params, "roq" ) ){ roq = 1; }

  if ( !varyphase && !roq ){
    /* get explicitly whitened versions of these sums */
    sumP = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumPWhite" );
    sumC = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumCWhite" );
    sumPC = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumPCWhite" );

    if ( nonGR ){
      sumX = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumXWhite" );
      sumY = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumYWhite" );
      sumB = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumBWhite" );
      sumL = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumLWhite" );

      sumPX = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumPXWhite" );
      sumPY = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumPYWhite" );
      sumPB = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumPBWhite" );
      sumPL = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumPLWhite" );
      sumCX = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumCXWhite" );
      sumCY = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumCYWhite" );
      sumCB = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumCBWhite" );
      sumCL = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumCLWhite" );
      sumXY = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumXYWhite" );
      sumXB = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumXBWhite" );
      sumXL = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumXLWhite" );
      sumYB = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumYBWhite" );
      sumYL = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumYLWhite" );
      sumBL = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "sumBLWhite" );
    }
  }

  length = data->compTimeData->data->length;

  for ( i = 0; i < length; i+=(INT4)chunkLength ){
    REAL8 snrcRe = 0., snrcIm = 0.;

    chunkLength = (REAL8)chunkLengths->data[count];

    /* skip section of data if its length is less than the minimum allowed chunk length */
    if( chunkLength < chunkMin ){
      count++;
      continue;
    }

    cl = i + (INT4)chunkLength;

    /* check whether the full time domain model exists, or the pre-summed model */
    if ( varyphase || roq ){
      for( j = i ; j < cl ; j++ ){
        /* calculate optimal signal power */
        snrcRe += SQUARE( creal(ifo_model->compTimeSignal->data->data[j]) ) / data->varTimeData->data->data[j];
        snrcIm += SQUARE( cimag(ifo_model->compTimeSignal->data->data[j]) ) / data->varTimeData->data->data[j];
      }
    }
    else{
      COMPLEX16 Mp, Mc;
      Mp = ifo_model->compTimeSignal->data->data[0];
      Mc = ifo_model->compTimeSignal->data->data[1];

      snrcRe += sumP->data[count]*(creal(Mp)*creal(Mp) + cimag(Mp)*cimag(Mp)) +
                sumC->data[count]*(creal(Mc)*creal(Mc) + cimag(Mc)*cimag(Mc)) +
                2.*sumPC->data[count]*(creal(Mp)*creal(Mc) + cimag(Mp)*cimag(Mc));
      snrcIm = 0.;

      if ( nonGR ){
        COMPLEX16 Mx, My, Mb, Ml;

        Mx = ifo_model->compTimeSignal->data->data[2];
        My = ifo_model->compTimeSignal->data->data[3];
        Mb = ifo_model->compTimeSignal->data->data[4];
        Ml = ifo_model->compTimeSignal->data->data[5];

        snrcRe += sumX->data[count]*(creal(Mx)*creal(Mx) + cimag(Mx)*cimag(Mx)) +
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
      }
    }

    /* add SNRs for each chunk in quadrature */
    snrval += ( snrcRe + snrcIm );

    count++;
  }

  snrval = sqrt( snrval );

  return snrval;
}


/**
 * \brief Get the signal-to-noise ratio of the maximum likelihood signal
 *
 * The function uses the signal with the highest likelihood (which will be the final point in the live points array) and
 * calculates the optimal signal-to-noise ratio (SNR) for it. This is output to a file based on the \c outfile value,
 * but with \c _SNR appended to it. For multiple detector, and/or models with multiple data sets, the individual
 * detector/data set SNR values will be output, with the final value being the multi-detector SNR. If a fake signal has
 * been injected into the data this file will already contain the optimal SNR of the true signal.
 *
 * \param runState [in] The analysis information structure
 *
 * \sa calculate_time_domain_snr
 */
void get_loudest_snr( LALInferenceRunState *runState ){
  INT4 ndats = 0;
  UINT4 roq = 0;
  INT4 Nlive = *(INT4 *)LALInferenceGetVariable( runState->algorithmParams, "Nlive" );
  REAL8 snrmulti = 0.;
  REAL8Vector *freqFactors = NULL;

  CHAR *snrfile = NULL;
  FILE *fpsnr = NULL;

  ProcessParamsTable *ppt;
  ProcessParamsTable *commandLine = runState->commandLine;

  LALInferenceVariables *loudestParams = NULL;
  LALInferenceIFOData *data = runState->data;
  LALInferenceIFOModel *ifo_model = runState->model->ifo;

  loudestParams = XLALCalloc( 1, sizeof(LALInferenceVariables) );

  /* max likelihood point should have been sorted to be the final value */
  LALInferenceCopyVariables( runState->livePoints[Nlive-1], loudestParams );

  /* if using ROQ we need to reinstate the full time stamp vector to compute the model for the SNR calculation */
  if ( LALInferenceGetProcParamVal( runState->commandLine, "--roq" ) ){
    roq = 1;

    while ( ifo_model ){
      UINT4 varyphase = 1; /* must set varyphase */

      REAL8Vector *sidtime = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "siderealDayFull" );
      LALInferenceRemoveVariable( ifo_model->params, "siderealDay" );
      LALInferenceAddVariable( ifo_model->params, "siderealDay", &sidtime, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );

      LIGOTimeGPSVector *timestamps = *(LIGOTimeGPSVector **)LALInferenceGetVariable( ifo_model->params, "timeStampVectorFull" );
      XLALDestroyTimestampVector( ifo_model->times );
      ifo_model->times = timestamps;

      REAL8TimeSeries *timedata = *(REAL8TimeSeries **)LALInferenceGetVariable( ifo_model->params, "timeDataFull" );
      XLALDestroyREAL8TimeSeries( ifo_model->timeData );
      ifo_model->timeData = timedata;

      ifo_model->compTimeSignal = XLALResizeCOMPLEX16TimeSeries( ifo_model->compTimeSignal, 0, ifo_model->times->length );

      /* remove the ROQ variable for calculating the likelihood */
      LALInferenceRemoveVariable( ifo_model->params, "roq" );

      if ( !LALInferenceCheckVariable( ifo_model->params, "varyphase" ) ){
        LALInferenceAddVariable( ifo_model->params, "varyphase", &varyphase, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED );
      }

      if ( LALInferenceCheckVariable( ifo_model->params, "ssb_delays" ) ){
        REAL8Vector *ssbdelays = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "ssb_delays_full" );
        LALInferenceRemoveVariable( ifo_model->params, "ssb_delays" );
        LALInferenceAddVariable( ifo_model->params, "ssb_delays", &ssbdelays, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
      }

      if ( LALInferenceCheckVariable( ifo_model->params, "bsb_delays" ) ){
        REAL8Vector *bsbdelays = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "bsb_delays_full" );
        LALInferenceRemoveVariable( ifo_model->params, "bsb_delays" );
        LALInferenceAddVariable( ifo_model->params, "bsb_delays", &bsbdelays, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
      }

      ifo_model = ifo_model->next;
    }
  }

  /* make sure that the signal model in runState->data is that of the loudest signal */
  REAL8 logLnew = runState->likelihood( loudestParams, runState->data, runState->model );

  if ( !roq ) {
    /* we don't expect exactly identical likelihoods if using ROQ */
    if ( logLnew != *(REAL8 *)LALInferenceGetVariable(
      runState->livePoints[Nlive-1], "logL" ) ){
      fprintf(stderr, "Error... maximum log likelihood problem!\n");
      exit(0);
    }
  }
  else{
    /* check likelihoods agree to within 0.1% */
    fprintf(stderr, "Max. ROQ likeihood = %.12le, Max. full likelihood = %.12le\n", *(REAL8 *)LALInferenceGetVariable(runState->livePoints[Nlive-1], "logL"), logLnew);
    if ( 100.*(1.-fabs(*(REAL8 *)LALInferenceGetVariable(runState->livePoints[Nlive-1], "logL")/logLnew)) > 0.1 ){
      fprintf(stderr, "WARNING... maximum likelihood using ROQ is greater than 0.1%% different that the full likelihood calculation.\n");
    }
  }

  LALInferenceClearVariables( loudestParams );

  /* setup output file */
  ppt = LALInferenceGetProcParamVal( commandLine, "--outfile" );
  if( !ppt ){
    fprintf(stderr, "Error... no output file specified!\n");
    exit(0);
  }

  snrfile = XLALStringDuplicate( ppt->value );
  snrfile = XLALStringAppend( snrfile, "_SNR" );

  /* append to previous injection SNR file if it exists */
  if( (fpsnr = fopen(snrfile, "a")) == NULL ){
    fprintf(stderr, "Error... cannot open output SNR file!\n");
    exit(0);
  }

  /* get SNR of loudest point and print out to file */
  data = runState->data;
  ifo_model = runState->model->ifo;

  //FILE *fp = NULL;
  //fp = fopen("max_like_signal.txt", "w");

  fprintf(fpsnr, "# Recovered SNR\n");

  freqFactors = *(REAL8Vector **)LALInferenceGetVariable( ifo_model->params, "freqfactors" );

  while( data ){
    if ( roq ) { LALInferenceAddVariable( ifo_model->params, "roq", &roq, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED ); }

    REAL8 snrval = calculate_time_domain_snr( data, ifo_model );

    //UINT4 length = ifo_model->compTimeData->data->length;

    /* print out maxlikelihood template */
    //for ( UINT4 j=0; j < length; j++ ){
    //  fprintf(fp, "%lf\t%le\t%le\n",
    //          XLALGPSGetREAL8( &ifo_model->times->data[j] ),
    //          ifo_model->compTimeSignal->data->data[j].re,
    //          ifo_model->compTimeSignal->data->data[j].im );
    //}

    snrmulti += SQUARE( snrval );

    /* print out SNR value */
    fprintf(fpsnr, "%s\t%.3lf\t%le\n", data->name, freqFactors->data[ndats%(INT4)freqFactors->length], snrval);

    ndats++;

    data = data->next;
    ifo_model = ifo_model->next;
  }

  //fclose(fp);

  /* print out multi-detector/multi-datastream SNR value */
  if ( ndats > 1 ) { fprintf(fpsnr, "Coherent\t%le\n", sqrt(snrmulti)); }

  fclose( fpsnr );
}

