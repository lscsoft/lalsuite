/*
*  Copyright (C) 2015, 2016 Matthew Pitkin
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
 * \brief Reduced order quadrature generation functions for use in parameter estimation
 * codes for targeted pulsar searches.
 */

#include "config.h"
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
  UINT4 ntraining = 0, nenrichmax = 100, verbose = 0, outputroq = 0, inputroq = 0, nconsec = 3;
  UINT4 counter = 0, i = 0, j = 0, conseccount = 0;
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
  UINT4 nquadtot = 0, nlineartot = 0, ndatatot = 0;

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

    /* check if a maximum number of enrichment steps is set */
    ppt = LALInferenceGetProcParamVal( runState->commandLine, "--enrich-max" );
    if ( ppt ){ nenrichmax = atoi(ppt->value); }
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
  UINT4 nfreqs = freqFactors->length;
  REAL8Vector *freqsCopy = XLALCreateREAL8Vector( nfreqs );
  memcpy(freqsCopy->data, freqFactors->data, sizeof(REAL8)*nfreqs);

  while ( ifo ){
    UINT4Vector *chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( ifo->params, "chunkLength" );
    REAL8 dt = *(REAL8*)LALInferenceGetVariable( ifo->params, "dt" );
    INT4 startidx = 0;

    UINT4 dlen = 0;
    UINT4Vector *nbases = XLALCreateUINT4Vector( chunkLengths->length );
    UINT4Vector *nbasesquad = XLALCreateUINT4Vector( chunkLengths->length );

    /* vectors to hold the vector and matrices with the interpolation weights */
    COMPLEX16 *dmweights = NULL;
    REAL8 *mmweights = NULL;
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
      UINT4 tlen = chunkLengths->data[i], enrichcounts = 0, nbases0 = 0, nbases0quad = 0;
      LALInferenceCOMPLEXROQInterpolant *interp = NULL;
      LALInferenceREALROQInterpolant *interpquad = NULL;

      if ( !inputroq ){
        COMPLEX16Array *ts = NULL; /* training set for linear part of interpolant */
        REAL8Array *tsquad = NULL; /* training set for quadratic part of interpolant */
        REAL8 maxprojerr = 0.;

        /* a temporary run state containing just the required data for a single detector */
        LALInferenceRunState *tmpRS = XLALMalloc(sizeof(LALInferenceRunState));
        LALInferenceIFOModel *ifotmp = XLALMalloc(sizeof(LALInferenceIFOModel));
        ifotmp->next = NULL;
        tmpRS->threads = LALInferenceInitThreads(1);
        tmpRS->threads[0]->model = XLALMalloc(sizeof(LALInferenceModel));
        tmpRS->threads[0]->model->templt = runState->threads[0]->model->templt;
        LALInferenceClearVariables(tmpRS->threads[0]->currentParams); /* free memory as LALInferenceInitThreads already has allocated memory */
        tmpRS->threads[0]->currentParams = runState->threads[0]->currentParams;
        tmpRS->threads[0]->model->ifo = ifotmp;
        tmpRS->GSLrandom = runState->GSLrandom;
        tmpRS->priorArgs = runState->priorArgs;
        tmpRS->commandLine = runState->commandLine;

        ifotmp->times = XLALCreateTimestampVector( tlen );
        /* create a new model vector */
        ifotmp->compTimeSignal = XLALCreateCOMPLEX16TimeSeries( "", &gpstime, 0., 1., &lalSecondUnit, tlen );
        ifotmp->timeData = XLALCreateREAL8TimeSeries( "", &gpstime, 0., 1., &lalSecondUnit, tlen );

        ifotmp->params = XLALCalloc(1, sizeof(LALInferenceVariables));
        LALInferenceCopyVariables(ifo->params, ifotmp->params); /* copy parameters */
        ifotmp->ephem = ifo->ephem;
        ifotmp->detector = ifo->detector;
        ifotmp->tdat = ifo->tdat;
        ifotmp->ttype = ifo->ttype;

        REAL8Vector *thissiddayfrac = NULL;
        thissiddayfrac = XLALCreateREAL8Vector( tlen );
        LALInferenceRemoveVariable( ifotmp->params, "siderealDay" );

        /* get chunk times */
        for ( j=0; j<tlen; j++ ){
          ifotmp->times->data[j] = ifo->times->data[startidx+j];
          ifotmp->timeData->data->data[j] = ifo->timeData->data->data[startidx+j];
          thissiddayfrac->data[j] = sidDayFrac->data[startidx+j];
        }

        LALInferenceAddVariable( ifotmp->params, "siderealDay", &thissiddayfrac, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
        if ( !LALInferenceCheckVariable( ifotmp->params, "varyskypos" ) ){
          UINT4 vsp = 1;
          LALInferenceAddVariable( ifotmp->params, "varyskypos", &vsp, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED );
        }

        REAL8Vector *deltas = XLALCreateREAL8Vector( 1 );
        deltas->data[0] = dt;

        if ( verbose ){ fprintf(stderr, "Data chunk no. %d of %d (length: %d)\n", i+1, chunkLengths->length, tlen); }

        /* GENERATE TRAINING SET FOR PART OF ROQ THAT IS LINEAR WRT THE MODEL <d|h> */
        ts = generate_training_set( tmpRS, ntraining ); /* generate the training set */

        COMPLEX16Array *RB = NULL; /* the reduced basis */
        maxprojerr = LALInferenceGenerateCOMPLEX16OrthonormalBasis(&RB, deltas, tolerance, ts); /* generate reduced basis */
        XLAL_CHECK_VOID( RB != NULL, XLAL_EFUNC, "Could not produce linear basis set" );
        nbases->data[i] = RB->dimLength->data[0];

        /* iterate over enrichment steps */
        if ( nenrichmax > 0 && nbases->data[i] < tlen ) {
          enrichcounts = 0;
          conseccount = 0; /* counter for number of consecutive enrichments that add no new bases */
          if( verbose ){ fprintf(stderr, "Enriching linear basis (no. bases): %d", nbases->data[i]); }

          do {
            /* create a new training set to try and "enrich" the original one */
            COMPLEX16Array *tsenrich = NULL;
            tsenrich = generate_training_set( tmpRS, ntraining );

            maxprojerr = LALInferenceEnrichCOMPLEX16Basis(deltas, tolerance, &RB, &ts, tsenrich); /* regenerate reduced basis */

            /* check if no new bases have been added */
            if ( nbases->data[i] == RB->dimLength->data[0] ){
              /* check if LALInferenceGenerateCOMPLEX16OrthonormalBasis exited before reaching the required tolerance */
              if ( maxprojerr < tolerance ){ conseccount++; }
            }
            else{ conseccount = 0; }

            nbases->data[i] = RB->dimLength->data[0];

            /* check if there are more bases than actual data points */
            if ( nbases->data[i] >= tlen ){
              XLALDestroyCOMPLEX16Array( tsenrich );
              nbases->data[i] = tlen;
              break;
            }

            if ( verbose ){ fprintf(stderr, "...%d", nbases->data[i]); }
            XLALDestroyCOMPLEX16Array( tsenrich );
            enrichcounts++;
          } while( enrichcounts < nenrichmax && conseccount < nconsec );
          if ( verbose ){ fprintf(stderr, "\n"); }
        }

        XLALDestroyCOMPLEX16Array( ts );

        if ( verbose ){
          fprintf(stderr, "Number of linear reduced bases for ROQ generation is %d, with a maximum projection error of %le\n", nbases->data[i], maxprojerr);
        }

        nbases0 = nbases->data[i];
        nlineartot += nbases0;

        /* generate the interpolants (pass the data noise variances for weighting the interpolants in the case of using a Gaussian likelihood) */
        if ( nbases0 <= tlen-1 ){
          interp = LALInferenceGenerateCOMPLEXROQInterpolant(RB);
        }
        else{
          /* if the number of bases is greater than the length just use all the data points in a chunk */
          interp = XLALMalloc(sizeof(LALInferenceCOMPLEXROQInterpolant));
          interp->B = NULL;
          interp->nodes = XLALMalloc(sizeof(UINT4)*tlen);

          for ( j=0; j<tlen; j++ ){ interp->nodes[j] = (UINT4)j; }
        }

        /* free the reduced basis set and training set */
        XLALDestroyCOMPLEX16Array( RB );

        /* GENERATE TRAINING SET FOR PART OF ROQ THAT IS QUADRATIC WRT THE MODEL <h|h> */
        tsquad = generate_training_set_quad( tmpRS, ntraining ); /* generate the training set */

        REAL8Array *RBquad = NULL; /* the reduced basis */
        maxprojerr = LALInferenceGenerateREAL8OrthonormalBasis(&RBquad, deltas, tolerance, tsquad); /* generate reduced basis */
        XLAL_CHECK_VOID( RBquad != NULL, XLAL_EFUNC, "Could not produce quadratic basis set" );
        nbasesquad->data[i] = RBquad->dimLength->data[0];

        /* iterate over enrichment steps */
        if ( nenrichmax > 0 && nbasesquad->data[i] < tlen ) {
          enrichcounts = 0;
          conseccount = 0; /* counter for number of consecutive enrichments that add no new bases */

          if( verbose ){ fprintf(stderr, "Enriching quadratic basis (no. bases): %d", nbasesquad->data[i]); }
          do {
            /* create a new training set to try and "enrich" the original one */
            REAL8Array *tsenrich = NULL;
            tsenrich = generate_training_set_quad( tmpRS, ntraining );

            maxprojerr = LALInferenceEnrichREAL8Basis(deltas, tolerance, &RBquad, &tsquad, tsenrich); /* regenerate reduced basis */

            /* check if no new bases have been added */
            if ( nbasesquad->data[i] == RBquad->dimLength->data[0] ){
              /* check if LALInferenceGenerateREAL8OrthonormalBasis exited before reaching the required tolerance */
              if ( maxprojerr < tolerance ){ conseccount++; }
            }
            else{ conseccount = 0; }

            nbasesquad->data[i] = RBquad->dimLength->data[0];

            /* check if there are more bases than actual data points */
            if ( nbasesquad->data[i] >= tlen ){
              XLALDestroyREAL8Array( tsenrich );
              nbasesquad->data[i] = tlen;
              break;
            }

            if ( verbose ){ fprintf(stderr, "...%d", nbasesquad->data[i]); }
            XLALDestroyREAL8Array( tsenrich );
            enrichcounts++;
          } while( enrichcounts < nenrichmax && conseccount < nconsec );
          if ( verbose ){ fprintf(stderr, "\n"); }
        }

        XLALDestroyREAL8Array(tsquad);

        if ( verbose ){ fprintf(stderr, "Number of quadratic reduced bases for ROQ generation is %d, with a maximum projection error of %le\n", nbasesquad->data[i], maxprojerr);}

        nbases0quad = nbasesquad->data[i];
        nquadtot += nbases0quad;

        /* generate the interpolants (pass the data noise variances for weighting the interpolants in the case of using a Gaussian likelihood) */
        if ( nbases0quad <= tlen-1 ){
          interpquad = LALInferenceGenerateREALROQInterpolant(RBquad);
        }
        else{
          /* if the number of bases is greater than the length just use all the data points in a chunk */
          interpquad = XLALMalloc(sizeof(LALInferenceREALROQInterpolant));
          interpquad->B = NULL;
          interpquad->nodes = XLALMalloc(sizeof(UINT4)*tlen);

          for ( j=0; j<tlen; j++ ){ interpquad->nodes[j] = (UINT4)j; }
        }

        /* free the reduced basis set and training set */
        XLALDestroyREAL8Array( RBquad );

        XLALDestroyREAL8Vector( deltas );
        XLALDestroyTimestampVector( ifotmp->times );
        XLALDestroyCOMPLEX16TimeSeries( ifotmp->compTimeSignal );
        LALInferenceClearVariables( ifotmp->params );
        XLALDestroyREAL8TimeSeries( ifotmp->timeData );
        XLALFree( tmpRS->threads[0]->model );
        XLALFree( tmpRS->threads );
        XLALFree( tmpRS );

        /* get view of data and variance */
        REAL8Vector vari, *varist = NULL;
        COMPLEX16Vector datasub;

        /* point to section of data */
        datasub.data = &data->compTimeData->data->data[startidx];
        datasub.length = tlen;

        if ( gaussianLike ){
          /* point to section of variance required */
          vari.data = &data->varTimeData->data->data[startidx];
          vari.length = tlen;
        }
        else{
          varist = XLALCreateREAL8Vector( 1 );
          varist->data[0] = 1.; /* using Students-t likelihood, so just set to 1 */
        }

        /* create the data/model and model/model inner product weights */
        COMPLEX16Vector *dmw = NULL;
        REAL8Vector *mmw = NULL;
        if ( nbases0 <= tlen-1 ){
          if ( gaussianLike ) {
            dmw = LALInferenceGenerateCOMPLEX16LinearWeights(interp->B, &datasub, &vari);
          }
          else{
            dmw = LALInferenceGenerateCOMPLEX16LinearWeights(interp->B, &datasub, varist);
          }
        }
        else{
          /* if the number of bases is longer than the data then fill in the data-model weights with just the data */
          dmw = XLALCreateCOMPLEX16Vector( tlen );

          for ( j=0; j<tlen; j++ ){
            if ( gaussianLike ){ dmw->data[j] = datasub.data[j]/vari.data[j]; }
            else{ dmw->data[j] = datasub.data[j]/varist->data[0]; }
          }
        }

        if ( nbases0quad <= tlen-1 ){
          if ( gaussianLike ) {
            mmw = LALInferenceGenerateQuadraticWeights(interpquad->B, &vari);
          }
          else{
            mmw = LALInferenceGenerateQuadraticWeights(interpquad->B, varist);
          }
        }
        else{
          /* if the number of bases is longer than the data then fill the model-model weights filled with one over the variances */
          mmw = XLALCreateREAL8Vector( tlen );

          for ( j=0; j<tlen; j++ ){
            if ( gaussianLike ){ mmw->data[j] = 1./vari.data[j]; }
            else{ mmw->data[j] = 1./varist->data[0]; }
          }
        }

        if ( !gaussianLike ) { XLALDestroyREAL8Vector( varist ); }

        /* put weights into a vector */
        mmlength += nbases0quad;
        dmlength += nbases0;

        dmweights = XLALRealloc( dmweights, sizeof(COMPLEX16)*dmlength );
        mmweights = XLALRealloc( mmweights, sizeof(REAL8)*mmlength );

        memcpy(&dmweights[dmlength-nbases0], dmw->data, sizeof(COMPLEX16)*nbases0);
        memcpy(&mmweights[mmlength-nbases0quad], mmw->data, sizeof(REAL8)*nbases0quad);

        /* output the node indices */
        if ( outputroq ){
          XLAL_CHECK_VOID( fwrite(&nbases0, sizeof(UINT4), 1, fpout) == 1, XLAL_EIO, "Could not output number of linear nodes to file." );
          XLAL_CHECK_VOID( fwrite(&interp->nodes[0], sizeof(UINT4), nbases0, fpout) == nbases0, XLAL_EIO, "Could not output linear interpolation node indices." );
          XLAL_CHECK_VOID( fwrite(&nbases0quad, sizeof(UINT4), 1, fpout) == 1, XLAL_EIO, "Could not output number of quadratic nodes to file." );
          XLAL_CHECK_VOID( fwrite(&interpquad->nodes[0], sizeof(UINT4), nbases0quad, fpout) == nbases0quad, XLAL_EIO, "Could not output quadratic interpolation node indices." );
        }

        XLALDestroyCOMPLEX16Vector( dmw );
        XLALDestroyREAL8Vector( mmw );
      }
      else{ /* reading in previously computed interpolant nodes */
        interp = XLALMalloc(sizeof(LALInferenceCOMPLEXROQInterpolant));
        interp->B = NULL;

        /* read in the number of linear interpolant nodes for the current chunk and then the node indices */
        XLAL_CHECK_VOID( fread((void*)&nbases0, sizeof(UINT4), 1, fpin) == 1, XLAL_EIO, "Could not read number of linear nodes" );

        if ( verbose ){ fprintf(stderr, "Number of linear reduced bases for ROQ generation is %d\n", nbases0); }

        /* read in the number of nodes for this chunk */
        interp->nodes = XLALMalloc(nbases0*sizeof(UINT4));
        XLAL_CHECK_VOID( fread((void*)interp->nodes, sizeof(UINT4), nbases0, fpin) == nbases0, XLAL_EIO, "Could not read in linear interpolation indices" );
        nbases->data[i] = nbases0;
        dmlength += nbases0;
        nlineartot += nbases0;

        interpquad = XLALMalloc(sizeof(LALInferenceREALROQInterpolant));
        interpquad->B = NULL;

        /* read in the number of quadratic interpolant nodes for the current chunk and then the node indices */
        XLAL_CHECK_VOID( fread((void*)&nbases0quad, sizeof(UINT4), 1, fpin) == 1, XLAL_EIO, "Could not read number of linear nodes" );

        if ( verbose ){ fprintf(stderr, "Number of quadratic reduced bases for ROQ generation is %d\n", nbases0quad); }

        /* read in the number of nodes for this chunk */
        interpquad->nodes = XLALMalloc(nbases0quad*sizeof(UINT4));
        XLAL_CHECK_VOID( fread((void*)interpquad->nodes, sizeof(UINT4), nbases0quad, fpin) == nbases0quad, XLAL_EIO, "Could not read in quadratic interpolation indices" );
        nbasesquad->data[i] = nbases0quad;
        mmlength += nbases0quad;
        nquadtot += nbases0quad;
      }

      /* create vectors of time stamps at the nodes, the vector holds a "linear" set of nodes and then a "quadratic" set of nodes */
      /* add times at the linear interpolation nodes */
      UINT4 tnlength = 0;
      if ( !timenodes ){ tnlength = nbases0; }
      else { tnlength = timenodes->length + nbases0; }

      timenodes = XLALResizeTimestampVector( timenodes, tnlength );
      sidtimenodes = XLALResizeREAL8Vector( sidtimenodes, tnlength );
      timedatanodes = XLALResizeREAL8TimeSeries( timedatanodes, 0, tnlength );

      if ( ssbdelays != NULL ){ ssbnodes = XLALResizeREAL8Vector( ssbnodes, tnlength ); }
      if ( bsbdelays != NULL ){ bsbnodes = XLALResizeREAL8Vector( bsbnodes, tnlength ); }

      for ( j=0; j<nbases0; j++ ){
        timenodes->data[dlen+j] = ifo->times->data[startidx+interp->nodes[j]];
        sidtimenodes->data[dlen+j] = sidDayFrac->data[startidx+interp->nodes[j]];
        timedatanodes->data->data[dlen+j] = ifo->timeData->data->data[startidx+interp->nodes[j]];

        if ( ssbdelays != NULL ){ ssbnodes->data[dlen+j] = ssbdelays->data[startidx+interp->nodes[j]]; }
        if ( bsbdelays != NULL ){ bsbnodes->data[dlen+j] = bsbdelays->data[startidx+interp->nodes[j]]; }
      }

      dlen += nbases0;
      tnlength = timenodes->length + nbases0quad;

      /* get times at the quadratic interpolation nodes */
      timenodes = XLALResizeTimestampVector( timenodes, tnlength );
      sidtimenodes = XLALResizeREAL8Vector( sidtimenodes, tnlength );
      timedatanodes = XLALResizeREAL8TimeSeries( timedatanodes, 0, tnlength );

      if ( ssbdelays != NULL ){ ssbnodes = XLALResizeREAL8Vector( ssbnodes, tnlength ); }
      if ( bsbdelays != NULL ){ bsbnodes = XLALResizeREAL8Vector( bsbnodes, tnlength ); }

      for ( j=0; j<nbases0quad; j++ ){
        timenodes->data[dlen+j] = ifo->times->data[startidx+interpquad->nodes[j]];
        sidtimenodes->data[dlen+j] = sidDayFrac->data[startidx+interpquad->nodes[j]];
        timedatanodes->data->data[dlen+j] = ifo->timeData->data->data[startidx+interpquad->nodes[j]];

        if ( ssbdelays != NULL ){ ssbnodes->data[dlen+j] = ssbdelays->data[startidx+interpquad->nodes[j]]; }
        if ( bsbdelays != NULL ){ bsbnodes->data[dlen+j] = bsbdelays->data[startidx+interpquad->nodes[j]]; }
      }

      LALInferenceRemoveCOMPLEXROQInterpolant( interp );
      LALInferenceRemoveREALROQInterpolant( interpquad );

      startidx += tlen;
      dlen += nbases0quad;
    }

    ndatatot += ifo->times->length;

    /* resize time vectors and model vectors, so they just contain values at the interpolation nodes */
    //REAL8Vector *timescopy = XLALCreateREAL8Vector( ifo->times->length );
    //for ( j=0; j<ifo->times->length; j++ ){ timescopy->data[j] = XLALGPSGetREAL8( &ifo->times->data[j] ); } /* save times as a REAL8Vector rather than a LIGOTimeGPSVector */
    //LALInferenceAddVariable( ifo->params, "timeStampVectorFull", &timescopy, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
    LIGOTimeGPSVector *timescopy = XLALCreateTimestampVector( ifo->times->length );
    memcpy(timescopy->data, ifo->times->data, sizeof(LIGOTimeGPS)*ifo->times->length);
    LALInferenceAddVariable( ifo->params, "timeStampVectorFull", &timescopy, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED );
    XLALDestroyTimestampVector(ifo->times);
    ifo->times = timenodes;

    /* make a copy of the original sidereal time vectors and time stamp vectors - this is needed when calculating the SNR of the maximum likelihood point */
    REAL8Vector *siddaycopy = XLALCreateREAL8Vector( sidDayFrac->length );
    memcpy(siddaycopy->data, sidDayFrac->data, sizeof(REAL8)*sidDayFrac->length);
    LALInferenceRemoveVariable( ifo->params, "siderealDay" );
    LALInferenceAddVariable( ifo->params, "siderealDay", &sidtimenodes, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
    LALInferenceAddVariable( ifo->params, "siderealDayFull", &siddaycopy, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );

    //REAL8Vector *timedatacopy = XLALCreateREAL8Vector( ifo->timeData->data->length ); /* save times as a REAL8Vector rather than a REAL8TimesSeries */
    //memcpy(timedatacopy->data, ifo->timeData->data->data, sizeof(REAL8)*ifo->timeData->data->length);
    //LALInferenceAddVariable( ifo->params, "timeDataFull", &timedatacopy, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
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

    ifo->compTimeSignal = XLALResizeCOMPLEX16TimeSeries( ifo->compTimeSignal, 0, dmlength+mmlength );

    if ( inputroq ){
      /* now read in all the weights */
      dmweights = XLALMalloc( sizeof(COMPLEX16)*dmlength );
      XLAL_CHECK_VOID( fread((void*)dmweights, sizeof(COMPLEX16), dmlength, fpin) == dmlength, XLAL_EIO, "Could not read in data-model product weights" );

      mmweights = XLALMalloc( sizeof(REAL8)*mmlength );
      XLAL_CHECK_VOID( fread((void*)mmweights, sizeof(REAL8), mmlength, fpin) == mmlength, XLAL_EIO, "Could not read in model-model product weights" );
    }

    if ( outputroq ){
      XLAL_CHECK_VOID( fwrite(dmweights, sizeof(COMPLEX16), dmlength, fpout) == dmlength, XLAL_EIO, "Could not output data-model product weights" );
      XLAL_CHECK_VOID( fwrite(mmweights, sizeof(REAL8), mmlength, fpout) == mmlength, XLAL_EIO, "Could not output model-model product weights" );
    }

    /* fill in data/model weights into roq->weights complex matrix */
    data->roq->weightsLinear = XLALMalloc( dmlength*sizeof(COMPLEX16) );
    memcpy(data->roq->weightsLinear, dmweights, dmlength*sizeof(COMPLEX16));
    XLALFree( dmweights );

    data->roq->weightsQuadratic = XLALMalloc( mmlength*sizeof(REAL8) );
    memcpy(data->roq->weightsQuadratic, mmweights, mmlength*sizeof(REAL8));
    XLALFree( mmweights );

    /* add interpolation nodes to a variable in runState->threads[0]->model->ifo->params */
    LALInferenceAddVariable( ifo->params, "numBases", &nbases, LALINFERENCE_UINT4Vector_t, LALINFERENCE_PARAM_FIXED );
    LALInferenceAddVariable( ifo->params, "numBasesQuad", &nbasesquad, LALINFERENCE_UINT4Vector_t, LALINFERENCE_PARAM_FIXED );

    ifo = ifo->next;
    data = data->next;
    counter++;
  }

  if( verbose ){ fprintf(stderr, "Total number of linear bases = %d, total number of quadratic bases = %d, total number of data points = %d\n", nlineartot, nquadtot, ndatatot); }

  /* reset freqfactors to the correct value */
  ifo = runState->threads[0]->model->ifo;
  while ( ifo ){
    check_and_add_fixed_variable( ifo->params, "freqfactors", &freqsCopy, LALINFERENCE_REAL8Vector_t );
    ifo = ifo->next;
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
    if ( LALInferenceCheckVariable( runState->algorithmParams, "timefile" ) ){
      fclose(*(FILE **)LALInferenceGetVariable( runState->algorithmParams, "timefile" ));
    }
    fclose(fpout);
    if ( verbose ){ fprintf(stdout, "ROQ weights have been written to file.\nExiting program.\n"); }
    exit(0); /* exit the programme succussfully */
  }
}


/**
 * \brief Generate a training set of waveforms for the basis generation
 *
 * This function will create a set of \a n waveforms randomly placed over the
 * prior parameter space.
 *
 * @param[in] rs A temporary run state
 * @param[in] n The number of training waveforms
 *
 * @return A complex matrix containing an array of training waveforms
 */
COMPLEX16Array *generate_training_set( LALInferenceRunState *rs, UINT4 n ){
  UINT4 j = 0;
  UINT4Vector *dims = XLALCreateUINT4Vector( 2 );
  dims->data[0] = n;
  dims->data[1] = rs->threads[0]->model->ifo->times->length;
  COMPLEX16Array *ts = XLALCreateCOMPLEX16Array( dims );
  gsl_matrix_complex_view tsview = gsl_matrix_complex_view_array( (double *)ts->data, dims->data[0], dims->data[1] );
  XLALDestroyUINT4Vector( dims );

  rs->threads[0]->model->params = XLALCalloc(1, sizeof(LALInferenceVariables));

  /* for parameters with Gaussian priors reset them to be flat, with ranges over +/-5 sigma for training set generation if required */
  LALInferenceVariables *priorcopy = NULL;

  if( LALInferenceGetProcParamVal( rs->commandLine, "--roq-uniform" ) ){
    priorcopy = XLALCalloc(1, sizeof(LALInferenceVariables));
    LALInferenceCopyVariables( rs->priorArgs, priorcopy );
    LALInferenceVariableItem *item = rs->threads[0]->currentParams->head;
    UINT4 hascorrelated = 0; /* set if there is a correlated prior */
    for( ; item; item = item->next ){
      REAL8 mu = 0., sigma = 0., minval = 0., maxval = 0.;
      if ( LALInferenceCheckGaussianPrior( rs->priorArgs, item->name ) ){
        LALInferenceGetGaussianPrior( rs->priorArgs, item->name, &mu, &sigma );
        minval = mu - 5.*sigma;
        maxval = mu + 5.*sigma;
        LALInferenceRemoveGaussianPrior( rs->priorArgs, item->name ); /* remove prior */
      }
      else if( LALInferenceCheckCorrelatedPrior( rs->priorArgs, item->name ) ){
        gsl_matrix *cor = NULL, *invcor = NULL;
        UINT4 idx = 0;
        LALInferenceGetCorrelatedPrior( rs->priorArgs, item->name, &cor, &invcor, &mu, &sigma, &idx );
        minval = mu - 5.*sigma;
        maxval = mu + 5.*sigma;
        hascorrelated = 1;
      }
      else{ continue; }

      /* re-add prior as a uniform prior */
      LALInferenceAddMinMaxPrior( rs->priorArgs, item->name, &minval, &maxval, LALINFERENCE_REAL8_t );
    }
    if ( hascorrelated ) { LALInferenceRemoveCorrelatedPrior( rs->priorArgs ); } /* remove correlated prior */
  }

  /* choose random variables values and fill in runState->threads[0]->model->params */
  for ( j=0; j<n; j++ ){
    /* copy values to model parameters */
    LALInferenceCopyVariables(rs->threads[0]->currentParams, rs->threads[0]->model->params);

    /* draw from the prior distributions */
    LALInferenceDrawFromPrior( rs->threads[0]->model->params, rs->priorArgs, rs->GSLrandom );

    /* generate model */
    rs->threads[0]->model->templt( rs->threads[0]->model );

    /* place model into an array */
    gsl_vector_complex_view cview;
    cview = gsl_vector_complex_view_array((double*)rs->threads[0]->model->ifo->compTimeSignal->data->data, rs->threads[0]->model->ifo->times->length);
    gsl_matrix_complex_set_row(&tsview.matrix, j, &cview.vector);
  }

  LALInferenceClearVariables(rs->threads[0]->model->params);

  /* reset priors if required */
  if( LALInferenceGetProcParamVal( rs->commandLine, "--roq-uniform" ) ){
    LALInferenceClearVariables( rs->priorArgs );
    LALInferenceCopyVariables( priorcopy, rs->priorArgs );
    LALInferenceClearVariables( priorcopy );
  }

  return ts;
}


/**
 * \brief Generate a training set of waveforms for the quadratic term basis generation
 *
 * This function will create a set of \a n waveforms randomly placed over the
 * prior parameter space. The returned array will contain the real part of the product
 * of the waveform with its complex conjugate.
 *
 * @param[in] rs A temporary run state
 * @param[in] n The number of training waveforms
 *
 * @return A real matrix containing an array of training waveforms
 */
REAL8Array *generate_training_set_quad( LALInferenceRunState *rs, UINT4 n ){
  UINT4 j = 0, k = 0;
  UINT4Vector *dims = XLALCreateUINT4Vector( 2 );
  dims->data[0] = n;
  dims->data[1] = rs->threads[0]->model->ifo->times->length;
  REAL8Array *ts = XLALCreateREAL8Array( dims );
  gsl_matrix_view tsview = gsl_matrix_view_array( ts->data, dims->data[0], dims->data[1] );
  XLALDestroyUINT4Vector( dims );

  rs->threads[0]->model->params = XLALCalloc(1, sizeof(LALInferenceVariables));

  /* for parameters with Gaussian priors reset them to be flat, with ranges over +/-5 sigma for training set generation if required */
  LALInferenceVariables *priorcopy = NULL;
  if( LALInferenceGetProcParamVal( rs->commandLine, "--roq-uniform" ) ){
    priorcopy = XLALCalloc(1, sizeof(LALInferenceVariables));
    LALInferenceCopyVariables( rs->priorArgs, priorcopy );
    LALInferenceVariableItem *item = rs->threads[0]->currentParams->head;
    UINT4 hascorrelated = 0; /* set if there is a correlated prior */
    for( ; item; item = item->next ){
      REAL8 mu = 0., sigma = 0., minval = 0., maxval = 0.;
      if ( LALInferenceCheckGaussianPrior( rs->priorArgs, item->name ) ){
        LALInferenceGetGaussianPrior( rs->priorArgs, item->name, &mu, &sigma );
        minval = mu - 5.*sigma;
        maxval = mu + 5.*sigma;
        LALInferenceRemoveGaussianPrior( rs->priorArgs, item->name ); /* remove prior */
      }
      else if( LALInferenceCheckCorrelatedPrior( rs->priorArgs, item->name ) ){
        gsl_matrix *cor = NULL, *invcor = NULL;
        UINT4 idx = 0;
        LALInferenceGetCorrelatedPrior( rs->priorArgs, item->name, &cor, &invcor, &mu, &sigma, &idx );
        minval = mu - 5.*sigma;
        maxval = mu + 5.*sigma;
        hascorrelated = 1;
      }
      else{ continue; }

      /* re-add prior as a uniform prior */
      LALInferenceAddMinMaxPrior( rs->priorArgs, item->name, &minval, &maxval, LALINFERENCE_REAL8_t );
    }
    if ( hascorrelated ) { LALInferenceRemoveCorrelatedPrior( rs->priorArgs ); } /* remove correlated prior */
  }

  /* choose random variables values and fill in runState->threads[0]->model->params */
  for ( j=0; j<n; j++ ){
    /* copy values to model parameters */
    LALInferenceCopyVariables(rs->threads[0]->currentParams, rs->threads[0]->model->params);

    /* draw from the prior distributions */
    LALInferenceDrawFromPrior( rs->threads[0]->model->params, rs->priorArgs, rs->GSLrandom );

    /* generate model */
    rs->threads[0]->model->templt( rs->threads[0]->model );

    /* place model*model^* into an array */
    for ( k = 0; k < rs->threads[0]->model->ifo->times->length; k++ ){
      COMPLEX16 cval = rs->threads[0]->model->ifo->compTimeSignal->data->data[k];
      gsl_matrix_set(&tsview.matrix, j, k, creal(cval*conj(cval)));
    }
  }

  LALInferenceClearVariables(rs->threads[0]->model->params);

  /* reset priors if required */
  if( LALInferenceGetProcParamVal( rs->commandLine, "--roq-uniform" ) ){
    LALInferenceClearVariables( rs->priorArgs );
    LALInferenceCopyVariables( priorcopy, rs->priorArgs );
    LALInferenceClearVariables( priorcopy );
  }

  return ts;
}


/*--------------- END OF ROQ FUNCTIONS ----------------------*/
