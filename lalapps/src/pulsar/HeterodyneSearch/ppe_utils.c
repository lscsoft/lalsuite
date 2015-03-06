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
/*                            HELPER FUNCTIONS                                */
/******************************************************************************/

#include "ppe_utils.h"

/** \brief Compute the noise variance for each data segment
 *
 * Once the data has been split into segments calculate the noise variance (using
 * both the real and imaginary parts) in each segment and fill in the associated
 * noise vector. To calculate the noise the running median should first be
 * subtracted from the data.
 *
 * \param data [in] the LALInferenceIFOData variable
 * \param model [in] the LALInferenceIFOModel variable
 */
void compute_variance( LALInferenceIFOData *data, LALInferenceIFOModel *model ){
  REAL8 chunkLength = 0.;

  INT4 i = 0, j = 0, length = 0, cl = 0, counter = 0;

  COMPLEX16Vector *meddata = NULL; /* data with running median removed */

  /* subtract a running median value from the data to remove any underlying
     trends (e.g. caused by a strong signal) */
  meddata = subtract_running_median( data->compTimeData->data );

  UINT4Vector *chunkLengths = NULL;
  chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( model->params, "chunkLength" );

  length = data->compTimeData->data->length;

  for ( i = 0, counter = 0; i < length; i+=chunkLength, counter++ ){
    REAL8 vari = 0., meani = 0.;

    chunkLength = (REAL8)chunkLengths->data[counter];
    cl = i + (INT4)chunkLength;

    /* get the mean (should be close to zero given the running median subtraction), but
     * probably worth doing anyway) */
    for ( j = i ; j < cl ; j++ ){
      meani += (creal(meddata->data[j]) + cimag(meddata->data[j]));
    }

    meani /= (2.*chunkLength);

    for ( j = i ; j < cl ; j++ ){
      vari += SQUARE( (creal(meddata->data[j]) - meani) );
      vari += SQUARE( (cimag(meddata->data[j]) - meani) );
    }

    vari /= (2.*chunkLength - 1.); /* unbiased sample variance */

    /* fill in variance vector */
    for ( j = i ; j < cl ; j++ ){ data->varTimeData->data->data[j] = vari; }
  }
}


/**
 * \brief Split the data into segments
 *
 * This function is deprecated to \c chop_n_merge, but gives the functionality of the old code.
 *
 * It cuts the data into as many contiguous segments of data as possible of length \c chunkMax. Where contiguous is
 * defined as containing consecutive point within 180 seconds of each other. The length of segments that do not fit into
 * a \c chunkMax length are also included.
 *
 * \param ifo [in] the LALInferenceIFOModel variable
 * \param chunkMax [in] the maximum length of a data chunk/segment
 *
 * \return A vector of chunk/segment lengths
 */
UINT4Vector *get_chunk_lengths( LALInferenceIFOModel *ifo, INT4 chunkMax ){
  INT4 i = 0, j = 0, count = 0;
  INT4 length;

  REAL8 t1, t2;

  UINT4Vector *chunkLengths = NULL;

  length = ifo->times->length;

  chunkLengths = XLALCreateUINT4Vector( length );

  REAL8 dt = *(REAL8*)LALInferenceGetVariable( ifo->params, "dt" );

  /* create vector of data segment length */
  while( 1 ){
    count++; /* counter */

    /* break clause */
    if( i > length - 2 ){
      /* set final value of chunkLength */
      chunkLengths->data[j] = count;
      j++;
      break;
    }

    i++;

    t1 = XLALGPSGetREAL8( &ifo->times->data[i-1 ] );
    t2 = XLALGPSGetREAL8( &ifo->times->data[i] );

    /* if consecutive points are within two sample times of each other count as in the same chunk */
    if( t2 - t1 > 2.*dt || count == chunkMax ){
      chunkLengths->data[j] = count;
      count = 0; /* reset counter */

      j++;
    }
  }

  chunkLengths = XLALResizeUINT4Vector( chunkLengths, j );

  return chunkLengths;
}


/* function to use change point analysis to chop up and remerge the data to find stationary chunks (i.e. lengths of data
 * which look like they have the same statistics e.g. the same standard deviation) */
/**
 * \brief Chops and remerges data into stationary segments
 *
 * This function finds segments of data that appear to be stationary (have the same standard deviation).
 *
 * The function first attempts to chop up the data into as many stationary segments as possible. The splitting may not
 * be optimal, so it then tries remerging consecutive segments to see if the merged segments show more evidence of
 * stationarity. <b>[NOTE: Remerging is currently turned off and will make very little difference to the algorithm]</b>.
 * It then, if necessary, chops the segments again to make sure there are none greater than the required \c chunkMax.
 * The default \c chunkMax is 0, so this rechopping will not normally happen.
 *
 * This is all performed on data that has had a running median subtracted, to try and removed any underlying trends in
 * the data (e.g. those caused by a strong signal), which might affect the calculations (which assume the data is
 * Gaussian with zero mean).
 *
 * If the \c verbose flag is set then a list of the segments will be output to a file called \c data_segment_list.txt,
 * with a prefix of the detector name.
 *
 * \param data [in] A data structure
 * \param chunkMin [in] The minimum length of a segment
 * \param chunkMax [in] The maximum length of a segment
 *
 * \return A vector of segment/chunk lengths
 *
 * \sa subtract_running_median
 * \sa chop_data
 * \sa merge_data
 * \sa rechop_data
 */
UINT4Vector *chop_n_merge( LALInferenceIFOData *data, INT4 chunkMin, INT4 chunkMax ){
  UINT4 j = 0;

  UINT4Vector *chunkLengths = NULL;
  UINT4Vector *chunkIndex = NULL;

  COMPLEX16Vector *meddata = NULL;

  /* subtract a running median value from the data to remove any underlying trends (e.g. caused by a string signal) that
   * might affect the chunk calculations (which can assume the data is Gaussian with zero mean). */
  meddata = subtract_running_median( data->compTimeData->data );

  /* pass chop data a gsl_vector_view, so that internally it can use vector views rather than having to create new vectors */
  gsl_vector_complex_view meddatagsl = gsl_vector_complex_view_array((double*)meddata->data, meddata->length);
  chunkIndex = chop_data( &meddatagsl.vector, chunkMin );

  /* DON'T BOTHER WITH THE MERGING AS IT WILL MAKE VERY LITTLE DIFFERENCE */
  /* merge_data( meddata, chunkIndex ); */

  /* if a maximum chunk length is defined then rechop up the data, to segment any chunks longer than this value */
  if ( chunkMax > chunkMin ) { rechop_data( chunkIndex, chunkMax, chunkMin ); }

  chunkLengths = XLALCreateUINT4Vector( chunkIndex->length );

  /* go through segments and turn into vector of chunk lengths */
  for ( j = 0; j < chunkIndex->length; j++ ){
    if ( j == 0 ) { chunkLengths->data[j] = chunkIndex->data[j]; }
    else { chunkLengths->data[j] = chunkIndex->data[j] - chunkIndex->data[j-1]; }
  }

  /* if verbose print out the segment end indices to a file */
  if ( verbose_output ){
    FILE *fpsegs = NULL;

    CHAR *outfile = NULL;

    /* set detector name as prefix */
    outfile = XLALStringDuplicate( data->detector->frDetector.prefix );

    outfile = XLALStringAppend( outfile, "data_segment_list.txt" );

    if ( (fpsegs = fopen(outfile, "w")) == NULL ){
      fprintf(stderr, "Non-fatal error open file to output segment list.\n");
      return chunkLengths;
    }

    for ( j = 0; j < chunkIndex->length; j++ ) { fprintf(fpsegs, "%u\n", chunkIndex->data[j]); }

    /* add space at the end so that you can separate lists from different detector data streams */
    fprintf(fpsegs, "\n");

    fclose( fpsegs );
  }

  return chunkLengths;
}


/**
 * \brief Subtract the running median from complex data
 *
 * This function uses \c gsl_stats_median_from_sorted_data to subtract a running median, calculated from the 30
 * consecutive point around a set point, from the data. At the start of the data running median is calculated from
 * 30-15+(i-1) points, and at the end it is calculated from 15+(N-i) points, where i is the point index and N is the
 * total number of data points.
 *
 * \param data [in] A complex data vector
 *
 * \return A complex vector containing data with the running median removed
 */
COMPLEX16Vector *subtract_running_median( COMPLEX16Vector *data ){
  COMPLEX16Vector *submed = NULL;
  UINT4 length = data->length, i = 0, j = 0, n = 0;
  UINT4 mrange = 0;
  UINT4 N = 0;
  INT4 sidx = 0;

  if ( length > 30 ){ mrange = 30; } /* perform running median with 30 data points */
  else { mrange = 2*floor((length-1)/2); } /* an even number less than length */
  N = (UINT4)floor(mrange/2);

  submed = XLALCreateCOMPLEX16Vector( length );

  for ( i = 1; i < length+1; i++ ){
    REAL8 *dre = NULL;
    REAL8 *dim = NULL;

    /* get median of data within RANGE */
    if ( i < N ){
      n = N+i;
      sidx = 0;
    }
    else{
      if ( i > length - N ){ n = length - i + N; }
      else{ n = mrange; }

      sidx = i-N;
    }

    dre = XLALCalloc( n, sizeof(REAL8) );
    dim = XLALCalloc( n, sizeof(REAL8) );

    for ( j = 0; j < n; j++ ){
      dre[j] = creal(data->data[j+sidx]);
      dim[j] = cimag(data->data[j+sidx]);
    }

    /* sort data */
    gsl_sort( dre, 1, n );
    gsl_sort( dim, 1, n );

    /* get median and subtract from data*/
    submed->data[i-1] = ( creal(data->data[i-1]) - gsl_stats_median_from_sorted_data( dre, 1, n ) )
      + I * ( cimag(data->data[i-1]) - gsl_stats_median_from_sorted_data( dim, 1, n ) );

    XLALFree( dre );
    XLALFree( dim );
  }

  return submed;
}


/**
 * \brief Chops the data into stationary segments based on Bayesian change point analysis
 *
 * This function splits data into two (and recursively runs on those two segments) if it is found that the odds ratio
 * for them being from two independent Gaussian distributions is greater than a certain threshold.
 *
 * The threshold for the natural logarithm of the odds ratio is empirically set to be
 * \f[
 * T = 4.07 + 1.33\log{}_{10}{N},
 * \f]
 * where \f$N\f$ is the length in samples of the dataset. This is based on Monte Carlo simulations of
 * many realisations of Gaussian noise for data of different lengths. The threshold comes from a linear
 * fit to the log odds ratios required to give a 1% chance of splitting Gaussian data (drawn from a single
 * distribution) for data of various lengths.  Note, however, that this relation is not good for stretches of data
 * with lengths of less than about 30 points, and in fact is rather consevative for such short stretches
 * of data, i.e. such short stretches of data will require relatively larger odds ratios for splitting than
 * longer stretches.
 *
 * \param data [in] A complex data vector
 * \param chunkMin [in] The minimum allowed segment length
 *
 * \return A vector of segment lengths
 *
 * \sa find_change_point
 */
UINT4Vector *chop_data( gsl_vector_complex *data, INT4 chunkMin ){
  UINT4Vector *chunkIndex = NULL;

  UINT4 length = (UINT4)data->size;

  REAL8 logodds = 0.;
  UINT4 changepoint = 0;

  REAL8 threshold = 0.; /* may need tuning or setting globally */

  chunkIndex = XLALCreateUINT4Vector( 1 );

  changepoint = find_change_point( data, &logodds, chunkMin );

  /* threshold scaling for a 0.5% false alarm probability of splitting Gaussian data */
  threshold = 4.07 + 1.33*log10((REAL8)length);

  if ( logodds > threshold ){
    UINT4Vector *cp1 = NULL;
    UINT4Vector *cp2 = NULL;

    gsl_vector_complex_view data1 = gsl_vector_complex_subvector( data, 0, changepoint );
    gsl_vector_complex_view data2 = gsl_vector_complex_subvector( data, changepoint, length-changepoint );

    UINT4 i = 0, l = 0;

    cp1 = chop_data( &data1.vector, chunkMin );
    cp2 = chop_data( &data2.vector, chunkMin );

    l = cp1->length + cp2->length;

    chunkIndex = XLALResizeUINT4Vector( chunkIndex, l );

    /* combine new chunks */
    for (i = 0; i < cp1->length; i++) { chunkIndex->data[i] = cp1->data[i]; }
    for (i = 0; i < cp2->length; i++) { chunkIndex->data[i+cp1->length] = cp2->data[i] + changepoint; }
  }
  else{ chunkIndex->data[0] = length; }

  return chunkIndex;
}


/**
 * \brief Find a change point in complex data
 *
 * This function is based in the Bayesian Blocks algorithm of \cite Scargle1998 that finds "change points" in data -
 * points at which the statistics of the data change. It is based on calculating evidence, or odds, ratios. The
 * function first computes the marginal likelihood (or evidence) that the whole of the data is described by a single
 * Gaussian (with mean of zero). This comes from taking a Gaussian likelihood function and analytically marginalising
 * over the standard deviation (using a prior on the standard deviation of \f$1/\sigma\f$), giving (see
 * [\cite DupuisWoan2005]) a Students-t distribution (see
 * <a href="https://wiki.ligo.org/foswiki/pub/CW/PulsarParameterEstimationNestedSampling/studentst.pdf">here</a>).
 * Following this the data is split into two segments (with lengths greater than, or equal to the minimum chunk length)
 * for all possible combinations, and the joint evidence for each of the two segments consisting of independent
 * Gaussian (basically multiplying the above equation calculated for each segment separately) is calculated and the
 * split point recorded. However, the value required for comparing to that for the whole data set, to give the odds
 * ratio, is the evidence that having any split is better than having no split, so the individual split data evidences
 * need to be added incoherently to give the total evidence for a split. The index at which the evidence for a single
 * split is maximum (i.e. the most favoured split point) that which is returned.
 *
 * \param data [in] a complex data vector
 * \param logodds [in] a pointer to return the natural logarithm of the odds ratio/Bayes factor
 * \param minlength [in] the minimum chunk length
 *
 * \return The position of the change point
 */
UINT4 find_change_point( gsl_vector_complex *data, REAL8 *logodds, INT4 minlength ){
  UINT4 changepoint = 0, i = 0;
  UINT4 length = (UINT4)data->size, lsum = 0;

  REAL8 datasum = 0.;

  REAL8 logsingle = 0., logtot = -INFINITY;
  REAL8 logdouble = 0., logdouble_min = -INFINITY;
  REAL8 logratio = 0.;

  REAL8 sumforward = 0., sumback = 0.;
  gsl_complex dval;

  /* check that data is at least twice the minimum length, if not return an odds ratio of zero (log odds = -inf [or close to that!]) */
  if ( length < (UINT4)(2*minlength) ){
    logratio = -INFINITY;
    memcpy(logodds, &logratio, sizeof(REAL8));
    return 0;
  }

  /* calculate the sum of the data squared */
  for (i = 0; i < length; i++) {
    dval = gsl_vector_complex_get( data, i );
    datasum += SQUARE( gsl_complex_abs( dval ) );
  }

  /* calculate the evidence that the data consists of a Gaussian data with a single standard deviation */
  logsingle = -LAL_LN2 - (REAL8)length*LAL_LNPI + gsl_sf_lnfact(length-1) - (REAL8)length * log( datasum );

  lsum = length - 2*minlength + 1;

  for ( i = 0; i < length; i++ ){
    dval = gsl_vector_complex_get( data, i );
    if ( i < (UINT4)minlength-1 ){ sumforward += SQUARE( gsl_complex_abs( dval ) ); }
    else{ sumback += SQUARE( gsl_complex_abs( dval ) ); }
  }

  /* go through each possible change point and calculate the evidence for the data consisting of two independent
   * Gaussian's either side of the change point. Also calculate the total evidence for any change point.
   * Don't allow single points, so start at the second data point. */
  for (i = 0; i < lsum; i++){
    UINT4 ln1 = i+minlength, ln2 = (length-i-minlength);

    REAL8 log_1 = 0., log_2 = 0.;

    dval = gsl_vector_complex_get( data, ln1-1 );
    REAL8 adval = SQUARE( gsl_complex_abs( dval ) );
    sumforward += adval;
    sumback -= adval;

    /* get log evidences for the individual segments */
    log_1 = -LAL_LN2 - (REAL8)ln1*LAL_LNPI + gsl_sf_lnfact(ln1-1) - (REAL8)ln1 * log( sumforward );
    log_2 = -LAL_LN2 - (REAL8)ln2*LAL_LNPI + gsl_sf_lnfact(ln2-1) - (REAL8)ln2 * log( sumback );

    /* get evidence for the two segments */
    logdouble = log_1 + log_2;

    /* add to total evidence for a change point */
    logtot = LOGPLUS(logtot, logdouble);

    /* find maximum value of logdouble and record that as the change point */
    if ( logdouble > logdouble_min ){
      changepoint = ln1;
      logdouble_min = logdouble;
    }
  }

  /* get the log odds ratio of segmented versus non-segmented model */
  logratio = logtot - logsingle;
  memcpy(logodds, &logratio, sizeof(REAL8));

  return changepoint;
}


/**
 * \brief Chop up the data into chunks smaller the the maximum allowed length
 *
 * This function chops any chunks that are greater than \c chunkMax into chunks smaller than, or equal to \c chunkMax,
 * and greater than \c chunkMin. On some occasions this might result in a segment smaller than \c chunkMin, but these
 * are ignored in the likelihood calculation anyway.
 *
 * \param chunkIndex [in] a vector of segment split positions
 * \param chunkMax [in] the maximum allowed segment/chunk length
 * \param chunkMin [in] the minimum allowed segment/chunk length
 */
void rechop_data( UINT4Vector *chunkIndex, INT4 chunkMax, INT4 chunkMin ){
  INT4 i = 0, j = 0, count = 0;
  INT4 length = chunkIndex->length;
  INT4 endIndex = (INT4)chunkIndex->data[length-1];
  UINT4 startindex = 0, chunklength = 0;

  UINT4Vector *newindex = NULL;
  newindex = XLALCreateUINT4Vector( ceil((REAL8)endIndex / (REAL8)chunkMax ) );

  /* chop any chunks that are greater than chunkMax into chunks smaller than, or equal to chunkMax, and greater than chunkMin */
  for ( i = 0; i < length; i++ ){
    if ( i == 0 ) { startindex = 0; }
    else { startindex = chunkIndex->data[i-1]+1; }

    chunklength = chunkIndex->data[i] - startindex;

    if ( chunklength > (UINT4)chunkMax ){
      INT4 remain = chunklength % chunkMax;

      /* cut segment into as many chunkMin chunks as possible */
      for ( j = 0; j < floor(chunklength / chunkMax); j++ ){
        newindex->data[count] = startindex + (j+1)*chunkMax;
        count++;
      }

      /* last chunk values */
      if ( remain != 0 ){
        /* set final value */
        newindex->data[count] = startindex + j*chunkMax + remain;

        if ( remain < chunkMin ){
          /* split the last two cells into one that is chunkMin long and one that is (chunkMax+remainder)-chunkMin long
           * - this may leave a cell shorter than chunkMin, but we'll have to live with that! */
          INT4 n1 = (chunkMax + remain) - chunkMin;

          /* reset second to last value two values */
          newindex->data[count-1] = newindex->data[count] - chunkMin;

          if ( n1 < chunkMin && verbose_output ){
            fprintf(stderr, "Non-fatal error... segment no. %d is %d long, which is less than chunkMin = %d.\n",
                    count, n1, chunkMin);
          }
        }

        count++;
      }
    }
    else{
      newindex->data[count] = chunkIndex->data[i];
      count++;
    }
  }

  chunkIndex = XLALResizeUINT4Vector( chunkIndex, count );

  for ( i = 0; i < count; i++ ) { chunkIndex->data[i] = newindex->data[i]; }

  XLALDestroyUINT4Vector( newindex );
}


/**
 * \brief Merge adjacent segments
 *
 * This function will attempt to remerge adjacent segments if statistically favourable (as calculated by the odds
 * ratio). For each pair of adjacent segments the joint likelihood of them being from two independent distributions is
 * compared to the likelihood that combined they are from one distribution. If the likelihood is highest for the
 * combined segments they are merged.
 *
 * \param data [in] A complex data vector
 * \param segs [in] A vector of split segment indexes
 */
void merge_data( COMPLEX16Vector *data, UINT4Vector *segs ){
  UINT4 j = 0;
  REAL8 threshold = 0.; /* may need to be passed to function in the future, or defined globally */

  /* loop until stopping criterion is reached */
  while( 1 ){
    UINT4 ncells = segs->length;

    UINT4 mergepoint = 0;
    REAL8 logodds = 0., minl = -LAL_REAL8_MAX;

    for (j = 1; j < ncells; j++){
      REAL8 summerged = 0., sum1 = 0., sum2 = 0.;
      UINT4 i = 0, n1 = 0, n2 = 0, nm = 0;
      UINT4 cellstarts1 = 0, cellends1 = 0, cellstarts2 = 0, cellends2 = 0;
      REAL8 log_merged = 0., log_individual = 0.;

      /* get the evidence for merged adjacent cells */
      if( j == 1 ) { cellstarts1 = 0; }
      else { cellstarts1 = segs->data[j-2]; }

      cellends1 = segs->data[j-1];

      cellstarts2 = segs->data[j-1];
      cellends2 = segs->data[j];

      n1 = cellends1 - cellstarts1;
      n2 = cellends2 - cellstarts2;
      nm = cellends2 - cellstarts1;

      for( i = cellstarts1; i < cellends1; i++ ) { sum1 += SQUARE( cabs(data->data[i]) ); }

      for( i = cellstarts2; i < cellends2; i++ ) { sum2 += SQUARE( cabs(data->data[i]) ); }

      summerged = sum1 + sum2;

      /* calculated evidences */
      log_merged = -2 + gsl_sf_lnfact(nm-1) - (REAL8)nm * log( summerged );

      log_individual = -2 + gsl_sf_lnfact(n1-1) - (REAL8)n1 * log( sum1 );
      log_individual += -2 + gsl_sf_lnfact(n2-1) - (REAL8)n2 * log( sum2 );

      logodds = log_merged - log_individual;

      if ( logodds > minl ){
        mergepoint = j - 1;
        minl = logodds;
      }
    }

    /* set break criterion */
    if ( minl < threshold ) { break; }
    else{ /* merge cells */
      /* remove the cell end value between the two being merged and shift */
      for( UINT4 i=0; i < ncells-(mergepoint+1); i++ ){
        segs->data[mergepoint+i] = segs->data[mergepoint+i+1];
      }

      segs = XLALResizeUINT4Vector( segs, ncells - 1 );
    }
  }
}


/**
 * \brief Gzip the nested sample files
 *
 * This function gzips the output nested sample files. It will also strip unneccesary parameters
 * from the header "_params.txt" file.
 *
 * \param runState [in] The analysis information structure
 */
void gzip_output( LALInferenceRunState *runState ){
  /* Open original output output file */
  CHAR *outfile = NULL;
  ProcessParamsTable *ppt1 = LALInferenceGetProcParamVal( runState->commandLine, "--outfile" );

  if( ppt1 ){
    CHAR outfilepars[256] = "", outfileparstmp[256] = "";
    FILE *fppars = NULL, *fpparstmp = NULL;
    UINT4 nonfixed = 0;

    outfile = ppt1->value;

    /* open file for printing out list of parameter names - this should already exist */
    sprintf(outfilepars, "%s_params.txt", outfile);
    if( (fppars = fopen(outfilepars, "r")) == NULL ){
      XLALPrintError("Error... cannot open parameter name output file %s.\n", outfilepars);
      XLAL_ERROR_VOID(XLAL_EIO);
    }
    /* read in the parameter names and remove the "model" value */
    sprintf(outfileparstmp, "%s_params.txt_tmp", outfile);
    if( (fpparstmp = fopen(outfileparstmp, "w")) == NULL ){
      XLALPrintError("Error... cannot open parameter name output file %s.\n", outfileparstmp);
      XLAL_ERROR_VOID(XLAL_EIO);
    }

    if ( LALInferenceGetProcParamVal( runState->commandLine, "--non-fixed-only" ) ){ nonfixed = 1; }

    CHAR v[128] = "";
    while( fscanf(fppars, "%s", v) != EOF ){
      /* if outputing only non-fixed values then only re-output names of those non-fixed things */
      if ( nonfixed ){
        if ( LALInferenceCheckVariable( runState->currentParams, v ) ){
          if ( LALInferenceGetVariableVaryType( runState->currentParams, v ) != LALINFERENCE_PARAM_FIXED ){
            fprintf(fpparstmp, "%s\t", v);
          }
        }
      }
      else{
        /* re-output everything to a temporary file */
        fprintf(fpparstmp, "%s\t", v);
      }
    }

    fclose(fppars);
    fclose(fpparstmp);

    /* move the temporary file name to the standard outfile_param name */
    rename( outfileparstmp, outfilepars );

    /* gzip the output file if required */
    if( LALInferenceGetProcParamVal( runState->commandLine, "--gzip" ) ){
      if ( XLALGzipTextFile( outfile ) != XLAL_SUCCESS ){
        XLAL_PRINT_ERROR( "Error... Could not gzip the output file!\n" );
        XLAL_ERROR_VOID( XLAL_EIO );
      }
    }
  }
/* if we have XML enabled */
#ifdef HAVE_LIBLALXML
  ProcessParamsTable *ppt2 = LALInferenceGetProcParamVal( runState->commandLine, "--outxml" );
  if(!ppt2){
    ppt2=LALInferenceGetProcParamVal(runState->commandLine,"--outXML");
  }
  CHAR *outVOTable = NULL;

  if ( !ppt2 && !ppt1 ){
    XLALPrintError("Must specify either --outfile or --outXML\n");
    XLAL_ERROR_VOID( XLAL_EIO );
  }

  if( ppt2 ){
    outVOTable = ppt2->value;

    /* gzip the output file if required */
    if( LALInferenceGetProcParamVal( runState->commandLine, "--gzip" ) ){
      if ( XLALGzipTextFile( outVOTable ) != XLAL_SUCCESS ){
        XLAL_PRINT_ERROR( "Error... Could not gzip the output file!\n" );
        XLAL_ERROR_VOID( XLAL_EIO );
      }
    }
  }

#else
  if ( !ppt1 ){
    XLALPrintError("Error... --outfile not defined!\n");
    XLAL_ERROR_VOID( XLAL_EIO );
  }
#endif

  return;
}


/**
 * \brief Counts the number of comma separated values in a string
 *
 * This function counts the number of comma separated values in a given input string.
 *
 * \param csvline [in] Any string
 *
 * \return The number of comma separated value in the input string
 */
INT4 count_csv( CHAR *csvline ){
  CHAR *inputstr = NULL;
  INT4 count = 0;

  inputstr = XLALStringDuplicate( csvline );

  /* count number of commas */
  while(1){
    if( XLALStringToken(&inputstr, ",", 0) == NULL ){
      XLALPrintError("Error... problem counting number of commas!\n");
      XLAL_ERROR( XLAL_EFUNC );
    }

    if ( inputstr == NULL ) { break; }

    count++;
  }

  return count+1;
}


/**
 * \brief Checks if a given parameter is recognised
 *
 * This function checks whether a given parameter is one of the defined amplitude (\c amppars), frequency (\c freqpars),
 * sky location (\c skypars) or binary system (\c binpars) parameters given in the header file.
 *
 * \param parname [in] The name of a parameter
 *
 * \return true (1) if the parameter is recognised and false (0) if not
 */
INT4 recognised_parameter( CHAR *parname ){
  INT4 i = 0;

  for( i = 0; i < NUMAMPPARS; i++ ) { if (!strcmp(parname, amppars[i])) { return 1; } }
  for( i = 0; i < NUMFREQPARS; i++ ) { if (!strcmp(parname, freqpars[i])) { return 1; } }
  for( i = 0; i < NUMSKYPARS; i++ ) { if (!strcmp(parname, skypars[i])) { return 1; } }
  for( i = 0; i < NUMBINPARS; i++ ) { if (!strcmp(parname, binpars[i])) { return 1; } }

  return 0;
}


/**
 * \brief Automatically set the solar system ephemeris file based on environment variables and data time span
 *
 * This function will attempt to find Earth and Sun ephemeris files based on LAL environment variables (as set up by
 * <code> lalpulsar-user-env.(c)sh </code>) and a given start and end GPS time (presumably taken from the data that is
 * to be analysed). It requires \c LALPULSAR is installed and the \c LALPULSAR_PREFIX variable is set, which should mean
 * that ephemeris files are installed in the directory \c ${LALPULSAR_PREFIX}/share/lalpulsar/.
 *
 * \param efile [in] a string that will return the Earth ephemeris file
 * \param sfile [in] a string that will return the Sun ephemeris file
 * \param tfile [in] a string that will return the time correction file
 * \param pulsar [in] the pulsar parameters read from a .par file
 * \param gpsstart [in] the GPS time of the start of the data
 * \param gpsend [in] the GPS time of the end of the data
 *
 * \return The TimeCorrectionType e.g. TDB or TCB
 */
TimeCorrectionType XLALAutoSetEphemerisFiles( CHAR *efile, CHAR *sfile, CHAR *tfile, BinaryPulsarParams  pulsar,
                                              INT4 gpsstart, INT4 gpsend ){
  /* set the times that the ephemeris files span */
  INT4 ephemstart = 630720013; /* GPS time of Jan 1, 2000, 00:00:00 UTC */
  INT4 ephemend = 1261872015; /* GPS time of Jan 1, 2020, 00:00:00 UTC */
  CHAR *eftmp = NULL, *sftmp = NULL, *tftmp = NULL;
  TimeCorrectionType ttype = TIMECORRECTION_NONE;

  CHAR *lalpath = NULL, *lalpulsarpath = NULL;

  if( gpsstart < ephemstart || gpsend < ephemstart || gpsstart > ephemend || gpsend > ephemend ){
    XLALPrintError("Start and end times are outside the ephemeris file ranges!\n");
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* first check that the path to the Ephemeris files is available in the
     environment variables */
  if((lalpath = getenv("LALPULSAR_PREFIX")) == NULL){
    XLALPrintError("LALPULSAR_PREFIX environment variable not set. Cannot automatically generate ephemeris files!\n");
    XLAL_ERROR(XLAL_EFUNC);
  }

  lalpulsarpath = XLALStringDuplicate( lalpath );

  if ( (lalpulsarpath = XLALStringAppend(lalpulsarpath, "/share/lalpulsar/")) == NULL ) { XLAL_ERROR(XLAL_EFUNC); }

  eftmp = XLALStringDuplicate(lalpulsarpath);
  sftmp = XLALStringDuplicate(lalpulsarpath);
  tftmp = XLALStringDuplicate(lalpulsarpath);

  eftmp = XLALStringAppend(eftmp, "earth00-19-");
  sftmp = XLALStringAppend(sftmp, "sun00-19-");

  if( pulsar.ephem == NULL ){
    /* default to use DE405 */
    eftmp = XLALStringAppend(eftmp, "DE405");
    sftmp = XLALStringAppend(sftmp, "DE405");
  }
  else{
    if( !strcmp(pulsar.ephem, "DE405") || !strcmp(pulsar.ephem, "DE200") ||
        !strcmp(pulsar.ephem, "DE414") ){
      eftmp = XLALStringAppend(eftmp, pulsar.ephem);
      sftmp = XLALStringAppend(sftmp, pulsar.ephem);
    }
    else{
      XLALPrintError("Unknown ephemeris %s in par file\n", pulsar.ephem);
      XLAL_ERROR(XLAL_EFUNC);
    }
  }

  /* add .dat extension */
  eftmp = XLALStringAppend(eftmp, ".dat.gz");
  sftmp = XLALStringAppend(sftmp, ".dat.gz");

  if ( eftmp == NULL || sftmp == NULL ) { XLAL_ERROR(XLAL_EFUNC); }

  XLALStringCopy( efile, eftmp, strlen(eftmp)+1 );
  XLALStringCopy( sfile, sftmp, strlen(sftmp)+1 );

  /* double check that the files exist */
  if( fopen(sfile, "r") == NULL || fopen(efile, "r") == NULL ){
    XLALPrintError("Error... ephemeris files not, or incorrectly, defined!\n");
    XLAL_ERROR(XLAL_EFUNC);
  }

  if( pulsar.units == NULL ){
    /* default to using TCB units */
    tftmp = XLALStringAppend(tftmp, "te405_2000-2019.dat.gz");
    ttype = TIMECORRECTION_TCB;
  }
  else{
    if ( !strcmp( pulsar.units, "TDB" ) ){
      tftmp = XLALStringAppend(tftmp, "tdb_2000-2019.dat.gz");
      ttype = TIMECORRECTION_TDB;
    }
    else{
      XLALPrintError("Error... unknown units %s in par file!\n", pulsar.units);
      XLAL_ERROR(XLAL_EFUNC);
    }
  }

  if ( tftmp == NULL ) { XLAL_ERROR(XLAL_EFUNC); }

  XLALStringCopy( tfile, tftmp, strlen(tftmp)+1 );

  if( fopen(tfile, "r") == NULL ){
    XLALPrintError("Error... time ephemeris files not, or incorrectly, defined!\n");
    XLAL_ERROR(XLAL_EFUNC);
  }

  return ttype;
}

/**
 * \brief Add a variable, checking first if it already exists and is of type \c LALINFERENCE_PARAM_FIXED
 * and if so removing it before re-adding it
 *
 * This function is for use as an alternative to \c LALInferenceAddVariable, which does not
 * allow \c LALINFERENCE_PARAM_FIXED variables to be changed. If the variable already exists
 * and is of type  \c LALINFERENCE_PARAM_FIXED, then it will first be removed and then re-added
 */
void check_and_add_fixed_variable( LALInferenceVariables *vars, const char *name, void *value, LALInferenceVariableType type ){
   if ( LALInferenceCheckVariable( vars, name ) ){
     if ( LALInferenceGetVariableVaryType( vars, name ) == LALINFERENCE_PARAM_FIXED ) { LALInferenceRemoveVariable( vars, name ); }
   }
   LALInferenceAddVariable( vars, name, value, type, LALINFERENCE_PARAM_FIXED );
}
