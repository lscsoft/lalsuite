#define LAL_USE_OLD_COMPLEX_STRUCTS
#include "coh_PTF.h"

/* gets the data, performs any injections, and conditions the data */
REAL4TimeSeries *coh_PTF_get_data(
              struct coh_PTF_params *params,
              const char *ifoChannel, 
              const char *dataCache,
              UINT4 ifoNumber  )
{
  int stripPad = 0;
  REAL4TimeSeries *channel = NULL;

  /* compute the start and duration needed to pad data */
  params->frameDataStartTime = params->startTime;
  XLALGPSAdd( &params->frameDataStartTime, -1.0 * params->padData );
  params->frameDataDuration = params->duration + 2.0 * params->padData;

  if ( params->getData )
  {
    if ( params->simData )
    {
      // First set the simulated PSD
      REAL8FrequencySeries *spectrum = NULL;
      /* Value of 1E-20 is used for white spectrum. It sets the PSD scale
         to give SNRs that are kind of comparable to iLIGO PSD */
      spectrum = generate_theoretical_psd(1./params->sampleRate,
          params->duration,params->simDataType, 1E-20);
      channel = get_simulated_data_new( ifoChannel, &params->startTime,
          params->duration, params->sampleRate,
          params->randomSeed+ 100*ifoNumber, spectrum );
      XLALDestroyREAL8FrequencySeries(spectrum);
    }
    else if ( params->zeroData )
    {
      channel = get_zero_data( ifoChannel, &params->startTime,
          params->duration, LALRINGDOWN_DATATYPE_ZERO, params->sampleRate );
    }
    else if ( (ifoNumber == LAL_IFO_V1 && params->virgoDoubleData)
              || (ifoNumber != LAL_IFO_V1 && params->ligoDoubleData) )
    {
      channel = get_frame_data_dbl_convert( dataCache, ifoChannel,
          &params->frameDataStartTime, params->frameDataDuration,
          LALRINGDOWN_DATATYPE_HT_REAL8,
          params->highpassFrequency);
      stripPad = 1;
    }
    else
    {
      channel = ring_get_frame_data( dataCache, ifoChannel,
          &params->frameDataStartTime, params->frameDataDuration,
          LALRINGDOWN_DATATYPE_HT_REAL4 );
      stripPad = 1;
    }
    if ( params->writeRawData ) /* write raw data */
      write_REAL4TimeSeries( channel );

    /* Function to put injections overhead */
    /*snprintf( channel->name, LALNameLength * sizeof(CHAR), "ZENITH" );*/

    /* inject signals */
    if ( params->injectFile )
      ring_inject_signal( channel, LALRINGDOWN_EOBNR_INJECT, params->injectFile,
          NULL, 1.0, NULL );
    if ( params->writeRawData )
       write_REAL4TimeSeries( channel );

    /* condition the data: resample and highpass */
    resample_REAL4TimeSeries( channel, params->sampleRate );
    if ( params->writeProcessedData ) /* write processed data */
      write_REAL4TimeSeries( channel );

    if (! params->zeroData )
    {
      highpass_REAL4TimeSeries( channel, params->highpassFrequency );
      if ( params->writeProcessedData ) /* write processed data */
        write_REAL4TimeSeries( channel );
    }

    if ( stripPad )
    {
      trimpad_REAL4TimeSeries( channel, params->padData );
      if ( params->writeProcessedData ) /* write data with padding removed */
        write_REAL4TimeSeries( channel );
    }
  }

  return channel;
}



/* This function is used to generate the null stream */
/* Be aware this is separate from the null SNR! */
int coh_PTF_get_null_stream(
    struct coh_PTF_params *params,
    REAL4TimeSeries *channel[LAL_NUM_IFO + 1],
    REAL4 *Fplus,
    REAL4 *Fcross,
    REAL4 *timeOffsets )
{
  UINT4 j,n[3],ifoNumber;
  INT4 i,t[3];
  INT4 timeOffsetPoints[LAL_NUM_IFO];
  REAL4 deltaT;
  REAL4TimeSeries *series;
  /* Determine how many detectors we are dealing with and assign 1-3 */
  i = -1;
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if ( params->haveTrig[ifoNumber] )
    {
      i++;
      if (i == 3)
      {
        fprintf(stderr,"Currently Null stream can only be constructed for three ifo combinations \n");
        return 1;
      }
      n[i] = ifoNumber;
    }
  }
  if (i < 2)
  {
    fprintf(stderr,"Null stream requires at least three detectors\n");
    return 1;
  }


  /* Determine time offset as number of data poitns */
  deltaT = channel[n[0]]->deltaT;
  for (i = 0; i < LAL_NUM_IFO; i++ )
  {
    timeOffsetPoints[i]=(int)(timeOffsets[i]/deltaT);
  }

  /* Initialise the null stream structures */
  series = LALCalloc( 1, sizeof( *series ) );
  series->data = XLALCreateREAL4Vector( channel[n[0]]->data->length );
  snprintf( series->name, sizeof( series->name ), "NULL_STREAM");
  series->epoch  = channel[n[0]]->epoch;
  series->deltaT = deltaT;
  series->sampleUnits = lalStrainUnit;
  for ( j = 0; j < series->data->length; ++j )
  {
    for ( i =0 ; i < 3; i++ )
    {
      t[i] = j + timeOffsetPoints[n[i]];
      if (t[i] < 0)
      {
        t[i] = t[i] + series->data->length;
      }
      if (t[i] > (INT4)(series->data->length-1))
      {
        t[i] = t[i] - series->data->length;
      }
    }
    series->data->data[j] = (Fplus[n[1]]*Fcross[n[2]] - Fplus[n[2]]*Fcross[n[1]]) * channel[n[0]]->data->data[t[0]];
    series->data->data[j] += (Fplus[n[2]]*Fcross[n[0]] - Fplus[n[0]]*Fcross[n[2]]) * channel[n[1]]->data->data[t[1]];
    series->data->data[j] += (Fplus[n[0]]*Fcross[n[1]] - Fplus[n[1]]*Fcross[n[0]]) * channel[n[2]]->data->data[t[2]];
  }
  if ( params->writeProcessedData ) /* write processed data */
    write_REAL4TimeSeries( series );
  channel[LAL_NUM_IFO] = series;

  return 0;
}

/* computes the inverse power spectrum */
REAL4FrequencySeries *coh_PTF_get_invspec(
    REAL4TimeSeries         *channel,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *revplan,
    struct coh_PTF_params   *params
    )
{
  REAL4FrequencySeries *invspec = NULL;
  if ( params->getSpectrum )
  {
    if ( ! params->whiteSpectrum)
    {
      /* compute raw average spectrum; store spectrum in invspec for now */
      invspec = compute_average_spectrum( channel,
          LALRINGDOWN_SPECTRUM_MEDIAN_MEAN, params->segmentDuration,
          params->strideDuration, fwdplan, 0 );
    }
    else
    {
      UINT4 k;
      REAL8FrequencySeries *spectrum;
      spectrum = generate_theoretical_psd(1./params->sampleRate,
          params->segmentDuration,params->simDataType, 1E-20);
  
      /* Need to convert to a REAL4 FrequencySeries */
      UINT4 segmentLength = floor( params->segmentDuration/channel->deltaT\
                                    + 0.5 );
      invspec = XLALCreateREAL4FrequencySeries("TEMP",&(channel->epoch),0,
          1.0/params->segmentDuration,&lalDimensionlessUnit,
          segmentLength/2 + 1);
      snprintf( invspec->name, sizeof( invspec->name),
          "%s_SPEC", channel->name);
      for ( k = 0; k < spectrum->data->length; ++k )
      {
        invspec->data->data[k] = (REAL4)(spectrum->data->data[k]*1E40);
      }
      XLALDestroyREAL8FrequencySeries(spectrum);
    }
    
    if ( params->writeInvSpectrum ) /* Write spectrum before inversion */
      write_REAL4FrequencySeries( invspec );

    /* invert spectrum */
    /* FIXME: What lower frequency should go here? */
    invert_spectrum( invspec, params->sampleRate, params->strideDuration,
        params->truncateDuration, params->lowTemplateFrequency, fwdplan,
        revplan );

    if ( params->writeInvSpectrum ) /* write inverse calibrated spectrum */
      write_REAL4FrequencySeries( invspec );
  }

  return invspec;
}

/* Function to rescale the data to avoid floating point errors*/
void coh_PTF_rescale_data (REAL4TimeSeries *channel,REAL8 rescaleFactor)
{
  UINT4 k;
  for ( k = 0; k < channel->data->length; ++k )
  {
    channel->data->data[k] *= rescaleFactor;
  }
}

/* creates the requested data segments (those in the list of segments to do) */
RingDataSegments *coh_PTF_get_segments(
    REAL4TimeSeries         *channel,
    REAL4FrequencySeries    *invspec,
    REAL4FFTPlan            *fwdplan,
    InterferometerNumber     NumberIFO,
    REAL4                   *timeSlideVectors,
    struct coh_PTF_params   *params
    )
{
  RingDataSegments *segments = NULL;
  COMPLEX8FrequencySeries  *response = NULL;
  UINT4  sgmnt,i,j, slidSegNum;
  UINT4  segListToDo[params->numOverlapSegments];

  segments = LALCalloc( 1, sizeof( *segments ) );

  if ( params->analyzeInjSegsOnly )
  {
    for ( i = 0 ; i < params->numOverlapSegments; i++)
      segListToDo[i] = 0;
    SimInspiralTable        *injectList = NULL;
    REAL8 deltaTime,segBoundDiff;
    INT4 segNumber, UNUSED segLoc, UNUSED ninj;
    segLoc = 0;
    if (! params->injectList)
    {
      ninj = SimInspiralTableFromLIGOLw( &injectList, params->injectFile, params->startTime.gpsSeconds, params->startTime.gpsSeconds + params->duration );
      params->injectList = injectList;
    }
    else
    {
      injectList = params->injectList;
    }
    while (injectList)
    {
      deltaTime = injectList->geocent_end_time.gpsSeconds;
      deltaTime += injectList->geocent_end_time.gpsNanoSeconds * 1E-9;
      deltaTime -= params->startTime.gpsSeconds;
      deltaTime -= params->startTime.gpsNanoSeconds * 1E-9;
      segNumber = floor(2*(deltaTime/params->segmentDuration) - 0.5);
      segListToDo[segNumber] = 1;
      /* Check if injection is near a segment boundary */
      for ( j = 0 ; j < params->numOverlapSegments; j++)
      {
        segBoundDiff = deltaTime - (j+0.5) * params->segmentDuration/2;
        if (segBoundDiff > 0 && segBoundDiff < params->injSearchWindow)
        {
          if (j != 0)
          {
            segListToDo[segNumber-1] = 1;
          }
        }
        if (segBoundDiff < 0 && segBoundDiff > -params->injSearchWindow)
        {
          if ((j+1) != params->numOverlapSegments)
          {
            segListToDo[segNumber+1] = 1;
          }
        }
      }
      injectList = injectList->next;
    }
  }
  else
  {
    params->injectList = NULL;
  }

 /* FIXME: For all sky mode trig start/end time needs to be implemented */
 /* these probably work in the ring code! */

  /* if todo list is empty then do them all */
  if ( (! params->segmentsToDoList || ! strlen( params->segmentsToDoList )) && (! params->analyzeInjSegsOnly) )
  {
    segments->numSgmnt = params->numOverlapSegments;
    segments->sgmnt = LALCalloc( segments->numSgmnt, sizeof(*segments->sgmnt) );
    for ( sgmnt = 0; sgmnt < params->numOverlapSegments; ++sgmnt )
    {
      slidSegNum = ( sgmnt + ( params->slideSegments[NumberIFO] ) ) % ( segments->numSgmnt );
      timeSlideVectors[NumberIFO*params->numOverlapSegments + sgmnt] = 
          (sgmnt-slidSegNum)*params->segmentDuration;
      compute_data_segment( &segments->sgmnt[sgmnt], slidSegNum, channel,
          invspec, response, params->segmentDuration, params->strideDuration,
          fwdplan );
    }
  }
  else  /* only do the segments in the todo list */
  {
    UINT4 count;

    /* first count the number of segments to do */
    count = 0;
    for ( sgmnt = 0; sgmnt < params->numOverlapSegments; ++sgmnt )
      if ( params->analyzeInjSegsOnly )
      {
        if ( segListToDo[sgmnt] == 1)
          ++count;
        else
          continue;
      }
      else
      {
        if ( is_in_list( sgmnt, params->segmentsToDoList ) )
          ++count;
        else
          continue; /* skip this segment: it is not in todo list */
      }

    if ( ! count ) /* no segments to do */
      return NULL;

    segments->numSgmnt = count;
    segments->sgmnt = LALCalloc( segments->numSgmnt, sizeof(*segments->sgmnt) );

    count = 0;
    for ( sgmnt = 0; sgmnt < params->numOverlapSegments; ++sgmnt )
    {
      if ( params->analyzeInjSegsOnly )
      {
        if ( segListToDo[sgmnt] == 1)
          compute_data_segment( &segments->sgmnt[count++], sgmnt, channel,
            invspec, response, params->segmentDuration, params->strideDuration,
            fwdplan );
      }
      else
      {
        if ( is_in_list( sgmnt, params->segmentsToDoList ) )
      /* we are sliding the names of segments here */
        {
          slidSegNum = ( sgmnt + ( params->slideSegments[NumberIFO] ) ) % ( segments->numSgmnt );
          timeSlideVectors[NumberIFO*params->numOverlapSegments + sgmnt] =
              ((INT4)slidSegNum-(INT4)sgmnt)*params->segmentDuration/2.;
          compute_data_segment( &segments->sgmnt[count++], slidSegNum, channel,
            invspec, response, params->segmentDuration, params->strideDuration,
            fwdplan );
        }
      }
    }
  }
  if ( params->writeSegment) /* write data segment */
  {
    for ( sgmnt = 0; sgmnt < segments->numSgmnt; ++sgmnt )
    {
      write_COMPLEX8FrequencySeries( &segments->sgmnt[sgmnt] ) ;
    }
  }

  return segments;
}


void coh_PTF_calculate_bmatrix(
  struct coh_PTF_params   *params,
  gsl_matrix *eigenvecs,
  gsl_vector *eigenvals,
  REAL4 Fplus[LAL_NUM_IFO],
  REAL4 Fcross[LAL_NUM_IFO],
  REAL8Array              *PTFM[LAL_NUM_IFO+1],
  UINT4 vecLength,
  UINT4 vecLengthTwo,
  UINT4 PTFMlen
  )
{
  // This function calculates the eigenvectors and eigenvalues of the
  // coherent "B" matrix. This is the matrix that appears as B^{-1} in the
  // original definition of SNR for spin.
  // For non spin, the eigenvalues can be used to rotate to the dominant
  // polarization frame.
  UINT4 i,j,k;
  UINT4 vecLengthSquare = vecLength*vecLength;
  REAL4 zh[vecLengthSquare],sh[vecLengthSquare],yu[vecLengthSquare];
  gsl_matrix *B2 = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_eigen_symmv_workspace *matTemp = gsl_eigen_symmv_alloc (vecLengthTwo);
  /* Create and invert the Bmatrix */
  /* Note for nonSpin PTFM contains one entry per ifo, this is loop is then
     a lot simpler than it looks! For PTF there are 25 entries in PTFM! 
     For non spin this stores (h_0|h_0)*/
  for (i = 0; i < vecLength; i++ )
  {
    for (j = 0; j < vecLength; j++ )
    {
      zh[i*vecLength+j] = 0;
      sh[i*vecLength+j] = 0;
      yu[i*vecLength+j] = 0;
      for( k = 0; k < LAL_NUM_IFO; k++)
      {
        if ( params->haveTrig[k] )
        {
          if ( params->faceOnStatistic )
          {
            zh[i*vecLength+j] += (Fplus[k]*Fplus[k] + Fcross[k] * Fcross[k])* PTFM[k]->data[i*PTFMlen+j];
          }
          else
          {
            zh[i*vecLength+j] += Fplus[k]*Fplus[k] * PTFM[k]->data[i*PTFMlen+j];
            sh[i*vecLength+j] += Fcross[k]*Fcross[k] * PTFM[k]->data[i*PTFMlen+j];
            yu[i*vecLength+j] += Fplus[k]*Fcross[k] * PTFM[k]->data[i*PTFMlen+j];
          }
        }
      }
    }
  }

  /* GSL is used to do the rotation */
  for (i = 0; i < vecLengthTwo; i++ )
  {
    for (j = 0; j < vecLengthTwo; j++ )
    {
      if ( i < vecLength && j < vecLength )
      {
        gsl_matrix_set(B2,i,j,zh[i*vecLength+j]);
      }
      else if ( i > (vecLength-1) && j > (vecLength-1))
      {
        gsl_matrix_set(B2,i,j,sh[(i-vecLength)*vecLength + (j-vecLength)]);
      }
      else if ( i < vecLength && j > (vecLength-1))
      {
        gsl_matrix_set(B2,i,j, yu[i*vecLength + (j-vecLength)]);
      }
      else if ( i > (vecLength-1) && j < vecLength)
      {
        gsl_matrix_set(B2,i,j,yu[j*vecLength + (i-vecLength)]);
      }
      else
        fprintf(stderr,"BUGGER! Something went wrong.");
    }
  }

  /* Here we compute the eigenvalues and eigenvectors of B2 */
  gsl_eigen_symmv (B2,eigenvals,eigenvecs,matTemp);
  gsl_eigen_symmv_free(matTemp);
  gsl_matrix_free(B2);
}


void coh_PTF_calculate_rotated_vectors(
    struct coh_PTF_params   *params,
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
    REAL4 *u1,
    REAL4 *u2,
    REAL4 Fplus[LAL_NUM_IFO],
    REAL4 Fcross[LAL_NUM_IFO],
    INT4  timeOffsetPoints[LAL_NUM_IFO],
    gsl_matrix *eigenvecs,
    gsl_vector *eigenvals,
    UINT4 numPoints,
    UINT4 position,
    UINT4 vecLength,
    UINT4 vecLengthTwo,
    UINT4 detectorNum)
{
  // This function calculates the coherent time series and rotates them into
  // the basis where the B matrix is the identity.
  // This is the dominant polarization frame with some normalization
  
  UINT4 j,k;  
  REAL4 v1[vecLengthTwo],v2[vecLengthTwo];

  for ( j = 0; j < vecLengthTwo ; j++ ) /* Construct the vi vectors */
  {
    v1[j] = 0.;
    v2[j] = 0.;
    if (detectorNum == LAL_NUM_IFO)
    {
      for( k = 0; k < LAL_NUM_IFO; k++)
      {
        if ( params->haveTrig[k] )
        {
          if ( params->faceOnStatistic)
          {
            /* Currently non-spin only! */
            v1[j] += Fplus[k] * PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]].re;
            v1[j] += Fcross[k] * PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]].im;
            if (params->faceOnStatistic == 1)
            {
              v2[j] += Fcross[k] * PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]].re;
              v2[j] -= Fplus[k] * PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]].im;
            }
            else if (params->faceOnStatistic == 2)
            {
              v2[j] -= Fcross[k] * PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]].re;
              v2[j] += Fplus[k] * PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]].im;
            }
            else
            {
              fprintf(stderr,"Face-on stat is not working!");
            }
          }
          else if (j < vecLength)
          {
            v1[j] += Fplus[k] * PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]].re;
            v2[j] += Fplus[k] * PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]].im;
          }
          else
          {
            v1[j] += Fcross[k] * PTFqVec[k]->data[(j-vecLength)*numPoints+position+timeOffsetPoints[k]].re;
            v2[j] += Fcross[k] * PTFqVec[k]->data[(j-vecLength)*numPoints+position+timeOffsetPoints[k]].im;
          }
        }
      }
    }
    else
    {
      v1[j] += PTFqVec[detectorNum]->data[j*numPoints+position].re;
      v2[j] += PTFqVec[detectorNum]->data[j*numPoints+position].im;
    }
  }

  /* Now we rotate the v1 and v2 to be in orthogonal basis */
  /* We can use gsl multiplication stuff to do this */
  /* BLAS stuff is complicated so we'll do it explicitly! */

  for ( j = 0 ; j < vecLengthTwo ; j++ )
  {
    u1[j] = 0.;
    u2[j] = 0.;
    for ( k = 0 ; k < vecLengthTwo ; k++ )
    {
      u1[j] += gsl_matrix_get(eigenvecs,k,j)*v1[k];
      u2[j] += gsl_matrix_get(eigenvecs,k,j)*v2[k];
    }
    u1[j] = u1[j] / (pow(gsl_vector_get(eigenvals,j),0.5));
    u2[j] = u2[j] / (pow(gsl_vector_get(eigenvals,j),0.5));
  }

}

void coh_PTF_cleanup(
    struct coh_PTF_params   *params,
    ProcessParamsTable      *procpar,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *revplan,
    COMPLEX8FFTPlan         *invPlan,
    REAL4TimeSeries         *channel[LAL_NUM_IFO+1],
    REAL4FrequencySeries    *invspec[LAL_NUM_IFO+1],
    RingDataSegments        *segments[LAL_NUM_IFO+1],
    MultiInspiralTable      *events,
    InspiralTemplate        *PTFbankhead,
    FindChirpTemplate       *fcTmplt,
    FindChirpTmpltParams    *fcTmpltParams,
    FindChirpInitParams     *fcInitParams,
    REAL8Array              *PTFM[LAL_NUM_IFO+1],
    REAL8Array              *PTFN[LAL_NUM_IFO+1],
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
    REAL4                   *timeOffsets,
    REAL4                   *Fplus,
    REAL4                   *Fcross,
    REAL4                   *Fplustrig,
    REAL4                   *Fcrosstrig
    )
{
  if ( params->injectList )
  {
    SimInspiralTable        *injectList = params->injectList;
    SimInspiralTable        *thisInject = NULL;
    while (injectList)
    {
      thisInject = injectList;
      injectList = injectList->next;
      LALFree(thisInject);
    }
  }

  /* Clean up memory usage */
  UINT4 ifoNumber;
  while ( events )
  {
    MultiInspiralTable *thisEvent;
    thisEvent = events;
    events = events->next;
    if ( thisEvent->event_id )
    {
      LALFree( thisEvent->event_id );
    }
    if ( thisEvent->time_slide_id )
    {
      LALFree( thisEvent->time_slide_id );
    }
    LALFree( thisEvent );
  }
  while ( PTFbankhead )
  {
    InspiralTemplate *thisTmplt;
    thisTmplt = PTFbankhead;
    PTFbankhead = PTFbankhead->next;
    if ( thisTmplt->event_id )
    {
      LALFree( thisTmplt->event_id );
    }
    LALFree( thisTmplt );
  }
  UINT4 sgmnt;
  for( ifoNumber = 0; ifoNumber < (LAL_NUM_IFO+1); ifoNumber++)
  {
    if ( segments[ifoNumber] )
    {
      for ( sgmnt = 0; sgmnt < segments[ifoNumber]->numSgmnt; sgmnt++ )
      {
        if (segments[ifoNumber]->sgmnt[sgmnt].data)
          XLALDestroyCOMPLEX8Vector(segments[ifoNumber]->sgmnt[sgmnt].data);
      }
      LALFree( segments[ifoNumber]->sgmnt );
      LALFree( segments[ifoNumber] );
    }
    if ( invspec[ifoNumber] )
    {
      XLALDestroyREAL4Vector( invspec[ifoNumber]->data );
      LALFree( invspec[ifoNumber] );
    }
    if ( channel[ifoNumber] )
    {
      XLALDestroyREAL4Vector( channel[ifoNumber]->data );
      LALFree( channel[ifoNumber] );
    }
    if ( PTFM[ifoNumber] )
      XLALDestroyREAL8Array( PTFM[ifoNumber] );
    if ( PTFN[ifoNumber] )
      XLALDestroyREAL8Array( PTFN[ifoNumber] );
    if ( PTFqVec[ifoNumber] )
      XLALDestroyCOMPLEX8VectorSequence( PTFqVec[ifoNumber] );
  }
  if ( revplan )
    XLALDestroyREAL4FFTPlan( revplan );
  if ( fwdplan )
    XLALDestroyREAL4FFTPlan( fwdplan );
  if ( invPlan )
    XLALDestroyCOMPLEX8FFTPlan( invPlan );
  while ( procpar )
  {
    ProcessParamsTable *thisParam;
    thisParam = procpar;
    procpar = procpar->next;
    LALFree( thisParam );
  }
  if (fcTmpltParams)
  {
    if ( fcTmpltParams->fwdPlan )
      XLALDestroyREAL4FFTPlan( fcTmpltParams->fwdPlan );
    if ( fcTmpltParams->PTFe1 )
      XLALDestroyVectorSequence( fcTmpltParams->PTFe1 );
    if ( fcTmpltParams->PTFe2 )
      XLALDestroyVectorSequence( fcTmpltParams->PTFe2 );
    if ( fcTmpltParams->PTFphi )
      XLALDestroyVector( fcTmpltParams->PTFphi );
    if ( fcTmpltParams->PTFomega_2_3 )
      XLALDestroyVector( fcTmpltParams->PTFomega_2_3 );
    if ( fcTmpltParams->xfacVec )
      XLALDestroyVector( fcTmpltParams->xfacVec );
    LALFree( fcTmpltParams );
  }
  if ( fcTmplt )
  {
    if ( fcTmplt->PTFQ )
      XLALDestroyVectorSequence( fcTmplt->PTFQ );
    if ( fcTmplt->PTFQtilde )
      XLALDestroyCOMPLEX8VectorSequence( fcTmplt->PTFQtilde );
    if ( fcTmplt->data )
      XLALDestroyCOMPLEX8Vector( fcTmplt->data );
    LALFree( fcTmplt );
  }
  if ( fcInitParams )
    LALFree( fcInitParams );
  if ( timeOffsets )
    LALFree( timeOffsets );
  if ( Fplus )
    LALFree( Fplus );
  if ( Fcross )
    LALFree( Fcross );
  if ( Fplustrig )
    LALFree( Fplustrig );
  if ( Fcrosstrig )
    LALFree( Fcrosstrig );
}


/* gets the forward fft plan */
REAL4FFTPlan *coh_PTF_get_fft_fwdplan( struct coh_PTF_params *params )
{
  REAL4FFTPlan *plan = NULL;
  if ( params->segmentDuration > 0.0 )
  {
    UINT4 segmentLength;
    segmentLength = floor( params->segmentDuration * params->sampleRate + 0.5 );
    plan = XLALCreateForwardREAL4FFTPlan( segmentLength, params->fftLevel );
  }
  return plan;
}


/* gets the reverse fft plan */
REAL4FFTPlan *coh_PTF_get_fft_revplan( struct coh_PTF_params *params )
{
  REAL4FFTPlan *plan = NULL;
  if ( params->segmentDuration > 0.0 )
  {
    UINT4 segmentLength;
    segmentLength = floor( params->segmentDuration * params->sampleRate + 0.5 );
    plan = XLALCreateReverseREAL4FFTPlan( segmentLength, params->fftLevel );
  }
  return plan;
}

/* routine to see if integer i is in a list of integers to do */
/* e.g., 2, 7, and 222 are in the list "1-3,5,7-" but 4 is not */
int is_in_list( int i, const char *list )
{
  char  buffer[BUFFER_SIZE];
  char *str = buffer;
  int   ans = 0;

  strncpy( buffer, list, sizeof( buffer ) - 1 );

  while ( str )
  {
    char *tok;  /* token in a delimited list */
    char *tok2; /* second part of token if it is a range */
    tok = str;

    if ( ( str = strchr( str, ',' ) ) ) /* look for next delimiter */
      *str++ = 0; /* nul terminate current token; str is remaining string */

    /* now see if this token is a range */
    if ( ( tok2 = strchr( tok, '-' ) ) )
      *tok2++ = 0; /* nul terminate first part of token; tok2 is second part */  
    if ( tok2 ) /* range */
    {
      int n1, n2;
      if ( strcmp( tok, "^" ) == 0 )
        n1 = INT_MIN;
      else
        n1 = atoi( tok );
      if ( strcmp( tok2, "$" ) == 0 )
        n2 = INT_MAX;
      else
        n2 = atoi( tok2 );
      if ( i >= n1 && i <= n2 ) /* see if i is in the range */
        ans = 1;
    }
    else if ( i == atoi( tok ) )
      ans = 1;

    if ( ans ) /* i is in the list */
      break;
  }

  return ans;
}

SnglInspiralTable *conv_insp_tmpl_to_sngl_table(
    InspiralTemplate        *template,
    UINT4                   eventNumber
    )
{
  SnglInspiralTable *cnvTemplate;
  cnvTemplate = (SnglInspiralTable *) LALCalloc(1,sizeof(SnglInspiralTable));
  cnvTemplate->event_id = (EventIDColumn *)
      LALCalloc(1, sizeof(EventIDColumn) );
  cnvTemplate->event_id->id=eventNumber;
  cnvTemplate->mass1 = template->mass1;
  cnvTemplate->mass2 = template->mass2;
  cnvTemplate->chi = template->chi;
  cnvTemplate->kappa = template->kappa;
  cnvTemplate->f_final = template->fCutoff;
  return cnvTemplate;
}

/*
 *
 * construct CohPTFSkyPositions structures for the different sky patching cases:
 *
 * if all sky:
 *   FAIL - not implemented yet
 * if two-detectors with patch:
 *   generate an arc of points perpendicular to the line of constant time-delay
 * if multiple-detector patch or single point:
 *   generate a grid of concentric circles around trigger point
 *
 */

CohPTFSkyPositions *coh_PTF_generate_sky_points( 
    struct coh_PTF_params *params
    )
{

  CohPTFSkyPositions *skyPoints = NULL;

  /* if given file */
  if ( params->skyPositionsFile != NULL )
  {
    UINT4 raColumn  = 0;
    UINT4 decColumn = 1;
    skyPoints = coh_PTF_read_grid_from_file(params->skyPositionsFile,\
                                            raColumn, decColumn);
  }

  /* if all sky */
  else if (params->skyLooping == TWO_DET_ALL_SKY)
  {
    verbose("Generating 2 detector all sky map...\n");
    skyPoints = coh_PTF_two_det_sky_grid(params);
  }

  else if ((params->skyLooping == ALL_SKY) && (params->numIFO==3))
  {
    verbose("Generating 3 detector all sky map...\n");
    skyPoints = coh_PTF_three_det_sky_grid(params);
    //error("all sky mode is not implemented yet, however, you can use --sky-positions-file\n");

    //skyPoints = coh_PTF_sky_grid()
  }

  else if (params->skyLooping == ALL_SKY)
  {
    error("all sky mode is not implemented yet, however, you can use --sky-positions-file\n");
  }

  /* if sky region */
  else if ( params->skyLooping == SKY_PATCH\
            || params->skyLooping == TWO_DET_SKY_PATCH\
            || params->skyLooping == SINGLE_SKY_POINT )
  {
    skyPoints = coh_PTF_generate_sky_grid(params);
  }

  /* if two-detectors, remove time-delay degeneracy */
  if (params->skyLooping == TWO_DET_SKY_PATCH)
  {
    verbose("Generated necessary sky grid with %d points, ",
            skyPoints->numPoints);
    verbose("parsing for time-delay degeneracy\n");
    CohPTFSkyPositions *parsedSkyPoints = NULL; 
    parsedSkyPoints = coh_PTF_parse_time_delays(skyPoints, params);
    if (skyPoints->data)
      LALFree(skyPoints->data);
    if (skyPoints)
      LALFree(skyPoints);
    return parsedSkyPoints;
  }
  else
  {
    return skyPoints;
  }
}

/*
 * Generate a grid of sky points based on the sinusoidal map method used by
 * the xpipeline.
 */

CohPTFSkyPositions *coh_PTF_generate_sky_grid(
    struct coh_PTF_params *params    
    )
{
  CohPTFSkyPositions *skyPoints = NULL;
  UINT4              ifoNumber,i,j;
  LALDetector        *detectors[LAL_NUM_IFO];
  REAL4              angle;  /* angle between IFO baseline & sky localisation */
  REAL4              lambdamin,lambdamax,lambda; /* angle closest to pi/2 */
  REAL4              alpha=0,detalpha;
  REAL4              angularResolution;
  double             baseline,lightTravelTime;
  REAL4              raNp  = 0.;                 /* north */
  REAL4              decNp = LAL_PI_2;           /* pole */

  gsl_vector         *axis;                       /* rotation axis vector */
  gsl_vector         *npPos;                      /* north pole position */
  gsl_vector         *trigPos;                    /* trigger position */

  /* get site coordinates */
  for(ifoNumber=0; ifoNumber<LAL_NUM_IFO; ifoNumber++)
  {
    detectors[ifoNumber] = LALCalloc(1, sizeof(*detectors[ifoNumber]));
    XLALReturnDetector(detectors[ifoNumber], ifoNumber);
  }

  /* find pair of detectors whose opening angle to the GRB ±error is closest */
  /* to 90 degrees */
  for (i=0; i<LAL_NUM_IFO; i++)
  {
    if (params->haveTrig[i])
    {
      for (j=i+1; j<LAL_NUM_IFO; j++)
      {
        if (params->haveTrig[j])
        {
          /* get dot product (time delay) between sites */
          baseline = XLALArrivalTimeDiff(detectors[i]->location,
                                         detectors[j]->location,
                                         params->rightAscension,
                                         params->declination,
                                         &params->trigTime);

          /* get light travel time */
          lightTravelTime = XLALLightTravelTime(detectors[i], detectors[j]);
          lightTravelTime *= 1e-9;

          /* calculate opening angle */
          angle  = acos(baseline/lightTravelTime);

          /* generate angular window with sky error */
          lambdamin = angle-params->skyError;
          lambdamax = angle+params->skyError;
   
          /* if pi/2 is in the range, choose that, 
           * otherwise get as close as possible */
          if (lambdamin < LAL_PI_2 && lambdamax > LAL_PI_2)
          {
            lambda = LAL_PI_2;
          }
          else
            if (fabs(LAL_PI_2-lambdamin) < fabs(LAL_PI_2-lambdamax))
              lambda = lambdamin;
            else
              lambda = lambdamax;

          /* calculate alpha */
          detalpha = lightTravelTime * sin(lambda);
          if (detalpha > alpha)
            alpha = detalpha;
    
        }
      }
    }
  }

  /* calculate angular resolution */
  if (! params->singlePolFlag && (! params->numIFO == 1))
  {
    angularResolution = 2. * params->timingAccuracy / alpha;
  }
  else
  {
    angularResolution = 1;
    params->skyError = 0;
  }

  /* generate sky grid using sinusoidal map */
  skyPoints = coh_PTF_circular_grid(angularResolution, params->skyError);

  /*
   * Rotate sky grid to centre on the given (ra,dec)
   */
  
  /* calculate angle between north pole and (ra,dec) */
  raNp  = 0.;
  decNp = LAL_PI_2;

  angle = acos ( sin(decNp)*sin(params->declination) +
                 cos(decNp)*cos(params->declination) * 
                 cos( raNp-params->rightAscension ) );

  /* calculate unit vector to rotate around */
  npPos = gsl_vector_alloc(3);
  trigPos = gsl_vector_alloc(3);
  axis = gsl_vector_alloc(3);

  gsl_vector_set(npPos, 0, 0.);
  gsl_vector_set(npPos, 1, 0.);
  gsl_vector_set(npPos, 2, 1.);

  gsl_vector_set(trigPos, 0,\
                 cos(params->declination)*cos(params->rightAscension));
  gsl_vector_set(trigPos, 1,\
                 cos(params->declination)*sin(params->rightAscension));
  gsl_vector_set(trigPos, 2, sin(params->declination));

  cross_product(axis, npPos, trigPos);
  normalise(axis);

  /* rotate sky points */
  coh_PTF_rotate_skyPoints(skyPoints, axis, angle);
  
  /* free memory */
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if ( detectors[ifoNumber] )
      LALFree(detectors[ifoNumber]);
  }
  FREE_GSL_VECTOR(npPos);
  FREE_GSL_VECTOR(trigPos);
  FREE_GSL_VECTOR(axis);

  return skyPoints;

}

/*
 * generate circular map of sky points centred on north pole:
 *  * place central point
 *  * step out in declination by angularResolution
 *  * place a ring of points separated by angularResolution
 *  * repeat until maximum skyError radius passed
 */

CohPTFSkyPositions *coh_PTF_circular_grid(
    REAL4        angularResolution,
    REAL4        skyError
    )
{

  UINT4                    p               = 0;
  UINT4                    i,j;
  UINT4                    numTheta;            /* number of rings */
  REAL4                    dPhi,phi,theta;      /* sky point parameters */
  UINT4                    numSkyPoints    = 0; /* total number of sky points */
  CohPTFSkyPositions       *skyPoints;

  /* set range of theta */
  numTheta = (int) ceil( skyError / angularResolution ) + 1;

  /* set number of sky points */
  UINT4 numPhi[numTheta];
  for ( i=0; i < numTheta; i++ )
  {
    theta     = angularResolution * i;
    numPhi[i] = (int) ceil( LAL_TWOPI * sin(theta) / angularResolution );
    if ( numPhi[i] < 1 )
      numPhi[i] = 1;
    numSkyPoints += numPhi[i];
  }

  /* assign memory for sky points */
  skyPoints = LALCalloc(1, sizeof(*skyPoints));
  skyPoints->numPoints = numSkyPoints;
  skyPoints->data      = LALCalloc(1, numSkyPoints*sizeof(SkyPosition));

  /* loop over rings on sky, and around each ring */
  for ( i=0; i < numTheta; i++ )
  {
    dPhi  = LAL_TWOPI / numPhi[i];
    theta = angularResolution * i;
    for ( j=0; j < numPhi[i]; j++ )
    {
      /* calculate phi */
      phi = ( -LAL_PI + dPhi / 2. ) + dPhi * j;
      /* assign sky point */
      skyPoints->data[p].longitude  = phi;
      skyPoints->data[p].latitude = LAL_PI_2 - theta;
      skyPoints->data[p].system = COORDINATESYSTEM_EQUATORIAL;
      XLALNormalizeSkyPosition(&skyPoints->data[p].longitude, &skyPoints->data[p].latitude);

      p++;
    }

  }

  return skyPoints;

}

CohPTFSkyPositions *coh_PTF_parse_time_delays(
    CohPTFSkyPositions    *skyPoints,
    struct coh_PTF_params *params
)
{

  CohPTFSkyPositions *parsedSkyPoints;
  UINT4              numSkyPoints = 0;
  UINT4              i, p, ifoNumber, appendPoint[skyPoints->numPoints];
  REAL8              timeDelay, timeDelays[skyPoints->numPoints];
  LALDetector        *detectors[params->numIFO];
  REAL8              dt = params->timingAccuracy;

  /* get site coordinates */
  i=0; 
  for(ifoNumber=0; ifoNumber<LAL_NUM_IFO; ifoNumber++)
  {
    if (params->haveTrig[ifoNumber])
    {
      detectors[i] = LALCalloc(1, sizeof(*detectors[i]));
      XLALReturnDetector(detectors[i], ifoNumber);
      i++;
    }
  }

  for (p = 0; p < skyPoints->numPoints; p++)
  {

    /* set default timeDelay to zero */
    timeDelays[p] = 0;
    /* set default stance to keep point */
    appendPoint[p] = 0;
    timeDelay = XLALArrivalTimeDiff(detectors[0]->location,
                                    detectors[1]->location,
                                    skyPoints->data[p].longitude,
                                    skyPoints->data[p].latitude,
                                    &params->trigTime);

    /* loop over timeDelays list */
    for (i = 0; i < p; i++)
    {
      /* if we haven't assigned this timeDelay, move on */
      if (timeDelays[i]==0)
      {
        continue;
      }

      /* if we already have a timeDelay within the timing accuracy,
       * we don't want another one */
      if (fabsf(timeDelay-timeDelays[i]) < dt)
      {
        appendPoint[p] = 1;
        break;
      }
    } 

    /* if we want to keep this point, save the time delay and increment the
     * counter */
    if (appendPoint[p] == 0)
    {
      numSkyPoints++;
      timeDelays[p] = timeDelay;
    }
  }

  /* assign memory for sky points */
  parsedSkyPoints = LALCalloc(1, sizeof(*parsedSkyPoints));
  parsedSkyPoints->numPoints = numSkyPoints;
  parsedSkyPoints->data =
      LALCalloc(1, numSkyPoints*sizeof(*parsedSkyPoints->data));

  /* save the new list of points */
  i = 0;
  for (p = 0; p < skyPoints->numPoints; p++)
  {
    if (appendPoint[p] == 0)
    {
      parsedSkyPoints->data[i] = skyPoints->data[p];
      i++;
    }
  }

  /* free memory */
  for(ifoNumber = 0; ifoNumber < params->numIFO; ifoNumber++)
  {
    if (detectors[ifoNumber])
      LALFree(detectors[ifoNumber]);
  }

  return parsedSkyPoints;
}

long int timeval_subtract(struct timeval *t1)
{
  struct timeval t2;
  gettimeofday(&t2,NULL);
  long int diff = (t2.tv_usec + 1000000 * t2.tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);

  return diff;
}

void timeval_print(struct timeval *tv)
{
  char buffer[30];
  time_t curtime;

  printf("%lu.%06lu", (unsigned long)tv->tv_sec, (unsigned long)tv->tv_usec);
  curtime = tv->tv_sec;
  strftime(buffer, 30, "%m-%d-%Y  %T", localtime(&curtime));
  printf(" = %s.%06lu\n", buffer, (unsigned long)tv->tv_usec);
}

/*
 *
 * rotate a vector about a given angle in three dimensions
 *
 */

void coh_PTF_rotate_skyPoints(
    CohPTFSkyPositions *skyPoints,
    gsl_vector *axis,
    REAL8 angle
)
{
  /* initialise variables */
  UINT4 i;
  gsl_matrix *matrix;
  matrix = gsl_matrix_alloc(3, 3);

  /* construct rotation matrix */
  rotation_matrix(matrix, axis, angle);

  /* loop over points rotating by angle around axis */
  for ( i=0; i < skyPoints->numPoints; i++ )
  {
    coh_PTF_rotate_SkyPosition(&skyPoints->data[i], matrix);
  }
}

void coh_PTF_rotate_SkyPosition(
    SkyPosition *skyPoint,
    gsl_matrix  *matrix
)
{
  /* initialise variables */
  REAL4 phi,theta;
  gsl_vector *pos;    /* original position vector */
  gsl_vector *rotPos; /* rotated position vector */

  
  phi   = skyPoint->longitude;
  theta = LAL_PI_2 - skyPoint->latitude;
  /* convert to cartesian */
  pos = gsl_vector_alloc(3);
  gsl_vector_set(pos, 0, sin(theta)*cos(phi));
  gsl_vector_set(pos, 1, sin(theta)*sin(phi));
  gsl_vector_set(pos, 2, cos(theta));

  /* rotate */
  rotPos = gsl_vector_alloc(3);
  gsl_blas_dgemv(CblasNoTrans, 1.0, matrix, pos, 0.0, rotPos);

  /* convert back to (phi,theta) */
  theta = acos(gsl_vector_get(rotPos, 2));
  phi   = atan2(gsl_vector_get(rotPos, 1), gsl_vector_get(rotPos, 0));
  //verbose("theta2 = %e, phi2 = %e\n", theta, phi);
  skyPoint->longitude = phi;
  skyPoint->latitude  = LAL_PI_2 - theta;
  XLALNormalizeSkyPosition(&skyPoint->longitude, &skyPoint->latitude);

  /* free memory */
  FREE_GSL_VECTOR(pos);
  FREE_GSL_VECTOR(rotPos);
}

void cross_product(
    gsl_vector *product,
    const gsl_vector *u,
    const gsl_vector *v
)
{
  double p1 = gsl_vector_get(u, 1)*gsl_vector_get(v, 2)
              - gsl_vector_get(u, 2)*gsl_vector_get(v, 1);

  double p2 = gsl_vector_get(u, 2)*gsl_vector_get(v, 0)
              - gsl_vector_get(u, 0)*gsl_vector_get(v, 2);

  double p3 = gsl_vector_get(u, 0)*gsl_vector_get(v, 1)
              - gsl_vector_get(u, 1)*gsl_vector_get(v, 0);

  gsl_vector_set(product, 0, p1);
  gsl_vector_set(product, 1, p2);
  gsl_vector_set(product, 2, p3);
}

void normalise(
    gsl_vector *vec
)
{
  double mag;
  /* calculate magnitude of vector */
  gsl_blas_ddot(vec, vec, &mag);
  mag = sqrt(mag);
  /* scale vector by inverse magnitude */
  gsl_vector_scale(vec, 1./mag);
}

void rotation_matrix(
    gsl_matrix *matrix,
    gsl_vector *axis,
    REAL8 angle
)
{

  gsl_matrix_set(matrix, 0, 0, cos(angle) +\
                               pow(gsl_vector_get(axis, 0),2)*(1-cos(angle)));
  gsl_matrix_set(matrix, 0, 1, gsl_vector_get(axis, 0)*gsl_vector_get(axis, 1)*\
                                   (1-cos(angle)) -\
                               gsl_vector_get(axis, 2)*sin(angle));
  gsl_matrix_set(matrix, 0, 2, gsl_vector_get(axis, 0)*gsl_vector_get(axis, 2)*\
                                   (1-cos(angle)) +\
                               gsl_vector_get(axis, 1)*sin(angle));

  gsl_matrix_set(matrix, 1, 0, gsl_vector_get(axis, 1)*gsl_vector_get(axis, 0)*\
                                   (1-cos(angle)) +\
                               gsl_vector_get(axis, 2)*sin(angle));
  gsl_matrix_set(matrix, 1, 1, cos(angle) +\
                               pow(gsl_vector_get(axis, 1),2)*(1-cos(angle)));
  gsl_matrix_set(matrix, 1, 2, gsl_vector_get(axis, 1)*gsl_vector_get(axis, 2)*\
                                  (1-cos(angle)) -\
                               gsl_vector_get(axis, 0)*sin(angle));

  gsl_matrix_set(matrix, 2, 0, gsl_vector_get(axis, 2)*gsl_vector_get(axis, 0)*\
                                   (1-cos(angle)) -\
                               gsl_vector_get(axis, 1)*sin(angle));
  gsl_matrix_set(matrix, 2, 1, gsl_vector_get(axis, 2)*gsl_vector_get(axis, 1)*\
                                   (1-cos(angle)) +\
                               gsl_vector_get(axis, 0)*sin(angle));
  gsl_matrix_set(matrix, 2, 2, cos(angle) +\
                               pow(gsl_vector_get(axis, 2),2)*(1-cos(angle)));

}

CohPTFSkyPositions *coh_PTF_read_grid_from_file(
    const char *fname,
    UINT4      raColumn,
    UINT4      decColumn
)
{

  UINT4              j, i=0, numSkyPoints=0; /* counters */
  CohPTFSkyPositions *skyPoints;             /* sky positions structure */
  char               *value, line[256]; /* string holders */
  FILE               *data;                  /* file object */

  /* read file */
  data = fopen(fname, "r");

  /* check file */
  if (data == NULL)
  {
    error("Error reading sky locations file %s. Please verify this path and try again\n", fname);
  }

  /* find number of lines */
  while (fgets(line, sizeof(line), data))
  {
    numSkyPoints += 1;
  }

  /* seek to start of file again */
  fseek(data, 0, SEEK_SET);  

  /* assign memory for sky points */
  skyPoints = LALCalloc(1, sizeof(*skyPoints));
  skyPoints->numPoints = numSkyPoints;
  skyPoints->data      = LALCalloc(1, numSkyPoints*sizeof(SkyPosition));

  /* find last column we need */
  UINT4 lastColumn = raColumn;
  if (decColumn > raColumn)
    lastColumn = decColumn;

  /* read data line by line */
  while (fgets(line, sizeof(line), data))
  {
    /* set counter */
    j = 0;

    /* extract first value */
    value = strtok(line, " ");

    /* loop over the columns and extract the correct ones */
    while (j <= lastColumn)
    {
      if (j==raColumn)
        skyPoints->data[i].longitude = (REAL4) atof(value) * LAL_PI_180;
      else if (j==decColumn)
        skyPoints->data[i].latitude  = (REAL4) atof(value) * LAL_PI_180;

      /* move on to next column */
      value = strtok(NULL, " ");
      j++;
    }

    i++;
  }

  return skyPoints;

}

CohPTFSkyPositions *coh_PTF_two_det_sky_grid(
    struct coh_PTF_params *params
)
{
  /* set up variables */
  CohPTFSkyPositions *skyPoints;
  CohPTFSkyPositions *geoSkyPoints;
  LALDetector        *detectors[params->numIFO];
  REAL8              lightTravelTime, timeDelay, angle;
  gsl_vector         *locations[params->numIFO], *normal, *northPole, *axis;
  UINT4              i, ifoNumber, numSkyPoints;

  /* get site coordinates */
  i=0;
  for(ifoNumber=0; ifoNumber<LAL_NUM_IFO; ifoNumber++)
  {
    if (params->haveTrig[ifoNumber])
    {
      detectors[i] = LALCalloc(1, sizeof(*detectors[i]));
      XLALReturnDetector(detectors[i], ifoNumber);
      locations[i] = gsl_vector_alloc(3);
      REALToGSLVector(detectors[i]->location, locations[i], 3);
      i++;
    }
  }

  /* get light travel time */
  lightTravelTime = XLALLightTravelTime(detectors[0], detectors[1]);
  lightTravelTime *= 1e-9;

  /* calculate number of skypoints */
  numSkyPoints = 2*(UINT4) floor(lightTravelTime/params->timingAccuracy) + 1;
   
  /* assign memory for sky points */
  geoSkyPoints = LALCalloc(1, sizeof(CohPTFSkyPositions));
  geoSkyPoints->numPoints = numSkyPoints;
  geoSkyPoints->data = LALCalloc(1,geoSkyPoints->numPoints*sizeof(SkyPosition));

  /* generate arc across equator */
  for (i=0; i < numSkyPoints; i++)
  {
    timeDelay = (INT4) (i - numSkyPoints/2) * params->timingAccuracy;
    //verbose("%f\n", timeDelay);
    geoSkyPoints->data[i].longitude = acos(-timeDelay/lightTravelTime);
    //verbose("%f\n", geoSkyPoints->data[i].longitude);
    geoSkyPoints->data[i].latitude  = 0.;
    geoSkyPoints->data[i].system    = COORDINATESYSTEM_GEOGRAPHIC;
  }

  /* calculate normal for time delay arc */
  normal = gsl_vector_alloc(3);
  cross_product(normal, locations[0], locations[1]);
  normalise(normal);

  /* calculate angle between this normal, and normal of equator */
  /* calculate unit vector to rotate around */
  northPole = gsl_vector_alloc(3);
  gsl_vector_set(northPole, 0, 0);
  gsl_vector_set(northPole, 1, 0);
  gsl_vector_set(northPole, 2, 1);
  gsl_blas_ddot(normal, northPole, &angle);
  angle = acos(angle);
 
  /* generate rotation vector */
  axis = gsl_vector_alloc(3);
  cross_product(axis, normal, northPole);

  /* rotate arc into plane of detectors and earth centre */
  coh_PTF_rotate_skyPoints(geoSkyPoints, axis, angle);
  //verbose("\n");

  /* convert from earth-fixed to sky-fixed */
  LALStatus status = blank_status;
  skyPoints = LALCalloc(1, sizeof(*skyPoints));
  skyPoints->numPoints = geoSkyPoints->numPoints;
  skyPoints->data      = LALCalloc(1, skyPoints->numPoints*sizeof(SkyPosition));
  for (i=0; i<numSkyPoints; i++)
  {
    //verbose("%f\n", XLALArrivalTimeDiff(detectors[0]->location, detectors[1]->location, geoSkyPoints->data[i].longitude, geoSkyPoints->data[i].latitude, &params->trigTime));
    LALGeographicToEquatorial(&status, &skyPoints->data[i],
                              &geoSkyPoints->data[i], &params->trigTime);
    XLALNormalizeSkyPosition(&skyPoints->data[i].longitude, &skyPoints->data[i].latitude);
    //verbose("%f\n", XLALArrivalTimeDiff(detectors[0]->location, detectors[1]->location, skyPoints->data[i].longitude, skyPoints->data[i].latitude, &params->trigTime));
  }

  /* free memory */
  for( ifoNumber = 0; ifoNumber < params->numIFO; ifoNumber++)
  {
    if (detectors[ifoNumber] )
      LALFree(detectors[ifoNumber]);
    if (locations[ifoNumber])
      FREE_GSL_VECTOR(locations[ifoNumber]);
  }
  FREE_GSL_VECTOR(normal);
  FREE_GSL_VECTOR(northPole);
  FREE_GSL_VECTOR(axis);

  return skyPoints;
}

/*
 * Set up three detector sky grid by tiling time-delay between detector 0 and
 * detector 1, and for each of those values tiling time-delay between detector
 * 0 and detector 2
 */

CohPTFSkyPositions *coh_PTF_three_det_sky_grid(
    struct coh_PTF_params *params
){

  /* set up variables */
  CohPTFSkyPositions *skyPoints;
  LALDetector        *detectors[params->numIFO];
  gsl_vector         *locations[params->numIFO], *baseline[params->numIFO],\
                     *northPole, *xaxis, *normal;
  gsl_matrix         *matrix;
  REAL8              T[params->numIFO], t2, t3, A, B, angle, xangle, zangle,\
                     condition, xphi, nphi, ntheta;
  UINT4              i, j, k, p, ifoNumber, numXPoints, numYPoints,\
                     numSkyPoints, totalPoints;

  /* get site coordinates */
  i=0;
  for(ifoNumber=0; ifoNumber<LAL_NUM_IFO; ifoNumber++)
  {
    if (params->haveTrig[ifoNumber])
    {
      detectors[i] = LALCalloc(1, sizeof(*detectors[i]));
      XLALReturnDetector(detectors[i], ifoNumber);
      locations[i] = gsl_vector_alloc(3);
      REALToGSLVector(detectors[i]->location, locations[i], 3);
      if (i>=1)
      {
        /* get light travel times relative to first detector */
        T[i-1] = XLALLightTravelTime(detectors[0], detectors[i]);
        T[i-1] *= 1e-9;
        /* get detector baselines */
        baseline[i-1] = gsl_vector_alloc(3);
        gsl_vector_memcpy(baseline[i-1], locations[i]);
        gsl_vector_sub(baseline[i-1], locations[0]);
        normalise(baseline[i-1]);
      } 
      i++;
    }
  }

  /* calculate angle between baselines */
  gsl_blas_ddot(baseline[0], baseline[1], &angle);
  angle = acos(angle);
  verbose("angle = %e\n", angle);

  /* calculate number of points spanning first two-detector time delay */
  numXPoints = 2*(UINT4) floor(T[0]/params->timingAccuracy) + 1;
  numYPoints = 2*(UINT4) floor(T[1]/params->timingAccuracy) + 1;
  numSkyPoints = 0;

  totalPoints = numXPoints * numYPoints;
  REAL4 valid[totalPoints];

  /* calculate number of points in y direction for each x */
  k=0;
  for (i=0; i<numXPoints; i++)
  {
    t2 = (INT4) (i - numXPoints/2) * params->timingAccuracy;
    for (j=0; j<numYPoints; j++)
    {
      t3 = (INT4) (j - numYPoints/2) * params->timingAccuracy;
      A = -T[1]/T[0] * t2 * cos(angle);
      B = pow(T[1], 2) * (pow(t2/T[0], 2) - pow(sin(angle), 2));
      condition = pow(t3, 2) + 2*A*t3 + B;
      if (condition <= 0)
      {
        valid[k] = 0;
        numSkyPoints++;
      }
      else
      {
        valid[k] = 1;
      }
      k++;
    }
  }

  /* assign memory for sky points */
  skyPoints = LALCalloc(1, sizeof(*skyPoints));
  skyPoints->numPoints = numSkyPoints;
  skyPoints->data      = LALCalloc(1, numSkyPoints*sizeof(SkyPosition));

  /*
   * Construct rotations from Rabaste network coordinates to geographical
   * coordinates
   */

  normal = gsl_vector_alloc(3);
  matrix = gsl_matrix_alloc(3, 3);
  xaxis  = gsl_vector_alloc(3);

  /* construct x-axis of network coordinates */
  xangle = LAL_PI_2 - angle;
  cross_product(normal, baseline[0], baseline[1]); 
  rotation_matrix(matrix, normal, xangle);
  /* apply rotation */
  gsl_blas_dgemv(CblasNoTrans, 1.0, matrix, baseline[1], 0.0, xaxis);

  /*
   * Rabaste network coordinate system has z-axis as the baseline between
   * detectors 1 and 2.
   * We need to rotate that onto the geographical north pole:
   */
  
  northPole = gsl_vector_alloc(3);

  /* construct rotation matrix */
  gsl_vector_set(northPole, 0, 0);
  gsl_vector_set(northPole, 1, 0);
  gsl_vector_set(northPole, 2, 1);
  gsl_blas_ddot(baseline[0], northPole, &zangle);
  zangle = acos(zangle);
  cross_product(normal, baseline[0], northPole);
  rotation_matrix(matrix, normal, zangle);

  verbose("xangle = %e\n", xangle);
  verbose("zangle = %e\n", zangle);

  /* rotate xaxis */
  gsl_vector *xaxis2 = gsl_vector_alloc(3);
  gsl_blas_dgemv(CblasNoTrans, 1.0, matrix, xaxis2, 0.0, xaxis);
  normalise(xaxis2);
  xphi = atan2(gsl_vector_get(xaxis2, 1), gsl_vector_get(xaxis2, 0));

  verbose("xphi = %e\n", xphi);

  /* assign sky points */
  /* calculate number of points in y direction for each x */
  k = 0;
  p = 0;
  for (i=0; i<numXPoints; i++)
  {
    t2 = (INT4) (i - numXPoints/2) * params->timingAccuracy;
    for (j=0; j<numYPoints; j++)
    {
      t3 = (INT4) (j - numYPoints/2) * params->timingAccuracy;
      if (valid[k]==0)
      {
        /* calculate (phi, theta) in network coordinates */
        ntheta = acos(-t2/T[0]);
        nphi   = acos(-(T[0]*t3-T[1]*t2*cos(angle))/\
                       (T[1]*sqrt(pow(T[0],2)-pow(t2,2))*sin(angle)));
        skyPoints->data[p].longitude = nphi;
        skyPoints->data[p].latitude  = ntheta-LAL_PI_2;
        skyPoints->data[p].system    = COORDINATESYSTEM_EQUATORIAL;
        XLALNormalizeSkyPosition(&skyPoints->data[p].longitude, &skyPoints->data[p].latitude);
        coh_PTF_rotate_SkyPosition(&skyPoints->data[i], matrix);
        skyPoints->data[p].longitude -= xphi;

        p++;
      }
      k++;
    }
  }

  return skyPoints;

}

void REALToGSLVector(
    const REAL8 *input,
    gsl_vector  *output,
    size_t      size
)
{
  UINT4 i;
  for (i=0; i<size; i++)
  {
    gsl_vector_set(output, i, input[i]);
  }
}

void findInjectionSegment(
    UINT4 *start,
    UINT4 *end,
    LIGOTimeGPS *epoch,
    struct coh_PTF_params *params
    )
{
    /* define variables */
    LIGOTimeGPS injTime, segmentStart, segmentEnd;
    UINT4 injSamplePoint, injWindow, numPoints;
    REAL8 injDiff;
    INT8 startDiff, endDiff;
    SimInspiralTable *thisInject = NULL;

    /* set variables */
    segmentStart = *epoch;
    segmentEnd   = *epoch;
    XLALGPSAdd(&segmentEnd, params->segmentDuration/2.0);
    thisInject = params->injectList;
    numPoints = floor(params->segmentDuration * params->sampleRate + 0.5);

    /* loop over injections */
    while (thisInject)
    {
        injTime = thisInject->geocent_end_time;
        startDiff = XLALGPSToINT8NS(&injTime) - XLALGPSToINT8NS(&segmentStart);
        endDiff = XLALGPSToINT8NS(&injTime) - XLALGPSToINT8NS(&segmentEnd);

        if ((startDiff > 0) && (endDiff < 0))
        {
            verbose("Generating analysis segment for injection at %d.\n",
                    injTime.gpsSeconds);
            if (*start)
            {
                verbose("warning: multiple injections in this segment.\n");
                *start = numPoints/4;
                *end = 3 * numPoints/4;
            }
            else
            {
                injDiff = (REAL8) ((XLALGPSToINT8NS(&injTime) - \
                                    XLALGPSToINT8NS(&segmentStart)) / 1E9);
                injSamplePoint = floor(injDiff * params->sampleRate + 0.5);
                injSamplePoint += numPoints/4;
                injWindow = floor(params->injSearchWindow * params->sampleRate
                                  + 1);
                *start = injSamplePoint - injWindow;
                if (*start < numPoints/4)
                    *start = numPoints/4;
                *end = injSamplePoint + injWindow + 1;
                if (*end > 3*numPoints/4)
                    *end = 3*numPoints/4;
                verbose("Found analysis segment at [%d,%d).\n", *start, *end);
            }
        }
        thisInject = thisInject->next;
    }
}
