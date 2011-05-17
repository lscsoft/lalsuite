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
      channel = get_simulated_data( ifoChannel, &params->startTime,
          params->duration, LALRINGDOWN_DATATYPE_SIM, params->sampleRate,
          params->randomSeed+ 100*ifoNumber, 1E-20 );
    else if ( params->zeroData )
    {
      channel = get_zero_data( ifoChannel, &params->startTime,
          params->duration, LALRINGDOWN_DATATYPE_ZERO, params->sampleRate );
    }
    else if ( params->doubleData )
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
    REAL8 *Fplus,
    REAL8 *Fcross,
    REAL8 *timeOffsets )
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
    /* compute raw average spectrum; store spectrum in invspec for now */
    invspec = compute_average_spectrum( channel,
        LALRINGDOWN_SPECTRUM_MEDIAN_MEAN, params->segmentDuration,
        params->strideDuration, fwdplan, params->whiteSpectrum );

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
    struct coh_PTF_params   *params
    )
{
  RingDataSegments *segments = NULL;
  COMPLEX8FrequencySeries  *response = NULL;
  UINT4  sgmnt,i;
  UINT4  segListToDo[params->numOverlapSegments];

  segments = LALCalloc( 1, sizeof( *segments ) );

  if ( params->analyzeInjSegsOnly )
  {
    segListToDo[2] = 1;
    for ( i = 0 ; i < params->numOverlapSegments; i++)
      segListToDo[i] = 0;
    SimInspiralTable        *injectList = NULL;
    SimInspiralTable        *thisInject = NULL;
    LIGOTimeGPS UNUSED injTime;
    REAL8 deltaTime;
    INT4 segNumber, UNUSED segLoc, UNUSED ninj;
    segLoc = 0;
    ninj = SimInspiralTableFromLIGOLw( &injectList, params->injectFile, params->startTime.gpsSeconds, params->startTime.gpsSeconds + params->duration );
    while (injectList)
    {
      injTime = injectList->geocent_end_time;
      deltaTime = injectList->geocent_end_time.gpsSeconds;
      deltaTime += injectList->geocent_end_time.gpsNanoSeconds * 1E-9;
      deltaTime -= params->startTime.gpsSeconds;
      deltaTime -= params->startTime.gpsNanoSeconds * 1E-9;
      segNumber = floor(2*(deltaTime/params->segmentDuration) - 0.5);
      segListToDo[segNumber] = 1;
      thisInject = injectList;
      injectList = injectList->next;
      LALFree( thisInject );
    }
  }

 /* FIXME: For all sky mode trig start/end time needs to be implemented */
 /* these probably work in the ring code! */

  /* if todo list is empty then do them all */
  if ( (! params->segmentsToDoList || ! strlen( params->segmentsToDoList )) && (! params->analyzeInjSegsOnly) )
  {
    segments->numSgmnt = params->numOverlapSegments;
    segments->sgmnt = LALCalloc( segments->numSgmnt, sizeof(*segments->sgmnt) );
    for ( sgmnt = 0; sgmnt < params->numOverlapSegments; ++sgmnt )
      compute_data_segment( &segments->sgmnt[sgmnt], sgmnt, channel, invspec,
          response, params->segmentDuration, params->strideDuration, fwdplan );
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
          compute_data_segment( &segments->sgmnt[count++], sgmnt, channel,
            invspec, response, params->segmentDuration, params->strideDuration,
            fwdplan );
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
  REAL4 a[LAL_NUM_IFO],
  REAL4 b[LAL_NUM_IFO],
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
          zh[i*vecLength+j] += a[k]*a[k] * PTFM[k]->data[i*PTFMlen+j];
          sh[i*vecLength+j] += b[k]*b[k] * PTFM[k]->data[i*PTFMlen+j];
          yu[i*vecLength+j] += a[k]*b[k] * PTFM[k]->data[i*PTFMlen+j];
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
      /*fprintf(stdout,"%f ",gsl_matrix_get(B2,i,j));*/
    }
    /*fprintf(stdout,"\n");*/
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
    REAL4 a[LAL_NUM_IFO],
    REAL4 b[LAL_NUM_IFO],
    INT4  timeOffsetPoints[LAL_NUM_IFO],
    gsl_matrix *eigenvecs,
    gsl_vector *eigenvals,
    UINT4 numPoints,
    UINT4 position,
    UINT4 vecLength,
    UINT4 vecLengthTwo)
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
    for( k = 0; k < LAL_NUM_IFO; k++)
    {
      if ( params->haveTrig[k] )
      {
        if (j < vecLength)
        {
          v1[j] += a[k] * PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]].re;
          v2[j] += a[k] * PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]].im;
        }
        else
        {
          v1[j] += b[k] * PTFqVec[k]->data[(j-vecLength)*numPoints+position+timeOffsetPoints[k]].re;
          v2[j] += b[k] * PTFqVec[k]->data[(j-vecLength)*numPoints+position+timeOffsetPoints[k]].im;
        }
      }
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
    REAL8                   *timeOffsets,
    REAL8                   *Fplus,
    REAL8                   *Fcross,
    REAL8                   *Fplustrig,
    REAL8                   *Fcrosstrig
    )
{
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
    if ( fcTmpltParams->PTFQ )
      XLALDestroyVectorSequence( fcTmpltParams->PTFQ );
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
    plan = XLALCreateForwardREAL4FFTPlan( segmentLength, 1 );
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
    plan = XLALCreateReverseREAL4FFTPlan( segmentLength, 1 );
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

/* construct list of of (ra, dec) pairs for multiple sky points */

void coh_PTF_generate_sky_points( 
    struct coh_PTF_skyPoints *skyPoints,
    struct coh_PTF_params *params
    )
{

  /* if all sky */
  if ( params->skyLooping == ALL_SKY || params->skyLooping == TWO_DET_ALL_SKY )
  {
    fprintf( stderr, "all sky mode is not implemented yet, please provide --declination" );
    exit(1);

    //skyPoints = coh_PTF_sky_grid()
  }

  /* if sky region */
  else if ( params->skyLooping == SKY_POINT_ERROR\
            || params->skyLooping == TWO_DET_SKY_POINT_ERROR
            || params->skyLooping == SINGLE_SKY_POINT )
  {
    coh_PTF_sky_grid( skyPoints, params );
  }

}

/*
 * Generate a grid of sky points based on the sinusoidal map method used by
 * the xpipeline.
 */

void coh_PTF_sky_grid(
    struct coh_PTF_skyPoints *skyPoints,
    struct coh_PTF_params *params    
    )
{
  UINT4 ifoNumber,i,j,k;
  LALDetector *detectors[LAL_NUM_IFO];
  REAL4 angle;   /* opening angle between 2 IFO baseline and sky localisation */
  REAL4 lambdamin,lambdamax,lambda; /* opening angle closest to pi/2 */
  REAL4 alpha,detalpha;
  alpha = 0;                         
  REAL4 angularResolution;          /* angular resolution of grid in radians */
  double baseline,lightTravelTime;
  REAL4 theta,phi;                  /* sky parameters */
  REAL4 raNp  = 0.;                 /* north */
  REAL4 decNp = LAL_PI_2;        /* pole */

  REAL4 rho;                        /* magnitude of rotation vector */
  REAL4 axis[3];                    /* rotation axis vector */
  REAL4 npPos[3];                   /* north pole position */
  REAL4 trigPos[3];                 /* trigger position */
  REAL4 pos[3],rotPos[3];           /* sky point position (orig & rotated) */

  /* get site coordinates */
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    detectors[ifoNumber] = LALCalloc( 1, sizeof( *detectors[ifoNumber] ) );
    XLALReturnDetector( detectors[ifoNumber], ifoNumber );
  }

  /* find pair of detectors whose opening angle to the GRB ±error is closest */
  /* to 90 degrees */
  for ( i = 0; i < LAL_NUM_IFO; i++ )
  {
    if ( params->haveTrig[i] )
    {
      for ( j = i+1; j < LAL_NUM_IFO; j++ )
      {
        if ( params->haveTrig[j] )
        {
          /* get dot product (time delay) between sites */
          baseline = XLALArrivalTimeDiff( detectors[i]->location,
                                          detectors[j]->location,
                                          params->rightAscension,
                                          params->declination,
                                          &params->trigTime );

          /* get light travel time */
          lightTravelTime = XLALLightTravelTime( detectors[i], detectors[j] );
          lightTravelTime *= 1e-9;

          /* calculate opening angle */
          angle  = acos( baseline/lightTravelTime );

          /* generate angular window with sky error */
          lambdamin = angle-params->skyError;
          lambdamax = angle+params->skyError;
   
          /* if pi/2 is in the range, choose that, 
           * otherwise get as close as possible */
          if ( lambdamin < LAL_PI_2 && lambdamax > LAL_PI_2 )
          {
            lambda = LAL_PI_2;
          }
          else
            if ( abs( LAL_PI_2 - lambdamin) < abs( LAL_PI_2 - lambdamax ) )
              lambda = lambdamin;
            else
              lambda = lambdamax;

          /* calculate alpha */
          detalpha = lightTravelTime * sin( lambda );
          if ( detalpha > alpha )
            alpha = detalpha;
    
        }
      }
    }
  }

  /* calculate angular resolution */
  angularResolution = 2. * params->timingAccuracy / alpha;

  /* generate sky grid using sinusoidal map */
  *skyPoints = coh_PTF_circular_grid( angularResolution, params->skyError );

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
  npPos[0] = cos(decNp) * cos(raNp);
  npPos[1] = cos(decNp) * sin(raNp);
  npPos[2] = sin(decNp);

  trigPos[0] = cos(params->declination) * cos(params->rightAscension);
  trigPos[1] = cos(params->declination) * sin(params->rightAscension);
  trigPos[2] = sin(params->declination);

  axis[0] = npPos[1]*trigPos[2] - npPos[2]*trigPos[1];
  axis[1] = npPos[2]*trigPos[0] - npPos[0]*trigPos[2];
  axis[2] = npPos[0]*trigPos[1] - npPos[1]*trigPos[0];

  rho = sqrt( pow(axis[0],2) + pow(axis[1],2) + pow(axis[2],2) );
  for ( i=0; i<3; i++ )
    axis[i] /= rho;

  /* generate rotation matrix */
  REAL4 R[3][3];

  R[0][0] = cos(angle) + pow(axis[0],2) * ( 1 - cos(angle) );
  R[0][1] = axis[0] * axis[1] * ( 1 - cos(angle) ) - axis[2] * sin(angle);
  R[0][2] = axis[0] * axis[2] * ( 1 - cos(angle) ) + axis[1] * sin(angle);

  R[1][0] = axis[1] * axis[0] * ( 1 - cos(angle) ) + axis[2] * sin(angle);
  R[1][1] = cos(angle) + pow(axis[1],2) * ( 1 - cos(angle) );
  R[1][2] = axis[1] * axis[2] * ( 1 - cos(angle) ) - axis[0] * sin(angle);

  R[2][0] = axis[2] * axis[0] * ( 1 - cos(angle) ) - axis[1] * sin(angle);
  R[2][1] = axis[2] * axis[1] * ( 1 - cos(angle) ) + axis[0] * sin(angle);
  R[2][2] = cos(angle) + pow(axis[2],2) * ( 1 - cos(angle) );

  /* loop over points rotating by angle around axis */
  for ( i=0; i < skyPoints->numPoints; i++ )
  {

    /* convert to spherical */
    phi   = skyPoints->rightAscension[i];
    theta = LAL_PI_2 - skyPoints->declination[i];
    
    /* convert to cartesian */
    pos[0] = sin(theta)*cos(phi);
    pos[1] = sin(theta)*sin(phi);
    pos[2] = cos(theta);


    /* rotate */
    for ( k=0; k<3; k++ )
    {
      rotPos[k] = R[k][0]*pos[0] + R[k][1]*pos[1] + R[k][2]*pos[2];
    }

    /* convert back to (phi,theta) */
    theta = acos(rotPos[2]);
    phi   = atan2(rotPos[1],rotPos[0]);
    skyPoints->rightAscension[i] = phi;
    skyPoints->declination[i]    = LAL_PI_2 - theta;

    if ( skyPoints->rightAscension[i] < 0.0 )
      skyPoints->rightAscension[i] += LAL_TWOPI;
    if ( skyPoints->rightAscension[i] >= LAL_TWOPI )
      skyPoints->rightAscension[i] -= LAL_TWOPI;

  }

  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if ( detectors[ifoNumber] )
      LALFree(detectors[ifoNumber]);
  }


}

/*
 * generate cicular map of sky points centred on north pole:
 *  * place central point
 *  * step out in declination by angularResolution
 *  * place a ring of points separated by angularResolution
 *  * repeat until maximum skyError radius passed
 */

struct coh_PTF_skyPoints coh_PTF_circular_grid(
    REAL4        angularResolution,
    REAL4        skyError
    )
{

  UINT4                    p               = 0;
  UINT4                    i,j;
  UINT4                    numTheta;            /* number of rings */
  REAL4                    dPhi,phi,theta;      /* sky point parameters */
  UINT4                    numSkyPoints    = 0; /* total number of sky points */
  struct coh_PTF_skyPoints skyPoints;

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
  skyPoints.rightAscension = LALCalloc(1, numSkyPoints*sizeof( REAL4 ));
  skyPoints.declination    = LALCalloc(1, numSkyPoints*sizeof( REAL4 ));
  skyPoints.numPoints      = numSkyPoints;

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
      skyPoints.rightAscension[p] = phi;
      skyPoints.declination[p]    = LAL_PI_2 - theta;

      if ( skyPoints.rightAscension[p] < 0.0 )
        skyPoints.rightAscension[p] += LAL_TWOPI;
      else if ( skyPoints.rightAscension[p] >= LAL_TWOPI )
        skyPoints.rightAscension[p] -= LAL_TWOPI;

      p++;
    }

  }

  return skyPoints;

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

  printf("%ld.%06ld", (long)tv->tv_sec, (long)tv->tv_usec);
  curtime = tv->tv_sec;
  strftime(buffer, 30, "%m-%d-%Y  %T", localtime(&curtime));
  printf(" = %s.%06ld\n", buffer, (long)tv->tv_usec);
}
