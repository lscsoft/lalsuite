#include "coh_PTF.h"

INT4 coh_PTF_data_condition(
              struct coh_PTF_params *params,
              REAL4TimeSeries          **channel,
              REAL4FrequencySeries     **invspec,
              RingDataSegments         **segments,
              REAL4FFTPlan             *fwdplan,
              REAL4FFTPlan             *psdplan,
              REAL4FFTPlan             *revplan,
              REAL4                    **timeSlideVectorsP,
              struct timeval           startTime
)
{
  REAL4 *timeSlideVectors;
  timeSlideVectors=LALCalloc(1, (LAL_NUM_IFO+1)*
      params->numOverlapSegments*sizeof(REAL4));
  memset(timeSlideVectors, 0,
      (LAL_NUM_IFO+1) * params->numOverlapSegments * sizeof(REAL4));


  UINT4 ifoNumber;
  INT4 numSegments = -1;
  /* loop over ifos */
  for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {

    /* if ifo is on: */
    if (params->haveTrig[ifoNumber])
    {
      /* Read in data from the various ifos */
      channel[ifoNumber] = coh_PTF_get_data(params, params->channel[ifoNumber],
                                            params->dataCache[ifoNumber],
                                            ifoNumber);
      coh_PTF_rescale_data(channel[ifoNumber], params->dynRangeFac);

      /* compute the spectrum */
      invspec[ifoNumber] = coh_PTF_get_invspec(channel[ifoNumber], fwdplan,
                                               revplan, psdplan, params);

      /* create the segments */

      segments[ifoNumber] = coh_PTF_get_segments(channel[ifoNumber],
                                invspec[ifoNumber],fwdplan, ifoNumber,
                                timeSlideVectors, params);

      if (numSegments < 0)
      {
        if (segments[ifoNumber])
        {
          numSegments = segments[ifoNumber]->numSgmnt;
        }
        else
        {
          numSegments = 0;
        }
      }
      else
      {
        if ( (! segments[ifoNumber]) )
        {
          if (numSegments != 0)
          {
            error("ERROR: Disagreement in number of segments in ifos.");
          }
        }
        else if (numSegments != (INT4)segments[ifoNumber]->numSgmnt)
        {
          error("ERROR: Disagreement in number of segments in ifos.");
        }
      }

      verbose("Created segments for one ifo %ld \n",
               timeval_subtract(&startTime));
    }
  }
  *timeSlideVectorsP = timeSlideVectors;
  return numSegments;
}




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

void coh_PTF_setup_null_stream(
    struct coh_PTF_params   *params,
    REAL4TimeSeries         **channel,
    REAL4FrequencySeries    **invspec,
    RingDataSegments        **segments,
    REAL4                   *Fplustrig,
    REAL4                   *Fcrosstrig,
    REAL4                   *timeOffsets,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *revplan,
    REAL4FFTPlan            *psdplan,
    REAL4                   *timeSlideVectors,
    struct timeval           startTime
)
{
  /* FIXME: The null stream probably will be broken by the timesliding stuff */
  /* It may be worth deprecating this entirely, but if not then a timeslid */
  /* null stream construction should be implemented. */
  /* My vote is to deprecate the null stream. */

  UINT4 ui;

  /* Read in data from the various ifos */
  if (coh_PTF_get_null_stream(params, channel, Fplustrig, Fcrosstrig,
                                timeOffsets))
  {
    error("Null stream construction failure\n");
  }

  /* compute the spectrum */
  invspec[LAL_NUM_IFO] = coh_PTF_get_invspec(channel[LAL_NUM_IFO], fwdplan,\
                                            revplan, psdplan, params);
  /* If white spectrum need to scale this. */
  if (params->whiteSpectrum)
  {
    for(ui=0 ; ui < invspec[LAL_NUM_IFO]->data->length; ui++)
    {
      /* FIXME: The factor here is hardcoded and wrong. */
      /* Should be dynamically calculated. */
      invspec[LAL_NUM_IFO]->data->data[ui] *= pow(1./0.3403324,2);
    }
  }

  /* create the segments */
  segments[LAL_NUM_IFO] = coh_PTF_get_segments(channel[LAL_NUM_IFO],\
                                        invspec[LAL_NUM_IFO],fwdplan,\
                                        LAL_NUM_IFO,timeSlideVectors, params);

  verbose("Created segments for null stream at %ld \n",
          timeval_subtract(&startTime));
}

/* This function is used to generate the null stream */
/* Be aware this is separate from the null SNR! */
int coh_PTF_get_null_stream(
    struct coh_PTF_params *params,
    REAL4TimeSeries **channel,
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
  coh_PTF_convert_time_offsets_to_points(params,timeOffsets,timeOffsetPoints);

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
    REAL4FFTPlan            *psdplan,
    struct coh_PTF_params   *params
    )
{
  REAL4FrequencySeries *invspec = NULL;
  if ( params->getSpectrum )
  {
    if ( ! params->whiteSpectrum)
    {
      UINT4 recordLength,psdSegmentLength,segmentStride,numDoublePsdSegs;
      UINT4 neededDataLength;
      /* compute raw average spectrum; store spectrum in invspec for now */
      REAL4FrequencySeries *invspecorig = NULL;
      /* If the data length is no appropriate, don't use all of it */
      recordLength = params->duration * params->sampleRate;
      psdSegmentLength = floor(params->psdSegmentDuration * \
          params->sampleRate + 0.5);
      segmentStride = floor(params->psdStrideDuration * params->sampleRate \
          + 0.5);
      sanity_check( segmentStride > 0 );
      numDoublePsdSegs = 1 + (recordLength - psdSegmentLength)/ \
          (segmentStride);
      numDoublePsdSegs /= 2;
      // Need enough psd segments for good estimation
      sanity_check(numDoublePsdSegs > 5);
      neededDataLength = ((2*numDoublePsdSegs)-1) * segmentStride + psdSegmentLength;
      verbose("Using %e of %e seconds for PSD estimation in %d segments.\n",\
          neededDataLength/params->sampleRate,params->duration,\
          2*numDoublePsdSegs);
      // Adjust channel length for PSD calculation
      channel->data->length = neededDataLength;
      invspecorig = compute_average_spectrum( channel,
          LALRINGDOWN_SPECTRUM_MEDIAN_MEAN, params->psdSegmentDuration,
          params->psdStrideDuration, psdplan, 0 );
      // Set the channel length back
      channel->data->length = recordLength;
      /* Resample if necessary */
      if (params->psdSegmentDuration != params->segmentDuration)
      {
        if ( params->writeInvSpectrum ) /* Write spectrum before resampling */
          write_REAL4FrequencySeries( invspecorig );
        verbose("Resampling PSD\n");
        invspec = resample_psd(invspecorig, params->sampleRate,
                               params->segmentDuration);
      }
      else
      {
        invspec = invspecorig;
      } 
    }
    else
    {
      UINT4 k;
      REAL8FrequencySeries *spectrum;
      spectrum = generate_theoretical_psd(1./params->sampleRate,
          params->segmentDuration,params->simDataType, 1E-20);
  
      /* Need to convert to a REAL4 FrequencySeries */
      invspec = XLALCreateREAL4FrequencySeries("TEMP",&(channel->epoch),0,
          1.0/params->segmentDuration,&lalDimensionlessUnit,
          params->numFreqPoints);
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
    invert_spectrum( invspec, params->sampleRate, params->strideDuration,
        params->truncateDuration, params->highpassFrequency, fwdplan,
        revplan );

    if ( params->writeInvSpectrum ) /* write inverse calibrated spectrum */
      write_REAL4FrequencySeries( invspec );
  }

  return invspec;
}

/* Function to rescale the data to avoid floating point errors*/
void coh_PTF_rescale_data (REAL4TimeSeries *channel, REAL8 rescaleFactor)
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
    /* NOTE: Using this flag will override anything given by the option
     * --only-analyse-segments, do not try and use both together */ 
    /* Figure out which segments the injections are in. If an injection is
     * within 1s of a segment boundary, analyse both segments. */ 
    for ( i = 0 ; i < params->numOverlapSegments; i++)
      segListToDo[i] = 0;
    SimInspiralTable        *injectList = NULL;
    REAL8 deltaTime,segBoundDiff;
    INT4 segNumber, UNUSED ninj;
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
      /* Calculate the epoch of the injection relative to the gps-start */
      deltaTime = injectList->geocent_end_time.gpsSeconds;
      deltaTime += injectList->geocent_end_time.gpsNanoSeconds * 1E-9;
      deltaTime -= params->startTime.gpsSeconds;
      deltaTime -= params->startTime.gpsNanoSeconds * 1E-9;

      /* Adjust to the epoch relative to the analysis start */
      deltaTime -= params->analStartTime;
    
      /* And figure out the segment number */
      segNumber = floor( deltaTime / params->strideDuration);
      if (segNumber >= (INT4) params->numOverlapSegments)
      {
        verbose("Injection at %" LAL_INT8_FORMAT " after analysis window.\n",\
                injectList->geocent_end_time.gpsSeconds);
      }
      else if (segNumber < 0)
      {
        verbose("Injection at %" LAL_INT8_FORMAT " before analysis window.\n",\
                injectList->geocent_end_time.gpsSeconds);
      }
      else
      {
        segListToDo[segNumber] = 1;
      }
      /* Check if injection is near a segment boundary */
      for ( j = 0 ; j < params->numOverlapSegments; j++)
      {
        segBoundDiff = deltaTime - j * params->strideDuration/2;
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
    /* Not convinved the code ever goes here, because empty segmentsToDoList
     * does not evaluate to false. The loop below still works though! */ 
    segments->numSgmnt = params->numOverlapSegments;
    segments->sgmnt = LALCalloc( segments->numSgmnt, sizeof(*segments->sgmnt) );
    for ( sgmnt = 0; sgmnt < params->numOverlapSegments; ++sgmnt )
    {
      slidSegNum = ( sgmnt + ( params->slideSegments[NumberIFO] ) ) % ( segments->numSgmnt );
      timeSlideVectors[NumberIFO*params->numOverlapSegments + sgmnt] = 
          ((INT4)slidSegNum-(INT4)sgmnt)*params->strideDuration;
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
    {
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
    }

    if ( ! count ) /* no segments to do */
    {
      LALFree(segments);
      return NULL;
    }

    segments->numSgmnt = count;
    segments->sgmnt = LALCalloc( segments->numSgmnt, sizeof(*segments->sgmnt) );

    count = 0;
    for ( sgmnt = 0; sgmnt < params->numOverlapSegments; ++sgmnt )
    {
      if ( params->analyzeInjSegsOnly )
      {
        /* NOTE: Time slide injection analysis is not currently possible */
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
          slidSegNum = ( sgmnt + ( params->slideSegments[NumberIFO] ) ) % ( params->numOverlapSegments );
          timeSlideVectors[NumberIFO*params->numOverlapSegments + count] =
              ((INT4)slidSegNum-(INT4)sgmnt)*params->strideDuration;
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

void coh_PTF_create_time_slide_table(
  struct coh_PTF_params   *params,
  INT8                    *slideIDList,
  RingDataSegments        **segments,
  TimeSlide               **time_slide_headP,
  TimeSlideSegmentMapTable **time_slide_map_headP,
  SegmentTable            **segment_table_headP,
  TimeSlideVectorList     **longTimeSlideListP,
  TimeSlideVectorList     **shortTimeSlideListP,
  REAL4                   *timeSlideVectors,
  INT4                    numSegments
)
{
  TimeSlide *time_slide_head = NULL;
  TimeSlide *curr_slide = NULL;
  TimeSlideSegmentMapTable *time_slide_map_head = NULL;
  TimeSlideSegmentMapTable *curr_time_slide_map = NULL;
  SegmentTable *segment_table_head = NULL;
  SegmentTable *curr_slide_segment = NULL;
  TimeSlideVectorList *longTimeSlideList = LALCalloc(\
      numSegments, sizeof(TimeSlideVectorList));
  CHAR ifoName[LIGOMETA_IFO_MAX];
  UINT4 longSlideCount = 0;
  UINT4 shortSlideCount = 0;
  INT4 i;
  UINT4 ui,uj,ifoNumber,ifoNum,lastStartPoint,wrapPoint,currId,ifoCount;
  UINT4 slideSegmentCount, calculatedFlag;
  REAL8 currBaseOffset,currBaseIfoOffset[LAL_NUM_IFO],wrapTime,currIfoOffset;
  LIGOTimeGPS slideStartTime,slideEndTime, tmpSlideStartTime, tmpSlideEndTime;

  /* First construct the list of long slides */
  for (i = 0 ; i < numSegments ; i++)
  {
    UINT4 slideDuplicate = 0;
    if (longSlideCount)
    {
      for (uj = 0; uj < longSlideCount; uj++)
      {
        UINT4 slideChecking = 1;
        /* Here we check if the time-slid segment i is the same slide as the
         * timeSlideList[uj] */  
        for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
        {
          if (params->haveTrig[ifoNumber])
          {
            if (timeSlideVectors[ifoNumber*params->numOverlapSegments+i] != \
                longTimeSlideList[uj].timeSlideVectors[ifoNumber])
            {
              slideChecking = 0;
            }
          }
        }
        if (slideChecking)
        {
          slideDuplicate = 1;
          slideIDList[i] = longTimeSlideList[uj].timeSlideID;
        }
      }
    }
    if (! slideDuplicate)
    {
      for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if (params->haveTrig[ifoNumber])
        {
          longTimeSlideList[longSlideCount].timeSlideVectors[ifoNumber] = 
              timeSlideVectors[ifoNumber*params->numOverlapSegments+i];
        }
      }
      longTimeSlideList[longSlideCount].timeSlideID = longSlideCount;
      slideIDList[i] = longTimeSlideList[longSlideCount].timeSlideID;
      longSlideCount++;
    }
  }

  TimeSlideVectorList *shortTimeSlideList = LALCalloc(\
      params->numShortSlides, sizeof(TimeSlideVectorList));

  /* Next construct the list of short time slides */
  /* Always have a zero-lag short slide */
  for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if (params->haveTrig[ifoNumber])
    {
      shortTimeSlideList[0].timeSlideVectors[ifoNumber] = 0.;
    }
  }
  shortTimeSlideList[0].timeSlideID = 0;
  shortTimeSlideList[0].analStartPoint = params->analStartPoint;
  shortTimeSlideList[0].analEndPoint = params->analEndPoint;
  shortSlideCount = 1;

  if (params->doShortSlides)
  {
    currBaseOffset = params->shortSlideOffset;
    while (1)
    {
      /* Calculate offsets for each detector */
      ifoNum = 0;
      for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        currBaseIfoOffset[ifoNumber] = 0;
        if (params->haveTrig[ifoNumber])
        {
          currBaseIfoOffset[ifoNumber] = currBaseOffset * ifoNum;
          /* If the offset is bigger than the stride then this cannot be done*/
          /* FIXME: We should check if the offset is not bigger than the */
          /* stride minus the BaseOffset, otherwise we can end up with short */
          /* slides whose offsets are very similar to the long slides */
          if (currBaseIfoOffset[ifoNumber] > \
                  (params->strideDuration) )
          {
            break;
          }
          ifoNum += 1;
        }
      }  
      if (ifoNumber != LAL_NUM_IFO)
      { /* If I did not finish the previous for loop due to break, break here*/
        /* THIS IS THE ONLY PLACE THAT I MIGHT BREAK FROM THE WHILE LOOP */
        break;
      }

      /* This avoids "WARNING: May be used unset" warning */
      lastStartPoint = 0;
      /* Now we need to consider the wrapping and set slides */
      for (ui = 0; ui < ifoNum; ui++)
      {
        /* Begin by initializing the slide */
        shortTimeSlideList[shortSlideCount].timeSlideID = shortSlideCount;

        /* Next we populate the slide's start and end times */
        /* We're actually going to loop backwards, so ui = 0 corresponds to the
         * end of the segment where are slides will have wrapped */  
        if (ui == 0)
        {
          /* If ui=0 we know the end point of slide will be analEndPoint */
          shortTimeSlideList[shortSlideCount].analEndPoint = \
               params->analEndPoint;
        }
        else
        {
          /* Otherwise end at where we started before */
          shortTimeSlideList[shortSlideCount].analEndPoint = lastStartPoint;
        } 
        if (ui == ifoNum - 1)
        {
          /* If ui = ifoNum we have to start at analStartPoint */
          shortTimeSlideList[shortSlideCount].analStartPoint = \
              params->analStartPoint;
        }
        else
        {
          /* Now the tricky one, we have to figure out where the next (previous
           * as we're looping backwards) wrap point will be. */  
          wrapTime = params->strideDuration - (ui + 1) * currBaseOffset;
          wrapPoint = floor(wrapTime * params->sampleRate + 0.5);
          /* Start point is then start point + wrap point */
          shortTimeSlideList[shortSlideCount].analStartPoint = \
              params->analStartPoint + wrapPoint;
        }
        lastStartPoint = shortTimeSlideList[shortSlideCount].analStartPoint;
        /* Now we construct the wrapped offsets and store*/
        ifoCount = 0;
        for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
        {
          /* We are looping backwards so invert this */
          if (params->haveTrig[ifoNumber])
          {
            currIfoOffset = currBaseIfoOffset[ifoNumber];
            if (ifoCount > ui)
            { /* Need to wrap */
              currIfoOffset -= params->strideDuration;
            }
            shortTimeSlideList[shortSlideCount].timeSlideVectors[ifoNumber] =
                currIfoOffset;
            ifoCount++;
          }
        }
        shortSlideCount++;
      } /* End loop over wrap points */
      currBaseOffset += params->shortSlideOffset;
    } /* End the while(1) loop over base offsets */
  } /* End if do short slides */
  sanity_check(shortSlideCount == params->numShortSlides);

  /* Time slide table will contain every long+short slide combination */
  for (ui = 0 ; ui < longSlideCount; ui++)
  { /* Loop over long slides */
    for (uj = 0 ; uj < shortSlideCount; uj++)
    { /* Loop over short slides */
      /* Construct the ID from the current long and short slide */
      currId = longTimeSlideList[ui].timeSlideID*shortSlideCount;
      currId += shortTimeSlideList[uj].timeSlideID;
      /* Create an entry in the time_slide_map */
      if (! time_slide_map_head)
      {
        time_slide_map_head=XLALCreateTimeSlideSegmentMapTableRow();
        curr_time_slide_map = time_slide_map_head;
      }
      else
      {
        curr_time_slide_map->next = XLALCreateTimeSlideSegmentMapTableRow();
        curr_time_slide_map = curr_time_slide_map->next;
      }
      /* Set the ID in the time_slide_table */
      curr_time_slide_map->time_slide_id = currId;
      curr_time_slide_map->segment_def_id = currId;

      for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if (params->haveTrig[ifoNumber])
        {
          /* We need an entry for every ifo and every slide */
          if (! time_slide_head)
          { /* Create table if it doesn't exist */
            time_slide_head=XLALCreateTimeSlide();
            curr_slide= time_slide_head;
          }
          else
          { /* Or setup the next entry in the linked list */
            curr_slide->next=XLALCreateTimeSlide();
            curr_slide = curr_slide->next;
          }
          /* Set the ID in the time_slide_table */
          curr_slide->time_slide_id = currId;
          /* Set the IFO  */
          XLALReturnIFO(ifoName,ifoNumber);
          strncpy(curr_slide->instrument,ifoName,\
                  sizeof(curr_slide->instrument)-1);
          /* Set the offset as the sum of the short and long offests */
          curr_slide->offset = \
              longTimeSlideList[ui].timeSlideVectors[ifoNumber];
          curr_slide->offset += \
              shortTimeSlideList[uj].timeSlideVectors[ifoNumber];
          /* Process IDs will be 0 */
          curr_slide->process_id=0;
        }
      } /* End ifo loop */
    } /* End short slide loop */
  } /* End long slide loop */

  /* Now we construct the list of analysed segments for each time slide_id
   * This will be a segment for *every* segment+short_slide combination */ 
  slideSegmentCount = 0;
  for (i = 0 ; i < numSegments ; i++)
  { /* Loop over segments */
    for (uj = 0 ; uj < shortSlideCount; uj++)
    { /* Loop over short slides */
      /* Contruct the ID form the current long+short slide */
      currId = slideIDList[i]*shortSlideCount;
      currId += shortTimeSlideList[uj].timeSlideID;
      /* Construct the segment start and end times */
      calculatedFlag = 0;
      for (ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if (params->haveTrig[ifoNumber])
        {
          tmpSlideStartTime = segments[ifoNumber]->sgmnt[i].epoch; 
          tmpSlideEndTime = segments[ifoNumber]->sgmnt[i].epoch;
          XLALGPSAdd(&tmpSlideStartTime, \
                     timeSlideVectors[ifoNumber*params->numOverlapSegments+i]);
          XLALGPSAdd(&tmpSlideEndTime, \
                     timeSlideVectors[ifoNumber*params->numOverlapSegments+i]);
          if (calculatedFlag)
          {
            if (! XLALGPSCmp(&tmpSlideStartTime, &slideStartTime))
            {
              error("Error long slide epochs not agreeing with the slid times. This suggests broken code, please contact a developer.");

            }
            if (! XLALGPSCmp(&tmpSlideEndTime, &slideEndTime))
            {
              error("Error long slide epochs not agreeing with the slid times. This suggests broken code, please contact a developer.");
            }
          }
          else
          {
            slideStartTime.gpsSeconds=tmpSlideStartTime.gpsSeconds;
            slideStartTime.gpsNanoSeconds=tmpSlideStartTime.gpsNanoSeconds;
            slideEndTime.gpsSeconds=tmpSlideEndTime.gpsSeconds;
            slideEndTime.gpsNanoSeconds=tmpSlideEndTime.gpsNanoSeconds;
            calculatedFlag = 1;
          }
        }
      }
      XLALGPSAdd(&slideStartTime, \
          ((REAL8)shortTimeSlideList[uj].analStartPoint) / params->sampleRate );
      XLALGPSAdd(&slideEndTime, \
          ((REAL8)shortTimeSlideList[uj].analEndPoint) / params->sampleRate );
      /* Add this to the segment table */
      if (! segment_table_head)
      { /* Create table if it doesn't exist */
        segment_table_head=XLALCreateSegmentTableRow(NULL);
        curr_slide_segment= segment_table_head;
      }
      else
      { /* Or setup the next entry in the linked list */
        curr_slide_segment->next=XLALCreateSegmentTableRow(NULL);
        curr_slide_segment = curr_slide_segment->next;
      }
      curr_slide_segment->segment_id = slideSegmentCount;
      curr_slide_segment->segment_def_id = currId;
      curr_slide_segment->start_time = slideStartTime;
      curr_slide_segment->end_time = slideEndTime;
      curr_slide_segment->process_id = 0;
      slideSegmentCount++;
    }
  }

  *time_slide_headP = time_slide_head;
  *time_slide_map_headP = time_slide_map_head;
  *segment_table_headP = segment_table_head;
  *shortTimeSlideListP = shortTimeSlideList;
  *longTimeSlideListP = longTimeSlideList;
}

void coh_PTF_calculate_det_stuff(
  struct coh_PTF_params   *params,
  LALDetector             **detectors,
  REAL4                   *timeOffsets,
  REAL4                   *Fplustrig,
  REAL4                   *Fcrosstrig,
  CohPTFSkyPositions      *skyPoints,
  UINT4                   skyPointNum
)
{
  UINT4 ifoNumber,ui;
  REAL8 FplusTmp,FcrossTmp;
  REAL8 detLoc[3];

  for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if (! detectors[ifoNumber])
    {
      detectors[ifoNumber] = LALCalloc(1, sizeof(*detectors[ifoNumber]));
      XLALReturnDetector(detectors[ifoNumber], ifoNumber);
    }
    /* get location in three dimensions */
    for (ui = 0; ui < 3; ui++)
    {
      detLoc[ui] = (double) detectors[ifoNumber]->location[ui];
    }
    /* calculate time offsets */
    timeOffsets[ifoNumber] = (REAL4)
        XLALTimeDelayFromEarthCenter(detLoc,\
                                     skyPoints->data[skyPointNum].longitude,
                                     skyPoints->data[skyPointNum].latitude,
                                     &(params->trigTime));
    /* calculate response functions for trigger */
    XLALComputeDetAMResponse(&FplusTmp, &FcrossTmp,
                            (const REAL4 (*)[3])detectors[ifoNumber]->response,
                            skyPoints->data[skyPointNum].longitude,
                            skyPoints->data[skyPointNum].latitude, 0.,
                            XLALGreenwichMeanSiderealTime(&(params->trigTime)));
    Fplustrig[ifoNumber] = (REAL4) FplusTmp;
    Fcrosstrig[ifoNumber] = (REAL4) FcrossTmp;
  }
}

void coh_PTF_initialize_structures(
  struct coh_PTF_params    *params,
  FindChirpInitParams      **fcInitParamsP,
  FindChirpTemplate        **fcTmpltP,
  FindChirpTmpltParams     **fcTmpltParamsP,
  REAL8Array               **PTFM,
  REAL8Array               **PTFN,
  COMPLEX8VectorSequence   **PTFqVec,
  REAL4FFTPlan             *fwdplan
)
{
  UINT4 ifoNumber;

  FindChirpInitParams *fcInitParams;
  FindChirpTemplate *fcTmplt;
  FindChirpTmpltParams *fcTmpltParams;

  coh_PTF_set_null_input_REAL8Array(PTFM,LAL_NUM_IFO+1);
  coh_PTF_set_null_input_REAL8Array(PTFN,LAL_NUM_IFO+1);
  coh_PTF_set_null_input_COMPLEX8VectorSequence(PTFqVec,LAL_NUM_IFO+1);

  /* finchirp parameters */
  fcInitParams               = LALCalloc(1, sizeof(*fcInitParams));
  fcTmplt                    = LALCalloc(1, sizeof(*fcTmplt));
  fcTmpltParams              = LALCalloc(1, sizeof(*fcTmpltParams));
  fcTmpltParams->approximant = params->approximant;
  fcTmpltParams->order       = params->order;

  /* Note that although non-spinning only uses Q1, the PTF
 *   generator still generates Q1-5, thus size of these vectors */
  if (params->approximant == FindChirpPTF)
  {
    fcTmplt->PTFQtilde          =
        XLALCreateCOMPLEX8VectorSequence(5, params->numFreqPoints);
    fcTmplt->PTFQ         = XLALCreateVectorSequence(5, params->numTimePoints);
    fcTmpltParams->PTFphi = XLALCreateVector(params->numTimePoints);
    fcTmpltParams->PTFomega_2_3 = XLALCreateVector(params->numTimePoints);
    fcTmpltParams->PTFe1  = XLALCreateVectorSequence(3, params->numTimePoints);
    fcTmpltParams->PTFe2  = XLALCreateVectorSequence(3, params->numTimePoints);
  }
  else if (params->approximant == FindChirpSP)
  {
    fcTmplt->PTFQtilde      =
        XLALCreateCOMPLEX8VectorSequence(1, params->numFreqPoints);
    fcTmplt->data           = XLALCreateCOMPLEX8Vector(params->numFreqPoints);
    fcTmpltParams->xfacVec  = XLALCreateVector(params->numFreqPoints);
    fcTmpltParams->PTFphi   = XLALCreateVector(params->numFreqPoints);
    /* Set the values of xfacVec  This is k^(-1/3) 
     * also PTFphi which is k^(-7/6). As these are expensive to compute it
     * is cheaper to do it once rather than many times. */
    const REAL4                   xfacExponent = -1.0/3.0;
    REAL4                        *xfac = NULL;
    UINT4 ui;
    xfac = fcTmpltParams->xfacVec->data;
    xfac[0] = 0;
    fcTmpltParams->PTFphi->data[0] = 0;
    for (ui = 1; ui < fcTmpltParams->xfacVec->length; ++ui)
    {
      xfac[ui] = pow((REAL4) ui, xfacExponent);
      fcTmpltParams->PTFphi->data[ui] = pow((REAL4) ui, -7.0/6.0);
    }
  }
  else
  {
    fcTmplt->PTFQtilde      =
        XLALCreateCOMPLEX8VectorSequence(1, params->numFreqPoints);
    fcTmplt->data           = XLALCreateCOMPLEX8Vector(params->numFreqPoints);
    fcTmpltParams->xfacVec  = XLALCreateVector(params->numTimePoints);
  }

  fcTmpltParams->fwdPlan      = fwdplan;
  fcTmpltParams->deltaT       = 1.0/params->sampleRate;
  if (params->dynTempLength)
  {
    fcTmpltParams->dynamicTmpltFlow = 1;
  }
  else
  {
    fcTmpltParams->dynamicTmpltFlow = 0;
    fcTmpltParams->fLow = params->lowTemplateFrequency;
  }
  /* This option holds 2x the length of data that is junk in each segment
   * because of conditioning and the PSD.*/
  fcTmpltParams->invSpecTrunc = params->truncateDuration;
  fcTmpltParams->maxTempLength = params->maxTempLength;
  fcTmpltParams->dynRange = params->dynRangeFac;

  for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if (params->haveTrig[ifoNumber])
    {
      if (params->approximant == FindChirpPTF)
      {
        PTFM[ifoNumber]    = XLALCreateREAL8ArrayL(2, 5, 5);
        PTFN[ifoNumber]    = XLALCreateREAL8ArrayL(2, 5, 5);
        PTFqVec[ifoNumber] = XLALCreateCOMPLEX8VectorSequence\
                             (5, params->numTimePoints);
      }
      else
      {
        PTFM[ifoNumber]    = XLALCreateREAL8ArrayL(2, 1, 1);
        PTFN[ifoNumber]    = XLALCreateREAL8ArrayL(2, 1, 1);
        PTFqVec[ifoNumber] = XLALCreateCOMPLEX8VectorSequence\
                             (1, params->numTimePoints);
      }
    }
  }

  if (params->doNullStream)
  {
    if (params->approximant == FindChirpPTF)
    {
      PTFM[LAL_NUM_IFO]    = XLALCreateREAL8ArrayL(2, 5, 5);
      PTFN[LAL_NUM_IFO]    = XLALCreateREAL8ArrayL(2, 5, 5);
      PTFqVec[LAL_NUM_IFO] = XLALCreateCOMPLEX8VectorSequence\
                             (5, params->numTimePoints);
    }
    else
    {
      PTFM[LAL_NUM_IFO]    = XLALCreateREAL8ArrayL(2, 1, 1);
      PTFN[LAL_NUM_IFO]    = XLALCreateREAL8ArrayL(2, 1, 1);
      PTFqVec[LAL_NUM_IFO] = XLALCreateCOMPLEX8VectorSequence\
                             (1, params->numTimePoints);
    }
  }  
  /* Send the output back to the parent function */
  *fcTmpltP = fcTmplt;
  *fcInitParamsP = fcInitParams;
  *fcTmpltParamsP = fcTmpltParams; 

}

void coh_PTF_initialize_time_series(
  struct coh_PTF_params    *params,
  LIGOTimeGPS              segStartTime,
  REAL8                    fLower,
  REAL4TimeSeries          **cohSNRP,
  REAL4TimeSeries          **nullSNRP,
  REAL4TimeSeries          **traceSNRP,
  REAL4TimeSeries          **bankVeto,
  REAL4TimeSeries          **autoVeto,
  REAL4TimeSeries          **chiSquare,
  REAL4TimeSeries          **snrComps,
  REAL4TimeSeries          **pValues,
  REAL4TimeSeries          **gammaBeta,
  UINT4                    spinTemplates
)
{
  UINT4 ifoNumber,vecLength,vecLengthTwo,ui;
  REAL4TimeSeries *cohSNR,*nullSNR,*traceSNR;
  char name[LALNameLength];
  CHAR ifoName[LIGOMETA_IFO_MAX];
  cohSNR = NULL;
  nullSNR = NULL;
  traceSNR = NULL;

  /* Generate the various time series as needed*/
  /* We only want to store data from middle half of segment */
  cohSNR = XLALCreateREAL4TimeSeries("cohSNR", &segStartTime,\
                                     fLower,(1.0/params->sampleRate),\
                                     &lalDimensionlessUnit,\
                                     params->numAnalPoints);

  if (params->doNullStream)
    nullSNR = XLALCreateREAL4TimeSeries("nullSNR", &segStartTime,\
                                        fLower,(1.0/params->sampleRate),\
                                        &lalDimensionlessUnit,\
                                        params->numAnalPoints);

  if (params->doTraceSNR)
    traceSNR = XLALCreateREAL4TimeSeries("traceSNR", &segStartTime,\
                                         fLower,(1.0/params->sampleRate),\
                                         &lalDimensionlessUnit,\
                                         params->numAnalPoints);

  for (ifoNumber = 0;ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if (params->haveTrig[ifoNumber])
    {
      XLALReturnIFO(ifoName,ifoNumber);
      snprintf( name, sizeof( name ), "%s_snr",ifoName);
      snrComps[ifoNumber] = XLALCreateREAL4TimeSeries(name,
                            &segStartTime,fLower,(1.0/params->sampleRate),
                            &lalDimensionlessUnit,params->numAnalPointsBuf);
    }
  }

  if (params->doBankVeto)
  {
    if (params->numIFO != 1)
    {
      bankVeto[LAL_NUM_IFO] = XLALCreateREAL4TimeSeries("bank_veto",\
                                         &segStartTime,\
                                         fLower,(1.0/params->sampleRate),\
                                         &lalDimensionlessUnit,\
                                         params->numAnalPoints);
    }
    if (params->doSnglChiSquared)
    {
      for (ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if (params->haveTrig[ifoNumber])
        {
          XLALReturnIFO(ifoName,ifoNumber);
          snprintf( name, sizeof( name ), "%s_bank_veto",ifoName);
          bankVeto[ifoNumber] = XLALCreateREAL4TimeSeries(name,&segStartTime,\
                                         fLower,(1.0/params->sampleRate),\
                                         &lalDimensionlessUnit,\
                                         params->numAnalPoints);
        }
      }
    }
  } 

  if (params->doAutoVeto)
  {
    if (params->numIFO != 1)
    {
      autoVeto[LAL_NUM_IFO] = XLALCreateREAL4TimeSeries("auto_veto",\
                                         &segStartTime,\
                                         fLower,(1.0/params->sampleRate),\
                                         &lalDimensionlessUnit,\
                                         params->numAnalPoints);
    }
    if (params->doSnglChiSquared)
    {
      for (ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if (params->haveTrig[ifoNumber])
        {
          XLALReturnIFO(ifoName,ifoNumber);
          snprintf( name, sizeof( name ), "%s_auto_veto",ifoName);
          autoVeto[ifoNumber] = XLALCreateREAL4TimeSeries(name,&segStartTime,\
                                         fLower,(1.0/params->sampleRate),\
                                         &lalDimensionlessUnit,\
                                         params->numAnalPoints);
        }
      }
    }
  }

  if (params->doChiSquare)
  {
    if (params->numIFO != 1)
    {
      chiSquare[LAL_NUM_IFO] = XLALCreateREAL4TimeSeries("chi_square",\
                                         &segStartTime,\
                                         fLower,(1.0/params->sampleRate),\
                                         &lalDimensionlessUnit,\
                                         params->numAnalPoints);
    }
    if (params->doSnglChiSquared)
    {
      for (ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if (params->haveTrig[ifoNumber])
        {
          XLALReturnIFO(ifoName,ifoNumber);
          snprintf( name, sizeof( name ), "%s_chi_square",ifoName);
          chiSquare[ifoNumber] = XLALCreateREAL4TimeSeries(name,&segStartTime,\
                                         fLower,(1.0/params->sampleRate),\
                                         &lalDimensionlessUnit,\
                                         params->numAnalPoints);
        }
      }
    }
  }

  /* Work out how many pValue arrays are needed */
  if (spinTemplates)
    vecLength = 5;
  else
    vecLength = 2;
  if (params->numIFO == 1 || params->singlePolFlag || params->faceOnStatistic)
    vecLengthTwo = vecLength;
  else
    vecLengthTwo = 2* vecLength;

  for (ui = 0 ; ui < vecLengthTwo ; ui++)
  {
    pValues[ui] = XLALCreateREAL4TimeSeries("Pvalue", &segStartTime,\
                                             fLower, (1.0/params->sampleRate),\
                                             &lalDimensionlessUnit,\
                                             params->numAnalPoints);
  }

  if (spinTemplates)
  {
    for (ui = 0 ; ui < 2 ; ui++)
    {
      gammaBeta[ui] = XLALCreateREAL4TimeSeries("Pvalue", &segStartTime,\
                                             fLower, (1.0/params->sampleRate),\
                                             &lalDimensionlessUnit,\
                                             params->numAnalPoints);
    }
  }

  *cohSNRP = cohSNR;
  *nullSNRP = nullSNR;
  *traceSNRP = traceSNR;

}

void coh_PTF_reset_time_series(
  struct coh_PTF_params    *params,
  LIGOTimeGPS              segStartTime,
  REAL4TimeSeries          *cohSNR,
  REAL4TimeSeries          *nullSNR,
  REAL4TimeSeries          *traceSNR,
  REAL4TimeSeries          **bankVeto,
  REAL4TimeSeries          **autoVeto,
  REAL4TimeSeries          **chiSquare,
  REAL4TimeSeries          **snrComps,
  REAL4TimeSeries          **pValues,
  REAL4TimeSeries          **gammaBeta,
  UINT4                    spinTemplates
)
{
  UINT4 ifoNumber,vecLength,vecLengthTwo,ui;

  /* For each of the time series we reset the epoch*/
  cohSNR->epoch.gpsSeconds = segStartTime.gpsSeconds;
  cohSNR->epoch.gpsNanoSeconds = segStartTime.gpsNanoSeconds;
  /* And memset to 0 all values */
  memset(cohSNR->data->data, 0, params->numAnalPoints * sizeof(REAL4));

  for (ifoNumber = 0;ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if (params->haveTrig[ifoNumber])
    {
      snrComps[ifoNumber]->epoch.gpsSeconds = segStartTime.gpsSeconds;
      snrComps[ifoNumber]->epoch.gpsNanoSeconds = \
              segStartTime.gpsNanoSeconds;
      memset(snrComps[ifoNumber]->data->data, 0,\
              params->numAnalPointsBuf * sizeof(REAL4));
    }
  }


  if (params->doNullStream)
  {
    nullSNR->epoch.gpsSeconds = segStartTime.gpsSeconds;
    nullSNR->epoch.gpsNanoSeconds = segStartTime.gpsNanoSeconds;
    memset(nullSNR->data->data, 0, params->numAnalPoints * sizeof(REAL4));
  }

  if (params->doTraceSNR)
  {
    traceSNR->epoch.gpsSeconds = segStartTime.gpsSeconds;
    traceSNR->epoch.gpsNanoSeconds = segStartTime.gpsNanoSeconds;
    memset(traceSNR->data->data, 0, params->numAnalPoints * sizeof(REAL4));
  }

  if (params->doBankVeto)
  {
    if (params->numIFO != 1)
    {
      bankVeto[LAL_NUM_IFO]->epoch.gpsSeconds = segStartTime.gpsSeconds;
      bankVeto[LAL_NUM_IFO]->epoch.gpsNanoSeconds = segStartTime.gpsNanoSeconds;
      memset(bankVeto[LAL_NUM_IFO]->data->data, 0,\
          params->numAnalPoints * sizeof(REAL4));
    }
    if (params->doSnglChiSquared)
    {
      for (ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if (params->haveTrig[ifoNumber])
        {
          bankVeto[ifoNumber]->epoch.gpsSeconds = segStartTime.gpsSeconds;
          bankVeto[ifoNumber]->epoch.gpsNanoSeconds = \
              segStartTime.gpsNanoSeconds;
          memset(bankVeto[ifoNumber]->data->data, 0,\
              params->numAnalPoints * sizeof(REAL4));
        }
      }
    }
  }

  if (params->doAutoVeto)
  {
    if (params->numIFO != 1)
    {
      autoVeto[LAL_NUM_IFO]->epoch.gpsSeconds = segStartTime.gpsSeconds;
      autoVeto[LAL_NUM_IFO]->epoch.gpsNanoSeconds = segStartTime.gpsNanoSeconds;
      memset(autoVeto[LAL_NUM_IFO]->data->data, 0,\
          params->numAnalPoints * sizeof(REAL4));
    }
    if (params->doSnglChiSquared)
    {
      for (ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if (params->haveTrig[ifoNumber])
        {
          autoVeto[ifoNumber]->epoch.gpsSeconds = segStartTime.gpsSeconds;
          autoVeto[ifoNumber]->epoch.gpsNanoSeconds = \
              segStartTime.gpsNanoSeconds;
          memset(autoVeto[ifoNumber]->data->data, 0,\
              params->numAnalPoints * sizeof(REAL4));
        }
      }
    }
  }

  if (params->doChiSquare)
  {
    if (params->numIFO != 1)
    {
      chiSquare[LAL_NUM_IFO]->epoch.gpsSeconds = segStartTime.gpsSeconds;
      chiSquare[LAL_NUM_IFO]->epoch.gpsNanoSeconds = \
          segStartTime.gpsNanoSeconds;
      memset(chiSquare[LAL_NUM_IFO]->data->data, 0,\
          params->numAnalPoints * sizeof(REAL4));
    }
    if (params->doSnglChiSquared)
    {
      for (ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if (params->haveTrig[ifoNumber])
        {
          chiSquare[ifoNumber]->epoch.gpsSeconds = segStartTime.gpsSeconds;
          chiSquare[ifoNumber]->epoch.gpsNanoSeconds = \
              segStartTime.gpsNanoSeconds;
          memset(chiSquare[ifoNumber]->data->data, 0,\
              params->numAnalPoints * sizeof(REAL4));
        }
      }
    }
  }

  /* Work out how many pValue arrays are needed */
  if (spinTemplates)
    vecLength = 5;
  else
    vecLength = 2;
  if (params->numIFO == 1 || params->singlePolFlag || params->faceOnStatistic)
    vecLengthTwo = vecLength;
  else
    vecLengthTwo = 2* vecLength;

  for (ui = 0 ; ui < vecLengthTwo ; ui++)
  {
    pValues[ui]->epoch.gpsSeconds = segStartTime.gpsSeconds;
    pValues[ui]->epoch.gpsNanoSeconds = segStartTime.gpsNanoSeconds;
    memset(pValues[ui]->data->data, 0,\
        params->numAnalPoints * sizeof(REAL4));
  }

  if (spinTemplates)
  {
    for (ui = 0 ; ui < 2 ; ui++)
    {
      gammaBeta[ui]->epoch.gpsSeconds = segStartTime.gpsSeconds;
      gammaBeta[ui]->epoch.gpsNanoSeconds = segStartTime.gpsNanoSeconds;
      memset(gammaBeta[ui]->data->data,0 ,\
          params->numAnalPoints * sizeof(REAL4));
    }
  }

}


void coh_PTF_destroy_time_series(
  REAL4TimeSeries          *cohSNR,
  REAL4TimeSeries          *nullSNR,
  REAL4TimeSeries          *traceSNR,
  REAL4TimeSeries          **bankVeto,
  REAL4TimeSeries          **autoVeto,
  REAL4TimeSeries          **chiSquare,
  REAL4TimeSeries          **pValues,
  REAL4TimeSeries          **gammaBeta,
  REAL4TimeSeries          **snrComps
)
{
  UINT4 k;

  if (cohSNR) XLALDestroyREAL4TimeSeries(cohSNR);
  if (nullSNR) XLALDestroyREAL4TimeSeries(nullSNR);
  if (traceSNR) XLALDestroyREAL4TimeSeries(traceSNR);

  for (k = 0; k < LAL_NUM_IFO+1; k++)
  {
    if (bankVeto[k])
    {
      XLALDestroyREAL4TimeSeries(bankVeto[k]);
      bankVeto[k]=NULL;
    }
    if (autoVeto[k])
    {
      XLALDestroyREAL4TimeSeries(autoVeto[k]);
      autoVeto[k]=NULL;
    }
    if (chiSquare[k])
    {
      XLALDestroyREAL4TimeSeries(chiSquare[k]);
      chiSquare[k]=NULL;
    }
  }


  for (k = 0 ; k < 10 ; k++)
  {
    if (pValues[k])
    {
        XLALDestroyREAL4TimeSeries(pValues[k]);
        pValues[k] = NULL;
    }
  }

  for (k = 0; k < LAL_NUM_IFO; k++)
  {
    if (snrComps[k])
    {
      XLALDestroyREAL4TimeSeries(snrComps[k]);
      snrComps[k] = NULL;
    }
  }
  for (k = 0; k < 2; k++)
  {
    if (gammaBeta[k])
    {
      XLALDestroyREAL4TimeSeries(gammaBeta[k]);
      gammaBeta[k] = NULL;
    }
  }
}

void coh_PTF_calculate_single_detector_filters(
  struct coh_PTF_params      *params,
  FindChirpTemplate          *fcTmplt,
  REAL4FrequencySeries       **invspec,
  REAL8Array                 **PTFM,
  COMPLEX8VectorSequence     **PTFqVec,
  REAL4TimeSeries            **snrComps,
  UINT4                      **snglAcceptPoints,
  UINT4                      *snglAcceptCount,
  RingDataSegments           **segments,
  COMPLEX8FFTPlan            *invPlan,
  UINT4                      spinTemplate,
  UINT4                      segNum
)
{
  UINT4 ifoNumber,ui,uj,localCount;
  REAL4 reSNRcomp,imSNRcomp;
  REAL4 *localSNRData;
  UINT4 *localAcceptPoints;

  REAL4 acceptThresh = params->snglSNRThreshold;

  for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if (params->haveTrig[ifoNumber])
    {
      localCount = 0;
      localSNRData = snrComps[ifoNumber]->data->data;
      localAcceptPoints = snglAcceptPoints[ifoNumber];
      /* Zero the storage vectors for the PTF filters */
      if (params->approximant == FindChirpPTF)
      {
        memset(PTFM[ifoNumber]->data, 0, 25 * sizeof(REAL8));
        memset(PTFqVec[ifoNumber]->data, 0,
              5 * params->numTimePoints * sizeof(COMPLEX8));
      }
      else
      {
        memset(PTFM[ifoNumber]->data, 0, 1 * sizeof(REAL8));
        memset(PTFqVec[ifoNumber]->data, 0,
              1 * params->numTimePoints * sizeof(COMPLEX8));
      }

      /* Here (h|s) and (h|h) are calculated */
      coh_PTF_normalize(params, fcTmplt, invspec[ifoNumber],
                        PTFM[ifoNumber], NULL, PTFqVec[ifoNumber],
                        &segments[ifoNumber]->sgmnt[segNum], invPlan,
                        spinTemplate);

      /* Here we calculate the single detector SNR */
      if (spinTemplate)
      {
        localCount = coh_PTF_calculate_single_det_spin_snr(params,PTFM,PTFqVec,\
                         snrComps,ifoNumber,localAcceptPoints);
      }
      else
      {
        for (ui = params->analStartPointBuf; ui < params->analEndPointBuf; ++ui)
        { /* Loop over time */
          uj = ui - params->analStartPointBuf;
          reSNRcomp = crealf(PTFqVec[ifoNumber]->data[ui]);
          imSNRcomp = cimagf(PTFqVec[ifoNumber]->data[ui]);
          localSNRData[uj] = sqrt((reSNRcomp*reSNRcomp + imSNRcomp*imSNRcomp)/
                   PTFM[ifoNumber]->data[0]);
          if (localSNRData[uj] > acceptThresh)
          {
            localAcceptPoints[localCount] = uj;
            localCount++;
          }
        }
      }
      snglAcceptCount[ifoNumber] = localCount;
      
    }
  }
}

UINT4 coh_PTF_calculate_single_det_spin_snr(
  struct coh_PTF_params      *params,
  REAL8Array                 **PTFM,
  COMPLEX8VectorSequence     **PTFqVec,
  REAL4TimeSeries            **snrComps,
  UINT4                      ifoNumber,
  UINT4                      *localAcceptPoints
)
{
  UINT4 ui,uj,uk,localCount;
  gsl_matrix *PTFmatrix;
  gsl_vector *eigenvalsSngl;
  gsl_matrix *eigenvecsSngl;
  gsl_eigen_symmv_workspace *matTemp = gsl_eigen_symmv_alloc(5);
  REAL4 v1_dot_u1, v1_dot_u2, v2_dot_u1, v2_dot_u2,max_eigen;
  REAL4 *snglv1p,*snglv2p;
  eigenvecsSngl = gsl_matrix_alloc(5,5);
  eigenvalsSngl = gsl_vector_alloc(5);
  snglv1p = LALCalloc(5 , sizeof(REAL4));
  snglv2p = LALCalloc(5 , sizeof(REAL4));
  REAL4 acceptThresh = params->snglSNRThreshold;

  /* convert PTFM to gsl_matrix */
  PTFmatrix = gsl_matrix_alloc(5,5);
  for (ui = 0; ui < 5; ui++ )
  {
    for (uj = 0; uj < 5; uj++ )
    {
      gsl_matrix_set(PTFmatrix, ui, uj, PTFM[ifoNumber]->data[ui*5+uj]);
    }
  }

  /* calculate eigenvectors and eigenvalues of (h|h)*/
  gsl_eigen_symmv(PTFmatrix, eigenvalsSngl, eigenvecsSngl,matTemp);
  gsl_eigen_symmv_free(matTemp);
  gsl_matrix_free(PTFmatrix);
  localCount = 0;
  for (ui = params->analStartPointBuf; ui < params->analEndPointBuf; ++ui)
  {  /* loop over time */
    uk = ui - params->analStartPointBuf;
    coh_PTF_calculate_rotated_vectors(params,PTFqVec,snglv1p,snglv2p,NULL,\
            NULL,NULL,eigenvecsSngl,eigenvalsSngl,\
            params->numTimePoints,ui,5,5,ifoNumber);

    /* Compute the dot products */
    v1_dot_u1 = v1_dot_u2 = v2_dot_u1 = v2_dot_u2 = max_eigen = 0.0;
    for (uj = 0; uj < 5; uj++)
    {
      v1_dot_u1 += snglv1p[uj] * snglv1p[uj];
      v1_dot_u2 += snglv1p[uj] * snglv2p[uj];
      v2_dot_u2 += snglv2p[uj] * snglv2p[uj];
    }
    max_eigen =0.5 * (v1_dot_u1 + v2_dot_u2 + sqrt((v1_dot_u1 - v2_dot_u2) * \
        (v1_dot_u1 - v2_dot_u2) + 4 * v1_dot_u2 * v1_dot_u2));
    snrComps[ifoNumber]->data->data[ui-params->analStartPointBuf] =\
        sqrt(max_eigen);
    if (snrComps[ifoNumber]->data->data[ui-params->analStartPointBuf] >\
         acceptThresh)
    {
      localAcceptPoints[localCount] = uk;
      localCount++;
    }

  }
  LALFree(snglv1p);
  LALFree(snglv2p);
  return localCount;
}

REAL4 coh_PTF_get_spin_SNR(
  REAL4 *v1p,
  REAL4 *v2p,
  UINT4 vecLengthTwo
)
{
  UINT4 ui;
  REAL4 v1_dot_u1,v1_dot_u2,v2_dot_u2,max_eigen;
  v1_dot_u1 = v1_dot_u2 = v2_dot_u2 = 0.0;
  for (ui =0; ui < vecLengthTwo; ui++)
  {
    v1_dot_u1 += v1p[ui] * v1p[ui];
    v2_dot_u2 += v2p[ui] * v2p[ui];
    v1_dot_u2 += v1p[ui] * v2p[ui];
  }
  max_eigen = 0.5 * (v1_dot_u1 + v2_dot_u2 + sqrt((v1_dot_u1 - v2_dot_u2)
            * (v1_dot_u1 - v2_dot_u2) + 4 * v1_dot_u2 * v1_dot_u2));
  return sqrt(max_eigen);
}

/* THIS FUNCTION IS COMMENTED OUT SO IT CAN BE VERIFIED AND FIXED*/
/* void coh_PTF_get_spin_amp_terms(
XXXXX
)
  coh_PTF_calculate_rotated_vectors(params,PTFqVec,v1p,v2p,Fplus,Fcross,timeOffsetPoints,
        eigenvecs,eigenvals,numPoints,i,vecLength,vecLengthTwo,LAL_NUM_IFO);
  v1_dot_u1 = v1_dot_u2 = v2_dot_u1 = v2_dot_u2 = 0;
  for (j = 0; j < vecLengthTwo; j++)
  {
    u1[j] = v1p[j] / (pow(gsl_vector_get(eigenvals,j),0.5));
    u2[j] = v2p[j] / (pow(gsl_vector_get(eigenvals,j),0.5));
    v1[j] = u1[j] * gsl_vector_get(eigenvals,j);
    v2[j] = u2[j] * gsl_vector_get(eigenvals,j);
    v1_dot_u1 += v1[j]*u1[j];
    v1_dot_u2 += v1[j]*u2[j];
    v2_dot_u2 += v2[j]*u2[j];
  }
  dCee = (max_eigen - v1_dot_u1) / v1_dot_u2;
  dAlpha = 1./(v1_dot_u1 + dCee * 2 * v1_dot_u2 + dCee*dCee*v2_dot_u2);
  dAlpha = pow(dAlpha,0.5);
  dBeta = dCee*dAlpha;
  // The p Values are calculated in the rotated frame
  for (j = 0 ; j < vecLengthTwo ; j++)
  {
    pValsTemp[j] = dAlpha*u1[j] + dBeta*u2[j];
    pValues[j]->data->data[i - numPoints/4] = 0.;
  }
  // This loop can be used to verify that the SNR obtained is as before
  recSNR = 0;
  for (j = 0 ; j < vecLengthTwo ; j++)
  {
    for (k = 0 ; k < vecLengthTwo ; k++)
    {
      recSNR += pValsTemp[j]*pValsTemp[k] * (v1[j]*v1[k]+v2[j]*v2[k]);
    }
  }
  // Then we calculate the two phase/amplitude terms beta and gamma
  // These are explained in Diego's thesis
  betaGammaTemp[0] = 0;
  betaGammaTemp[1] = 0;
  for (j = 0 ; j < vecLengthTwo ; j++)
  {
    betaGammaTemp[0] += pValsTemp[j]*v1[j];
    betaGammaTemp[1] += pValsTemp[j]*v2[j];
  }
  gammaBeta[0]->data->data[i - numPoints/4] = betaGammaTemp[0];
  gammaBeta[1]->data->data[i - numPoints/4] = betaGammaTemp[1];

  // The p Values need to be rotated back into the original frame.
  // Currently we are recording values in rotated frame
  for (j = 0 ; j < vecLengthTwo ; j++)
  {
    for (k = 0 ; k < vecLengthTwo ; k++)
    {
      pValues[j]->data->data[i-numPoints/4]+=gsl_matrix_get(eigenvecs,j,k)*pValsTemp[k];
    }
  }

  // And we check that this still gives the expected SNR in the
  // unrotated basis.
  for (j = 0; j < vecLengthTwo ; j++) // Construct the vi vectors
  {
    v1[j] = 0.;
    v2[j] = 0.;
    for(k = 0; k < LAL_NUM_IFO; k++)
    {
      if (params->haveTrig[k])
      {
        if (j < vecLength)
        {
          v1[j] += Fplus[k] * crealf(PTFqVec[k]->data[j*numPoints+i+timeOffsetPoints[k]]);
          v2[j] += Fplus[k] * cimagf(PTFqVec[k]->data[j*numPoints+i+timeOffsetPoints[k]]);
        }
        else
        {
          v1[j] += Fcross[k] * crealf(PTFqVec[k]->data[(j-vecLength)*numPoints+i+timeOffsetPoints[k]]);
          v2[j] += Fcross[k] * cimagf(PTFqVec[k]->data[(j-vecLength)*numPoints+i+timeOffsetPoints[k]]);
        }
      }
    }
  }
  recSNR = 0;
  for (j = 0 ; j < vecLengthTwo ; j++)
  {
    for (k = 0 ; k < vecLengthTwo ; k++)
    {
      recSNR += pValues[j]->data->data[i-numPoints/4]*pValues[k]->data->data[i-numPoints/4] * (v1[j]*v1[k]+v2[j]*v2[k]);
    }
  }
}
*/

void coh_PTF_calculate_coherent_SNR(
  struct coh_PTF_params      *params,
  REAL4                      *snrData,
  REAL4TimeSeries            **pValues,
  REAL4TimeSeries            **snrComps,
  INT4                       *timeOffsetPoints,
  COMPLEX8VectorSequence     **PTFqVec,
  REAL4                      *Fplus,
  REAL4                      *Fcross,
  gsl_matrix                 *eigenvecs,
  gsl_vector                 *eigenvals,
  UINT4                      segStartPoint,
  UINT4                      segEndPoint,
  UINT4                      vecLength,
  UINT4                      vecLengthTwo,
  UINT4                      spinTemplate,
  UINT4                      **snglAcceptPoints,
  UINT4                      *snglAcceptCount
)
{
  REAL4 snglSNRthresh = params->snglSNRThreshold;
  REAL4 *v1p = LALCalloc(vecLengthTwo , sizeof(REAL4));
  REAL4 *v2p = LALCalloc(vecLengthTwo , sizeof(REAL4));

  UINT4 i,j,k,ifoNumber,ifoNumber2,currPointLoc,ifoNum1,ifoNum2;
  UINT4 localCount,*localAcceptPoints,localOffset;
  INT4 tOffset1,tOffset2;
  REAL4 v1_dot_u1,v2_dot_u2,max_eigen,coincSNR;
  REAL4 cohSNRThresholdSq = params->threshold * params->threshold;

  /* If only two detectors & standard analysis identify the 2 detectors
   * up front for speed
   */
  ifoNum1 = ifoNum2 = tOffset1 = tOffset2 = 0;
  if (params->numIFO == 2 && (! params->singlePolFlag) &&\
      (!params->faceOnStatistic) )
  {
    for (ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
    {
      if (params->haveTrig[ifoNumber])
      {
        if (ifoNum1 == 0)
          ifoNum1 = ifoNumber;
        else if (ifoNum2 == 0)
          ifoNum2 = ifoNumber;
      }
    }
    tOffset1 = timeOffsetPoints[ifoNum1] - params->analStartPointBuf;
    tOffset2 = timeOffsetPoints[ifoNum2] - params->analStartPointBuf;
  }

  for (ifoNumber2 = 0; ifoNumber2 < LAL_NUM_IFO; ifoNumber2++)
  {
    if (! params->haveTrig[ifoNumber2])
    {
      continue;
    }
    localCount = snglAcceptCount[ifoNumber2];
    localAcceptPoints = snglAcceptPoints[ifoNumber2];
    localOffset = params->analStartPointBuf-timeOffsetPoints[ifoNumber2];
    for (k = 0; k < localCount; ++k)
    {
      i = localAcceptPoints[k] + localOffset;
      currPointLoc = i-params->analStartPoint;
      /* Continue if not in this analysis segment */
      if ((i < segStartPoint) || (i > segEndPoint))
      {
        continue;
      }

      /* Don't bother calculating coherent SNR if all ifo's SNR is less than
         some value */
      for (ifoNumber = 0;ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if (params->haveTrig[ifoNumber])
        {
          if (snrComps[ifoNumber]->data->
                data[i-params->analStartPointBuf+timeOffsetPoints[ifoNumber]]
                > snglSNRthresh)
          {
            break;
          }
        }
      }
      if (ifoNumber == LAL_NUM_IFO)
      {
        snrData[currPointLoc] = 0;
        fprintf(stderr,"I shouldn't be here now!!\n");
        continue;
      }

      if (params->numIFO == 2 && (! params->singlePolFlag) &&\
          (!params->faceOnStatistic) )
      { /*If only 2 detectors cohSNR = coincident SNR. SO just use that */
        max_eigen = snrComps[ifoNum1]->data->data[i+tOffset1] *
                            snrComps[ifoNum1]->data->data[i+tOffset1] +
                            snrComps[ifoNum2]->data->data[i+tOffset2] *
                            snrComps[ifoNum2]->data->data[i+tOffset2];
        if (max_eigen < cohSNRThresholdSq)
        {
          snrData[currPointLoc] = 0;
          continue;
        }
        snrData[currPointLoc] = sqrt(max_eigen);
        if (params->storeAmpParams)
        {
         coh_PTF_calculate_rotated_vectors(params,PTFqVec,v1p,v2p,Fplus,Fcross,
            timeOffsetPoints,eigenvecs,eigenvals,
            params->numTimePoints,i,vecLength,vecLengthTwo,LAL_NUM_IFO);
          for (j = 0 ; j < vecLengthTwo ; j++)
          {
            pValues[j]->data->data[currPointLoc] = \
                           v1p[j] / (pow(gsl_vector_get(eigenvals,j),0.5));
            pValues[j+vecLengthTwo]->data->data[currPointLoc] = \
                           v2p[j] / (pow(gsl_vector_get(eigenvals,j),0.5));
          }
        }
      }
      else
      {
        coincSNR = 0;
        /* Calculate the coincident SNR at this time point*/
        for (ifoNumber = 0;ifoNumber < LAL_NUM_IFO; ifoNumber++)
        {
          if (params->haveTrig[ifoNumber])
          {
            coincSNR += snrComps[ifoNumber]->data->data[\
                i - params->analStartPointBuf + timeOffsetPoints[ifoNumber]]
                * snrComps[ifoNumber]->data->data[\
                i - params->analStartPointBuf + timeOffsetPoints[ifoNumber]];
          }
        }
        /* Do not need to calculate coherent SNR if coinc SNR < threshold */
        /* NOTE: Cheaper to compare coincSNRSq than use pow(x,0.5) */
        if (coincSNR < cohSNRThresholdSq)
        {
          snrData[currPointLoc] = 0;
          continue;
        }

        /* This function combines the various (Q_i | s) and rotates them into
         * the orthonormal basis, using the eigen[vector,value]s.
         */
        coh_PTF_calculate_rotated_vectors(params,PTFqVec,v1p,v2p,Fplus,Fcross,
          timeOffsetPoints,eigenvecs,eigenvals,
          params->numTimePoints,i,vecLength,vecLengthTwo,LAL_NUM_IFO);

        /* And SNR is calculated
         * For non-spin+multi-site+coherent: 
         *               v1p[0] * v1p[0] = (\bf{F}_+\bf{h}_0 | \bf{s})^2
         *               v1p[1] * v1p[1] = (\bf{F}_x\bf{h}_0 | \bf{s})^2
         *               v2p[0] * v2p[0] = (\bf{F}_+\bf{h}_{\pi/2} | \bf{s})^2
         *               v2p[1] * v2p[1] = (\bf{F}_x\bf{h}_{\pi/2} | \bf{s})^2
         * For non-spin+single-site/face-on there will only be two components
         * in this calculation, but otherwise the same.
         */
        if (spinTemplate == 0)
        {
          v1_dot_u1 = v2_dot_u2 = 0.0;
          for (j = 0; j < vecLengthTwo; j++)
          {
            v1_dot_u1 += v1p[j] * v1p[j];
            v2_dot_u2 += v2p[j] * v2p[j];
          }
          max_eigen = (v1_dot_u1 + v2_dot_u2);
          if (max_eigen < cohSNRThresholdSq)
          {
            snrData[currPointLoc] = 0;
            continue;
          }
          else
          {
            snrData[currPointLoc] = sqrt(max_eigen);
            if (params->storeAmpParams)
            {
              for (j = 0 ; j < vecLengthTwo ; j++)
              {
                pValues[j]->data->data[currPointLoc] = \
                               v1p[j] / (pow(gsl_vector_get(eigenvals,j),0.5));
                pValues[j+vecLengthTwo]->data->data[currPointLoc] = \
                               v2p[j] / (pow(gsl_vector_get(eigenvals,j),0.5));
              }
            }
          }
        }
        else
        { /* Spinning case follow PTF notation to get SNR */
          snrData[i-params->analStartPoint] = \
              coh_PTF_get_spin_SNR(v1p,v2p,vecLengthTwo);
          if (params->storeAmpParams)
          {
            fprintf(stderr,"Spinning amplitude stuff is currently disabled\n");
            /* coh_PTF_get_spin_amp_terms(.......) */
          }
        }
      }
    }
  }
}

UINT4 coh_PTF_template_time_series_cluster(
  struct coh_PTF_params      *params,
  REAL4TimeSeries *cohSNR,
  UINT4 *acceptPoints,
  INT4 *timeOffsetPoints,
  INT4 numPointCheck,
  UINT4 startPoint,
  UINT4 endPoint,
  UINT4 **snglAcceptPoints,
  UINT4 *snglAcceptCount
)
{
  UINT4 ui,ifoNumber2,k,i;
  UINT4 count = 0;
  UINT4 logicArray[cohSNR->data->length];
  INT4 j,tempPoint;
  /* FIXME: Bizarrely, the mac clang compiler will "optimize away" this variable
   * if I do not use the volatile keyword, which causes seg faults as it should
   * not be "optimized away". Is this a bug in clang??? */
  volatile INT4 localOffset;
  UINT4 localCount,*localAcceptPoints;
  /* Have to cast from UINT4 to INT4 to avoid warning */
  INT4 startPointI = (INT4) startPoint;
  INT4 endPointI = (INT4) endPoint;

  REAL4 *snrData = cohSNR->data->data;

  for (ifoNumber2 = 0; ifoNumber2 < LAL_NUM_IFO; ifoNumber2++)
  {
    if (! params->haveTrig[ifoNumber2])
    {
      continue;
    }
    localCount = snglAcceptCount[ifoNumber2];
    localAcceptPoints = snglAcceptPoints[ifoNumber2];
    localOffset = params->analStartPointBuf-timeOffsetPoints[ifoNumber2];
    for (k = 0; k < localCount; ++k)
    {
      i = localAcceptPoints[k] + localOffset;
      ui = i-params->analStartPoint;
      /* Continue if not in this analysis segment */
      if ((ui < startPoint) || (ui > endPoint))
      {
        continue;
      }
      logicArray[ui] = 1;
      if (snrData[ui])
      {
        for (j = -numPointCheck; j < numPointCheck; j++)
        {
          tempPoint = ui + j;
          if (tempPoint < startPointI)
          {
            continue;
          }
          if (tempPoint >= endPointI)
          {
            continue;
          }
          if (snrData[tempPoint] > snrData[ui])
          {
            logicArray[ui] = 0;
            break;
          }
        }
      }
    }
  }

  for (ifoNumber2 = 0; ifoNumber2 < LAL_NUM_IFO; ifoNumber2++)
  {
    if (! params->haveTrig[ifoNumber2])
    {
      continue;
    }
    localCount = snglAcceptCount[ifoNumber2];
    localAcceptPoints = snglAcceptPoints[ifoNumber2];
    localOffset = params->analStartPointBuf-timeOffsetPoints[ifoNumber2];
    for (k = 0; k < localCount; ++k)
    {
      i = localAcceptPoints[k] + localOffset;
      ui = i-params->analStartPoint;
      /* Continue if not in this analysis segment */
      if ((ui < startPoint) || (ui > endPoint))
      {
        continue;
      }

      if (! logicArray[ui])
      {
        snrData[ui] = 0.;
      }
      else if (logicArray[ui] == 1)
      {
        if (snrData[ui])
        {
          acceptPoints[count] = ui;
          count++;  
          logicArray[ui] = 2;
        }
      }
    }
  }
  return count;
}

UINT4 coh_PTF_test_veto_vals(
  struct coh_PTF_params    *params,
  REAL4TimeSeries          *cohSNR,
  REAL4TimeSeries          *nullSNR,
  REAL4TimeSeries          **bankVeto,
  REAL4TimeSeries          **autoVeto,
  UINT4                    currPointLoc
)
{
  /* This function checks whether chisq needs to be calculated. It does this
   * by checking if the trigger will fail the bank/auto/null vetoes */ 

  /* NOTE: The code currently uses the null stream SNR. It should instead/also
   * use the nullSNR calculated from coincSNR - cohSNR */ 

  UINT4 numDOF;
  REAL4 bestNR,currSNR;

  numDOF = 4.;
  if (params->singlePolFlag || params->faceOnStatistic || params->numIFO == 1)
    numDOF = 2.;

  currSNR = cohSNR->data->data[currPointLoc];

  /* IS the null stream too large? */
  if (params->doNullStream)
  {
    if (nullSNR->data->data[currPointLoc] > params->nullStatThreshold \
        && currSNR < params->nullStatGradOn)
    {
      return 1;
    }
    if (currSNR > params->nullStatGradOn)
    {
      if (nullSNR->data->data[currPointLoc] > (params->nullStatThreshold \
           + (currSNR - params->nullStatGradOn)*params->nullStatGradient))
      {
        return 1;
      }
    }
  }

  /* Is bank new SNR too large? */
  if (params->doBankVeto)
  {
    if (bankVeto[LAL_NUM_IFO]->data->data[currPointLoc] > \
        params->BVsubBankSize*numDOF)
    {
      bestNR = currSNR / pow( ( 1. + pow( \
          bankVeto[LAL_NUM_IFO]->data->data[currPointLoc] / \
          ((REAL4)params->BVsubBankSize*numDOF ) , \
          params->bankVetoq/params->bankVeton) )/2. , 1./params->bankVetoq);
      if (bestNR < params->chiSquareCalcThreshold)
      {
        return 1;
      }
    }
  }

  /* Is auto new SNR too large */
  if (params->doAutoVeto)
  {
    if (autoVeto[LAL_NUM_IFO]->data->data[currPointLoc] > \
        params->numAutoPoints*numDOF)
    {
      bestNR = currSNR / pow( ( 1. + pow( \
          autoVeto[LAL_NUM_IFO]->data->data[currPointLoc] / \
          ((REAL4)params->numAutoPoints*numDOF ) , \
          params->autoVetoq/params->autoVeton) )/2. , 1./params->autoVetoq);
      if (bestNR < params->chiSquareCalcThreshold)
      {
        return 1;
      }
    }
  }
  return 0;
}

void coh_PTF_calculate_null_stream_filters(
  struct coh_PTF_params      *params,
  FindChirpTemplate          *fcTmplt,
  REAL4FrequencySeries       **invspec,
  REAL8Array                 **PTFM,
  COMPLEX8VectorSequence     **PTFqVec,
  RingDataSegments           **segments,
  COMPLEX8FFTPlan            *invPlan,
  UINT4                      spinTemplate,
  UINT4                      segNum
)
{
  if (params-> approximant == FindChirpPTF)
  {
    memset(PTFM[LAL_NUM_IFO]->data, 0, 25*sizeof(REAL8));
    memset(PTFqVec[LAL_NUM_IFO]->data, 0,\
            5*params->numTimePoints*sizeof(COMPLEX8));
  }
  else
  {
    memset(PTFM[LAL_NUM_IFO]->data, 0, 1*sizeof(REAL8));
    memset(PTFqVec[LAL_NUM_IFO]->data, 0,\
            1*params->numTimePoints*sizeof(COMPLEX8));
  }
  coh_PTF_normalize(params, fcTmplt, invspec[LAL_NUM_IFO],
                    PTFM[LAL_NUM_IFO], NULL, PTFqVec[LAL_NUM_IFO],
                    &segments[LAL_NUM_IFO]->sgmnt[segNum], invPlan,
                    spinTemplate);
}

void coh_PTF_calculate_null_stream_norms(
  UINT4 vecLength,
  gsl_matrix *eigenvecsNull,
  gsl_vector *eigenvalsNull,
  REAL8Array *PTFM[LAL_NUM_IFO+1]
)
{
  /* This function calculates the eigenvectors/values needed to orthonormalize
   * the filters when using the null stream
   */

  /* For non-spinning case this is trivial */
  if (vecLength == 1)
  {
    gsl_matrix_set(eigenvecsNull, 0 , 0, 1);
    gsl_vector_set(eigenvalsNull, 0, PTFM[LAL_NUM_IFO]->data[0]);
  }
  /* For spinning case there is more to this */
  else
  {
    UINT4 i,j;
    gsl_eigen_symmv_workspace *matTempNull = gsl_eigen_symmv_alloc(vecLength);
    gsl_matrix *B2Null = gsl_matrix_alloc(vecLength, vecLength);
    for (i = 0; i < vecLength; i++)
    {
      for (j = 0; j < vecLength; j++)
      {
        gsl_matrix_set(B2Null, i, j, PTFM[LAL_NUM_IFO]->data[i*5+j]);
      }
    }
    gsl_eigen_symmv(B2Null, eigenvalsNull, eigenvecsNull, matTempNull);
  }
}

void coh_PTF_calculate_null_stream_snr(
  struct coh_PTF_params   *params,
  REAL4TimeSeries         *nullSNR,
  COMPLEX8VectorSequence  **PTFqVec,
  gsl_matrix              *eigenvecsNull,
  gsl_vector              *eigenvalsNull,
  UINT4                   spinTemplate,
  UINT4                   vecLength,
  UINT4                   vecLoc,
  UINT4                   snrLoc
)
{
  UINT4 j,k;
  REAL4 v1_dot_u1,v1_dot_u2,v2_dot_u2,max_eigen;
  REAL4 v1N[vecLength],v2N[vecLength],u1N[vecLength],u2N[vecLength];
  /* Begin by rotating to the preferred vector */
  /* NOTE: For non-spin vecLength=1 and some of this is trivial */
  /* NOTE: This code could be optimized, but is rarely used so has not been */
  for (j = 0; j < vecLength; j++)
  {
    v1N[j] = crealf(PTFqVec[LAL_NUM_IFO]->data[j*params->numTimePoints+vecLoc]);
    v2N[j] = cimagf(PTFqVec[LAL_NUM_IFO]->data[j*params->numTimePoints+vecLoc]);
  }

  for (j = 0 ; j < vecLength ; j++)
  {
    u1N[j] = 0.;
    u2N[j] = 0.;
    for (k = 0 ; k < vecLength ; k++)
    {
      u1N[j] += gsl_matrix_get(eigenvecsNull,k,j)*v1N[k];
      u2N[j] += gsl_matrix_get(eigenvecsNull,k,j)*v2N[k];
    }
    u1N[j] = u1N[j] / (pow(gsl_vector_get(eigenvalsNull,j),0.5));
    u2N[j] = u2N[j] / (pow(gsl_vector_get(eigenvalsNull,j),0.5));
  }
  /* Compute the dot products */
  v1_dot_u1 = v1_dot_u2 = v2_dot_u2 = 0.0;
  for (j = 0; j < vecLength; j++)
  {
    v1_dot_u1 += u1N[j] * u1N[j];
    v1_dot_u2 += u1N[j] * u2N[j];
    v2_dot_u2 += u2N[j] * u2N[j];
  }
  if (spinTemplate == 0)
  {
    max_eigen = 0.5 * (v1_dot_u1 + v2_dot_u2);
  }
  else
  {
    max_eigen = 0.5*(v1_dot_u1+v2_dot_u2+sqrt((v1_dot_u1-v2_dot_u2)
        * (v1_dot_u1 - v2_dot_u2) + 4 * v1_dot_u2 * v1_dot_u2));
  }
  nullSNR->data->data[snrLoc] = sqrt(max_eigen);
}      

void coh_PTF_calculate_trace_snr(
  struct coh_PTF_params   *params,
  REAL4TimeSeries         *traceSNR,
  COMPLEX8VectorSequence  **PTFqVec,
  gsl_matrix              *eigenvecs,
  gsl_vector              *eigenvals,
  REAL4                   *Fplus,
  REAL4                   *Fcross,
  INT4                    *timeOffsetPoints,
  UINT4                   spinTemplate,
  UINT4                   vecLength,
  UINT4                   vecLengthTwo,
  UINT4                   vecLoc,
  UINT4                   snrLoc
)
{
  UINT4 j,k,m;
  REAL4 v1_dot_u1,v1_dot_u2,v2_dot_u2,max_eigen,traceSNRsq;
  REAL4 v1[vecLength],v2[vecLength],u1[vecLength],u2[vecLength];

  /* Trace SNR is the coherent SNR with no cross-detector terms */
  traceSNRsq = 0;
  for(k = 0; k < LAL_NUM_IFO; k++)
  {
    if (params->haveTrig[k])
    {
      for (j = 0; j < vecLengthTwo ; j++)
      {
        if (j < vecLength)
        {
          v1[j] = Fplus[k] * crealf(PTFqVec[k]->data[\
                         j*params->numTimePoints+vecLoc+timeOffsetPoints[k]]);
          v2[j] = Fplus[k] * cimagf(PTFqVec[k]->data[\
                         j*params->numTimePoints+vecLoc+timeOffsetPoints[k]]);
        }
        else
        {
          v1[j] = Fcross[k] * crealf(PTFqVec[k]->data[ (j-vecLength) * \
                           params->numTimePoints+vecLoc+timeOffsetPoints[k]]);
          v2[j] = Fcross[k] * cimagf(PTFqVec[k]->data[ (j-vecLength) * \
                           params->numTimePoints+vecLoc+timeOffsetPoints[k]]);
        }
      }
      for (j = 0 ; j < vecLengthTwo ; j++)
      {
        u1[j] = 0.;
        u2[j] = 0.;
        for (m = 0 ; m < vecLengthTwo ; m++)
        {
          u1[j] += gsl_matrix_get(eigenvecs,m,j)*v1[m];
          u2[j] += gsl_matrix_get(eigenvecs,m,j)*v2[m];
        }
        u1[j] = u1[j] / (pow(gsl_vector_get(eigenvals,j),0.5));
        u2[j] = u2[j] / (pow(gsl_vector_get(eigenvals,j),0.5));
      }
      /* Compute the dot products */
      v1_dot_u1 = v1_dot_u2 = v2_dot_u2 = max_eigen = 0.0;
      for (j = 0; j < vecLengthTwo; j++)
      {
        v1_dot_u1 += u1[j] * u1[j];
        v1_dot_u2 += u1[j] * u2[j];
        v2_dot_u2 += u2[j] * u2[j];
      }
      if (spinTemplate == 0)
      {
        max_eigen = (v1_dot_u1 + v2_dot_u2);
      }
      else
      {
        max_eigen = 0.5 * (v1_dot_u1 + v2_dot_u2 + sqrt((v1_dot_u1 - v2_dot_u2)
            * (v1_dot_u1 - v2_dot_u2) + 4 * v1_dot_u2 * v1_dot_u2));
      }
      traceSNRsq += max_eigen;
    }
  }
  traceSNR->data->data[snrLoc] = sqrt(traceSNRsq);
}

void coh_PTF_convert_time_offsets_to_points(
  struct coh_PTF_params   *params,
  REAL4                   *timeOffsets,
  INT4                    *timeOffsetPoints
)
{
  UINT4 i;
  REAL4 deltaT =  (1.0/params->sampleRate);
  for (i = 0; i < LAL_NUM_IFO; i++)
  {
    timeOffsetPoints[i] = (int) floor(timeOffsets[i]/deltaT + 0.5);
  }
  /* The following is useful for debugging */

  /*
  verbose("Time offsets (s): ");
  for (i = 0; i < LAL_NUM_IFO; i++)
  {
    if (params->haveTrig[i])
    {
      verbose("%f ", timeOffsets[i]);
    }
  }
  verbose("\n");
  verbose("Time offsets (points): ");
  for (i = 0; i < LAL_NUM_IFO; i++)
  {
    if (params->haveTrig[i])
    {
      verbose("%d " , timeOffsetPoints[i]);
    }
  }
  verbose("\n");
  */
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
    COMPLEX8VectorSequence  **PTFqVec,
    REAL4 *u1,
    REAL4 *u2,
    REAL4 *Fplus,
    REAL4 *Fcross,
    INT4  *timeOffsetPoints,
    gsl_matrix *eigenvecs,
    gsl_vector *eigenvals,
    UINT4 numPoints,
    UINT4 position,
    UINT4 vecLength,
    UINT4 vecLengthTwo,
    UINT4 detectorNum)
{
  /* IMPORTANT!!! This function is probably the dominant computational cost
   * when running in lots-of-sky-points mode. Optimizing this is important! */ 

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
            if (params->faceOnStatistic == 1)
            {
              v1[j] += Fplus[k] * crealf(PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]]);
              v1[j] += Fcross[k] * cimagf(PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]]);
              v2[j] += Fcross[k] * crealf(PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]]);
              v2[j] -= Fplus[k] * cimagf(PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]]);
            }
            else if (params->faceOnStatistic == 2)
            {
              v1[j] += Fplus[k] * crealf(PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]]);
              v1[j] -= Fcross[k] * cimagf(PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]]);
              v2[j] += Fcross[k] * crealf(PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]]);
              v2[j] += Fplus[k] * cimagf(PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]]);
            }
            else
            {
              fprintf(stderr,"Face-on stat is not working!");
            }
          }
          else if (j < vecLength)
          {
            v1[j] += Fplus[k] * crealf(PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]]);
            v2[j] += Fplus[k] * cimagf(PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]]);
          }
          else
          {
            v1[j] += Fcross[k] * crealf(PTFqVec[k]->data[(j-vecLength)*numPoints+position+timeOffsetPoints[k]]);
            v2[j] += Fcross[k] * cimagf(PTFqVec[k]->data[(j-vecLength)*numPoints+position+timeOffsetPoints[k]]);
          }
        }
      }
    }
    else
    {
      v1[j] += crealf(PTFqVec[detectorNum]->data[j*numPoints+position]);
      v2[j] += cimagf(PTFqVec[detectorNum]->data[j*numPoints+position]);
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

MultiInspiralTable* coh_PTF_create_multi_event(
    struct coh_PTF_params   *params,
    REAL4TimeSeries         *cohSNR,
    FindChirpTemplate       *fcTmplt,
    InspiralTemplate        PTFTemplate,
    UINT8                   *eventId,
    UINT4                   spinTrigger,
    REAL4TimeSeries         **pValues,
    REAL4TimeSeries         **gammaBeta,
    REAL4TimeSeries         **snrComps,
    REAL4TimeSeries         *nullSNR,
    REAL4TimeSeries         *traceSNR,
    REAL4TimeSeries         **bankVeto,
    REAL4TimeSeries         **autoVeto,
    REAL4TimeSeries         **chiSquare,
    REAL8Array              **PTFM,
    REAL4                   rightAscension,
    REAL4                   declination,
    INT8                    slideId,
    INT4                   *timeOffsetPoints,
    UINT4                   currPos
)
{
  LIGOTimeGPS trigTime;
  UINT4 numDOF = 4;
  if ( params->singlePolFlag || params->faceOnStatistic )
    numDOF = 2;

  MultiInspiralTable *currEvent;
  currEvent = (MultiInspiralTable *)
      LALCalloc(1, sizeof(MultiInspiralTable));
  currEvent->event_id = (EventIDColumn *)
      LALCalloc(1, sizeof(EventIDColumn));
  currEvent->event_id->id=*eventId;
  currEvent->time_slide_id = (EventIDColumn *)
      LALCalloc(1, sizeof(EventIDColumn));
  currEvent->time_slide_id->id=slideId;
  (*eventId)++;
  trigTime = cohSNR->epoch;
  XLALGPSAdd(&trigTime,currPos*cohSNR->deltaT);
  currEvent->snr = cohSNR->data->data[currPos];
  currEvent->mass1 = PTFTemplate.mass1;
  currEvent->mass2 = PTFTemplate.mass2;
  currEvent->chi = PTFTemplate.chi;
  currEvent->kappa = PTFTemplate.kappa;
  currEvent->mchirp = PTFTemplate.totalMass*pow(PTFTemplate.eta,3.0/5.0);
  currEvent->eta = PTFTemplate.eta;
  /* FIXME: Add spins to MultiInspiralTable */
  currEvent->chi = PTFTemplate.spin1[2];
  currEvent->kappa = PTFTemplate.spin2[2];
  currEvent->end_time = trigTime;
  /* add sky position, but need to track back to sky fixed sky position */
  currEvent->ra = rightAscension -
                  XLALGreenwichMeanSiderealTime(&params->trigTime) +
                  XLALGreenwichMeanSiderealTime(&currEvent->end_time);
  currEvent->dec = declination;
  if (params->doNullStream)
    currEvent->null_statistic = nullSNR->data->data[currPos];
  if (params->doTraceSNR)
    currEvent->trace_snr = traceSNR->data->data[currPos];
  if (params->doBankVeto)
  {
    if (params->numIFO != 1)
    {
      currEvent->bank_chisq = bankVeto[LAL_NUM_IFO]->data->data[currPos];
      currEvent->bank_chisq_dof = numDOF * params->BVsubBankSize;
    }
    if (params->doSnglChiSquared)
    {
      if (bankVeto[LAL_IFO_G1])
      {
        currEvent->bank_chisq_g = bankVeto[LAL_IFO_G1]->data->data[currPos];
      }
      if (bankVeto[LAL_IFO_H1])
      {
        currEvent->bank_chisq_h1 = bankVeto[LAL_IFO_H1]->data->data[currPos];
      }
      if (bankVeto[LAL_IFO_H2])
      {
        currEvent->bank_chisq_h2 = bankVeto[LAL_IFO_H2]->data->data[currPos];
      }
      if (bankVeto[LAL_IFO_L1])
      {
        currEvent->bank_chisq_l = bankVeto[LAL_IFO_L1]->data->data[currPos];
      }
      if (bankVeto[LAL_IFO_T1])
      {
        currEvent->bank_chisq_t = bankVeto[LAL_IFO_T1]->data->data[currPos];
      }
      if (bankVeto[LAL_IFO_V1])
      {
        currEvent->bank_chisq_v = bankVeto[LAL_IFO_V1]->data->data[currPos];
      }
    }
  }
  if (params->doAutoVeto)
  {
    if (params->numIFO != 1)
    {
      currEvent->cont_chisq = autoVeto[LAL_NUM_IFO]->data->data[currPos];
      currEvent->cont_chisq_dof = numDOF * params->numAutoPoints;
    }
    if (params->doSnglChiSquared)
    {
      if (autoVeto[LAL_IFO_G1])
      {
        currEvent->cont_chisq_g = autoVeto[LAL_IFO_G1]->data->data[currPos];
      }
      if (autoVeto[LAL_IFO_H1])
      {
        currEvent->cont_chisq_h1 = autoVeto[LAL_IFO_H1]->data->data[currPos];
      }
      if (autoVeto[LAL_IFO_H2])
      {
        currEvent->cont_chisq_h2 = autoVeto[LAL_IFO_H2]->data->data[currPos];
      }
      if (autoVeto[LAL_IFO_L1])
      {
        currEvent->cont_chisq_l = autoVeto[LAL_IFO_L1]->data->data[currPos];
      }
      if (autoVeto[LAL_IFO_T1])
      {
        currEvent->cont_chisq_t = autoVeto[LAL_IFO_T1]->data->data[currPos];
      }
      if (autoVeto[LAL_IFO_V1])
      {
        currEvent->cont_chisq_v = autoVeto[LAL_IFO_V1]->data->data[currPos];
      }
    }
  }
  if (params->doChiSquare)
  {
    if (params->numIFO != 1)
    {
      currEvent->chisq = chiSquare[LAL_NUM_IFO]->data->data[currPos];
      currEvent->chisq_dof = numDOF * (params->numChiSquareBins - 1);
    }
    if (params->doSnglChiSquared)
    {
      if (chiSquare[LAL_IFO_G1])
      {
        currEvent->chisq_g = chiSquare[LAL_IFO_G1]->data->data[currPos];
      }
      if (chiSquare[LAL_IFO_H1])
      {
        currEvent->chisq_h1 = chiSquare[LAL_IFO_H1]->data->data[currPos];
      }
      if (chiSquare[LAL_IFO_H2])
      {
        currEvent->chisq_h2 = chiSquare[LAL_IFO_H2]->data->data[currPos];
      }
      if (chiSquare[LAL_IFO_L1])
      {
        currEvent->chisq_l = chiSquare[LAL_IFO_L1]->data->data[currPos];
      }
      if (chiSquare[LAL_IFO_T1])
      {
        currEvent->chisq_t = chiSquare[LAL_IFO_T1]->data->data[currPos];
      }
      if (chiSquare[LAL_IFO_V1])
      {
        currEvent->chisq_v = chiSquare[LAL_IFO_V1]->data->data[currPos];
      }
    }
  }

  REAL8 sigmasqCorrFac, invSigmaCorrFac;
  sigmasqCorrFac = ( (REAL8)fcTmplt->tmpltNorm);
  sigmasqCorrFac *= ( (REAL8)params->tempCorrFac);
  invSigmaCorrFac = 1. / pow(sigmasqCorrFac, 0.5);

  if (pValues[0])
    currEvent->amp_term_1 = pValues[0]->data->data[currPos] * invSigmaCorrFac;
  if (pValues[1])
    currEvent->amp_term_2 = pValues[1]->data->data[currPos] * invSigmaCorrFac;
  if (pValues[2])
    currEvent->amp_term_3 = pValues[2]->data->data[currPos] * invSigmaCorrFac;
  if (pValues[3])
    currEvent->amp_term_4 = pValues[3]->data->data[currPos] * invSigmaCorrFac;
    /* If here, also populate distance and other parameters */
    currEvent->distance = ( currEvent->amp_term_1*currEvent->amp_term_1 + \
                            currEvent->amp_term_2*currEvent->amp_term_2 + \
                            currEvent->amp_term_3*currEvent->amp_term_3 + \
                            currEvent->amp_term_4*currEvent->amp_term_4 );
    currEvent->distance = pow(2. / currEvent->distance, 0.5);
  if (pValues[4])
    currEvent->amp_term_5 = pValues[4]->data->data[currPos] * invSigmaCorrFac;
  if (pValues[5])
    currEvent->amp_term_6 = pValues[5]->data->data[currPos] * invSigmaCorrFac;
  if (pValues[6])
    currEvent->amp_term_7 = pValues[6]->data->data[currPos] * invSigmaCorrFac;
  if (pValues[7])
    currEvent->amp_term_8 = pValues[7]->data->data[currPos] * invSigmaCorrFac;
  if (pValues[8])
    currEvent->amp_term_9 = pValues[8]->data->data[currPos] * invSigmaCorrFac;
  if (pValues[9])
    currEvent->amp_term_10 = pValues[9]->data->data[currPos] * invSigmaCorrFac;
  /* Note that these two terms are only used for debugging
 *    * at the moment. When they are used properly they will be
 *       * moved into sane columns! For spin they give Amp*cos(Phi_0) and
 *          * Amp*sin(Phi_0). For non spinning the second is 0 and the
 *             * first is some arbitrary amplitude. */
  if (gammaBeta[0])
  {
    currEvent->g1quad = crectf( gammaBeta[0]->data->data[currPos], gammaBeta[1]->data->data[currPos] );
  }

  if (snrComps[LAL_IFO_G1])
  {
    currEvent->snr_g = snrComps[LAL_IFO_G1]->data->
        data[currPos+params->numBufferPoints+timeOffsetPoints[LAL_IFO_G1]];
    currEvent->sigmasq_g = PTFM[LAL_IFO_G1]->data[0] * sigmasqCorrFac;
  }
  if (snrComps[LAL_IFO_H1])
  {
    currEvent->snr_h1 = snrComps[LAL_IFO_H1]->data->
        data[currPos+params->numBufferPoints+timeOffsetPoints[LAL_IFO_H1]];
    currEvent->sigmasq_h1 = PTFM[LAL_IFO_H1]->data[0] * sigmasqCorrFac;
  }
  if (snrComps[LAL_IFO_H2])
  {
    currEvent->snr_h2 = snrComps[LAL_IFO_H2]->data->
        data[currPos+params->numBufferPoints+timeOffsetPoints[LAL_IFO_H2]];
    currEvent->sigmasq_h2 = PTFM[LAL_IFO_H2]->data[0] * sigmasqCorrFac;
  }
  if (snrComps[LAL_IFO_L1])
  {
    currEvent->snr_l = snrComps[LAL_IFO_L1]->data->
        data[currPos+params->numBufferPoints+timeOffsetPoints[LAL_IFO_L1]];
    currEvent->sigmasq_l = PTFM[LAL_IFO_L1]->data[0] * sigmasqCorrFac;
  }
  if (snrComps[LAL_IFO_T1])
  {
    currEvent->snr_t = snrComps[LAL_IFO_T1]->data->
        data[currPos+params->numBufferPoints+timeOffsetPoints[LAL_IFO_T1]];
    currEvent->sigmasq_t = PTFM[LAL_IFO_T1]->data[0] * sigmasqCorrFac;
  }
  if (snrComps[LAL_IFO_V1])
  {
    currEvent->snr_v = snrComps[LAL_IFO_V1]->data->
        data[currPos+params->numBufferPoints+timeOffsetPoints[LAL_IFO_V1]];
    currEvent->sigmasq_v = PTFM[LAL_IFO_V1]->data[0] * sigmasqCorrFac;
  }
  if (spinTrigger == 1)
  {
    if (params->numIFO == 1)
      currEvent->snr_dof = 6;
    else
      currEvent->snr_dof = 12;
  }
  else
  {
    currEvent->snr_dof = numDOF;
  }

  /* store ifos */
  if (params->numIFO == 1)
  {
    snprintf(currEvent->ifos, LIGOMETA_IFOS_MAX,\
              "%s", params->ifoName[0]);
  }
  else if(params->numIFO == 2)
  {
    snprintf(currEvent->ifos, LIGOMETA_IFOS_MAX, "%s%s",
             params->ifoName[0], params->ifoName[1]);
  }
  else if (params->numIFO == 3)
  {
    snprintf(currEvent->ifos, LIGOMETA_IFOS_MAX, "%s%s%s",
             params->ifoName[0], params->ifoName[1], params->ifoName[2]);
  }
  else if (params->numIFO == 4)
  {
    snprintf(currEvent->ifos, LIGOMETA_IFOS_MAX, "%s%s%s%s",
             params->ifoName[0], params->ifoName[1], params->ifoName[2],
             params->ifoName[3]);
  }
  if (params->faceOnStatistic == 2)
  {
    currEvent->inclination = LAL_PI/2.;
  }

  return currEvent;
}

UINT8 coh_PTF_add_sngl_triggers(
    struct coh_PTF_params   *params,
    SnglInspiralTable       **eventList,
    SnglInspiralTable       **thisEvent,
    REAL4TimeSeries         *cohSNR,
    FindChirpTemplate       *fcTmplt,
    InspiralTemplate        PTFTemplate,
    UINT8                   eventId,
    REAL4TimeSeries         **pValues,
    REAL4TimeSeries         **bankVeto,
    REAL4TimeSeries         **autoVeto,
    REAL4TimeSeries         **chiSquare,
    REAL8Array              **PTFM,
    UINT4                   startPoint,
    UINT4                   endPoint
)
{
  /* This function adds SnglInspiral events to the event list */
  
  UINT4 i;
  SnglInspiralTable *lastEvent = *thisEvent;
  SnglInspiralTable *currEvent = NULL;

  for (i = startPoint ; i < endPoint ; i++)
  {
    if (cohSNR->data->data[i])
    {
      currEvent = coh_PTF_create_sngl_event(params,cohSNR,fcTmplt,PTFTemplate,\
          &eventId,pValues,bankVeto,autoVeto,chiSquare,PTFM,i);

      /* Check trigger against trig times */
      if (coh_PTF_trig_time_check(params,currEvent->end_time,\
                                         currEvent->end_time))
      {
        if (currEvent->event_id)
        {
          LALFree(currEvent->event_id);
        }
        LALFree(currEvent);
        continue;
      }
      /* And add the trigger to the lists. IF it passes clustering! */
      if (!*eventList)
      {
        *eventList = currEvent;
        lastEvent = currEvent;
      }
      else
      {
        if (! params->clusterFlag)
        {
          lastEvent->next = currEvent;
          lastEvent = currEvent;
        }
        else if (coh_PTF_accept_sngl_trig_check(params,eventList,*currEvent) )
        {
          lastEvent->next = currEvent;
          lastEvent = currEvent;
        }
        else
        {
          if (currEvent->event_id)
          {
            LALFree(currEvent->event_id);
          }
          LALFree(currEvent);
        }
      }
    }
  }
  *thisEvent = lastEvent;
  return eventId;
}

SnglInspiralTable* coh_PTF_create_sngl_event(
    struct coh_PTF_params   *params,
    REAL4TimeSeries         *cohSNR,
    FindChirpTemplate       *fcTmplt,
    InspiralTemplate        PTFTemplate,
    UINT8                   *eventId,
    REAL4TimeSeries         **pValues,
    REAL4TimeSeries         **bankVeto,
    REAL4TimeSeries         **autoVeto,
    REAL4TimeSeries         **chiSquare,
    REAL8Array              **PTFM,
    UINT4                   currPos
)
{
  LIGOTimeGPS trigTime;
  UINT4 numDOF = 2;
  UINT4 ifoNumber;
  CHAR searchName[LIGOMETA_SEARCH_MAX];

  SnglInspiralTable *thisEvent;
  thisEvent = (SnglInspiralTable *)
      LALCalloc(1, sizeof(SnglInspiralTable));
  thisEvent->event_id = (EventIDColumn *)
      LALCalloc(1, sizeof(EventIDColumn));
  thisEvent->event_id->id=*eventId;
  (*eventId)++;
  /* Set end times */
  trigTime = cohSNR->epoch;
  XLALGPSAdd(&trigTime,currPos*cohSNR->deltaT);
  thisEvent->end_time = trigTime;
  thisEvent->end_time_gmst = fmod(XLALGreenwichMeanSiderealTime(
      &thisEvent->end_time), LAL_TWOPI) * 24.0 / LAL_TWOPI;     /* hours */

  /* Set SNR, chisqs, sigmasq, eff_distance */
  REAL8 sigmasqCorrFac;
  sigmasqCorrFac = ( (REAL8)fcTmplt->tmpltNorm);
  sigmasqCorrFac *= ( (REAL8)params->tempCorrFac);

  thisEvent->snr = cohSNR->data->data[currPos];
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if ( params->haveTrig[ifoNumber] )
    {
      thisEvent->chisq = chiSquare[ifoNumber]->data->data[currPos];
      thisEvent->bank_chisq = bankVeto[ifoNumber]->data->data[currPos];
      thisEvent->cont_chisq = autoVeto[ifoNumber]->data->data[currPos];
      thisEvent->sigmasq = PTFM[ifoNumber]->data[0] * sigmasqCorrFac;
    }
  }
  /* FIXME: Should be fixed so that this actually stores the DOF! */
  thisEvent->chisq_dof = params->numChiSquareBins;
  thisEvent->bank_chisq_dof = numDOF * params->BVsubBankSize;
  thisEvent->cont_chisq_dof = numDOF * params->numAutoPoints;
  thisEvent->eff_distance = sqrt( thisEvent->sigmasq ) / thisEvent->snr;
  /* FIXME: What is this? What is it used for? */
  thisEvent->event_duration = 0;
  /* FIXME: What does coa_phase actually mean here? */
  thisEvent->coa_phase = (REAL4)
      atan2( pValues[0]->data->data[currPos], pValues[1]->data->data[currPos]); 

  /* set the impulse time for the event FIXME:Check this! */
  thisEvent->template_duration = (REAL8) PTFTemplate.tC;
  
  /* record the ifo name for the event */
  snprintf(thisEvent->ifo, LIGOMETA_IFO_MAX, "%s", params->ifoName[0]);
  /* Set the channel name */
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if ( params->haveTrig[ifoNumber] )
    {
      strncpy( thisEvent->channel, (params->channel[ifoNumber]) + 3,
      (LALNameLength - 3) * sizeof(CHAR) );
    }
  }
  /* FIXME: NOt sure about this one either, inspiral just copies end_time */
  thisEvent->impulse_time = thisEvent->end_time;

  /* copy the template into the event */
  thisEvent->mass1   = (REAL4) PTFTemplate.mass1;
  thisEvent->mass2   = (REAL4) PTFTemplate.mass2;
  thisEvent->mtotal  = (REAL4) PTFTemplate.totalMass;
  thisEvent->mchirp  = (REAL4) PTFTemplate.chirpMass;
  thisEvent->eta     = (REAL4) PTFTemplate.eta;
  thisEvent->kappa   = (REAL4) PTFTemplate.kappa;
  thisEvent->chi     = (REAL4) PTFTemplate.chi;
  thisEvent->tau0    = (REAL4) PTFTemplate.t0;
  thisEvent->tau2    = (REAL4) PTFTemplate.t2;
  thisEvent->tau3    = (REAL4) PTFTemplate.t3;
  thisEvent->tau4    = (REAL4) PTFTemplate.t4;
  thisEvent->tau5    = (REAL4) PTFTemplate.t5;
  thisEvent->ttotal  = (REAL4) PTFTemplate.tC;
  thisEvent->f_final = (REAL4) PTFTemplate.fFinal;
  thisEvent->spin1z  = (REAL4) PTFTemplate.spin1[2];
  thisEvent->spin2z  = (REAL4) PTFTemplate.spin2[2];

  /* We can now memcpy the 10 metric co-efficients */
  memcpy (thisEvent->Gamma, PTFTemplate.Gamma, 10*sizeof(REAL4));

  /* Store the approximant */
  XLALInspiralGetApproximantString( searchName, LIGOMETA_SEARCH_MAX,
                                    params->approximant, params->order );
  memcpy( thisEvent->search, searchName,
      LIGOMETA_SEARCH_MAX * sizeof(CHAR) );

  return thisEvent;
}

UINT4 coh_PTF_accept_sngl_trig_check(
    struct coh_PTF_params   *params,
    SnglInspiralTable      **eventList,
    SnglInspiralTable      thisEvent
)
{
  SnglInspiralTable *currEvent = *eventList;
  LIGOTimeGPS time1,time2;
  UINT4 loudTrigBefore=0,loudTrigAfter=0;

  currEvent = *eventList;

  /* for each trigger, find out whether a louder trigger is within the
 *    * clustering time */
  time1.gpsSeconds=thisEvent.end_time.gpsSeconds;
  time1.gpsNanoSeconds = thisEvent.end_time.gpsNanoSeconds;
  while (currEvent)
  {
    time2.gpsSeconds=currEvent->end_time.gpsSeconds;
    time2.gpsNanoSeconds=currEvent->end_time.gpsNanoSeconds;
    if (fabs(XLALGPSDiff(&time1,&time2)) < params->clusterWindow)
    {
      if (thisEvent.snr < currEvent->snr\
          && (thisEvent.event_id->id != currEvent->event_id->id))
      {
        if ( XLALGPSDiff(&time1,&time2) < 0 )
          loudTrigBefore = 1;
        else
          loudTrigAfter = 1;
        if (loudTrigBefore && loudTrigAfter)
        {
          return 0;
        }
      }
    }
    currEvent = currEvent->next;
  }

  return 1;
}

void coh_PTF_cluster_sngl_triggers(
    struct coh_PTF_params   *params,
    SnglInspiralTable      **eventList,
    SnglInspiralTable      **thisEvent
)
{

  SnglInspiralTable *currEvent = *eventList;
  SnglInspiralTable *currEvent2 = NULL;
  SnglInspiralTable *newEvent = NULL;
  SnglInspiralTable *newEventHead = NULL;
  UINT4 triggerNum = 0;
  UINT4 lenTriggers = 0;
  UINT4 numRemovedTriggers = 0;

  /* find number of triggers */
  while (currEvent)
  {
    lenTriggers+=1;
    currEvent = currEvent->next;
  }

  currEvent = *eventList;
  UINT4 rejectTriggers[lenTriggers];

  /* for each trigger, find out whether a louder trigger is within the
 *    * clustering time */
  while (currEvent)
  {
    if (coh_PTF_accept_sngl_trig_check(params,eventList,*currEvent) )
    {
      rejectTriggers[triggerNum] = 0;
      triggerNum += 1;
    }
    else
    {
      rejectTriggers[triggerNum] = 1;
      triggerNum += 1;
      numRemovedTriggers += 1;
    }
    currEvent = currEvent->next;
  }

  currEvent = *eventList;
  triggerNum = 0;

  /* construct new event table with triggers to keep */
  while (currEvent)
  {
    if (! rejectTriggers[triggerNum])
    {
      if (! newEventHead)
      {
        newEventHead = currEvent;
        newEvent = currEvent;
      }
      else
      {
        newEvent->next = currEvent;
        newEvent = currEvent;
      }
      currEvent = currEvent->next;
    }
    else
    {
      if (currEvent->event_id)
      {
        LALFree(currEvent->event_id);
      }
      currEvent2 = currEvent->next;
      LALFree(currEvent);
      currEvent = currEvent2;
    }
    triggerNum+=1;
  }

  /* write new table over old one */
  if (newEvent)
  {
    newEvent->next = NULL;
    *eventList = newEventHead;
    *thisEvent = newEvent;
  }
}

void coh_PTF_cleanup(
    struct coh_PTF_params   *params,
    ProcessParamsTable      *procpar,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *psdplan,
    REAL4FFTPlan            *revplan,
    COMPLEX8FFTPlan         *invPlan,
    REAL4TimeSeries         **channel,
    REAL4FrequencySeries    **invspec,
    RingDataSegments        **segments,
    MultiInspiralTable      *events,
    SnglInspiralTable       *snglEvents,
    InspiralTemplate        *PTFbankhead,
    FindChirpTemplate       *fcTmplt,
    FindChirpTmpltParams    *fcTmpltParams,
    FindChirpInitParams     *fcInitParams,
    REAL8Array              **PTFM,
    REAL8Array              **PTFN,
    COMPLEX8VectorSequence  **PTFqVec,
    REAL4                   *timeOffsets,
    REAL4                   *slidTimeOffsets,
    REAL4                   *Fplus,
    REAL4                   *Fcross,
    REAL4                   *Fplustrig,
    REAL4                   *Fcrosstrig,
    CohPTFSkyPositions      *skyPoints,
    TimeSlide               *time_slide_head,
    TimeSlideVectorList     *longTimeSlideList,
    TimeSlideVectorList     *shortTimeSlideList,
    REAL4                   *timeSlideVectors,
    LALDetector             **detectors,
    INT8                    *slideIDList,
    TimeSlideSegmentMapTable *time_slide_map_head,
    SegmentTable            *segment_table_head
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
  while ( snglEvents )
  {
    SnglInspiralTable *thisSnglEvent;
    thisSnglEvent = snglEvents;
    snglEvents = snglEvents->next;
    if ( thisSnglEvent->event_id )
    {
      LALFree( thisSnglEvent->event_id );
    }
    LALFree( thisSnglEvent );
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
  if ( psdplan )
    XLALDestroyREAL4FFTPlan( psdplan );
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
  if ( slidTimeOffsets)
    LALFree( slidTimeOffsets );
  if ( Fplus )
    LALFree( Fplus );
  if ( Fcross )
    LALFree( Fcross );
  if ( Fplustrig )
    LALFree( Fplustrig );
  if ( Fcrosstrig )
    LALFree( Fcrosstrig );

  if (skyPoints)
  {
    if (skyPoints->data)
      LALFree(skyPoints->data);
    LALFree(skyPoints);
  }
  if (time_slide_head)
    XLALDestroyTimeSlideTable(time_slide_head);
  if (time_slide_map_head)
    XLALDestroyTimeSlideSegmentMapTable(time_slide_map_head);
  if (segment_table_head)
    XLALDestroySegmentTable(segment_table_head);
  if (longTimeSlideList)
    LALFree(longTimeSlideList);
  if (shortTimeSlideList)
    LALFree(shortTimeSlideList);
  if (timeSlideVectors)
    LALFree(timeSlideVectors);

  for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if (detectors[ifoNumber])
      LALFree(detectors[ifoNumber]);
  }

  if (slideIDList)
    LALFree(slideIDList);

}


/* gets the forward fft plan */
REAL4FFTPlan *coh_PTF_get_fft_fwdplan( struct coh_PTF_params *params )
{
  REAL4FFTPlan *plan = NULL;
  if ( params->segmentDuration > 0.0 )
  {
    UINT4 segmentLength;
    segmentLength = params->numTimePoints;
    plan = XLALCreateForwardREAL4FFTPlan( segmentLength, params->fftLevel );
  }
  return plan;
}

/* gets the psd forward fft plan */
REAL4FFTPlan *coh_PTF_get_fft_psdplan( struct coh_PTF_params *params )
{
  REAL4FFTPlan *plan = NULL;
  if ( params->psdSegmentDuration > 0.0 )
  {
    UINT4 segmentLength;
    segmentLength = floor( params->psdSegmentDuration * params->sampleRate
                           + 0.5 );
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
    segmentLength = params->numTimePoints;
    plan = XLALCreateReverseREAL4FFTPlan( segmentLength, params->fftLevel );
  }
  return plan;
}

/* gets the inverse fft plan */
COMPLEX8FFTPlan *coh_PTF_get_fft_invplan( struct coh_PTF_params *params )
{
  COMPLEX8FFTPlan *plan = NULL;
  if ( params->segmentDuration > 0.0 )
  {
    UINT4 segmentLength;
    segmentLength = params->numTimePoints;
    plan = XLALCreateReverseCOMPLEX8FFTPlan( segmentLength, params->fftLevel );
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
    verbose("Generated full sky grid with %d points, ",
            skyPoints->numPoints);
    verbose("parsing for time-delay degeneracy\n");
    CohPTFSkyPositions *parsedSkyPoints = NULL; 
    parsedSkyPoints = coh_PTF_parse_time_delays(skyPoints, params);
    if (skyPoints->data)
      LALFree(skyPoints->data);
    if (skyPoints)
      LALFree(skyPoints);
    skyPoints = parsedSkyPoints;
  }
  verbose("Generated final sky grid with %d points, ",
          skyPoints->numPoints);
  return skyPoints;
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
  if ((! params->singlePolFlag) && (params->numIFO != 1))
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
  gsl_vector_free(npPos);
  gsl_vector_free(trigPos);
  gsl_vector_free(axis);

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
      if (fabs(timeDelay-timeDelays[i]) < dt)
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
  gsl_vector_free(pos);
  gsl_vector_free(rotPos);
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
      gsl_vector_free(locations[ifoNumber]);
  }
  gsl_vector_free(normal);
  gsl_vector_free(northPole);
  gsl_vector_free(axis);

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
    /* WARNING: THIS FUNCTION WILL NOT WORK WITH SHORT SLIDES. DO NOT ATTEMPT
     * TO SLIDE INJECTION TRIGGERS WITHOUT FIXING THIS FUNCTION FIRST! */ 

    /* define variables */
    LIGOTimeGPS injTime, segmentStart, segmentEnd;
    UINT4 injSamplePoint, injWindow, tmpStart, tmpEnd;
    REAL8 injDiff;
    INT8 startDiff, endDiff;
    SimInspiralTable *thisInject = NULL;

    /* set variables */
    segmentStart = *epoch;
    segmentEnd   = *epoch;
    XLALGPSAdd(&segmentEnd, params->strideDuration);
    thisInject = params->injectList;

    tmpStart = tmpEnd = 0;

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
            if (tmpStart)
            {
                verbose("warning: multiple injections in this segment.\n");
                tmpStart = params->analStartPoint;
                tmpEnd = params->analEndPoint;
            }
            else
            {
                injDiff = (REAL8) ((XLALGPSToINT8NS(&injTime) - \
                                    XLALGPSToINT8NS(&segmentStart)) / 1E9);
                injSamplePoint = floor(injDiff * params->sampleRate + 0.5);
                injSamplePoint += params->analStartPoint;
                injWindow = floor(params->injSearchWindow * params->sampleRate
                                  + 1);
                tmpStart = injSamplePoint - injWindow;
                if (tmpStart < params->analStartPoint)
                    tmpStart = params->analStartPoint;
                tmpEnd = injSamplePoint + injWindow + 1;
                if (tmpEnd > params->analEndPoint)
                    tmpEnd = params->analEndPoint;
                verbose("Found analysis segment at [%d,%d).\n", tmpStart, tmpEnd);
            }
        }
        thisInject = thisInject->next;
    }
    *start = tmpStart;
    *end = tmpEnd;  
}

UINT4 coh_PTF_trig_time_check(
    struct coh_PTF_params *params,
    LIGOTimeGPS segStartTime,
    LIGOTimeGPS segEndTime)
{
  INT8 currTimeNS;
  /* Does the segment end before the trigStartTime */
  if (params->trigStartTimeNS)
  {
    currTimeNS = segEndTime.gpsSeconds * 1E9 + segEndTime.gpsNanoSeconds;
    if (params->trigStartTimeNS > currTimeNS)
    {
      return 1;
    }
  }
  /* Does the segment start before the trigEndTime */
  if (params->trigEndTimeNS)
  {
    currTimeNS = segStartTime.gpsSeconds * 1E9 + segStartTime.gpsNanoSeconds;
    if (params->trigStartTimeNS > currTimeNS)
    {
      return 1;
    }
  }
  /* Otherwise pass */
  return 0;
}

UINT4 checkInjectionMchirp(
    struct coh_PTF_params *params,
    InspiralTemplate *tmplt,
    LIGOTimeGPS *epoch
    )
{
  /* define variables */
  LIGOTimeGPS injTime, segmentStart, segmentEnd;
  REAL8 tmpltMchirp,injMchirp,mchirpDiff,mchirpWin;
  INT8 startDiff, endDiff;
  SimInspiralTable *thisInject = NULL;
  UINT4 passMchirpCheck;

  /* set variables */
  segmentStart = *epoch;
  segmentEnd   = *epoch;
  XLALGPSAdd(&segmentEnd, params->strideDuration);
  passMchirpCheck = 2;
  thisInject = params->injectList;
  
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
      if (passMchirpCheck != 2)
      {
        verbose("warning: multiple injections in this segment.\n");
        passMchirpCheck = 1;
        break;
      }
      injMchirp = thisInject->mchirp;
      tmpltMchirp = tmplt->chirpMass;
      mchirpDiff = (injMchirp - tmpltMchirp)/tmpltMchirp;
      /* The mchirp window is increased with mchirp */
      if (injMchirp < 2)
        mchirpWin = params->injMchirpWindow;
      else if (injMchirp < 3)
        mchirpWin = params->injMchirpWindow * 2.5;
      else if (injMchirp < 4)
        mchirpWin = params->injMchirpWindow * 5;
      else
        // Note that I haven't tuned this above Mchirp = 6
        mchirpWin = params->injMchirpWindow * 10;

      if (fabs(mchirpDiff) > mchirpWin)
        passMchirpCheck = 0;
      else
        passMchirpCheck = 1;
    }
    thisInject = thisInject->next;
  }

  if (passMchirpCheck == 2)
  {
    verbose("WARNING: No injections found in this segment!? Not analysing\n");
    passMchirpCheck = 0;
  }
  return passMchirpCheck;
}

void coh_PTF_set_null_input_REAL4TimeSeries(
  REAL4TimeSeries** timeSeries,
  UINT4 length
)
{
  UINT4 i;
  for (i = 0 ; i < length ; i++)
  {
    timeSeries[i] = NULL;
  }
}

void coh_PTF_set_null_input_REAL4FrequencySeries(
  REAL4FrequencySeries** freqSeries,
  UINT4 length
)
{
  UINT4 i;
  for (i = 0 ; i < length ; i++)
  {
    freqSeries[i] = NULL;
  }
}

void coh_PTF_set_null_input_RingDataSegments(
  RingDataSegments** segment,
  UINT4 length
)
{
  UINT4 i;
  for (i = 0 ; i < length ; i++)
  {
    segment[i] = NULL;
  }
}

void coh_PTF_set_null_input_REAL8Array(
  REAL8Array** array,
  UINT4 length
)
{
  UINT4 i;
  for (i = 0 ; i < length ; i++)
  {
    array[i] = NULL;
  }
}

void coh_PTF_set_null_input_COMPLEX8VectorSequence(
  COMPLEX8VectorSequence** vecSeq,
  UINT4 length
)
{
  UINT4 i;
  for (i = 0 ; i < length ; i++)
  {
    vecSeq[i] = NULL;
  }
}

void coh_PTF_set_null_input_REAL4(
  REAL4** array,
  UINT4 length
)
{
  UINT4 i;
  for (i = 0 ; i < length ; i++)
  {
    array[i] = NULL;
  }
}

void coh_PTF_set_null_input_UINT4(
  UINT4** array,
  UINT4 length
)
{
  UINT4 i;
  for (i = 0 ; i < length ; i++)
  {
    array[i] = NULL;
  }
}

void coh_PTF_set_null_input_LALDetector(
  LALDetector** detector,
  UINT4 length
)
{
  UINT4 i;
  for (i = 0 ; i < length ; i++)
  {
    detector[i] = NULL;
  }
}

