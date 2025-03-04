#include "config.h"
#include "coh_PTF.h"

#define PROGRAM_NAME "lalapps_coh_PTF_inspiral"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE   "$Source$"
#define CVS_DATE     "$Date$"

/* This function should be migrated to option.c */
/* warning: returns a pointer to a static variable... not reenterant */
/* only call this routine once to initialize params! */
/* also do not attempt to free this pointer! */
static struct coh_PTF_params *coh_PTF_get_params(int argc, char **argv)
{
  static struct coh_PTF_params params;
  static char   programName[]  = PROGRAM_NAME;
  static char   cvsRevision[]  = CVS_REVISION;
  static char   cvsSource[]    = CVS_SOURCE;
  static char   cvsDate[]      = CVS_DATE;
  coh_PTF_parse_options(&params, argc, argv);
  coh_PTF_params_sanity_check(&params); /* this also sets various params */
  coh_PTF_params_inspiral_sanity_check(&params);
  params.programName = programName;
  params.cvsRevision = cvsRevision;
  params.cvsSource   = cvsSource;
  params.cvsDate     = cvsDate;
  return &params;
}

static int XLALCountMultiInspiralTable(MultiInspiralTable **head, UINT4 len)
{
	UINT4 ui;
        int length = 0;
        MultiInspiralTable *temp;
	/* count the number of events in the list */
        for (ui=0; ui < len; ui++)
        {
          for(temp = head[ui]; temp; temp = temp->next)
		length++;
        }

	return(length);
}

int main(int argc, char **argv)
{

  /* Declarations of parameters */
  INT4  i,j,k;
  UINT4 uj,sp,slideNum;

  /* process structures */
  struct coh_PTF_params    *params                  = NULL;
  ProcessParamsTable       *procpar                 = NULL;

  /* sky position+time slide structures */
  UINT4                    numSkyPoints,currAnalStart,currAnalEnd;
  CohPTFSkyPositions       *skyPoints               = NULL;
  INT8  *slideIDList,currSlideID;
  TimeSlide *time_slide_head = NULL;
  TimeSlideVectorList *longTimeSlideList = NULL;
  TimeSlideVectorList * shortTimeSlideList = NULL;
  SegmentTable *segment_table_head = NULL;
  TimeSlideSegmentMapTable *time_slide_map_head = NULL;

  /* FFT structures */
  REAL4FFTPlan             *fwdplan                 = NULL;
  REAL4FFTPlan             *psdplan                 = NULL;
  REAL4FFTPlan             *revplan                 = NULL;
  COMPLEX8FFTPlan          *invplan                 = NULL;

  /* input data and spectrum storage */
  REAL4TimeSeries          *channel[LAL_NUM_IFO+1];
  REAL4FrequencySeries     *invspec[LAL_NUM_IFO+1];
  RingDataSegments         *segments[LAL_NUM_IFO+1];
  INT4                     numSegments              = 0;

  /* template counters */
  INT4                     numTmplts                = 0;
  INT4                     numSpinTmplts            = 0;
  INT4                     numNoSpinTmplts          = 0;

  /* template indexes */
  INT4                     startTemplate            = -1;
  INT4                     stopTemplate             = -1;

  /* template and findchirp data structures */
  InspiralTemplate         *PTFSpinTemplate         = NULL;
  InspiralTemplate         *PTFNoSpinTemplate       = NULL;
  InspiralTemplate         *PTFtemplate             = NULL;
  InspiralTemplate         *PTFbankhead             = NULL;
  FindChirpTemplate        *fcTmplt                 = NULL;
  FindChirpTemplate        *bankFcTmplts            = NULL;
  FindChirpTmpltParams     *fcTmpltParams           = NULL;
  UINT4                    ifoNumber,spinTemplate;
  REAL8Array               *PTFM[LAL_NUM_IFO+1];
  REAL8Array               *PTFN[LAL_NUM_IFO+1];
  COMPLEX8VectorSequence   *PTFqVec[LAL_NUM_IFO+1];
  UINT4                    *acceptPointList         = NULL;
  UINT4                    numAcceptPoints;

  /* triggered sky position and sensitivity structures */
  LIGOTimeGPS              segStartTime;
  LIGOTimeGPS              segEndTime;
  segStartTime.gpsSeconds = 0;
  segStartTime.gpsNanoSeconds = 0;
  struct timeval           startTime;
  LALDetector              *detectors[LAL_NUM_IFO+1];
  REAL4                    *timeOffsets,*slidTimeOffsets;
  REAL4                    *Fplus;
  REAL4                    *Fcross;
  REAL4                    *Fplustrig;
  REAL4                    *Fcrosstrig;
  REAL4                    *timeSlideVectors;

  /* coherent statistic structures */
  REAL4TimeSeries          *cohSNR                  = NULL;
  REAL4TimeSeries          *pValues[10];
  REAL4TimeSeries          *snrComps[LAL_NUM_IFO];
  UINT4                    *snglAcceptPoints[LAL_NUM_IFO];
  UINT4                    snglAcceptCount[LAL_NUM_IFO];
  REAL4TimeSeries          *gammaBeta[2];
  REAL4TimeSeries          *nullSNR                 = NULL;
  REAL4TimeSeries          *traceSNR                = NULL;

  /* consistency test structures */
  REAL4TimeSeries          *bankVeto[LAL_NUM_IFO+1];
  REAL4TimeSeries          *autoVeto[LAL_NUM_IFO+1];
  REAL4TimeSeries          *chiSquare[LAL_NUM_IFO+1];
  struct bankDataOverlaps  *chisqOverlaps           = NULL;
  struct bankDataOverlaps  *chisqSnglOverlaps       = NULL;
  REAL4                    *frequencyRangesPlus[LAL_NUM_IFO+1];
  REAL4                    *frequencyRangesCross[LAL_NUM_IFO+1];
  UINT4 subBankSize = 0;
  struct bankTemplateOverlaps *bankNormOverlaps = NULL;
  struct bankComplexTemplateOverlaps *bankOverlaps = NULL;
  struct bankDataOverlaps *dataOverlaps = NULL;
  UINT4 timeStepPoints = 0;
  struct bankComplexTemplateOverlaps *autoTempOverlaps = NULL;
  REAL4                   **overlapCont = NULL;
  REAL4                   **snglOverlapCont = NULL;

  /* output event structures */
  MultiInspiralTable       **eventList               = NULL;
  MultiInspiralTable       **thisEvent               = NULL;
  MultiInspiralTable       *finalEvents             = NULL;
  SnglInspiralTable        *snglEventList           = NULL;
  SnglInspiralTable        *snglThisEvent           = NULL;
  UINT8                    eventId                  = 0;
  INT4                     timeDiff;

  /*------------------------------------------------------------------------*
   * initialise                                                             *
   *------------------------------------------------------------------------*/

  gettimeofday(&startTime, NULL);

  /* set error handlers to abort on error */
  set_abrt_on_error();

  /*
   * no lal mallocs before this! *
   */

  /* options are parsed and debug level is set here... */
  params  = coh_PTF_get_params(argc, argv);

  /* create process params */
  procpar = create_process_params(argc, argv, PROGRAM_NAME);

  verbose("Read input params %ld \n", timeval_subtract(&startTime));

  /* create forward and reverse fft plans */
  fwdplan = coh_PTF_get_fft_fwdplan(params);
  psdplan = coh_PTF_get_fft_psdplan(params);
  revplan = coh_PTF_get_fft_revplan(params);
  invplan = coh_PTF_get_fft_invplan(params);

  verbose("Made fft plans %ld \n", timeval_subtract(&startTime));

  /* NULL out pointers where necessary */
  coh_PTF_set_null_input_REAL4TimeSeries(channel,LAL_NUM_IFO+1);
  coh_PTF_set_null_input_REAL4FrequencySeries(invspec,LAL_NUM_IFO+1);
  coh_PTF_set_null_input_RingDataSegments(segments,LAL_NUM_IFO+1);
  coh_PTF_set_null_input_REAL4TimeSeries(pValues,10);
  coh_PTF_set_null_input_REAL4TimeSeries(snrComps,LAL_NUM_IFO);
  coh_PTF_set_null_input_UINT4(snglAcceptPoints,LAL_NUM_IFO);
  coh_PTF_set_null_input_REAL4TimeSeries(bankVeto,LAL_NUM_IFO+1);
  coh_PTF_set_null_input_REAL4TimeSeries(autoVeto,LAL_NUM_IFO+1);
  coh_PTF_set_null_input_REAL4TimeSeries(chiSquare,LAL_NUM_IFO+1);
  coh_PTF_set_null_input_REAL4TimeSeries(gammaBeta,2);
  coh_PTF_set_null_input_REAL4(frequencyRangesPlus,LAL_NUM_IFO+1);
  coh_PTF_set_null_input_REAL4(frequencyRangesCross,LAL_NUM_IFO+1);
  coh_PTF_set_null_input_LALDetector(detectors,LAL_NUM_IFO+1);

  /* allocate memory for time offsets and detector responses */
  timeOffsets = LALCalloc(1, LAL_NUM_IFO*sizeof(REAL4));
  slidTimeOffsets = LALCalloc(1, LAL_NUM_IFO*sizeof(REAL4));
  Fplus       = LALCalloc(1, LAL_NUM_IFO*sizeof(REAL4));
  Fcross      = LALCalloc(1, LAL_NUM_IFO*sizeof(REAL4));
  Fplustrig   = LALCalloc(1, LAL_NUM_IFO*sizeof(REAL4));
  Fcrosstrig  = LALCalloc(1, LAL_NUM_IFO*sizeof(REAL4));

  /* Initialize template and filtering structures */
  coh_PTF_initialize_structures(params,&fcTmplt,&fcTmpltParams,\
                                PTFM,PTFN,PTFqVec,fwdplan);

  /*------------------------------------------------------------------------*
   * read the data, generate segments and the PSD                           *
   *------------------------------------------------------------------------*/


  numSegments = coh_PTF_data_condition(params,channel,invspec,segments,\
                         fwdplan,psdplan,revplan,&timeSlideVectors,startTime);


  /*------------------------------------------------------------------------*
   * Create a list of time slide ids for each segment and create time slide *
   * table.                                                                 *
   *------------------------------------------------------------------------*/

  slideIDList = LALCalloc(1, numSegments*sizeof(INT8));
  coh_PTF_create_time_slide_table(params,slideIDList,segments,&time_slide_head,\
                                  &time_slide_map_head,&segment_table_head,\
                                  &longTimeSlideList,&shortTimeSlideList,\
                                  timeSlideVectors,numSegments);

  /*------------------------------------------------------------------------*
   * Determine the list of sky points.                                      *
   * Determine time delays and response functions for central point         *
   * This is computed for all detectors, even if not being analyzed         *
   *------------------------------------------------------------------------*/

  /* generate sky points array */
  skyPoints = coh_PTF_generate_sky_points(params);
  numSkyPoints = skyPoints->numPoints;

  /*------------------------------------------------------------------------*
   * Initialize the trigger storage structures                              *
   *------------------------------------------------------------------------*/

  timeDiff = params->endTime.gpsSeconds - params->startTime.gpsSeconds + 1;
  eventList = LALCalloc(1, timeDiff*params->numShortSlides*sizeof(MultiInspiralTable*));
  thisEvent = LALCalloc(1, timeDiff*params->numShortSlides*sizeof(MultiInspiralTable*));
  for (i = 0; (i < (INT4) (timeDiff*params->numShortSlides)); i++)
  {
    eventList[i] = NULL;
    thisEvent[i] = NULL;
  }


  /* loop over ifos if doing triggered search and determine the time-offset */
  /* and detector responses for the "preferred" sky location. For GRBs where */
  /* a central sky point is provided, this will be that point. For cases */
  /* where a list of points is given this will be the first point. */
  /* For all sky search there is no "preferred" sky location. */
  /* This preferred location is used to center the chi-squared */
  /* (I THINK!!!) */
  /* The preferred location is also used for the null stream, if active */
  if ((params->skyLooping != ALL_SKY) &&
       (params->skyLooping != TWO_DET_ALL_SKY))
  {
    coh_PTF_calculate_det_stuff(params,detectors,timeOffsets,Fplustrig,\
                                Fcrosstrig,skyPoints,0);
  }

  /*------------------------------------------------------------------------*
   * Construct the null stream, its segments and its PSD                    *
   *------------------------------------------------------------------------*/

  if (params->doNullStream)
  {
    coh_PTF_setup_null_stream(params,channel,invspec,\
            segments,Fplustrig,Fcrosstrig,timeOffsets,\
            fwdplan,revplan,psdplan,timeSlideVectors,startTime);
  }

  /*------------------------------------------------------------------------*
   * At this point we can discard the calibrated data, only the segments    *
   * and spectrum are needed now                                            *
   *------------------------------------------------------------------------*/

  for(ifoNumber = 0; ifoNumber < (LAL_NUM_IFO+1); ifoNumber++)
  {
    if (channel[ifoNumber])
    {
      XLALDestroyREAL4Vector(channel[ifoNumber]->data);
      LALFree(channel[ifoNumber]);
      channel[ifoNumber] = NULL;
    }
  }

  /*------------------------------------------------------------------------*
   * Read in the tmpltbank xml files                                        *
   *------------------------------------------------------------------------*/

  /* read spinning bank */
  if (params->spinBank)
  {
    numSpinTmplts = InspiralTmpltBankFromLIGOLw(&PTFSpinTemplate,
      params->spinBankName,startTemplate, stopTemplate);
    if (numSpinTmplts != 0)
    {
      PTFtemplate = PTFSpinTemplate;
      numTmplts = numSpinTmplts;
    }
    else
      params->spinBank = 0;
  }

  /* read non-spinning bank */
  if (params->noSpinBank)
  {
    numNoSpinTmplts = InspiralTmpltBankFromLIGOLw(&PTFNoSpinTemplate,
      params->noSpinBankName,startTemplate, stopTemplate);
    if (numNoSpinTmplts != 0)
    {
      PTFtemplate = PTFNoSpinTemplate;
      numTmplts = numNoSpinTmplts;
    }
    else
      params->noSpinBank = 0;
  }

  /* If both banks present combine them and mark where to swap over */
  if (params->spinBank && params->noSpinBank)
  {
    for (i = 0; (i < numNoSpinTmplts); PTFtemplate = PTFtemplate->next, i++)
    {
      if (i == (numNoSpinTmplts - 1))
      {
        PTFtemplate->next = PTFSpinTemplate;
        numTmplts = numSpinTmplts + numNoSpinTmplts;
      }
    }
    PTFtemplate = PTFNoSpinTemplate;
  }
  PTFbankhead = PTFtemplate;

  /*------------------------------------------------------------------------*
   * Allocate RAM for the various time series that get calculated           *
   *------------------------------------------------------------------------*/

  coh_PTF_initialize_time_series(params,segStartTime,\
          params->lowTemplateFrequency,&cohSNR,&nullSNR,&traceSNR,bankVeto,\
          autoVeto,chiSquare,snrComps,pValues,gammaBeta,numSpinTmplts);
  /* FIXME: Move into function above if this works */
  overlapCont = LALCalloc(1, LAL_NUM_IFO*sizeof(*overlapCont));
  snglOverlapCont = LALCalloc(1, LAL_NUM_IFO*sizeof(*overlapCont));
  acceptPointList = LALCalloc(params->numAnalPoints, sizeof(UINT4));
  for (ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    overlapCont[ifoNumber] = NULL;
    snglOverlapCont[ifoNumber] = NULL;
    if (params->haveTrig[ifoNumber])
    {
      snglAcceptPoints[ifoNumber] = \
          LALCalloc(params->numAnalPointsBuf, sizeof(UINT4));
    }
  }

  verbose("Initialized storage arrays at %ld\n",timeval_subtract(&startTime));


  /*------------------------------------------------------------------------*
   * Initialise bank veto - This function does the following:
   *  - Generate the set of bank veto templates and store \tilde{h}
   *  - Calculate the overlaps between each pair of templates
   *------------------------------------------------------------------------*/

  if (params->doBankVeto)
  {
    subBankSize = coh_PTF_initialize_bank_veto(params,&bankNormOverlaps,\
            &bankOverlaps,&dataOverlaps,&bankFcTmplts,fcTmplt,fcTmpltParams,\
            invspec,startTime);
  }

  /*------------------------------------------------------------------------*
   * initialise auto veto
   * - Create the structures needed for the auto veto, if necessary
   *------------------------------------------------------------------------*/

  if (params->doAutoVeto)
  {
    timeStepPoints = coh_PTF_initialize_auto_veto(params,&autoTempOverlaps,\
                                                  startTime);
  }

  /*------------------------------------------------------------------------*
   * find gravitational waves
   *------------------------------------------------------------------------*/

  UINT4 numPoints = params->numTimePoints;

  /* This is the primary loop over segments */
  for (j = 0; j < numSegments; ++j)
  {
    /* Reset the template list to the first one */
    PTFtemplate = PTFbankhead;

    /* Determine the epoch of this segment */
    for (ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
    {
      if (params->haveTrig[ifoNumber])
      {
        segStartTime = segments[ifoNumber]->sgmnt[j].epoch;
        segEndTime = segments[ifoNumber]->sgmnt[j].epoch;
        break;
      }
    }
    /* We only analyse middle half so add duration/4 to epoch */
    XLALGPSAdd(&segStartTime, params->analStartTime);
    XLALGPSAdd(&segEndTime, params->analEndTime);

    /* Test if trig-start and trig-end options overlap this segment at all */
    if (coh_PTF_trig_time_check(params,segStartTime,segEndTime))
      continue;

    if (params->doBankVeto)
    {
      /* For every segment we need to calculate the overlap between bank veto
       * templates and the data for use in bank veto calculation */
      coh_PTF_bank_veto_segment_setup(params,dataOverlaps,bankFcTmplts,\
                                      segments,PTFqVec,invplan,j,startTime);
    }

    /* This is the primary loop over templates in the bank */
    for (i = 0; (i < numTmplts); PTFtemplate = PTFtemplate->next, i++)
    {
      /* If running injections, check whether to analyse this template*/
      if ( params->injectFile && params->injMchirpWindow )
      {
        if (! checkInjectionMchirp(params,PTFtemplate,&segStartTime))
        {
          verbose("Injection not within mchirp window for segment %d, template %d at %ld \n", j, i, timeval_subtract(&startTime));
          continue;
        }
      }

      /* Determine if this template is non-spinning */
      if (i >= numNoSpinTmplts)
        spinTemplate = 1;
      else
        spinTemplate = 0;

      /* This value is used for template generation */
      PTFtemplate->fLower = params->lowTemplateFrequency;

      /* This function generates the template */
      coh_PTF_template(fcTmplt,PTFtemplate,fcTmpltParams);

      /* Put the template in the array used by the coh_PTF filtering */
      if (params->approximant != FindChirpPTF)
      {
        for (uj = 0 ; uj < (numPoints/2 +1) ; uj++)
        {
          fcTmplt->PTFQtilde->data[uj] = fcTmplt->data->data[uj];
        }
      }

      if (spinTemplate)
        verbose("Generated spin template %d at %ld\n",\
                i,timeval_subtract(&startTime));
      else
        verbose("Generated nospin template %d at %ld\n",\
                i,timeval_subtract(&startTime));

      /* Reset the epoch here and memset to 0 all entries */
      coh_PTF_reset_time_series(params,segStartTime,\
          cohSNR,nullSNR,traceSNR,bankVeto,\
          autoVeto,chiSquare,snrComps,pValues,gammaBeta,numSpinTmplts);
      verbose("Initialized storage arrays for segment %d at %ld\n",\
          j, timeval_subtract(&startTime));

      /* Calculate single detector filters */
      coh_PTF_calculate_single_detector_filters(params,fcTmplt,invspec,PTFM,\
              PTFqVec,snrComps,snglAcceptPoints,snglAcceptCount,segments,\
              invplan,spinTemplate,j);
      verbose("Calculated sngl filters for segment %d template %d at %ld\n",\
          j, i, timeval_subtract(&startTime));

      /* Calculate single detector bank veto filters for this template */
      if (params->doBankVeto)
      {
        coh_PTF_calculate_bank_veto_template_filters(params,bankFcTmplts,\
                fcTmplt,invspec,bankOverlaps);
        verbose("Calculated bank-veto filters for template %d at %ld\n",\
            i, timeval_subtract(&startTime));
      }

      /* Calculate single detector auto veto filters */
      if (params->doAutoVeto)
      {
        coh_PTF_calculate_auto_veto_template_filters(params,fcTmplt,\
                autoTempOverlaps,invspec,invplan,timeStepPoints);
        verbose("Calculated auto-veto filters for template %d at %ld\n",\
            i, timeval_subtract(&startTime));
      }

      /* Calculate null stream filters */
      if (params->doNullStream)
      {
        coh_PTF_calculate_null_stream_filters(params,fcTmplt,invspec,PTFM,\
              PTFqVec,segments,invplan,spinTemplate,j);
        verbose(\
          "Calculated null stream filters for segment %d template %d at %ld\n",\
          j, i, timeval_subtract(&startTime));
      }

      verbose("Begin loop over sky points at %ld \n",
              timeval_subtract(&startTime));

      /* Primary loop over sky points */
      for (sp = 0; sp < numSkyPoints ; sp++)
      {
/*        if (! ((sp == 0) || (sp == 3285)))
        {
          continue;
        } */
        /* Calculate offsets and responses for this sky point */
        coh_PTF_calculate_det_stuff(params,detectors,timeOffsets,Fplus,\
                            Fcross,skyPoints,sp);
        /* Loop over short slides */
        for (slideNum = 0; slideNum < params->numShortSlides ; slideNum++)
        {
/*          if (! ((slideNum == 30) || (slideNum == 21)))
          {
            continue;
          }*/
          /* Update the offsets */
          for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
          {
            if (params->haveTrig[ifoNumber])
            {
              slidTimeOffsets[ifoNumber] = timeOffsets[ifoNumber] + \
                  shortTimeSlideList[slideNum].timeSlideVectors[ifoNumber];
            }
            else
            {
              slidTimeOffsets[ifoNumber] = 0;
            }
          }
          /* Determine slide ID */
          currSlideID = slideIDList[j]*params->numShortSlides;
          currSlideID += shortTimeSlideList[slideNum].timeSlideID;

          /* FIXME: This is a hack when doing faceOn+faceAway analysis */
          if (params->faceAwayAnalysis && params->faceOnAnalysis)
          {
            params->faceOnStatistic = 1;
          }

          // This function calculates the cohSNR time series and all of the
          // signal based vetoes as appropriate
          numAcceptPoints = coh_PTF_statistic(\
                         cohSNR, PTFM, PTFqVec, params, spinTemplate,
                         slidTimeOffsets, Fplus, Fcross,
                         j, pValues, gammaBeta, snrComps, nullSNR,
                         traceSNR, bankVeto, autoVeto,
                         chiSquare, subBankSize, bankOverlaps,
                         bankNormOverlaps, dataOverlaps, autoTempOverlaps,
                         fcTmplt, invspec, segments, invplan,
                         &chisqOverlaps,&chisqSnglOverlaps, frequencyRangesPlus,
                         frequencyRangesCross, overlapCont, snglOverlapCont,
                         startTime,
                         shortTimeSlideList[slideNum].analStartPoint,
                         shortTimeSlideList[slideNum].analEndPoint,\
                         snglAcceptPoints,snglAcceptCount,acceptPointList);

          verbose("Made coherent statistic for segment %d, template %d, "
                  "sky point %d at %ld \n", j, i, sp,
                  timeval_subtract(&startTime));

          currAnalStart = shortTimeSlideList[slideNum].analStartPoint - \
                      params->analStartPoint;
          currAnalEnd = shortTimeSlideList[slideNum].analEndPoint - \
                      params->analStartPoint;

          /* This function construct triggers from loud events */
          if (! params->writeSnglInspiralTable)
          {
            eventId = coh_PTF_add_triggers(params, eventList, thisEvent,
                                         cohSNR, fcTmplt, *PTFtemplate, eventId,
                                         spinTemplate,
                                         pValues, gammaBeta, snrComps,
                                         nullSNR, traceSNR, bankVeto,
                                         autoVeto, chiSquare, PTFM,
                                         skyPoints->data[sp].longitude,
                                         skyPoints->data[sp].latitude,
                                         currSlideID, slidTimeOffsets,
                                         acceptPointList,numAcceptPoints,
                                         slideNum, timeDiff,
                                         params->startTime.gpsSeconds);
          }
          else
          {
            eventId = coh_PTF_add_sngl_triggers(params, &snglEventList,\
                           &snglThisEvent,cohSNR,fcTmplt,*PTFtemplate,eventId,\
                           pValues,bankVeto,autoVeto,chiSquare,PTFM,\
                           currAnalStart,currAnalEnd);
          }

          /* FIXME: Also part of the faceAway + faceOn hack */
          if (params->faceAwayAnalysis && params->faceOnAnalysis)
          {
            params->faceOnStatistic = 2;

            numAcceptPoints = coh_PTF_statistic(\
                         cohSNR, PTFM, PTFqVec, params, spinTemplate,
                         slidTimeOffsets, Fplus, Fcross,
                         j, pValues, gammaBeta, snrComps, nullSNR,
                         traceSNR, bankVeto, autoVeto,
                         chiSquare, subBankSize, bankOverlaps,
                         bankNormOverlaps, dataOverlaps, autoTempOverlaps,
                         fcTmplt, invspec, segments, invplan,
                         &chisqOverlaps,&chisqSnglOverlaps, frequencyRangesPlus,
                         frequencyRangesCross, overlapCont, snglOverlapCont,
                         startTime,
                         shortTimeSlideList[slideNum].analStartPoint,
                         shortTimeSlideList[slideNum].analEndPoint,\
                         snglAcceptPoints,snglAcceptCount,acceptPointList);

            verbose("Made coherent statistic for segment %d, template %d, "
                    "sky point %d at %ld \n", j, i, sp,
                    timeval_subtract(&startTime));

            eventId = coh_PTF_add_triggers(params, eventList, thisEvent,
                                         cohSNR, fcTmplt, *PTFtemplate, eventId,
                                         spinTemplate,
                                         pValues, gammaBeta, snrComps,
                                         nullSNR, traceSNR, bankVeto,
                                         autoVeto, chiSquare, PTFM,
                                         skyPoints->data[sp].longitude,
                                         skyPoints->data[sp].latitude,
                                         currSlideID, slidTimeOffsets,
                                         acceptPointList,numAcceptPoints,
                                         slideNum, timeDiff,
                                         params->startTime.gpsSeconds);
          } /* End of if faceaway and faceon block */

          if (vrbflg)
          {
            if (! params->writeSnglInspiralTable)
            {
              params->numEvents = XLALCountMultiInspiralTable(eventList,
                                              timeDiff*params->numShortSlides);
            }
            else
            {
              params->numEvents = XLALCountSnglInspiral(snglEventList);
            }
            verbose("There are currently %d triggers.\n", params->numEvents);
            verbose("Generated triggers for segment %d, template %d, sky point %d, short slide %d at %ld \n", j, i, sp, slideNum, timeval_subtract(&startTime));
          }
        }/* End loop over time slides */
      }/* End loop over sky points*/

      /* Free memory for temporary chisq products */

      if (chisqOverlaps)
      {
        for(uj = 0; uj < 2*params->numChiSquareBins; uj++)
        {
          for(k = 0; k < LAL_NUM_IFO; k++)
          {
            if (chisqOverlaps[uj].PTFqVec[k])
            {
              XLALDestroyCOMPLEX8VectorSequence(chisqOverlaps[uj].PTFqVec[k]);
            }
          }
        }
        LALFree(chisqOverlaps);
        chisqOverlaps = NULL;
      }
      if (chisqSnglOverlaps)
      {
        for(uj = 0; uj < params->numChiSquareBins; uj++)
        {
          for(k = 0; k < LAL_NUM_IFO; k++)
          {
            if (chisqSnglOverlaps[uj].PTFqVec[k])
            {
              XLALDestroyCOMPLEX8VectorSequence(chisqSnglOverlaps[uj].PTFqVec[k]);
            }
          }
        }
        LALFree(chisqSnglOverlaps);
        chisqSnglOverlaps = NULL;
      }
      for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO+1; ifoNumber++)
      {
        if (frequencyRangesPlus[ifoNumber])
        {
          LALFree(frequencyRangesPlus[ifoNumber]);
          frequencyRangesPlus[ifoNumber] = NULL;
        }
        if (frequencyRangesCross[ifoNumber])
        {
          LALFree(frequencyRangesCross[ifoNumber]);
          frequencyRangesCross[ifoNumber] = NULL;
        }
      }
      for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if (overlapCont[ifoNumber])
        {
          LALFree(overlapCont[ifoNumber]);
          overlapCont[ifoNumber] = NULL;
        }
        if (snglOverlapCont[ifoNumber])
        {
          LALFree(snglOverlapCont[ifoNumber]);
          snglOverlapCont[ifoNumber] = NULL;
        }
      }

    } /* End of loop over templates */

  } /* End of loop over segments */

  /* calulate number of events and cluster if needed */
  if (! params->writeSnglInspiralTable)
  {
    params->numEvents = XLALCountMultiInspiralTable(eventList,
                                              timeDiff*params->numShortSlides);
  }
  else
  {
    params->numEvents = XLALCountSnglInspiral(snglEventList);
  }

  verbose("There are %d total triggers before cluster.\n", params->numEvents);
  if ( params->clusterFlag )
  {
    if (! params->writeSnglInspiralTable)
    {
      coh_PTF_cluster_triggers(params,eventList,&finalEvents,
                               params->numShortSlides, timeDiff);
      params->numEvents = XLALCountMultiInspiralTable(&finalEvents, 1);
    }
    else
    {
      coh_PTF_cluster_sngl_triggers(params,&snglEventList,&snglThisEvent);
      params->numEvents = XLALCountSnglInspiral(snglEventList);
    }
    verbose("There are %d total triggers after cluster.\n", params->numEvents);
  }

  /* Output events to xml */
  coh_PTF_output_events_xml(params->outputFile, finalEvents, snglEventList,\
                            params->injectList, procpar, time_slide_head,\
                            time_slide_map_head, segment_table_head, params);

  /* Everything that follows is memory cleanup */
  coh_PTF_destroy_time_series(cohSNR,nullSNR,traceSNR,bankVeto,autoVeto,\
          chiSquare,pValues,gammaBeta,snrComps);

  coh_PTF_cleanup(params,procpar,fwdplan,psdplan,revplan,invplan,channel,
      invspec,segments,finalEvents,snglEventList,\
      PTFbankhead,fcTmplt,fcTmpltParams,
      PTFM,PTFN,PTFqVec,timeOffsets,slidTimeOffsets,Fplus,Fcross,\
      Fplustrig,Fcrosstrig,skyPoints,time_slide_head,longTimeSlideList,
      shortTimeSlideList,timeSlideVectors,detectors, slideIDList,\
      time_slide_map_head,segment_table_head);
  LALFree(eventList);

  coh_PTF_free_veto_memory(params,bankNormOverlaps,bankFcTmplts,bankOverlaps,\
      dataOverlaps,autoTempOverlaps);

  /* FIXME: Move into above if this works */
  for (ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if (params->haveTrig[ifoNumber])
    {
      LALFree(snglAcceptPoints[ifoNumber]);
    }
  }
  LALFree(acceptPointList);

  verbose("Generated output xml file, cleaning up and exiting at %ld \n",
      timeval_subtract(&startTime));

  /* And check for memory leaks*/
  LALCheckMemoryLeaks();
  return 0;
}

UINT4 coh_PTF_statistic(
    REAL4TimeSeries         *cohSNR,
    REAL8Array              *PTFM[LAL_NUM_IFO+1],
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
    struct coh_PTF_params   *params,
    UINT4                   spinTemplate,
    REAL4                   *timeOffsets,
    REAL4                   *Fplus,
    REAL4                   *Fcross,
    INT4                    segmentNumber,
    REAL4TimeSeries         *pValues[10],
    UNUSED REAL4TimeSeries         *gammaBeta[2],
/* NOTE: This is unused because the spin record extrinsic parameters function
 * is broken. When fixed, this will be used, so DO NOT DELETE */
    REAL4TimeSeries         *snrComps[LAL_NUM_IFO],
    REAL4TimeSeries         *nullSNR,
    REAL4TimeSeries         *traceSNR,
    REAL4TimeSeries         *bankVeto[LAL_NUM_IFO+1],
    REAL4TimeSeries         *autoVeto[LAL_NUM_IFO+1],
    REAL4TimeSeries         *chiSquare[LAL_NUM_IFO+1],
    UINT4                   subBankSize,
    struct bankComplexTemplateOverlaps *bankOverlaps,
    struct bankTemplateOverlaps *bankNormOverlaps,
    struct bankDataOverlaps *dataOverlaps,
    struct bankComplexTemplateOverlaps *autoTempOverlaps,
    FindChirpTemplate       *fcTmplt,
    REAL4FrequencySeries    *invspec[LAL_NUM_IFO+1],
    RingDataSegments        *segments[LAL_NUM_IFO+1],
    COMPLEX8FFTPlan         *invPlan,
    struct bankDataOverlaps **chisqOverlapsP,
    struct bankDataOverlaps **chisqSnglOverlapsP,
    REAL4 *frequencyRangesPlus[LAL_NUM_IFO+1],
    REAL4 *frequencyRangesCross[LAL_NUM_IFO+1],
    REAL4                   **overlapCont,
    REAL4                   **snglOverlapCont,
    struct timeval          startTime,
    UINT4                   segStartPoint,
    UINT4                   segEndPoint,
    UINT4                   **snglAcceptPoints,
    UINT4                   *snglAcceptCount,
    UINT4                   *acceptPointList
)

{
  verbose("Entering the statistic loop at %ld \n",
          timeval_subtract(&startTime));

  /* This function generates the SNR for every point in time and, where
   * appropriate calculates the desired signal based vetoes. */

  /* Begin with all the declarations */
  UINT4  csVecLength,csVecLengthTwo;
  UINT4  i, j, k, vecLength, vecLengthTwo;
  INT4   timeOffsetPoints[LAL_NUM_IFO],numPointCheck;
  REAL4 *v1p,*v2p; /*snglSNRthresh;*/
  /*REAL4 cohSNRThreshold, cohSNRThresholdSq;*/
  REAL4 *powerBinsPlus[LAL_NUM_IFO+1],*powerBinsCross[LAL_NUM_IFO+1];
  REAL4 *snrData;
  UINT4 currPointLoc,numAcceptPoints;
  struct bankCohTemplateOverlaps *bankCohOverlaps,*autoCohOverlaps;
  gsl_matrix *eigenvecs,*eigenvecsNull,*Autoeigenvecs;
  gsl_vector *eigenvals,*eigenvalsNull,*Autoeigenvals;
  eigenvecs = NULL;
  eigenvecsNull = NULL;
  Autoeigenvecs = NULL;
  eigenvals = NULL;
  eigenvalsNull = NULL;
  Autoeigenvals = NULL;
  /* FIXME: the 50s below seem to hardcode a limit on the number of templates
   * this should not be hardcoded. Note that this value is hardcoded in some
   * function declarations as well as here! Double pointers will fix this*/
  gsl_matrix *Bankeigenvecs[50];
  gsl_vector *Bankeigenvals[50];

  struct bankDataOverlaps *chisqOverlaps = *chisqOverlapsP;
  struct bankDataOverlaps *chisqSnglOverlaps = *chisqSnglOverlapsP;

  /* Code works slightly differently if spin/non spin and single/coherent */
  /* First we set the various vector lengths */
  if (spinTemplate)
    vecLength = 5;
  else
    vecLength = 1;
  if (params->numIFO == 1 || params->singlePolFlag || params->faceOnStatistic)
    vecLengthTwo = vecLength;
  else
    vecLengthTwo = 2* vecLength;

  csVecLength = 1;
  csVecLengthTwo = 2;
  if (params->numIFO == 1 || params->singlePolFlag || params->faceOnStatistic)
    csVecLengthTwo = 1;

  /* Initialize pointers, some initialize to NULL */
  Autoeigenvals = NULL;
  bankCohOverlaps = NULL;
  autoCohOverlaps = NULL;
  for (i = 0; i < 50; i++)
  {
    Bankeigenvecs[i] = NULL;
    Bankeigenvals[i] = NULL;
  }
  for(i = 0; i < LAL_NUM_IFO+1; i++)
  {
    powerBinsPlus[i] = NULL;
    powerBinsCross[i] = NULL;
  }
  eigenvecs   = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  eigenvals   = gsl_vector_alloc(vecLengthTwo);
  eigenvecsNull = gsl_matrix_alloc(vecLength,vecLength);
  eigenvalsNull = gsl_vector_alloc(vecLength);
  v1p = LALCalloc(vecLengthTwo , sizeof(REAL4));
  v2p = LALCalloc(vecLengthTwo , sizeof(REAL4));

  /* Pick the relevant SNR threshold */
  /*cohSNRThreshold = params->threshold;
  if (spinTemplate)
    cohSNRThreshold = params->spinThreshold;
  snglSNRthresh = params->snglSNRThreshold;
  cohSNRThresholdSq = cohSNRThreshold*cohSNRThreshold;*/

  /* This function takes the (Q_i|Q_j) matrices, combines it across the ifos
   * and returns the eigenvalues and eigenvectors of this new matrix.
   * We later rotate and rescale the (Q_i|s) values such that in the new basis
   * this matrix will be the identity matrix.
   * For non-spin this describes the rotation into the dominant polarization
   */
  coh_PTF_calculate_bmatrix(params,eigenvecs,eigenvals,Fplus,Fcross,PTFM,\
                            vecLength,vecLengthTwo,vecLength);

  /* If required also calculate these eigenvalues/vectors for the null stream*/
  if (params->doNullStream)
  {
    coh_PTF_calculate_null_stream_norms(vecLength,eigenvecsNull,eigenvalsNull,\
                                        PTFM);
  }

  /* This function takes the time offset in seconds and converts to time offset
   * in data points (rounded) */
  coh_PTF_convert_time_offsets_to_points(params,timeOffsets,timeOffsetPoints);

  /* If we have injections we only want to analyse the time around the
   * injection. Here we figure out what that time shou/rd be. No injections
   * equates to analyse the whole segment.
   */

  if ( params->injectFile && params->analyzeInjSegsOnly )
  {
    findInjectionSegment(&segStartPoint, &segEndPoint, &cohSNR->epoch, params);
  }

  if (! segStartPoint)
  {
    segStartPoint = params->analStartPoint;
    segEndPoint = params->analEndPoint;
  }

  verbose("-->Begin SNR calculation at %ld \n",timeval_subtract(&startTime));

  snrData = cohSNR->data->data;

  coh_PTF_calculate_coherent_SNR(params,snrData,pValues,snrComps,\
                                 timeOffsetPoints,PTFqVec,Fplus,Fcross,\
                                 eigenvecs,eigenvals,segStartPoint,segEndPoint,\
                                 vecLength,vecLengthTwo,spinTemplate,\
                                 snglAcceptPoints,snglAcceptCount);

  verbose("-->Calculated all SNRs at %ld \n",timeval_subtract(&startTime));

  numPointCheck = floor(params->timeWindow/cohSNR->deltaT + 0.5);

  numAcceptPoints = coh_PTF_template_time_series_cluster(params,
                                       cohSNR,acceptPointList,timeOffsetPoints,\
                                       numPointCheck,\
                                       segStartPoint - params->analStartPoint,\
                                       segEndPoint - params->analStartPoint,\
                                       snglAcceptPoints,snglAcceptCount);

  verbose("-->Done template clustering at %ld \n",timeval_subtract(&startTime));

  /* Now we calculate all the extrinsic parameters and signal based vetoes
   * Only calculated if this will be a trigger
   */

  for (j = 0; j < numAcceptPoints; ++j) /* loop over time */
  { /* We only loop over points that are not already rejected for speed */
    currPointLoc = acceptPointList[j];
    i = currPointLoc + params->analStartPoint;
    /* Check if point is going to be rejected */
    if (snrData[currPointLoc])
    {
      /* First sbv to be calculated is the null stream SNR. */
      if (params->doNullStream)
      {
        coh_PTF_calculate_null_stream_snr(params,nullSNR,PTFqVec,\
            eigenvecsNull,eigenvalsNull,spinTemplate,vecLength,i,\
            currPointLoc);
      }

      /* Next up is Trace SNR */
      if (params->doTraceSNR)
      {
        coh_PTF_calculate_trace_snr(params,traceSNR,PTFqVec,eigenvecs,\
            eigenvals,Fplus,Fcross,timeOffsetPoints,spinTemplate,vecLength,\
            vecLengthTwo,i,currPointLoc);
      }

      /* Next is the bank veto */
      if (params->doBankVeto)
      {
        if (params->numIFO != 1)
        {
          /* Begin by calculating variouse overlaps that are needed.
           * This only needs to doing once (and is only done once), but might
           * be better put outside of a loop over time. However, do not want
           * to calculate this if *no* points will be stored in this instance
           */
          coh_PTF_bank_veto_coh_setup(params,Bankeigenvecs,Bankeigenvals,\
              &bankCohOverlaps,bankOverlaps,Fplus,Fcross,PTFM,\
              bankNormOverlaps,csVecLength,csVecLengthTwo,vecLength);
          /* In this function all the filters are combined to produce the
           * value of the bank veto. */
          bankVeto[LAL_NUM_IFO]->data->data[currPointLoc] = \
              coh_PTF_calculate_bank_veto(params->numTimePoints,i,subBankSize,\
                  Fplus,Fcross,params,bankCohOverlaps,NULL,dataOverlaps,NULL,\
                  PTFqVec,NULL,timeOffsetPoints,Bankeigenvecs,Bankeigenvals,\
                  LAL_NUM_IFO,csVecLength,csVecLengthTwo);
        }
        /* As well as the coherent bank veto calculated above, we calculate
         * the single detector bank veto */
        if (params->doSnglChiSquared)
        {
          for(k = 0; k < LAL_NUM_IFO; k++)
          {
            if (params->haveTrig[k])
            {
              bankVeto[k]->data->data[currPointLoc] = \
                  coh_PTF_calculate_bank_veto(params->numTimePoints,i,\
                      subBankSize,Fplus,Fcross,params,NULL,bankOverlaps,\
                      dataOverlaps,bankNormOverlaps,PTFqVec,PTFM,\
                      timeOffsetPoints,NULL,NULL,k,1,1);
            }
          }
        }
      }

      /* Now we do the auto veto */
      if (params->doAutoVeto)
      {
        if (params->numIFO!=1)
        {
          /* As with bank_veto, we begin by calculating the various coherent
           * overlaps that are needed. Same caveats as with bank veto */
          coh_PTF_auto_veto_coh_setup(params,&Autoeigenvecs,&Autoeigenvals,\
              &autoCohOverlaps,autoTempOverlaps,Fplus,Fcross,PTFM,\
              csVecLength,csVecLengthTwo,vecLength);
          /* Auto veto is calculated */
          autoVeto[LAL_NUM_IFO]->data->data[currPointLoc] = \
              coh_PTF_calculate_auto_veto(params->numTimePoints,i,Fplus,Fcross,\
                  params,autoCohOverlaps,NULL,PTFqVec,NULL,timeOffsetPoints,\
                  Autoeigenvecs,Autoeigenvals,LAL_NUM_IFO,csVecLength,\
                  csVecLengthTwo);
        }
        /* Also do single detector auto veto */
        if (params->doSnglChiSquared)
        {
          for(k = 0; k < LAL_NUM_IFO; k++)
          {
            if (params->haveTrig[k])
            {
              autoVeto[k]->data->data[currPointLoc] = \
                  coh_PTF_calculate_auto_veto(params->numTimePoints,i,Fplus,\
                      Fcross,params,NULL,autoTempOverlaps,PTFqVec,PTFM,\
                      timeOffsetPoints,NULL,NULL,k,1,1);
            }
          }
        }

      }
    }
  }

  verbose("-->Calculated most vetoes at %ld \n",timeval_subtract(&startTime));

  /* And do the loop again to calculate chi square */

  if (params->doChiSquare)
  {
    for (j = 0; j < numAcceptPoints; ++j) /* loop over time */
    { /* We only loop over points that are not already rejected for speed */
      currPointLoc = acceptPointList[j];
      i = currPointLoc + params->analStartPoint;
      if (snrData[currPointLoc])
      {
        /* Test whether to do chi^2 */
        if (params->chiSquareCalcThreshold)
        {
          /* FIXME: Does not work for single detector runs */
          if (params->numIFO!=1)
          {
            if ( coh_PTF_test_veto_vals(params,cohSNR,nullSNR,bankVeto,\
                                        autoVeto,currPointLoc ) )
            {
              chiSquare[LAL_NUM_IFO]->data->data[currPointLoc] = 0;
              continue;
            }
          }
        }
        /* If no problems then calculate chi squared */
        if (params->numIFO != 1)
        {
          /* As with bank_veto, we begin by calculating the various coherent
           * overlaps that are needed. Same caveats as with bank veto
           * For chi square there are two types of variables here some of the
           * variables are only calculated once for all sky points, at the first
           * sky point they are needed. This includes the computationally
           * expensive filters, which means we must also use the same frequency
           * ranges for all sky points. Because of this we recalculate the
           * unequal power bins each time.
           */
          coh_PTF_chi_square_coh_setup(params,&Autoeigenvecs,\
              &Autoeigenvals,frequencyRangesPlus,frequencyRangesCross,\
              powerBinsPlus,powerBinsCross,overlapCont,&chisqOverlaps,fcTmplt,\
              invspec,segments,Fplus,Fcross,PTFM,invPlan,segmentNumber,\
              csVecLength,csVecLengthTwo,vecLength);

          /* Calculate chi square here */
          chiSquare[LAL_NUM_IFO]->data->data[currPointLoc] = \
              coh_PTF_calculate_chi_square(params,i,chisqOverlaps,\
                  PTFqVec,NULL,Fplus,Fcross,timeOffsetPoints,Autoeigenvecs,\
                  Autoeigenvals,powerBinsPlus[LAL_NUM_IFO],\
                  powerBinsCross[LAL_NUM_IFO],LAL_NUM_IFO,csVecLength,\
                  csVecLengthTwo);
        }
        if (params->doSnglChiSquared)
        {
          /* Begin with the setup, this is only done once. As with the
           * coherent chi squared, the filters are reused for every sky point
           */
          coh_PTF_chi_square_sngl_setup(params,frequencyRangesPlus,\
              frequencyRangesCross,powerBinsPlus,powerBinsCross,\
              snglOverlapCont,&chisqSnglOverlaps,fcTmplt,invspec,segments,\
              PTFM,invPlan,segmentNumber);
          for(k = 0; k < LAL_NUM_IFO; k++)
          {
            if (params->haveTrig[k])
            {
              /* Calculate chi squared */
              chiSquare[k]->data->data[currPointLoc] = \
                  coh_PTF_calculate_chi_square(params,i,\
                      chisqSnglOverlaps,PTFqVec,PTFM,Fplus,Fcross,\
                      timeOffsetPoints,NULL,NULL,powerBinsPlus[k],\
                      powerBinsCross[k],k,1,1);
            }
          } /* End loop over ifos */
        }
      }
    } /* End loop over time */
  }

  verbose("-->Calculated chi squared at %ld \n",timeval_subtract(&startTime));

  /* Free memory of stuff that will not be reused. */
  if (params->doBankVeto)
  {
    for (j = 0 ; j < subBankSize+1 ; j++)
    {
      if (Bankeigenvecs[j])
        gsl_matrix_free(Bankeigenvecs[j]);
      if (Bankeigenvals[j])
        gsl_vector_free(Bankeigenvals[j]);
    }
    if (bankCohOverlaps)
    {
      for (j = 0 ; j < subBankSize ; j++)
      {
        gsl_matrix_free(bankCohOverlaps[j].rotReOverlaps);
        gsl_matrix_free(bankCohOverlaps[j].rotImOverlaps);
      }
      LALFree(bankCohOverlaps);
    }
  }

  if (params->doAutoVeto)
  {
    if (Autoeigenvecs)
      gsl_matrix_free(Autoeigenvecs);
      Autoeigenvecs = NULL;
    if (Autoeigenvals)
      gsl_vector_free(Autoeigenvals);
      Autoeigenvals = NULL;
    if (autoCohOverlaps)
    {
      for (j = 0 ; j < params->numAutoPoints ; j++)
      {
        gsl_matrix_free(autoCohOverlaps[j].rotReOverlaps);
        gsl_matrix_free(autoCohOverlaps[j].rotImOverlaps);
      }
      LALFree(autoCohOverlaps);
    }
  }

  if (params->doChiSquare)
  {
    if (Autoeigenvecs)
      gsl_matrix_free(Autoeigenvecs);
    if (Autoeigenvals)
      gsl_vector_free(Autoeigenvals);
    for(k = 0; k < LAL_NUM_IFO+1; k++)
    {
      if (powerBinsPlus[k])
        LALFree(powerBinsPlus[k]);
      if (powerBinsCross[k])
        LALFree(powerBinsCross[k]);
    }
  }

  LALFree(v1p);
  LALFree(v2p);
  gsl_matrix_free(eigenvecs);
  gsl_vector_free(eigenvals);
  gsl_matrix_free(eigenvecsNull);
  gsl_vector_free(eigenvalsNull);

  *chisqOverlapsP = chisqOverlaps;
  *chisqSnglOverlapsP = chisqSnglOverlaps;

  verbose("-->Completed memory cleanup and exiting loop at %ld \n",\
          timeval_subtract(&startTime));

  return numAcceptPoints;
}

UINT8 coh_PTF_add_triggers(
    struct coh_PTF_params   *params,
    MultiInspiralTable      **eventList,
    MultiInspiralTable      **thisEvent,
    REAL4TimeSeries         *cohSNR,
    FindChirpTemplate       *fcTmplt,
    InspiralTemplate        PTFTemplate,
    UINT8                   eventId,
    UINT4                   spinTrigger,
    REAL4TimeSeries         *pValues[10],
    REAL4TimeSeries         *gammaBeta[2],
    REAL4TimeSeries         *snrComps[LAL_NUM_IFO],
    REAL4TimeSeries         *nullSNR,
    REAL4TimeSeries         *traceSNR,
    REAL4TimeSeries         *bankVeto[LAL_NUM_IFO+1],
    REAL4TimeSeries         *autoVeto[LAL_NUM_IFO+1],
    REAL4TimeSeries         *chiSquare[LAL_NUM_IFO+1],
    REAL8Array              *PTFM[LAL_NUM_IFO+1],
    REAL4                   rightAscension,
    REAL4                   declination,
    INT8                    slideId,
    REAL4                   *timeOffsets,
    UINT4                   *acceptPointList,
    UINT4                   numAcceptPoints,
    UINT4                   slideNum,
    INT4                    timeDiff,
    INT4                    startTime
)
{
  // This function adds a trigger to the event list

  UINT4 i,j;
  UINT4 currTimeDiff, currStorageID;
  INT4   timeOffsetPoints[LAL_NUM_IFO];
  MultiInspiralTable *lastEvent = NULL;
  MultiInspiralTable *currEvent = NULL;

  coh_PTF_convert_time_offsets_to_points(params,timeOffsets,timeOffsetPoints);

  for (j = 0; j < numAcceptPoints; ++j)
  {
    i = acceptPointList[j];
    if (cohSNR->data->data[i])
    {
      currEvent = coh_PTF_create_multi_event(params,cohSNR,fcTmplt,PTFTemplate,\
          &eventId,spinTrigger,pValues,gammaBeta,snrComps,nullSNR,traceSNR,\
          bankVeto,autoVeto,chiSquare,PTFM,rightAscension,declination,slideId,\
          timeOffsetPoints,i);
      /* Important to zero out the cohSNR array as we go */
      cohSNR->data->data[i] = 0.;

      /* Check trigger against trig times */
      if (coh_PTF_trig_time_check(params,currEvent->end_time,\
                                         currEvent->end_time))
      {
        LALFree(currEvent);
        continue;
      }
      currTimeDiff = (currEvent->end_time.gpsSeconds - startTime);
      currStorageID = timeDiff * slideNum + currTimeDiff;
      /* And add the trigger to the lists. IF it passes clustering! */
      if (!eventList[currStorageID])
      {
        eventList[currStorageID] = currEvent;
        thisEvent[currStorageID] = currEvent;
      }
      else
      {
        lastEvent = thisEvent[currStorageID];
        if (! params->clusterFlag)
        {
          lastEvent->next = currEvent;
          thisEvent[currStorageID] = currEvent;
        }
        else if (coh_PTF_accept_trig_check(params,eventList,*currEvent,
                                           timeDiff, currTimeDiff,
                                           currStorageID) )
        {
          lastEvent->next = currEvent;
          thisEvent[currStorageID] = currEvent;
        }
        else
        {
          LALFree(currEvent);
        }
      }
    }
  }
  return eventId;
}

void coh_PTF_cluster_triggers(
    struct coh_PTF_params   *params,
    MultiInspiralTable      **eventList,
    MultiInspiralTable      **newEventHead,
    UINT4                   numSlides,
    INT4                    timeDiff
)
{
  UINT4 slideNum, currTimeDiff, currStorageID;
  UINT4 triggerNum, lenTriggers;
  UINT4 **rejectTriggers;
  rejectTriggers = LALCalloc(1, numSlides*timeDiff*sizeof(UINT4*));

  MultiInspiralTable *currEvent=NULL;
  MultiInspiralTable *currEvent2=NULL;
  MultiInspiralTable *newEvent=NULL;

  for (slideNum = 0; slideNum < numSlides; slideNum++)
  {
    for (currTimeDiff = 0; currTimeDiff < (UINT4)timeDiff; currTimeDiff++)
    {
      currStorageID = timeDiff * slideNum + currTimeDiff;

      currEvent = eventList[currStorageID];
      currEvent2 = NULL;
      triggerNum = 0;
      lenTriggers = 0;

      /* find number of triggers */
      while (currEvent)
      {
        lenTriggers+=1;
        currEvent = currEvent->next;
      }

      currEvent = eventList[currStorageID];
      rejectTriggers[currStorageID] = LALCalloc(1, lenTriggers*sizeof(UINT4));

      /* for each trigger, find out whether a louder trigger is within the
       * clustering time */
      while (currEvent)
      {
        if (coh_PTF_accept_trig_check(params,eventList,*currEvent,
                                      timeDiff, currTimeDiff,
                                      currStorageID) )
        {
          rejectTriggers[currStorageID][triggerNum] = 0;
          triggerNum += 1;
        }
        else
        {
          rejectTriggers[currStorageID][triggerNum] = 1;
          triggerNum += 1;
        }
        currEvent = currEvent->next;
      }

    }
  }

  for (slideNum = 0; slideNum < numSlides; slideNum++)
  {
    for (currTimeDiff = 0; currTimeDiff < (UINT4)timeDiff; currTimeDiff++)
    {
      triggerNum = 0;
      currStorageID = timeDiff * slideNum + currTimeDiff;
      currEvent = eventList[currStorageID];
      /* construct new event table with triggers to keep */
      while (currEvent)
      {
        if (! rejectTriggers[currStorageID][triggerNum])
        {
          if (! *newEventHead)
          {
            *newEventHead = currEvent;
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
          currEvent2 = currEvent->next;
          LALFree(currEvent);
          currEvent = currEvent2;
        }
        triggerNum+=1;
      }
      LALFree(rejectTriggers[currStorageID]);
    }
  }
  LALFree(rejectTriggers);

  /* write new table over old one */
  if (newEvent)
  {
    newEvent->next = NULL;
  }
}

UINT4 coh_PTF_accept_trig_check(
    struct coh_PTF_params   *params,
    MultiInspiralTable      **eventList,
    MultiInspiralTable      thisEvent,
    INT4                    timeDiff,
    UINT4                   currTimeDiff,
    UINT4                   currStorageID
)
{
  MultiInspiralTable *currEvent = *eventList;
  LIGOTimeGPS time1,time2;
  UINT4 loudTrigBefore=0,loudTrigAfter=0;
  INT4 checkOffset;
  REAL8 GPSDiff;

  time1.gpsSeconds=thisEvent.end_time.gpsSeconds;
  time1.gpsNanoSeconds = thisEvent.end_time.gpsNanoSeconds;

  for (checkOffset=-1; checkOffset < 2; checkOffset++)
  {
    if ((currTimeDiff + checkOffset) == 0)
    {
      continue;
    }
    if ((INT4)(currTimeDiff + checkOffset) == timeDiff)
    {
      continue;
    }

    currEvent = eventList[currStorageID+checkOffset];

    /* for each trigger, find out whether a louder trigger is within the
     * clustering time */
    while (currEvent)
    {
      time2.gpsSeconds=currEvent->end_time.gpsSeconds;
      time2.gpsNanoSeconds=currEvent->end_time.gpsNanoSeconds;
      if (thisEvent.time_slide_id == currEvent->time_slide_id)
      {
        GPSDiff = XLALGPSDiff(&time1,&time2);
        if (fabs(GPSDiff) < params->clusterWindow)
        {
          if (thisEvent.snr_dof == currEvent->snr_dof)
          {
            if (thisEvent.snr < currEvent->snr\
                && (thisEvent.event_id != currEvent->event_id))
            {
              /* If at identical time, return 0 */
              if (GPSDiff == 0)
                return 0;
              else if ( GPSDiff < 0 )
                loudTrigBefore = 1;
              else
                loudTrigAfter = 1;

              if (loudTrigBefore && loudTrigAfter)
              {
                return 0;
              }
            }
          }
        }
      }
      currEvent = currEvent->next;
    }
  }

  return 1;
}
