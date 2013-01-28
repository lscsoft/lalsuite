#define LAL_USE_OLD_COMPLEX_STRUCTS
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

int main(int argc, char **argv)
{

  /* Declarations of parameters */
  INT4  i,j,k;
  UINT4 ui,uj,sp;
  char name[LALNameLength];
  CHAR ifoName[LIGOMETA_IFO_MAX];

  /* process structures */
  struct coh_PTF_params    *params                  = NULL;
  ProcessParamsTable       *procpar                 = NULL;

  /* sky position structures */
  UINT4                    numSkyPoints;
  CohPTFSkyPositions       *skyPoints               = NULL;

  /* FFT structures */
  REAL4FFTPlan             *fwdplan                 = NULL;
  REAL4FFTPlan             *psdplan                 = NULL;
  REAL4FFTPlan             *revplan                 = NULL;
  COMPLEX8FFTPlan          *invPlan                 = NULL;

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

  /* bank structures */
  UINT4                    UNUSED spinBank          = 0;
  char                     spinFileName[256];
  char                     noSpinFileName[256];

  /* template and findchirp data structures */
  InspiralTemplate         *PTFSpinTemplate         = NULL;
  InspiralTemplate         *PTFNoSpinTemplate       = NULL;
  InspiralTemplate         *PTFtemplate             = NULL;
  InspiralTemplate         *PTFbankhead             = NULL;
  FindChirpTemplate        *fcTmplt                 = NULL;
  InspiralTemplate         *PTFBankTemplates        = NULL;
  InspiralTemplate         *PTFBankvetoHead         = NULL;
  FindChirpTemplate        *bankFcTmplts            = NULL;
  FindChirpTmpltParams     *fcTmpltParams           = NULL;
  FindChirpInitParams      *fcInitParams            = NULL;
  UINT4                    numPoints,ifoNumber,spinTemplate;
  REAL8Array               *PTFM[LAL_NUM_IFO+1];
  REAL8Array               *PTFN[LAL_NUM_IFO+1];
  COMPLEX8VectorSequence   *PTFqVec[LAL_NUM_IFO+1];

  /* triggered sky position and sensitivity structures */
  LIGOTimeGPS              segStartTime;
  struct timeval           startTime;
  LALDetector              *detectors[LAL_NUM_IFO+1];
  REAL4                    *timeOffsets;
  REAL4                    *Fplus;
  REAL4                    *Fcross;
  REAL4                    *Fplustrig;
  REAL4                    *Fcrosstrig;
  REAL8                    FplusTmp;
  REAL8                    FcrossTmp;
  REAL8                    detLoc[3];
  REAL4                    *timeSlideVectors;

  /* coherent statistic structures */
  REAL4TimeSeries          *cohSNR                  = NULL;
  REAL4TimeSeries          *pValues[10];
  REAL4TimeSeries          *snrComps[LAL_NUM_IFO];
  REAL4TimeSeries          *gammaBeta[2];
  REAL4TimeSeries          *nullSNR                 = NULL;
  REAL4TimeSeries          *traceSNR                = NULL;

  /* consistency test structures */
  REAL4TimeSeries          *bankVeto[LAL_NUM_IFO+1];
  REAL4TimeSeries          *autoVeto[LAL_NUM_IFO+1];
  REAL4TimeSeries          *chiSquare[LAL_NUM_IFO+1];

  /* output event structures */
  MultiInspiralTable       *eventList               = NULL;
  MultiInspiralTable       *thisEvent               = NULL;
  UINT8                    eventId                  = 0;
  
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

  verbose("Made fft plans %ld \n", timeval_subtract(&startTime));

  /* NULL out pointers where necessary */
  for (i = 0 ; i < 10 ; i++)
  {
    pValues[i] = NULL;
  }   
  for (i = 0 ; i < LAL_NUM_IFO ; i++)
  {
    snrComps[i] = NULL;
    bankVeto[i] = NULL;
    autoVeto[i] = NULL;
    chiSquare[i] = NULL;
  }
  bankVeto[LAL_NUM_IFO] = NULL;
  autoVeto[LAL_NUM_IFO] = NULL;
  chiSquare[LAL_NUM_IFO] = NULL;
  gammaBeta[0] = NULL;
  gammaBeta[1] = NULL;

  /* Initialise some of the input file names */
  if (params->spinBank)
  {
    spinBank = 1;
    strncpy(spinFileName, params->spinBank, sizeof(spinFileName)-1);
  }
  if (params->noSpinBank)
    strncpy(noSpinFileName, params->noSpinBank, sizeof(noSpinFileName)-1);

  if (params->numIFO == 0)
  {
    fprintf(stderr, "You have not specified any detectors to analyse");
    return 1;
  }
  else if (params->numIFO == 1)
  {
    fprintf(stdout, "You have only specified one detector, "
                    "why are you using the coherent code? \n");
  }

  /*------------------------------------------------------------------------*
   * read the data, generate segments and the PSD                           *
   *------------------------------------------------------------------------*/

  timeSlideVectors=LALCalloc(1, (LAL_NUM_IFO+1)*
                                params->numOverlapSegments*sizeof(REAL4));
  memset(timeSlideVectors, 0,
         (LAL_NUM_IFO+1) * params->numOverlapSegments * sizeof(REAL4));

  /* loop over ifos */ 
  for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    /* Initialize some of the structures */
    channel[ifoNumber]  = NULL;
    invspec[ifoNumber]  = NULL;
    segments[ifoNumber] = NULL;
    PTFM[ifoNumber]     = NULL;
    PTFN[ifoNumber]     = NULL;
    PTFqVec[ifoNumber]  = NULL;

    /* if ifo is on: */
    if (params->haveTrig[ifoNumber])
    {
      /* Read in data from the various ifos */
      channel[ifoNumber] = coh_PTF_get_data(params, params->channel[ifoNumber],
                                            params->dataCache[ifoNumber],
                                            ifoNumber);
      coh_PTF_rescale_data(channel[ifoNumber], 1E20);

      /* compute the spectrum */
      invspec[ifoNumber] = coh_PTF_get_invspec(channel[ifoNumber], fwdplan,
                                               revplan, psdplan, params);

      /* create the segments */

      segments[ifoNumber] = coh_PTF_get_segments(channel[ifoNumber],
                                invspec[ifoNumber],fwdplan, ifoNumber,
                                timeSlideVectors, params);
      
      numSegments = segments[ifoNumber]->numSgmnt;

      verbose("Created segments for one ifo %ld \n",
               timeval_subtract(&startTime));
    }
  }

  /* Create a list of time slide ids for each segment and create time slide
     table. */

  TimeSlideVectorList timeSlideList[numSegments];
  UINT4 slideCount = 0;
  INT8  slideIDList[numSegments];
  
  for (i = 0 ; i < numSegments ; i++)
  {
    UINT4 slideDuplicate = 0;
    if (slideCount)
    {
      for (uj = 0; uj < slideCount; uj++)
      {
        UINT4 slideChecking = 1;
        for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
        {
          if (params->haveTrig[ifoNumber])
          {
            if (timeSlideVectors[ifoNumber*params->numOverlapSegments+i] != \
                timeSlideList[uj].timeSlideVectors[ifoNumber])
            {
              slideChecking = 0;
            }
          }
        }
        if (slideChecking)
        {
          slideDuplicate = 1;
          slideIDList[i] = timeSlideList[uj].timeSlideID;
        }
      }
    }
    if (! slideDuplicate)
    {
      for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if (params->haveTrig[ifoNumber])
        {
          timeSlideList[slideCount].timeSlideVectors[ifoNumber] = \
              timeSlideVectors[ifoNumber*params->numOverlapSegments+i];
        }
      }
      timeSlideList[slideCount].timeSlideID = slideCount;
      slideIDList[i] = timeSlideList[slideCount].timeSlideID;
      slideCount++;
    }
  }

  TimeSlide *time_slide_head=NULL;
  TimeSlide *curr_slide = NULL;

  for (ui = 0 ; ui < slideCount; ui++)
  {
    for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
    {
      if (params->haveTrig[ifoNumber])
      {
        if (! time_slide_head)
        {
          time_slide_head=XLALCreateTimeSlide();
          curr_slide= time_slide_head;
        }
        else
        {
          curr_slide->next=XLALCreateTimeSlide();
          curr_slide = curr_slide->next;
        }
        curr_slide->time_slide_id = timeSlideList[ui].timeSlideID;
        /* FIXME */
        XLALReturnIFO(ifoName,ifoNumber);
        strncpy(curr_slide->instrument,ifoName,sizeof(curr_slide->instrument)-1);
        curr_slide->offset = timeSlideList[ui].timeSlideVectors[ifoNumber];
        curr_slide->process_id=0;
      }
    }
  }

  /*------------------------------------------------------------------------*
   * Determine time delays and response functions                           *
   * This is computed for all detectors, even if not being analyzed         *
   *------------------------------------------------------------------------*/

  /* generate sky points array */
  skyPoints = coh_PTF_generate_sky_points(params);
  numSkyPoints = skyPoints->numPoints;

  verbose("Generated necessary sky grid with %d points %ld \n",\
          numSkyPoints, timeval_subtract(&startTime));

  for (sp=0; sp<numSkyPoints; sp++)
  {
    verbose("ra = %f dec = %f\n", skyPoints->data[sp].longitude,
            skyPoints->data[sp].latitude);
  }

  /* allocate memory */ 
  timeOffsets = LALCalloc(1, LAL_NUM_IFO*sizeof(REAL4));
  Fplus       = LALCalloc(1, LAL_NUM_IFO*sizeof(REAL4));
  Fcross      = LALCalloc(1, LAL_NUM_IFO*sizeof(REAL4));
  Fplustrig   = LALCalloc(1, LAL_NUM_IFO*sizeof(REAL4));
  Fcrosstrig  = LALCalloc(1, LAL_NUM_IFO*sizeof(REAL4));

  /* loop over ifos if doing triggered search */
  if ((params->skyLooping != ALL_SKY) &&
       (params->skyLooping != TWO_DET_ALL_SKY))
  {
    for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
    {
      detectors[ifoNumber] = LALCalloc(1, sizeof(*detectors[ifoNumber]));
      XLALReturnDetector(detectors[ifoNumber], ifoNumber);
      /* get location in three dimensions */
      for (i = 0; i < 3; i++)
      {
        detLoc[i] = (double) detectors[ifoNumber]->location[i];
      }
      /* set 'segStartTime' to trigger time */
      segStartTime = params->trigTime;
      /* calculate time offsets */
      timeOffsets[ifoNumber] = (REAL4)
          XLALTimeDelayFromEarthCenter(detLoc, skyPoints->data[0].longitude,
                                       skyPoints->data[0].latitude,
                                       &segStartTime);
      /* calculate response functions for trigger */
      XLALComputeDetAMResponse(&FplusTmp, &FcrossTmp,
                               detectors[ifoNumber]->response,
                               skyPoints->data[0].longitude,
                               skyPoints->data[0].latitude, 0.,
                               XLALGreenwichMeanSiderealTime(&segStartTime));
      Fplustrig[ifoNumber] = (REAL4) FplusTmp;
      Fcrosstrig[ifoNumber] = (REAL4) FcrossTmp;
    }
  }

  numPoints = floor(params->segmentDuration * params->sampleRate + 0.5);

  struct bankDataOverlaps *chisqOverlaps = NULL;
  struct bankDataOverlaps *chisqSnglOverlaps = NULL;
  REAL4                   *frequencyRangesPlus[LAL_NUM_IFO+1];
  REAL4                   *frequencyRangesCross[LAL_NUM_IFO+1];

  for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO+1; ifoNumber++)
  {
    frequencyRangesPlus[ifoNumber] = NULL;
    frequencyRangesCross[ifoNumber] = NULL;
  }

  /* Initialize some of the structures */
  ifoNumber           = LAL_NUM_IFO;
  channel[ifoNumber]  = NULL;
  invspec[ifoNumber]  = NULL;
  segments[ifoNumber] = NULL; 
  PTFM[ifoNumber]     = NULL;
  PTFN[ifoNumber]     = NULL;
  PTFqVec[ifoNumber]  = NULL;


  /*------------------------------------------------------------------------*
   * Construct the null stream, its segments and its PSD                    *
   *------------------------------------------------------------------------*/

  if (params->doNullStream)
  {
    /* Read in data from the various ifos */
    if (coh_PTF_get_null_stream(params, channel, Fplustrig, Fcrosstrig,
                                  timeOffsets))
    {
      fprintf(stderr,"Null stream construction failure\n");
      return 1;
    }

    /* compute the spectrum */
    invspec[ifoNumber] = coh_PTF_get_invspec(channel[ifoNumber], fwdplan,\
                                              revplan, psdplan, params);
    /* If white spectrum need to scale this. FIX ME!!! */
    if (params->whiteSpectrum)
    {
      for(ui=0 ; ui < invspec[ifoNumber]->data->length; ui++)
      {
        invspec[ifoNumber]->data->data[ui] *= pow(1./0.3403324,2);
      }
    }

    /* create the segments */
    segments[ifoNumber] = coh_PTF_get_segments(channel[ifoNumber],
                                               invspec[ifoNumber],fwdplan,
                                               ifoNumber, timeSlideVectors,
                                               params);

    numSegments = segments[ifoNumber]->numSgmnt;

    verbose("Created segments for null stream at %ld \n",
            timeval_subtract(&startTime));
    if (params->approximant == FindChirpPTF)
    {
      PTFM[ifoNumber]    = XLALCreateREAL8ArrayL(2, 5, 5);
      PTFN[ifoNumber]    = XLALCreateREAL8ArrayL(2, 5, 5);
      PTFqVec[ifoNumber] = XLALCreateCOMPLEX8VectorSequence (5, numPoints);
    }
    else
    {
      PTFM[ifoNumber]    = XLALCreateREAL8ArrayL(2, 1, 1);
      PTFN[ifoNumber]    = XLALCreateREAL8ArrayL(2, 1, 1);
      PTFqVec[ifoNumber] = XLALCreateCOMPLEX8VectorSequence (1, numPoints);
    }
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
   * Create the relevant structures that will be needed                     *
   *------------------------------------------------------------------------*/

  /* finchirp parameters */
  fcInitParams               = LALCalloc(1, sizeof(*fcInitParams));
  fcTmplt                    = LALCalloc(1, sizeof(*fcTmplt));
  fcTmpltParams              = LALCalloc(1, sizeof(*fcTmpltParams));
  fcTmpltParams->approximant = params->approximant;
  fcTmpltParams->order       = params->order;

  /* Note that although non-spinning only uses Q1, the PTF
  generator still generates Q1-5, thus size of these vectors */
  if (params->approximant == FindChirpPTF)
  {
    fcTmplt->PTFQtilde          = 
        XLALCreateCOMPLEX8VectorSequence(5, numPoints / 2 + 1); 
    fcTmplt->PTFQ         = XLALCreateVectorSequence(5, numPoints);
    fcTmpltParams->PTFphi       = XLALCreateVector(numPoints);
    fcTmpltParams->PTFomega_2_3 = XLALCreateVector(numPoints);
    fcTmpltParams->PTFe1        = XLALCreateVectorSequence(3, numPoints);
    fcTmpltParams->PTFe2        = XLALCreateVectorSequence(3, numPoints);
  }
  else if (params->approximant == FindChirpSP)
  {
    fcTmplt->PTFQtilde          =
        XLALCreateCOMPLEX8VectorSequence(1, numPoints / 2 + 1);
    fcTmplt->data               = XLALCreateCOMPLEX8Vector(numPoints / 2 + 1);
    fcTmpltParams->xfacVec      = XLALCreateVector(numPoints / 2 + 1);
    fcTmpltParams->PTFphi       = XLALCreateVector(numPoints / 2 + 1);
    /* Set the values of xfacVec  This is k^(-1/3) 
       also PTFphi which is k^(-7/6). As these are expensive to compute it
       is cheaper to do it once rather than many times. */
    const REAL4                   xfacExponent = -1.0/3.0;
    REAL4                        *xfac = NULL;
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
    fcTmplt->PTFQtilde          =
        XLALCreateCOMPLEX8VectorSequence(1, numPoints / 2 + 1);
    fcTmplt->data               = XLALCreateCOMPLEX8Vector(numPoints / 2 + 1);
    fcTmpltParams->xfacVec      = XLALCreateVector(numPoints);
  }


  fcTmpltParams->fwdPlan      = XLALCreateForwardREAL4FFTPlan(numPoints, 
                                    params->fftLevel);
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
  // This option holds 2x the length of data that is junk in each segment
  // because of conditioning and the PSD.
  fcTmpltParams->invSpecTrunc = params->truncateDuration;

  for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if (params->haveTrig[ifoNumber])
    {
      if (params->approximant == FindChirpPTF)
      {
        PTFM[ifoNumber]    = XLALCreateREAL8ArrayL(2, 5, 5);
        PTFN[ifoNumber]    = XLALCreateREAL8ArrayL(2, 5, 5);
        PTFqVec[ifoNumber] = XLALCreateCOMPLEX8VectorSequence (5, numPoints);
      }
      else
      {
        PTFM[ifoNumber]    = XLALCreateREAL8ArrayL(2, 1, 1);
        PTFN[ifoNumber]    = XLALCreateREAL8ArrayL(2, 1, 1);
        PTFqVec[ifoNumber] = XLALCreateCOMPLEX8VectorSequence (1, numPoints);
      }
    }
  }

  /* Create an inverse FFT plan */
  invPlan = XLALCreateReverseCOMPLEX8FFTPlan(numPoints, params->fftLevel);

  /*------------------------------------------------------------------------*
   * Read in the tmpltbank xml files                                        *
   *------------------------------------------------------------------------*/

  /* read spinning bank */
  if (params->spinBank)
  {
    numSpinTmplts = InspiralTmpltBankFromLIGOLw(&PTFSpinTemplate,
      spinFileName,startTemplate, stopTemplate);
    if (numSpinTmplts != 0)
    {
      PTFtemplate = PTFSpinTemplate;
      numTmplts = numSpinTmplts;
    }
    else
      params->spinBank = NULL;
      spinBank = 0;
  }

  /* read non-spinning bank */
  if (params->noSpinBank)
  {
    numNoSpinTmplts = InspiralTmpltBankFromLIGOLw(&PTFNoSpinTemplate,
      noSpinFileName,startTemplate, stopTemplate);
    if (numNoSpinTmplts != 0)
    {
      PTFtemplate = PTFNoSpinTemplate;
      numTmplts = numNoSpinTmplts;
    }
    else
      params->noSpinBank = NULL;
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

  /*------------------------------------------------------------------------*
   * initialise bank veto
   *------------------------------------------------------------------------*/

  /* Create the templates needed for the bank veto, if necessary */
  UINT4 subBankSize = 0;
  struct bankTemplateOverlaps *bankNormOverlaps = NULL;
  struct bankComplexTemplateOverlaps *bankOverlaps = NULL;
  struct bankDataOverlaps *dataOverlaps = NULL;

  if (params->doBankVeto)
  {
    verbose("Initializing bank veto filters at %ld \n", timeval_subtract(&startTime));
    /* Reads in and initializes the bank veto sub bank */
    subBankSize = coh_PTF_read_sub_bank(params,&PTFBankTemplates);
    params->BVsubBankSize = subBankSize;
    bankNormOverlaps = LALCalloc(subBankSize,sizeof(*bankNormOverlaps));
    bankOverlaps = LALCalloc(subBankSize,sizeof(*bankOverlaps));
    dataOverlaps = LALCalloc(subBankSize,sizeof(*dataOverlaps));
    bankFcTmplts = LALCalloc(subBankSize, sizeof(*bankFcTmplts));
    /* Create necessary structure to hold Q(f) */
    for (ui =0 ; ui < subBankSize; ui++)
    {
      bankFcTmplts[ui].PTFQtilde = 
          XLALCreateCOMPLEX8VectorSequence(1, numPoints / 2 + 1);
    }
    PTFBankvetoHead = PTFBankTemplates;
    
    for (ui=0 ; ui < subBankSize ; ui++)
    {
      coh_PTF_template(fcTmplt,PTFBankTemplates,
          fcTmpltParams);
      PTFBankTemplates = PTFBankTemplates->next;
      /* Only store Q1. Structures used in fcTmpltParams will be overwritten */
      for (uj = 0 ; uj < (numPoints/2 +1) ; uj++)
      {
        if (params->approximant == FindChirpPTF)
          bankFcTmplts[ui].PTFQtilde->data[uj] = fcTmplt->PTFQtilde->data[uj];
        else
          bankFcTmplts[ui].PTFQtilde->data[uj] = fcTmplt->data->data[uj];
      }
    }
    /* Calculate the overlap between templates for bank veto */
    for (ui = 0 ; ui < subBankSize; ui++)
    {
      for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if (params->haveTrig[ifoNumber])
        {
          bankNormOverlaps[ui].PTFM[ifoNumber]=
              XLALCreateREAL8ArrayL(2, 1, 1);
          memset(bankNormOverlaps[ui].PTFM[ifoNumber]->data,
              0, 1 * sizeof(REAL8));
          /* This function calculates the overlaps between templates */
          /* This returns a REAL4 as the overlap between identical templates*/
          /* must be real. */
          coh_PTF_template_overlaps(params,&(bankFcTmplts[ui]),
              &(bankFcTmplts[ui]),invspec[ifoNumber],0,
              bankNormOverlaps[ui].PTFM[ifoNumber]);
        }
      }
    }
 
    verbose("Generated bank veto filters at %ld \n", timeval_subtract(&startTime));
        
  }

  /*------------------------------------------------------------------------*
   * initialise auto veto
   *------------------------------------------------------------------------*/

  /* Create the structures needed for the auto veto, if necessary */

  UINT4 timeStepPoints = 0;
  struct bankComplexTemplateOverlaps *autoTempOverlaps = NULL;

  if (params->doAutoVeto)
  {
    verbose("Initializing auto veto filters at %ld \n",
            timeval_subtract(&startTime));
    /* Initializations */
    autoTempOverlaps = LALCalloc(params->numAutoPoints,
        sizeof(*autoTempOverlaps));
    timeStepPoints = params->autoVetoTimeStep*params->sampleRate;
    for (uj = 0; uj < params->numAutoPoints; uj++)
    {
      for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if (params->haveTrig[ifoNumber])
        {
          /* If it will be used initialize and zero out the overlap structure*/
          autoTempOverlaps[uj].PTFM[ifoNumber] = XLALCreateCOMPLEX8ArrayL(2, 1,
                                                                          1);
          memset(autoTempOverlaps[uj].PTFM[ifoNumber]->data, 0,
                 1*sizeof(COMPLEX8));
        }
        else
          autoTempOverlaps[uj].PTFM[ifoNumber] = NULL;
      }
    }
    verbose("Generated auto veto filters at %ld \n", 
            timeval_subtract(&startTime));
  }
   
  /*------------------------------------------------------------------------*
   * find gravitational waves
   *------------------------------------------------------------------------*/

  PTFbankhead = PTFtemplate;

  /* loop over segments */
  for (j = 0; j < numSegments; ++j)
  {
    if (params->doBankVeto)
    {
      /* Calculate overlap between bank templates and data for bank veto */
      /* loop over bank veto templates */
      for (ui = 0 ; ui < subBankSize ; ui++)
      {
        for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
        {
          if (params->haveTrig[ifoNumber])
          {
            dataOverlaps[ui].PTFqVec[ifoNumber] =
                XLALCreateCOMPLEX8VectorSequence(1, 3*numPoints/4-numPoints/4+
                                                    10000);
            bankOverlaps[ui].PTFM[ifoNumber]=XLALCreateCOMPLEX8ArrayL(2,1,1);
            /* This function calculates the overlap */
            coh_PTF_bank_filters(params, &(bankFcTmplts[ui]), 0,
                                 &segments[ifoNumber]->sgmnt[j], invPlan,
                                 PTFqVec[ifoNumber],
                                 dataOverlaps[ui].PTFqVec[ifoNumber], 0, 0);
          }
        }
      }
      verbose("Generated bank veto filters for segment %d at %ld \n", j,
              timeval_subtract(&startTime));
    }
    PTFtemplate = PTFbankhead;

    /* Loop over templates in the bank */
    for (i = 0; (i < numTmplts); PTFtemplate = PTFtemplate->next, i++)
    {
      for (ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if (params->haveTrig[ifoNumber])
        {
          segStartTime = segments[ifoNumber]->sgmnt[j].epoch;
          break;
        }
      }

      /* We only analyse middle half so add duration/4 to epoch */
      XLALGPSAdd(&segStartTime, params->segmentDuration/4.0);

      /* If running injections, check whether to analyse */
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

      /* Generate the template */
      /* PTF generator called here. For non spin Q1-5 is generated and stored
         only Q1 will be used. This would need some alteration if we planned to
         use other templates. */
       
      coh_PTF_template(fcTmplt,PTFtemplate,fcTmpltParams);

      if (params->approximant != FindChirpPTF)
      {
        for (uj = 0 ; uj < (numPoints/2 +1) ; uj++)
        {
          fcTmplt->PTFQtilde->data[uj] = fcTmplt->data->data[uj];
        }
      }

      if (spinTemplate)
        verbose("Generated spin template %d at %ld \n", i,
                timeval_subtract(&startTime));
      else
        verbose("Generated no spin template %d at %ld \n", i,
                timeval_subtract(&startTime));

      /* Generate the various time series as needed*/
      /* Need to zero these out */
      /* We only want to store data from middle half of segment */
      cohSNR = XLALCreateREAL4TimeSeries("cohSNR", &segStartTime,
                                         PTFtemplate->fLower,
                                         (1.0/params->sampleRate),
                                         &lalDimensionlessUnit,
                                         3*numPoints/4 - numPoints/4);
      if (params->doNullStream)
        nullSNR = XLALCreateREAL4TimeSeries("nullSNR", &segStartTime,
                                            PTFtemplate->fLower,
                                            (1.0/params->sampleRate),
                                            &lalDimensionlessUnit,
                                            3*numPoints/4 - numPoints/4);
      if (params->doTraceSNR)
        traceSNR = XLALCreateREAL4TimeSeries("traceSNR", &segStartTime,
                                             PTFtemplate->fLower,
                                             (1.0/params->sampleRate),
                                             &lalDimensionlessUnit,
                                             3*numPoints/4 - numPoints/4);
      if (params->doBankVeto)
      {
        if (params->numIFO != 1)
        {
          bankVeto[LAL_NUM_IFO] = XLALCreateREAL4TimeSeries("bank_veto",
                                             &segStartTime,
                                             PTFtemplate->fLower,
                                             (1.0/params->sampleRate),
                                             &lalDimensionlessUnit,
                                             3*numPoints/4 - numPoints/4);
        }
        if (params->doSnglChiSquared)
        {
          for (ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
          {
            if (params->haveTrig[ifoNumber])
            {
              XLALReturnIFO(ifoName,ifoNumber);
              snprintf( name, sizeof( name ), "%s_bank_veto",ifoName);
              bankVeto[ifoNumber] = XLALCreateREAL4TimeSeries(name,
                  &segStartTime,PTFtemplate->fLower,(1.0/params->sampleRate),
                  &lalDimensionlessUnit,3*numPoints/4 - numPoints/4);
            }
          }
        }
      }
      if (params->doAutoVeto)
      {
        if (params->numIFO != 1)
        {
          autoVeto[LAL_NUM_IFO] = XLALCreateREAL4TimeSeries("auto_veto",
                                             &segStartTime,
                                             PTFtemplate->fLower,
                                             (1.0/params->sampleRate),
                                             &lalDimensionlessUnit,
                                             3*numPoints/4 - numPoints/4);
        }
        if (params->doSnglChiSquared)
        {
          for (ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
          {
            if (params->haveTrig[ifoNumber])
            {
              XLALReturnIFO(ifoName,ifoNumber);
              snprintf( name, sizeof( name ), "%s_auto_veto",ifoName);
              autoVeto[ifoNumber] = XLALCreateREAL4TimeSeries(name,
                  &segStartTime,PTFtemplate->fLower,(1.0/params->sampleRate),
                  &lalDimensionlessUnit,3*numPoints/4 - numPoints/4);
            }
          }
        }

      }
      if (params->doChiSquare)
      {
        if (params->numIFO != 1)
        {
          chiSquare[LAL_NUM_IFO] = XLALCreateREAL4TimeSeries("chi_square",
                                                &segStartTime,
                                                PTFtemplate->fLower,
                                                (1.0/params->sampleRate),
                                                &lalDimensionlessUnit,
                                                3*numPoints/4 - numPoints/4);
        }
        if (params->doSnglChiSquared)
        {
          for (ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
          {
            if (params->haveTrig[ifoNumber])
            {
              XLALReturnIFO(ifoName,ifoNumber);
              snprintf( name, sizeof( name ), "%s_chi_square",ifoName);
              chiSquare[ifoNumber] = XLALCreateREAL4TimeSeries(name,
                  &segStartTime,PTFtemplate->fLower,(1.0/params->sampleRate),
                  &lalDimensionlessUnit,3*numPoints/4 - numPoints/4);
            }
          }
        }
      }

      /* Loop over ifos */
      for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if (params->haveTrig[ifoNumber])
        {
          /* Zero the storage vectors for the PTF filters */
          if (params->approximant == FindChirpPTF)
          {
            memset(PTFM[ifoNumber]->data, 0, 25 * sizeof(REAL8));
            memset(PTFqVec[ifoNumber]->data, 0, 
                  5 * numPoints * sizeof(COMPLEX8));
          }
          else
          {
            memset(PTFM[ifoNumber]->data, 0, 1 * sizeof(REAL8));
            memset(PTFqVec[ifoNumber]->data, 0,
                  1 * numPoints * sizeof(COMPLEX8));
          }

          /* Here (h|s) and (h|h) are calculated */
          coh_PTF_normalize(params, fcTmplt, invspec[ifoNumber],
                            PTFM[ifoNumber], NULL, PTFqVec[ifoNumber],
                            &segments[ifoNumber]->sgmnt[j], invPlan,
                            spinTemplate);

          // In this subroutine the overlap between template h and the 
          // bank veto template is calculated
          if (params->doBankVeto)
          {
          for (ui = 0 ; ui < subBankSize ; ui++)
          {
            memset(bankOverlaps[ui].PTFM[ifoNumber]->data,0,1*sizeof(COMPLEX8));
            coh_PTF_complex_template_overlaps(params, &(bankFcTmplts[ui]),
                                              fcTmplt, invspec[ifoNumber], 0,
                                              bankOverlaps[ui].PTFM[ifoNumber]);
          }
          }
          // And if necessary the overlap between h(deltaT) and h
          if (params->doAutoVeto)
          {
            coh_PTF_auto_veto_overlaps(params,fcTmplt,autoTempOverlaps,
                invspec[ifoNumber],invPlan,0,params->numAutoPoints,
                timeStepPoints,ifoNumber);
          }

          verbose("Made filters for ifo %d,segment %d, template %d at %ld \n", 
                   ifoNumber, j, i, timeval_subtract(&startTime));
        }
      }
      /* If necessary calculate the null stream filters */
      if (params->doNullStream)
      {
        if (params-> approximant == FindChirpPTF)
        {
          memset(PTFM[LAL_NUM_IFO]->data, 0, 25*sizeof(REAL8));
          memset(PTFqVec[LAL_NUM_IFO]->data, 0, 5*numPoints*sizeof(COMPLEX8));
        }
        else
        {
          memset(PTFM[LAL_NUM_IFO]->data, 0, 1*sizeof(REAL8));
          memset(PTFqVec[LAL_NUM_IFO]->data, 0, 1*numPoints*sizeof(COMPLEX8));
        }
        coh_PTF_normalize(params, fcTmplt, invspec[LAL_NUM_IFO],
                          PTFM[LAL_NUM_IFO], NULL, PTFqVec[LAL_NUM_IFO],
                          &segments[LAL_NUM_IFO]->sgmnt[j], invPlan,
                          spinTemplate);
        verbose("Made filters for NULL stream,segmen %d, template %d at %ld\n",
                j, i, timeval_subtract(&startTime));
      }
      
      switch(params->skyLooping)
      {
        case SINGLE_SKY_POINT:
        case TWO_DET_SKY_PATCH:
        case SKY_PATCH:
          verbose("Begin loop over sky points at %ld \n",
                  timeval_subtract(&startTime));

          /* set 'segStartTime' to trigger time */
          segStartTime = params->trigTime;

          /* loop over sky points */
          for (sp = 0; sp < numSkyPoints ; sp++)
          {
            for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
            {
              /* get location in three dimensions */
              for (ui = 0; ui < 3; ui++)
              {
                 detLoc[ui] = (double) detectors[ifoNumber]->location[ui];
              }
              /* calculate time offsets */
              timeOffsets[ifoNumber] = (REAL4)
                  XLALTimeDelayFromEarthCenter(detLoc,
                                               skyPoints->data[sp].longitude,
                                               skyPoints->data[sp].latitude,
                                               &segStartTime);
              /* calculate response functions */
              XLALComputeDetAMResponse(&FplusTmp, &FcrossTmp,
                                       detectors[ifoNumber]->response,
                                       skyPoints->data[sp].longitude,
                                       skyPoints->data[sp].latitude,0.,
                                       XLALGreenwichMeanSiderealTime(
                                           &segStartTime));
              Fplus[ifoNumber] = (REAL4) FplusTmp;
              Fcross[ifoNumber] = (REAL4) FcrossTmp;
            }
            if (params->faceAwayAnalysis && params->faceOnAnalysis)
            {
              params->faceOnStatistic = 1;
            }

            // This function calculates the cohSNR time series and all of the
            // signal based vetoes as appropriate
            coh_PTF_statistic(cohSNR, PTFM, PTFqVec, params, spinTemplate,
                              timeOffsets, Fplus, Fcross,
                              j, pValues, gammaBeta, snrComps, nullSNR,
                              traceSNR, bankVeto, autoVeto,
                              chiSquare, subBankSize, bankOverlaps, 
                              bankNormOverlaps, dataOverlaps, autoTempOverlaps,
                              fcTmplt, invspec, segments, invPlan, 
                              &chisqOverlaps,&chisqSnglOverlaps, frequencyRangesPlus,
                              frequencyRangesCross, startTime);
     
            verbose("Made coherent statistic for segment %d, template %d, "
                    "sky point %d at %ld \n", j, i, sp,
                    timeval_subtract(&startTime));      

            /* This function construct triggers from loud events */
            eventId = coh_PTF_add_triggers(params, &eventList, &thisEvent,
                                           cohSNR, *PTFtemplate, eventId,
                                           spinTemplate,
                                           pValues, gammaBeta, snrComps,
                                           nullSNR, traceSNR, bankVeto,
                                           autoVeto, chiSquare, PTFM,
                                           skyPoints->data[sp].longitude,
                                           skyPoints->data[sp].latitude,
                                           slideIDList[j], timeOffsets);

            if (params->faceAwayAnalysis && params->faceOnAnalysis)
            {
              params->faceOnStatistic = 2; 

              coh_PTF_statistic(cohSNR, PTFM, PTFqVec, params, spinTemplate,
                              timeOffsets, Fplus, Fcross,
                              j, pValues, gammaBeta, snrComps, nullSNR,
                              traceSNR, bankVeto, autoVeto,
                              chiSquare, subBankSize, bankOverlaps,
                              bankNormOverlaps, dataOverlaps, autoTempOverlaps,
                              fcTmplt, invspec, segments, invPlan,
                              &chisqOverlaps,&chisqSnglOverlaps, frequencyRangesPlus,
                              frequencyRangesCross, startTime);

              verbose("Made coherent statistic for segment %d, template %d, "
                      "sky point %d at %ld \n", j, i, sp,
                      timeval_subtract(&startTime));

              eventId = coh_PTF_add_triggers(params, &eventList, &thisEvent,
                                           cohSNR, *PTFtemplate, eventId,
                                           spinTemplate,
                                           pValues, gammaBeta, snrComps,
                                           nullSNR, traceSNR, bankVeto,
                                           autoVeto, chiSquare, PTFM,
                                           skyPoints->data[sp].longitude,
                                           skyPoints->data[sp].latitude,
                                           slideIDList[j], timeOffsets);
            }

            params->numEvents = XLALCountMultiInspiral(eventList);
            verbose("There are currently %d triggers.\n", params->numEvents);
            verbose("Generated triggers for segment %d, template %d, sky point %d at %ld \n", j, i, sp, timeval_subtract(&startTime));


          }
          break;
        case ALL_SKY:
        case TWO_DET_ALL_SKY:
          error("All sky code not implemented yet\n"); 
          break;
        default:
          error("Oops. No sky type set, Ian done broke the code.\n");
          break;

      } 

      /* cluster triggers */
      /* XLALClusterMultiInspiralTable(); */
      /* verbose("Clustered triggers for segment %d, template %d at %ld \n",
               j, i, time(NULL)-startTime); */

      /* Then we get a bunch of memory freeing statements */
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
      if (nullSNR) XLALDestroyREAL4TimeSeries(nullSNR);
      if (traceSNR) XLALDestroyREAL4TimeSeries(traceSNR);
      if (cohSNR) XLALDestroyREAL4TimeSeries(cohSNR);

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
    }
    if (params->doBankVeto)
    {
      for (ui = 0 ; ui < subBankSize ; ui++)
      {
        for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
        {
          if (dataOverlaps[ui].PTFqVec[ifoNumber])
            XLALDestroyCOMPLEX8VectorSequence(dataOverlaps[ui].PTFqVec[ifoNumber]);
          if (bankOverlaps[ui].PTFM[ifoNumber])
            XLALDestroyCOMPLEX8Array(bankOverlaps[ui].PTFM[ifoNumber]);
        }
      }
    }
  } // Main loop is ended here
  /* calulate number of events */
  params->numEvents = XLALCountMultiInspiral(eventList);
  fprintf(stderr,"There are %d total triggers before cluster.\n", params->numEvents);
  if ( params->clusterFlag )
  {
    coh_PTF_cluster_triggers(params,&eventList,&thisEvent);
    params->numEvents = XLALCountMultiInspiral(eventList);
    fprintf(stderr,"There are %d total triggers after cluster.\n", params->numEvents);
  }

  coh_PTF_output_events_xml(params->outputFile, eventList, params->injectList,\
                            procpar, time_slide_head, params);

  if (skyPoints->data)
    LALFree(skyPoints->data);
  if (skyPoints)
    LALFree(skyPoints);

  // This function cleans up memory usage
  XLALDestroyTimeSlideTable(time_slide_head);
  LALFree(timeSlideVectors);
  coh_PTF_cleanup(params,procpar,fwdplan,psdplan,revplan,invPlan,channel,
      invspec,segments,eventList,PTFbankhead,fcTmplt,fcTmpltParams,
      fcInitParams,PTFM,PTFN,PTFqVec,timeOffsets,Fplus,Fcross,Fplustrig,Fcrosstrig);
  
  while (PTFBankvetoHead)
  {
    InspiralTemplate *thisTmplt;
    thisTmplt = PTFBankvetoHead;
    PTFBankvetoHead = PTFBankvetoHead->next;
    if (thisTmplt->event_id)
    {
      LALFree(thisTmplt->event_id);
    }
    LALFree(thisTmplt);
  }

  coh_PTF_free_bank_veto_memory(bankNormOverlaps,PTFBankTemplates,bankFcTmplts,subBankSize,bankOverlaps,dataOverlaps);

  if (autoTempOverlaps)
  {
    for (uj = 0; uj < params->numAutoPoints; uj++)
    {
      for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if (params->haveTrig[ifoNumber])
        {
          if (autoTempOverlaps[uj].PTFM[ifoNumber])
          {
            XLALDestroyCOMPLEX8Array(autoTempOverlaps[uj].PTFM[ifoNumber]);
          }
        }
      }
    }
    LALFree(autoTempOverlaps);
  }

  for(ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    LALFree(detectors[ifoNumber]);
  }
  verbose("Generated output xml file, cleaning up and exiting at %ld \n",
      timeval_subtract(&startTime));
  LALCheckMemoryLeaks();
  return 0;
}

void coh_PTF_statistic(
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
    REAL4TimeSeries         *gammaBeta[2],
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
    struct timeval          startTime
)

{
  verbose("Entering the statistic loop at %ld \n",
          timeval_subtract(&startTime));

  /* This function generates the SNR for every point in time and, where
   * appropriate calculates the desired signal based vetoes. */

  CHAR ifoName[LIGOMETA_IFO_MAX];
  char name[LALNameLength];
  UINT4  check;
  UINT4  i, j, k, m, ifoNumber, vecLength, vecLengthTwo;
  INT4   l, timeOffsetPoints[LAL_NUM_IFO];
  REAL4  deltaT    = cohSNR->deltaT;
  UINT4  numPoints = floor(params->segmentDuration*params->sampleRate+0.5);

  struct bankDataOverlaps *chisqOverlaps = *chisqOverlapsP;
  struct bankDataOverlaps *chisqSnglOverlaps = *chisqSnglOverlapsP;

  /* Code works slightly differently if spin/non spin and single/coherent */
  if (spinTemplate)
    vecLength = 5;
  else
    vecLength = 1;
  if (params->numIFO == 1 || params->singlePolFlag || params->faceOnStatistic)
    vecLengthTwo = vecLength;
  else
    vecLengthTwo = 2* vecLength;

  REAL4 cohSNRThreshold = params->threshold;
  if (spinTemplate)
    cohSNRThreshold = params->spinThreshold;

  UINT4 csVecLength = 1;
  UINT4 csVecLengthTwo = 2;
  if (params->numIFO == 1 || params->singlePolFlag || params->faceOnStatistic)
    csVecLengthTwo = 1;

  /* These arrays are used to store the maximized quantities
   * For non spin these are the 4 F-stat parameters (only 2 for one detector)
   * For spin these are the P values */
  for (i = 0 ; i < vecLengthTwo ; i++)
  {
    if (!pValues[i])
    {
      pValues[i] = XLALCreateREAL4TimeSeries("Pvalue", &cohSNR->epoch, 
                                             cohSNR->f0, cohSNR->deltaT,
                                             &lalDimensionlessUnit,
                                             cohSNR->data->length);
    }
  }
  if (! spinTemplate)
  {
    for (i = vecLengthTwo ; i < 2*vecLengthTwo ; i++)
    {
      if (!pValues[i])
      {
        pValues[i] = XLALCreateREAL4TimeSeries("Pvalue", &cohSNR->epoch,
                                               cohSNR->f0, cohSNR->deltaT,
                                               &lalDimensionlessUnit,
                                               cohSNR->data->length);
      }
    }
  }
  if (spinTemplate)
  {
    if (!gammaBeta[0])
    {
      /* These store a amplitude and phase information for PTF search */
      gammaBeta[0] = XLALCreateREAL4TimeSeries("Gamma", &cohSNR->epoch,
                                               cohSNR->f0, cohSNR->deltaT,
                                               &lalDimensionlessUnit,
                                               cohSNR->data->length);
    }
    if (!gammaBeta[1])
    {
      gammaBeta[1] = XLALCreateREAL4TimeSeries("Beta", &cohSNR->epoch,
                                               cohSNR->f0, cohSNR->deltaT,
                                               &lalDimensionlessUnit,
                                               cohSNR->data->length);
    }
  }

  /* FIXME: All the time series should be outputtable */
  REAL4 u1[vecLengthTwo], u2[vecLengthTwo], v1[vecLengthTwo], v2[vecLengthTwo];
  REAL4 *v1p,*v2p,*snglv1p,*snglv2p;
  REAL4 u1N[vecLength], u2N[vecLength], v1N[vecLength], v2N[vecLength];
  REAL4 v1_dot_u1, v1_dot_u2, v2_dot_u1, v2_dot_u2,max_eigen;
  REAL4 recSNR, traceSNRsq;
  REAL4 dAlpha, dBeta, dCee;
  REAL4 pValsTemp[vecLengthTwo];
  REAL4 betaGammaTemp[2];

  gsl_matrix *B2Null = gsl_matrix_alloc(vecLength, vecLength);
  /* FIXME: the 50s below seem to hardcode a limit on the number of templates
     this should not be hardcoded. Note that this value is hardcoded in some
     function declarations as well as here! Double pointers will fix this*/
  gsl_matrix *Bankeigenvecs[50];
  gsl_vector *Bankeigenvals[50];
  gsl_matrix *Autoeigenvecs = NULL;
  gsl_vector *Autoeigenvals = NULL;
  for (i = 0; i < 50; i++)
  {
    Bankeigenvecs[i] = NULL;
    Bankeigenvals[i] = NULL;  
  }

  gsl_eigen_symmv_workspace *matTempNull = gsl_eigen_symmv_alloc(vecLength);
  gsl_matrix                *eigenvecs   = gsl_matrix_alloc(vecLengthTwo,
                                                             vecLengthTwo);
  gsl_vector                *eigenvals   = gsl_vector_alloc(vecLengthTwo);
  gsl_matrix                *eigenvecsNull = gsl_matrix_alloc(vecLength,
                                                              vecLength);
  gsl_vector                *eigenvalsNull = gsl_vector_alloc(vecLength);
  gsl_matrix                *eigenvecsSngl = gsl_matrix_alloc(vecLength,
                                                              vecLength);
  gsl_vector                *eigenvalsSngl = gsl_vector_alloc(vecLength);

  // This function takes the (Q_i|Q_j) matrices, combines it across the ifos
  // and returns the eigenvalues and eigenvectors of this new matrix.
  // We later rotate and rescale the (Q_i|s) values such that in the new basis
  // this matrix will be the identity matrix.
  // For non-spin this describes the rotation into the dominant polarization
  coh_PTF_calculate_bmatrix(params,eigenvecs,eigenvals,Fplus,Fcross,PTFM,vecLength,vecLengthTwo,vecLength);

  // If required also calculate these eigenvalues/vectors for the null stream 
  if (params->doNullStream)
  {
    for (i = 0; i < vecLength; i++)
    {
      for (j = 0; j < vecLength; j++)
      {
        gsl_matrix_set(B2Null, i, j, PTFM[LAL_NUM_IFO]->data[i*5+j]);
      }
    }
    gsl_eigen_symmv(B2Null, eigenvalsNull, eigenvecsNull, matTempNull); 
  }

  /* This loop takes the time offset in seconds and converts to time offset
  * in data points (rounded) */
  verbose("Time offsets (s): ");
  for (i = 0; i < LAL_NUM_IFO; i++)
    if (params->haveTrig[i])
      verbose("%f ", timeOffsets[i]);
  verbose("\n");
  verbose("Time offsets (points): ");
  for (i = 0; i < LAL_NUM_IFO; i++)
  {
    timeOffsetPoints[i] = (int) floor(timeOffsets[i]/deltaT + 0.5);
    if (params->haveTrig[i])
      verbose("%d " , timeOffsetPoints[i]);
  }
  verbose("\n");

  v1p = LALCalloc(vecLengthTwo , sizeof(REAL4));
  v2p = LALCalloc(vecLengthTwo , sizeof(REAL4));
  snglv1p = LALCalloc(vecLength , sizeof(REAL4));
  snglv2p = LALCalloc(vecLength , sizeof(REAL4));

  /* If we have injections we only want to analyse the time around the
   * injection. Here we figure out what that time should be. No injections
   * equates to analyse the whole segment.
   */

  UINT4 segStartPoint = 0;
  UINT4 segEndPoint = 0;

  if ( params->injectFile )
  {
    findInjectionSegment(&segStartPoint, &segEndPoint, &cohSNR->epoch, params);
  }

  if (! segStartPoint)
  {
    segStartPoint = numPoints/4;
    segEndPoint = 3*numPoints/4;
  }

  /* First, calculate the single detector SNR time series. Only do this for
     the first sky point. */
  /* NOT CORRECT FOR SPINNING FIXME! */

  REAL4 reSNRcomp,imSNRcomp;
  UINT4 ifoNum1,ifoNum2;

  ifoNum1 = 0;
  ifoNum2 = 0;

  /* FIXME: For >2 detectors this should choose the 2 most sensitive ifos */
  /* NOTE: This is only used in two det case where the 2 most sensitive ifos
     is obvious! */
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

  for (ifoNumber = 0;ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if (params->haveTrig[ifoNumber] && (! snrComps[ifoNumber]))
    {
      verbose("Calculating single detector SNR for ifo %d at %ld.\n",
              ifoNumber, timeval_subtract(&startTime));
      XLALReturnIFO(ifoName,ifoNumber);
      snprintf( name, sizeof( name ), "%s_snr",ifoName);
      snrComps[ifoNumber] = XLALCreateREAL4TimeSeries(name,
                                                      &cohSNR->epoch,
                                                      cohSNR->f0,
                                                      cohSNR->deltaT,
                                                      &lalDimensionlessUnit,
                                                      3*numPoints/4
                                                          -numPoints/4+10000);

      if (spinTemplate)
      {
        /* convert PTFM to gsl_matrix */
        gsl_matrix_view PTFmatrix;
        PTFmatrix = gsl_matrix_view_array(PTFM[ifoNumber]->data,\
                                          PTFM[ifoNumber]->dimLength->data[0],
                                          PTFM[ifoNumber]->dimLength->data[1]);
        /* calculate eigenvectors and eigenvalues of (h|h)*/
        gsl_eigen_symmv_workspace *matTemp = gsl_eigen_symmv_alloc(5);
        gsl_eigen_symmv(&(PTFmatrix.matrix), eigenvalsSngl, eigenvecsSngl,
                        matTemp);
        gsl_eigen_symmv_free(matTemp);
        for (i = segStartPoint-5000; i < segEndPoint+5000; ++i)
        {  /* loop over time */ 
          // This function combines the various (Q_i | s) and rotates them into
          // the basis as discussed above.
          coh_PTF_calculate_rotated_vectors(params,PTFqVec,snglv1p,snglv2p,Fplus,Fcross,
            timeOffsetPoints,eigenvecsSngl,eigenvalsSngl,
            numPoints,i,vecLength,vecLength,ifoNumber);

          /* Compute the dot products */
          v1_dot_u1 = v1_dot_u2 = v2_dot_u1 = v2_dot_u2 = max_eigen = 0.0;
          for (j = 0; j < vecLength; j++)
          {
            v1_dot_u1 += snglv1p[j] * snglv1p[j];
            v1_dot_u2 += snglv1p[j] * snglv2p[j];
            v2_dot_u2 += snglv2p[j] * snglv2p[j];
          }
          // And SNR is calculated
          // For spin this follows Diego's notation
          max_eigen =0.5 * (v1_dot_u1 + v2_dot_u2 + sqrt((v1_dot_u1 - v2_dot_u2)
              * (v1_dot_u1 - v2_dot_u2) + 4 * v1_dot_u2 * v1_dot_u2));
          snrComps[ifoNumber]->data->data[i-((numPoints/4)-5000)] =
              sqrt(max_eigen);
        }
      }
      else
      {
        for (i = segStartPoint-5000; i < segEndPoint+5000; ++i)
        {  /* loop over time */
          reSNRcomp = PTFqVec[ifoNumber]->data[i].re;
          imSNRcomp = PTFqVec[ifoNumber]->data[i].im;
          snrComps[ifoNumber]->data->data[i - ((numPoints/4)-5000)] =
              sqrt((reSNRcomp*reSNRcomp + imSNRcomp*imSNRcomp)/
                   PTFM[ifoNumber]->data[0]);
        }
      }
    }
  }


  verbose("Begin loop over time at %ld \n",timeval_subtract(&startTime));

  /* These are declared to optimize what comes below */
  REAL4 coincSNR;
  INT4  tOffset0 = 5000 - numPoints/4;
  INT4  tOffset1 = tOffset0 + timeOffsetPoints[ifoNum1];
  INT4  tOffset2 = tOffset0 + timeOffsetPoints[ifoNum2];
  INT4  sOffset  = numPoints/4;
  REAL4 SNRthresh = params->snglSNRThreshold;

  for (i = numPoints/4; i < 3*numPoints/4; ++i) /* Main loop over time */
  {
    cohSNR->data->data[i-sOffset] = 0;
  }

  for (i = segStartPoint; i < segEndPoint; ++i) /* Main loop over time */
  {
    /* Don't bother calculating coherent SNR if all ifo's SNR is less than
       some value */
    for (ifoNumber = 0;ifoNumber < LAL_NUM_IFO; ifoNumber++)
    {
      if (params->haveTrig[ifoNumber])
      {
        if (snrComps[ifoNumber]->data->
                data[i+tOffset0+timeOffsetPoints[ifoNumber]] > SNRthresh)
        {
          break;
        }
      }
    }
    if (ifoNumber == LAL_NUM_IFO)
    {
      cohSNR->data->data[i-sOffset] = 0;
      continue;
    }

    if (params->numIFO == 2 && (! params->singlePolFlag) && (!params->faceOnStatistic) )
    {
      cohSNR->data->data[i-sOffset] = sqrt(snrComps[ifoNum1]->data->
                                               data[i+tOffset1] *
                                           snrComps[ifoNum1]->data->
                                               data[i+tOffset1] +
                                           snrComps[ifoNum2]->data->
                                               data[i+tOffset2] *
                                           snrComps[ifoNum2]->data->
                                               data[i+tOffset2]);
    }
    else
    {
      coincSNR = 0;
      // We do not need to calculate coherent SNR if coinc SNR < threshold
      for (ifoNumber = 0;ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if (params->haveTrig[ifoNumber])
        {
          coincSNR += snrComps[ifoNumber]->data->data[\
              i + tOffset0 + timeOffsetPoints[ifoNumber]]
              * snrComps[ifoNumber]->data->data[\
              i + tOffset0 + timeOffsetPoints[ifoNumber]];
        }
      }
      if (coincSNR < cohSNRThreshold)
      {
        cohSNR->data->data[i-sOffset] = 0;
        continue;
      }

      // This function combines the various (Q_i | s) and rotates them into
      // the basis as discussed above.
      coh_PTF_calculate_rotated_vectors(params,PTFqVec,v1p,v2p,Fplus,Fcross,
        timeOffsetPoints,eigenvecs,eigenvals,
        numPoints,i,vecLength,vecLengthTwo,LAL_NUM_IFO);

      /* Compute the dot products */
      v1_dot_u1 = v1_dot_u2 = v2_dot_u1 = v2_dot_u2 = max_eigen = 0.0;
      for (j = 0; j < vecLengthTwo; j++)
      {
        v1_dot_u1 += v1p[j] * v1p[j];
        v2_dot_u2 += v2p[j] * v2p[j];
      }
      // And SNR is calculated
      // For non spin: v1p[0] * v1p[0] = (\bf{F}_+\bf{h}_0 | \bf{s})^2
      //               v1p[1] * v1p[1] = (\bf{F}_x\bf{h}_0 | \bf{s})^2
      //               v2p[0] * v2p[0] = (\bf{F}_+\bf{h}_{\pi/2} | \bf{s})^2
      //               v2p[1] * v2p[1] = (\bf{F}_x\bf{h}_{\pi/2} | \bf{s})^2
      //
      // For spin this follows Diego's notation
      if (spinTemplate == 0)
      {
        max_eigen = (v1_dot_u1 + v2_dot_u2);
      }
      else
      {
        for (j = 0; j < vecLengthTwo; j++)
        {
          v1_dot_u2 += v1p[j] * v2p[j];
        }
        max_eigen = 0.5 * (v1_dot_u1 + v2_dot_u2 + sqrt((v1_dot_u1 - v2_dot_u2)
            * (v1_dot_u1 - v2_dot_u2) + 4 * v1_dot_u2 * v1_dot_u2));
      }
      cohSNR->data->data[i-sOffset] = sqrt(max_eigen);
    }
  }

  verbose("Calculated all SNRs at %ld \n",timeval_subtract(&startTime));

  UINT4 chisqCheck = 0;
  REAL4 bestNR;
  INT4 numPointCheck = floor(params->timeWindow/cohSNR->deltaT + 0.5);
  struct bankCohTemplateOverlaps *bankCohOverlaps = NULL;
  struct bankCohTemplateOverlaps *autoCohOverlaps = NULL;
  COMPLEX8VectorSequence *tempqVec = NULL;
  REAL4 *powerBinsPlus[LAL_NUM_IFO+1];
  REAL4 *powerBinsCross[LAL_NUM_IFO+1];
  REAL4 fLowPlus,fHighPlus,fLowCross,fHighCross;

  for(k = 0; k < LAL_NUM_IFO+1; k++)
  {
    powerBinsPlus[k] = NULL;
    powerBinsCross[k] = NULL;
  }

  // Now we calculate all the extrinsic parameters and signal based vetoes
  // Only calculated if this will be a trigger

  for (i = segStartPoint; i < segEndPoint; ++i) /* loop over time */
  {
    if (cohSNR->data->data[i-numPoints/4] > cohSNRThreshold)
    {
      check = 1;
      for (l = (INT4)(i-numPoints/4)-numPointCheck; l < (INT4)(i-numPoints/4)+numPointCheck; l++)
      {
        if (l < 0)
          l = 0;
        if (l > (INT4)(cohSNR->data->length-1))
          break;
        if (cohSNR->data->data[l] > cohSNR->data->data[i-numPoints/4])
        {
          check = 0;
          break;
        }
      }
      if (check)
      {
        // The following block extracts the values of extrinsic parameters.
        // This follows the method set out in Diego's thesis to do this.
        /* FIXME: This is probably broke for coherent case now! */
        /* FIXME: First block should be optional! */
        if (0)
        {
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
          if (! spinTemplate)
          {
            /* This is a lot easier when there is no spin 
               Note that we output values in the dominant polarization frame
               for the case of non-spin */
            for (j = 0 ; j < vecLengthTwo ; j++)
            {
              pValues[j]->data->data[i - numPoints/4] = sqrt(v1[j]*u1[j]);
              pValues[j+vecLengthTwo]->data->data[i - numPoints/4] =sqrt(v2[j]*u2[j]);
            }
          }
            
          /* For PTF it is a bit more tricksy. Here we use the methods given
             in Diego's thesis.
             Values are outputted in the original basis */
          if (spinTemplate)
          {
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
            for (j = 0; j < vecLengthTwo ; j++) /* Construct the vi vectors */
            {
              v1[j] = 0.;
              v2[j] = 0.;
              for(k = 0; k < LAL_NUM_IFO; k++)
              {
                if (params->haveTrig[k])
                {
                  if (j < vecLength)
                  {
                    v1[j] += Fplus[k] * PTFqVec[k]->data[j*numPoints+i+timeOffsetPoints[k]].re;
                    v2[j] += Fplus[k] * PTFqVec[k]->data[j*numPoints+i+timeOffsetPoints[k]].im;
                  }
                  else
                  {
                    v1[j] += Fcross[k] * PTFqVec[k]->data[(j-vecLength)*numPoints+i+timeOffsetPoints[k]].re;
                    v2[j] += Fcross[k] * PTFqVec[k]->data[(j-vecLength)*numPoints+i+timeOffsetPoints[k]].im;
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
        }

        // First sbv to be calculated is the null SNR.
        if (params->doNullStream)
        {
          // As with SNR do the rotation and calculate SNR.
          for (j = 0; j < vecLength; j++) /* Construct the vi vectors */
          {
            v1N[j] = PTFqVec[LAL_NUM_IFO]->data[j*numPoints+i].re;
            v2N[j] = PTFqVec[LAL_NUM_IFO]->data[j*numPoints+i].im;
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
          v1_dot_u1 = v1_dot_u2 = v2_dot_u1 = v2_dot_u2 = max_eigen = 0.0;
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
          nullSNR->data->data[i-numPoints/4] = sqrt(max_eigen);
        }

        // Next up is Trace SNR and the SNR components
        /* FIXME: Sngl detector SNRs are only correct for coherent non spin*/
        /* FIXME: This loop should never be run in single detector mode */
        if (params->doTraceSNR)
        {
          // traceSNR is calculated as normal SNR but cross terms are not added
          traceSNRsq = 0;
          for(k = 0; k < LAL_NUM_IFO; k++)
          {
            if (params->haveTrig[k])
            {
              for (j = 0; j < vecLengthTwo ; j++)
              {
                if (j < vecLength)
                {
                  v1[j] = Fplus[k] * PTFqVec[k]->data[j*numPoints+i+timeOffsetPoints[k]].re;
                  v2[j] = Fplus[k] * PTFqVec[k]->data[j*numPoints+i+timeOffsetPoints[k]].im;
                }
                else
                {
                  v1[j] = Fcross[k] * PTFqVec[k]->data[(j-vecLength)*numPoints+i+timeOffsetPoints[k]].re;
                  v2[j] = Fcross[k] * PTFqVec[k]->data[(j-vecLength)*numPoints+i+timeOffsetPoints[k]].im;
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
              v1_dot_u1 = v1_dot_u2 = v2_dot_u1 = v2_dot_u2 = max_eigen = 0.0;
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
              // This needs to be converted for spinning case!
             /* snglSNRsq = v1[0]*v1[0] + v2[0]*v2[0];
              snglSNRsq = snglSNRsq/(a[k]*a[k]*PTFM[k]->data[0]);
              traceSNRsq += max_eigen;
              snrComps[k]->data->data[i-numPoints/4]=sqrt(snglSNRsq);*/
            }
          }
          traceSNR->data->data[i-numPoints/4] = sqrt(traceSNRsq);
        }

        // Next is the bank veto
        if (params->doBankVeto)
        {
          if (params->numIFO != 1)
          {
            for (j = 0 ; j < subBankSize+1 ; j++)
            {
              if (! Bankeigenvecs[j])
              {
                // FIXME: Lots of hardcoded vector lengths under here
                Bankeigenvecs[j] = gsl_matrix_alloc(csVecLengthTwo,
                                                    csVecLengthTwo);
                Bankeigenvals[j] = gsl_vector_alloc(csVecLengthTwo);
                // Here we calculate the eigenvectors for each bank template
                if (j == subBankSize)
                {
                  coh_PTF_calculate_bmatrix(params,Bankeigenvecs[j],
                      Bankeigenvals[j],Fplus,Fcross,PTFM,csVecLength,csVecLengthTwo,
                      vecLength);
                }
                else
                {
                  coh_PTF_calculate_bmatrix(params,Bankeigenvecs[j],
                      Bankeigenvals[j],Fplus,Fcross,bankNormOverlaps[j].PTFM,csVecLength,
                      csVecLengthTwo,vecLength);
                }
              }
            }

            if (! bankCohOverlaps)
            {
              bankCohOverlaps = LALCalloc(subBankSize,sizeof(*bankCohOverlaps));
              for (j = 0 ; j < subBankSize; j++)
              {
                bankCohOverlaps[j].rotReOverlaps = gsl_matrix_alloc(
                    csVecLengthTwo,csVecLengthTwo);
                bankCohOverlaps[j].rotImOverlaps = gsl_matrix_alloc(
                    csVecLengthTwo,csVecLengthTwo);
                // We calculate the coherent overlaps in this function
                coh_PTF_calculate_coherent_bank_overlaps(params,bankOverlaps[j],
                    bankCohOverlaps[j],Fplus,Fcross,Bankeigenvecs[subBankSize],
                    Bankeigenvals[subBankSize],Bankeigenvecs[j],
                    Bankeigenvals[j],csVecLength,csVecLengthTwo);
              }
            }
            // In this function all the filters are combined to produce the
            // value of the bank veto.
            bankVeto[LAL_NUM_IFO]->data->data[i-numPoints/4] = coh_PTF_calculate_bank_veto(numPoints,i,subBankSize,Fplus,Fcross,params,bankCohOverlaps,NULL,dataOverlaps,NULL,PTFqVec,NULL,timeOffsetPoints,Bankeigenvecs,Bankeigenvals,LAL_NUM_IFO,csVecLength,csVecLengthTwo);
          }
          // Now, as well as the coherent bank veto calculated above, we can calculate
          // the single detector bank veto
          if (params->doSnglChiSquared)
          {
            for(k = 0; k < LAL_NUM_IFO; k++)
            {
              if (params->haveTrig[k])
              {
                bankVeto[k]->data->data[i-numPoints/4] = coh_PTF_calculate_bank_veto(numPoints,i,subBankSize,Fplus,Fcross,params,NULL,bankOverlaps,dataOverlaps,bankNormOverlaps,PTFqVec,PTFM,timeOffsetPoints,NULL,NULL,k,1,1);
              }
            }            
          }
        }   

        // Now we do the auto veto
        if (params->doAutoVeto)
        {
          if (params->numIFO!=1)
          {
            if (! Autoeigenvecs)
            {
              Autoeigenvecs = gsl_matrix_alloc(csVecLengthTwo,csVecLengthTwo);
              Autoeigenvals = gsl_vector_alloc(csVecLengthTwo);
              // Again the eigenvectors/values are calculated
              coh_PTF_calculate_bmatrix(params,Autoeigenvecs,Autoeigenvals,
                  Fplus,Fcross,PTFM,csVecLength,csVecLengthTwo,vecLength);
            }

            if (! autoCohOverlaps)
            {
              autoCohOverlaps = LALCalloc(params->numAutoPoints,sizeof(*autoCohOverlaps));
              for (j = 0 ; j < params->numAutoPoints; j++)
              {
                autoCohOverlaps[j].rotReOverlaps = gsl_matrix_alloc(
                    csVecLengthTwo,csVecLengthTwo);
                autoCohOverlaps[j].rotImOverlaps = gsl_matrix_alloc(
                    csVecLengthTwo,csVecLengthTwo);
                // The coherent rotated overlaps are calculated
                coh_PTF_calculate_coherent_bank_overlaps(
                    params,autoTempOverlaps[j],
                    autoCohOverlaps[j],Fplus,Fcross,Autoeigenvecs,Autoeigenvals,
                    Autoeigenvecs,Autoeigenvals,csVecLength,csVecLengthTwo);
              }
            }
            // Auto veto is calculated
            autoVeto[LAL_NUM_IFO]->data->data[i-numPoints/4] = coh_PTF_calculate_auto_veto(numPoints,i,Fplus,Fcross,params,autoCohOverlaps,NULL,PTFqVec,NULL,timeOffsetPoints,Autoeigenvecs,Autoeigenvals,LAL_NUM_IFO,csVecLength,csVecLengthTwo);
          }
          if (params->doSnglChiSquared)
          {
            for(k = 0; k < LAL_NUM_IFO; k++)
            {
              if (params->haveTrig[k])
              {
                autoVeto[k]->data->data[i-numPoints/4] = coh_PTF_calculate_auto_veto(numPoints,i,Fplus,Fcross,params,NULL,autoTempOverlaps,PTFqVec,PTFM,timeOffsetPoints,NULL,NULL,k,1,1);
              }
            }
          }

        }
      }
    }
  }
  /* To save memory we cut the loop here, clean the memory before calculating
  chi square */
  LALFree(v1p);
  LALFree(v2p);
  LALFree(snglv1p);
  LALFree(snglv2p);

  verbose("Calculated most vetoes at %ld \n",timeval_subtract(&startTime));

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

  /* And do the loop again to calculate chi square */
  REAL4 numDOF = 4.;
  if (params->singlePolFlag)
    numDOF = 2.;

  for (i = segStartPoint; i < segEndPoint; ++i) /* loop over time */
  {
    if (cohSNR->data->data[i-numPoints/4] > cohSNRThreshold)
    {
      check = 1;
      for (l = (INT4)(i-numPoints/4)-numPointCheck; l < (INT4)(i-numPoints/4)+numPointCheck; l++)
      {
        if (l < 0)
          l = 0;
        if (l > (INT4)(cohSNR->data->length-1))
          break;
        if (cohSNR->data->data[l] > cohSNR->data->data[i-numPoints/4])
        {
          check = 0;
          break;
        }
      }
      if (check)
      {
        /* Test whether to do chi^2 */
        if (params->chiSquareCalcThreshold)
        {
          chisqCheck = 1;
          
          bestNR = cohSNR->data->data[i-numPoints/4];

          /* IS the null stream too large? */
          if (params->doNullStream)
          {
            if (nullSNR->data->data[i-numPoints/4] > params->nullStatThreshold \
                && bestNR < params->nullStatGradOn)
            {
              chisqCheck = 0;
            }
            else if (bestNR > params->nullStatGradOn)
            {
              if (nullSNR->data->data[i-numPoints/4] > (params->nullStatThreshold + (bestNR - params->nullStatGradOn)*params->nullStatGradient))
              {
                chisqCheck = 0;
              }
            }
          }
  
          /* Is bank new SNR too large? */
          if (params->doBankVeto)
          {
            if (bankVeto[LAL_NUM_IFO]->data->data[i-numPoints/4] > 
                                              subBankSize*numDOF)
              bestNR = bestNR/pow((1 + pow(bankVeto[LAL_NUM_IFO]->data->data[i-numPoints/4]/((REAL4)subBankSize*numDOF),params->bankVetoq/params->bankVeton))/2.,1./params->bankVetoq);
            if (bestNR < params->chiSquareCalcThreshold)
              chisqCheck = 0;
          }

          bestNR = cohSNR->data->data[i-numPoints/4];

          /* Is auto new SNR too large */
          if (params->doAutoVeto)
          {
            if (autoVeto[LAL_NUM_IFO]->data->data[i-numPoints/4] > 
                                              params->numAutoPoints*numDOF)
              bestNR = bestNR/pow((1 + pow(autoVeto[LAL_NUM_IFO]->data->data[i-numPoints/4]/((REAL4)params->numAutoPoints*numDOF),params->autoVetoq/params->autoVeton))/2.,1./params->autoVetoq);
            if (bestNR < params->chiSquareCalcThreshold)
              chisqCheck = 0;
          } 
        }
        else
          chisqCheck = 1;

        /* If no problems then calculate chi squared */
        if (params->doChiSquare && params->numIFO != 1)
        {
          if (chisqCheck)
          {
            if (! Autoeigenvecs)
            {
              Autoeigenvecs = gsl_matrix_alloc(csVecLengthTwo,csVecLengthTwo);
              Autoeigenvals = gsl_vector_alloc(csVecLengthTwo);
              // Again the eigenvectors/values are calculated
              coh_PTF_calculate_bmatrix(params,Autoeigenvecs,Autoeigenvals,
                  Fplus,Fcross,PTFM,csVecLength,csVecLengthTwo,vecLength);
            }
            if (! frequencyRangesPlus[LAL_NUM_IFO])
            {
              frequencyRangesPlus[LAL_NUM_IFO] = (REAL4 *)
                LALCalloc(params->numChiSquareBins-1, sizeof(REAL4));
              frequencyRangesCross[LAL_NUM_IFO] = (REAL4 *)
                LALCalloc(params->numChiSquareBins-1, sizeof(REAL4));
              coh_PTF_calculate_standard_chisq_freq_ranges(params,fcTmplt,invspec,PTFM,Fplus,Fcross,frequencyRangesPlus[LAL_NUM_IFO],frequencyRangesCross[LAL_NUM_IFO],Autoeigenvecs,LAL_NUM_IFO,params->singlePolFlag);
            }
            if (! powerBinsPlus[LAL_NUM_IFO])
            {
              powerBinsPlus[LAL_NUM_IFO] = (REAL4 *) 
                LALCalloc(params->numChiSquareBins, sizeof(REAL4));
              powerBinsCross[LAL_NUM_IFO] = (REAL4 *)
                LALCalloc(params->numChiSquareBins, sizeof(REAL4));
              coh_PTF_calculate_standard_chisq_power_bins(params,fcTmplt,invspec,PTFM,Fplus,Fcross,frequencyRangesPlus[LAL_NUM_IFO],frequencyRangesCross[LAL_NUM_IFO],powerBinsPlus[LAL_NUM_IFO],powerBinsCross[LAL_NUM_IFO],Autoeigenvecs,LAL_NUM_IFO,params->singlePolFlag);
            }
            if (! tempqVec)
              tempqVec = XLALCreateCOMPLEX8VectorSequence (1, numPoints);
            if (! chisqOverlaps)
            {
              chisqOverlaps = LALCalloc(2*params->numChiSquareBins,sizeof(*chisqOverlaps));
              for(j = 0; j < params->numChiSquareBins; j++)
              {
                if (params->numChiSquareBins == 1)
                {
                  fLowPlus = 0;
                  fHighPlus = 0;
                  fLowCross = 0;
                  fHighCross = 0;
                }
                else if (j == 0)
                {
                  fLowPlus = 0;
                  fHighPlus = frequencyRangesPlus[LAL_NUM_IFO][0];
                  fLowCross = 0;
                  fHighCross = frequencyRangesCross[LAL_NUM_IFO][0];
                }
                else if (j == params->numChiSquareBins-1)
                {
                  fLowPlus = frequencyRangesPlus[LAL_NUM_IFO][params->numChiSquareBins-2];
                  fHighPlus = 0;
                  fLowCross = frequencyRangesCross[LAL_NUM_IFO][params->numChiSquareBins-2];
                  fHighCross = 0;
                }
                else
                {
                  fLowPlus = frequencyRangesPlus[LAL_NUM_IFO][j-1];
                  fHighPlus = frequencyRangesPlus[LAL_NUM_IFO][j];
                  fLowCross = frequencyRangesCross[LAL_NUM_IFO][j-1];
                  fHighCross = frequencyRangesCross[LAL_NUM_IFO][j];
                }                 
                for(k = 0; k < LAL_NUM_IFO; k++)
                {
                  if (params->haveTrig[k])
                  {
                    chisqOverlaps[j].PTFqVec[k] =
                        XLALCreateCOMPLEX8VectorSequence (1,
                        3*numPoints/4 - numPoints/4 + 10000);
                    coh_PTF_bank_filters(params,fcTmplt,0,
                    &segments[k]->sgmnt[segmentNumber],invPlan,tempqVec,
                    chisqOverlaps[j].PTFqVec[k],fLowPlus,fHighPlus);

                    chisqOverlaps[j+params->numChiSquareBins].PTFqVec[k] =
                        XLALCreateCOMPLEX8VectorSequence (1,
                        3*numPoints/4 - numPoints/4 + 10000);
                    coh_PTF_bank_filters(params,fcTmplt,0,
                    &segments[k]->sgmnt[segmentNumber],invPlan,tempqVec,
                    chisqOverlaps[j+params->numChiSquareBins].PTFqVec[k],
                    fLowCross,fHighCross);
                  }
                  else
                  {
                    chisqOverlaps[j].PTFqVec[k] = NULL;
                    chisqOverlaps[j+params->numChiSquareBins].PTFqVec[k] = NULL;
                  }
                }
              }
            }
            /* Calculate chi square here */
            chiSquare[LAL_NUM_IFO]->data->data[i-numPoints/4] = coh_PTF_calculate_chi_square(params,numPoints,i,chisqOverlaps,PTFqVec,NULL,Fplus,Fcross,timeOffsetPoints,Autoeigenvecs,Autoeigenvals,powerBinsPlus[LAL_NUM_IFO],powerBinsCross[LAL_NUM_IFO],LAL_NUM_IFO,csVecLength,csVecLengthTwo);
          }
          else if (params->doChiSquare)
            chiSquare[LAL_NUM_IFO]->data->data[i-numPoints/4] = 0;
        }
        if (params->doChiSquare && params->doSnglChiSquared)
        {
          //FIXME: Add a check on calculate single detector chi squared!!!!
          if (! chisqSnglOverlaps)
          {
            chisqSnglOverlaps = LALCalloc(params->numChiSquareBins,sizeof(*chisqSnglOverlaps));
            for (k = 0; k < LAL_NUM_IFO; k++)
            {
              for(j = 0; j < params->numChiSquareBins; j++)
              {
                chisqSnglOverlaps[j].PTFqVec[k] = NULL;
              }
            }
          }
          for(k = 0; k < LAL_NUM_IFO; k++)
          {
            if (params->haveTrig[k])
            {
              if (! frequencyRangesPlus[k])
              {
                frequencyRangesPlus[k] = (REAL4 *)
                    LALCalloc(params->numChiSquareBins-1, sizeof(REAL4));
                frequencyRangesCross[k] = (REAL4 *)
                    LALCalloc(params->numChiSquareBins-1, sizeof(REAL4));
                coh_PTF_calculate_standard_chisq_freq_ranges(params,fcTmplt,invspec,PTFM,Fplus,Fcross,frequencyRangesPlus[k],frequencyRangesCross[k],NULL,k,0);
              }
              if (! powerBinsPlus[k])
              {
                powerBinsPlus[k] = (REAL4 *) 
                  LALCalloc(params->numChiSquareBins, sizeof(REAL4));
                powerBinsCross[k] = (REAL4 *)
                  LALCalloc(params->numChiSquareBins, sizeof(REAL4));
                coh_PTF_calculate_standard_chisq_power_bins(params,fcTmplt,invspec,PTFM,Fplus,Fcross,frequencyRangesPlus[k],frequencyRangesCross[k],powerBinsPlus[k],powerBinsCross[k],NULL,k,0);
              }
              if (! tempqVec)
                tempqVec = XLALCreateCOMPLEX8VectorSequence (1, numPoints);
             
              for(j = 0; j < params->numChiSquareBins; j++)
              {
                if (! chisqSnglOverlaps[j].PTFqVec[k])
                {
                  if (params->numChiSquareBins == 1)
                  {
                    fLowPlus = 0;
                    fHighPlus = 0;
                    fLowCross = 0;
                    fHighCross = 0;
                  }
                  else if (j == 0)
                  {
                    fLowPlus = 0;
                    fHighPlus = frequencyRangesPlus[k][0];
                    fLowCross = 0;
                    fHighCross = frequencyRangesCross[k][0];
                  }
                  else if (j == params->numChiSquareBins-1)
                  {
                    fLowPlus = frequencyRangesPlus[k][params->numChiSquareBins-2];
                    fHighPlus = 0;
                    fLowCross = frequencyRangesCross[k][params->numChiSquareBins-2];
                    fHighCross = 0;
                  }
                  else
                  {
                    fLowPlus = frequencyRangesPlus[k][j-1];
                    fHighPlus = frequencyRangesPlus[k][j];
                    fLowCross = frequencyRangesCross[k][j-1];
                    fHighCross = frequencyRangesCross[k][j];
                  }
                  // Set the overlaps, for single detector plus and cross must
                  // be the same
                  // FIXME: Do not need to calculate cross filters above!
                  chisqSnglOverlaps[j].PTFqVec[k] =
                        XLALCreateCOMPLEX8VectorSequence (1,
                        3*numPoints/4 - numPoints/4 + 10000);
                  coh_PTF_bank_filters(params,fcTmplt,0,
                  &segments[k]->sgmnt[segmentNumber],invPlan,tempqVec,
                  chisqSnglOverlaps[j].PTFqVec[k],fLowPlus,fHighPlus);
                }
              }
              // Calculate chi squared
              chiSquare[k]->data->data[i-numPoints/4] = coh_PTF_calculate_chi_square(params,numPoints,i,chisqSnglOverlaps,PTFqVec,PTFM,Fplus,Fcross,timeOffsetPoints,NULL,NULL,powerBinsPlus[k],powerBinsCross[k],k,1,1);
 
            }
          }
        }
              

      }
    }
  }

  verbose("Calculated chi squared at %ld \n",timeval_subtract(&startTime));


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
    if (tempqVec)
      XLALDestroyCOMPLEX8VectorSequence(tempqVec);
  }

  gsl_matrix_free(B2Null);
  gsl_eigen_symmv_free(matTempNull);
  gsl_matrix_free(eigenvecs);
  gsl_vector_free(eigenvals);
  gsl_matrix_free(eigenvecsSngl);
  gsl_vector_free(eigenvalsSngl);
  gsl_matrix_free(eigenvecsNull);
  gsl_vector_free(eigenvalsNull);

  *chisqOverlapsP = chisqOverlaps;
  *chisqSnglOverlapsP = chisqSnglOverlaps;

}

UINT8 coh_PTF_add_triggers(
    struct coh_PTF_params   *params,
    MultiInspiralTable      **eventList,
    MultiInspiralTable      **thisEvent,
    REAL4TimeSeries         *cohSNR,
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
    REAL4                   *timeOffsets
)
{
  // This function adds a trigger to the event list

  UINT4 i;
  INT4 j;
  UINT4 check;
  REAL4  deltaT    = cohSNR->deltaT;
  INT4 numPointCheck = floor(params->timeWindow/cohSNR->deltaT + 0.5);
  INT4   timeOffsetPoints[LAL_NUM_IFO];
  LIGOTimeGPS trigTime;
  MultiInspiralTable *lastEvent = *thisEvent;
  MultiInspiralTable *currEvent = NULL;
  UINT4 numDOF = 4;
  if ( params->singlePolFlag || params->faceOnStatistic )
    numDOF = 2;

  REAL4 cohSNRThreshold = params->threshold;
  if (spinTrigger)
    cohSNRThreshold = params->spinThreshold;

  for (i = 0; i < LAL_NUM_IFO; i++)
  {
    timeOffsetPoints[i] = (int) floor(timeOffsets[i]/deltaT + 0.5);
  }

  for (i = 0 ; i < cohSNR->data->length ; i++)
  {
    if (cohSNR->data->data[i] > cohSNRThreshold)
    {
      check = 1;
      for (j = ((INT4)i)-numPointCheck; j < ((INT4)i)+numPointCheck; j++)
      {
        if (j < 0)
          j = 0;
        if (j > (INT4)(cohSNR->data->length-1))
          break;
        if (cohSNR->data->data[j] > cohSNR->data->data[i])
        {
          check = 0;
          break;
        }
      }
      if (check) /* Add trigger to event list */
      {
        currEvent = (MultiInspiralTable *) 
                          LALCalloc(1, sizeof(MultiInspiralTable));

        currEvent->event_id = (EventIDColumn *) 
                                  LALCalloc(1, sizeof(EventIDColumn));
        currEvent->event_id->id=eventId;
        currEvent->time_slide_id = (EventIDColumn *)
                                       LALCalloc(1, sizeof(EventIDColumn));
        currEvent->time_slide_id->id=slideId;
        eventId++;
        trigTime = cohSNR->epoch;
        XLALGPSAdd(&trigTime,i*cohSNR->deltaT);
        currEvent->snr = cohSNR->data->data[i];
        currEvent->mass1 = PTFTemplate.mass1;
        currEvent->mass2 = PTFTemplate.mass2;
        currEvent->chi = PTFTemplate.chi;
        currEvent->kappa = PTFTemplate.kappa;
        currEvent->mchirp = PTFTemplate.totalMass*pow(PTFTemplate.eta,3.0/5.0);
        currEvent->eta = PTFTemplate.eta;
        currEvent->end_time = trigTime;
        /* add sky position, but need to track back to sky fixed sky position */
        currEvent->ra = rightAscension - 
                            XLALGreenwichMeanSiderealTime(&params->trigTime) + 
                            XLALGreenwichMeanSiderealTime(&currEvent->end_time);
        currEvent->dec = declination;
        if (params->doNullStream)
          currEvent->null_statistic = nullSNR->data->data[i];
        if (params->doTraceSNR)
          currEvent->trace_snr = traceSNR->data->data[i];
        if (params->doBankVeto)
        {
          if (params->numIFO != 1)
          {
            currEvent->bank_chisq = bankVeto[LAL_NUM_IFO]->data->data[i];
            currEvent->bank_chisq_dof = numDOF * params->BVsubBankSize;
          }
          if (params->doSnglChiSquared)
          {
            if (bankVeto[LAL_IFO_G1])
            {
              currEvent->bank_chisq_g = bankVeto[LAL_IFO_G1]->data->data[i];
            }
            if (bankVeto[LAL_IFO_H1])
            {
              currEvent->bank_chisq_h1 = bankVeto[LAL_IFO_H1]->data->data[i];
            }
            if (bankVeto[LAL_IFO_H2])
            {
              currEvent->bank_chisq_h2 = bankVeto[LAL_IFO_H2]->data->data[i];
            }
            if (bankVeto[LAL_IFO_L1])
            {
              currEvent->bank_chisq_l = bankVeto[LAL_IFO_L1]->data->data[i];
            }
            if (bankVeto[LAL_IFO_T1])
            {
              currEvent->bank_chisq_t = bankVeto[LAL_IFO_T1]->data->data[i];
            }
            if (bankVeto[LAL_IFO_V1])
            {
              currEvent->bank_chisq_v = bankVeto[LAL_IFO_V1]->data->data[i];
            }
          }
        }
        if (params->doAutoVeto)
        {
          if (params->numIFO != 1)
          {
            currEvent->cont_chisq = autoVeto[LAL_NUM_IFO]->data->data[i];
            currEvent->cont_chisq_dof = numDOF * params->numAutoPoints;
          }
          if (params->doSnglChiSquared)
          {
            if (autoVeto[LAL_IFO_G1])
            {
              currEvent->cont_chisq_g = autoVeto[LAL_IFO_G1]->data->data[i];
            }
            if (autoVeto[LAL_IFO_H1])
            {
              currEvent->cont_chisq_h1 = autoVeto[LAL_IFO_H1]->data->data[i];
            }
            if (autoVeto[LAL_IFO_H2])
            {
              currEvent->cont_chisq_h2 = autoVeto[LAL_IFO_H2]->data->data[i];
            }
            if (autoVeto[LAL_IFO_L1])
            {
              currEvent->cont_chisq_l = autoVeto[LAL_IFO_L1]->data->data[i];
            }
            if (autoVeto[LAL_IFO_T1])
            {
              currEvent->cont_chisq_t = autoVeto[LAL_IFO_T1]->data->data[i];
            }
            if (autoVeto[LAL_IFO_V1])
            {
              currEvent->cont_chisq_v = autoVeto[LAL_IFO_V1]->data->data[i];
            }
          }
        }
        if (params->doChiSquare)
        {
          if (params->numIFO != 1)
          {
            currEvent->chisq = chiSquare[LAL_NUM_IFO]->data->data[i];
            currEvent->chisq_dof = numDOF * (params->numChiSquareBins - 1);
          }
          if (params->doSnglChiSquared)
          {
            if (chiSquare[LAL_IFO_G1])
            {
              currEvent->chisq_g = chiSquare[LAL_IFO_G1]->data->data[i];
            }
            if (chiSquare[LAL_IFO_H1])
            {
              currEvent->chisq_h1 = chiSquare[LAL_IFO_H1]->data->data[i];
            }
            if (chiSquare[LAL_IFO_H2])
            {
              currEvent->chisq_h2 = chiSquare[LAL_IFO_H2]->data->data[i];
            }
            if (chiSquare[LAL_IFO_L1])
            {
              currEvent->chisq_l = chiSquare[LAL_IFO_L1]->data->data[i];
            }
            if (chiSquare[LAL_IFO_T1])
            {
              currEvent->chisq_t = chiSquare[LAL_IFO_T1]->data->data[i];
            }
            if (chiSquare[LAL_IFO_V1])
            {
              currEvent->chisq_v = chiSquare[LAL_IFO_V1]->data->data[i];
            }
          }
        }
        if (pValues[0])
          currEvent->amp_term_1 = pValues[0]->data->data[i];
        if (pValues[1]) 
          currEvent->amp_term_2 = pValues[1]->data->data[i];
        if (pValues[2]) 
          currEvent->amp_term_3 = pValues[2]->data->data[i];
        if (pValues[3]) 
          currEvent->amp_term_4 = pValues[3]->data->data[i];
        if (pValues[4]) 
          currEvent->amp_term_5 = pValues[4]->data->data[i];
        if (pValues[5]) 
          currEvent->amp_term_6 = pValues[5]->data->data[i];
        if (pValues[6]) 
          currEvent->amp_term_7 = pValues[6]->data->data[i];
        if (pValues[7]) 
          currEvent->amp_term_8 = pValues[7]->data->data[i];
        if (pValues[8]) 
          currEvent->amp_term_9 = pValues[8]->data->data[i];
        if (pValues[9]) 
          currEvent->amp_term_10 = pValues[9]->data->data[i];
        /* Note that these two terms are only used for debugging
        at the moment. When they are used properly they will be
        moved into sane columns! For spin they give Amp*cos(Phi_0) and
        Amp*sin(Phi_0). For non spinning the second is 0 and the
        first is some arbitrary amplitude. */
        if (gammaBeta[0])
        {
          currEvent->g1quad.re = gammaBeta[0]->data->data[i];
          currEvent->g1quad.im = gammaBeta[1]->data->data[i];
        }
        if (snrComps[LAL_IFO_G1])
        {
          currEvent->snr_g = snrComps[LAL_IFO_G1]->data->
                                 data[i+5000+timeOffsetPoints[LAL_IFO_G1]];
          currEvent->sigmasq_g = PTFM[LAL_IFO_G1]->data[0];
        }
        if (snrComps[LAL_IFO_H1])
        {
          currEvent->snr_h1 = snrComps[LAL_IFO_H1]->data->
                                  data[i+5000+timeOffsetPoints[LAL_IFO_H1]];
          currEvent->sigmasq_h1 = PTFM[LAL_IFO_H1]->data[0];
        }
        if (snrComps[LAL_IFO_H2])
        {
          currEvent->snr_h2 = snrComps[LAL_IFO_H2]->data->
                                  data[i+5000+timeOffsetPoints[LAL_IFO_H2]];
          currEvent->sigmasq_h2 = PTFM[LAL_IFO_H2]->data[0];
        }
        if (snrComps[LAL_IFO_L1])
        {
          currEvent->snr_l = snrComps[LAL_IFO_L1]->data->
                                 data[i+5000+timeOffsetPoints[LAL_IFO_L1]];
          currEvent->sigmasq_l = PTFM[LAL_IFO_L1]->data[0];
        }
        if (snrComps[LAL_IFO_T1])
        {
          currEvent->snr_t = snrComps[LAL_IFO_T1]->data->
                                 data[i+5000+timeOffsetPoints[LAL_IFO_T1]];
          currEvent->sigmasq_t = PTFM[LAL_IFO_T1]->data[0];
        }
        if (snrComps[LAL_IFO_V1])
        {
          currEvent->snr_v = snrComps[LAL_IFO_V1]->data->
                                 data[i+5000+timeOffsetPoints[LAL_IFO_V1]];
          currEvent->sigmasq_v = PTFM[LAL_IFO_V1]->data[0];
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
          if (params->numIFO == 1)
            currEvent->snr_dof = 2;
          else
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
          else if (coh_PTF_accept_trig_check(params,eventList,*currEvent) )
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
            if (currEvent->time_slide_id)
            {
              LALFree(currEvent->time_slide_id);
            }
            LALFree(currEvent);
          }
        }

      }
    }
  }
  *thisEvent = lastEvent;
  return eventId;
}

void coh_PTF_cluster_triggers(
    struct coh_PTF_params   *params,
    MultiInspiralTable      **eventList,
    MultiInspiralTable      **thisEvent
)
{

  /* This clustering function is currently unused. Currently clustering is
     done in post-processing, though this may need to be changed. */
  MultiInspiralTable *currEvent = *eventList;
  MultiInspiralTable *currEvent2 = NULL;
  MultiInspiralTable *newEvent = NULL;
  MultiInspiralTable *newEventHead = NULL;
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
   * clustering time */
  while (currEvent)
  {
    if (coh_PTF_accept_trig_check(params,eventList,*currEvent) )
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
      if (currEvent->time_slide_id)
      {
        LALFree(currEvent->time_slide_id);
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

UINT4 coh_PTF_accept_trig_check(
    struct coh_PTF_params   *params,
    MultiInspiralTable      **eventList,
    MultiInspiralTable      thisEvent
)
{
  /* This clustering function is currently unused. Currently clustering is
     done in post-processing, though this may need to be changed. */
  MultiInspiralTable *currEvent = *eventList;
  LIGOTimeGPS time1,time2;
  UINT4 loudTrigBefore=0,loudTrigAfter=0;

  currEvent = *eventList;

  /* for each trigger, find out whether a louder trigger is within the
   * clustering time */
  time1.gpsSeconds=thisEvent.end_time.gpsSeconds;
  time1.gpsNanoSeconds = thisEvent.end_time.gpsNanoSeconds;
  while (currEvent)
  {
    time2.gpsSeconds=currEvent->end_time.gpsSeconds;
    time2.gpsNanoSeconds=currEvent->end_time.gpsNanoSeconds;
    if (fabs(XLALGPSDiff(&time1,&time2)) < params->clusterWindow)
    {
      if (thisEvent.snr_dof == currEvent->snr_dof)
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
    }
    currEvent = currEvent->next;
  }

  return 1;
}

