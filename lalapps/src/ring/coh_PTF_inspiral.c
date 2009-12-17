#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include "coh_PTF.h"

#include "lalapps.h"
#include "getdata.h"
#include "injsgnl.h"
#include "getresp.h"
#include "spectrm.h"
#include "segment.h"
#include "errutil.h"

RCSID( "$Id$" );
#define PROGRAM_NAME "lalapps_coh_PTF_inspiral"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE   "$Source$"
#define CVS_DATE     "$Date$"

static struct coh_PTF_params *coh_PTF_get_params( int argc, char **argv );
static REAL4FFTPlan *coh_PTF_get_fft_fwdplan( struct coh_PTF_params *params );
static REAL4FFTPlan *coh_PTF_get_fft_revplan( struct coh_PTF_params *params );
static REAL4TimeSeries *coh_PTF_get_data( struct coh_PTF_params *params,\
                       char *ifoChannel, char *dataCache );
static REAL4FrequencySeries *coh_PTF_get_invspec(
    REAL4TimeSeries         *channel,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *revplan,
    struct coh_PTF_params   *params
    );
void rescale_data (REAL4TimeSeries *channel,REAL8 rescaleFactor);
static RingDataSegments *coh_PTF_get_segments(
    REAL4TimeSeries         *channel,
    REAL4FrequencySeries    *invspec,
    REAL4FFTPlan            *fwdplan,
    struct coh_PTF_params      *params
    );
static int is_in_list( int i, const char *list );
void fake_template (InspiralTemplate *template);
void generate_PTF_template(
    InspiralTemplate         *PTFtemplate,
    FindChirpTemplate        *fcTmplt,
    FindChirpTmpltParams     *fcTmpltParams);
void cohPTFTemplate (
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *InspTmplt,
    FindChirpTmpltParams       *params    );
void
cohPTFNormalize(
    FindChirpTemplate          *fcTmplt,
    REAL4FrequencySeries       *invspec,
    REAL8Array                 *PTFM,
    COMPLEX8VectorSequence     *PTFqVec,
    COMPLEX8FrequencySeries    *sgmnt,
    COMPLEX8FFTPlan            *invPlan,
    UINT4                      spinTemplate
    );
void cohPTFunconstrainedStatistic(
    REAL4TimeSeries         *cohSNR,
    REAL8Array              *PTFM[LAL_NUM_IFO],
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO],
    struct coh_PTF_params   *params,
    UINT4                   spinTemplate,
    UINT4                   singleDetector,
    REAL8                   *timeOffsets,
    REAL8                   *Fplus,
    REAL8                   *Fcross,
    INT4                    segmentNumber
    );
void cohPTFmodBasesUnconstrainedStatistic(
    REAL4TimeSeries         *cohSNR,
    REAL8Array              *PTFM[LAL_NUM_IFO],
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO],
    struct coh_PTF_params   *params,
    UINT4                   spinTemplate,
    UINT4                   singleDetector,
    REAL8                   *timeOffsets,
    REAL8                   *Fplus,
    REAL8                   *Fcross,
    INT4                    segmentNumber,
    REAL4TimeSeries         *pValues[10],
    REAL4TimeSeries         *gammaBeta[2]
    );
void cohPTFaddTriggers(
    struct coh_PTF_params   *params,
    MultiInspiralTable      **eventList,
    MultiInspiralTable      **thisEvent,
    REAL4TimeSeries         *cohSNR,
    InspiralTemplate        PTFTemplate,
    UINT4                   *eventId,
    UINT4                   spinTrigger,
    REAL4TimeSeries         *pValues[10],
    REAL4TimeSeries         *gammaBeta[2]
    );
static void coh_PTF_cleanup(
    struct coh_PTF_params   *params,
    ProcessParamsTable      *procpar,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *revplan,
    COMPLEX8FFTPlan         *invPlan,
    REAL4TimeSeries         *channel[LAL_NUM_IFO],
    REAL4FrequencySeries    *invspec[LAL_NUM_IFO],
    RingDataSegments        *segments[LAL_NUM_IFO],
    MultiInspiralTable      *events,
    InspiralTemplate        *PTFbankhead,
    FindChirpTemplate       *fcTmplt,
    FindChirpTmpltParams    *fcTmpltParams,
    FindChirpInitParams     *fcInitParams,
    REAL8Array              *PTFM[LAL_NUM_IFO],
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO]
    );

int main( int argc, char **argv )
{
  INT4 i,j,k;
  static LALStatus      status;
  struct coh_PTF_params      *params    = NULL;
  ProcessParamsTable      *procpar   = NULL;
  REAL4FFTPlan            *fwdplan   = NULL;
  REAL4FFTPlan            *revplan   = NULL;
  COMPLEX8FFTPlan         *invPlan   = NULL;
  REAL4TimeSeries         *channel[LAL_NUM_IFO];
  REAL4FrequencySeries    *invspec[LAL_NUM_IFO];
  RingDataSegments        *segments[LAL_NUM_IFO];
  INT4                    numTmplts = 0;
  INT4  startTemplate     = -1;           /* index of first template      */
  INT4  stopTemplate      = -1;           /* index of last template       */
  INT4 numSegments        = 0;
  InspiralTemplate        *PTFtemplate = NULL;
  InspiralTemplate        *PTFbankhead = NULL;
  FindChirpTemplate       *fcTmplt     = NULL;
  FindChirpTmpltParams    *fcTmpltParams      = NULL;
  FindChirpInitParams     *fcInitParams = NULL;
  UINT4                   numPoints,ifoNumber,spinTemplate;
  REAL8Array              *PTFM[LAL_NUM_IFO];
  COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO];
  time_t                  startTime;
  LALDetector             *detectors[LAL_NUM_IFO];
  REAL8                   *timeOffsets;
  REAL8                   *Fplus;
  REAL8                   *Fcross;
  REAL8                   detLoc[3];
  REAL4TimeSeries         *cohSNR = NULL;
  REAL4TimeSeries         *pValues[10];
  REAL4TimeSeries         *gammaBeta[2];
  LIGOTimeGPS             segStartTime;
  MultiInspiralTable      *eventList = NULL;
  MultiInspiralTable      *thisEvent = NULL;
  UINT4                   eventId = 0;
  UINT4                   numDetectors = 0;
  UINT4                   singleDetector = 0;
  
  startTime = time(NULL);

  /* set error handlers to abort on error */
  set_abrt_on_error();

  /* options are parsed and debug level is set here... */
  

  /* no lal mallocs before this! */
  params = coh_PTF_get_params( argc, argv );

  /* create process params */
  procpar = create_process_params( argc, argv, PROGRAM_NAME );

  verbose("Read input params %ld \n", time(NULL)-startTime);

  /* create forward and reverse fft plans */
  fwdplan = coh_PTF_get_fft_fwdplan( params );
  revplan = coh_PTF_get_fft_revplan( params );

  verbose("Made fft plans %ld \n", time(NULL)-startTime);

  /* Determine if we are analyzing single or multiple ifo data */

  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if ( params->haveTrig[ifoNumber] )
    {
      numDetectors++;
    }
  }
  /* NULL out the pValues pointer array */
  for ( i = 0 ; i < 10 ; i++ )
  {
    pValues[i] = NULL;
  }   
  gammaBeta[0] = NULL;
  gammaBeta[1] = NULL;

  if (numDetectors == 0 )
  {
    fprintf(stderr,"You have not specified any detectors to analyse");
    return 1;
  }
  else if (numDetectors == 1 )
  {
    fprintf(stdout,"You have only specified one detector, why are you using the coherent code? \n");
    singleDetector = 1;
  }

  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if ( params->haveTrig[ifoNumber] )
    {
      /* Initialize some of the structures */
      channel[ifoNumber] = NULL;
      invspec[ifoNumber] = NULL;
      segments[ifoNumber] = NULL;
      PTFM[ifoNumber] = NULL;
      PTFqVec[ifoNumber] = NULL;
      /* Read in data from the various ifos */
      params->doubleData = 1;
      if ( params->simData )
          params->doubleData = 0;
      else if ( ifoNumber == LAL_IFO_V1 )
          params->doubleData = 0;
      channel[ifoNumber] = coh_PTF_get_data(params,params->channel[ifoNumber],\
                               params->dataCache[ifoNumber] );
      rescale_data (channel[ifoNumber],1E20);

      /* compute the spectrum */
      invspec[ifoNumber] = coh_PTF_get_invspec( channel[ifoNumber], fwdplan,\
                               revplan, params );

      /* create the segments */
      segments[ifoNumber] = coh_PTF_get_segments( channel[ifoNumber],\
           invspec[ifoNumber], fwdplan, params );
      
      numSegments = segments[ifoNumber]->numSgmnt;

      verbose("Created segments for one ifo %ld \n", time(NULL)-startTime);
    }
  }

  /* Determine time delays and response functions */ 

  timeOffsets = LALCalloc(1, numSegments*LAL_NUM_IFO*sizeof( REAL8 ));
  Fplus = LALCalloc(1, numSegments*LAL_NUM_IFO*sizeof( REAL8 ));
  Fcross = LALCalloc(1, numSegments*LAL_NUM_IFO*sizeof( REAL8 ));
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    detectors[ifoNumber] = LALCalloc( 1, sizeof( *detectors[ifoNumber] ));
    XLALReturnDetector(detectors[ifoNumber] ,ifoNumber);
    for ( i = 0; i < 3; i++ )
    {
      detLoc[i] = (double) detectors[ifoNumber]->location[i];
    }
    for ( j = 0; j < numSegments; ++j )
    {
      /* Despite being called segStartTime we use the time at the middle 
      * of a segment */
      segStartTime = params->startTime;
      
      /*XLALGPSAdd(&segStartTime,(j+1)*params->segmentDuration/2.0);*/
      XLALGPSAdd(&segStartTime,8.5*params->segmentDuration/2.0);
      /*XLALGPSMultiply(&segStartTime,0.);
      XLALGPSAdd(&segStartTime,874610713.072549154);*/
      timeOffsets[j*LAL_NUM_IFO+ifoNumber] = 
          XLALTimeDelayFromEarthCenter(detLoc,params->rightAscension,
          params->declination,&segStartTime);
      XLALComputeDetAMResponse(&Fplus[j*LAL_NUM_IFO+ifoNumber],
         &Fcross[j*LAL_NUM_IFO+ifoNumber],
         detectors[ifoNumber]->response,params->rightAscension,
         params->declination,0.,XLALGreenwichMeanSiderealTime(&segStartTime));
    }
  }

  /* Create the relevant structures that will be needed */
  numPoints = floor( params->segmentDuration * params->sampleRate + 0.5 );
  fcInitParams = LALCalloc( 1, sizeof( *fcInitParams ));
  fcTmplt = LALCalloc( 1, sizeof( *fcTmplt ) );
  fcTmpltParams = LALCalloc ( 1, sizeof( *fcTmpltParams ) );
  fcTmpltParams->approximant = FindChirpPTF;
  fcTmplt->PTFQtilde =
      XLALCreateCOMPLEX8VectorSequence( 5, numPoints / 2 + 1 );
/*  fcTmplt->PTFBinverse = XLALCreateArrayL( 2, 5, 5 );
  fcTmplt->PTFB = XLALCreateArrayL( 2, 5, 5 );*/
  fcTmpltParams->PTFQ = XLALCreateVectorSequence( 5, numPoints );
  fcTmpltParams->PTFphi = XLALCreateVector( numPoints );
  fcTmpltParams->PTFomega_2_3 = XLALCreateVector( numPoints );
  fcTmpltParams->PTFe1 = XLALCreateVectorSequence( 3, numPoints );
  fcTmpltParams->PTFe2 = XLALCreateVectorSequence( 3, numPoints );
  fcTmpltParams->fwdPlan =
        XLALCreateForwardREAL4FFTPlan( numPoints, 0 );
  fcTmpltParams->deltaT = 1.0/params->sampleRate;
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if ( params->haveTrig[ifoNumber] )
    {
      PTFM[ifoNumber] = XLALCreateREAL8ArrayL( 2, 5, 5 );
      PTFqVec[ifoNumber] = XLALCreateCOMPLEX8VectorSequence ( 5, numPoints );
    }
  }
  /* Create an inverser FFT plan */
  invPlan = XLALCreateReverseCOMPLEX8FFTPlan( numPoints, 0 );

  /* Read in the tmpltbank xml file */
  numTmplts = InspiralTmpltBankFromLIGOLw( &PTFtemplate, params->bankFile,
      startTemplate, stopTemplate );
  PTFbankhead = PTFtemplate;
  /*fake_template (PTFtemplate);*/

  for (i = 0; (i < numTmplts); PTFtemplate = PTFtemplate->next, i++)
  {
    /* Determine if we can model this template as non-spinning */
    spinTemplate = 1;
    if (PTFtemplate->chi < 0.1) 
      spinTemplate = 0;
    PTFtemplate->approximant = FindChirpPTF;
    PTFtemplate->order = LAL_PNORDER_TWO;
    PTFtemplate->fLower = 30.;
    /* Generate the Q freq series of the template */
    generate_PTF_template(PTFtemplate,fcTmplt,fcTmpltParams);

    verbose("Generated template %d at %ld \n", i, time(NULL)-startTime);

    for ( j = 0; j < numSegments; ++j ) /* Loop over segments */
    {
      for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if ( params->haveTrig[ifoNumber] )
        {
          segStartTime = segments[ifoNumber]->sgmnt[j].epoch;
          break;
        }
      }
      XLALGPSAdd(&segStartTime,params->segmentDuration/4.0);
      cohSNR = XLALCreateREAL4TimeSeries("cohSNR",
          &segStartTime,PTFtemplate->fLower,
          (1.0/params->sampleRate),&lalDimensionlessUnit,
          3*numPoints/4 - numPoints/4);
      for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if ( params->haveTrig[ifoNumber] )
        {
          /* Zero the storage vectors for the PTF filters */
          memset( PTFM[ifoNumber]->data, 0, 25 * sizeof(REAL8) );
          memset( PTFqVec[ifoNumber]->data, 0, 
                  5 * numPoints * sizeof(COMPLEX8) );

          /* And calculate A^I B^I and M^IJ */
          cohPTFNormalize(fcTmplt,invspec[ifoNumber],PTFM[ifoNumber],
              PTFqVec[ifoNumber],&segments[ifoNumber]->sgmnt[j],invPlan,
              spinTemplate);

          verbose("Made filters for ifo %d,segment %d, template %d at %ld \n", 
              ifoNumber,j,i,time(NULL)-startTime);
        }
      }
      
      /* Calculate the cohSNR time series */
      cohPTFmodBasesUnconstrainedStatistic(cohSNR,PTFM,PTFqVec,params,
                                 spinTemplate,singleDetector,timeOffsets,
                                 Fplus,Fcross,j,pValues,gammaBeta);
     
      verbose("Made coherent statistic for segment %d, template %d at %ld \n",
          j,i,time(NULL)-startTime);      

      /* From this we want to construct triggers */
      cohPTFaddTriggers(params,&eventList,&thisEvent,cohSNR,*PTFtemplate,&eventId,spinTemplate,pValues,gammaBeta);
      for ( k = 0 ; k < 10 ; k++ )
      {
        if (pValues[k])
            XLALDestroyREAL4TimeSeries(pValues[k]);
      }
      if (gammaBeta[0]) XLALDestroyREAL4TimeSeries(gammaBeta[0]);
      if (gammaBeta[1]) XLALDestroyREAL4TimeSeries(gammaBeta[1]);
      verbose("Generated triggers for segment %d, template %d at %ld \n",
          j,i,time(NULL)-startTime);
      XLALDestroyREAL4TimeSeries(cohSNR);
    }
  }
  cohPTF_output_events_xml( params->outputFile, eventList, procpar, params );

  verbose("Generated output xml file, cleaning up and exiting at %ld \n",
      time(NULL)-startTime);

  coh_PTF_cleanup(params,procpar,fwdplan,revplan,invPlan,channel,
      invspec,segments,eventList,PTFbankhead,fcTmplt,fcTmpltParams,
      fcInitParams,PTFM,PTFqVec);
  LALCheckMemoryLeaks();
  return 0;
}

/* warning: returns a pointer to a static variable... not reenterant */
/* only call this routine once to initialize params! */
/* also do not attempt to free this pointer! */
static struct coh_PTF_params *coh_PTF_get_params( int argc, char **argv )
{
  static struct coh_PTF_params params;
  static char programName[] = PROGRAM_NAME;
  static char cvsRevision[] = CVS_REVISION;
  static char cvsSource[]   = CVS_SOURCE;
  static char cvsDate[]     = CVS_DATE;
  coh_PTF_parse_options( &params, argc, argv );
  coh_PTF_params_sanity_check( &params ); /* this also sets various params */
  params.programName = programName;
  params.cvsRevision = cvsRevision;
  params.cvsSource   = cvsSource;
  params.cvsDate     = cvsDate;
  return &params;
}

/* gets the data, performs any injections, and conditions the data */
static REAL4TimeSeries *coh_PTF_get_data( struct coh_PTF_params *params,\
                       char *ifoChannel, char *dataCache  )
{
  int stripPad = 0;
  REAL4TimeSeries *channel = NULL;
  UINT4 j;

  /* compute the start and duration needed to pad data */
  params->frameDataStartTime = params->startTime;
  XLALGPSAdd( &params->frameDataStartTime, -1.0 * params->padData );
  params->frameDataDuration = params->duration + 2.0 * params->padData;

  if ( params->getData )
  {
    if ( params->simData )
      channel = get_simulated_data( ifoChannel, &params->startTime,
          params->duration, params->strainData, params->sampleRate,
          params->randomSeed, 1E-20 );
    else if ( params->zeroData )
    {
      channel = get_zero_data( ifoChannel, &params->startTime,
          params->duration, params->strainData, params->sampleRate );
    }
    else if ( params->doubleData )
    {
      channel = get_frame_data_dbl_convert( dataCache, ifoChannel,
          &params->frameDataStartTime, params->frameDataDuration,
          params->strainData,
          params->highpassFrequency, 1. );
      stripPad = 1;
    }
    else
    {
      channel = get_frame_data( dataCache, ifoChannel,
          &params->frameDataStartTime, params->frameDataDuration,
          params->strainData );
      stripPad = 1;
    }
    if ( params->writeRawData ) /* write raw data */
      write_REAL4TimeSeries( channel );

    /* Function to put injections overhead */
    snprintf( channel->name, LALNameLength * sizeof(CHAR), "ZENITH" );

    /* inject signals */
    if ( params->injectFile )
      inject_signal( channel, EOBNR_inject, params->injectFile,
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

/* gets the forward fft plan */
static REAL4FFTPlan *coh_PTF_get_fft_fwdplan( struct coh_PTF_params *params )
{
  REAL4FFTPlan *plan = NULL;
  if ( params->segmentDuration > 0.0 )
  {
    UINT4 segmentLength;
    segmentLength = floor( params->segmentDuration * params->sampleRate + 0.5 );
    plan = XLALCreateForwardREAL4FFTPlan( segmentLength, 0 );
  }
  return plan;
}


/* gets the reverse fft plan */
static REAL4FFTPlan *coh_PTF_get_fft_revplan( struct coh_PTF_params *params )
{
  REAL4FFTPlan *plan = NULL;
  if ( params->segmentDuration > 0.0 )
  {
    UINT4 segmentLength;
    segmentLength = floor( params->segmentDuration * params->sampleRate + 0.5 );
    plan = XLALCreateReverseREAL4FFTPlan( segmentLength, 0 );
  }
  return plan;
}

/* computes the inverse power spectrum */
static REAL4FrequencySeries *coh_PTF_get_invspec(
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
    invspec = compute_average_spectrum( channel, params->segmentDuration,
        params->strideDuration, fwdplan, params->whiteSpectrum );

    if ( params->writeInvSpectrum ) /* Write spectrum before inversion */
      write_REAL4FrequencySeries( invspec );

    /* invert spectrum */
    invert_spectrum( invspec, params->sampleRate, params->strideDuration,
        params->truncateDuration, params->lowCutoffFrequency, fwdplan,
        revplan );

    if ( params->writeInvSpectrum ) /* write inverse calibrated spectrum */
      write_REAL4FrequencySeries( invspec );
  }

  return invspec;
}

void rescale_data (REAL4TimeSeries *channel,REAL8 rescaleFactor)
{
  /* Function to dynamically rescale the data */
  UINT4 k;
  for ( k = 0; k < channel->data->length; ++k )
  {
    channel->data->data[k] *= rescaleFactor;
  }
}

/* creates the requested data segments (those in the list of segments to do) */
static RingDataSegments *coh_PTF_get_segments(
    REAL4TimeSeries         *channel,
    REAL4FrequencySeries    *invspec,
    REAL4FFTPlan            *fwdplan,
    struct coh_PTF_params   *params
    )
{
  RingDataSegments *segments = NULL;
  COMPLEX8FrequencySeries  *response = NULL;
  UINT4  sgmnt;

  segments = LALCalloc( 1, sizeof( *segments ) );

  /* TODO: trig start/end time condition */

  /* if todo list is empty then do them all */
  if ( ! params->segmentsToDoList || ! strlen( params->segmentsToDoList ) )
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
      if ( is_in_list( sgmnt, params->segmentsToDoList ) )
        ++count;
      else
        continue; /* skip this segment: it is not in todo list */

    if ( ! count ) /* no segments to do */
      return NULL;

    segments->numSgmnt = count;
    segments->sgmnt = LALCalloc( segments->numSgmnt, sizeof(*segments->sgmnt) );

    count = 0;
    for ( sgmnt = 0; sgmnt < params->numOverlapSegments; ++sgmnt )
      if ( is_in_list( sgmnt, params->segmentsToDoList ) )
        compute_data_segment( &segments->sgmnt[count++], sgmnt, channel,
            invspec, response, params->segmentDuration, params->strideDuration,
            fwdplan );
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
 
/* routine to see if integer i is in a list of integers to do */
/* e.g., 2, 7, and 222 are in the list "1-3,5,7-" but 4 is not */
#define BUFFER_SIZE 256
static int is_in_list( int i, const char *list )
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

void fake_template (InspiralTemplate *template)
{
  /* Define the various options that a template needs */
  template->approximant = FindChirpPTF;
  template->order = LAL_PNORDER_TWO;
  template->mass1 = 14.;
  template->mass2 = 1.;
  template->fLower = 30.;
  template->chi = 0.9;
  template->kappa = 0.1;
/*  template->t0 = 6.090556;
  template->t2 = 0.854636;
  template->t3 = 1.136940;
  template->t4 = 0.07991391;
  template->tC = 5.888166;
  template->fFinal = 2048;*/
}

void generate_PTF_template(
    InspiralTemplate         *PTFtemplate,
    FindChirpTemplate        *fcTmplt,
    FindChirpTmpltParams     *fcTmpltParams)
{
  cohPTFTemplate( fcTmplt,PTFtemplate, fcTmpltParams );
}

void cohPTFunconstrainedStatistic(
    REAL4TimeSeries         *cohSNR,
    REAL8Array              *PTFM[LAL_NUM_IFO],
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO],
    struct coh_PTF_params   *params, 
    UINT4                   spinTemplate,
    UINT4                   singleDetector,
    REAL8                   *timeOffsets,
    REAL8                   *Fplus,
    REAL8                   *Fcross,
    INT4                    segmentNumber
)

{
  UINT4 i,j,k,vecLength,vecLengthTwo,vecLengthSquare,vecLengthTwoSquare;
  INT4 timeOffsetPoints[LAL_NUM_IFO];
  REAL4 deltaT = cohSNR->deltaT;
  UINT4 numPoints = floor( params->segmentDuration * params->sampleRate + 0.5 );
  vecLength = 5;
  vecLengthTwo = 10;
  vecLengthSquare = 25;
  vecLengthTwoSquare = 100;
  if (spinTemplate == 0 && singleDetector == 1)
  {
    vecLength = 2;
    vecLengthTwo = 2;
    vecLengthSquare = 4;
    vecLengthTwoSquare = 4;
  }
  else if (spinTemplate == 0)
  {
    vecLength = 2;
    vecLengthTwo = 4;
    vecLengthSquare = 4;
    vecLengthTwoSquare = 16; 
  }
  else if (singleDetector == 1)
  {
    vecLengthTwo = 5;
    vecLengthTwoSquare = 25;
  }
  REAL4 count;
  FILE *outfile;
  REAL8        *det         = NULL;
/*  REAL8Array  *B, *Binv;*/
  REAL4 u1[vecLengthTwo],u2[vecLengthTwo],v1[vecLengthTwo],v2[vecLengthTwo];
  REAL4 v1_dot_u1, v1_dot_u2, v2_dot_u1, v2_dot_u2,max_eigen;
  REAL4 a[LAL_NUM_IFO], b[LAL_NUM_IFO],theta[LAL_NUM_IFO],phi[LAL_NUM_IFO];
  REAL4 zh[vecLengthSquare],sh[vecLengthSquare],yu[vecLengthSquare];
  REAL4 beta,lamda;

  gsl_matrix *B2 = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_matrix *Binv2 = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_permutation *p = gsl_permutation_alloc(vecLengthTwo);
  gsl_vector *u1vec = gsl_vector_alloc(vecLengthTwo);
  gsl_vector *u2vec = gsl_vector_alloc(vecLengthTwo);
  gsl_vector *v1vec = gsl_vector_alloc(vecLengthTwo);
  gsl_vector *v2vec = gsl_vector_alloc(vecLengthTwo);

/*  B = XLALCreateREAL8ArrayL( 2, 10, 10 );
  Binv = XLALCreateREAL8ArrayL( 2, 10, 10 );
  memset( B->data, 0, 100 * sizeof(REAL8) );
  memset( Binv->data, 0, 100 * sizeof(REAL8) );*/

  for (i = 0; i < LAL_NUM_IFO; i++)
  {
    a[i] = Fplus[segmentNumber*LAL_NUM_IFO+i];
    b[i] = Fcross[segmentNumber*LAL_NUM_IFO+i];
  }
  
  /* Create and invert the Bmatrix */
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
          /* Note that PTFM is always 5x5 even for non-spin */
          zh[i*vecLength+j] += a[k]*a[k] * PTFM[k]->data[i*5+j];
          sh[i*vecLength+j] += b[k]*b[k] * PTFM[k]->data[i*5+j];
          yu[i*vecLength+j] += a[k]*b[k] * PTFM[k]->data[i*5+j];
        }
      }
    }
  }

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
      gsl_matrix_set(Binv2,i,j,gsl_matrix_get(B2,i,j));
    }
  }

/*  for (i = 0; i < vecLengthTwo; i++ )
  {
    for (j = 0; j < vecLengthTwo; j++ )
    {
      fprintf(stderr,"%f ",gsl_matrix_get(B2,i,j));
    }
    fprintf(stderr,"\n");
  }*/ 


  /* This is the LU decomposition of the B matrix */
  int signum;
  gsl_linalg_LU_decomp(Binv2,p,&signum);
/*  gsl_linalg_LU_invert(Binv2,p,Binv2a);*/

  /* This loop takes the time offset in seconds and converts to time offset
  * in data points */
  for (i = 0; i < LAL_NUM_IFO; i++ )
  {
    timeOffsetPoints[i]=(int)(timeOffsets[segmentNumber*LAL_NUM_IFO+i]/deltaT);
  }

  for ( i = numPoints/4; i < 3*numPoints/4; ++i ) /* Main loop over time */
  {
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
            v1[j] += a[k] * PTFqVec[k]->data[j*numPoints+i+timeOffsetPoints[k]].re;
            v2[j] += a[k] * PTFqVec[k]->data[j*numPoints+i+timeOffsetPoints[k]].im;
          }
          else
          {
            v1[j] += b[k] * PTFqVec[k]->data[(j-vecLength)*numPoints+i+timeOffsetPoints[k]].re;
            v2[j] += b[k] * PTFqVec[k]->data[(j-vecLength)*numPoints+i+timeOffsetPoints[k]].im;
          }
        }
      }
      gsl_vector_set(v1vec,j,v1[j]);
      gsl_vector_set(v2vec,j,v2[j]);
    }
    gsl_linalg_LU_solve(Binv2,p,v1vec,u1vec);
    gsl_linalg_LU_solve(Binv2,p,v2vec,u2vec);
    for ( j = 0 ; j < vecLengthTwo ; j ++ ) 
    {
      u1[j] = gsl_vector_get(u1vec,j);
      u2[j] = gsl_vector_get(u2vec,j);
    }
    
    /* Compute the dot products */
    v1_dot_u1 = v1_dot_u2 = v2_dot_u1 = v2_dot_u2 = max_eigen = 0.0;
    for (j = 0; j < vecLengthTwo; j++)
    {
      v1_dot_u1 += v1[j] * u1[j];
      v1_dot_u2 += v1[j] * u2[j];
      v2_dot_u1 += v2[j] * u1[j];
      v2_dot_u2 += v2[j] * u2[j];
    }
    if (spinTemplate == 0)
    {
      max_eigen = 0.5 * ( v1_dot_u1 + v2_dot_u2 );
    }
    else
    {
      max_eigen = 0.5 * ( v1_dot_u1 + v2_dot_u2 + sqrt( (v1_dot_u1 - v2_dot_u2)
          * (v1_dot_u1 - v2_dot_u2) + 4 * v1_dot_u2 * v2_dot_u1 ));
    }
    /*fprintf(stdout,"%f %f %f %f\n",v1_dot_u1,v2_dot_u2,v1_dot_u2,v2_dot_u1);*/
    cohSNR->data->data[i-numPoints/4] = sqrt(max_eigen);
  }
    
  /*outfile = fopen("cohSNR_timeseries.dat","w");
  for ( i = 0; i < cohSNR->data->length; ++i)
  {
    fprintf (outfile,"%f %f \n",deltaT*i,cohSNR->data->data[i]);
  }
  fclose(outfile);*/
}

void cohPTFmodBasesUnconstrainedStatistic(
    REAL4TimeSeries         *cohSNR,
    REAL8Array              *PTFM[LAL_NUM_IFO],
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO],
    struct coh_PTF_params   *params,
    UINT4                   spinTemplate,
    UINT4                   singleDetector,
    REAL8                   *timeOffsets,
    REAL8                   *Fplus,
    REAL8                   *Fcross,
    INT4                    segmentNumber,
    REAL4TimeSeries         *pValues[10],
    REAL4TimeSeries         *gammaBeta[2]
)

{
  UINT4 i,j,k,vecLength,vecLengthTwo,vecLengthSquare,vecLengthTwoSquare;
  INT4 timeOffsetPoints[LAL_NUM_IFO];
  REAL4 deltaT = cohSNR->deltaT;
  UINT4 numPoints = floor( params->segmentDuration * params->sampleRate + 0.5 );
  vecLength = 5;
  vecLengthTwo = 10;
  vecLengthSquare = 25;
  vecLengthTwoSquare = 100;
  if (spinTemplate == 0 && singleDetector == 1)
  {
    vecLength = 2;
    vecLengthTwo = 2;
    vecLengthSquare = 4;
    vecLengthTwoSquare = 4;
  }
  else if (spinTemplate == 0)
  {
    vecLength = 2;
    vecLengthTwo = 4;
    vecLengthSquare = 4;
    vecLengthTwoSquare = 16;
  }
  else if (singleDetector == 1)
  {
    vecLengthTwo = 5;
    vecLengthTwoSquare = 25;
  }
  for ( i = 0 ; i < vecLengthTwo ; i++ )
  {
    pValues[i] = XLALCreateREAL4TimeSeries("Pvalue",
          &cohSNR->epoch,cohSNR->f0,cohSNR->deltaT,
          &lalDimensionlessUnit,cohSNR->data->length);
  }
  gammaBeta[0] = XLALCreateREAL4TimeSeries("Gamma",
          &cohSNR->epoch,cohSNR->f0,cohSNR->deltaT,
          &lalDimensionlessUnit,cohSNR->data->length);
  gammaBeta[1] = XLALCreateREAL4TimeSeries("Beta",
          &cohSNR->epoch,cohSNR->f0,cohSNR->deltaT,
          &lalDimensionlessUnit,cohSNR->data->length);

  REAL4 count;
  FILE *outfile;
  REAL8        *det         = NULL;
/*  REAL8Array  *B, *Binv;*/
  REAL4 u1[vecLengthTwo],u2[vecLengthTwo],v1[vecLengthTwo],v2[vecLengthTwo];
  REAL4 v1_dot_u1, v1_dot_u2, v2_dot_u1, v2_dot_u2,max_eigen;
  REAL4 recSNR;
  REAL4 dAlpha,dBeta,dCee;
  REAL4 pValsTemp[vecLengthTwo];
  REAL4 betaGammaTemp[2];
  REAL4 a[LAL_NUM_IFO], b[LAL_NUM_IFO],theta[LAL_NUM_IFO],phi[LAL_NUM_IFO];
  REAL4 zh[vecLengthSquare],sh[vecLengthSquare],yu[vecLengthSquare];
  REAL4 beta,lamda;

  gsl_matrix *B2 = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_matrix *Binv2 = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_permutation *p = gsl_permutation_alloc(vecLengthTwo);
  gsl_vector *u1vec = gsl_vector_alloc(vecLengthTwo);
  gsl_vector *u2vec = gsl_vector_alloc(vecLengthTwo);
  gsl_vector *v1vec = gsl_vector_alloc(vecLengthTwo);
  gsl_vector *v2vec = gsl_vector_alloc(vecLengthTwo);
  gsl_eigen_symmv_workspace *matTemp = gsl_eigen_symmv_alloc (vecLengthTwo);
  gsl_matrix *eigenvecs = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_vector *eigenvals = gsl_vector_alloc(vecLengthTwo);

  for (i = 0; i < LAL_NUM_IFO; i++)
  {
    a[i] = Fplus[segmentNumber*LAL_NUM_IFO+i];
    b[i] = Fcross[segmentNumber*LAL_NUM_IFO+i];
  }

  /* Create and invert the Bmatrix */
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
          /* Note that PTFM is always 5x5 even for non-spin */
          zh[i*vecLength+j] += a[k]*a[k] * PTFM[k]->data[i*5+j];
          sh[i*vecLength+j] += b[k]*b[k] * PTFM[k]->data[i*5+j];
          yu[i*vecLength+j] += a[k]*b[k] * PTFM[k]->data[i*5+j];
        }
      }
    }
  }

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
      gsl_matrix_set(Binv2,i,j,gsl_matrix_get(B2,i,j));
      /*fprintf(stdout,"%f ",gsl_matrix_get(B2,i,j));*/
    }
    /*fprintf(stdout,"\n");*/
  }

  /*fprintf(stdout,"\n \n");*/

  /* Here we compute the eigenvalues and eigenvectors of B2 */
  gsl_eigen_symmv (Binv2,eigenvals,eigenvecs,matTemp);

  for (i = 0; i < vecLengthTwo; i++ )
  {
    for (j = 0; j < vecLengthTwo; j++ )
    {
      fprintf(stdout,"%f ",gsl_matrix_get(eigenvecs,i,j));
    }
    fprintf(stdout,"\n");
  }

  fprintf(stdout,"\n \n");

  for (i = 0; i < vecLengthTwo; i++ )
  {
    fprintf(stdout,"%f ",gsl_vector_get(eigenvals,i));
  }

  fprintf(stdout,"\n \n");

  /* This loop takes the time offset in seconds and converts to time offset
  * in data points */
  for (i = 0; i < LAL_NUM_IFO; i++ )
  {
    timeOffsetPoints[i]=(int)(timeOffsets[segmentNumber*LAL_NUM_IFO+i]/deltaT);
  }

  for ( i = numPoints/4; i < 3*numPoints/4; ++i ) /* Main loop over time */
  {
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
            v1[j] += a[k] * PTFqVec[k]->data[j*numPoints+i+timeOffsetPoints[k]].re;
            v2[j] += a[k] * PTFqVec[k]->data[j*numPoints+i+timeOffsetPoints[k]].im;
          }
          else
          {
            v1[j] += b[k] * PTFqVec[k]->data[(j-vecLength)*numPoints+i+timeOffsetPoints[k]].re;
            v2[j] += b[k] * PTFqVec[k]->data[(j-vecLength)*numPoints+i+timeOffsetPoints[k]].im;
          }
        }
      }
    }

    /* Now we rotate the v1 and v2 to be in orthogonal basis */
    /* We can use gsl multiplication stuff to do this */
    /* BLAS stuff is stupid so we'll do it explicitly! */
   
    for ( j = 0 ; j < vecLengthTwo ; j++ )
    {
      u1[j] = 0.;
      u2[j] = 0.;
      for ( k = 0 ; k < vecLengthTwo ; k++ )
      {
        /* NOTE: not too sure whether its j,k or k,j */
        u1[j] += gsl_matrix_get(eigenvecs,k,j)*v1[k];
        u2[j] += gsl_matrix_get(eigenvecs,k,j)*v2[k];
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
      max_eigen = 0.5 * ( v1_dot_u1 + v2_dot_u2 );
    }
    else
    {
      max_eigen = 0.5 * ( v1_dot_u1 + v2_dot_u2 + sqrt( (v1_dot_u1 - v2_dot_u2)
          * (v1_dot_u1 - v2_dot_u2) + 4 * v1_dot_u2 * v1_dot_u2 ));
    }
    /*fprintf(stdout,"%f %f %f %f\n",v1_dot_u1,v2_dot_u2,v1_dot_u2,v2_dot_u1);*/
    cohSNR->data->data[i-numPoints/4] = sqrt(max_eigen);
    if (cohSNR->data->data[i-numPoints/4] > 27.00 )
    {
      /* IF louder than threshold calculate maximized quantities */
      v1_dot_u1 = v1_dot_u2 = v2_dot_u1 = v2_dot_u2 = 0;
      for (j = 0; j < vecLengthTwo; j++)
      {
        u1[j] = u1[j] / (pow(gsl_vector_get(eigenvals,j),0.5));
        u2[j] = u2[j] / (pow(gsl_vector_get(eigenvals,j),0.5));
        v1[j] = u1[j] * gsl_vector_get(eigenvals,j);
        v2[j] = u2[j] * gsl_vector_get(eigenvals,j);
        v1_dot_u1 += v1[j]*u1[j];
        v1_dot_u2 += v1[j]*u2[j];
        v2_dot_u2 += v2[j]*u2[j];
      }
      if ( spinTemplate == 1 )
      {
        dCee = (max_eigen - v1_dot_u1) / v1_dot_u2;
        dCee = 1.04707590;
        fprintf(stdout,"%f %f %f \n",max_eigen,v1_dot_u1,v1_dot_u2);
      }
      else
        dCee = 0;
      dAlpha = 1./(v1_dot_u1 + dCee * 2 * v1_dot_u2 + dCee*dCee*v2_dot_u2);
      dAlpha = pow(dAlpha,0.5);
      dBeta = dCee*dAlpha;
      fprintf(stdout,"dAlpha %f dBeta %f \n",dAlpha,dBeta);
      for ( j = 0 ; j < vecLengthTwo ; j++ )
      {
        pValsTemp[j] = dAlpha*u1[j] + dBeta*u2[j];  
        pValues[j]->data->data[i - numPoints/4] = 0;
        fprintf(stdout,"Rot v1: %f  v2: %f  u1:%e u2:%e P %e\n",v1[j],v2[j],u1[j],u2[j],pValsTemp[j]);
      } 
      recSNR = 0;
      for ( j = 0 ; j < vecLengthTwo ; j++ )
      {
        for ( k = 0 ; k < vecLengthTwo ; k++ )
        {
          recSNR += pValsTemp[j]*pValsTemp[k] * (v1[j]*v1[k]+v2[j]*v2[k]);
        }
      }
      /*fprintf(stdout,"%e %e \n",max_eigen,recSNR);*/
      betaGammaTemp[0] = 0;
      betaGammaTemp[1] = 0;
      for ( j = 0 ; j < vecLengthTwo ; j++ )
      {
        betaGammaTemp[0] += pValsTemp[j]*v1[j];
        betaGammaTemp[1] += pValsTemp[j]*v2[j];
      }
      gammaBeta[0]->data->data[i - numPoints/4] = betaGammaTemp[0];
      gammaBeta[1]->data->data[i - numPoints/4] = betaGammaTemp[1];

      for ( j = 0 ; j < vecLengthTwo ; j++ )
      {
        for ( k = 0 ; k < vecLengthTwo ; k++ )
        {
          pValues[j]->data->data[i-numPoints/4]+=gsl_matrix_get(eigenvecs,j,k)*pValsTemp[k];
        }
      }

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
              v1[j] += a[k] * PTFqVec[k]->data[j*numPoints+i+timeOffsetPoints[k]].re;
              v2[j] += a[k] * PTFqVec[k]->data[j*numPoints+i+timeOffsetPoints[k]].im;
            }
            else
            {
              v1[j] += b[k] * PTFqVec[k]->data[(j-vecLength)*numPoints+i+timeOffsetPoints[k]].re;
              v2[j] += b[k] * PTFqVec[k]->data[(j-vecLength)*numPoints+i+timeOffsetPoints[k]].im;
            }
          }
        }
      }
      recSNR = 0;
      for ( j = 0 ; j < vecLengthTwo ; j++ )
      {
        for ( k = 0 ; k < vecLengthTwo ; k++ )
        {
          recSNR += pValues[j]->data->data[i-numPoints/4]*pValues[k]->data->data[i-numPoints/4] * (v1[j]*v1[k]+v2[j]*v2[k]);
        }
        fprintf(stdout,"true  v1:%f v2:%f P %f\n",v1[j],v2[j],pValues[j]->data->data[i-numPoints/4]);
        
      }

      /*fprintf(stdout,"%e %e %e %e \n",v1[0],v1[1],v1[2],v1[3]);
      fprintf(stdout,"%e %e %e %e \n",v2[0],v2[1],v2[2],v2[3]);*/
      /*fprintf(stdout,"%e %e \n",betaGammaTemp[0],betaGammaTemp[1]);*/

    }
  }

  
  /*outfile = fopen("cohSNR_timeseries.dat","w");
  for ( i = 0; i < cohSNR->data->length; ++i)
  {
    fprintf (outfile,"%f %f \n",deltaT*i,cohSNR->data->data[i]);
  }
  fclose(outfile);
  */
  /*
  outfile = fopen("rebased_timeseries.dat","w");
  if (spinTemplate == 1 && singleDetector == 1 )
  {
    for ( i = 0; i < cohSNR->data->length; ++i)
    {
      fprintf (outfile,"%f %f %f %f %f %f %f %f %f %f %f\n",deltaT*i,testing1[i],testing2[i],testing3[i],testing4[i],testing5[i],testing6[i],testing7[i],testing8[i],testing9[i],testing10[i]);
    }
  }
  else if (singleDetector == 1 )
  {
    for ( i = 0; i < cohSNR->data->length; ++i)
    {
      fprintf (outfile,"%f %f %f %f %f \n",deltaT*i,testing1[i],testing2[i],testing6[i],testing7[i]);
    }
  }
  fclose(outfile);
  */
}

void cohPTFaddTriggers(
    struct coh_PTF_params   *params,
    MultiInspiralTable      **eventList,
    MultiInspiralTable      **thisEvent,
    REAL4TimeSeries         *cohSNR,
    InspiralTemplate        PTFTemplate,
    UINT4                   *eventId,
    UINT4                   spinTrigger,
    REAL4TimeSeries         *pValues[10],
    REAL4TimeSeries         *gammaBeta[2])
{
  INT4 i,j;
  UINT4 check;
  INT4 numPointCheck = floor(params->timeWindow/cohSNR->deltaT + 0.5);
  LIGOTimeGPS trigTime;
  MultiInspiralTable *lastEvent = NULL;
  MultiInspiralTable *currEvent = *thisEvent;

  for (i = 0 ; i < cohSNR->data->length ; i++)
  {
    if (cohSNR->data->data[i] > params->threshold)
    {
      check = 1;
      for (j = i-numPointCheck; j < i+numPointCheck; j++)
      {
        if (j < 0)
          j = 0;
        if (j > (cohSNR->data->length-1))
          break;
        if (cohSNR->data->data[j] > cohSNR->data->data[i])
        {
          check = 0;
          break;
        }
      }
      if (check) /* Add trigger to event list */
      {
        if ( !*eventList ) 
        {
          *eventList = (MultiInspiralTable *) 
              LALCalloc( 1, sizeof(MultiInspiralTable) );
          currEvent = *eventList;
        }
        else
        {
          lastEvent = currEvent;
          currEvent = (MultiInspiralTable *) 
              LALCalloc( 1, sizeof(MultiInspiralTable) );
          lastEvent->next = currEvent;
        }
        currEvent->event_id = (EventIDColumn *) 
            LALCalloc(1, sizeof(EventIDColumn) );
        currEvent->event_id->id=eventId;
        *eventId++;
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
        if (pValues[0])
          currEvent->h1quad.re = pValues[0]->data->data[i];
        if (pValues[1]) 
          currEvent->h1quad.im = pValues[1]->data->data[i];
        if (pValues[2]) 
          currEvent->h2quad.re = pValues[2]->data->data[i];
        if (pValues[3]) 
          currEvent->h2quad.im = pValues[3]->data->data[i];
        if (pValues[4]) 
          currEvent->l1quad.re = pValues[4]->data->data[i];
        if (pValues[5]) 
          currEvent->l1quad.im = pValues[5]->data->data[i];
        if (pValues[6]) 
          currEvent->v1quad.re = pValues[6]->data->data[i];
        if (pValues[7]) 
          currEvent->v1quad.im = pValues[7]->data->data[i];
        if (pValues[8]) 
          currEvent->t1quad.re = pValues[8]->data->data[i];
        if (pValues[9]) 
          currEvent->t1quad.im = pValues[9]->data->data[i];
        currEvent->g1quad.re = gammaBeta[0]->data->data[i];
        currEvent->g1quad.im = gammaBeta[1]->data->data[i];
        if (spinTrigger == 1)
          currEvent->snr_dof = 6;
        else
          currEvent->snr_dof = 2;
      }
    }
  }
  *thisEvent = currEvent;
}
        
static void coh_PTF_cleanup(
    struct coh_PTF_params   *params,
    ProcessParamsTable      *procpar,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *revplan,
    COMPLEX8FFTPlan         *invPlan,
    REAL4TimeSeries         *channel[LAL_NUM_IFO],
    REAL4FrequencySeries    *invspec[LAL_NUM_IFO],
    RingDataSegments        *segments[LAL_NUM_IFO],
    MultiInspiralTable      *events,
    InspiralTemplate        *PTFbankhead,
    FindChirpTemplate       *fcTmplt,
    FindChirpTmpltParams    *fcTmpltParams,
    FindChirpInitParams     *fcInitParams,
    REAL8Array              *PTFM[LAL_NUM_IFO],
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO]
    )
{
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
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if ( params->haveTrig[ifoNumber] )
    {
      if ( segments[LAL_NUM_IFO] )
      {
        for ( sgmnt = 0; sgmnt < segments[ifoNumber]->numSgmnt; ++sgmnt )
          if (segments[ifoNumber]->sgmnt[sgmnt].data)
            XLALDestroyCOMPLEX8Vector(segments[ifoNumber]->sgmnt[sgmnt].data);
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
      XLALDestroyREAL8Array( PTFM[ifoNumber] );
      XLALDestroyCOMPLEX8VectorSequence( PTFqVec[ifoNumber] );
    }
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
    LALFree( fcTmpltParams );
  }
  if ( fcTmplt )
  {
    if ( fcTmplt->PTFQtilde )
      XLALDestroyCOMPLEX8VectorSequence( fcTmplt->PTFQtilde );
    LALFree( fcTmplt );
  }
  if ( fcInitParams )
    LALFree( fcInitParams );

}

  

