#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include <lal/LALStdlib.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/Date.h>
#include <lal/RealFFT.h>
#include <lal/Ring.h>
#include <lal/FrameStream.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirpDatatypes.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpPTF.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/DetectorSite.h>
#include <lal/TimeDelay.h>
#include <lal/DetResponse.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#include "lalapps.h"
#include "getdata.h"
#include "injsgnl.h"
#include "getresp.h"
#include "spectrm.h"
#include "segment.h"
#include "errutil.h"

#include "coh_PTF.h"

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
    INT4                    segmentNumber
    );
void cohPTFaddTriggers(
    struct coh_PTF_params   *params,
    MultiInspiralTable      **eventList,
    REAL4TimeSeries         *cohSNR,
    InspiralTemplate        PTFTemplate,
    UINT4                   *eventId);


int main( int argc, char **argv )
{
  INT4 i,j,k;
  static LALStatus      status;
  struct coh_PTF_params      *params    = NULL;
  ProcessParamsTable      *procpar   = NULL;
  REAL4FFTPlan            *fwdplan   = NULL;
  REAL4FFTPlan            *revplan   = NULL;
  COMPLEX8FFTPlan          *invPlan   = NULL;
  REAL4TimeSeries         *channel[LAL_NUM_IFO];
  REAL4FrequencySeries    *invspec[LAL_NUM_IFO];
  RingDataSegments        *segments[LAL_NUM_IFO];
  INT4                    numTmplts = 0;
  INT4  startTemplate     = -1;           /* index of first template      */
  INT4  stopTemplate      = -1;           /* index of last template       */
  INT4 numSegments        = 0;
  InspiralTemplate        *PTFtemplate = NULL;
  FindChirpTemplate       *fcTmplt     = NULL;
  FindChirpTmpltParams     *fcTmpltParams      = NULL;
  FindChirpInitParams     *fcInitParams = NULL;
  UINT4                   numPoints,ifoNumber,spinTemplate;
  REAL8Array              *PTFM[LAL_NUM_IFO];
  COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO];
  time_t                  startTime;
  LALDetector             *detectors[LAL_NUM_IFO];
  REAL8                   *timeOffsets;
  REAL8                   *Fplus;
  REAL8                   *Fcross;
  REAL8                   rA = 1.1;
  REAL8                   dec = 1.1;
  REAL8                   detLoc[3];
  REAL4TimeSeries         *cohSNR = NULL;
  LIGOTimeGPS             segStartTime;
  MultiInspiralTable      *eventList = NULL;
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

  fprintf(stdout,"Read input params %ld \n", time(NULL)-startTime);

  /* create forward and reverse fft plans */
  fwdplan = coh_PTF_get_fft_fwdplan( params );
  revplan = coh_PTF_get_fft_revplan( params );

  fprintf(stdout,"Made fft plans %ld \n", time(NULL)-startTime);

  /* Determine if we are analyzing single or multiple ifo data */

  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if ( params->haveTrig[ifoNumber] )
    {
      numDetectors++;
    }
  }
  
  if (numDetectors == 0 )
  {
    fprintf(stderr,"You have not specified any detectors to analyse");
    return 1;
  }
  else if (numDetectors == 1 )
  {
    fprintf(stdout,"You have only specified one detector, why are you using the coherent code?");
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

      fprintf(stdout,"Created segments for one ifo %ld \n", time(NULL)-startTime);
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
      XLALGPSAdd(&segStartTime,(j+1)*params->segmentDuration/2.0);
      timeOffsets[j*LAL_NUM_IFO+ifoNumber] = 
          XLALTimeDelayFromEarthCenter(detLoc,rA,dec,&segStartTime);
      XLALComputeDetAMResponse(&Fplus[j*LAL_NUM_IFO+ifoNumber],
         &Fcross[j*LAL_NUM_IFO+ifoNumber],
         detectors[ifoNumber]->response,rA,dec,0.,
         XLALGreenwichMeanSiderealTime(&segStartTime));
    }
  }

  /* Create the relevant structures that will be needed */
  numPoints = floor( params->segmentDuration * params->sampleRate + 0.5 );
  fcInitParams = LALCalloc( 1, sizeof( *fcInitParams ));
  fcTmplt = LALCalloc( 1, sizeof( *fcTmplt ) );
  fcTmpltParams = LALCalloc ( 1, sizeof( *fcTmpltParams ) );
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
  /*fake_template (PTFtemplate);*/

  for (i = 0; (i < numTmplts); PTFtemplate = PTFtemplate->next, i++)
  {
    /* Determine if we can model this template as non-spinning */
    spinTemplate = 1;
    if (PTFtemplate->chi < 0.1) 
      spinTemplate = 0;
    PTFtemplate->approximant = FindChirpPTF;
    PTFtemplate->order = LAL_PNORDER_TWO;
    PTFtemplate->fLower = 40.;
    /* Generate the Q freq series of the template */
    generate_PTF_template(PTFtemplate,fcTmplt,fcTmpltParams);

    fprintf(stdout,"Generated template %d at %ld \n", i, time(NULL)-startTime);

    for ( j = 0; j < 1; ++j ) /* Loop over segments */
    {
      segStartTime = params->startTime;
      XLALGPSAdd(&segStartTime,(j+1)*params->segmentDuration/2.0);
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

          fprintf(stdout,
              "Made filters for ifo %d,segment %d, template %d at %ld \n", 
              ifoNumber,j,i,time(NULL)-startTime);
        }
      }
      
      /* Calculate the cohSNR time series */
      cohPTFmodBasesUnconstrainedStatistic(cohSNR,PTFM,PTFqVec,params,
                                 spinTemplate,singleDetector,timeOffsets,
                                 Fplus,Fcross,j);
     
      fprintf(stdout,
          "Made coherent statistic for segment %d, template %d at %ld \n",
          j,i,time(NULL)-startTime);      

      /* From this we want to construct triggers */
      cohPTFaddTriggers(params,&eventList,cohSNR,*PTFtemplate,&eventId);

      fprintf(stdout,
          "Generated triggers for segment %d, template %d at %ld \n",
          j,i,time(NULL)-startTime);
      XLALDestroyREAL4TimeSeries(cohSNR);
    }
  }
  cohPTF_output_events_xml( params->outputFile, eventList, procpar, params );
  exit(0);

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

    /* inject signals */
    if ( params->injectFile )
      inject_signal( channel, 3, params->injectFile,
          NULL, 1.0, NULL );
    if ( params->writeRawData )
       write_REAL4TimeSeries( channel );

    /* condition the data: resample and highpass */
    resample_REAL4TimeSeries( channel, params->sampleRate );
    if ( params->writeProcessedData ) /* write processed data */
      write_REAL4TimeSeries( channel );

    highpass_REAL4TimeSeries( channel, params->highpassFrequency );
    if ( params->writeProcessedData ) /* write processed data */
      write_REAL4TimeSeries( channel );

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
  template->fLower = 40.;
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
    
  outfile = fopen("cohSNR_timeseries.dat","w");
  for ( i = 0; i < cohSNR->data->length; ++i)
  {
    fprintf (outfile,"%f %f \n",deltaT*i,cohSNR->data->data[i]);
  }
  fclose(outfile);
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
  REAL4 *testing1;
  testing1  = LALCalloc(1, cohSNR->data->length*sizeof( REAL4 ));
  REAL4 *testing2;
  testing2  = LALCalloc(1, cohSNR->data->length*sizeof( REAL4 ));
  REAL4 *testing3;
  testing3  = LALCalloc(1, cohSNR->data->length*sizeof( REAL4 ));
  REAL4 *testing4;
  testing4  = LALCalloc(1, cohSNR->data->length*sizeof( REAL4 ));
  REAL4 *testing5;
  testing5  = LALCalloc(1, cohSNR->data->length*sizeof( REAL4 ));
  REAL4 *testing6;
  testing6  = LALCalloc(1, cohSNR->data->length*sizeof( REAL4 ));
  REAL4 *testing7;
  testing7  = LALCalloc(1, cohSNR->data->length*sizeof( REAL4 ));
  REAL4 *testing8;
  testing8  = LALCalloc(1, cohSNR->data->length*sizeof( REAL4 ));
  REAL4 *testing9;
  testing9  = LALCalloc(1, cohSNR->data->length*sizeof( REAL4 ));
  REAL4 *testing10;
  testing10  = LALCalloc(1, cohSNR->data->length*sizeof( REAL4 ));

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
      fprintf(stdout,"%f ",gsl_matrix_get(B2,i,j));
    }
    fprintf(stdout,"\n");
  }

  fprintf(stdout,"\n \n");

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
    testing1[i-numPoints/4] = u1[0];
    testing2[i-numPoints/4] = u1[1];
    testing6[i-numPoints/4] = u2[0];
    testing7[i-numPoints/4] = u2[1];
    if (spinTemplate == 1)
    {
      testing3[i-numPoints/4] = u1[2];
      testing4[i-numPoints/4] = u1[3];
      testing5[i-numPoints/4] = u1[4];
      testing8[i-numPoints/4] = u2[2];
      testing9[i-numPoints/4] = u2[3];
      testing10[i-numPoints/4] = u2[4];
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
  }

  outfile = fopen("cohSNR_timeseries.dat","w");
  for ( i = 0; i < cohSNR->data->length; ++i)
  {
    fprintf (outfile,"%f %f \n",deltaT*i,cohSNR->data->data[i]);
  }
  fclose(outfile);
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

}

void cohPTFaddTriggers(
    struct coh_PTF_params   *params,
    MultiInspiralTable      **eventList,
    REAL4TimeSeries         *cohSNR,
    InspiralTemplate        PTFTemplate,
    UINT4                   *eventId)
{
  UINT4 i,j;
  UINT4 check;
  INT4 numPointCheck = floor(params->timeWindow/cohSNR->deltaT + 0.5);
  LIGOTimeGPS trigTime;
  MultiInspiralTable *lastEvent = NULL;
  MultiInspiralTable *thisEvent = NULL;

  for (i = 0 ; i < cohSNR->data->length ; i++)
  {
    if (cohSNR->data->data[i] > params->threshold)
    {
      check = 1;
      for (j = i-numPointCheck; j < i+numPointCheck; j++)
      {
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
          thisEvent = *eventList;
        }
        else
        {
          lastEvent = thisEvent;
          thisEvent = (MultiInspiralTable *) 
              LALCalloc( 1, sizeof(MultiInspiralTable) );
          lastEvent->next = thisEvent;
        }
        thisEvent->event_id = (EventIDColumn *) 
            LALCalloc(1, sizeof(EventIDColumn) );
        thisEvent->event_id->id=eventId;
        *eventId++;
        trigTime = cohSNR->epoch;
        XLALGPSAdd(&trigTime,i*cohSNR->deltaT);
        thisEvent->snr = cohSNR->data->data[i];
        thisEvent->mass1 = PTFTemplate.mass1;
        thisEvent->mass2 = PTFTemplate.mass2;
/*        thisEvent->chi = PTFTemplate.chi;
        thisEvent->kappa = PTFTemplate.kappa;*/
        thisEvent->mchirp = PTFTemplate.totalMass*pow(PTFTemplate.eta,3.0/5.0);
        thisEvent->eta = PTFTemplate.eta;
        thisEvent->end_time = trigTime;
      }
    }
  }
}
        

  

