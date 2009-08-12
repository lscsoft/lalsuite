#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

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

static struct ring_params *ring_get_params( int argc, char **argv );
static REAL4FFTPlan *ring_get_fft_fwdplan( struct ring_params *params );
static REAL4FFTPlan *ring_get_fft_revplan( struct ring_params *params );
static REAL4TimeSeries *ring_get_data( struct ring_params *params,\
                       char *ifoChannel, char *dataCache );
static REAL4FrequencySeries *ring_get_invspec(
    REAL4TimeSeries         *channel,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *revplan,
    struct ring_params      *params
    );
void rescale_data (REAL4TimeSeries *channel,REAL8 rescaleFactor);
static RingDataSegments *ring_get_segments(
    REAL4TimeSeries         *channel,
    REAL4FrequencySeries    *invspec,
    REAL4FFTPlan            *fwdplan,
    struct ring_params      *params
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
    REAL4Array                 *PTFM,
    COMPLEX8VectorSequence     *PTFqVec,
    COMPLEX8FrequencySeries    *sgmnt,
    COMPLEX8FFTPlan            *invPlan
    );

int main( int argc, char **argv )
{
  static LALStatus      status;
  struct ring_params      *params    = NULL;
  ProcessParamsTable      *procpar   = NULL;
  REAL4FFTPlan            *fwdplan   = NULL;
  REAL4FFTPlan            *revplan   = NULL;
  COMPLEX8FFTPlan          *invPlan   = NULL;
  REAL4TimeSeries         *h1channel   = NULL;
  REAL4TimeSeries         *l1channel   = NULL;
  REAL4TimeSeries         *v1channel   = NULL;
  REAL4FrequencySeries    *h1invspec   = NULL;
  REAL4FrequencySeries    *l1invspec   = NULL;
  REAL4FrequencySeries    *v1invspec   = NULL;
  RingDataSegments        *h1segments  = NULL;
  RingDataSegments        *l1segments  = NULL;
  RingDataSegments        *v1segments  = NULL;
  InspiralTemplate        *PTFtemplate = NULL;
  FindChirpTemplate       *fcTmplt     = NULL;
  FindChirpTmpltParams     *fcTmpltParams      = NULL;
  FindChirpInitParams     *fcInitParams = NULL;
  UINT4                   numPoints;
  REAL4Array              *h1PTFM = NULL;
  COMPLEX8VectorSequence  *h1PTFqVec = NULL;
  REAL4Array              *l1PTFM = NULL;
  COMPLEX8VectorSequence  *l1PTFqVec = NULL;
  REAL4Array              *v1PTFM = NULL;
  COMPLEX8VectorSequence  *v1PTFqVec = NULL;

  /* set error handlers to abort on error */
  set_abrt_on_error();

  /* options are parsed and debug level is set here... */
  /* no lal mallocs before this! */
  params = ring_get_params( argc, argv );

  /* create process params */
/*  procpar = create_process_params( argc, argv, PROGRAM_NAME );*/

  /* create forward and reverse fft plans */
  fwdplan = ring_get_fft_fwdplan( params );
  revplan = ring_get_fft_revplan( params );

  /* Read in data from the various ifos */
  params->doubleData = 1;
  h1channel = ring_get_data( params,params->h1_channel,params->h1_dataCache );
  rescale_data (h1channel,1E20);
  l1channel = ring_get_data( params,params->l1_channel,params->l1_dataCache );
  rescale_data (l1channel,1E20);
  params->doubleData = 0;
  v1channel = ring_get_data( params,params->v1_channel,params->v1_dataCache );
  rescale_data (v1channel,1E20);

  /* compute the spectrum */
  h1invspec = ring_get_invspec( h1channel, fwdplan, revplan, params );
  l1invspec = ring_get_invspec( l1channel, fwdplan, revplan, params );
  v1invspec = ring_get_invspec( v1channel, fwdplan, revplan, params );

  /* create the segments */
  h1segments = ring_get_segments( h1channel,  h1invspec, fwdplan, params );
  l1segments = ring_get_segments( l1channel,  l1invspec, fwdplan, params );
  v1segments = ring_get_segments( v1channel,  v1invspec, fwdplan, params );

  /* Create the relevant structures that will be needed */
  numPoints = floor( params->segmentDuration * params->sampleRate + 0.5 );
  fcInitParams = LALCalloc( 1, sizeof( *fcInitParams ));
  PTFtemplate = LALCalloc( 1, sizeof( *PTFtemplate ) );
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

  /* Create an inverser FFT plan */
  invPlan = XLALCreateReverseCOMPLEX8FFTPlan( numPoints, 0 );

  /* Create storage vectors for the PTF filters */
  h1PTFM = XLALCreateArrayL( 2, 5, 5 );
  h1PTFqVec = XLALCreateCOMPLEX8VectorSequence ( 5, numPoints );
  l1PTFM = XLALCreateArrayL( 2, 5, 5 );
  l1PTFqVec = XLALCreateCOMPLEX8VectorSequence ( 5, numPoints );
  v1PTFM = XLALCreateArrayL( 2, 5, 5 );
  v1PTFqVec = XLALCreateCOMPLEX8VectorSequence ( 5, numPoints );
  memset( h1PTFM->data, 0, 25 * sizeof(REAL4) );
  memset( h1PTFqVec->data, 0, 5 * numPoints * sizeof(COMPLEX8) );
  memset( l1PTFM->data, 0, 25 * sizeof(REAL4) );
  memset( l1PTFqVec->data, 0, 5 * numPoints * sizeof(COMPLEX8) );
  memset( v1PTFM->data, 0, 25 * sizeof(REAL4) );
  memset( v1PTFqVec->data, 0, 5 * numPoints * sizeof(COMPLEX8) );

  /* A temporary call to create a template with specified values */
  fake_template (PTFtemplate);

  /* Generate the Q freq series of the template */
  generate_PTF_template(PTFtemplate,fcTmplt,fcTmpltParams);

  /* And calculate A^I B^I and M^IJ for every IFO */
  cohPTFNormalize(fcTmplt,h1invspec,h1PTFM,h1PTFqVec,
                  &h1segments->sgmnt[13],invPlan);
/*  cohPTFNormalize(fcTmplt,l1invspec,l1PTFM,l1PTFqVec,
                  &l1segments->sgmnt[0],invPlan);
  cohPTFNormalize(fcTmplt,h1invspec,v1PTFM,v1PTFqVec,
                  &v1segments->sgmnt[0],invPlan);*/
 
  exit(0);

}

/* warning: returns a pointer to a static variable... not reenterant */
/* only call this routine once to initialize params! */
/* also do not attempt to free this pointer! */
static struct ring_params *ring_get_params( int argc, char **argv )
{
  static struct ring_params params;
  static char programName[] = PROGRAM_NAME;
  static char cvsRevision[] = CVS_REVISION;
  static char cvsSource[]   = CVS_SOURCE;
  static char cvsDate[]     = CVS_DATE;
  ring_parse_options( &params, argc, argv );
  ring_params_sanity_check( &params ); /* this also sets various params */
  params.programName = programName;
  params.cvsRevision = cvsRevision;
  params.cvsSource   = cvsSource;
  params.cvsDate     = cvsDate;
  return &params;
}

/* gets the data, performs any injections, and conditions the data */
static REAL4TimeSeries *ring_get_data( struct ring_params *params,\
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
    if ( params->doubleData )
    {
      channel = get_frame_data_dbl_convert( dataCache, ifoChannel,
          &params->frameDataStartTime, params->frameDataDuration,
          params->strainData,
          params->geoHighpassFrequency, params->geoScale );
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

    /* inject ring signals */
/*    if ( params->injectFile )
      inject_signal( channel, params->injectType, params->injectFile,
          params->calibCache, 1.0, ifoChannel );
    if ( params->writeRawData )
       write_REAL4TimeSeries( channel );
    if( params->geoData ){
      for (j=0; j<channel->data->length; j++){
        channel->data->data[j] *= params->geoScale;
      }
    } */

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
static REAL4FFTPlan *ring_get_fft_fwdplan( struct ring_params *params )
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
static REAL4FFTPlan *ring_get_fft_revplan( struct ring_params *params )
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
static REAL4FrequencySeries *ring_get_invspec(
    REAL4TimeSeries         *channel,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *revplan,
    struct ring_params      *params
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
static RingDataSegments *ring_get_segments(
    REAL4TimeSeries         *channel,
    REAL4FrequencySeries    *invspec,
    REAL4FFTPlan            *fwdplan,
    struct ring_params      *params
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
