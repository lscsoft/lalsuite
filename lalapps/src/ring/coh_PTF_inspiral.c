#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

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
#include <lal/MatrixUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>

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
    REAL4Array                 *PTFM,
    COMPLEX8VectorSequence     *PTFqVec,
    COMPLEX8FrequencySeries    *sgmnt,
    COMPLEX8FFTPlan            *invPlan
    );
void cohPTFStatistic(
    REAL4Array                 *h1PTFM,
    COMPLEX8VectorSequence     *h1PTFqVec,
    REAL4Array                 *l1PTFM,
    COMPLEX8VectorSequence     *l1PTFqVec,
    REAL4Array                 *v1PTFM,
    COMPLEX8VectorSequence     *v1PTFqVec
    );
void cohPTFunconstrainedStatistic(
    REAL4Array              *PTFM[LAL_NUM_IFO],
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO],
    struct coh_PTF_params   *params,
    REAL4                   deltaT,
    UINT4                   numPoints
    );
void cohPTFoldStatistic(
    REAL4Array                 *h1PTFM,
    COMPLEX8VectorSequence     *h1PTFqVec,
    REAL4Array                 *l1PTFM,
    COMPLEX8VectorSequence     *l1PTFqVec,
    REAL4Array                 *v1PTFM,
    COMPLEX8VectorSequence     *v1PTFqVec
    );



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
  UINT4                   numPoints,ifoNumber;
  REAL4Array              *PTFM[LAL_NUM_IFO];
  COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO];
  time_t                  startTime;
  
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
      PTFM[ifoNumber] = XLALCreateArrayL( 2, 5, 5 );
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
    PTFtemplate->approximant = FindChirpPTF;
    PTFtemplate->order = LAL_PNORDER_TWO;
    PTFtemplate->fLower = 40.;
    /* Generate the Q freq series of the template */
    generate_PTF_template(PTFtemplate,fcTmplt,fcTmpltParams);

    fprintf(stdout,"Generated template %d at %ld \n", i, time(NULL)-startTime);

    for ( j = 0; j < numSegments; ++j )
    {
      for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if ( params->haveTrig[ifoNumber] )
        {
          /* Zero the storage vectors for the PTF filters */
          memset( PTFM[ifoNumber]->data, 0, 25 * sizeof(REAL4) );
          memset( PTFqVec[ifoNumber]->data, 0, 
                  5 * numPoints * sizeof(COMPLEX8) );

          /* And calculate A^I B^I and M^IJ */
          cohPTFNormalize(fcTmplt,invspec[ifoNumber],PTFM[ifoNumber],
              PTFqVec[ifoNumber],&segments[ifoNumber]->sgmnt[j],invPlan);
          fprintf(stdout,
              "Made filters for ifo %d,segment %d, template %d at %ld \n", 
              ifoNumber,j,i,time(NULL)-startTime);
        }
      }
    
      cohPTFunconstrainedStatistic(PTFM,PTFqVec,params,
                                 fcTmpltParams->deltaT,numPoints);
      fprintf(stdout,
          "Made coherent statistic for segment %d, template %d at %ld \n",
          j,i,time(NULL)-startTime);
    }
  }
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

void cohPTFStatistic(
    REAL4Array                 *h1PTFM,
    COMPLEX8VectorSequence     *h1PTFqVec,
    REAL4Array                 *l1PTFM,
    COMPLEX8VectorSequence     *l1PTFqVec,
    REAL4Array                 *v1PTFM,
    COMPLEX8VectorSequence     *v1PTFqVec)
{
  UINT4 numPoints,i,j,k,l,m;
  numPoints = h1PTFqVec->vectorLength;
  REAL4Array *PTFM[3];
  COMPLEX8VectorSequence *PTFqVec[3];
  PTFqVec[0] = h1PTFqVec;
  PTFqVec[1] = l1PTFqVec;
  PTFqVec[2] = v1PTFqVec;
  PTFM[0] = h1PTFM;
  PTFM[1] = l1PTFM;
  PTFM[2] = v1PTFM;
  REAL4 P[10];
  REAL4 a[3], b[3],theta[3],phi[3];
  REAL8 A[12];
  REAL8 CEmBF,AFmCD,BDmAE;
  REAL4 tempSNR[2],cohSNR[numPoints];
  REAL4 N,Ntmp[4];
  REAL4 beta,lamda,Theta,varphi;
  REAL4 ZH[25],SH[25],YU[25],zh[25],sh[25],yu[25];
  FILE *outfile;
  REAL4 deltaT = 1./4076.;
  REAL4 numpi = 3.141592654;
  UINT4 numvarphi = 10.;
  UINT4 numTheta = 10.;

  beta = 0.5;
  lamda = 2.13;

  /* For a given beta,lamda,Theta,varphi calculate the SNR */
  /* We need to calculate a[3] and b[3]. This function needs to be written
     properly (calculating theta and phi from beta and lamda for all ifos */
  theta[0] = 0.5;
  theta[1] = 1.1;
  theta[2] = 0.2;
  phi[0] = 1.5;
  phi[1] = 2.7;
  phi[0] = 0.03;
  for (i = 0; i < 3; i++)
  {
    a[i] = 0.5 * (1. + pow(cos(theta[i]),2.))*cos(2.*phi[i]);
    b[i] = cos(theta[i]) * sin(2.*phi[i]);
  }

  for (i = 0; i < numPoints; i++)
  {
    cohSNR[i] = 0.;
    /* Calculate the cyrillics */
    for (j = 0; j < 5; j++)
    {
      for (k = 0; k < 5; k++)
      {
        ZH[j*5+k] = 0;
        SH[j*5+k] = 0;
        YU[j*5+k] = 0;
        zh[j*5+k] = 0;
        sh[j*5+k] = 0;
        yu[j*5+k] = 0;
        for (l=0; l < 3; l++)
        {
          for (m=0; m < 3; m++)
          {
            Ntmp[0] = PTFqVec[l]->data[j*numPoints+i].im;  
            Ntmp[1] = PTFqVec[l]->data[k*numPoints+i].im;
            Ntmp[2] = PTFqVec[l]->data[j*numPoints+i].re; 
            Ntmp[3] = PTFqVec[l]->data[k*numPoints+i].re;
            N = Ntmp[0]*Ntmp[1] + Ntmp[2]*Ntmp[3];
            ZH[j*5+k] += a[l]*a[m] * N;
            SH[j*5+k] += b[l]*b[m] * N;
            YU[j*5+k] += 2.* a[l]*b[m] * N;
          }
          zh[j*5+k] += a[l]*a[l] * PTFM[l]->data[j*5+k];
          sh[j*5+k] += b[l]*b[l] * PTFM[l]->data[j*5+k];
          yu[j*5+k] += 2.*a[l]*b[l] * PTFM[l]->data[j*5+k];
        }
      }
    }
    /* Now we loop over varphi and Theta */
    for (l=0; l < numTheta; l++)
    {
      for (m=0; m < numvarphi; m++)
      {
        A[0] = 0;
        A[1] = 0;
        A[2] = 0;
        A[3] = 0;
        A[4] = 0;
        A[5] = 0;
        /* Define varphi and Theta */
        varphi = l/(REAL4) numvarphi * 2.*numpi;
        Theta = l/(REAL4) numTheta * numpi;
        /* Calculate the Ps */
        P[0] = -0.25 * cos(2.*varphi) * ( 3. + cos(2*Theta));
        P[1] = cos(Theta) * sin(2.*varphi);
        P[2] = -0.25 * sin(2.*varphi) * ( 3. + cos(2*Theta));
        P[3] = cos(Theta) * cos(2.*varphi);
        P[4] = 0.5 * cos(varphi) * sin(2.*Theta);
        P[5] = - sin(Theta) * sin(varphi);
        P[6] = 0.5 * sin(varphi) * sin(2.*Theta);
        P[7] = sin(Theta) * cos(varphi);
        P[8] = 0.5 * (1 - cos(2.*Theta));
        P[9] = 0;

        /* Calculate the A[6] */
        for (j = 0; j < 5; j++)
        {
          for (k = 0; k < 5; k++)
          {
            A[0]+=0.5*(P[2*j]*P[2*k] + P[2*j+1]*P[2*k+1])*(ZH[j*5+k]+SH[j*5+k]);
            A[1]+=0.5*(P[2*j]*P[2*k] - P[2*j+1]*P[2*k+1])*(ZH[j*5+k]+SH[j*5+k]);
            A[1]+=P[2*j]*P[2*k+1]*YU[j*5+k];
            A[2]+=P[2*j]*P[2*k+1]*(ZH[j*5+k]-SH[j*5+k]);
            A[2]+=0.5*(P[2*j+1]*P[2*k+1] - P[2*j]*P[2*k])*YU[j*5+k];
            A[3]+=0.5*(P[2*j]*P[2*k] + P[2*j+1]*P[2*k+1])*(zh[j*5+k]+sh[j*5+k]);
            A[4]+=0.5*(P[2*j]*P[2*k] - P[2*j+1]*P[2*k+1])*(zh[j*5+k]+sh[j*5+k]);
            A[4]+=P[2*j]*P[2*k+1]*yu[j*5+k];
            A[5]+=P[2*j]*P[2*k+1]*(zh[j*5+k]-sh[j*5+k]);
            A[5]+=0.5*(P[2*j+1]*P[2*k+1] - P[2*j]*P[2*k])*yu[j*5+k];       
          }
        }
        /* Now we calculate the SNR maximized over psi */
        BDmAE = A[1]*A[3]-A[0]*A[4];
        CEmBF = A[2]*A[4]-A[1]*A[5];
        AFmCD = A[0]*A[5] - A[2]*A[3];
        A[6] = pow(pow(AFmCD,2.) + pow(BDmAE,2.),0.5); 
        A[7] = pow(pow(A[6],2.) - pow(CEmBF,2.),0.5);
        A[8] = (A[7]*BDmAE + CEmBF*AFmCD) / pow(A[6],2.);
        A[9] = (CEmBF*BDmAE - A[7]*AFmCD) / pow(A[6],2.);
        A[10] = (-A[7]*BDmAE + CEmBF*AFmCD) / pow(A[6],2.);
        A[11] = (CEmBF*BDmAE + A[7]*AFmCD) / pow(A[6],2.);
        if ( A[8] > (1. + 1E-7) )
          fprintf(stderr,"Error cos psi is greater than one %lf",A[8]);
        if ( A[9] > (1. + 1E-7) )
          fprintf(stderr,"Error sin psi is greater than one %lf",A[9]);
        if ( A[10] > (1. + 1E-7) )
          fprintf(stderr,"Error cos psi 2 is greater than one %lf",A[10]);
        if ( A[11] > (1. + 1E-7) )
          fprintf(stderr,"Error sin psi 2 is greater than one %lf",A[11]);
        tempSNR[0] = (A[0]+A[1]*A[8]+A[2]*A[9])/(A[3]+A[4]*A[8]+A[5]*A[9]);
        tempSNR[1] = (A[0]+A[1]*A[10]+A[2]*A[11])/(A[3]+A[4]*A[10]+A[5]*A[11]);
        if ( tempSNR[0] > cohSNR[i])
          cohSNR[i] = tempSNR[0];
        if ( tempSNR[1] > cohSNR[i])
          cohSNR[i] = tempSNR[1];
      }
    }
  }
  outfile = fopen("cohSNR_timeseries.dat","w");
  for ( i = 0; i < numPoints; ++i)
  {
    fprintf (outfile,"%f %f \n",deltaT*i,cohSNR[i]);
  }
  fclose(outfile);

}

void cohPTFunconstrainedStatistic(
    REAL4Array              *PTFM[LAL_NUM_IFO],
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO],
    struct coh_PTF_params   *params, 
    REAL4                   deltaT,
    UINT4                   numPoints
)

{
  LALStatus status = blank_status;
  UINT4 i,j,k;
  REAL4 count;
  FILE *outfile;
  REAL8        *det         = NULL;
  REAL4 u1[10], u2[10], v1[10], v2[10];
/*  REAL8Array  *B, *Binv;*/
  REAL4 v1_dot_u1, v1_dot_u2, v2_dot_u1, v2_dot_u2,max_eigen;
  REAL4 a[LAL_NUM_IFO], b[LAL_NUM_IFO],theta[LAL_NUM_IFO],phi[LAL_NUM_IFO];
  REAL4 zh[25],sh[25],yu[25];
  REAL4 beta,lamda;
  REAL4 cohSNR[numPoints];

  gsl_matrix *B2 = gsl_matrix_alloc(10,10);
  gsl_matrix *Binv2a = gsl_matrix_alloc(10,10);
  gsl_matrix *Binv2 = gsl_matrix_alloc(10,10);
  gsl_permutation *p = gsl_permutation_alloc(10);

/*  B = XLALCreateREAL8ArrayL( 2, 10, 10 );
  Binv = XLALCreateREAL8ArrayL( 2, 10, 10 );
  memset( B->data, 0, 100 * sizeof(REAL8) );
  memset( Binv->data, 0, 100 * sizeof(REAL8) );*/

  beta = 0.5;
  lamda = 2.13;
  /* We need to calculate a[3] and b[3]. This function needs to be written
     properly (calculating theta and phi from beta and lamda for all ifos */
  theta[1] = 0.5;
  theta[3] = 1.1;
  theta[5] = 0.2;
  phi[1] = 1.5;
  phi[3] = 2.7;
  phi[5] = 0.03;
  for (i = 1; i < 6; i = i+2)
  {
    a[i] = 0.5 * (1. + pow(cos(theta[i]),2.))*cos(2.*phi[i]);
    b[i] = cos(theta[i]) * sin(2.*phi[i]);
  }
  
  /* Create and invert the Bmatrix */
  for (i = 0; i < 5; i++ )
  {
    for (j = 0; j < 5; j++ )
    {
      zh[i*5+j] = 0;
      sh[i*5+j] = 0;
      yu[i*5+j] = 0;
      for( k = 0; k < LAL_NUM_IFO; k++)
      {
        if ( params->haveTrig[k] )
        {
          zh[i*5+j] += a[k]*a[k] * PTFM[k]->data[i*5+j];
          sh[i*5+j] += b[k]*b[k] * PTFM[k]->data[i*5+j];
          yu[i*5+j] += a[k]*b[k] * PTFM[k]->data[i*5+j];
        }
      }
    }
  }

  for (i = 0; i < 10; i++ )
  {
    for (j = 0; j < 10; j++ )
    {
      if ( i < 5 && j < 5 )
      {
        gsl_matrix_set(B2,i,j,zh[i*5+j]);
      }
      else if ( i > 4 && j > 4)
      {
        gsl_matrix_set(B2,i,j,sh[(i-5)*5 + (j-5)]);
      }
      else if ( i < 5 && j > 4)
      {
        gsl_matrix_set(B2,i,j, yu[i*5 + (j-5)]);
      }
      else if ( i > 4 && j < 5)
      {
        gsl_matrix_set(B2,i,j,yu[j*5 + (i-5)]);
      }
      else
        fprintf(stderr,"BUGGER! Something went wrong.");
      gsl_matrix_set(Binv2,i,j,gsl_matrix_get(B2,i,j));
    }
  }

  /* This is the inversion of the B matrix */
  int signum;
  gsl_linalg_LU_decomp(Binv2,p,&signum);
  gsl_linalg_LU_invert(Binv2,p,Binv2a);


  for ( i = 0; i < numPoints; ++i ) /* Main loop over time */
  {
    for ( j = 0; j < 10 ; j++ ) /* Construct the vi vectors */
    {
      v1[j] = 0;
      v2[j] = 0;
      for( k = 0; k < LAL_NUM_IFO; k++)
      {
        if ( params->haveTrig[k] )
        {
          if (j < 5)
          {
            v1[j] += a[k] * PTFqVec[k]->data[j*numPoints+i].re;
            v2[j] += a[k] * PTFqVec[k]->data[j*numPoints+i].im;
          }
          else
          {
            v1[j] += b[k] * PTFqVec[k]->data[(j-5)*numPoints+i].re;
            v2[j] += b[k] * PTFqVec[k]->data[(j-5)*numPoints+i].im;
          }
        }
      }
    }
    for ( j = 0 ; j < 10 ; j ++ ) 
    {
      u1[j] = 0.0;
      u2[j] = 0.0;
      for ( k = 0; k < 10; k++ )
      {
        u1[j] += gsl_matrix_get(Binv2a,j,k) * v1[k];
        u2[j] += gsl_matrix_get(Binv2a,j,k) * v2[k];
      }
    }

    
    /* Compute the dot products */
    v1_dot_u1 = v1_dot_u2 = v2_dot_u1 = v2_dot_u2 = max_eigen = 0.0;
    for (j = 0; j < 10; j++)
    {
      v1_dot_u1 += v1[j] * u1[j];
      v1_dot_u2 += v1[j] * u2[j];
      v2_dot_u1 += v2[j] * u1[j];
      v2_dot_u2 += v2[j] * u2[j];
    }
    max_eigen = 0.5 * ( v1_dot_u1 + v2_dot_u2 + sqrt( (v1_dot_u1 - v2_dot_u2)
          * (v1_dot_u1 - v2_dot_u2) + 4 * v1_dot_u2 * v2_dot_u1 ));
/*    fprintf(stderr,"%f \n",max_eigen);*/
    cohSNR[i] = sqrt(max_eigen);
  }
    
  outfile = fopen("cohSNR_timeseries.dat","w");
  for ( i = 0; i < numPoints; ++i)
  {
    fprintf (outfile,"%f %f \n",deltaT*i,cohSNR[i]);
  }
  fclose(outfile);
}


void cohPTFoldStatistic(
    REAL4Array                 *h1PTFM,
    COMPLEX8VectorSequence     *h1PTFqVec,
    REAL4Array                 *l1PTFM,
    COMPLEX8VectorSequence     *l1PTFqVec,
    REAL4Array                 *v1PTFM,
    COMPLEX8VectorSequence     *v1PTFqVec)
{
  UINT4 numPoints,i,j,k;
  REAL4 MsumH,MsumL,MsumV,Msum;
  REAL4Vector *Asum,*Bsum,*SNR;
  REAL4 P[15];
  FILE *outfile;
  REAL4 deltaT = 1./4076.;
  numPoints = h1PTFqVec->vectorLength;

  Asum = XLALCreateREAL4Vector( numPoints );
  Bsum = XLALCreateREAL4Vector( numPoints );
  SNR = XLALCreateREAL4Vector( numPoints );
  P[0] = 1.;
  P[1] = 1.11503773;
  P[2] = 0.61949214;
  P[3] = 0.65593259;
  P[4] = 0.27879013;
  P[5] = 0.37759221;
  P[6] = 0.47075501;
  P[7] = 0.14902839;
  P[8] = 0.20966414;
  P[9] = 0.47096755;
  P[10] = 1.100864;
  P[11] = 1.2352822;
  P[12] = 0.66869988;
  P[13] = 0.7161475;
  P[14] = 0.36410693;
/*  P[0] = 0.;
  P[1] = 0.;
  P[2] = 0.;
  P[3] = 0.;
  P[4] = 0.;
  P[5] = 0.;
  P[6] = 0.;
  P[7] = 0.;
  P[8] = 0.;
  P[9] = 0.;
  P[10] = 0.;
  P[11] = 0.;
  P[12] = 0.;
  P[13] = 0.;
  P[14] = 0.; */
  

  memset( Asum->data, 0, Asum->length * sizeof(REAL4) );
  memset( Bsum->data, 0, Bsum->length * sizeof(REAL4) );
  memset( SNR->data, 0, Bsum->length * sizeof(REAL4) );

  for ( i = 0; i < numPoints ; i++ ) 
  {
    for ( j = 0; j < 5; j++ )
    {
      Asum->data[i] += P[j] * h1PTFqVec->data[i + j*numPoints].re;
      Asum->data[i] += P[j+5] * l1PTFqVec->data[i + j*numPoints].re;
      Asum->data[i] += P[j+10] * v1PTFqVec->data[i + j*numPoints].re;
      Bsum->data[i] += P[j] * h1PTFqVec->data[i + j*numPoints].im;
      Bsum->data[i] += P[j+5] * l1PTFqVec->data[i + j*numPoints].im;
      Bsum->data[i] += P[j+10] * v1PTFqVec->data[i + j*numPoints].im;
    }
  }
  
  MsumH = 0;
  MsumL = 0;
  MsumV = 0;
  Msum  = 0;
  for ( i = 0; i < 5 ; i++ )
  {
    for ( j = 0; j < 5; j++ )
    {
      MsumH += P[i]* P[j]* h1PTFM->data[i + j*5];
      MsumL += P[i+5]* P[j+5]* l1PTFM->data[i + j*5];
      MsumV += P[i+10]* P[j+10]* v1PTFM->data[i + j*5];
    }
  }

  Msum = MsumH + MsumL + MsumV;  
/*  Msum = pow(MsumH*MsumH+MsumL*MsumL+MsumV+MsumV,0.5);*/
  for ( i = 0; i < numPoints ; i++ )
  {
    SNR->data[i] = (pow(Asum->data[i],2.) + pow(Bsum->data[i],2))/Msum;
    SNR->data[i] = pow(SNR->data[i],0.5);
  }

  outfile = fopen("cohSNR_timeseries.dat","w");
  for ( i = 0; i < numPoints; ++i)
  {
    fprintf (outfile,"%f %f \n",deltaT*i,SNR->data[i]);
  }
  fclose(outfile);

}
  
    

  

