#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>

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

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

static struct coh_PTF_params *coh_PTF_get_params( int argc, char **argv );
static REAL4FFTPlan *coh_PTF_get_fft_fwdplan( struct coh_PTF_params *params );
static REAL4FFTPlan *coh_PTF_get_fft_revplan( struct coh_PTF_params *params );
static REAL4TimeSeries *coh_PTF_get_data( struct coh_PTF_params *params,\
               const char *ifoChannel, const char *dataCache, UINT4 ifoNumber );
int coh_PTF_get_null_stream(
    struct coh_PTF_params *params,
    REAL4TimeSeries *channel[LAL_NUM_IFO + 1],
    REAL8 *Fplus,
    REAL8 *Fcross,
    REAL8 *timeOffsets );
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
    REAL8Array                 *PTFN,
    COMPLEX8VectorSequence     *PTFqVec,
    COMPLEX8FrequencySeries    *sgmnt,
    COMPLEX8FFTPlan            *invPlan,
    UINT4                      spinTemplate
    );
void cohPTFmodBasesUnconstrainedStatistic(
    REAL4TimeSeries         *cohSNR,
    REAL8Array              *PTFM[LAL_NUM_IFO+1],
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
    struct coh_PTF_params   *params,
    UINT4                   spinTemplate,
    UINT4                   singleDetector,
    REAL8                   *timeOffsets,
    REAL8                   *Fplus,
    REAL8                   *Fcross,
    INT4                    segmentNumber,
    REAL4TimeSeries         *pValues[10],
    REAL4TimeSeries         *gammaBeta[2],
    REAL4TimeSeries         *snrComps[LAL_NUM_IFO],
    REAL4TimeSeries         *nullSNR,
    REAL4TimeSeries         *traceSNR,
    REAL4TimeSeries         *bankVeto,
    REAL4TimeSeries         *autoVeto,
    REAL4TimeSeries         *chiSquare,
    UINT4                   subBankSize,
    struct bankComplexTemplateOverlaps *bankOverlaps,
    struct bankTemplateOverlaps *bankNormOverlaps,
    struct bankDataOverlaps *dataOverlaps,
    struct bankComplexTemplateOverlaps *autoTempOverlaps,
    FindChirpTemplate       *fcTmplt,
    REAL4FrequencySeries    *invspec[LAL_NUM_IFO+1],
    RingDataSegments        *segment[LAL_NUM_IFO+1],
    COMPLEX8FFTPlan         *invPlan
);
int cohPTFspinChecker(
    REAL8Array              *PTFM[LAL_NUM_IFO+1],
    REAL8Array              *PTFN[LAL_NUM_IFO+1],
    struct coh_PTF_params   *params,
    UINT4                   singleDetector,
    REAL8                   *Fplus,
    REAL8                   *Fcross,
    INT4                    segmentNumber
);
UINT8 cohPTFaddTriggers(
    struct coh_PTF_params   *params,
    MultiInspiralTable      **eventList,
    MultiInspiralTable      **thisEvent,
    REAL4TimeSeries         *cohSNR,
    InspiralTemplate        PTFTemplate,
    UINT8                   eventId,
    UINT4                   spinTrigger,
    UINT4                   singleDetector,
    REAL4TimeSeries         *pValues[10],
    REAL4TimeSeries         *gammaBeta[2],
    REAL4TimeSeries         *snrComps[LAL_NUM_IFO],
    REAL4TimeSeries         *nullSNR,
    REAL4TimeSeries         *traceSNR,
    REAL4TimeSeries         *bankVeto,
    REAL4TimeSeries         *autoVeto,
    REAL4TimeSeries         *chiSquare,
    REAL8Array              *PTFM[LAL_NUM_IFO+1]
);

void cohPTFclusterTriggers(
  MultiInspiralTable      **eventList,
  MultiInspiralTable      **thisEvent
);

static void coh_PTF_cleanup(
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
    REAL8                   *Fcross
    );

int main( int argc, char **argv )
{
// Declarations of parameters
  INT4 i,j,k;
  UINT4 ui,uj;
  struct coh_PTF_params      *params    = NULL;
  ProcessParamsTable      *procpar   = NULL;
  REAL4FFTPlan            *fwdplan   = NULL;
  REAL4FFTPlan            *revplan   = NULL;
  COMPLEX8FFTPlan         *invPlan   = NULL;
  REAL4TimeSeries         *channel[LAL_NUM_IFO+1];
  REAL4FrequencySeries    *invspec[LAL_NUM_IFO+1];
  RingDataSegments        *segments[LAL_NUM_IFO+1];
  INT4                    numTmplts = 0;
  INT4                    numSpinTmplts = 0;
  INT4                    numNoSpinTmplts = 0;
  INT4  startTemplate     = -1;           /* index of first template      */
  INT4  stopTemplate      = -1;           /* index of last template       */
  INT4 numSegments        = 0;
  InspiralTemplate        *PTFSpinTemplate = NULL;
  InspiralTemplate        *PTFNoSpinTemplate = NULL;
  InspiralTemplate        *PTFtemplate = NULL;
  InspiralTemplate        *PTFbankhead = NULL;
  FindChirpTemplate       *fcTmplt     = NULL;
  InspiralTemplate        *PTFBankTemplates = NULL;
  InspiralTemplate        *PTFBankvetoHead = NULL;
  FindChirpTemplate       *bankFcTmplts = NULL;
  FindChirpTmpltParams    *fcTmpltParams      = NULL;
  FindChirpInitParams     *fcInitParams = NULL;
  UINT4                   numPoints,ifoNumber,spinTemplate;
  REAL8Array              *PTFM[LAL_NUM_IFO+1];
  REAL8Array              *PTFN[LAL_NUM_IFO+1];
  COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1];
  time_t                  startTime;
  LALDetector             *detectors[LAL_NUM_IFO+1];
  REAL8                   *timeOffsets;
  REAL8                   *Fplus;
  REAL8                   *Fcross;
  REAL8                   detLoc[3];
  REAL4TimeSeries         *cohSNR = NULL;
  REAL4TimeSeries         *pValues[10];
  REAL4TimeSeries         *snrComps[LAL_NUM_IFO];
  REAL4TimeSeries         *gammaBeta[2];
  REAL4TimeSeries         *nullSNR = NULL;
  REAL4TimeSeries         *traceSNR = NULL;
  REAL4TimeSeries         *bankVeto = NULL;
  REAL4TimeSeries         *autoVeto = NULL;
  REAL4TimeSeries         *chiSquare = NULL;
  LIGOTimeGPS             segStartTime;
  MultiInspiralTable      *eventList = NULL;
  MultiInspiralTable      *thisEvent = NULL;
  UINT8                   eventId = 0;
  UINT4                   numDetectors = 0;
  UINT4                   singleDetector = 0;
  UINT4                   UNUSED spinBank = 0;
  char                    spinFileName[256];
  char                    noSpinFileName[256];
  
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
  /* NULL out pointers where necessary */
  for ( i = 0 ; i < 10 ; i++ )
  {
    pValues[i] = NULL;
  }   
  for ( i = 0 ; i < LAL_NUM_IFO ; i++ )
  {
    snrComps[i] = NULL;
  }
  gammaBeta[0] = NULL;
  gammaBeta[1] = NULL;

  /* Initialise some of the input file names */
  if ( params->spinBank )
  {
    spinBank = 1;
    strncpy(spinFileName,params->spinBank,sizeof(spinFileName)-1);
  }
  if ( params->noSpinBank )
    strncpy(noSpinFileName,params->noSpinBank,sizeof(noSpinFileName)-1);

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

  /* In this loop we read the data, generate segments and the PSD */
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    /* Initialize some of the structures */
    channel[ifoNumber] = NULL;
    invspec[ifoNumber] = NULL;
    segments[ifoNumber] = NULL;
    PTFM[ifoNumber] = NULL;
    PTFN[ifoNumber] = NULL;
    PTFqVec[ifoNumber] = NULL;
    if ( params->haveTrig[ifoNumber] )
    {
      /* Read in data from the various ifos */
      params->doubleData = 1;
      if ( params->simData )
          params->doubleData = 0;
      else if ( ifoNumber == LAL_IFO_V1 )
          params->doubleData = 0;
      channel[ifoNumber] = coh_PTF_get_data(params,params->channel[ifoNumber],\
                               params->dataCache[ifoNumber],ifoNumber );
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
  /* This is computed for all detectors, even if not being analyzed */ 
  timeOffsets = LALCalloc(1, LAL_NUM_IFO*sizeof( REAL8 ));
  Fplus = LALCalloc(1, LAL_NUM_IFO*sizeof( REAL8 ));
  Fcross = LALCalloc(1, LAL_NUM_IFO*sizeof( REAL8 ));
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
      segStartTime = params->trigTime;
      
      /*XLALGPSAdd(&segStartTime,(j+1)*params->segmentDuration/2.0);*/
//      XLALGPSAdd(&segStartTime,8.5*params->segmentDuration/2.0);
      /*XLALGPSMultiply(&segStartTime,0.);
      XLALGPSAdd(&segStartTime,874610713.072549154);*/
      timeOffsets[ifoNumber] = 
          XLALTimeDelayFromEarthCenter(detLoc,params->rightAscension,
          params->declination,&segStartTime);
      XLALComputeDetAMResponse(&Fplus[ifoNumber],
         &Fcross[ifoNumber],
         detectors[ifoNumber]->response,params->rightAscension,
         params->declination,0.,XLALGreenwichMeanSiderealTime(&segStartTime));
    }
    LALFree(detectors[ifoNumber]);
  }
  

  numPoints = floor( params->segmentDuration * params->sampleRate + 0.5 );

  /* Initialize some of the structures */
  ifoNumber = LAL_NUM_IFO;
  channel[ifoNumber] = NULL;
  invspec[ifoNumber] = NULL;
  segments[ifoNumber] = NULL;
  PTFM[ifoNumber] = NULL;
  PTFN[ifoNumber] = NULL;
  PTFqVec[ifoNumber] = NULL;


  /* Construct the null stream, its segments and its PSD */
  if ( params->doNullStream )
  {
    /* Read in data from the various ifos */
    if (coh_PTF_get_null_stream(params,channel,Fplus,Fcross,timeOffsets ))
    {
      fprintf(stderr,"Null stream construction failure\n");
      return 1;
    }

    /* compute the spectrum */
    invspec[ifoNumber] = coh_PTF_get_invspec( channel[ifoNumber], fwdplan,\
                             revplan, params );
    /* If white spectrum need to scale this. FIX ME!!! */
    if (params->whiteSpectrum)
    {
      for( ui=0 ; ui < invspec[ifoNumber]->data->length; ui++)
      {
        invspec[ifoNumber]->data->data[ui] *= pow(1./0.3403324,2);
      }
    }

    /* create the segments */
    segments[ifoNumber] = coh_PTF_get_segments( channel[ifoNumber],\
         invspec[ifoNumber], fwdplan, params );

    numSegments = segments[ifoNumber]->numSgmnt;

    verbose("Created segments for null stream at %ld \n", time(NULL)-startTime);
    PTFM[ifoNumber] = XLALCreateREAL8ArrayL( 2, 5, 5 );
    PTFN[ifoNumber] = XLALCreateREAL8ArrayL( 2, 5, 5 );
    PTFqVec[ifoNumber] = XLALCreateCOMPLEX8VectorSequence ( 5, numPoints );
  }

  /* At this point we can discard the calibrated data, only the segments
     and spectrum is needed now */

  for( ifoNumber = 0; ifoNumber < (LAL_NUM_IFO+1); ifoNumber++)
  {
    if ( channel[ifoNumber] )
    {
      XLALDestroyREAL4Vector( channel[ifoNumber]->data );
      LALFree( channel[ifoNumber] );
      channel[ifoNumber] = NULL;
    }
  }

  /* Create the relevant structures that will be needed */
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
      PTFN[ifoNumber] = XLALCreateREAL8ArrayL( 2, 5, 5 );
      PTFqVec[ifoNumber] = XLALCreateCOMPLEX8VectorSequence ( 5, numPoints );
    }
  }
  /* Create an inverse FFT plan */
  invPlan = XLALCreateReverseCOMPLEX8FFTPlan( numPoints, 0 );

  /* Read in the tmpltbank xml files */
  if ( params->spinBank )
  {
    numSpinTmplts = InspiralTmpltBankFromLIGOLw( &PTFSpinTemplate,
      spinFileName,startTemplate, stopTemplate );
    if (numSpinTmplts != 0 )
    {
      PTFtemplate = PTFSpinTemplate;
      numTmplts = numSpinTmplts;
    }
    else
      params->spinBank = NULL;
      spinBank = 0;
  }
  if ( params->noSpinBank )
  {
    numNoSpinTmplts = InspiralTmpltBankFromLIGOLw( &PTFNoSpinTemplate,
      noSpinFileName,startTemplate, stopTemplate );
    if ( numNoSpinTmplts != 0 )
    {
      PTFtemplate = PTFNoSpinTemplate;
      numTmplts = numNoSpinTmplts;
    }
    else
      params->noSpinBank = NULL;
  }
  /* If both banks present combine them and mark where to swap over */
  if ( params->spinBank && params->noSpinBank )
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

  /* Create the templates needed for the bank veto, if necessary */
  UINT4 subBankSize = 0;
  struct bankTemplateOverlaps *bankNormOverlaps = NULL;
  struct bankComplexTemplateOverlaps *bankOverlaps = NULL;
  struct bankDataOverlaps *dataOverlaps = NULL;
  
  if ( params->doBankVeto )
  {
    /* Reads in and initializes the bank veto sub bank */
    subBankSize = read_sub_bank(params,&PTFBankTemplates);
    bankNormOverlaps = LALCalloc( subBankSize,sizeof( *bankNormOverlaps));
    bankOverlaps = LALCalloc( subBankSize,sizeof( *bankOverlaps));
    dataOverlaps = LALCalloc(subBankSize,sizeof( *dataOverlaps));
    bankFcTmplts = LALCalloc( subBankSize, sizeof( *bankFcTmplts ));
    /* Create necessary structure to hold Q(f) */
    for (ui =0 ; ui < subBankSize; ui++)
    {
      bankFcTmplts[ui].PTFQtilde = 
          XLALCreateCOMPLEX8VectorSequence( 1, numPoints / 2 + 1 );
    }
    PTFBankvetoHead = PTFBankTemplates;
    
    for ( ui=0 ; ui < subBankSize ; ui++ )
    {
/*      if (! spinBank)
      {*/
      /* Generate the bank veto templates */
      generate_PTF_template(PTFBankTemplates,fcTmplt,
          fcTmpltParams);
      PTFBankTemplates = PTFBankTemplates->next;
      /* Only store Q1 */
      for ( uj = 0 ; uj < (numPoints/2 +1) ; uj++ )
      {
        bankFcTmplts[ui].PTFQtilde->data[uj] = fcTmplt->PTFQtilde->data[uj];
      }
/*      } 
      else
        generate_PTF_template(&(PTFBankTemplates[ui]),&(bankFcTmplts[ui]),
            fcTmpltParams);*/
    }
    /* Calculate the overlap between templates for bank veto */
    for ( ui = 0 ; ui < subBankSize; ui++ )
    {
      for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if ( params->haveTrig[ifoNumber] )
        {
/*          if ( spinBank )
          {
            bankNormOverlaps[ui].PTFM[ifoNumber]=
                XLALCreateREAL8ArrayL( 2, 5, 5 );
            memset( bankNormOverlaps[ui].PTFM[ifoNumber]->data,
                0, 25 * sizeof(REAL8) );
          }
          else
          {*/
          /* Initialize and zero out the structures */
          bankNormOverlaps[ui].PTFM[ifoNumber]=
              XLALCreateREAL8ArrayL( 2, 1, 1 );
          memset( bankNormOverlaps[ui].PTFM[ifoNumber]->data,
              0, 1 * sizeof(REAL8) );
//          }
          /* This function calculates the overlaps between templates */
          /* This returns a REAL4 as the overlap between identical templates*/
          /* must be real. */
          cohPTFTemplateOverlaps(&(bankFcTmplts[ui]),&(bankFcTmplts[ui]),
              invspec[ifoNumber],0,
              bankNormOverlaps[ui].PTFM[ifoNumber]);
        }
      }
    }
 
    verbose("Generated bank veto filters at %ld \n", time(NULL)-startTime);
        
  }

  /* Create the structures needed for the auto veto, if necessary */
  UINT4 timeStepPoints = 0;
  struct bankComplexTemplateOverlaps *autoTempOverlaps = NULL;

  if ( params->doAutoVeto )
  {
    /* Initializations */
    autoTempOverlaps = LALCalloc( params->numAutoPoints,
        sizeof( *autoTempOverlaps));
    timeStepPoints = params->autoVetoTimeStep*params->sampleRate;
    for (uj = 0; uj < params->numAutoPoints; uj++ )
    {
      for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if ( params->haveTrig[ifoNumber] )
        {
          /* If it will be used initialize and zero out the overlap structure*/
          autoTempOverlaps[uj].PTFM[ifoNumber] =
              XLALCreateCOMPLEX8ArrayL( 2, 1, 1 );
          memset( autoTempOverlaps[uj].PTFM[ifoNumber]->data,
              0, 1 * sizeof(COMPLEX8) );
        }
        else
          autoTempOverlaps[uj].PTFM[ifoNumber] = NULL;
      }
    }
  }
    

  PTFbankhead = PTFtemplate;
  /*fake_template (PTFtemplate);*/
  for ( j = 0; j < numSegments; ++j ) /* Loop over segments */
  {
    if ( params->doBankVeto )
    {
      /* Calculate overlap between bank templates and data for bank veto */
      for ( ui = 0 ; ui < subBankSize ; ui++ ) /* Loop over bank veto temps*/
      {
        for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
        {
          if ( params->haveTrig[ifoNumber] )
          {
/*            if (spinBank)
            {
              dataOverlaps[ui].PTFqVec[ifoNumber] =
                  XLALCreateCOMPLEX8VectorSequence ( 5, 3*numPoints/4 - numPoints/4 + 10000);
              bankOverlaps[ui].PTFM[ifoNumber]=XLALCreateCOMPLEX8ArrayL(2,5,5);
            }
            else
            {*/
            dataOverlaps[ui].PTFqVec[ifoNumber] =
                XLALCreateCOMPLEX8VectorSequence ( 1, 
                3*numPoints/4 - numPoints/4 + 10000);
            bankOverlaps[ui].PTFM[ifoNumber]=XLALCreateCOMPLEX8ArrayL(2,1,1);
//            }
            /* This function calculates the overlap */
            cohPTFBankFilters(&(bankFcTmplts[ui]),0,
                &segments[ifoNumber]->sgmnt[j],invPlan,PTFqVec[ifoNumber],
                dataOverlaps[ui].PTFqVec[ifoNumber],0,0);
          }
        }
      }
      verbose("Generated bank veto filters for segment %d at %ld \n",j, time(NULL)-startTime);
    }
    PTFtemplate = PTFbankhead;

    /* Loop over templates in the bank */
    for (i = 0; (i < numTmplts); PTFtemplate = PTFtemplate->next, i++)
    {
      /* Determine if this template is non-spinning */
      if (i >= numNoSpinTmplts)
        spinTemplate = 1;
      else
        spinTemplate = 0;
      PTFtemplate->approximant = FindChirpPTF;
      PTFtemplate->order = LAL_PNORDER_TWO;
      PTFtemplate->fLower = 38.;
      /* Generate the template */
      generate_PTF_template(PTFtemplate,fcTmplt,fcTmpltParams);

      if (spinTemplate)
        verbose("Generated spin template %d at %ld \n",i,time(NULL)-startTime);
      else
        verbose("Generated no spin template %d at %ld \n",i,time(NULL)-startTime);
      for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if ( params->haveTrig[ifoNumber] )
        {
          segStartTime = segments[ifoNumber]->sgmnt[j].epoch;
          break;
        }
      }
      /* Here we figure out the start time for the segment */
      XLALGPSAdd(&segStartTime,params->segmentDuration/4.0);

      /* Generate the various time series as needed*/
      /* Need to zero these out */
      cohSNR = XLALCreateREAL4TimeSeries("cohSNR",
          &segStartTime,PTFtemplate->fLower,
          (1.0/params->sampleRate),&lalDimensionlessUnit,
          3*numPoints/4 - numPoints/4);
      if (params->doTraceSNR)
      {
        for (ifoNumber = 0;ifoNumber < LAL_NUM_IFO; ifoNumber++)
        {
          if ( params->haveTrig[ifoNumber] )
          {
            snrComps[ifoNumber] = XLALCreateREAL4TimeSeries("snrComps",
                &cohSNR->epoch,cohSNR->f0,cohSNR->deltaT,
                &lalDimensionlessUnit,cohSNR->data->length);
          }
        }
      }
      if (params->doNullStream)
        nullSNR = XLALCreateREAL4TimeSeries("nullSNR",
            &segStartTime,PTFtemplate->fLower,
            (1.0/params->sampleRate),&lalDimensionlessUnit,
            3*numPoints/4 - numPoints/4);
      if (params->doTraceSNR)
        traceSNR = XLALCreateREAL4TimeSeries("traceSNR",
            &segStartTime,PTFtemplate->fLower,
            (1.0/params->sampleRate),&lalDimensionlessUnit,
            3*numPoints/4 - numPoints/4);
      if ( params->doBankVeto )
      {
        bankVeto = XLALCreateREAL4TimeSeries("bankVeto",
            &segStartTime,PTFtemplate->fLower,
            (1.0/params->sampleRate),&lalDimensionlessUnit,
            3*numPoints/4 - numPoints/4);
      }
      if ( params->doAutoVeto )
      {
        autoVeto = XLALCreateREAL4TimeSeries("autoVeto",
            &segStartTime,PTFtemplate->fLower,
            (1.0/params->sampleRate),&lalDimensionlessUnit,
            3*numPoints/4 - numPoints/4);
      }
      if ( params->doChiSquare )
      {
        chiSquare = XLALCreateREAL4TimeSeries("chiSquare",
            &segStartTime,PTFtemplate->fLower,
            (1.0/params->sampleRate),&lalDimensionlessUnit,
            3*numPoints/4 - numPoints/4);
      }

      /* Loop over ifos */
      for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if ( params->haveTrig[ifoNumber] )
        {
          /* Zero the storage vectors for the PTF filters */
          memset( PTFM[ifoNumber]->data, 0, 25 * sizeof(REAL8) );
          memset( PTFqVec[ifoNumber]->data, 0, 
                  5 * numPoints * sizeof(COMPLEX8) );

          /* Here (h|s) and (h|h) are calculated */
          cohPTFNormalize(fcTmplt,invspec[ifoNumber],PTFM[ifoNumber],NULL,
              PTFqVec[ifoNumber],&segments[ifoNumber]->sgmnt[j],invPlan,
              spinTemplate);

          // In this subroutine the overlap between template h and the 
          // bank veto template is calculated
          if ( params->doBankVeto )
          {
          for ( ui = 0 ; ui < subBankSize ; ui++ )
          {
            /*if ( spinBank )
            {
              memset(bankOverlaps[ui].PTFM[ifoNumber]->data,0,25*sizeof(COMPLEX8));
            }
            else
            {*/
            memset(bankOverlaps[ui].PTFM[ifoNumber]->data,0,1*sizeof(COMPLEX8));
//            }
            cohPTFComplexTemplateOverlaps(&(bankFcTmplts[ui]),fcTmplt,
                invspec[ifoNumber],0,
                bankOverlaps[ui].PTFM[ifoNumber]);
          }
          }
          // And if necessary the overlap between h(deltaT) and h
          if ( params->doAutoVeto )
          {
            autoVetoOverlaps(fcTmplt,autoTempOverlaps,invspec[ifoNumber],
                invPlan,0,params->numAutoPoints,timeStepPoints,ifoNumber);
          }

          verbose("Made filters for ifo %d,segment %d, template %d at %ld \n", 
              ifoNumber,j,i,time(NULL)-startTime);
        }
      }
      /* If necessary calculate the null stream filters */
      if ( params->doNullStream)
      {
        memset( PTFM[LAL_NUM_IFO]->data, 0, 25 * sizeof(REAL8) );
        memset( PTFqVec[LAL_NUM_IFO]->data, 0,
                5 * numPoints * sizeof(COMPLEX8) );
        cohPTFNormalize(fcTmplt,invspec[LAL_NUM_IFO],PTFM[LAL_NUM_IFO],NULL,
              PTFqVec[LAL_NUM_IFO],&segments[LAL_NUM_IFO]->sgmnt[j],invPlan,
              spinTemplate);
        verbose("Made filters for NULL stream,segmen %d, template %d at %ld\n",
              j,i,time(NULL)-startTime);
      }
      
      // This function calculates the cohSNR time series and all of the
      // signal based vetoes as appropriate
      cohPTFmodBasesUnconstrainedStatistic(cohSNR,PTFM,PTFqVec,params,
          spinTemplate,singleDetector,timeOffsets,Fplus,Fcross,j,pValues,
          gammaBeta,snrComps,nullSNR,traceSNR,bankVeto,autoVeto,chiSquare,
          subBankSize,bankOverlaps,bankNormOverlaps,dataOverlaps,
          autoTempOverlaps,fcTmplt,invspec,segments,invPlan);
     
      verbose("Made coherent statistic for segment %d, template %d at %ld \n",
          j,i,time(NULL)-startTime);      

      // This function adds any loud events to the list of triggers 
      eventId = cohPTFaddTriggers(params,&eventList,&thisEvent,cohSNR,*PTFtemplate,eventId,spinTemplate,singleDetector,pValues,gammaBeta,snrComps,nullSNR,traceSNR,bankVeto,autoVeto,chiSquare,PTFM);
      verbose("Generated triggers for segment %d, template %d at %ld \n",
          j,i,time(NULL)-startTime);
//      cohPTFclusterTriggers(params,&eventList,&thisEvent);
//      verbose("Clustered triggers for segment %d, template %d at %ld \n",
 //         j,i,time(NULL)-startTime);
      // Then we get a bunch of memory freeing statements
      for ( k = 0 ; k < 10 ; k++ )
      {
        if (pValues[k])
        {
            XLALDestroyREAL4TimeSeries(pValues[k]);
            pValues[k] = NULL;
        }
      }
      for ( k = 0; k < LAL_NUM_IFO; k++ )
      {
        if (snrComps[k])
        {
          XLALDestroyREAL4TimeSeries(snrComps[k]);
          snrComps[k] = NULL;
        }
      }
      if (gammaBeta[0]) XLALDestroyREAL4TimeSeries(gammaBeta[0]);
      if (gammaBeta[1]) XLALDestroyREAL4TimeSeries(gammaBeta[1]);
      if (nullSNR) XLALDestroyREAL4TimeSeries(nullSNR);
      if (traceSNR) XLALDestroyREAL4TimeSeries(traceSNR);
      if (bankVeto) XLALDestroyREAL4TimeSeries(bankVeto);
      if (autoVeto) XLALDestroyREAL4TimeSeries(autoVeto);
      if (chiSquare) XLALDestroyREAL4TimeSeries(chiSquare);
      XLALDestroyREAL4TimeSeries(cohSNR);
    }
    if ( params->doBankVeto )
    {
      for ( ui = 0 ; ui < subBankSize ; ui++ )
      {
        for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
        {
          if ( dataOverlaps[ui].PTFqVec[ifoNumber] )
            XLALDestroyCOMPLEX8VectorSequence( dataOverlaps[ui].PTFqVec[ifoNumber]);
          if ( bankOverlaps[ui].PTFM[ifoNumber] )
            XLALDestroyCOMPLEX8Array( bankOverlaps[ui].PTFM[ifoNumber]);
        }
      }
    }
  } // Main loop is ended here
  cohPTF_output_events_xml( params->outputFile, eventList, procpar, params );

  verbose("Generated output xml file, cleaning up and exiting at %ld \n",
      time(NULL)-startTime);

  // This function cleans up memory usage
  coh_PTF_cleanup(procpar,fwdplan,revplan,invPlan,channel,
      invspec,segments,eventList,PTFbankhead,fcTmplt,fcTmpltParams,
      fcInitParams,PTFM,PTFN,PTFqVec,Fplus,Fcross,timeOffsets);
  
  while ( PTFBankvetoHead )
  {
    InspiralTemplate *thisTmplt;
    thisTmplt = PTFBankvetoHead;
    PTFBankvetoHead = PTFBankvetoHead->next;
    if ( thisTmplt->event_id )
    {
      LALFree( thisTmplt->event_id );
    }
    LALFree( thisTmplt );
  }

  free_bank_veto_memory(bankNormOverlaps,PTFBankTemplates,bankFcTmplts,subBankSize,bankOverlaps,dataOverlaps);

  if ( autoTempOverlaps )
  {
    for (uj = 0; uj < params->numAutoPoints; uj++ )
    {
      for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if ( params->haveTrig[ifoNumber] )
        {
          if ( autoTempOverlaps[uj].PTFM[ifoNumber] )
          {
            XLALDestroyCOMPLEX8Array( autoTempOverlaps[uj].PTFM[ifoNumber] );
          }
        }
      }
    }
    LALFree( autoTempOverlaps );
  }


  verbose("Generated output xml file, cleaning up and exiting at %ld \n",
      time(NULL)-startTime);
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
  coh_PTF_params_inspiral_sanity_check( &params );
  params.programName = programName;
  params.cvsRevision = cvsRevision;
  params.cvsSource   = cvsSource;
  params.cvsDate     = cvsDate;
  return &params;
}

/* gets the data, performs any injections, and conditions the data */
static REAL4TimeSeries *coh_PTF_get_data( struct coh_PTF_params *params,\
               const char *ifoChannel, const char *dataCache, UINT4 ifoNumber  )
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
    invspec = compute_average_spectrum( channel,
        LALRINGDOWN_SPECTRUM_MEDIAN_MEAN, params->segmentDuration,
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

 /* TODO: trig start/end time condition */

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
  template->fLower = 38.;
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

void cohPTFmodBasesUnconstrainedStatistic(
    REAL4TimeSeries         *cohSNR,
    REAL8Array              *PTFM[LAL_NUM_IFO+1],
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
    struct coh_PTF_params   *params,
    UINT4                   spinTemplate,
    UINT4                   singleDetector,
    REAL8                   *timeOffsets,
    REAL8                   *Fplus,
    REAL8                   *Fcross,
    INT4                    segmentNumber,
    REAL4TimeSeries         *pValues[10],
    REAL4TimeSeries         *gammaBeta[2],
    REAL4TimeSeries         *snrComps[LAL_NUM_IFO],
    REAL4TimeSeries         *nullSNR,
    REAL4TimeSeries         *traceSNR,
    REAL4TimeSeries         *bankVeto,
    REAL4TimeSeries         *autoVeto,
    REAL4TimeSeries         *chiSquare,
    UINT4                   subBankSize,
    struct bankComplexTemplateOverlaps *bankOverlaps,
    struct bankTemplateOverlaps *bankNormOverlaps,
    struct bankDataOverlaps *dataOverlaps,
    struct bankComplexTemplateOverlaps *autoTempOverlaps,
    FindChirpTemplate       *fcTmplt,
    REAL4FrequencySeries    *invspec[LAL_NUM_IFO+1],
    RingDataSegments        *segments[LAL_NUM_IFO+1],
    COMPLEX8FFTPlan         *invPlan
)

{
// This function generates the SNR for every point in time and, where
// appropriate calculates the desired signal based vetoes.

  UINT4 i, j, k, m, vecLength, vecLengthTwo, UNUSED vecLengthSquare, UNUSED vecLengthTwoSquare;
  INT4 l;
  INT4 timeOffsetPoints[LAL_NUM_IFO];
  REAL4 deltaT = cohSNR->deltaT;
  UINT4 numPoints = floor( params->segmentDuration * params->sampleRate + 0.5 );

  // Code works slightly differently if spin/non spin and single/coherent
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

  // These arrays would store the maximized quantities if necessary
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

  FILE UNUSED *outfile;
  outfile = NULL;
/*  REAL8Array  *B, *Binv;*/
  REAL4 u1[vecLengthTwo],u2[vecLengthTwo],v1[vecLengthTwo],v2[vecLengthTwo];
  REAL4 *v1p,*v2p;
  REAL4 u1N[vecLength],u2N[vecLength],v1N[vecLength],v2N[vecLength];
  REAL4 v1_dot_u1, v1_dot_u2, v2_dot_u1, v2_dot_u2,max_eigen;
  REAL4 recSNR,traceSNRsq;
  REAL4 dAlpha,dBeta,dCee;
  REAL4 pValsTemp[vecLengthTwo];
  REAL4 betaGammaTemp[2];
  REAL4 a[LAL_NUM_IFO], b[LAL_NUM_IFO];

  gsl_matrix *BNull = gsl_matrix_alloc(vecLength,vecLength);
  gsl_matrix *B2Null = gsl_matrix_alloc(vecLength,vecLength);
  gsl_matrix *Bankeigenvecs[50];
  gsl_vector *Bankeigenvals[50];
  gsl_matrix *Autoeigenvecs = NULL;
  gsl_vector *Autoeigenvals = NULL;
  for (i = 0; i < 50; i++)
  {
    Bankeigenvecs[i] = NULL;
    Bankeigenvals[i] = NULL;  
  }
  gsl_permutation *p = gsl_permutation_alloc(vecLengthTwo);
  gsl_permutation *pNull = gsl_permutation_alloc(vecLength);
  gsl_eigen_symmv_workspace *matTempNull = gsl_eigen_symmv_alloc (vecLength);
  gsl_matrix *eigenvecs = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_vector *eigenvals = gsl_vector_alloc(vecLengthTwo);
  gsl_matrix *eigenvecsNull = gsl_matrix_alloc(vecLength,vecLength);
  gsl_vector *eigenvalsNull = gsl_vector_alloc(vecLength);

  // a = Fplus , b = Fcross
  for (i = 0; i < LAL_NUM_IFO; i++)
  {
    a[i] = Fplus[i];
    b[i] = Fcross[i];
  }

  // This function takes the (Q_i|Q_j) matrices, combines it across the ifos
  // and returns the eigenvalues and eigenvectors of this new matrix.
  // We later rotate and rescale the (Q_i|s) values such that in the new basis
  // this matrix will be the identity matrix.
  calculate_bmatrix(params,eigenvecs,eigenvals,a,b,PTFM,vecLength,vecLengthTwo,5);

  // If required also calculate these eigenvalues/vectors for the null stream 
  if ( params->doNullStream )
  {
    for (i = 0; i < vecLength; i++ )
    {
      for (j = 0; j < vecLength; j++ )
      {
        gsl_matrix_set(BNull,i,j,PTFM[LAL_NUM_IFO]->data[i*5+j]);
        gsl_matrix_set(B2Null,i,j,PTFM[LAL_NUM_IFO]->data[i*5+j]);
      }
    }
    gsl_eigen_symmv(B2Null,eigenvalsNull,eigenvecsNull,matTempNull); 
  }

  /*fprintf(stdout,"\n \n");*/

  /*for (i = 0; i < vecLengthTwo; i++ )
  {
    for (j = 0; j < vecLengthTwo; j++ )
    {
      fprintf(stdout,"%f ",gsl_matrix_get(eigenvecs,i,j));
    }
    fprintf(stdout,"\n");
  }

  fprintf(stdout,"\n \n");*/

  /*for (i = 0; i < vecLengthTwo; i++ )
  {
    fprintf(stdout,"%f ",gsl_vector_get(eigenvals,i));
  }

  fprintf(stdout,"\n \n");*/

  /* This loop takes the time offset in seconds and converts to time offset
  * in data points */
  for (i = 0; i < LAL_NUM_IFO; i++ )
  {
    timeOffsetPoints[i]=(int)(timeOffsets[i]/deltaT);
  }

  v1p = LALCalloc(vecLengthTwo , sizeof(REAL4));
  v2p = LALCalloc(vecLengthTwo , sizeof(REAL4));


  for ( i = numPoints/4; i < 3*numPoints/4; ++i ) /* Main loop over time */
  {
    // This function combines the various (Q_i | s) and rotates them into
    // the basis as discussed above.
    calculate_rotated_vectors(params,PTFqVec,v1p,v2p,a,b,timeOffsetPoints,
        eigenvecs,eigenvals,numPoints,i,vecLength,vecLengthTwo);

    /* Compute the dot products */
    v1_dot_u1 = v1_dot_u2 = v2_dot_u1 = v2_dot_u2 = max_eigen = 0.0;
    for (j = 0; j < vecLengthTwo; j++)
    {
      v1_dot_u1 += v1p[j] * v1p[j];
      v1_dot_u2 += v1p[j] * v2p[j];
      v2_dot_u2 += v2p[j] * v2p[j];
    }
    // And SNR is calculated
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

  UINT4 check;
  UINT4 chisqCheck = 0;
  REAL4 bestNR,snglSNRsq;
  INT4 numPointCheck = floor(params->timeWindow/cohSNR->deltaT + 0.5);
  struct bankCohTemplateOverlaps *bankCohOverlaps = NULL;
  struct bankCohTemplateOverlaps *autoCohOverlaps = NULL;
  struct bankDataOverlaps *chisqOverlaps = NULL;
  COMPLEX8VectorSequence *tempqVec = NULL;
  REAL4 *frequencyRangesPlus = NULL;
  REAL4 *frequencyRangesCross = NULL;
  REAL4 *powerBinsPlus = NULL;
  REAL4 *powerBinsCross = NULL;
  REAL4 fLow,fHigh;

  // Now we calculate all the extrinsic parameters and signal based vetoes
  // Only calculated if this will be a trigger

  for ( i = numPoints/4; i < 3*numPoints/4; ++i ) /* loop over time */
  {
    if (cohSNR->data->data[i-numPoints/4] > params->threshold)
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
      check = 1;
      if (check)
      {
        // The following block extracts the values of extrinsic parameters.
        // This follows the method set out in Diego's thesis to do this.
        calculate_rotated_vectors(params,PTFqVec,v1p,v2p,a,b,timeOffsetPoints,
        eigenvecs,eigenvals,numPoints,i,vecLength,vecLengthTwo);
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
        if ( spinTemplate == 1 )
        {
          dCee = (max_eigen - v1_dot_u1) / v1_dot_u2;
        }
        else
          dCee = 0;
        dAlpha = 1./(v1_dot_u1 + dCee * 2 * v1_dot_u2 + dCee*dCee*v2_dot_u2);
        dAlpha = pow(dAlpha,0.5);
        dBeta = dCee*dAlpha;
        // The p Values are calculated in the rotated frame
        for ( j = 0 ; j < vecLengthTwo ; j++ )
        {
          pValsTemp[j] = dAlpha*u1[j] + dBeta*u2[j];  
          pValues[j]->data->data[i - numPoints/4] = 0.;
        } 
        // This loop verifies that the SNR obtained is as before
        recSNR = 0;
        for ( j = 0 ; j < vecLengthTwo ; j++ )
        {
          for ( k = 0 ; k < vecLengthTwo ; k++ )
          {
            recSNR += pValsTemp[j]*pValsTemp[k] * (v1[j]*v1[k]+v2[j]*v2[k]);
          }
        }
        /*fprintf(stdout,"%e %e \n",max_eigen,recSNR);*/
        // Then we calculate the two phase/amplitude terms beta and gamma
        // These are explained in my PDF document
        betaGammaTemp[0] = 0;
        betaGammaTemp[1] = 0;
        for ( j = 0 ; j < vecLengthTwo ; j++ )
        {
          betaGammaTemp[0] += pValsTemp[j]*v1[j];
          betaGammaTemp[1] += pValsTemp[j]*v2[j];
        }
        gammaBeta[0]->data->data[i - numPoints/4] = betaGammaTemp[0];
        gammaBeta[1]->data->data[i - numPoints/4] = betaGammaTemp[1];
  
        // The p Values need to be rotated back into the original frame.
        // Currently we are recording values in rotated frame
        for ( j = 0 ; j < vecLengthTwo ; j++ )
        {
          for ( k = 0 ; k < vecLengthTwo ; k++ )
          {
            pValues[j]->data->data[i-numPoints/4]+=gsl_matrix_get(eigenvecs,j,k)*pValsTemp[k];
          }
        }
 
        // And we check that this still gives the expected SNR in the
        // unrotated basis.
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
          
        }

        // First sbv to be calculated is the null SNR.
        if (params->doNullStream)
        {
          // As with SNR do the rotation and calculate SNR.
          for ( j = 0; j < vecLength; j++ ) /* Construct the vi vectors */
          {
            v1N[j] = PTFqVec[LAL_NUM_IFO]->data[j*numPoints+i].re;
            v2N[j] = PTFqVec[LAL_NUM_IFO]->data[j*numPoints+i].im;
          }
     
          for ( j = 0 ; j < vecLength ; j++ )
          {
            u1N[j] = 0.;
            u2N[j] = 0.;
            for ( k = 0 ; k < vecLength ; k++ )
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
            max_eigen = 0.5 * ( v1_dot_u1 + v2_dot_u2 );
          }
          else
          {
            max_eigen = 0.5*(v1_dot_u1+v2_dot_u2+sqrt((v1_dot_u1-v2_dot_u2)
                * (v1_dot_u1 - v2_dot_u2) + 4 * v1_dot_u2 * v1_dot_u2 ));
          }
          nullSNR->data->data[i-numPoints/4] = sqrt(max_eigen);
        }

        // Next up is Trace SNR and the SNR components
        if (params->doTraceSNR)
        {
          // traceSNR is calculated as normal SNR but cross terms are not added
          traceSNRsq = 0;
          for( k = 0; k < LAL_NUM_IFO; k++)
          {
            if ( params->haveTrig[k] )
            {
              for ( j = 0; j < vecLengthTwo ; j++ )
              {
                if (j < vecLength)
                {
                  v1[j] = a[k] * PTFqVec[k]->data[j*numPoints+i+timeOffsetPoints[k]].re;
                  v2[j] = a[k] * PTFqVec[k]->data[j*numPoints+i+timeOffsetPoints[k]].im;
                }
                else
                {
                  v1[j] = b[k] * PTFqVec[k]->data[(j-vecLength)*numPoints+i+timeOffsetPoints[k]].re;
                  v2[j] = b[k] * PTFqVec[k]->data[(j-vecLength)*numPoints+i+timeOffsetPoints[k]].im;
                }
              }
              for ( j = 0 ; j < vecLengthTwo ; j++ )
              {
                u1[j] = 0.;
                u2[j] = 0.;
                for ( m = 0 ; m < vecLengthTwo ; m++ )
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
                max_eigen = 0.5 * ( v1_dot_u1 + v2_dot_u2 );
              }
              else
              {
                max_eigen = 0.5 * ( v1_dot_u1 + v2_dot_u2 + sqrt( (v1_dot_u1 - v2_dot_u2)
                * (v1_dot_u1 - v2_dot_u2) + 4 * v1_dot_u2 * v1_dot_u2 ));
              }
              // This needs to be converted for spinning case!
              snglSNRsq = v1[0]*v1[0] + v1[1]*v1[1];
              snglSNRsq = snglSNRsq/(a[k]*a[k]*PTFM[k]->data[0]);
              traceSNRsq += max_eigen;
              snrComps[k]->data->data[i-numPoints/4]=sqrt(snglSNRsq);
            }
          }
          traceSNR->data->data[i-numPoints/4] = sqrt(traceSNRsq);
        }

        // Next is the bank veto
        if ( params->doBankVeto )
        {
          if (! singleDetector)
          {
            for ( j = 0 ; j < subBankSize+1 ; j++ )
            {
              if (! Bankeigenvecs[j] )
              {
                Bankeigenvecs[j] = gsl_matrix_alloc(2,2);
                Bankeigenvals[j] = gsl_vector_alloc(2);
                // Here we calculate the eigenvectors for each bank template
                if (j == subBankSize)
                {
                  calculate_bmatrix(params,Bankeigenvecs[j],Bankeigenvals[j],
                      a,b,PTFM,1,2,5);
                }
                else
                {
                  calculate_bmatrix(params,Bankeigenvecs[j],Bankeigenvals[j],
                      a,b,bankNormOverlaps[j].PTFM,1,2,1);
                }
              }
            }

            if (! bankCohOverlaps)
            {
              bankCohOverlaps = LALCalloc(subBankSize,sizeof(*bankCohOverlaps));
              for ( j = 0 ; j < subBankSize; j++ )
              {
                bankCohOverlaps[j].rotReOverlaps = gsl_matrix_alloc(2,2);
                bankCohOverlaps[j].rotImOverlaps = gsl_matrix_alloc(2,2);
                // We calculate the coherent overlaps in this function
                calculate_coherent_bank_overlaps(params,bankOverlaps[j],
                    bankCohOverlaps[j],a,b,Bankeigenvecs[subBankSize],
                    Bankeigenvals[subBankSize],Bankeigenvecs[j],
                    Bankeigenvals[j]);
              }
            }
            // In this function all the filters are combined to produce the
            // value of the bank veto.
            bankVeto->data->data[i-numPoints/4] = calculate_bank_veto_max_phase_coherent(numPoints,i,subBankSize,a,b,params,bankCohOverlaps,dataOverlaps,PTFqVec,timeOffsetPoints,Bankeigenvecs,Bankeigenvals);
          }
//          autoVeto->data->data[i-numPoints/4] =calculate_bank_veto(numPoints,i,subBankSize,vecLength,a,b,cohSNR->data->data[i-numPoints/4],PTFM,params,bankOverlaps,bankNormOverlaps,dataOverlaps,pValues,gammaBeta,PTFqVec,timeOffsetPoints,singleDetector);
          // The single detector function is a little messy at the moment
          if (singleDetector)
            bankVeto->data->data[i-numPoints/4] =calculate_bank_veto_max_phase(numPoints,i,subBankSize,PTFM,params,bankOverlaps,bankNormOverlaps,dataOverlaps,PTFqVec,timeOffsetPoints);
        }   

        // Now we do the auto veto
        if ( params->doAutoVeto )
        {
          if (! singleDetector)
          {
            if (! Autoeigenvecs )
            {
              Autoeigenvecs = gsl_matrix_alloc(2,2);
              Autoeigenvals = gsl_vector_alloc(2);
              // Again the eigenvectors/values are calculated
              calculate_bmatrix(params,Autoeigenvecs,Autoeigenvals,
                  a,b,PTFM,1,2,5);
//              fprintf(stderr,"Eigenvecs %e %e %e %e\n",gsl_matrix_get(Autoeigenvecs,0,0),gsl_matrix_get(Autoeigenvecs,1,0),gsl_matrix_get(Autoeigenvecs,0,1),gsl_matrix_get(Autoeigenvecs,1,1));
//              fprintf(stderr,"Eigenvals %e %e \n",gsl_vector_get(Autoeigenvals,0),gsl_vector_get(Autoeigenvals,1));
            }

            if (! autoCohOverlaps)
            {
              autoCohOverlaps = LALCalloc(params->numAutoPoints,sizeof(*autoCohOverlaps));
              for ( j = 0 ; j < params->numAutoPoints; j++ )
              {
                autoCohOverlaps[j].rotReOverlaps = gsl_matrix_alloc(2,2);
                autoCohOverlaps[j].rotImOverlaps = gsl_matrix_alloc(2,2);
                // The coherent rotated overlaps are calculated
                calculate_coherent_bank_overlaps(params,autoTempOverlaps[j],
                    autoCohOverlaps[j],a,b,Autoeigenvecs,
                    Autoeigenvals,Autoeigenvecs,Autoeigenvals);
              }
            }
          }
          // Auto veto is calculated
          autoVeto->data->data[i-numPoints/4] = calculate_auto_veto_max_phase_coherent(numPoints,i,a,b,params,autoCohOverlaps,PTFqVec,timeOffsetPoints,Autoeigenvecs,Autoeigenvals);
        }
      }
    }
  }
  /* To save memory we cut the loop here, clean the memory before calculating
  chi square */
  LALFree(v1p);
  LALFree(v2p);

  if (params->doBankVeto)
  {
    for ( j = 0 ; j < subBankSize+1 ; j++ )
    {
      if (Bankeigenvecs[j])
        gsl_matrix_free(Bankeigenvecs[j]);
      if (Bankeigenvals[j])
        gsl_vector_free(Bankeigenvals[j]);
    }
    if (bankCohOverlaps)
    {
      for ( j = 0 ; j < subBankSize ; j++ )
      {
        gsl_matrix_free(bankCohOverlaps[j].rotReOverlaps);
        gsl_matrix_free(bankCohOverlaps[j].rotImOverlaps);
      }
      LALFree(bankCohOverlaps);
    }
  }

  if (params->doAutoVeto)
  {
    if ( Autoeigenvecs )
      gsl_matrix_free( Autoeigenvecs );
      Autoeigenvecs = NULL;
    if ( Autoeigenvals )
      gsl_vector_free( Autoeigenvals );
      Autoeigenvals = NULL;
    if (autoCohOverlaps)
    {
      for ( j = 0 ; j < params->numAutoPoints ; j++ )
      {
        gsl_matrix_free(autoCohOverlaps[j].rotReOverlaps);
        gsl_matrix_free(autoCohOverlaps[j].rotImOverlaps);
      }
      LALFree(autoCohOverlaps);
    }
  }

  /* And do the loop again to calculate chi square */

  for ( i = numPoints/4; i < 3*numPoints/4; ++i ) /* loop over time */
  {
    if (cohSNR->data->data[i-numPoints/4] > params->threshold)
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
        if ( params->chiSquareCalcThreshold )
        {
          chisqCheck = 1;
          
          bestNR = cohSNR->data->data[i-numPoints/4];

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
  
          if (params->doBankVeto)
          {
            if (bankVeto->data->data[i-numPoints/4] > 40)
              bestNR = bestNR/pow(( 1 + pow(bankVeto->data->data[i-numPoints/4]/((REAL4)subBankSize*4.),params->bankVetoq/params->bankVeton))/2.,1./params->bankVetoq);
            if (bestNR < params->chiSquareCalcThreshold)
              chisqCheck = 0;
          }

          bestNR = cohSNR->data->data[i-numPoints/4];

          if (params->doAutoVeto)
          {
            if (autoVeto->data->data[i-numPoints/4] > 40)
              bestNR = bestNR/pow(( 1 + pow(autoVeto->data->data[i-numPoints/4]/((REAL4)params->numAutoPoints*4.),params->autoVetoq/params->autoVeton))/2.,1./params->autoVetoq);
            if (bestNR < params->chiSquareCalcThreshold)
              chisqCheck = 0;
          } 
        }
        else
          chisqCheck = 1;

        if (params->doChiSquare && chisqCheck )
        {
          if (! Autoeigenvecs )
          {
            Autoeigenvecs = gsl_matrix_alloc(2,2);
            Autoeigenvals = gsl_vector_alloc(2);
            // Again the eigenvectors/values are calculated
            calculate_bmatrix(params,Autoeigenvecs,Autoeigenvals,
                a,b,PTFM,1,2,5);
          }
          if (! frequencyRangesPlus)
          {
            frequencyRangesPlus = (REAL4 *)
              LALCalloc( params->numChiSquareBins-1, sizeof(REAL4) );
            frequencyRangesCross = (REAL4 *)
              LALCalloc( params->numChiSquareBins-1, sizeof(REAL4) );
            powerBinsPlus = (REAL4 *) 
              LALCalloc( params->numChiSquareBins, sizeof(REAL4) );
            powerBinsCross = (REAL4 *)
              LALCalloc( params->numChiSquareBins, sizeof(REAL4) );
            calculate_standard_chisq_freq_ranges(params,fcTmplt,invspec,PTFM,a,b,frequencyRangesPlus,frequencyRangesCross,Autoeigenvecs);
            calculate_standard_chisq_power_bins(params,fcTmplt,invspec,PTFM,a,b,frequencyRangesPlus,powerBinsPlus,powerBinsCross,Autoeigenvecs);
          }
          if (! tempqVec)
            tempqVec = XLALCreateCOMPLEX8VectorSequence ( 1, numPoints );
          if (! chisqOverlaps)
          {
            chisqOverlaps = LALCalloc(params->numChiSquareBins,sizeof( *chisqOverlaps));
            for( j = 0; j < params->numChiSquareBins; j++)
            {
              if (params->numChiSquareBins == 1)
              {
                fLow = 0;
                fHigh = 0;
              }
              else if (j == 0)
              {
                fLow = 0;
                fHigh = frequencyRangesPlus[0];
              }
              else if (j == params->numChiSquareBins-1)
              {
                fLow = frequencyRangesPlus[params->numChiSquareBins-2];
                fHigh = 0;
              }
              else
              {
                fLow = frequencyRangesPlus[j-1];
                fHigh = frequencyRangesPlus[j];
              }                 
              for( k = 0; k < LAL_NUM_IFO; k++)
              {
                if ( params->haveTrig[k] )
                {
                  chisqOverlaps[j].PTFqVec[k] =
                      XLALCreateCOMPLEX8VectorSequence ( 1,
                      3*numPoints/4 - numPoints/4 + 10000);
                  cohPTFBankFilters(fcTmplt,0,
                  &segments[k]->sgmnt[segmentNumber],invPlan,tempqVec,
                  chisqOverlaps[j].PTFqVec[k],fLow,fHigh);
                }
                else
                  chisqOverlaps[j].PTFqVec[k] = NULL;
              }
            }
          }
          /* Calculate chi square here */
          chiSquare->data->data[i-numPoints/4] = calculate_chi_square(params,numPoints,i,chisqOverlaps,PTFqVec,a,b,timeOffsetPoints,Autoeigenvecs,Autoeigenvals,powerBinsPlus,powerBinsCross);
        
        }
        else if (params->doChiSquare)
          chiSquare->data->data[i-numPoints/4] = 0;
      }
    }
  }

  if (params->doChiSquare)
  {
    if ( Autoeigenvecs )
      gsl_matrix_free( Autoeigenvecs );
    if ( Autoeigenvals )
      gsl_vector_free( Autoeigenvals );
    if (frequencyRangesPlus)
      LALFree(frequencyRangesPlus);
    if (frequencyRangesCross)
      LALFree(frequencyRangesCross);
    if (tempqVec)
      XLALDestroyCOMPLEX8VectorSequence( tempqVec );
    if (chisqOverlaps)
    {
      for( j = 0; j < params->numChiSquareBins; j++)
      {
        for( k = 0; k < LAL_NUM_IFO; k++)
        {
          if (chisqOverlaps[j].PTFqVec[k])
          {
            XLALDestroyCOMPLEX8VectorSequence(chisqOverlaps[j].PTFqVec[k]);
          }
        }
      }
      LALFree(chisqOverlaps);      
    }

  }

/*  outfile = fopen("cohSNR_timeseries.dat","w");
  for ( i = 0; i < cohSNR->data->length; ++i)
  {
    fprintf (outfile,"%f %f \n",deltaT*i,cohSNR->data->data[i]);
  }
  fclose(outfile);*/


/*  outfile = fopen("bank_veto_timeseries.dat","w");
  for ( i = 0; i < bankVeto->data->length; ++i)
  {
    fprintf (outfile,"%f %f \n",deltaT*i,bankVeto->data->data[i]);
  }
  fclose(outfile);*/

/*  outfile = fopen("auto_veto_timeseries.dat","w");
  for ( i = 0; i < autoVeto->data->length; ++i)
  {
    fprintf (outfile,"%f %f \n",deltaT*i,autoVeto->data->data[i]);
  }
  fclose(outfile);*/

/*  outfile = fopen("chi_square_timeseries.dat","w");
  for ( i = 0; i < chiSquare->data->length; ++i)
  {
    fprintf (outfile,"%f %f \n",deltaT*i,chiSquare->data->data[i]);
  }
  fclose(outfile);*/

  
/*  outfile = fopen("nullSNR_timeseries.dat","w");
  for ( i = 0; i < nullSNR->data->length; ++i)
  {
    fprintf (outfile,"%f %f \n",deltaT*i,nullSNR->data->data[i]);
  }
  fclose(outfile);*/
  

  
/*  outfile = fopen("traceSNR_timeseries.dat","w");
  for ( i = 0; i < traceSNR->data->length; ++i)
  {
    fprintf (outfile,"%f %f \n",deltaT*i,traceSNR->data->data[i]);
  }
  fclose(outfile);*/
  
/*  outfile = fopen("h1SNR_timeseries.dat","w");
  for ( i = 0; i < snrComps[1]->data->length; ++i)
  {
    fprintf (outfile,"%f %f \n",deltaT*i,snrComps[1]->data->data[i]);
  }
  fclose(outfile);*/

/*  outfile = fopen("l1SNR_timeseries.dat","w");
  for ( i = 0; i < snrComps[3]->data->length; ++i)
  {
    fprintf (outfile,"%f %f \n",deltaT*i,snrComps[3]->data->data[i]);
  }
  fclose(outfile);*/

/*  outfile = fopen("v1SNR_timeseries.dat","w");
  for ( i = 0; i < snrComps[5]->data->length; ++i)
  {
    fprintf (outfile,"%f %f \n",deltaT*i,snrComps[5]->data->data[i]);
  }
  fclose(outfile);*/

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

  gsl_matrix_free(BNull);
  gsl_matrix_free(B2Null);
  gsl_permutation_free(p);
  gsl_permutation_free(pNull);
  gsl_eigen_symmv_free(matTempNull);
  gsl_matrix_free(eigenvecs);
  gsl_vector_free(eigenvals);
  gsl_matrix_free(eigenvecsNull);
  gsl_vector_free(eigenvalsNull);


}

UINT8 cohPTFaddTriggers(
    struct coh_PTF_params   *params,
    MultiInspiralTable      **eventList,
    MultiInspiralTable      **thisEvent,
    REAL4TimeSeries         *cohSNR,
    InspiralTemplate        PTFTemplate,
    UINT8                   eventId,
    UINT4                   spinTrigger,
    UINT4                   singleDetector,
    REAL4TimeSeries         *pValues[10],
    REAL4TimeSeries         *gammaBeta[2],
    REAL4TimeSeries         *snrComps[LAL_NUM_IFO],
    REAL4TimeSeries         *nullSNR,
    REAL4TimeSeries         *traceSNR,
    REAL4TimeSeries         *bankVeto,
    REAL4TimeSeries         *autoVeto,
    REAL4TimeSeries         *chiSquare,
    REAL8Array              *PTFM[LAL_NUM_IFO+1]
)
{
  // This function adds a trigger to the event list

  UINT4 i;
  INT4 j;
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
        if (params->doNullStream)
          currEvent->null_statistic = nullSNR->data->data[i];
        if (params->doTraceSNR)
          currEvent->null_stat_degen = traceSNR->data->data[i];
        if (params->doBankVeto)
        {
          currEvent->bank_chisq = bankVeto->data->data[i];
          currEvent->bank_chisq_dof = params->BVsubBankSize;
        }
        if (params->doAutoVeto)
        {
          currEvent->cont_chisq = autoVeto->data->data[i];
//          fprintf(stderr, "Auto Veto %e \n",currEvent->cont_chisq);
          currEvent->cont_chisq_dof = params->numAutoPoints;
        }
        if (params->doChiSquare)
        {
          currEvent->chisq = chiSquare->data->data[i];
        }
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
        if (snrComps[LAL_IFO_G1])
        {
          currEvent->chisq_g = snrComps[LAL_IFO_G1]->data->data[i];
          currEvent->sigmasq_g = PTFM[LAL_IFO_G1]->data[0];
        }
        if (snrComps[LAL_IFO_H1])
        {
          currEvent->chisq_h1 = snrComps[LAL_IFO_H1]->data->data[i];
          currEvent->sigmasq_h1 = PTFM[LAL_IFO_H1]->data[0];
        }
        if (snrComps[LAL_IFO_H2])
        {
          currEvent->chisq_h2 = snrComps[LAL_IFO_H2]->data->data[i];
          currEvent->sigmasq_h2 = PTFM[LAL_IFO_H2]->data[0];
        }
        if (snrComps[LAL_IFO_L1])
        {
          currEvent->chisq_l = snrComps[LAL_IFO_L1]->data->data[i];
          currEvent->sigmasq_l = PTFM[LAL_IFO_L1]->data[0];
        }
        if (snrComps[LAL_IFO_T1])
        {
          currEvent->chisq_t = snrComps[LAL_IFO_T1]->data->data[i];
          currEvent->sigmasq_t = PTFM[LAL_IFO_T1]->data[0];
        }
        if (snrComps[LAL_IFO_V1])
        {
          currEvent->chisq_v = snrComps[LAL_IFO_V1]->data->data[i];
          currEvent->sigmasq_v = PTFM[LAL_IFO_V1]->data[0];
        }
        if (spinTrigger == 1)
        {
          if (singleDetector == 1)
            currEvent->snr_dof = 6;
          else
            currEvent->snr_dof = 12;
        }
        else
        {
          if (singleDetector == 1)
            currEvent->snr_dof = 2;
          else
            currEvent->snr_dof = 4;
        }
      }
    }
  }
  *thisEvent = currEvent;
  return eventId;
}

void cohPTFclusterTriggers(
  MultiInspiralTable      **eventList,
  MultiInspiralTable      **thisEvent
)
{
  MultiInspiralTable *currEvent = *eventList;
  MultiInspiralTable *currEvent2 = NULL;
  MultiInspiralTable *newEvent = NULL;
  MultiInspiralTable *newEventHead = NULL;
  LIGOTimeGPS time1,time2;
  UINT4 rejectTrigger;
  UINT4 triggerNum = 0;
  UINT4 lenTriggers = 0;
  UINT4 numRemovedTriggers = 0;

  while (currEvent)
  {
    lenTriggers+=1;
    currEvent = currEvent->next;
  }

  currEvent = *eventList;
  UINT4 rejectTriggers[lenTriggers];

  while (currEvent)
  {
    rejectTrigger = 0;
    time1.gpsSeconds=currEvent->end_time.gpsSeconds;
    time1.gpsNanoSeconds = currEvent->end_time.gpsNanoSeconds;
    currEvent2 = *eventList;
    while (currEvent2)
    {
      time2.gpsSeconds=currEvent2->end_time.gpsSeconds;
      time2.gpsNanoSeconds=currEvent2->end_time.gpsNanoSeconds;
      if (fabs(XLALGPSDiff(&time1,&time2)) < 0.1)
      {
        if (currEvent->snr < currEvent2->snr && (currEvent->event_id->id != currEvent2->event_id->id))
        {
          rejectTrigger = 1;
          numRemovedTriggers +=1;
          break;
        }
        else
        {
          currEvent2 = currEvent2->next;
        }
      }
      else
      {
        currEvent2 = currEvent2->next;
      }
    }
    rejectTriggers[triggerNum] = rejectTrigger;
    triggerNum += 1;
    currEvent = currEvent->next;
  }

  currEvent = *eventList;
  triggerNum = 0;

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
      if ( currEvent->event_id )
      {
        LALFree( currEvent->event_id );
      }
      currEvent2 = currEvent->next;
      LALFree( currEvent );  
      currEvent = currEvent2;
    }
    triggerNum+=1;
  }
  if (newEvent)
  {
    newEvent->next = NULL;
    *eventList = newEventHead;
    *thisEvent = newEvent;
  }
} 
        
static void coh_PTF_cleanup(
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
    REAL8                   *Fcross
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
  if ( timeOffsets )
    LALFree( timeOffsets );
  if ( Fplus )
    LALFree( Fplus );
  if ( Fcross )
    LALFree( Fcross );
}

  

