/********************************** <lalVerbatim file="PowerInitSearchCV">
Author: Brady, P and Brown, D
$Id$
**************************************************** </lalVerbatim> */

#include <lal/LALStdlib.h>
#include <lal/ExcessPower.h>
#include <lal/TFTransform.h>
#include <lal/EPData.h>
#include "InitSearch.h"

NRCSID( INITSEARCHC, "power $Id$");

void LALInitSearch(
    LALStatus             *status,
    void                 **searchParams,
    LALInitSearchParams   *initSearchParams
    )
{
  
  EPSearchParams *params;


  INITSTATUS (status, "LALInitSearch", INITSEARCHC);
  ATTATCHSTATUSPTR (status);

  /*
   *
   * sleep to give me time to attach ddd
   *
   */

#ifdef POWER_SO_ATTACH_DDD
  {
   unsigned int doze = 10;
   pid_t myPID;
   myPID = getpid( );
   fprintf( stdout, "pid %d sleeping for %d seconds\n", myPID, doze );
   fflush( stdout );
   sleep( doze );
   fprintf( stdout, "pid %d awake\n", myPID );
   fflush( stdout );
  }
#endif

  /*
   *
   * Sanity checks on the input, output, params
   *
   */

  if ( !(initSearchParams) ){
    ABORT( status, INITSEARCHH_ENULL, INITSEARCHH_MSGENULL);
  }
  if ( (*searchParams) ){
    ABORT( status, INITSEARCHH_ENNUL, INITSEARCHH_MSGENNUL);
  }

  if ( !(initSearchParams->argv) ){
    ABORT( status, INITSEARCHH_ENULL, INITSEARCHH_MSGENULL);
  }
  if ( initSearchParams->argc < 16 || initSearchParams->argc > 17 ){
        ABORT( status, INITSEARCHH_EARGS, INITSEARCHH_MSGEARGS ); 
  }
 
  if (strcmp( initSearchParams->argv[0], "-filterparams" )){
    ABORT( status, INITSEARCHH_EARGS, INITSEARCHH_MSGEARGS);
  }
  if ( atoi( initSearchParams->argv[1] ) <= 0 ){
    ABORT( status, INITSEARCHH_ENUMZ, INITSEARCHH_MSGENUMZ);
  }
  if ( atoi( initSearchParams->argv[2] ) <= 0 ){
    ABORT( status, INITSEARCHH_ESEGZ, INITSEARCHH_MSGESEGZ);
  }
  if ( atoi( initSearchParams->argv[3] ) < 0 ){
    ABORT( status, INITSEARCHH_EOVLP, INITSEARCHH_MSGEOVLP);
  }
  if ( atoi( initSearchParams->argv[4] ) < 0 ){
    ABORT( status, INITSEARCHH_EOVPF, INITSEARCHH_MSGEOVPF);
  }
  if ( atoi( initSearchParams->argv[5] ) <= 0 ){
    ABORT( status, INITSEARCHH_EMFBZ, INITSEARCHH_MSGEMFBZ);
  }
  if ( atoi( initSearchParams->argv[6] ) <= 0 ){
    ABORT( status, INITSEARCHH_EMTBZ, INITSEARCHH_MSGEMTBZ);
  }
  if ( atof( initSearchParams->argv[7] ) <= 0.0 ){
    ABORT( status, INITSEARCHH_EFLOW, INITSEARCHH_MSGEFLOW);
  }
  if ( atof( initSearchParams->argv[8] ) <= 0.0 ){
    ABORT( status, INITSEARCHH_EDELF, INITSEARCHH_MSGEDELF);
  }
  if ( atoi( initSearchParams->argv[9] ) <= 0 ){
    ABORT( status, INITSEARCHH_ELTFZ, INITSEARCHH_MSGELTFZ);
  }
  if ( atof( initSearchParams->argv[10] ) <= 1.0 ){
    ABORT( status, INITSEARCHH_ESIGM, INITSEARCHH_MSGESIGM);
  }
  if ( atof( initSearchParams->argv[11] ) <= 0.0 || 
       atof( initSearchParams->argv[11] ) >= 1.0 ){
    ABORT( status, INITSEARCHH_EALPH, INITSEARCHH_MSGEALPH);
  }
  if ( atoi( initSearchParams->argv[12] ) < 1 ){
    ABORT( status, INITSEARCHH_EDUTY, INITSEARCHH_MSGEDUTY);
  }
  if ( atof( initSearchParams->argv[13] ) < 0 ){
    ABORT( status, INITSEARCHH_EAMAX, INITSEARCHH_MSGEAMAX);
  }
  /* EK - modified, argv[14] is the number of communicated events (integer) */
  if ( atoi( initSearchParams->argv[14] ) < 1 ||
      atoi( initSearchParams->argv[14] ) > 9999){
    ABORT( status, INITSEARCHH_EE2MS, INITSEARCHH_MSGEE2MS);
  }
  /* EK - modified, argv[15] is the string representing the channel to apply */
  if ( strlen( initSearchParams->argv[15] ) == 0 ){
    ABORT( status, INITSEARCHH_ECHNL, INITSEARCHH_MSGECHNL);
  }
  
  /*
   *
   * allocate memory 
   *
   */

  *searchParams = params = LALMalloc (sizeof( EPSearchParams )); 
  if ( !params )
  {
    ABORT (status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
  }

  params->tfTilingInput = 
    (CreateTFTilingIn *) LALMalloc (sizeof(CreateTFTilingIn));
  if ( !params->tfTilingInput )
  {
    LALFree (*searchParams); *searchParams = NULL;
    ABORT (status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
  }

  params->initParams = (EPInitParams *) LALMalloc (sizeof(EPInitParams));
  if ( !params-> initParams)
  {
    LALFree (params->tfTilingInput); params->tfTilingInput = NULL;
    LALFree (*searchParams); *searchParams = NULL;
    ABORT (status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
  }

  params->compEPInput = 
    (ComputeExcessPowerIn *) LALMalloc (sizeof(ComputeExcessPowerIn));
  if ( !params-> compEPInput)
  {
    LALFree (params->initParams); params->initParams = NULL;
    LALFree (params->tfTilingInput); params->tfTilingInput = NULL;
    LALFree (*searchParams); *searchParams = NULL;
    ABORT (status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
  }

   /* EK - Channel name to process (e.g. "ifodmro") */
   params->channelName = (CHAR *) LALMalloc(strlen( initSearchParams->argv[15] )+1 );
   if ( !params->channelName )
   {
     LALFree (params->compEPInput); params->compEPInput = NULL;
     LALFree (params->initParams); params->initParams = NULL;
     LALFree (params->tfTilingInput); params->tfTilingInput = NULL;
    ABORT (status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
  }     


  /*
   * 
   * assign the initial parameter values
   *
   */

  /* Number of data points in a segment */
  params->initParams->numPoints         = atoi( initSearchParams->argv[1] );
  /* Number of overlapping data segments */
  params->initParams->numSegments       = atoi( initSearchParams->argv[2] );
  /* Number of segments sent to slave */
  params->initParams->segDutyCycle      = atoi( initSearchParams->argv[12] ); 
  /* Overlap betweeen segments (# of points) */
  params->ovrlap                        = atoi( initSearchParams->argv[3] );  
  /* Identify events with alpha less that this value */
  params->alphaThreshold                = atof( initSearchParams->argv[13] ); 
  /* Number of data points in a segment */
  params->ntotT                         = params->initParams->numPoints;      
  
  /* Amount of overlap between neighboring TF tiles */
  params->tfTilingInput->overlapFactor  = atoi( initSearchParams->argv[4] );  
  /* Smallest extent in freq of TF tiles to search */
  params->tfTilingInput->minFreqBins    = atoi( initSearchParams->argv[5] );  
  /* Smallest extent in time of TF tiles to search */
  params->tfTilingInput->minTimeBins    = atoi( initSearchParams->argv[6] );  
  /* Lowest frequency in Hz to be searched */
  params->tfTilingInput->flow           = atof( initSearchParams->argv[7] );  
  /* Frequency resolution of first TF plane */
  params->tfTilingInput->deltaF         = atof( initSearchParams->argv[8] );  
  /* Length (N_F) of first TF plane (with N_T = 1) */
  params->tfTilingInput->length         = atoi( initSearchParams->argv[9] );  
  /* threshold number of sigma */
  params->compEPInput->numSigmaMin      = atof( initSearchParams->argv[10] ); 
  /* default alpha value for tiles with sigma < numSigmaMin */
  params->compEPInput->alphaDefault     = atof( initSearchParams->argv[11] ); 
  /* EK - Max. number of events to communicate to master (default: 16) */
  params->events2Master                 = atoi( initSearchParams->argv[14] );
  /* EK - Channel name to process (e.g. "ifodmro") */
  strcpy( params->channelName, initSearchParams->argv[15] );
  /* Simulation type */
  if ( initSearchParams->argc > 16 ){
    if ( atoi( initSearchParams->argv[16] ) < 0 ||
         atoi( initSearchParams->argv[16] ) > 2 ){
      ABORT( status, INITSEARCHH_ESIM, INITSEARCHH_MSGESIM);
    }
    else
    {
      params->simType                       =  atoi( initSearchParams->argv[16] );
    }
  }
  else
  {
    /* default is no simulation if not provided in filterparams */
    params->simType = 0;
  }

  /* initialize parameter structures */
  params->tfTiling     = NULL;
  params->epSegVec     = NULL;
  params->numSlaves    = NULL;

  /* allocate memory for the conditioned segments */
  LALCreateEPDataSegmentVector (status->statusPtr, &(params->epSegVec), 
      params->initParams);  
  CHECKSTATUSPTR (status);

  /* initialize parameters */
  params->searchMaster    = initSearchParams->nodeClass;
  params->haveData        = 0;
  params->currentSegment  = 0;
  params->winParams.type = Bartlett;
  params->initParams->method = useMean;

  DETATCHSTATUSPTR (status);
  RETURN (status);
}

