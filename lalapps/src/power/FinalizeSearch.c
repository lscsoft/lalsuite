/********************************** <lalVerbatim file="PowerFinalizeSearchCV">
Author: Brady, P
$Id$
**************************************************** </lalVerbatim> */

#include <lal/LALStdlib.h>
#include <lal/ExcessPower.h>
#include <lal/TFTransform.h>
#include <lal/EPData.h>
#include "FinalizeSearch.h"

NRCSID( FINALIZESEARCHC, "power $Id$" );

void
LALFinalizeSearch(
    LALStatus             *status,
    void                 **searchParams
    )
{
  EPSearchParams *params;

  INITSTATUS (status, "LALFinalizeSearch", FINALIZESEARCHC);
  ATTATCHSTATUSPTR (status);

  /* check searchParams exists and cast approriately */
  if ( !(*searchParams) ){
    ABORT( status, FINALIZESEARCHH_ENULL, FINALIZESEARCHH_MSGENULL );
  }
  params = (EPSearchParams *) *searchParams;

  /* destroy memory for conditioned segments */
  LALDestroyEPDataSegmentVector(status->statusPtr, &(params->epSegVec));
  
  if ( params->searchMaster )
  {
    /* free numSlaves */
    LALFree (params->numSlaves);
  }
  
  /* free compEPInput */
  LALFree (params->compEPInput);

  /* free EPInitParams */
  LALFree (params->initParams); 

  /* free TFTiling initialisation parameters */
  LALFree (params->tfTilingInput);
  
  /* EK - modified; added channel name */
  LALFree(params->channelName);

  /* free searchParams */
  LALFree (*searchParams);
  *searchParams = NULL;

  DETATCHSTATUSPTR (status);
  RETURN (status);
}

