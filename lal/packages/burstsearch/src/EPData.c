/*----------------------------------------------------------------------- 
 * 
 * File Name: EPData.c
 *
 * Author: Brady, P R
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include "EPData.h"

NRCSID (EPDATAC, "$Id$");

void
LALCreateEPDataSegmentVector (
    LALStatus                  *status,
    EPDataSegmentVector       **vector,
    EPInitParams               *params
    )
{
  INT4                  i;
  EPDataSegmentVector    *vectorPtr;
  EPDataSegment          *segPtr;

  INITSTATUS (status, "LALCreateEPDataSegmentVector", EPDATAC);
  ATTATCHSTATUSPTR (status);


  ASSERT (!*vector, status, EPDATA_ENNUL, EPDATA_MSGENNUL);

  ASSERT (params, status, EPDATA_ENULL, EPDATA_MSGENULL);
  ASSERT (params->numSegments > 0, status, EPDATA_ESEGZ, EPDATA_MSGESEGZ);
  ASSERT (params->numPoints > 0, status, EPDATA_ENUMZ, EPDATA_MSGENUMZ);

  vectorPtr = *vector = (EPDataSegmentVector *) LALMalloc (sizeof(EPDataSegmentVector));
  vectorPtr->length = params->numSegments;

  segPtr = vectorPtr->data = (EPDataSegment *) 
    LALMalloc (params->numSegments*sizeof(EPDataSegment));
  if ( !segPtr )
  {
    LALFree( vectorPtr );
    vectorPtr = NULL;
    ABORT( status, EPDATA_EALOC, EPDATA_MSGEALOC );
  }

  for (i = 0; i < vectorPtr->length ; ++i)
  {
    segPtr[i].data = (REAL4TimeSeries *) 
      LALMalloc (sizeof(REAL4TimeSeries));
    segPtr[i].spec = (REAL4FrequencySeries *) 
      LALMalloc (sizeof(REAL4FrequencySeries));
    segPtr[i].resp = (COMPLEX8FrequencySeries *) 
      LALMalloc (sizeof(COMPLEX8FrequencySeries));

    strncpy( segPtr[i].data->name, "anonymous", LALNameLength * sizeof(CHAR) );
    segPtr[i].data->data        = NULL;
    strncpy( segPtr[i].spec->name, "anonymous", LALNameLength * sizeof(CHAR) );
    segPtr[i].spec->data        = NULL;
    strncpy( segPtr[i].resp->name, "anonymous", LALNameLength * sizeof(CHAR) );
    segPtr[i].resp->data        = NULL;

    LALCreateVector (status->statusPtr, 
        &segPtr[i].data->data, params->numPoints);
    CHECKSTATUSPTR (status);

    LALCreateVector (status->statusPtr, 
        &segPtr[i].spec->data, params->numPoints/2 + 1);
    CHECKSTATUSPTR (status);

    LALCCreateVector (status->statusPtr, 
        &segPtr[i].resp->data, params->numPoints/2 + 1);
    CHECKSTATUSPTR (status);
    
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}



void
LALDestroyEPDataSegmentVector (
    LALStatus                  *status,
    EPDataSegmentVector       **vector
    )
{
  INT4                  i;
  EPDataSegment        *segPtr;

  INITSTATUS (status, "LALDestroyEPDataSegmentVector", EPDATAC);
  ATTATCHSTATUSPTR (status);


  ASSERT (*vector, status, 
      EPDATA_ENULL, EPDATA_MSGENULL);

  segPtr = (*vector)->data;
  
 for (i = 0; i < (*vector)->length; ++i)
  {
    LALCDestroyVector (status->statusPtr, &segPtr[i].resp->data);
    CHECKSTATUSPTR (status);

    LALDestroyVector (status->statusPtr, &segPtr[i].spec->data);
    CHECKSTATUSPTR (status);

    LALDestroyVector (status->statusPtr, &segPtr[i].data->data);
    CHECKSTATUSPTR (status);

    LALFree (segPtr[i].resp);
    LALFree (segPtr[i].spec);
    LALFree (segPtr[i].data);
  }

  LALFree (segPtr);

  LALFree (*vector);
  *vector = NULL;

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


