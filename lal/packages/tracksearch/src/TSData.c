/*----------------------------------------------------------------------- 
 * 
 * File Name: TSData.c
 *
 * Author: Torres C
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include "TSData.h"

NRCSID (TSDATAC, "$Id$");

void
LALCreateTSDataSegmentVector (
    LALStatus                    *status,
    TSSegmentVector             **vector,
    TSCreateParams               *params
    )
{
  INT4                           i;
  TSSegmentVector       *vectorPtr;
  /* REAL4TimeSeries          *segPtr; */

  INITSTATUS (status, "LALCreateTSSegmentVector", TSDATAC);
  ATTATCHSTATUSPTR (status);

  ASSERT (!*vector, status, TSDATA_ENNUL, TSDATA_MSGENNUL);
  ASSERT (params, status, TSDATA_ENULL, TSDATA_MSGENULL);
  ASSERT (params->numberDataSegments > 0, status, TSDATA_ESEGZ, TSDATA_MSGESEGZ);
  ASSERT (params->dataSegmentPoints > 0, status, TSDATA_ENUMZ, TSDATA_MSGENUMZ);

  vectorPtr = *vector = (TSSegmentVector *) LALMalloc (sizeof(TSSegmentVector));
  vectorPtr->length = params->numberDataSegments;

  vectorPtr->dataSeg  = (REAL4TimeSeries *) 
    LALMalloc (params->numberDataSegments*sizeof(REAL4TimeSeries));
  if ( !(vectorPtr->dataSeg) )
  {
    LALFree( vectorPtr );
    vectorPtr = NULL;
    ABORT( status, TSDATA_EALOC, TSDATA_MSGEALOC );
  }
  for (i = 0; i < (INT4)(vectorPtr->length) ; ++i)
  {
    strncpy( vectorPtr->dataSeg->name, "anonymous", LALNameLength * sizeof(CHAR) );
    vectorPtr->dataSeg[i].data = NULL;
    LALCreateVector (status->statusPtr, 
        &(vectorPtr->dataSeg[i].data), params->dataSegmentPoints);
    CHECKSTATUSPTR (status);
  }
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



void
LALDestroyTSDataSegmentVector (
    LALStatus                  *status,
    TSSegmentVector           **vector
    )
{
  INT4                    i;
  TSSegmentVector       *vectorPtr;

  INITSTATUS (status, "LALDestroyTSDataSegmentVector", TSDATAC);
  ATTATCHSTATUSPTR (status);


  ASSERT (*vector, status, 
      TSDATA_ENULL, TSDATA_MSGENULL);

  vectorPtr = (*vector);
  for (i = 0; i < (INT4)(vectorPtr->length) ; ++i)
  {
    LALDestroyVector (status->statusPtr, &(vectorPtr->dataSeg[i].data));
    CHECKSTATUSPTR (status);
  }
  LALFree (vectorPtr->dataSeg);
  LALFree (*vector);
  /*  *vector = NULL; */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


