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

  INITSTATUS (status, "LALCreateEPDataSegmentVector", EPDATAC);
  ATTATCHSTATUSPTR (status);


  ASSERT (!*vector, status, EPDATA_ENNUL, EPDATA_MSGENNUL);

  ASSERT (params, status, EPDATA_ENULL, EPDATA_MSGENULL);
  ASSERT (params->numSegments > 0, status, EPDATA_ESEGZ, EPDATA_MSGESEGZ);
  ASSERT (params->numPoints > 0, status, EPDATA_ENUMZ, EPDATA_MSGENUMZ);

  vectorPtr = *vector = (EPDataSegmentVector *) LALMalloc (sizeof(EPDataSegmentVector));
  vectorPtr->length = params->numSegments;

  vectorPtr->series = LALMalloc(params->numSegments * sizeof(*vectorPtr->series));
  if ( !vectorPtr->series )
  {
    LALFree( vectorPtr );
    vectorPtr = NULL;
    ABORT( status, EPDATA_EALOC, EPDATA_MSGEALOC );
  }

  for (i = 0; i < (INT4)vectorPtr->length ; ++i)
  {
    strncpy( vectorPtr->series[i].name, "anonymous", LALNameLength * sizeof(CHAR) );
    vectorPtr->series[i].data = NULL;

    LALCreateVector (status->statusPtr, &vectorPtr->series[i].data, 2 * params->numPoints);
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

  INITSTATUS (status, "LALDestroyEPDataSegmentVector", EPDATAC);
  ATTATCHSTATUSPTR (status);


  ASSERT (*vector, status, 
      EPDATA_ENULL, EPDATA_MSGENULL);

 for (i = 0; i < (INT4)(*vector)->length; ++i)
  {
    LALDestroyVector (status->statusPtr, &(*vector)->series[i].data);
    CHECKSTATUSPTR (status);
  }

  LALFree ((*vector)->series);

  LALFree (*vector);
  *vector = NULL;

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


