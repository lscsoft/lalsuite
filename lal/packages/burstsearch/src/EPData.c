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

  LALFree (*vector);
  *vector = NULL;

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


