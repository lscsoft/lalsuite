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
#include "EPData.h"

NRCSID (EPDATAC, "$Id$");

void
LALCreateEPDataSegmentVector (
    LALStatus                  *status,
    EPDataSegmentVector       **vector,
    EPInitParams               *params
    )
{
  INITSTATUS (status, "LALCreateEPDataSegmentVector", EPDATAC);
  ATTATCHSTATUSPTR (status);

  ASSERT (vector, status, EPDATA_ENUL, EPDATA_MSGENUL);

  *vector = NULL;

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
LALDestroyEPDataSegmentVector (
    LALStatus                  *status,
    EPDataSegmentVector       **vector
    )
{
  INITSTATUS (status, "LALDestroyEPDataSegmentVector", EPDATAC);
  ATTATCHSTATUSPTR (status);

  ASSERT (vector, status, EPDATA_ENULL, EPDATA_MSGENULL);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}
