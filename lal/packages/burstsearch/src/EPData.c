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

  ASSERT (vector, status, EPDATA_ENUL, EPDATA_MSGENUL);

  *vector = NULL;

  RETURN (status);
}


void
LALDestroyEPDataSegmentVector (
    LALStatus                  *status,
    EPDataSegmentVector       **vector
    )
{
  INITSTATUS (status, "LALDestroyEPDataSegmentVector", EPDATAC);

  ASSERT (vector, status, EPDATA_ENULL, EPDATA_MSGENULL);

  *vector = NULL;

  RETURN (status);
}
