#include <lal/LALStdio.h>
#include <lal/LIGOLwXMLRead.h>

NRCSID( LIGOLWXMLREADC, "$Id$" );

void LALLIGOLwXMLRead( LALStatus *status )
{
  INITSTATUS( status, "LALLIGOLwXMLRead", LIGOLWXMLREADC );
  RETURN( status );
}
