#include <math.h>
#include <lal/LALStdlib.h>
#include "lalapps.h"
#include "gpstime.h"

RCSID( "$Id$" );

/* one second in nano seconds as an INT8 */
#define I8_SEC_NS LAL_INT8_C(1000000000)

/* one second in nano seconds as a REAL8 */
#define D_SEC_NS (1e9)

/* one nano seconds in seconds as a REAL8 */
#define D_NS_SEC (1e-9)

INT8 epoch_to_ns( LIGOTimeGPS *epoch )
{
  INT8 ns;
  ns  = (INT8)(epoch->gpsSeconds) * I8_SEC_NS;
  ns += (INT8)(epoch->gpsNanoSeconds);
  return ns;
}

LIGOTimeGPS * ns_to_epoch( LIGOTimeGPS *epoch, INT8 ns )
{
  epoch->gpsSeconds     = (INT4)( ns / I8_SEC_NS );
  epoch->gpsNanoSeconds = (INT4)( ns % I8_SEC_NS );
  return epoch;
}

INT8 sec_to_ns( REAL8 sec )
{
  INT8 ns;
  ns = (INT8) floor( D_SEC_NS * sec );
  return ns;
}

REAL8 ns_to_sec( INT8 ns )
{
  REAL8 sec;
  sec = D_NS_SEC * ns;
  return sec;
}
