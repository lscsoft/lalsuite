#ifndef GPSTIME_H
#define GPSTIME_H

#include <lal/LALDatatypes.h>

INT8 epoch_to_ns( LIGOTimeGPS *epoch );
LIGOTimeGPS * ns_to_epoch( LIGOTimeGPS *epoch, INT8 ns );
INT8 sec_to_ns( REAL8 sec );
REAL8 ns_to_sec( INT8 ns );

#endif /* GPSTIME_H */
