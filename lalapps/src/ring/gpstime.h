#ifndef GPSTIME_H
#define GPSTIME_H

/*
 *
 * Routines to convert times between LIGOTimeGPS epochs, INT8 nanoseconds,
 * and REAL8 seconds.
 *
 */

#include <lal/LALDatatypes.h>

/* convert from LIGOTimeGPS epoch to INT8 nanoseconds */
INT8 epoch_to_ns( LIGOTimeGPS *epoch );

/* convert from INT8 nanoseconds to LIGOTimeGPS epoch */
LIGOTimeGPS * ns_to_epoch( LIGOTimeGPS *epoch, INT8 ns );

/* convert from REAL8 seconds to INT8 nanoseconds */
INT8 sec_to_ns( REAL8 sec );

/* convert from INT8 nanoseconds to REAL8 seconds */
REAL8 ns_to_sec( INT8 ns );

#endif /* GPSTIME_H */
