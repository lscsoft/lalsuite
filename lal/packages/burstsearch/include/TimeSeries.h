/* <lalVerbatim file="TimeSeriesHV">
Author: Cannon, K. C.
$Id$
</lalVerbatim>
 */

#ifndef _TIMESERIES_H
#define _TIMESERIES_H

#include <lal/LALDatatypes.h>

void XLALDestroyCOMPLEX8TimeSeries(
	COMPLEX8TimeSeries *series
);

void LALDestroyCOMPLEX8TimeSeries(
	LALStatus *status,
	COMPLEX8TimeSeries *series
);

void XLALDestroyREAL4TimeSeries(
	REAL4TimeSeries *series
);

void LALDestroyREAL4TimeSeries(
	LALStatus *status,
	REAL4TimeSeries *series
);

COMPLEX8TimeSeries *XLALCreateCOMPLEX8TimeSeries(
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaT,
	LALUnit sampleUnits,
	size_t length
);

void LALCreateCOMPLEX8TimeSeries(
	LALStatus *status,
	COMPLEX8TimeSeries **output,
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaT,
	LALUnit sampleUnits,
	size_t length
);

REAL4TimeSeries *XLALCreateREAL4TimeSeries(
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaT,
	LALUnit sampleUnits,
	size_t length
);

void LALCreateREAL4TimeSeries(
	LALStatus *status,
	REAL4TimeSeries **output,
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaT,
	LALUnit sampleUnits,
	size_t length
);

REAL4TimeSeries *XLALCutREAL4TimeSeries(
	REAL4TimeSeries *series,
	size_t first_sample,
	size_t num_samples
);

void LALCutREAL4TimeSeries(
	LALStatus *status,
	REAL4TimeSeries **output,
	REAL4TimeSeries *input,
	size_t first_sample,
	size_t num_samples
);

#endif  /* _TIMESERIES_H */
