/* <lalVerbatim file="TimeSeriesHV">
Author: Cannon, K. C.
$Id$
</lalVerbatim>
 */

#ifndef _TIMESERIES_H
#define _TIMESERIES_H

#include <lal/LALDatatypes.h>

#define LAL_FAIL_ERR	1
#define LAL_FAIL_MSG	"operation failed"
#define LAL_NULL_ERR	2
#define LAL_NULL_MSG	"unexpected NULL pointer"
#define LAL_NOMEM_ERR	3
#define LAL_NOMEM_MSG	"out of memory"
#define LAL_RANGE_ERR	4
#define LAL_RANGE_MSG	"parameter out of range"

void XLALDestroyREAL4TimeSeries(
	REAL4TimeSeries *series
);

void LALDestroyREAL4TimeSeries(
	LALStatus *status,
	REAL4TimeSeries *series
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
