/* <lalVerbatim file="FrequencySeriesHV">
Author: Cannon, K. C.
$Id$
</lalVerbatim>
 */

#ifndef _FREQUENCYSERIES_H
#define _FREQUENCYSERIES_H

#include <lal/LALDatatypes.h>

#define LAL_FAIL_ERR	1
#define LAL_FAIL_MSG	"operation failed"
#define LAL_NULL_ERR	2
#define LAL_NULL_MSG	"unexpected NULL pointer"
#define LAL_NOMEM_ERR	3
#define LAL_NOMEM_MSG	"out of memory"
#define LAL_RANGE_ERR	4
#define LAL_RANGE_MSG	"parameter out of range"

void XLALDestroyREAL4FrequencySeries(
	REAL4FrequencySeries *series
);

void LALDestroyREAL4FrequencySeries(
	LALStatus *status,
	REAL4FrequencySeries *series
);

REAL4FrequencySeries *XLALCreateREAL4FrequencySeries(
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaF,
	LALUnit sampleUnits,
	size_t length
);

void LALCreateREAL4FrequencySeries(
	LALStatus *status,
	REAL4FrequencySeries **output,
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaF,
	LALUnit sampleUnits,
	size_t length
);

#endif  /* _FREQUENCYSERIES_H */
