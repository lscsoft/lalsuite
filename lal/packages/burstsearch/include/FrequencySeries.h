/* <lalVerbatim file="FrequencySeriesHV">
Author: Kipp Cannon <kipp@gravity.phys.uwm.edu>
$Id$
</lalVerbatim>
 */

#ifndef _FREQUENCYSERIES_H
#define _FREQUENCYSERIES_H

#include <lal/LALDatatypes.h>

void XLALDestroyCOMPLEX8FrequencySeries(
	COMPLEX8FrequencySeries *series
);

void LALDestroyCOMPLEX8FrequencySeries(
	LALStatus *status,
	COMPLEX8FrequencySeries *series
);

void XLALDestroyREAL4FrequencySeries(
	REAL4FrequencySeries *series
);

void LALDestroyREAL4FrequencySeries(
	LALStatus *status,
	REAL4FrequencySeries *series
);

COMPLEX8FrequencySeries *XLALCreateCOMPLEX8FrequencySeries(
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaF,
	LALUnit sampleUnits,
	size_t length
);

void LALCreateCOMPLEX8FrequencySeries(
	LALStatus *status,
	COMPLEX8FrequencySeries **output,
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaF,
	LALUnit sampleUnits,
	size_t length
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
