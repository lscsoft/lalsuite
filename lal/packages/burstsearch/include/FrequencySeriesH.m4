dnl $Id$
define(SERIESTYPE,DATATYPE`FrequencySeries')
void `XLALDestroy'SERIESTYPE (
	SERIESTYPE *series
);

void `LALDestroy'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE *series
);

SERIESTYPE *`XLALCreate'SERIESTYPE (
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaF,
	LALUnit sampleUnits,
	size_t length
);

void `LALCreate'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE **output,
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaF,
	LALUnit sampleUnits,
	size_t length
);

SERIESTYPE *`XLALCut'SERIESTYPE (
	SERIESTYPE *series,
	size_t first,
	size_t length
);

void `LALCut'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE **output,
	SERIESTYPE *input,
	size_t first,
	size_t length
);
