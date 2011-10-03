dnl $Id$
define(`SERIESTYPE',DATATYPE`TimeSeries')
void `XLALDestroy'SERIESTYPE (
	SERIESTYPE *series
);


SERIESTYPE *`XLALCreate'SERIESTYPE (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaT,
	const LALUnit *sampleUnits,
	size_t length
);


SERIESTYPE *`XLALCut'SERIESTYPE (
	const SERIESTYPE *series,
	size_t first,
	size_t length
);


SERIESTYPE *`XLALResize'SERIESTYPE (
	SERIESTYPE *series,
	int first,
	size_t length
);


SERIESTYPE *`XLALShrink'SERIESTYPE (
	SERIESTYPE *series,
	size_t first,
	size_t length
);


SERIESTYPE *`XLALAdd'SERIESTYPE (
	SERIESTYPE *arg1,
	const SERIESTYPE *arg2
);
