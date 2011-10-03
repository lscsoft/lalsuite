dnl $Id$
define(`SERIESTYPE',DATATYPE`FrequencySeries')
void `XLALDestroy'SERIESTYPE (
	SERIESTYPE *series
);


SERIESTYPE *`XLALCreate'SERIESTYPE (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaF,
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


ifelse(DATATYPE, COMPLEX8,
SERIESTYPE *`XLALConjugate'SERIESTYPE (
        SERIESTYPE *series
);

, DATATYPE, COMPLEX16,
SERIESTYPE *`XLALConjugate'SERIESTYPE (
        SERIESTYPE *series
);

,)

SERIESTYPE *`XLALMultiply'SERIESTYPE (
	SERIESTYPE *arg1,
	const SERIESTYPE *arg2
);
