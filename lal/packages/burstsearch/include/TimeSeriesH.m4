dnl $Id$
define(`SERIESTYPE',DATATYPE`TimeSeries')
/* <lalVerbatim file="TimeSeriesDestroyP"> */
void `XLALDestroy'SERIESTYPE (
	SERIESTYPE *series
);
/* </lalVerbatim> */

/* <lalVerbatim file="TimeSeriesDestroyP"> */
void `LALDestroy'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE *series
);
/* </lalVerbatim> */

/* <lalVerbatim file="TimeSeriesCreateP"> */
SERIESTYPE *`XLALCreate'SERIESTYPE (
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaT,
	LALUnit sampleUnits,
	size_t length
);
/* </lalVerbatim> */

/* <lalVerbatim file="TimeSeriesCreateP"> */
void `LALCreate'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE **output,
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaT,
	LALUnit sampleUnits,
	size_t length
);
/* </lalVerbatim> */

/* <lalVerbatim file="TimeSeriesCutP"> */
SERIESTYPE *`XLALCut'SERIESTYPE (
	SERIESTYPE *series,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalVerbatim file="TimeSeriesCutP"> */
void `LALCut'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE **output,
	SERIESTYPE *input,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalVerbatim file="TimeSeriesShrinkP"> */
size_t `XLALShrink'SERIESTYPE (
	SERIESTYPE *series,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalVerbatim file="TimeSeriesShrinkP"> */
void `LALShrink'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE *series,
	size_t first,
	size_t length
);
/* </lalVerbatim> */
