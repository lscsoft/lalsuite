dnl $Id$
define(`SERIESTYPE',DATATYPE`FrequencySeries')
/* <lalVerbatim file="FrequencySeriesDestroyP"> */
void `XLALDestroy'SERIESTYPE (
	SERIESTYPE *series
);
/* </lalVerbatim> */

/* <lalVerbatim file="FrequencySeriesDestroyP"> */
void `LALDestroy'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE *series
);
/* </lalVerbatim> */

/* <lalVerbatim file="FrequencySeriesCreateP"> */
SERIESTYPE *`XLALCreate'SERIESTYPE (
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaF,
	LALUnit sampleUnits,
	size_t length
);
/* </lalVerbatim> */

/* <lalVerbatim file="FrequencySeriesCreateP"> */
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
/* </lalVerbatim> */

/* <lalVerbatim file="FrequencySeriesCutP"> */
SERIESTYPE *`XLALCut'SERIESTYPE (
	SERIESTYPE *series,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalVerbatim file="FrequencySeriesCutP"> */
void `LALCut'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE **output,
	SERIESTYPE *input,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalVerbatim file="FrequencySeriesShrinkP"> */
size_t `XLALShrink'SERIESTYPE (
	SERIESTYPE *series,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalVerbatim file="FrequencySeriesShrinkP"> */
void `LALShrink'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE *series,
	size_t first,
	size_t length
);
/* </lalVerbatim> */
