dnl $Id$
define(`SERIESTYPE',DATATYPE`TimeSeries')
/* <lalLaTeX file="TimeSeriesDestroyP">
\idx{`XLALDestroy'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesDestroyP"> */
void `XLALDestroy'SERIESTYPE (
	SERIESTYPE *series
);
/* </lalVerbatim> */

/* <lalLaTeX file="TimeSeriesDestroyP">
\idx{`LALDestroy'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesDestroyP"> */
void `LALDestroy'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE *series
);
/* </lalVerbatim> */

/* <lalLaTeX file="TimeSeriesCreateP">
\idx{`XLALCreate'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesCreateP"> */
SERIESTYPE *`XLALCreate'SERIESTYPE (
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaT,
	LALUnit sampleUnits,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="TimeSeriesCreateP">
\idx{`LALCreate'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesCreateP"> */
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

/* <lalLaTeX file="TimeSeriesCutP">
\idx{`XLALCut'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesCutP"> */
SERIESTYPE *`XLALCut'SERIESTYPE (
	SERIESTYPE *series,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="TimeSeriesCutP">
\idx{`LALCut'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesCutP"> */
void `LALCut'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE **output,
	SERIESTYPE *input,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="TimeSeriesShrinkP">
\idx{`XLALShrink'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesShrinkP"> */
size_t `XLALShrink'SERIESTYPE (
	SERIESTYPE *series,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="TimeSeriesShrinkP">
\idx{`LALShrink'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesShrinkP"> */
void `LALShrink'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE *series,
	size_t first,
	size_t length
);
/* </lalVerbatim> */
