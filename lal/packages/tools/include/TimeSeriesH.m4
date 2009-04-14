dnl $Id$
define(`SERIESTYPE',DATATYPE`TimeSeries')
/* <lalLaTeX file="TimeSeriesDestroyP">
\idx{`XLALDestroy'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesDestroyP"> */
void `XLALDestroy'SERIESTYPE (
	SERIESTYPE *series
);
/* </lalVerbatim> */

/* <lalLaTeX file="TimeSeriesCreateP">
\idx{`XLALCreate'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesCreateP"> */
SERIESTYPE *`XLALCreate'SERIESTYPE (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaT,
	const LALUnit *sampleUnits,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="TimeSeriesCutP">
\idx{`XLALCut'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesCutP"> */
SERIESTYPE *`XLALCut'SERIESTYPE (
	const SERIESTYPE *series,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="TimeSeriesResizeP">
\idx{`XLALResize'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesResizeP"> */
SERIESTYPE *`XLALResize'SERIESTYPE (
	SERIESTYPE *series,
	int first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="TimeSeriesResizeP">
\idx{`XLALShrink'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesResizeP"> */
SERIESTYPE *`XLALShrink'SERIESTYPE (
	SERIESTYPE *series,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="TimeSeriesAddP">
\idx{`XLALAdd'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesAddP"> */
SERIESTYPE *`XLALAdd'SERIESTYPE (
	SERIESTYPE *arg1,
	const SERIESTYPE *arg2
);
/* </lalVerbatim> */

