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
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaT,
	const LALUnit *sampleUnits,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="TimeSeriesCreateP">
\idx{`LALCreate'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesCreateP"> */
void `LALCreate'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE **output,
	const CHAR *name,
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
	const SERIESTYPE *series,
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
	const SERIESTYPE *input,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="TimeSeriesResizeP">
\idx{`XLALResize'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesResizeP"> */
size_t `XLALResize'SERIESTYPE (
	SERIESTYPE *series,
	int first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="TimeSeriesResizeP">
\idx{`LALResize'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesResizeP"> */
void `LALResize'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE *series,
	int first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="TimeSeriesResizeP">
\idx{`XLALShrink'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesResizeP"> */
size_t `XLALShrink'SERIESTYPE (
	SERIESTYPE *series,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="TimeSeriesResizeP">
\idx{`LALShrink'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="TimeSeriesResizeP"> */
void `LALShrink'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE *series,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

