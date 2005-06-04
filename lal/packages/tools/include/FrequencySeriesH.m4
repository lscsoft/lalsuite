dnl $Id$
define(`SERIESTYPE',DATATYPE`FrequencySeries')
/* <lalLaTeX file="FrequencySeriesDestroyP">
\idx{`XLALDestroy'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="FrequencySeriesDestroyP"> */
void `XLALDestroy'SERIESTYPE (
	SERIESTYPE *series
);
/* </lalVerbatim> */

/* <lalLaTeX file="FrequencySeriesDestroyP">
\idx{`LALDestroy'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="FrequencySeriesDestroyP"> */
void `LALDestroy'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE *series
);
/* </lalVerbatim> */

/* <lalLaTeX file="FrequencySeriesCreateP">
\idx{`XLALCreate'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="FrequencySeriesCreateP"> */
SERIESTYPE *`XLALCreate'SERIESTYPE (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaF,
	const LALUnit *sampleUnits,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="FrequencySeriesCreateP">
\idx{`LALCreate'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="FrequencySeriesCreateP"> */
void `LALCreate'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE **output,
	const CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaF,
	LALUnit sampleUnits,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="FrequencySeriesCutP">
\idx{`XLALCut'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="FrequencySeriesCutP"> */
SERIESTYPE *`XLALCut'SERIESTYPE (
	SERIESTYPE *series,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="FrequencySeriesCutP">
\idx{`LALCut'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="FrequencySeriesCutP"> */
void `LALCut'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE **output,
	SERIESTYPE *input,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="FrequencySeriesShrinkP">
\idx{`XLALShrink'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="FrequencySeriesShrinkP"> */
size_t `XLALShrink'SERIESTYPE (
	SERIESTYPE *series,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="FrequencySeriesShrinkP">
\idx{`LALShrink'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="FrequencySeriesShrinkP"> */
void `LALShrink'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE *series,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="FrequencySeriesAddP">
\idx{`XLALAdd'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="FrequencySeriesAddP"> */
SERIESTYPE *`XLALAdd'SERIESTYPE (
	SERIESTYPE *arg1,
	SERIESTYPE *arg2
);
/* </lalVerbatim> */
