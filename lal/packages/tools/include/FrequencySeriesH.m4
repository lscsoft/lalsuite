dnl $Id$
define(`SERIESTYPE',DATATYPE`FrequencySeries')
/* <lalLaTeX file="FrequencySeriesDestroyP">
\idx{`XLALDestroy'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="FrequencySeriesDestroyP"> */
void `XLALDestroy'SERIESTYPE (
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

/* <lalLaTeX file="FrequencySeriesCutP">
\idx{`XLALCut'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="FrequencySeriesCutP"> */
SERIESTYPE *`XLALCut'SERIESTYPE (
	const SERIESTYPE *series,
	size_t first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="FrequencySeriesShrinkP">
\idx{`XLALResize'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="FrequencySeriesShrinkP"> */
SERIESTYPE *`XLALResize'SERIESTYPE (
	SERIESTYPE *series,
	int first,
	size_t length
);
/* </lalVerbatim> */

/* <lalLaTeX file="FrequencySeriesShrinkP">
\idx{`XLALShrink'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="FrequencySeriesShrinkP"> */
SERIESTYPE *`XLALShrink'SERIESTYPE (
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
	const SERIESTYPE *arg2
);
/* </lalVerbatim> */

ifelse(DATATYPE, COMPLEX8,
/* <lalLaTeX file="FrequencySeriesConjugateP">
\idx{`XLALConjugate'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="FrequencySeriesConjugateP"> */
SERIESTYPE *`XLALConjugate'SERIESTYPE (
        SERIESTYPE *series
);
/* </lalVerbatim> */
, DATATYPE, COMPLEX16,
/* <lalLaTeX file="FrequencySeriesConjugateP">
\idx{`XLALConjugate'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="FrequencySeriesConjugateP"> */
SERIESTYPE *`XLALConjugate'SERIESTYPE (
        SERIESTYPE *series
);
/* </lalVerbatim> */
,)

/* <lalLaTeX file="FrequencySeriesMultiplyP">
\idx{`XLALMultiply'SERIESTYPE ()}
</lalLaTeX> <lalVerbatim file="FrequencySeriesMultiplyP"> */
SERIESTYPE *`XLALMultiply'SERIESTYPE (
	SERIESTYPE *arg1,
	const SERIESTYPE *arg2
);
/* </lalVerbatim> */
