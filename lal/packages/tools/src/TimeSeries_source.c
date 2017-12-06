#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define SERIESTYPE CONCAT2(DATATYPE,TimeSeries)
#define SEQUENCETYPE CONCAT2(DATATYPE,Sequence)

#define DSERIES CONCAT2(XLALDestroy,SERIESTYPE)
#define CSERIES CONCAT2(XLALCreate,SERIESTYPE)
#define XSERIES CONCAT2(XLALCut,SERIESTYPE)
#define RSERIES CONCAT2(XLALResize,SERIESTYPE)
#define SSERIES CONCAT2(XLALShrink,SERIESTYPE)
#define ASERIES CONCAT2(XLALAdd,SERIESTYPE)

#define DSEQUENCE CONCAT2(XLALDestroy,SEQUENCETYPE)
#define CSEQUENCE CONCAT2(XLALCreate,SEQUENCETYPE)
#define XSEQUENCE CONCAT2(XLALCut,SEQUENCETYPE)
#define RSEQUENCE CONCAT2(XLALResize,SEQUENCETYPE)

void DSERIES (
	SERIESTYPE *series
)
{
	if(series)
		DSEQUENCE (series->data);
	XLALFree(series);
}


SERIESTYPE *CSERIES (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaT,
	const LALUnit *sampleUnits,
	size_t length
)
{
	SERIESTYPE *new;
	SEQUENCETYPE *sequence;

	new = XLALMalloc(sizeof(*new));
	sequence = CSEQUENCE (length);
	if(!new || !sequence) {
		XLALFree(new);
		DSEQUENCE (sequence);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	if(name)
		strncpy(new->name, name, LALNameLength);
	else
		new->name[0] = '\0';
	new->epoch = *epoch;
	new->f0 = f0;
	new->deltaT = deltaT;
	new->sampleUnits = *sampleUnits;
	new->data = sequence;

	return new;
}


SERIESTYPE *XSERIES (
	const SERIESTYPE *series,
	size_t first,
	size_t length
)
{
	SERIESTYPE *new;
	SEQUENCETYPE *sequence;

	new = XLALMalloc(sizeof(*new));
	sequence = XSEQUENCE (series->data, first, length);
	if(!new || !sequence) {
		XLALFree(new);
		DSEQUENCE (sequence);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	*new = *series;
	new->data = sequence;
	XLALGPSAdd(&new->epoch, first * new->deltaT);

	return new;
}


SERIESTYPE *RSERIES (
	SERIESTYPE *series,
	int first,
	size_t length
)
{

	if(!RSEQUENCE (series->data, first, length))
		XLAL_ERROR_NULL(XLAL_EFUNC);
	XLALGPSAdd(&series->epoch, first * series->deltaT);

	return series;
}


SERIESTYPE *SSERIES (
	SERIESTYPE *series,
	size_t first,
	size_t length
)
{

	if(!RSERIES (series, first, length))
		XLAL_ERROR_NULL(XLAL_EFUNC);

	return series;
}


SERIESTYPE *ASERIES (
	SERIESTYPE *arg1,
	const SERIESTYPE *arg2
)
{
	REAL8 Delta_epoch = XLALGPSDiff(&arg2->epoch, &arg1->epoch);
	/* number of arg1 units per arg2 unit.  XLALUnitRatio() returns the
	 * number one obtains when one divides 1 of the first argument by 1
	 * of the second argument, for example if arg2 is in m and arg1 is
	 * in cm then unit_ratio = 100.0 */
	REAL8 unit_ratio = XLALUnitRatio(&arg2->sampleUnits, &arg1->sampleUnits);
	unsigned i, j;

	/* make sure arguments are compatible */
	if(XLALIsREAL8FailNaN(unit_ratio)) {
		XLALPrintError("%s(): incompatible units\n", __func__);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}
	if(arg1->f0 != arg2->f0) {
		XLALPrintError("%s(): incompatible heterodyne frequencies\n", __func__);
		XLAL_ERROR_NULL(XLAL_EDATA);
	}
	if(arg1->deltaT != arg2->deltaT) {
		XLALPrintError("%s(): incompatible sample periods\n", __func__);
		XLAL_ERROR_NULL(XLAL_EDATA);
	}

	/* set start indexes */
	if(Delta_epoch >= 0) {
		i = floor(Delta_epoch / arg1->deltaT + 0.5);
		j = 0;
	} else {
		i = 0;
		j = floor(-Delta_epoch / arg2->deltaT + 0.5);
	}

	/* add arg2 to arg1, adjusting the units */
	for(; i < arg1->data->length && j < arg2->data->length; i++, j++)
	{
		arg1->data->data[i] += arg2->data->data[j] * unit_ratio;
	}

	return arg1;
}

#undef SERIESTYPE
#undef SEQUENCETYPE

#undef DSERIES
#undef CSERIES
#undef XSERIES
#undef RSERIES
#undef SSERIES
#undef ASERIES
#undef MSERIES

#undef DSEQUENCE
#undef CSEQUENCE
#undef XSEQUENCE
#undef RSEQUENCE
