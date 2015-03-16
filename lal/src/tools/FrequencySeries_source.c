#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)

#define SERIESTYPE CONCAT2(DATATYPE,FrequencySeries)
#define SEQUENCETYPE CONCAT2(DATATYPE,Sequence)

#define DSERIES CONCAT2(XLALDestroy,SERIESTYPE)
#define CSERIES CONCAT2(XLALCreate,SERIESTYPE)
#define XSERIES CONCAT2(XLALCut,SERIESTYPE)
#define RSERIES CONCAT2(XLALResize,SERIESTYPE)
#define SSERIES CONCAT2(XLALShrink,SERIESTYPE)
#define ASERIES CONCAT2(XLALAdd,SERIESTYPE)
#define MSERIES CONCAT2(XLALMultiply,SERIESTYPE)

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
	REAL8 deltaF,
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
	new->deltaF = deltaF;
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
	new->f0 += first * new->deltaF;

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
	series->f0 += first * series->deltaF;

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
	REAL8 Delta_f0 = arg2->f0 - arg1->f0;
	/* number of arg1 units per arg2 unit.  XLALUnitRatio() returns the
	 * number one obtains when one divides 1 of the first argument by 1
	 * of the second argument, for example if arg2 is in m and arg1 is
	 * in cm then unit_ratio = 100.0 */
	REAL8 unit_ratio = XLALUnitRatio(&arg2->sampleUnits, &arg1->sampleUnits);
	unsigned i, j;

	/* make sure arguments are compatible */
	if(XLALIsREAL8FailNaN(unit_ratio)) {
		XLALPrintError("%s(): incompatible sample units\n", __func__);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}
	if(XLALGPSCmp(&arg1->epoch, &arg2->epoch)) {
		XLALPrintError("%s(): incompatible epochs\n", __func__);
		XLAL_ERROR_NULL(XLAL_EDATA);
	}
	if(arg1->deltaF != arg2->deltaF) {
		XLALPrintError("%s(): incompatible frequency resolution\n", __func__);
		XLAL_ERROR_NULL(XLAL_EDATA);
	}

	/* set start indexes */
	if(Delta_f0 >= 0) {
		i = floor(Delta_f0 / arg1->deltaF + 0.5);
		j = 0;
	} else {
		i = 0;
		j = floor(-Delta_f0 / arg2->deltaF + 0.5);
	}

	/* add arg2 to arg1, adjusting the units */
	for(; i < arg1->data->length && j < arg2->data->length; i++, j++)
	{
		arg1->data->data[i] += arg2->data->data[j] * unit_ratio;
	}

	return arg1;
}


SERIESTYPE *MSERIES (
	SERIESTYPE *arg1,
	const SERIESTYPE *arg2
)
{
	REAL8 Delta_f0 = arg2->f0 - arg1->f0;
	/* number of arg1 units per arg2 unit.  XLALUnitRatio() returns the
	 * number one obtains when one divides 1 of the first argument by 1
	 * of the second argument, for example if arg2 is in m and arg1 is
	 * in cm then unit_ratio = 100.0 */
	REAL8 unit_ratio = XLALUnitRatio(&arg2->sampleUnits, &arg1->sampleUnits);
	unsigned i, j;

	/* make sure arguments are compatible */

	if(XLALIsREAL8FailNaN(unit_ratio)) {
		XLALPrintError("%s(): incompatible sample units\n", __func__);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}
	if(XLALGPSCmp(&arg1->epoch, &arg2->epoch)) {
		XLALPrintError("%s(): incompatible epochs\n", __func__);
		XLAL_ERROR_NULL(XLAL_EDATA);
	}
	if(arg1->deltaF != arg2->deltaF) {
		XLALPrintError("%s(): incompatible frequency resolution\n", __func__);
		XLAL_ERROR_NULL(XLAL_EDATA);
	}

	/* set start indexes */

	if(Delta_f0 >= 0) {
		unsigned start = floor(Delta_f0 / arg1->deltaF + 0.5);

		/* zero the first part of arg1 */

		for(i = 0; i < start && i < arg1->data->length; i++)
		{
			arg1->data->data[i] = 0.0;
		}
		/* now i = min(start, arg1->data->length) */

		j = 0;
	} else {
		i = 0;
		j = floor(-Delta_f0 / arg2->deltaF + 0.5);
	}

	/* multiply arg1 by arg2 */

	for(; i < arg1->data->length && j < arg2->data->length; i++, j++)
	{
		arg1->data->data[i] *= arg2->data->data[j] * unit_ratio;
	}

	/* zero the last part of arg1 */

	for(; i < arg1->data->length; i++)
	{
		arg1->data->data[i] = 0.0;
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
