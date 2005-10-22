dnl "$Id"
define(`SERIESTYPE',DATATYPE`TimeSeries')
define(`SEQUENCETYPE',DATATYPE`Sequence')
void `XLALDestroy'SERIESTYPE (
	SERIESTYPE *series
)
{
	if(series)
		`XLALDestroy'SEQUENCETYPE (series->data);
	LALFree(series);
}


void `LALDestroy'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE *series
)
{
	INITSTATUS(status, "`LALDestroy'SERIESTYPE", TIMESERIESC);
	`XLALDestroy'SERIESTYPE (series);
	RETURN(status);
}


SERIESTYPE *`XLALCreate'SERIESTYPE (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaT,
	const LALUnit *sampleUnits,
	size_t length
)
{
	static const char *func = "`XLALCreate'SERIESTYPE";
	SERIESTYPE *new;
	SEQUENCETYPE *sequence;

	new = LALMalloc(sizeof(*new));
	sequence = `XLALCreate'SEQUENCETYPE (length);
	if(!new || !sequence) {
		LALFree(new);
		`XLALDestroy'SEQUENCETYPE (sequence);
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
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

	return(new);
}


void `LALCreate'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE **output,
	const CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaT,
	LALUnit sampleUnits,
	size_t length
)
{
	INITSTATUS(status, "`LALCreate'SERIESTYPE", TIMESERIESC);
	ASSERT(output != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	*output = `XLALCreate'SERIESTYPE (name, &epoch, f0, deltaT, &sampleUnits, length);
	if(*output == NULL) {
		XLALClearErrno();
		ABORT(status, LAL_NOMEM_ERR, LAL_NOMEM_MSG);
	}
	RETURN(status);
}


SERIESTYPE *`XLALCut'SERIESTYPE (
	const SERIESTYPE *series,
	size_t first,
	size_t length
)
{
	static const char *func = "`XLALCut'SERIESTYPE";
	SERIESTYPE *new;
	SEQUENCETYPE *sequence;

	if(!series || !series->data)
		return(NULL);

	new = LALMalloc(sizeof(*new));
	sequence = `XLALCut'SEQUENCETYPE (series->data, first, length);
	if(!new || !sequence) {
		LALFree(new);
		`XLALDestroy'SEQUENCETYPE (sequence);
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
	}

	*new = *series;
	new->data = sequence;
	XLALGPSAdd(&new->epoch, first * new->deltaT);

	return(new);
}


void `LALCut'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE **output,
	const SERIESTYPE *input,
	size_t first,
	size_t length
)
{
	INITSTATUS(status, "`LALCut'SERIESTYPE", TIMESERIESC);

	ASSERT(output != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(input != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(input->data != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(first + length <= input->data->length, status, LAL_RANGE_ERR, LAL_RANGE_MSG);
	*output = `XLALCut'SERIESTYPE (input, first, length);
	if(*output == NULL) {
		XLALClearErrno();
		ABORT(status, LAL_NOMEM_ERR, LAL_NOMEM_MSG);
	}

	RETURN(status);
}


size_t `XLALResize'SERIESTYPE (
	SERIESTYPE *series,
	int first,
	size_t length
)
{
	if(!series)
		return(0);

	XLALGPSAdd(&series->epoch, first * series->deltaT);
	return(`XLALResize'SEQUENCETYPE (series->data, first, length));
}


void `LALResize'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE *series,
	int first,
	size_t length
)
{
	size_t new_length;

	INITSTATUS(status, "`LALResize'SERIESTYPE", TIMESERIESC);

	ASSERT(series != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(series->data != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	new_length = `XLALResize'SERIESTYPE (series, first, length);
	if(new_length != length) {
		XLALClearErrno();
		ABORT(status, LAL_FAIL_ERR, LAL_FAIL_MSG);
	}
	RETURN(status);
}


size_t `XLALShrink'SERIESTYPE (
	SERIESTYPE *series,
	size_t first,
	size_t length
)
{
	return(`XLALResize'SERIESTYPE (series, first, length));
}


void `LALShrink'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE *series,
	size_t first,
	size_t length
)
{
	size_t new_length;

	INITSTATUS(status, "`LALShrink'SERIESTYPE", TIMESERIESC);

	ASSERT(series != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(series->data != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(first + length <= series->data->length, status, LAL_RANGE_ERR, LAL_RANGE_MSG);
	new_length = `XLALShrink'SERIESTYPE (series, first, length);
	if(new_length != length) {
		XLALClearErrno();
		ABORT(status, LAL_FAIL_ERR, LAL_FAIL_MSG);
	}
	RETURN(status);
}

