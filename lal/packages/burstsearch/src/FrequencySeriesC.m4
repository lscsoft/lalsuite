dnl $Id$
define(`SERIESTYPE',DATATYPE`FrequencySeries')
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
	INITSTATUS(status, "`LALDestroy'SERIESTYPE", FREQUENCYSERIESC);
	`XLALDestroy'SERIESTYPE (series);
	RETURN(status);
}


SERIESTYPE *`XLALCreate'SERIESTYPE (
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaF,
	LALUnit sampleUnits,
	size_t length
)
{
	SERIESTYPE *new;
	SEQUENCETYPE *sequence;

	new = LALMalloc(sizeof(*new));
	sequence = `XLALCreate'SEQUENCETYPE (length);
	if(!new || !sequence) {
		LALFree(new);
		`XLALDestroy'SEQUENCETYPE (sequence);
		return(NULL);
	}

	if(new->name)
		strncpy(new->name, name, LALNameLength);
	else
		new->name[0] = '\0';
	new->epoch = epoch;
	new->f0 = f0;
	new->deltaF = deltaF;
	new->sampleUnits = sampleUnits;
	new->data = sequence;

	return(new);
}


void `LALCreate'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE **output,
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaF,
	LALUnit sampleUnits,
	size_t length
)
{
	INITSTATUS(status, "`LALCreate'SERIESTYPE", FREQUENCYSERIESC);
	ASSERT(output != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	*output = `XLALCreate'SERIESTYPE (name, epoch, f0, deltaF, sampleUnits, length);
	ASSERT(*output != NULL, status, LAL_NOMEM_ERR, LAL_NOMEM_MSG);
	RETURN(status);
}


SERIESTYPE *`XLALCut'SERIESTYPE (
	SERIESTYPE *series,
	size_t first,
	size_t length
)
{
	SERIESTYPE *new;
	SEQUENCETYPE *sequence;

	if(!series || !series->data)
		return(NULL);

	new = LALMalloc(sizeof(*new));
	sequence = `XLALCut'SEQUENCETYPE (series->data, first, length);
	if(!new || !sequence) {
		LALFree(new);
		`XLALDestroy'SEQUENCETYPE (sequence);
		return(NULL);
	}

	*new = *series;
	new->data = sequence;
	new->f0 += first * new->deltaF;

	return(new);
}


void `LALCut'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE **output,
	SERIESTYPE *input,
	size_t first,
	size_t length
)
{
	INITSTATUS(status, "`LALCut'SERIESTYPE", FREQUENCYSERIESC);

	ASSERT(output != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(input != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(input->data != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(first + length <= input->data->length, status, LAL_RANGE_ERR, LAL_RANGE_MSG);
	*output = `XLALCut'SERIESTYPE (input, first, length);
	ASSERT(*output != NULL, status, LAL_NOMEM_ERR, LAL_NOMEM_MSG);

	RETURN(status);
}


size_t `XLALShrink'SERIESTYPE (
	SERIESTYPE *series,
	size_t first,
	size_t length
)
{
	if(!series)
		return(0);

	series->f0 += first * series->deltaF;
	return(`XLALShrink'SEQUENCETYPE (series->data, first, length));
}


void `LALShrink'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE *series,
	size_t first,
	size_t length
)
{
	size_t new_length;

	INITSTATUS(status, "`LALShrink'SERIESTYPE", FREQUENCYSERIESC);

	ASSERT(series != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(series->data != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(first + length <= series->data->length, status, LAL_RANGE_ERR, LAL_RANGE_MSG);
	new_length = `XLALShrink'SERIESTYPE (series, first, length);
	ASSERT(new_length == length, status, LAL_FAIL_ERR, LAL_FAIL_MSG);
	RETURN(status);
}
