dnl $Id$
define(`SERIESTYPE',DATATYPE`FrequencySeries')
define(`SEQUENCETYPE',DATATYPE`Sequence')
/* <lalVerbatim> */
void `XLALDestroy'SERIESTYPE (
	SERIESTYPE *series
)
/* </lalVerbatim> */
{
	if(series)
		`XLALDestroy'SEQUENCETYPE (series->data);
	LALFree(series);
}


/* <lalVerbatim> */
void `LALDestroy'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE *series
)
/* </lalVerbatim> */
{
	INITSTATUS(status, "`LALDestroy'SERIESTYPE", FREQUENCYSERIESC);
	`XLALDestroy'SERIESTYPE (series);
	RETURN(status);
}


/* <lalVerbatim> */
SERIESTYPE *`XLALCreate'SERIESTYPE (
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaF,
	LALUnit sampleUnits,
	size_t length
)
/* </lalVerbatim> */
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


/* <lalVerbatim> */
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
/* </lalVerbatim> */
{
	INITSTATUS(status, "`LALCreate'SERIESTYPE", FREQUENCYSERIESC);
	ASSERT(output != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	*output = `XLALCreate'SERIESTYPE (name, epoch, f0, deltaF, sampleUnits, length);
	ASSERT(*output != NULL, status, LAL_FAIL_ERR, LAL_FAIL_MSG);
	RETURN(status);
}


/* <lalVerbatim> */
SERIESTYPE *`XLALCut'SERIESTYPE (
	SERIESTYPE *series,
	size_t first,
	size_t length
)
/* </lalVerbatim> */
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


/* <lalVerbatim> */
void `LALCut'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE **output,
	SERIESTYPE *input,
	size_t first,
	size_t length
)
/* </lalVerbatim> */
{
	INITSTATUS(status, "`LALCut'SERIESTYPE", FREQUENCYSERIESC);

	ASSERT(output != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(input != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(input->data != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	*output = `XLALCut'SERIESTYPE (input, first, length);
	ASSERT(*output != NULL, status, LAL_FAIL_ERR, LAL_FAIL_MSG);

	RETURN(status);
}


/* <lalVerbatim> */
size_t `XLALShrink'SERIESTYPE (
	SERIESTYPE *series,
	size_t first,
	size_t length
)
/* </lalVerbatim> */
{
	if(!series)
		return(0);

	series->f0 += first * series->deltaF;
	return(`XLALShrink'SEQUENCETYPE (series->data, first, length));
}


/* <lalVerbatim> */
void `LALShrink'SERIESTYPE (
	LALStatus *status,
	SERIESTYPE *series,
	size_t first,
	size_t length
)
/* </lalVerbatim> */
{
	INITSTATUS(status, "`LALShrink'SERIESTYPE", FREQUENCYSERIESC);

	ASSERT(series != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(series->data != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	`XLALShrink'SERIESTYPE (series, first, length);
	RETURN(status);
}
