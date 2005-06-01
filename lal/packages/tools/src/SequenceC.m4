dnl $Id$
define(`SEQUENCETYPE',DATATYPE`Sequence')
void `XLALDestroy'SEQUENCETYPE (
	SEQUENCETYPE *sequence
)
{
	if(sequence)
		LALFree(sequence->data);
	LALFree(sequence);
}


void `LALDestroy'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE *sequence
)
{
	INITSTATUS(status, "`LALDestroy'SEQUENCETYPE", SEQUENCEC);
	`XLALDestroy'SEQUENCETYPE (sequence);
	RETURN(status);
}


SEQUENCETYPE *`XLALCreate'SEQUENCETYPE (
	size_t length
)
{
	static const char *func = "`XLALCreate'SEQUENCETYPE";
	SEQUENCETYPE *new;
	DATATYPE *data;

	new = LALMalloc(sizeof(*new));
	data = LALMalloc(length * sizeof(*data));
	if(!new || !data) {
		LALFree(new);
		LALFree(data);
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
	}

	new->data = data;
	new->length = length;

	return(new);
}


void `LALCreate'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE **output,
	size_t length
)
{
	INITSTATUS(status, "`LALCreate'SEQUENCETYPE", SEQUENCEC);
	ASSERT(output != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	*output = `XLALCreate'SEQUENCETYPE (length);
	if(*output == NULL) {
		XLALClearErrno();
		ABORT(status, LAL_NOMEM_ERR, LAL_NOMEM_MSG);
	}
	RETURN(status);
}


SEQUENCETYPE *`XLALCut'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
)
{
	static const char *func = "`XLALCut'SEQUENCETYPE";
	SEQUENCETYPE *new = NULL;

	if(sequence && sequence->data) {
		new = `XLALCreate'SEQUENCETYPE (length);
		if(!new)
			XLAL_ERROR_NULL(func, XLAL_EFUNC);
		memcpy(new->data, sequence->data + first, length * sizeof(*new->data));
	}

	return(new);
}


void `LALCut'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE **output,
	SEQUENCETYPE *input,
	size_t first,
	size_t length
)
{
	INITSTATUS(status, "`LALCut'SEQUENCETYPE", SEQUENCEC);
	ASSERT(output != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(input != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(first + length <= input->length, status, LAL_RANGE_ERR, LAL_RANGE_MSG);
	*output = `XLALCut'SEQUENCETYPE (input, first, length);
	if(*output == NULL) {
		XLALClearErrno();
		ABORT(status, LAL_NOMEM_ERR, LAL_NOMEM_MSG);
	}
	RETURN(status);
}


SEQUENCETYPE *`XLALCopy'SEQUENCETYPE (
	SEQUENCETYPE *sequence
)
{
	return(`XLALCut'SEQUENCETYPE (sequence, 0, sequence->length));
}


void `LALCopy'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE **output,
	SEQUENCETYPE *input
)
{
	INITSTATUS(status, "`LALCopy'SEQUENCETYPE", SEQUENCEC);
	ASSERT(output != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(input != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	*output = `XLALCopy'SEQUENCETYPE (input);
	if(*output == NULL) {
		XLALClearErrno();
		ABORT(status, LAL_NOMEM_ERR, LAL_NOMEM_MSG);
	}
	RETURN(status);
}


static void `Shift'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	int count
)
{
	memshift(sequence->data, sequence->length * sizeof(*sequence->data), count * (int) sizeof(*sequence->data));
}


void `XLALShift'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	int count
)
{
	if(!sequence || !sequence->data || !count)
		return;

	`Shift'SEQUENCETYPE (sequence, count);
	if((size_t) labs(count) >= sequence->length)
		memset(sequence->data, 0, sequence->length * sizeof(*sequence->data));
	else if(count > 0)
		memset(sequence->data, 0, count * sizeof(*sequence->data));
	else
		memset(sequence->data + sequence->length + count, 0, -count * sizeof(*sequence->data));
}


void `LALShift'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE *sequence,
	int count
)
{
	INITSTATUS(status, "`LALShift'SEQUENCETYPE", SEQUENCEC);
	ASSERT(sequence != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	`XLALShift'SEQUENCETYPE (sequence, count);
	RETURN(status);
}


/* FIXME: this function does not take care to move the least number of bytes
 * possible.  A performance gain would be realized by being more careful. */
/* FIXME: this function does not conform to the XLAL error reporting
 * convention. */
size_t `XLALResize'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	int first,
	size_t length
)
{
	if(!sequence || !sequence->data)
		return(0);

	if(length > sequence->length) {
		/* need to increase memory */
		sequence->data = LALRealloc(sequence->data, length * sizeof(*sequence->data));

		if(sequence->data)
			`Shift'SEQUENCETYPE (sequence, -first);
		else
			length = 0;
	} else {
		/* do not need to increase memory */
		`Shift'SEQUENCETYPE (sequence, -first);
		sequence->data = LALRealloc(sequence->data, length * sizeof(*sequence->data));
		if(!sequence->data)
			length = 0;
	}

	return(sequence->length = length);
}


void `LALResize'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE *sequence,
	int first,
	size_t length
)
{
	size_t new_length;

	INITSTATUS(status, "`LALResize'SEQUENCETYPE", SEQUENCEC);
	ASSERT(sequence != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	new_length = `XLALResize'SEQUENCETYPE (sequence, first, length);
	if(new_length != length) {
		XLALClearErrno();
		ABORT(status, LAL_FAIL_ERR, LAL_FAIL_MSG);
	}
	RETURN(status);
}


size_t `XLALShrink'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
)
{
	return(`XLALResize'SEQUENCETYPE (sequence, first, length));
}


void `LALShrink'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
)
{
	size_t new_length;

	INITSTATUS(status, "`LALShrink'SEQUENCETYPE", SEQUENCEC);
	ASSERT(sequence != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(first + length <= sequence->length, status, LAL_RANGE_ERR, LAL_RANGE_MSG);
	new_length = `XLALShrink'SEQUENCETYPE (sequence, first, length);
	if(new_length != length) {
		XLALClearErrno();
		ABORT(status, LAL_FAIL_ERR, LAL_FAIL_MSG);
	}
	RETURN(status);
}
