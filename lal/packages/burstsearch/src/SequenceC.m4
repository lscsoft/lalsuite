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
	SEQUENCETYPE *new;
	DATATYPE *data;

	new = LALMalloc(sizeof(*new));
	data = LALMalloc(length * sizeof(*data));
	if(!new || !data) {
		LALFree(new);
		LALFree(data);
		return(NULL);
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
	ASSERT(*output != NULL, status, LAL_FAIL_ERR, LAL_FAIL_MSG);
	RETURN(status);
}


SEQUENCETYPE *`XLALCut'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
)
{
	SEQUENCETYPE *new = NULL;

	if(sequence && sequence->data) {
		if(first >= sequence->length)
			length = 0;
		else if(first + length > sequence->length)
			length = sequence->length - first;
		new = `XLALCreate'SEQUENCETYPE (length);
		if(new && length)
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
	*output = `XLALCut'SEQUENCETYPE (input, first, length);
	ASSERT(*output != NULL, status, LAL_FAIL_ERR, LAL_FAIL_MSG);
	RETURN(status);
}


size_t `XLALShrink'SEQUENCETYPE (
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
)
{
	if(!sequence || !sequence->data)
		return(0);

	if(first >= sequence->length)
		length = 0;
	else {
		if(first + length > sequence->length)
			length = sequence->length - first;
		memmove(sequence->data, sequence->data + first, length * sizeof(*sequence->data));
	}

	sequence->data = LALRealloc(sequence->data, length * sizeof(*sequence->data));
	if(sequence->data)
		sequence->length = length;
	else
		sequence->length = 0;

	return(sequence->length);
}


void `LALShrink'SEQUENCETYPE (
	LALStatus *status,
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
)
{
	INITSTATUS(status, "`LALShrink'SEQUENCETYPE", SEQUENCEC);
	ASSERT(sequence != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	`XLALShrink'SEQUENCETYPE (sequence, first, length);
	RETURN(status);
}
