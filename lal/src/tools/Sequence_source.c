#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define SEQUENCETYPE CONCAT2(DATATYPE,Sequence)

#define DFUNC CONCAT2(XLALDestroy,SEQUENCETYPE)
#define CFUNC CONCAT2(XLALCreate,SEQUENCETYPE)
#define XFUNC CONCAT2(XLALCut,SEQUENCETYPE)
#define CPFUNC CONCAT2(XLALCopy,SEQUENCETYPE)
#define SFUNC CONCAT2(XLALShift,SEQUENCETYPE)
#define RFUNC CONCAT2(XLALResize,SEQUENCETYPE)
#define SHFUNC CONCAT2(XLALShrink,SEQUENCETYPE)

#define TYPESUM CONCAT3(XLAL,DATATYPE,Sum)
#define TYPESUMSQ CONCAT3(XLAL,DATATYPE,SumSquares)

#define SEQUENCESUM CONCAT3(XLAL,SEQUENCETYPE,Sum)
#define SEQUENCESUMSQ CONCAT3(XLAL,SEQUENCETYPE,SumSquares)

void DFUNC (
	SEQUENCETYPE *sequence
)
{
#ifdef USE_ALIGNED_MEMORY_ROUTINES
	if(sequence)
		XLALFreeAligned(sequence->data);
#else
	if(sequence)
		XLALFree(sequence->data);
#endif
	XLALFree(sequence);
}


SEQUENCETYPE *CFUNC (
	size_t length
)
{
	SEQUENCETYPE *new;
	DATATYPE *data;

	new = XLALMalloc(sizeof(*new));

#ifdef USE_ALIGNED_MEMORY_ROUTINES
	data = XLALMallocAligned(length * sizeof(*data));
#else
	data = XLALMalloc(length * sizeof(*data));
#endif /*  USE_ALIGNED_MEMORY_ROUTINES */

	/* data == NULL is OK if length == 0 */
	if(!new || (length && !data)) {
		XLALFree(new);
		XLALFree(data);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	new->data = data;
	new->length = length;

	return new;
}


SEQUENCETYPE *XFUNC (
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
)
{
	SEQUENCETYPE *new = CFUNC (length);
	if(!new)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	memcpy(new->data, sequence->data + first, length * sizeof(*new->data));

	return new;
}


SEQUENCETYPE *CPFUNC (
	SEQUENCETYPE *sequence
)
{
	return XFUNC (sequence, 0, sequence->length);
}


void SFUNC (
	SEQUENCETYPE *sequence,
	int count
)
{
	if(!count)
		return;

	memshift(sequence->data, sequence->length * sizeof(*sequence->data), count * (int) sizeof(*sequence->data));
	if((size_t) labs(count) >= sequence->length)
		memset(sequence->data, 0, sequence->length * sizeof(*sequence->data));
	else if(count > 0)
		memset(sequence->data, 0, count * sizeof(*sequence->data));
	else
		memset(sequence->data + sequence->length + count, 0, -count * sizeof(*sequence->data));
}


/* FIXME: this function does not take care to move and zero the least
 * number of bytes possible.  A performance gain would be realized by being
 * more careful. */
SEQUENCETYPE *RFUNC (
	SEQUENCETYPE *sequence,
	int first,
	size_t length
)
{
	DATATYPE *new_data;

	if(length > sequence->length) {
		/* need to increase memory */

#ifdef USE_ALIGNED_MEMORY_ROUTINES
		new_data = XLALReallocAligned(sequence->data, length * sizeof(*sequence->data));
#else
		new_data = XLALRealloc(sequence->data, length * sizeof(*sequence->data));
#endif /* USE_ALIGNED_MEMORY_ROUTINES */

		if(new_data) {
			sequence->data = new_data;
			memset(sequence->data + sequence->length, 0, (length - sequence->length) * sizeof(*sequence->data));
			sequence->length = length;
			SFUNC (sequence, -first);
		} else
			XLAL_ERROR_NULL(XLAL_EFUNC);
	} else {
		/* do not need to increase memory */
		SFUNC (sequence, -first);
#ifdef USE_ALIGNED_MEMORY_ROUTINES
		new_data = XLALReallocAligned(sequence->data, length * sizeof(*sequence->data));
#else
		new_data = XLALRealloc(sequence->data, length * sizeof(*sequence->data));
#endif /* USE_ALIGNED_MEMORY_ROUTINES */
		if(new_data) {
			sequence->data = new_data;
			sequence->length = length;
		} else
			XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	return sequence;
}


SEQUENCETYPE *SHFUNC (
	SEQUENCETYPE *sequence,
	size_t first,
	size_t length
)
{
	return RFUNC (sequence, first, length);
}


DATATYPE TYPESUM (
	const DATATYPE *data,
	size_t first,
	size_t count
)
{
	DATATYPE sum = 0;

	for(data += first; count-- > 0; data++)
		sum += *data;

	return sum;
}


SQUAREDATATYPE TYPESUMSQ (
	const DATATYPE *data,
	size_t first,
	size_t count
)
{
	SQUAREDATATYPE sum = 0;

	for(data += first; count-- > 0; data++)
		/* clang cannot compile complex compound assignments yet
		 * radr://11224126 */
		/* sum += *data * *data; */
		sum = sum + (*data * *data);

	return sum;
}


DATATYPE SEQUENCESUM (
	const SEQUENCETYPE *sequence,
	size_t first,
	size_t count
)
{
	if(first >= sequence->length)
		return 0;
	if(first + count > sequence->length)
		count = sequence->length - first;
	return TYPESUM (sequence->data, first, count);
}


SQUAREDATATYPE SEQUENCESUMSQ (
	const SEQUENCETYPE *sequence,
	size_t first,
	size_t count
)
{
	if(first >= sequence->length)
		return 0;
	if(first + count > sequence->length)
		count = sequence->length - first;
	return TYPESUMSQ (sequence->data, first, count);
}

#undef SEQUENCETYPE
#undef DFUNC
#undef CFUNC
#undef XFUNC
#undef CPFUNC
#undef SFUNC
#undef RFUNC
#undef SHFUNC
#undef TYPESUM
#undef TYPESUMSQ
#undef SEQUENCESUM
#undef SEQUENCESUMSQ
