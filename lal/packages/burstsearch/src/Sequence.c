/* <lalVerbatim file="SequenceCV">
Author: Kipp Cannon <kipp@gravity.phys.uwm.edu>
$Id$
</lalVerbatim>
 */

#include <string.h>
#include <lal/Date.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALErrno.h>
#include <lal/LALError.h>
#include <lal/Sequence.h>

NRCSID(SEQUENCEC, "$Id$");


void XLALDestroyCOMPLEX8Sequence(
	COMPLEX8Sequence *sequence
)
{
	if(sequence) {
		LALFree(sequence->data);
	}
	LALFree(sequence);
}


void LALDestroyCOMPLEX8Sequence(
	LALStatus *status,
	COMPLEX8Sequence *sequence
)
{
	INITSTATUS(status, "LALDestroyCOMPLEX8Sequence", SEQUENCEC);

	XLALDestroyCOMPLEX8Sequence(sequence);

	RETURN(status);
}


void XLALDestroyREAL4Sequence(
	REAL4Sequence *sequence
)
{
	if(sequence) {
		LALFree(sequence->data);
	}
	LALFree(sequence);
}


void LALDestroyREAL4Sequence(
	LALStatus *status,
	REAL4Sequence *sequence
)
{
	INITSTATUS(status, "LALDestroyREAL4Sequence", SEQUENCEC);

	XLALDestroyREAL4Sequence(sequence);

	RETURN(status);
}


COMPLEX8Sequence *XLALCreateCOMPLEX8Sequence(
	size_t length
)
{
	COMPLEX8Sequence *new;
	COMPLEX8 *data;

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


void LALCreateCOMPLEX8Sequence(
	LALStatus *status,
	COMPLEX8Sequence **output,
	size_t length
)
{
	INITSTATUS(status, "LALCreateCOMPLEX8Sequence", SEQUENCEC);

	ASSERT(output != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);

	*output = XLALCreateCOMPLEX8Sequence(length);

	ASSERT(*output != NULL, status, LAL_FAIL_ERR, LAL_FAIL_MSG);

	RETURN(status);
}


REAL4Sequence *XLALCreateREAL4Sequence(
	size_t length
)
{
	REAL4Sequence *new;
	REAL4 *data;

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


void LALCreateREAL4Sequence(
	LALStatus *status,
	REAL4Sequence **output,
	size_t length
)
{
	INITSTATUS(status, "LALCreateREAL4Sequence", SEQUENCEC);

	ASSERT(output != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);

	*output = XLALCreateREAL4Sequence(length);

	ASSERT(*output != NULL, status, LAL_FAIL_ERR, LAL_FAIL_MSG);

	RETURN(status);
}


REAL4Sequence *XLALCutREAL4Sequence(
	REAL4Sequence *sequence,
	size_t first,
	size_t length
)
{
	REAL4Sequence *new;

	new = XLALCreateREAL4Sequence(length);
	if(!sequence || !sequence->data || !new) {
		XLALDestroyREAL4Sequence(new);
		return(NULL);
	}

	memcpy(new->data, sequence->data + first, length * sizeof(*new->data));

	return(new);
}


void LALCutREAL4Sequence(
	LALStatus *status,
	REAL4Sequence **output,
	REAL4Sequence *input,
	size_t first,
	size_t length
)
{
	INITSTATUS(status, "LALCutREAL4Sequence", SEQUENCEC);

	ASSERT(output != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);

	*output = XLALCutREAL4Sequence(input, first, length);

	ASSERT(*output != NULL, status, LAL_FAIL_ERR, LAL_FAIL_MSG);

	RETURN(status);
}


REAL4Sequence *XLALShrinkREAL4Sequence(
	REAL4Sequence *sequence,
	size_t first,
	size_t length
)
{
	if(sequence && sequence->data) {
		memmove(sequence->data, sequence->data + first, length * sizeof(*sequence->data));
		realloc(sequence->data, length * sizeof(*sequence->data));
		sequence->length = length;
	}

	return(sequence);
}


void LALShrinkREAL4Sequence(
	LALStatus *status,
	REAL4Sequence **sequence,
	size_t first,
	size_t length
)
{
	INITSTATUS(status, "LALCutREAL4Sequence", SEQUENCEC);

	ASSERT(sequence != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);

	*sequence = XLALShrinkREAL4Sequence(*sequence, first, length);

	ASSERT(*sequence != NULL, status, LAL_FAIL_ERR, LAL_FAIL_MSG);

	RETURN(status);
}
