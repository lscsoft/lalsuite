/* <lalVerbatim file="SequenceCV">
Author: Cannon, K. C.
$Id$
</lalVerbatim>
 */

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

	*output = XLALCreateCOMPLEX8Sequence(length);

	if(!*output)
		ABORT(status, LAL_FAIL_ERR, LAL_FAIL_MSG);

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

	*output = XLALCreateREAL4Sequence(length);

	if(!*output)
		ABORT(status, LAL_FAIL_ERR, LAL_FAIL_MSG);

	RETURN(status);
}


REAL4Sequence *XLALCutREAL4Sequence(
	REAL4Sequence *sequence,
	size_t first,
	size_t length
)
{
	REAL4Sequence *new;
	REAL4 *data;

	new = LALMalloc(sizeof(*new));
	data = LALMalloc(length * sizeof(*data));
	if(!sequence || !sequence->data || !new || !data) {
		LALFree(new);
		LALFree(data);
		return(NULL);
	}

	new->data = data;
	new->length = length;
	memcpy(new->data, sequence->data + first, length * sizeof(*data));

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

	*output = XLALCutREAL4Sequence(input, first, length);

	if(!*output)
		ABORT(status, LAL_FAIL_ERR, LAL_FAIL_MSG);

	RETURN(status);
}
