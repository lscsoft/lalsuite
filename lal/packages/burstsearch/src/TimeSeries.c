#include <lal/Date.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALError.h>
#include <lal/TimeSeries.h>

NRCSID(TIMESERIESC, "$Id$");


void XLALDestroyREAL4TimeSeries(
	REAL4TimeSeries *series
)
{
	if(series) {
		if(series->data)
			LALFree(series->data->data);
		LALFree(series->data);
	}
	LALFree(series);
}


void LALDestroyREAL4TimeSeries(
	LALStatus *status,
	REAL4TimeSeries *series
)
{
	INITSTATUS(status, "LALDestroyREAL4TimeSeries", TIMESERIESC);

	XLALDestroyREAL4TimeSeries(series);

	RETURN(status);
}


REAL4TimeSeries *XLALCutREAL4TimeSeries(
	REAL4TimeSeries *series,
	size_t first_sample,
	size_t num_samples
)
{
	REAL4TimeSeries *new;
	REAL4Sequence *sequence;
	REAL4 *data;

	new = LALMalloc(sizeof(*new));
	sequence = LALMalloc(sizeof(*sequence));
	data = LALMalloc(num_samples * sizeof(*data));
	if(!series || !series->data || !series->data->data || !new || !sequence || !data) {
		LALFree(new);
		LALFree(sequence);
		LALFree(data);
		return(NULL);
	}

	*new = *series;
	new->data = sequence;
	sequence->data = data;
	sequence->length = num_samples;
	memcpy(sequence->data, series->data->data + first_sample, num_samples * sizeof(*data));

	new->epoch = XLALAddFloatToGPS(new->epoch, first_sample * new->deltaT);

	return(new);
}


void LALCutREAL4TimeSeries(
	LALStatus *status,
	REAL4TimeSeries **output,
	REAL4TimeSeries *input,
	size_t first_sample,
	size_t num_samples
)
{
	INITSTATUS(status, "LALCutREAL4TimeSeries", TIMESERIESC);

	*output = XLALCutREAL4TimeSeries(input, first_sample, num_samples);

	if(!*output)
		ABORT(status, LAL_FAIL_ERR, LAL_FAIL_MSG);

	RETURN(status);
}
