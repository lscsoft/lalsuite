#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALError.h>
#include <lal/TimeSeries.h>

NRCSID(TIMESERIESC, "$Id$");


void LALDestroyREAL4TimeSeries(
	LALStatus *status,
	REAL4TimeSeries *series
)
{
	INITSTATUS(status, "LALDestroyREAL4TimeSeries", TIMESERIESC);

	if(series) {
		if(series->data)
			LALFree(series->data->data);
		LALFree(series->data);
	}
	LALFree(series);

	RETURN(status);
}


void LALCutREAL4TimeSeries(
	LALStatus *status,
	REAL4TimeSeries **output,
	REAL4TimeSeries *input,
	size_t first_sample,
	size_t num_samples
)
{
	REAL4Sequence *sequence;
	REAL4 *data;

	INITSTATUS(status, "LALCutREAL4TimeSeries", TIMESERIESC);
	ATTATCHSTATUSPTR(status);

	*output = LALMalloc(sizeof(**output));
	sequence = LALMalloc(sizeof(*sequence));
	data = LALMalloc(num_samples * sizeof(*data));
	if(!input || !input->data->data || !*output || !sequence || !data) {
		LALFree(*output);
		LALFree(sequence);
		LALFree(data);
		*output = NULL;
		ABORT(status, LAL_FAIL_ERR, LAL_FAIL_MSG);
	}

	**output = *input;
	(*output)->data = sequence;
	sequence->data = data;
	sequence->length = num_samples;
	memcpy(sequence->data, input->data->data + first_sample, num_samples * sizeof(*data));

	LALAddFloatToGPS(status->statusPtr, &(*output)->epoch, &(*output)->epoch, first_sample * (*output)->deltaT);

	DETATCHSTATUSPTR(status);
	RETURN(status);
}
