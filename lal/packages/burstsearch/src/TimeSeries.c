#include <lal/Date.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALError.h>
#include <lal/TimeSeries.h>

NRCSID(TIMESERIESC, "$Id$");


void XLALDestroyCOMPLEX8TimeSeries(
	COMPLEX8TimeSeries *series
)
{
	if(series) {
		if(series->data)
			LALFree(series->data->data);
		LALFree(series->data);
	}
	LALFree(series);
}


void LALDestroyCOMPLEX8TimeSeries(
	LALStatus *status,
	COMPLEX8TimeSeries *series
)
{
	INITSTATUS(status, "LALDestroyCOMPLEX8TimeSeries", TIMESERIESC);

	XLALDestroyCOMPLEX8TimeSeries(series);

	RETURN(status);
}


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


COMPLEX8TimeSeries *XLALCreateCOMPLEX8TimeSeries(
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaT,
	LALUnit sampleUnits,
	size_t length
)
{
	COMPLEX8TimeSeries *new;
	COMPLEX8Sequence *sequence;
	COMPLEX8 *data;

	new = LALMalloc(sizeof(*new));
	sequence = LALMalloc(sizeof(*sequence));
	data = LALMalloc(length * sizeof(*data));
	if(!new || !sequence || !data) {
		LALFree(new);
		LALFree(sequence);
		LALFree(data);
		return(NULL);
	}

	strncpy(new->name, name, LALNameLength);
	new->epoch = epoch;
	new->f0 = f0;
	new->deltaT = deltaT;
	new->sampleUnits = sampleUnits;
	new->data = sequence;
	sequence->data = data;
	sequence->length = length;

	return(new);
}


void LALCreateCOMPLEX8TimeSeries(
	LALStatus *status,
	COMPLEX8TimeSeries **output,
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaT,
	LALUnit sampleUnits,
	size_t length
)
{
	INITSTATUS(status, "LALCreateCOMPLEX8TimeSeries", TIMESERIESC);

	*output = XLALCreateCOMPLEX8TimeSeries(name, epoch, f0, deltaT, sampleUnits, length);

	if(!*output)
		ABORT(status, LAL_FAIL_ERR, LAL_FAIL_MSG);

	RETURN(status);
}


REAL4TimeSeries *XLALCreateREAL4TimeSeries(
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaT,
	LALUnit sampleUnits,
	size_t length
)
{
	REAL4TimeSeries *new;
	REAL4Sequence *sequence;
	REAL4 *data;

	new = LALMalloc(sizeof(*new));
	sequence = LALMalloc(sizeof(*sequence));
	data = LALMalloc(length * sizeof(*data));
	if(!new || !sequence || !data) {
		LALFree(new);
		LALFree(sequence);
		LALFree(data);
		return(NULL);
	}

	strncpy(new->name, name, LALNameLength);
	new->epoch = epoch;
	new->f0 = f0;
	new->deltaT = deltaT;
	new->sampleUnits = sampleUnits;
	new->data = sequence;
	sequence->data = data;
	sequence->length = length;

	return(new);
}


void LALCreateREAL4TimeSeries(
	LALStatus *status,
	REAL4TimeSeries **output,
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaT,
	LALUnit sampleUnits,
	size_t length
)
{
	INITSTATUS(status, "LALCreateREAL4TimeSeries", TIMESERIESC);

	*output = XLALCreateREAL4TimeSeries(name, epoch, f0, deltaT, sampleUnits, length);

	if(!*output)
		ABORT(status, LAL_FAIL_ERR, LAL_FAIL_MSG);

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
