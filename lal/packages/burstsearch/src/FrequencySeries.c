#include <string.h>

#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALErrno.h>
#include <lal/LALError.h>
#include <lal/FrequencySeries.h>

NRCSID(FREQUENCYSERIESC, "$Id$");


void XLALDestroyCOMPLEX8FrequencySeries(
	COMPLEX8FrequencySeries *series
)
{
	if(series) {
		if(series->data)
			LALFree(series->data->data);
		LALFree(series->data);
	}
	LALFree(series);
}


void LALDestroyCOMPLEX8FrequencySeries(
	LALStatus *status,
	COMPLEX8FrequencySeries *series
)
{
	INITSTATUS(status, "LALDestroyCOMPLEX8FrequencySeries", FREQUENCYSERIESC);

	XLALDestroyCOMPLEX8FrequencySeries(series);

	RETURN(status);
}


void XLALDestroyREAL4FrequencySeries(
	REAL4FrequencySeries *series
)
{
	if(series) {
		if(series->data)
			LALFree(series->data->data);
		LALFree(series->data);
	}
	LALFree(series);
}


void LALDestroyREAL4FrequencySeries(
	LALStatus *status,
	REAL4FrequencySeries *series
)
{
	INITSTATUS(status, "LALDestroyREAL4FrequencySeries", FREQUENCYSERIESC);

	XLALDestroyREAL4FrequencySeries(series);

	RETURN(status);
}


COMPLEX8FrequencySeries *XLALCreateCOMPLEX8FrequencySeries(
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaF,
	LALUnit sampleUnits,
	size_t length
)
{
	COMPLEX8FrequencySeries *new;
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
	new->deltaF = deltaF;
	new->sampleUnits = sampleUnits;
	new->data = sequence;
	sequence->data = data;
	sequence->length = length;

	return(new);
}


void LALCreateCOMPLEX8FrequencySeries(
	LALStatus *status,
	COMPLEX8FrequencySeries **output,
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaF,
	LALUnit sampleUnits,
	size_t length
)
{
	INITSTATUS(status, "LALCreateCOMPLEX8FrequencySeries", FREQUENCYSERIESC);

	*output = XLALCreateCOMPLEX8FrequencySeries(name, epoch, f0, deltaF, sampleUnits, length);

	if(!*output)
		ABORT(status, LAL_FAIL_ERR, LAL_FAIL_MSG);

	RETURN(status);
}


REAL4FrequencySeries *XLALCreateREAL4FrequencySeries(
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaF,
	LALUnit sampleUnits,
	size_t length
)
{
	REAL4FrequencySeries *new;
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
	new->deltaF = deltaF;
	new->sampleUnits = sampleUnits;
	new->data = sequence;
	sequence->data = data;
	sequence->length = length;

	return(new);
}


void LALCreateREAL4FrequencySeries(
	LALStatus *status,
	REAL4FrequencySeries **output,
	CHAR *name,
	LIGOTimeGPS epoch,
	REAL8 f0,
	REAL8 deltaF,
	LALUnit sampleUnits,
	size_t length
)
{
	INITSTATUS(status, "LALCreateREAL4FrequencySeries", FREQUENCYSERIESC);

	*output = XLALCreateREAL4FrequencySeries(name, epoch, f0, deltaF, sampleUnits, length);

	if(!*output)
		ABORT(status, LAL_FAIL_ERR, LAL_FAIL_MSG);

	RETURN(status);
}
