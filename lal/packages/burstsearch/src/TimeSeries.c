/* <lalVerbatim file="TimeSeriesCV">
Author: Cannon, K. C.
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
#include <lal/TimeSeries.h>

NRCSID(TIMESERIESC, "$Id$");


void XLALDestroyCOMPLEX8TimeSeries(
	COMPLEX8TimeSeries *series
)
{
	if(series)
		XLALDestroyCOMPLEX8Sequence(series->data);
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
	if(series)
		XLALDestroyREAL4Sequence(series->data);
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

	new = LALMalloc(sizeof(*new));
	sequence = XLALCreateCOMPLEX8Sequence(length);
	if(!new || !sequence) {
		LALFree(new);
		XLALDestroyCOMPLEX8Sequence(sequence);
		return(NULL);
	}

	if(name)
		strncpy(new->name, name, LALNameLength);
	else
		*new->name = '\0';
	new->epoch = epoch;
	new->f0 = f0;
	new->deltaT = deltaT;
	new->sampleUnits = sampleUnits;
	new->data = sequence;

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

	ASSERT(output != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);

	*output = XLALCreateCOMPLEX8TimeSeries(name, epoch, f0, deltaT, sampleUnits, length);

	ASSERT(*output != NULL, status, LAL_FAIL_ERR, LAL_FAIL_MSG);

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

	new = LALMalloc(sizeof(*new));
	sequence = XLALCreateREAL4Sequence(length);
	if(!new || !sequence) {
		LALFree(new);
		XLALDestroyREAL4Sequence(sequence);
		return(NULL);
	}

	if(name)
		strncpy(new->name, name, LALNameLength);
	else
		*new->name = '\0';
	new->epoch = epoch;
	new->f0 = f0;
	new->deltaT = deltaT;
	new->sampleUnits = sampleUnits;
	new->data = sequence;

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

	ASSERT(output != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);

	*output = XLALCreateREAL4TimeSeries(name, epoch, f0, deltaT, sampleUnits, length);

	ASSERT(*output != NULL, status, LAL_FAIL_ERR, LAL_FAIL_MSG);

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

	if(!series || !series->data)
		return(NULL);

	new = LALMalloc(sizeof(*new));
	sequence = XLALCutREAL4Sequence(series->data, first_sample, num_samples);
	if(!new || !sequence) {
		LALFree(new);
		XLALDestroyREAL4Sequence(sequence);
		return(NULL);
	}

	*new = *series;
	new->data = sequence;
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

	ASSERT(output != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);

	*output = XLALCutREAL4TimeSeries(input, first_sample, num_samples);

	ASSERT(*output != NULL, status, LAL_FAIL_ERR, LAL_FAIL_MSG);

	RETURN(status);
}


REAL4TimeSeries *XLALShrinkREAL4TimeSeries(
	REAL4TimeSeries *series,
	size_t first_sample,
	size_t num_samples
)
{
	if(series) {
		series->data = XLALShrinkREAL4Sequence(series->data, first_sample, num_samples);
		series->epoch = XLALAddFloatToGPS(series->epoch, first_sample * series->deltaT);
	}

	return(series);
}


void LALShrinkREAL4TimeSeries(
	LALStatus *status,
	REAL4TimeSeries **series,
	size_t first_sample,
	size_t num_samples
)
{
	INITSTATUS(status, "LALCutREAL4TimeSeries", TIMESERIESC);

	ASSERT(series != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);

	*series = XLALShrinkREAL4TimeSeries(*series, first_sample, num_samples);

	ASSERT(*series != NULL, status, LAL_FAIL_ERR, LAL_FAIL_MSG);

	RETURN(status);
}
