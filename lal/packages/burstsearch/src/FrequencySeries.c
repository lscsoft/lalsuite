/* <lalVerbatim file="FrequencySeriesCV">
Author: Cannon, K. C.
$Id$
</lalVerbatim>
 */

#include <string.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALErrno.h>
#include <lal/LALError.h>
#include <lal/FrequencySeries.h>
#include <lal/Sequence.h>

NRCSID(FREQUENCYSERIESC, "$Id$");


void XLALDestroyCOMPLEX8FrequencySeries(
	COMPLEX8FrequencySeries *series
)
{
	if(series)
		XLALDestroyCOMPLEX8Sequence(series->data);
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
	if(series)
		XLALDestroyREAL4Sequence(series->data);
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

	new = LALMalloc(sizeof(*new));
	sequence = XLALCreateCOMPLEX8Sequence(length);
	if(!new || !sequence) {
		LALFree(new);
		XLALDestroyCOMPLEX8Sequence(sequence);
		return(NULL);
	}

	if(new->name)
		strncpy(new->name, name, LALNameLength);
	else
		*new->name = '\0';
	new->epoch = epoch;
	new->f0 = f0;
	new->deltaF = deltaF;
	new->sampleUnits = sampleUnits;
	new->data = sequence;

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

	ASSERT(output != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);

	*output = XLALCreateCOMPLEX8FrequencySeries(name, epoch, f0, deltaF, sampleUnits, length);

	ASSERT(*output != NULL, status, LAL_FAIL_ERR, LAL_FAIL_MSG);

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

	new = LALMalloc(sizeof(*new));
	sequence = XLALCreateREAL4Sequence(length);
	if(!new || !sequence) {
		LALFree(new);
		XLALDestroyREAL4Sequence(sequence);
		return(NULL);
	}

	if(new->name)
		strncpy(new->name, name, LALNameLength);
	else
		*new->name = '\0';
	new->epoch = epoch;
	new->f0 = f0;
	new->deltaF = deltaF;
	new->sampleUnits = sampleUnits;
	new->data = sequence;

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

	ASSERT(output != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);

	*output = XLALCreateREAL4FrequencySeries(name, epoch, f0, deltaF, sampleUnits, length);

	ASSERT(*output != NULL, status, LAL_FAIL_ERR, LAL_FAIL_MSG);

	RETURN(status);
}
