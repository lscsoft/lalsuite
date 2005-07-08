/******** <lalVerbatim file="CreateTFPlaneCV"> ********
Author: Flanagan, E. and Cannon, K.
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>

NRCSID(CREATETFPLANEC, "$Id$");

#include <lal/LALStdlib.h>
#include <lal/TFTransform.h>
#include <lal/XLALError.h>


/******** <lalVerbatim file="CreateTFPlaneCP"> ********/
REAL4TimeFrequencyPlane *XLALCreateTFPlane(
	TFPlaneParams *params
)
/******** </lalVerbatim> ********/
{
	static const char *func = "XLALCreateTFPlane";
	REAL4TimeFrequencyPlane *plane;
	REAL4 *data;

	/* Make sure that input parameters are reasonable */
	if((params->flow < 0) ||
	   (params->timeBins <= 0) ||
	   (params->freqBins <= 0) ||
	   (params->deltaF <= 0.0) ||
	   (params->deltaT <= 0.0))
		XLAL_ERROR_NULL(func, XLAL_EDATA);

	/* Allocate memory */
	plane = LALMalloc(sizeof(*plane));
	data = LALMalloc(params->timeBins * params->freqBins * sizeof(*data));
	if(!plane || !data) {
		LALFree(plane);
		LALFree(data);
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
	}

	/* 
	 * Fill some of the fields with nominal values.
	 */

	plane->epoch.gpsSeconds = 0;
	plane->epoch.gpsNanoSeconds = 0;
	plane->params = *params;
	plane->data = data;
	memset(data, 0, params->timeBins * params->freqBins * sizeof(*data));

	/* Normal exit */
	return(plane);
}
