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
COMPLEX8TimeFrequencyPlane *XLALCreateTFPlane(
	TFPlaneParams *params,
	INT4 minFreqBins
)
/******** </lalVerbatim> ********/
{
	static const char *func = "XLALCreateTFPlane";
	COMPLEX8TimeFrequencyPlane *plane;
	COMPLEX8 *data;

	/* compute some extra TF plane params */
	params->timeBins = params->timeDuration / params->deltaT;
	params->freqBins = (params->fhigh - params->flow) / params->deltaF;

	/* Make sure that input parameters are reasonable */
	if((params->timeDuration <= 0) ||
	   (params->flow < 0) ||
	   (params->fhigh < params->flow) ||
	   (params->freqBins <= 0) ||
	   (params->deltaT <= 0.0) ||
	   (params->freqBins < minFreqBins))
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
	 * Fill some of the fields with nominal values pending the
	 * allocation of correct values for these fields.
	 */

	plane->epoch.gpsSeconds = 0;
	plane->epoch.gpsNanoSeconds = 0;
	plane->params = *params;
	plane->data = data;
	memset(data, 0, params->timeBins * params->freqBins * sizeof(*data));

	/* Normal exit */
	return(plane);
}
