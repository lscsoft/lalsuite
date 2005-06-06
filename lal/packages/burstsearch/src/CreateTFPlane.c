/******** <lalVerbatim file="CreateTFPlaneCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/


#include <lal/LALRCSID.h>

NRCSID(CREATETFPLANEC, "$Id$");

#include <lal/LALStdlib.h>
#include <lal/TFTransform.h>
#include <lal/XLALError.h>


/******** <lalVerbatim file="CreateTFPlaneCP"> ********/
COMPLEX8TimeFrequencyPlane * XLALCreateTFPlane(
	TFPlaneParams * input
)
/******** </lalVerbatim> ********/
{
	static const char *func = "XLALCreateTFPlane";
	COMPLEX8TimeFrequencyPlane *tfp;
	TFPlaneParams *params;
	COMPLEX8 *data;

	/* Make sure that input parameters are reasonable */
	if((input->timeBins <= 0) || (input->freqBins <= 0) || (input->deltaT <= 0.0))
		XLAL_ERROR_NULL(func, XLAL_EDATA);

	/*  Assign memory */
	tfp = LALMalloc(sizeof(*tfp));
	params = LALMalloc(sizeof(*params));
	data = LALMalloc(input->timeBins * input->freqBins * sizeof(*data));
	if(!tfp || !params || !data) {
		LALFree(tfp);
		LALFree(params);
		LALFree(data);
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
	}

	/* 
	 *  Fill some of the fields with nominal values pending the allocation
	 *  of correct values for these fields.
	 */

	tfp->epoch.gpsSeconds = 0;
	tfp->epoch.gpsNanoSeconds = 0;
	tfp->params = params;
	tfp->data = data;
	*(tfp->params) = *input;

	/* Normal exit */
	return(tfp);
}
