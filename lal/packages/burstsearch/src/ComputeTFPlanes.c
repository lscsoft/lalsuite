/******** <lalVerbatim file="ComputeTFPlanesCV">
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>

NRCSID(COMPUTETFPLANESC, "$Id$");

#include <lal/ExcessPower.h>
#include <lal/LALErrno.h>
#include <lal/XLALError.h>

#define TRUE 1
#define FALSE 0

/******** <lalVerbatim file="ComputeTFPlanesCP"> ********/
int XLALComputeTFPlanes(
	TFTiling * tfTiling,
	COMPLEX8FrequencySeries * freqSeries,
	UINT4 windowShift,
	REAL4 * norm,
	REAL4FrequencySeries * psd
)
/******** </lalVerbatim> ********/
{
	static const char *func = "XLALComputeTFPlanes";
	HorizontalTFTransformIn transformparams;

	/* check ranges */
	if (tfTiling->numPlanes <= 0)
		XLAL_ERROR(func, XLAL_EDOM);

	/* setup input structure for computing TF transform */
	transformparams.startT = 0;	/* not used for horizontal transforms */
	transformparams.dftParams = *tfTiling->dftParams;
	transformparams.windowShift = windowShift;

	/* Compute TF transform */
	if (XLALFreqSeriesToTFPlane(*tfTiling->tfp, freqSeries, &transformparams, norm, psd))
		XLAL_ERROR(func, XLAL_EFUNC);

	/* set flags saying TF planes have been computed, but not EP or sorted */
	tfTiling->planesComputed = TRUE;
	tfTiling->excessPowerComputed = FALSE;
	tfTiling->tilesSorted = FALSE;

	/* normal exit */
	return(0);
}
