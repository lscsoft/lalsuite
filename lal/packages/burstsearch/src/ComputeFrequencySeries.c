/******** <lalVerbatim file="ComputeFrequencySeriesCV"> ********
Author: Flanagan, E. and Cannon, K.
$Id$
********* </lalVerbatim> ********/


#include <math.h>
#include <lal/LALRCSID.h>

NRCSID (COMPUTEFREQUENCYSERIESC, "$Id$");

#include <lal/FrequencySeries.h>
#include <lal/RealFFT.h>
#include <lal/Sequence.h>
#include <lal/TFTransform.h>
#include <lal/Units.h>
#include <lal/Window.h>
#include <lal/XLALError.h>

/******** <lalVerbatim file="ComputeFrequencySeriesCP"> ********/
COMPLEX8FrequencySeries *
XLALComputeFrequencySeries(
	const REAL4TimeSeries *tseries,
	const REAL4Window *window,
	const REAL4FFTPlan *plan
)
/******** </lalVerbatim> ********/
{
	/*
	 * This function accepts a time series and a window function, and
	 * computes and returns the Fourier transform of the windowed time
	 * series.  The output is equal to the left-hand-side of (3) in
	 * LIGO-T010095-00-Z with the normalization factor (5) applied.
	 */

	const char func[] = "XLALComputeFrequencySeries";
	COMPLEX8FrequencySeries *fseries;
	REAL4Sequence *tmp;
	REAL4 A;
	size_t i;
	size_t length = tseries->data->length;

	/* make sure sizes are reasonable and agree */
	if((window->data->length != length) ||
	   (tseries->deltaT <= 0.0) ||
	   (window->sumofsquares <= 0.0))
		XLAL_ERROR_NULL(func, XLAL_EINVAL);

	/* create the frequency series, and a copy of the time series data */
	fseries = XLALCreateCOMPLEX8FrequencySeries(tseries->name, &tseries->epoch, tseries->f0, 1.0 / (length * tseries->deltaT), &lalDimensionlessUnit, length / 2 + 1);
	tmp = XLALCutREAL4Sequence(tseries->data, 0, length);
	if(!fseries || !tmp)
		XLAL_ERROR_NULL(func, XLAL_EFUNC);

	/* compute normalization factor */
	A = sqrt(length / window->sumofsquares) * tseries->deltaT;

	/* compute windowed version of time series data */
	for(i = 0; i < length; i++)
		tmp->data[i] *= A * window->data->data[i];

	/* compute the DFT */
	XLALREAL4ForwardFFT(fseries->data, tmp, plan);

	/* clean up */
	XLALDestroyREAL4Sequence(tmp);
	return(fseries);
}
