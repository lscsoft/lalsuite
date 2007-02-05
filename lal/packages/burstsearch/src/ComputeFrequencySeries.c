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
	unsigned i;

	/* validate input */
	if(window->data->length != tseries->data->length)
		XLAL_ERROR_NULL(func, XLAL_EBADLEN);
	if(window->sumofsquares == 0.0)
		XLAL_ERROR_NULL(func, XLAL_EFPDIV0);

	/* create the frequency series, and a copy of the time series data */
	fseries = XLALCreateCOMPLEX8FrequencySeries(tseries->name, &tseries->epoch, tseries->f0, 1.0 / (tseries->data->length * tseries->deltaT), &lalDimensionlessUnit, tseries->data->length / 2 + 1);
	tmp = XLALCutREAL4Sequence(tseries->data, 0, tseries->data->length);
	if(!fseries || !tmp) {
		XLALDestroyCOMPLEX8FrequencySeries(fseries);
		XLALDestroyREAL4Sequence(tmp);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/* compute normalization factor */
	A = sqrt(tseries->data->length / window->sumofsquares) * tseries->deltaT;

	/* compute windowed version of time series data */
	for(i = 0; i < tseries->data->length; i++)
		tmp->data[i] *= A * window->data->data[i];

	/* compute the DFT */
	if(XLALREAL4ForwardFFT(fseries->data, tmp, plan)) {
		XLALDestroyCOMPLEX8FrequencySeries(fseries);
		XLALDestroyREAL4Sequence(tmp);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/* clean up */
	XLALDestroyREAL4Sequence(tmp);
	return(fseries);
}
