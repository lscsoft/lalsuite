/******** <lalVerbatim file="ComputeFrequencySeriesCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/


#include <math.h>
#include <lal/LALRCSID.h>

NRCSID (COMPUTEFREQUENCYSERIESC, "$Id$");

#include <lal/RealFFT.h>
#include <lal/Sequence.h>
#include <lal/TFTransform.h>
#include <lal/Window.h>
#include <lal/XLALError.h>

/******** <lalVerbatim file="ComputeFrequencySeriesCP"> ********/
int
XLALComputeFrequencySeries(
	COMPLEX8FrequencySeries *freqSeries,
	REAL4TimeSeries *timeSeries,
	REAL4Window *window,
	REAL4FFTPlan *plan
)
/******** </lalVerbatim> ********/
{
	static const char *func = "XLALComputeFrequencySeries";
	REAL4Sequence *tmp;
	REAL4 fac;
	size_t i;
	size_t length = timeSeries->data->length;

	/* make sure sizes are reasonable and agree */
	if((freqSeries->data->length != length/2 + 1) ||
	   (window->data->length != length) ||
	   (timeSeries->deltaT <= 0.0) ||
	   (window->sumofsquares <= 0.0))
		XLAL_ERROR(func, XLAL_EINVAL);

	/* copy metadata into frequency series structure */
	freqSeries->epoch = timeSeries->epoch;
	freqSeries->f0 = timeSeries->f0;
	freqSeries->deltaF = 1.0 / (length * timeSeries->deltaT);

	/* create working copy of time series */
	tmp = XLALCutREAL4Sequence(timeSeries->data, 0, length);
	if(!tmp)
		XLAL_ERROR(func, XLAL_EFUNC);

	/* compute normalization factor */
	fac = sqrt(length / window->sumofsquares) * timeSeries->deltaT;

	/* compute windowed version of time series data */
	for(i = 0; i < length; i++)
		tmp->data[i] *= fac * window->data->data[i];

	/* compute the DFT */
	XLALREAL4ForwardFFT(freqSeries->data, tmp, plan);

	/* clean up */
	XLALDestroyREAL4Sequence(tmp);
	return(0);
}
