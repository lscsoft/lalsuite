/*
 * $Id$
 *
 * Copyright (C) 2007  Kipp Cannon and Flanagan, E
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#include <math.h>


#include <lal/FrequencySeries.h>
#include <lal/RealFFT.h>
#include <lal/Sequence.h>
#include <lal/EPSearch.h>	/* for our own prototypes */
#include <lal/Units.h>
#include <lal/Window.h>
#include <lal/XLALError.h>


/**
 * This function accepts a time series and a window function, and computes
 * and returns the Fourier transform of the windowed time series.  The
 * output is equal to the left-hand-side of (7) in LIGO-T010095-00-Z with
 * the normalization factor (9) applied.  If window is NULL, then a
 * rectangular window (all 1s) is used.
 */


COMPLEX8FrequencySeries *XLALWindowedREAL4ForwardFFT(
	const REAL4TimeSeries *tseries,
	const REAL4Window *window,
	const REAL4FFTPlan *plan
)
{
	static const char func[] = "XLALWindowedREAL4ForwardFFT";
	COMPLEX8FrequencySeries *fseries;
	REAL4Sequence *tmp;
	unsigned i;

	/*
	 * validate input
	 */

	if(((window != NULL) && (window->data->length != tseries->data->length)) ||
	   (tseries->data->length == 0))
		XLAL_ERROR_NULL(func, XLAL_EBADLEN);
	if(((window != NULL) && (window->sumofsquares == 0.0)) ||
	   (tseries->deltaT == 0.0))
		XLAL_ERROR_NULL(func, XLAL_EFPDIV0);

	/*
	 * create the frequency series, and a copy of the time series data
	 */

	fseries = XLALCreateCOMPLEX8FrequencySeries(tseries->name, &tseries->epoch, tseries->f0, 1.0 / (tseries->data->length * tseries->deltaT), &tseries->sampleUnits, tseries->data->length / 2 + 1);
	tmp = XLALCutREAL4Sequence(tseries->data, 0, tseries->data->length);
	if(!fseries || !tmp) {
		XLALDestroyCOMPLEX8FrequencySeries(fseries);
		XLALDestroyREAL4Sequence(tmp);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/*
	 * set the frequency series' units
	 */

	XLALUnitMultiply(&fseries->sampleUnits, &fseries->sampleUnits, &lalSecondUnit);

	/*
	 * apply normalized window to time series data;  emulate a
	 * rectangular window (all 1s) if none was supplied
	 */

	if(window) {
		const float A = tseries->deltaT / sqrt(window->sumofsquares / window->data->length);
		for(i = 0; i < tmp->length; i++)
			tmp->data[i] *= A * window->data->data[i];
	} else {
		const float A = tseries->deltaT;
		for(i = 0; i < tmp->length; i++)
			tmp->data[i] *= A;
	}

	/*
	 * compute the DFT
	 */

	if(XLALREAL4ForwardFFT(fseries->data, tmp, plan)) {
		XLALDestroyCOMPLEX8FrequencySeries(fseries);
		XLALDestroyREAL4Sequence(tmp);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/*
	 * clean up
	 */

	XLALDestroyREAL4Sequence(tmp);

	/*
	 * done
	 */

	return fseries;
}


/**
 * This function accepts a time series and a window function, and computes
 * and returns the Fourier transform of the windowed time series.  The
 * output is equal to the left-hand-side of (7) in LIGO-T010095-00-Z with
 * the normalization factor (9) applied.  If window is NULL then a
 * rectanguar window (all 1s) is used.
 */


COMPLEX16FrequencySeries *XLALWindowedREAL8ForwardFFT(
	const REAL8TimeSeries *tseries,
	const REAL8Window *window,
	const REAL8FFTPlan *plan
)
{
	static const char func[] = "XLALWindowedREAL8ForwardFFT";
	COMPLEX16FrequencySeries *fseries;
	REAL8Sequence *tmp;
	unsigned i;

	/*
	 * validate input
	 */

	if(((window != NULL) && (window->data->length != tseries->data->length)) ||
	   (tseries->data->length == 0))
		XLAL_ERROR_NULL(func, XLAL_EBADLEN);
	if(((window != NULL) && (window->sumofsquares == 0.0)) ||
	   (tseries->deltaT == 0.0))
		XLAL_ERROR_NULL(func, XLAL_EFPDIV0);

	/*
	 * create the frequency series, and a copy of the time series data
	 */

	fseries = XLALCreateCOMPLEX16FrequencySeries(tseries->name, &tseries->epoch, tseries->f0, 1.0 / (tseries->data->length * tseries->deltaT), &tseries->sampleUnits, tseries->data->length / 2 + 1);
	tmp = XLALCutREAL8Sequence(tseries->data, 0, tseries->data->length);
	if(!fseries || !tmp) {
		XLALDestroyCOMPLEX16FrequencySeries(fseries);
		XLALDestroyREAL8Sequence(tmp);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/*
	 * set the frequency series' units
	 */

	XLALUnitMultiply(&fseries->sampleUnits, &fseries->sampleUnits, &lalSecondUnit);

	/*
	 * apply normalized window to time series data;  emulate a
	 * rectangular window (all 1s) if none was supplied
	 */

	if(window) {
		const double A = sqrt(tseries->data->length / window->sumofsquares) * tseries->deltaT;
		for(i = 0; i < tmp->length; i++)
			tmp->data[i] *= A * window->data->data[i];
	} else {
		const double A = tseries->deltaT;
		for(i = 0; i < tmp->length; i++)
			tmp->data[i] *= A;
	}

	/*
	 * compute the DFT
	 */

	if(XLALREAL8ForwardFFT(fseries->data, tmp, plan)) {
		XLALDestroyCOMPLEX16FrequencySeries(fseries);
		XLALDestroyREAL8Sequence(tmp);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/*
	 * clean up
	 */

	XLALDestroyREAL8Sequence(tmp);

	/*
	 * done
	 */

	return fseries;
}
