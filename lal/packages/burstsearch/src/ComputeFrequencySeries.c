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


#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/EPSearch.h>	/* for our own prototypes */
#include <lal/Units.h>
#include <lal/Window.h>
#include <lal/XLALError.h>


#include <lal/LALRCSID.h>
NRCSID(COMPUTEFREQUENCYSERIESC, "$Id$");


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
	REAL4TimeSeries *tmp;
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
	 * the Fourier transform routine will set the frequency series
	 * metadata itself, so we don't have to get it correct but we do
	 * anyway.
	 */

	fseries = XLALCreateCOMPLEX8FrequencySeries(tseries->name, &tseries->epoch, tseries->f0, 1.0 / (tseries->data->length * tseries->deltaT), &tseries->sampleUnits, tseries->data->length / 2 + 1);
	tmp = XLALCutREAL4TimeSeries(tseries, 0, tseries->data->length);
	if(!fseries || !tmp) {
		XLALDestroyCOMPLEX8FrequencySeries(fseries);
		XLALDestroyREAL4TimeSeries(tmp);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/*
	 * apply normalized window to time series data;  emulate a
	 * rectangular window (all 1s) if none was supplied
	 */

	if(window) {
		const float A = sqrt(window->data->length / window->sumofsquares);
		for(i = 0; i < window->data->length; i++)
			tmp->data->data[i] *= A * window->data->data[i];
	}

	/*
	 * compute the DFT
	 */

	if(XLALREAL4TimeFreqFFT(fseries, tmp, plan)) {
		XLALDestroyCOMPLEX8FrequencySeries(fseries);
		XLALDestroyREAL4TimeSeries(tmp);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/*
	 * clean up
	 */

	XLALDestroyREAL4TimeSeries(tmp);

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
	REAL8TimeSeries *tmp;
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
	 * the Fourier transform routine will set the frequency series
	 * metadata itself, so we don't have to get it correct but we do
	 * anyway.
	 */

	fseries = XLALCreateCOMPLEX16FrequencySeries(tseries->name, &tseries->epoch, tseries->f0, 1.0 / (tseries->data->length * tseries->deltaT), &tseries->sampleUnits, tseries->data->length / 2 + 1);
	tmp = XLALCutREAL8TimeSeries(tseries, 0, tseries->data->length);
	if(!fseries || !tmp) {
		XLALDestroyCOMPLEX16FrequencySeries(fseries);
		XLALDestroyREAL8TimeSeries(tmp);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/*
	 * apply normalized window to time series data;  emulate a
	 * rectangular window (all 1s) if none was supplied
	 */

	if(window) {
		const double A = sqrt(window->data->length / window->sumofsquares);
		for(i = 0; i < window->data->length; i++)
			tmp->data->data[i] *= A * window->data->data[i];
	}

	/*
	 * compute the DFT
	 */

	if(XLALREAL8TimeFreqFFT(fseries, tmp, plan)) {
		XLALDestroyCOMPLEX16FrequencySeries(fseries);
		XLALDestroyREAL8TimeSeries(tmp);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/*
	 * clean up
	 */

	XLALDestroyREAL8TimeSeries(tmp);

	/*
	 * done
	 */

	return fseries;
}
