/*
 * Copyright (C) 2009  Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <Python.h>
#include <lal/ComplexFFT.h>
#include <lal/RealFFT.h>
#include <lal/TimeFreqFFT.h>
#include <misc.h>
#include <datatypes/complex16fftplan.h>
#include <datatypes/complex16frequencyseries.h>
#include <datatypes/complex16timeseries.h>
#include <datatypes/real8fftplan.h>
#include <datatypes/real8frequencyseries.h>
#include <datatypes/real8timeseries.h>


#define MODULE_NAME "pylal.xlal.fft"


/*
 * ============================================================================
 *
 *                             Function Wrappers
 *
 * ============================================================================
 */


static PyObject *pylal_XLALCOMPLEX16FreqTimeFFT(PyObject *self, PyObject *args)
{
	pylal_COMPLEX16TimeSeries *tser;
	pylal_COMPLEX16FrequencySeries *fser;
	pylal_COMPLEX16FFTPlan *plan;

	if(!PyArg_ParseTuple(args, "O!O!O!", &pylal_COMPLEX16TimeSeries_Type, &tser, &pylal_COMPLEX16FrequencySeries_Type, &fser, &pylal_COMPLEX16FFTPlan_Type, &plan))
		return NULL;

	if(XLALCOMPLEX16FreqTimeFFT(tser->series, fser->series, plan->plan)) {
		pylal_set_exception_from_xlalerrno();
		return NULL;
	}

	Py_INCREF(Py_None);
	return Py_None;
}


static PyObject *pylal_XLALREAL8TimeFreqFFT(PyObject *self, PyObject *args)
{
	pylal_COMPLEX16FrequencySeries *fser;
	pylal_REAL8TimeSeries *tser;
	pylal_REAL8FFTPlan *plan;

	if(!PyArg_ParseTuple(args, "O!O!O!", &pylal_COMPLEX16FrequencySeries_Type, &fser, &pylal_REAL8TimeSeries_Type, &tser, &pylal_REAL8FFTPlan_Type, &plan))
		return NULL;

	if(XLALREAL8TimeFreqFFT(fser->series, tser->series, plan->plan)) {
		pylal_set_exception_from_xlalerrno();
		return NULL;
	}

	Py_INCREF(Py_None);
	return Py_None;
}


static PyObject *pylal_XLALREAL8FreqTimeFFT(PyObject *self, PyObject *args)
{
	pylal_REAL8TimeSeries *tser;
	pylal_COMPLEX16FrequencySeries *fser;
	pylal_REAL8FFTPlan *plan;

	if(!PyArg_ParseTuple(args, "O!O!O!", &pylal_REAL8TimeSeries_Type, &tser, &pylal_COMPLEX16FrequencySeries_Type, &fser, &pylal_REAL8FFTPlan_Type, &plan))
		return NULL;

	if(XLALREAL8FreqTimeFFT(tser->series, fser->series, plan->plan)) {
		pylal_set_exception_from_xlalerrno();
		return NULL;
	}

	Py_INCREF(Py_None);
	return Py_None;
}


static PyObject *pylal_XLALWhitenCOMPLEX16FrequencySeries(PyObject *self, PyObject *args)
{
	pylal_COMPLEX16FrequencySeries *fseries;
	pylal_REAL8FrequencySeries *psd;

	if(!PyArg_ParseTuple(args, "O!O!:XLALWhitenCOMPLEX16FrequencySeries", &pylal_COMPLEX16FrequencySeries_Type, &fseries, &pylal_REAL8FrequencySeries_Type, &psd))
		return NULL;

	if(!XLALWhitenCOMPLEX16FrequencySeries(fseries->series, psd->series)) {
		pylal_set_exception_from_xlalerrno();
		return NULL;
	}

	Py_INCREF(Py_None);
	return Py_None;
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static struct PyMethodDef module_methods[] = {
	{"XLALCOMPLEX16FreqTimeFFT", pylal_XLALCOMPLEX16FreqTimeFFT, METH_VARARGS, NULL},
	{"XLALREAL8TimeFreqFFT", pylal_XLALREAL8TimeFreqFFT, METH_VARARGS, NULL},
	{"XLALREAL8FreqTimeFFT", pylal_XLALREAL8FreqTimeFFT, METH_VARARGS, NULL},
	{"XLALWhitenCOMPLEX16FrequencySeries", pylal_XLALWhitenCOMPLEX16FrequencySeries, METH_VARARGS, NULL},
	{NULL,}
};


PyMODINIT_FUNC initfft(void)
{
	/* commented out to silence warning */
	/*PyObject *module = */Py_InitModule3(MODULE_NAME, module_methods, "Wrapper for LAL's fft package.");

	pylal_complex16fftplan_import();
	pylal_complex16frequencyseries_import();
	pylal_complex16timeseries_import();
	pylal_real8fftplan_import();
	pylal_real8frequencyseries_import();
	pylal_real8timeseries_import();
}
