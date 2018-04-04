/*
 * Copyright (C) 2010  Kipp Cannon
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
 *                        Python Wrapper For LALBurst
 *
 * ============================================================================
 */


#include <math.h>


#include <Python.h>
#include <numpy/arrayobject.h>


#include <lal/EPSearch.h>
#include <lal/LALDatatypes.h>
#include <lal/Sequence.h>


#include <misc.h>
#include <datatypes/complex16frequencyseries.h>
#include <datatypes/real8frequencyseries.h>


#define MODULE_NAME "pylal.xlal.lalburst"


/*
 * ============================================================================
 *
 *                                 Functions
 *
 * ============================================================================
 */


/*
 * pylal_XLALExcessPowerFilterInnerProduct()
 */


static PyObject *pylal_XLALExcessPowerFilterInnerProduct(PyObject *self, PyObject *args)
{
	pylal_COMPLEX16FrequencySeries *filter1;
	pylal_COMPLEX16FrequencySeries *filter2;
	PyArrayObject *correlation;
	REAL8Sequence corr;
	pylal_REAL8FrequencySeries *psd = NULL;
	double result;

	if(!PyArg_ParseTuple(args, "O!O!O!|O", &pylal_COMPLEX16FrequencySeries_Type, &filter1, &pylal_COMPLEX16FrequencySeries_Type, &filter2, &PyArray_Type, &correlation, &psd))
		return NULL;
	correlation = PyArray_GETCONTIGUOUS(correlation);
	if(!correlation)
		return NULL;
	corr.length = PyArray_DIM(correlation, 0);
	corr.data = PyArray_DATA(correlation);
	if((PyObject *) psd == Py_None)
		psd = NULL;
	else if(psd && !PyObject_TypeCheck((PyObject *) psd, &pylal_REAL8FrequencySeries_Type)) {
		Py_DECREF(correlation);
		PyErr_SetObject(PyExc_TypeError, (PyObject *) psd);
		return NULL;
	}

	result = XLALExcessPowerFilterInnerProduct(filter1->series, filter2->series, &corr, psd ? psd->series : NULL);
	if(XLAL_IS_REAL8_FAIL_NAN(result)) {
		Py_DECREF(correlation);
		pylal_set_exception_from_xlalerrno();
		return NULL;
	}

	Py_DECREF(correlation);
	return PyFloat_FromDouble(result);
}


/*
 * pylal_XLALCreateExcessPowerFilter()
 */


static PyObject *pylal_XLALCreateExcessPowerFilter(PyObject *self, PyObject *args)
{
	double flow;
	double width;
	pylal_REAL8FrequencySeries *psd = NULL;
	PyArrayObject *correlation = NULL;
	REAL8Sequence corr;
	COMPLEX16FrequencySeries *filter;

	if(!PyArg_ParseTuple(args, "ddO!O!", &flow, &width, &pylal_REAL8FrequencySeries_Type, &psd, &PyArray_Type, &correlation))
		return NULL;
	correlation = PyArray_GETCONTIGUOUS(correlation);
	if(!correlation)
		return NULL;
	corr.length = PyArray_DIM(correlation, 0);
	corr.data = PyArray_DATA(correlation);

	filter = XLALCreateExcessPowerFilter(flow, width, psd->series, &corr);
	Py_DECREF(correlation);
	if(!filter) {
		pylal_set_exception_from_xlalerrno();
		return NULL;
	}

	return pylal_COMPLEX16FrequencySeries_new(filter, NULL);
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static struct PyMethodDef methods[] = {
	{"XLALExcessPowerFilterInnerProduct", pylal_XLALExcessPowerFilterInnerProduct, METH_VARARGS, NULL},
	{"XLALCreateExcessPowerFilter", pylal_XLALCreateExcessPowerFilter, METH_VARARGS, NULL},
	{NULL,}
};


PyMODINIT_FUNC initlalburst(void)
{
	Py_InitModule3(MODULE_NAME, methods, "Wrapper for LALBurst package.");

	import_array()

	pylal_complex16frequencyseries_import();
	pylal_real8frequencyseries_import();
}
