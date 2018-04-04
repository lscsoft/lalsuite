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


#ifndef _PYLAL_XLAL_DATATYPES_REAL8TIMESERIES_H_
#define _PYLAL_XLAL_DATATYPES_REAL8TIMESERIES_H_


#include <Python.h>
#include <lal/LALDatatypes.h>
#include <lal/TimeSeries.h>


#define PYLAL_REAL8TIMESERIES_MODULE_NAME "pylal.xlal.datatypes.real8timeseries"


/*
 * ============================================================================
 *
 *                                    Type
 *
 * ============================================================================
 */


static PyTypeObject *_pylal_REAL8TimeSeries_Type = NULL;
#define pylal_REAL8TimeSeries_Type (*_pylal_REAL8TimeSeries_Type)


typedef struct {
	PyObject_HEAD
	PyObject *owner;
	REAL8TimeSeries *series;
} pylal_REAL8TimeSeries;


static PyObject *pylal_real8timeseries_import(void)
{
	PyObject *name = PyString_FromString(PYLAL_REAL8TIMESERIES_MODULE_NAME);
	PyObject *module = PyImport_Import(name);
	Py_DECREF(name);

	name = PyString_FromString("REAL8TimeSeries");
	_pylal_REAL8TimeSeries_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);
	Py_INCREF(&pylal_REAL8TimeSeries_Type);
	Py_DECREF(name);

	return module;
}


PyObject *pylal_REAL8TimeSeries_new(REAL8TimeSeries *series, PyObject *owner)
{
	PyObject *empty_tuple = PyTuple_New(0);
	pylal_REAL8TimeSeries *obj = (pylal_REAL8TimeSeries *) PyType_GenericNew(&pylal_REAL8TimeSeries_Type, empty_tuple, NULL);
	Py_DECREF(empty_tuple);
	if(!obj) {
		if(!owner)
			XLALDestroyREAL8TimeSeries(series);
		return NULL;
	}
	if(owner)
		Py_INCREF(owner);
	obj->owner = owner;
	XLALDestroyREAL8TimeSeries(obj->series);
	obj->series = series;
	return (PyObject *) obj;
}


#endif /* _PYLAL_XLAL_DATATYPES_REAL8TIMESERIES_H_ */
