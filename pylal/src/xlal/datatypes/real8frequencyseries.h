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


#ifndef _PYLAL_XLAL_DATATYPES_REAL8FREQUENCYSERIES_H_
#define _PYLAL_XLAL_DATATYPES_REAL8FREQUENCYSERIES_H_


#include <Python.h>
#include <lal/LALDatatypes.h>
#include <lal/FrequencySeries.h>


#define PYLAL_REAL8FREQUENCYSERIES_MODULE_NAME "pylal.xlal.datatypes.real8frequencyseries"


/*
 * ============================================================================
 *
 *                                    Type
 *
 * ============================================================================
 */


static PyTypeObject *_pylal_REAL8FrequencySeries_Type = NULL;
#define pylal_REAL8FrequencySeries_Type (*_pylal_REAL8FrequencySeries_Type)


typedef struct {
	PyObject_HEAD
	PyObject *owner;
	REAL8FrequencySeries *series;
} pylal_REAL8FrequencySeries;


static PyObject *pylal_real8frequencyseries_import(void)
{
	PyObject *name = PyString_FromString(PYLAL_REAL8FREQUENCYSERIES_MODULE_NAME);
	PyObject *module = PyImport_Import(name);
	Py_DECREF(name);

	name = PyString_FromString("REAL8FrequencySeries");
	_pylal_REAL8FrequencySeries_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);
	Py_INCREF(&pylal_REAL8FrequencySeries_Type);
	Py_DECREF(name);

	return module;
}


PyObject *pylal_REAL8FrequencySeries_new(REAL8FrequencySeries *series, PyObject *owner)
{
	PyObject *empty_tuple = PyTuple_New(0);
	pylal_REAL8FrequencySeries *obj = (pylal_REAL8FrequencySeries *) PyType_GenericNew(&pylal_REAL8FrequencySeries_Type, empty_tuple, NULL);
	Py_DECREF(empty_tuple);
	if(!obj) {
		if(!owner)
			XLALDestroyREAL8FrequencySeries(series);
		return NULL;
	}
	if(owner)
		Py_INCREF(owner);
	obj->owner = owner;
	XLALDestroyREAL8FrequencySeries(obj->series);
	obj->series = series;
	return (PyObject *) obj;
}


#endif /* _PYLAL_XLAL_DATATYPES_REAL8FREQUENCYSERIES_H_ */
