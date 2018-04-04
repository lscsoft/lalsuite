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


#ifndef _PYLAL_XLAL_DATATYPES_COMPLEX16FREQUENCYSERIES_H_
#define _PYLAL_XLAL_DATATYPES_COMPLEX16FREQUENCYSERIES_H_


#include <Python.h>
#include <lal/LALDatatypes.h>
#include <lal/FrequencySeries.h>


#define PYLAL_COMPLEX16FREQUENCYSERIES_MODULE_NAME "pylal.xlal.datatypes.complex16frequencyseries"


/*
 * ============================================================================
 *
 *                                    Type
 *
 * ============================================================================
 */


static PyTypeObject *_pylal_COMPLEX16FrequencySeries_Type = NULL;
#define pylal_COMPLEX16FrequencySeries_Type (*_pylal_COMPLEX16FrequencySeries_Type)


typedef struct {
	PyObject_HEAD
	PyObject *owner;
	COMPLEX16FrequencySeries *series;
} pylal_COMPLEX16FrequencySeries;


static PyObject *pylal_complex16frequencyseries_import(void)
{
	PyObject *name = PyString_FromString(PYLAL_COMPLEX16FREQUENCYSERIES_MODULE_NAME);
	PyObject *module = PyImport_Import(name);
	Py_DECREF(name);

	name = PyString_FromString("COMPLEX16FrequencySeries");
	_pylal_COMPLEX16FrequencySeries_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);
	Py_INCREF(&pylal_COMPLEX16FrequencySeries_Type);
	Py_DECREF(name);

	return module;
}


PyObject *pylal_COMPLEX16FrequencySeries_new(COMPLEX16FrequencySeries *series, PyObject *owner)
{
	PyObject *empty_tuple = PyTuple_New(0);
	pylal_COMPLEX16FrequencySeries *obj = (pylal_COMPLEX16FrequencySeries *) PyType_GenericNew(&pylal_COMPLEX16FrequencySeries_Type, empty_tuple, NULL);
	Py_DECREF(empty_tuple);
	if(!obj) {
		if(!owner)
			XLALDestroyCOMPLEX16FrequencySeries(series);
		return NULL;
	}
	if(owner)
		Py_INCREF(owner);
	obj->owner = owner;
	XLALDestroyCOMPLEX16FrequencySeries(obj->series);
	obj->series = series;
	return (PyObject *) obj;
}


#endif /* _PYLAL_XLAL_DATATYPES_COMPLEX16FREQUENCYSERIES_H_ */
